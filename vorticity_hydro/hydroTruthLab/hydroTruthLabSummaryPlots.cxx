// hydroTruthLabSummaryPlots.cxx
//
// Compares the runs of hydroTruthLab across configurations. For every comparison of the runs registry: the integrals, lobe features and
// curves of each run against the x of the comparison, and each run minus the first one, jackknifed jointly over the events the two share
// when they start from the same initial conditions.
//
// Usage:
//   ./hydroTruthLabSummaryPlots.exe <outputFile> --compare=<name>,<x>,<perConfigFile>,<perConfigFile>[,...] [--compare=...]
//   Normally compiled and launched by runHydroTruthLab.sh after every run, with one --compare per compare line of hydroTruthLab_runs.conf.
//
// Needs: ROOT, the PerConfig files of the runs, and the producer files they were made from (paths in their Config/).
//
// Definitions and the output reference: README.md in this folder ("SummaryPlots"). What it shares with PerConfig lives in
// hydroTruthLabCommon.h.

#include <TDirectory.h>
#include <TFile.h>
#include <TKey.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TStyle.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "hydroTruthLabCommon.h"

// ======================================================================================================================================
// Constants
// ======================================================================================================================================

constexpr int kNRingEstimators = 4; // R (and R_z) by its four estimators, the first (second) block of the integral tables
const std::array<std::string, kNCurveSamples> kSampleTitles = {"bullet", "no_bullet", "bullet - no_bullet"};

// Sets compared run against run. The delta-scan is left out: its integrals average over delta, which vanishes identically (half-turn
// identity), whatever the run
const std::array<int, 3> kDifferenceSets = {kSetTrue, kSetMixing, kSetSmudge};
const std::array<std::string, 3> kDifferenceSetKeys = {"True", "Mixed", "Smudged"};

// Producer options that fix a trigger set event by event. Two runs that agree on them give identical per-event values, so their difference
// over the events they share vanishes exactly: a closure of the whole chain. Mixing also depends on the event range (partner rule), and
// smudging on its seed. The index is the trigger set
const std::vector<std::string> kTrueTriggerKeys = {"sampleDir", "jetFile", "trigger", "inputName", "ptMin", "ptMax", "absYMax", "lambdaMass", "grid", "columnToPolarization"};
const std::array<std::vector<std::string>, kNTriggerSets> kSetExtraKeys = {{{}, {"deltaScanSteps"}, {"mixDepth", "firstEvent", "nEvents"}, {"smudgeSigma", "smudgeMode", "smudgeDraws", "seed"}}};

// Seed spread: the producer entries that may differ between runs that differ only by their seed, and the PerConfig options that must agree
const std::vector<std::string> kSeedFreeEntries = {"outputFile", "seed", "warnings"};
const std::vector<std::string> kPerConfigOptionKeys = {"paired", "alphaLambda", "nLambdaPrecision", "centralityLabel", "initialConditions"};

// ======================================================================================================================================
// Runs
// ======================================================================================================================================

void printUsage()
{
  std::cout << "Usage: ./hydroTruthLabSummaryPlots.exe <outputFile> --compare=<name>,<x>,<perConfigFile>,<perConfigFile>[,...] [--compare=...]\n"
            << "  name   folder of the comparison in the output file\n"
            << "  x      what labels the runs: an entry of the PerConfig Config/ (centralityLabel, runLabel, ...) or of its bullet producer's\n"
            << "         (smudgeSigma, seed, ...)\n"
            << "See README.md for the definitions." << std::endl;
}

/// \brief One comparison of the runs registry.
struct ComparisonSpec {
  std::string name;
  std::string x;
  std::vector<std::string> perConfigFiles;
};

/// \brief One run of a comparison: the metadata of its PerConfig file, the configurations of its producers, their event trees and the
///        estimates, computed with the code of PerConfig (so the numbers are those of its PerConfig file).
struct RunData {
  std::string runLabel;
  std::string xLabel;
  std::string perConfigFile;
  std::string initialConditions;
  bool hasNoBullet = false;
  bool paired = true;
  double alphaLambda = 0.;
  double nLambdaPrecision = 0.;
  std::map<std::string, std::string> perConfigEntries;                     // Config/ of the PerConfig file
  std::array<std::map<std::string, std::string>, kNSamples> producerConfig; // its Config/ProducerBullet/ and Config/ProducerNoBullet/
  ProducerSetup setup;
  std::array<SampleData, kNSamples> samples;
  std::array<std::map<int, int>, kNSamples> rowOfEvent;
  std::array<std::vector<int>, kNSamples> pairRows;
  std::array<SetResult, kNTriggerSets> results;
  std::array<std::array<EstimateSet, kNCurveSamples>, kNTriggerSets> leakage;
};

/// \brief Every TNamed of a directory, as name -> title (empty when the directory is missing).
std::map<std::string, std::string> readNamedEntries(TFile* file, const std::string& path)
{
  std::map<std::string, std::string> entries;
  TDirectory* directory = file->GetDirectory(path.c_str());
  if (directory == nullptr)
    return entries;
  TIter nextKey(directory->GetListOfKeys());
  while (auto* key = static_cast<TKey*>(nextKey())) {
    if (std::string(key->GetClassName()) != "TNamed")
      continue;
    std::unique_ptr<TNamed> named(key->ReadObject<TNamed>());
    entries[named->GetName()] = named->GetTitle();
  }
  return entries;
}

/// \brief Value of an entry, with a hard stop when it is missing (a PerConfig file older than this macro, for instance).
std::string requiredEntry(const std::map<std::string, std::string>& entries, const std::string& name, const std::string& where)
{
  const auto found = entries.find(name);
  if (found == entries.end())
    throw std::runtime_error(where + " has no " + name + " (rerun PerConfig)");
  return found->second;
}

/// \brief Value of an entry, empty when missing.
std::string entryOf(const std::map<std::string, std::string>& entries, const std::string& name)
{
  const auto found = entries.find(name);
  return found == entries.end() ? "" : found->second;
}

/// \brief Path of a producer file of a run: the one PerConfig read or, when the campaign folder has moved since, the file of the same name
///        next to the PerConfig file.
std::string producerPath(const std::string& storedPath, const std::string& perConfigFile)
{
  if (std::filesystem::is_regular_file(storedPath))
    return storedPath;
  const std::filesystem::path nextToPerConfig = std::filesystem::path(perConfigFile).parent_path() / std::filesystem::path(storedPath).filename();
  if (std::filesystem::is_regular_file(nextToPerConfig))
    return nextToPerConfig.string();
  throw std::runtime_error("producer file " + storedPath + " not found, neither there nor next to " + perConfigFile);
}

/// \brief Reads a run: its PerConfig metadata, the event trees of its producers, and its estimates. The warnings of the estimates repeat
///        those of the run's PerConfig log, so only their number is reported here.
RunData loadRun(const std::string& perConfigFile, const std::string& prefix, std::vector<std::string>& warnings)
{
  RunData run;
  run.perConfigFile = perConfigFile;
  std::unique_ptr<TFile> perConfig(TFile::Open(perConfigFile.c_str(), "READ"));
  if (!perConfig || perConfig->IsZombie())
    throw std::runtime_error("cannot open " + perConfigFile);
  run.perConfigEntries = readNamedEntries(perConfig.get(), "Config");
  run.producerConfig[kSampleBullet] = readNamedEntries(perConfig.get(), "Config/ProducerBullet");
  run.producerConfig[kSampleNoBullet] = readNamedEntries(perConfig.get(), "Config/ProducerNoBullet");
  run.runLabel = requiredEntry(run.perConfigEntries, "runLabel", perConfigFile);
  run.initialConditions = requiredEntry(run.perConfigEntries, "initialConditions", perConfigFile);
  run.paired = requiredEntry(run.perConfigEntries, "paired", perConfigFile) == "1";
  run.alphaLambda = std::stod(requiredEntry(run.perConfigEntries, "alphaLambda", perConfigFile));
  run.nLambdaPrecision = std::stod(requiredEntry(run.perConfigEntries, "nLambdaPrecision", perConfigFile));
  const std::string noBulletEntry = requiredEntry(run.perConfigEntries, "noBulletFile", perConfigFile);
  run.hasNoBullet = noBulletEntry != "none";
  if (!run.hasNoBullet)
    warnings.push_back(prefix + run.runLabel + " has no no_bullet sample, so it is missing from every no_bullet and difference plot");

  const std::array<std::string, kNSamples> storedPaths = {requiredEntry(run.perConfigEntries, "bulletFile", perConfigFile), noBulletEntry};
  std::array<std::unique_ptr<TFile>, kNSamples> producerFiles;
  for (int sample = 0; sample < kNSamples; ++sample) {
    if (sample == kSampleNoBullet && !run.hasNoBullet)
      continue;
    const std::string path = producerPath(storedPaths[sample], perConfigFile);
    producerFiles[sample].reset(TFile::Open(path.c_str(), "READ"));
    if (!producerFiles[sample] || producerFiles[sample]->IsZombie())
      throw std::runtime_error("cannot open " + path);
  }
  run.setup = readProducerSetup(producerFiles[kSampleBullet].get());
  for (int sample = 0; sample < kNSamples; ++sample) {
    if (!producerFiles[sample])
      continue;
    run.samples[sample] = readSample(producerFiles[sample].get(), run.runLabel + ", " + kSampleLabels[sample], run.setup.deltaPhiRange.nBins, run.setup.nDelta, run.setup.scanZeroBin, run.setup.mixDepth, run.setup.smudgeSigma > 0.);
    run.rowOfEvent[sample] = rowsByEvent(run.samples[sample], run.runLabel + ", " + kSampleLabels[sample]);
  }
  if (run.hasNoBullet)
    run.pairRows = pairRowsOf(run.samples[kSampleBullet], run.samples[kSampleNoBullet], run.runLabel + ", ");
  std::vector<std::string> estimateWarnings;
  run.results = computeSetResults(run.samples, run.hasNoBullet, run.paired, run.pairRows, run.setup.deltaPhiCentres, run.runLabel + ", ", estimateWarnings);
  run.leakage = computeLeakage(run.results, run.hasNoBullet, run.paired, run.runLabel + ", ", estimateWarnings);
  if (!estimateWarnings.empty())
    warnings.push_back(prefix + run.runLabel + ": " + std::to_string(estimateWarnings.size()) + " warnings of its estimates, the same as in its PerConfig log");
  return run;
}

/// \brief Value of the x of a comparison in a run: an entry of its PerConfig Config/ or, failing that, of its bullet producer's. Numbers are
///        printed short ("0.100000" -> "0.1").
std::string xValueOf(const RunData& run, const std::string& x)
{
  std::string value = entryOf(run.perConfigEntries, x);
  if (value.empty())
    value = entryOf(run.producerConfig[kSampleBullet], x);
  if (value.empty())
    throw std::runtime_error(run.runLabel + ": no entry " + x + " in the Config/ of its PerConfig file nor of its bullet producer");
  char* parseEnd = nullptr;
  const double number = std::strtod(value.c_str(), &parseEnd);
  if (parseEnd == value.c_str() || *parseEnd != '\0')
    return value;
  std::ostringstream shortForm;
  shortForm << number;
  return shortForm.str();
}

bool hasSample(const RunData& run, int sample) { return sample == kSampleBullet || run.hasNoBullet; }

bool setInAllRuns(const std::vector<RunData>& runs, int set)
{
  return std::all_of(runs.begin(), runs.end(), [set](const RunData& run) { return run.results[set].present; });
}

bool setInSomeRun(const std::vector<RunData>& runs, int set)
{
  return std::any_of(runs.begin(), runs.end(), [set](const RunData& run) { return run.results[set].present; });
}

// ======================================================================================================================================
// Estimates across runs
// ======================================================================================================================================

EstimatePoint pointAt(const EstimateSet& estimates, int index) { return {estimates.value[index], estimates.error[index]}; }

/// \brief One estimate of a run, undefined when the run lacks the set or the sample.
EstimatePoint pointOf(const RunData& run, int set, int sample, int index)
{
  if (!run.results[set].present || !hasSample(run, sample))
    return {};
  return pointAt(run.results[set].estimates[sample], index);
}

/// \brief Event IDs of the events two runs share, ascending: for one sample, the events both of its trees have; for the difference, the
///        pairs both runs have.
std::vector<int> sharedEventIds(const RunData& first, const RunData& other, int sample)
{
  std::array<std::vector<int>, 2> ids;
  for (int which = 0; which < 2; ++which) {
    const RunData& run = which == 0 ? first : other;
    if (sample != kSampleDifference) {
      ids[which] = run.samples[sample].eventId;
    } else {
      for (const int row : run.pairRows[kSampleBullet])
        ids[which].push_back(run.samples[kSampleBullet].eventId[row]);
    }
    std::sort(ids[which].begin(), ids[which].end());
  }
  std::vector<int> shared;
  std::set_intersection(ids[0].begin(), ids[0].end(), ids[1].begin(), ids[1].end(), std::back_inserter(shared));
  return shared;
}

/// \brief Rows of the given events in the tree of a sample, in the order of the events.
std::vector<int> rowsOfEvents(const std::map<int, int>& rowOfEvent, const std::vector<int>& eventIds)
{
  std::vector<int> rows;
  rows.reserve(eventIds.size());
  for (const int eventId : eventIds)
    rows.push_back(rowOfEvent.at(eventId));
  return rows;
}

/// \brief Stage-1 estimates of one sample of a run over the given events, in their order, with one replica per event. For the difference,
///        the paired bullet - no_bullet vectors, whose replica e drops pair e.
Replicated restrictedStage1(const RunData& run, int set, int sample, const std::vector<int>& eventIds, const EstimateLayout& layout)
{
  if (sample != kSampleDifference)
    return replicateStage1(selectEvents(run.samples[sample].sets[set], rowsOfEvents(run.rowOfEvent[sample], eventIds)));
  const Replicated bullet = replicateStage1(selectEvents(run.samples[kSampleBullet].sets[set], rowsOfEvents(run.rowOfEvent[kSampleBullet], eventIds)));
  const Replicated noBullet = replicateStage1(selectEvents(run.samples[kSampleNoBullet].sets[set], rowsOfEvents(run.rowOfEvent[kSampleNoBullet], eventIds)));
  const DifferenceReplicas pairs = differenceReplicas(bullet, noBullet, true, layout);
  return {pairs.full, pairs.paired};
}

/// \brief A run minus the first run of the comparison, for one set and one sample: a stage-1 vector with its errors.
struct RunDifference {
  bool computed = false;
  bool joint = false;    // jackknifed jointly over shared events, or independent runs
  int nSharedEvents = 0; // joint only: events, or pairs for the difference sample
  EstimateSet estimates;
};

/// \brief Run minus the first run. With the same initial conditions, both runs are restricted to the events (pairs) they share, in the same
///        order, and replica e drops event e from both: the fluctuation of the initial conditions cancels in the difference, which a sum of
///        the two errors in quadrature would count twice. Otherwise the runs are independent: the difference of their estimates over all
///        their events, with the two jackknife errors in quadrature. The difference sample is paired across runs only when both runs pair
///        their own samples (--paired=1 in PerConfig).
RunDifference differenceToFirst(const RunData& first, const RunData& other, int set, int sample, const std::string& label, std::vector<std::string>& warnings)
{
  RunDifference difference;
  if (!hasSample(first, sample) || !hasSample(other, sample))
    return difference;
  difference.computed = true;
  const EstimateLayout& layout = first.results[set].layout;
  const bool sameInitialConditions = first.initialConditions != "none" && first.initialConditions == other.initialConditions;
  difference.joint = sameInitialConditions && (sample != kSampleDifference || (first.paired && other.paired));
  if (!difference.joint) {
    const EstimateSet& firstEstimates = first.results[set].estimates[sample];
    const EstimateSet& otherEstimates = other.results[set].estimates[sample];
    difference.estimates.value = differenceOfVectors(otherEstimates.value, firstEstimates.value);
    difference.estimates.error.assign(difference.estimates.value.size(), kNoValue);
    for (std::size_t index = 0; index < difference.estimates.value.size(); ++index)
      if (firstEstimates.error[index] != kNoValue && otherEstimates.error[index] != kNoValue)
        difference.estimates.error[index] = std::hypot(firstEstimates.error[index], otherEstimates.error[index]);
    return difference;
  }
  const std::vector<int> eventIds = sharedEventIds(first, other, sample);
  difference.nSharedEvents = eventIds.size();
  if (eventIds.size() < 2) {
    warnings.push_back(label + ": " + std::to_string(eventIds.size()) + " shared events, too few for the jackknife, so no difference");
    difference.computed = false;
    return difference;
  }
  const Replicated firstStage1 = restrictedStage1(first, set, sample, eventIds, layout);
  const Replicated otherStage1 = restrictedStage1(other, set, sample, eventIds, layout);
  Replicated delta;
  delta.full = differenceOfVectors(otherStage1.full, firstStage1.full);
  for (std::size_t event = 0; event < eventIds.size(); ++event)
    delta.replicas.push_back(differenceOfVectors(otherStage1.replicas[event], firstStage1.replicas[event]));
  difference.estimates = estimatesFromReplicas(delta.full, delta.replicas, label, warnings);
  return difference;
}

/// \brief True when two runs define a trigger set identically, event by event, for one sample (both producers for the difference).
bool sameSetDefinition(const RunData& first, const RunData& other, int set, int sample)
{
  const std::vector<int> producers = sample == kSampleDifference ? std::vector<int>{kSampleBullet, kSampleNoBullet} : std::vector<int>{sample};
  for (const int producer : producers) {
    for (const std::vector<std::string>* keys : {&kTrueTriggerKeys, &kSetExtraKeys[set]})
      for (const std::string& key : *keys)
        if (entryOf(first.producerConfig[producer], key) != entryOf(other.producerConfig[producer], key))
          return false;
  }
  return true;
}

/// \brief True when two runs differ only by their seed: same producer entries (but kSeedFreeEntries) and same PerConfig options.
bool differOnlyBySeed(const RunData& first, const RunData& other)
{
  if (first.hasNoBullet != other.hasNoBullet)
    return false;
  for (int sample = 0; sample < kNSamples; ++sample) {
    std::map<std::string, std::string> firstEntries = first.producerConfig[sample];
    std::map<std::string, std::string> otherEntries = other.producerConfig[sample];
    for (const std::string& key : kSeedFreeEntries) {
      firstEntries.erase(key);
      otherEntries.erase(key);
    }
    if (firstEntries != otherEntries)
      return false;
  }
  for (const std::string& key : kPerConfigOptionKeys)
    if (entryOf(first.perConfigEntries, key) != entryOf(other.perConfigEntries, key))
      return false;
  return true;
}

// ======================================================================================================================================
// Drawing
// ======================================================================================================================================

/// \brief A table over the runs (or the runs but the first): one labelled bin per run, empty where the estimate is undefined.
Part runTable(TDirectory* directory, const std::string& name, const std::string& title, const std::vector<std::string>& runLabels, const std::vector<EstimatePoint>& points, double scale, std::vector<std::string>& warnings)
{
  std::vector<double> values, errors;
  for (const EstimatePoint& point : points) {
    values.push_back(point.value);
    errors.push_back(point.error);
  }
  return tableHistogram(directory, name, title, runLabels, values, errors, blockIndices(0, runLabels.size()), scale, warnings);
}

/// \brief Canvas of several estimates over the runs, one table per estimate ([estimate][run]), side by side in each bin.
void drawRunTables(const OutputFolder& folder, const std::string& name, const std::string& title, const std::string& axisTitles, const std::vector<std::string>& runLabels, const std::vector<std::string>& estimateLabels, const std::vector<std::vector<EstimatePoint>>& points, double scale, double referenceY, std::vector<std::string>& warnings)
{
  std::vector<DrawEntry> entries;
  for (std::size_t estimate = 0; estimate < estimateLabels.size(); ++estimate) {
    const Part part = runTable(folder.parts, "h" + name.substr(1) + "_" + std::to_string(estimate), estimateLabels[estimate] + axisTitles, runLabels, points[estimate], scale, warnings);
    entries.push_back({part, estimateLabels[estimate], kPaletteLooks[estimate], kStyleMarker});
  }
  drawOverlay(folder.canvases, name, title, entries, referenceY);
}

/// \brief Canvas of one curve of one set and sample, one band per run that has it.
void drawRunCurves(const OutputFolder& folder, const std::vector<RunData>& runs, int set, int curve, int sample, const std::string& name, const std::string& title, std::vector<std::string>& warnings)
{
  std::vector<DrawEntry> entries;
  const AxisRange& range = runs.front().setup.deltaPhiRange;
  const std::string symbol = kFeatureCurveSymbols[curve];
  for (std::size_t run = 0; run < runs.size(); ++run) {
    if (!runs[run].results[set].present || !hasSample(runs[run], sample)) // logged when the run was read
      continue;
    const EstimateSet& estimates = runs[run].results[set].estimates[sample];
    const Part part = curveHistogram(folder.parts, "h" + name.substr(1) + "_run" + std::to_string(run), "<" + symbol + "> (" + kSetLabels[set] + ", " + kSampleTitles[sample] + ", " + runs[run].xLabel + ", jackknife errors);" + range.title + ";<" + symbol + ">", range, estimates.value, estimates.error, runs[run].results[set].layout.curve(curve, 0), 1., warnings);
    entries.push_back({part, runs[run].xLabel, kPaletteLooks[run], kStyleBand});
  }
  drawOverlay(folder.canvases, name, title, entries, 0.);
}

// ======================================================================================================================================
// main
// ======================================================================================================================================

int main(int argc, char** argv)
{
  try {
    // ==================================================================================================================================
    // 1. Options
    // ==================================================================================================================================
    std::string outputFileName;
    std::vector<ComparisonSpec> comparisons;
    for (int iArg = 1; iArg < argc; ++iArg) {
      const std::string argument = argv[iArg];
      if (argument == "--help" || argument == "-h") {
        printUsage();
        return 0;
      }
      if (argument.rfind("--compare=", 0) == 0) {
        std::vector<std::string> fields;
        std::istringstream specStream(argument.substr(std::string("--compare=").size()));
        std::string field;
        while (std::getline(specStream, field, ','))
          fields.push_back(field);
        if (fields.size() < 4)
          throw std::runtime_error("expected --compare=<name>,<x>,<perConfigFile>,<perConfigFile>[,...], got " + argument);
        comparisons.push_back({fields[0], fields[1], std::vector<std::string>(fields.begin() + 2, fields.end())});
      } else if (argument.rfind("--", 0) == 0) {
        throw std::runtime_error("unknown option " + argument + " (see --help)");
      } else if (outputFileName.empty()) {
        outputFileName = argument;
      } else {
        throw std::runtime_error("expected one positional argument, <outputFile>, got a second: " + argument);
      }
    }
    if (outputFileName.empty() || comparisons.empty()) {
      printUsage();
      throw std::runtime_error("expected <outputFile> and at least one --compare");
    }
    for (std::size_t index = 0; index < comparisons.size(); ++index)
      for (std::size_t other = 0; other < index; ++other)
        if (comparisons[index].name == comparisons[other].name)
          throw std::runtime_error("comparison " + comparisons[index].name + " given twice");

    std::cout << "========================================================" << std::endl;
    std::cout << "  hydroTruthLabSummaryPlots" << std::endl;
    std::cout << "  outputFile = " << outputFileName << std::endl;
    for (const ComparisonSpec& spec : comparisons)
      std::cout << "  " << spec.name << ": x = " << spec.x << ", " << spec.perConfigFiles.size() << " runs" << std::endl;
    std::cout << "========================================================" << std::endl;

    std::vector<std::string> warnings; // printed at the end and stored in Config/: nothing is skipped silently

    // ==================================================================================================================================
    // 2. Output file
    // ==================================================================================================================================
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLegendBorderSize(0);
    std::unique_ptr<TFile> outputFile(TFile::Open(outputFileName.c_str(), "RECREATE"));
    if (!outputFile || outputFile->IsZombie())
      throw std::runtime_error("cannot create " + outputFileName);
    TDirectory* configDir = outputFile->mkdir("Config");

    // ==================================================================================================================================
    // 3. Comparisons, one folder each
    // ==================================================================================================================================
    for (const ComparisonSpec& spec : comparisons) {
      const std::string prefix = spec.name + ": ";
      std::cout << std::endl << "Comparison " << spec.name << std::endl;

      // --------------------------------------------------------------------------------------------------------------------------------
      // 3a. Runs. A run without a PerConfig file (not made yet, or failed) or with unreadable inputs is left out
      // --------------------------------------------------------------------------------------------------------------------------------
      std::vector<RunData> runs;
      for (const std::string& file : spec.perConfigFiles) {
        if (!std::filesystem::is_regular_file(file)) {
          warnings.push_back(prefix + "no PerConfig file " + file + " (run not made, or failed), so that run is left out");
          continue;
        }
        try {
          runs.push_back(loadRun(file, prefix, warnings));
          runs.back().xLabel = xValueOf(runs.back(), spec.x);
        } catch (const std::exception& error) {
          warnings.push_back(prefix + file + " left out: " + error.what());
        }
      }
      if (runs.size() < 2 || runs.size() > kPaletteLooks.size()) {
        warnings.push_back(prefix + std::to_string(runs.size()) + " runs, while a comparison needs 2 to " + std::to_string(kPaletteLooks.size()) + ": skipped");
        continue;
      }
      // The runs are compared bin by bin: every run must have the Delta phi axis of the first
      const AxisRange& firstRange = runs.front().setup.deltaPhiRange;
      bool sameAxes = true;
      for (const RunData& run : runs)
        if (run.setup.deltaPhiRange.nBins != firstRange.nBins || run.setup.deltaPhiRange.low != firstRange.low || run.setup.deltaPhiRange.high != firstRange.high)
          sameAxes = false;
      if (!sameAxes) {
        warnings.push_back(prefix + "the runs have different Delta phi axes, so they cannot be compared bin by bin: skipped");
        continue;
      }
      // Run labels are <dataset>_<variation>, "_" being reserved by the coordinator for that split. As x, they keep only the part that
      // tells the runs apart
      if (spec.x == "runLabel") {
        std::vector<std::string> datasets, variations;
        for (const RunData& run : runs) {
          const std::size_t split = run.runLabel.find('_');
          datasets.push_back(run.runLabel.substr(0, split));
          variations.push_back(split == std::string::npos ? "" : run.runLabel.substr(split + 1));
        }
        const bool sameDataset = std::set<std::string>(datasets.begin(), datasets.end()).size() == 1;
        const bool sameVariation = std::set<std::string>(variations.begin(), variations.end()).size() == 1;
        for (std::size_t run = 0; run < runs.size(); ++run)
          runs[run].xLabel = sameDataset ? variations[run] : sameVariation ? datasets[run] : runs[run].runLabel;
      }
      // Two runs with the same x would share a bin label: the run label tells them apart
      std::map<std::string, int> runsOfXLabel;
      for (const RunData& run : runs)
        ++runsOfXLabel[run.xLabel];
      for (RunData& run : runs)
        if (runsOfXLabel[run.xLabel] > 1)
          run.xLabel += " (" + run.runLabel + ")";
      std::vector<std::string> runLabels;
      for (const RunData& run : runs)
        runLabels.push_back(run.xLabel);
      const std::vector<std::string> otherRunLabels(runLabels.begin() + 1, runLabels.end());
      const std::string xTitle = spec.x == "runLabel" ? "run" : spec.x == "centralityLabel" ? "centrality" : spec.x;

      // A seed group shares its seed within one invocation of the coordinator only (README, "Seeds"): when some runs share a seed and
      // others do not, part of a group was probably rerun alone
      const bool allSmudged = std::all_of(runs.begin(), runs.end(), [](const RunData& run) { return run.setup.smudgeSigma > 0.; });
      std::map<std::string, std::vector<std::string>> runsOfSeed;
      std::string seedSummary;
      for (const RunData& run : runs) {
        const std::string seed = entryOf(run.producerConfig[kSampleBullet], "seed");
        runsOfSeed[seed].push_back(run.runLabel);
        seedSummary += run.runLabel + ": " + seed + "\n";
      }
      if (allSmudged && runsOfSeed.size() > 1 && runsOfSeed.size() < runs.size())
        warnings.push_back(prefix + "some runs share a seed and others do not, which a seed group partly rerun with --only leaves behind (seeds in Config/" + spec.name + "/seeds)");

      const OutputFolder folder = makeFolder(outputFile.get(), spec.name);
      const OutputFolder differenceFolder = makeFolder(folder.canvases, "DifferenceToFirst");
      TDirectory* comparisonConfig = configDir->mkdir(spec.name.c_str());
      std::vector<int> presentSamples;
      for (int sample = 0; sample < kNCurveSamples; ++sample)
        if (std::any_of(runs.begin(), runs.end(), [sample](const RunData& run) { return hasSample(run, sample); }))
          presentSamples.push_back(sample);
      const EstimateLayout& trueLayout = runs.front().results[kSetTrue].layout;

      // --------------------------------------------------------------------------------------------------------------------------------
      // 3b. Each run against x. Every point carries the jackknife error of its own run, and points that share initial conditions are
      //     correlated: DifferenceToFirst/ is where they are compared
      // --------------------------------------------------------------------------------------------------------------------------------
      // Fig. 11 analogue: the integrated R of the bullet sample by three estimators, in percent
      {
        std::vector<std::vector<EstimatePoint>> points(kFig11Scalars.size());
        for (std::size_t estimator = 0; estimator < kFig11Scalars.size(); ++estimator)
          for (const RunData& run : runs)
            points[estimator].push_back(pointOf(run, kSetTrue, kSampleBullet, trueLayout.scalar(kFig11Scalars[estimator])));
        drawRunTables(folder, "cFig11", "Fig. 11 analogue: integrated #it{R} of the bullet sample (%), jackknife errors per run", ";" + xTitle + ";#it{R} (%)", runLabels, std::vector<std::string>(kFig11Labels.begin(), kFig11Labels.end()), points, kPercent, kNoValue, warnings);
      }
      // Integrals of the true trigger and its lobe features, per sample
      const std::array<int, 3> integralBlocks = {0, kNRingEstimators, kNRingIntegrals}; // first scalar of R, R_z and the other projections
      const std::array<int, 3> integralBlockSizes = {kNRingEstimators, kNRingEstimators, kScalarNzSqD3N - kNRingIntegrals};
      const std::array<std::string, 3> integralBlockKeys = {"Ring", "RingZ", "Projections"};
      const std::array<std::string, 3> integralBlockTitles = {"#it{R} by its four estimators", "#it{R}_{z} by its four estimators", "the other projections"};
      for (const int sample : presentSamples) {
        const std::string key = kSampleKeys[sample];
        for (int block = 0; block < 3; ++block) {
          std::vector<std::string> estimateLabels;
          std::vector<std::vector<EstimatePoint>> points;
          for (int scalar = integralBlocks[block]; scalar < integralBlocks[block] + integralBlockSizes[block]; ++scalar) {
            estimateLabels.push_back(kScalarLatex[scalar]);
            points.emplace_back();
            for (const RunData& run : runs)
              points.back().push_back(pointOf(run, kSetTrue, sample, trueLayout.scalar(scalar)));
          }
          drawRunTables(folder, "c" + integralBlockKeys[block] + "_" + key, "Integrals of the true trigger, " + integralBlockTitles[block] + ", " + kSampleTitles[sample] + ", jackknife errors per run", ";" + xTitle + ";units of P", runLabels, estimateLabels, points, 1., 0., warnings);
        }
        // Lobes: the bullet sample and the jet-induced difference (the no_bullet lobes are noise with fluctuating initial conditions)
        if (sample == kSampleNoBullet)
          continue;
        for (int curve = 0; curve < kNFeatureCurves; ++curve) {
          for (const bool positions : {true, false}) {
            const int firstFeature = positions ? 0 : kNLobePositions;
            const int lastFeature = positions ? kNLobePositions : kNLobeFeatures;
            std::vector<std::string> estimateLabels;
            std::vector<std::vector<EstimatePoint>> points;
            for (int feature = firstFeature; feature < lastFeature; ++feature) {
              estimateLabels.push_back(kLobeLabels[feature]);
              points.emplace_back();
              for (const RunData& run : runs)
                points.back().push_back(pointOf(run, kSetTrue, sample, trueLayout.lobe(curve, feature)));
            }
            const std::string kind = positions ? "Positions" : "Heights";
            drawRunTables(folder, "cLobe" + kind + kFeatureCurveKeys[curve] + "_" + key, "Lobe " + std::string(positions ? "positions" : "heights") + " of " + kFeatureCurveSymbols[curve] + "(#Delta#varphi), " + kSampleTitles[sample] + ", jackknife errors per run", ";" + xTitle + ";" + (positions ? "#Delta#varphi (rad)" : kFeatureCurveSymbols[curve]), runLabels, estimateLabels, points, 1., 0., warnings);
          }
        }
      }
      // Curves of the true trigger, one band per run
      for (const int sample : presentSamples)
        for (int curve = 0; curve < kNFeatureCurves; ++curve)
          drawRunCurves(folder, runs, kSetTrue, curve, sample, "c" + kFeatureCurveKeys[curve] + "VsDeltaPhi_" + kSampleKeys[sample], "<" + kFeatureCurveSymbols[curve] + ">(#Delta#varphi) of the true trigger, " + kSampleTitles[sample] + ", jackknife bands per run", warnings);

      // Smudging: the smudged curves, the smudged/true integrals (jackknifed jointly within each run) and the leakage residual
      if (setInAllRuns(runs, kSetSmudge)) {
        const int smudgeBins = runs.front().results[kSetSmudge].layout.nBins;
        for (const int sample : presentSamples) {
          const std::string key = kSampleKeys[sample];
          for (int curve = 0; curve < kNFeatureCurves; ++curve)
            drawRunCurves(folder, runs, kSetSmudge, curve, sample, "cSmudged" + kFeatureCurveKeys[curve] + "VsDeltaPhi_" + key, "<" + kFeatureCurveSymbols[curve] + ">(#Delta#varphi) of the smudged triggers, " + kSampleTitles[sample] + ", jackknife bands per run", warnings);
          std::vector<std::vector<EstimatePoint>> ratioPoints(2);
          const std::array<int, 2> ratioScalars = {kScalarRingD3N, kScalarRingZD3N};
          for (int ratio = 0; ratio < 2; ++ratio) {
            for (const RunData& run : runs) {
              const int index = run.results[kSetTrue].layout.scalar(ratioScalars[ratio]);
              if (!hasSample(run, sample)) {
                ratioPoints[ratio].emplace_back();
              } else if (sample != kSampleDifference) {
                ratioPoints[ratio].push_back(crossSetRatio(run.results[kSetSmudge].stage1[sample].full, run.results[kSetSmudge].stage1[sample].replicas, run.results[kSetTrue].stage1[sample].full, run.results[kSetTrue].stage1[sample].replicas, index));
              } else if (run.paired) {
                ratioPoints[ratio].push_back(crossSetRatio(run.results[kSetSmudge].differenceStage1.full, run.results[kSetSmudge].differenceStage1.paired, run.results[kSetTrue].differenceStage1.full, run.results[kSetTrue].differenceStage1.paired, index));
              } else {
                warnings.push_back(prefix + run.runLabel + ": --paired=0, so no smudged/true ratio of the difference");
                ratioPoints[ratio].emplace_back();
              }
            }
          }
          drawRunTables(folder, "cSmudgedOverTrue_" + key, "Integrals, smudged over true trigger, " + kSampleTitles[sample] + ", jackknife errors per run (joint within the run)", ";" + xTitle + ";smudged/true", runLabels, {kScalarLatex[kScalarRingD3N], kScalarLatex[kScalarRingZD3N]}, ratioPoints, 1., 1., warnings);
        }
        // Leakage residual of the integral, R' - (R_true <cos zeta> + P_e <sin zeta>), per sample (README, "Ring overlap")
        const LeakageLayout leakageLayout{smudgeBins};
        std::vector<std::string> leakageLabels;
        std::vector<std::vector<EstimatePoint>> leakagePoints;
        for (const int sample : presentSamples) {
          leakageLabels.push_back(kSampleTitles[sample]);
          leakagePoints.emplace_back();
          for (const RunData& run : runs) {
            const EstimateSet& leakage = run.leakage[kSetSmudge][sample];
            leakagePoints.back().push_back(hasSample(run, sample) && !leakage.value.empty() ? pointAt(leakage, leakageLayout.residual(smudgeBins)) : EstimatePoint{});
          }
        }
        drawRunTables(folder, "cLeakageResidualSmudged", "Leakage ansatz of the smudged triggers: #it{R}' - prediction of the integral, jackknife errors per run (joint within the run)", ";" + xTitle + ";#it{R}' - prediction", runLabels, leakageLabels, leakagePoints, 1., 0., warnings);
      } else if (setInSomeRun(runs, kSetSmudge)) {
        warnings.push_back(prefix + "only some runs have smudged triggers, so there are no smudging plots");
      }

      // Expected significance of the bullet - no_bullet integrals (README, "Expected precision"), with the N and alpha of each run
      {
        std::vector<std::vector<EstimatePoint>> points(kPrecisionKeys.size());
        std::string nLambdaSummary;
        for (const RunData& run : runs) {
          nLambdaSummary += run.runLabel + ": " + std::to_string(run.nLambdaPrecision) + "\n";
          const std::array<double, 3> eSquared = integratedESquared(run.results[kSetTrue].estimates[kSampleBullet].value, run.results[kSetTrue].layout);
          for (std::size_t observable = 0; observable < kPrecisionKeys.size(); ++observable) {
            EstimatePoint significance;
            const double signal = run.hasNoBullet ? run.results[kSetTrue].estimates[kSampleDifference].value[run.results[kSetTrue].layout.scalar(kPrecisionDifferenceScalars[observable])] : kNoValue;
            if (signal != kNoValue && eSquared[observable] != kNoValue && run.nLambdaPrecision > 0.) {
              significance.value = signal / (sigmaSqrtNOf(eSquared[observable], run.alphaLambda) / std::sqrt(run.nLambdaPrecision));
              significance.error = 0.; // a prediction, no statistical error
            }
            points[observable].push_back(significance);
          }
        }
        drawRunTables(folder, "cSignificance", "Expected significance of the bullet - no_bullet integrals (N of each run in Config/" + spec.name + "/nLambdaPrecision)", ";" + xTitle + ";signal/#sigma", runLabels, std::vector<std::string>(kPrecisionSymbols.begin(), kPrecisionSymbols.end()), points, 1., 0., warnings);
        comparisonConfig->cd();
        TNamed("nLambdaPrecision", nLambdaSummary.c_str()).Write();
      }

      // --------------------------------------------------------------------------------------------------------------------------------
      // 3c. DifferenceToFirst/: every run minus the first, which is the reference only because it comes first in the compare line. Joint
      //     jackknife over the shared events when the two runs share initial conditions, independent errors otherwise. The runs that
      //     define a set identically must coincide exactly over the shared events: those differences are closures
      // --------------------------------------------------------------------------------------------------------------------------------
      std::string differenceSummary, closureSummary;
      for (std::size_t index = 0; index < kDifferenceSets.size(); ++index) {
        const int set = kDifferenceSets[index];
        if (!setInAllRuns(runs, set)) {
          if (setInSomeRun(runs, set))
            warnings.push_back(prefix + "only some runs have the " + kSetLabels[set] + ", so there are no differences of them");
          continue;
        }
        const EstimateLayout& layout = runs.front().results[set].layout;
        for (const int sample : presentSamples) {
          if (!hasSample(runs.front(), sample)) {
            warnings.push_back(prefix + "the first run, " + runs.front().runLabel + ", has no " + kSampleTitles[sample] + ", so there are no differences of the " + kSetLabels[set] + " for it");
            continue;
          }
          std::vector<RunDifference> differences;
          for (std::size_t run = 1; run < runs.size(); ++run) {
            const std::string label = prefix + runs[run].runLabel + " - " + runs.front().runLabel + ", " + kSetLabels[set] + ", " + kSampleTitles[sample];
            differences.push_back(differenceToFirst(runs.front(), runs[run], set, sample, label, warnings));
            const RunDifference& difference = differences.back();
            if (!difference.computed)
              continue;
            differenceSummary += kSetLabels[set] + ", " + kSampleTitles[sample] + ", " + runs[run].runLabel + ": " + (difference.joint ? "joint jackknife over " + std::to_string(difference.nSharedEvents) + " shared " + (sample == kSampleDifference ? "pairs" : "events") : "independent runs, errors in quadrature") + "\n";
            if (!difference.joint || !sameSetDefinition(runs.front(), runs[run], set, sample))
              continue;
            // Closure: same definition, same events, so the same per-event values and an exactly vanishing difference
            double maxDifference = 0., maxValue = 0.;
            for (int component = 0; component < layout.stage1Size(); ++component) {
              const double value = runs.front().results[set].estimates[sample].value[component];
              if (value != kNoValue)
                maxValue = std::max(maxValue, std::abs(value));
              if (difference.estimates.value[component] != kNoValue)
                maxDifference = std::max(maxDifference, std::abs(difference.estimates.value[component]));
            }
            const bool closes = maxDifference <= kSumCheckTolerance * maxValue;
            closureSummary += kSetLabels[set] + ", " + kSampleTitles[sample] + ", " + runs[run].runLabel + " vs " + runs.front().runLabel + ": defined identically, max |difference| = " + std::to_string(maxDifference) + (closes ? " (identical, as expected)" : " (DIFFERS)") + "\n";
            if (!closes)
              warnings.push_back(label + ": defined identically by the producer options, yet differ by up to " + std::to_string(maxDifference) + " over the shared events");
          }
          const bool allJoint = std::all_of(differences.begin(), differences.end(), [](const RunDifference& difference) { return !difference.computed || difference.joint; });
          const bool noneJoint = std::none_of(differences.begin(), differences.end(), [](const RunDifference& difference) { return difference.computed && difference.joint; });
          const std::string errorTitle = allJoint ? "joint jackknife over the shared events" : noneJoint ? "independent runs, errors in quadrature" : "joint or independent errors (Config/" + spec.name + "/differenceErrors)";
          for (const bool ringZ : {false, true}) {
            std::vector<std::string> estimateLabels;
            std::vector<std::vector<EstimatePoint>> points;
            for (int scalar = ringZ ? kNRingEstimators : 0; scalar < (ringZ ? kNRingIntegrals : kNRingEstimators); ++scalar) {
              estimateLabels.push_back(kScalarLatex[scalar]);
              points.emplace_back();
              for (const RunDifference& difference : differences)
                points.back().push_back(difference.computed ? pointAt(difference.estimates, layout.scalar(scalar)) : EstimatePoint{});
            }
            const std::string symbol = ringZ ? "#it{R}_{z}" : "#it{R}";
            drawRunTables(differenceFolder, "c" + std::string(ringZ ? "RingZ" : "Ring") + kDifferenceSetKeys[index] + "_" + kSampleKeys[sample], symbol + " of the " + kSetLabels[set] + ", " + kSampleTitles[sample] + ": each run minus " + runs.front().xLabel + " (the first of the line), " + errorTitle, ";" + xTitle + ";" + symbol + " - " + symbol + "_{first}", otherRunLabels, estimateLabels, points, 1., 0., warnings);
          }
        }
      }

      // --------------------------------------------------------------------------------------------------------------------------------
      // 3d. Seed spread: when the runs differ only by their seed, the spread of their smudged integrals across runs is the part of the error
      //     that comes from the smudge draws alone, since the events are the same. Mixing has no seed (partner rule), and the closures of
      //     DifferenceToFirst/ show that its results do not change
      // --------------------------------------------------------------------------------------------------------------------------------
      std::string seedSpreadSummary = "not applicable: the runs differ by more than their seed, or have no smudging";
      const bool onlySeedsDiffer = allSmudged && runsOfSeed.size() == runs.size() && std::all_of(runs.begin() + 1, runs.end(), [&runs](const RunData& run) { return differOnlyBySeed(runs.front(), run); });
      if (onlySeedsDiffer) {
        seedSpreadSummary.clear();
        const EstimateLayout& smudgeLayout = runs.front().results[kSetSmudge].layout;
        const std::vector<std::string> scalarLabels(kScalarLatex.begin(), kScalarLatex.begin() + kNRingIntegrals);
        std::vector<DrawEntry> entries;
        for (const int sample : presentSamples) {
          std::vector<double> ratios(kNRingIntegrals, kNoValue);
          for (int scalar = 0; scalar < kNRingIntegrals; ++scalar) {
            std::vector<double> values;
            double sumErrors = 0.;
            for (const RunData& run : runs) {
              const EstimatePoint point = pointOf(run, kSetSmudge, sample, smudgeLayout.scalar(scalar));
              if (point.value == kNoValue || point.error == kNoValue)
                break;
              values.push_back(point.value);
              sumErrors += point.error;
            }
            if (values.size() != runs.size() || !(sumErrors > 0.))
              continue;
            const double mean = std::accumulate(values.begin(), values.end(), 0.) / values.size();
            double sumSquares = 0.;
            for (const double value : values)
              sumSquares += (value - mean) * (value - mean);
            const double spread = std::sqrt(sumSquares / (values.size() - 1.));
            ratios[scalar] = spread / (sumErrors / values.size());
            seedSpreadSummary += kSampleTitles[sample] + ", " + kScalarText[scalar] + ": spread across seeds " + std::to_string(spread) + ", mean jackknife error " + std::to_string(sumErrors / values.size()) + ", ratio " + std::to_string(ratios[scalar]) + "\n";
          }
          const Part part = tableHistogram(folder.parts, "hSeedSpread_" + kSampleKeys[sample], "Smudged integrals: spread across seeds over the mean jackknife error, " + kSampleTitles[sample] + ";;spread/#sigma_{jackknife}", scalarLabels, ratios, std::vector<double>(kNRingIntegrals, 0.), blockIndices(0, kNRingIntegrals), 1., warnings);
          entries.push_back({part, kSampleTitles[sample], kSampleLooks[sample], kStyleMarker});
        }
        drawOverlay(folder.canvases, "cSeedSpread", "Smudged integrals over " + std::to_string(runs.size()) + " seeds: spread across seeds over the mean jackknife error (the share of the error from the smudge draws)", entries, 1., kRingIntegralSeparators);
        std::cout << seedSpreadSummary;
      }

      // --------------------------------------------------------------------------------------------------------------------------------
      // 3e. Config/<comparison>/
      // --------------------------------------------------------------------------------------------------------------------------------
      std::string runSummary, fileSummary, initialConditionsSummary;
      for (const RunData& run : runs) {
        runSummary += run.runLabel + " (" + run.xLabel + ")\n";
        fileSummary += run.perConfigFile + "\n";
        initialConditionsSummary += run.runLabel + ": " + run.initialConditions + "\n";
      }
      comparisonConfig->cd();
      TNamed("x", spec.x.c_str()).Write();
      TNamed("runs", runSummary.c_str()).Write();
      TNamed("files", fileSummary.c_str()).Write();
      TNamed("seeds", seedSummary.c_str()).Write();
      TNamed("initialConditions", initialConditionsSummary.c_str()).Write();
      TNamed("differenceErrors", differenceSummary.empty() ? "none" : differenceSummary.c_str()).Write();
      TNamed("closures", closureSummary.empty() ? "none: no two runs define a set identically" : closureSummary.c_str()).Write();
      TNamed("seedSpread", seedSpreadSummary.c_str()).Write();
      std::cout << differenceSummary << closureSummary;
    }

    // ==================================================================================================================================
    // 4. Config/ and warnings
    // ==================================================================================================================================
    std::string comparisonSummary;
    for (const ComparisonSpec& spec : comparisons) {
      comparisonSummary += spec.name + ", x = " + spec.x + ":";
      for (const std::string& file : spec.perConfigFiles)
        comparisonSummary += " " + file;
      comparisonSummary += "\n";
    }
    std::string warningSummary;
    for (const std::string& warning : warnings)
      warningSummary += warning + "\n";
    configDir->cd();
    TNamed("outputFile", outputFileName.c_str()).Write();
    TNamed("comparisons", comparisonSummary.c_str()).Write();
    TNamed("warnings", warningSummary.empty() ? "none" : warningSummary.c_str()).Write();
    outputFile->Close();

    std::cout << std::endl;
    for (const std::string& warning : warnings)
      std::cout << "[WARNING] " << warning << std::endl;
    std::cout << "Wrote " << outputFileName << std::endl;
  } catch (const std::exception& error) {
    std::cerr << "[ERROR] " << error.what() << std::endl;
    return 1;
  }
  return 0;
}