/*
* ============================================================
* compareConsumerOutputs.cxx
* The purpose of this code is to check if some renames, simplifications and condensing I've made in a
* large commit (close to the twelfth of October, 2026) didn't break my histograms. Better than checking by eye!
* Not sure if this will ever come in hand later, but it is the main review tool of that commit.
* - Some technical description:
* Object-by-object regression check between two consumer outputs: one produced before a code change,
* one after. Every histogram and profile is compared bin by bin, flows included. Plain histograms are
* compared on contents and errors, profiles on their raw accumulators (sum of w*y, sum of w*y^2, bin
* entries, sum of w^2). Titles, axis titles and bin labels are compared too, but reported apart: a title
* fix is not a content change.
*
* Usage:
*   root -l -b -q 'compareConsumerOutputs.cxx("before.root", "after.root")'
*   root -l -b -q 'compareConsumerOutputs.cxx("before.root", "after.root", "nameMap.txt", 1e-12)'
* A concrete example: ~/RingPol/LHC25ae_pass2/ITSandTPC/results_consumer $> root -l -b -q \ 
* '/home/users/cicerodm/PhD_codes/O2PhysicsTools/compareConsumerOutputs.cxx("ConsumerResults_JustLambda_useRingZ_analysisCuts.root", "ConsumerResults_JustLambda_useRingZ_analysisCuts.root")'
* Or using the actual map in full:
* ~/RingPol/LHC25ae_pass2/ITSandTPC/results_consumer $> root -l -b -q '/home/users/cicerodm/PhD_codes/O2PhysicsTools/compareConsumerOutputs.cxx("/home/users/cicerodm/RingPol/LHC25ae_pass2/ITSandTPC/results_consumer-BeforeMajorCondensingAndSimplification/ConsumerResults_JustLambda.root", "/home/users/cicerodm/RingPol/LHC25ae_pass2/ITSandTPC/results_consumer/ConsumerResults_JustLambda.root", "/home/users/cicerodm/PhD_codes/O2PhysicsTools/compareConsumerOutputs-nameMap_S3.txt")'
* Adding some tolerance to that:
* ~/RingPol/LHC25ae_pass2/ITSandTPC/results_consumer $> root -l -b -q '/home/users/cicerodm/PhD_codes/O2PhysicsTools/compareConsumerOutputs.cxx("/home/users/cicerodm/RingPol/LHC25ae_pass2/ITSandTPC/results_consumer-BeforeMajorCondensingAndSimplification/ConsumerResults_JustLambda.root", "/home/users/cicerodm/RingPol/LHC25ae_pass2/ITSandTPC/results_consumer/ConsumerResults_JustLambda.root", "/home/users/cicerodm/PhD_codes/O2PhysicsTools/compareConsumerOutputs-nameMap_S3.txt", 1e-14)' 
* nameMap (optional): one "oldPath newPath" pair per line, paths relative to the file root. A newPath
*   of "-" marks an object that is expected to disappear. Lines starting with "#" are comments.
*   Objects not listed are expected at the same path in both files.
*   One map serves every config: entries whose old path is absent (a family switched off, an object booked
*   only in some modes) are listed for information but do not decide the verdict. A typo cannot hide there,
*   because the misspelt object then shows up as MISSING (old-path typo) or as MISSING plus EXTRA (new-path typo).
* relTol (optional, default 0): relative tolerance per bin. 0 means bit-exact.
* (This nameMap is *the* list for the original comparison this was made for in 10/2026)
*
* Needs only ROOT. Smoke test: compare a file against itself, which must print IDENTICAL.
* ============================================================
*/

#include <TArrayD.h>
#include <TAxis.h>
#include <TClass.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// Long lists are cut at this many lines per category, with the total count printed after them:
// a wrong name map can otherwise flood the terminal with hundreds of "missing" lines.
constexpr size_t kMaxPrintedPerCategory = 100;

/// @brief Everything the comparison found, sorted into categories that mean different things.
struct RegressionReport {
  std::vector<std::string> contentDiffs;        // Same object, same binning, different numbers: the category that matters
  std::vector<std::string> structureDiffs;      // Different class, dimension or binning: the bins cannot even be compared
  std::vector<std::string> labelDiffs;          // Titles, axis titles or bin labels only
  std::vector<std::string> missing;             // Expected in "after" (same or mapped path), not found there
  std::vector<std::string> unexpectedSurvivors; // Mapped to "-" but still present in "after" at the same path
  std::vector<std::string> extra;               // In "after" only, with nothing in "before" mapping to them
  std::vector<std::string> unusedMapEntries;    // Name-map lines whose old path does not exist in "before" (informational)
  int nIdentical = 0;
  int nRemovedAsExpected = 0;
  double maxRelDiff = 0.;
  std::string maxRelDiffPath;
};

// ------------------------------------------------------------
// Input
// ------------------------------------------------------------

/// @brief Collects the path of every TH1-derived object below dir, recursively.
/// Only the newest cycle of each key is taken, as a plain Get() would.
void CollectPaths(TDirectory* dir, const std::string& prefix, std::set<std::string>& paths)
{
  if (!dir) return;
  std::set<std::string> seenNames; // GetListOfKeys() lists the newest cycle first
  TIter next(dir->GetListOfKeys());
  while (TKey* key = static_cast<TKey*>(next())) {
    const std::string name = key->GetName();
    if (!seenNames.insert(name).second) continue;
    TClass* cls = TClass::GetClass(key->GetClassName());
    if (!cls) continue;
    const std::string path = prefix.empty() ? name : prefix + "/" + name;
    if (cls->InheritsFrom(TDirectory::Class())) {
      CollectPaths(dir->GetDirectory(name.c_str()), path, paths);
    } else if (cls->InheritsFrom(TH1::Class())) { // TH1, TH2, TH3 and every profile
      paths.insert(path);
    }
  }
}

/// @brief Reads the optional "oldPath newPath" map. Sets ok = false on malformed or duplicated lines.
std::map<std::string, std::string> LoadNameMap(const char* mapPath, bool& ok)
{
  std::map<std::string, std::string> nameMap;
  ok = true;
  if (!mapPath || mapPath[0] == '\0') return nameMap;

  std::ifstream in(mapPath);
  if (!in) {
    std::cerr << " ERROR: cannot open name map " << mapPath << std::endl;
    ok = false;
    return nameMap;
  }
  std::string line;
  int lineNumber = 0;
  while (std::getline(in, line)) {
    ++lineNumber;
    std::istringstream fields(line);
    std::string oldPath, newPath, trailing;
    if (!(fields >> oldPath) || oldPath[0] == '#') continue; // Blank line or comment
    if (!(fields >> newPath) || (fields >> trailing)) {
      std::cerr << " ERROR: name map line " << lineNumber << " is not \"oldPath newPath\": " << line << std::endl;
      ok = false;
      continue;
    }
    if (!nameMap.emplace(oldPath, newPath).second) {
      std::cerr << " ERROR: name map line " << lineNumber << " repeats the old path " << oldPath << std::endl;
      ok = false;
    }
  }
  return nameMap;
}

// ------------------------------------------------------------
// Comparison helpers
// ------------------------------------------------------------

/// @brief Relative difference, symmetric in a and b. Exactly zero when bit-identical, two NaNs included.
/// NaN against a number returns NaN, which every "<= tolerance" test then rejects.
double RelDiff(double a, double b)
{
  if (a == b) return 0.;
  if (std::isnan(a) && std::isnan(b)) return 0.;
  return std::fabs(a - b) / std::max(std::fabs(a), std::fabs(b)); // The scale is > 0 here, since a != b
}

/// @brief True for TProfile, TProfile2D and TProfile3D.
bool IsProfile(const TH1* h)
{
  return h->InheritsFrom(TProfile::Class()) || h->InheritsFrom(TProfile2D::Class()) || h->InheritsFrom(TProfile3D::Class());
}

/// @brief Bin entries of a profile cell (sum of weights), whichever the profile dimension.
double ProfileEntries(TH1* h, int cell)
{
  if (auto* p1 = dynamic_cast<TProfile*>(h)) return p1->GetBinEntries(cell);
  if (auto* p2 = dynamic_cast<TProfile2D*>(h)) return p2->GetBinEntries(cell);
  if (auto* p3 = dynamic_cast<TProfile3D*>(h)) return p3->GetBinEntries(cell);
  return 0.;
}

/// @brief Sum of w^2 of a profile cell. Zero when the profile was never filled with weights (array not allocated).
double ProfileSumW2(TH1* h, int cell)
{
  const TArrayD* arr = nullptr;
  if (auto* p1 = dynamic_cast<TProfile*>(h)) arr = p1->GetBinSumw2();
  else if (auto* p2 = dynamic_cast<TProfile2D*>(h)) arr = p2->GetBinSumw2();
  else if (auto* p3 = dynamic_cast<TProfile3D*>(h)) arr = p3->GetBinSumw2();
  return (arr && cell < arr->GetSize()) ? arr->At(cell) : 0.;
}

/// @brief Empty when both axes have the same bins and bit-identical edges, else a short description.
std::string CompareAxis(const TAxis* a, const TAxis* b, const char* axisName)
{
  if (a->GetNbins() != b->GetNbins())
    return Form("%s axis: %d vs %d bins", axisName, a->GetNbins(), b->GetNbins());
  for (int i = 1; i <= a->GetNbins() + 1; ++i) { // Low edge of bin N+1 is the upper edge of the axis
    if (a->GetBinLowEdge(i) != b->GetBinLowEdge(i))
      return Form("%s axis: edge %d is %.17g vs %.17g", axisName, i, a->GetBinLowEdge(i), b->GetBinLowEdge(i));
  }
  return "";
}

/// @brief Empty when titles, axis titles and bin labels all agree, else the first difference found.
std::string CompareLabels(TH1* a, TH1* b)
{
  if (std::string(a->GetTitle()) != b->GetTitle())
    return Form("title \"%s\" vs \"%s\"", a->GetTitle(), b->GetTitle());
  const TAxis* axesA[3] = {a->GetXaxis(), a->GetYaxis(), a->GetZaxis()};
  const TAxis* axesB[3] = {b->GetXaxis(), b->GetYaxis(), b->GetZaxis()};
  const char* axisNames[3] = {"x", "y", "z"};
  for (int k = 0; k < a->GetDimension(); ++k) {
    if (std::string(axesA[k]->GetTitle()) != axesB[k]->GetTitle())
      return Form("%s axis title \"%s\" vs \"%s\"", axisNames[k], axesA[k]->GetTitle(), axesB[k]->GetTitle());
    for (int i = 1; i <= axesA[k]->GetNbins(); ++i) {
      if (std::string(axesA[k]->GetBinLabel(i)) != axesB[k]->GetBinLabel(i))
        return Form("%s axis bin %d label \"%s\" vs \"%s\"", axisNames[k], i, axesA[k]->GetBinLabel(i), axesB[k]->GetBinLabel(i));
    }
  }
  return "";
}

/// @brief Compares one object pair and files the outcome in the report.
void CompareObjects(TH1* a, TH1* b, const std::string& pathA, const std::string& pathB,
                    double relTol, RegressionReport& report)
{
  const std::string tag = (pathA == pathB) ? pathA : pathA + " --> " + pathB;

  // --- Structure: nothing below makes sense if the cells do not correspond ---
  if (std::string(a->ClassName()) != b->ClassName()) {
    report.structureDiffs.push_back(tag + ": class " + a->ClassName() + " vs " + b->ClassName());
    return;
  }
  const TAxis* axesA[3] = {a->GetXaxis(), a->GetYaxis(), a->GetZaxis()};
  const TAxis* axesB[3] = {b->GetXaxis(), b->GetYaxis(), b->GetZaxis()};
  const char* axisNames[3] = {"x", "y", "z"};
  for (int k = 0; k < a->GetDimension(); ++k) {
    const std::string axisDiff = CompareAxis(axesA[k], axesB[k], axisNames[k]);
    if (!axisDiff.empty()) {
      report.structureDiffs.push_back(tag + ": " + axisDiff);
      return;
    }
  }

  // --- Content: every cell, flows included ---
  // Profiles are compared on their raw accumulators rather than on mean and error, since equal means
  // can hide different entries, and a refactor that changes the entries is exactly what this must catch.
  const bool isProfile = IsProfile(a);
  const TArrayD* sumWyA = isProfile ? dynamic_cast<const TArrayD*>(a) : nullptr; // Profiles derive from TArrayD
  const TArrayD* sumWyB = isProfile ? dynamic_cast<const TArrayD*>(b) : nullptr;
  const TArrayD* sumWy2A = a->GetSumw2(); // For a profile, fSumw2 holds the sum of w*y^2
  const TArrayD* sumWy2B = b->GetSumw2();

  double objectMaxRel = 0.;
  std::string firstDiff;
  std::string maxRelCell; // Where the largest relative difference sits (usually not the first differing cell)
  // Per quantity: the largest absolute difference and the largest |value| of the object, so a difference can also be judged
  // against the scale of the whole object (a bin whose sum is close to zero makes the per-cell relative difference blow up):
  double maxAbsDiff[4] = {0., 0., 0., 0.};
  double maxAbsValue[4] = {0., 0., 0., 0.};
  const int nCells = a->GetNcells();
  for (int cell = 0; cell < nCells; ++cell) {
    double valuesA[4], valuesB[4];
    const char* quantityNames[4];
    int nQuantities = 0;
    if (isProfile) {
      valuesA[0] = sumWyA->At(cell);
      valuesB[0] = sumWyB->At(cell);
      quantityNames[0] = "sum of w*y";
      valuesA[1] = (cell < sumWy2A->GetSize()) ? sumWy2A->At(cell) : 0.;
      valuesB[1] = (cell < sumWy2B->GetSize()) ? sumWy2B->At(cell) : 0.;
      quantityNames[1] = "sum of w*y^2";
      valuesA[2] = ProfileEntries(a, cell);
      valuesB[2] = ProfileEntries(b, cell);
      quantityNames[2] = "bin entries";
      valuesA[3] = ProfileSumW2(a, cell);
      valuesB[3] = ProfileSumW2(b, cell);
      quantityNames[3] = "sum of w^2";
      nQuantities = 4;
    } else {
      valuesA[0] = a->GetBinContent(cell);
      valuesB[0] = b->GetBinContent(cell);
      quantityNames[0] = "content";
      valuesA[1] = a->GetBinError(cell);
      valuesB[1] = b->GetBinError(cell);
      quantityNames[1] = "error";
      nQuantities = 2;
    }

    for (int q = 0; q < nQuantities; ++q) {
      const double rel = RelDiff(valuesA[q], valuesB[q]);
      if (std::isnan(rel) || rel > objectMaxRel) {
        objectMaxRel = std::isnan(rel) ? INFINITY : rel;
        std::ostringstream msg;
        msg << std::setprecision(17) << "cell " << cell << ", " << quantityNames[q] << ": " << valuesA[q] << " vs " << valuesB[q];
        maxRelCell = msg.str();
      }
      maxAbsDiff[q] = std::max(maxAbsDiff[q], std::fabs(valuesA[q] - valuesB[q]));
      maxAbsValue[q] = std::max(maxAbsValue[q], std::max(std::fabs(valuesA[q]), std::fabs(valuesB[q])));
      if (firstDiff.empty() && !(rel <= relTol)) {
        std::ostringstream msg;
        msg << std::setprecision(17) << "cell " << cell << ", " << quantityNames[q] << ": " << valuesA[q] << " vs " << valuesB[q];
        firstDiff = msg.str();
      }
    }
  }

  // Total entries is a separate counter in ROOT, so it is checked on its own:
  if (firstDiff.empty() && !(RelDiff(a->GetEntries(), b->GetEntries()) <= relTol)) {
    std::ostringstream msg;
    msg << std::setprecision(17) << "GetEntries(): " << a->GetEntries() << " vs " << b->GetEntries();
    firstDiff = msg.str();
  }

  if (objectMaxRel > report.maxRelDiff) {
    report.maxRelDiff = objectMaxRel;
    report.maxRelDiffPath = tag;
  }

  const std::string labelDiff = CompareLabels(a, b);
  if (!labelDiff.empty()) report.labelDiffs.push_back(tag + ": " + labelDiff);

  // Largest difference against the scale of the object, over the compared quantities:
  double objectScaleRel = 0.;
  for (int q = 0; q < 4; ++q) {
    if (maxAbsValue[q] > 0.) objectScaleRel = std::max(objectScaleRel, maxAbsDiff[q] / maxAbsValue[q]);
  }

  if (!firstDiff.empty()) {
    std::ostringstream msg;
    msg << tag << ": first difference at " << firstDiff << "\n       max relative difference " << std::setprecision(3) << objectMaxRel;
    if (maxRelCell != firstDiff) msg << " (at " << maxRelCell << ")";
    msg << std::setprecision(3) << "\n       max difference relative to the object's largest value: " << objectScaleRel;
    report.contentDiffs.push_back(msg.str());
  } else if (labelDiff.empty()) {
    report.nIdentical++;
  }
}

// ------------------------------------------------------------
// Output
// ------------------------------------------------------------

/// @brief Prints one category, cut at kMaxPrintedPerCategory lines.
void PrintCategory(const std::string& heading, const std::vector<std::string>& lines)
{
  if (lines.empty()) return;
  std::cout << "\n " << heading << " (" << lines.size() << "):\n";
  for (size_t i = 0; i < lines.size() && i < kMaxPrintedPerCategory; ++i)
    std::cout << "   " << lines[i] << "\n";
  if (lines.size() > kMaxPrintedPerCategory)
    std::cout << "   ... and " << lines.size() - kMaxPrintedPerCategory << " more\n";
}

// ------------------------------------------------------------
// Main macro
// ------------------------------------------------------------
void compareConsumerOutputs(const char* beforePath, const char* afterPath, const char* nameMapPath = "", double relTol = 0.)
{
  std::cout << "\n=======================================================\n";
  std::cout << " Consumer output regression check\n";
  std::cout << "=======================================================\n";
  std::cout << " Before : " << beforePath << "\n";
  std::cout << " After  : " << afterPath << "\n";
  std::cout << " Map    : " << ((nameMapPath && nameMapPath[0]) ? nameMapPath : "(none)") << "\n";
  std::cout << " relTol : " << relTol << (relTol > 0. ? "" : " (bit-exact)") << "\n";

  bool mapOk = true;
  const std::map<std::string, std::string> nameMap = LoadNameMap(nameMapPath, mapOk);
  if (!mapOk) {
    std::cerr << " ERROR: fix the name map first, since a half-read map would report wrong differences." << std::endl;
    return;
  }

  TFile* fileA = TFile::Open(beforePath, "READ");
  TFile* fileB = TFile::Open(afterPath, "READ");
  if (!fileA || fileA->IsZombie() || !fileB || fileB->IsZombie()) {
    std::cerr << " ERROR: cannot open one of the input files." << std::endl;
    return;
  }

  std::set<std::string> pathsA, pathsB;
  CollectPaths(fileA, "", pathsA);
  CollectPaths(fileB, "", pathsB);
  std::cout << " Objects: " << pathsA.size() << " before, " << pathsB.size() << " after\n";

  RegressionReport report;
  std::set<std::string> matchedB; // "after" objects some "before" object was compared against

  for (const auto& entry : nameMap) {
    if (!pathsA.count(entry.first)) report.unusedMapEntries.push_back(entry.first);
  }

  for (const std::string& pathA : pathsA) {
    const auto mapped = nameMap.find(pathA);
    const std::string pathB = (mapped != nameMap.end()) ? mapped->second : pathA;

    if (pathB == "-") { // Expected removal
      if (pathsB.count(pathA)) {
        report.unexpectedSurvivors.push_back(pathA);
        matchedB.insert(pathA); // Reported once, here, and not again as "extra"
      } else {
        report.nRemovedAsExpected++;
      }
      continue;
    }
    if (!pathsB.count(pathB)) {
      report.missing.push_back((pathA == pathB) ? pathA : pathA + " --> " + pathB);
      continue;
    }
    matchedB.insert(pathB);

    TH1* hA = fileA->Get<TH1>(pathA.c_str());
    TH1* hB = fileB->Get<TH1>(pathB.c_str());
    if (!hA || !hB) { // Listed as TH1-derived by its key, so this should not happen
      report.structureDiffs.push_back(pathA + ": could not be read as TH1");
    } else {
      CompareObjects(hA, hB, pathA, pathB, relTol, report);
    }
    delete hA; // Read objects belong to us, and deleting them keeps memory flat over thousands of objects
    delete hB;
  }

  for (const std::string& pathB : pathsB) {
    if (!matchedB.count(pathB)) report.extra.push_back(pathB);
  }

  // --- Report ---
  PrintCategory("CONTENT DIFFERENCES", report.contentDiffs);
  PrintCategory("STRUCTURE DIFFERENCES (class or binning)", report.structureDiffs);
  PrintCategory("MISSING in after", report.missing);
  PrintCategory("EXPECTED REMOVALS still present", report.unexpectedSurvivors);
  PrintCategory("EXTRA in after (new objects)", report.extra);
  PrintCategory("NAME-MAP entries not found in before (informational: families off, mode-dependent bookings)", report.unusedMapEntries);
  PrintCategory("LABEL DIFFERENCES only (titles, axis titles, bin labels)", report.labelDiffs);

  std::cout << "\n-------------------------------------------------------\n";
  std::cout << " Identical objects    : " << report.nIdentical << "\n";
  std::cout << " Removed as expected  : " << report.nRemovedAsExpected << "\n";
  std::cout << " Max relative diff    : " << report.maxRelDiff
            << (report.maxRelDiffPath.empty() ? "" : "  (" + report.maxRelDiffPath + ")") << "\n";

  // Labels and unused map entries do not decide the verdict: a title fix is intended to show up there and nowhere else,
  // and a map entry for a family that is switched off has nothing to match.
  const bool identical = report.contentDiffs.empty() && report.structureDiffs.empty() && report.missing.empty() &&
                         report.unexpectedSurvivors.empty() && report.extra.empty();
  std::cout << " RESULT: " << (identical ? "IDENTICAL" : "DIFFERENCES FOUND")
            << (report.labelDiffs.empty() ? "" : Form(" (plus %zu label-only differences)", report.labelDiffs.size()))
            << "\n=======================================================\n\n";

  fileA->Close();
  fileB->Close();
  delete fileA;
  delete fileB;
}
