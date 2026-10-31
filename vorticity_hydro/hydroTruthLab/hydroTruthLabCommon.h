// hydroTruthLabCommon.h
//
// What hydroTruthLabPerConfig.cxx and hydroTruthLabSummaryPlots.cxx share: the layout of the producer files, the estimates with their
// jackknife, and the drawing of the output files. Not a macro of its own: each macro includes it from this folder, so the compilation of
// runHydroTruthLab.sh is unchanged.
//
// Every function is inline, so that the header could be included by two files of one program without breaking the one-definition rule.
//
// Estimators, errors and the output reference: README.md in this folder ("Errors", "PerConfig" and "SummaryPlots").

#ifndef HYDROTRUTHLABCOMMON_H
#define HYDROTRUTHLABCOMMON_H

#include <TAxis.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TNamed.h>
#include <TProfile.h>
#include <TTree.h>
#include <TVirtualPad.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// ======================================================================================================================================
// Constants
// ======================================================================================================================================

constexpr double kNoValue = -999.; // undefined estimate. Same guard value as the producer
// <n_z^2> of a bin below this is rounding noise, not geometry: for t_z = 0, n_z vanishes at Delta phi = 0 and pi, where the code gets
// sin(pi) = 1.2e-16 instead of 0. Half a cell away from those points <n_z^2> is already ~1e-3, so the cut separates the two cleanly
constexpr double kMinNzSqForAmplitude = 1e-12;
// |R| below this is rounding noise too: the delta-scan integrated over delta vanishes identically (half-turn identity), and R_z/R would
// divide by it
constexpr double kMinAbsRingForRatio = 1e-12;
constexpr double kPercent = 100.;            // Fig. 11 and the polarization maps only: the paper and the consumer's maps use percent
constexpr double kSumCheckTolerance = 1e-12; // relative. Sum of sumWDeltaPhi over bins vs sumWInCuts: the same numbers, summed in the same order
constexpr double kAxisTolerance = 1e-9;      // in units of the bin width, for the bins centred at Delta phi = 0 and delta = 0
constexpr double kPi = TMath::Pi();

// Drawing
constexpr double kLabelSize = 0.04;
constexpr double kTitleSize = 0.045;
constexpr double kYTitleOffset = 1.25;
constexpr double kLegendRowHeight = 0.05;   // NDC per legend entry
constexpr int kMaxHorizontalLabels = 4;     // tables with more bin labels than this get them vertical
constexpr double kTableBottomMargin = 0.32; // room for the vertical labels

// Trigger sets of a producer file, each read from the event tree as per-event, per-bin sums (README, "Event tree")
enum TriggerSet { kSetTrue = 0,
                  kSetScan,
                  kSetMixing,
                  kSetSmudge,
                  kNTriggerSets };
// Producer folder and per-bin profiles of each set: their pooled-master errors are compared with the jackknife in Jackknife/
const std::array<std::string, kNTriggerSets> kSetProducerFolder = {"Ring", "DeltaScan", "Mixing", "Smudge"};
const std::array<std::string, kNTriggerSets> kSetRingProfile = {"pRingVsDeltaPhi", "pRingD3NVsDelta", "pRingVsDeltaPhi", "pRingVsDeltaPhi"};
const std::array<std::string, kNTriggerSets> kSetRingZProfile = {"pRingZVsDeltaPhi", "pRingZD3NVsDelta", "pRingZVsDeltaPhi", "pRingZVsDeltaPhi"};
// Names used in the PerConfig file
const std::array<std::string, kNTriggerSets> kSetKeys = {"true", "deltaScan", "mixing", "smudge"};
const std::array<std::string, kNTriggerSets> kSetLabels = {"true trigger", "delta-scan", "mixed triggers", "smudged triggers"};
const std::array<std::string, kNTriggerSets> kSetAxisKeys = {"DeltaPhi", "Delta", "DeltaPhi", "DeltaPhi"};

// Per-bin sums of one set: the weights B_b and the numerators A_b of Ring, RingZ, n_z^2 and of the ring overlap cos(zeta) and sin(zeta)
enum BinArray { kBinW = 0,
                kBinARing,
                kBinARingZ,
                kBinANzSq,
                kBinACosZeta,
                kBinASinZeta,
                kNBinArrays };
// An empty branch name means that the set has no such sums: the true trigger has no ring overlap with itself
const std::array<std::string, kNBinArrays> kTrueBinBranches = {"sumWDeltaPhi", "sumARingDeltaPhi", "sumARingZDeltaPhi", "sumANzSqDeltaPhi", "", ""};
const std::array<std::string, kNBinArrays> kMixBinBranches = {"mixSumWDeltaPhi", "mixSumARingDeltaPhi", "mixSumARingZDeltaPhi", "mixSumANzSqDeltaPhi", "mixSumACosZetaDeltaPhi", "mixSumASinZetaDeltaPhi"};
const std::array<std::string, kNBinArrays> kSmudgeBinBranches = {"smudgeSumWDeltaPhi", "smudgeSumARingDeltaPhi", "smudgeSumARingZDeltaPhi", "smudgeSumANzSqDeltaPhi", "smudgeSumACosZetaDeltaPhi", "smudgeSumASinZetaDeltaPhi"};
// The scan stores integrals per delta, not sums: they are read back as A = value x sumWInCuts, with sumWInCuts as the weight (no kBinW branch)
const std::array<std::string, kNBinArrays> kScanBranches = {"", "scanRingD3N", "scanRingZD3N", "scanNzSqD3N", "scanCosZetaD3N", "scanSinZetaD3N"};

// Integrated-only observables of the true trigger (no per-bin sums in the tree): per-event integrals, weighted with sumWInCuts
enum WeightedExtra { kExtraRingZT = 0,
                     kExtraHelicityTerm,
                     kExtraRingZLab,
                     kExtraRingPerp,
                     kExtraNzSqSinSqTheta,
                     kExtraPolE,
                     kNWeightedExtras };
const std::array<std::string, kNWeightedExtras> kExtraBranches = {"ringZTD3N", "helicityTermD3N", "ringZLabD3N", "ringPerpD3N", "nzSqSinSqThetaD3N", "polED3N"};

// Per-event values whose plain event mean is an estimator of its own (the _EvtMean estimators, Fig. 11)
enum EventValue { kEvtRingD3N = 0,
                  kEvtRingFlat,
                  kEvtRingZD3N,
                  kEvtRingZFlat,
                  kNEventValues };
const std::array<std::string, kNEventValues> kTrueValueBranches = {"ringD3N", "ringFlat", "ringZD3N", "ringZFlat"};
const std::array<std::string, kNEventValues> kMixValueBranches = {"mixRingD3N", "mixRingFlat", "mixRingZD3N", "mixRingZFlat"};
const std::array<std::string, kNEventValues> kSmudgeValueBranches = {"smudgeRingD3N", "smudgeRingFlat", "smudgeRingZD3N", "smudgeRingZFlat"};

// Curves of an estimate vector, one value per bin. Ring and RingZ come first: they are the curves with lobe features, fold and asymmetry
enum CurveEstimate { kCurveRing = 0,
                     kCurveRingZ,
                     kCurveNzSq,
                     kCurveAmplitudeZ,
                     kCurveYieldFraction,
                     kCurveCosZeta, // variant sets only
                     kCurveSinZeta, // variant sets only
                     kNCurveEstimates };
constexpr int kNFeatureCurves = 2; // kCurveRing and kCurveRingZ
const std::array<std::string, kNFeatureCurves> kFeatureCurveKeys = {"Ring", "RingZ"};
const std::array<std::string, kNFeatureCurves> kFeatureCurveSymbols = {"#it{R}", "#it{R}_{z}"};

// Scalars of an estimate vector, in the order of the integral tables: R by four estimators, R_z by the same four, the other projections,
// all in units of the polarization, then the dimensionless ones. Everything is d3N-weighted and pooled unless labelled flat or evt. mean
enum ScalarEstimate { kScalarRingD3N = 0,
                      kScalarRingFlat, // flat average of the pooled curve
                      kScalarRingD3NEvtMean,
                      kScalarRingFlatEvtMean,
                      kScalarRingZD3N,
                      kScalarRingZFlat,
                      kScalarRingZD3NEvtMean,
                      kScalarRingZFlatEvtMean,
                      kScalarRingZTD3N,
                      kScalarHelicityTermD3N,
                      kScalarRingZLabD3N,
                      kScalarRingPerpD3N,
                      kScalarAmplitudeZD3N,
                      kScalarPolED3N,
                      kScalarNzSqD3N,
                      kScalarNzSqSinSqThetaD3N,
                      kScalarRingZOverRing,
                      kScalarCosZetaD3N, // variant sets only
                      kScalarSinZetaD3N, // variant sets only
                      kNScalarEstimates };
constexpr int kNPolarizationScalars = kScalarNzSqD3N;
constexpr int kNGeometryScalars = kScalarCosZetaD3N - kScalarNzSqD3N;
constexpr int kNRingIntegrals = kScalarRingZTD3N; // R and R_z by the four estimators, the integrals every set has
// Dashed separators of the integral tables (bin edges): after the four R estimators and after the four R_z estimators
const std::vector<double> kRingIntegralSeparators = {4., 8.};
// Labels: LaTeX for the plots, ASCII for the log
const std::array<std::string, kNScalarEstimates> kScalarLatex = {"#it{R}_{d^{3}N weight}", "#it{R}_{flat}", "#it{R}_{d^{3}N weight}, evt. mean", "#it{R}_{flat}, evt. mean", "#it{R}_{z, d^{3}N weight}", "#it{R}_{z, flat}", "#it{R}_{z, d^{3}N weight}, evt. mean", "#it{R}_{z, flat}, evt. mean", "#it{R}_{z}^{T}", "#it{H}", "#it{R}_{z}^{lab}", "#it{R}_{#perp}", "#it{A}_{z}", "#it{P}_{e}", "<n_{z}^{2}>", "<n_{z}^{2} sin^{2}#theta>", "#it{R}_{z}/#it{R}", "<cos#zeta>", "<sin#zeta>"};
const std::array<std::string, kNScalarEstimates> kScalarText = {"R_d3N", "R_flat", "R_d3N (evt. mean)", "R_flat (evt. mean)", "R_z_d3N", "R_z_flat", "R_z_d3N (evt. mean)", "R_z_flat (evt. mean)", "R_z^T_d3N", "H_d3N", "R_z^lab_d3N", "R_perp_d3N", "A_z_d3N", "P_e_d3N", "<n_z^2>_d3N", "<n_z^2 sin^2 theta>_d3N", "R_z/R (d3N)", "<cos zeta>_d3N", "<sin zeta>_d3N"};
const std::array<int, kNWeightedExtras> kExtraScalar = {kScalarRingZTD3N, kScalarHelicityTermD3N, kScalarRingZLabD3N, kScalarRingPerpD3N, kScalarNzSqSinSqThetaD3N, kScalarPolED3N};
const std::array<int, kNEventValues> kEventValueScalar = {kScalarRingD3NEvtMean, kScalarRingFlatEvtMean, kScalarRingZD3NEvtMean, kScalarRingZFlatEvtMean};

// Lobe features of a curve (README, "Errors"), in the order of the tables: the positions (rad), then the amplitudes (units of the curve)
enum LobeFeature { kLobeCentroidNeg = 0,
                   kLobeSeparation,
                   kLobeCentroidPos,
                   kLobeHeightNeg,
                   kLobeDip,
                   kLobeHeightPos,
                   kNLobeFeatures };
constexpr int kNLobePositions = kLobeHeightNeg;
const std::array<std::string, kNLobeFeatures> kLobeLabels = {"centroid (#Delta#varphi < 0)", "separation", "centroid (#Delta#varphi > 0)", "height (#Delta#varphi < 0)", "dip (#Delta#varphi = 0)", "height (#Delta#varphi > 0)"};
enum Lobe { kLobeNegative = 0,
            kLobePositive,
            kNLobes };
const std::array<int, kNLobes> kLobeCentroid = {kLobeCentroidNeg, kLobeCentroidPos};
const std::array<int, kNLobes> kLobeHeight = {kLobeHeightNeg, kLobeHeightPos};
// R_z against R, from the lobe features of both curves: the shifts R_z - R of the positions (same order as the positions), then the height
// ratios R_z/R of the two lobes
enum LobeComparison { kShiftCentroidNeg = 0,
                      kShiftSeparation,
                      kShiftCentroidPos,
                      kHeightRatioNeg,
                      kHeightRatioPos,
                      kNLobeComparisons };
constexpr int kNLobeShifts = kHeightRatioNeg;
const std::array<std::string, kNLobeComparisons> kLobeComparisonLabels = {"centroid (#Delta#varphi < 0)", "separation", "centroid (#Delta#varphi > 0)", "height (#Delta#varphi < 0)", "height (#Delta#varphi > 0)"};

// Curves of every set: each sample over all its events, and the bullet - no_bullet difference
enum CurveSample { kSampleBullet = 0,
                   kSampleNoBullet,
                   kSampleDifference,
                   kNCurveSamples };
constexpr int kNSamples = 2; // bullet and no_bullet, the first two CurveSample entries (one producer file each)
const std::array<std::string, kNCurveSamples> kSampleKeys = {"bullet", "noBullet", "difference"};
const std::array<std::string, kNSamples> kSampleLabels = {"bullet", "no_bullet"}; // in the logs

// Fig. 11 point: the integrated R of the bullet sample by three estimators
const std::array<int, 3> kFig11Scalars = {kScalarRingFlatEvtMean, kScalarRingD3NEvtMean, kScalarRingD3N};
const std::array<std::string, 3> kFig11Keys = {"FlatEvtMean", "D3NEvtMean", "D3N"};
const std::array<std::string, 3> kFig11Labels = {"flat, evt. mean #pm SEM (Vitor's estimator)", "d^{3}N weight, evt. mean #pm SEM", "d^{3}N weight, pooled #pm jackknife"};

// Expected precision (README, "Expected precision"): the observables, and the scalar of the bullet - no_bullet difference that is the signal
// of each
const std::array<std::string, 3> kPrecisionKeys = {"Ring", "RingZ", "RingZT"};
const std::array<std::string, 3> kPrecisionSymbols = {"#it{R}", "#it{R}_{z}", "#it{R}_{z}^{T}"};
const std::array<int, 3> kPrecisionDifferenceScalars = {kScalarRingD3N, kScalarRingZD3N, kScalarRingZTD3N};

// ======================================================================================================================================
// Drawing styles
// ======================================================================================================================================

/// \brief Colours and marker of one drawn entry. The bands are opaque light shades of the line colour: JSROOT (the vscode viewer) does not
///        render the E3 bands of a TH1, and the transparent colours of SetFillColorAlpha are session-dependent indices, wrong when a canvas is
///        reopened, so every band is a TGraphErrors drawn "3" with one of ROOT's predefined light colours.
struct Look {
  int color = kBlack;    // line, markers and error bars
  int bandColor = kGray; // band fill
  int marker = kFullCircle;
};
const std::array<Look, kNCurveSamples> kSampleLooks = {{{kRed + 1, kRed - 9, kFullCircle}, {kBlue + 1, kBlue - 9, kFullSquare}, {kBlack, kGray, kOpenCircle}}};
// Overlays of several observables of one sample, or of the runs of a comparison (SummaryPlots: one look per run, so at most eight runs)
const std::array<Look, 8> kPaletteLooks = {{{kBlack, kGray, kFullCircle}, {kRed + 1, kRed - 9, kFullSquare}, {kBlue + 1, kBlue - 9, kFullTriangleUp}, {kGreen + 2, kGreen - 9, kOpenCircle}, {kMagenta + 1, kMagenta - 9, kOpenSquare}, {kOrange + 7, kOrange - 9, kOpenTriangleUp}, {kCyan + 2, kCyan - 9, kFullDiamond}, {kViolet + 1, kViolet - 9, kFullCross}}};

enum DrawStyle { kStyleBand = 0, // error band and central line: every curve with errors
                 kStyleBandOnly, // error band without its line: the wider band of nested bands
                 kStyleLine,     // central line only: curves without errors
                 kStyleMarker,   // markers with error bars: tables
                 kStyleHistogram // the histogram itself: distributions
};

// ======================================================================================================================================
// Options
// ======================================================================================================================================

// The std::sto* functions throw with messages like "stoi", which say nothing about the option: these wrappers name it (as the producer does)
inline int toInt(const std::string& key, const std::string& value)
{
  std::size_t parsedLength = 0;
  int result = 0;
  try {
    result = std::stoi(value, &parsedLength);
  } catch (const std::exception&) {
    parsedLength = 0;
  }
  if (parsedLength == 0 || parsedLength != value.size())
    throw std::runtime_error("--" + key + " expects an integer, got '" + value + "'");
  return result;
}

inline double toDouble(const std::string& key, const std::string& value)
{
  std::size_t parsedLength = 0;
  double result = 0.;
  try {
    result = std::stod(value, &parsedLength);
  } catch (const std::exception&) {
    parsedLength = 0;
  }
  if (parsedLength == 0 || parsedLength != value.size())
    throw std::runtime_error("--" + key + " expects a number, got '" + value + "'");
  return result;
}

// ======================================================================================================================================
// Producer files
// ======================================================================================================================================

/// \brief Object at a path of a producer file, with a hard stop if it is missing.
template <typename T>
T* getObject(TFile* file, const std::string& path)
{
  T* object = file->Get<T>(path.c_str());
  if (object == nullptr)
    throw std::runtime_error(std::string(file->GetName()) + " has no " + path + " (not a hydroTruthLab producer file, or an older one?)");
  return object;
}

/// \brief Value of one option of a producer (stored in Config/ as name = option, title = value).
inline std::string readConfigValue(TFile* file, const std::string& name) { return getObject<TNamed>(file, "Config/" + name)->GetTitle(); }

/// \brief SetBranchAddress with a hard stop on a missing branch (a producer file older than this macro, for instance).
inline void attachBranch(TTree* tree, const std::string& name, void* address)
{
  if (tree->GetBranch(name.c_str()) == nullptr)
    throw std::runtime_error("the event tree has no branch " + name + ": rerun the producer");
  tree->SetBranchAddress(name.c_str(), address);
}

/// \brief Per-event sums of one trigger set of one sample, event-major: perBin[array][iEvent * nBins + bin]. An array the set does not have
///        (the ring overlap of the true trigger) stays empty.
struct SetSums {
  bool present = false;
  int nBins = 0;
  int nEvents = 0;
  std::array<std::vector<double>, kNBinArrays> perBin;
  bool hasWeightedExtras = false;                                  // true trigger only
  std::array<std::vector<double>, kNWeightedExtras> weightedExtra; // per event, A = integral x sumWInCuts
  bool hasEventValues = false;                                     // every set but the scan
  std::array<std::vector<double>, kNEventValues> eventValue;       // per event
};

/// \brief Everything PerConfig needs from one event tree.
struct SampleData {
  std::vector<int> eventId;
  std::vector<double> phiTrigger;
  std::vector<double> sumWInCuts;
  std::vector<int> mixDeltaEffBin; // event-major, mixDepth per event: delta_eff bin of each mixing partner, -1 outside the axis
  std::array<SetSums, kNTriggerSets> sets;
  // Delta-scan at delta = 0 against the true trigger, event by event: max |difference| and max |value|, indexed like BinArray (kBinW unused)
  std::array<double, kNBinArrays> scanZeroMaxResidual{};
  std::array<double, kNBinArrays> scanZeroMaxValue{};
};

/// \brief Reads the event tree of one producer file into the per-event sums of every trigger set, with entry-by-entry consistency checks.
inline SampleData readSample(TFile* file, const std::string& label, int nDeltaPhiBins, int nDelta, int scanZeroBin, int mixDepth, bool hasSmudge)
{
  const bool hasMixing = mixDepth > 0;
  TTree* tree = getObject<TTree>(file, "events");
  const Long64_t nEntries = tree->GetEntries();
  if (nEntries < 2)
    throw std::runtime_error(label + ": the jackknife needs at least two events, the tree has " + std::to_string(nEntries));

  // Array buffers sized by the largest count anywhere in the tree, so that a malformed entry is caught by the check below instead of
  // overflowing them
  const int trueSize = std::max<int>(nDeltaPhiBins, tree->GetMaximum("nDeltaPhiBins"));
  const int scanSize = std::max<int>({1, nDelta, static_cast<int>(tree->GetMaximum("nDelta"))});
  const int mixSize = std::max<int>(nDeltaPhiBins, tree->GetMaximum("nMixDeltaPhiBins"));
  const int smudgeSize = std::max<int>(nDeltaPhiBins, tree->GetMaximum("nSmudgeDeltaPhiBins"));

  Int_t eventId = 0;
  Double_t phiTrigger = 0.;
  Double_t sumWInCuts = 0.;
  Double_t nzSqD3N = 0.;
  Int_t nTrueBins = 0, nScanBins = 0, nMixBins = 0, nSmudgeBins = 0, nMixPartners = 0;
  std::array<Double_t, kNWeightedExtras> extraD3N{};
  std::array<Double_t, kNEventValues> trueValues{}, mixValues{}, smudgeValues{};
  std::array<std::vector<Double_t>, kNBinArrays> trueBins, scanValues, mixBins, smudgeBins;
  std::vector<Int_t> mixDeltaEffBin;

  // Event and true trigger
  attachBranch(tree, "eventId", &eventId);
  attachBranch(tree, "phiTrigger", &phiTrigger);
  attachBranch(tree, "sumWInCuts", &sumWInCuts);
  attachBranch(tree, "nzSqD3N", &nzSqD3N);
  attachBranch(tree, "nDeltaPhiBins", &nTrueBins);
  for (int extra = 0; extra < kNWeightedExtras; ++extra)
    attachBranch(tree, kExtraBranches[extra], &extraD3N[extra]);
  for (int value = 0; value < kNEventValues; ++value)
    attachBranch(tree, kTrueValueBranches[value], &trueValues[value]);
  for (int array = 0; array < kNBinArrays; ++array) {
    if (kTrueBinBranches[array].empty())
      continue;
    trueBins[array].assign(trueSize, 0.);
    attachBranch(tree, kTrueBinBranches[array], trueBins[array].data());
  }
  // Delta-scan. Its weight is sumWInCuts, attached above: kBinW has no array of its own
  attachBranch(tree, "nDelta", &nScanBins);
  if (nDelta > 0) {
    for (int array = kBinARing; array < kNBinArrays; ++array) {
      scanValues[array].assign(scanSize, 0.);
      attachBranch(tree, kScanBranches[array], scanValues[array].data());
    }
  }
  // Mixing
  attachBranch(tree, "nMixDeltaPhiBins", &nMixBins);
  attachBranch(tree, "nMixPartners", &nMixPartners);
  if (hasMixing) {
    mixDeltaEffBin.assign(std::max<int>(mixDepth, tree->GetMaximum("nMixPartners")), -1);
    attachBranch(tree, "mixDeltaEffBin", mixDeltaEffBin.data());
    for (int value = 0; value < kNEventValues; ++value)
      attachBranch(tree, kMixValueBranches[value], &mixValues[value]);
    for (int array = 0; array < kNBinArrays; ++array) {
      mixBins[array].assign(mixSize, 0.);
      attachBranch(tree, kMixBinBranches[array], mixBins[array].data());
    }
  }
  // Smudging
  attachBranch(tree, "nSmudgeDeltaPhiBins", &nSmudgeBins);
  if (hasSmudge) {
    for (int value = 0; value < kNEventValues; ++value)
      attachBranch(tree, kSmudgeValueBranches[value], &smudgeValues[value]);
    for (int array = 0; array < kNBinArrays; ++array) {
      smudgeBins[array].assign(smudgeSize, 0.);
      attachBranch(tree, kSmudgeBinBranches[array], smudgeBins[array].data());
    }
  }

  SampleData sample;
  SetSums& trueSet = sample.sets[kSetTrue];
  SetSums& scanSet = sample.sets[kSetScan];
  SetSums& mixSet = sample.sets[kSetMixing];
  SetSums& smudgeSet = sample.sets[kSetSmudge];
  trueSet.present = true;
  trueSet.nBins = nDeltaPhiBins;
  trueSet.hasWeightedExtras = true;
  trueSet.hasEventValues = true;
  scanSet.present = nDelta > 0;
  scanSet.nBins = nDelta;
  mixSet.present = hasMixing;
  mixSet.nBins = nDeltaPhiBins;
  mixSet.hasEventValues = true;
  smudgeSet.present = hasSmudge;
  smudgeSet.nBins = nDeltaPhiBins;
  smudgeSet.hasEventValues = true;

  for (Long64_t entry = 0; entry < nEntries; ++entry) {
    tree->GetEntry(entry);
    if (nTrueBins != nDeltaPhiBins || nScanBins != nDelta || (hasMixing && (nMixBins != nDeltaPhiBins || nMixPartners != mixDepth)) || (hasSmudge && nSmudgeBins != nDeltaPhiBins))
      throw std::runtime_error(label + ", event " + std::to_string(eventId) + ": the array sizes in the tree disagree with the producer's configuration");
    sample.eventId.push_back(eventId);
    sample.phiTrigger.push_back(phiTrigger);
    sample.sumWInCuts.push_back(sumWInCuts);

    // True trigger. Its integrated-only observables are weighted with sumWInCuts, which must equal the sum of the per-bin weights
    double sumWOverBins = 0.;
    for (int bin = 0; bin < nDeltaPhiBins; ++bin)
      sumWOverBins += trueBins[kBinW][bin];
    if (std::abs(sumWOverBins - sumWInCuts) > kSumCheckTolerance * std::abs(sumWInCuts))
      throw std::runtime_error(label + ", event " + std::to_string(eventId) + ": sumWInCuts differs from the sum of sumWDeltaPhi");
    for (int array = 0; array < kNBinArrays; ++array)
      if (!kTrueBinBranches[array].empty())
        trueSet.perBin[array].insert(trueSet.perBin[array].end(), trueBins[array].begin(), trueBins[array].begin() + nDeltaPhiBins);
    for (int extra = 0; extra < kNWeightedExtras; ++extra)
      trueSet.weightedExtra[extra].push_back(extraD3N[extra] * sumWInCuts);
    for (int value = 0; value < kNEventValues; ++value)
      trueSet.eventValue[value].push_back(trueValues[value]);

    // Delta-scan: integrals per delta, turned back into sums with the event weight. Exact for transverse triggers, which see every in-cut cell
    if (nDelta > 0) {
      for (int step = 0; step < nDelta; ++step) {
        scanSet.perBin[kBinW].push_back(sumWInCuts);
        for (int array = kBinARing; array < kNBinArrays; ++array)
          scanSet.perBin[array].push_back(scanValues[array][step] * sumWInCuts);
      }
      // delta = 0 is the true trigger itself (same cells, same trigger, summed in another order): the residual is float noise. Its ring
      // direction is the true one, so cos(zeta) = 1 and sin(zeta) = 0 there
      const std::array<double, kNBinArrays> trueIntegrals = {0., trueValues[kEvtRingD3N], trueValues[kEvtRingZD3N], nzSqD3N, 1., 0.};
      for (int array = kBinARing; array < kNBinArrays; ++array) {
        sample.scanZeroMaxResidual[array] = std::max(sample.scanZeroMaxResidual[array], std::abs(scanValues[array][scanZeroBin] - trueIntegrals[array]));
        sample.scanZeroMaxValue[array] = std::max(sample.scanZeroMaxValue[array], std::abs(trueIntegrals[array]));
      }
    }

    // Mixing and smudging: per-bin sums pooled over the partners (draws) of the event, and its per-event values
    if (hasMixing) {
      sample.mixDeltaEffBin.insert(sample.mixDeltaEffBin.end(), mixDeltaEffBin.begin(), mixDeltaEffBin.begin() + mixDepth);
      for (int array = 0; array < kNBinArrays; ++array)
        mixSet.perBin[array].insert(mixSet.perBin[array].end(), mixBins[array].begin(), mixBins[array].begin() + nDeltaPhiBins);
      for (int value = 0; value < kNEventValues; ++value)
        mixSet.eventValue[value].push_back(mixValues[value]);
    }
    if (hasSmudge) {
      for (int array = 0; array < kNBinArrays; ++array)
        smudgeSet.perBin[array].insert(smudgeSet.perBin[array].end(), smudgeBins[array].begin(), smudgeBins[array].begin() + nDeltaPhiBins);
      for (int value = 0; value < kNEventValues; ++value)
        smudgeSet.eventValue[value].push_back(smudgeValues[value]);
    }
  }
  for (SetSums& sums : sample.sets)
    if (sums.present)
      sums.nEvents = nEntries;
  return sample;
}

/// \brief Copy of a set restricted to some events (rows), in the given order: the pairs of the difference.
inline SetSums selectEvents(const SetSums& sums, const std::vector<int>& rows)
{
  SetSums selected;
  selected.present = sums.present;
  selected.nBins = sums.nBins;
  selected.nEvents = rows.size();
  selected.hasWeightedExtras = sums.hasWeightedExtras;
  selected.hasEventValues = sums.hasEventValues;
  for (const int row : rows) {
    for (int array = 0; array < kNBinArrays; ++array) {
      if (sums.perBin[array].empty())
        continue;
      const auto rowBegin = sums.perBin[array].begin() + row * sums.nBins;
      selected.perBin[array].insert(selected.perBin[array].end(), rowBegin, rowBegin + sums.nBins);
    }
    if (sums.hasWeightedExtras)
      for (int extra = 0; extra < kNWeightedExtras; ++extra)
        selected.weightedExtra[extra].push_back(sums.weightedExtra[extra][row]);
    if (sums.hasEventValues)
      for (int value = 0; value < kNEventValues; ++value)
        selected.eventValue[value].push_back(sums.eventValue[value][row]);
  }
  return selected;
}

struct AxisRange {
  int nBins = 0;
  double low = 0.;
  double high = 0.;
  std::string title;
};

/// \brief What every macro downstream of the producer needs from its configuration and axes, with the checks the estimates rely on.
struct ProducerSetup {
  std::string triggerMode;
  int nDelta = 0;
  int mixDepth = 0;
  double smudgeSigma = 0.;
  AxisRange deltaPhiRange;
  std::vector<double> deltaPhiCentres;
  AxisRange deltaRange; // empty without the delta-scan
  int scanZeroBin = -1; // bin of delta = 0, the true trigger itself
};

inline ProducerSetup readProducerSetup(TFile* file)
{
  ProducerSetup setup;
  setup.triggerMode = readConfigValue(file, "trigger");
  setup.nDelta = std::stoi(readConfigValue(file, "deltaScanSteps"));
  setup.mixDepth = std::stoi(readConfigValue(file, "mixDepth"));
  setup.smudgeSigma = std::stod(readConfigValue(file, "smudgeSigma"));
  // Delta phi axis, taken from the producer. The lobes, the fold and the asymmetry pair bin middle + j with bin middle - j, so the axis must
  // span 2 pi with an even number of bins and a bin centred at Delta phi = 0 (the producer's grid binning does)
  const TAxis* deltaPhiAxis = getObject<TProfile>(file, "Ring/pRingVsDeltaPhi")->GetXaxis();
  setup.deltaPhiRange = {deltaPhiAxis->GetNbins(), deltaPhiAxis->GetXmin(), deltaPhiAxis->GetXmax(), deltaPhiAxis->GetTitle()};
  const int nDeltaPhiBins = setup.deltaPhiRange.nBins;
  const double deltaPhiWidth = (setup.deltaPhiRange.high - setup.deltaPhiRange.low) / nDeltaPhiBins;
  setup.deltaPhiCentres.resize(nDeltaPhiBins);
  for (int bin = 0; bin < nDeltaPhiBins; ++bin)
    setup.deltaPhiCentres[bin] = deltaPhiAxis->GetBinCenter(bin + 1);
  if (nDeltaPhiBins % 2 != 0 || std::abs(setup.deltaPhiCentres[nDeltaPhiBins / 2]) > kAxisTolerance * deltaPhiWidth || std::abs(setup.deltaPhiRange.high - setup.deltaPhiRange.low - 2. * kPi) > kAxisTolerance * deltaPhiWidth)
    throw std::runtime_error(std::string(file->GetName()) + ": the Delta phi axis does not span 2 pi with an even number of bins centred on Delta phi = 0");
  // Delta axis of the scan (and of mixing vs delta_eff), with the bin of delta = 0, the true trigger itself
  if (setup.nDelta > 0) {
    const TAxis* deltaAxis = getObject<TProfile>(file, "DeltaScan/pRingD3NVsDelta")->GetXaxis();
    setup.deltaRange = {deltaAxis->GetNbins(), deltaAxis->GetXmin(), deltaAxis->GetXmax(), deltaAxis->GetTitle()};
    setup.scanZeroBin = deltaAxis->FindFixBin(0.) - 1;
    if (setup.deltaRange.nBins != setup.nDelta || std::abs(deltaAxis->GetBinCenter(setup.scanZeroBin + 1)) > kAxisTolerance * deltaAxis->GetBinWidth(1))
      throw std::runtime_error(std::string(file->GetName()) + ": the delta axis does not have deltaScanSteps bins with one centred on delta = 0");
  }
  return setup;
}

/// \brief Row of each event in the tree of a sample, with a hard stop on an event that appears twice.
inline std::map<int, int> rowsByEvent(const SampleData& sample, const std::string& label)
{
  std::map<int, int> rows;
  for (std::size_t row = 0; row < sample.eventId.size(); ++row)
    if (!rows.emplace(sample.eventId[row], row).second)
      throw std::runtime_error(label + ": event " + std::to_string(sample.eventId[row]) + " appears twice in the tree");
  return rows;
}

/// \brief Pairs of a run: the rows of the events both samples have, in the order of the bullet tree, with the trigger of each pair checked
///        (README, key hypothesis).
/// \param label  prefix of the error messages (the run, in SummaryPlots)
inline std::array<std::vector<int>, kNSamples> pairRowsOf(const SampleData& bullet, const SampleData& noBullet, const std::string& label)
{
  std::array<std::vector<int>, kNSamples> pairRows;
  rowsByEvent(bullet, label + "bullet"); // for its duplicate check only: the pairs follow the bullet tree
  const std::map<int, int> noBulletRows = rowsByEvent(noBullet, label + "no_bullet");
  for (std::size_t row = 0; row < bullet.eventId.size(); ++row) {
    const int eventId = bullet.eventId[row];
    const auto partner = noBulletRows.find(eventId);
    if (partner == noBulletRows.end())
      continue;
    // Both producers read phi_t from the same text with the same parser (or set it to 0), so the triggers of a pair are bit-identical
    if (bullet.phiTrigger[row] != noBullet.phiTrigger[partner->second])
      throw std::runtime_error(label + "event " + std::to_string(eventId) + ": phiTrigger differs between bullet and no_bullet, so the samples are not paired");
    pairRows[kSampleBullet].push_back(row);
    pairRows[kSampleNoBullet].push_back(partner->second);
  }
  return pairRows;
}

// ======================================================================================================================================
// Estimates and their jackknife (README, "Errors")
// ======================================================================================================================================

/// \brief Positions of the blocks of an estimate vector. Stage 1 (curves and scalars) exists for every set, the curve features (lobes of Ring
///        and RingZ, their comparison, fold and asymmetry) for the true trigger only.
struct EstimateLayout {
  int nBins = 0;
  int curve(int curveIndex, int bin) const { return curveIndex * nBins + bin; }
  int scalar(int scalarIndex) const { return kNCurveEstimates * nBins + scalarIndex; }
  int stage1Size() const { return kNCurveEstimates * nBins + kNScalarEstimates; }
  int middleBin() const { return nBins / 2; }          // centred at Delta phi = 0 (checked in main)
  int nFoldBins() const { return nBins / 2 + 1; }      // |Delta phi| = 0, w, ..., pi
  int nAsymmetryBins() const { return nBins / 2 - 1; } // |Delta phi| = w, ..., pi - w
  int lobe(int curveIndex, int feature) const { return stage1Size() + curveIndex * kNLobeFeatures + feature; }
  int lobeComparison(int comparison) const { return stage1Size() + kNFeatureCurves * kNLobeFeatures + comparison; }
  int fold(int curveIndex, int bin) const { return lobeComparison(kNLobeComparisons) + curveIndex * nFoldBins() + bin; }
  int asymmetry(int curveIndex, int bin) const { return fold(kNFeatureCurves, 0) + curveIndex * nAsymmetryBins() + bin; }
  int fullSize() const { return asymmetry(kNFeatureCurves, 0); }
};

/// \brief Positions in a leakage vector (README, "Ring overlap"): the prediction and the residual at every bin of the variant's curves, then
///        for its integral (point nBins).
struct LeakageLayout {
  int nBins = 0;
  int prediction(int point) const { return point; }
  int residual(int point) const { return nBins + 1 + point; }
  int size() const { return 2 * (nBins + 1); }
};

/// \brief Totals over the events of one set: the full sample, or one jackknife replica.
struct SetTotals {
  int nBins = 0;
  int nEvents = 0;
  std::array<std::vector<double>, kNBinArrays> perBin; // empty where the set has no such sums
  bool hasWeightedExtras = false;
  std::array<double, kNWeightedExtras> weightedExtra{};
  bool hasEventValues = false;
  std::array<double, kNEventValues> eventValueSum{};
};

/// \brief Estimates of the full sample and of its N leave-one-out replicas.
struct Replicated {
  std::vector<double> full;
  std::vector<std::vector<double>> replicas;
};

/// \brief Replicas of a bullet - no_bullet difference: leave-one-pair-out, and the two groups of the independent-samples jackknife.
struct DifferenceReplicas {
  std::vector<double> full;
  std::vector<std::vector<double>> paired;          // replica k drops pair k from both samples. Empty with --paired=0
  std::vector<std::vector<double>> bulletDropped;   // one bullet event dropped, no_bullet whole
  std::vector<std::vector<double>> noBulletDropped; // one no_bullet event dropped, bullet whole
};

/// \brief Value and errors of every component of an estimate vector. kNoValue marks an undefined estimate.
struct EstimateSet {
  std::vector<double> value;
  std::vector<double> error;         // jackknife. Leave-one-pair-out for a paired difference
  std::vector<double> bias;          // jackknife bias estimate, (N - 1)(mean of the replicas - value)
  std::vector<double> unpairedError; // differences only: independent-samples jackknife on the same events
};

struct JackknifeMoments {
  std::vector<double> variance; // kNoValue where the estimate, or one of its replicas, is undefined
  std::vector<double> bias;
  int nUndefinedInReplicas = 0; // components defined on the full sample but undefined in some replica
};

/// \brief Bins of each lobe of one curve, fixed once from the full-sample curve so that every replica uses the same definition.
struct LobeSpec {
  std::array<std::vector<int>, kNLobes> bins; // bins of the lobe where the full-sample curve is positive
  std::array<int, kNLobes> maxBin{{-1, -1}};  // bin of the full-sample maximum of the lobe
};
using LobeSpecs = std::array<LobeSpec, kNFeatureCurves>;

/// \brief Totals over every event of a set.
inline SetTotals sumOverEvents(const SetSums& sums)
{
  SetTotals totals;
  totals.nBins = sums.nBins;
  totals.nEvents = sums.nEvents;
  totals.hasWeightedExtras = sums.hasWeightedExtras;
  totals.hasEventValues = sums.hasEventValues;
  for (int array = 0; array < kNBinArrays; ++array) {
    if (sums.perBin[array].empty())
      continue;
    totals.perBin[array].assign(sums.nBins, 0.);
    for (int iEvent = 0; iEvent < sums.nEvents; ++iEvent)
      for (int bin = 0; bin < sums.nBins; ++bin)
        totals.perBin[array][bin] += sums.perBin[array][iEvent * sums.nBins + bin];
  }
  if (sums.hasWeightedExtras)
    for (int extra = 0; extra < kNWeightedExtras; ++extra)
      totals.weightedExtra[extra] = std::accumulate(sums.weightedExtra[extra].begin(), sums.weightedExtra[extra].end(), 0.);
  if (sums.hasEventValues)
    for (int value = 0; value < kNEventValues; ++value)
      totals.eventValueSum[value] = std::accumulate(sums.eventValue[value].begin(), sums.eventValue[value].end(), 0.);
  return totals;
}

/// \brief Totals of the jackknife replica without event iEvent: the full totals minus that event.
inline void removeEvent(const SetTotals& totals, const SetSums& sums, int iEvent, SetTotals& replica)
{
  replica.nBins = totals.nBins;
  replica.nEvents = totals.nEvents - 1;
  replica.hasWeightedExtras = totals.hasWeightedExtras;
  replica.hasEventValues = totals.hasEventValues;
  for (int array = 0; array < kNBinArrays; ++array) {
    replica.perBin[array].resize(totals.perBin[array].size());
    for (std::size_t bin = 0; bin < totals.perBin[array].size(); ++bin)
      replica.perBin[array][bin] = totals.perBin[array][bin] - sums.perBin[array][iEvent * totals.nBins + bin];
  }
  if (totals.hasWeightedExtras)
    for (int extra = 0; extra < kNWeightedExtras; ++extra)
      replica.weightedExtra[extra] = totals.weightedExtra[extra] - sums.weightedExtra[extra][iEvent];
  if (totals.hasEventValues)
    for (int value = 0; value < kNEventValues; ++value)
      replica.eventValueSum[value] = totals.eventValueSum[value] - sums.eventValue[value][iEvent];
}

/// \brief Stage-1 estimates of one set from its totals: the pooled curves (sum A / sum B per bin), their integrals, and the event means.
inline std::vector<double> computeStage1(const SetTotals& totals)
{
  const EstimateLayout layout{totals.nBins};
  std::vector<double> estimates(layout.stage1Size(), kNoValue);
  const bool hasOverlap = !totals.perBin[kBinACosZeta].empty();
  double sumW = 0., sumARing = 0., sumARingZ = 0., sumANzSq = 0., sumACosZeta = 0., sumASinZeta = 0.;
  double sumRingCurve = 0., sumRingZCurve = 0.;
  int nDefinedBins = 0;
  for (int bin = 0; bin < totals.nBins; ++bin) {
    const double binW = totals.perBin[kBinW][bin];
    const double binARing = totals.perBin[kBinARing][bin];
    const double binARingZ = totals.perBin[kBinARingZ][bin];
    const double binANzSq = totals.perBin[kBinANzSq][bin];
    sumW += binW;
    sumARing += binARing;
    sumARingZ += binARingZ;
    sumANzSq += binANzSq;
    if (hasOverlap) {
      sumACosZeta += totals.perBin[kBinACosZeta][bin];
      sumASinZeta += totals.perBin[kBinASinZeta][bin];
    }
    if (!(binW > 0.)) // no cell in this bin in any event: every estimate of the bin stays undefined
      continue;
    estimates[layout.curve(kCurveRing, bin)] = binARing / binW;
    estimates[layout.curve(kCurveRingZ, bin)] = binARingZ / binW;
    estimates[layout.curve(kCurveNzSq, bin)] = binANzSq / binW;
    if (binANzSq > kMinNzSqForAmplitude * binW)
      estimates[layout.curve(kCurveAmplitudeZ, bin)] = binARingZ / binANzSq;
    if (hasOverlap) {
      estimates[layout.curve(kCurveCosZeta, bin)] = totals.perBin[kBinACosZeta][bin] / binW;
      estimates[layout.curve(kCurveSinZeta, bin)] = totals.perBin[kBinASinZeta][bin] / binW;
    }
    sumRingCurve += binARing / binW;
    sumRingZCurve += binARingZ / binW;
    ++nDefinedBins;
  }
  if (!(sumW > 0.)) // a set without weight: nothing is defined
    return estimates;
  for (int bin = 0; bin < totals.nBins; ++bin)
    estimates[layout.curve(kCurveYieldFraction, bin)] = totals.perBin[kBinW][bin] / sumW;
  estimates[layout.scalar(kScalarRingD3N)] = sumARing / sumW;
  estimates[layout.scalar(kScalarRingZD3N)] = sumARingZ / sumW;
  estimates[layout.scalar(kScalarNzSqD3N)] = sumANzSq / sumW;
  if (sumANzSq > kMinNzSqForAmplitude * sumW)
    estimates[layout.scalar(kScalarAmplitudeZD3N)] = sumARingZ / sumANzSq;
  if (std::abs(sumARing) > kMinAbsRingForRatio * sumW)
    estimates[layout.scalar(kScalarRingZOverRing)] = sumARingZ / sumARing;
  if (hasOverlap) {
    estimates[layout.scalar(kScalarCosZetaD3N)] = sumACosZeta / sumW;
    estimates[layout.scalar(kScalarSinZetaD3N)] = sumASinZeta / sumW;
  }
  // Flat average of the pooled curve over its defined bins (every bin, on the grid binning), as the producer's Flat of one event
  estimates[layout.scalar(kScalarRingFlat)] = sumRingCurve / nDefinedBins;
  estimates[layout.scalar(kScalarRingZFlat)] = sumRingZCurve / nDefinedBins;
  // The integrated-only observables share the weight of the curves, since sumWInCuts is the sum of sumWDeltaPhi (checked when reading)
  if (totals.hasWeightedExtras)
    for (int extra = 0; extra < kNWeightedExtras; ++extra)
      estimates[layout.scalar(kExtraScalar[extra])] = totals.weightedExtra[extra] / sumW;
  if (totals.hasEventValues && totals.nEvents > 0)
    for (int value = 0; value < kNEventValues; ++value)
      estimates[layout.scalar(kEventValueScalar[value])] = totals.eventValueSum[value] / totals.nEvents;
  return estimates;
}

/// \brief Stage-1 estimates of a set and of its leave-one-out replicas. Each replica subtracts one event from the totals, so the whole
///        jackknife costs one pass per event over the bins.
inline Replicated replicateStage1(const SetSums& sums)
{
  Replicated result;
  const SetTotals totals = sumOverEvents(sums);
  result.full = computeStage1(totals);
  result.replicas.reserve(sums.nEvents);
  SetTotals replicaTotals;
  for (int iEvent = 0; iEvent < sums.nEvents; ++iEvent) {
    removeEvent(totals, sums, iEvent, replicaTotals);
    result.replicas.push_back(computeStage1(replicaTotals));
  }
  return result;
}

/// \brief Component-by-component difference of two estimate vectors, undefined where either input is.
inline std::vector<double> differenceOfVectors(const std::vector<double>& minuend, const std::vector<double>& subtrahend)
{
  std::vector<double> difference(minuend.size(), kNoValue);
  for (std::size_t index = 0; index < minuend.size(); ++index)
    if (minuend[index] != kNoValue && subtrahend[index] != kNoValue)
      difference[index] = minuend[index] - subtrahend[index];
  return difference;
}

/// \brief Bullet - no_bullet difference of two stage-1 vectors, component by component, except R_z/R. That one becomes the ratio of the
///        differences, the fraction of the jet-induced ring that the z projection sees: the difference of the two ratios would divide by the
///        no_bullet ring, which is ~0 for the jet trigger.
inline std::vector<double> differenceOfStage1(const std::vector<double>& bullet, const std::vector<double>& noBullet, const EstimateLayout& layout)
{
  std::vector<double> difference = differenceOfVectors(bullet, noBullet);
  const double ringDifference = difference[layout.scalar(kScalarRingD3N)];
  const double ringZDifference = difference[layout.scalar(kScalarRingZD3N)];
  difference[layout.scalar(kScalarRingZOverRing)] = kNoValue;
  if (ringDifference != kNoValue && ringZDifference != kNoValue && std::abs(ringDifference) > kMinAbsRingForRatio)
    difference[layout.scalar(kScalarRingZOverRing)] = ringZDifference / ringDifference;
  return difference;
}

/// \brief Replicas of the bullet - no_bullet difference. With paired = true the two inputs must hold the same pairs in the same order.
inline DifferenceReplicas differenceReplicas(const Replicated& bullet, const Replicated& noBullet, bool paired, const EstimateLayout& layout)
{
  DifferenceReplicas difference;
  difference.full = differenceOfStage1(bullet.full, noBullet.full, layout);
  if (paired) {
    if (bullet.replicas.size() != noBullet.replicas.size())
      throw std::runtime_error("paired difference with different numbers of events (internal error)");
    for (std::size_t pair = 0; pair < bullet.replicas.size(); ++pair)
      difference.paired.push_back(differenceOfStage1(bullet.replicas[pair], noBullet.replicas[pair], layout));
  }
  for (const std::vector<double>& replica : bullet.replicas)
    difference.bulletDropped.push_back(differenceOfStage1(replica, noBullet.full, layout));
  for (const std::vector<double>& replica : noBullet.replicas)
    difference.noBulletDropped.push_back(differenceOfStage1(bullet.full, replica, layout));
  return difference;
}

/// \brief Jackknife variance and bias of every component: sigma^2 = (N - 1)/N sum_i (theta_(i) - mean of theta_(i))^2, with the full-sample
///        value as the estimate (no bias correction).
inline JackknifeMoments jackknife(const std::vector<double>& full, const std::vector<std::vector<double>>& replicas)
{
  const int nReplicas = replicas.size();
  JackknifeMoments moments;
  moments.variance.assign(full.size(), kNoValue);
  moments.bias.assign(full.size(), kNoValue);
  for (std::size_t index = 0; index < full.size(); ++index) {
    if (full[index] == kNoValue)
      continue;
    double sumReplicas = 0.;
    bool allDefined = true;
    for (const std::vector<double>& replica : replicas) {
      if (replica[index] == kNoValue) {
        allDefined = false;
        break;
      }
      sumReplicas += replica[index];
    }
    if (!allDefined) {
      ++moments.nUndefinedInReplicas;
      continue;
    }
    const double meanReplica = sumReplicas / nReplicas;
    double sumSquares = 0.;
    for (const std::vector<double>& replica : replicas) {
      const double deviation = replica[index] - meanReplica;
      sumSquares += deviation * deviation;
    }
    moments.variance[index] = (nReplicas - 1.) / nReplicas * sumSquares;
    moments.bias[index] = (nReplicas - 1.) * (meanReplica - full[index]);
  }
  return moments;
}

inline std::vector<double> errorsFromVariances(const std::vector<double>& variance)
{
  std::vector<double> errors(variance.size(), kNoValue);
  for (std::size_t index = 0; index < variance.size(); ++index)
    if (variance[index] != kNoValue)
      errors[index] = std::sqrt(variance[index]);
  return errors;
}

/// \brief Values and jackknife errors of an estimate vector from its full-sample value and its replicas.
inline EstimateSet estimatesFromReplicas(const std::vector<double>& full, const std::vector<std::vector<double>>& replicas, const std::string& label, std::vector<std::string>& warnings)
{
  const JackknifeMoments moments = jackknife(full, replicas);
  if (moments.nUndefinedInReplicas > 0)
    warnings.push_back(label + ": " + std::to_string(moments.nUndefinedInReplicas) + " estimates are undefined in some jackknife replica and carry no error");
  EstimateSet estimates;
  estimates.value = full;
  estimates.error = errorsFromVariances(moments.variance);
  estimates.bias = moments.bias;
  return estimates;
}

/// \brief Lobes of Ring and RingZ from a full-sample estimate vector: the bins of each half (Delta phi < 0 and > 0, without the bins at 0 and
///        pi) where the curve is positive, and the bin of its maximum. Fixed once, so that the features are smooth functions of the replicas.
inline LobeSpecs buildLobeSpecs(const std::vector<double>& estimates, const EstimateLayout& layout, const std::string& label, std::vector<std::string>& warnings)
{
  LobeSpecs specs;
  const int middleBin = layout.middleBin();
  for (int curve = 0; curve < kNFeatureCurves; ++curve) {
    for (int lobe = 0; lobe < kNLobes; ++lobe) {
      const int firstBin = lobe == kLobeNegative ? 1 : middleBin + 1;
      const int lastBin = lobe == kLobeNegative ? middleBin - 1 : layout.nBins - 1;
      double maxValue = 0.;
      for (int bin = firstBin; bin <= lastBin; ++bin) {
        const double value = estimates[layout.curve(curve, bin)];
        if (!(value > 0.)) // also drops undefined bins (kNoValue < 0)
          continue;
        specs[curve].bins[lobe].push_back(bin);
        if (value > maxValue) {
          maxValue = value;
          specs[curve].maxBin[lobe] = bin;
        }
      }
      if (specs[curve].bins[lobe].empty())
        warnings.push_back(label + ": " + kFeatureCurveKeys[curve] + " has no positive bin for " + (lobe == kLobeNegative ? "Delta phi < 0" : "Delta phi > 0") + ", so the features of that lobe are undefined");
    }
  }
  return specs;
}

/// \brief Appends the curve features of Ring and RingZ to a stage-1 vector: lobe centroids, their separation, the heights at the fixed
///        maximum bins and the dip at Delta phi = 0 (README, "Errors"), the comparison of the RingZ lobes with the Ring lobes, then the folded
///        curve and the asymmetry R(dphi) - R(-dphi) (QA only).
inline void appendCurveFeatures(std::vector<double>& estimates, const LobeSpecs& specs, const EstimateLayout& layout, const std::vector<double>& binCentres)
{
  estimates.resize(layout.fullSize(), kNoValue);
  const int middleBin = layout.middleBin();
  for (int curve = 0; curve < kNFeatureCurves; ++curve) {
    // Lobes: centroid sum(dphi R)/sum(R) over the fixed bins, and the height at the fixed maximum bin
    for (int lobe = 0; lobe < kNLobes; ++lobe) {
      double sumValues = 0.;
      double sumPositionValues = 0.;
      bool isDefined = !specs[curve].bins[lobe].empty();
      for (const int bin : specs[curve].bins[lobe]) {
        const double value = estimates[layout.curve(curve, bin)];
        if (value == kNoValue) {
          isDefined = false;
          break;
        }
        sumValues += value;
        sumPositionValues += binCentres[bin] * value;
      }
      if (isDefined && sumValues > 0.)
        estimates[layout.lobe(curve, kLobeCentroid[lobe])] = sumPositionValues / sumValues;
      const int maxBin = specs[curve].maxBin[lobe];
      if (maxBin >= 0)
        estimates[layout.lobe(curve, kLobeHeight[lobe])] = estimates[layout.curve(curve, maxBin)];
    }
    const double centroidNeg = estimates[layout.lobe(curve, kLobeCentroidNeg)];
    const double centroidPos = estimates[layout.lobe(curve, kLobeCentroidPos)];
    if (centroidNeg != kNoValue && centroidPos != kNoValue)
      estimates[layout.lobe(curve, kLobeSeparation)] = centroidPos - centroidNeg;
    estimates[layout.lobe(curve, kLobeDip)] = estimates[layout.curve(curve, middleBin)];

    // Fold at |Delta phi| = j w: the bin at 0 and the bin at pi (bin 0, which is also -pi) stand alone, the others average with their mirror
    estimates[layout.fold(curve, 0)] = estimates[layout.curve(curve, middleBin)];
    estimates[layout.fold(curve, middleBin)] = estimates[layout.curve(curve, 0)];
    for (int step = 1; step < middleBin; ++step) {
      const double valuePos = estimates[layout.curve(curve, middleBin + step)];
      const double valueNeg = estimates[layout.curve(curve, middleBin - step)];
      if (valuePos == kNoValue || valueNeg == kNoValue)
        continue;
      estimates[layout.fold(curve, step)] = 0.5 * (valuePos + valueNeg);
      estimates[layout.asymmetry(curve, step - 1)] = valuePos - valueNeg;
    }
  }

  // RingZ against Ring: shifts of the positions (the shift of the separation is the difference of the two centroid shifts) and height ratios
  for (int position = 0; position < kNLobePositions; ++position) {
    const double ring = estimates[layout.lobe(kCurveRing, position)];
    const double ringZ = estimates[layout.lobe(kCurveRingZ, position)];
    if (ring != kNoValue && ringZ != kNoValue)
      estimates[layout.lobeComparison(kShiftCentroidNeg + position)] = ringZ - ring;
  }
  for (int lobe = 0; lobe < kNLobes; ++lobe) {
    const double ring = estimates[layout.lobe(kCurveRing, kLobeHeight[lobe])];
    const double ringZ = estimates[layout.lobe(kCurveRingZ, kLobeHeight[lobe])];
    if (ring != kNoValue && ringZ != kNoValue && ring != 0.)
      estimates[layout.lobeComparison(kHeightRatioNeg + lobe)] = ringZ / ring;
  }
}

/// \brief Estimates of one sample with their jackknife errors, the curve features included for the true trigger.
inline EstimateSet singleSampleEstimates(Replicated stage1, bool withFeatures, const EstimateLayout& layout, const std::vector<double>& binCentres, const std::string& label, std::vector<std::string>& warnings)
{
  if (withFeatures) {
    const LobeSpecs specs = buildLobeSpecs(stage1.full, layout, label, warnings);
    appendCurveFeatures(stage1.full, specs, layout, binCentres);
    for (std::vector<double>& replica : stage1.replicas)
      appendCurveFeatures(replica, specs, layout, binCentres);
  }
  return estimatesFromReplicas(stage1.full, stage1.replicas, label, warnings);
}

/// \brief Estimates of the bullet - no_bullet difference. The error is the leave-one-pair-out jackknife when the replicas are paired, and the
///        independent-samples jackknife otherwise. The latter is always kept as unpairedError, so that the gain of the pairing can be read off.
inline EstimateSet differenceEstimates(DifferenceReplicas difference, bool withFeatures, const EstimateLayout& layout, const std::vector<double>& binCentres, const std::string& label, std::vector<std::string>& warnings)
{
  if (withFeatures) {
    // The lobes of the difference curve itself, not the differences of the lobe features of the two samples
    const LobeSpecs specs = buildLobeSpecs(difference.full, layout, label, warnings);
    appendCurveFeatures(difference.full, specs, layout, binCentres);
    for (std::vector<double>& replica : difference.paired)
      appendCurveFeatures(replica, specs, layout, binCentres);
    for (std::vector<double>& replica : difference.bulletDropped)
      appendCurveFeatures(replica, specs, layout, binCentres);
    for (std::vector<double>& replica : difference.noBulletDropped)
      appendCurveFeatures(replica, specs, layout, binCentres);
  }
  // Independent samples: the two groups of replicas vary independent events, so their variances add
  const JackknifeMoments bulletDropped = jackknife(difference.full, difference.bulletDropped);
  const JackknifeMoments noBulletDropped = jackknife(difference.full, difference.noBulletDropped);
  std::vector<double> unpairedVariance(difference.full.size(), kNoValue);
  std::vector<double> unpairedBias(difference.full.size(), kNoValue);
  for (std::size_t index = 0; index < difference.full.size(); ++index) {
    if (bulletDropped.variance[index] == kNoValue || noBulletDropped.variance[index] == kNoValue)
      continue;
    unpairedVariance[index] = bulletDropped.variance[index] + noBulletDropped.variance[index];
    unpairedBias[index] = bulletDropped.bias[index] + noBulletDropped.bias[index];
  }
  EstimateSet estimates;
  if (difference.paired.empty()) {
    estimates.value = difference.full;
    estimates.error = errorsFromVariances(unpairedVariance);
    estimates.bias = unpairedBias;
  } else {
    estimates = estimatesFromReplicas(difference.full, difference.paired, label, warnings);
  }
  estimates.unpairedError = errorsFromVariances(unpairedVariance);
  return estimates;
}

/// \brief Leakage ansatz of one variant trigger set against the true trigger (README, "Ring overlap"): the prediction
///        L = R_true <cos zeta> + P_e <sin zeta>, with R_true and P_e the d3N-weighted integrals of the true trigger, and the residual R' - L,
///        at every bin of the variant's curves and for its integral.
inline std::vector<double> leakageVector(const std::vector<double>& trueEstimates, const EstimateLayout& trueLayout, const std::vector<double>& variantEstimates, const EstimateLayout& variantLayout)
{
  const LeakageLayout leakageLayout{variantLayout.nBins};
  std::vector<double> leakage(leakageLayout.size(), kNoValue);
  const double ringTrue = trueEstimates[trueLayout.scalar(kScalarRingD3N)];
  const double polE = trueEstimates[trueLayout.scalar(kScalarPolED3N)];
  if (ringTrue == kNoValue || polE == kNoValue)
    return leakage;
  for (int point = 0; point <= variantLayout.nBins; ++point) {
    const bool isIntegral = point == variantLayout.nBins;
    const double ring = variantEstimates[isIntegral ? variantLayout.scalar(kScalarRingD3N) : variantLayout.curve(kCurveRing, point)];
    const double cosZeta = variantEstimates[isIntegral ? variantLayout.scalar(kScalarCosZetaD3N) : variantLayout.curve(kCurveCosZeta, point)];
    const double sinZeta = variantEstimates[isIntegral ? variantLayout.scalar(kScalarSinZetaD3N) : variantLayout.curve(kCurveSinZeta, point)];
    if (ring == kNoValue || cosZeta == kNoValue || sinZeta == kNoValue)
      continue;
    const double prediction = ringTrue * cosZeta + polE * sinZeta;
    leakage[leakageLayout.prediction(point)] = prediction;
    leakage[leakageLayout.residual(point)] = ring - prediction;
  }
  return leakage;
}

/// \brief Leakage vectors of the full sample and of every replica. Both inputs come from the same events in the same order, so replica k of
///        the true trigger and of the variant drop the same event: the jackknife of the leakage is joint.
inline Replicated replicateLeakage(const Replicated& trueStage1, const EstimateLayout& trueLayout, const Replicated& variantStage1, const EstimateLayout& variantLayout)
{
  if (trueStage1.replicas.size() != variantStage1.replicas.size())
    throw std::runtime_error("leakage of two sets with different numbers of events (internal error)");
  Replicated leakage;
  leakage.full = leakageVector(trueStage1.full, trueLayout, variantStage1.full, variantLayout);
  for (std::size_t replica = 0; replica < trueStage1.replicas.size(); ++replica)
    leakage.replicas.push_back(leakageVector(trueStage1.replicas[replica], trueLayout, variantStage1.replicas[replica], variantLayout));
  return leakage;
}

/// \brief Results of one trigger set.
struct SetResult {
  bool present = false;
  EstimateLayout layout;
  std::array<Replicated, kNSamples> stage1;     // all events of each sample
  std::array<Replicated, kNSamples> pairStage1; // the pairs only, in pair order (--paired=1): the per-sample replicas of the paired difference
  DifferenceReplicas differenceStage1;          // the events entering the difference (the pairs, or all of them with --paired=0)
  std::array<EstimateSet, kNCurveSamples> estimates;
};

struct EstimatePoint {
  double value = kNoValue;
  double error = kNoValue;
};

/// \brief Jackknife of the ratio of one estimate between two sets of the same events (smudged over true, for instance): replica k of both
///        sets drops the same event (or pair), so the ratio is formed replica by replica.
inline EstimatePoint crossSetRatio(const std::vector<double>& numeratorFull, const std::vector<std::vector<double>>& numeratorReplicas, const std::vector<double>& denominatorFull, const std::vector<std::vector<double>>& denominatorReplicas, int index)
{
  EstimatePoint point;
  if (numeratorReplicas.size() != denominatorReplicas.size())
    throw std::runtime_error("ratio of two sets with different numbers of events (internal error)");
  if (numeratorFull[index] == kNoValue || denominatorFull[index] == kNoValue || denominatorFull[index] == 0.)
    return point;
  const std::vector<double> full = {numeratorFull[index] / denominatorFull[index]};
  std::vector<std::vector<double>> replicas;
  for (std::size_t replica = 0; replica < numeratorReplicas.size(); ++replica) {
    const double numerator = numeratorReplicas[replica][index];
    const double denominator = denominatorReplicas[replica][index];
    replicas.push_back({numerator == kNoValue || denominator == kNoValue || denominator == 0. ? kNoValue : numerator / denominator});
  }
  const JackknifeMoments moments = jackknife(full, replicas);
  point.value = full[0];
  if (moments.variance[0] != kNoValue)
    point.error = std::sqrt(moments.variance[0]);
  return point;
}

/// \brief Estimates of every trigger set of a run: each sample over all its events, and the bullet - no_bullet difference over the pairs
///        (over all events with paired = false). The true trigger carries the curve features.
/// \param label  prefix of the warnings (the run, in SummaryPlots)
inline std::array<SetResult, kNTriggerSets> computeSetResults(const std::array<SampleData, kNSamples>& samples, bool hasNoBullet, bool paired, const std::array<std::vector<int>, kNSamples>& pairRows, const std::vector<double>& deltaPhiCentres, const std::string& label, std::vector<std::string>& warnings)
{
  std::array<SetResult, kNTriggerSets> results;
  for (int set = 0; set < kNTriggerSets; ++set) {
    if (!samples[kSampleBullet].sets[set].present) // the configurations of the samples agree, so the no_bullet file has the same sets
      continue;
    SetResult& result = results[set];
    result.present = true;
    result.layout.nBins = samples[kSampleBullet].sets[set].nBins;
    const bool withFeatures = set == kSetTrue;
    for (int sample = 0; sample < kNSamples; ++sample) {
      if (sample == kSampleNoBullet && !hasNoBullet)
        continue;
      result.stage1[sample] = replicateStage1(samples[sample].sets[set]);
      result.estimates[sample] = singleSampleEstimates(result.stage1[sample], withFeatures, result.layout, deltaPhiCentres, label + kSampleLabels[sample] + ", " + kSetLabels[set], warnings);
    }
    if (!hasNoBullet)
      continue;
    if (paired) {
      for (int sample = 0; sample < kNSamples; ++sample)
        result.pairStage1[sample] = replicateStage1(selectEvents(samples[sample].sets[set], pairRows[sample]));
      result.differenceStage1 = differenceReplicas(result.pairStage1[kSampleBullet], result.pairStage1[kSampleNoBullet], true, result.layout);
    } else {
      result.differenceStage1 = differenceReplicas(result.stage1[kSampleBullet], result.stage1[kSampleNoBullet], false, result.layout);
    }
    result.estimates[kSampleDifference] = differenceEstimates(result.differenceStage1, withFeatures, result.layout, deltaPhiCentres, label + "bullet - no_bullet, " + kSetLabels[set], warnings);
  }
  return results;
}

/// \brief Leakage ansatz of every variant set (README, "Ring overlap"), [set][sample], jackknifed jointly with the true trigger. Empty where
///        not computed: the difference needs the per-sample replicas of the pairs, so paired = true.
inline std::array<std::array<EstimateSet, kNCurveSamples>, kNTriggerSets> computeLeakage(const std::array<SetResult, kNTriggerSets>& results, bool hasNoBullet, bool paired, const std::string& label, std::vector<std::string>& warnings)
{
  std::array<std::array<EstimateSet, kNCurveSamples>, kNTriggerSets> leakage;
  const EstimateLayout& trueLayout = results[kSetTrue].layout;
  for (const int set : {kSetScan, kSetMixing, kSetSmudge}) {
    if (!results[set].present)
      continue;
    const EstimateLayout& variantLayout = results[set].layout;
    for (int sample = 0; sample < kNSamples; ++sample) {
      if (sample == kSampleNoBullet && !hasNoBullet)
        continue;
      const Replicated sampleLeakage = replicateLeakage(results[kSetTrue].stage1[sample], trueLayout, results[set].stage1[sample], variantLayout);
      leakage[set][sample] = estimatesFromReplicas(sampleLeakage.full, sampleLeakage.replicas, label + kSampleLabels[sample] + ", leakage of the " + kSetLabels[set], warnings);
    }
    if (!hasNoBullet)
      continue;
    if (!paired) {
      warnings.push_back(label + "RingOverlap/: --paired=0, so the leakage test of the bullet - no_bullet difference (" + kSetLabels[set] + ") is not computed");
      continue;
    }
    const Replicated bulletLeakage = replicateLeakage(results[kSetTrue].pairStage1[kSampleBullet], trueLayout, results[set].pairStage1[kSampleBullet], variantLayout);
    const Replicated noBulletLeakage = replicateLeakage(results[kSetTrue].pairStage1[kSampleNoBullet], trueLayout, results[set].pairStage1[kSampleNoBullet], variantLayout);
    std::vector<std::vector<double>> differenceReplicaList;
    for (std::size_t pair = 0; pair < bulletLeakage.replicas.size(); ++pair)
      differenceReplicaList.push_back(differenceOfVectors(bulletLeakage.replicas[pair], noBulletLeakage.replicas[pair]));
    leakage[set][kSampleDifference] = estimatesFromReplicas(differenceOfVectors(bulletLeakage.full, noBulletLeakage.full), differenceReplicaList, label + "bullet - no_bullet, leakage of the " + kSetLabels[set], warnings);
  }
  return leakage;
}

/// \brief sigma x sqrt(N) of one estimate for an ideal detector (README, "Expected precision"): sqrt(3 <|e|^2>)/|alpha|. For one bin, pass
///        <|e|^2> over the yield fraction of the bin.
inline double sigmaSqrtNOf(double meanESquared, double alphaLambda) { return std::sqrt(3. * meanESquared) / std::abs(alphaLambda); }

/// \brief <|e|^2> of the integrals of kPrecisionKeys (1, n_z^2 and n_z^2 sin^2 theta), from the true-trigger estimates of the bullet sample.
inline std::array<double, 3> integratedESquared(const std::vector<double>& bulletTrueValues, const EstimateLayout& layout)
{
  return {1., bulletTrueValues[layout.scalar(kScalarNzSqD3N)], bulletTrueValues[layout.scalar(kScalarNzSqSinSqThetaD3N)]};
}

// ======================================================================================================================================
// Output objects
// ======================================================================================================================================

/// \brief Output folder: the canvases at its top, for a first skim, and every object they draw in its Parts/ subfolder.
struct OutputFolder {
  TDirectory* canvases = nullptr;
  TDirectory* parts = nullptr;
};

inline OutputFolder makeFolder(TDirectory* parent, const std::string& name)
{
  OutputFolder folder;
  folder.canvases = parent->mkdir(name.c_str());
  folder.parts = folder.canvases->mkdir("Parts");
  return folder;
}

/// \brief A drawable quantity: the histogram stored in Parts/, and which of its bins hold a defined value. An undefined estimate is an empty
///        bin of the histogram (content and error 0), and is left out of every graph drawn from it.
struct Part {
  TH1* hist = nullptr;
  std::vector<bool> defined; // per bin, from 0
};

struct DrawEntry {
  Part part;
  std::string label;
  Look look;
  DrawStyle style = kStyleMarker;
};

/// \brief Indices offset, offset + 1, ..., offset + count - 1 of an estimate vector (one curve, or one block).
inline std::vector<int> blockIndices(int offset, int count)
{
  std::vector<int> indices(count);
  std::iota(indices.begin(), indices.end(), offset);
  return indices;
}

/// \brief Fills bins 1, 2, ... of a histogram with the given components of an estimate vector, scaled, and returns which are defined. An
///        undefined estimate leaves its bin empty, and is logged.
inline std::vector<bool> fillFromEstimates(TH1* hist, const std::vector<double>& values, const std::vector<double>& errors, const std::vector<int>& indices, double scale, std::vector<std::string>& warnings)
{
  std::vector<bool> defined(indices.size(), false);
  int nUndefined = 0;
  for (std::size_t position = 0; position < indices.size(); ++position) {
    const double value = values[indices[position]];
    const double error = errors[indices[position]];
    if (value == kNoValue) {
      ++nUndefined;
      continue;
    }
    defined[position] = true;
    hist->SetBinContent(position + 1, scale * value);
    hist->SetBinError(position + 1, error == kNoValue ? 0. : scale * error); // an undefined error is already reported by the jackknife
  }
  if (nUndefined > 0)
    warnings.push_back(std::string(hist->GetDirectory()->GetPath()) + "/" + hist->GetName() + ": " + std::to_string(nUndefined) + " undefined estimates, left as empty bins and out of the graphs");
  return defined;
}

/// \brief One curve of an estimate vector on the given axis, stored in the directory.
inline Part curveHistogram(TDirectory* directory, const std::string& name, const std::string& title, const AxisRange& range, const std::vector<double>& values, const std::vector<double>& errors, int offset, double scale, std::vector<std::string>& warnings)
{
  directory->cd();
  Part part;
  part.hist = new TH1D(name.c_str(), title.c_str(), range.nBins, range.low, range.high);
  part.defined = fillFromEstimates(part.hist, values, errors, blockIndices(offset, range.nBins), scale, warnings);
  return part;
}

/// \brief A table: one labelled bin per estimate, stored in the directory.
inline Part tableHistogram(TDirectory* directory, const std::string& name, const std::string& title, const std::vector<std::string>& labels, const std::vector<double>& values, const std::vector<double>& errors, const std::vector<int>& indices, double scale, std::vector<std::string>& warnings)
{
  directory->cd();
  const int nLabels = labels.size();
  Part part;
  part.hist = new TH1D(name.c_str(), title.c_str(), nLabels, 0., nLabels);
  for (int bin = 0; bin < nLabels; ++bin)
    part.hist->GetXaxis()->SetBinLabel(bin + 1, labels[bin].c_str());
  part.defined = fillFromEstimates(part.hist, values, errors, indices, scale, warnings);
  return part;
}

// ======================================================================================================================================
// Drawing
// ======================================================================================================================================

/// \brief Runs of consecutive defined bins of a part, as [first, last] bins from 0. A band or a line is drawn per run, so that it never
///        bridges an undefined bin.
inline std::vector<std::pair<int, int>> definedRuns(const Part& part)
{
  std::vector<std::pair<int, int>> runs;
  const int nBins = part.defined.size();
  for (int bin = 0; bin < nBins; ++bin) {
    if (!part.defined[bin])
      continue;
    if (runs.empty() || runs.back().second != bin - 1)
      runs.emplace_back(bin, bin);
    else
      runs.back().second = bin;
  }
  return runs;
}

/// \brief Graph of the bins firstBin to lastBin of a part, at the bin centres shifted by xOffset bin widths, owned by the pad it is drawn on.
inline TGraphErrors* graphOfRun(const Part& part, int firstBin, int lastBin, double xOffset)
{
  const TAxis* axis = part.hist->GetXaxis();
  auto* graph = new TGraphErrors(lastBin - firstBin + 1);
  graph->SetBit(kCanDelete);
  for (int bin = firstBin; bin <= lastBin; ++bin) {
    graph->SetPoint(bin - firstBin, axis->GetBinCenter(bin + 1) + xOffset * axis->GetBinWidth(bin + 1), part.hist->GetBinContent(bin + 1));
    graph->SetPointError(bin - firstBin, 0., part.hist->GetBinError(bin + 1));
  }
  return graph;
}

inline void styleMarkers(TGraphErrors* graph, const Look& look)
{
  graph->SetMarkerStyle(look.marker);
  graph->SetMarkerColor(look.color);
  graph->SetMarkerSize(1.);
  graph->SetLineColor(look.color);
  graph->SetLineWidth(2);
}

/// \brief Draws the entries on the current pad: a frame, every band first, the reference and separator lines over them, then the lines,
///        markers and histograms, the legend, and the axes once more (the opaque bands cover their ticks). The y range covers every
///        defined point with its error and leaves room for the legend. Every drawn object is a copy owned by the pad: the stored parts carry
///        no drawing style, and one part can appear in several canvases with different looks.
/// \param referenceY   dashed horizontal line (0 for curves, 1 for ratios), or kNoValue for none
/// \param separators   dashed vertical lines at these x: the groups of a table
/// \param textScale    axis labels, titles and legend x this: the small pad of a two-pad canvas
/// \param showXLabels  false for the upper pad of a two-pad canvas, which shares the x axis of the lower one
inline void drawEntriesOnPad(const std::vector<DrawEntry>& entries, const std::string& title, double referenceY, const std::vector<double>& separators, double textScale, bool showXLabels)
{
  double yMin = std::numeric_limits<double>::max();
  double yMax = std::numeric_limits<double>::lowest();
  bool onlyHistograms = true;
  for (const DrawEntry& entry : entries) {
    if (entry.style != kStyleHistogram)
      onlyHistograms = false;
    const double errorFactor = entry.style == kStyleHistogram ? 0. : 1.;
    for (std::size_t bin = 0; bin < entry.part.defined.size(); ++bin) {
      if (!entry.part.defined[bin])
        continue;
      const double content = entry.part.hist->GetBinContent(bin + 1);
      const double error = errorFactor * entry.part.hist->GetBinError(bin + 1);
      yMin = std::min(yMin, content - error);
      yMax = std::max(yMax, content + error);
    }
  }
  if (referenceY != kNoValue) {
    yMin = std::min(yMin, referenceY);
    yMax = std::max(yMax, referenceY);
  }
  if (yMin > yMax) { // nothing defined (already logged): an empty frame
    yMin = 0.;
    yMax = 1.;
  }
  double span = yMax - yMin;
  if (!(span > 0.))
    span = yMax != 0. ? std::abs(yMax) : 1.;
  const double bottom = onlyHistograms && yMin >= 0. ? 0. : yMin - 0.1 * span;

  // Frame: an empty copy of the first histogram (axis, bin labels and axis titles) drawn with an invisible line, so only the axes and the
  // title show. Clone() registers copies in the current directory, so every copy is detached at once and owned by the pad instead
  auto* frame = static_cast<TH1*>(entries.front().part.hist->Clone("frame"));
  frame->SetDirectory(nullptr);
  frame->SetBit(kCanDelete);
  frame->Reset();
  frame->SetTitle(title.c_str());
  frame->SetLineWidth(0);
  frame->SetStats(false);
  TAxis* xAxis = frame->GetXaxis();
  TAxis* yAxis = frame->GetYaxis();
  xAxis->SetLabelSize(showXLabels ? kLabelSize * textScale : 0.);
  xAxis->SetTitleSize(showXLabels ? kTitleSize * textScale : 0.);
  yAxis->SetLabelSize(kLabelSize * textScale);
  yAxis->SetTitleSize(kTitleSize * textScale);
  yAxis->SetTitleOffset(kYTitleOffset / textScale);
  if (textScale > 1.) // a small pad: fewer, larger labels, and no room above the axis for a common exponent
    yAxis->SetNdivisions(505);
  else
    yAxis->SetMaxDigits(3);
  if (xAxis->GetLabels() != nullptr && xAxis->GetNbins() > kMaxHorizontalLabels) {
    frame->LabelsOption("v", "X");
    gPad->SetBottomMargin(kTableBottomMargin);
  }
  // The legend hangs from the top of the frame, one row per entry: the y range keeps the data below it
  const double frameTop = 1. - gPad->GetTopMargin();
  const double frameHeight = frameTop - gPad->GetBottomMargin();
  const double legendTop = frameTop - 0.01;
  const double legendHeight = std::min(0.5 * frameHeight, kLegendRowHeight * textScale * entries.size());
  const double dataFraction = (frameHeight - legendHeight - 0.03) / frameHeight;
  const double top = bottom + (yMax - bottom) / dataFraction;
  frame->SetMinimum(bottom);
  frame->SetMaximum(top);
  frame->Draw("HIST");

  // 1. Bands. "3" draws no outline, so the line colour of a band graph shows in its legend entry only
  std::vector<TObject*> legendObjects(entries.size(), nullptr);
  for (std::size_t index = 0; index < entries.size(); ++index) {
    const DrawEntry& entry = entries[index];
    if (entry.style != kStyleBand && entry.style != kStyleBandOnly)
      continue;
    for (const auto& [firstBin, lastBin] : definedRuns(entry.part)) {
      TGraphErrors* band = graphOfRun(entry.part, firstBin, lastBin, 0.);
      band->SetFillColor(entry.look.bandColor);
      band->SetFillStyle(1001);
      band->SetLineColor(entry.style == kStyleBand ? entry.look.color : entry.look.bandColor);
      band->SetLineWidth(entry.style == kStyleBand ? 2 : 1);
      if (firstBin == lastBin) { // a single point has no band: a marker with its error bar
        styleMarkers(band, entry.look);
        band->Draw("P");
      } else {
        band->Draw("3");
      }
      if (legendObjects[index] == nullptr)
        legendObjects[index] = band;
    }
  }

  // 2. Reference line and separators, over the bands
  if (referenceY != kNoValue) {
    auto* line = new TLine(xAxis->GetXmin(), referenceY, xAxis->GetXmax(), referenceY);
    line->SetBit(kCanDelete);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kGray + 2);
    line->Draw();
  }
  for (const double separatorX : separators) {
    auto* line = new TLine(separatorX, bottom, separatorX, top);
    line->SetBit(kCanDelete);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kGray + 2);
    line->Draw();
  }

  // 3. Lines, markers (side by side within a bin when several tables share the canvas) and histograms
  int nMarkerEntries = 0;
  for (const DrawEntry& entry : entries)
    if (entry.style == kStyleMarker)
      ++nMarkerEntries;
  const double markerStep = nMarkerEntries > 1 ? std::min(0.15, 0.6 / nMarkerEntries) : 0.;
  int markerIndex = 0;
  for (std::size_t index = 0; index < entries.size(); ++index) {
    const DrawEntry& entry = entries[index];
    if (entry.style == kStyleHistogram) {
      auto* copy = static_cast<TH1*>(entry.part.hist->Clone());
      copy->SetDirectory(nullptr);
      copy->SetBit(kCanDelete);
      copy->SetLineColor(entry.look.color);
      copy->SetLineWidth(2);
      copy->SetFillStyle(0);
      copy->Draw("HIST SAME");
      legendObjects[index] = copy;
      continue;
    }
    if (entry.style == kStyleBandOnly)
      continue;
    const double xOffset = entry.style == kStyleMarker ? (markerIndex - 0.5 * (nMarkerEntries - 1)) * markerStep : 0.;
    if (entry.style == kStyleMarker)
      ++markerIndex;
    for (const auto& [firstBin, lastBin] : definedRuns(entry.part)) {
      TGraphErrors* graph = graphOfRun(entry.part, firstBin, lastBin, xOffset);
      if (entry.style == kStyleMarker || firstBin == lastBin) {
        styleMarkers(graph, entry.look);
        graph->Draw("P");
      } else {
        graph->SetLineColor(entry.look.color);
        graph->SetLineWidth(2);
        graph->Draw("LX"); // X: no error bars, the band shows them
      }
      if (legendObjects[index] == nullptr)
        legendObjects[index] = graph;
    }
  }

  // 4. Legend
  auto* legend = new TLegend(0.55, legendTop - legendHeight, 0.89, legendTop);
  legend->SetBit(kCanDelete);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.03 * textScale);
  for (std::size_t index = 0; index < entries.size(); ++index) {
    const DrawEntry& entry = entries[index];
    if (legendObjects[index] == nullptr) { // every bin undefined (already logged)
      legend->AddEntry(static_cast<TObject*>(nullptr), (entry.label + " (undefined)").c_str(), "");
      continue;
    }
    const std::array<std::string, 5> legendOptions = {"lf", "f", "l", "pe", "l"};
    legend->AddEntry(legendObjects[index], entry.label.c_str(), legendOptions[entry.style].c_str());
  }
  legend->Draw();
  gPad->RedrawAxis();
}

/// \brief Canvas of overlaid entries (drawEntriesOnPad), written into the directory.
inline void drawOverlay(TDirectory* directory, const std::string& name, const std::string& title, const std::vector<DrawEntry>& entries, double referenceY, const std::vector<double>& separators = {})
{
  if (entries.empty())
    return;
  directory->cd();
  TCanvas canvas(name.c_str(), title.c_str(), 900, 650);
  canvas.SetLeftMargin(0.12);
  drawEntriesOnPad(entries, title, referenceY, separators, 1., true);
  directory->WriteTObject(&canvas);
}

#endif // HYDROTRUTHLABCOMMON_H