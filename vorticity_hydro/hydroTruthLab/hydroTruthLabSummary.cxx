// hydroTruthLabSummary.cxx
//
// Combines the producer files of hydroTruthLab.cxx (bullet and no_bullet) into the results of the study: every pooled quantity with its
// jackknife error from the event trees, the bullet - no_bullet differences over paired initial conditions, the checks of the trigger variants
// (the leakage ansatz of the ring overlap included), the polarization maps as vector fields, and the analogues of the figures of
// PhysRevC.109.014905 for this dataset.
//
// Usage:
//   ./hydroTruthLabSummary.exe <bulletFile> <noBulletFile|none> <outputFile> [--option=value ...]      (--help lists every option)
//   Normally compiled and launched by runHydroTruthLab.sh, which holds every configurable.
//
// Needs: ROOT (Hist, Tree, RIO, Gpad, Graf) and producer files written with the same configuration.
//
// Estimators, errors and the output reference: README.md in this folder ("Errors" and "Summary").

#include <TArrow.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TNamed.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
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
constexpr double kPercent = 100.;            // Paper/ and PolarizationMaps/ only: the paper and the consumer's maps use percent
constexpr double kSumCheckTolerance = 1e-12; // relative. Sum of sumWDeltaPhi over bins vs sumWInCuts: the same numbers, summed in the same order
constexpr double kAxisTolerance = 1e-9;      // in units of the bin width, for the bins centred at Delta phi = 0 and delta = 0
constexpr double kRebuildTolerance = 1e-9;   // relative. Mixed weights rebuilt from the tree vs those the producer filled: same numbers, other order
constexpr int kPullBins = 50;
constexpr double kPullRange = 5.;
constexpr int kArrowBlockSize = 4;             // PolarizationMaps/: each arrow averages 4 x 4 bins (0.4 GeV/c), as the consumer's maps
constexpr double kArrowScalePercentile = 0.95; // arrows are scaled to this quantile of the tile magnitudes, and capped there
constexpr double kPi = TMath::Pi();

// Drawing
constexpr double kLabelSize = 0.04;
constexpr double kTitleSize = 0.045;
constexpr double kYTitleOffset = 1.25;
constexpr double kLegendRowHeight = 0.05;   // NDC per legend entry
constexpr int kMaxHorizontalLabels = 4;     // tables with more bin labels than this get them vertical
constexpr double kTableBottomMargin = 0.32; // room for the vertical labels
constexpr double kLowerPadHeight = 0.3;     // cAmplitudeZ: the difference pad, as a fraction of the canvas height

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
// Names used in the summary file
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
// cAmplitudeZ: A_z next to the R of the same sample
const std::array<Look, kNCurveSamples> kAmplitudeZLooks = {{{kOrange + 7, kOrange - 9, kOpenCircle}, {kCyan + 2, kCyan - 9, kOpenSquare}, {kGreen + 2, kGreen - 9, kOpenTriangleUp}}};
// Overlays of several observables of one sample
const std::array<Look, 6> kPaletteLooks = {{{kBlack, kGray, kFullCircle}, {kRed + 1, kRed - 9, kFullSquare}, {kBlue + 1, kBlue - 9, kFullTriangleUp}, {kGreen + 2, kGreen - 9, kOpenCircle}, {kMagenta + 1, kMagenta - 9, kOpenSquare}, {kOrange + 7, kOrange - 9, kOpenTriangleUp}}};
// Difference/: the wide independent-samples band behind the paired one
const Look kUnpairedBandLook = {kBlack, kAzure - 9, kOpenCircle};

enum DrawStyle { kStyleBand = 0, // error band and central line: every curve with errors
                 kStyleBandOnly, // error band without its line: the wider band of nested bands
                 kStyleLine,     // central line only: curves without errors
                 kStyleMarker,   // markers with error bars: tables
                 kStyleHistogram // the histogram itself: distributions
};

// ======================================================================================================================================
// Configuration
// ======================================================================================================================================

struct SummaryConfig {
  std::string bulletFile;
  std::string noBulletFile; // "none" skips everything that needs the no_bullet sample
  std::string outputFile;
  bool paired = true;
  double alphaLambda = 0.749;
  double nLambdaPrecision = 0.; // 0: Precision/ holds sigma x sqrt(N) only
  std::string centralityLabel = "40-50%";
};

void printUsage()
{
  std::cout << "Usage: ./hydroTruthLabSummary.exe <bulletFile> <noBulletFile|none> <outputFile> [--option=value ...]\n"
            << "  --paired=1|0             paired differences (key hypothesis), or differences of means         (default: 1)\n"
            << "  --alphaLambda=<value>    Lambda decay parameter, used only in Precision/                      (default: 0.749)\n"
            << "  --nLambdaPrecision=<N>   number of Lambdas for sigma and the significance in Precision/, 0 = off (default: 0)\n"
            << "  --centralityLabel=<text> label of the Fig. 11 point                                           (default: 40-50%)\n"
            << "With none as the no_bullet file, everything that needs it (Difference/, Fig. 3 left, significances) is skipped.\n"
            << "See README.md for the definitions." << std::endl;
}

// The std::sto* functions throw with messages like "stoi", which say nothing about the option: these wrappers name it (same as the producer)
int toInt(const std::string& key, const std::string& value)
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

double toDouble(const std::string& key, const std::string& value)
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
std::string readConfigValue(TFile* file, const std::string& name) { return getObject<TNamed>(file, "Config/" + name)->GetTitle(); }

/// \brief SetBranchAddress with a hard stop on a missing branch (a producer file older than the summary, for instance).
void attachBranch(TTree* tree, const std::string& name, void* address)
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

/// \brief Everything the summary needs from one event tree.
struct SampleData {
  std::vector<int> eventId;
  std::vector<double> phiTrigger;
  std::vector<double> sumWInCuts;
  std::array<SetSums, kNTriggerSets> sets;
  // Delta-scan at delta = 0 against the true trigger, event by event: max |difference| and max |value|, indexed like BinArray (kBinW unused)
  std::array<double, kNBinArrays> scanZeroMaxResidual{};
  std::array<double, kNBinArrays> scanZeroMaxValue{};
};

/// \brief Reads the event tree of one producer file into the per-event sums of every trigger set, with entry-by-entry consistency checks.
SampleData readSample(TFile* file, const std::string& label, int nDeltaPhiBins, int nDelta, int scanZeroBin, bool hasMixing, bool hasSmudge)
{
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
  Int_t nTrueBins = 0, nScanBins = 0, nMixBins = 0, nSmudgeBins = 0;
  std::array<Double_t, kNWeightedExtras> extraD3N{};
  std::array<Double_t, kNEventValues> trueValues{}, mixValues{}, smudgeValues{};
  std::array<std::vector<Double_t>, kNBinArrays> trueBins, scanValues, mixBins, smudgeBins;

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
  if (hasMixing) {
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
    if (nTrueBins != nDeltaPhiBins || nScanBins != nDelta || (hasMixing && nMixBins != nDeltaPhiBins) || (hasSmudge && nSmudgeBins != nDeltaPhiBins))
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
SetSums selectEvents(const SetSums& sums, const std::vector<int>& rows)
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
SetTotals sumOverEvents(const SetSums& sums)
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
void removeEvent(const SetTotals& totals, const SetSums& sums, int iEvent, SetTotals& replica)
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
std::vector<double> computeStage1(const SetTotals& totals)
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
Replicated replicateStage1(const SetSums& sums)
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
std::vector<double> differenceOfVectors(const std::vector<double>& minuend, const std::vector<double>& subtrahend)
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
std::vector<double> differenceOfStage1(const std::vector<double>& bullet, const std::vector<double>& noBullet, const EstimateLayout& layout)
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
DifferenceReplicas differenceReplicas(const Replicated& bullet, const Replicated& noBullet, bool paired, const EstimateLayout& layout)
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
JackknifeMoments jackknife(const std::vector<double>& full, const std::vector<std::vector<double>>& replicas)
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

std::vector<double> errorsFromVariances(const std::vector<double>& variance)
{
  std::vector<double> errors(variance.size(), kNoValue);
  for (std::size_t index = 0; index < variance.size(); ++index)
    if (variance[index] != kNoValue)
      errors[index] = std::sqrt(variance[index]);
  return errors;
}

/// \brief Values and jackknife errors of an estimate vector from its full-sample value and its replicas.
EstimateSet estimatesFromReplicas(const std::vector<double>& full, const std::vector<std::vector<double>>& replicas, const std::string& label, std::vector<std::string>& warnings)
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
LobeSpecs buildLobeSpecs(const std::vector<double>& estimates, const EstimateLayout& layout, const std::string& label, std::vector<std::string>& warnings)
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
void appendCurveFeatures(std::vector<double>& estimates, const LobeSpecs& specs, const EstimateLayout& layout, const std::vector<double>& binCentres)
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
EstimateSet singleSampleEstimates(Replicated stage1, bool withFeatures, const EstimateLayout& layout, const std::vector<double>& binCentres, const std::string& label, std::vector<std::string>& warnings)
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
EstimateSet differenceEstimates(DifferenceReplicas difference, bool withFeatures, const EstimateLayout& layout, const std::vector<double>& binCentres, const std::string& label, std::vector<std::string>& warnings)
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
std::vector<double> leakageVector(const std::vector<double>& trueEstimates, const EstimateLayout& trueLayout, const std::vector<double>& variantEstimates, const EstimateLayout& variantLayout)
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
Replicated replicateLeakage(const Replicated& trueStage1, const EstimateLayout& trueLayout, const Replicated& variantStage1, const EstimateLayout& variantLayout)
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
EstimatePoint crossSetRatio(const std::vector<double>& numeratorFull, const std::vector<std::vector<double>>& numeratorReplicas, const std::vector<double>& denominatorFull, const std::vector<std::vector<double>>& denominatorReplicas, int index)
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

// ======================================================================================================================================
// Output objects
// ======================================================================================================================================

/// \brief Output folder: the canvases at its top, for a first skim, and every object they draw in its Parts/ subfolder.
struct OutputFolder {
  TDirectory* canvases = nullptr;
  TDirectory* parts = nullptr;
};

OutputFolder makeFolder(TDirectory* parent, const std::string& name)
{
  OutputFolder folder;
  folder.canvases = parent->mkdir(name.c_str());
  folder.parts = folder.canvases->mkdir("Parts");
  return folder;
}

struct AxisRange {
  int nBins = 0;
  double low = 0.;
  double high = 0.;
  std::string title;
};

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
std::vector<int> blockIndices(int offset, int count)
{
  std::vector<int> indices(count);
  std::iota(indices.begin(), indices.end(), offset);
  return indices;
}

/// \brief Fills bins 1, 2, ... of a histogram with the given components of an estimate vector, scaled, and returns which are defined. An
///        undefined estimate leaves its bin empty, and is logged.
std::vector<bool> fillFromEstimates(TH1* hist, const std::vector<double>& values, const std::vector<double>& errors, const std::vector<int>& indices, double scale, std::vector<std::string>& warnings)
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
Part curveHistogram(TDirectory* directory, const std::string& name, const std::string& title, const AxisRange& range, const std::vector<double>& values, const std::vector<double>& errors, int offset, double scale, std::vector<std::string>& warnings)
{
  directory->cd();
  Part part;
  part.hist = new TH1D(name.c_str(), title.c_str(), range.nBins, range.low, range.high);
  part.defined = fillFromEstimates(part.hist, values, errors, blockIndices(offset, range.nBins), scale, warnings);
  return part;
}

/// \brief A table: one labelled bin per estimate, stored in the directory.
Part tableHistogram(TDirectory* directory, const std::string& name, const std::string& title, const std::vector<std::string>& labels, const std::vector<double>& values, const std::vector<double>& errors, const std::vector<int>& indices, double scale, std::vector<std::string>& warnings)
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

/// \brief Copy of a producer profile in the directory: the curves without per-bin sums in the tree, with their pooled-master errors. A bin
///        without entries is undefined.
Part profilePart(TDirectory* directory, TFile* file, const std::string& path, const std::string& newName)
{
  const TProfile* source = getObject<TProfile>(file, path);
  directory->cd(); // Clone() registers the copy in the current directory before SetDirectory() moves it: keep that in the target
  auto* copy = static_cast<TProfile*>(source->Clone(newName.c_str()));
  copy->SetDirectory(directory);
  Part part;
  part.hist = copy;
  part.defined.assign(copy->GetNbinsX(), false);
  for (int bin = 0; bin < copy->GetNbinsX(); ++bin)
    part.defined[bin] = copy->GetBinEntries(bin + 1) > 0.;
  return part;
}

/// \brief A histogram drawn as it is (a distribution): every bin is defined.
Part histogramPart(TH1* hist)
{
  Part part;
  part.hist = hist;
  part.defined.assign(hist->GetNbinsX(), true);
  return part;
}

/// \brief TH2D of the means of a producer TProfile2D (errors kept), scaled, stored in the directory.
TH2D* profileMap(TDirectory* directory, const std::string& name, const std::string& title, const TProfile2D* profile, double scale)
{
  directory->cd();
  TH2D* map = profile->ProjectionXY(name.c_str(), "E");
  map->SetDirectory(directory);
  map->SetTitle(title.c_str());
  map->Scale(scale);
  return map;
}

/// \brief Copy of a producer TH2D in the directory.
TH2D* copyMap(TDirectory* directory, TFile* file, const std::string& path, const std::string& newName, const std::string& titleSuffix)
{
  const TH2D* source = getObject<TH2D>(file, path);
  directory->cd();
  auto* copy = static_cast<TH2D*>(source->Clone(newName.c_str()));
  copy->SetDirectory(directory);
  copy->SetTitle((std::string(source->GetTitle()) + titleSuffix).c_str());
  return copy;
}

/// \brief Canvas of one 2D map (COLZ), written into the directory. The map itself is stored in Parts/.
void drawMap(TDirectory* directory, const std::string& name, TH2D* map)
{
  directory->cd();
  TCanvas canvas(name.c_str(), map->GetTitle(), 900, 700);
  canvas.SetRightMargin(0.16);
  map->Draw("COLZ");
  directory->WriteTObject(&canvas);
}

// ======================================================================================================================================
// Drawing
// ======================================================================================================================================

/// \brief Runs of consecutive defined bins of a part, as [first, last] bins from 0. A band or a line is drawn per run, so that it never
///        bridges an undefined bin.
std::vector<std::pair<int, int>> definedRuns(const Part& part)
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
TGraphErrors* graphOfRun(const Part& part, int firstBin, int lastBin, double xOffset)
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

void styleMarkers(TGraphErrors* graph, const Look& look)
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
void drawEntriesOnPad(const std::vector<DrawEntry>& entries, const std::string& title, double referenceY, const std::vector<double>& separators, double textScale, bool showXLabels)
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
void drawOverlay(TDirectory* directory, const std::string& name, const std::string& title, const std::vector<DrawEntry>& entries, double referenceY, const std::vector<double>& separators = {})
{
  if (entries.empty())
    return;
  directory->cd();
  TCanvas canvas(name.c_str(), title.c_str(), 900, 650);
  canvas.SetLeftMargin(0.12);
  drawEntriesOnPad(entries, title, referenceY, separators, 1., true);
  directory->WriteTObject(&canvas);
}

/// \brief Canvas with the entries on an upper pad and their bullet - no_bullet differences on a lower pad (30% of the height) sharing its
///        x axis, written into the directory.
void drawWithDifferencePad(TDirectory* directory, const std::string& name, const std::string& title, const std::vector<DrawEntry>& upperEntries, const std::vector<DrawEntry>& lowerEntries, const std::string& lowerYTitle)
{
  directory->cd();
  TCanvas canvas(name.c_str(), title.c_str(), 900, 850);
  auto* upperPad = new TPad("upperPad", "", 0., kLowerPadHeight, 1., 1.);
  auto* lowerPad = new TPad("lowerPad", "", 0., 0., 1., kLowerPadHeight);
  upperPad->SetBit(kCanDelete);
  lowerPad->SetBit(kCanDelete);
  upperPad->SetLeftMargin(0.12);
  upperPad->SetBottomMargin(0.02);
  lowerPad->SetLeftMargin(0.12);
  lowerPad->SetTopMargin(0.04);
  lowerPad->SetBottomMargin(0.3);
  upperPad->Draw();
  lowerPad->Draw();
  upperPad->cd();
  drawEntriesOnPad(upperEntries, title, 0., {}, 1., false);
  lowerPad->cd();
  // The x title is the lower pad's own: its frame title sets it explicitly, since an empty field would clear it
  const std::string xTitle = lowerEntries.front().part.hist->GetXaxis()->GetTitle();
  drawEntriesOnPad(lowerEntries, ";" + xTitle + ";" + lowerYTitle, 0., {}, (1. - kLowerPadHeight) / kLowerPadHeight, true);
  directory->WriteTObject(&canvas);
}

// --------------------------------------------------------------------------------------------------------------------------------------
// Polarization maps. The vector field is ported from DrawVectorFieldPanel() of auxiliaryPerConfigPlots.cxx (itself from
// plotHelicityEfficiency.cxx), with three changes for hydro: the maps are TH2Ds in percent (a difference of two samples is no TProfile2D), a bin
// enters where the yield map is positive (hydro has no counting statistics to cut on), and the arrow of a tile is the yield-weighted mean of its
// bins, i.e. the d3N-weighted mean polarization of the Lambdas in the tile.
// --------------------------------------------------------------------------------------------------------------------------------------

/// \brief One arrow of a vector-field panel.
struct ArrowTile {
  double xCentre = 0.;
  double yCentre = 0.;
  double horizontal = 0.; // in-plane components, along the horizontal and vertical axes
  double vertical = 0.;
  double magnitude = 0.;
  double magnitudeError = 0.;
};

bool hasSmallerMagnitude(const ArrowTile& first, const ArrowTile& second) { return first.magnitude < second.magnitude; }

/// \brief Text in NDC of the current pad, centred at (x, y), owned by the pad.
void addLabel(double x, double y, const std::string& text, double size)
{
  auto* label = new TLatex(x, y, text.c_str());
  label->SetBit(kCanDelete);
  label->SetNDC();
  label->SetTextSize(size);
  label->SetTextAlign(22);
  label->Draw();
}

/// \brief A signed map as COLZ on the current pad, on a range symmetric about 0 set by the populated bins (positive yield).
void drawSymmetricMap(const TH2D* map, const TH2D* yield, const std::string& axes)
{
  auto* copy = static_cast<TH2D*>(map->Clone((std::string(map->GetName()) + "_drawn").c_str()));
  copy->SetDirectory(nullptr);
  copy->SetBit(kCanDelete);
  copy->SetStats(false);
  copy->SetTitle(axes.c_str());
  double maxAbs = 0.;
  for (int binX = 1; binX <= copy->GetNbinsX(); ++binX)
    for (int binY = 1; binY <= copy->GetNbinsY(); ++binY)
      if (yield->GetBinContent(binX, binY) > 0.)
        maxAbs = std::max(maxAbs, std::abs(copy->GetBinContent(binX, binY)));
  maxAbs = maxAbs > 0. ? 1.1 * maxAbs : 1.;
  // Bins without yield go below the minimum, which COLZ leaves unpainted: a 0 there would show the colour of the middle of the scale
  for (int binX = 1; binX <= copy->GetNbinsX(); ++binX)
    for (int binY = 1; binY <= copy->GetNbinsY(); ++binY)
      if (!(yield->GetBinContent(binX, binY) > 0.))
        copy->SetBinContent(binX, binY, -2. * maxAbs);
  copy->SetMinimum(-maxAbs);
  copy->SetMaximum(maxAbs);
  copy->GetZaxis()->SetTitleOffset(1.5);
  copy->Draw("COLZ");
}

/// \brief A yield map as COLZ on the current pad, from 0.
void drawYieldMap(const TH2D* yield, const std::string& axes)
{
  auto* copy = static_cast<TH2D*>(yield->Clone((std::string(yield->GetName()) + "_drawn").c_str()));
  copy->SetDirectory(nullptr);
  copy->SetBit(kCanDelete);
  copy->SetStats(false);
  copy->SetTitle(axes.c_str());
  copy->SetMinimum(0.);
  copy->GetZaxis()->SetTitleOffset(1.5);
  copy->Draw("COLZ");
}

/// \brief Vector-field panel on the current pad: the out-of-plane component as COLZ and arrows of the in-plane components, each the
///        yield-weighted mean over a tile of kArrowBlockSize x kArrowBlockSize bins, drawn at the true angle on screen.
void drawVectorField(const TH2D* horizontal, const TH2D* vertical, const TH2D* background, const TH2D* yield, const std::string& axes, const std::string& magnitudeLabel)
{
  drawSymmetricMap(background, yield, axes);
  gPad->Update(); // the pad geometry read below is valid only once the pad is laid out

  // The arrow direction lives in polarization space, where both components share their units, but its tip is placed in data coordinates.
  // Both axes rarely map the same data units per pixel, so the vertical displacement is converted through the pad geometry: the arrow keeps
  // the true angle on screen (the aspect correction of the original)
  const double frameWidth = gPad->GetWw() * gPad->GetAbsWNDC() * (1. - gPad->GetLeftMargin() - gPad->GetRightMargin());
  const double frameHeight = gPad->GetWh() * gPad->GetAbsHNDC() * (1. - gPad->GetTopMargin() - gPad->GetBottomMargin());
  const TAxis* xAxis = background->GetXaxis();
  const TAxis* yAxis = background->GetYaxis();
  const double xSpan = xAxis->GetXmax() - xAxis->GetXmin();
  const double ySpan = yAxis->GetXmax() - yAxis->GetXmin();
  double aspect = 1.;
  if (frameWidth > 0. && frameHeight > 0. && xSpan > 0. && ySpan > 0.)
    aspect = (ySpan / frameHeight) / (xSpan / frameWidth);

  // Tiles: yield-weighted means of the components, with their errors propagated as if the bins were independent (as in the original: the
  // bins share events, so the error of the scale reference is indicative only)
  std::vector<ArrowTile> tiles;
  const int nTilesX = background->GetNbinsX() / kArrowBlockSize;
  const int nTilesY = background->GetNbinsY() / kArrowBlockSize;
  for (int tileX = 0; tileX < nTilesX; ++tileX) {
    for (int tileY = 0; tileY < nTilesY; ++tileY) {
      double sumYield = 0., sumHorizontal = 0., sumVertical = 0., sumHorizontalErrorSq = 0., sumVerticalErrorSq = 0.;
      for (int stepX = 0; stepX < kArrowBlockSize; ++stepX) {
        for (int stepY = 0; stepY < kArrowBlockSize; ++stepY) {
          const int binX = tileX * kArrowBlockSize + stepX + 1;
          const int binY = tileY * kArrowBlockSize + stepY + 1;
          const double binYield = yield->GetBinContent(binX, binY);
          if (!(binYield > 0.))
            continue;
          sumYield += binYield;
          sumHorizontal += binYield * horizontal->GetBinContent(binX, binY);
          sumVertical += binYield * vertical->GetBinContent(binX, binY);
          sumHorizontalErrorSq += std::pow(binYield * horizontal->GetBinError(binX, binY), 2);
          sumVerticalErrorSq += std::pow(binYield * vertical->GetBinError(binX, binY), 2);
        }
      }
      if (!(sumYield > 0.))
        continue;
      ArrowTile tile;
      tile.xCentre = 0.5 * (xAxis->GetBinLowEdge(tileX * kArrowBlockSize + 1) + xAxis->GetBinUpEdge((tileX + 1) * kArrowBlockSize));
      tile.yCentre = 0.5 * (yAxis->GetBinLowEdge(tileY * kArrowBlockSize + 1) + yAxis->GetBinUpEdge((tileY + 1) * kArrowBlockSize));
      tile.horizontal = sumHorizontal / sumYield;
      tile.vertical = sumVertical / sumYield;
      const double horizontalError = std::sqrt(sumHorizontalErrorSq) / sumYield;
      const double verticalError = std::sqrt(sumVerticalErrorSq) / sumYield;
      tile.magnitude = std::hypot(tile.horizontal, tile.vertical);
      if (tile.magnitude > 0.)
        tile.magnitudeError = std::hypot(tile.horizontal * horizontalError, tile.vertical * verticalError) / tile.magnitude;
      tiles.push_back(tile);
    }
  }
  if (tiles.empty())
    return;

  // Scale reference: the kArrowScalePercentile quantile of the magnitudes (the first nonzero one above it if that is 0). Its arrow is 0.6
  // tile widths long, and longer arrows are capped there
  std::vector<ArrowTile> sortedTiles = tiles;
  std::sort(sortedTiles.begin(), sortedTiles.end(), hasSmallerMagnitude);
  std::size_t referenceIndex = kArrowScalePercentile * (sortedTiles.size() - 1);
  while (referenceIndex + 1 < sortedTiles.size() && !(sortedTiles[referenceIndex].magnitude > 0.))
    ++referenceIndex;
  const ArrowTile reference = sortedTiles[referenceIndex];
  if (!(reference.magnitude > 0.))
    return;
  const double lengthPerUnit = 0.6 * kArrowBlockSize * xAxis->GetBinWidth(1) / reference.magnitude;
  for (const ArrowTile& tile : tiles) {
    if (tile.magnitude < 0.05 * reference.magnitude) // near-zero arrows would be dots
      continue;
    const double length = std::min(tile.magnitude, reference.magnitude) * lengthPerUnit;
    auto* arrow = new TArrow(tile.xCentre, tile.yCentre, tile.xCentre + tile.horizontal / tile.magnitude * length, tile.yCentre + tile.vertical / tile.magnitude * length * aspect, 0.012, ">");
    arrow->SetBit(kCanDelete);
    arrow->SetLineColor(kBlack);
    arrow->SetFillColor(kBlack);
    arrow->SetLineWidth(2);
    arrow->Draw();
  }

  // Scale box, top left
  auto* box = new TPaveText(xAxis->GetXmin() + 0.05 * xSpan, yAxis->GetXmax() - 0.16 * ySpan, xAxis->GetXmin() + 0.65 * xSpan, yAxis->GetXmax() - 0.05 * ySpan, "arc");
  box->SetBit(kCanDelete);
  box->SetCornerRadius(0.15);
  box->SetFillColor(kWhite);
  box->SetFillStyle(1001);
  box->SetBorderSize(0);
  box->SetMargin(0.02);
  box->SetTextAlign(12);
  box->SetTextFont(63);
  box->SetTextSize(18);
  box->AddText(Form("%s_{%.0fpct} = (%.2f #pm %.2f)%%", magnitudeLabel.c_str(), 100. * kArrowScalePercentile, reference.magnitude, reference.magnitudeError));
  box->Draw();
}

enum PanelKind { kPanelVectorField = 0,
                 kPanelSymmetric,
                 kPanelYield };

/// \brief One panel of a map canvas.
struct MapPanel {
  PanelKind kind = kPanelVectorField;
  const TH2D* map = nullptr;        // the colour map: the out-of-plane component of a vector field, a ring scalar, or a yield
  const TH2D* horizontal = nullptr; // vector field: in-plane component along the horizontal axis
  const TH2D* vertical = nullptr;   // vector field: in-plane component along the vertical axis
  const TH2D* yield = nullptr;      // populated bins and tile weights (vector field and symmetric)
  std::string banner;               // short text above the panel
  std::string magnitudeLabel;       // vector field: symbol of the scale reference
};

/// \brief Axis titles of a map for its panel, with no main title (the banner and the canvas title replace it).
std::string panelAxes(const TH2D* map)
{
  return std::string(";") + map->GetXaxis()->GetTitle() + ";" + map->GetYaxis()->GetTitle() + ";" + map->GetZaxis()->GetTitle();
}

/// \brief Canvas of map panels side by side, 700 px each, with a banner per panel and the title across the top (the layout of
///        MakePanelCanvas() of auxiliaryPerConfigPlots.cxx), written into the directory.
void drawPanelCanvas(TDirectory* directory, const std::string& name, const std::string& title, const std::vector<MapPanel>& panels)
{
  directory->cd();
  const int nPanels = panels.size();
  TCanvas canvas(name.c_str(), title.c_str(), 700 * nPanels, 650);
  canvas.Divide(nPanels, 1, 0.005, 0.005);
  for (int index = 0; index < nPanels; ++index) {
    const MapPanel& panel = panels[index];
    canvas.cd(index + 1);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.16); // the palette and its title
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.10);
    if (panel.kind == kPanelVectorField)
      drawVectorField(panel.horizontal, panel.vertical, panel.map, panel.yield, panelAxes(panel.map), panel.magnitudeLabel);
    else if (panel.kind == kPanelSymmetric)
      drawSymmetricMap(panel.map, panel.yield, panelAxes(panel.map));
    else
      drawYieldMap(panel.map, panelAxes(panel.map));
    addLabel(0.5, 0.95, panel.banner, 0.045);
  }
  canvas.cd(0);
  addLabel(0.5, 0.98, title, 0.032);
  directory->WriteTObject(&canvas);
}

/// \brief Maps in percent of one producer TProfile2D for every sample, and the difference of the two pooled maps (errors added in
///        quadrature), stored in the directory: [sample], nullptr where absent.
std::array<TH2D*, kNCurveSamples> sampleMaps(TDirectory* directory, const std::array<TFile*, kNSamples>& files, const std::string& path, const std::string& name, const std::array<std::string, kNCurveSamples>& legends)
{
  std::array<TH2D*, kNCurveSamples> maps{};
  std::string producerTitle;
  for (int sample = 0; sample < kNSamples; ++sample) {
    if (files[sample] == nullptr)
      continue;
    const TProfile2D* profile = getObject<TProfile2D>(files[sample], path);
    producerTitle = profile->GetTitle();
    const std::string axes = std::string(";") + profile->GetXaxis()->GetTitle() + ";" + profile->GetYaxis()->GetTitle() + ";" + profile->GetZaxis()->GetTitle() + " (%)";
    maps[sample] = profileMap(directory, name + "_" + kSampleKeys[sample], producerTitle + ", " + legends[sample] + axes, profile, kPercent);
  }
  if (files[kSampleNoBullet] != nullptr) {
    directory->cd();
    maps[kSampleDifference] = static_cast<TH2D*>(maps[kSampleBullet]->Clone((name + "_" + kSampleKeys[kSampleDifference]).c_str()));
    maps[kSampleDifference]->SetDirectory(directory);
    maps[kSampleDifference]->Add(maps[kSampleNoBullet], -1.);
    maps[kSampleDifference]->SetTitle((producerTitle + ", " + legends[kSampleDifference]).c_str());
  }
  return maps;
}

/// \brief Integrals of R and R_z by the four estimators, true trigger against a variant set, one canvas per sample (Mixing/ and Smudge/).
void drawTrueVsVariant(const OutputFolder& folder, int set, const std::string& variantKey, const std::array<SetResult, kNTriggerSets>& results, const std::vector<int>& presentSamples, const std::array<std::string, kNCurveSamples>& sampleLegends, std::vector<std::string>& warnings)
{
  const std::vector<std::string> labels(kScalarLatex.begin(), kScalarLatex.begin() + kNRingIntegrals);
  for (const int sample : presentSamples) {
    const std::string key = kSampleKeys[sample];
    const EstimateSet& trueEstimates = results[kSetTrue].estimates[sample];
    const EstimateSet& variantEstimates = results[set].estimates[sample];
    const Part truePart = tableHistogram(folder.parts, "hIntegratedTrue_" + key, "Integrals of the true trigger, " + sampleLegends[sample] + ", jackknife errors", labels, trueEstimates.value, trueEstimates.error, blockIndices(results[kSetTrue].layout.scalar(0), kNRingIntegrals), 1., warnings);
    const Part variantPart = tableHistogram(folder.parts, "hIntegrated" + variantKey + "_" + key, "Integrals of the " + kSetLabels[set] + ", " + sampleLegends[sample] + ", jackknife errors", labels, variantEstimates.value, variantEstimates.error, blockIndices(results[set].layout.scalar(0), kNRingIntegrals), 1., warnings);
    drawOverlay(folder.canvases, "cIntegratedTrueVs" + variantKey + "_" + key, "Integrals, true against " + kSetLabels[set] + ", " + sampleLegends[sample] + ", jackknife errors", {{truePart, "true trigger", kPaletteLooks[0], kStyleMarker}, {variantPart, kSetLabels[set], kPaletteLooks[1], kStyleMarker}}, 0., kRingIntegralSeparators);
  }
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
    SummaryConfig config;
    std::vector<std::string> positionalArgs;
    std::map<std::string, std::string> options;
    for (int iArg = 1; iArg < argc; ++iArg) {
      const std::string argument = argv[iArg];
      if (argument == "--help" || argument == "-h") {
        printUsage();
        return 0;
      }
      if (argument.rfind("--", 0) == 0) {
        const std::size_t equalPos = argument.find('=');
        if (equalPos == std::string::npos)
          throw std::runtime_error("option without value: " + argument + " (use --option=value)");
        options[argument.substr(2, equalPos - 2)] = argument.substr(equalPos + 1);
      } else {
        positionalArgs.push_back(argument);
      }
    }
    if (positionalArgs.size() != 3) {
      printUsage();
      throw std::runtime_error("expected exactly three positional arguments, <bulletFile>, <noBulletFile|none> and <outputFile>");
    }
    config.bulletFile = positionalArgs[0];
    config.noBulletFile = positionalArgs[1];
    config.outputFile = positionalArgs[2];
    for (const auto& [key, value] : options) {
      if (key == "paired") {
        const int pairedFlag = toInt(key, value);
        if (pairedFlag != 0 && pairedFlag != 1)
          throw std::runtime_error("--paired must be 0 or 1, got '" + value + "'");
        config.paired = pairedFlag == 1;
      } else if (key == "alphaLambda") {
        config.alphaLambda = toDouble(key, value);
      } else if (key == "nLambdaPrecision") {
        config.nLambdaPrecision = toDouble(key, value);
      } else if (key == "centralityLabel") {
        config.centralityLabel = value;
      } else {
        throw std::runtime_error("unknown option --" + key + " (see --help)");
      }
    }
    if (config.alphaLambda == 0.)
      throw std::runtime_error("--alphaLambda must be nonzero");
    if (config.nLambdaPrecision < 0.)
      throw std::runtime_error("--nLambdaPrecision must be >= 0");
    const bool hasNoBullet = config.noBulletFile != "none";

    std::vector<std::pair<std::string, std::string>> configEntries;
    configEntries.emplace_back("bulletFile", config.bulletFile);
    configEntries.emplace_back("noBulletFile", config.noBulletFile);
    configEntries.emplace_back("outputFile", config.outputFile);
    configEntries.emplace_back("paired", config.paired ? "1" : "0");
    configEntries.emplace_back("alphaLambda", std::to_string(config.alphaLambda));
    configEntries.emplace_back("nLambdaPrecision", std::to_string(config.nLambdaPrecision));
    configEntries.emplace_back("centralityLabel", config.centralityLabel);
    std::cout << "========================================================" << std::endl;
    std::cout << "  hydroTruthLabSummary" << std::endl;
    for (const auto& [name, value] : configEntries)
      std::cout << "  " << name << " = " << value << std::endl;
    std::cout << "========================================================" << std::endl;

    std::vector<std::string> warnings; // printed at the end and stored in Config/: nothing is skipped silently

    // ==================================================================================================================================
    // 2. Producer files: configurations that must agree, the producers' own warnings, and the axes
    // ==================================================================================================================================
    std::unique_ptr<TFile> bulletFile(TFile::Open(config.bulletFile.c_str(), "READ"));
    if (!bulletFile || bulletFile->IsZombie())
      throw std::runtime_error("cannot open " + config.bulletFile);
    std::unique_ptr<TFile> noBulletFile;
    if (hasNoBullet) {
      noBulletFile.reset(TFile::Open(config.noBulletFile.c_str(), "READ"));
      if (!noBulletFile || noBulletFile->IsZombie())
        throw std::runtime_error("cannot open " + config.noBulletFile);
    }
    const std::array<TFile*, kNSamples> sampleFiles = {bulletFile.get(), noBulletFile.get()}; // nullptr without no_bullet
    const std::array<std::string, kNSamples> sampleLabels = {"bullet", "no_bullet"};

    const std::string triggerMode = readConfigValue(bulletFile.get(), "trigger");
    const int nDelta = std::stoi(readConfigValue(bulletFile.get(), "deltaScanSteps"));
    const int mixDepth = std::stoi(readConfigValue(bulletFile.get(), "mixDepth"));
    const double smudgeSigma = std::stod(readConfigValue(bulletFile.get(), "smudgeSigma"));
    if (hasNoBullet) {
      // Settings that define the observables must agree whatever the pairing (the producers print them in the same format)
      const std::vector<std::string> sharedKeys = {"trigger", "inputName", "ptMin", "ptMax", "absYMax", "lambdaMass", "deltaScanSteps", "mixDepth", "smudgeSigma", "smudgeMode", "smudgeDraws", "grid", "columnToPolarization", "mapSubdivisions"};
      // Paired trigger variants (README, key hypothesis): the smudged triggers need the same seed, the mixing partners the same event range
      std::vector<std::string> pairingKeys;
      if (config.paired && smudgeSigma > 0.)
        pairingKeys.push_back("seed");
      if (config.paired && mixDepth > 0) {
        pairingKeys.push_back("firstEvent");
        pairingKeys.push_back("nEvents");
      }
      for (const std::string& key : sharedKeys) {
        const std::string bulletValue = readConfigValue(bulletFile.get(), key);
        const std::string noBulletValue = readConfigValue(noBulletFile.get(), key);
        if (bulletValue != noBulletValue)
          throw std::runtime_error("Config/" + key + " differs between the samples ('" + bulletValue + "' vs '" + noBulletValue + "')");
      }
      for (const std::string& key : pairingKeys) {
        const std::string bulletValue = readConfigValue(bulletFile.get(), key);
        const std::string noBulletValue = readConfigValue(noBulletFile.get(), key);
        if (bulletValue != noBulletValue)
          throw std::runtime_error("Config/" + key + " differs between the samples ('" + bulletValue + "' vs '" + noBulletValue + "'), which breaks the pairing of the trigger variants: rerun both producers with the coordinator, or use --paired=0");
      }
    }
    // The producers' warnings, repeated so that the summary log shows them
    for (int sample = 0; sample < kNSamples; ++sample) {
      if (sampleFiles[sample] == nullptr)
        continue;
      const std::string producerWarnings = readConfigValue(sampleFiles[sample], "warnings");
      if (producerWarnings != "none")
        warnings.push_back("producer " + sampleLabels[sample] + ": " + producerWarnings);
    }

    // Delta phi axis, taken from the producer. The lobes, the fold and the asymmetry pair bin middle + j with bin middle - j, so the axis must
    // span 2 pi with an even number of bins and a bin centred at Delta phi = 0 (the producer's grid binning does)
    const TAxis* producerDeltaPhiAxis = getObject<TProfile>(bulletFile.get(), "Ring/pRingVsDeltaPhi")->GetXaxis();
    const AxisRange deltaPhiRange{producerDeltaPhiAxis->GetNbins(), producerDeltaPhiAxis->GetXmin(), producerDeltaPhiAxis->GetXmax(), producerDeltaPhiAxis->GetTitle()};
    const int nDeltaPhiBins = deltaPhiRange.nBins;
    const double deltaPhiWidth = (deltaPhiRange.high - deltaPhiRange.low) / nDeltaPhiBins;
    std::vector<double> deltaPhiCentres(nDeltaPhiBins);
    for (int bin = 0; bin < nDeltaPhiBins; ++bin)
      deltaPhiCentres[bin] = producerDeltaPhiAxis->GetBinCenter(bin + 1);
    if (nDeltaPhiBins % 2 != 0 || std::abs(deltaPhiCentres[nDeltaPhiBins / 2]) > kAxisTolerance * deltaPhiWidth || std::abs(deltaPhiRange.high - deltaPhiRange.low - 2. * kPi) > kAxisTolerance * deltaPhiWidth)
      throw std::runtime_error("the Delta phi axis of the producer does not span 2 pi with an even number of bins centred on Delta phi = 0");
    // Delta axis of the scan (and of mixing vs delta_eff), with the bin of delta = 0, the true trigger itself
    AxisRange deltaRange;
    int scanZeroBin = -1;
    if (nDelta > 0) {
      const TAxis* producerDeltaAxis = getObject<TProfile>(bulletFile.get(), "DeltaScan/pRingD3NVsDelta")->GetXaxis();
      deltaRange = {producerDeltaAxis->GetNbins(), producerDeltaAxis->GetXmin(), producerDeltaAxis->GetXmax(), producerDeltaAxis->GetTitle()};
      scanZeroBin = producerDeltaAxis->FindFixBin(0.) - 1;
      if (deltaRange.nBins != nDelta || std::abs(producerDeltaAxis->GetBinCenter(scanZeroBin + 1)) > kAxisTolerance * producerDeltaAxis->GetBinWidth(1))
        throw std::runtime_error("the delta axis of the producer does not have deltaScanSteps bins with one centred on delta = 0");
    }

    // ==================================================================================================================================
    // 3. Event trees
    // ==================================================================================================================================
    std::array<SampleData, kNSamples> samples;
    for (int sample = 0; sample < kNSamples; ++sample) {
      if (sampleFiles[sample] == nullptr)
        continue;
      samples[sample] = readSample(sampleFiles[sample], sampleLabels[sample], nDeltaPhiBins, nDelta, scanZeroBin, mixDepth > 0, smudgeSigma > 0.);
      std::cout << sampleLabels[sample] << ": " << samples[sample].eventId.size() << " events" << std::endl;
    }
    const int nBulletEvents = samples[kSampleBullet].eventId.size();
    const int nNoBulletEvents = samples[kSampleNoBullet].eventId.size();

    // ==================================================================================================================================
    // 4. Pairs: the event IDs both samples have, with the trigger of each pair checked (README, key hypothesis)
    // ==================================================================================================================================
    std::array<std::vector<int>, kNSamples> pairRows; // row of each pair in the tree of each sample, in pair order
    if (hasNoBullet) {
      std::map<int, int> noBulletRowOfEvent, bulletRowOfEvent;
      for (int row = 0; row < nNoBulletEvents; ++row)
        if (!noBulletRowOfEvent.emplace(samples[kSampleNoBullet].eventId[row], row).second)
          throw std::runtime_error("no_bullet: event " + std::to_string(samples[kSampleNoBullet].eventId[row]) + " appears twice in the tree");
      for (int row = 0; row < nBulletEvents; ++row) {
        const int eventId = samples[kSampleBullet].eventId[row];
        if (!bulletRowOfEvent.emplace(eventId, row).second)
          throw std::runtime_error("bullet: event " + std::to_string(eventId) + " appears twice in the tree");
        const auto partner = noBulletRowOfEvent.find(eventId);
        if (partner == noBulletRowOfEvent.end())
          continue;
        // Both producers read phi_t from the same text with the same parser (or set it to 0), so the triggers of a pair are bit-identical
        if (samples[kSampleBullet].phiTrigger[row] != samples[kSampleNoBullet].phiTrigger[partner->second])
          throw std::runtime_error("event " + std::to_string(eventId) + ": phiTrigger differs between bullet and no_bullet, so the samples are not paired");
        pairRows[kSampleBullet].push_back(row);
        pairRows[kSampleNoBullet].push_back(partner->second);
      }
    }
    const int nPairs = pairRows[kSampleBullet].size();
    if (hasNoBullet) {
      std::cout << nPairs << " pairs, " << nBulletEvents - nPairs << " bullet events without a partner, " << nNoBulletEvents - nPairs << " no_bullet events without a partner" << std::endl;
      if (config.paired && nPairs < 2)
        throw std::runtime_error("the paired jackknife needs at least two pairs");
      // Every bullet event has a no_bullet partner in this dataset. The extra no_bullet events (rp trigger) enter the no_bullet curves only
      if (config.paired && nBulletEvents > nPairs)
        warnings.push_back(std::to_string(nBulletEvents - nPairs) + " bullet events have no no_bullet partner: they enter the bullet curves but not the difference");
      if (config.paired && nNoBulletEvents > nPairs)
        std::cout << "[NOTE] " << nNoBulletEvents - nPairs << " no_bullet events have no bullet partner: they enter the no_bullet curves but not the difference" << std::endl;
    }

    // ==================================================================================================================================
    // 5. Estimates of every trigger set: each sample over all its events, and the difference over the pairs (all events with --paired=0).
    //    Then the leakage ansatz of every variant set, jackknifed jointly with the true trigger
    // ==================================================================================================================================
    std::vector<int> presentSamples = {kSampleBullet};
    if (hasNoBullet) {
      presentSamples.push_back(kSampleNoBullet);
      presentSamples.push_back(kSampleDifference);
    }
    std::array<SetResult, kNTriggerSets> results;
    for (int set = 0; set < kNTriggerSets; ++set) {
      if (!samples[kSampleBullet].sets[set].present) // the configurations agree, so the no_bullet file has the same sets
        continue;
      SetResult& result = results[set];
      result.present = true;
      result.layout.nBins = samples[kSampleBullet].sets[set].nBins;
      const bool withFeatures = set == kSetTrue;
      for (int sample = 0; sample < kNSamples; ++sample) {
        if (sampleFiles[sample] == nullptr)
          continue;
        result.stage1[sample] = replicateStage1(samples[sample].sets[set]);
        result.estimates[sample] = singleSampleEstimates(result.stage1[sample], withFeatures, result.layout, deltaPhiCentres, sampleLabels[sample] + ", " + kSetLabels[set], warnings);
      }
      if (!hasNoBullet)
        continue;
      if (config.paired) {
        for (int sample = 0; sample < kNSamples; ++sample)
          result.pairStage1[sample] = replicateStage1(selectEvents(samples[sample].sets[set], pairRows[sample]));
        result.differenceStage1 = differenceReplicas(result.pairStage1[kSampleBullet], result.pairStage1[kSampleNoBullet], true, result.layout);
      } else {
        result.differenceStage1 = differenceReplicas(result.stage1[kSampleBullet], result.stage1[kSampleNoBullet], false, result.layout);
      }
      result.estimates[kSampleDifference] = differenceEstimates(result.differenceStage1, withFeatures, result.layout, deltaPhiCentres, "bullet - no_bullet, " + kSetLabels[set], warnings);
    }
    const EstimateLayout& trueLayout = results[kSetTrue].layout;

    // Leakage: [set][sample], empty where not computed. The difference needs the per-sample replicas of the pairs, so --paired=1
    std::array<std::array<EstimateSet, kNCurveSamples>, kNTriggerSets> leakage;
    for (const int set : {kSetScan, kSetMixing, kSetSmudge}) {
      if (!results[set].present)
        continue;
      const EstimateLayout& variantLayout = results[set].layout;
      for (int sample = 0; sample < kNSamples; ++sample) {
        if (sampleFiles[sample] == nullptr)
          continue;
        const Replicated sampleLeakage = replicateLeakage(results[kSetTrue].stage1[sample], trueLayout, results[set].stage1[sample], variantLayout);
        leakage[set][sample] = estimatesFromReplicas(sampleLeakage.full, sampleLeakage.replicas, sampleLabels[sample] + ", leakage of the " + kSetLabels[set], warnings);
      }
      if (!hasNoBullet)
        continue;
      if (!config.paired) {
        warnings.push_back("RingOverlap/: --paired=0, so the leakage test of the bullet - no_bullet difference (" + kSetLabels[set] + ") is not computed");
        continue;
      }
      const Replicated bulletLeakage = replicateLeakage(results[kSetTrue].pairStage1[kSampleBullet], trueLayout, results[set].pairStage1[kSampleBullet], variantLayout);
      const Replicated noBulletLeakage = replicateLeakage(results[kSetTrue].pairStage1[kSampleNoBullet], trueLayout, results[set].pairStage1[kSampleNoBullet], variantLayout);
      std::vector<std::vector<double>> differenceReplicaList;
      for (std::size_t pair = 0; pair < bulletLeakage.replicas.size(); ++pair)
        differenceReplicaList.push_back(differenceOfVectors(bulletLeakage.replicas[pair], noBulletLeakage.replicas[pair]));
      leakage[set][kSampleDifference] = estimatesFromReplicas(differenceOfVectors(bulletLeakage.full, noBulletLeakage.full), differenceReplicaList, "bullet - no_bullet, leakage of the " + kSetLabels[set], warnings);
    }

    std::array<std::string, kNCurveSamples> sampleLegends;
    sampleLegends[kSampleBullet] = "bullet (" + std::to_string(nBulletEvents) + " evt.)";
    sampleLegends[kSampleNoBullet] = "no_bullet (" + std::to_string(nNoBulletEvents) + " evt.)";
    sampleLegends[kSampleDifference] = config.paired ? "bullet - no_bullet (" + std::to_string(nPairs) + " pairs)" : "bullet - no_bullet (means of all evt.)";

    // Integrated values of the true trigger in the log (it has no ring overlap)
    for (const int sample : presentSamples) {
      const EstimateSet& estimates = results[kSetTrue].estimates[sample];
      std::cout << "Integrated, true trigger, " << sampleLegends[sample] << ":" << std::endl;
      for (int scalar = 0; scalar < kScalarCosZetaD3N; ++scalar)
        std::cout << "  " << kScalarText[scalar] << " = " << estimates.value[trueLayout.scalar(scalar)] << " +- " << estimates.error[trueLayout.scalar(scalar)] << std::endl;
    }
    // Closure of the jet trigger: no_bullet carries no jet-induced ring, so its integrals must be consistent with zero
    std::string noBulletClosure = "not applicable (no no_bullet file, or the rp trigger)";
    if (hasNoBullet && triggerMode == "jet") {
      noBulletClosure.clear();
      for (const int scalar : {kScalarRingD3N, kScalarRingZD3N}) {
        const double value = results[kSetTrue].estimates[kSampleNoBullet].value[trueLayout.scalar(scalar)];
        const double error = results[kSetTrue].estimates[kSampleNoBullet].error[trueLayout.scalar(scalar)];
        noBulletClosure += kScalarText[scalar] + " = " + std::to_string(value) + " +- " + std::to_string(error) + " (" + std::to_string(error > 0. ? value / error : kNoValue) + " sigma). ";
      }
      std::cout << "no_bullet closure (jet trigger, expected consistent with 0): " << noBulletClosure << std::endl;
    }

    // ==================================================================================================================================
    // 6. Output file and folders: canvases at the top of each folder, the objects they draw in its Parts/
    // ==================================================================================================================================
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLegendBorderSize(0);
    std::unique_ptr<TFile> outputFile(TFile::Open(config.outputFile.c_str(), "RECREATE"));
    if (!outputFile || outputFile->IsZombie())
      throw std::runtime_error("cannot create " + config.outputFile);
    TDirectory* configDir = outputFile->mkdir("Config");
    const OutputFolder paper = makeFolder(outputFile.get(), "Paper");
    const OutputFolder integrated = makeFolder(outputFile.get(), "Integrated");
    const OutputFolder differentialRing = makeFolder(outputFile.get(), "DifferentialRing");
    const OutputFolder difference = hasNoBullet ? makeFolder(outputFile.get(), "Difference") : OutputFolder{};
    const OutputFolder scan = nDelta > 0 ? makeFolder(outputFile.get(), "DeltaScan") : OutputFolder{};
    const OutputFolder mixing = mixDepth > 0 ? makeFolder(outputFile.get(), "Mixing") : OutputFolder{};
    const OutputFolder smudge = smudgeSigma > 0. ? makeFolder(outputFile.get(), "Smudge") : OutputFolder{};
    const OutputFolder ringOverlap = nDelta > 0 || mixDepth > 0 || smudgeSigma > 0. ? makeFolder(outputFile.get(), "RingOverlap") : OutputFolder{};
    TDirectory* polarizationMapsDir = outputFile->mkdir("PolarizationMaps");
    const OutputFolder precision = makeFolder(outputFile.get(), "Precision");
    const OutputFolder jackknifeFolder = makeFolder(outputFile.get(), "Jackknife");
    const OutputFolder qa = makeFolder(outputFile.get(), "QA");
    // Folder of the Ring and RingZ curves of each set (the true trigger's live in DifferentialRing/, next to the other projections)
    const std::array<OutputFolder, kNTriggerSets> setFolders = {differentialRing, scan, mixing, smudge};

    // ==================================================================================================================================
    // 7. Ring and RingZ curves of every set, with jackknife errors: [set][curve][sample], drawn again by the later sections
    // ==================================================================================================================================
    std::array<std::array<std::array<Part, kNCurveSamples>, kNFeatureCurves>, kNTriggerSets> curveParts;
    for (int set = 0; set < kNTriggerSets; ++set) {
      if (!results[set].present)
        continue;
      const AxisRange& range = set == kSetScan ? deltaRange : deltaPhiRange;
      for (int curve = 0; curve < kNFeatureCurves; ++curve) {
        const std::string symbol = kFeatureCurveSymbols[curve];
        std::vector<DrawEntry> entries;
        for (const int sample : presentSamples) {
          const EstimateSet& estimates = results[set].estimates[sample];
          const std::string name = "h" + kFeatureCurveKeys[curve] + "Vs" + kSetAxisKeys[set] + "_" + kSampleKeys[sample];
          const std::string title = "<" + symbol + "> (" + kSetLabels[set] + ", " + sampleLegends[sample] + ", jackknife errors);" + range.title + ";<" + symbol + ">";
          curveParts[set][curve][sample] = curveHistogram(setFolders[set].parts, name, title, range, estimates.value, estimates.error, results[set].layout.curve(curve, 0), 1., warnings);
          entries.push_back({curveParts[set][curve][sample], sampleLegends[sample], kSampleLooks[sample], kStyleBand});
        }
        drawOverlay(setFolders[set].canvases, "c" + kFeatureCurveKeys[curve] + "Vs" + kSetAxisKeys[set], "<" + symbol + "> (" + kSetLabels[set] + ", pooled, jackknife bands)", entries, 0.);
      }
    }

    // ==================================================================================================================================
    // 8. Paper/: the figures of the paper for this dataset, in percent
    // ==================================================================================================================================
    // Fig. 2 (right): lab P_z vs (Delta phi, pT) with the jet
    const TProfile2D* fig2Profile = getObject<TProfile2D>(bulletFile.get(), "Paper/p2dPzLabVsDeltaPhiVsPt");
    TH2D* fig2Map = profileMap(paper.parts, "hPzLabVsDeltaPhiVsPt", "Fig. 2 (right): lab P_{z} (%), bullet, y cut only;" + std::string(fig2Profile->GetXaxis()->GetTitle()) + ";" + fig2Profile->GetYaxis()->GetTitle() + ";P_{z} (%)", fig2Profile, kPercent);
    drawMap(paper.canvases, "cFig2Right", fig2Map);
    // Fig. 3 (left): lab P_z vs (phi - Psi_RP, pT) without the jet
    if (hasNoBullet) {
      const TProfile2D* fig3Profile = getObject<TProfile2D>(noBulletFile.get(), "Paper/p2dPzLabVsPhiVsPt");
      TH2D* fig3Map = profileMap(paper.parts, "hPzLabVsPhiVsPt", "Fig. 3 (left): lab P_{z} (%), no_bullet, y cut only;" + std::string(fig3Profile->GetXaxis()->GetTitle()) + ";" + fig3Profile->GetYaxis()->GetTitle() + ";P_{z} (%)", fig3Profile, kPercent);
      drawMap(paper.canvases, "cFig3Left", fig3Map);
    } else {
      warnings.push_back("Paper/: no no_bullet file, so no Fig. 3 (left)");
    }
    // Fig. 4 analogue: R(Delta phi) of both samples and of the difference, with jackknife bands. The no_bullet curve uses all its events, the
    // difference only the pairs (the legends say how many)
    {
      std::vector<DrawEntry> fig4Entries;
      for (const int sample : presentSamples) {
        const EstimateSet& estimates = results[kSetTrue].estimates[sample];
        const Part part = curveHistogram(paper.parts, "hRingVsDeltaPhi_" + kSampleKeys[sample], "#it{R} (%), " + sampleLegends[sample] + ", pooled, jackknife errors;" + deltaPhiRange.title + ";#it{R} (%)", deltaPhiRange, estimates.value, estimates.error, trueLayout.curve(kCurveRing, 0), kPercent, warnings);
        fig4Entries.push_back({part, sampleLegends[sample], kSampleLooks[sample], kStyleBand});
      }
      drawOverlay(paper.canvases, "cFig4", "Fig. 4 analogue (" + config.centralityLabel + "): #it{R}(#Delta#varphi), pooled, jackknife bands", fig4Entries, 0.);
    }
    // Fig. 11 point: the integrated R of the bullet sample by three estimators, side by side in one centrality bin
    {
      const std::array<int, 3> fig11Scalars = {kScalarRingFlatEvtMean, kScalarRingD3NEvtMean, kScalarRingD3N};
      const std::array<std::string, 3> fig11Keys = {"FlatEvtMean", "D3NEvtMean", "D3N"};
      const std::array<std::string, 3> fig11Labels = {"flat, evt. mean #pm SEM (Vitor's estimator)", "d^{3}N weight, evt. mean #pm SEM", "d^{3}N weight, pooled #pm jackknife"};
      const EstimateSet& bulletEstimates = results[kSetTrue].estimates[kSampleBullet];
      std::vector<DrawEntry> fig11Entries;
      for (int estimator = 0; estimator < 3; ++estimator) {
        const Part part = tableHistogram(paper.parts, "hFig11" + fig11Keys[estimator], "Fig. 11 point, " + fig11Labels[estimator] + ";;#it{R} (%)", {config.centralityLabel}, bulletEstimates.value, bulletEstimates.error, {trueLayout.scalar(fig11Scalars[estimator])}, kPercent, warnings);
        fig11Entries.push_back({part, fig11Labels[estimator], kPaletteLooks[estimator], kStyleMarker});
      }
      drawOverlay(paper.canvases, "cFig11", "Fig. 11 point: integrated #it{R} of the bullet sample", fig11Entries, kNoValue);
    }

    // ==================================================================================================================================
    // 9. Integrated/: every integral of the true trigger, the lobe features of Ring and RingZ, and the comparison of their lobes
    // ==================================================================================================================================
    const std::vector<std::string> polarizationLabels(kScalarLatex.begin(), kScalarLatex.begin() + kNPolarizationScalars);
    const std::vector<std::string> geometryLabels(kScalarLatex.begin() + kScalarNzSqD3N, kScalarLatex.begin() + kScalarCosZetaD3N);
    const std::vector<std::string> lobePositionLabels(kLobeLabels.begin(), kLobeLabels.begin() + kNLobePositions);
    const std::vector<std::string> lobeHeightLabels(kLobeLabels.begin() + kNLobePositions, kLobeLabels.end());
    const std::vector<std::string> lobeShiftLabels(kLobeComparisonLabels.begin(), kLobeComparisonLabels.begin() + kNLobeShifts);
    const std::vector<std::string> lobeRatioLabels(kLobeComparisonLabels.begin() + kNLobeShifts, kLobeComparisonLabels.end());
    {
      std::vector<DrawEntry> integratedEntries, geometryEntries, shiftEntries, ratioEntries;
      std::array<std::vector<DrawEntry>, kNFeatureCurves> lobePositionEntries, lobeHeightEntries;
      for (const int sample : presentSamples) {
        const EstimateSet& estimates = results[kSetTrue].estimates[sample];
        const std::string key = kSampleKeys[sample];
        const Look& look = kSampleLooks[sample];
        const Part integratedPart = tableHistogram(integrated.parts, "hIntegrated_" + key, "Integrals in units of P, " + sampleLegends[sample] + ", jackknife errors", polarizationLabels, estimates.value, estimates.error, blockIndices(trueLayout.scalar(0), kNPolarizationScalars), 1., warnings);
        const Part geometryPart = tableHistogram(integrated.parts, "hIntegratedGeometry_" + key, "Dimensionless integrals, " + sampleLegends[sample] + ", jackknife errors", geometryLabels, estimates.value, estimates.error, blockIndices(trueLayout.scalar(kScalarNzSqD3N), kNGeometryScalars), 1., warnings);
        integratedEntries.push_back({integratedPart, sampleLegends[sample], look, kStyleMarker});
        geometryEntries.push_back({geometryPart, sampleLegends[sample], look, kStyleMarker});
        if (sample == kSampleDifference) {
          tableHistogram(integrated.parts, "hIntegrated_" + key + "_unpaired", "Integrals in units of P, " + sampleLegends[sample] + ", independent-samples errors", polarizationLabels, estimates.value, estimates.unpairedError, blockIndices(trueLayout.scalar(0), kNPolarizationScalars), 1., warnings);
          tableHistogram(integrated.parts, "hIntegratedGeometry_" + key + "_unpaired", "Dimensionless integrals, " + sampleLegends[sample] + ", independent-samples errors", geometryLabels, estimates.value, estimates.unpairedError, blockIndices(trueLayout.scalar(kScalarNzSqD3N), kNGeometryScalars), 1., warnings);
        }
        std::array<Part, kNFeatureCurves> positionParts;
        for (int curve = 0; curve < kNFeatureCurves; ++curve) {
          positionParts[curve] = tableHistogram(integrated.parts, "hLobePositions" + kFeatureCurveKeys[curve] + "_" + key, "Lobe positions of " + kFeatureCurveSymbols[curve] + "(#Delta#varphi), " + sampleLegends[sample] + ";;#Delta#varphi (rad)", lobePositionLabels, estimates.value, estimates.error, blockIndices(trueLayout.lobe(curve, 0), kNLobePositions), 1., warnings);
          const Part heights = tableHistogram(integrated.parts, "hLobeHeights" + kFeatureCurveKeys[curve] + "_" + key, "Lobe heights of " + kFeatureCurveSymbols[curve] + "(#Delta#varphi), " + sampleLegends[sample] + ";;" + kFeatureCurveSymbols[curve], lobeHeightLabels, estimates.value, estimates.error, blockIndices(trueLayout.lobe(curve, kNLobePositions), kNLobeFeatures - kNLobePositions), 1., warnings);
          lobePositionEntries[curve].push_back({positionParts[curve], sampleLegends[sample], look, kStyleMarker});
          lobeHeightEntries[curve].push_back({heights, sampleLegends[sample], look, kStyleMarker});
        }
        // R against R_z for this sample: where their lobes sit, and the jackknifed shifts and height ratios
        drawOverlay(integrated.canvases, "cLobePositionsRingVsRingZ_" + key, "Lobe positions of #it{R} and #it{R}_{z}, " + sampleLegends[sample] + ", jackknife errors", {{positionParts[kCurveRing], "#it{R}", kPaletteLooks[0], kStyleMarker}, {positionParts[kCurveRingZ], "#it{R}_{z}", kPaletteLooks[1], kStyleMarker}}, 0.);
        const Part shifts = tableHistogram(integrated.parts, "hLobeShifts_" + key, "Lobe positions, #it{R}_{z} - #it{R}, " + sampleLegends[sample] + ", jackknife errors;;#Delta#varphi shift (rad)", lobeShiftLabels, estimates.value, estimates.error, blockIndices(trueLayout.lobeComparison(0), kNLobeShifts), 1., warnings);
        const Part ratios = tableHistogram(integrated.parts, "hLobeHeightRatios_" + key, "Lobe heights, #it{R}_{z}/#it{R}, " + sampleLegends[sample] + ", jackknife errors;;#it{R}_{z}/#it{R}", lobeRatioLabels, estimates.value, estimates.error, blockIndices(trueLayout.lobeComparison(kNLobeShifts), kNLobeComparisons - kNLobeShifts), 1., warnings);
        shiftEntries.push_back({shifts, sampleLegends[sample], look, kStyleMarker});
        ratioEntries.push_back({ratios, sampleLegends[sample], look, kStyleMarker});
      }
      drawOverlay(integrated.canvases, "cIntegrated", "Integrals of the true trigger in units of P (d^{3}N-weighted and pooled unless labelled), jackknife errors", integratedEntries, 0., kRingIntegralSeparators);
      drawOverlay(integrated.canvases, "cIntegratedGeometry", "Dimensionless integrals of the true trigger, jackknife errors", geometryEntries, kNoValue);
      for (int curve = 0; curve < kNFeatureCurves; ++curve) {
        drawOverlay(integrated.canvases, "cLobePositions" + kFeatureCurveKeys[curve], "Lobe positions of " + kFeatureCurveSymbols[curve] + "(#Delta#varphi), jackknife errors", lobePositionEntries[curve], 0.);
        drawOverlay(integrated.canvases, "cLobeHeights" + kFeatureCurveKeys[curve], "Lobe heights of " + kFeatureCurveSymbols[curve] + "(#Delta#varphi), jackknife errors", lobeHeightEntries[curve], 0.);
      }
      drawOverlay(integrated.canvases, "cLobeShifts", "Lobe positions of #it{R}_{z} minus those of #it{R}, jackknife errors", shiftEntries, 0.);
      drawOverlay(integrated.canvases, "cLobeHeightRatios", "Lobe heights of #it{R}_{z} over those of #it{R}, jackknife errors", ratioEntries, kNoValue);
    }

    // ==================================================================================================================================
    // 10. DifferentialRing/: the projections of the true trigger vs Delta phi, <n_z^2>, and A_z = <R_z>/<n_z^2> against R (a diagnostic:
    //     README, "Observables")
    // ==================================================================================================================================
    {
      std::array<Part, kNCurveSamples> amplitudeZParts;
      std::vector<DrawEntry> nzSqEntries;
      for (const int sample : presentSamples) {
        const EstimateSet& estimates = results[kSetTrue].estimates[sample];
        const std::string key = kSampleKeys[sample];
        amplitudeZParts[sample] = curveHistogram(differentialRing.parts, "hAmplitudeZVsDeltaPhi_" + key, "#it{A}_{z} = <#it{R}_{z}>/<n_{z}^{2}> (" + sampleLegends[sample] + ", jackknife errors);" + deltaPhiRange.title + ";#it{A}_{z}", deltaPhiRange, estimates.value, estimates.error, trueLayout.curve(kCurveAmplitudeZ, 0), 1., warnings);
        std::vector<DrawEntry> projectionEntries = {{curveParts[kSetTrue][kCurveRing][sample], "#it{R} (jackknife)", kPaletteLooks[0], kStyleBand}, {curveParts[kSetTrue][kCurveRingZ][sample], "#it{R}_{z} (jackknife)", kPaletteLooks[1], kStyleBand}};
        if (sample != kSampleDifference) {
          const Part nzSq = curveHistogram(differentialRing.parts, "hNzSqVsDeltaPhi_" + key, "<n_{z}^{2}> (" + sampleLegends[sample] + ", jackknife errors);" + deltaPhiRange.title + ";<n_{z}^{2}>", deltaPhiRange, estimates.value, estimates.error, trueLayout.curve(kCurveNzSq, 0), 1., warnings);
          nzSqEntries.push_back({nzSq, sampleLegends[sample], kSampleLooks[sample], kStyleBand});
          // R_z^T, H and R_z^lab have no per-bin sums in the tree: the producer's pooled curves, with their pooled-master errors
          projectionEntries.push_back({profilePart(differentialRing.parts, sampleFiles[sample], "Ring/pRingZTVsDeltaPhi", "pRingZTVsDeltaPhi_" + key), "#it{R}_{z}^{T} (pooled-master errors)", kPaletteLooks[2], kStyleBand});
          projectionEntries.push_back({profilePart(differentialRing.parts, sampleFiles[sample], "Ring/pHelicityTermVsDeltaPhi", "pHelicityTermVsDeltaPhi_" + key), "#it{H} (pooled-master errors)", kPaletteLooks[3], kStyleBand});
          projectionEntries.push_back({profilePart(differentialRing.parts, sampleFiles[sample], "Ring/pRingZLabVsDeltaPhi", "pRingZLabVsDeltaPhi_" + key), "#it{R}_{z}^{lab} (pooled-master errors)", kPaletteLooks[4], kStyleBand});
        }
        drawOverlay(differentialRing.canvases, "cProjections_" + key, "Projections of the true trigger vs #Delta#varphi, " + sampleLegends[sample], projectionEntries, 0.);
      }
      drawOverlay(differentialRing.canvases, "cNzSqVsDeltaPhi", "<n_{z}^{2}> vs #Delta#varphi (true trigger), jackknife bands", nzSqEntries, kNoValue);
      // R and A_z of both samples, and their differences on a lower pad
      std::vector<DrawEntry> upperEntries, lowerEntries;
      for (int sample = 0; sample < kNSamples; ++sample) {
        if (sampleFiles[sample] == nullptr)
          continue;
        upperEntries.push_back({curveParts[kSetTrue][kCurveRing][sample], "#it{R}, " + sampleLabels[sample], kSampleLooks[sample], kStyleBand});
        upperEntries.push_back({amplitudeZParts[sample], "#it{A}_{z}, " + sampleLabels[sample], kAmplitudeZLooks[sample], kStyleBand});
      }
      const std::string amplitudeTitle = "#it{R} and #it{A}_{z} = <#it{R}_{z}>/<n_{z}^{2}> vs #Delta#varphi (true trigger), jackknife bands";
      if (hasNoBullet) {
        lowerEntries.push_back({curveParts[kSetTrue][kCurveRing][kSampleDifference], "#it{R}, " + sampleLegends[kSampleDifference], kSampleLooks[kSampleDifference], kStyleBand});
        lowerEntries.push_back({amplitudeZParts[kSampleDifference], "#it{A}_{z}, " + sampleLegends[kSampleDifference], kAmplitudeZLooks[kSampleDifference], kStyleBand});
        drawWithDifferencePad(differentialRing.canvases, "cAmplitudeZ", amplitudeTitle, upperEntries, lowerEntries, "bullet - no_bullet");
      } else {
        drawOverlay(differentialRing.canvases, "cAmplitudeZ", amplitudeTitle, upperEntries, 0.);
      }
    }

    // ==================================================================================================================================
    // 11. Difference/: paired against independent-samples errors, i.e. the gain of the key hypothesis
    // ==================================================================================================================================
    if (hasNoBullet && config.paired) {
      for (int set = 0; set < kNTriggerSets; ++set) {
        if (!results[set].present)
          continue;
        const EstimateSet& differenceEstimate = results[set].estimates[kSampleDifference];
        const AxisRange& range = set == kSetScan ? deltaRange : deltaPhiRange;
        // sigma_unpaired / sigma_paired for every component
        std::vector<double> gain(differenceEstimate.value.size(), kNoValue);
        for (std::size_t index = 0; index < gain.size(); ++index)
          if (differenceEstimate.error[index] != kNoValue && differenceEstimate.unpairedError[index] != kNoValue && differenceEstimate.error[index] > 0.)
            gain[index] = differenceEstimate.unpairedError[index] / differenceEstimate.error[index];
        const std::vector<double> noErrors(gain.size(), 0.);
        std::vector<DrawEntry> gainEntries;
        for (int curve = 0; curve < kNFeatureCurves; ++curve) {
          const std::string symbol = kFeatureCurveSymbols[curve];
          const std::string suffix = kFeatureCurveKeys[curve] + "Vs" + kSetAxisKeys[set] + "_" + kSetKeys[set];
          const Part unpaired = curveHistogram(difference.parts, "h" + suffix + "_unpaired", "<" + symbol + ">, bullet - no_bullet (" + kSetLabels[set] + ", independent-samples errors);" + range.title + ";<" + symbol + ">", range, differenceEstimate.value, differenceEstimate.unpairedError, results[set].layout.curve(curve, 0), 1., warnings);
          const Part gainPart = curveHistogram(difference.parts, "hPairingGain" + suffix, "#sigma_{unpaired}/#sigma_{paired} of " + symbol + " (" + kSetLabels[set] + ");" + range.title + ";#sigma_{unpaired}/#sigma_{paired}", range, gain, noErrors, results[set].layout.curve(curve, 0), 1., warnings);
          gainEntries.push_back({gainPart, symbol, kPaletteLooks[curve], kStyleLine});
          // One central line (the paired and unpaired values are the same numbers), the wide independent-samples band behind the paired one
          drawOverlay(difference.canvases, "c" + suffix, "<" + symbol + ">, bullet - no_bullet (" + kSetLabels[set] + "): paired against independent-samples errors", {{unpaired, "independent samples", kUnpairedBandLook, kStyleBandOnly}, {curveParts[set][curve][kSampleDifference], "paired (leave-one-pair-out)", kSampleLooks[kSampleDifference], kStyleBand}}, 0.);
        }
        drawOverlay(difference.canvases, "cPairingGain_" + kSetKeys[set], "Gain of the pairing, #sigma_{unpaired}/#sigma_{paired} (" + kSetLabels[set] + ")", gainEntries, 1.);
        if (set == kSetTrue) {
          const Part gainIntegrated = tableHistogram(difference.parts, "hPairingGainIntegrated", "Gain of the pairing for the integrals of the true trigger;;#sigma_{unpaired}/#sigma_{paired}", polarizationLabels, gain, noErrors, blockIndices(trueLayout.scalar(0), kNPolarizationScalars), 1., warnings);
          drawOverlay(difference.canvases, "cPairingGainIntegrated", "Gain of the pairing for the integrals of the true trigger", {{gainIntegrated, "#sigma_{unpaired}/#sigma_{paired}", kPaletteLooks[0], kStyleMarker}}, 1., kRingIntegralSeparators);
        }
      }
    } else if (hasNoBullet) {
      std::cout << "[NOTE] --paired=0: Difference/ holds no pairing gain, and every difference carries independent-samples errors" << std::endl;
    }

    // ==================================================================================================================================
    // 12. DeltaScan/: closures and the (Delta phi', delta) maps
    // ==================================================================================================================================
    std::string scanClosureSummary = "delta-scan off";
    if (nDelta > 0) {
      std::vector<std::string> closureLabels;
      std::vector<double> closureValues;
      scanClosureSummary.clear();
      for (int sample = 0; sample < kNSamples; ++sample) {
        if (sampleFiles[sample] == nullptr)
          continue;
        const SampleData& data = samples[sample];
        // delta = 0 against the true trigger, per event: float noise expected
        const std::array<std::string, kNBinArrays> integralNames = {"", "R", "R_z", "<n_z^2>", "<cos zeta> - 1", "<sin zeta>"};
        for (int array = kBinARing; array < kNBinArrays; ++array) {
          const double relative = data.scanZeroMaxValue[array] > 0. ? data.scanZeroMaxResidual[array] / data.scanZeroMaxValue[array] : data.scanZeroMaxResidual[array];
          closureLabels.push_back("delta = 0 vs true, " + integralNames[array] + ", " + sampleLabels[sample]);
          closureValues.push_back(relative);
        }
        // Half-turn on the pooled curve, max_k |R(delta_k + pi) + R(delta_k)| / max_k |R(delta_k)|, and the uniform average over delta. Both
        // vanish up to float noise, since R(delta + pi) = -R(delta) holds event by event for transverse triggers
        const EstimateSet& scanEstimates = results[kSetScan].estimates[sample];
        double maxResidual = 0.;
        double maxAbsValue = 0.;
        for (int step = 0; step < nDelta / 2; ++step) {
          const double value = scanEstimates.value[results[kSetScan].layout.curve(kCurveRing, step)];
          const double opposite = scanEstimates.value[results[kSetScan].layout.curve(kCurveRing, step + nDelta / 2)];
          maxResidual = std::max(maxResidual, std::abs(value + opposite));
          maxAbsValue = std::max({maxAbsValue, std::abs(value), std::abs(opposite)});
        }
        const double uniformAverage = scanEstimates.value[results[kSetScan].layout.scalar(kScalarRingFlat)];
        closureLabels.push_back("half-turn, pooled R, " + sampleLabels[sample]);
        closureValues.push_back(maxAbsValue > 0. ? maxResidual / maxAbsValue : maxResidual);
        closureLabels.push_back("|average over delta| / max |R|, " + sampleLabels[sample]);
        closureValues.push_back(maxAbsValue > 0. ? std::abs(uniformAverage) / maxAbsValue : std::abs(uniformAverage));
      }
      for (std::size_t check = 0; check < closureLabels.size(); ++check) {
        std::cout << "DeltaScan closure, " << closureLabels[check] << ": " << closureValues[check] << " (relative, expected at float noise)" << std::endl;
        scanClosureSummary += closureLabels[check] + ": " + std::to_string(closureValues[check]) + "\n";
      }
      const std::vector<double> noErrors(closureValues.size(), 0.);
      const Part closureTable = tableHistogram(scan.parts, "hClosure", "Closures of the delta-scan (relative residuals, float noise expected)", closureLabels, closureValues, noErrors, blockIndices(0, closureLabels.size()), 1., warnings);
      drawOverlay(scan.canvases, "cClosure", "Closures of the delta-scan (relative residuals, float noise expected)", {{closureTable, "relative residual", kPaletteLooks[0], kStyleMarker}}, 0.);
      // Per-event half-turn closure of the producers, summed over the samples
      scan.parts->cd();
      auto* halfTurn = static_cast<TH1D*>(getObject<TH1D>(bulletFile.get(), "DeltaScan/hClosureHalfTurn")->Clone("hClosureHalfTurn"));
      halfTurn->SetDirectory(scan.parts);
      if (hasNoBullet)
        halfTurn->Add(getObject<TH1D>(noBulletFile.get(), "DeltaScan/hClosureHalfTurn"));
      halfTurn->SetTitle("Per-event half-turn closure, both samples;log_{10}(relative residual);Events");
      drawOverlay(scan.canvases, "cClosureHalfTurn", "Per-event half-turn closure, both samples", {{histogramPart(halfTurn), "events", kPaletteLooks[0], kStyleHistogram}}, kNoValue);
      // (Delta phi', delta) maps of the producers (pooled: no per-event sums of these)
      for (int sample = 0; sample < kNSamples; ++sample) {
        if (sampleFiles[sample] == nullptr)
          continue;
        for (int curve = 0; curve < kNFeatureCurves; ++curve) {
          const TProfile2D* scanProfile = getObject<TProfile2D>(sampleFiles[sample], "DeltaScan/p2d" + kFeatureCurveKeys[curve] + "VsDeltaPhiVsDelta");
          TH2D* scanMap = profileMap(scan.parts, "h" + kFeatureCurveKeys[curve] + "VsDeltaPhiPrimeVsDelta_" + kSampleKeys[sample], "<" + kFeatureCurveSymbols[curve] + "> vs (#Delta#varphi', #delta), " + sampleLegends[sample] + ";" + scanProfile->GetXaxis()->GetTitle() + ";" + scanProfile->GetYaxis()->GetTitle() + ";<" + kFeatureCurveSymbols[curve] + ">", scanProfile, 1.);
          drawMap(scan.canvases, "c" + kFeatureCurveKeys[curve] + "VsDeltaPhiPrimeVsDelta_" + kSampleKeys[sample], scanMap);
        }
      }
    }

    // ==================================================================================================================================
    // 13. Mixing/: true against mixed triggers, and the pulls of mixing vs delta_eff against the delta-scan vs delta
    // ==================================================================================================================================
    std::string pullSummary = "mixing off";
    if (mixDepth > 0) {
      for (const int sample : presentSamples)
        for (int curve = 0; curve < kNFeatureCurves; ++curve)
          drawOverlay(mixing.canvases, "c" + kFeatureCurveKeys[curve] + "TrueVsMixed_" + kSampleKeys[sample], "<" + kFeatureCurveSymbols[curve] + "> vs #Delta#varphi, true against mixed triggers, " + sampleLegends[sample], {{curveParts[kSetTrue][curve][sample], "true trigger", kPaletteLooks[0], kStyleBand}, {curveParts[kSetMixing][curve][sample], "mixed triggers", kPaletteLooks[1], kStyleBand}}, 0.);
      drawTrueVsVariant(mixing, kSetMixing, "Mixed", results, presentSamples, sampleLegends, warnings);
      if (nDelta > 0) {
        // The mixed integrals in one delta_eff bin come from a subsample of the events the scan uses at that delta. With v_e the per-event scan
        // value, s^2 its d3N-weighted spread around the pooled scan value, W_e the scan weights (total T) and c_e the summed weights of the
        // mixed entries of event e in the bin (total C), Var(mix - scan) = s^2 sum_e (c_e/C - W_e/T)^2. The sum accounts for the events the
        // two share, and s comes from all events, not from the few in the mixing bin. Events, not entries, are the units: two partners of one
        // event in the same bin carry nearly the same value (README, "Delta-scan, mixing and smudging")
        pullSummary.clear();
        std::vector<std::string> pullLabels;
        std::vector<double> pullMeans, pullRmsValues;
        std::vector<DrawEntry> pullVsDeltaEntries, pullEntries;
        int nSkippedBins = 0;
        for (int sample = 0; sample < kNSamples; ++sample) {
          if (sampleFiles[sample] == nullptr)
            continue;
          const SampleData& data = samples[sample];
          const SetSums& scanSums = data.sets[kSetScan];
          const EstimateSet& scanEstimates = results[kSetScan].estimates[sample];
          // c_e per delta_eff bin, event-major, with the partners rebuilt by the producer's rule: event i takes the phi_J of the events
          // i + 1, ..., i + mixDepth (mod nEvents), and phi_J is the phiTrigger of the jet trigger. Each mixed entry weighs the in-cut yield
          const int firstEvent = std::stoi(readConfigValue(sampleFiles[sample], "firstEvent"));
          const int nConfiguredEvents = std::stoi(readConfigValue(sampleFiles[sample], "nEvents"));
          std::map<int, int> rowOfEvent;
          for (int row = 0; row < scanSums.nEvents; ++row)
            rowOfEvent.emplace(data.eventId[row], row);
          const TAxis* deltaEffAxis = getObject<TProfile>(sampleFiles[sample], "Mixing/pRingD3NVsDeltaEff")->GetXaxis();
          std::vector<double> mixedWeight(scanSums.nEvents * nDelta, 0.);
          std::string rebuildProblem;
          for (int row = 0; row < scanSums.nEvents && rebuildProblem.empty(); ++row) {
            for (int partner = 1; partner <= mixDepth; ++partner) {
              const int partnerId = firstEvent + (data.eventId[row] - firstEvent + partner) % nConfiguredEvents;
              const auto partnerRow = rowOfEvent.find(partnerId);
              if (partnerRow == rowOfEvent.end()) {
                rebuildProblem = "partner event " + std::to_string(partnerId) + " is not in the event tree";
                break;
              }
              // The producer's wrap and the profile's own binning, on the same numbers: the bins agree exactly. An entry the rounding puts
              // beyond the axis went to the overflow of the producer's profile as well
              const double angle = data.phiTrigger[partnerRow->second] - data.phiTrigger[row];
              const double deltaEff = angle - TMath::TwoPi() * std::floor((angle - deltaRange.low) / TMath::TwoPi());
              const int bin = deltaEffAxis->FindFixBin(deltaEff) - 1;
              if (bin >= 0 && bin < nDelta)
                mixedWeight[row * nDelta + bin] += data.sumWInCuts[row];
            }
          }
          for (int curve = 0; curve < kNFeatureCurves; ++curve) {
            const std::string suffix = kFeatureCurveKeys[curve] + "_" + kSampleKeys[sample];
            const Part mixPart = profilePart(mixing.parts, sampleFiles[sample], "Mixing/p" + kFeatureCurveKeys[curve] + "D3NVsDeltaEff", "p" + kFeatureCurveKeys[curve] + "D3NVsDeltaEff_" + kSampleKeys[sample]);
            const auto* mixProfile = static_cast<const TProfile*>(mixPart.hist);
            if (mixProfile->GetNbinsX() != nDelta)
              throw std::runtime_error("Mixing/p" + kFeatureCurveKeys[curve] + "D3NVsDeltaEff does not use the delta-scan binning");
            // The two curves the pulls compare
            drawOverlay(mixing.canvases, "c" + kFeatureCurveKeys[curve] + "D3NMixedVsScan_" + kSampleKeys[sample], "Integrated <" + kFeatureCurveSymbols[curve] + ">: mixing vs #delta_{eff} against the delta-scan vs #delta, " + sampleLegends[sample], {{curveParts[kSetScan][curve][sample], "delta-scan vs #delta (jackknife)", kPaletteLooks[0], kStyleBand}, {mixPart, "mixed vs #delta_{eff} (pooled-master errors)", kPaletteLooks[1], kStyleMarker}}, 0.);
            // The rebuilt entries must carry the weights the producer filled, bin by bin
            for (int step = 0; step < nDelta && rebuildProblem.empty(); ++step) {
              double rebuiltWeight = 0.;
              for (int iEvent = 0; iEvent < scanSums.nEvents; ++iEvent)
                rebuiltWeight += mixedWeight[iEvent * nDelta + step];
              const double filledWeight = mixProfile->GetBinEntries(step + 1);
              if (std::abs(rebuiltWeight - filledWeight) > kRebuildTolerance * std::abs(filledWeight))
                rebuildProblem = "the rebuilt mixed entries do not match the weights of " + std::string(mixProfile->GetName()) + " in delta bin " + std::to_string(step);
            }
            if (!rebuildProblem.empty()) {
              warnings.push_back("Mixing/, " + sampleLabels[sample] + ", " + kFeatureCurveKeys[curve] + ": " + rebuildProblem + ", so no pulls (has the partner rule of the producer changed?)");
              continue;
            }
            const int array = curve == kCurveRing ? kBinARing : kBinARingZ;
            mixing.parts->cd();
            auto* pullVsDelta = new TH1D(("hPullVsDelta" + suffix).c_str(), ("Pull of mixed vs #delta_{eff} against the scan vs #delta, " + kFeatureCurveSymbols[curve] + ", " + sampleLegends[sample] + ";" + deltaRange.title + ";pull").c_str(), deltaRange.nBins, deltaRange.low, deltaRange.high);
            auto* pullDistribution = new TH1D(("hPull" + suffix).c_str(), ("Pulls of mixed against scan, " + kFeatureCurveSymbols[curve] + ", " + sampleLegends[sample] + ";pull;#delta bins").c_str(), kPullBins, -kPullRange, kPullRange);
            Part pullVsDeltaPart = histogramPart(pullVsDelta);
            std::vector<double> pulls;
            for (int step = 0; step < nDelta; ++step) {
              pullVsDeltaPart.defined[step] = false;
              const double scanValue = scanEstimates.value[results[kSetScan].layout.curve(curve, step)];
              double sumW = 0., sumMixedWeight = 0., sumWDeviationSq = 0.;
              for (int iEvent = 0; iEvent < scanSums.nEvents; ++iEvent) {
                sumW += scanSums.perBin[kBinW][iEvent * nDelta + step];
                sumMixedWeight += mixedWeight[iEvent * nDelta + step];
              }
              if (scanValue == kNoValue || !(sumMixedWeight > 0.) || !(sumW > 0.)) {
                ++nSkippedBins;
                continue;
              }
              double overlapFactor = 0.;
              for (int iEvent = 0; iEvent < scanSums.nEvents; ++iEvent) {
                const double weight = scanSums.perBin[kBinW][iEvent * nDelta + step];
                overlapFactor += std::pow(mixedWeight[iEvent * nDelta + step] / sumMixedWeight - weight / sumW, 2);
                if (weight > 0.)
                  sumWDeviationSq += weight * std::pow(scanSums.perBin[array][iEvent * nDelta + step] / weight - scanValue, 2);
              }
              const double spread = std::sqrt(sumWDeviationSq / sumW);
              if (!(spread > 0.) || !(overlapFactor > 0.)) {
                ++nSkippedBins;
                continue;
              }
              const double pull = (mixProfile->GetBinContent(step + 1) - scanValue) / (spread * std::sqrt(overlapFactor));
              pullVsDelta->SetBinContent(step + 1, pull);
              pullVsDeltaPart.defined[step] = true;
              pullDistribution->Fill(pull);
              pulls.push_back(pull);
            }
            // Mean and RMS from the pulls themselves: the histogram statistics would drop any pull beyond its range
            double pullMean = kNoValue;
            double pullRms = kNoValue;
            if (!pulls.empty()) {
              pullMean = std::accumulate(pulls.begin(), pulls.end(), 0.) / pulls.size();
              double sumSquares = 0.;
              for (const double pull : pulls)
                sumSquares += (pull - pullMean) * (pull - pullMean);
              pullRms = std::sqrt(sumSquares / pulls.size());
            }
            pullLabels.push_back(kFeatureCurveKeys[curve] + ", " + sampleLabels[sample]);
            pullMeans.push_back(pullMean);
            pullRmsValues.push_back(pullRms);
            pullSummary += kFeatureCurveKeys[curve] + ", " + sampleLabels[sample] + ": mean " + std::to_string(pullMean) + ", RMS " + std::to_string(pullRms) + " (" + std::to_string(pulls.size()) + " bins)\n";
            const Look& look = kPaletteLooks[pullEntries.size()];
            pullVsDeltaEntries.push_back({pullVsDeltaPart, kFeatureCurveSymbols[curve] + ", " + sampleLabels[sample], look, kStyleLine});
            pullEntries.push_back({histogramPart(pullDistribution), kFeatureCurveSymbols[curve] + ", " + sampleLabels[sample], look, kStyleHistogram});
          }
        }
        std::cout << "Mixing pulls against the delta-scan (mean ~0, RMS ~1 expected):\n" << pullSummary;
        if (nSkippedBins > 0)
          warnings.push_back("Mixing/: " + std::to_string(nSkippedBins) + " delta bins left out of the pulls (no mixed entry, or no spread of the scan)");
        const std::vector<double> noErrors(pullLabels.size(), 0.);
        const Part meanTable = tableHistogram(mixing.parts, "hPullMean", "Mean of the pulls of mixed against scan (0 expected);;mean", pullLabels, pullMeans, noErrors, blockIndices(0, pullLabels.size()), 1., warnings);
        const Part rmsTable = tableHistogram(mixing.parts, "hPullRms", "RMS of the pulls of mixed against scan (1 expected);;RMS", pullLabels, pullRmsValues, noErrors, blockIndices(0, pullLabels.size()), 1., warnings);
        drawOverlay(mixing.canvases, "cPullMean", "Mean of the pulls of mixed against scan (0 expected)", {{meanTable, "mean", kPaletteLooks[0], kStyleMarker}}, 0.);
        drawOverlay(mixing.canvases, "cPullRms", "RMS of the pulls of mixed against scan (1 expected)", {{rmsTable, "RMS", kPaletteLooks[0], kStyleMarker}}, 1.);
        drawOverlay(mixing.canvases, "cPullVsDelta", "Pulls of mixing vs #delta_{eff} against the delta-scan vs #delta", pullVsDeltaEntries, 0.);
        drawOverlay(mixing.canvases, "cPulls", "Pulls of mixing vs #delta_{eff} against the delta-scan vs #delta", pullEntries, kNoValue);
      } else {
        warnings.push_back("Mixing/: the delta-scan is off, so mixing is not compared with it");
      }
    }

    // ==================================================================================================================================
    // 14. Smudge/: true against smudged triggers, the smudged/true ratio of the integrals, and the integrals vs smudge angle
    // ==================================================================================================================================
    if (smudgeSigma > 0.) {
      for (const int sample : presentSamples)
        for (int curve = 0; curve < kNFeatureCurves; ++curve)
          drawOverlay(smudge.canvases, "c" + kFeatureCurveKeys[curve] + "TrueVsSmudged_" + kSampleKeys[sample], "<" + kFeatureCurveSymbols[curve] + "> vs #Delta#varphi, true against smudged triggers, " + sampleLegends[sample], {{curveParts[kSetTrue][curve][sample], "true trigger", kPaletteLooks[0], kStyleBand}, {curveParts[kSetSmudge][curve][sample], "smudged triggers", kPaletteLooks[1], kStyleBand}}, 0.);
      drawTrueVsVariant(smudge, kSetSmudge, "Smudged", results, presentSamples, sampleLegends, warnings);
      // Smudged over true, replica by replica over the same events. For the difference, the ratio of the paired differences: without pairs its
      // replicas would come in two independent groups for numerator and denominator alike, which this ratio does not handle
      if (hasNoBullet && !config.paired)
        warnings.push_back("Smudge/: --paired=0, so the smudged/true ratio of the difference is not computed");
      const std::vector<std::string> ratioLabels = {kScalarLatex[kScalarRingD3N], kScalarLatex[kScalarRingZD3N]};
      std::vector<DrawEntry> ratioEntries;
      for (const int sample : presentSamples) {
        if (sample == kSampleDifference && !config.paired)
          continue;
        std::vector<double> ratioValues, ratioErrors;
        for (const int scalar : {kScalarRingD3N, kScalarRingZD3N}) {
          EstimatePoint ratio;
          if (sample != kSampleDifference)
            ratio = crossSetRatio(results[kSetSmudge].stage1[sample].full, results[kSetSmudge].stage1[sample].replicas, results[kSetTrue].stage1[sample].full, results[kSetTrue].stage1[sample].replicas, trueLayout.scalar(scalar));
          else
            ratio = crossSetRatio(results[kSetSmudge].differenceStage1.full, results[kSetSmudge].differenceStage1.paired, results[kSetTrue].differenceStage1.full, results[kSetTrue].differenceStage1.paired, trueLayout.scalar(scalar));
          ratioValues.push_back(ratio.value);
          ratioErrors.push_back(ratio.error);
        }
        const Part ratioTable = tableHistogram(smudge.parts, "hSmudgedOverTrue_" + kSampleKeys[sample], "Integrals, smudged over true trigger, " + sampleLegends[sample] + ", jackknife errors;;smudged/true", ratioLabels, ratioValues, ratioErrors, blockIndices(0, 2), 1., warnings);
        ratioEntries.push_back({ratioTable, sampleLegends[sample], kSampleLooks[sample], kStyleMarker});
      }
      drawOverlay(smudge.canvases, "cSmudgedOverTrue", "Integrals, smudged over true trigger, jackknife errors", ratioEntries, 1.);
      // Integrals vs smudge angle: per draw in the producer, so only pooled (pooled-master errors)
      for (int curve = 0; curve < kNFeatureCurves; ++curve) {
        std::vector<DrawEntry> angleEntries;
        for (int sample = 0; sample < kNSamples; ++sample) {
          if (sampleFiles[sample] == nullptr)
            continue;
          const Part anglePart = profilePart(smudge.parts, sampleFiles[sample], "Smudge/p" + kFeatureCurveKeys[curve] + "D3NVsSmudgeAngle", "p" + kFeatureCurveKeys[curve] + "D3NVsSmudgeAngle_" + kSampleKeys[sample]);
          angleEntries.push_back({anglePart, sampleLegends[sample] + " (pooled-master errors)", kSampleLooks[sample], kStyleBand});
        }
        drawOverlay(smudge.canvases, "c" + kFeatureCurveKeys[curve] + "D3NVsSmudgeAngle", "Integrated <" + kFeatureCurveSymbols[curve] + "> vs smudge angle", angleEntries, 0.);
      }
    }

    // ==================================================================================================================================
    // 15. RingOverlap/: <cos zeta> and <sin zeta> of every variant set, their distributions, and the leakage ansatz (README, "Ring overlap")
    // ==================================================================================================================================
    std::string leakageSummary = "no variant trigger";
    if (ringOverlap.canvases != nullptr) {
      leakageSummary.clear();
      const std::array<int, 2> overlapCurves = {kCurveCosZeta, kCurveSinZeta};
      const std::array<int, 2> overlapScalars = {kScalarCosZetaD3N, kScalarSinZetaD3N};
      const std::array<std::string, 2> overlapKeys = {"CosZeta", "SinZeta"};
      const std::array<std::string, 2> overlapSymbols = {"cos#zeta", "sin#zeta"};
      const std::array<std::string, 2> overlapProducerHists = {"hCosZeta", "hSinZeta"};
      std::vector<std::string> integratedLabels;
      std::array<std::vector<double>, kNSamples> integratedValues, integratedErrors;
      for (const int set : {kSetScan, kSetMixing, kSetSmudge}) {
        if (!results[set].present)
          continue;
        const AxisRange& range = set == kSetScan ? deltaRange : deltaPhiRange;
        for (int overlap = 0; overlap < 2; ++overlap) {
          // Curves vs Delta phi' (vs delta for the scan), bullet and no_bullet: the overlap is geometry, the samples differ by their weights only
          std::vector<DrawEntry> curveEntries, distributionEntries;
          for (int sample = 0; sample < kNSamples; ++sample) {
            if (sampleFiles[sample] == nullptr)
              continue;
            const EstimateSet& estimates = results[set].estimates[sample];
            const Part part = curveHistogram(ringOverlap.parts, "h" + overlapKeys[overlap] + "Vs" + kSetAxisKeys[set] + "_" + kSetKeys[set] + "_" + kSampleKeys[sample], "<" + overlapSymbols[overlap] + "> (" + kSetLabels[set] + ", " + sampleLegends[sample] + ", jackknife errors);" + range.title + ";<" + overlapSymbols[overlap] + ">", range, estimates.value, estimates.error, results[set].layout.curve(overlapCurves[overlap], 0), 1., warnings);
            curveEntries.push_back({part, sampleLegends[sample], kSampleLooks[sample], kStyleBand});
            // Distributions over the cells (d3N-weighted, all events), normalized to unit area: the weights have no absolute scale
            if (set == kSetScan) {
              const TH2D* producerMap = getObject<TH2D>(sampleFiles[sample], "DeltaScan/h2d" + overlapKeys[overlap] + "VsDelta");
              ringOverlap.parts->cd();
              auto* map = static_cast<TH2D*>(producerMap->Clone(("h2d" + overlapKeys[overlap] + "VsDelta_" + kSampleKeys[sample]).c_str()));
              map->SetDirectory(ringOverlap.parts);
              map->SetTitle((std::string(producerMap->GetTitle()) + ", " + sampleLegends[sample]).c_str());
              drawMap(ringOverlap.canvases, "c" + overlapKeys[overlap] + "VsDeltaMap_" + kSampleKeys[sample], map);
            } else {
              const TH1D* producerHist = getObject<TH1D>(sampleFiles[sample], kSetProducerFolder[set] + "/" + overlapProducerHists[overlap]);
              ringOverlap.parts->cd();
              auto* distribution = static_cast<TH1D*>(producerHist->Clone((overlapProducerHists[overlap] + "_" + kSetKeys[set] + "_" + kSampleKeys[sample]).c_str()));
              distribution->SetDirectory(ringOverlap.parts);
              if (distribution->Integral() > 0.)
                distribution->Scale(1. / distribution->Integral());
              distribution->SetTitle((overlapSymbols[overlap] + " of the " + kSetLabels[set] + " (d^{3}N-weighted, normalized), " + sampleLegends[sample] + ";" + overlapSymbols[overlap] + ";fraction").c_str());
              distributionEntries.push_back({histogramPart(distribution), sampleLegends[sample], kSampleLooks[sample], kStyleHistogram});
            }
          }
          drawOverlay(ringOverlap.canvases, "c" + overlapKeys[overlap] + "Vs" + kSetAxisKeys[set] + "_" + kSetKeys[set], "<" + overlapSymbols[overlap] + "> (" + kSetLabels[set] + "), jackknife bands", curveEntries, 0.);
          drawOverlay(ringOverlap.canvases, "c" + overlapKeys[overlap] + "Distribution_" + kSetKeys[set], overlapSymbols[overlap] + " of the " + kSetLabels[set] + " (d^{3}N-weighted, normalized)", distributionEntries, kNoValue);
          // Integrals. The scan's average over delta vanishes by the half-turn identity and is not shown
          if (set == kSetScan)
            continue;
          integratedLabels.push_back("<" + overlapSymbols[overlap] + ">, " + kSetLabels[set]);
          for (int sample = 0; sample < kNSamples; ++sample) {
            if (sampleFiles[sample] == nullptr)
              continue;
            integratedValues[sample].push_back(results[set].estimates[sample].value[results[set].layout.scalar(overlapScalars[overlap])]);
            integratedErrors[sample].push_back(results[set].estimates[sample].error[results[set].layout.scalar(overlapScalars[overlap])]);
          }
        }
      }
      if (!integratedLabels.empty()) {
        std::vector<DrawEntry> integratedEntries;
        for (int sample = 0; sample < kNSamples; ++sample) {
          if (sampleFiles[sample] == nullptr)
            continue;
          const Part part = tableHistogram(ringOverlap.parts, "hOverlapIntegrated_" + kSampleKeys[sample], "Integrated ring overlap (d^{3}N-weighted), " + sampleLegends[sample] + ", jackknife errors", integratedLabels, integratedValues[sample], integratedErrors[sample], blockIndices(0, integratedLabels.size()), 1., warnings);
          integratedEntries.push_back({part, sampleLegends[sample], kSampleLooks[sample], kStyleMarker});
        }
        drawOverlay(ringOverlap.canvases, "cOverlapIntegrated", "Integrated ring overlap with the true trigger, jackknife errors", integratedEntries, 0.);
      }

      // Leakage ansatz. Scan: R'(delta) against the prediction, and the residual. Mixing and smudging: integrated only (with constant A and B
      // the ansatz cannot follow R(Delta phi') bin by bin)
      const std::string predictionLabel = "#it{R}_{true}<cos#zeta> + <#it{P}_{e}><sin#zeta>";
      if (results[kSetScan].present) {
        const LeakageLayout leakageLayout{results[kSetScan].layout.nBins};
        std::vector<DrawEntry> residualEntries;
        for (const int sample : presentSamples) {
          const EstimateSet& sampleLeakage = leakage[kSetScan][sample];
          if (sampleLeakage.value.empty())
            continue;
          const std::string key = kSampleKeys[sample];
          const Part prediction = curveHistogram(ringOverlap.parts, "hLeakagePredictionVsDelta_" + key, "Leakage ansatz " + predictionLabel + " (delta-scan, " + sampleLegends[sample] + ", jackknife errors);" + deltaRange.title + ";prediction", deltaRange, sampleLeakage.value, sampleLeakage.error, leakageLayout.prediction(0), 1., warnings);
          const Part residual = curveHistogram(ringOverlap.parts, "hLeakageResidualVsDelta_" + key, "#it{R}' - prediction (delta-scan, " + sampleLegends[sample] + ", joint jackknife errors);" + deltaRange.title + ";#it{R}' - prediction", deltaRange, sampleLeakage.value, sampleLeakage.error, leakageLayout.residual(0), 1., warnings);
          drawOverlay(ringOverlap.canvases, "cLeakageVsDelta_" + key, "Leakage ansatz against the delta-scan, " + sampleLegends[sample] + ", jackknife bands", {{curveParts[kSetScan][kCurveRing][sample], "#it{R}'(#delta), delta-scan", kPaletteLooks[0], kStyleBand}, {prediction, predictionLabel, kPaletteLooks[1], kStyleBand}}, 0.);
          residualEntries.push_back({residual, sampleLegends[sample], kSampleLooks[sample], kStyleBand});
        }
        drawOverlay(ringOverlap.canvases, "cLeakageResidualVsDelta", "#it{R}' - (" + predictionLabel + ") vs #delta, joint jackknife bands", residualEntries, 0.);
      }
      for (const int set : {kSetMixing, kSetSmudge}) {
        if (!results[set].present)
          continue;
        const LeakageLayout leakageLayout{results[set].layout.nBins};
        const int integral = results[set].layout.nBins;
        const std::vector<std::string> leakageLabels = {"#it{R}'_{d^{3}N weight}", "prediction", "#it{R}' - prediction"};
        std::vector<DrawEntry> leakageEntries;
        for (const int sample : presentSamples) {
          const EstimateSet& sampleLeakage = leakage[set][sample];
          if (sampleLeakage.value.empty())
            continue;
          const EstimateSet& variant = results[set].estimates[sample];
          const std::vector<double> values = {variant.value[results[set].layout.scalar(kScalarRingD3N)], sampleLeakage.value[leakageLayout.prediction(integral)], sampleLeakage.value[leakageLayout.residual(integral)]};
          const std::vector<double> errors = {variant.error[results[set].layout.scalar(kScalarRingD3N)], sampleLeakage.error[leakageLayout.prediction(integral)], sampleLeakage.error[leakageLayout.residual(integral)]};
          const Part part = tableHistogram(ringOverlap.parts, "hLeakageIntegrated_" + kSetKeys[set] + "_" + kSampleKeys[sample], "Leakage ansatz (" + predictionLabel + "), " + kSetLabels[set] + ", " + sampleLegends[sample] + ", jackknife errors", leakageLabels, values, errors, blockIndices(0, 3), 1., warnings);
          leakageEntries.push_back({part, sampleLegends[sample], kSampleLooks[sample], kStyleMarker});
          const std::string line = kSetLabels[set] + ", " + sampleLabels[sample == kSampleDifference ? kSampleBullet : sample] + (sample == kSampleDifference ? " - no_bullet" : "") + ": R' = " + std::to_string(values[0]) + " +- " + std::to_string(errors[0]) + ", prediction = " + std::to_string(values[1]) + " +- " + std::to_string(errors[1]) + ", residual = " + std::to_string(values[2]) + " +- " + std::to_string(errors[2]) + "\n";
          leakageSummary += line;
        }
        drawOverlay(ringOverlap.canvases, "cLeakageIntegrated_" + kSetKeys[set], "Leakage ansatz, " + kSetLabels[set] + ": #it{R}' against " + predictionLabel + ", jackknife errors", leakageEntries, 0.);
      }
      std::cout << "Leakage ansatz, integrated (residual = Cov(P_R, cos zeta) + Cov(P_e, sin zeta)):\n" << leakageSummary;
    }

    // ==================================================================================================================================
    // 16. PolarizationMaps/: the producer's maps of the true trigger, in percent, as vector fields (README, "Polarization maps")
    // ==================================================================================================================================
    {
      const std::array<std::string, kNCurveSamples> mapLegends = {sampleLegends[kSampleBullet], sampleLegends[kSampleNoBullet], "bullet - no_bullet (pooled maps, errors in quadrature)"};
      const std::string inPlaneLabel = "|<P*_{in-plane}>|";
      // Lab and PrimeJet: three planes each. Producer names p2d<component>Star<tag>_vs<plane><tag> and h2dYieldVs<plane><tag>
      const std::array<std::string, 2> systemKeys = {"Lab", "PrimeJet"};
      const std::array<std::array<std::string, 3>, 2> planeBanners = {{{"X-Y plane", "Z-X plane", "Y-Z plane"}, {"X'-Y' plane (transverse to the jet)", "Z'-X' plane (the beam-jet plane)", "Y'-Z' plane"}}};
      const std::array<std::array<std::string, 3>, 2> magnitudeLabels = {{{inPlaneLabel, inPlaneLabel, inPlaneLabel}, {"|<P*_{#perp Jet}>|", inPlaneLabel, inPlaneLabel}}};
      const std::array<std::string, 3> planeKeys = {"PxPy", "PzPx", "PyPz"};
      const std::array<std::string, 3> componentKeys = {"Px", "Py", "Pz"};
      // Arrow components of each plane (along its horizontal and vertical axes), and the out-of-plane component as colour
      const std::array<std::array<int, 3>, 3> planeComponents = {{{0, 1, 2}, {2, 0, 1}, {1, 2, 0}}};
      for (int system = 0; system < 2; ++system) {
        const OutputFolder folder = makeFolder(polarizationMapsDir, systemKeys[system]);
        const std::string tag = system == 0 ? "" : "PrimeJet";
        const std::string producerFolder = "PolarizationMaps/" + systemKeys[system] + "/";
        // [plane][component][sample], [plane][sample] for the yields, [plane][curve][sample] for the ring scalars (PrimeJet only)
        std::array<std::array<std::array<TH2D*, kNCurveSamples>, 3>, 3> maps{};
        std::array<std::array<TH2D*, kNSamples>, 3> yields{};
        std::array<std::array<std::array<TH2D*, kNCurveSamples>, kNFeatureCurves>, 3> ringMaps{};
        for (int plane = 0; plane < 3; ++plane) {
          for (int component = 0; component < 3; ++component) {
            const std::string producerName = "p2d" + componentKeys[component] + "Star" + tag + "_vs" + planeKeys[plane] + tag;
            maps[plane][component] = sampleMaps(folder.parts, sampleFiles, producerFolder + producerName, "h2d" + producerName.substr(3), mapLegends);
          }
          if (system == 1)
            for (int curve = 0; curve < kNFeatureCurves; ++curve)
              ringMaps[plane][curve] = sampleMaps(folder.parts, sampleFiles, producerFolder + "p2d" + kFeatureCurveKeys[curve] + "Vs" + planeKeys[plane] + tag, "h2d" + kFeatureCurveKeys[curve] + "Vs" + planeKeys[plane] + tag, mapLegends);
          for (int sample = 0; sample < kNSamples; ++sample)
            if (sampleFiles[sample] != nullptr)
              yields[plane][sample] = copyMap(folder.parts, sampleFiles[sample], producerFolder + "h2dYieldVs" + planeKeys[plane] + tag, "h2dYieldVs" + planeKeys[plane] + tag + "_" + kSampleKeys[sample], ", " + sampleLegends[sample]);
        }
        for (const int sample : presentSamples) {
          // The difference maps use the bullet yield for the populated bins and the tile weights
          const int yieldSample = sample == kSampleDifference ? kSampleBullet : sample;
          const std::string key = kSampleKeys[sample];
          std::vector<MapPanel> vectorPanels, yieldPanels;
          std::array<std::vector<MapPanel>, kNFeatureCurves> ringPanels;
          for (int plane = 0; plane < 3; ++plane) {
            const std::array<int, 3>& components = planeComponents[plane];
            vectorPanels.push_back({kPanelVectorField, maps[plane][components[2]][sample], maps[plane][components[0]][sample], maps[plane][components[1]][sample], yields[plane][yieldSample], planeBanners[system][plane], magnitudeLabels[system][plane]});
            if (sample != kSampleDifference)
              yieldPanels.push_back({kPanelYield, yields[plane][sample], nullptr, nullptr, nullptr, planeBanners[system][plane], ""});
            if (system == 1)
              for (int curve = 0; curve < kNFeatureCurves; ++curve)
                ringPanels[curve].push_back({kPanelSymmetric, ringMaps[plane][curve][sample], nullptr, nullptr, yields[plane][yieldSample], planeBanners[system][plane], ""});
          }
          drawPanelCanvas(folder.canvases, "cVectorField" + systemKeys[system] + "_" + key, "<P*> in " + systemKeys[system] + " coordinates (%), " + mapLegends[sample], vectorPanels);
          if (sample != kSampleDifference)
            drawPanelCanvas(folder.canvases, "cYield" + systemKeys[system] + "_" + key, "Yields in " + systemKeys[system] + " coordinates, " + mapLegends[sample], yieldPanels);
          if (system == 1)
            for (int curve = 0; curve < kNFeatureCurves; ++curve)
              drawPanelCanvas(folder.canvases, "c" + kFeatureCurveKeys[curve] + "PrimeJet_" + key, "<" + kFeatureCurveSymbols[curve] + "> in PrimeJet coordinates (%), " + mapLegends[sample], ringPanels[curve]);
        }
      }
      // PrimeV0: the production plane (p_z, p_T), whose own in-plane axes are z and x'V0. The colour is P*_y'V0, the ring observable with the
      // beam as the trigger
      const OutputFolder folder = makeFolder(polarizationMapsDir, "PrimeV0");
      const std::array<TH2D*, kNCurveSamples> horizontalMaps = sampleMaps(folder.parts, sampleFiles, "PolarizationMaps/PrimeV0/p2dPzStar_vsPzPt", "h2dPzStar_vsPzPt", mapLegends);
      const std::array<TH2D*, kNCurveSamples> verticalMaps = sampleMaps(folder.parts, sampleFiles, "PolarizationMaps/PrimeV0/p2dPxStarPrimeV0_vsPzPt", "h2dPxStarPrimeV0_vsPzPt", mapLegends);
      const std::array<TH2D*, kNCurveSamples> backgroundMaps = sampleMaps(folder.parts, sampleFiles, "PolarizationMaps/PrimeV0/p2dPyStarPrimeV0_vsPzPt", "h2dPyStarPrimeV0_vsPzPt", mapLegends);
      std::array<TH2D*, kNSamples> primeV0Yields{};
      for (int sample = 0; sample < kNSamples; ++sample)
        if (sampleFiles[sample] != nullptr)
          primeV0Yields[sample] = copyMap(folder.parts, sampleFiles[sample], "PolarizationMaps/PrimeV0/h2dYieldVsPzPt", "h2dYieldVsPzPt_" + kSampleKeys[sample], ", " + sampleLegends[sample]);
      for (const int sample : presentSamples) {
        const int yieldSample = sample == kSampleDifference ? kSampleBullet : sample;
        std::vector<MapPanel> panels = {{kPanelVectorField, backgroundMaps[sample], horizontalMaps[sample], verticalMaps[sample], primeV0Yields[yieldSample], "#Lambda production plane", inPlaneLabel}};
        if (sample != kSampleDifference)
          panels.push_back({kPanelYield, primeV0Yields[sample], nullptr, nullptr, nullptr, "Yield", ""});
        drawPanelCanvas(folder.canvases, "cVectorFieldPrimeV0_" + kSampleKeys[sample], "<P*> in PrimeV0 coordinates (%), " + mapLegends[sample], panels);
      }
    }

    // ==================================================================================================================================
    // 17. Precision/: expected sigma per Delta phi bin for an ideal detector (README, "Expected precision"), from the bullet sample
    // ==================================================================================================================================
    {
      const EstimateSet& bulletTrue = results[kSetTrue].estimates[kSampleBullet];
      const TProfile* nzSqSinSqThetaProfile = getObject<TProfile>(bulletFile.get(), "Ring/pNzSqSinSqThetaVsDeltaPhi");
      const std::array<std::string, 3> precisionKeys = {"Ring", "RingZ", "RingZT"};
      const std::array<std::string, 3> precisionSymbols = {"#it{R}", "#it{R}_{z}", "#it{R}_{z}^{T}"};
      // sigma_b sqrt(N) = sqrt(3 <|e|^2>_b / f_b) / |alpha|, with |e|^2 = 1, n_z^2 and n_z^2 sin^2 theta, and f_b the yield fraction of the bin
      std::array<std::vector<double>, 3> sigmaSqrtN;
      for (std::vector<double>& values : sigmaSqrtN)
        values.assign(nDeltaPhiBins, kNoValue);
      for (int bin = 0; bin < nDeltaPhiBins; ++bin) {
        const double yieldFraction = bulletTrue.value[trueLayout.curve(kCurveYieldFraction, bin)];
        const double nzSq = bulletTrue.value[trueLayout.curve(kCurveNzSq, bin)];
        const std::array<double, 3> eSquared = {1., nzSq, nzSqSinSqThetaProfile->GetBinContent(bin + 1)};
        if (!(yieldFraction > 0.) || nzSq == kNoValue)
          continue;
        for (int observable = 0; observable < 3; ++observable)
          sigmaSqrtN[observable][bin] = std::sqrt(3. * eSquared[observable] / yieldFraction) / std::abs(config.alphaLambda);
      }
      const std::array<double, 3> eSquaredIntegrated = {1., bulletTrue.value[trueLayout.scalar(kScalarNzSqD3N)], bulletTrue.value[trueLayout.scalar(kScalarNzSqSinSqThetaD3N)]};
      std::vector<double> sigmaSqrtNIntegrated(3, kNoValue);
      for (int observable = 0; observable < 3; ++observable)
        if (eSquaredIntegrated[observable] != kNoValue)
          sigmaSqrtNIntegrated[observable] = std::sqrt(3. * eSquaredIntegrated[observable]) / std::abs(config.alphaLambda);
      const std::vector<double> noErrors(nDeltaPhiBins, 0.);
      const std::vector<std::string> precisionLabels(precisionSymbols.begin(), precisionSymbols.end());
      std::vector<DrawEntry> sigmaEntries;
      for (int observable = 0; observable < 3; ++observable) {
        const Part part = curveHistogram(precision.parts, "hSigmaSqrtN_" + precisionKeys[observable], "#sigma #sqrt{N} of " + precisionSymbols[observable] + " per bin (ideal detector, bullet yields);" + deltaPhiRange.title + ";#sigma #sqrt{N}", deltaPhiRange, sigmaSqrtN[observable], noErrors, 0, 1., warnings);
        sigmaEntries.push_back({part, precisionSymbols[observable], kPaletteLooks[observable], kStyleLine});
      }
      drawOverlay(precision.canvases, "cSigmaSqrtN", "Expected #sigma #sqrt{N} per #Delta#varphi bin (ideal detector)", sigmaEntries, kNoValue);
      const Part sigmaSqrtNTable = tableHistogram(precision.parts, "hSigmaSqrtNIntegrated", "Expected #sigma #sqrt{N} of the integrals (ideal detector);;#sigma #sqrt{N}", precisionLabels, sigmaSqrtNIntegrated, std::vector<double>(3, 0.), blockIndices(0, 3), 1., warnings);
      drawOverlay(precision.canvases, "cSigmaSqrtNIntegrated", "Expected #sigma #sqrt{N} of the integrals (ideal detector)", {{sigmaSqrtNTable, "#sigma #sqrt{N}", kPaletteLooks[0], kStyleMarker}}, kNoValue);
      if (config.nLambdaPrecision > 0.) {
        const double sqrtN = std::sqrt(config.nLambdaPrecision);
        std::vector<DrawEntry> significanceEntries;
        for (int observable = 0; observable < 3; ++observable) {
          std::vector<double> sigma(nDeltaPhiBins, kNoValue);
          for (int bin = 0; bin < nDeltaPhiBins; ++bin)
            if (sigmaSqrtN[observable][bin] != kNoValue)
              sigma[bin] = sigmaSqrtN[observable][bin] / sqrtN;
          curveHistogram(precision.parts, "hSigma_" + precisionKeys[observable], "Expected #sigma of " + precisionSymbols[observable] + " per bin, N = " + std::to_string(config.nLambdaPrecision) + ";" + deltaPhiRange.title + ";#sigma", deltaPhiRange, sigma, noErrors, 0, 1., warnings);
          // Expected significance of the jet-induced signal: the difference over the error of one measurement of R (or R_z) with N Lambdas.
          // R_z^T has no per-bin difference (no per-bin sums of it in the tree)
          if (!hasNoBullet || observable == 2)
            continue;
          const EstimateSet& differenceEstimate = results[kSetTrue].estimates[kSampleDifference];
          std::vector<double> significance(nDeltaPhiBins, kNoValue);
          for (int bin = 0; bin < nDeltaPhiBins; ++bin) {
            const double signal = differenceEstimate.value[trueLayout.curve(observable == 0 ? kCurveRing : kCurveRingZ, bin)];
            if (sigma[bin] != kNoValue && signal != kNoValue && sigma[bin] > 0.)
              significance[bin] = signal / sigma[bin];
          }
          const Part part = curveHistogram(precision.parts, "hSignificance_" + precisionKeys[observable], "Expected significance of the bullet - no_bullet " + precisionSymbols[observable] + " per bin, N = " + std::to_string(config.nLambdaPrecision) + ";" + deltaPhiRange.title + ";signal/#sigma", deltaPhiRange, significance, noErrors, 0, 1., warnings);
          significanceEntries.push_back({part, precisionSymbols[observable], kPaletteLooks[observable], kStyleLine});
        }
        drawOverlay(precision.canvases, "cSignificance", "Expected significance of the jet-induced signal per #Delta#varphi bin", significanceEntries, 0.);
        std::vector<double> sigmaIntegrated(3, kNoValue), significanceIntegrated(3, kNoValue);
        const std::array<int, 3> differenceScalars = {kScalarRingD3N, kScalarRingZD3N, kScalarRingZTD3N};
        for (int observable = 0; observable < 3; ++observable) {
          if (sigmaSqrtNIntegrated[observable] == kNoValue)
            continue;
          sigmaIntegrated[observable] = sigmaSqrtNIntegrated[observable] / sqrtN;
          if (hasNoBullet) {
            const double signal = results[kSetTrue].estimates[kSampleDifference].value[trueLayout.scalar(differenceScalars[observable])];
            if (signal != kNoValue)
              significanceIntegrated[observable] = signal / sigmaIntegrated[observable];
          }
        }
        tableHistogram(precision.parts, "hSigmaIntegrated", "Expected #sigma of the integrals, N = " + std::to_string(config.nLambdaPrecision) + ";;#sigma", precisionLabels, sigmaIntegrated, std::vector<double>(3, 0.), blockIndices(0, 3), 1., warnings);
        if (hasNoBullet) {
          const Part significanceTable = tableHistogram(precision.parts, "hSignificanceIntegrated", "Expected significance of the bullet - no_bullet integrals, N = " + std::to_string(config.nLambdaPrecision) + ";;signal/#sigma", precisionLabels, significanceIntegrated, std::vector<double>(3, 0.), blockIndices(0, 3), 1., warnings);
          drawOverlay(precision.canvases, "cSignificanceIntegrated", "Expected significance of the bullet - no_bullet integrals, N = " + std::to_string(config.nLambdaPrecision), {{significanceTable, "signal/#sigma", kPaletteLooks[0], kStyleMarker}}, 0.);
        }
      }
    }

    // ==================================================================================================================================
    // 18. Jackknife/: validation of the error treatment
    // ==================================================================================================================================
    // Pooled-master error of the producer profiles over the jackknife error, per bin: close to 1 expected (both estimate the error of a ratio
    // of sums, the profile through the effective number of events)
    for (int set = 0; set < kNTriggerSets; ++set) {
      if (!results[set].present)
        continue;
      const AxisRange& range = set == kSetScan ? deltaRange : deltaPhiRange;
      std::vector<DrawEntry> ratioEntries;
      for (int sample = 0; sample < kNSamples; ++sample) {
        if (sampleFiles[sample] == nullptr)
          continue;
        const EstimateSet& estimates = results[set].estimates[sample];
        for (int curve = 0; curve < kNFeatureCurves; ++curve) {
          const std::string profilePath = kSetProducerFolder[set] + "/" + (curve == kCurveRing ? kSetRingProfile[set] : kSetRingZProfile[set]);
          const TProfile* profile = getObject<TProfile>(sampleFiles[sample], profilePath);
          if (profile->GetNbinsX() != range.nBins)
            throw std::runtime_error(profilePath + " does not have the binning of the event tree");
          std::vector<double> errorRatio(range.nBins, kNoValue);
          for (int bin = 0; bin < range.nBins; ++bin) {
            const double jackknifeError = estimates.error[results[set].layout.curve(curve, bin)];
            if (jackknifeError != kNoValue && jackknifeError > 0.)
              errorRatio[bin] = profile->GetBinError(bin + 1) / jackknifeError;
          }
          const Part part = curveHistogram(jackknifeFolder.parts, "hErrorRatio" + kFeatureCurveKeys[curve] + "_" + kSetKeys[set] + "_" + kSampleKeys[sample], "#sigma_{pooled master}/#sigma_{jackknife} of " + kFeatureCurveSymbols[curve] + " (" + kSetLabels[set] + ", " + sampleLabels[sample] + ");" + range.title + ";#sigma_{pooled master}/#sigma_{jackknife}", range, errorRatio, std::vector<double>(range.nBins, 0.), 0, 1., warnings);
          ratioEntries.push_back({part, kFeatureCurveSymbols[curve] + ", " + sampleLabels[sample], kPaletteLooks[ratioEntries.size()], kStyleLine});
        }
      }
      drawOverlay(jackknifeFolder.canvases, "cErrorRatio_" + kSetKeys[set], "Pooled-master over jackknife error (" + kSetLabels[set] + ")", ratioEntries, 1.);
    }
    // The same for the single-bin profiles. The _EvtMean ones are plain means, whose jackknife is the SEM with the N - 1 spread, while TProfile
    // uses the population spread: their ratio is sqrt((N - 1)/N) exactly, a check of the jackknife itself
    {
      const std::array<std::pair<std::string, int>, 6> integratedProfiles = {{{"pRingD3N", kScalarRingD3N}, {"pRingZD3N", kScalarRingZD3N}, {"pRingD3N_EvtMean", kScalarRingD3NEvtMean}, {"pRingFlat_EvtMean", kScalarRingFlatEvtMean}, {"pRingZD3N_EvtMean", kScalarRingZD3NEvtMean}, {"pRingZFlat_EvtMean", kScalarRingZFlatEvtMean}}};
      std::vector<std::string> ratioLabels;
      std::vector<double> ratioValues;
      for (const int set : {kSetTrue, kSetMixing, kSetSmudge}) {
        if (!results[set].present)
          continue;
        for (int sample = 0; sample < kNSamples; ++sample) {
          if (sampleFiles[sample] == nullptr)
            continue;
          const EstimateSet& estimates = results[set].estimates[sample];
          for (const auto& [profileName, scalar] : integratedProfiles) {
            const double jackknifeError = estimates.error[results[set].layout.scalar(scalar)];
            const double profileError = getObject<TProfile>(sampleFiles[sample], kSetProducerFolder[set] + "/" + profileName)->GetBinError(1);
            ratioLabels.push_back(profileName + ", " + kSetKeys[set] + ", " + sampleLabels[sample]);
            ratioValues.push_back(jackknifeError != kNoValue && jackknifeError > 0. ? profileError / jackknifeError : kNoValue);
          }
        }
      }
      const Part ratioTable = tableHistogram(jackknifeFolder.parts, "hErrorRatioIntegrated", "Single-bin profiles: #sigma_{profile}/#sigma_{jackknife} (_EvtMean: #sqrt{(N-1)/N} expected);;#sigma_{profile}/#sigma_{jackknife}", ratioLabels, ratioValues, std::vector<double>(ratioValues.size(), 0.), blockIndices(0, ratioLabels.size()), 1., warnings);
      drawOverlay(jackknifeFolder.canvases, "cErrorRatioIntegrated", "Single-bin profiles: profile over jackknife error", {{ratioTable, "#sigma_{profile}/#sigma_{jackknife}", kPaletteLooks[0], kStyleMarker}}, 1.);
    }
    // Jackknife bias over error for the true trigger: well below 1 means the full-sample values need no bias correction
    for (int curve = 0; curve < kNFeatureCurves; ++curve) {
      std::vector<DrawEntry> biasEntries;
      for (const int sample : presentSamples) {
        const EstimateSet& estimates = results[kSetTrue].estimates[sample];
        std::vector<double> biasOverError(estimates.value.size(), kNoValue);
        for (std::size_t index = 0; index < biasOverError.size(); ++index)
          if (estimates.bias[index] != kNoValue && estimates.error[index] != kNoValue && estimates.error[index] > 0.)
            biasOverError[index] = estimates.bias[index] / estimates.error[index];
        const Part part = curveHistogram(jackknifeFolder.parts, "hBiasOverError" + kFeatureCurveKeys[curve] + "_" + kSampleKeys[sample], "Jackknife bias/error of " + kFeatureCurveSymbols[curve] + " (true trigger, " + sampleLegends[sample] + ");" + deltaPhiRange.title + ";bias/#sigma", deltaPhiRange, biasOverError, std::vector<double>(biasOverError.size(), 0.), trueLayout.curve(curve, 0), 1., warnings);
        biasEntries.push_back({part, sampleLegends[sample], kSampleLooks[sample], kStyleLine});
      }
      drawOverlay(jackknifeFolder.canvases, "cBiasOverError" + kFeatureCurveKeys[curve], "Jackknife bias over error of " + kFeatureCurveSymbols[curve] + " (true trigger)", biasEntries, 0.);
    }

    // ==================================================================================================================================
    // 19. QA/: folded curves and asymmetries R(dphi) - R(-dphi). Symmetric on average only (README, "Errors"), so never folded elsewhere
    // ==================================================================================================================================
    {
      const AxisRange foldRange{trueLayout.nFoldBins(), -0.5 * deltaPhiWidth, (trueLayout.nFoldBins() - 0.5) * deltaPhiWidth, "|#Delta#varphi| (rad)"};
      const AxisRange asymmetryRange{trueLayout.nAsymmetryBins(), 0.5 * deltaPhiWidth, (trueLayout.nAsymmetryBins() + 0.5) * deltaPhiWidth, "|#Delta#varphi| (rad)"};
      for (int curve = 0; curve < kNFeatureCurves; ++curve) {
        const std::string symbol = kFeatureCurveSymbols[curve];
        std::vector<DrawEntry> foldEntries, asymmetryEntries;
        for (const int sample : presentSamples) {
          const EstimateSet& estimates = results[kSetTrue].estimates[sample];
          const Part fold = curveHistogram(qa.parts, "h" + kFeatureCurveKeys[curve] + "Folded_" + kSampleKeys[sample], "[" + symbol + "(#Delta#varphi) + " + symbol + "(-#Delta#varphi)]/2 (" + sampleLegends[sample] + ", jackknife errors);" + foldRange.title + ";" + symbol + ", folded", foldRange, estimates.value, estimates.error, trueLayout.fold(curve, 0), 1., warnings);
          const Part asymmetry = curveHistogram(qa.parts, "h" + kFeatureCurveKeys[curve] + "Asymmetry_" + kSampleKeys[sample], symbol + "(#Delta#varphi) - " + symbol + "(-#Delta#varphi) (" + sampleLegends[sample] + ", jackknife errors);" + asymmetryRange.title + ";" + symbol + "(#Delta#varphi) - " + symbol + "(-#Delta#varphi)", asymmetryRange, estimates.value, estimates.error, trueLayout.asymmetry(curve, 0), 1., warnings);
          foldEntries.push_back({fold, sampleLegends[sample], kSampleLooks[sample], kStyleBand});
          asymmetryEntries.push_back({asymmetry, sampleLegends[sample], kSampleLooks[sample], kStyleBand});
        }
        drawOverlay(qa.canvases, "c" + kFeatureCurveKeys[curve] + "Folded", "Folded " + symbol + " (QA only: the curve is symmetric on average, not event by event)", foldEntries, 0.);
        drawOverlay(qa.canvases, "c" + kFeatureCurveKeys[curve] + "Asymmetry", "Asymmetry " + symbol + "(#Delta#varphi) - " + symbol + "(-#Delta#varphi), consistent with 0 expected", asymmetryEntries, 0.);
      }
    }

    // ==================================================================================================================================
    // 20. Config/, warnings, write
    // ==================================================================================================================================
    std::string warningSummary;
    for (const std::string& warningLine : warnings) {
      std::cout << "[WARNING] " << warningLine << std::endl;
      warningSummary += warningLine + "\n";
    }
    configDir->cd();
    for (const auto& [name, value] : configEntries)
      TNamed(name.c_str(), value.c_str()).Write();
    TNamed("trigger", triggerMode.c_str()).Write();
    TNamed("nEventsBullet", std::to_string(nBulletEvents).c_str()).Write();
    TNamed("nEventsNoBullet", std::to_string(nNoBulletEvents).c_str()).Write();
    TNamed("nPairs", std::to_string(nPairs).c_str()).Write();
    TNamed("noBulletClosure", noBulletClosure.c_str()).Write();
    TNamed("deltaScanClosures", scanClosureSummary.c_str()).Write();
    TNamed("mixingPulls", pullSummary.c_str()).Write();
    TNamed("leakageTest", leakageSummary.c_str()).Write();
    TNamed("warnings", warningSummary.empty() ? "none" : warningSummary.c_str()).Write();

    outputFile->Write();
    outputFile->Close();
    std::cout << "Output written to " << config.outputFile << std::endl;
  } catch (const std::exception& error) {
    std::cerr << "[ERROR] " << error.what() << std::endl;
    return 1;
  }
  return 0;
}