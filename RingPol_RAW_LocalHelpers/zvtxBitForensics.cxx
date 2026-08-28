// ============================================================================
// zvtxBitForensics.cxx
//
// Bit-level forensics and integrity QA on the derived AO2D tables that feed
// the event-mixing machinery of lambdaJetPolarizationIonsDerived.cxx.
//
// The physics motivation, the interpretation of every histogram and the
// birthday-paradox model behind the "expected" counters are documented in
// RingPol_RAW_LocalHelpers/README.md. This file documents only the logic it
// needs in order to run.
//
// Two modes:
//
//   Worker mode (one process per batch of AO2D files):
//     zvtxBitForensics.exe <manifest.txt> <output.root>
//   where manifest.txt holds one absolute AO2D path per line.
//
//   Finalize mode (once, after hadd has merged every batch output):
//     zvtxBitForensics.exe --finalize <merged.root>
//   which divides the additive numerator/denominator counters written by the
//   workers and stores the resulting ratio histograms back into the same file.
//
// A ratio is not additive under hadd, so every ratio in this tool is stored as
// a (numerator, denominator) pair by the workers and only divided in finalize
// mode. Do not add a histogram here that is itself a ratio.
// ============================================================================

#include <TAxis.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TKey.h>
#include <TList.h>
#include <TProfile.h>
#include <TString.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

// ============================================================================
// SECTION 1 -- IEEE-754 helpers
// ============================================================================

/// \brief Reinterprets the four raw bytes of a float as an unsigned 32-bit word.
/// \note  memcpy rather than a pointer cast, so no strict-aliasing assumption is made.
static uint32_t floatBits(float value)
{
  uint32_t bits = 0;
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

/// \brief Counts how many low-order bits of the 23-bit mantissa field are zero.
/// \return -1 when the mantissa is entirely zero (exact zero or exact power of
///         two), which carries no truncation information and must be excluded
///         from any minimum taken over the sample.
static int countTrailingZeroMantissaBits(float value)
{
  const uint32_t mantissa = floatBits(value) & 0x007FFFFFu;
  if (mantissa == 0u) {
    return -1;
  }
  int trailingZeros = 0;
  uint32_t work = mantissa;
  while ((work & 1u) == 0u) {
    work >>= 1;
    ++trailingZeros;
  }
  return trailingZeros;
}

/// \brief Spacing between adjacent representable values at a given magnitude,
///        for a float carrying keptBits informative mantissa bits.
/// \note  For 2^e <= |x| < 2^(e+1) the float32 spacing is 2^(e-23); clearing the
///        low (23 - keptBits) bits scales that up to 2^(e - keptBits).
static double localGridStep(double value, int keptBits)
{
  if (value == 0.0 || !std::isfinite(value)) {
    return 0.0;
  }
  return std::ldexp(1.0, std::ilogb(std::fabs(value)) - keptBits);
}

// ============================================================================
// SECTION 2 -- Configuration
// ============================================================================

/// \brief One float column to be analysed, and where to find it.
/// \note  axisLow/axisHigh fix the profile binning at book time. Auto-ranging
///        would give each worker a different axis and make hadd unmergeable.
struct ColumnSpec {
  const char* label;      //! short name used in histogram names
  const char* treeName;   //! tree inside each DF_* directory
  const char* branchName; //! branch inside that tree
  double axisLow;         //! lower edge of the value axis
  double axisHigh;        //! upper edge of the value axis
};

// Columns analysed for truncation and duplicate structure. The collision
// columns come first because the fingerprint test consumes them in this order.
static const std::vector<ColumnSpec> kColumns = {
  {"Zvtx", "O2ringcollision", "fZvtx", -15.0, 15.0},
  {"CentFT0M", "O2ringcollision", "fCentFT0M", -5.0, 105.0},
  {"CentFT0C", "O2ringcollision", "fCentFT0C", -5.0, 105.0},
  {"CentFV0A", "O2ringcollision", "fCentFV0A", -5.0, 105.0},
  {"IntRate", "O2ringcollision", "fInteractionRate", 0.0, 1.0e6},
  {"JetPt", "O2ringjet", "fJetPt", -20.0, 200.0},
  {"LeadPPt", "O2ringleadp", "fLeadParticlePt", 0.0, 200.0}};

// The collision columns participate in the cumulative fingerprint test, in the
// order above. Five columns means five cumulative levels.
static const int kNFingerprintColumns = 5;

// Index column shared by every non-collision derived table.
static const char* kIndexBranch = "fIndexRingCollisions";

// Tables carrying that index, checked by the integrity module.
static const std::vector<const char*> kIndexedTrees = {"O2ringjet", "O2ringleadp", "O2ringlav0"};

// Mixing-pool axes. These MUST stay in sync with axisConfigurations in
// lambdaJetPolarizationIonsDerived.cxx, otherwise the occupancy prediction is
// not comparable to the consumer's own mixing QA.
// TODO: read these from the dpl-config JSON instead of hardcoding them, and
// extend this module into a full reproduction of the SameKindPair windowing so
// the predicted candidate distribution can be compared bin by bin against
// EventMixingQA/CollLoopOutcome/hMixedEventLeadPCandidates.
static const double kAxisPVzLow = -10.0;
static const double kAxisPVzHigh = 10.0;
static const int kAxisPVzBins = 60;
static const std::vector<double> kAxisCentralityEdges = {0.0, 20.0, 50.0, 100.0};
static const std::vector<double> kAxisJetPtEdges = {0, 2, 4, 6, 8, 10, 15, 20, 30, 40, 60, 80, 120, 160, 200};

/// \brief Bin index for a value on a variable-width axis, or -1 when outside.
static int findVariableBin(const std::vector<double>& edges, double value)
{
  if (value < edges.front() || value >= edges.back()) {
    return -1;
  }
  const auto upper = std::upper_bound(edges.begin(), edges.end(), value);
  return static_cast<int>(upper - edges.begin()) - 1;
}

/// \brief Bin index for a value on a uniform axis, or -1 when outside.
static int findUniformBin(double low, double high, int nBins, double value)
{
  if (value < low || value >= high) {
    return -1;
  }
  return static_cast<int>((value - low) / (high - low) * nBins);
}

// ============================================================================
// SECTION 3 -- Histogram registry
// ============================================================================

/// \brief Every histogram written by a worker. Booked once, filled per dataframe.
/// \note  Counter histograms come in numerator/denominator pairs because hadd
///        can only sum; finalize mode turns each pair into a ratio.
struct ForensicsHistograms {
  // --- Truncation ---
  std::vector<TH1D*> hTrailingZeros;         // one per column
  std::vector<TH2D*> hTrailingZerosVsBinade; // one per column

  // --- Duplicates ---
  std::vector<TH1D*> hMultiplicity;
  std::vector<TH1D*> hLogSpacing;
  std::vector<TH1D*> hSpacingOverGridStep;
  std::vector<TH1D*> hDuplicateRowGap;
  std::vector<TProfile*> pMultiplicityVsValue;

  // --- Additive counters, per column (bin i = column i) ---
  TH1D* hEntriesPerColumn = nullptr;
  TH1D* hDistinctPerColumn = nullptr;
  TH1D* hObservedDuplicatePairs = nullptr;
  TH1D* hExpectedDuplicatePairs = nullptr;
  TH1D* hNonFinitePerColumn = nullptr;
  TH1D* hConstantColumnDataframes = nullptr;

  // --- Birthday-model scaling test: counters binned in dataframe size ---
  TH1D* hObservedPairsVsDfSize = nullptr;
  TH1D* hExpectedPairsVsDfSize = nullptr;

  // --- Fingerprints (bin k = match on the first k collision columns) ---
  TH1D* hObservedFingerprintPairs = nullptr;
  TH1D* hExpectedFingerprintPairs = nullptr;
  TH1D* hFingerprintGroupSize = nullptr;
  TH1D* hFingerprintRowGap = nullptr;

  // --- Integrity ---
  TH1D* hIntegrityViolations = nullptr;
  TH1D* hRowsPerCollisionLeadP = nullptr;
  TH1D* hRowsPerCollisionJet = nullptr;
  TH1D* hRowsPerCollisionV0 = nullptr;
  TH1D* hCollisionContentPattern = nullptr;
  TH2D* hIndexRangeExcess = nullptr;

  // --- Mixing pool ---
  TH1D* hCollisionsPerDataframe = nullptr;
  TH1D* hMixBinOccupancyLeadP = nullptr;
  TH1D* hMixBinOccupancyLeadJet = nullptr;
  TH1D* hMixPoolOutcomeLeadP = nullptr;
  TH1D* hMixPoolOutcomeLeadJet = nullptr;
  TH1D* hZvtxAcceptance = nullptr;

  void book(TFile* outputFile);
};

// Labels for the integrity violation counter. Order defines the bin index.
static const std::vector<const char*> kIntegrityLabels = {
  "index < 0",
  "index >= nCollisions",
  "index not sorted",
  "leadP rows per collision > 1",
  "collision tree missing",
  "indexed tree missing",
  "non-finite value seen"};

// Labels for the collision-content pattern counter (bitmask: 1=jet, 2=leadP, 4=V0).
static const std::vector<const char*> kContentLabels = {
  "empty", "jet only", "leadP only", "jet+leadP",
  "V0 only", "jet+V0", "leadP+V0", "jet+leadP+V0"};

void ForensicsHistograms::book(TFile* outputFile)
{
  const int nColumns = static_cast<int>(kColumns.size());

  TDirectory* dTruncation = outputFile->mkdir("Truncation");
  TDirectory* dDuplicates = outputFile->mkdir("Duplicates");
  TDirectory* dCounters = outputFile->mkdir("Counters");
  TDirectory* dFingerprints = outputFile->mkdir("Fingerprints");
  TDirectory* dIntegrity = outputFile->mkdir("Integrity");
  TDirectory* dMixingPool = outputFile->mkdir("MixingPool");

  // --- Per-column histograms ---
  for (const ColumnSpec& column : kColumns) {
    dTruncation->cd();
    hTrailingZeros.push_back(new TH1D(
      Form("hTrailingZeros_%s", column.label),
      Form("%s: trailing zero mantissa bits;N_{trailing zeros};entries", column.label),
      25, -0.5, 24.5));
    hTrailingZerosVsBinade.push_back(new TH2D(
      Form("hTrailingZerosVsBinade_%s", column.label),
      Form("%s: truncation vs magnitude;#lfloor log_{2}|value| #rfloor;N_{trailing zeros}", column.label),
      60, -30.5, 29.5, 25, -0.5, 24.5));

    dDuplicates->cd();
    hMultiplicity.push_back(new TH1D(
      Form("hMultiplicity_%s", column.label),
      Form("%s: entries sharing one stored value;multiplicity;distinct values", column.label),
      200, 0.5, 200.5));
    hLogSpacing.push_back(new TH1D(
      Form("hLogSpacing_%s", column.label),
      Form("%s: gap between adjacent distinct values;log_{10}(#Delta);gaps", column.label),
      280, -9.0, 5.0));
    hSpacingOverGridStep.push_back(new TH1D(
      Form("hSpacingOverGridStep_%s", column.label),
      Form("%s: gap in units of the local grid step;log_{10}(#Delta / q);gaps", column.label),
      300, -0.5, 5.5));
    hDuplicateRowGap.push_back(new TH1D(
      Form("hDuplicateRowGap_%s", column.label),
      Form("%s: row separation of value-sharing pairs;log_{10}(|#Delta row|);pairs", column.label),
      140, 0.0, 7.0));
    pMultiplicityVsValue.push_back(new TProfile(
      Form("pMultiplicityVsValue_%s", column.label),
      Form("%s: mean multiplicity vs value;value;<multiplicity>", column.label),
      200, column.axisLow, column.axisHigh));
  }

  // --- Additive per-column counters ---
  dCounters->cd();
  auto bookColumnCounter = [nColumns](const char* name, const char* title) {
    TH1D* histogram = new TH1D(name, title, nColumns, -0.5, nColumns - 0.5);
    for (int i = 0; i < nColumns; ++i) {
      histogram->GetXaxis()->SetBinLabel(i + 1, kColumns[i].label);
    }
    return histogram;
  };
  hEntriesPerColumn = bookColumnCounter("hEntriesPerColumn", "Rows read;column;rows");
  hDistinctPerColumn = bookColumnCounter("hDistinctPerColumn", "Distinct values, summed over dataframes;column;distinct values");
  hObservedDuplicatePairs = bookColumnCounter("hObservedDuplicatePairs", "Observed value-sharing pairs;column;pairs");
  hExpectedDuplicatePairs = bookColumnCounter("hExpectedDuplicatePairs", "Expected value-sharing pairs;column;pairs");
  hNonFinitePerColumn = bookColumnCounter("hNonFinitePerColumn", "Non-finite values;column;entries");
  hConstantColumnDataframes = bookColumnCounter("hConstantColumnDataframes", "Dataframes where the column was constant;column;dataframes");

  // Dataframe-size binning shared by the two scaling counters. Log-spaced so a
  // single set of bins covers small and large dataframes.
  const int nSizeBins = 30;
  hObservedPairsVsDfSize = new TH1D("hObservedPairsVsDfSize",
                                    "Observed Zvtx pairs vs dataframe size;log_{10}(collisions in dataframe);pairs",
                                    nSizeBins, 0.0, 6.0);
  hExpectedPairsVsDfSize = new TH1D("hExpectedPairsVsDfSize",
                                    "Expected Zvtx pairs vs dataframe size;log_{10}(collisions in dataframe);pairs",
                                    nSizeBins, 0.0, 6.0);

  // --- Fingerprints ---
  dFingerprints->cd();
  auto bookFingerprintCounter = [](const char* name, const char* title) {
    TH1D* histogram = new TH1D(name, title, kNFingerprintColumns, 0.5, kNFingerprintColumns + 0.5);
    for (int k = 0; k < kNFingerprintColumns; ++k) {
      histogram->GetXaxis()->SetBinLabel(k + 1, Form("+%s", kColumns[k].label));
    }
    return histogram;
  };
  hObservedFingerprintPairs = bookFingerprintCounter(
    "hObservedFingerprintPairs", "Collision pairs matching the first k columns;columns used;pairs");
  hExpectedFingerprintPairs = bookFingerprintCounter(
    "hExpectedFingerprintPairs", "Chance expectation under column independence (lower bound);columns used;pairs");
  hFingerprintGroupSize = new TH1D("hFingerprintGroupSize",
                                   "Collisions sharing a full fingerprint;group size;groups",
                                   50, 1.5, 51.5);
  hFingerprintRowGap = new TH1D("hFingerprintRowGap",
                                "Row separation of full-fingerprint pairs;|#Delta row| within dataframe;pairs",
                                200, -0.5, 199.5);

  // --- Integrity ---
  dIntegrity->cd();
  hIntegrityViolations = new TH1D("hIntegrityViolations",
                                  "Integrity violations;check;occurrences",
                                  static_cast<int>(kIntegrityLabels.size()), -0.5,
                                  kIntegrityLabels.size() - 0.5);
  for (size_t i = 0; i < kIntegrityLabels.size(); ++i) {
    hIntegrityViolations->GetXaxis()->SetBinLabel(static_cast<int>(i) + 1, kIntegrityLabels[i]);
  }
  hRowsPerCollisionLeadP = new TH1D("hRowsPerCollisionLeadP",
                                    "Leading-particle rows per collision;rows;collisions", 11, -0.5, 10.5);
  hRowsPerCollisionJet = new TH1D("hRowsPerCollisionJet",
                                  "Jet rows per collision;rows;collisions", 101, -0.5, 100.5);
  hRowsPerCollisionV0 = new TH1D("hRowsPerCollisionV0",
                                 "V0 rows per collision;rows;collisions", 101, -0.5, 100.5);
  hCollisionContentPattern = new TH1D("hCollisionContentPattern",
                                      "Which tables reference each collision;content;collisions",
                                      static_cast<int>(kContentLabels.size()), -0.5,
                                      kContentLabels.size() - 0.5);
  for (size_t i = 0; i < kContentLabels.size(); ++i) {
    hCollisionContentPattern->GetXaxis()->SetBinLabel(static_cast<int>(i) + 1, kContentLabels[i]);
  }
  hIndexRangeExcess = new TH2D("hIndexRangeExcess",
                               "Out-of-range index values;table (0=jet, 1=leadP, 2=V0);signed distance outside [0, nCollisions)",
                               3, -0.5, 2.5, 201, -100.5, 100.5);

  // --- Mixing pool ---
  dMixingPool->cd();
  hCollisionsPerDataframe = new TH1D("hCollisionsPerDataframe",
                                     "Collisions per dataframe;collisions;dataframes", 200, 0.0, 20000.0);
  hMixBinOccupancyLeadP = new TH1D("hMixBinOccupancyLeadP",
                                   "Collisions per (Z_{vtx}, p_{T}, cent) bin, leadP proxy;occupancy;bins",
                                   200, 0.5, 200.5);
  hMixBinOccupancyLeadJet = new TH1D("hMixBinOccupancyLeadJet",
                                     "Collisions per (Z_{vtx}, p_{T}, cent) bin, leading-jet proxy;occupancy;bins",
                                     200, 0.5, 200.5);
  auto bookPoolOutcome = [](const char* name, const char* title) {
    TH1D* histogram = new TH1D(name, title, 3, -0.5, 2.5);
    histogram->GetXaxis()->SetBinLabel(1, "no proxy");
    histogram->GetXaxis()->SetBinLabel(2, "proxy, alone in bin");
    histogram->GetXaxis()->SetBinLabel(3, "proxy, partner available");
    return histogram;
  };
  hMixPoolOutcomeLeadP = bookPoolOutcome("hMixPoolOutcomeLeadP",
                                         "Mixing prospects per collision, leadP proxy;outcome;collisions");
  hMixPoolOutcomeLeadJet = bookPoolOutcome("hMixPoolOutcomeLeadJet",
                                           "Mixing prospects per collision, leading-jet proxy;outcome;collisions");
  hZvtxAcceptance = new TH1D("hZvtxAcceptance",
                             "Collisions inside the |Z_{vtx}| < 10 cm filter;0 = rejected, 1 = accepted;collisions",
                             2, -0.5, 1.5);

  outputFile->cd();
}

// ============================================================================
// SECTION 4 -- Per-dataframe column analysis
// ============================================================================

/// \brief Truncation, spacing and duplicate structure of one column in one dataframe.
/// \param values     the column, in table row order
/// \param columnIndex position of this column in kColumns, used for counter bins
/// \param histos     filled in place
static void analyseColumn(const std::vector<float>& values, int columnIndex, ForensicsHistograms& histos)
{
  const int nValues = static_cast<int>(values.size());
  if (nValues == 0) {
    return;
  }

  // --- Truncation depth. Measured per dataframe: with thousands of entries the
  // minimum is saturated, and every later step needs the grid step it implies.
  int minTrailingZeros = 24;
  int nNonFinite = 0;
  std::vector<std::pair<float, int>> sorted; // (value, original row index)
  sorted.reserve(nValues);

  for (int row = 0; row < nValues; ++row) {
    const float value = values[row];
    if (!std::isfinite(value)) {
      ++nNonFinite;
      continue;
    }
    sorted.emplace_back(value, row);
    const int trailingZeros = countTrailingZeroMantissaBits(value);
    if (trailingZeros < 0) {
      continue; // zero mantissa carries no truncation information
    }
    minTrailingZeros = std::min(minTrailingZeros, trailingZeros);
    histos.hTrailingZeros[columnIndex]->Fill(trailingZeros);
    histos.hTrailingZerosVsBinade[columnIndex]->Fill(std::ilogb(std::fabs(value)), trailingZeros);
  }

  histos.hNonFinitePerColumn->Fill(columnIndex, nNonFinite);
  if (nNonFinite > 0) {
    histos.hIntegrityViolations->Fill(6, nNonFinite); // "non-finite value seen"
  }
  if (sorted.empty()) {
    return;
  }
  const int keptBits = 23 - minTrailingZeros;

  std::sort(sorted.begin(), sorted.end(),
            [](const std::pair<float, int>& a, const std::pair<float, int>& b) { return a.first < b.first; });

  const int nFinite = static_cast<int>(sorted.size());
  histos.hEntriesPerColumn->Fill(columnIndex, nFinite);

  // --- Empirical density, used only for the birthday expectation. Bins are
  // many orders of magnitude wider than the grid step, so the density is
  // locally flat on the scale that matters.
  const double valueLow = sorted.front().first;
  const double valueHigh = sorted.back().first;
  const bool columnIsConstant = !(valueHigh > valueLow);
  if (columnIsConstant) {
    histos.hConstantColumnDataframes->Fill(columnIndex);
  }

  const int nDensityBins = 256;
  std::vector<double> densityCounts(nDensityBins, 0.0);
  const double densitySpan = columnIsConstant ? 1.0 : (valueHigh - valueLow);
  const double densityBinWidth = densitySpan / nDensityBins;
  if (!columnIsConstant) {
    for (const auto& entry : sorted) {
      int bin = static_cast<int>((entry.first - valueLow) / densitySpan * nDensityBins);
      bin = std::min(std::max(bin, 0), nDensityBins - 1);
      densityCounts[bin] += 1.0;
    }
  }

  // --- Walk the sorted values, grouping exact repeats.
  long long observedPairs = 0;
  long long nDistinct = 0;
  double expectedPairs = 0.0;

  for (int i = 0; i < nFinite;) {
    int j = i + 1;
    while (j < nFinite && sorted[j].first == sorted[i].first) {
      ++j;
    }
    const int multiplicity = j - i;
    ++nDistinct;
    observedPairs += static_cast<long long>(multiplicity) * (multiplicity - 1) / 2;

    histos.hMultiplicity[columnIndex]->Fill(multiplicity);
    histos.pMultiplicityVsValue[columnIndex]->Fill(sorted[i].first, multiplicity);

    // Row separation of every value-sharing pair inside this group.
    if (multiplicity > 1) {
      for (int a = i; a < j; ++a) {
        for (int b = a + 1; b < j; ++b) {
          const int rowGap = std::abs(sorted[a].second - sorted[b].second);
          if (rowGap > 0) {
            histos.hDuplicateRowGap[columnIndex]->Fill(std::log10(static_cast<double>(rowGap)));
          }
        }
      }
    }

    // Birthday contribution of this group:
    //   E[pairs] = C(N,2) * integral p^2 q, estimated as the sample mean of
    //   p(x_i) q(x_i) times C(N,2)/N, accumulated entry by entry.
    if (!columnIsConstant) {
      int densityBin = static_cast<int>((sorted[i].first - valueLow) / densitySpan * nDensityBins);
      densityBin = std::min(std::max(densityBin, 0), nDensityBins - 1);
      const double density = densityCounts[densityBin] / (nFinite * densityBinWidth);
      const double step = localGridStep(sorted[i].first, keptBits);
      expectedPairs += 0.5 * (nFinite - 1) * density * step * multiplicity;
    }

    // Gap to the next distinct value, raw and in units of the local grid step.
    if (j < nFinite) {
      const double gap = static_cast<double>(sorted[j].first) - static_cast<double>(sorted[i].first);
      if (gap > 0.0) {
        histos.hLogSpacing[columnIndex]->Fill(std::log10(gap));
        const double step = localGridStep(sorted[i].first, keptBits);
        if (step > 0.0) {
          histos.hSpacingOverGridStep[columnIndex]->Fill(std::log10(gap / step));
        }
      }
    }
    i = j;
  }

  // A constant column shares every pair by construction; report it as such
  // rather than letting the density model produce a meaningless number.
  if (columnIsConstant) {
    expectedPairs = 0.5 * static_cast<double>(nFinite) * (nFinite - 1);
  }

  histos.hDistinctPerColumn->Fill(columnIndex, nDistinct);
  histos.hObservedDuplicatePairs->Fill(columnIndex, observedPairs);
  histos.hExpectedDuplicatePairs->Fill(columnIndex, expectedPairs);

  // Scaling test, Zvtx only: both counters share a dataframe-size binning, so
  // the finalized ratio shows whether the quadratic birthday law holds.
  if (columnIndex == 0 && nFinite > 1) {
    const double sizeBin = std::log10(static_cast<double>(nFinite));
    histos.hObservedPairsVsDfSize->Fill(sizeBin, observedPairs);
    histos.hExpectedPairsVsDfSize->Fill(sizeBin, expectedPairs);
  }
}

// ============================================================================
// SECTION 5 -- Cumulative fingerprint test
// ============================================================================

/// \brief Counts collision pairs agreeing bit-for-bit on the first k stored
///        columns, for every k from 1 to kNFingerprintColumns.
/// \param columnValues one vector per collision column, all of equal length
/// \note  Restricted to a single dataframe on purpose: SameKindPair never
///        crosses a dataframe, so only within-dataframe repeats can reach the
///        mixing pool.
static void analyseFingerprints(const std::vector<std::vector<float>>& columnValues,
                                const std::vector<bool>& columnPresent,
                                ForensicsHistograms& histos)
{
  if (columnValues.size() < static_cast<size_t>(kNFingerprintColumns)) {
    return;
  }
  const int nCollisions = static_cast<int>(columnValues[0].size());
  if (nCollisions < 2) {
    return;
  }
  // A missing branch would enter the key as a constant and inflate every deeper
  // level, so the cumulative scan stops at the last column actually present.
  int nUsableColumns = 0;
  while (nUsableColumns < kNFingerprintColumns && columnPresent[nUsableColumns]) {
    ++nUsableColumns;
  }
  if (nUsableColumns == 0) {
    return;
  }

  // Independence-based chance expectation. Accumulated per level as the product
  // of the per-column match probabilities measured in this same dataframe.
  // The centrality columns are mutually correlated, so this product is a LOWER
  // BOUND on the true chance rate, never an upper bound.
  std::vector<double> perColumnMatchProbability(kNFingerprintColumns, 0.0);

  for (int k = 0; k < nUsableColumns; ++k) {
    std::unordered_map<uint32_t, int> valueCounts;
    valueCounts.reserve(nCollisions * 2);
    for (const float value : columnValues[k]) {
      ++valueCounts[floatBits(value)];
    }
    long long matchingPairs = 0;
    for (const auto& entry : valueCounts) {
      matchingPairs += static_cast<long long>(entry.second) * (entry.second - 1) / 2;
    }
    const double totalPairs = 0.5 * static_cast<double>(nCollisions) * (nCollisions - 1);
    perColumnMatchProbability[k] = (totalPairs > 0.0) ? matchingPairs / totalPairs : 0.0;
  }

  // Cumulative match counting. The key at level k is the concatenation of the
  // first k+1 column bit patterns.
  std::map<std::vector<uint32_t>, std::vector<int>> groupsByFingerprint;
  double independenceProbability = 1.0;

  for (int k = 0; k < nUsableColumns; ++k) {
    groupsByFingerprint.clear();
    for (int collision = 0; collision < nCollisions; ++collision) {
      std::vector<uint32_t> fingerprint;
      fingerprint.reserve(k + 1);
      for (int column = 0; column <= k; ++column) {
        fingerprint.push_back(floatBits(columnValues[column][collision]));
      }
      groupsByFingerprint[fingerprint].push_back(collision);
    }

    long long matchingPairs = 0;
    for (const auto& group : groupsByFingerprint) {
      const int groupSize = static_cast<int>(group.second.size());
      if (groupSize < 2) {
        continue;
      }
      matchingPairs += static_cast<long long>(groupSize) * (groupSize - 1) / 2;

      // Only the deepest level is worth detailing: those are the candidate
      // duplicated rows.
      if (k == nUsableColumns - 1) {
        histos.hFingerprintGroupSize->Fill(groupSize);
        for (int a = 0; a < groupSize; ++a) {
          for (int b = a + 1; b < groupSize; ++b) {
            histos.hFingerprintRowGap->Fill(std::abs(group.second[b] - group.second[a]));
          }
        }
      }
    }

    independenceProbability *= perColumnMatchProbability[k];
    const double totalPairs = 0.5 * static_cast<double>(nCollisions) * (nCollisions - 1);
    histos.hObservedFingerprintPairs->Fill(k + 1, static_cast<double>(matchingPairs));
    histos.hExpectedFingerprintPairs->Fill(k + 1, independenceProbability * totalPairs);
  }
}

// ============================================================================
// SECTION 6 -- Integrity and mixing pool
// ============================================================================

/// \brief Per-collision row counts collected from one indexed table.
struct IndexScanResult {
  std::vector<int> rowsPerCollision;
  std::vector<float> firstValuePerCollision; //! first proxy value seen, for the mixing pool
  bool treePresent = false;
};

/// \brief Reads an indexed table, validating its index column against the
///        collision count and recording per-collision row multiplicities.
/// \param proxyBranch optional branch whose maximum per collision is kept
///                    (jet or leading-particle pT); pass nullptr to skip.
static IndexScanResult scanIndexedTree(TDirectory* dataframe, const char* treeName,
                                       int nCollisions, int tableIndex,
                                       const char* proxyBranch, ForensicsHistograms& histos)
{
  IndexScanResult result;
  result.rowsPerCollision.assign(nCollisions, 0);
  result.firstValuePerCollision.assign(nCollisions, -1.f);

  TTree* tree = nullptr;
  dataframe->GetObject(treeName, tree);
  if (tree == nullptr) {
    histos.hIntegrityViolations->Fill(5); // "indexed tree missing"
    return result;
  }
  result.treePresent = true;

  tree->SetBranchStatus("*", 0);
  if (tree->GetBranch(kIndexBranch) == nullptr) {
    histos.hIntegrityViolations->Fill(5);
    return result;
  }
  tree->SetBranchStatus(kIndexBranch, 1);
  int collisionIndex = 0;
  tree->SetBranchAddress(kIndexBranch, &collisionIndex);

  float proxyValue = 0.f;
  const bool haveProxy = (proxyBranch != nullptr) && (tree->GetBranch(proxyBranch) != nullptr);
  if (haveProxy) {
    tree->SetBranchStatus(proxyBranch, 1);
    tree->SetBranchAddress(proxyBranch, &proxyValue);
  }

  int previousIndex = -1;
  const Long64_t nRows = tree->GetEntries();
  for (Long64_t row = 0; row < nRows; ++row) {
    tree->GetEntry(row);

    if (collisionIndex < 0) {
      histos.hIntegrityViolations->Fill(0);
      histos.hIndexRangeExcess->Fill(tableIndex, collisionIndex);
      continue;
    }
    if (collisionIndex >= nCollisions) {
      histos.hIntegrityViolations->Fill(1);
      histos.hIndexRangeExcess->Fill(tableIndex, collisionIndex - (nCollisions - 1));
      continue;
    }
    if (collisionIndex < previousIndex) {
      histos.hIntegrityViolations->Fill(2);
    }
    previousIndex = collisionIndex;

    ++result.rowsPerCollision[collisionIndex];
    if (haveProxy && proxyValue > result.firstValuePerCollision[collisionIndex]) {
      result.firstValuePerCollision[collisionIndex] = proxyValue;
    }
  }
  tree->ResetBranchAddresses();
  return result;
}

/// \brief Bin occupancy of the (Zvtx, proxy pT, centrality) mixing grid, and the
///        resulting per-collision mixing prospects.
static void analyseMixingPool(const std::vector<float>& zvtx,
                              const std::vector<float>& centrality,
                              const std::vector<float>& proxyPt,
                              TH1D* occupancyHistogram, TH1D* outcomeHistogram)
{
  const int nCollisions = static_cast<int>(zvtx.size());
  std::unordered_map<long long, int> occupancy;
  std::vector<long long> binOfCollision(nCollisions, -1);

  for (int collision = 0; collision < nCollisions; ++collision) {
    if (proxyPt[collision] < 0.f) {
      continue; // no proxy in this collision: it never enters the pool
    }
    const int binZ = findUniformBin(kAxisPVzLow, kAxisPVzHigh, kAxisPVzBins, zvtx[collision]);
    const int binPt = findVariableBin(kAxisJetPtEdges, proxyPt[collision]);
    const int binCent = findVariableBin(kAxisCentralityEdges, centrality[collision]);
    if (binZ < 0 || binPt < 0 || binCent < 0) {
      continue; // outside the mixing grid, exactly as the binning policy would
    }
    const long long nPtBins = static_cast<long long>(kAxisJetPtEdges.size());
    const long long nCentBins = static_cast<long long>(kAxisCentralityEdges.size());
    const long long globalBin = (static_cast<long long>(binZ) * nPtBins + binPt) * nCentBins + binCent;
    binOfCollision[collision] = globalBin;
    ++occupancy[globalBin];
  }

  for (const auto& entry : occupancy) {
    occupancyHistogram->Fill(entry.second);
  }
  for (int collision = 0; collision < nCollisions; ++collision) {
    if (binOfCollision[collision] < 0) {
      outcomeHistogram->Fill(0); // no proxy, or off-grid
    } else if (occupancy[binOfCollision[collision]] < 2) {
      outcomeHistogram->Fill(1); // proxy exists but nothing to mix with
    } else {
      outcomeHistogram->Fill(2); // at least one partner in the same bin
    }
  }
}

// ============================================================================
// SECTION 7 -- Dataframe and file drivers
// ============================================================================

/// \brief Runs every module on one DF_* directory.
static void processDataframe(TDirectory* dataframe, ForensicsHistograms& histos)
{
  TTree* collisionTree = nullptr;
  dataframe->GetObject("O2ringcollision", collisionTree);
  if (collisionTree == nullptr) {
    histos.hIntegrityViolations->Fill(4); // "collision tree missing"
    return;
  }

  // --- Read every collision column once.
  const int nCollisionColumns = kNFingerprintColumns;
  std::vector<std::vector<float>> collisionColumns(nCollisionColumns);
  std::vector<float> readBuffer(nCollisionColumns, 0.f);

  collisionTree->SetBranchStatus("*", 0);
  std::vector<bool> columnPresent(nCollisionColumns, false);
  for (int column = 0; column < nCollisionColumns; ++column) {
    if (collisionTree->GetBranch(kColumns[column].branchName) != nullptr) {
      collisionTree->SetBranchStatus(kColumns[column].branchName, 1);
      collisionTree->SetBranchAddress(kColumns[column].branchName, &readBuffer[column]);
      columnPresent[column] = true;
    }
  }

  const Long64_t nCollisions = collisionTree->GetEntries();
  for (int column = 0; column < nCollisionColumns; ++column) {
    collisionColumns[column].reserve(nCollisions);
  }
  for (Long64_t row = 0; row < nCollisions; ++row) {
    collisionTree->GetEntry(row);
    for (int column = 0; column < nCollisionColumns; ++column) {
      collisionColumns[column].push_back(columnPresent[column] ? readBuffer[column] : 0.f);
    }
  }
  collisionTree->ResetBranchAddresses();

  histos.hCollisionsPerDataframe->Fill(static_cast<double>(nCollisions));
  for (const float z : collisionColumns[0]) {
    histos.hZvtxAcceptance->Fill(std::fabs(z) < 10.0f ? 1 : 0);
  }

  // --- Column-wise truncation and duplicate structure.
  for (int column = 0; column < nCollisionColumns; ++column) {
    if (columnPresent[column]) {
      analyseColumn(collisionColumns[column], column, histos);
    }
  }

  // --- Cumulative fingerprint test.
  analyseFingerprints(collisionColumns, columnPresent, histos);

  // --- Indexed tables: integrity plus the proxy values the mixing pool needs.
  const int nCollisionsInt = static_cast<int>(nCollisions);
  const IndexScanResult jetScan =
    scanIndexedTree(dataframe, "O2ringjet", nCollisionsInt, 0, "fJetPt", histos);
  const IndexScanResult leadPScan =
    scanIndexedTree(dataframe, "O2ringleadp", nCollisionsInt, 1, "fLeadParticlePt", histos);
  const IndexScanResult v0Scan =
    scanIndexedTree(dataframe, "O2ringlav0", nCollisionsInt, 2, nullptr, histos);

  for (int collision = 0; collision < nCollisionsInt; ++collision) {
    const int nJets = jetScan.rowsPerCollision[collision];
    const int nLeadP = leadPScan.rowsPerCollision[collision];
    const int nV0s = v0Scan.rowsPerCollision[collision];

    histos.hRowsPerCollisionJet->Fill(nJets);
    histos.hRowsPerCollisionLeadP->Fill(nLeadP);
    histos.hRowsPerCollisionV0->Fill(nV0s);

    // The producer writes at most one leading-particle row per collision, and
    // the consumer's mixing pool relies on that by taking rows.begin().
    if (nLeadP > 1) {
      histos.hIntegrityViolations->Fill(3);
    }

    const int contentPattern = (nJets > 0 ? 1 : 0) | (nLeadP > 0 ? 2 : 0) | (nV0s > 0 ? 4 : 0);
    histos.hCollisionContentPattern->Fill(contentPattern);
  }

  // --- Jet and leading-particle columns, analysed as their own populations.
  for (int column = nCollisionColumns; column < static_cast<int>(kColumns.size()); ++column) {
    TTree* tree = nullptr;
    dataframe->GetObject(kColumns[column].treeName, tree);
    if (tree == nullptr) {
      continue;
    }
    tree->SetBranchStatus("*", 0);
    if (tree->GetBranch(kColumns[column].branchName) == nullptr) {
      continue;
    }
    tree->SetBranchStatus(kColumns[column].branchName, 1);
    float value = 0.f;
    tree->SetBranchAddress(kColumns[column].branchName, &value);

    const Long64_t nRows = tree->GetEntries();
    std::vector<float> values;
    values.reserve(nRows);
    for (Long64_t row = 0; row < nRows; ++row) {
      tree->GetEntry(row);
      values.push_back(value);
    }
    tree->ResetBranchAddresses();
    analyseColumn(values, column, histos);
  }

  // --- Mixing-pool occupancy, once per proxy type.
  analyseMixingPool(collisionColumns[0], collisionColumns[1], leadPScan.firstValuePerCollision,
                    histos.hMixBinOccupancyLeadP, histos.hMixPoolOutcomeLeadP);
  analyseMixingPool(collisionColumns[0], collisionColumns[1], jetScan.firstValuePerCollision,
                    histos.hMixBinOccupancyLeadJet, histos.hMixPoolOutcomeLeadJet);
}

/// \brief Opens one AO2D and runs every DF_* directory through processDataframe.
/// \return number of dataframes processed, or -1 if the file could not be opened.
static int processFile(const std::string& path, ForensicsHistograms& histos)
{
  TFile* input = TFile::Open(path.c_str(), "READ");
  if (input == nullptr || input->IsZombie()) {
    std::cerr << "  could not open " << path << std::endl;
    delete input;
    return -1;
  }

  int nDataframes = 0;
  TIter nextKey(input->GetListOfKeys());
  TKey* key = nullptr;
  while ((key = static_cast<TKey*>(nextKey())) != nullptr) {
    const TString name = key->GetName();
    if (!name.BeginsWith("DF_")) {
      continue;
    }
    TDirectory* dataframe = nullptr;
    input->GetObject(name, dataframe);
    if (dataframe == nullptr) {
      continue;
    }
    processDataframe(dataframe, histos);
    ++nDataframes;
  }

  input->Close();
  delete input;
  return nDataframes;
}

// ============================================================================
// SECTION 8 -- Worker and finalize entry points
// ============================================================================

static int runWorker(const std::string& manifestPath, const std::string& outputPath)
{
  std::ifstream manifest(manifestPath);
  if (!manifest.is_open()) {
    std::cerr << "Could not open manifest: " << manifestPath << std::endl;
    return 1;
  }
  std::vector<std::string> files;
  std::string line;
  while (std::getline(manifest, line)) {
    if (!line.empty() && line[0] != '#') {
      files.push_back(line);
    }
  }
  manifest.close();

  if (files.empty()) {
    std::cerr << "Manifest is empty: " << manifestPath << std::endl;
    return 1;
  }

  TFile* output = TFile::Open(outputPath.c_str(), "RECREATE");
  if (output == nullptr || output->IsZombie()) {
    std::cerr << "Could not create output: " << outputPath << std::endl;
    delete output;
    return 1;
  }

  ForensicsHistograms histos;
  histos.book(output);

  int nDataframes = 0;
  int nFailures = 0;
  for (const std::string& file : files) {
    const int result = processFile(file, histos);
    if (result < 0) {
      ++nFailures;
    } else {
      nDataframes += result;
    }
  }

  output->Write();
  output->Close();
  delete output;

  std::cout << "Processed " << files.size() << " files, " << nDataframes
            << " dataframes, " << nFailures << " unreadable." << std::endl;
  return (nFailures == static_cast<int>(files.size())) ? 1 : 0;
}

/// \brief Divides one additive counter pair into a ratio histogram.
/// \note  Binomial errors are wrong here (numerator and denominator are not
///        counts of the same trials), so uncertainties are left off deliberately.
static void writeRatio(TDirectory* directory, const char* numeratorName,
                       const char* denominatorName, const char* ratioName,
                       const char* ratioTitle)
{
  TH1D* numerator = nullptr;
  TH1D* denominator = nullptr;
  directory->GetObject(numeratorName, numerator);
  directory->GetObject(denominatorName, denominator);
  if (numerator == nullptr || denominator == nullptr) {
    std::cerr << "  missing pair for " << ratioName << ", skipped." << std::endl;
    return;
  }

  TH1D* ratio = static_cast<TH1D*>(numerator->Clone(ratioName));
  ratio->SetTitle(ratioTitle);
  ratio->SetDirectory(directory);
  ratio->Reset();

  for (int bin = 1; bin <= numerator->GetNbinsX(); ++bin) {
    const double denominatorValue = denominator->GetBinContent(bin);
    if (denominatorValue > 0.0) {
      ratio->SetBinContent(bin, numerator->GetBinContent(bin) / denominatorValue);
    }
    ratio->SetBinError(bin, 0.0);
  }
  directory->cd();
  ratio->Write(ratioName, TObject::kOverwrite);
}

static int runFinalize(const std::string& mergedPath)
{
  TFile* merged = TFile::Open(mergedPath.c_str(), "UPDATE");
  if (merged == nullptr || merged->IsZombie()) {
    std::cerr << "Could not open merged file: " << mergedPath << std::endl;
    delete merged;
    return 1;
  }

  TDirectory* counters = merged->GetDirectory("Counters");
  if (counters != nullptr) {
    writeRatio(counters, "hObservedDuplicatePairs", "hExpectedDuplicatePairs",
               "hDuplicatePairRatio",
               "Observed / expected value-sharing pairs;column;ratio");
    writeRatio(counters, "hDistinctPerColumn", "hEntriesPerColumn",
               "hDistinctFraction",
               "Fraction of rows carrying a distinct value;column;distinct / rows");
    writeRatio(counters, "hObservedPairsVsDfSize", "hExpectedPairsVsDfSize",
               "hDuplicatePairRatioVsDfSize",
               "Observed / expected Zvtx pairs vs dataframe size;log_{10}(collisions);ratio");
  }

  TDirectory* fingerprints = merged->GetDirectory("Fingerprints");
  if (fingerprints != nullptr) {
    writeRatio(fingerprints, "hObservedFingerprintPairs", "hExpectedFingerprintPairs",
               "hFingerprintPairRatio",
               "Observed / chance-expected matches;columns used;ratio");
  }

  merged->Close();
  delete merged;
  std::cout << "Finalized " << mergedPath << std::endl;
  return 0;
}

// ============================================================================
// SECTION 9 -- main
// ============================================================================

int main(int argc, char** argv)
{
  if (argc == 3 && std::strcmp(argv[1], "--finalize") == 0) {
    return runFinalize(argv[2]);
  }
  if (argc == 3) {
    return runWorker(argv[1], argv[2]);
  }

  std::cerr << "Usage:\n"
            << "  " << argv[0] << " <manifest.txt> <output.root>\n"
            << "  " << argv[0] << " --finalize <merged.root>\n";
  return 1;
}