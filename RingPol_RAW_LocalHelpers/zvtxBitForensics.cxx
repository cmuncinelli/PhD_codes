// ============================================================================
// zvtxBitForensics.cxx
//
// Bit-level forensics on float columns stored in AO2D derived tables, with
// the Zvtx column of O2ringcollision as the default target.
//
// What this answers:
//   (1) How many mantissa bits actually survive in the stored floats, i.e.
//       recover the O2 storage truncation empirically instead of trusting a
//       remembered hex mask.
//   (2) Whether the observed number of bit-identical repeats is consistent
//       with the birthday-paradox expectation on that measured grid, or
//       whether there is a genuine excess pointing at duplicated rows.
//   (3) Whether duplicate-valued entries sit at adjacent row indices, which
//       separates a harmless numerical coincidence from a split vertex or a
//       double-written collision.
//   (4) Whether any two collisions coincide on ALL stored float columns at
//       once, which cannot happen by chance and is therefore a definitive
//       duplicated-row test.
//
// Usage from the ROOT CLI (no includes or O2 headers needed on the user side):
//   root -l AO2D_1.root
//   .x zvtxBitForensics.cxx
//   .x zvtxBitForensics.cxx("AO2D_1.root")
//   .x zvtxBitForensics.cxx("AO2D_1.root", "O2ringjet", "fJetPt")
//
// ASCII only throughout.
// ============================================================================

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TKey.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TString.h>
#include <TTree.h>
#include <TVirtualPad.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

// ----------------------------------------------------------------------------
// Low-level IEEE-754 helpers
// ----------------------------------------------------------------------------

// Reinterpret the four raw bytes of a float as an unsigned 32-bit word.
// memcpy is used rather than a pointer cast so that no strict-aliasing
// assumption is made; this is the only portable way to inspect the bits.
static uint32_t floatBits(float value)
{
  uint32_t bits = 0;
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

// Count how many low-order bits of the 23-bit mantissa field are exactly zero.
// If the AO2D writer truncated the value by ANDing with a mask that clears the
// N lowest mantissa bits, then EVERY stored value has at least N trailing
// zeros, and the true N is the minimum of this count over the whole sample
// (attained by whichever event happened to use all of its surviving precision).
// Returns -1 for the uninformative case of a zero mantissa (exact zero, or an
// exact power of two), which must be excluded from the minimum.
static int countTrailingZeroMantissaBits(float value)
{
  const uint32_t bits = floatBits(value);
  const uint32_t mantissa = bits & 0x007FFFFFu; // low 23 bits = mantissa field
  if (mantissa == 0u) {
    return -1; // no information: cannot distinguish truncation from a clean value
  }
  int trailingZeros = 0;
  uint32_t work = mantissa;
  while ((work & 1u) == 0u) {
    work >>= 1;
    ++trailingZeros;
  }
  return trailingZeros;
}

// Local grid step (in the units of the variable) for a float that carries
// keptBits informative mantissa bits.
// For 2^e <= |x| < 2^(e+1) the spacing of representable float32 values is
// 2^(e-23); clearing the low (23 - keptBits) bits multiplies that spacing by
// 2^(23-keptBits), giving a step of 2^(e - keptBits).
static double localGridStep(double value, int keptBits)
{
  if (value == 0.0 || !std::isfinite(value)) {
    return 0.0;
  }
  const int exponent = std::ilogb(std::fabs(value));
  return std::ldexp(1.0, exponent - keptBits);
}

// ----------------------------------------------------------------------------
// Per-sample analysis result
// ----------------------------------------------------------------------------

struct ColumnStats {
  Long64_t nEntries = 0;      // total rows seen
  Long64_t nDistinct = 0;     // number of distinct stored values
  Long64_t nDuplicatePairs = 0; // sum over value groups of C(multiplicity, 2)
  double expectedPairs = 0.0; // birthday-paradox expectation on the measured grid
  int minTrailingZeros = 24;  // measured truncation depth (bits cleared)
  int keptBits = 23;          // 23 - minTrailingZeros
  double minSpacing = 0.0;    // smallest gap between distinct values
  Long64_t nNonFinite = 0;    // NaN / inf guard
};

// ----------------------------------------------------------------------------
// Core analysis of one float sample.
//
// values : the stored floats
// rows   : the row index each value came from, parallel to "values". Used for
//          the duplicate-adjacency test. Pass an empty vector to skip that test.
// label  : prefix used for histogram names and printout
// verbose: print the per-value table for the first few distinct values
// ----------------------------------------------------------------------------
static ColumnStats analyseFloatColumn(const std::vector<float>& values,
                                      const std::vector<Long64_t>& rows,
                                      const char* label,
                                      bool makeHistograms,
                                      bool verbose)
{
  ColumnStats stats;
  stats.nEntries = static_cast<Long64_t>(values.size());
  if (stats.nEntries == 0) {
    return stats;
  }
  const bool haveRows = (rows.size() == values.size());

  // --- Pass 1: bit statistics -----------------------------------------------
  // Measure the truncation depth before anything else, because the grid step
  // used by every later step depends on it.
  TH1I* hTrailingZeros = nullptr;
  TH2D* hTrailingZerosVsBinade = nullptr;
  if (makeHistograms) {
    hTrailingZeros = new TH1I(Form("hTrailingZeros_%s", label),
                              Form("%s: trailing zero mantissa bits;N_{trailing zeros};entries", label),
                              25, -0.5, 24.5);
    // A flat band in this 2D plot means a relative (mantissa-mask) truncation:
    // constant number of significant bits, absolute step growing with |value|.
    // A band that slopes with the binade means an absolute fixed-step rounding.
    hTrailingZerosVsBinade = new TH2D(Form("hTrailingZerosVsBinade_%s", label),
                                      Form("%s: truncation vs magnitude;#lfloor log_{2}|value| #rfloor;N_{trailing zeros}", label),
                                      40, -20.5, 19.5, 25, -0.5, 24.5);
  }

  for (const float value : values) {
    if (!std::isfinite(value)) {
      ++stats.nNonFinite;
      continue;
    }
    const int trailingZeros = countTrailingZeroMantissaBits(value);
    if (trailingZeros < 0) {
      continue; // zero mantissa: uninformative, must not drag the minimum down
    }
    if (trailingZeros < stats.minTrailingZeros) {
      stats.minTrailingZeros = trailingZeros;
    }
    if (makeHistograms) {
      hTrailingZeros->Fill(trailingZeros);
      if (value != 0.f) {
        hTrailingZerosVsBinade->Fill(std::ilogb(std::fabs(value)), trailingZeros);
      }
    }
  }
  stats.keptBits = 23 - stats.minTrailingZeros;

  // --- Pass 2: sort, keeping the row index attached -------------------------
  // Sorting by value groups identical values together; the row index rides
  // along so that the adjacency test below can ask where the duplicates came
  // from in the original table ordering.
  std::vector<std::pair<float, Long64_t>> sorted;
  sorted.reserve(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    if (!std::isfinite(values[i])) {
      continue;
    }
    sorted.emplace_back(values[i], haveRows ? rows[i] : static_cast<Long64_t>(i));
  }
  std::sort(sorted.begin(), sorted.end(),
            [](const std::pair<float, Long64_t>& a, const std::pair<float, Long64_t>& b) {
              return a.first < b.first;
            });

  // --- Empirical density, needed for the birthday-paradox expectation -------
  // A coarse histogram (bins far wider than the grid step) is a perfectly good
  // estimator of the smooth density p(z); it only has to be locally flat on
  // the scale of the grid step, which it is by many orders of magnitude.
  const double densityBinWidth = 0.25; // cm for Zvtx; adjust if reusing for pT
  const double densityLow = sorted.front().first - densityBinWidth;
  const double densityHigh = sorted.back().first + densityBinWidth;
  const int nDensityBins = std::max(1, static_cast<int>((densityHigh - densityLow) / densityBinWidth));
  TH1D hDensity(Form("hDensity_%s", label), "", nDensityBins, densityLow, densityHigh);
  hDensity.SetDirectory(nullptr); // scratch object: keep it out of the ROOT file
  for (const auto& entry : sorted) {
    hDensity.Fill(entry.first);
  }
  const double densityNorm = static_cast<double>(sorted.size()) * hDensity.GetXaxis()->GetBinWidth(1);

  // --- Spacing / multiplicity / duplicate histograms ------------------------
  TH1I* hMultiplicity = nullptr;
  TH1D* hLogSpacing = nullptr;
  TH1D* hSpacingOverStep = nullptr;
  TH1D* hDuplicateRowGap = nullptr;
  TProfile* pMultiplicityVsValue = nullptr;
  if (makeHistograms) {
    hMultiplicity = new TH1I(Form("hMultiplicity_%s", label),
                             Form("%s: multiplicity of identical values;entries per distinct value;number of distinct values", label),
                             200, 0.5, 200.5);
    // Log-binned in the sense that the FILLED quantity is log10(dz): this keeps
    // uniform bins while covering the many decades produced by the sparse tails,
    // where consecutive distinct values are separated by many empty grid cells.
    hLogSpacing = new TH1D(Form("hLogSpacing_%s", label),
                           Form("%s: spacing between distinct values;log_{10}(#Delta z / cm);number of spacings", label),
                           280, -9.0, 5.0);
    // The physically meaningful version: spacing divided by the LOCAL grid step.
    // On a true quantisation grid this piles up at integers, so the picket-fence
    // of peaks at 1, 2, 3, ... is far stronger evidence of a grid than any single
    // minimum-spacing scalar could be.
    hSpacingOverStep = new TH1D(Form("hSpacingOverStep_%s", label),
                                Form("%s: spacing in units of the local grid step;log_{10}(#Delta z / q(z));number of spacings", label),
                                300, -0.5, 5.5);
    // Row-index separation of duplicate-valued entries. Uniform-ish means
    // coincidence; a spike at small gaps means split or double-written rows.
    hDuplicateRowGap = new TH1D(Form("hDuplicateRowGap_%s", label),
                                Form("%s: row gap between duplicate-valued entries;log_{10}(|#Delta row|);pairs", label),
                                140, 0.0, 7.0);
    pMultiplicityVsValue = new TProfile(Form("pMultiplicityVsValue_%s", label),
                                        Form("%s: mean multiplicity vs value;value;<multiplicity>", label),
                                        120, densityLow, densityHigh);
  }

  stats.minSpacing = std::numeric_limits<double>::max();

  if (verbose) {
    std::cout << "\nFirst 30 distinct values of " << label << ":\n";
    std::cout << "            value       multiplicity           Delta          Delta/step\n";
    std::cout << "-------------------------------------------------------------------------\n";
  }

  for (size_t i = 0; i < sorted.size();) {
    // Find the end of this run of identical values.
    size_t j = i + 1;
    while (j < sorted.size() && sorted[j].first == sorted[i].first) {
      ++j;
    }
    const size_t multiplicity = j - i;
    ++stats.nDistinct;
    stats.nDuplicatePairs += static_cast<Long64_t>(multiplicity) * static_cast<Long64_t>(multiplicity - 1) / 2;

    if (makeHistograms) {
      hMultiplicity->Fill(static_cast<double>(multiplicity));
      pMultiplicityVsValue->Fill(sorted[i].first, static_cast<double>(multiplicity));
      // Row-adjacency test: for every duplicate pair inside this group, record
      // how far apart the two rows were in the original table.
      if (multiplicity > 1) {
        for (size_t a = i; a < j; ++a) {
          for (size_t b = a + 1; b < j; ++b) {
            const Long64_t gap = std::llabs(sorted[a].second - sorted[b].second);
            if (gap > 0) {
              hDuplicateRowGap->Fill(std::log10(static_cast<double>(gap)));
            }
          }
        }
      }
    }

    // Birthday-paradox contribution of every entry in this group.
    // E[pairs] = C(N,2) * integral p^2 q dz, and the integral is estimated as
    // the sample mean of p(z_i) q(z_i), so each entry contributes
    // (N-1)/2 * p(z_i) q(z_i) / 1, summed over all N entries.
    const double density = hDensity.GetBinContent(hDensity.FindBin(sorted[i].first)) / densityNorm;
    const double step = localGridStep(sorted[i].first, stats.keptBits);
    stats.expectedPairs += 0.5 * static_cast<double>(sorted.size() - 1) * density * step * static_cast<double>(multiplicity);

    // Spacing to the next distinct value.
    double delta = 0.0;
    if (j < sorted.size()) {
      delta = static_cast<double>(sorted[j].first) - static_cast<double>(sorted[i].first);
      if (delta > 0.0) {
        if (makeHistograms) {
          hLogSpacing->Fill(std::log10(delta));
          if (step > 0.0) {
            hSpacingOverStep->Fill(std::log10(delta / step));
          }
        }
        if (delta < stats.minSpacing) {
          stats.minSpacing = delta;
        }
      }
    }

    if (verbose && stats.nDistinct <= 30) {
      std::cout << std::fixed << std::setprecision(9)
                << std::setw(17) << sorted[i].first
                << std::setw(15) << multiplicity;
      if (j < sorted.size()) {
        std::cout << std::setw(20) << delta;
        if (step > 0.0) {
          std::cout << std::setw(20) << std::setprecision(3) << (delta / step);
        }
      }
      std::cout << std::endl;
    }

    i = j;
  }

  if (stats.minSpacing == std::numeric_limits<double>::max()) {
    stats.minSpacing = 0.0;
  }
  return stats;
}

// ----------------------------------------------------------------------------
// Print one stats block.
// ----------------------------------------------------------------------------
static void printStats(const ColumnStats& stats, const char* label)
{
  std::cout << "\n=========================================================\n";
  std::cout << "SUMMARY: " << label << "\n";
  std::cout << "=========================================================\n";
  std::cout << std::defaultfloat;
  std::cout << "Entries                    : " << stats.nEntries << "\n";
  if (stats.nNonFinite > 0) {
    std::cout << "NON-FINITE entries         : " << stats.nNonFinite << "   <-- investigate\n";
  }
  std::cout << "Distinct values            : " << stats.nDistinct << "\n";
  std::cout << "Mantissa bits cleared      : " << stats.minTrailingZeros
            << "   (kept: " << stats.keptBits << " of 23)\n";
  std::cout << "Relative grid step         : " << std::scientific << std::setprecision(3)
            << std::ldexp(1.0, -stats.keptBits) << "\n";
  std::cout << "Smallest observed spacing  : " << stats.minSpacing << "\n";
  std::cout << std::defaultfloat;
  std::cout << "Duplicate pairs observed   : " << stats.nDuplicatePairs << "\n";
  std::cout << "Duplicate pairs expected   : " << std::fixed << std::setprecision(1)
            << stats.expectedPairs << "   (birthday paradox on the measured grid)\n";
  if (stats.expectedPairs > 0.0) {
    const double ratio = static_cast<double>(stats.nDuplicatePairs) / stats.expectedPairs;
    std::cout << "Observed / expected        : " << std::setprecision(3) << ratio;
    if (ratio > 3.0) {
      std::cout << "   <-- EXCESS: suspect duplicated rows";
    } else if (ratio < 0.33) {
      std::cout << "   <-- DEFICIT: density model or truncation misidentified";
    } else {
      std::cout << "   <-- consistent with pure numerical coincidence";
    }
    std::cout << "\n";
  }
  std::cout << "=========================================================\n";
  std::cout << std::defaultfloat;
}

// ----------------------------------------------------------------------------
// Definitive duplicated-row test.
//
// Two distinct collisions coinciding on a single float column is common (see
// the birthday-paradox numbers above). Coinciding SIMULTANEOUSLY on several
// independent float columns is not: the probability is the product of the
// individual ones, so any hit here is a genuinely duplicated row rather than a
// coincidence. This needs no modelling assumption whatsoever.
// ----------------------------------------------------------------------------
static void reportTupleDuplicates(const std::vector<std::vector<uint32_t>>& tuples,
                                  const std::vector<Long64_t>& rows,
                                  const char* label)
{
  if (tuples.empty()) {
    return;
  }
  std::map<std::vector<uint32_t>, std::vector<Long64_t>> seen;
  for (size_t i = 0; i < tuples.size(); ++i) {
    seen[tuples[i]].push_back(rows[i]);
  }
  Long64_t nDuplicateGroups = 0;
  Long64_t nDuplicateRows = 0;
  std::cout << "\n--- Full-tuple duplicate test (" << label << ") ---\n";
  for (const auto& entry : seen) {
    if (entry.second.size() > 1) {
      ++nDuplicateGroups;
      nDuplicateRows += static_cast<Long64_t>(entry.second.size());
      if (nDuplicateGroups <= 20) {
        std::cout << "  identical on all columns, rows:";
        for (const Long64_t row : entry.second) {
          std::cout << " " << row;
        }
        std::cout << "\n";
      }
    }
  }
  if (nDuplicateGroups == 0) {
    std::cout << "  none. No two rows agree on every stored float column.\n";
    std::cout << "  => the single-column repeats are numerical coincidence, not duplicated rows.\n";
  } else {
    std::cout << "  " << nDuplicateGroups << " duplicated groups covering "
              << nDuplicateRows << " rows.\n";
    std::cout << "  => these are GENUINE duplicated collision rows. Event mixing can pair\n";
    std::cout << "     them with each other and produce self-correlation despite distinct indices.\n";
  }
}

// ----------------------------------------------------------------------------
// Main entry point.
// ----------------------------------------------------------------------------
void zvtxBitForensics(const char* fileName = "",
                      const char* treeName = "O2ringcollision",
                      const char* branchName = "fZvtx",
                      bool doTupleTest = true)
{
  // Use the already-open file when the macro is launched as "root AO2D_1.root".
  TFile* file = nullptr;
  bool ownFile = false;
  if (fileName != nullptr && std::strlen(fileName) > 0) {
    file = TFile::Open(fileName, "READ");
    ownFile = true;
  } else {
    file = gFile;
  }
  if (file == nullptr || file->IsZombie()) {
    std::cerr << "Could not open input file." << std::endl;
    return;
  }

  // Every AO2D stores its tables inside one directory per dataframe. Event
  // mixing operates WITHIN a dataframe, so the per-DF numbers are the ones
  // that matter for the mixing question; the pooled numbers are the ones that
  // matter for characterising the storage truncation.
  std::vector<float> allValues;
  std::vector<Long64_t> allRows;
  std::vector<std::vector<uint32_t>> allTuples;

  // Extra collision columns used by the full-tuple test. Missing branches are
  // skipped silently so the macro still works on other tables.
  const char* extraBranches[] = {"fCentFT0M", "fCentFT0C", "fCentFV0A", "fInteractionRate"};
  const int nExtraBranches = 4;

  int nDataframes = 0;
  Long64_t globalRow = 0;

  TIter nextKey(file->GetListOfKeys());
  TKey* key = nullptr;
  while ((key = static_cast<TKey*>(nextKey())) != nullptr) {
    const TString dirName = key->GetName();
    if (!dirName.BeginsWith("DF_")) {
      continue;
    }
    TDirectory* dir = nullptr;
    file->GetObject(dirName, dir);
    if (dir == nullptr) {
      continue;
    }
    TTree* tree = nullptr;
    dir->GetObject(treeName, tree);
    if (tree == nullptr) {
      continue;
    }
    ++nDataframes;

    // Read only what is needed: with many columns this is a large I/O saving.
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus(branchName, 1);
    float value = 0.f;
    tree->SetBranchAddress(branchName, &value);

    // Bind the extra columns when they exist in this tree.
    std::vector<float> extraValues(nExtraBranches, 0.f);
    std::vector<bool> extraPresent(nExtraBranches, false);
    if (doTupleTest) {
      for (int b = 0; b < nExtraBranches; ++b) {
        if (tree->GetBranch(extraBranches[b]) != nullptr) {
          tree->SetBranchStatus(extraBranches[b], 1);
          tree->SetBranchAddress(extraBranches[b], &extraValues[b]);
          extraPresent[b] = true;
        }
      }
    }

    const Long64_t nEntries = tree->GetEntries();
    std::vector<float> dfValues;
    std::vector<Long64_t> dfRows;
    dfValues.reserve(nEntries);
    dfRows.reserve(nEntries);

    for (Long64_t i = 0; i < nEntries; ++i) {
      tree->GetEntry(i);
      dfValues.push_back(value);
      dfRows.push_back(i); // row index WITHIN this dataframe
      allValues.push_back(value);
      allRows.push_back(globalRow + i);
      if (doTupleTest) {
        std::vector<uint32_t> tuple;
        tuple.push_back(floatBits(value));
        for (int b = 0; b < nExtraBranches; ++b) {
          if (extraPresent[b]) {
            tuple.push_back(floatBits(extraValues[b]));
          }
        }
        if (tuple.size() > 1) { // only meaningful with at least two columns
          allTuples.push_back(tuple);
        }
      }
    }
    globalRow += nEntries;

    // Per-dataframe report, no histograms (they would multiply without limit).
    const ColumnStats dfStats = analyseFloatColumn(dfValues, dfRows, dirName.Data(), false, false);
    std::cout << "DF " << std::setw(3) << nDataframes
              << "  " << dirName
              << "  entries " << std::setw(8) << dfStats.nEntries
              << "  distinct " << std::setw(8) << dfStats.nDistinct
              << "  dup pairs obs " << std::setw(7) << dfStats.nDuplicatePairs
              << "  exp " << std::fixed << std::setprecision(1) << std::setw(9) << dfStats.expectedPairs
              << std::defaultfloat << std::endl;
  }

  if (nDataframes == 0) {
    std::cerr << "No DF_* directories containing tree '" << treeName << "' were found." << std::endl;
    if (ownFile) {
      file->Close();
    }
    return;
  }

  std::cout << "\nScanned " << nDataframes << " dataframes." << std::endl;

  // Pooled analysis with the full histogram set.
  const ColumnStats pooled = analyseFloatColumn(allValues, allRows, "pooled", true, true);
  printStats(pooled, branchName);

  // Definitive duplicated-row test on the collision tuple.
  if (doTupleTest && !allTuples.empty()) {
    reportTupleDuplicates(allTuples, allRows, treeName);
  }

  // ------------------------------------------------------------------------
  // Draw. Four canvases, each answering one of the questions above.
  // ------------------------------------------------------------------------
  TCanvas* cBits = new TCanvas("cBits", "truncation depth", 1200, 500);
  cBits->Divide(2, 1);
  cBits->cd(1);
  gPad->SetLogy();
  if (gDirectory->Get("hTrailingZeros_pooled") != nullptr) {
    gDirectory->Get("hTrailingZeros_pooled")->Draw();
  }
  cBits->cd(2);
  gPad->SetLogz();
  if (gDirectory->Get("hTrailingZerosVsBinade_pooled") != nullptr) {
    gDirectory->Get("hTrailingZerosVsBinade_pooled")->Draw("COLZ");
  }

  TCanvas* cSpacing = new TCanvas("cSpacing", "spacing structure", 1200, 500);
  cSpacing->Divide(2, 1);
  cSpacing->cd(1);
  gPad->SetLogy();
  if (gDirectory->Get("hLogSpacing_pooled") != nullptr) {
    gDirectory->Get("hLogSpacing_pooled")->Draw();
  }
  cSpacing->cd(2);
  gPad->SetLogy();
  if (gDirectory->Get("hSpacingOverStep_pooled") != nullptr) {
    gDirectory->Get("hSpacingOverStep_pooled")->Draw();
  }

  TCanvas* cDuplicates = new TCanvas("cDuplicates", "duplicate structure", 1200, 500);
  cDuplicates->Divide(2, 1);
  cDuplicates->cd(1);
  gPad->SetLogy();
  if (gDirectory->Get("hMultiplicity_pooled") != nullptr) {
    gDirectory->Get("hMultiplicity_pooled")->Draw();
  }
  cDuplicates->cd(2);
  gPad->SetLogy();
  if (gDirectory->Get("hDuplicateRowGap_pooled") != nullptr) {
    gDirectory->Get("hDuplicateRowGap_pooled")->Draw();
  }

  TCanvas* cShape = new TCanvas("cShape", "multiplicity vs value", 700, 500);
  if (gDirectory->Get("pMultiplicityVsValue_pooled") != nullptr) {
    gDirectory->Get("pMultiplicityVsValue_pooled")->Draw();
  }

  // Deliberately not closing the file when it was already open in the CLI.
  if (ownFile) {
    // Histograms live in gDirectory; detach before closing so they survive.
    std::cout << "\n(Input file was opened by the macro and is left open so the plots stay valid.)\n";
  }
}
