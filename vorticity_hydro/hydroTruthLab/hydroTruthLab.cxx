// hydroTruthLab.cxx
//
// Truth-level ring observable and its longitudinal projection R_z, computed directly from the per-cell Lambda polarization that iSS writes for
// each hydro event: no sampling, no decays (way simpler than lambda_pol_toy_model.cxx's implementation and way more stable).
// One output file per sample (bullet or no_bullet).
//
// Usage:
//   ./hydroTruthLab.exe <sampleDir> <outputFile> [--option=value ...]      (--help lists every option and its default)
//   Normally compiled and launched by runHydroTruthLab.sh, which holds every configurable.
//
// Needs: ROOT (Hist, Tree, RIO, GenVector), the iSS files <sampleDir>/<i>/results/<inputName>, and the jet file for --trigger=jet and mixing.
//
// Design, definitions, error treatment and the output reference: README.md in this folder.

#include <Math/GenVector/VectorUtil.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <TAxis.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TNamed.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using ROOT::Math::PxPyPzEVector;
using ROOT::Math::XYZVector;

// ======================================================================================================================================
// !!! READ THIS BEFORE TOUCHING THE POLARIZATION COLUMNS !!!
// The iSS columns labelled "S^t S^x S^y S^z" (and their "_LRF" twins) already hold the POLARIZATION P^mu = S^mu/(1/2), NOT the spin vector S^mu.
// iSS commit a8e1b80 (2021-03-27) moved the factor 2 into the code but kept the S labels, and Vitor's build (/home/vribeiro/vorticity/iSS-pol in
// jarvis4) has that normalization: prefactor = 1./(2.*mass), with no 0.5 in front of the vorticity terms. Multiplying by 2 here doubles every
// result, which is exactly what the old lambda_pol_toy_model.cxx did. Evidence, and the check to repeat on regenerated files: README.md, "IMPORTANT".
// ======================================================================================================================================
constexpr double kColumnToPolarization = 1.0;

constexpr double kNoValue = -999.;            // "not filled" guard value. Same as ALICE standard values
constexpr int kNumberOfColumnsISSOutput = 12; // y, pT, phi, dN/(pTdpTdphidy), P^mu (lab), P^mu (LRF)
constexpr double kNodeTolerance = 1e-5;       // grid nodes are printed with 7 significant digits, so they sit far closer than this to an exactly uniform grid
constexpr double kJetAngleTolerance = 1e-4;   // rad. Same printing argument, applied to the jet file
constexpr double kEdgeTolerance = 1e-6;       // in units of the cell width, for the cuts-on-cell-edges check
constexpr int kProgressEvery = 25;            // convenience counter for prints
constexpr int kDistributionBins = 100;        // per-event distributions, whose range is taken from the values themselves
constexpr int kZetaBins = 100;                // distributions of cos(zeta) and sin(zeta), over [-1, 1]

// Closure residuals are histogrammed as log10 of a relative residual: their scale (the print precision) spans decades and is not known beforehand
constexpr int kLogResidualBins = 160;
constexpr double kLogResidualLow = -16.;
constexpr double kLogResidualHigh = 0.;

constexpr double kPi = TMath::Pi();
constexpr double kTwoPi = TMath::TwoPi();

// Column order of the iSS per-cell files: "# y pT[GeV] phi dN/(pTdpTdphidy)[GeV^-2] S^t S^x S^y S^z S^t_LRF S^x_LRF S^y_LRF S^z_LRF"
enum CellColumn { kColY = 0,
                  kColPt,
                  kColPhi,
                  kColMult,
                  kColPolT,
                  kColPolX,
                  kColPolY,
                  kColPolZ,
                  kColPolStarT,
                  kColPolStarX,
                  kColPolStarY,
                  kColPolStarZ };
static_assert(kColPolStarZ + 1 == kNumberOfColumnsISSOutput, "CellColumn and kNumberOfColumnsISSOutput disagree");

// Per-cell observables of the true trigger (README, "Observables"). The index is only used to fill: every object is booked by name.
enum RingObservable { kRing = 0,
                      kRingZ,
                      kRingZT,
                      kHelicityTerm,
                      kRingZLab,
                      kRingPerp,
                      kNzSq,
                      kNzSqSinSqTheta,
                      kNRingObservables };

// Polarization components of PolVsAngles/: lab P and rest-frame P*. In the trigger-aligned coordinate system there, X and Y stand for x' (along t)
// and y' = z x x' (not the consumer's PrimeJet system, which PolarizationMaps/ uses)
// (For a better reference of what these maps mean, just take a look at my consumer's bookings on O2Physics' lambdaJetPolarizationIonsDerived.cxx)
enum PolComponent { kPolLabX = 0,
                    kPolLabY,
                    kPolLabZ,
                    kPolStarX,
                    kPolStarY,
                    kPolStarZ,
                    kNPolComponents };

// PolarizationMaps/, as RawPolarization/PolMaps/ in the consumer (lambdaJetPolarizationIonsDerived.cxx): the three P* components of a
// coordinate system (x, y, z of Lab, x', y', z' of PrimeJet and PrimeV0), and the two ring scalars mapped on the PrimeJet planes
enum MapComponent { kMapX = 0,
                    kMapY,
                    kMapZ,
                    kNMapComponents };
enum MapRing { kMapRing = 0,
               kMapRingZ,
               kNMapRings };

// Pick the jet direction from random_parameters.dat or use phi_RP = 0 in Trento's initialization
// (following what Vitor seems to have used for the PRC paper)
enum class TriggerMode { kJet,
                         kReactionPlane };
// Either we do a full 3D smudging of the jet direction, or a transverse plane only smudge:
enum class SmudgeMode { kCone,
                        kPhi };

// Axis labels shared by many titles
const std::string kDeltaPhiLabel = "#Delta#varphi = #varphi - #varphi_{t} (rad)";
const std::string kPhiLabel = "#varphi - #Psi_{RP} (rad)";
const std::string kPtLabel = "#it{p}_{T} (GeV/#it{c})";
const std::string kYLabel = "#it{y}";
const std::string kDeltaLabel = "#delta (rad)";
const std::string kPxLabel = "#it{p}_{x} (GeV/#it{c})";
const std::string kPyLabel = "#it{p}_{y} (GeV/#it{c})";
const std::string kPzLabel = "#it{p}_{z} (GeV/#it{c})";
const std::string kPxPrimeJetLabel = "#it{p}_{x'Jet} (GeV/#it{c})";
const std::string kPyPrimeJetLabel = "#it{p}_{y'Jet} (GeV/#it{c})";
const std::string kPzPrimeJetLabel = "#it{p}_{z'Jet} (GeV/#it{c})";

// ======================================================================================================================================
// Configuration
// ======================================================================================================================================

struct AnalysisConfig {
  std::string sampleDir;
  std::string outputFile;
  std::string jetFile; // empty when not given
  TriggerMode triggerMode = TriggerMode::kJet;
  int firstEvent = 0;
  int nEvents = 250;
  std::string inputName = "Smu_dpTdphidy_Thermal_rapidity_3122.dat";
  // Kinematic cuts of Vitor's phenomenology paper (PhysRevC.109.014905), applied on cell centres. There, 0.5 < pT < 1.5 GeV/c keeps the region
  // most affected by the thermalized jet: the upper limit reduces the contamination from the anisotropic expansion (v2), and the lower one follows
  // the experimental acceptance and the small Lambda yield below it. |y| < 0.5 is the usual experimental acceptance. We repeat these cuts as a
  // closure test of this code: does it recover Vitor's results?
  double ptMin = 0.5;
  double ptMax = 1.5;
  double absYMax = 0.5;
  double lambdaMass = 1.115683;
  int deltaScanSteps = 100;
  int mixDepth = 0;
  double smudgeSigma = 0.;
  SmudgeMode smudgeMode = SmudgeMode::kCone;
  int smudgeDraws = 1;
  std::uint32_t seed = 0; // drawn from the current time unless --seed is given (main, section 1). Recorded in Config/ either way
  int mapSubdivisions = 3; // PolarizationMaps/: each cell spread over n^3 points of its momentum volume, 1 = nodes only (README, "Polarization maps")
};

void printUsage()
{
  std::cout << "Usage: ./hydroTruthLab.exe <sampleDir> <outputFile> [--option=value ...]\n"
            << "  --jetFile=<path>         jet file, required for --trigger=jet and for mixing      (default: none)\n"
            << "  --trigger=jet|rp         phi_t = phi_J, or phi_t = 0 (reaction plane)             (default: jet)\n"
            << "  --firstEvent=<int>       first event folder                                       (default: 0)\n"
            << "  --nEvents=<int>          number of events                                         (default: 250)\n"
            << "  --inputName=<name>       file inside <sampleDir>/<i>/results/                     (default: Smu_dpTdphidy_Thermal_rapidity_3122.dat)\n"
            << "  --ptMin=<GeV/c>          lower pT cut, on cell centres                            (default: 0.5)\n"
            << "  --ptMax=<GeV/c>          upper pT cut, on cell centres                            (default: 1.5)\n"
            << "  --absYMax=<value>        |y| cut, on cell centres                                 (default: 0.5)\n"
            << "  --lambdaMass=<GeV/c^2>   mass used to build E and p_z from (y, pT)                (default: 1.115683)\n"
            << "  --deltaScanSteps=<int>   deltas of the scan, 0 = off. Even, divides the phi nodes (default: 100)\n"
            << "  --mixDepth=<int>         mixing partners per event, 0 = off (jet trigger only)    (default: 0)\n"
            << "  --smudgeSigma=<rad>      smudge width, 0 = off                                    (default: 0)\n"
            << "  --smudgeMode=cone|phi    isotropic cone around t, or azimuth only                 (default: cone)\n"
            << "  --smudgeDraws=<int>      smudged triggers per event                               (default: 1)\n"
            << "  --seed=<uint32>          base seed, each event seeded from (seed, eventId)        (default: from the current time)\n"
            << "                           Paired smudging needs the same seed for both samples: runHydroTruthLab.sh passes one per run\n"
            << "  --mapSubdivisions=<int>  PolarizationMaps/: cells spread over n^3 points, 1 = nodes  (default: 3)\n"
            << "See README.md for the definitions." << std::endl;
}

// The std::sto* functions throw with messages like "stoi", which say nothing about the option: these wrappers name it
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

std::uint32_t toUInt32(const std::string& key, const std::string& value)
{
  std::size_t parsedLength = 0;
  unsigned long long result = 0;
  try {
    result = std::stoull(value, &parsedLength);
  } catch (const std::exception&) {
    parsedLength = 0;
  }
  if (parsedLength == 0 || parsedLength != value.size() || value[0] == '-' || result > std::numeric_limits<std::uint32_t>::max())
    throw std::runtime_error("--" + key + " expects an unsigned 32-bit integer, got '" + value + "'");
  return result;
}

// ======================================================================================================================================
// Small helpers
// ======================================================================================================================================

double squared(double value) { return value * value; }

/// \brief Wraps an angle into [low, low + 2 pi), up to rounding: an angle one ulp below low + 2 pi can come out one ulp below low. Fine for
///        ROOT fills, which put such a value in the underflow. Array indices use AxisSpec::periodicBin instead.
double wrapAngle(double angle, double low) { return angle - kTwoPi * std::floor((angle - low) / kTwoPi); }

/// \brief Input file of one event: <sampleDir>/<eventId>/results/<inputName> (README, "Data layout").
std::string eventInputPath(const AnalysisConfig& config, int eventId) { return config.sampleDir + "/" + std::to_string(eventId) + "/results/" + config.inputName; }

/// \brief Sorted event IDs as ranges, "28-29, 33-34, 49-199", for the log and Config/.
std::string formatEventRanges(const std::vector<int>& eventIds)
{
  std::string text;
  for (std::size_t first = 0; first < eventIds.size();) {
    std::size_t last = first;
    while (last + 1 < eventIds.size() && eventIds[last + 1] == eventIds[last] + 1)
      ++last;
    text += (text.empty() ? "" : ", ") + std::to_string(eventIds[first]) + (last > first ? "-" + std::to_string(eventIds[last]) : "");
    first = last + 1;
  }
  return text;
}

/// \brief Residual over its natural scale, for the closure histograms. A vanishing scale with a vanishing residual counts as exact.
double relativeResidual(double residual, double scale)
{
  if (scale > 0.)
    return std::abs(residual) / scale;
  return residual == 0. ? 0. : std::numeric_limits<double>::infinity();
}

/// \brief ROOT bin (0 = underflow, which also takes exact zeros) of log10(residual) on the closure axis.
int logResidualBin(double residual)
{
  if (!(residual > 0.))
    return 0;
  const double logResidual = std::log10(residual);
  if (logResidual < kLogResidualLow)
    return 0;
  if (logResidual >= kLogResidualHigh)
    return kLogResidualBins + 1;
  return 1 + static_cast<int>((logResidual - kLogResidualLow) / (kLogResidualHigh - kLogResidualLow) * kLogResidualBins);
}

/// \brief 0-based bin of cos(zeta) or sin(zeta) on [-1, 1]. Rounding can push a product of unit vectors just past +-1: such values go to the edge bins.
int zetaBin(double value)
{
  const int bin = std::floor((value + 1.) * 0.5 * kZetaBins);
  return std::clamp(bin, 0, kZetaBins - 1);
}

/// \brief Nodes of one grid axis: the distinct values of its column, sorted (30000 y entries --> 10 y nodes). Their number and spacing give the
///        shape of the grid without hardcoding it, and without relying on the order in which iSS writes the rows.
/// \note  Every repeat of a node is printed from the same double, so identical text gives bit-identical values: exact comparison is intended.
std::vector<double> determineGridShape(const std::vector<double>& columnValues)
{
  std::vector<double> nodes(columnValues);
  std::sort(nodes.begin(), nodes.end());
  nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
  return nodes;
}

/// \brief Index of a value on the uniform grid first + index * width, with a hard stop if it is not a grid node.
int gridIndex(double value, double first, double width, const std::string& columnName)
{
  const long index = std::lround((value - first) / width);
  if (std::abs(value - (first + index * width)) > kNodeTolerance)
    throw std::runtime_error("grid column " + columnName + " is not uniform (value " + std::to_string(value) + ")");
  return index;
}

/// \brief Whether a cut value sits on a cell edge of the uniform grid firstNode + i * width (nodes are cell centres).
bool isOnCellEdge(double cut, double firstNode, double width)
{
  const double edgeIndex = (cut - (firstNode - 0.5 * width)) / width;
  return std::abs(edgeIndex - std::round(edgeIndex)) < kEdgeTolerance;
}

struct AxisSpec {
  int nBins = 0;
  double low = 0.;
  double high = 0.;
  double width() const { return (high - low) / nBins; }
  double binCentre(int bin) const { return low + (bin + 0.5) * width(); } // 0-based bin
  /// 0-based bin, -1 outside [low, high)
  int findBin(double value) const
  {
    if (!(value >= low && value < high))
      return -1;
    return std::min(static_cast<int>((value - low) / width()), nBins - 1); // min() guards the rounding as value --> high
  }
  /// 0-based bin of an angle on an axis that spans 2 pi, wrapped as an integer: always a valid bin. Wrapping the angle first and then calling
  /// findBin() is not safe, since wrapAngle() can return a value one ulp outside the axis, and findBin() then returns -1
  int periodicBin(double angle) const
  {
    const long bin = std::floor((angle - low) / width());
    return ((bin % nBins) + nBins) % nBins;
  }
};

/// \brief Fills one event's per-bin ratios A_b/B_b into a pooled master (weight B_b, so the mean is sum(A)/sum(B)) and, if given, into an
///        equal-weight master (weight 1, so the mean is the event mean). README, "Errors".
/// \note  A bin with no cells in this event carries no information and is not filled.
void flushPerEventBins(TProfile* pooledMaster, TProfile* evtMeanMaster, const std::vector<double>& sumW, const std::vector<double>& sumWX)
{
  const TAxis* xAxis = pooledMaster->GetXaxis();
  for (std::size_t bin = 0; bin < sumW.size(); ++bin) {
    if (sumW[bin] <= 0.)
      continue;
    const double binMean = sumWX[bin] / sumW[bin];
    const double binCentre = xAxis->GetBinCenter(bin + 1);
    pooledMaster->Fill(binCentre, binMean, sumW[bin]);
    if (evtMeanMaster != nullptr)
      evtMeanMaster->Fill(binCentre, binMean);
  }
}

/// \brief One per-event value at abscissa x: weight eventWeight (the event's summed cell weight) in the pooled master, weight 1 in the equal-weight one.
void fillEventPoint(TProfile* pooledMaster, TProfile* evtMeanMaster, double x, double value, double eventWeight)
{
  pooledMaster->Fill(x, value, eventWeight);
  if (evtMeanMaster != nullptr)
    evtMeanMaster->Fill(x, value);
}

/// \brief Flat average over bins of the per-bin ratio, (1/N_filled) sum_b A_b/B_b. The grid binning fills every bin (see countEmptyBins).
double flatAverage(const std::vector<double>& sumW, const std::vector<double>& sumWX)
{
  double sumOfRatios = 0.;
  int nFilledBins = 0;
  for (std::size_t bin = 0; bin < sumW.size(); ++bin) {
    if (sumW[bin] > 0.) {
      sumOfRatios += sumWX[bin] / sumW[bin];
      ++nFilledBins;
    }
  }
  return nFilledBins > 0 ? sumOfRatios / nFilledBins : kNoValue;
}

int countEmptyBins(const std::vector<double>& sumW)
{
  return std::count_if(sumW.begin(), sumW.end(), [](double binSumW) { return !(binSumW > 0.); });
}

double sumOf(const std::vector<double>& values) { return std::accumulate(values.begin(), values.end(), 0.); }

/// \brief Histogram of per-event values, with its range taken from the values themselves (their scale is not known beforehand).
void bookDistribution(TDirectory* directory, const std::string& name, const std::string& title, const std::vector<double>& values)
{
  if (values.empty())
    return;
  const auto [minIt, maxIt] = std::minmax_element(values.begin(), values.end());
  const double padding = *maxIt > *minIt ? 0.05 * (*maxIt - *minIt) : std::max(1e-12, 0.05 * std::abs(*minIt));
  directory->cd();
  TH1D* distribution = new TH1D(name.c_str(), title.c_str(), kDistributionBins, *minIt - padding, *maxIt + padding);
  for (const double value : values)
    distribution->Fill(value);
}

/// \brief Clones an object into another directory under a new name (the Paper/ copies).
void cloneInto(TDirectory* directory, const TH1* source, const std::string& newName)
{
  directory->cd(); // Clone() may register the copy in the current directory before SetDirectory() moves it: keep that in the target
  auto* clone = static_cast<TH1*>(source->Clone(newName.c_str()));
  clone->SetDirectory(directory);
}

// ======================================================================================================================================
// Input file handling:
// ======================================================================================================================================

/// \brief All 12 columns of one iSS per-cell file, in file order.
struct CellTable {
  std::array<std::vector<double>, kNumberOfColumnsISSOutput> column;
  int size() const { return column[kColY].size(); }
};

/// \brief Reads one iSS per-cell file, with strict format checks (hard stop on any malformed line).
/// \note  The file is read into memory at once and parsed with strtod, several times faster than ifstream >> double on these ~5 MB files.
void readCellTable(const std::string& path, CellTable& table)
{
  std::ifstream inputStream(path, std::ios::binary);
  if (!inputStream)
    throw std::runtime_error("cannot open " + path);
  const std::string buffer((std::istreambuf_iterator<char>(inputStream)), std::istreambuf_iterator<char>());

  const int previousSize = table.size();
  for (std::vector<double>& columnValues : table.column) {
    columnValues.clear();
    columnValues.reserve(previousSize);
  }

  const char* cursor = buffer.c_str();
  const char* const bufferEnd = cursor + buffer.size();
  int lineNumber = 0;
  while (cursor < bufferEnd) {
    ++lineNumber;
    const char* lineEnd = static_cast<const char*>(std::memchr(cursor, '\n', bufferEnd - cursor));
    if (lineEnd == nullptr)
      lineEnd = bufferEnd;
    const char* firstChar = cursor;
    while (firstChar < lineEnd && std::isspace(static_cast<unsigned char>(*firstChar)))
      ++firstChar;
    if (firstChar == lineEnd || *firstChar == '#') { // header or blank line
      cursor = lineEnd + 1;
      continue;
    }
    const char* parsePos = firstChar;
    for (int iColumn = 0; iColumn < kNumberOfColumnsISSOutput; ++iColumn) {
      char* afterNumber = nullptr;
      const double value = std::strtod(parsePos, &afterNumber);
      // strtod also skips newlines, so a short line would silently borrow numbers from the next one: hence the lineEnd bound
      if (afterNumber == parsePos || afterNumber > lineEnd)
        throw std::runtime_error(path + ", line " + std::to_string(lineNumber) + ": expected " + std::to_string(kNumberOfColumnsISSOutput) + " numeric columns");
      table.column[iColumn].push_back(value);
      parsePos = afterNumber;
    }
    for (; parsePos < lineEnd; ++parsePos)
      if (!std::isspace(static_cast<unsigned char>(*parsePos)))
        throw std::runtime_error(path + ", line " + std::to_string(lineNumber) + ": more than " + std::to_string(kNumberOfColumnsISSOutput) + " columns");
    cursor = lineEnd + 1;
  }
  if (table.size() == 0)
    throw std::runtime_error(path + " has no data lines");
}

struct JetEntry {
  bool present = false;
  double phi = kNoValue;
};

/// \brief Event-independent grid: node indices, Lambda kinematics, cell measure and cut flags, per cell in file order.
struct CellGrid {
  int nY = 0, nPt = 0, nPhi = 0, nCells = 0;
  double yFirst = 0., ptFirst = 0., phiFirst = 0.;
  double yWidth = 0., ptWidth = 0., phiWidth = 0.; // phiWidth = 2 pi / nPhi exactly
  std::vector<double> y, pt, phi;                  // as printed, the reference for the event-by-event check
  std::vector<int> iY, iPt, iPhi;
  std::vector<double> px, py, pz, energy, pMag, pHatX, pHatY, pHatZ, sinSqTheta, cellMeasure;
  std::vector<XYZVector> boostToRest;
  std::vector<char> inPtCut, inYCut; // inside the paper's cuts (see AnalysisConfig). char instead of vector<bool>, whose proxy bits are slow in the cell loop
};

/// \brief Compact copy of the cells inside both of the paper's cuts (see AnalysisConfig), for the variant loops (delta-scan, mixing, smudging):
///        kinematics once, fields refreshed every event.
struct InCutCells {
  std::vector<int> cellIndex, phiIndex;
  std::vector<double> px, py, pz, pHatX, pHatY, pHatZ;      // event independent
  std::vector<double> weight, polStarX, polStarY, polStarZ; // refreshed every event
  std::vector<double> nTrueX, nTrueY, nTrueZ;               // ring direction of the true trigger, the reference of cos(zeta): refreshed every event
  int size() const { return cellIndex.size(); }
};

/// \brief Points of the smoothed momentum maps (PolarizationMaps/): every cell spread over n^3 points of its own (y, pT, phi) volume, each carrying
///        the field of the node (README, "Polarization maps"). Event independent, so built once, indexed [iCell * perCell + iPoint].
struct SubCellPoints {
  int perCell = 0;
  std::vector<double> px, py;
  std::vector<double> weightFraction; // d3N share of the point, (pT_point/pT_node)/n^3: the shares of one cell add up to 1
  std::vector<int> pxBin, pyBin, pzBin; // bins that do not depend on the trigger. pT needs none: every point stays in the pT bin of its cell
};

// ======================================================================================================================================
// Per-event accumulators (README, "Errors")
// ======================================================================================================================================

/// \brief One observable on one axis: cells are summed per event (A_b, B_b), then flushed once per event into the masters.
struct PerEventProfile1D {
  TProfile* pooledMaster = nullptr;
  TProfile* evtMeanMaster = nullptr; // optional
  std::vector<double> sumW;          // B_b, summed cell weights
  std::vector<double> sumWX;         // A_b, summed weight times observable

  void attach(TProfile* master)
  {
    pooledMaster = master;
    sumW.assign(master->GetNbinsX(), 0.);
    sumWX.assign(master->GetNbinsX(), 0.);
  }
  void attachEvtMean(TProfile* master) { evtMeanMaster = master; }
  void add(int bin, double weight, double value)
  {
    sumW[bin] += weight;
    sumWX[bin] += weight * value;
  }
  void flushAndReset()
  {
    flushPerEventBins(pooledMaster, evtMeanMaster, sumW, sumWX);
    std::fill(sumW.begin(), sumW.end(), 0.);
    std::fill(sumWX.begin(), sumWX.end(), 0.);
  }
};

/// \brief Same as PerEventProfile1D, on two axes (pooled master only).
struct PerEventProfile2D {
  TProfile2D* pooledMaster = nullptr;
  int nBinsX = 0;
  std::vector<double> sumW;
  std::vector<double> sumWX;

  void attach(TProfile2D* master)
  {
    pooledMaster = master;
    nBinsX = master->GetNbinsX();
    sumW.assign(nBinsX * master->GetNbinsY(), 0.);
    sumWX.assign(nBinsX * master->GetNbinsY(), 0.);
  }
  void add(int binX, int binY, double weight, double value)
  {
    const int index = binX + nBinsX * binY;
    sumW[index] += weight;
    sumWX[index] += weight * value;
  }
  void flushAndReset()
  {
    const TAxis* xAxis = pooledMaster->GetXaxis();
    const TAxis* yAxis = pooledMaster->GetYaxis();
    for (std::size_t index = 0; index < sumW.size(); ++index) {
      if (sumW[index] <= 0.) // no cells in this bin in this event: no information
        continue;
      const int binX = index % nBinsX;
      const int binY = index / nBinsX;
      pooledMaster->Fill(xAxis->GetBinCenter(binX + 1), yAxis->GetBinCenter(binY + 1), sumWX[index] / sumW[index], sumW[index]);
    }
    std::fill(sumW.begin(), sumW.end(), 0.);
    std::fill(sumWX.begin(), sumWX.end(), 0.);
  }
};

/// \brief A TH1D filled through a plain array (ROOT bin numbering, under- and overflow included) and flushed once at the end: keeps the
///        virtual Fill() out of the cell loop. Only the contents are meaningful (yields and closure counts): with the default Sumw2 their errors read 0.
/// \note  This type of optimization is very close to what we did for the PhD_codes/ToyModels/helicityEfficiencyToyModel.cxx macro. A lot of computing time
///        can be wasted just because of Fill() calls that could have been buffered into such a struct!
struct BufferedHist1D {
  TH1D* hist = nullptr;
  std::vector<double> content;

  void attach(TH1D* target)
  {
    hist = target;
    content.assign(target->GetNbinsX() + 2, 0.);
  }
  void add(int rootBin, double weight) { content[rootBin] += weight; }
  void flush() const
  {
    for (std::size_t rootBin = 0; rootBin < content.size(); ++rootBin)
      hist->SetBinContent(rootBin, content[rootBin]);
    hist->ResetStats();
  }
};

/// \brief Same as BufferedHist1D, on two axes (0-based in-range bins only).
struct BufferedHist2D {
  TH2D* hist = nullptr;
  int nBinsX = 0;
  std::vector<double> content;

  void attach(TH2D* target)
  {
    hist = target;
    nBinsX = target->GetNbinsX();
    content.assign(nBinsX * target->GetNbinsY(), 0.);
  }
  void add(int binX, int binY, double weight) { content[binX + nBinsX * binY] += weight; }
  void flush() const
  {
    for (std::size_t index = 0; index < content.size(); ++index)
      hist->SetBinContent(index % nBinsX + 1, index / nBinsX + 1, content[index]);
    hist->ResetStats();
  }
};

/// \brief Per-bin sums of the reduced set (README, "Delta-scan, mixing and smudging") for one trigger variant within one event, plus the
///        d3N-weighted distributions of cos(zeta) and sin(zeta) (bins of zetaBin(), not Delta phi bins).
struct ReducedAccumulators {
  std::vector<double> sumW, sumWRing, sumWRingZ, sumWNzSq, sumWCosZeta, sumWSinZeta;
  std::vector<double> cosZetaDistribution, sinZetaDistribution;

  void resize(int nBins)
  {
    sumW.assign(nBins, 0.);
    sumWRing.assign(nBins, 0.);
    sumWRingZ.assign(nBins, 0.);
    sumWNzSq.assign(nBins, 0.);
    sumWCosZeta.assign(nBins, 0.);
    sumWSinZeta.assign(nBins, 0.);
    cosZetaDistribution.assign(kZetaBins, 0.);
    sinZetaDistribution.assign(kZetaBins, 0.);
  }
  void reset()
  {
    std::fill(sumW.begin(), sumW.end(), 0.);
    std::fill(sumWRing.begin(), sumWRing.end(), 0.);
    std::fill(sumWRingZ.begin(), sumWRingZ.end(), 0.);
    std::fill(sumWNzSq.begin(), sumWNzSq.end(), 0.);
    std::fill(sumWCosZeta.begin(), sumWCosZeta.end(), 0.);
    std::fill(sumWSinZeta.begin(), sumWSinZeta.end(), 0.);
    std::fill(cosZetaDistribution.begin(), cosZetaDistribution.end(), 0.);
    std::fill(sinZetaDistribution.begin(), sinZetaDistribution.end(), 0.);
  }
};

/// \brief Integrated sums of one trigger axis alone.
struct IntegratedSums {
  double sumW = 0.;
  double sumWRing = 0.;
  double sumWRingZ = 0.;
  double sumWNzSq = 0.;
  double sumWCosZeta = 0.;
  double sumWSinZeta = 0.;
};

struct ReducedValues {
  double ringD3N = kNoValue;
  double ringFlat = kNoValue;
  double ringZD3N = kNoValue;
  double ringZFlat = kNoValue;
  double cosZetaD3N = kNoValue;
  double sinZetaD3N = kNoValue;
};

// ======================================================================================================================================
// Output objects
// ======================================================================================================================================

/// \brief Everything filled with the true trigger.
struct TrueTriggerOutputs {
  // Inputs/
  BufferedHist1D yieldVsPhi, yieldVsPt, yieldVsY;
  BufferedHist2D yieldVsDeltaPhiVsPt;
  TH1D* phiTrigger = nullptr;
  // Ring/ (cuts: pT and y)
  std::array<PerEventProfile1D, kNRingObservables> obsVsDeltaPhi;
  std::array<TProfile*, kNRingObservables> obsD3N{};
  std::array<TProfile*, kNRingObservables> obsD3NEvtMean{}; // Ring and RingZ only
  TProfile* ringFlat = nullptr;
  TProfile* ringFlatEvtMean = nullptr;
  TProfile* ringZFlat = nullptr;
  TProfile* ringZFlatEvtMean = nullptr;
  PerEventProfile2D ringVsDeltaPhiVsPt, ringZVsDeltaPhiVsPt; // y cut only, all pT
  PerEventProfile2D ringVsDeltaPhiVsY, ringZVsDeltaPhiVsY;   // pT cut only
  PerEventProfile1D ringVsPt, ringZVsPt;                     // y cut only
  PerEventProfile1D ringVsY, ringZVsY;                       // pT cut only
  PerEventProfile1D ringVsDeltaPhiNoCuts, ringZVsDeltaPhiNoCuts;
  // PolVsAngles/
  std::array<PerEventProfile1D, kNPolComponents> polVsPhi;             // cuts
  std::array<PerEventProfile2D, kNPolComponents> polVsPhiVsPt;         // y cut only
  std::array<PerEventProfile1D, kNPolComponents> polJetVsDeltaPhi;     // cuts
  std::array<PerEventProfile2D, kNPolComponents> polJetVsDeltaPhiVsPt; // y cut only
  // PolarizationMaps/ (full grid, smoothed): P* components, yields, and on the PrimeJet planes the ring scalars
  std::array<PerEventProfile2D, kNMapComponents> labMapPxPy, labMapPzPx, labMapPyPz;
  BufferedHist2D labYieldPxPy, labYieldPzPx, labYieldPyPz;
  std::array<PerEventProfile2D, kNMapComponents> primeJetMapPxPy, primeJetMapPzPx, primeJetMapPyPz;
  BufferedHist2D primeJetYieldPxPy, primeJetYieldPzPx, primeJetYieldPyPz;
  std::array<PerEventProfile2D, kNMapRings> primeJetRingPxPy, primeJetRingPzPx, primeJetRingPyPz;
  std::array<PerEventProfile2D, kNMapComponents> primeV0MapPzPt;
  BufferedHist2D primeV0YieldPzPt;
  // QA/
  PerEventProfile1D polEVsDeltaPhi, polParVsDeltaPhi, alignmentVsDeltaPhi; // cuts. P_R is Ring/pRingVsDeltaPhi itself
  PerEventProfile1D ringVsChi, ringZVsChi;
  BufferedHist1D chiYield;
  BufferedHist1D closurePdotP, closurePStarTime, closureBoost, closureRingLabMinusStar;
  // Per-event values, for the distributions built at the end
  std::vector<double> ringD3NPerEvent, ringFlatPerEvent, ringZD3NPerEvent, ringZFlatPerEvent, eventYieldPerEvent;

  // Naming scheme may be confusing, but notice that flushAndReset is a method of TrueTriggerOutputs,
  // and it calls same-named methods from, for instance, PerEventProfile1D:
  // (no recursive functions nor anything as complex as that!)
  void flushAndReset()
  {
    for (PerEventProfile1D& profile : obsVsDeltaPhi)
      profile.flushAndReset();
    ringVsDeltaPhiVsPt.flushAndReset();
    ringZVsDeltaPhiVsPt.flushAndReset();
    ringVsDeltaPhiVsY.flushAndReset();
    ringZVsDeltaPhiVsY.flushAndReset();
    ringVsPt.flushAndReset();
    ringZVsPt.flushAndReset();
    ringVsY.flushAndReset();
    ringZVsY.flushAndReset();
    ringVsDeltaPhiNoCuts.flushAndReset();
    ringZVsDeltaPhiNoCuts.flushAndReset();
    for (int iComp = 0; iComp < kNPolComponents; ++iComp) {
      polVsPhi[iComp].flushAndReset();
      polVsPhiVsPt[iComp].flushAndReset();
      polJetVsDeltaPhi[iComp].flushAndReset();
      polJetVsDeltaPhiVsPt[iComp].flushAndReset();
    }
    for (int iComp = 0; iComp < kNMapComponents; ++iComp) {
      labMapPxPy[iComp].flushAndReset();
      labMapPzPx[iComp].flushAndReset();
      labMapPyPz[iComp].flushAndReset();
      primeJetMapPxPy[iComp].flushAndReset();
      primeJetMapPzPx[iComp].flushAndReset();
      primeJetMapPyPz[iComp].flushAndReset();
      primeV0MapPzPt[iComp].flushAndReset();
    }
    for (int iRing = 0; iRing < kNMapRings; ++iRing) {
      primeJetRingPxPy[iRing].flushAndReset();
      primeJetRingPzPx[iRing].flushAndReset();
      primeJetRingPyPz[iRing].flushAndReset();
    }
    polEVsDeltaPhi.flushAndReset();
    polParVsDeltaPhi.flushAndReset();
    alignmentVsDeltaPhi.flushAndReset();
    ringVsChi.flushAndReset();
    ringZVsChi.flushAndReset();
  }
};

/// \brief The reduced set of one trigger variant (Mixing/ or Smudge/).
struct ReducedSetOutputs {
  TProfile* ringVsDeltaPhi = nullptr;
  TProfile* ringVsDeltaPhiEvtMean = nullptr;
  TProfile* ringZVsDeltaPhi = nullptr;
  TProfile* ringZVsDeltaPhiEvtMean = nullptr;
  TProfile* nzSqVsDeltaPhi = nullptr;
  TProfile* ringD3N = nullptr;
  TProfile* ringD3NEvtMean = nullptr;
  TProfile* ringFlat = nullptr;
  TProfile* ringFlatEvtMean = nullptr;
  TProfile* ringZD3N = nullptr;
  TProfile* ringZD3NEvtMean = nullptr;
  TProfile* ringZFlat = nullptr;
  TProfile* ringZFlatEvtMean = nullptr;
  // Overlap with the true ring direction (README, "Ring overlap")
  TProfile* cosZetaVsDeltaPhi = nullptr;
  TProfile* sinZetaVsDeltaPhi = nullptr;
  TProfile* cosZetaD3N = nullptr;
  TProfile* sinZetaD3N = nullptr;
  BufferedHist1D cosZetaDistribution, sinZetaDistribution;
  std::vector<double> ringD3NPerEvent, ringFlatPerEvent, ringZD3NPerEvent, ringZFlatPerEvent;
};

struct DeltaScanOutputs {
  TProfile* ringD3NVsDelta = nullptr;
  TProfile* ringD3NVsDeltaEvtMean = nullptr;
  TProfile* ringFlatVsDelta = nullptr;
  TProfile* ringFlatVsDeltaEvtMean = nullptr;
  TProfile* ringZD3NVsDelta = nullptr;
  TProfile* ringZD3NVsDeltaEvtMean = nullptr;
  TProfile* ringZFlatVsDelta = nullptr;
  TProfile* ringZFlatVsDeltaEvtMean = nullptr;
  TProfile* nzSqD3NVsDelta = nullptr;
  TProfile* nzSqD3NVsDeltaEvtMean = nullptr;
  TProfile2D* ringVsDeltaPhiVsDelta = nullptr;
  TProfile2D* ringZVsDeltaPhiVsDelta = nullptr;
  TH1D* closureHalfTurn = nullptr;
  // Overlap with the true ring direction (README, "Ring overlap")
  TProfile* cosZetaD3NVsDelta = nullptr;
  TProfile* sinZetaD3NVsDelta = nullptr;
  BufferedHist2D cosZetaVsDelta, sinZetaVsDelta;
};

struct MixingExtras {
  TProfile* ringD3NVsDeltaEff = nullptr;
  TProfile* ringD3NVsDeltaEffEvtMean = nullptr;
  TProfile* ringZD3NVsDeltaEff = nullptr;
  TProfile* ringZD3NVsDeltaEffEvtMean = nullptr;
};

struct SmudgeExtras {
  TH1D* smudgeAngle = nullptr;
  TH1D* triggerEta = nullptr;
  TProfile* ringD3NVsSmudgeAngle = nullptr;
  TProfile* ringD3NVsSmudgeAngleEvtMean = nullptr;
  TProfile* ringZD3NVsSmudgeAngle = nullptr;
  TProfile* ringZD3NVsSmudgeAngleEvtMean = nullptr;
};

/// \brief Branch buffers of the event tree (README, "Event tree").
struct EventRecord {
  Int_t eventId = 0;
  Double_t phiJ = kNoValue;
  Double_t phiTrigger = kNoValue;
  Int_t nCellsInCuts = 0;
  Double_t sumWInCuts = 0.;
  std::array<Double_t, kNRingObservables> observableD3N{};
  Double_t ringFlat = kNoValue;
  Double_t ringZFlat = kNoValue;
  Double_t polED3N = kNoValue; // P_e = P*.(p_hat x n) of the true trigger, the B of the leakage ansatz (README, "Ring overlap")
  Int_t nDeltaPhiBins = 0;
  std::vector<Double_t> sumWDeltaPhi, sumARingDeltaPhi, sumARingZDeltaPhi, sumANzSqDeltaPhi;
  Int_t nDelta = 0;
  std::vector<Double_t> scanRingD3N, scanRingFlat, scanRingZD3N, scanRingZFlat, scanNzSqD3N, scanCosZetaD3N, scanSinZetaD3N;
  Double_t mixRingD3N = kNoValue;
  Double_t mixRingFlat = kNoValue;
  Double_t mixRingZD3N = kNoValue;
  Double_t mixRingZFlat = kNoValue;
  Int_t nMixPartners = 0;     // mixDepth
  std::vector<Int_t> mixDeltaEffBin; // per partner: its delta_eff bin of Mixing/p*VsDeltaEff (from 0), -1 outside the axis
  Int_t nMixDeltaPhiBins = 0; // nDeltaPhiBins when mixing is on, 0 otherwise
  std::vector<Double_t> mixSumWDeltaPhi, mixSumARingDeltaPhi, mixSumARingZDeltaPhi, mixSumANzSqDeltaPhi, mixSumACosZetaDeltaPhi, mixSumASinZetaDeltaPhi;
  Double_t smudgeRingD3N = kNoValue;
  Double_t smudgeRingFlat = kNoValue;
  Double_t smudgeRingZD3N = kNoValue;
  Double_t smudgeRingZFlat = kNoValue;
  Int_t nSmudgeDeltaPhiBins = 0; // nDeltaPhiBins when smudging is on, 0 otherwise
  std::vector<Double_t> smudgeSumWDeltaPhi, smudgeSumARingDeltaPhi, smudgeSumARingZDeltaPhi, smudgeSumANzSqDeltaPhi, smudgeSumACosZetaDeltaPhi, smudgeSumASinZetaDeltaPhi;
};

/// \brief Counters of everything skipped or suspicious, reported at the end of the run and stored in Config/ (nothing is skipped silently).
struct WarningCounters {
  long long invalidCells = 0;            // non-finite or negative columns: the cell is skipped in that event
  long long emptyDeltaPhiBins = 0;       // the grid binning should make this impossible
  long long pointsOutsideMaps = 0;       // PolarizationMaps/ points outside the map ranges (the ranges are built to hold every cell)
  long long degenerateRingDirection = 0; // p parallel to t, so n is undefined
  long long zeroPolStarCells = 0;        // |P*| = 0, skipped in the alignment P_R/|P*| only
  std::vector<std::string> startupMessages;
};

// ======================================================================================================================================
// Pieces used by several sections of main()
// ======================================================================================================================================

/// \brief Integer shift that maps the phi index of a cell to its Delta phi bin for a trigger at phiTrigger: bin = (iPhi + shift) mod nPhi.
/// \note  The Delta phi axis is [-pi - w/2, pi - w/2) with w = 2 pi/nPhi, so its bin centres sit on the phi grid. A cell sits at phiFirst + iPhi * w,
///        so floor((Delta phi + pi + w/2)/w) is iPhi plus an event-wide integer, and a grid of spacing w puts exactly one node in every bin,
///        whatever the trigger angle.
int deltaPhiBinShift(double phiTrigger, const CellGrid& grid)
{
  const long shift = std::floor((grid.phiFirst - phiTrigger + kPi + 0.5 * grid.phiWidth) / grid.phiWidth);
  return ((shift % grid.nPhi) + grid.nPhi) % grid.nPhi;
}

/// \brief Accumulates the reduced set of one trigger axis over the in-cut cells, binned in Delta phi from the trigger azimuth, together with
///        the overlap of its ring direction n with the true one n_t (inCutCells.nTrue*, filled for the event beforehand).
/// \param triggerUnitVec trigger direction, unit norm. t_z may be nonzero (cone smudging), so the full cross product is used.
/// \return the integrated sums of this trigger alone (the accumulators keep summing across the triggers of one event)
IntegratedSums accumulateReducedSet(const XYZVector& triggerUnitVec, const InCutCells& inCutCells, const CellGrid& grid, ReducedAccumulators& accumulators, WarningCounters& warnings)
{
  const double triggerX = triggerUnitVec.X();
  const double triggerY = triggerUnitVec.Y();
  const double triggerZ = triggerUnitVec.Z();
  const int binShift = deltaPhiBinShift(std::atan2(triggerY, triggerX), grid);
  const int nPhiBins = grid.nPhi;
  const int nInCutCells = inCutCells.size();
  IntegratedSums sums;
  for (int iCut = 0; iCut < nInCutCells; ++iCut) {
    const double crossX = triggerY * inCutCells.pz[iCut] - triggerZ * inCutCells.py[iCut];
    const double crossY = triggerZ * inCutCells.px[iCut] - triggerX * inCutCells.pz[iCut];
    const double crossZ = triggerX * inCutCells.py[iCut] - triggerY * inCutCells.px[iCut];
    const double crossNormSq = crossX * crossX + crossY * crossY + crossZ * crossZ;
    if (!(crossNormSq > 0.)) { // p parallel to t (only conceivable for a smudged trigger): n undefined
      ++warnings.degenerateRingDirection;
      continue;
    }
    const double invCrossNorm = 1. / std::sqrt(crossNormSq);
    const double nX = crossX * invCrossNorm;
    const double nY = crossY * invCrossNorm;
    const double nZ = crossZ * invCrossNorm;
    const double ring = inCutCells.polStarX[iCut] * nX + inCutCells.polStarY[iCut] * nY + inCutCells.polStarZ[iCut] * nZ;
    const double ringZ = inCutCells.polStarZ[iCut] * nZ;
    // Overlap with the true ring direction, as computeRingOverlap in the consumer: zeta turns n_t into n, right-handed about p_hat. Both are
    // already unit vectors here, so the direct products replace the consumer's Binet-Cauchy form
    const double nTrueX = inCutCells.nTrueX[iCut];
    const double nTrueY = inCutCells.nTrueY[iCut];
    const double nTrueZ = inCutCells.nTrueZ[iCut];
    const double cosZeta = nTrueX * nX + nTrueY * nY + nTrueZ * nZ;
    const double sinZeta = (nTrueY * nZ - nTrueZ * nY) * inCutCells.pHatX[iCut] + (nTrueZ * nX - nTrueX * nZ) * inCutCells.pHatY[iCut] + (nTrueX * nY - nTrueY * nX) * inCutCells.pHatZ[iCut];
    int bin = inCutCells.phiIndex[iCut] + binShift;
    if (bin >= nPhiBins)
      bin -= nPhiBins;
    const double weight = inCutCells.weight[iCut];
    accumulators.sumW[bin] += weight;
    accumulators.sumWRing[bin] += weight * ring;
    accumulators.sumWRingZ[bin] += weight * ringZ;
    accumulators.sumWNzSq[bin] += weight * nZ * nZ;
    accumulators.sumWCosZeta[bin] += weight * cosZeta;
    accumulators.sumWSinZeta[bin] += weight * sinZeta;
    accumulators.cosZetaDistribution[zetaBin(cosZeta)] += weight;
    accumulators.sinZetaDistribution[zetaBin(sinZeta)] += weight;
    sums.sumW += weight;
    sums.sumWRing += weight * ring;
    sums.sumWRingZ += weight * ringZ;
    sums.sumWNzSq += weight * nZ * nZ;
    sums.sumWCosZeta += weight * cosZeta;
    sums.sumWSinZeta += weight * sinZeta;
  }
  return sums;
}

/// \brief Per-event integrals of a reduced set, its flush into the masters, and the reset of its accumulators.
ReducedValues flushReducedSet(ReducedAccumulators& accumulators, ReducedSetOutputs& outputs, WarningCounters& warnings)
{
  ReducedValues values;
  const double sumW = sumOf(accumulators.sumW);
  values.ringD3N = sumOf(accumulators.sumWRing) / sumW;
  values.ringZD3N = sumOf(accumulators.sumWRingZ) / sumW;
  warnings.emptyDeltaPhiBins += countEmptyBins(accumulators.sumW);
  values.ringFlat = flatAverage(accumulators.sumW, accumulators.sumWRing);
  values.ringZFlat = flatAverage(accumulators.sumW, accumulators.sumWRingZ);
  values.cosZetaD3N = sumOf(accumulators.sumWCosZeta) / sumW;
  values.sinZetaD3N = sumOf(accumulators.sumWSinZeta) / sumW;

  flushPerEventBins(outputs.ringVsDeltaPhi, outputs.ringVsDeltaPhiEvtMean, accumulators.sumW, accumulators.sumWRing);
  flushPerEventBins(outputs.ringZVsDeltaPhi, outputs.ringZVsDeltaPhiEvtMean, accumulators.sumW, accumulators.sumWRingZ);
  flushPerEventBins(outputs.nzSqVsDeltaPhi, nullptr, accumulators.sumW, accumulators.sumWNzSq);
  flushPerEventBins(outputs.cosZetaVsDeltaPhi, nullptr, accumulators.sumW, accumulators.sumWCosZeta);
  flushPerEventBins(outputs.sinZetaVsDeltaPhi, nullptr, accumulators.sumW, accumulators.sumWSinZeta);
  fillEventPoint(outputs.ringD3N, outputs.ringD3NEvtMean, 0.5, values.ringD3N, sumW);
  fillEventPoint(outputs.ringFlat, outputs.ringFlatEvtMean, 0.5, values.ringFlat, sumW);
  fillEventPoint(outputs.ringZD3N, outputs.ringZD3NEvtMean, 0.5, values.ringZD3N, sumW);
  fillEventPoint(outputs.ringZFlat, outputs.ringZFlatEvtMean, 0.5, values.ringZFlat, sumW);
  fillEventPoint(outputs.cosZetaD3N, nullptr, 0.5, values.cosZetaD3N, sumW);
  fillEventPoint(outputs.sinZetaD3N, nullptr, 0.5, values.sinZetaD3N, sumW);
  for (int bin = 0; bin < kZetaBins; ++bin) {
    outputs.cosZetaDistribution.add(bin + 1, accumulators.cosZetaDistribution[bin]);
    outputs.sinZetaDistribution.add(bin + 1, accumulators.sinZetaDistribution[bin]);
  }
  outputs.ringD3NPerEvent.push_back(values.ringD3N);
  outputs.ringFlatPerEvent.push_back(values.ringFlat);
  outputs.ringZD3NPerEvent.push_back(values.ringZD3N);
  outputs.ringZFlatPerEvent.push_back(values.ringZFlat);
  accumulators.reset();
  return values;
}

/// \brief Books the reduced set of one trigger variant (used by Mixing/ and Smudge/).
void bookReducedSet(TDirectory* directory, ReducedSetOutputs& outputs, const AxisSpec& deltaPhiAxis, const AxisSpec& integratedAxis, const AxisSpec& zetaAxis)
{
  directory->cd();
  outputs.ringVsDeltaPhi = new TProfile("pRingVsDeltaPhi", ("<#it{R}> vs #Delta#varphi (pooled);" + kDeltaPhiLabel + ";<#it{R}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high);
  outputs.ringVsDeltaPhiEvtMean = new TProfile("pRingVsDeltaPhi_EvtMean", ("<#it{R}> vs #Delta#varphi (event mean);" + kDeltaPhiLabel + ";<#it{R}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high);
  outputs.ringZVsDeltaPhi = new TProfile("pRingZVsDeltaPhi", ("<#it{R}_{z}> vs #Delta#varphi (pooled);" + kDeltaPhiLabel + ";<#it{R}_{z}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high);
  outputs.ringZVsDeltaPhiEvtMean = new TProfile("pRingZVsDeltaPhi_EvtMean", ("<#it{R}_{z}> vs #Delta#varphi (event mean);" + kDeltaPhiLabel + ";<#it{R}_{z}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high);
  outputs.nzSqVsDeltaPhi = new TProfile("pNzSqVsDeltaPhi", ("<n_{z}^{2}> vs #Delta#varphi (pooled);" + kDeltaPhiLabel + ";<n_{z}^{2}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high);
  outputs.ringD3N = new TProfile("pRingD3N", "Integrated <#it{R}> (d^{3}N-weighted, pooled);;<#it{R}>", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
  outputs.ringD3NEvtMean = new TProfile("pRingD3N_EvtMean", "Integrated <#it{R}> (d^{3}N-weighted, event mean);;<#it{R}>", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
  outputs.ringFlat = new TProfile("pRingFlat", "Integrated #it{R}, flat in #Delta#varphi (events weighted by their in-cut yield);;#it{R}_{flat}", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
  outputs.ringFlatEvtMean = new TProfile("pRingFlat_EvtMean", "Integrated #it{R}, flat in #Delta#varphi (event mean);;#it{R}_{flat}", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
  outputs.ringZD3N = new TProfile("pRingZD3N", "Integrated <#it{R}_{z}> (d^{3}N-weighted, pooled);;<#it{R}_{z}>", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
  outputs.ringZD3NEvtMean = new TProfile("pRingZD3N_EvtMean", "Integrated <#it{R}_{z}> (d^{3}N-weighted, event mean);;<#it{R}_{z}>", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
  outputs.ringZFlat = new TProfile("pRingZFlat", "Integrated #it{R}_{z}, flat in #Delta#varphi (events weighted by their in-cut yield);;#it{R}_{z,flat}", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
  outputs.ringZFlatEvtMean = new TProfile("pRingZFlat_EvtMean", "Integrated #it{R}_{z}, flat in #Delta#varphi (event mean);;#it{R}_{z,flat}", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
  outputs.cosZetaVsDeltaPhi = new TProfile("pCosZetaVsDeltaPhi", ("<cos#zeta> = <#hat{n}_{t}#upoint#hat{n}> vs #Delta#varphi (pooled);" + kDeltaPhiLabel + ";<cos#zeta>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high);
  outputs.sinZetaVsDeltaPhi = new TProfile("pSinZetaVsDeltaPhi", ("<sin#zeta> = <(#hat{n}_{t} #times #hat{n})#upoint#hat{p}> vs #Delta#varphi (pooled);" + kDeltaPhiLabel + ";<sin#zeta>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high);
  outputs.cosZetaD3N = new TProfile("pCosZetaD3N", "Integrated <cos#zeta> (d^{3}N-weighted, pooled);;<cos#zeta>", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
  outputs.sinZetaD3N = new TProfile("pSinZetaD3N", "Integrated <sin#zeta> (d^{3}N-weighted, pooled);;<sin#zeta>", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
  outputs.cosZetaDistribution.attach(new TH1D("hCosZeta", "cos#zeta = #hat{n}_{t}#upoint#hat{n} (d^{3}N-weighted, all events);cos#zeta;#sum w", zetaAxis.nBins, zetaAxis.low, zetaAxis.high));
  outputs.sinZetaDistribution.attach(new TH1D("hSinZeta", "sin#zeta = (#hat{n}_{t} #times #hat{n})#upoint#hat{p} (d^{3}N-weighted, all events);sin#zeta;#sum w", zetaAxis.nBins, zetaAxis.low, zetaAxis.high));
}

// ======================================================================================================================================
// main
// ======================================================================================================================================

int main(int argc, char** argv)
{
  try {
    // ==================================================================================================================================
    // 1. Options: parse, validate, and keep the (name, value) list for the run header and for Config/
    // ==================================================================================================================================
    AnalysisConfig config;
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
    if (positionalArgs.size() != 2) {
      printUsage();
      throw std::runtime_error("expected exactly two positional arguments, <sampleDir> and <outputFile>");
    }
    config.sampleDir = positionalArgs[0];
    config.outputFile = positionalArgs[1];

    bool seedGiven = false;
    for (const auto& [key, value] : options) {
      if (key == "jetFile") {
        config.jetFile = value;
      } else if (key == "trigger") {
        if (value == "jet")
          config.triggerMode = TriggerMode::kJet;
        else if (value == "rp")
          config.triggerMode = TriggerMode::kReactionPlane;
        else
          throw std::runtime_error("--trigger must be jet or rp, got '" + value + "'");
      } else if (key == "firstEvent") {
        config.firstEvent = toInt(key, value);
      } else if (key == "nEvents") {
        config.nEvents = toInt(key, value);
      } else if (key == "inputName") {
        config.inputName = value;
      } else if (key == "ptMin") {
        config.ptMin = toDouble(key, value);
      } else if (key == "ptMax") {
        config.ptMax = toDouble(key, value);
      } else if (key == "absYMax") {
        config.absYMax = toDouble(key, value);
      } else if (key == "lambdaMass") {
        config.lambdaMass = toDouble(key, value);
      } else if (key == "deltaScanSteps") {
        config.deltaScanSteps = toInt(key, value);
      } else if (key == "mixDepth") {
        config.mixDepth = toInt(key, value);
      } else if (key == "smudgeSigma") {
        config.smudgeSigma = toDouble(key, value);
      } else if (key == "smudgeMode") {
        if (value == "cone")
          config.smudgeMode = SmudgeMode::kCone;
        else if (value == "phi")
          config.smudgeMode = SmudgeMode::kPhi;
        else
          throw std::runtime_error("--smudgeMode must be cone or phi, got '" + value + "'");
      } else if (key == "smudgeDraws") {
        config.smudgeDraws = toInt(key, value);
      } else if (key == "seed") {
        config.seed = toUInt32(key, value);
        seedGiven = true;
      } else if (key == "mapSubdivisions") {
        config.mapSubdivisions = toInt(key, value);
      } else {
        throw std::runtime_error("unknown option --" + key + " (see --help)");
      }
    }

    // Checks that need no input file (the grid-dependent ones run once the first event is read)
    if (config.firstEvent < 0 || config.nEvents < 1)
      throw std::runtime_error("need --firstEvent >= 0 and --nEvents >= 1");
    if (!(config.ptMin < config.ptMax) || !(config.absYMax > 0.))
      throw std::runtime_error("need --ptMin < --ptMax and --absYMax > 0");
    if (config.triggerMode == TriggerMode::kJet && config.jetFile.empty())
      throw std::runtime_error("--trigger=jet needs --jetFile");
    if (config.deltaScanSteps < 0 || config.mixDepth < 0 || config.smudgeSigma < 0. || config.smudgeDraws < 1 || config.mapSubdivisions < 1)
      throw std::runtime_error("need --deltaScanSteps >= 0, --mixDepth >= 0, --smudgeSigma >= 0, --smudgeDraws >= 1 and --mapSubdivisions >= 1");
    // With phi_t = 0 for every event, all mixed triggers coincide with the true one
    if (config.mixDepth > 0 && config.triggerMode != TriggerMode::kJet)
      throw std::runtime_error("mixing needs --trigger=jet");
    if (config.mixDepth >= config.nEvents)
      throw std::runtime_error("--mixDepth must be smaller than --nEvents (an event cannot be its own partner)");
    // Without --seed, the low 32 bits of the current time in nanoseconds. Two producers started separately then get different seeds, which
    // breaks the pairing of the smudged triggers: the coordinator always passes one seed to both samples
    if (!seedGiven)
      config.seed = std::chrono::system_clock::now().time_since_epoch().count();

    // Every effective option as (name, value) strings, for the run header and for Config/
    std::vector<std::pair<std::string, std::string>> configEntries;
    configEntries.emplace_back("sampleDir", config.sampleDir);
    configEntries.emplace_back("outputFile", config.outputFile);
    configEntries.emplace_back("jetFile", config.jetFile.empty() ? "none" : config.jetFile);
    configEntries.emplace_back("trigger", config.triggerMode == TriggerMode::kJet ? "jet" : "rp");
    configEntries.emplace_back("firstEvent", std::to_string(config.firstEvent));
    configEntries.emplace_back("nEvents", std::to_string(config.nEvents));
    configEntries.emplace_back("inputName", config.inputName);
    configEntries.emplace_back("ptMin", std::to_string(config.ptMin));
    configEntries.emplace_back("ptMax", std::to_string(config.ptMax));
    configEntries.emplace_back("absYMax", std::to_string(config.absYMax));
    configEntries.emplace_back("lambdaMass", std::to_string(config.lambdaMass));
    configEntries.emplace_back("deltaScanSteps", std::to_string(config.deltaScanSteps));
    configEntries.emplace_back("mixDepth", std::to_string(config.mixDepth));
    configEntries.emplace_back("smudgeSigma", std::to_string(config.smudgeSigma));
    configEntries.emplace_back("smudgeMode", config.smudgeMode == SmudgeMode::kCone ? "cone" : "phi");
    configEntries.emplace_back("smudgeDraws", std::to_string(config.smudgeDraws));
    configEntries.emplace_back("seed", std::to_string(config.seed));
    configEntries.emplace_back("mapSubdivisions", std::to_string(config.mapSubdivisions));

    std::cout << "========================================================" << std::endl;
    std::cout << "  hydroTruthLab" << std::endl;
    for (const auto& [name, value] : configEntries)
      std::cout << "  " << name << " = " << value << std::endl;
    std::cout << "  column-to-polarization factor = " << kColumnToPolarization << " (the iSS S^mu columns already hold P^mu, see README)" << std::endl;
    std::cout << "========================================================" << std::endl;

    // ==================================================================================================================================
    // 2. Jet file, with atan2(momentum_y, momentum_x) checked against phi. Rows have five fields (event n_random phi(rad) momentum_x
    //    momentum_y), or four without n_random: the files of lumpy_events_0_5_fine and break_symmetry keep the five-column header over rows of
    //    four, so the layout is taken from each row, not from the header
    // ==================================================================================================================================
    std::vector<JetEntry> jets;
    if (!config.jetFile.empty()) {
      std::ifstream jetStream(config.jetFile);
      if (!jetStream)
        throw std::runtime_error("cannot open the jet file " + config.jetFile);
      std::string line;
      int lineNumber = 0;
      while (std::getline(jetStream, line)) {
        ++lineNumber;
        const std::size_t firstChar = line.find_first_not_of(" \t\r");
        if (firstChar == std::string::npos || line[firstChar] == '#')
          continue;
        std::istringstream lineStream(line);
        std::vector<double> fields;
        double field = 0.;
        while (lineStream >> field)
          fields.push_back(field);
        if (!lineStream.eof() || (fields.size() != 5 && fields.size() != 4))
          throw std::runtime_error(config.jetFile + ", line " + std::to_string(lineNumber) + ": expected 5 numeric fields (event n_random phi momentum_x momentum_y) or 4 (without n_random)");
        const std::size_t phiField = fields.size() - 3;
        const double eventColumn = fields[0];
        const double phiColumn = fields[phiField];
        const double momentumX = fields[phiField + 1];
        const double momentumY = fields[phiField + 2];
        const long jetEventId = std::lround(eventColumn); // read as a number, since the event column may be printed in floating-point format
        if (std::abs(eventColumn - jetEventId) > 1e-9 || jetEventId < 0)
          throw std::runtime_error(config.jetFile + ", line " + std::to_string(lineNumber) + ": the event column is not a non-negative integer");
        const double angleMismatch = std::abs(std::remainder(std::atan2(momentumY, momentumX) - phiColumn, kTwoPi));
        if (angleMismatch > kJetAngleTolerance)
          throw std::runtime_error(config.jetFile + ", line " + std::to_string(lineNumber) + ": atan2(momentum_y, momentum_x) differs from phi by " + std::to_string(angleMismatch) + " rad");
        if (jetEventId >= static_cast<long>(jets.size()))
          jets.resize(jetEventId + 1);
        if (jets[jetEventId].present)
          throw std::runtime_error(config.jetFile + ": event " + std::to_string(jetEventId) + " appears twice");
        jets[jetEventId] = {true, phiColumn};
      }
    }
    // The jet trigger and mixing need an entry for every event of the run: stop now rather than halfway
    if (config.triggerMode == TriggerMode::kJet || config.mixDepth > 0) {
      for (int eventId = config.firstEvent; eventId < config.firstEvent + config.nEvents; ++eventId)
        if (eventId >= static_cast<int>(jets.size()) || !jets[eventId].present)
          throw std::runtime_error("event " + std::to_string(eventId) + " has no entry in " + config.jetFile);
    }

    // ==================================================================================================================================
    // 3. Output file
    // ==================================================================================================================================
    std::unique_ptr<TFile> outputFile(TFile::Open(config.outputFile.c_str(), "RECREATE"));
    if (!outputFile || outputFile->IsZombie())
      throw std::runtime_error("cannot create " + config.outputFile);
    // Every histogram and profile created from here on gets Sumw2(). The pooled masters are filled with weights, and without it their effective
    // entries, hence their errors, are wrong
    TH1::SetDefaultSumw2(kTRUE);

    // ==================================================================================================================================
    // 4. First event: the grid it defines, which every later event must reproduce (README, "File format and grid"), and the axes built on it
    // ==================================================================================================================================
    // Events of the range that have their input file. A missing one (an iSS run that never finished, in some datasets) is skipped with a
    // warning instead of stopping the run: its mixing partners still take their phi_J from the jet file, and PerConfig pairs the samples by
    // eventId, so the pair is simply absent
    std::vector<int> eventIds, skippedEventIds;
    for (int eventId = config.firstEvent; eventId < config.firstEvent + config.nEvents; ++eventId) {
      if (std::filesystem::is_regular_file(eventInputPath(config, eventId)))
        eventIds.push_back(eventId);
      else
        skippedEventIds.push_back(eventId);
    }
    if (eventIds.empty())
      throw std::runtime_error("no event of " + std::to_string(config.firstEvent) + "-" + std::to_string(config.firstEvent + config.nEvents - 1) + " has " + config.inputName);
    const std::string skippedEventsText = skippedEventIds.empty() ? "none" : formatEventRanges(skippedEventIds);
    CellTable table;
    readCellTable(eventInputPath(config, eventIds.front()), table);

    CellGrid grid;
    grid.nCells = table.size();
    const std::vector<double> yNodes = determineGridShape(table.column[kColY]);
    const std::vector<double> ptNodes = determineGridShape(table.column[kColPt]);
    const std::vector<double> phiNodes = determineGridShape(table.column[kColPhi]);
    grid.nY = yNodes.size();
    grid.nPt = ptNodes.size();
    grid.nPhi = phiNodes.size();
    if (grid.nY < 2 || grid.nPt < 2 || grid.nPhi < 2)
      throw std::runtime_error("the grid needs at least two nodes in y, pT and phi");
    if (grid.nY * grid.nPt * grid.nPhi != grid.nCells)
      throw std::runtime_error("the grid is not a full (y, pT, phi) product: " + std::to_string(grid.nY) + " x " + std::to_string(grid.nPt) + " x " + std::to_string(grid.nPhi) + " != " + std::to_string(grid.nCells) + " cells");

    // Widths from the nodes themselves, never hardcoded
    grid.yFirst = yNodes.front();
    grid.yWidth = (yNodes.back() - yNodes.front()) / (grid.nY - 1);
    grid.ptFirst = ptNodes.front();
    grid.ptWidth = (ptNodes.back() - ptNodes.front()) / (grid.nPt - 1);
    grid.phiFirst = phiNodes.front();
    grid.phiWidth = kTwoPi / grid.nPhi;
    // phi is periodic: the nodes must cover the full circle, otherwise the one-node-per-bin property of the Delta phi binning breaks
    const double phiWidthFromNodes = (phiNodes.back() - phiNodes.front()) / (grid.nPhi - 1);
    if (std::abs(phiWidthFromNodes - grid.phiWidth) * grid.nPhi > 1e-4)
      throw std::runtime_error("the phi nodes do not cover the full circle uniformly");

    const int nCells = grid.nCells;
    grid.y.resize(nCells);
    grid.pt.resize(nCells);
    grid.phi.resize(nCells);
    grid.iY.resize(nCells);
    grid.iPt.resize(nCells);
    grid.iPhi.resize(nCells);
    grid.px.resize(nCells);
    grid.py.resize(nCells);
    grid.pz.resize(nCells);
    grid.energy.resize(nCells);
    grid.pMag.resize(nCells);
    grid.pHatX.resize(nCells);
    grid.pHatY.resize(nCells);
    grid.pHatZ.resize(nCells);
    grid.sinSqTheta.resize(nCells);
    grid.cellMeasure.resize(nCells);
    grid.boostToRest.resize(nCells);
    grid.inPtCut.resize(nCells);
    grid.inYCut.resize(nCells);

    InCutCells inCutCells;
    std::vector<char> cellSeen(nCells, 0);
    const double lambdaMassSq = squared(config.lambdaMass);
    for (int iCell = 0; iCell < nCells; ++iCell) {
      const double cellY = table.column[kColY][iCell];
      const double cellPt = table.column[kColPt][iCell];
      const double cellPhi = table.column[kColPhi][iCell];
      grid.y[iCell] = cellY;
      grid.pt[iCell] = cellPt;
      grid.phi[iCell] = cellPhi;
      grid.iY[iCell] = gridIndex(cellY, grid.yFirst, grid.yWidth, "y");
      grid.iPt[iCell] = gridIndex(cellPt, grid.ptFirst, grid.ptWidth, "pT");
      grid.iPhi[iCell] = gridIndex(cellPhi, grid.phiFirst, grid.phiWidth, "phi");
      const int flatIndex = (grid.iY[iCell] * grid.nPt + grid.iPt[iCell]) * grid.nPhi + grid.iPhi[iCell];
      if (cellSeen[flatIndex]++)
        throw std::runtime_error("grid node (y, pT, phi) = (" + std::to_string(cellY) + ", " + std::to_string(cellPt) + ", " + std::to_string(cellPhi) + ") appears twice");

      // Lambda kinematics from (y, pT, phi): p_z = m_T sinh(y), E = m_T cosh(y)
      const double transverseMass = std::sqrt(lambdaMassSq + squared(cellPt));
      grid.px[iCell] = cellPt * std::cos(cellPhi);
      grid.py[iCell] = cellPt * std::sin(cellPhi);
      grid.pz[iCell] = transverseMass * std::sinh(cellY);
      grid.energy[iCell] = transverseMass * std::cosh(cellY);
      grid.pMag[iCell] = std::sqrt(squared(cellPt) + squared(grid.pz[iCell]));
      grid.pHatX[iCell] = grid.px[iCell] / grid.pMag[iCell];
      grid.pHatY[iCell] = grid.py[iCell] / grid.pMag[iCell];
      grid.pHatZ[iCell] = grid.pz[iCell] / grid.pMag[iCell];
      grid.sinSqTheta[iCell] = squared(cellPt / grid.pMag[iCell]);
      // Boost trivector from the lab to the Lambda rest frame, as in the consumer (only used by the boost closure)
      grid.boostToRest[iCell] = PxPyPzEVector(grid.px[iCell], grid.py[iCell], grid.pz[iCell], grid.energy[iCell]).BoostToCM();
      // Nodes are cell centres (midpoint rule), so this is the momentum-space volume of the cell, times the pT Jacobian
      grid.cellMeasure[iCell] = cellPt * grid.ptWidth * grid.yWidth * grid.phiWidth;
      // The paper's cuts (see AnalysisConfig), on cell centres
      grid.inPtCut[iCell] = cellPt > config.ptMin && cellPt < config.ptMax;
      grid.inYCut[iCell] = std::abs(cellY) < config.absYMax;
      if (grid.inPtCut[iCell] && grid.inYCut[iCell]) {
        inCutCells.cellIndex.push_back(iCell);
        inCutCells.phiIndex.push_back(grid.iPhi[iCell]);
        inCutCells.px.push_back(grid.px[iCell]);
        inCutCells.py.push_back(grid.py[iCell]);
        inCutCells.pz.push_back(grid.pz[iCell]);
        inCutCells.pHatX.push_back(grid.pHatX[iCell]);
        inCutCells.pHatY.push_back(grid.pHatY[iCell]);
        inCutCells.pHatZ.push_back(grid.pHatZ[iCell]);
      }
    }
    const int nInCutCells = inCutCells.size();
    if (nInCutCells == 0)
      throw std::runtime_error("no cell passes the pT and y cuts");
    inCutCells.weight.resize(nInCutCells);
    inCutCells.polStarX.resize(nInCutCells);
    inCutCells.polStarY.resize(nInCutCells);
    inCutCells.polStarZ.resize(nInCutCells);
    inCutCells.nTrueX.resize(nInCutCells);
    inCutCells.nTrueY.resize(nInCutCells);
    inCutCells.nTrueZ.resize(nInCutCells);
    std::cout << "Grid: " << grid.nY << " y nodes from " << grid.yFirst << " (step " << grid.yWidth << "), " << grid.nPt << " pT nodes from " << grid.ptFirst << " (step " << grid.ptWidth << "), " << grid.nPhi << " phi nodes from " << grid.phiFirst << " (step 2pi/" << grid.nPhi << "). " << nInCutCells << " of " << nCells << " cells pass the cuts." << std::endl;

    // Options that depend on the grid
    WarningCounters warnings;
    if (!skippedEventIds.empty()) {
      const std::string message = std::to_string(skippedEventIds.size()) + " events of the range have no " + config.inputName + " and were skipped: " + skippedEventsText;
      std::cout << "[WARNING] " << message << std::endl;
      warnings.startupMessages.push_back(message);
    }
    // The delta grid must sit on the phi grid (rotated cells then land on bin centres), and be even (half-turn closure)
    if (config.deltaScanSteps > 0 && (config.deltaScanSteps % 2 != 0 || grid.nPhi % config.deltaScanSteps != 0))
      throw std::runtime_error("--deltaScanSteps must be even and divide the number of phi nodes (" + std::to_string(grid.nPhi) + ")");
    // Cuts act on cell centres, so a cut inside a cell keeps or drops the whole cell: legal, but the effective cut is not the requested one
    const bool ptMinOnEdge = isOnCellEdge(config.ptMin, grid.ptFirst, grid.ptWidth);
    const bool ptMaxOnEdge = isOnCellEdge(config.ptMax, grid.ptFirst, grid.ptWidth);
    const bool yCutOnEdges = isOnCellEdge(config.absYMax, grid.yFirst, grid.yWidth) && isOnCellEdge(-config.absYMax, grid.yFirst, grid.yWidth);
    if (!ptMinOnEdge || !ptMaxOnEdge || !yCutOnEdges) {
      const std::string message = std::string("a kinematic cut falls inside a cell (ptMin on an edge: ") + (ptMinOnEdge ? "yes" : "no") + ", ptMax on an edge: " + (ptMaxOnEdge ? "yes" : "no") + ", |y| cut on edges: " + (yCutOnEdges ? "yes" : "no") + "): whole cells are kept or dropped by their centres";
      std::cout << "[WARNING] " << message << std::endl;
      warnings.startupMessages.push_back(message);
    }

    // Axes (README, "Binning"). Delta phi (and phi - Psi_RP): bin centres on the phi grid, [-pi - w/2, pi - w/2) with w = 2 pi / nPhi
    const AxisSpec deltaPhiAxis{grid.nPhi, -kPi - 0.5 * grid.phiWidth, kPi - 0.5 * grid.phiWidth};
    const AxisSpec ptAxis{grid.nPt, grid.ptFirst - 0.5 * grid.ptWidth, grid.ptFirst + (grid.nPt - 0.5) * grid.ptWidth};
    const AxisSpec yAxis{grid.nY, grid.yFirst - 0.5 * grid.yWidth, grid.yFirst + (grid.nY - 0.5) * grid.yWidth};
    const AxisSpec integratedAxis{1, 0., 1.};
    // Delta grid: bin centres on delta_k = 2 pi k / nSteps. Mixing reuses it for delta_eff (with the phi binning if the scan is off)
    const int nDeltaBins = config.deltaScanSteps > 0 ? config.deltaScanSteps : grid.nPhi;
    const AxisSpec deltaAxis{nDeltaBins, -kPi / nDeltaBins, kTwoPi - kPi / nDeltaBins};
    // Momentum axes of PolarizationMaps/, with the pT cell width as bin width and ranges that just hold the momentum volume of the grid (the
    // smoothed maps spread every cell over its own volume): transverse components up to the last pT edge, |p_z| up to m_T sinh|y| at the outer
    // pT and y edges, both rounded up to a whole bin
    const int nTransverseHalfBins = std::ceil(ptAxis.high / grid.ptWidth);
    const double maxAbsYEdge = std::max(std::abs(yAxis.low), std::abs(yAxis.high));
    const int nPzHalfBins = std::ceil(std::sqrt(lambdaMassSq + squared(ptAxis.high)) * std::sinh(maxAbsYEdge) / grid.ptWidth);
    const AxisSpec transverseMomentumAxis{2 * nTransverseHalfBins, -nTransverseHalfBins * grid.ptWidth, nTransverseHalfBins * grid.ptWidth};
    const AxisSpec longitudinalMomentumAxis{2 * nPzHalfBins, -nPzHalfBins * grid.ptWidth, nPzHalfBins * grid.ptWidth};
    const AxisSpec chiAxis{grid.nPhi, -kPi, kPi};
    const AxisSpec phiTriggerAxis{grid.nPhi, 0., kTwoPi};
    const AxisSpec logResidualAxis{kLogResidualBins, kLogResidualLow, kLogResidualHigh};
    const AxisSpec zetaAxis{kZetaBins, -1., 1.};
    const AxisSpec smudgeAngleAxis{50, 0., 5. * config.smudgeSigma};
    const AxisSpec triggerEtaAxis{50, -5. * config.smudgeSigma, 5. * config.smudgeSigma};
    const int reactionPlaneBinShift = deltaPhiBinShift(0., grid); // phi wrapped like Delta phi with phi_t = 0, i.e. phi - Psi_RP
    const int nPhiBins = grid.nPhi;

    // Sub-cell points of PolarizationMaps/ (README, "Polarization maps"): offsets (k + 1/2)/n - 1/2 of a cell width along y, pT and phi, so
    // n = 1 is the node itself. Within a cell, d3N/(pT dpT dy dphi) and the field are those of the node (the midpoint rule of every other output),
    // so a point carries the d3N share (pT_point/pT_node)/n^3 of its cell
    SubCellPoints subCellPoints;
    const int nSubdivisions = config.mapSubdivisions;
    subCellPoints.perCell = nSubdivisions * nSubdivisions * nSubdivisions;
    const int nSubCellPoints = nCells * subCellPoints.perCell;
    subCellPoints.px.resize(nSubCellPoints);
    subCellPoints.py.resize(nSubCellPoints);
    subCellPoints.weightFraction.resize(nSubCellPoints);
    subCellPoints.pxBin.resize(nSubCellPoints);
    subCellPoints.pyBin.resize(nSubCellPoints);
    subCellPoints.pzBin.resize(nSubCellPoints);
    for (int iCell = 0; iCell < nCells; ++iCell) {
      for (int iSubY = 0; iSubY < nSubdivisions; ++iSubY) {
        for (int iSubPt = 0; iSubPt < nSubdivisions; ++iSubPt) {
          for (int iSubPhi = 0; iSubPhi < nSubdivisions; ++iSubPhi) {
            const int iPoint = iCell * subCellPoints.perCell + (iSubY * nSubdivisions + iSubPt) * nSubdivisions + iSubPhi;
            const double pointY = grid.y[iCell] + ((iSubY + 0.5) / nSubdivisions - 0.5) * grid.yWidth;
            const double pointPt = grid.pt[iCell] + ((iSubPt + 0.5) / nSubdivisions - 0.5) * grid.ptWidth;
            const double pointPhi = grid.phi[iCell] + ((iSubPhi + 0.5) / nSubdivisions - 0.5) * grid.phiWidth;
            const double pointPz = std::sqrt(lambdaMassSq + squared(pointPt)) * std::sinh(pointY);
            subCellPoints.px[iPoint] = pointPt * std::cos(pointPhi);
            subCellPoints.py[iPoint] = pointPt * std::sin(pointPhi);
            subCellPoints.weightFraction[iPoint] = pointPt / grid.pt[iCell] / subCellPoints.perCell;
            subCellPoints.pxBin[iPoint] = transverseMomentumAxis.findBin(subCellPoints.px[iPoint]);
            subCellPoints.pyBin[iPoint] = transverseMomentumAxis.findBin(subCellPoints.py[iPoint]);
            subCellPoints.pzBin[iPoint] = longitudinalMomentumAxis.findBin(pointPz);
            if (subCellPoints.pxBin[iPoint] < 0 || subCellPoints.pyBin[iPoint] < 0 || subCellPoints.pzBin[iPoint] < 0)
              throw std::runtime_error("a sub-cell point falls outside the PolarizationMaps/ axes, which are built to hold the whole grid (internal error)");
          }
        }
      }
    }
    // Title tag of every map: smoothed unless the cells stay on their nodes
    const std::string mapTag = nSubdivisions > 1 ? "full grid, smoothed map" : "full grid";

    // ==================================================================================================================================
    // 5. Bookings, one block per folder
    // ==================================================================================================================================
    TDirectory* configDir = outputFile->mkdir("Config");
    TDirectory* inputsDir = outputFile->mkdir("Inputs");
    TDirectory* paperDir = outputFile->mkdir("Paper");
    TDirectory* ringDir = outputFile->mkdir("Ring");
    TDirectory* noCutsDir = ringDir->mkdir("NoCuts");
    TDirectory* polVsAnglesDir = outputFile->mkdir("PolVsAngles");
    TDirectory* polMapsDir = outputFile->mkdir("PolarizationMaps");
    TDirectory* polMapsLabDir = polMapsDir->mkdir("Lab");
    TDirectory* polMapsPrimeJetDir = polMapsDir->mkdir("PrimeJet");
    TDirectory* polMapsPrimeV0Dir = polMapsDir->mkdir("PrimeV0");
    TDirectory* qaDir = outputFile->mkdir("QA");
    TDirectory* decompositionDir = qaDir->mkdir("Decomposition");
    TDirectory* chiDir = qaDir->mkdir("Chi");
    TDirectory* closureDir = qaDir->mkdir("Closure");
    TDirectory* scanDir = config.deltaScanSteps > 0 ? outputFile->mkdir("DeltaScan") : nullptr;
    TDirectory* mixingDir = config.mixDepth > 0 ? outputFile->mkdir("Mixing") : nullptr;
    TDirectory* smudgeDir = config.smudgeSigma > 0. ? outputFile->mkdir("Smudge") : nullptr;

    TrueTriggerOutputs trueTriggerOutputs;
    DeltaScanOutputs deltaScanOutputs;
    ReducedSetOutputs mixingOutputs, smudgeOutputs;
    MixingExtras mixingExtras;
    SmudgeExtras smudgeExtras;

    // ---------------- Inputs/ ----------------
    inputsDir->cd();
    trueTriggerOutputs.yieldVsPhi.attach(new TH1D("hYieldVsPhi", ("Yield vs #varphi (all cells, all events);" + kPhiLabel + ";#sum w").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.yieldVsPt.attach(new TH1D("hYieldVsPt", ("Yield vs #it{p}_{T} (all cells, all events);" + kPtLabel + ";#sum w").c_str(), ptAxis.nBins, ptAxis.low, ptAxis.high));
    trueTriggerOutputs.yieldVsY.attach(new TH1D("hYieldVsY", ("Yield vs #it{y} (all cells, all events);" + kYLabel + ";#sum w").c_str(), yAxis.nBins, yAxis.low, yAxis.high));
    trueTriggerOutputs.yieldVsDeltaPhiVsPt.attach(new TH2D("h2dYieldVsDeltaPhiVsPt", ("Yield vs #Delta#varphi vs #it{p}_{T} (y cut, all events);" + kDeltaPhiLabel + ";" + kPtLabel).c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high, ptAxis.nBins, ptAxis.low, ptAxis.high));
    trueTriggerOutputs.phiTrigger = new TH1D("hPhiTrigger", "Trigger azimuth;#varphi_{t} (rad);Events", phiTriggerAxis.nBins, phiTriggerAxis.low, phiTriggerAxis.high);

    // ---------------- Ring/ (true trigger) ----------------
    ringDir->cd();
    trueTriggerOutputs.obsVsDeltaPhi[kRing].attach(new TProfile("pRingVsDeltaPhi", ("<#it{R}> vs #Delta#varphi (pooled);" + kDeltaPhiLabel + ";<#it{R}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.obsVsDeltaPhi[kRing].attachEvtMean(new TProfile("pRingVsDeltaPhi_EvtMean", ("<#it{R}> vs #Delta#varphi (event mean);" + kDeltaPhiLabel + ";<#it{R}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.obsVsDeltaPhi[kRingZ].attach(new TProfile("pRingZVsDeltaPhi", ("<#it{R}_{z}> vs #Delta#varphi (pooled);" + kDeltaPhiLabel + ";<#it{R}_{z}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.obsVsDeltaPhi[kRingZ].attachEvtMean(new TProfile("pRingZVsDeltaPhi_EvtMean", ("<#it{R}_{z}> vs #Delta#varphi (event mean);" + kDeltaPhiLabel + ";<#it{R}_{z}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.obsVsDeltaPhi[kRingZT].attach(new TProfile("pRingZTVsDeltaPhi", ("<#it{R}_{z}^{T}> vs #Delta#varphi (pooled);" + kDeltaPhiLabel + ";<#it{R}_{z}^{T}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.obsVsDeltaPhi[kHelicityTerm].attach(new TProfile("pHelicityTermVsDeltaPhi", ("<#it{H}> = <#it{R}_{z} - #it{R}_{z}^{T}> vs #Delta#varphi (pooled);" + kDeltaPhiLabel + ";<#it{H}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.obsVsDeltaPhi[kRingZLab].attach(new TProfile("pRingZLabVsDeltaPhi", ("<#it{R}_{z}^{lab}> vs #Delta#varphi (pooled);" + kDeltaPhiLabel + ";<#it{R}_{z}^{lab}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.obsVsDeltaPhi[kRingPerp].attach(new TProfile("pRingPerpVsDeltaPhi", ("<#it{R}_{#perp}> = <#it{R} - #it{R}_{z}> vs #Delta#varphi (pooled);" + kDeltaPhiLabel + ";<#it{R}_{#perp}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.obsVsDeltaPhi[kNzSq].attach(new TProfile("pNzSqVsDeltaPhi", ("<n_{z}^{2}> vs #Delta#varphi (pooled);" + kDeltaPhiLabel + ";<n_{z}^{2}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.obsVsDeltaPhi[kNzSqSinSqTheta].attach(new TProfile("pNzSqSinSqThetaVsDeltaPhi", ("<n_{z}^{2} sin^{2}#theta> vs #Delta#varphi (pooled);" + kDeltaPhiLabel + ";<n_{z}^{2} sin^{2}#theta>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.obsD3N[kRing] = new TProfile("pRingD3N", "Integrated <#it{R}> (d^{3}N-weighted, pooled);;<#it{R}>", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
    trueTriggerOutputs.obsD3NEvtMean[kRing] = new TProfile("pRingD3N_EvtMean", "Integrated <#it{R}> (d^{3}N-weighted, event mean);;<#it{R}>", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
    trueTriggerOutputs.obsD3N[kRingZ] = new TProfile("pRingZD3N", "Integrated <#it{R}_{z}> (d^{3}N-weighted, pooled);;<#it{R}_{z}>", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
    trueTriggerOutputs.obsD3NEvtMean[kRingZ] = new TProfile("pRingZD3N_EvtMean", "Integrated <#it{R}_{z}> (d^{3}N-weighted, event mean);;<#it{R}_{z}>", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
    trueTriggerOutputs.obsD3N[kRingZT] = new TProfile("pRingZTD3N", "Integrated <#it{R}_{z}^{T}> (d^{3}N-weighted, pooled);;<#it{R}_{z}^{T}>", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
    trueTriggerOutputs.obsD3N[kHelicityTerm] = new TProfile("pHelicityTermD3N", "Integrated <#it{H}> (d^{3}N-weighted, pooled);;<#it{H}>", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
    trueTriggerOutputs.obsD3N[kRingZLab] = new TProfile("pRingZLabD3N", "Integrated <#it{R}_{z}^{lab}> (d^{3}N-weighted, pooled);;<#it{R}_{z}^{lab}>", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
    trueTriggerOutputs.obsD3N[kRingPerp] = new TProfile("pRingPerpD3N", "Integrated <#it{R}_{#perp}> (d^{3}N-weighted, pooled);;<#it{R}_{#perp}>", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
    trueTriggerOutputs.obsD3N[kNzSq] = new TProfile("pNzSqD3N", "Integrated <n_{z}^{2}> (d^{3}N-weighted, pooled);;<n_{z}^{2}>", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
    trueTriggerOutputs.obsD3N[kNzSqSinSqTheta] = new TProfile("pNzSqSinSqThetaD3N", "Integrated <n_{z}^{2} sin^{2}#theta> (d^{3}N-weighted, pooled);;<n_{z}^{2} sin^{2}#theta>", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
    trueTriggerOutputs.ringFlat = new TProfile("pRingFlat", "Integrated #it{R}, flat in #Delta#varphi (events weighted by their in-cut yield);;#it{R}_{flat}", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
    trueTriggerOutputs.ringFlatEvtMean = new TProfile("pRingFlat_EvtMean", "Integrated #it{R}, flat in #Delta#varphi (event mean);;#it{R}_{flat}", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
    trueTriggerOutputs.ringZFlat = new TProfile("pRingZFlat", "Integrated #it{R}_{z}, flat in #Delta#varphi (events weighted by their in-cut yield);;#it{R}_{z,flat}", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
    trueTriggerOutputs.ringZFlatEvtMean = new TProfile("pRingZFlat_EvtMean", "Integrated #it{R}_{z}, flat in #Delta#varphi (event mean);;#it{R}_{z,flat}", integratedAxis.nBins, integratedAxis.low, integratedAxis.high);
    trueTriggerOutputs.ringVsDeltaPhiVsPt.attach(new TProfile2D("p2dRingVsDeltaPhiVsPt", ("<#it{R}> vs #Delta#varphi vs #it{p}_{T} (y cut only);" + kDeltaPhiLabel + ";" + kPtLabel + ";<#it{R}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high, ptAxis.nBins, ptAxis.low, ptAxis.high));
    trueTriggerOutputs.ringZVsDeltaPhiVsPt.attach(new TProfile2D("p2dRingZVsDeltaPhiVsPt", ("<#it{R}_{z}> vs #Delta#varphi vs #it{p}_{T} (y cut only);" + kDeltaPhiLabel + ";" + kPtLabel + ";<#it{R}_{z}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high, ptAxis.nBins, ptAxis.low, ptAxis.high));
    trueTriggerOutputs.ringVsDeltaPhiVsY.attach(new TProfile2D("p2dRingVsDeltaPhiVsY", ("<#it{R}> vs #Delta#varphi vs #it{y} (pT cut only);" + kDeltaPhiLabel + ";" + kYLabel + ";<#it{R}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high, yAxis.nBins, yAxis.low, yAxis.high));
    trueTriggerOutputs.ringZVsDeltaPhiVsY.attach(new TProfile2D("p2dRingZVsDeltaPhiVsY", ("<#it{R}_{z}> vs #Delta#varphi vs #it{y} (pT cut only);" + kDeltaPhiLabel + ";" + kYLabel + ";<#it{R}_{z}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high, yAxis.nBins, yAxis.low, yAxis.high));
    trueTriggerOutputs.ringVsPt.attach(new TProfile("pRingVsPt", ("<#it{R}> vs #it{p}_{T} (y cut only);" + kPtLabel + ";<#it{R}>").c_str(), ptAxis.nBins, ptAxis.low, ptAxis.high));
    trueTriggerOutputs.ringZVsPt.attach(new TProfile("pRingZVsPt", ("<#it{R}_{z}> vs #it{p}_{T} (y cut only);" + kPtLabel + ";<#it{R}_{z}>").c_str(), ptAxis.nBins, ptAxis.low, ptAxis.high));
    trueTriggerOutputs.ringVsY.attach(new TProfile("pRingVsY", ("<#it{R}> vs #it{y} (pT cut only);" + kYLabel + ";<#it{R}>").c_str(), yAxis.nBins, yAxis.low, yAxis.high));
    trueTriggerOutputs.ringZVsY.attach(new TProfile("pRingZVsY", ("<#it{R}_{z}> vs #it{y} (pT cut only);" + kYLabel + ";<#it{R}_{z}>").c_str(), yAxis.nBins, yAxis.low, yAxis.high));
    noCutsDir->cd();
    trueTriggerOutputs.ringVsDeltaPhiNoCuts.attach(new TProfile("pRingVsDeltaPhi", ("<#it{R}> vs #Delta#varphi (full grid, pooled);" + kDeltaPhiLabel + ";<#it{R}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.ringZVsDeltaPhiNoCuts.attach(new TProfile("pRingZVsDeltaPhi", ("<#it{R}_{z}> vs #Delta#varphi (full grid, pooled);" + kDeltaPhiLabel + ";<#it{R}_{z}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));

    // ---------------- PolVsAngles/ ----------------
    polVsAnglesDir->cd();
    trueTriggerOutputs.polVsPhi[kPolLabX].attach(new TProfile("pPolXLabVsPhi", ("<P_{x}> (lab) vs #varphi;" + kPhiLabel + ";<P_{x}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.polVsPhi[kPolLabY].attach(new TProfile("pPolYLabVsPhi", ("<P_{y}> (lab) vs #varphi;" + kPhiLabel + ";<P_{y}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.polVsPhi[kPolLabZ].attach(new TProfile("pPolZLabVsPhi", ("<P_{z}> (lab) vs #varphi;" + kPhiLabel + ";<P_{z}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.polVsPhi[kPolStarX].attach(new TProfile("pPolXStarVsPhi", ("<P*_{x}> vs #varphi;" + kPhiLabel + ";<P*_{x}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.polVsPhi[kPolStarY].attach(new TProfile("pPolYStarVsPhi", ("<P*_{y}> vs #varphi;" + kPhiLabel + ";<P*_{y}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.polVsPhi[kPolStarZ].attach(new TProfile("pPolZStarVsPhi", ("<P*_{z}> vs #varphi;" + kPhiLabel + ";<P*_{z}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.polVsPhiVsPt[kPolLabX].attach(new TProfile2D("p2dPolXLabVsPhiVsPt", ("<P_{x}> (lab) vs #varphi vs #it{p}_{T};" + kPhiLabel + ";" + kPtLabel + ";<P_{x}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high, ptAxis.nBins, ptAxis.low, ptAxis.high));
    trueTriggerOutputs.polVsPhiVsPt[kPolLabY].attach(new TProfile2D("p2dPolYLabVsPhiVsPt", ("<P_{y}> (lab) vs #varphi vs #it{p}_{T};" + kPhiLabel + ";" + kPtLabel + ";<P_{y}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high, ptAxis.nBins, ptAxis.low, ptAxis.high));
    trueTriggerOutputs.polVsPhiVsPt[kPolLabZ].attach(new TProfile2D("p2dPolZLabVsPhiVsPt", ("<P_{z}> (lab) vs #varphi vs #it{p}_{T};" + kPhiLabel + ";" + kPtLabel + ";<P_{z}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high, ptAxis.nBins, ptAxis.low, ptAxis.high));
    trueTriggerOutputs.polVsPhiVsPt[kPolStarX].attach(new TProfile2D("p2dPolXStarVsPhiVsPt", ("<P*_{x}> vs #varphi vs #it{p}_{T};" + kPhiLabel + ";" + kPtLabel + ";<P*_{x}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high, ptAxis.nBins, ptAxis.low, ptAxis.high));
    trueTriggerOutputs.polVsPhiVsPt[kPolStarY].attach(new TProfile2D("p2dPolYStarVsPhiVsPt", ("<P*_{y}> vs #varphi vs #it{p}_{T};" + kPhiLabel + ";" + kPtLabel + ";<P*_{y}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high, ptAxis.nBins, ptAxis.low, ptAxis.high));
    trueTriggerOutputs.polVsPhiVsPt[kPolStarZ].attach(new TProfile2D("p2dPolZStarVsPhiVsPt", ("<P*_{z}> vs #varphi vs #it{p}_{T};" + kPhiLabel + ";" + kPtLabel + ";<P*_{z}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high, ptAxis.nBins, ptAxis.low, ptAxis.high));
    trueTriggerOutputs.polJetVsDeltaPhi[kPolLabX].attach(new TProfile("pPolXJetLabVsDeltaPhi", ("<P_{x'}> (lab, x' along t) vs #Delta#varphi;" + kDeltaPhiLabel + ";<P_{x'}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.polJetVsDeltaPhi[kPolLabY].attach(new TProfile("pPolYJetLabVsDeltaPhi", ("<P_{y'}> (lab, y' = z #times x') vs #Delta#varphi;" + kDeltaPhiLabel + ";<P_{y'}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.polJetVsDeltaPhi[kPolLabZ].attach(new TProfile("pPolZLabVsDeltaPhi", ("<P_{z}> (lab) vs #Delta#varphi;" + kDeltaPhiLabel + ";<P_{z}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.polJetVsDeltaPhi[kPolStarX].attach(new TProfile("pPolXJetStarVsDeltaPhi", ("<P*_{x'}> (x' along t) vs #Delta#varphi;" + kDeltaPhiLabel + ";<P*_{x'}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.polJetVsDeltaPhi[kPolStarY].attach(new TProfile("pPolYJetStarVsDeltaPhi", ("<P*_{y'}> (y' = z #times x') vs #Delta#varphi;" + kDeltaPhiLabel + ";<P*_{y'}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.polJetVsDeltaPhi[kPolStarZ].attach(new TProfile("pPolZStarVsDeltaPhi", ("<P*_{z}> vs #Delta#varphi;" + kDeltaPhiLabel + ";<P*_{z}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.polJetVsDeltaPhiVsPt[kPolLabX].attach(new TProfile2D("p2dPolXJetLabVsDeltaPhiVsPt", ("<P_{x'}> (lab, x' along t) vs #Delta#varphi vs #it{p}_{T};" + kDeltaPhiLabel + ";" + kPtLabel + ";<P_{x'}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high, ptAxis.nBins, ptAxis.low, ptAxis.high));
    trueTriggerOutputs.polJetVsDeltaPhiVsPt[kPolLabY].attach(new TProfile2D("p2dPolYJetLabVsDeltaPhiVsPt", ("<P_{y'}> (lab, y' = z #times x') vs #Delta#varphi vs #it{p}_{T};" + kDeltaPhiLabel + ";" + kPtLabel + ";<P_{y'}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high, ptAxis.nBins, ptAxis.low, ptAxis.high));
    trueTriggerOutputs.polJetVsDeltaPhiVsPt[kPolLabZ].attach(new TProfile2D("p2dPolZLabVsDeltaPhiVsPt", ("<P_{z}> (lab) vs #Delta#varphi vs #it{p}_{T};" + kDeltaPhiLabel + ";" + kPtLabel + ";<P_{z}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high, ptAxis.nBins, ptAxis.low, ptAxis.high));
    trueTriggerOutputs.polJetVsDeltaPhiVsPt[kPolStarX].attach(new TProfile2D("p2dPolXJetStarVsDeltaPhiVsPt", ("<P*_{x'}> (x' along t) vs #Delta#varphi vs #it{p}_{T};" + kDeltaPhiLabel + ";" + kPtLabel + ";<P*_{x'}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high, ptAxis.nBins, ptAxis.low, ptAxis.high));
    trueTriggerOutputs.polJetVsDeltaPhiVsPt[kPolStarY].attach(new TProfile2D("p2dPolYJetStarVsDeltaPhiVsPt", ("<P*_{y'}> (y' = z #times x') vs #Delta#varphi vs #it{p}_{T};" + kDeltaPhiLabel + ";" + kPtLabel + ";<P*_{y'}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high, ptAxis.nBins, ptAxis.low, ptAxis.high));
    trueTriggerOutputs.polJetVsDeltaPhiVsPt[kPolStarZ].attach(new TProfile2D("p2dPolZStarVsDeltaPhiVsPt", ("<P*_{z}> vs #Delta#varphi vs #it{p}_{T};" + kDeltaPhiLabel + ";" + kPtLabel + ";<P*_{z}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high, ptAxis.nBins, ptAxis.low, ptAxis.high));

    // ---------------- PolarizationMaps/ (README, "Polarization maps") ----------------
    // Lab coordinate system
    polMapsLabDir->cd();
    trueTriggerOutputs.labMapPxPy[kMapX].attach(new TProfile2D("p2dPxStar_vsPxPy", ("<P*_{x}> vs (#it{p}_{x}, #it{p}_{y}) (" + mapTag + ");" + kPxLabel + ";" + kPyLabel + ";<P*_{x}>").c_str(), transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high, transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high));
    trueTriggerOutputs.labMapPxPy[kMapY].attach(new TProfile2D("p2dPyStar_vsPxPy", ("<P*_{y}> vs (#it{p}_{x}, #it{p}_{y}) (" + mapTag + ");" + kPxLabel + ";" + kPyLabel + ";<P*_{y}>").c_str(), transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high, transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high));
    trueTriggerOutputs.labMapPxPy[kMapZ].attach(new TProfile2D("p2dPzStar_vsPxPy", ("<P*_{z}> vs (#it{p}_{x}, #it{p}_{y}) (" + mapTag + ");" + kPxLabel + ";" + kPyLabel + ";<P*_{z}>").c_str(), transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high, transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high));
    trueTriggerOutputs.labYieldPxPy.attach(new TH2D("h2dYieldVsPxPy", ("Yield vs (#it{p}_{x}, #it{p}_{y}) (" + mapTag + ");" + kPxLabel + ";" + kPyLabel + ";#sum w").c_str(), transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high, transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high));
    trueTriggerOutputs.labMapPzPx[kMapX].attach(new TProfile2D("p2dPxStar_vsPzPx", ("<P*_{x}> vs (#it{p}_{z}, #it{p}_{x}) (" + mapTag + ");" + kPzLabel + ";" + kPxLabel + ";<P*_{x}>").c_str(), longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high, transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high));
    trueTriggerOutputs.labMapPzPx[kMapY].attach(new TProfile2D("p2dPyStar_vsPzPx", ("<P*_{y}> vs (#it{p}_{z}, #it{p}_{x}) (" + mapTag + ");" + kPzLabel + ";" + kPxLabel + ";<P*_{y}>").c_str(), longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high, transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high));
    trueTriggerOutputs.labMapPzPx[kMapZ].attach(new TProfile2D("p2dPzStar_vsPzPx", ("<P*_{z}> vs (#it{p}_{z}, #it{p}_{x}) (" + mapTag + ");" + kPzLabel + ";" + kPxLabel + ";<P*_{z}>").c_str(), longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high, transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high));
    trueTriggerOutputs.labYieldPzPx.attach(new TH2D("h2dYieldVsPzPx", ("Yield vs (#it{p}_{z}, #it{p}_{x}) (" + mapTag + ");" + kPzLabel + ";" + kPxLabel + ";#sum w").c_str(), longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high, transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high));
    trueTriggerOutputs.labMapPyPz[kMapX].attach(new TProfile2D("p2dPxStar_vsPyPz", ("<P*_{x}> vs (#it{p}_{y}, #it{p}_{z}) (" + mapTag + ");" + kPyLabel + ";" + kPzLabel + ";<P*_{x}>").c_str(), transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high, longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high));
    trueTriggerOutputs.labMapPyPz[kMapY].attach(new TProfile2D("p2dPyStar_vsPyPz", ("<P*_{y}> vs (#it{p}_{y}, #it{p}_{z}) (" + mapTag + ");" + kPyLabel + ";" + kPzLabel + ";<P*_{y}>").c_str(), transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high, longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high));
    trueTriggerOutputs.labMapPyPz[kMapZ].attach(new TProfile2D("p2dPzStar_vsPyPz", ("<P*_{z}> vs (#it{p}_{y}, #it{p}_{z}) (" + mapTag + ");" + kPyLabel + ";" + kPzLabel + ";<P*_{z}>").c_str(), transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high, longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high));
    trueTriggerOutputs.labYieldPyPz.attach(new TH2D("h2dYieldVsPyPz", ("Yield vs (#it{p}_{y}, #it{p}_{z}) (" + mapTag + ");" + kPyLabel + ";" + kPzLabel + ";#sum w").c_str(), transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high, longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high));
    // PrimeJet: the consumer's jet coordinate system, z' = t, x' = beam projected perpendicular to t, y' = z' x x'
    polMapsPrimeJetDir->cd();
    trueTriggerOutputs.primeJetMapPxPy[kMapX].attach(new TProfile2D("p2dPxStarPrimeJet_vsPxPyPrimeJet", ("<P*_{x'Jet}> vs (#it{p}_{x'Jet}, #it{p}_{y'Jet}) (" + mapTag + ");" + kPxPrimeJetLabel + ";" + kPyPrimeJetLabel + ";<P*_{x'Jet}>").c_str(), longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high, transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high));
    trueTriggerOutputs.primeJetMapPxPy[kMapY].attach(new TProfile2D("p2dPyStarPrimeJet_vsPxPyPrimeJet", ("<P*_{y'Jet}> vs (#it{p}_{x'Jet}, #it{p}_{y'Jet}) (" + mapTag + ");" + kPxPrimeJetLabel + ";" + kPyPrimeJetLabel + ";<P*_{y'Jet}>").c_str(), longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high, transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high));
    trueTriggerOutputs.primeJetMapPxPy[kMapZ].attach(new TProfile2D("p2dPzStarPrimeJet_vsPxPyPrimeJet", ("<P*_{z'Jet}> vs (#it{p}_{x'Jet}, #it{p}_{y'Jet}) (" + mapTag + ");" + kPxPrimeJetLabel + ";" + kPyPrimeJetLabel + ";<P*_{z'Jet}>").c_str(), longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high, transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high));
    trueTriggerOutputs.primeJetYieldPxPy.attach(new TH2D("h2dYieldVsPxPyPrimeJet", ("Yield vs (#it{p}_{x'Jet}, #it{p}_{y'Jet}) (" + mapTag + ");" + kPxPrimeJetLabel + ";" + kPyPrimeJetLabel + ";#sum w").c_str(), longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high, transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high));
    trueTriggerOutputs.primeJetRingPxPy[kMapRing].attach(new TProfile2D("p2dRingVsPxPyPrimeJet", ("<#it{R}> vs (#it{p}_{x'Jet}, #it{p}_{y'Jet}) (" + mapTag + ");" + kPxPrimeJetLabel + ";" + kPyPrimeJetLabel + ";<#it{R}>").c_str(), longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high, transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high));
    trueTriggerOutputs.primeJetRingPxPy[kMapRingZ].attach(new TProfile2D("p2dRingZVsPxPyPrimeJet", ("<#it{R}_{z}> vs (#it{p}_{x'Jet}, #it{p}_{y'Jet}) (" + mapTag + ");" + kPxPrimeJetLabel + ";" + kPyPrimeJetLabel + ";<#it{R}_{z}>").c_str(), longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high, transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high));
    trueTriggerOutputs.primeJetMapPzPx[kMapX].attach(new TProfile2D("p2dPxStarPrimeJet_vsPzPxPrimeJet", ("<P*_{x'Jet}> vs (#it{p}_{z'Jet}, #it{p}_{x'Jet}) (" + mapTag + ");" + kPzPrimeJetLabel + ";" + kPxPrimeJetLabel + ";<P*_{x'Jet}>").c_str(), transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high, longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high));
    trueTriggerOutputs.primeJetMapPzPx[kMapY].attach(new TProfile2D("p2dPyStarPrimeJet_vsPzPxPrimeJet", ("<P*_{y'Jet}> vs (#it{p}_{z'Jet}, #it{p}_{x'Jet}) (" + mapTag + ");" + kPzPrimeJetLabel + ";" + kPxPrimeJetLabel + ";<P*_{y'Jet}>").c_str(), transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high, longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high));
    trueTriggerOutputs.primeJetMapPzPx[kMapZ].attach(new TProfile2D("p2dPzStarPrimeJet_vsPzPxPrimeJet", ("<P*_{z'Jet}> vs (#it{p}_{z'Jet}, #it{p}_{x'Jet}) (" + mapTag + ");" + kPzPrimeJetLabel + ";" + kPxPrimeJetLabel + ";<P*_{z'Jet}>").c_str(), transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high, longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high));
    trueTriggerOutputs.primeJetYieldPzPx.attach(new TH2D("h2dYieldVsPzPxPrimeJet", ("Yield vs (#it{p}_{z'Jet}, #it{p}_{x'Jet}) (" + mapTag + ");" + kPzPrimeJetLabel + ";" + kPxPrimeJetLabel + ";#sum w").c_str(), transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high, longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high));
    trueTriggerOutputs.primeJetRingPzPx[kMapRing].attach(new TProfile2D("p2dRingVsPzPxPrimeJet", ("<#it{R}> vs (#it{p}_{z'Jet}, #it{p}_{x'Jet}) (" + mapTag + ");" + kPzPrimeJetLabel + ";" + kPxPrimeJetLabel + ";<#it{R}>").c_str(), transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high, longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high));
    trueTriggerOutputs.primeJetRingPzPx[kMapRingZ].attach(new TProfile2D("p2dRingZVsPzPxPrimeJet", ("<#it{R}_{z}> vs (#it{p}_{z'Jet}, #it{p}_{x'Jet}) (" + mapTag + ");" + kPzPrimeJetLabel + ";" + kPxPrimeJetLabel + ";<#it{R}_{z}>").c_str(), transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high, longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high));
    trueTriggerOutputs.primeJetMapPyPz[kMapX].attach(new TProfile2D("p2dPxStarPrimeJet_vsPyPzPrimeJet", ("<P*_{x'Jet}> vs (#it{p}_{y'Jet}, #it{p}_{z'Jet}) (" + mapTag + ");" + kPyPrimeJetLabel + ";" + kPzPrimeJetLabel + ";<P*_{x'Jet}>").c_str(), transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high, transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high));
    trueTriggerOutputs.primeJetMapPyPz[kMapY].attach(new TProfile2D("p2dPyStarPrimeJet_vsPyPzPrimeJet", ("<P*_{y'Jet}> vs (#it{p}_{y'Jet}, #it{p}_{z'Jet}) (" + mapTag + ");" + kPyPrimeJetLabel + ";" + kPzPrimeJetLabel + ";<P*_{y'Jet}>").c_str(), transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high, transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high));
    trueTriggerOutputs.primeJetMapPyPz[kMapZ].attach(new TProfile2D("p2dPzStarPrimeJet_vsPyPzPrimeJet", ("<P*_{z'Jet}> vs (#it{p}_{y'Jet}, #it{p}_{z'Jet}) (" + mapTag + ");" + kPyPrimeJetLabel + ";" + kPzPrimeJetLabel + ";<P*_{z'Jet}>").c_str(), transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high, transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high));
    trueTriggerOutputs.primeJetYieldPyPz.attach(new TH2D("h2dYieldVsPyPzPrimeJet", ("Yield vs (#it{p}_{y'Jet}, #it{p}_{z'Jet}) (" + mapTag + ");" + kPyPrimeJetLabel + ";" + kPzPrimeJetLabel + ";#sum w").c_str(), transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high, transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high));
    trueTriggerOutputs.primeJetRingPyPz[kMapRing].attach(new TProfile2D("p2dRingVsPyPzPrimeJet", ("<#it{R}> vs (#it{p}_{y'Jet}, #it{p}_{z'Jet}) (" + mapTag + ");" + kPyPrimeJetLabel + ";" + kPzPrimeJetLabel + ";<#it{R}>").c_str(), transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high, transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high));
    trueTriggerOutputs.primeJetRingPyPz[kMapRingZ].attach(new TProfile2D("p2dRingZVsPyPzPrimeJet", ("<#it{R}_{z}> vs (#it{p}_{y'Jet}, #it{p}_{z'Jet}) (" + mapTag + ");" + kPyPrimeJetLabel + ";" + kPzPrimeJetLabel + ";<#it{R}_{z}>").c_str(), transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high, transverseMomentumAxis.nBins, transverseMomentumAxis.low, transverseMomentumAxis.high));
    // PrimeV0: the Lambda production plane (p_z, p_T), x'V0 along p_T and y'V0 = z x x'V0, the azimuthal direction
    polMapsPrimeV0Dir->cd();
    trueTriggerOutputs.primeV0MapPzPt[kMapX].attach(new TProfile2D("p2dPxStarPrimeV0_vsPzPt", ("<P*_{x'V0}> vs (#it{p}_{z}, #it{p}_{T}) (" + mapTag + ");" + kPzLabel + ";" + kPtLabel + ";<P*_{x'V0}>").c_str(), longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high, ptAxis.nBins, ptAxis.low, ptAxis.high));
    trueTriggerOutputs.primeV0MapPzPt[kMapY].attach(new TProfile2D("p2dPyStarPrimeV0_vsPzPt", ("<P*_{y'V0}> vs (#it{p}_{z}, #it{p}_{T}) (" + mapTag + ");" + kPzLabel + ";" + kPtLabel + ";<P*_{y'V0}>").c_str(), longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high, ptAxis.nBins, ptAxis.low, ptAxis.high));
    trueTriggerOutputs.primeV0MapPzPt[kMapZ].attach(new TProfile2D("p2dPzStar_vsPzPt", ("<P*_{z}> vs (#it{p}_{z}, #it{p}_{T}) (" + mapTag + ");" + kPzLabel + ";" + kPtLabel + ";<P*_{z}>").c_str(), longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high, ptAxis.nBins, ptAxis.low, ptAxis.high));
    trueTriggerOutputs.primeV0YieldPzPt.attach(new TH2D("h2dYieldVsPzPt", ("Yield vs (#it{p}_{z}, #it{p}_{T}) (" + mapTag + ");" + kPzLabel + ";" + kPtLabel + ";#sum w").c_str(), longitudinalMomentumAxis.nBins, longitudinalMomentumAxis.low, longitudinalMomentumAxis.high, ptAxis.nBins, ptAxis.low, ptAxis.high));

    // ---------------- QA/ ----------------
    decompositionDir->cd();
    trueTriggerOutputs.polEVsDeltaPhi.attach(new TProfile("pPolEVsDeltaPhi", ("<P_{e}> = <#vec{P}*#upoint(#hat{p} #times #hat{n})> vs #Delta#varphi;" + kDeltaPhiLabel + ";<P_{e}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.polParVsDeltaPhi.attach(new TProfile("pPolParVsDeltaPhi", ("<P_{#parallel}> = <#vec{P}*#upoint#hat{p}> vs #Delta#varphi;" + kDeltaPhiLabel + ";<P_{#parallel}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    trueTriggerOutputs.alignmentVsDeltaPhi.attach(new TProfile("pAlignmentVsDeltaPhi", ("<P_{R}/|#vec{P}*|> vs #Delta#varphi;" + kDeltaPhiLabel + ";<P_{R}/|#vec{P}*|>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high));
    chiDir->cd();
    trueTriggerOutputs.ringVsChi.attach(new TProfile("pRingVsChi", "<#it{R}> vs #chi (#hat{n} = cos#chi #hat{#varphi} + sin#chi #hat{#theta});#chi (rad);<#it{R}>", chiAxis.nBins, chiAxis.low, chiAxis.high));
    trueTriggerOutputs.ringZVsChi.attach(new TProfile("pRingZVsChi", "<#it{R}_{z}> vs #chi (#hat{n} = cos#chi #hat{#varphi} + sin#chi #hat{#theta});#chi (rad);<#it{R}_{z}>", chiAxis.nBins, chiAxis.low, chiAxis.high));
    trueTriggerOutputs.chiYield.attach(new TH1D("hChi", "Yield vs #chi;#chi (rad);#sum w", chiAxis.nBins, chiAxis.low, chiAxis.high));
    closureDir->cd();
    trueTriggerOutputs.closurePdotP.attach(new TH1D("hPdotP", "Lab P#upoint p / (E |#vec{P}|) per cell;log_{10}(relative residual);Cells", logResidualAxis.nBins, logResidualAxis.low, logResidualAxis.high));
    trueTriggerOutputs.closurePStarTime.attach(new TH1D("hPStarTime", "P*^{0} / |#vec{P}*| per cell (LRF block);log_{10}(relative residual);Cells", logResidualAxis.nBins, logResidualAxis.low, logResidualAxis.high));
    trueTriggerOutputs.closureBoost.attach(new TH1D("hBoostClosure", "|#vec{P}*_{our boost} - #vec{P}*_{LRF}| / |#vec{P}*_{LRF}| per cell;log_{10}(relative residual);Cells", logResidualAxis.nBins, logResidualAxis.low, logResidualAxis.high));
    trueTriggerOutputs.closureRingLabMinusStar.attach(new TH1D("hRingLabMinusStar", "|#vec{P}#upoint#hat{n} - #vec{P}*#upoint#hat{n}| / |#vec{P}*| per cell;log_{10}(relative residual);Cells", logResidualAxis.nBins, logResidualAxis.low, logResidualAxis.high));

    // ---------------- DeltaScan/ ----------------
    ReducedAccumulators scanAccumulators, mixingAccumulators, smudgeAccumulators;
    if (config.deltaScanSteps > 0) {
      scanDir->cd();
      deltaScanOutputs.ringD3NVsDelta = new TProfile("pRingD3NVsDelta", ("Integrated <#it{R}> vs #delta (d^{3}N-weighted, pooled);" + kDeltaLabel + ";<#it{R}>").c_str(), deltaAxis.nBins, deltaAxis.low, deltaAxis.high);
      deltaScanOutputs.ringD3NVsDeltaEvtMean = new TProfile("pRingD3NVsDelta_EvtMean", ("Integrated <#it{R}> vs #delta (d^{3}N-weighted, event mean);" + kDeltaLabel + ";<#it{R}>").c_str(), deltaAxis.nBins, deltaAxis.low, deltaAxis.high);
      deltaScanOutputs.ringFlatVsDelta = new TProfile("pRingFlatVsDelta", ("Integrated #it{R} vs #delta (flat, events weighted by their in-cut yield);" + kDeltaLabel + ";#it{R}_{flat}").c_str(), deltaAxis.nBins, deltaAxis.low, deltaAxis.high);
      deltaScanOutputs.ringFlatVsDeltaEvtMean = new TProfile("pRingFlatVsDelta_EvtMean", ("Integrated #it{R} vs #delta (flat, event mean);" + kDeltaLabel + ";#it{R}_{flat}").c_str(), deltaAxis.nBins, deltaAxis.low, deltaAxis.high);
      deltaScanOutputs.ringZD3NVsDelta = new TProfile("pRingZD3NVsDelta", ("Integrated <#it{R}_{z}> vs #delta (d^{3}N-weighted, pooled);" + kDeltaLabel + ";<#it{R}_{z}>").c_str(), deltaAxis.nBins, deltaAxis.low, deltaAxis.high);
      deltaScanOutputs.ringZD3NVsDeltaEvtMean = new TProfile("pRingZD3NVsDelta_EvtMean", ("Integrated <#it{R}_{z}> vs #delta (d^{3}N-weighted, event mean);" + kDeltaLabel + ";<#it{R}_{z}>").c_str(), deltaAxis.nBins, deltaAxis.low, deltaAxis.high);
      deltaScanOutputs.ringZFlatVsDelta = new TProfile("pRingZFlatVsDelta", ("Integrated #it{R}_{z} vs #delta (flat, events weighted by their in-cut yield);" + kDeltaLabel + ";#it{R}_{z,flat}").c_str(), deltaAxis.nBins, deltaAxis.low, deltaAxis.high);
      deltaScanOutputs.ringZFlatVsDeltaEvtMean = new TProfile("pRingZFlatVsDelta_EvtMean", ("Integrated #it{R}_{z} vs #delta (flat, event mean);" + kDeltaLabel + ";#it{R}_{z,flat}").c_str(), deltaAxis.nBins, deltaAxis.low, deltaAxis.high);
      deltaScanOutputs.nzSqD3NVsDelta = new TProfile("pNzSqD3NVsDelta", ("Integrated <n_{z}^{2}> vs #delta (d^{3}N-weighted, pooled);" + kDeltaLabel + ";<n_{z}^{2}>").c_str(), deltaAxis.nBins, deltaAxis.low, deltaAxis.high);
      deltaScanOutputs.nzSqD3NVsDeltaEvtMean = new TProfile("pNzSqD3NVsDelta_EvtMean", ("Integrated <n_{z}^{2}> vs #delta (d^{3}N-weighted, event mean);" + kDeltaLabel + ";<n_{z}^{2}>").c_str(), deltaAxis.nBins, deltaAxis.low, deltaAxis.high);
      deltaScanOutputs.ringVsDeltaPhiVsDelta = new TProfile2D("p2dRingVsDeltaPhiVsDelta", ("<#it{R}> vs #Delta#varphi' vs #delta (pooled);#Delta#varphi' from the rotated trigger (rad);" + kDeltaLabel + ";<#it{R}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high, deltaAxis.nBins, deltaAxis.low, deltaAxis.high);
      deltaScanOutputs.ringZVsDeltaPhiVsDelta = new TProfile2D("p2dRingZVsDeltaPhiVsDelta", ("<#it{R}_{z}> vs #Delta#varphi' vs #delta (pooled);#Delta#varphi' from the rotated trigger (rad);" + kDeltaLabel + ";<#it{R}_{z}>").c_str(), deltaPhiAxis.nBins, deltaPhiAxis.low, deltaPhiAxis.high, deltaAxis.nBins, deltaAxis.low, deltaAxis.high);
      deltaScanOutputs.closureHalfTurn = new TH1D("hClosureHalfTurn", "max_{k} |R(#delta_{k} + #pi) + R(#delta_{k})| / max_{k} |R(#delta_{k})| per event;log_{10}(relative residual);Events", logResidualAxis.nBins, logResidualAxis.low, logResidualAxis.high);
      deltaScanOutputs.cosZetaD3NVsDelta = new TProfile("pCosZetaD3NVsDelta", ("Integrated <cos#zeta> = <#hat{n}_{t}#upoint#hat{n}> vs #delta (d^{3}N-weighted, pooled);" + kDeltaLabel + ";<cos#zeta>").c_str(), deltaAxis.nBins, deltaAxis.low, deltaAxis.high);
      deltaScanOutputs.sinZetaD3NVsDelta = new TProfile("pSinZetaD3NVsDelta", ("Integrated <sin#zeta> = <(#hat{n}_{t} #times #hat{n})#upoint#hat{p}> vs #delta (d^{3}N-weighted, pooled);" + kDeltaLabel + ";<sin#zeta>").c_str(), deltaAxis.nBins, deltaAxis.low, deltaAxis.high);
      deltaScanOutputs.cosZetaVsDelta.attach(new TH2D("h2dCosZetaVsDelta", ("cos#zeta vs #delta (d^{3}N-weighted, all events);" + kDeltaLabel + ";cos#zeta").c_str(), deltaAxis.nBins, deltaAxis.low, deltaAxis.high, zetaAxis.nBins, zetaAxis.low, zetaAxis.high));
      deltaScanOutputs.sinZetaVsDelta.attach(new TH2D("h2dSinZetaVsDelta", ("sin#zeta vs #delta (d^{3}N-weighted, all events);" + kDeltaLabel + ";sin#zeta").c_str(), deltaAxis.nBins, deltaAxis.low, deltaAxis.high, zetaAxis.nBins, zetaAxis.low, zetaAxis.high));
      scanAccumulators.resize(grid.nPhi);
    }

    // ---------------- Mixing/ ----------------
    if (config.mixDepth > 0) {
      bookReducedSet(mixingDir, mixingOutputs, deltaPhiAxis, integratedAxis, zetaAxis);
      mixingDir->cd();
      mixingExtras.ringD3NVsDeltaEff = new TProfile("pRingD3NVsDeltaEff", "Integrated <#it{R}> vs #delta_{eff} = #varphi_{J}(partner) - #varphi_{J} (pooled);#delta_{eff} (rad);<#it{R}>", deltaAxis.nBins, deltaAxis.low, deltaAxis.high);
      mixingExtras.ringD3NVsDeltaEffEvtMean = new TProfile("pRingD3NVsDeltaEff_EvtMean", "Integrated <#it{R}> vs #delta_{eff} (event mean);#delta_{eff} (rad);<#it{R}>", deltaAxis.nBins, deltaAxis.low, deltaAxis.high);
      mixingExtras.ringZD3NVsDeltaEff = new TProfile("pRingZD3NVsDeltaEff", "Integrated <#it{R}_{z}> vs #delta_{eff} = #varphi_{J}(partner) - #varphi_{J} (pooled);#delta_{eff} (rad);<#it{R}_{z}>", deltaAxis.nBins, deltaAxis.low, deltaAxis.high);
      mixingExtras.ringZD3NVsDeltaEffEvtMean = new TProfile("pRingZD3NVsDeltaEff_EvtMean", "Integrated <#it{R}_{z}> vs #delta_{eff} (event mean);#delta_{eff} (rad);<#it{R}_{z}>", deltaAxis.nBins, deltaAxis.low, deltaAxis.high);
      mixingAccumulators.resize(grid.nPhi);
    }

    // ---------------- Smudge/ ----------------
    if (config.smudgeSigma > 0.) {
      bookReducedSet(smudgeDir, smudgeOutputs, deltaPhiAxis, integratedAxis, zetaAxis);
      smudgeDir->cd();
      smudgeExtras.smudgeAngle = new TH1D("hSmudgeAngle", "Smudge angle;angle between smudged and true trigger (rad);Draws", smudgeAngleAxis.nBins, smudgeAngleAxis.low, smudgeAngleAxis.high);
      smudgeExtras.triggerEta = new TH1D("hTriggerEta", "Pseudorapidity of the smudged trigger;#eta_{t};Draws", triggerEtaAxis.nBins, triggerEtaAxis.low, triggerEtaAxis.high);
      smudgeExtras.ringD3NVsSmudgeAngle = new TProfile("pRingD3NVsSmudgeAngle", "Integrated <#it{R}> vs smudge angle (pooled);smudge angle (rad);<#it{R}>", smudgeAngleAxis.nBins, smudgeAngleAxis.low, smudgeAngleAxis.high);
      smudgeExtras.ringD3NVsSmudgeAngleEvtMean = new TProfile("pRingD3NVsSmudgeAngle_EvtMean", "Integrated <#it{R}> vs smudge angle (event mean);smudge angle (rad);<#it{R}>", smudgeAngleAxis.nBins, smudgeAngleAxis.low, smudgeAngleAxis.high);
      smudgeExtras.ringZD3NVsSmudgeAngle = new TProfile("pRingZD3NVsSmudgeAngle", "Integrated <#it{R}_{z}> vs smudge angle (pooled);smudge angle (rad);<#it{R}_{z}>", smudgeAngleAxis.nBins, smudgeAngleAxis.low, smudgeAngleAxis.high);
      smudgeExtras.ringZD3NVsSmudgeAngleEvtMean = new TProfile("pRingZD3NVsSmudgeAngle_EvtMean", "Integrated <#it{R}_{z}> vs smudge angle (event mean);smudge angle (rad);<#it{R}_{z}>", smudgeAngleAxis.nBins, smudgeAngleAxis.low, smudgeAngleAxis.high);
      smudgeAccumulators.resize(grid.nPhi);
    }

    // ==================================================================================================================================
    // 6. Event tree (README, "Event tree")
    // ==================================================================================================================================
    outputFile->cd();
    TTree* eventTree = new TTree("events", "One entry per event (README, Event tree)");
    EventRecord eventRecord;
    eventRecord.nDeltaPhiBins = grid.nPhi;
    eventRecord.nDelta = config.deltaScanSteps;
    // Sized once, never reallocated afterwards: the branches keep these addresses (at least one element, so that data() is never null)
    const int nScanSlots = std::max(config.deltaScanSteps, 1);
    eventRecord.sumWDeltaPhi.assign(grid.nPhi, 0.);
    eventRecord.sumARingDeltaPhi.assign(grid.nPhi, 0.);
    eventRecord.sumARingZDeltaPhi.assign(grid.nPhi, 0.);
    eventRecord.sumANzSqDeltaPhi.assign(grid.nPhi, 0.);
    eventRecord.scanRingD3N.assign(nScanSlots, kNoValue);
    eventRecord.scanRingFlat.assign(nScanSlots, kNoValue);
    eventRecord.scanRingZD3N.assign(nScanSlots, kNoValue);
    eventRecord.scanRingZFlat.assign(nScanSlots, kNoValue);
    eventRecord.scanNzSqD3N.assign(nScanSlots, kNoValue);
    eventRecord.scanCosZetaD3N.assign(nScanSlots, kNoValue);
    eventRecord.scanSinZetaD3N.assign(nScanSlots, kNoValue);
    eventRecord.nMixPartners = config.mixDepth;
    eventRecord.mixDeltaEffBin.assign(std::max(1, config.mixDepth), -1); // one slot at least, so that the branch has an address
    eventRecord.nMixDeltaPhiBins = config.mixDepth > 0 ? grid.nPhi : 0;
    eventRecord.mixSumWDeltaPhi.assign(grid.nPhi, 0.);
    eventRecord.mixSumARingDeltaPhi.assign(grid.nPhi, 0.);
    eventRecord.mixSumARingZDeltaPhi.assign(grid.nPhi, 0.);
    eventRecord.mixSumANzSqDeltaPhi.assign(grid.nPhi, 0.);
    eventRecord.mixSumACosZetaDeltaPhi.assign(grid.nPhi, 0.);
    eventRecord.mixSumASinZetaDeltaPhi.assign(grid.nPhi, 0.);
    eventRecord.nSmudgeDeltaPhiBins = config.smudgeSigma > 0. ? grid.nPhi : 0;
    eventRecord.smudgeSumWDeltaPhi.assign(grid.nPhi, 0.);
    eventRecord.smudgeSumARingDeltaPhi.assign(grid.nPhi, 0.);
    eventRecord.smudgeSumARingZDeltaPhi.assign(grid.nPhi, 0.);
    eventRecord.smudgeSumANzSqDeltaPhi.assign(grid.nPhi, 0.);
    eventRecord.smudgeSumACosZetaDeltaPhi.assign(grid.nPhi, 0.);
    eventRecord.smudgeSumASinZetaDeltaPhi.assign(grid.nPhi, 0.);

    // Event and trigger
    eventTree->Branch("eventId", &eventRecord.eventId, "eventId/I");
    eventTree->Branch("phiJ", &eventRecord.phiJ, "phiJ/D");
    eventTree->Branch("phiTrigger", &eventRecord.phiTrigger, "phiTrigger/D");
    eventTree->Branch("nCellsInCuts", &eventRecord.nCellsInCuts, "nCellsInCuts/I");
    eventTree->Branch("sumWInCuts", &eventRecord.sumWInCuts, "sumWInCuts/D");
    // Ring
    eventTree->Branch("ringD3N", &eventRecord.observableD3N[kRing], "ringD3N/D");
    eventTree->Branch("ringFlat", &eventRecord.ringFlat, "ringFlat/D");
    eventTree->Branch("ringPerpD3N", &eventRecord.observableD3N[kRingPerp], "ringPerpD3N/D");
    // Ring Z
    eventTree->Branch("ringZD3N", &eventRecord.observableD3N[kRingZ], "ringZD3N/D");
    eventTree->Branch("ringZFlat", &eventRecord.ringZFlat, "ringZFlat/D");
    eventTree->Branch("ringZLabD3N", &eventRecord.observableD3N[kRingZLab], "ringZLabD3N/D");
    // Ring ZT and the helicity term
    eventTree->Branch("ringZTD3N", &eventRecord.observableD3N[kRingZT], "ringZTD3N/D");
    eventTree->Branch("helicityTermD3N", &eventRecord.observableD3N[kHelicityTerm], "helicityTermD3N/D");
    // Normalizations of A_z and of the expected precision
    eventTree->Branch("nzSqD3N", &eventRecord.observableD3N[kNzSq], "nzSqD3N/D");
    eventTree->Branch("nzSqSinSqThetaD3N", &eventRecord.observableD3N[kNzSqSinSqTheta], "nzSqSinSqThetaD3N/D");
    // P_e, for the leakage ansatz of the trigger variants
    eventTree->Branch("polED3N", &eventRecord.polED3N, "polED3N/D");
    // Per-bin sums A_b and B_b of the true trigger, for the jackknife of any pooled quantity downstream
    eventTree->Branch("nDeltaPhiBins", &eventRecord.nDeltaPhiBins, "nDeltaPhiBins/I");
    eventTree->Branch("sumWDeltaPhi", eventRecord.sumWDeltaPhi.data(), "sumWDeltaPhi[nDeltaPhiBins]/D");
    eventTree->Branch("sumARingDeltaPhi", eventRecord.sumARingDeltaPhi.data(), "sumARingDeltaPhi[nDeltaPhiBins]/D");
    eventTree->Branch("sumARingZDeltaPhi", eventRecord.sumARingZDeltaPhi.data(), "sumARingZDeltaPhi[nDeltaPhiBins]/D");
    eventTree->Branch("sumANzSqDeltaPhi", eventRecord.sumANzSqDeltaPhi.data(), "sumANzSqDeltaPhi[nDeltaPhiBins]/D");
    // Delta scan
    eventTree->Branch("nDelta", &eventRecord.nDelta, "nDelta/I");
    eventTree->Branch("scanRingD3N", eventRecord.scanRingD3N.data(), "scanRingD3N[nDelta]/D");
    eventTree->Branch("scanRingFlat", eventRecord.scanRingFlat.data(), "scanRingFlat[nDelta]/D");
    eventTree->Branch("scanRingZD3N", eventRecord.scanRingZD3N.data(), "scanRingZD3N[nDelta]/D");
    eventTree->Branch("scanRingZFlat", eventRecord.scanRingZFlat.data(), "scanRingZFlat[nDelta]/D");
    eventTree->Branch("scanNzSqD3N", eventRecord.scanNzSqD3N.data(), "scanNzSqD3N[nDelta]/D");
    eventTree->Branch("scanCosZetaD3N", eventRecord.scanCosZetaD3N.data(), "scanCosZetaD3N[nDelta]/D");
    eventTree->Branch("scanSinZetaD3N", eventRecord.scanSinZetaD3N.data(), "scanSinZetaD3N[nDelta]/D");
    // Mixing (pooled over the partners of each event)
    eventTree->Branch("mixRingD3N", &eventRecord.mixRingD3N, "mixRingD3N/D");
    eventTree->Branch("mixRingFlat", &eventRecord.mixRingFlat, "mixRingFlat/D");
    eventTree->Branch("mixRingZD3N", &eventRecord.mixRingZD3N, "mixRingZD3N/D");
    eventTree->Branch("mixRingZFlat", &eventRecord.mixRingZFlat, "mixRingZFlat/D");
    eventTree->Branch("nMixPartners", &eventRecord.nMixPartners, "nMixPartners/I");
    eventTree->Branch("mixDeltaEffBin", eventRecord.mixDeltaEffBin.data(), "mixDeltaEffBin[nMixPartners]/I");
    eventTree->Branch("nMixDeltaPhiBins", &eventRecord.nMixDeltaPhiBins, "nMixDeltaPhiBins/I");
    eventTree->Branch("mixSumWDeltaPhi", eventRecord.mixSumWDeltaPhi.data(), "mixSumWDeltaPhi[nMixDeltaPhiBins]/D");
    eventTree->Branch("mixSumARingDeltaPhi", eventRecord.mixSumARingDeltaPhi.data(), "mixSumARingDeltaPhi[nMixDeltaPhiBins]/D");
    eventTree->Branch("mixSumARingZDeltaPhi", eventRecord.mixSumARingZDeltaPhi.data(), "mixSumARingZDeltaPhi[nMixDeltaPhiBins]/D");
    eventTree->Branch("mixSumANzSqDeltaPhi", eventRecord.mixSumANzSqDeltaPhi.data(), "mixSumANzSqDeltaPhi[nMixDeltaPhiBins]/D");
    eventTree->Branch("mixSumACosZetaDeltaPhi", eventRecord.mixSumACosZetaDeltaPhi.data(), "mixSumACosZetaDeltaPhi[nMixDeltaPhiBins]/D");
    eventTree->Branch("mixSumASinZetaDeltaPhi", eventRecord.mixSumASinZetaDeltaPhi.data(), "mixSumASinZetaDeltaPhi[nMixDeltaPhiBins]/D");
    // Smudging (pooled over the draws of each event)
    eventTree->Branch("smudgeRingD3N", &eventRecord.smudgeRingD3N, "smudgeRingD3N/D");
    eventTree->Branch("smudgeRingFlat", &eventRecord.smudgeRingFlat, "smudgeRingFlat/D");
    eventTree->Branch("smudgeRingZD3N", &eventRecord.smudgeRingZD3N, "smudgeRingZD3N/D");
    eventTree->Branch("smudgeRingZFlat", &eventRecord.smudgeRingZFlat, "smudgeRingZFlat/D");
    eventTree->Branch("nSmudgeDeltaPhiBins", &eventRecord.nSmudgeDeltaPhiBins, "nSmudgeDeltaPhiBins/I");
    eventTree->Branch("smudgeSumWDeltaPhi", eventRecord.smudgeSumWDeltaPhi.data(), "smudgeSumWDeltaPhi[nSmudgeDeltaPhiBins]/D");
    eventTree->Branch("smudgeSumARingDeltaPhi", eventRecord.smudgeSumARingDeltaPhi.data(), "smudgeSumARingDeltaPhi[nSmudgeDeltaPhiBins]/D");
    eventTree->Branch("smudgeSumARingZDeltaPhi", eventRecord.smudgeSumARingZDeltaPhi.data(), "smudgeSumARingZDeltaPhi[nSmudgeDeltaPhiBins]/D");
    eventTree->Branch("smudgeSumANzSqDeltaPhi", eventRecord.smudgeSumANzSqDeltaPhi.data(), "smudgeSumANzSqDeltaPhi[nSmudgeDeltaPhiBins]/D");
    eventTree->Branch("smudgeSumACosZetaDeltaPhi", eventRecord.smudgeSumACosZetaDeltaPhi.data(), "smudgeSumACosZetaDeltaPhi[nSmudgeDeltaPhiBins]/D");
    eventTree->Branch("smudgeSumASinZetaDeltaPhi", eventRecord.smudgeSumASinZetaDeltaPhi.data(), "smudgeSumASinZetaDeltaPhi[nSmudgeDeltaPhiBins]/D");

    // ==================================================================================================================================
    // 7. Event loop
    // ==================================================================================================================================
    // Per-event fields, already converted to polarization (kColumnToPolarization)
    std::vector<double> cellWeight(nCells), polLabT(nCells), polLabX(nCells), polLabY(nCells), polLabZ(nCells);
    std::vector<double> polStarT(nCells), polStarX(nCells), polStarY(nCells), polStarZ(nCells);
    std::vector<char> cellValid(nCells);
    std::array<double, kNRingObservables> observable{};
    std::array<double, kNPolComponents> polComponent{};
    std::array<double, kNPolComponents> polJetComponent{};
    int eventsProcessed = 0;

    for (const int eventId : eventIds) {
      // ------------------------------------------------------------------------------------------------------------------------------
      // 7a. Read the event and check its grid (the first event is already in memory). The files are written from the same iSS grid arrays
      //     with the same format, so identical nodes give identical text, hence identical doubles: the exact comparison is intended
      // ------------------------------------------------------------------------------------------------------------------------------
      if (eventId != eventIds.front()) {
        readCellTable(eventInputPath(config, eventId), table);
        if (table.size() != grid.nCells)
          throw std::runtime_error("event " + std::to_string(eventId) + " has " + std::to_string(table.size()) + " cells instead of " + std::to_string(grid.nCells));
        for (int iCell = 0; iCell < grid.nCells; ++iCell)
          if (table.column[kColY][iCell] != grid.y[iCell] || table.column[kColPt][iCell] != grid.pt[iCell] || table.column[kColPhi][iCell] != grid.phi[iCell])
            throw std::runtime_error("event " + std::to_string(eventId) + ", data line " + std::to_string(iCell + 1) + ": (y, pT, phi) differs from the first event");
      }

      // ------------------------------------------------------------------------------------------------------------------------------
      // 7b. Cell weights and polarization of this event (columns converted with kColumnToPolarization), plus the compact in-cut copy
      // ------------------------------------------------------------------------------------------------------------------------------
      int nInvalidCells = 0;
      for (int iCell = 0; iCell < nCells; ++iCell) {
        const double multiplicity = table.column[kColMult][iCell];
        bool isValid = std::isfinite(multiplicity) && multiplicity >= 0.;
        for (int iColumn = kColPolT; iColumn <= kColPolStarZ; ++iColumn)
          isValid = isValid && std::isfinite(table.column[iColumn][iCell]);
        cellValid[iCell] = isValid;
        if (!isValid) {
          // Zeroed as well as flagged, so that the compact copy (weight 0) cannot carry a NaN into the variant sums
          ++nInvalidCells;
          cellWeight[iCell] = polLabT[iCell] = polLabX[iCell] = polLabY[iCell] = polLabZ[iCell] = 0.;
          polStarT[iCell] = polStarX[iCell] = polStarY[iCell] = polStarZ[iCell] = 0.;
          continue;
        }
        cellWeight[iCell] = multiplicity * grid.cellMeasure[iCell];
        polLabT[iCell] = kColumnToPolarization * table.column[kColPolT][iCell];
        polLabX[iCell] = kColumnToPolarization * table.column[kColPolX][iCell];
        polLabY[iCell] = kColumnToPolarization * table.column[kColPolY][iCell];
        polLabZ[iCell] = kColumnToPolarization * table.column[kColPolZ][iCell];
        polStarT[iCell] = kColumnToPolarization * table.column[kColPolStarT][iCell];
        polStarX[iCell] = kColumnToPolarization * table.column[kColPolStarX][iCell];
        polStarY[iCell] = kColumnToPolarization * table.column[kColPolStarY][iCell];
        polStarZ[iCell] = kColumnToPolarization * table.column[kColPolStarZ][iCell];
      }
      if (nInvalidCells > 0) {
        std::cout << "[WARNING] event " << eventId << ": " << nInvalidCells << " cells with non-finite or negative columns are skipped" << std::endl;
        warnings.invalidCells += nInvalidCells;
      }
      for (int iCut = 0; iCut < nInCutCells; ++iCut) {
        const int iCell = inCutCells.cellIndex[iCut];
        inCutCells.weight[iCut] = cellWeight[iCell];
        inCutCells.polStarX[iCut] = polStarX[iCell];
        inCutCells.polStarY[iCut] = polStarY[iCell];
        inCutCells.polStarZ[iCut] = polStarZ[iCell];
      }

      const bool hasJetEntry = eventId < static_cast<int>(jets.size()) && jets[eventId].present;
      const double phiTrigger = config.triggerMode == TriggerMode::kJet ? jets[eventId].phi : 0.;
      eventRecord.eventId = eventId;
      eventRecord.phiJ = hasJetEntry ? jets[eventId].phi : kNoValue;
      eventRecord.phiTrigger = phiTrigger;
      eventRecord.nCellsInCuts = nInCutCells;
      trueTriggerOutputs.phiTrigger->Fill(wrapAngle(phiTrigger, 0.));

      // ------------------------------------------------------------------------------------------------------------------------------
      // 7c. True trigger, in one pass over the cells: observables, maps, QA and closures (README, "Output reference").
      //     The true trigger is always transverse (t_z = 0). Regions: 1D curves use the pT and y cuts, 2D maps with a pT axis use the y cut
      //     only, the y curves use the pT cut only, and the momentum maps, NoCuts/, Inputs/ and the closures use every cell.
      // ------------------------------------------------------------------------------------------------------------------------------
      const double triggerX = std::cos(phiTrigger);
      const double triggerY = std::sin(phiTrigger);
      const int trueBinShift = deltaPhiBinShift(phiTrigger, grid);
      // Ring direction of the true trigger for the in-cut cells: the reference n_t of cos(zeta) and sin(zeta) in the variant loops (7e-7g).
      // A cell with p parallel to t (impossible without a y = 0 node) is counted by the loop below, and gets n_t = 0, hence zero overlap
      for (int iCut = 0; iCut < nInCutCells; ++iCut) {
        const double crossX = triggerY * inCutCells.pz[iCut];
        const double crossY = -triggerX * inCutCells.pz[iCut];
        const double crossZ = triggerX * inCutCells.py[iCut] - triggerY * inCutCells.px[iCut];
        const double crossNorm = std::sqrt(squared(crossX) + squared(crossY) + squared(crossZ));
        const double invCrossNorm = crossNorm > 0. ? 1. / crossNorm : 0.;
        inCutCells.nTrueX[iCut] = crossX * invCrossNorm;
        inCutCells.nTrueY[iCut] = crossY * invCrossNorm;
        inCutCells.nTrueZ[iCut] = crossZ * invCrossNorm;
      }
      for (int iCell = 0; iCell < nCells; ++iCell) {
        if (!cellValid[iCell]) // already counted and reported above
          continue;
        const double weight = cellWeight[iCell];
        const double px = grid.px[iCell];
        const double py = grid.py[iCell];
        const double pz = grid.pz[iCell];
        const double pHatX = grid.pHatX[iCell];
        const double pHatY = grid.pHatY[iCell];
        const double pHatZ = grid.pHatZ[iCell];
        const int iPt = grid.iPt[iCell];
        const int iY = grid.iY[iCell];
        int deltaPhiBin = grid.iPhi[iCell] + trueBinShift;
        if (deltaPhiBin >= nPhiBins)
          deltaPhiBin -= nPhiBins;
        int phiBin = grid.iPhi[iCell] + reactionPlaneBinShift; // phi - Psi_RP, wrapped like Delta phi
        if (phiBin >= nPhiBins)
          phiBin -= nPhiBins;

        // Closures, one count per cell, as relative residuals so that they read directly against the print precision (README, "Closure checks")
        const double polLabMag = std::sqrt(squared(polLabX[iCell]) + squared(polLabY[iCell]) + squared(polLabZ[iCell]));
        const double polStarMag = std::sqrt(squared(polStarX[iCell]) + squared(polStarY[iCell]) + squared(polStarZ[iCell]));
        const double minkowskiPolDotP = polLabT[iCell] * grid.energy[iCell] - (polLabX[iCell] * px + polLabY[iCell] * py + polLabZ[iCell] * pz);
        trueTriggerOutputs.closurePdotP.add(logResidualBin(relativeResidual(minkowskiPolDotP, grid.energy[iCell] * polLabMag)), 1.);
        trueTriggerOutputs.closurePStarTime.add(logResidualBin(relativeResidual(polStarT[iCell], polStarMag)), 1.);
        const PxPyPzEVector polLab4Vec(polLabX[iCell], polLabY[iCell], polLabZ[iCell], polLabT[iCell]);
        const PxPyPzEVector polStarOurs4Vec = ROOT::Math::VectorUtil::boost(polLab4Vec, grid.boostToRest[iCell]);
        const double boostResidual = std::sqrt(squared(polStarOurs4Vec.X() - polStarX[iCell]) + squared(polStarOurs4Vec.Y() - polStarY[iCell]) + squared(polStarOurs4Vec.Z() - polStarZ[iCell]));
        trueTriggerOutputs.closureBoost.add(logResidualBin(relativeResidual(boostResidual, polStarMag)), 1.);

        // Inputs/
        trueTriggerOutputs.yieldVsPhi.add(phiBin + 1, weight);
        trueTriggerOutputs.yieldVsPt.add(iPt + 1, weight);
        trueTriggerOutputs.yieldVsY.add(iY + 1, weight);

        // Ring direction n = (t x p)/|t x p| in lab components, with t transverse
        const double crossX = triggerY * pz;
        const double crossY = -triggerX * pz;
        const double crossZ = triggerX * py - triggerY * px;
        const double crossNormSq = squared(crossX) + squared(crossY) + squared(crossZ);
        if (!(crossNormSq > 0.)) { // p parallel to t: n undefined. Impossible on a grid without a y = 0 node, but never assumed
          ++warnings.degenerateRingDirection;
          continue;
        }
        const double invCrossNorm = 1. / std::sqrt(crossNormSq);
        const double nX = crossX * invCrossNorm;
        const double nY = crossY * invCrossNorm;
        const double nZ = crossZ * invCrossNorm;

        // P* from the LRF block. The lab projection P.n is identical because n is orthogonal to p: its difference is a closure only
        const double ringStar = polStarX[iCell] * nX + polStarY[iCell] * nY + polStarZ[iCell] * nZ;
        const double ringLab = polLabX[iCell] * nX + polLabY[iCell] * nY + polLabZ[iCell] * nZ;
        trueTriggerOutputs.closureRingLabMinusStar.add(logResidualBin(relativeResidual(ringLab - ringStar, polStarMag)), 1.);
        const double polStarAlongP = polStarX[iCell] * pHatX + polStarY[iCell] * pHatY + polStarZ[iCell] * pHatZ; // P*_parallel
        observable[kRing] = ringStar;
        observable[kRingZ] = polStarZ[iCell] * nZ;
        observable[kHelicityTerm] = polStarAlongP * pHatZ * nZ;
        observable[kRingZT] = observable[kRingZ] - observable[kHelicityTerm]; // (P*_z - P*_parallel p_hat_z) n_z
        observable[kRingZLab] = polLabZ[iCell] * nZ;
        observable[kRingPerp] = ringStar - observable[kRingZ];
        observable[kNzSq] = nZ * nZ;
        observable[kNzSqSinSqTheta] = nZ * nZ * grid.sinSqTheta[iCell];

        // Trigger-aligned coordinate system of PolVsAngles/: x' along t, y' = z x x', z' = z
        polComponent = {polLabX[iCell], polLabY[iCell], polLabZ[iCell], polStarX[iCell], polStarY[iCell], polStarZ[iCell]};
        polJetComponent = {polLabX[iCell] * triggerX + polLabY[iCell] * triggerY, -polLabX[iCell] * triggerY + polLabY[iCell] * triggerX, polLabZ[iCell],
                           polStarX[iCell] * triggerX + polStarY[iCell] * triggerY, -polStarX[iCell] * triggerY + polStarY[iCell] * triggerX, polStarZ[iCell]};

        // Full grid
        trueTriggerOutputs.ringVsDeltaPhiNoCuts.add(deltaPhiBin, weight, ringStar);
        trueTriggerOutputs.ringZVsDeltaPhiNoCuts.add(deltaPhiBin, weight, observable[kRingZ]);

        // PolarizationMaps/: P* in the three coordinate systems of the consumer, with the field of the node spread over the sub-cell points
        // (README, "Polarization maps"). The true trigger is transverse (t_T = 1), so the PrimeJet axes are x' = z, y' = t x z and z' = t
        const std::array<double, kNMapComponents> polStarLab = {polStarX[iCell], polStarY[iCell], polStarZ[iCell]};
        const std::array<double, kNMapComponents> polStarPrimeJet = {polStarZ[iCell], polStarX[iCell] * triggerY - polStarY[iCell] * triggerX, polStarX[iCell] * triggerX + polStarY[iCell] * triggerY};
        const double invPt = 1. / grid.pt[iCell];
        const std::array<double, kNMapComponents> polStarPrimeV0 = {(polStarX[iCell] * px + polStarY[iCell] * py) * invPt, (polStarY[iCell] * px - polStarX[iCell] * py) * invPt, polStarZ[iCell]};
        const std::array<double, kNMapRings> ringScalars = {ringStar, observable[kRingZ]};
        const int firstPoint = iCell * subCellPoints.perCell;
        for (int iPoint = firstPoint; iPoint < firstPoint + subCellPoints.perCell; ++iPoint) {
          const double pointWeight = weight * subCellPoints.weightFraction[iPoint];
          const int pxBin = subCellPoints.pxBin[iPoint];
          const int pyBin = subCellPoints.pyBin[iPoint];
          const int pzBin = subCellPoints.pzBin[iPoint]; // p_z is also p_x'Jet
          const int pyPrimeJetBin = transverseMomentumAxis.findBin(subCellPoints.px[iPoint] * triggerY - subCellPoints.py[iPoint] * triggerX);
          const int pzPrimeJetBin = transverseMomentumAxis.findBin(subCellPoints.px[iPoint] * triggerX + subCellPoints.py[iPoint] * triggerY);
          if (pyPrimeJetBin < 0 || pzPrimeJetBin < 0) { // the transverse axis holds |p_T| of every point, so this should never happen
            ++warnings.pointsOutsideMaps;
            continue;
          }
          for (int iComp = 0; iComp < kNMapComponents; ++iComp) {
            trueTriggerOutputs.labMapPxPy[iComp].add(pxBin, pyBin, pointWeight, polStarLab[iComp]);
            trueTriggerOutputs.labMapPzPx[iComp].add(pzBin, pxBin, pointWeight, polStarLab[iComp]);
            trueTriggerOutputs.labMapPyPz[iComp].add(pyBin, pzBin, pointWeight, polStarLab[iComp]);
            trueTriggerOutputs.primeJetMapPxPy[iComp].add(pzBin, pyPrimeJetBin, pointWeight, polStarPrimeJet[iComp]);
            trueTriggerOutputs.primeJetMapPzPx[iComp].add(pzPrimeJetBin, pzBin, pointWeight, polStarPrimeJet[iComp]);
            trueTriggerOutputs.primeJetMapPyPz[iComp].add(pyPrimeJetBin, pzPrimeJetBin, pointWeight, polStarPrimeJet[iComp]);
            trueTriggerOutputs.primeV0MapPzPt[iComp].add(pzBin, iPt, pointWeight, polStarPrimeV0[iComp]);
          }
          for (int iRing = 0; iRing < kNMapRings; ++iRing) {
            trueTriggerOutputs.primeJetRingPxPy[iRing].add(pzBin, pyPrimeJetBin, pointWeight, ringScalars[iRing]);
            trueTriggerOutputs.primeJetRingPzPx[iRing].add(pzPrimeJetBin, pzBin, pointWeight, ringScalars[iRing]);
            trueTriggerOutputs.primeJetRingPyPz[iRing].add(pyPrimeJetBin, pzPrimeJetBin, pointWeight, ringScalars[iRing]);
          }
          trueTriggerOutputs.labYieldPxPy.add(pxBin, pyBin, pointWeight);
          trueTriggerOutputs.labYieldPzPx.add(pzBin, pxBin, pointWeight);
          trueTriggerOutputs.labYieldPyPz.add(pyBin, pzBin, pointWeight);
          trueTriggerOutputs.primeJetYieldPxPy.add(pzBin, pyPrimeJetBin, pointWeight);
          trueTriggerOutputs.primeJetYieldPzPx.add(pzPrimeJetBin, pzBin, pointWeight);
          trueTriggerOutputs.primeJetYieldPyPz.add(pyPrimeJetBin, pzPrimeJetBin, pointWeight);
          trueTriggerOutputs.primeV0YieldPzPt.add(pzBin, iPt, pointWeight);
        }

        const bool inPtCut = grid.inPtCut[iCell];
        const bool inYCut = grid.inYCut[iCell];

        // y cut only (all pT)
        if (inYCut) {
          trueTriggerOutputs.yieldVsDeltaPhiVsPt.add(deltaPhiBin, iPt, weight);
          trueTriggerOutputs.ringVsDeltaPhiVsPt.add(deltaPhiBin, iPt, weight, ringStar);
          trueTriggerOutputs.ringZVsDeltaPhiVsPt.add(deltaPhiBin, iPt, weight, observable[kRingZ]);
          trueTriggerOutputs.ringVsPt.add(iPt, weight, ringStar);
          trueTriggerOutputs.ringZVsPt.add(iPt, weight, observable[kRingZ]);
          for (int iComp = 0; iComp < kNPolComponents; ++iComp) {
            trueTriggerOutputs.polVsPhiVsPt[iComp].add(phiBin, iPt, weight, polComponent[iComp]);
            trueTriggerOutputs.polJetVsDeltaPhiVsPt[iComp].add(deltaPhiBin, iPt, weight, polJetComponent[iComp]);
          }
        }

        // pT cut only (all y)
        if (inPtCut) {
          trueTriggerOutputs.ringVsDeltaPhiVsY.add(deltaPhiBin, iY, weight, ringStar);
          trueTriggerOutputs.ringZVsDeltaPhiVsY.add(deltaPhiBin, iY, weight, observable[kRingZ]);
          trueTriggerOutputs.ringVsY.add(iY, weight, ringStar);
          trueTriggerOutputs.ringZVsY.add(iY, weight, observable[kRingZ]);
        }

        if (!(inPtCut && inYCut))
          continue;

        // Both cuts
        for (int iObs = 0; iObs < kNRingObservables; ++iObs)
          trueTriggerOutputs.obsVsDeltaPhi[iObs].add(deltaPhiBin, weight, observable[iObs]);
        for (int iComp = 0; iComp < kNPolComponents; ++iComp) {
          trueTriggerOutputs.polVsPhi[iComp].add(phiBin, weight, polComponent[iComp]);
          trueTriggerOutputs.polJetVsDeltaPhi[iComp].add(deltaPhiBin, weight, polJetComponent[iComp]);
        }

        // Ring-basis decomposition of P*: along n (that is Ring itself), along e = p_hat x n, and along p_hat
        const double eX = pHatY * nZ - pHatZ * nY;
        const double eY = pHatZ * nX - pHatX * nZ;
        const double eZ = pHatX * nY - pHatY * nX;
        trueTriggerOutputs.polEVsDeltaPhi.add(deltaPhiBin, weight, polStarX[iCell] * eX + polStarY[iCell] * eY + polStarZ[iCell] * eZ);
        trueTriggerOutputs.polParVsDeltaPhi.add(deltaPhiBin, weight, polStarAlongP);
        if (polStarMag > 0.)
          trueTriggerOutputs.alignmentVsDeltaPhi.add(deltaPhiBin, weight, ringStar / polStarMag);
        else
          ++warnings.zeroPolStarCells;

        // chi exactly as in computeProxyRing of the consumer (n = cos(chi) phi_hat + sin(chi) theta_hat), from the unnormalized cross product
        const double chi = std::atan2(-crossZ * grid.pMag[iCell], crossY * px - crossX * py);
        // atan2 returns [-pi, pi], and chi within an ulp of +pi does occur: when phi_J and a phi node are the same printed number (the regular
        // jet grids of break_symmetry and lumpy_events_0_5_fine), Delta phi = 0 exactly and only rounding is left in crossZ
        const int chiBin = chiAxis.periodicBin(chi);
        trueTriggerOutputs.ringVsChi.add(chiBin, weight, ringStar);
        trueTriggerOutputs.ringZVsChi.add(chiBin, weight, observable[kRingZ]);
        trueTriggerOutputs.chiYield.add(chiBin + 1, weight);
      }

      // ------------------------------------------------------------------------------------------------------------------------------
      // 7d. Per-event integrals of the true trigger, the tree buffers, and the flush of every true-trigger accumulator
      // ------------------------------------------------------------------------------------------------------------------------------
      const PerEventProfile1D& ringProfile = trueTriggerOutputs.obsVsDeltaPhi[kRing];
      const double sumWInCuts = sumOf(ringProfile.sumW); // same cells, hence same sum, for every observable
      if (!(sumWInCuts > 0.))                            // every integral of this event would be 0/0
        throw std::runtime_error("event " + std::to_string(eventId) + " has no valid cell (or zero yield) inside the cuts");
      eventRecord.sumWInCuts = sumWInCuts;
      for (int iObs = 0; iObs < kNRingObservables; ++iObs) {
        eventRecord.observableD3N[iObs] = sumOf(trueTriggerOutputs.obsVsDeltaPhi[iObs].sumWX) / sumWInCuts;
        fillEventPoint(trueTriggerOutputs.obsD3N[iObs], trueTriggerOutputs.obsD3NEvtMean[iObs], 0.5, eventRecord.observableD3N[iObs], sumWInCuts);
      }
      warnings.emptyDeltaPhiBins += countEmptyBins(ringProfile.sumW);
      eventRecord.ringFlat = flatAverage(ringProfile.sumW, ringProfile.sumWX);
      eventRecord.ringZFlat = flatAverage(ringProfile.sumW, trueTriggerOutputs.obsVsDeltaPhi[kRingZ].sumWX);
      fillEventPoint(trueTriggerOutputs.ringFlat, trueTriggerOutputs.ringFlatEvtMean, 0.5, eventRecord.ringFlat, sumWInCuts);
      fillEventPoint(trueTriggerOutputs.ringZFlat, trueTriggerOutputs.ringZFlatEvtMean, 0.5, eventRecord.ringZFlat, sumWInCuts);
      // P_e is filled for the same cells as the observables (both cuts, n defined), hence the same weight
      eventRecord.polED3N = sumOf(trueTriggerOutputs.polEVsDeltaPhi.sumWX) / sumWInCuts;

      std::copy(ringProfile.sumW.begin(), ringProfile.sumW.end(), eventRecord.sumWDeltaPhi.begin());
      std::copy(ringProfile.sumWX.begin(), ringProfile.sumWX.end(), eventRecord.sumARingDeltaPhi.begin());
      std::copy(trueTriggerOutputs.obsVsDeltaPhi[kRingZ].sumWX.begin(), trueTriggerOutputs.obsVsDeltaPhi[kRingZ].sumWX.end(), eventRecord.sumARingZDeltaPhi.begin());
      std::copy(trueTriggerOutputs.obsVsDeltaPhi[kNzSq].sumWX.begin(), trueTriggerOutputs.obsVsDeltaPhi[kNzSq].sumWX.end(), eventRecord.sumANzSqDeltaPhi.begin());

      trueTriggerOutputs.ringD3NPerEvent.push_back(eventRecord.observableD3N[kRing]);
      trueTriggerOutputs.ringFlatPerEvent.push_back(eventRecord.ringFlat);
      trueTriggerOutputs.ringZD3NPerEvent.push_back(eventRecord.observableD3N[kRingZ]);
      trueTriggerOutputs.ringZFlatPerEvent.push_back(eventRecord.ringZFlat);
      trueTriggerOutputs.eventYieldPerEvent.push_back(sumWInCuts);

      trueTriggerOutputs.flushAndReset();

      // ------------------------------------------------------------------------------------------------------------------------------
      // 7e. Delta-scan: deterministic event mixing, with the true trigger rotated by delta_k = 2 pi k / nSteps (README, "Delta-scan, mixing
      //     and smudging")
      // ------------------------------------------------------------------------------------------------------------------------------
      if (config.deltaScanSteps > 0) {
        const int nSteps = config.deltaScanSteps;
        for (int iStep = 0; iStep < nSteps; ++iStep) {
          const double delta = deltaAxis.binCentre(iStep);
          const double rotatedPhi = phiTrigger + delta;
          const IntegratedSums sums = accumulateReducedSet(XYZVector(std::cos(rotatedPhi), std::sin(rotatedPhi), 0.), inCutCells, grid, scanAccumulators, warnings);
          warnings.emptyDeltaPhiBins += countEmptyBins(scanAccumulators.sumW);
          eventRecord.scanRingD3N[iStep] = sums.sumWRing / sums.sumW;
          eventRecord.scanRingZD3N[iStep] = sums.sumWRingZ / sums.sumW;
          eventRecord.scanNzSqD3N[iStep] = sums.sumWNzSq / sums.sumW;
          eventRecord.scanRingFlat[iStep] = flatAverage(scanAccumulators.sumW, scanAccumulators.sumWRing);
          eventRecord.scanRingZFlat[iStep] = flatAverage(scanAccumulators.sumW, scanAccumulators.sumWRingZ);

          fillEventPoint(deltaScanOutputs.ringD3NVsDelta, deltaScanOutputs.ringD3NVsDeltaEvtMean, delta, eventRecord.scanRingD3N[iStep], sums.sumW);
          fillEventPoint(deltaScanOutputs.ringFlatVsDelta, deltaScanOutputs.ringFlatVsDeltaEvtMean, delta, eventRecord.scanRingFlat[iStep], sums.sumW);
          fillEventPoint(deltaScanOutputs.ringZD3NVsDelta, deltaScanOutputs.ringZD3NVsDeltaEvtMean, delta, eventRecord.scanRingZD3N[iStep], sums.sumW);
          fillEventPoint(deltaScanOutputs.ringZFlatVsDelta, deltaScanOutputs.ringZFlatVsDeltaEvtMean, delta, eventRecord.scanRingZFlat[iStep], sums.sumW);
          fillEventPoint(deltaScanOutputs.nzSqD3NVsDelta, deltaScanOutputs.nzSqD3NVsDeltaEvtMean, delta, eventRecord.scanNzSqD3N[iStep], sums.sumW);
          eventRecord.scanCosZetaD3N[iStep] = sums.sumWCosZeta / sums.sumW;
          eventRecord.scanSinZetaD3N[iStep] = sums.sumWSinZeta / sums.sumW;
          fillEventPoint(deltaScanOutputs.cosZetaD3NVsDelta, nullptr, delta, eventRecord.scanCosZetaD3N[iStep], sums.sumW);
          fillEventPoint(deltaScanOutputs.sinZetaD3NVsDelta, nullptr, delta, eventRecord.scanSinZetaD3N[iStep], sums.sumW);
          for (int bin = 0; bin < kZetaBins; ++bin) {
            deltaScanOutputs.cosZetaVsDelta.add(iStep, bin, scanAccumulators.cosZetaDistribution[bin]);
            deltaScanOutputs.sinZetaVsDelta.add(iStep, bin, scanAccumulators.sinZetaDistribution[bin]);
          }
          for (int bin = 0; bin < nPhiBins; ++bin) {
            if (scanAccumulators.sumW[bin] <= 0.) // already counted as an empty bin
              continue;
            const double deltaPhiPrime = deltaPhiAxis.binCentre(bin);
            deltaScanOutputs.ringVsDeltaPhiVsDelta->Fill(deltaPhiPrime, delta, scanAccumulators.sumWRing[bin] / scanAccumulators.sumW[bin], scanAccumulators.sumW[bin]);
            deltaScanOutputs.ringZVsDeltaPhiVsDelta->Fill(deltaPhiPrime, delta, scanAccumulators.sumWRingZ[bin] / scanAccumulators.sumW[bin], scanAccumulators.sumW[bin]);
          }
          scanAccumulators.reset();
        }

        // Half-turn closure: for t_z = 0, n(-t) = -n(t), so the integral at delta + pi is exactly minus the one at delta in every event
        double maxResidual = 0.;
        double maxAbsIntegral = 0.;
        for (int iStep = 0; iStep < nSteps / 2; ++iStep) {
          maxResidual = std::max(maxResidual, std::abs(eventRecord.scanRingD3N[iStep] + eventRecord.scanRingD3N[iStep + nSteps / 2]));
          maxAbsIntegral = std::max({maxAbsIntegral, std::abs(eventRecord.scanRingD3N[iStep]), std::abs(eventRecord.scanRingD3N[iStep + nSteps / 2])});
        }
        const double halfTurnResidual = relativeResidual(maxResidual, maxAbsIntegral);
        // The closure axis is in log10, so exact zeros are drawn at its lowest edge instead of the underflow
        deltaScanOutputs.closureHalfTurn->Fill(halfTurnResidual > 0. ? std::log10(halfTurnResidual) : kLogResidualLow);
      }

      // ------------------------------------------------------------------------------------------------------------------------------
      // 7f. Uniform event mixing: this event takes the phi_J of the next mixDepth events (mod nEvents) as trigger. Partners are chosen by
      //     event index, so the bullet and no_bullet events of one pair always mix with the same jets (README, key hypothesis)
      // ------------------------------------------------------------------------------------------------------------------------------
      if (config.mixDepth > 0) {
        for (int iPartner = 1; iPartner <= config.mixDepth; ++iPartner) {
          const int partnerId = config.firstEvent + (eventId - config.firstEvent + iPartner) % config.nEvents;
          const double partnerPhi = jets[partnerId].phi;
          const IntegratedSums sums = accumulateReducedSet(XYZVector(std::cos(partnerPhi), std::sin(partnerPhi), 0.), inCutCells, grid, mixingAccumulators, warnings);
          const double deltaEff = wrapAngle(partnerPhi - phiTrigger, deltaAxis.low);
          // The bin of the profile itself, so that the pulls downstream group the partners exactly as the profile does (-1: the rounding put
          // delta_eff beyond the axis, into the profile's overflow)
          const int deltaEffBin = mixingExtras.ringD3NVsDeltaEff->GetXaxis()->FindFixBin(deltaEff) - 1;
          eventRecord.mixDeltaEffBin[iPartner - 1] = deltaEffBin < deltaAxis.nBins ? deltaEffBin : -1;
          fillEventPoint(mixingExtras.ringD3NVsDeltaEff, mixingExtras.ringD3NVsDeltaEffEvtMean, deltaEff, sums.sumWRing / sums.sumW, sums.sumW);
          fillEventPoint(mixingExtras.ringZD3NVsDeltaEff, mixingExtras.ringZD3NVsDeltaEffEvtMean, deltaEff, sums.sumWRingZ / sums.sumW, sums.sumW);
        }
        // Per-bin sums for the tree (the jackknife and the pairing of the mixed curves downstream), before the flush resets them
        std::copy(mixingAccumulators.sumW.begin(), mixingAccumulators.sumW.end(), eventRecord.mixSumWDeltaPhi.begin());
        std::copy(mixingAccumulators.sumWRing.begin(), mixingAccumulators.sumWRing.end(), eventRecord.mixSumARingDeltaPhi.begin());
        std::copy(mixingAccumulators.sumWRingZ.begin(), mixingAccumulators.sumWRingZ.end(), eventRecord.mixSumARingZDeltaPhi.begin());
        std::copy(mixingAccumulators.sumWNzSq.begin(), mixingAccumulators.sumWNzSq.end(), eventRecord.mixSumANzSqDeltaPhi.begin());
        std::copy(mixingAccumulators.sumWCosZeta.begin(), mixingAccumulators.sumWCosZeta.end(), eventRecord.mixSumACosZetaDeltaPhi.begin());
        std::copy(mixingAccumulators.sumWSinZeta.begin(), mixingAccumulators.sumWSinZeta.end(), eventRecord.mixSumASinZetaDeltaPhi.begin());
        const ReducedValues mixingValues = flushReducedSet(mixingAccumulators, mixingOutputs, warnings); // pooled over the partners
        eventRecord.mixRingD3N = mixingValues.ringD3N;
        eventRecord.mixRingFlat = mixingValues.ringFlat;
        eventRecord.mixRingZD3N = mixingValues.ringZD3N;
        eventRecord.mixRingZFlat = mixingValues.ringZFlat;
      }

      // ------------------------------------------------------------------------------------------------------------------------------
      // 7g. Smudged triggers around the true one. One generator per event, seeded from (seed, eventId) rather than one running stream: the
      //     bullet and no_bullet events of one pair draw identical smudged triggers (README, key hypothesis)
      // ------------------------------------------------------------------------------------------------------------------------------
      if (config.smudgeSigma > 0.) {
        std::seed_seq seedSequence{config.seed, static_cast<std::uint32_t>(eventId)};
        std::mt19937 generator(seedSequence);
        std::normal_distribution<double> smudgeGaus(0., config.smudgeSigma);
        std::uniform_real_distribution<double> azimuthUniform(0., kTwoPi);
        const XYZVector trueTrigger(triggerX, triggerY, 0.);
        for (int iDraw = 0; iDraw < config.smudgeDraws; ++iDraw) {
          XYZVector smudgedTrigger;
          double smudgeAngle = 0.; // angle between the smudged and the true trigger
          if (config.smudgeMode == SmudgeMode::kPhi) {
            // phi: a Gaussian shift of the azimuth only (t_z stays 0)
            const double phiShift = smudgeGaus(generator);
            smudgeAngle = std::abs(phiShift);
            smudgedTrigger = XYZVector(std::cos(phiTrigger + phiShift), std::sin(phiTrigger + phiShift), 0.);
          } else {
            // cone: as forceJetDirectionSmudge in the consumer, a rotation by a half-normal angle around a uniformly random axis perpendicular
            // to t (isotropic around t, unlike smearing eta and phi separately).
            // 1) A uniformly random axis perpendicular to t (same construction as the consumer)
            const XYZVector perpVec = std::abs(trueTrigger.X()) > 0.99 ? XYZVector(-trueTrigger.Z(), 0., trueTrigger.X()).Unit() : XYZVector(0., trueTrigger.Z(), -trueTrigger.Y()).Unit();
            const double smearAzimuth = azimuthUniform(generator);
            const XYZVector smearAxis = perpVec * std::cos(smearAzimuth) + trueTrigger.Cross(perpVec) * std::sin(smearAzimuth);
            // 2) Half-normal rotation angle, 3) Rodrigues rotation (the (k.v)(1 - cos) term vanishes, since the axis is perpendicular to t)
            smudgeAngle = std::abs(smudgeGaus(generator));
            smudgedTrigger = trueTrigger * std::cos(smudgeAngle) + smearAxis.Cross(trueTrigger) * std::sin(smudgeAngle);
          }
          const IntegratedSums sums = accumulateReducedSet(smudgedTrigger, inCutCells, grid, smudgeAccumulators, warnings);
          smudgeExtras.smudgeAngle->Fill(smudgeAngle);
          smudgeExtras.triggerEta->Fill(smudgedTrigger.Eta());
          fillEventPoint(smudgeExtras.ringD3NVsSmudgeAngle, smudgeExtras.ringD3NVsSmudgeAngleEvtMean, smudgeAngle, sums.sumWRing / sums.sumW, sums.sumW);
          fillEventPoint(smudgeExtras.ringZD3NVsSmudgeAngle, smudgeExtras.ringZD3NVsSmudgeAngleEvtMean, smudgeAngle, sums.sumWRingZ / sums.sumW, sums.sumW);
        }
        // Per-bin sums for the tree (the jackknife and the pairing of the smudged curves downstream), before the flush resets them
        std::copy(smudgeAccumulators.sumW.begin(), smudgeAccumulators.sumW.end(), eventRecord.smudgeSumWDeltaPhi.begin());
        std::copy(smudgeAccumulators.sumWRing.begin(), smudgeAccumulators.sumWRing.end(), eventRecord.smudgeSumARingDeltaPhi.begin());
        std::copy(smudgeAccumulators.sumWRingZ.begin(), smudgeAccumulators.sumWRingZ.end(), eventRecord.smudgeSumARingZDeltaPhi.begin());
        std::copy(smudgeAccumulators.sumWNzSq.begin(), smudgeAccumulators.sumWNzSq.end(), eventRecord.smudgeSumANzSqDeltaPhi.begin());
        std::copy(smudgeAccumulators.sumWCosZeta.begin(), smudgeAccumulators.sumWCosZeta.end(), eventRecord.smudgeSumACosZetaDeltaPhi.begin());
        std::copy(smudgeAccumulators.sumWSinZeta.begin(), smudgeAccumulators.sumWSinZeta.end(), eventRecord.smudgeSumASinZetaDeltaPhi.begin());
        const ReducedValues smudgeValues = flushReducedSet(smudgeAccumulators, smudgeOutputs, warnings); // pooled over the draws
        eventRecord.smudgeRingD3N = smudgeValues.ringD3N;
        eventRecord.smudgeRingFlat = smudgeValues.ringFlat;
        eventRecord.smudgeRingZD3N = smudgeValues.ringZD3N;
        eventRecord.smudgeRingZFlat = smudgeValues.ringZFlat;
      }

      // ------------------------------------------------------------------------------------------------------------------------------
      // 7h. Tree and progress
      // ------------------------------------------------------------------------------------------------------------------------------
      eventTree->Fill();
      ++eventsProcessed;
      if (eventsProcessed % kProgressEvery == 0 || eventId == eventIds.back())
        std::cout << "Processed " << eventsProcessed << "/" << eventIds.size() << " events (last: " << eventId << ")" << std::endl;
    }

    // ==================================================================================================================================
    // 8. Finalize: buffered histograms, per-event distributions, Paper/ copies, warnings, Config/, write
    // ==================================================================================================================================
    trueTriggerOutputs.yieldVsPhi.flush();
    trueTriggerOutputs.yieldVsPt.flush();
    trueTriggerOutputs.yieldVsY.flush();
    trueTriggerOutputs.yieldVsDeltaPhiVsPt.flush();
    trueTriggerOutputs.chiYield.flush();
    trueTriggerOutputs.labYieldPxPy.flush();
    trueTriggerOutputs.labYieldPzPx.flush();
    trueTriggerOutputs.labYieldPyPz.flush();
    trueTriggerOutputs.primeJetYieldPxPy.flush();
    trueTriggerOutputs.primeJetYieldPzPx.flush();
    trueTriggerOutputs.primeJetYieldPyPz.flush();
    trueTriggerOutputs.primeV0YieldPzPt.flush();
    if (config.deltaScanSteps > 0) {
      deltaScanOutputs.cosZetaVsDelta.flush();
      deltaScanOutputs.sinZetaVsDelta.flush();
    }
    if (config.mixDepth > 0) {
      mixingOutputs.cosZetaDistribution.flush();
      mixingOutputs.sinZetaDistribution.flush();
    }
    if (config.smudgeSigma > 0.) {
      smudgeOutputs.cosZetaDistribution.flush();
      smudgeOutputs.sinZetaDistribution.flush();
    }
    trueTriggerOutputs.closurePdotP.flush();
    trueTriggerOutputs.closurePStarTime.flush();
    trueTriggerOutputs.closureBoost.flush();
    trueTriggerOutputs.closureRingLabMinusStar.flush();

    // Per-event distributions (README: Paper/ holds the Ring ones, as the band of Fig. 11)
    bookDistribution(paperDir, "hRingD3NPerEvent", "Per-event integrated #it{R} (d^{3}N-weighted);#it{R};Events", trueTriggerOutputs.ringD3NPerEvent);
    bookDistribution(paperDir, "hRingFlatPerEvent", "Per-event integrated #it{R} (flat in #Delta#varphi);#it{R}_{flat};Events", trueTriggerOutputs.ringFlatPerEvent);
    bookDistribution(ringDir, "hRingZD3NPerEvent", "Per-event integrated #it{R}_{z} (d^{3}N-weighted);#it{R}_{z};Events", trueTriggerOutputs.ringZD3NPerEvent);
    bookDistribution(ringDir, "hRingZFlatPerEvent", "Per-event integrated #it{R}_{z} (flat in #Delta#varphi);#it{R}_{z,flat};Events", trueTriggerOutputs.ringZFlatPerEvent);
    bookDistribution(inputsDir, "hEventYield", "Per-event yield inside the cuts;#sum w;Events", trueTriggerOutputs.eventYieldPerEvent);
    if (config.mixDepth > 0) {
      bookDistribution(mixingDir, "hRingD3NPerEvent", "Per-event integrated #it{R} (d^{3}N-weighted, mixed);#it{R};Events", mixingOutputs.ringD3NPerEvent);
      bookDistribution(mixingDir, "hRingFlatPerEvent", "Per-event integrated #it{R} (flat, mixed);#it{R}_{flat};Events", mixingOutputs.ringFlatPerEvent);
      bookDistribution(mixingDir, "hRingZD3NPerEvent", "Per-event integrated #it{R}_{z} (d^{3}N-weighted, mixed);#it{R}_{z};Events", mixingOutputs.ringZD3NPerEvent);
      bookDistribution(mixingDir, "hRingZFlatPerEvent", "Per-event integrated #it{R}_{z} (flat, mixed);#it{R}_{z,flat};Events", mixingOutputs.ringZFlatPerEvent);
    }
    if (config.smudgeSigma > 0.) {
      bookDistribution(smudgeDir, "hRingD3NPerEvent", "Per-event integrated #it{R} (d^{3}N-weighted, smudged);#it{R};Events", smudgeOutputs.ringD3NPerEvent);
      bookDistribution(smudgeDir, "hRingFlatPerEvent", "Per-event integrated #it{R} (flat, smudged);#it{R}_{flat};Events", smudgeOutputs.ringFlatPerEvent);
      bookDistribution(smudgeDir, "hRingZD3NPerEvent", "Per-event integrated #it{R}_{z} (d^{3}N-weighted, smudged);#it{R}_{z};Events", smudgeOutputs.ringZD3NPerEvent);
      bookDistribution(smudgeDir, "hRingZFlatPerEvent", "Per-event integrated #it{R}_{z} (flat, smudged);#it{R}_{z,flat};Events", smudgeOutputs.ringZFlatPerEvent);
    }

    // Paper/: copies of the objects that make the figures (README, "Output reference"), filled once in their own folders
    cloneInto(paperDir, trueTriggerOutputs.polJetVsDeltaPhiVsPt[kPolLabZ].pooledMaster, "p2dPzLabVsDeltaPhiVsPt");
    cloneInto(paperDir, trueTriggerOutputs.polVsPhiVsPt[kPolLabZ].pooledMaster, "p2dPzLabVsPhiVsPt");
    cloneInto(paperDir, trueTriggerOutputs.obsVsDeltaPhi[kRing].pooledMaster, "pRingVsDeltaPhi");
    cloneInto(paperDir, trueTriggerOutputs.obsVsDeltaPhi[kRing].evtMeanMaster, "pRingVsDeltaPhi_EvtMean");
    cloneInto(paperDir, trueTriggerOutputs.ringFlat, "pRingFlat");
    cloneInto(paperDir, trueTriggerOutputs.ringFlatEvtMean, "pRingFlat_EvtMean");
    cloneInto(paperDir, trueTriggerOutputs.obsD3N[kRing], "pRingD3N");
    cloneInto(paperDir, trueTriggerOutputs.obsD3NEvtMean[kRing], "pRingD3N_EvtMean");

    // Warnings: printed and stored, nothing is skipped silently
    std::vector<std::string> warningLines = warnings.startupMessages;
    if (warnings.invalidCells > 0)
      warningLines.push_back(std::to_string(warnings.invalidCells) + " cell-events with non-finite or negative columns were skipped");
    if (warnings.emptyDeltaPhiBins > 0)
      warningLines.push_back(std::to_string(warnings.emptyDeltaPhiBins) + " empty Delta phi bins inside the cuts (the grid binning should make this impossible: check the grid)");
    if (warnings.pointsOutsideMaps > 0)
      warningLines.push_back(std::to_string(warnings.pointsOutsideMaps) + " map points (one per sub-cell point and event) fell outside PolarizationMaps/ and were left out of it");
    if (warnings.degenerateRingDirection > 0)
      warningLines.push_back(std::to_string(warnings.degenerateRingDirection) + " cell-triggers with p parallel to t (n undefined) were skipped");
    if (warnings.zeroPolStarCells > 0)
      warningLines.push_back(std::to_string(warnings.zeroPolStarCells) + " cell-events with |P*| = 0 were left out of the alignment P_R/|P*|");
    std::string warningSummary;
    for (const std::string& warningLine : warningLines) {
      std::cout << "[WARNING] " << warningLine << std::endl;
      warningSummary += warningLine + "\n";
    }

    // Config/: every option as (name = option, title = value), plus provenance
    configDir->cd();
    for (const auto& [name, value] : configEntries)
      TNamed(name.c_str(), value.c_str()).Write();
    TNamed("columnToPolarization", (std::to_string(kColumnToPolarization) + " -- the iSS S^mu columns already hold P^mu = S^mu/(1/2): see README, IMPORTANT").c_str()).Write();
    const std::string gridSummary = std::to_string(grid.nY) + " y nodes from " + std::to_string(grid.yFirst) + " (step " + std::to_string(grid.yWidth) + "), " + std::to_string(grid.nPt) + " pT nodes from " + std::to_string(grid.ptFirst) + " (step " + std::to_string(grid.ptWidth) + "), " + std::to_string(grid.nPhi) + " phi nodes from " + std::to_string(grid.phiFirst) + " (step 2pi/" + std::to_string(grid.nPhi) + "), " + std::to_string(nInCutCells) + " cells in the cuts";
    TNamed("grid", gridSummary.c_str()).Write();
    TNamed("eventsProcessed", std::to_string(eventsProcessed).c_str()).Write();
    TNamed("skippedEvents", skippedEventsText.c_str()).Write();
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