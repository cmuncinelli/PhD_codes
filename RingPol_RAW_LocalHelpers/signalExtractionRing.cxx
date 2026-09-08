/**************************************************************************************************
 *
 * signalExtractionRing.cxx
 *
 * Authors: Cicero & Gemini (for coding speed, obviously. Everything was checked and edited when needed!)
 *
 * Description:
 * ----------------------------------------------------------------------------------------------
 * Signal extraction framework for the Ring Observable polarization analysis in Lambda--jet systems.
 * Expected usage: root -l -b -q 'signalExtractionRing.cxx("/home/cicero/RingPol/LHC25ae_pass2_small/results_consumer/ConsumerResults_JustTOFBase-FullQA_permissivePt-BothHyperons.root",
  "/home/cicero/RingPol/LHC25ae_pass2_small/results_consumer/")'
 *
 * This code performs invariant mass fits in bins of angular and/or kinematic observables in order
 * to extract the true signal polarization <R>_signal from a mixture of signal and combinatorial
 * background Lambda candidates.
 *
 * The necessity of signal extraction arises from the fact that the measured polarization inside
 * the Lambda invariant mass peak is not purely signal:
 * <R>_measured = f_S <R>_signal + f_B <R>_background
 * where background polarization is NOT guaranteed to be zero.
 *
 * In heavy-ion collisions, local angular momentum can generate polarization.
 * However, combinatorial background may:
 * - Dilute the signal
 * - Be polarized differently
 * - Have distorted angular structure
 *
 * Therefore, proper signal extraction is mandatory.
 *
 * IMPORTANT NON-TRIVIAL ASSUMPTIONS:
 * ----------------------------------------------------------------------------------------------
 * 1) Background polarization is assumed smooth in invariant mass.
 * 2) No strong mass--observable correlation exists (must be validated).
 * 3) Sidebands reliably represent the background under the peak.
 * If these assumptions fail, a simultaneous mass--polarization fit will need to be implemented:
 * Numerator(m)   = S(m) R_S + B(m) R_B
 * Denominator(m) = S(m) + B(m)
 * which avoids sideband assumptions entirely.
 *
 * INVOCATION:
 * ----------------------------------------------------------------------------------------------
 *   signalExtractionRing <inputFilePath> <outputFolderPath> [--key=value ...]
 * Run with --help for the full list of options and their current defaults. The extraction windows
 * (signal half-width, sideband inner and outer edges) and the background polynomial order are all
 * set from there; nothing about a region is hardcoded anywhere below, including in the plot labels.
 * With no options the historical windows are reproduced, so a default run stays comparable with
 * older output. See the TODO: block on MakeSharedSidebandConfig for where the defaults are headed.
 **************************************************************************************************/
#include <iostream>
#include <string>
#include <vector>
#include <utility>  // std::pair, used for the selection-flow chain segments
#include <algorithm> // std::min, std::max
#include <cmath>
#include <stdexcept> // std::exception, thrown by std::stod on a malformed --key=value
#include <filesystem> // To create output folder whenever it is not present

// ROOT Includes
#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TMatrixDSym.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TString.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>

#include <TBox.h>
#include <TLine.h>
#include <TLegend.h>
#include <TStyle.h>      // gStyle, which must be set before anything is painted
#include <TPaveStats.h>  // The stats box, edited on the object so it survives serialisation
#include <TPaveText.h>   // The formula / parameter annotations on the QA canvases
#include <TArrayD.h>     // Variable bin edges, when cloning a profile's axis

// Extra includes for simultaneous signal extraction fit to Ring Observable numerator and Counts denominator:
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

namespace fs = std::filesystem;

constexpr double lambdaPDGMassApprox = 1.11568; // Just for some fit initial guesses

// ------------------------------------------------------------------------------------------------
// V0 selection flow layout -- mirrors V0SelectionFlowCounter in the producer task.
// CAUTION! If the selection order changes in the task, these must change here too.
//
//   x = 0 .. 30    generic V0 cuts, shared by both mass hypotheses  -> ROOT bins 1 .. 31
//   x = 31 .. 41   Lambda-specific cuts                             -> ROOT bins 32 .. 42
//   x = 42 .. 52   AntiLambda-specific cuts                         -> ROOT bins 43 .. 53
//
// Note that the producer's fill() writes BOTH h2dSelectionLambdaMass and h2dSelectionAntiLambdaMass
// at every step, so each 2D histogram physically contains all 53 bins. But the two hypothesis
// blocks are NOT nested in one another: a candidate goes down one branch or the other, and only the
// generic block is a common ancestor of both. Reading the AntiLambda block out of the Lambda-mass
// histogram (or vice versa) means looking at candidates selected under the wrong hypothesis, and
// chaining across the branch boundary would break the nesting assumption that the retention
// uncertainties rely on. Each species is therefore processed along its own chain:
//   Lambda     : bins 1..31 then 32..42
//   AntiLambda : bins 1..31 then 43..53   (the Lambda block is skipped entirely)
// ------------------------------------------------------------------------------------------------
constexpr int kNGenericV0Cuts = 31;  // Shared cuts, occupying ROOT bins 1 .. kNGenericV0Cuts
constexpr int kNHypothesisCuts = 11; // Cuts per mass hypothesis, in one contiguous block each


// ================================================================================================
// SHARED SIDEBAND SIGNAL EXTRACTION
// ================================================================================================
// Four places in this file perform the same bin-counting + sideband-subtraction measurement:
//   1. ExtractObservable2D, per angular bin              (numerator available -> yields R_S)
//   2. ExtractObservable2D, angle-integrated             (numerator available -> yields R_S)
//   3. PerformDenominatorQA                              (denominator only)
//   4. PerformSelectionCutFlowExtraction, per cut         (denominator only)
//
// Each used to carry its own copy of the algorithm, which is how the same error-propagation bug
// ended up needing the same fix in three separate places. The three pieces that are genuinely
// identical now live here exactly once:
//   FitMassPeak()             -- locate the peak (mu, sigma) with a gaus + pol2 fit
//   BuildSidebandGraph()      -- collect the counts-density sideband points
//   BuildRingSidebandGraph()  -- collect the <R>(m) sideband points
//   ComputePeakWindowYields() -- count, subtract the background, and propagate every uncertainty
//
// What is NOT unified is the orchestration around them. The four workflows really do differ in
// what they are for, and forcing them into one entry point would have required a pile of mode
// switches that hide the differences rather than document them. Each caller therefore keeps a
// thin, readable sequence of calls into the helpers above.
//
// THE NUMERATOR IS OPTIONAL. Pass hNumSum = nullptr (and ringBkgFit = nullptr) to run in
// "denominator only" mode: the mass spectrum is fitted, the background is subtracted, and the
// yields, purity and significance come out exactly as usual -- only the ring-observable quantities
// (R_peak, R_B, R_S) are left untouched, with hasNumerator = false to say so. That is what lets
// PerformDenominatorQA and PerformSelectionCutFlowExtraction share this code with the ring
// extraction: they measure the same yields, they simply have no <R> to carry along.
// ================================================================================================

// How the initial fit parameters are guessed before the gaus + pol2 mass fit.
enum class PeakInitGuessMode {
    SimpleMaximum,           // amplitude = histogram maximum; pol2 coefficients left at ROOT defaults
    SimpleMaximumZeroedPol2, // same, but the three pol2 coefficients are explicitly set to zero
    ShapeEstimated           // amplitude, mean and pol2 offset estimated from the histogram shape
};

// What to do when the mass fit does not converge cleanly.
enum class PeakFitFallback {
    Reject,        // treat the bin as unusable
    SalvageWindow, // keep (mu, sigma) if they land inside a hand-set "physically sane" window
    PdgDefaults    // fall back to PDG-approximate values and carry on (used for drawing-only QA)
};

// ------------------------------------------------------------------------------------------------
// SidebandConfig -- every tuned number of the extraction, in one place.
// ------------------------------------------------------------------------------------------------
// HISTORY: this used to be four independently tuned presets (per-bin, integrated, denominator QA,
// cut flow), which had drifted apart on things that are not workflow-specific at all: one fitted
// raw counts while the others fitted densities, one kept empty sideband bins, and the sigma limits
// disagreed between two of them. The differences were originally meant to encode different levels
// of available statistics -- the TH3D projections are far thinner than the integrated spectrum --
// but a looser fit is not what thin statistics need. They now share ONE set of fitting and
// stability numbers, chosen as the most numerically defensive of each pair, and differ only where
// the difference is semantic rather than numerical (see MakeExtractionConfigs below).
struct SidebandConfig {
    // --- Stability gates --------------------------------------------------------------------
    double minEntries;        // Require hMass->GetEntries() >= minEntries
    double minIntegral;       // Require hMass->Integral()   >  minIntegral
    int    minSidebandPoints; // Minimum points in the sideband graph before attempting the pol fit
    double minSidebandCounts; // Minimum RAW counts in the sidebands (<= 0 disables the check)

    // --- Peak-finding fit -------------------------------------------------------------------
    PeakInitGuessMode initGuessMode;
    PeakFitFallback   fallback;
    const char*       fitOptions;                 // ROOT fit option string
    bool              constrainAmplitudePositive; // SetParLimits(0, 0, 1e18)
    double            muInitGuess;                // Starting guess for the peak position
    double            sigmaInitGuess;             // Starting guess for the Gaussian width
    double            muLimitWindow;              // SetParLimits(1, muInitGuess -/+ muLimitWindow)
    double            sigmaLimitMin, sigmaLimitMax; // SetParLimits(2, ...)

    // --- Salvage window (only consulted when fallback == SalvageWindow) ----------------------
    double salvageMuMin, salvageMuMax;
    double salvageSigmaMin, salvageSigmaMax;

    // --- Region definitions, in units of the fitted peak sigma -------------------------------
    // Signal window   : [mu - nSigmaPeak*sigma, mu + nSigmaPeak*sigma]
    // Sideband band   : [mu -/+ nSigmaExclusion*sigma, mu -/+ nSigmaExclusionOuter*sigma]
    double nSigmaPeak;           // Half-width of the signal counting window
    double nSigmaExclusion;      // Inner edge of the sidebands (half-width of the excluded zone)
    double nSigmaExclusionOuter; // Outer edge of the sidebands.
                                 // <= 0 means "consider the edges of the mass axis"

    // --- Background model -------------------------------------------------------------------
    // Order of the "polN" fitted to the sideband graphs. The two spectra are NOT the same shape
    // and have no business sharing an order: the COUNTS sideband is a smooth combinatorial mass
    // spectrum with real curvature, while the NUMERATOR is Sum_R against mass, whose background is
    // close to featureless. A pol2 there fits noise and then extrapolates it under the peak, which
    // is the worst place for a spurious quadratic to land.
    int bkgPolOrder;     // Counts sideband, fitted as a density in (m - muInitGuess)
    int ringBkgPolOrder; // <R>_bkg(m) in the sidebands. NOT a density: <R> is intensive.

    // --- Sideband collection ----------------------------------------------------------------
    bool requirePositiveSidebandBins; // Skip empty / zero-error bins when building the graph
};

/// @brief How ExtractObservable2D reduces its angular bins to one integrated number.
enum class IntegralMode {
    ProjectThenExtract, // Project the 2D inputs onto mass, extract once. Legacy path.
    CombinePerBin       // Skip the projection entirely; combine the per-bin extractions instead.
};

// ------------------------------------------------------------------------------------------------
// The one shared set of numbers. Everything here is a physics or numerical-safety choice that has
// no reason to differ between the four workflows.
// ------------------------------------------------------------------------------------------------
SidebandConfig MakeSharedSidebandConfig()
{
    SidebandConfig c;

    // Stability gates: the strictest of the values the four old presets used, so a bin that would
    // have been fitted by only the loosest of them is now skipped everywhere.
    c.minEntries = 50.0;
    c.minIntegral = 30.0;
    c.minSidebandCounts = 8.0; // Below this, the matrix inversion floods Minuit with errors
    c.minSidebandPoints = 4;   // 0 means "derive from bkgPolOrder"; resolved in ResolveDerivedFields

    // Peak-finding fit. ShapeEstimated gives the Gaussian only the excess above an estimated
    // background level instead of the raw maximum, which is the behaviour that survives a
    // background-dominated spectrum -- exactly the regime the thin TH3D slices live in.
    c.initGuessMode = PeakInitGuessMode::ShapeEstimated;
    c.fallback = PeakFitFallback::Reject; // Overridden per workflow below
    c.fitOptions = "Q 0 R S";
    c.constrainAmplitudePositive = true;
    c.muInitGuess = lambdaPDGMassApprox;
    c.sigmaInitGuess = 0.002; // Typical for Lambda, from an earlier analysis
    c.muLimitWindow = 0.01;

    // Sigma limits: the two old presets disagreed (0.1-5 MeV vs 0.5-8 MeV). The wider one is kept
    // because it is the only one under which the salvage window below is actually reachable; with
    // a 5 MeV cap, SetParLimits already forbids everything the salvage test was written to accept.
    c.sigmaLimitMin = 0.0005; // 0.5 MeV
    c.sigmaLimitMax = 0.008;  // 8 MeV
    c.salvageMuMin = 1.105;
    c.salvageMuMax = 1.125;
    c.salvageSigmaMin = 0.0005;
    c.salvageSigmaMax = 0.008;

    // Region definitions.
    c.nSigmaPeak = 2.0; // Historic value: 4.0
    c.nSigmaExclusion = 5.0; // Historic value: 6.0
    c.nSigmaExclusionOuter = 7.0; // Historic value: -1.0

    // These two order parameters should be kept at the same values for the hNumExtBinDensity fits
    // For the density fit, we are fitting <R>*dN/dm (i.e., <R>_meas * (N_sig + N_bkg) = (<R>_sig * N_sig + <R>_bkg * N_bkg)), from <R>_meas = (<R>_sig * N_sig + <R>_bkg * N_bkg)/(N_sig + N_bkg)
    // Thus, the shape is conditioned to the dN/dm format and we should have the same pol2 estimating both of them.
    c.bkgPolOrder = 2;    // The Lambda combinatorial background is almost linear, with a little curvature
    // pol0 says exactly what sideband subtraction has always assumed -- that <R> in the sidebands
    // does not depend on mass. Raising it to 1 is the natural systematic variation, and the QA
    // canvas shows directly whether the flat assumption holds.
    c.ringBkgPolOrder = 0;
    c.requirePositiveSidebandBins = true;

    return c;
}

// Fields that are derived from other fields rather than set directly.
void ResolveDerivedFields(SidebandConfig& c)
{
    // A polN has (N+1) parameters. Asking for (N+1)+2 points leaves 2 degrees of freedom, which is
    // the minimum at which the covariance matrix TF1::IntegralError needs is worth anything.
    // Driven by the HIGHER of the two orders, since both fits run on the same set of points.
    if (c.minSidebandPoints <= 0)
        c.minSidebandPoints = std::max(c.bkgPolOrder, c.ringBkgPolOrder) + 3;
}

// ------------------------------------------------------------------------------------------------
// The four workflow configurations. They now differ ONLY in their failure policy, which is genuinely
// semantic: it says what the caller wants done with a bin whose peak fit did not converge, not how
// carefully the fit should have been attempted.
// ------------------------------------------------------------------------------------------------
struct ExtractionConfigSet {
    SidebandConfig perBin;     // ExtractObservable2D, per angular bin
    SidebandConfig integrated; // ExtractObservable2D, angle-integrated
    SidebandConfig denomQA;    // PerformDenominatorQA
    SidebandConfig cutFlow;    // PerformSelectionCutFlowExtraction
};

ExtractionConfigSet MakeExtractionConfigs()
{
    SidebandConfig base = MakeSharedSidebandConfig();
    ExtractionConfigSet s;

    // Per angular bin: a non-converged fit whose (mu, sigma) still land somewhere physically sane
    // is worth keeping, because only the central values are consumed downstream and a non-zero
    // Minuit status is usually just a complaint about the covariance matrix.
    s.perBin = base;
    s.perBin.fallback = PeakFitFallback::SalvageWindow;

    // Angle-integrated. Salvaged rather than rejected: Minuit status 3 means the error matrix had
    // to be forced positive-definite, which happens routinely on this spectrum and says nothing
    // about whether (mu, sigma) are sensible. Rejecting on it produced zero-filled results that
    // were indistinguishable from a genuine measurement of zero -- see the status reporting below,
    // which now makes any salvage or failure visible in the log and in hExtractionStatus.
    s.integrated = base;
    s.integrated.fallback = PeakFitFallback::SalvageWindow;

    // Denominator QA: drawing-only. It must always produce a canvas, so a failed fit falls back to
    // PDG-approximate values purely to keep the drawn regions meaningful. Yields are reported only
    // when the fit actually converged -- see the guard in PerformDenominatorQA.
    s.denomQA = base;
    s.denomQA.fallback = PeakFitFallback::PdgDefaults;

    // Cut flow: the early selection steps sit under an enormous background, so salvaging is the
    // right call for the same reason as the per-bin path.
    s.cutFlow = base;
    s.cutFlow.fallback = PeakFitFallback::SalvageWindow;

    ResolveDerivedFields(s.perBin);
    ResolveDerivedFields(s.integrated);
    ResolveDerivedFields(s.denomQA);
    ResolveDerivedFields(s.cutFlow);
    return s;
}

// ------------------------------------------------------------------------------------------------
// Gaussian coverage of an arbitrary (possibly asymmetric) window around a fitted peak.
// ------------------------------------------------------------------------------------------------
// Returns the fraction of a Gaussian(mu, sigma) that lies inside [xLow, xHigh]:
//   0.5 * ( erf((xHigh-mu)/(sigma*sqrt2)) - erf((xLow-mu)/(sigma*sqrt2)) )
// For a symmetric window at n sigma this reduces to erf(n/sqrt2), which is where the old hardcoded
// 0.9999366 came from (it is erf(4/sqrt2)). That constant silently assumed a 4-sigma window, so it
// had to go the moment the window became configurable: at 2 sigma the true coverage is 0.9545 and
// the yield would have been biased by 4.7%. Using the realized bin-snapped edges rather than
// mu +/- n*sigma also handles the asymmetry that finite binning introduces.
double GaussianWindowCoverage(double xLow, double xHigh, double mu, double sigma)
{
    if (sigma <= 0.0) return 0.0;
    const double invSqrt2 = 1.0 / std::sqrt(2.0);
    return 0.5 * (std::erf((xHigh - mu) / sigma * invSqrt2) - std::erf((xLow - mu) / sigma * invSqrt2));
}

// ------------------------------------------------------------------------------------------------
// EnsureDir -- fetch a subdirectory, creating it only if it is not there yet.
// ------------------------------------------------------------------------------------------------
// ExtractObservable2D is called many times with the same parent directory (once per angular
// observable, and once per kinematic window inside each 3D slice folder), and every one of them
// wants to drop something into a shared IntegratedSummary/ folder. TDirectory::mkdir on an
// existing name returns nullptr, so a plain mkdir would work for the first caller and silently
// fail for all the others.
TDirectory* EnsureDir(TDirectory* parent, const char* name)
{
    if (!parent) return nullptr;
    TDirectory* d = parent->GetDirectory(name);
    return d ? d : parent->mkdir(name);
}

// ------------------------------------------------------------------------------------------------
// Label helpers -- so that no plot can ever claim a window the extraction did not use.
// ------------------------------------------------------------------------------------------------
// Every canvas title, axis label and legend entry that mentions a region is built from these, which
// is the whole point: the strings used to be typed out as "#pm4#sigma" and "6#sigma excl." next to
// code that read the number from a variable, so changing the variable quietly made the plots lie.
/// @brief A polN written in (x - shift) rather than in x.
///
/// @note The {1, x, x^2} basis is nearly degenerate over a 70 MeV range centred on 1.115: the three
/// columns are almost parallel, the fitted coefficients come out around 1e12 with alternating
/// signs, and the value at the peak is a cancellation down to about 1e9. Double precision absorbs
/// the value, but the PARAMETER COVARIANCE does not survive it -- TF1::IntegralError integrates a
/// gradient built from that matrix and was returning "cannot reach tolerance because of roundoff error"
/// with an integral error of several thousand counts. Var(bkgCounts) feeds every downstream
/// uncertainty, so this was a correctness problem and not merely cosmetic.
///
/// Recentring on the peak makes the basis {1, dm, dm^2} with dm of order 1e-2, the coefficients
/// come out at the scale of the density itself, and the cancellation disappears.
TString PolShifted(int order, double shift)
{
    TString f = "[0]";
    for (int k = 1; k <= order; ++k) {
        TString term = Form("[%d]", k);
        for (int p = 0; p < k; ++p) term += Form("*(x-%.6f)", shift);
        f += "+" + term;
    }
    return f;
}

/// @brief Sideband graph of the ring observable itself, <R>(m), rather than of Sum_R.
///
/// WHY NOT Sum_R. In the sidebands d(Sum_R)/dm = <R>_bkg(m) * dN_bkg/dm. With dN_bkg/dm modelled as
/// a pol2 and <R>_bkg roughly constant, that product is a pol2 as well -- so fitting Sum_R with a
/// pol1 asserted that a quadratic times a constant is linear, and the fit absorbed the mismatch by
/// forcing <R>_bkg(m) = pol1(m)/pol2(m). That ratio has a pole wherever the denominator crosses
/// zero, which is what made the drawn background look like nothing physical.
///
/// Fitting <R> directly removes the misspecification entirely. <R> is INTENSIVE, so unlike counts
/// and unlike Sum_R it must NOT be divided by the bin width: it does not scale with binning.
/// A pol0 then says exactly what sideband subtraction always assumed, and the assumption becomes
/// visible on the QA canvas instead of being buried in a ratio of two fits.
///
/// @param hCounts Candidate counts vs mass (raw, not a density).
/// @param hNumSum Sum_R vs mass (raw, not a density). <R> = hNumSum/hCounts bin by bin.
TGraphErrors* BuildRingSidebandGraph(TH1D* hCounts, TH1D* hNumSum, double mu, double sigma,
                                     double massMin, double massMax, const SidebandConfig& cfg,
                                     TString grName)
{
    if (!hCounts || !hNumSum) return nullptr;

    const bool boundedOuter = (cfg.nSigmaExclusionOuter > 0.0);
    const double xInnerLow = mu - cfg.nSigmaExclusion * sigma;
    const double xInnerHigh = mu + cfg.nSigmaExclusion * sigma;
    const double xOuterLow = boundedOuter ? std::max(massMin, mu - cfg.nSigmaExclusionOuter * sigma) : massMin;
    const double xOuterHigh = boundedOuter ? std::min(massMax, mu + cfg.nSigmaExclusionOuter * sigma) : massMax;

    TGraphErrors* gr = new TGraphErrors();
    gr->SetName(grName);

    int ptIdx = 0;
    for (int jBin = 1; jBin <= hCounts->GetNbinsX(); ++jBin) {
        const double x = hCounts->GetBinCenter(jBin);
        const bool inLeft = (x >= xOuterLow && x <= xInnerLow);
        const bool inRight = (x >= xInnerHigh && x <= xOuterHigh);
        if (!inLeft && !inRight) continue;

        const double n = hCounts->GetBinContent(jBin);
        const double errSum = hNumSum->GetBinError(jBin);
        if (!(n > 0.0) || !(errSum > 0.0)) continue; // Content may be zero or negative; the error may not

        gr->SetPoint(ptIdx, x, hNumSum->GetBinContent(jBin) / n);
        gr->SetPointError(ptIdx, 0.0, errSum / n); // errSum = sigma_R*sqrt(N), so this is the SEM
        ptIdx++;
    }
    return gr;
}

TString SignalWindowLabel(const SidebandConfig& cfg)
{
    return TString(Form("#mu#pm%.3g#sigma", cfg.nSigmaPeak));
}

TString SidebandBandLabel(const SidebandConfig& cfg)
{
    if (cfg.nSigmaExclusionOuter > 0.0)
        return TString(Form("%.3g-%.3g#sigma", cfg.nSigmaExclusion, cfg.nSigmaExclusionOuter));
    return TString(Form("%.3g#sigma excl.", cfg.nSigmaExclusion));
}

// ------------------------------------------------------------------------------------------------
// FitMassPeak -- locate the Lambda peak with a combined gaus(0) + pol2(3) fit.
// ------------------------------------------------------------------------------------------------
// Only (mu, sigma) are consumed downstream: the pol2 half of this fit is thrown away, and the
// background is refitted on the sidebands alone by BuildSidebandGraph + a separate pol2. Fitting
// both together here is purely a robust way of finding the peak in the presence of background.
//
// RETURNS: outcome.fit is a NEW TF1 owned by the caller. Attach it to a histogram's list of
// functions (which transfers ownership) or delete it.
struct PeakFitOutcome {
    bool   valid = false;
    bool   converged = false; // Fit itself converged, independently of any fallback that was used
    double mu = 0.0;
    double sigma = 0.0;
    int    status = -1;
    TF1*   fit = nullptr;
};

PeakFitOutcome FitMassPeak(TH1D* hForFit, double massMin, double massMax,
                           const SidebandConfig& cfg, TString fitName)
{
    PeakFitOutcome out;
    if (!hForFit) return out;

    TF1* fitFunc = new TF1(fitName, "gaus(0) + pol2(3)", massMin, massMax);

    // --- Initial guesses ---
    double maxVal = hForFit->GetMaximum();
    switch (cfg.initGuessMode) {
        case PeakInitGuessMode::SimpleMaximum:
            fitFunc->SetParameter(0, maxVal);            // Amplitude
            fitFunc->SetParameter(1, cfg.muInitGuess);   // Mean close to the PDG value
            fitFunc->SetParameter(2, cfg.sigmaInitGuess);
            break;
        case PeakInitGuessMode::SimpleMaximumZeroedPol2:
            // Explicitly sets the pol2 coefficients to zero as well
            fitFunc->SetParameters(maxVal, cfg.muInitGuess, cfg.sigmaInitGuess, 0.0, 0.0, 0.0);
            break;
        case PeakInitGuessMode::ShapeEstimated: {
            // Estimate the background level from the two extreme bins and give the Gaussian only
            // the excess above it, which behaves much better on a very background-dominated
            // spectrum than starting the amplitude at the raw maximum.
            double estMu = hForFit->GetBinCenter(hForFit->GetMaximumBin());
            double firstBinC = hForFit->GetBinContent(1);
            double lastBinC = hForFit->GetBinContent(hForFit->GetNbinsX());
            double estBkg = (firstBinC + lastBinC) / 2.0;
            if (estBkg < 0.0) estBkg = 0.0;
            double estSigAmp = maxVal - estBkg;
            if (estSigAmp < 0.0) estSigAmp = maxVal * 0.5; // Fallback for oddly-shaped histograms
            fitFunc->SetParameter(0, estSigAmp);
            fitFunc->SetParameter(1, estMu);
            fitFunc->SetParameter(2, cfg.sigmaInitGuess);
            fitFunc->SetParameter(3, estBkg);
            fitFunc->SetParameter(4, 0.0);
            fitFunc->SetParameter(5, 0.0);
            break;
        }
    }

    // --- Stability constraints ---
    if (cfg.constrainAmplitudePositive) fitFunc->SetParLimits(0, 0.0, 1e18);
    fitFunc->SetParLimits(1, cfg.muInitGuess - cfg.muLimitWindow,
                             cfg.muInitGuess + cfg.muLimitWindow);
    fitFunc->SetParLimits(2, cfg.sigmaLimitMin, cfg.sigmaLimitMax);

    TFitResultPtr r = hForFit->Fit(fitFunc, cfg.fitOptions);

    out.fit = fitFunc;
    out.status = int(r);
    out.mu = fitFunc->GetParameter(1);
    out.sigma = fitFunc->GetParameter(2);

    // --- Acceptance ---
    switch (cfg.fallback) {
        case PeakFitFallback::SalvageWindow:
            // Status 0 means a clean convergence. Anything else is usually just Minuit complaining
            // about the covariance (error) matrix, which does not matter here because only the
            // central values of mu and sigma are used. So a non-zero status is salvaged as long as
            // the parameters themselves are physically sensible.
            out.converged = (out.status == 0);
            out.valid = out.converged ||
                        (out.mu > cfg.salvageMuMin && out.mu < cfg.salvageMuMax &&
                         out.sigma > cfg.salvageSigmaMin && out.sigma < cfg.salvageSigmaMax);
            if (!out.valid) { out.mu = 0.0; out.sigma = 0.0; }
            break;

        case PeakFitFallback::Reject:
            out.converged = r->IsValid();
            out.valid = out.converged;
            break;

        case PeakFitFallback::PdgDefaults:
            // Never fails: the caller still wants to draw something. Region definitions simply fall
            // back to PDG-approximate values so the sidebands remain meaningful.
            out.converged = r->IsValid();
            out.valid = true;
            if (!out.converged) {
                out.mu = cfg.muInitGuess;
                out.sigma = cfg.sigmaInitGuess;
            }
            break;
    }
    return out;
}

// ------------------------------------------------------------------------------------------------
// BuildSidebandGraph -- collect the mass bins inside the sideband band.
// ------------------------------------------------------------------------------------------------
// The points are stored as DENSITIES (content / bin width). This matters because the mass axis has
// variable bin widths: fitting densities makes the polN bin-agnostic, so TF1::Integral over the
// signal window returns counts directly, with no bin-width bookkeeping afterwards.
//
// A TGraphErrors is used rather than a TH1 with a masked range because ROOT histograms handle
// discontinuous fit ranges badly.
//
// The band is [mu -/+ nSigmaExclusion*sigma, mu -/+ nSigmaExclusionOuter*sigma], clipped to
// [massMin, massMax]. When nSigmaExclusionOuter <= 0 the outer edge is the mass axis edge itself,
// which is the historical behaviour: an earlier version used a fixed [4 sigma, 6 sigma] band and
// was widened because that band left only about 10 bins per side.
//
// kind:
//   Counts     -- the denominator spectrum. Empty and zero-error bins are skipped when the config
//                 asks for it, and the raw (unscaled) counts are accumulated into span.rawCounts
//                 for the caller's stability check.
//   Numerator  -- the Sum_R_i spectrum. Only the error is checked, never the content: <R> may
//                 legitimately be zero or negative, so a content > 0 cut would bias the sidebands.
//
// spanOut (optional) reports what was actually taken, so that a "too few points" failure can be
// diagnosed instead of silently invalidating a bin.
//
// RETURNS: a NEW TGraphErrors owned by the caller.
enum class SidebandGraphKind { Counts, Numerator };

struct SidebandSpan {
    int    nPointsLeft = 0, nPointsRight = 0;
    double rawCounts = 0.0;      // Raw (undivided) counts, Counts kind only
    // Requested band edges in mass units, before any bin snapping. Zero when the side is empty.
    double xInnerLow = 0.0, xOuterLow = 0.0;   // Left sideband:  [xOuterLow, xInnerLow]
    double xInnerHigh = 0.0, xOuterHigh = 0.0; // Right sideband: [xInnerHigh, xOuterHigh]
    // Achieved band edges in units of sigma, measured from the fitted mu on each side. These are
    // the bin centres of the extreme points actually collected, so they say what the polN really
    // saw rather than what was asked for.
    double nSigmaInnerLowAchieved = 0.0, nSigmaOuterLowAchieved = 0.0;
    double nSigmaInnerHighAchieved = 0.0, nSigmaOuterHighAchieved = 0.0;
};

TGraphErrors* BuildSidebandGraph(TH1D* hSource, double mu, double sigma,
                                 double massMin, double massMax,
                                 const SidebandConfig& cfg, SidebandGraphKind kind,
                                 SidebandSpan* spanOut, TString grName)
{
    SidebandSpan span;
    if (spanOut) *spanOut = span;
    if (!hSource) return nullptr;

    // Resolve the band. A non-positive outer half-width means "run to the mass axis edge".
    const bool boundedOuter = (cfg.nSigmaExclusionOuter > 0.0);
    span.xInnerLow = mu - cfg.nSigmaExclusion * sigma;
    span.xInnerHigh = mu + cfg.nSigmaExclusion * sigma;
    span.xOuterLow = boundedOuter ? std::max(massMin, mu - cfg.nSigmaExclusionOuter * sigma) : massMin;
    span.xOuterHigh = boundedOuter ? std::min(massMax, mu + cfg.nSigmaExclusionOuter * sigma) : massMax;

    TGraphErrors* gr = new TGraphErrors();
    gr->SetName(grName);

    int ptIdx = 0;
    for (int jBin = 1; jBin <= hSource->GetNbinsX(); ++jBin) {
        double x = hSource->GetBinCenter(jBin);

        bool inLeftSideband = (x >= span.xOuterLow && x <= span.xInnerLow);
        bool inRightSideband = (x >= span.xInnerHigh && x <= span.xOuterHigh);
        if (!inLeftSideband && !inRightSideband) continue;

        double binContent = hSource->GetBinContent(jBin);
        double binErr = hSource->GetBinError(jBin);
        double bw = hSource->GetBinWidth(jBin);

        if (cfg.requirePositiveSidebandBins) {
            if (kind == SidebandGraphKind::Counts) {
                if (!(binContent > 0.0 && binErr > 0.0)) continue;
            } else {
                if (!(binErr > 0.0)) continue; // Content may be zero or negative for <R>
            }
        }

        if (kind == SidebandGraphKind::Counts) span.rawCounts += binContent;

        // Track the extreme points actually taken on each side, in sigma units from mu.
        if (sigma > 0.0) {
            if (inLeftSideband) {
                double nSig = (mu - x) / sigma;
                if (span.nPointsLeft == 0 || nSig < span.nSigmaInnerLowAchieved) span.nSigmaInnerLowAchieved = nSig;
                if (nSig > span.nSigmaOuterLowAchieved) span.nSigmaOuterLowAchieved = nSig;
                span.nPointsLeft++;
            } else {
                double nSig = (x - mu) / sigma;
                if (span.nPointsRight == 0 || nSig < span.nSigmaInnerHighAchieved) span.nSigmaInnerHighAchieved = nSig;
                if (nSig > span.nSigmaOuterHighAchieved) span.nSigmaOuterHighAchieved = nSig;
                span.nPointsRight++;
            }
        }

        gr->SetPoint(ptIdx, x, binContent / bw);
        gr->SetPointError(ptIdx, 0.0, binErr / bw);
        ptIdx++;
    }

    if (spanOut) *spanOut = span;
    return gr;
}

// ------------------------------------------------------------------------------------------------
// ComputePeakWindowYields -- count in the signal window, subtract the fitted background, and
// propagate every uncertainty.
// ------------------------------------------------------------------------------------------------
// This is the single place in the file where the extraction arithmetic lives. Every formula below
// is derived once, here, and every caller inherits it.
//
// The numerator arguments (hNumSum, ringBkgFit, rRingBkg) are OPTIONAL: pass nullptr for hNumSum
// and ringBkgFit to run in denominator-only mode. Yields, purity and significance are produced either
// way; R_peak, R_B and R_S are filled only when a numerator was supplied, flagged by hasNumerator.
struct PeakWindowYields {
    bool valid = false;
    int  failStage = 0; // 0 = ok, 1 = non-positive totCounts, 2 = non-positive sigCounts

    // Signal window actually used, snapped to bin edges
    double xLow = 0.0, xHigh = 0.0;
    int firstBin = 0, lastBin = 0;
    // What that window is worth in sigma units, measured from the fitted mu on each side. The two
    // are NOT equal in general: FindBin snaps to bin edges, and the mass axis has variable bin
    // widths, so a symmetric request comes back asymmetric. These are the numbers to quote when
    // saying how wide the signal window really was.
    double nSigmaAchievedLow = 0.0, nSigmaAchievedHigh = 0.0;
    double coverage = 0.0; // Gaussian fraction inside [xLow, xHigh]; erf(n/sqrt2) if symmetric

    // Yields
    double totCounts = 0.0, totCountsErrSq = 0.0;
    double bkgCounts = 0.0, errBkgCounts = 0.0;
    double sigCounts = 0.0, errSigCounts = 0.0;

    // Fractions and QA
    double fB = 0.0, fS = 0.0, var_fB = 0.0;
    double purity = 0.0, errPurity = 0.0;
    double significance = 0.0, errSignificance = 0.0;

    // Ring observable (only when a numerator was supplied)
    bool hasNumerator = false;
    double totNum = 0.0, totNumErrSq = 0.0;
    double bkgNum = 0.0, errBkgNum = 0.0;
    double R_peak = 0.0, errR_peak = 0.0;
    double R_B = 0.0, varR_B = 0.0, errR_B = 0.0;
    double R_S = 0.0, errR_S = 0.0;

    // Differences between the three ring observables, with the correlations kept.
    // These CANNOT be recomputed downstream by adding the errors above in quadrature: the three
    // are built from the same four primitives and are strongly correlated. See the derivation at
    // the bottom of ComputePeakWindowYields.
    double diffPeakMinusSig = 0.0, errDiffPeakMinusSig = 0.0;
    double diffPeakMinusBkg = 0.0, errDiffPeakMinusBkg = 0.0;
    double diffSigMinusBkg = 0.0, errDiffSigMinusBkg = 0.0;
};

// ------------------------------------------------------------------------------------------------
// FinalizeDerivedQuantities -- everything that follows from the four primitives.
// ------------------------------------------------------------------------------------------------
// Takes a PeakWindowYields whose primitives (totCounts, totNum, bkgCounts, bkgNum and their
// variances) are already filled, and computes every derived quantity from them: yields, purity,
// significance, the three ring observables and the three correlated differences.
//
// It was split out of ComputePeakWindowYields so that the ANGLE-COMBINED result can reuse it. Two
// independent angular bins have independent primitives, so summing the primitives and running this
// on the totals is algebraically identical to a single extraction over the union of the bins --
// which means the combined result is guaranteed to use exactly the same error propagation as the
// per-bin one, rather than a second implementation that could drift away from it.
//
// Sets y.valid on success, or y.failStage on rejection.
void FinalizeDerivedQuantities(PeakWindowYields& y)
{
    // --- Signal yield ---
    // totCounts and bkgCounts are independent: bkgCounts comes from a fit to the SIDEBANDS only,
    // while totCounts counts the peak window, and the two mass regions are disjoint (hence disjoint
    // candidates, hence independent Poisson populations).
    y.sigCounts = y.totCounts - y.bkgCounts;
    y.errSigCounts = std::sqrt(y.totCountsErrSq + y.errBkgCounts * y.errBkgCounts);

    if (y.sigCounts <= 0) { y.failStage = 2; return; }

    y.fB = y.bkgCounts / y.totCounts;
    y.fS = y.sigCounts / y.totCounts; // Equivalent to 1 - fB

    // Variance of f_B: a ratio of the two independent quantities discussed just above.
    y.var_fB = std::pow(y.fB, 2) * (std::pow(y.errBkgCounts / y.bkgCounts, 2)
                                    + (y.totCountsErrSq / std::pow(y.totCounts, 2)));

    y.purity = y.fS;
    y.errPurity = std::sqrt(y.var_fB); // Since fS = 1 - fB, the variance is identical

    // --- Significance ---
    // Error propagation for Z = S/sqrt(N). Here N is the MEASURED totCounts, not S + B: the
    // denominator is counted directly and does not move when the background estimate moves.
    // Writing Z in the two independent primitives, Z = (totCounts - bkgCounts)/sqrt(totCounts),
    //   dZ/d(totCounts) = (totCounts + bkgCounts) / (2 * totCounts^(3/2)),
    //   dZ/d(bkgCounts) = -1 / sqrt(totCounts),
    // hence
    //   Var(Z) = (totCounts + bkgCounts)^2 * Var(totCounts)/(4*totCounts^3)
    //            + Var(bkgCounts)/totCounts.
    // An earlier version used sigCounts in place of (totCounts + bkgCounts) AND added
    // Var(sigCounts)/totCounts, which double counts the totCounts fluctuation because sigCounts and
    // totCounts are not independent (Cov = Var(totCounts)). In the zero-background limit the
    // correct expression reproduces the textbook Z = sqrt(N) with sigma_Z = 1/2, while the old one
    // gave 1.12.
    y.significance = y.sigCounts / std::sqrt(y.totCounts);
    y.errSignificance = std::sqrt(std::pow(y.totCounts + y.bkgCounts, 2) * y.totCountsErrSq
                                     / (4.0 * std::pow(y.totCounts, 3))
                                 + y.errBkgCounts * y.errBkgCounts / y.totCounts);

    if (!y.hasNumerator) { y.valid = true; return; }

    // --- Ring observable ---
    y.R_peak = y.totNum / y.totCounts;
    y.R_B = (y.bkgCounts > 0) ? (y.bkgNum / y.bkgCounts) : 0.0;
    // True signal polarization: R_S = (R_peak - f_B * R_B) / f_S
    y.R_S = (y.R_peak - y.fB * y.R_B) / y.fS;

    // --- Error propagation for R_peak (kept as documentation of where totNumErrSq comes from) ---
    // An OLD method used totSqNum (from a separate squared-observable TH2D) to compute
    //   varR_peak = (totSqNum/totCounts - R_peak^2) / (totCounts - 1)   [Bessel-corrected]
    // The current method instead uses the TProfile-derived errors already stored in hNumSum:
    //   sigma^2(Sum_R_i) = totNumErrSq       (quadrature sum of sigma_R*sqrt(N_bin) over the window)
    //   sigma^2(<R>_peak) = totNumErrSq / totCounts^2
    // The two are mathematically equivalent, but the second is consistent with the corrected
    // sideband errors and needs no separate squared-observable histogram.
    //
    // IMPORTANT (this is what makes everything below consistent): totNumErrSq is a CONDITIONAL
    // variance, i.e. the spread of Sum_R_i at FIXED candidate counts. That is exactly what a
    // TProfile error gives, and it is why no covariance correction between totNum and totCounts is
    // needed for R_peak itself. It also fixes the convention used below:
    //   Cov(totNum, totCounts) = R_peak * Var(totCounts),
    //   Cov(bkgNum, bkgCounts) = R_B    * Var(bkgCounts),
    // since fluctuating the candidate population scales the sum and the count together, with the
    // mean as the proportionality factor.

    // --- Error propagation for R_B ---
    // R_B = bkgNum / bkgCounts is a mean, exactly like R_peak, so it obeys the same logic.
    // errBkgNum is the error on a SUM of R values at fixed sideband populations, while errBkgCounts
    // is the error on those populations. The two are correlated: fluctuating the sideband population
    // scales bkgNum and bkgCounts together and leaves their ratio untouched. Inserting
    // Cov(bkgNum, bkgCounts) = R_B * Var(bkgCounts) into the ratio expansion,
    //   Var(R_B)/R_B^2 = var(bkgNum)/bkgNum^2 + var(bkgCounts)/bkgCounts^2
    //                    - 2*Cov(bkgNum,bkgCounts)/(bkgNum*bkgCounts),
    // the last two terms cancel exactly and the result collapses to the same "error of the mean"
    // form already used for R_peak:
    //   Var(R_B) = errBkgNum^2 / bkgCounts^2.
    // An earlier version treated bkgNum and bkgCounts as independent and therefore kept a spurious
    // (errBkgCounts/bkgCounts)^2 term, overestimating the error on R_B.
    y.varR_B = std::pow(y.errBkgNum / y.bkgCounts, 2);
    y.errR_B = std::sqrt(y.varR_B);
    y.errR_peak = std::sqrt(y.totNumErrSq) / y.totCounts; // Same "error of the mean" form

    // --- Error propagation for R_S ---
    // The three-step expression R_S = (R_peak - f_B*R_B)/f_S is, algebraically, nothing but
    //   R_S = (totNum - bkgNum) / (totCounts - bkgCounts),
    // i.e. signal numerator over signal counts. The central value is identical either way, but
    // propagating in this compact form is the correct thing to do: R_peak, R_B and f_B all share the
    // denominators totCounts and bkgCounts, so they are strongly correlated, and the older formula
    // summed their variances as if they were not.
    //
    // Propagating the compact form in the four primitives (totNum, bkgNum, totCounts, bkgCounts)
    // with the two covariance relations stated above, every cross term regroups into a perfect
    // square and the result is simply:
    //
    //   Var(R_S) = [ totNumErrSq + errBkgNum^2
    //                + (R_peak - R_S)^2 * Var(totCounts)
    //                + (R_B    - R_S)^2 * Var(bkgCounts) ] / sigCounts^2
    //
    // The two squared factors are physically transparent: a yield uncertainty only matters to the
    // extent that moving candidates between the signal and background pools actually changes the
    // mean. In the limit R_B -> R_S the background subtraction contributes no uncertainty at all,
    // which the old formula could not reproduce.
    //
    // ASSUMPTIONS: (i) peak window and sidebands are disjoint in mass, so totCounts/totNum are
    // independent of bkgCounts/bkgNum; (ii) the TProfile-derived errors are conditional on the
    // candidate counts, as discussed above; (iii) linear (first-order) propagation, which requires
    // sigCounts to be comfortably away from zero. Assumption (iii) is the binding one: in bins where
    // sigCounts is only a couple of sigma from zero the ratio distribution develops heavy tails and
    // NO first-order formula is reliable.
    double varR_S = (y.totNumErrSq
                     + y.errBkgNum * y.errBkgNum
                     + std::pow(y.R_peak - y.R_S, 2) * y.totCountsErrSq
                     + std::pow(y.R_B - y.R_S, 2) * y.errBkgCounts * y.errBkgCounts)
                    / std::pow(y.sigCounts, 2);
    y.errR_S = std::sqrt(varR_S);

    // --- Differences between the three ring observables ---------------------------------------
    // R_peak, R_S and R_B are three functions of the SAME four primitives, so they are strongly
    // correlated and their differences must be propagated jointly. Adding the errors above in
    // quadrature is not a conservative approximation here: it is wrong in both directions.
    //
    // Working in the unconditional covariance of (totCounts, totNum, bkgCounts, bkgNum),
    //   Var(totNum) = totNumErrSq + R_peak^2 * Var(totCounts),  Cov = R_peak * Var(totCounts)
    //   Var(bkgNum) = errBkgNum^2 + R_B^2    * Var(bkgCounts),  Cov = R_B    * Var(bkgCounts)
    // (the extra terms convert the count-conditional TProfile variances into unconditional ones),
    // the three gradients regroup exactly as the R_S variance above does, and give:
    //
    //   Var(R_peak - R_S) = (1/T - 1/S)^2 * totNumErrSq
    //                     + [ errBkgNum^2 + (R_peak-R_S)^2 Var(T) + (R_B-R_S)^2 Var(B) ] / S^2
    //
    //   Var(R_peak - R_B) = totNumErrSq/T^2 + errBkgNum^2/B^2      (exactly Var(R_peak)+Var(R_B))
    //
    //   Var(R_S - R_B)    = [ totNumErrSq + (R_peak-R_S)^2 Var(T) + (R_S-R_B)^2 Var(B) ] / S^2
    //                     + (1/S + 1/B)^2 * errBkgNum^2
    //
    // with T = totCounts, B = bkgCounts, S = sigCounts. Two sanity checks are built into these:
    //   - The first is the R_S variance with totNumErrSq/S^2 replaced by (B/(T*S))^2*totNumErrSq,
    //     i.e. suppressed by f_B^2. The peak numerator fluctuation largely cancels between R_peak
    //     and R_S, which it must: they are built from the same candidates.
    //   - The second comes out as pure quadrature, because R_peak and R_B are measured in disjoint
    //     mass regions. That is the one pair for which the naive treatment happens to be right.
    // Same binding assumption as everywhere above: first-order propagation, unreliable once
    // sigCounts approaches zero.
    const double T = y.totCounts, B = y.bkgCounts, S = y.sigCounts;
    const double varT = y.totCountsErrSq, varB = y.errBkgCounts * y.errBkgCounts;
    const double varNumPeak = y.totNumErrSq, varNumBkg = y.errBkgNum * y.errBkgNum;
    const double sharedYieldTerm = std::pow(y.R_peak - y.R_S, 2) * varT;

    y.diffPeakMinusSig = y.R_peak - y.R_S;
    y.errDiffPeakMinusSig = std::sqrt(std::pow(1.0 / T - 1.0 / S, 2) * varNumPeak
                                      + (varNumBkg + sharedYieldTerm
                                         + std::pow(y.R_B - y.R_S, 2) * varB) / (S * S));

    y.diffPeakMinusBkg = y.R_peak - y.R_B;
    y.errDiffPeakMinusBkg = std::sqrt(varNumPeak / (T * T) + varNumBkg / (B * B));

    y.diffSigMinusBkg = y.R_S - y.R_B;
    y.errDiffSigMinusBkg = std::sqrt((varNumPeak + sharedYieldTerm
                                      + std::pow(y.R_S - y.R_B, 2) * varB) / (S * S)
                                     + std::pow(1.0 / S + 1.0 / B, 2) * varNumBkg);

    y.valid = true;
}

PeakWindowYields ComputePeakWindowYields(TH1D* hMassCounts, TH1D* hNumSum,
                                         double mu, double sigma, double nSigmaPeak,
                                         TF1* bkgFit, TFitResultPtr rBkg,
                                         TF1* ringBkgFit, TFitResultPtr rRingBkg)
{
    PeakWindowYields y;
    if (!hMassCounts || !bkgFit) return y;

    y.hasNumerator = (hNumSum != nullptr && ringBkgFit != nullptr);

    // --- Signal window, snapped to bin edges ---
    y.firstBin = hMassCounts->FindBin(mu - nSigmaPeak * sigma);
    y.lastBin = hMassCounts->FindBin(mu + nSigmaPeak * sigma);
    y.xLow = hMassCounts->GetBinLowEdge(y.firstBin);
    y.xHigh = hMassCounts->GetBinLowEdge(y.lastBin) + hMassCounts->GetBinWidth(y.lastBin);
    if (sigma > 0.0) {
        y.nSigmaAchievedLow = (mu - y.xLow) / sigma;
        y.nSigmaAchievedHigh = (y.xHigh - mu) / sigma;
    }
    y.coverage = GaussianWindowCoverage(y.xLow, y.xHigh, mu, sigma);

    // --- Raw sums over the window ---
    // Bins inside the peak are independent, so their errors add in quadrature.
    for (int jBin = y.firstBin; jBin <= y.lastBin; ++jBin) {
        y.totCounts += hMassCounts->GetBinContent(jBin);
        y.totCountsErrSq += std::pow(hMassCounts->GetBinError(jBin), 2);
        if (y.hasNumerator) {
            y.totNum += hNumSum->GetBinContent(jBin); // Accumulate Sum_R_i
            // hNumSum->GetBinError is sigma_R * sqrt(N_bin), the correct error on the sum
            y.totNumErrSq += std::pow(hNumSum->GetBinError(jBin), 2);
        }
    }

    if (y.totCounts <= 0) { y.failStage = 1; return y; }

    /*
    * MEDIUM COMMENT: BACKGROUND INTEGRATION AND COVARIANCE
    * ----------------------------------------------------------------------------------
    * The background must be integrated in one step using the full covariance matrix
    * because polynomial parameters are highly correlated.
    * Performing per-bin integrations and summing the errors in quadrature incorrectly
    * assumes independent uncertainties between bins. This artificial inflation of the
    * background error severely overestimates the final signal uncertainty.
    * TF1::IntegralError correctly uses the Jacobian of the integral with respect to
    * the parameters and the full parameter covariance matrix.
    *
    * The pol2 was fitted to a DENSITY (counts per unit mass), so TF1::Integral over the window
    * already returns counts -- exactly as integrating dN/dpT over a pT window returns counts.
    * No division by bin width is needed anywhere below.
    */
    y.bkgCounts = bkgFit->Integral(y.xLow, y.xHigh);
    y.errBkgCounts = bkgFit->IntegralError(y.xLow, y.xHigh, rBkg->GetParams(),
                                           rBkg->GetCovarianceMatrix().GetMatrixArray());

    if (y.hasNumerator) {
        // --- Background numerator, from the fitted <R>_bkg(m) ---------------------------------
        // M = Integral over the window of R(m) * b(m) dm, where b is the counts density fitted
        // above. Writing it as the b-WEIGHTED AVERAGE of R over the window,
        //     R_B_eff = Integral R(m) b(m) dm / Integral b(m) dm  =  M / B,
        // makes M = R_B_eff * B by construction, which is exactly the relation the covariance
        // structure downstream already assumes (Cov(M,B) = R_B * Var(B)). It also makes the
        // conditional error trivial to get right: R_B_eff is a LINEAR functional of the fitted
        // coefficients,
        //     R_B_eff = sum_k p_k * w_k,   w_k = Integral (m-shift)^k b(m) dm / B,
        // so its variance is w^T Cov(p) w exactly, with no gradient integration and therefore none
        // of the roundoff failures that TF1::IntegralError was reporting.
        const int nPar = ringBkgFit->GetNpar();

        // Weights, by Simpson quadrature on a smooth polynomial: cheap, and exact enough that the
        // quadrature error is orders of magnitude below the fit error.
        const int nSteps = 400; // even, for Simpson
        const double h = (y.xHigh - y.xLow) / nSteps;
        std::vector<double> w(nPar, 0.0);
        double norm = 0.0;
        for (int i = 0; i <= nSteps; ++i) {
            const double m = y.xLow + i * h;
            const double coef = (i == 0 || i == nSteps) ? 1.0 : ((i % 2) ? 4.0 : 2.0);
            const double bVal = bkgFit->Eval(m);
            norm += coef * bVal;
            // The k-th basis function of ringBkgFit evaluated at m, obtained by unit-poking the
            // parameters. This keeps the weights correct whatever basis PolShifted used.
            for (int k = 0; k < nPar; ++k) {
                std::vector<double> unit(nPar, 0.0);
                unit[k] = 1.0;
                w[k] += coef * bVal * ringBkgFit->EvalPar(&m, unit.data());
            }
        }
        if (norm > 0.0) {
            for (int k = 0; k < nPar; ++k) w[k] /= norm;

            double rEff = 0.0;
            for (int k = 0; k < nPar; ++k) rEff += ringBkgFit->GetParameter(k) * w[k];

            double varREff = 0.0;
            const TMatrixDSym& cov = rRingBkg->GetCovarianceMatrix();
            for (int j = 0; j < nPar; ++j)
                for (int k = 0; k < nPar; ++k) varREff += w[j] * w[k] * cov(j, k);

            y.bkgNum = rEff * y.bkgCounts;
            // Conditional on bkgCounts, exactly as the propagation downstream expects.
            y.errBkgNum = y.bkgCounts * std::sqrt(std::max(varREff, 0.0));
        } else {
            y.hasNumerator = false;
        }
    }

    FinalizeDerivedQuantities(y);
    return y;
}

// ================================================================================================
// COMBINING ANGULAR BINS
// ================================================================================================
// Two genuinely different ways to reduce a set of per-bin extractions to one number, and the
// DIFFERENCE between them is itself the measurement this folder exists to make.
//
//   SIGNAL-WEIGHTED. Sum the primitives across bins and extract from the totals:
//       <R>_S = sum_i (N_i - M_i) / sum_i (T_i - B_i)
//   which is the signal-yield weighted mean, sum_i S_i R_i / sum_i S_i. This is what "the mean of R
//   over all candidates" means, and it is therefore the ACCEPTANCE-WEIGHTED answer: the weights S_i
//   are exactly the non-uniform occupancies that the azimuthal efficiency effect produces.
//
//   FLAT-ACCEPTANCE. The unweighted bin average, (1/n) sum_i R_i. Bins are equal width in the
//   angular variable and hold disjoint candidates, so this is a plain independent average and its
//   error really is quadrature over n. This is what would be measured with uniform acceptance.
//
// Their difference is the effect. Neither number alone shows it.
//
// WHY NOT PROJECT AND EXTRACT ONCE. Projecting the 2D histograms onto the mass axis and running a
// single extraction gives the signal-weighted answer too, but it fits ONE sideband polynomial to
// the angle-summed numerator. Where <R> changes sign across the angular variable -- which is the
// entire reason for splitting on phi_Lambda - phi_p* -- that sum is a mixture of opposite-sign
// contributions with different slopes, and it is under no obligation to be linear. Splitting first
// and combining afterwards keeps every fit in the regime where a low-order sideband model is
// defensible. Combining is also legitimate: different angular bins hold disjoint candidates, so
// their primitives are independent and simply add.
// ================================================================================================

/// @brief Signal-weighted combination: sums the primitives, then reuses the standard algebra.
/// @note Independent bins, so variances add. Passing the totals through FinalizeDerivedQuantities
///       means the combined errors come from the same code as the per-bin ones, correlations and
///       all, instead of a parallel derivation that could quietly diverge from it.
PeakWindowYields CombineYieldsSignalWeighted(const std::vector<PeakWindowYields>& perBin)
{
    PeakWindowYields agg;
    if (perBin.empty()) return agg;

    agg.hasNumerator = perBin.front().hasNumerator;
    for (const auto& y : perBin) {
        if (!y.valid) continue;
        agg.totCounts += y.totCounts;
        agg.totCountsErrSq += y.totCountsErrSq;
        agg.totNum += y.totNum;
        agg.totNumErrSq += y.totNumErrSq;
        agg.bkgCounts += y.bkgCounts;
        agg.errBkgCounts += y.errBkgCounts * y.errBkgCounts; // Accumulate variance, root it below
        agg.bkgNum += y.bkgNum;
        agg.errBkgNum += y.errBkgNum * y.errBkgNum;
        agg.hasNumerator = agg.hasNumerator && y.hasNumerator;
    }
    if (agg.totCounts <= 0.0) return agg;

    agg.errBkgCounts = std::sqrt(agg.errBkgCounts);
    agg.errBkgNum = std::sqrt(agg.errBkgNum);

    FinalizeDerivedQuantities(agg);
    return agg;
}

/// @brief One flat-acceptance average and its error.
struct FlatAverage {
    double value = 0.0, error = 0.0;
    int    nBins = 0;
    bool   valid = false;
};

/// @brief Unweighted mean over independent angular bins, error = sqrt(sum sigma_i^2) / n.
FlatAverage AverageOverBins(const std::vector<double>& vals, const std::vector<double>& errs)
{
    FlatAverage out;
    if (vals.empty() || vals.size() != errs.size()) return out;
    double sum = 0.0, varSum = 0.0;
    for (size_t i = 0; i < vals.size(); ++i) { sum += vals[i]; varSum += errs[i] * errs[i]; }
    out.nBins = static_cast<int>(vals.size());
    out.value = sum / out.nBins;
    out.error = std::sqrt(varSum) / out.nBins;
    out.valid = true;
    return out;
}

/// @brief Signal-weighted minus flat-acceptance, with the shared candidates accounted for.
///
/// Both are linear combinations of the SAME independent per-bin values R_i, so quadrature between
/// them is wrong -- they are strongly positively correlated and it would badly overstate the error.
/// Writing D = sum_i (w_i - 1/n) R_i with w_i = S_i / sum_j S_j, the R_i are independent across
/// bins and the propagation collapses to
///     Var(D) = sum_i (w_i - 1/n)^2 * Var(R_i).
///
/// @note The weights are treated as fixed. They are ratios of signal yields, determined far more
///       precisely than the R_i themselves, and the neglected term is suppressed by
///       (R_i - R_weighted) / sum_j S_j. Worth revisiting only if the bin yields become comparable
///       to their own uncertainties.
/// @param vals   Per-bin R_i.
/// @param errs   Per-bin sigma(R_i).
/// @param sigCts Per-bin signal yields S_i, used for the weights.
FlatAverage WeightedMinusFlat(const std::vector<double>& vals, const std::vector<double>& errs,
                              const std::vector<double>& sigCts)
{
    FlatAverage out;
    const size_t n = vals.size();
    if (n == 0 || errs.size() != n || sigCts.size() != n) return out;

    double sumS = 0.0;
    for (double s : sigCts) sumS += s;
    if (sumS <= 0.0) return out;

    double diff = 0.0, var = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const double coeff = sigCts[i] / sumS - 1.0 / static_cast<double>(n);
        diff += coeff * vals[i];
        var += coeff * coeff * errs[i] * errs[i];
    }
    out.value = diff;
    out.error = std::sqrt(var);
    out.nBins = static_cast<int>(n);
    out.valid = true;
    return out;
}

// ================================================================================================
// ================================================================================================
// Helper: BuildNumFromProfile
// ================================================================================================
// PURPOSE:
//   Takes a TProfile2D (p2d, axes: angle x mass) and the corresponding TH2D counts histogram,
//   and constructs a corrected TH2D "hNum" where:
//     bin content = profile_mean * N_counts   (= Sum R_i, the absolute sum)
//     bin error   = profile_error * N_counts  (= sigma_R * sqrt(N), the correct error on the sum)
//
// ARGUMENTS:
//   p2d          -- the TProfile2D filled as Fill(angle, mass, R_i)
//   h2dCounts    -- the corresponding TH2D filled with unweighted counts (denominator)
//   histoName    -- a unique TString name for the output TH2D
//
// RETURNS:
//   A new TH2D* (caller owns it -- delete when done). Returns nullptr on bad input.
// NOTE: The binning of p2d and h2dCounts must match exactly.
// ================================================================================================
TH2D* BuildNumFromProfile(TProfile2D* p2d, TH2D* h2dCounts, TString histoName)
{
    // --- Input validation ---
    if (!p2d || !h2dCounts) {
        std::cerr << "[BuildNumFromProfile] ERROR: null input pointer.\n";
        return nullptr;
    }

    int nX = p2d->GetNbinsX();
    int nY = p2d->GetNbinsY();

    if (nX != h2dCounts->GetNbinsX() || nY != h2dCounts->GetNbinsY()) {
        std::cerr << "[BuildNumFromProfile] ERROR: bin mismatch between TProfile2D and TH2D.\n";
        return nullptr;
    }

    // --- Create the output TH2D with identical axis layout ---
    TH2D* hNum = (TH2D*)h2dCounts->Clone(histoName);
    hNum->SetDirectory(nullptr); // Sole ownership is the caller's.
                                 // Without this, Clone() also registers hNum in current
                                 // gDirectory, so the later "delete" by the caller would
                                 // leave a remaining entry in that directory's in-memory list.
    hNum->Reset(); // Clear all content and errors
    if (hNum->GetSumw2N() == 0) hNum->Sumw2(); // Ensure Sumw2 is on so errors are stored correctly

    // --- Fill bin by bin ---
    for (int ix = 1; ix <= nX; ++ix) {
        for (int iy = 1; iy <= nY; ++iy) {

            // Profile mean = <R> in this (angle, mass) cell
            double mean  = p2d->GetBinContent(ix, iy);

            // Profile error = sigma_R / sqrt(N) = standard error of the mean
            double seom  = p2d->GetBinError(ix, iy);

            // Number of Lambda candidates (unweighted count) in this cell
            double nCand = h2dCounts->GetBinContent(ix, iy);

            if (nCand <= 0) continue;  // Empty bin: leave at zero

            // Sum = mean * N
            double sumR     = mean * nCand;

            // Error on sum = (sigma_R / sqrt(N)) * N = sigma_R * sqrt(N)
            double errSumR  = seom * nCand;

            hNum->SetBinContent(ix, iy, sumR);
            hNum->SetBinError(ix, iy, errSumR);
        }
    }

    return hNum;
}

// ================================================================================================
// Helper:  BuildNumFromProfile3D
// ================================================================================================
// PURPOSE:
//   Constructs a corrected TH3D from a TProfile3D + counts TH3D such that:
//     bin content = Sum_R_i  (profile mean * N counts in that cell)
//     bin error   = sigma_R * sqrt(N)  (correct error on the sum)
//
//   The output TH3D is a drop-in replacement for the original h3dRingObservable* histograms.
//   It can be projected with Project3D("... e") exactly as before, and the projected
//   TH2D/TH1D will carry correct errors automatically.
//
// ARGUMENTS:
//   p3d       -- TProfile3D filled as Fill(x, y, z, R_i).
//                Axes must match h3dCounts exactly.
//   h3dCounts -- TH3D with unweighted Lambda candidate counts (the denominator histogram).
//                Essentially just a reference for the new format!
//   histoName -- unique name for the output TH3D.
//
// RETURNS:
//   A new TH3D* owned by the caller. Returns nullptr on bad input.
//
// NOTE ON PROJECT3D + ERROR PROPAGATION:
//   When TH3D::Project3D("yx e") is called after SetRange on the Z axis, ROOT sums bin
//   contents over Z and adds errors in quadrature. This is exactly what we want: the
//   projected bin content is Sum_{iz} Sum_R_i(ix,iy,iz) and the projected bin error is
//   sqrt(Sum_{iz} sigma_R_iz^2 * N_iz). Both are the correct statistics for the
//   angle-mass slice in that kinematic window.
// ================================================================================================
TH3D* BuildNumFromProfile3D(TProfile3D* p3d, TH3D* h3dCounts, TString histoName)
{
    // --- Input validation ---
    if (!p3d || !h3dCounts) {
        std::cerr << "[BuildNumFromProfile3D] ERROR: null input pointer.\n";
        return nullptr;
    }

    int nX = p3d->GetNbinsX();
    int nY = p3d->GetNbinsY();
    int nZ = p3d->GetNbinsZ();

    // Axis consistency check -- binning must match between profile and counts
    if (nX != h3dCounts->GetNbinsX() ||
        nY != h3dCounts->GetNbinsY() ||
        nZ != h3dCounts->GetNbinsZ()) {
        std::cerr << "[BuildNumFromProfile3D] ERROR: bin mismatch between TProfile3D and TH3D.\n";
        std::cerr << "  Profile: (" << nX << ", " << nY << ", " << nZ << ")\n";
        std::cerr << "  Counts:  (" << h3dCounts->GetNbinsX() << ", "
                                    << h3dCounts->GetNbinsY() << ", "
                                    << h3dCounts->GetNbinsZ() << ")\n";
        return nullptr;
    }

    // --- Create output TH3D with identical axis layout ---
    // Clone from h3dCounts to inherit axis labels, titles, and bin edges exactly.
    TH3D* hOut = (TH3D*)h3dCounts->Clone(histoName);
    hOut->SetDirectory(nullptr); // Sole ownership is the caller's. Same as the TProfile2D built by ConvertToProfile2D.
    hOut->Reset();    // Clear content; keep axis structure
    if (hOut->GetSumw2N() == 0) hOut->Sumw2();

    // --- Fill bin by bin ---
    for (int ix = 1; ix <= nX; ++ix) {
        for (int iy = 1; iy <= nY; ++iy) {
            for (int iz = 1; iz <= nZ; ++iz) {

                // Profile mean in this cell: <R> = Sum_R_i / N
                double mean = p3d->GetBinContent(ix, iy, iz);

                // Profile standard error of the mean: sigma_R / sqrt(N)
                // (ROOT TProfile default error option "" gives the standard error of the mean)
                double seom = p3d->GetBinError(ix, iy, iz);

                // Unweighted Lambda count in this cell
                double nCand = h3dCounts->GetBinContent(ix, iy, iz);

                // Skip empty cells: leave content and error at zero
                if (nCand <= 0) continue;

                // Reconstruct the absolute sum: Sum_R_i = <R> * N
                double sumR = mean * nCand;

                // Reconstruct the correct error on the sum:
                //   sigma(Sum_R_i) = (sigma_R / sqrt(N)) * N = sigma_R * sqrt(N)
                double errSumR = seom * nCand;

                hOut->SetBinContent(ix, iy, iz, sumR);
                hOut->SetBinError(ix, iy, iz, errSumR);
            }
        }
    }

    return hOut;
}

// ================================================================================================
// Helper: ConvertToProfile2D
// ================================================================================================
// PURPOSE:
//   Converts a corrected error TH2D (content = Sum_R_i, error = sigma_R * sqrt(N)) and
//   its matching counts TH2D into a TProfile2D (content = <R>, error = sigma_R / sqrt(N)).
//
//   The resulting TProfile2D is a valid input to BuildNumFromProfile inside
//   ExtractObservable2D, completing the DRY (don't repeat yourself) round-trip without a 
//   second overload.
//
// ARGUMENTS:
//   h2dNumCorr  -- corrected TH2D: content = Sum_R_i, error = sigma_R * sqrt(N).
//                  Built by BuildNumFromProfile or by projecting a corrected TH3D.
//   h2dCounts   -- matching counts TH2D: content = N (unweighted). Axes must match h2dNumCorr.
//   profileName -- unique name for the output TProfile2D.
//
// RETURNS:
//   A new TProfile2D* owned by the caller. Returns nullptr on bad input or bin mismatch.
//
// USAGE (inside Step 7 of the main macro):
//   TH2D* h2dCountsPhi   = (TH2D*)h3dDeltaPhiVsMassVsLambdaPt->Project3D("yx e" ...);
//   TH2D* h2dNumCorrPhi  = (TH2D*)h3dNumCorrDeltaPhiVsMassVsLambdaPt->Project3D("yx e" ...);
//   TProfile2D* p2dPhi   = ConvertToProfile2D(h2dNumCorrPhi, h2dCountsPhi, "p2d_phi_name");
//   ExtractObservable2D(h2dCountsPhi, p2dPhi, ...);   // single function, no overload needed
//   delete p2dPhi;        // caller owns it
//   // h2dCountsPhi and h2dNumCorrPhi are owned by ROOT (Project3D output); do NOT delete.
// ================================================================================================
TProfile2D* ConvertToProfile2D(TH2D* h2dNumCorr, TH2D* h2dCounts, TString profileName)
{
    // --- Input validation ---
    if (!h2dNumCorr || !h2dCounts) {
        std::cerr << "[ConvertToProfile2D] ERROR: null input pointer.\n";
        return nullptr;
    }

    int nX = h2dNumCorr->GetNbinsX();
    int nY = h2dNumCorr->GetNbinsY();

    if (nX != h2dCounts->GetNbinsX() || nY != h2dCounts->GetNbinsY()) {
        std::cerr << "[ConvertToProfile2D] ERROR: bin mismatch between h2dNumCorr and h2dCounts.\n";
        std::cerr << "  NumCorr: (" << nX << ", " << nY << ")\n";
        std::cerr << "  Counts:  (" << h2dCounts->GetNbinsX() << ", "
                                    << h2dCounts->GetNbinsY() << ")\n";
        return nullptr;
    }

    // --- Create TProfile2D with matching axis layout ---
    // We use the counts histogram axis edges to initialise the profile, since we cannot
    // pass a variable-bin-edge array to TProfile2D's standard constructor directly.
    // The cleanest approach: build from the TAxis objects of h2dCounts.
    TAxis* xAx = h2dCounts->GetXaxis();
    TAxis* yAx = h2dCounts->GetYaxis();

    // Fallback for uniform-bin axes (GetXbins()->GetArray() returns nullptr if uniform)
    // For each axis: use bin-edge array if variable, xmin/xmax if uniform
    const double* xEdges = xAx->GetXbins()->GetArray();
    const double* yEdges = yAx->GetXbins()->GetArray(); // This is the mass axis, usually, so there is no problem with the X axis being uniform
    TProfile2D* prof = nullptr;
    if (xEdges && yEdges) {
        prof = new TProfile2D(profileName, Form("%s;%s;%s;<R>", profileName.Data(),xAx->GetTitle(), yAx->GetTitle()), nX, xEdges, nY, yEdges);
    } else if (xEdges) {
        prof = new TProfile2D(profileName, Form("%s;%s;%s;<R>", profileName.Data(),xAx->GetTitle(), yAx->GetTitle()), nX, xEdges, nY, yAx->GetXmin(), yAx->GetXmax());
    } else if (yEdges) {
        prof = new TProfile2D(profileName, Form("%s;%s;%s;<R>", profileName.Data(),xAx->GetTitle(), yAx->GetTitle()), nX, xAx->GetXmin(), xAx->GetXmax(), nY, yEdges);
    } else {
        prof = new TProfile2D(profileName, Form("%s;%s;%s;<R>", profileName.Data(),xAx->GetTitle(), yAx->GetTitle()), nX, xAx->GetXmin(), xAx->GetXmax(), nY, yAx->GetXmin(), yAx->GetXmax());
    }

    prof->SetDirectory(nullptr); // To enable our own management of this TProfile in memory (easier in our weird application)

    // // Set the error option to "" (standard error of the mean), which is also the ROOT default.
    // // This is explicit documentation of intent; it matches the TProfile2D used upstream.
    // prof->SetErrorOption("");

    // --- Fill the TProfile2D bin by bin using SetBinContent / SetBinEntries ---
    // We bypass the Fill(x, y, val) interface entirely because we do not have the
    // individual R_i values -- only their sum and error. Instead we directly set the
    // internal TProfile2D accumulators:
    //
    //   TProfile2D internally stores (per bin):
    //     fArray[bin]    = Sum_W  = Sum_R_i  (for unit weights, = Sum_R_i)
    //     fSumw2[bin]    = Sum_W2 = Sum_R_i^2  (needed to compute the spread)
    //     fBinEntries[bin] = N   (the count)
    //   From our inputs:
    //     Sum_R_i  = h2dNumCorr->GetBinContent(ix, iy)
    //     N        = h2dCounts->GetBinContent(ix, iy)
    //     errSumR  = h2dNumCorr->GetBinError(ix, iy)  = sigma_R * sqrt(N)
    //   We need Sum_R_i^2 to reproduce the correct standard error of the mean:
    //     seom = sqrt( (Sum_R_i^2/N - (Sum_R_i/N)^2) / N )
    //          = sqrt( Var(R) / N )
    //          = sigma_R / sqrt(N)
    //   Since errSumR = sigma_R * sqrt(N), we have sigma_R = errSumR / sqrt(N), so:
    //     Var(R) = sigma_R^2 = errSumR^2 / N
    //   And:
    //     Sum_R_i^2 = N * (Var(R) + <R>^2)
    //               = N * (errSumR^2/N + (Sum_R_i/N)^2)
    //               = errSumR^2 + Sum_R_i^2 / N
    //   Which gives:
    //     Sum_R_i^2 = N * errSumR^2 / (N - 1)   [Bessel-corrected]
    //   or for large N:
    //     Sum_R_i^2 ~ errSumR^2 + (Sum_R_i)^2 / N
    //   We use the exact relation:
    //     Sum_R_i^2 = errSumR^2 * N + Sum_R_i^2 / N
    //               = N * (errSumR^2 + mean^2)
    //   since Sum_W2 = N * (Var + mean^2) = N * E[R^2].
    for (int ix = 1; ix <= nX; ++ix) {
        for (int iy = 1; iy <= nY; ++iy) {
            double sumR    = h2dNumCorr->GetBinContent(ix, iy);  // Sum_R_i
            double errSumR = h2dNumCorr->GetBinError(ix, iy);    // sigma_R * sqrt(N)
            double nCand   = h2dCounts->GetBinContent(ix, iy);   // N

            // Skip empty or undefined cells
            if (nCand <= 0) continue;

            double mean  = sumR / nCand;                          // <R>
            double seom  = errSumR / nCand;                       // sigma_R / N
            // sigma_R^2 = seom^2 * N  (variance of individual R values)
            double varR  = seom * seom * nCand;

            // Sum_R_i^2 = N * (Var(R) + <R>^2) = N * E[R^2]
            double sumR2 = nCand * (varR + mean * mean);

            // Global bin index (TProfile2D uses the same global bin as TH2)
            int gBin = prof->GetBin(ix, iy);

            // Directly set the internal TProfile2D accumulators
            prof->SetBinContent(gBin, sumR);            // fArray[gBin] = Sum_W = Sum_R_i
            prof->SetBinError(gBin, seom);              // sets fSumw2 correctly via ROOT internals
            prof->SetBinEntries(gBin, nCand);           // fBinEntries[gBin] = N

            // SetBinError on a TProfile2D sets fSumw2[gBin] = seom^2 * nCand^2 (ROOT convention).
            // But we need fSumw2[gBin] = Sum_R_i^2 for the variance formula.
            // ROOT's TProfile::GetBinError() computes:
            //   seom = sqrt(max(0, fSumw2/fBinEntries - (fArray/fBinEntries)^2) / fBinEntries)
            //        = sqrt(max(0, Sum_R_i^2/N - <R>^2) / N)
            // So we must set fSumw2 = sumR2 directly via the Sumw2 array:
            prof->GetSumw2()->SetAt(sumR2, gBin);
        }
    }
    return prof;
}

// LEGACY CODE! TF1 actually performs the exact \int f(m) * dm integration we need,
// which transforms the f(m) = dN/dm (a density of counts as function of mass, properly
// normalized to the bin size of each bin to be an actual density) into an actual number
// of counts after integrating! We DO NOT need this code!
// // NOTE: the polynomial is HARD CODED to be of order 2! Careful when changing this!
// std::pair<double,double> IntegratePolynomialOverBins(
//     TF1* func, TFitResultPtr& fitResult,
//     TH1D* histo, int firstBin, int lastBin)
// {
//     // Build the Jacobian vector d_k = sum_j (m_j^k * delta_m_j)
//     double d[3] = {0.0, 0.0, 0.0};
//     double central = 0.0;

//     for (int jBin = firstBin; jBin <= lastBin; ++jBin) {
//         double m   = histo->GetBinCenter(jBin);
//         double dm  = histo->GetBinWidth(jBin);
//         double fval = func->Eval(m);      // f(m) = c0 + c1*m + c2*m^2 [density]

//         central += fval * dm;             // Riemann sum for count

//         // Jacobian components: d(N_B)/d(c_k) = sum_j m_j^k * delta_m_j
//         d[0] += dm;
//         d[1] += m * dm;
//         d[2] += m * m * dm;
//     }

//     // Propagate through covariance matrix: sigma^2 = d^T V d
//     TMatrixDSym cov = fitResult->GetCovarianceMatrix();
//     double var = 0.0;
//     for (int i = 0; i < 3; ++i)
//         for (int j = 0; j < 3; ++j)
//             var += d[i] * cov(i, j) * d[j];

//     return {central, std::sqrt(var)};
// }
// ================================================================================================
// ================================================================================================

// QA additional plot for signal extraction using Simultaneous Fit strategy:
// (uses previously obtained values as initial guesses from the regular signal extraction strategy)
    // Struct to hold the QA fit results:
struct SimFitResult {
    double R_S, err_R_S;
    double R_B, err_R_B;
    double SigYield, err_SigYield;
    double BkgYield, err_BkgYield;
    double Purity, err_Purity;
    double Significance, err_Significance;
    int status; 
};
    // The standalone simultaneous fit function:
// (TODO: implement a more sophisticated method on PerformSimultaneousFitQA that actually does bin counting in the peak, whilst also using a covariance-matrix like entity for full error propagation)
/**
 * @brief Performs a simultaneous Chi2 fit on mass-density and numerator-density spectra
 * to extract signal/background yields and ring observable values (R_S, R_B).
 *
 * IMPORTANT INPUT EXPECTATIONS:
 * -----------------------------
 * @param hMassDensity   Mass spectrum scaled by bin width (Y-axis = dN/dM, counts/GeV).
 *                       Obtained via hMassRaw->Clone(...); hMassDensity->Scale(1.0,"width").
 *                       Variable binning is natively supported.
 * @param hNumDensity    Numerator spectrum scaled by bin width (Y-axis = Sum_R_i / GeV).
 *                       Must have the exact same binning structure as hMassDensity.
 *                       Obtained the same way from the raw hNum histogram.
 * @param massMin        Lower bound of the mass fit range.
 * @param massMax        Upper bound of the mass fit range.
 * @param init_sigYield  Initial guess for the total signal yield (absolute counts,
 *                       i.e. the Gaussian area = amplitude * sigma * sqrt(2*pi)).
 * @param init_mu        Initial guess for the Gaussian mean.
 * @param init_sigma     Initial guess for the Gaussian width.
 * @param bkgFitFunc     A TF1 (pol2) previously fitted on a MASS DENSITY histogram.
 *                       par[0..2] = c0, c1, c2 in counts/GeV units.
 *                       Now consistently required to be a density-space pol2 for
 *                       both the per-bin calls (from grBkg sideband fit) and the
 *                       integrated call (from preFit par[3..5] transplanted into
 *                       a fresh pol2). No unit conversion is needed at the call site.
 * @param init_RS        Initial guess for the signal ring observable value.
 * @param init_RB        Initial guess for the background ring observable value.
 *
 * INTERNAL MECHANICS:
 * -------------------
 * - ALL histograms are in density convention (counts/GeV) on entry. The chi2
 *   compares density model vs density data directly, with no delta_m_i weighting
 *   of data. Signal and background yields are recovered by analytical integration
 *   of the fitted density functions over the signal window.
 * - Fits are performed in "normalized density" space (divided by hMassDensity maximum)
 *   to keep Minuit parameters O(1) and ensure convergence. Yields are recovered by
 *   integrating the normalized density and rescaling by normFactor * signal_window_width.
 * - Background yields within the +/- 4 sigma window are calculated via exact analytical
 *   integration, avoiding discrete bin-width estimation errors.
 * - The numerator density model is R_S * S_density(m) + R_B * B_density(m), compared
 *   bin-by-bin against hNumDensity. This is consistent because both hMassDensity and
 *   hNumDensity share the same bin-width scaling.
 * @return SimFitResult struct containing absolute yields, ratios, purity, significance,
 * and their analytically propagated errors using the Minuit covariance matrix.
 */
SimFitResult PerformSimultaneousFitQA(TH1D* hMassDensity, TH1D* hNumDensity, double massMin, double massMax, double init_sigYield, double init_mu, double init_sigma, TF1* bkgFitFunc, double init_RS, double init_RB, double nSigmaPeak) {
    SimFitResult result = { // Initialized as null result
        0.0, 0.0, // R_S, err_R_S
        0.0, 0.0, // R_B, err_R_B
        0.0, 0.0, // SigYield, err_SigYield
        0.0, 0.0, // BkgYield, err_BkgYield
        0.0, 0.0, // Purity, err_Purity
        0.0, 0.0, // Significance, err_Significance
        -1        // status (Default to failed)
    };
    if (!hMassDensity || !hNumDensity || !bkgFitFunc) return result;

    // All bin widths are fetched per-bin inside the loop where needed.
    // In the density convention the chi2 compares f(m) vs data directly --
    // no delta_m_i multiplication of the data point is required here.

    // --- COMBINED STABILITY CHECKS 1 & 1.5 (START) ---
    // Find the bins corresponding to the requested mass range.
    int binMin = hMassDensity->FindBin(massMin);
    int binMax = hMassDensity->FindBin(massMax);
    double sumMassEstimate = 0.0; // Actually counts estimator via bin width
    int nPointsMass = 0;
    // A single, efficient loop to compute both the integrals and non-empty bin counts
    for (int i = binMin; i <= binMax; ++i) {
        double valMassCounts = hMassDensity->GetBinContent(i) * hMassDensity->GetBinWidth(i);
        sumMassEstimate += valMassCounts;
        if (valMassCounts > 0.0) nPointsMass++; // Count non-empty bins to ensure we have degrees of freedom
    }
    // CHECK 1: Enough non-empty bins for a simultaneous fit with 8 parameters?
    if (nPointsMass < 8) return result; // Returns the default struct with status = -1
    // CHECK 1.5: Are the histograms basically empty?
    if (sumMassEstimate <= 0.0) return result;
    // --- COMBINED STABILITY CHECKS 1 & 1.5 (END!) ---

    // --- SCALING FACTOR ---
    // We normalize the mass density spectrum so the peak density is 1.0 (counts/GeV / peak_density).
    // This keeps Minuit parameters O(1) regardless of the absolute luminosity,
    // fixing the "Status 1" convergence issues that arose from O(1e9) density values.
    // MAJOR CHANGE: normFactor is now the peak of the DENSITY histogram, not the raw counts peak.
    // All internal model quantities are in normalized-density space throughout.
    double normFactor = hMassDensity->GetMaximum();
    if (normFactor <= 0) normFactor = 1.0;

    // 1. Define the Joint Chi2
    // MAJOR CHANGE: chi2 now compares density model vs density data directly.
    // The mass component: model = (S_density_norm + B_density_norm), data = hMassDensity / normFactor.
    // The numerator component: model = R_S * S_density_abs + R_B * B_density_abs,
    //   data = hNumDensity (already in density units, not rescaled further).
    // Both components are bin-point comparisons of density values -- no delta_m_i
    // multiplication of data is needed. delta_m_i only enters when converting
    // a density integral to counts (done analytically in the result extraction below).
    auto globalChi2 = [&](const double *par) {
        double chi2 = 0;
        for (int i = 1; i <= hMassDensity->GetNbinsX(); ++i) {
            double m = hMassDensity->GetBinCenter(i);
            if (m < massMin || m > massMax) continue;

            // Get density data and errors
            double y_massDensity      = hMassDensity->GetBinContent(i);
            double e_massDensity      = hMassDensity->GetBinError(i);
            double y_numDensity       = hNumDensity->GetBinContent(i);
            double e_numDensity       = hNumDensity->GetBinError(i);

            if (e_massDensity <= 0) e_massDensity = 1.0;
            if (e_numDensity  <= 0) e_numDensity  = 1.0;

            // --- MASS DENSITY COMPONENT (NORMALIZED) ---
            // par[0] = Normalized signal amplitude (density at peak / normFactor, approx 1.0).
            // par[1] = mu, par[2] = sigma.
            // par[3,4,5] = Normalized background density coefficients (counts/GeV / normFactor).
            //
            // S_density_norm(m) = par[0] * Gaus(m; mu, sigma)  [normalized Gaussian = 1/GeV]
            // B_density_norm(m) = par[3] + par[4]*m + par[5]*m^2
            // Both are already densities -- no delta_m_i needed for the chi2 comparison.
            double S_density_norm = par[0] * TMath::Gaus(m, par[1], par[2], true);
            double B_density_norm = par[3] + par[4]*m + par[5]*m*m;

            double M_density_model_norm = S_density_norm + B_density_norm;

            // Normalize data to match model scale
            double y_massDensity_norm = y_massDensity / normFactor;
            double e_massDensity_norm = e_massDensity / normFactor;

            chi2 += std::pow((y_massDensity_norm - M_density_model_norm) / e_massDensity_norm, 2);

            // --- NUMERATOR DENSITY COMPONENT (ABSOLUTE) ---
            // The numerator density model is: R_S * S_density_abs(m) + R_B * B_density_abs(m)
            // where S_density_abs and B_density_abs are the unscaled (absolute) density values.
            // We recover absolute densities by multiplying the normalized model by normFactor.
            // This is consistent because hNumDensity was NOT divided by normFactor.
            double S_density_abs = S_density_norm * normFactor;
            double B_density_abs = B_density_norm * normFactor;

            // par[6] = R_S, par[7] = R_B
            double Num_density_model = par[6] * S_density_abs + par[7] * B_density_abs;

            chi2 += std::pow((y_numDensity - Num_density_model) / e_numDensity, 2);
        }
        return chi2;
    };

    // 2. Setup Minimizer
    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    if (!min) return result; // Fallback if Minuit2 isn't loaded

    // Increased tolerance and calls slightly for stability with the new scaling
    min->SetMaxFunctionCalls(1000000);
    min->SetTolerance(0.01);
    min->SetPrintLevel(0); // Keep it quiet so it doesn't flood the terminal

    ROOT::Math::Functor f(globalChi2, 8);
    min->SetFunction(f);

    // 3. Initialize Parameters (NORMALIZED DENSITY SPACE)
    // MAJOR CHANGE: all initial guesses are now divided by normFactor (the density peak),
    // consistent with the normalized-density chi2 above.
    //
    // par[0]: Normalized signal amplitude = (init_sigYield / (sigma * sqrt(2*pi))) / normFactor
    //   init_sigYield is total counts; dividing by (sigma*sqrt(2*pi)) converts to peak amplitude
    //   in density units; dividing by normFactor puts it in normalized space (~1.0 at peak).
    double init_sigAmplitude_norm = (init_sigYield / (init_sigma * std::sqrt(2.0 * TMath::Pi()))) / normFactor;
    double step_S = std::max(0.01, init_sigAmplitude_norm * 0.05);
    min->SetVariable(0, "SigAmplitude_Norm", init_sigAmplitude_norm, step_S);
    min->SetVariable(1, "Mu",                init_mu,                0.0001);
    min->SetVariable(2, "Sigma",             init_sigma,             0.0001);

    // Use the actual parameter errors from the density-fitted bkgFitFunc as step sizes.
    // bkgFitFunc is a pol2 fitted to a density histogram, so par[0..2] are already
    // in counts/GeV units. Dividing by normFactor puts them in normalized-density space.
    double err_c0 = bkgFitFunc->GetParError(0) / normFactor;
    double err_c1 = bkgFitFunc->GetParError(1) / normFactor;
    double err_c2 = bkgFitFunc->GetParError(2) / normFactor;
    // (We keep a 0.1 fallback just in case the error matrix was empty)
    min->SetVariable(3, "Bkg_c0_Norm", bkgFitFunc->GetParameter(0) / normFactor, err_c0 > 0 ? err_c0 : 0.1);
    min->SetVariable(4, "Bkg_c1_Norm", bkgFitFunc->GetParameter(1) / normFactor, err_c1 > 0 ? err_c1 : 0.1);
    min->SetVariable(5, "Bkg_c2_Norm", bkgFitFunc->GetParameter(2) / normFactor, err_c2 > 0 ? err_c2 : 0.1);
    min->SetVariable(6, "R_S",         init_RS,                                  0.01);
    min->SetVariable(7, "R_B",         init_RB,                                  0.01);

    // Lock the kinematics tight to prevent the fit from wandering
    min->SetVariableLimits(1, init_mu - 0.003, init_mu + 0.003); // Slightly wider window
    min->SetVariableLimits(2, init_sigma * 0.5, init_sigma * 1.5);
    // Don't let the signal amplitude go negative
    min->SetVariableLowerLimit(0, 0.0);

    // 4. Minimize
    min->Minimize();
    result.status = min->Status();

    if (result.status == 0) {
        const double *fitVals = min->X();
        const double *fitErrs = min->Errors();

        // R_S and R_B are pure ratios -- scale-invariant and density-invariant.
        // They are extracted directly from the fit without any rescaling.
        result.R_S     = fitVals[6];
        result.err_R_S = fitErrs[6];
        result.R_B     = fitVals[7];
        result.err_R_B = fitErrs[7];

        // =====================================================================
        // Recover absolute signal yield from the normalized density amplitude.
        // MAJOR CHANGE: par[0] is now the normalized Gaussian amplitude (counts/GeV
        // divided by normFactor), not a normalized total count.
        // The total Gaussian area (= absolute signal yield in counts) is:
        //   SigYield = amplitude_abs * sigma * sqrt(2*pi)
        //            = (par[0] * normFactor) * par[2] * sqrt(2*pi)
        // This is the exact inverse of how init_sigAmplitude_norm was constructed.
        // =====================================================================
        double amplitude_abs   = fitVals[0] * normFactor; // density amplitude, counts/GeV
        double err_amplitude_abs = fitErrs[0] * normFactor;

        double muVal    = fitVals[1];
        double sigmaVal = fitVals[2];

        // Full Gaussian integral = amplitude * sigma * sqrt(2*pi), of which only the fraction
        // inside the signal window is counted. This used to be the constant 0.9999366, which is
        // erf(4/sqrt2) and therefore silently assumed a 4-sigma window; at 2 sigma the true
        // coverage is 0.9545 and the yield would have come out 4.7% high. The window here is
        // analytic (mu +/- nSigmaPeak*sigma, no bin snapping), so the symmetric form is exact.
        double xLow  = muVal - nSigmaPeak * sigmaVal;
        double xHigh = muVal + nSigmaPeak * sigmaVal;
        double windowCoverage = GaussianWindowCoverage(xLow, xHigh, muVal, sigmaVal);

        double raw_SigYield     = amplitude_abs * sigmaVal * std::sqrt(2.0 * TMath::Pi());
        double raw_errSigYield  = err_amplitude_abs * sigmaVal * std::sqrt(2.0 * TMath::Pi());
        result.SigYield    = raw_SigYield    * windowCoverage;
        result.err_SigYield = raw_errSigYield * windowCoverage;

        // =====================================================================
        // Recover absolute background yield by analytically integrating the
        // fitted background density over the same signal window used just above.
        // MAJOR CHANGE: the integration is now entirely in density space.
        // integral of (c0 + c1*m + c2*m^2) dm from xLow to xHigh gives counts
        // directly (density [counts/GeV] * mass interval [GeV] = counts).
        // We integrate the ABSOLUTE density (par_norm * normFactor), so the
        // final result is already in counts with no further rescaling needed
        // beyond the normFactor multiplication applied to c0_norm, c1_norm, c2_norm.
        // =====================================================================
        // Extract normalized background polynomial coefficients
        double c0_norm = fitVals[3];
        double c1_norm = fitVals[4];
        double c2_norm = fitVals[5];

        // Analytical integration terms: d(integral)/d(c_k) for a pol2
        // (polynomial integration is simple enough to do by hand)
        double dx1 = (xHigh - xLow);
        double dx2 = (xHigh * xHigh - xLow * xLow) / 2.0;
        double dx3 = (xHigh * xHigh * xHigh - xLow * xLow * xLow) / 3.0;

        // Integrate the normalized background density over the signal window.
        // Multiplying by normFactor converts from normalized-density counts to absolute counts.
        double integral_norm = (c0_norm * dx1) + (c1_norm * dx2) + (c2_norm * dx3);
        result.BkgYield = integral_norm * normFactor;

        // Analytically propagate the background yield error using the Minuit covariance matrix.
        // Var(integral_norm) = sum_{i,j} (dI/dc_i_norm)(dI/dc_j_norm) Cov(c_i_norm, c_j_norm)
        // then Var(BkgYield) = Var(integral_norm) * normFactor^2
        double dI_dc[3] = {dx1, dx2, dx3};
        double var_integral_norm = 0.0;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                // CovMatrix indices for c0, c1, c2 are 3, 4, 5
                double cov_ij = min->CovMatrix(3 + i, 3 + j);
                var_integral_norm += dI_dc[i] * dI_dc[j] * cov_ij;
            }
        }
        result.err_BkgYield = std::sqrt(var_integral_norm) * normFactor;

        // ==========================================================
        // Calculate Purity, Significance, and their propagated errors
        // using the full Minuit covariance matrix, including Cov(S, B).
        // ==========================================================
        // Cov(S_abs, B_abs) = normFactor^2 * sum_k [ dI/dc_k_norm * Cov(par[0], c_k_norm) ]
        // Note: par[0] is SigAmplitude_Norm; SigYield = par[0]*normFactor*sigmaVal*sqrt(2pi)*coverage,
        // so Var(SigYield) carries an extra factor of (sigmaVal*sqrt(2pi)*coverage)^2 relative
        // to Var(par[0]). For Cov(S,B) the same factor applies to the par[0] row/column.
        double gaussFactor = sigmaVal * std::sqrt(2.0 * TMath::Pi()) * windowCoverage;

        double cov_S_B_norm = 0.0;
        // --- Cov(SigAmplitude_Norm, integral_norm) in normalized parameter space ---
        for (int k = 0; k < 3; ++k) {
            double cov_0k = min->CovMatrix(0, 3 + k); // Cov(par[0], c_k) in normalized space
            cov_S_B_norm += dI_dc[k] * cov_0k;
        }
        // Scale to absolute counts space:
        // Cov(SigYield, BkgYield) = gaussFactor * normFactor * integral_factor * normFactor * cov_S_B_norm
        double cov_SigYield_BkgYield = gaussFactor * normFactor * normFactor * cov_S_B_norm;

        // --- Variances in absolute counts ---
        double Var_SigYield = min->CovMatrix(0, 0) * std::pow(gaussFactor * normFactor, 2);
        double Var_BkgYield = var_integral_norm * (normFactor * normFactor);

        // =======================
        // Purity and Significance
        // =======================
        double S = result.SigYield;
        double B = result.BkgYield;
        double N = S + B;
        if (N > 0) {
            // Central values
            result.Purity      = S / N;
            result.Significance = S / std::sqrt(N);

            // Propagated error for Purity: dP/dS = B/N^2, dP/dB = -S/N^2
            // Var(P) = (dP/dS)^2 Var(S) + (dP/dB)^2 Var(B) + 2(dP/dS)(dP/dB) Cov(S,B)
            double dP_dS = B / (N * N);
            double dP_dB = -S / (N * N);
            double VarP  = dP_dS*dP_dS * Var_SigYield
                         + dP_dB*dP_dB * Var_BkgYield
                         + 2.0 * dP_dS * dP_dB * cov_SigYield_BkgYield;
            result.err_Purity = (VarP > 0.0) ? std::sqrt(VarP) : 0.0;

            // Propagated error for Significance = S/sqrt(N):
            // dSig/dS = (B + S/2) / N^(3/2),  dSig/dB = -(S/2) / N^(3/2)
            double dSig_dS = (B + S / 2.0) / std::pow(N, 1.5);
            double dSig_dB = -(S / 2.0)    / std::pow(N, 1.5);
            double VarSig  = dSig_dS*dSig_dS * Var_SigYield
                           + dSig_dB*dSig_dB * Var_BkgYield
                           + 2.0 * dSig_dS * dSig_dB * cov_SigYield_BkgYield;
            result.err_Significance = (VarSig > 0.0) ? std::sqrt(VarSig) : 0.0;
        }
        else {
            result.Purity = 0; result.err_Purity = 0;
            result.Significance = 0; result.err_Significance = 0;
        }
    }

    delete min;
    return result;
}

// (TODO: make helper function numerically stable for full experimental statistics, i.e.,
// normalize particle counters in a smart way that does not break signal extraction)
// =================================================================================================
// HELPER FUNCTION: Full 2D to 1D Signal Extraction Engine
// =================================================================================================
void ExtractObservable2D(TH2D* h2dCounts, TProfile2D* p2dRingObs, TDirectory* parentDir,
                         TString extractionName, TString axisTitle, double massMin, double massMax,
                         const SidebandConfig& cfg, const SidebandConfig& cfgInt,
                         IntegralMode integralMode = IntegralMode::ProjectThenExtract,
                         bool printHeader = true){
    if (printHeader) std::cout << "\n[ExtractObservable2D] Starting extraction: " << extractionName << std::endl;
    // Create subdirectories for organized output
    TDirectory* dirBase = parentDir->mkdir(extractionName);
    TDirectory* dirFits = dirBase->mkdir("MassFits");
    TDirectory* dirRingMass = dirBase->mkdir("RingObservable_vs_Mass");
    // TDirectory* dirBkgFits = dirBase->mkdir("BackgroundFits");
    TDirectory* dirResults = dirBase->mkdir("Results");
    TDirectory* dirResultsSim = dirBase->mkdir("ResultsCombinedFit");
    TDirectory* dirDiagnostics = dirBase->mkdir("Diagnostics");

    // Build the corrected hNum TH2D from the TProfile2D
        // (that is, a TH2D with the correct "error of the mean" error bars of <R> = \sum_i R_i / N_\Lambda, instead of the errors ROOT assigned via "sqrt(\sum_i R_i^2)")
        // From this point on, h2dNum is used exactly as before, but now carries correct errors.
    TH2D* h2dNum = BuildNumFromProfile(p2dRingObs, h2dCounts, Form("hNumCorrectedError_%s", extractionName.Data()));
    if (!h2dNum) {
        std::cerr << "[ExtractObservable2D] ERROR: BuildNumFromProfile returned null. Skipping.\n";
        return;
    }

    // // DEBUG: confirm h2dNum is fresh and correct
    // if (h2dNum) {
    //     // Project the full Y axis to check the integral
    //     TH1D* dbgProj = h2dNum->ProjectionY("dbgProj_tmp", 1, h2dNum->GetNbinsX(), "e");
    //     std::cout << "[DEBUG " << extractionName << "] h2dNum=" << h2dNum
    //             << " name=" << h2dNum->GetName()
    //             << " hNumProj integral=" << dbgProj->Integral() << std::endl;
    //     delete dbgProj;
    //     std::cout << "[DEBUG " << extractionName << "] p2dRingObsIntegral: " << p2dRingObs->Integral() << std::endl;
    //     std::cout << "[DEBUG " << extractionName << "] h2dCountsIntegral: " << h2dCounts->Integral() << std::endl;
    // }

    int nBins = h2dCounts->GetNbinsX();

    // Switch to the results directory before booking anything below:
    // Without this, the objects would land in whatever directory was current,
    // which is harmless for the numbers (they are explicitly Write()-n further down) but makes
    // ownership harder to reason about.
    dirResults->cd();

    // Prepare output histograms (Cloning ProjectionX to keep exact angular binning)
    TH1D* hSigYield = (TH1D*)h2dCounts->ProjectionX(Form("hSigYield_%s", extractionName.Data()));
    hSigYield->Reset();
    hSigYield->SetTitle(Form("Signal Yield vs %s;%s;Counts", axisTitle.Data(), axisTitle.Data()));
    
    TH1D* hBkgYield = (TH1D*)hSigYield->Clone(Form("hBkgYield_%s", extractionName.Data()));
    hBkgYield->SetTitle(Form("Background Yield vs %s;%s;Counts", axisTitle.Data(), axisTitle.Data()));
    
    TH1D* hRSig = (TH1D*)hSigYield->Clone(Form("hRSig_%s", extractionName.Data()));
    hRSig->SetTitle(Form("Extracted <R>_{signal} vs %s;%s;<R>_{S}", axisTitle.Data(), axisTitle.Data()));
    
    TH1D* hRBkg = (TH1D*)hSigYield->Clone(Form("hRBkg_%s", extractionName.Data()));
    hRBkg->SetTitle(Form("Background <R>_{bkg} vs %s;%s;<R>_{B}", axisTitle.Data(), axisTitle.Data()));

    // --- QA HISTOGRAMS: statistical significance and purity ---
    TH1D* hPurity = (TH1D*)hSigYield->Clone(Form("hPurity_%s", extractionName.Data()));
    hPurity->SetTitle(Form("Purity S/(S+B) vs %s;%s;Purity", axisTitle.Data(), axisTitle.Data()));
    
    TH1D* hSigStat = (TH1D*)hSigYield->Clone(Form("hSigStat_%s", extractionName.Data()));
    hSigStat->SetTitle(Form("Significance S/#sqrt{S+B} vs %s;%s;Significance", axisTitle.Data(), axisTitle.Data()));

    // --- DIAGNOSTICS: what the configured windows actually became, bin by bin ---
    // The mass axis is finite and variable-width, so a request for "mu +/- N sigma" comes back
    // snapped to bin edges, asymmetric, and different in every angular bin (mu and sigma are
    // refitted each time). These record what the extraction really used.
    dirDiagnostics->cd();
    TH1D* hAchievedNSigmaLow = (TH1D*)hSigYield->Clone(Form("hAchievedNSigmaLow_%s", extractionName.Data()));
    hAchievedNSigmaLow->SetTitle(Form("Achieved signal window, low side (requested %.2f) vs %s;%s;(#mu - x_{low}) / #sigma",
                                      cfg.nSigmaPeak, axisTitle.Data(), axisTitle.Data()));
    TH1D* hAchievedNSigmaHigh = (TH1D*)hSigYield->Clone(Form("hAchievedNSigmaHigh_%s", extractionName.Data()));
    hAchievedNSigmaHigh->SetTitle(Form("Achieved signal window, high side (requested %.2f) vs %s;%s;(x_{high} - #mu) / #sigma",
                                       cfg.nSigmaPeak, axisTitle.Data(), axisTitle.Data()));
    TH1D* hAchievedCoverage = (TH1D*)hSigYield->Clone(Form("hAchievedCoverage_%s", extractionName.Data()));
    hAchievedCoverage->SetTitle(Form("Gaussian coverage of the achieved signal window vs %s;%s;Coverage",
                                     axisTitle.Data(), axisTitle.Data()));
    TH1D* hSidebandPointsLeft = (TH1D*)hSigYield->Clone(Form("hSidebandPointsLeft_%s", extractionName.Data()));
    hSidebandPointsLeft->SetTitle(Form("Sideband points used, left of the peak vs %s;%s;Points",
                                       axisTitle.Data(), axisTitle.Data()));
    TH1D* hSidebandPointsRight = (TH1D*)hSigYield->Clone(Form("hSidebandPointsRight_%s", extractionName.Data()));
    hSidebandPointsRight->SetTitle(Form("Sideband points used, right of the peak vs %s;%s;Points",
                                        axisTitle.Data(), axisTitle.Data()));
    dirResults->cd();

    // <R>_S - <R>_B in each angular bin. Stored rather than left to be subtracted downstream,
    // because the two share the sideband primitives (M and B) and are strongly correlated: naive
    // quadrature UNDERSTATES this particular difference, which is the dangerous direction.
    dirResults->cd();
    TH1D* hRSigMinusRBkg = (TH1D*)hSigYield->Clone(Form("hRSigMinusRBkg_%s", extractionName.Data()));
    hRSigMinusRBkg->SetTitle(Form("<R>_{S} - <R>_{B} vs %s;%s;<R>_{S} - <R>_{B}",
                                  axisTitle.Data(), axisTitle.Data()));

    // Per-bin results kept for the angle combination below. Only filled when the extraction
    // succeeded, so a rejected bin simply does not contribute rather than contributing a zero.
    std::vector<PeakWindowYields> perBinYields;
    std::vector<double> perBinRSig, perBinRSigErr, perBinRBkg, perBinRBkgErr;
    std::vector<double> perBinRMeas, perBinRMeasErr, perBinSigCounts;

    // Running summary for the one log line printed after the loop.
    int    nBinsExtracted = 0, nSidebandStarved = 0;
    double achievedLowMin = 1e9, achievedLowMax = -1e9;
    double achievedHighMin = 1e9, achievedHighMax = -1e9;
    int    sidebandPointsMin = 1000000, sidebandPointsMax = -1;

    // --- SIMULTANEOUS FIT QA HISTOGRAMS ---
    TH1D* hSigYield_Sim = (TH1D*)h2dCounts->ProjectionX(Form("hSigYield_Sim_%s", extractionName.Data()));
    hSigYield_Sim->Reset();
    hSigYield_Sim->SetTitle(Form("SimFit Signal Yield vs %s;%s;Counts", axisTitle.Data(), axisTitle.Data()));
    
    TH1D* hBkgYield_Sim = (TH1D*)hSigYield_Sim->Clone(Form("hBkgYield_Sim_%s", extractionName.Data()));
    hBkgYield_Sim->SetTitle(Form("SimFit Background Yield vs %s;%s;Counts", axisTitle.Data(), axisTitle.Data()));
    
    TH1D* hRSig_Sim = (TH1D*)hSigYield_Sim->Clone(Form("hRSig_Sim_%s", extractionName.Data()));
    hRSig_Sim->SetTitle(Form("SimFit Extracted <R>_{signal} vs %s;%s;<R>_{S}", axisTitle.Data(), axisTitle.Data()));
    
    TH1D* hRBkg_Sim = (TH1D*)hSigYield_Sim->Clone(Form("hRBkg_Sim_%s", extractionName.Data()));
    hRBkg_Sim->SetTitle(Form("SimFit Background <R>_{bkg} vs %s;%s;<R>_{B}", axisTitle.Data(), axisTitle.Data()));

    TH1D* hPurity_Sim = (TH1D*)hSigYield_Sim->Clone(Form("hPurity_Sim_%s", extractionName.Data()));
    hPurity_Sim->SetTitle(Form("SimFit Purity S/(S+B) vs %s;%s;Purity", axisTitle.Data(), axisTitle.Data()));
    
    TH1D* hSigStat_Sim = (TH1D*)hSigYield_Sim->Clone(Form("hSigStat_Sim_%s", extractionName.Data()));
    hSigStat_Sim->SetTitle(Form("SimFit Significance S/#sqrt{S+B} vs %s;%s;Significance", axisTitle.Data(), axisTitle.Data()));

        // Integrated fit result:
    TH1D* hIntegratedRSig_Sim = new TH1D(Form("hIntegratedRSig_Sim_%s", extractionName.Data()), 
                                         Form("Integrated R_{S} (SimFit) - %s; ;R_{S}", extractionName.Data()), 
                                         1, 0, 1);
    TH1D* hIntegratedRBkg_Sim = new TH1D(Form("hIntegratedRBkg_Sim_%s", extractionName.Data()), 
                                         Form("Integrated R_{B} (SimFit) - %s; ;R_{B}", extractionName.Data()), 
                                         1, 0, 1);

    // =========================================================================================
    // Steps 4/5/6/6.5: Per-bin mass fit, sideband fit, signal extraction and SimFit QA
    // =========================================================================================
    // These used to be two separate loops: one that fitted every mass projection and stashed
    // (mu, sigma) in a vector, and a second one that re-fetched the projections from dirFits by
    // name to do the sideband work. That round-trip through a ROOT directory was the source of a
    // silent bug (the first bin was registered in whichever directory happened to be current, so
    // dirFits->Get() returned nullptr and bin 1 was dropped). Merging the two loops removes the
    // round-trip entirely: every projection now stays in a local pointer for its whole lifetime.
    //
    // All tuned numbers come from the caller-supplied config; the shared helpers at the top of
    // this file do the fitting, the sideband collection and the error propagation.
    if (printHeader) std::cout << "  -> Steps 4/5/6/6.5: Mass fits, sideband fits and extraction..." << std::endl;

    for (int iBin = 1; iBin <= nBins; ++iBin){
        // 1. Project the mass spectrum (Counts/Denominator) using "e" for error propagation
        TH1D* hMassProj = h2dCounts->ProjectionY(Form("hMass_Bin%d", iBin), iBin, iBin, "e");
        hMassProj->SetDirectory(nullptr); // Local to this iteration; deleted before leaving it

        TAxis* xAxis = h2dCounts->GetXaxis();
        double xLow  = xAxis->GetBinLowEdge(iBin);
        double xHigh = xAxis->GetBinUpEdge(iBin);
        hMassProj->SetTitle(Form("Angle bin %d [%.4f,%.4f)];M_{p#pi} (GeV/c^{2});Counts", iBin, xLow, xHigh));

        // --- STABILITY CHECK: Does this bin even have data? ---
        // If there are fewer than ~30 counts in the entire mass window, a good fit is probably
        // impossible (and if possible, useless).
        if (hMassProj->GetEntries() < cfg.minEntries || hMassProj->Integral() <= cfg.minIntegral) {
            if (printHeader) std::cout << "    Bin " << iBin << ": invalid mu/sigma, skipped.\n";
            // We don't even want to see these problematic histograms, so just skip right away!
            delete hMassProj;
            continue;
        }

        // 2. Project the Numerator. This carries proper sigma(Sum_R_i) errors from the TProfile.
        TH1D* hNumProj = h2dNum->ProjectionY(Form("hNumExtBin,Angle[%.2f,%.2f),Bin%d", xLow, xHigh, iBin), iBin, iBin, "e");
        hNumProj->SetDirectory(nullptr); // Local to this iteration

        // 4. Density versions of both spectra.
        // It is actually just way more natural to fit the density histogram! Easier to perform QA
        // later on too, and it makes the fit insensitive to the variable bin widths.
        TH1D* hMassProjDensity = (TH1D*)hMassProj->Clone(Form("hMassDensity_Bin%d", iBin));
        hMassProjDensity->SetDirectory(nullptr);
        hMassProjDensity->SetTitle(Form("Angle bin %d [%.4f,%.4f);M_{p#pi} (GeV/c^{2});dN/dM", iBin, xLow, xHigh));
        hMassProjDensity->Scale(1.0, "width");

        // (changed to a density for proper normalization over variable bin widths. Notice that
        //  because <R>(m) = sum_i R_i(m) / N_Lambda(m) is a ratio, we can do this bin width
        //  scaling without much worry)
        // Named for what it is: d(Sum_R)/dm. Sum_R is EXTENSIVE, like counts, so scaling it by the
        // bin width is meaningful. It exists only to feed the simultaneous-fit QA; the sideband
        // background is no longer fitted to it.
        TH1D* hNumProjDensity = (TH1D*)hNumProj->Clone(Form("hSumRDensity_Bin%d", iBin));
        hNumProjDensity->SetDirectory(nullptr);
        hNumProjDensity->Scale(1.0, "width");

        // 5. Locate the peak. Only (mu, sigma) are kept; the pol2 half of this fit is discarded and
        //    the background is refitted on the sidebands alone just below.
        PeakFitOutcome peak = FitMassPeak(hMassProjDensity, massMin, massMax, cfg, Form("fit_Bin%d", iBin));

        dirFits->cd();
        hMassProjDensity->GetListOfFunctions()->Add(peak.fit); // Histogram takes ownership of the TF1
            // Save the histogram only if the projection is not empty!
        if (hMassProj->GetEntries() > 0) hMassProjDensity->Write(); // For QA

        if (!peak.valid) {
            if (printHeader) std::cout << "    Bin " << iBin << ": invalid mu/sigma, skipped.\n";
            delete hMassProjDensity; delete hNumProjDensity; delete hNumProj; delete hMassProj;
            continue;
        }
        double mu = peak.mu;
        double sigma = peak.sigma;

        // 6. Discontinuous sideband background fits, using TGraphErrors (TH1D's are bad at that!).
        //    The counts and numerator graphs are collected independently on purpose: they are only
        //    used to obtain a background shape, so there is no need for them to share points, and
        //    keeping them independent avoids discarding useful counts points whenever the numerator
        //    is invalid in that bin (and vice versa).
        SidebandSpan spanBkg;
        TGraphErrors* grBkg = BuildSidebandGraph(hMassProj, mu, sigma, massMin, massMax, cfg,
                                                 SidebandGraphKind::Counts, &spanBkg,
                                                 Form("grBkg_Bin%d", iBin));
        grBkg->SetTitle(Form("Sideband Bkg Bin %d;M_{p#pi} (GeV/c^{2});Counts/BinWidth", iBin));

        TGraphErrors* grRingBkg = BuildRingSidebandGraph(hMassProj, hNumProj, mu, sigma,
                                                         massMin, massMax, cfg,
                                                         Form("grRingBkg_Bin%d", iBin));
        grRingBkg->SetTitle(Form("<R> sideband points, bin %d;M_{p#pi} (GeV/c^{2});<R>_{bkg}", iBin));

        // --- STABILITY CHECK 1: Enough points for a pol2? ---
        // Adding some stability checks because the 3D histogram projections have much smaller
        // statistics to be handled! A pol2 has 3 parameters, so 5 points leaves 2 degrees of freedom.
        // --- STABILITY CHECK 1.5: Are the sidebands just empty zeros? ---
        // If the background counts sum to (nearly) zero the matrix inversion fails and we get a
        // looooot of Minuit errors. spanBkg.rawCounts holds the RAW counts, which is what this
        // check needs: the graph's own Y values are densities now, so they cannot be summed here.
        if (grBkg->GetN() < cfg.minSidebandPoints || grRingBkg->GetN() < cfg.minSidebandPoints ||
            (cfg.minSidebandCounts > 0.0 && spanBkg.rawCounts <= cfg.minSidebandCounts)) {
            nSidebandStarved++;
            delete grBkg; delete grRingBkg;
            delete hMassProjDensity; delete hNumProjDensity; delete hNumProj; delete hMassProj;
            continue;
        }

        // Fit the discontinuous graphs with a polN, to cover the small curvature in Lambda background.
        // From QA plots you can see that the background is almost linear, maybe a bit quadratic,
        // when studying <R>(m_Lambda). Thus, regular signal extraction should work!
        TF1* bkgFitFunc = new TF1(Form("bkgFit_Bin%d", iBin),
                                  PolShifted(cfg.bkgPolOrder, cfg.muInitGuess),
                                  massMin, massMax);
        TFitResultPtr rBkg = grBkg->Fit(bkgFitFunc, "Q 0 S"); // "Q" = quiet, "0" = don't draw, "S" = return TFitResultPtr

        TF1* ringBkgFitFunc = new TF1(Form("ringBkgFit_Bin%d", iBin),
                                     PolShifted(cfg.ringBkgPolOrder, cfg.muInitGuess),
                                     massMin, massMax);
        TFitResultPtr rRingBkg = grRingBkg->Fit(ringBkgFitFunc, "Q 0 S");

        // --- STABILITY CHECK 2: Did the fits converge properly? ---
        if (!rBkg->IsValid() || !rRingBkg->IsValid()) {
            // Here you still need to delete the TF1s manually, as they have not yet been taken
            // ownership of by ROOT through an "Add()" call!
            delete bkgFitFunc; delete ringBkgFitFunc;
            delete grBkg; delete grRingBkg;
            delete hMassProjDensity; delete hNumProjDensity; delete hNumProj; delete hMassProj;
            continue;
        }

        dirFits->cd();
        grBkg->GetListOfFunctions()->Add(bkgFitFunc); // Attach fit for viewing in TBrowser
        grBkg->Write();

        // --- <R>(m) with the fitted sideband background drawn on it ---------------------------
        // Built here rather than earlier so the fitted function can be attached to it: this is the
        // one plot that shows whether the flat-<R> assumption behind sideband subtraction actually
        // holds in this angular bin. The fit is CLONED because the graph below takes ownership of
        // the original.
        TH1D* hRingVsMass = (TH1D*)hNumProj->Clone(Form("hRingVsMass_Bin%d", iBin));
        hRingVsMass->SetDirectory(nullptr);
        hRingVsMass->SetTitle(Form("<R> vs mass, bin %d [%.4f,%.4f);M_{p#pi} (GeV/c^{2});<R>",
                                   iBin, xLow, xHigh));
        hRingVsMass->Divide(hNumProj, hMassProj, 1.0, 1.0, "B");
        hRingVsMass->GetListOfFunctions()->Add((TF1*)ringBkgFitFunc->Clone());

        dirRingMass->cd();
        hRingVsMass->Write();
        grRingBkg->GetListOfFunctions()->Add(ringBkgFitFunc); // Graph takes ownership
        grRingBkg->Write();
        delete hRingVsMass;

        // 7. Signal extraction. Every count, subtraction and uncertainty lives in this one call.
        PeakWindowYields y = ComputePeakWindowYields(hMassProj, hNumProj, mu, sigma, cfg.nSigmaPeak,
                                                     bkgFitFunc, rBkg, ringBkgFitFunc, rRingBkg);

        if (!y.valid) {
            // failStage 1 = no counts in the window at all, 2 = background exceeded the peak
            if (y.failStage == 2) std::cout << "    Bin " << iBin << ": non-positive signal, skipped.\n";
            delete grBkg; delete grRingBkg;
            delete hMassProjDensity; delete hNumProjDensity; delete hNumProj; delete hMassProj;
            continue;
        }

        // Fill Histograms
        hSigYield->SetBinContent(iBin, y.sigCounts);
        hSigYield->SetBinError(iBin, y.errSigCounts);

        hBkgYield->SetBinContent(iBin, y.bkgCounts);
        hBkgYield->SetBinError(iBin, y.errBkgCounts);

        hRSig->SetBinContent(iBin, y.R_S);
        hRSig->SetBinError(iBin, y.errR_S);

        hRBkg->SetBinContent(iBin, y.R_B);
        hRBkg->SetBinError(iBin, std::sqrt(y.varR_B));

        hPurity->SetBinContent(iBin, y.purity);
        hPurity->SetBinError(iBin, y.errPurity);

        hSigStat->SetBinContent(iBin, y.significance);
        hSigStat->SetBinError(iBin, y.errSignificance);

        // Diagnostics: what the windows became in this bin. Errors are set to zero because these
        // are exact statements about the binning, not measurements.
        hAchievedNSigmaLow->SetBinContent(iBin, y.nSigmaAchievedLow);
        hAchievedNSigmaLow->SetBinError(iBin, 0.0);
        hAchievedNSigmaHigh->SetBinContent(iBin, y.nSigmaAchievedHigh);
        hAchievedNSigmaHigh->SetBinError(iBin, 0.0);
        hAchievedCoverage->SetBinContent(iBin, y.coverage);
        hAchievedCoverage->SetBinError(iBin, 0.0);
        hSidebandPointsLeft->SetBinContent(iBin, spanBkg.nPointsLeft);
        hSidebandPointsLeft->SetBinError(iBin, 0.0);
        hSidebandPointsRight->SetBinContent(iBin, spanBkg.nPointsRight);
        hSidebandPointsRight->SetBinError(iBin, 0.0);

        hRSigMinusRBkg->SetBinContent(iBin, y.diffSigMinusBkg);
        hRSigMinusRBkg->SetBinError(iBin, y.errDiffSigMinusBkg);

        perBinYields.push_back(y);
        perBinRSig.push_back(y.R_S);       perBinRSigErr.push_back(y.errR_S);
        perBinRBkg.push_back(y.R_B);       perBinRBkgErr.push_back(y.errR_B);
        perBinRMeas.push_back(y.R_peak);   perBinRMeasErr.push_back(y.errR_peak);
        perBinSigCounts.push_back(y.sigCounts);

        nBinsExtracted++;
        achievedLowMin = std::min(achievedLowMin, y.nSigmaAchievedLow);
        achievedLowMax = std::max(achievedLowMax, y.nSigmaAchievedLow);
        achievedHighMin = std::min(achievedHighMin, y.nSigmaAchievedHigh);
        achievedHighMax = std::max(achievedHighMax, y.nSigmaAchievedHigh);
        int nSideband = spanBkg.nPointsLeft + spanBkg.nPointsRight;
        sidebandPointsMin = std::min(sidebandPointsMin, nSideband);
        sidebandPointsMax = std::max(sidebandPointsMax, nSideband);

        // =========================================================================
        // Step 6.5: Simultaneous Fit QA
        // =========================================================================
        // PerformSimultaneousFitQA receives the corrected hNumProjDensity (the TProfile fix is
        // entirely in what the numerator projection contains).
        // init_sigYield (= sigCounts from sideband subtraction) remains a valid initial guess: it is still
        // an absolute count, which PerformSimultaneousFitQA converts internally to a density
        // amplitude via sigCounts / (sigma * sqrt(2*pi)). We can still pass sigCounts because the
        // normalization is done against hMassProjDensity, independently of bin widths.
        SimFitResult qaResult = PerformSimultaneousFitQA(hMassProjDensity, hNumProjDensity, massMin, massMax, y.sigCounts, mu, sigma, bkgFitFunc, y.R_S, y.R_B, cfg.nSigmaPeak);
        if (qaResult.status == 0) {
            hSigYield_Sim->SetBinContent(iBin, qaResult.SigYield);
            hSigYield_Sim->SetBinError(iBin, qaResult.err_SigYield);

            hBkgYield_Sim->SetBinContent(iBin, qaResult.BkgYield);
            hBkgYield_Sim->SetBinError(iBin, qaResult.err_BkgYield);

            hRSig_Sim->SetBinContent(iBin, qaResult.R_S);
            hRSig_Sim->SetBinError(iBin, qaResult.err_R_S);

            hRBkg_Sim->SetBinContent(iBin, qaResult.R_B);
            hRBkg_Sim->SetBinError(iBin, qaResult.err_R_B);

            hPurity_Sim->SetBinContent(iBin, qaResult.Purity);
            hPurity_Sim->SetBinError(iBin, qaResult.err_Purity);

            hSigStat_Sim->SetBinContent(iBin, qaResult.Significance);
            hSigStat_Sim->SetBinError(iBin, qaResult.err_Significance);
        }

        // Note: we do NOT delete bkgFitFunc or ringBkgFitFunc here, because ROOT took ownership of
        // them when they were added to the graphs' lists of functions. Deleting them would crash!
        // Likewise, peak.fit is owned by hMassProjDensity.
        delete grBkg;
        delete grRingBkg;
        delete hNumProjDensity;
        delete hMassProjDensity;
        delete hNumProj;
        delete hMassProj;
    }

    // -----------------------------------------------------------------------------------------
    // Achieved-window report.
    // -----------------------------------------------------------------------------------------
    // One line per extraction rather than one per bin: the configured value is a single number,
    // but what it becomes after bin snapping varies from bin to bin, so the range is the honest
    // thing to quote. This is also the place where a starved sideband band announces itself --
    // silently invalidating every bin is exactly the failure mode that is hard to notice.
    if (printHeader) {
        if (nBinsExtracted > 0) {
            std::cout << Form("  -> Windows [%s]: signal requested %.2f sigma, achieved low %.2f-%.2f, "
                              "high %.2f-%.2f (coverage varies with it); sidebands %.2f sigma to %s, "
                              "%d-%d points; %d/%d bins extracted",
                              extractionName.Data(), cfg.nSigmaPeak,
                              achievedLowMin, achievedLowMax, achievedHighMin, achievedHighMax,
                              cfg.nSigmaExclusion,
                              cfg.nSigmaExclusionOuter > 0.0 ? Form("%.2f sigma", cfg.nSigmaExclusionOuter) : "the axis edge",
                              sidebandPointsMin, sidebandPointsMax, nBinsExtracted, nBins)
                      << std::endl;
        }
        if (nSidebandStarved > 0) {
            std::cout << Form("  -> WARNING [%s]: %d/%d bins had too few sideband points or counts "
                              "(need >= %d points and > %.1f counts). Widen the band or coarsen the "
                              "mass axis.", extractionName.Data(), nSidebandStarved, nBins,
                              cfg.minSidebandPoints, cfg.minSidebandCounts)
                      << std::endl;
        }
    }

    dirDiagnostics->cd();
    hAchievedNSigmaLow->Write();
    hAchievedNSigmaHigh->Write();
    hAchievedCoverage->Write();
    hSidebandPointsLeft->Write();
    hSidebandPointsRight->Write();

    dirResults->cd();
    hRSigMinusRBkg->Write();
    delete hRSigMinusRBkg;

    // Resetting stats boxes to get an estimate of Signal vs Background Lambda counts:
    hSigYield->ResetStats();
    hBkgYield->ResetStats();

    hSigYield->Write();
    hBkgYield->Write();
    hRSig->Write();
    hRBkg->Write();
    hPurity->Write();
    hSigStat->Write();

    // Saving histograms from the Simultaneous Ring+Mass fit:
    dirResultsSim->cd();
    hSigYield_Sim->ResetStats();
    hBkgYield_Sim->ResetStats();

    hSigYield_Sim->Write();
    hBkgYield_Sim->Write();
    hRSig_Sim->Write();
    hRBkg_Sim->Write();
    hPurity_Sim->Write();
    hSigStat_Sim->Write();

    // =========================================================================================
    // Step 7.9: Angle-COMBINED results, built from the per-bin extractions
    // =========================================================================================
    // Produced whenever there is more than one usable angular bin. See the header block on
    // CombineYieldsSignalWeighted for why combining beats projecting, and for what the two flavours
    // mean. In CombinePerBin mode these REPLACE Steps 8 and 9 rather than sitting beside them.
    if (perBinYields.size() >= 2) {
        TDirectory* dirCombined = EnsureDir(dirBase, "IntegratedCombined");

        PeakWindowYields aggSW = CombineYieldsSignalWeighted(perBinYields);

        FlatAverage flatSig = AverageOverBins(perBinRSig, perBinRSigErr);
        FlatAverage flatBkg = AverageOverBins(perBinRBkg, perBinRBkgErr);
        FlatAverage flatMeas = AverageOverBins(perBinRMeas, perBinRMeasErr);

        FlatAverage aeeSig = WeightedMinusFlat(perBinRSig, perBinRSigErr, perBinSigCounts);
        FlatAverage aeeBkg = WeightedMinusFlat(perBinRBkg, perBinRBkgErr, perBinSigCounts);
        FlatAverage aeeMeas = WeightedMinusFlat(perBinRMeas, perBinRMeasErr, perBinSigCounts);

        auto scalarOut = [&](const char* stem, const char* yTitle, double v, double e) {
            TH1D* h = new TH1D(Form("%s_%s", stem, extractionName.Data()),
                               Form("%s; ;%s", extractionName.Data(), yTitle), 1, 0, 1);
            h->SetDirectory(nullptr);
            h->SetBinContent(1, v);
            h->SetBinError(1, e);
            dirCombined->cd();
            h->Write();
            delete h;
        };

        if (aggSW.valid) {
            scalarOut("hCombinedRMeas_SignalWeighted", "<R>_{measured}", aggSW.R_peak, aggSW.errR_peak);
            scalarOut("hCombinedRSig_SignalWeighted", "<R>_{S}", aggSW.R_S, aggSW.errR_S);
            scalarOut("hCombinedRBkg_SignalWeighted", "<R>_{B}", aggSW.R_B, aggSW.errR_B);
            scalarOut("hCombinedDiffSigMinusBkg_SignalWeighted", "<R>_{S} - <R>_{B}",
                      aggSW.diffSigMinusBkg, aggSW.errDiffSigMinusBkg);
            scalarOut("hCombinedDiffMeasMinusSig_SignalWeighted", "<R>_{measured} - <R>_{S}",
                      aggSW.diffPeakMinusSig, aggSW.errDiffPeakMinusSig);
            scalarOut("hCombinedDiffMeasMinusBkg_SignalWeighted", "<R>_{measured} - <R>_{B}",
                      aggSW.diffPeakMinusBkg, aggSW.errDiffPeakMinusBkg);
            scalarOut("hCombinedPurity_SignalWeighted", "S/(S+B)", aggSW.purity, aggSW.errPurity);
            scalarOut("hCombinedSignificance_SignalWeighted", "S/#sqrt{S+B}",
                      aggSW.significance, aggSW.errSignificance);
        }
        if (flatMeas.valid) scalarOut("hCombinedRMeas_Flat", "<R>_{measured}", flatMeas.value, flatMeas.error);
        if (flatSig.valid)  scalarOut("hCombinedRSig_Flat", "<R>_{S}", flatSig.value, flatSig.error);
        if (flatBkg.valid)  scalarOut("hCombinedRBkg_Flat", "<R>_{B}", flatBkg.value, flatBkg.error);

        // Azimuthal Efficiency Effect: signal-weighted minus flat, per quantity.
        if (aeeMeas.valid) scalarOut("hAEE_RMeas", "<R>_{measured}: weighted - flat", aeeMeas.value, aeeMeas.error);
        if (aeeSig.valid)  scalarOut("hAEE_RSig", "<R>_{S}: weighted - flat", aeeSig.value, aeeSig.error);
        if (aeeBkg.valid)  scalarOut("hAEE_RBkg", "<R>_{B}: weighted - flat", aeeBkg.value, aeeBkg.error);

        // Bookkeeping, so a combined number can always be traced back to how many bins fed it.
        TH1D* hNUsed = new TH1D(Form("hCombinedBinsUsed_%s", extractionName.Data()),
                                Form("%s; ;Angular bins", extractionName.Data()), 2, 0, 2);
        hNUsed->SetDirectory(nullptr);
        hNUsed->GetXaxis()->SetBinLabel(1, "used");
        hNUsed->GetXaxis()->SetBinLabel(2, "total");
        hNUsed->SetBinContent(1, static_cast<double>(perBinYields.size()));
        hNUsed->SetBinContent(2, static_cast<double>(nBins));
        dirCombined->cd();
        hNUsed->Write();
        delete hNUsed;

        if (printHeader && aggSW.valid) {
            std::cout << Form("  -> Combined [%s] over %d/%d angular bins:",
                              extractionName.Data(), (int)perBinYields.size(), nBins) << std::endl;
            std::cout << Form("       signal-weighted: <R>_S = %+.6f +/- %.6f, <R>_B = %+.6f +/- %.6f",
                              aggSW.R_S, aggSW.errR_S, aggSW.R_B, aggSW.errR_B) << std::endl;
            if (flatSig.valid && flatBkg.valid)
                std::cout << Form("       flat acceptance: <R>_S = %+.6f +/- %.6f, <R>_B = %+.6f +/- %.6f",
                                  flatSig.value, flatSig.error, flatBkg.value, flatBkg.error) << std::endl;
            if (aeeSig.valid)
                std::cout << Form("       AEE (weighted - flat): <R>_S %+.6f +/- %.6f, <R>_B %+.6f +/- %.6f",
                                  aeeSig.value, aeeSig.error, aeeBkg.value, aeeBkg.error) << std::endl;
        }
    }

    if (integralMode == IntegralMode::CombinePerBin) {
        // Steps 8 and 9 both collapse every angular bin onto the mass axis before extracting, which
        // is precisely what the angular split exists to avoid. Where <R> changes sign across the
        // angular variable, that projection mixes opposite-sign numerator backgrounds into one
        // spectrum and then asks a low-order polynomial to describe it. The combined results above
        // are the same measurement done in the order that keeps every fit well conditioned, so the
        // projected versions are not merely redundant here -- they would be worse.
        if (printHeader)
            std::cout << "  -> Steps 8 and 9 skipped: the angle-integrated result comes from "
                         "IntegratedCombined/, not from a projection over the angular bins."
                      << std::endl;
        delete h2dNum;
        return;
    }

    // =========================================================================================
    // Step 8: Calculate the Angle-Integrated Ring Observable
    // =========================================================================================
    // Repeating the same procedure, but for an integrated observable!
    if (printHeader) std::cout << "  -> Step 8: Calculating Integrated Observable..." << std::endl;
    
    // Project the entire TH2D onto the Y-axis (Mass) to integrate over all angles
    TH1D* hMassInt = h2dCounts->ProjectionY(Form("hMassInt_%s", extractionName.Data()), 1, nBins, "e");
        // Project the corrected h2dNum (built from TProfile) to get proper Sum and errors:
    TH1D* hNumInt = h2dNum->ProjectionY(Form("hNumInt_%s", extractionName.Data()), 1, nBins, "e");

    // Density versions, for the fits. BEHAVIOUR CHANGE: this path used to fit the RAW counts while
    // the per-bin path fitted the density spectrum. On a variable-width mass axis those are not the
    // same fit, and the density one is the correct treatment (it is what makes TF1::Integral over
    // the signal window return counts with no bin-width bookkeeping). The two paths were never
    // meant to disagree about this; the per-bin path simply got the fix first.
    TH1D* hMassIntDensity = (TH1D*)hMassInt->Clone(Form("hMassIntDensity_%s", extractionName.Data()));
    hMassIntDensity->SetDirectory(nullptr);
    hMassIntDensity->Scale(1.0, "width");

    // -----------------------------------------------------------------------------------------
    // Integrated results. One single-bin histogram per quantity, all written into
    // IntegratedSummary/<extraction>/ so that downstream plotting code has one flat, predictable
    // place to read from instead of digging through the per-bin folders.
    // -----------------------------------------------------------------------------------------
    // One single-bin histogram per quantity. The x axis is empty on purpose: these are scalars, and
    // giving them a real axis would only invite someone to plot them against it.
    auto makeScalarHisto = [&](const char* stem, const char* yTitle) {
        TH1D* h = new TH1D(Form("%s_%s", stem, extractionName.Data()),
                           Form("%s; ;%s", extractionName.Data(), yTitle), 1, 0, 1);
        h->SetDirectory(nullptr);
        return h;
    };

    TH1D* hIntegratedRSig  = makeScalarHisto("hIntegratedRSig",  "<R>_{S}");
    TH1D* hIntegratedRBkg  = makeScalarHisto("hIntegratedRBkg",  "<R>_{B}");
    TH1D* hIntegratedRMeas = makeScalarHisto("hIntegratedRMeas", "<R>_{measured}");
    TH1D* hIntegratedDiffMeasMinusSig = makeScalarHisto("hIntegratedDiffMeasMinusSig",
                                                        "<R>_{measured} - <R>_{S}");
    TH1D* hIntegratedDiffMeasMinusBkg = makeScalarHisto("hIntegratedDiffMeasMinusBkg",
                                                        "<R>_{measured} - <R>_{B}");
    TH1D* hIntegratedDiffSigMinusBkg  = makeScalarHisto("hIntegratedDiffSigMinusBkg",
                                                        "<R>_{S} - <R>_{B}");
    TH1D* hIntegratedPurity       = makeScalarHisto("hIntegratedPurity", "S/(S+B)");
    TH1D* hIntegratedSignificance = makeScalarHisto("hIntegratedSignificance", "S/#sqrt{S+B}");

    // Two labelled bins rather than two histograms: signal and background yields are always read
    // together, and their correlation is not carried anywhere, so splitting them buys nothing.
    TH1D* hIntegratedYields = new TH1D(Form("hIntegratedYields_%s", extractionName.Data()),
                                       Form("%s; ;Counts in the signal window", extractionName.Data()),
                                       2, 0, 2);
    hIntegratedYields->SetDirectory(nullptr);
    hIntegratedYields->GetXaxis()->SetBinLabel(1, "Signal");
    hIntegratedYields->GetXaxis()->SetBinLabel(2, "Background");

    // What the configured windows actually became on this mass axis, so that a plot made from the
    // numbers above can state the window it really used rather than the one that was requested.
    TH1D* hAchievedWindow = new TH1D(Form("hAchievedWindow_%s", extractionName.Data()),
                                     Form("%s: achieved extraction windows; ;Value", extractionName.Data()),
                                     6, 0, 6);
    hAchievedWindow->SetDirectory(nullptr);
    hAchievedWindow->GetXaxis()->SetBinLabel(1, "nSigmaLow");
    hAchievedWindow->GetXaxis()->SetBinLabel(2, "nSigmaHigh");
    hAchievedWindow->GetXaxis()->SetBinLabel(3, "coverage");
    hAchievedWindow->GetXaxis()->SetBinLabel(4, "sidebandPtsLeft");
    hAchievedWindow->GetXaxis()->SetBinLabel(5, "sidebandPtsRight");
    hAchievedWindow->GetXaxis()->SetBinLabel(6, "bkgPolOrder");

    // Only proceed if we have enough global statistics to fit a peak.
    if (hMassInt->GetEntries() >= cfgInt.minEntries && hMassInt->Integral() > cfgInt.minIntegral) {
        // 1. Global Mass Fit, on the density spectrum (see the note above).
        PeakFitOutcome peakInt = FitMassPeak(hMassIntDensity, massMin, massMax, cfgInt,
                                             Form("fitMassInt_%s", extractionName.Data()));

        if (peakInt.valid) {
            double muInt = peakInt.mu;
            double sigmaInt = peakInt.sigma;

            // 2. Global Sideband Extraction.
            SidebandSpan spanInt;
            TGraphErrors* grBkgInt = BuildSidebandGraph(hMassInt, muInt, sigmaInt, massMin, massMax,
                                                        cfgInt, SidebandGraphKind::Counts, &spanInt,
                                                        Form("grBkgInt_%s", extractionName.Data()));
            TGraphErrors* grRingBkgInt = BuildRingSidebandGraph(hMassInt, hNumInt, muInt, sigmaInt,
                                                               massMin, massMax, cfgInt,
                                                               Form("grRingBkgInt_%s", extractionName.Data()));

            if (grBkgInt->GetN() >= cfgInt.minSidebandPoints && grRingBkgInt->GetN() >= cfgInt.minSidebandPoints) {
                TF1* bkgFitInt = new TF1(Form("bkgFitInt_%s", extractionName.Data()),
                                         PolShifted(cfgInt.bkgPolOrder, cfgInt.muInitGuess), massMin, massMax);
                TF1* ringBkgFitInt = new TF1(Form("ringBkgFitInt_%s", extractionName.Data()),
                                            PolShifted(cfgInt.ringBkgPolOrder, cfgInt.muInitGuess), massMin, massMax);

                TFitResultPtr rBkgInt = grBkgInt->Fit(bkgFitInt, "Q 0 S");
                TFitResultPtr rRingBkgInt = grRingBkgInt->Fit(ringBkgFitInt, "Q 0 S");

                if (rBkgInt->IsValid() && rRingBkgInt->IsValid()) {
                    // 3. Global integration and error propagation -- the very same helper used by
                    //    the per-bin loop, so the two paths can never drift apart again.
                    PeakWindowYields yInt = ComputePeakWindowYields(hMassInt, hNumInt, muInt, sigmaInt,
                                                                    cfgInt.nSigmaPeak,
                                                                    bkgFitInt, rBkgInt,
                                                                    ringBkgFitInt, rRingBkgInt);
                    if (yInt.valid) {
                        hIntegratedRSig->SetBinContent(1, yInt.R_S);
                        hIntegratedRSig->SetBinError(1, yInt.errR_S);

                        hIntegratedRBkg->SetBinContent(1, yInt.R_B);
                        hIntegratedRBkg->SetBinError(1, yInt.errR_B);

                        // <R>_measured is the raw mean inside the signal window, background and
                        // all. It is the signal-extraction counterpart of the "Data" column in the
                        // summary plots, but NOT the same number: the summary plots average over
                        // whatever mass range the consumer profile covers, while this one is
                        // restricted to the achieved peak window recorded below.
                        hIntegratedRMeas->SetBinContent(1, yInt.R_peak);
                        hIntegratedRMeas->SetBinError(1, yInt.errR_peak);

                        // The three differences carry the correlations between the quantities
                        // above, which is why they are stored rather than left to be recomputed.
                        hIntegratedDiffMeasMinusSig->SetBinContent(1, yInt.diffPeakMinusSig);
                        hIntegratedDiffMeasMinusSig->SetBinError(1, yInt.errDiffPeakMinusSig);
                        hIntegratedDiffMeasMinusBkg->SetBinContent(1, yInt.diffPeakMinusBkg);
                        hIntegratedDiffMeasMinusBkg->SetBinError(1, yInt.errDiffPeakMinusBkg);
                        hIntegratedDiffSigMinusBkg->SetBinContent(1, yInt.diffSigMinusBkg);
                        hIntegratedDiffSigMinusBkg->SetBinError(1, yInt.errDiffSigMinusBkg);

                        hIntegratedPurity->SetBinContent(1, yInt.purity);
                        hIntegratedPurity->SetBinError(1, yInt.errPurity);
                        hIntegratedSignificance->SetBinContent(1, yInt.significance);
                        hIntegratedSignificance->SetBinError(1, yInt.errSignificance);

                        hIntegratedYields->SetBinContent(1, yInt.sigCounts);
                        hIntegratedYields->SetBinError(1, yInt.errSigCounts);
                        hIntegratedYields->SetBinContent(2, yInt.bkgCounts);
                        hIntegratedYields->SetBinError(2, yInt.errBkgCounts);

                        hAchievedWindow->SetBinContent(1, yInt.nSigmaAchievedLow);
                        hAchievedWindow->SetBinContent(2, yInt.nSigmaAchievedHigh);
                        hAchievedWindow->SetBinContent(3, yInt.coverage);
                        hAchievedWindow->SetBinContent(4, spanInt.nPointsLeft);
                        hAchievedWindow->SetBinContent(5, spanInt.nPointsRight);
                        hAchievedWindow->SetBinContent(6, cfgInt.bkgPolOrder);

                        if (printHeader) {
                            std::cout << Form("  -> Integrated [%s]: <R>_meas = %.5f +/- %.5f, "
                                              "<R>_S = %.5f +/- %.5f, <R>_B = %.5f +/- %.5f",
                                              extractionName.Data(),
                                              yInt.R_peak, yInt.errR_peak,
                                              yInt.R_S, yInt.errR_S,
                                              yInt.R_B, yInt.errR_B) << std::endl;
                            std::cout << Form("     window: requested %.2f sigma, achieved [-%.2f, +%.2f] "
                                              "(mass %.4f to %.4f, coverage %.4f); sidebands %d + %d points",
                                              cfgInt.nSigmaPeak,
                                              yInt.nSigmaAchievedLow, yInt.nSigmaAchievedHigh,
                                              yInt.xLow, yInt.xHigh, yInt.coverage,
                                              spanInt.nPointsLeft, spanInt.nPointsRight) << std::endl;
                        }
                    }
                }
                delete bkgFitInt; delete ringBkgFitInt;
            }
            else if (printHeader) {
                std::cout << Form("  -> WARNING [%s]: integrated sidebands have only %d/%d points "
                                  "(need >= %d). No integrated result produced.",
                                  extractionName.Data(), grBkgInt->GetN(), grRingBkgInt->GetN(),
                                  cfgInt.minSidebandPoints) << std::endl;
            }
            delete grBkgInt; delete grRingBkgInt;
        }
        delete peakInt.fit; // Never attached to a histogram here, so we own it
    } // end of step 8 conditionals

    // Written into Results/ as a per-extraction cross-check ONLY.
    // These are NOT the analysis integrated numbers and nothing downstream should read them: every
    // angular extraction integrates over the same candidates, so DeltaPhi and DeltaTheta give
    // identical values by construction and publishing both would be duplication dressed up as two
    // measurements. The integrated result now comes from ExtractIntegratedFromProfile, once per
    // proxy, in <variation>/IntegratedSummary/.
    dirResults->cd();
    hIntegratedRSig->Write();
    hIntegratedRBkg->Write();
    hIntegratedRMeas->Write();
    hIntegratedDiffMeasMinusSig->Write();
    hIntegratedDiffMeasMinusBkg->Write();
    hIntegratedDiffSigMinusBkg->Write();
    hIntegratedPurity->Write();
    hIntegratedSignificance->Write();
    hIntegratedYields->Write();
    hAchievedWindow->Write();

    delete hIntegratedRSig; delete hIntegratedRBkg; delete hIntegratedRMeas;
    delete hIntegratedDiffMeasMinusSig; delete hIntegratedDiffMeasMinusBkg;
    delete hIntegratedDiffSigMinusBkg;
    delete hIntegratedPurity; delete hIntegratedSignificance;
    delete hIntegratedYields; delete hAchievedWindow;

        // Clean up memory for the 1D projections:
    delete hMassIntDensity; delete hMassInt; delete hNumInt;

    // ==========================================================================================
    // Simultaneous Fit - Integrated results (integrated over all angles)
    // ==========================================================================================
    // Instead of averaging the R_S values from individual bins (which requires complex error propagation 
    // of weighted means), we sum the raw histograms first and fit ONCE.
    // 1. It provides a "Global Anchor" with high statistics.
    // 2. The error bar comes directly from the Global Hessian Matrix, automatically accounting 
    //    for correlations between Signal Yield and Polarization.
    if (printHeader) std::cout << "  -> Step 9: Integrated Simultaneous Fit..." << std::endl;

    // 1. Project Total Statistics (Summing over all angular bins 1 to nBins)
    TH1D* hTotalMass_Sim = h2dCounts->ProjectionY(Form("hTotalMass_Sim_%s", extractionName.Data()), 1, nBins, "e");
    TH1D* hTotalNum_Sim  = h2dNum->ProjectionY(Form("hTotalNum_Sim_%s", extractionName.Data()), 1, nBins, "e");
    // Should represent the same object for Delta Phi and Delta Theta, as it is integrated in deltaTheta and deltaPhi.
    // The integral of the TPrifle2Ds is different, though, as it retrieves an unweighted sum of per-bin means, with
    // no regard whatsoever for how many candidates are in each bin.

    // -----------------------------------------------------------------------------------------
    // CONDENSED STABILITY CHECKS
    // -----------------------------------------------------------------------------------------
    bool doIntegratedFit = true;
    // CHECK A: Data Existence & Matrix Safety
    // If there are fewer than 30 counts in the entire mass window, a good fit is probably impossible.
    // Also, if sumBkg is too low, the covariance matrix inversion for pol2 will fail.
    double totalInteg = hTotalMass_Sim->Integral();
    if (hTotalMass_Sim->GetEntries() < 30 || totalInteg <= 30) {
        if (printHeader) std::cout << "    [SimFit Integrated] Skipped: Too few entries (" << totalInteg << ")" << std::endl;
        doIntegratedFit = false;
    }
    // CHECK B: Ensure we have data in the sidebands (approximation of the "Enough points" check)
    // If the histogram is empty in the wings, the background parabola cannot be constrained.
    if (doIntegratedFit) {
        int nBinsH = hTotalMass_Sim->GetNbinsX();
        double leftWing  = hTotalMass_Sim->Integral(1, nBinsH/5); // First 20%
        double rightWing = hTotalMass_Sim->Integral(nBinsH - nBinsH/5, nBinsH); // Last 20%
        if (leftWing + rightWing < 5) { // Threshold for "empty zeros"
            if (printHeader) std::cout << "    [SimFit Integrated] Skipped: Sidebands are empty (Matrix Inversion Risk)." << std::endl;
            doIntegratedFit = false;
        }
    }
    // -----------------------------------------------------------------------------------------
    // PRE-FIT & SIMULTANEOUS EXTRACTION
    // -----------------------------------------------------------------------------------------
    // We only proceed if the histogram passed the stability checks
    if (doIntegratedFit) {
        // 2. Obtain robust initial guesses for the simultaneous fit using a standard pre-fit
        TF1* preFit = new TF1(Form("preFit_%s", extractionName.Data()), "gaus(0) + pol2(3)", massMin, massMax);

        // Transforming the mass histogram into a density histogram (easier QA):
        TH1D* hTotalMass_SimDensity = (TH1D*)hTotalMass_Sim->Clone(Form("hTotalMass_SimDensity_%s", extractionName.Data()));
            // Scale by bin width to get density
        hTotalMass_SimDensity->Scale(1.0, "width");

        // Same for numerator:
        TH1D* hTotalNum_SimDensity = (TH1D*)hTotalNum_Sim->Clone(Form("hTotalNum_SimDensity_%s", extractionName.Data()));
        hTotalNum_SimDensity->Scale(1.0, "width");
        
        // double maxVal = hTotalMass_Sim->GetMaximum();
        double maxVal = hTotalMass_SimDensity->GetMaximum();
        double estimatedMu = hTotalMass_SimDensity->GetBinCenter(hTotalMass_SimDensity->GetMaximumBin());
        
        // Estimate background from the edges (average of first and last bin)
        // (we NEED this type of initial guess to be better for the integrated fit because the statistics is just so much higher!)
        double firstBinC = hTotalMass_SimDensity->GetBinContent(1);
        double lastBinC  = hTotalMass_SimDensity->GetBinContent(hTotalMass_SimDensity->GetNbinsX());
        double estBkgLevel = (firstBinC + lastBinC) / 2.0;
        if (estBkgLevel < 0) estBkgLevel = 0;

        // Estimate Signal Amplitude (Total Height - Background Level)
        double estSigAmp = maxVal - estBkgLevel;
        if (estSigAmp < 0) estSigAmp = maxVal * 0.5; // Fallback if shape is weird

        // Set the Parameters
        preFit->SetParameter(0, estSigAmp);   // Signal Amplitude
        preFit->SetParameter(1, estimatedMu); // Mean
        preFit->SetParameter(2, 0.002);       // Sigma
        preFit->SetParameter(3, estBkgLevel); // c0 (Background constant)
        preFit->SetParameter(4, 0);           // c1
        preFit->SetParameter(5, 0);           // c2
        
        // Relax the sigma limit slightly (High stats might pick up slight resolution broadening)
        preFit->SetParLimits(2, 0.0005, 0.008); 

        // 3. Fit
        TFitResultPtr rPre = hTotalMass_SimDensity->Fit(preFit, "Q N 0 R S");

        // We check IsValid() instead of strictly == 0, as high stats sometimes gives status 4000 (Converged but non-pos-def covariance)
        // which is perfectly fine as initial guesses.
        if (rPre->IsValid()) { 
            double initGuess_sigma = preFit->GetParameter(2);
            // Conversion: Area = Amplitude * Sigma * sqrt(2*pi)
                // PerformSimultaneousFitQA expects to see an actual counts value for the initial guess, so we integrate (analytically) to actually get the counts estimate:
            double initGuess_Yield = preFit->GetParameter(0) * initGuess_sigma * std::sqrt(2 * TMath::Pi());
            double initGuess_mu    = preFit->GetParameter(1);

            // // DEBUG: print initial guesses
            // std::cout << "  [DEBUG " << extractionName << "] initGuess_Yield=" << initGuess_Yield
            //         << " initGuess_mu=" << initGuess_mu << " initGuess_sigma=" << initGuess_sigma
            //         << " hTotalMass entries=" << hTotalMass_Sim->GetEntries()
            //         << " hTotalMass integral=" << hTotalMass_Sim->Integral()
            //         << " hTotalNum Integral=" << hTotalNum_Sim->Integral() << std::endl;
            // // DEBUG: histogram addresses
            // std::cout << "  [DEBUG " << extractionName << "] h2dCounts=" << h2dCounts << " p2dRing=" << p2dRingObs << std::endl;
            // std::cout << "  [DEBUG " << extractionName << "] preFit=" << preFit  << " c0 " << preFit->GetParameter(3)
            //           << " c1 " << preFit->GetParameter(4) << " c2 " << preFit->GetParameter(5) << std::endl;

            // Extract the background polynomial from the full preFit (gaus(0)+pol2(3)).
            // par[3..5] of preFit are c0, c1, c2 of the background -- exactly what
            // PerformSimultaneousFitQA expects from a pure pol2 at par[0..2].
            // (PerformSimultaneousFitQA takes a pol2 function as input!)
            TF1* bkgInitGuessFunc = new TF1(Form("bkgInitGuess_%s", extractionName.Data()),
                                            PolShifted(cfgInt.bkgPolOrder, cfgInt.muInitGuess), massMin, massMax);
            bkgInitGuessFunc->SetParameter(0, preFit->GetParameter(3)); // c0
            bkgInitGuessFunc->SetParameter(1, preFit->GetParameter(4)); // c1
            bkgInitGuessFunc->SetParameter(2, preFit->GetParameter(5)); // c2
            // Transplant errors too, so PerformSimultaneousFitQA gets sensible step sizes
            bkgInitGuessFunc->SetParError(0, preFit->GetParError(3));
            bkgInitGuessFunc->SetParError(1, preFit->GetParError(4));
            bkgInitGuessFunc->SetParError(2, preFit->GetParError(5));

            // 4. Perform the Simultaneous Fit
                // The Simultaneous Fit expects a 
            SimFitResult simResTotal = PerformSimultaneousFitQA(
                hTotalMass_SimDensity, hTotalNum_SimDensity, massMin, massMax, // hTotalMass_SimDensity carries correct errors from h2dNum (which was built from TProfile2D)
                initGuess_Yield, initGuess_mu, initGuess_sigma, bkgInitGuessFunc,
                1e-3, 1e-3, // Guesses for background and signal
                cfgInt.nSigmaPeak
            );

            if (simResTotal.status == 0) {
                hIntegratedRSig_Sim->SetBinContent(1, simResTotal.R_S);
                hIntegratedRSig_Sim->SetBinError(1, simResTotal.err_R_S); 

                hIntegratedRBkg_Sim->SetBinContent(1, simResTotal.R_B);
                hIntegratedRBkg_Sim->SetBinError(1, simResTotal.err_R_B);

                if (printHeader) std::cout << "    [SimFit Integrated] Success. R_S: " << simResTotal.R_S << " +/- " << simResTotal.err_R_S << std::endl;
            }
            else if (printHeader) std::cout << "    [SimFit Integrated] Failed (Status " << simResTotal.status << ")" << std::endl;
        }
        else if (printHeader) std::cout << "    [SimFit Integrated] Skipped: Pre-fit failed to converge." << std::endl;
        delete preFit;
        delete hTotalMass_SimDensity;
        delete hTotalNum_SimDensity;
    }

    // 5. Write to Disk.
    // The simultaneous fit stays in ResultsCombinedFit/ and is not mirrored anywhere. It is a
    // cross-check, not a result: the joint fit assumes a CONSTANT <R>_bkg across the mass window
    // while the sideband method fits its mass dependence, so a disagreement between the two can
    // just as well be that modelling difference.
    dirResultsSim->cd();
    hIntegratedRSig_Sim->Write();
    hIntegratedRBkg_Sim->Write();

    delete hTotalMass_Sim;
    delete hTotalNum_Sim;

    // h2dNum was allocated by BuildNumFromProfile in this function. We remove it:
    delete h2dNum;
}

// A small helper function that calculates the NSigma with which each bin deviates from zero (useful to easily spot signal in the ring observable plots!)
TGraphErrors* makeSignificanceGraph(TProfile* prof, const char* name) {
    if (!prof) return nullptr;
    int nBins = prof->GetNbinsX();
    TGraphErrors* gr = new TGraphErrors();
    gr->SetName(name);
    gr->SetTitle(prof->GetTitle());
    int point = 0;
    for (int i = 1; i <= nBins; ++i) {
        double mean  = prof->GetBinContent(i);
        double error = prof->GetBinError(i);
        if (error <= 0) continue;  // skip empty bins
        double x  = prof->GetBinCenter(i);
        double ex = prof->GetBinWidth(i)/2.0; // Or 0.0, if you'd like them to have zero width
        double y  = mean / error;

        gr->SetPoint(point, x, y);
        gr->SetPointError(point, ex, 0.0);
        point++;
    }
    gr->GetXaxis()->SetTitle(prof->GetXaxis()->GetTitle());
    gr->GetYaxis()->SetTitle("<R>/#sigma");

    // Style options:
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.1);
    gr->SetLineStyle(0);   // disable lines
    gr->SetLineWidth(0);   // no line drawing
    return gr;
}

void drawSigAndDashedLines(TGraphErrors* gr, TDirectory* outDir, const std::string& canvasName) {
    if (!gr) return;
    outDir->cd();
    TCanvas* c = new TCanvas(canvasName.c_str(), "", 800, 600);
    c->SetMargin(0.12, 0.04, 0.12, 0.04);

    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.1);
    gr->SetLineWidth(2);
    gr->Draw("AP"); // Draw first to get the axes

    double xmin = gr->GetXaxis()->GetXmin();
    double xmax = gr->GetXaxis()->GetXmax();

    // ----- 1-sigma compatibility band -----
    TBox* band = new TBox(xmin, -1.0, xmax, 1.0);
    // band->SetFillColorAlpha(kGray+1, 0.25);
    band->SetFillColor(kGray);
    band->SetFillStyle(3001);   // hatched or light fill
    band->SetLineColor(0);
    band->Draw("same");

    // ----- dashed +/- 1-sigma lines -----
    TLine* lineUp = new TLine(xmin, 1.0, xmax, 1.0);
    TLine* lineDown = new TLine(xmin, -1.0, xmax, -1.0);

    lineUp->SetLineStyle(2);
    lineDown->SetLineStyle(2);
    lineUp->SetLineWidth(2);
    lineDown->SetLineWidth(2);

    lineUp->Draw("same");
    lineDown->Draw("same");
    gr->Draw("P same"); // redraw points on top

    // ----- legend -----
    TLegend* leg = new TLegend(0.60, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    leg->AddEntry(gr, "<R>/#sigma", "p");
    leg->AddEntry(band, "Compatible with 0 within 1#sigma", "f");

    leg->Draw();

    c->Write();
}

// PURPOSE:
//   For each analysis variation (Ring, RingKinematicCuts, etc.), perform signal extraction
//   on the denominator only (the invariant mass spectrum), integrated over ALL angular bins.
//   This is a broad QA that does NOT analyze any differential observable over angles.
//
//   Stored outputs:
//     1) hDenomQA_MassDensity_* : Full density histogram with a gaus+pol2 combined fit attached
//                                  (analogous to hMassDensity_BinX, but angle-integrated).
//     2) grDenomQA_Sidebands_*  : TGraphErrors of sideband points with a pol2 background fit
//                                  attached to it (visible in TBrowser as a curve through the
//                                  sideband markers). Color: green.
//     3) hDenomQA_PeakRegion_*  : Clone of the density histogram that is non-zero ONLY in
//                                  [mu - 4sigma, mu + 4sigma] (the bin-counting window).
//                                  Displayed as red hatched columns.
//     4) cDenomQA_*             : Summary TCanvas with drawing order:
//                                    a) Full density histogram (black, error bars, base layer)
//                                    b) Peak-region histogram  (red hatched columns)
//                                    c) Sideband graph markers (green open circles)
//                                    d) Background pol2 fit    (green dashed, extrapolated)
//                                    e) Combined gaus+pol2 fit (blue solid)
//
// ARGUMENTS:
//   hMassSigExtract -- TH1D* with the raw Lambda candidate counts on the mass axis,
//                      already integrated over all angular bins. This is exactly
//                      hMassSigExtract from the input file.
//   outDir          -- TDirectory where all QA output will be written.
//   qaName          -- Label string used in all histogram/canvas names to avoid ROOT name
//                      clashes. Typically the variation name ("Ring", "RingKinematicCuts"...).
//   massMin         -- Lower fit/display limit on the mass axis (GeV/c^2).
//   massMax         -- Upper fit/display limit on the mass axis (GeV/c^2).
// =================================================================================================

// ================================================================================================
// INTEGRATED EXTRACTION FROM THE 1D <R> vs MASS PROFILES
// ================================================================================================
// This is where the integrated <R>_measured, <R>_S and <R>_B actually come from.
//
// WHY NOT FROM ExtractObservable2D's STEP 8. Every angular extraction integrates over exactly the
// same candidates, so "integrated DeltaPhi" and "integrated DeltaTheta" are the same measurement
// written twice.. The three quantities that genuinely differ are the three jet-side proxies,
// because a Lambda with no subleading jet never enters the SubJet profile at all. So the integrated
// result is produced once per proxy,
// straight from the 1D profiles the consumer already writes:
//     <variation>/pRingObservableMass          -> LeadJet
//     <variation>/pRingObservableLeadPMass     -> LeadP
//     <variation>/pRingObservable2ndJetMass    -> SubJet
//
// WHY A TProfile IS THE RIGHT INPUT. A TProfile carries its own denominator: GetBinEntries(i) is
// the candidate count and GetBinContent(i)*GetBinEntries(i) is Sum_R_i over the very same
// candidates. Pairing the profile with a separate counts histogram (hMassSigExtract) would be
// wrong for LeadP and SubJet, whose candidate sets are subsets of it -- R_measured = N/T would
// then divide a numerator and a denominator built from different candidates.
//
// OPTIONAL FINE-AXIS PEAK REFERENCE. (mu, sigma) may be taken from a much finer mass spectrum
// (QA/hMass, 450 bins) rather than from the profile's own 24-bin axis. Only the peak POSITION and
// WIDTH are borrowed; every count, every sideband point and both background fits stay on the
// profile's own axis, so the candidate sets never mix. The two axes are not aligned and do not
// even cover the same range, which is exactly why nothing may be counted across them.
// ================================================================================================

struct IntegratedProxySpec {
    const char* profilePath; // Path inside the variation directory
    const char* name;        // Output folder name
    const char* label;       // For titles and logs
};

// Number of entries in a TProfile bin, promoted to a TH1D. Errors are Poisson on the count.
TH1D* CountsFromProfile(TProfile* prof, const char* name)
{
    if (!prof) return nullptr;
    const int nb = prof->GetNbinsX();
    TH1D* h = new TH1D(name, "", nb, prof->GetXaxis()->GetXbins()->GetArray()
                                       ? prof->GetXaxis()->GetXbins()->GetArray()
                                       : nullptr);
    if (!prof->GetXaxis()->GetXbins()->GetArray()) {
        delete h;
        h = new TH1D(name, "", nb, prof->GetXaxis()->GetXmin(), prof->GetXaxis()->GetXmax());
    }
    h->SetDirectory(nullptr);
    for (int i = 1; i <= nb; ++i) {
        double n = prof->GetBinEntries(i);
        h->SetBinContent(i, n);
        h->SetBinError(i, std::sqrt(std::max(n, 0.0)));
    }
    // ResetStats rather than SetEntries: it recomputes the statistics from the bin contents,
    // so GetEntries() reports the candidate count the extraction gates on, not the profile's fill count.
    h->ResetStats();
    return h;
}

// Candidate counts from a TProfile2D's own entry array, as a TH2D.
// This is what lets ExtractObservable2D run on a TProfile2D for which the consumer books no counts
// histogram at all. It is also the safer pairing in general: GetBinEntries counts exactly the
// candidates that went into the profile, so the numerator and the denominator can never come from
// different selections, which is a real risk when a separately booked TH2D is filled under a
// slightly different condition.
TH2D* CountsFromProfile2D(TProfile2D* prof, const char* name)
{
    if (!prof) return nullptr;
    const TAxis* ax = prof->GetXaxis();
    const TAxis* ay = prof->GetYaxis();
    const int nx = ax->GetNbins(), ny = ay->GetNbins();

    TH2D* h = nullptr;
    const bool varX = (ax->GetXbins()->GetSize() > 0);
    const bool varY = (ay->GetXbins()->GetSize() > 0);
    if (varX && varY)
        h = new TH2D(name, "", nx, ax->GetXbins()->GetArray(), ny, ay->GetXbins()->GetArray());
    else if (varY)
        h = new TH2D(name, "", nx, ax->GetXmin(), ax->GetXmax(), ny, ay->GetXbins()->GetArray());
    else if (varX)
        h = new TH2D(name, "", nx, ax->GetXbins()->GetArray(), ny, ay->GetXmin(), ay->GetXmax());
    else
        h = new TH2D(name, "", nx, ax->GetXmin(), ax->GetXmax(), ny, ay->GetXmin(), ay->GetXmax());
    h->SetDirectory(nullptr);

    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            double n = prof->GetBinEntries(prof->GetBin(ix, iy));
            h->SetBinContent(ix, iy, n);
            h->SetBinError(ix, iy, std::sqrt(std::max(n, 0.0)));
        }
    }
    // h->SetEntries(prof->GetEntries());
    h->ResetStats();
    return h;
}

// Sum_R_i in each bin, with the error of the SUM rather than of the mean.
// A TProfile bin error is the standard error of the mean, s/sqrt(n), so Var(Sum) = n^2 * SEM^2 at
// fixed n. That is a CONDITIONAL variance -- the count fluctuation is accounted for separately in
// ComputePeakWindowYields, which is why it must not be folded in here.
TH1D* NumeratorFromProfile(TProfile* prof, const char* name)
{
    TH1D* h = CountsFromProfile(prof, name);
    if (!h) return nullptr;
    for (int i = 1; i <= prof->GetNbinsX(); ++i) {
        double n = prof->GetBinEntries(i);
        h->SetBinContent(i, n * prof->GetBinContent(i));
        h->SetBinError(i, n * prof->GetBinError(i));
    }
    return h;
}

// Draws one annotated QA canvas: a spectrum, its fitted function, the shaded signal window, a
// stats box carrying the fit parameters, and a TPaveText carrying the literal formula.
//
// The two ordering details below are the ones that bite. gStyle has to be set BEFORE anything is
// painted, because the TPaveStats is built at paint time and reads gStyle then. And c->Update()
// has to be forced before Write(), because TCanvas::Write() does not paint: without it the canvas
// is stored with no stats object at all and ROOT rebuilds a default one when the file is reopened,
// which is how a canvas ends up with no fit information on it.
void DrawFitQACanvas(TH1* h, TF1* fit, const char* canvasName, const char* canvasTitle,
                     const char* formulaText, const std::vector<TString>& extraLines,
                     double shadeLow, double shadeHigh, TDirectory* outDir)
{
    if (!h || !outDir) return;

    // OptStat is deliberately OFF here. TPaveStats lives in the histogram's function list and is
    // carried across by Clone(), so a cloned histogram inherits the box of whatever it was cloned
    // from and c->Update() reuses it rather than rebuilding. That is how these canvases came to
    // show an entry count belonging to a different object -- billions where the profile held
    // hundreds of millions. The count is now written into the TPaveText below, from a source that
    // cannot be inherited by accident.
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(111);  // chi2/ndf, values and errors; no probability line

    TCanvas* c = new TCanvas(canvasName, canvasTitle, 900, 650);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);
    c->SetGridx();
    c->SetGridy();

    h->SetStats(1);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.8);
    h->SetLineWidth(2);
    h->Draw("PE");

    // Signal window, drawn underneath the points so it never hides them.
    if (shadeHigh > shadeLow) {
        TH1* hShade = (TH1*)h->Clone(Form("%s_shade", canvasName));
        hShade->SetDirectory(nullptr);
        for (int i = 1; i <= hShade->GetNbinsX(); ++i) {
            double x = hShade->GetBinCenter(i);
            if (x < shadeLow || x > shadeHigh) { hShade->SetBinContent(i, 0); hShade->SetBinError(i, 0); }
        }
        hShade->SetFillColorAlpha(kAzure - 9, 0.35);
        hShade->SetFillStyle(1001);
        hShade->SetLineColor(kAzure - 9);
        hShade->SetMarkerSize(0);
        hShade->Draw("HIST SAME");
        h->Draw("PE SAME");
    }

    if (fit) { fit->SetLineColor(kRed + 1); fit->SetLineWidth(2); fit->SetNpx(600); fit->Draw("SAME"); }

    c->Update(); // Force the paint pass that actually builds the TPaveStats

    TPaveStats* stats = dynamic_cast<TPaveStats*>(h->FindObject("stats"));
    if (stats) {
        // Set the options ON THE OBJECT, not only on gStyle, so they are serialised with the
        // canvas and survive whatever global style the reader happens to have active.
        stats->SetOptStat(0);
        stats->SetOptFit(111);
        stats->SetX1NDC(0.60); stats->SetX2NDC(0.90);
        stats->SetY1NDC(0.58); stats->SetY2NDC(0.90);
        stats->SetFillColor(kWhite);
        stats->SetFillStyle(1001); // Opaque, so the grid does not run through the numbers
        stats->SetBorderSize(1);
        stats->SetTextSize(0.026);
    }

    double paveTop = 0.56;
    double paveHeight = 0.045 * (1.0 + extraLines.size());
    TPaveText* pave = new TPaveText(0.42, paveTop - paveHeight, 0.90, paveTop, "NDC");
    pave->SetFillColor(kWhite);
    pave->SetFillStyle(1001);
    pave->SetBorderSize(1);
    pave->SetTextSize(0.025);
    pave->SetTextAlign(12);
    pave->AddText(formulaText);
    for (const auto& line : extraLines) pave->AddText(line.Data());
    pave->Draw();

    c->Modified();
    c->Update(); // Second pass, so the repositioned box and the label are part of what is written

    outDir->cd();
    c->Write();
    delete pave;
    delete c;
}

// ------------------------------------------------------------------------------------------------
// ExtractIntegratedFromProfile -- the integrated measurement, once per proxy.
// ------------------------------------------------------------------------------------------------
// hFinePeakRef is optional: when supplied, (mu, sigma) are taken from it instead of from the
// profile's own axis. Pass nullptr to fit the profile's axis directly.
void ExtractIntegratedFromProfile(TProfile* prof, TH1D* hFinePeakRef,
                                  TDirectory* parentDir, const IntegratedProxySpec& spec,
                                  const SidebandConfig& cfg)
{
    std::cout << Form("  [Integrated %s] ----------------------------------------------------------",
                      spec.label) << std::endl;

    TDirectory* dirProxy = EnsureDir(EnsureDir(parentDir, "IntegratedSummary"), spec.name);
    TDirectory* dirQA = EnsureDir(dirProxy, "QA");

    // hExtractionStatus is written unconditionally, so that "the extraction failed" is always
    // distinguishable from "the extraction never ran". The value histograms below are written ONLY
    // on success: a zero-filled TH1D looks exactly like a genuine measurement of zero, and that
    // ambiguity is precisely the failure mode this reporting exists to remove.
    TH1D* hStatus = new TH1D(Form("hExtractionStatus_%s", spec.name),
                             Form("%s: extraction status; ;Value", spec.label), 6, 0, 6);
    hStatus->SetDirectory(nullptr);
    hStatus->GetXaxis()->SetBinLabel(1, "valid");
    hStatus->GetXaxis()->SetBinLabel(2, "peakFitConverged");
    hStatus->GetXaxis()->SetBinLabel(3, "peakFitSalvaged");
    hStatus->GetXaxis()->SetBinLabel(4, "sidebandPoints");
    hStatus->GetXaxis()->SetBinLabel(5, "bkgPolOrder");
    hStatus->SetBinContent(5, cfg.bkgPolOrder);
    hStatus->GetXaxis()->SetBinLabel(6, "ringBkgPolOrder");
    hStatus->SetBinContent(6, cfg.ringBkgPolOrder);

    auto bail = [&](const char* why) {
        std::cout << Form("  [Integrated %s] NO RESULT: %s", spec.label, why) << std::endl;
        dirProxy->cd();
        hStatus->Write();
        delete hStatus;
    };

    if (!prof) { bail("the profile is missing from this file"); return; }

    const double massMin = prof->GetXaxis()->GetXmin();
    const double massMax = prof->GetXaxis()->GetXmax();

    TH1D* hCounts = CountsFromProfile(prof, Form("hCountsInt_%s", spec.name));
    TH1D* hNum = NumeratorFromProfile(prof, Form("hNumInt_%s", spec.name));

    std::cout << Form("    entries = %.0f over %d mass bins in [%.5f, %.5f]",
                      prof->GetEntries(), prof->GetNbinsX(), massMin, massMax) << std::endl;

    if (hCounts->GetEntries() < cfg.minEntries || hCounts->Integral() <= cfg.minIntegral) {
        bail(Form("too few candidates (%.0f entries, %.0f integral; need %.0f and %.0f)",
                  hCounts->GetEntries(), hCounts->Integral(), cfg.minEntries, cfg.minIntegral));
        delete hCounts; delete hNum;
        return;
    }

    // --- Peak position and width -------------------------------------------------------------
    // Densities, because the mass axis is variable-width: fitting counts there would fit the
    // binning as much as the physics.
    TH1D* hPeakSource = hFinePeakRef ? (TH1D*)hFinePeakRef->Clone(Form("hPeakSrc_%s", spec.name))
                                     : (TH1D*)hCounts->Clone(Form("hPeakSrc_%s", spec.name));
    hPeakSource->SetDirectory(nullptr);
    hPeakSource->Scale(1.0, "width");

    PeakFitOutcome peak = FitMassPeak(hPeakSource, hPeakSource->GetXaxis()->GetXmin(),
                                      hPeakSource->GetXaxis()->GetXmax(), cfg,
                                      Form("fitPeakInt_%s", spec.name));

    hStatus->SetBinContent(2, peak.converged ? 1 : 0);
    hStatus->SetBinContent(3, (peak.valid && !peak.converged) ? 1 : 0);

    if (!peak.valid) {
        bail(Form("the mass peak fit failed (Minuit status %d, mu = %.5f, sigma = %.5f); "
                  "neither convergence nor the salvage window accepted it",
                  peak.status, peak.mu, peak.sigma));
        delete hPeakSource; delete hCounts; delete hNum;
        return;
    }
    if (!peak.converged) {
        std::cout << Form("    NOTE: peak fit did not converge (Minuit status %d) but (mu, sigma) "
                          "landed inside the salvage window, so the result is kept.", peak.status)
                  << std::endl;
    }
    std::cout << Form("    peak: mu = %.5f, sigma = %.5f (%.3f MeV), taken from %s",
                      peak.mu, peak.sigma, peak.sigma * 1000.0,
                      hFinePeakRef ? "the fine QA mass axis" : "the profile's own axis") << std::endl;

    // --- Sidebands, always on the profile's own axis ------------------------------------------
    SidebandSpan span;
    TGraphErrors* grBkg = BuildSidebandGraph(hCounts, peak.mu, peak.sigma, massMin, massMax, cfg,
                                             SidebandGraphKind::Counts, &span,
                                             Form("grBkgInt_%s", spec.name));
    TGraphErrors* grRingBkg = BuildRingSidebandGraph(hCounts, hNum, peak.mu, peak.sigma,
                                                    massMin, massMax, cfg,
                                                    Form("grRingBkg_%s", spec.name));
    hStatus->SetBinContent(4, span.nPointsLeft + span.nPointsRight);

    std::cout << Form("    sidebands: %s -> %d left + %d right = %d points (need >= %d; counts pol%d, "
                      "<R>_bkg pol%d)",
                      SidebandBandLabel(cfg).Data(), span.nPointsLeft, span.nPointsRight,
                      grBkg->GetN(), cfg.minSidebandPoints, cfg.bkgPolOrder, cfg.ringBkgPolOrder) << std::endl;

    if (grBkg->GetN() < cfg.minSidebandPoints || grRingBkg->GetN() < cfg.minSidebandPoints) {
        bail("too few sideband points for the background polynomial");
        delete grBkg; delete grRingBkg; delete hPeakSource; delete hCounts; delete hNum;
        return;
    }

    TF1* bkgFit = new TF1(Form("bkgFitInt_%s", spec.name),
                          PolShifted(cfg.bkgPolOrder, cfg.muInitGuess), massMin, massMax);
    TF1* ringBkgFit = new TF1(Form("ringBkgFitInt_%s", spec.name),
                             PolShifted(cfg.ringBkgPolOrder, cfg.muInitGuess), massMin, massMax);
    TFitResultPtr rBkg = grBkg->Fit(bkgFit, "Q 0 S");
    TFitResultPtr rRingBkg = grRingBkg->Fit(ringBkgFit, "Q 0 S");

    if (!rBkg->IsValid() || !rRingBkg->IsValid()) {
        bail(Form("a sideband fit failed (counts status %d, numerator status %d)",
                  rBkg->Status(), rRingBkg->Status()));
        delete bkgFit; delete ringBkgFit; delete grBkg; delete grRingBkg;
        delete hPeakSource; delete hCounts; delete hNum;
        return;
    }

    PeakWindowYields y = ComputePeakWindowYields(hCounts, hNum, peak.mu, peak.sigma, cfg.nSigmaPeak,
                                                 bkgFit, rBkg, ringBkgFit, rRingBkg);
    if (!y.valid) {
        bail("the yield computation rejected the result (signal counts are non-positive)");
        delete bkgFit; delete ringBkgFit; delete grBkg; delete grRingBkg;
        delete hPeakSource; delete hCounts; delete hNum;
        return;
    }

    hStatus->SetBinContent(1, 1);

    // --- Results -------------------------------------------------------------------------------
    auto scalar = [&](const char* stem, const char* yTitle, double v, double e) {
        TH1D* h = new TH1D(Form("%s_%s", stem, spec.name),
                           Form("%s; ;%s", spec.label, yTitle), 1, 0, 1);
        h->SetDirectory(nullptr);
        h->SetBinContent(1, v);
        h->SetBinError(1, e);
        return h;
    };

    std::vector<TH1D*> outputs = {
        scalar("hIntegratedRMeas", "<R>_{measured}", y.R_peak, y.errR_peak),
        scalar("hIntegratedRSig", "<R>_{S}", y.R_S, y.errR_S),
        scalar("hIntegratedRBkg", "<R>_{B}", y.R_B, y.errR_B),
        scalar("hIntegratedDiffMeasMinusSig", "<R>_{measured} - <R>_{S}",
               y.diffPeakMinusSig, y.errDiffPeakMinusSig),
        scalar("hIntegratedDiffMeasMinusBkg", "<R>_{measured} - <R>_{B}",
               y.diffPeakMinusBkg, y.errDiffPeakMinusBkg),
        scalar("hIntegratedDiffSigMinusBkg", "<R>_{S} - <R>_{B}",
               y.diffSigMinusBkg, y.errDiffSigMinusBkg),
        scalar("hIntegratedPurity", "S/(S+B)", y.purity, y.errPurity),
        scalar("hIntegratedSignificance", "S/#sqrt{S+B}", y.significance, y.errSignificance)
    };

    TH1D* hYields = new TH1D(Form("hIntegratedYields_%s", spec.name),
                             Form("%s; ;Counts in the signal window", spec.label), 2, 0, 2);
    hYields->SetDirectory(nullptr);
    hYields->GetXaxis()->SetBinLabel(1, "Signal");
    hYields->GetXaxis()->SetBinLabel(2, "Background");
    hYields->SetBinContent(1, y.sigCounts); hYields->SetBinError(1, y.errSigCounts);
    hYields->SetBinContent(2, y.bkgCounts); hYields->SetBinError(2, y.errBkgCounts);

    TH1D* hWindow = new TH1D(Form("hAchievedWindow_%s", spec.name),
                             Form("%s: achieved extraction windows; ;Value", spec.label), 6, 0, 6);
    hWindow->SetDirectory(nullptr);
    const char* wLabels[6] = {"nSigmaLow", "nSigmaHigh", "coverage", "massLow", "massHigh", "requested"};
    double wVals[6] = {y.nSigmaAchievedLow, y.nSigmaAchievedHigh, y.coverage,
                       y.xLow, y.xHigh, cfg.nSigmaPeak};
    for (int i = 0; i < 6; ++i) { hWindow->GetXaxis()->SetBinLabel(i + 1, wLabels[i]); hWindow->SetBinContent(i + 1, wVals[i]); }

    std::cout << Form("    window: requested %.2f sigma -> achieved [-%.2f, +%.2f] sigma "
                      "(mass %.5f to %.5f), Gaussian coverage %.4f",
                      cfg.nSigmaPeak, y.nSigmaAchievedLow, y.nSigmaAchievedHigh,
                      y.xLow, y.xHigh, y.coverage) << std::endl;
    std::cout << Form("    yields: S = %.0f +/- %.0f, B = %.0f +/- %.0f, purity = %.4f, "
                      "significance = %.2f",
                      y.sigCounts, y.errSigCounts, y.bkgCounts, y.errBkgCounts,
                      y.purity, y.significance) << std::endl;
    std::cout << Form("    <R>_measured = %+.6f +/- %.6f", y.R_peak, y.errR_peak) << std::endl;
    std::cout << Form("    <R>_S        = %+.6f +/- %.6f", y.R_S, y.errR_S) << std::endl;
    std::cout << Form("    <R>_B        = %+.6f +/- %.6f", y.R_B, y.errR_B) << std::endl;
    std::cout << Form("    <R>_S - <R>_B = %+.6f +/- %.6f  (naive quadrature would give %.6f)",
                      y.diffSigMinusBkg, y.errDiffSigMinusBkg,
                      std::sqrt(y.errR_S * y.errR_S + y.errR_B * y.errR_B)) << std::endl;

    // --- QA canvases ---------------------------------------------------------------------------
    TF1* fPeakDraw = peak.fit;
    // Candidate counts stated explicitly rather than left to a TPaveStats. The sum of the profile's
    // own bin entries is the number that actually fed the extraction; TH1::GetEntries() is a fill
    // counter that can drift from it, and on a clone it may not even belong to this histogram.
    double profCandidates = 0.0;
    for (int i = 1; i <= prof->GetNbinsX(); ++i) profCandidates += prof->GetBinEntries(i);
    double peakSrcCandidates = hPeakSource->Integral();

    std::vector<TString> massLines = {
        Form("candidates in the spectrum: %.0f", peakSrcCandidates),
        Form("#mu = %.5f, #sigma = %.5f GeV/c^{2} (%.3f MeV)", peak.mu, peak.sigma, peak.sigma * 1000.0),
        Form("signal window [%.5f, %.5f] = [-%.2f, +%.2f] #sigma", y.xLow, y.xHigh,
             y.nSigmaAchievedLow, y.nSigmaAchievedHigh),
        Form("S = %.0f #pm %.0f, B = %.0f #pm %.0f, S/(S+B) = %.4f",
             y.sigCounts, y.errSigCounts, y.bkgCounts, y.errBkgCounts, y.purity)
    };
    DrawFitQACanvas(hPeakSource, fPeakDraw,
                    Form("cMassFit_%s", spec.name),
                    Form("%s: mass spectrum and peak fit", spec.label),
                    Form("[p0]*exp(-0.5*((x-[p1])/[p2])^{2}) + pol%d(3)", cfg.bkgPolOrder),
                    massLines, y.xLow, y.xHigh, dirQA);

    std::vector<TString> ringLines;
    ringLines.push_back(Form("candidates in this profile: %.0f", profCandidates));
    ringLines.push_back(Form("sideband band: %s, %d + %d points",
                             SidebandBandLabel(cfg).Data(), span.nPointsLeft, span.nPointsRight));
    for (int ip = 0; ip <= cfg.ringBkgPolOrder; ++ip)
        ringLines.push_back(Form("p%d = %+.6e #pm %.3e", ip, ringBkgFit->GetParameter(ip), ringBkgFit->GetParError(ip)));
    ringLines.push_back(Form("<R>_{meas} = %+.6f #pm %.6f", y.R_peak, y.errR_peak));
    ringLines.push_back(Form("<R>_{S} = %+.6f #pm %.6f", y.R_S, y.errR_S));
    ringLines.push_back(Form("<R>_{B} = %+.6f #pm %.6f", y.R_B, y.errR_B));

    TProfile* profDraw = (TProfile*)prof->Clone(Form("pRingVsMass_%s", spec.name));
    profDraw->SetDirectory(nullptr);
    profDraw->SetTitle(Form("%s: <R> vs mass;m_{p#pi} (GeV/c^{2});<R>", spec.label));
    // The numerator background is fitted as a DENSITY of Sum_R, so dividing by the counts density
    // is what turns it into the <R> the profile shows. Drawn over the profile's own range only.
    // <R>_bkg(m) is now fitted DIRECTLY, so this is simply that function -- no longer a ratio of
    // two polynomials, which is what used to give it poles and a shape with no physical meaning.
    TF1* fRingBkg = (TF1*)ringBkgFit->Clone(Form("fRingBkg_%s", spec.name));
    DrawFitQACanvas(profDraw, fRingBkg,
                    Form("cRingVsMass_%s", spec.name),
                    Form("%s: <R> vs mass with the sideband background", spec.label),
                    Form("<R>_{bkg}(m) = pol%d in (m - %.5f), fitted to the sidebands",
                         cfg.ringBkgPolOrder, cfg.muInitGuess),
                    ringLines, y.xLow, y.xHigh, dirQA);

    // Keep the sideband graphs and the fitted functions, so the QA canvases can be re-derived.
    dirQA->cd();
    grBkg->Write(Form("grSidebandCounts_%s", spec.name));
    grRingBkg->Write(Form("grSidebandRingObservable_%s", spec.name));
    bkgFit->Write(Form("fBkgCounts_%s", spec.name));
    ringBkgFit->Write(Form("fRingBkgFit_%s", spec.name));

    dirProxy->cd();
    for (auto* h : outputs) { h->Write(); delete h; }
    hYields->Write();
    hWindow->Write();
    hStatus->Write();

    delete hYields; delete hWindow; delete hStatus;
    delete fRingBkg; delete profDraw;
    delete bkgFit; delete ringBkgFit; delete grBkg; delete grRingBkg;
    delete hPeakSource; delete hCounts; delete hNum;
    delete peak.fit;
}

void PerformDenominatorQA(TH1D* hMassSigExtract, TDirectory* outDir,
                           const TString& qaName,
                           double massMin, double massMax,
                           const SidebandConfig& cfg)
{
    // -----------------------------------------------------------------------------------------
    // Input guard
    // -----------------------------------------------------------------------------------------
    if (!hMassSigExtract || !outDir) {
        std::cerr << "[PerformDenominatorQA] ERROR: null input pointer for " << qaName << ".\n";
        return;
    }
    std::cout << "  [DenomQA] Running denominator QA for variation: " << qaName << std::endl;
 
    // Basic statistics guard: do not attempt a fit on an empty or near-empty histogram.
    if (hMassSigExtract->GetEntries() < cfg.minEntries || hMassSigExtract->Integral() <= cfg.minIntegral) {
        std::cout << "  [DenomQA] Too few entries for " << qaName
                  << " (" << hMassSigExtract->GetEntries() << " entries). Skipping.\n";
        return;
    }
 
    // -----------------------------------------------------------------------------------------
    // Step 1: Convert raw counts to density (counts / GeV) for fitting
    // -----------------------------------------------------------------------------------------
    // hMassSigExtract already carries the full angle-integrated Lambda candidate counts --
    // it is exactly the denominator of the ring observable. We clone it and scale by bin
    // width to get dN/dM, consistent with the density convention used throughout this code.
    // The original histogram is NOT modified (it is input-file-owned).
    TH1D* hMassDensity = (TH1D*)hMassSigExtract->Clone(
        Form("hDenomQA_MassDensity_%s", qaName.Data())
    );
    hMassDensity->SetTitle(
        Form("Invariant Mass Density - %s;M_{p#pi} (GeV/c^{2});dN/dM (GeV^{-1}c^{2})", qaName.Data())
    );
    hMassDensity->Scale(1.0, "width"); // each bin content /= bin_width  (density convention)
 
    // -----------------------------------------------------------------------------------------
    // Step 2: Combined gaus + pol2 fit on the density histogram
    // -----------------------------------------------------------------------------------------
    // Initial guesses are estimated from the histogram shape (PeakInitGuessMode::ShapeEstimated),
    // which behaves much better than a raw maximum on this very background-dominated spectrum.
    // With PeakFitFallback::PdgDefaults this never "fails": if the fit does not converge, the
    // region definitions fall back to PDG-approximate values so the canvas is still meaningful.
    PeakFitOutcome peak = FitMassPeak(hMassDensity, massMin, massMax, cfg,
                                      Form("fDenomQA_Combined_%s", qaName.Data()));
    TF1* fitCombined = peak.fit;
    fitCombined->SetLineColor(kBlue);
    fitCombined->SetLineWidth(2);
    fitCombined->SetLineStyle(1); // solid

    double fitMu = peak.mu;
    double fitSigma = peak.sigma;
    bool fitOK = peak.converged;
    if (fitOK) {
        std::cout << Form("  [DenomQA] %s: mu=%.5f GeV/c^2, sigma=%.5f GeV/c^2, status=%d\n",
                          qaName.Data(), fitMu, fitSigma, peak.status);
    } else {
        std::cout << "  [DenomQA] WARNING: Combined fit did not fully converge for "
                  << qaName << " (status=" << peak.status
                  << "). Using PDG-approximate values for region definitions.\n";
    }
 
    // -----------------------------------------------------------------------------------------
    // Step 3: Sideband TGraphErrors + pol2 background fit
    // -----------------------------------------------------------------------------------------
    // Excludes [mu - 6sigma, mu + 6sigma] from the sideband fit, consistent with the 6-sigma
    // exclusion zone used by ExtractObservable2D. hMassSigExtract is passed rather than
    // hMassDensity because BuildSidebandGraph performs the bin-width scaling itself.
    TGraphErrors* grSidebands = BuildSidebandGraph(hMassSigExtract, fitMu, fitSigma,
                                                   massMin, massMax, cfg,
                                                   SidebandGraphKind::Counts, nullptr,
                                                   Form("grDenomQA_Sidebands_%s", qaName.Data()));
    grSidebands->SetTitle(
        Form("Sideband Points (%s) - %s;M_{p#pi} (GeV/c^{2});dN/dM (GeV^{-1}c^{2})",
             SidebandBandLabel(cfg).Data(), qaName.Data())
    );
 
    // Style: green open circles
    grSidebands->SetMarkerStyle(24);
    grSidebands->SetMarkerSize(1.0);
    grSidebands->SetMarkerColor(kGreen + 2);
    grSidebands->SetLineColor(kGreen + 2);
 
    // pol2 background fit (same degree as in ExtractObservable2D)
    TF1* fitBkg = new TF1(
        Form("fDenomQA_Bkg_%s", qaName.Data()),
        PolShifted(cfg.bkgPolOrder, cfg.muInitGuess), massMin, massMax
    );
    fitBkg->SetLineColor(kGreen + 2);
    fitBkg->SetLineWidth(2);
    fitBkg->SetLineStyle(2); // dashed: visually distinguishes background from combined fit
 
    bool bkgFitOK = false;
    TFitResultPtr rBkg; // Kept in scope: needed by ComputePeakWindowYields below
    if (grSidebands->GetN() >= cfg.minSidebandPoints) {
        rBkg = grSidebands->Fit(fitBkg, "Q 0 S");
        bkgFitOK = rBkg->IsValid();
        if (!bkgFitOK)
            std::cout << "  [DenomQA] WARNING: Sideband pol2 fit did not converge for "
                      << qaName << ".\n";
    } else {
        std::cout << "  [DenomQA] WARNING: Only " << grSidebands->GetN()
                  << " sideband points for " << qaName
                  << " -- skipping background fit (need >= " << cfg.minSidebandPoints << ").\n";
    }

    // -----------------------------------------------------------------------------------------
    // Step 3.5: Yields, purity and significance
    // -----------------------------------------------------------------------------------------
    // BEHAVIOUR CHANGE: this routine used to be drawing-only and threw away every number it
    // computed. Now that the extraction arithmetic lives in a shared helper, the same fits also
    // give the denominator yields for free, at no extra fitting cost.
    //
    // These numbers are reported ONLY when BOTH the peak fit and the sideband fit converged. The
    // canvas below is still drawn either way, exactly as before: with PeakFitFallback::PdgDefaults
    // a non-converged peak fit falls back to PDG-approximate values, which is good enough to
    // salvage a readable plot but NOT good enough to quote a yield from, since mu and sigma would
    // then be assumptions rather than measurements.
    if (fitOK && bkgFitOK) {
        PeakWindowYields y = ComputePeakWindowYields(hMassSigExtract, nullptr,
                                                     fitMu, fitSigma, cfg.nSigmaPeak,
                                                     fitBkg, rBkg, nullptr, TFitResultPtr());
        if (y.valid) {
            std::cout << Form("  [DenomQA] %s: S = %.1f +/- %.1f, B = %.1f +/- %.1f, "
                              "purity = %.4f +/- %.4f, S/sqrt(S+B) = %.2f +/- %.2f\n",
                              qaName.Data(), y.sigCounts, y.errSigCounts,
                              y.bkgCounts, y.errBkgCounts,
                              y.purity, y.errPurity,
                              y.significance, y.errSignificance);
        } else {
            std::cout << "  [DenomQA] WARNING: non-positive yield for " << qaName
                      << " -- no numbers reported (the QA canvas is still produced).\n";
        }
    }
 
    // Attach the background fit to the graph so it is drawn automatically in TBrowser
    grSidebands->GetListOfFunctions()->Add(fitBkg);
 
    // -----------------------------------------------------------------------------------------
    // Step 4: Peak-region histogram (non-zero only in [mu - 4sigma, mu + 4sigma])
    // -----------------------------------------------------------------------------------------
    // This gives a visual representation of the exact bin-counting integration window
    // used in Step 6 of ExtractObservable2D. Bins outside the window are zeroed out.
    TH1D* hPeakRegion = (TH1D*)hMassDensity->Clone(
        Form("hDenomQA_PeakRegion_%s", qaName.Data())
    );
    hPeakRegion->SetTitle(
        Form("Signal Region [%s] - %s;M_{p#pi} (GeV/c^{2});dN/dM (GeV^{-1}c^{2})",
             SignalWindowLabel(cfg).Data(), qaName.Data())
    );
    hPeakRegion->Reset(); // zero all bins; axis structure preserved by Clone
    if (hPeakRegion->GetSumw2N() == 0) hPeakRegion->Sumw2();
 
    int peakBinLo = hMassDensity->FindBin(fitMu - cfg.nSigmaPeak * fitSigma);
    int peakBinHi = hMassDensity->FindBin(fitMu + cfg.nSigmaPeak * fitSigma);
    for (int jBin = peakBinLo; jBin <= peakBinHi; ++jBin) {
        hPeakRegion->SetBinContent(jBin, hMassDensity->GetBinContent(jBin));
        hPeakRegion->SetBinError(jBin,   hMassDensity->GetBinError(jBin));
    }
 
    // Style: red hatched columns
    hPeakRegion->SetFillColor(kRed);
    hPeakRegion->SetFillStyle(3354); // forward-diagonal hatching
    hPeakRegion->SetLineColor(kRed);
    hPeakRegion->SetLineWidth(1);
 
    // -----------------------------------------------------------------------------------------
    // Step 5: Write individual objects to the output directory BEFORE building the canvas
    // -----------------------------------------------------------------------------------------
    // Writing before canvas creation avoids ROOT ownership ambiguity: once written, the
    // objects are safe to draw on the canvas and to delete at the end of this function.
    outDir->cd();
 
    // Attach the combined fit to the histogram so TBrowser draws it automatically on open.
    hMassDensity->GetListOfFunctions()->Add(fitCombined); // fitCombined is now ROOT-owned by hMassDensity
    hMassDensity->Write();  // Full density histogram + combined fit
    hPeakRegion->Write();   // Peak-only histogram
    grSidebands->Write();   // Sideband graph + pol2 background fit
 
    // -----------------------------------------------------------------------------------------
    // Step 6: Build summary canvas
    // -----------------------------------------------------------------------------------------
    // Drawing order chosen to avoid visual information loss when superposing with "same":
    //   a) Full density histogram (base layer, sets axis range, black with error bars)
    //   b) Peak region histogram  (red filled columns -- sits ON TOP of base in peak region)
    //   c) Sideband graph markers (green open circles -- visible in sideband regions)
    //   d) Background pol2 extrapolation (green dashed line -- runs through sidebands into peak)
    //   e) Combined gaus+pol2 fit (blue solid line -- drawn last so it is always visible)
 
    TCanvas* cQA = new TCanvas(
        Form("cDenomQA_%s", qaName.Data()),
        Form("Denominator QA - %s", qaName.Data()),
        1000, 700
    );
    cQA->SetLeftMargin(0.13);
    cQA->SetBottomMargin(0.13);
    cQA->SetRightMargin(0.05);
    cQA->SetTopMargin(0.07);
 
    // -- a) Full density histogram --
    // SetStats(0): suppress the statistics box (it clutters the QA canvas)
    // The fit function was just attached above; drawing the histogram will also draw the fit.
    // We suppress that here so we can control Z-order manually with explicit Draw calls below.
    hMassDensity->SetStats(0);
    hMassDensity->SetLineColor(kBlack);
    hMassDensity->SetLineWidth(2);
    hMassDensity->SetMarkerStyle(20);
    hMassDensity->SetMarkerSize(0.7);
    hMassDensity->SetMarkerColor(kBlack);
    // Draw without the attached fit function to maintain our explicit Z-order.
    // We use "HIST E" to show error bars but suppress the function list auto-draw by
    // temporarily disabling the function painting; then redraw manually below.
    // Simpler approach: draw with "E FUNC" at the end -- ROOT will draw data then functions.
    // We want functions drawn AFTER the peak histogram and sideband graph, so we use "E" first.
    hMassDensity->Draw("E"); // This sets axis range; fit function is drawn automatically AFTER
 
    // Temporarily remove the fit from the function list so ROOT doesn't auto-draw it yet.
    // We will draw it manually after the other objects, then re-add it.
    hMassDensity->GetListOfFunctions()->Remove(fitCombined);
 
    // Re-draw the histogram without attached functions (clean base layer)
    hMassDensity->Draw("E");
 
    // -- b) Peak-region histogram (red filled columns) --
    hPeakRegion->Draw("HIST same"); // "HIST" forces column drawing; "same" overlays on base
 
    // -- c) Sideband graph (green open circle markers) --
    if (grSidebands->GetN() > 0) {
        // "P same" = draw markers only on the same pad; attached fit function draws automatically
        // after the markers. To preserve Z-order, we remove the fit from the graph's list
        // here too and draw it explicitly below.
        grSidebands->GetListOfFunctions()->Remove(fitBkg);
        grSidebands->Draw("P same");
 
        // -- d) Background pol2 extrapolation (green dashed) --
        if (bkgFitOK) fitBkg->Draw("same");
 
        // Re-attach for TBrowser (doesn't affect the already-rendered canvas)
        grSidebands->GetListOfFunctions()->Add(fitBkg);
    }
 
    // -- e) Combined gaus+pol2 fit (blue solid, topmost layer) --
    if (fitOK) fitCombined->Draw("same");
 
    // Re-attach the combined fit to the histogram (for TBrowser; Write() already captured it)
    hMassDensity->GetListOfFunctions()->Add(fitCombined);
 
    // -----------------------------------------------------------------------------------------
    // Legend
    // -----------------------------------------------------------------------------------------
    TLegend* leg = new TLegend(0.54, 0.57, 0.92, 0.91);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.033);
    leg->AddEntry(hMassDensity, "Invariant mass", "lep");
    leg->AddEntry(hPeakRegion, Form("Signal region [%s]", SignalWindowLabel(cfg).Data()), "f");
    leg->AddEntry(grSidebands, Form("Sideband points (%s)", SidebandBandLabel(cfg).Data()), "p");
    leg->AddEntry(fitBkg, "Background pol2 fit", "l");
    leg->AddEntry(fitCombined, "Signal + bkg fit (gaus+pol2)", "l");
    leg->Draw();
 
    // Write the canvas to the output directory
    outDir->cd();
    cQA->Write();
 
    // -----------------------------------------------------------------------------------------
    // Cleanup
    // -----------------------------------------------------------------------------------------
    // fitCombined: owned by hMassDensity via GetListOfFunctions() -- deleted with hMassDensity.
    // fitBkg:      owned by grSidebands via GetListOfFunctions() -- deleted with grSidebands.
    // cQA: does not own the drawn objects (only holds pointers); safe to delete independently.
    delete hMassDensity;  // also deletes fitCombined
    delete hPeakRegion;
    delete grSidebands;   // also deletes fitBkg
    delete cQA;
 
    std::cout << "  [DenomQA] Done for " << qaName << ". Objects written to "
              << outDir->GetPath() << std::endl;
}

// ================================================================================================
// Helper: PerformSelectionCutFlowExtraction
// ================================================================================================
// PURPOSE:
//   Runs the standard bin-counting + sideband-subtraction signal extraction ONCE PER SELECTION
//   STEP of the V0 selection flow, using the "selection flow vs invariant mass" TH2D produced by
//   the TableProducer task (lambdajetpolarizationions/GeneralQA/h2dSelectionLambdaMass and its AntiLambda counterpart).
//
//   The point is to see how the estimated signal and background evolve along the cuts: a cut that
//   removes a lot of background and little signal is doing its job; a cut that removes both in the
//   same proportion is only costing statistics.
//
//   The per-step procedure is deliberately identical to Steps 4-6 of ExtractObservable2D:
//     1. Project the mass spectrum of that selection step.
//     2. Fit gaus + pol2 on the density spectrum to locate the peak (mu, sigma).
//     3. Refit the background alone on the sidebands (6 sigma exclusion zone) with a pol2.
//     4. Count in the [mu - 4 sigma, mu + 4 sigma] window and subtract the integrated background.
//
// ARGUMENTS:
//   h2dSelectionMass -- TH2D with the labelled selection axis on X and the invariant mass on Y.
//   outDir           -- TDirectory in which the results are stored.
//   tag              -- short label used to build unique ROOT names ("Lambda", "AntiLambda").
//   hypothesisTitle  -- mass axis title for the QA plots.
//   chainSegments    -- inclusive ROOT bin ranges, listed in nesting order, that make up this
//                       species' selection chain. Bins outside these ranges belong to the other
//                       hypothesis' branch and are skipped entirely: they are neither extracted
//                       nor used as a reference for any retention ratio. Consecutive segments are
//                       treated as consecutive links of one nested chain, so the step-to-step
//                       ratio at the start of a segment correctly refers back to the last bin of
//                       the previous segment rather than to the adjacent (unrelated) bin.
//
// NOTE ON BINNING AND LABELS:
//   Nothing about the selection axis is hardcoded here. The number of steps, the bin edges and
//   the custom cut labels are all inherited from the input histogram: ProjectionX copies the
//   alphanumeric label list, so the cut names stay automatically in sync with the task, including
//   the greyed-out "(off)" markers for selections that were disabled in that run.
// ================================================================================================
void PerformSelectionCutFlowExtraction(TH2D* h2dSelectionMass,
                                       TDirectory* outDir,
                                       TString tag,
                                       TString hypothesisTitle,
                                       const std::vector<std::pair<int, int>>& chainSegments,
                                       const SidebandConfig& cfg)
{
    if (!h2dSelectionMass || !outDir) {
        std::cerr << "[CutFlow] ERROR: null input for tag '" << tag << "'. Skipping.\n";
        return;
    }

    // Mass window and selection binning both come straight from the input histogram
    const double massMin = h2dSelectionMass->GetYaxis()->GetXmin();
    const double massMax = h2dSelectionMass->GetYaxis()->GetXmax();
    const int nCuts = h2dSelectionMass->GetNbinsX();

    std::cout << "  -> [CutFlow " << tag << "] " << nCuts << " selection steps, mass window ["
              << massMin << ", " << massMax << "] GeV/c^2" << std::endl;

    // Per-step mass fits live in their own subfolder, to keep the results folder readable
    TDirectory* dirCutFits = outDir->mkdir(Form("MassFits_%s", tag.Data()));

    // ---------------------------------------------------------------------------------------
    // Output histograms
    // ---------------------------------------------------------------------------------------
    // The first one is built by projecting the selection axis (which carries the labels) and
    // then reset; all the others are clones of it, so every output shares the same labelled axis.
    TH1D* hTot = (TH1D*)h2dSelectionMass->ProjectionX(Form("hCutFlowTotal_%s", tag.Data()));
    hTot->SetDirectory(nullptr); // Owned by this function until it is written below
    hTot->Reset();
    if (hTot->GetSumw2N() == 0) hTot->Sumw2();
    hTot->SetTitle(Form("%s: counts in the [%s] window vs selection step; ;Counts",
                        tag.Data(), SignalWindowLabel(cfg).Data()));

    // Small helper so the clone + reset + retitle sequence is written only once
    auto bookLike = [&](const char* name, const char* title) {
        TH1D* h = (TH1D*)hTot->Clone(name);
        h->SetDirectory(nullptr);
        h->Reset();
        if (h->GetSumw2N() == 0) h->Sumw2();
        h->SetTitle(title);
        return h;
    };

    TH1D* hSig = bookLike(Form("hCutFlowSignal_%s", tag.Data()),
                          Form("%s: extracted signal vs selection step; ;Signal counts", tag.Data()));
    TH1D* hBkg = bookLike(Form("hCutFlowBackground_%s", tag.Data()),
                          Form("%s: estimated background vs selection step; ;Background counts", tag.Data()));
    TH1D* hSoverB = bookLike(Form("hCutFlowSoverB_%s", tag.Data()),
                             Form("%s: S/B vs selection step; ;S/B", tag.Data()));
    TH1D* hPurity = bookLike(Form("hCutFlowPurity_%s", tag.Data()),
                             Form("%s: purity S/(S+B) vs selection step; ;Purity", tag.Data()));
    TH1D* hSignificance = bookLike(Form("hCutFlowSignificance_%s", tag.Data()),
                                   Form("%s: significance S/#sqrt{S+B} vs selection step; ;Significance", tag.Data()));
    // Cumulative retention, relative to the first successfully extracted step
    TH1D* hSigCum = bookLike(Form("hCutFlowSignalRetention_%s", tag.Data()),
                             Form("%s: cumulative signal retention; ;S_{i} / S_{ref}", tag.Data()));
    TH1D* hBkgCum = bookLike(Form("hCutFlowBackgroundRetention_%s", tag.Data()),
                             Form("%s: cumulative background retention; ;B_{i} / B_{ref}", tag.Data()));
    // Step-to-step retention: the marginal cost/benefit of each individual cut
    TH1D* hSigStep = bookLike(Form("hCutFlowSignalStep_%s", tag.Data()),
                              Form("%s: step-to-step signal retention; ;S_{i} / S_{i-1}", tag.Data()));
    TH1D* hBkgStep = bookLike(Form("hCutFlowBackgroundStep_%s", tag.Data()),
                              Form("%s: step-to-step background retention; ;B_{i} / B_{i-1}", tag.Data()));

    // Buffers, so the ratio histograms can be filled after the extraction loop
    std::vector<double> sigVal(nCuts + 1, 0.0), sigErr(nCuts + 1, 0.0);
    std::vector<double> bkgVal(nCuts + 1, 0.0), bkgErr(nCuts + 1, 0.0);
    std::vector<bool> stepOk(nCuts + 1, false);

    // Flatten the segments into a single ordered list of bins to visit. This list IS the nesting
    // chain: consecutive entries are nested samples even where the bin numbers jump across the
    // other hypothesis' block.
    std::vector<int> chainBins;
    for (const auto& seg : chainSegments) {
        for (int b = std::max(1, seg.first); b <= std::min(nCuts, seg.second); ++b) chainBins.push_back(b);
    }
    if (chainBins.empty()) {
        std::cerr << "[CutFlow] ERROR: empty selection chain for tag '" << tag << "'. Skipping.\n";
        return;
    }
    std::cout << "  -> [CutFlow " << tag << "] chain covers " << chainBins.size()
              << " of the " << nCuts << " bins (bins outside this species' branch are skipped)."
              << std::endl;

    // ---------------------------------------------------------------------------------------
    // Extraction loop, one selection step at a time (in chain order)
    // ---------------------------------------------------------------------------------------
    int nGood = 0;
    for (int iCut : chainBins) {

        const char* cutLabel = h2dSelectionMass->GetXaxis()->GetBinLabel(iCut);

        TH1D* hMassProj = h2dSelectionMass->ProjectionY(Form("hMass_%s_Cut%02d", tag.Data(), iCut), iCut, iCut, "e");
        hMassProj->SetDirectory(nullptr); // Owned here, deleted before leaving this iteration
        hMassProj->SetTitle(Form("%s step %d: %s;%s;Counts", tag.Data(), iCut, cutLabel, hypothesisTitle.Data()));

        // --- STABILITY CHECK: is there anything to fit? ---
        // Disabled selections are typically empty bins, so they fall out here naturally and stay
        // at zero in the outputs.
        if (hMassProj->GetEntries() < cfg.minEntries || hMassProj->Integral() <= cfg.minIntegral) {
            delete hMassProj;
            continue;
        }

        // --- Peak location from a gaus + pol2 fit on the density spectrum ---
        TH1D* hMassDensity = (TH1D*)hMassProj->Clone(Form("hMassDensity_%s_Cut%02d", tag.Data(), iCut));
        hMassDensity->SetDirectory(nullptr);
        hMassDensity->Scale(1.0, "width");
        hMassDensity->SetTitle(Form("%s step %d: %s;%s;dN/dM", tag.Data(), iCut, cutLabel, hypothesisTitle.Data()));

        PeakFitOutcome peak = FitMassPeak(hMassDensity, massMin, massMax, cfg,
                                          Form("fit_%s_Cut%02d", tag.Data(), iCut));

        dirCutFits->cd();
        hMassDensity->GetListOfFunctions()->Add(peak.fit); // Histogram now owns the TF1
        hMassDensity->Write();

        if (!peak.valid) {
            delete hMassDensity; // Also deletes peak.fit, which it owns
            delete hMassProj;
            continue;
        }
        double mu = peak.mu;
        double sigma = peak.sigma;

        // --- Sideband background fit, 6 sigma exclusion zone ---
        SidebandSpan spanBkg; // Reports the raw counts used by the stability check below
        TGraphErrors* grBkg = BuildSidebandGraph(hMassProj, mu, sigma, massMin, massMax, cfg,
                                                 SidebandGraphKind::Counts, &spanBkg,
                                                 Form("grBkg_%s_Cut%02d", tag.Data(), iCut));
        grBkg->SetTitle(Form("%s step %d sidebands: %s;%s;Counts/BinWidth", tag.Data(), iCut, cutLabel, hypothesisTitle.Data()));

        // --- STABILITY CHECKS: enough points for a pol2, and non-empty sidebands ---
        if (grBkg->GetN() < cfg.minSidebandPoints ||
            (cfg.minSidebandCounts > 0.0 && spanBkg.rawCounts <= cfg.minSidebandCounts)) {
            delete grBkg;
            delete hMassDensity;
            delete hMassProj;
            continue;
        }

        TF1* bkgFitFunc = new TF1(Form("bkgFit_%s_Cut%02d", tag.Data(), iCut),
                                  PolShifted(cfg.bkgPolOrder, cfg.muInitGuess), massMin, massMax);
        TFitResultPtr rBkg = grBkg->Fit(bkgFitFunc, "Q 0 S");

        if (!rBkg->IsValid()) {
            delete bkgFitFunc; // Not yet owned by the graph, so we still delete it ourselves
            delete grBkg;
            delete hMassDensity;
            delete hMassProj;
            continue;
        }

        dirCutFits->cd();
        grBkg->GetListOfFunctions()->Add(bkgFitFunc); // Graph now owns bkgFitFunc
        grBkg->Write();

        // --- Signal extraction (denominator only: there is no ring observable here) ---
        PeakWindowYields y = ComputePeakWindowYields(hMassProj, nullptr, mu, sigma, cfg.nSigmaPeak,
                                                     bkgFitFunc, rBkg, nullptr, TFitResultPtr());

        if (!y.valid) {
            std::cout << "    [CutFlow " << tag << "] step " << iCut << ": non-positive yield, skipped.\n";
            delete grBkg;
            delete hMassDensity;
            delete hMassProj;
            continue;
        }

        // --- Fill ---
        hTot->SetBinContent(iCut, y.totCounts);
        hTot->SetBinError(iCut, std::sqrt(y.totCountsErrSq));

        hSig->SetBinContent(iCut, y.sigCounts);
        hSig->SetBinError(iCut, y.errSigCounts);

        hBkg->SetBinContent(iCut, y.bkgCounts);
        hBkg->SetBinError(iCut, y.errBkgCounts);

        // S/B is deliberately propagated as T/B - 1 rather than as a ratio of S and B.
        // S = T - B, so S and B are strongly anti-correlated and treating them as independent
        // would be wrong. In terms of T (peak-region counts) and B (from a fit to the sidebands
        // only, hence independent of T) the transformation is exact:
        //   S/B = T/B - 1   =>   var(S/B) = (T/B)^2 * [ var(T)/T^2 + var(B)/B^2 ]
        double sOverB = y.sigCounts / y.bkgCounts;
        double errSOverB = (y.totCounts / y.bkgCounts) * std::sqrt(y.totCountsErrSq / std::pow(y.totCounts, 2)
                                                                   + std::pow(y.errBkgCounts / y.bkgCounts, 2));
        hSoverB->SetBinContent(iCut, sOverB);
        hSoverB->SetBinError(iCut, errSOverB);

        hPurity->SetBinContent(iCut, y.purity);
        hPurity->SetBinError(iCut, y.errPurity);

        hSignificance->SetBinContent(iCut, y.significance);
        hSignificance->SetBinError(iCut, y.errSignificance);

        sigVal[iCut] = y.sigCounts;  sigErr[iCut] = y.errSigCounts;
        bkgVal[iCut] = y.bkgCounts;  bkgErr[iCut] = y.errBkgCounts;
        stepOk[iCut] = true;
        nGood++;

        delete grBkg;        // Also deletes bkgFitFunc, which it owns
        delete hMassDensity; // Also deletes peak.fit, which it owns
        delete hMassProj;
    }

    // ---------------------------------------------------------------------------------------
    // Retention fractions
    // ---------------------------------------------------------------------------------------
    // These are ratios of NESTED yields: the candidates surviving the later step are a subset of
    // those surviving the earlier one. Propagating them as if numerator and denominator were
    // independent badly overestimates the error -- a cut that removes nothing would come out as
    // 1.0 +/- sqrt(2)*rel instead of the correct 1.0 +/- 0.0.
    //
    // DERIVATION. Write each yield as S = T - B, with T the peak-window counts and B the
    // integrated sideband background. Because the later sample is a subset of the earlier one,
    // the earlier one splits into two DISJOINT pieces,
    //     T_den = T_num + T_rejected,
    // and disjoint Poisson populations are independent, so Cov(T_num, T_den) = Var(T_num).
    // The same argument applies to the sideband counts that drive B. The peak region and the
    // sidebands are disjoint in mass, so the cross terms Cov(T, B) between the two steps vanish.
    // Adding up:
    //     Cov(S_num, S_den) = Var(T_num) + Var(B_num) = Var(S_num).
    // The usual ratio expansion for r = num/den then gives
    //     var(r)/r^2 = var_num/num^2 + var_den/den^2 - 2*var_num/(num*den).
    //
    // CONSISTENCY CHECK. In the pure counting limit (no background, var = S) this collapses to
    //     var(r) = eps*(1 - eps)/S_den,
    // which is exactly the binomial efficiency variance. The expression below is therefore the
    // background-subtracted generalization of the standard binomial treatment. No bootstrap is
    // needed: a bootstrap would require per-candidate cut-survival information, which the 2D
    // selection-flow histogram does not carry, and under the same assumptions it would converge
    // to this same result anyway.
    //
    // CAVEAT. The peak window is redefined at every step, since mu and sigma are refitted. The
    // nesting argument strictly holds for a fixed mass window, so this is exact only to the
    // extent that the fitted peak position and width are stable along the selection flow.
    auto nestedRatioWithError = [](double num, double errNum, double den, double errDen,
                                   double& val, double& err) {
        val = num / den;
        double varNum = errNum * errNum;
        double varDen = errDen * errDen;
        double relVar = varNum / (num * num) + varDen / (den * den) - 2.0 * varNum / (num * den);
        // relVar is non-negative whenever var_den >= var_num, which is what nesting implies.
        // Background-fit fluctuations can marginally break that, so clamp rather than take the
        // square root of a negative number.
        if (relVar < 0.0) relVar = 0.0;
        err = std::fabs(val) * std::sqrt(relVar);
    };

    int refIdx = -1;
    for (int b : chainBins) {
        if (stepOk[b]) { refIdx = b; break; }
    }

    int prevIdx = -1;
    for (int i : chainBins) {
        if (!stepOk[i]) continue;

        double val = 0.0, err = 0.0;

        // Cumulative, relative to the first step that could be extracted.
        // For i == refIdx the formula above returns exactly 1.0 +/- 0.0 on its own, so there is
        // no need to special-case the reference step.
        if (refIdx > 0) {
            nestedRatioWithError(sigVal[i], sigErr[i], sigVal[refIdx], sigErr[refIdx], val, err);
            hSigCum->SetBinContent(i, val); hSigCum->SetBinError(i, err);
            nestedRatioWithError(bkgVal[i], bkgErr[i], bkgVal[refIdx], bkgErr[refIdx], val, err);
            hBkgCum->SetBinContent(i, val); hBkgCum->SetBinError(i, err);
        }

        // Step to step, relative to the previous step that could be extracted
        if (prevIdx > 0) {
            nestedRatioWithError(sigVal[i], sigErr[i], sigVal[prevIdx], sigErr[prevIdx], val, err);
            hSigStep->SetBinContent(i, val); hSigStep->SetBinError(i, err);
            nestedRatioWithError(bkgVal[i], bkgErr[i], bkgVal[prevIdx], bkgErr[prevIdx], val, err);
            hBkgStep->SetBinContent(i, val); hBkgStep->SetBinError(i, err);
        }

        prevIdx = i;
    }

    // ---------------------------------------------------------------------------------------
    // Write everything out
    // ---------------------------------------------------------------------------------------
    std::vector<TH1D*> outputs = {hTot, hSig, hBkg, hSoverB, hPurity, hSignificance,
                                  hSigCum, hBkgCum, hSigStep, hBkgStep};

    outDir->cd();
    for (auto* h : outputs) {
        h->GetXaxis()->LabelsOption("v"); // Vertical labels: the selection names are long
        h->SetStats(0);
        h->Write();
    }

    // --- Summary canvas: signal and background side by side along the selection flow ---
    TCanvas* cYields = new TCanvas(Form("cCutFlowYields_%s", tag.Data()),
                                   Form("Signal and background vs selection step (%s)", tag.Data()), 1400, 800);
    cYields->SetBottomMargin(0.42); // The selection labels are long, so they need the room
    cYields->SetLogy();
    cYields->SetGridx(); // NOTE: the two branches disagreed here (SetGridy in the other one)

    hBkg->SetLineColor(kRed + 1);
    hBkg->SetMarkerColor(kRed + 1);
    hBkg->SetMarkerStyle(21);
    hSig->SetLineColor(kBlue + 1);
    hSig->SetMarkerColor(kBlue + 1);
    hSig->SetMarkerStyle(20);

    // With "SAME", the frame is defined ENTIRELY by the first histogram drawn. Since the background
    // is far above the signal, letting hBkg set the range pushed the signal points below the bottom
    // of the log-scale pad -- they were drawn, just outside the visible frame. So compute the range
    // over BOTH histograms explicitly. Only strictly positive values are considered, because a log
    // axis cannot show zeros (skipped steps sit at exactly zero and are simply not drawn).
    double yLow = 0.0, yHigh = 0.0;
    for (const TH1D* h : {hSig, hBkg}) {
        for (int b = 1; b <= h->GetNbinsX(); ++b) {
            double c = h->GetBinContent(b);
            if (c <= 0.0) continue;
            double e = h->GetBinError(b);
            double lo = (c - e > 0.0) ? (c - e) : c; // Keep it positive for the log axis
            if (yLow == 0.0 || lo < yLow) yLow = lo;
            if (c + e > yHigh) yHigh = c + e;
        }
    }
    if (yLow <= 0.0 || yHigh <= 0.0) { yLow = 0.5; yHigh = 10.0; } // Degenerate fallback
    hBkg->SetMinimum(0.5 * yLow);  // A little breathing room below
    hBkg->SetMaximum(5.0 * yHigh); // Extra headroom above, so the legend does not sit on the points

    hBkg->SetTitle(Form("%s: signal and background vs selection step; ;Counts in the [%s] window",
                        tag.Data(), SignalWindowLabel(cfg).Data()));
    hBkg->Draw("E1");
    hSig->Draw("E1 SAME");

    TLegend* legYields = new TLegend(0.68, 0.78, 0.89, 0.89);
    legYields->AddEntry(hSig, "Extracted signal", "lep");
    legYields->AddEntry(hBkg, "Estimated background", "lep");
    legYields->Draw();

    outDir->cd();
    cYields->Write();

    std::cout << "  -> [CutFlow " << tag << "] extraction succeeded for " << nGood
              << " of " << nCuts << " selection steps." << std::endl;

    // Cleanup: the canvas goes first, so nothing it still points at is deleted underneath it
    delete cYields;
    delete legYields;
    for (auto* h : outputs) delete h;
}


// ------------------------------------------------------------------------------------------------
// Main Macro
// ------------------------------------------------------------------------------------------------
// ================================================================================================
// ================================================================================================
// COMMAND LINE OPTIONS
// ================================================================================================
// Two positional arguments as before, so existing callers (run_all_wagons.sh) keep working
// untouched, plus optional --key=value flags for the extraction windows and the background model.
// Everything defaults to the historical behaviour, so a run with no flags is comparable with older
// output; see the TODO: block on MakeSharedSidebandConfig for where the defaults are headed.
// ================================================================================================
void PrintUsage(const char* exeName)
{
    SidebandConfig d = MakeSharedSidebandConfig();
    ResolveDerivedFields(d);
    std::cout
        << "Usage: " << exeName << " <inputFilePath> <outputFolderPath> [options]\n\n"
        << "Options (all optional; the value shown is the current default):\n"
        << "  --nSigmaPeak=<x>            Signal window half-width, in fitted sigmas   [" << d.nSigmaPeak << "]\n"
        << "  --nSigmaExclusion=<x>       Inner sideband edge, in fitted sigmas        [" << d.nSigmaExclusion << "]\n"
        << "  --nSigmaExclusionOuter=<x>  Outer sideband edge, in fitted sigmas;\n"
        << "                              <= 0 means the edges of the mass axis        [" << d.nSigmaExclusionOuter << "]\n"
        << "  --bkgPolOrder=<n>           Order of the polN fitted to the counts\n"
        << "                              sideband                                     [" << d.bkgPolOrder << "]\n"
        << "  --ringBkgPolOrder=<n>       Order of the polN fitted to <R>_bkg(m) in\n"
        << "                              the sidebands (0 = flat)                     [" << d.ringBkgPolOrder << "]\n"
        << "  --minSidebandPoints=<n>     Points required before fitting a sideband;\n"
        << "                              <= 0 derives it as bkgPolOrder + 3           [" << d.minSidebandPoints << "]\n"
        << "  --minSidebandCounts=<x>     Raw counts required in the sidebands         [" << d.minSidebandCounts << "]\n"
        << "  --minEntries=<x>            Entries required before fitting a spectrum   [" << d.minEntries << "]\n"
        << "  --minIntegral=<x>           Integral required before fitting a spectrum  [" << d.minIntegral << "]\n"
        << "  --muInitGuess=<x>           Initial guess for the peak position, GeV/c^2 [" << d.muInitGuess << "]\n"
        << "  --sigmaInitGuess=<x>        Initial guess for the peak width, GeV/c^2    [" << d.sigmaInitGuess << "]\n"
        << "  --help                      Print this message and exit\n\n"
        << "The windows are applied identically to all four workflows (per-bin, integrated,\n"
        << "denominator QA and cut flow). Those four differ only in what they do with a mass fit\n"
        << "that fails to converge, which is not something a command line should be choosing.\n";
}

// Applies the parsed overrides to every workflow config at once. The window definitions are
// physics and have no business differing between workflows; only the failure policy does.
bool ParseCommandLine(int argc, char** argv, ExtractionConfigSet& configs, bool& wantsHelp)
{
    wantsHelp = false;

    auto applyToAll = [&configs](void (*setter)(SidebandConfig&, double), double v) {
        setter(configs.perBin, v);
        setter(configs.integrated, v);
        setter(configs.denomQA, v);
        setter(configs.cutFlow, v);
    };

    for (int i = 3; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") { wantsHelp = true; return true; }

        size_t eq = arg.find('=');
        if (arg.rfind("--", 0) != 0 || eq == std::string::npos) {
            std::cerr << "  Error: unrecognized argument '" << arg << "'. Expected --key=value.\n";
            return false;
        }
        std::string key = arg.substr(2, eq - 2);
        std::string valStr = arg.substr(eq + 1);

        double val = 0.0;
        try { val = std::stod(valStr); }
        catch (const std::exception&) {
            std::cerr << "  Error: could not read a number from '" << arg << "'.\n";
            return false;
        }

        if      (key == "nSigmaPeak")           applyToAll([](SidebandConfig& c, double v){ c.nSigmaPeak = v; }, val);
        else if (key == "nSigmaExclusion")      applyToAll([](SidebandConfig& c, double v){ c.nSigmaExclusion = v; }, val);
        else if (key == "nSigmaExclusionOuter") applyToAll([](SidebandConfig& c, double v){ c.nSigmaExclusionOuter = v; }, val);
        else if (key == "bkgPolOrder")          applyToAll([](SidebandConfig& c, double v){ c.bkgPolOrder = int(v); c.minSidebandPoints = 0; }, val);
        else if (key == "ringBkgPolOrder")      applyToAll([](SidebandConfig& c, double v){ c.ringBkgPolOrder = int(v); c.minSidebandPoints = 0; }, val);
        else if (key == "minSidebandPoints")    applyToAll([](SidebandConfig& c, double v){ c.minSidebandPoints = int(v); }, val);
        else if (key == "minSidebandCounts")    applyToAll([](SidebandConfig& c, double v){ c.minSidebandCounts = v; }, val);
        else if (key == "minEntries")           applyToAll([](SidebandConfig& c, double v){ c.minEntries = v; }, val);
        else if (key == "minIntegral")          applyToAll([](SidebandConfig& c, double v){ c.minIntegral = v; }, val);
        else if (key == "muInitGuess")          applyToAll([](SidebandConfig& c, double v){ c.muInitGuess = v; }, val);
        else if (key == "sigmaInitGuess")       applyToAll([](SidebandConfig& c, double v){ c.sigmaInitGuess = v; }, val);
        else {
            std::cerr << "  Error: unknown option '--" << key << "'. Run with --help for the list.\n";
            return false;
        }
    }

    // bkgPolOrder resets minSidebandPoints to the derived value, so resolve after the whole parse.
    ResolveDerivedFields(configs.perBin);
    ResolveDerivedFields(configs.integrated);
    ResolveDerivedFields(configs.denomQA);
    ResolveDerivedFields(configs.cutFlow);

    if (configs.perBin.nSigmaExclusion <= configs.perBin.nSigmaPeak) {
        std::cerr << Form("  Error: nSigmaExclusion (%.3g) must exceed nSigmaPeak (%.3g), otherwise "
                          "the sidebands overlap the signal window.\n",
                          configs.perBin.nSigmaExclusion, configs.perBin.nSigmaPeak);
        return false;
    }
    if (configs.perBin.nSigmaExclusionOuter > 0.0 &&
        configs.perBin.nSigmaExclusionOuter <= configs.perBin.nSigmaExclusion) {
        std::cerr << Form("  Error: nSigmaExclusionOuter (%.3g) must exceed nSigmaExclusion (%.3g).\n",
                          configs.perBin.nSigmaExclusionOuter, configs.perBin.nSigmaExclusion);
        return false;
    }
    return true;
}

// ------------------------------------------------------------------------------------------------
// CheckSidebandBandIsFittable -- refuse to start a run that cannot possibly produce a result.
// ------------------------------------------------------------------------------------------------
// The requested band is resolved once against the reference mass spectrum, using the PDG mass and
// the sigma initial guess in place of a per-bin fit. That is approximate, but it is more than
// enough to catch the failure this guards against: on a coarse mass axis a bounded sideband band
// can contain fewer points than the polN needs, and the symptom is 144 extraction calls that all
// silently invalidate every bin and write out empty histograms. One line at startup is cheaper
// than finding that in the output.
bool CheckSidebandBandIsFittable(TH1D* hRefMass, const SidebandConfig& cfg)
{
    if (!hRefMass) return true; // Nothing to check against; let the run proceed

    double mu = cfg.muInitGuess;
    double sigma = cfg.sigmaInitGuess;
    double massMin = hRefMass->GetXaxis()->GetXmin();
    double massMax = hRefMass->GetXaxis()->GetXmax();

    SidebandSpan span;
    TGraphErrors* gr = BuildSidebandGraph(hRefMass, mu, sigma, massMin, massMax, cfg,
                                          SidebandGraphKind::Counts, &span, "grSidebandPrecheck");
    int nPoints = gr ? gr->GetN() : 0;
    delete gr;

    std::cout << Form("  Sideband precheck: band %s around mu = %.5f with sigma = %.5f gives "
                      "%d + %d = %d points (need >= %d for a pol%d).",
                      SidebandBandLabel(cfg).Data(), mu, sigma,
                      span.nPointsLeft, span.nPointsRight, nPoints,
                      cfg.minSidebandPoints, cfg.bkgPolOrder)
              << std::endl;

    if (nPoints < cfg.minSidebandPoints) {
        std::cerr << "  Error: the requested sideband band cannot be fitted on this mass axis. "
                     "Every bin would be invalidated.\n"
                  << "         Either widen the band (raise --nSigmaExclusionOuter, or set it <= 0 "
                     "to use the axis edges),\n"
                  << "         lower --bkgPolOrder, or refine axisLambdaMassSigExtract in the "
                     "consumer and re-run it.\n";
        return false;
    }
    return true;
}

int main(int argc, char** argv) {
    if (argc >= 2 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) {
        PrintUsage(argv[0]);
        return 0;
    }
    if (argc < 3) {
        PrintUsage(argv[0]);
        return 1;
    }

    std::string inputFilePath = argv[1];
    std::string outputFolderPath = argv[2];

    ExtractionConfigSet configs = MakeExtractionConfigs();
    bool wantsHelp = false;
    if (!ParseCommandLine(argc, argv, configs, wantsHelp)) return 1;
    if (wantsHelp) { PrintUsage(argv[0]); return 0; }

    std::cout << "\n========================================================================" << std::endl;
    std::cout << " Starting Signal Extraction for the Ring Observable" << std::endl;
    std::cout << "========================================================================\n" << std::endl;

    // Echo the configuration that is actually in force. Without this, a run with non-default
    // windows is indistinguishable from a default one in the logs, and the output file carries no
    // record of how it was made.
    {
        const SidebandConfig& c = configs.perBin;
        std::cout << "  Extraction windows:" << std::endl;
        std::cout << Form("    signal   : %s", SignalWindowLabel(c).Data()) << std::endl;
        std::cout << Form("    sidebands: %s", SidebandBandLabel(c).Data()) << std::endl;
        std::cout << Form("    bkg model: counts pol%d, <R>_bkg pol%d, needing >= %d points "
                          "and > %.1f raw counts",
                          c.bkgPolOrder, c.ringBkgPolOrder, c.minSidebandPoints, c.minSidebandCounts) << std::endl;
        std::cout << Form("    gates    : >= %.0f entries and > %.0f integral per spectrum",
                          c.minEntries, c.minIntegral) << std::endl;
        std::cout << Form("    peak fit : mu init %.5f +/- %.3f, sigma init %.5f in [%.4f, %.4f]",
                          c.muInitGuess, c.muLimitWindow, c.sigmaInitGuess,
                          c.sigmaLimitMin, c.sigmaLimitMax) << std::endl;
    }

    // =========================================================================================
    // Step 2: Parse Input Path and Set Up Output File
    // =========================================================================================
    std::cout << "\n[Step 2] Parsing input file path and setting up I/O..." << std::endl;
    
    // Extract suffix after "ConsumerResults_"
    std::string suffix = "";
    std::string searchStr = "ConsumerResults_";
    size_t pos = inputFilePath.find(searchStr);
    
    if (pos != std::string::npos){
        size_t startPos = pos + searchStr.length();
        size_t endPos = inputFilePath.find(".root", startPos);
        if (endPos != std::string::npos){
            suffix = inputFilePath.substr(startPos, endPos - startPos);
        }
        else {
            suffix = inputFilePath.substr(startPos); // Just in case .root is missing
        }
    }
    else {
        suffix = "UnknownSuffix";
        std::cout << "  Warning: 'ConsumerResults_' not found in input path. Using 'UnknownSuffix'." << std::endl;
    }

    // Ensure output directory exists (like mkdir -p)
    try {fs::create_directories(outputFolderPath);} 
    catch (const fs::filesystem_error& e) {
        std::cerr << "Error creating output directory: " << e.what() << std::endl;
        return 1;
    }

    std::string outFileName = outputFolderPath + "signalExtractionRing_" + suffix + ".root";
    std::cout << "  Input file: " << inputFilePath << std::endl;
    std::cout << "  Output file will be: " << outFileName << std::endl;

    // Open input file (Read-Only)
    TFile* inFile = TFile::Open(inputFilePath.c_str(), "READ");
    if (!inFile || inFile->IsZombie()){
        std::cerr << "  Error: Could not open input file!" << std::endl;
        return 1;
    }

    // Create output file (Recreate to avoid appending to old runs)
    TFile* outFile = TFile::Open(outFileName.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()){
        std::cerr << "  Error: Could not create output file!" << std::endl;
        if(inFile) inFile->Close();
        return 1;
    }

    // =========================================================================================
    // Step 2.5: V0 selection cut-flow signal extraction
    // =========================================================================================
    // The selection-flow QA histograms live in AnalysisResults_merged.root, which sits one level
    // ABOVE the folder holding the ConsumerResults files. For example:
    //   input:  .../ITSandTPC_min3ITS/results_consumer/ConsumerResults_BothHyperons.root
    //   QA:     .../ITSandTPC_min3ITS/AnalysisResults_merged.root
    // This file is the same for all consumer variations, so the extraction is done exactly once,
    // before the variation loop. It is cheap compared to the ring signal extraction itself.
    std::cout << "\n[Step 2.5] Running V0 selection cut-flow signal extraction..." << std::endl;

    fs::path consumerPath(inputFilePath);
    fs::path analysisResultsPath = consumerPath.parent_path().parent_path() / "AnalysisResults_merged.root";
    std::cout << "  Looking for: " << analysisResultsPath.string() << std::endl;

    TFile* qaFile = TFile::Open(analysisResultsPath.c_str(), "READ");
    if (!qaFile || qaFile->IsZombie()) {
        std::cerr << "  Warning: could not open " << analysisResultsPath.string()
                  << ". Skipping the cut-flow extraction and continuing.\n";
        if (qaFile) { qaFile->Close(); qaFile = nullptr; }
    }
    else {
        // The producer's histogram registry is stored under "lambdajetpolarizationions/".
        TDirectory* qaTaskDir = (TDirectory*)qaFile->Get("lambdajetpolarizationions/GeneralQA");
        if (qaTaskDir) std::cout << "  Found selection QA histograms in: lambdajetpolarizationions/GeneralQA" << std::endl;

        // ALTERNATIVE (kept for reference): if the task is ever renamed, the folder can be found
        // without hardcoding the O2Physics device name by looking for the top-level directory
        // that actually owns GeneralQA/h2dSelectionLambdaMass.
        //
        // TDirectory* qaTaskDir = nullptr;
        // TIter nextTopKey(qaFile->GetListOfKeys());
        // TKey* topKey;
        // while ((topKey = (TKey*)nextTopKey())) {
        //     TObject* topObj = topKey->ReadObj();
        //     TDirectory* candidateDir = dynamic_cast<TDirectory*>(topObj);
        //     if (!candidateDir) {
        //         delete topObj; // Not a directory: safe to drop (TDirectory lifetimes belong to the file)
        //         continue;
        //     }
        //     TDirectory* generalQADir = (TDirectory*)candidateDir->Get("GeneralQA");
        //     if (generalQADir && generalQADir->Get("h2dSelectionLambdaMass")) {
        //         qaTaskDir = generalQADir;
        //         std::cout << "  Found selection QA histograms in: " << candidateDir->GetName() << "/GeneralQA" << std::endl;
        //         break;
        //     }
        // }

        if (!qaTaskDir) {
            std::cerr << "  Warning: could not open lambdajetpolarizationions/GeneralQA in "
                      << analysisResultsPath.string() << ". Skipping the cut-flow extraction.\n";
        }
        else {
            TH2D* h2dSelLambdaMass = (TH2D*)qaTaskDir->Get("h2dSelectionLambdaMass");
            TH2D* h2dSelAntiLambdaMass = (TH2D*)qaTaskDir->Get("h2dSelectionAntiLambdaMass");
            TH1* hSelectionV0s = (TH1*)qaTaskDir->Get("hSelectionV0s");

            TDirectory* dirCutFlow = outFile->mkdir("V0SelectionCutFlow");

            // Nested chains for each species, as inclusive ROOT bin ranges in nesting order.
            // See the kNGenericV0Cuts / kNHypothesisCuts block at the top of this file for why the
            // two hypothesis branches must not be chained into one another.
            const std::vector<std::pair<int, int>> lambdaChain = {
                {1, kNGenericV0Cuts},                                                  // shared V0 cuts
                {kNGenericV0Cuts + 1, kNGenericV0Cuts + kNHypothesisCuts}               // Lambda block
            };
            const std::vector<std::pair<int, int>> antiLambdaChain = {
                {1, kNGenericV0Cuts},                                                  // shared V0 cuts
                {kNGenericV0Cuts + kNHypothesisCuts + 1, kNGenericV0Cuts + 2 * kNHypothesisCuts} // AntiLambda block
            };

            // Keep the raw selection flow next to the extracted yields, for traceability
            if (hSelectionV0s) {
                dirCutFlow->cd();
                hSelectionV0s->Write();
            }

            if (h2dSelLambdaMass) {
                dirCutFlow->cd();
                h2dSelLambdaMass->Write(); // Input histogram, saved for traceability
                PerformSelectionCutFlowExtraction(h2dSelLambdaMass, dirCutFlow, "Lambda",
                                                  "M_{p#pi^{-}} (GeV/#it{c}^{2})", lambdaChain,
                                                  configs.cutFlow);
            }
            else std::cerr << "  Warning: h2dSelectionLambdaMass missing. Skipping Lambda cut flow.\n";

            if (h2dSelAntiLambdaMass) {
                dirCutFlow->cd();
                h2dSelAntiLambdaMass->Write(); // Input histogram, saved for traceability
                PerformSelectionCutFlowExtraction(h2dSelAntiLambdaMass, dirCutFlow, "AntiLambda",
                                                  "M_{#bar{p}#pi^{+}} (GeV/#it{c}^{2})", antiLambdaChain,
                                                  configs.cutFlow);
            }
            else std::cerr << "  Warning: h2dSelectionAntiLambdaMass missing. Skipping AntiLambda cut flow.\n";
        }

        qaFile->Close();
        qaFile = nullptr;
    }

    outFile->cd(); // Restore a well-defined current directory before the variation loop starts

    // All results of the derived data consumer live inside a folder named after the task.
    const std::string baseFolderForAEE = "lambdajetpolarizationionsderived";

    // =========================================================================================
    // Step 2.75: Azimuthal Efficiency Effect probe -- signal extraction vs phi_Lambda - phi_p*
    // =========================================================================================
    // WHY THIS SPLIT IS WORTH EXTRACTING. <R> changes sign across phi_Lambda - phi_p*, so within a
    // single bin of it the background's dependence on mass is far closer to linear than it is in
    // the angle-integrated spectrum, where opposite-sign contributions pile up and leave a shape
    // the sideband polynomial has to absorb. That is precisely the regime a low-order sideband fit
    // is trustworthy in, which makes this the cleanest place to read <R>_S and <R>_B off.
    //
    // WHERE THEY LIVE. Unlike every other observable in this macro, these profiles are booked with
    // a bare path rather than (folder + "/..."), so they sit at the task root ALONGSIDE the
    // variation directories rather than once inside each of them. They are therefore extracted
    // once per file. If that placement was not deliberate, this is the line of code that shows it.
    //
    // NO COUNTS HISTOGRAM IS NEEDED. The consumer books no TH2D counterpart for these, and none is
    // required: CountsFromProfile2D reads the denominator out of the profile's own entry array,
    // which also guarantees the numerator and denominator describe the same candidates.
    {
        const std::string aeeFolder = baseFolderForAEE + "/HelicityEfficiencyQA/PhiLambdaPhiProtonStar";
        TDirectory* aeeIn = (TDirectory*)inFile->Get(aeeFolder.c_str());

        if (!aeeIn) {
            std::cout << "\n[Step 2.75] No HelicityEfficiencyQA/PhiLambdaPhiProtonStar folder in this "
                         "file; skipping the AEE probe extraction. This is expected for consumer "
                         "output produced before those profiles were added." << std::endl;
        } else {
            std::cout << "\n[Step 2.75] Extracting <R> vs phi_Lambda - phi_p* per mass bin..." << std::endl;
            TDirectory* aeeOut = outFile->mkdir("HelicityEfficiencyQA_PhiLambdaPhiProtonStar");

            // Three proxies x three species selections. The inclusive "LambdaLike" set is QA for the
            // competition between the two species: with equal Lambda and antiLambda yields the
            // dependence should cancel, so a residual there is itself the measurement of imbalance.
            struct AEESpec { const char* profileName; const char* outName; const char* axisTitle; };
            const std::vector<AEESpec> aeeSpecs = {
                {"p2dRingObservableLeadJetVsPhiLambdaLikePhiProtonStarVsMass", "LeadJet_LambdaLike",
                 "#phi_{#Lambda-like}-#phi_{p-like}^{*}"},
                {"p2dRingObservableLeadPVsPhiLambdaLikePhiProtonStarVsMass", "LeadP_LambdaLike",
                 "#phi_{#Lambda-like}-#phi_{p-like}^{*}"},
                {"p2dRingObservable2ndJetVsPhiLambdaLikePhiProtonStarVsMass", "SubJet_LambdaLike",
                 "#phi_{#Lambda-like}-#phi_{p-like}^{*}"},
                {"p2dRingObservableLeadJetVsPhiLambdaPhiProtonStarVsMass", "LeadJet_Lambda",
                 "#phi_{#Lambda}-#phi_{p}^{*}"},
                {"p2dRingObservableLeadPVsPhiLambdaPhiProtonStarVsMass", "LeadP_Lambda",
                 "#phi_{#Lambda}-#phi_{p}^{*}"},
                {"p2dRingObservable2ndJetVsPhiLambdaPhiProtonStarVsMass", "SubJet_Lambda",
                 "#phi_{#Lambda}-#phi_{p}^{*}"},
                {"p2dRingObservableLeadJetVsPhiAntiLambdaPhiProtonStarVsMass", "LeadJet_AntiLambda",
                 "#phi_{#bar{#Lambda}}-#phi_{#bar{p}}^{*}"},
                {"p2dRingObservableLeadPVsPhiAntiLambdaPhiProtonStarVsMass", "LeadP_AntiLambda",
                 "#phi_{#bar{#Lambda}}-#phi_{#bar{p}}^{*}"},
                {"p2dRingObservable2ndJetVsPhiAntiLambdaPhiProtonStarVsMass", "SubJet_AntiLambda",
                 "#phi_{#bar{#Lambda}}-#phi_{#bar{p}}^{*}"}
            };

            int nDone = 0, nSkipped = 0;
            for (const auto& spec : aeeSpecs) {
                TProfile2D* prof = (TProfile2D*)aeeIn->Get(spec.profileName);
                if (!prof) {
                    std::cout << "    -> Skipping " << spec.outName << ": " << spec.profileName
                              << " not found." << std::endl;
                    nSkipped++;
                    continue;
                }

                TH2D* counts = CountsFromProfile2D(prof, Form("h2dCountsAEE_%s", spec.outName));
                double aeeMassMin = prof->GetYaxis()->GetXmin();
                double aeeMassMax = prof->GetYaxis()->GetXmax();

                std::cout << "    -> Processing " << spec.outName << " ("
                          << prof->GetXaxis()->GetNbins() << " angular bins, "
                          << prof->GetEntries() << " entries)..." << std::endl;

                // CombinePerBin: no projection over phi_Lambda - phi_p*. See the note at the top
                // of the combination helpers -- projecting first would undo the entire reason this
                // folder splits on an angular variable across which <R> changes sign.
                ExtractObservable2D(counts, prof, aeeOut, spec.outName, spec.axisTitle,
                                    aeeMassMin, aeeMassMax,
                                    configs.perBin, configs.integrated,
                                    IntegralMode::CombinePerBin);
                delete counts;
                nDone++;
            }
            std::cout << Form("  [AEE probe] %d of %d extractions produced; %d profiles missing.",
                              nDone, (int)aeeSpecs.size(), nSkipped) << std::endl;
            outFile->cd();
        }
    }

    // =========================================================================================
    // Step 3: Loop Over Variations & Fetch Histograms
    // =========================================================================================
    std::cout << "\n[Step 3] Looping over variations and fetching histograms..." << std::endl;
    
    std::string baseFolder = baseFolderForAEE; // Same task folder the AEE block above used
    // Only "Ring" is generally booked by the consumer.
    // The other three families have their switches turned off in most cases.
    struct VariationSpec { std::string name; bool mandatory; };
    const std::vector<VariationSpec> variations = {
        {"Ring",                      true },
        {"RingKinematicCuts",         false},
        {"JetKinematicCuts",          false},
        {"JetAndLambdaKinematicCuts", false}
    };

    // --- Resolve which variations are actually present, once, before any work starts ---
    // Doing this up front buys two things: the absent ones are reported on a single line instead of
    // one warning per variation deep inside the loop, and no empty <Variation>/ directory is written
    // to the output file for a variation that was never produced. Previously mkdir() ran before the
    // existence check, so a trimmed run and a complete one looked identical in a TBrowser.
    std::vector<std::string> presentVariations;
    std::string presentList, absentList;
    for (const auto& v : variations) {
        if (inFile->Get((baseFolder + "/" + v.name).c_str())) {
            presentVariations.push_back(v.name);
            if (!presentList.empty()) presentList += ", ";
            presentList += v.name;
        } else if (v.mandatory) {
            std::cerr << "  Error: mandatory variation '" << v.name << "' not found in the input file.\n";
            outFile->Close();
            inFile->Close();
            return 1;
        } else {
            if (!absentList.empty()) absentList += ", ";
            absentList += v.name;
        }
    }

    bool sidebandPrecheckDone = false; // The mass axis is the same for every variation

    std::cout << "  Variations present : " << (presentList.empty() ? "(none)" : presentList) << std::endl;
    if (!absentList.empty())
        std::cout << "  Variations absent  : " << absentList << "  (optional, not enabled in this config)" << std::endl;

    for (const auto& var : presentVariations){
        std::cout << "\n#########################################" << std::endl;
        std::cout << "-> Processing variation: " << var << std::endl;
        std::cout << "#########################################\n" << std::endl;

        std::string fullFolderPath = baseFolder + "/" + var;
        TDirectory* inDir = (TDirectory*)inFile->Get(fullFolderPath.c_str());

        if (!inDir){ // Defensive only: the presence scan above should already have resolved this!
            std::cerr << "  Warning: Directory " << fullFolderPath << " not found in input file! Skipping." << std::endl;
            continue;
        }

        // Step 3.1: Output file gets one subfolder per PRESENT variation
        TDirectory* outDirVar = outFile->mkdir(var.c_str());
        if (!outDirVar){
            std::cerr << "  Warning: Could not create output directory '" << var << "'. Skipping." << std::endl;
            continue;
        }

        // Step 3.2: Fetching Histograms
        // 1D QA & Mass
        TH1D* hMassSigExtract = (TH1D*)inDir->Get("hMassSigExtract");
        
        // 2D: Observable vs Invariant Mass (Now using TProfiles to get the proper errors!)
        TProfile2D* p2dRingObservableDeltaPhiVsMass = (TProfile2D*)inDir->Get("p2dRingObservableDeltaPhiVsMass");
        TProfile2D* p2dRingObservableDeltaThetaVsMass = (TProfile2D*)inDir->Get("p2dRingObservableDeltaThetaVsMass");
        
            // Counters
        TH2D* h2dDeltaPhiVsMass = (TH2D*)inDir->Get("QA/h2dDeltaPhiVsMass");
        TH2D* h2dDeltaThetaVsMass = (TH2D*)inDir->Get("QA/h2dDeltaThetaVsMass");
        
        // 3D: Observable vs Mass vs Lambda pT (Now using TProfiles to get the proper errors!)
        TProfile3D* p3dRingObservableDeltaPhiVsMassVsLambdaPt = (TProfile3D*)inDir->Get("p3dRingObservableDeltaPhiVsMassVsLambdaPt");
        TProfile3D* p3dRingObservableDeltaThetaVsMassVsLambdaPt = (TProfile3D*)inDir->Get("p3dRingObservableDeltaThetaVsMassVsLambdaPt");
            // Counters
        TH3D* h3dDeltaPhiVsMassVsLambdaPt = (TH3D*)inDir->Get("QA/h3dDeltaPhiVsMassVsLambdaPt");
        TH3D* h3dDeltaThetaVsMassVsLambdaPt = (TH3D*)inDir->Get("QA/h3dDeltaThetaVsMassVsLambdaPt");

        // 3D: Observable vs Mass vs Lead Jet pT (Now using TProfiles to get the proper errors!)
        TProfile3D* p3dRingObservableDeltaPhiVsMassVsLeadJetPt = (TProfile3D*)inDir->Get("p3dRingObservableDeltaPhiVsMassVsLeadJetPt");
        TProfile3D* p3dRingObservableDeltaThetaVsMassVsLeadJetPt = (TProfile3D*)inDir->Get("p3dRingObservableDeltaThetaVsMassVsLeadJetPt"); 
            // Counters
        TH3D* h3dDeltaPhiVsMassVsLeadJetPt = (TH3D*)inDir->Get("QA/h3dDeltaPhiVsMassVsLeadJetPt");
        TH3D* h3dDeltaThetaVsMassVsLeadJetPt = (TH3D*)inDir->Get("QA/h3dDeltaThetaVsMassVsLeadJetPt");

        // 2D: Observable vs eta of the jet-side proxy vs Invariant Mass.
        // OPTIONAL: these were added to the consumer well after the angular ones, so any
        // ConsumerResults file produced before that simply does not contain them. Missing objects
        // are reported once and skipped -- they must not abort a run over older files, which is why
        // they are deliberately kept out of the essential-histogram check below.
        // Note the counters live at folder level rather than under QA/: they are inputs to the
        // extraction, not QA plots.
        TProfile2D* p2dRingObsEtaLeadJetVsMass = (TProfile2D*)inDir->Get("p2dRingObservableEtaLeadJetVsMass");
        TProfile2D* p2dRingObsEtaLeadPVsMass = (TProfile2D*)inDir->Get("p2dRingObservableLeadPEtaLeadPVsMass");
        TProfile2D* p2dRingObsEtaSubJetVsMass = (TProfile2D*)inDir->Get("p2dRingObservable2ndJetEta2ndJetVsMass");
            // Counters
        TH2D* h2dCounterEtaLeadJetVsMass = (TH2D*)inDir->Get("h2dCounterEtaLeadJetVsMass");
        TH2D* h2dCounterEtaLeadPVsMass = (TH2D*)inDir->Get("h2dCounterLeadPEtaLeadPVsMass");
        TH2D* h2dCounterEtaSubJetVsMass = (TH2D*)inDir->Get("h2dCounter2ndJetEta2ndJetVsMass");

        // Basic check to ensure critical histograms were loaded correctly
        if (!h2dDeltaPhiVsMass || !hMassSigExtract || !p2dRingObservableDeltaPhiVsMass || !p2dRingObservableDeltaThetaVsMass || !p3dRingObservableDeltaPhiVsMassVsLambdaPt ||
            !p3dRingObservableDeltaThetaVsMassVsLambdaPt || !p3dRingObservableDeltaPhiVsMassVsLeadJetPt || !p3dRingObservableDeltaThetaVsMassVsLeadJetPt) {
            std::cerr << "  Error: Missing essential histograms for " << var << ". Skipping.\n";
            continue;
        }

        // =========================================================================================
        // Steps 4, 5, and 6: Mass Projections, Fits, QA, and Signal Extraction (using Helper)
        // =========================================================================================
        std::cout << "[Steps 4-6] Projecting mass bins, fitting backgrounds, and extracting <R>..." << std::endl;

        // Get the mass limits from your standard reference histogram
        double massMin = hMassSigExtract->GetXaxis()->GetXmin();
        double massMax = hMassSigExtract->GetXaxis()->GetXmax();

        // Refuse to start if the requested sideband band cannot be fitted on this mass axis.
        // Done once per variation, against the reference spectrum, before any of the 39 extraction
        // calls below. Cheap, and it turns a silent 39-call no-op into one line.
        if (!sidebandPrecheckDone) {
            if (!CheckSidebandBandIsFittable(hMassSigExtract, configs.perBin)) {
                outFile->Close();
                inFile->Close();
                return 1;
            }
            sidebandPrecheckDone = true;
        }

        // -----------------------------------------------------------------------------------------
        // Extraction for Delta Phi Bins
        // -----------------------------------------------------------------------------------------
        std::cout << "    -> Processing Delta Phi observables..." << std::endl;
        ExtractObservable2D(h2dDeltaPhiVsMass, 
                            p2dRingObservableDeltaPhiVsMass, 
                            outDirVar,       // The parent TDirectory to save everything inside
                            "DeltaPhi",      // This creates a "DeltaPhi" subfolder for organized output
                            "#Delta#phi",    // Axis title for generated histograms
                            massMin, 
                            massMax,
                            configs.perBin, configs.integrated);
        // std::cout<< "DEBUG! Integral of p2dRingObservableDeltaPhiVsMass: " << p2dRingObservableDeltaPhiVsMass->Integral() << std::endl;
        // std::cout<< "DEBUG! Integral of h2dDeltaPhiVsMass: " << h2dDeltaPhiVsMass->Integral() << std::endl;
        // -----------------------------------------------------------------------------------------
        // Extraction for Delta Theta Bins
        // -----------------------------------------------------------------------------------------
        std::cout << "    -> Processing Delta Theta observables..." << std::endl;
        ExtractObservable2D(h2dDeltaThetaVsMass, 
                            p2dRingObservableDeltaThetaVsMass, 
                            outDirVar,       // The parent TDirectory to save everything inside
                            "DeltaTheta",    // This creates a "DeltaTheta" subfolder for organized output
                            "#Delta(#theta)",   // Axis title for generated histograms
                            massMin, 
                            massMax,
                            configs.perBin, configs.integrated);
        // std::cout<< "DEBUG! Integral of p2dRingObservableDeltaThetaVsMass: " << p2dRingObservableDeltaThetaVsMass->Integral() << std::endl;
        // std::cout<< "DEBUG! Integral of h2dDeltaThetaVsMass: " << h2dDeltaThetaVsMass->Integral() << std::endl;

        // -----------------------------------------------------------------------------------------
        // Integrated extraction, once per jet-side proxy
        // -----------------------------------------------------------------------------------------
        // The 1D <R> vs mass profiles are written by every consumer that has ever run, so this
        // needs no re-run to produce. QA/hMass is a 450-bin spectrum over the same mass region and
        // is used, when present, to pin down (mu, sigma) far better than the 24-bin extraction axis
        // can. Nothing is COUNTED across the two axes -- their edges do not align and they do not
        // even cover the same range -- only the peak position and width are borrowed.
        TH1D* hMassFine = (TH1D*)inDir->Get("QA/hMass");
        std::cout << "\n[Integrated] Extracting <R>_measured, <R>_S and <R>_B per proxy"
                  << (hMassFine ? ", with (mu, sigma) from the fine QA mass axis"
                                : " (no fine QA mass spectrum found; fitting the extraction axis)")
                  << "..." << std::endl;

        const std::vector<IntegratedProxySpec> proxySpecs = {
            {"pRingObservableMass",        "LeadJet", "Leading jet"},
            {"pRingObservableLeadPMass",   "LeadP",   "Leading particle"},
            {"pRingObservable2ndJetMass",  "SubJet",  "Subleading jet"}
        };
        for (const auto& ps : proxySpecs) {
            TProfile* profMass = (TProfile*)inDir->Get(ps.profilePath);
            ExtractIntegratedFromProfile(profMass, hMassFine, outDirVar, ps, configs.integrated);
        }

        // -----------------------------------------------------------------------------------------
        // Extraction vs the eta of the jet-side proxy
        // -----------------------------------------------------------------------------------------
        // Same engine, different pair of input objects: ExtractObservable2D only ever needed a
        // (counts TH2D, TProfile2D) pair with mass on the Y axis, so nothing about the eta
        // dependence needs new machinery. The eta binning is read from the histogram itself
        // (nBins = h2dCounts->GetNbinsX(), and every output is cloned from ProjectionX), so
        // changing axisEtaSigExtract in the consumer requires no edit here.
        //
        // Unlike DeltaPhi and DeltaTheta -- which are two angular views of the SAME candidates and
        // therefore give identical angle-integrated numbers -- these three genuinely differ, because
        // a Lambda with no subleading jet never enters the SubJet profile at all. Their integrated
        // values are three separate measurements, not three views of one.
        struct EtaExtractionSpec {
            TH2D*        counts;
            TProfile2D*  profile;
            const char*  extractionName;
            const char*  axisTitle;
            const char*  countsPath;  // Only for the "missing input" message
            const char*  profilePath;
        };
        const std::vector<EtaExtractionSpec> etaExtractions = {
            {h2dCounterEtaLeadJetVsMass, p2dRingObsEtaLeadJetVsMass, "EtaLeadJet", "#eta_{Jet}",
             "h2dCounterEtaLeadJetVsMass", "p2dRingObservableEtaLeadJetVsMass"},
            {h2dCounterEtaLeadPVsMass, p2dRingObsEtaLeadPVsMass, "EtaLeadP", "#eta_{LeadP}",
             "h2dCounterLeadPEtaLeadPVsMass", "p2dRingObservableLeadPEtaLeadPVsMass"},
            {h2dCounterEtaSubJetVsMass, p2dRingObsEtaSubJetVsMass, "EtaSubJet", "#eta_{SubJet}",
             "h2dCounter2ndJetEta2ndJetVsMass", "p2dRingObservable2ndJetEta2ndJetVsMass"}
        };

        for (const auto& spec : etaExtractions) {
            if (!spec.counts || !spec.profile) {
                std::cout << "    -> Skipping " << spec.extractionName << ": "
                          << (spec.counts ? spec.profilePath : spec.countsPath)
                          << " not found in this file. This is expected for ConsumerResults files "
                             "produced before the eta-vs-mass histograms were added." << std::endl;
                continue;
            }
            std::cout << "    -> Processing " << spec.extractionName << " observables..." << std::endl;
            ExtractObservable2D(spec.counts,
                                spec.profile,
                                outDirVar,
                                spec.extractionName,
                                spec.axisTitle,
                                massMin,
                                massMax,
                                configs.perBin, configs.integrated);
        }
        // =========================================================================================
        // Step 7: 3D Projections and Slicing (Moving Kinematic Windows)
        // =========================================================================================
        std::cout << "\n[Step 7] Slicing 3D histograms and extracting kinematic windows..." << std::endl;

        // Define the kinematic windows for Lambda pT and Leading Jet pT
        std::vector<std::pair<double, double>> lambdaPtWindows = {{0, 0.5}, {0.5, 1.5}, {1.5, 3.0}, {3.0, 6.0}, {6.0, 8.0}, {8.0, 15.0}, {15.0, 30.0}, {30.0, 50.0}};
        std::vector<std::pair<double, double>> leadingJetPtWindows = {{0, 5}, {5, 10}, {10, 20}, {20, 30}, {30, 40}, {40, 60}, {60, 80}, {80, 100}, {100, 200}};

        // Helper lambda to format the Pt strings (e.g., "0.5" -> "05", "3.0" -> "3")
        auto fmtWindow = [](double val) {
            TString s = Form("%.1f", val);
            s.ReplaceAll(".0", ""); // 3.0 -> 3
            s.ReplaceAll(".", "");  // 0.5 -> 05, 1.5 -> 15
            return s;
        };

        // Before entering the loops, we actually build TH3Ds based on the TProfile3Ds that have been imported from the derived data consumer output.
        // This ensures proper error propagation of <R> (just use error of the mean formula!), and keeps the same TH3D structure from before!
            // When you call p3dRingObs->Project3D("yx") on a TProfile3D, ROOT produces a TProfile2D where each (ix, iy) cell contains the grand mean
            // of all R values that fell into any iz bin within the Z range. That is NOT what we want here, so we need to do the dirty work of converting
            // the TProfile3Ds into TH3Ds before doing any projections. Error propagation would be wrong otherwise!
        TH3D* h3dNumCorrErrDeltaPhiVsMassVsLambdaPt = BuildNumFromProfile3D(p3dRingObservableDeltaPhiVsMassVsLambdaPt, h3dDeltaPhiVsMassVsLambdaPt, Form("h3dNumCorrErr_DeltaPhi_LambdaPt_%s", var.c_str()));
        TH3D* h3dNumCorrErrDeltaThetaVsMassVsLambdaPt = BuildNumFromProfile3D(p3dRingObservableDeltaThetaVsMassVsLambdaPt, h3dDeltaThetaVsMassVsLambdaPt, Form("h3dNumCorrErr_DeltaTheta_LambdaPt_%s", var.c_str()));
        TH3D* h3dNumCorrErrDeltaPhiVsMassVsLeadJetPt = BuildNumFromProfile3D(p3dRingObservableDeltaPhiVsMassVsLeadJetPt, h3dDeltaPhiVsMassVsLeadJetPt, Form("h3dNumCorrErr_DeltaPhi_LeadJetPt_%s", var.c_str()));
        TH3D* h3dNumCorrErrDeltaThetaVsMassVsLeadJetPt = BuildNumFromProfile3D( p3dRingObservableDeltaThetaVsMassVsLeadJetPt, h3dDeltaThetaVsMassVsLeadJetPt, Form("h3dNumCorrErr_DeltaTheta_LeadJetPt_%s", var.c_str()));
        // Null check for corrected histograms (BuildNumFromProfile3D returns nullptr on failure)
        if (!h3dNumCorrErrDeltaPhiVsMassVsLambdaPt || !h3dNumCorrErrDeltaThetaVsMassVsLambdaPt || !h3dNumCorrErrDeltaPhiVsMassVsLeadJetPt || !h3dNumCorrErrDeltaThetaVsMassVsLeadJetPt) {
            std::cerr << "  Error: BuildNumFromProfile3D failed for " << var << ". Skipping.\n";
            continue;
        }

        // -----------------------------------------------------------------------------------------
        // 7.A: Lambda pT Projections
        // -----------------------------------------------------------------------------------------
        std::cout << "    -> Processing Lambda pT 3D windows..." << std::endl;
        TDirectory* dir3D_LambdaPt = outDirVar->mkdir("3D_Projections_LambdaPt");

        for (const auto& window : lambdaPtWindows) {
            double ptMin = window.first;
            double ptMax = window.second;
            TString ptStr = Form("Pt_%s_%s", fmtWindow(ptMin).Data(), fmtWindow(ptMax).Data()); // Formats pt values such as 0.5 to 05 or 3.0 to 3 in the folder name!
            
            // Find bins on the Z-axis (using a tiny epsilon to ensure we are safely inside the bin)
            int zBinMin = h3dDeltaPhiVsMassVsLambdaPt->GetZaxis()->FindBin(ptMin + 1e-4);
            int zBinMax = h3dDeltaPhiVsMassVsLambdaPt->GetZaxis()->FindBin(ptMax - 1e-4);

            // --- Delta Phi ---
            // Set the Z-axis ranges for the 3D histograms (Lambda Pt ranges)
            h3dDeltaPhiVsMassVsLambdaPt->GetZaxis()->SetRange(zBinMin, zBinMax);
            h3dNumCorrErrDeltaPhiVsMassVsLambdaPt->GetZaxis()->SetRange(zBinMin, zBinMax);

            // Project to 2D using "yx e" (Y vs X, which places Mass on the Y-axis and Angle on the X-axis)
            // ROOT's Project3D allows naming the output histogram by putting the name followed by a space before the option
                // This projection now carries correct errors because h3dNumCorr* also does!
            TH2D* h2dCountsPhi_Lpt = (TH2D*)h3dDeltaPhiVsMassVsLambdaPt->Project3D(Form("h2dCountsPhi_Lpt_%s yx e", ptStr.Data()));
            TH2D* h2dCorrErrNumPhi_Lpt = (TH2D*)h3dNumCorrErrDeltaPhiVsMassVsLambdaPt->Project3D(Form("h2dNumPhi_Lpt_%s yx e", ptStr.Data()));
            // These projections are pure intermediates: nothing downstream ever fetches them back
            // from a directory, they are only passed along by pointer. Detaching them makes this
            // scope their sole owner, so they can be deleted deterministically below.
            if (h2dCountsPhi_Lpt) h2dCountsPhi_Lpt->SetDirectory(nullptr);
            if (h2dCorrErrNumPhi_Lpt) h2dCorrErrNumPhi_Lpt->SetDirectory(nullptr);

            // Convert the TH2D with corrected error bars into a TProfile2D so ExtractObservable2D can keep its TProfile2D signature:
            TProfile2D* p2dPhi_Lpt = ConvertToProfile2D(h2dCorrErrNumPhi_Lpt, h2dCountsPhi_Lpt, Form("p2dPhi_Lpt_%s", ptStr.Data()));

            // Send to helper!
            ExtractObservable2D(h2dCountsPhi_Lpt, p2dPhi_Lpt, dir3D_LambdaPt, Form("DeltaPhi_%s", ptStr.Data()), "#Delta#phi", massMin, massMax, configs.perBin, configs.integrated,
                            IntegralMode::ProjectThenExtract, false);
            delete p2dPhi_Lpt;          // We own all three now
            delete h2dCountsPhi_Lpt;
            delete h2dCorrErrNumPhi_Lpt;

            // --- Cos Theta ---
            // Set the Z-axis ranges (Lambda Pt ranges)
            h3dDeltaThetaVsMassVsLambdaPt->GetZaxis()->SetRange(zBinMin, zBinMax);
            h3dNumCorrErrDeltaThetaVsMassVsLambdaPt->GetZaxis()->SetRange(zBinMin, zBinMax);

            TH2D* h2dCountsTheta_Lpt = (TH2D*)h3dDeltaThetaVsMassVsLambdaPt->Project3D(Form("h2dCountsTheta_Lpt_%s yx e", ptStr.Data()));
            TH2D* h2dNumCorrErrTheta_Lpt    = (TH2D*)h3dNumCorrErrDeltaThetaVsMassVsLambdaPt->Project3D(Form("h2dNumTheta_Lpt_%s yx e", ptStr.Data()));
            if (h2dCountsTheta_Lpt) h2dCountsTheta_Lpt->SetDirectory(nullptr);
            if (h2dNumCorrErrTheta_Lpt) h2dNumCorrErrTheta_Lpt->SetDirectory(nullptr);

            // Convert the TH2D with corrected error bars into a TProfile2D so ExtractObservable2D can keep its TProfile2D signature:
            TProfile2D* p2dTheta_Lpt = ConvertToProfile2D(h2dNumCorrErrTheta_Lpt, h2dCountsTheta_Lpt, Form("p2dTheta_Lpt_%s", ptStr.Data()));

            // Send to helper!
            ExtractObservable2D(h2dCountsTheta_Lpt, p2dTheta_Lpt, dir3D_LambdaPt, Form("DeltaTheta_%s", ptStr.Data()), "#Delta#theta", massMin, massMax, configs.perBin, configs.integrated,
                            IntegralMode::ProjectThenExtract, false);
            delete p2dTheta_Lpt;
            delete h2dCountsTheta_Lpt;
            delete h2dNumCorrErrTheta_Lpt;
        }

        // -----------------------------------------------------------------------------------------
        // 7.B: Leading Jet pT Projections
        // -----------------------------------------------------------------------------------------
        std::cout << "    -> Processing Leading Jet pT 3D windows..." << std::endl;
        TDirectory* dir3D_LeadJetPt = outDirVar->mkdir("3D_Projections_LeadJetPt");

        for (const auto& window : leadingJetPtWindows) {
            double ptMin = window.first;
            double ptMax = window.second;
            TString ptStr = Form("Pt_%s_%s", fmtWindow(ptMin).Data(), fmtWindow(ptMax).Data());
            
            int zBinMin = h3dDeltaPhiVsMassVsLeadJetPt->GetZaxis()->FindBin(ptMin + 1e-4);
            int zBinMax = h3dDeltaPhiVsMassVsLeadJetPt->GetZaxis()->FindBin(ptMax - 1e-4);

            // --- Delta Phi ---
            h3dDeltaPhiVsMassVsLeadJetPt->GetZaxis()->SetRange(zBinMin, zBinMax); // Leading jet Pt projection
            h3dNumCorrErrDeltaPhiVsMassVsLeadJetPt->GetZaxis()->SetRange(zBinMin, zBinMax);

            TH2D* h2dCountsPhi_Jpt = (TH2D*)h3dDeltaPhiVsMassVsLeadJetPt->Project3D(Form("h2dCountsPhi_Jpt_%s yx e", ptStr.Data()));
            TH2D* h2dNumCorrErrPhi_Jpt = (TH2D*)h3dNumCorrErrDeltaPhiVsMassVsLeadJetPt->Project3D(Form("h2dNumPhi_Jpt_%s yx e", ptStr.Data()));
            if (h2dCountsPhi_Jpt) h2dCountsPhi_Jpt->SetDirectory(nullptr);
            if (h2dNumCorrErrPhi_Jpt) h2dNumCorrErrPhi_Jpt->SetDirectory(nullptr);

            // Convert the TH2D with corrected error bars into a TProfile2D so ExtractObservable2D can keep its TProfile2D signature:
            TProfile2D* p2dPhi_Jpt = ConvertToProfile2D(h2dNumCorrErrPhi_Jpt, h2dCountsPhi_Jpt, Form("p2dPhi_Jpt_%s", ptStr.Data()));
            ExtractObservable2D(h2dCountsPhi_Jpt, p2dPhi_Jpt, dir3D_LeadJetPt, Form("DeltaPhi_%s", ptStr.Data()), "#Delta#phi", massMin, massMax, configs.perBin, configs.integrated,
                            IntegralMode::ProjectThenExtract, false);
            delete p2dPhi_Jpt;
            delete h2dCountsPhi_Jpt;
            delete h2dNumCorrErrPhi_Jpt;

            // --- Cos Theta ---
            h3dDeltaThetaVsMassVsLeadJetPt->GetZaxis()->SetRange(zBinMin, zBinMax);
            h3dNumCorrErrDeltaThetaVsMassVsLeadJetPt->GetZaxis()->SetRange(zBinMin, zBinMax);

            TH2D* h2dCountsTheta_Jpt = (TH2D*)h3dDeltaThetaVsMassVsLeadJetPt->Project3D(Form("h2dCountsTheta_Jpt_%s yx e", ptStr.Data()));
            TH2D* h2dNumCorrErrTheta_Jpt = (TH2D*)h3dNumCorrErrDeltaThetaVsMassVsLeadJetPt->Project3D(Form("h2dNumTheta_Jpt_%s yx e", ptStr.Data()));
            if (h2dCountsTheta_Jpt) h2dCountsTheta_Jpt->SetDirectory(nullptr);
            if (h2dNumCorrErrTheta_Jpt) h2dNumCorrErrTheta_Jpt->SetDirectory(nullptr);

            // Convert the TH2D with corrected error bars into a TProfile2D so ExtractObservable2D can keep its TProfile2D signature:
            TProfile2D* p2dTheta_Jpt = ConvertToProfile2D(h2dNumCorrErrTheta_Jpt, h2dCountsTheta_Jpt, Form("p2dTheta_Jpt_%s", ptStr.Data()));
            ExtractObservable2D(h2dCountsTheta_Jpt, p2dTheta_Jpt, dir3D_LeadJetPt, Form("DeltaTheta_%s", ptStr.Data()), "#Delta#theta", massMin, massMax, configs.perBin, configs.integrated,
                            IntegralMode::ProjectThenExtract, false);
            delete p2dTheta_Jpt;
            delete h2dCountsTheta_Jpt;
            delete h2dNumCorrErrTheta_Jpt;
        }

        // Cleanup -- Delete the corrected error bar TH3D objects we just created!
            // The projected TH2D pointers (h2dCountsPhi_Lpt etc.) are now detached from any
            // directory and deleted inside the loops above, right after each extraction.
        delete h3dNumCorrErrDeltaPhiVsMassVsLambdaPt;
        delete h3dNumCorrErrDeltaThetaVsMassVsLambdaPt;
        delete h3dNumCorrErrDeltaPhiVsMassVsLeadJetPt;
        delete h3dNumCorrErrDeltaThetaVsMassVsLeadJetPt;

        // =========================================================================================
        // Step 8: Save Input Histograms for Traceability
        // =========================================================================================
        std::cout << "  [Step 8] Saving all input histograms for " << var << "..." << std::endl;

        TDirectory* dirInput = outDirVar->mkdir("Input_Histograms");
        
        // Iterate over all keys in the current variation's input directory
        TIter nextKey(inDir->GetListOfKeys());
        TKey *key;
        while ((key = (TKey*)nextKey())) {
            TObject *obj = key->ReadObj();
            if (!obj) continue;
            if (obj->InheritsFrom(TH1::Class())) {
                dirInput->cd();
                obj->Write();
                delete obj;
            }
            else if (obj->InheritsFrom(TDirectory::Class())) {
                // Do NOT recurse or delete TDirectory objects here:
                // their lifetime is managed by the file. Simply skip.
                // (If you want to copy the QA subfolder too, handle it separately.)
                continue;
            }
            else {
                delete obj; // Safe for TTree, TNamed, etc.
            }
        }

        inFile->cd(fullFolderPath.c_str()); // Restore gDirectory to the input variation folder
                                            // so Step 9's inDir->Get() calls associate objects
                                            // with the input file, not the output file.

        // =========================================================================================
        // Step 9: Calculating significance plots for quick 1D QA
        // =========================================================================================
        std::cout << "  [Step 9] Calculating significance plots for quick 1D QA, in variation: " << var << "..." << std::endl;
        // Significance plots based on all 1D TProfiles (NSigma with which each bin deviates from zero):
        // (This has to be done AFTER processing the derived data, because these cannot be merged when code is pipelined!)
        TDirectory* outDirSig = outDirVar->mkdir("1D_Significance");

        // Fetching 1D TProfiles:
        TProfile* pDeltaPhi = (TProfile*)inDir->Get("pRingObservableDeltaPhi");
        TProfile* pDeltaTheta = (TProfile*)inDir->Get("pRingObservableDeltaTheta");
        TProfile* pIntegrated = (TProfile*)inDir->Get("pRingObservableIntegrated");
        TProfile* pLambdaPt = (TProfile*)inDir->Get("pRingObservableLambdaPt");
        TProfile* pMass = (TProfile*)inDir->Get("pRingObservableMass"); 

        outDirSig->cd();

        TGraphErrors* gSigDeltaPhi   = makeSignificanceGraph(pDeltaPhi,   "gRingSignificanceDeltaPhi");
        TGraphErrors* gSigDeltaTheta = makeSignificanceGraph(pDeltaTheta, "gRingSignificanceDeltaTheta");
        TGraphErrors* gSigIntegrated = makeSignificanceGraph(pIntegrated, "gRingSignificanceIntegrated");
        TGraphErrors* gSigLambdaPt   = makeSignificanceGraph(pLambdaPt,   "gRingSignificanceLambdaPt");
        TGraphErrors* gSigMass       = makeSignificanceGraph(pMass,       "gRingSignificanceMass");

        gSigDeltaPhi->Write();
        gSigDeltaTheta->Write();
        gSigIntegrated->Write();
        gSigLambdaPt->Write();
        gSigMass->Write();

        // Draw and save canvases
        drawSigAndDashedLines(gSigDeltaPhi,   outDirSig, "cSigDeltaPhi_" + var);
        drawSigAndDashedLines(gSigDeltaTheta, outDirSig, "cSigDeltaTheta_" + var);
        drawSigAndDashedLines(gSigIntegrated, outDirSig, "cSigIntegrated_" + var);
        drawSigAndDashedLines(gSigLambdaPt,   outDirSig, "cSigLambdaPt_" + var);
        drawSigAndDashedLines(gSigMass,       outDirSig, "cSigMass_" + var);

        // =========================================================================================
        // Step 10 (QA): Denominator-only signal extraction, integrated over all angles
        // =========================================================================================
        // Fit the invariant mass spectrum of ALL selected Lambda candidates with a gaus+pol2,
        // define the signal and sideband regions, and save a summary canvas.
        // hMassSigExtract is the dedicated denominator histogram from the input file -- it
        // already contains the full angle-integrated counts, so no projection is needed.
        std::cout << "\n[Step 10 DenomQA] Running denominator QA for " << var << "..." << std::endl;
        // Create a dedicated subdirectory to keep the QA output self-contained
        TDirectory* dirDenomQA = outDirVar->mkdir("DenominatorQA");

        // Call function:
        PerformDenominatorQA(hMassSigExtract, dirDenomQA,
                                TString(var),  // unique label for ROOT name generation
                                massMin, massMax, configs.denomQA);
    } // <--- This closes the variation loop ("Ring", "RingKinematicCuts", etc.)

    // Clean up files
    std::cout << "\n#############################################" << std::endl;
    std::cout << "Writing everything to disk and closing output files!" << std::endl;
    std::cout << "#############################################\n" << std::endl;
    outFile->Close();
    inFile->Close();
    
    std::cout << "\n========================================================================" << std::endl;
    std::cout << " Signal Extraction Complete. File saved as: " << outFileName << std::endl;
    std::cout << "========================================================================\n" << std::endl;
    
    return 0; // Exit successfully
}