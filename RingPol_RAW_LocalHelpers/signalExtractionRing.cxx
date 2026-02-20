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
 **************************************************************************************************/
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

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

// Extra includes for simultaneous signal extraction fit to Ring Observable numerator and Counts denominator:
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"


constexpr double lambdaPDGMassApprox = 1.11568; // Just for some fit initial guesses

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
    hNum->Reset();                     // Clear all content and errors
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
SimFitResult PerformSimultaneousFitQA(TH1D* hMassDensity, TH1D* hNumDensity, double massMin, double massMax, double init_sigYield, double init_mu, double init_sigma, TF1* bkgFitFunc, double init_RS, double init_RB) {
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
    // MAJOR CHANGE: all seeds are now divided by normFactor (the density peak),
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

        // Full Gaussian integral = amplitude * sigma * sqrt(2*pi)
        // The signal yield in the +/- 4 sigma window is ~99.99366% of the total integral.
        double raw_SigYield     = amplitude_abs * sigmaVal * std::sqrt(2.0 * TMath::Pi());
        double raw_errSigYield  = err_amplitude_abs * sigmaVal * std::sqrt(2.0 * TMath::Pi());
        result.SigYield    = raw_SigYield    * 0.9999366;
        result.err_SigYield = raw_errSigYield * 0.9999366;

        // =====================================================================
        // Recover absolute background yield by analytically integrating the
        // fitted background density over the +/- 4 sigma signal window.
        // MAJOR CHANGE: the integration is now entirely in density space.
        // integral of (c0 + c1*m + c2*m^2) dm from xLow to xHigh gives counts
        // directly (density [counts/GeV] * mass interval [GeV] = counts).
        // We integrate the ABSOLUTE density (par_norm * normFactor), so the
        // final result is already in counts with no further rescaling needed
        // beyond the normFactor multiplication applied to c0_norm, c1_norm, c2_norm.
        // =====================================================================
        double xLow  = muVal - 4.0 * sigmaVal; // Wrote as 4.0*sigmaVal to be easy to change
        double xHigh = muVal + 4.0 * sigmaVal; // (TODO: update if signal region definition changes)

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
        // Note: par[0] is SigAmplitude_Norm; SigYield = par[0]*normFactor*sigmaVal*sqrt(2pi)*0.9999366,
        // so Var(SigYield) carries an extra factor of (sigmaVal*sqrt(2pi)*0.9999366)^2 relative
        // to Var(par[0]). For Cov(S,B) the same factor applies to the par[0] row/column.
        double gaussFactor = sigmaVal * std::sqrt(2.0 * TMath::Pi()) * 0.9999366;

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
void ExtractObservable2D(TH2D* h2dCounts, TProfile2D* p2dRingObs, TDirectory* parentDir, TString extractionName, TString axisTitle, double massMin, double massMax, bool printHeader = true){
    if (printHeader) std::cout << "\n[ExtractObservable2D] Starting extraction: " << extractionName << std::endl;
    // Create subdirectories for organized output
    TDirectory* dirBase = parentDir->mkdir(extractionName);
    TDirectory* dirFits = dirBase->mkdir("MassFits");
    TDirectory* dirQARingMass = dirBase->mkdir("QA_RingObservable_vs_Mass");
    // TDirectory* dirBkgFits = dirBase->mkdir("BackgroundFits");
    TDirectory* dirBkgNumFits = dirBase->mkdir("BkgNumFits");
    TDirectory* dirResults = dirBase->mkdir("Results");
    TDirectory* dirResultsSim = dirBase->mkdir("ResultsCombinedFit");

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

    // We will need to store the extracted mu and sigma for Step 5 & 6
    // Using a simple struct to hold the fit results for each bin
    struct FitResult { double mu; double sigma; bool valid; };
    std::vector<FitResult> fitResults(nBins + 1);

    // =========================================================================================
    // Step 4: Mass Projections, Fits, and QA Plots for each Bin
    // =========================================================================================
    if (printHeader) std::cout << "  -> Step 4: Mass fits and QA..." << std::endl;
    for (int iBin = 1; iBin <= nBins; ++iBin){
        // 1. Project the mass spectrum (Counts/Denominator) using "e" for error propagation
        TString projName = Form("hMass_Bin%d", iBin);
        TH1D* hMassProj = h2dCounts->ProjectionY(projName, iBin, iBin, "e");

        TAxis* xAxis = h2dCounts->GetXaxis();
        double xLow  = xAxis->GetBinLowEdge(iBin);
        double xHigh = xAxis->GetBinUpEdge(iBin);
        hMassProj->SetTitle(Form("Angle bin %d [%.4f,%.4f)];M_{p#pi} (GeV/c^{2});Counts", iBin, xLow, xHigh));

        // --- STABILITY CHECK: Does this bin even have data? ---
        // If there are fewer than 30 counts in the entire mass window, a good fit is probably
        // impossible (and if possible, useless).
        if (hMassProj->GetEntries() < 30 || hMassProj->Integral() <= 30) {
            fitResults[iBin] = {0, 0, false};
            // We don't even want to see these problematic histograms, so just skip right away!
            continue;
        }

        // 2. Project the Numerator to calculate <R>(m) QA plot
        TString numProjName = Form("hNum_Bin%d", iBin);
        TH1D* hNumProj = h2dNum->ProjectionY(numProjName, iBin, iBin, "e");
        
        TString qaName = Form("hQARingVsMass_Bin%d", iBin);
        TH1D* hQARing = (TH1D*)hNumProj->Clone(qaName);
        hQARing->SetTitle(Form("<R> vs Mass for Bin %d;M_{p#pi} (GeV/c^{2});<R>", iBin));
        hQARing->Divide(hNumProj, hMassProj, 1.0, 1.0, "B"); // Binomial errors might be tricky here, standard divide is okay for now as it's weighted
        
        dirQARingMass->cd();
        hQARing->Write();

        // 3. Setup the Fit Function (gaus + pol2)
        TString fitName = Form("fit_Bin%d", iBin);
        TF1* fitFunc = new TF1(fitName, "gaus(0) + pol2(3)", massMin, massMax); // Short-hand initialization of fit function

            // Scale the histogram to densities (easier to QA the fit):

        // Project the corrected hNum -- this now carries proper sigma(Sum_R_i) errors:
        // Now using a clean, retrievable name (title stays descriptive separately).
        TH1D* hMassProjDensity = (TH1D*)hMassProj->Clone(Form("hMassDensity_Bin%d", iBin));
        hMassProjDensity->SetTitle(Form("Angle bin %d [%.4f,%.4f);M_{p#pi} (GeV/c^{2});dN/dM", iBin, xLow, xHigh));
        hMassProjDensity->Scale(1.0, "width");
        
        // Smart initialization (starting parameters should be somewhat physical)
        // double maxCount = hMassProj->GetMaximum();
        double maxDensity = hMassProjDensity->GetMaximum();
        fitFunc->SetParameter(0, maxDensity);        // Constant
        fitFunc->SetParameter(1, lambdaPDGMassApprox);         // Mean close to PDG value
        fitFunc->SetParameter(2, 0.002);           // Sigma (typical for Lambda, derived from some earlier analysis of mine)

        // 4. Apply Stability Constraints
        fitFunc->SetParameter(0, maxDensity);        // Constant (Height)
        // - Constrain mean within mass axis limits
        fitFunc->SetParLimits(1, lambdaPDGMassApprox-0.01, lambdaPDGMassApprox+0.01); // Constraints close to mass (1.15 to 1.125)
        // Sigma: Restrict to physical limits (e.g., 0.1 to 5 MeV)
        fitFunc->SetParLimits(2, 0.0001, 0.005);

        // 5. Fit using Binned Poisson Log-Likelihood ("L"), Quiet ("Q"), No Draw ("0"), Respect Range ("R")
        // We use "S" to get the TFitResultPtr to check if the fit converged
        // TFitResultPtr r = hMassProj->Fit(fitFunc, "Q 0 R S"); // Removed log-likelihood to make fits behave better in high-statistics 1D environment!
            // To fit with variable bin widths, you must use the "I" option, which does not compare f(x_center), but to \int_{bin start}^{bin end} f(x)dx
        // TFitResultPtr r = hMassProj->Fit(fitFunc, "Q 0 R S I");

        // It is actually just way more natural to just fit the density histogram! Easier to perform QA later on too!
        TFitResultPtr r = hMassProjDensity->Fit(fitFunc, "Q 0 R S");

        dirFits->cd();
        hMassProjDensity->GetListOfFunctions()->Add(fitFunc); // Attach fit to histogram for saving later
            // Save the histogram only if the projection is not empty!
        if (hMassProj->GetEntries() > 0) hMassProjDensity->Write(); // For QA

        // Store results if fit is somewhat reasonable
        int fitStatus = int(r); // Extract the integer status code from the fit result
        
        // 0 = perfectly converged. 
        if (fitStatus == 0) {
            fitResults[iBin] = {fitFunc->GetParameter(1), fitFunc->GetParameter(2), true};
        } else {
            // Most of the other errors in this stage are simply error code 1, which means Minuit had problems with
            // the covariance matrix (the parameter "error" matrix). As we don't care about errors right now, this
            // is perfectly fine and we can just proceed. Thus, added a "salvage" block that just checks if the
            // fit parameters are somewhat OK before proceeding.

            // --- Salvage Logic ---
            // If the mean and sigma are within a physically reasonable window, keep them!
            double mu = fitFunc->GetParameter(1);
            double sig = fitFunc->GetParameter(2);
            
            // We can adjust these "sane" bounds later on, if needed
            if (mu > 1.105 && mu < 1.125 && sig > 0.0005 && sig < 0.008) {
                fitResults[iBin] = {mu, sig, true};
            } else {
                fitResults[iBin] = {0, 0, false};
            }
        }

        delete hNumProj; // Temporary projection; will be re-projected below in Step 5
        delete hQARing;
    }

    // =========================================================================================
    // Step 5: Discontinuous Sideband Fits
    // The key change here after new usage of TProfile2D: varR_peak is computed from hNumProj errors, NOT from totSqNum.
    // =========================================================================================
    if (printHeader) std::cout << "  -> Steps 5/6/6.5: Discontinuous sideband fits..." << std::endl; // Single print to avoid flooding screen
    for (int iBin = 1; iBin <= nBins; ++iBin){
        if (!fitResults[iBin].valid) {
            if (printHeader) std::cout << "    Bin " << iBin << ": invalid mu/sigma, skipped.\n";
            continue;
        }

        double mu = fitResults[iBin].mu;
        double sigma = fitResults[iBin].sigma;

        TH1D* hMassProj = (TH1D*)dirFits->Get(Form("hMass_Bin%d", iBin)); // Name is simple, title is descriptive though

        // MAJOR CHANGE: also retrieve the density version (written to dirFits in Step 4)
        // for use as input to PerformSimultaneousFitQA in Step 6.5.
        // hMassProjDensity = hMassProj scaled by 1/binWidth, in units of counts/GeV.
        TH1D* hMassProjDensity = (TH1D*)dirFits->Get(Form("hMassDensity_Bin%d", iBin));

        TAxis* xAxis = h2dNum->GetXaxis();
        double xLow  = xAxis->GetBinLowEdge(iBin);
        double xHigh = xAxis->GetBinUpEdge(iBin);
        // Project the corrected hNum -- this now carries proper sigma(Sum_R_i) errors:
        TH1D* hNumProj = h2dNum->ProjectionY(Form("hNumExtBin,Angle[%.2f,%.2f),Bin%d", xLow, xHigh, iBin), iBin, iBin, "e");
            // Creating a density version of this histogram (easier to save and to use):
        TH1D* hNumProjDensity = (TH1D*)hNumProj->Clone(Form("hNumExtBinDensity,Angle[%.2f,%.2f),Bin%d", xLow, xHigh, iBin));
            // Scale by bin width to get density
        hNumProjDensity->Scale(1.0, "width");
        
        if (!hMassProj || !hMassProjDensity || !hNumProj) continue;

        // --- Step 5: Discontinuous Sideband Background Fit using TGraphErrors (TH1D's are bad at that!) ---
        
        // Create a TGraphErrors to hold ONLY the sideband points
        TGraphErrors* grBkg = new TGraphErrors();
        grBkg->SetName(Form("grBkg_Bin%d", iBin));
        grBkg->SetTitle(Form("Sideband Bkg Bin %d;M_{p#pi} (GeV/c^{2});Counts/BinWidth", iBin));

        TGraphErrors* grBkgNum = new TGraphErrors();
        grBkgNum->SetName(Form("grBkgNum_Bin%d", iBin));
        grBkgNum->SetTitle(Form("Numerator Sidebands Bin %d;M_{p#pi} (GeV/c^{2});#Sigma r_{i}/BinWidth", iBin));

        int ptIdxMass = 0; // Counter for mass graph points
        int ptIdxNum = 0; // Counter for numerator graph points
            // Tracking raw counts for stability check 1.5:
        // (The new variable-width histogram generalization cannot depend on grBkg->GetY()[k]; for counts! (that is a density value in the Y axis now!))
        double rawBkgCounts = 0;
        // Loop over histogram bins and select only those in the sidebands
        for (int jBin = 1; jBin <= hMassProj->GetNbinsX(); ++jBin){
            double x = hMassProj->GetBinCenter(jBin);

            // UPDATE! Fixed to use larger bands. With 175 bins from 1.08 to 1.15 GeV/c^2, [4sigma,6sigma]
            // would mean to use only 10 bins on the right side! (same for other side!)
            // Now considering all bins in the left and right sides!
            // UPDATE 2! Now using a 6 sigma exclusion zone -- the peak was getting into the fit!
            // This uses ALL available statistics outside the 6-sigma exclusion zone:
            bool inLeftSideband = (x >= massMin && x <= (mu - 6.0 * sigma));
            bool inRightSideband = (x >= (mu + 6.0 * sigma) && x <= massMax);

            if (inLeftSideband || inRightSideband){
                double binContent = hMassProj->GetBinContent(jBin);
                double binErr = hMassProj->GetBinError(jBin);
                double bw = hMassProj->GetBinWidth(jBin);

                // 1. Process Mass Background
                if (binContent > 0.0 && binErr > 0.0) { // Safety check to avoid including empty points, or points with zero error
                    rawBkgCounts += binContent; // Add actual counts before they are converted to densities

                    // To account for variable bin widths, will now make the sideband plots actually a density instead of raw counts.
                    // This will impact the fits and should now be a correct approximation of the area under the peak.
                    grBkg->SetPoint(ptIdxMass, x, binContent / bw);
                    grBkg->SetPointError(ptIdxMass, 0.0, binErr / bw); // Standard x-error is 0 for fits, y-error is standard bin error
                    ptIdxMass++;
                }

                // 2. Process Numerator Background
                // (we can process this TGraphErrors independently and avoid discarding useful "counts hist" points
                //  if the numerator is invalid, because these TGraphErrors will only be used to get a fit function
                //  for the background! There is no need for them to have the same exact number of points!)
                double numContent = hNumProjDensity->GetBinContent(jBin);
                double numErr = hNumProjDensity->GetBinError(jBin); // hNumProj error is now sigma_R * sqrt(N) -- correct for the mean!
                if (numErr > 0.0) {
                    // Only add to the numerator graph if it has a valid, non-zero error!
                    // We do NOT check if numContent > 0 because <R> can be negative or zero.
                    grBkgNum->SetPoint(ptIdxNum, x, numContent);
                    grBkgNum->SetPointError(ptIdxNum, 0.0, numErr); 
                    ptIdxNum++;
                }
            }
        }

        // --- STABILITY CHECK 1: Enough points for a pol2? ---
        // Adding some stability checks because the 3D histogram projections have much smaller statistics to be handled!
        // A pol2 has 3 parameters. We want at least 4 points to have 1 degree of freedom.
        // Demanded 6 for stricter stability.
        if (grBkg->GetN() < 5 || grBkgNum->GetN() < 5) {
            fitResults[iBin].valid = false; // Invalidate this bin
            delete grBkg; delete grBkgNum;
            continue; // Skip fitting
        }

        // --- STABILITY CHECK 1.5: Are the sidebands just empty zeros? ---
        // If the sum of the background counts is 0, the matrix inversion will fail and we get a looooot of Minuit errors:
            // Old code to estimate the counts in the background (new code uses rawBkgCounts to account for variable bin width adaptations)
        // double sumBkg = 0;
        // for (int k = 0; k < grBkg->GetN(); ++k) {
        //     sumBkg += grBkg->GetY()[k];
        // }
        if (rawBkgCounts <= 8.0) { // 8 counts is still way to little though
            fitResults[iBin].valid = false;
            delete grBkg; delete grBkgNum;
            continue;
        }
            // A print only for the 3D histogram slicing attempts:
        // if (!printHeader) std::cout << "DEBUG! Passed stability checks and will now fit sidebands" << std::endl;

        // Fit the discontinuous graph with a pol2
        // We use pol2 to cover the small curvature in Lambda background
        TF1* bkgFitFunc = new TF1(Form("bkgFit_Bin%d", iBin), "pol2", massMin, massMax);
        TFitResultPtr rBkg = grBkg->Fit(bkgFitFunc, "Q 0 S"); // "Q" = quiet, "0" = don't draw, "S" = return TFitResultPtr

        // if (!printHeader) std::cout << "DEBUG! Mass background fit done" << std::endl; // This print revealed we lacked a stability check on grBkgNum errors on the new version!
        
        // Fit the Numerator Sidebands to get Background Polarization
            // From QA plots, you can see that the background is almost linear, maybe a bit quadratic,
            // when studying <R>(m_\Lambda). Thus, regular signal extraction should work!
        TF1* bkgNumFitFunc = new TF1(Form("bkgNumFit_Bin%d", iBin), "pol2", massMin, massMax);
        TFitResultPtr rBkgNum = grBkgNum->Fit(bkgNumFitFunc, "Q 0 S");
            // A print only for the 3D histogram slicing attempts:
        // if (!printHeader) std::cout << "DEBUG! Print after fits" << std::endl;

        // --- STABILITY CHECK 2: Did the fits converge properly? ---
        if (!rBkg->IsValid() || !rBkgNum->IsValid()) {
            fitResults[iBin].valid = false;
            delete bkgFitFunc; delete bkgNumFitFunc; // Here you still need to delete them manually, as they have not yet been taken ownership of by ROOT due to an "Add()" call!
            delete grBkg; delete grBkgNum;
            continue;
        }

        dirFits->cd();
        grBkg->GetListOfFunctions()->Add(bkgFitFunc); // Attach fit for viewing in TBrowser
        grBkg->Write();

        // Writing bkgNum histograms in a specific folder for QA:
        dirBkgNumFits->cd();
        grBkgNum->GetListOfFunctions()->Add(bkgNumFitFunc); // Attach the fit function to the graph for easy viewing in TBrowser
        grBkgNum->Write(); // Write both the sideband graph (with fit) and the raw projected histogram
        // hNumProj->Write(); // Save the corrected-error numerator for visual QA
        hNumProjDensity->Write(); // Save density to actually cross-check the fit!

        // =========================================================================================
        // Step 6: Define Signal Region and Exact Limits
        // =========================================================================================
        // if (printHeader) std::cout << "  -> Step 6: Signal extraction..." << std::endl; // Removed to avoid print flood
        int firstBin = hMassProj->FindBin(mu - 4.0 * sigma); // Now defining signal in the [mu-4sigma, mu+4sigma] region, to match Gianni's Sigma0 analysis codes
        int lastBin  = hMassProj->FindBin(mu + 4.0 * sigma);
        double x_low  = hMassProj->GetBinLowEdge(firstBin);
        double x_high = hMassProj->GetBinLowEdge(lastBin) + hMassProj->GetBinWidth(lastBin);
        // double binWidth = hMassProj->GetBinWidth(1); // No longer assume uniform binning!

        // Raw Signal Counts (Peak Region)
        // (As peak was not rescaled, this is actually counts, not a density!)
        double totCounts = 0;
        double totCountsErrSq = 0;
        double totNum = 0;
        double totNumErrSq = 0; // replaces the old totSqNum accumulation

        for (int jBin = firstBin; jBin <= lastBin; ++jBin){
            totCounts += hMassProj->GetBinContent(jBin);
            totCountsErrSq += std::pow(hMassProj->GetBinError(jBin), 2);
            
            totNum += hNumProj->GetBinContent(jBin); // Accumulate Sum_R_i
            
            // hNumProj->GetBinError(jBin) is sigma_R * sqrt(N_bin) -- correct error on the mean, <R>.
                // Bins in the peak are independent, so errors add in quadrature:
            totNumErrSq += std::pow(hNumProj->GetBinError(jBin), 2);
        }
        double errTotCounts = std::sqrt(totCountsErrSq);

        // --- STABILITY CHECK 3: Ensure there are actually any counts to divide by ---
        if (totCounts <= 0) {
            delete grBkg; delete grBkgNum; delete hNumProj;
            continue;
        }

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
        */
        
        // Background Estimation in one step
            // pol2 is fitted to counts density per mass bin, so estimating the background in the peak
            // should actually be easier now, as we have a bin-agnostic fit for the background!
            // --> TF1::Integral will give us f(m) * dm, which is an integral of a density of counts
            // per mass, over a mass interval, so the value retrieved is actually counts! (same thing
            // as dN/dpT, where integrating over a pT window will give us an area which corresponds to
            // counts, with no further efforts!)
        double bkgCounts    = bkgFitFunc->Integral(x_low, x_high);
        double errBkgCounts = bkgFitFunc->IntegralError(x_low, x_high, rBkg->GetParams(), rBkg->GetCovarianceMatrix().GetMatrixArray());

        double bkgNum = bkgNumFitFunc->Integral(x_low, x_high);
        double errBkgNum = bkgNumFitFunc->IntegralError(x_low, x_high, rBkgNum->GetParams(), rBkgNum->GetCovarianceMatrix().GetMatrixArray());

        // Final Signal Yield
        double sigCounts = totCounts - bkgCounts;
        double errSigCounts = std::sqrt(totCountsErrSq + errBkgCounts * errBkgCounts);

        if (sigCounts <= 0 || totCounts <= 0){
            std::cout << "    Bin " << iBin << ": non-positive signal, skipped.\n";
            delete grBkg; delete grBkgNum; delete hNumProj;
            continue;
        }

        // Fractions
        double fB = bkgCounts / totCounts;
        double fS = sigCounts / totCounts; // Equivalent to 1 - fB
        
        // Calculate Observables
        double R_peak = totNum / totCounts;
        double R_B = bkgNum / bkgCounts;
        
        // True Signal Polarization: R_S = (R_peak - f_B * R_B) / f_S
        double R_S = (R_peak - fB * R_B) / fS;

        // --- Error Propagation for R_peak ---
        // OLD method used totSqNum (from a separate squared-observable TH2D) to compute:
        //   varR_peak = (totSqNum/totCounts - R_peak^2) / (totCounts - 1)  [Bessel-corrected sample variance / N]
        //
        // NEW method uses the TProfile-derived errors already stored in hNumProj.
        // The uncertainty on Sum_R_i in the peak window is:
        //   sigma^2(Sum_R_i) = totNumErrSq  (sum in quadrature of sigma_R*sqrt(N_bin) over peak bins)
        // The uncertainty on the mean <R>_peak = Sum_R_i / totCounts is therefore:
        //   sigma^2(<R>_peak) = totNumErrSq / totCounts^2
        //
        // This is mathematically equivalent to the Bessel-corrected estimator but uses the
        // TProfile's properly propagated errors instead of the raw squared sum, which is
        // consistent with the corrected sideband errors.
        double varR_peak = totNumErrSq / std::pow(totCounts, 2);

        // Variance of R_B (simple error propagation for ratio B_num / B_counts)
        double varR_B = std::pow(R_B, 2) * (std::pow(errBkgNum/bkgNum, 2) + std::pow(errBkgCounts/bkgCounts, 2));
        
        // Variance of f_B
        double var_fB = std::pow(fB, 2) * (std::pow(errBkgCounts/bkgCounts, 2) + (totCountsErrSq/std::pow(totCounts, 2)));

        // Final error on R_S (formula unchanged, inputs improved)
        // sigma^2(R_S) = (1/f_S^2)*sigma^2(R_peak) + (f_B^2/f_S^2)*sigma^2(R_B) + ((R_peak - R_B)/f_S^2)^2 * sigma^2(f_B)
        double errR_S = std::sqrt(
            varR_peak / std::pow(fS, 2) +
            (std::pow(fB, 2) / std::pow(fS, 2)) * varR_B +
            std::pow((R_peak - R_B) / std::pow(fS, 2), 2) * var_fB);
        
        // --- Calculate Purity & Significance ---
        double purity = fS; 
        double errPurity = std::sqrt(var_fB); // Since fS = 1 - fB, the variance is identical

        double significance = sigCounts / std::sqrt(totCounts);
        // Standard error propagation for Z = S/sqrt(N)
        double errSignificance = std::sqrt( (errSigCounts*errSigCounts)/totCounts + (sigCounts*sigCounts*totCountsErrSq)/(4.0 * std::pow(totCounts, 3)));

        // Fill Histograms
        hSigYield->SetBinContent(iBin, sigCounts);
        hSigYield->SetBinError(iBin, errSigCounts);

        hBkgYield->SetBinContent(iBin, bkgCounts);
        hBkgYield->SetBinError(iBin, errBkgCounts);

        hRSig->SetBinContent(iBin, R_S);
        hRSig->SetBinError(iBin, errR_S);

        hRBkg->SetBinContent(iBin, R_B);
        hRBkg->SetBinError(iBin, std::sqrt(varR_B));

        hPurity->SetBinContent(iBin, purity);
        hPurity->SetBinError(iBin, errPurity);

        hSigStat->SetBinContent(iBin, significance);
        hSigStat->SetBinError(iBin, errSignificance);

        // =========================================================================
        // Step 6.5: Simultaneous Fit QA
        // =========================================================================
        // PerformSimultaneousFitQA now receives the corrected hNumProj (the TProfile fix is entirely in what hNumProj contains)
        // (changed to hNumProjDensity for proper normalization over variable bin widths. Notice that because <R>(m) = \sum_i R_i (m)/N_\Lambda(m)
        //  is a ratio, we can just do this bin width scaling without much worry)
        // init_sigYield (= sigCounts from sideband subtraction) remains a valid seed:
        // it is still an absolute count, which PerformSimultaneousFitQA converts
        // internally to a density amplitude via sigCounts / (sigma * sqrt(2*pi)).

        // if (printHeader) std::cout << "  -> Step 6.5: Simultaneous fit QA..." << std::endl; // Removed to avoid print flood
        SimFitResult qaResult = PerformSimultaneousFitQA(hMassProjDensity, hNumProjDensity, massMin, massMax, sigCounts, mu, sigma, bkgFitFunc, R_S, R_B);
        // (Notice we can still pass sigCounts to the SimultaneousFit because we will normalize the counts with hMassProj. That is independent
        // of the bin width normalizations!)
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

        // Note: We do NOT delete bkgFitFunc or bkgNumFitFunc here anymore because ROOT takes ownership 
        // when we add it to the graph's list of functions. Deleting it will cause a crash!
        delete grBkg; 
        delete grBkgNum;
        delete hNumProj;        // We own this projection
        delete hNumProjDensity; // We own this clone
    }

    dirResults->cd();
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
    // Step 8: Calculate the Angle-Integrated Ring Observable
    // =========================================================================================
    // Repeating the same procedure, but for an integrated observable!
    if (printHeader) std::cout << "  -> Step 8: Calculating Integrated Observable..." << std::endl;
    
    // Project the entire TH2D onto the Y-axis (Mass) to integrate over all angles
    TH1D* hMassInt = h2dCounts->ProjectionY(Form("hMassInt_%s", extractionName.Data()), 1, nBins, "e");
        // Project the corrected h2dNum (built from TProfile) to get proper Sum and errors:
    TH1D* hNumInt = h2dNum->ProjectionY(Form("hNumInt_%s", extractionName.Data()), 1, nBins, "e");

    // We need a place to save the final integrated value
    TH1D* hIntegratedRSig = new TH1D(Form("hIntegratedRSig_%s", extractionName.Data()), 
                                        Form("Integrated R_{S} - %s; ;R_{S}", extractionName.Data()), 
                                        1, 0, 1);

    // Only proceed if we have enough global statistics to fit a peak
    if (hMassInt->GetEntries() > 50 && hMassInt->Integral() > 0) {
        // 1. Global Mass Fit
        TF1* fitFuncInt = new TF1(Form("fitMassInt_%s", extractionName.Data()), "gaus(0) + pol2(3)", massMin, massMax);
        double maxValInt = hMassInt->GetMaximum();
        fitFuncInt->SetParameters(maxValInt, lambdaPDGMassApprox, 0.002, 0, 0, 0);
        fitFuncInt->SetParLimits(1, lambdaPDGMassApprox - 0.01, lambdaPDGMassApprox + 0.01);
        fitFuncInt->SetParLimits(2, 0.0001, 0.005);
        
        TFitResultPtr rInt = hMassInt->Fit(fitFuncInt, "Q 0 S");
        
        if (rInt->IsValid()) {
            double muInt = fitFuncInt->GetParameter(1);
            double sigmaInt = fitFuncInt->GetParameter(2);

            // 2. Global Sideband Extraction
            TGraphErrors* grBkgInt = new TGraphErrors();
            TGraphErrors* grBkgNumInt = new TGraphErrors();
            int ptIdx = 0;
            
            for (int jBin = 1; jBin <= hMassInt->GetNbinsX(); ++jBin) {
                double x = hMassInt->GetBinCenter(jBin);
                if ((x >= massMin && x <= (muInt - 6.0 * sigmaInt)) || (x >= (muInt + 6.0 * sigmaInt) && x <= massMax)) {
                    // To account for variable bin widths, will now make the sideband plots actually a density instead of raw counts.
                    // This will impact the fits and should now be a correct approximation of the area under the peak.
                    double bw = hMassInt->GetBinWidth(jBin);
                    grBkgInt->SetPoint(ptIdx, x, hMassInt->GetBinContent(jBin) / bw);
                    grBkgInt->SetPointError(ptIdx, 0.0, hMassInt->GetBinError(jBin) / bw);
                    
                    grBkgNumInt->SetPoint(ptIdx, x, hNumInt->GetBinContent(jBin) / bw);
                    grBkgNumInt->SetPointError(ptIdx, 0.0, hNumInt->GetBinError(jBin) / bw);
                    ptIdx++;
                }
            }

            if (grBkgInt->GetN() >= 4 && grBkgNumInt->GetN() >= 4) {
                TF1* bkgFitInt = new TF1(Form("bkgFitInt_%s", extractionName.Data()), "pol2", massMin, massMax);
                TF1* bkgNumFitInt = new TF1(Form("bkgNumFitInt_%s", extractionName.Data()), "pol2", massMin, massMax);
                
                TFitResultPtr rBkgInt = grBkgInt->Fit(bkgFitInt, "Q 0 S");
                TFitResultPtr rBkgNumInt = grBkgNumInt->Fit(bkgNumFitInt, "Q 0 S");

                if (rBkgInt->IsValid() && rBkgNumInt->IsValid()) {
                    
                    // 3. Global Integration
                    int firstBin = hMassInt->FindBin(muInt - 4.0 * sigmaInt);
                    int lastBin  = hMassInt->FindBin(muInt + 4.0 * sigmaInt);
                    double x_low  = hMassInt->GetBinLowEdge(firstBin);
                    double x_high = hMassInt->GetBinLowEdge(lastBin) + hMassInt->GetBinWidth(lastBin);
                    
                    double totCountsInt = 0, totCountsErrSqInt = 0;
                    double totNumInt = 0, totNumErrSqInt = 0; // Last one replaces totSqNumInt

                    for (int jBin = firstBin; jBin <= lastBin; ++jBin){
                        totCountsInt += hMassInt->GetBinContent(jBin);
                        totCountsErrSqInt += std::pow(hMassInt->GetBinError(jBin), 2);
                        totNumInt += hNumInt->GetBinContent(jBin);
                        totNumErrSqInt += std::pow(hNumInt->GetBinError(jBin), 2);
                    }

                    if (totCountsInt > 0) {
                        // Integrate the background density over the signal window.
                        // TF1::Integral returns counts directly (density * mass interval),
                        // because the new grBkgInt and grBkgNumInt give us densities instead
                        // of raw counts, so the fit is now properly bin-agnostic and the
                        // new density-fitted functions can be integrated into any mass
                        // interval and will give the proper counts (or integrated <R>_Bkg)
                        // with no further procedures other than TF1::IntegralError()!
                        // TF1::IntegralError propagates the full covariance matrix of the
                        // polynomial parameters, correctly handling their correlations.
                        // No division by binWidth: the integral is already in count units
                        // with the new definition of the TGraphErrors points as densities!
                        double bkgCountsInt = bkgFitInt->Integral(x_low, x_high);
                        double errBkgCountsInt = bkgFitInt->IntegralError(x_low, x_high, rBkgInt->GetParams(), rBkgInt->GetCovarianceMatrix().GetMatrixArray());

                        double bkgNumInt_val = bkgNumFitInt->Integral(x_low, x_high);
                        double errBkgNumInt = bkgNumFitInt->IntegralError(x_low, x_high, rBkgNumInt->GetParams(), rBkgNumInt->GetCovarianceMatrix().GetMatrixArray());

                        double sigCountsInt = totCountsInt - bkgCountsInt;

                        if (sigCountsInt > 0) {
                            double fB_Int = bkgCountsInt / totCountsInt;
                            double fS_Int = sigCountsInt / totCountsInt; 
                            
                            double R_peak_Int = totNumInt / totCountsInt;
                            double R_B_Int = (bkgCountsInt > 0) ? (bkgNumInt_val / bkgCountsInt) : 0;
                            
                            double R_S_Int = (R_peak_Int - fB_Int * R_B_Int) / fS_Int;

                            // Error Propagation (Using Unbiased Estimator!)
                            // varR_peak_Int: uses TProfile-derived errors (totNumErrSqInt), same logic as per-bin loop
                            // OLD METHOD: double varR_peak_Int = (totCountsInt > 1) ? ((totSqNumInt / totCountsInt) - (R_peak_Int * R_peak_Int)) / (totCountsInt - 1.0) : 0;
                            double varR_peak_Int = totNumErrSqInt / std::pow(totCountsInt, 2);
                            
                            double varR_B_Int = 0;
                            if (bkgCountsInt > 0 && bkgNumInt_val != 0) {
                                varR_B_Int = std::pow(R_B_Int, 2) * (std::pow(errBkgNumInt / bkgNumInt_val, 2) + std::pow(errBkgCountsInt / bkgCountsInt, 2));
                            }

                            double var_fB_Int = std::pow(fB_Int, 2) * (std::pow(errBkgCountsInt / bkgCountsInt, 2) + (totCountsErrSqInt / std::pow(totCountsInt, 2)));

                            double errR_S_Int = std::sqrt(varR_peak_Int / std::pow(fS_Int, 2) + (std::pow(fB_Int, 2) / std::pow(fS_Int, 2)) * varR_B_Int + std::pow((R_peak_Int - R_B_Int) / std::pow(fS_Int, 2), 2) * var_fB_Int);

                            // Save to histogram
                            hIntegratedRSig->SetBinContent(1, R_S_Int);
                            hIntegratedRSig->SetBinError(1, errR_S_Int);
                            
                            // std::cout << Form("    [Global] Integrated R_S = %.5f +/- %.5f", R_S_Int, errR_S_Int) << std::endl;
                        }
                    }
                }
                delete bkgFitInt; delete bkgNumFitInt;
            }
            delete grBkgInt; delete grBkgNumInt;
        }
        delete fitFuncInt;
    } // end of step 8 conditionals
    
    // Write the integrated histogram into the main directory
    dirResults->cd();
    hIntegratedRSig->Write();
        // Clean up memory for the 1D projections:
    delete hMassInt; delete hNumInt;

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
        // 2. Obtain robust seeds for the simultaneous fit using a standard pre-fit
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
        // which is perfectly fine for seeding values.
        if (rPre->IsValid()) { 
            double seed_sigma = preFit->GetParameter(2);
            // Conversion: Area = Amplitude * Sigma * sqrt(2*pi)
                // PerformSimultaneousFitQA expects to see an actual counts value for the initial seed, so we integrate (analytically) to actually get the counts estimate:
            double seed_Yield = preFit->GetParameter(0) * seed_sigma * std::sqrt(2 * TMath::Pi());
            double seed_mu    = preFit->GetParameter(1);

            // // DEBUG: print seeds
            // std::cout << "  [DEBUG " << extractionName << "] seed_Yield=" << seed_Yield
            //         << " seed_mu=" << seed_mu << " seed_sigma=" << seed_sigma
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
            TF1* bkgSeedFunc = new TF1(Form("bkgSeed_%s", extractionName.Data()), "pol2", massMin, massMax);
            bkgSeedFunc->SetParameter(0, preFit->GetParameter(3)); // c0
            bkgSeedFunc->SetParameter(1, preFit->GetParameter(4)); // c1
            bkgSeedFunc->SetParameter(2, preFit->GetParameter(5)); // c2
            // Transplant errors too, so PerformSimultaneousFitQA gets sensible step sizes
            bkgSeedFunc->SetParError(0, preFit->GetParError(3));
            bkgSeedFunc->SetParError(1, preFit->GetParError(4));
            bkgSeedFunc->SetParError(2, preFit->GetParError(5));

            // 4. Perform the Simultaneous Fit
                // The Simultaneous Fit expects a 
            SimFitResult simResTotal = PerformSimultaneousFitQA(
                hTotalMass_SimDensity, hTotalNum_SimDensity, massMin, massMax, // hTotalMass_SimDensity carries correct errors from h2dNum (which was built from TProfile2D)
                seed_Yield, seed_mu, seed_sigma, bkgSeedFunc,
                1e-3, 1e-3 // Guesses for background and signal
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

    // 5. Write to Disk
    dirResultsSim->cd(); 
    hIntegratedRSig_Sim->Write();
    hIntegratedRBkg_Sim->Write();

    delete hTotalMass_Sim;
    delete hTotalNum_Sim;

    // h2dNum was allocated by BuildNumFromProfile in this function. We remove it:
    delete h2dNum;
}

// ------------------------------------------------------------------------------------------------
// Main Macro
// ------------------------------------------------------------------------------------------------
void signalExtractionRing(const std::string& inputFilePath, const std::string& outputFolderPath){
    std::cout << "\n========================================================================" << std::endl;
    std::cout << " Starting Signal Extraction for the Ring Observable" << std::endl;
    std::cout << "========================================================================\n" << std::endl;

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

    std::string outFileName = outputFolderPath + "signalExtractionRing_" + suffix + ".root";
    std::cout << "  Input file: " << inputFilePath << std::endl;
    std::cout << "  Output file will be: " << outFileName << std::endl;

    // Open input file (Read-Only)
    TFile* inFile = TFile::Open(inputFilePath.c_str(), "READ");
    if (!inFile || inFile->IsZombie()){
        std::cerr << "  Error: Could not open input file!" << std::endl;
        return;
    }

    // Create output file (Recreate to avoid appending to old runs)
    TFile* outFile = TFile::Open(outFileName.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()){
        std::cerr << "  Error: Could not create output file!" << std::endl;
        if(inFile) inFile->Close();
        return;
    }

    // =========================================================================================
    // Step 3: Loop Over Variations & Fetch Histograms
    // =========================================================================================
    std::cout << "\n[Step 3] Looping over variations and fetching histograms..." << std::endl;
    
    std::string baseFolder = "lambdajetpolarizationionsderived"; // All results of the derived data consumer are stored inside a folder related to the task name
    std::vector<std::string> variations = {
        "Ring", 
        "RingKinematicCuts", 
        "JetKinematicCuts", 
        "JetAndLambdaKinematicCuts"
    };

    for (const auto& var : variations){
        std::cout << "\n#########################################" << std::endl;
        std::cout << "-> Processing variation: " << var << std::endl;
        std::cout << "#########################################\n" << std::endl;
        
        // Step 3.1: Output file should have these 4 cases as subfolders
        TDirectory* outDirVar = outFile->mkdir(var.c_str());
        
        std::string fullFolderPath = baseFolder + "/" + var;
        TDirectory* inDir = (TDirectory*)inFile->Get(fullFolderPath.c_str());
        
        if (!inDir){
            std::cerr << "  Warning: Directory " << fullFolderPath << " not found in input file! Skipping." << std::endl;
            continue;
        }

        // Step 3.2: Fetching Histograms
        // 1D QA & Mass
        TH1D* hMass = (TH1D*)inDir->Get("hMass");
        TH1D* hMassSigExtract = (TH1D*)inDir->Get("hMassSigExtract"); // Loaded just to get 
        // TH1D* hRingObservableMass = (TH1D*)inDir->Get("hRingObservableMass");
        
        // 2D: Observable vs Invariant Mass (Now using TProfiles to get the proper errors!)
        TProfile2D* p2dRingObservableDeltaPhiVsMass = (TProfile2D*)inDir->Get("p2dRingObservableDeltaPhiVsMass");
        TProfile2D* p2dRingObservableDeltaThetaVsMass = (TProfile2D*)inDir->Get("p2dRingObservableDeltaThetaVsMass");
        
        // Counters (denominators)
        TH2D* h2dDeltaPhiVsMass = (TH2D*)inDir->Get("h2dDeltaPhiVsMass");
        TH2D* h2dDeltaThetaVsMass = (TH2D*)inDir->Get("h2dDeltaThetaVsMass");
        
        // 3D: Observable vs Mass vs Lambda pT (Now using TProfiles to get the proper errors!)
        TProfile3D* p3dRingObservableDeltaPhiVsMassVsLambdaPt = (TProfile3D*)inDir->Get("p3dRingObservableDeltaPhiVsMassVsLambdaPt");
        TProfile3D* p3dRingObservableDeltaThetaVsMassVsLambdaPt = (TProfile3D*)inDir->Get("p3dRingObservableDeltaThetaVsMassVsLambdaPt");
            // Counters:
        TH3D* h3dDeltaPhiVsMassVsLambdaPt = (TH3D*)inDir->Get("h3dDeltaPhiVsMassVsLambdaPt");
        TH3D* h3dDeltaThetaVsMassVsLambdaPt = (TH3D*)inDir->Get("h3dDeltaThetaVsMassVsLambdaPt");

        // 3D: Observable vs Mass vs Lead Jet pT (Now using TProfiles to get the proper errors!)
        TProfile3D* p3dRingObservableDeltaPhiVsMassVsLeadJetPt = (TProfile3D*)inDir->Get("p3dRingObservableDeltaPhiVsMassVsLeadJetPt");
        TProfile3D* p3dRingObservableDeltaThetaVsMassVsLeadJetPt = (TProfile3D*)inDir->Get("p3dRingObservableDeltaThetaVsMassVsLeadJetPt"); 
            // Counters:
        TH3D* h3dDeltaPhiVsMassVsLeadJetPt = (TH3D*)inDir->Get("h3dDeltaPhiVsMassVsLeadJetPt");
        TH3D* h3dDeltaThetaVsMassVsLeadJetPt = (TH3D*)inDir->Get("h3dDeltaThetaVsMassVsLeadJetPt");

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
                            massMax);
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
                            massMax);
        // std::cout<< "DEBUG! Integral of p2dRingObservableDeltaThetaVsMass: " << p2dRingObservableDeltaThetaVsMass->Integral() << std::endl;
        // std::cout<< "DEBUG! Integral of h2dDeltaThetaVsMass: " << h2dDeltaThetaVsMass->Integral() << std::endl;
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

            // Convert the TH2D with corrected error bars into a TProfile2D so ExtractObservable2D can keep its TProfile2D signature:
            TProfile2D* p2dPhi_Lpt = ConvertToProfile2D(h2dCorrErrNumPhi_Lpt, h2dCountsPhi_Lpt, Form("p2dPhi_Lpt_%s", ptStr.Data()));

            // Send to helper!
            ExtractObservable2D(h2dCountsPhi_Lpt, p2dPhi_Lpt, dir3D_LambdaPt, Form("DeltaPhi_%s", ptStr.Data()), "#Delta#phi", massMin, massMax, false);
            delete p2dPhi_Lpt; // We own this. h2dCounts* and h2dNumCorr* are ROOT-owned

            // --- Cos Theta ---
            // Set the Z-axis ranges (Lambda Pt ranges)
            h3dDeltaThetaVsMassVsLambdaPt->GetZaxis()->SetRange(zBinMin, zBinMax);
            h3dNumCorrErrDeltaThetaVsMassVsLambdaPt->GetZaxis()->SetRange(zBinMin, zBinMax);

            TH2D* h2dCountsTheta_Lpt = (TH2D*)h3dDeltaThetaVsMassVsLambdaPt->Project3D(Form("h2dCountsTheta_Lpt_%s yx e", ptStr.Data()));
            TH2D* h2dNumCorrErrTheta_Lpt    = (TH2D*)h3dNumCorrErrDeltaThetaVsMassVsLambdaPt->Project3D(Form("h2dNumTheta_Lpt_%s yx e", ptStr.Data()));

            // Convert the TH2D with corrected error bars into a TProfile2D so ExtractObservable2D can keep its TProfile2D signature:
            TProfile2D* p2dTheta_Lpt = ConvertToProfile2D(h2dNumCorrErrTheta_Lpt, h2dCountsTheta_Lpt, Form("p2dTheta_Lpt_%s", ptStr.Data()));

            // Send to helper!
            ExtractObservable2D(h2dCountsTheta_Lpt, p2dTheta_Lpt, dir3D_LambdaPt, Form("DeltaTheta_%s", ptStr.Data()), "#Delta#theta", massMin, massMax, false);
            delete p2dTheta_Lpt;
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

            // Convert the TH2D with corrected error bars into a TProfile2D so ExtractObservable2D can keep its TProfile2D signature:
            TProfile2D* p2dPhi_Jpt = ConvertToProfile2D(h2dNumCorrErrPhi_Jpt, h2dCountsPhi_Jpt, Form("p2dPhi_Jpt_%s", ptStr.Data()));
            ExtractObservable2D(h2dCountsPhi_Jpt, p2dPhi_Jpt, dir3D_LeadJetPt, Form("DeltaPhi_%s", ptStr.Data()), "#Delta#phi", massMin, massMax, false);
            delete p2dPhi_Jpt;

            // --- Cos Theta ---
            h3dDeltaThetaVsMassVsLeadJetPt->GetZaxis()->SetRange(zBinMin, zBinMax);
            h3dNumCorrErrDeltaThetaVsMassVsLeadJetPt->GetZaxis()->SetRange(zBinMin, zBinMax);

            TH2D* h2dCountsTheta_Jpt = (TH2D*)h3dDeltaThetaVsMassVsLeadJetPt->Project3D(Form("h2dCountsTheta_Jpt_%s yx e", ptStr.Data()));
            TH2D* h2dNumCorrErrTheta_Jpt = (TH2D*)h3dNumCorrErrDeltaThetaVsMassVsLeadJetPt->Project3D(Form("h2dNumTheta_Jpt_%s yx e", ptStr.Data()));

            // Convert the TH2D with corrected error bars into a TProfile2D so ExtractObservable2D can keep its TProfile2D signature:
            TProfile2D* p2dTheta_Jpt = ConvertToProfile2D(h2dNumCorrErrTheta_Jpt, h2dCountsTheta_Jpt, Form("p2dTheta_Jpt_%s", ptStr.Data()));
            ExtractObservable2D(h2dCountsTheta_Jpt, p2dTheta_Jpt, dir3D_LeadJetPt, Form("DeltaTheta_%s", ptStr.Data()), "#Delta#theta", massMin, massMax, false);
            delete p2dTheta_Jpt;
        }

        // Cleanup -- Delete the corrected error bar TH3D objects we just created!
            // The projected TH2D pointers (h2dCountsPhi_Lpt etc.) are owned by ROOT's
            // current directory after Project3D, so we do not delete them manually.
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
        while ((key = (TKey*)nextKey())){
            TObject *obj = key->ReadObj();
            // Check if the object is a histogram (1D, 2D, or 3D)
            if (obj && obj->InheritsFrom(TH1::Class())){
                dirInput->cd();
                obj->Write(); // Write to our new output subfolder
            }
            delete obj; // Clean up memory after writing
        }
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
    
}