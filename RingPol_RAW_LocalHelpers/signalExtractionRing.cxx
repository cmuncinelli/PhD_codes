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

// Extra includes for simultaneous signal extraction fit to Ring Observable numerator and Counts denominator:
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"


constexpr double lambdaPDGMassApprox = 1.11568; // Just for some fit initial guesses

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
SimFitResult PerformSimultaneousFitQA(TH1D* hMass, TH1D* hNum, double massMin, double massMax, double init_sigYield, double init_mu, double init_sigma, TF1* bkgFitFunc, double init_RS, double init_RB){
    SimFitResult result = { // Initialized as null result
        0.0, 0.0, // R_S, err_R_S
        0.0, 0.0, // R_B, err_R_B
        0.0, 0.0, // SigYield, err_SigYield
        0.0, 0.0, // BkgYield, err_BkgYield
        0.0, 0.0, // Purity, err_Purity
        0.0, 0.0, // Significance, err_Significance
        -1        // status (Default to failed)
    };
    if (!hMass || !hNum || !bkgFitFunc) return result;

    double binW = hMass->GetBinWidth(1);

    // --- SCALING FACTOR ---
    // We normalize the mass spectrum so the peak is at 1.0.
    // This keeps Minuit parameters ~O(1) instead of O(10^5), fixing the "Status 1" convergence issues.
    double normFactor = hMass->GetMaximum();
    if (normFactor <= 0) normFactor = 1.0; 

    // 1. Define the Joint Chi2
    auto globalChi2 = [&](const double *par) {
        double chi2 = 0;
        
        for (int i = 1; i <= hMass->GetNbinsX(); ++i) {
            double m = hMass->GetBinCenter(i);
            if (m < massMin || m > massMax) continue;

            // Get Raw Data
            double y_mass_raw = hMass->GetBinContent(i);
            double e_mass_raw = hMass->GetBinError(i);
            double y_num  = hNum->GetBinContent(i); // We do NOT scale numerator data
            double e_num  = hNum->GetBinError(i);

            if (e_mass_raw <= 0) e_mass_raw = 1.0; 
            if (e_num <= 0) e_num = 1.0;

            // --- MASS COMPONENT (SCALED) ---
            // par[0] = Normalized Signal Yield (approx 1.0).
            // par[1] = mu, par[2] = sigma
            // par[3,4,5] = Normalized background coefficients
            
            double S_shape = TMath::Gaus(m, par[1], par[2], true); // Normalized Gaussian
            double S_norm  = par[0] * binW * S_shape; 
            
            double B_norm  = par[3] + par[4]*m + par[5]*m*m;
            
            double M_model_norm = S_norm + B_norm;

            // Scale data down to compare with scaled model
            double y_mass_norm = y_mass_raw / normFactor;
            double e_mass_norm = e_mass_raw / normFactor;

            chi2 += std::pow((y_mass_norm - M_model_norm) / e_mass_norm, 2);

            // --- NUMERATOR COMPONENT (UNSCALED) ---
            // The Numerator is "Weight * AbsoluteCount". 
            // We must un-scale our model by multiplying by normFactor to match the raw hNum data.
            
            // Absolute Yields = Normalized Yields * normFactor
            double S_abs = S_norm * normFactor; 
            double B_abs = B_norm * normFactor;

            // par[6] = R_S, par[7] = R_B
            double Num_model = par[6] * S_abs + par[7] * B_abs;

            chi2 += std::pow((y_num - Num_model) / e_num, 2);
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

    // 3. Initialize Parameters (SCALED)
    // We must divide the yield and background guesses by the normFactor
    min->SetVariable(0, "SigYield_Norm", init_sigYield / normFactor, 0.1);
    min->SetVariable(1, "Mu",            init_mu,       0.0001);
    min->SetVariable(2, "Sigma",         init_sigma,    0.0001);
    min->SetVariable(3, "Bkg_c0_Norm",   bkgFitFunc->GetParameter(0) / normFactor, 0.1);
    min->SetVariable(4, "Bkg_c1_Norm",   bkgFitFunc->GetParameter(1) / normFactor, 0.1);
    min->SetVariable(5, "Bkg_c2_Norm",   bkgFitFunc->GetParameter(2) / normFactor, 0.1);
    min->SetVariable(6, "R_S",           init_RS,       0.01); 
    min->SetVariable(7, "R_B",           init_RB,       0.01); 

    // Lock the kinematics tight to prevent the fit from wandering
    min->SetVariableLimits(1, init_mu - 0.003, init_mu + 0.003); // Slightly wider window
    min->SetVariableLimits(2, init_sigma * 0.5, init_sigma * 1.5);
    // Don't let yields go negative
    min->SetVariableLowerLimit(0, 0.0);

    // 4. Minimize
    min->Minimize();

    result.status = min->Status();
    if (result.status == 0) {
        const double *fitVals = min->X();
        const double *fitErrs = min->Errors();
        
        // R_S and R_B are scale-invariant!
        result.R_S = fitVals[6];
        result.err_R_S = fitErrs[6];
        result.R_B = fitVals[7];
        result.err_R_B = fitErrs[7];

        // Yields must be RESCALED back to absolute counts
        // The fit returned "Normalized Yield". Real Yield = Fit * normFactor.
        double raw_SigYield = fitVals[0] * normFactor; 
        double raw_errSigYield = fitErrs[0] * normFactor;

        // The signal yield in the +/- 4 sigma window is ~99.99366% of the total Gaussian integral
        // Using that as an estimate of the total signal yield from the fits
        result.SigYield = raw_SigYield * 0.9999366; 
        result.err_SigYield = raw_errSigYield * 0.9999366;

        // Reconstruct background polynomial to integrate it in the +/- 4 sigma window
        // (another neat way of recovering the integrated background in this region, for comparison with regular signal extraction method)
        double muVal = fitVals[1];
        double sigma = fitVals[2];
        double xLow = muVal - 4.0 * sigma; // Wrote as 4.0 * sigma to be easier to replace if I change the signal extraction region
        double xHigh = muVal + 4.0 * sigma;

        // Extract SCALED polynomial parameters
        double c0_norm = fitVals[3];
        double c1_norm = fitVals[4];
        double c2_norm = fitVals[5];

        // 1. Calculate the integration terms (derivatives wrt parameters)
        // (polynomial integration is simple enough to do by hand!)
        double dx1 = (xHigh - xLow);
        double dx2 = (xHigh * xHigh - xLow * xLow) / 2.0;
        double dx3 = (xHigh * xHigh * xHigh - xLow * xLow * xLow) / 3.0;
        
        // Analytically integrate c0 + c1*x + c2*x^2 from xLow to xHigh (Normalized Space)
        double integral_norm = (c0_norm * dx1) + (c1_norm * dx2) + (c2_norm * dx3);
        
        // Rescale Background Result to Absolute counts
        result.BkgYield = (integral_norm / binW) * normFactor;
        
        // 2. Analytically propagate the error using the Minuit Covariance Matrix
        // We calculate variance in Normalized space, then scale the final error.
        double dI_dc[3] = {dx1, dx2, dx3};
        double var_integral_norm = 0.0;
        
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                // CovMatrix indices for c0, c1, c2 are 3, 4, 5
                double cov_ij = min->CovMatrix(3 + i, 3 + j); 
                var_integral_norm += dI_dc[i] * dI_dc[j] * cov_ij;
            }
        }
        
        double err_integral_norm = std::sqrt(var_integral_norm);
        result.err_BkgYield = (err_integral_norm / binW) * normFactor;
        
        // ==========================================================
        // Calculate Purity, Significance, and their propagated errors
        // ==========================================================
        if ((result.SigYield + result.BkgYield) > 0) {
            double S = result.SigYield;
            double errS = result.err_SigYield;
            
            double B = result.BkgYield;
            double errB = result.err_BkgYield;
            
            double N = S + B;
            
            // Central values
            result.Purity = S / N;
            result.Significance = S / std::sqrt(N);
            
            // Propagated Error for Purity
            double dP_dS = B / (N * N);
            double dP_dB = -S / (N * N);
            result.err_Purity = std::sqrt( (dP_dS * dP_dS * errS * errS) + (dP_dB * dP_dB * errB * errB) );
            
            // Propagated Error for Significance
            double dSig_dS = (B + S/2.0) / std::pow(N, 1.5);
            double dSig_dB = -(S/2.0) / std::pow(N, 1.5);
            result.err_Significance = std::sqrt( (dSig_dS * dSig_dS * errS * errS) + (dSig_dB * dSig_dB * errB * errB) );
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
void ExtractObservable2D(TH2D* h2dCounts, TH2D* h2dNum, TH2D* h2dSqNum, TDirectory* parentDir, TString extractionName, TString axisTitle, double massMin, double massMax, bool printHeader = true){
    if (printHeader) std::cout << "\n[ExtractObservable2D] Starting extraction: " << extractionName << std::endl;
    // Create subdirectories for organized output
    TDirectory* dirBase = parentDir->mkdir(extractionName);
    TDirectory* dirFits = dirBase->mkdir("MassFits");
    TDirectory* dirQARingMass = dirBase->mkdir("QA_RingObservable_vs_Mass");
    TDirectory* dirBkgFits = dirBase->mkdir("BackgroundFits");
    TDirectory* dirBkgNumFits = dirBase->mkdir("BkgNumFits");
    TDirectory* dirResults = dirBase->mkdir("Results");
    TDirectory* dirResultsSim = dirBase->mkdir("ResultsCombinedFit");

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
        hMassProj->SetTitle(Form("Mass Projection Bin %d;M_{p#pi} (GeV/c^{2});Counts", iBin));

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
        
        // Smart initialization (starting parameters should be somewhat physical)
        double maxCount = hMassProj->GetMaximum();
        fitFunc->SetParameter(0, maxCount);        // Constant
        fitFunc->SetParameter(1, lambdaPDGMassApprox);         // Mean close to PDG value
        fitFunc->SetParameter(2, 0.002);           // Sigma (typical for Lambda, derived from some earlier analysis of mine)

        // 4. Apply Stability Constraints
        fitFunc->SetParameter(0, maxCount);        // Constant (Height)
        // - Constrain mean within mass axis limits
        fitFunc->SetParLimits(1, lambdaPDGMassApprox-0.01, lambdaPDGMassApprox+0.01); // Constraints close to mass (1.15 to 1.125)
        // Sigma: Restrict to physical limits (e.g., 0.1 to 5 MeV)
        fitFunc->SetParLimits(2, 0.0001, 0.005);

        // 5. Fit using Binned Poisson Log-Likelihood ("L"), Quiet ("Q"), No Draw ("0"), Respect Range ("R")
        // We use "S" to get the TFitResultPtr to check if the fit converged
        TFitResultPtr r = hMassProj->Fit(fitFunc, "Q 0 R S"); // Removed log-likelihood to make fits behave better in high-statistics 1D environment!

        dirFits->cd();
        hMassProj->GetListOfFunctions()->Add(fitFunc); // Attach fit to histogram for saving later
            // Save the histogram only if the projection is not empty!
        if (hMassProj->GetEntries() > 0) hMassProj->Write(); // For QA

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
    }

    // =========================================================================================
    // Step 5: Discontinuous Sideband Fits
    // =========================================================================================
    if (printHeader) std::cout << "  -> Steps 5/6/6.5: Discontinuous sideband fits..." << std::endl; // Single print to avoid flooding screen
    for (int iBin = 1; iBin <= nBins; ++iBin){
        if (!fitResults[iBin].valid) {
            if (printHeader) std::cout << "    Bin " << iBin << ": invalid mu/sigma, skipped.\n";
            continue;
        }

        double mu = fitResults[iBin].mu;
        double sigma = fitResults[iBin].sigma;

        // Re-fetch the mass projections
        TH1D* hMassProj = (TH1D*)dirFits->Get(Form("hMass_Bin%d", iBin));
        TH1D* hNumProj = h2dNum->ProjectionY(Form("hNum_Ext_Bin%d", iBin), iBin, iBin, "e");
        // Also fetch the squared observable projection to get the variance of R_peak
        TH1D* hSqNumProj = h2dSqNum->ProjectionY(Form("hSqNum_Ext_Bin%d", iBin), iBin, iBin, "e");
        
        if (!hMassProj || !hNumProj || !hSqNumProj) continue;

        // --- Step 5: Discontinuous Sideband Background Fit using TGraphErrors (TH1D's are bad at that!) ---
        
        // Create a TGraphErrors to hold ONLY the sideband points
        TGraphErrors* grBkg = new TGraphErrors();
        grBkg->SetName(Form("grBkg_Bin%d", iBin));
        grBkg->SetTitle(Form("Sideband Bkg Bin %d;M_{p#pi} (GeV/c^{2});Counts", iBin));

        TGraphErrors* grBkgNum = new TGraphErrors();
        grBkgNum->SetName(Form("grBkgNum_Bin%d", iBin));
        grBkgNum->SetTitle(Form("Numerator Sidebands Bin %d;M_{p#pi} (GeV/c^{2});#Sigma r_{i}", iBin));

        int ptIdx = 0;
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
                grBkg->SetPoint(ptIdx, x, hMassProj->GetBinContent(jBin));
                grBkg->SetPointError(ptIdx, 0.0, hMassProj->GetBinError(jBin)); // Standard x-error is 0 for fits, y-error is standard bin error
                
                grBkgNum->SetPoint(ptIdx, x, hNumProj->GetBinContent(jBin));
                grBkgNum->SetPointError(ptIdx, 0.0, hNumProj->GetBinError(jBin));
                ptIdx++;
            }
        }

        // --- STABILITY CHECK 1: Enough points for a pol2? ---
        // NEW: Adding some stability checks because the 3D histogram projections have much smaller statistics to be handled!
        // A pol2 has 3 parameters. We want at least 4 points to have 1 degree of freedom.
        if (grBkg->GetN() < 4 || grBkgNum->GetN() < 4) {
            fitResults[iBin].valid = false; // Invalidate this bin
            delete grBkg; delete grBkgNum;
            continue; // Skip fitting
        }

        // --- STABILITY CHECK 1.5: Are the sidebands just empty zeros? ---
        // If the sum of the background counts is 0, the matrix inversion will fail and we get a looooot of Minuit errors:
        double sumBkg = 0;
        for (int k = 0; k < grBkg->GetN(); ++k) {
            sumBkg += grBkg->GetY()[k];
        }
        if (sumBkg <= 5) { // 5 counts is still way to little
            fitResults[iBin].valid = false;
            delete grBkg; delete grBkgNum;
            continue;
        }

        // Fit the discontinuous graph with a pol2
        // We use pol2 to cover the small curvature in Lambda background
        TF1* bkgFitFunc = new TF1(Form("bkgFit_Bin%d", iBin), "pol2", massMin, massMax);
        TFitResultPtr rBkg = grBkg->Fit(bkgFitFunc, "Q 0 S"); // "Q" = quiet, "0" = don't draw, "S" = return TFitResultPtr
        
        // Fit the Numerator Sidebands to get Background Polarization
            // From QA plots, you can see that the background is almost linear, maybe a bit quadratic,
            // when studying <R>(m_\Lambda). Thus, regular signal extraction should work!
        TF1* bkgNumFitFunc = new TF1(Form("bkgNumFit_Bin%d", iBin), "pol2", massMin, massMax);
        TFitResultPtr rBkgNum = grBkgNum->Fit(bkgNumFitFunc, "Q 0 S");

        // --- STABILITY CHECK 2: Did the fits converge properly? ---
        if (!rBkg->IsValid() || !rBkgNum->IsValid()) {
            fitResults[iBin].valid = false;
            delete bkgFitFunc; delete bkgNumFitFunc; // Here you still need to delete them manually, as they have not yet been taken ownership of by ROOT due to an "Add()" call!
            delete grBkg; delete grBkgNum;
            continue;
        }

        dirBkgFits->cd();
        grBkg->GetListOfFunctions()->Add(bkgFitFunc); // Attach fit for viewing in TBrowser
        grBkg->Write();

        // Writing bkgNum histograms in a specific folder for QA:
        dirBkgNumFits->cd();
        grBkgNum->GetListOfFunctions()->Add(bkgNumFitFunc); // Attach the fit function to the graph for easy viewing in TBrowser
        grBkgNum->Write(); // Write both the sideband graph (with fit) and the raw projected histogram
        hNumProj->Write();

        // =========================================================================================
        // Step 6: Define Signal Region and Exact Limits
        // =========================================================================================
        // if (printHeader) std::cout << "  -> Step 6: Signal extraction..." << std::endl; // Removed to avoid print flood
        int firstBin = hMassProj->FindBin(mu - 4.0 * sigma); // Now defining signal in the [mu-4sigma, mu+4sigma] region, to match Gianni's Sigma0 analysis codes
        int lastBin  = hMassProj->FindBin(mu + 4.0 * sigma);
        double x_low  = hMassProj->GetBinLowEdge(firstBin);
        double x_high = hMassProj->GetBinLowEdge(lastBin) + hMassProj->GetBinWidth(lastBin);
        double binWidth = hMassProj->GetBinWidth(1); // Assuming uniform binning

        // Raw Signal Counts (Peak Region)
        double totCounts = 0;
        double totCountsErrSq = 0;
        double totNum = 0;
        double totNumErrSq = 0;
        double totSqNum = 0; // For variance of the observable

        for (int jBin = firstBin; jBin <= lastBin; ++jBin){
            totCounts += hMassProj->GetBinContent(jBin);
            totCountsErrSq += std::pow(hMassProj->GetBinError(jBin), 2);
            
            totNum += hNumProj->GetBinContent(jBin);
            totNumErrSq += std::pow(hNumProj->GetBinError(jBin), 2);
            
            totSqNum += hSqNumProj->GetBinContent(jBin);
        }
        double errTotCounts = std::sqrt(totCountsErrSq);

        // --- STABILITY CHECK 3: Ensure there are actually any counts to divide by ---
        if (totCounts <= 0) {
            delete grBkg; delete grBkgNum;
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
        double bkgCounts = bkgFitFunc->Integral(x_low, x_high) / binWidth;
        double errBkgCounts = bkgFitFunc->IntegralError(x_low, x_high, rBkg->GetParams(), rBkg->GetCovarianceMatrix().GetMatrixArray()) / binWidth;

        double bkgNum = bkgNumFitFunc->Integral(x_low, x_high) / binWidth;
        double errBkgNum = bkgNumFitFunc->IntegralError(x_low, x_high, rBkgNum->GetParams(), rBkgNum->GetCovarianceMatrix().GetMatrixArray()) / binWidth;

        // Final Signal Yield
        double sigCounts = totCounts - bkgCounts;
        double errSigCounts = std::sqrt(totCountsErrSq + errBkgCounts * errBkgCounts);

        if (sigCounts <= 0 || totCounts <= 0){
            std::cout << "    Bin " << iBin << ": non-positive signal, skipped.\n";
            delete grBkg; delete grBkgNum;
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

        // --- Error Propagation for Polarization ---
        // Variance of R_peak (using squared observable data)
        // double varR_peak = ((totSqNum / totCounts) - (R_peak * R_peak)) / totCounts; // Standard error of the mean
           // Variance of R_peak (using the unbiased estimator with Bessel's correction for unbiased estimators)
        double varR_peak = ((totSqNum / totCounts) - (R_peak * R_peak)) / (totCounts - 1.0);

        // Variance of R_B (simple error propagation for ratio B_num / B_counts)
        double varR_B = std::pow(R_B, 2) * (std::pow(errBkgNum/bkgNum, 2) + std::pow(errBkgCounts/bkgCounts, 2));
        
        // Variance of f_B
        double var_fB = std::pow(fB, 2) * (std::pow(errBkgCounts/bkgCounts, 2) + (totCountsErrSq/std::pow(totCounts, 2)));

        // Final Variance Formula:
        // sigma^2(R_S) = (1/f_S^2)*sigma^2(R_peak) + (f_B^2/f_S^2)*sigma^2(R_B) + ((R_peak - R_B)/f_S^2)^2 * sigma^2(f_B)
        double errR_S = std::sqrt(varR_peak/std::pow(fS, 2) + (std::pow(fB, 2)/std::pow(fS, 2))*varR_B + std::pow((R_peak - R_B)/std::pow(fS, 2), 2)*var_fB);
        
        // --- Calculate Purity & Significance ---
        double purity = fS; 
        double errPurity = std::sqrt(var_fB); // Since fS = 1 - fB, the variance is identical

        double significance = sigCounts / std::sqrt(totCounts);
        // Standard error propagation for Z = S/sqrt(N)
        double errSignificance = std::sqrt( (errSigCounts*errSigCounts)/totCounts + (sigCounts*sigCounts*totCountsErrSq)/(4.0*std::pow(totCounts, 3)) );

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
        // if (printHeader) std::cout << "  -> Step 6.5: Simultaneous fit QA..." << std::endl; // Removed to avoid print flood
        SimFitResult qaResult = PerformSimultaneousFitQA(hMassProj, hNumProj, massMin, massMax, sigCounts, mu, sigma, bkgFitFunc, R_S, R_B);
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
    TH1D* hMassInt  = h2dCounts->ProjectionY(Form("hMassInt_%s", extractionName.Data()), 1, nBins, "e");
    TH1D* hNumInt   = h2dNum->ProjectionY(Form("hNumInt_%s", extractionName.Data()), 1, nBins, "e");
    TH1D* hSqNumInt = h2dSqNum->ProjectionY(Form("hSqNumInt_%s", extractionName.Data()), 1, nBins, "e");

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
                    grBkgInt->SetPoint(ptIdx, x, hMassInt->GetBinContent(jBin));
                    grBkgInt->SetPointError(ptIdx, 0.0, hMassInt->GetBinError(jBin));
                    
                    grBkgNumInt->SetPoint(ptIdx, x, hNumInt->GetBinContent(jBin));
                    grBkgNumInt->SetPointError(ptIdx, 0.0, hNumInt->GetBinError(jBin));
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
                    double binWidth = hMassInt->GetBinWidth(1);
                    
                    double totCountsInt = 0, totCountsErrSqInt = 0;
                    double totNumInt = 0, totSqNumInt = 0;

                    for (int jBin = firstBin; jBin <= lastBin; ++jBin){
                        totCountsInt += hMassInt->GetBinContent(jBin);
                        totCountsErrSqInt += std::pow(hMassInt->GetBinError(jBin), 2);
                        totNumInt += hNumInt->GetBinContent(jBin);
                        totSqNumInt += hSqNumInt->GetBinContent(jBin);
                    }

                    if (totCountsInt > 0 && binWidth > 0) {
                        double bkgCountsInt = bkgFitInt->Integral(x_low, x_high) / binWidth;
                        double errBkgCountsInt = bkgFitInt->IntegralError(x_low, x_high, rBkgInt->GetParams(), rBkgInt->GetCovarianceMatrix().GetMatrixArray()) / binWidth;

                        double bkgNumInt_val = bkgNumFitInt->Integral(x_low, x_high) / binWidth;
                        double errBkgNumInt = bkgNumFitInt->IntegralError(x_low, x_high, rBkgNumInt->GetParams(), rBkgNumInt->GetCovarianceMatrix().GetMatrixArray()) / binWidth;

                        double sigCountsInt = totCountsInt - bkgCountsInt;

                        if (sigCountsInt > 0) {
                            double fB_Int = bkgCountsInt / totCountsInt;
                            double fS_Int = sigCountsInt / totCountsInt; 
                            
                            double R_peak_Int = totNumInt / totCountsInt;
                            double R_B_Int = (bkgCountsInt > 0) ? (bkgNumInt_val / bkgCountsInt) : 0;
                            
                            double R_S_Int = (R_peak_Int - fB_Int * R_B_Int) / fS_Int;

                            // Error Propagation (Using Unbiased Estimator!)
                            double varR_peak_Int = (totCountsInt > 1) ? ((totSqNumInt / totCountsInt) - (R_peak_Int * R_peak_Int)) / (totCountsInt - 1.0) : 0;
                            
                            double varR_B_Int = 0;
                            if (bkgCountsInt > 0) {
                                varR_B_Int = (std::pow(errBkgNumInt, 2) / std::pow(bkgCountsInt, 2)) + 
                                                (std::pow(bkgNumInt_val, 2) * std::pow(errBkgCountsInt, 2) / std::pow(bkgCountsInt, 4));
                            }
                            
                            double var_fB_Int = (std::pow(errBkgCountsInt, 2) / std::pow(totCountsInt, 2)) + 
                                                (std::pow(bkgCountsInt, 2) * totCountsErrSqInt / std::pow(totCountsInt, 4));

                            double errR_S_Int = std::sqrt(varR_peak_Int/std::pow(fS_Int, 2) + (std::pow(fB_Int, 2)/std::pow(fS_Int, 2))*varR_B_Int + std::pow((R_peak_Int - R_B_Int)/std::pow(fS_Int, 2), 2)*var_fB_Int);

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

    // (Clean up memory for the 1D projections)
    delete hMassInt; delete hNumInt; delete hSqNumInt;

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
        
        double maxVal = hTotalMass_Sim->GetMaximum();
        double estimatedMu = hTotalMass_Sim->GetBinCenter(hTotalMass_Sim->GetMaximumBin());
        
        // Estimate background from the edges (average of first and last bin)
        // (we NEED this type of initial guess to be better for the integrated fit because the statistics is just so much higher!)
        double firstBinC = hTotalMass_Sim->GetBinContent(1);
        double lastBinC  = hTotalMass_Sim->GetBinContent(hTotalMass_Sim->GetNbinsX());
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
        TFitResultPtr rPre = hTotalMass_Sim->Fit(preFit, "Q N 0 R S");

        // We check IsValid() instead of strictly == 0, as high stats sometimes gives status 4000 (Converged but non-pos-def covariance)
        // which is perfectly fine for seeding values.
        if (rPre->IsValid()) { 
            double seed_sigma = preFit->GetParameter(2);
            // Conversion: Area = Amplitude * Sigma * sqrt(2*pi)
            double seed_Yield = preFit->GetParameter(0) * seed_sigma * std::sqrt(2 * TMath::Pi());
            double seed_mu    = preFit->GetParameter(1);

            // 4. Perform the Simultaneous Fit
            SimFitResult simResTotal = PerformSimultaneousFitQA(
                hTotalMass_Sim, hTotalNum_Sim, massMin, massMax,
                seed_Yield, seed_mu, seed_sigma, preFit,
                0.0, 0.0
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
    }

    // 5. Write to Disk
    dirResultsSim->cd(); 
    hIntegratedRSig_Sim->Write();
    hIntegratedRBkg_Sim->Write();

    delete hTotalMass_Sim;
    delete hTotalNum_Sim;
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
        TH1D* hMassSigExtract = (TH1D*)inDir->Get("hMassSigExtract");
        TH1D* hRingObservableMass = (TH1D*)inDir->Get("hRingObservableMass");
        
        // 2D: Observable vs Invariant Mass
        TH2D* h2dRingObservableDeltaPhiVsMass = (TH2D*)inDir->Get("h2dRingObservableDeltaPhiVsMass");
        TH2D* h2dRingObservableDeltaThetaVsMass = (TH2D*)inDir->Get("h2dRingObservableDeltaThetaVsMass");
        
        // Squared observable (for variance propagation)
        TH2D* h2dRingObservableSquaredDeltaPhiVsMass = (TH2D*)inDir->Get("h2dRingObservableSquaredDeltaPhiVsMass");
        TH2D* h2dRingObservableSquaredDeltaThetaVsMass = (TH2D*)inDir->Get("h2dRingObservableSquaredDeltaThetaVsMass");
        
        // Counters (denominators)
        TH2D* h2dDeltaPhiVsMass = (TH2D*)inDir->Get("h2dDeltaPhiVsMass");
        TH2D* h2dDeltaThetaVsMass = (TH2D*)inDir->Get("h2dDeltaThetaVsMass");
        
        // 3D: Observable vs Mass vs Lambda pT
        TH3D* h3dRingObservableDeltaPhiVsMassVsLambdaPt = (TH3D*)inDir->Get("h3dRingObservableDeltaPhiVsMassVsLambdaPt");
        TH3D* h3dRingObservableDeltaThetaVsMassVsLambdaPt = (TH3D*)inDir->Get("h3dRingObservableDeltaThetaVsMassVsLambdaPt");
        
        TH3D* h3dRingObservableSquaredDeltaPhiVsMassVsLambdaPt = (TH3D*)inDir->Get("h3dRingObservableSquaredDeltaPhiVsMassVsLambdaPt");
        TH3D* h3dRingObservableSquaredDeltaThetaVsMassVsLambdaPt = (TH3D*)inDir->Get("h3dRingObservableSquaredDeltaThetaVsMassVsLambdaPt");
        
        TH3D* h3dDeltaPhiVsMassVsLambdaPt = (TH3D*)inDir->Get("h3dDeltaPhiVsMassVsLambdaPt");
        TH3D* h3dDeltaThetaVsMassVsLambdaPt = (TH3D*)inDir->Get("h3dDeltaThetaVsMassVsLambdaPt");

        // 3D: Observable vs Mass vs Lead Jet pT
        TH3D* h3dRingObservableDeltaPhiVsMassVsLeadJetPt = (TH3D*)inDir->Get("h3dRingObservableDeltaPhiVsMassVsLeadJetPt");
        TH3D* h3dRingObservableDeltaThetaVsMassVsLeadJetPt = (TH3D*)inDir->Get("h3dRingObservableDeltaThetaVsMassVsLeadJetPt");
        // --- Squared version ---
        TH3D* h3dRingObservableSquaredDeltaPhiVsMassVsLeadJetPt = (TH3D*)inDir->Get("h3dRingObservableSquaredDeltaPhiVsMassVsLeadJetPt");
        TH3D* h3dRingObservableSquaredDeltaThetaVsMassVsLeadJetPt = (TH3D*)inDir->Get("h3dRingObservableSquaredDeltaThetaVsMassVsLeadJetPt");
        // --- Counters ---
        TH3D* h3dDeltaPhiVsMassVsLeadJetPt = (TH3D*)inDir->Get("h3dDeltaPhiVsMassVsLeadJetPt");
        TH3D* h3dDeltaThetaVsMassVsLeadJetPt = (TH3D*)inDir->Get("h3dDeltaThetaVsMassVsLeadJetPt");

        // Basic check to ensure critical histograms loaded correctly
        if (!h2dRingObservableDeltaPhiVsMass || !hMassSigExtract){
            std::cerr << "  Error: Missing essential histograms for " << var << ". Skipping this variation." << std::endl;
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
                            h2dRingObservableDeltaPhiVsMass, 
                            h2dRingObservableSquaredDeltaPhiVsMass, 
                            outDirVar,       // The parent TDirectory to save everything inside
                            "DeltaPhi",      // This creates a "DeltaPhi" subfolder for organized output
                            "#Delta#phi",    // Axis title for generated histograms
                            massMin, 
                            massMax);

        // -----------------------------------------------------------------------------------------
        // Extraction for Cos Theta Bins
        // -----------------------------------------------------------------------------------------
        std::cout << "    -> Processing Cos Theta observables..." << std::endl;
        ExtractObservable2D(h2dDeltaThetaVsMass, 
                            h2dRingObservableDeltaThetaVsMass, 
                            h2dRingObservableSquaredDeltaThetaVsMass, 
                            outDirVar,       // The parent TDirectory to save everything inside
                            "DeltaTheta",    // This creates a "DeltaTheta" subfolder for organized output
                            "cos(#theta)",   // Axis title for generated histograms
                            massMin, 
                            massMax);

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
            h3dRingObservableDeltaPhiVsMassVsLambdaPt->GetZaxis()->SetRange(zBinMin, zBinMax);
            h3dRingObservableSquaredDeltaPhiVsMassVsLambdaPt->GetZaxis()->SetRange(zBinMin, zBinMax);

            // Project to 2D using "yx e" (Y vs X, which places Mass on the Y-axis and Angle on the X-axis)
            // ROOT's Project3D allows naming the output histogram by putting the name followed by a space before the option
            TH2D* h2dCountsPhi_Lpt = (TH2D*)h3dDeltaPhiVsMassVsLambdaPt->Project3D(Form("h2dCountsPhi_Lpt_%s yx e", ptStr.Data()));
            TH2D* h2dNumPhi_Lpt    = (TH2D*)h3dRingObservableDeltaPhiVsMassVsLambdaPt->Project3D(Form("h2dNumPhi_Lpt_%s yx e", ptStr.Data()));
            TH2D* h2dSqNumPhi_Lpt  = (TH2D*)h3dRingObservableSquaredDeltaPhiVsMassVsLambdaPt->Project3D(Form("h2dSqNumPhi_Lpt_%s yx e", ptStr.Data()));

            // Send to helper!
            ExtractObservable2D(h2dCountsPhi_Lpt, h2dNumPhi_Lpt, h2dSqNumPhi_Lpt, dir3D_LambdaPt, Form("DeltaPhi_%s", ptStr.Data()), "#Delta#phi", massMin, massMax, false);

            // --- Cos Theta ---
            // Set the Z-axis ranges (Lambda Pt ranges)
            h3dDeltaThetaVsMassVsLambdaPt->GetZaxis()->SetRange(zBinMin, zBinMax);
            h3dRingObservableDeltaThetaVsMassVsLambdaPt->GetZaxis()->SetRange(zBinMin, zBinMax);
            h3dRingObservableSquaredDeltaThetaVsMassVsLambdaPt->GetZaxis()->SetRange(zBinMin, zBinMax);

            TH2D* h2dCountsTheta_Lpt = (TH2D*)h3dDeltaThetaVsMassVsLambdaPt->Project3D(Form("h2dCountsTheta_Lpt_%s yx e", ptStr.Data()));
            TH2D* h2dNumTheta_Lpt    = (TH2D*)h3dRingObservableDeltaThetaVsMassVsLambdaPt->Project3D(Form("h2dNumTheta_Lpt_%s yx e", ptStr.Data()));
            TH2D* h2dSqNumTheta_Lpt  = (TH2D*)h3dRingObservableSquaredDeltaThetaVsMassVsLambdaPt->Project3D(Form("h2dSqNumTheta_Lpt_%s yx e", ptStr.Data()));

            // Send to helper!
            ExtractObservable2D(h2dCountsTheta_Lpt, h2dNumTheta_Lpt, h2dSqNumTheta_Lpt, dir3D_LambdaPt, Form("DeltaTheta_%s", ptStr.Data()), "#Delta#theta", massMin, massMax, false);
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
            h3dRingObservableDeltaPhiVsMassVsLeadJetPt->GetZaxis()->SetRange(zBinMin, zBinMax);
            h3dRingObservableSquaredDeltaPhiVsMassVsLeadJetPt->GetZaxis()->SetRange(zBinMin, zBinMax);

            TH2D* h2dCountsPhi_Jpt = (TH2D*)h3dDeltaPhiVsMassVsLeadJetPt->Project3D(Form("h2dCountsPhi_Jpt_%s yx e", ptStr.Data()));
            TH2D* h2dNumPhi_Jpt    = (TH2D*)h3dRingObservableDeltaPhiVsMassVsLeadJetPt->Project3D(Form("h2dNumPhi_Jpt_%s yx e", ptStr.Data()));
            TH2D* h2dSqNumPhi_Jpt  = (TH2D*)h3dRingObservableSquaredDeltaPhiVsMassVsLeadJetPt->Project3D(Form("h2dSqNumPhi_Jpt_%s yx e", ptStr.Data()));

            ExtractObservable2D(h2dCountsPhi_Jpt, h2dNumPhi_Jpt, h2dSqNumPhi_Jpt, dir3D_LeadJetPt, Form("DeltaPhi_%s", ptStr.Data()), "#Delta#phi", massMin, massMax, false);

            // --- Cos Theta ---
            h3dDeltaThetaVsMassVsLeadJetPt->GetZaxis()->SetRange(zBinMin, zBinMax);
            h3dRingObservableDeltaThetaVsMassVsLeadJetPt->GetZaxis()->SetRange(zBinMin, zBinMax);
            h3dRingObservableSquaredDeltaThetaVsMassVsLeadJetPt->GetZaxis()->SetRange(zBinMin, zBinMax);

            TH2D* h2dCountsTheta_Jpt = (TH2D*)h3dDeltaThetaVsMassVsLeadJetPt->Project3D(Form("h2dCountsTheta_Jpt_%s yx e", ptStr.Data()));
            TH2D* h2dNumTheta_Jpt    = (TH2D*)h3dRingObservableDeltaThetaVsMassVsLeadJetPt->Project3D(Form("h2dNumTheta_Jpt_%s yx e", ptStr.Data()));
            TH2D* h2dSqNumTheta_Jpt  = (TH2D*)h3dRingObservableSquaredDeltaThetaVsMassVsLeadJetPt->Project3D(Form("h2dSqNumTheta_Jpt_%s yx e", ptStr.Data()));

            ExtractObservable2D(h2dCountsTheta_Jpt, h2dNumTheta_Jpt, h2dSqNumTheta_Jpt, dir3D_LeadJetPt, Form("DeltaTheta_%s", ptStr.Data()), "#Delta#theta", massMin, massMax, false);
        }

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