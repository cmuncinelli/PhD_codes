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

constexpr double lambdaPDGMassApprox = 1.11568;

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
        std::cout << "\n  -> Processing variation: " << var << std::endl;
        
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

        // Basic check to ensure critical histograms loaded correctly
        if (!h2dRingObservableDeltaPhiVsMass || !hMassSigExtract){
            std::cerr << "  Error: Missing essential histograms for " << var << ". Skipping this variation." << std::endl;
            continue;
        }

        // =========================================================================================
        // Step 4: Mass Projections, Fits, and QA Plots for each Angular Bin
        // =========================================================================================
        std::cout << "  [Step 4] Projecting mass bins, fitting mu/sigma, and building QA plots..." << std::endl;

        // Create subdirectories for organized output
        TDirectory* dirFitsPhi = outDirVar->mkdir("MassFits_DeltaPhi");
        TDirectory* dirFitsTheta = outDirVar->mkdir("MassFits_DeltaTheta");
        TDirectory* dirQARingMass = outDirVar->mkdir("QA_RingObservable_vs_Mass");

        // We will need to store the extracted mu and sigma for Step 5
        // Using a simple struct to hold the fit results for each bin
        struct FitResult { double mu; double sigma; bool valid; };
        std::vector<FitResult> fitResultsPhi(h2dDeltaPhiVsMass->GetNbinsX() + 1);
        std::vector<FitResult> fitResultsTheta(h2dDeltaThetaVsMass->GetNbinsX() + 1);

        double massMin = hMassSigExtract->GetXaxis()->GetXmin();
        double massMax = hMassSigExtract->GetXaxis()->GetXmax();
        double massRange = massMax - massMin;

        // -----------------------------------------------------------------------------------------
        // 4.A: Delta Phi Bins
        // -----------------------------------------------------------------------------------------
        for (int iBin = 1; iBin <= h2dDeltaPhiVsMass->GetNbinsX(); ++iBin){
            // 1. Project the mass spectrum (Counts/Denominator) using "e" for error propagation
            TString projName = Form("hMass_DeltaPhi_Bin%d", iBin);
            TH1D* hMassProj = h2dDeltaPhiVsMass->ProjectionY(projName, iBin, iBin, "e");
            hMassProj->SetTitle(Form("Mass Projection #Delta#phi Bin %d;M_{p#pi} (GeV/c^{2});Counts", iBin));

            // 2. Project the Numerator to calculate <R>(m) QA plot
            TString numProjName = Form("hNum_DeltaPhi_Bin%d", iBin);
            TH1D* hNumProj = h2dRingObservableDeltaPhiVsMass->ProjectionY(numProjName, iBin, iBin, "e");
            
            TString qaName = Form("hQARingVsMass_DeltaPhi_Bin%d", iBin);
            TH1D* hQARing = (TH1D*)hNumProj->Clone(qaName);
            hQARing->SetTitle(Form("<R> vs Mass for #Delta#phi Bin %d;M_{p#pi} (GeV/c^{2});<R>", iBin));
            hQARing->Divide(hNumProj, hMassProj, 1.0, 1.0, "B"); // Binomial errors might be tricky here, standard divide is okay for now as it's weighted
            
            dirQARingMass->cd();
            hQARing->Write();

            // 3. Setup the Fit Function (gaus + pol2)
            TString fitName = Form("fit_DeltaPhi_Bin%d", iBin);
            TF1* fitFunc = new TF1(fitName, "gaus(0) + pol2(3)", massMin, massMax); // Short-hand initialization of fit function
            
            // Smart initialization (starting parameters should be somewhat physicaç)
            double maxBinCenter = hMassProj->GetBinCenter(hMassProj->GetMaximumBin());
            double maxCount = hMassProj->GetMaximum();
            fitFunc->SetParameter(0, maxCount);        // Constant
            fitFunc->SetParameter(1, lambdaPDGMassApprox);         // Mean close to PDG value
            fitFunc->SetParameter(2, 0.002);           // Sigma (typical for Lambda, derived from some earlier analysis of mine)

            // 4. Apply Stability Constraints
            fitFunc->SetParameter(0, maxCount);        // Constant (Height)
            // - Constrain mean within mass axis limits
            fitFunc->SetParLimits(1, lambdaPDGMassApprox-0.01, massMax+0.01); // Constraints close to mass (1.15 to 1.125)
            // Sigma: Restrict to physical limits (e.g., 0.1 to 5 MeV)
            fitFunc->SetParLimits(2, 0.0001, 0.005);

            // 5. Fit using Binned Poisson Log-Likelihood ("L"), Quiet ("Q"), No Draw ("0"), Respect Range ("R")
            // We use "S" to get the TFitResultPtr to check if the fit converged
            TFitResultPtr r = hMassProj->Fit(fitFunc, "Q 0 R S"); // Removed log-likelihood to make fits behave better in high-statistics 1D environment!

            dirFitsPhi->cd();
            hMassProj->GetListOfFunctions()->Add(fitFunc); // Attach fit to histogram for saving later
            hMassProj->Write(); // For QA

            // Store results if fit is somewhat reasonable
            int fitStatus = int(r); // Extract the integer status code from the fit result
            
            // 0 = perfectly converged. 
            if (fitStatus == 0) {
                fitResultsPhi[iBin] = {fitFunc->GetParameter(1), fitFunc->GetParameter(2), true};
            }
            // Most of the other errors in this stage are simply error code 1, which means Minuit had problems with
            // the covariance matrix (the parameter "error" matrix). As we don't care about errors right now, this
            // is perfectly fine and we can just proceed. Thus, added a "salvage" block that just checks if the
            // fit parameters are somewhat OK before proceeding.
            else {
                // std::cerr << "    [DEBUG] Fit 'failed' for DeltaPhi Bin " << iBin << " (Status: " << fitStatus << ")" << std::endl;
                // std::cerr << "      -> Valid? " << (r->IsValid() ? "Yes" : "No") << std::endl;
                
                // // Loop through parameters and print their values vs limits
                // for (int p = 0; p < fitFunc->GetNpar(); ++p) {
                //     double val = fitFunc->GetParameter(p);
                //     double err = fitFunc->GetParError(p);
                //     double minLim, maxLim;
                //     fitFunc->GetParLimits(p, minLim, maxLim);
                    
                //     std::cerr << "      Par " << p << " (" << fitFunc->GetParName(p) << "): " << val << " +/- " << err;
                    
                //     // Check if limits were applied
                //     if (minLim != maxLim) {
                //         std::cerr << "  [Limits: " << minLim << ", " << maxLim << "]";
                //         // Check if the value is pegged to the boundary (within a tiny tolerance)
                //         if (std::abs(val - minLim) < 1e-6 || std::abs(val - maxLim) < 1e-6) {
                //             std::cerr << "  <--- PEGGED TO LIMIT!";
                //         }
                //     }
                //     std::cerr << std::endl;
                // }

                // --- Salvage Logic ---
                // If the mean and sigma are within a physically reasonable window, keep them!
                double mu = fitFunc->GetParameter(1);
                double sig = fitFunc->GetParameter(2);
                
                // We can adjust these "sane" bounds later on, if needed
                if (mu > 1.105 && mu < 1.125 && sig > 0.0005 && sig < 0.008) {
                    // std::cerr << "      -> Salvaging fit: parameters are physically reasonable." << std::endl;
                    fitResultsPhi[iBin] = {mu, sig, true};
                } else {
                    // std::cerr << "      -> Discarding fit: parameters are too wild." << std::endl;
                    fitResultsPhi[iBin] = {0, 0, false};
                }
            }
        }

        // -----------------------------------------------------------------------------------------
        // 4.B: Cos Theta Bins (same code, different input histogram)
        // -----------------------------------------------------------------------------------------
        for (int iBin = 1; iBin <= h2dDeltaThetaVsMass->GetNbinsX(); ++iBin){
            // Project Counts
            TString projName = Form("hMass_DeltaTheta_Bin%d", iBin);
            TH1D* hMassProj = h2dDeltaThetaVsMass->ProjectionY(projName, iBin, iBin, "e");
            hMassProj->SetTitle(Form("Mass Projection cos(#theta) Bin %d;M_{p#pi} (GeV/c^{2});Counts", iBin));

            // Project Numerator and create <R>(m) QA plot
            TString numProjName = Form("hNum_DeltaTheta_Bin%d", iBin);
            TH1D* hNumProj = h2dRingObservableDeltaThetaVsMass->ProjectionY(numProjName, iBin, iBin, "e");
            
            TString qaName = Form("hQARingVsMass_DeltaTheta_Bin%d", iBin);
            TH1D* hQARing = (TH1D*)hNumProj->Clone(qaName);
            hQARing->SetTitle(Form("<R> vs Mass for cos(#theta) Bin %d;M_{p#pi} (GeV/c^{2});<R>", iBin));
            hQARing->Divide(hNumProj, hMassProj, 1.0, 1.0, "B"); 
            
            dirQARingMass->cd();
            hQARing->Write();

            // Fit Function
            TString fitName = Form("fit_DeltaTheta_Bin%d", iBin);
            TF1* fitFunc = new TF1(fitName, "gaus(0) + pol2(3)", massMin, massMax);
            
            double maxBinCenter = hMassProj->GetBinCenter(hMassProj->GetMaximumBin());
            double maxCount = hMassProj->GetMaximum();
            fitFunc->SetParameter(0, maxCount);
            fitFunc->SetParameter(1, lambdaPDGMassApprox);
            fitFunc->SetParameter(2, 0.002);

            fitFunc->SetParameter(0, maxCount);
            fitFunc->SetParLimits(1, lambdaPDGMassApprox-0.01, massMax+0.01); // Constraints close to mass (1.15 to 1.125)
            fitFunc->SetParLimits(2, 0.0001, 0.005);

            TFitResultPtr r = hMassProj->Fit(fitFunc, "Q 0 R S"); // Removed log-likelihood to make fits behave better in high-statistics 1D environment!

            dirFitsTheta->cd();
            hMassProj->GetListOfFunctions()->Add(fitFunc);
            hMassProj->Write();

            // Store results if fit is somewhat reasonable
            int fitStatus = int(r); // Extract the integer status code from the fit result
            
            // 0 = perfectly converged. 
            if (fitStatus == 0) {
                fitResultsTheta[iBin] = {fitFunc->GetParameter(1), fitFunc->GetParameter(2), true};
            }
            else {
                // --- Salvage Logic ---
                // If the mean and sigma are within a physically reasonable window, keep them!
                double muCheck = fitFunc->GetParameter(1);
                double sigmaCheck = fitFunc->GetParameter(2);
                
                // We can adjust these "sane" bounds later on, if needed
                if (muCheck > 1.105 && muCheck < 1.125 && sigmaCheck > 0.0005 && sigmaCheck < 0.008) {
                    fitResultsTheta[iBin] = {muCheck, sigmaCheck, true};
                } else {
                    fitResultsTheta[iBin] = {0, 0, false};
                }
            }
        }

        // =========================================================================================
        // Step 5: Discontinuous Sideband Background Fit using TGraphErrors (TH1D's are bad at that!)
        // =========================================================================================
        std::cout << "  [Step 5] Performing discontinuous sideband fits..." << std::endl;

        // Subdirectories to store the background graphs and their fits
        TDirectory* dirBkgFitsPhi = outDirVar->mkdir("BackgroundFits_DeltaPhi");
        TDirectory* dirBkgFitsTheta = outDirVar->mkdir("BackgroundFits_DeltaTheta");

        // -----------------------------------------------------------------------------------------
        // 5.A: Background Fits for Delta Phi Bins
        // -----------------------------------------------------------------------------------------
        for (int iBin = 1; iBin <= h2dDeltaPhiVsMass->GetNbinsX(); ++iBin){
            if (!fitResultsPhi[iBin].valid) continue; // Skip if step 4's fit failed

            double mu = fitResultsPhi[iBin].mu;
            double sigma = fitResultsPhi[iBin].sigma;

            // Re-fetch the mass projection (we just need the bin contents)
            TString projName = Form("hMass_DeltaPhi_Bin%d", iBin);
            TH1D* hMassProj = (TH1D*)dirFitsPhi->Get(projName);
            if (!hMassProj) continue;

            // Create a TGraphErrors to hold ONLY the sideband points
            TString graphName = Form("grBkg_DeltaPhi_Bin%d", iBin);
            TGraphErrors* grBkg = new TGraphErrors();
            grBkg->SetName(graphName);
            grBkg->SetTitle(Form("Sideband Bkg #Delta#phi Bin %d;M_{p#pi} (GeV/c^{2});Counts", iBin));

            int ptIdx = 0;
            // Loop over histogram bins and select only those in the sidebands
            for (int jBin = 1; jBin <= hMassProj->GetNbinsX(); ++jBin){
                double x = hMassProj->GetBinCenter(jBin);
                double y = hMassProj->GetBinContent(jBin);
                double ey = hMassProj->GetBinError(jBin);

                // Define intervals: [mu - 6*sigma, mu - 4*sigma] U [mu + 4*sigma, mu + 6*sigma]
                // The gap [mu - 4*sigma, mu - 3*sigma] and [mu + 3*sigma, mu + 4*sigma] decouples the fit
                // bool inLeftSideband = (x >= (mu - 6.0 * sigma) && x <= (mu - 4.0 * sigma));
                // bool inRightSideband = (x >= (mu + 4.0 * sigma) && x <= (mu + 6.0 * sigma));
                // UPDATE! Fixed to use larger bands. With 175 bins from 1.08 to 1.15 GeV/c^2, [4sigma,6sigma]
                // would mean to use only 10 bins on the right side! (same for other side!)
                // Now considering all bins in the left and right sides!
                // UPDATE 2! Now using a 6 sigma exclusion zone -- the peak was getting into the fit!
                    // This uses ALL available statistics outside the 6-sigma exclusion zone:
                bool inLeftSideband = (x >= massMin && x <= (mu - 6.0 * sigma));
                bool inRightSideband = (x >= (mu + 6.0 * sigma) && x <= massMax);

                if (inLeftSideband || inRightSideband){
                    grBkg->SetPoint(ptIdx, x, y);
                    // Standard x-error is 0 for fits, y-error is standard bin error
                    grBkg->SetPointError(ptIdx, 0.0, ey); 
                    ptIdx++;
                }
            }

            // Fit the discontinuous graph with a pol2
            // We use pol2 to cover the small curvature in Lambda background
            TString bkgFitName = Form("bkgFit_DeltaPhi_Bin%d", iBin);
            TF1* bkgFitFunc = new TF1(bkgFitName, "pol2", massMin, massMax);
            
            // "Q" = quiet, "0" = don't draw, "S" = return TFitResultPtr
            TFitResultPtr rBkg = grBkg->Fit(bkgFitFunc, "Q 0 S");

            dirBkgFitsPhi->cd();
            grBkg->GetListOfFunctions()->Add(bkgFitFunc); // Attach fit for viewing in TBrowser
            grBkg->Write();
        }

        // -----------------------------------------------------------------------------------------
        // 5.B: Background Fits for Cos Theta Bins
        // -----------------------------------------------------------------------------------------
        for (int iBin = 1; iBin <= h2dDeltaThetaVsMass->GetNbinsX(); ++iBin){
            if (!fitResultsTheta[iBin].valid) continue;

            double mu = fitResultsTheta[iBin].mu;
            double sigma = fitResultsTheta[iBin].sigma;

            TString projName = Form("hMass_DeltaTheta_Bin%d", iBin);
            TH1D* hMassProj = (TH1D*)dirFitsTheta->Get(projName);
            if (!hMassProj) continue;

            TString graphName = Form("grBkg_DeltaTheta_Bin%d", iBin);
            TGraphErrors* grBkg = new TGraphErrors();
            grBkg->SetName(graphName);
            grBkg->SetTitle(Form("Sideband Bkg cos(#theta) Bin %d;M_{p#pi} (GeV/c^{2});Counts", iBin));

            int ptIdx = 0;
            for (int jBin = 1; jBin <= hMassProj->GetNbinsX(); ++jBin){
                double x = hMassProj->GetBinCenter(jBin);
                double y = hMassProj->GetBinContent(jBin);
                double ey = hMassProj->GetBinError(jBin);

                bool inLeftSideband = (x >= massMin && x <= (mu - 6.0 * sigma));
                bool inRightSideband = (x >= (mu + 6.0 * sigma) && x <= massMax);

                if (inLeftSideband || inRightSideband){
                    grBkg->SetPoint(ptIdx, x, y);
                    grBkg->SetPointError(ptIdx, 0.0, ey);
                    ptIdx++;
                }
            }

            TString bkgFitName = Form("bkgFit_DeltaTheta_Bin%d", iBin);
            TF1* bkgFitFunc = new TF1(bkgFitName, "pol2", massMin, massMax);
            
            TFitResultPtr rBkg = grBkg->Fit(bkgFitFunc, "Q 0 S");

            dirBkgFitsTheta->cd();
            grBkg->GetListOfFunctions()->Add(bkgFitFunc);
            grBkg->Write();
        }

        // =========================================================================================
        // Step 6: Signal Yield and Polarization Extraction
        // =========================================================================================
        std::cout << "  [Step 6] Extracting Signal Yields and <R>_signal..." << std::endl;

        TDirectory* dirResultsPhi = outDirVar->mkdir("Results_DeltaPhi");
        TDirectory* dirResultsTheta = outDirVar->mkdir("Results_DeltaTheta");
        
        // NEW: Directories for Numerator Fit QA
        TDirectory* dirBkgNumFitsPhi = outDirVar->mkdir("BkgNumFits_DeltaPhi");
        TDirectory* dirBkgNumFitsTheta = outDirVar->mkdir("BkgNumFits_DeltaTheta");

        // Prepare output histograms (Cloning ProjectionX to keep exact angular binning)
        TH1D* hSigYieldPhi = (TH1D*)h2dDeltaPhiVsMass->ProjectionX("hSigYieldPhi"); hSigYieldPhi->Reset();
        hSigYieldPhi->SetTitle("Signal Yield vs #Delta#phi;#Delta#phi;Counts");
        
        TH1D* hBkgYieldPhi = (TH1D*)hSigYieldPhi->Clone("hBkgYieldPhi");
        hBkgYieldPhi->SetTitle("Background Yield vs #Delta#phi;#Delta#phi;Counts");
        
        TH1D* hRSigPhi = (TH1D*)hSigYieldPhi->Clone("hRSigPhi");
        hRSigPhi->SetTitle("Extracted <R>_{signal} vs #Delta#phi;#Delta#phi;<R>_{S}");
        
        TH1D* hRBkgPhi = (TH1D*)hSigYieldPhi->Clone("hRBkgPhi");
        hRBkgPhi->SetTitle("Background <R>_{bkg} vs #Delta#phi;#Delta#phi;<R>_{B}");

        // --- NEW QA HISTOGRAMS: Delta Phi statistical significance and purity ---
        TH1D* hPurityPhi = (TH1D*)hSigYieldPhi->Clone("hPurityPhi");
        hPurityPhi->SetTitle("Purity S/(S+B) vs #Delta#phi;#Delta#phi;Purity");
        
        TH1D* hSigStatPhi = (TH1D*)hSigYieldPhi->Clone("hSigStatPhi");
        hSigStatPhi->SetTitle("Significance S/#sqrt{S+B} vs #Delta#phi;#Delta#phi;Significance");

        // -----------------------------------------------------------------------------------------
        // 6.A: Extraction for Delta Phi Bins
        // -----------------------------------------------------------------------------------------
        for (int iBin = 1; iBin <= h2dDeltaPhiVsMass->GetNbinsX(); ++iBin){
            if (!fitResultsPhi[iBin].valid) continue;

            double mu = fitResultsPhi[iBin].mu;
            double sigma = fitResultsPhi[iBin].sigma;

            // Fetch denominator (counts) and numerator (observable) mass projections
            TString projName = Form("hMass_DeltaPhi_Bin%d", iBin);
            TH1D* hMassProj = (TH1D*)dirFitsPhi->Get(projName);
            
            TString numProjName = Form("hNum_DeltaPhi_Bin%d", iBin);
            TH1D* hNumProj = h2dRingObservableDeltaPhiVsMass->ProjectionY(numProjName, iBin, iBin, "e");
            
            // Also fetch the squared observable projection to get the variance of R_peak
            TString sqProjName = Form("hSqNum_DeltaPhi_Bin%d", iBin);
            TH1D* hSqNumProj = h2dRingObservableSquaredDeltaPhiVsMass->ProjectionY(sqProjName, iBin, iBin, "e");

            // Fetch the background counts fit and its result pointer
            TString bkgFitName = Form("bkgFit_DeltaPhi_Bin%d", iBin);
            TGraphErrors* grBkg = (TGraphErrors*)dirBkgFitsPhi->Get(Form("grBkg_DeltaPhi_Bin%d", iBin));
            if (!grBkg) continue;
            TF1* bkgFitFunc = grBkg->GetFunction(bkgFitName);
            TFitResultPtr rBkg = grBkg->Fit(bkgFitFunc, "Q 0 S"); // Re-fit quickly to grab the TFitResultPtr covariance

            // --- Fit the Numerator Sidebands to get Background Polarization ---
            TGraphErrors* grBkgNum = new TGraphErrors();
            int ptIdx = 0;
            for (int jBin = 1; jBin <= hNumProj->GetNbinsX(); ++jBin){
                double x = hNumProj->GetBinCenter(jBin);
                bool inLeft = (x >= massMin && x <= (mu - 6.0 * sigma)); // Updated to use larger sidebands
                bool inRight = (x >= (mu + 6.0 * sigma) && x <= massMax);
                if (inLeft || inRight){
                    grBkgNum->SetPoint(ptIdx, x, hNumProj->GetBinContent(jBin));
                    grBkgNum->SetPointError(ptIdx, 0.0, hNumProj->GetBinError(jBin));
                    ptIdx++;
                }
            }
            TF1* bkgNumFitFunc = new TF1(Form("bkgNumFit_Phi_%d", iBin), "pol2", massMin, massMax);
            TFitResultPtr rBkgNum = grBkgNum->Fit(bkgNumFitFunc, "Q 0 S");

            // --- 6.1 & 6.3: Define Signal Region and Exact Limits ---
            int firstBin = hMassProj->FindBin(mu - 3.0 * sigma);
            int lastBin  = hMassProj->FindBin(mu + 3.0 * sigma);
            double x_low  = hMassProj->GetBinLowEdge(firstBin);
            double x_high = hMassProj->GetBinLowEdge(lastBin) + hMassProj->GetBinWidth(lastBin);
            double binWidth = hMassProj->GetBinWidth(1); // Assuming uniform binning

            // --- 6.2: Raw Signal Counts (Peak Region) ---
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
            
            // --- 6.4: Background Estimation in ONE STEP ---
            double bkgCounts = bkgFitFunc->Integral(x_low, x_high) / binWidth;
            double errBkgCounts = bkgFitFunc->IntegralError(x_low, x_high, rBkg->GetParams(), rBkg->GetCovarianceMatrix().GetMatrixArray()) / binWidth;

            double bkgNum = bkgNumFitFunc->Integral(x_low, x_high) / binWidth;
            double errBkgNum = bkgNumFitFunc->IntegralError(x_low, x_high, rBkgNum->GetParams(), rBkgNum->GetCovarianceMatrix().GetMatrixArray()) / binWidth;

            // --- 6.5: Final Signal Yield ---
            double sigCounts = totCounts - bkgCounts;
            double errSigCounts = std::sqrt(totCountsErrSq + errBkgCounts * errBkgCounts);

            if (sigCounts <= 0 || totCounts <= 0) continue; // Safety check

            // Fractions
            double fB = bkgCounts / totCounts;
            double fS = sigCounts / totCounts; // Equivalent to 1 - fB

            // --- Calculate Observables ---
            double R_peak = totNum / totCounts;
            double R_B = bkgNum / bkgCounts;
            
            // True Signal Polarization: R_S = (R_peak - f_B * R_B) / f_S
            double R_S = (R_peak - fB * R_B) / fS;

            // --- 6.6: Error Propagation for Polarization ---
            // Variance of R_peak (using squared observable data)
            double varR_peak = (totSqNum / totCounts) - (R_peak * R_peak);
            varR_peak = varR_peak / totCounts; // Standard error of the mean

            // Variance of R_B (simple error propagation for ratio B_num / B_counts)
            double varR_B = std::pow(R_B, 2) * ( std::pow(errBkgNum/bkgNum, 2) + std::pow(errBkgCounts/bkgCounts, 2) );

            // Variance of f_B
            double var_fB = std::pow(fB, 2) * ( std::pow(errBkgCounts/bkgCounts, 2) + std::pow(errTotCounts/totCounts, 2) );

            // Final Variance Formula:
            // sigma^2(R_S) = (1/f_S^2)*sigma^2(R_peak) + (f_B^2/f_S^2)*sigma^2(R_B) + ((R_peak - R_B)/f_S^2)^2 * sigma^2(f_B)
            double term1 = varR_peak / std::pow(fS, 2);
            double term2 = (std::pow(fB, 2) / std::pow(fS, 2)) * varR_B;
            double term3 = std::pow((R_peak - R_B) / std::pow(fS, 2), 2) * var_fB;
            
            double errR_S = std::sqrt(term1 + term2 + term3);

            // Fill Histograms
            hSigYieldPhi->SetBinContent(iBin, sigCounts);
            hSigYieldPhi->SetBinError(iBin, errSigCounts);
            
            hBkgYieldPhi->SetBinContent(iBin, bkgCounts);
            hBkgYieldPhi->SetBinError(iBin, errBkgCounts);

            hRSigPhi->SetBinContent(iBin, R_S);
            hRSigPhi->SetBinError(iBin, errR_S);

            hRBkgPhi->SetBinContent(iBin, R_B);
            hRBkgPhi->SetBinError(iBin, std::sqrt(varR_B));

            // --- Calculate Purity & Significance ---
            double purity = fS; 
            double errPurity = std::sqrt(var_fB); // Since fS = 1 - fB, the variance is identical

            double significance = sigCounts / std::sqrt(totCounts);
            // Standard error propagation for Z = S/sqrt(N)
            double errSignificance = std::sqrt( (errSigCounts*errSigCounts)/totCounts + 
                                                (sigCounts*sigCounts*totCountsErrSq)/(4.0*std::pow(totCounts, 3)) );

            hPurityPhi->SetBinContent(iBin, purity);
            hPurityPhi->SetBinError(iBin, errPurity);

            hSigStatPhi->SetBinContent(iBin, significance);
            hSigStatPhi->SetBinError(iBin, errSignificance);

            // Writing bkgNum histograms in a specific folder for QA:
            dirBkgNumFitsPhi->cd();
            
            // Give the graph a unique name so it doesn't overwrite
            grBkgNum->SetName(Form("grBkgNum_DeltaPhi_Bin%d", iBin));
            grBkgNum->SetTitle(Form("Numerator Sidebands DeltaPhi Bin %d;Invariant Mass (GeV/c^{2});#Sigma r_{i}", iBin));
            
            // Attach the fit function to the graph for easy viewing in TBrowser
            grBkgNum->GetListOfFunctions()->Add(bkgNumFitFunc);
            
            // Write both the sideband graph (with fit) and the raw projected histogram
            grBkgNum->Write();
            hNumProj->Write(); 
            
            delete grBkgNum;
            // Note: We do NOT delete bkgNumFitFunc here anymore because ROOT takes ownership 
            // when we add it to the graph's list of functions. Deleting it will cause a crash!
            // delete bkgNumFitFunc;
        }
        
        dirResultsPhi->cd();
            // Resetting stats boxes to get an estimate of Signal vs Background Lambda counts:
        hSigYieldPhi->ResetStats();
        hBkgYieldPhi->ResetStats();
        hSigYieldPhi->Write();
        hBkgYieldPhi->Write();
        hPurityPhi->Write();
        hSigStatPhi->Write();
        hRSigPhi->Write();
        hRBkgPhi->Write();

        // -----------------------------------------------------------------------------------------
        // 6.B: Extraction for Cos Theta Bins
        // -----------------------------------------------------------------------------------------
        
        // Prepare output histograms
        TH1D* hSigYieldTheta = (TH1D*)h2dDeltaThetaVsMass->ProjectionX("hSigYieldTheta"); hSigYieldTheta->Reset();
        hSigYieldTheta->SetTitle("Signal Yield vs cos(#theta);cos(#theta);Counts");
        
        TH1D* hBkgYieldTheta = (TH1D*)hSigYieldTheta->Clone("hBkgYieldTheta");
        hBkgYieldTheta->SetTitle("Background Yield vs cos(#theta);cos(#theta);Counts");
        
        TH1D* hRSigTheta = (TH1D*)hSigYieldTheta->Clone("hRSigTheta");
        hRSigTheta->SetTitle("Extracted <R>_{signal} vs cos(#theta);cos(#theta);<R>_{S}");
        
        TH1D* hRBkgTheta = (TH1D*)hSigYieldTheta->Clone("hRBkgTheta");
        hRBkgTheta->SetTitle("Background <R>_{bkg} vs cos(#theta);cos(#theta);<R>_{B}");

        // --- NEW QA HISTOGRAMS: Cos Theta ---
        TH1D* hPurityTheta = (TH1D*)hSigYieldTheta->Clone("hPurityTheta");
        hPurityTheta->SetTitle("Purity S/(S+B) vs cos(#theta);cos(#theta);Purity");
        
        TH1D* hSigStatTheta = (TH1D*)hSigYieldTheta->Clone("hSigStatTheta");
        hSigStatTheta->SetTitle("Significance S/#sqrt{S+B} vs cos(#theta);cos(#theta);Significance");

        for (int iBin = 1; iBin <= h2dDeltaThetaVsMass->GetNbinsX(); ++iBin){
            if (!fitResultsTheta[iBin].valid) continue;

            double mu = fitResultsTheta[iBin].mu;
            double sigma = fitResultsTheta[iBin].sigma;

            // Fetch denominator (counts) and numerator (observable) mass projections
            TString projName = Form("hMass_DeltaTheta_Bin%d", iBin);
            TH1D* hMassProj = (TH1D*)dirFitsTheta->Get(projName);
            
            TString numProjName = Form("hNum_DeltaTheta_Bin%d", iBin);
            TH1D* hNumProj = h2dRingObservableDeltaThetaVsMass->ProjectionY(numProjName, iBin, iBin, "e");
            
            // Fetch squared observable projection
            TString sqProjName = Form("hSqNum_DeltaTheta_Bin%d", iBin);
            TH1D* hSqNumProj = h2dRingObservableSquaredDeltaThetaVsMass->ProjectionY(sqProjName, iBin, iBin, "e");

            // Fetch the background counts fit and its result pointer
            TString bkgFitName = Form("bkgFit_DeltaTheta_Bin%d", iBin);
            TGraphErrors* grBkg = (TGraphErrors*)dirBkgFitsTheta->Get(Form("grBkg_DeltaTheta_Bin%d", iBin));
            if (!grBkg) continue;
            TF1* bkgFitFunc = grBkg->GetFunction(bkgFitName);
            TFitResultPtr rBkg = grBkg->Fit(bkgFitFunc, "Q 0 S"); 

            // --- Fit the Numerator Sidebands to get Background Polarization ---
            TGraphErrors* grBkgNum = new TGraphErrors();
            int ptIdx = 0;
            for (int jBin = 1; jBin <= hNumProj->GetNbinsX(); ++jBin){
                double x = hNumProj->GetBinCenter(jBin);
                bool inLeft = (x >= massMin && x <= (mu - 6.0 * sigma)); // Updated to use larger sidebands
                bool inRight = (x >= (mu + 6.0 * sigma) && x <= massMax);
                if (inLeft || inRight){
                    grBkgNum->SetPoint(ptIdx, x, hNumProj->GetBinContent(jBin));
                    grBkgNum->SetPointError(ptIdx, 0.0, hNumProj->GetBinError(jBin));
                    ptIdx++;
                }
            }
            TF1* bkgNumFitFunc = new TF1(Form("bkgNumFit_Theta_%d", iBin), "pol2", massMin, massMax);
            TFitResultPtr rBkgNum = grBkgNum->Fit(bkgNumFitFunc, "Q 0 S");

            // --- Define Signal Region and Exact Limits ---
            int firstBin = hMassProj->FindBin(mu - 3.0 * sigma);
            int lastBin  = hMassProj->FindBin(mu + 3.0 * sigma);
            double x_low  = hMassProj->GetBinLowEdge(firstBin);
            double x_high = hMassProj->GetBinLowEdge(lastBin) + hMassProj->GetBinWidth(lastBin);
            double binWidth = hMassProj->GetBinWidth(1); 

            // --- Raw Signal Counts (Peak Region) ---
            double totCounts = 0;
            double totCountsErrSq = 0;
            double totNum = 0;
            double totNumErrSq = 0;
            double totSqNum = 0; 

            for (int jBin = firstBin; jBin <= lastBin; ++jBin){
                totCounts += hMassProj->GetBinContent(jBin);
                totCountsErrSq += std::pow(hMassProj->GetBinError(jBin), 2);
                
                totNum += hNumProj->GetBinContent(jBin);
                totNumErrSq += std::pow(hNumProj->GetBinError(jBin), 2);
                
                totSqNum += hSqNumProj->GetBinContent(jBin);
            }
            double errTotCounts = std::sqrt(totCountsErrSq);
            
            // --- Background Estimation in ONE STEP (Covariance Integration) ---
            double bkgCounts = bkgFitFunc->Integral(x_low, x_high) / binWidth;
            double errBkgCounts = bkgFitFunc->IntegralError(x_low, x_high, rBkg->GetParams(), rBkg->GetCovarianceMatrix().GetMatrixArray()) / binWidth;

            double bkgNum = bkgNumFitFunc->Integral(x_low, x_high) / binWidth;
            double errBkgNum = bkgNumFitFunc->IntegralError(x_low, x_high, rBkgNum->GetParams(), rBkgNum->GetCovarianceMatrix().GetMatrixArray()) / binWidth;

            // --- Final Signal Yield ---
            double sigCounts = totCounts - bkgCounts;
            double errSigCounts = std::sqrt(totCountsErrSq + errBkgCounts * errBkgCounts);

            if (sigCounts <= 0 || totCounts <= 0) continue; 

            // Fractions
            double fB = bkgCounts / totCounts;
            double fS = sigCounts / totCounts; 

            // --- Calculate Observables ---
            double R_peak = totNum / totCounts;
            double R_B = bkgNum / bkgCounts;
            double R_S = (R_peak - fB * R_B) / fS;

            // --- Error Propagation for Polarization ---
            double varR_peak = (totSqNum / totCounts) - (R_peak * R_peak);
            varR_peak = varR_peak / totCounts; 

            double varR_B = std::pow(R_B, 2) * ( std::pow(errBkgNum/bkgNum, 2) + std::pow(errBkgCounts/bkgCounts, 2) );
            double var_fB = std::pow(fB, 2) * ( std::pow(errBkgCounts/bkgCounts, 2) + std::pow(errTotCounts/totCounts, 2) );

            double term1 = varR_peak / std::pow(fS, 2);
            double term2 = (std::pow(fB, 2) / std::pow(fS, 2)) * varR_B;
            double term3 = std::pow((R_peak - R_B) / std::pow(fS, 2), 2) * var_fB;
            
            double errR_S = std::sqrt(term1 + term2 + term3);

            // Fill Histograms
            hSigYieldTheta->SetBinContent(iBin, sigCounts);
            hSigYieldTheta->SetBinError(iBin, errSigCounts);
            
            hBkgYieldTheta->SetBinContent(iBin, bkgCounts);
            hBkgYieldTheta->SetBinError(iBin, errBkgCounts);

            hRSigTheta->SetBinContent(iBin, R_S);
            hRSigTheta->SetBinError(iBin, errR_S);

            hRBkgTheta->SetBinContent(iBin, R_B);
            hRBkgTheta->SetBinError(iBin, std::sqrt(varR_B));

            // --- Calculate Purity & Significance ---
            double purity = fS; 
            double errPurity = std::sqrt(var_fB); // Since fS = 1 - fB, the variance is identical

            double significance = sigCounts / std::sqrt(totCounts);
            // Standard error propagation for Z = S/sqrt(N)
            double errSignificance = std::sqrt( (errSigCounts*errSigCounts)/totCounts + 
                                                (sigCounts*sigCounts*totCountsErrSq)/(4.0*std::pow(totCounts, 3)) );

            hPurityTheta->SetBinContent(iBin, purity);
            hPurityTheta->SetBinError(iBin, errPurity);
            
            hSigStatTheta->SetBinContent(iBin, significance);
            hSigStatTheta->SetBinError(iBin, errSignificance);
            
            // Writing bkgNum histograms in a specific folder for QA:
            dirBkgNumFitsTheta->cd();
            
            // Give the graph a unique name so it doesn't overwrite
            grBkgNum->SetName(Form("grBkgNum_DeltaTheta_Bin%d", iBin));
            grBkgNum->SetTitle(Form("Numerator Sidebands DeltaTheta Bin %d;Invariant Mass (GeV/c^{2});#Sigma r_{i}", iBin));
            
            // Attach the fit function to the graph for easy viewing in TBrowser
            grBkgNum->GetListOfFunctions()->Add(bkgNumFitFunc);
            
            // Write both the sideband graph (with fit) and the raw projected histogram
            grBkgNum->Write();
            hNumProj->Write(); 
            
            delete grBkgNum;
        }
        
        dirResultsTheta->cd();
            // Resetting stats boxes to get an estimate of Signal vs Background Lambda counts:
        hSigYieldTheta->ResetStats();
        hBkgYieldTheta->ResetStats();
        hSigYieldTheta->Write();
        hBkgYieldTheta->Write();
        hSigStatTheta->Write();
        hPurityTheta->Write();
        hRSigTheta->Write();
        hRBkgTheta->Write();

        // =========================================================================================
        // Step 7: Save Input Histograms for Traceability
        // =========================================================================================
        std::cout << "  [Step 7] Saving all input histograms for " << var << "..." << std::endl;

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
    outFile->Close();
    inFile->Close();
    
    std::cout << "\n========================================================================" << std::endl;
    std::cout << " Signal Extraction Complete. File saved as: " << outFileName << std::endl;
    std::cout << "========================================================================\n" << std::endl;
    
}