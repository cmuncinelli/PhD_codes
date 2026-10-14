/**************************************************************************************************
 *
 *  signalExtractionRing.cxx
 *
 *  Description:
 *  ----------------------------------------------------------------------------------------------
 *  Signal extraction framework for the Ring Observable polarization analysis in Λ–jet systems.
 *
 *  This code performs invariant mass fits in bins of angular and/or kinematic observables in order
 *  to extract the true signal polarization ⟨R⟩_signal from a mixture of signal and combinatorial
 *  background Λ candidates.
 *
 *  The necessity of signal extraction arises from the fact that the measured polarization inside
 *  the Λ invariant mass peak is not purely signal:
 *
 *      ⟨R⟩_measured = f_S ⟨R⟩_signal + f_B ⟨R⟩_background
 *
 *  where background polarization is NOT guaranteed to be zero.
 *
 *  In heavy-ion collisions, local angular momentum can generate polarization.
 *  However, combinatorial background may:
 *      • Dilute the signal
 *      • Be polarized differently
 *      • Have distorted angular structure
 *
 *  Therefore, proper signal extraction is mandatory.
 *
 *
 *  IMPORTANT NON-TRIVIAL ASSUMPTIONS:
 *  ----------------------------------------------------------------------------------------------
 *  1) Background polarization is assumed smooth in invariant mass.
 *  2) No strong mass–observable correlation exists.
 *  3) Sidebands represent the background under the peak.
 *
 *  These assumptions MUST be validated.
 *
 *  If they fail, a simultaneous mass–polarization fit must be implemented:
 *
 *      Numerator(m)   = S(m) R_S + B(m) R_B
 *      Denominator(m) = S(m) + B(m)
 *
 *  which avoids sideband assumptions entirely.
 *
 *
 *  Authors:
 *      Cicero (Analysis Lead)
 *      ChatGPT (AI-assisted implementation, 2026)
 *
 *  Compilation:
 *      g++ signalExtractionRing.cxx `root-config --cflags --libs` -O2 -o signalExtractionRing.exe
 *  Execution:
 *      ./signalExtractionRing /path/to/ConsumerResults_XXXX.root
 **************************************************************************************************/

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cmath>

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TMatrixDSym.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TCanvas.h"
#include "TMath.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Helper function: Extract suffix from input file name
////////////////////////////////////////////////////////////////////////////////////////////////////
std::string extractSuffix(const std::string& inputPath)
{
    size_t pos = inputPath.find("ConsumerResults_");
    if (pos == std::string::npos)
        return "Unknown";

    std::string suffix = inputPath.substr(pos + std::string("ConsumerResults_").length());

    size_t dot = suffix.find(".root");
    if (dot != std::string::npos)
        suffix = suffix.substr(0, dot);

    return suffix;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Main Analysis Function for One Observable Variation
////////////////////////////////////////////////////////////////////////////////////////////////////
void processVariation(TDirectory* inputDir,
                      TDirectory* outputDir,
                      const std::string& variationName)
{
    std::cout << "\n============================================================";
    std::cout << "\nProcessing variation: " << variationName;
    std::cout << "\n============================================================\n";

    // Create subdirectory in output
    TDirectory* varOutDir = outputDir->mkdir(variationName.c_str());
    varOutDir->cd();

    // ---------------------------------------------------------------------
    // Load required histograms
    // ---------------------------------------------------------------------

    std::cout << "Loading histograms...\n";

    TH2D* h2dDeltaPhiVsMass =
        dynamic_cast<TH2D*>(inputDir->Get(
            (variationName + "/h2dDeltaPhiVsMass").c_str()));

    TH2D* h2dRingObservableDeltaPhiVsMass =
        dynamic_cast<TH2D*>(inputDir->Get(
            (variationName + "/h2dRingObservableDeltaPhiVsMass").c_str()));

    TH2D* h2dRingObservableSquaredDeltaPhiVsMass =
        dynamic_cast<TH2D*>(inputDir->Get(
            (variationName + "/h2dRingObservableSquaredDeltaPhiVsMass").c_str()));

    if (!h2dDeltaPhiVsMass || !h2dRingObservableDeltaPhiVsMass)
    {
        std::cerr << "ERROR: Required histograms not found for "
                  << variationName << std::endl;
        return;
    }

    std::cout << "Histograms loaded successfully.\n";

    // ---------------------------------------------------------------------
    // Prepare output subfolders
    // ---------------------------------------------------------------------

    TDirectory* fitsDir = varOutDir->mkdir("MassFits");
    fitsDir->cd();

    // ---------------------------------------------------------------------
    // Get axis information
    // ---------------------------------------------------------------------

    int nPhiBins = h2dDeltaPhiVsMass->GetXaxis()->GetNbins();
    double massMin = h2dDeltaPhiVsMass->GetYaxis()->GetXmin();
    double massMax = h2dDeltaPhiVsMass->GetYaxis()->GetXmax();
    double massRange = massMax - massMin;

    std::cout << "Number of DeltaPhi bins: " << nPhiBins << std::endl;
    std::cout << "Mass range: [" << massMin << ", " << massMax << "]\n";

    // Containers to store extracted mean and sigma per bin
    std::vector<double> meanPhi(nPhiBins+1, 0.0);
    std::vector<double> sigmaPhi(nPhiBins+1, 0.0);

    // ---------------------------------------------------------------------
    // Loop over DeltaPhi bins
    // ---------------------------------------------------------------------

    for (int iPhi = 1; iPhi <= nPhiBins; ++iPhi)
    {
        std::cout << "Processing DeltaPhi bin " << iPhi
                  << " / " << nPhiBins << std::endl;

        // Restrict to one DeltaPhi bin
        h2dDeltaPhiVsMass->GetXaxis()->SetRange(iPhi, iPhi);

        // Project mass (Y axis)
        TH1D* hMassProj =
            dynamic_cast<TH1D*>(h2dDeltaPhiVsMass->ProjectionY(
                Form("hMassProj_bin%d", iPhi),
                iPhi, iPhi,
                "e"));  // propagate errors

        if (!hMassProj || hMassProj->GetEntries() < 10)
        {
            std::cout << "Skipping bin (low statistics)\n";
            continue;
        }

        // -----------------------------------------------------------------
        // Define Gaussian + quadratic function
        // -----------------------------------------------------------------

        TF1* fitFunc = new TF1(
            Form("fit_bin%d", iPhi),
            "[0]*TMath::Gaus(x,[1],[2],true) + [3] + [4]*x + [5]*x*x",
            massMin, massMax);

        // Initial parameter guesses
        double peakHeight = hMassProj->GetMaximum();
        double peakPos    = hMassProj->GetBinCenter(hMassProj->GetMaximumBin());
        double sigmaGuess = 0.005;

        fitFunc->SetParameters(
            peakHeight,  // Gaussian amplitude
            peakPos,     // mean
            sigmaGuess,  // sigma
            1.0, 0.0, 0.0);  // polynomial params

        // -----------------------------------------------------------------
        // Stability constraints
        // -----------------------------------------------------------------

        fitFunc->SetParLimits(1, massMin, massMax);       // mean bounds
        fitFunc->SetParLimits(2, 1e-4, massRange);        // sigma bounds

        // -----------------------------------------------------------------
        // Perform Poisson likelihood fit
        // -----------------------------------------------------------------

        TFitResultPtr fitResult =
            hMassProj->Fit(fitFunc, "QLS");

        if (int(fitResult) != 0)
        {
            std::cout << "Fit failed in bin " << iPhi << "\n";
            delete hMassProj;
            delete fitFunc;
            continue;
        }

        double mean  = fitFunc->GetParameter(1);
        double sigma = fitFunc->GetParameter(2);

        std::cout << "  Extracted mean  = " << mean  << "\n";
        std::cout << "  Extracted sigma = " << sigma << "\n";

        meanPhi[iPhi]  = mean;
        sigmaPhi[iPhi] = sigma;

        // Save projection + fit
        fitsDir->cd();
        hMassProj->Write();
        fitFunc->Write();

        delete hMassProj;
        delete fitFunc;
    }

    // Reset axis range
    h2dDeltaPhiVsMass->GetXaxis()->SetRange(0, -1);

    // ============================================================
    // 2D CASE: DeltaTheta vs Mass
    // ============================================================

    std::cout << "\nStarting DeltaTheta fits...\n";

    TH2D* h2dDeltaThetaVsMass =
        dynamic_cast<TH2D*>(inputDir->Get(
            (variationName + "/h2dDeltaThetaVsMass").c_str()));

    if (!h2dDeltaThetaVsMass)
    {
        std::cerr << "WARNING: Missing h2dDeltaThetaVsMass for "
                  << variationName << std::endl;
    }
    else
    {
        TDirectory* fitsThetaDir = varOutDir->mkdir("MassFits_DeltaTheta");
        fitsThetaDir->cd();

        int nThetaBins = h2dDeltaThetaVsMass->GetXaxis()->GetNbins();

        std::vector<double> meanTheta(nThetaBins+1, 0.0);
        std::vector<double> sigmaTheta(nThetaBins+1, 0.0);

        for (int iTheta = 1; iTheta <= nThetaBins; ++iTheta)
        {
            std::cout << "Processing DeltaTheta bin "
                      << iTheta << " / " << nThetaBins << std::endl;

            h2dDeltaThetaVsMass->GetXaxis()->SetRange(iTheta, iTheta);

            TH1D* hMassProj =
                dynamic_cast<TH1D*>(h2dDeltaThetaVsMass->ProjectionY(
                    Form("hMassProj_theta_bin%d", iTheta),
                    iTheta, iTheta,
                    "e"));

            if (!hMassProj || hMassProj->GetEntries() < 10)
            {
                std::cout << "Skipping bin (low statistics)\n";
                continue;
            }

            TF1* fitFunc = new TF1(
                Form("fit_theta_bin%d", iTheta),
                "[0]*TMath::Gaus(x,[1],[2],true) + [3] + [4]*x + [5]*x*x",
                massMin, massMax);

            double peakHeight = hMassProj->GetMaximum();
            double peakPos    = hMassProj->GetBinCenter(hMassProj->GetMaximumBin());
            double sigmaGuess = 0.005;

            fitFunc->SetParameters(
                peakHeight, peakPos, sigmaGuess,
                1.0, 0.0, 0.0);

            fitFunc->SetParLimits(1, massMin, massMax);
            fitFunc->SetParLimits(2, 1e-4, massRange);

            TFitResultPtr fitResult =
                hMassProj->Fit(fitFunc, "QLS");

            if (int(fitResult) != 0)
            {
                std::cout << "Fit failed in bin " << iTheta << "\n";
                delete hMassProj;
                delete fitFunc;
                continue;
            }

            meanTheta[iTheta]  = fitFunc->GetParameter(1);
            sigmaTheta[iTheta] = fitFunc->GetParameter(2);

            fitsThetaDir->cd();
            hMassProj->Write();
            fitFunc->Write();

            delete hMassProj;
            delete fitFunc;
        }

        h2dDeltaThetaVsMass->GetXaxis()->SetRange(0, -1);
    }

    // ============================================================
    // 3D CASE: DeltaPhi vs Mass vs LambdaPt
    // ============================================================

    std::cout << "\nStarting 3D DeltaPhi–Mass–LambdaPt fits...\n";

    TH3D* h3dDeltaPhiVsMassVsLambdaPt =
        dynamic_cast<TH3D*>(inputDir->Get(
            (variationName + "/h3dDeltaPhiVsMassVsLambdaPt").c_str()));

    if (h3dDeltaPhiVsMassVsLambdaPt)
    {
        TDirectory* fits3DPhiDir =
            varOutDir->mkdir("MassFits_3D_DeltaPhi_LambdaPt");

        int nPtBins =
            h3dDeltaPhiVsMassVsLambdaPt->GetZaxis()->GetNbins();

        for (int iPt = 1; iPt <= nPtBins; ++iPt)
        {
            std::cout << "Processing pT bin "
                      << iPt << " / " << nPtBins << std::endl;

            // Restrict pT
            h3dDeltaPhiVsMassVsLambdaPt->GetZaxis()->SetRange(iPt, iPt);

            // Project to 2D: DeltaPhi vs Mass
            TH2D* h2dProj =
                dynamic_cast<TH2D*>(
                    h3dDeltaPhiVsMassVsLambdaPt->Project3D("yx"));

            if (!h2dProj)
                continue;

            int nPhiBins3D = h2dProj->GetXaxis()->GetNbins();

            for (int iPhi = 1; iPhi <= nPhiBins3D; ++iPhi)
            {
                h2dProj->GetXaxis()->SetRange(iPhi, iPhi);

                TH1D* hMassProj =
                    dynamic_cast<TH1D*>(h2dProj->ProjectionY(
                        Form("hMassProj_3D_phi%d_pt%d", iPhi, iPt),
                        iPhi, iPhi,
                        "e"));

                if (!hMassProj || hMassProj->GetEntries() < 10)
                    continue;

                TF1* fitFunc = new TF1(
                    Form("fit_3D_phi%d_pt%d", iPhi, iPt),
                    "[0]*TMath::Gaus(x,[1],[2],true) + [3] + [4]*x + [5]*x*x",
                    massMin, massMax);

                fitFunc->SetParameters(
                    hMassProj->GetMaximum(),
                    hMassProj->GetBinCenter(hMassProj->GetMaximumBin()),
                    0.005,
                    1.0, 0.0, 0.0);

                fitFunc->SetParLimits(1, massMin, massMax);
                fitFunc->SetParLimits(2, 1e-4, massRange);

                hMassProj->Fit(fitFunc, "QLS");

                fits3DPhiDir->cd();
                hMassProj->Write();
                fitFunc->Write();

                delete hMassProj;
                delete fitFunc;
            }

            delete h2dProj;
        }

        h3dDeltaPhiVsMassVsLambdaPt->GetZaxis()->SetRange(0, -1);
    }

        // ============================================================
    // 3D CASE: DeltaTheta vs Mass vs LambdaPt
    // ============================================================

    std::cout << "\nStarting 3D DeltaTheta–Mass–LambdaPt fits...\n";

    TH3D* h3dDeltaThetaVsMassVsLambdaPt =
        dynamic_cast<TH3D*>(inputDir->Get(
            (variationName + "/h3dDeltaThetaVsMassVsLambdaPt").c_str()));

    if (h3dDeltaThetaVsMassVsLambdaPt)
    {
        TDirectory* fits3DThetaDir =
            varOutDir->mkdir("MassFits_3D_DeltaTheta_LambdaPt");

        int nPtBins =
            h3dDeltaThetaVsMassVsLambdaPt->GetZaxis()->GetNbins();

        for (int iPt = 1; iPt <= nPtBins; ++iPt)
        {
            h3dDeltaThetaVsMassVsLambdaPt->GetZaxis()->SetRange(iPt, iPt);

            TH2D* h2dProj =
                dynamic_cast<TH2D*>(
                    h3dDeltaThetaVsMassVsLambdaPt->Project3D("yx"));

            if (!h2dProj)
                continue;

            int nThetaBins3D = h2dProj->GetXaxis()->GetNbins();

            for (int iTheta = 1; iTheta <= nThetaBins3D; ++iTheta)
            {
                h2dProj->GetXaxis()->SetRange(iTheta, iTheta);

                TH1D* hMassProj =
                    dynamic_cast<TH1D*>(h2dProj->ProjectionY(
                        Form("hMassProj_3D_theta%d_pt%d", iTheta, iPt),
                        iTheta, iTheta,
                        "e"));

                if (!hMassProj || hMassProj->GetEntries() < 10)
                    continue;

                TF1* fitFunc = new TF1(
                    Form("fit_3D_theta%d_pt%d", iTheta, iPt),
                    "[0]*TMath::Gaus(x,[1],[2],true) + [3] + [4]*x + [5]*x*x",
                    massMin, massMax);

                fitFunc->SetParameters(
                    hMassProj->GetMaximum(),
                    hMassProj->GetBinCenter(hMassProj->GetMaximumBin()),
                    0.005,
                    1.0, 0.0, 0.0);

                fitFunc->SetParLimits(1, massMin, massMax);
                fitFunc->SetParLimits(2, 1e-4, massRange);

                hMassProj->Fit(fitFunc, "QLS");

                fits3DThetaDir->cd();
                hMassProj->Write();
                fitFunc->Write();

                delete hMassProj;
                delete fitFunc;
            }

            delete h2dProj;
        }

        h3dDeltaThetaVsMassVsLambdaPt->GetZaxis()->SetRange(0, -1);
    }


    // ===============================
    // STEP 3: BACKGROUND SUBTRACTION
    // ===============================

    // Inputs required:
    //   hSignal_RingObservable              (Σ r_i in signal window)
    //   hSignal_RingObservableSquared       (Σ r_i^2 in signal window)
    //   hBkg_RingObservable                 (Σ r_i in sideband)
    //   hBkg_RingObservableSquared          (Σ r_i^2 in sideband)
    //   double alpha                        (background scaling factor)
    //
    // alpha = N_signal_sideband_norm / N_background_sideband
    // (computed previously from invariant mass fit or sideband integrals)


    // ---------------------------------------
    // Create output histogram
    // ---------------------------------------

    TH1D* hRingObservable_Subtracted =
        (TH1D*)hSignal_RingObservable->Clone("hRingObservable_Subtracted");

    hRingObservable_Subtracted->Reset();


    // ---------------------------------------
    // Loop over bins and subtract background
    // ---------------------------------------

    for (int iBin = 1; iBin <= hSignal_RingObservable->GetNbinsX(); ++iBin)
    {
        // ----- Signal -----
        double S  = hSignal_RingObservable->GetBinContent(iBin);          // Σ r_i
        double S2 = hSignal_RingObservableSquared->GetBinContent(iBin);   // Σ r_i^2

        // Statistical variance of weighted signal
        double varS = S2;   // variance = Σ r_i^2


        // ----- Background -----
        double B  = hBkg_RingObservable->GetBinContent(iBin);
        double B2 = hBkg_RingObservableSquared->GetBinContent(iBin);

        double varB = B2;   // variance = Σ r_i^2


        // ---------------------------------------
        // Background subtraction
        // ---------------------------------------

        double value = S - alpha * B;

        // ---------------------------------------
        // Error propagation
        //
        // Var(S - αB) = Var(S) + α² Var(B)
        //
        // (assuming S and B statistically independent)
        // ---------------------------------------

        double variance = varS + alpha * alpha * varB;

        double error = (variance > 0.0) ? std::sqrt(variance) : 0.0;


        // ---------------------------------------
        // Store result
        // ---------------------------------------

        hRingObservable_Subtracted->SetBinContent(iBin, value);
        hRingObservable_Subtracted->SetBinError(iBin, error);
    }



    std::cout << "Finished fits for variation " << variationName << "\n";
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
    std::cout << "\nStarting signalExtractionRing...\n";

    if (argc < 2)
    {
        std::cerr << "Usage: ./signalExtractionRing inputFile.root\n";
        return 1;
    }

    std::string inputFileName = argv[1];
    std::string suffix = extractSuffix(inputFileName);

    std::cout << "Input file: " << inputFileName << std::endl;
    std::cout << "Output suffix: " << suffix << std::endl;

    // ---------------------------------------------------------------------
    // Open input file
    // ---------------------------------------------------------------------

    TFile* inputFile = TFile::Open(inputFileName.c_str(), "READ");
    if (!inputFile || inputFile->IsZombie())
    {
        std::cerr << "ERROR: Could not open input file.\n";
        return 1;
    }

    TDirectory* mainDir =
        dynamic_cast<TDirectory*>(inputFile->Get("lambdajetpolarizationionsderived"));

    if (!mainDir)
    {
        std::cerr << "ERROR: Could not find folder lambdajetpolarizationionsderived\n";
        return 1;
    }

    // ---------------------------------------------------------------------
    // Create output file
    // ---------------------------------------------------------------------

    std::string outputFileName = "SignalExtraction_" + suffix + ".root";
    TFile* outputFile = TFile::Open(outputFileName.c_str(), "RECREATE");

    std::cout << "Output file: " << outputFileName << std::endl;

    // ---------------------------------------------------------------------
    // Observable variations
    // ---------------------------------------------------------------------

    std::vector<std::string> variations = {
        "Ring",
        "RingKinematicCuts",
        "JetKinematicCuts",
        "JetAndLambdaKinematicCuts"
    };

    for (const auto& var : variations)
    {
        processVariation(mainDir, outputFile, var);
    }

    outputFile->Write();
    outputFile->Close();
    inputFile->Close();

    std::cout << "\nSignal extraction completed successfully.\n";

    return 0;
}
