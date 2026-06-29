/*
* ============================================================
* makeCumulativeDCAdauProfile.cxx
* Convert differential TProfile2Ds (phiStar vs DCA) into
* cumulative ones where each (i, j) bin holds <R> for all
* entries with DCA > lower edge of bin j.
* A very short, very simple macro, that simply tells us what is the behavior of
* experimental signal with increasingly higher minDCAdau cuts, to test for robustness
* of the ring observable wrt the AEE effect's main source.
*
* Usage: root -l -b -q 'makeCumulativeDCAdauProfile.cxx("path/to/ConsumerResults_SUFFIX.root")'
* ============================================================
*/

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TString.h"

// Declaring helper functions to process all TProfiles and THDs that could benefit from accumulation to draw minDCAcut plots:

// ---------------------------------------------------------
// Helper 1: Core TProfile2D Processing
// ---------------------------------------------------------
TH2F* MakeCumulativeProfile2D(TProfile2D* pDiff, const char* outName, const char* outTitle) {
    int nPhi = pDiff->GetNbinsX();
    int nDCA = pDiff->GetNbinsY();

    TH2F* hCumul = new TH2F(outName, outTitle,
        nPhi, pDiff->GetXaxis()->GetXmin(), pDiff->GetXaxis()->GetXmax(),
        nDCA, pDiff->GetYaxis()->GetXmin(), pDiff->GetYaxis()->GetXmax());

    // --- Build cumulative from high DCA downward ---
    // At each step j (decreasing), we add the entries in differential
    // bin (i, j) to the running sum. The result at (i, j) represents
    // <R> for all candidates with DCA > lower_edge(bin j).
    for (int i = 1; i <= nPhi; i++) {
        // Some accumulators for error propagation:
        double sum_N   = 0.0; // Running total: sum of entry counts
        double sum_NR  = 0.0; // Running total: sum of n_{ik} * mean_R_{ik}
        double sum_Nv  = 0.0; // Running total: sum of n_{ik}^2 * SEM_{ik}^2

        for (int j = nDCA; j >= 1; j--) {
            double n_ij = pDiff->GetBinEntries(i, j);
            double R_ij = pDiff->GetBinContent(i, j);
            double sem_ij = pDiff->GetBinError(i, j);

            sum_N  += n_ij;
            sum_NR += n_ij * R_ij;
            // Variance of cumulative mean contribution from this bin:
            sum_Nv += n_ij * n_ij * sem_ij * sem_ij;

            if (sum_N > 0.0) {
                double cumul_mean = sum_NR / sum_N;
                double cumul_err  = std::sqrt(sum_Nv) / sum_N;

                // Condensing into single output histogram:
                hCumul->SetBinContent(i, j, cumul_mean);
                hCumul->SetBinError(i, j, cumul_err);
            }
        }
    }
    return hCumul;
}

// ---------------------------------------------------------
// Helper 2: Core TH2D (Counts) Processing
// ---------------------------------------------------------
TH2D* MakeCumulativeTH2D(TH2D* hDiff, const char* outName, const char* outTitle) {
    int nPhi = hDiff->GetNbinsX();
    int nDCA = hDiff->GetNbinsY();

    TH2D* hCumul = (TH2D*)hDiff->Clone(outName);
    hCumul->SetTitle(outTitle);
    hCumul->Reset(); // Clear old contents

    for (int i = 1; i <= nPhi; i++) {
        double sum_Counts = 0.0;
        double sum_Err2   = 0.0;
        for (int j = nDCA; j >= 1; j--) {
            sum_Counts += hDiff->GetBinContent(i, j);
            double err = hDiff->GetBinError(i, j);
            sum_Err2   += err * err;
            
            hCumul->SetBinContent(i, j, sum_Counts);
            hCumul->SetBinError(i, j, std::sqrt(sum_Err2));
        }
    }
    return hCumul;
}

// ---------------------------------------------------------
// Helper 3: TProfile3D Projection & Processing
// ---------------------------------------------------------
// Notice that the TProfile3D is actually just a set of two TProfile2Ds, one for each Jet or Lambda eta:
// Bin 1 is -0.9 to 0 (Negative Eta), Bin 2 is 0 to 0.9 (Positive Eta)
void ProcessProfile3D(TProfile3D* p3D, TFile* fOut, const std::string& baseName) {
    // Process Negative Eta (Z bin 1)
    p3D->GetZaxis()->SetRange(1, 1);
    TProfile2D* p2D_neg = p3D->Project3DProfile("xy_neg");
    TH2F* hCumul_neg = MakeCumulativeProfile2D(p2D_neg, 
        Form("hCumul_%s_NegEta", baseName.c_str()), 
        Form("Cumul. %s (Eta < 0); #phi^{*}; min DCA", baseName.c_str()));
    
    // Process Positive Eta (Z bin 2)
    p3D->GetZaxis()->SetRange(2, 2);
    TProfile2D* p2D_pos = p3D->Project3DProfile("xy_pos");
    TH2F* hCumul_pos = MakeCumulativeProfile2D(p2D_pos, 
        Form("hCumul_%s_PosEta", baseName.c_str()), 
        Form("Cumul. %s (Eta > 0); #phi^{*}; min DCA", baseName.c_str()));

    fOut->cd();
    hCumul_neg->Write();
    hCumul_pos->Write();

    delete p2D_neg; delete p2D_pos;
}


// ---------------------------------------------------------
// Main Function
// ---------------------------------------------------------
void makeCumulativeDCAdauProfile(const char* inFileStr)
{
    // --- Dynamic Output Filename ---
    std::string inPath(inFileStr);
    std::string directory = "";
    std::string filename = inPath;

    size_t lastSlash = inPath.find_last_of('/');
    if (lastSlash != std::string::npos) {
        directory = inPath.substr(0, lastSlash + 1);
        filename = inPath.substr(lastSlash + 1);
    }

    std::string prefixToReplace = "ConsumerResults_";
    size_t pos = filename.find(prefixToReplace);
    if (pos != std::string::npos) {
        filename.replace(pos, prefixToReplace.length(), "CumulativeProfiles_");
    } else {
        filename = "CumulativeProfiles_" + filename;
    }

    std::string outFileStr = directory + filename;

    std::cout << "\n=======================================================\n";
    std::cout << " Starting Cumulative Profile Drawing\n";
    std::cout << "=======================================================\n";
    std::cout << " Input File:  " << inPath << "\n";
    std::cout << " Output File: " << outFileStr << "\n";

    // --- Open Input File ---
    TFile* fIn = TFile::Open(inPath.c_str(), "READ");
    if (!fIn || fIn->IsZombie()) {
        std::cerr << " ERROR: cannot open " << inPath << std::endl;
        return;
    }

    TFile* fOut = TFile::Open(outFileStr.c_str(), "RECREATE");

    std::string baseDir = "lambdajetpolarizationionsderived/HelicityEfficiencyQA/";

    // 1. Process TH2D (Counts)
    std::vector<std::string> countsList = {
        "hFakePolCountsJet_PhiStarVsDCAdau",
        "hFakePolCountsJet_PhiStarVsDCAProLike",
        "hFakePolCountsJet_PhiStarVsDCAPiLike"
    };

    std::cout << " -> Processing TH2D Count Histograms...\n";
    for(const auto& name : countsList) {
        TH2D* hIn = (TH2D*)fIn->Get((baseDir + name).c_str());
        if(hIn) {
            TH2D* hCumul = MakeCumulativeTH2D(hIn, ("hCumul_" + name).c_str(), hIn->GetTitle());
            fOut->cd(); hCumul->Write();
        } else {
             std::cerr << "    [Warning] Missing: " << name << std::endl;
        }
    }

    // 2. Process TProfile2D (Signals)
    std::vector<std::string> prof2DList = {
        "pFakePolSignalJet_PhiStarVsDCAdau",
        "pFakePolSignalJet_PhiStarVsDCAProLike",
        "pFakePolSignalJet_PhiStarVsDCAPiLike"
    };

    std::cout << " -> Processing TProfile2D Signal Histograms...\n";
    for(const auto& name : prof2DList) {
        TProfile2D* pIn = (TProfile2D*)fIn->Get((baseDir + name).c_str());
        if(pIn) {
            TH2F* pCumul = MakeCumulativeProfile2D(pIn, ("hCumul_" + name).c_str(), pIn->GetTitle());
            fOut->cd(); pCumul->Write();
        } else {
             std::cerr << "    [Warning] Missing: " << name << std::endl;
        }
    }

    // 3. Process TProfile3D (Signals with Eta splits)
    std::vector<std::string> prof3DList = {
        "pFakePolSignalJet_PhiStarVsDCAdauVsEtaJet",
        "pFakePolSignalJet_PhiStarVsDCAdauVsEtaLambda",
        "pFakePolSignalJet_PhiStarVsDCAProLikeVsEtaJet",
        "pFakePolSignalJet_PhiStarVsDCAProLikeVsEtaLambda",
        "pFakePolSignalJet_PhiStarVsDCAPiLikeVsEtaJet",
        "pFakePolSignalJet_PhiStarVsDCAPiLikeVsEtaLambda"
    };

    std::cout << " -> Processing TProfile3D Eta-Split Histograms...\n";
    for(const auto& name : prof3DList) {
        TProfile3D* p3DIn = (TProfile3D*)fIn->Get((baseDir + name).c_str());
        if(p3DIn) {
            ProcessProfile3D(p3DIn, fOut, name);
        } else {
             std::cerr << "    [Warning] Missing: " << name << std::endl;
        }
    }

    std::cout << " Successfully generated and saved all cumulative profiles.\n Thank you for waiting!\n";
    std::cout << "=======================================================\n\n";

    fOut->Close();
    fIn->Close();
}