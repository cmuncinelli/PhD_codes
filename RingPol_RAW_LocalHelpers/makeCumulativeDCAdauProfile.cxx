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
#include "TH1F.h"
#include "TString.h"

// Declaring helper functions to process all TProfiles and THDs that could benefit from accumulation to draw minDCAcut plots:

// ---------------------------------------------------------
// Helper 1: Core TProfile2D Processing (Outputs 2D)
// ---------------------------------------------------------
TH2F* MakeCumulativeProfile2D(TProfile2D* pDiff, const char* outName, const char* outTitle) {
    int nPhi = pDiff->GetNbinsX();
    int nDCA = pDiff->GetNbinsY();

    // Transpose the output: X-axis becomes DCA (from input Y), Y-axis becomes phi* (from input X)
    TH2F* hCumul = new TH2F(outName, outTitle,
        nDCA, pDiff->GetYaxis()->GetXmin(), pDiff->GetYaxis()->GetXmax(),
        nPhi, pDiff->GetXaxis()->GetXmin(), pDiff->GetXaxis()->GetXmax());

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
            // Get the global bin index first!
            int global_bin = pDiff->GetBin(i, j);

            double n_ij = pDiff->GetBinEntries(global_bin);
            double R_ij = pDiff->GetBinContent(global_bin);
            double sem_ij = pDiff->GetBinError(global_bin);

            sum_N  += n_ij;
            sum_NR += n_ij * R_ij;
            // Variance of cumulative mean contribution from this bin:
            sum_Nv += n_ij * n_ij * sem_ij * sem_ij;

            if (sum_N > 0.0) {
                double cumul_mean = sum_NR / sum_N;
                double cumul_err  = std::sqrt(sum_Nv) / sum_N;

                // SetBinContent(X, Y) -> We put DCA (j) on X and phi* (i) on Y
                hCumul->SetBinContent(j, i, cumul_mean);
                hCumul->SetBinError(j, i, cumul_err);
            }
        }
    }
    return hCumul;
}

// ---------------------------------------------------------
// Helper 2: Core TProfile2D Processing to 1D (Integrated over phi*)
// ---------------------------------------------------------
TH1F* MakeCumulativeProfile1D(TProfile2D* pDiff, const char* outName, const char* outTitle) {
    int nPhi = pDiff->GetNbinsX();
    int nDCA = pDiff->GetNbinsY();

    // The X-axis of this 1D histogram corresponds to the Y-axis (DCA) of the 2D profile
    TH1F* hCumul1D = new TH1F(outName, outTitle,
        nDCA, pDiff->GetYaxis()->GetXmin(), pDiff->GetYaxis()->GetXmax());
    
    hCumul1D->GetXaxis()->SetTitle("min DCA cut");
    hCumul1D->GetYaxis()->SetTitle("Integrated <R>");

    // Accumulators spanning across ALL phi* bins
    double total_sum_N  = 0.0; 
    double total_sum_NR = 0.0; 
    double total_sum_Nv = 0.0; 

    // Build cumulative from high DCA downward
    for (int j = nDCA; j >= 1; j--) {
        // Integrate over all phi* bins (i) for this specific DCA bin (j)
        for (int i = 1; i <= nPhi; i++) {
            // Get the global bin index first!
            int global_bin = pDiff->GetBin(i, j);
            
            double n_ij = pDiff->GetBinEntries(global_bin);
            double R_ij = pDiff->GetBinContent(global_bin);
            double sem_ij = pDiff->GetBinError(global_bin);

            total_sum_N  += n_ij;
            total_sum_NR += n_ij * R_ij;
            total_sum_Nv += n_ij * n_ij * sem_ij * sem_ij;
        }

        // Calculate and set the integrated bin content for this min DCA cut
        if (total_sum_N > 0.0) {
            double cumul_mean = total_sum_NR / total_sum_N;
            double cumul_err  = std::sqrt(total_sum_Nv) / total_sum_N;

            hCumul1D->SetBinContent(j, cumul_mean);
            hCumul1D->SetBinError(j, cumul_err);
        }
    }
    return hCumul1D;
}


// ---------------------------------------------------------
// Helper 3: Core TH2D (Counts) Processing
// ---------------------------------------------------------
TH2D* MakeCumulativeTH2D(TH2D* hDiff, const char* outName, const char* outTitle) {
    int nPhi = hDiff->GetNbinsX();
    int nDCA = hDiff->GetNbinsY();

    // Transpose output to get X = DCA, Y = phi*
    // (phi* looks better in the Y axis. It is usually a supporting variable, not the main one here)
    TH2D* hCumul = new TH2D(outName, outTitle,
        nDCA, hDiff->GetYaxis()->GetXmin(), hDiff->GetYaxis()->GetXmax(),
        nPhi, hDiff->GetXaxis()->GetXmin(), hDiff->GetXaxis()->GetXmax());

    for (int i = 1; i <= nPhi; i++) {
        double sum_Counts = 0.0;
        double sum_Err2   = 0.0;
        for (int j = nDCA; j >= 1; j--) {
            int global_bin = hDiff->GetBin(i, j);
            sum_Counts += hDiff->GetBinContent(global_bin);
            double err = hDiff->GetBinError(global_bin);
            sum_Err2   += err * err;
            
            // Populate transposed coordinates
            hCumul->SetBinContent(j, i, sum_Counts);
            hCumul->SetBinError(j, i, std::sqrt(sum_Err2));
        }
    }
    return hCumul;
}

// ---------------------------------------------------------
// Helper 4: TProfile3D Projection & Processing
// ---------------------------------------------------------
// Notice that the TProfile3D is actually just a set of two TProfile2Ds, one for each Jet or Lambda eta:
// Bin 1 is -0.9 to 0 (Negative Eta), Bin 2 is 0 to 0.9 (Positive Eta)
// Extracts All, Positive, and Negative Eta variations and generates both 2D and 1D cumulatives
void ProcessProfile3D(TProfile3D* p3D, TDirectory* dir2D, TDirectory* dir1D, const std::string& baseName) {
    
    // Pass "yx" so the 2D projection maintains native mapping (X = phi*, Y = DCA)
    // We assign names manually afterward.
    
    // 1. All Eta (Project across all Z bins)
    p3D->GetZaxis()->SetRange(1, p3D->GetNbinsZ());
    TProfile2D* p2D_all = p3D->Project3DProfile("yx");
    p2D_all->SetName(Form("%s_xy_all", baseName.c_str()));
    
    // 2. Negative Eta (Z bin 1: -0.9 to 0)
    p3D->GetZaxis()->SetRange(1, 1);
    TProfile2D* p2D_neg = p3D->Project3DProfile("yx");
    p2D_neg->SetName(Form("%s_xy_neg", baseName.c_str()));
    
    // 3. Positive Eta (Z bin 2: 0 to 0.9)
    p3D->GetZaxis()->SetRange(2, 2);
    TProfile2D* p2D_pos = p3D->Project3DProfile("yx");
    p2D_pos->SetName(Form("%s_xy_pos", baseName.c_str()));

    // Resolve Precise DCA Label
    std::string dcaLabel = "min DCA";
    if (baseName.find("DCAdau") != std::string::npos) dcaLabel = "min DCA_{V0 Daughters}";
    else if (baseName.find("DCAProLike") != std::string::npos) dcaLabel = "min DCA_{PosPV}";
    else if (baseName.find("DCAPiLike") != std::string::npos) dcaLabel = "min DCA_{NegPV}";

    // Generate 2D Cumulatives (Transposed now: X = min DCA, Y = phi*)
    TH2F* hCumul2D_all = MakeCumulativeProfile2D(p2D_all, Form("hCumul2D_%s_AllEta", baseName.c_str()), Form("Cumul. %s (All Eta); %s; #phi^{*}", baseName.c_str(), dcaLabel.c_str()));
    TH2F* hCumul2D_neg = MakeCumulativeProfile2D(p2D_neg, Form("hCumul2D_%s_NegEta", baseName.c_str()), Form("Cumul. %s (Eta < 0); %s; #phi^{*}", baseName.c_str(), dcaLabel.c_str()));
    TH2F* hCumul2D_pos = MakeCumulativeProfile2D(p2D_pos, Form("hCumul2D_%s_PosEta", baseName.c_str()), Form("Cumul. %s (Eta > 0); %s; #phi^{*}", baseName.c_str(), dcaLabel.c_str()));

    // Generate 1D Integrated Cumulatives
    TH1F* hCumul1D_all = MakeCumulativeProfile1D(p2D_all, Form("hCumul1D_%s_AllEta", baseName.c_str()), Form("Integrated Cumul. %s (All Eta); %s; Integrated <R>", baseName.c_str(), dcaLabel.c_str()));
    TH1F* hCumul1D_neg = MakeCumulativeProfile1D(p2D_neg, Form("hCumul1D_%s_NegEta", baseName.c_str()), Form("Integrated Cumul. %s (Eta < 0); %s; Integrated <R>", baseName.c_str(), dcaLabel.c_str()));
    TH1F* hCumul1D_pos = MakeCumulativeProfile1D(p2D_pos, Form("hCumul1D_%s_PosEta", baseName.c_str()), Form("Integrated Cumul. %s (Eta > 0); %s; Integrated <R>", baseName.c_str(), dcaLabel.c_str()));

    if (dir2D) { dir2D->cd(); hCumul2D_all->Write(); hCumul2D_neg->Write(); hCumul2D_pos->Write(); }
    if (dir1D) { dir1D->cd(); hCumul1D_all->Write(); hCumul1D_neg->Write(); hCumul1D_pos->Write(); }

    delete p2D_all; delete p2D_neg; delete p2D_pos;
    delete hCumul2D_all; delete hCumul2D_neg; delete hCumul2D_pos;
    delete hCumul1D_all; delete hCumul1D_neg; delete hCumul1D_pos;
}


// ---------------------------------------------------------
// Main Function
// ---------------------------------------------------------
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <inputFilePath>\n";
        return 1;
    }
    
    const char* inFileStr = argv[1];

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
        return 1;
    }

    TFile* fOut = TFile::Open(outFileStr.c_str(), "RECREATE");
    
    // --- Create Organized Subfolders ---
    TDirectory* dirCounts = fOut->mkdir("Cumulative_Counts");
    TDirectory* dir2D = fOut->mkdir("Cumulative_2D");
    TDirectory* dir1D = fOut->mkdir("Cumulative_1D");

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
            dirCounts->cd(); hCumul->Write();
        } else {
             std::cerr << "    [Warning] Missing: " << name << std::endl;
        }
    }

    // 2. Process TProfile2D (Signals) -> Generates both 2D and 1D versions
    std::vector<std::string> prof2DList = {
        "pFakePolSignalJet_PhiStarVsDCAdau",
        "pFakePolSignalJet_PhiStarVsDCAProLike",
        "pFakePolSignalJet_PhiStarVsDCAPiLike"
    };

    std::cout << " -> Processing TProfile2D Signal Histograms...\n";
    for(const auto& name : prof2DList) {
        TProfile2D* pIn = (TProfile2D*)fIn->Get((baseDir + name).c_str());
        if(pIn) {
            // Smart labeling
            std::string dcaLabel = "min DCA";
            if (name.find("DCAdau") != std::string::npos) dcaLabel = "min DCA_{V0 Daughters}";
            else if (name.find("DCAProLike") != std::string::npos) dcaLabel = "min DCA_{PosPV}";
            else if (name.find("DCAPiLike") != std::string::npos) dcaLabel = "min DCA_{NegPV}";

            TH2F* pCumul2D = MakeCumulativeProfile2D(pIn, ("hCumul2D_" + name).c_str(), ("Cumul. " + name + "; " + dcaLabel + "; #phi^{*}").c_str());
            dir2D->cd(); pCumul2D->Write();

            TH1F* pCumul1D = MakeCumulativeProfile1D(pIn, ("hCumul1D_" + name).c_str(), ("Integrated Cumul. " + name + "; " + dcaLabel + "; Integrated <R>").c_str());
            dir1D->cd(); pCumul1D->Write();

            // Deleting to avoid leaks:
            delete pCumul2D; delete pCumul1D;
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
            ProcessProfile3D(p3DIn, dir2D, dir1D, name);
        } else {
             std::cerr << "    [Warning] Missing: " << name << std::endl;
        }
    }

    std::cout << " Successfully generated and saved all cumulative profiles.\n Thank you for waiting!\n";
    std::cout << "=======================================================\n\n";

    fOut->Write();
    fOut->Close();

    fIn->Close();

    delete fOut;
    delete fIn;

    return 0;
}