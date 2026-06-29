// ============================================================
// makeCumulativeProfile.cxx
// Convert a differential TProfile2D (phiStar vs DCA) into
// a cumulative one where each (i, j) bin holds <R> for all
// entries with DCA > lower edge of bin j.
// A very short, very simple macro, that simply tells us what is the behavior of
// experimental signal with increasingly higher minDCAdau cuts, to test for robustness
// of the ring observable wrt the AEE effect's main source.
//
// Usage: root -l -b -q 'makeCumulativeDCAdauProfile.cxx("path/to/ConsumerResults_SUFFIX.root")'
// ============================================================

#include <iostream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TProfile2D.h"
#include "TH2F.h"

void makeCumulativeDCAdauProfile(const char* inFileStr)
{
    // --- Dynamic Output Filename Drawing ---
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

    // --- Retrieve the differential TProfile2D ---
    std::string histoPath = "lambdajetpolarizationionsderived/HelicityEfficiencyQA/pFakePolSignalJet_PhiStarVsDCAProLike";
    TProfile2D* pDiff = (TProfile2D*)fIn->Get(histoPath.c_str());
    if (!pDiff) {
        std::cerr << " ERROR: Target histogram not found: " << histoPath << std::endl;
        fIn->Close();
        return;
    }

    int nPhi = pDiff->GetNbinsX();
    int nDCA = pDiff->GetNbinsY();

    // --- Create a single condensed output histogram ---
    // Stores Mean in BinContent and Error in BinError
    TH2F* hCumul = new TH2F("hCumul_PhiStarVsMinDCAProLike",
        "Cumul. ring signal (DCA_{pos}>cut); #phi^{*}; min DCA_{pos} (cm)",
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
            // Number of entries in differential bin (i, j)
            double n_ij = pDiff->GetBinEntries(i, j);

            // Mean ring observable in differential bin (i, j)
            double R_ij = pDiff->GetBinContent(i, j);

            // Standard error of the mean in differential bin (i, j)
            // GetBinError returns SEM = sigma/sqrt(n) for TProfile
            double sem_ij = pDiff->GetBinError(i, j);

            // Accumulate
            sum_N  += n_ij;
            sum_NR += n_ij * R_ij;
            // Variance of cumulative mean contribution from this bin:
            // Var(n*R_bar / N_tot) = n^2 * Var(R_bar) / N_tot^2
            //                      = n^2 * sem^2 / N_tot^2
            // We store n^2 * sem^2 here and divide by N_tot^2 at the end.
            sum_Nv += n_ij * n_ij * sem_ij * sem_ij;

            if (sum_N > 0.0) {
                double cumul_mean = sum_NR / sum_N;
                // Error: sqrt(sum_Nv) / sum_N
                double cumul_err  = std::sqrt(sum_Nv) / sum_N;

                // Condensing into single output histogram:
                hCumul->SetBinContent(i, j, cumul_mean);
                hCumul->SetBinError(i, j, cumul_err);
            }
        }
    }

    // --- Write output ---
    TFile* fOut = TFile::Open(outFileStr.c_str(), "RECREATE");
    hCumul->Write();
    
    std::cout << " Successfully generated and saved cumulative profile.\n Thank you for waiting!\n";
    std::cout << "=======================================================\n\n";

    fOut->Close();
    fIn->Close();
}