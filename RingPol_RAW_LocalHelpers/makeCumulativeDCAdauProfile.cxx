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
* Input: the consumer's FakePolDrivers/PhiStar/DCACorrelations/LeadJet/ folder, booked only with
*        doFakePolDiagnosticsQA. Without it there is nothing to do: the macro says so and writes no output file.
*
* Usage: ./makeCumulativeDCAdauProfile.exe <inputFilePath> <outputDirectory>
* (outputDirectory is created recursively if it does not already exist -- this
*  keeps post-processing output out of results_consumer/, see run_all_wagons.sh)
* ============================================================
*/

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TH2D.h"
#include "TArrayD.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"

// ---------------------------------------------------------
// Input location and naming tables
// ---------------------------------------------------------
// The consumer fills the (phi*, DCA) correlations for the leading jet only. The proxy folder name also goes into every
// title (<R>_{LeadJet}, #eta_{LeadJet}), so the "Proxy" tokens in the object names are never ambiguous.
const std::string kProxyFolder = "LeadJet";
const std::string kInputDir = "lambdajetpolarizationionsderived/FakePolDrivers/PhiStar/DCACorrelations/" + kProxyFolder;
const std::string kRingLabel = "<#it{R}>_{" + kProxyFolder + "}"; // Same convention as the consumer titles

// One entry per DCA definition. The tag is the consumer's token, used for the input and the output names alike.
// The daughters are told apart by species (proton-like, pion-like), not by charge: for an AntiLambda the proton-like one is negative.
struct DCASpec { std::string tag; std::string axisTitle; std::string legendLabel; };
const std::vector<DCASpec> kDCASpecs = {
    { "DCAdau",     "min DCA_{V0 daughters} (cm)", "DCA_{V0 daughters}" },
    { "DCAProLike", "min DCA_{p-like, PV} (cm)",   "DCA_{p-like, PV}"   },
    { "DCAPiLike",  "min DCA_{#pi-like, PV} (cm)", "DCA_{#pi-like, PV}" }
};

// One entry per eta classifier of the TProfile3Ds (the third axis carries the sign of this eta):
struct EtaSplitSpec { std::string tag; std::string etaLabel; };
const std::vector<EtaSplitSpec> kEtaSplitSpecs = {
    { "ProxyEta",  "#eta_{" + kProxyFolder + "}" },
    { "EtaLambda", "#eta_{#Lambda}" }
};

// Output names follow the consumer convention: dimension prefix, then the quantity (RingObservable for <R>, Counts for
// entries), then <X>Vs<Y>. "Cumul" marks the cumulative versions, whose X axis is the minimum DCA cut.
// Shared by the producers below and by DrawComparisonCanvases, which reads the 1D ones back.
// sliceSuffix is empty for the plain (phi*, DCA) profiles, or "_<eta tag><All|Neg|Pos>" for the eta slices.
std::string CumulCountsName(const DCASpec& dca) { return "h2dCumulCountsMin" + dca.tag + "VsPhiStar"; }
std::string Cumul2DName(const DCASpec& dca, const std::string& sliceSuffix = "") { return "p2dCumulRingObservableMin" + dca.tag + "VsPhiStar" + sliceSuffix; }
std::string Cumul1DName(const DCASpec& dca, const std::string& sliceSuffix = "") { return "pCumulRingObservableMin" + dca.tag + sliceSuffix; }
std::string SliceSuffix(const EtaSplitSpec& split, const std::string& slice) { return "_" + split.tag + slice; }

// Declaring helper functions to process all TProfiles and THDs that could benefit from accumulation to draw minDCAcut plots:


// ---------------------------------------------------------
// Helper 1: Core TProfile2D Processing (Outputs 2D)
// ---------------------------------------------------------
TProfile2D* MakeCumulativeProfile2D(TProfile2D* pDiff, const char* outName, const char* outTitle) {  
    // Null input is a legitimate outcome upstream (a projection that could not be built), so this
    // returns nullptr rather than dereferencing it. Callers must check the result before writing.
    if (!pDiff) return nullptr;

    int nPhi = pDiff->GetNbinsX();
    int nDCA = pDiff->GetNbinsY();

    // Transpose the output: X-axis becomes DCA (from input Y), Y-axis becomes phi* (from input X)
    // Kept as a TProfile2D (not a TH2F) so the per-bin sum/entries/sum-of-squares survive,
    // which is what lets TBrowser rebin this histogram correctly on the fly.
    TProfile2D* hCumul = new TProfile2D(outName, outTitle,
        nDCA, pDiff->GetYaxis()->GetXmin(), pDiff->GetYaxis()->GetXmax(),
        nPhi, pDiff->GetXaxis()->GetXmin(), pDiff->GetXaxis()->GetXmax());

    // --- Build cumulative from high DCA downward ---
    // At each step j (decreasing), we add the entries in differential
    // bin (i, j) to the running sum. The result at (i, j) represents
    // <R> for all candidates with DCA > lower_edge(bin j).
    for (int i = 1; i <= nPhi; i++) {
        // Accumulators for error propagation now hold the RAW profile statistics (entries, sum of Y,
        // sum of Y^2) rather than a precomputed mean/error, so they can be
        // written straight into the TProfile2D's own per-bin storage below:
        double sum_N  = 0.0; // Running total: sum of entry counts         -> bin entries
        double sum_Y  = 0.0; // Running total: sum of n_{ik} * mean_R_{ik} -> bin content (raw sum of Y)
        double sum_Y2 = 0.0; // Running total: sum of Y^2 across bins      -> bin sum-of-squares

        for (int j = nDCA; j >= 1; j--) {
            // Get the global bin index first!
            int global_bin = pDiff->GetBin(i, j);

            double n_ij  = pDiff->GetBinEntries(global_bin);
            double R_ij  = pDiff->GetBinContent(global_bin);
            double Y2_ij = pDiff->GetSumw2()->At(global_bin); // exact sum(Y^2) stored by the source profile

            sum_N  += n_ij;
            sum_Y  += n_ij * R_ij;
            sum_Y2 += Y2_ij;

            if (sum_N > 0.0) {
                // SetBinContent(X, Y) -> We put DCA (j) on X and phi* (i) on Y
                // Write the raw statistics (NOT the mean/error!) - TProfile2D
                // derives <R> and its error from these automatically.
                int dest_bin = hCumul->GetBin(j, i);
                hCumul->SetBinEntries(dest_bin, sum_N);
                hCumul->SetBinContent(dest_bin, sum_Y);
                hCumul->GetSumw2()->SetAt(sum_Y2, dest_bin);
            }
        }
    }
    return hCumul;
}

// ---------------------------------------------------------
// Helper 2: Core TProfile2D Processing to 1D (Integrated over phi*)
// ---------------------------------------------------------
TProfile* MakeCumulativeProfile1D(TProfile2D* pDiff, const char* outName, const char* outTitle) {
    // Same contract as MakeCumulativeProfile2D: if nullptr in, then nullptr out.
    if (!pDiff) return nullptr;

    int nPhi = pDiff->GetNbinsX();
    int nDCA = pDiff->GetNbinsY();

    // The X-axis of this 1D histogram corresponds to the Y-axis (DCA) of the 2D profile
    // Kept as a TProfile (not a TH1F) for the same reason as MakeCumulativeProfile2D:
    // preserving entries/sum/sum-of-squares per bin allows correct TBrowser rebinning.
    TProfile* hCumul1D = new TProfile(outName, outTitle,
        nDCA, pDiff->GetYaxis()->GetXmin(), pDiff->GetYaxis()->GetXmax());

    // Accumulators spanning across ALL phi* bins, holding RAW profile statistics
    double total_sum_N  = 0.0; // -> bin entries
    double total_sum_Y  = 0.0; // -> bin content (raw sum of Y)
    double total_sum_Y2 = 0.0; // -> bin sum-of-squares

    // Build cumulative from high DCA downward
    for (int j = nDCA; j >= 1; j--) {
        // Integrate over all phi* bins (i) for this specific DCA bin (j)
        for (int i = 1; i <= nPhi; i++) {
            // Get the global bin index first!
            int global_bin = pDiff->GetBin(i, j);
            
            double n_ij  = pDiff->GetBinEntries(global_bin);
            double R_ij  = pDiff->GetBinContent(global_bin);
            double Y2_ij = pDiff->GetSumw2()->At(global_bin); // exact sum(Y^2) stored by the source profile

            total_sum_N  += n_ij;
            total_sum_Y  += n_ij * R_ij;
            total_sum_Y2 += Y2_ij;
        }

        // Write the raw statistics (NOT the mean/error!) for this min DCA cut
        if (total_sum_N > 0.0) {
            hCumul1D->SetBinEntries(j, total_sum_N);
            hCumul1D->SetBinContent(j, total_sum_Y);
            hCumul1D->GetSumw2()->SetAt(total_sum_Y2, j);
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
// Notice that the TProfile3D is actually just a set of two TProfile2Ds, one for each sign of the split eta:
// Bin 1 is -0.9 to 0 (Negative Eta), Bin 2 is 0 to 0.9 (Positive Eta)
// Extracts All, Positive, and Negative Eta variations and generates both 2D and 1D cumulatives
// -- Note: "All" is NOT a copy of the plain (phi*, DCA) profile. The projection excludes the under- and overflow of the
//    eta axis, so "All" drops the candidates with |eta| >= 0.9 (jets in particular can sit beyond it).
void ProcessProfile3D(TProfile3D* p3D, TDirectory* dir2D, TDirectory* dir1D, const DCASpec& dca, const EtaSplitSpec& split) {
    // Pass "yx" so the 2D projection maintains native mapping (X = phi*, Y = DCA)
    // We assign names manually afterward.
    //
    // Project3DProfile returns nullptr on a rejected option string or an unusable axis range, so each result is checked before it is named:
    // (Dereferencing it blind turns a recoverable ROOT complaint into a segfault, which the coordinator can only report as an opaque crash)
        // Thus, we create a small lambda to encapsulate everything:
    auto project = [&](int zLo, int zHi, const std::string& slice) -> TProfile2D* {
        p3D->GetZaxis()->SetRange(zLo, zHi);
        TProfile2D* p2D = p3D->Project3DProfile("yx");
        if (!p2D) {
            std::cerr << "    [Warning] Project3DProfile(\"yx\") failed for " << p3D->GetName()
                      << " (" << slice << "). Skipping this eta slice.\n";
            return nullptr;
        }
        p2D->SetName((std::string(p3D->GetName()) + "_yx_" + slice).c_str()); // Temporary: deleted below
        return p2D;
    };

    // The three slices: name suffix, Z bin range and title tag
    struct EtaSlice { std::string slice; int zLo; int zHi; std::string titleTag; };
    const std::vector<EtaSlice> slices = {
        { "All", 1, p3D->GetNbinsZ(), "all " + split.etaLabel }, // 1. All Eta (Project across all Z bins)
        { "Neg", 1, 1,                split.etaLabel + " < 0" },  // 2. Negative Eta (Z bin 1: -0.9 to 0)
        { "Pos", 2, 2,                split.etaLabel + " #geq 0" } // 3. Positive Eta (Z bin 2: 0 to 0.9)
    };

    for (const auto& sl : slices) {
        TProfile2D* p2D = project(sl.zLo, sl.zHi, sl.slice);
        const std::string suffix = SliceSuffix(split, sl.slice);

        // Generate the 2D (transposed now: X = min DCA, Y = phi*) and the 1D integrated cumulatives
        TProfile2D* hCumul2D = MakeCumulativeProfile2D(p2D, Cumul2DName(dca, suffix).c_str(),
            ("Cumulative " + kRingLabel + ", " + sl.titleTag + ";" + dca.axisTitle + ";#phi^{*}").c_str());
        TProfile* hCumul1D = MakeCumulativeProfile1D(p2D, Cumul1DName(dca, suffix).c_str(),
            ("Cumulative " + kRingLabel + ", #phi^{*}-integrated, " + sl.titleTag + ";" + dca.axisTitle + ";" + kRingLabel).c_str());

        // Written individually: a failed projection leaves its cumulative null, and the other two slices
        // are still perfectly good output.
        if (dir2D && hCumul2D) { dir2D->cd(); hCumul2D->Write(); }
        if (dir1D && hCumul1D) { dir1D->cd(); hCumul1D->Write(); }

        delete p2D; delete hCumul2D; delete hCumul1D;
    }

    // Restore the full Z range: the caller's TProfile3D is shared, and leaving it clipped to bin 2
    // would silently bias anything that reads the same histogram afterwards (if we ever intend to use it afterwards).
    p3D->GetZaxis()->SetRange(1, p3D->GetNbinsZ());
}



// ---------------------------------------------------------
// Helper 5: Draw comparison TCanvases from stored 1D cumulatives
// ---------------------------------------------------------
// Produces 11 canvases in total, saved to dirCanvases:
//   A (1): phi*-integrated, all three DCA types overlaid
//   B (3): leading-jet (proxy) eta splits (All/Neg/Pos), one canvas per DCA type
//   C (3): Lambda eta splits (All/Neg/Pos), one canvas per DCA type
//   D (4): fixed eta hemisphere (ProxyEta/EtaLambda x Pos/Neg), all three DCAs overlaid
void DrawComparisonCanvases(TDirectory* dir1D, TDirectory* dirCanvases) {

    // Retrieve a 1D cumulative from dir1D by name; warn if absent.
    auto Get1D = [&](const std::string& name) -> TProfile* {
        TProfile* h = (TProfile*)dir1D->Get(name.c_str());
        if (!h) std::cerr << "    [Warning] Canvas: missing " << name << "\n";
        return h;
    };

    // Apply line/marker style to a histogram
    auto Style = [](TProfile* h, int color, int marker) {
        if (!h) return;
        h->SetLineColor(color); h->SetMarkerColor(color);
        h->SetMarkerStyle(marker); h->SetMarkerSize(0.9);
    };

    // Draw a set of 1D histograms onto one canvas and write it to dirCanvases.
    // Automatically scales the Y-axis with a 15% margin around the global min/max.
    auto WriteCanvas = [&](const char* cName, const char* frameTitle,
                           const std::vector<TProfile*>& hv,
                           const std::vector<const char*>& lv) {
        TCanvas* c = new TCanvas(cName, cName, 800, 600);
        c->SetLeftMargin(0.14); c->SetBottomMargin(0.13);

        TLegend* leg = new TLegend(0.55, 0.70, 0.88, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);

        // 1. Find the global minimum and maximum (including error bars!)
        double yMin = 1e9;
        double yMax = -1e9;
        TProfile* firstValidHist = nullptr;

        for (size_t k = 0; k < hv.size(); k++) {
            TProfile* h = hv[k];
            if (!h) continue;
            
            // Keep track of the first valid histogram to use as our frame setter
            if (!firstValidHist) firstValidHist = h;
            leg->AddEntry(h, lv[k], "lp");

            // Scan bins to find the true min and max
            for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
                double val = h->GetBinContent(bin);
                double err = h->GetBinError(bin);
                
                // Optional: skip completely empty bins if they distort the scaling
                if (val == 0 && err == 0) continue; 
                
                if (val + err > yMax) yMax = val + err;
                if (val - err < yMin) yMin = val - err;
            }
        }

        // 2. Apply the dynamic limits if we found valid data
        if (firstValidHist) {
            firstValidHist->SetTitle(frameTitle);

            // yMin/yMax only get updated inside the scan above if at least one
            // bin had nonzero content or error; if every histogram was empty,
            // they are still at their +-1e9 sentinel values (yMax < yMin), and
            // SetRangeUser with an inverted range would silently render nothing.
            // In that case, skip it and let ROOT fall back to its own auto-range.
            if (yMax > yMin) {
                double range = yMax - yMin;
                if (range == 0) range = 1.0; // Fallback in case of flat lines

                // Add a 15% margin to top and bottom
                double margin = range * 0.15;
                firstValidHist->GetYaxis()->SetRangeUser(yMin - margin, yMax + margin);
            }

            // 3. Draw everything
            bool firstDrawn = false;
            for (TProfile* h : hv) {
                if (!h) continue;
                if (!firstDrawn) { 
                    h->Draw("E1"); 
                    firstDrawn = true; 
                } else {
                    h->Draw("E1 SAME");
                }
            }
            leg->Draw();
        }

        dirCanvases->cd(); 
        c->Write();
        delete leg; 
        delete c;
    };

    // Color/marker scheme: DCA type (Canvases A and D)
    const int kColDau = kRed+1,   kMkrDau = 20;
    const int kColPro = kAzure+2, kMkrPro = 21;
    const int kColPi  = kGreen+2, kMkrPi  = 22;
    // Color/marker scheme: eta slice (Canvases B and C)
    const int kColAll = kBlack,   kMkrAll = 20;
    const int kColNeg = kRed+1,   kMkrNeg = 25;
    const int kColPos = kAzure+2, kMkrPos = 24;

    // Markers and colours of the DCA-overlay canvases (A and D) follow kDCASpecs order: V0 daughters, p-like, pi-like
    const std::vector<int> kColDCA = { kColDau, kColPro, kColPi };
    const std::vector<int> kMkrDCA = { kMkrDau, kMkrPro, kMkrPi };
    std::vector<const char*> dcaLegends;
    for (const auto& d : kDCASpecs) dcaLegends.push_back(d.legendLabel.c_str());

    // Fetches and styles the three DCA types for one name suffix (shared by canvases A and D)
    auto GetAllDCA = [&](const std::string& sliceSuffix) {
        std::vector<TProfile*> hv;
        for (size_t k = 0; k < kDCASpecs.size(); ++k) {
            TProfile* h = Get1D(Cumul1DName(kDCASpecs[k], sliceSuffix));
            Style(h, kColDCA[k], kMkrDCA[k]);
            hv.push_back(h);
        }
        return hv;
    };

    // A. phi*-integrated, all three DCA types overlaid
    WriteCanvas("cCumulRingObservable_AllDCA",
                ("Cumulative " + kRingLabel + " vs min DCA (#phi^{*}-integrated); min DCA cut (cm); " + kRingLabel).c_str(),
                GetAllDCA(""), dcaLegends);

    // B/C. Eta-slice variations per DCA type (2 eta classifiers x 3 DCA types = 6 canvases)
    for (const auto& d : kDCASpecs) {
        for (const auto& e : kEtaSplitSpecs) {
            TProfile* hAll = Get1D(Cumul1DName(d, SliceSuffix(e, "All")));
            TProfile* hNeg = Get1D(Cumul1DName(d, SliceSuffix(e, "Neg")));
            TProfile* hPos = Get1D(Cumul1DName(d, SliceSuffix(e, "Pos")));
            Style(hAll, kColAll, kMkrAll); Style(hNeg, kColNeg, kMkrNeg); Style(hPos, kColPos, kMkrPos);

            const std::string cName  = "cCumulRingObservable_" + d.tag + "_" + e.tag + "Split";
            const std::string cTitle = "Cumulative " + kRingLabel + " vs " + d.axisTitle + " (" + e.etaLabel
                                       + " splits);" + d.axisTitle + ";" + kRingLabel;
            const std::string legAll = "All " + e.etaLabel, legNeg = e.etaLabel + " < 0", legPos = e.etaLabel + " #geq 0";
            WriteCanvas(cName.c_str(), cTitle.c_str(),
                        { hAll, hNeg, hPos },
                        { legAll.c_str(), legNeg.c_str(), legPos.c_str() });
        }
    }

    // D. Fixed eta hemisphere, all three DCA types overlaid (4 canvases)
    // Isolates which DCA type is responsible for the AEE signal in each eta hemisphere.
    for (const auto& e : kEtaSplitSpecs) {
        for (const std::string slice : { "Pos", "Neg" }) {
            const std::string etaRange = e.etaLabel + (slice == "Pos" ? " #geq 0" : " < 0");
            const std::string cName  = "cCumulRingObservable_AllDCA_" + e.tag + slice;
            const std::string cTitle = "Cumulative " + kRingLabel + " vs min DCA cut (" + etaRange + "); min DCA cut (cm); " + kRingLabel;
            WriteCanvas(cName.c_str(), cTitle.c_str(), GetAllDCA(SliceSuffix(e, slice)), dcaLegends);
        }
    }
}


// ---------------------------------------------------------
// Main Function
// ---------------------------------------------------------
int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <inputFilePath> <outputDirectory>\n";
        return 1;
    }
    
    const char* inFileStr = argv[1];
    const char* outDirStr = argv[2];

    // --- Dynamic Output Filename ---
    // The output directory is now taken explicitly from argv[2] (instead of being
    // inferred from the input file's directory), so the filename is derived from
    // the input's basename only.
    std::string inPath(inFileStr);
    std::string filename = inPath;

    size_t lastSlash = inPath.find_last_of('/');
    if (lastSlash != std::string::npos) {
        filename = inPath.substr(lastSlash + 1);
    }

    std::string prefixToReplace = "ConsumerResults_";
    size_t pos = filename.find(prefixToReplace);
    if (pos != std::string::npos) {
        filename.replace(pos, prefixToReplace.length(), "CumulativeProfiles_");
    } else {
        filename = "CumulativeProfiles_" + filename;
    }

    // Normalize the output directory to always end with exactly one trailing slash
    std::string outDir(outDirStr);
    if (outDir.empty() || outDir.back() != '/') outDir += "/";

    // Create the output directory (and any missing parents) if it does not exist yet.
    // kTRUE requests recursive creation, mirroring signalExtractionRing.cxx's convention.
    gSystem->mkdir(outDir.c_str(), kTRUE);

    std::string outFileStr = outDir + filename;

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

    // Check the input folder before anything is written. It is booked only with doFakePolDiagnosticsQA: without it this
    // config has nothing to accumulate, which is reported once and is not a failure (no output file is created).
    if (!fIn->Get(kInputDir.c_str())) {
        std::cout << " No " << kInputDir << " folder in this file (doFakePolDiagnosticsQA off in this config).\n"
                  << " Nothing to do: no output file written.\n";
        std::cout << "=======================================================\n\n";
        fIn->Close();
        delete fIn;
        return 0;
    }

    TFile* fOut = TFile::Open(outFileStr.c_str(), "RECREATE");
    // Checking output file before closing:
    if (!fOut || fOut->IsZombie()) {
        std::cerr << " ERROR: cannot create " << outFileStr << std::endl;
        if (fOut) { fOut->Close(); delete fOut; }
        fIn->Close();
        delete fIn;
        return 1;
    }

    // --- Create Organized Subfolders ---
    TDirectory* dirCounts   = fOut->mkdir("Cumulative_Counts");
    TDirectory* dir2D       = fOut->mkdir("Cumulative_2D");
    TDirectory* dir1D       = fOut->mkdir("Cumulative_1D");
    TDirectory* dirCanvases = fOut->mkdir("Comparison_Canvases");

    const std::string inDir = kInputDir + "/";

    // 1. Process TH2D (Counts)
    std::cout << " -> Processing TH2D Count Histograms...\n";
    for (const auto& dca : kDCASpecs) {
        const std::string name = "h2dFakePolCounts_PhiStarVs" + dca.tag;
        TH2D* hIn = (TH2D*)fIn->Get((inDir + name).c_str());
        if(hIn) {
            // Title built here rather than copied: the output axes are transposed (X = min DCA, Y = phi*)
            TH2D* hCumul = MakeCumulativeTH2D(hIn, CumulCountsName(dca).c_str(),
                ("Cumulative counts (" + kProxyFolder + ");" + dca.axisTitle + ";#phi^{*}").c_str());
            dirCounts->cd(); hCumul->Write();
            delete hCumul;
        } else {
             std::cerr << "    [Warning] Missing: " << name << std::endl;
        }
    }

    // 2. Process TProfile2D (Signals) -> Generates both 2D and 1D versions
    std::cout << " -> Processing TProfile2D Signal Histograms...\n";
    for (const auto& dca : kDCASpecs) {
        const std::string name = "p2dFakePolSignal_PhiStarVs" + dca.tag;
        TProfile2D* pIn = (TProfile2D*)fIn->Get((inDir + name).c_str());
        if(pIn) {
            TProfile2D* pCumul2D = MakeCumulativeProfile2D(pIn, Cumul2DName(dca).c_str(),
                ("Cumulative " + kRingLabel + ";" + dca.axisTitle + ";#phi^{*}").c_str());
            dir2D->cd(); pCumul2D->Write();

            TProfile* pCumul1D = MakeCumulativeProfile1D(pIn, Cumul1DName(dca).c_str(),
                ("Cumulative " + kRingLabel + ", #phi^{*}-integrated;" + dca.axisTitle + ";" + kRingLabel).c_str());
            dir1D->cd(); pCumul1D->Write();

            // Deleting to avoid leaks:
            delete pCumul2D; delete pCumul1D;
        } else {
             std::cerr << "    [Warning] Missing: " << name << std::endl;
        }
    }

    // 3. Process TProfile3D (Signals with Eta splits)
    std::cout << " -> Processing TProfile3D Eta-Split Histograms...\n";
    for (const auto& dca : kDCASpecs) {
        for (const auto& split : kEtaSplitSpecs) {
            const std::string name = "p3dFakePolSignal_PhiStarVs" + dca.tag + "Vs" + split.tag;
            TProfile3D* p3DIn = (TProfile3D*)fIn->Get((inDir + name).c_str());
            if(p3DIn) {
                ProcessProfile3D(p3DIn, dir2D, dir1D, dca, split);
            } else {
                 std::cerr << "    [Warning] Missing: " << name << std::endl;
            }
        }
    }

    // 4. Draw comparison canvases from the stored 1D cumulatives
    std::cout << " -> Drawing comparison canvases...\n";
    DrawComparisonCanvases(dir1D, dirCanvases);

    std::cout << " Successfully generated and saved all cumulative profiles.\n Thank you for waiting!\n";
    std::cout << "=======================================================\n\n";

    fOut->Write("", TObject::kOverwrite); // Just to make sure we wrote everything in memory, and overwriting whatever was already written
    fOut->Close();

    fIn->Close();

    delete fOut;
    delete fIn;

    return 0;
}