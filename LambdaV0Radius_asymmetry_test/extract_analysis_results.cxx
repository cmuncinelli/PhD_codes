// ----------------------------------------------------------------------
// processLambda.C
//
// Usage:
//   root -l -q 'processLambda.C("input.root", "output/")'
//
// This script:
//   1) Loads a ROOT file and retrieves the histogram:
//           asymmetric_rapidity_test/Lambda/hV0RadiusVsY (TH2D)
//
//   2) For EACH X-bin (V0Radius): project onto Y-axis -> TH1D in rapidity
//     Stored in: YProjections/
//
//   3) Each projection has Sumw2 enabled and is titled as:
//        "dN/dy, bin = N (V0Radius = Xlow-Xhigh)"
//
//   4) Then performs **Y-axis symmetric projection ratios**:
//        If Y has 50 bins (−1 to 1):
//           ratio = projection(26) / projection(25)
//           ratio = projection(27) / projection(24)
//           ...
//      with proper error propagation: sigma_R = R * sqrt( (sigmaA/A)^2 + (sigmaB/B)^2 )
//
//   5) Create MULTIPLOT canvas with all ratio histograms
//     Stored in: XProjectionRatios/RatiosMultiplot
//
//   6) Writes everything to a new ROOT file inside output folder.
//
// ----------------------------------------------------------------------

#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TDirectory.h>
#include <TString.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>
#include <filesystem>
#include <cmath>

using namespace std;

void extract_analysis_results(const char* inputPath, const char* outputFolder)
{
    // ------------------------------
    // Ensure output folder exists
    // ------------------------------
    std::filesystem::create_directories(outputFolder);

    TString outFileName = TString::Format("%sProcessedProjections.root", outputFolder); // Notice we already expect the folder to have the trailing "/"
    TFile* fout = new TFile(outFileName, "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "ERROR: Cannot create output ROOT file in folder: " << outputFolder << std::endl;
        return;
    }

    // ------------------------------
    // Open input ROOT file
    // ------------------------------
    TFile* fin = TFile::Open(inputPath, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "ERROR: Cannot open input file: " << inputPath << std::endl;
        return;
    }

    // ------------------------------
    // Retrieve histogram
    // ------------------------------
    TH2D* h2 = dynamic_cast<TH2D*>(
        fin->Get("asymmetric_rapidity_test/Lambda/hV0RadiusVsY")
    );

    if (!h2) {
        std::cerr << "ERROR: Histogram hV0RadiusVsY not found in file." << std::endl;
        return;
    }

    std::cout << "Histogram loaded: hV0RadiusVsY" << std::endl;

    // Write a copy of the original into output file
    fout->cd();
    h2->Write("hV0RadiusVsY_original");

    // Bin counts
    const int nBinsX = h2->GetNbinsX();
    const int nBinsY = h2->GetNbinsY();

    // Make directories inside output file
    fout->mkdir("YProjections");
    fout->mkdir("XProjectionRatios");

    // ==================================================================
    // (1) PROJECT ONTO Y FOR EACH X-BIN
    // ==================================================================
    fout->cd("YProjections");

    vector<TH1D*> yProjs;
    yProjs.reserve(nBinsX);

    for (int xbin = 1; xbin <= nBinsX; xbin++) {
        TString name = TString::Format("projY_xbin%d", xbin);

        TH1D* hproj = h2->ProjectionY(name, xbin, xbin);
        hproj->Sumw2();

        double xlow  = h2->GetXaxis()->GetBinLowEdge(xbin);
        double xhigh = h2->GetXaxis()->GetBinUpEdge(xbin);

        hproj->SetTitle(
            TString::Format("dN/dy, X-bin=%d (V0Radius = %.1f-%.1f)", 
                            xbin, xlow, xhigh)
        );

        hproj->GetXaxis()->SetTitle("y");
        hproj->GetYaxis()->SetTitle("Counts");

        hproj->Write();
        yProjs.push_back(hproj);
    }

    std::cout << "Created Y projections for each X bin." << std::endl;

    // -----------------------------------------------------------
    // (2) BUILD RATIO OF Y-SYMMETRIC PROJECTIONS
    // -----------------------------------------------------------

    // Y bins must be even
    if (nBinsY % 2 != 0) {
        std::cerr << "ERROR: Y-axis does not have an even number of bins!" << std::endl;
        return;
    }

    const int mid = nBinsY/2;   // example: 50 -> mid=25

    fout->cd("XProjectionRatios");
    vector<TH1D*> ratioHists;

    // Loop over symmetric Y bins:
    //   (mid+1)/(mid), (mid+2)/(mid-1), ...
    for (int i = 1; i <= mid; i++) {
        // Weird way of doing it, but notice that because i starts as 1 the Plus bins always have a higher index, as they should
        int yPlus  = mid + i;       // positive rapidity bin
        int yMinus = mid + 1 - i;   // symmetric negative rapidity bin

        TString nPlus  = TString::Format("temp_plus_%d",  yPlus);
        TString nMinus = TString::Format("temp_minus_%d", yMinus);

        // Project onto Y slices
        TH1D* A = h2->ProjectionX(nPlus, yPlus,  yPlus);
        TH1D* B = h2->ProjectionX(nMinus, yMinus, yMinus);

        double A_ylow  = h2->GetYaxis()->GetBinLowEdge(yPlus);
        double A_yhigh = h2->GetYaxis()->GetBinUpEdge(yPlus);
        double B_ylow  = h2->GetYaxis()->GetBinLowEdge(yMinus);
        double B_yhigh = h2->GetYaxis()->GetBinUpEdge(yMinus);

        A->Sumw2();
        B->Sumw2();

        TString rname = TString::Format("ratio_y%d_over_y%d", yPlus, yMinus);

        TH1D* R = (TH1D*)A->Clone(rname);
        R->Reset();
        R->SetTitle(
            TString::Format("ratio y=%.3f-%.3f over y = %.3f-%.3f (bin%d_over_y%d)", A_ylow, A_yhigh, B_ylow, B_yhigh, yPlus, yMinus)
        );

        // Adding some axis titles:
        R->GetXaxis()->SetTitle("#Lambda V0 Radius (cm)");
        R->GetYaxis()->SetTitle("yPlus / yMinus");
        // R->GetYaxis()->SetTitleOffset( 1.2);

        // Compute ratio with errors
        for (int by = 1; by <= nBinsY; by++) {
            double a  = A->GetBinContent(by);
            double b  = B->GetBinContent(by);
            double ea = A->GetBinError(by);
            double eb = B->GetBinError(by);

            if (b == 0) {
                R->SetBinContent(by, 0);
                R->SetBinError(by, 0);
                continue;
            }

            double ratio = a / b;
            double err = 0;

            if (a > 0 && b > 0) {
                err = ratio * std::sqrt( (ea*ea)/(a*a) + (eb*eb)/(b*b) );
            }

            R->SetBinContent(by, ratio);
            R->SetBinError(by, err);
        }

        R->Write();
        ratioHists.push_back(R);

        delete A;
        delete B;
    }

    std::cout << "Created symmetric Y-projection ratios." << std::endl;

    // ==================================================================
    // (3) MOSAIC MULTIPLOT CANVAS
    // ==================================================================
    int N = ratioHists.size();
    int cols = ceil(sqrt(N));
    int rows = ceil(double(N) / cols);

    fout->cd("XProjectionRatios");

    TCanvas* cMosaic = new TCanvas("RatiosMultiplot_Mosaic", "All Ratios - Mosaic", 1400, 1000);
    cMosaic->Divide(cols, rows);

    for (int i = 0; i < N; i++) {
        cMosaic->cd(i+1);
        ratioHists[i]->SetLineWidth(2);
        ratioHists[i]->SetStats(0);
        ratioHists[i]->GetYaxis()->SetTitle("Ratio");
        ratioHists[i]->Draw("E");

        TLegend* leg = new TLegend(0.1, 0.75, 0.9, 0.93);
        leg->SetHeader(ratioHists[i]->GetTitle(), "C");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetNColumns(2);
        leg->Draw();
    }

    cMosaic->Write();
    std::cout << "Created mosaic multiplot of all ratio histograms." << std::endl;

    // ==================================================================
    // (4) SINGLE OVERLAY MULTIPLOT (all ratios in one frame)
    // ==================================================================

    TCanvas* cOverlay = new TCanvas("RatiosMultiplot_Overlay", "All Ratios - Overlay", 1200, 800);

    double globalMax = 0.0;
    double globalMin = +1e9;

    int colorIndex = 2;

    for (auto* ratio : ratioHists) {
        for (int bx = 1; bx <= ratio->GetNbinsX(); bx++) {
            double val = ratio->GetBinContent(bx);
            if (val > globalMax) globalMax = val;
            if (val < globalMin && val > 0) globalMin = val;
        }
    }

    bool first = true;

    TLegend* L = new TLegend(0.15, 0.15, 0.85, 0.35);
    L->SetBorderSize(0);
    L->SetFillStyle(0);
    L->SetNColumns(4);

    // for (auto* ratio : ratioHists) {
    //     ratio->SetStats(0);
    //     ratio->SetLineWidth(2);
    //     ratio->SetLineColor(colorIndex);

    //     ratio->GetYaxis()->SetRangeUser(globalMin * 0.8, globalMax * 1.2);
    //     ratio->GetYaxis()->SetTitle("Ratio");

    //     if (first) {
    //         ratio->Draw("E");
    //         first = false;
    //     } else {
    //         ratio->Draw("E SAME");
    //     }

    //     L->AddEntry(ratio, ratio->GetTitle(), "l");
    //     colorIndex++;
    // }

    // Actually, rewriting to make it so we only plot 6, evenly spaced in y, curves to see the evolution with y:
    // int mid = ratioHists.size();       // or ratioHists.size()/2 depending on your logic
    int Nsel = 6;

    // Select 6 evenly spaced indices between 1 and mid
    std::vector<int> idxSel;
    idxSel.reserve(Nsel);
    for (int i = 0; i < Nsel; i++) {
        int idx = 1 + std::round(i * (mid - 1) / double(Nsel)); 
        if (idx >= mid) idx = mid - 1;
        idxSel.push_back(idx);
    }

    // Now plot only these
    // bool first = true;
    // int colorIndex = 2;

    for (int idx : idxSel) {

        TH1D* ratio = ratioHists[idx];
        if (!ratio) continue;

        ratio->SetStats(0);
        ratio->SetLineWidth(2);
        ratio->SetLineColor(colorIndex++);

        ratio->GetYaxis()->SetTitle("yPlus / yMinus");
        ratio->GetXaxis()->SetTitle("#Lambda V0 Radius (cm)");
        ratio->GetYaxis()->SetRangeUser(globalMin * 0.8, globalMax * 1.2);

        if (first) {
            ratio->Draw("E");
            first = false;
        } else {
            ratio->Draw("E SAME");
        }

        L->AddEntry(ratio, ratio->GetTitle(), "l");
    }

    L->Draw();
    cOverlay->Write();
    std::cout << "Created overlay multiplot of all ratio histograms." << std::endl;

    // -----------------------------------------------------------
    // Export everything
    // -----------------------------------------------------------

    fout->Close();
    fin->Close();

    std::cout << "\nAll results written to:\n  " << outFileName << "\n" << std::endl;
}
