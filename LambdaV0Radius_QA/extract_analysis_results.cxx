// ----------------------------------------------------------------------
// extract_analysis_results.cxx
//
// Usage:
//   root -l -q 'extract_analysis_results.cxx("input.root", "output/")'
// Example:
//   root -l -q 'extract_analysis_results.cxx("/home/cicero/Datasets/asymmetric_rapidity_tests_LHC23_pass4_Thin_small/results/AnalysisResults-FinerBins-1file.root", "/home/cicero/Datasets/asymmetric_rapidity_tests_LHC23_pass4_Thin_small/results/")'
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
#include <TH3D.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TF1.h>
#include <TDirectory.h>
#include <TString.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <iostream>
#include <filesystem>
#include <cmath>
#include <vector>

using namespace std;

void extract_analysis_results(const char* inputPath, const char* outputFolder)
{
    // Find where "AnalysisResults" occurs
    std::string key = "AnalysisResults";
    std::string inputPath_string = (std::string) inputPath;
    size_t pos = inputPath_string.find(key);

    std::string suffix = "";
    if (pos != std::string::npos) {
        pos += key.length(); // jump to after "AnalysisResults"
        size_t end = inputPath_string.rfind(".root");
        if (end != std::string::npos && end > pos) {
            suffix = inputPath_string.substr(pos, end - pos);  // e.g. "-hasTPCnoITS"
        }
    }

    std::string outFileName = (std::string) outputFolder + "ProcessedProjections" + suffix + ".root";
    TFile* fout = new TFile(outFileName.c_str(), "RECREATE");
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

        hproj->GetXaxis()->SetTitle("#Lambda y");
        hproj->GetYaxis()->SetTitle("#Lambda Counts");

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
        R->GetYaxis()->SetTitle("Counts yPlus / Counts yMinus");
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

        ratio->GetYaxis()->SetTitle("Counts yPlus / Counts yMinus");
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

    // ==================================================================
    // PART B: Invariant-mass filtering using hMass and hV0RadiusVsYVsMass
    // ==================================================================

    // --- (B1) Fit hMass with gaussian + linear to get mean and sigma ----
    TH1D* hMass = dynamic_cast<TH1D*>( fin->Get("asymmetric_rapidity_test/Lambda/hMass") );
    if (!hMass) {
        std::cerr << "ERROR: hMass not found in file; skipping InvMassFiltered workflow." << std::endl;
        fout->Close();
        fin->Close();
        return;
    }
    fout->cd();
    hMass->Sumw2();
    hMass->Write("hMass_original");

    // // Fit with gaussian + linear background:
    // // f(x) = G + (p0 + p1*x)
    // // Use a TF1 with 5 parameters: [0]*exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x
    // double massMin = hMass->GetXaxis()->GetXmin();
    // double massMax = hMass->GetXaxis()->GetXmax();
    
    // TF1* fitFunc = new TF1("fitFunc", "[0]*exp(-0.5*((x-[1])/[2])**2) + [3] + [4]*x", massMin, massMax);
    
    // Actually best to fit with Gaussian + quadratic, after looking at the fit:
    // f(x) = [0]*exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x + [5]*x^2
    double massMin = hMass->GetXaxis()->GetXmin();
    double massMax = hMass->GetXaxis()->GetXmax();
    TF1* fitFunc = new TF1("fitFunc", "[0]*exp(-0.5*((x-[1])/[2])**2) + [3] + [4]*x + [5]*x*x", massMin, massMax);

    // Reasonable initial guesses: peak height, mean at bin center of max bin, sigma ~ bin width * few
    // [0] Amp, [1] Mean, [2] Sigma, [3] BG_const, [4] BG_slope, [5] BG_quad
    int maxBin = hMass->GetMaximumBin(); // This is a simple bin, so int should work
    double peak = hMass->GetBinContent(maxBin);
    double mean0 = hMass->GetBinCenter(maxBin);
    double sigma0 = hMass->GetRMS();
    if (sigma0 <= 0) sigma0 = (massMax - massMin)/100.0;
    hMass->Scale(1./peak);
    TH1D* hMassScaled = (TH1D*)hMass->Clone("hMassScaled");
    hMassScaled->Write(); // For comparison with the fit.

    // Constraints on parameters:
        // With too much data, fitting hMass without normalizing is numerically unstable.
        // Have to stabilize this as much as possible.
    // fitFunc->SetParLimits(0, 0.0, 2*peak); // Amplitude should be positive, and should not be nearly as high as the peak
    fitFunc->SetParLimits(0, 0.0, 2); // New value for the peak is 1. by definition of the rescaling.
    fitFunc->SetParLimits(1, massMin, massMax); // Mean is within the minimum and maximum of mass
    fitFunc->SetParLimits(2, 0.0, 1); // Sigma should always be positive
    
    fitFunc->SetParameters(peak, mean0, sigma0, 0.1, 0.1, 0.1);
    fitFunc->SetParNames("Amp","Mean","Sigma","BG_const","BG_slope","BG_quad");
    // fitFunc->SetParameters(peak, mean0, sigma0, 0.1, 0.1);
    // fitFunc->SetParNames("Amp","Mean","Sigma","BG_const","BG_slope");

    // Fit in a window around the peak to avoid far tails, say +/- 5*sigma0 initially (conservative)
    double fit_low = mean0 - 5.0 * sigma0;
    double fit_high = mean0 + 5.0 * sigma0;
    if (fit_low < massMin) fit_low = massMin;
    if (fit_high > massMax) fit_high = massMax;

    // Perform fit (quiet)
        // Could use RQL for log-likelihood fits too.
    int fitStatus = hMass->Fit(fitFunc, "RQ", "", fit_low, fit_high);

    double mean = fitFunc->GetParameter(1);
    double sigma = fabs(fitFunc->GetParameter(2));
    if (sigma <= 0) sigma = sigma0;

    // ---- Print full fit summary ----
    std::cout << "\n\n================ Mass Fit Results (Gaussian + Quadratic) ================\n";

    // Print mean and sigma explicitly
    std::cout << "Peak position (mean)        = " << mean  << " GeV/c^2\n";
    std::cout << "Gaussian sigma              = " << sigma << " GeV/c^2\n";

    // Fit status
    std::cout << "Fit status code             = " << fitStatus << "\n";

    // Parameter table
    std::cout << "\nFit parameters:\n";
    int npar = fitFunc->GetNpar();
    for (int i = 0; i < npar; i++) {
        std::cout << "  " << std::setw(12) << fitFunc->GetParName(i)
                << " = " << std::setw(12) << fitFunc->GetParameter(i)
                << " ± " << fitFunc->GetParError(i) << "\n";
    }

    // Goodness of fit
    double chi2 = fitFunc->GetChisquare();
    int ndf = fitFunc->GetNDF();
    double chi2_ndf = (ndf > 0 ? chi2 / ndf : -1);

    std::cout << "\nChi-square                  = " << chi2 << "\n";
    std::cout << "NDF                         = " << ndf  << "\n";
    std::cout << "Reduced chi-square          = " << chi2_ndf << "\n";
    // Gives off a humongous chi2 because the error bars are waaaaay to small
    std::cout << "=======================================================================\n\n";

    // Save fitted function
    fout->cd();
    // fitFunc->SetName("hMass_gausPlusLinear");
    fitFunc->SetName("hMass_gausPlusQuadratic");
    fitFunc->Write();

    // --- Define mass window: mean +/- 2.5*sigma
    double mLow = mean - 2.5 * sigma;
    double mHigh = mean + 2.5 * sigma;
    std::cout << "Using invariant mass window: [" << mLow << " , " << mHigh << "]" << std::endl;

    // --- (B2) Slice TH3D in mass axis to obtain TH2D (mass filtered)
    TH3D* h3 = dynamic_cast<TH3D*>( fin->Get("asymmetric_rapidity_test/Lambda/hV0RadiusVsYVsMass") );
    if (!h3) {
        std::cerr << "ERROR: TH3D hV0RadiusVsYVsMass not found; skipping InvMassFiltered workflow." << std::endl;
        fout->Close();
        fin->Close();
        return;
    }

    // Find z-bins corresponding to mLow, mHigh
    int zbinLow  = h3->GetZaxis()->FindBin(mLow);
    int zbinHigh = h3->GetZaxis()->FindBin(mHigh);

    // Ensure bin ordering
    if (zbinLow > zbinHigh) std::swap(zbinLow, zbinHigh);

    // --- Recommended method in ROOT 6.36 ---
    // Restrict Z range to the desired mass window -- Then we can integrate over it to get the TH2D in the desired mass interval
    // ROOT 6.38 has a fancier function (below), but this is not available for the current O2env's ROOT 6.36:
    // TH2D* h2 = h3->ProjectionXY("name", zbinLow, zbinHigh);
    h3->GetZaxis()->SetRange(zbinLow, zbinHigh);

    // Projection: produce TH2D over X-Y axes integrating over Z in [zbinLow, zbinHigh]
    // Project onto X-Y plane
    TH2D* h2_massFiltered = (TH2D*)h3->Project3D("yx"); // Had to project as yx to make the dimensions the same as the
                                                        // previous TH2D (V0Radius on X, rapidity on Y).
                                                        // Be careful with this!!!
    h2_massFiltered->SetName("hV0RadiusVsY_InvMassFiltered");

    // Reset Z-axis range back to full histogram
    h3->GetZaxis()->SetRange(0, 0);

    if (!h2_massFiltered) {
        std::cerr << "ERROR: Project3D(xy) failed; skipping InvMassFiltered workflow." << std::endl;
        fout->Close();
        fin->Close();
        return;
    }
    h2_massFiltered->Sumw2();
    fout->cd();
    h2_massFiltered->Write("hV0RadiusVsY_InvMassFiltered_original");

    std::cout << "Created mass-filtered TH2D: hV0RadiusVsY_InvMassFiltered" << std::endl;

    // ============================
    // Now reproduce the entire workflow for the mass-filtered TH2D
    // ============================

    const int nBinsX_f = h2_massFiltered->GetNbinsX();
    const int nBinsY_f = h2_massFiltered->GetNbinsY();

    // Create directories for the filtered workflow
    fout->mkdir("YProjections_InvMassFiltered");
    fout->mkdir("XProjectionRatios_InvMassFiltered");

    // ==================================================================
    // (1) PROJECT ONTO Y FOR EACH X-BIN (mass-filtered)
    // ==================================================================
    fout->cd("YProjections_InvMassFiltered");

    vector<TH1D*> yProjs_f;
    yProjs_f.reserve(nBinsX_f);

    for (int xbin = 1; xbin <= nBinsX_f; xbin++) {
        TString name = TString::Format("InvProjY_xbin%d", xbin);

        TH1D* hproj = h2_massFiltered->ProjectionY(name, xbin, xbin);
        // hproj->Sumw2();
        if (!hproj->GetSumw2N()) hproj->Sumw2(); // Checks if there is an error storage already attributed to this projection
                                                 // (ProjectionY would create it, as h2 is already a projection with Sumw2()
                                                 // activated!)

        double xlow  = h2_massFiltered->GetXaxis()->GetBinLowEdge(xbin);
        double xhigh = h2_massFiltered->GetXaxis()->GetBinUpEdge(xbin);

        hproj->SetTitle( TString::Format("InvMass filtered dN/dy, X-bin=%d (V0Radius = %.3f-%.3f)", xbin, xlow, xhigh) );
        hproj->GetXaxis()->SetTitle("#Lambda y");
        hproj->GetYaxis()->SetTitle("#Lambda Counts");

        hproj->Write();
        yProjs_f.push_back(hproj);
    }

    std::cout << "Created Y projections for each X bin (InvMassFiltered TH2D)." << std::endl;

    // -----------------------------------------------------------
    // (2) BUILD RATIO OF Y-SYMMETRIC PROJECTIONS (InvMassFiltered)
    // -----------------------------------------------------------
    if (nBinsY_f % 2 != 0) {
        std::cerr << "ERROR: Filtered Y-axis does not have an even number of bins!" << std::endl;
        // proceed but warn (or return). Here we will return to avoid mis-indexing
        fout->Close();
        fin->Close();
        return;
    }

    const int mid_f = nBinsY_f / 2;
    fout->cd("XProjectionRatios_InvMassFiltered");
    vector<TH1D*> ratioHists_f;
    ratioHists_f.reserve(mid_f);

    for (int i = 1; i <= mid_f; i++) {
        int yPlus  = mid_f + i;
        int yMinus = mid_f + 1 - i;

        TString namePlus  = TString::Format("Inv_plus_%d",  yPlus);
        TString nameMinus = TString::Format("Inv_minus_%d", yMinus);

        TH1D* A = h2_massFiltered->ProjectionX(namePlus,  yPlus,  yPlus);
        TH1D* B = h2_massFiltered->ProjectionX(nameMinus, yMinus, yMinus);

        if (A->GetSumw2N() == 0) A->Sumw2(); // Same check as earlier, before activating Sumw2()
        if (B->GetSumw2N() == 0) B->Sumw2();

        double A_ylow  = h2_massFiltered->GetYaxis()->GetBinLowEdge(yPlus);
        double A_yhigh = h2_massFiltered->GetYaxis()->GetBinUpEdge(yPlus);
        double B_ylow  = h2_massFiltered->GetYaxis()->GetBinLowEdge(yMinus);
        double B_yhigh = h2_massFiltered->GetYaxis()->GetBinUpEdge(yMinus);

        TString rname = TString::Format("Inv_ratio_y%d_over_y%d", yPlus, yMinus);

        TH1D* R = (TH1D*)A->Clone(rname);
        R->Reset();
        R->SetTitle( TString::Format("InvMass ratio y=%.3f-%.3f over y=%.3f-%.3f (bin%d_over_y%d)", A_ylow, A_yhigh, B_ylow, B_yhigh, yPlus, yMinus) );

        // Axis titles
        R->GetXaxis()->SetTitle("#Lambda V0 Radius (cm)");
        R->GetYaxis()->SetTitle("Counts yPlus / Counts yMinus");

        // Compute ratio with errors
        for (int bx = 1; bx <= R->GetNbinsX(); bx++) {
            double a  = A->GetBinContent(bx);
            double b  = B->GetBinContent(bx);
            double ea = A->GetBinError(bx);
            double eb = B->GetBinError(bx);

            if (b == 0) {
                R->SetBinContent(bx, 0);
                R->SetBinError(bx, 0);
                continue;
            }

            double ratio = a / b;
            double err = 0;
            if (a > 0 && b > 0) {
                err = ratio * std::sqrt( (ea*ea)/(a*a) + (eb*eb)/(b*b) );
            }

            R->SetBinContent(bx, ratio);
            R->SetBinError(bx, err);
        }

        R->Write();
        ratioHists_f.push_back(R);

        delete A;
        delete B;
    }

    std::cout << "Created symmetric Y-projection ratios (InvMassFiltered)." << std::endl;

    // ==================================================================
    // (3) MOSAIC MULTIPLOT CANVAS (InvMassFiltered)
    // ==================================================================
    {
        int Nf = ratioHists_f.size();
        int cols = ceil(sqrt(Nf));
        int rows = ceil(double(Nf) / cols);

        fout->cd("XProjectionRatios_InvMassFiltered");
        TCanvas* cMosaicF = new TCanvas("RatiosMultiplot_Mosaic_InvMassFiltered", "All Ratios - Mosaic (InvMassFiltered)", 1400, 1000);
        cMosaicF->Divide(cols, rows);

        for (int i = 0; i < Nf; i++) {
            cMosaicF->cd(i+1);
            ratioHists_f[i]->SetLineWidth(2);
            ratioHists_f[i]->SetStats(0);
            ratioHists_f[i]->GetYaxis()->SetTitle("Ratio");
            ratioHists_f[i]->Draw("E");

            TLegend* leg = new TLegend(0.1, 0.75, 0.9, 0.93);
            leg->SetHeader(ratioHists_f[i]->GetTitle(), "C");
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetNColumns(2);
            leg->Draw();
        }

        cMosaicF->Write();
        std::cout << "Created mosaic multiplot of all ratio histograms (InvMassFiltered)." << std::endl;
    }

    // ==================================================================
    // (4) SINGLE OVERLAY MULTIPLOT (InvMassFiltered) - select 6 indices
    // ==================================================================
    {
        TCanvas* cOverlayF = new TCanvas("RatiosMultiplot_Overlay_InvMassFiltered", "All Ratios - Overlay (InvMassFiltered)", 1200, 800);

        double globalMax = 0.0;
        double globalMin = +1e9;
        for (auto* ratio : ratioHists_f) {
            for (int bx = 1; bx <= ratio->GetNbinsX(); bx++) {
                double val = ratio->GetBinContent(bx);
                if (val > globalMax) globalMax = val;
                if (val < globalMin && val > 0) globalMin = val;
            }
        }

        int Nsel = 6;
        std::vector<int> idxSel;
        idxSel.reserve(Nsel);
        for (int i = 0; i < Nsel; i++) {
            int idx = 1 + std::round(i * (mid_f - 1) / double(Nsel));
            if (idx >= mid_f) idx = mid_f - 1;
            idxSel.push_back(idx);
        }

        bool first = true;
        int colorIndex = 2;
        TLegend* Lf = new TLegend(0.15, 0.15, 0.85, 0.35);
        Lf->SetBorderSize(0);
        Lf->SetFillStyle(0);
        Lf->SetNColumns(3);

        for (int idx : idxSel) {
            if (idx < 0 || idx >= (int)ratioHists_f.size()) continue;
            TH1D* ratio = ratioHists_f[idx];
            if (!ratio) continue;

            ratio->SetStats(0);
            ratio->SetLineWidth(2);
            ratio->SetLineColor(colorIndex++);

            ratio->GetYaxis()->SetTitle("Counts yPlus / Counts yMinus");
            ratio->GetXaxis()->SetTitle("#Lambda V0 Radius (cm)");
            if (globalMin < globalMax) ratio->GetYaxis()->SetRangeUser(globalMin * 0.8, globalMax * 1.2);

            if (first) {
                ratio->Draw("E");
                first = false;
            } else {
                ratio->Draw("E SAME");
            }

            Lf->AddEntry(ratio, ratio->GetTitle(), "l");
        }

        Lf->Draw();
        cOverlayF->Write();
        std::cout << "Created overlay multiplot of selected ratio histograms (InvMassFiltered)." << std::endl;
    }

    // ======================================================================
    // PART C: Coarse Y-bin ratios (InvMassFiltered)
    //         (8 coarse bins of width Δy = 0.25 -> 4 symmetric ratios)
    //         (just assumed limits go from -1 to 1 and did some things
    //          in a sloppier, faster way, just to compare quickly)
    // ======================================================================
    std::cout << "\n=== Building COARSE Y-binning ratios (InvMassFiltered) ===\n";

    fout->mkdir("XProjectionRatios_InvMassFiltered_CoarseY");
    fout->cd("XProjectionRatios_InvMassFiltered_CoarseY");

    // ---------------------------------------------------------
    // Determine coarse bins
    // ---------------------------------------------------------
    double yMin  = h2_massFiltered->GetYaxis()->GetXmin();   // should be -1
    double yMax  = h2_massFiltered->GetYaxis()->GetXmax();   // should be +1. Not hard coded, but whatever
    double width = (yMax - yMin) / 8.0;                      // 0.25

    // Coarse bin boundaries
    std::vector<double> coarseEdges;
    coarseEdges.reserve(9);

    for (int i = 0; i <= 8; i++)
        coarseEdges.push_back(yMin + i * width);

    // Helper function to get fine-bin ranges for a given coarse interval
    auto getFineBinRange = [&](double low, double high) {
        int b1 = h2_massFiltered->GetYaxis()->FindBin(low  + 1e-6);
        int b2 = h2_massFiltered->GetYaxis()->FindBin(high - 1e-6);
        if (b1 < 1) b1 = 1;
        if (b2 > h2_massFiltered->GetNbinsY()) b2 = h2_massFiltered->GetNbinsY();
        return std::make_pair(b1, b2);
    };

    // ---------------------------------------------------------
    // Build the 4 symmetric ratios
    // ---------------------------------------------------------
    std::vector<TH1D*> coarseRatios;

    for (int i = 0; i < 4; i++) {

        double pos_low  = coarseEdges[4 + i];
        double pos_high = coarseEdges[5 + i];
        double neg_low  = coarseEdges[3 - i];
        double neg_high = coarseEdges[4 - i];

        auto [ypL, ypH] = getFineBinRange(pos_low, pos_high);  // positive slice
        auto [ymL, ymH] = getFineBinRange(neg_low, neg_high);  // negative slice

        TString nameA = TString::Format("CoarsePos_%d", i);
        TString nameB = TString::Format("CoarseNeg_%d", i);

        // Integrate over many Y bins --> produce a TH1D for X
        TH1D* A = h2_massFiltered->ProjectionX(nameA, ypL, ypH);
        TH1D* B = h2_massFiltered->ProjectionX(nameB, ymL, ymH);

        if (A->GetSumw2N() == 0) A->Sumw2();
        if (B->GetSumw2N() == 0) B->Sumw2();

        TString rname = TString::Format("CoarseRatio_%d", i);
        TH1D* R = (TH1D*)A->Clone(rname);
        R->Reset();

        R->SetTitle(
            TString::Format("Ratio y=[%.3f,%.3f] / y=[%.3f,%.3f]",
                            pos_low,pos_high, neg_low,neg_high)
        );

        R->GetXaxis()->SetTitle("#Lambda V0 Radius (cm)");
        R->GetYaxis()->SetTitle("Counts yPlus / Counts yMinus");

        // Compute ratio
        int nXB = A->GetNbinsX();
        for (int bx = 1; bx <= nXB; bx++) {

            double a  = A->GetBinContent(bx);
            double b  = B->GetBinContent(bx);
            double ea = A->GetBinError(bx);
            double eb = B->GetBinError(bx);

            if (b == 0) {
                R->SetBinContent(bx, 0);
                R->SetBinError(bx, 0);
                continue;
            }

            double ratio = a / b;
            double err = 0;
            if (a > 0 && b > 0)
                err = ratio * sqrt( (ea*ea)/(a*a) + (eb*eb)/(b*b) );

            R->SetBinContent(bx, ratio);
            R->SetBinError(bx, err);
        }

        R->Write();
        coarseRatios.push_back(R);

        delete A;
        delete B;
    }

    std::cout << "Created 4 coarse Y-bin symmetric ratios.\n";

    // ---------------------------------------------------------
    // Overlay plot
    // ---------------------------------------------------------
    {
        TCanvas* cO = new TCanvas("CoarseRatios_Overlay", "Coarse Y ratios — Overlay", 1100, 800);

        double gmax = 0, gmin = 1e9;
        for (auto* R : coarseRatios) {
            for (int b=1; b<=R->GetNbinsX(); b++) {
                double v = R->GetBinContent(b);
                if (v > gmax) gmax = v;
                if (v > 0 && v < gmin) gmin = v;
            }
        }

        bool first = true;
        int ci = 2;
        TLegend* L = new TLegend(0.15, 0.15, 0.85, 0.35);
        L->SetBorderSize(0);
        L->SetFillStyle(0);
        L->SetNColumns(2);

        for (int i=0; i<4; i++) {
            TH1D* R = coarseRatios[i];
            R->SetLineColor(ci++);
            R->SetLineWidth(3);
            R->SetStats(0);
            R->GetYaxis()->SetRangeUser(gmin*0.8, gmax*1.2);

            if (first) { R->Draw("E"); first = false; }
            else       { R->Draw("E SAME"); }

            L->AddEntry(R, R->GetTitle(), "l");
        }

        L->Draw();
        cO->Write();
        std::cout << "Created coarse-ratio overlay plot.\n";
    }

    // ======================================================================
    // PART C extension: Build NEW normalized coarse ratios
    //         RECOMPUTE A(x) and B(x) using the original TH2D
    // ======================================================================

    std::cout << "Building normalized coarse ratios from scratch...\n";

    // 1. Compute total positive and total negative rapidity yields
    double totalPos = 0.0;
    double totalNeg = 0.0;

    int yMid = h2_massFiltered->GetYaxis()->FindBin(0.0);

    for (int i = 1; i <= mid_f; i++) {
        int yPlus  = mid_f + i;
        int yMinus = mid_f + 1 - i;

        TString namePlus  = TString::Format("Inv_plus_%d",  yPlus);
        TString nameMinus = TString::Format("Inv_minus_%d", yMinus);

        TH1D* projPos = h2_massFiltered->ProjectionX(namePlus,  yPlus,  yPlus);
        TH1D* projNeg = h2_massFiltered->ProjectionX(nameMinus, yMinus, yMinus);
        totalPos += projPos->GetEntries();
        totalNeg += projNeg->GetEntries();
    }

    std::cout << "  Total +y Λ = " << totalPos << "\n"
            << "  Total -y Λ = " << totalNeg << "\n";

    if (totalPos <= 0 || totalNeg <= 0) {
        std::cerr << "ERROR: Zero totalPos or totalNeg; cannot normalize.\n";
    } else {

        fout->cd("XProjectionRatios_InvMassFiltered_CoarseY");

        std::vector<TH1D*> coarseRatiosNorm;

        // Loop over each coarse pair
        for (int i = 0; i < 4; i++) {

            // Determine coarse Y ranges
            double pos_low  = coarseEdges[4 + i];
            double pos_high = coarseEdges[5 + i];
            double neg_low  = coarseEdges[3 - i];
            double neg_high = coarseEdges[4 - i];

            auto [ypL, ypH] = getFineBinRange(pos_low,  pos_high);
            auto [ymL, ymH] = getFineBinRange(neg_low,  neg_high);

            // Recompute A(x) and B(x) from the original TH2D
            TH1D* Araw = h2_massFiltered->ProjectionX(
                TString::Format("CoarseNorm_PosRaw_%d", i), ypL, ypH);

            TH1D* Braw = h2_massFiltered->ProjectionX(
                TString::Format("CoarseNorm_NegRaw_%d", i), ymL, ymH);

            Araw->Sumw2();
            Braw->Sumw2();

            int nXB = Araw->GetNbinsX();

            // Create normalized ratio histogram
            TH1D* R = new TH1D(
                TString::Format("CoarseRatio_%d_Norm", i),
                TString::Format("Normalized ratio y=[%.2f,%.2f] / y=[%.2f,%.2f]",
                                pos_low,pos_high, neg_low,neg_high),
                nXB,
                Araw->GetXaxis()->GetXmin(),
                Araw->GetXaxis()->GetXmax()
            );

            R->Sumw2();
            R->GetXaxis()->SetTitle("#Lambda V0 Radius (cm)");
            R->GetYaxis()->SetTitle("(N_pos(r,y)/int N_pos(r,y) dr) / (N_neg(r,y)/int N_neg(r,y) dr)");

            // Fill R(x) = [A(x)/totalPos] / [B(x)/totalNeg]
            for (int bx = 1; bx <= nXB; bx++) {

                double A = Araw->GetBinContent(bx);
                double B = Braw->GetBinContent(bx);
                double eA = Araw->GetBinError(bx);
                double eB = Braw->GetBinError(bx);

                if (A < 0) A = 0;
                if (B < 0) B = 0;

                if (B == 0) {
                    R->SetBinContent(bx, 0);
                    R->SetBinError(bx, 0);
                    continue;
                }

                double A_norm = A / totalPos;
                double B_norm = B / totalNeg;

                double ratio = A_norm / B_norm;

                // Propagate uncertainty correctly
                double err = 0;
                if (A > 0 && B > 0) {
                    err = ratio * std::sqrt(
                        (eA*eA)/(A*A) + (eB*eB)/(B*B)
                    );
                }

                R->SetBinContent(bx, ratio);
                R->SetBinError(bx, err);
            }

            R->Write();
            coarseRatiosNorm.push_back(R);

            delete Araw;
            delete Braw;
        }

        std::cout << "Normalized coarse ratios successfully created.\n";
    }


    std::cout << "=== Finished Part C (coarse Y ratios). ===\n\n";


    // ==================================================================
    // PART D: Repeat workflow, but for Z vertex position instead of y
    // ==================================================================

    std::cout << "\n\n========== PART D: Z-vertex projection analysis ==========\n";

    // Load TH3D for Z coordinate
    TH3D* h3_ZPos = dynamic_cast<TH3D*>( fin->Get("asymmetric_rapidity_test/Lambda/hV0RadiusVsZVsMass") );
    if (!h3_ZPos) {
        std::cerr << "ERROR: TH3D hV0RadiusVsZVsMass not found! Skipping Part D." << std::endl;
    } 
    // Convert the mass window (mLow, mHigh) to Z-bins
    int zbinLow2  = h3_ZPos->GetZaxis()->FindBin(mLow);
    int zbinHigh2 = h3_ZPos->GetZaxis()->FindBin(mHigh);
    if (zbinLow2 > zbinHigh2) std::swap(zbinLow2, zbinHigh2);

    // Restrict Z range for projection
    h3_ZPos->GetZaxis()->SetRange(zbinLow2, zbinHigh2);

    // Project 3D --> 2D (X = radius, Y = Z-vertex)
    TH2D* h2_ZPos = (TH2D*) h3_ZPos->Project3D("yx");
    h2_ZPos->SetName("hV0RadiusVsZ_InvMassFiltered");

    // Restore original range
    h3_ZPos->GetZaxis()->SetRange(0, 0);

    // Write original filtered 2D
    fout->cd();
    h2_ZPos->Sumw2();
    h2_ZPos->SetOption("COLZ");
    h2_ZPos->Write("hV0RadiusVsZ_InvMassFiltered_original");

    std::cout << "Created mass-filtered TH2D: hV0RadiusVsZ_InvMassFiltered\n";


    // ============================
    // Now perform the Z-based workflow
    // ============================

    const int nBinsX_Z = h2_ZPos->GetNbinsX();
    const int nBinsZ   = h2_ZPos->GetNbinsY();   // IMPORTANT: second coordinate is Z, not y

    // Directories for Z projections
    fout->mkdir("ZProjections_InvMassFiltered");
    fout->mkdir("XProjectionRatios_Z_InvMassFiltered");

    // ==================================================================
    // (1) PROJECT ONTO Z FOR EACH X-BIN
    // ==================================================================

    fout->cd("ZProjections_InvMassFiltered");

    std::vector<TH1D*> ZProjs;
    ZProjs.reserve(nBinsX_Z);

    for (int xbin = 1; xbin <= nBinsX_Z; xbin++) {

        TString name = TString::Format("InvProjZ_xbin%d", xbin);

        TH1D* hproj = h2_ZPos->ProjectionY(name, xbin, xbin);

        if (!hproj->GetSumw2N()) hproj->Sumw2();

        double xlow  = h2_ZPos->GetXaxis()->GetBinLowEdge(xbin);
        double xhigh = h2_ZPos->GetXaxis()->GetBinUpEdge(xbin);

        hproj->SetTitle(
            TString::Format("InvMass filtered Z counts, X-bin=%d (Radius=%.3f-%.3f)",
                            xbin, xlow, xhigh)
        );

        hproj->GetXaxis()->SetTitle("V0 Z vertex (cm)");
        hproj->GetYaxis()->SetTitle("#Lambda Counts");

        hproj->Write();
        ZProjs.push_back(hproj);
    }

    std::cout << "Created Z projections for each X bin.\n";


    // ==================================================================
    // (2) BUILD RATIO OF Z-SYMMETRIC PROJECTIONS
    // ==================================================================

    if (nBinsZ % 2 != 0) {
        std::cerr << "ERROR: Z-axis does not have an even number of bins! Aborting Z-ratios.\n";
    } else {

        const int midZ = nBinsZ / 2;

        fout->cd("XProjectionRatios_Z_InvMassFiltered");

        std::vector<TH1D*> ZRatios;
        ZRatios.reserve(midZ);

        for (int i = 1; i <= midZ; i++) {

            int zPlus  = midZ + i;
            int zMinus = midZ + 1 - i;

            TString namePlus  = TString::Format("InvZ_plus_%d",  zPlus);
            TString nameMinus = TString::Format("InvZ_minus_%d", zMinus);

            TH1D* A = h2_ZPos->ProjectionX(namePlus,  zPlus,  zPlus);
            TH1D* B = h2_ZPos->ProjectionX(nameMinus, zMinus, zMinus);

            if (A->GetSumw2N() == 0) A->Sumw2();
            if (B->GetSumw2N() == 0) B->Sumw2();

            double A_zlow  = h2_ZPos->GetYaxis()->GetBinLowEdge(zPlus);
            double A_zhigh = h2_ZPos->GetYaxis()->GetBinUpEdge(zPlus);
            double B_zlow  = h2_ZPos->GetYaxis()->GetBinLowEdge(zMinus);
            double B_zhigh = h2_ZPos->GetYaxis()->GetBinUpEdge(zMinus);

            TString rname =
                TString::Format("InvZ_ratio_z%d_over_z%d", zPlus, zMinus);

            TH1D* R = (TH1D*) A->Clone(rname);
            R->Reset();
            R->SetTitle(
                TString::Format("InvMass Z-ratio Z=[%.2f,%.2f] / Z=[%.2f,%.2f]",
                                A_zlow,A_zhigh, B_zlow,B_zhigh)
            );

            R->GetXaxis()->SetTitle("#Lambda V0 Radius (cm)");
            R->GetYaxis()->SetTitle("Counts zPlus / Counts zMinus");

            // Compute ratio with errors
            for (int bx = 1; bx <= R->GetNbinsX(); bx++) {

                double a  = A->GetBinContent(bx);
                double b  = B->GetBinContent(bx);
                double ea = A->GetBinError(bx);
                double eb = B->GetBinError(bx);

                if (b <= 0) {
                    R->SetBinContent(bx, 0);
                    R->SetBinError(bx, 0);
                    continue;
                }

                double ratio = a / b;
                double err   = 0;

                if (a > 0 && b > 0) {
                    err = ratio * std::sqrt(
                        (ea*ea)/(a*a) + (eb*eb)/(b*b)
                    );
                }

                R->SetBinContent(bx, ratio);
                R->SetBinError(bx, err);
            }

            R->Write();
            ZRatios.push_back(R);

            delete A;
            delete B;
        }

        std::cout << "Created symmetric Z-projection ratios (InvMassFiltered).\n";
    }

    // ======================================================================
    // PART D-extra: Additional coarse-bin ratios, differences, and full 2D asymmetries
    // ======================================================================

    std::cout << "\n\n========== PART D-extra: Additional Z-vertex asymmetry analyses ==========\n";

    // Safety check
    if (!h2_ZPos) {
        std::cerr << "ERROR: h2_ZPos not available — cannot run Part D-extra.\n";
    } else {

        const int nBinsZ_d = h2_ZPos->GetNbinsY();
        if (nBinsZ_d % 2 != 0) {
            std::cerr << "ERROR: Z-axis has an odd number of bins — cannot build symmetric structures.\n";
        } else {

            const int midZ_d = nBinsZ_d / 2;

            // Build coarse binning identical to Part E
            const int nCoarseZ = 8;  // → 4 symmetric pairs
            double coarseEdgesZ[nCoarseZ + 1];

            double zMin = h2_ZPos->GetYaxis()->GetXmin();
            double zMax = h2_ZPos->GetYaxis()->GetXmax();
            double zStep = (zMax - zMin) / nCoarseZ;

            for (int i = 0; i <= nCoarseZ; i++)
                coarseEdgesZ[i] = zMin + i * zStep;

            auto getZBins = [&](double z1, double z2) {
                int b1 = h2_ZPos->GetYaxis()->FindBin(z1);
                int b2 = h2_ZPos->GetYaxis()->FindBin(z2);
                if (b1 > b2) std::swap(b1, b2);
                return std::pair<int,int>(b1, b2);
            };

            // ------------------------------------------------------------------
            // (1) RAW coarse ratios (A/B) — NO normalization
            // ------------------------------------------------------------------

            fout->mkdir("XProjectionRatios_Z_InvMassFiltered_CoarseZ_raw");
            fout->cd("XProjectionRatios_Z_InvMassFiltered_CoarseZ_raw");

            for (int i = 0; i < 4; i++) {

                double pos_low  = coarseEdgesZ[4 + i];
                double pos_high = coarseEdgesZ[5 + i];
                double neg_low  = coarseEdgesZ[3 - i];
                double neg_high = coarseEdgesZ[4 - i];

                auto [zpL, zpH] = getZBins(pos_low, pos_high);
                auto [zmL, zmH] = getZBins(neg_low, neg_high);

                TH1D* A = h2_ZPos->ProjectionX(
                    TString::Format("CoarseZ_PosRawNoNorm_%d", i), zpL, zpH);
                TH1D* B = h2_ZPos->ProjectionX(
                    TString::Format("CoarseZ_NegRawNoNorm_%d", i), zmL, zmH);

                A->Sumw2();
                B->Sumw2();

                int nXB = A->GetNbinsX();
                TH1D* R = new TH1D(
                    TString::Format("CoarseZRatio_%d_Raw", i),
                    TString::Format("Z ratio Z=[%.1f,%.1f] / Z=[%.1f,%.1f]",
                                    pos_low, pos_high, neg_low, neg_high),
                    nXB, A->GetXaxis()->GetXmin(), A->GetXaxis()->GetXmax()
                );

                R->Sumw2();
                R->GetXaxis()->SetTitle("#Lambda V0 Radius (cm)");
                R->GetYaxis()->SetTitle("Counts zPlus / Counts zMinus");

                for (int bx = 1; bx <= nXB; bx++) {

                    double a  = A->GetBinContent(bx);
                    double b  = B->GetBinContent(bx);
                    double ea = A->GetBinError(bx);
                    double eb = B->GetBinError(bx);

                    if (b == 0) {
                        R->SetBinContent(bx, 0);
                        R->SetBinError(bx, 0);
                        continue;
                    }

                    double ratio = a / b;
                    double err = 0;
                    if (a > 0 && b > 0)
                        err = ratio * std::sqrt((ea*ea)/(a*a) + (eb*eb)/(b*b));

                    R->SetBinContent(bx, ratio);
                    R->SetBinError(bx, err);
                }

                R->Write();
                delete A;
                delete B;
            }

            // ------------------------------------------------------------------
            // (2) RAW coarse differences: A(x) - B(x)
            // ------------------------------------------------------------------

            fout->mkdir("XProjectionDiff_Z_InvMassFiltered_CoarseZ_raw");
            fout->cd("XProjectionDiff_Z_InvMassFiltered_CoarseZ_raw");

            for (int i = 0; i < 4; i++) {

                double pos_low  = coarseEdgesZ[4 + i];
                double pos_high = coarseEdgesZ[5 + i];
                double neg_low  = coarseEdgesZ[3 - i];
                double neg_high = coarseEdgesZ[4 - i];

                auto [zpL, zpH] = getZBins(pos_low, pos_high);
                auto [zmL, zmH] = getZBins(neg_low, neg_high);

                TH1D* A = h2_ZPos->ProjectionX(
                    TString::Format("CoarseZ_PosRawDiff_%d", i), zpL, zpH);
                TH1D* B = h2_ZPos->ProjectionX(
                    TString::Format("CoarseZ_NegRawDiff_%d", i), zmL, zmH);

                A->Sumw2();
                B->Sumw2();

                TH1D* D = (TH1D*) A->Clone(TString::Format("CoarseZDiff_%d", i));
                D->SetTitle(
                    TString::Format("Coarse Z difference: Z=[%.1f,%.1f] - Z=[%.1f,%.1f]",
                                    pos_low,pos_high, neg_low,neg_high)
                );
                D->Add(B, -1.0);

                D->GetXaxis()->SetTitle("#Lambda V0 Radius (cm)");
                D->GetYaxis()->SetTitle("A - B");

                D->Write();

                delete A;
                delete B;
            }

            // ------------------------------------------------------------------
            // (3) Full TH2D raw difference: Z>0 minus Z<0, but KEEP ONLY Z>0 range
            // ------------------------------------------------------------------

            fout->mkdir("ZFull2D_Differences_PosOnly");
            fout->cd("ZFull2D_Differences_PosOnly");

            int nX = h2_ZPos->GetNbinsX();
            int nZ = h2_ZPos->GetNbinsY();
            int midZcut = nZ / 2;   // index separating -Z and +Z (Z=0 boundary)

            // Build output TH2D with *only positive Z half*
            int nZpos = nZ - midZcut;

            // Z range edges
            double zPosMin = h2_ZPos->GetYaxis()->GetBinLowEdge(midZcut+1);
            double zPosMax = h2_ZPos->GetYaxis()->GetBinUpEdge(nZ);

            // Create final differential histogram
            TH2D* h2_Z_diff_raw_posOnly = new TH2D(
                "h2_ZFull_Diff_raw_PosMinusNeg_PosOnly",
                "Raw 2D Difference: (Z>0) - (Z<0)  [Only +Z bins]",
                nX,
                h2_ZPos->GetXaxis()->GetXmin(),
                h2_ZPos->GetXaxis()->GetXmax(),
                nZpos,
                zPosMin,
                zPosMax
            );

            h2_Z_diff_raw_posOnly->Sumw2();
            h2_Z_diff_raw_posOnly->GetXaxis()->SetTitle("V0 Radius (cm)");
            h2_Z_diff_raw_posOnly->GetYaxis()->SetTitle("Z vertex (cm)");
            h2_Z_diff_raw_posOnly->SetTitle("Raw Z Difference D(x,z) = N(+z) - N(-z)");

            for (int ix = 1; ix <= nX; ix++) {
                for (int izp = 1; izp <= nZpos; izp++) {

                    int zPlusBin  = midZcut + izp;
                    int zMinusBin = midZcut + 1 - izp;

                    double plus  = h2_ZPos->GetBinContent(ix, zPlusBin);
                    double minus = h2_ZPos->GetBinContent(ix, zMinusBin);

                    double diff = plus - minus;

                    h2_Z_diff_raw_posOnly->SetBinContent(ix, izp, diff);
                    // For simplicity, use sqrt(sum of variances)
                    double eplus  = h2_ZPos->GetBinError(ix, zPlusBin);
                    double eminus = h2_ZPos->GetBinError(ix, zMinusBin);
                    h2_Z_diff_raw_posOnly->SetBinError(
                        ix, izp,
                        std::sqrt(eplus*eplus + eminus*eminus)
                    );
                }
            }

            h2_Z_diff_raw_posOnly->SetOption("COLZ");
            h2_Z_diff_raw_posOnly->Write();


            // ------------------------------------------------------------------
            // (4) Normalized full 2D difference (Z>0-only version)
            // ------------------------------------------------------------------

            fout->mkdir("ZFull2D_Normalized_PosOnly");
            fout->cd("ZFull2D_Normalized_PosOnly");

            // Compute total yields for normalization
            double totPosZ_full = 0.0;
            double totNegZ_full = 0.0;

            for (int ix = 1; ix <= nX; ix++) {
                for (int iz = 1; iz <= nZ; iz++) {

                    double zc = h2_ZPos->GetYaxis()->GetBinCenter(iz);
                    double v  = h2_ZPos->GetBinContent(ix, iz);

                    if (zc > 0) totPosZ_full += v;
                    if (zc < 0) totNegZ_full += v;
                }
            }

            if (totPosZ_full <= 0 || totNegZ_full <= 0) {
                std::cerr << "ERROR: zero positive or negative Z totals! Cannot normalize.\n";
            } else {

                TH2D* h2_Z_diff_norm_posOnly = (TH2D*) h2_Z_diff_raw_posOnly->Clone(
                    "h2_ZFull_Diff_norm_PosMinusNeg_PosOnly"
                );
                h2_Z_diff_norm_posOnly->SetTitle("Normalized 2D Asymmetry: f(+z) - f(-z) [Only +Z]");

                for (int ix = 1; ix <= nX; ix++) {
                    for (int izp = 1; izp <= nZpos; izp++) {

                        int zPlusBin  = midZcut + izp;
                        int zMinusBin = midZcut + 1 - izp;

                        double A = h2_ZPos->GetBinContent(ix, zPlusBin)  / totPosZ_full;
                        double B = h2_ZPos->GetBinContent(ix, zMinusBin) / totNegZ_full;

                        double eA = h2_ZPos->GetBinError(ix, zPlusBin)  / totPosZ_full;
                        double eB = h2_ZPos->GetBinError(ix, zMinusBin) / totNegZ_full;

                        double diff = A - B;
                        double err  = std::sqrt(eA*eA + eB*eB);

                        h2_Z_diff_norm_posOnly->SetBinContent(ix, izp, diff);
                        h2_Z_diff_norm_posOnly->SetBinError(ix, izp, err);
                    }
                }
                h2_Z_diff_norm_posOnly->SetOption("COLZ");
                h2_Z_diff_norm_posOnly->Write();
            }
            std::cout << "Saved positive-Z-only difference histograms.\n";
        }
    }

    std::cout << "=== Finished Part D-extra. ===\n\n";

    // ======================================================================
    // PART E: Normalized coarse Z-ratios
    //         (recompute A(x) and B(x) from the original TH2D)
    // ======================================================================

    std::cout << "\n\n========== PART E: Normalized coarse Z-ratios ==========\n";

    if (!h2_ZPos) {
        std::cerr << "ERROR: h2_ZPos (InvMassFiltered Z TH2D) not found! Skipping Part E.\n";
    } else {

    // Ensure even number of Z bins
    const int nBinsZ_f = h2_ZPos->GetNbinsY();
    if (nBinsZ_f % 2 != 0) {
        std::cerr << "ERROR: Z-axis has odd number of bins (" 
                << nBinsZ_f << ")! Cannot form symmetric ratios.\n";
    }
    const int midZ_f = nBinsZ_f / 2;

    // Create output directory
    fout->mkdir("XProjectionRatios_Z_InvMassFiltered_CoarseZ");
    fout->cd("XProjectionRatios_Z_InvMassFiltered_CoarseZ");

    std::cout << "Building normalized coarse Z-ratios from scratch...\n";

    // --------------------------------------------------------
    // 1. Compute total positive-Z and negative-Z Λ counts
    // --------------------------------------------------------
    double totalPosZ = 0.0;
    double totalNegZ = 0.0;

    for (int i = 1; i <= midZ_f; i++) {

        int zPlus  = midZ_f + i;       // positive side
        int zMinus = midZ_f + 1 - i;   // negative side

        TH1D* projPos = h2_ZPos->ProjectionX(
            TString::Format("ZNorm_PosRaw_tmp_%d", i), zPlus, zPlus);

        TH1D* projNeg = h2_ZPos->ProjectionX(
            TString::Format("ZNorm_NegRaw_tmp_%d", i), zMinus, zMinus);

        totalPosZ += projPos->GetEntries();
        totalNegZ += projNeg->GetEntries();

        delete projPos;
        delete projNeg;
    }

    std::cout << "  Total +Z Λ = " << totalPosZ << "\n";
    std::cout << "  Total -Z Λ = " << totalNegZ << "\n";

    if (totalPosZ <= 0 || totalNegZ <= 0) {
        std::cerr << "ERROR: Zero totalPosZ or totalNegZ — cannot normalize.\n";
    } else {

        // Number of coarse bins we want: Z in [-100,100], 60 bins --> coarse as you choose.
        // Following your Part C choice: 8 coarse bins --> 4 symmetric pairs.
        // Let’s define coarse Z edges explicitly:
        const int nCoarse = 8;
        double coarseEdgesZ[nCoarse + 1];
        double zMin = h2_ZPos->GetYaxis()->GetXmin(); // should be -100
        double zMax = h2_ZPos->GetYaxis()->GetXmax(); // should be +100
        double zStep = (zMax - zMin) / nCoarse;

        for (int i = 0; i <= nCoarse; i++)
            coarseEdgesZ[i] = zMin + i * zStep;

        auto getZFineBins = [&](double z1, double z2) {
            int b1 = h2_ZPos->GetYaxis()->FindBin(z1);
            int b2 = h2_ZPos->GetYaxis()->FindBin(z2);
            if (b1 > b2) std::swap(b1, b2);
            return std::pair<int,int>(b1, b2);
        };

        std::vector<TH1D*> coarseZNormRatios;

        // --------------------------------------------------------
        // 2. Build 4 symmetric coarse Z-ratio histograms
        // --------------------------------------------------------
        for (int i = 0; i < 4; i++) {

            double pos_low  = coarseEdgesZ[4 + i];
            double pos_high = coarseEdgesZ[5 + i];
            double neg_low  = coarseEdgesZ[3 - i];
            double neg_high = coarseEdgesZ[4 - i];

            auto [zpL, zpH] = getZFineBins(pos_low, pos_high);
            auto [zmL, zmH] = getZFineBins(neg_low, neg_high);

            TH1D* Araw = h2_ZPos->ProjectionX(
                TString::Format("CoarseZ_PosRaw_%d", i), zpL, zpH);

            TH1D* Braw = h2_ZPos->ProjectionX(
                TString::Format("CoarseZ_NegRaw_%d", i), zmL, zmH);

            Araw->Sumw2();
            Braw->Sumw2();

            int nXB = Araw->GetNbinsX();

            TH1D* R = new TH1D(
                TString::Format("CoarseZRatio_%d_Norm", i),
                TString::Format("Normalized coarse Z ratio Z=[%.1f,%.1f] / Z=[%.1f,%.1f]",
                                pos_low,pos_high, neg_low,neg_high),
                nXB,
                Araw->GetXaxis()->GetXmin(),
                Araw->GetXaxis()->GetXmax()
            );

            R->Sumw2();
            R->GetXaxis()->SetTitle("#Lambda V0 Radius (cm)");
            R->GetYaxis()->SetTitle("(N_pos/N_totPosZ) / (N_neg/N_totNegZ)");

            // Fill ratio
            for (int bx = 1; bx <= nXB; bx++) {

                double A = Araw->GetBinContent(bx);
                double B = Braw->GetBinContent(bx);
                double eA = Araw->GetBinError(bx);
                double eB = Braw->GetBinError(bx);

                if (A < 0) A = 0;
                if (B < 0) B = 0;

                if (B == 0) {
                    R->SetBinContent(bx, 0);
                    R->SetBinError(bx, 0);
                    continue;
                }

                double A_norm = A / totalPosZ;
                double B_norm = B / totalNegZ;
                double ratio  = A_norm / B_norm;

                double err = 0;
                if (A > 0 && B > 0) {
                    err = ratio * std::sqrt(
                        (eA*eA)/(A*A) + (eB*eB)/(B*B)
                    );
                }

                R->SetBinContent(bx, ratio);
                R->SetBinError(bx, err);
            }

            R->Write();
            coarseZNormRatios.push_back(R);

            delete Araw;
            delete Braw;
        }
        std::cout << "Normalized coarse Z ratios created.\n";
    }
    } // end even number of Z bins check

    std::cout << "=== Finished Part E (coarse Z normalized ratios). ===\n\n";

    // ======================================================================
    // PART F: V0 XY vertex asymmetry (posZ vs negZ)
    // ======================================================================

    std::cout << "\n\n========== PART F: XY vertex asymmetry (posZ vs negZ) ==========\n";

    // ---------------------------------------------------------------
    // Load the TH3Ds for XY vertex vs invariant mass, split by Z sign
    // ---------------------------------------------------------------
    TH3D* h3_XY_posZ = dynamic_cast<TH3D*>( fin->Get("asymmetric_rapidity_test/Lambda/hV0XYvsMass_posZ") );
    TH3D* h3_XY_negZ = dynamic_cast<TH3D*>( fin->Get("asymmetric_rapidity_test/Lambda/hV0XYvsMass_negZ") );

    if (!h3_XY_posZ || !h3_XY_negZ) {
        std::cerr << "ERROR: XY vs Mass histograms not found! Skipping Part F.\n";
    } else {

    // ==================================================================
    // (1) Save original histograms
    // ==================================================================
    fout->mkdir("XY_original");
    fout->cd("XY_original");

    h3_XY_posZ->Write("hV0XYvsMass_posZ_original");
    h3_XY_negZ->Write("hV0XYvsMass_negZ_original");

    // ==================================================================
    // (2) Define folders we will need
    // ==================================================================
    fout->mkdir("XY_MassFiltered");
    fout->mkdir("XY_Differences");
    fout->mkdir("XY_Normalized");
    fout->mkdir("XY_Normalized_Differences");

    // ==================================================================
    // (3) Mass filtering of both XY histograms using mLow, mHigh
    // ==================================================================
    fout->cd("XY_MassFiltered");

    int zbinLow_p  = h3_XY_posZ->GetZaxis()->FindBin(mLow);
    int zbinHigh_p = h3_XY_posZ->GetZaxis()->FindBin(mHigh);
    if (zbinLow_p > zbinHigh_p) std::swap(zbinLow_p, zbinHigh_p);

    int zbinLow_n  = h3_XY_negZ->GetZaxis()->FindBin(mLow);
    int zbinHigh_n = h3_XY_negZ->GetZaxis()->FindBin(mHigh);
    if (zbinLow_n > zbinHigh_n) std::swap(zbinLow_n, zbinHigh_n);

    // Restrict range and project XY plane
    h3_XY_posZ->GetZaxis()->SetRange(zbinLow_p, zbinHigh_p);
    TH2D* h2_XY_posZ = (TH2D*) h3_XY_posZ->Project3D("yx");
    h2_XY_posZ->SetName("hV0XY_posZ_InvMassFiltered");

    h3_XY_negZ->GetZaxis()->SetRange(zbinLow_n, zbinHigh_n);
    TH2D* h2_XY_negZ = (TH2D*) h3_XY_negZ->Project3D("yx");
    h2_XY_negZ->SetName("hV0XY_negZ_InvMassFiltered");

    // Restore full ranges
    h3_XY_posZ->GetZaxis()->SetRange(0, 0);
    h3_XY_negZ->GetZaxis()->SetRange(0, 0);

    // --- VISUALIZATION SETTINGS (Standard plots) ---
    h2_XY_posZ->SetOption("COLZ");
    h2_XY_negZ->SetOption("COLZ");
    // -----------------------------------------------

    // Save filtered TH2Ds
    h2_XY_posZ->Sumw2();
    h2_XY_negZ->Sumw2();
    h2_XY_posZ->Write();
    h2_XY_negZ->Write();

    std::cout << "Mass-filtered XY histograms written.\n";

    // ==================================================================
    // (3.5) Build UNNORMALIZED DIFFERENCE: posZ - negZ
    // ==================================================================
    fout->cd("XY_Differences");

    TH2D* h2_XY_diff_raw =
        (TH2D*) h2_XY_posZ->Clone("hV0XY_Diff_raw_posZ_minus_negZ");

    h2_XY_diff_raw->Add(h2_XY_negZ, -1.0);
    h2_XY_diff_raw->SetTitle("Raw difference: XY(posZ) - XY(negZ)");
    // --- VISUALIZATION SETTINGS ---
    h2_XY_diff_raw->SetOption("COLZ");

    h2_XY_diff_raw->Write();

    std::cout << "Raw XY difference histogram written.\n";

    // ==================================================================
    // (3.6) Build RELATIVE ASYMMETRY: (Pos - Neg) / (Pos + Neg)
    // ==================================================================
    
    // 1. Prepare Numerator: (Pos - Neg)
    // We can just clone the diff we calculated in 3.5, or recalculate to be safe/clean
    TH2D* h2_XY_rel_num = (TH2D*) h2_XY_posZ->Clone("hV0XY_RelDiff_Numerator");
    h2_XY_rel_num->Add(h2_XY_negZ, -1.0);

    // 2. Prepare Denominator: (Pos + Neg) / 2
    TH2D* h2_XY_rel_den = (TH2D*) h2_XY_posZ->Clone("hV0XY_RelDiff_Denominator");
    h2_XY_rel_den->Add(h2_XY_negZ, 1.0); // Add Neg (Sum)
    // h2_XY_rel_den->Scale(0.5);           // Divide by 2 (Average)
    // No need to take the average with that factor of 2!

    // 3. Divide: Num / Den
    // Note: We use the numerator as the base for the final histogram
    TH2D* h2_XY_rel_asym = (TH2D*) h2_XY_rel_num->Clone("hV0XY_Relative_Asymmetry");
    h2_XY_rel_asym->Divide(h2_XY_rel_den);

    // 4. Clean up Labels and Titles
    h2_XY_rel_asym->SetTitle("Relative Asymmetry: (XY(+Z) - XY(-Z))/(XY(+Z) + XY(-Z))");
    // h2_XY_rel_asym->GetZaxis()->SetTitle("(N_{+} - N_{-}) / <N>");

    // 5. Visualization Settings
    h2_XY_rel_asym->SetOption("COLZ");
    h2_XY_rel_asym->Write();
    
    // Clean up temporary histograms to free memory
    delete h2_XY_rel_num;
    delete h2_XY_rel_den;

    std::cout << "Relative Asymmetry histogram written.\n";

    // ==================================================================
    // (4) Create NORMALIZED versions (divide each by total yield)
    // ==================================================================

    fout->cd("XY_Normalized");

    // Compute total positive-Z and total negative-Z counts
    double totalPosXY = h2_XY_posZ->Integral();
    double totalNegXY = h2_XY_negZ->Integral();

    std::cout << "Total +Z mass-filtered XY counts = " << totalPosXY << "\n";
    std::cout << "Total -Z mass-filtered XY counts = " << totalNegXY << "\n";

    if (totalPosXY <= 0 || totalNegXY <= 0) {
        std::cerr << "ERROR: Zero totalPosXY or totalNegXY — cannot normalize Part F.\n";
    } else {

        // Normalize copies
        TH2D* h2_XY_posZ_norm =
            (TH2D*) h2_XY_posZ->Clone("hV0XY_posZ_InvMassFiltered_normalized");
        TH2D* h2_XY_negZ_norm =
            (TH2D*) h2_XY_negZ->Clone("hV0XY_negZ_InvMassFiltered_normalized");

        h2_XY_posZ_norm->Scale(1.0 / totalPosXY);
        h2_XY_negZ_norm->Scale(1.0 / totalNegXY);

        // --- VISUALIZATION SETTINGS (Standard plots) ---
        h2_XY_posZ_norm->SetOption("COLZ");
        h2_XY_negZ_norm->SetOption("COLZ");
        // -----------------------------------------------

        h2_XY_posZ_norm->SetTitle("Normalized XY(+Z): f(x,y | Z>0)");
        h2_XY_negZ_norm->SetTitle("Normalized XY(-Z): f(x,y | Z<0)");

        h2_XY_posZ_norm->Write();
        h2_XY_negZ_norm->Write();

        std::cout << "Normalized XY histograms written.\n";

        // --------------------------------------------------------------
        // (4.5) NORMALIZED DIFFERENCE: posZ_norm - negZ_norm
        // --------------------------------------------------------------

        fout->cd("XY_Normalized_Differences");

        TH2D* h2_XY_diff_norm =
            (TH2D*) h2_XY_posZ_norm->Clone("hV0XY_Diff_norm_posZ_minus_negZ");

        h2_XY_diff_norm->Add(h2_XY_negZ_norm, -1.0);
        h2_XY_diff_norm->SetTitle(
            "Normalized difference: XY(+Z)/N_{+Z} - XY(-Z)/N_{-Z}"
        );
        h2_XY_diff_norm->SetOption("COLZ"); // Show color bar

        h2_XY_diff_norm->Write();

        std::cout << "Normalized XY difference histogram written.\n";

        // --------------------------------------------------------------
        // (4.6) NORMALIZED RELATIVE ASYMMETRY: DiffNorm / SumNorm (SumNorm = AverageNorm*2, if you'd prefer)
        // Formula: (PosNorm - NegNorm) / (PosNorm + NegNorm)
        // --------------------------------------------------------------

        // 1. Numerator: We already have (PosNorm - NegNorm) in h2_XY_diff_norm.
        // We clone it so we don't modify the original difference plot.
        TH2D* h2_XY_rel_norm_num = 
            (TH2D*) h2_XY_diff_norm->Clone("hV0XY_RelDiff_Norm_Num");

        // 2. Denominator: (PosNorm + NegNorm)
        TH2D* h2_XY_rel_norm_den = 
            (TH2D*) h2_XY_posZ_norm->Clone("hV0XY_RelDiff_Norm_Den");
        h2_XY_rel_norm_den->Add(h2_XY_negZ_norm, 1.0); // Add (Sum)
        // h2_XY_rel_norm_den->Scale(0.5);                // Divide by 2 (Average)
        // Actuallu this factor of 2 is not that useful. Just remove it!

        // 3. Divide: Num / Den
        TH2D* h2_XY_rel_norm_asym = 
            (TH2D*) h2_XY_rel_norm_num->Clone("hV0XY_RelNorm_Asymmetry");
        h2_XY_rel_norm_asym->Divide(h2_XY_rel_norm_den);

        // 4. Set Titles and Option
        h2_XY_rel_norm_asym->SetTitle(
            "((XY(+Z)/N_{+Z}) - (XY(-Z)/N_{-Z})) / (XY(+Z)/N_{+Z} + XY(-Z)/N_{-Z})"
        );
        h2_XY_rel_norm_asym->SetOption("COLZ");

        h2_XY_rel_norm_asym->Write();

        // Clean up temporary helper histograms
        delete h2_XY_rel_norm_num;
        delete h2_XY_rel_norm_den;

        std::cout << "Normalized Relative Asymmetry histogram written.\n";
    } // end successful normalization case
    } // end if both XY TH3Ds found

    std::cout << "=== Finished Part F (XY vertex asymmetry) ===\n\n";

    

    // -----------------------------------------------------------------
    // Finalize: write and close files
    // -----------------------------------------------------------------
    fout->Close();
    fin->Close();

    std::cout << "\nAll results written to:\n  " << outFileName << "\n" << std::endl;
}
