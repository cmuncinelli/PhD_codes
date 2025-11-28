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

        // Integrate over many Y bins → produce a TH1D for X
        TH1D* A = h2_massFiltered->ProjectionX(nameA, ypL, ypH);
        TH1D* B = h2_massFiltered->ProjectionX(nameB, ymL, ymH);

        if (A->GetSumw2N() == 0) A->Sumw2();
        if (B->GetSumw2N() == 0) B->Sumw2();

        TString rname = TString::Format("CoarseRatio_%d", i);
        TH1D* R = (TH1D*)A->Clone(rname);
        R->Reset();

        R->SetTitle(
            TString::Format("Coarse ratio y=[%.3f,%.3f] / y=[%.3f,%.3f]",
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

    std::cout << "=== Finished Part C (coarse Y ratios). ===\n\n";
    


    // -----------------------------------------------------------------
    // Finalize: write and close files
    // -----------------------------------------------------------------
    fout->Close();
    fin->Close();

    std::cout << "\nAll results written to:\n  " << outFileName << "\n" << std::endl;
}
