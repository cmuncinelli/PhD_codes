// ==========================================================================
// plotHelicityEfficiency.cxx
// ==========================================================================
//
// PURPOSE
// -------
// Companion plotting macro for helicityEfficiencyToyModel.cxx.
// Reads the output ROOT file produced by that macro and generates a
// comprehensive set of diagnostic figures as PDF files (one PDF per
// figure, saved to a user-specified output directory).
//
// WHAT THIS MACRO SHOWS
// ----------------------
// For each of the four cut scenarios (NoCuts, pTCutOnly, DCACutOnly,
// BothCuts) and each eta half (EtaPos, EtaNeg, All), this macro produces:
//
//   Figure 1: 2D (cos theta*, phi*) maps -- the PRIMARY diagnostic.
//             Four-panel layout: NoCuts / pTCut / DCACut / BothCuts
//             for a chosen eta selection.
//             Expected textures:
//               NoCuts    -> uniform (flat) -- bug check
//               pTCutOnly -> rows in cos(theta*) suppressed near +1
//                            (helicity / forward-backward effect)
//               DCACutOnly-> columns modulated as ~sin(phi*)
//                            (left-right / magnetic-field effect)
//               BothCuts  -> combination of both
//
//   Figure 2: cos(theta*) projections, all four scenarios overlaid,
//             separately for EtaPos, EtaNeg, All.
//             Normalised to area 1 for shape comparison.
//
//   Figure 3: phi* projections, all four scenarios overlaid,
//             separately for EtaPos, EtaNeg, All.
//             A fit to A * (1 + B * sin(phi*)) is drawn on DCACutOnly.
//
//   Figure 4: Ring proxy <R_proxy> as a function of Lambda eta,
//             all four scenarios on one panel.
//             Key feature: sign flip between eta > 0 and eta < 0
//             for DCA-sensitive cuts -- this is the observed antisymmetry.
//
//   Figure 5: Ring proxy <R_proxy> as a function of Lambda pT,
//             all four scenarios on one panel.
//
//   Figure 6: Integrated <R_proxy> bar chart (one value per scenario),
//             separately for EtaPos, EtaNeg, All.
//             Includes a dashed line at zero for reference.
//
//   Figure 7: Daughter pT distributions (proton and pion), pre-cut and
//             post-pT-cut, to visualise which daughter population is
//             being removed by the threshold.
//
//   Figure 8: Daughter DCA_xy distributions (proton and pion), pre-cut
//             and post-DCA-cut, showing the displaced-vertex selection.
//
//   Figure 9: Transverse decay radius distribution (all cuts combined)
//             and Lambda pT / eta spectra from generation.
//
//   Figure 10: eta>0 minus eta<0 difference plots for cos(theta*) and
//              phi*, showing the antisymmetric component explicitly.
//              This directly mimics what is visible in the O-O data.
//
// USAGE
// -----
// From the ROOT command line:
//
//   root -l -b -q 'plotHelicityEfficiency.cxx'
//       (reads helicityEffOutput.root, saves plots to ./plots/)
//
//   root -l -b -q 'plotHelicityEfficiency.cxx("myout.root","myplots")'
//       (custom input file and output directory)
//
// The output directory is created if it does not exist.
//
// DEPENDENCIES
// ------------
// Only standard ROOT headers are needed.  No external libraries required.
// Tested with ROOT 6.28+; should work with ROOT 6.14+.
//
//   Reference analysis code: lambdaJetPolarizationIonsDerived.cxx (ALICE O2)
// ==========================================================================

#include <TFile.h>
#include <TDirectory.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TF1.h>
#include <TStyle.h>
#include <TMath.h>
#include <TSystem.h>
#include <TROOT.h>

#include <cstdio>
#include <cstring>


// ==========================================================================
// STYLE CONSTANTS
// Colors and line widths used consistently across all figures.
// ==========================================================================

// One color per scenario (ROOT color indices)
static const int kColNoCuts  = kBlack;           // no cuts: black
static const int kColPtCut   = kRed + 1;         // pT cut:  red
static const int kColDcaCut  = kBlue + 1;        // DCA cut: blue
static const int kColBoth    = kGreen + 2;       // both:    dark green

static const int kLineWidth  = 2;

// Marker styles and sizes for TProfile overlays
static const int kMarkerNC = 20;  // full circle
static const int kMarkerPT = 21;  // full square
static const int kMarkerDC = 22;  // full triangle up
static const int kMarkerBC = 23;  // full triangle down
static const double kMarkerSz = 0.8;


// ==========================================================================
// HELPER: SafeGet
// --------------------------------------------------------------------------
// Retrieves a named object from a TDirectory with error checking.
// Returns nullptr and prints a warning if the object is not found.
// Cast to the desired type with static_cast<T*>() at the call site.
// ==========================================================================
static TObject* SafeGet(TDirectory* dir, const char* name)
{
    if (!dir) {
        printf("WARNING: SafeGet called with null directory for '%s'\n", name);
        return nullptr;
    }
    TObject* obj = dir->Get(name);
    if (!obj) {
        printf("WARNING: Object '%s' not found in directory '%s'\n",
               name, dir->GetName());
    }
    return obj;
}


// ==========================================================================
// HELPER: GetDir
// --------------------------------------------------------------------------
// Gets a subdirectory from a TDirectory with error checking.
// Prints a warning and returns nullptr if not found.
// ==========================================================================
static TDirectory* GetDir(TDirectory* parent, const char* name)
{
    if (!parent) {
        printf("WARNING: GetDir called with null parent for '%s'\n", name);
        return nullptr;
    }
    TDirectory* dir = (TDirectory*)parent->Get(name);
    if (!dir) {
        printf("WARNING: Directory '%s' not found in '%s'\n",
               name, parent->GetName());
    }
    return dir;
}


// ==========================================================================
// HELPER: SetHistStyle
// --------------------------------------------------------------------------
// Applies a consistent visual style to a TH1 or TProfile.
// ==========================================================================
static void SetHistStyle(TH1* h, int color, int marker, int lineStyle = 1)
{
    if (!h) return;
    h->SetLineColor(color);
    h->SetLineWidth(kLineWidth);
    h->SetLineStyle(lineStyle);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(kMarkerSz);
    h->SetStats(0);  // Turn off the stat box; we add our own annotations
}


// ==========================================================================
// HELPER: SaveCanvas
// --------------------------------------------------------------------------
// Saves a TCanvas as a PDF to the given output directory.
// Creates the directory if it does not exist (using gSystem->Exec).
// ==========================================================================
static void SaveCanvas(TCanvas* c, const char* outDir, const char* filename, TFile *fout)
{
    // Ensure output directory exists
    gSystem->Exec(Form("mkdir -p %s", outDir));

    // Build full path and save
    TString path = Form("%s/%s.pdf", outDir, filename);
    c->SaveAs(path.Data());
    fout->cd();
    c->Write();
    printf("  Saved: %s\n", path.Data());
}


// ==========================================================================
// HELPER: AddLabel
// --------------------------------------------------------------------------
// Draws a TLatex label on the current pad at normalized (x, y) coordinates.
// ==========================================================================
static void AddLabel(double x, double y, const char* text,
                     double size = 0.04, int align = 12)
{
    TLatex* lat = new TLatex(x, y, text);
    lat->SetNDC();
    lat->SetTextSize(size);
    lat->SetTextAlign(align);
    lat->Draw();
}


// ==========================================================================
// HELPER: GetScenarioDir
// --------------------------------------------------------------------------
// Convenience wrapper: returns the directory for a given scenario name and
// eta selection string (one of "EtaPos", "EtaNeg", "All").
// ==========================================================================
static TDirectory* GetScenarioDir(TFile* f,
                                   const char* scenario,
                                   const char* etaSel)
{
    TDirectory* dirScen = GetDir(f, scenario);
    if (!dirScen) return nullptr;
    return GetDir(dirScen, etaSel);
}


// ==========================================================================
// STRUCT: ScenDirs
// --------------------------------------------------------------------------
// Bundles the four scenario directories for a given eta selection so that
// they can be passed together to plotting functions.
// ==========================================================================
struct ScenDirs {
    TDirectory* nc;   // NoCuts
    TDirectory* pt;   // pTCutOnly
    TDirectory* dc;   // DCACutOnly
    TDirectory* bc;   // BothCuts
};

static ScenDirs GetScenDirs(TFile* f, const char* etaSel)
{
    ScenDirs sd;
    sd.nc = GetScenarioDir(f, "NoCuts",     etaSel);
    sd.pt = GetScenarioDir(f, "pTCutOnly",  etaSel);
    sd.dc = GetScenarioDir(f, "DCACutOnly", etaSel);
    sd.bc = GetScenarioDir(f, "BothCuts",   etaSel);
    return sd;
}


// ==========================================================================
// HELPER: MakeLegend
// --------------------------------------------------------------------------
// Creates a 4-entry TLegend for the four scenarios.
// The caller is responsible for placing (x1,y1,x2,y2) appropriately.
// Returns a pointer owned by the current pad (no manual delete needed).
// ==========================================================================
static TLegend* MakeLegend(double x1, double y1, double x2, double y2)
{
    TLegend* leg = new TLegend(x1, y1, x2, y2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    return leg;
}


// ==========================================================================
// FIGURE 1: 2D (cos theta*, phi*) maps
// --------------------------------------------------------------------------
// Four-panel canvas: NoCuts | pTCut | DCACut | BothCuts for one eta half.
// Drawn with COLZ (color palette map).
// ==========================================================================
static void MakeFig1_2DMaps(TFile* f, TFile* fout,
                              const char* etaSel,
                              const char* outDir)
{
    // Names and titles for panel labels
    const char* scenNames[4]  = {"NoCuts", "pTCutOnly", "DCACutOnly", "BothCuts"};
    const char* panelLabels[4]= {"No cuts (bug-check)",
                                  "p_{T} cut only (helicity effect)",
                                  "DCA cut only (left-right effect)",
                                  "Both cuts (combined)"};

    TCanvas* c = new TCanvas(Form("c_2Dmap_%s", etaSel),
                              Form("2D (cos#theta*, #phi*) -- %s", etaSel),
                              1400, 350);
    c->Divide(4, 1, 0.002, 0.002);

    // Global style: use palette "kRainBow" (ROOT 6) for readability
    gStyle->SetPalette(kRainBow);
    gStyle->SetNumberContours(64);

    for (int i = 0; i < 4; ++i) {
        c->cd(i + 1);
        gPad->SetRightMargin(0.15);  // Room for the Z-axis (COLZ)
        gPad->SetLeftMargin(0.12);

        TDirectory* dir = GetScenarioDir(f, scenNames[i], etaSel);
        if (!dir) continue;

        TH2D* h = static_cast<TH2D*>(SafeGet(dir, "h2d_cosTheta_phi"));
        if (!h) continue;

        // Draw with color palette
        h->SetStats(0);
        h->GetXaxis()->SetTitleSize(0.055);
        h->GetYaxis()->SetTitleSize(0.055);
        h->Draw("COLZ");

        // Panel label at top
        AddLabel(0.5, 0.94, panelLabels[i], 0.042, 22);  // centered, top
    }

    c->cd(0);  // Return to the canvas (not a pad) for the super-title
    AddLabel(0.5, 0.995,
             Form("Proton emission angles in #Lambda rest frame -- %s", etaSel),
             0.038, 22);

    SaveCanvas(c, outDir, Form("fig1_2Dmap_%s", etaSel), fout);
    delete c;
}


// ==========================================================================
// FIGURE 2: cos(theta*) projections, all scenarios overlaid
// --------------------------------------------------------------------------
// Normalised to area 1 so shape differences are visible.
// One canvas per eta selection (EtaPos, EtaNeg, All).
// ==========================================================================
static void MakeFig2_CosTheta(TFile* f, TFile* fout,
                                const char* etaSel,
                                const char* outDir)
{
    ScenDirs sd = GetScenDirs(f, etaSel);

    // Retrieve histograms; clone so we can normalise without modifying originals
    TH1D* hNC = sd.nc ? static_cast<TH1D*>(
        SafeGet(sd.nc, "h1d_cosTheta")) : nullptr;
    TH1D* hPT = sd.pt ? static_cast<TH1D*>(
        SafeGet(sd.pt, "h1d_cosTheta")) : nullptr;
    TH1D* hDC = sd.dc ? static_cast<TH1D*>(
        SafeGet(sd.dc, "h1d_cosTheta")) : nullptr;
    TH1D* hBC = sd.bc ? static_cast<TH1D*>(
        SafeGet(sd.bc, "h1d_cosTheta")) : nullptr;

    if (!hNC) return;  // Need at least the no-cuts histogram to set axis range

    // Clone and normalise each to unit area for shape comparison
    auto NormClone = [](TH1D* h, const char* suffix) -> TH1D* {
        if (!h) return nullptr;
        TH1D* hc = static_cast<TH1D*>(h->Clone(Form("%s_%s", h->GetName(), suffix)));
        hc->SetDirectory(nullptr);  // Detach from any TDirectory
        double integral = hc->Integral();
        if (integral > 0.) hc->Scale(1. / integral);
        return hc;
    };

    TH1D* hNCn = NormClone(hNC, "normNC");
    TH1D* hPTn = NormClone(hPT, "normPT");
    TH1D* hDCn = NormClone(hDC, "normDC");
    TH1D* hBCn = NormClone(hBC, "normBC");

    // Apply styles
    SetHistStyle(hNCn, kColNoCuts, kMarkerNC);
    SetHistStyle(hPTn, kColPtCut,  kMarkerPT);
    SetHistStyle(hDCn, kColDcaCut, kMarkerDC);
    SetHistStyle(hBCn, kColBoth,   kMarkerBC);

    TCanvas* c = new TCanvas(Form("c_cosTheta_%s", etaSel),
                              Form("cos(theta*) -- %s", etaSel), 700, 600);
    c->SetLeftMargin(0.14);
    c->SetBottomMargin(0.13);

    // Determine y-axis range (10% above maximum across all histograms)
    double ymax = 0.;
    for (TH1D* h : {hNCn, hPTn, hDCn, hBCn}) {
        if (h && h->GetMaximum() > ymax) ymax = h->GetMaximum();
    }
    ymax *= 1.25;

    hNCn->SetTitle(Form("cos(#theta*) -- %s;cos(#theta*);Normalised counts", etaSel));
    hNCn->GetYaxis()->SetRangeUser(0., ymax);
    hNCn->Draw("HIST");
    if (hPTn) hPTn->Draw("HIST SAME");
    if (hDCn) hDCn->Draw("HIST SAME");
    if (hBCn) hBCn->Draw("HIST SAME");

    // Reference line at the flat (no-cut) value = 1 / N_bins
    // (visual guide for where suppression appears)
    double flatVal = (hNCn->GetNbinsX() > 0) ? 1. / hNCn->GetNbinsX() : 0.;
    TLine* flatLine = new TLine(-1., flatVal, 1., flatVal);
    flatLine->SetLineColor(kGray + 2);
    flatLine->SetLineStyle(2);
    flatLine->SetLineWidth(1);
    flatLine->Draw();

    // Legend
    TLegend* leg = MakeLegend(0.60, 0.68, 0.88, 0.88);
    leg->AddEntry(hNCn, "No cuts (bug-check)", "l");
    if (hPTn) leg->AddEntry(hPTn, "p_{T} cut only", "l");
    if (hDCn) leg->AddEntry(hDCn, "DCA cut only",   "l");
    if (hBCn) leg->AddEntry(hBCn, "Both cuts",      "l");
    leg->Draw();

    // Annotation about the expected asymmetry direction
    AddLabel(0.16, 0.87,
             "Expect: suppression near cos#theta* ~ +1 after p_{T} cut",
             0.033, 12);
    AddLabel(0.16, 0.82,
             "(proton emitted backward in #Lambda frame is easier to reconstruct)",
             0.028, 12);

    SaveCanvas(c, outDir, Form("fig2_cosTheta_%s", etaSel), fout);
    // Clean up clones
    delete hNCn; delete hPTn; delete hDCn; delete hBCn;
    delete c;
}


// ==========================================================================
// FIGURE 3: phi* projections, all scenarios overlaid
// --------------------------------------------------------------------------
// A fit to A*(1 + B*sin(phi*)) is drawn on the DCACutOnly histogram
// to measure the left-right modulation amplitude.
// ==========================================================================
static void MakeFig3_PhiStar(TFile* f, TFile* fout,
                               const char* etaSel,
                               const char* outDir)
{
    ScenDirs sd = GetScenDirs(f, etaSel);

    TH1D* hNC = sd.nc ? static_cast<TH1D*>(SafeGet(sd.nc, "h1d_phi")) : nullptr;
    TH1D* hPT = sd.pt ? static_cast<TH1D*>(SafeGet(sd.pt, "h1d_phi")) : nullptr;
    TH1D* hDC = sd.dc ? static_cast<TH1D*>(SafeGet(sd.dc, "h1d_phi")) : nullptr;
    TH1D* hBC = sd.bc ? static_cast<TH1D*>(SafeGet(sd.bc, "h1d_phi")) : nullptr;

    if (!hNC) return;

    // Clone and normalise
    auto NormClone = [](TH1D* h, const char* suf) -> TH1D* {
        if (!h) return nullptr;
        TH1D* hc = static_cast<TH1D*>(h->Clone(Form("%s_%s", h->GetName(), suf)));
        hc->SetDirectory(nullptr);
        double ig = hc->Integral();
        if (ig > 0.) hc->Scale(1. / ig);
        return hc;
    };

    TH1D* hNCn = NormClone(hNC, "normNC");
    TH1D* hPTn = NormClone(hPT, "normPT");
    TH1D* hDCn = NormClone(hDC, "normDC");
    TH1D* hBCn = NormClone(hBC, "normBC");

    SetHistStyle(hNCn, kColNoCuts, kMarkerNC);
    SetHistStyle(hPTn, kColPtCut,  kMarkerPT);
    SetHistStyle(hDCn, kColDcaCut, kMarkerDC);
    SetHistStyle(hBCn, kColBoth,   kMarkerBC);

    TCanvas* c = new TCanvas(Form("c_phiStar_%s", etaSel),
                              Form("phi* -- %s", etaSel), 700, 600);
    c->SetLeftMargin(0.14);
    c->SetBottomMargin(0.13);

    double ymax = 0.;
    for (TH1D* h : {hNCn, hPTn, hDCn, hBCn}) {
        if (h && h->GetMaximum() > ymax) ymax = h->GetMaximum();
    }
    ymax *= 1.35;

    hNCn->SetTitle(Form("#phi* -- %s;#phi* [rad];Normalised counts", etaSel));
    hNCn->GetYaxis()->SetRangeUser(0., ymax);
    hNCn->Draw("HIST");
    if (hPTn) hPTn->Draw("HIST SAME");
    if (hDCn) hDCn->Draw("HIST SAME");
    if (hBCn) hBCn->Draw("HIST SAME");

    // --- Fit sin modulation to the DCA-cut histogram ---
    // Model: f(phi) = A * (1 + B * sin(phi))
    // where A = 1/N_bins (flat component) and B is the asymmetry amplitude.
    // This is the left-right asymmetry expected from the DCA-magnetic-field effect.
    if (hDCn && hDCn->Integral() > 0.) {
        // The normalised histogram has bin width = 2pi/N_bins, so the
        // flat value per bin is 1/N_bins.  We fit amplitude * (1 + B*sin(phi)).
        double nbins  = hDCn->GetNbinsX();
        double flatV  = 1. / nbins;

        // Fit function: p0 * (1 + p1 * sin(x + p2))
        // p0 ~ flatV, p1 ~ asymmetry magnitude, p2 ~ phase offset
        TF1* fSin = new TF1("fSinDCA",
                             "[0] * (1. + [1]*TMath::Sin(x + [2]))",
                             -TMath::Pi(), TMath::Pi());
        fSin->SetParameters(flatV, 0.1, 0.);    // Initial guesses
        fSin->SetParLimits(0, 0.5*flatV, 2.*flatV);   // A must be positive
        fSin->SetParLimits(1, -1., 1.);          // |B| <= 1
        fSin->SetLineColor(kBlue - 3);
        fSin->SetLineWidth(2);
        fSin->SetLineStyle(2);
        hDCn->Fit(fSin, "RQ");  // "Q" = quiet (no printout per fit step)
        fSin->Draw("SAME");

        // Report fit result
        printf("  Fig 3 [%s] DCA-cut sin-fit: A=%.4f, B=%.4f +/- %.4f, phase=%.4f\n",
               etaSel,
               fSin->GetParameter(0),
               fSin->GetParameter(1), fSin->GetParError(1),
               fSin->GetParameter(2));
    }

    // Legend
    TLegend* leg = MakeLegend(0.13, 0.73, 0.55, 0.90);
    leg->AddEntry(hNCn, "No cuts (expected: flat)",    "l");
    if (hPTn) leg->AddEntry(hPTn, "p_{T} cut (expected: flat)",   "l");
    if (hDCn) leg->AddEntry(hDCn, "DCA cut (expected: sin modulation)", "l");
    if (hBCn) leg->AddEntry(hBCn, "Both cuts",                          "l");
    leg->Draw();

    AddLabel(0.16, 0.67,
             "DCA cut -> sin(#phi*) modulation = left-right asymmetry",
             0.033, 12);
    AddLabel(0.16, 0.62,
             "Dashed blue line: A(1 + B sin(#phi* + #delta)) fit to DCA-cut",
             0.028, 12);

    SaveCanvas(c, outDir, Form("fig3_phiStar_%s", etaSel), fout);
    delete hNCn; delete hPTn; delete hDCn; delete hBCn;
    delete c;
}


// ==========================================================================
// FIGURE 4: <R_proxy> vs Lambda eta -- the eta antisymmetry plot
// --------------------------------------------------------------------------
// This is the single most important plot for understanding the O-O data.
// The DCA-cut curve should be approximately antisymmetric: positive at
// eta > 0, negative at eta < 0 (or vice versa, depending on B-field sign).
// The pT-cut curve should be approximately zero everywhere.
// ==========================================================================
static void MakeFig4_RingVsEta(TFile* f, TFile* fout, const char* outDir)
{
    // For this figure we use the "All" (combined eta) directory, since
    // pRingProxyVsEta spans both halves on one profile.
    ScenDirs sd = GetScenDirs(f, "All");

    TProfile* pNC = sd.nc ? static_cast<TProfile*>(
        SafeGet(sd.nc, "pRingProxyVsEta")) : nullptr;
    TProfile* pPT = sd.pt ? static_cast<TProfile*>(
        SafeGet(sd.pt, "pRingProxyVsEta")) : nullptr;
    TProfile* pDC = sd.dc ? static_cast<TProfile*>(
        SafeGet(sd.dc, "pRingProxyVsEta")) : nullptr;
    TProfile* pBC = sd.bc ? static_cast<TProfile*>(
        SafeGet(sd.bc, "pRingProxyVsEta")) : nullptr;

    if (!pNC) return;

    SetHistStyle(pNC, kColNoCuts, kMarkerNC);
    SetHistStyle(pPT, kColPtCut,  kMarkerPT);
    SetHistStyle(pDC, kColDcaCut, kMarkerDC);
    SetHistStyle(pBC, kColBoth,   kMarkerBC);

    TCanvas* c = new TCanvas("c_RingVsEta",
                              "<R_{proxy}> vs #Lambda #eta", 750, 600);
    c->SetLeftMargin(0.14);
    c->SetBottomMargin(0.13);

    // Find symmetric y-axis range
    double ymax = 0.;
    for (TProfile* p : {pNC, pPT, pDC, pBC}) {
        if (!p) continue;
        for (int ib = 1; ib <= p->GetNbinsX(); ++ib) {
            double val = std::fabs(p->GetBinContent(ib) +
                                    p->GetBinError(ib));
            if (val > ymax) ymax = val;
        }
    }
    ymax = (ymax < 1.e-6) ? 0.05 : ymax * 1.35;

    pNC->SetTitle("<R_{proxy}> vs #Lambda pseudorapidity;"
                  "#eta_{#Lambda};<R_{proxy}>");
    pNC->GetYaxis()->SetRangeUser(-ymax, ymax);
    pNC->Draw("EP");
    if (pPT) pPT->Draw("EP SAME");
    if (pDC) pDC->Draw("EP SAME");
    if (pBC) pBC->Draw("EP SAME");

    // Zero reference line
    TLine* zeroline = new TLine(-0.9, 0., 0.9, 0.);
    zeroline->SetLineColor(kGray + 2);
    zeroline->SetLineStyle(2);
    zeroline->SetLineWidth(1);
    zeroline->Draw();

    // Vertical line at eta = 0 to highlight the antisymmetry axis
    TLine* etazeroline = new TLine(0., -ymax, 0., ymax);
    etazeroline->SetLineColor(kGray + 1);
    etazeroline->SetLineStyle(3);
    etazeroline->SetLineWidth(1);
    etazeroline->Draw();

    TLegend* leg = MakeLegend(0.15, 0.72, 0.50, 0.89);
    leg->AddEntry(pNC, "No cuts",      "ep");
    if (pPT) leg->AddEntry(pPT, "p_{T} cut only", "ep");
    if (pDC) leg->AddEntry(pDC, "DCA cut only",   "ep");
    if (pBC) leg->AddEntry(pBC, "Both cuts",       "ep");
    leg->Draw();

    AddLabel(0.16, 0.67,
             "Key: DCA-cut curve should be ANTISYMMETRIC in #eta",
             0.033, 12);
    AddLabel(0.16, 0.62,
             "(sign flip at #eta = 0 mirrors the fake signal in O-O data)",
             0.028, 12);

    SaveCanvas(c, outDir, "fig4_ringProxyVsEta", fout);
    delete c;
}


// ==========================================================================
// FIGURE 5: <R_proxy> vs Lambda pT
// ==========================================================================
static void MakeFig5_RingVsPt(TFile* f, TFile* fout, const char* outDir)
{
    ScenDirs sd = GetScenDirs(f, "All");

    TProfile* pNC = sd.nc ? static_cast<TProfile*>(
        SafeGet(sd.nc, "pRingProxyVsPt")) : nullptr;
    TProfile* pPT = sd.pt ? static_cast<TProfile*>(
        SafeGet(sd.pt, "pRingProxyVsPt")) : nullptr;
    TProfile* pDC = sd.dc ? static_cast<TProfile*>(
        SafeGet(sd.dc, "pRingProxyVsPt")) : nullptr;
    TProfile* pBC = sd.bc ? static_cast<TProfile*>(
        SafeGet(sd.bc, "pRingProxyVsPt")) : nullptr;

    if (!pNC) return;

    SetHistStyle(pNC, kColNoCuts, kMarkerNC);
    SetHistStyle(pPT, kColPtCut,  kMarkerPT);
    SetHistStyle(pDC, kColDcaCut, kMarkerDC);
    SetHistStyle(pBC, kColBoth,   kMarkerBC);

    double ymax = 0.;
    for (TProfile* p : {pNC, pPT, pDC, pBC}) {
        if (!p) continue;
        for (int ib = 1; ib <= p->GetNbinsX(); ++ib) {
            double v = std::fabs(p->GetBinContent(ib) + p->GetBinError(ib));
            if (v > ymax) ymax = v;
        }
    }
    ymax = (ymax < 1.e-6) ? 0.05 : ymax * 1.35;

    TCanvas* c = new TCanvas("c_RingVsPt",
                              "<R_{proxy}> vs #Lambda p_{T}", 750, 600);
    c->SetLeftMargin(0.14);
    c->SetBottomMargin(0.13);

    pNC->SetTitle("<R_{proxy}> vs #Lambda p_{T};"
                  "p_{T}^{#Lambda} [GeV/c];<R_{proxy}>");
    pNC->GetYaxis()->SetRangeUser(-ymax, ymax);
    pNC->Draw("EP");
    if (pPT) pPT->Draw("EP SAME");
    if (pDC) pDC->Draw("EP SAME");
    if (pBC) pBC->Draw("EP SAME");

    TLine* zeroline = new TLine(0., 0., 5., 0.);
    zeroline->SetLineColor(kGray + 2);
    zeroline->SetLineStyle(2);
    zeroline->Draw();

    TLegend* leg = MakeLegend(0.55, 0.72, 0.88, 0.89);
    leg->AddEntry(pNC, "No cuts",      "ep");
    if (pPT) leg->AddEntry(pPT, "p_{T} cut only", "ep");
    if (pDC) leg->AddEntry(pDC, "DCA cut only",   "ep");
    if (pBC) leg->AddEntry(pBC, "Both cuts",       "ep");
    leg->Draw();

    AddLabel(0.16, 0.22,
             "Expect: DCA-cut signal may have p_{T} dependence from helix geometry",
             0.030, 12);

    SaveCanvas(c, outDir, "fig5_ringProxyVsPt", fout);
    delete c;
}


// ==========================================================================
// FIGURE 6: Integrated <R_proxy> bar-chart per scenario
// --------------------------------------------------------------------------
// Uses the single-bin TProfile pRingProxy from each scenario+eta directory.
// Three panels side by side: EtaPos | EtaNeg | All.
// ==========================================================================
static void MakeFig6_IntegratedRing(TFile* f, TFile* fout, const char* outDir)
{
    // We will manually build a TH1D with one bin per scenario
    // to represent the bar chart, since ROOT's bar chart requires TH1.
    // Scenarios: NoCuts, pTCutOnly, DCACutOnly, BothCuts (4 entries)
    // Three versions: EtaPos, EtaNeg, All

    const char* scenNames[4]  = {"NoCuts", "pTCutOnly", "DCACutOnly", "BothCuts"};
    const char* etaSels[3]    = {"EtaPos", "EtaNeg", "All"};
    const char* etaLabels[3]  = {"#eta > 0", "#eta < 0", "All #eta"};
    const int   scenColors[4] = {kColNoCuts, kColPtCut, kColDcaCut, kColBoth};
    const char* scenLabels[4] = {"No cuts", "p_{T} cut", "DCA cut", "Both"};

    TCanvas* c = new TCanvas("c_intRing",
                              "Integrated <R_{proxy}>", 1050, 500);
    c->Divide(3, 1, 0.004, 0.002);

    // Find global y-range across all panels for consistent axes
    double globalMax = 0.;
    for (int ie = 0; ie < 3; ++ie) {
        for (int is = 0; is < 4; ++is) {
            TDirectory* dir = GetScenarioDir(f, scenNames[is], etaSels[ie]);
            if (!dir) continue;
            TProfile* p = static_cast<TProfile*>(SafeGet(dir, "pRingProxy"));
            if (!p) continue;
            double v = std::fabs(p->GetBinContent(1)) + p->GetBinError(1);
            if (v > globalMax) globalMax = v;
        }
    }
    globalMax = (globalMax < 1.e-6) ? 0.05 : globalMax * 1.5;

    for (int ie = 0; ie < 3; ++ie) {
        c->cd(ie + 1);
        gPad->SetLeftMargin(0.18);
        gPad->SetBottomMargin(0.16);

        // Build a TH1D of width 1 per scenario (bin label on x-axis)
        TH1D* hBar = new TH1D(Form("hBar_%s", etaSels[ie]),
                               Form("<R_{proxy}> -- %s;"
                                    "Scenario;<R_{proxy}>", etaLabels[ie]),
                               4, 0., 4.);
        hBar->GetYaxis()->SetRangeUser(-globalMax, globalMax);
        hBar->SetStats(0);

        // Draw individual points with error bars using a temporary holder
        // per scenario, overlaid on the same pad
        hBar->Draw();  // Draws axes; bars added below via individual draws

        // Zero reference line
        TLine* zl = new TLine(0., 0., 4., 0.);
        zl->SetLineColor(kGray + 2);
        zl->SetLineStyle(2);
        zl->SetLineWidth(2);
        zl->Draw();

        // For each scenario, place a marker at the bin center
        for (int is = 0; is < 4; ++is) {
            TDirectory* dir = GetScenarioDir(f, scenNames[is], etaSels[ie]);
            if (!dir) continue;
            TProfile* p = static_cast<TProfile*>(SafeGet(dir, "pRingProxy"));
            if (!p) continue;

            double val = p->GetBinContent(1);
            double err = p->GetBinError(1);

            // Draw as a vertical error bar + horizontal tick at the center
            // x position = is + 0.5 (center of bin is+1 in the hBar axis)
            double xc = is + 0.5;
            TLine* barV = new TLine(xc, val - err, xc, val + err);
            barV->SetLineColor(scenColors[is]);
            barV->SetLineWidth(3);
            barV->Draw();

            // Horizontal caps on the error bar
            double capW = 0.12;
            TLine* capU = new TLine(xc - capW, val + err,
                                     xc + capW, val + err);
            TLine* capL = new TLine(xc - capW, val - err,
                                     xc + capW, val - err);
            capU->SetLineColor(scenColors[is]);
            capU->SetLineWidth(2);
            capL->SetLineColor(scenColors[is]);
            capL->SetLineWidth(2);
            capU->Draw(); capL->Draw();

            // Central value marker (filled circle)
            TLine* cent = new TLine(xc - capW, val, xc + capW, val);
            cent->SetLineColor(scenColors[is]);
            cent->SetLineWidth(4);
            cent->Draw();

            // x-axis label below the tick
            hBar->GetXaxis()->SetBinLabel(is + 1, scenLabels[is]);
        }
        hBar->GetXaxis()->SetLabelSize(0.06);

        AddLabel(0.5, 0.96, etaLabels[ie], 0.045, 22);
    }

    c->cd(0);
    AddLabel(0.5, 0.995,
             "Integrated <R_{proxy}> per cut scenario and #eta selection",
             0.038, 22);

    SaveCanvas(c, outDir, "fig6_integratedRingProxy", fout);
    delete c;
}


// ==========================================================================
// FIGURE 7: Daughter pT distributions
// --------------------------------------------------------------------------
// Shows proton and pion pT before any cut and after the pT cut.
// Two-panel canvas: proton (left), pion (right).
// The threshold line shows where the minimum pT cut falls.
// ==========================================================================
static void MakeFig7_DaughterPt(TFile* f, TFile* fout, const char* outDir)
{
    // Pre-cut histograms are in the Kinematics directory
    TDirectory* dirKin = GetDir(f, "Kinematics");
    if (!dirKin) return;

    // Post-pT-cut histograms from the NoCuts/All and pTCutOnly/All dirs
    // For DCA-only-cut scenario, pT is NOT cut, so we read from pTCutOnly
    // which already has pT cut applied.
    TDirectory* dirPT_All = GetScenarioDir(f, "pTCutOnly", "All");

    TH1D* hKin_p  = static_cast<TH1D*>(SafeGet(dirKin, "hKin_pT_proton"));
    TH1D* hKin_pi = static_cast<TH1D*>(SafeGet(dirKin, "hKin_pT_pion"));
    TH1D* hCut_p  = dirPT_All ?
        static_cast<TH1D*>(SafeGet(dirPT_All, "h1d_pT_proton")) : nullptr;
    TH1D* hCut_pi = dirPT_All ?
        static_cast<TH1D*>(SafeGet(dirPT_All, "h1d_pT_pion")) : nullptr;

    if (!hKin_p || !hKin_pi) return;

    TCanvas* c = new TCanvas("c_daughterPt",
                              "Daughter p_{T} distributions", 1000, 500);
    c->Divide(2, 1, 0.004, 0.002);

    for (int id = 0; id < 2; ++id) {
        c->cd(id + 1);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.14);
        gPad->SetBottomMargin(0.13);

        TH1D* hPre = (id == 0) ? hKin_p  : hKin_pi;
        TH1D* hPost = (id == 0) ? hCut_p : hCut_pi;

        // Clone, normalise pre-cut to 1, scale post-cut by same factor
        TH1D* hPreN  = static_cast<TH1D*>(hPre->Clone(Form("hPreN_%d",  id)));
        TH1D* hPostN = hPost ?
            static_cast<TH1D*>(hPost->Clone(Form("hPostN_%d", id))) : nullptr;
        hPreN->SetDirectory(nullptr);
        if (hPostN) hPostN->SetDirectory(nullptr);

        double norm = hPreN->Integral();
        if (norm > 0.) {
            hPreN->Scale(1. / norm);
            if (hPostN) hPostN->Scale(1. / norm);
        }

        SetHistStyle(hPreN, kBlack, kMarkerNC);
        if (hPostN) SetHistStyle(hPostN, kRed + 1, kMarkerPT);

        const char* pname = (id == 0) ? "Proton" : "Pion";
        hPreN->SetTitle(Form("%s p_{T};"
                             "p_{T}^{%s} [GeV/c];Normalised counts",
                             pname, pname));
        hPreN->GetYaxis()->SetRangeUser(1.e-5, 1.0);
        hPreN->Draw("HIST");
        if (hPostN) hPostN->Draw("HIST SAME");

        // Draw vertical line at the cut threshold
        // We cannot read it from the ROOT file, so we label it generically
        // (the user knows their cut value from the run parameters printout).
        AddLabel(0.50, 0.88, Form("%s p_{T} dist.", pname), 0.04, 22);
        AddLabel(0.50, 0.82, "Black: pre-cut   Red: after p_{T} cut", 0.032, 22);
        AddLabel(0.50, 0.76,
                 "(log scale -- look for low-pT depletion after cut)",
                 0.028, 22);

        delete hPreN;
        if (hPostN) delete hPostN;
    }

    SaveCanvas(c, outDir, "fig7_daughterPt", fout);
    delete c;
}


// ==========================================================================
// FIGURE 8: Daughter DCA_xy distributions
// ==========================================================================
static void MakeFig8_DaughterDCA(TFile* f, TFile* fout, const char* outDir)
{
    TDirectory* dirKin = GetDir(f, "Kinematics");
    if (!dirKin) return;
    TDirectory* dirDCA_All = GetScenarioDir(f, "DCACutOnly", "All");

    TH1D* hKin_p  = static_cast<TH1D*>(SafeGet(dirKin, "hKin_DCA_proton"));
    TH1D* hKin_pi = static_cast<TH1D*>(SafeGet(dirKin, "hKin_DCA_pion"));
    TH1D* hCut_p  = dirDCA_All ?
        static_cast<TH1D*>(SafeGet(dirDCA_All, "h1d_DCA_proton")) : nullptr;
    TH1D* hCut_pi = dirDCA_All ?
        static_cast<TH1D*>(SafeGet(dirDCA_All, "h1d_DCA_pion")) : nullptr;

    if (!hKin_p || !hKin_pi) return;

    TCanvas* c = new TCanvas("c_daughterDCA",
                              "Daughter DCA_{xy} distributions", 1000, 500);
    c->Divide(2, 1, 0.004, 0.002);

    for (int id = 0; id < 2; ++id) {
        c->cd(id + 1);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.14);
        gPad->SetBottomMargin(0.13);

        TH1D* hPre  = (id == 0) ? hKin_p  : hKin_pi;
        TH1D* hPost = (id == 0) ? hCut_p  : hCut_pi;

        TH1D* hPreN  = static_cast<TH1D*>(
            hPre->Clone(Form("hDCApreN_%d", id)));
        TH1D* hPostN = hPost ?
            static_cast<TH1D*>(hPost->Clone(Form("hDCApostN_%d", id))) : nullptr;
        hPreN->SetDirectory(nullptr);
        if (hPostN) hPostN->SetDirectory(nullptr);

        double norm = hPreN->Integral();
        if (norm > 0.) {
            hPreN->Scale(1. / norm);
            if (hPostN) hPostN->Scale(1. / norm);
        }

        SetHistStyle(hPreN, kBlack, kMarkerNC);
        if (hPostN) SetHistStyle(hPostN, kBlue + 1, kMarkerDC);

        const char* pname = (id == 0) ? "Proton" : "Pion";
        hPreN->SetTitle(Form("%s DCA_{xy};"
                             "DCA_{xy}^{%s} [cm];Normalised counts",
                             pname, pname));
        hPreN->Draw("HIST");
        if (hPostN) hPostN->Draw("HIST SAME");

        AddLabel(0.50, 0.88, Form("%s DCA_{xy}", pname), 0.04, 22);
        AddLabel(0.50, 0.82, "Black: pre-cut   Blue: after DCA cut", 0.032, 22);
        AddLabel(0.50, 0.76,
                 "(note depletion at small DCA after cut)",
                 0.028, 22);

        delete hPreN;
        if (hPostN) delete hPostN;
    }

    SaveCanvas(c, outDir, "fig8_daughterDCA", fout);
    delete c;
}


// ==========================================================================
// FIGURE 9: Lambda kinematic spectra and decay radius
// --------------------------------------------------------------------------
// Two-row, three-column layout:
//   Row 1: Lambda pT | Lambda eta | Lambda phi
//   Row 2: Decay radius (all) | Proton pT pre-cut | Pion pT pre-cut
// ==========================================================================
static void MakeFig9_LambdaKin(TFile* f, TFile* fout, const char* outDir)
{
    TDirectory* dirKin = GetDir(f, "Kinematics");
    if (!dirKin) return;

    TH1D* hpT  = static_cast<TH1D*>(SafeGet(dirKin, "hKin_pT_lambda"));
    TH1D* heta = static_cast<TH1D*>(SafeGet(dirKin, "hKin_eta_lambda"));
    TH1D* hphi = static_cast<TH1D*>(SafeGet(dirKin, "hKin_phi_lambda"));
    TH1D* hR   = static_cast<TH1D*>(SafeGet(dirKin, "hKin_decayR"));
    TH1D* hpP  = static_cast<TH1D*>(SafeGet(dirKin, "hKin_pT_proton"));
    TH1D* hpPi = static_cast<TH1D*>(SafeGet(dirKin, "hKin_pT_pion"));

    TCanvas* c = new TCanvas("c_lamKin",
                              "#Lambda kinematics", 1200, 700);
    c->Divide(3, 2, 0.003, 0.003);

    struct PadInfo { TH1D* h; bool logy; const char* note; };
    PadInfo pads[6] = {
        {hpT,  true,
         "Input #Lambda p_{T} spectrum (Boltzmann mT)"},
        {heta, false,
         "Generated #Lambda #eta (flat in rapidity -> slight dip at |#eta|~0.9)"},
        {hphi, false,
         "Generated #Lambda #phi (should be uniform)"},
        {hR,   true,
         "Transverse decay radius (exponential in proper length)"},
        {hpP,  true,
         "Proton p_{T} before any cut"},
        {hpPi, true,
         "Pion p_{T} before any cut"}
    };

    for (int i = 0; i < 6; ++i) {
        c->cd(i + 1);
        if (pads[i].logy) gPad->SetLogy();
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        if (!pads[i].h) continue;
        pads[i].h->SetStats(0);
        pads[i].h->SetLineColor(kBlack);
        pads[i].h->SetLineWidth(2);
        pads[i].h->Draw("HIST");
        AddLabel(0.5, 0.90, pads[i].note, 0.032, 22);
    }

    SaveCanvas(c, outDir, "fig9_lambdaKinematics", fout);
    delete c;
}


// ==========================================================================
// FIGURE 10: eta>0 minus eta<0 difference plots
// --------------------------------------------------------------------------
// This directly shows the antisymmetric component of each distribution,
// which is what appears in the data as a non-zero ring observable when
// averaging over positive and negative eta separately.
//
// For the phi* histogram:  difference = EtaPos - EtaNeg (normalised)
// For the cos(theta*) histogram: same.
// For <R_proxy> vs eta: the profile already shows this when plotted as Fig 4.
//
// Here we produce the histogram-level antisymmetric component.
// ==========================================================================
static void MakeFig10_EtaDifference(TFile* f, TFile* fout, const char* outDir)
{
    const char* scenNames[4]  = {"NoCuts", "pTCutOnly", "DCACutOnly", "BothCuts"};
    const int   scenColors[4] = {kColNoCuts, kColPtCut, kColDcaCut, kColBoth};
    const char* scenLabels[4] = {"No cuts", "p_{T} cut", "DCA cut", "Both"};

    TCanvas* c = new TCanvas("c_etaDiff",
                              "#eta>0 minus #eta<0 difference", 1200, 500);
    c->Divide(3, 1, 0.003, 0.002);  // phi*, cos(theta*), ring proxy vs eta diff

    // ---- Panel 1: phi* difference ----
    c->cd(1);
    gPad->SetLeftMargin(0.14);
    bool firstDrawn1 = false;
    double ymax1 = 0.;
    TLegend* leg1 = MakeLegend(0.15, 0.70, 0.60, 0.90);

    for (int is = 0; is < 4; ++is) {
        TDirectory* dirP = GetScenarioDir(f, scenNames[is], "EtaPos");
        TDirectory* dirN = GetScenarioDir(f, scenNames[is], "EtaNeg");
        if (!dirP || !dirN) continue;

        TH1D* hP = static_cast<TH1D*>(SafeGet(dirP, "h1d_phi"));
        TH1D* hN = static_cast<TH1D*>(SafeGet(dirN, "h1d_phi"));
        if (!hP || !hN) continue;

        // Normalise each to its integral so counts differences don't dominate
        TH1D* hPn = static_cast<TH1D*>(hP->Clone(Form("hPhi_Pn_%d", is)));
        TH1D* hNn = static_cast<TH1D*>(hN->Clone(Form("hPhi_Nn_%d", is)));
        hPn->SetDirectory(nullptr);
        hNn->SetDirectory(nullptr);
        double iP = hPn->Integral(); if (iP > 0.) hPn->Scale(1./iP);
        double iN = hNn->Integral(); if (iN > 0.) hNn->Scale(1./iN);

        // Difference: positive eta minus negative eta
        hPn->Add(hNn, -1.);

        SetHistStyle(hPn, scenColors[is], kMarkerNC + is);

        double locmax = std::max(std::fabs(hPn->GetMinimum()),
                                  hPn->GetMaximum());
        if (locmax > ymax1) ymax1 = locmax;

        if (!firstDrawn1) {
            hPn->SetTitle("#phi* (#eta>0) - (#eta<0) difference;"
                           "#phi* [rad];#Delta(normalised counts)");
            hPn->GetYaxis()->SetRangeUser(-0.001, 0.001); // placeholder
            hPn->Draw("HIST");
            firstDrawn1 = true;
        } else {
            hPn->Draw("HIST SAME");
        }
        leg1->AddEntry(hPn, scenLabels[is], "l");
    }
    if (firstDrawn1) {
        // Update axis range using actual data
        gPad->GetListOfPrimitives()->At(0);  // The first histogram
        ymax1 = (ymax1 < 1.e-8) ? 0.002 : ymax1 * 1.35;
        // Note: changing the range after Draw requires GetHistogram()
        // We add a zero reference line instead
        TLine* zl1 = new TLine(-TMath::Pi(), 0., TMath::Pi(), 0.);
        zl1->SetLineColor(kGray + 2); zl1->SetLineStyle(2); zl1->Draw();
    }
    leg1->Draw();
    AddLabel(0.5, 0.94,
             "Expect: sin-like modulation for DCA cut",
             0.030, 22);

    // ---- Panel 2: cos(theta*) difference ----
    c->cd(2);
    gPad->SetLeftMargin(0.14);
    bool firstDrawn2 = false;
    TLegend* leg2 = MakeLegend(0.15, 0.70, 0.60, 0.90);

    for (int is = 0; is < 4; ++is) {
        TDirectory* dirP = GetScenarioDir(f, scenNames[is], "EtaPos");
        TDirectory* dirN = GetScenarioDir(f, scenNames[is], "EtaNeg");
        if (!dirP || !dirN) continue;

        TH1D* hP = static_cast<TH1D*>(SafeGet(dirP, "h1d_cosTheta"));
        TH1D* hN = static_cast<TH1D*>(SafeGet(dirN, "h1d_cosTheta"));
        if (!hP || !hN) continue;

        TH1D* hPn = static_cast<TH1D*>(hP->Clone(Form("hCos_Pn_%d", is)));
        TH1D* hNn = static_cast<TH1D*>(hN->Clone(Form("hCos_Nn_%d", is)));
        hPn->SetDirectory(nullptr);
        hNn->SetDirectory(nullptr);
        double iP = hPn->Integral(); if (iP > 0.) hPn->Scale(1./iP);
        double iN = hNn->Integral(); if (iN > 0.) hNn->Scale(1./iN);
        hPn->Add(hNn, -1.);

        SetHistStyle(hPn, scenColors[is], kMarkerNC + is);

        if (!firstDrawn2) {
            hPn->SetTitle("cos(#theta*) (#eta>0) - (#eta<0) difference;"
                           "cos(#theta*);#Delta(normalised counts)");
            hPn->Draw("HIST");
            firstDrawn2 = true;
        } else {
            hPn->Draw("HIST SAME");
        }
        leg2->AddEntry(hPn, scenLabels[is], "l");
    }
    if (firstDrawn2) {
        TLine* zl2 = new TLine(-1., 0., 1., 0.);
        zl2->SetLineColor(kGray + 2); zl2->SetLineStyle(2); zl2->Draw();
    }
    leg2->Draw();
    AddLabel(0.5, 0.94,
             "Expect: near-flat difference (helicity effect is eta-symmetric)",
             0.030, 22);

    // ---- Panel 3: Integrated ring proxy EtaPos minus EtaNeg bar chart ----
    c->cd(3);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.16);

    // Build values manually
    double vals[4] = {0.}, errs[4] = {0.};
    for (int is = 0; is < 4; ++is) {
        TDirectory* dirP = GetScenarioDir(f, scenNames[is], "EtaPos");
        TDirectory* dirN = GetScenarioDir(f, scenNames[is], "EtaNeg");
        if (!dirP || !dirN) continue;
        TProfile* pP = static_cast<TProfile*>(SafeGet(dirP, "pRingProxy"));
        TProfile* pN = static_cast<TProfile*>(SafeGet(dirN, "pRingProxy"));
        if (!pP || !pN) continue;
        vals[is] = pP->GetBinContent(1) - pN->GetBinContent(1);
        errs[is] = std::sqrt(pP->GetBinError(1) * pP->GetBinError(1) +
                              pN->GetBinError(1) * pN->GetBinError(1));
    }

    double ymax3 = 0.;
    for (int is = 0; is < 4; ++is)
        if (std::fabs(vals[is]) + errs[is] > ymax3)
            ymax3 = std::fabs(vals[is]) + errs[is];
    ymax3 = (ymax3 < 1.e-8) ? 0.01 : ymax3 * 1.5;

    TH1D* hDummy = new TH1D("hDummy3",
                              "<R_{proxy}>_{#eta>0} - <R_{proxy}>_{#eta<0};"
                              "Scenario;#Delta<R_{proxy}>",
                              4, 0., 4.);
    hDummy->SetStats(0);
    hDummy->GetYaxis()->SetRangeUser(-ymax3, ymax3);
    for (int is = 0; is < 4; ++is) hDummy->GetXaxis()->SetBinLabel(is+1, scenLabels[is]);
    hDummy->GetXaxis()->SetLabelSize(0.06);
    hDummy->Draw();

    TLine* zl3 = new TLine(0., 0., 4., 0.);
    zl3->SetLineColor(kGray + 2); zl3->SetLineStyle(2); zl3->SetLineWidth(2);
    zl3->Draw();

    for (int is = 0; is < 4; ++is) {
        double xc = is + 0.5;
        double capW = 0.12;
        // Central value bar
        TLine* barC = new TLine(xc - capW, vals[is], xc + capW, vals[is]);
        barC->SetLineColor(scenColors[is]); barC->SetLineWidth(4); barC->Draw();
        // Error bar
        TLine* barE = new TLine(xc, vals[is]-errs[is], xc, vals[is]+errs[is]);
        barE->SetLineColor(scenColors[is]); barE->SetLineWidth(2); barE->Draw();
    }

    AddLabel(0.5, 0.94,
             "Antisymmetric component: non-zero = fake ring signal",
             0.030, 22);
    AddLabel(0.5, 0.88,
             "Expect large value for DCA cut; ~0 for pT cut",
             0.026, 22);

    c->cd(0);
    AddLabel(0.5, 0.997,
             "#eta>0 minus #eta<0 difference plots: isolating the antisymmetric fake signal",
             0.036, 22);

    SaveCanvas(c, outDir, "fig10_etaDifference", fout);
    delete c;
}


// ==========================================================================
// MAIN FUNCTION: plotHelicityEfficiency
// --------------------------------------------------------------------------
// Reads the output ROOT file from helicityEfficiencyToyModel.cxx and
// generates all 10 diagnostic figure PDFs.
//
// Parameters:
//   inputFile  -- path to the ROOT file produced by the toy model
//                 (default: "helicityEffOutput.root")
//   outDir     -- directory where PDF figures are saved
//                 (default: "plots")
// ==========================================================================
void plotHelicityEfficiency(
    const char* inputFile = "helicityEffOutput.root",
    const char* outDir    = "plots")
{
    // Suppress ROOT info messages to keep terminal output clean
    gROOT->ProcessLine("gErrorIgnoreLevel = kWarning;");

    printf("\n");
    printf("========================================================\n");
    printf("  plotHelicityEfficiency\n");
    printf("  Input : %s\n", inputFile);
    printf("  Output: %s/\n", outDir);
    printf("========================================================\n\n");

    TString rootOutName = (TString) outDir + "/../plots.root"; // To see it all in a single .root file. Saves on the same level as the helicity_baseline.root file
    TFile *fout = new TFile(rootOutName, "RECREATE");

    // Open the input file in READ mode
    TFile* f = TFile::Open(inputFile, "READ");
    if (!f || f->IsZombie()) {
        printf("ERROR: Cannot open file '%s'. Aborting.\n", inputFile);
        return;
    }

    // ---- Global style settings ----
    gStyle->SetOptStat(0);          // No stat boxes
    gStyle->SetOptTitle(0);         // No automatic title boxes (we draw manually)
    gStyle->SetPadTickX(1);         // Ticks on both sides of X axis
    gStyle->SetPadTickY(1);         // Ticks on both sides of Y axis
    gStyle->SetFrameLineWidth(1);
    gStyle->SetTitleSize(0.05, "XYZ");
    gStyle->SetLabelSize(0.045, "XYZ");
    gStyle->SetTitleOffset(1.1, "Y");


    // ---- Figure 1: 2D maps for each eta selection ----
    printf("Producing Figure 1: 2D (cos theta*, phi*) maps...\n");
    for (const char* eta : {"EtaPos", "EtaNeg", "All"})
        MakeFig1_2DMaps(f, fout, eta, outDir);

    // ---- Figure 2: cos(theta*) projections ----
    printf("Producing Figure 2: cos(theta*) projections...\n");
    for (const char* eta : {"EtaPos", "EtaNeg", "All"})
        MakeFig2_CosTheta(f, fout, eta, outDir);

    // ---- Figure 3: phi* projections with sin fit ----
    printf("Producing Figure 3: phi* projections...\n");
    for (const char* eta : {"EtaPos", "EtaNeg", "All"})
        MakeFig3_PhiStar(f, fout, eta, outDir);

    // ---- Figure 4: <R_proxy> vs Lambda eta ----
    printf("Producing Figure 4: <R_proxy> vs eta...\n");
    MakeFig4_RingVsEta(f, fout, outDir);

    // ---- Figure 5: <R_proxy> vs Lambda pT ----
    printf("Producing Figure 5: <R_proxy> vs pT...\n");
    MakeFig5_RingVsPt(f, fout, outDir);

    // ---- Figure 6: Integrated ring proxy bar chart ----
    printf("Producing Figure 6: Integrated <R_proxy>...\n");
    MakeFig6_IntegratedRing(f, fout, outDir);

    // ---- Figure 7: Daughter pT distributions ----
    printf("Producing Figure 7: Daughter pT distributions...\n");
    MakeFig7_DaughterPt(f, fout, outDir);

    // ---- Figure 8: Daughter DCA distributions ----
    printf("Producing Figure 8: Daughter DCA distributions...\n");
    MakeFig8_DaughterDCA(f, fout, outDir);

    // ---- Figure 9: Lambda kinematics ----
    printf("Producing Figure 9: Lambda kinematics...\n");
    MakeFig9_LambdaKin(f, fout, outDir);

    // ---- Figure 10: eta>0 minus eta<0 antisymmetry ----
    printf("Producing Figure 10: eta antisymmetry difference plots...\n");
    MakeFig10_EtaDifference(f, fout, outDir);

    f->Close();
    fout->Close();

    printf("\n========================================================\n");
    printf("  All figures saved to: %s/\n", outDir);
    printf("  Total: 10 multi-panel figures (13 PDF files).\n");
    printf("========================================================\n\n");
}
// end of plotHelicityEfficiency.cxx
