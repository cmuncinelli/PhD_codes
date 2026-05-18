// ==========================================================================
// plotHelicityEfficiency.cxx
// ==========================================================================
//
// PURPOSE
// -------
// Companion plotting macro for helicityEfficiencyToyModel.cxx.
// Reads the ROOT file produced by that macro and writes a comprehensive set
// of diagnostic TCanvas objects to a single plots.root file placed next to
// the input file.  No PDFs are produced.
//
// The macro processes both histogram families produced by the generator:
//   WithEtaGate     -- daughter |eta| < etaMaxDetector required (consistent set)
//   WithoutEtaGate  -- no daughter eta requirement (kept for comparison)
//
// For each family the following figures are produced:
//
//   fig1_2Dmap_<eta>     -- 4-panel COLZ: cos(theta*) vs phi* for each scenario
//   fig2_cosTheta_<eta>  -- cos(theta*) projections, all 4 scenarios overlaid
//   fig3_phiStar_<eta>   -- phi* projections, all 4 scenarios overlaid
//   fig4_ringProxy       -- 2-panel: R_proxy distributions + <R_proxy> vs Lambda eta
//   fig5_ringVsPt        -- <R_proxy> vs Lambda pT
//   fig6_intRing         -- integrated <R_proxy> bar chart (EtaPos | EtaNeg | All)
//   fig7_daughterPt      -- proton/pion pT before and after the pT cut
//   fig8_daughterDCA     -- proton/pion DCA_xy before and after the DCA cut
//   fig9_etaDiff         -- eta>0 minus eta<0 antisymmetric component
//
// Additionally, one shared figure is produced at the top level:
//   figKin               -- Lambda kinematics and decay radius (generation-level)
//
// OUTPUT
// ------
// The output file is plots.root, placed in the same directory as inputFile.
// Canvas objects are written under subdirectories WithEtaGate/ and
// WithoutEtaGate/ inside plots.root.  The kinematics canvas is at the top level.
//
// USAGE
// -----
//   root -l -b -q 'plotHelicityEfficiency.cxx'
//       (reads helicityEffOutput.root; writes plots.root next to it)
//
//   root -l -b -q 'plotHelicityEfficiency.cxx("myout.root","unused_dir")'
//       (second argument kept for run-script compatibility; only input path matters)
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
#include <TStyle.h>
#include <TMath.h>
#include <TSystem.h>
#include <TROOT.h>

#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>


// ==========================================================================
// STYLE CONSTANTS
// ==========================================================================

static const int    kColNoCuts  = kBlack;
static const int    kColPtCut   = kRed + 1;
static const int    kColDcaCut  = kBlue + 1;
static const int    kColBoth    = kGreen + 2;
static const int    kLineWidth  = 2;
static const int    kMarkerNC   = 20;   // full circle
static const int    kMarkerPT   = 21;   // full square
static const int    kMarkerDC   = 22;   // full triangle up
static const int    kMarkerBC   = 23;   // full triangle down
static const double kMarkerSz   = 0.8;

// Scenario names, labels, colors and markers in a fixed order
// (NoCuts, pTCutOnly, DCACutOnly, BothCuts)
static const char* kScenNames[4]  = {"NoCuts", "pTCutOnly", "DCACutOnly", "BothCuts"};
static const char* kScenLabels[4] = {"No cuts", "p_{T} cut only", "DCA cut only", "Both cuts"};
static const int   kScenColors[4] = {kColNoCuts, kColPtCut, kColDcaCut, kColBoth};
static const int   kScenMarkers[4]= {kMarkerNC, kMarkerPT, kMarkerDC, kMarkerBC};


// ==========================================================================
/**
 * @brief Retrieves a named object from a TDirectory with a null-safety check.
 * @param dir   Source directory; a null pointer is silently tolerated.
 * @param name  Name of the object to retrieve.
 * @return      Pointer to the object, or nullptr if not found.
 */
// ==========================================================================
static TObject* SafeGet(TDirectory* dir, const char* name)
{
    if (!dir) { printf("WARNING: SafeGet called with null directory for '%s'\n", name); return nullptr; }
    TObject* obj = dir->Get(name);
    if (!obj) printf("WARNING: '%s' not found in directory '%s'\n", name, dir->GetName());
    return obj;
}


// ==========================================================================
/**
 * @brief Gets a named sub-directory from a TDirectory with a null-safety check.
 * @param parent  Parent directory; null is tolerated.
 * @param name    Sub-directory name.
 * @return        Pointer to the sub-directory, or nullptr if not found.
 */
// ==========================================================================
static TDirectory* GetDir(TDirectory* parent, const char* name)
{
    if (!parent) { printf("WARNING: GetDir called with null parent for '%s'\n", name); return nullptr; }
    TDirectory* dir = static_cast<TDirectory*>(parent->Get(name));
    if (!dir) printf("WARNING: directory '%s' not found in '%s'\n", name, parent->GetName());
    return dir;
}


// ==========================================================================
/**
 * @brief Applies a consistent line/marker/color style to a TH1 or TProfile.
 * @param h          Histogram to style; a null pointer is silently tolerated.
 * @param color      ROOT color index.
 * @param marker     ROOT marker style.
 * @param lineStyle  ROOT line style (default: solid).
 */
// ==========================================================================
static void SetHistStyle(TH1* h, int color, int marker, int lineStyle = 1)
{
    if (!h) return;
    h->SetLineColor(color);  h->SetLineWidth(kLineWidth);  h->SetLineStyle(lineStyle);
    h->SetMarkerColor(color); h->SetMarkerStyle(marker);   h->SetMarkerSize(kMarkerSz);
    h->SetStats(0);
}


// ==========================================================================
/**
 * @brief Draws a TLatex string on the current pad at normalized coordinates.
 * @param x      NDC x position.
 * @param y      NDC y position.
 * @param text   LaTeX string.
 * @param size   Text size (default 0.04).
 * @param align  Text alignment code (default 12 = left-middle).
 */
// ==========================================================================
static void AddLabel(double x, double y, const char* text, double size = 0.04, int align = 12)
{
    TLatex* lat = new TLatex(x, y, text);
    lat->SetNDC();
    lat->SetTextSize(size);
    lat->SetTextAlign(align);
    lat->Draw();
}


// ==========================================================================
/**
 * @brief Writes a TCanvas to a TDirectory and returns to the file root.
 *
 * Replaces the old SaveCanvas/PDF approach.  The canvas is written directly
 * to @p outDir inside the output ROOT file; no disk directory or PDF is
 * created.
 *
 * @param c       Canvas to write.
 * @param outDir  Target TDirectory in the output ROOT file.
 */
// ==========================================================================
static void WriteCanvas(TCanvas* c, TDirectory* outDir)
{
    if (!c || !outDir) return;
    outDir->cd();
    c->Write();
}


// ==========================================================================
/**
 * @brief Creates a TLegend with the house style (no border, transparent fill).
 * @param x1,y1,x2,y2  NDC corner coordinates.
 * @return              Pointer to the new TLegend (owned by the current pad).
 */
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
/**
 * @brief Returns the leaf TDirectory for a given scenario and eta selection
 *        inside a family directory.
 *
 * Path: familyDir / scenarioName / etaSel
 *
 * @param familyDir    Top-level family directory (WithEtaGate or WithoutEtaGate).
 * @param scenarioName One of "NoCuts", "pTCutOnly", "DCACutOnly", "BothCuts".
 * @param etaSel       One of "EtaPos", "EtaNeg", "All".
 * @return             Leaf TDirectory, or nullptr if any step fails.
 */
// ==========================================================================
static TDirectory* GetScenarioDir(TDirectory* familyDir, const char* scenarioName, const char* etaSel)
{
    TDirectory* dirScen = GetDir(familyDir, scenarioName);
    if (!dirScen) return nullptr;
    return GetDir(dirScen, etaSel);
}


// ==========================================================================
/**
 * @brief Bundles the four scenario leaf directories for one (family, etaSel)
 *        combination into a single struct for convenience.
 */
// ==========================================================================
struct ScenDirs {
    TDirectory* nc = nullptr;   ///< NoCuts
    TDirectory* pt = nullptr;   ///< pTCutOnly
    TDirectory* dc = nullptr;   ///< DCACutOnly
    TDirectory* bc = nullptr;   ///< BothCuts
};

static ScenDirs GetScenDirs(TDirectory* famDir, const char* etaSel)
{
    ScenDirs sd;
    sd.nc = GetScenarioDir(famDir, "NoCuts",     etaSel);
    sd.pt = GetScenarioDir(famDir, "pTCutOnly",  etaSel);
    sd.dc = GetScenarioDir(famDir, "DCACutOnly", etaSel);
    sd.bc = GetScenarioDir(famDir, "BothCuts",   etaSel);
    return sd;
}

static int gCloneIdx = 0; // Global counter to ensure unique clone names
static TH1D* SafeClone(TH1D* h) {
    if (!h) return nullptr;
    TH1D* hc = static_cast<TH1D*>(h->Clone(Form("%s_clone_%d", h->GetName(), gCloneIdx++)));
    
    // CRITICAL FIX: Detaches the clone from the output ROOT file, preventing 
    // histogram pollution and ensuring the TCanvas renders the overlays correctly!
    hc->SetDirectory(nullptr); 
    return hc;
}

// ==========================================================================
/**
 * @brief Fig 1 -- 4-panel COLZ map of cos(theta*) vs phi* for each scenario.
 *
 * One canvas per (family, etaSel) combination.  The four panels show NoCuts,
 * pTCutOnly, DCACutOnly, and BothCuts side by side.
 *
 * @param famDir    Family directory in the input file.
 * @param famOut    Corresponding sub-directory in the output ROOT file.
 * @param famLabel  Short label used in canvas and object names.
 * @param etaSel    One of "EtaPos", "EtaNeg", "All".
 */
// ==========================================================================
static void MakeFig1_2DMaps(TDirectory* famDir, TDirectory* famOut,
                              const char* famLabel, const char* etaSel)
{
    TCanvas* c = new TCanvas(Form("c_%s_2Dmap_%s", famLabel, etaSel), Form("cos#theta* vs #phi* -- %s -- %s", famLabel, etaSel), 1400, 350);
    c->Divide(4, 1, 0.002, 0.002);

    for (int i = 0; i < 4; ++i) {
        c->cd(i + 1);
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.12);
        TDirectory* dir = GetScenarioDir(famDir, kScenNames[i], etaSel);
        if (!dir) continue;
        TH2D* h = static_cast<TH2D*>(SafeGet(dir, "h2d_cosTheta_phi"));
        if (!h) continue;
        h->SetStats(0);
        h->GetXaxis()->SetTitleSize(0.055);
        h->GetYaxis()->SetTitleSize(0.055);
        // Removed "SAME" so the COLZ palette renders properly for all frames
        h->Draw("COLZ");
        AddLabel(0.5, 0.94, kScenLabels[i], 0.042, 22);
    }

    c->cd(0);
    AddLabel(0.5, 0.995, Form("Proton emission angles in #Lambda rest frame -- %s -- %s", famLabel, etaSel), 0.036, 22);
    WriteCanvas(c, famOut);
    delete c;
}


// ==========================================================================
/**
 * @brief Fig 2 -- cos(theta*) projections, all four scenarios normalised to
 *        unit area and overlaid on one canvas.
 *
 * A dashed reference line marks the flat (1/N_bins) level.
 *
 * @param famDir    Family directory in the input file.
 * @param famOut    Output sub-directory.
 * @param famLabel  Short label used in object names.
 * @param etaSel    One of "EtaPos", "EtaNeg", "All".
 */
// ==========================================================================
static void MakeFig2_CosTheta(TDirectory* famDir, TDirectory* famOut,
                                const char* famLabel, const char* etaSel)
{
    ScenDirs sd = GetScenDirs(famDir, etaSel);

    TH1D* hNC = sd.nc ? static_cast<TH1D*>(SafeGet(sd.nc, "h1d_cosTheta")) : nullptr;
    TH1D* hPT = sd.pt ? static_cast<TH1D*>(SafeGet(sd.pt, "h1d_cosTheta")) : nullptr;
    TH1D* hDC = sd.dc ? static_cast<TH1D*>(SafeGet(sd.dc, "h1d_cosTheta")) : nullptr;
    TH1D* hBC = sd.bc ? static_cast<TH1D*>(SafeGet(sd.bc, "h1d_cosTheta")) : nullptr;
    if (!hNC) return;

    TH1D* hNCn = SafeClone(hNC);
    TH1D* hPTn = SafeClone(hPT);
    TH1D* hDCn = SafeClone(hDC);
    TH1D* hBCn = SafeClone(hBC);

    SetHistStyle(hNCn, kColNoCuts, kMarkerNC);  SetHistStyle(hPTn, kColPtCut,  kMarkerPT);
    SetHistStyle(hDCn, kColDcaCut, kMarkerDC);  SetHistStyle(hBCn, kColBoth,   kMarkerBC);

    double ymax = 0.;
    for (TH1D* h : {hNCn, hPTn, hDCn, hBCn}) { if (h && h->GetMaximum() > ymax) ymax = h->GetMaximum(); }
    ymax *= 1.25;

    TCanvas* c = new TCanvas(Form("c_%s_cosTheta_%s", famLabel, etaSel), "", 700, 600);
    c->SetLeftMargin(0.14);  c->SetBottomMargin(0.13);

    hNCn->SetTitle(Form("cos(#theta*) -- %s -- %s;cos(#theta*);Counts", famLabel, etaSel));
    hNCn->GetYaxis()->SetRangeUser(0., ymax);
    hNCn->Draw("HIST");
    if (hPTn) hPTn->Draw("HIST SAME");
    if (hDCn) hDCn->Draw("HIST SAME");
    if (hBCn) hBCn->Draw("HIST SAME");

    TLegend* leg = MakeLegend(0.60, 0.68, 0.88, 0.88);
    leg->AddEntry(hNCn, kScenLabels[0], "l");
    if (hPTn) leg->AddEntry(hPTn, kScenLabels[1], "l");
    if (hDCn) leg->AddEntry(hDCn, kScenLabels[2], "l");
    if (hBCn) leg->AddEntry(hBCn, kScenLabels[3], "l");
    leg->Draw("SAME");

    WriteCanvas(c, famOut);
    // Removed deletes so the canvas keeps the objects!
    delete c;
}


// ==========================================================================
/**
 * @brief Fig 3 -- phi* projections, all four scenarios normalised to unit
 *        area and overlaid on one canvas.
 *
 * @param famDir    Family directory in the input file.
 * @param famOut    Output sub-directory.
 * @param famLabel  Short label used in object names.
 * @param etaSel    One of "EtaPos", "EtaNeg", "All".
 */
// ==========================================================================
static void MakeFig3_PhiStar(TDirectory* famDir, TDirectory* famOut,
                               const char* famLabel, const char* etaSel)
{
    ScenDirs sd = GetScenDirs(famDir, etaSel);

    TH1D* hNC = sd.nc ? static_cast<TH1D*>(SafeGet(sd.nc, "h1d_phi")) : nullptr;
    TH1D* hPT = sd.pt ? static_cast<TH1D*>(SafeGet(sd.pt, "h1d_phi")) : nullptr;
    TH1D* hDC = sd.dc ? static_cast<TH1D*>(SafeGet(sd.dc, "h1d_phi")) : nullptr;
    TH1D* hBC = sd.bc ? static_cast<TH1D*>(SafeGet(sd.bc, "h1d_phi")) : nullptr;
    if (!hNC) return;

    TH1D* hNCn = SafeClone(hNC);
    TH1D* hPTn = SafeClone(hPT);
    TH1D* hDCn = SafeClone(hDC);
    TH1D* hBCn = SafeClone(hBC);

    SetHistStyle(hNCn, kColNoCuts, kMarkerNC);  SetHistStyle(hPTn, kColPtCut,  kMarkerPT);
    SetHistStyle(hDCn, kColDcaCut, kMarkerDC);  SetHistStyle(hBCn, kColBoth,   kMarkerBC);

    double ymax = 0.;
    for (TH1D* h : {hNCn, hPTn, hDCn, hBCn}) { if (h && h->GetMaximum() > ymax) ymax = h->GetMaximum(); }
    ymax *= 1.35;

    TCanvas* c = new TCanvas(Form("c_%s_phiStar_%s", famLabel, etaSel), "", 700, 600);
    c->SetLeftMargin(0.14);  c->SetBottomMargin(0.13);

    hNCn->SetTitle(Form("#phi* -- %s -- %s;#phi* [rad];Counts", famLabel, etaSel));
    hNCn->GetYaxis()->SetRangeUser(0., ymax);
    hNCn->Draw("HIST");
    if (hPTn) hPTn->Draw("HIST SAME");
    if (hDCn) hDCn->Draw("HIST SAME");
    if (hBCn) hBCn->Draw("HIST SAME");

    TLegend* leg = MakeLegend(0.13, 0.73, 0.55, 0.90);
    leg->AddEntry(hNCn, kScenLabels[0], "l");
    if (hPTn) leg->AddEntry(hPTn, kScenLabels[1], "l");
    if (hDCn) leg->AddEntry(hDCn, kScenLabels[2], "l");
    if (hBCn) leg->AddEntry(hBCn, kScenLabels[3], "l");
    leg->Draw("SAME");

    WriteCanvas(c, famOut);
    delete c;
}


// ==========================================================================
/**
 * @brief Fig 4 -- 2-panel ring observable proxy figure.
 *
 * Left panel: distributions of R_proxy for all four scenarios ("All" eta),
 * normalised to unit area.  Right panel: profile <R_proxy> vs Lambda eta for
 * all four scenarios.
 *
 * @param famDir    Family directory in the input file.
 * @param famOut    Output sub-directory.
 * @param famLabel  Short label used in object names.
 */
// ==========================================================================
static void MakeFig4_RingProxy(TDirectory* famDir, TDirectory* famOut, const char* famLabel)
{
    // --- Retrieve distributions (h1d_ringProxy) from "All" directories ---
    TH1D* hDistNC = nullptr;  TH1D* hDistPT = nullptr;
    TH1D* hDistDC = nullptr;  TH1D* hDistBC = nullptr;
    {
        TDirectory* dNC = GetScenarioDir(famDir, "NoCuts",     "All");
        TDirectory* dPT = GetScenarioDir(famDir, "pTCutOnly",  "All");
        TDirectory* dDC = GetScenarioDir(famDir, "DCACutOnly", "All");
        TDirectory* dBC = GetScenarioDir(famDir, "BothCuts",   "All");
        if (dNC) hDistNC = static_cast<TH1D*>(SafeGet(dNC, "h1d_ringProxy"));
        if (dPT) hDistPT = static_cast<TH1D*>(SafeGet(dPT, "h1d_ringProxy"));
        if (dDC) hDistDC = static_cast<TH1D*>(SafeGet(dDC, "h1d_ringProxy"));
        if (dBC) hDistBC = static_cast<TH1D*>(SafeGet(dBC, "h1d_ringProxy"));
    }

    // --- Retrieve <R_proxy> vs eta profiles from "All" ---
    ScenDirs sd = GetScenDirs(famDir, "All");
    TProfile* pNC = sd.nc ? static_cast<TProfile*>(SafeGet(sd.nc, "pRingProxyVsEta")) : nullptr;
    TProfile* pPT = sd.pt ? static_cast<TProfile*>(SafeGet(sd.pt, "pRingProxyVsEta")) : nullptr;
    TProfile* pDC = sd.dc ? static_cast<TProfile*>(SafeGet(sd.dc, "pRingProxyVsEta")) : nullptr;
    TProfile* pBC = sd.bc ? static_cast<TProfile*>(SafeGet(sd.bc, "pRingProxyVsEta")) : nullptr;
    if (!pNC && !hDistNC) return;

    TH1D* hDnc = SafeClone(hDistNC);
    TH1D* hDpt = SafeClone(hDistPT);
    TH1D* hDdc = SafeClone(hDistDC);
    TH1D* hDbc = SafeClone(hDistBC);
    
    SetHistStyle(hDnc, kColNoCuts, kMarkerNC);  SetHistStyle(hDpt, kColPtCut,  kMarkerPT);
    SetHistStyle(hDdc, kColDcaCut, kMarkerDC);  SetHistStyle(hDbc, kColBoth,   kMarkerBC);

    SetHistStyle(pNC, kColNoCuts, kMarkerNC);   SetHistStyle(pPT, kColPtCut,  kMarkerPT);
    SetHistStyle(pDC, kColDcaCut, kMarkerDC);   SetHistStyle(pBC, kColBoth,   kMarkerBC);

    TCanvas* c = new TCanvas(Form("c_%s_ringProxy", famLabel), "", 1200, 550);
    c->Divide(2, 1, 0.004, 0.002);

    // ---- Left panel: R_proxy distribution ----
    c->cd(1);
    gPad->SetLeftMargin(0.14);  gPad->SetBottomMargin(0.13);
    {
        double ymax = 0.;
        for (TH1D* h : {hDnc, hDpt, hDdc, hDbc}) { if (h && h->GetMaximum() > ymax) ymax = h->GetMaximum(); }
        ymax *= 1.3;
        if (hDnc) {
            hDnc->SetTitle(Form("R_{proxy} distribution -- %s;R_{proxy};Counts", famLabel));
            hDnc->GetYaxis()->SetRangeUser(0., ymax);
            hDnc->Draw("HIST");
        }
        if (hDpt) hDpt->Draw("HIST SAME");
        if (hDdc) hDdc->Draw("HIST SAME");
        if (hDbc) hDbc->Draw("HIST SAME");
        TLine* zl = new TLine(0., 0., 0., ymax);
        zl->SetLineColor(kGray + 2);  zl->SetLineStyle(3);  zl->Draw("SAME");
        TLegend* leg = MakeLegend(0.15, 0.70, 0.60, 0.88);
        if (hDnc) leg->AddEntry(hDnc, kScenLabels[0], "l");
        if (hDpt) leg->AddEntry(hDpt, kScenLabels[1], "l");
        if (hDdc) leg->AddEntry(hDdc, kScenLabels[2], "l");
        if (hDbc) leg->AddEntry(hDbc, kScenLabels[3], "l");
        leg->Draw("SAME");
    }

    // ---- Right panel: <R_proxy> vs eta ----
    c->cd(2);
    gPad->SetLeftMargin(0.14);  gPad->SetBottomMargin(0.13);
    {
        double ymax = 0.;
        for (TProfile* p : {pNC, pPT, pDC, pBC}) {
            if (!p) continue;
            for (int ib = 1; ib <= p->GetNbinsX(); ++ib) {
                double v = std::fabs(p->GetBinContent(ib)) + p->GetBinError(ib);
                if (v > ymax) ymax = v;
            }
        }
        ymax = (ymax < 1.e-6) ? 0.05 : ymax * 1.35;
        if (pNC) {
            pNC->SetTitle(Form("<R_{proxy}> vs #eta_{#Lambda} -- %s;#eta_{#Lambda};<R_{proxy}>", famLabel));
            pNC->GetYaxis()->SetRangeUser(-ymax, ymax);
            pNC->Draw("EP");
        }
        if (pPT) pPT->Draw("EP SAME");
        if (pDC) pDC->Draw("EP SAME");
        if (pBC) pBC->Draw("EP SAME");
        double etaMax = pNC ? pNC->GetXaxis()->GetXmax() : 0.9;
        TLine* zl  = new TLine(-etaMax, 0., etaMax, 0.);
        zl->SetLineColor(kGray + 2);   zl->SetLineStyle(2);  zl->Draw("SAME");
        TLine* eta0 = new TLine(0., -ymax, 0., ymax);
        eta0->SetLineColor(kGray + 1); eta0->SetLineStyle(3); eta0->Draw("SAME");
        TLegend* leg = MakeLegend(0.15, 0.72, 0.50, 0.89);
        if (pNC) leg->AddEntry(pNC, kScenLabels[0], "ep");
        if (pPT) leg->AddEntry(pPT, kScenLabels[1], "ep");
        if (pDC) leg->AddEntry(pDC, kScenLabels[2], "ep");
        if (pBC) leg->AddEntry(pBC, kScenLabels[3], "ep");
        leg->Draw("SAME");
    }

    c->cd(0);
    AddLabel(0.5, 0.995, Form("Ring observable proxy -- %s", famLabel), 0.036, 22);
    WriteCanvas(c, famOut);
    delete c;
}


// ==========================================================================
/**
 * @brief Fig 5 -- <R_proxy> vs Lambda pT for all four scenarios.
 *
 * Uses the "All" eta combined TProfile from each scenario.
 *
 * @param famDir    Family directory in the input file.
 * @param famOut    Output sub-directory.
 * @param famLabel  Short label used in object names.
 */
// ==========================================================================
static void MakeFig5_RingVsPt(TDirectory* famDir, TDirectory* famOut, const char* famLabel)
{
    ScenDirs sd = GetScenDirs(famDir, "All");
    TProfile* pNC = sd.nc ? static_cast<TProfile*>(SafeGet(sd.nc, "pRingProxyVsPt")) : nullptr;
    TProfile* pPT = sd.pt ? static_cast<TProfile*>(SafeGet(sd.pt, "pRingProxyVsPt")) : nullptr;
    TProfile* pDC = sd.dc ? static_cast<TProfile*>(SafeGet(sd.dc, "pRingProxyVsPt")) : nullptr;
    TProfile* pBC = sd.bc ? static_cast<TProfile*>(SafeGet(sd.bc, "pRingProxyVsPt")) : nullptr;
    if (!pNC) return;

    SetHistStyle(pNC, kColNoCuts, kMarkerNC);  SetHistStyle(pPT, kColPtCut,  kMarkerPT);
    SetHistStyle(pDC, kColDcaCut, kMarkerDC);  SetHistStyle(pBC, kColBoth,   kMarkerBC);

    double ymax = 0.;
    for (TProfile* p : {pNC, pPT, pDC, pBC}) {
        if (!p) continue;
        for (int ib = 1; ib <= p->GetNbinsX(); ++ib) {
            double v = std::fabs(p->GetBinContent(ib)) + p->GetBinError(ib);
            if (v > ymax) ymax = v;
        }
    }
    ymax = (ymax < 1.e-6) ? 0.05 : ymax * 1.35;

    TCanvas* c = new TCanvas(Form("c_%s_ringVsPt", famLabel), "", 750, 600);
    c->SetLeftMargin(0.14);  c->SetBottomMargin(0.13);

    pNC->SetTitle(Form("<R_{proxy}> vs #Lambda p_{T} -- %s;p_{T}^{#Lambda} [GeV/c];<R_{proxy}>", famLabel));
    pNC->GetYaxis()->SetRangeUser(-ymax, ymax);
    pNC->Draw("EP");
    if (pPT) pPT->Draw("EP SAME");
    if (pDC) pDC->Draw("EP SAME");
    if (pBC) pBC->Draw("EP SAME");

    double ptMax = pNC->GetXaxis()->GetXmax();
    TLine* zl = new TLine(0., 0., ptMax, 0.);
    zl->SetLineColor(kGray + 2);  zl->SetLineStyle(2);  zl->Draw("SAME");

    TLegend* leg = MakeLegend(0.55, 0.72, 0.88, 0.89);
    leg->AddEntry(pNC, kScenLabels[0], "ep");
    if (pPT) leg->AddEntry(pPT, kScenLabels[1], "ep");
    if (pDC) leg->AddEntry(pDC, kScenLabels[2], "ep");
    if (pBC) leg->AddEntry(pBC, kScenLabels[3], "ep");
    leg->Draw("SAME");

    WriteCanvas(c, famOut);
    delete c;
}


// ==========================================================================
/**
 * @brief Fig 6 -- integrated <R_proxy> bar chart, one panel per eta selection.
 *
 * Three panels side by side: EtaPos | EtaNeg | All.  Each panel shows the
 * single-bin mean <R_proxy> from the pRingProxy TProfile for all four
 * scenarios, drawn as a central-value bar with error caps.
 *
 * @param famDir    Family directory in the input file.
 * @param famOut    Output sub-directory.
 * @param famLabel  Short label used in object names.
 */
// ==========================================================================
static void MakeFig6_IntegratedRing(TDirectory* famDir, TDirectory* famOut, const char* famLabel)
{
    const char* etaSels[3]   = {"EtaPos", "EtaNeg", "All"};
    const char* etaLabels[3] = {"#eta_{#Lambda} > 0", "#eta_{#Lambda} < 0", "All #eta"};

    // Find the global y range across all panels for consistent axes
    double globalMax = 0.;
    for (int ie = 0; ie < 3; ++ie) {
        for (int is = 0; is < 4; ++is) {
            TDirectory* dir = GetScenarioDir(famDir, kScenNames[is], etaSels[ie]);
            if (!dir) continue;
            TProfile* p = static_cast<TProfile*>(SafeGet(dir, "pRingProxy"));
            if (!p) continue;
            double v = std::fabs(p->GetBinContent(1)) + p->GetBinError(1);
            if (v > globalMax) globalMax = v;
        }
    }
    
    // // Enforce a hard minimum scale of 0.02 so we don't zoom in on microscopic 
    // // statistical noise, and increase the margin multiplier to zoom out.
    // globalMax = std::max(0.02, globalMax * 2.5);
    globalMax = globalMax * 1.5;

    TCanvas* c = new TCanvas(Form("c_%s_intRing", famLabel), "", 1050, 500);
    c->Divide(3, 1, 0.004, 0.002);

    for (int ie = 0; ie < 3; ++ie) {
        c->cd(ie + 1);
        gPad->SetLeftMargin(0.18);  gPad->SetBottomMargin(0.16);

        TH1D* hBar = new TH1D(Form("hBar_%s_%s", famLabel, etaSels[ie]),
                               Form("<R_{proxy}> -- %s;Scenario;<R_{proxy}>", etaLabels[ie]),
                               4, 0., 4.);
        hBar->SetDirectory(nullptr); // Critical to make the bars not be owned by the current TDirectory!
                                     // Otherwise, they would dump into the output folder, even if they are
                                     // just temporary histograms!
        hBar->GetYaxis()->SetRangeUser(-globalMax, globalMax);
        hBar->SetStats(0);
        for (int is = 0; is < 4; ++is) hBar->GetXaxis()->SetBinLabel(is + 1, kScenLabels[is]);
        hBar->GetXaxis()->SetLabelSize(0.055);
        hBar->Draw("AXIS"); // Draw frame

        TLine* zl = new TLine(0., 0., 4., 0.);
        zl->SetLineColor(kGray + 2);  zl->SetLineStyle(2);  zl->SetLineWidth(2);  zl->Draw("SAME");

        for (int is = 0; is < 4; ++is) {
            TDirectory* dir = GetScenarioDir(famDir, kScenNames[is], etaSels[ie]);
            if (!dir) continue;
            TProfile* p = static_cast<TProfile*>(SafeGet(dir, "pRingProxy"));
            if (!p) continue;
            
            // Replaces manual lines: fill a single bin and use standard ROOT E1
            TH1D* hPt = new TH1D(Form("hPt6_%s_%d_%d", famLabel, ie, is), "", 4, 0., 4.);
            hPt->SetDirectory(nullptr); // Same fix as above
            hPt->SetBinContent(is + 1, p->GetBinContent(1));
            hPt->SetBinError(is + 1, p->GetBinError(1));
            SetHistStyle(hPt, kScenColors[is], kScenMarkers[is]);
            hPt->SetLineWidth(3);
            hPt->Draw("E1 SAME");
        }
        AddLabel(0.5, 0.96, etaLabels[ie], 0.045, 22);
    }

    c->cd(0);
    AddLabel(0.5, 0.995, Form("Integrated <R_{proxy}> per scenario -- %s", famLabel), 0.036, 22);
    WriteCanvas(c, famOut);
    delete c;
}


// ==========================================================================
/**
 * @brief Fig 7 -- daughter pT distributions before and after the pT cut.
 *
 * Two panels: proton (left) and pion (right).  Pre-cut spectra come from the
 * top-level Kinematics directory; post-cut spectra from pTCutOnly/All inside
 * the family directory.  Both are normalised to the pre-cut integral.
 *
 * @param f         Input ROOT file (for Kinematics).
 * @param famDir    Family directory (for post-cut scenario).
 * @param famOut    Output sub-directory.
 * @param famLabel  Short label used in object names.
 */
// ==========================================================================
static void MakeFig7_DaughterPt(TFile* f, TDirectory* famDir,
                                  TDirectory* famOut, const char* famLabel)
{
    TDirectory* dirKin   = GetDir(f, "Kinematics");
    TDirectory* dirPtAll = GetScenarioDir(famDir, "pTCutOnly", "All");
    if (!dirKin) return;

    TH1D* hKin_p  = static_cast<TH1D*>(SafeGet(dirKin,   "hKin_pT_proton"));
    TH1D* hKin_pi = static_cast<TH1D*>(SafeGet(dirKin,   "hKin_pT_pion"));
    TH1D* hCut_p  = dirPtAll ? static_cast<TH1D*>(SafeGet(dirPtAll, "h1d_pT_proton")) : nullptr;
    TH1D* hCut_pi = dirPtAll ? static_cast<TH1D*>(SafeGet(dirPtAll, "h1d_pT_pion"))   : nullptr;
    if (!hKin_p || !hKin_pi) return;

    TCanvas* c = new TCanvas(Form("c_%s_daughterPt", famLabel), "", 1000, 500);
    c->Divide(2, 1, 0.004, 0.002);

    for (int id = 0; id < 2; ++id) {
        c->cd(id + 1);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.14);  gPad->SetBottomMargin(0.13);

        TH1D* hPre  = (id == 0) ? hKin_p  : hKin_pi;
        TH1D* hPost = (id == 0) ? hCut_p  : hCut_pi;

        TH1D* hPreN = SafeClone(hPre);
        TH1D* hPostN = SafeClone(hPost);

        SetHistStyle(hPreN, kBlack, kMarkerNC);
        if (hPostN) SetHistStyle(hPostN, kColPtCut, kMarkerPT);

        const char* pname = (id == 0) ? "Proton" : "Pion";
        hPreN->SetTitle(Form("%s p_{T} -- %s;p_{T} [GeV/c];Counts", pname, famLabel));
        hPreN->Draw("HIST");
        if (hPostN) hPostN->Draw("HIST SAME");

        TLegend* leg = MakeLegend(0.45, 0.72, 0.88, 0.88);
        leg->AddEntry(hPreN, "Pre-cut (no #eta cut)", "l");
        if (hPostN) leg->AddEntry(hPostN, "After p_{T} cut (with #eta cut)", "l");
        leg->Draw("SAME");

        AddLabel(0.5, 0.93, Form("%s p_{T}", pname), 0.04, 22);
    }

    WriteCanvas(c, famOut);
    delete c;
}


// ==========================================================================
/**
 * @brief Fig 8 -- daughter DCA_xy distributions before and after the DCA cut.
 *
 * Two panels: proton (left) and pion (right).  Structure mirrors Fig 7.
 *
 * @param f         Input ROOT file (for Kinematics).
 * @param famDir    Family directory (for post-cut scenario).
 * @param famOut    Output sub-directory.
 * @param famLabel  Short label used in object names.
 */
// ==========================================================================
static void MakeFig8_DaughterDCA(TFile* f, TDirectory* famDir,
                                   TDirectory* famOut, const char* famLabel)
{
    TDirectory* dirKin    = GetDir(f, "Kinematics");
    TDirectory* dirDcaAll = GetScenarioDir(famDir, "DCACutOnly", "All");
    if (!dirKin) return;

    TH1D* hKin_p  = static_cast<TH1D*>(SafeGet(dirKin,    "hKin_DCA_proton"));
    TH1D* hKin_pi = static_cast<TH1D*>(SafeGet(dirKin,    "hKin_DCA_pion"));
    TH1D* hCut_p  = dirDcaAll ? static_cast<TH1D*>(SafeGet(dirDcaAll, "h1d_DCA_proton")) : nullptr;
    TH1D* hCut_pi = dirDcaAll ? static_cast<TH1D*>(SafeGet(dirDcaAll, "h1d_DCA_pion"))   : nullptr;
    if (!hKin_p || !hKin_pi) return;

    TCanvas* c = new TCanvas(Form("c_%s_daughterDCA", famLabel), "", 1000, 500);
    c->Divide(2, 1, 0.004, 0.002);

    for (int id = 0; id < 2; ++id) {
        c->cd(id + 1);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.14);  gPad->SetBottomMargin(0.13);

        TH1D* hPre  = (id == 0) ? hKin_p  : hKin_pi;
        TH1D* hPost = (id == 0) ? hCut_p  : hCut_pi;

        TH1D* hPreN = SafeClone(hPre);
        TH1D* hPostN = SafeClone(hPost);

        SetHistStyle(hPreN, kBlack, kMarkerNC);
        if (hPostN) SetHistStyle(hPostN, kColDcaCut, kMarkerDC);

        const char* pname = (id == 0) ? "Proton" : "Pion";
        hPreN->SetTitle(Form("%s DCA_{xy} -- %s;DCA_{xy} [cm];Counts", pname, famLabel));
        hPreN->Draw("HIST");
        if (hPostN) hPostN->Draw("HIST SAME");

        TLegend* leg = MakeLegend(0.45, 0.72, 0.88, 0.88);
        leg->AddEntry(hPreN, "Pre-cut (no #eta cut)", "l");
        if (hPostN) leg->AddEntry(hPostN, "After DCA cut (with #eta cut)", "l");
        leg->Draw("SAME");

        AddLabel(0.5, 0.93, Form("%s DCA_{xy}", pname), 0.04, 22);
    }

    WriteCanvas(c, famOut);
    delete c;
}


// ==========================================================================
/**
 * @brief Fig 9 -- eta>0 minus eta<0 antisymmetric component plots.
 *
 * Three panels:
 *   Panel 1: phi* difference (EtaPos - EtaNeg, normalised) for all 4 scenarios.
 *   Panel 2: cos(theta*) difference for all 4 scenarios.
 *   Panel 3: <R_proxy>(EtaPos) - <R_proxy>(EtaNeg) bar chart for all 4 scenarios.
 *
 * Bug fix relative to the old macro: y-axis ranges are computed from the data
 * before the first Draw call, so no placeholder range is needed.
 *
 * @param famDir    Family directory in the input file.
 * @param famOut    Output sub-directory.
 * @param famLabel  Short label used in object names.
 */
// ==========================================================================
static void MakeFig9_EtaDiff(TDirectory* famDir, TDirectory* famOut, const char* famLabel)
{
    TCanvas* c = new TCanvas(Form("c_%s_etaDiff", famLabel), "", 1200, 500);
    c->Divide(3, 1, 0.003, 0.002);

    auto MakeDiff = [&](const char* scenName, const char* histName) -> TH1D* {
        TDirectory* dP = GetScenarioDir(famDir, scenName, "EtaPos");
        TDirectory* dN = GetScenarioDir(famDir, scenName, "EtaNeg");
        if (!dP || !dN) return nullptr;
        TH1D* hP = static_cast<TH1D*>(SafeGet(dP, histName));
        TH1D* hN = static_cast<TH1D*>(SafeGet(dN, histName));
        if (!hP || !hN) return nullptr;
        
        TH1D* hPn = static_cast<TH1D*>(hP->Clone(Form("%s_diff_%d", histName, gCloneIdx++)));
        hPn->SetDirectory(nullptr); // Same fix as above
        hPn->Add(hN, -1.);
        return hPn;
    };

    // Helper: draw a set of difference histograms with a shared y range
    auto DrawDiffPanel = [&](const char* histName, const char* title, const char* xTitle) {
        std::vector<TH1D*> diffs;
        for (int is = 0; is < 4; ++is) {
            TH1D* h = MakeDiff(kScenNames[is], histName);
            if (h) { SetHistStyle(h, kScenColors[is], kScenMarkers[is]); }
            diffs.push_back(h);
        }
        
        double ymax = 0.;
        for (auto* h : diffs) {
            if (!h) continue;
            double lm = std::max(std::fabs(h->GetMinimum()), h->GetMaximum());
            if (lm > ymax) ymax = lm;
        }
        ymax = (ymax < 1.e-8) ? 0.005 : ymax * 1.4;

        bool first = true;
        for (int is = 0; is < 4; ++is) {
            if (!diffs[is]) continue;
            if (first) {
                diffs[is]->SetTitle(Form("%s -- %s;%s;#Delta(Counts)", title, famLabel, xTitle));
                diffs[is]->GetYaxis()->SetRangeUser(-ymax, ymax);
                diffs[is]->Draw("HIST");
                first = false;
            } else {
                diffs[is]->Draw("HIST SAME");
            }
        }
        TLine* zl = new TLine(diffs[0] ? diffs[0]->GetXaxis()->GetXmin() : -1., 0.,
                               diffs[0] ? diffs[0]->GetXaxis()->GetXmax() :  1., 0.);
        zl->SetLineColor(kGray + 2);  zl->SetLineStyle(2);  zl->Draw("SAME");

        TLegend* leg = MakeLegend(0.15, 0.70, 0.60, 0.90);
        for (int is = 0; is < 4; ++is) {
            if (diffs[is]) leg->AddEntry(diffs[is], kScenLabels[is], "l");
        }
        leg->Draw("SAME");
        // Removed histogram deletion here so the canvas successfully writes!
    };

    // ---- Panel 1: phi* difference ----
    c->cd(1);  gPad->SetLeftMargin(0.14);
    DrawDiffPanel("h1d_phi", "#phi* (#eta>0) - (#eta<0)", "#phi* [rad]");

    // ---- Panel 2: cos(theta*) difference ----
    c->cd(2);  gPad->SetLeftMargin(0.14);
    DrawDiffPanel("h1d_cosTheta", "cos(#theta*) (#eta>0) - (#eta<0)", "cos(#theta*)");

    // ---- Panel 3: <R_proxy> EtaPos minus EtaNeg bar chart ----
    c->cd(3);
    gPad->SetLeftMargin(0.14);  gPad->SetBottomMargin(0.16);

    double vals[4] = {0.}, errs[4] = {0.};
    for (int is = 0; is < 4; ++is) {
        TDirectory* dP = GetScenarioDir(famDir, kScenNames[is], "EtaPos");
        TDirectory* dN = GetScenarioDir(famDir, kScenNames[is], "EtaNeg");
        if (!dP || !dN) continue;
        TProfile* pP = static_cast<TProfile*>(SafeGet(dP, "pRingProxy"));
        TProfile* pN = static_cast<TProfile*>(SafeGet(dN, "pRingProxy"));
        if (!pP || !pN) continue;
        vals[is] = pP->GetBinContent(1) - pN->GetBinContent(1);
        errs[is] = std::sqrt(pP->GetBinError(1)*pP->GetBinError(1) +
                              pN->GetBinError(1)*pN->GetBinError(1));
    }

    double ymax3 = 0.;
    for (int is = 0; is < 4; ++is) {
        double v = std::fabs(vals[is]) + errs[is];
        if (v > ymax3) ymax3 = v;
    }
    ymax3 = (ymax3 < 1.e-8) ? 0.01 : ymax3 * 1.5;

    TH1D* hBar = new TH1D(Form("hBar3_%s", famLabel),
                            Form("<R_{proxy}>_{#eta>0} - <R_{proxy}>_{#eta<0} -- %s;"
                                 "Scenario;#Delta<R_{proxy}>", famLabel),
                            4, 0., 4.);
    hBar->SetDirectory(nullptr); // Same fix as above
    hBar->SetStats(0);
    hBar->GetYaxis()->SetRangeUser(-ymax3, ymax3);
    for (int is = 0; is < 4; ++is) hBar->GetXaxis()->SetBinLabel(is + 1, kScenLabels[is]);
    hBar->GetXaxis()->SetLabelSize(0.055);
    hBar->Draw("AXIS"); // Draw frame

    TLine* zl3 = new TLine(0., 0., 4., 0.);
    zl3->SetLineColor(kGray + 2);  zl3->SetLineStyle(2);  zl3->SetLineWidth(2);  zl3->Draw("SAME");

    for (int is = 0; is < 4; ++is) {
        // Replaces manual lines: fill a single bin and use standard ROOT E1
        TH1D* hPt = new TH1D(Form("hPt9_%s_%d", famLabel, is), "", 4, 0., 4.);
        hPt->SetDirectory(nullptr); // Same fix as above
        hPt->SetBinContent(is + 1, vals[is]);
        hPt->SetBinError(is + 1, errs[is]);
        SetHistStyle(hPt, kScenColors[is], kScenMarkers[is]);
        hPt->SetLineWidth(3);
        hPt->Draw("E1 SAME");
    }

    TLegend* leg3 = MakeLegend(0.15, 0.70, 0.60, 0.90);
    for (int is = 0; is < 4; ++is) {
        TLine* dummy = new TLine(); dummy->SetLineColor(kScenColors[is]); dummy->SetLineWidth(3);
        leg3->AddEntry(dummy, kScenLabels[is], "l");
    }
    leg3->Draw("SAME");

    c->cd(0);
    AddLabel(0.5, 0.995, Form("#eta>0 minus #eta<0 antisymmetric component -- %s", famLabel), 0.034, 22);
    WriteCanvas(c, famOut);
    delete c;
}

// ==========================================================================
/**
 * @brief Fig 10 -- <R_proxy_jet> vs Lambda eta and Jet eta.
 *
 * 2 panels: 
 * Left:  <R_proxy_jet> vs \eta_\Lambda
 * Right: <R_proxy_jet> vs \eta_{jet}
 *
 * @param famDir    Family directory in the input file.
 * @param famOut    Output sub-directory.
 * @param famLabel  Short label used in object names.
 */
// ==========================================================================
static void MakeFig10_RingProxyJetVsEta(TDirectory* famDir, TDirectory* famOut, const char* famLabel)
{
    ScenDirs sd = GetScenDirs(famDir, "All");
    
    // --- Retrieve <R_proxy_jet> vs Lambda eta ---
    TProfile* pNCL = sd.nc ? static_cast<TProfile*>(SafeGet(sd.nc, "pRingProxyJetVsEta")) : nullptr;
    TProfile* pPTL = sd.pt ? static_cast<TProfile*>(SafeGet(sd.pt, "pRingProxyJetVsEta")) : nullptr;
    TProfile* pDCL = sd.dc ? static_cast<TProfile*>(SafeGet(sd.dc, "pRingProxyJetVsEta")) : nullptr;
    TProfile* pBCL = sd.bc ? static_cast<TProfile*>(SafeGet(sd.bc, "pRingProxyJetVsEta")) : nullptr;

    // --- Retrieve <R_proxy_jet> vs Jet eta ---
    TProfile* pNCJ = sd.nc ? static_cast<TProfile*>(SafeGet(sd.nc, "pRingProxyJetVsEtaJet")) : nullptr;
    TProfile* pPTJ = sd.pt ? static_cast<TProfile*>(SafeGet(sd.pt, "pRingProxyJetVsEtaJet")) : nullptr;
    TProfile* pDCJ = sd.dc ? static_cast<TProfile*>(SafeGet(sd.dc, "pRingProxyJetVsEtaJet")) : nullptr;
    TProfile* pBCJ = sd.bc ? static_cast<TProfile*>(SafeGet(sd.bc, "pRingProxyJetVsEtaJet")) : nullptr;

    if (!pNCL && !pNCJ) return;

    SetHistStyle(pNCL, kColNoCuts, kMarkerNC);  SetHistStyle(pPTL, kColPtCut,  kMarkerPT);
    SetHistStyle(pDCL, kColDcaCut, kMarkerDC);  SetHistStyle(pBCL, kColBoth,   kMarkerBC);

    SetHistStyle(pNCJ, kColNoCuts, kMarkerNC);  SetHistStyle(pPTJ, kColPtCut,  kMarkerPT);
    SetHistStyle(pDCJ, kColDcaCut, kMarkerDC);  SetHistStyle(pBCJ, kColBoth,   kMarkerBC);

    TCanvas* c = new TCanvas(Form("c_%s_ringProxyJetVsEta", famLabel), "", 1200, 550);
    c->Divide(2, 1, 0.004, 0.002);

    auto drawPanel = [](TProfile* pNC, TProfile* pPT, TProfile* pDC, TProfile* pBC, const char* xTitle) {
        double ymax = 0.;
        for (TProfile* p : {pNC, pPT, pDC, pBC}) {
            if (!p) continue;
            for (int ib = 1; ib <= p->GetNbinsX(); ++ib) {
                double v = std::fabs(p->GetBinContent(ib)) + p->GetBinError(ib);
                if (v > ymax) ymax = v;
            }
        }
        ymax = (ymax < 1.e-6) ? 0.05 : ymax * 1.35;

        if (pNC) {
            pNC->SetTitle(Form(";<R_{proxy}^{jet}>"));
            pNC->GetXaxis()->SetTitle(xTitle);
            pNC->GetYaxis()->SetRangeUser(-ymax, ymax);
            pNC->Draw("EP");
        }
        if (pPT) pPT->Draw("EP SAME");
        if (pDC) pDC->Draw("EP SAME");
        if (pBC) pBC->Draw("EP SAME");

        double xMax = pNC ? pNC->GetXaxis()->GetXmax() : 1.0;
        TLine* zl  = new TLine(-xMax, 0., xMax, 0.);
        zl->SetLineColor(kGray + 2);   zl->SetLineStyle(2);  zl->Draw("SAME");
        TLine* eta0 = new TLine(0., -ymax, 0., ymax);
        eta0->SetLineColor(kGray + 1); eta0->SetLineStyle(3); eta0->Draw("SAME");

        TLegend* leg = MakeLegend(0.15, 0.72, 0.50, 0.89);
        if (pNC) leg->AddEntry(pNC, kScenLabels[0], "ep");
        if (pPT) leg->AddEntry(pPT, kScenLabels[1], "ep");
        if (pDC) leg->AddEntry(pDC, kScenLabels[2], "ep");
        if (pBC) leg->AddEntry(pBC, kScenLabels[3], "ep");
        leg->Draw("SAME");
    };

    c->cd(1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
    drawPanel(pNCL, pPTL, pDCL, pBCL, "#eta_{#Lambda}");
    AddLabel(0.5, 0.94, "Dependence on #Lambda pseudorapidity", 0.045, 22);

    c->cd(2); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
    drawPanel(pNCJ, pPTJ, pDCJ, pBCJ, "#eta_{jet}");
    AddLabel(0.5, 0.94, "Dependence on Jet pseudorapidity", 0.045, 22);

    c->cd(0);
    AddLabel(0.5, 0.995, Form("Jet Ring observable proxy vs Eta -- %s", famLabel), 0.036, 22);
    WriteCanvas(c, famOut);
    delete c;
}


// ==========================================================================
/**
 * @brief Fig 11 -- <R_proxy_jet> vs Lambda pT for all four scenarios.
 */
// ==========================================================================
static void MakeFig11_ringJetVsPt(TDirectory* famDir, TDirectory* famOut, const char* famLabel)
{
    ScenDirs sd = GetScenDirs(famDir, "All");
    TProfile* pNC = sd.nc ? static_cast<TProfile*>(SafeGet(sd.nc, "pRingProxyJetVsPt")) : nullptr;
    TProfile* pPT = sd.pt ? static_cast<TProfile*>(SafeGet(sd.pt, "pRingProxyJetVsPt")) : nullptr;
    TProfile* pDC = sd.dc ? static_cast<TProfile*>(SafeGet(sd.dc, "pRingProxyJetVsPt")) : nullptr;
    TProfile* pBC = sd.bc ? static_cast<TProfile*>(SafeGet(sd.bc, "pRingProxyJetVsPt")) : nullptr;
    if (!pNC) return;

    SetHistStyle(pNC, kColNoCuts, kMarkerNC);  SetHistStyle(pPT, kColPtCut,  kMarkerPT);
    SetHistStyle(pDC, kColDcaCut, kMarkerDC);  SetHistStyle(pBC, kColBoth,   kMarkerBC);

    double ymax = 0.;
    for (TProfile* p : {pNC, pPT, pDC, pBC}) {
        if (!p) continue;
        for (int ib = 1; ib <= p->GetNbinsX(); ++ib) {
            double v = std::fabs(p->GetBinContent(ib)) + p->GetBinError(ib);
            if (v > ymax) ymax = v;
        }
    }
    ymax = (ymax < 1.e-6) ? 0.05 : ymax * 1.35;

    TCanvas* c = new TCanvas(Form("c_%s_ringJetVsPt", famLabel), "", 750, 600);
    c->SetLeftMargin(0.14);  c->SetBottomMargin(0.13);

    pNC->SetTitle(Form("<R_{proxy}^{jet}> vs #Lambda p_{T} -- %s;p_{T}^{#Lambda} [GeV/c];<R_{proxy}^{jet}>", famLabel));
    pNC->GetYaxis()->SetRangeUser(-ymax, ymax);
    pNC->Draw("EP");
    if (pPT) pPT->Draw("EP SAME");
    if (pDC) pDC->Draw("EP SAME");
    if (pBC) pBC->Draw("EP SAME");

    double ptMax = pNC->GetXaxis()->GetXmax();
    TLine* zl = new TLine(0., 0., ptMax, 0.);
    zl->SetLineColor(kGray + 2);  zl->SetLineStyle(2);  zl->Draw("SAME");

    TLegend* leg = MakeLegend(0.55, 0.72, 0.88, 0.89);
    leg->AddEntry(pNC, kScenLabels[0], "ep");
    if (pPT) leg->AddEntry(pPT, kScenLabels[1], "ep");
    if (pDC) leg->AddEntry(pDC, kScenLabels[2], "ep");
    if (pBC) leg->AddEntry(pBC, kScenLabels[3], "ep");
    leg->Draw("SAME");

    WriteCanvas(c, famOut);
    delete c;
}


// ==========================================================================
/**
 * @brief Fig 12 -- integrated <R_proxy_jet> bar chart, one panel per eta selection.
 */
// ==========================================================================
static void MakeFig12_IntegratedRingJet(TDirectory* famDir, TDirectory* famOut, const char* famLabel)
{
    const char* etaSels[3]   = {"EtaPos", "EtaNeg", "All"};
    const char* etaLabels[3] = {"#eta_{#Lambda} > 0", "#eta_{#Lambda} < 0", "All #eta"};

    // Find the global y range across all panels for consistent axes
    double globalMax = 0.;
    for (int ie = 0; ie < 3; ++ie) {
        for (int is = 0; is < 4; ++is) {
            TDirectory* dir = GetScenarioDir(famDir, kScenNames[is], etaSels[ie]);
            if (!dir) continue;
            TProfile* p = static_cast<TProfile*>(SafeGet(dir, "pRingProxyJet"));
            if (!p) continue;
            double v = std::fabs(p->GetBinContent(1)) + p->GetBinError(1);
            if (v > globalMax) globalMax = v;
        }
    }
    
    // // Enforce a hard minimum scale of 0.02 so we don't zoom in on microscopic 
    // // statistical noise, and increase the margin multiplier to zoom out.
    // globalMax = std::max(0.02, globalMax * 2.5);
    globalMax = globalMax * 1.2;

    TCanvas* c = new TCanvas(Form("c_%s_intRingJet", famLabel), "", 1050, 500);
    c->Divide(3, 1, 0.004, 0.002);

    for (int ie = 0; ie < 3; ++ie) {
        c->cd(ie + 1);
        gPad->SetLeftMargin(0.18);  gPad->SetBottomMargin(0.16);

        TH1D* hBar = new TH1D(Form("hBarJet_%s_%s", famLabel, etaSels[ie]),
                               Form("<R_{proxy}^{jet}> -- %s;Scenario;<R_{proxy}^{jet}>", etaLabels[ie]),
                               4, 0., 4.);
        hBar->SetDirectory(nullptr);
        hBar->GetYaxis()->SetRangeUser(-globalMax, globalMax);
        hBar->SetStats(0);
        for (int is = 0; is < 4; ++is) hBar->GetXaxis()->SetBinLabel(is + 1, kScenLabels[is]);
        hBar->GetXaxis()->SetLabelSize(0.055);
        hBar->Draw("AXIS"); // Draw frame

        TLine* zl = new TLine(0., 0., 4., 0.);
        zl->SetLineColor(kGray + 2);  zl->SetLineStyle(2);  zl->SetLineWidth(2);  zl->Draw("SAME");

        for (int is = 0; is < 4; ++is) {
            TDirectory* dir = GetScenarioDir(famDir, kScenNames[is], etaSels[ie]);
            if (!dir) continue;
            TProfile* p = static_cast<TProfile*>(SafeGet(dir, "pRingProxyJet"));
            if (!p) continue;
            
            TH1D* hPt = new TH1D(Form("hPtJet12_%s_%d_%d", famLabel, ie, is), "", 4, 0., 4.);
            hPt->SetDirectory(nullptr); 
            hPt->SetBinContent(is + 1, p->GetBinContent(1));
            hPt->SetBinError(is + 1, p->GetBinError(1));
            SetHistStyle(hPt, kScenColors[is], kScenMarkers[is]);
            hPt->SetLineWidth(3);
            hPt->Draw("E1 SAME");
        }
        AddLabel(0.5, 0.96, etaLabels[ie], 0.045, 22);
    }

    c->cd(0);
    AddLabel(0.5, 0.995, Form("Integrated <R_{proxy}^{jet}> per scenario -- %s", famLabel), 0.036, 22);
    WriteCanvas(c, famOut);
    delete c;
}

// ==========================================================================
/**
 * @brief Fig 13 -- 2-panel random jet ring observable proxy figure.
 *
 * Left panel: distributions of R_proxyJet for all four scenarios ("All" eta),
 * normalised to unit area.  Right panel: profile <R_proxyJet> vs Lambda eta for
 * all four scenarios.
 *
 * @param famDir    Family directory in the input file.
 * @param famOut    Output sub-directory.
 * @param famLabel  Short label used in object names.
 */
// ==========================================================================
static void MakeFig13_RingProxyJet(TDirectory* famDir, TDirectory* famOut, const char* famLabel)
{
    // --- Retrieve distributions (h1d_ringProxyJet) from "All" directories ---
    TH1D* hDistNC = nullptr;  TH1D* hDistPT = nullptr;
    TH1D* hDistDC = nullptr;  TH1D* hDistBC = nullptr;
    {
        TDirectory* dNC = GetScenarioDir(famDir, "NoCuts",     "All");
        TDirectory* dPT = GetScenarioDir(famDir, "pTCutOnly",  "All");
        TDirectory* dDC = GetScenarioDir(famDir, "DCACutOnly", "All");
        TDirectory* dBC = GetScenarioDir(famDir, "BothCuts",   "All");
        if (dNC) hDistNC = static_cast<TH1D*>(SafeGet(dNC, "h1d_ringProxyJet"));
        if (dPT) hDistPT = static_cast<TH1D*>(SafeGet(dPT, "h1d_ringProxyJet"));
        if (dDC) hDistDC = static_cast<TH1D*>(SafeGet(dDC, "h1d_ringProxyJet"));
        if (dBC) hDistBC = static_cast<TH1D*>(SafeGet(dBC, "h1d_ringProxyJet"));
    }

    // --- Retrieve <R_proxyJet> vs eta profiles from "All" ---
    ScenDirs sd = GetScenDirs(famDir, "All");
    TProfile* pNC = sd.nc ? static_cast<TProfile*>(SafeGet(sd.nc, "pRingProxyJetVsEta")) : nullptr;
    TProfile* pPT = sd.pt ? static_cast<TProfile*>(SafeGet(sd.pt, "pRingProxyJetVsEta")) : nullptr;
    TProfile* pDC = sd.dc ? static_cast<TProfile*>(SafeGet(sd.dc, "pRingProxyJetVsEta")) : nullptr;
    TProfile* pBC = sd.bc ? static_cast<TProfile*>(SafeGet(sd.bc, "pRingProxyJetVsEta")) : nullptr;
    if (!pNC && !hDistNC) return;

    TH1D* hDnc = SafeClone(hDistNC);
    TH1D* hDpt = SafeClone(hDistPT);
    TH1D* hDdc = SafeClone(hDistDC);
    TH1D* hDbc = SafeClone(hDistBC);
    
    SetHistStyle(hDnc, kColNoCuts, kMarkerNC);  SetHistStyle(hDpt, kColPtCut,  kMarkerPT);
    SetHistStyle(hDdc, kColDcaCut, kMarkerDC);  SetHistStyle(hDbc, kColBoth,   kMarkerBC);

    SetHistStyle(pNC, kColNoCuts, kMarkerNC);   SetHistStyle(pPT, kColPtCut,  kMarkerPT);
    SetHistStyle(pDC, kColDcaCut, kMarkerDC);   SetHistStyle(pBC, kColBoth,   kMarkerBC);

    TCanvas* c = new TCanvas(Form("c_%s_ringProxyJetDist", famLabel), "", 1200, 550);
    c->Divide(2, 1, 0.004, 0.002);

    // ---- Left panel: R_proxyJet distribution ----
    c->cd(1);
    gPad->SetLeftMargin(0.14);  gPad->SetBottomMargin(0.13);
    {
        double ymax = 0.;
        for (TH1D* h : {hDnc, hDpt, hDdc, hDbc}) { if (h && h->GetMaximum() > ymax) ymax = h->GetMaximum(); }
        ymax *= 1.3;
        if (hDnc) {
            hDnc->SetTitle(Form("R_{proxy}^{jet} distribution -- %s;R_{proxy}^{jet};Counts", famLabel));
            hDnc->GetYaxis()->SetRangeUser(0., ymax);
            hDnc->Draw("HIST");
        }
        if (hDpt) hDpt->Draw("HIST SAME");
        if (hDdc) hDdc->Draw("HIST SAME");
        if (hDbc) hDbc->Draw("HIST SAME");
        TLine* zl = new TLine(0., 0., 0., ymax);
        zl->SetLineColor(kGray + 2);  zl->SetLineStyle(3);  zl->Draw("SAME");
        TLegend* leg = MakeLegend(0.15, 0.70, 0.60, 0.88);
        if (hDnc) leg->AddEntry(hDnc, kScenLabels[0], "l");
        if (hDpt) leg->AddEntry(hDpt, kScenLabels[1], "l");
        if (hDdc) leg->AddEntry(hDdc, kScenLabels[2], "l");
        if (hDbc) leg->AddEntry(hDbc, kScenLabels[3], "l");
        leg->Draw("SAME");
    }

    // ---- Right panel: <R_proxyJet> vs eta ----
    c->cd(2);
    gPad->SetLeftMargin(0.14);  gPad->SetBottomMargin(0.13);
    {
        double ymax = 0.;
        for (TProfile* p : {pNC, pPT, pDC, pBC}) {
            if (!p) continue;
            for (int ib = 1; ib <= p->GetNbinsX(); ++ib) {
                double v = std::fabs(p->GetBinContent(ib)) + p->GetBinError(ib);
                if (v > ymax) ymax = v;
            }
        }
        ymax = (ymax < 1.e-6) ? 0.05 : ymax * 1.35;
        if (pNC) {
            pNC->SetTitle(Form("<R_{proxy}^{jet}> vs #eta_{#Lambda} -- %s;#eta_{#Lambda};<R_{proxy}^{jet}>", famLabel));
            pNC->GetYaxis()->SetRangeUser(-ymax, ymax);
            pNC->Draw("EP");
        }
        if (pPT) pPT->Draw("EP SAME");
        if (pDC) pDC->Draw("EP SAME");
        if (pBC) pBC->Draw("EP SAME");
        double etaMax = pNC ? pNC->GetXaxis()->GetXmax() : 0.9;
        TLine* zl  = new TLine(-etaMax, 0., etaMax, 0.);
        zl->SetLineColor(kGray + 2);   zl->SetLineStyle(2);  zl->Draw("SAME");
        TLine* eta0 = new TLine(0., -ymax, 0., ymax);
        eta0->SetLineColor(kGray + 1); eta0->SetLineStyle(3); eta0->Draw("SAME");
        TLegend* leg = MakeLegend(0.15, 0.72, 0.50, 0.89);
        if (pNC) leg->AddEntry(pNC, kScenLabels[0], "ep");
        if (pPT) leg->AddEntry(pPT, kScenLabels[1], "ep");
        if (pDC) leg->AddEntry(pDC, kScenLabels[2], "ep");
        if (pBC) leg->AddEntry(pBC, kScenLabels[3], "ep");
        leg->Draw("SAME");
    }

    c->cd(0);
    AddLabel(0.5, 0.995, Form("Jet Ring observable proxy -- %s", famLabel), 0.036, 22);
    WriteCanvas(c, famOut);
    delete c;
}


// ==========================================================================
/**
 * @brief Lambda kinematics figure -- drawn once, shared across families.
 *
 * Six panels (2 rows x 3 columns):
 *   Row 1: Lambda pT | Lambda eta | Lambda phi
 *   Row 2: Transverse decay radius | Proton pT (pre-cut) | Pion pT (pre-cut)
 *
 * All histograms come from the top-level Kinematics directory.
 *
 * @param f     Input ROOT file.
 * @param fout  Output ROOT file (canvas written at top level).
 */
// ==========================================================================
static void MakeFigKin(TFile* f, TFile* fout)
{
    TDirectory* dirKin = GetDir(f, "Kinematics");
    if (!dirKin) return;

    TH1D* hpT  = static_cast<TH1D*>(SafeGet(dirKin, "hKin_pT_lambda"));
    TH1D* heta = static_cast<TH1D*>(SafeGet(dirKin, "hKin_eta_lambda"));
    TH1D* hphi = static_cast<TH1D*>(SafeGet(dirKin, "hKin_phi_lambda"));
    TH1D* hR   = static_cast<TH1D*>(SafeGet(dirKin, "hKin_decayR"));
    TH1D* hpP  = static_cast<TH1D*>(SafeGet(dirKin, "hKin_pT_proton"));
    TH1D* hpPi = static_cast<TH1D*>(SafeGet(dirKin, "hKin_pT_pion"));

    TCanvas* c = new TCanvas("c_lamKin", "#Lambda kinematics (generation level)", 1200, 700);
    c->Divide(3, 2, 0.003, 0.003);

    struct PadInfo { TH1D* h; bool logy; };
    PadInfo pads[6] = {{hpT, true}, {heta, false}, {hphi, false},
                       {hR,  true}, {hpP,  true},  {hpPi, true}};

    for (int i = 0; i < 6; ++i) {
        c->cd(i + 1);
        if (pads[i].logy) gPad->SetLogy();
        gPad->SetLeftMargin(0.15);  gPad->SetBottomMargin(0.15);
        if (!pads[i].h) continue;
        pads[i].h->SetStats(0);
        pads[i].h->SetLineColor(kBlack);  pads[i].h->SetLineWidth(2);
        pads[i].h->Draw("HIST");
    }

    fout->cd();
    c->Write();
    delete c;
}


// ==========================================================================
/**
 * @brief Orchestrates all figure production for helicityEfficiencyToyModel output.
 *
 * @details
 * Opens @p inputFile (a ROOT file produced by helicityEfficiencyToyModel.cxx),
 * creates a companion plots.root in the same directory, and writes all
 * diagnostic TCanvas objects into it under per-family sub-directories.
 *
 * @param inputFile  Path to the generator output ROOT file. Output will be
                     written to the folder where this lies in.
 */
// ==========================================================================
void plotHelicityEfficiency(const char* inputFile = "helicityEffOutput.root")
{
    gROOT->ProcessLine("gErrorIgnoreLevel = kWarning;");

    TString outDir = gSystem->DirName(inputFile);

    printf("\n========================================================\n");
    printf("  plotHelicityEfficiency\n");
    printf("  Input : %s\n", inputFile);
    printf("  Output folder : %s\n", outDir.Data());
    printf("========================================================\n\n");

    // Output ROOT file placed next to the input file
    TString outRootPath = outDir + "/plots.root";
    TFile* fout = new TFile(outRootPath, "RECREATE");
    if (!fout || fout->IsZombie()) {
        printf("ERROR: Cannot create output file '%s'. Aborting.\n", outRootPath.Data());
        return;
    }

    TFile* f = TFile::Open(inputFile, "READ");
    if (!f || f->IsZombie()) {
        printf("ERROR: Cannot open '%s'. Aborting.\n", inputFile);
        fout->Close();
        return;
    }

    // Global style
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetTitleSize(0.05, "XYZ");
    gStyle->SetLabelSize(0.045, "XYZ");
    gStyle->SetTitleOffset(1.1, "Y");

    // Per-family figures
    const char* famNames[2] = {"WithEtaGate", "WithoutEtaGate"};
    for (const char* famName : famNames) {
        TDirectory* famDir = GetDir(f, famName);
        if (!famDir) { printf("WARNING: family '%s' not found -- skipping.\n", famName); continue; }

        // Create the output sub-directory for this family
        TDirectory* famOut = fout->mkdir(famName);

        printf("Producing figures for family: %s\n", famName);

        // Figures 1-3: one canvas per eta selection
        for (const char* eta : {"EtaPos", "EtaNeg", "All"}) {
            MakeFig1_2DMaps  (famDir, famOut, famName, eta);
            MakeFig2_CosTheta(famDir, famOut, famName, eta);
            MakeFig3_PhiStar (famDir, famOut, famName, eta);
        }

        // Figures 4-9: one canvas each (use "All" eta combined internally)
        MakeFig4_RingProxy    (famDir, famOut, famName);
        MakeFig5_RingVsPt     (famDir, famOut, famName);
        MakeFig6_IntegratedRing(famDir, famOut, famName);
        MakeFig7_DaughterPt   (f, famDir, famOut, famName);
        MakeFig8_DaughterDCA  (f, famDir, famOut, famName);
        MakeFig9_EtaDiff      (famDir, famOut, famName);

        // Random jet direction proxy figures:
        // (numeration got confusing, but the plots are now saved in the same order for the \hat z and \hat jet proxies)
        MakeFig13_RingProxyJet      (famDir, famOut, famName);
        MakeFig11_ringJetVsPt       (famDir, famOut, famName);
        MakeFig10_RingProxyJetVsEta (famDir, famOut, famName);
        MakeFig12_IntegratedRingJet (famDir, famOut, famName);
    }

    // Lambda kinematics -- shared, written at the fout top level
    printf("Producing kinematics figure...\n");
    MakeFigKin(f, fout);

    f->Close();
    fout->Write("", TObject::kOverwrite);
    fout->Close();

    printf("\n========================================================\n");
    printf("  Output written to: %s\n", outRootPath.Data());
    printf("========================================================\n\n");
}
// end of plotHelicityEfficiency.cxx