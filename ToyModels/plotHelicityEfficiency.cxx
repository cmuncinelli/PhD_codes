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
// New for the non-SEM error estimators of integrated ring observables:
#include <TGraphErrors.h>
#include <TMarker.h>

#include <TArrow.h>
#include <TProfile2D.h>

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
 * @brief Draws one vector-field panel: pPstarZ as a COLZ background with
 *        block-averaged transverse (px, py) arrows overlaid.
 *
 * @details
 * The colour background encodes <p*_z> on a diverging palette (set by the
 * caller via gStyle->SetPalette) with a symmetric Z range derived from the
 * data.  Arrows represent the mean transverse proton direction in the Lambda
 * rest frame, averaged over arrowBlockSize x arrowBlockSize tiles of bins.
 *
 * Arrow scaling is robust against edge outliers: all valid tile magnitudes
 * are collected, sorted, and the scalePercentile-th percentile is used as
 * the normalisation reference (scaleRef).  The longest "typical" arrow spans
 * 0.45 * tileWidth.  Tiles whose averaged magnitude exceeds scaleRef have
 * their arrow length capped at 0.45 * tileWidth, so one noisy edge tile
 * cannot compress the rest of the field into invisibility.
 *
 * Tiles with fewer than minEntries occupied bins inside them, or whose
 * averaged transverse magnitude is below 5% of scaleRef, are suppressed.
 *
 * Must be called while the target TPad is current (gPad).
 *
 * @param hX              TProfile2D of <p*_x> vs (px_lam, py_lam).
 * @param hY              TProfile2D of <p*_y> vs (px_lam, py_lam).
 * @param hZ              TProfile2D of <p*_z> vs (px_lam, py_lam) [colormap].
 * @param title           Histogram title string; pass nullptr to suppress.
 * @param minEntries      Minimum per-bin entry count to include in a tile.
 *                        Default: 50.
 * @param arrowBlockSize  Side length of the averaging tile in bins.
 *                        Default: 4 (4x4 tiles, 25 arrows/axis for 20 bins).
 *                        Valid choices for 20 bins: 1, 2, 4, 5, 10, 20.
 * @param scalePercentile Fraction of tile magnitudes used as the length
 *                        normalisation reference.  0.90 (default) means the
 *                        90th-percentile magnitude maps to the full arrow
 *                        length; the top 10% are drawn capped at that length.
 *                        Tighten toward 0.75 if outliers are still visible;
 *                        loosen toward 1.0 to recover the old max-based behaviour.
 */
// ==========================================================================
static void DrawVectorFieldPanel(TProfile2D* hX, TProfile2D* hY, TProfile2D* hZ,
                                  const char* title,
                                  double minEntries    = 50.,
                                  int    arrowBlockSize = 4,
                                  double scalePercentile = 0.95)
{
    if (!hX || !hY || !hZ) return;

    // ---- COLZ background ----
    TProfile2D* hZd = static_cast<TProfile2D*>(
        hZ->Clone(Form("hVFZ_tmp_%d", gCloneIdx++)));
    hZd->SetDirectory(nullptr);
    hZd->SetStats(0);

    // Scale the clone by 100 to convert the colormap to a percentage
    hZd->Scale(100.0);

    if (title) hZd->SetTitle(title);

    double zmax = 0.;
    for (int ix = 1; ix <= hZd->GetNbinsX(); ++ix)
        for (int iy = 1; iy <= hZd->GetNbinsY(); ++iy) {
            int gb = hZd->GetBin(ix, iy);
            if (hZd->GetBinEntries(gb) < minEntries) continue;
            double v = std::fabs(hZd->GetBinContent(ix, iy));
            if (v > zmax) zmax = v;
        }
    zmax = (zmax < 1.e-10) ? 1.e-4 : zmax * 1.1;
    hZd->SetMinimum(-zmax);
    hZd->SetMaximum( zmax);
    hZd->Draw("COLZ");
    gPad->Update();

    // ---- Tile grid parameters ----
    int    bs      = arrowBlockSize;
    int    nBinsX  = hX->GetNbinsX();
    int    nBinsY  = hX->GetNbinsY();
    double binW    = hX->GetXaxis()->GetBinWidth(1);
    double tileW   = bs * binW;
    int    nTilesX = nBinsX / bs;
    int    nTilesY = nBinsY / bs;

    // ---- Pass 1: accumulate tile averages & propagate errors ----
    struct Tile { 
        double xc, yc;
        double bx, by, err_bx, err_by;
        double mag, err_mag;
        bool valid = false; 
    };
    std::vector<Tile> tiles(nTilesX * nTilesY);

    for (int itx = 0; itx < nTilesX; ++itx) {
        for (int ity = 0; ity < nTilesY; ++ity) {
            double sumBx = 0., sumBy = 0.;
            double sumErrBx2 = 0., sumErrBy2 = 0.;
            int    nUsed = 0;
            
            for (int dix = 0; dix < bs; ++dix) {
                for (int diy = 0; diy < bs; ++diy) {
                    int ix = itx * bs + dix + 1;
                    int iy = ity * bs + diy + 1;
                    int gb = hX->GetBin(ix, iy);
                    if (hX->GetBinEntries(gb) < minEntries) continue;
                    
                    sumBx += hX->GetBinContent(ix, iy);
                    sumBy += hY->GetBinContent(ix, iy);
                    
                    double errX = hX->GetBinError(ix, iy);
                    double errY = hY->GetBinError(ix, iy);
                    sumErrBx2 += errX * errX;
                    sumErrBy2 += errY * errY;
                    
                    nUsed++;
                }
            }

            int    firstBinX = itx * bs + 1, lastBinX = firstBinX + bs - 1;
            int    firstBinY = ity * bs + 1, lastBinY = firstBinY + bs - 1;
            double xc = 0.5 * (hX->GetXaxis()->GetBinLowEdge(firstBinX) +
                                hX->GetXaxis()->GetBinUpEdge (lastBinX));
            double yc = 0.5 * (hX->GetYaxis()->GetBinLowEdge(firstBinY) +
                                hX->GetYaxis()->GetBinUpEdge (lastBinY));

            Tile& t = tiles[itx * nTilesY + ity];
            t.xc = xc;  t.yc = yc;
            if (nUsed > 0) {
                t.bx     = sumBx / nUsed;
                t.by     = sumBy / nUsed;
                t.err_bx = std::sqrt(sumErrBx2) / nUsed;
                t.err_by = std::sqrt(sumErrBy2) / nUsed;
                
                t.mag = std::sqrt(t.bx * t.bx + t.by * t.by);
                
                // Magnitude error propagation
                if (t.mag > 1.e-12) 
                    t.err_mag = std::sqrt(t.bx * t.bx * t.err_bx * t.err_bx + t.by * t.by * t.err_by * t.err_by) / t.mag;
                else
                    t.err_mag = 0.;
                
                t.valid = true;
            }
        }
    }

    // ---- Collect valid tiles & sort by magnitude ----
    std::vector<Tile*> validTiles;
    validTiles.reserve(tiles.size());
    for (Tile& t : tiles) {
        if (t.valid) validTiles.push_back(&t);
    }
    if (validTiles.empty()) return;

    std::sort(validTiles.begin(), validTiles.end(),
              [](const Tile* a, const Tile* b) { return a->mag < b->mag; });

    // ---- Percentile index ----
    int pIdx = static_cast<int>(scalePercentile * static_cast<double>(validTiles.size() - 1)); // Notice the -1 subtraction to not bias upwards.
                                                                                               // For instance, if we have 100 elements, the 95th percentile
                                                                                               // is actually on index 94, as indices go from 0 to 99
    pIdx = std::max(0, std::min(pIdx, static_cast<int>(validTiles.size()) - 1));

    // ---- Percentile reference ----
    double scaleRef    = validTiles[pIdx]->mag;
    double scaleRefErr = validTiles[pIdx]->err_mag;

    // Fallback: if the chosen percentile lands on zero (e.g. many empty tiles),
    // walk up to the first non-zero magnitude so we always draw something
    if (scaleRef < 1.e-12) {
        for (Tile* t : validTiles) {
            if (t->mag > 1.e-12) {
                scaleRef    = t->mag;
                scaleRefErr = t->err_mag; // Ensures the error tracks the new fallback magnitude
                break;
            }
        }
    }
    if (scaleRef < 1.e-12) return;

    // full-reference arrow = 0.6 tile widths
    double scale = 0.6 * tileW / scaleRef;

    // ---- Pass 2: draw arrows, capping length at the scale reference ----
    for (const Tile& t : tiles) {
        if (!t.valid) continue;
        if (t.mag < 0.05 * scaleRef) continue; // suppress near-zero noise

        // Cap outlier arrows at scaleRef length; direction is always correct
        double drawLen = std::min(t.mag, scaleRef) * scale;
        double x2 = t.xc + (t.bx / t.mag) * drawLen;
        double y2 = t.yc + (t.by / t.mag) * drawLen;

        TArrow* arr = new TArrow(t.xc, t.yc, x2, y2, 0.012, ">");
        arr->SetLineColor(kBlack);
        arr->SetFillColor(kBlack);
        arr->SetLineWidth(2);
        arr->Draw();
    }

    // Reporting the maximum scale of magnitudes -- this must be on top of the arrows, for readability:
    // ----- Compute plot-coordinate placement -----
    double xMin = hZd->GetXaxis()->GetXmin();
    double xMax = hZd->GetXaxis()->GetXmax();
    double yMin = hZd->GetYaxis()->GetXmin();
    double yMax = hZd->GetYaxis()->GetXmax();
    double x1 = xMin + 0.05 * (xMax - xMin);
    double x2 = xMin + 0.65 * (xMax - xMin); // Increased from 0.48 to 0.58 to accommodate the error bar
    double y1 = yMax - 0.16 * (yMax - yMin);
    double y2 = yMax - 0.05 * (yMax - yMin);

    // ----- Background box -----
    TPaveText* pave = new TPaveText(x1, y1, x2, y2, "arc");
    pave->SetCornerRadius(0.15);
    pave->SetFillColor(kWhite);
    pave->SetFillStyle(1001);
    pave->SetBorderSize(0);
    pave->SetMargin(0.02);
    pave->SetTextAlign(12); // left-center in the pave
    pave->SetTextFont(63);
    pave->SetTextSize(18);
    
    double percentile = scalePercentile * 100.;
    pave->AddText(Form("|#LTp*_{T}#GT|_{%.0fpct} = (%.2f #pm %.2f)%%", percentile, scaleRef * 100., scaleRefErr * 100.));
    pave->Draw();
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
 *
 * Single panel figure displaying the profile of the jet ring observable 
 * proxy (<R_proxy_jet>) as a function of the Lambda transverse momentum (pT).
 * Overlays all four cut scenarios (No Cuts, pT Cut Only, DCA Cut Only, Both Cuts)
 * using the inclusive ("All") pseudorapidity selection.
 *
 * @param famDir    Family directory in the input file.
 * @param famOut    Output sub-directory.
 * @param famLabel  Short label used in object names.
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
 * @brief Fig 12 -- Integrated <R_proxy_jet> bar chart per scenario.
 *
 * 3 panels:
 * Left:   EtaPos (#eta_{#Lambda} > 0)
 * Middle: EtaNeg (#eta_{#Lambda} < 0)
 * Right:  All Eta
 *
 * Each panel displays a 1D chart of the integrated jet ring observable proxy 
 * (<R_proxy_jet>) extracted from the single-bin TProfiles for the four 
 * different cut scenarios.
 *
 * @param famDir    Family directory in the input file.
 * @param famOut    Output sub-directory.
 * @param famLabel  Short label used in object names.
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
 * @brief Fig 12 (event-mean variant) -- integrated <R_proxy_jet> bar chart
 *        using the event-mean uncertainty estimator, one panel per eta selection.
 *
 * @details
 * Mirrors MakeFig12_IntegratedRingJet exactly in layout, but replaces the
 * TProfile SEM with the event-mean-based uncertainty that correctly accounts
 * for intra-event correlations arising from Lambdas sharing the same jet
 * direction.
 *
 * Central value: hEventMeanRingProxyJet->GetMean()
 *   -- mean of the per-event averages R_e = (1/N_e) sum_i R_i
 * Uncertainty:   hEventMeanRingProxyJet->GetMeanError() = StdDev(R_e) / sqrt(N_events)
 *   -- properly treats jet groups as the independent sampling unit
 *
 * @param famDir    Family directory in the input file.
 * @param famOut    Output sub-directory in the plots ROOT file.
 * @param famLabel  Short label used in canvas and object names.
 */
// ==========================================================================
static void MakeFig12_IntegratedRingJetEvt(TDirectory* famDir, TDirectory* famOut, const char* famLabel)
{
    const char* etaSels[3]   = {"EtaPos", "EtaNeg", "All"};
    const char* etaLabels[3] = {"#eta_{#Lambda} > 0", "#eta_{#Lambda} < 0", "All #eta"};

    // Global y range from event-mean errors across all panels
    double globalMax = 0.;
    for (int ie = 0; ie < 3; ++ie) {
        for (int is = 0; is < 4; ++is) {
            TDirectory* dir = GetScenarioDir(famDir, kScenNames[is], etaSels[ie]);
            if (!dir) continue;
            TH1D* h = static_cast<TH1D*>(SafeGet(dir, "hEventMeanRingProxyJet"));
            if (!h || h->GetEntries() == 0) continue;
            double v = std::fabs(h->GetMean()) + h->GetMeanError();
            if (v > globalMax) globalMax = v;
        }
    }
    globalMax = (globalMax < 1.e-8) ? 0.01 : globalMax * 1.5;

    TCanvas* c = new TCanvas(Form("c_%s_intRingJetEvt", famLabel), "", 1050, 500);
    c->Divide(3, 1, 0.004, 0.002);

    for (int ie = 0; ie < 3; ++ie) {
        c->cd(ie + 1);
        gPad->SetLeftMargin(0.18);  gPad->SetBottomMargin(0.16);

        TH1D* hBar = new TH1D(Form("hBarJetEvt_%s_%s", famLabel, etaSels[ie]),
                               Form("<R_{proxy}^{jet}> (event-mean) -- %s;Scenario;<R_{proxy}^{jet}>", etaLabels[ie]),
                               4, 0., 4.);
        hBar->SetDirectory(nullptr);
        hBar->GetYaxis()->SetRangeUser(-globalMax, globalMax);
        hBar->SetStats(0);
        for (int is = 0; is < 4; ++is) hBar->GetXaxis()->SetBinLabel(is + 1, kScenLabels[is]);
        hBar->GetXaxis()->SetLabelSize(0.055);
        hBar->Draw("AXIS");

        TLine* zl = new TLine(0., 0., 4., 0.);
        zl->SetLineColor(kGray + 2);  zl->SetLineStyle(2);  zl->SetLineWidth(2);  zl->Draw("SAME");

        for (int is = 0; is < 4; ++is) {
            TDirectory* dir = GetScenarioDir(famDir, kScenNames[is], etaSels[ie]);
            if (!dir) continue;
            TH1D* h = static_cast<TH1D*>(SafeGet(dir, "hEventMeanRingProxyJet"));
            if (!h || h->GetEntries() == 0) continue;

            TH1D* hPt = new TH1D(Form("hPtJetEvt12_%s_%d_%d", famLabel, ie, is), "", 4, 0., 4.);
            hPt->SetDirectory(nullptr);
            hPt->SetBinContent(is + 1, h->GetMean());
            hPt->SetBinError  (is + 1, h->GetMeanError());
            SetHistStyle(hPt, kScenColors[is], kScenMarkers[is]);
            hPt->SetLineWidth(3);
            hPt->Draw("E1 SAME");
        }
        AddLabel(0.5, 0.96, etaLabels[ie], 0.045, 22);
    }

    c->cd(0);
    AddLabel(0.5, 0.995, Form("Integrated <R_{proxy}^{jet}> (event-mean estimator) -- %s", famLabel), 0.036, 22);
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
 * @brief Fig 14 -- three-way overlay comparison of TProfile SEM, event-mean
 *        uncertainty, and data-chunking uncertainty for the integrated
 *        <R_proxy_jet>, one panel per eta selection.
 *
 * @details
 * For each (scenario, eta) combination, three points are drawn at slightly
 * offset x positions so all estimators are visible simultaneously:
 *
 *   Filled marker  (x - 0.18): naive TProfile SEM; assumes all N_Lambda entries
 *                               are independent -- underestimates if intra-event
 *                               correlations from the shared jet axis are present.
 *
 *   Open marker    (x + 0.00): event-mean estimator; sigma(R_e)/sqrt(N_events),
 *                               where R_e is the per-jet-group mean.  Correctly
 *                               treats the jet group as the independent unit, but
 *                               uses a weighted TH1D internally, which can still
 *                               carry a ROOT-level weighting artefact.
 *
 *   Star marker    (x + 0.18): data-chunking estimator; stddev(mu_k)/sqrt(K_filled),
 *                               where mu_k = raw_sum_k / raw_count_k for each of
 *                               the kChunks sequential chunks of jet-group events.
 *                               Entirely unweighted; immune to ROOT's effective-
 *                               degrees-of-freedom trap.
 *
 * Hierarchy of trust (highest to lowest): chunking > event-mean > TProfile SEM.
 * If all three agree, the intra-event correlation is negligible at the current
 * jet-group size.  If the chunking bar is taller than the TProfile bar, the SEM
 * was underestimating.
 *
 * @param famDir    Family directory in the input file.
 * @param famOut    Output sub-directory in the plots ROOT file.
 * @param famLabel  Short label used in canvas and object names.
 */
// ==========================================================================
static void MakeFig14_RingJetErrComparison(TDirectory* famDir, TDirectory* famOut, const char* famLabel)
{
    const char* etaSels[3]   = {"EtaPos", "EtaNeg", "All"};
    const char* etaLabels[3] = {"#eta_{#Lambda} > 0", "#eta_{#Lambda} < 0", "All #eta"};

    // Open marker counterparts for kScenMarkers {20,21,22,23}
    static const int kOpenMarkers[4] = {24, 25, 26, 32};

    // Global y range covering all three estimators
    double globalMax = 0.;
    for (int ie = 0; ie < 3; ++ie) {
        for (int is = 0; is < 4; ++is) {
            TDirectory* dir = GetScenarioDir(famDir, kScenNames[is], etaSels[ie]);
            if (!dir) continue;
            TProfile* p  = static_cast<TProfile*>(SafeGet(dir, "pRingProxyJet"));
            TH1D*     he = static_cast<TH1D*>(SafeGet(dir, "hEventMeanRingProxyJet"));
            TH1D*     hc = static_cast<TH1D*>(SafeGet(dir, "hChunkMeansRingProxyJet"));
            if (p) {
                double v = std::fabs(p->GetBinContent(1)) + p->GetBinError(1);
                if (v > globalMax) globalMax = v;
            }
            if (he && he->GetEntries() > 0) {
                double v = std::fabs(he->GetMean()) + he->GetMeanError();
                if (v > globalMax) globalMax = v;
            }
            if (hc && hc->GetEntries() > 0) {
                double v = std::fabs(hc->GetMean()) + hc->GetMeanError();
                if (v > globalMax) globalMax = v;
            }
        }
    }
    globalMax = (globalMax < 1.e-8) ? 0.01 : globalMax * 1.6;

    TCanvas* c = new TCanvas(Form("c_%s_ringJetErrComp", famLabel), "", 1050, 500);
    c->Divide(3, 1, 0.004, 0.002);

    for (int ie = 0; ie < 3; ++ie) {
        c->cd(ie + 1);
        gPad->SetLeftMargin(0.18);  gPad->SetBottomMargin(0.16);

        // Axis frame with scenario bin labels
        TH1D* hBar = new TH1D(Form("hBarErrComp_%s_%s", famLabel, etaSels[ie]),
                               Form("<R_{proxy}^{jet}> -- %s;Scenario;<R_{proxy}^{jet}>", etaLabels[ie]),
                               4, 0., 4.);
        hBar->SetDirectory(nullptr);
        hBar->GetYaxis()->SetRangeUser(-globalMax, globalMax);
        hBar->SetStats(0);
        for (int is = 0; is < 4; ++is) hBar->GetXaxis()->SetBinLabel(is + 1, kScenLabels[is]);
        hBar->GetXaxis()->SetLabelSize(0.055);
        hBar->Draw("AXIS");

        TLine* zl = new TLine(0., 0., 4., 0.);
        zl->SetLineColor(kGray + 2);  zl->SetLineStyle(2);  zl->SetLineWidth(2);  zl->Draw("SAME");

        for (int is = 0; is < 4; ++is) {
            TDirectory* dir = GetScenarioDir(famDir, kScenNames[is], etaSels[ie]);
            if (!dir) continue;
            TProfile* p  = static_cast<TProfile*>(SafeGet(dir, "pRingProxyJet"));
            TH1D*     he = static_cast<TH1D*>(SafeGet(dir, "hEventMeanRingProxyJet"));
            TH1D*     hc = static_cast<TH1D*>(SafeGet(dir, "hChunkMeansRingProxyJet"));

            // Bin centre for scenario 'is' in a 4-bin axis spanning [0,4]
            double xCtr = is + 0.5;

            // -- TProfile SEM: filled marker, shifted left --
            if (p) {
                TGraphErrors* grS = new TGraphErrors(1);
                grS->SetPoint    (0, xCtr - 0.18, p->GetBinContent(1));
                grS->SetPointError(0, 0.,          p->GetBinError(1));
                grS->SetMarkerStyle(kScenMarkers[is]);
                grS->SetMarkerColor(kScenColors[is]);
                grS->SetLineColor  (kScenColors[is]);
                grS->SetLineWidth(2);
                grS->SetMarkerSize(1.1);
                grS->Draw("P SAME");
            }

            // -- Event-mean sigma: open marker, centred --
            if (he && he->GetEntries() > 0) {
                TGraphErrors* grE = new TGraphErrors(1);
                grE->SetPoint    (0, xCtr + 0.00, he->GetMean());
                grE->SetPointError(0, 0.,          he->GetMeanError());
                grE->SetMarkerStyle(kOpenMarkers[is]);
                grE->SetMarkerColor(kScenColors[is]);
                grE->SetLineColor  (kScenColors[is]);
                grE->SetLineWidth(2);
                grE->SetMarkerSize(1.1);
                grE->Draw("P SAME");
            }

            // -- Chunking estimator: full-star marker, shifted right --
            if (hc && hc->GetEntries() > 0) {
                TGraphErrors* grC = new TGraphErrors(1);
                grC->SetPoint    (0, xCtr + 0.18, hc->GetMean());
                grC->SetPointError(0, 0.,          hc->GetMeanError());
                grC->SetMarkerStyle(29); // full 5-pointed star -- visually distinct from filled/open sets
                grC->SetMarkerColor(kScenColors[is]);
                grC->SetLineColor  (kScenColors[is]);
                grC->SetLineWidth(2);
                grC->SetMarkerSize(1.4); // stars look small at 1.1; bump up slightly
                grC->Draw("P SAME");
            }
        }
        AddLabel(0.5, 0.96, etaLabels[ie], 0.045, 22);

        // Legend on the rightmost (All) panel only
        if (ie == 2) {
            TLegend* leg = MakeLegend(0.18, 0.68, 0.95, 0.88);
            for (int is = 0; is < 4; ++is) {
                TLine* dummy = new TLine();
                dummy->SetLineColor(kScenColors[is]);  dummy->SetLineWidth(2);
                leg->AddEntry(dummy, kScenLabels[is], "l");
            }
            TMarker* mFull = new TMarker(0., 0., 20);
            mFull->SetMarkerColor(kBlack);  mFull->SetMarkerSize(1.1);
            leg->AddEntry(mFull, "TProfile SEM (filled, left)", "p");
            TMarker* mOpen = new TMarker(0., 0., 24);
            mOpen->SetMarkerColor(kBlack);  mOpen->SetMarkerSize(1.1);
            leg->AddEntry(mOpen, "Event-mean #sigma (open, centre)", "p");
            TMarker* mStar = new TMarker(0., 0., 29);
            mStar->SetMarkerColor(kBlack);  mStar->SetMarkerSize(1.4);
            leg->AddEntry(mStar, "Chunking #sigma (star, right)", "p");
            leg->Draw("SAME");
        }
    }

    c->cd(0);
    AddLabel(0.5, 0.995, Form("SEM vs event-mean vs chunking uncertainty -- %s", famLabel), 0.036, 22);
    WriteCanvas(c, famOut);
    delete c;
}

// ==========================================================================
/**
 * @brief Fig 15 -- distribution of per-event R_proxy_jet means, all scenarios.
 *
 * @details
 * Analog of Fig13's left panel, but drawing hEventMeanRingProxyJet (the
 * distribution of R_e = mean ringProxyJet within each jet group) rather than
 * the per-Lambda h1d_ringProxyJet.
 *
 * This distribution is expected to be:
 *   - Much narrower than h1d_ringProxyJet by a factor ~1/sqrt(N_per_event),
 *   - Approximately Gaussian (central limit theorem acting on the group means),
 *   - Centred near zero for all scenarios (any shift would indicate a bias).
 *
 * Its standard deviation divided by sqrt(N_events) is exactly the event-mean
 * uncertainty plotted in MakeFig12_IntegratedRingJetEvt and MakeFig14.
 *
 * @param famDir    Family directory in the input file.
 * @param famOut    Output sub-directory in the plots ROOT file.
 * @param famLabel  Short label used in canvas and object names.
 */
// ==========================================================================
static void MakeFig15_EventMeanDist(TDirectory* famDir, TDirectory* famOut, const char* famLabel)
{
    TH1D* hNC = nullptr;  TH1D* hPT = nullptr;
    TH1D* hDC = nullptr;  TH1D* hBC = nullptr;
    {
        TDirectory* dNC = GetScenarioDir(famDir, "NoCuts",     "All");
        TDirectory* dPT = GetScenarioDir(famDir, "pTCutOnly",  "All");
        TDirectory* dDC = GetScenarioDir(famDir, "DCACutOnly", "All");
        TDirectory* dBC = GetScenarioDir(famDir, "BothCuts",   "All");
        if (dNC) hNC = static_cast<TH1D*>(SafeGet(dNC, "hEventMeanRingProxyJet"));
        if (dPT) hPT = static_cast<TH1D*>(SafeGet(dPT, "hEventMeanRingProxyJet"));
        if (dDC) hDC = static_cast<TH1D*>(SafeGet(dDC, "hEventMeanRingProxyJet"));
        if (dBC) hBC = static_cast<TH1D*>(SafeGet(dBC, "hEventMeanRingProxyJet"));
    }
    if (!hNC) return;

    TH1D* hNCc = SafeClone(hNC);
    TH1D* hPTc = SafeClone(hPT);
    TH1D* hDCc = SafeClone(hDC);
    TH1D* hBCc = SafeClone(hBC);

    SetHistStyle(hNCc, kColNoCuts, kMarkerNC);  SetHistStyle(hPTc, kColPtCut,  kMarkerPT);
    SetHistStyle(hDCc, kColDcaCut, kMarkerDC);  SetHistStyle(hBCc, kColBoth,   kMarkerBC);

    double ymax = 0.;
    for (TH1D* h : {hNCc, hPTc, hDCc, hBCc}) {
        if (h && h->GetMaximum() > ymax) ymax = h->GetMaximum();
    }
    ymax *= 1.3;

    TCanvas* c = new TCanvas(Form("c_%s_evtMeanDist", famLabel), "", 750, 600);
    c->SetLeftMargin(0.14);  c->SetBottomMargin(0.13);

    if (hNCc) {
        hNCc->SetTitle(Form("Per-event mean R_{proxy}^{jet} distribution -- %s;"
                            "#bar{R}_{proxy}^{jet} per jet group;Events", famLabel));
        hNCc->GetYaxis()->SetRangeUser(0., ymax);
        hNCc->Draw("HIST");
    }
    if (hPTc) hPTc->Draw("HIST SAME");
    if (hDCc) hDCc->Draw("HIST SAME");
    if (hBCc) hBCc->Draw("HIST SAME");

    TLine* zl = new TLine(0., 0., 0., ymax);
    zl->SetLineColor(kGray + 2);  zl->SetLineStyle(3);  zl->Draw("SAME");

    TLegend* leg = MakeLegend(0.15, 0.70, 0.60, 0.88);
    if (hNCc) leg->AddEntry(hNCc, kScenLabels[0], "l");
    if (hPTc) leg->AddEntry(hPTc, kScenLabels[1], "l");
    if (hDCc) leg->AddEntry(hDCc, kScenLabels[2], "l");
    if (hBCc) leg->AddEntry(hBCc, kScenLabels[3], "l");
    leg->Draw("SAME");

    WriteCanvas(c, famOut);
    delete c;
}

// ==========================================================================
/**
 * @brief Fig 16 -- main proton rest-frame polarisation vector field.
 *
 * @details
 * A 2 x 3 canvas (cut x eta-selection):
 *   Row 1: NoCuts    x { All eta | eta>0 | eta<0 }
 *   Row 2: BothCuts  x { All eta | eta>0 | eta<0 }
 *
 * Each panel shows the average proton direction in the Lambda rest frame as
 * a 2D vector field in the Lambda transverse momentum plane (px, py):
 *   - Background colour:  <p*_z>  (out-of-plane component, diverging palette,
 *                          centred at zero; encodes any azimuthal z-polarisation
 *                          pattern induced by the cuts).
 *   - Arrow direction:   (<p*_x>, <p*_y>) per (px_lam, py_lam) cell.
 *   - Arrow length:       normalised to the global maximum transverse magnitude
 *                          within the panel so relative strengths are preserved.
 *
 * The expected result in the NoCuts panel is a uniform near-zero field; after
 * DCA cuts a ring-like pattern (arrows tangent to circles of constant pT)
 * should emerge, driven by the magnetic-field-dependent acceptance.
 *
 * The supplemental figure MakeFig16s_PstarVectorFieldSupp shows the two
 * intermediate cut scenarios (pTCutOnly, DCACutOnly).
 *
 * @param famDir    Family directory in the input file.
 * @param famOut    Output sub-directory in the plots ROOT file.
 * @param famLabel  Short label used in canvas and object names.
 */
// ==========================================================================
static void MakeFig16_PstarVectorField(TDirectory* famDir, TDirectory* famOut,
                                        const char* famLabel)
{
    const char* etaSels[3]    = {"All",       "EtaPos",            "EtaNeg"};
    const char* etaLabels[3]  = {"All #eta",  "#eta_{#Lambda}>0",  "#eta_{#Lambda}<0"};
    const char* cutRows[2]    = {"NoCuts",    "BothCuts"};
    const char* cutLabels[2]  = {"No Cuts",   "Both Cuts"};

    // gStyle->SetPalette(kTemperatureMap); // diverging blue-white-red, zero = white

    TCanvas* c = new TCanvas(Form("c_%s_pstarVF", famLabel), "", 1500, 950);
    c->Divide(3, 2, 0.002, 0.002);

    int panel = 0;
    for (int irow = 0; irow < 2; ++irow) {
        for (int icol = 0; icol < 3; ++icol) {
            c->cd(++panel);
            gPad->SetTopMargin(0.05);
            if (icol == 0)
                gPad->SetLeftMargin(0.08);
            else
                gPad->SetLeftMargin(0.09);
            gPad->SetRightMargin(0.16);

            if (irow == 0)
                gPad->SetTopMargin(0.12);
            else
                gPad->SetTopMargin(0.05);

            if (irow == 1)
                gPad->SetBottomMargin(0.10);
            else
                gPad->SetBottomMargin(0.08);

            TDirectory* dir = GetScenarioDir(famDir, cutRows[irow], etaSels[icol]);
            if (!dir) continue;

            TProfile2D* hX = static_cast<TProfile2D*>(SafeGet(dir, "pPstarX_vsPxPy"));
            TProfile2D* hY = static_cast<TProfile2D*>(SafeGet(dir, "pPstarY_vsPxPy"));
            TProfile2D* hZ = static_cast<TProfile2D*>(SafeGet(dir, "pPstarZ_vsPxPy"));
            if (!hX || !hY || !hZ) continue;

            // Giving a better offset for the Z axis label -- This is done outside of the function to avoid gPad updates:
            if (icol == 2)
                hZ->GetZaxis()->SetTitleOffset(1.5);
            else
                hZ->GetZaxis()->SetTitleOffset(1.6);
            DrawVectorFieldPanel(hX, hY, hZ, Form(" ;p_{x}^{#Lambda} [GeV/c];p_{y}^{#Lambda} [GeV/c];<p*_{z}> [%]")); // No title works best here, so don't use Form("<p*_{z}> + #vec{<p*_{T}>};"
                                                                                                                      // The Z axis will be converted to a percentage only inside the function call,
                                                                                                                      // to preserve the original object and only alter the function's inner clone.
            if (irow == 0)
                AddLabel(0.50, 0.905, Form("%s  --  %s", cutLabels[irow], etaLabels[icol]), 0.040, 22);
            else
                AddLabel(0.50, 0.975, Form("%s  --  %s", cutLabels[irow], etaLabels[icol]), 0.040, 22);
        }
    }

    c->cd(0);
    AddLabel(0.5, 0.98, Form("Spurious #LTp*#GT polarisation vector field  --  %s", famLabel), 0.033, 22);
    WriteCanvas(c, famOut);
    delete c;
}

// ==========================================================================
/**
 * @brief Fig 16s -- supplemental proton rest-frame polarisation vector field.
 *
 * @details
 * Identical layout to MakeFig16_PstarVectorField but showing the two
 * intermediate cut scenarios:
 *   Row 1: pTCutOnly  x { All eta | eta>0 | eta<0 }
 *   Row 2: DCACutOnly x { All eta | eta>0 | eta<0 }
 *
 * Comparing with the main figure allows the pT-cut and DCA-cut contributions
 * to the fake ring pattern to be separated.  The DCA-only row is expected to
 * show the clearest ring-like structure, since it is the DCA cut that couples
 * directly to the helix geometry and thus to the Lambda azimuthal angle.
 *
 * @param famDir    Family directory in the input file.
 * @param famOut    Output sub-directory in the plots ROOT file.
 * @param famLabel  Short label used in canvas and object names.
 */
// ==========================================================================
static void MakeFig16s_PstarVectorFieldSupp(TDirectory* famDir, TDirectory* famOut,
                                             const char* famLabel)
{
    const char* etaSels[3]    = {"All",       "EtaPos",            "EtaNeg"};
    const char* etaLabels[3]  = {"All #eta",  "#eta_{#Lambda}>0",  "#eta_{#Lambda}<0"};
    const char* cutRows[2]    = {"pTCutOnly",     "DCACutOnly"};
    const char* cutLabels[2]  = {"p_{T} Cut Only", "DCA Cut Only"};

    // gStyle->SetPalette(kTemperatureMap);

    TCanvas* c = new TCanvas(Form("c_%s_pstarVFsupp", famLabel), "", 1500, 950);
    c->Divide(3, 2, 0.002, 0.002);

    int panel = 0;
    for (int irow = 0; irow < 2; ++irow) {
        for (int icol = 0; icol < 3; ++icol) {
            c->cd(++panel);
            gPad->SetTopMargin(0.05);
            if (icol == 0)
                gPad->SetLeftMargin(0.08);
            else
                gPad->SetLeftMargin(0.09);
            gPad->SetRightMargin(0.16);

            if (irow == 0)
                gPad->SetTopMargin(0.12);
            else
                gPad->SetTopMargin(0.05);

            if (irow == 1)
                gPad->SetBottomMargin(0.10);
            else
                gPad->SetBottomMargin(0.08);

            TDirectory* dir = GetScenarioDir(famDir, cutRows[irow], etaSels[icol]);
            if (!dir) continue;

            TProfile2D* hX = static_cast<TProfile2D*>(SafeGet(dir, "pPstarX_vsPxPy"));
            TProfile2D* hY = static_cast<TProfile2D*>(SafeGet(dir, "pPstarY_vsPxPy"));
            TProfile2D* hZ = static_cast<TProfile2D*>(SafeGet(dir, "pPstarZ_vsPxPy"));
            if (!hX || !hY || !hZ) continue;

            // Giving a better offset for the Z axis label -- This is done outside of the function to avoid gPad updates:
            if (icol == 2)
                hZ->GetZaxis()->SetTitleOffset(1.5);
            else
                hZ->GetZaxis()->SetTitleOffset(1.6);
            DrawVectorFieldPanel(hX, hY, hZ, Form(" ;p_{x}^{#Lambda} [GeV/c];p_{y}^{#Lambda} [GeV/c];<p*_{z}>")); // No title works best here, so don't use Form("<p*_{z}> + #vec{<p*_{T}>};"
            if (irow == 0)
                AddLabel(0.50, 0.905, Form("%s  --  %s", cutLabels[irow], etaLabels[icol]), 0.040, 22);
            else
                AddLabel(0.50, 0.975, Form("%s  --  %s", cutLabels[irow], etaLabels[icol]), 0.040, 22);
        }
    }

    c->cd(0);
    AddLabel(0.5, 0.98, Form("Spurious #LTp*#GT vector field (intermediate cuts)  --  %s", famLabel), 0.033, 22);
    WriteCanvas(c, famOut);
    delete c;
}

// ==========================================================================
/**
 * @brief Fig 17 -- <p*> components vs Lambda azimuthal angle phi_lam.
 *
 * @details
 * Three canvases (one per eta selection: All, EtaPos, EtaNeg).  Each canvas
 * has three panels showing <p*_x>, <p*_y>, <p*_z> vs phi_lam, with all four
 * cut scenarios overlaid in every panel.
 *
 * This is the compact 1D projection of the vector field: a ring-like
 * polarisation bias from DCA cuts will appear as a sinusoidal modulation
 * in <p*_x> and <p*_y> with period 2pi and a phase determined by the
 * magnetic field orientation.  <p*_z> should remain flat if there is no
 * out-of-plane component to the fake signal.
 *
 * The y-axis range is set globally across all three component panels and
 * both positive and negative extremes so the relative amplitudes are
 * directly comparable between panels.
 *
 * @param famDir    Family directory in the input file.
 * @param famOut    Output sub-directory in the plots ROOT file.
 * @param famLabel  Short label used in canvas and object names.
 */
// ==========================================================================
static void MakeFig17_PstarVsPhiLam(TDirectory* famDir, TDirectory* famOut,
                                     const char* famLabel)
{
    const char* etaSels[3]   = {"All",      "EtaPos",            "EtaNeg"};
    const char* etaLabels[3] = {"All #eta", "#eta_{#Lambda}>0",  "#eta_{#Lambda}<0"};

    // Profile names and axis labels for the three components
    const char* profNames[3]  = {"pPstarX_vsPhiLam",
                                  "pPstarY_vsPhiLam",
                                  "pPstarZ_vsPhiLam"};
    const char* compLabels[3] = {"<p*_{x}>", "<p*_{y}>", "<p*_{z}>"};

    for (int ie = 0; ie < 3; ++ie) {

        // Collect all 12 profiles (4 scenarios x 3 components) and find global y range
        TProfile* profs[4][3]; // [scenario][component]
        double ymax = 0.;
        for (int is = 0; is < 4; ++is) {
            TDirectory* dir = GetScenarioDir(famDir, kScenNames[is], etaSels[ie]);
            for (int ic = 0; ic < 3; ++ic) {
                profs[is][ic] = (dir)
                    ? static_cast<TProfile*>(SafeGet(dir, profNames[ic]))
                    : nullptr;
                TProfile* p = profs[is][ic];
                if (!p) continue;
                for (int ib = 1; ib <= p->GetNbinsX(); ++ib) {
                    double v = std::fabs(p->GetBinContent(ib)) + p->GetBinError(ib);
                    if (v > ymax) ymax = v;
                }
            }
        }
        ymax = (ymax < 1.e-10) ? 0.01 : ymax * 1.45;

        TCanvas* c = new TCanvas(
            Form("c_%s_pstarPhi_%s", famLabel, etaSels[ie]), "", 1400, 500);
        c->Divide(3, 1, 0.004, 0.002);

        for (int ic = 0; ic < 3; ++ic) {
            c->cd(ic + 1);
            gPad->SetLeftMargin(0.15);
            gPad->SetBottomMargin(0.15);

            // Invisible axis frame so the y range and title are set before any Draw
            TH1D* hFrame = new TH1D(
                Form("hVFframe_%s_%s_%d", famLabel, etaSels[ie], ic),
                Form("%s vs #phi_{#Lambda}  --  %s;"
                     "#phi_{#Lambda} [rad];%s",
                     compLabels[ic], etaLabels[ie], compLabels[ic]),
                1, -TMath::Pi(), TMath::Pi());
            hFrame->SetDirectory(nullptr);
            hFrame->SetStats(0);
            hFrame->GetYaxis()->SetRangeUser(-ymax, ymax);
            hFrame->Draw("AXIS");

            TLine* zl = new TLine(-TMath::Pi(), 0., TMath::Pi(), 0.);
            zl->SetLineColor(kGray + 2);
            zl->SetLineStyle(2);
            zl->Draw("SAME");

            for (int is = 0; is < 4; ++is) {
                TProfile* p = profs[is][ic];
                if (!p) continue;
                SetHistStyle(p, kScenColors[is], kScenMarkers[is]);
                p->Draw("EP SAME");
            }

            // Legend on the rightmost panel only
            if (ic == 2) {
                TLegend* leg = MakeLegend(0.16, 0.68, 0.62, 0.88);
                for (int is = 0; is < 4; ++is)
                    if (profs[is][ic])
                        leg->AddEntry(profs[is][ic], kScenLabels[is], "ep");
                leg->Draw("SAME");
            }
        }

        c->cd(0);
        AddLabel(0.5, 0.997,
            Form("#LTp*#GT components vs #phi_{#Lambda}  --  %s  --  %s",
                 famLabel, etaLabels[ie]),
            0.034, 22);
        WriteCanvas(c, famOut);
        delete c;
    }
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

        // ---- Output sub-directories (keeps each family folder navigable) ----
        TDirectory* outEmis = famOut->mkdir("EmissionAngles");    // figs 1-3
        TDirectory* outRingZ= famOut->mkdir("RingZHat");          // figs 4-6
        TDirectory* outDaug = famOut->mkdir("Daughters");         // figs 7-8
        // TDirectory* outEta  = famOut->mkdir("EtaAsymmetry");      // fig  9
        TDirectory* outJet  = famOut->mkdir("RingJet");           // figs 10-15
        TDirectory* outVF   = famOut->mkdir("PolarizationVectorField"); // figs 16-17

        printf("Producing figures for family: %s\n", famName);

        // Figs 1-3: one canvas per eta selection
        for (const char* eta : {"EtaPos", "EtaNeg", "All"}) {
            MakeFig1_2DMaps  (famDir, outEmis, famName, eta);
            MakeFig2_CosTheta(famDir, outEmis, famName, eta);
            MakeFig3_PhiStar (famDir, outEmis, famName, eta);
        }

        // Figs 4-6: z-hat ring proxy
        MakeFig4_RingProxy    (famDir, outRingZ, famName);
        MakeFig5_RingVsPt     (famDir, outRingZ, famName);
        MakeFig6_IntegratedRing(famDir, outRingZ, famName);

        // Figs 7-8: daughter kinematics
        MakeFig7_DaughterPt  (f, famDir, outDaug, famName);
        MakeFig8_DaughterDCA (f, famDir, outDaug, famName);

        // Fig 9: eta asymmetry
        MakeFig9_EtaDiff     (famDir, famOut, famName);

        // Figs 10-15: random jet-direction ring proxy and error estimators
        MakeFig10_RingProxyJetVsEta   (famDir, outJet, famName);
        MakeFig11_ringJetVsPt         (famDir, outJet, famName);
        MakeFig12_IntegratedRingJet   (famDir, outJet, famName);
        MakeFig12_IntegratedRingJetEvt(famDir, outJet, famName);
        MakeFig13_RingProxyJet        (famDir, outJet, famName);
        MakeFig14_RingJetErrComparison(famDir, outJet, famName);
        MakeFig15_EventMeanDist       (famDir, outJet, famName);

        // Figs 16-17: polarisation vector field
        MakeFig16_PstarVectorField    (famDir, outVF, famName);
        MakeFig16s_PstarVectorFieldSupp(famDir, outVF, famName);
        MakeFig17_PstarVsPhiLam       (famDir, outVF, famName);
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