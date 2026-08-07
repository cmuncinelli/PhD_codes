/*
* ============================================================
* auxiliaryPerConfigPlots.cxx
* ============================================================
*
* PURPOSE
* -------
* General-purpose, per-config post-processing macro for the lambdaJetPolarizationIonsDerived.cxx consumer output from O2Physics.
*
* Like makeCumulativeDCAdauProfile.cxx, this operates on a single ConsumerResults_<SUFFIX>.root file
* (one wagon and one consumer config) and is meant to be invoked once per (wagon, config) pair by
* run_all_wagons.sh.
* Unlike makeCumulativeDCAdauProfile.cxx, this file is NOT scoped to one specific derived quantity:
* it is meant to grow into the general home for "derivative" plots that are cheap/easy to
* build in a post-processing stage instead of inside the O2Physics consumer itself
* This code starts with the polarization vector-field plot.
*
* SECTIONS
* -----------------
*   1. Polarization vector-field plots ("pretty arrows" to compare with the Helicity Efficiency Toy Model):
*      <P*_x,y,z> vs the Lambda momentum plane, both (px,py) and (pz,px), analogous to the pPstarX_vsPxPy family
*      in helicityEfficiencyToyModel.cxx / plotHelicityEfficiency.cxx, but reading the real p2dP{x,y,z}Star_vsPxPy
*      / p2dP{x,y,z}Star_vsPzPx TProfile2Ds booked by the consumer.
*   2. Ring observable 2D scalar maps: <R> vs the same two momentum planes (p2dRingObservableVsPxPy / p2dRingObservableVsPzPx).
*      No arrow overlay is needed here (wouldn't make sense to do so, actually) -- the COLZ background alone is the physical result.
*
* Both sections are produced once per kinematic-cut folder (see kFolders below): Ring, RingKinematicCuts, JetKinematicCuts, JetAndLambdaKinematicCuts.
*
* ADDING A NEW SECTION: see the "ADD MORE POST-PROCESSING SECTIONS HERE" comment inside main() below.
*
* ON THE DUPLICATED DrawVectorFieldPanel():
* DrawVectorFieldPanel() below is copied close to verbatim from
* plotHelicityEfficiency.cxx (the toy-model plotter). This is a deliberate
* choice: the toy-model plotter and this per-config post-processor are two independent workflows,
* compiled and ran separately, and the function itself is small, self-contained, and should not
* change that often. Sharing it via a common header would couple two otherwise-unrelated build targets
* for little benefit. If it ever needs to diverge meaningfully or grows non-trivially, that is the
* trigger to revisit and actually factor it out.
*
* Usage mirrors makeCumulativeDCAdauProfile.cxx; compiled ahead-of-time
* by run_all_wagons.sh, not invoked as a ROOT interpreted macro.
* In case you do want to run it by hand, command would look like this:
*   ./auxiliaryPerConfigPlots.exe path/to/ConsumerResults_SUFFIX.root
* ============================================================
*/

#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TArrow.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

// ==========================================================================
// Small navigation / drawing helpers:
// (Signatures mirror the same-named helpers in plotHelicityEfficiency.cxx)
// ==========================================================================

static TObject* SafeGet(TDirectory* dir, const char* name)
{
    if (!dir) { printf("WARNING: SafeGet called with null directory for '%s'\n", name); return nullptr; }
    TObject* obj = dir->Get(name);
    if (!obj) printf("WARNING: '%s' not found in directory '%s'\n", name, dir->GetName());
    return obj;
}

static TDirectory* GetDir(TDirectory* parent, const char* name)
{
    if (!parent) { printf("WARNING: GetDir called with null parent for '%s'\n", name); return nullptr; }
    TDirectory* dir = static_cast<TDirectory*>(parent->Get(name));
    if (!dir) printf("WARNING: directory '%s' not found in '%s'\n", name, parent->GetName());
    return dir;
}

static void AddLabel(double x, double y, const char* text, double size = 0.04, int align = 22)
{
    TLatex* lat = new TLatex(x, y, text);
    lat->SetNDC();
    lat->SetTextSize(size);
    lat->SetTextAlign(align);
    lat->Draw();
}

static void WriteCanvas(TCanvas* c, TDirectory* outDir)
{
    if (!c || !outDir) return;
    outDir->cd();
    c->Write();
}

static int gCloneIdx = 0; // Unique-name counter for internal clones (mirrors plotHelicityEfficiency.cxx)
// ==========================================================================
/**
 * @brief Draws one vector-field panel: a COLZ background with block-averaged
 * transverse/longitudinal arrows overlaid.
 *
 * DUPLICATED FROM plotHelicityEfficiency.cxx
 * Behavior is intentionally identical to that copy, and the full parameter
 * documentation lives there.
 *
 * @param hX              TProfile2D for the arrow's horizontal component.
 * @param hY              TProfile2D for the arrow's vertical component.
 * @param hZ              TProfile2D for the COLZ background (out-of-plane component).
 * @param title           Histogram title string; pass nullptr to suppress.
 * @param magLabel        Label used in the pave text for the vector magnitude.
 * @param minEntries      Minimum per-bin entry count to include in a tile.
 * @param arrowBlockSize  Side length of the averaging tile in bins (how many tiles in the grid to use for averaging the polarization arrows).
 * @param scalePercentile Fraction of tile magnitudes used as the length-normalisation reference.
 */
// ==========================================================================
static void DrawVectorFieldPanel(TProfile2D* hX, TProfile2D* hY, TProfile2D* hZ,
                                  const char* title,
                                  const char* magLabel   = "|#LTp*_{T}#GT|",
                                  double minEntries      = 50.,
                                  int    arrowBlockSize  = 4,
                                  double scalePercentile = 0.95)
{
    if (!hX || !hY || !hZ) return;

    // ---- COLZ background ----
    TProfile2D* hZd = static_cast<TProfile2D*>(hZ->Clone(Form("hVFZ_tmp_%d", gCloneIdx++)));
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
    int pIdx = static_cast<int>(scalePercentile * static_cast<double>(validTiles.size() - 1));
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
                scaleRefErr = t->err_mag;
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

        double drawLen = std::min(t.mag, scaleRef) * scale;
        double x2 = t.xc + (t.bx / t.mag) * drawLen;
        double y2 = t.yc + (t.by / t.mag) * drawLen;

        TArrow* arr = new TArrow(t.xc, t.yc, x2, y2, 0.012, ">");
        arr->SetLineColor(kBlack);
        arr->SetFillColor(kBlack);
        arr->SetLineWidth(2);
        arr->Draw();
    }

    // ----- Compute plot-coordinate placement -----
    double xMin = hZd->GetXaxis()->GetXmin();
    double xMax = hZd->GetXaxis()->GetXmax();
    double yMin = hZd->GetYaxis()->GetXmin();
    double yMax = hZd->GetYaxis()->GetXmax();
    double x1 = xMin + 0.05 * (xMax - xMin);
    double x2 = xMin + 0.65 * (xMax - xMin);
    double y1 = yMax - 0.16 * (yMax - yMin);
    double y2 = yMax - 0.05 * (yMax - yMin);

    // ----- Background box -----
    TPaveText* pave = new TPaveText(x1, y1, x2, y2, "arc");
    pave->SetCornerRadius(0.15);
    pave->SetFillColor(kWhite);
    pave->SetFillStyle(1001);
    pave->SetBorderSize(0);
    pave->SetMargin(0.02);
    pave->SetTextAlign(12);
    pave->SetTextFont(63);
    pave->SetTextSize(18);

    double percentile = scalePercentile * 100.;
    pave->AddText(Form("%s_{%.0fpct} = (%.2f #pm %.2f)%%", magLabel, percentile, scaleRef * 100., scaleRefErr * 100.));
    pave->Draw();
}

// ==========================================================================
/**
 * @brief Draws a single TProfile2D as a COLZ map with a symmetric
 * (diverging) Z range derived from the data -- the same background-map
 * logic used inside DrawVectorFieldPanel() above, kept as a separate small
 * helper since the ring-observable panels need only the background, no
 * arrows.
 *
 * @param hSrc            Source TProfile2D; left untouched (a scaled clone is drawn).
 * @param title           Full ROOT title string (";xTitle;yTitle;zTitle"); nullptr keeps the source title.
 * @param minEntries      Minimum per-bin entry count for a bin to count toward the Z-range scan.
 * @param scaleToPercent  If true, multiplies bin content by 100 before drawing (fraction --> %).
 */
// ==========================================================================
static void DrawSymmetricColz(TProfile2D* hSrc, const char* title, double minEntries = 50., bool scaleToPercent = true)
{
    if (!hSrc) return;
    TProfile2D* hd = static_cast<TProfile2D*>(hSrc->Clone(Form("hColz_tmp_%d", gCloneIdx++)));
    hd->SetDirectory(nullptr);
    hd->SetStats(0);
    if (scaleToPercent) hd->Scale(100.0);
    if (title) hd->SetTitle(title);

    double zmax = 0.;
    for (int ix = 1; ix <= hd->GetNbinsX(); ++ix)
        for (int iy = 1; iy <= hd->GetNbinsY(); ++iy) {
            int gb = hd->GetBin(ix, iy);
            if (hd->GetBinEntries(gb) < minEntries) continue;
            double v = std::fabs(hd->GetBinContent(ix, iy));
            if (v > zmax) zmax = v;
        }
    zmax = (zmax < 1.e-10) ? 1.e-4 : zmax * 1.1;
    hd->SetMinimum(-zmax);
    hd->SetMaximum(zmax);
    hd->Draw("COLZ");
}

// ==========================================================================
// Folder registry: the four kinematic-cut scenarios booked by lambdaJetPolarizationIonsDerived.cxx's
// addRingObservableFamily() lambda (Ring, RingKinematicCuts, JetKinematicCuts, JetAndLambdaKinematicCuts).
// If the consumer ever receives another folder, you can just add it here! At the same time,
// if the consumer switch is off for any of these three families in the source .root file, safeGet() and
// GetDir() will make the program fallback into a non-fatal crash (i.e., it just ignores these families).
// ==========================================================================
struct FolderSpec {
    const char* name;  ///< O2 histogram-registry folder name (matches addRingObservableFamily(...) argument)
    const char* label; ///< Human-readable label used in canvas titles
};

static const std::vector<FolderSpec> kFolders = {
    {"Ring",                      "Ring (no kinematic cuts)"},
    {"RingKinematicCuts",         "Ring, #Lambda kinematic cuts"},
    {"JetKinematicCuts",          "Ring, jet kinematic cuts"},
    {"JetAndLambdaKinematicCuts", "Ring, jet & #Lambda kinematic cuts"},
};

static const char* kTaskDir = "lambdajetpolarizationionsderived"; // My O2Physics task name
// ==========================================================================
/**
 * @brief One canvas per folder: X-Y plane and Z-X plane polarization vector
 * fields side by side, reading the p2dP{x,y,z}Star_vsPxPy /
 * p2dP{x,y,z}Star_vsPzPx TProfile2Ds booked under <folder>/QA/ by the
 * consumer.
 *
 * Left panel (X-Y plane):  colormap = <P*_z> (out-of-plane), arrows = (<P*_x>, <P*_y>).
 * Right panel (Z-X plane): colormap = <P*_y> (out-of-plane), arrows sourced as
 * (hPxZX, hPzZX, hPyZX) -- this argument order matches the established
 * convention in plotHelicityEfficiency.cxx's own (pz,px)-plane calls; see
 * that file's DrawVectorFieldPanel() call sites.
 *
 * @param taskDir  Top-level task TDirectory ("lambdajetpolarizationionsderived").
 * @param outDir   Output sub-directory to write the canvas into.
 * @param f        Folder to process (see kFolders).
 */
// ==========================================================================
static void MakeVectorFieldCanvas(TDirectory* taskDir, TDirectory* outDir, const FolderSpec& f)
{
    TDirectory* qaDir = GetDir(GetDir(taskDir, f.name), "QA");
    if (!qaDir) { printf("WARNING: skipping vector-field canvas for '%s' (QA dir missing)\n", f.name); return; }

    TProfile2D* hPxXY = static_cast<TProfile2D*>(SafeGet(qaDir, "p2dPxStar_vsPxPy"));
    TProfile2D* hPyXY = static_cast<TProfile2D*>(SafeGet(qaDir, "p2dPyStar_vsPxPy"));
    TProfile2D* hPzXY = static_cast<TProfile2D*>(SafeGet(qaDir, "p2dPzStar_vsPxPy"));
    TProfile2D* hPxZX = static_cast<TProfile2D*>(SafeGet(qaDir, "p2dPxStar_vsPzPx"));
    TProfile2D* hPyZX = static_cast<TProfile2D*>(SafeGet(qaDir, "p2dPyStar_vsPzPx"));
    TProfile2D* hPzZX = static_cast<TProfile2D*>(SafeGet(qaDir, "p2dPzStar_vsPzPx"));
    if (!hPxXY || !hPyXY || !hPzXY || !hPxZX || !hPyZX || !hPzZX) {
        printf("WARNING: skipping vector-field canvas for '%s' (one or more histograms missing)\n", f.name);
        return;
    }

    TCanvas* c = new TCanvas(Form("cVectorField_%s", f.name), "", 1400, 650);
    c->Divide(2, 1, 0.005, 0.005);

    c->cd(1);
    gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.16); gPad->SetTopMargin(0.08); gPad->SetBottomMargin(0.10);
    hPzXY->GetZaxis()->SetTitleOffset(1.5);
    DrawVectorFieldPanel(hPxXY, hPyXY, hPzXY, " ;p_{x}^{#Lambda} [GeV/c];p_{y}^{#Lambda} [GeV/c];<P*_{z}> [%]");
    AddLabel(0.5, 0.95, "X-Y plane", 0.045, 22);

    c->cd(2);
    gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.16); gPad->SetTopMargin(0.08); gPad->SetBottomMargin(0.10);
    hPyZX->GetZaxis()->SetTitleOffset(1.5);
    DrawVectorFieldPanel(hPzZX, hPxZX, hPyZX, " ;p_{z}^{#Lambda} [GeV/c];p_{x}^{#Lambda} [GeV/c];<P*_{y}> [%]");
    AddLabel(0.5, 0.95, "Z-X plane", 0.045, 22);

    c->cd(0);
    AddLabel(0.5, 0.98, Form("<P*> polarization vector field  --  %s", f.label), 0.032, 22);

    WriteCanvas(c, outDir);
    delete c;
}

// ==========================================================================
/**
 * @brief One canvas per folder: X-Y plane and Z-X plane ring-observable
 * COLZ maps, reading p2dRingObservableVsPxPy / p2dRingObservableVsPzPx
 * (booked directly under <folder>/, NOT under QA/). No arrow overlay is
 * needed here -- the ring observable is already a single scalar per
 * candidate, so the COLZ background alone is the full result.
 *
 * @param taskDir  Top-level task TDirectory.
 * @param outDir   Output sub-directory to write the canvas into.
 * @param f        Folder to process (see kFolders).
 */
// ==========================================================================
static void MakeRingObservable2DCanvas(TDirectory* taskDir, TDirectory* outDir, const FolderSpec& f)
{
    TDirectory* dir = GetDir(taskDir, f.name);
    if (!dir) { printf("WARNING: skipping ring-observable canvas for '%s' (folder missing)\n", f.name); return; }

    TProfile2D* hXY = static_cast<TProfile2D*>(SafeGet(dir, "p2dRingObservableVsPxPy"));
    TProfile2D* hZX = static_cast<TProfile2D*>(SafeGet(dir, "p2dRingObservableVsPzPx"));
    if (!hXY || !hZX) {
        printf("WARNING: skipping ring-observable canvas for '%s' (one or both histograms missing)\n", f.name);
        return;
    }

    TCanvas* c = new TCanvas(Form("cRingObservable2D_%s", f.name), "", 1400, 650);
    c->Divide(2, 1, 0.005, 0.005);

    c->cd(1);
    gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.16); gPad->SetTopMargin(0.08); gPad->SetBottomMargin(0.10);
    DrawSymmetricColz(hXY, " ;p_{x}^{#Lambda} [GeV/c];p_{y}^{#Lambda} [GeV/c];<#it{R}> [%]");
    AddLabel(0.5, 0.95, "X-Y plane", 0.045, 22);

    c->cd(2);
    gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.16); gPad->SetTopMargin(0.08); gPad->SetBottomMargin(0.10);
    DrawSymmetricColz(hZX, " ;p_{z}^{#Lambda} [GeV/c];p_{x}^{#Lambda} [GeV/c];<#it{R}> [%]");
    AddLabel(0.5, 0.95, "Z-X plane", 0.045, 22);

    c->cd(0);
    AddLabel(0.5, 0.98, Form("Ring observable <#it{R}>  --  %s", f.label), 0.032, 22);

    WriteCanvas(c, outDir);
    delete c;
}

// ---------------------------------------------------------
// Main Function
// ---------------------------------------------------------
int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <inputFilePath (ConsumerResults_*.root)>\n";
        return 1;
    }

    const char* inFileStr = argv[1];

    // --- Dynamic Output Filename (mirrors makeCumulativeDCAdauProfile.cxx) ---
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
        filename.replace(pos, prefixToReplace.length(), "AuxiliaryPerConfigPlots_");
    } else {
        filename = "AuxiliaryPerConfigPlots_" + filename;
    }

    std::string outFileStr = directory + filename;

    std::cout << "\n=======================================================\n";
    std::cout << " Starting Auxiliary Per-Config Plots\n";
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
    if (!fOut || fOut->IsZombie()) {
        std::cerr << " ERROR: cannot create " << outFileStr << std::endl;
        fIn->Close();
        return 1;
    }

    // --- Global style, matching plotHelicityEfficiency.cxx ---
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetTitleSize(0.05, "XYZ");
    gStyle->SetLabelSize(0.045, "XYZ");
    gStyle->SetTitleOffset(1.1, "Y");
    // gStyle->SetPalette(kTemperatureMap); // diverging blue-white-red, zero = white -- flip on if preferred over the default palette

    TDirectory* taskDir = GetDir(fIn, kTaskDir);
    if (!taskDir) {
        std::cerr << " ERROR: '" << kTaskDir << "' directory not found in " << inPath << std::endl;
        fOut->Close(); fIn->Close();
        return 1;
    }

    // --- Create organized output subfolders ---
    TDirectory* dirVectorField = fOut->mkdir("VectorField_Canvases");
    TDirectory* dirRingObs2D   = fOut->mkdir("RingObservable2D_Canvases");

    // 1. Polarization vector-field canvases (one per folder)
    std::cout << " -> Drawing vector-field canvases...\n";
    for (const auto& f : kFolders)
        MakeVectorFieldCanvas(taskDir, dirVectorField, f);

    // 2. Ring-observable 2D scalar canvases (one per folder)
    std::cout << " -> Drawing ring-observable 2D canvases...\n";
    for (const auto& f : kFolders)
        MakeRingObservable2DCanvas(taskDir, dirRingObs2D, f);

    // ===========================================================
    // ADD MORE POST-PROCESSING SECTIONS HERE.
    //
    // This file is meant to be the general home for per-config
    // derivative plots that live outside O2Physics -- i.e. anything
    // that reads a finished ConsumerResults_*.root and produces
    // additional canvases/histograms from it, but does not need the
    // cross-config/wagon-level aggregation that auxiliarySummaryPlots.cxx
    // handles. Follow the same pattern as sections 1-2 above:
    //   - a FolderSpec-style loop if the new plot is booked per
    //     kinematic-cut folder (reuse kFolders if so),
    //   - a dedicated TDirectory in fOut via fOut->mkdir(...),
    //   - a Make...Canvas(taskDir, outDir, ...) function defined above
    //     main(), called from a loop right here.
    // ===========================================================

    std::cout << " Successfully generated and saved auxiliary per-config plots.\n";
    std::cout << "=======================================================\n\n";

    fOut->Write("", TObject::kOverwrite);
    fOut->Close();
    fIn->Close();

    delete fOut;
    delete fIn;

    return 0;
}