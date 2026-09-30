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
*   1. Polarization vector fields (the "pretty arrows" plots), one canvas per coordinate system:
*        Lab       -- the detector frame (three planes: X-Y, Z-X, Y-Z)
*        PrimeJet  -- rotated onto the jet axis (three planes). The X'-Y' panel is the headline plot:
*                     the ring observable is exactly the azimuthal component of the arrows there.
*        PrimeV0   -- the Lambda production plane, (p_z, p_T), plus its counts map
*   2. Ring observable 2D scalar maps, one canvas per coordinate system and per jet proxy.
*      No arrow overlay: the ring observable is already one scalar per candidate.
*   3. AEE acceptance maps: counts, <P*_T> and <P*_z> over the three Aee planes.
*   4. Candidate counts over the Lab and PrimeJet planes.
*   5. AEE fold: the odd-in-Phi_AEE component of the acceptance -- which IS the AEE --
*      together with the three parity nulls that accompany it.
*   6. Ring projection kernel: per-slice linear fits of <R> against the proxy direction cosine,
*      which isolate the acceptance moment M0 from both the A/C yield asymmetry and the signal.
*   7. Kernel moments: B_phi recovered differentially in cos(theta_Lambda) from the 3D kernel
*      profiles, and the M0 / M1 moments built from it -- including the eta-symmetrised M1.
*   8. Kernel symmetry decomposition: the 3D kernel split into the four sectors of the group generated
*      by the z-reflection and the antipodal reference axis, which separates detector-induced from
*      genuine contributions and leaves one sector as a null test.
*   9. KappaEff: the response coefficient kappa per mass bin, per proxy (both ring definitions).
*  10. R_z only (useRingZ files, "_useRingZ" in the name): the Delta phi fold (east/west), the chi symmetry
*      sectors, and the Phi_AEE flatness test. In R_z files Sections 6-8 are null tests.
*
* All sections are produced once per kinematic-cut folder (see kFolders below): Ring, RingKinematicCuts,
* JetKinematicCuts, JetAndLambdaKinematicCuts.
*
* The geometry of the three coordinate systems, the reason each one exists, and the parity checks the
* resulting maps must satisfy are all in the folder README ("Coordinate systems used by the
* polarization maps"). They are NOT repeated here.
*
* ADDING A NEW SECTION: see the "ADD MORE POST-PROCESSING SECTIONS HERE" comment inside main() below.
*
* ON THE DUPLICATED DrawVectorFieldPanel():
* DrawVectorFieldPanel() below is copied close to verbatim from plotHelicityEfficiency.cxx (the toy-model plotter).
* This is a deliberate choice: the toy-model plotter and this per-config post-processor are two independent workflows,
* compiled and ran separately, and the function itself is small, self-contained, and should not change that often.
* Sharing it via a common header would couple two otherwise-unrelated build targets for little benefit.
* If it ever needs to diverge meaningfully or grows non-trivially, that is the trigger to revisit and actually factor it out.
*
* Usage mirrors makeCumulativeDCAdauProfile.cxx; compiled ahead-of-time by run_all_wagons.sh, not invoked as a ROOT interpreted macro.
* In case you do want to run it by hand anyways, command would look like this:
*   ./auxiliaryPerConfigPlots.exe path/to/ConsumerResults_SUFFIX.root path/to/outputDirectory/
* (outputDirectory is created recursively if it does not already exist -- this keeps
*  post-processing output out of results_consumer/, see run_all_wagons.sh)
* ============================================================
*/

#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TArrow.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TProfile3D.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TLine.h>
#include <TMath.h>

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

// warnIfMissing is a functionality to avoid printing warnings in some cases. It is not always an anomaly to have some plots
// missing: the consumer may book only some of the kinematic cut families at once, so this is expected behavior in most cases.
// Pass false when the caller is doing the probing and will report the outcome itself.
// Leave it true when the object is genuinely expected to be there.
static TObject* SafeGet(TDirectory* dir, const char* name, bool warnIfMissing = true)
{
    if (!dir) {
        if (warnIfMissing) printf("WARNING: SafeGet called with null directory for '%s'\n", name);
        return nullptr;
    }
    TObject* obj = dir->Get(name);
    if (!obj && warnIfMissing) printf("WARNING: '%s' not found in directory '%s'\n", name, dir->GetName());
    return obj;
}

static TDirectory* GetDir(TDirectory* parent, const char* name, bool warnIfMissing = true)
{
    if (!parent) {
        if (warnIfMissing) printf("WARNING: GetDir called with null parent for '%s'\n", name);
        return nullptr;
    }
    TDirectory* dir = static_cast<TDirectory*>(parent->Get(name));
    if (!dir && warnIfMissing) printf("WARNING: directory '%s' not found in '%s'\n", name, parent->GetName());
    return dir;
}

// Like GetDir, but resolves a multi-level path ("PolMaps/PrimeJet") in one call. An empty path returns
// the parent itself, which is what the ring-observable maps need: they live at the folder root.
static TDirectory* GetDirPath(TDirectory* parent, const std::string& path, bool warnIfMissing = true)
{
    if (!parent) {
        if (warnIfMissing) printf("WARNING: GetDirPath called with null parent for '%s'\n", path.c_str());
        return nullptr;
    }
    if (path.empty()) return parent;
    TDirectory* dir = parent->GetDirectory(path.c_str());
    if (!dir && warnIfMissing) printf("WARNING: directory '%s' not found in '%s'\n", path.c_str(), parent->GetName());
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

// Every panel in this file uses the same margins. The right margin is the wide one: it has to fit the
// COLZ palette and its axis title.
static void SetupPadMargins()
{
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.16);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.10);
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
 * documentation lives there -- EXCEPT for the axis-aspect correction below.
 *
 * THE ASPECT CORRECTION, and why it matters:
 * The arrow direction (bx, by) lives in polarization-component space, where both
 * components carry the same units, so the drawn angle is meaningful. The arrow
 * tip, however, is placed in DATA coordinates. Displacing by the same number of
 * data units along both axes therefore only renders the correct angle when both
 * axes map the same number of data units per pixel. They generally do not: a
 * (p_z, p_x) panel spans [-4,4] horizontally and [-3,3] vertically on a pad that
 * is not square, so a 45-degree polarization drew at roughly 34 degrees. Since
 * the whole point of these panels is reading directions off the page, that was a
 * silent, systematic misreading of every non-square panel.
 * The fix converts the displacement through the pad's pixel geometry, so the
 * arrow has a fixed pixel LENGTH and the true polarization ANGLE on screen.
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
    gPad->Update(); // Required: the pad geometry read below is only valid once the pad has been laid out

    // ---- Tile grid parameters ----
    int    bs      = arrowBlockSize;
    int    nBinsX  = hX->GetNbinsX();
    int    nBinsY  = hX->GetNbinsY();
    double binW    = hX->GetXaxis()->GetBinWidth(1);
    double tileW   = bs * binW; // Tile width, in data units of the HORIZONTAL axis
    int    nTilesX = nBinsX / bs;
    int    nTilesY = nBinsY / bs;

    // ---- Axis aspect: data units per pixel, vertical over horizontal ----
    // The frame occupies the pad minus its margins, so the pixel extents are taken from the pad size
    // times the margin-free fractions. aspect = 1 recovers the old (incorrect) behavior exactly, which
    // is what a square panel with equal axis ranges would give anyway.
    const double frameWpx = gPad->GetWw() * gPad->GetAbsWNDC() * (1.0 - gPad->GetLeftMargin() - gPad->GetRightMargin());
    const double frameHpx = gPad->GetWh() * gPad->GetAbsHNDC() * (1.0 - gPad->GetTopMargin() - gPad->GetBottomMargin());
    const double xRange   = hZd->GetXaxis()->GetXmax() - hZd->GetXaxis()->GetXmin();
    const double yRange   = hZd->GetYaxis()->GetXmax() - hZd->GetYaxis()->GetXmin();

    double aspect = 1.0;
    if (frameWpx > 0. && frameHpx > 0. && xRange > 0. && yRange > 0.) {
        const double dataPerPxX = xRange / frameWpx;
        const double dataPerPxY = yRange / frameHpx;
        aspect = dataPerPxY / dataPerPxX;
    }

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

    // full-reference arrow = 0.6 tile widths, measured along the HORIZONTAL axis
    double scale = 0.6 * tileW / scaleRef;

    // ---- Pass 2: draw arrows, capping length at the scale reference ----
    for (const Tile& t : tiles) {
        if (!t.valid) continue;
        if (t.mag < 0.05 * scaleRef) continue; // suppress near-zero noise

        // drawLen is the tip displacement expressed in HORIZONTAL data units; the vertical displacement
        // is the same physical length on screen, which is what the aspect factor converts it into.
        double drawLen = std::min(t.mag, scaleRef) * scale;
        double x2 = t.xc + (t.bx / t.mag) * drawLen;
        double y2 = t.yc + (t.by / t.mag) * drawLen * aspect;

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
/**
 * @brief Draws a TH2D occupancy map as a plain COLZ, floor pinned at zero.
 *
 * Counts are a one-sided quantity, so the diverging/symmetric treatment used for
 * the polarization maps would waste half the palette and put the empty regions in
 * mid-scale. This is the plot that carries the azimuthal efficiency modulation on
 * the AEE planes, so it deserves the full dynamic range.
 *
 * @param hSrc  Source TH2D; left untouched (a clone is drawn).
 * @param title Full ROOT title string (";xTitle;yTitle;zTitle"); nullptr keeps the source title.
 */
// ==========================================================================
static void DrawCountsColz(TH2* hSrc, const char* title)
{
    if (!hSrc) return;
    TH2* hd = static_cast<TH2*>(hSrc->Clone(Form("hCounts_tmp_%d", gCloneIdx++)));
    hd->SetDirectory(nullptr);
    hd->SetStats(0);
    if (title) hd->SetTitle(title);
    hd->SetMinimum(0.);
    hd->Draw("COLZ");
}


// ==========================================================================
// The AEE fold.
//
// On the (x_AEE, y_AEE) plane the radius is p_T^Lambda and the azimuth is Phi_AEE, so the reflection
// y_AEE --> -y_AEE is exactly Phi_AEE --> -Phi_AEE. Parity fixes how each map must behave under it
// (the full argument is in the folder README):
//
//     counts, <sin(theta*)>  -->  EVEN        <cos(theta*)>, i.e. <P*_z>  -->  ODD        <R>  -->  EVEN
//
// and the magnetic field is the one thing that breaks the relevant mirror. So the odd part of the
// counts map is not merely a diagnostic of the AEE -- it IS the AEE, isolated from acceptance that is
// non-uniform for ordinary geometric reasons. The even/odd parts of the remaining maps are then free
// null tests: each one should vanish, and a non-vanishing residue localises a problem.
//
// Folding requires the vertical axis to pair bin iy with bin N+1-iy, i.e. to be symmetric about zero
// with an even bin count. axisLambdaPRot is, but it is a ConfigurableAxis, so it is checked at runtime
// rather than assumed.
// ==========================================================================

/// @brief True when bin b of an axis is the mirror image of bin N+1-b about zero, for every b.
/// Checked edge by edge rather than from the range alone, so a variable-width ConfigurableAxis that
/// merely happens to span a symmetric range is still caught.
static bool IsMirrorPairable(TAxis* ax)
{
    if (!ax) return false;
    const int n = ax->GetNbins();
    if (n < 2 || n % 2 != 0) return false;
    const double span = ax->GetXmax() - ax->GetXmin();
    if (span <= 0.) return false;
    for (int b = 1; b <= n; ++b)
        if (std::fabs(ax->GetBinLowEdge(b) + ax->GetBinUpEdge(n + 1 - b)) > 1.e-6 * span) return false;
    return true;
}

/// @brief True when the Y axis pairs bins exactly under y --> -y (see IsMirrorPairable()).
static bool HasFoldableYAxis(TH1* h)
{
    return h && IsMirrorPairable(h->GetYaxis());
}

// ==========================================================================
/**
 * @brief Normalised count asymmetry under y --> -y.
 *
 * Returns \f$ A(x,y) = [N(x,y) - N(x,-y)] / [N(x,y) + N(x,-y)] \f$, which is dimensionless, bounded in
 * [-1, 1], and reads directly as the fractional efficiency modulation. For independent Poisson counts
 * \f$a\f$ and \f$b\f$ the propagated variance is exactly \f$ 4ab/(a+b)^3 \f$.
 *
 * The result is filled over the FULL plane rather than a half-plane: it is antisymmetric by
 * construction, so drawing both halves makes that manifest and costs nothing but redundancy.
 * Bin pairs with no entries at all are left empty (content and error zero).
 */
// ==========================================================================
static TH2D* MakeCountsAsymmetryY(TH2* src, const char* newName)
{
    if (!src || !HasFoldableYAxis(src)) return nullptr;

    TH2D* out = static_cast<TH2D*>(src->Clone(newName));
    out->SetDirectory(nullptr);
    out->SetStats(0);
    out->Reset();

    const int nx = src->GetNbinsX();
    const int ny = src->GetNbinsY();

    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            const int iyMirror = ny + 1 - iy;
            const double a = src->GetBinContent(ix, iy);
            const double b = src->GetBinContent(ix, iyMirror);
            const double sum = a + b;
            if (sum <= 0.) continue;

            out->SetBinContent(ix, iy, (a - b) / sum);
            out->SetBinError(ix, iy, std::sqrt(4.0 * a * b / (sum * sum * sum)));
        }
    }
    return out;
}

// ==========================================================================
/**
 * @brief Even or odd part of a TProfile2D under y --> -y.
 *
 * Returns \f$ \tfrac{1}{2}[P(x,y) \pm P(x,-y)] \f$ with error \f$ \tfrac{1}{2}\sqrt{e_+^2 + e_-^2} \f$.
 * A bin pair contributes only when BOTH partners clear minEntries, so a half-populated pair cannot
 * masquerade as a parity violation -- which is the failure mode this plot exists to detect.
 *
 * @param odd true for the antisymmetric part, false for the symmetric one.
 */
// ==========================================================================
static TH2D* MakeProfileFoldY(TProfile2D* src, bool odd, double minEntries, const char* newName)
{
    if (!src || !HasFoldableYAxis(src)) return nullptr;

    const int nx = src->GetNbinsX();
    const int ny = src->GetNbinsY();

    // Built as a plain TH2D: the folded quantity is a combination of two means, not a profile any more.
    TH2D* out = new TH2D(newName, "",
                         nx, src->GetXaxis()->GetXmin(), src->GetXaxis()->GetXmax(),
                         ny, src->GetYaxis()->GetXmin(), src->GetYaxis()->GetXmax());
    out->SetDirectory(nullptr);
    out->SetStats(0);

    const double sign = odd ? -1.0 : +1.0;

    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            const int iyMirror = ny + 1 - iy;
            if (src->GetBinEntries(src->GetBin(ix, iy)) < minEntries) continue;
            if (src->GetBinEntries(src->GetBin(ix, iyMirror)) < minEntries) continue;

            const double p = src->GetBinContent(ix, iy);
            const double m = src->GetBinContent(ix, iyMirror);
            const double ep = src->GetBinError(ix, iy);
            const double em = src->GetBinError(ix, iyMirror);

            out->SetBinContent(ix, iy, 0.5 * (p + sign * m));
            out->SetBinError(ix, iy, 0.5 * std::sqrt(ep * ep + em * em));
        }
    }
    return out;
}

// ==========================================================================
/**
 * @brief COLZ with a symmetric Z range, for a plain TH2 rather than a TProfile2D.
 *
 * Deliberately separate from DrawSymmetricColz(): that one gates the Z-range scan on
 * TProfile2D::GetBinEntries(), which a TH2D does not have. Here an unfilled bin is identified by a
 * zero error instead, which is exactly what the fold helpers above leave behind.
 */
// ==========================================================================
static void DrawSymmetricColzH2(TH2* hSrc, const char* title, bool scaleToPercent = true)
{
    if (!hSrc) return;
    TH2* hd = static_cast<TH2*>(hSrc->Clone(Form("hFoldColz_tmp_%d", gCloneIdx++)));
    hd->SetDirectory(nullptr);
    hd->SetStats(0);
    if (scaleToPercent) hd->Scale(100.0);
    if (title) hd->SetTitle(title);

    double zmax = 0.;
    for (int ix = 1; ix <= hd->GetNbinsX(); ++ix)
        for (int iy = 1; iy <= hd->GetNbinsY(); ++iy) {
            if (hd->GetBinError(ix, iy) <= 0.) continue; // never filled by the fold
            const double v = std::fabs(hd->GetBinContent(ix, iy));
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
// If the consumer ever receives another folder, you can just add it here!
//
// Only "Ring" is generally booked in all executions. 
// The absence of the other three is a design feature, so it is resolved ONCE by ScanPresentFolders() below and reported as a single line,
// rather than being rediscovered (and re-warned about) by every drawing section in turn (this was polluting the logs way too much!).
// ==========================================================================
struct FolderSpec {
    const char* name;      ///< O2 histogram-registry folder name (matches addRingObservableFamily(...) argument in the consumer)
    const char* label;     ///< Human-readable label used in canvas titles
    bool        mandatory; ///< true -> its absence means the input file is broken, not merely trimmed
};

static const std::vector<FolderSpec> kFolders = {
    {"Ring",                      "Ring (no kinematic cuts)",           true },
    {"RingKinematicCuts",         "Ring, #Lambda kinematic cuts",       false},
    {"JetKinematicCuts",          "Ring, jet kinematic cuts",           false},
    {"JetAndLambdaKinematicCuts", "Ring, jet & #Lambda kinematic cuts", false},
};

static const char* kTaskDir = "lambdajetpolarizationionsderived"; // My O2Physics task name

// Ring definition of the input file. A useRingZ consumer output is identified by its file name alone: the
// config names carry "_useRingZ" right after the family (BothHyperons_useRingZ_MixedEventProxies), so the tag
// is searched anywhere in the basename (see the README, "The longitudinal ring"). main() sets both globals.
static const char* kRingZTag = "_useRingZ";
static std::string gRingSym  = "#it{R}"; // "#it{R}_{z}" for R_z files
static std::string gNullTag  = "";       // Prepended to the kernel titles, which are null tests for R_z

/// @brief Swaps the full-ring symbol for the one of the current file (a no-op for full-ring files).
static std::string RingLabel(std::string text)
{
    static const std::string kFull = "#it{R}";
    for (size_t pos = text.find(kFull); pos != std::string::npos; pos = text.find(kFull, pos + gRingSym.size()))
        text.replace(pos, kFull.size(), gRingSym);
    return text;
}

// ==========================================================================
// Panel registry.
//
// Every canvas in this file is a horizontal strip of panels, and every panel is one of three kinds.
// Describing them as data instead of as code means a new coordinate system costs one table, not one
// function: the builder below handles fetching, missing-object reporting, layout and labelling.
//
// Axis-title fragments are named once here and reused by every table, so a relabelled axis cannot
// drift between the Lab, Aee, PrimeV0 and PrimeJet versions of the same plane.
// ==========================================================================
enum class PanelKind {
    kVectorField, ///< hA, hB = arrow components (horizontal, vertical); hC = COLZ background. All TProfile2D.
    kProfileColz, ///< hA = TProfile2D drawn as a diverging COLZ map. hB, hC unused.
    kCountsColz   ///< hA = TH2D drawn as a one-sided COLZ occupancy map. hB, hC unused.
};

struct PanelSpec {
    PanelKind   kind;
    std::string hA;
    std::string hB;
    std::string hC;
    std::string axes;     ///< Full ROOT title string " ;xTitle;yTitle;zTitle"
    std::string banner;   ///< Short text drawn above the panel
    std::string magLabel; ///< kVectorField only; empty falls back to the DrawVectorFieldPanel default
};

// --- Axis-title fragments ---
static const std::string kPx    = "p_{x}^{#Lambda} [GeV/c]";
static const std::string kPy    = "p_{y}^{#Lambda} [GeV/c]";
static const std::string kPz    = "p_{z}^{#Lambda} [GeV/c]";
static const std::string kPt    = "p_{T}^{#Lambda} [GeV/c]";
static const std::string kPxA   = "p_{x,AEE}^{#Lambda} [GeV/c]";
static const std::string kPyA   = "p_{y,AEE}^{#Lambda} [GeV/c]";
static const std::string kPxJ   = "p_{x'Jet}^{#Lambda} [GeV/c]";
static const std::string kPyJ   = "p_{y'Jet}^{#Lambda} [GeV/c]";
static const std::string kPzJ   = "p_{z'Jet}^{#Lambda} [GeV/c]";
static const std::string kZCnt  = "Counts";

// Builds the " ;x;y;z" title string the ROOT drawing helpers expect.
static std::string Axes(const std::string& x, const std::string& y, const std::string& z)
{
    return " ;" + x + ";" + y + ";" + z;
}

// ==========================================================================
/**
 * @brief Builds one canvas from a panel table: fetches, checks, lays out, draws and writes.
 *
 * Fetching is all-or-nothing on purpose. f is already known to exist (main() only iterates the
 * folders ScanPresentFolders() resolved), so a missing histogram here is a genuine anomaly --
 * most likely a consumer/macro name drift -- and half a canvas would hide it rather than show it.
 * The warning names the first offender so the drift is one grep away.
 *
 * @param taskDir    Top-level task TDirectory ("lambdajetpolarizationionsderived").
 * @param outDir     Output sub-directory to write the canvas into.
 * @param f          Folder to process (see kFolders).
 * @param subPath    Directory inside the folder holding the histograms; "" means the folder root.
 * @param canvasName Canvas name prefix; the folder name is appended.
 * @param mainTitle  Title drawn across the top of the canvas.
 * @param panels     The panel table.
 */
// ==========================================================================
static void MakePanelCanvas(TDirectory* taskDir, TDirectory* outDir, const FolderSpec& f,
                            const std::string& subPath, const std::string& canvasName,
                            const std::string& mainTitle, const std::vector<PanelSpec>& panels)
{
    if (panels.empty()) return;

    TDirectory* folderDir = GetDir(taskDir, f.name);
    if (!folderDir) { printf("WARNING: skipping '%s' for '%s' (folder missing)\n", canvasName.c_str(), f.name); return; }

    TDirectory* dir = GetDirPath(folderDir, subPath);
    if (!dir) { printf("WARNING: skipping '%s' for '%s' (subdirectory '%s' missing)\n", canvasName.c_str(), f.name, subPath.c_str()); return; }

    // ---- Resolve every object first, so a partial canvas is never written ----
    struct Resolved { TProfile2D* a; TProfile2D* b; TProfile2D* c; TH2* counts; };
    std::vector<Resolved> objs(panels.size(), Resolved{nullptr, nullptr, nullptr, nullptr});

    for (size_t i = 0; i < panels.size(); ++i) {
        const PanelSpec& p = panels[i];
        const char* missing = nullptr;

        if (p.kind == PanelKind::kCountsColz) {
            objs[i].counts = static_cast<TH2*>(SafeGet(dir, p.hA.c_str(), false));
            if (!objs[i].counts) missing = p.hA.c_str();
        } else {
            objs[i].a = static_cast<TProfile2D*>(SafeGet(dir, p.hA.c_str(), false));
            if (!objs[i].a) missing = p.hA.c_str();
            if (p.kind == PanelKind::kVectorField) {
                objs[i].b = static_cast<TProfile2D*>(SafeGet(dir, p.hB.c_str(), false));
                objs[i].c = static_cast<TProfile2D*>(SafeGet(dir, p.hC.c_str(), false));
                if (!missing && !objs[i].b) missing = p.hB.c_str();
                if (!missing && !objs[i].c) missing = p.hC.c_str();
            }
        }

        if (missing) {
            printf("WARNING: skipping '%s' for '%s' (missing '%s' in '%s')\n",
                   canvasName.c_str(), f.name, missing, subPath.empty() ? f.name : subPath.c_str());
            return;
        }
    }

    // ---- Layout: one column per panel, 700 px each ----
    const int nPanels = static_cast<int>(panels.size());
    TCanvas* c = new TCanvas(Form("%s_%s", canvasName.c_str(), f.name), "", 700 * nPanels, 650);
    c->Divide(nPanels, 1, 0.005, 0.005);

    for (int i = 0; i < nPanels; ++i) {
        const PanelSpec& p = panels[i];
        c->cd(i + 1);
        SetupPadMargins();

        switch (p.kind) {
            case PanelKind::kVectorField:
                objs[i].c->GetZaxis()->SetTitleOffset(1.5);
                DrawVectorFieldPanel(objs[i].a, objs[i].b, objs[i].c, p.axes.c_str(),
                                     p.magLabel.empty() ? "|#LTp*_{T}#GT|" : p.magLabel.c_str());
                break;
            case PanelKind::kProfileColz:
                objs[i].a->GetZaxis()->SetTitleOffset(1.5);
                DrawSymmetricColz(objs[i].a, p.axes.c_str());
                break;
            case PanelKind::kCountsColz:
                objs[i].counts->GetZaxis()->SetTitleOffset(1.5);
                DrawCountsColz(objs[i].counts, p.axes.c_str());
                break;
        }
        AddLabel(0.5, 0.95, p.banner.c_str(), 0.045, 22);
    }

    c->cd(0);
    AddLabel(0.5, 0.98, Form("%s  --  %s", mainTitle.c_str(), f.label), 0.032, 22);

    WriteCanvas(c, outDir);
    delete c;
}

// ==========================================================================
// Panel tables.
// ==========================================================================

// --- Section 1a: lab-frame polarization vector fields ---
static const std::vector<PanelSpec> kPanelsVectorLab = {
    {PanelKind::kVectorField, "p2dPxStar_vsPxPy", "p2dPyStar_vsPxPy", "p2dPzStar_vsPxPy",
     Axes(kPx, kPy, "<P*_{z}> [%]"), "X-Y plane", ""},
    {PanelKind::kVectorField, "p2dPzStar_vsPzPx", "p2dPxStar_vsPzPx", "p2dPyStar_vsPzPx",
     Axes(kPz, kPx, "<P*_{y}> [%]"), "Z-X plane", ""},
    {PanelKind::kVectorField, "p2dPyStar_vsPyPz", "p2dPzStar_vsPyPz", "p2dPxStar_vsPyPz",
     Axes(kPy, kPz, "<P*_{x}> [%]"), "Y-Z plane", ""},
};

// --- Section 1b: jet-frame polarization vector fields ---
// The X'-Y' panel is the measurement itself: the ring observable is exactly the azimuthal component of
// these arrows about the jet axis, so a ring signal appears as arrows CIRCULATING around the origin
// while a radial pattern is something else. The plane radius is |p^Lambda| sin(DeltaTheta_jet).
static const std::vector<PanelSpec> kPanelsVectorPrimeJet = {
    {PanelKind::kVectorField, "p2dPxStarPrimeJet_vsPxPyPrimeJet", "p2dPyStarPrimeJet_vsPxPyPrimeJet", "p2dPzStarPrimeJet_vsPxPyPrimeJet",
     Axes(kPxJ, kPyJ, "<P*_{z'Jet}> [%]"), "X'-Y' plane (transverse to the jet)", "|#LTP*_{#perp Jet}#GT|"},
    {PanelKind::kVectorField, "p2dPzStarPrimeJet_vsPzPxPrimeJet", "p2dPxStarPrimeJet_vsPzPxPrimeJet", "p2dPyStarPrimeJet_vsPzPxPrimeJet",
     Axes(kPzJ, kPxJ, "<P*_{y'Jet}> [%]"), "Z'-X' plane (the beam-jet plane)", "|#LTP*_{#perp Jet}#GT|"},
    {PanelKind::kVectorField, "p2dPyStarPrimeJet_vsPyPzPrimeJet", "p2dPzStarPrimeJet_vsPyPzPrimeJet", "p2dPxStarPrimeJet_vsPyPzPrimeJet",
     Axes(kPyJ, kPzJ, "<P*_{x'Jet}> [%]"), "Y'-Z' plane", "|#LTP*_{#perp Jet}#GT|"},
};

// --- Section 1c: the Lambda production plane ---
// (p_z, p_T) IS the production plane, and (z, x'V0) are its own in-plane unit vectors, so axes and
// arrows share one basis here. The COLZ is the out-of-plane component, which is identically the ring
// observable computed with the beam as the jet proxy.
static const std::vector<PanelSpec> kPanelsVectorPrimeV0 = {
    {PanelKind::kVectorField, "p2dPzStar_vsPzPt", "p2dPxStarPrimeV0_vsPzPt", "p2dPyStarPrimeV0_vsPzPt",
     Axes(kPz, kPt, "<P*_{y'V0}> [%]"), "#Lambda production plane", "|#LTP*_{in-plane}#GT|"},
    {PanelKind::kCountsColz, "h2dCountsVsPzPt", "", "",
     Axes(kPz, kPt, kZCnt), "Candidate counts", ""},
};

// --- Section 3: AEE acceptance maps, one canvas per plane ---
// Radius on the first plane is p_T^Lambda and its azimuth is PhiAEE, so the counts panel there IS the
// azimuthal efficiency modulation resolved in p_T. Only two polarization observables exist in this
// frame, and that is not an omission: see the README.
static const std::vector<PanelSpec> kPanelsAeeXY = {
    {PanelKind::kCountsColz, "h2dCountsVsPxAeePyAee", "", "",
     Axes(kPxA, kPyA, kZCnt), "Candidate counts", ""},
    {PanelKind::kProfileColz, "p2dPtStar_vsPxAeePyAee", "", "",
     Axes(kPxA, kPyA, "<P*_{T}> [%]"), "<P*_{T}>  (i.e. sin#theta*)", ""},
    {PanelKind::kProfileColz, "p2dPzStar_vsPxAeePyAee", "", "",
     Axes(kPxA, kPyA, "<P*_{z}> [%]"), "<P*_{z}>  (i.e. cos#theta*)", ""},
};

static const std::vector<PanelSpec> kPanelsAeeZX = {
    {PanelKind::kCountsColz, "h2dCountsVsPzPxAee", "", "",
     Axes(kPz, kPxA, kZCnt), "Candidate counts", ""},
    {PanelKind::kProfileColz, "p2dPtStar_vsPzPxAee", "", "",
     Axes(kPz, kPxA, "<P*_{T}> [%]"), "<P*_{T}>  (i.e. sin#theta*)", ""},
    {PanelKind::kProfileColz, "p2dPzStar_vsPzPxAee", "", "",
     Axes(kPz, kPxA, "<P*_{z}> [%]"), "<P*_{z}>  (i.e. cos#theta*)", ""},
};

static const std::vector<PanelSpec> kPanelsAeeYZ = {
    {PanelKind::kCountsColz, "h2dCountsVsPyAeePz", "", "",
     Axes(kPyA, kPz, kZCnt), "Candidate counts", ""},
    {PanelKind::kProfileColz, "p2dPtStar_vsPyAeePz", "", "",
     Axes(kPyA, kPz, "<P*_{T}> [%]"), "<P*_{T}>  (i.e. sin#theta*)", ""},
    {PanelKind::kProfileColz, "p2dPzStar_vsPyAeePz", "", "",
     Axes(kPyA, kPz, "<P*_{z}> [%]"), "<P*_{z}>  (i.e. cos#theta*)", ""},
};

// --- Section 4: candidate counts over the lab and jet-frame planes ---
static const std::vector<PanelSpec> kPanelsCountsLab = {
    {PanelKind::kCountsColz, "h2dCountsVsPxPy", "", "", Axes(kPx, kPy, kZCnt), "X-Y plane", ""},
    {PanelKind::kCountsColz, "h2dCountsVsPzPx", "", "", Axes(kPz, kPx, kZCnt), "Z-X plane", ""},
    {PanelKind::kCountsColz, "h2dCountsVsPyPz", "", "", Axes(kPy, kPz, kZCnt), "Y-Z plane", ""},
};

static const std::vector<PanelSpec> kPanelsCountsPrimeJet = {
    {PanelKind::kCountsColz, "h2dCountsVsPxPyPrimeJet", "", "", Axes(kPxJ, kPyJ, kZCnt), "X'-Y' plane", ""},
    {PanelKind::kCountsColz, "h2dCountsVsPzPxPrimeJet", "", "", Axes(kPzJ, kPxJ, kZCnt), "Z'-X' plane", ""},
    {PanelKind::kCountsColz, "h2dCountsVsPyPzPrimeJet", "", "", Axes(kPyJ, kPzJ, kZCnt), "Y'-Z' plane", ""},
};

// ==========================================================================
/**
 * @brief Builds the three ring-observable panels for one coordinate system.
 *
 * The ring maps differ between jet proxies only by a histogram-name prefix and between coordinate
 * systems only by the plane suffixes and the axis labels, so they are generated rather than tabulated:
 * four near-identical hand-written tables would be four places for a typo to hide.
 *
 * @param prefix   Histogram-name prefix, e.g. "p2dRingObservable" or "p2dRingObservableLeadP".
 * @param suffixes The three plane suffixes, e.g. {"PxPy", "PzPx", "PyPz"}.
 * @param xTitles  Horizontal axis title of each plane, in the same order.
 * @param yTitles  Vertical axis title of each plane, in the same order.
 * @param banners  Panel banner of each plane, in the same order.
 */
// ==========================================================================
static std::vector<PanelSpec> MakeRingPanels(const std::string& prefix,
                                             const std::vector<std::string>& suffixes,
                                             const std::vector<std::string>& xTitles,
                                             const std::vector<std::string>& yTitles,
                                             const std::vector<std::string>& banners)
{
    std::vector<PanelSpec> panels;
    panels.reserve(suffixes.size());
    for (size_t i = 0; i < suffixes.size(); ++i) {
        panels.push_back({PanelKind::kProfileColz, prefix + "Vs" + suffixes[i], "", "",
                          Axes(xTitles[i], yTitles[i], RingLabel("<#it{R}> [%]")), banners[i], ""});
    }
    return panels;
}


// ==========================================================================
/**
 * @brief One canvas per folder: the AEE isolated, plus the three parity nulls that come with it.
 *
 * This one does not go through MakePanelCanvas(): it computes its panels instead of fetching them, and
 * it reads from two directories at once (the maps live under PolMaps/Aee, the ring map at the folder
 * root). The folded histograms are written alongside the canvas so the result can be used
 * quantitatively -- projected, integrated, fitted -- and not only looked at.
 *
 * | Panel | Content | Expectation |
 * |---|---|---|
 * | 1 | odd part of the counts map, normalised | **this is the AEE**; non-zero if present |
 * | 2 | odd part of <P*_T> | null: <sin(theta*)> is even in Phi_AEE |
 * | 3 | even part of <P*_z> | null: <cos(theta*)> is odd in Phi_AEE |
 * | 4 | odd part of <R> | null: R is even in Phi_AEE |
 *
 * A residue in panels 2-4 is not automatically a bug: the same magnetic field that produces panel 1
 * breaks the mirror these nulls rest on. Read them together, not in isolation.
 *
 * @param taskDir    Top-level task TDirectory.
 * @param outDir     Output sub-directory for the canvas and the folded histograms.
 * @param f          Folder to process (see kFolders).
 * @param minEntries Minimum per-bin entry count required of BOTH partners of a folded pair.
 */
// ==========================================================================
static void MakeAeeFoldCanvas(TDirectory* taskDir, TDirectory* outDir, const FolderSpec& f, double minEntries = 50.)
{
    TDirectory* folderDir = GetDir(taskDir, f.name);
    if (!folderDir) { printf("WARNING: skipping AEE fold for '%s' (folder missing)\n", f.name); return; }

    TDirectory* aeeDir = GetDirPath(folderDir, "PolMaps/Aee");
    if (!aeeDir) { printf("WARNING: skipping AEE fold for '%s' (PolMaps/Aee missing)\n", f.name); return; }

    TH2*        hCounts = static_cast<TH2*>(SafeGet(aeeDir, "h2dCountsVsPxAeePyAee", false));
    TProfile2D* hPt     = static_cast<TProfile2D*>(SafeGet(aeeDir, "p2dPtStar_vsPxAeePyAee", false));
    TProfile2D* hPz     = static_cast<TProfile2D*>(SafeGet(aeeDir, "p2dPzStar_vsPxAeePyAee", false));
    TProfile2D* hRing   = static_cast<TProfile2D*>(SafeGet(folderDir, "p2dRingObservableVsPxAeePyAee", false));

    if (!hCounts || !hPt || !hPz || !hRing) {
        printf("WARNING: skipping AEE fold for '%s' (one or more (x_AEE, y_AEE) maps missing)\n", f.name);
        return;
    }

    if (!HasFoldableYAxis(hCounts)) {
        printf("WARNING: skipping AEE fold for '%s': the y_AEE axis is not symmetric about zero with an "
               "even bin count, so the fold cannot pair bins. Check axisLambdaPRot.\n", f.name);
        return;
    }

    TH2D* foldCounts = MakeCountsAsymmetryY(hCounts, Form("hAeeFoldCounts_%s", f.name));
    TH2D* foldPt     = MakeProfileFoldY(hPt,   true,  minEntries, Form("hAeeFoldPtOdd_%s", f.name));
    TH2D* foldPz     = MakeProfileFoldY(hPz,   false, minEntries, Form("hAeeFoldPzEven_%s", f.name));
    TH2D* foldRing   = MakeProfileFoldY(hRing, true,  minEntries, Form("hAeeFoldRingOdd_%s", f.name));

    if (!foldCounts || !foldPt || !foldPz || !foldRing) {
        printf("WARNING: skipping AEE fold for '%s' (fold construction failed)\n", f.name);
        return;
    }

    const std::string axPair = ";" + kPxA + ";" + kPyA + ";";

    TCanvas* c = new TCanvas(Form("cAeeFold_%s", f.name), "", 2800, 650);
    c->Divide(4, 1, 0.005, 0.005);

    c->cd(1);
    SetupPadMargins();
    // Already a ratio, so it is shown as a percentage of the pair total rather than rescaled again.
    DrawSymmetricColzH2(foldCounts, (" " + axPair + "counts asymmetry [%]").c_str(), true);
    AddLabel(0.5, 0.95, "AEE: odd part of the counts", 0.045, 22);

    c->cd(2);
    SetupPadMargins();
    DrawSymmetricColzH2(foldPt, (" " + axPair + "odd part of <P*_{T}> [%]").c_str(), true);
    AddLabel(0.5, 0.95, "Null: <P*_{T}> should be even", 0.045, 22);

    c->cd(3);
    SetupPadMargins();
    DrawSymmetricColzH2(foldPz, (" " + axPair + "even part of <P*_{z}> [%]").c_str(), true);
    AddLabel(0.5, 0.95, "Null: <P*_{z}> should be odd", 0.045, 22);

    c->cd(4);
    SetupPadMargins();
    DrawSymmetricColzH2(foldRing, (" " + axPair + RingLabel("odd part of <#it{R}> [%]")).c_str(), true);
    AddLabel(0.5, 0.95, RingLabel("Null: <#it{R}> should be even").c_str(), 0.045, 22);

    c->cd(0);
    AddLabel(0.5, 0.98, Form("AEE fold about #Phi_{AEE} = 0  --  %s", f.label), 0.032, 22);

    WriteCanvas(c, outDir);

    // The folded maps themselves, so the asymmetry can be projected or integrated later:
    outDir->cd();
    foldCounts->Write();
    foldPt->Write();
    foldPz->Write();
    foldRing->Write();

    delete c;
}


// ==========================================================================
// Section 6: the ring projection kernel.
//
// Writing the acceptance moment as \vec B = <\hat p*>_epsilon and decomposing it in the Lambda's own
// triad (\hat p_Lambda, \hat phi, \hat theta), the ring projection closes in elementary form. At fixed
// DeltaTheta_jet there is no azimuthal freedom left -- cos(dphi) is determined -- and the result is
//
//     <R> = (3/alpha) B_phi (t_z - lambda_z cos dTheta) / (sin(theta_Lambda) sin dTheta)
//
// with t_z = cos(theta_proxy) = tanh(eta_proxy) and lambda_z = cos(theta_Lambda). Averaging over the
// accepted Lambdas at fixed (cos dTheta, t_z) leaves only two numbers:
//
//     <R>(c, t_z) = [ M0 t_z - M1 c ] / sqrt(1 - c^2)  +  S(c),     c = cos dTheta
//     M0 = <B_phi cosh(eta_Lambda)>,   M1 = <B_phi sinh(eta_Lambda)>
//
// (M0 and M1 absorb the 3/alpha, which the consumer already applied to R.)
//
// The kernel is EXACTLY affine in t_z only at fixed (c, lambda_z). This surface integrates over
// lambda_z, and at fixed (c, t_z) a Lambda enters with phase-space weight 1/sqrt(G), G the Gram
// determinant of (beam, reference axis, Lambda), which depends on t_z. M0 and M1 are therefore
// cell-local, functions of (c, t_z) themselves, and the slope in t_z at fixed c is
//
//     slope * sqrt(1 - c^2) = M0 + t_z dM0/dt_z - c dM1/dt_z   (+ an east-west branch term)
//
// rather than M0 alone. The geometric part of M1 is odd in both t_z and c, so -c dM1/dt_z is EVEN in
// c. The predictions are therefore (derivations in the ring-geometry note):
//
//   * slope * sqrt(1 - c^2) is a BOWL, even in c, whose minimum at c = 0 is the cell-local M0 there.
//     Its rise is geometry, present for a constant polarization and a perfect detector. Outside the
//     fully allowed window (dashed lines, see FullyAllowedCosWindow()) kinematic truncation steepens
//     it further and makes it binning-dependent;
//   * the jet and leading-particle bowls need NOT coincide: a bowl depends on the proxy's own t_z
//     range and on its correlation with the Lambda. Proxy independence is exact only at fixed
//     lambda_z, i.e. in Section 7;
//   * the fit range is symmetric in t_z, so the intercept is the value at t_z = 0, where the geometric
//     M1 vanishes: intercept * sqrt(1 - c^2) = -M1(0, c) c + S(c) sqrt(1 - c^2), odd in c whenever the
//     signal (and the branch term) vanish.
//
// Deliberately NOT done here: separating M1 from S by parity in c. That would need S to be even about
// c = 0, and it is not -- the recoil jet sits near pi - dTheta. Section 8 does that separation at
// fixed lambda_z, where it is exact.
// ==========================================================================

/// @brief One (cos dTheta) slice: the straight-line fit of <R> against the proxy direction cosine.
struct KernelSliceFit {
    double cosCenter = 0.0;   ///< Bin-group centre in cos(dTheta)
    double sinTheta  = 0.0;   ///< sqrt(1 - cosCenter^2), the kernel's denominator
    double entries   = 0.0;
    double slope = 0.0, slopeErr = 0.0;         ///< d<R>/dt_z at fixed cos(dTheta)
    double intercept = 0.0, interceptErr = 0.0; ///< <R> at t_z = 0
    double chi2 = 0.0;
    int    ndf = 0;
    bool   valid = false;
};

// ==========================================================================
/**
 * @brief Fits <R> against the proxy direction cosine in slices of cos(dTheta).
 *
 * Fit configuration, and how it differs from the tanh fits elsewhere in this analysis:
 *  - A straight line, not a curve, and that is the point. The kernel is exactly linear in t_z at fixed
 *    cos(dTheta), so there is no model choice to make and no degeneracy to canonicalise.
 *  - Plain chi2 ("QS"), no Minos. The tanh fits need Minos because p0 and p1 collapse onto a curved
 *    valley at small p1*eta; a straight line has no such valley, and with a t_z axis symmetric about
 *    zero the slope and intercept are very nearly uncorrelated. Parabolic errors are adequate here,
 *    and asking for Minos would only cost time.
 *  - Bins are means carrying SEM errors, so chi2 is the correct method. A likelihood option would be
 *    wrong for the same reason it is wrong for the tanh fits: these are not Poisson counts.
 *
 * @param h2         Source TProfile2D, x = cos(dTheta), y = proxy direction cosine.
 * @param cosRebin   Number of x bins merged per slice. The per-slice fit needs the statistics.
 * @param maxAbsCos  Slices centred beyond this |cos| are skipped: the kernel diverges as
 *                   1/sqrt(1 - c^2) there, and the acceptance has already emptied those bins.
 * @param minEntries Minimum entries in a slice for it to be fitted at all.
 * @param fitDir     Directory for the per-slice canvases; pass nullptr to skip drawing them.
 * @param tag        Name prefix for those canvases.
 * @return One entry per attempted slice, in increasing cos(dTheta). Failed fits carry valid = false.
 */
// ==========================================================================
static std::vector<KernelSliceFit> FitRingKernelSlices(TProfile2D* h2, int cosRebin, double maxAbsCos,
                                                       double minEntries, TDirectory* fitDir,
                                                       const std::string& tag)
{
    std::vector<KernelSliceFit> out;
    if (!h2 || cosRebin < 1) return out;

    const int nx = h2->GetNbinsX();
    const double yLo = h2->GetYaxis()->GetXmin();
    const double yHi = h2->GetYaxis()->GetXmax();

    for (int ix = 1; ix + cosRebin - 1 <= nx; ix += cosRebin) {
        const int ixHi = ix + cosRebin - 1;
        const double cLo = h2->GetXaxis()->GetBinLowEdge(ix);
        const double cHi = h2->GetXaxis()->GetBinUpEdge(ixHi);

        KernelSliceFit k;
        k.cosCenter = 0.5 * (cLo + cHi);
        if (std::fabs(k.cosCenter) > maxAbsCos) continue;
        k.sinTheta = std::sqrt(std::max(0.0, 1.0 - k.cosCenter * k.cosCenter));

        // ProfileY collapses the x range into one profile along y, which is exactly the slice wanted.
        TProfile* slice = h2->ProfileY(Form("%s_slice_%d_%d", tag.c_str(), ix, gCloneIdx++), ix, ixHi);
        if (!slice) continue;
        slice->SetDirectory(nullptr);
        k.entries = slice->GetEntries();

        if (k.entries < minEntries) { out.push_back(k); delete slice; continue; }

        TF1* fLin = new TF1(Form("%s_fn_%d", tag.c_str(), gCloneIdx++), "[p0]+[p1]*x", yLo, yHi);
        fLin->SetParameters(0.0, 0.0);
        fLin->SetLineWidth(2);
        fLin->SetNpx(300);

        TFitResultPtr res = slice->Fit(fLin, "QS", "", yLo, yHi);
        if (res.Get()) {
            k.chi2 = res->Chi2();
            k.ndf  = res->Ndf();
            k.intercept    = res->Parameter(0);
            k.interceptErr = res->ParError(0);
            k.slope        = res->Parameter(1);
            k.slopeErr     = res->ParError(1);
            k.valid = res->IsValid() && (k.ndf > 0);
        }

        if (fitDir && k.valid) {
            // Style set before the first paint: the stats box is built at paint time and reads gStyle
            // then, so setting it afterwards has no effect.
            gStyle->SetOptStat(10);
            gStyle->SetOptFit(111);

            TCanvas* c = new TCanvas(Form("%s_cos%+.2f", tag.c_str(), k.cosCenter), "", 800, 600);
            c->SetLeftMargin(0.13);
            c->SetBottomMargin(0.12);
            c->SetGridx();
            c->SetGridy();

            slice->SetStats(1);
            slice->SetMarkerStyle(20);
            slice->SetLineWidth(2);
            slice->GetXaxis()->SetTitle("#hat{t}_{z}");
            slice->GetYaxis()->SetTitle(RingLabel("<#it{R}>").c_str());
            slice->SetTitle(Form("%.2f < cos#Delta#theta < %.2f", cLo, cHi));
            slice->Draw("PE");

            // Forcing a paint pass so the TPaveStats is actually built: TCanvas::Write() does not
            // paint, and without this the canvas is stored with no stats object at all.
            c->Update();
            TPaveStats* st = dynamic_cast<TPaveStats*>(slice->FindObject("stats"));
            if (st) {
                st->SetOptStat(10);
                st->SetOptFit(111);
                st->SetX1NDC(0.15); st->SetX2NDC(0.47);
                st->SetY1NDC(0.72); st->SetY2NDC(0.90);
                st->SetFillColor(kWhite);
                st->SetFillStyle(1001);
                st->SetBorderSize(1);
                st->SetTextSize(0.030);
            }
            c->Modified();
            c->Update();

            fitDir->cd();
            c->Write();
            delete c;
        }

        delete fLin;
        delete slice;
        out.push_back(k);
    }
    return out;
}

// ==========================================================================
/**
 * @brief Builds the kernel graph from a set of slice fits.
 *
 * @param fits    Slice fits, as returned by FitRingKernelSlices().
 * @param useSlope true -> slope * sin(dTheta) (the M0 estimator);
 *                 false -> intercept * sin(dTheta) (the M1 + signal combination).
 * @param name    Object name for the graph.
 * @return A graph owned by the caller, or nullptr if no slice was usable.
 */
// ==========================================================================
static TGraphErrors* MakeKernelGraph(const std::vector<KernelSliceFit>& fits, bool useSlope, const char* name)
{
    std::vector<double> x, y, ex, ey;
    for (const KernelSliceFit& k : fits) {
        if (!k.valid) continue;
        const double v = useSlope ? k.slope    : k.intercept;
        const double e = useSlope ? k.slopeErr : k.interceptErr;
        x.push_back(k.cosCenter);
        y.push_back(v * k.sinTheta);
        ex.push_back(0.0);
        ey.push_back(e * k.sinTheta);
    }
    if (x.empty()) return nullptr;

    TGraphErrors* g = new TGraphErrors(static_cast<int>(x.size()), x.data(), y.data(), ex.data(), ey.data());
    g->SetName(name);
    g->SetMarkerStyle(20);
    g->SetLineWidth(2);
    return g;
}

// ==========================================================================
/**
 * @brief Largest |bin edge| along one axis of a TProfile3D, among bins holding a non-negligible share
 * of the entries: the acceptance edge along that axis, read off the data rather than repeated here.
 *
 * Deliberately conservative: it returns the OUTER edge of the outermost populated bin, so the true
 * acceptance edge is at or inside it, and any window built from it is at or inside the true window.
 *
 * @param h3           Source profile.
 * @param axisIndex    0 = x, 1 = y, 2 = z.
 * @param relThreshold A bin counts as populated above this fraction of the most populated one.
 * @return The extent, or 0 if the profile is empty.
 */
// ==========================================================================
static double OccupiedExtent(TProfile3D* h3, int axisIndex, double relThreshold = 1.e-3)
{
    if (!h3) return 0.;
    TAxis* ax = (axisIndex == 0) ? h3->GetXaxis() : (axisIndex == 1) ? h3->GetYaxis() : h3->GetZaxis();
    const int nx = h3->GetNbinsX(), ny = h3->GetNbinsY(), nz = h3->GetNbinsZ();
    std::vector<double> proj(ax->GetNbins() + 2, 0.);

    for (int ix = 1; ix <= nx; ++ix)
        for (int iy = 1; iy <= ny; ++iy)
            for (int iz = 1; iz <= nz; ++iz) {
                const int b = (axisIndex == 0) ? ix : (axisIndex == 1) ? iy : iz;
                proj[b] += h3->GetBinEntries(h3->GetBin(ix, iy, iz));
            }

    const double peak = *std::max_element(proj.begin(), proj.end());
    if (peak <= 0.) return 0.;
    double extent = 0.;
    for (int b = 1; b <= ax->GetNbins(); ++b)
        if (proj[b] > relThreshold * peak)
            extent = std::max(extent, std::max(std::fabs(ax->GetBinLowEdge(b)), std::fabs(ax->GetBinUpEdge(b))));
    return extent;
}

// ==========================================================================
/**
 * @brief Half-width in cos(dTheta) of the fully allowed window: the range of opening angles inside
 * which every (reference axis, Lambda) pair of the acceptance is kinematically allowed.
 *
 * From the spherical triangle inequality, with vartheta the smallest polar angle in each acceptance
 * (vartheta = acos of the largest direction cosine), the window is
 *     pi - vartheta_t - vartheta_L <= dTheta <= vartheta_t + vartheta_L,
 * i.e. |cos dTheta| <= -cos(vartheta_t + vartheta_L), and it is empty when vartheta_t + vartheta_L < pi/2.
 * Both acceptances are read off the 3D kernel profile: its y axis is t_z, its z axis lambda_z.
 *
 * @return The half-width, or -1 when the window is empty or the profile is missing.
 */
// ==========================================================================
static double FullyAllowedCosWindow(TProfile3D* h3)
{
    if (!h3) return -1.;
    const double tMax = std::min(OccupiedExtent(h3, 1), 1.0);
    const double lMax = std::min(OccupiedExtent(h3, 2), 1.0);
    if (tMax <= 0. || lMax <= 0.) return -1.;
    const double sum = std::acos(tMax) + std::acos(lMax);
    const double halfPi = 0.5 * std::acos(-1.0);
    return (sum < halfPi) ? -1. : -std::cos(sum);
}

/// @brief Dashed vertical lines at cos(dTheta) = +-halfWidth across the current pad (no-op if <= 0).
static void DrawWindowLines(double halfWidth, int color)
{
    if (halfWidth <= 0.) return;
    gPad->Update(); // The user y range is only valid once the pad has been laid out
    for (double x : {-halfWidth, halfWidth}) {
        TLine* l = new TLine(x, gPad->GetUymin(), x, gPad->GetUymax());
        l->SetLineStyle(2);
        l->SetLineColor(color);
        l->SetLineWidth(2);
        l->Draw();
    }
}

// ==========================================================================
/**
 * @brief One canvas per folder: the kernel surfaces, the M0 estimator, and the intercept combination.
 *
 * | Panel | Content | Expectation |
 * |---|---|---|
 * | 1 | <R> over (cos dTheta, t_z), leading jet | a tilted surface, steeper towards \|cos\| -> 1 |
 * | 2 | slope * sin(dTheta), jet and LeadP overlaid | a bowl, even in cos, minimum near M0 at cos = 0 |
 * | 3 | intercept * sin(dTheta) | -M1(0,c) c + S sin(dTheta); odd in cos if S = 0 |
 *
 * The dashed lines mark each proxy's fully allowed window (FullyAllowedCosWindow()).
 *
 * Panel 2 is NOT a flatness test: the bowl is geometry (see the Section 6 header). What is
 * informative is its symmetry about cos = 0 -- an odd component points at a near-side correlation or
 * at an efficiency that knows where the reference axis is -- and its minimum, which is close to, but
 * not identical with, Section 7's sample-wide M0. The jet and LeadP bowls need not coincide.
 *
 * @param taskDir    Top-level task TDirectory.
 * @param outDir     Output sub-directory for the canvas and the graphs.
 * @param f          Folder to process (see kFolders).
 * @param cosRebin   cos(dTheta) bins merged per slice.
 * @param maxAbsCos  |cos(dTheta)| ceiling for a slice to be fitted.
 * @param minEntries Minimum entries per slice.
 * @param drawSlices Whether to write the individual slice-fit canvases into a Fits/ subdirectory.
 */
// ==========================================================================
static void MakeRingKernelCanvas(TDirectory* taskDir, TDirectory* outDir, const FolderSpec& f,
                                 int cosRebin = 5, double maxAbsCos = 0.95,
                                 double minEntries = 200., bool drawSlices = true)
{
    TDirectory* folderDir = GetDir(taskDir, f.name);
    if (!folderDir) { printf("WARNING: skipping ring kernel for '%s' (folder missing)\n", f.name); return; }

    TDirectory* kDir = GetDirPath(folderDir, "RingKernel");
    if (!kDir) { printf("WARNING: skipping ring kernel for '%s' (RingKernel missing)\n", f.name); return; }

    TProfile2D* hJet = static_cast<TProfile2D*>(SafeGet(kDir, "p2dRingObservableCosDeltaThetaVsJetZ", false));
    TProfile2D* hLeadP = static_cast<TProfile2D*>(SafeGet(kDir, "p2dRingObservableLeadPCosDeltaThetaVsLeadPZ", false));
    if (!hJet) { printf("WARNING: skipping ring kernel for '%s' (jet surface missing)\n", f.name); return; }

    TDirectory* fitDir = nullptr;
    if (drawSlices) {
        fitDir = outDir->mkdir(Form("Fits_%s", f.name));
    }

    // Fully allowed windows, from each proxy's own 3D occupancy (absent profiles just skip the lines)
    const double cwJet = FullyAllowedCosWindow(static_cast<TProfile3D*>(
        SafeGet(kDir, "p3dRingObservableCosDeltaThetaVsJetZVsLambdaZ", false)));
    const double cwLeadP = FullyAllowedCosWindow(static_cast<TProfile3D*>(
        SafeGet(kDir, "p3dRingObservableLeadPCosDeltaThetaVsLeadPZVsLambdaZ", false)));

    const std::vector<KernelSliceFit> fitJet =
        FitRingKernelSlices(hJet, cosRebin, maxAbsCos, minEntries, fitDir, Form("kJet_%s", f.name));
    const std::vector<KernelSliceFit> fitLeadP =
        hLeadP ? FitRingKernelSlices(hLeadP, cosRebin, maxAbsCos, minEntries, fitDir, Form("kLeadP_%s", f.name))
               : std::vector<KernelSliceFit>{};

    TGraphErrors* gSlopeJet   = MakeKernelGraph(fitJet,   true,  Form("gM0Jet_%s", f.name));
    TGraphErrors* gSlopeLeadP = MakeKernelGraph(fitLeadP, true,  Form("gM0LeadP_%s", f.name));
    TGraphErrors* gIntJet     = MakeKernelGraph(fitJet,   false, Form("gInterceptJet_%s", f.name));
    TGraphErrors* gIntLeadP   = MakeKernelGraph(fitLeadP, false, Form("gInterceptLeadP_%s", f.name));

    if (!gSlopeJet) {
        printf("WARNING: skipping ring kernel for '%s' (no slice reached %.0f entries)\n", f.name, minEntries);
        return;
    }

    TCanvas* c = new TCanvas(Form("cRingKernel_%s", f.name), "", 2100, 650);
    c->Divide(3, 1, 0.005, 0.005);

    c->cd(1);
    SetupPadMargins();
    hJet->GetZaxis()->SetTitleOffset(1.5);
    DrawSymmetricColz(hJet, RingLabel(" ;cos#Delta#theta_{jet};#hat{t}_{z};<#it{R}> [%]").c_str());
    AddLabel(0.5, 0.95, "Kernel surface, leading jet", 0.045, 22);

    // Panels 2 and 3 share their styling; only the quantity and the expectation differ.
    for (int panel = 2; panel <= 3; ++panel) {
        const bool isSlope = (panel == 2);
        TGraphErrors* gJ = isSlope ? gSlopeJet   : gIntJet;
        TGraphErrors* gL = isSlope ? gSlopeLeadP : gIntLeadP;
        if (!gJ) continue;

        c->cd(panel);
        SetupPadMargins();
        gPad->SetGridx();
        gPad->SetGridy();

        gJ->SetMarkerColor(kAzure + 2);
        gJ->SetLineColor(kAzure + 2);
        gJ->SetTitle(isSlope ? " ;cos#Delta#theta;slope #times sin#Delta#theta"
                             : " ;cos#Delta#theta;intercept #times sin#Delta#theta");
        gJ->Draw("AP");
        if (gL) {
            gL->SetMarkerColor(kRed + 1);
            gL->SetLineColor(kRed + 1);
            gL->SetMarkerStyle(21);
            gL->Draw("P SAME");
        }

        DrawWindowLines(cwJet, kAzure + 2);
        if (gL) DrawWindowLines(cwLeadP, kRed + 1);

        TLegend* leg = new TLegend(0.58, 0.76, 0.84, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(gJ, "Leading jet", "lp");
        if (gL) leg->AddEntry(gL, "Leading particle", "lp");
        leg->Draw();

        AddLabel(0.5, 0.95, isSlope ? "Geometric bowl, minimum #approx M_{0} at cos#Delta#theta = 0"
                                    : "Odd in cos#Delta#theta unless S #neq 0", 0.045, 22);
    }

    c->cd(0);
    AddLabel(0.5, 0.98, Form("%sRing projection kernel  --  %s", gNullTag.c_str(), f.label), 0.032, 22);

    WriteCanvas(c, outDir);

    // The graphs themselves, so M0 can be averaged, compared across configs, or propagated:
    outDir->cd();
    gSlopeJet->Write();
    if (gSlopeLeadP) gSlopeLeadP->Write();
    if (gIntJet) gIntJet->Write();
    if (gIntLeadP) gIntLeadP->Write();

    delete c;
}


// ==========================================================================
// Section 7: B_phi, and the M0 / M1 moments.
//
// Section 6's bowl brackets M0 without pinning it down. This one goes behind it, recovering the
// axis-independent azimuthal component differentially in cos(theta_Lambda), where the kernel is
// exactly affine in t_z, and building both moments from it.
//
// Naming: the code calls that component B_phi. It is the WHOLE axis-independent azimuthal component --
// detector-induced plus any genuine polarization normal to the production plane (P_N) -- and the two
// are not separated here. Its moments are sample-wide (whole-slab weights), so the geometric, t_z-odd
// part of Section 6's cell-local M1 cancels in them. What remains in M1 is the eta_Lambda-odd content
// of the sample: an A/C yield asymmetry, an A/C acceptance asymmetry, or a genuine P_N.
//
// Predicting it instead is straightforward. Writing lambda_z = cos(theta_Lambda) and
// lambda_T = sin(theta_Lambda) = sqrt(1 - lambda_z^2), the closed form reads
//
//     <R>(c, t_z, lambda_z) = B_phi(lambda_z) (t_z - lambda_z c) / (lambda_T sqrt(1 - c^2))
//
// so at fixed (lambda_z, c) the slope in t_z is m = K / sqrt(1 - c^2), with K = B_phi / lambda_T.
// That gives the extraction, and the identities that follow it:
//
//     m * sqrt(1 - c^2) = K,   CONSTANT in c            <-- the model test, per lambda_z slice
//     cosh(eta_Lambda) = 1 / lambda_T                   <-- so M0 = <K>_w
//     sinh(eta_Lambda) = lambda_z / lambda_T            <-- so M1 = <lambda_z K>_w
//
// with w(lambda_z) the Lambda occupancy, read off the profile's own per-bin entry counts. Nothing
// else is needed: the 3D profile alone carries the moments, the weights and the model test.
//
// Every step here inherits Section 6's immunity to the signal, because K comes from a slope in t_z
// and the signal is flat in t_z. Verified on a synthetic profile carrying a known B_phi, a 35% A/C
// yield asymmetry and a non-trivial S(c): B_phi recovered to 0.1%, M0 and M1 to better than 1%.
//
// Symmetrising w(lambda_z) removes the yield asymmetry, and nothing else. The residual M1 after it is
// non-zero when B_phi itself has a part odd in eta_Lambda: an A/C ACCEPTANCE asymmetry, or a genuine
// P_N. Field reversal tells them apart -- a detector-induced B_phi is odd in the field, P_N is not.
// ==========================================================================

/// @brief Weighted straight-line fit, done in closed form rather than through Minuit.
struct LinFit {
    double slope = 0.0, slopeErr = 0.0;
    double intercept = 0.0, interceptErr = 0.0;
    int    n = 0;
    bool   valid = false;
};

// ==========================================================================
/**
 * @brief Inverse-variance weighted least squares for y = intercept + slope * x.
 *
 * Done by hand on purpose. This is called once per (lambda_z, cos dTheta) cell -- hundreds of times
 * per folder -- and a two-parameter linear fit has a closed-form solution, so routing each one
 * through TF1 and Minuit would buy nothing but overhead and a pile of transient ROOT objects.
 * Section 6's fits go through ROOT because their canvases are meant to be looked at; these are not.
 *
 * @param x,y,e Points and their errors. Entries with e <= 0 must already be filtered out.
 * @return The fit; valid = false when fewer than three points or a degenerate normal matrix.
 */
// ==========================================================================
static LinFit WeightedLineFit(const std::vector<double>& x, const std::vector<double>& y,
                              const std::vector<double>& e)
{
    LinFit out;
    out.n = static_cast<int>(x.size());
    if (out.n < 3) return out;

    double s = 0., sx = 0., sy = 0., sxx = 0., sxy = 0.;
    for (int i = 0; i < out.n; ++i) {
        const double w = 1.0 / (e[i] * e[i]);
        s += w; sx += w * x[i]; sy += w * y[i];
        sxx += w * x[i] * x[i]; sxy += w * x[i] * y[i];
    }
    const double det = s * sxx - sx * sx;
    if (std::fabs(det) < 1.e-300) return out;

    out.slope        = (s * sxy - sx * sy) / det;
    out.slopeErr     = std::sqrt(s / det);
    out.intercept    = (sxx * sy - sx * sxy) / det;
    out.interceptErr = std::sqrt(sxx / det);
    out.valid = true;
    return out;
}

/// @brief B_phi and its supporting numbers for one cos(theta_Lambda) slab.
struct KernelBphiPoint {
    double lambdaZ = 0.0;   ///< Slab centre, cos(theta_Lambda) = tanh(eta_Lambda)
    double lambdaT = 0.0;   ///< sin(theta_Lambda) = sqrt(1 - lambdaZ^2)
    double weight  = 0.0;   ///< Lambda occupancy of the slab, i.e. w(lambda_z)
    double kappa = 0.0, kappaErr = 0.0; ///< K = B_phi / lambdaT = B_phi cosh(eta_Lambda)
    double bPhi  = 0.0, bPhiErr  = 0.0; ///< K * lambdaT
    double spreadChi2 = 0.0;            ///< Scatter of m*sqrt(1-c^2) about K, across cos(dTheta)
    int    spreadNdf = 0;               ///< Its ndf. A large chi2/ndf means the closed form is failing
    bool   valid = false;
};

// ==========================================================================
/**
 * @brief Recovers B_phi per cos(theta_Lambda) slab from a 3D kernel profile.
 *
 * For each slab, every cos(dTheta) column supplies an independent estimate of the same K through the
 * slope of <R> against t_z. They are combined with inverse-variance weights, and their scatter about
 * the combination is kept as spreadChi2: that scatter is the model test, since the closed form says
 * m * sqrt(1 - c^2) must not depend on c at all.
 *
 * @param h3              Kernel TProfile3D: x = cos(dTheta), y = proxy direction cosine, z = lambda_z.
 * @param maxAbsCos       Columns centred beyond this |cos| are skipped (the kernel diverges there).
 * @param minEntriesCell  Minimum entries for a single (x,y,z) cell to enter a column fit.
 * @return One entry per slab, in increasing lambda_z. Unusable slabs carry valid = false.
 */
// ==========================================================================
static std::vector<KernelBphiPoint> ExtractBphiFromKernel3D(TProfile3D* h3, double maxAbsCos,
                                                             double minEntriesCell)
{
    std::vector<KernelBphiPoint> out;
    if (!h3) return out;

    const int nc = h3->GetNbinsX();
    const int nt = h3->GetNbinsY();
    const int nl = h3->GetNbinsZ();

    for (int iz = 1; iz <= nl; ++iz) {
        KernelBphiPoint p;
        p.lambdaZ = h3->GetZaxis()->GetBinCenter(iz);
        p.lambdaT = std::sqrt(std::max(0.0, 1.0 - p.lambdaZ * p.lambdaZ));
        if (p.lambdaT < 1.e-6) { out.push_back(p); continue; }

        // Combination accumulators over the cos(dTheta) columns of this slab:
        double sumW = 0., sumWV = 0.;
        std::vector<double> colVal, colErr;

        for (int ix = 1; ix <= nc; ++ix) {
            const double c = h3->GetXaxis()->GetBinCenter(ix);
            if (std::fabs(c) > maxAbsCos) continue;
            const double sinC = std::sqrt(std::max(0.0, 1.0 - c * c));
            if (sinC < 1.e-6) continue;

            std::vector<double> x, y, e;
            for (int iy = 1; iy <= nt; ++iy) {
                const int gb = h3->GetBin(ix, iy, iz);
                if (h3->GetBinEntries(gb) < minEntriesCell) continue;
                const double err = h3->GetBinError(ix, iy, iz);
                if (err <= 0.) continue;
                x.push_back(h3->GetYaxis()->GetBinCenter(iy));
                y.push_back(h3->GetBinContent(ix, iy, iz));
                e.push_back(err);
                p.weight += h3->GetBinEntries(gb);
            }

            const LinFit lf = WeightedLineFit(x, y, e);
            if (!lf.valid) continue;

            // m * sqrt(1 - c^2) estimates K, and must do so independently of c:
            const double v = lf.slope * sinC;
            const double ev = lf.slopeErr * sinC;
            if (ev <= 0.) continue;

            const double w = 1.0 / (ev * ev);
            sumW += w; sumWV += w * v;
            colVal.push_back(v); colErr.push_back(ev);
        }

        if (sumW <= 0. || colVal.size() < 2) { out.push_back(p); continue; }

        p.kappa    = sumWV / sumW;
        p.kappaErr = 1.0 / std::sqrt(sumW);
        p.bPhi     = p.kappa * p.lambdaT;
        p.bPhiErr  = p.kappaErr * p.lambdaT;

        for (size_t i = 0; i < colVal.size(); ++i) {
            const double d = (colVal[i] - p.kappa) / colErr[i];
            p.spreadChi2 += d * d;
        }
        p.spreadNdf = static_cast<int>(colVal.size()) - 1;
        p.valid = true;
        out.push_back(p);
    }
    return out;
}

/// @brief The two moments of the kernel, as built from a set of KernelBphiPoint.
struct KernelMoments {
    double m0 = 0.0, m0Err = 0.0;       ///< <B_phi cosh(eta_Lambda)>, the field-driven term
    double m1 = 0.0, m1Err = 0.0;       ///< <B_phi sinh(eta_Lambda)>, the A/C yield-asymmetry term
    double m0Sym = 0.0, m0SymErr = 0.0; ///< Same, with the lambda_z weights symmetrised
    double m1Sym = 0.0, m1SymErr = 0.0;
    double yieldAsym = 0.0;             ///< (N(+) - N(-)) / (N(+) + N(-)) over the whole slab set
    bool   symmetrisable = false;       ///< False when the lambda_z axis cannot be mirror-paired
    bool   valid = false;
};

// ==========================================================================
/**
 * @brief Builds M0 and M1 from the per-slab B_phi, both as measured and eta-symmetrised.
 *
 * Symmetrising replaces w(lambda_z) by min(w(lambda_z), w(-lambda_z)), the largest weight set that is
 * mirror-symmetric and never up-weights a slab beyond what was actually recorded. It requires the
 * lambda_z axis to pair bin iz with bin N+1-iz, which axisLambdaZ does -- but it is a ConfigurableAxis,
 * so the pairing is verified rather than assumed, exactly as the AEE fold verifies its own.
 *
 * @param pts Slab results from ExtractBphiFromKernel3D(), in increasing lambda_z.
 */
// ==========================================================================
static KernelMoments ComputeKernelMoments(const std::vector<KernelBphiPoint>& pts)
{
    KernelMoments out;
    const int n = static_cast<int>(pts.size());
    if (n < 2) return out;

    // The axis pairs iz with n+1-iz only if it is symmetric about zero with an even bin count:
    out.symmetrisable = (n % 2 == 0);
    if (out.symmetrisable) {
        for (int i = 0; i < n && out.symmetrisable; ++i) {
            const double span = std::fabs(pts[n - 1].lambdaZ - pts[0].lambdaZ);
            if (std::fabs(pts[i].lambdaZ + pts[n - 1 - i].lambdaZ) > 1.e-6 * std::max(span, 1.e-6))
                out.symmetrisable = false;
        }
    }

    double sw = 0., swSym = 0., nPos = 0., nNeg = 0.;
    for (int i = 0; i < n; ++i) {
        if (!pts[i].valid) continue;
        sw += pts[i].weight;
        (pts[i].lambdaZ >= 0. ? nPos : nNeg) += pts[i].weight;
        if (out.symmetrisable && pts[n - 1 - i].valid)
            swSym += std::min(pts[i].weight, pts[n - 1 - i].weight);
    }
    if (sw <= 0.) return out;

    if (nPos + nNeg > 0.) out.yieldAsym = (nPos - nNeg) / (nPos + nNeg);

    double v0 = 0., v1 = 0., e0 = 0., e1 = 0.;
    double v0s = 0., v1s = 0., e0s = 0., e1s = 0.;
    for (int i = 0; i < n; ++i) {
        const KernelBphiPoint& p = pts[i];
        if (!p.valid) continue;

        const double w = p.weight / sw;
        v0 += w * p.kappa;
        v1 += w * p.lambdaZ * p.kappa;
        e0 += (w * p.kappaErr) * (w * p.kappaErr);
        e1 += (w * p.lambdaZ * p.kappaErr) * (w * p.lambdaZ * p.kappaErr);

        if (out.symmetrisable && swSym > 0. && pts[n - 1 - i].valid) {
            const double ws = std::min(p.weight, pts[n - 1 - i].weight) / swSym;
            v0s += ws * p.kappa;
            v1s += ws * p.lambdaZ * p.kappa;
            e0s += (ws * p.kappaErr) * (ws * p.kappaErr);
            e1s += (ws * p.lambdaZ * p.kappaErr) * (ws * p.lambdaZ * p.kappaErr);
        }
    }

    out.m0 = v0; out.m0Err = std::sqrt(e0);
    out.m1 = v1; out.m1Err = std::sqrt(e1);
    out.m0Sym = v0s; out.m0SymErr = std::sqrt(e0s);
    out.m1Sym = v1s; out.m1SymErr = std::sqrt(e1s);
    out.valid = true;
    return out;
}

/// @brief Graph of B_phi against lambda_z, optionally mirrored so the evenness test is a visual overlay.
static TGraphErrors* MakeBphiGraph(const std::vector<KernelBphiPoint>& pts, bool mirror, const char* name)
{
    std::vector<double> x, y, ex, ey;
    for (const KernelBphiPoint& p : pts) {
        if (!p.valid) continue;
        x.push_back(mirror ? -p.lambdaZ : p.lambdaZ);
        y.push_back(p.bPhi * 100.0); // percent, matching every other polarization plot in this file
        ex.push_back(0.0);
        ey.push_back(p.bPhiErr * 100.0);
    }
    if (x.empty()) return nullptr;

    TGraphErrors* g = new TGraphErrors(static_cast<int>(x.size()), x.data(), y.data(), ex.data(), ey.data());
    g->SetName(name);
    g->SetMarkerStyle(mirror ? 24 : 20);
    g->SetLineWidth(2);
    return g;
}

// ==========================================================================
/**
 * @brief One canvas per folder: B_phi(lambda_z), the A/C yield asymmetry, and the moments.
 *
 * | Panel | Content | Expectation |
 * |---|---|---|
 * | 1 | B_phi vs lambda_z, jet, with its own mirror image overlaid | the two must lie on top of each other |
 * | 2 | B_phi vs lambda_z, jet and leading particle | must coincide: B_phi does not know the proxy |
 * | 3 | K = B_phi cosh(eta), with M0, M1 and the symmetrised M1 in a pave | M1sym consistent with zero |
 *
 * Panel 1: a detector-induced B_phi from an A/C-symmetric acceptance is even in eta_Lambda, so open and
 * filled markers should overlie. An odd part is either an A/C difference in the acceptance itself (not
 * in the yield, which panel 3 reports separately) or a genuine P_N; compare field polarities to tell.
 *
 * The numbers also go to stdout, since they are what the next stage of the chain consumes.
 *
 * @param taskDir        Top-level task TDirectory.
 * @param outDir         Output sub-directory for the canvas and the graphs.
 * @param f              Folder to process (see kFolders).
 * @param maxAbsCos      |cos(dTheta)| ceiling for a column to be used.
 * @param minEntriesCell Minimum entries for one 3D cell to enter a column fit.
 */
// ==========================================================================
static void MakeKernelMomentsCanvas(TDirectory* taskDir, TDirectory* outDir, const FolderSpec& f,
                                    double maxAbsCos = 0.95, double minEntriesCell = 30.)
{
    TDirectory* folderDir = GetDir(taskDir, f.name);
    if (!folderDir) { printf("WARNING: skipping kernel moments for '%s' (folder missing)\n", f.name); return; }

    TDirectory* kDir = GetDirPath(folderDir, "RingKernel");
    if (!kDir) { printf("WARNING: skipping kernel moments for '%s' (RingKernel missing)\n", f.name); return; }

    TProfile3D* h3Jet = static_cast<TProfile3D*>(
        SafeGet(kDir, "p3dRingObservableCosDeltaThetaVsJetZVsLambdaZ", false));
    TProfile3D* h3LeadP = static_cast<TProfile3D*>(
        SafeGet(kDir, "p3dRingObservableLeadPCosDeltaThetaVsLeadPZVsLambdaZ", false));
    if (!h3Jet) { printf("WARNING: skipping kernel moments for '%s' (jet 3D profile missing)\n", f.name); return; }

    const std::vector<KernelBphiPoint> ptsJet = ExtractBphiFromKernel3D(h3Jet, maxAbsCos, minEntriesCell);
    const std::vector<KernelBphiPoint> ptsLeadP =
        h3LeadP ? ExtractBphiFromKernel3D(h3LeadP, maxAbsCos, minEntriesCell)
                : std::vector<KernelBphiPoint>{};

    const KernelMoments mJet = ComputeKernelMoments(ptsJet);
    const KernelMoments mLeadP = ComputeKernelMoments(ptsLeadP);

    if (!mJet.valid) {
        printf("WARNING: skipping kernel moments for '%s' (no usable cos#theta_{#Lambda} slab)\n", f.name);
        return;
    }

    TGraphErrors* gJet    = MakeBphiGraph(ptsJet,   false, Form("gBphiJet_%s", f.name));
    TGraphErrors* gJetMir = MakeBphiGraph(ptsJet,   true,  Form("gBphiJetMirror_%s", f.name));
    TGraphErrors* gLeadP  = MakeBphiGraph(ptsLeadP, false, Form("gBphiLeadP_%s", f.name));
    if (!gJet) return;

    TCanvas* c = new TCanvas(Form("cKernelMoments_%s", f.name), "", 2100, 650);
    c->Divide(3, 1, 0.005, 0.005);

    // --- Panel 1: the Mz-symmetry test ---
    c->cd(1);
    SetupPadMargins();
    gPad->SetGridx(); gPad->SetGridy();
    gJet->SetMarkerColor(kAzure + 2); gJet->SetLineColor(kAzure + 2);
    gJet->SetTitle(" ;cos#theta_{#Lambda};B_{#varphi} [%]");
    gJet->Draw("AP");
    if (gJetMir) {
        gJetMir->SetMarkerColor(kAzure + 2); gJetMir->SetLineColor(kAzure + 2);
        gJetMir->Draw("P SAME");
    }
    {
        TLegend* leg = new TLegend(0.55, 0.78, 0.84, 0.90);
        leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->AddEntry(gJet, "B_{#varphi}", "lp");
        if (gJetMir) leg->AddEntry(gJetMir, "mirrored", "lp");
        leg->Draw();
    }
    AddLabel(0.5, 0.95, "Odd part: A/C acceptance asymmetry or genuine P_{N}", 0.045, 22);

    // --- Panel 2: the proxy-independence test ---
    c->cd(2);
    SetupPadMargins();
    gPad->SetGridx(); gPad->SetGridy();
    gJet->Draw("AP");
    if (gLeadP) {
        gLeadP->SetMarkerColor(kRed + 1); gLeadP->SetLineColor(kRed + 1);
        gLeadP->SetMarkerStyle(21);
        gLeadP->Draw("P SAME");
    }
    {
        TLegend* leg = new TLegend(0.55, 0.78, 0.84, 0.90);
        leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->AddEntry(gJet, "Leading jet", "lp");
        if (gLeadP) leg->AddEntry(gLeadP, "Leading particle", "lp");
        leg->Draw();
    }
    AddLabel(0.5, 0.95, "B_{#varphi} does not know the proxy", 0.045, 22);

    // --- Panel 3: the moments ---
    c->cd(3);
    SetupPadMargins();
    gPad->SetGridx(); gPad->SetGridy();
    gJet->Draw("AP");
    {
        TPaveText* pave = new TPaveText(0.15, 0.58, 0.88, 0.90, "NDC");
        pave->SetFillColor(kWhite);
        pave->SetFillStyle(1001);
        pave->SetBorderSize(1);
        pave->SetTextAlign(12);
        pave->SetTextSize(0.035);
        pave->AddText(Form("M_{0} = (%.4f #pm %.4f)%%", mJet.m0 * 100., mJet.m0Err * 100.));
        pave->AddText(Form("M_{1} = (%.4f #pm %.4f)%%", mJet.m1 * 100., mJet.m1Err * 100.));
        if (mJet.symmetrisable)
            pave->AddText(Form("M_{1}^{sym} = (%.4f #pm %.4f)%%", mJet.m1Sym * 100., mJet.m1SymErr * 100.));
        else
            pave->AddText("M_{1}^{sym}: cos#theta_{#Lambda} axis not mirror-pairable");
        pave->AddText(Form("#Lambda yield asymmetry = %+.3f", mJet.yieldAsym));
        if (mLeadP.valid)
            pave->AddText(Form("M_{0}^{LeadP} = (%.4f #pm %.4f)%%", mLeadP.m0 * 100., mLeadP.m0Err * 100.));
        pave->Draw();
    }
    AddLabel(0.5, 0.95, "Kernel moments", 0.045, 22);

    c->cd(0);
    AddLabel(0.5, 0.98, Form("%sB_{#varphi} and the kernel moments  --  %s", gNullTag.c_str(), f.label), 0.032, 22);

    WriteCanvas(c, outDir);

    outDir->cd();
    gJet->Write();
    if (gLeadP) gLeadP->Write();

    // These numbers are what the cross-config stage consumes, so they belong in the run log too:
    printf("    [%s] M0 = %.5f +- %.5f   M1 = %.5f +- %.5f   M1sym = %.5f +- %.5f   yieldAsym = %+.4f\n",
           f.name, mJet.m0, mJet.m0Err, mJet.m1, mJet.m1Err, mJet.m1Sym, mJet.m1SymErr, mJet.yieldAsym);
    if (mLeadP.valid)
        printf("    [%s] M0(LeadP) = %.5f +- %.5f   (must agree with M0 above)\n",
               f.name, mLeadP.m0, mLeadP.m0Err);

    delete c;
}


// ==========================================================================
// Section 8: symmetry decomposition of the 3D kernel profile.
//
// Two reflections act on the (c, t_z, lambda_z) cells of the 3D kernel profile, c = cos dTheta:
//   Mz  : reflection of the event through z = 0          (c, t_z, lambda_z) --> ( c, -t_z, -lambda_z)
//   A   : reference axis replaced by its antipode, t-->-t (c, t_z, lambda_z) --> (-c, -t_z,  lambda_z)
//   MzA : both                                            (c, t_z, lambda_z) --> (-c,  t_z, -lambda_z)
// With the identity they form a four-element group, so any cell function F splits exactly into four
// sectors labelled by its sign under Mz and under A:
//
//     F_ab = ( F + a F.Mz + b F.A + ab F.MzA ) / 4,      a, b = +-1.
//
// What lands where (derivations in the ring-geometry note, "Mirror Symmetries" and "The Antipodal
// Reflection of the Reference Axis"):
//
//   (+,+)  a genuine ring, averaged between dTheta and pi - dTheta
//   (+,-)  a genuine P_N; the near/away-odd part of a genuine ring; detector effects leaking through
//          an A/C-asymmetric acceptance
//   (-,-)  detector-induced contributions that do not know where the reference axis is
//   (-,+)  nothing: a null test
//
// Mz needs a barrel with equal acceptance on both sides (for the detector part); A needs an
// efficiency that does not know where the reference axis is, and nothing at all of the detector.
// Each tolerates exactly what the other forbids, which is what makes the null sector informative:
// a non-zero (-,+) is a detector effect that knows the reference axis, or a finite-bin residual.
// Validated on a synthetic profile carrying a detector-induced component, a genuine P_N and a genuine
// ring with different near- and away-side magnitudes: each landed in its sector, and the null stayed at
// the finite-bin level.
//
// Everything is done at fixed lambda_z, where the statements are exact per cell; the display then
// averages the lambda_z < 0 slabs only (see DecomposeKernel3D() for why that loses nothing). Integrating over lambda_z FIRST would not be exact:
// the phase-space weights of a cell and of its antipodal image differ as soon as the reference axis is
// correlated in azimuth with the Lambda.
// ==========================================================================

/// @brief Sector order used throughout Section 8: (+,+), (+,-), (-,-), (-,+) -- signs under (Mz, A).
static const int kSectorSign[4][2] = {{+1, +1}, {+1, -1}, {-1, -1}, {-1, +1}};

/// @brief The four sectors of one 3D kernel profile, averaged over lambda_z, plus per-sector chi2 vs 0.
struct KernelSectors {
    TH2D*  map[4]  = {nullptr, nullptr, nullptr, nullptr}; ///< Over (cos dTheta, t_z), order as kSectorSign
    double chi2[4] = {0., 0., 0., 0.};                     ///< Against zero, one term per group orbit
    int    ndf[4]  = {0, 0, 0, 0};
    bool   valid   = false;                                ///< At least one complete orbit was found
};

// ==========================================================================
/**
 * @brief Splits a 3D kernel profile into the four symmetry sectors, cell by cell at fixed lambda_z.
 *
 * Each cell is combined with its three images; the sector value is the signed quarter-sum, with
 * error a quarter of the quadrature sum of the four cell errors (the four members are distinct cells,
 * since mirror-pairable axes have no bin centred on zero). The chi2 against zero counts each orbit
 * once: the four members of an orbit carry the same sector values up to sign, so counting all four
 * would quadruple the ndf.
 *
 * The display maps average over the lambda_z < 0 slabs ONLY. Averaging over all of lambda_z would
 * cancel any sector content odd in lambda_z -- the -lambda_z c term of a detector-induced ring, for
 * one -- and hide it from the picture. Nothing is lost by using one half: Mz maps a cell at lambda_z to
 * one at -lambda_z, and a sector value there is exactly (sign under Mz) times the value here.
 *
 * @param h3             Kernel TProfile3D: x = cos(dTheta), y = proxy direction cosine, z = lambda_z.
 * @param minEntriesCell Every member of an orbit must clear this, or the orbit is skipped: a missing
 *                       mirror would otherwise enter a sector as a spurious signal.
 * @param tag            Name tag for the output histograms.
 */
// ==========================================================================
static KernelSectors DecomposeKernel3D(TProfile3D* h3, double minEntriesCell, const std::string& tag)
{
    KernelSectors out;
    if (!h3) return out;
    TAxis* ac = h3->GetXaxis();
    TAxis* at = h3->GetYaxis();
    TAxis* al = h3->GetZaxis();
    if (!IsMirrorPairable(ac) || !IsMirrorPairable(at) || !IsMirrorPairable(al)) return out;

    const int nc = ac->GetNbins(), nt = at->GetNbins(), nl = al->GetNbins();

    // Display maps share the profile's own (c, t_z) binning, edges copied so variable widths survive
    std::vector<double> edgesC(nc + 1), edgesT(nt + 1);
    for (int b = 0; b <= nc; ++b) edgesC[b] = (b < nc) ? ac->GetBinLowEdge(b + 1) : ac->GetBinUpEdge(nc);
    for (int b = 0; b <= nt; ++b) edgesT[b] = (b < nt) ? at->GetBinLowEdge(b + 1) : at->GetBinUpEdge(nt);

    static const char* kSectorName[4] = {"PP", "PM", "MM", "MP"};
    std::vector<double> sw(4 * nc * nt, 0.), swv(4 * nc * nt, 0.), sw2e2(4 * nc * nt, 0.);
    auto slot = [nc, nt](int s, int j, int i) { return (s * nc + (j - 1)) * nt + (i - 1); };

    for (int j = 1; j <= nc; ++j) {
        const int jc = nc + 1 - j;
        for (int i = 1; i <= nt; ++i) {
            const int it = nt + 1 - i;
            for (int k = 1; k <= nl; ++k) {
                const int kl = nl + 1 - k;
                // The orbit: the cell, its Mz image, its A image, its MzA image
                const int gb[4] = {h3->GetBin(j, i, k), h3->GetBin(j, it, kl),
                                   h3->GetBin(jc, it, k), h3->GetBin(jc, i, kl)};
                const int ix[4] = {j, j, jc, jc}, iy[4] = {i, it, it, i}, iz[4] = {k, kl, k, kl};

                bool complete = true;
                double F[4], E2 = 0.;
                for (int m = 0; m < 4 && complete; ++m) {
                    if (h3->GetBinEntries(gb[m]) < minEntriesCell) { complete = false; break; }
                    F[m] = h3->GetBinContent(ix[m], iy[m], iz[m]);
                    const double e = h3->GetBinError(ix[m], iy[m], iz[m]);
                    E2 += e * e;
                }
                if (!complete) continue;
                out.valid = true;

                const double err = 0.25 * std::sqrt(E2);
                const double n0 = h3->GetBinEntries(gb[0]);
                // Orbit representative: first half in both c and lambda_z (exactly one member qualifies)
                const bool representative = (2 * j <= nc) && (2 * k <= nl);

                const bool displayHalf = (2 * k <= nl); // lambda_z < 0 slabs; see the doc comment
                for (int s = 0; s < 4; ++s) {
                    const int a = kSectorSign[s][0], b = kSectorSign[s][1];
                    const double v = 0.25 * (F[0] + a * F[1] + b * F[2] + a * b * F[3]);
                    if (displayHalf) {
                        const int q = slot(s, j, i);
                        sw[q] += n0;
                        swv[q] += n0 * v;
                        sw2e2[q] += n0 * n0 * err * err;
                    }
                    if (representative && err > 0.) {
                        out.chi2[s] += (v / err) * (v / err);
                        out.ndf[s] += 1;
                    }
                }
            }
        }
    }
    if (!out.valid) return out;

    for (int s = 0; s < 4; ++s) {
        TH2D* h = new TH2D(Form("hKernelSector%s_%s", kSectorName[s], tag.c_str()), "",
                           nc, edgesC.data(), nt, edgesT.data());
        h->SetDirectory(nullptr);
        h->SetStats(0);
        for (int j = 1; j <= nc; ++j)
            for (int i = 1; i <= nt; ++i) {
                const int q = slot(s, j, i);
                if (sw[q] <= 0.) continue; // Never populated by a complete orbit: left empty (error 0)
                h->SetBinContent(j, i, swv[q] / sw[q]);
                h->SetBinError(j, i, std::sqrt(sw2e2[q]) / sw[q]);
            }
        out.map[s] = h;
    }
    return out;
}

// ==========================================================================
/**
 * @brief One canvas per folder and proxy: the four symmetry sectors of the 3D kernel profile.
 *
 * | Panel | Sector | Content |
 * |---|---|---|
 * | 1 | (+,+) | a genuine ring, averaged between dTheta and pi - dTheta |
 * | 2 | (+,-) | genuine P_N, near/away-odd ring, A/C-asymmetry leakage |
 * | 3 | (-,-) | detector-induced, axis-independent |
 * | 4 | (-,+) | null: must vanish |
 *
 * The sector maps are written next to the canvas, and each sector's chi2 against zero goes to stdout.
 *
 * @param taskDir        Top-level task TDirectory.
 * @param outDir         Output sub-directory.
 * @param f              Folder to process (see kFolders).
 * @param h3Name         Name of the 3D kernel profile inside <folder>/RingKernel/.
 * @param tag            Short proxy tag for names ("Jet", "LeadP").
 * @param proxyLabel     Human-readable proxy name for the canvas title.
 * @param minEntriesCell Minimum entries required of every member of an orbit.
 */
// ==========================================================================
static void MakeKernelSymmetryCanvas(TDirectory* taskDir, TDirectory* outDir, const FolderSpec& f,
                                     const char* h3Name, const char* tag, const char* proxyLabel,
                                     double minEntriesCell = 30.)
{
    TDirectory* folderDir = GetDir(taskDir, f.name);
    if (!folderDir) { printf("WARNING: skipping kernel symmetry for '%s' (folder missing)\n", f.name); return; }
    TDirectory* kDir = GetDirPath(folderDir, "RingKernel");
    if (!kDir) { printf("WARNING: skipping kernel symmetry for '%s' (RingKernel missing)\n", f.name); return; }

    TProfile3D* h3 = static_cast<TProfile3D*>(SafeGet(kDir, h3Name, false));
    if (!h3) { printf("WARNING: skipping kernel symmetry %s for '%s' (%s missing)\n", tag, f.name, h3Name); return; }

    const KernelSectors ks = DecomposeKernel3D(h3, minEntriesCell, Form("%s_%s", tag, f.name));
    if (!ks.valid) {
        printf("WARNING: skipping kernel symmetry %s for '%s': an axis does not mirror-pair, or no orbit "
               "has all four cells above %.0f entries\n", tag, f.name, minEntriesCell);
        return;
    }

    static const char* kBanner[4] = {
        "(+,+)  genuine ring, near/away averaged",
        "(+,-)  genuine P_{N}, A/C leakage",
        "(-,-)  detector-induced",
        "(-,+)  null: must vanish"};

    TCanvas* c = new TCanvas(Form("cKernelSymmetry%s_%s", tag, f.name), "", 2800, 650);
    c->Divide(4, 1, 0.005, 0.005);
    for (int s = 0; s < 4; ++s) {
        c->cd(s + 1);
        SetupPadMargins();
        DrawSymmetricColzH2(ks.map[s], " ;cos#Delta#theta;#hat{t}_{z};sector value, cos#theta_{#Lambda} < 0 [%]", true);
        AddLabel(0.5, 0.95, kBanner[s], 0.045, 22);
    }
    c->cd(0);
    AddLabel(0.5, 0.98, Form("%sKernel symmetry sectors (M_{z}, antipodal), %s  --  %s", gNullTag.c_str(), proxyLabel, f.label),
             0.032, 22);

    WriteCanvas(c, outDir);
    outDir->cd();
    for (int s = 0; s < 4; ++s) ks.map[s]->Write();

    // Consumed by the cross-config stage, so also in the run log:
    printf("    [%s, %s] sector chi2/ndf vs 0:  (+,+) %.1f/%d   (+,-) %.1f/%d   (-,-) %.1f/%d   (-,+) null %.1f/%d\n",
           f.name, tag, ks.chi2[0], ks.ndf[0], ks.chi2[1], ks.ndf[1], ks.chi2[2], ks.ndf[2], ks.chi2[3], ks.ndf[3]);

    delete c;
}

// ==========================================================================
// Section 9: KappaEff, the response coefficient per mass bin.
//
// kappa = 3 <u^2>/<w>, with u = R/prefactor and w = 1 (full ring) or n_z^2 (R_z). The four consumer
// profiles per proxy give <u^2>, <w>, <u^2 w> and <R w> per mass bin; the spreads come from the
// profiles themselves. What is drawn here is the RAW kappa per mass bin -- signal and background mixed --
// as a QA of its mass dependence. kappa_S belongs to the signal extraction, not to this file.
// ==========================================================================

/// @brief Per-bin spread sigma_y of a TProfile ("s" error option), restoring the caller's option.
static double ProfileSpread(TProfile* p, int bin)
{
    const TString opt = p->GetErrorOption();
    p->SetErrorOption("s");
    const double spread = p->GetBinError(bin);
    p->SetErrorOption(opt);
    return spread;
}

// ==========================================================================
/**
 * @brief kappa = 3 <num>/<den> per bin, with its first-order error.
 *
 * Both means come from the same candidates, so the error keeps the covariance:
 *   Var(kappa) = (9/n) [ s_U^2 - 2 (U/W) C + (U/W)^2 s_W^2 ] / W^2,   C = <u^2 w> - U W.
 * For the full ring w == 1, so s_W = C = 0 and this reduces to the standard error of 3 <u^2>.
 *
 * @param pNum        <u^2> vs mass.
 * @param pDen        <w> vs mass.
 * @param pNumDen     <u^2 w> vs mass.
 * @param newName     Name of the returned histogram.
 * @param minEntries  Bins with fewer candidates are left empty.
 * @return A TH1D on the mass axis, or nullptr if the three profiles are not aligned.
 */
// ==========================================================================
static TH1D* MakeKappaVsMass(TProfile* pNum, TProfile* pDen, TProfile* pNumDen, const char* newName, double minEntries = 50.)
{
    if (!pNum || !pDen || !pNumDen) return nullptr;
    const int n = pNum->GetNbinsX();
    if (pDen->GetNbinsX() != n || pNumDen->GetNbinsX() != n) return nullptr;

    // Same binning as the profiles, variable-width or not:
    const TAxis* ax = pNum->GetXaxis();
    TH1D* out = ax->GetXbins()->GetSize() ? new TH1D(newName, "", n, ax->GetXbins()->GetArray())
                                          : new TH1D(newName, "", n, ax->GetXmin(), ax->GetXmax());
    out->SetDirectory(nullptr);
    out->SetStats(0);

    for (int b = 1; b <= n; ++b) {
        const double entries = pNum->GetBinEntries(b);
        if (entries < minEntries) continue;
        const double U = pNum->GetBinContent(b);
        const double W = pDen->GetBinContent(b);
        if (W <= 0.) continue;
        const double sU = ProfileSpread(pNum, b);
        const double sW = ProfileSpread(pDen, b);
        const double C  = pNumDen->GetBinContent(b) - U * W;
        const double r  = U / W;
        const double var = 9. / entries * (sU * sU - 2. * r * C + r * r * sW * sW) / (W * W);
        out->SetBinContent(b, 3. * r);
        out->SetBinError(b, var > 0. ? std::sqrt(var) : 0.);
    }
    return out;
}

// ==========================================================================
/**
 * @brief One canvas per folder: raw kappa vs mass, one panel per proxy. Written for both ring modes.
 *
 * The kappa histograms are written next to the canvas, as hKappaVsMass<Proxy>_<folder>.
 */
// ==========================================================================
static void MakeKappaEffCanvas(TDirectory* taskDir, TDirectory* outDir, const FolderSpec& f)
{
    TDirectory* folderDir = GetDir(taskDir, f.name);
    if (!folderDir) { printf("WARNING: skipping KappaEff for '%s' (folder missing)\n", f.name); return; }
    TDirectory* kDir = GetDirPath(folderDir, "KappaEff", false);
    if (!kDir) { printf("WARNING: skipping KappaEff for '%s' (KappaEff/ missing: consumer older than the KappaEff moments?)\n", f.name); return; }

    static const char* kProxies[3] = {"LeadJet", "LeadP", "SubJet"};
    const std::string kappaSym = (gRingSym == "#it{R}") ? "#kappa_{eff}" : "#kappa_{z}";

    TCanvas* c = new TCanvas(Form("cKappaEff_%s", f.name), "", 2100, 650);
    c->Divide(3, 1, 0.005, 0.005);
    std::vector<TH1D*> written;

    for (int i = 0; i < 3; ++i) {
        TProfile* pNum    = static_cast<TProfile*>(SafeGet(kDir, Form("pKappaNum%sVsMass", kProxies[i]), false));
        TProfile* pDen    = static_cast<TProfile*>(SafeGet(kDir, Form("pKappaDen%sVsMass", kProxies[i]), false));
        TProfile* pNumDen = static_cast<TProfile*>(SafeGet(kDir, Form("pKappaNumTimesDen%sVsMass", kProxies[i]), false));
        TH1D* h = MakeKappaVsMass(pNum, pDen, pNumDen, Form("hKappaVsMass%s_%s", kProxies[i], f.name));
        c->cd(i + 1);
        SetupPadMargins();
        if (!h) {
            printf("WARNING: KappaEff %s for '%s': moment profiles missing or misaligned\n", kProxies[i], f.name);
            continue;
        }
        h->SetTitle((" ;m_{p#pi} (GeV/c^{2});raw " + kappaSym + " per mass bin").c_str());
        h->SetMarkerStyle(20);
        h->SetMarkerSize(0.6);
        h->Draw("E1");
        AddLabel(0.5, 0.95, kProxies[i], 0.045, 22);
        written.push_back(h);
    }
    c->cd(0);
    AddLabel(0.5, 0.98, Form("Raw %s vs mass (signal + background)  --  %s", kappaSym.c_str(), f.label), 0.032, 22);

    WriteCanvas(c, outDir);
    outDir->cd();
    for (auto* h : written) h->Write();
    delete c;
}

// ==========================================================================
// Section 10: R_z-only diagnostics (the README's "The longitudinal ring" has the derivations).
//
// Exactly, n_z = sin(Delta phi) sin(theta_t) sin(theta_Lambda) / sin(Delta theta), with
// Delta phi = phi_Lambda - phi_proxy as the consumer defines it. So the two halves of the Delta phi axis
// are the two sides of the Lambda ("east", "west"), and under R_z:
//   - the signal, weighted by n_z^2, is EVEN in Delta phi;
//   - the Lambda-kinematic fakes, weighted by n_z, are ODD in Delta phi.
// The odd part of <R_z>(Delta phi) is therefore pure leak, and it reaches the integrated value only
// through the east/west imbalance A_EW: <R_z> - symmetrized = A_EW x (half-difference).
// ==========================================================================

// ==========================================================================
/**
 * @brief The fold of a 1D profile about x = 0: 0.5 (f(x) +/- f(-x)) on the positive half-axis.
 *
 * The two partners hold disjoint candidates, so their errors add in quadrature.
 *
 * @param src        Profile with a mirror-pairable x axis (see IsMirrorPairable()).
 * @param odd        true: 0.5 (f(x) - f(-x)); false: 0.5 (f(x) + f(-x)).
 * @param minEntries Minimum entries required of BOTH partners.
 * @param newName    Name of the returned histogram.
 */
// ==========================================================================
static TH1D* MakeProfileFoldX1D(TProfile* src, bool odd, double minEntries, const char* newName)
{
    if (!src || !IsMirrorPairable(src->GetXaxis())) return nullptr;
    const int n = src->GetNbinsX();
    const int half = n / 2;

    // Edges of the positive half, which also covers variable-width axes:
    std::vector<double> edges(half + 1);
    for (int k = 0; k <= half; ++k) edges[k] = src->GetXaxis()->GetBinLowEdge(half + 1 + k);
    edges[half] = src->GetXaxis()->GetXmax();

    TH1D* out = new TH1D(newName, "", half, edges.data());
    out->SetDirectory(nullptr);
    out->SetStats(0);

    const double sign = odd ? -1.0 : +1.0;
    for (int k = 1; k <= half; ++k) {
        const int bPos = half + k;
        const int bNeg = half + 1 - k;
        if (src->GetBinEntries(bPos) < minEntries || src->GetBinEntries(bNeg) < minEntries) continue;
        const double ep = src->GetBinError(bPos);
        const double em = src->GetBinError(bNeg);
        out->SetBinContent(k, 0.5 * (src->GetBinContent(bPos) + sign * src->GetBinContent(bNeg)));
        out->SetBinError(k, 0.5 * std::sqrt(ep * ep + em * em));
    }
    return out;
}

// ==========================================================================
/**
 * @brief Integrated east/west summary of one <R_z>(Delta phi) profile.
 *
 * The halves are merged with TProfile::Rebin, which is exact (sums, not averages of means), so each half
 * carries its own candidates' mean and error. Bins:
 *   1 all, 2 Delta phi > 0, 3 Delta phi < 0, 4 symmetrized 0.5 (R+ + R-), 5 half-difference D = 0.5 (R+ - R-),
 *   6 leak in the integrated value, A_EW x D.
 * The weights in bins 4-6 are treated as fixed, which is safe while the counts are far better known than
 * the means. A_EW itself is returned through aEW.
 */
// ==========================================================================
static TH1D* MakeEastWestSummary(TProfile* src, const char* newName, double& aEW)
{
    aEW = 0.;
    if (!src || !IsMirrorPairable(src->GetXaxis())) return nullptr;
    const int n = src->GetNbinsX();

    TProfile* halves = static_cast<TProfile*>(src->Clone(Form("%s_halves_%d", src->GetName(), gCloneIdx++)));
    halves->SetDirectory(nullptr);
    halves->Rebin(n / 2);
    TProfile* all = static_cast<TProfile*>(src->Clone(Form("%s_all_%d", src->GetName(), gCloneIdx++)));
    all->SetDirectory(nullptr);
    all->Rebin(n);

    const double nNeg = halves->GetBinEntries(1), nPos = halves->GetBinEntries(2);
    if (nNeg <= 0. || nPos <= 0.) { delete halves; delete all; return nullptr; }
    aEW = (nPos - nNeg) / (nPos + nNeg);

    const double rPos = halves->GetBinContent(2), ePos = halves->GetBinError(2);
    const double rNeg = halves->GetBinContent(1), eNeg = halves->GetBinError(1);
    const double eHalf = 0.5 * std::sqrt(ePos * ePos + eNeg * eNeg);
    const double d = 0.5 * (rPos - rNeg);

    static const char* kLabels[6] = {"All", "#Delta#varphi > 0", "#Delta#varphi < 0", "Symmetrized",
                                     "Half-difference D", "A_{EW} #times D (leak)"};
    const double vals[6] = {all->GetBinContent(1), rPos, rNeg, 0.5 * (rPos + rNeg), d, aEW * d};
    const double errs[6] = {all->GetBinError(1), ePos, eNeg, eHalf, eHalf, std::fabs(aEW) * eHalf};

    TH1D* out = new TH1D(newName, "", 6, 0., 6.);
    out->SetDirectory(nullptr);
    out->SetStats(0);
    for (int b = 0; b < 6; ++b) {
        out->GetXaxis()->SetBinLabel(b + 1, kLabels[b]);
        out->SetBinContent(b + 1, vals[b]);
        out->SetBinError(b + 1, errs[b]);
    }
    delete halves;
    delete all;
    return out;
}

// ==========================================================================
/**
 * @brief One canvas per folder: the Delta phi fold of <R_z>, one column per proxy.
 *
 * Top pad: even (signal) and odd (leak) parts vs |Delta phi|. Bottom pad: the integrated east/west
 * summary, with A_EW in the banner. Fold and summary histograms are written next to the canvas.
 */
// ==========================================================================
static void MakeDeltaPhiFoldCanvas(TDirectory* taskDir, TDirectory* outDir, const FolderSpec& f, double minEntries = 50.)
{
    TDirectory* folderDir = GetDir(taskDir, f.name);
    if (!folderDir) { printf("WARNING: skipping Delta phi fold for '%s' (folder missing)\n", f.name); return; }

    static const char* kProfiles[3] = {"pRingObservableDeltaPhi", "pRingObservableLeadPDeltaPhi", "pRingObservable2ndJetDeltaPhi"};
    static const char* kProxies[3]  = {"LeadJet", "LeadP", "SubJet"};

    TCanvas* c = new TCanvas(Form("cDeltaPhiFold_%s", f.name), "", 2100, 1100);
    c->Divide(3, 2, 0.005, 0.005);
    std::vector<TH1*> written;

    for (int i = 0; i < 3; ++i) {
        TProfile* src = static_cast<TProfile*>(SafeGet(folderDir, kProfiles[i], false));
        if (!src || !IsMirrorPairable(src->GetXaxis())) {
            printf("WARNING: Delta phi fold %s for '%s': '%s' missing or its axis does not mirror-pair\n",
                   kProxies[i], f.name, kProfiles[i]);
            continue;
        }
        TH1D* even = MakeProfileFoldX1D(src, false, minEntries, Form("hRzDeltaPhiEven%s_%s", kProxies[i], f.name));
        TH1D* odd  = MakeProfileFoldX1D(src, true,  minEntries, Form("hRzDeltaPhiOdd%s_%s", kProxies[i], f.name));
        double aEW = 0.;
        TH1D* summary = MakeEastWestSummary(src, Form("hRzEastWest%s_%s", kProxies[i], f.name), aEW);
        if (!even || !odd || !summary) continue;

        c->cd(i + 1);
        SetupPadMargins();
        even->SetTitle(RingLabel(" ;|#Delta#varphi|;fold of <#it{R}>").c_str());
        even->SetMarkerStyle(20);
        even->SetMarkerColor(kBlue + 1);
        even->SetLineColor(kBlue + 1);
        odd->SetMarkerStyle(24);
        odd->SetMarkerColor(kRed + 1);
        odd->SetLineColor(kRed + 1);
        const double yMax = 1.3 * std::max({std::fabs(even->GetMaximum()), std::fabs(even->GetMinimum()),
                                            std::fabs(odd->GetMaximum()), std::fabs(odd->GetMinimum())});
        even->GetYaxis()->SetRangeUser(-yMax, yMax);
        even->Draw("E1");
        odd->Draw("E1 SAME");
        TLegend* leg = new TLegend(0.14, 0.76, 0.60, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(even, "even: signal", "lp");
        leg->AddEntry(odd, "odd: leak", "lp");
        leg->Draw();
        AddLabel(0.5, 0.95, kProxies[i], 0.045, 22);

        c->cd(i + 4);
        SetupPadMargins();
        summary->SetTitle(RingLabel(" ; ;<#it{R}>").c_str());
        summary->SetMarkerStyle(20);
        summary->Draw("E1");
        AddLabel(0.5, 0.95, Form("A_{EW} = %.4f", aEW), 0.045, 22);

        written.push_back(even);
        written.push_back(odd);
        written.push_back(summary);
        printf("    [%s, %s] A_EW = %+.4f   symmetrized = %+.3e +/- %.1e   leak = %+.3e +/- %.1e\n", f.name, kProxies[i],
               aEW, summary->GetBinContent(4), summary->GetBinError(4), summary->GetBinContent(6), summary->GetBinError(6));
    }
    c->cd(0);
    AddLabel(0.5, 0.99, RingLabel(Form("<#it{R}> folded about #Delta#varphi = 0 (east/west)  --  %s", f.label)).c_str(), 0.025, 22);

    WriteCanvas(c, outDir);
    outDir->cd();
    for (auto* h : written) h->Write();
    delete c;
}

// ==========================================================================
/**
 * @brief The four sectors of <R_z>(chi) under a: chi --> -chi and b: chi --> pi - chi.
 *
 *   f_(sa,sb)(chi) = 1/4 [ f(chi) + sa f(-chi) + sb f(pi - chi) + sa sb f(chi - pi) ],   chi in [0, pi/2].
 *
 * Sector order: (+,+) signal, (-,+) Lambda-kinematic fakes x A_EW, (-,-) towards-the-jet P_e, (+,-) null.
 * The four partners of an orbit are disjoint bins, so their errors add in quadrature. Requires a uniform
 * axis on [-pi, pi] with a bin count divisible by 4, so that every partner is a whole bin.
 *
 * @return false (and no histograms) if the axis does not qualify.
 */
// ==========================================================================
static bool DecomposeChi(TProfile* src, double minEntries, const std::string& tag, TH1D* out[4])
{
    const TAxis* ax = src ? src->GetXaxis() : nullptr;
    if (!ax || ax->GetXbins()->GetSize() || ax->GetNbins() % 4 != 0) return false;
    if (std::fabs(ax->GetXmin() + TMath::Pi()) > 1.e-4 || std::fabs(ax->GetXmax() - TMath::Pi()) > 1.e-4) return false;

    const int n = ax->GetNbins();
    const int quarter = n / 4;
    static const char* kSectorName[4] = {"PP", "MP", "MM", "PM"};
    static const int kSa[4] = {+1, -1, -1, +1};
    static const int kSb[4] = {+1, +1, -1, -1};

    for (int s = 0; s < 4; ++s) {
        out[s] = new TH1D(Form("hRzChiSector%s_%s", kSectorName[s], tag.c_str()), "", quarter, 0., 0.5 * TMath::Pi());
        out[s]->SetDirectory(nullptr);
        out[s]->SetStats(0);
    }
    for (int k = 1; k <= quarter; ++k) {
        const double chi = out[0]->GetXaxis()->GetBinCenter(k);
        // Partner bins of the orbit {chi, -chi, pi - chi, chi - pi}; centres map onto centres exactly:
        const int bins[4] = {src->FindFixBin(chi), src->FindFixBin(-chi), src->FindFixBin(TMath::Pi() - chi), src->FindFixBin(chi - TMath::Pi())};
        bool enough = true;
        for (int j = 0; j < 4; ++j) enough = enough && (src->GetBinEntries(bins[j]) >= minEntries);
        if (!enough) continue;

        double err2 = 0.;
        for (int j = 0; j < 4; ++j) err2 += src->GetBinError(bins[j]) * src->GetBinError(bins[j]);
        for (int s = 0; s < 4; ++s) {
            const double v = src->GetBinContent(bins[0]) + kSa[s] * src->GetBinContent(bins[1])
                           + kSb[s] * src->GetBinContent(bins[2]) + kSa[s] * kSb[s] * src->GetBinContent(bins[3]);
            out[s]->SetBinContent(k, 0.25 * v);
            out[s]->SetBinError(k, 0.25 * std::sqrt(err2));
        }
    }
    return true;
}

// ==========================================================================
/**
 * @brief One canvas per proxy: the four chi sectors of <R_z>, each fitted with its expected shape.
 *
 * | Panel | Sector | Shape fitted | Content |
 * |---|---|---|---|
 * | 1 | (+,+) | A sin^2(chi)  | signal (plus field-odd instrumental parts only) |
 * | 2 | (-,+) | B sin(chi)    | Lambda-kinematic fakes (B_theta, HEE) x A_EW |
 * | 3 | (-,-) | C sin(2 chi)  | towards-the-jet P_e, parity-odd, invisible to mixing |
 * | 4 | (+,-) | 0             | null |
 *
 * The shapes are those of the chi weights with theta_Lambda averaged; the amplitudes are summaries,
 * not measurements. Reads RzDiagnostics/ at task level, so there is no folder loop.
 */
// ==========================================================================
static void MakeChiSectorCanvas(TDirectory* taskDir, TDirectory* outDir, const char* proxy, double minEntries = 50.)
{
    TDirectory* rzDir = GetDir(taskDir, "RzDiagnostics", false);
    TProfile* src = rzDir ? static_cast<TProfile*>(SafeGet(rzDir, Form("pRingVsChi%s", proxy), false)) : nullptr;
    if (!src) { printf("WARNING: skipping chi sectors %s (RzDiagnostics/pRingVsChi%s missing)\n", proxy, proxy); return; }

    TH1D* sec[4] = {nullptr, nullptr, nullptr, nullptr};
    if (!DecomposeChi(src, minEntries, proxy, sec)) {
        printf("WARNING: skipping chi sectors %s: the chi axis must be uniform on [-pi, pi] with a multiple of 4 bins\n", proxy);
        return;
    }

    static const char* kBanner[4] = {"(+,+)  signal", "(-,+)  #Lambda-kinematic fakes #times A_{EW}",
                                     "(-,-)  towards-the-jet P_{e}", "(+,-)  null"};
    static const char* kShape[4]  = {"[0]*sin(x)*sin(x)", "[0]*sin(x)", "[0]*sin(2*x)", "[0]"};

    TCanvas* c = new TCanvas(Form("cRzChiSectors%s", proxy), "", 2800, 650);
    c->Divide(4, 1, 0.005, 0.005);
    for (int s = 0; s < 4; ++s) {
        c->cd(s + 1);
        SetupPadMargins();
        sec[s]->SetTitle(RingLabel(" ;#chi;sector of <#it{R}>").c_str());
        sec[s]->SetMarkerStyle(20);
        sec[s]->Draw("E1");
        // The null sector is fitted with a constant, whose compatibility with zero is the test:
        TF1* fit = new TF1(Form("fRzChi%s_%d", proxy, s), kShape[s], 0., 0.5 * TMath::Pi());
        sec[s]->Fit(fit, "Q0R");
        fit->SetLineColor(kRed + 1);
        fit->Draw("SAME");
        AddLabel(0.5, 0.95, kBanner[s], 0.040, 22);
        AddLabel(0.55, 0.86, Form("amp = %.2e #pm %.1e,  #chi^{2}/ndf = %.1f/%d", fit->GetParameter(0), fit->GetParError(0),
                                  fit->GetChisquare(), fit->GetNDF()), 0.035, 22);
        printf("    [%s] chi sector %d: amp = %+.3e +/- %.1e   chi2/ndf = %.1f/%d\n", proxy, s,
               fit->GetParameter(0), fit->GetParError(0), fit->GetChisquare(), fit->GetNDF());
    }
    c->cd(0);
    AddLabel(0.5, 0.98, RingLabel(Form("<#it{R}> vs #chi: symmetry sectors under (#chi #rightarrow -#chi, #chi #rightarrow #pi - #chi)  --  %s", proxy)).c_str(), 0.032, 22);

    WriteCanvas(c, outDir);
    outDir->cd();
    for (auto* h : sec) h->Write();
    delete c;
}

// ==========================================================================
/**
 * @brief One canvas: <R_z> and <R_perp> vs Phi_AEE, per proxy. The AEE should live entirely in R_perp.
 *
 * <R_z>(Phi_AEE) is a genuine null: p*_z is independent of phi*_p for isotropic decays, so any
 * structure is instrumental and E-W imbalanced. A constant is fitted to it and its chi2 printed.
 * Both profiles fill on the same candidates (valid proxy, fake-pol diagnostics QA on). The R_z one
 * lives under HelicityEfficiencyQA/, so the comparison is skipped quietly when that QA was off.
 */
// ==========================================================================
static void MakePhiAeeFlatnessCanvas(TDirectory* taskDir, TDirectory* outDir)
{
    TDirectory* rzDir  = GetDir(taskDir, "RzDiagnostics", false);
    TDirectory* aeeDir = GetDirPath(taskDir, "HelicityEfficiencyQA/PhiLambdaPhiProtonStar", false);
    if (!rzDir || !aeeDir) { printf(" (Phi_AEE flatness skipped: HelicityEfficiencyQA or RzDiagnostics absent in this config)\n"); return; }

    static const char* kProxies[2] = {"LeadJet", "LeadP"};
    TCanvas* c = new TCanvas("cRzPhiAeeFlatness", "", 1400, 650);
    c->Divide(2, 1, 0.005, 0.005);

    for (int i = 0; i < 2; ++i) {
        TProfile* pZ    = static_cast<TProfile*>(SafeGet(aeeDir, Form("pRingObservable%sVsPhiLambdaLikePhiProtonStar", kProxies[i]), false));
        TProfile* pPerp = static_cast<TProfile*>(SafeGet(rzDir, Form("pRingPerp%sVsPhiAEE", kProxies[i]), false));
        if (!pZ || !pPerp) { printf("WARNING: Phi_AEE flatness %s: a profile is missing\n", kProxies[i]); continue; }

        c->cd(i + 1);
        SetupPadMargins();
        TProfile* hPerp = static_cast<TProfile*>(pPerp->Clone(Form("pPerpClone%s_%d", kProxies[i], gCloneIdx++)));
        TProfile* hZ    = static_cast<TProfile*>(pZ->Clone(Form("pZClone%s_%d", kProxies[i], gCloneIdx++)));
        hPerp->SetDirectory(nullptr);
        hZ->SetDirectory(nullptr);
        hPerp->SetTitle(" ;#phi_{#Lambda-like}-#phi_{p-like}^{*};< #it{R}_{#perp} >,  < #it{R}_{z} >");
        hPerp->SetMarkerStyle(24);
        hPerp->SetMarkerColor(kGray + 2);
        hPerp->SetLineColor(kGray + 2);
        hZ->SetMarkerStyle(20);
        hZ->SetMarkerColor(kBlue + 1);
        hZ->SetLineColor(kBlue + 1);
        hPerp->Draw("E1");
        hZ->Draw("E1 SAME");

        TF1* flat = new TF1(Form("fRzFlat%s", kProxies[i]), "[0]", hZ->GetXaxis()->GetXmin(), hZ->GetXaxis()->GetXmax());
        hZ->Fit(flat, "Q0");
        const double pVal = flat->GetNDF() > 0 ? TMath::Prob(flat->GetChisquare(), flat->GetNDF()) : 0.;
        TLegend* leg = new TLegend(0.14, 0.76, 0.60, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(hPerp, "< #it{R}_{#perp} >: carries the AEE", "lp");
        leg->AddEntry(hZ, Form("< #it{R}_{z} >: flat, #chi^{2}/ndf = %.1f/%d (p = %.2f)", flat->GetChisquare(), flat->GetNDF(), pVal), "lp");
        leg->Draw();
        AddLabel(0.5, 0.95, kProxies[i], 0.045, 22);
        printf("    [%s] R_z vs Phi_AEE constant: %+.3e +/- %.1e   chi2/ndf = %.1f/%d  (p = %.3f)\n", kProxies[i],
               flat->GetParameter(0), flat->GetParError(0), flat->GetChisquare(), flat->GetNDF(), pVal);
    }
    c->cd(0);
    AddLabel(0.5, 0.98, "#Phi_{AEE} flatness: the AEE lives in #it{R}_{#perp}, not in #it{R}_{z}", 0.032, 22);
    WriteCanvas(c, outDir);
    delete c;
}

// ==========================================================================
/**
 * @brief Resolves which of the kFolders members actually exist in the input file.
 *
 * Probes quietly (warnIfMissing = false) and prints one "present / absent" summary line, so that a
 * config enabling only some families costs two lines of log instead of one warning cascade per
 * drawing section. Every section afterwards iterates the returned vector and can therefore treat any
 * further missing object as a genuine anomaly worth a loud WARNING if needed.
 *
 * @param taskDir Top-level task TDirectory.
 * @param ok      Set to false when a mandatory folder is absent -- left untouched otherwise.
 * @return The subset of kFolders present in the file, in registry order, inside a vector.
 */
// ==========================================================================
static std::vector<FolderSpec> ScanPresentFolders(TDirectory* taskDir, bool& ok)
{
    std::vector<FolderSpec> present;
    std::string absentList;

    for (const auto& f : kFolders) {
        if (GetDir(taskDir, f.name, false)) {
            present.push_back(f);
            continue;
        }
        if (f.mandatory) {
            printf("ERROR: mandatory folder '%s' not found in '%s'\n", f.name, kTaskDir);
            ok = false;
            continue;
        }
        if (!absentList.empty()) absentList += ", ";
        absentList += f.name;
    }

    // Built from the vector rather than from kFolders so the two lines cannot disagree:
    std::string presentList;
    for (const auto& f : present) {
        if (!presentList.empty()) presentList += ", ";
        presentList += f.name;
    }

    printf(" Families present : %s\n", presentList.empty() ? "(none)" : presentList.c_str());
    if (!absentList.empty())
        printf(" Families absent  : %s  (optional, not enabled in this config)\n", absentList.c_str());

    return present;
}

// ---------------------------------------------------------
// Main Function
// ---------------------------------------------------------
int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <inputFilePath (ConsumerResults_*.root)> <outputDirectory>\n";
        return 1;
    }

    const char* inFileStr = argv[1];
    const char* outDirStr = argv[2];

    // --- Dynamic Output Filename (mirrors makeCumulativeDCAdauProfile.cxx) ---
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
        filename.replace(pos, prefixToReplace.length(), "AuxiliaryPerConfigPlots_");
    } else {
        filename = "AuxiliaryPerConfigPlots_" + filename;
    }

    // Normalize the output directory to always end with exactly one trailing slash
    std::string outDir(outDirStr);
    if (outDir.empty() || outDir.back() != '/') outDir += "/";

    // Create the output directory (and any missing parents) if it does not exist yet.
    // kTRUE requests recursive creation, mirroring signalExtractionRing.cxx's convention.
    gSystem->mkdir(outDir.c_str(), kTRUE);

    std::string outFileStr = outDir + filename;

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

    // --- Ring definition: full ring or R_z, from the file suffix (see the README) ---
    // RzDiagnostics/ is booked only by useRingZ, so it cross-checks the suffix: a mislabelled file would be
    // drawn, and then compared downstream, as the wrong observable. That is fatal, not a warning.
    const bool ringZMode = filename.find(kRingZTag) != std::string::npos;
    const bool hasRzDiagnostics = GetDir(taskDir, "RzDiagnostics", false) != nullptr;
    if (ringZMode != hasRzDiagnostics) {
        std::cerr << " ERROR: the file suffix says " << (ringZMode ? "R_z" : "full ring") << ", but RzDiagnostics/ is "
                  << (hasRzDiagnostics ? "present" : "absent") << ". Check the useRingZ / suffix pairing of " << inPath << std::endl;
        fOut->Close(); fIn->Close();
        return 1;
    }
    if (ringZMode) {
        gRingSym = "#it{R}_{z}";
        gNullTag = "[null test for #it{R}_{z}]  ";
    }
    std::cout << " Ring definition: " << (ringZMode ? "R_z = P_z n_z (useRingZ)" : "full ring") << "\n";

    // --- Resolve which kinematic-cut families this config actually produced ---
    // Done once, before any drawing, so an absent optional family costs one summary line instead of
    // one warning per drawing section. A missing MANDATORY family is fatal: the input file is not the
    // consumer output it claims to be, and drawing from it would be meaningless rather than merely partial.
    bool foldersOk = true;
    const std::vector<FolderSpec> folders = ScanPresentFolders(taskDir, foldersOk);
    if (!foldersOk) {
        std::cerr << " ERROR: mandatory kinematic-cut family missing in " << inPath << std::endl;
        fOut->Close(); fIn->Close();
        return 1;
    }

    // --- Ring-observable panel tables, one per coordinate system ---
    // Built here rather than at namespace scope because they depend on the axis-title fragments and
    // are cheap to construct; the plane descriptions are written once and shared by both jet proxies.
    const std::vector<std::string> labSuffixes  = {"PxPy", "PzPx", "PyPz"};
    const std::vector<std::string> labX         = {kPx, kPz, kPy};
    const std::vector<std::string> labY         = {kPy, kPx, kPz};
    const std::vector<std::string> labBanners   = {"X-Y plane", "Z-X plane", "Y-Z plane"};

    const std::vector<std::string> aeeSuffixes  = {"PxAeePyAee", "PzPxAee", "PyAeePz"};
    const std::vector<std::string> aeeX         = {kPxA, kPz, kPyA};
    const std::vector<std::string> aeeY         = {kPyA, kPxA, kPz};
    const std::vector<std::string> aeeBanners   = {"AEE plane (radius p_{T}, azimuth #Phi_{AEE})", "Z-X_{AEE} plane", "Y_{AEE}-Z plane"};

    const std::vector<std::string> jetSuffixes  = {"PxPyPrimeJet", "PzPxPrimeJet", "PyPzPrimeJet"};
    const std::vector<std::string> jetX         = {kPxJ, kPzJ, kPyJ};
    const std::vector<std::string> jetY         = {kPyJ, kPxJ, kPzJ};
    const std::vector<std::string> jetBanners   = {"X'-Y' plane (transverse to the jet)", "Z'-X' plane", "Y'-Z' plane"};

    // --- Create organized output subfolders ---
    TDirectory* dirVectorField = fOut->mkdir("VectorField_Canvases");
    TDirectory* dirRingObs2D   = fOut->mkdir("RingObservable2D_Canvases");
    TDirectory* dirAeeMaps     = fOut->mkdir("AeeAcceptance_Canvases");
    TDirectory* dirCounts      = fOut->mkdir("Counts_Canvases");
    TDirectory* dirAeeFold     = fOut->mkdir("AeeFold_Canvases");
    TDirectory* dirRingKernel  = fOut->mkdir("RingKernel_Canvases");
    TDirectory* dirKernelMoments = fOut->mkdir("KernelMoments_Canvases");
    TDirectory* dirKernelSymmetry = fOut->mkdir("KernelSymmetry_Canvases");
    TDirectory* dirKappaEff    = fOut->mkdir("KappaEff_Canvases");
    TDirectory* dirRzDiag      = ringZMode ? fOut->mkdir("RzDiagnostics_Canvases") : nullptr; // Never an empty folder in full-ring files

    // 1. Polarization vector-field canvases, one per coordinate system
    std::cout << " -> Drawing vector-field canvases...\n";
    for (const auto& f : folders) {
        MakePanelCanvas(taskDir, dirVectorField, f, "PolMaps/Lab", "cVectorFieldLab",
                        "<P*> polarization vector field, lab frame", kPanelsVectorLab);
        MakePanelCanvas(taskDir, dirVectorField, f, "PolMaps/PrimeJet", "cVectorFieldPrimeJet",
                        "<P*> polarization vector field, jet frame", kPanelsVectorPrimeJet);
        MakePanelCanvas(taskDir, dirVectorField, f, "PolMaps/PrimeV0", "cVectorFieldPrimeV0",
                        "<P*> polarization in the #Lambda production plane", kPanelsVectorPrimeV0);
    }

    // 2. Ring-observable 2D scalar canvases, per proxy and per coordinate system
    std::cout << " -> Drawing ring-observable 2D canvases...\n";
    for (const auto& f : folders) {
        MakePanelCanvas(taskDir, dirRingObs2D, f, "", "cRingObservable2D", RingLabel("Ring observable <#it{R}>"),
                        MakeRingPanels("p2dRingObservable", labSuffixes, labX, labY, labBanners));
        MakePanelCanvas(taskDir, dirRingObs2D, f, "", "cRingObservableAee2D", RingLabel("Ring observable <#it{R}>, AEE planes"),
                        MakeRingPanels("p2dRingObservable", aeeSuffixes, aeeX, aeeY, aeeBanners));
        MakePanelCanvas(taskDir, dirRingObs2D, f, "", "cRingObservablePrimeJet2D", RingLabel("Ring observable <#it{R}>, jet frame"),
                        MakeRingPanels("p2dRingObservable", jetSuffixes, jetX, jetY, jetBanners));
    }

    // 3. Ring-observable 2D scalar canvases (Leading Particle Proxy)
    std::cout << " -> Drawing ring-observable (Leading Particle) 2D canvases...\n";
    for (const auto& f : folders)
        MakePanelCanvas(taskDir, dirRingObs2D, f, "", "cRingObservableLeadP2D", RingLabel("Ring observable <#it{R}>_{LeadP}"),
                        MakeRingPanels("p2dRingObservableLeadP", labSuffixes, labX, labY, labBanners));

    // 4. AEE acceptance maps, one canvas per plane
    std::cout << " -> Drawing AEE acceptance canvases...\n";
    for (const auto& f : folders) {
        MakePanelCanvas(taskDir, dirAeeMaps, f, "PolMaps/Aee", "cAeeMapsXY",
                        "AEE acceptance, (x_{AEE}, y_{AEE}) plane", kPanelsAeeXY);
        MakePanelCanvas(taskDir, dirAeeMaps, f, "PolMaps/Aee", "cAeeMapsZX",
                        "AEE acceptance, (z, x_{AEE}) plane", kPanelsAeeZX);
        MakePanelCanvas(taskDir, dirAeeMaps, f, "PolMaps/Aee", "cAeeMapsYZ",
                        "AEE acceptance, (y_{AEE}, z) plane", kPanelsAeeYZ);
    }

    // 5. Candidate counts over the lab and jet-frame momentum planes
    std::cout << " -> Drawing candidate-count canvases...\n";
    for (const auto& f : folders) {
        MakePanelCanvas(taskDir, dirCounts, f, "PolMaps/Lab", "cCountsLab",
                        "Candidate counts, lab frame", kPanelsCountsLab);
        MakePanelCanvas(taskDir, dirCounts, f, "PolMaps/PrimeJet", "cCountsPrimeJet",
                        "Candidate counts, jet frame", kPanelsCountsPrimeJet);
    }

    // 6. AEE fold: the odd-in-Phi_AEE part of the acceptance, plus its parity nulls
    std::cout << " -> Drawing AEE fold canvases...\n";
    for (const auto& f : folders)
        MakeAeeFoldCanvas(taskDir, dirAeeFold, f);

    // 7. Ring projection kernel: per-slice linear fits that isolate M0
    std::cout << " -> Fitting ring-kernel slices...\n";
    for (const auto& f : folders)
        MakeRingKernelCanvas(taskDir, dirRingKernel, f);

    // 8. Kernel moments: B_phi per cos(theta_Lambda) slab, and the M0 / M1 it implies
    std::cout << " -> Extracting kernel moments...\n";
    for (const auto& f : folders)
        MakeKernelMomentsCanvas(taskDir, dirKernelMoments, f);

    // 9. Kernel symmetry sectors, per proxy: detector-induced vs genuine, and one null
    std::cout << " -> Decomposing kernel symmetry sectors...\n";
    for (const auto& f : folders) {
        MakeKernelSymmetryCanvas(taskDir, dirKernelSymmetry, f,
                                 "p3dRingObservableCosDeltaThetaVsJetZVsLambdaZ", "Jet", "leading jet");
        MakeKernelSymmetryCanvas(taskDir, dirKernelSymmetry, f,
                                 "p3dRingObservableLeadPCosDeltaThetaVsLeadPZVsLambdaZ", "LeadP", "leading particle");
    }

    // 10. KappaEff: raw kappa per mass bin, both ring definitions
    std::cout << " -> Drawing KappaEff canvases...\n";
    for (const auto& f : folders)
        MakeKappaEffCanvas(taskDir, dirKappaEff, f);

    // 11. R_z only: east/west fold per folder, then the task-level chi sectors and Phi_AEE flatness
    if (ringZMode) {
        std::cout << " -> Drawing R_z diagnostics...\n";
        for (const auto& f : folders)
            MakeDeltaPhiFoldCanvas(taskDir, dirRzDiag, f);
        for (const char* proxy : {"LeadJet", "LeadP", "SubJet"})
            MakeChiSectorCanvas(taskDir, dirRzDiag, proxy);
        MakePhiAeeFlatnessCanvas(taskDir, dirRzDiag);
    }

    // ===========================================================
    // ADD MORE POST-PROCESSING SECTIONS HERE.
    //
    // This file is meant to be the general home for per-config
    // derivative plots that live outside O2Physics -- i.e. anything
    // that reads a finished ConsumerResults_*.root and produces
    // additional canvases/histograms from it, but does not need the
    // cross-config/wagon-level aggregation that auxiliarySummaryPlots.cxx
    // handles. Follow the same pattern as the sections above:
    //   - a FolderSpec-style loop if the new plot is booked per
    //     kinematic-cut folder -- iterate the resolved `folders` vector,
    //     NOT kFolders, so absent optional families stay silent,
    //   - a dedicated TDirectory in fOut via fOut->mkdir(...),
    //   - a PanelSpec table plus a MakePanelCanvas(...) call, if the new
    //     plot is a strip of momentum-plane maps. Anything that needs a
    //     different layout gets its own Make...Canvas() above main().
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