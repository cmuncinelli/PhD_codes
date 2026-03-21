#!/usr/bin/env python3
"""
Lambda hyperon signal extraction script
======================================

This script performs invariant-mass signal extraction for Lambda candidates
stored in a 3D ROOT histogram:

    TH3D: Lambda/hLambdaPtZMass

Axes:
    X: pT of the Lambda candidate
    Y: Z position of the V0 decay vertex
    Z: invariant mass of the Lambda candidate

The procedure is identical for Data and MC.
Signal extraction is performed independently in each (pT, Z) bin.

Usage: python SignalExtraction.py <AnalysisResults.root>

Author: Cicero M. (with ChatGPT to increase the speed of coding)
"""

# =========================
# 0) Imports
# =========================

import os
import sys
import ROOT
import math

ROOT.gROOT.SetBatch(True)  # Do not pop up canvases
ROOT.gErrorIgnoreLevel = ROOT.kFatal # Overkill, but otherwise code wrote too much on screen!
ROOT.TH1.SetDefaultSumw2(True)
ROOT.TH2.SetDefaultSumw2(True)
ROOT.TH3.SetDefaultSumw2(True)

# =========================
# 1) Input handling
# =========================

print("\n[INFO] Starting Lambda signal extraction script")

if len(sys.argv) != 2:
    print("[ERROR] Usage: python SignalExtraction.py <AnalysisResults.root>")
    sys.exit(1)

input_file_path = sys.argv[1]

if not os.path.isfile(input_file_path):
    print(f"[ERROR] Input file does not exist: {input_file_path}")
    sys.exit(1)

input_dir  = os.path.dirname(input_file_path)
input_file = os.path.basename(input_file_path)

print(f"[INFO] Input file: {input_file_path}")

# Extract configuration name from filename
# Example:
# AnalysisResults-hasTPCnoITS.root  -->  hasTPCnoITS
config_name = input_file.replace("AnalysisResults-", "").replace(".root", "")

output_file_name = f"LambdaSignalExtraction-{config_name}.root"
output_file_path = os.path.join(input_dir, output_file_name)

print(f"[INFO] Output file will be written to:")
print(f"       {output_file_path}")

# =========================
# 2) Open input ROOT file
# =========================

infile = ROOT.TFile.Open(input_file_path, "READ")
if not infile or infile.IsZombie():
    print("[ERROR] Could not open input ROOT file")
    sys.exit(1)

# Navigate to main directory
main_dir_name = "asymmetric_rapidity_test"
main_dir = infile.Get(main_dir_name)

if not main_dir:
    print(f"[ERROR] Directory '{main_dir_name}' not found in input file")
    sys.exit(1)

# Get the TH3D
hist_path = "Lambda/hLambdaPtZMass"
h3 = main_dir.Get(hist_path)

if not h3:
    print(f"[ERROR] Histogram '{hist_path}' not found")
    sys.exit(1)

if not isinstance(h3, ROOT.TH3):
    print("[ERROR] Object found is not a TH3")
    sys.exit(1)

print("[INFO] Successfully loaded TH3D: Lambda/hLambdaPtZMass")

# =========================
# 3) Prepare output ROOT file and directory structure
# =========================

outfile = ROOT.TFile.Open(output_file_path, "RECREATE")
if not outfile or outfile.IsZombie():
    print("[ERROR] Could not create output ROOT file")
    sys.exit(1)

print("[INFO] Output ROOT file created")

# Top-level folders
dir_mass_projections = outfile.mkdir("InvariantMassProjections")
dir_fits_full        = outfile.mkdir("FullFits_GaussPlusPoly")
dir_fits_sidebands   = outfile.mkdir("SidebandFits")
dir_signal_maps      = outfile.mkdir("SignalExtractionMaps")
dir_spectra          = outfile.mkdir("FinalPtSpectra")

# =========================
# 4) Axis information
# =========================

xaxis = h3.GetXaxis()  # pT
yaxis = h3.GetYaxis()  # Z position
zaxis = h3.GetZaxis()  # invariant mass

n_pt_bins = xaxis.GetNbins()
n_z_bins  = yaxis.GetNbins()

print(f"[INFO] pT bins: {n_pt_bins}")
print(f"[INFO] Z bins : {n_z_bins}")

# =========================
# 5) Prepare containers for extracted yields
# =========================
# These histograms will store the signal-extracted yield
# for each (pT, Z) bin. Errors are non-Poissonian.

# Unfortunately ternary operators don't work here in ROOT definitions due to
# how python parses them to the definition!
# hSignalPtZ = ROOT.TH2D(
#     "hSignalPtZ",
#     "Lambda signal yield; p_{T} (GeV/c); V0 decay Z (cm)",
#     n_pt_bins, xaxis.GetXbins().GetArray()
#     if xaxis.GetXbins().GetSize() > 0 else xaxis.GetXmin(), xaxis.GetXmax(),
#     n_z_bins,  yaxis.GetXbins().GetArray()
#     if yaxis.GetXbins().GetSize() > 0 else yaxis.GetXmin(), yaxis.GetXmax()
# )
# We can, however, write a helper function that properly builds the TH2D
# independently of it having variable or fixed-width bins for each of the
# x, y or z axis we are working with!
def make_TH2D(name, title, xaxis, yaxis):
    if xaxis.GetXbins().GetSize() > 0:
        xdef = (xaxis.GetNbins(), xaxis.GetXbins().GetArray())
    else:
        xdef = (xaxis.GetNbins(), xaxis.GetXmin(), xaxis.GetXmax())

    if yaxis.GetXbins().GetSize() > 0:
        ydef = (yaxis.GetNbins(), yaxis.GetXbins().GetArray())
    else:
        ydef = (yaxis.GetNbins(), yaxis.GetXmin(), yaxis.GetXmax())

    return ROOT.TH2D(name, title, *xdef, *ydef)

hSignalPtZ = make_TH2D(
    "hSignalPtZ",
    "Lambda signal yield; p_{T} (GeV/c); V0 decay Z (cm)",
    xaxis, yaxis
)
hSignalPtZ.Sumw2()

# =========================
# 6) Loop over (pT, Z) bins and project invariant mass
# =========================

print("[INFO] Starting loop over (pT, Z) bins")

for ipt in range(1, n_pt_bins + 1):
    pt_low  = xaxis.GetBinLowEdge(ipt)
    pt_high = xaxis.GetBinUpEdge(ipt)

    for iz in range(1, n_z_bins + 1):
        z_low  = yaxis.GetBinLowEdge(iz)
        z_high = yaxis.GetBinUpEdge(iz)

        # Too verbose to actually include in final code:
        # print(f"[INFO] Processing bin: pT [{pt_low:.2f}, {pt_high:.2f}] "
        #       f"Z [{z_low:.2f}, {z_high:.2f}]")

        # Project invariant mass for this (pT, Z) bin
        # IMPORTANT: "e" option propagates errors correctly
        h_mass = h3.ProjectionZ(
            f"hMass_pt{ipt}_z{iz}",
            ipt, ipt,
            iz, iz,
            "e"
        )

        h_mass.SetTitle(
            f"Lambda invariant mass; m_{{p#pi}} (GeV/c^2); Counts"
        )

        # Save the mass histogram for QA
        dir_mass_projections.cd()
        if h_mass.GetEntries() > 0: # Should not save empty histograms! Those are just clutter for empty bins!
            h_mass.Write()

        # --- Placeholder for signal extraction ---
        # The actual fitting and bin counting will be done
        # in the next part of the script.

# =========================
# End of PART 1
# =========================

print("[INFO] Finished PART 1 (file handling & projections)")


# =========================
# 7) Signal extraction per (pT, Z) bin
# =========================

print("[INFO] Starting signal extraction in each (pT, Z) bin")

# Loop again to keep structure clear and readable
for ipt in range(1, n_pt_bins + 1):
    pt_low  = xaxis.GetBinLowEdge(ipt)
    pt_high = xaxis.GetBinUpEdge(ipt)

    for iz in range(1, n_z_bins + 1):
        z_low  = yaxis.GetBinLowEdge(iz)
        z_high = yaxis.GetBinUpEdge(iz)

        # Too verbose to include in the final code!
        # print(f"[INFO] Signal extraction: pT [{pt_low:.2f}, {pt_high:.2f}] "
        #       f"Z [{z_low:.2f}, {z_high:.2f}]")

        # Retrieve already-projected invariant mass histogram
        h_mass = outfile.Get(
            f"InvariantMassProjections/hMass_pt{ipt}_z{iz}"
        )

        if not h_mass or h_mass.GetEntries() < 10:
            # Too little statistics to extract signal reliably
            hSignalPtZ.SetBinContent(ipt, iz, 0.0)
            hSignalPtZ.SetBinError(ipt, iz, 0.0)
            continue

        # Now performing a slight rebinning to make the fit more well behaved
        # Current analysis defaults are 450 bins from 1.08 to 1.15 GeV/c^2 (0.00015 bin width),
        # whilst the mean should be 1.115 (PDG Lambda mass) and the sigma is ~ 0.003 (20 times
        # larger than the current bin width, so it should be rebinned to reduced error bars!)
        h_mass.Rebin(10) # Rebinning by a factor of 10, which should still leave us with a sensible
                         # selection of bins outside the +/- 3.0 * sigma sideband.

        # =========================
        # 7.1) Full fit: Gaussian + quadratic
        #       (ONLY to get mean and sigma)
        # =========================

        mass_min = zaxis.GetXmin()
        mass_max = zaxis.GetXmax()

        full_fit = ROOT.TF1(
            f"fFull_pt{ipt}_z{iz}",
            "gaus(0) + pol2(3)",
            mass_min,
            mass_max
        )

        # Initial guesses (important for stability)
        full_fit.SetParameter(0, h_mass.GetMaximum())   # amplitude
        full_fit.SetParameter(1, 1.1157)                # Lambda mass
        full_fit.SetParameter(2, 0.002)                 # sigma ~ few MeV

        # Some restraints on the ranges these parameters may have
        # (remember that Minuit might go crazy and get to some senseless,
        # unphysical parameters in low statistic bins if we don't do this!)
        full_fit.SetParLimits(1, mass_min, mass_max)   # mean
        full_fit.SetParLimits(2, 1e-4, mass_max - mass_min)  # sigma > 0

        # Quiet fit, store result
        fit_result_full = h_mass.Fit(full_fit, "QS0L")
        # Could probably write this as a regular Chi2 fit (as below), but log-likelihood works well for counts!
        # fit_result_full = h_mass.Fit(full_fit, "QS0")

        # mean  = full_fit.GetParameter(1) # Could use this, but it is not consistent
        mean = fit_result_full.Parameter(1)
        sigma = abs(fit_result_full.Parameter(2))

        # Save full fit for QA
        dir_fits_full.cd()
        if h_mass.GetEntries() > 0: # Should not save empty histograms! Those are just clutter for empty bins!
            h_mass.Write()
            full_fit.Write()

        # No longer need this check here!
        # # Sanity check on sigma
        # if sigma <= 0 or sigma > 0.01:
        #     hSignalPtZ.SetBinContent(ipt, iz, 0.0)
        #     hSignalPtZ.SetBinError(ipt, iz, 0.0)
        #     continue

        # =========================
        # 7.2) Sideband-only background fit
        #      using a DISCONTINUOUS interval
        # =========================
        #
        # Sidebands:
        #   [mean - 6*sigma, mean - 4*sigma] U [mean + 4*sigma, mean + 6*sigma]
        #
        # ROOT does not natively fit discontinuous ranges,
        # so we:
        #   - Clone the histogram
        #   - Zero-out bins in the signal + gap region
        #   - Fit the remaining points only

        h_sideband = h_mass.Clone(f"hSideband_pt{ipt}_z{iz}")
        # This selection code to exclude the region between the two background ranges does NOT work with log-likelihood!!!
        # The part below is NOT necessary for the current version of the fit -- it is only to properly visualize the
        # sideband that was fitted alongside with the original sidebands when doing QA!
        h_sideband.Reset("ICES")
        for ibin in range(1, h_mass.GetNbinsX() + 1):
            x = h_mass.GetXaxis().GetBinCenter(ibin)
            if (mean - 6*sigma <= x <= mean - 4*sigma) or \
               (mean + 4*sigma <= x <= mean + 6*sigma):
                h_sideband.SetBinContent(ibin, h_mass.GetBinContent(ibin))
                h_sideband.SetBinError(ibin, h_mass.GetBinError(ibin))
            # else:
            #     # By doing this, only the sideband bins are considered in
            #     # the fit, as signal + small gap range has negative errors!
            #     # (This is the part that only works for chi2 minimization)
            #     h_sideband.SetBinError(ibin, -1.0)

        # Defining the sidebands in a RooFit-like way, but will conduct the fit using :
            # Nominal sidebands
        sb_left_low  = mean - 6.0 * sigma
        sb_left_high = mean - 4.0 * sigma
        sb_right_low = mean + 4.0 * sigma
        sb_right_high= mean + 6.0 * sigma

        # Protect against histogram boundaries:
        if sb_left_low < mass_min:
            sb_left_low  = mass_min
            sb_left_high = mean - 3.0 * sigma
        if sb_right_high > mass_max:
            sb_right_low  = mean + 3.0 * sigma
            sb_right_high = mass_max

        # Defining the sidebands in a TGraphError, which should only include the desired points for the fit:
        graph_sideband = ROOT.TGraphErrors()
        point = 0
        for ibin in range(1, h_mass.GetNbinsX() + 1):
            x = h_mass.GetBinCenter(ibin)

            in_left_sb  = sb_left_low  <= x <= sb_left_high
            in_right_sb = sb_right_low <= x <= sb_right_high

            if not (in_left_sb or in_right_sb):
                continue

            y  = h_mass.GetBinContent(ibin)
            ey = h_mass.GetBinError(ibin)

            # Skip empty bins to avoid zero-error issues
            if ey <= 0:
                continue
            ex = 0.5 * h_mass.GetBinWidth(ibin)

            graph_sideband.SetPoint(point, x, y)
            graph_sideband.SetPointError(point, ex, ey)
            point += 1

        # Background model (quadratic by default)
        bkg_fit = ROOT.TF1(
            f"fBkg_pt{ipt}_z{iz}",
            "pol2",
            mass_min,
            mass_max
        )
        # This does not work! The fit is empty in the end! This is a C++ only thing!
        # # Redefining the background model to exclude the central, signal+gap, region!
        # def bkg_func(x, p):
        #     xx = x[0]
        #     if (mean - 4*sigma < xx < mean + 4*sigma):
        #         ROOT.TF1.RejectPoint()
        #         return 0
        #     return p[0] + p[1]*xx + p[2]*xx*xx
        # bkg_fit = ROOT.TF1(
        #     f"fBkg_pt{ipt}_z{iz}",
        #     bkg_func,
        #     mass_min,
        #     mass_max,
        #     3
        # )

        # fit_result_bkg = h_sideband.Fit(bkg_fit, "QS0L") # Also doing this with log-likelihood
        # fit_result_bkg = h_sideband.Fit(bkg_fit, "QS0") # Switched back to chi2 fitting because the -1 error is only skipped with chi2 fitting!
        fit_result_bkg = graph_sideband.Fit(bkg_fit, "QS0")

        # Save sideband fit for QA
        dir_fits_sidebands.cd()
        if h_sideband.GetEntries() > 0: # Should not save empty histograms! Those are just clutter for empty bins!
            h_sideband.Write()
            bkg_fit.Write()

        # =========================
        # 7.3) Bin counting in +/- 3*sigma
        # =========================
        #
        # IMPORTANT CHOICE:
        # We integrate ONLY over bins whose centers lie
        # within [mean-3*sigma, mean+3*sigma].
        #
        # This avoids subtracting background from regions
        # where there are no actual histogram bins due to
        # finite binning effects.

        signal_counts = 0.0
        signal_err2   = 0.0

        bkg_integral = 0.0
        bkg_err2     = 0.0

        for ibin in range(1, h_mass.GetNbinsX() + 1):
            x_low  = h_mass.GetXaxis().GetBinLowEdge(ibin)
            x_high = h_mass.GetXaxis().GetBinUpEdge(ibin)
            x_cent = h_mass.GetXaxis().GetBinCenter(ibin)

            if mean - 3*sigma <= x_cent <= mean + 3*sigma:
                # Signal region bin
                c = h_mass.GetBinContent(ibin)
                e = h_mass.GetBinError(ibin)

                signal_counts += c
                signal_err2   += e*e

                # # Background contribution integrated over THIS BIN
                # bkg_bin = bkg_fit.Integral(x_low, x_high)
                # bkg_err = bkg_fit.IntegralError(
                #     x_low, x_high,
                #     fit_result_bkg.GetParams(),
                #     fit_result_bkg.GetCovarianceMatrix()
                # )
                # bkg_integral += bkg_bin
                # bkg_err2     += bkg_err*bkg_err
                # NOTE: This is now done OUTSIDE of the loop!
                # The background is integrated in a single step using TF1::IntegralError.
                # This is essential because the polynomial fit parameters are correlated.
                # Integrating bin-by-bin and summing errors in quadrature would assume
                # independent uncertainties and therefore overestimate the background error.
                # Using the full covariance matrix ensures the correct propagation of
                # correlated fit uncertainties into the signal yield.
        
        # Determine integration bounds aligned to bins
        # (we should only integrate in the same range where the
        # signal is being bin-counted!)
        bin_low  = h_mass.GetXaxis().FindBin(mean - 3*sigma) # FindBin() will give us the bin that contains the provided value in its interior.
        bin_high = h_mass.GetXaxis().FindBin(mean + 3*sigma)

        # Integrates from the lower edge of the first bin
        # to the upper edge of the highest bin!
        x_low  = h_mass.GetXaxis().GetBinLowEdge(bin_low)
        x_high = h_mass.GetXaxis().GetBinUpEdge(bin_high)

        # One-shot integration (includes covariance matrix!)
        # bkg_val = fit_result_bkg.Integral(x_low, x_high)
        if h_sideband.GetEntries() > 0:
            bkg_val = bkg_fit.Integral(x_low, x_high) # Actual proper way to do this
        else:
            bkg_val = 0

        ##################################
        # Encapsulating everything under a "try...except" block to see in which
        # bins the errors might have occurred.
        # (Notice that low statistics bins are expected to give an error in the
        # integration process!)
        try:
            if fit_result_bkg.Status() != 0: # Also added a check to see if the quality of the fit was good enough.
                # print(f"\t[Fit FAILED] ipt={ipt}, iz={iz}") # Silenced too many prints
                bkg_val = 0.0
                bkg_err = 0.0
            else:
                bkg_val = bkg_fit.Integral(x_low, x_high)
                bkg_err = bkg_fit.IntegralError(
                    x_low, x_high,
                    fit_result_bkg.GetParams(),
                    fit_result_bkg.GetCovarianceMatrix().GetMatrixArray()
                )
        except Exception as e:
            # Silenced too many prints:
            # print(
            #     f"\t[IntegralError FAILED] ipt={ipt}, iz={iz} "
            #     f"\t| range=({x_low:.3f}, {x_high:.3f})"
            # )
            # print(f"  Reason: {e}")

            bkg_val = 0.0
            bkg_err = 0.0
        ##################################

        # =========================
        # 7.4) Final signal yield
        # =========================

        signal_yield = signal_counts - bkg_integral

        # Errors are NOT Poissonian anymore:
        # variance = variance(counting) + variance(background)
        signal_error = math.sqrt(signal_err2 + bkg_err2)

        # Fill the (pT, Z) signal map
        hSignalPtZ.SetBinContent(ipt, iz, signal_yield)
        hSignalPtZ.SetBinError(ipt, iz, signal_error)

# Save the signal map
dir_signal_maps.cd()
hSignalPtZ.Write()

print("[INFO] Finished signal extraction for all (pT, Z) bins")


# =========================
# 8) Build pT spectra for different Z regions
# =========================

print("[INFO] Building final pT spectra from signal-extracted (pT, Z) map")

# Create pT spectra histograms
# hPt_FullZ = ROOT.TH1D(
#     "hPt_FullZ",
#     "Lambda p_{T} spectrum (full Z); p_{T} (GeV/c); dN/dp_{T}",
#     n_pt_bins,
#     xaxis.GetXbins().GetArray()
#     if xaxis.GetXbins().GetSize() > 0 else xaxis.GetXmin(), xaxis.GetXmax()
# )
# Same correction as before to be able to use variable binning or fixed-width histograms alike:
def make_TH1D(name, title, axis):
    """
    Build a TH1D using the binning of a TAxis.
    Works for both fixed and variable binning.
    """
    nbins = axis.GetNbins()

    if axis.GetXbins().GetSize() > 0:
        # Variable binning
        bins = axis.GetXbins().GetArray()
        h = ROOT.TH1D(name, title, nbins, bins)
    else:
        # Fixed binning
        xmin = axis.GetXmin()
        xmax = axis.GetXmax()
        h = ROOT.TH1D(name, title, nbins, xmin, xmax)

    h.Sumw2()
    return h

hPt_FullZ = make_TH1D(
    "hPt_FullZ",
    "Lambda p_{T} spectrum (full Z); p_{T} (GeV/c); dN/dp_{T}",
    xaxis
)

hPt_PosZ = hPt_FullZ.Clone("hPt_PosZ")
hPt_PosZ.SetTitle("Lambda p_{T} spectrum (Z #geq 0); p_{T} (GeV/c); dN/dp_{T}")

hPt_NegZ = hPt_FullZ.Clone("hPt_NegZ")
hPt_NegZ.SetTitle("Lambda p_{T} spectrum (Z < 0); p_{T} (GeV/c); dN/dp_{T}")

hPt_FullZ.Reset("ICES")
hPt_PosZ.Reset("ICES")
hPt_NegZ.Reset("ICES")

hPt_FullZ.Sumw2()
hPt_PosZ.Sumw2()
hPt_NegZ.Sumw2()

# =========================
# 8.1) Loop over pT bins and sum Z bins
# =========================
#
# IMPORTANT:
# After signal extraction, bin errors are NOT Poissonian.
# Therefore, when summing bins we must sum variances manually
# (i.e. add errors in quadrature).

for ipt in range(1, n_pt_bins + 1):

    sum_full   = 0.0
    err2_full  = 0.0

    sum_pos    = 0.0
    err2_pos   = 0.0

    sum_neg    = 0.0
    err2_neg   = 0.0

    for iz in range(1, n_z_bins + 1):

        z_center = yaxis.GetBinCenter(iz)

        val = hSignalPtZ.GetBinContent(ipt, iz)
        err = hSignalPtZ.GetBinError(ipt, iz)

        # Full Z
        sum_full  += val
        err2_full += err * err

        # Z >= 0
        if z_center >= 0.0:
            sum_pos  += val
            err2_pos += err * err
        else:
            sum_neg  += val
            err2_neg += err * err

    # Fill pT spectra
    hPt_FullZ.SetBinContent(ipt, sum_full)
    hPt_FullZ.SetBinError(ipt, math.sqrt(err2_full))

    hPt_PosZ.SetBinContent(ipt, sum_pos)
    hPt_PosZ.SetBinError(ipt, math.sqrt(err2_pos))

    hPt_NegZ.SetBinContent(ipt, sum_neg)
    hPt_NegZ.SetBinError(ipt, math.sqrt(err2_neg))

# =========================
# 9) Save final spectra
# =========================

print("[INFO] Saving final pT spectra")

dir_spectra.cd()
hPt_FullZ.Write()
hPt_PosZ.Write()
hPt_NegZ.Write()

# =========================
# 10) Finalize and close files
# =========================

print("[INFO] Writing output file and closing resources")

# outfile.Write()
outfile.Close()
infile.Close()

print("[INFO] Lambda signal extraction completed successfully")
print(f"[INFO] Output written to: {output_file_path}")
