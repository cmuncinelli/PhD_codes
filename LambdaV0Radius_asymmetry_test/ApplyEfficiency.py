#!/usr/bin/env python3
"""
ApplyEfficiency.py

Simplified, Lambda-only version of calculateCorrectedpTSpectra.py.

Inputs:
 - AnalysisResults_root : contains "Lambda/hLambdaPtZMass" (TH3D)
 - EfficiencyCorrections_root : contains hEfficiencyXAcceptance, hSignalLoss, hEventSplitting, hEventLoss
 - ProcessedProjections_root : contains TF1 "hMass_gausPlusQuadratic" used to determine mass mean & sigma

Output:
 - CorrectedpTSpectra_Lambda.root (saved in same folder as AnalysisResults_root)
   Contains:
     - TH2D with Lambda pT vs Z projected in +/- 3sigma mass window
     - Three corrected TH1D pT spectra: full Z, Z >= 0, Z < 0
     - Copies of raw and correction histograms for provenance
"""

import sys, os, array
import ROOT, math
ROOT.gErrorIgnoreLevel = 6001  # silence some ROOT warnings

# --------------------------- Helper functions (preserved/adapted) ---------------------------
def getBinEdges(hInputSpectra):
    """
    Return an array.array('d') with bin edges from a 1D histogram.
    Used to rebin generated histograms to the same binning as raw spectra.
    """
    nBinsX = hInputSpectra.GetNbinsX()
    XAxis = hInputSpectra.GetXaxis()
    pTbins = [XAxis.GetBinLowEdge(i + 1) for i in range(nBinsX + 1)]
    return array.array('d', pTbins)

def getErrorHist(CorrSpectra, RawSpectra, EffAcc):
    """
    Returns three TH1D histograms storing relative errors for corr / raw / eff.
    (Preserves the behaviour of the original helper.)
    """
    nbins = EffAcc.GetNbinsX()
    edges = [EffAcc.GetXaxis().GetBinLowEdge(b) for b in range(1, nbins+1)]
    edges.append(EffAcc.GetXaxis().GetBinUpEdge(nbins))
    edges_array = array.array('d', edges)

    hRawSpectraErrors = ROOT.TH1D("RawpTSpectraErrors", "", nbins, edges_array)
    hCorrSpectraErrors = ROOT.TH1D("CorrpTSpectraErrors", "", nbins, edges_array) 
    hEffAccErrors = ROOT.TH1D("EffAccErrors", "", nbins, edges_array)

    for b in range(1, nbins+1):
        CorrValue = CorrSpectra.GetBinContent(b)
        CorrError = CorrSpectra.GetBinError(b)
        RawValue = RawSpectra.GetBinContent(b)
        RawError = RawSpectra.GetBinError(b)
        EffValue = EffAcc.GetBinContent(b)
        EffError = EffAcc.GetBinError(b)

        hCorrSpectraErrors.SetBinContent(b, CorrError/CorrValue if CorrValue != 0 else 0)
        hCorrSpectraErrors.SetBinError(b, 1e-20)

        hRawSpectraErrors.SetBinContent(b, RawError/RawValue if RawValue != 0 else 0)
        hRawSpectraErrors.SetBinError(b, 1e-20)

        hEffAccErrors.SetBinContent(b, EffError/EffValue if EffValue != 0 else 0)
        hEffAccErrors.SetBinError(b, 1e-20)

    return hCorrSpectraErrors, hRawSpectraErrors, hEffAccErrors

# --------------------------- Core routine --------------------------------
def apply_efficiency(analysis_results_path, corrections_path, processed_projections_path, main_dir="asymmetric_rapidity_test"):
    """
    Main function.
    - analysis_results_path: AnalysisResults root with Lambda/hLambdaPtZMass
    - corrections_path: root with MC corrections (hEfficiencyXAcceptance, hSignalLoss, hEventSplitting, hEventLoss)
    - processed_projections_path: root with TF1 fit named "hMass_gausPlusQuadratic"
    """
    # 1) basic checks & open files
    for p in (analysis_results_path, corrections_path, processed_projections_path):
        if not os.path.isfile(p):
            raise FileNotFoundError(f"Input file not found: {p}")

    fin = ROOT.TFile.Open(analysis_results_path, "READ")
    fcor = ROOT.TFile.Open(corrections_path, "READ")
    fproc = ROOT.TFile.Open(processed_projections_path, "READ")

    if not fin or fin.IsZombie():
        raise RuntimeError(f"Could not open {analysis_results_path}")
    if not fcor or fcor.IsZombie():
        raise RuntimeError(f"Could not open {corrections_path}")
    if not fproc or fproc.IsZombie():
        raise RuntimeError(f"Could not open {processed_projections_path}")

    # 2) retrieve TF1 fit from ProcessedProjections to get mean and sigma
    fit_name = "hMass_gausPlusQuadratic"
    fitFunc = fproc.Get(fit_name)
    if not fitFunc or not isinstance(fitFunc, ROOT.TF1):
        # Try to search for the function in the file (some tasks store TF1s as primitives)
        print(f"Warning: TF1 '{fit_name}' not found directly in {processed_projections_path}. Trying file search...")
        fitFunc = None
        for key in fproc.GetListOfKeys():
            obj = key.ReadObj()
            if isinstance(obj, ROOT.TF1) and obj.GetName() == fit_name:
                fitFunc = obj
                break
        if not fitFunc:
            raise RuntimeError(f"TF1 '{fit_name}' not found in {processed_projections_path}. Needed to determine +/- 3sigma mass window.")

    mean = fitFunc.GetParameter(1)  # gaussian mean typically at par 1 for gaus
    sigma = abs(fitFunc.GetParameter(2))  # sigma typically at par 2
    mass_min = mean - 3.0 * sigma
    mass_max = mean + 3.0 * sigma

    print(f"Using mass window: mean = {mean:.6f}, sigma = {sigma:.6f} -> [{mass_min:.6f}, {mass_max:.6f}]")

    # 3) get raw 3D histogram: Lambda/hLambdaPtZMass
    th3_path = f"{main_dir}/Lambda/hLambdaPtZMass"
    th3 = fin.Get(th3_path)
    if not th3:
        # try without nested Lambda dir
        th3 = fin.Get("Lambda/hLambdaPtZMass") or fin.Get("hLambdaPtZMass")
    if not th3:
        raise RuntimeError(f"TH3D 'Lambda/hLambdaPtZMass' not found in {analysis_results_path} under '{main_dir}' or root. Make sure the path is correct.")

    if not isinstance(th3, ROOT.TH3):
        raise RuntimeError("Found object is not a TH3. Aborting.")

    # 4) determine mass axis bin indices and project TH3 -> TH2 (pT vs Z) integrating over mass window
    mass_axis = th3.GetZaxis()  # per your definition: (axisPt, axisZPos, axisLambdaMass) => mass is axis 3 (z)
    bin_min = mass_axis.FindBin(mass_min)
    bin_max = mass_axis.FindBin(mass_max)
    # clamp bin ranges
    if bin_min < 1: bin_min = 1
    if bin_max > mass_axis.GetNbins(): bin_max = mass_axis.GetNbins()

    print(f"Mass axis bins used: {bin_min} .. {bin_max} of {mass_axis.GetNbins()}")

    # ProjectionXY integrates over z-axis bins (3rd axis) -> gives TH2: x=pt, y=zpos
    # h2_pt_z_masswindow = th3.ProjectionXY("hLambdaPtZ_MassWindow", bin_min, bin_max)
    # h2_pt_z_masswindow.SetTitle("Lambda p_{T} vs Z (mass window +/- 3sigma)")
    # h2_pt_z_masswindow.Sumw2()
    # Had to actually create a ROOT 6.36 version equivalent -- Same problem as in extract_analysis_results.cxx!
    # This version does not have ProjectionXY yet!
        # Set the Z-axis to the desired mass window
    th3.GetZaxis().SetRange(bin_min, bin_max)

    # Equivalent to ProjectionXY: integrates over z, keeps X=pt, Y=zpos
    h2_pt_z_masswindow = th3.Project3D("yx")
    h2_pt_z_masswindow.SetName("hLambdaPtZ_MassWindow")
    h2_pt_z_masswindow.SetTitle("Lambda p_{T} vs Z (mass window +/- 3sigma)")
    # Restore full range of Z
    th3.GetZaxis().SetRange(0, 0)

    h2_pt_z_masswindow.Sumw2()


    # ------------------------------
    # 5) produce 1D pT spectra (projection on X) for:
    #   - full Z (all y bins)
    #   - Z >= 0  (y bins with center >= 0)
    #   - Z < 0   (y bins with center < 0)
    # ------------------------------
    yaxis = h2_pt_z_masswindow.GetYaxis()
    nbins_y = yaxis.GetNbins()
    # find first bin with center >= 0
    first_nonneg_bin = None
    for by in range(1, nbins_y + 1):
        center = yaxis.GetBinCenter(by)
        if center >= 0:
            first_nonneg_bin = by
            break
    if first_nonneg_bin is None:
        first_nonneg_bin = nbins_y + 1  # no non-negative bins

    # full Z projection
    h_pt_fullZ = h2_pt_z_masswindow.ProjectionX("hRaw_pT_fullZ", 1, nbins_y)
    h_pt_fullZ.SetTitle("Raw Lambda p_{T} (full Z) within mass window")
    h_pt_fullZ.Sumw2()

    # Z >= 0 projection
    if first_nonneg_bin <= nbins_y:
        h_pt_Zpos = h2_pt_z_masswindow.ProjectionX("hRaw_pT_Zpos", first_nonneg_bin, nbins_y)
    else:
        h_pt_Zpos = h_pt_fullZ.Clone("hRaw_pT_Zpos")
        for i in range(1, h_pt_Zpos.GetNbinsX() + 1):
            h_pt_Zpos.SetBinContent(i, 0.0); h_pt_Zpos.SetBinError(i, 0.0)
    h_pt_Zpos.SetTitle("Raw Lambda p_{T} (Z >= 0) within mass window")
    h_pt_Zpos.Sumw2()

    # Z < 0 projection
    if first_nonneg_bin > 1:
        h_pt_Zneg = h2_pt_z_masswindow.ProjectionX("hRaw_pT_Zneg", 1, first_nonneg_bin - 1)
    else:
        h_pt_Zneg = h_pt_fullZ.Clone("hRaw_pT_Zneg")
        for i in range(1, h_pt_Zneg.GetNbinsX() + 1):
            h_pt_Zneg.SetBinContent(i, 0.0); h_pt_Zneg.SetBinError(i, 0.0)
    h_pt_Zneg.SetTitle("Raw Lambda p_{T} (Z < 0) within mass window")
    h_pt_Zneg.Sumw2()

    # ------------------------------
    # 6) load corrections from corrections file (TH2D: pT vs Z)
    # ------------------------------
    # Load 2D correction histograms written by CalculateEfficiency.py
    hEffAcc = fcor.Get("hEfficiencyXAcceptance2D") # Efficiency × Acceptance (TH2D pT×Z)
    hSignalLoss = fcor.Get("hSignalLoss2D")     # Signal Loss (TH2D pT×Z)

    # Load scalar event-level corrections
    hEvtSplitting = fcor.Get("hEventLoss")
    hEvtLoss      = fcor.Get("hEventSplitting")


    # Basic checks
    if not hEffAcc or not hSignalLoss or not hEvtSplitting or not hEvtLoss:
        raise RuntimeError("One or more correction histograms not found in corrections file. Expected: hEffAcc, hSignalLoss, hEventSplitting, hEventLoss")

    # Ensure Sumw2 on TH2D corr histograms
    hEffAcc.Sumw2()
    hSignalLoss.Sumw2()
    hEvtSplitting.Sumw2()
    hEvtLoss.Sumw2()

    # enable Sumw2 on raw spectra
    h_pt_fullZ.Sumw2()
    h_pt_Zpos.Sumw2()
    h_pt_Zneg.Sumw2()

    # get event-level scalars
    EvtSplitValue = hEvtSplitting.GetBinContent(1)
    EvtSplitError = hEvtSplitting.GetBinError(1)
    EvtLossValue = hEvtLoss.GetBinContent(1)
    EvtLossError = hEvtLoss.GetBinError(1)

    if EvtSplitValue == 0:
        raise RuntimeError("Event splitting value is zero -> division by zero in scaling.")

    print("EvtSplitValue: ", EvtSplitValue, " EvtLossValue: ", EvtLossValue)

    # ------------------------------
    # 7) rebin/align correction histograms to raw pT binning (in X), keep Z bins untouched
    #    We'll produce 1D correction histograms for each Z-region by projecting the 2D corrections.
    # ------------------------------

    # get raw pT bin edges (for variable binning)
    raw_edges = getBinEdges(h_pt_fullZ)
    # array.array('d') of edges, as required by RebinX
    edges_array = raw_edges

    # Rebin the TH2D correction histograms in X (pT) to the raw pT binning if necessary.
    # Use RebinX to preserve Z axis.
    # def rebin2d_in_x_to_edges(h2_in, edges_array, new_name):
    #     # if same nbins in X already, clone
    #     if h2_in.GetNbinsX() == len(edges_array) - 1:
    #         return h2_in.Clone(new_name)
    #     # RebinX returns a new histogram with specified x-binning, keeping Y unchanged
    #     reb = h2_in.RebinX(len(edges_array) - 1, new_name, edges_array)
    #     return reb
    # Had to reimplement as RebinX was a method for TH1Ds!!!
    def rebin2d_in_x_to_edges(h2_in, edges_array, new_name):
        """
        Rebin a TH2D in X using variable bin edges.
        Keeps Y-axis binning unchanged.
        ROOT lacks TH2::RebinX with variable edges, so we do it manually.
        """

        old_nx = h2_in.GetNbinsX()
        old_ny = h2_in.GetNbinsY()

        new_nx = len(edges_array) - 1

        # If no rebinning needed:
        if old_nx == new_nx:
            return h2_in.Clone(new_name)

        # Get original Y binning (variable or fixed)
        yaxis = h2_in.GetYaxis()
        y_edges = yaxis.GetXbins().GetArray()
        if not y_edges:
            y_edges = array.array('d', [yaxis.GetBinLowEdge(i+1) for i in range(old_ny+1)])

        # Create new TH2D with desired X binning and original Y binning
        h2_out = ROOT.TH2D(
            new_name,
            h2_in.GetTitle(),
            new_nx, edges_array,
            old_ny, y_edges
        )
        h2_out.Sumw2()

        # Loop over all original bins and transfer their contents
        for ix in range(1, old_nx + 1):
            x_center = h2_in.GetXaxis().GetBinCenter(ix)

            for iy in range(1, old_ny + 1):
                y_center = h2_in.GetYaxis().GetBinCenter(iy)

                content = h2_in.GetBinContent(ix, iy)
                error   = h2_in.GetBinError(ix, iy)

                if content != 0 or error != 0:
                    # Add content to new bin
                    h2_out.Fill(x_center, y_center, content)

                    # ROOT's Fill() does NOT preserve errors, so set error manually:
                    new_ix = h2_out.GetXaxis().FindBin(x_center)
                    new_iy = h2_out.GetYaxis().FindBin(y_center)
                    current_err = h2_out.GetBinError(new_ix, new_iy)

                    # Add errors in quadrature
                    h2_out.SetBinError(new_ix, new_iy, (current_err**2 + error**2)**0.5)

        return h2_out


    hEffAcc_reb2d = rebin2d_in_x_to_edges(hEffAcc, edges_array, "hEffAcc_reb2d")
    hSignalLoss_reb2d = rebin2d_in_x_to_edges(hSignalLoss, edges_array, "hSignalLoss_reb2d")

    hEffAcc_reb2d.Sumw2()
    hSignalLoss_reb2d.Sumw2()

    # ------------------------------
    # Create 1D correction histograms for each Z-region by projecting the TH2D along Z
    # ------------------------------
    nbins_y_eff = hEffAcc_reb2d.GetNbinsY()

    # helper to project TH2D->TH1D over a range of Y (Z) bins and name/label properly
    def project_corr_2d_to_1d(h2_corr, y_bin_min, y_bin_max, out_name, out_title):
        # ProjectX takes (name, y_bin_min, y_bin_max)
        h1 = h2_corr.ProjectionX(out_name, y_bin_min, y_bin_max)
        h1.SetTitle(out_title)
        h1.Sumw2()
        return h1

    # full Z: project all Y bins
    hEff_fullZ = project_corr_2d_to_1d(hEffAcc_reb2d, 1, nbins_y_eff, "hEff_fullZ", "EffAcc (full Z)")
    hSig_fullZ = project_corr_2d_to_1d(hSignalLoss_reb2d, 1, nbins_y_eff, "hSig_fullZ", "SignalLoss (full Z)")

    # Z >= 0:
    if first_nonneg_bin <= nbins_y_eff:
        hEff_Zpos = project_corr_2d_to_1d(hEffAcc_reb2d, first_nonneg_bin, nbins_y_eff, "hEff_Zpos", "EffAcc (Z>=0)")
        hSig_Zpos = project_corr_2d_to_1d(hSignalLoss_reb2d, first_nonneg_bin, nbins_y_eff, "hSig_Zpos", "SignalLoss (Z>=0)")
    else:
        # create zeroed clones with same binning
        hEff_Zpos = hEff_fullZ.Clone("hEff_Zpos")
        hSig_Zpos = hSig_fullZ.Clone("hSig_Zpos")
        for ib in range(1, hEff_Zpos.GetNbinsX() + 1):
            hEff_Zpos.SetBinContent(ib, 0.0); hEff_Zpos.SetBinError(ib, 0.0)
            hSig_Zpos.SetBinContent(ib, 0.0); hSig_Zpos.SetBinError(ib, 0.0)

    # Z < 0:
    if first_nonneg_bin > 1:
        hEff_Zneg = project_corr_2d_to_1d(hEffAcc_reb2d, 1, first_nonneg_bin - 1, "hEff_Zneg", "EffAcc (Z<0)")
        hSig_Zneg = project_corr_2d_to_1d(hSignalLoss_reb2d, 1, first_nonneg_bin - 1, "hSig_Zneg", "SignalLoss (Z<0)")
    else:
        hEff_Zneg = hEff_fullZ.Clone("hEff_Zneg")
        hSig_Zneg = hSig_fullZ.Clone("hSig_Zneg")
        for ib in range(1, hEff_Zneg.GetNbinsX() + 1):
            hEff_Zneg.SetBinContent(ib, 0.0); hEff_Zneg.SetBinError(ib, 0.0)
            hSig_Zneg.SetBinContent(ib, 0.0); hSig_Zneg.SetBinError(ib, 0.0)

    # ------------------------------
    # 8) Calculate corrected spectra per region:
    # corrected = raw / EffAcc / SignalLoss * (EvtLoss / EvtSplit)
    # ------------------------------
    def make_corrected(raw_hist, eff1d, sig1d, label):
        raw_clone = raw_hist.Clone(f"hRawClone_{label}")
        raw_clone.Sumw2()
        corr = raw_clone.Clone(f"hCorrected_{label}")
        # protect against zero bins in eff/sig histograms
        # We'll perform a safe bin-wise division: corr_bin = raw_bin / (eff_bin * sig_bin)
        # ROOT's Divide already handles bin-by-bin division if bins match.
        # Ensure bin counts match:
        if corr.GetNbinsX() != eff1d.GetNbinsX():
            raise RuntimeError(f"Bin mismatch between raw and correction histograms for {label}")
        corr.Divide(eff1d)      # divide by EffAcc
        corr.Divide(sig1d)      # divide by SignalLoss
        corr.Scale(EvtLossValue / EvtSplitValue)  # Scale by event-level factor
        corr.SetTitle(f"Corrected Lambda pT ({label})")
        return corr

    hCorr_fullZ = make_corrected(h_pt_fullZ, hEff_fullZ, hSig_fullZ, "fullZ")
    hCorr_Zpos = make_corrected(h_pt_Zpos, hEff_Zpos, hSig_Zpos, "Zpos")
    hCorr_Zneg = make_corrected(h_pt_Zneg, hEff_Zneg, hSig_Zneg, "Zneg")

    # ------------------------------
    # 9) Create relative-error histograms (preserve original helper behavior)
    # ------------------------------
    hCorrErrors, hRawErrors, hEffErrors = getErrorHist(hCorr_fullZ, h_pt_fullZ, hEff_fullZ)



    # =====================================================================
    # 9.5) Diagnostic plots: efficiencies, signal-loss, raw and corrected spectra
    # =====================================================================

    # We will create all canvases and write them inside a "Diagnostics" dir in the output ROOT file.
    diag_canvases = []

    def make_ratio_hist(h_num, h_den, name):
        """Create ratio histogram with error propagation."""
        h_ratio = h_num.Clone(name)
        h_ratio.Reset("ICESM")  # clear but preserve binning

        for ib in range(1, h_ratio.GetNbinsX() + 1):
            A = h_num.GetBinContent(ib)
            B = h_den.GetBinContent(ib)
            eA = h_num.GetBinError(ib)
            eB = h_den.GetBinError(ib)

            if B > 0:
                R = A / B
                # propagate uncertainties
                err = R * math.sqrt((eA / A)**2 + (eB / B)**2) if A > 0 else 0
            else:
                R = 0
                err = 0

            h_ratio.SetBinContent(ib, R)
            h_ratio.SetBinError(ib, err)
        return h_ratio

    def normalize_spectrum(h, Nev, label):
        """Convert histogram to (1/Nev)*(1/ΔpT)*counts."""
        out = h.Clone(label)
        for ib in range(1, out.GetNbinsX() + 1):
            bw = out.GetBinWidth(ib)
            val = out.GetBinContent(ib)
            err = out.GetBinError(ib)
            out.SetBinContent(ib, val / (Nev * bw))
            out.SetBinError(ib, err / (Nev * bw))
        return out

    # Event normalization = total reconstructed events (same as in corrections)
    NEV = hEvtSplitting.GetBinContent(1) * hEvtLoss.GetBinContent(1)
    if NEV <= 0:
        raise RuntimeError("Got a null number of events!")

    # ------------------------------
    # PLOT 1: Signal Loss vs pT
    # ------------------------------

    c_sig = ROOT.TCanvas("cSignalLoss", "Signal Loss vs pT", 1200, 900)
    c_sig.SetMargin(0.12, 0.05, 0.10, 0.08)

    # format histos
    for h, col, name in [
        (hSig_fullZ, ROOT.kBlack, "Full Z"),
        (hSig_Zpos,  ROOT.kRed+1, "Z >= 0"),
        (hSig_Zneg,  ROOT.kBlue+1, "Z < 0")
    ]:
        h.SetLineColor(col)
        h.SetMarkerColor(col)
        h.SetMarkerStyle(20)
        h.SetTitle(";p_{T} (GeV/c);Signal Loss")
        h.GetYaxis().SetTitleOffset(1.3)

    # draw in order
    hSig_fullZ.Draw("E")
    hSig_Zpos.Draw("E SAME")
    hSig_Zneg.Draw("E SAME")

    leg = ROOT.TLegend(0.60, 0.70, 0.88, 0.88)
    leg.AddEntry(hSig_fullZ, "#Lambda full Z", "lep")
    leg.AddEntry(hSig_Zpos,  "#Lambda Z >= 0", "lep")
    leg.AddEntry(hSig_Zneg,  "#Lambda Z < 0", "lep")
    leg.Draw()

    text = ROOT.TLatex()
    text.SetNDC(True)
    text.SetTextSize(0.04)
    text.DrawLatex(0.18, 0.86, "pp, 13.6 TeV, #Lambda")

    diag_canvases.append(c_sig)

    # ------------------------------
    # PLOT 2: Efficiency * Acceptance
    # ------------------------------

    c_eff = ROOT.TCanvas("cEffAcc", "Efficiency vs pT", 1200, 900)
    c_eff.SetMargin(0.12, 0.05, 0.10, 0.08)

    for h, col, name in [
        (hEff_fullZ, ROOT.kBlack, "Full Z"),
        (hEff_Zpos,  ROOT.kRed+1, "Z >= 0"),
        (hEff_Zneg,  ROOT.kBlue+1, "Z < 0")
    ]:
        h.SetLineColor(col)
        h.SetMarkerColor(col)
        h.SetMarkerStyle(20)
        h.SetTitle(";p_{T} (GeV/c);Efficiency * Acceptance")
        h.GetYaxis().SetTitleOffset(1.3)

    hEff_fullZ.Draw("E")
    hEff_Zpos.Draw("E SAME")
    hEff_Zneg.Draw("E SAME")

    leg2 = ROOT.TLegend(0.60, 0.70, 0.88, 0.88)
    leg2.AddEntry(hEff_fullZ, "#Lambda full Z", "lep")
    leg2.AddEntry(hEff_Zpos,  "#Lambda Z >= 0", "lep")
    leg2.AddEntry(hEff_Zneg,  "#Lambda Z < 0", "lep")
    leg2.Draw()

    text.DrawLatex(0.18, 0.86, "pp, 13.6 TeV, #Lambda")

    diag_canvases.append(c_eff)

    # ------------------------------
    # PLOT 3: Raw spectra + Ratio subplot
    # ------------------------------

    # Normalize the raw spectra
    hRaw_fullZ_norm = normalize_spectrum(h_pt_fullZ, NEV, "hRaw_fullZ_norm")
    hRaw_Zpos_norm  = normalize_spectrum(h_pt_Zpos,  NEV, "hRaw_Zpos_norm")
    hRaw_Zneg_norm  = normalize_spectrum(h_pt_Zneg,  NEV, "hRaw_Zneg_norm")

    # Ratio histos
    hRaw_ratio_Zpos = make_ratio_hist(hRaw_Zpos_norm, hRaw_fullZ_norm, "hRawRatio_ZposFull")
    hRaw_ratio_Zneg = make_ratio_hist(hRaw_Zneg_norm, hRaw_fullZ_norm, "hRawRatio_ZnegFull")

    c_raw = ROOT.TCanvas("cRawSpectra", "Raw pT spectra", 1200, 1000)
    # Larger bottom panel (35%)
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.35, 1, 1.0)
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.00, 1, 0.35)

    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.12)
    pad1.SetRightMargin(0.05)

    pad2.SetTopMargin(0.03)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(0.12)
    pad2.SetRightMargin(0.05)

    pad1.Draw()
    pad2.Draw()

    # Upper panel
    pad1.cd()
    hRaw_fullZ_norm.GetXaxis().SetRangeUser(0, 20)
    hRaw_Zpos_norm.GetXaxis().SetRangeUser(0, 20)
    hRaw_Zneg_norm.GetXaxis().SetRangeUser(0, 20)

    for h, col in [
        (hRaw_fullZ_norm, ROOT.kBlack),
        (hRaw_Zpos_norm,  ROOT.kRed+1),
        (hRaw_Zneg_norm,  ROOT.kBlue+1)
    ]:
        h.SetLineColor(col)
        h.SetMarkerColor(col)
        h.SetMarkerStyle(20)
        h.SetTitle(";p_{T} (GeV/c);(1/N_{evt}) dN/dp_{T}")
        h.GetYaxis().SetTitleOffset(1.3)

    hRaw_fullZ_norm.Draw("E")
    hRaw_Zpos_norm.Draw("E SAME")
    hRaw_Zneg_norm.Draw("E SAME")

    leg3 = ROOT.TLegend(0.60, 0.70, 0.88, 0.88)
    leg3.AddEntry(hRaw_fullZ_norm, "#Lambda full Z", "lep")
    leg3.AddEntry(hRaw_Zpos_norm,  "#Lambda Z >= 0", "lep")
    leg3.AddEntry(hRaw_Zneg_norm,  "#Lambda Z < 0", "lep")
    leg3.Draw()

    text.DrawLatex(0.18, 0.86, "pp, 13.6 TeV, #Lambda")

    # Lower ratio subplot
    pad2.cd()
    hRaw_ratio_Zpos.GetXaxis().SetRangeUser(0, 20)
    hRaw_ratio_Zneg.GetXaxis().SetRangeUser(0, 20)
    hRaw_ratio_Zpos.SetLineColor(ROOT.kRed+1)
    hRaw_ratio_Zpos.SetMarkerColor(ROOT.kRed+1)
    hRaw_ratio_Zpos.SetMarkerStyle(20)
    hRaw_ratio_Zpos.SetTitle(";p_{T} (GeV/c);Ratio to full Z")

    hRaw_ratio_Zneg.SetLineColor(ROOT.kBlue+1)
    hRaw_ratio_Zneg.SetMarkerColor(ROOT.kBlue+1)
    hRaw_ratio_Zneg.SetMarkerStyle(20)

    hRaw_ratio_Zpos.Draw("E")
    hRaw_ratio_Zneg.Draw("E SAME")

    diag_canvases.append(c_raw)

    # ------------------------------
    # PLOT 4: Corrected spectra + Ratios
    # ------------------------------

    hCorr_fullZ_norm = normalize_spectrum(hCorr_fullZ, NEV, "hCorr_fullZ_norm")
    hCorr_Zpos_norm  = normalize_spectrum(hCorr_Zpos, NEV, "hCorr_Zpos_norm")
    hCorr_Zneg_norm  = normalize_spectrum(hCorr_Zneg, NEV, "hCorr_Zneg_norm")

    hCorr_ratio_Zpos = make_ratio_hist(hCorr_Zpos_norm, hCorr_fullZ_norm, "hCorrRatio_ZposFull")
    hCorr_ratio_Zneg = make_ratio_hist(hCorr_Zneg_norm, hCorr_fullZ_norm, "hCorrRatio_ZnegFull")

    c_corr = ROOT.TCanvas("cCorrectedSpectra", "Corrected pT spectra", 1200, 1000)
    pad1c = ROOT.TPad("pad1c", "pad1c", 0, 0.35, 1, 1.0)
    pad2c = ROOT.TPad("pad2c", "pad2c", 0, 0.00, 1, 0.35)

    pad1c.SetBottomMargin(0.02)
    pad1c.SetLeftMargin(0.12)
    pad1c.SetRightMargin(0.05)

    pad2c.SetTopMargin(0.03)
    pad2c.SetBottomMargin(0.35)
    pad2c.SetLeftMargin(0.12)
    pad2c.SetRightMargin(0.05)

    pad1c.Draw()
    pad2c.Draw()

    # Upper panel
    pad1c.cd()
    hCorr_fullZ_norm.GetXaxis().SetRangeUser(0, 20)
    hCorr_Zpos_norm.GetXaxis().SetRangeUser(0, 20)
    hCorr_Zneg_norm.GetXaxis().SetRangeUser(0, 20)
    for h, col in [
        (hCorr_fullZ_norm, ROOT.kBlack),
        (hCorr_Zpos_norm,  ROOT.kRed+1),
        (hCorr_Zneg_norm,  ROOT.kBlue+1)
    ]:
        h.SetLineColor(col)
        h.SetMarkerColor(col)
        h.SetMarkerStyle(20)
        h.SetTitle(";p_{T} (GeV/c);(1/N_{evt}) dN/dp_{T}")
        h.GetYaxis().SetTitleOffset(1.3)

    hCorr_fullZ_norm.Draw("E")
    hCorr_Zpos_norm.Draw("E SAME")
    hCorr_Zneg_norm.Draw("E SAME")

    leg4 = ROOT.TLegend(0.60, 0.70, 0.88, 0.88)
    leg4.AddEntry(hCorr_fullZ_norm, "#Lambda full Z", "lep")
    leg4.AddEntry(hCorr_Zpos_norm,  "#Lambda Z >= 0", "lep")
    leg4.AddEntry(hCorr_Zneg_norm,  "#Lambda Z < 0", "lep")
    leg4.Draw()

    text.DrawLatex(0.18, 0.86, "pp, 13.6 TeV, #Lambda")

    # Lower panel
    pad2c.cd()
    hCorr_ratio_Zpos.GetXaxis().SetRangeUser(0, 20)
    hCorr_ratio_Zneg.GetXaxis().SetRangeUser(0, 20)
    hCorr_ratio_Zpos.SetLineColor(ROOT.kRed+1)
    hCorr_ratio_Zpos.SetMarkerColor(ROOT.kRed+1)
    hCorr_ratio_Zpos.SetMarkerStyle(20)
    hCorr_ratio_Zpos.SetTitle(";p_{T} (GeV/c);Ratio to full Z")

    hCorr_ratio_Zneg.SetLineColor(ROOT.kBlue+1)
    hCorr_ratio_Zneg.SetMarkerColor(ROOT.kBlue+1)
    hCorr_ratio_Zneg.SetMarkerStyle(20)

    hCorr_ratio_Zpos.Draw("E")
    hCorr_ratio_Zneg.Draw("E SAME")

    diag_canvases.append(c_corr)



    # ------------------------------
    # 10) prepare output file (in same folder as analysis_results_path)
    # ------------------------------
    input_dir = os.path.dirname(os.path.abspath(analysis_results_path))
    output_name = os.path.join(input_dir, "CorrectedpTSpectra_Lambda.root")
    fout = ROOT.TFile.Open(output_name, "RECREATE")

    # Create subdirectories
    dir_main = fout.mkdir("ProcessedSpectra")
    dir_corr = fout.mkdir("CorrectionsUsed")
    dir_diag = fout.mkdir("Diagnostics")

    # Write spectra and processed items
    dir_main.cd()
    h2_pt_z_masswindow.Write()
    h_pt_fullZ.Write()
    h_pt_Zpos.Write()
    h_pt_Zneg.Write()
    hCorr_fullZ.Write("hCorrected_pT_fullZ")
    hCorr_Zpos.Write("hCorrected_pT_Zpos")
    hCorr_Zneg.Write("hCorrected_pT_Zneg")
    hRawErrors.Write()
    hCorrErrors.Write()
    hEffErrors.Write()

    # Write correction inputs
    dir_corr.cd()
    hEffAcc_reb2d.Write()
    hSignalLoss_reb2d.Write()
    hEff_fullZ.Write()
    hEff_Zpos.Write()
    hEff_Zneg.Write()
    hSig_fullZ.Write()
    hSig_Zpos.Write()
    hSig_Zneg.Write()
    hEvtSplitting.Write()
    hEvtLoss.Write()

    # Write Diagnostic canvases
    dir_diag.cd()
    for c in diag_canvases:
        c.Write()

    fout.Close()
    fin.Close()
    fcor.Close()
    fproc.Close()

    print(f"Done. Output saved to: {output_name}")


# --------------------------- CLI --------------------------------
def usage_and_exit():
    print("Usage: python ApplyEfficiency.py <AnalysisResults.root> <EfficiencyCorrections.root> <ProcessedProjections.root>")
    sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        usage_and_exit()
    analysis_results_path = sys.argv[1]
    corrections_path = sys.argv[2]
    processed_projections_path = sys.argv[3]
    apply_efficiency(analysis_results_path, corrections_path, processed_projections_path)
