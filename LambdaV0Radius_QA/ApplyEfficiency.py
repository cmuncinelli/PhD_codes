#!/usr/bin/env python3
"""
ApplyEfficiency.py

Simplified, Lambda-only version of calculateCorrectedpTSpectra.py.

Inputs:
 - AnalysisResults_root : contains "Lambda/hLambdaPtZMass" (TH3D)
 - EfficiencyCorrections_root : contains hEfficiencyXAcceptance, hSignalLoss, hEventSplitting, hEventLoss
#  - ProcessedProjections_root : contains TF1 "hMass_gausPlusQuadratic" used to determine mass mean & sigma --> Deprecated!
 - LambdaSignalExtraction_root : contains the signal-extracted pT spectra.

Output:
 - CorrectedpTSpectra_Lambda.root (saved in same folder as AnalysisResults_root)
   Contains:
     - TH2D with Lambda pT vs Z projected in +/- 3sigma mass window
     - Three corrected TH1D pT spectra: full Z, Z >= 0, Z < 0
     - Copies of raw and correction histograms for provenance
"""

import sys, os, re, array
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
def apply_efficiency(analysis_results_path, corrections_path, main_dir="asymmetric_rapidity_test"):
    """
    Main function.
    - analysis_results_path: AnalysisResults root with Lambda/hLambdaPtZMass
    - corrections_path: root with MC corrections (hEfficiencyXAcceptance, hSignalLoss, hEventSplitting, hEventLoss)
    """
    # - processed_projections_path: root with TF1 fit named "hMass_gausPlusQuadratic" --> NO LONGER NEED THIS. Now using Signal Extraction codes!
    # 1) basic checks & open files
    print("\nStarting code -- Basic checks and opening files")
    for p in (analysis_results_path, corrections_path):
        if not os.path.isfile(p):
            raise FileNotFoundError(f"Input file not found: {p}")

    fin = ROOT.TFile.Open(analysis_results_path, "READ")
    fcor = ROOT.TFile.Open(corrections_path, "READ")
    # fproc = ROOT.TFile.Open(processed_projections_path, "READ")

    # Extracting the file path from the current directory of the AnalysisResults -- Makes things a lot simpler!
        # Can't do the same about the MC efficiency calculation, as the
        # names of the folders of anchored datasets are never obvious!
    input_dir = os.path.dirname(os.path.abspath(analysis_results_path))
    
    running_OO = False
    if "OO" in input_dir or "LF_LHC25h3_pass2_Strangeness" in input_dir:
        print("\t\tHandling OO collisions from run 3.")
        running_OO = True
    running_MC = False
    if "LF_LHC" in input_dir:
        print("\t\tRunning in MC -- Closure tests")
        running_MC = True
    
        # Extract the base name without extension
    base_input_name = os.path.splitext(os.path.basename(analysis_results_path))[0]
        # Remove leading "AnalysisResults-" if present
    clean_name = base_input_name.replace("AnalysisResults-", "", 1)
    signal_extracted_root_path = f"LambdaSignalExtraction-{clean_name}.root"
    signal_extracted_root_path = os.path.join(input_dir, signal_extracted_root_path) # This is in the same folder as the Data/MC analysis results!
    fsigExtractedData = ROOT.TFile.Open(signal_extracted_root_path, "READ")

    if not fin or fin.IsZombie():
        raise RuntimeError(f"Could not open {analysis_results_path}")
    if not fcor or fcor.IsZombie():
        raise RuntimeError(f"Could not open {corrections_path}")
    if not fsigExtractedData or fsigExtractedData.IsZombie():
        raise RuntimeError(f"Could not open {signal_extracted_root_path}")

    # ------------------------------
    # 5) Fetch 1D pT spectra (projection on X) for:
    #   - full Z (all y bins)
    #   - Z >= 0  (y bins with center >= 0)
    #   - Z < 0   (y bins with center < 0)
    # ------------------------------
    h_pt_fullZ_2D = fsigExtractedData.Get("SignalExtractionMaps/hSignalPtZ")
    h_pt_fullZ = fsigExtractedData.Get("FinalPtSpectra/hPt_FullZ") # full Z projection
    h_pt_Zpos = fsigExtractedData.Get("FinalPtSpectra/hPt_PosZ") # Z >= 0 projection
    h_pt_Zneg = fsigExtractedData.Get("FinalPtSpectra/hPt_NegZ") # Z < 0 projection

    # ------------------------------
    # 6) load corrections from corrections file (TH2D: pT vs Z)
    # ------------------------------
    # Load 2D correction histograms written by CalculateEfficiency.py
    hEffAcc = fcor.Get("hEfficiencyXAcceptance2D") # Efficiency × Acceptance (TH2D pT×Z)
    hSignalLoss = fcor.Get("hSignalLoss2D")     # Signal Loss (TH2D pT×Z)

    # Load scalar event-level corrections
    hEvtSplitting = fcor.Get("hEventSplitting")
    hEvtLoss      = fcor.Get("hEventLoss")

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
    # EvtSplitError = hEvtSplitting.GetBinError(1)
    EvtLossValue = hEvtLoss.GetBinContent(1)
    # EvtLossError = hEvtLoss.GetBinError(1)

    if EvtSplitValue == 0:
        raise RuntimeError("Event splitting value is zero -> division by zero in scaling.")

    print("EvtSplitValue: ", EvtSplitValue, " EvtLossValue: ", EvtLossValue)

    # ------------------------------
    # 7) rebin/align correction histograms to raw pT binning (in X), keep Z bins untouched
    #    We'll produce 1D correction histograms for each Z-region by projecting the 2D corrections.
    # ------------------------------

    # get raw pT bin edges (for variable binning)
    # raw_edges = getBinEdges(h_pt_fullZ)
    # edges_array = raw_edges

    # --- Load 1D Corrections directly ---
    # Better than calculating this here, we can just use the result directly!
    # These names must match what you wrote in CalculateEfficiency.py
    hEff1D_Full_in = fcor.Get("hEfficiency1D_Full")
    hEff1D_Pos_in  = fcor.Get("hEfficiency1D_Pos")
    hEff1D_Neg_in  = fcor.Get("hEfficiency1D_Neg")
    
    hSig1D_Full_in = fcor.Get("hSignalLoss1D_Full")
    hSig1D_Pos_in  = fcor.Get("hSignalLoss1D_Pos")
    hSig1D_Neg_in  = fcor.Get("hSignalLoss1D_Neg")

    # Check
    if not hEff1D_Full_in or not hSig1D_Full_in:
        raise RuntimeError("1D Correction histograms not found in file!")

    # # Reworked this function to avoid the N_entries not being set appropriately!
    # def rebin2d_in_x_to_edges(h2_in, edges_array, new_name):
    #     """
    #     Rebin a TH2D in X using variable bin edges.
    #     accumulates Content and Error manually to preserve N_eff.
    #     """
    #     old_nx = h2_in.GetNbinsX()
    #     old_ny = h2_in.GetNbinsY()
    #     new_nx = len(edges_array) - 1

    #     # If no rebinning needed (and simple clone is desired):
    #     if old_nx == new_nx:
    #         # Note: If edges don't match exactly, this might be risky, 
    #         # but usually implies same binning.
    #         return h2_in.Clone(new_name)

    #     yaxis = h2_in.GetYaxis()
    #     y_edges = yaxis.GetXbins().GetArray()
    #     if not y_edges:
    #         y_edges = array.array('d', [yaxis.GetBinLowEdge(i+1) for i in range(old_ny+1)])

    #     h2_out = ROOT.TH2D(new_name, h2_in.GetTitle(), new_nx, edges_array, old_ny, y_edges)
    #     h2_out.Sumw2()

    #     for ix in range(1, old_nx + 1):
    #         x_center = h2_in.GetXaxis().GetBinCenter(ix)
    #         for iy in range(1, old_ny + 1):
    #             content = h2_in.GetBinContent(ix, iy)
    #             error   = h2_in.GetBinError(ix, iy)
                
    #             if content != 0 or error != 0:
    #                 y_center = h2_in.GetYaxis().GetBinCenter(iy)
                    
    #                 # Find new bin
    #                 new_ix = h2_out.GetXaxis().FindBin(x_center)
    #                 new_iy = h2_out.GetYaxis().FindBin(y_center)
                    
    #                 # Manual Accumulation (Fixes "Fill" weight bug)
    #                 cur_val = h2_out.GetBinContent(new_ix, new_iy)
    #                 cur_err = h2_out.GetBinError(new_ix, new_iy)
                    
    #                 h2_out.SetBinContent(new_ix, new_iy, cur_val + content)
    #                 h2_out.SetBinError(new_ix, new_iy, math.sqrt(cur_err**2 + error**2))

    #     return h2_out

    # def rebin1d_to_edges(h1_in, edges_array, new_name):
    #     """
    #     Rebin a TH1D using variable bin edges with manual error accumulation.
    #     """
    #     old_nx = h1_in.GetNbinsX()
    #     new_nx = len(edges_array) - 1
        
    #     if old_nx == new_nx:
    #         return h1_in.Clone(new_name)
            
    #     h1_out = ROOT.TH1D(new_name, h1_in.GetTitle(), new_nx, edges_array)
    #     h1_out.Sumw2()
        
    #     for ix in range(1, old_nx + 1):
    #         content = h1_in.GetBinContent(ix)
    #         error   = h1_in.GetBinError(ix)
            
    #         if content != 0 or error != 0:
    #             x_center = h1_in.GetXaxis().GetBinCenter(ix)
    #             new_ix = h1_out.GetXaxis().FindBin(x_center)
                
    #             # Manual Accumulation
    #             cur_val = h1_out.GetBinContent(new_ix)
    #             cur_err = h1_out.GetBinError(new_ix)
                
    #             h1_out.SetBinContent(new_ix, cur_val + content)
    #             h1_out.SetBinError(new_ix, math.sqrt(cur_err**2 + error**2))
                
    #     return h1_out

    ###################################################
    ### Can't rebin efficiency, goddamnit! The data is already made to be rebinned!
    ###################################################
    # # A) Rebin 2D maps (Just for saving to "CorrectionsUsed", not used for calculation anymore)
    # hEffAcc_reb2d = rebin2d_in_x_to_edges(hEffAcc, edges_array, "hEffAcc_reb2d")
    # hSignalLoss_reb2d = rebin2d_in_x_to_edges(hSignalLoss, edges_array, "hSignalLoss_reb2d")
    #
    # # B) Rebin 1D Efficiency & Signal Loss (Used for calculation)
    # # We use the new rebin1d helper to ensure errors are preserved
    # hEff_fullZ = rebin1d_to_edges(hEff1D_Full_in, edges_array, "hEff_fullZ")
    # hEff_Zpos  = rebin1d_to_edges(hEff1D_Pos_in,  edges_array, "hEff_Zpos")
    # hEff_Zneg  = rebin1d_to_edges(hEff1D_Neg_in,  edges_array, "hEff_Zneg")
    #
    # hSig_fullZ = rebin1d_to_edges(hSig1D_Full_in, edges_array, "hSig_fullZ")
    # hSig_Zpos  = rebin1d_to_edges(hSig1D_Pos_in,  edges_array, "hSig_Zpos")
    # hSig_Zneg  = rebin1d_to_edges(hSig1D_Neg_in,  edges_array, "hSig_Zneg")
    #
    # hEffAcc_reb2d.Sumw2()
    # hSignalLoss_reb2d.Sumw2()
        # These are the right corrections:
    hEff_fullZ = hEff1D_Full_in.Clone("hEff_fullZ")
    hEff_Zpos = hEff1D_Pos_in.Clone("hEff_Zpos")
    hEff_Zneg = hEff1D_Neg_in.Clone("hEff_Zneg")

    hSig_fullZ = hSig1D_Full_in.Clone("hSig_fullZ")
    hSig_Zpos = hSig1D_Pos_in.Clone("hSig_Zpos")
    hSig_Zneg = hSig1D_Neg_in.Clone("hSig_Zneg")

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
            # # It makes more sense to divide the Full Z value by two,
            # # as we are comparing two halves of it:
            # B /= 2.0
            # eB /= 2.0

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
        """Convert histogram to (1/Nev)*(1/dpT)*counts."""
        out = h.Clone(label)
        out.Scale(1, "width") # Shorthands
        out.Scale(1./Nev) # Shorthands
        return out

    # Event normalization = total reconstructed events (same as in corrections)
    # NEV = hEvtSplitting.GetBinContent(1) * hEvtLoss.GetBinContent(1)
    # if NEV <= 0:
    #     raise RuntimeError("Got a null number of events!")
    # Getting the actual number of events from the DATA part -- We don't want the Nev from MC!!!
    hCentrality = fin.Get(f"{main_dir}/hEventCentrality") # Just to get the number of entries!
    NEV = hCentrality.Integral()
    print("\tNEV from hCentrality:", NEV)

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
    if not running_OO:
        text.DrawLatex(0.18, 0.86, "pp, 13.6 TeV, #Lambda")
    else:
        text.DrawLatex(0.18, 0.86, "OO, 5.36 TeV, #Lambda")

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

    if not running_OO:
        text.DrawLatex(0.18, 0.86, "pp, 13.6 TeV, #Lambda")
    else:
        text.DrawLatex(0.18, 0.86, "OO, 5.36 TeV, #Lambda")

    diag_canvases.append(c_eff)

    # ------------------------------
    # PLOT 3: Raw spectra + Ratio subplot
    # ------------------------------

    # Normalize the raw spectra
    hRaw_fullZ_norm = normalize_spectrum(h_pt_fullZ, NEV, "hRaw_fullZ_norm")
    hRaw_Zpos_norm  = normalize_spectrum(h_pt_Zpos,  NEV, "hRaw_Zpos_norm")
    hRaw_Zneg_norm  = normalize_spectrum(h_pt_Zneg,  NEV, "hRaw_Zneg_norm")
    # Also normalizing by the size of the dz projection -- These should be d^{2}N/dp_{T}dz plots!
        # Getting the range of Z values used:
    yaxis = h_pt_fullZ_2D.GetYaxis()  # Z position
    full_z_range_size = abs(yaxis.GetXmax() - yaxis.GetXmin())
    hRaw_fullZ_norm.Scale(1./full_z_range_size)
    hRaw_Zpos_norm.Scale(1./(full_z_range_size/2))
    hRaw_Zneg_norm.Scale(1./(full_z_range_size/2))

    # Ratio histos
    # No need for halves when we already have d^2N/dpT dz!
    # # Halving the FullZ histogram for direct comparison with the corrected histograms:
    # hRaw_fullZ_norm_halved = hRaw_fullZ_norm.Clone("hRaw_fullZ_norm_halved")
    # hRaw_fullZ_norm_halved.Scale(0.5)
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
        h.SetTitle(";p_{T} (GeV/c);(1/N_{evt}) d^{2}N/dp_{T}dz")
        h.GetYaxis().SetTitleOffset(1.3)

    # hRaw_fullZ_norm.Draw("E")
    hRaw_fullZ_norm.Draw("E")
    hRaw_Zpos_norm.Draw("E SAME")
    hRaw_Zneg_norm.Draw("E SAME")

    leg3 = ROOT.TLegend(0.60, 0.70, 0.88, 0.88)
    leg3.AddEntry(hRaw_fullZ_norm, "#Lambda full Z", "lep")
    leg3.AddEntry(hRaw_Zpos_norm,  "#Lambda Z >= 0", "lep")
    leg3.AddEntry(hRaw_Zneg_norm,  "#Lambda Z < 0", "lep")
    leg3.Draw()

    if not running_OO:
        text.DrawLatex(0.18, 0.86, "pp, 13.6 TeV, #Lambda")
    else:
        text.DrawLatex(0.18, 0.86, "OO, 5.36 TeV, #Lambda")

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
    # Also normalizing by the size of the dz projection -- These should be d^{2}N/dp_{T}dz plots!
    hCorr_fullZ_norm.Scale(1./full_z_range_size)
    hCorr_Zpos_norm.Scale(1./(full_z_range_size/2))
    hCorr_Zneg_norm.Scale(1./(full_z_range_size/2))

    # Checking to see the number of particles per event in those bins:
    print("Integrated d^2N/dpTdz spectra, Full Z:", hCorr_fullZ_norm.Integral("width") * full_z_range_size)
    print("Integrated d^2N/dpTdz spectra, Z>=0:", hCorr_Zpos_norm.Integral("width")* (full_z_range_size/2))
    print("Integrated d^2N/dpTdz spectra, Z<0:", hCorr_Zneg_norm.Integral("width")* (full_z_range_size/2))

    # No need for halves when we already have d^2N/dpT dz!
    # # Also halving the corrected spectra:
    # hCorr_fullZ_norm_halved = hCorr_fullZ_norm.Clone("hCorr_fullZ_norm_halved")
    # hCorr_fullZ_norm_halved.Scale(0.5)
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
        h.SetTitle(";p_{T} (GeV/c);(1/N_{evt}) d^{2}N/dp_{T}dz")
        h.GetYaxis().SetTitleOffset(1.3)

    # hCorr_fullZ_norm.Draw("E")
    hCorr_fullZ_norm.Draw("E")
    hCorr_Zpos_norm.Draw("E SAME")
    hCorr_Zneg_norm.Draw("E SAME")

    leg4 = ROOT.TLegend(0.60, 0.70, 0.88, 0.88)
    leg4.AddEntry(hCorr_fullZ_norm, "#Lambda full Z", "lep")
    leg4.AddEntry(hCorr_Zpos_norm,  "#Lambda Z >= 0", "lep")
    leg4.AddEntry(hCorr_Zneg_norm,  "#Lambda Z < 0", "lep")
    leg4.Draw()

    if not running_OO:
        text.DrawLatex(0.18, 0.86, "pp, 13.6 TeV, #Lambda")
    else:
        text.DrawLatex(0.18, 0.86, "OO, 5.36 TeV, #Lambda")

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
    # 9.6) MC Closure test
    # ------------------------------
    if running_MC:
        closure_canvases = []
        # # h2d_PtVsZ_GenMC = fcor.Get("h3dGenLambdaVsZVsMultMC_RecoedEvt_yz") # Despite the name, this is the TH2D of pT vs Z
        #     # The actual proper version should not be Recoed!!! It should be the actual Generated data!
        # # h2d_PtVsZ_GenMC = fcor.Get("hGenLambdaPtVsZ")
        # h2d_PtVsZ_GenMC = fcor.Get("hGenLambdaPtVsZ_rebinned") # Now using the rebinned version to make sure it matches Reco MC
        
        # # Creating TH1D structures as in SignalExtraction, with the same binning as the TH2D at first:
        # def make_TH1D(name, title, axis):
        #     """
        #     Build a TH1D using the binning of a TAxis.
        #     Works for both fixed and variable binning.
        #     """
        #     nbins = axis.GetNbins()

        #     if axis.GetXbins().GetSize() > 0:
        #         # Variable binning
        #         bins = axis.GetXbins().GetArray()
        #         h = ROOT.TH1D(name, title, nbins, bins)
        #     else:
        #         # Fixed binning
        #         xmin = axis.GetXmin()
        #         xmax = axis.GetXmax()
        #         h = ROOT.TH1D(name, title, nbins, xmin, xmax)

        #     h.Sumw2()
        #     return h
        
        # xaxis = h2d_PtVsZ_GenMC.GetXaxis() # pT
        # yaxis = h2d_PtVsZ_GenMC.GetYaxis() # Z position

        # n_pt_bins = xaxis.GetNbins()
        # n_z_bins  = yaxis.GetNbins()

        # h_pt_fullZ_GenMC = make_TH1D(
        #     "h_pt_fullZ_GenMC",
        #     "#Lambda p_{T} spectrum (full Z); p_{T} (GeV/c); dN/dp_{T}",
        #     xaxis
        # )
        # # h_pt_fullZ_GenMC = h_pt_fullZ.Clone("h_pt_fullZ_GenMC")
        # hPt_NegZ_GenMC = h_pt_fullZ_GenMC.Clone("hPt_NegZ_GenMC")
        # hPt_PosZ_GenMC = h_pt_fullZ_GenMC.Clone("hPt_PosZ_GenMC")
        
        # h_pt_fullZ_GenMC.Reset("ICES")
        # hPt_NegZ_GenMC.Reset("ICES")
        # hPt_PosZ_GenMC.Reset("ICES")
        
        # # =========================
        # # 8.1) Loop over pT bins and summing Z bins
        # # =========================
        # for ipt in range(1, n_pt_bins + 1):

        #     sum_full   = 0.0
        #     err2_full  = 0.0

        #     sum_pos    = 0.0
        #     err2_pos   = 0.0

        #     sum_neg    = 0.0
        #     err2_neg   = 0.0

        #     for iz in range(1, n_z_bins + 1):

        #         z_center = yaxis.GetBinCenter(iz)

        #         val = h2d_PtVsZ_GenMC.GetBinContent(ipt, iz)
        #         err = h2d_PtVsZ_GenMC.GetBinError(ipt, iz)

        #         # Full Z
        #         sum_full  += val
        #         err2_full += err * err

        #         # Z >= 0
        #         if z_center >= 0.0:
        #             sum_pos  += val
        #             err2_pos += err * err
        #         else:
        #             sum_neg  += val
        #             err2_neg += err * err

        #     # Fill pT spectra
        #     h_pt_fullZ_GenMC.SetBinContent(ipt, sum_full)
        #     h_pt_fullZ_GenMC.SetBinError(ipt, math.sqrt(err2_full))

        #     hPt_PosZ_GenMC.SetBinContent(ipt, sum_pos)
        #     hPt_PosZ_GenMC.SetBinError(ipt, math.sqrt(err2_pos))

        #     hPt_NegZ_GenMC.SetBinContent(ipt, sum_neg)
        #     hPt_NegZ_GenMC.SetBinError(ipt, math.sqrt(err2_neg))
            
        # # Applying some rebins to make this the same as the Reco MC -- Binning should actually already be the same! No need to do so!
        # # h_pt_fullZ_GenMC = rebin1d_to_edges(h_pt_fullZ_GenMC, edges_array, "h_pt_fullZ_GenMC_rebinned")
        # # hPt_NegZ_GenMC  = rebin1d_to_edges(hPt_NegZ_GenMC,  edges_array, "hPt_NegZ_GenMC_rebinned")
        # # hPt_PosZ_GenMC  = rebin1d_to_edges(hPt_PosZ_GenMC,  edges_array, "hPt_PosZ_GenMC_rebinned")
        
        # Getting the spectra already in a TH1D form:
            # Could use the full block above, but this is way shorter!
        h_pt_fullZ_GenMC = fcor.Get("hLambdaGeneratorLevelPt_full")
        hPt_NegZ_GenMC = fcor.Get("hLambdaGeneratorLevelPt_pos")
        hPt_PosZ_GenMC = fcor.Get("hLambdaGeneratorLevelPt_neg")
        
        # Now normalizing these plots:
        hGenEvent = fcor.Get("hGenEvent")
        NEV_Gen = hGenEvent.GetBinContent(1)
        print("\tNumber of events generated is {}, while NEV reconstructed is {}".format(NEV_Gen, NEV))
        h_pt_fullZ_GenMC_norm = normalize_spectrum(h_pt_fullZ_GenMC, NEV_Gen, "h_pt_fullZ_GenMC_norm")
        hPt_NegZ_GenMC_norm  = normalize_spectrum(hPt_NegZ_GenMC, NEV_Gen, "hPt_NegZ_GenMC_norm")
        hPt_PosZ_GenMC_norm  = normalize_spectrum(hPt_PosZ_GenMC, NEV_Gen, "hPt_PosZ_GenMC_norm")
        # Also normalizing by the size of the dz projection -- These should be d^{2}N/dp_{T}dz plots!
        h_pt_fullZ_GenMC_norm.Scale(1./full_z_range_size)
        hPt_NegZ_GenMC_norm.Scale(1./(full_z_range_size/2))
        hPt_PosZ_GenMC_norm.Scale(1./(full_z_range_size/2))

        # Checking to see the number of particles per event in those bins:
        print("\tIntegrated d^2N/dpTdz spectra, MC_Gen, Full Z:", h_pt_fullZ_GenMC_norm.Integral("width") * full_z_range_size)
        print("\tIntegrated d^2N/dpTdz spectra, MC_Gen, Z>=0:", hPt_NegZ_GenMC_norm.Integral("width")* (full_z_range_size/2))
        print("\tIntegrated d^2N/dpTdz spectra, MC_Gen, Z<0:", hPt_PosZ_GenMC_norm.Integral("width")* (full_z_range_size/2))
        
        # Now drawing these plots along with the MC corrected spectra -- Things should match exactly!
            # This is a copy of the related section of the diagnostics plots
        # ------------------------------
        # Corrected spectra + Ratios with the Generated spectra
        # ------------------------------
        # hCorr_fullZ_norm_rebinned = rebin1d_to_edges(hCorr_fullZ_norm, edges_array, "hCorr_fullZ_norm_rebinned")
        # hCorr_Zpos_norm_rebinned  = rebin1d_to_edges(hCorr_Zpos_norm,  edges_array, "hCorr_Zpos_norm_rebinned")
        # hCorr_Zneg_norm_rebinned  = rebin1d_to_edges(hCorr_Zneg_norm,  edges_array, "hCorr_Zneg_norm_rebinned")

        hCorrectedReco_over_Gen_ratio_FullZ = make_ratio_hist(hCorr_fullZ_norm, h_pt_fullZ_GenMC_norm, "hCorrectedReco_over_Gen_ratio_FullZ")
        hCorrectedReco_over_Gen_ratio_Zpos = make_ratio_hist(hCorr_Zpos_norm, hPt_NegZ_GenMC_norm, "hCorrectedReco_over_Gen_ratio_Zpos")
        hCorrectedReco_over_Gen_ratio_Zneg = make_ratio_hist(hCorr_Zneg_norm, hPt_PosZ_GenMC_norm, "hCorrectedReco_over_Gen_ratio_Zneg")

        c_corr = ROOT.TCanvas("cCorrectedSpectra_RatioWithGenMC", "Corrected pT spectra - RatioWithGenMC", 1200, 1000)
        pad1cGenVsReco = ROOT.TPad("pad1cGenVsReco", "pad1cGenVsReco", 0, 0.35, 1, 1.0)
        pad2cGenVsReco = ROOT.TPad("pad2cGenVsReco", "pad2cGenVsReco", 0, 0.00, 1, 0.35)

        pad1cGenVsReco.SetBottomMargin(0.02)
        pad1cGenVsReco.SetLeftMargin(0.12)
        pad1cGenVsReco.SetRightMargin(0.05)

        pad2cGenVsReco.SetTopMargin(0.03)
        pad2cGenVsReco.SetBottomMargin(0.35)
        pad2cGenVsReco.SetLeftMargin(0.12)
        pad2cGenVsReco.SetRightMargin(0.05)

        pad1cGenVsReco.Draw()
        pad2cGenVsReco.Draw()

        # Upper panel
        pad1cGenVsReco.cd()
        hCorr_fullZ_norm.GetXaxis().SetRangeUser(0, 20)
        hCorr_Zpos_norm.GetXaxis().SetRangeUser(0, 20)
        hCorr_Zneg_norm.GetXaxis().SetRangeUser(0, 20)
        h_pt_fullZ_GenMC_norm.GetXaxis().SetRangeUser(0, 20)
        hPt_NegZ_GenMC_norm.GetXaxis().SetRangeUser(0, 20)
        hPt_PosZ_GenMC_norm.GetXaxis().SetRangeUser(0, 20)
        for h, col in [
            (hCorr_fullZ_norm, ROOT.kBlack),
            (hCorr_Zpos_norm, ROOT.kRed+1),
            (hCorr_Zneg_norm, ROOT.kBlue+1),
            (h_pt_fullZ_GenMC_norm, ROOT.kGreen+2),
            (hPt_NegZ_GenMC_norm, ROOT.kMagenta+1),
            (hPt_PosZ_GenMC_norm, ROOT.kOrange+7),
        ]:
            h.SetLineColor(col)
            h.SetMarkerColor(col)
            h.SetMarkerStyle(20)

            h.SetTitle(";p_{T} (GeV/c);(1/N_{evt}) d^{2}N/dp_{T}dz")
            h.GetYaxis().SetTitleOffset(1.3)

        # hCorr_fullZ_norm.Draw("E")
        hCorr_fullZ_norm.Draw("E")
        hCorr_Zpos_norm.Draw("E SAME")
        hCorr_Zneg_norm.Draw("E SAME")
        h_pt_fullZ_GenMC_norm.Draw("E SAME")
        hPt_NegZ_GenMC_norm.Draw("E SAME")
        hPt_PosZ_GenMC_norm.Draw("E SAME")

        leg4 = ROOT.TLegend(0.60, 0.52, 0.88, 0.88)
        leg4.AddEntry(hCorr_fullZ_norm, "#Lambda MC corr. full Z", "lep")
        leg4.AddEntry(hCorr_Zpos_norm,  "#Lambda MC corr. Z >= 0", "lep")
        leg4.AddEntry(hCorr_Zneg_norm,  "#Lambda MC corr. Z < 0", "lep")
        leg4.AddEntry(h_pt_fullZ_GenMC_norm, "#Lambda Gen full Z", "lep")
        leg4.AddEntry(hPt_NegZ_GenMC_norm,  "#Lambda Gen Z >= 0", "lep")
        leg4.AddEntry(hPt_PosZ_GenMC_norm,  "#Lambda Gen Z < 0", "lep")
        leg4.Draw()

        if not running_OO:
            text.DrawLatex(0.18, 0.86, "pp, 13.6 TeV, #Lambda")
        else:
            text.DrawLatex(0.18, 0.86, "OO, 5.36 TeV, #Lambda")

        # Lower panel
        pad2cGenVsReco.cd()
        hCorrectedReco_over_Gen_ratio_FullZ.GetXaxis().SetRangeUser(0, 20)
        hCorrectedReco_over_Gen_ratio_Zpos.GetXaxis().SetRangeUser(0, 20)
        hCorrectedReco_over_Gen_ratio_Zneg.GetXaxis().SetRangeUser(0, 20)
        
        hCorrectedReco_over_Gen_ratio_FullZ.SetLineColor(ROOT.kBlack)
        hCorrectedReco_over_Gen_ratio_Zpos.SetLineColor(ROOT.kRed+1)
        hCorrectedReco_over_Gen_ratio_Zneg.SetLineColor(ROOT.kBlue+1)
        
        hCorrectedReco_over_Gen_ratio_FullZ.SetMarkerColor(ROOT.kBlack)
        hCorrectedReco_over_Gen_ratio_Zpos.SetMarkerColor(ROOT.kRed+1)
        hCorrectedReco_over_Gen_ratio_Zneg.SetMarkerColor(ROOT.kBlue+1)
        
        hCorrectedReco_over_Gen_ratio_FullZ.SetMarkerStyle(20)
        hCorrectedReco_over_Gen_ratio_Zpos.SetMarkerStyle(20)
        hCorrectedReco_over_Gen_ratio_Zneg.SetMarkerStyle(20)

        hCorrectedReco_over_Gen_ratio_FullZ.SetTitle(";p_{T} (GeV/c);Corrected Reco MC / Gen MC")

        hCorrectedReco_over_Gen_ratio_FullZ.Draw("E")
        hCorrectedReco_over_Gen_ratio_Zpos.Draw("E SAME")
        hCorrectedReco_over_Gen_ratio_Zneg.Draw("E SAME")

        closure_canvases.append(c_corr)

    # ------------------------------
    # 10) prepare output file (in same folder as analysis_results_path)
    # ------------------------------
    # input_dir = os.path.dirname(os.path.abspath(analysis_results_path))
    # output_name = os.path.join(input_dir, "CorrectedpTSpectra.root")
    
    # Determine config tag from input file name:
    # Expected pattern: "AnalysisResults-<config>.root"
    base = os.path.basename(analysis_results_path)  # e.g. AnalysisResults-hasTPCnoITS.root

    match = re.match(r"AnalysisResults-(.+)\.root", base)
    if match:
        config_tag = match.group(1)  # e.g. "hasTPCnoITS"
    else:
        print("Warning: could not extract config tag from input filename. Using 'UnknownConfig'.")
        config_tag = "UnknownConfig"

    # Prepare output file name
    input_dir = os.path.dirname(os.path.abspath(analysis_results_path))
    output_name = os.path.join(input_dir, f"CorrectedpTSpectra-{config_tag}.root")

    print(f"Writing corrected spectra to: {output_name}")
    fout = ROOT.TFile.Open(output_name, "RECREATE")

    # Create subdirectories
    dir_main = fout.mkdir("ProcessedSpectra")
    dir_corr = fout.mkdir("CorrectionsUsed")
    dir_diag = fout.mkdir("Diagnostics")
    if running_MC:
        dir_closure = fout.mkdir("ClosureTest")

    # Write spectra and processed items
    dir_main.cd()
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
    # hEffAcc_reb2d.Write()
    # hSignalLoss_reb2d.Write()
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
    
    # Writ MC closure canvases:
    dir_closure.cd()
    for c in closure_canvases:
        c.Write()

    fout.Close()
    fin.Close()
    fcor.Close()
    fsigExtractedData.Close()

    print(f"Done. Output saved to: {output_name}")


# --------------------------- CLI --------------------------------
def usage_and_exit():
    print("Usage: python ApplyEfficiency.py <AnalysisResults.root> <EfficiencyCorrections.root>")
    sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        usage_and_exit()
    analysis_results_path = sys.argv[1]
    corrections_path = sys.argv[2]
    apply_efficiency(analysis_results_path, corrections_path)
