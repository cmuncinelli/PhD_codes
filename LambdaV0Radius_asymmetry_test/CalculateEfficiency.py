#!/usr/bin/env python3
"""
CalculateEfficiency_Lambda.py

Cleaned and focused version of Gianni's CalculateEfficiency.py adapted for
Lambda (Λ) efficiency calculations.

This script reads a single AnalysisResults ROOT file and writes an
EfficiencyCorrections_*.root file in the same folder as the input file.

Usage:
    python CalculateEfficiency.py /path/to/MCAnalysisResults_<Dataset>.root \
        --main-dir MyMainDir --species Lambda

Author: adapted by ChatGPT for the user (keeps original code comments and helper functions)
Based on: Gianni Shigeru Setoue Liveraro's CalculateEfficiency.py
"""

# -------------------------------  DEFINITIONS --------------------------------
# EFFICIENCY X ACC X B.R:
# N. Reco in the selected reco colls / N. Gen V0s with at least one reco coll
#
# SIGNAL LOSS:
# N. Gen V0s with at least one reco coll / N. Gen V0s in ALL Generated colls
#
# EVENT LOSS:
# N. gen evts with least one reco coll (after evt sel) / N. of ALL gen evts
#
# EVENT SPLITTING:
# N. gen evts with least one reco coll (after evt sel) / N. of reco colls
#
# FULL pT SPECTRA CORRECTION:
# (1/EFFICIENCY) * (EVENT LOSS / SIGNAL LOSS * EVENT SPLITTING)
# -----------------------------------------------------------------------------

# Imports
import os
import sys
import argparse
import array

import ROOT
import numpy as np

ROOT.gErrorIgnoreLevel = 6001 # Ignore "Error in <TList::Clear>:...""

# ------------------- Helper Functions --------------------------------
def getBinEdges(hInputSpectra):
    """
    Return an array.array('d') with bin edges from a 1D histogram.
    Used to rebin to the same binning as raw spectra.
    """
    nBinsX = hInputSpectra.GetNbinsX()
    XAxis = hInputSpectra.GetXaxis()
    pTbins = [XAxis.GetBinLowEdge(i + 1) for i in range(nBinsX + 1)]
    pTbinsArray = array.array('d', pTbins)
    return pTbinsArray

def project_th2d(hist2D_0, axis, x_range=None, y_range=None, rebin_factor=None, number=-1):
    """
    Project a TH2D into TH1D along given axis with optional axis ranges and rebin.
    axis: 'X' or 'Y' (which axis to keep after projection)
    x_range/y_range: tuples (low, high) in axis units (not bin indices)
    rebin_factor: integer rebin factor (>1 means Rebin)
    number: identifier used in projection name
    """
    if not hist2D_0:
        raise ValueError("Invalid histogram provided to project_th2d.")
    
    x_bin_min, x_bin_max = 1, hist2D_0.GetNbinsX()
    y_bin_min, y_bin_max = 1, hist2D_0.GetNbinsY()
    
    if x_range:
        x_bin_min = hist2D_0.GetXaxis().FindBin(x_range[0])
        x_bin_max = hist2D_0.GetXaxis().FindBin(x_range[1])
    
    if y_range:
        y_bin_min = hist2D_0.GetYaxis().FindBin(y_range[0])
        y_bin_max = hist2D_0.GetYaxis().FindBin(y_range[1])
    
    if axis == 'X':
        hist1D = hist2D_0.ProjectionX(f"proj_{number}", y_bin_min, y_bin_max)
    else:
        hist1D = hist2D_0.ProjectionY(f"proj_{number}", x_bin_min, x_bin_max)
    
    if rebin_factor and rebin_factor > 1:
        hist1D.Rebin(rebin_factor)
    
    return hist1D

# Actually need 3D projections now that we are working with:
# h3dGenLambdaVsZ                        # (centrality, Z, pT) - all generated Lambdas
# h3dGenLambdaVsZVsMultMC_RecoedEvt     # (Nch, Z, pT) - generated Lambda with >=1 reco track
# h3dGenLambdaVsZVsMultMC               # (Nch, Z, pT) - all generated Lambda (same as above but maybe different binnings)
def project_th3d(hist3D, axis, 
                 event_range=None, z_range=None, pt_range=None,
                 number=-1):
    if not hist3D:
        raise ValueError("Invalid TH3D provided.")

    # Determine bin limits
    e_min, e_max = 1, hist3D.GetNbinsX()
    z_min, z_max = 1, hist3D.GetNbinsY()
    p_min, p_max = 1, hist3D.GetNbinsZ()

    if event_range:
        e_min = hist3D.GetXaxis().FindBin(event_range[0])
        e_max = hist3D.GetXaxis().FindBin(event_range[1])

    if z_range:
        z_min = hist3D.GetYaxis().FindBin(z_range[0])
        z_max = hist3D.GetYaxis().FindBin(z_range[1])

    if pt_range:
        p_min = hist3D.GetZaxis().FindBin(pt_range[0])
        p_max = hist3D.GetZaxis().FindBin(pt_range[1])

    # Axis selection
    if axis == "pt":
        return hist3D.ProjectionZ(f"projZ_{number}", e_min, e_max, z_min, z_max)
    if axis == "Z":
        return hist3D.ProjectionY(f"projY_{number}", e_min, e_max, p_min, p_max)
    if axis == "evt":
        return hist3D.ProjectionX(f"projX_{number}", z_min, z_max, p_min, p_max)

    raise ValueError("axis must be one of: 'pt','Z','evt'")



# ------------------- Settings (can be adapted or exposed to CLI) --------------------------------
plot_xmin, plot_xmax = 0, 5.0
plot_ymin, plot_ymax = 0.0, 0.06
rebin_factor = 1
x_range = None
y_range = None
z_range = None

# ------------------- Main calculation --------------------------------
def calculate_efficiency(input_root_path, main_dir="asymmetric_rapidity_test", species="Lambda"):
    """
    input_root_path : path to MCAnalysisResults_*.root
    main_dir : path inside the root file where Gen histograms are stored
    species : 'Lambda', 'ALambda', or 'LambdaALambda'
    """
    # Basic checks
    if not os.path.isfile(input_root_path):
        raise FileNotFoundError(f"Input file not found: {input_root_path}")
    
    # Output filename in same folder as input file
    input_dir = os.path.dirname(os.path.abspath(input_root_path))

    # Extract the base name without extension
    base_input_name = os.path.splitext(os.path.basename(input_root_path))[0]

    # Remove leading "AnalysisResults-" if present
    clean_name = base_input_name.replace("AnalysisResults-", "", 1)

    output_name = f"EfficiencyCorrections-{clean_name}.root"
    output_path = os.path.join(input_dir, output_name)


    print(f"Input file : {input_root_path}")
    print(f"Main dir   : {main_dir}")
    print(f"Species    : {species}")
    print(f"Output file: {output_path}")

    # Open file
    fin = ROOT.TFile.Open(input_root_path, "READ")
    if not fin or fin.IsZombie():
        raise RuntimeError(f"Could not open input ROOT file: {input_root_path}")

    # Determine histogram names for Lambda
    # The user wanted replacements:
    # h2dGenSigma0VsMultMC -> h2dGenLambdaVsMultMC
    # h2dGenSigma0VsMultMC_RecoedEvt -> h2dGenLambdaVsMultMC_RecoedEvt
    # For anti-Lambda we assume suffix "AntiLambda" following original pattern.
    # Combined option sums Lambda + AntiLambda.
    # def get_gen_hist_names2d(spec):
    #     if spec == "Lambda":
    #         h2d_gen = f"{main_dir}/h2dGenLambdaVsMultMC"
    #         h2d_gen_recoed = f"{main_dir}/h2dGenLambdaVsMultMC_RecoedEvt"
    #         return h2d_gen_recoed, h2d_gen
    #     elif spec == "ALambda":
    #         h2d_gen = f"{main_dir}/h2dGenAntiLambdaVsMultMC"
    #         h2d_gen_recoed = f"{main_dir}/h2dGenAntiLambdaVsMultMC_RecoedEvt"
    #         return h2d_gen_recoed, h2d_gen
    #     elif spec == "LambdaALambda":
    #         # will fetch both and sum
    #         return ("COMBINED", "COMBINED")
    #     else:
    #         raise ValueError("Unsupported species. Use Lambda, ALambda or LambdaALambda.")
    
    # Version for 3D:
    def get_gen_hist_names(spec):
        if spec == "Lambda":
            h3d_gen = f"{main_dir}/h3dGenLambdaVsZVsMultMC"
            h3d_gen_recoed = f"{main_dir}/h3dGenLambdaVsZVsMultMC_RecoedEvt"
            return h3d_gen_recoed, h3d_gen
        elif spec == "ALambda":
            h3d_gen = f"{main_dir}/h2dGenAntiLambdaVsZVsMultMC"
            h3d_gen_recoed = f"{main_dir}/h2dGenAntiLambdaVsZVsMultMC_RecoedEvt"
            return h3d_gen_recoed, h3d_gen
        elif spec == "LambdaALambda":
            # will fetch both and sum
            return ("COMBINED", "COMBINED")
        else:
            raise ValueError("Unsupported species. Use Lambda, ALambda or LambdaALambda.")

    # Fetch histograms depending on species
    # if species in ("Lambda", "ALambda"):
        # ------------------- TH3D MC LOADING + PROJECTION -------------------

    h3d_gen_recoed_name, h3d_gen_name = get_gen_hist_names(species)

    # === Load TH3D: Generated Λ with ≥1 reco event ===
    h3dGen_RecoedEvt = fin.Get(h3d_gen_recoed_name)
    if not h3dGen_RecoedEvt:
        print(f"Warning: could not find {h3d_gen_recoed_name} in file.")
    else:
        # Project TH3D(event, Z, pT) → TH2D(pT, Z)
        # ROOT syntax: Project3D("zy") → x'=Z, y'=Y? No:
        #  For TH3D( X=event , Y=Z , Z=pT ):
        #    Project3D("zy") → Z→X', Y→Y' = (pT , Z)
        hGenRecoEvt = h3dGen_RecoedEvt.Project3D("zy")
        hGenRecoEvt.SetName(f"Gen{species}_RecoEvts")
        hGenRecoEvt.SetTitle(f"Gen{species}_RecoEvts")

    # === Load TH3D: All generated Λ ===
    h3dGen = fin.Get(h3d_gen_name)
    if not h3dGen:
        print(f"Warning: could not find {h3d_gen_name} in file.")
    else:
        # Same projection: TH3D(event, Z, pT) → TH2D(pT, Z)
        hGen = h3dGen.Project3D("zy")
        hGen.SetName(f"Gen{species}_GenEvts")
        hGen.SetTitle(f"Gen{species}_GenEvts")
    
    # else:  # LambdaALambda combined
    #     # Lambda
    #     h2dGen_RecoedEvt_L = fin.Get(f"{main_dir}/h2dGenLambdaVsMultMC_RecoedEvt")
    #     hGenReco_L = project_th2d(h2dGen_RecoedEvt_L, "Y", x_range, y_range, rebin_factor, "GenLambda_RecoEvts")

    #     h2dGen_L = fin.Get(f"{main_dir}/h2dGenLambdaVsMultMC")
    #     hGen_L = project_th2d(h2dGen_L, "Y", x_range, y_range, rebin_factor, "GenLambda_GenEvts")

    #     # Anti-Lambda
    #     h2dGen_RecoedEvt_A = fin.Get(f"{main_dir}/h2dGenAntiLambdaVsMultMC_RecoedEvt")
    #     hGenReco_A = project_th2d(h2dGen_RecoedEvt_A, "Y", x_range, y_range, rebin_factor, "GenALambda_RecoEvts")

    #     h2dGen_A = fin.Get(f"{main_dir}/h2dGenAntiLambdaVsMultMC")
    #     hGen_A = project_th2d(h2dGen_A, "Y", x_range, y_range, rebin_factor, "GenALambda_GenEvts")

    #     # Sum
    #     hGenRecoEvt = hGenReco_L.Clone("GenLambdaALambda_RecoEvts")
    #     hGenRecoEvt.Add(hGenReco_A)

    #     hGen = hGen_L.Clone("GenLambdaALambda_GenEvts")
    #     hGen.Add(hGen_A)

    # =====================================================================
    # N. Reco in the selected reco colls : NOW USING hLambdaPtZ (TH2D)
    # =====================================================================

    # raw MC spectra: TH2D(pT,Z)
    rawMCSpectra = fin.Get(f"{main_dir}/Lambda/hLambdaPtZ")
    if not rawMCSpectra:
        raise RuntimeError(f"Could not find '{main_dir}/Lambda/hLambdaPtZ' in the input file.")
    rawMCSpectra.Sumw2()

    # =====================================================================
    # Event histograms
    # =====================================================================

    h2dGenAndRecoEvents = fin.Get(f"{main_dir}/hGenEvents")
    if not h2dGenAndRecoEvents:
        raise RuntimeError(f"Could not find {main_dir}/hGenEvents in input file.")

    # For event counts we still want 1D, even though MC histos are 3D.
    # Project X (event-axis) → 1D two-bin hist:
    hGenAndRecoEvts = h2dGenAndRecoEvents.ProjectionX("EventGenReco")
    hGenAndRecoEvts.Sumw2()

    # All gen events = bin 1
    hGenEvent = ROOT.TH1D("hGenEvent", "hGenEvent", 1, 0, 1)
    hGenEvent.SetBinContent(1, hGenAndRecoEvts.GetBinContent(1))
    hGenEvent.SetBinError(1, hGenAndRecoEvts.GetBinError(1))

    # Gen events with ≥1 reco = bin 2
    hGenRecoEvt_events = ROOT.TH1D("hGenRecoEvt_events", "hGenRecoEvt_events", 1, 0, 1)
    hGenRecoEvt_events.SetBinContent(1, hGenAndRecoEvts.GetBinContent(2))
    hGenRecoEvt_events.SetBinError(1, hGenAndRecoEvts.GetBinError(2))

    print("Gen Events:", hGenAndRecoEvts.GetBinContent(1))
    print("Gen Reco Events:", hGenAndRecoEvts.GetBinContent(2))

    # =====================================================================
    # N reco colls: hEventCentrality
    # =====================================================================

    hCentrality = fin.Get(f"{main_dir}/hEventCentrality")
    if not hCentrality:
        print("Warning: hEventCentrality not found, creating fallback from event histogram.")
        hCentrality = ROOT.TH1D("hEventCentrality", "hEventCentrality", 1, 0, 1)
        hCentrality.SetBinContent(1, hGenAndRecoEvts.GetBinContent(2))

    try:
        hCentrality.Rebin(12)
    except Exception:
        pass

    hRecoEvents = ROOT.TH1D("hRecoEvents", "hRecoEvents", 1, 0, 1)
    hRecoEvents.SetBinContent(1, hCentrality.Integral())
    hRecoEvents.Sumw2()

    # =====================================================================
    # PROPAGATE ERRORS
    # =====================================================================

    hGenEvent.Sumw2()
    hGenRecoEvt_events.Sumw2()
    hGenAndRecoEvts.Sumw2()
    hCentrality.Sumw2()
    hRecoEvents.Sumw2()

    # =====================================================================
    # Now compute efficiency using 2D histograms (pT,Z)
    # =====================================================================

    # --- Load TH3D MC histograms computed earlier ---
    # (This code block must already have created:)
    #   hGenRecoEvt = TH2D(pT,Z)    # from Project3D("zy")
    #   hGen        = TH2D(pT,Z)

    hGen.Sumw2()
    hGenRecoEvt.Sumw2()

    # =====================================================================
    # Match pT binning of rawMCSpectra, keep original Z binning
    # =====================================================================
    nPtBins = rawMCSpectra.GetNbinsX()
    ptAxis = rawMCSpectra.GetXaxis()

    # Build the variable pT bin array
    ptBins = [ptAxis.GetBinLowEdge(i + 1) for i in range(nPtBins + 1)]
    ptBinsArray = array.array('d', ptBins)

    # Original Z binning stays the same
    nZBins = hGenRecoEvt.GetNbinsY()
    zBinsArray = hGenRecoEvt.GetYaxis().GetXbins().GetArray()
    # If uniform, zBinsArray may be None → handle that:
    if not zBinsArray:
        zBinsArray = array.array(
            'd',
            [hGenRecoEvt.GetYaxis().GetBinLowEdge(i+1) for i in range(nZBins+1)]
        )

    # Create new histograms with desired X binning (pt) and original Y binning (z)
    hGenRecoEvt_rebinned = ROOT.TH2D(
        "hGenRecoEvt_rebinned",
        hGenRecoEvt.GetTitle(),
        nPtBins, ptBinsArray,
        nZBins, zBinsArray
    )

    hGen_rebinned = ROOT.TH2D(
        "hGen_rebinned",
        hGen.GetTitle(),
        nPtBins, ptBinsArray,
        nZBins, zBinsArray
    )

    # Enable storage of sum of squares of weights
    hGenRecoEvt_rebinned.Sumw2()
    hGen_rebinned.Sumw2()

    # Fill new histograms by mapping old bin centers
    for ix in range(1, hGenRecoEvt.GetNbinsX() + 1):
        for iy in range(1, nZBins + 1):

            # get old values
            content_evt = hGenRecoEvt.GetBinContent(ix, iy)
            err_evt = hGenRecoEvt.GetBinError(ix, iy)

            content_gen = hGen.GetBinContent(ix, iy)
            err_gen = hGen.GetBinError(ix, iy)

            # find coordinates
            x = hGenRecoEvt.GetXaxis().GetBinCenter(ix)
            y = hGenRecoEvt.GetYaxis().GetBinCenter(iy)

            # fill new histograms
            if content_evt != 0:
                hGenRecoEvt_rebinned.Fill(x, y, content_evt)

            if content_gen != 0:
                hGen_rebinned.Fill(x, y, content_gen)


    # =====================================================================
    # EFF × ACC (2D: pT vs Z)
    # =====================================================================
    # Calculate Efficiency x Acc x B.R
    # N. Reco in the selected reco colls / N. Gen V0s with at least one reco coll
    hEffAcc = rawMCSpectra.Clone("hEfficiencyXAcceptance2D")
    hEffAcc.Divide(hGenRecoEvt_rebinned)

    # =====================================================================
    # SIGNAL LOSS (2D)
    # =====================================================================
    # Calculate Signal loss
    # Signal Loss = N. Gen V0s with at least one reco coll / N. Gen V0s in ALL Generated colls

    hSignalLoss = hGenRecoEvt_rebinned.Clone("hSignalLoss2D")
    hSignalLoss.Divide(hGen_rebinned)

    # =====================================================================
    # EVENT LOSS (scalar)
    # =====================================================================
    # Calculate event loss
    # Event Loss = N. gen evts with least one reco coll (after evt sel) / N. of ALL gen evts
    hEventLoss = hGenRecoEvt_events.Clone("hEventLoss")
    hEventLoss.Divide(hGenEvent)

    # =====================================================================
    # EVENT SPLITTING (scalar)
    # =====================================================================
    # Calculate event splitting
    # Event Splitting = N. gen evts with least one reco coll (after evt sel) / N. of reco colls
    hEventSplitting = hGenRecoEvt_events.Clone("hEventSplitting")
    hEventSplitting.Divide(hRecoEvents)

    # =====================================================================
    # Write output file
    # =====================================================================

    fout = ROOT.TFile(output_path, "RECREATE")
    fout.cd()

    hGenAndRecoEvts.Write() # Generated and reco events
    hGen.Write() # Generated Lambda (all generated)
    hGenRecoEvt.Write() # Generated Lambda with at least one reco
    rawMCSpectra.Write() # Raw MC Spectra

    hEffAcc.Write()          # 2D pT×Z efficiency × acceptance # Efficiency x Acc x B.R. Output name is hEfficiencyXAcceptance2D
    hRecoEvents.Write()       # Reco events
    hSignalLoss.Write()      # 2D pT×Z signal loss. Output name is hSignalLoss2D
    hEventLoss.Write()       # scalar. Output name is hEventLoss
    hEventSplitting.Write()  # scalar. Output name is hEventSplitting

    hGenRecoEvt_rebinned.Write()
    hGen_rebinned.Write()

    hCentrality.Write()
    hGenEvent.Write()
    hGenRecoEvt_events.Write()

    fout.Close()
    fin.Close()

    print(f"Output written to: {output_path}")


# ------------------- Command-line interface --------------------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Calculate efficiency corrections for Lambda from AnalysisResults root file.")
    parser.add_argument("input_root", help="Path to MCAnalysisResults_*.root (single file).")
    parser.add_argument("--main-dir", default="Main", help="Main directory inside the ROOT file where Gen/ histograms live (default: 'Main').")
    parser.add_argument("--species", default="Lambda", choices=["Lambda", "ALambda", "LambdaALambda"],
                        help="Species to compute for: Lambda, ALambda (anti-Lambda) or LambdaALambda combined.")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    calculate_efficiency(args.input_root, main_dir=args.main_dir, species=args.species)