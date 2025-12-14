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

# # Actually need 3D projections now that we are working with:
# # h3dGenLambdaVsZ                        # (centrality, Z, pT) - all generated Lambdas
# # h3dGenLambdaVsZVsMultMC_RecoedEvt     # (Nch, Z, pT) - generated Lambda with >=1 reco track
# # h3dGenLambdaVsZVsMultMC               # (Nch, Z, pT) - all generated Lambda (same as above but maybe different binnings)
# def project_th3d(hist3D, axis, 
#                  event_range=None, z_range=None, pt_range=None,
#                  number=-1):
#     if not hist3D:
#         raise ValueError("Invalid TH3D provided.")

#     # Determine bin limits
#     e_min, e_max = 1, hist3D.GetNbinsX()
#     z_min, z_max = 1, hist3D.GetNbinsY()
#     p_min, p_max = 1, hist3D.GetNbinsZ()

#     if event_range:
#         e_min = hist3D.GetXaxis().FindBin(event_range[0])
#         e_max = hist3D.GetXaxis().FindBin(event_range[1])

#     if z_range:
#         z_min = hist3D.GetYaxis().FindBin(z_range[0])
#         z_max = hist3D.GetYaxis().FindBin(z_range[1])

#     if pt_range:
#         p_min = hist3D.GetZaxis().FindBin(pt_range[0])
#         p_max = hist3D.GetZaxis().FindBin(pt_range[1])

#     # Axis selection
#     if axis == "pt":
#         return hist3D.ProjectionZ(f"projZ_{number}", e_min, e_max, z_min, z_max)
#     if axis == "Z":
#         return hist3D.ProjectionY(f"projY_{number}", e_min, e_max, p_min, p_max)
#     if axis == "evt":
#         return hist3D.ProjectionX(f"projX_{number}", z_min, z_max, p_min, p_max)

#     raise ValueError("axis must be one of: 'pt','Z','evt'")

def rebinTH2D(src, dst):
    src.Sumw2()
    dst.Sumw2()

    for ix in range(1, src.GetNbinsX() + 1):
        x = src.GetXaxis().GetBinCenter(ix)

        for iy in range(1, src.GetNbinsY() + 1):
            y = src.GetYaxis().GetBinCenter(iy)

            content = src.GetBinContent(ix, iy)
            error   = src.GetBinError(ix, iy)

            if content == 0:
                continue

            # Find destination bin
            bin_dst = dst.FindBin(x, y)

            # Accumulate content
            dst.AddBinContent(bin_dst, content)

            # Proper error propagation
            old_err = dst.GetBinError(bin_dst)
            dst.SetBinError(
                bin_dst,
                (old_err**2 + error**2)**0.5
            )


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
    print("\nStarting main execution.")
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
            h3d_gen = f"{main_dir}/h3dGenAntiLambdaVsZVsMultMC"
            h3d_gen_recoed = f"{main_dir}/h3dGenAntiLambdaVsZVsMultMC_RecoedEvt"
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

    # === Load TH3D: Generated Λ with >=1 reco event ===
    h3dGenLambda_RecoedEvt = fin.Get(h3d_gen_recoed_name)
    if not h3dGenLambda_RecoedEvt:
        print(f"Warning: could not find {h3d_gen_recoed_name} in file.")
    else:
        # hGenLambdaRecoEvt = h3dGenLambda_RecoedEvt.Project3D("zy") # This was giving the wrong ordering!
        # For a TH3D( X=event , Y=Z , Z=pT ), the projection below will give (X = pT , Y = Z)
        hGenLambdaRecoEvtPtVsZ = h3dGenLambda_RecoedEvt.Project3D("yz")
        hGenLambdaRecoEvtPtVsZ.SetTitle(f"Gen{species}_RecoEvts")

    # === Load TH3D: All generated Λ ===
    h3dGenLambda = fin.Get(h3d_gen_name)
    if not h3dGenLambda:
        print(f"Warning: could not find {h3d_gen_name} in file.")
    else:
        hGenLambdaPtVsZ = h3dGenLambda.Project3D("yz") # Same fix as previous block
        hGenLambdaPtVsZ.SetTitle(f"Gen{species}_GenEvts")

    hGenLambdaRecoEvtPtVsZ.Sumw2()
    hGenLambdaPtVsZ.Sumw2()

    # =====================================================================
    # N. Reco in the selected reco colls :
    # =====================================================================
    # raw MC spectra: TH2D(pT,Z)
    rawMCSpectraPtVsZ = fin.Get(f"{main_dir}/Lambda/hLambdaPtZ") # --> This is the guy you will need to change for the signal extraction code!!!
    if not rawMCSpectraPtVsZ:
        raise RuntimeError(f"Could not find '{main_dir}/Lambda/hLambdaPtZ' in the input file.")
    rawMCSpectraPtVsZ.Sumw2()

        # Rebinning other histograms (hGenLambdaRecoEvtPtVsZ) to match this rawMCSpectraPtVsZ:
    # First, get the desired geometry via a clone + error elimination:
    hGenLambdaRecoEvtPtVsZ_rebinned = rawMCSpectraPtVsZ.Clone("hGenLambdaRecoEvtPtVsZ_rebinned")
    hGenLambdaRecoEvtPtVsZ_rebinned.Reset("ICES")  # reset contents + errors
    hGenLambdaRecoEvtPtVsZ_rebinned.Sumw2()
    rebinTH2D(src = hGenLambdaRecoEvtPtVsZ, dst = hGenLambdaRecoEvtPtVsZ_rebinned)

    # Same for hGenLambdaPtVsZ:
    hGenLambdaPtVsZ_rebinned = rawMCSpectraPtVsZ.Clone("hGenLambdaPtVsZ_rebinned")
    hGenLambdaPtVsZ_rebinned.Reset("ICES")  # reset contents + errors
    hGenLambdaPtVsZ_rebinned.Sumw2()
    rebinTH2D(src = hGenLambdaPtVsZ, dst = hGenLambdaPtVsZ_rebinned)

    # =====================================================================
    # Event histograms
    # =====================================================================

    h2dGenAndRecoEvents = fin.Get(f"{main_dir}/hGenEvents")
    if not h2dGenAndRecoEvents:
        raise RuntimeError(f"Could not find {main_dir}/hGenEvents in input file.")

    # For event counts we still want 1D, even though MC histos are 3D.
    hGenAndRecoEvts = project_th2d(h2dGenAndRecoEvents, "Y", x_range, y_range, rebin_factor)

    ## All Gen evts
    hGenEvent = ROOT.TH1D("hGenEvent", "hGenEvent", 1, 0, 1)    
    hGenEvent.SetBinContent(1, hGenAndRecoEvts.GetBinContent(1))
    hGenEvent.SetBinError(1, hGenAndRecoEvts.GetBinError(1))

    ## Gen evt with at least one reco coll (after evt sel)
    hGenRecoEvt = ROOT.TH1D("hGenRecoEvt", "hGenRecoEvt", 1, 0, 1)    
    hGenRecoEvt.SetBinContent(1, hGenAndRecoEvts.GetBinContent(2))
    hGenRecoEvt.SetBinError(1, hGenAndRecoEvts.GetBinError(2))

    print("Gen Events:", hGenAndRecoEvts.GetBinContent(1))
    print("Gen Reco Events:", hGenAndRecoEvts.GetBinContent(2))
    print("Evt Loss:", hGenAndRecoEvts.GetBinContent(2)/hGenAndRecoEvts.GetBinContent(1))

    # Propagate errors:
    hGenRecoEvt.Sumw2()

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
    hRecoEvents.SetBinContent(1, hCentrality.Integral()) # Number of reconstructed collisions!
    hRecoEvents.SetBinError(1, 0)
    hRecoEvents.Sumw2()
    print("hCentrality.GetBinContent(1): ", hCentrality.GetBinContent(1))
    print("hCentrality.Integral()", hCentrality.Integral())

    # =====================================================================
    # PROPAGATE ERRORS
    # =====================================================================

    hGenEvent.Sumw2()
    hGenAndRecoEvts.Sumw2()
    hCentrality.Sumw2()
    hRecoEvents.Sumw2()

    # =====================================================================
    # EFF × ACC (2D: pT vs Z)
    # =====================================================================
    # Calculate Efficiency x Acc x B.R
    # N. Reco in the selected reco colls / N. Gen V0s with at least one reco coll
    hEffAcc = rawMCSpectraPtVsZ.Clone("hEfficiencyXAcceptance2D")
    # An explicitly binomial error propagation -- These two histograms are of course correlated!!!
    # The error bars add differently!
    hEffAcc.Divide(rawMCSpectraPtVsZ, hGenLambdaRecoEvtPtVsZ_rebinned, 1, 1, "B")

    # =====================================================================
    # SIGNAL LOSS (2D)
    # =====================================================================
    # Calculate Signal loss
    # Signal Loss = N. Gen V0s with at least one reco coll / N. Gen V0s in ALL Generated colls

    hSignalLoss2D = hGenLambdaRecoEvtPtVsZ_rebinned.Clone("hSignalLoss2D")
    # Same correction as above -- binomial error bars!
    hSignalLoss2D.Divide(hGenLambdaRecoEvtPtVsZ_rebinned, hGenLambdaPtVsZ_rebinned, 1, 1, "B")

    # =====================================================================
    # EVENT LOSS (scalar)
    # =====================================================================
    # Calculate event loss
    # Event Loss = N. gen evts with least one reco coll (after evt sel) / N. of ALL gen evts
    hEventLoss = hGenRecoEvt.Clone(f"hEventLoss")
    hEventLoss.Divide(hGenEvent)

    # =====================================================================
    # EVENT SPLITTING (scalar)
    # =====================================================================
    # Calculate event splitting
    # Event Splitting = N. gen evts with least one reco coll (after evt sel) / N. of reco colls
    hEventSplitting = hGenRecoEvt.Clone(f"hEventSplitting")
    hEventSplitting.Divide(hRecoEvents)

    # =====================================================================
    # 1D Projections: Efficiency & Signal Loss (Z-integrated, Z>0, Z<0)
    # =====================================================================
    print("Computing 1D Efficiency and Signal Loss for Z ranges (Pos, Neg, Full)...")

    def get_z_projections(h2d, hist_name_base):
        """
        Helper to project 2D (pT, Z) into 1D (pT) for Full, Pos (Z>=0), and Neg (Z<0).
        """
        yaxis = h2d.GetYaxis()
        n_bins_y = h2d.GetNbinsY()
        
        # Find the bin index corresponding to Z=0
        bin_zero = yaxis.FindBin(0.0)
        bin_center = yaxis.GetBinCenter(bin_zero)
        
        # Determine split indices based on where 0 falls
        if bin_center < 0:
            bin_neg_end = bin_zero
            bin_pos_start = bin_zero + 1
        else:
            # If center >= 0, include this bin in Positive
            bin_neg_end = bin_zero - 1
            bin_pos_start = bin_zero
            
        # 1. Full Projection (Bin 1 to N)
        h1d_full = h2d.ProjectionX(f"{hist_name_base}_Full", 1, n_bins_y, "e")
        
        # 2. Negative Z Projection
        if bin_neg_end >= 1:
            h1d_neg = h2d.ProjectionX(f"{hist_name_base}_Neg", 1, bin_neg_end, "e")
        else:
            # Handle edge case where no bins are negative
            h1d_neg = h2d.ProjectionX(f"{hist_name_base}_Neg", 0, 0, "e") 
            h1d_neg.Reset()

        # 3. Positive Z Projection
        if bin_pos_start <= n_bins_y:
            h1d_pos = h2d.ProjectionX(f"{hist_name_base}_Pos", bin_pos_start, n_bins_y, "e")
        else:
            h1d_pos = h2d.ProjectionX(f"{hist_name_base}_Pos", 0, 0, "e")
            h1d_pos.Reset()
            
        return h1d_full, h1d_pos, h1d_neg

    # --- Step A: Project Numerators and Denominators ---

    # 1. Reco (Numerator for Efficiency) -> from rawMCSpectraPtVsZ
    rec_full, rec_pos, rec_neg = get_z_projections(rawMCSpectraPtVsZ, "hRec1D")
    
    # 2. GenReco (Denominator for Efficiency / Numerator for Signal Loss) -> from hGenLambdaRecoEvtPtVsZ_rebinned
    genrec_full, genrec_pos, genrec_neg = get_z_projections(hGenLambdaRecoEvtPtVsZ_rebinned, "hGenRec1D")
    
    # 3. GenAll (Denominator for Signal Loss) -> from hGenLambdaPtVsZ_rebinned
    genall_full, genall_pos, genall_neg = get_z_projections(hGenLambdaPtVsZ_rebinned, "hGenAll1D")

    # --- Step B: Calculate Ratios (Efficiency & Signal Loss) with "B" Option ---

    # Efficiency * Acceptance = Reco / GenReco
    hEff_Full = rec_full.Clone("hEfficiency1D_Full")
    hEff_Full.SetTitle("Efficiency x Acceptance (Full Z)")
    hEff_Full.Divide(rec_full, genrec_full, 1, 1, "B") # <--- Added "B" for binomial errors

    hEff_Pos = rec_pos.Clone("hEfficiency1D_Pos")
    hEff_Pos.SetTitle("Efficiency x Acceptance (Z >= 0)")
    hEff_Pos.Divide(rec_pos, genrec_pos, 1, 1, "B")    # <--- Added "B"

    hEff_Neg = rec_neg.Clone("hEfficiency1D_Neg")
    hEff_Neg.SetTitle("Efficiency x Acceptance (Z < 0)")
    hEff_Neg.Divide(rec_neg, genrec_neg, 1, 1, "B")    # <--- Added "B"

    # Signal Loss = GenReco / GenAll
    hSigLoss_Full = genrec_full.Clone("hSignalLoss1D_Full")
    hSigLoss_Full.SetTitle("Signal Loss (Full Z)")
    hSigLoss_Full.Divide(genrec_full, genall_full, 1, 1, "B") # <--- Added "B"

    hSigLoss_Pos = genrec_pos.Clone("hSignalLoss1D_Pos")
    hSigLoss_Pos.SetTitle("Signal Loss (Z > 0)")
    hSigLoss_Pos.Divide(genrec_pos, genall_pos, 1, 1, "B")    # <--- Added "B"

    hSigLoss_Neg = genrec_neg.Clone("hSignalLoss1D_Neg")
    hSigLoss_Neg.SetTitle("Signal Loss (Z < 0)")
    hSigLoss_Neg.Divide(genrec_neg, genall_neg, 1, 1, "B")    # <--- Added "B"

    # =====================================================================
    # Write output file
    # =====================================================================

    fout = ROOT.TFile(output_path, "RECREATE")
    fout.cd()

    hGenAndRecoEvts.Write() # Generated and reco events
    hGenLambdaPtVsZ.Write() # Generated Lambda (all generated)
    hGenLambdaRecoEvtPtVsZ.Write() # Generated Lambda with at least one reco
    rawMCSpectraPtVsZ.Write() # Raw MC Spectra

    hEffAcc.Write()          # 2D pT×Z efficiency × acceptance # Efficiency x Acc x B.R. Output name is hEfficiencyXAcceptance2D
    hRecoEvents.Write()       # Reco events
    hSignalLoss2D.Write()      # 2D pT×Z signal loss. Output name is hSignalLoss2D
    hEventLoss.Write()       # scalar. Output name is hEventLoss
    hEventSplitting.Write()  # scalar. Output name is hEventSplitting

    # Rebins:
    hGenLambdaRecoEvtPtVsZ_rebinned.Write()
    hGenLambdaPtVsZ_rebinned.Write()

    hCentrality.Write()
    hGenEvent.Write()
    hGenRecoEvt.Write()

    # 1D Projections:
    hEff_Full.Write()
    hEff_Pos.Write()
    hEff_Neg.Write()
    
    hSigLoss_Full.Write()
    hSigLoss_Pos.Write()
    hSigLoss_Neg.Write()

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