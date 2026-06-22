// ============================================================================
// compareHelicityFamilies.cxx
// ============================================================================
//
// PURPOSE
// -------
// Post-generator cross-family and within-family comparison script for the
// helicity efficiency toy model study.  Reads every ROOT file produced by
// the generator (via runHelicityToyModel.sh) and writes a single output ROOT
// file containing TCanvas objects organized in a clean directory tree.
//
// Three comparison types are produced, all using the SEM error estimator
// (TProfile SEM as primary markers; Kahan SEM as small diamond overlays):
//
//   (1) WITHIN-FAMILY  -- For each family, x-axis = generator configurations
//       in registration order (typically increasing cut value), 4 scenario
//       series (NoCuts / pTCutOnly / DCACutOnly / BothCuts) as colored markers.
//       Also produces a NET-EFFECT canvas: (scenario - NoCuts) for each config,
//       isolating the fake polarization introduced by the cuts.
//
//   (2) CROSS-FAMILY (representatives) -- One entry per family on the x-axis,
//       using a pre-defined representative configuration that reflects the
//       typical experimental working point.  See repIdx comments in BuildFamilies.
//       Also produces a BothCuts-only version for easier reading.
//
//   (3) BRUTE-FORCE -- Every configuration from every family in family order,
//       all on one wide canvas.  Vertical lines + family labels separate groups.
//       Also produces a BothCuts-only version.
//
// OUTPUT DIRECTORY STRUCTURE (inside compare_families.root)
// ---------------------------------------------------------
//   Family_N_<ShortName>/
//     WithEtaGate/    WithoutEtaGate/
//       EtaLambda/    EtaJet/      <-- within-family (all 4 scenarios)
//       NetEffect/
//         EtaLambda/  EtaJet/      <-- scenario - NoCuts
//   CrossFamily/
//     WithEtaGate/    WithoutEtaGate/
//       EtaLambda/    EtaJet/      <-- representative per family, all 4 scenarios
//       BothCutsOnly/
//         EtaLambda/  EtaJet/      <-- BothCuts only, cleaner view
//       NetEffect/
//         EtaLambda/  EtaJet/      <-- (scenario - NoCuts) for representatives
//   BruteForce/
//     WithEtaGate/    WithoutEtaGate/
//       EtaLambda/    EtaJet/      <-- all configs, all families, all scenarios
//       BothCutsOnly/
//         EtaLambda/  EtaJet/
//       NetEffect/
//         EtaLambda/  EtaJet/
//
// HISTOGRAM SOURCES
// -----------------
// EtaLambda selections read pRingProxyJet / hRingProxyJet_Kahan from the
//   etaSel sub-directory (EtaPos | EtaNeg | All) of each scenario.
// EtaJet selections read pRingProxyJet_JetEtaPos(Neg) / their Kahan
//   counterparts from the "All" sub-directory -- those conditional TProfiles
//   are only stored in the "All" Lambda-eta directory.
//
// USAGE
// -----
//   root -l -b -q 'compareHelicityFamilies.cxx'
//       (uses the default BASE_DIR hardcoded near the bottom)
//
//   root -l -b -q 'compareHelicityFamilies.cxx("/path/to/BaseDir")'
//   root -l -b -q 'compareHelicityFamilies.cxx("/path/to/BaseDir","out.root")'
//
//   Called automatically from runHelicityToyModel.sh after all per-run jobs
//   finish.  Missing ROOT files are silently skipped -- it is safe to call
//   this after a partial family run (e.g., --family 9).
//
// ============================================================================

#include "TFile.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystem.h"
#include "TROOT.h"

#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdio>


// ============================================================================
// STYLE CONSTANTS  (mirror plotHelicityEfficiency.cxx for visual consistency)
// ============================================================================

static const int    kColNoCuts = kBlack;
static const int    kColPtCut  = kRed + 1;
static const int    kColDcaCut = kBlue + 1;
static const int    kColBoth   = kGreen + 2;
static const int    kLineW     = 2;
static const double kMkSz      = 0.95;

static const char* kScenNames [4] = {"NoCuts","pTCutOnly","DCACutOnly","BothCuts"};
static const char* kScenLabels[4] = {"No cuts","p_{T} cut only","DCA cut only","Both cuts"};
static const int   kScenCols  [4] = {kColNoCuts, kColPtCut, kColDcaCut, kColBoth};
static const int   kScenMkF   [4] = {20, 21, 22, 23};  // filled: circle/square/triU/triD
static const int   kScenMkO   [4] = {24, 25, 26, 32};  // open counterparts

static const char* kEtaGates [2] = {"WithEtaGate","WithoutEtaGate"};
static const char* kEtaGLbls [2] = {"w/ #eta gate","w/o #eta gate"};

// ============================================================================
// RING SOURCE DESCRIPTORS  (5 integrated ring selections)
// ============================================================================

struct RingSource {
    const char* etaSel;    ///< Sub-directory under each scenario (for EtaLambda srcs)
    const char* profName;  ///< TProfile histogram name
    const char* kahanName; ///< Kahan TH1D histogram name
    const char* axisLbl;   ///< TLatex string for canvas annotations
    const char* subDir;    ///< Output sub-directory ("EtaLambda" or "EtaJet")
    const char* sufx;      ///< Short suffix for canvas/directory names
};

// Sources 0-2: Lambda-eta selections.  Sources 3-4: jet-eta conditional selections.
// JetEtaPos/Neg live in the "All" Lambda-eta directory; etaSel is set to "All"
// but profName / kahanName differ.
static const RingSource kSrcs[5] = {
    {"EtaPos", "pRingProxyJet",         "hRingProxyJet_Kahan",
     "#eta_{#Lambda} > 0",  "EtaLambda", "EtaPos"},
    {"EtaNeg", "pRingProxyJet",         "hRingProxyJet_Kahan",
     "#eta_{#Lambda} < 0",  "EtaLambda", "EtaNeg"},
    {"All",    "pRingProxyJet",         "hRingProxyJet_Kahan",
     "All #eta_{#Lambda}",  "EtaLambda", "All"},
    {"All",    "pRingProxyJet_JetEtaPos","hRingProxyJet_JetEtaPos_Kahan",
     "#eta_{jet} #geq 0",   "EtaJet",    "JetEtaPos"},
    {"All",    "pRingProxyJet_JetEtaNeg","hRingProxyJet_JetEtaNeg_Kahan",
     "#eta_{jet} < 0",      "EtaJet",    "JetEtaNeg"},
};
static const int kNSrc = 5;


// ============================================================================
// DATA STRUCTURES
// ============================================================================

struct GenConfig {
    std::string name;    ///< generator run name, e.g. "dca_sym_005"
    std::string subdir;  ///< relative path under BASE_DIR, e.g. "2_SymDCA/dca_sym_005"
    std::string label;   ///< short ASCII label used as x-axis bin label
};

struct FamilyInfo {
    int         num;
    std::string shortName;   ///< "SymDCA" -- used in directory/canvas names
    std::string longName;    ///< human-readable description
    std::vector<GenConfig> cfgs; ///< in display order (bash registration / increasing cut)
    int         repIdx;      ///< index into cfgs for the "representative" config
};

/// All integrated ring values pre-loaded from one ROOT file.
/// Indices: [etaGate 0..1][scenario 0..3][source 0..4]
struct RunVals {
    double val  [2][4][5];  ///< TProfile mean (SEM)
    double err  [2][4][5];  ///< TProfile SEM
    double kval [2][4][5];  ///< Kahan mean
    double kerr [2][4][5];  ///< Kahan SEM
    bool   ok   [2][4][5];  ///< whether the TProfile was found and non-empty
    bool   kok  [2][4][5];  ///< whether the Kahan TH1D was found
    bool   found;           ///< whether the ROOT file existed and was readable

    RunVals() : found(false) {
        memset(val,   0, sizeof(val));   memset(err,  0, sizeof(err));
        memset(kval,  0, sizeof(kval));  memset(kerr, 0, sizeof(kerr));
        memset(ok,    0, sizeof(ok));    memset(kok,  0, sizeof(kok));
    }
};


// ============================================================================
// FAMILY REGISTRY
// Encode the exact same families/configs/order as runHelicityToyModel.sh.
// repIdx points to the config best representing the "standard" working point.
// ============================================================================

static std::vector<FamilyInfo> BuildFamilies()
{
    std::vector<FamilyInfo> F;
    F.reserve(13);

    // -- F0: Baseline (no cuts) -----------------------------------------------
    { FamilyInfo f; f.num=0; f.repIdx=0;
      f.shortName="Baseline"; f.longName="Baseline -- no cuts";
      f.cfgs={{"baseline","0_Baseline/baseline","base"}};
      F.push_back(std::move(f)); }

    // -- F1: Asymmetric DCA cuts [AEE daughter-level asymmetry] ---------------
    // repIdx=3 --> dca_asym_std: the standard ALICE-like (0.05/0.10 cm) asymmetric config
    { FamilyInfo f; f.num=1; f.repIdx=3;
      f.shortName="AsymDCA"; f.longName="Asymmetric DCA cuts [AEE probe]";
      f.cfgs={
        {"dca_p_only_005",  "1_AsymDCA/dca_p_only_005",  "p:0.05"},
        {"dca_pi_only_005", "1_AsymDCA/dca_pi_only_005", "pi:0.05"},
        {"dca_pi_only_010", "1_AsymDCA/dca_pi_only_010", "pi:0.10"},
        {"dca_asym_std",    "1_AsymDCA/dca_asym_std",    "std"},
        {"dca_asym_wide",   "1_AsymDCA/dca_asym_wide",   "wide"},
        {"dca_asym_rev",    "1_AsymDCA/dca_asym_rev",    "rev"},
      };
      F.push_back(std::move(f)); }

    // -- F2: Symmetric DCA cuts [AEE strength scan] ---------------------------
    // repIdx=1 --> dca_sym_010: a clean, single-parameter reference
    { FamilyInfo f; f.num=2; f.repIdx=1;
      f.shortName="SymDCA"; f.longName="Symmetric DCA cuts [AEE strength]";
      f.cfgs={
        {"dca_sym_005","2_SymDCA/dca_sym_005","0.05cm"},
        {"dca_sym_010","2_SymDCA/dca_sym_010","0.10cm"},
        {"dca_sym_015","2_SymDCA/dca_sym_015","0.15cm"},
        {"dca_sym_020","2_SymDCA/dca_sym_020","0.20cm"},
        {"dca_sym_030","2_SymDCA/dca_sym_030","0.30cm"},
      };
      F.push_back(std::move(f)); }

    // -- F3: Magnetic field scan, DCA only [AEE field dependence] -------------
    // repIdx=4 --> field_b050: standard B = +0.5 T
    { FamilyInfo f; f.num=3; f.repIdx=4;
      f.shortName="BField"; f.longName="B-field scan (DCA only) [AEE field dep.]";
      f.cfgs={
        {"field_vlow",    "3_BField/field_vlow",    "B=0.0001T"},
        {"field_b010",    "3_BField/field_b010",    "B=0.10T"},
        {"field_b020",    "3_BField/field_b020",    "B=0.20T"},
        {"field_b030",    "3_BField/field_b030",    "B=0.30T"},
        {"field_b050",    "3_BField/field_b050",    "B=0.50T"},
        {"field_b075",    "3_BField/field_b075",    "B=0.75T"},
        {"field_b100",    "3_BField/field_b100",    "B=1.00T"},
        {"field_b050neg", "3_BField/field_b050neg", "B=-0.50T"},
      };
      F.push_back(std::move(f)); }

    // -- F4: Lambda pT minimum [kinematic context, standard ALICE cuts] --------
    // repIdx=0 --> lam_ptmin_000: fully inclusive spectrum
    { FamilyInfo f; f.num=4; f.repIdx=0;
      f.shortName="LamPtMin"; f.longName="Lambda pT minimum [kinematic context]";
      f.cfgs={
        {"lam_ptmin_000","4_LamPtMin/lam_ptmin_000","pT>0.0"},
        {"lam_ptmin_030","4_LamPtMin/lam_ptmin_030","pT>0.3"},
        {"lam_ptmin_050","4_LamPtMin/lam_ptmin_050","pT>0.5"},
        {"lam_ptmin_100","4_LamPtMin/lam_ptmin_100","pT>1.0"},
        {"lam_ptmin_150","4_LamPtMin/lam_ptmin_150","pT>1.5"},
        {"lam_ptmin_200","4_LamPtMin/lam_ptmin_200","pT>2.0"},
      };
      F.push_back(std::move(f)); }

    // -- F5: Daughter pT cuts [HEE primary probe] ------------------------------
    // repIdx=2 --> pt_sym_150: the ALICE standard-ish symmetric cut
    { FamilyInfo f; f.num=5; f.repIdx=2;
      f.shortName="DaughterPt"; f.longName="Daughter pT cuts [HEE probe]";
      f.cfgs={
        {"pt_sym_005","5_DaughterPt/pt_sym_005","sym5"},
        {"pt_sym_010","5_DaughterPt/pt_sym_010","sym10"},
        {"pt_sym_150","5_DaughterPt/pt_sym_150","sym15"},
        {"pt_sym_200","5_DaughterPt/pt_sym_200","sym20"},
        {"pt_sym_300","5_DaughterPt/pt_sym_300","sym30"},
        {"pt_pi_010", "5_DaughterPt/pt_pi_010", "pi10"},
        {"pt_pi_150", "5_DaughterPt/pt_pi_150", "pi15"},
        {"pt_pi_200", "5_DaughterPt/pt_pi_200", "pi20"},
        {"pt_p_150",  "5_DaughterPt/pt_p_150",  "p15"},
      };
      F.push_back(std::move(f)); }

    // -- F6: Eta acceptance window [HEE+AEE vs detector acceptance] -----------
    // repIdx=2 --> eta_090: the ALICE inner barrel
    { FamilyInfo f; f.num=6; f.repIdx=2;
      f.shortName="EtaSym"; f.longName="Eta acceptance window [HEE+AEE vs acceptance]";
      f.cfgs={
        {"eta_050","6_EtaSym/eta_050","|eta|<0.50"},
        {"eta_070","6_EtaSym/eta_070","|eta|<0.70"},
        {"eta_090","6_EtaSym/eta_090","|eta|<0.90"},
        {"eta_120","6_EtaSym/eta_120","|eta|<1.20"},
        {"eta_150","6_EtaSym/eta_150","|eta|<1.50"},
        {"eta_200","6_EtaSym/eta_200","|eta|<2.00"},
        {"eta_300","6_EtaSym/eta_300","|eta|<3.00"},
      };
      F.push_back(std::move(f)); }

    // -- F7: Boltzmann temperature [spectrum shape] ---------------------------
    // repIdx=2 --> temp_030: T = 0.30 GeV, nominal LHC Lambda inverse slope
    { FamilyInfo f; f.num=7; f.repIdx=2;
      f.shortName="Temperature"; f.longName="Boltzmann temperature [spectrum shape]";
      f.cfgs={
        {"temp_020","7_Temperature/temp_020","T=0.20"},
        {"temp_025","7_Temperature/temp_025","T=0.25"},
        {"temp_030","7_Temperature/temp_030","T=0.30"},
        {"temp_035","7_Temperature/temp_035","T=0.35"},
        {"temp_040","7_Temperature/temp_040","T=0.40"},
        {"temp_050","7_Temperature/temp_050","T=0.50"},
      };
      F.push_back(std::move(f)); }

    // -- F8: Ring kinematic pT windows ----------------------------------------
    // repIdx=3 --> win_ring: the ring analysis window [0.5, 1.5] GeV/c
    { FamilyInfo f; f.num=8; f.repIdx=3;
      f.shortName="KinWindow"; f.longName="Lambda pT windows [fake signal per band]";
      f.cfgs={
        {"win_inclusive","8_KinWindow/win_inclusive","incl"},
        {"win_vsoft",    "8_KinWindow/win_vsoft",    "vsoft"},
        {"win_soft",     "8_KinWindow/win_soft",     "soft"},
        {"win_ring",     "8_KinWindow/win_ring",     "ring"},
        {"win_mid",      "8_KinWindow/win_mid",      "mid"},
        {"win_hard",     "8_KinWindow/win_hard",     "hard"},
        {"win_vhard",    "8_KinWindow/win_vhard",    "vhard"},
      };
      F.push_back(std::move(f)); }

    // -- F9: Realistic ALICE cuts [combined HEE+AEE estimate] -----------------
    // repIdx=1 --> alice_std: the standard working point (the one most relevant to data)
    { FamilyInfo f; f.num=9; f.repIdx=1;
      f.shortName="RealisticAlice"; f.longName="Realistic ALICE cuts [combined HEE+AEE]";
      f.cfgs={
        {"alice_loose",   "9_RealisticAlice/alice_loose",   "loose"},
        {"alice_std",     "9_RealisticAlice/alice_std",     "std"},
        {"alice_std_neg", "9_RealisticAlice/alice_std_neg", "std/neg"},
        {"alice_tight",   "9_RealisticAlice/alice_tight",   "tight"},
        {"alice_ring",    "9_RealisticAlice/alice_ring",    "ring"},
      };
      F.push_back(std::move(f)); }

    // -- F10: B-field with pT+DCA cuts [combined AEE+HEE field dep.] ----------
    // repIdx=1 --> fieldCuts_b050: standard +0.5 T with both cut types active
    { FamilyInfo f; f.num=10; f.repIdx=1;
      f.shortName="BFieldCuts"; f.longName="B-field scan (pT+DCA cuts) [AEE+HEE field dep.]";
      f.cfgs={
        {"fieldCuts_vlow",    "10_BField/fieldCuts_vlow",    "B=0.0001T"},
        {"fieldCuts_b050",    "10_BField/fieldCuts_b050",    "B=0.50T"},
        {"fieldCuts_b100",    "10_BField/fieldCuts_b100",    "B=1.00T"},
        {"fieldCuts_b050neg", "10_BField/fieldCuts_b050neg", "B=-0.50T"},
      };
      F.push_back(std::move(f)); }

    // -- F11: Asymmetric eta acceptance [Zvtx shift / asymmetric window] ------
    // repIdx=0 --> pos10cmZvtx: simplest emulated Zvtx shift; clear baseline
    { FamilyInfo f; f.num=11; f.repIdx=0;
      f.shortName="EtaAsym"; f.longName="Asymmetric eta acceptance [Zvtx shift]";
      f.cfgs={
        {"pos10cmZvtx",        "11_EtaAsym/pos10cmZvtx",        "+10cmZvtx"},
        {"neg10cmZvtx",        "11_EtaAsym/neg10cmZvtx",        "-10cmZvtx"},
        {"pos10cmZvtxNegField","11_EtaAsym/pos10cmZvtxNegField", "+10cm/B-"},
        {"EtaMin130etaMax090", "11_EtaAsym/EtaMin130etaMax090",  "[-1.3,0.9]"},
        {"EtaMin110etaMax090", "11_EtaAsym/EtaMin110etaMax090",  "[-1.1,0.9]"},
        {"EtaMin070etaMax090", "11_EtaAsym/EtaMin070etaMax090",  "[-0.7,0.9]"},
        {"EtaMin050etaMax090", "11_EtaAsym/EtaMin050etaMax090",  "[-0.5,0.9]"},
        {"EtaMin030etaMax090", "11_EtaAsym/EtaMin030etaMax090",  "[-0.3,0.9]"},
      };
      F.push_back(std::move(f)); }

    // -- F12: Lambda max rapidity at generator level --------------------------
    // repIdx=3 --> yLam4: the standard rapidity range used throughout the study
    { FamilyInfo f; f.num=12; f.repIdx=3;
      f.shortName="LambdaRap"; f.longName="Lambda max rapidity [generator acceptance]";
      f.cfgs={
        {"yLam0_5GuardTest","12_LambdaRap/yLam0_5GuardTest","y<0.5"},
        {"yLam2",           "12_LambdaRap/yLam2",           "y<2"},
        {"yLam3",           "12_LambdaRap/yLam3",           "y<3"},
        {"yLam4",           "12_LambdaRap/yLam4",           "y<4"},
        {"yLam5",           "12_LambdaRap/yLam5",           "y<5"},
        {"yLam6",           "12_LambdaRap/yLam6",           "y<6"},
        {"yLam10",          "12_LambdaRap/yLam10",          "y<10"},
      };
      F.push_back(std::move(f)); }

    return F;
}


// ============================================================================
// SHARED HELPERS
// ============================================================================

static void SetupStyle()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetTitleSize(0.05, "XYZ");
    gStyle->SetLabelSize(0.045, "XYZ");
    gStyle->SetTitleOffset(1.2, "Y");
}

static void AddLabel(double x, double y, const char* txt, double sz=0.032, int align=22)
{
    TLatex* l = new TLatex(x, y, txt);
    l->SetNDC(); l->SetTextSize(sz); l->SetTextAlign(align);
    l->Draw();
}

static TLegend* MakeLeg(double x1, double y1, double x2, double y2)
{
    TLegend* l = new TLegend(x1, y1, x2, y2);
    l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.030);
    return l;
}

/// Write canvas to directory, then delete it.
static void SaveTo(TCanvas* c, TDirectory* d)
{
    if (!c || !d) return;
    d->cd(); c->Write(c->GetName(), TObject::kOverwrite);
    delete c;
}

/// Get-or-create a named sub-directory.
static TDirectory* Mkdir(TDirectory* parent, const char* name)
{
    if (!parent) return nullptr;
    TDirectory* d = dynamic_cast<TDirectory*>(parent->Get(name));
    if (!d) d = parent->mkdir(name);
    return d;
}

// Global counter for unique canvas names (avoids ROOT's duplicate-name warning)
static int gCI = 0;


// ============================================================================
// DATA LOADING
// ============================================================================

/// Open one generator ROOT file and extract all integrated ring values.
/// Missing or unreadable files return an empty RunVals with found=false.
static RunVals LoadRunVals(const std::string& baseDir, const GenConfig& cfg)
{
    RunVals rv;
    std::string path = baseDir + "/" + cfg.subdir + "/helicity_" + cfg.name + ".root";
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) return rv;  // silently skip
    rv.found = true;

    for (int ieg = 0; ieg < 2; ++ieg) {
        TDirectory* egDir = dynamic_cast<TDirectory*>(f->Get(kEtaGates[ieg]));
        if (!egDir) continue;

        for (int is = 0; is < 4; ++is) {
            TDirectory* dScen = dynamic_cast<TDirectory*>(egDir->Get(kScenNames[is]));
            if (!dScen) continue;

            // Sources 0-2: EtaLambda selections (EtaPos, EtaNeg, All)
            for (int isrc = 0; isrc < 3; ++isrc) {
                TDirectory* d = dynamic_cast<TDirectory*>(dScen->Get(kSrcs[isrc].etaSel));
                if (!d) continue;
                TProfile* p = dynamic_cast<TProfile*>(d->Get(kSrcs[isrc].profName));
                if (p && p->GetEntries() > 0) {
                    rv.val[ieg][is][isrc] = p->GetBinContent(1);
                    rv.err[ieg][is][isrc] = p->GetBinError(1);
                    rv.ok [ieg][is][isrc] = true;
                }
                TH1D* hK = dynamic_cast<TH1D*>(d->Get(kSrcs[isrc].kahanName));
                if (hK) {
                    rv.kval[ieg][is][isrc] = hK->GetBinContent(1);
                    rv.kerr[ieg][is][isrc] = hK->GetBinError(1);
                    rv.kok [ieg][is][isrc] = true;
                }
            }

            // Sources 3-4: EtaJet selections (always from "All" Lambda-eta directory)
            TDirectory* dAll = dynamic_cast<TDirectory*>(dScen->Get("All"));
            if (dAll) {
                for (int isrc = 3; isrc < 5; ++isrc) {
                    TProfile* p = dynamic_cast<TProfile*>(dAll->Get(kSrcs[isrc].profName));
                    if (p && p->GetEntries() > 0) {
                        rv.val[ieg][is][isrc] = p->GetBinContent(1);
                        rv.err[ieg][is][isrc] = p->GetBinError(1);
                        rv.ok [ieg][is][isrc] = true;
                    }
                    TH1D* hK = dynamic_cast<TH1D*>(dAll->Get(kSrcs[isrc].kahanName));
                    if (hK) {
                        rv.kval[ieg][is][isrc] = hK->GetBinContent(1);
                        rv.kerr[ieg][is][isrc] = hK->GetBinError(1);
                        rv.kok [ieg][is][isrc] = true;
                    }
                }
            }
        }
    }
    f->Close();
    return rv;
}

/// Build a "net-effect" RunVals vector: every scenario minus a reference scenario.
/// The reference scenario (typically NoCuts = index 0) is subtracted from all others.
/// Errors are propagated in quadrature (conservative: assumes statistical independence,
/// which over-estimates the error since NoCuts and cut scenarios share the same events).
static std::vector<RunVals> ComputeNetEffect(const std::vector<RunVals>& data, int ref = 0)
{
    std::vector<RunVals> out(data.size());
    for (int ic = 0; ic < (int)data.size(); ++ic) {
        out[ic].found = data[ic].found;
        for (int ieg = 0; ieg < 2; ++ieg) {
            for (int is = 0; is < 4; ++is) {
                for (int isrc = 0; isrc < kNSrc; ++isrc) {
                    // TProfile SEM difference
                    if (data[ic].ok[ieg][is][isrc] && data[ic].ok[ieg][ref][isrc]) {
                        out[ic].val[ieg][is][isrc] =
                            data[ic].val[ieg][is][isrc] - data[ic].val[ieg][ref][isrc];
                        out[ic].err[ieg][is][isrc] = std::sqrt(
                            data[ic].err[ieg][is][isrc] * data[ic].err[ieg][is][isrc] +
                            data[ic].err[ieg][ref][isrc] * data[ic].err[ieg][ref][isrc]);
                        out[ic].ok[ieg][is][isrc] = true;
                    }
                    // Kahan difference
                    if (data[ic].kok[ieg][is][isrc] && data[ic].kok[ieg][ref][isrc]) {
                        out[ic].kval[ieg][is][isrc] =
                            data[ic].kval[ieg][is][isrc] - data[ic].kval[ieg][ref][isrc];
                        out[ic].kerr[ieg][is][isrc] = std::sqrt(
                            data[ic].kerr[ieg][is][isrc] * data[ic].kerr[ieg][is][isrc] +
                            data[ic].kerr[ieg][ref][isrc] * data[ic].kerr[ieg][ref][isrc]);
                        out[ic].kok[ieg][is][isrc] = true;
                    }
                }
            }
        }
    }
    return out;
}


// ============================================================================
// CORE CANVAS BUILDER
// ============================================================================

/// Build one integrated-ring comparison canvas.
///
/// @param data       One RunVals per x-bin (nullptr/!found entries are skipped).
/// @param labels     x-axis bin labels (one per entry in data).
/// @param title      Canvas super-title (TLatex rendered at top).
/// @param ieg        0 = WithEtaGate, 1 = WithoutEtaGate.
/// @param isrc       Ring source index 0..4 (see kSrcs).
/// @param scenMask   Bitmask of which scenarios to draw (bit k = scenario k).
///                   Default 0xF = all four.  Use 0x8 for BothCuts only.
///                   Use 0xE to suppress NoCuts (handy for net-effect plots).
/// @param famBdries  x-positions of vertical separator lines (brute-force only).
/// @param famLabels  Annotations above each separator (same length as famBdries).
static TCanvas* MakeComparisonCanvas(
    const std::vector<const RunVals*>& data,
    const std::vector<std::string>&    labels,
    const char* title,
    int  ieg,
    int  isrc,
    int  scenMask   = 0xF,
    const std::vector<double>&      famBdries = {},
    const std::vector<std::string>& famLabels = {}
)
{
    int N = (int)data.size();
    if (N == 0) return nullptr;

    // Canvas width: narrow for small families, wide for brute-force
    int pixPerBin = (N > 40) ? 55 : (N > 20) ? 80 : (N > 8) ? 120 : 180;
    int cW = std::max(700, pixPerBin * N + 280);
    int cH = 520;

    TCanvas* c = new TCanvas(Form("cCmp_%05d", gCI++), "", cW, cH);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.03);
    gPad->SetBottomMargin((N > 8) ? 0.26 : 0.18);
    gPad->SetTopMargin(0.13);

    // Global Y range (all active scenarios, both estimators)
    double yMax = 0.;
    for (int ic = 0; ic < N; ++ic) {
        if (!data[ic] || !data[ic]->found) continue;
        for (int is = 0; is < 4; ++is) {
            if (!(scenMask & (1 << is))) continue;
            if (data[ic]->ok[ieg][is][isrc]) {
                double v = std::fabs(data[ic]->val[ieg][is][isrc])
                         + data[ic]->err[ieg][is][isrc];
                if (v > yMax) yMax = v;
            }
        }
    }
    yMax = (yMax < 1.e-9) ? 1.e-4 : yMax * 1.55;

    // Axis frame (N bins, one per x-entry)
    TH1D* hA = new TH1D(Form("hFr%05d", gCI),
                         ";<Configuration>;<R_{proxy}^{jet}>", N, 0., double(N));
    hA->SetDirectory(nullptr);
    hA->SetStats(0);
    hA->GetYaxis()->SetRangeUser(-yMax, yMax);
    hA->GetYaxis()->SetTitleSize(0.042);
    hA->GetYaxis()->SetLabelSize(0.036);
    double lblSz = (N > 40) ? 0.020 : (N > 20) ? 0.026 : (N > 8) ? 0.034 : 0.040;
    hA->GetXaxis()->SetLabelSize(lblSz);
    for (int ic = 0; ic < N; ++ic)
        hA->GetXaxis()->SetBinLabel(ic + 1, labels[ic].c_str());
    if (N > 5) hA->GetXaxis()->LabelsOption("v");
    hA->Draw("AXIS");

    // Zero reference line
    TLine* zl = new TLine(0., 0., double(N), 0.);
    zl->SetLineColor(kGray + 2); zl->SetLineStyle(2); zl->SetLineWidth(1);
    zl->Draw("SAME");

    // Family boundary separators (brute-force)
    for (int ib = 0; ib < (int)famBdries.size(); ++ib) {
        double xb = famBdries[ib];
        TLine* vl = new TLine(xb, -yMax * 0.97, xb, yMax * 0.97);
        vl->SetLineColor(kGray); vl->SetLineStyle(3); vl->SetLineWidth(1);
        vl->Draw("SAME");
        if (ib < (int)famLabels.size() && !famLabels[ib].empty()) {
            TLatex* ann = new TLatex(xb + 1.3, yMax * 0.45, famLabels[ib].c_str());
            ann->SetTextSize((N > 40) ? 0.018 : 0.024);
            ann->SetTextAlign(12);
            ann->SetTextAngle(90); // Rotate 90 degrees counterclockwise to make the new version of famLabels fit with the shortName string attached
            ann->Draw();
        }
    }

    // X-offsets for 4 scenarios (centred pattern, mirrored from Fig14)
    static const double dX[4] = {-0.27, -0.09, +0.09, +0.27};

    // Draw one TGraphErrors per scenario (filled marker = TProfile SEM;
    // diamond overlay = Kahan SEM for cross-check)
    for (int is = 0; is < 4; ++is) {
        if (!(scenMask & (1 << is))) continue;

        std::vector<double> xs, ys, ex, ey;
        std::vector<double> xK, yK, eKx, eKy;

        for (int ic = 0; ic < N; ++ic) {
            if (!data[ic] || !data[ic]->found) continue;
            double xc = ic + 0.5 + dX[is];
            if (data[ic]->ok[ieg][is][isrc]) {
                xs.push_back(xc);
                ys.push_back(data[ic]->val[ieg][is][isrc]);
                ex.push_back(0.);
                ey.push_back(data[ic]->err[ieg][is][isrc]);
            }
            if (data[ic]->kok[ieg][is][isrc]) {
                xK.push_back(xc);
                yK.push_back(data[ic]->kval[ieg][is][isrc]);
                eKx.push_back(0.);
                eKy.push_back(data[ic]->kerr[ieg][is][isrc]);
            }
        }

        if (!xs.empty()) {
            TGraphErrors* gr = new TGraphErrors((int)xs.size(),
                xs.data(), ys.data(), ex.data(), ey.data());
            gr->SetMarkerStyle(kScenMkF[is]); gr->SetMarkerSize(kMkSz);
            gr->SetMarkerColor(kScenCols[is]); gr->SetLineColor(kScenCols[is]);
            gr->SetLineWidth(kLineW);
            gr->Draw("P SAME");
        }
        // Kahan diamond: same color, slightly larger, no error bar drawn
        // (kerr should be virtually identical to err; shown as a cross-check)
        if (!xK.empty()) {
            TGraphErrors* grK = new TGraphErrors((int)xK.size(),
                xK.data(), yK.data(), eKx.data(), eKy.data());
            grK->SetMarkerStyle(33); grK->SetMarkerSize(kMkSz * 1.30);
            grK->SetMarkerColor(kScenCols[is]); grK->SetLineColor(kScenCols[is]);
            grK->SetLineWidth(1);
            grK->Draw("P SAME");
        }
    }

    // Canvas title
    AddLabel(0.5, 0.955, title, 0.028, 22);

    // Legend -- placed in bottom-left
    double ly1 = 0.33;
    double ly2 = 0.48;
    TLegend* leg = MakeLeg(0.13, ly1, 0.43, ly2);
    for (int is = 0; is < 4; ++is) {
        if (!(scenMask & (1 << is))) continue;
        TLine* dl = new TLine(); dl->SetLineColor(kScenCols[is]); dl->SetLineWidth(kLineW);
        leg->AddEntry(dl, kScenLabels[is], "l");
    }
    TMarker* mF = new TMarker(0.,0.,20); mF->SetMarkerColor(kBlack); mF->SetMarkerSize(kMkSz);
    TMarker* mD = new TMarker(0.,0.,33); mD->SetMarkerColor(kBlack); mD->SetMarkerSize(kMkSz*1.30);
    leg->AddEntry(mF, "TProfile SEM", "p");
    leg->AddEntry(mD, "Kahan SEM", "p");
    leg->Draw();

    return c;
}


// ============================================================================
// WITHIN-FAMILY PLOTS
// ============================================================================

static void MakeWithinFamilyPlots(
    const FamilyInfo&         fam,
    const std::vector<RunVals>& data,   ///< one per config, same order as fam.cfgs
    TDirectory*               famOut
)
{
    int N = (int)fam.cfgs.size();
    std::vector<std::string>    labels(N);
    std::vector<const RunVals*> ptrs  (N);
    for (int i = 0; i < N; ++i) { labels[i] = fam.cfgs[i].label; ptrs[i] = &data[i]; }

    // Net-effect dataset: (scenario) - NoCuts for each config
    std::vector<RunVals> netData = ComputeNetEffect(data, /*ref=*/0);
    std::vector<const RunVals*> netPtrs(N);
    for (int i = 0; i < N; ++i) netPtrs[i] = &netData[i];

    for (int ieg = 0; ieg < 2; ++ieg) {
        TDirectory* egDir = Mkdir(famOut, kEtaGates[ieg]);

        // ---- Main comparison (all 4 scenarios) ----
        for (int isrc = 0; isrc < kNSrc; ++isrc) {
            TDirectory* srcDir = Mkdir(egDir, kSrcs[isrc].subDir);
            std::string tit = Form("Within-family F%d %s | %s | %s",
                fam.num, fam.shortName.c_str(),
                kSrcs[isrc].axisLbl, kEtaGLbls[ieg]);
            TCanvas* c = MakeComparisonCanvas(ptrs, labels, tit.c_str(), ieg, isrc);
            if (!c) continue;
            c->SetName(Form("cWF_F%d_%s_%s", fam.num, kEtaGates[ieg], kSrcs[isrc].sufx));
            SaveTo(c, srcDir);
        }

        // ---- Net-effect: (scenario - NoCuts), NoCuts trivially 0 at reference ----
        // scenMask=0xE skips NoCuts in the drawing (it would plot at 0 ± δ, distracting)
        TDirectory* netDir = Mkdir(egDir, "NetEffect");
        for (int isrc = 0; isrc < kNSrc; ++isrc) {
            TDirectory* srcDir = Mkdir(netDir, kSrcs[isrc].subDir);
            std::string tit = Form("Net effect (scen #minus NoCuts) F%d %s | %s | %s",
                fam.num, fam.shortName.c_str(),
                kSrcs[isrc].axisLbl, kEtaGLbls[ieg]);
            // Draw pT/DCA/Both but omit NoCuts (which is exactly 0 after subtraction)
            TCanvas* c = MakeComparisonCanvas(netPtrs, labels, tit.c_str(), ieg, isrc,
                                              /*scenMask=*/0xE);
            if (!c) continue;
            c->SetName(Form("cNet_F%d_%s_%s", fam.num, kEtaGates[ieg], kSrcs[isrc].sufx));
            SaveTo(c, srcDir);
        }
    }
}


// ============================================================================
// CROSS-FAMILY REPRESENTATIVE PLOTS
// ============================================================================

static void MakeCrossFamilyPlots(
    const std::vector<FamilyInfo>&             fams,
    const std::vector<std::vector<RunVals>>&   allData,  ///< [famIdx][cfgIdx]
    TDirectory*                                outDir
)
{
    // Collect representative data + labels from each family
    std::vector<const RunVals*> repPtrs;
    std::vector<std::string>    repLabels;
    for (int if_ = 0; if_ < (int)fams.size(); ++if_) {
        const FamilyInfo& f = fams[if_];
        int ri = f.repIdx;
        if (ri < 0 || ri >= (int)allData[if_].size()) continue;
        repPtrs  .push_back(&allData[if_][ri]);
        repLabels.push_back(Form("F%d:%s", f.num, f.cfgs[ri].label.c_str()));
    }
    if (repPtrs.empty()) return;

    // Net-effect for representatives
    std::vector<RunVals> repNetStore(repPtrs.size());
    for (int i = 0; i < (int)repPtrs.size(); ++i) {
        if (!repPtrs[i]) continue;
        // Build a single-element vector, compute net, extract result
        std::vector<RunVals> tmp(1, *repPtrs[i]);
        std::vector<RunVals> ntmp = ComputeNetEffect(tmp, 0);
        repNetStore[i] = ntmp[0];
    }
    std::vector<const RunVals*> repNetPtrs(repNetStore.size());
    for (int i = 0; i < (int)repNetStore.size(); ++i) repNetPtrs[i] = &repNetStore[i];

    TDirectory* crsDir = Mkdir(outDir, "CrossFamily");

    for (int ieg = 0; ieg < 2; ++ieg) {
        TDirectory* egDir = Mkdir(crsDir, kEtaGates[ieg]);

        // All 4 scenarios
        for (int isrc = 0; isrc < kNSrc; ++isrc) {
            TDirectory* srcDir = Mkdir(egDir, kSrcs[isrc].subDir);
            std::string tit = Form("Cross-family representatives | %s | %s",
                kSrcs[isrc].axisLbl, kEtaGLbls[ieg]);
            TCanvas* c = MakeComparisonCanvas(repPtrs, repLabels, tit.c_str(), ieg, isrc);
            if (!c) continue;
            c->SetName(Form("cCrossRep_%s_%s", kEtaGates[ieg], kSrcs[isrc].sufx));
            SaveTo(c, srcDir);
        }

        // BothCuts-only (cleaner for a direct cross-family physics comparison)
        TDirectory* bcDir = Mkdir(egDir, "BothCutsOnly");
        for (int isrc = 0; isrc < kNSrc; ++isrc) {
            TDirectory* srcDir = Mkdir(bcDir, kSrcs[isrc].subDir);
            std::string tit = Form("Cross-family reps (BothCuts only) | %s | %s",
                kSrcs[isrc].axisLbl, kEtaGLbls[ieg]);
            TCanvas* c = MakeComparisonCanvas(repPtrs, repLabels, tit.c_str(), ieg, isrc,
                                              /*scenMask=*/0x8); // BothCuts = bit 3
            if (!c) continue;
            c->SetName(Form("cCrossRepBC_%s_%s", kEtaGates[ieg], kSrcs[isrc].sufx));
            SaveTo(c, srcDir);
        }

        // Net-effect for representatives (BothCuts-NoCuts per family)
        TDirectory* netDir = Mkdir(egDir, "NetEffect");
        for (int isrc = 0; isrc < kNSrc; ++isrc) {
            TDirectory* srcDir = Mkdir(netDir, kSrcs[isrc].subDir);
            std::string tit = Form("Cross-family net effect | %s | %s",
                kSrcs[isrc].axisLbl, kEtaGLbls[ieg]);
            TCanvas* c = MakeComparisonCanvas(repNetPtrs, repLabels, tit.c_str(), ieg, isrc,
                                              /*scenMask=*/0xE);
            if (!c) continue;
            c->SetName(Form("cCrossNet_%s_%s", kEtaGates[ieg], kSrcs[isrc].sufx));
            SaveTo(c, srcDir);
        }
    }
}


// ============================================================================
// BRUTE-FORCE ALL-CONFIGS PLOTS
// ============================================================================

static void MakeBruteForcePlots(
    const std::vector<FamilyInfo>&           fams,
    const std::vector<std::vector<RunVals>>& allData,
    TDirectory*                              outDir
)
{
    // Flatten all configs in family order
    std::vector<const RunVals*> bfPtrs;
    std::vector<std::string>    bfLabels;
    std::vector<double>         bfBdries;   // x-positions of family separators
    std::vector<std::string>    bfFamLbls;  // family annotations

    for (int if_ = 0; if_ < (int)fams.size(); ++if_) {
        const FamilyInfo& f = fams[if_];
        int nCfg = (int)f.cfgs.size();
        if (nCfg == 0) continue;

        // Family separator before each family except the first
        double startX = (double)bfPtrs.size();
        if (if_ > 0) {
            bfBdries.push_back(startX);
            bfFamLbls.push_back(Form("F%d (%s)", f.num, f.shortName.c_str()));
        }

        for (int ic = 0; ic < nCfg; ++ic) {
            bfPtrs.push_back(&allData[if_][ic]);
            // Labels: "F{N}:{idx}" -- compact enough for a wide canvas
            bfLabels.push_back(Form("F%d:%d", f.num, ic));
        }
    }
    if (bfPtrs.empty()) return;

    // Net-effect parallel dataset
    std::vector<RunVals> bfNetStore(bfPtrs.size());
    {
        int idx = 0;
        for (int if_ = 0; if_ < (int)fams.size(); ++if_) {
            const FamilyInfo& f = fams[if_];
            int nCfg = (int)f.cfgs.size();
            if (nCfg == 0) continue;
            std::vector<RunVals> tmp(allData[if_].begin(), allData[if_].end());
            std::vector<RunVals> ntmp = ComputeNetEffect(tmp, 0);
            for (int ic = 0; ic < nCfg; ++ic) bfNetStore[idx++] = ntmp[ic];
        }
    }
    std::vector<const RunVals*> bfNetPtrs(bfNetStore.size());
    for (int i = 0; i < (int)bfNetStore.size(); ++i) bfNetPtrs[i] = &bfNetStore[i];

    TDirectory* bfDir = Mkdir(outDir, "BruteForce");

    for (int ieg = 0; ieg < 2; ++ieg) {
        TDirectory* egDir = Mkdir(bfDir, kEtaGates[ieg]);

        // All 4 scenarios, all configs
        for (int isrc = 0; isrc < kNSrc; ++isrc) {
            TDirectory* srcDir = Mkdir(egDir, kSrcs[isrc].subDir);
            std::string tit = Form("All configurations -- all families | %s | %s",
                kSrcs[isrc].axisLbl, kEtaGLbls[ieg]);
            TCanvas* c = MakeComparisonCanvas(bfPtrs, bfLabels, tit.c_str(), ieg, isrc,
                                              /*scenMask=*/0xF, bfBdries, bfFamLbls);
            if (!c) continue;
            c->SetName(Form("cBF_%s_%s", kEtaGates[ieg], kSrcs[isrc].sufx));
            SaveTo(c, srcDir);
        }

        // BothCuts only (for a condensed single-series overview)
        TDirectory* bcDir = Mkdir(egDir, "BothCutsOnly");
        for (int isrc = 0; isrc < kNSrc; ++isrc) {
            TDirectory* srcDir = Mkdir(bcDir, kSrcs[isrc].subDir);
            std::string tit = Form("All configurations (BothCuts only) | %s | %s",
                kSrcs[isrc].axisLbl, kEtaGLbls[ieg]);
            TCanvas* c = MakeComparisonCanvas(bfPtrs, bfLabels, tit.c_str(), ieg, isrc,
                                              /*scenMask=*/0x8, bfBdries, bfFamLbls);
            if (!c) continue;
            c->SetName(Form("cBFbc_%s_%s", kEtaGates[ieg], kSrcs[isrc].sufx));
            SaveTo(c, srcDir);
        }

        // Net-effect brute-force
        TDirectory* netDir = Mkdir(egDir, "NetEffect");
        for (int isrc = 0; isrc < kNSrc; ++isrc) {
            TDirectory* srcDir = Mkdir(netDir, kSrcs[isrc].subDir);
            std::string tit = Form("All configurations -- net effect | %s | %s",
                kSrcs[isrc].axisLbl, kEtaGLbls[ieg]);
            TCanvas* c = MakeComparisonCanvas(bfNetPtrs, bfLabels, tit.c_str(), ieg, isrc,
                                              /*scenMask=*/0xE, bfBdries, bfFamLbls);
            if (!c) continue;
            c->SetName(Form("cBFnet_%s_%s", kEtaGates[ieg], kSrcs[isrc].sufx));
            SaveTo(c, srcDir);
        }
    }
}


// ============================================================================
// MAIN ENTRY POINT
// ============================================================================

void compareHelicityFamilies(
    const char* baseDir = "/home/users/cicerodm/RingPol/HelicityToyModelOptRefactTest",
    const char* outFile = ""  ///< empty = auto: BASE_DIR/compare_families.root
)
{
    gROOT->ProcessLine("gErrorIgnoreLevel = kWarning;");
    SetupStyle();

    // Resolve output path
    std::string outPath;
    if (outFile && outFile[0] != '\0') {
        outPath = outFile;
    } else {
        outPath = std::string(baseDir) + "/compare_families.root";
    }

    printf("\n=================================================================\n");
    printf("  compareHelicityFamilies\n");
    printf("  Input  dir : %s\n", baseDir);
    printf("  Output file: %s\n", outPath.c_str());
    printf("=================================================================\n\n");

    // Build family registry
    std::vector<FamilyInfo> fams = BuildFamilies();
    int nFams = (int)fams.size();

    // ---- Data loading phase -------------------------------------------------
    // Pre-load ALL values into memory before opening the output file.
    // This avoids the root file ownership hazard (TFile::Open changes gDirectory).
    printf("  Loading generator ROOT files...\n");

    std::vector<std::vector<RunVals>> allData(nFams);
    int nFound = 0, nMissing = 0;
    for (int if_ = 0; if_ < nFams; ++if_) {
        const FamilyInfo& f = fams[if_];
        int nCfg = (int)f.cfgs.size();
        allData[if_].resize(nCfg);
        for (int ic = 0; ic < nCfg; ++ic) {
            allData[if_][ic] = LoadRunVals(std::string(baseDir), f.cfgs[ic]);
            if (allData[if_][ic].found) ++nFound; else ++nMissing;
        }
    }
    printf("  Loaded: %d found, %d missing/skipped\n\n", nFound, nMissing);

    if (nFound == 0) {
        printf("  ERROR: No ROOT files found under '%s'.\n", baseDir);
        printf("         Run the generator first, or check BASE_DIR.\n");
        return;
    }

    // ---- Output file --------------------------------------------------------
    TFile* fout = new TFile(outPath.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        printf("  ERROR: Cannot create '%s'. Aborting.\n", outPath.c_str());
        return;
    }

    // ---- Per-family within-family plots ------------------------------------
    printf("  Producing within-family plots...\n");
    for (int if_ = 0; if_ < nFams; ++if_) {
        const FamilyInfo& f = fams[if_];
        // Count how many configs in this family actually have data
        int nPresent = 0;
        for (const auto& rv : allData[if_]) if (rv.found) ++nPresent;
        if (nPresent == 0) {
            printf("    F%d %-14s -- no data found, skipping\n",
                   f.num, f.shortName.c_str());
            continue;
        }
        printf("    F%d %-14s -- %d/%d configs\n",
               f.num, f.shortName.c_str(), nPresent, (int)f.cfgs.size());

        std::string famDirName = Form("Family_%02d_%s", f.num, f.shortName.c_str());
        TDirectory* famOut = Mkdir(fout, famDirName.c_str());
        MakeWithinFamilyPlots(f, allData[if_], famOut);
    }

    // ---- Cross-family representative plots ---------------------------------
    printf("\n  Producing cross-family comparison plots...\n");
    MakeCrossFamilyPlots(fams, allData, fout);

    // ---- Brute-force all-configs plots -------------------------------------
    printf("  Producing brute-force all-configs plots...\n");
    MakeBruteForcePlots(fams, allData, fout);

    // ---- Finalize ----------------------------------------------------------
    fout->Write("", TObject::kOverwrite);
    fout->Close();

    printf("\n=================================================================\n");
    printf("  Output written: %s\n", outPath.c_str());
    printf("  Canvas count  : ~%d total (within-family + cross + brute-force)\n",
           nFound * 10 + 60);  // rough estimate
    printf("  To browse: root -l '%s'\n", outPath.c_str());
    printf("=================================================================\n\n");
}


// ============================================================================
// STANDALONE ENTRY POINT (for compiled execution outside ROOT CLING)
// ============================================================================

#ifndef __CLING__
int main(int argc, char** argv)
{
    const char* baseDir = (argc > 1) ? argv[1]
        : "/home/users/cicerodm/RingPol/HelicityToyModelOptRefactTest";
    const char* outFile = (argc > 2) ? argv[2] : "";
    compareHelicityFamilies(baseDir, outFile);
    return 0;
}
#endif
// end of compareHelicityFamilies.cxx
