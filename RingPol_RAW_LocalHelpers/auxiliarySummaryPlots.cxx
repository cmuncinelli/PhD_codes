/*
* ============================================================
* auxiliarySummaryPlots.cxx
* Extracts and compares various 1D TProfiles (e.g., Ring Observable 
* vs Eta of different reference particles) across multiple Consumer 
* configuration outputs within a single wagon's working directory.
*   Generates combined canvases for different analysis families 
* (Lambda, AntiLambda, BothHyperons) and their respective systematic 
* variations (Data-like Jet, Rand Jet, etc.), appropriately scaled.
*
* Observables are declared in a single table (see kObservableGroups in the main
* function): each entry carries its own in-file directory, so profiles living in
* different folders of the consumer output (EtaDependence/, ProxyPtDependence/,
* the cut-folder root, or task-level folders such as EtaStudy/) all flow through
* the same machinery.
*
* Usage: root -l -b -q 'auxiliarySummaryPlots.cxx("path/to/wagon/dir", "[optional_MC_ref_path]")'
* (or just use the bash coordinator script, which properly compiles this into production-level code)
* ============================================================
*/

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <utility>   // std::swap, used when canonicalising the tanh fit parameter signs

#include "TFile.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TBox.h"
#include "TLegend.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TGraphAsymmErrors.h"
#include <functional> // std::function, for the subset-selection predicates
#include "TStyle.h"
#include "TMath.h"
#include "TLine.h"
#include "TString.h"
#include "TSystem.h"
#include "TDirectory.h"

// ---------------------------------------------------------
// Configuration Data Structures
// ---------------------------------------------------------

/// @brief Styling and file suffix for a single systematic variation.
/// @note isData marks the reference curve: it is drawn last, thicker, and acts as the minuend
///       of every subtraction. Exactly one entry per canvas should set it.
struct VariationConfig {
    std::string suffix;      // e.g., "JustLambda_forceDatalikeJet"
    std::string legendLabel; // e.g., "Data-like Jet"
    int color;
    int lineStyle;
    int markerStyle;
    bool isData;             // Flag to indicate if this is the main, thick black line

    // --- Which canvases this variation appears on -------------------------------------------
    // The full systematics canvas carries everything, which had become unreadable: nine curves plus
    // data, several of them gated/ungated pairs that differ from each other and from nothing else.
    // These flags split that into one complete canvas plus focused ones, without duplicating the
    // table. All have defaults, so an ordinary entry need not spell any of them out.

    bool inRedux = true; // On the reduced main canvas. Set false for variations whose interest is a
                         // comparison against a specific partner rather than against the data.

    // "Does a minimum-pT gate on artificial proxies change anything?" Both members of each
    // gated/ungated pair carry this, so the pairs can be read without unrelated curves on top.
    bool inPtGateStudy = false;

    // "Does removing the eta acceptance gate on artificial proxies change anything?" Note the
    // polarity is the OPPOSITE of the pT study: eta gating is now the DEFAULT, and the variation is
    // the one that switches it off. So the pair is (ungated variation, gated parent), not the other
    // way round, and the difference is read as the cost of NOT gating.
    bool inEtaGateStudy = false;
};

/// @brief A hyperon selection (Lambda, AntiLambda, Both), i.e. one set of consumer output files.
/// @note familyName names the output folder; dataSuffix is what actually appears in the filenames.
struct FamilyConfig {
    std::string familyName;                 // e.g., "Lambda"
    std::string dataSuffix;                 // e.g., "JustLambda"
};

/// @brief One 1D TProfile to extract, together with where it lives and how to label it.
struct ProfileConfig {
    std::string subDir;       // In-file directory, WITHOUT the task dir (e.g. "EtaDependence/")
    std::string profileName;  // e.g., "pRingObservableEtaLeadP"
    std::string xAxisTitle;   // e.g., "#eta_{LeadP}"
    std::string yAxisTitle;   // e.g., "R"
    bool underCutFolder;      // true  -> <task>/<cutFolder>/<subDir><profileName> (per kinematic-cut family)
                              // false -> <task>/<subDir><profileName>            (task-level, e.g. EtaStudy/)

    // One-vs-one canvases against every single variation:
    // (By far the largest contributor to the ~20 extra canvases per observable)
    // Defaults true only for the LeadP proxy.
    bool doIndividualComparisons = false;

    // Name of the matching Toy Model profile inside TOY_SCENARIO_DIR. Empty (the default) means the toy
    // has no counterpart for this observable, and the toy is then skipped for it rather than guessed at.
    std::string toyProfile = "";

    // Fit R(eta) = p0 * tanh(p1 * eta) to this observable, once per variation, plus the cross-variation
    // parameter summaries. Only meaningful for the fine-binned eta profiles: the coarse ones do not have
    // enough points across the turn-over for the two parameters to separate.
    bool doTanhFit = false;

    // Which proxy defines the ring for this observable ("LeadP", "LeadJet", "SubJet"). Used to group the
    // brute-force axis, so observables sharing an entry set sit together instead of being interleaved.
    // Empty means "keep off the brute-force axis entirely": that is what the eta-sign split observables
    // get, since they cover different subsets and would stretch the scale for everything else.
    std::string proxyTag = "";

    // Hyperon species this observable is defined for: "Lambda", "AntiLambda", "LambdaLike", or ""
    // for any. Some consumer objects are booked per species, and the family already selects the
    // species through its input file -- so drawing the AntiLambda-specific profile inside the
    // Lambda family gives an empty or meaningless curve. Empty (the default) keeps an observable on
    // every family, which is what almost all of them want.
    // README: "LambdaLike" is the mixed-species set, so it is only meaningful for BothHyperons;
    // in a single-species family it would just duplicate that species' own profile.
    std::string speciesTag = "";

    // Folder and label name for this observable. Empty means "derive it" -- see
    // ObservableDisplayName. Set it only where the derived name is still unwieldy.
    std::string displayName = "";
};

/**
 * @brief Folder name for one observable: the explicit displayName if set, else a derived one.
 *
 * The consumer's object names encode their ROOT type and a long "what is this" prefix, which made
 * output folders like Lambda/EtaProxy/pRingObservableEtaLeadP -- almost all of it repeated from the
 * group folder one level up. Stripping mechanically handles the bulk; the cases where the result is
 * still unreadable carry an explicit displayName in the table, where the judgment is visible.
 *
 * The rule, in order:
 *   1. Strip the ROOT type marker, "p2d" or "p".
 *   2. Strip a leading "RingObservable", but ONLY if something remains and it does not start with
 *      "Vs". That guard is what keeps pRingVsCentrality readable: "RingVsCentrality" is a good name
 *      and stripping further would leave the meaningless "VsCentrality".
 *
 * README: prefer fixing a bad derived name with a displayName override rather than by extending the
 * rule. The rule is meant to stay small enough to hold in your head.
 */
std::string DeriveObservableName(const std::string& objectName, const std::string& displayName) {
    if (!displayName.empty()) return displayName;

    std::string s = objectName;
    if (s.rfind("p2d", 0) == 0)      s = s.substr(3);
    else if (s.rfind("p", 0) == 0)   s = s.substr(1);

    const std::string kPrefix = "RingObservable";
    if (s.rfind(kPrefix, 0) == 0) {
        const std::string rest = s.substr(kPrefix.size());
        if (!rest.empty() && rest.rfind("Vs", 0) != 0) s = rest;
    }
    return s.empty() ? objectName : s;
}

/// @brief Folder name for a ProfileConfig.
std::string ObservableDisplayName(const ProfileConfig& profConfig) {
    return DeriveObservableName(profConfig.profileName, profConfig.displayName);
}

/// @brief Is this observable meaningful for the family currently being drawn?
bool SpeciesMatchesFamily(const std::string& speciesTag, const std::string& familyName) {
    if (speciesTag.empty()) return true;
    if (familyName == "BothHyperons") return true; // All three sets say something distinct here
    return speciesTag == familyName;               // "Lambda" / "AntiLambda"
}

/// @brief A set of profiles that share an output folder and the same set of valid transforms.
/// @note Adding a new observable to the macro means adding an entry here and nothing else.
/// @note Whether the Toy Model is overlaid is decided per profile, via ProfileConfig::toyProfile,
///       since within one group some observables have a toy counterpart and others do not.
struct ObservableGroup {
    std::string groupName;                  // Output folder name, e.g. "EtaProxy"
    std::vector<ProfileConfig> profiles;
    bool doFold;                            // Fold about x = 0. Only meaningful for axes that are
                                            // symmetric around zero (eta, PVz). Must stay false for
                                            // pT, mass, DeltaPhi, ... where folding is nonsense.
};

/**
 * @brief An external reference dataset drawn alongside the in-wagon systematics
 *        (MC baseline, pp baseline, Toy Model, ...).
 *
 * Collapsing these into one list removes the previously triplicated
 * "if (pMC) ... if (pPP) ... if (pToy) ..." blocks, and with them the manual index-offset
 * arithmetic that had to be kept in sync by hand every time a reference was added.
 *
 * @note An empty basePath means "not requested", and the reference is skipped everywhere with
 *       no further bookkeeping.
 */
struct ExternalRef {
    enum Kind {
        kConsumerDir, // basePath is a directory holding ConsumerResults_<suffix>.root
        kToyModel     // basePath is the full path to a single Toy Model file
    };

    VariationConfig config;
    std::string     basePath; // Empty means "not requested", and the reference is skipped entirely
    Kind            kind;
};

// ---------------------------------------------------------
// Global Definitions
// ---------------------------------------------------------
// The consumer writes everything below this single top-level task directory
const std::string TASK_DIR_IN_FILE = "lambdajetpolarizationionsderived/";

// Default kinematic-cut folder. The consumer books up to four of these
// (Ring, RingKinematicCuts, JetKinematicCuts, JetAndLambdaKinematicCuts), but only
// "Ring" is guaranteed to be enabled, so the others are opt-in via the cutFolder argument.
const std::string DEFAULT_CUT_FOLDER = "Ring";

// Toy Model scenario to compare against. The toy books a 2 x 4 x 3 grid of
// {WithEtaGate, WithoutEtaGate} x {NoCuts, pTCutOnly, DCACutOnly, BothCuts} x {EtaPos, EtaNeg, All}.
// WithEtaGate/BothCuts is the internally consistent one (the toy's own header notes that
// WithoutEtaGate exists only to demonstrate why omitting the daughter acceptance gate makes the
// physics cuts inconsistent), and "All" is the unsplit sample that matches how the consumer fills.
const std::string TOY_SCENARIO_DIR = "WithEtaGate/BothCuts/All/";

// Toy Model integrated <R>, split by the sign of its (random) jet direction. These three single-bin
// profiles are the toy counterparts of the first three bins of the EtaStudy cut profiles
// ("All Lambda", "eta_proxy >= 0", "eta_proxy < 0"), and are the only part of those 9- and 15-bin
// plots the toy can speak to: it has no hyperon-eta or jet-radius splits.
const std::string TOY_RING_INTEGRATED = "pRingProxyJet";
const std::string TOY_RING_ETA_POS    = "pRingProxyJet_JetEtaPos";
const std::string TOY_RING_ETA_NEG    = "pRingProxyJet_JetEtaNeg";

// ---------------------------------------------------------
// Toy Model coverage audit (HelicityToyModel, current version)
// ---------------------------------------------------------
// Ring-observable TProfiles the toy books per scenario, and what this macro does with them:
//
//   pRingProxyJetVsEtaJet   <R> vs eta_jet          -> used by EtaProxy
//   pRingProxyJetVsEta      <R> vs eta_Lambda       -> used by EtaV0
//   pRingProxyJetVsPt       <R> vs pT_Lambda        -> UNUSED: this macro has no Lambda-pT group.
//                                                      Adding one costs nothing on the toy side, since
//                                                      the profile already exists and the consumer books
//                                                      pRingObservable[LeadP|2ndJet]LambdaPt.
//   pRingProxyJet_JetEtaPos <R> integrated, eta_jet >= 0 -> UNUSED here, but these two are the natural
//   pRingProxyJet_JetEtaNeg <R> integrated, eta_jet <  0    toy counterparts of the EtaStudy eta-split
//                                                           closure test, so they become relevant later.
//   pRingProxy, pRingProxyVsEta, pRingProxyVsPt           -> the z-hat-as-jet variant of the above.
//                                                            A second reference rather than a counterpart.
//
// Observables with NO toy counterpart, and how expensive one would be:
//   DeltaPhi, DeltaTheta, phi_Lambda -> cheap. All the needed angles are already in scope where
//                                       ringProxyJet is computed; each costs a booking and one Fill,
//                                       with no extra sampling and no extra runtime. See the TODOs
//                                       on the corresponding groups below.
//   ProxyPt (and its eta splits)     -> not cheap. The toy's jet is a unit direction with no magnitude.
//   PVz, Centrality, NV0s            -> not applicable. The toy models neither a primary vertex nor
//                                       event activity.
// ---------------------------------------------------------

// ---------------------------------------------------------
// Helper 0: Small ROOT navigation utilities
// ---------------------------------------------------------

/**
 * @brief Turns a human-readable legend label into something usable as a ROOT object name.
 *
 * @note Legend labels carry spaces, dots and parentheses ("OO MC Base", "Perp to Jet (pT gated)"), and the
 * canvas names are built from them. ROOT tolerates such names but they are painful to retrieve with
 * dir->Get(), so every run of non-alphanumeric characters collapses to a single underscore.
 *
 * @return The sanitised label, with no leading or trailing underscore.
 */
std::string SanitizeName(const std::string& label) {
    std::string out;
    bool lastWasSep = false;
    for (char ch : label) {
        const bool keep = (ch >= '0' && ch <= '9') || (ch >= 'a' && ch <= 'z') ||
                          (ch >= 'A' && ch <= 'Z');
        if (keep) {
            out += ch;
            lastWasSep = false;
        } else if (!out.empty() && !lastWasSep) {
            out += '_';
            lastWasSep = true;
        }
    }
    if (!out.empty() && out.back() == '_') out.pop_back();
    return out;
}

/**
 * @brief Returns an existing subdirectory, creating it only if absent.
 *
 * Used to keep directory creation lazy: nothing is created until we actually have a canvas to
 * write into it, so consumer files that do not enable a given cut folder do not leave a trail
 * of empty folders behind.
 *
 * @return The subdirectory, owned by the parent, or nullptr if the parent itself is null.
 */
/**
 * @brief Returns an existing subdirectory, creating it only if absent.
 *
 * Used to keep directory creation lazy: nothing is created until we actually have a canvas to
 * write into it, so consumer files that do not enable a given cut folder do not leave a trail
 * of empty folders behind.
 *
 * @param name May be a slash-separated path ("ProxyPt/EtaSplits"), in which case every level is
 *             created or reused in turn. That is what lets an ObservableGroup declare its position
 *             in the output tree purely through its groupName.
 * @return The deepest subdirectory, owned by its parent, or nullptr if the parent itself is null.
 */
TDirectory* EnsureDir(TDirectory* parent, const std::string& name) {
    if (!parent) return nullptr;

    TDirectory* current = parent;
    size_t start = 0;
    while (start <= name.size()) {
        size_t slash = name.find('/', start);
        std::string part = (slash == std::string::npos) ? name.substr(start)
                                                        : name.substr(start, slash - start);
        if (!part.empty()) {
            TDirectory* next = current->GetDirectory(part.c_str());
            if (!next) next = current->mkdir(part.c_str());
            if (!next) return nullptr;
            current = next;
        }
        if (slash == std::string::npos) break;
        start = slash + 1;
    }
    return current;
}

/**
 * @brief Builds the full in-file path of a profile, honouring the selected cut folder.
 * @param cutFolder Kinematic-cut folder to read from (Ring, JetKinematicCuts, ...). Ignored
 *                  for profiles whose ProfileConfig marks them as task-level.
 * @return Path relative to the file root, ready to hand to TFile::Get.
 */
std::string BuildProfilePath(const ProfileConfig& profConfig, const std::string& cutFolder) {
    if (profConfig.underCutFolder) {
        return TASK_DIR_IN_FILE + cutFolder + "/" + profConfig.subDir + profConfig.profileName;
    }
    return TASK_DIR_IN_FILE + profConfig.subDir + profConfig.profileName;
}

// ---------------------------------------------------------
// Helper 1: Cached file access
// ---------------------------------------------------------
/**
 * @brief Opens each input file at most once and keeps it open for the whole run.
 *
 * Every fetch used to open and close the ROOT file again, so the number of opens scaled as
 * (families x observables x variations). With the observable list growing this became a real
 * cost for the macro, and it is pure overhead: the same handful of files is reopened
 * hundreds of times.  Each file is now opened once, kept open, and closed at the very end.
 * Failed opens are cached too, so a missing systematic warns once instead of once per observable.
 *
 * @note Objects handed out by FetchClone are detached clones, so they stay valid after CloseAll().
 */
class FileCache {
public:
    ~FileCache() { CloseAll(); }

    /**
     * @brief Opens a file, or returns the already-open handle.
     * @return The open file, or nullptr if it could not be opened. A failed path is remembered
     *         and will keep returning nullptr without retrying or re-warning.
     */
    TFile* Open(const std::string& path) {
        auto it = fFiles.find(path);
        if (it != fFiles.end()) return it->second; // Includes cached nullptr for known-bad paths

        TFile* f = TFile::Open(path.c_str(), "READ");
        if (f && f->IsZombie()) {
            f->Close();
            delete f;
            f = nullptr;
        }
        if (!f) {
            std::cerr << "  [Warning] Could not open file: " << path << std::endl;
        }
        fFiles[path] = f;
        return f;
    }

    /**
     * @brief Fetches an object and returns a detached clone that survives CloseAll().
     * @param objPath        Path of the object inside the file, as built by BuildProfilePath.
     * @param warnIfMissing  Set false for genuinely optional objects, to keep the log readable.
     * @return A clone owned by the caller (delete it), or nullptr if the file or object is missing.
     */
    template <typename T>
    T* FetchClone(const std::string& filePath, const std::string& objPath, bool warnIfMissing = true) {
        TFile* f = Open(filePath);
        if (!f) return nullptr;

        T* obj = dynamic_cast<T*>(f->Get(objPath.c_str()));
        if (!obj) {
            if (warnIfMissing) {
                std::cerr << "  [Warning] Object not found: " << objPath << " in " << filePath << std::endl;
            }
            return nullptr;
        }

        T* clone = static_cast<T*>(obj->Clone());
        clone->SetDirectory(nullptr); // Detach immediately: Clone() otherwise attaches to gDirectory,
                                      // which at this point is the output file
        return clone;
    }

    /// @brief Closes and releases every cached file. Safe to call more than once.
    void CloseAll() {
        for (auto& entry : fFiles) {
            if (entry.second) {
                entry.second->Close();
                delete entry.second;
            }
        }
        fFiles.clear();
    }

private:
    std::map<std::string, TFile*> fFiles;
};

// ---------------------------------------------------------
// Helper 1.1: Fold Profile (Sum Positive and Negative Eta)
// ---------------------------------------------------------
/**
 * @brief Folds a profile about x = 0, summing each bin with its mirror image.
 * @param newName Name for the returned histogram. Must be unique among the objects alive at the
 *                same time, since several folded copies of the same source coexist.
 * @return A detached TH1D on the positive half-axis, owned by the caller.
 * @note Only meaningful for axes symmetric about zero (eta, PVz). Callers must gate this on
 *       ObservableGroup::doFold, otherwise folding a pT or mass axis produces nonsense.
 */
TH1D* FoldProfile(TProfile* pIn, const std::string& newName) {
    if (!pIn) return nullptr;
    
    TAxis* ax = pIn->GetXaxis();
    int nBins = ax->GetNbins();
    
    // Dynamically build the positive side bin edges
    std::vector<double> posEdges;
    posEdges.push_back(0.0);
    for(int i = 1; i <= nBins; ++i) {
        if (ax->GetBinLowEdge(i) >= -1e-7) { // Identify bins that start at or above 0
            posEdges.push_back(ax->GetBinUpEdge(i));
        }
    }
    
    TH1D* hFolded = new TH1D(newName.c_str(), pIn->GetTitle(), posEdges.size() - 1, posEdges.data());
    hFolded->SetDirectory(nullptr); // Detach up front so the histogram never lands in the output file
    
    // Map and sum the corresponding bins
    for (int i = 1; i <= hFolded->GetNbinsX(); ++i) {
        double center = hFolded->GetBinCenter(i);
        int binPos = pIn->FindBin(center);
        int binNeg = pIn->FindBin(-center);

        double val = 0, err2 = 0;
        
        if (binPos >= 1 && binPos <= nBins) {
            val += pIn->GetBinContent(binPos);
            err2 += std::pow(pIn->GetBinError(binPos), 2);
        }
        
        // If the central bin spans exactly zero, binPos == binNeg, so we only add it once.
        // Otherwise, add the negative counterpart.
        if (binNeg >= 1 && binNeg <= nBins && binNeg != binPos) {
            val += pIn->GetBinContent(binNeg);
            err2 += std::pow(pIn->GetBinError(binNeg), 2);
        }
        
        hFolded->SetBinContent(i, val);
        hFolded->SetBinError(i, std::sqrt(err2));
    }
    
    return hFolded;
}

// ---------------------------------------------------------
// Helper 1.2: Subtract Profiles (Data - Systematics)
// ---------------------------------------------------------
/**
 * @brief Bin-by-bin difference pData - pSys, with the errors added in quadrature.
 * @param pData Reference (the minuend). Its binning defines the output binning.
 * @param pSys  Variation to subtract. May have a different binning: see the note below.
 * @param newName Name for the returned histogram, unique among simultaneously live objects.
 * @return A detached TH1D holding the difference, owned by the caller.
 * @note Bins where either input has no entries are set to exactly zero rather than left as
 *       spurious differences, so empty regions do not masquerade as a signal.
 * @note If the binnings differ (the Toy Model case) the variation is evaluated by interpolation
 *       and its error taken from the nearest bin, which is approximate by construction.
 */
// (Also handles binning mismatches appropriately)
//
// TODO (or at least a warning):
// the error bars below assume Data and the variation are statistically independent, so the
// combined error is sqrt(eData^2 + eSys^2). That is NOT true here: every variation is computed on the
// same reconstructed Lambda candidates as the data, differing only in the proxy direction, so the two
// are strongly positively correlated. Since Var(A-B) = sA^2 + sB^2 - 2*rho*sA*sB, ignoring the
// covariance OVERESTIMATES the uncertainty, i.e. this is a conservative bound rather than the correct
// error. Recovering rho is not possible after the fact from two independently-filled histograms: it
// needs joint per-event accumulators (cross terms sum(r_A*r_B), sum(r_A*n_B), ...) filled in the
// consumer, or a delete-1 jackknife / bootstrap over events applied directly to the difference.
// See the "Ring uncertainty estimators (covariance included)" note for the full derivation.
TH1D* SubtractProfiles(TH1* pData, TH1* pSys, const std::string& newName) {
    if (!pData || !pSys) return nullptr;
    
    // Safely clone the data binning structure
    TH1D* hSub = nullptr;
    if (auto pDataProf = dynamic_cast<TProfile*>(pData)) {
        hSub = pDataProf->ProjectionX(newName.c_str());
    } else {
        hSub = (TH1D*)pData->Clone(newName.c_str());
    }
    hSub->SetDirectory(nullptr); // Detach before anything else, ProjectionX/Clone attach to gDirectory
    hSub->Reset(); // Clear original contents and errors

    // Carry any categorical bin labels across. Reset() leaves them in place on most ROOT versions, but
    // ProjectionX does not always reproduce them, and the EtaStudy difference/pull panels are unreadable
    // without them. Copying explicitly costs nothing and makes the behaviour version-independent.
    for (int i = 1; i <= hSub->GetNbinsX(); ++i) {
        const char* label = pData->GetXaxis()->GetBinLabel(i);
        if (label && label[0] != '\0') hSub->GetXaxis()->SetBinLabel(i, label);
    }
    
    // Verify if the binning structures are identical
    bool binningMatches = (pData->GetNbinsX() == pSys->GetNbinsX() &&
                           std::abs(pData->GetXaxis()->GetXmin() - pSys->GetXaxis()->GetXmin()) < 1e-6 &&
                           std::abs(pData->GetXaxis()->GetXmax() - pSys->GetXaxis()->GetXmax()) < 1e-6);
    
    for (int i = 1; i <= hSub->GetNbinsX(); ++i) {
        double xCenter = pData->GetXaxis()->GetBinCenter(i);
        double vData = pData->GetBinContent(i);
        double eData = pData->GetBinError(i);
        
        double vSys = 0.0;
        double eSys = 0.0;
        bool sysHasInfo = false;

        if (binningMatches) { // Safe bin-by-bin alignment case
            vSys = pSys->GetBinContent(i);
            eSys = pSys->GetBinError(i);
            sysHasInfo = (eSys > 1e-9 || std::abs(vSys) > 1e-9 || (dynamic_cast<TProfile*>(pSys) && ((TProfile*)pSys)->GetBinEntries(i) > 0));
        } else {
            // Binning mismatch (e.g., Toy Model)! 
            // Protect boundaries: Only evaluate if the data bin center falls within the systematic axis range
            if (xCenter >= pSys->GetXaxis()->GetXmin() && xCenter <= pSys->GetXaxis()->GetXmax()) {
                vSys = pSys->Interpolate(xCenter); // A very convenient function in TProfiles! This allows for 
                
                // For error propagation on mismatched grids, proxy with the closest bin's error
                // (could interpolate or do something fancier, but come on! Interpolation is already somewhat wrong,
                // in the sense that we could have binned everything consistently, but this works as a good approximation)
                int closestBin = pSys->FindBin(xCenter);
                eSys = pSys->GetBinError(closestBin);
                sysHasInfo = true; 
            }
        }
        
        bool dataHasInfo = (eData > 1e-9 || std::abs(vData) > 1e-9 || (dynamic_cast<TProfile*>(pData) && ((TProfile*)pData)->GetBinEntries(i) > 0));
        
        if (dataHasInfo && sysHasInfo) {
            hSub->SetBinContent(i, vData - vSys);
            hSub->SetBinError(i, std::sqrt(eData*eData + eSys*eSys));
        } else {
            hSub->SetBinContent(i, 0.0);
            hSub->SetBinError(i, 0.0);
        }
    }
    
    return hSub;
}

// ---------------------------------------------------------
// Helper 1.3: Integrate Profile (Merge all bins)
// ---------------------------------------------------------
/**
 * @brief Collapses a profile to its single entry-weighted mean over the whole axis.
 * @return {value, error}. Both are zero for a null input.
 */
// TODO: this should go away. The consumer already writes the integrated value, once, in bin 1 of
//       the TProfile1Ds IntegratedCuts/pRingCuts, pRingCutsSubLeadingJet and pRingCutsLeadingP.
//       Re-integrating a multi-bin profile here recomputes a number that already exists, by a
//       different route, which is exactly how two values start disagreeing after an innocent
//       rebin. Read the consumer's value instead and keep this only as a cross-check.
//       (Companion note: the two-dimensional IntegratedCuts/p2dRingCuts*V0MassPeak carry the same
//       X axis crossed with a mass flag, where Y bin 1 is strictly OUT of the mass peak and Y bin 2
//       is strictly IN it.)
//
// KNOWN LIMITATION, and the one case where the two routes legitimately differ: this integrates the
// bins of the histogram it is given, so anything in the UNDERFLOW or OVERFLOW is silently excluded,
// while the consumer's integrated profile counts every candidate. Angular observables are safe --
// their axes cover the full interval by construction -- but the proxy-pT and EtaV0 axes do overflow,
// confirmed, with a small shift in some fake-ring estimators (RandJet without the eta gate among
// them). So a disagreement on those two is expected rather than a bug.
// README: this limitation is the reason the consumer's value is authoritative, not merely tidier.
std::pair<double, double> GetIntegratedProfile(TProfile* pIn) {
    if (!pIn) return {0.0, 0.0};
    
    // TProfile::Rebin properly merges bins, weighting them by their respective 
    // entries (fBinEntries) and correctly propagating the variance.
    // This seemed quicker (and simpler!) than any other solution...
    TProfile* pClone = (TProfile*)pIn->Clone("tmp_integrate");
    pClone->SetDirectory(nullptr); // Keep the scratch clone out of the output file
    pClone->Rebin(pClone->GetNbinsX()); 
    
    double val = pClone->GetBinContent(1);
    double err = pClone->GetBinError(1);
    
    delete pClone;
    return {val, err};
}

// ---------------------------------------------------------
// Helper 1.4: Single categorical point
// ---------------------------------------------------------
/**
 * @brief Builds a histogram holding a single point, for the labelled-axis summary canvases.
 *
 * The integrated summaries draw one point per variation on a shared labelled axis, which means
 * one single-filled histogram per variation, all sharing the same binning so they can be
 * overlaid. This used to be written out inline once per reference type; it is identical every
 * time apart from the bin index, so it lives here now.
 *
 * @param name  Object name. Only needs to be unique among simultaneously live objects, since the
 *              histogram is detached from any directory.
 * @param idx   Zero-based position of this variation on the categorical axis (bin idx + 1 is filled).
 * @param nCats Total number of categories on the axis, i.e. the bin count shared by every point
 *              drawn on the same canvas.
 * @return A detached TH1D owned by the caller.
 */
TH1D* MakeCategoricalPoint(const std::string& name, int idx, int nCats, double value, double error) {
    TH1D* h = new TH1D(name.c_str(), "", nCats, 0, nCats);
    h->SetDirectory(nullptr); // Also silences the "Replacing existing TH1" warnings these used to raise,
                              // since the same names recur for every observable and family
    h->SetBinContent(idx + 1, value);
    h->SetBinError(idx + 1, error);
    return h;
}

// ---------------------------------------------------------
// Helper 1.5: Toy Model Profile Extraction
// ---------------------------------------------------------
/**
 * @brief Fetches the Toy Model counterpart of a consumer profile.
 * @param filePath    Full path to the single Toy Model file. An empty string skips the fetch.
 * @param toyProfile  Profile name inside TOY_SCENARIO_DIR, taken from ProfileConfig::toyProfile.
 *                    An empty string means the toy has no counterpart, and nothing is fetched.
 * @return A detached clone owned by the caller, or nullptr if unavailable.
 * @note The Toy Model jet direction is always random and the hyperons are generated as Lambdas, so a
 *       single toy curve is compared against the Lambda, AntiLambda and BothHyperons families alike.
 *       It is a reference for the geometric/efficiency effect, not a per-species prediction.
 * @note Toy binnings do not match the consumer's, which SubtractProfiles handles by interpolation.
 */
TProfile* GetToyModelProfile(FileCache& cache, const std::string& filePath, const std::string& toyProfile) {
    if (filePath.empty() || toyProfile.empty()) return nullptr;
    return cache.FetchClone<TProfile>(filePath, TOY_SCENARIO_DIR + toyProfile);
}

// ---------------------------------------------------------
// Helper 2: Draw Comparison Canvas
// ---------------------------------------------------------

/// @brief Pairs a fetched curve with the styling it should be drawn in.
/// @note Held as TH1* rather than TProfile* because the folded and subtracted transforms produce
///       plain TH1Ds, and all three flavours have to flow through the same drawing helpers.
struct ProfileBundle {
    TH1* profile; // Converted to TH1* which is also the parent type of TProfile* in order to add the R(eta_pos) + R(eta_pos) plots into this structure
    VariationConfig config;
};

/**
 * @brief What to show in the lower pad of a comparison canvas, underneath the overlay.
 *
 * @note There is deliberately no ratio mode. The ring observable is signed and crosses zero (that is
 *       the whole point of the eta-dependence plots), so Data/Variation diverges wherever the
 *       denominator passes through zero, right in the middle of the plot!
 */
enum LowerPadMode {
    kNoLowerPad,
    kDifference, // Data - Variation, errors added in quadrature
    kPull,       // (Data - Variation) / sqrt(sigmaData^2 + sigmaVar^2)
    kSelfPull    // Each curve divided by its OWN error, for canvases that already show a difference
};

/// @brief Canvas-name suffix identifying a lower-pad variant. Empty for kNoLowerPad.
std::string LowerPadSuffix(LowerPadMode mode) {
    if (mode == kDifference) return "_Diff";
    if (mode == kPull) return "_Pull";
    if (mode == kSelfPull) return "_Pull";
    return "";
}

/// @brief Y-axis title for the lower pad.
std::string LowerPadYTitle(LowerPadMode mode) {
    if (mode == kDifference) return "D - V";
    if (mode == kPull) return "(D-V)/#sigma";
    if (mode == kSelfPull) return "value/#sigma";
    return "";
}

/**
 * @brief Builds the lower-pad curve for one variation against the data.
 *
 * Both modes are derived from the same subtraction, which is what keeps the binning-mismatch
 * handling in one place: SubtractProfiles already stores the difference as the bin content and
 * sqrt(sigmaData^2 + sigmaVar^2) as the bin error, so the pull is just content divided by error.
 *
 * @param pData Reference curve (the minuend), whose binning defines the output.
 * @param pVar  Variation to compare against.
 * @param name  Object name, unique among simultaneously live objects.
 * @return A detached TH1D owned by the caller, or nullptr if either input is missing.
 * @note Pull points carry no error bar: a pull is already expressed in units of its own sigma, so
 *       drawing one would be double counting. Bins with no combined error are set to zero.
 * @note The pull inherits the independence assumption flagged on SubtractProfiles. Because the
 *       denominator is overestimated, these pulls are biased TOWARDS zero, i.e. a variation can look
 *       more compatible with the data than it really is. Erring that way is the safe direction for a
 *       null test, but it is the reason a "passing" closure here is necessary and NOT SUFFICIENT (!).
 */
/**
 * @brief Divides a curve by its own error bars, bin by bin.
 *
 * Used on canvases that already show a difference, where "how far is each point from zero, in units of
 * its own uncertainty" is the question being asked. On a Data - Variation curve the bin error is already
 * sqrt(sigmaData^2 + sigmaVar^2), so this is numerically the same quantity as kPull -- it just does not
 * need a separate reference curve to compute it from.
 *
 * @return A detached TH1D owned by the caller, or nullptr if the input is missing.
 * @note Carries the same conservative-error caveat as SubtractProfiles: the denominator is overestimated,
 *       so these values are biased towards zero.
 */
TH1D* MakeSelfPullHist(TH1* src, const std::string& name) {
    if (!src) return nullptr;

    TH1D* h = nullptr;
    if (auto prof = dynamic_cast<TProfile*>(src)) h = prof->ProjectionX(name.c_str());
    else h = (TH1D*)src->Clone(name.c_str());
    h->SetDirectory(nullptr);

    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        const double v = src->GetBinContent(i);
        const double e = src->GetBinError(i);
        h->SetBinContent(i, (e > 1e-12) ? v / e : 0.0);
        h->SetBinError(i, 0.0); // A pull is already in units of its own sigma
        const char* lab = src->GetXaxis()->GetBinLabel(i);
        if (lab && lab[0] != '\0') h->GetXaxis()->SetBinLabel(i, lab);
    }
    return h;
}

TH1D* MakeLowerPadHist(TH1* pData, TH1* pVar, LowerPadMode mode, const std::string& name) {
    TH1D* h = SubtractProfiles(pData, pVar, name);
    if (!h) return nullptr;

    if (mode == kPull) {
        for (int i = 1; i <= h->GetNbinsX(); ++i) {
            double diff = h->GetBinContent(i);
            double comb = h->GetBinError(i);
            h->SetBinContent(i, (comb > 1e-12) ? diff / comb : 0.0);
            h->SetBinError(i, 0.0);
        }
    }
    return h;
}

/**
 * @brief Builds an empty framing histogram that reproduces a source's categorical bin labels.
 *
 * TPad::DrawFrame produces a plain numeric axis, which silently discards the bin labels the consumer
 * set on the EtaStudy cut profiles. Those labels ARE the axis there, so those canvases get this
 * instead: same binning, same labels, no content.
 *
 * @param src Source whose binning and labels are copied. Unlabelled bins simply stay numeric.
 * @return A detached, stats-free TH1D owned by the caller.
 */
TH1D* MakeLabelledFrame(TH1* src, const std::string& name, double yMin, double yMax) {
    const int nBins = src->GetNbinsX();
    TH1D* frame = new TH1D(name.c_str(), "", nBins,
                           src->GetXaxis()->GetXmin(), src->GetXaxis()->GetXmax());
    frame->SetDirectory(nullptr);
    frame->SetStats(0);
    for (int i = 1; i <= nBins; ++i) {
        const char* label = src->GetXaxis()->GetBinLabel(i);
        if (label && label[0] != '\0') frame->GetXaxis()->SetBinLabel(i, label);
    }

    // Long selection names ("#eta_{Jet} > R, #eta_{#Lambda} < 0") overlap badly when laid flat across
    // 9 or 15 bins, so anything past a handful of bins gets angled. ROOT's angle options are "d" and "u"
    // for 45 degrees ("u" tilts the way we want here, "d" tilts the other way) and "v" for vertical.
    // Short axes (the three-bin eta-sign digest) stay horizontal.
    if (nBins > 4) frame->GetXaxis()->LabelsOption("u");
    frame->SetMinimum(yMin);
    frame->SetMaximum(yMax);
    return frame;
}

/**
 * @brief A labelled-bin cut-scan profile, whose x-axis is a list of named selections.
 *
 * These live at task level (EtaStudy/) rather than inside a kinematic-cut folder, and are booked only
 * when the consumer's doFakePolDiagnosticsQA switch is on, so they are always fetched defensively.
 */
struct CategoricalObservable {
    std::string subDir;       // In-file directory below the task dir, e.g. "EtaStudy/"
    std::string profileName;
    std::string yAxisTitle;

    // Build the three-bin "All / eta_proxy >= 0 / eta_proxy < 0" digest of bins 1-3, which is the only
    // slice of these plots the Toy Model can be overlaid on. See kToyEtaSign* above.
    bool doEtaSignSummary = false;

    // Companion TProfile2D of shape (same bins) x (out of peak, in peak), projected into two curves.
    // Empty means there is none. Only the LeadingP scan has one.
    std::string massSplitProfile = "";

    // Same convention as ProfileConfig::displayName. Empty means derive it.
    // Declared here with a default so the existing table entries need no change.
    std::string displayName = "";
};

/// @brief Folder name for a CategoricalObservable, by the same rule as ProfileConfig.
std::string ObservableDisplayName(const CategoricalObservable& catObs) {
    return DeriveObservableName(catObs.profileName, catObs.displayName);
}

/**
 * @brief Copies the first three bins of a cut-scan into a standalone three-bin histogram.
 *
 * Bins 1-3 of every EtaStudy scan are "All Lambda", "eta_proxy >= 0" and "eta_proxy < 0". Pulling them
 * out lets the Toy Model, which only knows those three numbers, be compared against them on an axis
 * where every source is fully populated, rather than leaving 6 or 12 empty bins on the full scan.
 *
 * @return A detached TH1D owned by the caller, or nullptr if the source has fewer than three bins.
 */
TH1D* ExtractEtaSignBins(TH1* src, const std::string& name) {
    if (!src || src->GetNbinsX() < 3) return nullptr;

    TH1D* out = new TH1D(name.c_str(), "", 3, 0, 3);
    out->SetDirectory(nullptr);
    out->SetStats(0);
    for (int i = 1; i <= 3; ++i) {
        out->SetBinContent(i, src->GetBinContent(i));
        out->SetBinError(i, src->GetBinError(i));
        const char* label = src->GetXaxis()->GetBinLabel(i);
        if (label && label[0] != '\0') out->GetXaxis()->SetBinLabel(i, label);
    }
    return out;
}

/**
 * @brief Assembles the Toy Model's three eta-sign numbers into a histogram matching ExtractEtaSignBins.
 *
 * @param labelSrc Supplies the bin labels, so the toy curve lines up with the data's own wording.
 * @return A detached TH1D owned by the caller, or nullptr if any of the three toy profiles is missing.
 */
TH1D* BuildToyEtaSignHist(FileCache& cache, const std::string& toyFilePath,
                          TH1* labelSrc, const std::string& name) {
    if (toyFilePath.empty() || !labelSrc) return nullptr;

    const std::string toyNames[3] = {TOY_RING_INTEGRATED, TOY_RING_ETA_POS, TOY_RING_ETA_NEG};

    TH1D* out = new TH1D(name.c_str(), "", 3, 0, 3);
    out->SetDirectory(nullptr);
    out->SetStats(0);

    for (int i = 0; i < 3; ++i) {
        // Each of these is a single-bin TProfile holding one integrated mean and its error
        TProfile* pToy = cache.FetchClone<TProfile>(toyFilePath, TOY_SCENARIO_DIR + toyNames[i]);
        if (!pToy) { delete out; return nullptr; } // Partial toy curves would be misleading: all or nothing
        out->SetBinContent(i + 1, pToy->GetBinContent(1));
        out->SetBinError(i + 1, pToy->GetBinError(1));
        delete pToy;

        const char* label = labelSrc->GetXaxis()->GetBinLabel(i + 1);
        if (label && label[0] != '\0') out->GetXaxis()->SetBinLabel(i + 1, label);
    }
    return out;
}

/**
 * @brief Converts one profile into a zero-x-width point graph, horizontally offset from its siblings.
 *
 * WHY THESE REPLACED SUPERPOSED HISTOGRAMS. Drawing eight TProfiles on one pad with "PE" gave every
 * curve a horizontal error bar half a bin wide (ROOT's default gStyle ErrorX), so eight sets of bars
 * overlapped inside every bin and the markers all sat at the same x. On the eta-split canvases that
 * was unreadable: you could see that curves disagreed, but not which was which or by how much.
 *
 * A point graph fixes both halves of that. The x error is set to zero, because the horizontal bar
 * never carried information -- it was the bin width, which the axis already shows. And each curve is
 * nudged by a fraction of a bin so that markers which would otherwise coincide are separated. The
 * offsets are centred on the bin, so the set stays anchored to the x value it belongs to.
 *
 * This is the same construction DrawOffsetSeriesCanvas already used for the cross-family points; it
 * is factored out here so the differential canvases share it rather than reimplementing it.
 *
 * @param h          Source histogram or profile. Not modified.
 * @param idx        Index of this curve among the overlaid set, used to pick its offset.
 * @param nCurves    Total number of curves sharing the pad.
 * @param spreadFrac Total horizontal spread, as a fraction of the bin width. Zero disables offsets,
 *                   which is what a single-curve canvas wants.
 * @param name       Object name for the returned graph.
 * @return A new TGraphAsymmErrors owned by the caller, or nullptr when the source has no usable bin.
 * @note Bins with neither a value nor an error are skipped rather than plotted at zero. Filling them
 *       would park a marker on the axis wherever a variation has no entry, which reads as a genuine
 *       measurement of zero -- the exact misreading these plots exist to prevent.
 */
/// Total horizontal spread of an overlaid set, as a fraction of the bin width. Small enough that a
/// point stays visibly inside its own bin, large enough to separate eight markers.
constexpr double kOverlaySpreadFrac = 0.45;

TGraphAsymmErrors* MakeOffsetPointGraph(TH1* h, int idx, int nCurves, double spreadFrac,
                                        const std::string& name) {
    if (!h) return nullptr;

    const double step = (nCurves > 1) ? spreadFrac / (nCurves - 1) : 0.0;
    const double frac = (nCurves > 1) ? (idx - 0.5 * (nCurves - 1)) * step : 0.0;

    TGraphAsymmErrors* g = new TGraphAsymmErrors(0);
    g->SetName(name.c_str());

    int np = 0;
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        const double v = h->GetBinContent(i);
        const double e = h->GetBinError(i);
        if (std::abs(v) < 1e-12 && e < 1e-12) continue;
        // The offset is a fraction of THIS bin's width, so variable-width axes stay anchored.
        const double dx = frac * h->GetBinWidth(i);
        g->SetPoint(np, h->GetBinCenter(i) + dx, v);
        g->SetPointError(np, 0.0, 0.0, e, e); // No x width: the offset alone separates the curves
        ++np;
    }
    if (np == 0) { delete g; return nullptr; }
    return g;
}

/// @brief Applies a variation's styling to a point graph, matching the histogram styling it replaces.
void StylePointGraph(TGraphAsymmErrors* g, const VariationConfig& cfg) {
    if (!g) return;
    g->SetLineColor(cfg.color);
    g->SetMarkerColor(cfg.color);
    g->SetLineStyle(cfg.lineStyle);
    g->SetMarkerStyle(cfg.markerStyle);
    g->SetMarkerSize(cfg.isData ? 1.4 : 1.1);
    g->SetLineWidth(cfg.isData ? 3 : 1);
}

/**
 * @brief Overlays several profiles on one canvas and writes it to the output file.
 *
 * Axis ranges are computed from the data actually present (error bars included), so profiles with
 * different binnings can safely share a canvas. Whichever bundle is flagged as data is drawn last,
 * on top of the systematics, and with a thicker line.
 *
 * @param bundles      Curves to draw, each with its own styling. Drawn in the given order, so pass
 *                     the data first for it to head the legend.
 * @param outDir       Destination directory. Nothing is drawn if this is null.
 * @param profConfig   Supplies the default axis titles.
 * @param customXTitle Overrides the x-axis title from profConfig; empty means "use profConfig".
 * @param customYTitle Overrides the y-axis title from profConfig; empty means "use profConfig".
 * @param lowerPad     Adds a lower panel comparing every variation against the data. Silently falls
 *                     back to a single pad when the bundle list has no data entry or nothing to
 *                     compare it against, so callers need not special-case single-curve canvases.
 * @param categorical  Draws a labelled-bin axis copied from the first bundle instead of a numeric one.
 *                     Needed for the EtaStudy cut profiles, whose x-axis is a list of named selections.
 * @param isDifference Marks the canvas as a Data - Variation plot. The data bundle is then dropped from
 *                     both the drawing and the legend (it is identically zero by construction, so its
 *                     marker carries no information), zero is forced inside the y-range, and a dashed
 *                     reference line is drawn there instead.
 */
void DrawComparisonCanvas(const std::vector<ProfileBundle>& inputBundles,
                          const std::string& canvasName,
                          const std::string& canvasTitle,
                          TDirectory* outDir,
                          const ProfileConfig& profConfig,
                          const std::string& customXTitle = "",
                          const std::string& customYTitle = "",
                          LowerPadMode lowerPad = kNoLowerPad,
                          bool categorical = false,
                          bool isDifference = false) {
    
    if (inputBundles.empty() || !outDir) return;

    // On a difference plot the data curve is identically zero by construction, so its markers and its
    // legend entry say nothing that the axis does not already say. Drop it here and draw a dashed zero
    // reference line instead (see below), rather than making every caller filter its own list.
    std::vector<ProfileBundle> bundles;
    if (isDifference) {
        for (const auto& b : inputBundles) {
            if (!b.config.isData) bundles.push_back(b);
        }
    } else {
        bundles = inputBundles;
    }
    if (bundles.empty()) return;

    // Find global Min and Max across all bins and all provided profiles
    double globalMin = 999999.;
    double globalMax = -999999.;
    
    // Dynamically calculate X-axis bounds to accommodate different binnings safely
    // (no longer assuming that all profiles share the same X-axis binning structure)
    double xMin = 999999.;
    double xMax = -999999.;

    for (const auto& bundle : bundles) {
        double currentXMin = bundle.profile->GetXaxis()->GetXmin();
        double currentXMax = bundle.profile->GetXaxis()->GetXmax();
        if (currentXMin < xMin) xMin = currentXMin;
        if (currentXMax > xMax) xMax = currentXMax;

        TH1* p = bundle.profile;
        for (int i = 1; i <= p->GetNbinsX(); ++i) {
            // Check if bin has data (TProfile tracks entries, TH1D we just proxy via error/content)
            // Only consider bins that actually have data - Could probably think of a more ROOT-esque
            // way of doing this using some get min-max, but in this way I guarantee the error bars are
            // also considered in the interval for the Y axis limits
            bool hasData = false;
            if (auto prof = dynamic_cast<TProfile*>(p)) {
                hasData = (prof->GetBinEntries(i) > 0);
            } else {
                hasData = (p->GetBinError(i) > 1e-9 || std::abs(p->GetBinContent(i)) > 1e-9);
            }

            if (hasData) { 
                double val = p->GetBinContent(i);
                double err = p->GetBinError(i);
                if (val - err < globalMin) globalMin = val - err;
                if (val + err > globalMax) globalMax = val + err;
            }
        }
    }

    // Apply 5% margin (with a fallback if the curve is perfectly flat)
    double margin = (globalMax - globalMin) * 0.05;
    if (margin < 1e-12) margin = 0.05; 
    double yMin = globalMin - margin;
    double yMax = globalMax + margin;

    // A difference plot is read against zero, so keep zero in frame even when every variation sits
    // comfortably to one side of it
    if (isDifference) {
        if (yMin > -margin) yMin = -margin;
        if (yMax <  margin) yMax =  margin;
    }

    // 2. Setup Canvas and Legend
    // A lower pad only makes sense when there is a data curve AND something to compare it against,
    // so single-curve canvases (DataOnly, the external-reference-only ones) quietly stay single-pad.
    int dataIdx = -1;
    for (size_t i = 0; i < bundles.size(); ++i) {
        if (bundles[i].config.isData) { dataIdx = static_cast<int>(i); break; }
    }
    // kSelfPull needs no reference curve: every bundle is divided by its own errors. The other modes
    // compare against the data, so they still require a data bundle and something to compare it with.
    const bool useLowerPad = (lowerPad == kSelfPull)
                             ? (!bundles.empty())
                             : ((lowerPad != kNoLowerPad) && (dataIdx >= 0) && (bundles.size() > 1));

    TCanvas* c = new TCanvas(canvasName.c_str(), canvasTitle.c_str(), 800, useLowerPad ? 800 : 600);

    // The lower pad takes 30% of the height. Everything drawn inside it therefore needs its text
    // scaled up by the inverse of that fraction to end up the same apparent size as in the upper pad.
    const double kLowerFraction = 0.30;
    const double textScale = (1.0 - kLowerFraction) / kLowerFraction;

    TPad* padUpper = nullptr;
    TPad* padLower = nullptr;

    if (useLowerPad) {
        padUpper = new TPad((canvasName + "_up").c_str(), "", 0.0, kLowerFraction, 1.0, 1.0);
        padUpper->SetLeftMargin(0.12);
        padUpper->SetBottomMargin(0.02); // Axis labels live on the lower pad instead
        padUpper->SetGridx();
        padUpper->SetGridy();
        padUpper->Draw();

        padLower = new TPad((canvasName + "_lo").c_str(), "", 0.0, 0.0, 1.0, kLowerFraction);
        padLower->SetLeftMargin(0.12);
        padLower->SetTopMargin(0.02);
        padLower->SetBottomMargin(0.32); // Room for the x-axis title, which now sits down here
        padLower->SetGridx();
        padLower->Draw();

        padUpper->cd();
    } else {
        c->SetLeftMargin(0.12); // todo: could probably mess with these margins a little after the macro finishes running
        c->SetBottomMargin(0.12);
        c->SetGridx();
        c->SetGridy();
    }

    // Adjust legend coordinates based on how many items we have, to keep it tidy
    double legBottom = 0.88 - (bundles.size() * 0.04);
    TLegend* leg = new TLegend(0.15, legBottom, 0.45, 0.88); // eyeballed some sizes! Verify this later on!
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);

    // Determine the Axes titles
    std::string xTitle = customXTitle.empty() ? profConfig.xAxisTitle : customXTitle;
    std::string yTitle = customYTitle.empty() ? profConfig.yAxisTitle : customYTitle;

    // 3. Draw an empty frame to guarantee proper axes scaling
    // In categorical mode the frame is built by hand so the bin labels survive; DrawFrame would
    // replace them with a numeric axis. Hand-built frames are owned here and deleted after the write.
    TH1* frame = nullptr;
    TH1D* ownedUpperFrame = nullptr;
    if (categorical) {
        if (useLowerPad) padUpper->cd(); else c->cd();
        ownedUpperFrame = MakeLabelledFrame(bundles[0].profile, canvasName + "_frameUp", yMin, yMax);
        ownedUpperFrame->Draw();
        frame = ownedUpperFrame;
    } else {
        frame = useLowerPad ? padUpper->DrawFrame(xMin, yMin, xMax, yMax)
                            : c->DrawFrame(xMin, yMin, xMax, yMax);
    }
    frame->SetTitle(canvasTitle.c_str());
    frame->GetXaxis()->SetTitle(useLowerPad ? "" : xTitle.c_str()); // With a lower pad the x-axis is labelled down there
    frame->GetYaxis()->SetTitle(yTitle.c_str());
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);
    if (useLowerPad) {
        frame->GetXaxis()->SetLabelSize(0.0); // Tick numbers would collide with the lower pad's own
    }

    // 4. Build one zero-x-width point graph per curve, horizontally offset from its siblings.
    //    See MakeOffsetPointGraph for why these replaced the superposed histograms.
    //    A single-curve canvas gets no offset: there is nothing to separate it from.
    const int nCurves = static_cast<int>(bundles.size());
    const double spread = (nCurves > 1) ? kOverlaySpreadFrac : 0.0;

    std::vector<TGraphAsymmErrors*> upperGraphs(nCurves, nullptr);
    TGraphAsymmErrors* dataGraph = nullptr;
    bool hasData = false;

    for (int i = 0; i < nCurves; ++i) {
        TGraphAsymmErrors* g = MakeOffsetPointGraph(bundles[i].profile, i, nCurves, spread,
                                                    canvasName + Form("_g%d", i));
        if (!g) continue; // A variation with no usable bin is dropped rather than drawn empty
        StylePointGraph(g, bundles[i].config);
        upperGraphs[i] = g;

        if (bundles[i].config.isData) { dataGraph = g; hasData = true; }
        // Data is passed first, so it heads the legend in the order the caller intended
        leg->AddEntry(g, bundles[i].config.legendLabel.c_str(), "pe");
    }

    // 5. Draw the systematics first, then the data on top of them
    for (int i = 0; i < nCurves; ++i) {
        if (upperGraphs[i] && !bundles[i].config.isData) upperGraphs[i]->Draw("P SAME");
    }
    if (hasData && dataGraph) dataGraph->Draw("P SAME");

    // Dashed reference at zero, replacing the flat data curve that difference plots used to carry
    TLine* zeroRef = nullptr;
    if (isDifference) {
        zeroRef = new TLine(xMin, 0.0, xMax, 0.0);
        zeroRef->SetLineColor(kBlack);
        zeroRef->SetLineStyle(2);
        zeroRef->SetLineWidth(2);
        zeroRef->Draw("SAME");
    }

    leg->Draw();

    // 5b. Lower pad: every variation compared against the data
    std::vector<TH1*> lowerToDelete;
    std::vector<TGraphAsymmErrors*> lowerGraphs; // Owned here, deleted after the canvas is written
    TH1D* lowerFrameOwned = nullptr; // Only set in categorical mode, where the frame is built by hand
    if (useLowerPad) {
        padLower->cd();

        // Build the comparison curves first, so the pad's y-range can be set from what is actually there
        std::vector<ProfileBundle> lowerBundles;
        std::vector<size_t> lowerUpperIdx;   // Upper-pad index of each lower curve, for a matching offset
        for (size_t i = 0; i < bundles.size(); ++i) {
            if (lowerPad != kSelfPull && static_cast<int>(i) == dataIdx) continue; // Data against itself is zero

            std::string lowName = canvasName + "_lower_" + std::to_string(i);
            TH1D* hLow = (lowerPad == kSelfPull)
                ? MakeSelfPullHist(bundles[i].profile, lowName)
                : MakeLowerPadHist(bundles[dataIdx].profile, bundles[i].profile, lowerPad, lowName);
            if (!hLow) continue;

            hLow->SetLineColor(bundles[i].config.color);
            hLow->SetMarkerColor(bundles[i].config.color);
            hLow->SetLineStyle(bundles[i].config.lineStyle);
            hLow->SetMarkerStyle(bundles[i].config.markerStyle);
            hLow->SetLineWidth(1);

            lowerBundles.push_back({hLow, bundles[i].config});
            lowerUpperIdx.push_back(i);
            lowerToDelete.push_back(hLow);
        }

        // Y-range: driven by the data for the difference, but pinned to a sensible fixed scale for the
        // pull, where the interesting question is only "how many sigma", and where a single wild bin
        // would otherwise compress everything else into a flat line.
        double lowMin = 0.0, lowMax = 0.0;
        if (lowerPad == kPull || lowerPad == kSelfPull) {
            lowMin = -4.0;
            lowMax = 4.0;
            for (const auto& lb : lowerBundles) { // Widen (never shrink) if points genuinely run off scale
                for (int i = 1; i <= lb.profile->GetNbinsX(); ++i) {
                    double v = lb.profile->GetBinContent(i);
                    if (std::abs(v) < 1e-12) continue;
                    if (v < lowMin) lowMin = std::max(v * 1.1, -12.0); // Clamped: outliers must not flatten the rest
                    if (v > lowMax) lowMax = std::min(v * 1.1, 12.0);
                }
            }
        } else {
            double gMin = 999999., gMax = -999999.;
            for (const auto& lb : lowerBundles) {
                for (int i = 1; i <= lb.profile->GetNbinsX(); ++i) {
                    double v = lb.profile->GetBinContent(i);
                    double e = lb.profile->GetBinError(i);
                    if (std::abs(v) < 1e-12 && e < 1e-12) continue;
                    if (v - e < gMin) gMin = v - e;
                    if (v + e > gMax) gMax = v + e;
                }
            }
            if (gMax < gMin) { gMin = -1.0; gMax = 1.0; } // Nothing to show: keep a valid frame anyway
            double m = (gMax - gMin) * 0.1;
            if (m < 1e-12) m = 0.05;
            lowMin = std::min(gMin - m, -m); // Always keep zero inside the frame
            lowMax = std::max(gMax + m, m);
        }

        TH1* lowFrame = nullptr;
        TH1D* ownedLowerFrame = nullptr;
        if (categorical) {
            ownedLowerFrame = MakeLabelledFrame(bundles[0].profile, canvasName + "_frameLo", lowMin, lowMax);
            ownedLowerFrame->Draw();
            lowFrame = ownedLowerFrame;
        } else {
            lowFrame = padLower->DrawFrame(xMin, lowMin, xMax, lowMax);
        }
        lowerFrameOwned = ownedLowerFrame;
        lowFrame->GetXaxis()->SetTitle(xTitle.c_str());
        lowFrame->GetYaxis()->SetTitle(LowerPadYTitle(lowerPad).c_str());
        // Text is scaled by the inverse pad fraction so it matches the upper pad's apparent size
        lowFrame->GetXaxis()->SetTitleSize(0.045 * textScale);
        lowFrame->GetXaxis()->SetLabelSize(0.04 * textScale);
        lowFrame->GetYaxis()->SetTitleSize(0.045 * textScale);
        lowFrame->GetYaxis()->SetLabelSize(0.04 * textScale);
        lowFrame->GetYaxis()->SetTitleOffset(0.45);
        lowFrame->GetYaxis()->SetNdivisions(505); // Fewer ticks: the pad is short

        // Sigma bands for the pull. These are added to the pad BEFORE the points, so in principle the
        // points paint on top of them -- but a solid fill still swallowed them in practice, so the bands
        // are hatched rather than solid. A hatched fill physically cannot hide a marker: the gaps between
        // the hatch lines are transparent, so every point stays visible whatever the paint order ends up
        // being. Alpha was tried first and abandoned, since ROOT's transparency is ignored by several
        // canvas backends.
        if (lowerPad == kPull || lowerPad == kSelfPull) {
            TBox* band2 = new TBox(xMin, std::max(-2.0, lowMin), xMax, std::min(2.0, lowMax));
            band2->SetFillColor(kYellow - 3);
            band2->SetFillStyle(3005); // Coarse diagonal hatching, +/- 2 sigma
            band2->SetLineColor(0);
            band2->Draw("SAME");

            TBox* band1 = new TBox(xMin, std::max(-1.0, lowMin), xMax, std::min(1.0, lowMax));
            band1->SetFillColor(kGreen - 3);
            band1->SetFillStyle(3004); // Opposite diagonal, +/- 1 sigma, so the two bands read apart
            band1->SetLineColor(0);
            band1->Draw("SAME");
        }

        TLine* zeroLine = new TLine(xMin, 0.0, xMax, 0.0);
        zeroLine->SetLineColor(kBlack);
        zeroLine->SetLineStyle(2);
        zeroLine->SetLineWidth(2);
        zeroLine->Draw("SAME");

        // Same point-graph treatment as the upper pad, and crucially the SAME offset index, so a
        // curve sits at the same x in both panels and the eye can track it between them.
        for (size_t k = 0; k < lowerBundles.size(); ++k) {
            const bool isPull = (lowerPad == kPull || lowerPad == kSelfPull);
            TGraphAsymmErrors* gLow = MakeOffsetPointGraph(
                lowerBundles[k].profile, static_cast<int>(lowerUpperIdx[k]), nCurves, spread,
                canvasName + Form("_glow%zu", k));
            if (!gLow) continue;
            StylePointGraph(gLow, lowerBundles[k].config);
            if (isPull) {
                // A pull is already in units of its own sigma, so an error bar on it is meaningless
                for (int ip = 0; ip < gLow->GetN(); ++ip) gLow->SetPointError(ip, 0.0, 0.0, 0.0, 0.0);
            }
            gLow->Draw("P SAME");
            lowerGraphs.push_back(gLow);
        }
        padLower->RedrawAxis(); // Bands are filled shapes, so put the tick marks and frame back on top
    }

    // 6. Save to the specific directory within our output file
    outDir->cd();
    c->Write();
    
    // Clean up to prevent memory leaks during the loop (object should already have been written)
    for (auto p : lowerToDelete) delete p;
    for (auto g : upperGraphs) delete g;   // nullptr entries are fine: a dropped curve made no graph
    for (auto g : lowerGraphs) delete g;
    delete ownedUpperFrame; // nullptr unless categorical, and deleting nullptr is fine
    delete lowerFrameOwned;
    delete leg;
    delete c; 
}

// ---------------------------------------------------------
// Helper 3.5: tanh fit of the eta dependence
// ---------------------------------------------------------

/// @brief Outcome of one R(eta) = p0*tanh(p1*eta) fit. Errors are magnitudes, already "sign-canonicalised".
struct TanhFitResult {
    bool   valid = false; // Fit converged with a usable covariance and at least one degree of freedom
    double entries = 0.0;
    double chi2 = 0.0;
    int    ndf = 0;
    double p0 = 0.0, p0Low = 0.0, p0High = 0.0;
    double p1 = 0.0, p1Low = 0.0, p1High = 0.0;
};

/**
 * @brief Fits p0*tanh(p1*eta) to one profile, draws the annotated canvas, and returns the parameters.
 *
 * Fit configuration, and why:
 *  - Unconstrained, then canonicalised. The tanh equation is exactly invariant under (p0,p1) -> (-p0,-p1), so
 *    leaving both free gives two equivalent minima and makes convergence and errors ill-defined. Hard
 *    limits would fix that but introduce their own pathologies (Minos behaves badly against a boundary,
 *    and a limit at p1 >= 0 is hit whenever the distribution is flat). Fitting free and then flipping
 *    both signs when p1 < 0 removes the redundancy exactly, costs nothing numerically, and leaves p0
 *    carrying the sign of the slope.
 *    Note the asymmetric errors swap sides under that flip, which is handled below.
 *
 *  - Minos ("E"), which is worth its cost here specifically. For small p1*eta the model degenerates to
 *    p0*p1*eta, so only the PRODUCT is constrained and the p0/p1 contour becomes a long curved valley.
 *    Parabolic (Hesse) errors describe that valley badly, which is exactly the regime the data sits in at
 *    current statistics. Minos maps the real contour instead, and the resulting asymmetry is itself the
 *    diagnostic that the two parameters have not separated.
 *
 *  - "M" (improve) to re-search after convergence, for the same valley reason, and "S" to get the full
 *    result object back rather than just a status code.
 *
 * @param prof       Profile to fit. Empty bins carry zero error and are skipped by the chi2 method.
 * @param canvasName Name of the canvas written into outDir.
 * @return The fit outcome. A non-converged fit still reports chi2/ndf and entries, with valid = false.
 * @note Bins are means with SEM errors, so a chi2 fit with those errors is the correct method here.
 *       A likelihood option would be wrong: these are not Poisson counts.
 */
TanhFitResult FitTanhProfile(TProfile* prof,
                             const ProfileConfig& profConfig,
                             const VariationConfig& varConfig,
                             const std::string& canvasName,
                             TDirectory* outDir) {
    TanhFitResult out;
    if (!prof || !outDir) return out;

    out.entries = prof->GetEntries();

    const double xLo = prof->GetXaxis()->GetXmin();
    const double xHi = prof->GetXaxis()->GetXmax();

    // Named parameters, so the fit box shows "p0" and "p1" rather than "p0" style indices
    TF1* fTanh = new TF1((canvasName + "_fn").c_str(), "[p0]*TMath::TanH([p1]*x)", xLo, xHi);
    fTanh->SetParameters(1.0e-2, 0.5); // Orders of magnitude taken from the Toy Model fits
    fTanh->SetLineColor(varConfig.color);
    fTanh->SetLineWidth(2);
    fTanh->SetNpx(500);

    // gStyle must be set BEFORE anything is painted: the stats box is built at paint time and reads
    // these values then. Setting them after Draw (as this did) had no effect at all.
    // OptStat 10 keeps only "Entries"; OptFit 111 keeps chi2/ndf, parameter values and errors while
    // dropping the fit-probability line.
    gStyle->SetOptStat(10);
    gStyle->SetOptFit(111);

    TFitResultPtr res = prof->Fit(fTanh, "QSEM", "", xLo, xHi);

    TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str(), 800, 600);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);
    c->SetGridx();
    c->SetGridy();

    prof->SetStats(1);
    prof->SetLineColor(varConfig.color);
    prof->SetMarkerColor(varConfig.color);
    prof->SetMarkerStyle(varConfig.markerStyle);
    prof->SetLineWidth(varConfig.isData ? 3 : 1);
    prof->GetXaxis()->SetTitle(profConfig.xAxisTitle.c_str());
    prof->GetYaxis()->SetTitle(profConfig.yAxisTitle.c_str());
    prof->Draw("PE");

    // Force a paint pass so the TPaveStats actually gets built. TCanvas::Write() does NOT paint, so
    // without this the canvas was stored with no stats object at all, and ROOT rebuilt a default one
    // from whatever gStyle happened to be active when the file was later opened. That is precisely why
    // the box came back with default contents and no fit information in it.
    c->Update();

    // Having got hold of the object, set the options ON IT rather than only on gStyle. Those settings
    // are then serialised with the canvas and survive independently of the reader's global style.
    TPaveStats* stats = dynamic_cast<TPaveStats*>(prof->FindObject("stats"));
    if (stats) {
        stats->SetOptStat(10);
        stats->SetOptFit(111);
        // Upper LEFT corner: the tanh rises to the right, so the top-left is the empty part of the pad
        stats->SetX1NDC(0.14);
        stats->SetX2NDC(0.46);
        stats->SetY1NDC(0.72);
        stats->SetY2NDC(0.90);
        stats->SetFillColor(kWhite);
        stats->SetFillStyle(1001); // Opaque, so the grid does not run through the numbers
        stats->SetBorderSize(1);
        stats->SetTextSize(0.030);
    }

    // Function definition, sitting immediately below the stat box
    TPaveText* fnLabel = new TPaveText(0.14, 0.645, 0.46, 0.712, "NDC");
    fnLabel->SetFillColor(kWhite);
    fnLabel->SetFillStyle(1001); // Opaque white, as requested: transparent was unreadable over the grid
    fnLabel->SetBorderSize(1);
    fnLabel->SetTextSize(0.028);
    fnLabel->AddText(("f(" + profConfig.xAxisTitle + ")=[p0]*TMath::TanH([p1]*x)").c_str());
    fnLabel->Draw();

    // Second pass, so the repositioned stats box and the new label are part of what gets written
    c->Modified();
    c->Update();

    outDir->cd();
    c->Write();

    // Collect the numbers before things get out of scope:
    if (res.Get()) {
        out.chi2 = res->Chi2();
        out.ndf  = res->Ndf();
        out.valid = res->IsValid() && (out.ndf > 0);

        out.p0 = res->Parameter(0);
        out.p1 = res->Parameter(1);

        // Minos gives asymmetric errors; fall back to the parabolic one where it did not run.
        // std::abs guards against the Minuit sign convention for the lower side.
        out.p0Low  = res->HasMinosError(0) ? std::abs(res->LowerError(0)) : res->ParError(0);
        out.p0High = res->HasMinosError(0) ? std::abs(res->UpperError(0)) : res->ParError(0);
        out.p1Low  = res->HasMinosError(1) ? std::abs(res->LowerError(1)) : res->ParError(1);
        out.p1High = res->HasMinosError(1) ? std::abs(res->UpperError(1)) : res->ParError(1);

        // Canonicalise onto p1 > 0. Flipping a parameter's sign mirrors its error interval, so the
        // low and high sides swap; the magnitudes themselves are unchanged.
        if (out.p1 < 0.0) {
            out.p0 = -out.p0;
            out.p1 = -out.p1;
            std::swap(out.p0Low, out.p0High);
            std::swap(out.p1Low, out.p1High);
        }
    }

    delete fnLabel;
    delete c;
    delete fTanh;
    return out;
}

/// @brief One variation's value on a cross-variation summary axis.
struct FitSummaryPoint {
    VariationConfig config;
    double value = 0.0;
    double errLow = 0.0;
    double errHigh = 0.0;
};

/**
 * @brief Draws one point per variation on a labelled axis, using asymmetric error bars.
 *
 * The counterpart of DrawIntegratedCanvas for quantities whose errors are not symmetric. One
 * single-point TGraphAsymmErrors per variation, so each keeps its own colour and marker, and no
 * line is ever drawn between points (they are unordered categories, not a series).
 *
 * @param labels Bin labels in axis order. Must be the same length as points.
 * @param points One entry per label. Pass errors of zero for quantities that have none.
 */
// ------------------------------------------------------------------------------------------------
// Mass-region styling, shared by every canvas that slices on the three-bin mass axis.
// ------------------------------------------------------------------------------------------------
// Canvas_MassSelection, Canvas_InOutOfPeak and Canvas_SidebandAgreement all sit in the same folder
// and all show some subset of the same four curves. If each picked its own colours, the reader would
// have to re-learn the legend three times over. Defined once here so they cannot drift apart.
// README: any new mass-region canvas must draw its styling from these.
struct MassRegionStyle { int color; int marker; };
constexpr MassRegionStyle kStyleInPeak    = {kTeal + 3,    34}; // Taken from the massVariations table so the two mass-selection routes agree on sight
constexpr MassRegionStyle kStyleOutOfPeak = {kMagenta + 2, 23};
constexpr MassRegionStyle kStyleLeftSideband  = {kAzure + 2, 21};
constexpr MassRegionStyle kStyleRightSideband = {kGreen + 3, 22};

VariationConfig MassRegionConfig(const char* label, MassRegionStyle st, bool isReference) {
    VariationConfig c;
    c.legendLabel = label;
    c.color = st.color;
    c.lineStyle = 1;
    c.markerStyle = st.marker;
    c.isData = isReference; // The reference drives the pull panel
    return c;
}

/**
 * @brief The reduced subset of a bundle list: data, externals, and the variations flagged inRedux.
 *
 * README: every canvas that draws the full variation set should also emit a redux twin, from the
 * SAME call site. Emitting them apart is how two loops end up disagreeing about which variations
 * they show -- exactly the drift that let an all-V0s canvas into the Lambda folder.
 *
 * External references keep the default inRedux = true: they are context rather than clutter.
 */
std::vector<ProfileBundle> ReduxSubset(const std::vector<ProfileBundle>& bundles) {
    std::vector<ProfileBundle> out;
    for (const auto& b : bundles)
        if (b.config.isData || b.config.inRedux) out.push_back(b);
    return out;
}

/**
 * @brief Draws <R>(x) separately for each mass slice of a TProfile2D, on one canvas.
 *
 * The cheap alternative to a full signal extraction: instead of fitting the peak and subtracting a
 * background, the mass axis is cut into a peak region and its two sidebands, and the observable is
 * simply read off in each. Three curves on one canvas then answer, by eye, whether the observable
 * behaves differently inside the mass window than beside it -- which is what a real signal would
 * look like -- with none of the fitting machinery and none of the QA that machinery needs.
 *
 * README: this is deliberately independent of signalExtractionRing. It is a cross-check on the
 * extraction, so it must not share the extraction's assumptions.
 *
 * Drawn through the shared comparison canvas, so the three curves get the same zero-x-width offset
 * points as everything else: the in-window slice is flagged as the reference, and the pull panel
 * then answers "is each sideband consistent with the peak region?" directly.
 *
 * @param p2d      Source, with the observable axis on x and mass on y. Not modified.
 * @param outDir   Destination. Nothing is drawn if this or p2d is null.
 * @return Nothing; the canvas is written to outDir.
 */
void DrawMassSliceCanvas(TProfile2D* p2d,
                         const std::string& canvasName,
                         const std::string& canvasTitle,
                         TDirectory* outDir,
                         const std::string& xTitle,
                         const std::string& yTitle) {
    if (!p2d || !outDir) return;

    const int nMass = p2d->GetYaxis()->GetNbins();
    if (nMass < 2) return; // A single mass slice is just the inclusive profile

    const MassRegionStyle sliceStyles[3] = {kStyleLeftSideband, kStyleInPeak, kStyleRightSideband};

    std::vector<ProfileBundle> bundles;
    std::vector<TH1*> toDelete;

    for (int m = 1; m <= nMass && m <= 3; ++m) {
        TProfile* slice = p2d->ProfileX(Form("%s_slice%d", canvasName.c_str(), m), m, m);
        if (!slice) continue;
        slice->SetDirectory(nullptr);
        if (slice->GetEntries() < 1) { delete slice; continue; }

        const double mLo = p2d->GetYaxis()->GetBinLowEdge(m);
        const double mHi = p2d->GetYaxis()->GetBinUpEdge(m);
        const char* region = (nMass == 3) ? (m == 1 ? "left sideband"
                                          : m == 2 ? "in mass window"
                                                   : "right sideband")
                                          : "slice";

        // The in-window slice is the reference the sidebands are read against
        VariationConfig cfg = MassRegionConfig(Form("%s [%.4f, %.4f]", region, mLo, mHi),
                                               sliceStyles[(m - 1) % 3], (m == 2));

        bundles.push_back({slice, cfg});
        toDelete.push_back(slice);
    }
    if (bundles.size() < 2) { for (auto p : toDelete) delete p; return; }

    ProfileConfig pc;
    pc.xAxisTitle = xTitle;
    pc.yAxisTitle = yTitle;

    DrawComparisonCanvas(bundles, canvasName, canvasTitle, outDir, pc, xTitle, yTitle);
    // Pull panel against the in-window slice: "is the sideband consistent with the peak region?"
    DrawComparisonCanvas(bundles, canvasName + "_PullOnly", canvasTitle, outDir, pc,
                         xTitle, yTitle, kPull);

    for (auto p : toDelete) delete p;
}

/**
 * @brief The two questions a three-bin mass axis can answer, as two canvases.
 *
 * The single three-curve canvas this replaces mixed them together. They are different questions and
 * a reader should not have to disentangle them:
 *
 *   Canvas_InOutOfPeak       in-window against the two sidebands COMBINED. This is the cheap signal
 *                            extraction: if the observable carries a real signal, the in-window
 *                            curve departs from the background estimate beside it.
 *   Canvas_SidebandAgreement left sideband against right sideband. This is the assumption CHECK: a
 *                            sideband subtraction only makes sense if the two sides agree, and if
 *                            they do not, the number the canvas above produces is not trustworthy
 *                            whatever it shows.
 *
 * README: the second canvas is a precondition for the first. Read it first.
 *
 * The combined sideband is built by cloning the left slice and Add()-ing the right, so it is the
 * exact entry-weighted merge and shares binning with everything else by construction.
 *
 * @param p2d Source with the observable on x and a three-bin mass axis on y. Bin 2 is the peak.
 */
void DrawMassRegionCanvases(TProfile2D* p2d,
                            const std::string& titlePrefix,
                            TDirectory* outDir,
                            const std::string& xTitle,
                            const std::string& yTitle) {
    if (!p2d || !outDir) return;
    if (p2d->GetYaxis()->GetNbins() < 3) return; // Needs both sidebands and a peak

    TProfile* sLeft  = p2d->ProfileX("mrc_left",  1, 1);
    TProfile* sPeak  = p2d->ProfileX("mrc_peak",  2, 2);
    TProfile* sRight = p2d->ProfileX("mrc_right", 3, 3);
    if (!sLeft || !sPeak || !sRight) { delete sLeft; delete sPeak; delete sRight; return; }
    sLeft->SetDirectory(nullptr); sPeak->SetDirectory(nullptr); sRight->SetDirectory(nullptr);

    ProfileConfig pc;
    pc.xAxisTitle = xTitle;
    pc.yAxisTitle = yTitle;

    // ONLY the sideband-agreement canvas is produced here.
    //
    // The in/out-of-peak comparison already exists in this same folder, as Canvas_MassSelection,
    // built from the massVariations table (_excludeOutOfPeak / _excludeInPeak). That version is
    // strictly better: it carries a third curve, the unselected data, which this one cannot. When
    // this function also emitted its own three-region overview and its own in/out canvas, the folder
    // ended up with two different objects named Canvas_MassSelection and three ways of asking one
    // question, in two colour schemes.
    //
    // What the massVariations route CANNOT do is compare the two sidebands against each other: it
    // only knows "in peak" and "out of peak". That is the one thing left for this function, and it
    // is a precondition for trusting either -- a sideband subtraction is only meaningful if the two
    // sides agree.
    // README: the two routes also read DIFFERENT inputs. Canvas_MassSelection comes from separate
    // consumer files produced with a mass cut applied; this one slices the three-bin mass axis of a
    // single nominal file. They should agree, and a disagreement is itself worth knowing about.
    if (sLeft->GetEntries() > 0 && sRight->GetEntries() > 0) {
        // The in-window curve rides along as the anchor, in the SAME style Canvas_MassSelection
        // uses for it, so the two canvases can be read side by side without re-learning the legend.
        std::vector<ProfileBundle> b;
        if (sPeak->GetEntries() > 0)
            b.push_back({sPeak, MassRegionConfig("In Mass Peak", kStyleInPeak, true)});
        b.push_back({sLeft,  MassRegionConfig("left sideband", kStyleLeftSideband, sPeak->GetEntries() <= 0)});
        b.push_back({sRight, MassRegionConfig("right sideband", kStyleRightSideband, false)});
        const std::string t = titlePrefix + ", left vs right sideband";
        DrawComparisonCanvas(b, "Canvas_SidebandAgreement", t, outDir, pc, xTitle, yTitle);
        DrawComparisonCanvas(b, "Canvas_SidebandAgreement_PullOnly", t, outDir, pc, xTitle, yTitle, kPull);
    }

    delete sLeft; delete sPeak; delete sRight;
}

void DrawFitSummaryCanvas(const std::vector<std::string>& labels,
                          const std::vector<FitSummaryPoint>& points,
                          const std::string& canvasName,
                          const std::string& canvasTitle,
                          TDirectory* outDir,
                          const std::string& yTitle) {
    if (points.empty() || !outDir || labels.size() != points.size()) return;

    double gMin = 999999., gMax = -999999.;
    for (const auto& pt : points) {
        if (pt.value - pt.errLow  < gMin) gMin = pt.value - pt.errLow;
        if (pt.value + pt.errHigh > gMax) gMax = pt.value + pt.errHigh;
    }
    double margin = (gMax - gMin) * 0.1;
    if (margin < 1e-12) margin = (std::abs(gMax) > 1e-12) ? std::abs(gMax) * 0.1 : 0.05;

    const int nBins = static_cast<int>(labels.size());
    TCanvas* c = new TCanvas(canvasName.c_str(), canvasTitle.c_str(), std::max(800, nBins * 100), 600);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.15);
    c->SetGridx();

    TH1D* frame = new TH1D((canvasName + "_frame").c_str(), canvasTitle.c_str(), nBins, 0, nBins);
    frame->SetDirectory(nullptr);
    for (int i = 0; i < nBins; ++i) frame->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
    frame->GetYaxis()->SetTitle(yTitle.c_str());
    frame->SetStats(0);
    // Long category names overlap when laid flat; angle them once the axis gets busy
    if (nBins > 4) frame->GetXaxis()->LabelsOption("u");
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->SetMinimum(gMin - margin);
    frame->SetMaximum(gMax + margin);
    frame->Draw();

    std::vector<TGraphAsymmErrors*> graphs;
    for (size_t i = 0; i < points.size(); ++i) {
        TGraphAsymmErrors* g = new TGraphAsymmErrors(1);
        g->SetPoint(0, static_cast<double>(i) + 0.5, points[i].value);
        g->SetPointError(0, 0.0, 0.0, points[i].errLow, points[i].errHigh);
        g->SetLineColor(points[i].config.color);
        g->SetMarkerColor(points[i].config.color);
        g->SetMarkerStyle(points[i].config.markerStyle);
        g->SetMarkerSize(1.5);
        g->SetLineWidth(2);
        g->Draw("P SAME"); // "P" only: never "L", these are categories rather than an ordered series
        graphs.push_back(g);
    }

    outDir->cd();
    c->Write();

    for (auto* g : graphs) delete g;
    delete frame;
    delete c;
}

/**
 * @brief Draws several multi-bin series on one labelled axis, each nudged sideways within its bin.
 *
 * The counterpart of DrawIntegratedCanvas for the case where the SERIES are what is being compared
 * (hyperon selections) rather than the bins. Drawn as TGraphAsymmErrors rather than histograms so the
 * points can be offset: stacked exactly on top of one another their error bars overlap and become
 * impossible to tell apart, which is the whole thing one wants to read here.
 *
 * @param labels     Bin labels in axis order; its size sets the bin count.
 * @param series     One entry per hyperon selection, each holding one value per bin.
 * @param spreadFrac Fraction of a bin width the offsets span, measured between the outermost two series.
 *                   0.5 puts three series at -1/4, 0 and +1/4 of a bin. Points carry no x error bar, so
 *                   this only shifts the marker and its vertical bar.
 * @param addSelfPullPad Adds a lower panel of each point divided by its own error, at the same
 *                   horizontal offsets, so a series lines up vertically between the two panels.
 */
void DrawOffsetSeriesCanvas(const std::vector<std::string>& labels,
                            const std::vector<ProfileBundle>& series,
                            const std::string& canvasName,
                            const std::string& canvasTitle,
                            TDirectory* outDir,
                            const std::string& yTitle,
                            double spreadFrac = 0.5,
                            bool addSelfPullPad = false) {
    if (series.empty() || labels.empty() || !outDir) return;

    const int nBins = static_cast<int>(labels.size());
    const int nSeries = static_cast<int>(series.size());

    double gMin = 999999., gMax = -999999.;
    for (const auto& sr : series) {
        for (int i = 1; i <= sr.profile->GetNbinsX(); ++i) {
            const double v = sr.profile->GetBinContent(i);
            const double e = sr.profile->GetBinError(i);
            if (std::abs(v) < 1e-12 && e < 1e-12) continue;
            if (v - e < gMin) gMin = v - e;
            if (v + e > gMax) gMax = v + e;
        }
    }
    if (gMax < gMin) { gMin = -1.0; gMax = 1.0; }
    double margin = (gMax - gMin) * 0.15;
    if (margin < 1e-12) margin = 0.05;

    TCanvas* c = new TCanvas(canvasName.c_str(), canvasTitle.c_str(),
                             std::max(800, nBins * 110), addSelfPullPad ? 800 : 600);

    const double kLowerFraction = 0.30;
    const double textScale = (1.0 - kLowerFraction) / kLowerFraction;

    TPad* padUpper = nullptr;
    TPad* padLower = nullptr;
    if (addSelfPullPad) {
        padUpper = new TPad((canvasName + "_up").c_str(), "", 0.0, kLowerFraction, 1.0, 1.0);
        padUpper->SetLeftMargin(0.12);
        padUpper->SetBottomMargin(0.02); // Labels live on the lower pad
        padUpper->SetGridy();
        padUpper->Draw();

        padLower = new TPad((canvasName + "_lo").c_str(), "", 0.0, 0.0, 1.0, kLowerFraction);
        padLower->SetLeftMargin(0.12);
        padLower->SetTopMargin(0.02);
        padLower->SetBottomMargin(0.42); // Angled category labels need the room
        padLower->SetGridx();
        padLower->Draw();

        padUpper->cd();
    } else {
        c->SetLeftMargin(0.12);
        c->SetBottomMargin(0.22); // Room for the angled labels
        c->SetGridy();
    }

    TH1D* frame = new TH1D((canvasName + "_frame").c_str(), canvasTitle.c_str(), nBins, 0, nBins);
    frame->SetDirectory(nullptr);
    frame->SetStats(0);
    for (int i = 0; i < nBins; ++i) frame->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
    frame->GetYaxis()->SetTitle(yTitle.c_str());
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(addSelfPullPad ? 0.0 : 0.040); // Hidden when the lower pad carries them
    frame->SetMinimum(gMin - margin);
    frame->SetMaximum(gMax + margin);
    if (!addSelfPullPad && nBins > 4) frame->GetXaxis()->LabelsOption("u");
    frame->Draw();

    TLegend* leg = new TLegend(0.14, 0.80, 0.44, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.032);

    // Offsets are centred on the bin, so the set stays visually anchored to its label
    const double step = (nSeries > 1) ? spreadFrac / (nSeries - 1) : 0.0;

    // Shared by both pads, so a series sits at the same x in the value panel and the pull panel
    auto seriesOffset = [&](int j) {
        return (nSeries > 1) ? (j - 0.5 * (nSeries - 1)) * step : 0.0;
    };

    std::vector<TGraphAsymmErrors*> graphs;
    for (int j = 0; j < nSeries; ++j) {
        const double dx = seriesOffset(j);

        // Points are added only for bins that actually hold a value. Filling every bin regardless would
        // park a marker at exactly zero wherever a series has no entry, which reads as a real measurement
        // sitting on the axis -- harmless in the value panel (usually clipped out of range) but badly
        // misleading in the pull panel, which is pinned around zero.
        TGraphAsymmErrors* g = new TGraphAsymmErrors(0);
        int np = 0;
        for (int i = 0; i < nBins; ++i) {
            const double v = series[j].profile->GetBinContent(i + 1);
            const double e = series[j].profile->GetBinError(i + 1);
            if (std::abs(v) < 1e-12 && e < 1e-12) continue;
            g->SetPoint(np, i + 0.5 + dx, v);
            g->SetPointError(np, 0.0, 0.0, e, e); // No x width: the offset alone separates the series
            ++np;
        }
        g->SetLineColor(series[j].config.color);
        g->SetMarkerColor(series[j].config.color);
        g->SetMarkerStyle(series[j].config.markerStyle);
        g->SetMarkerSize(1.3);
        g->SetLineWidth(2);
        g->Draw("P SAME"); // "P" only: these are categories, never a connected series
        leg->AddEntry(g, series[j].config.legendLabel.c_str(), "PE");
        graphs.push_back(g);
    }
    leg->Draw();

    // Lower panel: each point divided by its own error, keeping the same horizontal offsets
    TH1D* lowFrame = nullptr;
    if (addSelfPullPad) {
        padLower->cd();

        lowFrame = new TH1D((canvasName + "_lowframe").c_str(), "", nBins, 0, nBins);
        lowFrame->SetDirectory(nullptr);
        lowFrame->SetStats(0);
        for (int i = 0; i < nBins; ++i) lowFrame->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
        if (nBins > 4) lowFrame->GetXaxis()->LabelsOption("u");
        lowFrame->GetYaxis()->SetTitle("value/#sigma");
        lowFrame->GetXaxis()->SetLabelSize(0.040 * textScale);
        lowFrame->GetYaxis()->SetTitleSize(0.045 * textScale);
        lowFrame->GetYaxis()->SetLabelSize(0.040 * textScale);
        lowFrame->GetYaxis()->SetTitleOffset(0.45);
        lowFrame->GetYaxis()->SetNdivisions(505);
        lowFrame->SetMinimum(-7.0);
        lowFrame->SetMaximum(7.0);
        lowFrame->Draw();

        TBox* band2 = new TBox(0.0, -2.0, (double)nBins, 2.0);
        band2->SetFillColor(kYellow - 3); band2->SetFillStyle(3005); band2->SetLineColor(0);
        band2->Draw("SAME");
        TBox* band1 = new TBox(0.0, -1.0, (double)nBins, 1.0);
        band1->SetFillColor(kGreen - 3); band1->SetFillStyle(3004); band1->SetLineColor(0);
        band1->Draw("SAME");

        TLine* zeroLine = new TLine(0.0, 0.0, (double)nBins, 0.0);
        zeroLine->SetLineColor(kBlack); zeroLine->SetLineStyle(2); zeroLine->SetLineWidth(2);
        zeroLine->Draw("SAME");

        for (int j = 0; j < nSeries; ++j) {
            const double dx = seriesOffset(j);

            TGraphAsymmErrors* gp = new TGraphAsymmErrors(0);
            int np = 0;
            for (int i = 0; i < nBins; ++i) {
                const double v = series[j].profile->GetBinContent(i + 1);
                const double e = series[j].profile->GetBinError(i + 1);
                if (std::abs(v) < 1e-12 && e < 1e-12) continue;
                gp->SetPoint(np, i + 0.5 + dx, (e > 1e-12) ? v / e : 0.0);
                gp->SetPointError(np, 0.0, 0.0, 0.0, 0.0); // A pull is already in units of its own sigma
                ++np;
            }
            gp->SetLineColor(series[j].config.color);
            gp->SetMarkerColor(series[j].config.color);
            gp->SetMarkerStyle(series[j].config.markerStyle);
            gp->SetMarkerSize(1.3);
            gp->Draw("P SAME");
            graphs.push_back(gp);
        }
        padLower->RedrawAxis(); // Bands are filled shapes: put the frame back on top
    }

    outDir->cd();
    c->Write();

    for (auto* g : graphs) delete g;
    delete leg;
    delete lowFrame;
    delete frame;
    delete c;
}


// ---------------------------------------------------------
// Helper 3: Draw Integrated Ring Observable Canvas
// ---------------------------------------------------------
/**
 * @brief Draws one point per variation on a labelled categorical axis.
 *
 * Each bundle is expected to be a single-point histogram from MakeCategoricalPoint, all sharing
 * the same bin count, so that overlaying them produces one point per category.
 *
 * @param bundles     One single-filled histogram per category, in axis order.
 * @param labels      Bin labels, one per category. Its size defines the axis bin count, so it must
 *                    match the nCats the bundles were built with.
 * @param outDir      Destination directory. Nothing is drawn if this is null.
 * @param isSubtracted Marks the canvas as a difference plot: forces zero to stay inside the y-range
 *                     and draws a dashed line there, so "consistent with zero" is readable at a glance.
 * @param drawLegend  Adds a legend naming each variation. Needed when the x-axis labels describe
 *                    something other than the variations themselves (the brute-force summary).
 * @param addSelfPullPad Adds a lower panel showing each point divided by its own error, so the value
 *                    and its significance can be read off the same canvas instead of two.
 */
void DrawIntegratedCanvas(const std::vector<ProfileBundle>& bundles,
                          const std::vector<std::string>& labels,
                          const std::string& canvasName,
                          const std::string& canvasTitle,
                          TDirectory* outDir,
                          const std::string& yTitle,
                          bool isSubtracted,
                          bool drawLegend = false,
                          bool addSelfPullPad = false) {
    if (bundles.empty() || !outDir) return;

    // 1. Find Min/Max to set Y-axis boundaries
    double globalMin = 999999., globalMax = -999999.;
    for (const auto& bundle : bundles) {
        TH1* p = bundle.profile;
        for (int i = 1; i <= p->GetNbinsX(); ++i) {
            if (p->GetBinError(i) > 1e-9 || std::abs(p->GetBinContent(i)) > 1e-9) {
                double val = p->GetBinContent(i);
                double err = p->GetBinError(i);
                if (val - err < globalMin) globalMin = val - err;
                if (val + err > globalMax) globalMax = val + err;
            }
        }
    }
    
    double margin = (globalMax - globalMin) * 0.1;
    if (margin < 1e-12) margin = 0.05;
    double yMin = globalMin - margin;
    double yMax = globalMax + margin;
    
    // Lock zero line bounds for subtracted plots
    if (isSubtracted && yMin > -margin) yMin = -margin;
    if (isSubtracted && yMax < margin) yMax = margin;

    int nBins = labels.size();
    TCanvas* c = new TCanvas(canvasName.c_str(), canvasTitle.c_str(),
                             std::max(800, nBins * 100), addSelfPullPad ? 800 : 600);

    // The lower panel takes 30% of the height, so its text needs scaling up by the inverse of that
    // fraction to end up the same apparent size as the upper pad's
    const double kLowerFraction = 0.30;
    const double textScale = (1.0 - kLowerFraction) / kLowerFraction;

    TPad* padUpper = nullptr;
    TPad* padLower = nullptr;
    if (addSelfPullPad) {
        padUpper = new TPad((canvasName + "_up").c_str(), "", 0.0, kLowerFraction, 1.0, 1.0);
        padUpper->SetLeftMargin(0.12);
        padUpper->SetBottomMargin(0.02); // Labels live on the lower pad instead
        padUpper->SetGridx();
        padUpper->Draw();

        padLower = new TPad((canvasName + "_lo").c_str(), "", 0.0, 0.0, 1.0, kLowerFraction);
        padLower->SetLeftMargin(0.12);
        padLower->SetTopMargin(0.02);
        padLower->SetBottomMargin(0.42); // Angled category labels need the room
        padLower->SetGridx();
        padLower->Draw();

        padUpper->cd();
    } else {
        c->SetLeftMargin(0.12);
        c->SetBottomMargin(0.15); // Larger bottom margin to fit text labels
        c->SetGridx();
    }

    // 2. Create a dummy histogram just to draw the custom axes and labels
    TH1D* frame = new TH1D((canvasName + "_frame").c_str(), canvasTitle.c_str(), nBins, 0, nBins);
    frame->SetDirectory(nullptr);
    for (int i = 0; i < nBins; ++i) {
        frame->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
    }
    if (addSelfPullPad) frame->GetXaxis()->SetLabelSize(0.0); // Labels would collide with the lower pad's
    frame->GetYaxis()->SetTitle(yTitle.c_str());
    frame->GetYaxis()->SetRangeUser(yMin, yMax);
    // Long category names overlap when laid flat; angle them once the axis gets busy
    if (nBins > 4) frame->GetXaxis()->LabelsOption("u");
    frame->SetStats(0);
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);
    
    // Drawing the frame by actually ensuring the minima as well:
    frame->SetMinimum(yMin);
    frame->SetMaximum(yMax);
    frame->Draw();

    // 3. Draw each integrated point
    for (auto& bundle : bundles) {
        TH1* p = bundle.profile;
        p->SetLineColor(bundle.config.color);
        p->SetMarkerColor(bundle.config.color);
        p->SetLineStyle(bundle.config.lineStyle);
        p->SetMarkerStyle(bundle.config.markerStyle);
        p->SetLineWidth(bundle.config.isData ? 3 : 2);
        p->SetMarkerSize(1.5);
        p->Draw("PE SAME");
    }

    if (isSubtracted) {
        TLine* line = new TLine(0, 0, nBins, 0);
        line->SetLineColor(kBlack); line->SetLineWidth(2); line->SetLineStyle(2);
        line->Draw("SAME");
    }

    // Drawing an extra legend that adds the systematic variation description as legend instead
    // of relying only on the x-axis:
        // This is not a problem for plots such as the Integrated_summary folder's, where the X
        // axis actually explains if we are dealing with randJets, datalike estimators, data
        // estimators, perpToJet, ..., but not for the "brute force" plot where the X axis actually
        // indicates if a given observable was calculated using the FastJet's output, the 2nd highest
        // pT jet from FastJet, or the leading particle! This information has to be in a legend!
    if (drawLegend) {
        // Positioned in the bottom left (it overlapped with the data points, otherwise, as they are
        // mostly positive and the top part is usually filled with data)
        TLegend* leg = new TLegend(0.12, 0.16, 0.38, 0.40);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0); // Transparent background
        leg->SetTextSize(0.035);
        for (const auto& bundle : bundles) {
            leg->AddEntry(bundle.profile, bundle.config.legendLabel.c_str(), "PE");
        }
        leg->Draw();
    }

    // 4. Lower panel: every point divided by its own error bars
    // Each bundle here is a single-filled-bin histogram, so the pull is drawn as a one-point graph at
    // that bin rather than as another histogram. Drawing a histogram would put a spurious marker at zero
    // in every bin the bundle does not own, which on a fixed +/-4 pull scale would be very visible.
    TH1D* lowFrame = nullptr;
    std::vector<TGraphAsymmErrors*> pullGraphs;
    if (addSelfPullPad) {
        padLower->cd();

        lowFrame = new TH1D((canvasName + "_lowframe").c_str(), "", nBins, 0, nBins);
        lowFrame->SetDirectory(nullptr);
        lowFrame->SetStats(0);
        for (int i = 0; i < nBins; ++i) lowFrame->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
        if (nBins > 4) lowFrame->GetXaxis()->LabelsOption("u");
        lowFrame->GetYaxis()->SetTitle("value/#sigma");
        lowFrame->GetXaxis()->SetLabelSize(0.040 * textScale);
        lowFrame->GetYaxis()->SetTitleSize(0.045 * textScale);
        lowFrame->GetYaxis()->SetLabelSize(0.040 * textScale);
        lowFrame->GetYaxis()->SetTitleOffset(0.45);
        lowFrame->GetYaxis()->SetNdivisions(505);
        lowFrame->SetMinimum(-4.0);
        lowFrame->SetMaximum(4.0);
        lowFrame->Draw();

        // Same hatched sigma bands as the differential pull panels, for consistency of reading
        TBox* band2 = new TBox(0.0, -2.0, (double)nBins, 2.0);
        band2->SetFillColor(kYellow - 3); band2->SetFillStyle(3005); band2->SetLineColor(0);
        band2->Draw("SAME");
        TBox* band1 = new TBox(0.0, -1.0, (double)nBins, 1.0);
        band1->SetFillColor(kGreen - 3); band1->SetFillStyle(3004); band1->SetLineColor(0);
        band1->Draw("SAME");

        TLine* zeroLine = new TLine(0.0, 0.0, (double)nBins, 0.0);
        zeroLine->SetLineColor(kBlack); zeroLine->SetLineStyle(2); zeroLine->SetLineWidth(2);
        zeroLine->Draw("SAME");

        for (const auto& bundle : bundles) {
            TH1* p = bundle.profile;
            for (int i = 1; i <= p->GetNbinsX(); ++i) {
                const double v = p->GetBinContent(i);
                const double e = p->GetBinError(i);
                if (e <= 1e-12 && std::abs(v) <= 1e-12) continue; // Bin this bundle does not own

                TGraphAsymmErrors* g = new TGraphAsymmErrors(1);
                g->SetPoint(0, p->GetXaxis()->GetBinCenter(i), (e > 1e-12) ? v / e : 0.0);
                g->SetPointError(0, 0.0, 0.0, 0.0, 0.0); // A pull is already in units of its own sigma
                g->SetMarkerColor(bundle.config.color);
                g->SetLineColor(bundle.config.color);
                g->SetMarkerStyle(bundle.config.markerStyle);
                g->SetMarkerSize(1.5);
                g->Draw("P SAME");
                pullGraphs.push_back(g);
            }
        }
        padLower->RedrawAxis(); // Bands are filled shapes: put the frame back on top
    }

    outDir->cd();
    c->Write();
    for (auto* g : pullGraphs) delete g;
    delete lowFrame;
    delete frame;
    delete c;
}

// ---------------------------------------------------------
// Main Macro 
// ---------------------------------------------------------
/**
 * @brief Builds the cross-configuration summary file for one wagon.
 *
 * Reads every ConsumerResults_<suffix>.root found in consumerDir (one per systematic variation),
 * plus the optional external references, and writes auxiliarySummaryPlots.root beside them.
 *
 * @param consumerDir  Directory holding this wagon's ConsumerResults files. Also the output location.
 * @param mcRefDir     Directory holding the MC baseline's ConsumerResults files. Empty disables the overlay.
 * @param ppRefDir     Directory holding the pp baseline's ConsumerResults files. Empty disables the overlay.
 * @param toyModelPath Full path to the single Toy Model file. Empty disables the overlay.
 * @param cutFolder    Which kinematic-cut folder to read. Only "Ring" is always booked by the consumer;
 *                     the others are opt-in and produce no output at all when absent from the files.
 * @param doIndividualComparisons Master switch for the one-vs-one canvases. Off by default: that block
 *                     alone is close to half of the entire output, and the overlays plus the pull panels
 *                     already carry the same comparisons. Per-observable opt-in still applies on top of
 *                     this, so turning it on enables it only for the profiles flagged in the table.
 */
// =================================================================================================
// SIGNAL-EXTRACTED INTEGRATED OBSERVABLES
// =================================================================================================
// Everything below reads the output of signalExtractionRing.cxx, which is a DIFFERENT file per
// consumer configuration:
//     <sigExtractDir>/signalExtractionRing_<dataSuffix><sysSuffix>.root
// and, inside it,
//     <cutFolder>/IntegratedSummary/<proxy>/<histogram>
//
// WHY THIS IS NOT JUST ANOTHER ObservableGroup. Every other group in this macro integrates a
// TProfile over its x-axis with GetIntegratedProfile. These quantities arrive ALREADY integrated,
// as single-bin histograms, and two of them (<R>_S and <R>_B) cannot be obtained from a consumer
// profile at all: they only exist once a sideband subtraction has been done. So they get their own
// small table and fetch path, and then rejoin the shared drawers -- MakeCategoricalPoint,
// DrawIntegratedCanvas and DrawOffsetSeriesCanvas do all the drawing, with no new plotting code.
//
// WHY THE DIFFERENCES ARE READ FROM FILE RATHER THAN COMPUTED HERE. <R>_measured, <R>_S and <R>_B
// are three functions of the same four primitives (peak counts, peak Sum_R, sideband counts,
// sideband Sum_R) and are strongly correlated. SubtractProfiles adds errors in quadrature, which
// for these is not a conservative approximation but simply wrong, in both directions: on a
// representative set of primitives it overstates the error on <R>_meas - <R>_S by about 60% and
// UNDERSTATES the one on <R>_S - <R>_B by about 20%. The extraction macro computes all three
// differences where the primitives and their covariance are still in scope, and this code plots the
// numbers it is handed. Do not be tempted to re-derive them here.
// =================================================================================================

/// @brief One signal-extracted quantity, i.e. one histogram name under IntegratedSummary/<proxy>/.
struct SigExtractQuantity {
    std::string histName;   // Without the "_<proxy>" suffix the extraction macro appends
    std::string legendLabel;
    std::string yAxisTitle;
    int color;
    int markerStyle;
    bool isData;            // Draws thick, like the Data column elsewhere in this macro
};

/// @brief One proxy folder inside IntegratedSummary/.
/// @note These are the three genuinely independent measurements. The angular extractions
///       (DeltaPhi, DeltaTheta) integrate over the SAME candidates and would only repeat one
///       another, which is why signalExtractionRing no longer publishes them here at all.
struct SigExtractProxy {
    std::string folderName; // "LeadJet", "LeadP", "SubJet"
    std::string legendLabel;
};

/// @brief A value with its uncertainty, plus whether it was actually found.
struct SigExtractPoint {
    double value = 0.0;
    double error = 0.0;
    bool   found = false;
};

/**
 * @brief The integrated <R> as the CONSUMER computed it, from bin 1 of the IntegratedCuts profiles.
 *
 * These TProfile1Ds carry the integrated value for each proxy, filled once over every candidate.
 * Reading them is strictly better than re-integrating a differential profile here: the number is
 * computed in one place, so it cannot drift after a rebin, and it counts candidates that fall in a
 * differential histogram's under- or overflow. See the note on GetIntegratedProfile for the two
 * axes where that difference is confirmed to matter.
 *
 * Bin 1 is "All Lambda". Bins 2-4 are the pT, rapidity and LeadP-eta cut scans, unused for now.
 * README: the two-dimensional companions, IntegratedCuts/p2dRingCuts*V0MassPeak, share this X axis
 * crossed with a mass flag -- Y bin 1 strictly OUT of the mass peak, Y bin 2 strictly IN it.
 *
 * @param proxy "LeadJet", "LeadP" or "SubJet".
 * @return {value, error}, and found = false when the object is absent.
 */
SigExtractPoint FetchIntegratedFromCuts(FileCache& cache, const std::string& filePath,
                                        const std::string& proxy) {
    SigExtractPoint out;
    std::string obj;
    if (proxy == "LeadJet")      obj = "pRingCuts";
    else if (proxy == "LeadP")   obj = "pRingCutsLeadingP";
    else if (proxy == "SubJet")  obj = "pRingCutsSubLeadingJet";
    else return out;

    TProfile* p = cache.FetchClone<TProfile>(filePath, TASK_DIR_IN_FILE + "IntegratedCuts/" + obj, false);
    if (!p) return out;
    if (p->GetNbinsX() >= 1 && p->GetBinEntries(1) > 0) {
        out.value = p->GetBinContent(1);
        out.error = p->GetBinError(1);
        out.found = true;
    }
    delete p;
    return out;
}

/**
 * @brief The integrated <R> split by mass window, from IntegratedCuts/p2dRingCuts*V0MassPeak.
 *
 * Same X axis as the one-dimensional pRingCuts profiles -- bin 1 is "All Lambda" -- crossed with a
 * two-state mass flag on Y, where BIN 1 is strictly OUT of the mass peak and BIN 2 strictly IN it.
 * (Bin numbers, not fill values: the consumer fills 0 and 1 on a {2, 0, 2} axis.)
 *
 * README: this is the cheap counterpart to signalExtractionRing and is deliberately independent of
 * it -- no fit, no sideband model, nothing to go wrong quietly. It answers "how far does the
 * background pull the integrated number?" rather than "what is the signal?", and a large gap
 * between the two canvases is the warning that the extraction has real work to do.
 *
 * @param inPeak true selects Y bin 2, false selects Y bin 1.
 */
SigExtractPoint FetchIntegratedByMassPeak(FileCache& cache, const std::string& filePath,
                                          const std::string& proxy, bool inPeak) {
    SigExtractPoint out;
    std::string obj;
    if (proxy == "LeadJet")     obj = "p2dRingCutsV0MassPeak";
    else if (proxy == "LeadP")  obj = "p2dRingCutsLeadingPV0MassPeak";
    else if (proxy == "SubJet") obj = "p2dRingCutsSubLeadingJetV0MassPeak";
    else return out;

    TProfile2D* p = cache.FetchClone<TProfile2D>(filePath,
                                                 TASK_DIR_IN_FILE + "IntegratedCuts/" + obj, false);
    if (!p) return out;
    const int yBin = inPeak ? 2 : 1;
    const int gBin = p->GetBin(1, yBin); // X bin 1 is "All Lambda"
    if (p->GetBinEntries(gBin) > 0) {
        out.value = p->GetBinContent(gBin);
        out.error = p->GetBinError(gBin);
        out.found = true;
    }
    delete p;
    return out;
}

/**
 * @brief Reads one single-bin histogram out of a signalExtractionRing output file.
 *
 * @param cache         Shared file cache; the extraction files are opened at most once each.
 * @param sigExtractDir Directory holding the signalExtractionRing_*.root files.
 * @param fileSuffix    <dataSuffix><sysSuffix>, exactly as in the consumer file names.
 * @param cutFolder     Top-level variation folder inside the extraction file ("Ring", ...).
 * @param proxy         Folder under IntegratedSummary/.
 * @param histName      Histogram stem; "_<proxy>" is appended, matching the extraction macro.
 * @param requireValid  When true, hExtractionStatus must report a valid extraction. Guards against
 *                      the failure mode where a fit did not converge: a zero-filled single-bin
 *                      histogram is visually indistinguishable from a genuine measurement of zero,
 *                      which is exactly why the extraction macro writes that status flag.
 * @return The point, with found = false when anything is missing or the extraction was invalid.
 */
SigExtractPoint FetchSigExtractPoint(FileCache& cache,
                                     const std::string& sigExtractDir,
                                     const std::string& fileSuffix,
                                     const std::string& cutFolder,
                                     const std::string& proxy,
                                     const std::string& histName,
                                     bool requireValid = true) {
    SigExtractPoint out;
    if (sigExtractDir.empty()) return out;

    const std::string filePath = sigExtractDir + "/signalExtractionRing_" + fileSuffix + ".root";
    const std::string dirPath = cutFolder + "/IntegratedSummary/" + proxy + "/";

    if (requireValid) {
        TH1D* hStatus = cache.FetchClone<TH1D>(filePath, dirPath + "hExtractionStatus_" + proxy, false);
        if (!hStatus) return out;           // No status object at all: treat as not extracted
        const bool valid = hStatus->GetBinContent(1) > 0.5; // Bin 1 is labelled "valid"
        delete hStatus;
        if (!valid) return out;
    }

    TH1D* h = cache.FetchClone<TH1D>(filePath, dirPath + histName + "_" + proxy, false);
    if (!h) return out;
    if (h->GetNbinsX() >= 1) {
        out.value = h->GetBinContent(1);
        out.error = h->GetBinError(1);
        out.found = true;
    }
    delete h;
    return out;
}

/**
 * @brief Builds the "components" canvases: <R>_measured, <R>_S and <R>_B side by side.
 *
 * One categorical axis with three entries, in the style of the Integrated_summary plots, plus a
 * self-pull lower panel so the value and its distance from zero are read off the same canvas.
 * The subtracted companion shows the two differences against <R>_measured, taken from file with
 * their correlations intact (see the header note above).
 */
void DrawSigExtractComponents(FileCache& cache,
                              const std::string& sigExtractDir,
                              const std::string& fileSuffix,
                              const std::string& cutFolder,
                              const SigExtractProxy& proxy,
                              const std::string& familyName,
                              TDirectory* outDir) {
    if (!outDir) return;

    // Order matters: measured first, so it occupies the "Data" slot the eye expects.
    const std::vector<SigExtractQuantity> components = {
        {"hIntegratedRMeas", "Measured", "Integrated <R>", kBlack,    8,  true },
        {"hIntegratedRSig",  "Signal",   "Integrated <R>", kRed + 1,  21, false},
        {"hIntegratedRBkg",  "Background","Integrated <R>", kBlue + 1, 22, false}
    };

    std::vector<std::string> labels;
    std::vector<ProfileBundle> bundles;
    std::vector<TH1*> toDelete;

    for (size_t i = 0; i < components.size(); ++i) {
        SigExtractPoint p = FetchSigExtractPoint(cache, sigExtractDir, fileSuffix, cutFolder,
                                                 proxy.folderName, components[i].histName);
        if (!p.found) return; // All three come from one extraction: if one is absent, so are the rest
        labels.push_back(components[i].legendLabel);
    }

    const int nCats = static_cast<int>(components.size());
    for (int i = 0; i < nCats; ++i) {
        SigExtractPoint p = FetchSigExtractPoint(cache, sigExtractDir, fileSuffix, cutFolder,
                                                 proxy.folderName, components[i].histName);
        TH1D* h = MakeCategoricalPoint(Form("SigComp_%s_%d", proxy.folderName.c_str(), i),
                                       i, nCats, p.value, p.error);
        VariationConfig vc{"", components[i].legendLabel, components[i].color, 1,
                           components[i].markerStyle, components[i].isData};
        bundles.push_back({h, vc});
        toDelete.push_back(h);
    }

    DrawIntegratedCanvas(bundles, labels, "Canvas_Components",
                         familyName + " " + proxy.legendLabel + ": measured, signal and background",
                         outDir, "Integrated <R>", false, true, true);

    // --- Subtracted companion ------------------------------------------------------------------
    // The measured bin is left empty rather than removed, so the axis and labels stay identical to
    // the canvas above it. Same convention as the existing Integrated_Subtracted plots.
    const std::vector<std::pair<std::string, int>> diffs = {
        {"hIntegratedDiffMeasMinusSig", 1},
        {"hIntegratedDiffMeasMinusBkg", 2}
    };

    std::vector<ProfileBundle> subBundles;
    for (const auto& d : diffs) {
        SigExtractPoint p = FetchSigExtractPoint(cache, sigExtractDir, fileSuffix, cutFolder,
                                                 proxy.folderName, d.first);
        if (!p.found) continue;
        TH1D* h = MakeCategoricalPoint(Form("SigCompSub_%s_%d", proxy.folderName.c_str(), d.second),
                                       d.second, nCats, p.value, p.error);
        subBundles.push_back({h, {"", components[d.second].legendLabel, components[d.second].color,
                                  1, components[d.second].markerStyle, false}});
        toDelete.push_back(h);
    }

    if (!subBundles.empty()) {
        DrawIntegratedCanvas(subBundles, labels, "Canvas_Components_Subtracted",
                             familyName + " " + proxy.legendLabel +
                                 ": measured minus each component",
                             outDir, "#Delta<R> (Measured - Component)", true, true, true);
    }

    for (auto q : toDelete) delete q;
}

/**
 * @brief Builds the systematics canvases for one signal-extracted quantity.
 *
 * Directly imitates Integrated_summary: one categorical axis whose entries are Data followed by
 * every systematic variation, drawn from the extraction file belonging to each variation.
 *
 * The subtracted companion is Y = <quantity>_data - <quantity>_variation, computed here rather
 * than read from file. That IS legitimate quadrature, unlike the component differences above:
 * these come from statistically independent consumer configurations, not from shared primitives.
 * The one caveat is that the proxy-distortion variations reuse the same underlying candidates, so
 * the difference errors are conservative for those -- the same approximation the rest of this
 * macro already makes for every other systematic difference.
 *
 * @return The Data value, so the caller can accumulate it for the cross-family comparison.
 */
SigExtractPoint DrawSigExtractSystematics(FileCache& cache,
                                          const std::string& sigExtractDir,
                                          const std::string& dataSuffix,
                                          const std::vector<VariationConfig>& sysVariations,
                                          const VariationConfig& dataConfig,
                                          const std::string& cutFolder,
                                          const SigExtractProxy& proxy,
                                          const SigExtractQuantity& quantity,
                                          const std::string& familyName,
                                          TDirectory* outDir) {
    SigExtractPoint dataPoint;
    if (!outDir) return dataPoint;

    dataPoint = FetchSigExtractPoint(cache, sigExtractDir, dataSuffix, cutFolder,
                                     proxy.folderName, quantity.histName);
    if (!dataPoint.found) return dataPoint;

    std::vector<std::string> labels{dataConfig.legendLabel};
    std::vector<VariationConfig> configs{dataConfig};
    std::vector<SigExtractPoint> points{dataPoint};

    for (const auto& sys : sysVariations) {
        SigExtractPoint p = FetchSigExtractPoint(cache, sigExtractDir, dataSuffix + sys.suffix,
                                                 cutFolder, proxy.folderName, quantity.histName);
        if (!p.found) continue; // A variation whose extraction failed drops off the axis entirely
        labels.push_back(sys.legendLabel);
        configs.push_back(sys);
        points.push_back(p);
    }

    if (labels.size() < 2) return dataPoint; // Data alone is not a systematics plot

    const int nCats = static_cast<int>(labels.size());
    std::vector<ProfileBundle> bVal, bSub;
    std::vector<TH1*> toDelete;

    for (int i = 0; i < nCats; ++i) {
        TH1D* hV = MakeCategoricalPoint(Form("SigSysVal_%s_%s_%d", proxy.folderName.c_str(),
                                             quantity.histName.c_str(), i),
                                        i, nCats, points[i].value, points[i].error);
        bVal.push_back({hV, configs[i]});
        toDelete.push_back(hV);

        if (i == 0) continue; // Data bin stays empty on the difference canvas, as elsewhere

        const double diff = points[0].value - points[i].value;
        const double ediff = std::sqrt(points[0].error * points[0].error +
                                       points[i].error * points[i].error);
        TH1D* hS = MakeCategoricalPoint(Form("SigSysSub_%s_%s_%d", proxy.folderName.c_str(),
                                             quantity.histName.c_str(), i),
                                        i, nCats, diff, ediff);
        bSub.push_back({hS, configs[i]});
        toDelete.push_back(hS);
    }

    DrawIntegratedCanvas(bVal, labels, "Canvas_Integrated",
                         familyName + " " + proxy.legendLabel + " " + quantity.legendLabel,
                         outDir, quantity.yAxisTitle, false, false, true);

    if (!bSub.empty()) {
        DrawIntegratedCanvas(bSub, labels, "Canvas_Integrated_Subtracted",
                             familyName + " " + proxy.legendLabel + " " + quantity.legendLabel +
                                 " difference",
                             outDir, "#Delta" + quantity.yAxisTitle + " (Data - Var)",
                             true, false, true);
    }

    for (auto q : toDelete) delete q;
    return dataPoint;
}

void auxiliarySummaryPlots(const std::string& consumerDir,
                           const std::string& mcRefDir = "",
                           const std::string& ppRefDir = "",
                           const std::string& toyModelPath = "",
                           const std::string& cutFolder = DEFAULT_CUT_FOLDER,
                           bool doIndividualComparisons = false,
                           const std::string& sigExtractDir = "") {
    
    // 1. Define Systematic Variations (the list of all useful variations I would like to track into this plot)
    // The data config (empty suffix) is handled separately in the logic to ensure it is always first
    // Field order after isData: inRedux, inPtGateStudy, inEtaGateStudy.
    std::vector<VariationConfig> sysVariations = {
        // suffix                          legend                      colour     ls mk  data   redux  ptGate etaGate
        {"_forceRandJet",                 "Rand Jet",                  kBlue,     1, 20, false, false,  true,  true}, // Full circle
        // Data-like is off the reduced canvas because it is the crudest of the mixing proxies: a
        // cheap event mixing, cheaper even than prevJet, which does not carry the correlations
        // between the angular distributions in full. It remains on the full systematics canvas and
        // in the pT gate study.
        {"_forceDatalikeJet",             "Data-like Jet",             kRed,      1, 21, false, false, true,  false}, // Square
        {"_forcePerpToJet",               "Perp to Jet",               kGreen+2,  1, 22, false, true,  true,  true},  // Triangle up
        {"_forcePreviousJet",             "Prev Jet",                  kOrange+1, 1, 34, false, false, false, false}, // Full cross (bold "+")
        {"_MixedEventProxies",            "MixedEv Jet",               kGray+2,   1, 24, false, true,  false, false}, // Open circle
        {"_MixedEventProxies_10ColGap",   "MixedEv Jet (10coll gap)",  kGray+1,   2, 25, false, false, false, false}, // Open circle, dashed line

        // Do the artificial proxies survive a minimum-pT gate? The pT of an invented direction is rebuilt from
        // the original |p| (see gatePtOnArtificialProxies in applyProxyDistortion), so gating on it is a genuine
        // selection rather than a no-op, and it may sculpt the sample. Same colour as the ungated parent, dashed
        // line and the open counterpart of its marker, so the pairs read together on a shared canvas.
        // The trailing "true" marks both members of each pair as part of the gate study, which gets its own
        // canvas (Canvas_PtGateStudy) free of the unrelated systematics.
        {"_forceRandJet_gatePtOnProxy",     "Rand Jet (pT gated)",      kBlue,     2, 24, false, false, true,  false}, // Open circle
        {"_forcePerpToJet_gatePtOnProxy",   "Perp to Jet (pT gated)",   kGreen+2,  2, 26, false, false, true,  false}, // Open triangle up
        {"_forceDatalikeJet_gatePtOnProxy", "Data-like Jet (pT gated)", kRed,      2, 25, false, false, true,  false}, // Open square

        // Do the artificial proxies still make sense WITHOUT the eta acceptance gate? Unlike the pT
        // study above, gating on eta is now the default: an invented direction that ignores the
        // experimental acceptance produces an observable the data could never have produced, so the
        // question is how much is lost by not gating. Only RandJet and PerpToJet exist in both forms;
        // the data-like proxy gates on eta by construction and has no ungated counterpart.
        // Same colour as the gated parent, dotted line and an open marker, so the pairs read together.
        {"_forceRandJet_noGateEtaOnProxy",   "Rand Jet (no #eta gate)",  kBlue,     3, 32, false, false, false, true}, // Open triangle down
        {"_forcePerpToJet_noGateEtaOnProxy", "Perp to Jet (no #eta gate)", kGreen+2, 3, 28, false, false, false, true}, // Open cross

        // Tighter primary-vertex selection. Unrelated to the proxy machinery, so it gets its own colour.
        {"_5cmZvtx",                      "|Zvtx| < 5 cm",             kBlack,    1, 24, false, true,  false, false} // Open circle, similar to dataConfig, as this is also data.

        // Kept for reference: the resampled variations were superseded by the gated ones above.
        // Uncomment (and re-run the corresponding consumer configs) if they are ever needed again.
        // {"_forceDatalikeJet_10resamples", "Data-like Jet (10 resam.)", kRed,     2, 25, false}, // Open square
        // {"_forcePerpToJet_10resamples",   "Perp to Jet (10 resam.)",   kGreen+2, 2, 26, false}  // Open triangle up
    };

    VariationConfig dataConfig = {"", "Data", kBlack,  1, 8, true}; // Thick black line. Thickness is controlled by the flag turned "true", essentially
    VariationConfig mcConfig   = {"", "OO MC Base", kCyan+1, 2, 29, false}; // Full stars for MC, with a darker cyan type of color (see TColor and TAttMarker)
                                                                             // Kept line as dashed because MC error bars are (VERY!) large with current statistics
    VariationConfig ppConfig   = {"", "pp Data", kAzure+2, 2, 47, false}; // X shape, something between blue and cyan
    VariationConfig toyConfig  = {"", "Toy Model", kViolet+1, 2, 33, false}; // Diamond shape, purple-ish

    // External reference datasets, in the order they should appear on the categorical axes.
    // An empty basePath simply drops that reference from every plot, with no index bookkeeping needed.
    std::vector<ExternalRef> externals = {
        {mcConfig,  mcRefDir,      ExternalRef::kConsumerDir},
        {ppConfig,  ppRefDir,      ExternalRef::kConsumerDir},
        {toyConfig, toyModelPath,  ExternalRef::kToyModel}
    };

    // Define Mass Selection Variations (treated as a separate standalone family to avoid polluting main plots)
    // (We are essentially checking to see if the eta dependency of the ring, the integrated observables and the folded spectra
    //  change whenever considering a ring observable calculated using particles strictly inside a mass window or strictly outside
    //  that window, defined by the v0InMassPeak bool in the consumer)
        // This is a naive alternative to signal extraction!
    std::vector<VariationConfig> massVariations = {
        {"_excludeOutOfPeak", "In Mass Peak",     kTeal+3,    1, 34, false}, // Cross
        {"_excludeInPeak",    "Out Of Mass Peak", kMagenta+2, 1, 23, false}  // Down triangle
    };

    // 2. Define the Families
    std::vector<FamilyConfig> families = {
        {"Lambda",       "JustLambda"},
        {"AntiLambda",   "JustAntiLambda"},
        {"BothHyperons", "BothHyperons"} // BothHyperons_2GeVLeadP was already the standard LeadP pT, so I modified this to follow the usual
                                         // naming scheme that was in place for other variations of BothHyperons configurables of the consumer
    };

    // 3. Define the Observable Groups
    // Split by what the x-axis actually is: the proxy's own eta, or the hyperon's eta.
    //
    // TODO: the consumer's profile naming is not uniform across observables. The leading jet has no infix
    // in some families (pRingObservableDeltaPhi) but does in others (pRingObservableLeadJetPVz), and the
    // subleading jet is "2ndJet" everywhere except PVz, where it is "SubLead". That is why every name below
    // is spelled out instead of being built from a proxy infix. Worth unifying in the consumer at some
    // point, but it would break other downstream code, so leaving it alone for now.
    std::vector<ObservableGroup> observableGroups = {
        {
            "EtaProxy", // x-axis is the eta of the jet/leading-particle proxy
            {
                {"EtaDependence/", "pRingObservableEtaLeadP", "#eta_{LeadP}", "R", true, true, "pRingProxyJetVsEtaJet", false, "LeadP"},
                {"EtaDependence/", "pRingObservableEtaJet", "#eta_{Jet}", "R", true, false, "pRingProxyJetVsEtaJet", false, "LeadJet"},
                {"EtaDependence/", "pRingObservableEta2ndJet", "#eta_{2ndJet}", "R", true, false, "pRingProxyJetVsEtaJet", false, "SubJet"},
                // Fine-binned counterparts of the two above. Same axis, finer granularity: these are the
                // ones that show the tanh-like shape cleanly in the Toy Model.
                {"EtaDependence/", "pRingObservableEtaJetHighEtaRes", "#eta_{Jet}", "R", true, false, "pRingProxyJetVsEtaJet", true, "LeadJet"},
                {"EtaDependence/", "pRingObservableEtaLeadPHighEtaRes", "#eta_{LeadP}", "R", true, true, "pRingProxyJetVsEtaJet", true, "LeadP"}
            },
            true // Eta axis is symmetric about zero, so folding is meaningful
        },
        {
            "EtaV0", // x-axis is the eta of the hyperon itself, split by which proxy defined the ring
            {
                {"EtaDependence/", "pRingObservableEtaLambda", "#eta_{#Lambda}", "R", true, false, "pRingProxyJetVsEta", false, "LeadJet"},
                {"EtaDependence/", "pRingObservableEtaLambdaLeadP", "#eta_{#Lambda(LeadP)}", "R", true, true, "pRingProxyJetVsEta", false, "LeadP"},
                {"EtaDependence/", "pRingObservableEtaLambda2ndJet", "#eta_{#Lambda(2ndJet)}", "R", true, false, "pRingProxyJetVsEta", false, "SubJet"}
            },
            true
        },
        {
            // <R> against the angle between the hyperon and its decay proton in the helicity frame.
            // Booked at TASK level rather than inside a cut folder, hence underCutFolder = false.
            // <R> changes sign across this axis, so folding it would cancel the very structure the
            // observable exists to show: doFold stays false.
            // README: the Lambda-like set is the QA for species competition -- with balanced Lambda
            // and AntiLambda yields its dependence should cancel, so a residual there measures the
            // imbalance rather than a physics asymmetry.
            "AEE",
            {
                {"HelicityEfficiencyQA/PhiLambdaPhiProtonStar/", "pRingObservableLeadJetVsPhiLambdaLikePhiProtonStar", "#phi_{#Lambda-like}-#phi_{p-like}^{*}", "R", false, false, "", false, "LeadJet", "LambdaLike", "RingLeadJetAEE_LambdaLike"},
                {"HelicityEfficiencyQA/PhiLambdaPhiProtonStar/", "pRingObservableLeadPVsPhiLambdaLikePhiProtonStar",   "#phi_{#Lambda-like}-#phi_{p-like}^{*}", "R", false, true,  "", false, "LeadP", "LambdaLike", "RingLeadPAEE_LambdaLike"},
                {"HelicityEfficiencyQA/PhiLambdaPhiProtonStar/", "pRingObservable2ndJetVsPhiLambdaLikePhiProtonStar",  "#phi_{#Lambda-like}-#phi_{p-like}^{*}", "R", false, false, "", false, "SubJet", "LambdaLike", "RingSubJetAEE_LambdaLike"},

                {"HelicityEfficiencyQA/PhiLambdaPhiProtonStar/", "pRingObservableLeadJetVsPhiLambdaPhiProtonStar", "#phi_{#Lambda}-#phi_{p}^{*}", "R", false, false, "", false, "LeadJet", "Lambda", "RingLeadJetAEE_Lambda"},
                {"HelicityEfficiencyQA/PhiLambdaPhiProtonStar/", "pRingObservableLeadPVsPhiLambdaPhiProtonStar",   "#phi_{#Lambda}-#phi_{p}^{*}", "R", false, true,  "", false, "LeadP", "Lambda", "RingLeadPAEE_Lambda"},
                {"HelicityEfficiencyQA/PhiLambdaPhiProtonStar/", "pRingObservable2ndJetVsPhiLambdaPhiProtonStar",  "#phi_{#Lambda}-#phi_{p}^{*}", "R", false, false, "", false, "SubJet", "Lambda", "RingSubJetAEE_Lambda"},

                {"HelicityEfficiencyQA/PhiLambdaPhiProtonStar/", "pRingObservableLeadJetVsPhiAntiLambdaPhiProtonStar", "#phi_{#bar{#Lambda}}-#phi_{#bar{p}}^{*}", "R", false, false, "", false, "LeadJet", "AntiLambda", "RingLeadJetAEE_AntiLambda"},
                {"HelicityEfficiencyQA/PhiLambdaPhiProtonStar/", "pRingObservableLeadPVsPhiAntiLambdaPhiProtonStar",   "#phi_{#bar{#Lambda}}-#phi_{#bar{p}}^{*}", "R", false, true,  "", false, "LeadP", "AntiLambda", "RingLeadPAEE_AntiLambda"},
                {"HelicityEfficiencyQA/PhiLambdaPhiProtonStar/", "pRingObservable2ndJetVsPhiAntiLambdaPhiProtonStar",  "#phi_{#bar{#Lambda}}-#phi_{#bar{p}}^{*}", "R", false, false, "", false, "SubJet", "AntiLambda", "RingSubJetAEE_AntiLambda"}
            },
            false
        },
        {
            "ProxyPt", // <R> against the pT of whichever proxy defined the ring
            //
            // No Toy Model counterpart, and not a cheap one to add: the toy's jet is a pure unit
            // direction built from (cosThetaJet, phiJet) and carries no magnitude at all, so there is
            // no "proxy pT" to bin against. Supplying one would mean inventing a jet momentum spectrum,
            // which is a modelling decision rather than a booking change. Deliberately left empty.
            {
                {"ProxyPtDependence/", "pRingVsPtJet", "#it{p}_{T}^{Jet} (GeV/c)", "R", true, false, "", false, "LeadJet"},
                {"ProxyPtDependence/", "pRingVsPtLeadP", "#it{p}_{T}^{LeadP} (GeV/c)", "R", true, true, "", false, "LeadP"},
                {"ProxyPtDependence/", "pRingVsPt2ndJet", "#it{p}_{T}^{SubJet} (GeV/c)", "R", true, false, "", false, "SubJet"}
            },
            false // pT is strictly positive: folding about zero is meaningless
        },
        {
            "ProxyPt/EtaSplits", // The LeadP pT dependence, split by the sign of the proxy's or the hyperon's eta.
                                // Kept to the four single-sign splits: the four combined (proxy AND V0) splits exist
                                // in the consumer too, but they multiply the plot count for little extra information.
                                // No toy counterpart, for the same reason as ProxyPt above.
            {
                {"ProxyPtDependence/", "pRingVsPtLeadP_PosEtaLeadP", "#it{p}_{T}^{LeadP} (GeV/c), #eta_{LeadP}>0", "R", true, true},
                {"ProxyPtDependence/", "pRingVsPtLeadP_NegEtaLeadP", "#it{p}_{T}^{LeadP} (GeV/c), #eta_{LeadP}<0", "R", true, true},
                {"ProxyPtDependence/", "pRingVsPtLeadP_PosEtaV0",    "#it{p}_{T}^{LeadP} (GeV/c), #eta_{V0}>0",    "R", true, true},
                {"ProxyPtDependence/", "pRingVsPtLeadP_NegEtaV0",    "#it{p}_{T}^{LeadP} (GeV/c), #eta_{V0}<0",    "R", true, true}
            },
            false
        },
        {
            "PVz", // Primary-vertex z dependence, a pure detector-geometry handle on the fake signal
                   // No toy counterpart by design: the toy has no primary vertex to displace.
            {
                {"", "pRingObservableLeadJetPVz", "PVz (cm), LeadJet", "R", true, false, "", false, "LeadJet"},
                {"", "pRingObservableLeadPPVz", "PVz (cm), LeadP", "R", true, true, "", false, "LeadP"},
                {"", "pRingObservableSubLeadPVz", "PVz (cm), SubJet", "R", true, false, "", false, "SubJet"}
            },
            true // PVz is symmetric about zero, so folding it is meaningful
        },
        {
            "QA_Mult", // Event-activity dependence. Ties directly to the ratio-estimator discussion
                       // in the "Ring uncertainty estimators" note: the ring denominator is itself
                       // a fluctuating per-event multiplicity.
                       // No toy counterpart: the toy has neither a centrality nor a multiplicity concept.
            {
                {"", "pRingVsCentrality", "Centrality (%)", "R", true, false, "", false, "LeadJet"},
                // Booked at task level rather than inside the kinematic-cut folders, hence underCutFolder = false
                {"IntegratedCuts/", "pRingVsNV0s", "N_{#Lambda}+N_{#bar{#Lambda}}", "R", false, false, "", false, "LeadJet"}
            },
            false
        },
        {
            "OtherAngDepndncs", // Lab-frame azimuth. Relevant to the acceptance/AEE story rather than to the ring itself
            //
            // TODO (Toy Model): <R> vs phi_Lambda would be the cheapest addition of the three flagged here.
            // phi_Lambda is already being used as an axis by pPstarX/Y/Z_vsPhiLam, so the binning and the
            // variable both exist; only a <R>-weighted TProfile against it is missing.
            {
                {"", "pRingObservablePhiJet", "#varphi_{Jet}", "R", true, false, "", false, "LeadJet"},
                {"", "pRingObservablePhiLambda", "#varphi_{#Lambda}", "R", true, false, "", false, "LeadJet"}
            },
            false // Lab-frame phi runs over [0, 2pi): there is no zero to fold about
        },
        {
            "OtherAngDepndncs/DeltaPhi", // Azimuthal separation between the hyperon and the proxy
            //
            // TODO (Toy Model): a <R> vs (phi_Lambda - phi_jet) TProfile would be a cheap addition. Both
            // angles are already in scope where ringProxyJet is computed (phiJet from the jet reshuffle
            // block, phi_Lambda from the hyperon kinematics), so it is a booking plus a single Fill, with
            // no new sampling and no extra runtime. Fill toyProfile below once it exists.
            {
                {"", "pRingObservableDeltaPhi", "#Delta#varphi_{jet}", "R", true, false, "", false, "LeadJet"},
                {"", "pRingObservableLeadPDeltaPhi", "#Delta#varphi_{LeadP}", "R", true, true, "", false, "LeadP"},
                {"", "pRingObservable2ndJetDeltaPhi", "#Delta#varphi_{SubJet}", "R", true, false, "", false, "SubJet"}
            },
            false // DeltaPhi runs over a full period, not symmetric about zero
        },
        {
            "OtherAngDepndncs/DeltaTheta", // Polar separation between the hyperon and the proxy
            //
            // TODO (Toy Model): same situation as DeltaPhi. cosThetaJet and the hyperon direction are both
            // already available at the ring computation, so <R> vs DeltaTheta is a booking plus one Fill.
            {
                {"", "pRingObservableDeltaTheta", "#Delta#theta_{jet}", "R", true, false, "", false, "LeadJet"},
                {"", "pRingObservableLeadPDeltaTheta", "#Delta#theta_{LeadP}", "R", true, true, "", false, "LeadP"},
                {"", "pRingObservable2ndJetDeltaTheta", "#Delta#theta_{SubJet}", "R", true, false, "", false, "SubJet"}
            },
            false
        }
    };

    // Which lower-pad variants to emit alongside every comparison canvas. Both are generated in the
    // same pass over the same bundles, so each overlay gains a "_Diff" and a "_Pull" sibling.
    // Emptying this vector reverts the macro to single-pad canvases only.
    std::vector<LowerPadMode> lowerPadModes = {kDifference, kPull};

    // Used where a full-size standalone difference canvas already exists for exactly the same bundles
    // (the Subtracted/ folder, and the per-variation Canvas_Subtracted* in Individual_Comparisons).
    // The pull has no standalone equivalent anywhere, so it is always worth generating.
    const std::vector<LowerPadMode> kPullOnly = {kPull};

    // 3b. Define the categorical (labelled-bin) cut scans
    // These are the eta-split closure tests: each bin is a named selection rather than a coordinate,
    // so they need the labelled-axis drawing path. Booked at task level and only when the consumer's
    // doFakePolDiagnosticsQA is on, hence fetched defensively like everything else.
    std::vector<CategoricalObservable> categoricalObservables = {
        {"EtaStudy/", "pRingEtaCuts",              "R", true, ""},
        {"EtaStudy/", "pRingEtaCutsSubLeadingJet", "R", true, ""},
        {"EtaStudy/", "pRingEtaCutsLeadingP",      "R", true, "pRingEtaCutsLeadingP_MassSignalVsBackground"}
    };

    // 4. Create the output file
    std::string outFilePath = consumerDir + "/auxiliarySummaryPlots.root";
    TFile* fOut = new TFile(outFilePath.c_str(), "RECREATE");
    if (!fOut || fOut->IsZombie()) {
        std::cerr << "[Error] Could not create output file: " << outFilePath << std::endl;
        return;
    }

    // One cache for the whole run: every consumer/MC/pp/toy file is opened at most once
    FileCache cache;

    std::cout << "========================================================\n";
    std::cout << " Generating Auxiliary Plots\n";
    std::cout << " Output:     " << outFilePath << "\n";
    std::cout << " Cut folder: " << cutFolder << "\n";
    std::cout << "========================================================\n";

    // One representative observable per proxy. Integrated values are a property of the entry set, not of
    // the binning, so a single observable per proxy captures everything the others would repeat.
    const std::vector<std::pair<std::string, std::string>> proxyRepresentatives = {
        {"LeadP",   "pRingObservableEtaLeadP"},
        {"LeadJet", "pRingObservableEtaJet"},
        {"SubJet",  "pRingObservableEta2ndJet"}
    };

    // Styling for the cross-family comparison, where the SERIES are the hyperon selections
    const std::vector<VariationConfig> familyStyles = {
        {"", "Lambda",       kBlue,  1, 20, false},
        {"", "AntiLambda",   kRed,   1, 21, false},
        {"", "BothHyperons", kBlack, 1, 8,  true}
    };

    // The three proxy folders written by signalExtractionRing, in the same order as
    // proxyRepresentatives above so the two cross-family blocks line up.
    const std::vector<SigExtractProxy> sigExtractProxies = {
        {"LeadP",   "Leading particle"},
        {"LeadJet", "Leading jet"},
        {"SubJet",  "Subleading jet"}
    };

    // The two signal-extracted quantities that get a full systematics treatment. <R>_measured is
    // deliberately NOT here: it already has a counterpart in the ordinary Integrated_summary plots,
    // and repeating it would invite the two to be read as independent measurements.
    const std::vector<SigExtractQuantity> sigExtractQuantities = {
        {"hIntegratedRSig", "Signal",     "Integrated <R>_{S}", kRed + 1,  21, false},
        {"hIntegratedRBkg", "Background", "Integrated <R>_{B}", kBlue + 1, 22, false}
    };

    // sigExtractCrossFamily[quantity][proxy][family] = the Data value from each extraction file
    std::vector<std::vector<std::vector<SigExtractPoint>>> sigExtractCrossFamily(
        sigExtractQuantities.size(),
        std::vector<std::vector<SigExtractPoint>>(sigExtractProxies.size(),
                                                  std::vector<SigExtractPoint>(families.size())));

    // crossFamilyVals[proxy][family] = one integrated value per variation, in allSystematics order
    std::vector<std::vector<std::vector<std::pair<double, double>>>>
        crossFamilyVals(proxyRepresentatives.size(),
                        std::vector<std::vector<std::pair<double, double>>>(families.size()));
    std::vector<std::vector<std::string>> crossFamilyLabels(proxyRepresentatives.size());

    // 5. Main Processing Loops
    for (size_t famIdx = 0; famIdx < families.size(); ++famIdx) {
        const auto& fam = families[famIdx];
        std::cout << " -> Processing Family: " << fam.familyName << std::endl;

        // Directories are created lazily (only once something is actually written into them),
        // so a consumer file that does not enable the requested cut folder leaves no empty folders.
        TDirectory* famDir = nullptr;

        // Variables to accumulate the integrated observable summary across all variations for this family
        // (a "grand summary" for short)
        std::vector<std::pair<double, double>> grandDataVals;
        std::vector<std::vector<std::pair<double, double>>> grandExtVals(externals.size());  // One row per external reference
        std::vector<std::vector<std::pair<double, double>>> grandSysVals(sysVariations.size()); // Track all systematic variations
        std::vector<std::vector<std::pair<double, double>>> grandMassVals(massVariations.size());  // Mass-window selections
        std::vector<std::string> grandLabels;
        std::vector<std::string> grandProxyTags; // Parallel to grandLabels, used to group the axis by proxy

        // EtaSplitStudies is emitted right after EtaProxy rather than at the end of the family, because
        // ROOT lists directories in creation order and this is the browsing order that was asked for:
        // EtaProxy, EtaSplitStudies, EtaV0, AEE, HEE, CheapSigExtract, ProxyPt, PVz, QA_Mult,
        // OtherAngDepndncs, BruteForce. The eta groups lead because they carry the main result; the
        // AEE and HEE blocks follow as the acceptance cross-checks on it.
        bool etaSplitStudiesDone = false;
        // Emitted as a lambda, like runEtaSplitStudies, so the AEE and HEE folders can be
        // written from inside the group loop. Directory order in a ROOT file is creation
        // order, so a block that runs after the loop lands after every group -- which is how
        // HEE ended up between OtherAngDepndncs and BruteForce instead of beside AEE.
        bool aeeHeeDone = false;
        auto runAeeHeeBlocks = [&]() {
            if (aeeHeeDone) return;
            aeeHeeDone = true;
            // -----------------------------------------------------------------------------------
            // Mass-slice cross-checks: <R> in the peak region against <R> in the two sidebands
            // -----------------------------------------------------------------------------------
            // Task-level objects, so they carry no cutFolder and are read once per family. Both groups
            // use the same three-bin mass axis and the same drawer; only the x-axis differs.
            // README: these exist to be independent of signalExtractionRing, so that a disagreement
            // between the two is informative rather than circular.
            {
                const std::string famFile = consumerDir + "/ConsumerResults_" + fam.dataSuffix + ".root";
                const std::string aeeDir = TASK_DIR_IN_FILE + "HelicityEfficiencyQA/";
            
                struct MassSliceSpec { std::string path, outName, xTitle, title, species; };
            
                // (a) <R> against phi_Lambda - phi_p*, three mass regions, per proxy and per species.
                std::vector<MassSliceSpec> phiSpecs;
                const std::pair<const char*, const char*> proxies[3] = {
                    {"LeadJet", "<#it{R}>_{LeadJet}"},
                    {"LeadP",   "<#it{R}>_{LeadP}"},
                    {"2ndJet",  "<#it{R}>_{SubJet}"}
                };
                const std::pair<const char*, const char*> species[3] = {
                    {"PhiLambdaLike",  "#phi_{#Lambda-like}-#phi_{p-like}^{*}"},
                    {"PhiLambda",      "#phi_{#Lambda}-#phi_{p}^{*}"},
                    {"PhiAntiLambda",  "#phi_{#bar{#Lambda}}-#phi_{#bar{p}}^{*}"}
                };
                for (const auto& pr : proxies) {
                    for (const auto& sp : species) {
                        // Same species restriction as the 1D group above
                        const std::string spTag = (std::string(sp.first) == "PhiLambda")     ? "Lambda"
                                                : (std::string(sp.first) == "PhiAntiLambda") ? "AntiLambda"
                                                                                             : "LambdaLike";
                        if (!SpeciesMatchesFamily(spTag, fam.familyName)) continue;
                        // Destination folder is the 1D observable's own folder, so the cheap
                        // check sits next to the distribution it is checking rather than in a
                        // separate tree the reader has to cross-reference.
                        // README: displayName here must match the AEE table entry exactly.
                        const std::string proxyDisp = (std::string(pr.first) == "2ndJet") ? "SubJet"
                                                                                          : pr.first;
                        phiSpecs.push_back({
                            aeeDir + "PhiLambdaPhiProtonStar/ThreeBinMass/p2dRingObservable" +
                                pr.first + "Vs" + sp.first + "PhiProtonStarVs3BinMass",
                            "Ring" + proxyDisp + "AEE_" + spTag,
                            sp.second,
                            std::string(pr.second) + " vs " + sp.second,
                            spTag});
                    }
                }
            
                for (const auto& spec : phiSpecs) {
                    TProfile2D* p2d = cache.FetchClone<TProfile2D>(famFile, spec.path, false);
                    if (!p2d) continue; // Absent for consumer output produced before these were booked
                    // spec.outName is the observable's displayName, so this lands in
                    // <Family>/AEE/<observable>/Mass_Selection/
                    TDirectory* obsMassDir = EnsureDir(
                        EnsureDir(EnsureDir(EnsureDir(fOut, fam.familyName), "AEE"), spec.outName),
                        "Mass_Selection");
                    DrawMassRegionCanvases(p2d, fam.familyName + " " + spec.title,
                                           obsMassDir, spec.xTitle, "R");
                    delete p2d;
                }
            
                // (b) <R> against cos(theta) in the helicity frame. Leading jet only, which is why the
                //     titles say so explicitly: there is no LeadP or SubJet counterpart to compare with.
                //     Motivated by the angular dependence reported in J. Adams' global-polarization
                //     thesis; the question is whether the ring observable shows anything similar.
                const std::vector<MassSliceSpec> heeSpecs = {
                    {aeeDir + "CosThetaHEEByMass/p2dRingLambdaCosThetaHEEVsMass", "Canvas_HEE_Lambda",
                     "cos(#theta)_{HEE} = #hat{p}^{*}_{D} . #vec{p}_{#Lambda}",
                     "<#it{R}>_{LeadJet} vs cos(#theta)_{HEE}, #Lambda, by mass region", "Lambda"},
                    {aeeDir + "CosThetaHEEByMass/p2dRingAntiLambdaCosThetaHEEVsMass", "Canvas_HEE_AntiLambda",
                     "cos(#theta)_{HEE} = #hat{p}^{*}_{D} . #vec{p}_{#bar{#Lambda}}",
                     "<#it{R}>_{LeadJet} vs cos(#theta)_{HEE}, #bar{#Lambda}, by mass region", "AntiLambda"},
                    {aeeDir + "CosThetaHEEByMass/p2dRingAllV0sCosThetaHEEVsMass", "Canvas_HEE_AllV0s",
                     "cos(#theta)_{HEE} = #hat{p}^{*}_{D} . #vec{p}_{#Lambda-like}",
                     "<#it{R}>_{LeadJet} vs cos(#theta)_{HEE}, all V0s, by mass region", "LambdaLike"}
                };
            
                TDirectory* heeDir = nullptr;
            
                // Is the HEE dependence a reflection of the AEE angle, or a separate effect?
                // Splitting on the sign of phi_Lambda - phi_p* answers it directly: if the cos(theta)
                // trend is driven by the AEE, the two halves separate and the inclusive curve sits
                // between them. If it survives with the same shape in both halves, it is its own effect.
                // README: this is the whole point of the AEE-split canvases below.
                struct HeeSplitSpec { std::string base, outName, xTitle, species; };
                const std::vector<HeeSplitSpec> heeSplits = {
                    {aeeDir + "CosThetaHEEByMass/p2dRingLambdaCosThetaHEEVsMass", "Lambda",
                     "cos(#theta)_{HEE} = #hat{p}^{*}_{D} . #vec{p}_{#Lambda}", "Lambda"},
                    {aeeDir + "CosThetaHEEByMass/p2dRingAntiLambdaCosThetaHEEVsMass", "AntiLambda",
                     "cos(#theta)_{HEE} = #hat{p}^{*}_{D} . #vec{p}_{#bar{#Lambda}}", "AntiLambda"},
                    {aeeDir + "CosThetaHEEByMass/p2dRingAllV0sCosThetaHEEVsMass", "AllV0s",
                     "cos(#theta)_{HEE} = #hat{p}^{*}_{D} . #vec{p}_{#Lambda-like}", "LambdaLike"}
                };
            
                for (const auto& hs : heeSplits) {
                    if (!SpeciesMatchesFamily(hs.species, fam.familyName)) continue;
            
                    TProfile2D* pPos = cache.FetchClone<TProfile2D>(famFile, hs.base + "PosDeltaPhiAEE", false);
                    TProfile2D* pNeg = cache.FetchClone<TProfile2D>(famFile, hs.base + "NegDeltaPhiAEE", false);
                    if (!pPos || !pNeg) { delete pPos; delete pNeg; continue; }
            
                    if (!heeDir) heeDir = EnsureDir(EnsureDir(fOut, fam.familyName), "HEE");
            
                    // The inclusive curve is REBUILT from the two halves rather than read from the
                    // unsplit profile. TProfile::Add merges sums and entry counts, so the result is the
                    // exact entry-weighted inclusive -- and it is guaranteed to share binning with the
                    // two halves, which a separately booked profile is not.
                    const int nMass = pPos->GetYaxis()->GetNbins();
                    const int peakBin = (nMass == 3) ? 2 : 1; // Middle slice is the in-window region
            
                    auto sliceOf = [&](TProfile2D* p, int m, const char* tag) {
                        TProfile* s = p->ProfileX(Form("hee_%s_%s_%d", hs.outName.c_str(), tag, m), m, m);
                        if (s) s->SetDirectory(nullptr);
                        return s;
                    };
            
                    // (a) In-window only: all, positive half, negative half
                    TProfile* sPos = sliceOf(pPos, peakBin, "pos");
                    TProfile* sNeg = sliceOf(pNeg, peakBin, "neg");
                    if (sPos && sNeg) {
                        TProfile* sAll = (TProfile*)sPos->Clone(Form("hee_%s_all", hs.outName.c_str()));
                        sAll->SetDirectory(nullptr);
                        sAll->Add(sNeg);
            
                        std::vector<ProfileBundle> hb;
                        VariationConfig cAll; cAll.legendLabel = "all #phi_{#Lambda}-#phi_{p}^{*}";
                        cAll.color = kBlack; cAll.lineStyle = 1; cAll.markerStyle = 20; cAll.isData = true;
                        VariationConfig cP; cP.legendLabel = "#phi_{#Lambda}-#phi_{p}^{*} > 0";
                        cP.color = kRed + 1; cP.lineStyle = 1; cP.markerStyle = 21; cP.isData = false;
                        VariationConfig cN; cN.legendLabel = "#phi_{#Lambda}-#phi_{p}^{*} < 0";
                        cN.color = kBlue + 1; cN.lineStyle = 1; cN.markerStyle = 22; cN.isData = false;
                        hb.push_back({sAll, cAll}); hb.push_back({sPos, cP}); hb.push_back({sNeg, cN});
            
                        ProfileConfig pc; pc.xAxisTitle = hs.xTitle; pc.yAxisTitle = "R";
                        const std::string nm = "Canvas_HEE_" + hs.outName + "_AEESplit_InPeak";
                        const std::string ttl = fam.familyName + " <#it{R}>_{LeadJet} vs cos(#theta)_{HEE}, "
                                              + hs.outName + ", in mass window, split by AEE angle sign";
                        DrawComparisonCanvas(hb, nm, ttl, heeDir, pc, hs.xTitle, "R");
                        DrawComparisonCanvas(hb, nm + "_PullOnly", ttl, heeDir, pc, hs.xTitle, "R", kPull);
                        delete sAll;
                    }
                    delete sPos; delete sNeg;
            
                    // (b) Everything: both halves across all three mass regions on one canvas
                    std::vector<ProfileBundle> hbAll;
                    std::vector<TH1*> hbDelete;
                    const int massColors[3] = {kAzure + 2, kRed + 1, kGreen + 3};
                    const char* massNames[3] = {"left sideband", "in mass window", "right sideband"};
                    for (int m = 1; m <= nMass && m <= 3; ++m) {
                        for (int side = 0; side < 2; ++side) {
                            TProfile2D* srcP = (side == 0) ? pPos : pNeg;
                            TProfile* s = sliceOf(srcP, m, side == 0 ? "P" : "N");
                            if (!s || s->GetEntries() < 1) { delete s; continue; }
                            VariationConfig c;
                            c.legendLabel = Form("%s, #Delta#phi_{AEE} %s", massNames[(m - 1) % 3],
                                                 side == 0 ? "> 0" : "< 0");
                            c.color = massColors[(m - 1) % 3];
                            c.lineStyle = (side == 0) ? 1 : 2;   // Solid positive, dashed negative
                            c.markerStyle = (side == 0) ? 20 : 24; // Full / open, same shape
                            c.isData = (m == peakBin && side == 0);
                            hbAll.push_back({s, c});
                            hbDelete.push_back(s);
                        }
                    }
                    if (hbAll.size() >= 2) {
                        ProfileConfig pc; pc.xAxisTitle = hs.xTitle; pc.yAxisTitle = "R";
                        DrawComparisonCanvas(hbAll, "Canvas_HEE_" + hs.outName + "_AEESplit_AllMass",
                                             fam.familyName + " <#it{R}>_{LeadJet} vs cos(#theta)_{HEE}, "
                                             + hs.outName + ", every mass region and AEE sign",
                                             heeDir, pc, hs.xTitle, "R");
                    }
                    for (auto p : hbDelete) delete p;
                    delete pPos; delete pNeg;
                }
            
                for (const auto& spec : heeSpecs) {
                    if (!SpeciesMatchesFamily(spec.species, fam.familyName)) continue;
                    TProfile2D* p2d = cache.FetchClone<TProfile2D>(famFile, spec.path, false);
                    if (!p2d) continue;
                    if (!heeDir) heeDir = EnsureDir(EnsureDir(fOut, fam.familyName), "HEE");
                    DrawMassSliceCanvas(p2d, spec.outName, fam.familyName + " " + spec.title,
                                        heeDir, spec.xTitle, "R");
                    delete p2d;
                }
            }
        };

        auto runEtaSplitStudies = [&]() {
            if (etaSplitStudiesDone) return;
            etaSplitStudiesDone = true;
            // ---------------------------------------------------------
            // Categorical cut scans (the eta-split closure test)
            // ---------------------------------------------------------
            // Same fetch/subtract machinery as above, but drawn on a labelled axis. The headline plot here is
            // the pull panel: after subtracting a fake-signal estimator from the data, a residual consistent
            // with zero across every eta category is what "the correction closes" actually looks like.
            for (const auto& catObs : categoricalObservables) {

                const std::string catPath = TASK_DIR_IN_FILE + catObs.subDir + catObs.profileName;
                const std::string dataFile = consumerDir + "/ConsumerResults_" + fam.dataSuffix + ".root";

                TProfile* pCatData = cache.FetchClone<TProfile>(dataFile, catPath);
                if (!pCatData) continue; // Not booked in this configuration: nothing is created on disk

                if (!famDir) famDir = EnsureDir(fOut, fam.familyName);
                TDirectory* etaStudyDir = EnsureDir(famDir, "EtaSplitStudies");
                TDirectory* catDir      = EnsureDir(etaStudyDir, ObservableDisplayName(catObs));
                TDirectory* catOverlay  = EnsureDir(catDir, "Overlays");
                TDirectory* catSubDir   = EnsureDir(catDir, "Subtracted");

                // The labelled axis carries the selection names, so the ProfileConfig titles stay blank
                ProfileConfig catConfig{catObs.subDir, catObs.profileName, "", catObs.yAxisTitle, false};

                std::vector<ProfileBundle> catBundles;
                std::vector<TH1*> catToDelete;
                catBundles.push_back({pCatData, dataConfig});
                catToDelete.push_back(pCatData);

                for (const auto& sys : sysVariations) {
                    std::string sysFile = consumerDir + "/ConsumerResults_" + fam.dataSuffix + sys.suffix + ".root";
                    TProfile* pSys = cache.FetchClone<TProfile>(sysFile, catPath);
                    if (pSys) { catBundles.push_back({pSys, sys}); catToDelete.push_back(pSys); }
                }
                // Consumer-style external references only. The Toy Model cannot fill a 9- or 15-bin cut scan,
                // so it is handled separately by the three-bin eta-sign summary further below.
                for (const auto& ext : externals) {
                    if (ext.basePath.empty() || ext.kind == ExternalRef::kToyModel) continue;
                    std::string extFile = ext.basePath + "/ConsumerResults_" + fam.dataSuffix + ".root";
                    TProfile* pExt = cache.FetchClone<TProfile>(extFile, catPath);
                    if (pExt) { catBundles.push_back({pExt, ext.config}); catToDelete.push_back(pExt); }
                }

                // Overlay, plus the pull closure panel underneath
                DrawComparisonCanvas(catBundles, "Canvas_AllInOne", fam.familyName + " " + ObservableDisplayName(catObs),
                                     catOverlay, catConfig, "", "", kNoLowerPad, true);
                DrawComparisonCanvas(catBundles, "Canvas_AllInOne_Pull", fam.familyName + " " + ObservableDisplayName(catObs),
                                     catOverlay, catConfig, "", "", kPull, true);

                // Explicit Data - Variation canvases, on the same labelled axis
                std::vector<ProfileBundle> catSubBundles;
                TH1D* pCatZero = SubtractProfiles(pCatData, pCatData, std::string(pCatData->GetName()) + "_catZero");
                for (int i = 1; i <= pCatZero->GetNbinsX(); ++i) pCatZero->SetBinError(i, 0.0);
                catSubBundles.push_back({pCatZero, dataConfig});
                catToDelete.push_back(pCatZero);

                for (size_t i = 1; i < catBundles.size(); ++i) {
                    std::string subName = std::string(pCatData->GetName()) + "_catSub_" + std::to_string(i);
                    TH1D* pSub = SubtractProfiles(pCatData, catBundles[i].profile, subName);
                    catSubBundles.push_back({pSub, catBundles[i].config});
                    catToDelete.push_back(pSub);
                }
                DrawComparisonCanvas(catSubBundles, "Canvas_Subtracted_AllInOne",
                                     fam.familyName + " " + ObservableDisplayName(catObs) + " difference",
                                     catSubDir, catConfig, "", "#Delta" + catObs.yAxisTitle + " (Data - Var)",
                                     kNoLowerPad, true, true);
                DrawComparisonCanvas(catSubBundles, "Canvas_Subtracted_AllInOne_Pull",
                                     fam.familyName + " " + ObservableDisplayName(catObs) + " difference",
                                     catSubDir, catConfig, "", "#Delta" + catObs.yAxisTitle + " (Data - Var)",
                                     kSelfPull, true, true);

                // --- Three-bin eta-sign digest, the one slice the Toy Model can speak to ---
                if (catObs.doEtaSignSummary) {
                    TDirectory* signDir = EnsureDir(catDir, "EtaSignSummary");

                    std::vector<ProfileBundle> signBundles;
                    std::vector<TH1*> signToDelete;
                    for (size_t i = 0; i < catBundles.size(); ++i) {
                        TH1D* h = ExtractEtaSignBins(catBundles[i].profile,
                                                     catObs.profileName + "_sign_" + std::to_string(i));
                        if (!h) continue;
                        signBundles.push_back({h, catBundles[i].config});
                        signToDelete.push_back(h);
                    }

                    // Toy overlay: all three numbers or none, so a partial curve never misleads
                    for (const auto& ext : externals) {
                        if (ext.kind != ExternalRef::kToyModel || ext.basePath.empty()) continue;
                        TH1D* hToy = BuildToyEtaSignHist(cache, ext.basePath, pCatData,
                                                         catObs.profileName + "_signToy");
                        if (hToy) {
                            signBundles.push_back({hToy, ext.config});
                            signToDelete.push_back(hToy);
                        }
                    }

                    if (!signBundles.empty()) {
                        const std::string signTitle =
                            fam.familyName + " " + ObservableDisplayName(catObs) + " eta-sign digest";

                        // Full and redux emitted together, so neither can be forgotten or drift
                        auto emitSign = [&](const std::vector<ProfileBundle>& bs,
                                            const std::string& suffix) {
                            if (bs.size() < 2) return;
                            DrawComparisonCanvas(bs, "Canvas_EtaSign" + suffix, signTitle,
                                                 signDir, catConfig, "", "", kNoLowerPad, true);
                            DrawComparisonCanvas(bs, "Canvas_EtaSign" + suffix + "_Pull", signTitle,
                                                 signDir, catConfig, "", "", kPull, true);
                        };
                        emitSign(signBundles, "");
                        emitSign(ReduxSubset(signBundles), "_Redux");
                    }
                    for (auto p : signToDelete) delete p;
                }

                // --- In-peak vs out-of-peak split of the same scan ---
                // If the fake signal is geometric it must be identical in the signal region and the sidebands,
                // which makes this an independent handle on whether what is being subtracted is really detector
                // response rather than physics. Kept in its own subfolder: useful, but not the headline.
                if (!catObs.massSplitProfile.empty()) {
                    TDirectory* massSplitDir = EnsureDir(catDir, "MassSignalVsBackground");
                    const std::string massPath = TASK_DIR_IN_FILE + catObs.subDir + catObs.massSplitProfile;

                    // Y bin 1 is "out of mass peak", y bin 2 is "in mass peak" (see the consumer's bin labels)
                    const int massYBin[2] = {1, 2};
                    const char* massTag[2] = {"OutOfPeak", "InPeak"};

                    for (int m = 0; m < 2; ++m) {
                        std::vector<ProfileBundle> massBundles;
                        std::vector<TH1*> massToDelete;

                        // Same variation list as above, projected onto the chosen mass slice
                        std::vector<std::pair<std::string, VariationConfig>> sources;
                        sources.push_back({dataFile, dataConfig});
                        for (const auto& sys : sysVariations)
                            sources.push_back({consumerDir + "/ConsumerResults_" + fam.dataSuffix + sys.suffix + ".root", sys});

                        for (const auto& src : sources) {
                            TProfile2D* p2d = cache.FetchClone<TProfile2D>(src.first, massPath, false);
                            if (!p2d) continue;
                            TProfile* proj = p2d->ProfileX((catObs.profileName + "_" + massTag[m] + "_" +
                                                            SanitizeName(src.second.legendLabel)).c_str(),
                                                           massYBin[m], massYBin[m]);
                            if (proj) {
                                proj->SetDirectory(nullptr);
                                // ProfileX builds a fresh x-axis and does NOT carry the parent's bin labels
                                // across, which is why these canvases came out with bare bin numbers instead
                                // of the selection names used everywhere else in this folder. Copy them over.
                                for (int ib = 1; ib <= proj->GetNbinsX(); ++ib) {
                                    const char* lab = p2d->GetXaxis()->GetBinLabel(ib);
                                    if (lab && lab[0] != '\0') proj->GetXaxis()->SetBinLabel(ib, lab);
                                }
                            }
                            delete p2d;
                            if (!proj) continue;
                            massBundles.push_back({proj, src.second});
                            massToDelete.push_back(proj);
                        }

                        if (massBundles.size() > 1) {
                            DrawComparisonCanvas(massBundles, std::string("Canvas_") + massTag[m],
                                                 fam.familyName + " " + ObservableDisplayName(catObs) + " (" + massTag[m] + ")",
                                                 massSplitDir, catConfig, "", "", kNoLowerPad, true);
                            DrawComparisonCanvas(massBundles, std::string("Canvas_") + massTag[m] + "_Pull",
                                                 fam.familyName + " " + ObservableDisplayName(catObs) + " (" + massTag[m] + ")",
                                                 massSplitDir, catConfig, "", "", kPull, true);
                        }
                        for (auto p : massToDelete) delete p;
                    }
                }

                for (auto p : catToDelete) delete p;
            } // end of categorical observables loop
        };

        for (const auto& group : observableGroups) {

            TDirectory* groupDir = nullptr; // Also created lazily

            for (const auto& profConfig : group.profiles) {
                // Species-specific observables only appear on the family they describe. See
                // SpeciesMatchesFamily: a per-species profile read from a single-species file is
                // either empty or a duplicate of that species' own curve.
                if (!SpeciesMatchesFamily(profConfig.speciesTag, fam.familyName)) continue;


                const std::string profilePath = BuildProfilePath(profConfig, cutFolder);

                // Fetch Data:
                std::string dataFile = consumerDir + "/ConsumerResults_" + fam.dataSuffix + ".root";
                TProfile* pData = cache.FetchClone<TProfile>(dataFile, profilePath);
                if (!pData) continue; // Skip to next observable if data is missing (nothing is created on disk)

                // Only now that we know there is something to draw do we materialize the directories
                if (!famDir)   famDir   = EnsureDir(fOut, fam.familyName);
                if (!groupDir) groupDir = EnsureDir(famDir, group.groupName);
                TDirectory* obsDir = EnsureDir(groupDir, ObservableDisplayName(profConfig));

                // Splitting the three canvas categories into their own subfolders, so the observable
                // folder stays navigable instead of holding a flat list of differently-flavoured canvases
                TDirectory* overlayDir = EnsureDir(obsDir, "Overlays");
                TDirectory* foldedDir  = group.doFold ? EnsureDir(obsDir, "Folded") : nullptr;
                TDirectory* subDir     = EnsureDir(obsDir, "Subtracted");

                // Fetch Systematics Data
                std::vector<ProfileBundle> allSystematics;
                allSystematics.push_back({pData, dataConfig}); // Always keep data at index 0 (it is the very first thing in this array)
                
                std::vector<TProfile*> profilesToDelete; // Track for memory cleanup
                profilesToDelete.push_back(pData);

                for (const auto& sys : sysVariations) {
                    std::string sysFile = consumerDir + "/ConsumerResults_" + fam.dataSuffix + sys.suffix + ".root";
                    TProfile* pSys = cache.FetchClone<TProfile>(sysFile, profilePath);
                    if (pSys) {
                        allSystematics.push_back({pSys, sys});
                        profilesToDelete.push_back(pSys);
                    }
                }

                // Fetch the external references (MC base, pp base, Toy Model, ...)
                // presentExternals keeps the original index so the grand summary rows stay aligned
                std::vector<std::pair<size_t, ProfileBundle>> presentExternals;
                for (size_t iExt = 0; iExt < externals.size(); ++iExt) {
                    const auto& ext = externals[iExt];
                    if (ext.basePath.empty()) continue;

                    TProfile* pExt = nullptr;
                    if (ext.kind == ExternalRef::kToyModel) {
                        // Skipped outright when this observable has no toy counterpart. Deciding it from an
                        // explicit name rather than from the consumer-side name means a missing counterpart
                        // produces no curve at all, instead of an unrelated profile that would look
                        // plausible overlaid on a pT, PVz or DeltaPhi axis and be entirely wrong.
                        if (profConfig.toyProfile.empty()) continue;
                        pExt = GetToyModelProfile(cache, ext.basePath, profConfig.toyProfile);
                    } else {
                        std::string extFile = ext.basePath + "/ConsumerResults_" + fam.dataSuffix + ".root";
                        pExt = cache.FetchClone<TProfile>(extFile, profilePath);
                    }

                    if (pExt) {
                        presentExternals.push_back({iExt, {pExt, ext.config}});
                        profilesToDelete.push_back(pExt);
                    }
                }

                // Fetch Mass Selection Data
                std::vector<ProfileBundle> allMassSystematics;
                allMassSystematics.push_back({pData, dataConfig}); // Keep data as baseline
                
                for (const auto& mass : massVariations) {
                    std::string massFile = consumerDir + "/ConsumerResults_" + fam.dataSuffix + mass.suffix + ".root";
                    TProfile* pMass = cache.FetchClone<TProfile>(massFile, profilePath);
                    if (pMass) {
                        allMassSystematics.push_back({pMass, mass});
                        profilesToDelete.push_back(pMass); // Track for standard cleanup
                    }
                }

                // ---------------------------------------------------------
                // FOLDING THE DATA (Summing Pos and Neg Eta)
                // ---------------------------------------------------------
                // Gated on the group flag: folding about zero is only meaningful for symmetric axes
                std::vector<ProfileBundle> allFoldedSystematics;
                std::vector<ProfileBundle> allFoldedExternals;
                std::vector<ProfileBundle> allFoldedMass;
                std::vector<TH1*> foldedToDelete;

                if (group.doFold) {
                    for (const auto& bundle : allSystematics) {
                        // Safely cast to TProfile since we know the original inputs were TProfiles
                        TProfile* pOrig = dynamic_cast<TProfile*>(bundle.profile);
                        std::string foldName = std::string(pOrig->GetName()) + "_Folded_" + bundle.config.suffix;
                        
                        TH1D* pFolded = FoldProfile(pOrig, foldName);
                        allFoldedSystematics.push_back({pFolded, bundle.config});
                        foldedToDelete.push_back(pFolded);
                    }

                    for (const auto& ext : presentExternals) {
                        TProfile* pOrig = dynamic_cast<TProfile*>(ext.second.profile);
                        std::string foldName = std::string(pOrig->GetName()) + "_Folded_Ext" + std::to_string(ext.first);
                        TH1D* pFolded = FoldProfile(pOrig, foldName);
                        allFoldedExternals.push_back({pFolded, ext.second.config});
                        foldedToDelete.push_back(pFolded);
                    }

                    // Fold Mass Selection Data
                    for (const auto& bundle : allMassSystematics) {
                        TProfile* pOrig = dynamic_cast<TProfile*>(bundle.profile);
                        std::string foldName = std::string(pOrig->GetName()) + "_FoldedMass_" + bundle.config.suffix;
                        TH1D* pFolded = FoldProfile(pOrig, foldName);
                        allFoldedMass.push_back({pFolded, bundle.config});
                        foldedToDelete.push_back(pFolded);
                    }
                }

                // Creating the Modulus X-Axis title (e.g. |#eta_{LeadP}|)
                std::string foldedXTitle = "|" + profConfig.xAxisTitle + "|";

                // ---------------------------------------------------------
                // SUBTRACTING THE DATA (Data - Variations)
                // ---------------------------------------------------------
                std::vector<ProfileBundle> allSubtractedSystematics;
                std::vector<ProfileBundle> allSubtractedExternals;
                std::vector<TH1*> subtractedToDelete;

                    // Create a zero-line reference for the Data bundle (could have made this into just a true vertical line,
                // but in this implementation we actually get some points in the plot, which could look better)
                TH1D* pDataZero = SubtractProfiles(pData, pData, std::string(pData->GetName()) + "_Subtracted_Data");
                // Force errors to 0 so the data acts purely as a flat reference line
                for (int i = 1; i <= pDataZero->GetNbinsX(); ++i) pDataZero->SetBinError(i, 0.0);
                
                allSubtractedSystematics.push_back({pDataZero, dataConfig});
                subtractedToDelete.push_back(pDataZero);

                    // Subtract systematics (start loop at 1 to skip the data itself)
                for (size_t i = 1; i < allSystematics.size(); ++i) {
                    const auto& bundle = allSystematics[i];
                    std::string subName = std::string(pData->GetName()) + "_Subtracted_" + bundle.config.suffix;
                    TH1D* pSub = SubtractProfiles(pData, bundle.profile, subName);
                    
                    allSubtractedSystematics.push_back({pSub, bundle.config});
                    subtractedToDelete.push_back(pSub);
                }

                    // Subtract the external references
                for (const auto& ext : presentExternals) {
                    std::string subName = std::string(pData->GetName()) + "_Subtracted_Ext" + std::to_string(ext.first);
                    TH1D* pSub = SubtractProfiles(pData, ext.second.profile, subName);
                    allSubtractedExternals.push_back({pSub, ext.second.config});
                    subtractedToDelete.push_back(pSub);
                }

                    // Subtract Mass Selection Data
                std::vector<ProfileBundle> allSubtractedMass;
                allSubtractedMass.push_back({pDataZero, dataConfig}); // Re-use the existing zero-line reference
                
                for (size_t i = 1; i < allMassSystematics.size(); ++i) {
                    const auto& bundle = allMassSystematics[i];
                    std::string subName = std::string(pData->GetName()) + "_SubtractedMass_" + bundle.config.suffix;
                    TH1D* pSub = SubtractProfiles(pData, bundle.profile, subName);
                    allSubtractedMass.push_back({pSub, bundle.config});
                    subtractedToDelete.push_back(pSub);
                }

                // Creating the Y-Axis title for the difference plots
                std::string subYTitle = "#Delta" + profConfig.yAxisTitle + " (Data - Var)";

                // Draws the plain overlay, then one extra canvas per requested lower-pad mode.
                // Canvases without a data curve to compare against fall back to a single pad inside
                // DrawComparisonCanvas, so this is safe to call on every overlay uniformly.
                // Callers that already emit a full-size standalone difference canvas for the same bundles
                // pass kPullOnly, since a squeezed copy of that difference in a 30%-height pad adds nothing.
                auto drawWithLowerPads = [&](const std::vector<ProfileBundle>& drawBundles,
                                             const std::string& name, const std::string& title,
                                             TDirectory* dir, const std::vector<LowerPadMode>& modes,
                                             const std::string& xT = "", const std::string& yT = "") {
                    DrawComparisonCanvas(drawBundles, name, title, dir, profConfig, xT, yT);
                    if (drawBundles.size() < 2) return; // Nothing to compare: the variants would duplicate the overlay
                    for (auto mode : modes) {
                        DrawComparisonCanvas(drawBundles, name + LowerPadSuffix(mode), title, dir,
                                             profConfig, xT, yT, mode);
                    }
                };

                // ---------------------------------------------------------
                // DRAWING THE VARIATIONS
                // ---------------------------------------------------------
                
                // --- 1. Standard (Unfolded) Plots ---
                    // Variation 1: Data Only
                DrawComparisonCanvas({{pData, dataConfig}}, "Canvas_DataOnly", fam.familyName + " Data", overlayDir, profConfig);
                    // Variation 2: Systematics (+ Data)
                drawWithLowerPads(allSystematics, "Canvas_Systematics", fam.familyName + " Systematics", overlayDir, kPullOnly);
                    // Variation 3+: one canvas per external reference, then the all-in-one
                std::vector<ProfileBundle> allInOne = allSystematics;
                for (const auto& ext : presentExternals) {
                    const std::string& label = ext.second.config.legendLabel;
                    DrawComparisonCanvas({ext.second}, "Canvas_Only_" + SanitizeName(externals[ext.first].config.legendLabel),
                                         fam.familyName + " " + label, overlayDir, profConfig);
                    allInOne.push_back(ext.second);
                }
                if (!presentExternals.empty()) {
                    drawWithLowerPads(allInOne, "Canvas_AllInOne", fam.familyName + " All Comparisons", overlayDir, kPullOnly);
                }

                    // Variation 2b: the reduced canvas. Same content, minus the variations whose
                    // interest is a comparison against a specific partner rather than against the
                    // data -- the gated proxies and the data-like one. Those still have their own
                    // canvases below, so nothing is lost, and this one stays legible.
                    // README: the inRedux / inPtGateStudy / inEtaGateStudy flags on VariationConfig
                    // decide membership; adding a variation is a one-line table edit.
                std::vector<ProfileBundle> reduxBundles = ReduxSubset(allSystematics);
                if (reduxBundles.size() > 1) {
                    drawWithLowerPads(reduxBundles, "Canvas_Systematics_Redux",
                                      fam.familyName + " Systematics (reduced)", overlayDir, kPullOnly);
                }

                    // The two gate studies, each on its own canvas, in their own folder. Data is kept
                    // as the reference in both.
                    // README: the pT study asks whether ADDING a gate changes an artificial proxy;
                    // the eta study asks what is lost by REMOVING one, since eta gating is now the
                    // default. Opposite polarity, same plot shape.
                TDirectory* gateDir = nullptr;
                auto drawGateStudy = [&](bool VariationConfig::*flag, const std::string& name,
                                         const std::string& title) {
                    std::vector<ProfileBundle> gb;
                    gb.push_back(allSystematics[0]); // Data heads the legend and drives the lower pad
                    for (size_t i = 1; i < allSystematics.size(); ++i)
                        if (allSystematics[i].config.*flag) gb.push_back(allSystematics[i]);
                    if (gb.size() < 2) return;
                    if (!gateDir) gateDir = EnsureDir(obsDir, "GateStudies");
                    drawWithLowerPads(gb, name, fam.familyName + " " + title, gateDir, kPullOnly);
                };
                drawGateStudy(&VariationConfig::inPtGateStudy, "Canvas_PtGateStudy",
                              "pT gate on artificial proxies");
                drawGateStudy(&VariationConfig::inEtaGateStudy, "Canvas_EtaGateStudy",
                              "#eta acceptance gate on artificial proxies");

                // --- 2. Folded Plots ---
                if (group.doFold) {
                    DrawComparisonCanvas({allFoldedSystematics[0]}, "Canvas_Folded_DataOnly", fam.familyName + " (R(#eta_{pos})+R(#eta_{neg})) Data", foldedDir, profConfig, foldedXTitle);
                    drawWithLowerPads(allFoldedSystematics, "Canvas_Folded_Systematics", fam.familyName + " (R(#eta_{pos})+R(#eta_{neg})) Systematics", foldedDir, lowerPadModes, foldedXTitle);

                    std::vector<ProfileBundle> allInOneFolded = allFoldedSystematics; // For folded plots tracking
                    for (size_t i = 0; i < allFoldedExternals.size(); ++i) {
                        const std::string& label = allFoldedExternals[i].config.legendLabel;
                        DrawComparisonCanvas({allFoldedExternals[i]}, "Canvas_Folded_Only_" + SanitizeName(label),
                                             fam.familyName + " (R(#eta_{pos})+R(#eta_{neg})) " + label, foldedDir, profConfig, foldedXTitle);
                        allInOneFolded.push_back(allFoldedExternals[i]);
                    }
                    if (!allFoldedExternals.empty()) {
                        drawWithLowerPads(allInOneFolded, "Canvas_Folded_AllInOne", fam.familyName + " (R(#eta_{pos})+R(#eta_{neg})) Comparisons", foldedDir, lowerPadModes, foldedXTitle);
                    }
                }

                // --- 3. Subtracted Data - Systematics ---
                    // Data minus Data Systematics
                DrawComparisonCanvas(allSubtractedSystematics, "Canvas_Subtracted_Systematics", fam.familyName + " Systematics Difference", subDir, profConfig, "", subYTitle, kNoLowerPad, false, true);
                DrawComparisonCanvas(allSubtractedSystematics, "Canvas_Subtracted_Systematics_Pull", fam.familyName + " Systematics Difference", subDir, profConfig, "", subYTitle, kSelfPull, false, true);
                    // Data minus each external reference
                std::vector<ProfileBundle> allInOneSubtracted = allSubtractedSystematics;
                for (size_t i = 0; i < allSubtractedExternals.size(); ++i) {
                    const std::string& label = allSubtractedExternals[i].config.legendLabel;
                    DrawComparisonCanvas({allSubtractedExternals[i]}, "Canvas_Subtracted_Only_" + SanitizeName(label),
                                         fam.familyName + " " + label + " Difference", subDir, profConfig, "", subYTitle, kNoLowerPad, false, true);
                    allInOneSubtracted.push_back(allSubtractedExternals[i]);
                }
                    // All-In-One
                if (!allSubtractedExternals.empty()) {
                    DrawComparisonCanvas(allInOneSubtracted, "Canvas_Subtracted_AllInOne", fam.familyName + " All Comparisons Difference", subDir, profConfig, "", subYTitle, kNoLowerPad, false, true);
                    DrawComparisonCanvas(allInOneSubtracted, "Canvas_Subtracted_AllInOne_Pull", fam.familyName + " All Comparisons Difference", subDir, profConfig, "", subYTitle, kSelfPull, false, true);
                }

                // --- Mass Selection Proxy Plots ---
                TDirectory* massDir = nullptr;
                if (allMassSystematics.size() > 1) { // Only draw if mass files were found
                    massDir = EnsureDir(obsDir, "Mass_Selection");
                    
                    // Standard
                    DrawComparisonCanvas(allMassSystematics, "Canvas_MassSelection", fam.familyName + " Mass Selection", massDir, profConfig);
                    // Folded
                    if (group.doFold) {
                        DrawComparisonCanvas(allFoldedMass, "Canvas_Folded_MassSelection", fam.familyName + " (R(#eta_{pos})+R(#eta_{neg})) Mass Selection", massDir, profConfig, foldedXTitle);
                    }
                    // Subtracted
                    DrawComparisonCanvas(allSubtractedMass, "Canvas_Subtracted_MassSelection", fam.familyName + " Mass Selection Difference", massDir, profConfig, "", subYTitle, kNoLowerPad, false, true);
                    DrawComparisonCanvas(allSubtractedMass, "Canvas_Subtracted_MassSelection_Pull", fam.familyName + " Mass Selection Difference", massDir, profConfig, "", subYTitle, kSelfPull, false, true);
                }

                // --- 4. Individual Comparisons (One-by-One) ---
                // (condensed all of this into a single block of code because the loops become cleaner!)
                // Gated per observable: this block alone is roughly 20 canvases per observable, and the
                // overlays plus the difference canvases already cover the non-LeadP proxies adequately.
                if (doIndividualComparisons && profConfig.doIndividualComparisons) {
                    // Create a sub-directory specifically for the individual comparisons to keep things organized
                    TDirectory* indivDir = EnsureDir(obsDir, "Individual_Comparisons");

                    // Loop through all variations (starting at i = 1 to skip the data comparing against itself)
                    for (size_t i = 1; i < allSystematics.size(); ++i) {
                        std::string sysLabel = allSystematics[i].config.legendLabel;
                        std::string sysSuffix = allSystematics[i].config.suffix; // e.g., "_forceDatalikeJet"

                        // Standard (now also carrying the difference and pull panels underneath)
                        drawWithLowerPads({allSystematics[0], allSystematics[i]}, "Canvas_Standard" + sysSuffix, fam.familyName + " Data vs " + sysLabel, indivDir, kPullOnly);
                        // Folded
                        if (group.doFold) {
                            drawWithLowerPads({allFoldedSystematics[0], allFoldedSystematics[i]}, "Canvas_Folded" + sysSuffix, fam.familyName + " Folded vs " + sysLabel, indivDir, lowerPadModes, foldedXTitle);
                        }
                        // Subtracted
                        DrawComparisonCanvas({allSubtractedSystematics[i]}, "Canvas_Subtracted" + sysSuffix, fam.familyName + " Difference vs " + sysLabel, indivDir, profConfig, "", subYTitle, kNoLowerPad, false, true);
                    }

                    // Do the same for every external reference, which are introduced separately
                    for (size_t i = 0; i < presentExternals.size(); ++i) {
                        const std::string& label = presentExternals[i].second.config.legendLabel;

                        drawWithLowerPads({allSystematics[0], presentExternals[i].second}, "Canvas_Standard_" + SanitizeName(label), fam.familyName + " Data vs " + label, indivDir, kPullOnly);
                        if (group.doFold) {
                            drawWithLowerPads({allFoldedSystematics[0], allFoldedExternals[i]}, "Canvas_Folded_" + SanitizeName(label), fam.familyName + " Folded vs " + label, indivDir, lowerPadModes, foldedXTitle);
                        }
                        DrawComparisonCanvas({allSubtractedExternals[i]}, "Canvas_Subtracted_" + SanitizeName(label), fam.familyName + " Difference vs " + label, indivDir, profConfig, "", subYTitle, kNoLowerPad, false, true);
                    }
                }

                // --- 5. Integrated accumulation (drawn once per family, not once per observable) ---
                // Integrating a profile just recovers the entry-weighted mean of whatever filled it, so
                // two observables filled from the SAME entries with different x-binnings give the SAME
                // number. Every LeadP observable therefore produced an identical Integrated_Summary, as
                // did every LeadJet one, and so on: 20 of the 28 folders per family were exact duplicates
                // of one of the other 8. Nothing is drawn here any more; the values are accumulated and
                // presented once, at family level, by the brute-force summaries below.
                auto dataInteg = GetIntegratedProfile(pData);

                // Observables with no proxy tag stay off the brute-force axis. That is the eta-sign pT
                // splits: they cover different subsets of the sample, so their integrals sit at a
                // different scale and would flatten everything else on a shared axis.
                if (!profConfig.proxyTag.empty()) {
                    grandProxyTags.push_back(profConfig.proxyTag);
                    grandDataVals.push_back(dataInteg);
                    grandLabels.push_back(ObservableDisplayName(profConfig));

                    // (allSystematics[0] is the data, so sysVariations[i] corresponds to allSystematics[i+1])
                    for (size_t i = 0; i < sysVariations.size(); ++i) {
                        if (i + 1 < allSystematics.size()) {
                            auto sysInteg = GetIntegratedProfile(dynamic_cast<TProfile*>(allSystematics[i+1].profile));
                            grandSysVals[i].push_back(sysInteg);
                        }
                    }
                    for (const auto& ext : presentExternals) {
                        grandExtVals[ext.first].push_back(GetIntegratedProfile(dynamic_cast<TProfile*>(ext.second.profile)));
                    }
                    for (size_t i = 0; i < massVariations.size(); ++i) {
                        if (i + 1 < allMassSystematics.size()) {
                            grandMassVals[i].push_back(GetIntegratedProfile(dynamic_cast<TProfile*>(allMassSystematics[i+1].profile)));
                        }
                    }
                }

                // NOTE: the per-observable IntegratedSummary block used to live here. It re-integrated
                // each differential profile with GetIntegratedProfile, which recomputed a number the
                // consumer already publishes in bin 1 of IntegratedCuts/pRingCuts*. That is now read
                // directly, once per family, in <Family>/IntegratedSummary/<proxy>/.
                // README: verified against the old values before removal. The only disagreements were
                // the two known overflow cases (proxy pT, EtaV0) plus a small one on the centrality
                // axis, still unexplained; in all three the consumer value is the one to report.

                // Cross-family comparison: only the three representative observables are recorded, one per
                // proxy, since anything else would just repeat one of them (same argument as above).
                for (size_t iRep = 0; iRep < proxyRepresentatives.size(); ++iRep) {
                    if (profConfig.profileName != proxyRepresentatives[iRep].second) continue;
                    auto& slot = crossFamilyVals[iRep][famIdx];
                    slot.push_back(dataInteg);
                    for (size_t i = 1; i < allSystematics.size(); ++i) {
                        slot.push_back(GetIntegratedProfile(dynamic_cast<TProfile*>(allSystematics[i].profile)));
                    }
                    if (crossFamilyLabels[iRep].empty()) {
                        for (const auto& b2 : allSystematics) crossFamilyLabels[iRep].push_back(b2.config.legendLabel);
                    }
                }

                // --- 6. tanh fit of the eta dependence ---
                // Only for the fine-binned eta profiles, where there are enough points across the
                // turn-over for p0 and p1 to separate at all.
                if (profConfig.doTanhFit) {
                    TDirectory* fitDir    = EnsureDir(obsDir, "TanhFit");
                    TDirectory* perVarDir = EnsureDir(fitDir, "PerVariation");

                    // Fit every variation, then every external reference, in the same axis order used
                    // by the integrated summaries so the four canvases below line up with them.
                    std::vector<std::string> fitLabels;
                    std::vector<TanhFitResult> fitResults;
                    std::vector<VariationConfig> fitConfigs;
                    std::vector<bool> fitIsToy; // Toy is excluded from the entries summary, see below

                    for (const auto& bundle : allSystematics) {
                        TProfile* pFit = dynamic_cast<TProfile*>(bundle.profile);
                        if (!pFit) continue;
                        fitLabels.push_back(bundle.config.legendLabel);
                        fitConfigs.push_back(bundle.config);
                        fitIsToy.push_back(false);
                        fitResults.push_back(FitTanhProfile(pFit, profConfig, bundle.config,
                                                            "Canvas_TanhFit_" + SanitizeName(bundle.config.legendLabel),
                                                            perVarDir));
                    }
                    for (const auto& ext : presentExternals) {
                        TProfile* pFit = dynamic_cast<TProfile*>(ext.second.profile);
                        if (!pFit) continue;
                        fitLabels.push_back(ext.second.config.legendLabel);
                        fitConfigs.push_back(ext.second.config);
                        fitIsToy.push_back(externals[ext.first].kind == ExternalRef::kToyModel);
                        fitResults.push_back(FitTanhProfile(pFit, profConfig, ext.second.config,
                                                            "Canvas_TanhFit_" + SanitizeName(ext.second.config.legendLabel),
                                                            perVarDir));
                    }

                    // Four cross-variation summaries. Entries and reduced chi2 keep EVERY variation,
                    // including the ones whose fit did not converge: a wild chi2/ndf is itself the
                    // quickest way to spot a failed fit. The parameter summaries keep only the
                    // converged ones, since a non-converged p0 or p1 is not a measurement of anything.
                    std::vector<std::string> allLabels, okLabels, entriesLabels;
                    std::vector<FitSummaryPoint> entriesPts, chi2Pts, p0Pts, p1Pts;

                    for (size_t i = 0; i < fitResults.size(); ++i) {
                        const auto& fr = fitResults[i];

                        allLabels.push_back(fitLabels[i]);
                        chi2Pts.push_back({fitConfigs[i], (fr.ndf > 0) ? fr.chi2 / fr.ndf : 0.0, 0.0, 0.0});

                        // The Toy Model is run with something like three orders of magnitude more
                        // statistics than any data or MC sample, so including it forces a y-range in
                        // which every real point collapses onto the axis. It is left out of the entries
                        // summary entirely; its statistics are a property of how long the toy was run
                        // for, not a result worth comparing.
                        if (!fitIsToy[i]) {
                            entriesLabels.push_back(fitLabels[i]);
                            entriesPts.push_back({fitConfigs[i], fr.entries, 0.0, 0.0});
                        }

                        if (!fr.valid) continue;
                        okLabels.push_back(fitLabels[i]);
                        p0Pts.push_back({fitConfigs[i], fr.p0, fr.p0Low, fr.p0High});
                        p1Pts.push_back({fitConfigs[i], fr.p1, fr.p1Low, fr.p1High});
                    }

                    DrawFitSummaryCanvas(entriesLabels, entriesPts, "Canvas_FitSummary_Entries",
                                         fam.familyName + " tanh fit: entries (Toy excluded)", fitDir, "Entries");
                    DrawFitSummaryCanvas(allLabels, chi2Pts, "Canvas_FitSummary_RedChi2",
                                         fam.familyName + " tanh fit: reduced #chi^{2}", fitDir, "#chi^{2}/ndf");
                    DrawFitSummaryCanvas(okLabels, p0Pts, "Canvas_FitSummary_p0",
                                         fam.familyName + " tanh fit: p0 (amplitude)", fitDir, "p0");
                    DrawFitSummaryCanvas(okLabels, p1Pts, "Canvas_FitSummary_p1",
                                         fam.familyName + " tanh fit: p1 (slope scale)", fitDir, "p1");
                }

                // Cleanup dynamically allocated profiles for this observable iteration
                for (auto p : profilesToDelete) delete p;
                for (auto p : foldedToDelete) delete p;
                for (auto p : subtractedToDelete) delete p;
            } // end of profiles loop, inside groups loop

            // Emit the eta-split cut scans immediately after EtaProxy, so they land second in the file
            if (group.groupName == "EtaProxy") runEtaSplitStudies();
            // Same reasoning: HEE belongs immediately after AEE, not at the end of the family
            if (group.groupName == "AEE") runAeeHeeBlocks();
        } // end of groups loop, inside families loop

        // Safety net: if EtaProxy produced nothing at all (its profiles absent from these files), the
        // scans have still not been written, so do it here rather than losing them entirely.
        runEtaSplitStudies();
        runAeeHeeBlocks(); // Safety net, as above


        // --- 6. Integrated summaries, drawn ONCE per family ---
        // Everything integrated now lives here rather than being repeated inside every observable folder.
        // The x-axis is the observable, the series are the variations, so this single canvas carries all
        // the information the 28 per-observable folders used to hold between them.
        // Comparing the eta_Jet and eta_Lambda entries is also a way of probing for overflow/underflow
        // problems (there were some! The eta axis range in the consumer was fixed after this was spotted).
        int nGrand = grandLabels.size();

        // -----------------------------------------------------------------------------------
        // Family-level integrated summary, from the consumer's own integrated profiles
        // -----------------------------------------------------------------------------------
        // One place per family where "the integrated value" lives, instead of one copy inside every
        // observable folder. Nothing is integrated here: each point is read from bin 1 of the
        // matching IntegratedCuts profile, so every variation reports the number its own consumer
        // run computed.
        // README: this supersedes the per-observable IntegratedSummary folders. Those recompute the
        // same quantity from differential profiles and are kept only until the cross-check below has
        // been looked at; see the TODO on GetIntegratedProfile.
        {
            const std::string dataFileInt = consumerDir + "/ConsumerResults_" + fam.dataSuffix + ".root";
            TDirectory* famIntDir = nullptr;

            for (const auto& proxy : sigExtractProxies) {
                std::vector<std::string> labels;
                std::vector<VariationConfig> cfgs;
                std::vector<std::pair<double, double>> vals;

                SigExtractPoint dataPt = FetchIntegratedFromCuts(cache, dataFileInt, proxy.folderName);
                if (!dataPt.found) continue;
                labels.push_back(dataConfig.legendLabel);
                cfgs.push_back(dataConfig);
                vals.push_back({dataPt.value, dataPt.error});

                for (const auto& sys : sysVariations) {
                    const std::string f = consumerDir + "/ConsumerResults_" + fam.dataSuffix + sys.suffix + ".root";
                    SigExtractPoint p = FetchIntegratedFromCuts(cache, f, proxy.folderName);
                    if (!p.found) continue;
                    labels.push_back(sys.legendLabel);
                    cfgs.push_back(sys);
                    vals.push_back({p.value, p.error});
                }

                // External references, in the same order they appear on every other categorical axis.
                //
                // MC and pp are consumer output, so they carry the same IntegratedCuts profiles and
                // are read exactly like the variations above.
                //
                // The Toy Model has its own integrated object, TOY_SCENARIO_DIR/pRingProxyJet, a
                // single-bin TProfile. Nothing is integrated here either: bin 1 is already the
                // answer. Note pRingProxyJet, NOT pRingProxy -- the latter uses z-hat as the axis
                // rather than a jet direction, and would be a different observable entirely.
                // The same toy value is placed on all three proxy axes, which is the convention the
                // differential plots already follow: the toy generates one random jet direction and
                // one hyperon species, so it is a reference for the geometric effect rather than a
                // per-proxy or per-species prediction.
                for (const auto& ext : externals) {
                    if (ext.basePath.empty()) continue;

                    SigExtractPoint p;
                    if (ext.kind == ExternalRef::kToyModel) {
                        TProfile* pToy = GetToyModelProfile(cache, ext.basePath, "pRingProxyJet");
                        if (pToy) {
                            if (pToy->GetNbinsX() >= 1 && pToy->GetBinEntries(1) > 0) {
                                p.value = pToy->GetBinContent(1);
                                p.error = pToy->GetBinError(1);
                                p.found = true;
                            }
                            delete pToy;
                        }
                    } else {
                        const std::string extFile = ext.basePath + "/ConsumerResults_" + fam.dataSuffix + ".root";
                        p = FetchIntegratedFromCuts(cache, extFile, proxy.folderName);
                    }

                    if (!p.found) continue;
                    labels.push_back(ext.config.legendLabel);
                    cfgs.push_back(ext.config);
                    vals.push_back({p.value, p.error});
                }

                // The signal-extracted value, as an extra column. Its error carries the full
                // propagation done in signalExtractionRing, so it is fetched, never recomputed.
                if (!sigExtractDir.empty()) {
                    SigExtractPoint sp = FetchSigExtractPoint(cache, sigExtractDir, fam.dataSuffix,
                                                              cutFolder, proxy.folderName,
                                                              "hIntegratedRSig");
                    if (sp.found) {
                        VariationConfig c;
                        c.legendLabel = "Sig. extracted";
                        c.color = kMagenta + 2;
                        c.lineStyle = 1;
                        c.markerStyle = 29;
                        c.isData = false;
                        c.inRedux = true;
                        labels.push_back(c.legendLabel);
                        cfgs.push_back(c);
                        vals.push_back({sp.value, sp.error});
                    }
                }

                if (labels.size() < 2) continue;
                if (!famIntDir) famIntDir = EnsureDir(EnsureDir(fOut, fam.familyName), "IntegratedSummary");
                TDirectory* pxDir = EnsureDir(famIntDir, proxy.folderName);

                // Full and redux emitted from the same call site, as everywhere else
                auto emitAxis = [&](std::function<bool(const VariationConfig&)> keep,
                                    const std::string& name, const std::string& title,
                                    TDirectory* dir) {
                    std::vector<int> pick;
                    for (size_t i = 0; i < cfgs.size(); ++i)
                        if (i == 0 || keep(cfgs[i])) pick.push_back(static_cast<int>(i));
                    if (pick.size() < 2) return;

                    const int n = static_cast<int>(pick.size());
                    std::vector<std::string> subLabels;
                    std::vector<ProfileBundle> bVal, bSub;
                    std::vector<TH1*> del;
                    for (int k = 0; k < n; ++k) {
                        const int i = pick[k];
                        subLabels.push_back(labels[i]);
                        TH1D* hV = MakeCategoricalPoint(Form("%s_%s_V_%d", name.c_str(),
                                                             proxy.folderName.c_str(), k),
                                                        k, n, vals[i].first, vals[i].second);
                        bVal.push_back({hV, cfgs[i]}); del.push_back(hV);
                        if (k == 0) continue;
                        const double d = vals[0].first - vals[i].first;
                        const double e = std::sqrt(vals[0].second * vals[0].second +
                                                   vals[i].second * vals[i].second);
                        TH1D* hS = MakeCategoricalPoint(Form("%s_%s_S_%d", name.c_str(),
                                                             proxy.folderName.c_str(), k),
                                                        k, n, d, e);
                        bSub.push_back({hS, cfgs[i]}); del.push_back(hS);
                    }
                    DrawIntegratedCanvas(bVal, subLabels, name,
                                         fam.familyName + " " + proxy.legendLabel + " " + title,
                                         dir, "Integrated <R>", false, false, true);
                    if (!bSub.empty())
                        DrawIntegratedCanvas(bSub, subLabels, name + "_Subtracted",
                                             fam.familyName + " " + proxy.legendLabel + " " + title +
                                                 ", difference",
                                             dir, "#Delta<R> (Data - Var)", true, false, true);
                    for (auto p : del) delete p;
                };

                emitAxis([](const VariationConfig&) { return true; },
                         "Canvas_Integrated", "integrated <R>", pxDir);
                emitAxis([](const VariationConfig& c) { return c.inRedux; },
                         "Canvas_Integrated_Redux", "integrated <R> (reduced)", pxDir);

                TDirectory* gDir = EnsureDir(pxDir, "GateStudies");
                emitAxis([](const VariationConfig& c) { return c.inPtGateStudy; },
                         "Canvas_Integrated_PtGate", "integrated, pT gate", gDir);
                emitAxis([](const VariationConfig& c) { return c.inEtaGateStudy; },
                         "Canvas_Integrated_EtaGate", "integrated, #eta acceptance gate", gDir);
            }
        }

        // -----------------------------------------------------------------------------------
        // CheapSigExtract: the same integrated axis, in and out of the mass peak
        // -----------------------------------------------------------------------------------
        // Two canvases of the same shape as Canvas_Integrated_Redux, one built from candidates
        // strictly inside the mass window and one from candidates strictly outside it. The gap
        // between them is how much the background moves the integrated number -- measured without
        // fitting anything, so it cannot be wrong in the way a sideband model can.
        {
            const std::string dataFileCheap = consumerDir + "/ConsumerResults_" + fam.dataSuffix + ".root";
            TDirectory* cheapDir = nullptr;

            for (const auto& proxy : sigExtractProxies) {
                for (int pass = 0; pass < 2; ++pass) {
                    const bool inPeak = (pass == 0);
                    std::vector<std::string> labels;
                    std::vector<VariationConfig> cfgs;
                    std::vector<std::pair<double, double>> vals;

                    SigExtractPoint d = FetchIntegratedByMassPeak(cache, dataFileCheap,
                                                                  proxy.folderName, inPeak);
                    if (!d.found) continue;
                    labels.push_back(dataConfig.legendLabel);
                    cfgs.push_back(dataConfig);
                    vals.push_back({d.value, d.error});

                    for (const auto& sys : sysVariations) {
                        if (!sys.inRedux) continue; // Redux membership, same flag as everywhere else
                        const std::string f = consumerDir + "/ConsumerResults_" + fam.dataSuffix +
                                              sys.suffix + ".root";
                        SigExtractPoint p = FetchIntegratedByMassPeak(cache, f, proxy.folderName, inPeak);
                        if (!p.found) continue;
                        labels.push_back(sys.legendLabel);
                        cfgs.push_back(sys);
                        vals.push_back({p.value, p.error});
                    }
                    if (labels.size() < 2) continue;

                    if (!cheapDir) cheapDir = EnsureDir(EnsureDir(fOut, fam.familyName), "CheapSigExtract");
                    TDirectory* pxDir = EnsureDir(cheapDir, proxy.folderName);

                    const int n = static_cast<int>(labels.size());
                    const std::string tag = inPeak ? "InMassPeak" : "OutOfMassPeak";
                    const std::string ttl = fam.familyName + " " + proxy.legendLabel + ", " +
                                            (inPeak ? "in mass peak" : "out of mass peak");

                    std::vector<ProfileBundle> bVal, bSub;
                    std::vector<TH1*> del;
                    for (int k = 0; k < n; ++k) {
                        TH1D* hV = MakeCategoricalPoint(Form("Cheap_%s_%s_V_%d", tag.c_str(),
                                                             proxy.folderName.c_str(), k),
                                                        k, n, vals[k].first, vals[k].second);
                        bVal.push_back({hV, cfgs[k]}); del.push_back(hV);
                        if (k == 0) continue;
                        const double diff = vals[0].first - vals[k].first;
                        const double ediff = std::sqrt(vals[0].second * vals[0].second +
                                                       vals[k].second * vals[k].second);
                        TH1D* hS = MakeCategoricalPoint(Form("Cheap_%s_%s_S_%d", tag.c_str(),
                                                             proxy.folderName.c_str(), k),
                                                        k, n, diff, ediff);
                        bSub.push_back({hS, cfgs[k]}); del.push_back(hS);
                    }

                    DrawIntegratedCanvas(bVal, labels, "Canvas_Integrated_Redux_" + tag, ttl,
                                         pxDir, "Integrated <R>", false, false, true);
                    if (!bSub.empty())
                        DrawIntegratedCanvas(bSub, labels, "Canvas_Integrated_Redux_" + tag + "_Subtracted",
                                             ttl + ", difference", pxDir,
                                             "#Delta<R> (Data - Var)", true, false, true);
                    for (auto p : del) delete p;
                }
            }
        }

        // -----------------------------------------------------------------------------------
        // Signal-extracted integrated observables
        // -----------------------------------------------------------------------------------
        // Fed by signalExtractionRing rather than by the consumer, so this whole block is skipped
        // when no extraction directory was given. Directories are created lazily, exactly as
        // everywhere else here, so a run without extraction output leaves no empty folders.
        if (!sigExtractDir.empty()) {
            TDirectory* sigExtractGroupDir = nullptr;

            for (size_t iProxy = 0; iProxy < sigExtractProxies.size(); ++iProxy) {
                const auto& proxy = sigExtractProxies[iProxy];

                // Is there anything at all for this proxy? Probing one histogram avoids creating a
                // folder for a proxy whose extraction never ran or did not converge.
                SigExtractPoint probe = FetchSigExtractPoint(cache, sigExtractDir, fam.dataSuffix,
                                                             cutFolder, proxy.folderName,
                                                             "hIntegratedRSig");
                if (!probe.found) continue;

                if (!famDir) famDir = EnsureDir(fOut, fam.familyName);
                if (!sigExtractGroupDir) sigExtractGroupDir = EnsureDir(famDir, "SignalExtraction");
                TDirectory* proxyDir = EnsureDir(sigExtractGroupDir, proxy.folderName);

                // (a) The three components side by side, plus their correlated differences
                DrawSigExtractComponents(cache, sigExtractDir, fam.dataSuffix, cutFolder,
                                         proxy, fam.familyName,
                                         EnsureDir(proxyDir, "Components"));

                // (b) One systematics canvas pair per quantity, imitating Integrated_summary
                for (size_t iQ = 0; iQ < sigExtractQuantities.size(); ++iQ) {
                    const auto& quantity = sigExtractQuantities[iQ];
                    SigExtractPoint dataVal = DrawSigExtractSystematics(
                        cache, sigExtractDir, fam.dataSuffix, sysVariations, dataConfig,
                        cutFolder, proxy, quantity, fam.familyName,
                        EnsureDir(proxyDir, SanitizeName(quantity.legendLabel)));
                    sigExtractCrossFamily[iQ][iProxy][famIdx] = dataVal;
                }

                // (c) Signal and background as two series over the systematics axis, so the two are
                //     read against each other directly instead of on facing canvases. This is the
                //     comparison that matters: an artificial-proxy variation that moves <R>_S and
                //     <R>_B together is saying something quite different from one that moves only
                //     the signal.
                std::vector<std::string> sysLabels{dataConfig.legendLabel};
                std::vector<std::vector<SigExtractPoint>> seriesPts(sigExtractQuantities.size());
                for (size_t iQ = 0; iQ < sigExtractQuantities.size(); ++iQ) {
                    seriesPts[iQ].push_back(FetchSigExtractPoint(cache, sigExtractDir,
                                                                 fam.dataSuffix, cutFolder,
                                                                 proxy.folderName,
                                                                 sigExtractQuantities[iQ].histName));
                }
                for (const auto& sys : sysVariations) {
                    std::vector<SigExtractPoint> row;
                    bool complete = true;
                    for (const auto& quantity : sigExtractQuantities) {
                        SigExtractPoint p = FetchSigExtractPoint(cache, sigExtractDir,
                                                                 fam.dataSuffix + sys.suffix,
                                                                 cutFolder, proxy.folderName,
                                                                 quantity.histName);
                        if (!p.found) { complete = false; break; }
                        row.push_back(p);
                    }
                    if (!complete) continue; // Keep the two series on a common axis
                    sysLabels.push_back(sys.legendLabel);
                    for (size_t iQ = 0; iQ < row.size(); ++iQ) seriesPts[iQ].push_back(row[iQ]);
                }

                if (sysLabels.size() >= 2) {
                    const int nSys = static_cast<int>(sysLabels.size());
                    std::vector<ProfileBundle> series;
                    std::vector<TH1*> seriesToDelete;
                    for (size_t iQ = 0; iQ < sigExtractQuantities.size(); ++iQ) {
                        TH1D* h = new TH1D(Form("SigSeries_%s_%zu", proxy.folderName.c_str(), iQ),
                                           "", nSys, 0, nSys);
                        h->SetDirectory(nullptr);
                        for (int i = 0; i < nSys; ++i) {
                            h->SetBinContent(i + 1, seriesPts[iQ][i].value);
                            h->SetBinError(i + 1, seriesPts[iQ][i].error);
                        }
                        VariationConfig vc{"", sigExtractQuantities[iQ].legendLabel,
                                           sigExtractQuantities[iQ].color, 1,
                                           sigExtractQuantities[iQ].markerStyle, false};
                        series.push_back({h, vc});
                        seriesToDelete.push_back(h);
                    }
                    DrawOffsetSeriesCanvas(sysLabels, series, "Canvas_SignalVsBackground",
                                           fam.familyName + " " + proxy.legendLabel +
                                               ": <R>_{S} and <R>_{B} across systematics",
                                           proxyDir, "Integrated <R>", 0.5, true);
                    for (auto q : seriesToDelete) delete q;
                }
            }
        }

        if (nGrand > 0 && famDir) {
            TDirectory* integDir = EnsureDir(famDir, "BruteForce");

            // Reorder the axis so observables sharing a proxy sit together (all LeadP, then all LeadJet,
            // then all SubJet) instead of being interleaved by whichever group happened to run first.
            // Observables sharing a proxy share an entry set, so their integrals are identical, and having
            // them adjacent turns each block into a flat plateau that is trivial to read.
            std::vector<size_t> order(grandLabels.size());
            for (size_t i = 0; i < order.size(); ++i) order[i] = i;
            const std::vector<std::string> proxyOrder = {"LeadP", "LeadJet", "SubJet"};
            auto proxyRank = [&](const std::string& tag) {
                for (size_t r = 0; r < proxyOrder.size(); ++r) if (proxyOrder[r] == tag) return r;
                return proxyOrder.size(); // Anything unrecognised sorts last rather than being dropped
            };
            std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
                return proxyRank(grandProxyTags[a]) < proxyRank(grandProxyTags[b]);
            });

            // Rebuild the labels in the new order, prefixed with the proxy so the blocks are self-evident
            std::vector<std::string> sortedLabels;
            for (size_t k : order) sortedLabels.push_back(grandProxyTags[k] + ": " + grandLabels[k]);

            auto reorder = [&](const std::vector<std::pair<double, double>>& v) {
                std::vector<std::pair<double, double>> out;
                if (v.size() != order.size()) return out; // Incomplete series: skipped by makeSeries below
                for (size_t k : order) out.push_back(v[k]);
                return out;
            };

            std::vector<TH1*> grandToDelete;

            // Small local helper: one multi-bin series across all observables, for a given variation
            auto makeSeries = [&](const std::vector<std::pair<double, double>>& vals,
                                  const std::string& name) -> TH1D* {
                if (vals.size() != static_cast<size_t>(nGrand)) return nullptr;
                TH1D* h = new TH1D(name.c_str(), "", nGrand, 0, nGrand);
                h->SetDirectory(nullptr);
                for (int i = 0; i < nGrand; ++i) {
                    h->SetBinContent(i + 1, vals[i].first);
                    h->SetBinError(i + 1, vals[i].second);
                }
                return h;
            };

            // (a) Values
            std::vector<ProfileBundle> grandBundles;
            TH1D* hGrandData = makeSeries(reorder(grandDataVals), "GrandData");
            if (hGrandData) { grandBundles.push_back({hGrandData, dataConfig}); grandToDelete.push_back(hGrandData); }

            for (size_t sysIdx = 0; sysIdx < sysVariations.size(); ++sysIdx) {
                TH1D* h = makeSeries(reorder(grandSysVals[sysIdx]), Form("GrandSys_%zu", sysIdx));
                if (h) { grandBundles.push_back({h, sysVariations[sysIdx]}); grandToDelete.push_back(h); }
            }
            for (size_t extIdx = 0; extIdx < externals.size(); ++extIdx) {
                if (externals[extIdx].basePath.empty()) continue;
                TH1D* h = makeSeries(reorder(grandExtVals[extIdx]), Form("GrandExt_%zu", extIdx));
                if (h) { grandBundles.push_back({h, externals[extIdx].config}); grandToDelete.push_back(h); }
            }
            DrawIntegratedCanvas(grandBundles, sortedLabels, "Canvas_Summary",
                                 fam.familyName + " integrated R across observables", integDir,
                                 "Integrated R", false, true);

            // (b) Differences against the data, which the per-observable folders used to show
            std::vector<ProfileBundle> grandSubBundles;
            for (size_t k = 1; k < grandBundles.size(); ++k) { // Skip the data: Data - Data is zero
                TH1D* hSub = (TH1D*)grandBundles[k].profile->Clone(Form("GrandSub_%zu", k));
                hSub->SetDirectory(nullptr);
                for (int i = 1; i <= nGrand; ++i) {
                    double d = grandDataVals[order[i-1]].first,  ed = grandDataVals[order[i-1]].second;
                    double v = grandBundles[k].profile->GetBinContent(i);
                    double ev = grandBundles[k].profile->GetBinError(i);
                    hSub->SetBinContent(i, d - v);
                    hSub->SetBinError(i, std::sqrt(ed*ed + ev*ev));
                }
                grandSubBundles.push_back({hSub, grandBundles[k].config});
                grandToDelete.push_back(hSub);
            }
            if (!grandSubBundles.empty()) {
                DrawIntegratedCanvas(grandSubBundles, sortedLabels, "Canvas_Subtracted",
                                     fam.familyName + " integrated difference across observables", integDir,
                                     "#DeltaR (Data - Var)", true, true);
            }

            // (c) Mass-window selections, hoisted out of every observable's Mass_Selection folder
            std::vector<ProfileBundle> grandMassBundles;
            if (hGrandData) grandMassBundles.push_back({hGrandData, dataConfig});
            for (size_t mIdx = 0; mIdx < massVariations.size(); ++mIdx) {
                TH1D* h = makeSeries(reorder(grandMassVals[mIdx]), Form("GrandMass_%zu", mIdx));
                if (h) { grandMassBundles.push_back({h, massVariations[mIdx]}); grandToDelete.push_back(h); }
            }
            if (grandMassBundles.size() > 1) {
                DrawIntegratedCanvas(grandMassBundles, sortedLabels, "Canvas_MassSelection",
                                     fam.familyName + " integrated R, mass windows", integDir,
                                     "Integrated R", false, true);
            }

            for (auto p : grandToDelete) delete p;
        }
    } // end of families loop

    // ---------------------------------------------------------
    // Cross-family comparison, one level above Lambda / AntiLambda / BothHyperons
    // ---------------------------------------------------------
    // Same integrated quantity, but now the SERIES are the hyperon selections and the axis is the
    // variation. One canvas per proxy, because integrated values are per entry set: a single canvas
    // mixing proxies would be comparing three different things on one axis.
    {
        TDirectory* xFamDir = nullptr;

        for (size_t iRep = 0; iRep < proxyRepresentatives.size(); ++iRep) {
            const auto& labels = crossFamilyLabels[iRep];
            if (labels.empty()) continue;

            std::vector<ProfileBundle> famBundles;
            std::vector<ProfileBundle> famSubBundles;
            std::vector<TH1*> famToDelete;

            for (size_t iFam = 0; iFam < families.size(); ++iFam) {
                const auto& vals = crossFamilyVals[iRep][iFam];
                if (vals.size() != labels.size()) continue; // Family missing this observable entirely

                // familyStyles is indexed in the same order as families, so the legend names the selection
                const VariationConfig& style = familyStyles[std::min(iFam, familyStyles.size() - 1)];

                TH1D* h = new TH1D(Form("XFam_%zu_%zu", iRep, iFam), "", (int)labels.size(), 0, (double)labels.size());
                h->SetDirectory(nullptr);
                for (size_t i = 0; i < vals.size(); ++i) {
                    h->SetBinContent((int)i + 1, vals[i].first);
                    h->SetBinError((int)i + 1, vals[i].second);
                }
                famBundles.push_back({h, style});
                famToDelete.push_back(h);

                // Difference against THIS family's own data, which is bin 1 of its own series. Subtracting
                // each selection from its own baseline is what makes the three comparable: an offset that
                // is common to a selection cancels, so what is left is how differently each systematic
                // moves Lambda, AntiLambda and BothHyperons.
                TH1D* hSub = new TH1D(Form("XFamSub_%zu_%zu", iRep, iFam), "", (int)labels.size(), 0, (double)labels.size());
                hSub->SetDirectory(nullptr);
                const double d  = vals[0].first;
                const double ed = vals[0].second;
                for (size_t i = 1; i < vals.size(); ++i) { // Bin 1 (the data itself) is left empty
                    hSub->SetBinContent((int)i + 1, d - vals[i].first);
                    hSub->SetBinError((int)i + 1, std::sqrt(ed * ed + vals[i].second * vals[i].second));
                }
                famSubBundles.push_back({hSub, style});
                famToDelete.push_back(hSub);
            }

            if (famBundles.size() > 1) {
                if (!xFamDir) xFamDir = EnsureDir(fOut, "CrossFamily_Integrated");
                const std::string proxyName = proxyRepresentatives[iRep].first;

                DrawOffsetSeriesCanvas(labels, famBundles,
                                       "Canvas_CrossFamily_" + proxyName,
                                       "Integrated R by hyperon selection (" + proxyName + " proxy)",
                                       xFamDir, "Integrated R", 0.5, true);

                DrawOffsetSeriesCanvas(labels, famSubBundles,
                                       "Canvas_CrossFamily_" + proxyName + "_Subtracted",
                                       "Integrated difference by hyperon selection (" + proxyName + " proxy)",
                                       xFamDir, "#DeltaR (Data - Var)", 0.5, true);
            }
            for (auto p : famToDelete) delete p;
        }
    }

    // ---------------------------------------------------------------------------------------
    // Cross-family comparison of the signal-extracted values
    // ---------------------------------------------------------------------------------------
    // Same shape as the block above, but the series are the two extracted quantities and the axis
    // is the hyperon selection. This is where a Lambda / AntiLambda asymmetry would show up in the
    // signal rather than in the raw measurement: the background is largely combinatorial and has
    // little reason to know which species it is sitting under, so <R>_B tracking <R>_S across the
    // families would point at something common to both rather than at a physics asymmetry.
    if (!sigExtractDir.empty()) {
        TDirectory* xFamSigDir = nullptr;

        for (size_t iProxy = 0; iProxy < sigExtractProxies.size(); ++iProxy) {
            std::vector<std::string> famLabels;
            std::vector<int> famIndices;
            for (size_t f = 0; f < families.size(); ++f) {
                bool anyFound = false;
                for (size_t iQ = 0; iQ < sigExtractQuantities.size(); ++iQ)
                    if (sigExtractCrossFamily[iQ][iProxy][f].found) anyFound = true;
                if (!anyFound) continue;
                famLabels.push_back(families[f].familyName);
                famIndices.push_back(static_cast<int>(f));
            }
            if (famLabels.size() < 2) continue; // A single family is not a cross-family plot

            const int nFam = static_cast<int>(famLabels.size());
            std::vector<ProfileBundle> series;
            std::vector<TH1*> toDelete;
            for (size_t iQ = 0; iQ < sigExtractQuantities.size(); ++iQ) {
                TH1D* h = new TH1D(Form("XFamSig_%zu_%zu", iQ, iProxy), "", nFam, 0, nFam);
                h->SetDirectory(nullptr);
                for (int i = 0; i < nFam; ++i) {
                    const auto& pt = sigExtractCrossFamily[iQ][iProxy][famIndices[i]];
                    h->SetBinContent(i + 1, pt.value);
                    h->SetBinError(i + 1, pt.error);
                }
                VariationConfig vc{"", sigExtractQuantities[iQ].legendLabel,
                                   sigExtractQuantities[iQ].color, 1,
                                   sigExtractQuantities[iQ].markerStyle, false};
                series.push_back({h, vc});
                toDelete.push_back(h);
            }

            if (!xFamSigDir)
                xFamSigDir = EnsureDir(EnsureDir(fOut, "CrossFamily_Integrated"), "SigExtracted");

            const std::string proxyName = sigExtractProxies[iProxy].folderName;
            DrawOffsetSeriesCanvas(famLabels, series,
                                   "Canvas_CrossFamily_SigExtracted_" + proxyName,
                                   "Signal-extracted integrated <R> by hyperon selection (" +
                                       proxyName + " proxy)",
                                   xFamSigDir, "Integrated <R>", 0.5, true);
            for (auto q : toDelete) delete q;
        }
    }

    cache.CloseAll(); // Release every input file before closing the output

    fOut->Close();
    delete fOut;
    std::cout << " Done!\n";
}

// To be able to run this both as a compiled C++ macro or a ROOT script:
#ifndef __CINT__
int main(int argc, char** argv) {
    if (argc < 2) { // Check is argc < 2 because this is the bare minimum. mcRefDir, ppRefDir, toyModelPath and cutFolder are all optionals
        std::cerr << "Usage: " << argv[0]
                  << " <consumerDir> [mcRefDir] [ppRefDir] [toyModelPath] [cutFolder]"
                     " [doIndividualComparisons] [sigExtractDir]\n"
                     "  sigExtractDir defaults to <consumerDir>/../results_SigExtract;"
                     " pass \"none\" to skip the signal-extraction plots.\n";
        return 1;
    }
    std::string consumerDir = argv[1];
    std::string mcRefDir = (argc > 2) ? argv[2] : "";
    std::string ppRefDir = (argc > 3) ? argv[3] : "";
    std::string toyModelPath = (argc > 4) ? argv[4] : "";
    std::string cutFolder = (argc > 5) ? argv[5] : DEFAULT_CUT_FOLDER;
    // signalExtractionRing writes next to the consumer output, so the default follows the pipeline
    // layout (<workDir>/results_consumer and <workDir>/results_SigExtract) and can be overridden.
    // Passing "none" disables the signal-extraction plots entirely.
    std::string sigExtractDir = (argc > 7) ? argv[7] : (consumerDir + "/../results_SigExtract");
    if (sigExtractDir == "none") sigExtractDir = "";
    // Accepts 1/0, true/false, yes/no. Anything else (including an absent argument) leaves it off.
    bool doIndividual = false;
    if (argc > 6) {
        std::string flag = argv[6];
        doIndividual = (flag == "1" || flag == "true" || flag == "yes" || flag == "on");
    }

    auxiliarySummaryPlots(consumerDir, mcRefDir, ppRefDir, toyModelPath, cutFolder, doIndividual,
                          sigExtractDir);
    return 0;
}
#endif