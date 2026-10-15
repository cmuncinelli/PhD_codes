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
* Observables are declared in a single table (see observableGroups in the main
* function): each entry carries its own in-file directory, so profiles living in
* different folders of the consumer output (Ring/<Proxy>/EtaDependence/,
* Ring/<Proxy>/PtDependence/, or task-level folders such as FakePolDrivers/) all
* flow through the same machinery.
*
* Usage: root -l -b -q 'auxiliarySummaryPlots.cxx("path/to/wagon/dir", "[optional_MC_ref_path]")'
* (or just use the bash coordinator script, which properly compiles this into production-level code)
* ============================================================
*/

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
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
#include <cctype>     // std::isalnum, for canvas-safe names
#include "TStyle.h"
#include "TMath.h"
#include "TLine.h"
#include "TString.h"
#include "TSystem.h"
#include <sstream> // IsLongitudinalSuffix
#include "TDirectory.h"
#include "TExec.h"  // Carries a map's palette inside its canvas; gStyle is not serialised
#include "TKey.h"       // Longitudinal relabel pass: walking the output file
#include "TClass.h"
#include "TGraph.h"
#include "TText.h"
#include "TLegendEntry.h"
#include "TObjString.h"
#include "THashList.h"

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
/// @brief A longitudinal ring definition of the consumer (README, "The longitudinal ring"). See LONGITUDINAL_RINGS.
struct LongitudinalRing {
    std::string name;         // "RingZ": the output folders are <Family>_<name>
    std::string fileTag;      // "_useRingZ": right after the family in the consumer config and file names
    std::string symbolSuffix; // "_{z}": appended to the ring symbol in display text
    std::string shortName;    // "Rz": used in canvas names
};

struct FamilyConfig {
    std::string familyName;                 // e.g., "Lambda"
    std::string dataSuffix;                 // e.g., "JustLambda"
    double alphaSq;                         // Squared decay constant of the selection (see AlphaSqBoth for the mixture)
    const LongitudinalRing* ring = nullptr; // Set for the longitudinal siblings of a full-ring family
};

/// @brief One 1D TProfile to extract, together with where it lives and how to label it.
struct ProfileConfig {
    // In-file path below the task dir, spelled exactly as the consumer books it (so it can be grepped in
    // both codes), e.g. "Ring/LeadP/EtaDependence/pRingObservableProxyEta". Its proxy folder (LeadJet, LeadP
    // or SubJet) names the proxy that defines the ring: see ProxyOfPath.
    std::string path;

    // One-vs-one canvases against every single variation:
    // (By far the largest contributor to the ~20 extra canvases per observable)
    // Set only on the LeadP entries.
    bool doIndividualComparisons = false;

    // Name of the matching Toy Model profile inside TOY_SCENARIO_DIR. Empty (the default) means the toy
    // has no counterpart for this observable, and the toy is then skipped for it rather than guessed at.
    std::string toyProfile = "";

    // Fit R(eta) = p0 * tanh(p1 * eta) to this observable, once per variation, plus the cross-variation
    // parameter summaries. Only meaningful for the fine-binned eta profiles: the coarse ones do not have
    // enough points across the turn-over for the two parameters to separate.
    bool doTanhFit = false;

    // Whether this observable joins the brute-force axis, where it is grouped by its proxy. False keeps it
    // off: that is what the eta-sign split observables get, since they cover different subsets and would
    // stretch the scale for everything else.
    bool inBruteForce = true;

    // Folder and label name for this observable. Empty means "derive it" -- see
    // ObservableDisplayName. Set it only where the derived name is still unwieldy.
    std::string displayName = "";

    // Appended to the x title read from the data profile. Only for observables whose selection the consumer
    // states in the histogram title rather than on the axis (the eta-sign splits of the proxy pT).
    std::string xTitleSuffix = "";

    // Axis titles. Left empty in the observable tables: the observable loop fills them from the data
    // profile, i.e. from the consumer's own booking, with the proxy as a subscript of the y title. Set
    // explicitly only where a ProfileConfig just carries titles for a drawer.
    std::string xAxisTitle = "";
    std::string yAxisTitle = "";
};

/**
 * @brief Folder name for one observable: the explicit displayName if set, else a derived one.
 *
 * The consumer's object names encode their ROOT type and a long "what is this" prefix, which made
 * output folders like Lambda/EtaProxy/pRingObservableProxyEta -- almost all of it repeated from the
 * group folder one level up. Stripping mechanically handles the bulk; the cases where the result is
 * still unreadable carry an explicit displayName in the table, where the judgment is visible.
 * The object names are the same in every proxy folder, so ObservableDisplayName then puts the proxy
 * in (WithProxy).
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

// The consumer books every ring observable under a single folder, Ring/<Proxy>/..., and the old
// kinematic-cut families (RingKinematicCuts, JetKinematicCuts, ...) no longer exist. Kept as a
// constant so that every path below spells it the same way.
const std::string RING_DIR_IN_FILE = "Ring/";
// signalExtractionRing mirrors that folder at the top of its own output (Ring/IntegratedSummary/,
// Ring/KappaEff/), so every read of an extraction file goes through this one.
const std::string SIG_EXTRACT_RING_DIR = "Ring/";

// The consumer's three proxy folders, in the order every per-proxy loop and axis of this macro uses.
const std::vector<std::string> PROXIES = {"LeadP", "LeadJet", "SubJet"};

// Toy Model scenario to compare against. The toy books a 2 x 4 x 3 grid of
// {WithEtaGate, WithoutEtaGate} x {NoCuts, pTCutOnly, DCACutOnly, BothCuts} x {EtaPos, EtaNeg, All}.
// WithEtaGate/BothCuts is the internally consistent one (the toy's own header notes that
// WithoutEtaGate exists only to demonstrate why omitting the daughter acceptance gate makes the
// physics cuts inconsistent), and "All" is the unsplit sample that matches how the consumer fills.
const std::string TOY_SCENARIO_DIR = "WithEtaGate/BothCuts/All/";

// Toy Model integrated <R>, split by the sign of its (random) jet direction. These three single-bin
// profiles are the toy counterparts of the first three bins of the ProxyEta cut scans
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
//                                                      Ring/<Proxy>/PtDependence/pRingObservableLambdaPt.
//   pRingProxyJet_JetEtaPos <R> integrated, eta_jet >= 0 -> UNUSED here, but these two are the natural
//   pRingProxyJet_JetEtaNeg <R> integrated, eta_jet <  0    toy counterparts of the ProxyEta eta-split
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

// ---------------------------------------------------------
// Skip log: nothing is dropped silently, and nothing floods the log
// ---------------------------------------------------------
// One consumer switch turned off can take the same object out of every configuration file, and a warning per
// read would bury everything else. So every skip is counted here, per reason and per object, together with the
// distinct files or canvases it happened in, and PrintSkipSummary lists each (reason, object) once, at the end.
//   CountSkip  expected absences (a file not run, an object behind a consumer switch): counted only.
//   WarnSkip   anomalies, and objects every file should hold: also printed the first time they are met, so they
//              show up in context. Their repeats are only counted.

/// @brief How often one (reason, object) was skipped, and in which distinct files or canvases.
struct SkipEntry { int count = 0; std::set<std::string> places; };
std::map<std::pair<std::string, std::string>, SkipEntry> gSkips; // (reason, object) --> its entry

/// @brief Counts one skip for the end-of-run summary. @return True the first time this (reason, object) is seen.
bool CountSkip(const std::string& what, const std::string& why, const std::string& place = "") {
    SkipEntry& e = gSkips[{why, what}];
    if (!place.empty()) e.places.insert(place);
    return ++e.count == 1;
}

/// @brief Counts one skip, and prints it the first time this (reason, object) is met.
void WarnSkip(const std::string& what, const std::string& why, const std::string& place = "") {
    if (!CountSkip(what, why, place)) return;
    std::cerr << "  [Warning] Skipped " << what << (place.empty() ? "" : " in " + place) << ": " << why
              << " (repeats: see the skip summary at the end)" << std::endl;
}

/// @brief The end-of-run summary: one line per (reason, object), with how often it was skipped and where.
void PrintSkipSummary() {
    std::cout << "========================================================\n";
    if (gSkips.empty()) {
        std::cout << " Skip summary: nothing was skipped\n";
        return;
    }
    std::cout << " Skip summary (one line per reason and object, over the whole run)\n";
    for (const auto& kv : gSkips) {
        const SkipEntry& e = kv.second;
        std::cout << "  " << kv.first.first << ": " << kv.first.second << " -- " << e.count
                  << (e.count == 1 ? " time" : " times");
        if (e.places.size() == 1) {
            std::cout << ", in " << *e.places.begin();
        } else if (e.places.size() > 1) {
            bool allFiles = true; // Files for the absences, canvases or folders for everything else
            for (const auto& p : e.places)
                allFiles = allFiles && p.size() > 5 && p.compare(p.size() - 5, 5, ".root") == 0;
            std::cout << ", in " << e.places.size() << (allFiles ? " files" : " places");
        }
        std::cout << "\n";
    }
    std::cout << "========================================================\n";
}

/**
 * @brief Returns an existing subdirectory, creating it only if absent.
 *
 * Used to keep directory creation lazy: nothing is created until we actually have a canvas to
 * write into it, so inputs that lack a given object do not leave a trail of empty folders behind.
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
 * Failed opens are cached too, so a missing file is reported once and then only counted (see the skip log).
 *
 * @note Objects handed out by FetchClone are detached clones, so they stay valid after CloseAll().
 */
class FileCache {
public:
    ~FileCache() { CloseAll(); }

    /**
     * @brief Opens a file, or returns the already-open handle.
     * @param required False for files whose absence is common (a configuration that was not run): their
     *                 absence is then only counted, never printed on the spot.
     * @return The open file, or nullptr if it could not be opened. A failed path is remembered and keeps
     *         returning nullptr without retrying; each such read is counted in the skip log, except for the
     *         paths marked absent on purpose (MarkAbsent), which the run banner already lists.
     */
    TFile* Open(const std::string& path, bool required = true) {
        auto it = fFiles.find(path);
        if (it != fFiles.end()) { // Includes cached nullptr for known-bad paths
            if (!it->second && !fExpectedAbsent.count(path)) CountSkip(path, "file absent or unreadable");
            return it->second;
        }

        // Checked first, so an absent file costs no ROOT error message (it goes to the skip log instead)
        TFile* f = gSystem->AccessPathName(path.c_str()) ? nullptr : TFile::Open(path.c_str(), "READ");
        if (f && f->IsZombie()) {
            f->Close();
            delete f;
            f = nullptr;
        }
        if (!f) {
            if (required) WarnSkip(path, "file absent or unreadable");
            else CountSkip(path, "file absent or unreadable");
        }
        fFiles[path] = f;
        return f;
    }

    /**
     * @brief Fetches an object and returns a detached clone that survives CloseAll().
     * @param objPath        Path of the object inside the file, task folder included.
     * @param warnIfMissing  Set false for genuinely optional objects (booked behind a consumer switch, or in
     *                       files that may not exist): their absence is then only counted for the skip
     *                       summary. True prints the first absence of each object as well.
     * @return A clone owned by the caller (delete it), or nullptr if the file or object is missing.
     */
    template <typename T>
    T* FetchClone(const std::string& filePath, const std::string& objPath, bool warnIfMissing = true) {
        TFile* f = Open(filePath, warnIfMissing);
        if (!f) return nullptr;

        T* obj = dynamic_cast<T*>(f->Get(objPath.c_str()));
        if (!obj) {
            // The task folder is the same in every path: left out of the log, which stays readable
            const std::string shown = (objPath.rfind(TASK_DIR_IN_FILE, 0) == 0) ? objPath.substr(TASK_DIR_IN_FILE.size())
                                                                                 : objPath;
            if (warnIfMissing) WarnSkip(shown, "object absent", filePath);
            else CountSkip(shown, "object absent (optional)", filePath);
            return nullptr;
        }

        T* clone = static_cast<T*>(obj->Clone());
        clone->SetDirectory(nullptr); // Detach immediately: Clone() otherwise attaches to gDirectory,
                                      // which at this point is the output file
        return clone;
    }

    /// @brief Records a path as absent, with no warning and no count: for files the run banner already lists.
    void MarkAbsent(const std::string& path) {
        fFiles.emplace(path, nullptr);
        fExpectedAbsent.insert(path);
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
    std::set<std::string> fExpectedAbsent; // MarkAbsent paths: listed in the banner, so not counted again
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
 * @note REQUIREMENT: "same binning" means the same bin count and the same axis range. The bin edges
 *       themselves are not compared. That holds for every current call (identical consumer bookings,
 *       slices of one 2D profile, and the Toy Model, which takes the interpolation branch), but a
 *       reference produced with a different variable-width axis of the same count and range would be
 *       subtracted bin by bin across different edges, without a warning.
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
    // ProjectionX does not always reproduce them, and the EtaSplitStudies difference/pull panels are unreadable
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
//       the TProfile1Ds IntegratedCuts/<Proxy>/pRingCuts, one per proxy.
//       Re-integrating a multi-bin profile here recomputes a number that already exists, by a
//       different route, which is exactly how two values start disagreeing after an innocent
//       rebin. Read the consumer's value instead and keep this only as a cross-check.
//       (Companion note: the two-dimensional IntegratedCuts/<Proxy>/p2dRingCutsV0MassPeak carry the same
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
 * set on the ProxyEta cut scans. Those labels ARE the axis there, so those canvases get this
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
 * These live in FakePolDrivers/ProxyEta/<Proxy>/ and are booked only when the consumer's
 * doFakePolDiagnosticsQA switch is on, so they are always fetched defensively.
 */
struct CategoricalObservable {
    std::string path;         // In-file path below the task dir, e.g. "FakePolDrivers/ProxyEta/LeadP/pRingEtaCuts"

    // Build the three-bin "All / eta_proxy >= 0 / eta_proxy < 0" digest of bins 1-3, which is the only
    // slice of these plots the Toy Model can be overlaid on. See kToyEtaSign* above.
    bool doEtaSignSummary = false;

    // Companion TProfile2D of shape (same bins) x (out of peak, in peak), projected into two curves. Its
    // object name only: it sits in the same folder. Empty means there is none. Only the LeadP scan has one.
    std::string massSplitProfile = "";

    // Same convention as ProfileConfig::displayName. Empty means derive it.
    std::string displayName = "";
};

/// @brief The object name at the end of an in-file path ("Ring/LeadP/EtaDependence/pRingObservableProxyEta" --> "pRingObservableProxyEta").
std::string ObjectNameOf(const std::string& path) {
    const size_t slash = path.rfind('/');
    return (slash == std::string::npos) ? path : path.substr(slash + 1);
}

/// @brief The folder part of an in-file path, trailing slash included ("" for a bare object name).
std::string FolderOf(const std::string& path) {
    const size_t slash = path.rfind('/');
    return (slash == std::string::npos) ? "" : path.substr(0, slash + 1);
}

/// @brief The proxy folder named in an in-file path ("Ring/LeadP/..." --> "LeadP"), or "" if there is none.
std::string ProxyOfPath(const std::string& path) {
    const std::string padded = "/" + path + "/"; // Whole folder names only
    for (const auto& proxy : PROXIES)
        if (padded.find("/" + proxy + "/") != std::string::npos) return proxy;
    return "";
}

/**
 * @brief Puts the proxy into a display name: "ProxyEta" --> "LeadPEta", "PVz" --> "LeadPPVz".
 *
 * Every proxy folder of the consumer holds the same object names, so a name derived from the object alone
 * would send the three proxies into one output folder, as new cycles of the same canvases (no error, and
 * the browser shows only the last). Every "Proxy" token is replaced by the proxy's name; a name without
 * one gets the proxy as a prefix, unless it already names it (the explicit AEE display names do).
 */
std::string WithProxy(std::string name, const std::string& proxy) {
    if (proxy.empty()) return name;
    const std::string kToken = "Proxy";
    if (name.find(kToken) == std::string::npos)
        return (name.find(proxy) != std::string::npos) ? name : proxy + name;
    for (size_t pos = name.find(kToken); pos != std::string::npos; pos = name.find(kToken, pos + proxy.size()))
        name.replace(pos, kToken.size(), proxy);
    return name;
}

/// @brief A y title that names its proxy: "<#it{R}>" --> "<#it{R}>_{LeadP}". Unchanged if it already does, or if either is empty.
std::string WithProxySubscript(const std::string& title, const std::string& proxy) {
    if (proxy.empty() || title.empty() || title.find("_{" + proxy + "}") != std::string::npos) return title;
    return title + "_{" + proxy + "}";
}

/// @brief Folder and label name for a ProfileConfig: derived (or explicit), then carrying its proxy.
std::string ObservableDisplayName(const ProfileConfig& profConfig) {
    return WithProxy(DeriveObservableName(ObjectNameOf(profConfig.path), profConfig.displayName),
                     ProxyOfPath(profConfig.path));
}

/// @brief Folder name for a CategoricalObservable, by the same rule as ProfileConfig.
std::string ObservableDisplayName(const CategoricalObservable& catObs) {
    return WithProxy(DeriveObservableName(ObjectNameOf(catObs.path), catObs.displayName),
                     ProxyOfPath(catObs.path));
}

/**
 * @brief Copies the first three bins of a cut-scan into a standalone three-bin histogram.
 *
 * Bins 1-3 of every ProxyEta scan are "All Lambda", "eta_proxy >= 0" and "eta_proxy < 0". Pulling them
 * out lets the Toy Model, which only knows those three numbers, be compared against them on an axis
 * where every source is fully populated, rather than leaving 6 or 12 empty bins on the full scan.
 *
 * @return A detached TH1D owned by the caller, or nullptr if the source has fewer than three bins.
 */
TH1D* ExtractEtaSignBins(TH1* src, const std::string& name) {
    if (!src) return nullptr;
    if (src->GetNbinsX() < 3) {
        WarnSkip(src->GetName(), "eta-sign digest needs at least three bins, the scan has " + std::to_string(src->GetNbinsX()));
        return nullptr;
    }

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
constexpr double kOverlaySpreadFrac = 0.4;

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
 *                     Needed for the ProxyEta cut scans, whose x-axis is a list of named selections.
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
        if (!g) { // A variation with no usable bin is dropped rather than drawn empty
            WarnSkip(bundles[i].config.legendLabel, "curve with no usable bin", std::string(outDir->GetPath()) + "/" + canvasName);
            continue;
        }
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

// ------------------------------------------------------------------------------------------------
// Mass regions: styling and bins, shared by every canvas that slices on the mass-region axis.
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

// The consumer's mass-region axis (axisMassRegions) has exactly five bins: left sideband, left gap, peak,
// right gap, right sideband. Masses outside the sidebands go to the under/overflow, which is never read.
// Only the sidebands and the peak are drawn for now. TODO: the gap bins (2 and 4): draw them, or not?
// An axis with any other number of bins is refused with a warning (HasMassRegionAxis), never sliced by guesswork.
constexpr int kNMassRegionBins = 5;
struct MassRegion { int bin; const char* label; MassRegionStyle style; };
constexpr MassRegion kMassLeftSideband  = {1, "left sideband",  kStyleLeftSideband};
constexpr MassRegion kMassPeak          = {3, "in mass window", kStyleInPeak};
constexpr MassRegion kMassRightSideband = {5, "right sideband", kStyleRightSideband};
constexpr MassRegion kMassRegions[3] = {kMassLeftSideband, kMassPeak, kMassRightSideband}; // Drawing order

/// @brief True for the consumer's five-bin mass-region axis (on y); any other axis is refused with a warning.
bool HasMassRegionAxis(TProfile2D* p2d, const std::string& what) {
    const int n = p2d->GetYaxis()->GetNbins();
    if (n == kNMassRegionBins) return true;
    WarnSkip(what, std::to_string(n) + " mass bins, the mass-region axis has " + std::to_string(kNMassRegionBins) +
                       " (left sideband, left gap, peak, right gap, right sideband)");
    return false;
}

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
 * @param p2d      Source, with the observable axis on x and the five-bin mass-region axis on y (see
 *                 kMassRegions). Not modified. Any other mass axis is refused with a warning.
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
    if (!HasMassRegionAxis(p2d, canvasTitle)) return;

    std::vector<ProfileBundle> bundles;
    std::vector<TH1*> toDelete;

    for (const MassRegion& region : kMassRegions) {
        const int m = region.bin;
        TProfile* slice = p2d->ProfileX(Form("%s_slice%d", canvasName.c_str(), m), m, m);
        if (!slice) continue;
        slice->SetDirectory(nullptr);
        if (slice->GetEntries() < 1) {
            WarnSkip(canvasTitle + ", " + region.label, "no entries");
            delete slice;
            continue;
        }

        const double mLo = p2d->GetYaxis()->GetBinLowEdge(m);
        const double mHi = p2d->GetYaxis()->GetBinUpEdge(m);

        // The in-window slice is the reference the sidebands are read against
        VariationConfig cfg = MassRegionConfig(Form("%s [%.4f, %.4f]", region.label, mLo, mHi),
                                               region.style, (m == kMassPeak.bin));

        bundles.push_back({slice, cfg});
        toDelete.push_back(slice);
    }
    if (bundles.size() < 2) {
        WarnSkip(canvasTitle, "fewer than two mass regions with entries");
        for (auto p : toDelete) delete p;
        return;
    }

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
 * @brief The two questions a mass-region axis can answer, as two canvases.
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
 * @param p2d Source with the observable on x and the five-bin mass-region axis on y (see kMassRegions).
 *            Any other mass axis is refused with a warning.
 */
void DrawMassRegionCanvases(TProfile2D* p2d,
                            const std::string& titlePrefix,
                            TDirectory* outDir,
                            const std::string& xTitle,
                            const std::string& yTitle) {
    if (!p2d || !outDir) return;
    if (!HasMassRegionAxis(p2d, titlePrefix)) return;

    TProfile* sLeft  = p2d->ProfileX("mrc_left",  kMassLeftSideband.bin,  kMassLeftSideband.bin);
    TProfile* sPeak  = p2d->ProfileX("mrc_peak",  kMassPeak.bin,          kMassPeak.bin);
    TProfile* sRight = p2d->ProfileX("mrc_right", kMassRightSideband.bin, kMassRightSideband.bin);
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
    // consumer files produced with a mass cut applied; this one slices the mass-region axis of a
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
    } else {
        WarnSkip(titlePrefix + ", left vs right sideband", "a sideband has no entries");
    }

    delete sLeft; delete sPeak; delete sRight;
}

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
 * @param bottomMargin Canvas bottom margin without the pull pad. The default leaves room for long
 *                   angled labels; short column labels can afford much less.
 */
void DrawOffsetSeriesCanvas(const std::vector<std::string>& labels,
                            const std::vector<ProfileBundle>& series,
                            const std::string& canvasName,
                            const std::string& canvasTitle,
                            TDirectory* outDir,
                            const std::string& yTitle,
                            double spreadFrac = 0.5,
                            bool addSelfPullPad = false,
                            double bottomMargin = 0.22) {
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
        c->SetBottomMargin(bottomMargin); // Room for the angled labels
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

    // Grows with the series count past three, so canvases with up to three series are unchanged
    const double legHeight = std::max(0.10, 0.034 * nSeries);
    TLegend* leg = new TLegend(0.14, 0.90 - legHeight, 0.44, 0.90);
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

// =================================================================================================
// SIGNAL-EXTRACTED INTEGRATED OBSERVABLES
// =================================================================================================
// Everything below reads the output of signalExtractionRing.cxx, which is a DIFFERENT file per
// consumer configuration:
//     <sigExtractDir>/signalExtractionRing_<dataSuffix><sysSuffix>.root
// and, inside it,
//     Ring/IntegratedSummary/<proxy>/<histogram>
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

/// @brief Colour for one eta_Proxy slice: blue below zero, red above, deeper the further out.
/// The sign is the hue, so an antisymmetry reads at a glance; |eta| is the depth, so neighbouring
/// slices stay apart. Used together with an open/filled marker split by sign, which survives a
/// colour-blind reader and a black-and-white print.
// Fixed ROOT colours rather than new TColor indices: the canvases do not store custom colour indices consistently.
int EtaSliceColor(double etaCentre, double etaMax)
{
    const double t = std::min(1.0, std::abs(etaCentre) / etaMax);

    // Four shades on each side, from lighter near eta = 0
    // to darker towards the edge.
    static const int blue[] = {
        kAzure + 1,
        kAzure + 2,
        kBlue + 1,
        kBlue + 2
    };

    static const int red[] = {
        kOrange + 1,
        kRed + 1,
        kRed + 2,
        kRed + 3
    };

    const int n = sizeof(blue) / sizeof(blue[0]);
    int i = std::min(n - 1, static_cast<int>(t * n));

    return (etaCentre < 0.0) ? blue[i] : red[i];
}

/**
 * @brief A per-cell map over (AEE angle, eta_Proxy), with the cells that produced no result hatched.
 *
 * A cell that was not extracted holds zero, which in a signed quantity paints exactly like a measured
 * zero. Only the status map tells the two apart, so it is overlaid here rather than left to a separate
 * canvas. Its code 1 means "extracted", as in signalExtractionRing's CellStatus.
 *
 * The palette travels with the canvas in a TExec: gStyle is not serialised, so without it a map would
 * reopen in whatever palette the reader has active. Signed quantities get a diverging palette on a
 * range symmetric about zero, so zero is the neutral colour and the sign reads from the hue.
 */
void DrawMapCanvas(TH2* hMap, TH2* hStatus, const std::string& canvasName,
                   const std::string& canvasTitle, TDirectory* outDir,
                   bool isSigned, bool drawText = false) {
    if (!hMap || !outDir) return;
    outDir->cd();

    TH2* h = (TH2*)hMap->Clone((canvasName + "_map").c_str());
    h->SetDirectory(nullptr);
    h->SetStats(0);
    h->GetZaxis()->SetTitleOffset(1.3);
    if (isSigned) {
        double aMax = 0.0;
        for (int ix = 1; ix <= h->GetNbinsX(); ++ix)
            for (int iy = 1; iy <= h->GetNbinsY(); ++iy)
                aMax = std::max(aMax, std::fabs(h->GetBinContent(ix, iy)));
        if (aMax > 0.0) { h->SetMinimum(-aMax); h->SetMaximum(aMax); }
    }

    TCanvas* c = new TCanvas(canvasName.c_str(), canvasTitle.c_str(), 1000, 650);
    c->SetLeftMargin(0.11);
    c->SetRightMargin(0.16);
    c->SetBottomMargin(0.12);

    const std::string opt = drawText ? "COLZ TEXT" : "COLZ";
    TExec* palette = new TExec("ex_palette", isSigned ? "gStyle->SetPalette(kLightTemperature);"
                                                      : "gStyle->SetPalette(kBird);");
    palette->SetBit(kCanDelete);
    // Draw, palette, draw again: the TExec only acts on what is painted after it
    h->Draw(opt.c_str());
    palette->Draw();
    h->Draw((opt + " SAME").c_str());

    if (hStatus) {
        for (int ix = 1; ix <= hStatus->GetNbinsX(); ++ix) {
            for (int iy = 1; iy <= hStatus->GetNbinsY(); ++iy) {
                if (std::lround(hStatus->GetBinContent(ix, iy)) == 1) continue;
                TBox* b = new TBox(hStatus->GetXaxis()->GetBinLowEdge(ix), hStatus->GetYaxis()->GetBinLowEdge(iy),
                                   hStatus->GetXaxis()->GetBinUpEdge(ix), hStatus->GetYaxis()->GetBinUpEdge(iy));
                b->SetFillColor(kGray + 1);
                b->SetFillStyle(3354);
                b->SetLineWidth(0);
                b->SetBit(kCanDelete);
                b->Draw();
            }
        }
    }

    c->Update(); // TCanvas::Write() does not paint
    c->Write();
    delete c;    // Before h, which the canvas still lists but does not own
    delete h;
    gStyle->SetPalette(kBird); // Leave nothing behind for the next canvas
}

/**
 * @brief One integrated number per variation, drawn categorically, in its full and redux subsets.
 *
 * Entry 0 is the reference and is kept in every subset. Each subset gets the values and their
 * Data - Variation differences, the latter in quadrature: every entry is a separate consumer run,
 * which is the one situation where that is right (the mixed-event wagons reuse the V0 sample, so
 * they over-cover slightly, as accepted throughout).
 */
void DrawCategoricalCrossSystem(const std::vector<std::string>& labels,
                                const std::vector<VariationConfig>& cfgs,
                                const std::vector<std::pair<double, double>>& vals,
                                const std::string& nameStem, const std::string& title,
                                TDirectory* outDir, const std::string& yTitle) {
    if (!outDir || labels.size() < 2) return;

    auto emit = [&](std::function<bool(const VariationConfig&)> keep, const std::string& suffix) {
        std::vector<int> pick;
        for (size_t i = 0; i < cfgs.size(); ++i)
            if (i == 0 || keep(cfgs[i])) pick.push_back(static_cast<int>(i));
        if (pick.size() < 2) return;

        const int n = static_cast<int>(pick.size());
        std::vector<std::string> lb;
        std::vector<ProfileBundle> bV, bS;
        std::vector<TH1*> del;
        for (int k = 0; k < n; ++k) {
            const int i = pick[k];
            lb.push_back(labels[i]);
            TH1D* hV = MakeCategoricalPoint(Form("%s%s_V_%d", nameStem.c_str(), suffix.c_str(), k),
                                            k, n, vals[i].first, vals[i].second);
            bV.push_back({hV, cfgs[i]});
            del.push_back(hV);
            if (k == 0) continue;
            const double dd = vals[0].first - vals[i].first;
            const double ee = std::sqrt(vals[0].second * vals[0].second + vals[i].second * vals[i].second);
            TH1D* hS = MakeCategoricalPoint(Form("%s%s_S_%d", nameStem.c_str(), suffix.c_str(), k),
                                            k, n, dd, ee);
            bS.push_back({hS, cfgs[i]});
            del.push_back(hS);
        }
        DrawIntegratedCanvas(bV, lb, nameStem + suffix, title, outDir, yTitle, false, false, true);
        if (!bS.empty())
            DrawIntegratedCanvas(bS, lb, nameStem + suffix + "_Subtracted", title + ", difference",
                                 outDir, "#Delta" + yTitle + " (Data - Var)", true, false, true);
        for (auto p : del) delete p;
    };
    // Full and redux from the same call, as everywhere else
    emit([](const VariationConfig&) { return true; }, "");
    emit([](const VariationConfig& c) { return c.inRedux; }, "_Redux");
}

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
 * Bin 1 is "All Lambda". Bins 2-4 hold the old kinematic-cut selections (Lambda cuts, proxy |eta| < 0.5,
 * both), unused for now.
 * README: the two-dimensional companion, IntegratedCuts/<Proxy>/p2dRingCutsV0MassPeak, shares this X axis
 * crossed with a mass flag -- Y bin 1 strictly OUT of the mass peak, Y bin 2 strictly IN it.
 *
 * @param proxy "LeadJet", "LeadP" or "SubJet": the folder under IntegratedCuts/.
 * @return {value, error}, and found = false when the object is absent.
 */
SigExtractPoint FetchIntegratedFromCuts(FileCache& cache, const std::string& filePath,
                                        const std::string& proxy) {
    SigExtractPoint out;
    // Every proxy folder holds the same object name, so the proxy only picks the folder
    TProfile* p = cache.FetchClone<TProfile>(filePath, TASK_DIR_IN_FILE + "IntegratedCuts/" + proxy + "/pRingCuts", false);
    if (!p) return out;
    if (p->GetNbinsX() >= 1 && p->GetBinEntries(1) > 0) {
        out.value = p->GetBinContent(1);
        out.error = p->GetBinError(1);
        out.found = true;
    } else {
        WarnSkip("IntegratedCuts/" + proxy + "/pRingCuts, bin 1", "no entries", filePath);
    }
    delete p;
    return out;
}

/**
 * @brief The integrated <R> split by mass window, from IntegratedCuts/<Proxy>/p2dRingCutsV0MassPeak.
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
    // Every proxy folder holds the same object name, so the proxy only picks the folder
    TProfile2D* p = cache.FetchClone<TProfile2D>(filePath,
                                                 TASK_DIR_IN_FILE + "IntegratedCuts/" + proxy + "/p2dRingCutsV0MassPeak", false);
    if (!p) return out;
    const int yBin = inPeak ? 2 : 1;
    const int gBin = p->GetBin(1, yBin); // X bin 1 is "All Lambda"
    if (p->GetBinEntries(gBin) > 0) {
        out.value = p->GetBinContent(gBin);
        out.error = p->GetBinError(gBin);
        out.found = true;
    } else {
        WarnSkip("IntegratedCuts/" + proxy + "/p2dRingCutsV0MassPeak, " + (inPeak ? "in peak" : "out of peak"), "no entries", filePath);
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
 * @param proxy         Folder under Ring/IntegratedSummary/ (or under Ring/KappaEff/, see summaryFolder).
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
                                     const std::string& proxy,
                                     const std::string& histName,
                                     bool requireValid = true,
                                     const std::string& summaryFolder = "IntegratedSummary") {
    SigExtractPoint out;
    if (sigExtractDir.empty()) return out;

    // summaryFolder is "KappaEff" for the moment extractions, whose "proxy" is <Proxy>_<Num|Den>
    const std::string filePath = sigExtractDir + "/signalExtractionRing_" + fileSuffix + ".root";
    const std::string dirPath = SIG_EXTRACT_RING_DIR + summaryFolder + "/" + proxy + "/";

    if (requireValid) {
        TH1D* hStatus = cache.FetchClone<TH1D>(filePath, dirPath + "hExtractionStatus_" + proxy, false);
        if (!hStatus) return out;           // No status object at all: treat as not extracted
        const bool valid = hStatus->GetBinContent(1) > 0.5; // Bin 1 is labelled "valid"
        delete hStatus;
        if (!valid) {
            WarnSkip(dirPath + histName + "_" + proxy, "extraction flagged invalid (hExtractionStatus)", filePath);
            return out;
        }
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
        SigExtractPoint p = FetchSigExtractPoint(cache, sigExtractDir, fileSuffix,
                                                 proxy.folderName, components[i].histName);
        if (!p.found) return; // All three come from one extraction: if one is absent, so are the rest
        labels.push_back(components[i].legendLabel);
    }

    const int nCats = static_cast<int>(components.size());
    for (int i = 0; i < nCats; ++i) {
        SigExtractPoint p = FetchSigExtractPoint(cache, sigExtractDir, fileSuffix,
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
        SigExtractPoint p = FetchSigExtractPoint(cache, sigExtractDir, fileSuffix,
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
                                          const SigExtractProxy& proxy,
                                          const SigExtractQuantity& quantity,
                                          const std::string& familyName,
                                          TDirectory* outDir) {
    SigExtractPoint dataPoint;
    if (!outDir) return dataPoint;

    dataPoint = FetchSigExtractPoint(cache, sigExtractDir, dataSuffix,
                                     proxy.folderName, quantity.histName);
    if (!dataPoint.found) return dataPoint;

    std::vector<std::string> labels{dataConfig.legendLabel};
    std::vector<VariationConfig> configs{dataConfig};
    std::vector<SigExtractPoint> points{dataPoint};

    for (const auto& sys : sysVariations) {
        SigExtractPoint p = FetchSigExtractPoint(cache, sigExtractDir, dataSuffix + sys.suffix,
                                                 proxy.folderName, quantity.histName);
        if (!p.found) continue; // A variation whose extraction failed drops off the axis entirely
        labels.push_back(sys.legendLabel);
        configs.push_back(sys);
        points.push_back(p);
    }

    if (labels.size() < 2) { // Data alone is not a systematics plot
        CountSkip("SignalExtraction, " + quantity.legendLabel, "Data alone, nothing to compare it with",
                  familyName + "/" + proxy.folderName);
        return dataPoint;
    }

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

// ---------------------------------------------------------
// Longitudinal families (useRingZ and useRingZT consumer outputs). See the README, "The longitudinal ring".
// ---------------------------------------------------------
// A longitudinal consumer output carries its tag right after the family suffix, e.g.
// ConsumerResults_BothHyperons_useRingZT_MixedEventProxies.root. A family whose dataSuffix ends in it
// therefore resolves every variation, reference and signal-extraction file exactly like its full-ring
// sibling, which is why the longitudinal rings enter this macro as extra families and not as a separate
// code path. Files are built from these tags, never matched against them, so "_useRingZ" being a prefix
// of "_useRingZT" is harmless here (see IsLongitudinalSuffix for the one place that matches).
const std::vector<LongitudinalRing> LONGITUDINAL_RINGS = {
    {"RingZ",  "_useRingZ",  "_{z}",     "Rz"},   // R_z = P_z n_z
    {"RingZT", "_useRingZT", "_{z}^{T}", "RzT"}}; // R_z^T = R_z minus its helicity term

// Display style of each longitudinal definition, shared by every cross-definition canvas so a definition looks the
// same everywhere. Indexed like LONGITUDINAL_RINGS (modulo 2); the full ring is always kGray + 2 with marker 24.
const Color_t LONG_RING_COLOR[2] = {kRed + 1, kGreen + 2};
const Style_t LONG_RING_MARKER[2] = {21, 22};      // Filled: the headline quantity of the definition
const Style_t LONG_RING_MARKER_OPEN[2] = {25, 26}; // Open counterpart, for the raw mean beside it

/// @brief True when a consumer file suffix is a longitudinal one: exact match on its "_"-separated tokens.
bool IsLongitudinalSuffix(const std::string& suffix) {
    std::stringstream tokens(suffix);
    for (std::string token; std::getline(tokens, token, '_');)
        for (const auto& lr : LONGITUDINAL_RINGS)
            if ("_" + token == lr.fileTag) return true;
    return false;
}

// Decay constants, as in the consumer. Only alpha^2 enters here (the exact full-ring inversion).
constexpr double LambdaWeakDecayConstant = 0.749; // PDG 2025 update, same as DerivedDataConsumer's
constexpr double AntiLambdaWeakDecayConstant = -0.758; // PDG 2025 update, same as DerivedDataConsumer's
constexpr double AlphaSqLambda = LambdaWeakDecayConstant * LambdaWeakDecayConstant;
constexpr double AlphaSqAntiLambda = AntiLambdaWeakDecayConstant * AntiLambdaWeakDecayConstant;
// BothHyperons mixes the two, and the exact term is per species. The plain mean is within 1.2% of either
// alpha^2, and the term it multiplies is itself a ~10% correction, so the choice moves P by <~0.2%.
constexpr double AlphaSqBoth = 0.5 * (AlphaSqLambda + AlphaSqAntiLambda);

/**
 * @brief Swaps the full-ring symbol for a longitudinal one (suffix sfx, e.g. "_{z}") in one piece of display text.
 *
 * The titles here are built from many literals ("<R>", "Integrated <R>_{S}", a bare "R" y title, and
 * "#Delta" + that title), so longitudinal families are relabelled in one pass at the end rather than at each of them.
 * Only whole-token forms are touched: a jet-radius "R = 0.4" is left alone.
 */
std::string SwapRingText(std::string t, const std::string& sfx) {
    auto swapAll = [&t](const std::string& from, const std::string& to) {
        for (size_t pos = t.find(from); pos != std::string::npos; pos = t.find(from, pos + to.size()))
            t.replace(pos, from.size(), to);
    };
    if (t == "R") return "R" + sfx;
    if (t.rfind("#DeltaR ", 0) == 0) t.replace(0, 8, "#DeltaR" + sfx + " ");
    swapAll("<R>", "<R" + sfx + ">");
    swapAll("#it{R}", "#it{R}" + sfx);
    return t;
}

/// @brief Relabels every text an object carries: titles, axis titles, bin labels, legend and pave lines.
void RelabelRingObject(TObject* obj, const std::string& sfx) {
    if (!obj) return;
    if (auto* pad = dynamic_cast<TPad*>(obj)) { // TCanvas included; recurse into everything drawn on it
        pad->SetTitle(SwapRingText(pad->GetTitle(), sfx).c_str());
        TIter next(pad->GetListOfPrimitives());
        while (TObject* prim = next()) RelabelRingObject(prim, sfx);
    } else if (auto* h = dynamic_cast<TH1*>(obj)) {
        h->SetTitle(SwapRingText(h->GetTitle(), sfx).c_str());
        for (TAxis* ax : {h->GetXaxis(), h->GetYaxis(), h->GetZaxis()}) {
            ax->SetTitle(SwapRingText(ax->GetTitle(), sfx).c_str());
            if (THashList* labels = ax->GetLabels()) {
                TIter nextLabel(labels);
                while (auto* lab = static_cast<TObjString*>(nextLabel())) lab->SetString(SwapRingText(lab->GetString().Data(), sfx).c_str());
            }
        }
    } else if (auto* g = dynamic_cast<TGraph*>(obj)) {
        g->SetTitle(SwapRingText(g->GetTitle(), sfx).c_str());
        if (TH1* frame = g->GetHistogram()) RelabelRingObject(frame, sfx);
    } else if (auto* leg = dynamic_cast<TLegend*>(obj)) {
        TIter nextEntry(leg->GetListOfPrimitives());
        while (auto* e = static_cast<TLegendEntry*>(nextEntry())) e->SetLabel(SwapRingText(e->GetLabel(), sfx).c_str());
    } else if (auto* pave = dynamic_cast<TPaveText*>(obj)) {
        TIter nextLine(pave->GetListOfLines());
        while (auto* line = dynamic_cast<TText*>(nextLine())) line->SetTitle(SwapRingText(line->GetTitle(), sfx).c_str());
    } else if (auto* txt = dynamic_cast<TText*>(obj)) { // TLatex included
        txt->SetTitle(SwapRingText(txt->GetTitle(), sfx).c_str());
    }
}

/**
 * @brief Relabels every canvas, histogram and graph under a directory, recursively, rewriting each in place.
 *
 * Key names are collected first: writing while iterating would add cycles to the very list being walked.
 */
void RelabelRingDirectory(TDirectory* dir, const std::string& sfx) {
    if (!dir) return;
    std::vector<std::pair<std::string, TClass*>> items;
    std::map<std::string, bool> seen; // Newest cycle only
    TIter nextKey(dir->GetListOfKeys());
    while (auto* key = static_cast<TKey*>(nextKey())) {
        if (seen[key->GetName()]) continue;
        seen[key->GetName()] = true;
        items.push_back({key->GetName(), TClass::GetClass(key->GetClassName())});
    }
    for (const auto& item : items) {
        const std::string& name = item.first;
        TClass* cls = item.second;
        if (!cls) continue;
        if (cls->InheritsFrom(TDirectory::Class())) { RelabelRingDirectory(dir->GetDirectory(name.c_str()), sfx); continue; }
        if (!cls->InheritsFrom(TPad::Class()) && !cls->InheritsFrom(TH1::Class()) && !cls->InheritsFrom(TGraph::Class())) continue;
        TObject* obj = dir->Get(name.c_str());
        RelabelRingObject(obj, sfx);
        dir->cd();
        obj->Write(name.c_str(), TObject::kOverwrite);
        if (auto* h = dynamic_cast<TH1*>(obj)) h->SetDirectory(nullptr);
        delete obj;
    }
}

/// @brief Sum_i y_i and Sum_i y_i^2 over every candidate of a profile, underflow/overflow bins included.
struct CandidateMoments {
    double n = 0., sum = 0., sum2 = 0.;
    double Mean() const { return n > 0. ? sum / n : 0.; }
    double Var() const { return n > 0. ? sum2 / n - Mean() * Mean() : 0.; } // Per-candidate variance
};
// TODO: revise all the instances of CandidateMoments and check if it is indeed best to have overflow/underflow bins included or not.
// It seemed like a good idea at the designing stage, but there may be some edge case where this is not the best idea.

/// @brief Rebuilds the candidate moments from the per-cell means and spreads ("s" error option).
template <typename P>
CandidateMoments MomentsOfProfile(P* p) {
    CandidateMoments m;
    const TString opt = p->GetErrorOption();
    p->SetErrorOption("s");
    for (int b = 0; b < p->GetNcells(); ++b) {
        const double n = p->GetBinEntries(b);
        if (n <= 0.) continue;
        const double mean = p->GetBinContent(b), s = p->GetBinError(b);
        m.n += n;
        m.sum += n * mean;
        m.sum2 += n * (s * s + mean * mean);
    }
    p->SetErrorOption(opt);
    return m;
}

/// @brief A simplified getter that wraps the call for either a TProfile2D or 1D.
CandidateMoments MomentsOf(TH1* h) {
    if (auto* p2 = dynamic_cast<TProfile2D*>(h)) return MomentsOfProfile(p2);
    if (auto* p1 = dynamic_cast<TProfile*>(h)) return MomentsOfProfile(p1);
    return {};
}

/**
 * @brief Candidate moments of one consumer file, for one proxy: the ring and the KappaEff profiles.
 *
 * Every profile is mass-binned and filled on exactly the proxy's candidates, and every bin and flow is
 * summed, so all five share one candidate set. The ring comes from each proxy's 1D pRingObservableMass:
 * the (eta_proxy, mass) maps would silently lose every candidate outside their eta axis. Nothing here is
 * signal-extracted: these are all candidates, signal and background mixed.
 */
struct RingSample {
    CandidateMoments ring;    // R (or R_z, R_z^T), with the prefactor
    CandidateMoments num;     // u v: u = R / prefactor, v = the unitless full ring (u^2 for the full ring)
    CandidateMoments den;     // w: 1 (full ring) or n_z^2 (R_z, R_z^T)
    CandidateMoments numDen;  // u v w
    CandidateMoments ringDen; // R w
};

bool LoadRingSample(FileCache& cache, const std::string& file, const std::string& proxy, RingSample& s) {
    if (!cache.Open(file, false)) return false; // Absent config: only counted, it is the common case
    // Every proxy folder holds the same object names, so the proxy only picks the folder
    const std::string base = TASK_DIR_IN_FILE + RING_DIR_IN_FILE + proxy + "/";
    const std::string kap  = base + "KappaEff/"; // Each proxy carries its own KappaEff/

    TH1* h[5] = {cache.FetchClone<TH1>(file, base + "MassDependence/pRingObservableMass", false),
                 cache.FetchClone<TH1>(file, kap + "pKappaNumVsMass", false),
                 cache.FetchClone<TH1>(file, kap + "pKappaDenVsMass", false),
                 cache.FetchClone<TH1>(file, kap + "pKappaNumTimesDenVsMass", false),
                 cache.FetchClone<TH1>(file, kap + "pRingTimesDenVsMass", false)};
    const bool ok = h[0] && h[1] && h[2] && h[3] && h[4];
    if (ok) {
        s.ring = MomentsOf(h[0]); s.num = MomentsOf(h[1]); s.den = MomentsOf(h[2]);
        s.numDen = MomentsOf(h[3]); s.ringDen = MomentsOf(h[4]);
    }
    for (TH1* p : h) delete p;
    if (!ok) return false; // Each absent object is already counted by FetchClone
    if (s.ring.n < 2.) {
        WarnSkip(proxy + " ring sample", "fewer than two candidates", file);
        return false;
    }
    if (std::fabs(s.ring.n - s.num.n) > 0.5) {
        WarnSkip(proxy + " ring sample", "ring and KappaEff candidate counts differ", file);
        return false;
    }
    return true;
}

/// @brief <R> of a sample with its standard error.
std::pair<double, double> MeanOf(const CandidateMoments& m) {
    return {m.Mean(), std::sqrt(std::max(m.Var(), 0.) / m.n)};
}

/**
 * @brief <a>/<b> over one candidate set, with the first-order error that keeps their covariance:
 *   Var(r) = [ s_a^2 - 2 r C + r^2 s_b^2 ] / (n <b>^2),   r = <a>/<b>,   C = <ab> - <a><b>.
 * For the full 3D ring, b = 1 (the full ring's "w" factor), and this becomes the plain standard error of <a>.
 */
std::pair<double, double> RatioOfMeans(const CandidateMoments& a, const CandidateMoments& b, double meanAB) {
    const double bm = b.Mean();
    if (bm <= 0. || a.n < 2.) return {0., -1.};
    const double r = a.Mean() / bm;
    const double cov = meanAB - a.Mean() * bm;
    const double var = (a.Var() - 2. * r * cov + r * r * b.Var()) / (a.n * bm * bm);
    return {r, std::sqrt(std::max(var, 0.))};
}

/// @brief kappa = 3 <u v>/<w>, the cross moment: the same expression for every ring definition.
std::pair<double, double> KappaOf(const RingSample& s) {
    const auto r = RatioOfMeans(s.num, s.den, s.numDen.Mean());
    return {3. * r.first, 3. * r.second};
}

/**
 * @brief Several series on one categorical axis, one bin per column, the series nudged sideways within it.
 * A thin adapter onto DrawOffsetSeriesCanvas: it only packs vals[series][column] into one histogram per
 * series. A value with a negative error is treated as missing and left empty, which the offset canvas
 * skips (no marker parked at zero). Its columns are sample labels, which are short, so the frame reaches
 * much lower than the offset canvas default.
 */
void DrawSeriesByColumn(const std::vector<std::string>& columnLabels, const std::vector<VariationConfig>& series,
                        const std::vector<std::vector<std::pair<double, double>>>& vals,
                        const std::string& canvasName, const std::string& title, TDirectory* dir, const std::string& yTitle) {
    const int nCols = static_cast<int>(columnLabels.size());
    if (!dir || nCols == 0 || series.empty()) return;

    std::vector<ProfileBundle> bundles;
    for (size_t sIdx = 0; sIdx < series.size(); ++sIdx) {
        TH1D* h = new TH1D(Form("%s_series%zu", canvasName.c_str(), sIdx), "", nCols, 0, nCols);
        h->SetDirectory(nullptr);
        for (int c = 0; c < nCols; ++c) {
            if (vals[sIdx][c].second < 0.) continue;
            h->SetBinContent(c + 1, vals[sIdx][c].first);
            h->SetBinError(c + 1, vals[sIdx][c].second);
        }
        bundles.push_back({h, series[sIdx]});
    }
    DrawOffsetSeriesCanvas(columnLabels, bundles, canvasName, title, dir, yTitle, 0.5, false, 0.10);
    for (auto& b : bundles) delete b.profile;
}

/**
 * @brief CheapSigExtract inputs of one consumer file, for one proxy: the four *VsMassRegion profiles.
 *
 * The consumer fills them from the same two flags as IntegratedCuts/<Proxy>/p2dRingCutsV0MassPeak:
 * v0InMassWindow (bin 1, the sideband) and v0InMassPeak (bin 2). Index 0 is the sideband, 1 the peak.
 */
struct RegionSample {
    CandidateMoments ring[2], num[2], den[2], numDen[2];
    CandidateMoments ringDen[2]; // R w: only in consumer outputs that book it (empty, n = 0, in older ones)
};

/// @brief The candidate moments of one bin of a TProfile, from its mean and spread ("s" error option).
CandidateMoments BinMoments(TProfile* p, int bin) {
    CandidateMoments m;
    const double n = p->GetBinEntries(bin);
    if (n <= 0.) return m;
    const TString opt = p->GetErrorOption();
    p->SetErrorOption("s");
    const double mean = p->GetBinContent(bin), sd = p->GetBinError(bin);
    p->SetErrorOption(opt);
    m.n = n;
    m.sum = n * mean;
    m.sum2 = n * (sd * sd + mean * mean);
    return m;
}

bool LoadRegionSample(FileCache& cache, const std::string& file, const std::string& proxy, RegionSample& r) {
    if (!cache.Open(file, false)) return false; // Absent config: only counted
    // Every proxy folder holds the same object names, so the proxy only picks the folder
    const std::string kap = TASK_DIR_IN_FILE + RING_DIR_IN_FILE + proxy + "/KappaEff/";
    TProfile* p[4] = {cache.FetchClone<TProfile>(file, kap + "pRingVsMassRegion", false),
                      cache.FetchClone<TProfile>(file, kap + "pKappaNumVsMassRegion", false),
                      cache.FetchClone<TProfile>(file, kap + "pKappaDenVsMassRegion", false),
                      cache.FetchClone<TProfile>(file, kap + "pKappaNumTimesDenVsMassRegion", false)};
    const bool ok = p[0] && p[1] && p[2] && p[3];
    if (ok)
        for (int i = 0; i < 2; ++i) {
            r.ring[i] = BinMoments(p[0], i + 1); r.num[i] = BinMoments(p[1], i + 1);
            r.den[i] = BinMoments(p[2], i + 1);  r.numDen[i] = BinMoments(p[3], i + 1);
        }
    for (TProfile* q : p) delete q;
    // <R w> only serves the <R_long>/<n_z^2> ratio, so its absence costs that ratio alone (left as a gap), never the sample:
    if (ok)
        if (TProfile* rw = cache.FetchClone<TProfile>(file, kap + "pRingTimesDenVsMassRegion", false)) {
            for (int i = 0; i < 2; ++i) r.ringDen[i] = BinMoments(rw, i + 1);
            delete rw;
        }
    if (!ok) return false; // Each absent object is already counted by FetchClone
    if (!(r.ring[1].n > r.ring[0].n)) { // A peak no richer than its sideband has no signal to estimate
        WarnSkip(proxy + " mass-region sample", "the peak holds no more candidates than the sideband", file);
        return false;
    }
    return true;
}

/**
 * @brief The cheap signal mean of y: <y>_S = (Sum_P y - Sum_B y) / (N_P - N_B).
 *
 * Single sideband of the same total width as the peak window (the consumer enforces it), so for a background
 * whose density and <y>_B are linear in mass, the sideband holds exactly the background under the peak.
 * A linear density times a linear <y>_B is quadratic, and that curvature is what this misses; splitting the
 * sideband into left and right would recover it (see the README).
 *
 * The two regions are disjoint, so they are independent. Treating each count as Poisson (a compound sum):
 *   Var(A) = N_P <y^2>_P + N_B <y^2>_B,  Var(D) = N_P + N_B,  Cov(A, D) = Sum_P y + Sum_B y,
 *   Var(<y>_S) = [Var(A) - 2 r Cov(A, D) + r^2 Var(D)] / D^2,   A = Sum_P - Sum_B,  D = N_P - N_B,  r = A/D.
 * With no background it reduces to the standard error of <y>_P; for y == 1 it is exactly 1 with no error.
 */
std::pair<double, double> CheapSignalMean(const CandidateMoments& peak, const CandidateMoments& side) {
    const double D = peak.n - side.n;
    if (!(D > 0.)) return {0., -1.};
    const double A = peak.sum - side.sum;
    const double r = A / D;
    const double varA = peak.sum2 + side.sum2;
    const double varD = peak.n + side.n;
    const double cov  = peak.sum + side.sum;
    const double var  = (varA - 2. * r * cov + r * r * varD) / (D * D);
    return {r, std::sqrt(std::max(var, 0.))};
}

/// @brief Adds one candidate set's moments to another's (disjoint sets: n, Sum y and Sum y^2 all add up).
void AddMoments(CandidateMoments& to, const CandidateMoments& m) {
    to.n += m.n; to.sum += m.sum; to.sum2 += m.sum2;
}

/// @brief The candidate moments of one cell of a TProfile2D (global bin), from its mean and spread ("s" error option).
CandidateMoments CellMoments(TProfile2D* p, int globalBin) {
    CandidateMoments m;
    const double n = p->GetBinEntries(globalBin);
    if (n <= 0.) return m;
    const TString opt = p->GetErrorOption();
    p->SetErrorOption("s");
    const double mean = p->GetBinContent(globalBin), sd = p->GetBinError(globalBin);
    p->SetErrorOption(opt);
    m.n = n;
    m.sum = n * mean;
    m.sum2 = n * (sd * sd + mean * mean);
    return m;
}

/**
 * @brief Per-bin candidate moments of one proxy, vs proxy pT or centrality (Ring/<Proxy>/KappaEff/Differential/).
 *
 * Per x bin: all candidates (every mass-region cell, underflow included, as in the integrated Raw) and the two
 * CheapSigExtract regions. Binning and title are read from the consumer's profiles, which share the event mixing
 * axes, so each bin holds the same Lambdas in data and in its mixed events.
 */
struct DifferentialSample {
    std::vector<double> edges;        // x bin edges, nBins + 1
    std::string xTitle;
    std::vector<RingSample> all;      // One per x bin
    std::vector<RegionSample> region; // One per x bin
};

/**
 * @brief Loads the differential moments of one file, proxy and variable ("ProxyPt" or "Centrality").
 *
 * Closure: summed over every cell, x flows included, they must hold exactly the integrated candidates (same fill
 * list, same gate). A mismatch is reported, not fatal: the five differential profiles stay consistent among
 * themselves, but they no longer describe the integrated sample.
 */
bool LoadDifferentialSample(FileCache& cache, const std::string& file, const std::string& proxy,
                            const std::string& var, DifferentialSample& d) {
    if (!cache.Open(file, false)) return false; // Absent config: only counted
    static const char* kMoment[5] = {"Ring", "KappaNum", "KappaDen", "KappaNumTimesDen", "RingTimesDen"};
    // Every proxy folder holds the same object names, so the proxy only picks the folder
    const std::string dir = TASK_DIR_IN_FILE + RING_DIR_IN_FILE + proxy + "/KappaEff/Differential/p2d";
    TProfile2D* p[5] = {nullptr, nullptr, nullptr, nullptr, nullptr};
    bool ok = true;
    for (int k = 0; k < 5; ++k) {
        p[k] = cache.FetchClone<TProfile2D>(file, dir + kMoment[k] + "Vs" + var + "VsMassRegion", false);
        ok = ok && p[k];
    }
    if (ok) {
        const TAxis* ax = p[0]->GetXaxis();
        const int nx = ax->GetNbins(), ny = p[0]->GetNbinsY();
        d.xTitle = ax->GetTitle();
        d.edges.clear();
        for (int ix = 1; ix <= nx + 1; ++ix) d.edges.push_back(ax->GetBinLowEdge(ix));
        d.all.assign(nx, RingSample{});
        d.region.assign(nx, RegionSample{});
        double nTotal = 0.;
        for (int ix = 0; ix <= nx + 1; ++ix)
            for (int iy = 0; iy <= ny + 1; ++iy) {
                const int g = p[0]->GetBin(ix, iy);
                CandidateMoments m[5];
                for (int k = 0; k < 5; ++k) m[k] = CellMoments(p[k], g);
                nTotal += m[0].n;
                if (ix < 1 || ix > nx) continue; // x flows only enter the closure
                RingSample& a = d.all[ix - 1];
                AddMoments(a.ring, m[0]); AddMoments(a.num, m[1]); AddMoments(a.den, m[2]);
                AddMoments(a.numDen, m[3]); AddMoments(a.ringDen, m[4]);
                if (iy == 1 || iy == 2) { // Mass-region bins: 1 sideband (index 0), 2 peak (index 1)
                    RegionSample& r = d.region[ix - 1];
                    const int i = iy - 1;
                    r.ring[i] = m[0]; r.num[i] = m[1]; r.den[i] = m[2]; r.numDen[i] = m[3]; r.ringDen[i] = m[4];
                }
            }
        static std::set<std::string> checked; // Every folder reloads the same files: check each one once
        RingSample integrated;
        if (checked.insert(file + "|" + proxy + "|" + var).second && LoadRingSample(cache, file, proxy, integrated) &&
            std::fabs(integrated.ring.n - nTotal) > 0.5)
            std::cerr << "  [Warning] " << proxy << " vs " << var << ": " << nTotal << " differential candidates against "
                      << integrated.ring.n << " integrated in " << file << std::endl;
    }
    for (TProfile2D* q : p) delete q;
    return ok;
}

/**
 * @brief The cheap signal ratio <a>_S/<b>_S = A_a/A_b, with A_y = Sum_P y - Sum_B y (the N_P - N_B cancels).
 *
 * Same compound-Poisson treatment as CheapSignalMean, with the a-b covariance from Sum ab of each region:
 *   Var(A_y) = Sum_P y^2 + Sum_B y^2,   Cov(A_a, A_b) = Sum_P ab + Sum_B ab,
 *   Var(rho) = [Var(A_a) - 2 rho Cov(A_a, A_b) + rho^2 Var(A_b)] / A_b^2,   rho = A_a/A_b.
 * This is exact to first order, so kappa_S = 3 rho(u v, w) no longer needs its two errors in quadrature.
 *
 * @param abPeak, abSide The moments of the product a b in each region (only their sums are read).
 */
std::pair<double, double> CheapSignalRatio(const CandidateMoments& aPeak, const CandidateMoments& aSide,
                                           const CandidateMoments& bPeak, const CandidateMoments& bSide,
                                           const CandidateMoments& abPeak, const CandidateMoments& abSide) {
    const double Ab = bPeak.sum - bSide.sum;
    if (!(Ab > 0.) || abPeak.n <= 0.) return {0., -1.};
    const double rho  = (aPeak.sum - aSide.sum) / Ab;
    const double varA = aPeak.sum2 + aSide.sum2;
    const double varB = bPeak.sum2 + bSide.sum2;
    const double cov  = abPeak.sum + abSide.sum;
    const double var  = (varA - 2. * rho * cov + rho * rho * varB) / (Ab * Ab);
    return {rho, std::sqrt(std::max(var, 0.))};
}

/// @brief One column of a full-ring vs longitudinal comparison: a config, read from both definitions.
struct RingPairColumn { VariationConfig cfg; std::string base, suffix; };

/**
 * @brief Data, every in-wagon variation and every consumer-style reference. The file of a column is
 *        <base>/ConsumerResults_<family><suffix>.root, and its longitudinal siblings the same with their tag after
 *        the family; a column is drawn only where both exist.
 */
std::vector<RingPairColumn> MakeRingPairColumns(const std::string& consumerDir, const VariationConfig& dataConfig,
                                                const std::vector<VariationConfig>& sysVariations,
                                                const std::vector<ExternalRef>& externals) {
    std::vector<RingPairColumn> columns = {{dataConfig, consumerDir, ""}};
    for (const auto& sys : sysVariations) columns.push_back({sys, consumerDir, sys.suffix});
    for (const auto& ext : externals)
        if (!ext.basePath.empty() && ext.kind == ExternalRef::kConsumerDir) columns.push_back({ext.config, ext.basePath, ""});
    return columns;
}

/**
 * @brief The series of every Full 3D Ring vs longitudinal canvas: <R> of the full ring, then <R_long> and
 *        <R_long>/<n_z^2> for each of LONGITUDINAL_RINGS. Series 1 + 2l and 2 + 2l belong to definition l.
 *
 * One colour per definition, as in Corrections/; the raw mean is open and the per-candidate amplitude filled,
 * the latter being the one emphasised.
 */
std::vector<VariationConfig> MakeRingVsLongitudinalSeries() {
    std::vector<VariationConfig> series = {{"", "<R> (Full 3D Ring)", kGray + 2, 1, 24, false}};
    for (size_t l = 0; l < LONGITUDINAL_RINGS.size(); ++l) {
        const std::string sym = "R" + LONGITUDINAL_RINGS[l].symbolSuffix;
        series.push_back({"", "<" + sym + ">", LONG_RING_COLOR[l % 2], 1, LONG_RING_MARKER_OPEN[l % 2], false});
        series.push_back({"", "<" + sym + ">/<n_{z}^{2}>", LONG_RING_COLOR[l % 2], 1, LONG_RING_MARKER[l % 2], true});
    }
    return series;
}

/**
 * @brief The headline check: <R> of the Full 3D Ring against <R_long> and <R_long>/<n_z^2> of every longitudinal
 *        definition, per family and proxy, across configs, all on one canvas.
 *
 * Columns are Data, then every in-wagon variation and consumer-style external reference that has its full-ring
 * file and at least one longitudinal one; a definition whose file is absent is left as a gap. If the longitudinal
 * rings do what they should, they are ~0 for MixedEv and survive only in data, while the full ring carries its
 * large fake everywhere. <R_long>/<n_z^2> is the per-candidate ring amplitude, comparable with <R> once the
 * latter's fake is removed (under ring alignment). RAW values: no signal extraction and no kappa (for those, see
 * Corrections/). Nothing is written when no such file exists. Each definition also gets its <R_long>/<n_z^2>
 * alone, with its Data - Var differences and redux subset, over the columns where that definition exists.
 * Written to RingVsLongitudinal/<Family>/<Proxy>/.
 *
 * @param families Untagged families: the longitudinal files of each are the same name plus each fileTag.
 */
void DrawRingVsLongitudinal(FileCache& cache, TDirectory* fileRoot, const std::string& consumerDir,
                            const std::vector<FamilyConfig>& families, const VariationConfig& dataConfig,
                            const std::vector<VariationConfig>& sysVariations, const std::vector<ExternalRef>& externals) {
    const std::vector<VariationConfig> series = MakeRingVsLongitudinalSeries();
    const std::pair<double, double> kMissing = {0., -1.};
    const std::vector<RingPairColumn> columns = MakeRingPairColumns(consumerDir, dataConfig, sysVariations, externals);

    for (const auto& fam : families) {
        for (const auto& proxy : PROXIES) {
            std::vector<std::string> labels;
            std::vector<VariationConfig> cfgs;
            std::vector<std::vector<std::pair<double, double>>> vals(series.size()); // [series][column]
            bool dataFirst = false;
            for (const auto& col : columns) {
                RingSample full;
                const std::string stem = col.base + "/ConsumerResults_" + fam.dataSuffix;
                if (!LoadRingSample(cache, stem + col.suffix + ".root", proxy, full)) continue;

                std::vector<std::pair<double, double>> row(series.size(), kMissing); // This column, every series
                row[0] = MeanOf(full.ring);
                bool anyLong = false;
                for (size_t l = 0; l < LONGITUDINAL_RINGS.size(); ++l) {
                    RingSample lng;
                    if (!LoadRingSample(cache, stem + LONGITUDINAL_RINGS[l].fileTag + col.suffix + ".root", proxy, lng))
                        continue;
                    row[1 + 2 * l] = MeanOf(lng.ring);
                    row[2 + 2 * l] = RatioOfMeans(lng.ring, lng.den, lng.ringDen.Mean());
                    anyLong = true;
                }
                if (!anyLong) continue;

                if (labels.empty()) dataFirst = col.cfg.isData;
                labels.push_back(col.cfg.legendLabel);
                cfgs.push_back(col.cfg);
                for (size_t s = 0; s < series.size(); ++s) vals[s].push_back(row[s]);
                printf("    [RingVsLongitudinal %s, %s, %s] <R> = %+.3e +/- %.1e", fam.familyName.c_str(), proxy.c_str(),
                       col.cfg.legendLabel.c_str(), row[0].first, row[0].second);
                for (size_t l = 0; l < LONGITUDINAL_RINGS.size(); ++l) {
                    const char* sn = LONGITUDINAL_RINGS[l].shortName.c_str();
                    if (row[1 + 2 * l].second < 0.) { printf("   <%s>: absent", sn); continue; }
                    printf("   <%s> = %+.3e +/- %.1e   <%s>/<n_z^2> = %+.3e +/- %.1e", sn, row[1 + 2 * l].first,
                           row[1 + 2 * l].second, sn, row[2 + 2 * l].first, row[2 + 2 * l].second);
                }
                printf("\n");
            }
            if (!dataFirst) { // Differences need Data as the first column
                CountSkip("RingVsLongitudinal", labels.empty() ? "no column holds both ring definitions"
                                                               : "no Data column, which the differences need",
                          fam.familyName + "/" + proxy);
                continue;
            }
            TDirectory* pxDir = EnsureDir(EnsureDir(EnsureDir(fileRoot, "RingVsLongitudinal"), fam.familyName), proxy);

            DrawSeriesByColumn(labels, series, vals, "Canvas_RingVsLongitudinal_" + proxy,
                               fam.familyName + ", " + proxy + ": Full 3D Ring vs longitudinal rings (raw, all candidates)",
                               pxDir, "integrated value");

            // <R_long>/<n_z^2> of each definition alone. The categorical canvas has no notion of a gap, so it gets
            // only the columns where this definition has a value, and still needs Data first for its differences.
            for (size_t l = 0; l < LONGITUDINAL_RINGS.size(); ++l) {
                const auto& lr = LONGITUDINAL_RINGS[l];
                const auto& ratio = vals[2 + 2 * l];
                std::vector<std::string> subLabels;
                std::vector<VariationConfig> subCfgs;
                std::vector<std::pair<double, double>> subVals;
                for (size_t c = 0; c < labels.size(); ++c) {
                    if (ratio[c].second < 0.) continue;
                    subLabels.push_back(labels[c]);
                    subCfgs.push_back(cfgs[c]);
                    subVals.push_back(ratio[c]);
                }
                if (subCfgs.empty() || !subCfgs.front().isData) {
                    CountSkip("RingVsLongitudinal, <R_long>/<n_z^2> alone", "no Data column for this definition",
                              fam.familyName + "/" + proxy + "/" + lr.name);
                    continue;
                }
                const std::string sym = "R" + lr.symbolSuffix;
                DrawCategoricalCrossSystem(subLabels, subCfgs, subVals, "Canvas_" + lr.shortName + "OverNz2_" + proxy,
                                           fam.familyName + ", " + proxy + ": <" + sym + ">/<n_{z}^{2}> (raw)", pxDir,
                                           "<" + sym + ">/<n_{z}^{2}>");
            }
        }
    }
}

/**
 * @brief CheapSigExtract of the same comparison: <R> and every <R_long>, in and out of the mass peak, and how far
 *        the background pulls each (in - out).
 *
 * Same columns and series as RingVsLongitudinal/, read from IntegratedCuts/<Proxy>/p2dRingCutsV0MassPeak
 * (FetchIntegratedByMassPeak): no fit and no sideband model. In-peak and out-of-peak candidates are disjoint, so
 * in - out is plain quadrature. The peak and sideband canvases add <R_long>/<n_z^2> of that region, from the
 * consumer's Ring/<Proxy>/KappaEff/p*VsMassRegion profiles: ring, n_z^2 and their product on the same candidates,
 * so the ratio keeps its covariance. Consumer outputs without pRingTimesDenVsMassRegion leave it as a gap.
 * Written to RingVsLongitudinal_CheapSigExtract/<Family>/<Proxy>/.
 */
void DrawCheapRingVsLongitudinal(FileCache& cache, TDirectory* fileRoot, const std::string& consumerDir,
                                 const std::vector<FamilyConfig>& families, const VariationConfig& dataConfig,
                                 const std::vector<VariationConfig>& sysVariations, const std::vector<ExternalRef>& externals) {
    const std::vector<VariationConfig> series = MakeRingVsLongitudinalSeries();
    const std::pair<double, double> kMissing = {0., -1.};
    const std::vector<RingPairColumn> columns = MakeRingPairColumns(consumerDir, dataConfig, sysVariations, externals);
    using StateSeries = std::vector<std::vector<std::pair<double, double>>>; // [state][series], state 0 in, 1 out, 2 in - out

    // Fills the three states of one series from a file's in-peak and out-of-peak points; false if either is absent
    auto fillInOut = [&](StateSeries& row, size_t s, const std::string& file, const std::string& proxy) {
        const SigExtractPoint in = FetchIntegratedByMassPeak(cache, file, proxy, true);
        const SigExtractPoint out = FetchIntegratedByMassPeak(cache, file, proxy, false);
        if (!in.found || !out.found) return false;
        row[0][s] = {in.value, in.error};
        row[1][s] = {out.value, out.error};
        row[2][s] = {in.value - out.value, std::hypot(in.error, out.error)};
        return true;
    };

    for (const auto& fam : families) {
        for (const auto& proxy : PROXIES) {
            std::vector<std::string> labels;
            std::vector<StateSeries> vals(3, StateSeries(series.size())); // [state][series][column]
            for (const auto& col : columns) {
                const std::string stem = col.base + "/ConsumerResults_" + fam.dataSuffix;
                const std::string fullFile = stem + col.suffix + ".root";
                if (!cache.Open(fullFile, false)) continue; // Absent config: only counted

                StateSeries row(3, std::vector<std::pair<double, double>>(series.size(), kMissing)); // This column
                if (!fillInOut(row, 0, fullFile, proxy)) continue;
                bool anyLong = false;
                for (size_t l = 0; l < LONGITUDINAL_RINGS.size(); ++l) {
                    const std::string lngFile = stem + LONGITUDINAL_RINGS[l].fileTag + col.suffix + ".root";
                    if (!cache.Open(lngFile, false)) continue; // Absent config: only counted
                    if (!fillInOut(row, 1 + 2 * l, lngFile, proxy)) continue;
                    anyLong = true;
                    // <R_long>/<n_z^2> per region (index 1 peak, 0 sideband); in - out of a ratio is not drawn
                    RegionSample region;
                    if (!LoadRegionSample(cache, lngFile, proxy, region)) continue;
                    for (int r = 0; r < 2; ++r) {
                        const int idx = (r == 0) ? 1 : 0; // State 0 is the peak (region index 1), state 1 the sideband
                        if (region.ringDen[idx].n > 0.)
                            row[r][2 + 2 * l] = RatioOfMeans(region.ring[idx], region.den[idx], region.ringDen[idx].Mean());
                    }
                }
                if (!anyLong) continue;

                labels.push_back(col.cfg.legendLabel);
                for (int st = 0; st < 3; ++st)
                    for (size_t s = 0; s < series.size(); ++s) vals[st][s].push_back(row[st][s]);
            }
            if (labels.empty()) {
                CountSkip("RingVsLongitudinal_CheapSigExtract", "no column holds both ring definitions", fam.familyName + "/" + proxy);
                continue;
            }
            TDirectory* pxDir = EnsureDir(EnsureDir(EnsureDir(fileRoot, "RingVsLongitudinal_CheapSigExtract"), fam.familyName), proxy);

            const std::string head = fam.familyName + ", " + proxy + ": Full 3D Ring vs longitudinal rings, ";
            const std::string stemName = "Canvas_RingVsLongitudinal";
            DrawSeriesByColumn(labels, series, vals[0], stemName + "_InMassPeak_" + proxy,
                               head + "in mass peak", pxDir, "integrated value");
            DrawSeriesByColumn(labels, series, vals[1], stemName + "_OutOfMassPeak_" + proxy,
                               head + "out of mass peak", pxDir, "integrated value");
            DrawSeriesByColumn(labels, series, vals[2], stemName + "_InMinusOut_" + proxy,
                               head + "in - out of mass peak", pxDir, "#Delta (in - out)");
        }
    }
}

/**
 * @brief Corrections: the fake subtracted and the response divided out, with or without signal extraction.
 *
 * Both ring definitions are reduced to the same quantity, the ring density P_R (under ring alignment):
 *   full ring  P = (R_d - R_f) / (kappa_eff - (alpha^2/3) R_d R_f)    exact inversion (eq:finalcorrection)
 *   R_z, R_z^T P = (R_z,d - R_z,f) / (3 <u v>) = R_z,true / <n_z^2>   linear, see below
 * with d the data sample, f its mixed-event partner, and kappa and <u v> measured on the data sample.
 * The exact R_z inversion is R_z,true = [R_z,d (1 + <delta>_0) - R_z,f] / kappa_z, where <delta>_0 ~
 * (alpha^2/3) R_true R_fake is the normalization of the WHOLE sample -- full-ring quantities, whatever R_z's
 * own fake is. It multiplies R_z,d, so the dropped term is <delta>_0 R_z,d ~ 1e-4 x R_z,d: negligible even
 * with R_z's fake as large as its signal. The full ring cannot drop it, because there it multiplies R_d ~ 0.8.
 * For
 * the full ring kappa_eff = 3 <u v> = 3 <u^2>, so all are the one expression P = (R_d - R_f) / (3 <u v> - k R_d R_f),
 * with k = alpha^2/3 for the full ring (FamilyConfig::alphaSq) and k = 0 for R_z and R_z^T.
 *
 * Errors are first order, with the full derivative. The data and fake rings are separate consumer runs,
 * added in quadrature, although the mixed events reuse the data's Lambdas (the accepted over-coverage).
 * The error of 3<u v> is propagated; its covariance with R_d is not -- there is no <R u v> profile. That is
 * safe: every term in it carries dP/dD = -P/D, so the dropped term is at most ~2 P sigma_D/sigma_R ~ 10^-3
 * of the leading one.
 */

/// @brief One sample's inputs: <R>, 3<u v> and kappa, each with its error.
struct CorrectionInputs { std::pair<double, double> ring, uvx3, kappa; };

/// @brief A data sample and its fake estimate, where their consumer and signal-extraction files live, and how the
///        corrected curve is drawn wherever pairings are series (the differential canvases).
struct FakePairing { std::string consumerBase, sigBase, dataSuffix, fakeSuffix, dataLabel, pairLabel; VariationConfig style; };

/// @brief Fills the inputs of one sample (file suffix includes family and variation); false if unavailable.
using CorrectionFetcher = std::function<bool(const FakePairing&, const std::string&, const std::string&, CorrectionInputs&)>;

/// @brief The inputs of one candidate set (all candidates, or one mass region): kappa keeps the exact <u v w> covariance.
bool InputsOfSample(const RingSample& smp, CorrectionInputs& in) {
    if (smp.ring.n < 2.) return false;
    in.ring  = MeanOf(smp.ring);
    in.uvx3  = {3. * smp.num.Mean(), 3. * std::sqrt(std::max(smp.num.Var(), 0.) / smp.num.n)};
    in.kappa = KappaOf(smp);
    return true;
}

/// @brief One mass region of a RegionSample as a sample (index 1 peak, 0 sideband); ringDen comes along when booked.
RingSample RegionAsSample(const RegionSample& r, int region) {
    RingSample smp;
    smp.ring = r.ring[region]; smp.num = r.num[region]; smp.den = r.den[region];
    smp.numDen = r.numDen[region]; smp.ringDen = r.ringDen[region];
    return smp;
}

/**
 * @brief The sideband-subtracted inputs: every moment as Sum_P - Sum_B, then the usual correction.
 *
 * kappa_S = 3 <u v>_S/<w>_S is a ratio of two subtracted sums, with the exact covariance of CheapSignalRatio.
 */
bool InputsOfSidebandSubtracted(const RegionSample& r, CorrectionInputs& in) {
    if (!(r.ring[1].n > r.ring[0].n)) return false; // A peak no richer than its sideband has no signal to estimate
    const auto ring = CheapSignalMean(r.ring[1], r.ring[0]);
    const auto u    = CheapSignalMean(r.num[1], r.num[0]);
    const auto k    = CheapSignalRatio(r.num[1], r.num[0], r.den[1], r.den[0], r.numDen[1], r.numDen[0]);
    if (ring.second < 0. || u.second < 0. || k.second < 0.) return false;
    in.ring  = ring;
    in.uvx3  = {3. * u.first, 3. * u.second};
    in.kappa = {3. * k.first, 3. * k.second};
    return true;
}

/// @brief P = (R_d - R_f) / (D0 - k R_d R_f), with its first-order error.
std::pair<double, double> CorrectRing(const std::pair<double, double>& rd, const std::pair<double, double>& rf,
                                      const std::pair<double, double>& d0, double k) {
    const double den = d0.first - k * rd.first * rf.first;
    if (!(den > 0.)) return {0., -1.};
    const double p   = (rd.first - rf.first) / den;
    const double dRd = 1. / den + p * k * rf.first / den;  // dP/dR_d
    const double dRf = -1. / den + p * k * rd.first / den; // dP/dR_f
    const double dD0 = -p / den;                           // dP/dD0
    return {p, std::sqrt(dRd * dRd * rd.second * rd.second + dRf * dRf * rf.second * rf.second +
                         dD0 * dD0 * d0.second * d0.second)};
}

/// @brief Tag of the resampled mixed-event wagons (nProxyResamples in the consumer), after the usual mixed-event suffix.
const std::string MIXED_RESAMPLES_TAG = "_10resamples";

/// @brief The count in a resampling tag ("10" from "_10resamples"), so labels take it from the tag and never disagree.
std::string ResampleCountOf(const std::string& tag) {
    const size_t first = tag.find_first_of("0123456789");
    if (first == std::string::npos) return "";
    return tag.substr(first, tag.find_first_not_of("0123456789", first) - first);
}

/// @brief Data / fake pairs: the in-wagon two, the same two against resampled mixed events, then each consumer-style
///        reference against its own mixed events.
/// @note A reference's signal extraction is looked for where the pipeline puts it: <base>/../results_SigExtract.
std::vector<FakePairing> MakeFakePairings(const std::string& consumerDir, const std::string& sigExtractDir,
                                          const std::vector<ExternalRef>& externals) {
    const std::string resampled = "," + ResampleCountOf(MIXED_RESAMPLES_TAG) + "Resamp.";
    // Styles: black for the single borrow, the resampled mixed events' dark grey and diamond; AN-Cuts dashed and open
    std::vector<FakePairing> pairings = {
        {consumerDir, sigExtractDir, "", "_MixedEventProxies", "Data", "Data - MixedEv",
         {"", "Data - MixedEv", kBlack, 1, 20, true}},
        {consumerDir, sigExtractDir, "_analysisCuts", "_MixedEventProxies_analysisCuts", "Data (AN-Cuts)", "AN-Cuts: Data - MixedEv",
         {"", "AN-Cuts: Data - MixedEv", kBlack, 2, 24, false}},
        // Same data, same Lambdas, each now against several distinct mixed jets: the fake with less jet-sampling noise
        {consumerDir, sigExtractDir, "", "_MixedEventProxies" + MIXED_RESAMPLES_TAG, "Data", "Data - MixedEv" + resampled,
         {"", "Data - MixedEv" + resampled, kGray + 3, 1, 33, false}},
        {consumerDir, sigExtractDir, "_analysisCuts", "_MixedEventProxies_analysisCuts" + MIXED_RESAMPLES_TAG, "Data (AN-Cuts)",
         "AN-Cuts: Data - MixedEv" + resampled, {"", "AN-Cuts: Data - MixedEv" + resampled, kGray + 3, 2, 27, false}}};
    for (const auto& ext : externals)
        if (!ext.basePath.empty() && ext.kind == ExternalRef::kConsumerDir) {
            VariationConfig style = ext.config;
            style.legendLabel = ext.config.legendLabel + " - MixedEv";
            style.isData = false;
            pairings.push_back({ext.basePath, ext.basePath + "/../results_SigExtract", "", "_MixedEventProxies",
                                ext.config.legendLabel, style.legendLabel, style});
        }
    return pairings;
}

/// @brief One series per ring definition: the full ring first, then LONGITUDINAL_RINGS in order (shared styles).
std::vector<VariationConfig> MakeDefinitionSeries(const std::string& fullLabel, const std::string& longSymbol) {
    std::vector<VariationConfig> series = {{"", fullLabel, kGray + 2, 1, 24, false}};
    for (size_t l = 0; l < LONGITUDINAL_RINGS.size(); ++l)
        series.push_back({"", longSymbol + LONGITUDINAL_RINGS[l].symbolSuffix, LONG_RING_COLOR[l % 2], 1, LONG_RING_MARKER[l % 2], true});
    return series;
}

/**
 * @brief Writes Corrections/<folder>/<Family>/: the kappa of every data-like sample and the corrected P per pair.
 *
 * The same code serves the raw and the signal-extracted inputs; only the fetcher differs. The full ring is
 * drawn wherever its inputs exist, and R_z joins it wherever the R_z files do.
 *
 * @param families Full-ring families: the longitudinal samples of each are the same suffix plus each tag.
 */
void DrawCorrections(TDirectory* fileRoot, const std::vector<FamilyConfig>& families,
                     const std::vector<FakePairing>& pairings, const std::string& folder,
                     const std::string& tag, const CorrectionFetcher& fetch) {
    // One series per ring definition: the full ring first, then LONGITUDINAL_RINGS in order
    const std::vector<VariationConfig> kappaSeries = MakeDefinitionSeries("#kappa_{eff} (Full 3D Ring)", "#kappa");
    const std::vector<VariationConfig> corrSeries = MakeDefinitionSeries("Full 3D Ring", "R");
    std::vector<std::string> fileTags = {""};
    for (const auto& lr : LONGITUDINAL_RINGS) fileTags.push_back(lr.fileTag);
    const size_t nDef = fileTags.size();
    const std::pair<double, double> kMissing = {0., -1.};
    // Corrected rows of every family, kept for the species combination at the end: [family + "|" + proxy] --> rows
    struct CorrectedRow { std::string pairLabel; std::vector<std::pair<double, double>> perDefinition; };
    std::map<std::string, std::vector<CorrectedRow>> rowsByFamilyProxy;

    for (const auto& fam : families) {
        TDirectory* famDir = nullptr;
        for (const auto& proxy : PROXIES) {
            std::vector<std::string> kLabels, cLabels;
            std::vector<std::vector<std::pair<double, double>>> kVals(nDef), cVals(nDef);
            std::set<std::string> kappaSamplesDrawn; // Pairings share data samples, and kappa is the data's: one column each

            for (const auto& pr : pairings) {
                std::vector<std::pair<double, double>> kappa(nDef, kMissing), corrected(nDef, kMissing);
                bool anyData = false, anyCorrected = false;
                for (size_t d = 0; d < nDef; ++d) {
                    CorrectionInputs data, fake;
                    if (!fetch(pr, fam.dataSuffix + fileTags[d] + pr.dataSuffix, proxy, data)) continue;
                    anyData = true;
                    kappa[d] = data.kappa;
                    if (!fetch(pr, fam.dataSuffix + fileTags[d] + pr.fakeSuffix, proxy, fake)) continue;
                    // Exact inversion for the full ring (k = alpha^2/3), linear for the longitudinal rings (k = 0)
                    corrected[d] = CorrectRing(data.ring, fake.ring, data.uvx3, d == 0 ? fam.alphaSq / 3. : 0.);
                    if (corrected[d].second < 0.)
                        WarnSkip(pr.pairLabel + ", " + proxy + (d == 0 ? std::string(", full ring") : ", " + LONGITUDINAL_RINGS[d - 1].name),
                                 "correction denominator 3<u v> - k R_d R_f is not positive", "Corrections/" + folder + "/" + fam.familyName);
                    anyCorrected = anyCorrected || corrected[d].second >= 0.;
                }
                if (!anyData) continue;
                if (kappaSamplesDrawn.insert(pr.consumerBase + "|" + pr.dataSuffix).second) {
                    kLabels.push_back(pr.dataLabel);
                    for (size_t d = 0; d < nDef; ++d) kVals[d].push_back(kappa[d]);
                }
                if (anyCorrected) {
                    cLabels.push_back(pr.pairLabel);
                    for (size_t d = 0; d < nDef; ++d) cVals[d].push_back(corrected[d]);
                    rowsByFamilyProxy[fam.familyName + "|" + proxy].push_back({pr.pairLabel, corrected});
                }
                printf("    [Corrections/%s %s, %s, %s]", folder.c_str(), fam.familyName.c_str(), proxy.c_str(), pr.pairLabel.c_str());
                for (size_t d = 0; d < nDef; ++d)
                    printf("  %s: kappa = %.4f, P = %+.3e +/- %.1e", d == 0 ? "R" : LONGITUDINAL_RINGS[d - 1].shortName.c_str(),
                           kappa[d].first, corrected[d].first, corrected[d].second);
                printf("\n");
            }
            if (kLabels.empty()) {
                CountSkip("Corrections/" + folder, "no data sample with inputs", fam.familyName + "/" + proxy);
                continue;
            }
            if (!famDir) famDir = EnsureDir(EnsureDir(EnsureDir(fileRoot, "Corrections"), folder), fam.familyName);
            DrawSeriesByColumn(kLabels, kappaSeries, kVals, "Canvas_Kappa_" + proxy,
                               fam.familyName + ", " + proxy + ": #kappa of the data sample (" + tag + ")", famDir, "#kappa");
            if (!cLabels.empty())
                DrawSeriesByColumn(cLabels, corrSeries, cVals, "Canvas_Corrected_" + proxy,
                                   fam.familyName + ", " + proxy + ": fake-subtracted ring amplitude (" + tag + ")", famDir, "P_{R}");
            else
                CountSkip("Corrections/" + folder + ", corrected P", "no pairing with both data and fake inputs",
                          fam.familyName + "/" + proxy);
        }
    }

    // ---- Species combination (LongitudinalRing.tex, Prop. 6.2) ----
    // A topological fake has the same unitless moment for Lambda and antiLambda, so it flips with alpha, while a
    // vorticity ring does not: the average of the two corrected results suppresses such fakes (to 0.6%, from
    // |alpha_Lambda| != |alpha_antiLambda|), and the half-difference holds them, with no vorticity signal. The two
    // species are disjoint samples, so their errors add in quadrature. Charge-driven fakes and a magnetic-field
    // polarization do NOT cancel; they land in the half-difference.
    for (const auto& proxy : PROXIES) {
        const auto itL = rowsByFamilyProxy.find("Lambda|" + proxy);
        const auto itA = rowsByFamilyProxy.find("AntiLambda|" + proxy);
        if (itL == rowsByFamilyProxy.end() || itA == rowsByFamilyProxy.end()) {
            CountSkip("Corrections/" + folder + "/SpeciesCombination", "no corrected rows for both Lambda and AntiLambda", proxy);
            continue;
        }
        std::vector<std::string> labels;
        std::vector<std::vector<std::pair<double, double>>> avg(nDef), half(nDef);
        for (const auto& rowL : itL->second)
            for (const auto& rowA : itA->second) {
                if (rowA.pairLabel != rowL.pairLabel) continue;
                labels.push_back(rowL.pairLabel);
                for (size_t d = 0; d < nDef; ++d) {
                    const auto& pl = rowL.perDefinition[d];
                    const auto& pa = rowA.perDefinition[d];
                    if (pl.second < 0. || pa.second < 0.) { avg[d].push_back(kMissing); half[d].push_back(kMissing); continue; }
                    const double err = 0.5 * std::hypot(pl.second, pa.second);
                    avg[d].push_back({0.5 * (pl.first + pa.first), err});
                    half[d].push_back({0.5 * (pl.first - pa.first), err});
                }
            }
        if (labels.empty()) {
            CountSkip("Corrections/" + folder + "/SpeciesCombination", "no pairing corrected in both Lambda and AntiLambda", proxy);
            continue;
        }
        TDirectory* spDir = EnsureDir(EnsureDir(EnsureDir(fileRoot, "Corrections"), folder), "SpeciesCombination");
        DrawSeriesByColumn(labels, corrSeries, avg, "Canvas_SpeciesAverage_" + proxy,
                           proxy + ": #frac{1}{2}(#Lambda + #bar{#Lambda}), topological fakes cancel (" + tag + ")", spDir, "P_{R}");
        DrawSeriesByColumn(labels, corrSeries, half, "Canvas_SpeciesHalfDifference_" + proxy,
                           proxy + ": #frac{1}{2}(#Lambda - #bar{#Lambda}), topological fakes only (" + tag + ")", spDir, "#Delta P_{R}");
    }
}

/**
 * @brief Corrections/CheapSigExtract/<Family>/Canvas_KappaByRegion_<Proxy>: the raw kappa of the peak and of the
 *        sideband, for every ring definition, per data-like sample.
 *
 * Each is 3 <u v>/<w> over that region's candidates, with the ratio-of-means error that keeps their covariance.
 * Where the two regions disagree, background candidates do not share the signal's response, which is exactly
 * why SidebandSubtracted/ divides by the signal kappa, and InPeakQA/ and OffPeakQA/ each by their own region's.
 */
void DrawKappaByRegion(FileCache& cache, TDirectory* fileRoot, const std::vector<FamilyConfig>& families,
                       const std::vector<FakePairing>& pairings) {
    // Two series (peak, sideband) per ring definition: the full ring first, then LONGITUDINAL_RINGS in order
    static const Color_t kColor[3] = {kGray + 2, kBlue + 1, kGreen + 2};
    std::vector<VariationConfig> series;
    std::vector<std::string> fileTags = {""}, symbols = {"_{eff}"};
    for (const auto& lr : LONGITUDINAL_RINGS) { fileTags.push_back(lr.fileTag); symbols.push_back(lr.symbolSuffix); }
    for (size_t d = 0; d < fileTags.size(); ++d) {
        series.push_back({"", "#kappa" + symbols[d] + ", peak",     kColor[d % 3], 1, 20, d > 0});
        series.push_back({"", "#kappa" + symbols[d] + ", sideband", kColor[d % 3], 1, 24, d > 0});
    }
    const std::pair<double, double> kMissing = {0., -1.};
    auto kappaOf = [](const RegionSample& r, int i) {
        const auto k = RatioOfMeans(r.num[i], r.den[i], r.numDen[i].Mean());
        return std::make_pair(3. * k.first, 3. * k.second);
    };

    for (const auto& fam : families) {
        for (const auto& proxy : PROXIES) {
            std::vector<std::string> labels;
            std::vector<std::vector<std::pair<double, double>>> vals(series.size());
            std::set<std::string> samplesDrawn; // Pairings share data samples: one column each
            for (const auto& pr : pairings) {
                if (!samplesDrawn.insert(pr.consumerBase + "|" + pr.dataSuffix).second) continue;
                const std::string stem = pr.consumerBase + "/ConsumerResults_" + fam.dataSuffix;
                std::vector<std::pair<double, double>> row(series.size(), kMissing);
                bool any = false;
                for (size_t d = 0; d < fileTags.size(); ++d) {
                    RegionSample r;
                    if (!LoadRegionSample(cache, stem + fileTags[d] + pr.dataSuffix + ".root", proxy, r)) continue;
                    any = true;
                    row[2 * d] = kappaOf(r, 1);     // Peak
                    row[2 * d + 1] = kappaOf(r, 0); // Sideband
                }
                if (!any) continue;
                labels.push_back(pr.dataLabel);
                for (size_t k = 0; k < series.size(); ++k) vals[k].push_back(row[k]);
            }
            if (labels.empty()) {
                CountSkip("Corrections/CheapSigExtract, kappa by region", "no data sample with mass-region inputs", fam.familyName + "/" + proxy);
                continue;
            }
            TDirectory* famDir = EnsureDir(EnsureDir(EnsureDir(fileRoot, "Corrections"), "CheapSigExtract"), fam.familyName);
            DrawSeriesByColumn(labels, series, vals, "Canvas_KappaByRegion_" + proxy,
                               fam.familyName + ", " + proxy + ": #kappa in the mass peak and in the sideband (raw)", famDir, "#kappa");
        }
    }
}

// ===============================
// Differential corrections and ring comparisons (vs proxy pT, vs centrality)
// ===============================

/// @brief The consumer's differential variables ("ProxyPt" is each proxy's own pT): binning and titles come from its files.
const std::vector<std::string> DIFFERENTIAL_VARIABLES = {"ProxyPt", "Centrality"};

/// @brief A canvas-safe name: letters and digits only ("Data (AN-Cuts)" --> "DataANCuts").
std::string CanvasSafe(const std::string& s) {
    std::string out;
    for (char c : s)
        if (std::isalnum(static_cast<unsigned char>(c))) out += c;
    return out;
}

/// @brief An empty TH1D on the consumer's own x binning. Bins left unset are skipped when drawn, never parked at zero.
TH1D* MakeVsXHist(const std::string& stem, const std::vector<double>& edges) {
    static int counter = 0; // Unique names: no "Replacing existing TH1" across families, proxies and samples
    TH1D* h = new TH1D((stem + "_" + std::to_string(counter++)).c_str(), "", static_cast<int>(edges.size()) - 1, edges.data());
    h->SetDirectory(nullptr);
    return h;
}

/// @brief Sets one bin from a (value, error) pair; a negative error marks a missing value and leaves the bin empty.
void SetVsXPoint(TH1D* h, int bin, const std::pair<double, double>& v) {
    if (!h || v.second < 0.) return;
    h->SetBinContent(bin, v.first);
    h->SetBinError(bin, v.second);
}

/// @brief The inputs of one x bin of a sample: the integrated fetchers' formulas, applied per bin.
using DifferentialBinInputs = std::function<bool(const DifferentialSample&, size_t, CorrectionInputs&)>;

/**
 * @brief Writes Corrections/<folder>/Differential/<Family>/<Proxy>/: the corrected P and kappa vs proxy pT and vs
 *        centrality, bin by bin with the same formulas as the integrated folder.
 *
 * Canvas_Corrected_Vs<Var>_<Definition>: one curve per pairing, so the estimates of the fake compare at a glance.
 * Canvas_Kappa_Vs<Var>_<Sample>: one curve per ring definition, for every data sample. Each bin holds the same Lambdas
 * in data and in its mixed events, as the consumer bins on the event mixing axes; a pairing whose two files disagree
 * on the binning is skipped with a warning rather than subtracted across different bins.
 */
void DrawDifferentialCorrections(FileCache& cache, TDirectory* fileRoot, const std::vector<FamilyConfig>& families,
                                 const std::vector<FakePairing>& pairings,
                                 const std::string& folder, const std::string& tag, const DifferentialBinInputs& binInputs) {
    const std::vector<VariationConfig> kappaSeries = MakeDefinitionSeries("#kappa_{eff} (Full 3D Ring)", "#kappa");
    const std::vector<VariationConfig> corrSeries = MakeDefinitionSeries("Full 3D Ring", "R");
    std::vector<std::string> fileTags = {""}, defNames = {"FullRing"};
    for (const auto& lr : LONGITUDINAL_RINGS) { fileTags.push_back(lr.fileTag); defNames.push_back(lr.name); }
    const size_t nDef = fileTags.size();
    const ProfileConfig pc{}; // Every title is passed explicitly

    for (const auto& fam : families)
        for (const auto& proxy : PROXIES)
            for (const auto& var : DIFFERENTIAL_VARIABLES) {
                TDirectory* pxDir = nullptr;
                auto dirOf = [&]() {
                    if (!pxDir) pxDir = EnsureDir(fileRoot, "Corrections/" + folder + "/Differential/" + fam.familyName + "/" + proxy);
                    return pxDir;
                };
                std::string xTitle;
                std::vector<std::vector<ProfileBundle>> corrBundles(nDef); // [definition] --> one curve per pairing
                std::vector<TH1*> owned;
                std::set<std::string> kappaSamplesDrawn; // Pairings share data samples: one kappa canvas each

                for (const auto& pr : pairings) {
                    std::vector<TH1D*> kappa(nDef, nullptr), corrected(nDef, nullptr);
                    for (size_t d = 0; d < nDef; ++d) {
                        const std::string stem = pr.consumerBase + "/ConsumerResults_" + fam.dataSuffix + fileTags[d];
                        DifferentialSample dataS, fakeS;
                        if (!LoadDifferentialSample(cache, stem + pr.dataSuffix + ".root", proxy, var, dataS)) continue;
                        xTitle = dataS.xTitle;
                        kappa[d] = MakeVsXHist("hKappa", dataS.edges);
                        owned.push_back(kappa[d]);
                        bool haveFake = LoadDifferentialSample(cache, stem + pr.fakeSuffix + ".root", proxy, var, fakeS);
                        if (haveFake && fakeS.edges != dataS.edges) {
                            WarnSkip(pr.pairLabel + ", " + proxy + " vs " + var, "data and fake binnings differ",
                                     "Corrections/" + folder + "/Differential/" + fam.familyName);
                            haveFake = false;
                        }
                        if (haveFake) {
                            corrected[d] = MakeVsXHist("hCorrected", dataS.edges);
                            owned.push_back(corrected[d]);
                        }
                        for (size_t i = 0; i < dataS.all.size(); ++i) {
                            CorrectionInputs data, fake;
                            if (!binInputs(dataS, i, data)) continue;
                            SetVsXPoint(kappa[d], static_cast<int>(i) + 1, data.kappa);
                            // Exact inversion for the full ring (k = alpha^2/3), linear for the longitudinal rings (k = 0)
                            if (haveFake && binInputs(fakeS, i, fake))
                                SetVsXPoint(corrected[d], static_cast<int>(i) + 1,
                                            CorrectRing(data.ring, fake.ring, data.uvx3, d == 0 ? fam.alphaSq / 3. : 0.));
                        }
                    }
                    for (size_t d = 0; d < nDef; ++d)
                        if (corrected[d]) corrBundles[d].push_back({corrected[d], pr.style});

                    // Kappa belongs to the data sample, so each is drawn once, all definitions together
                    std::vector<ProfileBundle> kb;
                    for (size_t d = 0; d < nDef; ++d)
                        if (kappa[d]) kb.push_back({kappa[d], kappaSeries[d]});
                    if (kb.empty() || !kappaSamplesDrawn.insert(pr.consumerBase + "|" + pr.dataSuffix).second) continue;
                    DrawComparisonCanvas(kb, "Canvas_Kappa_Vs" + var + "_" + CanvasSafe(pr.dataLabel),
                                         fam.familyName + ", " + proxy + ": #kappa of " + pr.dataLabel + " (" + tag + ")",
                                         dirOf(), pc, xTitle, "#kappa");
                }

                for (size_t d = 0; d < nDef; ++d) {
                    if (corrBundles[d].empty()) {
                        CountSkip("Corrections/" + folder + "/Differential, corrected P", "no pairing with both data and fake inputs",
                                  fam.familyName + "/" + proxy + "/" + var + "/" + defNames[d]);
                        continue;
                    }
                    DrawComparisonCanvas(corrBundles[d], "Canvas_Corrected_Vs" + var + "_" + defNames[d],
                                         fam.familyName + ", " + proxy + ": " + corrSeries[d].legendLabel +
                                             ", fake-subtracted ring amplitude (" + tag + ")",
                                         dirOf(), pc, xTitle, "P_{R}");
                }
                for (TH1* h : owned) delete h;
            }
}

/**
 * @brief RingVsLongitudinal{,_CheapSigExtract}/Differential/<Family>/<Proxy>/: the five series of the integrated
 *        comparison (<R>, then <R_long> and <R_long>/<n_z^2> per definition) vs proxy pT and vs centrality.
 *
 * One canvas per sample, so only the samples of the Data - MixedEv pairings are drawn (data and mixed events alike),
 * as every systematic would multiply them. Raw uses all candidates; CheapSigExtract the peak, the sideband and the
 * in - out of the means, exactly as the integrated version (the in - out of a ratio is not drawn).
 * A _DataOnly canvas then overlays the in-wagon data samples (Data, AN-Cuts), each in its own line style, so the
 * selections compare by eye on every series at once.
 */
void DrawDifferentialRingVsLongitudinal(FileCache& cache, TDirectory* fileRoot, const std::string& consumerDir,
                                        const std::vector<FamilyConfig>& families, const VariationConfig& dataConfig,
                                        const std::vector<VariationConfig>& sysVariations, const std::vector<ExternalRef>& externals,
                                        const std::vector<FakePairing>& pairings) {
    const std::vector<VariationConfig> series = MakeRingVsLongitudinalSeries();
    const ProfileConfig pc{}; // Every title is passed explicitly

    // The samples of the pairings, styled and labelled by the column table (one source for both)
    std::set<std::string> wanted;
    for (const auto& pr : pairings) { wanted.insert(pr.consumerBase + "|" + pr.dataSuffix); wanted.insert(pr.consumerBase + "|" + pr.fakeSuffix); }
    std::vector<RingPairColumn> columns;
    for (const auto& col : MakeRingPairColumns(consumerDir, dataConfig, sysVariations, externals))
        if (wanted.count(col.base + "|" + col.suffix)) columns.push_back(col);
    // The in-wagon data samples, which also go to the _DataOnly canvas
    std::set<std::string> dataOnly;
    for (const auto& pr : pairings)
        if (pr.consumerBase == consumerDir) dataOnly.insert(pr.consumerBase + "|" + pr.dataSuffix);

    // States: all candidates, then the CheapSigExtract peak, sideband and in - out (series 0 is <R>, 1 + 2l <R_long>, 2 + 2l the ratio)
    struct State { std::string folder, nameTag, titleTag; };
    const std::vector<State> states = {{"RingVsLongitudinal", "", "raw, all candidates"},
                                       {"RingVsLongitudinal_CheapSigExtract", "_InMassPeak", "in mass peak"},
                                       {"RingVsLongitudinal_CheapSigExtract", "_OutOfMassPeak", "out of mass peak"},
                                       {"RingVsLongitudinal_CheapSigExtract", "_InMinusOut", "in - out of mass peak"}};
    // <R> (or <R_long>) and the ratio of one candidate set; the ratio keeps the R-w covariance
    auto meanAndRatio = [](const RingSample& s, std::pair<double, double>& mean, std::pair<double, double>& ratio) {
        mean = (s.ring.n >= 2.) ? MeanOf(s.ring) : std::make_pair(0., -1.);
        ratio = RatioOfMeans(s.ring, s.den, s.ringDen.Mean());
    };

    for (const auto& fam : families)
        for (const auto& proxy : PROXIES)
            for (const auto& var : DIFFERENTIAL_VARIABLES) {
                // Curves of the data samples, kept past their own canvases for the _DataOnly one: [sample] --> [state][series]
                struct KeptCurves { std::string label; std::vector<std::vector<TH1D*>> h; };
                std::vector<KeptCurves> dataCurves;
                std::string xTitle;
                for (const auto& col : columns) {
                    const std::string stem = col.base + "/ConsumerResults_" + fam.dataSuffix;
                    DifferentialSample full;
                    if (!LoadDifferentialSample(cache, stem + col.suffix + ".root", proxy, var, full)) continue;
                    std::vector<DifferentialSample> lng(LONGITUDINAL_RINGS.size());
                    std::vector<bool> haveLong(LONGITUDINAL_RINGS.size(), false);
                    bool anyLong = false;
                    for (size_t l = 0; l < LONGITUDINAL_RINGS.size(); ++l) {
                        haveLong[l] = LoadDifferentialSample(cache, stem + LONGITUDINAL_RINGS[l].fileTag + col.suffix + ".root",
                                                             proxy, var, lng[l]);
                        if (haveLong[l] && lng[l].edges != full.edges) {
                            WarnSkip(LONGITUDINAL_RINGS[l].name + " vs " + var + ", " + proxy + ", " + col.cfg.legendLabel,
                                     "binning differs from the full ring's", fam.familyName);
                            haveLong[l] = false;
                        }
                        anyLong = anyLong || haveLong[l];
                    }
                    if (!anyLong) continue;

                    // [state][series] --> one curve over x
                    std::vector<std::vector<TH1D*>> h(states.size(), std::vector<TH1D*>(series.size(), nullptr));
                    for (auto& row : h)
                        for (auto& p : row) p = MakeVsXHist("hRingVsLong", full.edges);
                    for (size_t i = 0; i < full.all.size(); ++i) {
                        const int bin = static_cast<int>(i) + 1;
                        // The three candidate sets of this bin, for one file: all, peak (region 1), sideband (region 0)
                        auto setsOf = [&](const DifferentialSample& d) {
                            return std::vector<RingSample>{d.all[i], RegionAsSample(d.region[i], 1), RegionAsSample(d.region[i], 0)};
                        };
                        const std::vector<RingSample> fullSets = setsOf(full);
                        std::pair<double, double> mean, ratio;
                        for (size_t st = 0; st < 3; ++st) {
                            meanAndRatio(fullSets[st], mean, ratio);
                            SetVsXPoint(h[st][0], bin, mean);
                        }
                        for (size_t l = 0; l < LONGITUDINAL_RINGS.size(); ++l) {
                            if (!haveLong[l]) continue;
                            const std::vector<RingSample> sets = setsOf(lng[l]);
                            for (size_t st = 0; st < 3; ++st) {
                                meanAndRatio(sets[st], mean, ratio);
                                SetVsXPoint(h[st][1 + 2 * l], bin, mean);
                                SetVsXPoint(h[st][2 + 2 * l], bin, ratio);
                            }
                        }
                        // In - out of the means: disjoint candidates, so plain quadrature
                        for (size_t s = 0; s < series.size(); s += (s == 0 ? 1 : 2)) {
                            const double eIn = h[1][s]->GetBinError(bin), eOut = h[2][s]->GetBinError(bin);
                            if (eIn > 0. && eOut > 0.)
                                SetVsXPoint(h[3][s], bin, {h[1][s]->GetBinContent(bin) - h[2][s]->GetBinContent(bin), std::hypot(eIn, eOut)});
                        }
                    }

                    for (size_t st = 0; st < states.size(); ++st) {
                        std::vector<ProfileBundle> bundles;
                        for (size_t s = 0; s < series.size(); ++s) bundles.push_back({h[st][s], series[s]});
                        TDirectory* pxDir = EnsureDir(fileRoot, states[st].folder + "/Differential/" + fam.familyName + "/" + proxy);
                        DrawComparisonCanvas(bundles, "Canvas_RingVsLongitudinal_Vs" + var + states[st].nameTag + "_" + CanvasSafe(col.cfg.legendLabel),
                                             fam.familyName + ", " + proxy + ", " + col.cfg.legendLabel + ": Full 3D Ring vs longitudinal rings (" +
                                                 states[st].titleTag + ")",
                                             pxDir, pc, full.xTitle, st == 3 ? "#Delta (in - out)" : "value");
                    }
                    if (dataOnly.count(col.base + "|" + col.suffix)) {
                        dataCurves.push_back({col.cfg.legendLabel, h});
                        xTitle = full.xTitle;
                        continue; // Deleted after the _DataOnly canvas
                    }
                    for (auto& row : h)
                        for (TH1D* p : row) delete p;
                }

                // Data samples only: the series' own styles, one line style per sample (solid for the first, i.e. Data)
                if (dataCurves.size() >= 2)
                    for (size_t st = 0; st < states.size(); ++st) {
                        std::vector<ProfileBundle> bundles;
                        for (size_t k = 0; k < dataCurves.size(); ++k)
                            for (size_t s = 0; s < series.size(); ++s) {
                                VariationConfig cfg = series[s];
                                cfg.legendLabel = series[s].legendLabel + ", " + dataCurves[k].label;
                                cfg.lineStyle = 1 + static_cast<int>(k);
                                cfg.isData = cfg.isData && k == 0; // Emphasis stays on the first sample
                                bundles.push_back({dataCurves[k].h[st][s], cfg});
                            }
                        TDirectory* pxDir = EnsureDir(fileRoot, states[st].folder + "/Differential/" + fam.familyName + "/" + proxy);
                        DrawComparisonCanvas(bundles, "Canvas_RingVsLongitudinal_Vs" + var + states[st].nameTag + "_DataOnly",
                                             fam.familyName + ", " + proxy + ", data samples: Full 3D Ring vs longitudinal rings (" +
                                                 states[st].titleTag + ")",
                                             pxDir, pc, xTitle, st == 3 ? "#Delta (in - out)" : "value");
                    }
                for (auto& kept : dataCurves)
                    for (auto& row : kept.h)
                        for (TH1D* p : row) delete p;
            }
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
 * @param doIndividualComparisons Master switch for the one-vs-one canvases. Off by default: that block
 *                     alone is close to half of the entire output, and the overlays plus the pull panels
 *                     already carry the same comparisons. Per-observable opt-in still applies on top of
 *                     this, so turning it on enables it only for the profiles flagged in the table.
 * @param sigExtractDir Directory holding the signalExtractionRing_*.root files. Empty disables every
 *                     signal-extracted block.
 */
void auxiliarySummaryPlots(const std::string& consumerDir,
                           const std::string& mcRefDir = "",
                           const std::string& ppRefDir = "",
                           const std::string& toyModelPath = "",
                           bool doIndividualComparisons = false,
                           const std::string& sigExtractDir = "") {
    
    // 1. Define Systematic Variations (the list of all useful variations I would like to track into this plot)
    // The data config (empty suffix) is handled separately in the logic to ensure it is always first
    // Field order after isData: inRedux, inPtGateStudy, inEtaGateStudy.
    const std::string resampCount = ResampleCountOf(MIXED_RESAMPLES_TAG); // Legends take the count from the file tag
    std::vector<VariationConfig> sysVariations = {
        // suffix                          legend                      colour     ls mk  data   redux  ptGate etaGate
        {"_forceRandJet",                 "Rand Jet",                  kBlue,     1, 20, false, false,  true,  true}, // Full circle
        // Data-like is off the reduced canvas because it is the crudest of the mixing proxies: a
        // cheap event mixing, cheaper even than prevJet, which does not carry the correlations
        // between the angular distributions in full. It remains on the full systematics canvas and
        // in the pT gate study.
        {"_forceDatalikeJet",               "Data-like Jet",             kRed,      1, 21, false, false, true,  false}, // Square
        {"_forcePerpToJet",                 "Perp to Jet",               kGreen+2,  1, 22, false, true,  true,  true},  // Triangle up
        {"_forcePreviousJet",               "Prev Jet",                  kOrange+1, 1, 34, false, false, false, false}, // Full cross (bold "+")
        {"_MixedEventProxies",              "MixedEv Jet",               kGray+2,   1, 24, false, true,  false, false}, // Open circle
        {"_MixedEventProxies_10ColGap",     "MixedEv Jet (10coll gap)",  kGray+1,   2, 25, false, false, false, false}, // Open circle, dashed line
        {"_analysisCuts",                   "Data (AN-Cuts)",            kBlack,    2, 24, false, true,  false, false}, // Open circle, dashed line
        {"_MixedEventProxies_analysisCuts", "MixedEv Jet (AN-Cuts)",     kGray+2,   2, 24, false, true,  false, false}, // Open circle
        // The same mixed events, each target now borrowing nProxyResamples distinct jets: darker grey and an open
        // diamond, solid or dashed exactly as its single-borrow parent (off the reduced canvases)
        {"_MixedEventProxies" + MIXED_RESAMPLES_TAG,               "MixedEv Jet (" + resampCount + " resamp.)",         kGray+3, 1, 27, false, false, false, false}, // Open diamond
        {"_MixedEventProxies_analysisCuts" + MIXED_RESAMPLES_TAG, "MixedEv Jet (AN-Cuts, " + resampCount + " resamp.)", kGray+3, 2, 27, false, false, false, false}, // Open diamond, dashed line

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
        {"_5cmZvtx",                      "Data, |Zvtx| < 5 cm",         kBlack,    1, 24, false, true,  false, false} // Open circle, similar to dataConfig, as this is also data.

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
    // Each family reads them through its own list (see externalsNoToy), so they are "All" here.
    std::vector<ExternalRef> externalsAll = {
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
        {"_excludeOutOfPeak", "InPeak (#mu#pm 1#sigma)",             kViolet+2, 1, 34, false}, // Cross
        {"_excludeInPeak",    "OutOfPeak (#mu#pm[6#sigma,7#sigma])", kViolet,   1, 28, false}  // Open cross
    };

    // 2. Define the Families
    std::vector<FamilyConfig> families = {
        {"Lambda",       "JustLambda",     AlphaSqLambda},
        {"AntiLambda",   "JustAntiLambda", AlphaSqAntiLambda},
        {"BothHyperons", "BothHyperons",   AlphaSqBoth} // BothHyperons_2GeVLeadP was already the standard LeadP pT, so I modified this to follow the usual
                                                        // naming scheme that was in place for other variations of BothHyperons configurables of the consumer
    };

    // Longitudinal families (README, "The longitudinal ring"). The consumer configs carry the definition's tag
    // right after the family, so a longitudinal output is just another family whose suffix ends in that tag:
    // every variation, reference and signal-extraction file then resolves exactly as for its full-ring sibling.
    // Added only where the data file exists; fullFamilies keeps the originals for the blocks that pair definitions.
    const std::vector<FamilyConfig> fullFamilies = families;
    std::vector<FamilyConfig> siblingsNotRun; // No data file: listed once in the banner, their files never read
    for (const auto& lr : LONGITUDINAL_RINGS)
        for (const auto& fam : fullFamilies) {
            FamilyConfig sibling{fam.familyName + "_" + lr.name, fam.dataSuffix + lr.fileTag, fam.alphaSq, &lr};
            if (!gSystem->AccessPathName((consumerDir + "/ConsumerResults_" + sibling.dataSuffix + ".root").c_str()))
                families.push_back(sibling);
            else
                siblingsNotRun.push_back(sibling);
        }

    // 3. Define the Observable Groups
    // Split by what the x-axis actually is: the proxy's own eta, or the hyperon's eta.
    //
    // Every entry is the consumer's in-file path, spelled as the consumer books it. Every proxy folder holds
    // the same object names, so the proxy folder in the path is what tells the entries apart, and
    // ObservableDisplayName carries it into the output names. The axis titles are read from the data profile
    // in the observable loop, so the table does not repeat them.
    // Field order: path, doIndividualComparisons, toyProfile, doTanhFit, inBruteForce, displayName, xTitleSuffix.
    std::vector<ObservableGroup> observableGroups = {
        {
            "EtaProxy", // x-axis is the eta of the jet/leading-particle proxy
            {
                {"Ring/LeadP/EtaDependence/pRingObservableProxyEta",   true,  "pRingProxyJetVsEtaJet"},
                {"Ring/LeadJet/EtaDependence/pRingObservableProxyEta", false, "pRingProxyJetVsEtaJet"},
                {"Ring/SubJet/EtaDependence/pRingObservableProxyEta",  false, "pRingProxyJetVsEtaJet"},
                // Fine-binned counterparts of the LeadJet and LeadP ones above. Same axis, finer granularity:
                // these are the ones that show the tanh-like shape cleanly in the Toy Model.
                {"Ring/LeadJet/EtaDependence/pRingObservableProxyEtaHighEtaRes", false, "pRingProxyJetVsEtaJet", true},
                {"Ring/LeadP/EtaDependence/pRingObservableProxyEtaHighEtaRes",   true,  "pRingProxyJetVsEtaJet", true}
            },
            true // Eta axis is symmetric about zero, so folding is meaningful
        },
        {
            "EtaV0", // x-axis is the eta of the hyperon itself, split by which proxy defined the ring
            {
                {"Ring/LeadJet/EtaDependence/pRingObservableEtaLambda", false, "pRingProxyJetVsEta"},
                {"Ring/LeadP/EtaDependence/pRingObservableEtaLambda",   true,  "pRingProxyJetVsEta"},
                {"Ring/SubJet/EtaDependence/pRingObservableEtaLambda",  false, "pRingProxyJetVsEta"}
            },
            true
        },
        {
            // <R> against the angle between the hyperon and its decay proton in the helicity frame.
            // Booked under FakePolDrivers/PhiAEE/<Proxy>/ (doFakePolDiagnosticsQA), one species-agnostic
            // profile per proxy: the family's file sets the species.
            // <R> changes sign across this axis, so folding it would cancel the very structure the
            // observable exists to show: doFold stays false.
            // README: in the BothHyperons family this is the Lambda-like set, the QA for species competition --
            // with balanced Lambda and AntiLambda yields its dependence should cancel, so a residual there
            // measures the imbalance rather than a physics asymmetry.
            // The display names are spelled out to keep the folder names the AEE outputs have always had, minus the
            // species tag the family folder now carries (the derived ones would read LeadJetRingObservableVsPhiAEE).
            // The AEE mass block (runAeeHeeBlocks) takes them from here.
            "AEE",
            {
                {"FakePolDrivers/PhiAEE/LeadJet/pRingObservableVsPhiAEE", false, "", false, true, "RingLeadJetAEE"},
                {"FakePolDrivers/PhiAEE/LeadP/pRingObservableVsPhiAEE",   true,  "", false, true, "RingLeadPAEE"},
                {"FakePolDrivers/PhiAEE/SubJet/pRingObservableVsPhiAEE",  false, "", false, true, "RingSubJetAEE"}
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
                {"Ring/LeadJet/PtDependence/pRingVsProxyPt", false},
                {"Ring/LeadP/PtDependence/pRingVsProxyPt",   true},
                {"Ring/SubJet/PtDependence/pRingVsProxyPt",  false}
            },
            false // pT is strictly positive: folding about zero is meaningless
        },
        {
            "ProxyPt/EtaSplits", // The LeadP pT dependence, split by the sign of the proxy's or the hyperon's eta.
                                // Kept to the four single-sign splits: the four combined (proxy AND V0) splits exist
                                // in the consumer too, but they multiply the plot count for little extra information.
                                // No toy counterpart, for the same reason as ProxyPt above.
                                // The consumer states the split in the histogram title only, hence the x title suffix.
            {
                {"Ring/LeadP/PtDependence/pRingVsProxyPt_PosProxyEta", true, "", false, false, "", ", #eta_{LeadP}>0"},
                {"Ring/LeadP/PtDependence/pRingVsProxyPt_NegProxyEta", true, "", false, false, "", ", #eta_{LeadP}<0"},
                {"Ring/LeadP/PtDependence/pRingVsProxyPt_PosEtaV0",    true, "", false, false, "", ", #eta_{V0}>0"},
                {"Ring/LeadP/PtDependence/pRingVsProxyPt_NegEtaV0",    true, "", false, false, "", ", #eta_{V0}<0"}
            },
            false
        },
        {
            "PVz", // Primary-vertex z dependence, a pure detector-geometry handle on the fake signal
                   // No toy counterpart by design: the toy has no primary vertex to displace.
            {
                {"Ring/LeadJet/EventLevelDependence/pRingObservablePVz", false},
                {"Ring/LeadP/EventLevelDependence/pRingObservablePVz",   true},
                {"Ring/SubJet/EventLevelDependence/pRingObservablePVz",  false}
            },
            true // PVz is symmetric about zero, so folding it is meaningful
        },
        {
            "QA_Mult", // Event-activity dependence. Ties directly to the ratio-estimator discussion
                       // in the "Ring uncertainty estimators" note: the ring denominator is itself
                       // a fluctuating per-event multiplicity.
                       // No toy counterpart: the toy has neither a centrality nor a multiplicity concept.
            {
                {"Ring/LeadJet/EventLevelDependence/pRingVsCentrality"},
                // Booked under IntegratedCuts/ rather than Ring/: the leading-jet ring against the number of
                // Lambda-like V0s in the collision
                {"IntegratedCuts/LeadJet/pRingVsNV0s"}
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
                {"Ring/LeadJet/AngularDependence/pRingObservableProxyPhi"},
                {"Ring/LeadJet/AngularDependence/pRingObservablePhiLambda"}
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
                {"Ring/LeadJet/AngularDependence/pRingObservableDeltaPhi", false},
                {"Ring/LeadP/AngularDependence/pRingObservableDeltaPhi",   true},
                {"Ring/SubJet/AngularDependence/pRingObservableDeltaPhi",  false}
            },
            false // DeltaPhi runs over a full period, not symmetric about zero
        },
        {
            "OtherAngDepndncs/DeltaTheta", // Polar separation between the hyperon and the proxy
            //
            // TODO (Toy Model): same situation as DeltaPhi. cosThetaJet and the hyperon direction are both
            // already available at the ring computation, so <R> vs DeltaTheta is a booking plus one Fill.
            {
                {"Ring/LeadJet/AngularDependence/pRingObservableDeltaTheta", false},
                {"Ring/LeadP/AngularDependence/pRingObservableDeltaTheta",   true},
                {"Ring/SubJet/AngularDependence/pRingObservableDeltaTheta",  false}
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
    // so they need the labelled-axis drawing path. Booked under FakePolDrivers/ProxyEta/<Proxy>/ and only
    // when the consumer's doFakePolDiagnosticsQA is on, hence fetched defensively like everything else.
    // Field order: path, doEtaSignSummary, massSplitProfile.
    std::vector<CategoricalObservable> categoricalObservables = {
        {"FakePolDrivers/ProxyEta/LeadJet/pRingEtaCuts", true, ""},
        {"FakePolDrivers/ProxyEta/SubJet/pRingEtaCuts",  true, ""},
        {"FakePolDrivers/ProxyEta/LeadP/pRingEtaCuts",   true, "p2dRingEtaCuts_MassSignalVsBackground"}
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

    // Only a subset of the variations is run in R_z. Their absence is expected, so it is listed once here,
    // per family, instead of being warned about per observable (and those files are then not counted again).
    // A longitudinal definition never run for a family gets one line, and none of its files is read.
    for (const auto& fam : families) {
        if (!fam.ring) continue;
        std::string absent;
        for (const auto* vars : {&sysVariations, &massVariations})
            for (const auto& v : *vars) {
                const std::string f = consumerDir + "/ConsumerResults_" + fam.dataSuffix + v.suffix + ".root";
                if (!gSystem->AccessPathName(f.c_str())) continue;
                cache.MarkAbsent(f);
                absent += (absent.empty() ? "" : ", ") + v.suffix;
            }
        if (!absent.empty()) std::cout << " " << fam.familyName << ": variations not run: " << absent << "\n";
    }
    for (const auto& sibling : siblingsNotRun) {
        cache.MarkAbsent(consumerDir + "/ConsumerResults_" + sibling.dataSuffix + ".root");
        for (const auto* vars : {&sysVariations, &massVariations})
            for (const auto& v : *vars) {
                const std::string f = consumerDir + "/ConsumerResults_" + sibling.dataSuffix + v.suffix + ".root";
                if (gSystem->AccessPathName(f.c_str())) cache.MarkAbsent(f);
            }
        std::cout << " " << sibling.familyName << ": not run (no data file), skipped entirely\n";
    }

    std::cout << "========================================================\n";
    std::cout << " Generating Auxiliary Plots\n";
    std::cout << " Output:            " << outFilePath << "\n";
    // The references and the signal extraction are choices of the run: an absent one is stated here, once
    for (const auto& ext : externalsAll)
        std::cout << " " << ext.config.legendLabel << ": " << std::string(std::max(0, 17 - static_cast<int>(ext.config.legendLabel.size())), ' ')
                  << (ext.basePath.empty() ? "not requested" : ext.basePath) << "\n";
    std::cout << " Signal extraction: " << (sigExtractDir.empty() ? "disabled" : sigExtractDir) << "\n";
    std::cout << "========================================================\n";

    // One representative observable per proxy. Integrated values are a property of the entry set, not of
    // the binning, so a single observable per proxy captures everything the others would repeat.
    // {proxy, path}, in PROXIES order. Matched by the full path: the object name is the same in every proxy
    // folder, so on its own it would match all three.
    std::vector<std::pair<std::string, std::string>> proxyRepresentatives;
    for (const auto& proxy : PROXIES)
        proxyRepresentatives.push_back({proxy, RING_DIR_IN_FILE + proxy + "/EtaDependence/pRingObservableProxyEta"});

    // Styling for the cross-family comparison, where the SERIES are the hyperon selections
    const std::vector<VariationConfig> familyStyles = {
        {"", "Lambda",       kBlue,  1, 20, false},
        {"", "AntiLambda",   kRed,   1, 21, false},
        {"", "BothHyperons", kBlack, 1, 8,  true}
    };

    // The three proxy folders written by signalExtractionRing, in PROXIES order like
    // proxyRepresentatives above, so the two cross-family blocks line up.
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

    // crossFamilyVals[proxy][family] = one integrated value per column of crossFamilyLabels; a column the family
    // lacks keeps a negative error (a gap)
    std::vector<std::vector<std::vector<std::pair<double, double>>>>
        crossFamilyVals(proxyRepresentatives.size(),
                        std::vector<std::vector<std::pair<double, double>>>(families.size()));
    // The cross-family columns: Data, then every variation of the table, in table order
    std::vector<std::string> crossFamilyLabels{dataConfig.legendLabel};
    for (const auto& sys : sysVariations) crossFamilyLabels.push_back(sys.legendLabel);

    // 5. Main Processing Loops
    // The Toy Model predicts the full ring only, so longitudinal families see the references without it
    std::vector<ExternalRef> externalsNoToy = externalsAll;
    for (auto& ext : externalsNoToy)
        if (ext.kind == ExternalRef::kToyModel) ext.basePath.clear(); // Kept in place: grand-summary rows stay aligned

    for (size_t famIdx = 0; famIdx < families.size(); ++famIdx) {
        const auto& fam = families[famIdx];
        std::cout << " -> Processing Family: " << fam.familyName << std::endl;
        // This family's references: everything below reads them through this one name, so it needs no R_z branch
        const std::vector<ExternalRef>& externals = fam.ring ? externalsNoToy : externalsAll;

        // Directories are created lazily (only once something is actually written into them),
        // so a consumer file that lacks an object leaves no empty folders.
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
            // Read once per family, from the family's data file, which also sets the species: every object
            // here is species-agnostic. Both groups use the same mass-region axis and the same drawers; only
            // the x-axis differs. The titles are read from the source profiles, as in the observable loop.
            // README: these exist to be independent of signalExtractionRing, so that a disagreement
            // between the two is informative rather than circular.
            {
                const std::string famFile = consumerDir + "/ConsumerResults_" + fam.dataSuffix + ".root";

                // (a) <R> against phi_Lambda - phi_p*, by mass region, one map per entry of the AEE group.
                //     Destination folder is the 1D observable's own folder, so the cheap check sits next to
                //     the distribution it is checking rather than in a separate tree the reader has to
                //     cross-reference. The folder name comes from the table entry itself, so the two cannot
                //     drift apart.
                for (const auto& g : observableGroups) {
                    if (g.groupName != "AEE") continue;
                    for (const auto& aee : g.profiles) {
                        const std::string proxy = ProxyOfPath(aee.path);
                        TProfile2D* p2d = cache.FetchClone<TProfile2D>(famFile,
                            TASK_DIR_IN_FILE + "FakePolDrivers/PhiAEE/FiveBinMass/" + proxy + "/p2dRingObservableVsPhiAEEVs5BinMass");
                        if (!p2d) continue; // Booked only with doFakePolDiagnosticsQA. FetchClone has already warned
                        const std::string xTitle = p2d->GetXaxis()->GetTitle();
                        const std::string yTitle = WithProxySubscript(p2d->GetZaxis()->GetTitle(), proxy);
                        // Lands in <Family>/AEE/<observable>/Mass_Selection/
                        TDirectory* obsMassDir = EnsureDir(
                            EnsureDir(EnsureDir(EnsureDir(fOut, fam.familyName), "AEE"), ObservableDisplayName(aee)),
                            "Mass_Selection");
                        DrawMassRegionCanvases(p2d, fam.familyName + " " + yTitle + " vs " + xTitle,
                                               obsMassDir, xTitle, yTitle);
                        delete p2d;
                    }
                }

                // (b) <R> against cos(theta) in the helicity frame. Leading jet only: the consumer fills these
                //     with the leading-jet ring (their folder names no proxy), which is why the titles say so
                //     explicitly: there is no LeadP or SubJet counterpart to compare with.
                //     Motivated by the angular dependence reported in J. Adams' global-polarization
                //     thesis; the question is whether the ring observable shows anything similar.
                //     The species is the family's (its folder), so the canvas names carry none.
                const std::string heeBase = TASK_DIR_IN_FILE + "FakePolDrivers/CosHEE/ByMass/p2dRingCosThetaHEEVsMass";
                const std::string heeName = "Canvas_HEE";
                TProfile2D* heeNominal = cache.FetchClone<TProfile2D>(famFile, heeBase);
                if (!heeNominal) return; // Booked only with doFakePolDiagnosticsQA. FetchClone has already warned, and HEE is the last block here
                const std::string heeXTitle = heeNominal->GetXaxis()->GetTitle();
                const std::string heeYTitle = WithProxySubscript(heeNominal->GetZaxis()->GetTitle(), "LeadJet");
                TDirectory* heeDir = nullptr;

                // Is the HEE dependence a reflection of the AEE angle, or a separate effect?
                // Splitting on the sign of phi_Lambda - phi_p* answers it directly: if the cos(theta)
                // trend is driven by the AEE, the two halves separate and the inclusive curve sits
                // between them. If it survives with the same shape in both halves, it is its own effect.
                // README: this is the whole point of the AEE-split canvases below.
                auto drawAeeSplit = [&]() {
                    TProfile2D* pPos = cache.FetchClone<TProfile2D>(famFile, heeBase + "PosPhiAEE");
                    TProfile2D* pNeg = cache.FetchClone<TProfile2D>(famFile, heeBase + "NegPhiAEE");
                    const std::string splitWhat = fam.familyName + " HEE, AEE split";
                    if (!pPos || !pNeg || !HasMassRegionAxis(pPos, splitWhat) || !HasMassRegionAxis(pNeg, splitWhat)) {
                        delete pPos; delete pNeg;
                        return; // Each check above has already warned
                    }

                    if (!heeDir) heeDir = EnsureDir(EnsureDir(fOut, fam.familyName), "HEE");

                    // The inclusive curve is REBUILT from the two halves rather than read from the
                    // unsplit profile. TProfile::Add merges sums and entry counts, so the result is the
                    // exact entry-weighted inclusive -- and it is guaranteed to share binning with the
                    // two halves, which a separately booked profile is not.

                    auto sliceOf = [&](TProfile2D* p, int m, const char* tag) {
                        TProfile* s = p->ProfileX(Form("hee_%s_%d", tag, m), m, m);
                        if (s) s->SetDirectory(nullptr);
                        return s;
                    };

                    // (a) In-window only: all, positive half, negative half
                    TProfile* sPos = sliceOf(pPos, kMassPeak.bin, "pos");
                    TProfile* sNeg = sliceOf(pNeg, kMassPeak.bin, "neg");
                    if (sPos && sNeg) {
                        TProfile* sAll = (TProfile*)sPos->Clone("hee_all");
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

                        ProfileConfig pc; pc.xAxisTitle = heeXTitle; pc.yAxisTitle = heeYTitle;
                        const std::string nm = heeName + "_AEESplit_InPeak";
                        const std::string ttl = fam.familyName + " " + heeYTitle + " vs cos(#theta)_{HEE}, "
                                              + "in mass window, split by AEE angle sign";
                        DrawComparisonCanvas(hb, nm, ttl, heeDir, pc, heeXTitle, heeYTitle);
                        DrawComparisonCanvas(hb, nm + "_PullOnly", ttl, heeDir, pc, heeXTitle, heeYTitle, kPull);
                        delete sAll;
                    }
                    delete sPos; delete sNeg;

                    // (b) Everything: both halves across the three drawn mass regions on one canvas,
                    //     coloured like every other mass-region canvas
                    std::vector<ProfileBundle> hbAll;
                    std::vector<TH1*> hbDelete;
                    for (const MassRegion& region : kMassRegions) {
                        for (int side = 0; side < 2; ++side) {
                            TProfile2D* srcP = (side == 0) ? pPos : pNeg;
                            const std::string label = Form("%s, #Delta#phi_{AEE} %s", region.label, side == 0 ? "> 0" : "< 0");
                            TProfile* s = sliceOf(srcP, region.bin, side == 0 ? "P" : "N");
                            if (!s || s->GetEntries() < 1) {
                                WarnSkip(splitWhat + ", " + label, "no entries");
                                delete s;
                                continue;
                            }
                            VariationConfig c;
                            c.legendLabel = label;
                            c.color = region.style.color;
                            c.lineStyle = (side == 0) ? 1 : 2;   // Solid positive, dashed negative
                            c.markerStyle = (side == 0) ? 20 : 24; // Full / open, same shape
                            c.isData = (region.bin == kMassPeak.bin && side == 0);
                            hbAll.push_back({s, c});
                            hbDelete.push_back(s);
                        }
                    }
                    if (hbAll.size() >= 2) {
                        ProfileConfig pc; pc.xAxisTitle = heeXTitle; pc.yAxisTitle = heeYTitle;
                        DrawComparisonCanvas(hbAll, heeName + "_AEESplit_AllMass",
                                             fam.familyName + " " + heeYTitle + " vs cos(#theta)_{HEE}, "
                                             + "every mass region and AEE sign",
                                             heeDir, pc, heeXTitle, heeYTitle);
                    }
                    for (auto p : hbDelete) delete p;
                    delete pPos; delete pNeg;
                };
                drawAeeSplit();

                // Every mass region of the nominal file on one canvas
                if (!heeDir) heeDir = EnsureDir(EnsureDir(fOut, fam.familyName), "HEE");
                DrawMassSliceCanvas(heeNominal, heeName, fam.familyName + " " + heeYTitle + " vs cos(#theta)_{HEE}, by mass region",
                                    heeDir, heeXTitle, heeYTitle);

                // --- CrossSystem: the same HEE dependence across every consumer variation --------
                // This block read only the nominal file until now, which wasted the one thing the
                // variation machinery exists for. Same nomenclature as AEE_SignalExtracted's
                // CrossSystem, so the two read alike.
                //
                // The curves here are DIFFERENTIAL, not integrated single points, so the variations
                // are overlaid through the comparison drawer rather than placed on a categorical
                // axis. Only the in-peak mass slice is used: the sidebands are already covered by
                // Canvas_MassSelection, and overlaying nine variations times three mass regions
                // would be unreadable.
                {
                    std::vector<ProfileBundle> xb;
                    std::vector<TH1*> xbDel;

                    auto addVariation = [&](const std::string& suffix, const VariationConfig& vc) {
                        const std::string f = consumerDir + "/ConsumerResults_" + fam.dataSuffix +
                                              suffix + ".root";
                        TProfile2D* p = cache.FetchClone<TProfile2D>(f, heeBase);
                        if (!p) return; // FetchClone has already warned (or the file is a known absent variation)
                        const std::string what = fam.familyName + " HEE CrossSystem, " + vc.legendLabel;
                        if (!HasMassRegionAxis(p, what)) { delete p; return; }
                        const std::string tag = suffix.empty() ? "data" : SanitizeName(suffix);
                        TProfile* s = p->ProfileX(Form("heeX_%s", tag.c_str()), kMassPeak.bin, kMassPeak.bin);
                        delete p;
                        if (!s) return;
                        s->SetDirectory(nullptr);
                        if (s->GetEntries() < 1) { WarnSkip(what, "no entries in the mass window"); delete s; return; }
                        xb.push_back({s, vc});
                        xbDel.push_back(s);
                    };

                    addVariation("", dataConfig);
                    for (const auto& sys : sysVariations) addVariation(sys.suffix, sys);

                    if (xb.size() >= 2) {
                        TDirectory* xDir = EnsureDir(heeDir, "CrossSystem");

                        ProfileConfig pcx;
                        pcx.xAxisTitle = heeXTitle;
                        pcx.yAxisTitle = heeYTitle;
                        const std::string xt = fam.familyName + " " + heeYTitle + " vs cos(#theta)_{HEE}" +
                                               ", in mass window, across variations";

                        // Full and redux from the same call site, as everywhere else
                        DrawComparisonCanvas(xb, heeName, xt, xDir, pcx,
                                             heeXTitle, heeYTitle, kPull);
                        std::vector<ProfileBundle> xr = ReduxSubset(xb);
                        if (xr.size() >= 2)
                            DrawComparisonCanvas(xr, heeName + "_Redux", xt + " (reduced)",
                                                 xDir, pcx, heeXTitle, heeYTitle, kPull);
                    }
                    for (auto p : xbDel) delete p;
                }
                delete heeNominal;
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

                const std::string catPath = TASK_DIR_IN_FILE + catObs.path;
                const std::string dataFile = consumerDir + "/ConsumerResults_" + fam.dataSuffix + ".root";

                TProfile* pCatData = cache.FetchClone<TProfile>(dataFile, catPath);
                if (!pCatData) continue; // Not booked in this configuration (FetchClone has reported it): nothing is created on disk

                if (!famDir) famDir = EnsureDir(fOut, fam.familyName);
                TDirectory* etaStudyDir = EnsureDir(famDir, "EtaSplitStudies");
                TDirectory* catDir      = EnsureDir(etaStudyDir, ObservableDisplayName(catObs));
                TDirectory* catOverlay  = EnsureDir(catDir, "Overlays");
                TDirectory* catSubDir   = EnsureDir(catDir, "Subtracted");

                // The labelled axis carries the selection names, so the x title stays blank. The y title is the
                // data profile's, with the proxy as a subscript, as for the observable groups.
                ProfileConfig catConfig;
                catConfig.yAxisTitle = WithProxySubscript(pCatData->GetYaxis()->GetTitle(), ProxyOfPath(catObs.path));

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
                                     catSubDir, catConfig, "", "#Delta" + catConfig.yAxisTitle + " (Data - Var)",
                                     kNoLowerPad, true, true);
                DrawComparisonCanvas(catSubBundles, "Canvas_Subtracted_AllInOne_Pull",
                                     fam.familyName + " " + ObservableDisplayName(catObs) + " difference",
                                     catSubDir, catConfig, "", "#Delta" + catConfig.yAxisTitle + " (Data - Var)",
                                     kSelfPull, true, true);

                // --- Three-bin eta-sign digest, the one slice the Toy Model can speak to ---
                if (catObs.doEtaSignSummary) {
                    TDirectory* signDir = EnsureDir(catDir, "EtaSignSummary");

                    std::vector<ProfileBundle> signBundles;
                    std::vector<TH1*> signToDelete;
                    for (size_t i = 0; i < catBundles.size(); ++i) {
                        TH1D* h = ExtractEtaSignBins(catBundles[i].profile,
                                                     ObservableDisplayName(catObs) + "_sign_" + std::to_string(i));
                        if (!h) continue;
                        signBundles.push_back({h, catBundles[i].config});
                        signToDelete.push_back(h);
                    }

                    // Toy overlay: all three numbers or none, so a partial curve never misleads
                    for (const auto& ext : externals) {
                        if (ext.kind != ExternalRef::kToyModel || ext.basePath.empty()) continue;
                        TH1D* hToy = BuildToyEtaSignHist(cache, ext.basePath, pCatData,
                                                         ObservableDisplayName(catObs) + "_signToy");
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
                    const std::string massPath = TASK_DIR_IN_FILE + FolderOf(catObs.path) + catObs.massSplitProfile;

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
                            TProfile* proj = p2d->ProfileX((ObservableDisplayName(catObs) + "_" + massTag[m] + "_" +
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

            for (const auto& tableConfig : group.profiles) {
                const std::string profilePath = TASK_DIR_IN_FILE + tableConfig.path;

                // Fetch Data:
                std::string dataFile = consumerDir + "/ConsumerResults_" + fam.dataSuffix + ".root";
                TProfile* pData = cache.FetchClone<TProfile>(dataFile, profilePath);
                if (!pData) continue; // Data missing (FetchClone has reported it): next observable, nothing is created on disk

                // The axis titles come from the data profile, i.e. from the consumer's own booking, so they
                // cannot drift from it. Only the proxy is added, as a subscript of the y title: several x
                // titles (EtaV0, PVz) do not name it. Everything below reads this resolved copy.
                const std::string obsProxy = ProxyOfPath(tableConfig.path);
                ProfileConfig profConfig = tableConfig;
                if (profConfig.xAxisTitle.empty())
                    profConfig.xAxisTitle = std::string(pData->GetXaxis()->GetTitle()) + profConfig.xTitleSuffix;
                if (profConfig.yAxisTitle.empty())
                    profConfig.yAxisTitle = WithProxySubscript(pData->GetYaxis()->GetTitle(), obsProxy);

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

                // Only the variations present are kept, so allSystematics is NOT aligned with sysVariations.
                // sysIndexOf[k - 1] records which variation allSystematics[k] came from.
                std::vector<size_t> sysIndexOf;
                for (size_t iSys = 0; iSys < sysVariations.size(); ++iSys) {
                    const auto& sys = sysVariations[iSys];
                    std::string sysFile = consumerDir + "/ConsumerResults_" + fam.dataSuffix + sys.suffix + ".root";
                    TProfile* pSys = cache.FetchClone<TProfile>(sysFile, profilePath);
                    if (pSys) {
                        allSystematics.push_back({pSys, sys});
                        sysIndexOf.push_back(iSys);
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
                
                std::vector<size_t> massIndexOf; // Same as sysIndexOf: allMassSystematics[k] <-- massVariations[massIndexOf[k - 1]]
                for (size_t iMass = 0; iMass < massVariations.size(); ++iMass) {
                    const auto& mass = massVariations[iMass];
                    std::string massFile = consumerDir + "/ConsumerResults_" + fam.dataSuffix + mass.suffix + ".root";
                    TProfile* pMass = cache.FetchClone<TProfile>(massFile, profilePath);
                    if (pMass) {
                        allMassSystematics.push_back({pMass, mass});
                        massIndexOf.push_back(iMass);
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

                // Observables with inBruteForce = false stay off the brute-force axis. That is the eta-sign pT
                // splits: they cover different subsets of the sample, so their integrals sit at a
                // different scale and would flatten everything else on a shared axis.
                if (profConfig.inBruteForce) {
                    grandProxyTags.push_back(obsProxy);
                    grandDataVals.push_back(dataInteg);
                    grandLabels.push_back(ObservableDisplayName(profConfig));

                    // Each value goes to its OWN variation's row (sysIndexOf), never by position: allSystematics
                    // skips absent variations, and a positional match hands their rows to whoever comes next.
                    // An absent variation's row then stays short and makeSeries drops it, as it should.
                    for (size_t k = 1; k < allSystematics.size(); ++k)
                        grandSysVals[sysIndexOf[k - 1]].push_back(GetIntegratedProfile(dynamic_cast<TProfile*>(allSystematics[k].profile)));
                    for (const auto& ext : presentExternals) {
                        grandExtVals[ext.first].push_back(GetIntegratedProfile(dynamic_cast<TProfile*>(ext.second.profile)));
                    }
                    for (size_t k = 1; k < allMassSystematics.size(); ++k) // By index, as above
                        grandMassVals[massIndexOf[k - 1]].push_back(GetIntegratedProfile(dynamic_cast<TProfile*>(allMassSystematics[k].profile)));
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
                    if (profConfig.path != proxyRepresentatives[iRep].second) continue;
                    // Each value goes to its OWN variation's column (sysIndexOf), never by position: families need not
                    // share the same absent variations, and a positional match would put values under the wrong label
                    auto& slot = crossFamilyVals[iRep][famIdx];
                    slot.assign(crossFamilyLabels.size(), {0., -1.});
                    slot[0] = dataInteg;
                    for (size_t k = 1; k < allSystematics.size(); ++k)
                        slot[1 + sysIndexOf[k - 1]] = GetIntegratedProfile(dynamic_cast<TProfile*>(allSystematics[k].profile));
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

                        if (!fr.valid) {
                            WarnSkip(fitLabels[i], "tanh fit did not converge, left out of the p0 / p1 summaries",
                                     std::string(fitDir->GetPath()));
                            continue;
                        }
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
                if (!dataPt.found) { // Already reported by FetchIntegratedFromCuts
                    CountSkip("IntegratedSummary", "no Data value", fam.familyName + "/" + proxy.folderName);
                    continue;
                }
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
                            } else {
                                WarnSkip("Toy Model pRingProxyJet, bin 1", "no entries", ext.basePath);
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
                                                              proxy.folderName,
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

                if (labels.size() < 2) {
                    CountSkip("IntegratedSummary", "Data alone, nothing to compare it with", fam.familyName + "/" + proxy.folderName);
                    continue;
                }
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
                    if (!d.found) { // Already reported by FetchIntegratedByMassPeak
                        CountSkip("CheapSigExtract", "no Data value", fam.familyName + "/" + proxy.folderName);
                        continue;
                    }
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
                    if (labels.size() < 2) {
                        CountSkip("CheapSigExtract", "Data alone, nothing to compare it with", fam.familyName + "/" + proxy.folderName);
                        continue;
                    }

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

                // --- In against out, and their difference -------------------------------------
                // The two canvases above put in-peak and out-of-peak on separate pads, which makes
                // the comparison that actually matters -- how far the background pulls the number --
                // something the reader has to do by eye across two pictures. These put both states
                // on one axis per variation, and their difference on the next.
                //
                // The difference IS plain quadrature here, unlike the differences inside the
                // signal extraction: in-peak and out-of-peak candidates occupy disjoint mass
                // regions, so the two measurements share no candidate and are independent.
                {
                    std::vector<std::string> labels;
                    std::vector<VariationConfig> cfgs;
                    std::vector<SigExtractPoint> in, out;

                    auto collect = [&](const std::string& file, const VariationConfig& c) {
                        SigExtractPoint pin = FetchIntegratedByMassPeak(cache, file, proxy.folderName, true);
                        SigExtractPoint pout = FetchIntegratedByMassPeak(cache, file, proxy.folderName, false);
                        if (!pin.found || !pout.found) return; // Both states or neither
                        labels.push_back(c.legendLabel);
                        cfgs.push_back(c);
                        in.push_back(pin);
                        out.push_back(pout);
                    };

                    collect(dataFileCheap, dataConfig);
                    for (const auto& sys : sysVariations) {
                        if (!sys.inRedux) continue;
                        collect(consumerDir + "/ConsumerResults_" + fam.dataSuffix + sys.suffix + ".root", sys);
                    }

                    if (labels.size() >= 1) {
                        if (!cheapDir) cheapDir = EnsureDir(EnsureDir(fOut, fam.familyName), "CheapSigExtract");
                        TDirectory* pxDir = EnsureDir(cheapDir, proxy.folderName);

                        const int n = static_cast<int>(labels.size());
                        std::vector<ProfileBundle> bBoth, bDiff;
                        std::vector<TH1*> del2;

                        // Two series on one categorical axis: the in-peak state is the reference, so
                        // the lower pad reads as "how far out of peak sits from in peak".
                        TH1D* hIn = new TH1D(Form("CheapBoth_In_%s", proxy.folderName.c_str()), "", n, 0, n);
                        TH1D* hOut = new TH1D(Form("CheapBoth_Out_%s", proxy.folderName.c_str()), "", n, 0, n);
                        hIn->SetDirectory(nullptr);
                        hOut->SetDirectory(nullptr);
                        for (int k = 0; k < n; ++k) {
                            hIn->SetBinContent(k + 1, in[k].value);
                            hIn->SetBinError(k + 1, in[k].error);
                            hOut->SetBinContent(k + 1, out[k].value);
                            hOut->SetBinError(k + 1, out[k].error);

                            const double d = in[k].value - out[k].value;
                            const double e = std::sqrt(in[k].error * in[k].error +
                                                       out[k].error * out[k].error);
                            TH1D* hD = MakeCategoricalPoint(Form("CheapDiff_%s_%d",
                                                                 proxy.folderName.c_str(), k),
                                                            k, n, d, e);
                            bDiff.push_back({hD, cfgs[k]});
                            del2.push_back(hD);
                        }

                        VariationConfig cIn;
                        cIn.legendLabel = "In mass peak";
                        cIn.color = kTeal + 3;
                        cIn.lineStyle = 1;
                        cIn.markerStyle = 34;
                        cIn.isData = true;
                        VariationConfig cOut;
                        cOut.legendLabel = "Out of mass peak";
                        cOut.color = kMagenta + 2;
                        cOut.lineStyle = 1;
                        cOut.markerStyle = 23;
                        cOut.isData = false;
                        bBoth.push_back({hIn, cIn});
                        bBoth.push_back({hOut, cOut});
                        del2.push_back(hIn);
                        del2.push_back(hOut);

                        const std::string ttl2 = fam.familyName + " " + proxy.legendLabel +
                                                 ", integrated <R> in and out of the mass peak";
                        DrawOffsetSeriesCanvas(labels, bBoth, "Canvas_InVsOutOfMassPeak", ttl2,
                                               pxDir, "Integrated <R>", 0.5, true);
                        DrawIntegratedCanvas(bDiff, labels, "Canvas_InMinusOutOfMassPeak",
                                             ttl2 + ", difference", pxDir,
                                             "<R>_{in peak} - <R>_{out of peak}", true, false, true);
                        for (auto p : del2) delete p;
                    }
                }
            }
        }

        // -----------------------------------------------------------------------------------
        // Signal extraction of the AEE angle: differential, angle-combined, and cross-system
        // -----------------------------------------------------------------------------------
        // These live at the TOP LEVEL of the extraction file, outside its Ring/ folder, because their
        // consumer source (FakePolDrivers/PhiAEE/) sits outside Ring/ too. That is the reason a path
        // built like the per-proxy IntegratedSummary ones would find nothing.
        //
        // The layout is proxy-major then split-major: everything about one proxy is under a single
        // folder, so a question about the leading jet never means opening three sibling folders in
        // turn. The eta split rides in the canvas NAME instead, which puts its three variants side by
        // side in one listing -- which is how they are actually read. The species is the family's: every
        // extraction file carries one set per proxy, named after the proxy alone, read in every family.
        //
        // hRSigMinusRBkg is read, never recomputed: <R>_S and <R>_B share the sideband primitives
        // and quadrature UNDERSTATES their difference by roughly a fifth. The same correlation is
        // why the decomposition canvases below carry no lower pad.
        if (!sigExtractDir.empty()) {
            const std::string sigFile = sigExtractDir + "/signalExtractionRing_" + fam.dataSuffix + ".root";
            const std::string aeeRoot = "HelicityEfficiencyQA_PhiLambdaPhiProtonStar/";
            // A consumer path, so it carries the task folder (the extraction paths above do not)
            const std::string aeeConsumerRoot = TASK_DIR_IN_FILE + "FakePolDrivers/PhiAEE/";

            // The extraction writes three eta-split variants, each filled identically. Mirrored here
            // one-for-one, so a name in the summary output maps to a folder in the extraction output
            // with no translation.
            struct EtaSplitView { const char* dir; const char* label; int color; int marker; };
            const std::vector<EtaSplitView> etaSplits = {
                {"NoEtaSplit",  "all #eta_{proxy}", kBlack,    20},
                {"EtaProxyPos", "#eta_{proxy} > 0", kRed + 1,  21},
                {"EtaProxyNeg", "#eta_{proxy} < 0", kBlue + 1, 22}
            };
            TDirectory* aeeSigDir = nullptr;

            // The full-range curve comes from the consumer's <R> vs AEE angle vs proxy-eta-sign profile. It has no mass
            // axis at all, so it is the one measured curve that cannot lose candidates to a mass
            // underflow/overflow bin, and it shares axisDeltaPhiCoarse with everything the extraction produced.
            // Safer if we want to look at the full, uncut, observable!
            // name is the proxy: the extraction's folder and object suffix, and this macro's output folder. The
            // labels name only the proxy too: the family name in front of them names the species.
            struct AeeSpec { std::string name, label; };
            const std::vector<AeeSpec> aeeSpecs = {
                {"LeadJet", "Leading jet"},
                {"LeadP",   "Leading particle"},
                {"SubJet",  "Subleading jet"}
            };

            // Axis titles are not written here as part of a DRY philosophy.
            // The consumer's booking string is the only place this axis should be named, minimizing the maintenance surface for me later on.
            // The extraction code repeats it onto everything it produces, so the first object that carries one "wins" and nothing in this macro can drift from the consumer.
            auto axisTitleOf = [](const std::vector<TH1*>& candidates) {
                for (TH1* h : candidates) {
                    if (!h || !h->GetXaxis()) continue;
                    const TString t = h->GetXaxis()->GetTitle();
                    if (t.Length() > 0) return std::string(t.Data());
                }
                return std::string("");
            };

            // The measured <R> over the FULL mass range, split by proxy eta, straight from the consumer and untouched by the extraction:
            // The Y axis is the proxy eta SIGN, two bins filled at -0.5 and +0.5 they are addressed by VALUE rather than by index so the
            // convention stays where the fill put it and a later reordering of that axis cannot silently swap the halves.
            auto measuredFullRange = [&](const std::string& consumerFile, const std::string& obj,
                                         const std::string& splitDir, const std::string& uniq) -> TH1D* {
                TProfile2D* p2 = cache.FetchClone<TProfile2D>(consumerFile, aeeConsumerRoot + obj, false);
                if (!p2) return nullptr;
                const int bNeg = p2->GetYaxis()->FindBin(-0.5);
                const int bPos = p2->GetYaxis()->FindBin(+0.5);
                int lo = std::min(bNeg, bPos), hi = std::max(bNeg, bPos);
                if (splitDir == "EtaProxyPos") { lo = bPos; hi = bPos; }
                else if (splitDir == "EtaProxyNeg") { lo = bNeg; hi = bNeg; }

                TH1D* h = nullptr;
                TProfile* pr = p2->ProfileX((uniq + "_prof").c_str(), lo, hi);
                if (pr) {
                    // ProjectionX of a TProfile carries the mean in the content and the error OF THE
                    // MEAN in the error, which is what every other <R> curve on these canvases is.
                    h = pr->ProjectionX(uniq.c_str(), "e");
                    if (h) {
                        h->SetDirectory(nullptr);
                        h->GetXaxis()->SetTitle(p2->GetXaxis()->GetTitle());
                    }
                    delete pr;
                }
                delete p2;
                return h;
            };

            for (const auto& spec : aeeSpecs) {
                // One folder per proxy:
                // (created lazily, on the first split that has anything in it, so a family that never produced this spec leaves no empty folder)
                TDirectory* specDir = nullptr;   // <Proxy>/
                TDirectory* angleDir = nullptr;  // <Proxy>/AEEangle/
                TDirectory* combDir = nullptr;   // <Proxy>/Combined/
                TDirectory* combQaDir = nullptr; // <Proxy>/Combined/QA/
                TDirectory* xAngDir = nullptr;   // <Proxy>/AEEangle/CrossSystem_Angular/
                TDirectory* xCatDir = nullptr;   // <Proxy>/Combined/CrossSystem_Categorical/

                auto ensureSpecDirs = [&]() {
                    if (specDir) return;
                    if (!aeeSigDir)
                        aeeSigDir = EnsureDir(EnsureDir(fOut, fam.familyName), "AEE_SignalExtracted");
                    specDir = EnsureDir(aeeSigDir, spec.name);
                    angleDir = EnsureDir(specDir, "AEEangle");
                    combDir = EnsureDir(specDir, "Combined");
                    combQaDir = EnsureDir(combDir, "QA");
                };

              for (const auto& esp : etaSplits) {
                const std::string split = esp.dir;
                const std::string base = aeeRoot + split + "/" + spec.name + "/";
                TH1D* hSig = cache.FetchClone<TH1D>(sigFile, base + "Results/hRSig_" + spec.name, false);
                TH1D* hBkg = cache.FetchClone<TH1D>(sigFile, base + "Results/hRBkg_" + spec.name, false);
                TH1D* hMeas = cache.FetchClone<TH1D>(sigFile, base + "Results/hRMeas_" + spec.name, false);
                TH1D* hDiff = cache.FetchClone<TH1D>(sigFile, base + "Results/hRSigMinusRBkg_" + spec.name, false);
                if (!hSig || !hBkg) { delete hSig; delete hBkg; delete hMeas; delete hDiff; continue; }

                TH1D* hMeasFull = measuredFullRange(
                    consumerDir + "/ConsumerResults_" + fam.dataSuffix + ".root",
                    spec.name + "/p2dRingObservableVsPhiAEEVsProxyEtaSign", split,
                    "measFull_" + fam.familyName + "_" + spec.name + "_" + split);

                ensureSpecDirs();

                // Colour carries the QUANTITY on these canvases, because the wagon is fixed. On the
                // cross-system canvases further down it is the other way round.
                VariationConfig cS;
                cS.legendLabel = "<R>_{S}";
                cS.color = kRed + 1;
                cS.lineStyle = 1;
                cS.markerStyle = 20;
                cS.isData = true;
                VariationConfig cB;
                cB.legendLabel = "<R>_{B}";
                cB.color = kBlue + 1;
                cB.lineStyle = 1;
                cB.markerStyle = 21;
                cB.isData = false;
                VariationConfig cM;
                cM.color = kBlack;
                cM.lineStyle = 1;
                cM.markerStyle = 24;
                cM.isData = false;

                ProfileConfig pc;
                pc.xAxisTitle = axisTitleOf({hSig, hBkg, hMeasFull});
                pc.yAxisTitle = "R";

                const std::string t = fam.familyName + " " + spec.label + ", signal extracted, " +
                                      esp.label;

                //// The decomposition, two mass treatments: ////
                // Both carry <R>_S and <R>_B unchanged and differ only in what "measured" means, so
                // the pair isolates the cost of the mass window and the binning rather than mixing
                // it into a single curve nobody can then unpick.
                //
                // NO LOWER PAD, on either. The three curves are functions of the same four sideband
                // primitives, so the quadrature that kDifference and kPull both perform would
                // overstate Var(<R>_meas - <R>_S) by roughly 60 percent. The one difference with a
                // correctly propagated error is hRSigMinusRBkg, which gets its own canvas below.
                {
                    VariationConfig cMp = cM;
                    cMp.legendLabel = "<R>_{meas}^{binnedInMass}";
                    std::vector<ProfileBundle> b = {{hSig, cS}, {hBkg, cB}};
                    if (hMeas) b.push_back({hMeas, cMp});
                    // The identity <R>_meas = p <R>_S + (1-p) <R>_B holds bin by bin in the peak
                    // window, so this canvas is also a closure test: a measured curve that does not
                    // sit between the other two, weighted by purity, is a defect and not a
                    // fluctuation.
                    DrawComparisonCanvas(b, "Canvas_Decomposition_PeakWindow_" + split,
                                         t + ", peak-window decomposition", angleDir, pc,
                                         pc.xAxisTitle, pc.yAxisTitle);
                }
                if (hMeasFull) {
                    VariationConfig cMf = cM;
                    cMf.legendLabel = "<R>_{meas}^{fullMass}";
                    // This measured curve obeys no identity with the other two: it counts every
                    // candidate at every mass, while they describe the peak window only. It is here
                    // because it is the one curve in the picture that no fit, window or bin
                    // rejection has touched.
                    std::vector<ProfileBundle> b = {{hSig, cS}, {hBkg, cB}, {hMeasFull, cMf}};
                    DrawComparisonCanvas(b, "Canvas_Decomposition_FullMass_" + split,
                                         t + ", full-mass-range measured", angleDir, pc,
                                         pc.xAxisTitle, pc.yAxisTitle);
                }

                if (hDiff) {
                    VariationConfig cD;
                    cD.legendLabel = "<R>_{S} - <R>_{B}";
                    cD.color = kGreen + 3;
                    cD.lineStyle = 1;
                    cD.markerStyle = 22;
                    cD.isData = true;
                    std::vector<ProfileBundle> bd = {{hDiff, cD}};
                    DrawComparisonCanvas(bd, "Canvas_SigMinusBkg_" + split, t + ", difference",
                                         angleDir, pc, pc.xAxisTitle, "<R>_{S} - <R>_{B}");
                }

                //// Cross-system, differential ////
                // Curve against curve across wagons, as opposed to the categorical block below,
                // which compares one integrated number per wagon. Both are wanted: a wagon can agree
                // on the integral and disagree everywhere along the angle.
                //
                // Here the lower pad IS valid. Different wagons are separate consumer runs, so the
                // quadrature kDifference performs is the right error for the comparison -- the same
                // argument the systematics canvases rest on. (The mixed-event estimator reuses the
                // same V0 sample, so its pad over-covers slightly; that is accepted here exactly as
                // it is accepted everywhere else in this macro rather than special-cased.)
                {
                    struct XQuantity { const char* obj; const char* canvas; const char* legend;
                                       const char* yTitle; int color; int marker; };
                    const std::vector<XQuantity> xq = {
                        {"hRSig", "RSig", "<R>_{S}", "<R>_{S}", kRed + 1, 20},
                        {"hRBkg", "RBkg", "<R>_{B}", "<R>_{B}", kBlue + 1, 21},
                        {"hRMeas", "RMeas", "<R>_{meas}^{binnedInMass}", "<R>_{meas}", kBlack, 24},
                        {"hRSigMinusRBkg", "RSigMinusRBkg", "<R>_{S} - <R>_{B}",
                         "<R>_{S} - <R>_{B}", kGreen + 3, 22}
                    };

                    // Not every systematic: the two pairings the eye actually compares. Adding a
                    // wagon here multiplies the canvas count by the number of quantities, which is
                    // why the exhaustive sweep stays in the categorical block.
                    struct XWagonSet { const char* dirTag; const char* title;
                                       std::vector<std::string> suffixes; };
                    const std::vector<XWagonSet> wagonSets = {
                        {"DataVsMixedEv", "data vs mixed-event", {"", "_MixedEventProxies", "_MixedEventProxies" + MIXED_RESAMPLES_TAG}},
                        {"ANCuts", "analysis-level cuts",
                         {"", "_analysisCuts", "_MixedEventProxies_analysisCuts", "_MixedEventProxies_analysisCuts" + MIXED_RESAMPLES_TAG}}
                    };

                    for (const auto& ws : wagonSets) {
                        // Fetched once and used twice: the per-quantity canvases and the overview
                        // are built from the same pointers, so the two views cannot end up showing
                        // different things under the same wagon label.
                        struct XWagon { VariationConfig cfg; std::map<std::string, TH1D*> byQty; };
                        std::vector<XWagon> wagons;
                        std::vector<TH1*> ownX;

                        for (const auto& sfx : ws.suffixes) {
                            // Styling comes from the variation table, so these canvases cannot drift
                            // from the colours used everywhere else in the file.
                            VariationConfig vc = dataConfig;
                            if (!sfx.empty()) {
                                bool known = false;
                                for (const auto& sys : sysVariations)
                                    if (sys.suffix == sfx) { vc = sys; known = true; break; }
                                if (!known) {
                                    WarnSkip("wagon " + sfx, "not in the variation table (sysVariations)", ws.dirTag);
                                    continue;
                                }
                            }
                            vc.isData = sfx.empty(); // The uncut data is the reference

                            XWagon w;
                            w.cfg = vc;
                            const std::string wf = sigExtractDir + "/signalExtractionRing_" +
                                                   fam.dataSuffix + sfx + ".root";
                            for (const auto& q : xq) {
                                TH1D* h = cache.FetchClone<TH1D>(
                                    wf, aeeRoot + split + "/" + spec.name + "/Results/" +
                                        std::string(q.obj) + "_" + spec.name, false);
                                if (!h) continue;
                                w.byQty[q.obj] = h;
                                ownX.push_back(h);
                            }
                            if (!w.byQty.empty()) wagons.push_back(w);
                        }

                        if (wagons.size() < 2) {
                            CountSkip("AEE_SignalExtracted, CrossSystem_Angular", "fewer than two wagons extracted",
                                      fam.familyName + "/" + spec.name + "/" + ws.dirTag + "/" + split);
                            for (auto p : ownX) delete p;
                            continue;
                        }
                        if (!xAngDir) xAngDir = EnsureDir(angleDir, "CrossSystem_Angular");
                        TDirectory* wsDir = EnsureDir(xAngDir, ws.dirTag);

                        // Per quantity: colour is the WAGON, since the quantity is fixed here.
                        for (const auto& q : xq) {
                            std::vector<ProfileBundle> perQ;
                            for (const auto& w : wagons) {
                                auto it = w.byQty.find(q.obj);
                                if (it == w.byQty.end()) continue;
                                perQ.push_back({it->second, w.cfg});
                            }
                            if (perQ.size() < 2) continue;
                            DrawComparisonCanvas(perQ,
                                                 "Canvas_" + std::string(q.canvas) + "_" + split,
                                                 t + ", " + q.legend + " across wagons (" +
                                                     ws.title + ")",
                                                 wsDir, pc, pc.xAxisTitle, q.yTitle, kDifference);
                        }

                        // Everything at once: colour is the QUANTITY and the line style is the
                        // wagon, so a reader follows one quantity across wagons by colour and one
                        // wagon across quantities by dash pattern. Busy on purpose -- it is the QA
                        // view for the canvases above, not the one to read a number from, and it
                        // carries no lower pad because its curves are different quantities with no
                        // common reference to subtract.
                        std::vector<ProfileBundle> overview;
                        for (size_t iw = 0; iw < wagons.size(); ++iw) {
                            for (const auto& q : xq) {
                                auto it = wagons[iw].byQty.find(q.obj);
                                if (it == wagons[iw].byQty.end()) continue;
                                VariationConfig ov;
                                ov.legendLabel = std::string(q.legend) + ", " +
                                                 wagons[iw].cfg.legendLabel;
                                ov.color = q.color;
                                ov.lineStyle = static_cast<int>(iw) + 1;
                                ov.markerStyle = q.marker;
                                ov.isData = false;
                                overview.push_back({it->second, ov});
                            }
                        }
                        if (overview.size() >= 2)
                            DrawComparisonCanvas(overview, "Canvas_All_" + split,
                                                 t + ", all quantities, " + ws.title, wsDir, pc,
                                                 pc.xAxisTitle, pc.yAxisTitle, kNoLowerPad);

                        for (auto p : ownX) delete p;
                    }
                }

                //// The angle-combined digests ////
                // Not differential anymore, so they sit in Combined/ rather than AEEangle/. 
                // Part of the digest is drawn on SEPARATE canvases on purpose: purity and significance are of order 0.1 and order 100~1000 while
                // every <R> is of order 1e-2, so a shared axis would flatten the ring observables onto zero.
                auto drawDigest = [&](const std::string& obj, const std::string& canvas,
                                      const std::string& title, const std::string& yTitle,
                                      TDirectory* dest, const std::string& note = "") {
                    TH1D* h = cache.FetchClone<TH1D>(sigFile, base + "IntegratedCombined/" + obj +
                                                              "_" + spec.name, false);
                    if (!h) return;
                    VariationConfig c;
                    c.legendLabel = "Angle-combined";
                    c.color = kBlack;
                    c.lineStyle = 1;
                    c.markerStyle = 20;
                    c.isData = true;
                    std::vector<ProfileBundle> bc = {{h, c}};
                    // The note travels in the title, so it survives into the canvas rather than
                    // living only in whatever documentation the reader did not open.
                    const std::string fullTitle = note.empty() ? title : (title + " [" + note + "]");
                    DrawComparisonCanvas(bc, canvas, fullTitle, dest, pc, "", yTitle,
                                         kNoLowerPad, true);
                    delete h;
                };
                drawDigest("hCombinedSummary", "Canvas_Combined_" + split,
                           t + ", angle-combined results", "Integrated <R>", combDir);
                drawDigest("hCombinedQA", "Canvas_CombinedQA_" + split,
                           t + ", cross-checks and acceptance diagnostics", "Value", combQaDir,
                           "\"weighted-flat\" = integrating WITHOUT counts as weights");
                drawDigest("hCombinedQuality", "Canvas_CombinedQuality_" + split,
                           t + ", extraction quality", "Value", combQaDir);

                //// Cross-system, categorical ////
                // One integrated number per wagon, swept over every systematic in the table. The
                // differential counterpart lives in AEEangle/CrossSystem_Angular/ and deliberately
                // covers fewer wagons; the two are named apart because a TBrowser listing of
                // "CrossSystem" holding both was a standing source of confusion.
                // README: this couples to the row order of hCombinedSummary in signalExtractionRing.
                {
                    struct XRow { int bin; const char* name; const char* title; const char* yTitle; };
                    // hCombinedSummary row order, which this is coupled to:
                    //   1 <R>_meas^FullMassRange   2 <R>_meas^PeakWin,AccBins
                    //   3 <R>_S                    4 <R>_B                    5 <R>_S - <R>_B
                    const std::vector<XRow> xrows = {
                        {3, "Canvas_CrossSystem_RSig",  "<R>_{S} across variations", "<R>_{S}"},
                        {4, "Canvas_CrossSystem_RBkg",  "<R>_{B} across variations", "<R>_{B}"},
                        {5, "Canvas_CrossSystem_RSigMinusRBkg",
                            "<R>_{S} - <R>_{B} across variations", "<R>_{S} - <R>_{B}"}
                    };

                    auto fetchRow = [&](const std::string& suffix, int bin) {
                        SigExtractPoint p;
                        const std::string f = sigExtractDir + "/signalExtractionRing_" + suffix + ".root";
                        TH1D* h = cache.FetchClone<TH1D>(
                            f, aeeRoot + split + "/" + spec.name +
                               "/IntegratedCombined/hCombinedSummary_" + spec.name, false);
                        if (!h) return p;
                        if (h->GetNbinsX() >= bin && h->GetBinError(bin) > 0.0) {
                            p.value = h->GetBinContent(bin);
                            p.error = h->GetBinError(bin);
                            p.found = true;
                        }
                        delete h;
                        return p;
                    };

                    for (const auto& xr : xrows) {
                        std::vector<std::string> xLabels;
                        std::vector<VariationConfig> xCfgs;
                        std::vector<std::pair<double, double>> xVals;

                        SigExtractPoint dp = fetchRow(fam.dataSuffix, xr.bin);
                        if (!dp.found) {
                            CountSkip("AEE_SignalExtracted, CrossSystem_Categorical", "no Data value",
                                      fam.familyName + "/" + spec.name + "/" + xr.name + "_" + split);
                            continue;
                        }
                        xLabels.push_back(dataConfig.legendLabel);
                        xCfgs.push_back(dataConfig);
                        xVals.push_back({dp.value, dp.error});

                        for (const auto& sys : sysVariations) {
                            SigExtractPoint p = fetchRow(fam.dataSuffix + sys.suffix, xr.bin);
                            if (!p.found) continue;
                            xLabels.push_back(sys.legendLabel);
                            xCfgs.push_back(sys);
                            xVals.push_back({p.value, p.error});
                        }
                        if (xLabels.size() < 2) {
                            CountSkip("AEE_SignalExtracted, CrossSystem_Categorical", "Data alone, nothing to compare it with",
                                      fam.familyName + "/" + spec.name + "/" + xr.name + "_" + split);
                            continue;
                        }

                        if (!xCatDir) xCatDir = EnsureDir(combDir, "CrossSystem_Categorical");

                        // The split is in the canvas name, not in a folder: the three variants then
                        // list side by side, which is how they get compared.
                        DrawCategoricalCrossSystem(xLabels, xCfgs, xVals, std::string(xr.name) + "_" + split,
                                                   fam.familyName + " " + spec.label + " " + xr.title + ", " + esp.label,
                                                   xCatDir, xr.yTitle);
                    }
                }

                delete hSig; delete hBkg; delete hMeas; delete hDiff; delete hMeasFull;
              } // eta-split variants

                // --- The eta splits on one canvas -------------------------------------------
                // Same content as each variant's Canvas_Combined, overlaid, so "does the result
                // depend on which side of the detector the proxy came from?" is one picture rather
                // than three canvases opened in turn. A split-dependent result points at acceptance
                // rather than at physics.
                //
                // Built through a lambda so the same four-curve set can be produced for any consumer
                // variation. Data alone, Data against MixedEv, and the analysis-cut wagons all reuse
                // it rather than each assembling their own and drifting apart.
                if (specDir) {
                    // Returns {NoEtaSplit, EtaProxyPos, EtaProxyNeg, recombined} for one wagon.
                    auto buildEtaSet = [&](const std::string& suffix, const std::string& tag,
                                           int lineStyle, std::vector<TH1*>& own) {
                        std::vector<ProfileBundle> out;
                        std::map<std::string, TH1*> byName;
                        const std::string f = sigExtractDir + "/signalExtractionRing_" +
                                              fam.dataSuffix + suffix + ".root";

                        for (const auto& esp : etaSplits) {
                            TH1D* h = cache.FetchClone<TH1D>(
                                f, aeeRoot + esp.dir + "/" + spec.name +
                                   "/IntegratedCombined/hCombinedSummary_" + spec.name, false);
                            if (!h) continue;
                            VariationConfig c;
                            c.legendLabel = tag.empty() ? esp.label
                                                        : (tag + ", " + std::string(esp.label));
                            c.color = esp.color;
                            c.lineStyle = lineStyle;
                            c.markerStyle = esp.marker;
                            c.isData = (std::string(esp.dir) == "NoEtaSplit" && lineStyle == 1);
                            out.push_back({h, c});
                            own.push_back(h);
                            byName[esp.dir] = h;
                        }

                        // The recombined point. See the per-row weighting note below.
                        // hCombinedQuality bins: 1 purity, 2 significance, 3 signal, 4 background,
                        // 5 bins used, 6 bins total, 7 full-range counts.
                        struct HalfYields { double sig = 1.0, bkg = 1.0, tot = 2.0, full = 2.0; };
                        auto yieldsOf = [&](const char* dir) {
                            HalfYields y;
                            TH1D* hq = cache.FetchClone<TH1D>(
                                f, aeeRoot + std::string(dir) + "/" + spec.name +
                                   "/IntegratedCombined/hCombinedQuality_" + spec.name, false);
                            if (hq) {
                                if (hq->GetNbinsX() >= 4) {
                                    if (hq->GetBinContent(3) > 0.0) y.sig = hq->GetBinContent(3);
                                    if (hq->GetBinContent(4) > 0.0) y.bkg = hq->GetBinContent(4);
                                    y.tot = y.sig + y.bkg;
                                    y.full = y.tot;
                                }
                                if (hq->GetNbinsX() >= 7 && hq->GetBinContent(7) > 0.0)
                                    y.full = hq->GetBinContent(7);
                                delete hq;
                            }
                            return y;
                        };

                        auto itP = byName.find("EtaProxyPos");
                        auto itN = byName.find("EtaProxyNeg");
                        TH1* hPos = (itP != byName.end()) ? itP->second : nullptr;
                        TH1* hNeg = (itN != byName.end()) ? itN->second : nullptr;
                        if (hPos && hNeg && hPos->GetNbinsX() == hNeg->GetNbinsX()) {
                            const HalfYields yP = yieldsOf("EtaProxyPos");
                            const HalfYields yN = yieldsOf("EtaProxyNeg");

                            TH1D* hAvg = (TH1D*)hPos->Clone(Form("EtaRecomb_%s_%s",
                                                                 spec.name.c_str(),
                                                                 suffix.empty() ? "data"
                                                                                : SanitizeName(suffix).c_str()));
                            hAvg->SetDirectory(nullptr);
                            for (int ib = 1; ib <= hAvg->GetNbinsX(); ++ib) {
                                // Every row of hCombinedSummary is a mean over a DIFFERENT population, so every row gets its own weight.
                                // Row order:
                                //   1 <R>_meas^FullMassRange  -> all candidates, all mass  -> full
                                //   2 <R>_meas^PeakWin,AccBins-> peak window, kept bins    -> T=S+B
                                //   3 <R>_S                                                -> S
                                //   4 <R>_B                                                -> B
                                //   5 <R>_S - <R>_B                                        -> S
                                //
                                // CLOSURE: row 1 involves no fit, no mass window and no bin
                                // rejection, so recombining it MUST reproduce NoEtaSplit. If it does
                                // not, the weights are wrong. Row 2 need not close, and the gap
                                // between rows 1 and 2 is exactly what it cost to restrict the mass
                                // and drop the failed bins.
                                double nP = yP.full, nN = yN.full;
                                if (ib == 2) { nP = yP.tot; nN = yN.tot; }
                                else if (ib == 3) { nP = yP.sig; nN = yN.sig; }
                                else if (ib == 4) { nP = yP.bkg; nN = yN.bkg; }
                                else if (ib == 5) {
                                    // Row 5 is <R>_S - <R>_B, and it is THE result of this whole
                                    // construction, not a leftover diagnostic.
                                    //
                                    // <R>_S is not the observable on its own: a true Lambda is
                                    // subject to the same acceptance, AEE and HEE distortions as a
                                    // combinatorial one, so <R>_S carries the fake signal too. The
                                    // background, by assumption, carries ONLY the fake signal --
                                    // there is no physics polarization in a combinatorial pair. The
                                    // difference is therefore what survives once the fake part is
                                    // removed from a sample that is genuinely made of Lambdas.
                                    //
                                    // The subtraction inside each half is already done in the
                                    // extraction, with the correlated error: <R>_S and <R>_B share
                                    // the sideband primitives there, and quadrature would understate
                                    // it by about a fifth. So this recombines the stored per-half
                                    // DIFFERENCES rather than subtracting two separately recombined
                                    // rows -- which would have been wrong twice over, once for
                                    // mixing weight sets and once for dropping the correlation.
                                    //
                                    // Weighted by SIGNAL yield: the quantity being combined is the
                                    // physics polarization of the signal, so a half with more
                                    // Lambdas carries more of it. The fake part is assumed COMMON to
                                    // both populations rather than proportional to either yield,
                                    // which is exactly why it cancels in the difference before this
                                    // weighting is ever applied. Across halves the candidates are
                                    // disjoint, so quadrature here is exact.
                                    nP = yP.sig;
                                    nN = yN.sig;
                                }
                                const double wP = nP / (nP + nN);
                                const double wN = nN / (nP + nN);
                                const double vP = hPos->GetBinContent(ib), eP = hPos->GetBinError(ib);
                                const double vN = hNeg->GetBinContent(ib), eN = hNeg->GetBinError(ib);
                                hAvg->SetBinContent(ib, wP * vP + wN * vN);
                                hAvg->SetBinError(ib, std::sqrt(wP * wP * eP * eP + wN * wN * eN * eN));
                            }
                            VariationConfig cA;
                            cA.legendLabel = (tag.empty() ? std::string("") : tag + ", ") +
                                             "Pos #oplus Neg (recombined)";
                            cA.color = kGreen + 3;
                            cA.lineStyle = lineStyle;
                            cA.markerStyle = 33;
                            cA.isData = false;
                            out.push_back({hAvg, cA});
                            own.push_back(hAvg);
                        }
                        return out;
                    };

                    ProfileConfig pcc;
                    pcc.xAxisTitle = "";
                    pcc.yAxisTitle = "Integrated <R>";
                    const std::string ctBase = fam.familyName + " " + spec.label +
                                               ", angle-combined results by proxy #eta sign";

                    // Reduced view: only the two curves that carry the message. The gap between the
                    // unsplit result and the recombined one is the quantity of interest, so the
                    // question this canvas asks is whether a second estimator shows the same gap.
                    auto reducedOf = [](const std::vector<ProfileBundle>& full) {
                        std::vector<ProfileBundle> r;
                        for (const auto& b : full)
                            if (b.config.legendLabel.find("all #eta_{proxy}") != std::string::npos ||
                                b.config.legendLabel.find("recombined") != std::string::npos)
                                r.push_back(b);
                        return r;
                    };

                    // One pair of canvases per wagon pairing: the uncut data against the mixed-event
                    // estimator, and the same two with the analysis-level cuts applied.
                    struct EtaPairing { std::string dataSuffix, mixedSuffix, canvasSuffix, dirName; };
                    const std::vector<EtaPairing> pairings = {
                        {"", "_MixedEventProxies", "", ""},
                        {"_analysisCuts", "_MixedEventProxies_analysisCuts", "_ANCuts", "ANCuts"}
                    };

                    for (const auto& pr : pairings) {
                        std::vector<TH1*> own;
                        std::vector<ProfileBundle> dataSet =
                            buildEtaSet(pr.dataSuffix, "", 1, own);
                        std::vector<ProfileBundle> mixedSet =
                            buildEtaSet(pr.mixedSuffix, "MixedEv", 2, own); // Dashed, same colours

                        // ANCuts/ holds two different comparisons and they are easy to confuse from
                        // a TBrowser listing, so they go in named subfolders rather than side by
                        // side: these compare the eta SPLITS for one wagon, while the per-split
                        // canvases further down compare the WAGONS for one split. Both now sit under
                        // this spec's folder rather than once at the top of AEE_SignalExtracted.
                        TDirectory* dest = pr.dirName.empty()
                                           ? EnsureDir(specDir, "EtaComparison")
                                           : EnsureDir(EnsureDir(specDir, "ANCuts"), "BySplit");

                        if (dataSet.size() >= 2) {
                            DrawComparisonCanvas(dataSet,
                                                 "Canvas_EtaComparison" + pr.canvasSuffix,
                                                 ctBase, dest, pcc, "", "Integrated <R>",
                                                 kNoLowerPad, true);
                        }

                        // Eight curves: every eta variant for both estimators. Busy on purpose --
                        // it is the QA view for the reduced canvas below, not the one to read from.
                        std::vector<ProfileBundle> both = dataSet;
                        both.insert(both.end(), mixedSet.begin(), mixedSet.end());
                        if (!mixedSet.empty() && both.size() >= 2) {
                            DrawComparisonCanvas(both,
                                                 "Canvas_EtaComparison_WithMixedEv" + pr.canvasSuffix,
                                                 ctBase + ", data and mixed-event", dest, pcc, "",
                                                 "Integrated <R>", kNoLowerPad, true);

                            std::vector<ProfileBundle> red = reducedOf(dataSet);
                            std::vector<ProfileBundle> redM = reducedOf(mixedSet);
                            red.insert(red.end(), redM.begin(), redM.end());
                            if (red.size() >= 2)
                                DrawComparisonCanvas(red,
                                                     "Canvas_EtaComparison_Reduced" + pr.canvasSuffix,
                                                     ctBase + ", unsplit vs recombined", dest, pcc, "",
                                                     "Integrated <R>", kNoLowerPad, true);
                        }
                        for (auto p : own) delete p;
                    }
                }

                // --- Analysis-level cuts, per eta-split mode -----------------------------------
                // The _analysisCuts wagons re-run the consumer with isV0Accepted enabled, so they
                // answer a different question from the proxy variations: not "does an artificial
                // proxy reproduce this", but "does the result survive a tighter selection", and
                // whether the background shrank relative to the signal. One canvas per eta-split
                // mode, so the comparison is against the SAME split rather than across splits.
                //
                // Complementary to the ANCuts canvases produced above: those compare the SPLITS for
                // a given wagon, these compare the WAGONS for a given split.
                if (specDir) {
                    const std::vector<std::string> cutVariations = {
                        "", "_analysisCuts", "_MixedEventProxies_analysisCuts", "_MixedEventProxies_analysisCuts" + MIXED_RESAMPLES_TAG
                    };

                    for (const auto& esp : etaSplits) {
                        std::vector<ProfileBundle> cb;
                        std::vector<TH1*> cbDel;

                        for (const auto& suffix : cutVariations) {
                            // Find the variation's styling in the table, so these canvases cannot
                            // drift from the colours used everywhere else.
                            VariationConfig vc = dataConfig;
                            if (!suffix.empty()) {
                                bool known = false;
                                for (const auto& sys : sysVariations)
                                    if (sys.suffix == suffix) { vc = sys; known = true; break; }
                                if (!known) {
                                    WarnSkip("wagon " + suffix, "not in the variation table (sysVariations)", "cutVariations");
                                    continue;
                                }
                            }

                            TH1D* h = cache.FetchClone<TH1D>(
                                sigExtractDir + "/signalExtractionRing_" + fam.dataSuffix + suffix + ".root",
                                aeeRoot + esp.dir + "/" + spec.name +
                                    "/IntegratedCombined/hCombinedSummary_" + spec.name, false);
                            if (!h) continue;
                            vc.isData = suffix.empty(); // The uncut data is the reference
                            cb.push_back({h, vc});
                            cbDel.push_back(h);
                        }

                        if (cb.size() >= 2) {
                            ProfileConfig pcc;
                            pcc.xAxisTitle = "";
                            pcc.yAxisTitle = "Integrated <R>";
                            DrawComparisonCanvas(cb, std::string("Canvas_") + esp.dir +
                                                     "_AnalysisCuts",
                                                 fam.familyName + " " + spec.label + ", " +
                                                     esp.label + ", analysis-level cuts",
                                                 EnsureDir(EnsureDir(specDir, "ANCuts"), "ByWagon"),
                                                 pcc, "", "Integrated <R>", kNoLowerPad, true);
                        }
                        for (auto p : cbDel) delete p;
                    }
                }

                //// The eta_Proxy-binned extraction ////
                // Everything here reads the Summary/ folder the extraction's Step 2.8 writes: each
                // quantity once, on the consumer's own axes, so no binning is rebuilt here and no
                // slice folder is walked.
                //
                // There is no full-mass measured curve in this set, on purpose. The only mass-free
                // source, p2d...VsProxyEta, is binned on an eta axis whose edges never meet the
                // slices', and a curve on another binning is not comparable cell by cell. The
                // FullMassRange row of each slice is used once, in the recombination closure, where it
                // is the one row that must close.
                {
                    const std::string etaSum = "HelicityEfficiencyQA_PhiLambdaPhiProtonStar_EtaBinned/" +
                                               spec.name + "/Summary/";
                    auto sigFileOf = [&](const std::string& sfx) {
                        return sigExtractDir + "/signalExtractionRing_" + fam.dataSuffix + sfx + ".root";
                    };
                    auto fetchMap = [&](const std::string& sfx, const std::string& q) {
                        return cache.FetchClone<TH2D>(sigFileOf(sfx), etaSum + "hMap_" + q, false);
                    };
                    auto fetchVsEta = [&](const std::string& sfx, const std::string& q) {
                        return cache.FetchClone<TH1D>(sigFileOf(sfx), etaSum + "hVsEta_" + q, false);
                    };

                    TH2D* hStatus = fetchMap("", "CellStatus");
                    if (hStatus) {
                        ensureSpecDirs();
                        TDirectory* etaDir = EnsureDir(specDir, "EtaBinned");
                        const TAxis* etaAxis = hStatus->GetYaxis();
                        const double etaMax = std::max(std::fabs(etaAxis->GetXmin()), std::fabs(etaAxis->GetXmax()));
                        const std::string tEta = fam.familyName + " " + spec.label + ", #eta_{Proxy}-binned";

                        //// Maps over (angle, eta_Proxy) ////
                        // The slope significance says whether THIS cell's background is measurably
                        // tilted. It is not a flatness score to compare across binnings: at a fixed
                        // true slope it scales with sqrt(N), so finer cells look flatter for purely
                        // statistical reasons.
                        TDirectory* mapDir = EnsureDir(etaDir, "Maps");
                        struct MapView { const char* q; bool isSigned; };
                        const std::vector<MapView> mapViews = {
                            {"RSig", true}, {"RBkg", true}, {"RMeas", true}, {"RSigMinusRBkg", true},
                            {"Purity", false}, {"SigStat", false}, {"SigYield", false},
                            {"RingBkgSlopeSignif", true}, {"RingBkgChi2Ndf", false}
                        };
                        for (const auto& mv : mapViews) {
                            TH2D* h = fetchMap("", mv.q);
                            if (!h) continue;
                            DrawMapCanvas(h, hStatus, std::string("Canvas_Map_") + mv.q,
                                          tEta + ", " + h->GetTitle() + " (hatched: not extracted)",
                                          mapDir, mv.isSigned);
                            delete h;
                        }
                        // The status map itself, codes printed: which cells failed, and why
                        DrawMapCanvas(hStatus, nullptr, "Canvas_Map_CellStatus",
                                      tEta + ", cell status (1 extracted, 2 few entries, 3 peak fit, "
                                             "4 sidebands, 5 bkg fit, 6 no signal; 0 slice skipped)",
                                      mapDir, false, true);

                        //// Every eta_Proxy bin against the angle, one quantity per canvas ////
                        // No lower pad: there is no reference slice, and the message is the pattern
                        // across all of them rather than any one difference.
                        TDirectory* angDir = EnsureDir(etaDir, "AEEangle");
                        struct AngView { const char* q; const char* yTitle; };
                        const std::vector<AngView> angViews = {
                            {"RSig", "<R>_{S}"}, {"RBkg", "<R>_{B}"},
                            {"RMeas", "<R>_{meas}^{binnedInMass}"}, {"RSigMinusRBkg", "<R>_{S} - <R>_{B}"}
                        };
                        for (const auto& av : angViews) {
                            TH2D* h = fetchMap("", av.q);
                            if (!h) continue;
                            std::vector<ProfileBundle> bs;
                            std::vector<TH1*> own;
                            for (int k = 1; k <= h->GetNbinsY(); ++k) {
                                TH1D* row = h->ProjectionX(Form("%s_%s_eta%02d", spec.name.c_str(), av.q, k), k, k, "e");
                                row->SetDirectory(nullptr);
                                own.push_back(row);
                                bool any = false;
                                for (int ix = 1; ix <= row->GetNbinsX() && !any; ++ix)
                                    any = (row->GetBinError(ix) > 0.0);
                                if (!any) { // A slice with no extracted cell (see Canvas_Map_CellStatus): no legend entry for nothing
                                    CountSkip("AEE_SignalExtracted, eta_Proxy slice", "no extracted cell, slice not drawn",
                                              fam.familyName + "/" + spec.name + "/" + av.q);
                                    continue;
                                }

                                const double lo = etaAxis->GetBinLowEdge(k), hi = etaAxis->GetBinUpEdge(k);
                                VariationConfig c;
                                c.legendLabel = Form("#eta_{Proxy} #in [%+.2f, %+.2f)", lo, hi);
                                c.color = EtaSliceColor(0.5 * (lo + hi), etaMax);
                                c.lineStyle = 1;
                                c.markerStyle = (0.5 * (lo + hi) < 0.0) ? 24 : 20;
                                c.isData = false;
                                bs.push_back({row, c});
                            }
                            if (!bs.empty()) {
                                ProfileConfig pca;
                                pca.xAxisTitle = h->GetXaxis()->GetTitle();
                                pca.yAxisTitle = av.yTitle;
                                DrawComparisonCanvas(bs, std::string("Canvas_") + av.q + "_AllEtaBins",
                                                     tEta + ", " + av.yTitle + " in every #eta_{Proxy} bin",
                                                     angDir, pca, pca.xAxisTitle, av.yTitle);
                            }
                            for (auto p : own) delete p;
                            delete h;
                        }

                        //// Angle-combined, against eta_Proxy ////
                        TDirectory* depDir = EnsureDir(etaDir, "EtaDependence");

                        // One row of the Step 2.75 angle-combined result for an eta-sign split
                        auto combinedRow = [&](const std::string& sfx, const std::string& split, int row) {
                            SigExtractPoint p;
                            TH1D* h = cache.FetchClone<TH1D>(sigFileOf(sfx), aeeRoot + split + "/" + spec.name +
                                                                 "/IntegratedCombined/hCombinedSummary_" + spec.name, false);
                            if (!h) return p;
                            if (h->GetNbinsX() >= row && h->GetBinError(row) > 0.0) {
                                p.value = h->GetBinContent(row);
                                p.error = h->GetBinError(row);
                                p.found = true;
                            }
                            delete h;
                            return p;
                        };
                        // A coarser result repeated over the eta_Proxy bins it covers, so it overlays the
                        // binned curve as a reference. Only where the binned curve has a point, so a jet
                        // proxy's structurally empty edge bins stay empty.
                        auto referenceOver = [&](TH1D* shape, const std::string& name, int row, bool bySign) -> TH1D* {
                            const SigExtractPoint all = combinedRow("", "NoEtaSplit", row);
                            const SigExtractPoint pos = combinedRow("", "EtaProxyPos", row);
                            const SigExtractPoint neg = combinedRow("", "EtaProxyNeg", row);
                            TH1D* r = (TH1D*)shape->Clone(name.c_str());
                            r->SetDirectory(nullptr);
                            r->Reset();
                            bool any = false;
                            for (int k = 1; k <= r->GetNbinsX(); ++k) {
                                if (shape->GetBinError(k) <= 0.0) continue;
                                const SigExtractPoint& p = !bySign ? all : (r->GetBinCenter(k) < 0.0 ? neg : pos);
                                if (!p.found) continue;
                                r->SetBinContent(k, p.value);
                                r->SetBinError(k, p.error);
                                any = true;
                            }
                            if (!any) { delete r; return nullptr; }
                            return r;
                        };

                        ProfileConfig pce;
                        pce.xAxisTitle = etaAxis->GetTitle();
                        pce.yAxisTitle = "Angle-combined <R>";

                        // The decomposition, as on the angular canvases. No lower pad, for the same
                        // reason: the three are functions of the same sideband primitives.
                        {
                            struct Part { const char* q; const char* legend; int color; int marker; bool isData; };
                            const std::vector<Part> parts = {
                                {"RSig", "<R>_{S}", kRed + 1, 20, true},
                                {"RBkg", "<R>_{B}", kBlue + 1, 21, false},
                                {"RMeasPeakWindow", "<R>_{meas}^{binnedInMass}", kBlack, 24, false}
                            };
                            std::vector<ProfileBundle> b;
                            std::vector<TH1*> own;
                            for (const auto& pt : parts) {
                                TH1D* h = fetchVsEta("", pt.q);
                                if (!h) continue;
                                VariationConfig c;
                                c.legendLabel = pt.legend;
                                c.color = pt.color;
                                c.lineStyle = 1;
                                c.markerStyle = pt.marker;
                                c.isData = pt.isData;
                                b.push_back({h, c});
                                own.push_back(h);
                            }
                            if (!b.empty())
                                DrawComparisonCanvas(b, "Canvas_Decomposition_VsEta",
                                                     tEta + ", angle-combined decomposition", depDir, pce,
                                                     pce.xAxisTitle, pce.yAxisTitle);
                            for (auto p : own) delete p;
                        }

                        // Each result against the coarser splits it refines. Correlated with them --
                        // the same candidates -- so no lower pad; the question is whether the binned
                        // curve averages to the sign-split points, not how far apart they are.
                        struct DepView { const char* q; int row; const char* yTitle; };
                        const std::vector<DepView> depViews = {
                            {"RSig", 3, "<R>_{S}"}, {"RBkg", 4, "<R>_{B}"}, {"RSigMinusRBkg", 5, "<R>_{S} - <R>_{B}"}
                        };
                        for (const auto& dv : depViews) {
                            TH1D* hB = fetchVsEta("", dv.q);
                            if (!hB) continue;
                            TH1D* hSgn = referenceOver(hB, spec.name + "_refSign_" + dv.q, dv.row, true);
                            TH1D* hAll = referenceOver(hB, spec.name + "_refAll_" + dv.q, dv.row, false);

                            VariationConfig cB;
                            cB.legendLabel = "#eta_{Proxy}-binned";
                            cB.color = kBlack;
                            cB.lineStyle = 1;
                            cB.markerStyle = 20;
                            cB.isData = true;
                            VariationConfig cS;
                            cS.legendLabel = "#eta_{Proxy} sign split, repeated per half";
                            cS.color = kAzure + 2;
                            cS.lineStyle = 2;
                            cS.markerStyle = 25;
                            cS.isData = false;
                            VariationConfig cA;
                            cA.legendLabel = "no #eta_{Proxy} split, repeated";
                            cA.color = kGray + 2;
                            cA.lineStyle = 3;
                            cA.markerStyle = 27;
                            cA.isData = false;

                            std::vector<ProfileBundle> b = {{hB, cB}};
                            if (hSgn) b.push_back({hSgn, cS});
                            if (hAll) b.push_back({hAll, cA});
                            DrawComparisonCanvas(b, std::string("Canvas_") + dv.q + "_VsEta",
                                                 tEta + ", " + dv.yTitle + " against the coarser splits",
                                                 depDir, pce, pce.xAxisTitle, dv.yTitle);
                            delete hB; delete hSgn; delete hAll;
                        }

                        //// The eta_Proxy slices recombined ////
                        // One weight per row, because each row of hCombinedSummary is a mean over a
                        // different population -- the same argument as the Pos/Neg recombination above,
                        // with N pieces instead of two. Row order, coupled to the extraction:
                        //   1 <R>_meas^FullMassRange  -> full-range counts
                        //   2 <R>_meas^PeakWin        -> S + B
                        //   3 <R>_S                   -> S
                        //   4 <R>_B                   -> B
                        //   5 <R>_S - <R>_B           -> S (the fake part cancels inside each slice, not between slices)
                        // The slices are disjoint, so quadrature across them is (should be) exact.
                        auto recombine = [&](const std::string& sfx) -> TH1D* {
                            // The unsplit result supplies the labelled five-row axis, and is the reference the closure (of this recombination) compares against
                            // (see hRef below)
                            TH1D* out = cache.FetchClone<TH1D>(sigFileOf(sfx), aeeRoot + "NoEtaSplit/" + spec.name +
                                                                   "/IntegratedCombined/hCombinedSummary_" + spec.name, false);
                            TH1D* hN = fetchVsEta(sfx, "FullRangeCounts");
                            TH1D* hS = fetchVsEta(sfx, "SigYield");
                            TH1D* hB = fetchVsEta(sfx, "BkgYield");
                            if (!out || !hN || !hS || !hB) { delete out; delete hN; delete hS; delete hB; return nullptr; }
                            out->SetName(Form("EtaBinnedRecomb_%s_%s", spec.name.c_str(),
                                              sfx.empty() ? "data" : SanitizeName(sfx).c_str()));
                            out->Reset();

                            const char* rowQ[5] = {"RMeasFullMass", "RMeasPeakWindow", "RSig", "RBkg", "RSigMinusRBkg"};
                            for (int r = 1; r <= std::min(5, out->GetNbinsX()); ++r) {
                                TH1D* hv = fetchVsEta(sfx, rowQ[r - 1]);
                                if (!hv) continue;
                                double sw = 0.0, swv = 0.0, sw2e2 = 0.0;
                                for (int k = 1; k <= hv->GetNbinsX(); ++k) {
                                    const double e = hv->GetBinError(k);
                                    if (!(e > 0.0)) continue; // Not valid in that slice
                                    const double w = (r == 1) ? hN->GetBinContent(k)
                                                   : (r == 2) ? hS->GetBinContent(k) + hB->GetBinContent(k)
                                                   : (r == 4) ? hB->GetBinContent(k)
                                                              : hS->GetBinContent(k);
                                    if (!(w > 0.0)) continue;
                                    sw += w;
                                    swv += w * hv->GetBinContent(k);
                                    sw2e2 += w * w * e * e;
                                }
                                if (sw > 0.0) {
                                    out->SetBinContent(r, swv / sw);
                                    out->SetBinError(r, std::sqrt(sw2e2) / sw);
                                }
                                delete hv;
                            }
                            delete hN; delete hS; delete hB;
                            return out;
                        };

                        // Closure against the unsplit extraction. Only row 1 must close, and only up to
                        // the candidates outside the eta_Proxy axis (reported by the extraction as
                        // hSliceClosure): it involves no fit, no window and no bin rejection. The other
                        // rows need not, and their gap is what the finer split changed.
                        {
                            TH1D* hRef = cache.FetchClone<TH1D>(sigFileOf(""), aeeRoot + "NoEtaSplit/" + spec.name +
                                                                    "/IntegratedCombined/hCombinedSummary_" + spec.name, false);
                            TH1D* hRec = recombine("");
                            if (hRef && hRec) {
                                VariationConfig cRef;
                                cRef.legendLabel = "no #eta_{Proxy} split";
                                cRef.color = kBlack;
                                cRef.lineStyle = 1;
                                cRef.markerStyle = 20;
                                cRef.isData = true;
                                VariationConfig cRec;
                                cRec.legendLabel = "#eta_{Proxy} bins, recombined";
                                cRec.color = kGreen + 3;
                                cRec.lineStyle = 1;
                                cRec.markerStyle = 33;
                                cRec.isData = false;
                                ProfileConfig pcc;
                                pcc.xAxisTitle = "";
                                pcc.yAxisTitle = "Integrated <R>";
                                std::vector<ProfileBundle> b = {{hRef, cRef}, {hRec, cRec}};
                                DrawComparisonCanvas(b, "Canvas_RecombinationClosure",
                                                     tEta + ", recombined against the unsplit result (row 1 must close)",
                                                     depDir, pcc, "", "Integrated <R>", kNoLowerPad, true);
                            }
                            delete hRef; delete hRec;
                        }

                        //// Data against the mixed-event estimator, differential in eta_Proxy ////
                        // Separate consumer runs, so here the lower pad is a genuine Data - MixedEv. Once per
                        // mixed-event wagon: single borrow, then resampled (canvas name and title tagged).
                        struct MixWagon { std::string suffix, nameTag, titleTag; };
                        const std::vector<MixWagon> mixWagons = {
                            {"_MixedEventProxies", "", ""},
                            {"_MixedEventProxies" + MIXED_RESAMPLES_TAG, "_" + resampCount + "Resamp", " (" + resampCount + " resamp.)"}};
                        for (const auto& mxw : mixWagons) {
                            const VariationConfig* mix = nullptr;
                            for (const auto& sys : sysVariations)
                                if (sys.suffix == mxw.suffix) { mix = &sys; break; }
                            if (mix) {
                                TDirectory* dmDir = nullptr;
                                struct MixView { const char* q; const char* yTitle; };
                                const std::vector<MixView> mixViews = {
                                    {"RSig", "<R>_{S}"}, {"RBkg", "<R>_{B}"},
                                    {"RMeasPeakWindow", "<R>_{meas}^{binnedInMass}"},
                                    {"RSigMinusRBkg", "<R>_{S} - <R>_{B}"}
                                };
                                for (const auto& mw : mixViews) {
                                    TH1D* hD = fetchVsEta("", mw.q);
                                    TH1D* hM = fetchVsEta(mix->suffix, mw.q);
                                    if (hD && hM) {
                                        if (!dmDir) dmDir = EnsureDir(depDir, "CrossSystem_DataVsMixedEv");
                                        VariationConfig cD = dataConfig;
                                        cD.isData = true;
                                        VariationConfig cM = *mix;
                                        cM.isData = false;
                                        std::vector<ProfileBundle> b = {{hD, cD}, {hM, cM}};
                                        DrawComparisonCanvas(b, std::string("Canvas_") + mw.q + "_VsEta" + mxw.nameTag,
                                                             tEta + ", " + mw.yTitle + ", data vs mixed-event" + mxw.titleTag,
                                                             dmDir, pce, pce.xAxisTitle, mw.yTitle, kDifference);
                                    }
                                    delete hD; delete hM;
                                }
                            }
                        }

                        //// Cross-system of the final result: the recombined numbers ////
                        // Every systematic in the table, each recombined with its own yields, so a
                        // variation that moves the signal between eta_Proxy bins is weighted as it
                        // should be rather than with the data's weights.
                        {
                            struct Recomb { std::string label; VariationConfig cfg; TH1D* h; };
                            std::vector<Recomb> rec;
                            if (TH1D* h = recombine("")) rec.push_back({dataConfig.legendLabel, dataConfig, h});
                            if (!rec.empty()) {
                                for (const auto& sys : sysVariations)
                                    if (TH1D* h = recombine(sys.suffix)) rec.push_back({sys.legendLabel, sys, h});

                                struct XRow { int bin; const char* name; const char* title; const char* yTitle; };
                                const std::vector<XRow> xrows = {
                                    {3, "Canvas_CrossSystem_RSig", "<R>_{S} across variations", "<R>_{S}"},
                                    {4, "Canvas_CrossSystem_RBkg", "<R>_{B} across variations", "<R>_{B}"},
                                    {5, "Canvas_CrossSystem_RSigMinusRBkg", "<R>_{S} - <R>_{B} across variations",
                                     "<R>_{S} - <R>_{B}"}
                                };
                                TDirectory* xDir = nullptr;
                                for (const auto& xr : xrows) {
                                    std::vector<std::string> labels;
                                    std::vector<VariationConfig> cfgs;
                                    std::vector<std::pair<double, double>> vals;
                                    for (const auto& r : rec) {
                                        if (r.h->GetNbinsX() < xr.bin || !(r.h->GetBinError(xr.bin) > 0.0)) {
                                            if (labels.empty()) break; // No data point, nothing to compare against
                                            continue;
                                        }
                                        labels.push_back(r.label);
                                        cfgs.push_back(r.cfg);
                                        vals.push_back({r.h->GetBinContent(xr.bin), r.h->GetBinError(xr.bin)});
                                    }
                                    if (labels.size() < 2) {
                                        CountSkip("AEE_SignalExtracted, recombined CrossSystem", "Data alone, nothing to compare it with",
                                                  fam.familyName + "/" + spec.name + "/" + xr.name);
                                        continue;
                                    }
                                    if (!xDir) xDir = EnsureDir(etaDir, "CrossSystem_Categorical");
                                    DrawCategoricalCrossSystem(labels, cfgs, vals, xr.name,
                                                               tEta + " recombined, " + xr.title, xDir, xr.yTitle);
                                }
                            }
                            for (auto& r : rec) delete r.h;
                        }
                        delete hStatus;
                    }
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
                                                             proxy.folderName,
                                                             "hIntegratedRSig");
                if (!probe.found) { // Already reported by FetchSigExtractPoint
                    CountSkip("SignalExtraction", "no extracted <R>_S", fam.familyName + "/" + proxy.folderName);
                    continue;
                }

                if (!famDir) famDir = EnsureDir(fOut, fam.familyName);
                if (!sigExtractGroupDir) sigExtractGroupDir = EnsureDir(famDir, "SignalExtraction");
                TDirectory* proxyDir = EnsureDir(sigExtractGroupDir, proxy.folderName);

                // (a) The three components side by side, plus their correlated differences
                DrawSigExtractComponents(cache, sigExtractDir, fam.dataSuffix,
                                         proxy, fam.familyName,
                                         EnsureDir(proxyDir, "Components"));

                // (b) One systematics canvas pair per quantity, imitating Integrated_summary
                for (size_t iQ = 0; iQ < sigExtractQuantities.size(); ++iQ) {
                    const auto& quantity = sigExtractQuantities[iQ];
                    SigExtractPoint dataVal = DrawSigExtractSystematics(
                        cache, sigExtractDir, fam.dataSuffix, sysVariations, dataConfig,
                        proxy, quantity, fam.familyName,
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
                                                                 fam.dataSuffix,
                                                                 proxy.folderName,
                                                                 sigExtractQuantities[iQ].histName));
                }
                for (const auto& sys : sysVariations) {
                    std::vector<SigExtractPoint> row;
                    bool complete = true;
                    for (const auto& quantity : sigExtractQuantities) {
                        SigExtractPoint p = FetchSigExtractPoint(cache, sigExtractDir,
                                                                 fam.dataSuffix + sys.suffix,
                                                                 proxy.folderName,
                                                                 quantity.histName);
                        if (!p.found) { complete = false; break; }
                        row.push_back(p);
                    }
                    if (!complete) { // Keep the two series on a common axis
                        CountSkip("SignalExtraction, <R>_S and <R>_B series", "variation without both quantities",
                                  fam.familyName + "/" + proxy.folderName + "/" + sys.legendLabel);
                        continue;
                    }
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
            auto proxyRank = [&](const std::string& tag) {
                for (size_t r = 0; r < PROXIES.size(); ++r) if (PROXIES[r] == tag) return r;
                return PROXIES.size(); // Anything unrecognised sorts last rather than being dropped
            };
            std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
                return proxyRank(grandProxyTags[a]) < proxyRank(grandProxyTags[b]);
            });

            // Rebuild the labels in the new order. The display names carry the proxy, so the blocks are self-evident
            std::vector<std::string> sortedLabels;
            for (size_t k : order) sortedLabels.push_back(grandLabels[k]);

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

            // One series per variation, or none: a variation absent for some observables would misalign the axis
            auto seriesOf = [&](const std::vector<std::pair<double, double>>& raw, const std::string& name,
                                const std::string& label) -> TH1D* {
                if (!raw.empty() && raw.size() != order.size())
                    CountSkip("BruteForce, " + label, "absent for some observables, so its series is left out", fam.familyName);
                return makeSeries(reorder(raw), name);
            };

            // (a) Values
            std::vector<ProfileBundle> grandBundles;
            TH1D* hGrandData = seriesOf(grandDataVals, "GrandData", dataConfig.legendLabel);
            if (hGrandData) { grandBundles.push_back({hGrandData, dataConfig}); grandToDelete.push_back(hGrandData); }

            for (size_t sysIdx = 0; sysIdx < sysVariations.size(); ++sysIdx) {
                TH1D* h = seriesOf(grandSysVals[sysIdx], Form("GrandSys_%zu", sysIdx), sysVariations[sysIdx].legendLabel);
                if (h) { grandBundles.push_back({h, sysVariations[sysIdx]}); grandToDelete.push_back(h); }
            }
            for (size_t extIdx = 0; extIdx < externals.size(); ++extIdx) {
                if (externals[extIdx].basePath.empty()) continue;
                TH1D* h = seriesOf(grandExtVals[extIdx], Form("GrandExt_%zu", extIdx), externals[extIdx].config.legendLabel);
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
                TH1D* h = seriesOf(grandMassVals[mIdx], Form("GrandMass_%zu", mIdx), massVariations[mIdx].legendLabel);
                if (h) { grandMassBundles.push_back({h, massVariations[mIdx]}); grandToDelete.push_back(h); }
            }
            if (grandMassBundles.size() > 1) {
                DrawIntegratedCanvas(grandMassBundles, sortedLabels, "Canvas_MassSelection",
                                     fam.familyName + " integrated R, mass windows", integDir,
                                     "Integrated R", false, true);
            }

            for (auto p : grandToDelete) delete p;
        }

        // A longitudinal family is relabelled once complete, with its own symbol; nothing else is touched
        if (fam.ring) RelabelRingDirectory(fOut->GetDirectory(fam.familyName.c_str()), fam.ring->symbolSuffix);
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
            const std::string proxyName = proxyRepresentatives[iRep].first;

            // The families holding this observable, with a Data value to subtract from
            std::vector<size_t> heldBy;
            for (size_t iFam = 0; iFam < families.size(); ++iFam) {
                if (families[iFam].ring) continue; // Full ring only: the longitudinal rings are compared in RingVsLongitudinal*/
                const auto& vals = crossFamilyVals[iRep][iFam];
                if (vals.empty() || vals[0].second < 0.) {
                    CountSkip("CrossFamily_Integrated", "family without the representative observable",
                              families[iFam].familyName + "/" + proxyName);
                    continue;
                }
                heldBy.push_back(iFam);
            }
            // The columns: Data, then every variation at least one of those families holds
            std::vector<size_t> cols;
            for (size_t c = 0; c < crossFamilyLabels.size(); ++c)
                for (size_t iFam : heldBy)
                    if (crossFamilyVals[iRep][iFam][c].second >= 0.) { cols.push_back(c); break; }
            std::vector<std::string> labels;
            for (size_t c : cols) labels.push_back(crossFamilyLabels[c]);
            const int nCols = static_cast<int>(cols.size());

            std::vector<ProfileBundle> famBundles;
            std::vector<ProfileBundle> famSubBundles;
            std::vector<TH1*> famToDelete;

            for (size_t iFam : heldBy) {
                const auto& vals = crossFamilyVals[iRep][iFam];

                // familyStyles is indexed in the same order as families, so the legend names the selection
                const VariationConfig& style = familyStyles[std::min(iFam, familyStyles.size() - 1)];

                // Difference against THIS family's own data, which is bin 1 of its own series. Subtracting
                // each selection from its own baseline is what makes the three comparable: an offset that
                // is common to a selection cancels, so what is left is how differently each systematic
                // moves Lambda, AntiLambda and BothHyperons.
                TH1D* h = new TH1D(Form("XFam_%zu_%zu", iRep, iFam), "", nCols, 0, nCols);
                TH1D* hSub = new TH1D(Form("XFamSub_%zu_%zu", iRep, iFam), "", nCols, 0, nCols);
                h->SetDirectory(nullptr);
                hSub->SetDirectory(nullptr);
                const double d  = vals[0].first;
                const double ed = vals[0].second;
                for (int i = 0; i < nCols; ++i) {
                    const auto& v = vals[cols[i]];
                    if (v.second < 0.) continue; // Not in this family: a gap, never a zero
                    h->SetBinContent(i + 1, v.first);
                    h->SetBinError(i + 1, v.second);
                    if (cols[i] == 0) continue; // Bin 1 (the data itself) is left empty on the difference
                    hSub->SetBinContent(i + 1, d - v.first);
                    hSub->SetBinError(i + 1, std::sqrt(ed * ed + v.second * v.second));
                }
                famBundles.push_back({h, style});
                famSubBundles.push_back({hSub, style});
                famToDelete.push_back(h);
                famToDelete.push_back(hSub);
            }

            if (famBundles.size() < 2)
                CountSkip("CrossFamily_Integrated", "fewer than two families hold the representative observable", proxyName);
            if (famBundles.size() > 1) {
                if (!xFamDir) xFamDir = EnsureDir(fOut, "CrossFamily_Integrated");

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
                if (families[f].ring) continue; // Full ring only, as above
                bool anyFound = false;
                for (size_t iQ = 0; iQ < sigExtractQuantities.size(); ++iQ)
                    if (sigExtractCrossFamily[iQ][iProxy][f].found) anyFound = true;
                if (!anyFound) continue;
                famLabels.push_back(families[f].familyName);
                famIndices.push_back(static_cast<int>(f));
            }
            if (famLabels.size() < 2) { // A single family is not a cross-family plot
                CountSkip("CrossFamily_Integrated/SigExtracted", "fewer than two families with an extraction", sigExtractProxies[iProxy].folderName);
                continue;
            }

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

    // Cross-definition folders, from the full-ring families: each pairs a file with its longitudinal siblings and uses
    // only the files it finds, so a wagon without longitudinal outputs still gets its full-ring corrections.
    std::cout << " -> Full ring vs the longitudinal rings...\n";
    DrawRingVsLongitudinal(cache, fOut, consumerDir, fullFamilies, dataConfig, sysVariations, externalsAll);
    DrawCheapRingVsLongitudinal(cache, fOut, consumerDir, fullFamilies, dataConfig, sysVariations, externalsAll);
    const std::vector<FakePairing> pairings = MakeFakePairings(consumerDir, sigExtractDir, externalsAll);
    std::cout << " -> Full ring vs the longitudinal rings, differential...\n";
    DrawDifferentialRingVsLongitudinal(cache, fOut, consumerDir, fullFamilies, dataConfig, sysVariations, externalsAll, pairings);

    // Raw: every candidate, from the consumer's mass-binned moments
    std::cout << " -> Corrections, raw (no signal extraction)...\n";
    const CorrectionFetcher fetchRaw = [&](const FakePairing& pr, const std::string& sfx, const std::string& proxy,
                                           CorrectionInputs& in) {
        RingSample smp;
        // LoadRingSample already refuses a sample with fewer than two candidates, the one way InputsOfSample can fail
        if (!LoadRingSample(cache, pr.consumerBase + "/ConsumerResults_" + sfx + ".root", proxy, smp)) return false;
        return InputsOfSample(smp, in);
    };
    DrawCorrections(fOut, fullFamilies, pairings, "Raw", "raw, all candidates", fetchRaw);
    DrawDifferentialCorrections(cache, fOut, fullFamilies, pairings, "Raw", "raw, all candidates",
                                [](const DifferentialSample& d, size_t i, CorrectionInputs& in) { return InputsOfSample(d.all[i], in); });

    // Signal-extracted: <R>_S from IntegratedSummary/, <u v>_S and <w>_S from KappaEff/ (signalExtractionRing).
    // <w>_S is extracted for longitudinal files only: for the full ring w == 1 identically. kappa = 3<u v>_S/<w>_S
    // takes its two errors in quadrature; they are positively correlated (both scale as n_z^2), so this
    // overstates it slightly, on a term that only enters P through dP/dkappa = -P/D anyway.
    if (!sigExtractDir.empty()) {
        std::cout << " -> Corrections, signal-extracted...\n";
        const CorrectionFetcher fetchSig = [&](const FakePairing& pr, const std::string& sfx, const std::string& proxy,
                                               CorrectionInputs& in) {
            const std::string sigFileOfSample = pr.sigBase + "/signalExtractionRing_" + sfx + ".root";
            if (pr.sigBase.empty() || !cache.Open(sigFileOfSample, false))
                return false; // Not extracted: only counted, it is the common case for references
            const SigExtractPoint r = FetchSigExtractPoint(cache, pr.sigBase, sfx, proxy, "hIntegratedRSig");
            const SigExtractPoint u = FetchSigExtractPoint(cache, pr.sigBase, sfx, proxy + "_Num", "hIntegratedRSig", true, "KappaEff");
            if (!r.found || !u.found) return false; // Already counted or warned by FetchSigExtractPoint
            if (!(u.value > 0.)) {
                WarnSkip(proxy + " signal-extracted <u v>", "not positive", sigFileOfSample);
                return false;
            }
            SigExtractPoint w; w.value = 1.; w.error = 0.; w.found = true;
            if (IsLongitudinalSuffix(sfx))
                w = FetchSigExtractPoint(cache, pr.sigBase, sfx, proxy + "_Den", "hIntegratedRSig", true, "KappaEff");
            if (!w.found) return false; // Already counted or warned by FetchSigExtractPoint
            if (!(w.value > 0.)) {
                WarnSkip(proxy + " signal-extracted <w>", "not positive", sigFileOfSample);
                return false;
            }
            in.ring = {r.value, r.error};
            in.uvx3 = {3. * u.value, 3. * u.error};
            const double kap = 3. * u.value / w.value;
            in.kappa = {kap, kap * std::sqrt(std::pow(u.error / u.value, 2) + std::pow(w.error / w.value, 2))};
            return true;
        };
        DrawCorrections(fOut, fullFamilies, pairings, "SignalExtracted", "signal-extracted", fetchSig);
    }

    // CheapSigExtract: the consumer's *VsMassRegion moments (no fit, no sideband model), three ways, each region corrected
    // with its OWN kappa -- the inversion is a statement about one population, so one region's response is never applied
    // to another's candidates:
    //   InPeakQA/, OffPeakQA/  one region alone, background left in. QA ONLY, never a result: P and kappa are mixes of
    //                          signal and background (InPeakQA against SidebandSubtracted shows how far the background moves it)
    //   SidebandSubtracted/    the signal estimate: every moment sideband-subtracted first, then corrected
    std::cout << " -> Corrections, CheapSigExtract (in peak, off peak, sideband-subtracted)...\n";

    // One region's moments as a sample (index 1 peak, 0 sideband), with kappa's exact <u v w> covariance, as in Raw
    auto makeRegionFetcher = [&](int region) -> CorrectionFetcher {
        return [&, region](const FakePairing& pr, const std::string& sfx, const std::string& proxy, CorrectionInputs& in) {
            RegionSample r;
            const std::string file = pr.consumerBase + "/ConsumerResults_" + sfx + ".root";
            if (!LoadRegionSample(cache, file, proxy, r)) return false;
            if (InputsOfSample(RegionAsSample(r, region), in)) return true;
            WarnSkip(proxy + (region == 1 ? " peak" : " sideband") + " sample", "fewer than two candidates", file);
            return false;
        };
    };
    DrawCorrections(fOut, fullFamilies, pairings, "CheapSigExtract/InPeakQA", "in peak only, QA: background not removed", makeRegionFetcher(1));
    DrawCorrections(fOut, fullFamilies, pairings, "CheapSigExtract/OffPeakQA", "sideband only, QA: background candidates", makeRegionFetcher(0));
    DrawDifferentialCorrections(cache, fOut, fullFamilies, pairings, "CheapSigExtract/InPeakQA", "in peak only, QA: background not removed",
                                [](const DifferentialSample& d, size_t i, CorrectionInputs& in) { return InputsOfSample(RegionAsSample(d.region[i], 1), in); });
    DrawDifferentialCorrections(cache, fOut, fullFamilies, pairings, "CheapSigExtract/OffPeakQA", "sideband only, QA: background candidates",
                                [](const DifferentialSample& d, size_t i, CorrectionInputs& in) { return InputsOfSample(RegionAsSample(d.region[i], 0), in); });

    // Sideband-subtracted: <y>_S = (Sum_P y - Sum_B y) / (N_P - N_B) for each of R, u v and w, then the usual correction
    // (see InputsOfSidebandSubtracted: kappa_S keeps its exact covariance)
    const CorrectionFetcher fetchCheap = [&](const FakePairing& pr, const std::string& sfx, const std::string& proxy,
                                             CorrectionInputs& in) {
        RegionSample r;
        const std::string file = pr.consumerBase + "/ConsumerResults_" + sfx + ".root";
        if (!LoadRegionSample(cache, file, proxy, r)) return false;
        if (InputsOfSidebandSubtracted(r, in)) return true;
        WarnSkip(proxy + " sideband-subtracted inputs", "no usable signal estimate (a subtracted count is not positive)", file);
        return false;
    };
    DrawCorrections(fOut, fullFamilies, pairings, "CheapSigExtract/SidebandSubtracted", "sideband-subtracted", fetchCheap);
    DrawDifferentialCorrections(cache, fOut, fullFamilies, pairings, "CheapSigExtract/SidebandSubtracted", "sideband-subtracted",
                                [](const DifferentialSample& d, size_t i, CorrectionInputs& in) { return InputsOfSidebandSubtracted(d.region[i], in); });
    DrawKappaByRegion(cache, fOut, fullFamilies, pairings);

    cache.CloseAll(); // Release every input file before closing the output

    fOut->Close();
    delete fOut;
    PrintSkipSummary(); // Last, so it closes the log
    std::cout << " Done!\n";
}

// To be able to run this both as a compiled C++ macro or a ROOT script:
#ifndef __CINT__
int main(int argc, char** argv) {
    if (argc < 2) { // Check is argc < 2 because this is the bare minimum. Everything after consumerDir is optional
        std::cerr << "Usage: " << argv[0]
                  << " <consumerDir> [mcRefDir] [ppRefDir] [toyModelPath]"
                     " [doIndividualComparisons] [sigExtractDir]\n"
                     "  sigExtractDir defaults to <consumerDir>/../results_SigExtract;"
                     " pass \"none\" to skip the signal-extraction plots.\n"
                     "  Longitudinal-ring outputs (_useRingZ, _useRingZT) in <consumerDir> are detected by name and"
                     " summarised as <family>_RingZ and <family>_RingZT.\n";
        return 1;
    }
    std::string consumerDir = argv[1];
    std::string mcRefDir = (argc > 2) ? argv[2] : "";
    std::string ppRefDir = (argc > 3) ? argv[3] : "";
    std::string toyModelPath = (argc > 4) ? argv[4] : "";
    // The old argv[5] (the kinematic-cut folder) is gone: the consumer books a single Ring/ folder.
    // signalExtractionRing writes next to the consumer output, so the default follows the pipeline
    // layout (<workDir>/results_consumer and <workDir>/results_SigExtract) and can be overridden.
    // Passing "none" disables the signal-extraction plots entirely.
    std::string sigExtractDir = (argc > 6) ? argv[6] : (consumerDir + "/../results_SigExtract");
    if (sigExtractDir == "none") sigExtractDir = "";
    // Accepts 1/0, true/false, yes/no. Anything else (including an absent argument) leaves it off.
    bool doIndividual = false;
    if (argc > 5) {
        std::string flag = argv[5];
        doIndividual = (flag == "1" || flag == "true" || flag == "yes" || flag == "on");
    }

    auxiliarySummaryPlots(consumerDir, mcRefDir, ppRefDir, toyModelPath, doIndividual,
                          sigExtractDir);
    return 0;
}
#endif