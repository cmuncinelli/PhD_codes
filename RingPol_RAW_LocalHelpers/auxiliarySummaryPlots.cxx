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

#include "TFile.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include "TSystem.h"
#include "TDirectory.h"

// ---------------------------------------------------------
// Configuration Data Structures
// ---------------------------------------------------------

// Holds styling and file suffix for a single systematic variation
struct VariationConfig {
    std::string suffix;      // e.g., "JustLambda_forceDatalikeJet"
    std::string legendLabel; // e.g., "Data-like Jet"
    int color;
    int lineStyle;
    int markerStyle;
    bool isData;             // Flag to indicate if this is the main, thick black line
};

// Represents a full family (Lambda, AntiLambda, Both) and its associated variations
struct FamilyConfig {
    std::string familyName;                 // e.g., "Lambda"
    std::string dataSuffix;                 // e.g., "JustLambda"
};

// Represents the different 1D TProfiles we want to extract from the files
struct ProfileConfig {
    std::string subDir;       // In-file directory, WITHOUT the task dir (e.g. "EtaDependence/")
    std::string profileName;  // e.g., "pRingObservableEtaLeadP"
    std::string xAxisTitle;   // e.g., "#eta_{LeadP}"
    std::string yAxisTitle;   // e.g., "R"
    bool underCutFolder;      // true  -> <task>/<cutFolder>/<subDir><profileName> (per kinematic-cut family)
                              // false -> <task>/<subDir><profileName>            (task-level, e.g. EtaStudy/)
};

// A set of profiles that share an output folder and the same set of valid transforms.
// Adding a new observable to the macro means adding an entry here and nothing else.
struct ObservableGroup {
    std::string groupName;                  // Output folder name, e.g. "EtaProxy"
    std::vector<ProfileConfig> profiles;
    bool doFold;                            // Fold about x = 0. Only meaningful for axes that are
                                            // symmetric around zero (eta, PVz). Must stay false for
                                            // pT, mass, DeltaPhi, ... where folding is nonsense.
};

// An external reference dataset drawn alongside the in-wagon systematics
// (MC baseline, pp baseline, Toy Model, ...). Collapsing these into one list removes
// the previously triplicated "if (pMC) ... if (pPP) ... if (pToy) ..." blocks, and with them
// the manual index-offset arithmetic that had to be kept in sync by hand.
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

// ---------------------------------------------------------
// Helper 0: Small ROOT navigation utilities
// ---------------------------------------------------------

// Returns an existing subdirectory or creates it. Used to keep directory creation lazy:
// nothing is created until we actually have a canvas to write into it, so consumer files
// that do not enable a given cut folder do not leave a trail of empty folders behind.
TDirectory* EnsureDir(TDirectory* parent, const std::string& name) {
    if (!parent) return nullptr;
    TDirectory* dir = parent->GetDirectory(name.c_str());
    if (!dir) dir = parent->mkdir(name.c_str());
    return dir;
}

// Builds the full in-file path of a profile, honouring the selected cut folder
std::string BuildProfilePath(const ProfileConfig& profConfig, const std::string& cutFolder) {
    if (profConfig.underCutFolder) {
        return TASK_DIR_IN_FILE + cutFolder + "/" + profConfig.subDir + profConfig.profileName;
    }
    return TASK_DIR_IN_FILE + profConfig.subDir + profConfig.profileName;
}

// ---------------------------------------------------------
// Helper 1: Cached file access
// ---------------------------------------------------------
// Every fetch used to open and close the ROOT file again, so the number of opens scaled as
// (families x observables x variations). With the observable list growing this became a real
// cost for the macro, and it is pure overhead: the same handful of files is reopened
// hundreds of times. Each file is now opened once, kept open, and closed at the very end.
// Failed opens are cached too, so a missing systematic warns once instead of once per observable.
class FileCache {
public:
    ~FileCache() { CloseAll(); }

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

    // Fetches an object and returns a detached clone that survives CloseAll()
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
// Only valid for axes symmetric about zero. Callers must gate this on ObservableGroup::doFold.
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
// The integrated summaries draw one point per variation on a labelled axis, which means one
// single-filled histogram per variation. This used to be written out inline once per reference
// type; it is identical every time apart from the bin index, so it lives here now.
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
TProfile* GetToyModelProfile(FileCache& cache, const std::string& filePath, const std::string& profileName) {
    if (filePath.empty()) return nullptr;

    // Remark: The Toy Model trigger/jet direction is always random and fundamentally for Lambdas (though it wouldn't
    // really make any meaningful distinction for antiLambdas or a case where both hyperons are considered).
    // We map it to LeadP, Jet, 2ndJet, etc., strictly for plotting/comparison purposes!
    //
    // TODO: this substring heuristic only covers the eta-dependence observables and will silently pick the
    // wrong profile for anything else (e.g. "pRingObservableEtaLambdaMinusEtaJet" also contains "Lambda").
    // Once the Toy Model's own histogram set is revisited, replace this with an explicit per-observable
    // toy path carried in ProfileConfig, where an empty path means "no toy equivalent, skip".
    std::string internalPath = "WithEtaGate/BothCuts/All/";
    if (profileName.find("Lambda") != std::string::npos) {
        internalPath += "pRingProxyJetVsEta";
    } else {
        internalPath += "pRingProxyJetVsEtaJet";
    }

    return cache.FetchClone<TProfile>(filePath, internalPath);
}

// ---------------------------------------------------------
// Helper 2: Draw Comparison Canvas
// ---------------------------------------------------------

// A temporary struct to bundle a fetched profile/histogram with its styling instructions
struct ProfileBundle {
    TH1* profile; // Converted to TH1* which is also the parent type of TProfile* in order to add the R(eta_pos) + R(eta_pos) plots into this structure
    VariationConfig config;
};

void DrawComparisonCanvas(const std::vector<ProfileBundle>& bundles,
                          const std::string& canvasName,
                          const std::string& canvasTitle,
                          TDirectory* outDir,
                          const ProfileConfig& profConfig,
                          const std::string& customXTitle = "",
                          const std::string& customYTitle = "") {
    
    if (bundles.empty() || !outDir) return;

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

    // 2. Setup Canvas and Legend
    TCanvas* c = new TCanvas(canvasName.c_str(), canvasTitle.c_str(), 800, 600);
    c->SetLeftMargin(0.12); // todo: could probably mess with these margins a little after the macro finishes running
    c->SetBottomMargin(0.12);
    c->SetGridx();
    c->SetGridy();

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
    TH1* frame = c->DrawFrame(xMin, yMin, xMax, yMax);
    frame->SetTitle(canvasTitle.c_str());
    frame->GetXaxis()->SetTitle(xTitle.c_str());
    frame->GetYaxis()->SetTitle(yTitle.c_str());
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);

    // 4. Style everything and populate the legend first
    ProfileBundle dataBundle;
    bool hasData = false;

    for (auto& bundle : bundles) {
        TH1* p = bundle.profile;
        
        p->SetLineColor(bundle.config.color);
        p->SetMarkerColor(bundle.config.color);
        p->SetLineStyle(bundle.config.lineStyle);
        p->SetMarkerStyle(bundle.config.markerStyle);

        if (bundle.config.isData) {
            p->SetLineWidth(3); // Thick line for the data
            dataBundle = bundle;
            hasData = true;
        } else {
            p->SetLineWidth(1); // Standard width for systematics
        }
        
        // Add to legend (Data will be added first assuming it's passed first in the vector)
        leg->AddEntry(p, bundle.config.legendLabel.c_str(), "pe");
    }

    // 5. Draw the profiles (Systematics first, then Data on top)
    for (auto& bundle : bundles) {
        if (!bundle.config.isData) {
            bundle.profile->Draw("PE SAME");
        }
    }
    
    if (hasData) { // Drawing the data after all others
        dataBundle.profile->Draw("PE SAME");
    }

    leg->Draw();

    // 6. Save to the specific directory within our output file
    outDir->cd();
    c->Write();
    
    // Clean up to prevent memory leaks during the loop (object should already have been written)
    delete leg;
    delete c; 
}

// ---------------------------------------------------------
// Helper 3: Draw Integrated Ring Observable Canvas
// ---------------------------------------------------------
void DrawIntegratedCanvas(const std::vector<ProfileBundle>& bundles,
                          const std::vector<std::string>& labels,
                          const std::string& canvasName,
                          const std::string& canvasTitle,
                          TDirectory* outDir,
                          const std::string& yTitle,
                          bool isSubtracted,
                          bool drawLegend = false) {
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
    TCanvas* c = new TCanvas(canvasName.c_str(), canvasTitle.c_str(), std::max(800, nBins * 100), 600);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.15); // Larger bottom margin to fit text labels
    c->SetGridx();

    // 2. Create a dummy histogram just to draw the custom axes and labels
    TH1D* frame = new TH1D((canvasName + "_frame").c_str(), canvasTitle.c_str(), nBins, 0, nBins);
    frame->SetDirectory(nullptr);
    for (int i = 0; i < nBins; ++i) {
        frame->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
    }
    frame->GetYaxis()->SetTitle(yTitle.c_str());
    frame->GetYaxis()->SetRangeUser(yMin, yMax);
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
    
    outDir->cd();
    c->Write();
    delete frame;
    delete c;
}

// ---------------------------------------------------------
// Main Macro 
// ---------------------------------------------------------
void auxiliarySummaryPlots(const std::string& consumerDir,
                           const std::string& mcRefDir = "",
                           const std::string& ppRefDir = "",
                           const std::string& toyModelPath = "",
                           const std::string& cutFolder = DEFAULT_CUT_FOLDER) {
    
    // 1. Define Systematic Variations (the list of all useful variations I would like to track into this plot)
    // The data config (empty suffix) is handled separately in the logic to ensure it is always first
    std::vector<VariationConfig> sysVariations = {
        {"_forceRandJet",                 "Rand Jet",                  kBlue,     1, 20, false}, // Full circle
        {"_forceDatalikeJet",             "Data-like Jet",             kRed,      1, 21, false}, // Square
        {"_forceDatalikeJet_10resamples", "Data-like Jet (10 resam.)", kRed,      2, 25, false}, // Open square
        {"_forcePerpToJet",               "Perp to Jet",               kGreen+2,  1, 22, false}, // Triangle up
        {"_forcePerpToJet_10resamples",   "Perp to Jet (10 resam.)",   kGreen+2,  2, 26, false}, // Open triangle up
        {"_forcePreviousJet",             "Prev Jet",                  kOrange+1, 1, 34, false}, // Full cross (bold "+")
        {"_MixedEventProxies",            "MixedEv Jet",               kGray+1,   1, 24, false}  // Open circle
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
                {"EtaDependence/", "pRingObservableEtaLeadP",  "#eta_{LeadP}",  "R", true},
                {"EtaDependence/", "pRingObservableEtaJet",    "#eta_{Jet}",    "R", true},
                {"EtaDependence/", "pRingObservableEta2ndJet", "#eta_{2ndJet}", "R", true}
            },
            true // Eta axis is symmetric about zero, so folding is meaningful
        },
        {
            "EtaV0", // x-axis is the eta of the hyperon itself, split by which proxy defined the ring
            {
                {"EtaDependence/", "pRingObservableEtaLambda",       "#eta_{#Lambda}",         "R", true},
                {"EtaDependence/", "pRingObservableEtaLambdaLeadP",  "#eta_{#Lambda(LeadP)}",  "R", true},
                {"EtaDependence/", "pRingObservableEtaLambda2ndJet", "#eta_{#Lambda(2ndJet)}", "R", true}
            },
            true
        }
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

    // 5. Main Processing Loops
    for (const auto& fam : families) {
        std::cout << " -> Processing Family: " << fam.familyName << std::endl;

        // Directories are created lazily (only once something is actually written into them),
        // so a consumer file that does not enable the requested cut folder leaves no empty folders.
        TDirectory* famDir = nullptr;

        // Variables to accumulate the integrated observable summary across all variations for this family
        // (a "grand summary" for short)
        std::vector<std::pair<double, double>> grandDataVals;
        std::vector<std::vector<std::pair<double, double>>> grandExtVals(externals.size());  // One row per external reference
        std::vector<std::vector<std::pair<double, double>>> grandSysVals(sysVariations.size()); // Track all systematic variations
        std::vector<std::string> grandLabels;

        for (const auto& group : observableGroups) {

            TDirectory* groupDir = nullptr; // Also created lazily

            for (const auto& profConfig : group.profiles) {

                const std::string profilePath = BuildProfilePath(profConfig, cutFolder);

                // Fetch Data:
                std::string dataFile = consumerDir + "/ConsumerResults_" + fam.dataSuffix + ".root";
                TProfile* pData = cache.FetchClone<TProfile>(dataFile, profilePath);
                if (!pData) continue; // Skip to next observable if data is missing (nothing is created on disk)

                // Only now that we know there is something to draw do we materialize the directories
                if (!famDir)   famDir   = EnsureDir(fOut, fam.familyName);
                if (!groupDir) groupDir = EnsureDir(famDir, group.groupName);
                TDirectory* obsDir = EnsureDir(groupDir, profConfig.profileName);

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
                        pExt = GetToyModelProfile(cache, ext.basePath, profConfig.profileName);
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

                // ---------------------------------------------------------
                // DRAWING THE VARIATIONS
                // ---------------------------------------------------------
                
                // --- 1. Standard (Unfolded) Plots ---
                    // Variation 1: Data Only
                DrawComparisonCanvas({{pData, dataConfig}}, "Canvas_DataOnly", fam.familyName + " Data", overlayDir, profConfig);
                    // Variation 2: Systematics (+ Data)
                DrawComparisonCanvas(allSystematics, "Canvas_Systematics", fam.familyName + " Systematics", overlayDir, profConfig);
                    // Variation 3+: one canvas per external reference, then the all-in-one
                std::vector<ProfileBundle> allInOne = allSystematics;
                for (const auto& ext : presentExternals) {
                    const std::string& label = ext.second.config.legendLabel;
                    DrawComparisonCanvas({ext.second}, "Canvas_Only_" + externals[ext.first].config.legendLabel,
                                         fam.familyName + " " + label, overlayDir, profConfig);
                    allInOne.push_back(ext.second);
                }
                if (!presentExternals.empty()) {
                    DrawComparisonCanvas(allInOne, "Canvas_AllInOne", fam.familyName + " All Comparisons", overlayDir, profConfig);
                }

                // --- 2. Folded Plots ---
                if (group.doFold) {
                    DrawComparisonCanvas({allFoldedSystematics[0]}, "Canvas_Folded_DataOnly", fam.familyName + " (R(#eta_{pos})+R(#eta_{neg})) Data", foldedDir, profConfig, foldedXTitle);
                    DrawComparisonCanvas(allFoldedSystematics, "Canvas_Folded_Systematics", fam.familyName + " (R(#eta_{pos})+R(#eta_{neg})) Systematics", foldedDir, profConfig, foldedXTitle);

                    std::vector<ProfileBundle> allInOneFolded = allFoldedSystematics; // For folded plots tracking
                    for (size_t i = 0; i < allFoldedExternals.size(); ++i) {
                        const std::string& label = allFoldedExternals[i].config.legendLabel;
                        DrawComparisonCanvas({allFoldedExternals[i]}, "Canvas_Folded_Only_" + label,
                                             fam.familyName + " (R(#eta_{pos})+R(#eta_{neg})) " + label, foldedDir, profConfig, foldedXTitle);
                        allInOneFolded.push_back(allFoldedExternals[i]);
                    }
                    if (!allFoldedExternals.empty()) {
                        DrawComparisonCanvas(allInOneFolded, "Canvas_Folded_AllInOne", fam.familyName + " (R(#eta_{pos})+R(#eta_{neg})) Comparisons", foldedDir, profConfig, foldedXTitle);
                    }
                }

                // --- 3. Subtracted Data - Systematics ---
                    // Data minus Data Systematics
                DrawComparisonCanvas(allSubtractedSystematics, "Canvas_Subtracted_Systematics", fam.familyName + " Systematics Difference", subDir, profConfig, "", subYTitle);
                    // Data minus each external reference
                std::vector<ProfileBundle> allInOneSubtracted = allSubtractedSystematics;
                for (size_t i = 0; i < allSubtractedExternals.size(); ++i) {
                    const std::string& label = allSubtractedExternals[i].config.legendLabel;
                    DrawComparisonCanvas({{pDataZero, dataConfig}, allSubtractedExternals[i]}, "Canvas_Subtracted_Only_" + label,
                                         fam.familyName + " " + label + " Difference", subDir, profConfig, "", subYTitle);
                    allInOneSubtracted.push_back(allSubtractedExternals[i]);
                }
                    // All-In-One
                if (!allSubtractedExternals.empty()) {
                    DrawComparisonCanvas(allInOneSubtracted, "Canvas_Subtracted_AllInOne", fam.familyName + " All Comparisons Difference", subDir, profConfig, "", subYTitle);
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
                    DrawComparisonCanvas(allSubtractedMass, "Canvas_Subtracted_MassSelection", fam.familyName + " Mass Selection Difference", massDir, profConfig, "", subYTitle);
                }

                // --- 4. Individual Comparisons (One-by-One) ---
                // (condensed all of this into a single block of code because the loops become cleaner!)
                    // Create a sub-directory specifically for the individual comparisons to keep things organized
                TDirectory* indivDir = EnsureDir(obsDir, "Individual_Comparisons");

                // Loop through all variations (starting at i = 1 to skip the data comparing against itself)
                for (size_t i = 1; i < allSystematics.size(); ++i) {
                    std::string sysLabel = allSystematics[i].config.legendLabel;
                    std::string sysSuffix = allSystematics[i].config.suffix; // e.g., "_forceDatalikeJet"

                    // Standard
                    DrawComparisonCanvas({allSystematics[0], allSystematics[i]}, "Canvas_Standard" + sysSuffix, fam.familyName + " Data vs " + sysLabel, indivDir, profConfig);
                    // Folded
                    if (group.doFold) {
                        DrawComparisonCanvas({allFoldedSystematics[0], allFoldedSystematics[i]}, "Canvas_Folded" + sysSuffix, fam.familyName + " Folded vs " + sysLabel, indivDir, profConfig, foldedXTitle);
                    }
                    // Subtracted
                    DrawComparisonCanvas({allSubtractedSystematics[0], allSubtractedSystematics[i]}, "Canvas_Subtracted" + sysSuffix, fam.familyName + " Difference vs " + sysLabel, indivDir, profConfig, "", subYTitle);
                }

                // Do the same for every external reference, which are introduced separately
                for (size_t i = 0; i < presentExternals.size(); ++i) {
                    const std::string& label = presentExternals[i].second.config.legendLabel;

                    DrawComparisonCanvas({allSystematics[0], presentExternals[i].second}, "Canvas_Standard_" + label, fam.familyName + " Data vs " + label, indivDir, profConfig);
                    if (group.doFold) {
                        DrawComparisonCanvas({allFoldedSystematics[0], allFoldedExternals[i]}, "Canvas_Folded_" + label, fam.familyName + " Folded vs " + label, indivDir, profConfig, foldedXTitle);
                    }
                    DrawComparisonCanvas({allSubtractedSystematics[0], allSubtractedExternals[i]}, "Canvas_Subtracted_" + label, fam.familyName + " Difference vs " + label, indivDir, profConfig, "", subYTitle);
                }

                // --- 5. Integrated Summaries ---
                TDirectory* integDir = EnsureDir(obsDir, "Integrated_Summary");

                std::vector<std::string> integLabels;
                std::vector<ProfileBundle> integBundles;
                std::vector<ProfileBundle> integSubBundles;
                std::vector<TH1*> integToDelete;

                int numCats = allSystematics.size() + presentExternals.size();
                
                    // Compute Integrated Data and save to Grand Summary
                auto dataInteg = GetIntegratedProfile(pData);
                grandDataVals.push_back(dataInteg);
                grandLabels.push_back(profConfig.xAxisTitle);

                // Populate the grand summary for systematics
                // (allSystematics[0] is the data, so sysVariations[i] corresponds to allSystematics[i+1])
                for (size_t i = 0; i < sysVariations.size(); ++i) {
                    if (i + 1 < allSystematics.size()) {
                        auto sysInteg = GetIntegratedProfile(dynamic_cast<TProfile*>(allSystematics[i+1].profile));
                        grandSysVals[i].push_back(sysInteg);
                    }
                }
                // ... and for the external references, indexed by their position in the externals list
                for (const auto& ext : presentExternals) {
                    grandExtVals[ext.first].push_back(GetIntegratedProfile(dynamic_cast<TProfile*>(ext.second.profile)));
                }

                    // Loop over variations and construct isolated categorical bins
                for (size_t i = 0; i < allSystematics.size(); ++i) {
                    integLabels.push_back(allSystematics[i].config.legendLabel);
                    
                    auto integ = GetIntegratedProfile(dynamic_cast<TProfile*>(allSystematics[i].profile));
                    
                    TH1D* hInt = MakeCategoricalPoint(Form("Integ_%zu", i), i, numCats, integ.first, integ.second);
                    
                    // Performing also the subtraction between integrated observables, variation by variation
                    // (error bars are also properly propagated)
                    TH1D* hIntSub = (i == 0)
                        ? MakeCategoricalPoint(Form("IntegSub_%zu", i), i, numCats, 0.0, 0.0)
                        : MakeCategoricalPoint(Form("IntegSub_%zu", i), i, numCats,
                                               dataInteg.first - integ.first,
                                               std::sqrt(dataInteg.second*dataInteg.second + integ.second*integ.second));

                    integBundles.push_back({hInt, allSystematics[i].config});
                    integSubBundles.push_back({hIntSub, allSystematics[i].config});
                    integToDelete.push_back(hInt);
                    integToDelete.push_back(hIntSub);
                }
                
                    // Add the external references to the integrated list, continuing the same bin numbering
                for (size_t i = 0; i < presentExternals.size(); ++i) {
                    const auto& extBundle = presentExternals[i].second;
                    integLabels.push_back(extBundle.config.legendLabel);

                    auto integExt = GetIntegratedProfile(dynamic_cast<TProfile*>(extBundle.profile));
                    int idx = static_cast<int>(allSystematics.size() + i); // No manual offset bookkeeping any more

                    TH1D* hIntExt = MakeCategoricalPoint(Form("Integ_Ext%zu", i), idx, numCats, integExt.first, integExt.second);
                    TH1D* hIntSubExt = MakeCategoricalPoint(Form("IntegSub_Ext%zu", i), idx, numCats,
                                                            dataInteg.first - integExt.first,
                                                            std::sqrt(dataInteg.second*dataInteg.second + integExt.second*integExt.second));

                    integBundles.push_back({hIntExt, extBundle.config});
                    integSubBundles.push_back({hIntSubExt, extBundle.config});
                    integToDelete.push_back(hIntExt);
                    integToDelete.push_back(hIntSubExt);
                }

                // 4. Draw them using the custom categorical plotter
                DrawIntegratedCanvas(integBundles, integLabels, "Canvas_Integrated", fam.familyName + " Integrated Summary", integDir, "Integrated " + profConfig.yAxisTitle, false);
                DrawIntegratedCanvas(integSubBundles, integLabels, "Canvas_Integrated_Subtracted", fam.familyName + " Integrated Differences", integDir, "#Delta" + profConfig.yAxisTitle + " (Data - Var)", true);

                // 5. Draw Mass Selection Integrated Summaries
                if (allMassSystematics.size() > 1 && massDir) {
                    std::vector<std::string> massIntegLabels;
                    std::vector<ProfileBundle> massIntegBundles;
                    std::vector<ProfileBundle> massIntegSubBundles;
                    int numMassCats = allMassSystematics.size();
                    
                    for (size_t i = 0; i < allMassSystematics.size(); ++i) {
                        massIntegLabels.push_back(allMassSystematics[i].config.legendLabel);
                        auto integ = GetIntegratedProfile(dynamic_cast<TProfile*>(allMassSystematics[i].profile));
                        
                        TH1D* hInt = MakeCategoricalPoint(Form("MassInteg_%zu", i), i, numMassCats, integ.first, integ.second);
                        
                        TH1D* hIntSub = (i == 0)
                            ? MakeCategoricalPoint(Form("MassIntegSub_%zu", i), i, numMassCats, 0.0, 0.0)
                            : MakeCategoricalPoint(Form("MassIntegSub_%zu", i), i, numMassCats,
                                                   dataInteg.first - integ.first,
                                                   std::sqrt(dataInteg.second*dataInteg.second + integ.second*integ.second));
                        
                        massIntegBundles.push_back({hInt, allMassSystematics[i].config});
                        massIntegSubBundles.push_back({hIntSub, allMassSystematics[i].config});
                        
                        integToDelete.push_back(hInt);     // Now safely in scope!
                        integToDelete.push_back(hIntSub);
                    }
                    
                    DrawIntegratedCanvas(massIntegBundles, massIntegLabels, "Canvas_Integrated_MassSelection", fam.familyName + " Mass Selection Integrated", massDir, "Integrated " + profConfig.yAxisTitle, false);
                    DrawIntegratedCanvas(massIntegSubBundles, massIntegLabels, "Canvas_Integrated_Subtracted_MassSelection", fam.familyName + " Mass Selection Integrated Differences", massDir, "#Delta" + profConfig.yAxisTitle + " (Data - Var)", true);
                }

                // Cleanup dynamically allocated profiles for this observable iteration
                for (auto p : integToDelete) delete p;
                for (auto p : profilesToDelete) delete p;
                for (auto p : foldedToDelete) delete p;
                for (auto p : subtractedToDelete) delete p;
            } // end of profiles loop, inside groups loop
        } // end of groups loop, inside families loop

        // --- 6. Integrated observable Summary (Observable in each variation vs Observable for Data) ---
        // Comparing eta_Jet and eta_Lambda integrations is a way of probing if there were any overflows/underflows
        // and the such (there were! I then modified the eta axis range in the consumer, so this was a good catch!
        // The single-bin integrated observables from other scripts were impervious to this problem though, so good news!),
        // or even other possible problems lingering in the code (none found after the aformentioned fix!)
        int nGrand = grandLabels.size();
        if (nGrand > 0 && famDir) {
            std::vector<ProfileBundle> grandBundles;
            std::vector<TH1*> grandToDelete; // For memory cleanup
            
            // Add Data
            TH1D* hGrandData = new TH1D("GrandData", "", nGrand, 0, nGrand);
            hGrandData->SetDirectory(nullptr);
            for(int i = 0; i < nGrand; ++i) {
                hGrandData->SetBinContent(i + 1, grandDataVals[i].first);
                hGrandData->SetBinError(i + 1, grandDataVals[i].second);
            }
            grandBundles.push_back({hGrandData, dataConfig});
            grandToDelete.push_back(hGrandData);

            // Add All Systematics
            for (size_t sysIdx = 0; sysIdx < sysVariations.size(); ++sysIdx) {
                if (grandSysVals[sysIdx].size() == (size_t)nGrand) {
                    TH1D* hGrandSys = new TH1D(Form("GrandSys_%zu", sysIdx), "", nGrand, 0, nGrand);
                    hGrandSys->SetDirectory(nullptr);
                    for (int i = 0; i < nGrand; ++i) {
                        hGrandSys->SetBinContent(i + 1, grandSysVals[sysIdx][i].first);
                        hGrandSys->SetBinError(i + 1, grandSysVals[sysIdx][i].second);
                    }
                    grandBundles.push_back({hGrandSys, sysVariations[sysIdx]});
                    grandToDelete.push_back(hGrandSys);
                }
            }

            // Add every external reference that was present for all observables
            for (size_t extIdx = 0; extIdx < externals.size(); ++extIdx) {
                if (externals[extIdx].basePath.empty()) continue;
                if (grandExtVals[extIdx].size() != (size_t)nGrand) continue;

                TH1D* hGrandExt = new TH1D(Form("GrandExt_%zu", extIdx), "", nGrand, 0, nGrand);
                hGrandExt->SetDirectory(nullptr);
                for (int i = 0; i < nGrand; ++i) {
                    hGrandExt->SetBinContent(i + 1, grandExtVals[extIdx][i].first);
                    hGrandExt->SetBinError(i + 1, grandExtVals[extIdx][i].second);
                }
                grandBundles.push_back({hGrandExt, externals[extIdx].config});
                grandToDelete.push_back(hGrandExt);
            }
            
            // Draw the Brute Force Canvas
            DrawIntegratedCanvas(grandBundles, grandLabels, "Canvas_BruteForce_Summary", fam.familyName + " Global Summary Across Kinematics", famDir, "Integrated R", false, true);
            
            for (auto p : grandToDelete) delete p;
        }
    } // end of families loop

    cache.CloseAll(); // Release every input file before closing the output

    fOut->Close();
    delete fOut;
    std::cout << " Done!\n";
}

// To be able to run this both as a compiled C++ macro or a ROOT script:
#ifndef __CINT__
int main(int argc, char** argv) {
    if (argc < 2) { // Check is argc < 2 because this is the bare minimum. mcRefDir, ppRefDir, toyModelPath and cutFolder are all optionals
        std::cerr << "Usage: " << argv[0] << " <consumerDir> [mcRefDir] [ppRefDir] [toyModelPath] [cutFolder]\n";
        return 1;
    }
    std::string consumerDir = argv[1];
    std::string mcRefDir = (argc > 2) ? argv[2] : "";
    std::string ppRefDir = (argc > 3) ? argv[3] : "";
    std::string toyModelPath = (argc > 4) ? argv[4] : "";
    std::string cutFolder = (argc > 5) ? argv[5] : DEFAULT_CUT_FOLDER;
    
    auxiliarySummaryPlots(consumerDir, mcRefDir, ppRefDir, toyModelPath, cutFolder);
    return 0;
}
#endif