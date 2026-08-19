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
    std::string profileName;  // e.g., "pRingObservableEtaLeadP"
    std::string xAxisTitle;   // e.g., "#eta_{LeadP}"
    std::string yAxisTitle;   // e.g., "R"
};

// ---------------------------------------------------------
// Global Definitions
// ---------------------------------------------------------
const std::string BASE_DIR_IN_FILE = "lambdajetpolarizationionsderived/Ring/EtaDependence/";

// ---------------------------------------------------------
// Helper 1: Safe Profile Extraction
// ---------------------------------------------------------
TProfile* GetProfile(const std::string& filePath, const std::string& profileName) {
    TFile* fIn = TFile::Open(filePath.c_str(), "READ");
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "  [Warning] Could not open file: " << filePath << std::endl;
        if (fIn) fIn->Close();
        return nullptr;
    }

    std::string fullPath = BASE_DIR_IN_FILE + profileName;
    TProfile* prof = (TProfile*)fIn->Get(fullPath.c_str());
    
    if (!prof) {
        std::cerr << "  [Warning] Profile not found: " << fullPath 
                  << " in " << filePath << std::endl;
        fIn->Close();
        return nullptr;
    }

    // Clone it so it persists in memory after we close the input file
    TProfile* pClone = (TProfile*)prof->Clone();
    pClone->SetDirectory(nullptr); 
    
    fIn->Close();
    return pClone;
}

// ---------------------------------------------------------
// Helper 1.1: Fold Profile (Sum Positive and Negative Eta)
// ---------------------------------------------------------
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
    
    hFolded->SetDirectory(nullptr); // Protect memory
    return hFolded;
}

// ---------------------------------------------------------
// Helper 1.2: Subtract Profiles (Data - Systematics)
// ---------------------------------------------------------
// (Also handles binning mismatches appropriately)
TH1D* SubtractProfiles(TH1* pData, TH1* pSys, const std::string& newName) {
    if (!pData || !pSys) return nullptr;
    
    // Safely clone the data binning structure
    TH1D* hSub = nullptr;
    if (auto pDataProf = dynamic_cast<TProfile*>(pData)) {
        hSub = pDataProf->ProjectionX(newName.c_str());
    } else {
        hSub = (TH1D*)pData->Clone(newName.c_str());
    }
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
    
    hSub->SetDirectory(nullptr);
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
    pClone->Rebin(pClone->GetNbinsX()); 
    
    double val = pClone->GetBinContent(1);
    double err = pClone->GetBinError(1);
    
    delete pClone;
    return {val, err};
}

// ---------------------------------------------------------
// Helper 1.4: Toy Model Profile Extraction
// ---------------------------------------------------------
TProfile* GetToyModelProfile(const std::string& filePath, const std::string& profileName) {
    if (filePath.empty()) return nullptr;
    
    TFile* fIn = TFile::Open(filePath.c_str(), "READ"); // Open as read-only, of course
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "  [Warning] Could not open Toy Model file: " << filePath << std::endl;
        if (fIn) fIn->Close();
        return nullptr;
    }

    // Remark: The Toy Model trigger/jet direction is always random and fundamentally for Lambdas (though it wouldn't
    // really make any meaningful distinction for antiLambdas or a case where both hyperons are considered).
    // We map it to LeadP, Jet, 2ndJet, etc., strictly for plotting/comparison purposes!
    std::string internalPath = "WithEtaGate/BothCuts/All/";
    if (profileName.find("Lambda") != std::string::npos) {
        internalPath += "pRingProxyJetVsEta";
    } else {
        internalPath += "pRingProxyJetVsEtaJet";
    }

    TProfile* prof = (TProfile*)fIn->Get(internalPath.c_str());
    if (!prof) {
        std::cerr << "  [Warning] Toy Model Profile not found: " << internalPath << std::endl;
        fIn->Close();
        return nullptr;
    }

    TProfile* pClone = (TProfile*)prof->Clone();
    pClone->SetDirectory(nullptr); 
    fIn->Close();
    return pClone;
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
    if (bundles.empty()) return;

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
void auxiliarySummaryPlots(const std::string& consumerDir, const std::string& mcRefDir = "", const std::string& ppRefDir = "", const std::string& toyModelPath = "") {
    
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

    // 3. Define the Profiles to Extract
    std::vector<ProfileConfig> profiles = {
        {"pRingObservableEtaLeadP",        "#eta_{LeadP}",           "R"},
        {"pRingObservableEtaJet",          "#eta_{Jet}",             "R"},
        {"pRingObservableEta2ndJet",       "#eta_{2ndJet}",          "R"},
        {"pRingObservableEtaLambda",       "#eta_{#Lambda}",         "R"},
        {"pRingObservableEtaLambda2ndJet", "#eta_{#Lambda(2ndJet)}", "R"},
        {"pRingObservableEtaLambdaLeadP",  "#eta_{#Lambda(LeadP)}",  "R"}
    };

    // 4. Create the output file
    std::string outFilePath = consumerDir + "/auxiliarySummaryPlots.root";
    TFile* fOut = new TFile(outFilePath.c_str(), "RECREATE");
    if (!fOut || fOut->IsZombie()) {
        std::cerr << "[Error] Could not create output file: " << outFilePath << std::endl;
        return;
    }

    std::cout << "========================================================\n";
    std::cout << " Generating Auxiliary Plots\n";
    std::cout << " Output: " << outFilePath << "\n";
    std::cout << "========================================================\n";

    // 5. Main Processing Loops
    for (const auto& fam : families) {
        std::cout << " -> Processing Family: " << fam.familyName << std::endl;
        
        // Create a directory for this family in the output ROOT file
        TDirectory* famDir = fOut->mkdir(fam.familyName.c_str());

        // Variables to accumulate the integrated observable summary across all variations for this family
        // (a "grand summary" for short)
        std::vector<std::pair<double, double>> grandDataVals;
        std::vector<std::pair<double, double>> grandMCVals;
        std::vector<std::pair<double, double>> grandPPVals;
        std::vector<std::pair<double, double>> grandToyVals;
        std::vector<std::vector<std::pair<double, double>>> grandSysVals(sysVariations.size()); // Track all systematic variations
        std::vector<std::string> grandLabels;

        for (const auto& profConfig : profiles) {
            
            // Create a subdirectory for the observable to keep things extremely tidy
            TDirectory* obsDir = famDir->mkdir(profConfig.profileName.c_str());

            // Fetch Data:
            std::string dataFile = consumerDir + "/ConsumerResults_" + fam.dataSuffix + ".root";
            TProfile* pData = GetProfile(dataFile, profConfig.profileName);
            if (!pData) continue; // Skip to next observable if data is missing

            // Fetch Systematics Data
            std::vector<ProfileBundle> allSystematics;
            allSystematics.push_back({pData, dataConfig}); // Always keep data at index 0 (it is the very first thing in this array)
            
            std::vector<TProfile*> profilesToDelete; // Track for memory cleanup
            profilesToDelete.push_back(pData);

            for (const auto& sys : sysVariations) {
                std::string sysFile = consumerDir + "/ConsumerResults_" + fam.dataSuffix + sys.suffix + ".root";
                TProfile* pSys = GetProfile(sysFile, profConfig.profileName);
                if (pSys) {
                    allSystematics.push_back({pSys, sys});
                    profilesToDelete.push_back(pSys);
                }
            }

            // Fetch MC Base Data (if requested via non-null path)
            TProfile* pMC = nullptr;
            if (!mcRefDir.empty()) {
                std::string mcFile = mcRefDir + "/ConsumerResults_" + fam.dataSuffix + ".root";
                pMC = GetProfile(mcFile, profConfig.profileName);
                if (pMC) profilesToDelete.push_back(pMC);
            }

            // Fetch PP Data
            TProfile* pPP = nullptr;
            if (!ppRefDir.empty()) {
                std::string ppFile = ppRefDir + "/ConsumerResults_" + fam.dataSuffix + ".root";
                pPP = GetProfile(ppFile, profConfig.profileName);
                if (pPP) profilesToDelete.push_back(pPP);
            }

            // Fetch Toy Model Data (if requested via non-null path)
            TProfile* pToy = nullptr;
            if (!toyModelPath.empty()) {
                pToy = GetToyModelProfile(toyModelPath, profConfig.profileName);
                if (pToy) profilesToDelete.push_back(pToy);
            }

            // Fetch Mass Selection Data
            std::vector<ProfileBundle> allMassSystematics;
            allMassSystematics.push_back({pData, dataConfig}); // Keep data as baseline
            
            for (const auto& mass : massVariations) {
                std::string massFile = consumerDir + "/ConsumerResults_" + fam.dataSuffix + mass.suffix + ".root";
                TProfile* pMass = GetProfile(massFile, profConfig.profileName);
                if (pMass) {
                    allMassSystematics.push_back({pMass, mass});
                    profilesToDelete.push_back(pMass); // Track for standard cleanup
                }
            }

            // ---------------------------------------------------------
            // FOLDING THE DATA (Summing Pos and Neg Eta)
            // ---------------------------------------------------------
            std::vector<ProfileBundle> allFoldedSystematics;
            std::vector<TH1*> foldedToDelete;
            
            for (const auto& bundle : allSystematics) {
                // Safely cast to TProfile since we know the original inputs were TProfiles
                TProfile* pOrig = dynamic_cast<TProfile*>(bundle.profile);
                std::string foldName = std::string(pOrig->GetName()) + "_Folded_" + bundle.config.suffix;
                
                TH1D* pFolded = FoldProfile(pOrig, foldName);
                allFoldedSystematics.push_back({pFolded, bundle.config});
                foldedToDelete.push_back(pFolded);
            }

            TH1D* pMCFolded = nullptr;
            if (pMC) {
                pMCFolded = FoldProfile(pMC, std::string(pMC->GetName()) + "_Folded_MC");
                foldedToDelete.push_back(pMCFolded);
            }

            TH1D* pPPFolded = nullptr;
            if (pPP) {
                pPPFolded = FoldProfile(pPP, std::string(pPP->GetName()) + "_Folded_PP");
                foldedToDelete.push_back(pPPFolded);
            }

            TH1D* pToyFolded = nullptr;
            if (pToy) {
                pToyFolded = FoldProfile(pToy, std::string(pToy->GetName()) + "_Folded_Toy");
                foldedToDelete.push_back(pToyFolded);
            }

            // Fold Mass Selection Data
            std::vector<ProfileBundle> allFoldedMass;
            for (const auto& bundle : allMassSystematics) {
                TProfile* pOrig = dynamic_cast<TProfile*>(bundle.profile);
                std::string foldName = std::string(pOrig->GetName()) + "_Folded_" + bundle.config.suffix;
                TH1D* pFolded = FoldProfile(pOrig, foldName);
                allFoldedMass.push_back({pFolded, bundle.config});
                foldedToDelete.push_back(pFolded);
            }

            // Creating the Modulus X-Axis title (e.g. |#eta_{LeadP}|)
            std::string foldedXTitle = "|" + profConfig.xAxisTitle + "|";

            // ---------------------------------------------------------
            // SUBTRACTING THE DATA (Data - Variations)
            // ---------------------------------------------------------
            std::vector<ProfileBundle> allSubtractedSystematics;
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

                // Subtract MC
            TH1D* pMCSubtracted = nullptr;
            if (pMC) {
                pMCSubtracted = SubtractProfiles(pData, pMC, std::string(pMC->GetName()) + "_Subtracted_MC");
                subtractedToDelete.push_back(pMCSubtracted);
            }

                // Same for pp data:
            TH1D* pPPSubtracted = nullptr;
            if (pPP) {
                pPPSubtracted = SubtractProfiles(pData, pPP, std::string(pPP->GetName()) + "_Subtracted_PP");
                subtractedToDelete.push_back(pPPSubtracted);
            }

                // Subtract for Toy Model
            TH1D* pToySubtracted = nullptr;
            if (pToy) {
                pToySubtracted = SubtractProfiles(pData, pToy, std::string(pToy->GetName()) + "_Subtracted_Toy");
                subtractedToDelete.push_back(pToySubtracted);
            }

                // Subtract Mass Selection Data
            std::vector<ProfileBundle> allSubtractedMass;
            allSubtractedMass.push_back({pDataZero, dataConfig}); // Re-use the existing zero-line reference
            
            for (size_t i = 1; i < allMassSystematics.size(); ++i) {
                const auto& bundle = allMassSystematics[i];
                std::string subName = std::string(pData->GetName()) + "_Subtracted_" + bundle.config.suffix;
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
            DrawComparisonCanvas({{pData, dataConfig}}, "Canvas_DataOnly", fam.familyName + " Data", obsDir, profConfig);
                // Variation 2: Systematics (+ Data)
            DrawComparisonCanvas(allSystematics, "Canvas_Systematics", fam.familyName + " Systematics", obsDir, profConfig);
                // Variation 3: MC Only (if available)
            std::vector<ProfileBundle> allInOne = allSystematics;
            if (pMC) {
                DrawComparisonCanvas({{pMC, mcConfig}}, "Canvas_MCOnly", fam.familyName + " MC", obsDir, profConfig);
                allInOne.push_back({pMC, mcConfig});
            }
                // Variation 4: pp data data (if available)
            if (pPP) {
                DrawComparisonCanvas({{pPP, ppConfig}}, "Canvas_PPOnly", fam.familyName + " PP Data", obsDir, profConfig);
                allInOne.push_back({pPP, ppConfig});
            }
                // Variation 5: Toy Model Only (if available)
            if (pToy) {
                DrawComparisonCanvas({{pToy, toyConfig}}, "Canvas_ToyOnly", fam.familyName + " Toy Model", obsDir, profConfig);
                allInOne.push_back({pToy, toyConfig});
            }
                // Variation 6: All-in-One (if MC and Toy Model available)
            if (pMC || pPP || pToy) {
                DrawComparisonCanvas(allInOne, "Canvas_AllInOne", fam.familyName + " All Comparisons", obsDir, profConfig);
            }

            // --- 2. Folded Plots ---
            DrawComparisonCanvas({allFoldedSystematics[0]}, "Canvas_Folded_DataOnly", fam.familyName + " (R(#eta_{pos})+R(#eta_{neg})) Data", obsDir, profConfig, foldedXTitle);
            DrawComparisonCanvas(allFoldedSystematics, "Canvas_Folded_Systematics", fam.familyName + " (R(#eta_{pos})+R(#eta_{neg})) Systematics", obsDir, profConfig, foldedXTitle);
            std::vector<ProfileBundle> allInOneFolded = allFoldedSystematics; // For folded plots tracking
            if (pMCFolded) {
                DrawComparisonCanvas({{pMCFolded, mcConfig}}, "Canvas_Folded_MCOnly", fam.familyName + " (R(#eta_{pos})+R(#eta_{neg})) MC", obsDir, profConfig, foldedXTitle);
                allInOneFolded.push_back({pMCFolded, mcConfig});
            }
            if (pPPFolded) {
                DrawComparisonCanvas({{pPPFolded, ppConfig}}, "Canvas_Folded_PPOnly", fam.familyName + " (R(#eta_{pos})+R(#eta_{neg})) PP", obsDir, profConfig, foldedXTitle);
                allInOneFolded.push_back({pPPFolded, ppConfig});
            }
            if (pToyFolded) {
                DrawComparisonCanvas({{pToyFolded, toyConfig}}, "Canvas_Folded_ToyOnly", fam.familyName + " (R(#eta_{pos})+R(#eta_{neg})) Toy Model", obsDir, profConfig, foldedXTitle);
                allInOneFolded.push_back({pToyFolded, toyConfig});
            }
            if (pMCFolded || pPPFolded || pToyFolded) {
                DrawComparisonCanvas(allInOneFolded, "Canvas_Folded_AllInOne", fam.familyName + " (R(#eta_{pos})+R(#eta_{neg})) Comparisons", obsDir, profConfig, foldedXTitle);
            }

            // --- 3. Subtracted Data - Systematics ---
                // Data minus Data Systematics
            DrawComparisonCanvas(allSubtractedSystematics, "Canvas_Subtracted_Systematics", fam.familyName + " Systematics Difference", obsDir, profConfig, "", subYTitle);
                // Data minus MC Only
            std::vector<ProfileBundle> allInOneSubtracted = allSubtractedSystematics;
            if (pMCSubtracted) {
                DrawComparisonCanvas({{pDataZero, dataConfig}, {pMCSubtracted, mcConfig}}, "Canvas_Subtracted_MCOnly", fam.familyName + " MC Difference", obsDir, profConfig, "", subYTitle);
                allInOneSubtracted.push_back({pMCSubtracted, mcConfig});
            }
                // Data minus pp data
            if (pPPSubtracted) {
                DrawComparisonCanvas({{pDataZero, dataConfig}, {pPPSubtracted, ppConfig}}, "Canvas_Subtracted_PPOnly", fam.familyName + " PP Difference", obsDir, profConfig, "", subYTitle);
                allInOneSubtracted.push_back({pPPSubtracted, ppConfig});
            }
                // Data minus Toy Model Only
            if (pToySubtracted) {
                DrawComparisonCanvas({{pDataZero, dataConfig}, {pToySubtracted, toyConfig}}, "Canvas_Subtracted_ToyOnly", fam.familyName + " Toy Difference", obsDir, profConfig, "", subYTitle);
                allInOneSubtracted.push_back({pToySubtracted, toyConfig});
            }
                // All-In-One
            if (pMCSubtracted || pPPSubtracted || pToySubtracted) {
                DrawComparisonCanvas(allInOneSubtracted, "Canvas_Subtracted_AllInOne", fam.familyName + " All Comparisons Difference", obsDir, profConfig, "", subYTitle);
            }

            // --- Mass Selection Proxy Plots ---
            if (allMassSystematics.size() > 1) { // Only draw if mass files were found
                TDirectory* massDir = obsDir->mkdir("Mass_Selection");
                
                // Standard
                DrawComparisonCanvas(allMassSystematics, "Canvas_MassSelection", fam.familyName + " Mass Selection", massDir, profConfig);
                // Folded
                DrawComparisonCanvas(allFoldedMass, "Canvas_Folded_MassSelection", fam.familyName + " (R(#eta_{pos})+R(#eta_{neg})) Mass Selection", massDir, profConfig, foldedXTitle);
                // Subtracted
                DrawComparisonCanvas(allSubtractedMass, "Canvas_Subtracted_MassSelection", fam.familyName + " Mass Selection Difference", massDir, profConfig, "", subYTitle);
            }

            // --- 4. Individual Comparisons (One-by-One) ---
            // (condensed all of this into a single block of code because the loops become cleaner!)
                // Create a sub-directory specifically for the individual comparisons to keep things organized
            TDirectory* indivDir = obsDir->mkdir("Individual_Comparisons");

            // Loop through all variations (starting at i = 1 to skip the data comparing against itself)
            for (size_t i = 1; i < allSystematics.size(); ++i) {
                std::string sysLabel = allSystematics[i].config.legendLabel;
                std::string sysSuffix = allSystematics[i].config.suffix; // e.g., "_forceDatalikeJet"

                // Standard
                DrawComparisonCanvas({allSystematics[0], allSystematics[i]}, "Canvas_Standard" + sysSuffix, fam.familyName + " Data vs " + sysLabel, indivDir, profConfig);
                // Folded
                DrawComparisonCanvas({allFoldedSystematics[0], allFoldedSystematics[i]}, "Canvas_Folded" + sysSuffix, fam.familyName + " Folded vs " + sysLabel, indivDir, profConfig, foldedXTitle);
                // Subtracted
                DrawComparisonCanvas({allSubtractedSystematics[0], allSubtractedSystematics[i]}, "Canvas_Subtracted" + sysSuffix, fam.familyName + " Difference vs " + sysLabel, indivDir, profConfig, "", subYTitle);
            }

            // Do the same for the MC reference, which was introduced separately
            if (pMC) {
                // Standard
                DrawComparisonCanvas({allSystematics[0], {pMC, mcConfig}}, "Canvas_Standard_MC", fam.familyName + " Data vs MC", indivDir, profConfig);
                // Folded
                DrawComparisonCanvas({allFoldedSystematics[0], {pMCFolded, mcConfig}}, "Canvas_Folded_MC", fam.familyName + " Folded vs MC", indivDir, profConfig, foldedXTitle);
                // Subtracted
                DrawComparisonCanvas({allSubtractedSystematics[0], {pMCSubtracted, mcConfig}}, "Canvas_Subtracted_MC", fam.familyName + " Difference vs MC", indivDir, profConfig, "", subYTitle);
            }
            if (pPP) {
                DrawComparisonCanvas({allSystematics[0], {pPP, ppConfig}}, "Canvas_Standard_PP", fam.familyName + " Data vs PP", indivDir, profConfig);
                DrawComparisonCanvas({allFoldedSystematics[0], {pPPFolded, ppConfig}}, "Canvas_Folded_PP", fam.familyName + " Folded vs PP", indivDir, profConfig, foldedXTitle);
                DrawComparisonCanvas({allSubtractedSystematics[0], {pPPSubtracted, ppConfig}}, "Canvas_Subtracted_PP", fam.familyName + " Difference vs PP", indivDir, profConfig, "", subYTitle);
            }


            // --- 5. Integrated Summaries ---
            TDirectory* integDir = obsDir->mkdir("Integrated_Summary");

            std::vector<std::string> integLabels;
            std::vector<ProfileBundle> integBundles;
            std::vector<ProfileBundle> integSubBundles;
            std::vector<TH1*> integToDelete;

            int numCats = allSystematics.size() + (pMC ? 1 : 0) + (pPP ? 1 : 0) + (pToy ? 1 : 0);
            
                // Compute Integrated Data and save to Grand Summary
            auto dataInteg = GetIntegratedProfile(pData);
            grandDataVals.push_back(dataInteg);
            grandLabels.push_back(profConfig.xAxisTitle);
            if (pMC) grandMCVals.push_back(GetIntegratedProfile(pMC));
            if (pPP) grandPPVals.push_back(GetIntegratedProfile(pPP));
            if (pToy) grandToyVals.push_back(GetIntegratedProfile(pToy));

            // Populate the grand summary for systematics
            // (allSystematics[0] is the data, so sysVariations[i] corresponds to allSystematics[i+1])
            for (size_t i = 0; i < sysVariations.size(); ++i) {
                if (i + 1 < allSystematics.size()) {
                    auto sysInteg = GetIntegratedProfile(dynamic_cast<TProfile*>(allSystematics[i+1].profile));
                    grandSysVals[i].push_back(sysInteg);
                }
            }

                // Loop over variations and construct isolated categorical bins
            for (size_t i = 0; i < allSystematics.size(); ++i) {
                integLabels.push_back(allSystematics[i].config.legendLabel);
                
                auto integ = GetIntegratedProfile(dynamic_cast<TProfile*>(allSystematics[i].profile));
                
                TH1D* hInt = new TH1D(Form("Integ_%zu", i), "", numCats, 0, numCats);
                hInt->SetBinContent(i + 1, integ.first);
                hInt->SetBinError(i + 1, integ.second);
                
                // Performing also the subtraction between integrated observables, variation by variation
                // (error bars are also properly propagated)
                TH1D* hIntSub = new TH1D(Form("IntegSub_%zu", i), "", numCats, 0, numCats);
                if (i == 0) {
                    hIntSub->SetBinContent(i + 1, 0.0);
                    hIntSub->SetBinError(i + 1, 0.0);
                } else {
                    hIntSub->SetBinContent(i + 1, dataInteg.first - integ.first);
                    hIntSub->SetBinError(i + 1, std::sqrt(dataInteg.second*dataInteg.second + integ.second*integ.second));
                }

                integBundles.push_back({hInt, allSystematics[i].config});
                integSubBundles.push_back({hIntSub, allSystematics[i].config});
                integToDelete.push_back(hInt);
                integToDelete.push_back(hIntSub);
            }
            
                // Add MC to integrated list (if available)
            if (pMC) {
                integLabels.push_back(mcConfig.legendLabel);
                auto integMC = GetIntegratedProfile(pMC);
                int idx = allSystematics.size();
                
                TH1D* hIntMC = new TH1D("Integ_MC", "", numCats, 0, numCats);
                hIntMC->SetBinContent(idx + 1, integMC.first);
                hIntMC->SetBinError(idx + 1, integMC.second);
                
                TH1D* hIntSubMC = new TH1D("IntegSub_MC", "", numCats, 0, numCats);
                hIntSubMC->SetBinContent(idx + 1, dataInteg.first - integMC.first);
                hIntSubMC->SetBinError(idx + 1, std::sqrt(dataInteg.second*dataInteg.second + integMC.second*integMC.second));
                
                integBundles.push_back({hIntMC, mcConfig});
                integSubBundles.push_back({hIntSubMC, mcConfig});
                integToDelete.push_back(hIntMC);
                integToDelete.push_back(hIntSubMC);
            }
                // Same for pp data and Toy Model signals
            if (pPP) {
                integLabels.push_back(ppConfig.legendLabel);
                auto integPP = GetIntegratedProfile(pPP);
                int idx = allSystematics.size() + (pMC ? 1 : 0); // Correctly offsets the X-axis placement
                
                TH1D* hIntPP = new TH1D("Integ_PP", "", numCats, 0, numCats);
                hIntPP->SetBinContent(idx + 1, integPP.first);
                hIntPP->SetBinError(idx + 1, integPP.second);
                
                TH1D* hIntSubPP = new TH1D("IntegSub_PP", "", numCats, 0, numCats);
                hIntSubPP->SetBinContent(idx + 1, dataInteg.first - integPP.first);
                hIntSubPP->SetBinError(idx + 1, std::sqrt(dataInteg.second*dataInteg.second + integPP.second*integPP.second));
                
                integBundles.push_back({hIntPP, ppConfig});
                integSubBundles.push_back({hIntSubPP, ppConfig});
                integToDelete.push_back(hIntPP);
                integToDelete.push_back(hIntSubPP);
            }

            if (pToy) {
                integLabels.push_back(toyConfig.legendLabel);
                auto integToy = GetIntegratedProfile(pToy);
                int idx = allSystematics.size() + (pMC ? 1 : 0) + (pPP ? 1 : 0);
                
                TH1D* hIntToy = new TH1D("Integ_Toy", "", numCats, 0, numCats);
                hIntToy->SetBinContent(idx + 1, integToy.first);
                hIntToy->SetBinError(idx + 1, integToy.second);
                
                TH1D* hIntSubToy = new TH1D("IntegSub_Toy", "", numCats, 0, numCats);
                hIntSubToy->SetBinContent(idx + 1, dataInteg.first - integToy.first);
                hIntSubToy->SetBinError(idx + 1, std::sqrt(dataInteg.second*dataInteg.second + integToy.second*integToy.second));
                
                integBundles.push_back({hIntToy, toyConfig});
                integSubBundles.push_back({hIntSubToy, toyConfig});
                integToDelete.push_back(hIntToy);
                integToDelete.push_back(hIntSubToy);
            }

            // 4. Draw them using the custom categorical plotter
            DrawIntegratedCanvas(integBundles, integLabels, "Canvas_Integrated", fam.familyName + " Integrated Summary", integDir, "Integrated " + profConfig.yAxisTitle, false);
            DrawIntegratedCanvas(integSubBundles, integLabels, "Canvas_Integrated_Subtracted", fam.familyName + " Integrated Differences", integDir, "#Delta" + profConfig.yAxisTitle + " (Data - Var)", true);

            // 5. Draw Mass Selection Integrated Summaries
            if (allMassSystematics.size() > 1) {
                std::vector<std::string> massIntegLabels;
                std::vector<ProfileBundle> massIntegBundles;
                std::vector<ProfileBundle> massIntegSubBundles;
                int numMassCats = allMassSystematics.size();
                
                for (size_t i = 0; i < allMassSystematics.size(); ++i) {
                    massIntegLabels.push_back(allMassSystematics[i].config.legendLabel);
                    auto integ = GetIntegratedProfile(dynamic_cast<TProfile*>(allMassSystematics[i].profile));
                    
                    TH1D* hInt = new TH1D(Form("MassInteg_%zu", i), "", numMassCats, 0, numMassCats);
                    hInt->SetBinContent(i + 1, integ.first);
                    hInt->SetBinError(i + 1, integ.second);
                    
                    TH1D* hIntSub = new TH1D(Form("MassIntegSub_%zu", i), "", numMassCats, 0, numMassCats);
                    if (i == 0) {
                        hIntSub->SetBinContent(i + 1, 0.0);
                        hIntSub->SetBinError(i + 1, 0.0);
                    } else {
                        hIntSub->SetBinContent(i + 1, dataInteg.first - integ.first);
                        hIntSub->SetBinError(i + 1, std::sqrt(dataInteg.second*dataInteg.second + integ.second*integ.second));
                    }
                    
                    massIntegBundles.push_back({hInt, allMassSystematics[i].config});
                    massIntegSubBundles.push_back({hIntSub, allMassSystematics[i].config});
                    
                    integToDelete.push_back(hInt);     // Now safely in scope!
                    integToDelete.push_back(hIntSub);
                }
                
                TDirectory* massDir = (TDirectory*)obsDir->Get("Mass_Selection");
                if (massDir) {
                    DrawIntegratedCanvas(massIntegBundles, massIntegLabels, "Canvas_Integrated_MassSelection", fam.familyName + " Mass Selection Integrated", massDir, "Integrated " + profConfig.yAxisTitle, false);
                    DrawIntegratedCanvas(massIntegSubBundles, massIntegLabels, "Canvas_Integrated_Subtracted_MassSelection", fam.familyName + " Mass Selection Integrated Differences", massDir, "#Delta" + profConfig.yAxisTitle + " (Data - Var)", true);
                }
            }

            // Cleanup dynamically allocated profiles for this observable iteration
            for (auto p : integToDelete) delete p;
            for (auto p : profilesToDelete) delete p;
            for (auto p : foldedToDelete) delete p;
            for (auto p : subtractedToDelete) delete p;
        } // end of profiles loop, inside families loop

        // --- 6. Integrated observable Summary (Observable in each variation vs Observable for Data) ---
        // Comparing eta_Jet and eta_Lambda integrations is a way of probing if there were any overflows/underflows
        // and the such (there were! I then modified the eta axis range in the consumer, so this was a good catch!
        // The single-bin integrated observables from other scripts were impervious to this problem though, so good news!),
        // or even other possible problems lingering in the code (none found after the aformentioned fix!)
        int nGrand = grandLabels.size();
        if (nGrand > 0) {
            std::vector<ProfileBundle> grandBundles;
            std::vector<TH1*> grandToDelete; // For memory cleanup
            
            // Add Data
            TH1D* hGrandData = new TH1D("GrandData", "", nGrand, 0, nGrand);
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
                    for (int i = 0; i < nGrand; ++i) {
                        hGrandSys->SetBinContent(i + 1, grandSysVals[sysIdx][i].first);
                        hGrandSys->SetBinError(i + 1, grandSysVals[sysIdx][i].second);
                    }
                    grandBundles.push_back({hGrandSys, sysVariations[sysIdx]});
                    grandToDelete.push_back(hGrandSys);
                }
            }

            // Add MC (if available)
            if (!mcRefDir.empty() && grandMCVals.size() == (size_t)nGrand) { 
                TH1D* hGrandMC = new TH1D("GrandMC", "", nGrand, 0, nGrand);
                for(int i = 0; i < nGrand; ++i) {
                    hGrandMC->SetBinContent(i + 1, grandMCVals[i].first);
                    hGrandMC->SetBinError(i + 1, grandMCVals[i].second);
                }
                grandBundles.push_back({hGrandMC, mcConfig});
                grandToDelete.push_back(hGrandMC);
            }

            // Add pp data data (if available)
            if (!ppRefDir.empty() && grandPPVals.size() == (size_t)nGrand) { 
                TH1D* hGrandPP = new TH1D("GrandPP", "", nGrand, 0, nGrand);
                for(int i = 0; i < nGrand; ++i) {
                    hGrandPP->SetBinContent(i + 1, grandPPVals[i].first);
                    hGrandPP->SetBinError(i + 1, grandPPVals[i].second);
                }
                grandBundles.push_back({hGrandPP, ppConfig});
                grandToDelete.push_back(hGrandPP);
            }

            // Add Toy Model (if available)
            if (!toyModelPath.empty() && grandToyVals.size() == (size_t)nGrand) { 
                TH1D* hGrandToy = new TH1D("GrandToy", "", nGrand, 0, nGrand);
                for(int i = 0; i < nGrand; ++i) {
                    hGrandToy->SetBinContent(i + 1, grandToyVals[i].first);
                    hGrandToy->SetBinError(i + 1, grandToyVals[i].second);
                }
                grandBundles.push_back({hGrandToy, toyConfig});
                grandToDelete.push_back(hGrandToy);
            }
            
            // Draw the Brute Force Canvas
            DrawIntegratedCanvas(grandBundles, grandLabels, "Canvas_BruteForce_Summary", fam.familyName + " Global Summary Across Kinematics", famDir, "Integrated R", false, true);
            
            for (auto p : grandToDelete) delete p;
        }
    } // end of families loop

    fOut->Close();
    delete fOut;
    std::cout << " Done!\n";
}

// To be able to run this both as a compiled C++ macro or a ROOT script:
#ifndef __CINT__
int main(int argc, char** argv) {
    if (argc < 2) { // Check is argc < 2 because this is the bare minimum. mcRefDir, ppRefDir and toyModelPath are all optionals
        std::cerr << "Usage: " << argv[0] << " <consumerDir> [mcRefDir] [ppRefDir] [toyModelPath]\n";
        return 1;
    }
    std::string consumerDir = argv[1];
    std::string mcRefDir = (argc > 2) ? argv[2] : "";
    std::string ppRefDir = (argc > 3) ? argv[3] : "";
    std::string toyModelPath = (argc > 4) ? argv[4] : "";
    
    auxiliarySummaryPlots(consumerDir, mcRefDir, ppRefDir, toyModelPath);
    return 0;
}
#endif