/*
* ============================================================
* auxiliaryPlots.cxx
* Extracts and compares various 1D TProfiles (e.g., Ring Observable 
* vs Eta of different reference particles) across multiple Consumer 
* configuration outputs within a single wagon's working directory.
*   Generates combined canvases for different analysis families 
* (Lambda, AntiLambda, BothHyperons) and their respective systematic 
* variations (Data-like Jet, Rand Jet, etc.), appropriately scaled.
*
* Usage: root -l -b -q 'auxiliaryPlots.cxx("path/to/wagon/dir", "[optional_MC_ref_path]")'
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
    bool isBase;             // Flag to indicate if this is the main, thick black line
};

// Represents a full family (Lambda, AntiLambda, Both) and its associated variations
struct FamilyConfig {
    std::string familyName;                 // e.g., "Lambda"
    std::string baseSuffix;                 // e.g., "JustLambda"
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
// Helper 1.5: Fold Profile (Sum Positive and Negative Eta)
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
                          const std::string& customXTitle = "") {
    
    if (bundles.empty()) return;

    // 1. Find global Min and Max across all bins and all provided profiles
    double globalMin = 999999.;
    double globalMax = -999999.;
    
    // We assume all profiles share the same X-axis binning structure
    double xMin = bundles[0].profile->GetXaxis()->GetXmin();
    double xMax = bundles[0].profile->GetXaxis()->GetXmax();

    for (const auto& bundle : bundles) {
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

    // Determine the X-Axis title
    std::string xTitle = customXTitle.empty() ? profConfig.xAxisTitle : customXTitle;

    // 3. Draw an empty frame to guarantee proper axes scaling
    TH1* frame = c->DrawFrame(xMin, yMin, xMax, yMax);
    frame->SetTitle(canvasTitle.c_str());
    frame->GetXaxis()->SetTitle(xTitle.c_str());
    frame->GetYaxis()->SetTitle(profConfig.yAxisTitle.c_str());
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);

    // 4. Style everything and populate the legend first
    ProfileBundle baseBundle;
    bool hasBase = false;

    for (auto& bundle : bundles) {
        TH1* p = bundle.profile;
        
        p->SetLineColor(bundle.config.color);
        p->SetMarkerColor(bundle.config.color);
        p->SetLineStyle(bundle.config.lineStyle);
        p->SetMarkerStyle(bundle.config.markerStyle);

        if (bundle.config.isBase) {
            p->SetLineWidth(3); // Thick line for the baseline
            baseBundle = bundle;
            hasBase = true;
        } else {
            p->SetLineWidth(1); // Standard width for systematics
        }
        
        // Add to legend (Baseline will be added first assuming it's passed first in the vector)
        leg->AddEntry(p, bundle.config.legendLabel.c_str(), "pe");
    }

    // 5. Draw the profiles (Systematics first, then Baseline on top)
    for (auto& bundle : bundles) {
        if (!bundle.config.isBase) {
            bundle.profile->Draw("PE SAME");
        }
    }
    
    if (hasBase) { // Drawing the baseline after all others
        baseBundle.profile->Draw("PE SAME");
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
// Main Macro 
// ---------------------------------------------------------

void auxiliaryPlots(const std::string& wagonDir, const std::string& mcRefDir = "") {
    
    // 1. Define Systematic Variations
    // The base config (empty suffix) is handled separately in the logic to ensure it is always first
    std::vector<VariationConfig> sysVariations = {
        {"_forceDatalikeJet",             "Data-like Jet",             kRed,      1, 20, false},
        {"_forceRandJet",                 "Rand Jet",                  kBlue,     1, 21, false},
        {"_forceDatalikeJet_10resamples", "Data-like Jet (10 resam.)", kGreen+2,  2, 24, false},
        {"_forcePerpToJet",               "Perp to Jet",               kMagenta,  1, 22, false},
        {"_forcePerpToJet_10resamples",   "Perp to Jet (10 resam.)",   kOrange+1, 2, 26, false}
    };

    VariationConfig baseConfig = {"", "Baseline", kBlack, 1, 8, true}; // Thick black line
    VariationConfig mcConfig   = {"", "MC Baseline", kGray+2, 2, 25, false}; // Dashed gray for MC

    // 2. Define the Families
    std::vector<FamilyConfig> families = {
        {"Lambda",       "JustLambda"},
        {"AntiLambda",   "JustAntiLambda"},
        {"BothHyperons", "BothHyperons_2GeVLeadP"}
    };

    // 3. Define the Profiles to Extract
    std::vector<ProfileConfig> profiles = {
        {"pRingObservableEtaLeadP",       "#eta_{LeadP}",            "R"},
        {"pRingObservableEtaJet",         "#eta_{Jet}",              "R"},
        {"pRingObservableEta2ndJet",      "#eta_{2ndJet}",           "R"},
        {"pRingObservableEtaLambda",      "#eta_{#Lambda}",          "R"},
        {"pRingObservableEtaLambda2ndJet","#eta_{#Lambda(2ndJet)}",  "R"},
        {"pRingObservableEtaLambdaLeadP", "#eta_{#Lambda(LeadP)}",   "R"}
    };

    // 4. Create the output file
    std::string outFilePath = wagonDir + "/AuxiliaryPlots.root";
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

        for (const auto& profConfig : profiles) {
            
            // Create a subdirectory for the observable to keep things extremely tidy
            TDirectory* obsDir = famDir->mkdir(profConfig.profileName.c_str());

            // A. Fetch Baseline Data
            std::string baseFile = wagonDir + "/ConsumerResults_" + fam.baseSuffix + ".root";
            TProfile* pBase = GetProfile(baseFile, profConfig.profileName);
            if (!pBase) continue; // Skip to next observable if baseline is missing

            // B. Fetch Systematics Data
            std::vector<ProfileBundle> allSystematics;
            allSystematics.push_back({pBase, baseConfig}); // Always keep baseline at index 0 (it is the very first thing in this array)
            
            std::vector<TProfile*> profilesToDelete; // Track for memory cleanup
            profilesToDelete.push_back(pBase);

            for (const auto& sys : sysVariations) {
                std::string sysFile = wagonDir + "/ConsumerResults_" + fam.baseSuffix + sys.suffix + ".root";
                TProfile* pSys = GetProfile(sysFile, profConfig.profileName);
                if (pSys) {
                    allSystematics.push_back({pSys, sys});
                    profilesToDelete.push_back(pSys);
                }
            }

            // C. Fetch MC Baseline Data (if requested)
            TProfile* pMC = nullptr;
            if (!mcRefDir.empty()) {
                std::string mcFile = mcRefDir + "/ConsumerResults_" + fam.baseSuffix + ".root";
                pMC = GetProfile(mcFile, profConfig.profileName);
                if (pMC) profilesToDelete.push_back(pMC);
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

            // Creating the Modulus X-Axis title (e.g. |#eta_{LeadP}|)
            std::string foldedXTitle = "|" + profConfig.xAxisTitle + "|";

            // ---------------------------------------------------------
            // DRAWING THE VARIATIONS
            // ---------------------------------------------------------
            
            // --- 1. Standard (Unfolded) Plots ---
            // Variation 1: Baseline Only
            DrawComparisonCanvas({{pBase, baseConfig}}, "Canvas_BaselineOnly", fam.familyName + " Baseline", obsDir, profConfig);
            // Variation 2: Systematics (+ Baseline)
            DrawComparisonCanvas(allSystematics, "Canvas_Systematics", fam.familyName + " Systematics", obsDir, profConfig);
            if (pMC) {
                // Variation 3: MC Only (if available)
                DrawComparisonCanvas({{pMC, mcConfig}}, "Canvas_MCOnly", fam.familyName + " MC Baseline", obsDir, profConfig);
                std::vector<ProfileBundle> allInOne = allSystematics;
                allInOne.push_back({pMC, mcConfig});
                // Variation 4: All-in-One (if MC available)
                DrawComparisonCanvas(allInOne, "Canvas_AllInOne", fam.familyName + " All Comparisons", obsDir, profConfig);
            }

            // --- 2. Folded Plots ---
            DrawComparisonCanvas({allFoldedSystematics[0]}, "Canvas_Folded_BaselineOnly", fam.familyName + " Folded Baseline", obsDir, profConfig, foldedXTitle);
            DrawComparisonCanvas(allFoldedSystematics, "Canvas_Folded_Systematics", fam.familyName + " Folded Systematics", obsDir, profConfig, foldedXTitle);
            if (pMCFolded) {
                DrawComparisonCanvas({{pMCFolded, mcConfig}}, "Canvas_Folded_MCOnly", fam.familyName + " Folded MC", obsDir, profConfig, foldedXTitle);
                std::vector<ProfileBundle> allInOneFolded = allFoldedSystematics;
                allInOneFolded.push_back({pMCFolded, mcConfig});
                DrawComparisonCanvas(allInOneFolded, "Canvas_Folded_AllInOne", fam.familyName + " Folded Comparisons", obsDir, profConfig, foldedXTitle);
            }

            // Cleanup dynamically allocated profiles for this observable iteration
            for (auto p : profilesToDelete) delete p;
            for (auto p : foldedToDelete) delete p;
        }
    }

    fOut->Close();
    delete fOut;
    std::cout << " Done!\n";
}

// To be able to run this both as a compiled C++ macro or a ROOT script:
#ifndef __CINT__
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <wagonDir> [mcRefDir]\n";
        return 1;
    }
    std::string wagonDir = argv[1];
    std::string mcRefDir = (argc > 2) ? argv[2] : "";
    
    auxiliaryPlots(wagonDir, mcRefDir);
    return 0;
}
#endif
