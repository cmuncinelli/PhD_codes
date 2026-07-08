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
// Helper 1.2: Subtract Profiles (Base - Systematics)
// ---------------------------------------------------------
TH1D* SubtractProfiles(TH1* pBase, TH1* pSys, const std::string& newName) {
    if (!pBase || !pSys) return nullptr;
    
    // Safely clone the base binning structure
    TH1D* hSub = nullptr;
    if (auto pBaseProf = dynamic_cast<TProfile*>(pBase)) {
        hSub = pBaseProf->ProjectionX(newName.c_str());
    } else {
        hSub = (TH1D*)pBase->Clone(newName.c_str());
    }
    hSub->Reset(); // Clear original contents and errors
    
    for (int i = 1; i <= hSub->GetNbinsX(); ++i) {
        double vBase = pBase->GetBinContent(i);
        double eBase = pBase->GetBinError(i);
        double vSys = pSys->GetBinContent(i);
        double eSys = pSys->GetBinError(i);
        
        // Ensure both bins contain data before subtracting
        bool baseHasData = (eBase > 1e-9 || std::abs(vBase) > 1e-9 || (dynamic_cast<TProfile*>(pBase) && ((TProfile*)pBase)->GetBinEntries(i) > 0));
        bool sysHasData = (eSys > 1e-9 || std::abs(vSys) > 1e-9 || (dynamic_cast<TProfile*>(pSys) && ((TProfile*)pSys)->GetBinEntries(i) > 0));
        
        if (baseHasData && sysHasData) {
            hSub->SetBinContent(i, vBase - vSys);
            hSub->SetBinError(i, std::sqrt(eBase*eBase + eSys*eSys));
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
        p->SetLineWidth(bundle.config.isBase ? 3 : 2);
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
        // axis actually explains if we are dealing with randJets, datalike estimators, baseline
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

void auxiliaryPlots(const std::string& consumerDir, const std::string& mcRefDir = "") {
    
    // 1. Define Systematic Variations (the list of all useful variations I would like to track into this plot)
    // The base config (empty suffix) is handled separately in the logic to ensure it is always first
    std::vector<VariationConfig> sysVariations = {
        {"_forceRandJet",                 "Rand Jet",                  kBlue,    1, 20, false}, // Full circle
        {"_forceDatalikeJet",             "Data-like Jet",             kRed,     1, 21, false}, // Square
        {"_forceDatalikeJet_10resamples", "Data-like Jet (10 resam.)", kRed,     2, 25, false}, // Open square
        {"_forcePerpToJet",               "Perp to Jet",               kGreen+2, 1, 22, false}, // Triangle up
        {"_forcePerpToJet_10resamples",   "Perp to Jet (10 resam.)",   kGreen+2, 2, 26, false}  // Open triangle up
    };

    VariationConfig baseConfig = {"", "Baseline",    kBlack,  1, 8, true}; // Thick black line. Thickness is controlled by the flag turned "true", essentially
    VariationConfig mcConfig   = {"", "MC Baseline", kCyan+1, 2, 29, false}; // Full stars for MC, with a darker cyan type of color (see TColor and TAttMarker)
                                                                             // Kept line as dashed because MC error bars are (VERY!) large with current statistics

    // 2. Define the Families
    std::vector<FamilyConfig> families = {
        {"Lambda",       "JustLambda"},
        {"AntiLambda",   "JustAntiLambda"},
        {"BothHyperons", "BothHyperons_2GeVLeadP"}
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
    std::string outFilePath = consumerDir + "/AuxiliaryPlots.root";
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
        std::vector<std::pair<double, double>> grandBaseVals;
        std::vector<std::pair<double, double>> grandMCVals;
        std::vector<std::vector<std::pair<double, double>>> grandSysVals(sysVariations.size()); // Track all systematic variations
        std::vector<std::string> grandLabels;

        for (const auto& profConfig : profiles) {
            
            // Create a subdirectory for the observable to keep things extremely tidy
            TDirectory* obsDir = famDir->mkdir(profConfig.profileName.c_str());

            // A. Fetch Baseline Data
            std::string baseFile = consumerDir + "/ConsumerResults_" + fam.baseSuffix + ".root";
            TProfile* pBase = GetProfile(baseFile, profConfig.profileName);
            if (!pBase) continue; // Skip to next observable if baseline is missing

            // B. Fetch Systematics Data
            std::vector<ProfileBundle> allSystematics;
            allSystematics.push_back({pBase, baseConfig}); // Always keep baseline at index 0 (it is the very first thing in this array)
            
            std::vector<TProfile*> profilesToDelete; // Track for memory cleanup
            profilesToDelete.push_back(pBase);

            for (const auto& sys : sysVariations) {
                std::string sysFile = consumerDir + "/ConsumerResults_" + fam.baseSuffix + sys.suffix + ".root";
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
            // SUBTRACTING THE DATA (Baseline - Variations)
            // ---------------------------------------------------------
            std::vector<ProfileBundle> allSubtractedSystematics;
            std::vector<TH1*> subtractedToDelete;

                // 1. Create a zero-line reference for the Baseline bundle (could have made this into just a true vertical line,
            // but in this implementation we actually get some points in the plot, which could look better)
            TH1D* pBaseZero = SubtractProfiles(pBase, pBase, std::string(pBase->GetName()) + "_Subtracted_Base");
            // Force errors to 0 so the baseline acts purely as a flat reference line
            for (int i = 1; i <= pBaseZero->GetNbinsX(); ++i) pBaseZero->SetBinError(i, 0.0);
            
            allSubtractedSystematics.push_back({pBaseZero, baseConfig});
            subtractedToDelete.push_back(pBaseZero);

                // 2. Subtract systematics (start loop at 1 to skip the baseline itself)
            for (size_t i = 1; i < allSystematics.size(); ++i) {
                const auto& bundle = allSystematics[i];
                std::string subName = std::string(pBase->GetName()) + "_Subtracted_" + bundle.config.suffix;
                TH1D* pSub = SubtractProfiles(pBase, bundle.profile, subName);
                
                allSubtractedSystematics.push_back({pSub, bundle.config});
                subtractedToDelete.push_back(pSub);
            }

                // 3. Subtract MC
            TH1D* pMCSubtracted = nullptr;
            if (pMC) {
                pMCSubtracted = SubtractProfiles(pBase, pMC, std::string(pMC->GetName()) + "_Subtracted_MC");
                subtractedToDelete.push_back(pMCSubtracted);
            }

            // Creating the Y-Axis title for the difference plots
            std::string subYTitle = "#Delta" + profConfig.yAxisTitle + " (Base - Var)";

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
            DrawComparisonCanvas({allFoldedSystematics[0]}, "Canvas_Folded_BaselineOnly", fam.familyName + " (R(#eta_{pos})+R(#eta_{neg})) Baseline", obsDir, profConfig, foldedXTitle);
            DrawComparisonCanvas(allFoldedSystematics, "Canvas_Folded_Systematics", fam.familyName + " (R(#eta_{pos})+R(#eta_{neg})) Systematics", obsDir, profConfig, foldedXTitle);
            if (pMCFolded) {
                DrawComparisonCanvas({{pMCFolded, mcConfig}}, "Canvas_Folded_MCOnly", fam.familyName + " (R(#eta_{pos})+R(#eta_{neg})) MC", obsDir, profConfig, foldedXTitle);
                std::vector<ProfileBundle> allInOneFolded = allFoldedSystematics;
                allInOneFolded.push_back({pMCFolded, mcConfig});
                DrawComparisonCanvas(allInOneFolded, "Canvas_Folded_AllInOne", fam.familyName + " (R(#eta_{pos})+R(#eta_{neg})) Comparisons", obsDir, profConfig, foldedXTitle);
            }

            // --- 3. Subtracted Baseline - Systematics ---
                // Variation 1: Baseline minus Data Systematics:
            DrawComparisonCanvas(allSubtractedSystematics, "Canvas_Subtracted_Systematics", fam.familyName + " Systematics Difference", obsDir, profConfig, "", subYTitle);

                // Variation 2: Baseline minus MC Only (and Variation 3: All-In-One):
            if (pMCSubtracted) {
                DrawComparisonCanvas({{pBaseZero, baseConfig}, {pMCSubtracted, mcConfig}}, "Canvas_Subtracted_MCOnly", fam.familyName + " MC Difference", obsDir, profConfig, "", subYTitle);

                std::vector<ProfileBundle> allInOneSubtracted = allSubtractedSystematics;
                allInOneSubtracted.push_back({pMCSubtracted, mcConfig});
                DrawComparisonCanvas(allInOneSubtracted, "Canvas_Subtracted_AllInOne", fam.familyName + " All Comparisons Difference", obsDir, profConfig, "", subYTitle);
            }

            // --- 4. Individual Comparisons (One-by-One) ---
            // (condensed all of this into a single block of code because the loops become cleaner!)
                // Create a sub-directory specifically for the individual comparisons to keep things organized
            TDirectory* indivDir = obsDir->mkdir("Individual_Comparisons");

            // Loop through all variations (starting at i = 1 to skip the baseline comparing against itself)
            for (size_t i = 1; i < allSystematics.size(); ++i) {
                std::string sysLabel = allSystematics[i].config.legendLabel;
                std::string sysSuffix = allSystematics[i].config.suffix; // e.g., "_forceDatalikeJet"

                // Standard
                DrawComparisonCanvas({allSystematics[0], allSystematics[i]}, "Canvas_Standard" + sysSuffix, fam.familyName + " Baseline vs " + sysLabel, indivDir, profConfig);
                // Folded
                DrawComparisonCanvas({allFoldedSystematics[0], allFoldedSystematics[i]}, "Canvas_Folded" + sysSuffix, fam.familyName + " Folded vs " + sysLabel, indivDir, profConfig, foldedXTitle);
                // Subtracted
                DrawComparisonCanvas({allSubtractedSystematics[0], allSubtractedSystematics[i]}, "Canvas_Subtracted" + sysSuffix, fam.familyName + " Difference vs " + sysLabel, indivDir, profConfig, "", subYTitle);
            }

            // Do the same for the MC reference, which was introduced separately
            if (pMC) {
                // Standard
                DrawComparisonCanvas({allSystematics[0], {pMC, mcConfig}}, "Canvas_Standard_MC", fam.familyName + " Baseline vs MC", indivDir, profConfig);
                // Folded
                DrawComparisonCanvas({allFoldedSystematics[0], {pMCFolded, mcConfig}}, "Canvas_Folded_MC", fam.familyName + " Folded vs MC", indivDir, profConfig, foldedXTitle);
                // Subtracted
                DrawComparisonCanvas({allSubtractedSystematics[0], {pMCSubtracted, mcConfig}}, "Canvas_Subtracted_MC", fam.familyName + " Difference vs MC", indivDir, profConfig, "", subYTitle);
            }


            // --- 5. Integrated Summaries ---
            TDirectory* integDir = obsDir->mkdir("Integrated_Summary");

            std::vector<std::string> integLabels;
            std::vector<ProfileBundle> integBundles;
            std::vector<ProfileBundle> integSubBundles;
            std::vector<TH1*> integToDelete;

            int numCats = allSystematics.size() + (pMC ? 1 : 0);
            
                // 1. Compute Integrated Baseline and save to Grand Summary
            auto baseInteg = GetIntegratedProfile(pBase);
            grandBaseVals.push_back(baseInteg);
            grandLabels.push_back(profConfig.xAxisTitle);
            if (pMC) grandMCVals.push_back(GetIntegratedProfile(pMC));

            // Populate the grand summary for systematics
            // (allSystematics[0] is the baseline, so sysVariations[i] corresponds to allSystematics[i+1])
            for (size_t i = 0; i < sysVariations.size(); ++i) {
                if (i + 1 < allSystematics.size()) {
                    auto sysInteg = GetIntegratedProfile(dynamic_cast<TProfile*>(allSystematics[i+1].profile));
                    grandSysVals[i].push_back(sysInteg);
                }
            }

                // 2. Loop over variations and construct isolated categorical bins
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
                    hIntSub->SetBinContent(i + 1, baseInteg.first - integ.first);
                    hIntSub->SetBinError(i + 1, std::sqrt(baseInteg.second*baseInteg.second + integ.second*integ.second));
                }

                integBundles.push_back({hInt, allSystematics[i].config});
                integSubBundles.push_back({hIntSub, allSystematics[i].config});
                integToDelete.push_back(hInt);
                integToDelete.push_back(hIntSub);
            }
            
                // 3. Add MC to integrated list (if available)
            if (pMC) {
                integLabels.push_back(mcConfig.legendLabel);
                auto integMC = GetIntegratedProfile(pMC);
                int idx = allSystematics.size();
                
                TH1D* hIntMC = new TH1D("Integ_MC", "", numCats, 0, numCats);
                hIntMC->SetBinContent(idx + 1, integMC.first);
                hIntMC->SetBinError(idx + 1, integMC.second);
                
                TH1D* hIntSubMC = new TH1D("IntegSub_MC", "", numCats, 0, numCats);
                hIntSubMC->SetBinContent(idx + 1, baseInteg.first - integMC.first);
                hIntSubMC->SetBinError(idx + 1, std::sqrt(baseInteg.second*baseInteg.second + integMC.second*integMC.second));
                
                integBundles.push_back({hIntMC, mcConfig});
                integSubBundles.push_back({hIntSubMC, mcConfig});
                integToDelete.push_back(hIntMC);
                integToDelete.push_back(hIntSubMC);
            }

            // 4. Draw them using the custom categorical plotter
            DrawIntegratedCanvas(integBundles, integLabels, "Canvas_Integrated", fam.familyName + " Integrated Summary", integDir, "Integrated " + profConfig.yAxisTitle, false);
            DrawIntegratedCanvas(integSubBundles, integLabels, "Canvas_Integrated_Subtracted", fam.familyName + " Integrated Differences", integDir, "#Delta" + profConfig.yAxisTitle + " (Base - Var)", true);
            
            for (auto p : integToDelete) delete p;

            // Cleanup dynamically allocated profiles for this observable iteration
            for (auto p : profilesToDelete) delete p;
            for (auto p : foldedToDelete) delete p;
            for (auto p : subtractedToDelete) delete p;
        } // end of profiles loop, inside families loop

        // --- 6. Integrated observable Summary (Observable in each variation vs Observable for Baseline) ---
        // Comparing eta_Jet and eta_Lambda integrations is a way of probing if there were any overflows/underflows
        // and the such (there were! I then modified the eta axis range in the consumer, so this was a good catch!
        // The single-bin integrated observables from other scripts were impervious to this problem though, so good news!),
        // or even other possible problems lingering in the code (none found after the aformentioned fix!)
        int nGrand = grandLabels.size();
        if (nGrand > 0) {
            std::vector<ProfileBundle> grandBundles;
            std::vector<TH1*> grandToDelete; // For memory cleanup
            
            // A. Add Baseline
            TH1D* hGrandBase = new TH1D("GrandBase", "", nGrand, 0, nGrand);
            for(int i = 0; i < nGrand; ++i) {
                hGrandBase->SetBinContent(i + 1, grandBaseVals[i].first);
                hGrandBase->SetBinError(i + 1, grandBaseVals[i].second);
            }
            grandBundles.push_back({hGrandBase, baseConfig});
            grandToDelete.push_back(hGrandBase);

            // B. Add All Systematics
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

            // C. Add MC (if available)
            if (!mcRefDir.empty() && grandMCVals.size() == (size_t)nGrand) { 
                TH1D* hGrandMC = new TH1D("GrandMC", "", nGrand, 0, nGrand);
                for(int i = 0; i < nGrand; ++i) {
                    hGrandMC->SetBinContent(i + 1, grandMCVals[i].first);
                    hGrandMC->SetBinError(i + 1, grandMCVals[i].second);
                }
                grandBundles.push_back({hGrandMC, mcConfig});
                grandToDelete.push_back(hGrandMC);
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
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <consumerDir> [mcRefDir]\n";
        return 1;
    }
    std::string consumerDir = argv[1];
    std::string mcRefDir = (argc > 2) ? argv[2] : "";
    
    auxiliaryPlots(consumerDir, mcRefDir);
    return 0;
}
#endif
