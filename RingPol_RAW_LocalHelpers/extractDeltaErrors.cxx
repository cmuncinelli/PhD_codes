/*
* =====================================================================================
* Script: extractDeltaErrors.cxx
* Author: Cicero Domenico Muncinelli
* * Summary:
* This standalone ROOT post-processing script calculates the correct statistical 
* uncertainties for the Ring Observable ratio <R> = r/n using the Delta Method. 
* It reads the raw accumulator components (r, n, r^2, n^2, r*n) populated by the 
* O2 analysis task, avoiding the mathematical pitfalls of the Standard Error of 
* the Mean (SEM) by preserving the 2*r*n correlation term.
* This file will be extended to not only calculate errors via the delta method, but
* also (possibly) bootstrapping and jackknifing.
*
* Key Operations:
* 1. Reads integrated and 2D differential accumulators from the O2 output.
* 2. Computes <R> and its Delta Method error for the integrated yields.
* 3. Projects 2D differential accumulators (vs pT, Mass, DeltaTheta) into clean 
* 1D histograms of <R> with proper Delta Method error bars.
* 4. Builds a final 4-bin summary histogram ("pRingCuts_Delta") comparing the 
* integrated observable across different kinematic cut families.
* 5. Copies the original SEM-based TProfiles into the output file so both 
* error calculation methods can be directly compared.
*
* Input:  "/path/to/results_consumer/ConsumerResults_<SUFFIX>.root" (Output from the O2 consumer task)
* Output: "<outputDirectory>/ErrorPropagation_<SUFFIX>.root" (outputDirectory is given explicitly
*         as the second argument, and is created recursively if it does not already exist --
*         this keeps post-processing output out of results_consumer/, see run_all_wagons.sh)
*
* Usage:    ./extractDeltaErrors.exe <inputFilePath> <outputDirectory>
* Example:  ./extractDeltaErrors.exe /home/users/cicerodm/RingPol/LHC25ae_pass2_small/UsingTOF_min3ITS/results_consumer/ConsumerResults_BothHyperons.root /home/users/cicerodm/RingPol/LHC25ae_pass2_small/UsingTOF_min3ITS/results_DeltaErr/
* =====================================================================================
*/

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include "TArrayD.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TSystem.h"

// -----------------------------------------------------------------------------
// Helper: Compute Delta Method R and its Error from the 5 sums
// -----------------------------------------------------------------------------
bool ComputeDeltaError(double sum_r, double sum_n, double sum_r2, double sum_n2, double sum_rn, 
                       double& R, double& errR) {
    if (sum_n <= 0) {
        R = 0; errR = 0;
        return false;
    }

    R = sum_r / sum_n;
    
    // Formula: Var(R) = (1 / sum_n^2) * [ sum_r2 + R^2 * sum_n2 - 2*R * sum_rn ]
    double variance = (1.0 / (sum_n * sum_n)) * (sum_r2 + (R * R * sum_n2) - (2.0 * R * sum_rn));
    
    // Protect against floating-point inaccuracies making variance slightly negative
    variance = std::max(0.0, variance);
    errR = std::sqrt(variance);
    
    return true;
}

// -----------------------------------------------------------------------------
// Helper: Convert a 2D accumulator into a 1D histogram of <R>
// -----------------------------------------------------------------------------
TH1D* Process2DTo1D(TH2D* h2d, const std::string& newName) {
    if (!h2d) return nullptr;

    // Create a TH1D using the X-axis of the 2D histogram.
    //
    // The variable- and uniform-binning cases MUST be split. For a uniform axis ROOT stores no edge
    // array at all, so GetXbins()->GetArray() is a null pointer -- and TH1's variable-binning
    // constructor silently accepts that, falling back to fXaxis.Set(nbins, 0, 1). The result is a
    // histogram whose X axis runs 0..1 regardless of the source axis, with no error reported.
    const TArrayD* xBins = h2d->GetXaxis()->GetXbins();
    TH1D* h1d = nullptr;
    if (xBins && xBins->GetSize() > 0) { // Variable binning: copy the edge array verbatim
        h1d = new TH1D(newName.c_str(), h2d->GetTitle(), h2d->GetNbinsX(), xBins->GetArray());
    } else { // Uniform binning: reproduce it from the axis limits
        h1d = new TH1D(newName.c_str(), h2d->GetTitle(), h2d->GetNbinsX(), h2d->GetXaxis()->GetXmin(), h2d->GetXaxis()->GetXmax());
    }
    h1d->GetXaxis()->SetTitle(h2d->GetXaxis()->GetTitle());
    h1d->GetYaxis()->SetTitle("<#it{R}>");

    for (int ix = 1; ix <= h2d->GetNbinsX(); ++ix) {
        // Y-axis bins: 1=r, 2=n, 3=r^2, 4=n^2, 5=rn
        double sum_r  = h2d->GetBinContent(ix, 1);
        double sum_n  = h2d->GetBinContent(ix, 2);
        double sum_r2 = h2d->GetBinContent(ix, 3);
        double sum_n2 = h2d->GetBinContent(ix, 4);
        double sum_rn = h2d->GetBinContent(ix, 5);

        double R = 0, errR = 0;
        if (ComputeDeltaError(sum_r, sum_n, sum_r2, sum_n2, sum_rn, R, errR)) {
            h1d->SetBinContent(ix, R);
            h1d->SetBinError(ix, errR);
        }
    }
    return h1d;
}

// -----------------------------------------------------------------------------
// Main Macro
// -----------------------------------------------------------------------------
int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <inputFilePath> <outputDirectory>\n";
        return 1;
    }
    
    const char* inFileStr = argv[1];
    const char* outDirStr = argv[2];

    // --- Dynamic Output Filename Generation ---
    // The output directory is now taken explicitly from argv[2] (instead of being
    // inferred from the input file's directory), so the filename is derived from
    // the input's basename only.
    std::string inPath(inFileStr);
    std::string filename = inPath;

    // Strip any leading directory components off the input path, keeping only the filename
    size_t lastSlash = inPath.find_last_of('/');
    if (lastSlash != std::string::npos) {
        filename = inPath.substr(lastSlash + 1);
    }

    // Replace "ConsumerResults_" with "ErrorPropagation_"
    std::string prefixToReplace = "ConsumerResults_";
    size_t pos = filename.find(prefixToReplace);
    if (pos != std::string::npos) {
        filename.replace(pos, prefixToReplace.length(), "ErrorPropagation_");
    } else {
        // Fallback: If the prefix isn't found, just prepend to whatever it is
        filename = "ErrorPropagation_" + filename;
    }

    // Normalize the output directory to always end with exactly one trailing slash
    std::string outDir(outDirStr);
    if (outDir.empty() || outDir.back() != '/') outDir += "/";

    // Create the output directory (and any missing parents) if it does not exist yet.
    // kTRUE requests recursive creation, mirroring signalExtractionRing.cxx's convention.
    gSystem->mkdir(outDir.c_str(), kTRUE);

    std::string outFileStr = outDir + filename;

    std::cout << "\n=======================================================\n";
    std::cout << " Starting Delta Method Post-Processing\n";
    std::cout << "=======================================================\n";
    std::cout << " Input File:  " << inPath << "\n";
    std::cout << " Output File: " << outFileStr << "\n";

    // --- Open Files ---
    TFile* inFile = TFile::Open(inPath.c_str(), "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: Could not open input file " << inPath << std::endl;
        return 1;
    }

    TFile* outFile = TFile::Open(outFileStr.c_str(), "RECREATE");
    // Checking the output file's existence before use (e.g., to prevent full-disk errors before they happen):
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error: Could not create output file " << outFileStr << std::endl;
        if (outFile) { outFile->Close(); delete outFile; }
        inFile->Close();
        delete inFile;
        return 1;
    }

    // Only "Ring" is generally booked by the consumer. The other three families have their switches turned off
    // in most cases. "bool mandatory" marks that distinction.
    struct FamilySpec { const char* name; bool mandatory; };
    const std::vector<FamilySpec> families = {
        {"Ring",                      true },
        {"RingKinematicCuts",         false},
        {"JetKinematicCuts",          false},
        {"JetAndLambdaKinematicCuts", false}
    };

    std::vector<std::string> vars2D = {
        "h2dDeltaThetaVsDeltaComp", "h2dLambdaPtVsDeltaComp", "h2dMassVsDeltaComp"
    };

    // Prepare the summary 1D histogram for integrated observables
    TH1D* hRingCutsDelta = new TH1D("pRingCuts_Delta", "Integrated Ring Observable (Delta Method); ;<#it{R}>", 4, 0, 4);
    hRingCutsDelta->GetXaxis()->SetBinLabel(1, "All #Lambda");
    hRingCutsDelta->GetXaxis()->SetBinLabel(2, "p_{T}^{#Lambda}@[0.5,1.5],|y_{#Lambda}|<0.5");
    hRingCutsDelta->GetXaxis()->SetBinLabel(3, "|Jet_{#eta}|<0.5");
    hRingCutsDelta->GetXaxis()->SetBinLabel(4, "#Lambda + Jet cuts");

    // The O2 task's top-level directory:
    const std::string taskDir = "lambdajetpolarizationionsderived/";

    // --- Resolve which families this config actually produced ---
    // Done once, up front, so an absent optional family is reported on a single line instead of once
    // per accumulator it fails to supply:
    std::vector<size_t> presentFamilies; // Indices into families, preserving the summary-bin mapping
    std::string presentList, absentList;
    for (size_t iFam = 0; iFam < families.size(); ++iFam) {
        const bool found = (inFile->Get((taskDir + families[iFam].name).c_str()) != nullptr);
        if (found) {
            presentFamilies.push_back(iFam);
            if (!presentList.empty()) presentList += ", ";
            presentList += families[iFam].name;
        } else if (families[iFam].mandatory) {
            std::cerr << "Error: mandatory family '" << families[iFam].name
                      << "' not found in " << inPath << std::endl;
            outFile->Close(); delete outFile;
            inFile->Close();  delete inFile;
            return 1;
        } else {
            if (!absentList.empty()) absentList += ", ";
            absentList += families[iFam].name;
        }
    }

    std::cout << "\n=======================================================\n";
    std::cout << " Starting Delta Method Post-Processing\n";
    std::cout << "=======================================================\n";
    std::cout << " Families present : " << (presentList.empty() ? "(none)" : presentList) << "\n";
    if (!absentList.empty())
        std::cout << " Families absent  : " << absentList << "  (optional, not enabled in this config)\n";

    for (size_t iPresent = 0; iPresent < presentFamilies.size(); ++iPresent) {
        // iFam indexes the ORIGINAL registry, so the summary histogram keeps its fixed bin labels:
        // bin (iFam + 1) always means the same cut family whether or not its neighbours were produced.
        const size_t iFam = presentFamilies[iPresent];
        const std::string fam = families[iFam].name;
        std::cout << "\n---> Processing Family: " << fam << std::endl;

        // Create directory structure in output file.
        // TDirectoryFile::mkdir returns nullptr when a key of that name already exists, so the chain
        // is checked rather than dereferenced blind.
        TDirectory* dirFam = outFile->mkdir(fam.c_str());
        if (!dirFam) {
            std::cerr << "  [Warning] Could not create output directory '" << fam << "'. Skipping family.\n";
            continue;
        }
        TDirectory* dirDelta = dirFam->mkdir("DeltaMethod_err");
        TDirectory* dirSEM = dirFam->mkdir("SEM_method_err");
        if (!dirDelta || !dirSEM) {
            std::cerr << "  [Warning] Could not create sub-directories under '" << fam << "'. Skipping family.\n";
            continue;
        }

        // ---------------------------------------------------------------------
        // 1. Process Integrated Observable
        // ---------------------------------------------------------------------
        std::string pathInt = taskDir + fam + "/DeltaMethod/hIntegrated";
        TH1D* hInt = (TH1D*)inFile->Get(pathInt.c_str());
        
        if (hInt) {
            // X-axis bins for 1D: 1=r, 2=n, 3=r^2, 4=n^2, 5=rn
            double sum_r  = hInt->GetBinContent(1);
            double sum_n  = hInt->GetBinContent(2);
            double sum_r2 = hInt->GetBinContent(3);
            double sum_n2 = hInt->GetBinContent(4);
            double sum_rn = hInt->GetBinContent(5);

            double R = 0, errR = 0;
            if (ComputeDeltaError(sum_r, sum_n, sum_r2, sum_n2, sum_rn, R, errR)) {
                std::cout << "  [Integrated] <R> = " << R << " +/- " << errR << std::endl;
                // Fill the summary histogram (Bins are 1-indexed)
                hRingCutsDelta->SetBinContent(iFam + 1, R);
                hRingCutsDelta->SetBinError(iFam + 1, errR);
            }
        } else {
            std::cerr << "  [Warning] Could not find " << pathInt << std::endl;
        }

        // ---------------------------------------------------------------------
        // 2. Process 2D Differentials to 1D
        // ---------------------------------------------------------------------
        dirDelta->cd();
        for (const auto& var : vars2D) {
            std::string path2D = taskDir + fam + "/DeltaMethod/" + var;
            TH2D* h2d = (TH2D*)inFile->Get(path2D.c_str());
            
            if (h2d) {
                // Remove "VsDeltaComp" from the new name for clarity
                std::string baseName = var.substr(0, var.find("VsDeltaComp"));
                TH1D* h1d = Process2DTo1D(h2d, baseName + "_DeltaErr");
                if (h1d) { // Adding a check before writing
                    h1d->Write();
                    std::cout << "  [Differential] Processed " << var << " -> " << h1d->GetName() << std::endl;
                } else {
                    std::cerr << "  [Warning] Conversion to 1D failed for " << var << std::endl;
                }
            }
        }

        // // ---------------------------------------------------------------------
        // // 3. Copy Over Original SEM Histograms
        // // ---------------------------------------------------------------------
        // dirSEM->cd();
        // // Look through the input file's QA or base directories for TProfiles to copy
        // // We will do a generic fetch of objects in the family directory
        // std::string famDirPath = taskDir + fam;
        // TDirectory* inFamDir = (TDirectory*)inFile->Get(famDirPath.c_str());
        // if (inFamDir) {
        //     TIter next(inFamDir->GetListOfKeys());
        //     TKey* key;
        //     while ((key = (TKey*)next())) {
        //         TObject* obj = key->ReadObj();
        //         // Check if it's a TProfile or TProfile2D
        //         if (obj->InheritsFrom(TProfile::Class()) || obj->InheritsFrom(TProfile2D::Class())) {
        //             obj->Write(); // writes into dirSEM
        //         }
        //     }
        //     std::cout << "  [SEM] Copied original TProfiles to " << fam << "/SEM_method_err/" << std::endl;
        // }

        // ---------------------------------------------------------------------
        // 3. Copy Over Original SEM Histograms
        // ---------------------------------------------------------------------
        dirSEM->cd();
        
        // Explicit list of the exact TProfiles we want to preserve
        std::vector<std::string> targetSEMProfiles = {
            "pRingObservableDeltaTheta",
            "pRingObservableLambdaPt",
            "pRingObservableMass"
        };
        for (const auto& semName : targetSEMProfiles) {
            std::string fullPath = taskDir + fam + "/" + semName;
            TObject* obj = inFile->Get(fullPath.c_str());
            
            if (obj) {
                obj->Write(); // Safely writes the object with its original name
            } else {
                std::cerr << "  [Warning] SEM profile not found: " << fullPath << std::endl;
            }
        }
        std::cout << "  [SEM] Copied explicit target TProfiles to " << fam << "/SEM_method_err/" << std::endl;
    }

    // ---------------------------------------------------------------------
    // 4. Save Final Summary Histogram
    // ---------------------------------------------------------------------
    outFile->cd();
    hRingCutsDelta->Write();
    std::cout << "\n=======================================================\n";
    std::cout << " Saved summary histogram 'pRingCuts_Delta' to root directory.\n";
    
    // Copy the original pRingCuts TProfile from input to output root dir for direct comparison
    std::string origRingPath = taskDir + "pRingCuts";
    TProfile* pRingCutsOrig = (TProfile*)inFile->Get(origRingPath.c_str());
    if (pRingCutsOrig) {
        pRingCutsOrig->Write("pRingCuts_SEM");
        std::cout << " Copied original 'pRingCuts' as 'pRingCuts_SEM' to root directory.\n";
    }

    std::cout << " Processing Complete! File closed.\n";
    std::cout << "=======================================================\n\n";

    outFile->Close();
    inFile->Close();

    return 0;
}