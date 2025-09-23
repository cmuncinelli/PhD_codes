// This code is based on the lambda_pol_toy_model.cxx file, but reading
// from PYTHIA8313, which should not have any polarization nor any signal
// whatsoever (not only because PYTHIA does not have vorticity, but specially
// because PYTHIA decays everything isotropically!!!)
// I am assuming C++17, at least.
// Start of code: 18/09/2025

// ROOT:
#include "TROOT.h"
#include "TMath.h"
// ROOT objects:
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TH3.h"
#include "TH3D.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <TChain.h> // Joining all trees from all cores as if one giant tree.
// C and C++ libraries:
#include <vector>
#include <array>
#include <cmath>
#include <random>
#include <chrono>
#include <stdexcept>
#include <iostream>
#include <string>
#include <fstream>
// Parallelization capabilities:
#include <omp.h>
// Search for files and patterns (get all output files from a simulation, independently of how many cores it used for running):
#include <dirent.h>
#include <regex>
#include <algorithm>



// Some constants for the code:
const double PI = TMath::Pi();
    // Constant from Phys. Rev. D 110, 030001 - Published 1 August, 2024, i.e.: https://journals.aps.org/prd/pdf/10.1103/PhysRevD.110.030001
const double alpha_H = 0.747; // Fixed as a constant, just for this toy model!
const double mass_Lambda = 1115.683/1000.; // Passed from MeV/c^2 (+/- 0.006) to GeV/c^2
const double mass_proton = 938.27208816/1000.; // Passed from MeV/c^2 (+/- 0.00000029) to GeV/c^2
const double mass_pi = 139.57039/1000.; // Passed from MeV/c^2 (+/- 0.00018), for both pi+ and pi-, to GeV/c^2
const double mL2 = mass_Lambda * mass_Lambda;

void find_jet_proxy();

    // Helper functions:
inline double wrapToInterval(double phi, double phi_min, double phi_max);
void CopyTH3D(const TH3D* source, TH3D* target);
void SqrtHist(TH1* h);
void AddScalar(TH3D* h, double scalar);
void AddScalarIfLargerThanOne(TH3D* h, double scalar);
void SetErrorsFromErrHist(TH1* data, const TH1* errors);
double IntegralAbsNoOverflow(const TH1* h);
int FindClosestBin(TH1D* hist, double targetValue);



int main(int argc, char *argv[]){
    if (argc != 4){
		std::cerr << "Usage:" << argv[0] << " input_card_name N_events target_centrality_bin" << std::endl;
        return 1; // Indicate an error.
    }

	const char *input_card_name = (const char*) argv[1];
    double N_ev_receiver = std::atof(argv[2]); // This receives input like 1e9 and converts it into a proper double.
	uint N_events = static_cast<int>(N_ev_receiver);
    int target_centrality_bin = std::atoi(argv[3]);

    std::cout << "\nReading results for input card " << input_card_name << " and running for " << N_events << " events" << std::endl;

    // Defining the path related variables
    std::string output_folder = "/home/cicero/results/hydro_vorticity/PythiaPol/";
    
    std::string input_folder = "/storage2/cicero/pythia_data/";
    input_folder += input_card_name;
    input_folder += "_" + std::to_string(N_events) + "ev/"; // Specific folder for this input card and its number of events!
    
    std::string centrality_conversion_file_path = input_folder + "Ntracks_to_centrality_conversion.root";

    /////////////////////////////////////////////////////////////////////

    // Getting all the different "output_coreX.root" files that may have been generated (X is the core identifier) for this event:
        // (This method makes the getter independent of the number of N_cores chosen for this specific simulation)
    DIR* dir = opendir(input_folder.c_str());
    if (!dir){
        perror("opendir");
        return 1;
    }

    std::regex pattern(R"(output_core(\d{1,3})\.root)");
    std::vector<std::pair<int,std::string>> file_pair; // A pair where the first number is the core identifier and the second is the full file name

    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr){
        std::string fname = entry->d_name;

        std::smatch match;
        if (std::regex_match(fname, match, pattern)){
            int core_id = std::stoi(match[1]);
            file_pair.emplace_back(core_id, input_folder + fname); // Notice that we are storing the full file path here!
        }
    }
    closedir(dir);

    // Sort by core ID:
    std::sort(file_pair.begin(), file_pair.end());

    // Finally, the loop over all files:
        // Build a TChain of all matching files:
    TChain chain("t3");
    for (auto& current_file_pair : file_pair){ // A loop in the "for i in [...]" style from python! (Yes, I am assuming at least a C++17 standard)
        std::string file_path = current_file_pair.second;
        
        std::cout << "Reading " << file_path << " (core " << current_file_pair.first << ")\n";
        chain.Add(file_path.c_str());
    }

    /////////////////////////////////////////////////////////////////////
    std::cout << "\nInitializing variables and declaring histograms" << std::endl;
    // Now declaring all of the necessary histograms for this operation:
        // Using the exact same binning as what we have from iSS in the lambda_pol_toy_model.cxx, to have an apples-to-apples comparison
    const int N_bins_pT = 30;
    const int N_bins_rap = 10;
    const int N_bins_phi = 100;
    int N_bins_total = N_bins_phi*N_bins_pT*N_bins_rap;

    const double pT_min = 0.;
    const double pT_max = 3.;
    const double rap_max = 0.5;

        // Building the phi intervals:
    double dphi = 2*PI * 1./N_bins_phi; // Bin width
    double phi_min = -dphi/2.; // The lower edge is defined in such a way that phi = 0 is the bin center.
    double phi_max = (2*PI - dphi) + dphi/2.;
    double phi_min_Ring = phi_min - PI;
    double phi_max_Ring = phi_max - PI;

        // Histogram declarations:
        // (3D - pT,y,phi)
    auto hLambdaCounter_pT_y_DeltaPhiJ = new TH3D("hLambdaCounter_pT_y_DeltaPhiJ", "hLambdaCounter_pT_y_DeltaPhiJ", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);
    auto hLambdaAvgDotX_pT_y_DeltaPhiJ = new TH3D("hLambdaAvgDotX_pT_y_DeltaPhiJ", "hLambdaAvgDotX_pT_y_DeltaPhiJ", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);
    auto hLambdaAvgDotY_pT_y_DeltaPhiJ = new TH3D("hLambdaAvgDotY_pT_y_DeltaPhiJ", "hLambdaAvgDotY_pT_y_DeltaPhiJ", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);
    auto hLambdaAvgDotZ_pT_y_DeltaPhiJ = new TH3D("hLambdaAvgDotZ_pT_y_DeltaPhiJ", "hLambdaAvgDotZ_pT_y_DeltaPhiJ", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);

    // auto hLambdaCounterSquared_pT_y_DeltaPhiJ = new TH3D("hLambdaCounterSquared_pT_y_DeltaPhiJ", "hLambdaCounterSquared_pT_y_DeltaPhiJ", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);
        // No need for the squared counter here because there is no weight to each individual particle!
    auto hLambdaAvgDotXSquared_pT_y_DeltaPhiJ = new TH3D("hLambdaAvgDotXSquared_pT_y_DeltaPhiJ", "hLambdaAvgDotXSquared_pT_y_DeltaPhiJ", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);
    auto hLambdaAvgDotYSquared_pT_y_DeltaPhiJ = new TH3D("hLambdaAvgDotYSquared_pT_y_DeltaPhiJ", "hLambdaAvgDotYSquared_pT_y_DeltaPhiJ", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);
    auto hLambdaAvgDotZSquared_pT_y_DeltaPhiJ = new TH3D("hLambdaAvgDotZSquared_pT_y_DeltaPhiJ", "hLambdaAvgDotZSquared_pT_y_DeltaPhiJ", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);

    auto hLambdaPolX_pT_y_DeltaPhiJReco = new TH3D("hLambdaPolX_pT_y_DeltaPhiJReco", "hLambdaPolX_pT_y_DeltaPhiJReco", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);
    auto hLambdaPolY_pT_y_DeltaPhiJReco = new TH3D("hLambdaPolY_pT_y_DeltaPhiJReco", "hLambdaPolY_pT_y_DeltaPhiJReco", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);
    auto hLambdaPolZ_pT_y_DeltaPhiJReco = new TH3D("hLambdaPolZ_pT_y_DeltaPhiJReco", "hLambdaPolZ_pT_y_DeltaPhiJReco", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);

    auto hLambdaPolStarX_pT_y_DeltaPhiJReco = new TH3D("hLambdaPolStarX_pT_y_DeltaPhiJReco", "hLambdaPolStarX_pT_y_DeltaPhiJReco", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);
    auto hLambdaPolStarY_pT_y_DeltaPhiJReco = new TH3D("hLambdaPolStarY_pT_y_DeltaPhiJReco", "hLambdaPolStarY_pT_y_DeltaPhiJReco", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);
    auto hLambdaPolStarZ_pT_y_DeltaPhiJReco = new TH3D("hLambdaPolStarZ_pT_y_DeltaPhiJReco", "hLambdaPolStarZ_pT_y_DeltaPhiJReco", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);

        // (2D - phi,pT)
    auto hLambdaCounter_DeltaPhiJ_pT = new TH2D("hLambdaCounter_DeltaPhiJ_pT", "hLambdaCounter_DeltaPhiJ_pT", N_bins_phi, phi_min, phi_max, N_bins_pT, pT_min, pT_max);
    auto hLambdaPolX_DeltaPhiJ_pTReco = new TH2D("hLambdaPolX_DeltaPhiJ_pTReco", "hLambdaPolX_DeltaPhiJ_pTReco", N_bins_phi, phi_min, phi_max, N_bins_pT, pT_min, pT_max);
    auto hLambdaPolY_DeltaPhiJ_pTReco = new TH2D("hLambdaPolY_DeltaPhiJ_pTReco", "hLambdaPolY_DeltaPhiJ_pTReco", N_bins_phi, phi_min, phi_max, N_bins_pT, pT_min, pT_max);
    auto hLambdaPolZ_DeltaPhiJ_pTReco = new TH2D("hLambdaPolZ_DeltaPhiJ_pTReco", "hLambdaPolZ_DeltaPhiJ_pTReco", N_bins_phi, phi_min, phi_max, N_bins_pT, pT_min, pT_max);

        // (1D - phi)
    auto hLambdaCounter_DeltaPhiJ = new TH1D("hLambdaCounter_DeltaPhiJ", "hLambdaCounter_DeltaPhiJ", N_bins_phi, phi_min, phi_max);
    auto hLambdaCounter_DeltaPhiJRingAngles = new TH1D("hLambdaCounter_DeltaPhiJRingAngles", "hLambdaCounter_DeltaPhiJRingAngles", N_bins_phi, phi_min_Ring, phi_max_Ring);
    auto hProtonStarCounter_DeltaPhiJRing = new TH1D("hProtonStarCounter_DeltaPhiJRing", "hProtonStarCounter_DeltaPhiJRing", N_bins_phi, phi_min_Ring, phi_max_Ring);
            // A variation that includes only the PtCuts counting:
    auto hLambdaCounter_DeltaPhiJ_PtCuts = new TH1D("hLambdaCounter_DeltaPhiJ_PtCuts", "hLambdaCounter_DeltaPhiJ_PtCuts", N_bins_phi, phi_min, phi_max);
    auto hLambdaCounter_DeltaPhiJRingAngles_PtCuts = new TH1D("hLambdaCounter_DeltaPhiJRingAngles_PtCuts", "hLambdaCounter_DeltaPhiJRingAngles_PtCuts", N_bins_phi, phi_min_Ring, phi_max_Ring);
    auto hProtonStarCounter_DeltaPhiJRing_PtCuts = new TH1D("hProtonStarCounter_DeltaPhiJRing_PtCuts", "hProtonStarCounter_DeltaPhiJRing_PtCuts", N_bins_phi, phi_min_Ring, phi_max_Ring);
    
    auto hLambdaPolX_DeltaPhiJReco = new TH1D("hLambdaPolX_DeltaPhiJReco", "hLambdaPolX_DeltaPhiJReco", N_bins_phi, phi_min, phi_max);
    auto hLambdaPolY_DeltaPhiJReco = new TH1D("hLambdaPolY_DeltaPhiJReco", "hLambdaPolY_DeltaPhiJReco", N_bins_phi, phi_min, phi_max);
    auto hLambdaPolZ_DeltaPhiJReco = new TH1D("hLambdaPolZ_DeltaPhiJReco", "hLambdaPolZ_DeltaPhiJReco", N_bins_phi, phi_min, phi_max);
    

    auto hRingObservable_PolStar = new TH1D("hRingObservable_PolStar", "hRingObservable_PolStar", N_bins_phi, phi_min_Ring, phi_max_Ring);
    auto hRingObservable_PolStar_PtCuts = new TH1D("hRingObservable_PolStar_PtCuts", "hRingObservable_PolStar_PtCuts", N_bins_phi, phi_min_Ring, phi_max_Ring);
    auto hRingObservable_PolStar_PtCuts_integrated = new TH1D("hRingObservable_PolStar_PtCuts_integrated", "hRingObservable_PolStar_PtCuts_integrated", 1, phi_min_Ring, phi_max_Ring);
                // A variation calculated with \hat{p}*_D instead of <\hat{p}*_D>:
                // (A second test, where we multiply p_D * p_\Lambda directly, without the averaging -- should make no
                // difference in this Pythia simulation where P_\Lambda = 0, but it is good to work out this structure)
    auto hRingObservable_daughter_no_avg = new TH1D("hRingObservable_daughter_no_avg", "hRingObservable_daughter_no_avg", N_bins_phi, phi_min_Ring, phi_max_Ring);
    auto hRingObservable_daughter_no_avg_PtCuts = new TH1D("hRingObservable_daughter_no_avg_PtCuts", "hRingObservable_daughter_no_avg_PtCuts", N_bins_phi, phi_min_Ring, phi_max_Ring);
    auto hRingObservable_daughter_no_avg_PtCuts_integrated = new TH1D("hRingObservable_daughter_no_avg_PtCuts_integrated", "hRingObservable_daughter_no_avg_PtCuts_integrated", 1, phi_min_Ring, phi_max_Ring);

            // Error propagation histograms:
    auto hLambdaAvgDotX_pT_y_DeltaPhiJ_Err = new TH3D("hLambdaAvgDotX_pT_y_DeltaPhiJ_Err", "hLambdaAvgDotX_pT_y_DeltaPhiJ_Err", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);
    auto hLambdaAvgDotY_pT_y_DeltaPhiJ_Err = new TH3D("hLambdaAvgDotY_pT_y_DeltaPhiJ_Err", "hLambdaAvgDotY_pT_y_DeltaPhiJ_Err", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);
    auto hLambdaAvgDotZ_pT_y_DeltaPhiJ_Err = new TH3D("hLambdaAvgDotZ_pT_y_DeltaPhiJ_Err", "hLambdaAvgDotZ_pT_y_DeltaPhiJ_Err", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);
    auto hLambdaCounter_minus1_pT_y_DeltaPhiJ = new TH3D("hLambdaCounter_minus1_pT_y_DeltaPhiJ", "hLambdaCounter_minus1_pT_y_DeltaPhiJ", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);
        // The scaled version, that has the error on polarization:
    auto hLambdaPolStarX_pT_y_DeltaPhiJReco_Err = new TH3D("hLambdaPolStarX_pT_y_DeltaPhiJReco_Err", "hLambdaPolStarX_pT_y_DeltaPhiJReco_Err", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);
    auto hLambdaPolStarY_pT_y_DeltaPhiJReco_Err = new TH3D("hLambdaPolStarY_pT_y_DeltaPhiJReco_Err", "hLambdaPolStarY_pT_y_DeltaPhiJReco_Err", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);
    auto hLambdaPolStarZ_pT_y_DeltaPhiJReco_Err = new TH3D("hLambdaPolStarZ_pT_y_DeltaPhiJReco_Err", "hLambdaPolStarZ_pT_y_DeltaPhiJReco_Err", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);

    auto hRingObservable_PolStar_Err = new TH1D("hRingObservable_PolStar_Err", "hRingObservable_PolStar_Err", N_bins_phi, phi_min_Ring, phi_max_Ring);
    auto hRingObservable_PolStar_PtCuts_Err = new TH1D("hRingObservable_PolStar_PtCuts_Err", "hRingObservable_PolStar_PtCuts_Err", N_bins_phi, phi_min_Ring, phi_max_Ring);
    auto hRingObservable_PolStar_PtCuts_integrated_Err = new TH1D("hRingObservable_PolStar_PtCuts_integrated_Err", "hRingObservable_PolStar_PtCuts_integrated_Err", 1, phi_min_Ring, phi_max_Ring);


    // Histograms to know about the overall kinematics and production rates of the events:
        // Event counters (in the 40-50% centrality class only):
    auto hNEv_InCentralityClass = new TH1D("hNEv_InCentralityClass", "hNEv_InCentralityClass", 1, 0, 1);
    auto hNEv_WithLambda = new TH1D("hNEv_WithLambda", "hNEv_WithLambda", 1, 0, 1);
    auto hNEv_proton_from_Lambda = new TH1D("hNEv_proton_from_Lambda", "hNEv_proton_from_Lambda", 1, 0, 1); // Should be about 64.1% (the branching ratio of Lambda->p+pi^-)
    auto hNEv_WithLambda_andJets_Pt_leq_2_unit = new TH1D("hNEv_WithLambda_andJets_Pt_leq_2_unit", "hNEv_WithLambda_andJets_Pt_leq_2_unit", 1, 0, 1); // Events that have a particle with pT > 2 GeV/c at least
    auto hNUsefulEv = new TH1D("hNUsefulEv", "hNUsefulEv", 1, 0, 1); // Just to know the number of events. This is not used anywhere in the code and there is no 1./Nev normalization anywhere

        // Jet estimators in the |y| < 2.5 region only:
            // A pT histogram that has all of the jets produced
    auto hHighestMomentum_all_jets = new TH1D("hHighestMomentum_all_jets", "hHighestMomentum_all_jets", 100, 0, 10);
    auto hHighestMomentumPzOverNorm_all_jets = new TH1D("hHighestMomentumPzOverNorm_all_jets", "hHighestMomentumPzOverNorm_all_jets", 100, -1, 1);

            // A pT histogram for the jets that come from an event with Lambda:
    auto hHighestMomentum_EvWithLambda = new TH1D("hHighestMomentum_EvWithLambda", "hHighestMomentum_EvWithLambda", 100, 0, 10);
    auto hHighestMomentumPzOverNorm_EvWithLambda = new TH1D("hHighestMomentumPzOverNorm_EvWithLambda", "hHighestMomentumPzOverNorm_EvWithLambda", 100, -1, 1);

            // A pT histogram for the jets that come from an event with Lambda and have a jet in midrapidity:
    auto hHighestMomentum_WithLambda_Midrap = new TH1D("hHighestMomentum_WithLambda_Midrap", "hHighestMomentum_WithLambda_Midrap", 100, 0, 10); // To know the distribution of the jet proxies' momenta
    auto hHighestMomentumPzOverNorm_WithLambda_Midrap = new TH1D("hHighestMomentumPzOverNorm_WithLambda_Midrap", "hHighestMomentumPzOverNorm_WithLambda_Midrap", 100, -1, 1); // Tells me how much of the momentum of the highest momentum particles is in Z (in normalized vector units)

        // Kinematic information of the produced jets, in the same three categories as above:
    auto hJetProxyY_all_jets = new TH1D("hJetProxyY_all_jets", "hJetProxyY_all_jets", 100, -5, 5);
    auto hJetProxyPt_all_jets = new TH1D("hJetProxyPt_all_jets", "hJetProxyPt_all_jets", 100, 0, 10);
    auto hJetProxyPhi_all_jets = new TH1D("hJetProxyPhi_all_jets", "hJetProxyPhi_all_jets", 100, -PI, PI); // Pythia gives us Phi from -PI to PI!

    auto hJetProxyY_ProtonFromLambda_events = new TH1D("hJetProxyY_ProtonFromLambda_events", "hJetProxyY_ProtonFromLambda_events", 100, -5, 5);
    auto hJetProxyPt_ProtonFromLambda_events = new TH1D("hJetProxyPt_ProtonFromLambda_events", "hJetProxyPt_ProtonFromLambda_events", 100, 0, 10);
    auto hJetProxyPhi_ProtonFromLambda_events = new TH1D("hJetProxyPhi_ProtonFromLambda_events", "hJetProxyPhi_ProtonFromLambda_events", 100, -PI, PI);

    auto hJetProxyY_UsefulEvent = new TH1D("hJetProxyY_UsefulEvent", "hJetProxyY_UsefulEvent", 100, -5, 5);
    auto hJetProxyPt_UsefulEvent = new TH1D("hJetProxyPt_UsefulEvent", "hJetProxyPt_UsefulEvent", 100, 0, 10);
    auto hJetProxyPhi_UsefulEvent = new TH1D("hJetProxyPhi_UsefulEvent", "hJetProxyPhi_UsefulEvent", 100, -PI, PI);

        // To know all of the produced Lambdas:
        // (Do notice that these DO NOT include the coordinate axis rotation that shifts phi=0 to phi_{jet}!)
    auto hLambdaPt = new TH1D("hLambdaPt", "hLambdaPt", 100, 0, 10);
    auto hLambdaY = new TH1D("hLambdaY", "hLambdaY", 100, -5, 5);
    auto hLambdaPhi = new TH1D("hLambdaPhi", "hLambdaPhi", 100, -PI, PI); // Pythia gives phi in the range -PI to PI!
    auto hLambdaDeltaBetaJ_WithProtonFromLambda = new TH1D("hLambdaDeltaBetaJ_WithProtonFromLambda", "hLambdaDeltaBetaJ_WithProtonFromLambda", 100, 0, PI); // A histogram that tells us the angular distribution of Lambdas along a jet (this is the angle from the dot product, so more general than Phi)
        // Another histogram, which includes all possible Lambda decay modes:
    auto hLambdaDeltaBetaJ_AllLambdaDecayModes = new TH1D("hLambdaDeltaBetaJ_AllLambdaDecayModes", "hLambdaDeltaBetaJ_AllLambdaDecayModes", 100, 0, PI);

        // To know all of Lambdas that are also proton mothers:
    auto hLambdaPt_proton_mothers = new TH1D("hLambdaPt_proton_mothers", "hLambdaPt_proton_mothers", 100, 0, 10);
    auto hLambdaY_proton_mothers = new TH1D("hLambdaY_proton_mothers", "hLambdaY_proton_mothers", 100, -5, 5);
    auto hLambdaPhi_proton_mothers = new TH1D("hLambdaPhi_proton_mothers", "hLambdaPhi_proton_mothers", 100, -PI, PI);

        // Turning on Sumw2() for the histograms that only store counts:
    hHighestMomentum_all_jets->Sumw2();
    hHighestMomentumPzOverNorm_all_jets->Sumw2();

    hHighestMomentum_EvWithLambda->Sumw2();
    hHighestMomentumPzOverNorm_EvWithLambda->Sumw2();

    hHighestMomentum_WithLambda_Midrap->Sumw2();
    hHighestMomentumPzOverNorm_WithLambda_Midrap->Sumw2();

    hJetProxyY_all_jets->Sumw2();
    hJetProxyPt_all_jets->Sumw2();
    hJetProxyPhi_all_jets->Sumw2();
    hJetProxyY_ProtonFromLambda_events->Sumw2();
    hJetProxyPt_ProtonFromLambda_events->Sumw2();
    hJetProxyPhi_ProtonFromLambda_events->Sumw2();
    hJetProxyY_UsefulEvent->Sumw2();
    hJetProxyPt_UsefulEvent->Sumw2();
    hJetProxyPhi_UsefulEvent->Sumw2();

    hLambdaPt->Sumw2();
    hLambdaY->Sumw2();
    hLambdaPhi->Sumw2();
    hLambdaDeltaBetaJ_WithProtonFromLambda->Sumw2();

    hLambdaDeltaBetaJ_AllLambdaDecayModes->Sumw2();

    hLambdaPt_proton_mothers->Sumw2();
    hLambdaY_proton_mothers->Sumw2();
    hLambdaPhi_proton_mothers->Sumw2();
    
    
    // In case you want to speed the code execution up:
    //     // Buffer vector declarations:
    // std::vector<double> hLambdaCounter_pT_y_DeltaPhiJ_bufferGlobal(N_bins_total, 0.0);
    // std::vector<double> hLambdaAvgDotX_pT_y_DeltaPhiJ_bufferGlobal(N_bins_total, 0.0);
    // std::vector<double> hLambdaAvgDotY_pT_y_DeltaPhiJ_bufferGlobal(N_bins_total, 0.0);
    // std::vector<double> hLambdaAvgDotZ_pT_y_DeltaPhiJ_bufferGlobal(N_bins_total, 0.0);
    // std::vector<double> hLambdaCounterSquared_pT_y_DeltaPhiJ_bufferGlobal(N_bins_total, 0.0);
    // std::vector<double> hLambdaAvgDotXSquared_pT_y_DeltaPhiJ_bufferGlobal(N_bins_total, 0.0);
    // std::vector<double> hLambdaAvgDotYSquared_pT_y_DeltaPhiJ_bufferGlobal(N_bins_total, 0.0);
    // std::vector<double> hLambdaAvgDotZSquared_pT_y_DeltaPhiJ_bufferGlobal(N_bins_total, 0.0);

    // // std::vector<double> hLambdaCounter_pT_y_bufferGlobal(N_bins_pT*N_bins_rap, 0.0); // Not useful for the current study
    // std::vector<double> hLambdaPolX_DeltaPhiJ_pT_bufferGlobal(N_bins_phi*N_bins_pT, 0.0);
    // std::vector<double> hLambdaPolY_DeltaPhiJ_pT_bufferGlobal(N_bins_phi*N_bins_pT, 0.0);
    // std::vector<double> hLambdaPolZ_DeltaPhiJ_pT_bufferGlobal(N_bins_phi*N_bins_pT, 0.0);
    // std::vector<double> hLambdaCounter_DeltaPhiJ_pT_bufferGlobal(N_bins_phi*N_bins_pT, 0.0);
    // std::vector<double> hLambdaCounter_DeltaPhiJ_bufferGlobal(N_bins_phi, 0.0);

    // std::vector<double> hLambdaCounter_DeltaPhiJRingAngles_bufferGlobal(N_bins_phi, 0.0);
    // std::vector<double> hLambdaCounter_DeltaPhiJ_PtCuts_bufferGlobal(N_bins_phi, 0.0);
    // std::vector<double> hLambdaCounter_DeltaPhiJRingAngles_PtCuts_bufferGlobal(N_bins_phi, 0.0);
    // std::vector<double> hProtonStarCounter_DeltaPhiJRing_PtCuts_bufferGlobal(N_bins_phi, 0.0);

    // std::vector<double> hLambdaPolX_DeltaPhiJ_bufferGlobal(N_bins_phi, 0.0);
    // std::vector<double> hLambdaPolY_DeltaPhiJ_bufferGlobal(N_bins_phi, 0.0);
    // std::vector<double> hLambdaPolZ_DeltaPhiJ_bufferGlobal(N_bins_phi, 0.0);
    // std::vector<double> hLambdaPolX_pT_y_DeltaPhiJ_bufferGlobal(N_bins_total, 0.0);
    // std::vector<double> hLambdaPolY_pT_y_DeltaPhiJ_bufferGlobal(N_bins_total, 0.0);
    // std::vector<double> hLambdaPolZ_pT_y_DeltaPhiJ_bufferGlobal(N_bins_total, 0.0);

    // std::vector<double> hRingObservable_PolStar_PtCuts_bufferGlobal(N_bins_phi, 0.0);
    // std::vector<double> hRingObservable_PolStar_PtCuts_integrated_bufferGlobal(N_bins_phi, 0.0);
    // std::vector<double> hRingObservable_daughter_no_avg_PtCuts_bufferGlobal(N_bins_phi, 0.0);
    // std::vector<double> hRingObservable_daughter_no_avg_PtCuts_integrated_bufferGlobal(N_bins_phi, 0.0);

    /////////////////////////////////////////////////////////////////////
    std::cout << "Estimating centrality" << std::endl;
    // Centrality class definitions (This code comes from root_to_csv_custom_PYTHIA.cxx from the "laguerre_project" repository):
        // (The target class for comparison with Hydro is 40-50%)
    std::array<double, 10> cent_bin_limits{0, 5, 10, 20, 30, 40, 50, 60, 70, 80}; // Edited from 10 to 9 centrality classes.
    std::array<double, 11> cent_bin_limits_full{0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90};
        // Declaring the approximator for the centrality class limits (the centrality binning is very precise for large N_events, but we won't always get 5.000%,
        // we can get 4.98% as the finest possible limit with Ntracks_final_charged being a discrete scale)
    double cent_bin_limits_non_approximate[10] = {}; // The first entry will always be zero, but the rest will change! (and should be close to cent_bin_limits)
    std::array<int, 9> multiplicity_borders{}; // Initialized all of them with zeroes! --> This is the vector that we can use to see in which class an event fits
    
    TFile centrality_conversion_file(centrality_conversion_file_path.c_str());
    TH1D *centrality_conversion_histogram = (TH1D *) centrality_conversion_file.Get("hNtracks_to_centrality_final_charged"); // Using the number of charged final particles as the estimator.

    // Now approximating the multiplicity limits to match the desired percentages from cent_bin_limits
        // I emphasize: APPROXIMATE! For example, it is possible that the Ntracks observable has only enough resolution to give you 0-9.7% instead of the desired 0-10%
    for (int cent_limit_idx = 1; cent_limit_idx < cent_bin_limits.size(); cent_limit_idx++){ // Starts from the second limit, as the first one is always 0% and is not needed!!
        int current_border = FindClosestBin(centrality_conversion_histogram, cent_bin_limits[cent_limit_idx]);
        multiplicity_borders[cent_limit_idx - 1] = current_border;

        // Also collecting the actual multiplicity percentage that was achieved, not the approximate value (This will be useful for the Nev scaling later):
        double non_approximate_multiplicity_percentage = centrality_conversion_histogram->GetBinContent(current_border);
        cent_bin_limits_non_approximate[cent_limit_idx] = non_approximate_multiplicity_percentage;

        std::cout << "The ideal multiplicity class limit was " << std::to_string(cent_bin_limits[cent_limit_idx]) << " while the actual value achieved was "
                    << std::to_string(non_approximate_multiplicity_percentage) << " (with a multiplicity of " << current_border << ")"
                    << std::endl;
        std::cout << "  The error for this multiplicity class definition is: (" 
                    << (-1) * (1 - non_approximate_multiplicity_percentage * 1./cent_bin_limits[cent_limit_idx]) * 100 << "%)" // Multiplied by -1 to have the (estimate - value)/value error, not the (value - estimate)/value error.
                    << std::endl;
    }

    std::cout << "\nThe Ntracks_final_charged multiplicity cuts for centrality definition in each class are:" << std::endl;
    for (int i = 0; i < 9; i++){
        std::cout << multiplicity_borders[i] << std::endl;
    }

    /////////////////////////////////////////////////////////////////////
    std::cout << "\n\nStarting event loop procedure" << std::endl;
    // Variables to store the useful quantities of each event:
        // Actually don't need them for calculating the polarization -- The averaging will already introduce binning and will not need the full information
    std::vector<TVector3> jet_trimomenta_array; // The trimomentum of the jet, to determine \hat{t}
    // std::vector<Float_t> jet_phi_array; // To avoid recalculating, the azimuthal angle of the jet
    // std::vector<std::vector<TVector3>> proton_trimomenta_matrix; // First dimension is the number of events, second one is the particle idx
    // std::vector<std::vector<TVector3>> proton_trimomenta_Star_matrix;
    // std::vector<std::vector<TVector3>> lambda_trimomenta_matrix;


    // Finally, doing the loop on each event of the TTree (similar to the loop in root_to_csv_custom_PYTHIA.cxx):
    Long64_t nEntries = chain.GetEntries();
    std::cout << "Total number of events in the TChain is: " << nEntries << "\n";

        // Indexes for looping over all particles of the event:
    // const Int_t kMaxTrack = 7000; // This is no longer necessary with the update on PythiaGenMin.cxx, which uses C++ vectors, with no specific number of entries!
    Int_t ntrack; // The number of particles generated in that event (including carbon-copies and non-final state particles)
        // Centrality selection:
    Int_t Ntracks_final_charged;
        // PID selection:
    // Int_t PID[kMaxTrack];
    std::vector<Int_t> *PID_pointer = nullptr; // Newer C++ implementation --> As ROOT gives me a pointer for the std::vector from the way PythiaGenMin is coded, this has to also be a pointer
    std::vector<Int_t> *mother1_Idx_pointer = nullptr; // Not quite the PID of the mother. It is the Idx of the mother in Pythia's list of particles
    std::vector<Int_t> *mother2_Idx_pointer = nullptr; // Also initializing the pointer as a nullptr. If you don't, you will get a segmentation error!
    std::vector<Bool_t> *IsFinal_pointer = nullptr;
    std::vector<Bool_t> *IsCarbonCopy_pointer = nullptr;
        // 4-moment information:
    std::vector<Float_t> *m_pointer = nullptr; // The actual simulated mass of the particle in the dynamic, so may be somewhat different from the PDG value.
    std::vector<Float_t> *pt_pointer = nullptr;
    std::vector<Float_t> *y_pointer = nullptr;
    std::vector<Float_t> *Phi_pointer = nullptr;
            // If needed, can reconstruct these variables from (m, pT, y, Phi):
    // std::vector<Float_t> *px_pointer;
    // std::vector<Float_t> *py_pointer;
    // std::vector<Float_t> *pz_pointer;
    // std::vector<Float_t> *E_pointer;


        // Connecting these variables to the TChain:
    chain.SetBranchAddress("ntrack", &ntrack);
    chain.SetBranchAddress("Ntracks_final_charged", &Ntracks_final_charged);

    // chain.SetBranchAddress("ID", &PID); // Do notice that with this new implementation, based on std::vector, you need to pass the address of the object: we are no longer dealing with C arrays (pointers!)
    chain.SetBranchAddress("ID", &PID_pointer); // Do notice that with this new implementation, based on std::vector, you need to pass the std::vector address, thus the previous usage of a pointer!
    chain.SetBranchAddress("Mother1_ID", &mother1_Idx_pointer); // Weird reading that still needs an extra & for the address.
    chain.SetBranchAddress("Mother2_ID", &mother2_Idx_pointer);
    chain.SetBranchAddress("IsFinal", &IsFinal_pointer);
    chain.SetBranchAddress("IsCarbonCopy", &IsCarbonCopy_pointer);

    chain.SetBranchAddress("m", &m_pointer);
    chain.SetBranchAddress("pt", &pt_pointer);
    chain.SetBranchAddress("y", &y_pointer);
    chain.SetBranchAddress("Phi", &Phi_pointer);

    // Int_t N_useful_events = 0; // This number should be updated if there actually was any proton with a Lambda mother in the event, and if the jet was in the |y| < 0.5 range.
    for (Long64_t iEvent = 0; iEvent < nEntries; iEvent++){
        // Keeping track of the proccess for every 5% of the events:
        int five_percent_step = (N_events * 0.05); // Using the same low-statistics-safe (N_events < 20 safe) approach, so that int (N_events * 0.05) does not give us an undefined "% 0" operation.
        if (five_percent_step > 0 ? (iEvent % five_percent_step == 0) : false){ 
            std::cout << "Now on event " << iEvent << " of " << N_events << " (" << (100.0 * iEvent / N_events) << "%)" << std::endl;
        }

        chain.GetEntry(iEvent);

        // Dereferencing for the inner loop -- This new C++-like implementation with std::vector* is kinda messy:
            // In other words, copying the contents of the pointer into local variables
        std::vector<Int_t> PID = *PID_pointer;
        std::vector<Int_t> mother1_Idx = *mother1_Idx_pointer;
        std::vector<Int_t> mother2_Idx = *mother2_Idx_pointer;
        std::vector<Bool_t> IsFinal = *IsFinal_pointer;
        std::vector<Bool_t> IsCarbonCopy = *IsCarbonCopy_pointer;
        std::vector<Float_t> m = *m_pointer;
        std::vector<Float_t> pt = *pt_pointer;
        std::vector<Float_t> y = *y_pointer;
        std::vector<Float_t> Phi = *Phi_pointer;

        // Will use something like this to identify to which centrality class the event belongs to:
        int centrality_bin = 1;
        while (Ntracks_final_charged < multiplicity_borders[centrality_bin - 1]){
            if(centrality_bin == 10){
                centrality_bin += 1; // Don't actually need this, do I?
                break;
            }
            centrality_bin += 1;
        }

        // For now will analyze only the events in the 40-50% centrality class, that lie in the 5th centrality bin:
        if (centrality_bin != target_centrality_bin){continue;} // target_centrality_bin should be 1 for 0-5%, 6 for 40-50%
        hNEv_InCentralityClass->Fill(0); // Just to get a counter in the same .root file that has the number of total events in the centrality class (should either be 10% of total or 5% of total)

        // Declaring vectors for the particle loop (they won't need to be used elsewhere!):
        std::vector<TLorentzVector> proton_4momenta;
        // std::vector<TLorentzVector> proton_4momenta_Star;
        std::vector<TLorentzVector> lambda_4momenta;
        std::vector<TLorentzVector> lambda_4momenta_all_decay_modes; // Another vector, that stores all possible decays, even if they didn't become the useful p+pi^- for the analysis
        std::vector<Float_t> lambda_phi_array; // Just for the calculations done in this event's loop
        std::vector<Float_t> lambda_y_array;
        std::vector<Float_t> lambda_pT_array;

        // Entering the particle loop:
        Float_t highest_particle_momentum = 0;
        Int_t highest_particle_momentum_Idx = -1; // Saving this to be able to recover
        
            // Checks for increasingly restrictive conditions necessary for the event to be useful:
        Bool_t contains_lambda = false;
        // Bool_t jet_Pt_leq_2 = false; // These two can be checked after the event
        // Bool_t jet_Pt_leq_2_contains_Lambda = false;
        Bool_t contains_particle = false; // Contains a proton that came from a Lambda
        for (int particle_idx = 0; particle_idx < ntrack; particle_idx++){
            // Excluding particles that are useless for this comparison:
            if (IsCarbonCopy[particle_idx]){continue;} // Used this exclusion first, because I also want to analyze Lambda's that did not decay into protons:
            // Getting useful information for the jet-proxy estimation and the boosts:
            Float_t particle_y = y[particle_idx];
            Float_t particle_pT = pt[particle_idx];
            Float_t particle_phi = Phi[particle_idx];
            Float_t particle_m = m[particle_idx];

            if (PID[particle_idx] == 3122){
                hLambdaY->Fill(particle_y);
                hLambdaPt->Fill(particle_pT);
                hLambdaPhi->Fill(particle_phi);

                // Saving the 4 vector for this case that didn't necessarily decay into p+pi^-:
                // (Yes, I am recalculating the px,py,pz,E, but there are so few Lambda that this shouldn't matter)
                Float_t particle_px = particle_pT*std::cos(particle_phi);
                Float_t particle_py = particle_pT*std::sin(particle_phi);
                Float_t particle_mT = std::sqrt(particle_m*particle_m + particle_pT*particle_pT);
                Float_t particle_pz = particle_mT*std::sinh(particle_y);
                Float_t particle_momentum_squared = particle_px*particle_px + particle_py*particle_py + particle_pz*particle_pz;
                Float_t particle_E = std::sqrt(particle_m*particle_m + particle_momentum_squared);

                TLorentzVector lambda_4_momentum_all_decay_modes(particle_px, particle_py, particle_pz, particle_E);
                lambda_4momenta_all_decay_modes.push_back(lambda_4_momentum_all_decay_modes);

                contains_lambda = true;
            }

                // Now continuing the selection to see only the final particles
            if (!IsFinal[particle_idx]){continue;}
            // std::cout<<"Actually useful particle!!!" << std::endl;
                // Recovering (px,py,pz) to calculate the norm of the trimomentum:
            Float_t particle_px = particle_pT*std::cos(particle_phi);
            Float_t particle_py = particle_pT*std::sin(particle_phi);
            Float_t particle_mT = std::sqrt(particle_m*particle_m + particle_pT*particle_pT);
            Float_t particle_pz = particle_mT*std::sinh(particle_y);
            Float_t particle_momentum_squared = particle_px*particle_px + particle_py*particle_py + particle_pz*particle_pz; // Defined here to avoid recalculating stuff for particle_E and particle_momentum
            Float_t particle_E = std::sqrt(particle_m*particle_m + particle_momentum_squared);
            
            // Updating the proxy to see if this particle has the highest momentum of the produced particles (jet direction!):
            // (Notice that this comes BEFORE excluding this particle if it is not a proton! It can be a jet proxy, yet not be the particle that came from the Lambda decay!)
            Float_t particle_momentum = std::sqrt(particle_momentum_squared);
            if (std::fabs(particle_y) < 2.5 && (highest_particle_momentum < particle_momentum)){
                highest_particle_momentum = particle_momentum;
                highest_particle_momentum_Idx = particle_idx;
                // std::cout << "Particle momentum is " << highest_particle_momentum << std::endl;
            }

            //////////////////////////////////////
            //// Calculating the polarization ////
            //////////////////////////////////////
            if (PID[particle_idx] != 2212){continue;} // Will only use protons for the end of this loop, which is the part that actually calculates the polarization
                // Do notice that only PROTONS are useful here: pbar can be used for lambda-bar polarization studies, but in the first version of this code
                // the goal is to use ONLY lambda, not lambda-bar.

            int mother_idx = mother1_Idx[particle_idx];
            while (IsCarbonCopy[mother_idx]){ // Should keep searching for the true mother until it is an actual particle, not a carbon-copy
                mother_idx = mother1_Idx[mother_idx]; // Will search for this mother's mother untill it is no longer a carbon-copy
            }
            if (PID[mother_idx] != 3122){continue;} // This means that the mother of this particle IS NOT a Lambda. So excludes it.
            contains_particle = true;
            // std::cout << "This proton has a Lambda mother!" << std::endl; // Very few events should have this property!

            // Fetching the mother's kinematic variables:
            // (Necessary for the Ring Observable calculation)
            Float_t lambda_y = y[mother_idx];
            Float_t lambda_pT = pt[mother_idx];
            Float_t lambda_phi = Phi[mother_idx];
            Float_t lambda_m = m[mother_idx];
                // For the boost:
            Float_t lambda_px = lambda_pT*std::cos(lambda_phi);
            Float_t lambda_py = lambda_pT*std::sin(lambda_phi);
            Float_t lambda_mT = std::sqrt(lambda_m*lambda_m + lambda_pT*lambda_pT);
            Float_t lambda_pz = lambda_mT*std::sinh(lambda_y);
            
            Float_t lambda_momentum_squared = lambda_px*lambda_px + lambda_py*lambda_py + lambda_pz*lambda_pz;
            Float_t lambda_E = std::sqrt(lambda_m*lambda_m + lambda_momentum_squared);

                // Filling some useful histograms about the Lambdas that decayed into protons:
            hLambdaY_proton_mothers->Fill(lambda_y);
            hLambdaPt_proton_mothers->Fill(lambda_pT);
            hLambdaPhi_proton_mothers->Fill(lambda_phi);

            // Constructing the 4-momenta for the boost:
                // (As I am later rotating the coordinate system to have phi_jet = 0, I will only boost the particles after rotating them in the laboratory frame!
                // Remember that rotation is not Lorentz-invariant, i.e., rotate in lab + Boost is different from Boost + rotate in Lambda rest frame!)
            TLorentzVector lambda_4momentum(lambda_px, lambda_py, lambda_pz, lambda_E);
            TLorentzVector proton_4momentum(particle_px, particle_py, particle_pz, particle_E);
            // TLorentzVector proton_4momentum_star(particle_px, particle_py, particle_pz, particle_E); // Defined from the 4-momentum in the lab frame, but the boost method acts in-place,
                                                                                                         // so had to define it from the laboratory's 4-momentum.

            lambda_4momenta.push_back(lambda_4momentum);
            proton_4momenta.push_back(proton_4momentum);
            lambda_phi_array.push_back(lambda_phi);
            lambda_y_array.push_back(lambda_y);
            lambda_pT_array.push_back(lambda_pT);
        }
        // Getting the jet variables:
        Float_t jet_phi = Phi[highest_particle_momentum_Idx];
        Float_t jet_px = pt[highest_particle_momentum_Idx]*std::cos(jet_phi);
        Float_t jet_py = pt[highest_particle_momentum_Idx]*std::sin(jet_phi);
        Float_t jet_mT = std::sqrt(m[highest_particle_momentum_Idx]*m[highest_particle_momentum_Idx]
                                   + pt[highest_particle_momentum_Idx]*pt[highest_particle_momentum_Idx]);
        Float_t jet_pz = jet_mT*std::sinh(y[highest_particle_momentum_Idx]);

            // Saving the jet distribution before the coordinate shift:
        hJetProxyY_all_jets->Fill(y[highest_particle_momentum_Idx]);
        hJetProxyPt_all_jets->Fill(pt[highest_particle_momentum_Idx]);
        hJetProxyPhi_all_jets->Fill(jet_phi);

        TVector3 jet_trimomentum(jet_px, jet_py, jet_pz);
        // Rotating the jet to a coordinate system where \hat{t}_x = 1:
        jet_trimomentum.SetPhi(0);
        // Getting the unit vector for Z axis projections:
        TVector3 jet_trimomentum_unit = jet_trimomentum.Unit();

        hHighestMomentum_all_jets->Fill(highest_particle_momentum);
        hHighestMomentumPzOverNorm_all_jets->Fill(jet_trimomentum_unit.Z()); // How much of the unit vector is in z. Includes all jets with |y|<2.5

        // Filling event counter:
        if (contains_lambda){
            hNEv_WithLambda->Fill(0);

            // Another DeltaBetaJ histogram, that works for events where the Lambda decayed into something different from a proton (the other decay mode that differs p+pi^-)
            for (int i = 0; i < lambda_4momenta_all_decay_modes.size(); i++){ // Notice the usage of lambda_4momenta_all_decay_modes, which is filled even if the Lambda didn't decay into a p+pi^-
                TLorentzVector lambda_4momentum = lambda_4momenta_all_decay_modes[i];
                TVector3 lambda_trivec = lambda_4momentum.Vect();
                Double_t DeltaBetaJ = jet_trimomentum.Angle(lambda_trivec);
                hLambdaDeltaBetaJ_AllLambdaDecayModes->Fill(DeltaBetaJ);
            }
        }

        if(!contains_particle){continue;} // Skip the whole event: it doesn't have a single proton that came from Lambda.
        hNEv_proton_from_Lambda->Fill(0); // This event passed the "contains_particle" check
        hHighestMomentum_EvWithLambda->Fill(highest_particle_momentum);
        hHighestMomentumPzOverNorm_EvWithLambda->Fill(jet_trimomentum_unit.Z()); // How much of the unit vector is in z. Includes all jets with |y|<2.5, but only those with a Lambda

        hJetProxyY_ProtonFromLambda_events->Fill(y[highest_particle_momentum_Idx]);
        hJetProxyPt_ProtonFromLambda_events->Fill(pt[highest_particle_momentum_Idx]);
        hJetProxyPhi_ProtonFromLambda_events->Fill(jet_phi);

        if (pt[highest_particle_momentum_Idx] >= 2.){
            hNEv_WithLambda_andJets_Pt_leq_2_unit->Fill(0);
        }
        
        // Need to know if this event has any jets in the mid-rapidity region. If not, can't use it!
        // (TODO: If needed, you can relax this demand a little bit, to increase statistics in your testing!)
        Bool_t contains_jet_in_midrapidity = (std::fabs(y[highest_particle_momentum_Idx]) < 0.5); // Using fabs for easyness of reading for older coders.
                                                                                                  // C++11 already overloaded std::abs() to work appropriately with floats
        // Added an additional loop here just to get the dN/d DeltaBetaJ distribution to see 
        // if there is any preferential production of Lambdas wrlt the jet.
        // Do notice that this Beta angle is the angular aperture of the Lambda and the jet,
        // so it is a more general angle than Phi or Theta and is valid vor all rapidities!
        // (in the |y| < 2.5 though, so that we actually get jets, not beam remnants)
        for (int i = 0; i < lambda_4momenta.size(); i++){
            TLorentzVector lambda_4momentum = lambda_4momenta[i];
            TVector3 lambda_trivec = lambda_4momentum.Vect();
            Double_t DeltaBetaJ = jet_trimomentum.Angle(lambda_trivec);
            hLambdaDeltaBetaJ_WithProtonFromLambda->Fill(DeltaBetaJ); // Not necessarily a useful event in the sense of |y|<0.5, but still a useful statistic to know
        }
        if (!contains_jet_in_midrapidity){continue;}
        // std::cout << "\tThis event contains a jet in midrapidity!" << std::endl;

        // N_useful_events += 1;
        hNUsefulEv->Fill(0);

        hJetProxyY_UsefulEvent->Fill(y[highest_particle_momentum_Idx]);
        hJetProxyPt_UsefulEvent->Fill(pt[highest_particle_momentum_Idx]);
        hJetProxyPhi_UsefulEvent->Fill(jet_phi);
        
            // Saving this useful jet for the Ring Observable calculation:
        jet_trimomenta_array.push_back(jet_trimomentum);
        hHighestMomentum_WithLambda_Midrap->Fill(highest_particle_momentum);
        hHighestMomentumPzOverNorm_WithLambda_Midrap->Fill(jet_trimomentum_unit.Z()); // How much of the unit vector is in z. Only for |y| < 0.5 jets

        // std::cout << "\t\tThis jet's momentum is " << highest_particle_momentum << std::endl;

        // Now summing the contribution to the polarization after we confirmed the event is usable:
        for (int i = 0; i < proton_4momenta.size(); i++){ // Goes over all protons and their trimomenta
            TLorentzVector proton_4momentum = proton_4momenta[i];
            TLorentzVector lambda_4momentum = lambda_4momenta[i];
            
            // Rotating the whole coordinate system so that the jet's phi coordinate is zero:
            // (this allows easier combination across different events)
            // (Notice that this has to be done in the laboratory frame! Boosts and rotations do not commute in general!)
            double rotated_proton_phi = wrapToInterval(proton_4momentum.Phi() - jet_phi, 0, 2*PI); // With this, you can go back to using Phi from 0 to 2PI, as opposed to the -PI to PI from Pythia
            proton_4momentum.SetPhi(rotated_proton_phi);
            double lambda_phi = lambda_phi_array[i];
            double rotated_lambda_phi = wrapToInterval(lambda_phi - jet_phi, 0, 2*PI);
            lambda_4momentum.SetPhi(rotated_lambda_phi);
            
            TVector3 beta_inverse = -lambda_4momentum.BoostVector(); // Boost trivector that goes from laboratory frame to the rest frame
            TLorentzVector proton_4momentum_star = proton_4momentum;
            proton_4momentum_star.Boost(beta_inverse); // Takes proton to Lambda rest frame.

            // Finally, using the boosted (and already rotated) vectors:
            TVector3 proton_star_unit_vector = (proton_4momentum_star.Vect()).Unit();
            double X_dot = proton_star_unit_vector.X();
            double Y_dot = proton_star_unit_vector.Y();
            double Z_dot = proton_star_unit_vector.Z();

            double X_dotSquared = X_dot*X_dot;
            double Y_dotSquared = Y_dot*Y_dot;
            double Z_dotSquared = Z_dot*Z_dot;

            // double delta_phi_J = lambda_phi - jet_phi;
            double delta_phi_J = rotated_lambda_phi; // Actually, jet_phi = 0 in this formulation already! So everything is a delta_phi_J angle!
            delta_phi_J = wrapToInterval(delta_phi_J, phi_min_Ring, phi_max_Ring); // Actually the same as rotated_lambda_phi, but in the -PI to PI range instead of 0 to 2PI.

            // Incrementing counters (have to increment them AFTER we selected the event based on the presence of a jet in |y|<0.5):
            // double lambda_pT = std::sqrt(lambda_trimomentum.X()*lambda_trimomentum.X() + lambda_trimomentum.Y()*lambda_trimomentum.Y());
            // double lambda_pT = lambda_trimomentum.Perp(); // Shorthand
            Float_t lambda_y = lambda_y_array[i];
            Float_t lambda_pT = lambda_pT_array[i];
            
            hLambdaCounter_pT_y_DeltaPhiJ->Fill(lambda_pT, lambda_y, rotated_lambda_phi);
            hLambdaCounter_DeltaPhiJ_pT->Fill(rotated_lambda_phi, lambda_pT);
            hLambdaCounter_DeltaPhiJ->Fill(rotated_lambda_phi);
            hLambdaCounter_DeltaPhiJRingAngles->Fill(delta_phi_J); // Should use the wrapped version, from ~-PI to +PI
            hProtonStarCounter_DeltaPhiJRing->Fill(wrapToInterval(proton_4momentum_star.Phi(), phi_min_Ring, phi_max_Ring)); // A debug counter to see if there is any anistropy in the decays (on the Lambda rest frame!)
                // Incrementing X_dot histograms:
            hLambdaAvgDotX_pT_y_DeltaPhiJ->Fill(lambda_pT, lambda_y, rotated_lambda_phi, X_dot);
            hLambdaAvgDotY_pT_y_DeltaPhiJ->Fill(lambda_pT, lambda_y, rotated_lambda_phi, X_dot);
            hLambdaAvgDotZ_pT_y_DeltaPhiJ->Fill(lambda_pT, lambda_y, rotated_lambda_phi, X_dot);

            hLambdaAvgDotXSquared_pT_y_DeltaPhiJ->Fill(lambda_pT, lambda_y, rotated_lambda_phi, X_dotSquared);
            hLambdaAvgDotYSquared_pT_y_DeltaPhiJ->Fill(lambda_pT, lambda_y, rotated_lambda_phi, X_dotSquared);
            hLambdaAvgDotZSquared_pT_y_DeltaPhiJ->Fill(lambda_pT, lambda_y, rotated_lambda_phi, X_dotSquared);

                // Calculating the ring observable in the non-averaged \hat{p}*_D:
                // (Have to use jet_trimomentum instead of the jet_px, jet_py, jet_pz components due to the rotation!)
            double cross_x = jet_trimomentum.Y()*lambda_4momentum.Z() - jet_trimomentum.Z()*lambda_4momentum.Y();
            double cross_y = jet_trimomentum.Z()*lambda_4momentum.X() - jet_trimomentum.X()*lambda_4momentum.Z();
            double cross_z = jet_trimomentum.X()*lambda_4momentum.Y() - jet_trimomentum.Y()*lambda_4momentum.X();
            double cross_product_norm = std::sqrt(cross_x*cross_x + cross_y*cross_y + cross_z*cross_z);
            double no_avg_ring_observable = 3./(alpha_H) * ((proton_star_unit_vector.X()*cross_x + proton_star_unit_vector.Y()*cross_y
                                                            + proton_star_unit_vector.Z()*cross_z)/cross_product_norm);
            hRingObservable_daughter_no_avg->Fill(delta_phi_J, no_avg_ring_observable);

            if (lambda_pT > 0.5 && lambda_pT < 1.5){
                hRingObservable_daughter_no_avg_PtCuts->Fill(no_avg_ring_observable);
                hRingObservable_daughter_no_avg_PtCuts_integrated->Fill(no_avg_ring_observable);

                hLambdaCounter_DeltaPhiJ_PtCuts->Fill(rotated_lambda_phi);
                hLambdaCounter_DeltaPhiJRingAngles_PtCuts->Fill(delta_phi_J);
                hProtonStarCounter_DeltaPhiJRing_PtCuts->Fill(delta_phi_J); // Should use delta_phi_J because the angle goes from ~-PI to PI here
            }
        }
    }
    std::cout << "Finished the event loop!" << std::endl;

    ////////////////////////////////////////////////////////////
    //// Finishing the averaging process for the histograms ////
    ////////////////////////////////////////////////////////////
    std::cout << "\nNow averaging histograms" << std::endl;
    // (And notice that with Pythia's data these can be simple averages, without any weights!)
    hLambdaAvgDotX_pT_y_DeltaPhiJ->Divide(hLambdaCounter_pT_y_DeltaPhiJ);
    hLambdaAvgDotY_pT_y_DeltaPhiJ->Divide(hLambdaCounter_pT_y_DeltaPhiJ);
    hLambdaAvgDotZ_pT_y_DeltaPhiJ->Divide(hLambdaCounter_pT_y_DeltaPhiJ);

    // Averaging the second momentum of this unit vector distribution too:
    hLambdaAvgDotXSquared_pT_y_DeltaPhiJ->Divide(hLambdaCounter_pT_y_DeltaPhiJ);
    hLambdaAvgDotYSquared_pT_y_DeltaPhiJ->Divide(hLambdaCounter_pT_y_DeltaPhiJ);
    hLambdaAvgDotZSquared_pT_y_DeltaPhiJ->Divide(hLambdaCounter_pT_y_DeltaPhiJ);

    // Cloning everything by hand (usually behaves better than ROOT's Clone()):
    CopyTH3D(hLambdaAvgDotX_pT_y_DeltaPhiJ, hLambdaPolStarX_pT_y_DeltaPhiJReco);
    CopyTH3D(hLambdaAvgDotY_pT_y_DeltaPhiJ, hLambdaPolStarY_pT_y_DeltaPhiJReco);
    CopyTH3D(hLambdaAvgDotZ_pT_y_DeltaPhiJ, hLambdaPolStarZ_pT_y_DeltaPhiJReco);
    hLambdaPolStarX_pT_y_DeltaPhiJReco->SetTitle("hLambdaPolStarX_pT_y_DeltaPhiJReco");
    hLambdaPolStarY_pT_y_DeltaPhiJReco->SetTitle("hLambdaPolStarY_pT_y_DeltaPhiJReco");
    hLambdaPolStarZ_pT_y_DeltaPhiJReco->SetTitle("hLambdaPolStarZ_pT_y_DeltaPhiJReco");
    // Applying the transform that turns them into polarizations -- In the lambda rest frame!
    hLambdaPolStarX_pT_y_DeltaPhiJReco->Scale(3.0 / alpha_H);
    hLambdaPolStarY_pT_y_DeltaPhiJReco->Scale(3.0 / alpha_H);
    hLambdaPolStarZ_pT_y_DeltaPhiJReco->Scale(3.0 / alpha_H);

    ///////////////////////////////////////////////////////////
    //// Calculating the error bars for an UNWEIGHTED case ////
    ///////////////////////////////////////////////////////////
    // (using the same approach as in lambda_pol_toy_model.cxx, or section 1 of the "Error_propagation_discussion_of_Ring.pdf" document)
    CopyTH3D(hLambdaAvgDotX_pT_y_DeltaPhiJ, hLambdaAvgDotX_pT_y_DeltaPhiJ_Err);
    CopyTH3D(hLambdaAvgDotY_pT_y_DeltaPhiJ, hLambdaAvgDotY_pT_y_DeltaPhiJ_Err);
    CopyTH3D(hLambdaAvgDotZ_pT_y_DeltaPhiJ, hLambdaAvgDotZ_pT_y_DeltaPhiJ_Err);
    // Now taking the square to get \mu^2:
    hLambdaAvgDotX_pT_y_DeltaPhiJ_Err->Multiply(hLambdaAvgDotX_pT_y_DeltaPhiJ_Err);
    hLambdaAvgDotY_pT_y_DeltaPhiJ_Err->Multiply(hLambdaAvgDotY_pT_y_DeltaPhiJ_Err);
    hLambdaAvgDotZ_pT_y_DeltaPhiJ_Err->Multiply(hLambdaAvgDotZ_pT_y_DeltaPhiJ_Err);

    // Now replacing the contents with <p_i^2> - \mu^2
    hLambdaAvgDotX_pT_y_DeltaPhiJ_Err->Add(hLambdaAvgDotXSquared_pT_y_DeltaPhiJ, hLambdaAvgDotX_pT_y_DeltaPhiJ_Err, 1., -1.);
    hLambdaAvgDotY_pT_y_DeltaPhiJ_Err->Add(hLambdaAvgDotYSquared_pT_y_DeltaPhiJ, hLambdaAvgDotY_pT_y_DeltaPhiJ_Err, 1., -1.);
    hLambdaAvgDotZ_pT_y_DeltaPhiJ_Err->Add(hLambdaAvgDotZSquared_pT_y_DeltaPhiJ, hLambdaAvgDotZ_pT_y_DeltaPhiJ_Err, 1., -1.);

    // Finally, calculating Var(<p_D>) = sqrt(1/(N_particles-1) * (<p_i^2> - \mu^2)) for each component of <\hat{p}_D>:
        // (This is different from the original lambda_pol_toy_model.cxx case. There, we only had N_resamples things that
        // behaved like particles, i.e., N_resamples protons with their own angles and that had been sampled from the
        // angular distribution of the decay. Here, though, each particle is an individual object, not an average as in iSS!)
            // First, calculating the N_particles-1 part, only for the bins that have a non-zero number of particles:
    CopyTH3D(hLambdaCounter_pT_y_DeltaPhiJ, hLambdaCounter_minus1_pT_y_DeltaPhiJ); // Necessary for an unbiased error estimator for low-statistics
    AddScalarIfLargerThanOne(hLambdaCounter_minus1_pT_y_DeltaPhiJ, -1); // Different from AddScalar: it only subtracts 1 if that bin is >1 (I can't go creating negative multiplicities in a bin, nor making it zero!)
    hLambdaAvgDotX_pT_y_DeltaPhiJ_Err->Divide(hLambdaCounter_minus1_pT_y_DeltaPhiJ); // Notice you have to scale BEFORE the square root takes place
    hLambdaAvgDotY_pT_y_DeltaPhiJ_Err->Divide(hLambdaCounter_minus1_pT_y_DeltaPhiJ);
    hLambdaAvgDotZ_pT_y_DeltaPhiJ_Err->Divide(hLambdaCounter_minus1_pT_y_DeltaPhiJ);
            // Taking the sqrt for each of them with an in-place method:
    SqrtHist(hLambdaAvgDotX_pT_y_DeltaPhiJ_Err);
    SqrtHist(hLambdaAvgDotY_pT_y_DeltaPhiJ_Err);
    SqrtHist(hLambdaAvgDotZ_pT_y_DeltaPhiJ_Err);
            // Scaling these errors to get the error on polarization:
    CopyTH3D(hLambdaAvgDotX_pT_y_DeltaPhiJ_Err, hLambdaPolStarX_pT_y_DeltaPhiJReco_Err);
    CopyTH3D(hLambdaAvgDotY_pT_y_DeltaPhiJ_Err, hLambdaPolStarY_pT_y_DeltaPhiJReco_Err);
    CopyTH3D(hLambdaAvgDotZ_pT_y_DeltaPhiJ_Err, hLambdaPolStarZ_pT_y_DeltaPhiJReco_Err);
    hLambdaPolStarX_pT_y_DeltaPhiJReco_Err->Scale(3.0 / alpha_H);
    hLambdaPolStarY_pT_y_DeltaPhiJReco_Err->Scale(3.0 / alpha_H);
    hLambdaPolStarZ_pT_y_DeltaPhiJReco_Err->Scale(3.0 / alpha_H);
    ///////////////////////////////////////////////////////////////
    std::cout << "Starting the ring observable calculation from polarization's loop" << std::endl;
    // Calculating the Ring Observable for each of the binned entries:
        // Will also calculate a boosted polarization to compare with the plots in Vitor's paper
        // (even though these won't be used for the Ring Observable)
        // This part is essentially the same implementation as in lambda_pol_toy_model.cxx.
        // The data is already binned and the uncertainties are essentially the same too.
    for (int pT_idx = 1; pT_idx <= N_bins_pT; pT_idx++){ // bin index starts at 1 (not 0)
        for (int y_idx = 1; y_idx <= N_bins_rap; y_idx++){
            for (int phi_idx = 1; phi_idx <= N_bins_phi; phi_idx++){
                // Getting all the needed values for all three spatial components:
                double pT_bin_center = hLambdaCounter_pT_y_DeltaPhiJ->GetXaxis()->GetBinCenter(pT_idx);
                double y_bin_center = hLambdaCounter_pT_y_DeltaPhiJ->GetYaxis()->GetBinCenter(y_idx);
                double phi_bin_center = hLambdaCounter_pT_y_DeltaPhiJ->GetZaxis()->GetBinCenter(phi_idx);

                // std::cout << "pT bin " << pT_idx << ", center (" << pT_bin_center << "), y: " << y_bin_center << " phi: " << phi_bin_center << std::endl;

                double PolX_star_pTy_DeltaPhiJ = hLambdaPolStarX_pT_y_DeltaPhiJReco->GetBinContent(pT_idx, y_idx, phi_idx);
                double PolY_star_pTy_DeltaPhiJ = hLambdaPolStarY_pT_y_DeltaPhiJReco->GetBinContent(pT_idx, y_idx, phi_idx);
                double PolZ_star_pTy_DeltaPhiJ = hLambdaPolStarZ_pT_y_DeltaPhiJReco->GetBinContent(pT_idx, y_idx, phi_idx);

                double PolX_star_pTy_DeltaPhiJ_Err = hLambdaPolStarX_pT_y_DeltaPhiJReco_Err->GetBinContent(pT_idx, y_idx, phi_idx); // Or could just use a GetBinError() from hLambdaPolStarX_pT_y_DeltaPhiJReco
                double PolY_star_pTy_DeltaPhiJ_Err = hLambdaPolStarY_pT_y_DeltaPhiJReco_Err->GetBinContent(pT_idx, y_idx, phi_idx);
                double PolZ_star_pTy_DeltaPhiJ_Err = hLambdaPolStarZ_pT_y_DeltaPhiJReco_Err->GetBinContent(pT_idx, y_idx, phi_idx);
                // std::cout << "DEBUG! PolX_star_pTy_DeltaPhiJ_Err: " << PolX_star_pTy_DeltaPhiJ_Err << std::endl;

                double current_DeltaPhiJ_bin_multiplicity = hLambdaCounter_DeltaPhiJ->GetBinContent(phi_idx); // A value integrated over the y and pT axes, also matching the phi definition (0 to 2*PI) of this particular phi_idx.

                // Now calculating a mean 4-momentum of the Lambda in each bin:
                double mean_px = pT_bin_center * std::cos(phi_bin_center);
                double mean_py = pT_bin_center * std::sin(phi_bin_center);
                double mean_mT = std::sqrt(mL2 + pT_bin_center * pT_bin_center);
                double mean_pz = mean_mT * std::sinh(y_bin_center);
                double mean_E = std::sqrt(mL2 + mean_px * mean_px + mean_py * mean_py + mean_pz * mean_pz);

                TLorentzVector mean_Lambda_4vec_lab(mean_px, mean_py, mean_pz, mean_E);
                    // Building boost vector for the polarization in this bin:
                TVector3 beta = mean_Lambda_4vec_lab.BoostVector(); // The boost that takes from the Lambda rest frame into the lab frame.

                // Building the mean polarization 4-momentum, with the fact that S^0 = 0 in the particle's own rest frame (no temporal component of spin in the rest frame!):
                TLorentzVector P_Lambda_star_mean_4vec(PolX_star_pTy_DeltaPhiJ, PolY_star_pTy_DeltaPhiJ, PolZ_star_pTy_DeltaPhiJ, 0);
                TLorentzVector P_Lambda_lab_mean_4vec = P_Lambda_star_mean_4vec;
                P_Lambda_lab_mean_4vec.Boost(beta);

                // For the boosted polarization:
                hLambdaPolX_pT_y_DeltaPhiJReco->Fill(pT_bin_center, y_bin_center, phi_bin_center, P_Lambda_lab_mean_4vec.X());
                hLambdaPolY_pT_y_DeltaPhiJReco->Fill(pT_bin_center, y_bin_center, phi_bin_center, P_Lambda_lab_mean_4vec.Y());
                hLambdaPolZ_pT_y_DeltaPhiJReco->Fill(pT_bin_center, y_bin_center, phi_bin_center, P_Lambda_lab_mean_4vec.Z());

                // Also doing a projection in 2D over pT and y:
                    // Correction -- When reducing dimensionality, you need to ponder the polarization by the multiplicity of the bin it belongs to!
                    // (summing over all bins as if they were never there requires you to sum over particles and their appropriate weights):
                double current_pTyphi_bin_multiplicity = hLambdaCounter_pT_y_DeltaPhiJ->GetBinContent(pT_idx, y_idx, phi_idx);
                // hLambdaPolX_pT_yReco->Fill(pT_bin_center, y_bin_center, P_Lambda_lab_mean_4vec.X() * current_pTyphi_bin_multiplicity); // Don't need to average over N_phi bins or anything like that!
                // hLambdaPolY_pT_yReco->Fill(pT_bin_center, y_bin_center, P_Lambda_lab_mean_4vec.Y() * current_pTyphi_bin_multiplicity);
                // hLambdaPolZ_pT_yReco->Fill(pT_bin_center, y_bin_center, P_Lambda_lab_mean_4vec.Z() * current_pTyphi_bin_multiplicity);
                    // Don't actually need these pT_y projections here!
                
                    // A second 2D projection, in the same way it was plotted on Vitor's paper:
                hLambdaPolX_DeltaPhiJ_pTReco->Fill(phi_bin_center, pT_bin_center, P_Lambda_lab_mean_4vec.X() * current_pTyphi_bin_multiplicity); // Also with the multiplicity correction
                hLambdaPolY_DeltaPhiJ_pTReco->Fill(phi_bin_center, pT_bin_center, P_Lambda_lab_mean_4vec.Y() * current_pTyphi_bin_multiplicity);
                hLambdaPolZ_DeltaPhiJ_pTReco->Fill(phi_bin_center, pT_bin_center, P_Lambda_lab_mean_4vec.Z() * current_pTyphi_bin_multiplicity);

                    // Finally, a 1D projection for the distribution in Phi only:
                hLambdaPolX_DeltaPhiJReco->Fill(phi_bin_center, P_Lambda_lab_mean_4vec.X() * current_pTyphi_bin_multiplicity); // Also with the multiplicity correction
                hLambdaPolY_DeltaPhiJReco->Fill(phi_bin_center, P_Lambda_lab_mean_4vec.Y() * current_pTyphi_bin_multiplicity);
                hLambdaPolZ_DeltaPhiJReco->Fill(phi_bin_center, P_Lambda_lab_mean_4vec.Z() * current_pTyphi_bin_multiplicity);

                // Calculating the ring observable:
                    // Cross product (trigger X p)
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    // DISCLAIMER: This below is an APPROXIMATION! The coordinate system is defined such that jet_py = 0, but jet_pz may still be non-zero!
                    // I am assuming, though, that this Z component will be very small in comparison to the other components, or equivalently that we average
                    // over so many possible jets that the average \hat{t} jet vector is only facing X in this coordinate system (equivalently, we assume that
                    // the average jet vector lies in the XY plane, with |y| = 0).
                    // This is the easiest way of joining different events as if they had physically equivalent jets (jets flat in XY), and whose jets can be
                    // described by a single angular variable.
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                const double x_trigger = 1.;
                const double y_trigger = 0.;
                const double z_trigger = 0.;
                double cross_x = y_trigger*mean_pz - z_trigger*mean_py;
                double cross_y = z_trigger*mean_px - x_trigger*mean_pz;
                double cross_z = x_trigger*mean_py - y_trigger*mean_px;

                double cross_product_norm = std::sqrt(cross_x*cross_x + cross_y*cross_y + cross_z*cross_z);

                double delta_phi_J = phi_bin_center; // In this coordinate definition, the transverse part of the jet is always facing \hat{x}, and phi is already a delta_phi_J
                delta_phi_J = wrapToInterval(delta_phi_J, phi_min_Ring, phi_max_Ring);
                int phiRing_idx = hLambdaCounter_DeltaPhiJRingAngles->FindBin(delta_phi_J);
                double current_DeltaPhiJRing_bin_multiplicity_PtCuts = hLambdaCounter_DeltaPhiJRingAngles_PtCuts->GetBinContent(phiRing_idx);

                double RP_star = (PolX_star_pTy_DeltaPhiJ*cross_x + PolY_star_pTy_DeltaPhiJ*cross_y + PolZ_star_pTy_DeltaPhiJ*cross_z)/cross_product_norm;
                hRingObservable_PolStar->Fill(delta_phi_J, RP_star * current_pTyphi_bin_multiplicity);

                // Calculating the current bin's error (following "Error_propagation_discussion_of_Ring.pdf")
                    // -- Say the ring observable is a (weighted) average R = <P_i * A_i>_i where i is each bin, so that sigma(R) = sqrt(sum_i A_i^2 * sigma(P_i)^2 * w_i/W^2)
                    // w_i is the current 3D bin's particle multiplicity, W is the current_DeltaPhiJRing_bin_multiplicity, P_i the polarization and A_i the rest of the ring observable.
                    // Actually, notice that you have to combine the three spatial components here via the cross_x, cross_y, cross_z products too!
                    // So instead of sigma(P_i), you will have sigma(P_i_k) * cross_k for each component, and then you just sum the 3 squares inside each i-th bin
                double err_x = PolX_star_pTy_DeltaPhiJ_Err * cross_x;
                double err_y = PolY_star_pTy_DeltaPhiJ_Err * cross_y;
                double err_z = PolZ_star_pTy_DeltaPhiJ_Err * cross_z;
                double current_RingObservableReco_Err_squared = (err_x*err_x + err_y*err_y + err_z*err_z)/(cross_product_norm*cross_product_norm) // Summing the errors of the three spatial components already weighted by their cross product component in err_i
                                                                 * current_pTyphi_bin_multiplicity/current_DeltaPhiJ_bin_multiplicity; // w_i --> Careful! This is not w_i^2! We are really just summing the square errors as if independent!
                                                                // Just for testing, placed this normalization outside the loop:
                                                                //  /(current_DeltaPhiJ_bin_multiplicity*current_DeltaPhiJ_bin_multiplicity); // 1/W^2
                hRingObservable_PolStar_Err->Fill(delta_phi_J, current_RingObservableReco_Err_squared); // This will need to receive a sqrt() later, after the averaging!

                if (pT_bin_center > 0.5 && pT_bin_center < 1.5){
                    hRingObservable_PolStar_PtCuts->Fill(delta_phi_J, RP_star * current_pTyphi_bin_multiplicity);
                    hRingObservable_PolStar_PtCuts_integrated->Fill(delta_phi_J, RP_star * current_pTyphi_bin_multiplicity);

                    // Now for the PtCuts variation of error:
                    hRingObservable_PolStar_PtCuts_Err->Fill(delta_phi_J, current_RingObservableReco_Err_squared); // This will need to receive a sqrt() later, after the averaging!
                    hRingObservable_PolStar_PtCuts_integrated_Err->Fill(delta_phi_J, current_RingObservableReco_Err_squared);
                }
            }
        }
    }
    std::cout << "Loop finished. Starting error histogram operations" << std::endl;
    // Taking the square root of the summed current_RingObservableReco_Err_squared components to properly get the error, and then setting that as the error of the ring TH1D's
    SqrtHist(hRingObservable_PolStar_Err);
    SqrtHist(hRingObservable_PolStar_PtCuts_Err);
    SqrtHist(hRingObservable_PolStar_PtCuts_integrated_Err);

    // Reintroducing the three removed normalizations (AFTER the square root is taken!):
    hRingObservable_PolStar_Err->Divide(hLambdaCounter_DeltaPhiJRingAngles);
    hRingObservable_PolStar_PtCuts_Err->Divide(hLambdaCounter_DeltaPhiJRingAngles_PtCuts);
    hRingObservable_PolStar_PtCuts_integrated_Err->Scale(1./hLambdaCounter_DeltaPhiJRingAngles_PtCuts->Integral());

    // The previous method of error propagation of the integral was not working. Will fix it by hand:
    double integrated_err = 0;
    for (int i = 1; i <= hRingObservable_PolStar_PtCuts_Err->GetNbinsX(); i++){
        double err_i = hRingObservable_PolStar_PtCuts_Err->GetBinContent(i);
        integrated_err += err_i*err_i;
    }
    integrated_err = std::sqrt(integrated_err);
    hRingObservable_PolStar_PtCuts_integrated_Err->SetBinContent(1, integrated_err);

    // Normalizing the polarization histograms that correspond to lower dimensional projections of the TH3D's of pT,y,phi polarization:
    hLambdaPolX_DeltaPhiJ_pTReco->Divide(hLambdaCounter_DeltaPhiJ_pT);
    hLambdaPolY_DeltaPhiJ_pTReco->Divide(hLambdaCounter_DeltaPhiJ_pT);
    hLambdaPolZ_DeltaPhiJ_pTReco->Divide(hLambdaCounter_DeltaPhiJ_pT);

    hLambdaPolX_DeltaPhiJReco->Divide(hLambdaCounter_DeltaPhiJ);
    hLambdaPolY_DeltaPhiJReco->Divide(hLambdaCounter_DeltaPhiJ);
    hLambdaPolZ_DeltaPhiJReco->Divide(hLambdaCounter_DeltaPhiJ);

    // Normalizing all ring observables by the particle weights to complete the averaging:
        // (Same normalization as done with the errors, but without the SqrtHist part)
    hRingObservable_PolStar->Divide(hLambdaCounter_DeltaPhiJRingAngles);
    hRingObservable_PolStar_PtCuts->Divide(hLambdaCounter_DeltaPhiJRingAngles_PtCuts);
    hRingObservable_PolStar_PtCuts_integrated->Scale(1./hLambdaCounter_DeltaPhiJRingAngles_PtCuts->Integral());

    SetErrorsFromErrHist(hRingObservable_PolStar, hRingObservable_PolStar_Err);
    SetErrorsFromErrHist(hRingObservable_PolStar_PtCuts, hRingObservable_PolStar_PtCuts_Err);
    SetErrorsFromErrHist(hRingObservable_PolStar_PtCuts_integrated, hRingObservable_PolStar_PtCuts_integrated_Err);

    /////////////////////////////////////
    //// Exporting into a .root file ////
    /////////////////////////////////////
    std::cout << "\nExporting into a .root" << std::endl;
    std::string output_file_path = output_folder + input_card_name;
    std::vector<std::string> centrality_limits = {"0-5", "5-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90"};
    output_file_path += (std::string) "_" + centrality_limits[target_centrality_bin-1] + "_only_" + std::to_string(N_events) + "evInAllClasses.root"; // Split in two to make the "string + const char*" operation work
    TFile f(output_file_path.c_str(), "RECREATE");

    // Overall event statistics:
    hNEv_InCentralityClass->Write();
    hNEv_WithLambda->Write();
    hNEv_proton_from_Lambda->Write();
    hNEv_WithLambda_andJets_Pt_leq_2_unit->Write();
    hNUsefulEv->Write();

    hLambdaPt->Write();
    hLambdaY->Write();
    hLambdaPhi->Write();
    hLambdaDeltaBetaJ_WithProtonFromLambda->Write();
    hLambdaDeltaBetaJ_AllLambdaDecayModes->Write();

    hLambdaPt_proton_mothers->Write();
    hLambdaY_proton_mothers->Write();
    hLambdaPhi_proton_mothers->Write();

    hHighestMomentum_all_jets->Write();
    hHighestMomentumPzOverNorm_all_jets->Write();
    hHighestMomentum_EvWithLambda->Write();
    hHighestMomentumPzOverNorm_EvWithLambda->Write();
    hHighestMomentum_WithLambda_Midrap->Write();
    hHighestMomentumPzOverNorm_WithLambda_Midrap->Write();

    // Jet kinematics:
    hJetProxyY_all_jets->Write();
    hJetProxyPt_all_jets->Write();
    hJetProxyPhi_all_jets->Write();

    hJetProxyY_ProtonFromLambda_events->Write();
    hJetProxyPt_ProtonFromLambda_events->Write();
    hJetProxyPhi_ProtonFromLambda_events->Write();

    hJetProxyY_UsefulEvent->Write();
    hJetProxyPt_UsefulEvent->Write();
    hJetProxyPhi_UsefulEvent->Write();

    hProtonStarCounter_DeltaPhiJRing_PtCuts->Write();

        // Polarization measurements:
    hLambdaCounter_pT_y_DeltaPhiJ->Write();
    hLambdaAvgDotX_pT_y_DeltaPhiJ->Write();
    hLambdaAvgDotY_pT_y_DeltaPhiJ->Write();
    hLambdaAvgDotZ_pT_y_DeltaPhiJ->Write();

    hLambdaAvgDotXSquared_pT_y_DeltaPhiJ->Write();
    hLambdaAvgDotYSquared_pT_y_DeltaPhiJ->Write();
    hLambdaAvgDotZSquared_pT_y_DeltaPhiJ->Write();
    hLambdaPolX_pT_y_DeltaPhiJReco->Write();
    hLambdaPolY_pT_y_DeltaPhiJReco->Write();
    hLambdaPolZ_pT_y_DeltaPhiJReco->Write();

    hLambdaCounter_DeltaPhiJ_pT->Write();
    hLambdaPolX_DeltaPhiJ_pTReco->Write();
    hLambdaPolY_DeltaPhiJ_pTReco->Write();
    hLambdaPolZ_DeltaPhiJ_pTReco->Write();

    hLambdaCounter_DeltaPhiJ->Write();
    hLambdaCounter_DeltaPhiJRingAngles->Write();
    hProtonStarCounter_DeltaPhiJRing->Write();

    hLambdaCounter_DeltaPhiJ_PtCuts->Write();
    hLambdaCounter_DeltaPhiJRingAngles_PtCuts->Write();

    hLambdaPolX_DeltaPhiJReco->Write();
    hLambdaPolY_DeltaPhiJReco->Write();
    hLambdaPolZ_DeltaPhiJReco->Write();
    hRingObservable_PolStar->Write();
    hRingObservable_PolStar_PtCuts->Write();
    hRingObservable_PolStar_PtCuts_integrated->Write();

    hRingObservable_daughter_no_avg->Write();
    hRingObservable_daughter_no_avg_PtCuts->Write();
    hRingObservable_daughter_no_avg_PtCuts_integrated->Write();
    // TODO: Implement error bars for these no_avg estimates!

    hLambdaAvgDotX_pT_y_DeltaPhiJ_Err->Write();
    hLambdaAvgDotY_pT_y_DeltaPhiJ_Err->Write();
    hLambdaAvgDotZ_pT_y_DeltaPhiJ_Err->Write();
    hLambdaCounter_minus1_pT_y_DeltaPhiJ->Write();

    hLambdaPolStarX_pT_y_DeltaPhiJReco_Err->Write();
    hLambdaPolStarY_pT_y_DeltaPhiJReco_Err->Write();
    hLambdaPolStarZ_pT_y_DeltaPhiJReco_Err->Write();
    hRingObservable_PolStar_Err->Write();
    hRingObservable_PolStar_PtCuts_Err->Write();
    hRingObservable_PolStar_PtCuts_integrated_Err->Write();

    f.Close();
    std::cout << "Done exporting!" << std::endl;

    std::cout << "\n\nCode execution finished with code 0. Thank you!" << std::endl;
    return 0;
}









///////////////////////////////////
///////// Helper functions ////////
///////////////////////////////////
// Wrap an angle into [phi_min, phi_max) for the coordinate shift that redefines the phi angle as having phi = 0 pointing towards the jet direction.
// This is useful to define a \hat{t} vector that is always pointing towards (1, 0, 0), making it easier to calculate the ring observable in the
// "with bullet" case:
// Fast wrapping into [phi_min, phi_max)
inline double wrapToInterval(double phi, double phi_min, double phi_max){
    double range = phi_max - phi_min;
    // Use fmod to reduce directly
    phi = std::fmod(phi - phi_min, range);
    if (phi < 0) phi += range;
    return phi + phi_min;
}

// Copy contents of a TH3D into a new TH3D -- Better than the Clone() function...
void CopyTH3D(const TH3D* source, TH3D* target){
    if (!source || !target) return;
    // std::cout << "Now copying " << source->GetTitle() << " into " << target->GetTitle() << std::endl;

    // Check binning matches
    if (source->GetNbinsX() != target->GetNbinsX() || source->GetNbinsY() != target->GetNbinsY() || source->GetNbinsZ() != target->GetNbinsZ()){
        printf("Error: Histograms have different binning!\n");
        return;
    }

    // Copy bin contents and errors
    for (int ix = 0; ix <= source->GetNbinsX() + 1; ++ix){
        for (int iy = 0; iy <= source->GetNbinsY() + 1; ++iy){
            for (int iz = 0; iz <= source->GetNbinsZ() + 1; ++iz){
                int bin = source->GetBin(ix, iy, iz);
                target->SetBinContent(bin, source->GetBinContent(bin));
                target->SetBinError(bin, source->GetBinError(bin));
            }
        }
    }
}

// ROOT lacks a sqrt function for its histograms, thus implemented it here by hand.
// This function modifies the histogram in-place, so no returns needed.
// An uglier rewrite of this function, that discriminates all three possible dimensions instead of using shorter syntax.
void SqrtHist(TH1* h){
    int dim = h->GetDimension();

    if (dim == 1){
        int nbx = h->GetNbinsX();
        for (int ix = 0; ix <= nbx + 1; ++ix){
            double val = h->GetBinContent(ix);
            if (val >= 0)
                h->SetBinContent(ix, std::sqrt(val));
            else {
                std::cerr << "Warning: negative bin content in bin " << ix
                          << ": " << val << std::endl;
                h->SetBinContent(ix, 0.);
            }
        }
    }
    else if (dim == 2){
        int nbx = h->GetNbinsX();
        int nby = h->GetNbinsY();
        for (int ix = 0; ix <= nbx + 1; ++ix){
            for (int iy = 0; iy <= nby + 1; ++iy){
                int bin = h->GetBin(ix, iy);
                double val = h->GetBinContent(bin);
                if (val >= 0)
                    h->SetBinContent(bin, std::sqrt(val));
                else {
                    std::cerr << "Warning: negative bin content in bin (" << ix
                              << "," << iy << "): " << val << std::endl;
                    h->SetBinContent(bin, 0.);
                }
            }
        }
    }
    else if (dim == 3){
        int nbx = h->GetNbinsX();
        int nby = h->GetNbinsY();
        int nbz = h->GetNbinsZ();
        for (int ix = 0; ix <= nbx + 1; ++ix){
            for (int iy = 0; iy <= nby + 1; ++iy){
                for (int iz = 0; iz <= nbz + 1; ++iz){
                    int bin = h->GetBin(ix, iy, iz);
                    double val = h->GetBinContent(bin);
                    if (val >= 0)
                        h->SetBinContent(bin, std::sqrt(val));
                    else {
                        std::cerr << "Warning: negative bin content in bin ("
                                  << ix << "," << iy << "," << iz << "): " << val << std::endl;
                        h->SetBinContent(bin, 0.);
                    }
                }
            }
        }
    }
    else {
        std::cerr << "Error: unsupported histogram dimension " << dim << std::endl;
    }
}

// Function that subtracts given scalar from all bins in the histogram:
void AddScalar(TH3D* h, double scalar){
  int nbx = h->GetNbinsX();
  int nby = h->GetNbinsY();
  int nbz = h->GetNbinsZ();

  for (int ix = 1; ix <= nbx; ++ix){
    for (int iy = 1; iy <= nby; ++iy){
      for (int iz = 1; iz <= nbz; ++iz){
        double val = h->GetBinContent(ix, iy, iz);
        h->SetBinContent(ix, iy, iz, val + scalar);
      }
    }
  }
}

void AddScalarIfLargerThanOne(TH3D* h, double scalar){
  int nbx = h->GetNbinsX();
  int nby = h->GetNbinsY();
  int nbz = h->GetNbinsZ();

  for (int ix = 1; ix <= nbx; ++ix){
    for (int iy = 1; iy <= nby; ++iy){
      for (int iz = 1; iz <= nbz; ++iz){
        double val = h->GetBinContent(ix, iy, iz);
        h->SetBinContent(ix, iy, iz, (val > 1) ? val + scalar : val); // Should only subtract if the bin count is larger than one
      }
    }
  }
}



// A quick helper to set the errors in one histogram from the values in an error-histogram:
// (safer rewrite with one check for each dimension. Seems I got something wrong in the previous
// loop, even when replacing stuff like iz <= nbz+1 with iz <= (dim>2) ? nbz+1 : 0)
void SetErrorsFromErrHist(TH1* data, const TH1* errors){
    if (data->GetDimension() != errors->GetDimension()){
        std::cerr << "Error: Histograms have different dimensions!" << std::endl;
        return;
    }

    int dim = data->GetDimension();

    if (dim == 1){
        int nbx = data->GetNbinsX();
        for (int ix = 0; ix <= nbx + 1; ++ix){
            double err = errors->GetBinContent(ix);
            data->SetBinError(ix, err);
        }
    } 
    else if (dim == 2){
        int nbx = data->GetNbinsX();
        int nby = data->GetNbinsY();
        for (int ix = 0; ix <= nbx + 1; ++ix){
            for (int iy = 0; iy <= nby + 1; ++iy){
                int bin = data->GetBin(ix, iy);
                double err = errors->GetBinContent(bin);
                data->SetBinError(bin, err);
            }
        }
    }
    else if (dim == 3){
        int nbx = data->GetNbinsX();
        int nby = data->GetNbinsY();
        int nbz = data->GetNbinsZ();
        for (int ix = 0; ix <= nbx + 1; ++ix){
            for (int iy = 0; iy <= nby + 1; ++iy){
                for (int iz = 0; iz <= nbz + 1; ++iz){
                    int bin = data->GetBin(ix, iy, iz);
                    double err = errors->GetBinContent(bin);
                    data->SetBinError(bin, err);
                }
            }
        }
    }
    else {
        std::cerr << "Error: unsupported histogram dimension " << dim << std::endl;
    }
}

double IntegralAbsNoOverflow(const TH1* h){
    double sum = 0.0;
    if (auto h1 = dynamic_cast<const TH1*>(h)){
        for (int ix = 1; ix <= h1->GetNbinsX(); ix++){
            for (int iy = 1; iy <= h1->GetNbinsY(); iy++){
                for (int iz = 1; iz <= h1->GetNbinsZ(); iz++){
                    sum += std::fabs(h->GetBinContent(ix, iy, iz));
                }
            }
        }
    }
    return sum;
}


// This function will find the bin with the closest y value with relation to the targetValue.
// It searches from right to left and stops when the difference between the current bin and
// targetValue starts to rise, instead of lowering. It only works if the difference is monotonically
// rising or falling.
// USED IN: estimating the closest MC centrality class limits given an experimental limits list.
int FindClosestBin(TH1D* hist, double targetValue){
    int nBins = hist->GetNbinsX();
    int closestBin_idx = nBins;  // Starts at the rightmost bin

    // std::cout << "Nbins: " << nBins << std::endl;
    // std::cout << "Desired value: " << targetValue << std::endl;
    // std::cout << "Starting value: " << hist->GetBinContent(closestBin_idx) << std::endl;
    // std::cout << "Some middle value: " << hist->GetBinContent(40) << std::endl;
    // std::cout << "The first three values from the right: " << hist->GetBinContent(nBins-1) << " " << hist->GetBinContent(nBins-2) << " " << hist->GetBinContent(nBins-3) << std::endl;

        // Starting the initial difference with the first bin from the right:
    double previousDiff = std::fabs(hist->GetBinContent(closestBin_idx) - targetValue);

    // std::cout << "Starting difference: " << previousDiff << std::endl;

    for (int i = nBins - 1; i >= 1; i--){
        double currentContent = hist->GetBinContent(i);
        double currentDiff = std::fabs(currentContent - targetValue);

        // std::cout << "Current centrality edge: " << currentContent << std::endl; // Debug tip: This will print even in the case that the loop is broken! Don't worry, it will not actually take the step that increases the difference!
        // std::cout << "Current difference from target: " << currentDiff << std::endl;

        // if (currentDiff < previousDiff || (currentContent - previousContent) < 1e-15){ // There is an even cleaner method!
        if (currentDiff <= previousDiff){
            // The first bins to the left will be very close to zero, and the code must not stop when the difference does not minimize significantly for those cases!
            // Also, if two bins vary too little, then the currentDiff will be almost the same as the previous one, thus you should continue searching for the
            // best bin.
            closestBin_idx = i;
            previousDiff = currentDiff;
        }
        else{
            // std::cout << "Breaking the loop!" << std::endl;
            break;
        } // diff increased: stop search
    }

    return closestBin_idx;
}