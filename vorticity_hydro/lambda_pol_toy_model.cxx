#include "TROOT.h"
#include "TFile.h"
// #include "TTree.h"

#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TH3.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRotation.h"

#include <vector>
#include <array>
#include <cmath>
#include <random>
#include <chrono>
#include <stdexcept> // For the crash in mapGlobalToCompact
// #include <algorithm> // Useful for std::fill() commands for the standard c++ array buffers
// #include "math.h" // Same as cmath!
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TMath.h"

#include <dirent.h>
#include <iostream>
#include <string>
#include <fstream>

// Parallelization capabilities -- Do notice the new style of ROOT headers used for the TThreadedObject class!
#include "ROOT/TThreadedObject.hxx"
#include <omp.h>

typedef std::vector<std::vector<double>> DoubleMatrix; // An alias

// Some constants for the code:
const double PI = TMath::Pi();
    // Constant from Phys. Rev. D 110, 030001 - Published 1 August, 2024, i.e.: https://journals.aps.org/prd/pdf/10.1103/PhysRevD.110.030001
const double alpha_H = 0.747; // Fixed as a constant, just for this toy model!
const double mass_Lambda = 1115.683; // In MeV/c^2 (+/- 0.006)
const double mass_proton = 938.27208816; // (+/- 0.00000029)
const double mass_pi = 139.57039; // (+/- 0.00018), for both pi+ and pi-
const double mL2 = mass_Lambda * mass_Lambda;
const double mp2 = mass_proton * mass_proton;
const double mpi2 = mass_pi * mass_pi;
const double mp_plus_mpi2 = std::pow(mass_proton + mass_pi, 2);
const double mp_minus_mpi2 = std::pow(mass_proton - mass_pi, 2);

// Function prototypes:
    // Getter functions:
void get_lambda(DoubleMatrix &y_matrix, DoubleMatrix &phi_matrix, DoubleMatrix &px_matrix, DoubleMatrix &py_matrix,
                DoubleMatrix &pz_matrix, DoubleMatrix &pT_matrix, DoubleMatrix &E_matrix, DoubleMatrix& mult_matrix,
                DoubleMatrix &St_matrix, DoubleMatrix &Sx_matrix, DoubleMatrix &Sy_matrix, DoubleMatrix &Sz_matrix,
                int N_events, bool with_bullet);
void get_jet(std::vector<double> &n_event, std::vector<double> &n_random, std::vector<double> &phi_random,
             std::vector<double> &momentum_x, std::vector<double> &momentum_y);

    // Physics part functions
std::pair<double, double> sample_P_angle_proton(double P_Lambda_star_mag, std::mt19937 &rng, std::uniform_real_distribution<double> dist_x,
                                                std::uniform_real_distribution<double> dist_y, std::uniform_real_distribution<double> &dist_azimuth);
std::pair<double, double> sample_P_angle_proton_from_cos_xi(double P_Lambda_star_mag, std::mt19937 &rng, std::uniform_real_distribution<double> &dist_unit,
                                                               std::uniform_real_distribution<double> &dist_azimuth);
std::pair<TLorentzVector, TLorentzVector> Lambda_decay(TLorentzVector Lambda_4_momentum_lab, TVector3 P_Lambda_star, double xi_star, double phi_star);
TVector3 boost_polarization_to_rest_frame(TLorentzVector Lambda_4_momentum_lab, TLorentzVector P_Lambda_lab_4vec);

    // Helper functions:
inline double wrapToInterval(double phi, double phi_min, double phi_max);
void CheckBins(const char* name, const int* bins, int N_bins_total, TH1* hist);
inline int mapGlobalToCompact(const TH1& h, int globalBin);
inline int mapCompactToGlobal(const TH1& h, int compactIndex); // The inverse transform
TH3D* CopyTH3D(const TH3D* h);



// This code assumes we are in Jarvis4! Also, the getter functions are based on Pol_Analysis_Random_hist_ebe.C from the HadrEx_Ph repository
// int lambda_pol_toy_model(){
int main(int argc, char *argv[]){ // Changed the code into a compiler-friendly way, that searches for a "main()" function instead of a main with the same name as the .cxx file that ROOT's .x compiler expects
    // Receiving useful values:
    if (argc != 5) {
        std::cerr << "Usage:" << argv[0] << " N_resamples N_threads with_bullet(int) force_DecayDist_mode(0=default, 1=AlongP, 2=PerpP)" << std::endl;
        return 1; // Indicate an error. (This warning code block is chat-gpt made, by the way!)
    }
    int N_resamples = atoi(argv[1]);
    int N_threads = atoi(argv[2]);
    bool with_bullet = atoi(argv[3]);
    int force_DecayDist_mode = atoi(argv[4]); // 0 for the usual, physical, decay. 1 for decays along Pol. 2 for decays perpendicular to Pol

    ////////////////////
    //// 1 - Initializing variables and histograms:
    ////////////////////
    // A counter for the execution time:
    auto start_time = std::chrono::high_resolution_clock::now();

    std::cout << "Initializing variables" << std::endl;
    int N_events = 250; // There are only 250 events in the desired folder (40_50 has 300 for the no_bullet case, for some reason). Possibly oversampled to get statistics.
    // bool with_bullet = false;

    ROOT::EnableThreadSafety(); // THIS IS MANDATORY TO MAKE THREADING WORK!!!
    auto hLambdaCounter = new TH1D("hLambdaCounter", "", 1, -1, 1);

    // Declaring histogram-related variables:
        // Corrected with the values used in /home/vribeiro/vorticity/iSS-pol/src/spin_polarization.h
    const int N_bins_pT = 30;
    const int N_bins_rap = 10;
    const int N_bins_phi = 100;

    const double pT_min = 0.;
    const double pT_max = 3.;
    const double rap_max = 0.5; // Am selecting particles only in this interval for easier comparison. 
                                // The original data is also clipped to this upper limit (see /home/vribeiro/vorticity/iSS-pol/src/spin_polarization.cpp),
                                // so this binning and this interval will properly reconstruct the original information.
                                // This simulation only has |y|<0.5 and pT between 0 and 3 GeV/c.
    double dpT = (pT_max - pT_min)* 1./N_bins_pT;
    double drap = (rap_max + rap_max) * 1./N_bins_rap;
    // Defining phi intervals in the same way as defined in dev/vribeiro/vorticity/analysis/Pol_Analysis/Pol_Analysis_Random_hist_ebe.C:
    double dphi = 2*PI * 1./N_bins_phi; // Bin width
    double phi_min = -dphi/2.; // The lower edge is defined in such a way that phi = 0 is the bin center.
    double phi_max = (2*PI - dphi) + dphi/2.; // The last step of the vector (2pi - bin_width, or the last bin center), summed with the size of the last bin to get the upper edge.
        // ALTERNATIVELY, YOU COULD JUST SUM dphi/2. TO EVERY PHI VALUE ON phi_matrix! It is probably better to do this than to use these weird phi_min--phi_max intervals.
        // todo: check if the angle convention is the same for the jet's phi angle!
    double phi_min_Ring = phi_min - PI;
    double phi_max_Ring = phi_max - PI;


    // For peace of mind, the same plots that appear on the paper, before any kind of reconstruction:
        // (Properly weighted by the multiplicity of each bin!)
    auto hLambdaCounter_phi_pT_Weighted = ROOT::TThreadedObject<TH2D>("hLambdaCounter_phi_pT_Weighted", "hLambdaCounter_phi_pT_Weighted", N_bins_phi, phi_min, phi_max, N_bins_pT, pT_min, pT_max);
    auto hLambdaPolX_phi_pT_Weighted = ROOT::TThreadedObject<TH2D>("hLambdaPolX_phi_pT_Weighted", "hLambdaPolX_phi_pT_Weighted", N_bins_phi, phi_min, phi_max, N_bins_pT, pT_min, pT_max);
    auto hLambdaPolY_phi_pT_Weighted = ROOT::TThreadedObject<TH2D>("hLambdaPolY_phi_pT_Weighted", "hLambdaPolY_phi_pT_Weighted", N_bins_phi, phi_min, phi_max, N_bins_pT, pT_min, pT_max);
    auto hLambdaPolZ_phi_pT_Weighted = ROOT::TThreadedObject<TH2D>("hLambdaPolZ_phi_pT_Weighted", "hLambdaPolZ_phi_pT_Weighted", N_bins_phi, phi_min, phi_max, N_bins_pT, pT_min, pT_max);

        // A third version of these plots, but now for 3D, to then convert into 2D -- Testing if the method of reducing the dimensionality through an averaging works, or if it is the source of the problem!
    auto hLambdaPolX_pT_y_phi_Weighted = ROOT::TThreadedObject<TH3D>("hLambdaPolX_pT_y_phi_Weighted", "hLambdaPolX_pT_y_phi_Weighted", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);
    auto hLambdaPolY_pT_y_phi_Weighted = ROOT::TThreadedObject<TH3D>("hLambdaPolY_pT_y_phi_Weighted", "hLambdaPolY_pT_y_phi_Weighted", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);
    auto hLambdaPolZ_pT_y_phi_Weighted = ROOT::TThreadedObject<TH3D>("hLambdaPolZ_pT_y_phi_Weighted", "hLambdaPolZ_pT_y_phi_Weighted", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);

            // A test to show that a projection from 3D to 2D works:
    TH2D *hLambdaPolX_phi_pT_Weighted_ProjTEST = new TH2D("hLambdaPolX_phi_pT_Weighted_ProjTEST", "hLambdaPolX_phi_pT_Weighted_ProjTEST", N_bins_phi, phi_min, phi_max, N_bins_pT, pT_min, pT_max);
    TH2D *hLambdaPolY_phi_pT_Weighted_ProjTEST = new TH2D("hLambdaPolY_phi_pT_Weighted_ProjTEST", "hLambdaPolY_phi_pT_Weighted_ProjTEST", N_bins_phi, phi_min, phi_max, N_bins_pT, pT_min, pT_max);
    TH2D *hLambdaPolZ_phi_pT_Weighted_ProjTEST = new TH2D("hLambdaPolZ_phi_pT_Weighted_ProjTEST", "hLambdaPolZ_phi_pT_Weighted_ProjTEST", N_bins_phi, phi_min, phi_max, N_bins_pT, pT_min, pT_max);

            // Declaring a counter histogram and some summation-storing histograms to define polarization in each (pT, y) bin:
    auto hLambdaCounter_pT_y_Weighted = ROOT::TThreadedObject<TH2D>("hLambdaCounter_pT_y_Weighted", "hLambdaCounter_pT_y_Weighted", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max); // Properly weighted counts

        // Declaring some histograms to calculate the polarization in the lab frame -- You will need the phi information to build the full 4-vector of the "mean Lambda" in that bin:
            // Weighting the observables by the multiplicity of each Lambda bin (the mult_matrix data):
    auto hLambdaCounter_pT_y_phi_Weighted = ROOT::TThreadedObject<TH3D>("hLambdaCounter_pT_y_phi_Weighted", "hLambdaCounter_pT_y_phi_Weighted", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);
    auto hLambdaAvgDotX_pT_y_phi_Weighted = ROOT::TThreadedObject<TH3D>("hLambdaAvgDotX_pT_y_phi_Weighted", "hLambdaAvgDotX_pT_y_phi_Weighted", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);
    auto hLambdaAvgDotY_pT_y_phi_Weighted = ROOT::TThreadedObject<TH3D>("hLambdaAvgDotY_pT_y_phi_Weighted", "hLambdaAvgDotY_pT_y_phi_Weighted", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);
    auto hLambdaAvgDotZ_pT_y_phi_Weighted = ROOT::TThreadedObject<TH3D>("hLambdaAvgDotZ_pT_y_phi_Weighted", "hLambdaAvgDotZ_pT_y_phi_Weighted", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);

        // A bunch of alternative counters that are useful from a debugging standpoint (other than the actually useful hLambdaCounter_DeltaphiJRing_Weighted normalizer):
    auto hLambdaCounter_phi_Weighted = ROOT::TThreadedObject<TH1D>("hLambdaCounter_phi_Weighted", "hLambdaCounter_phi_Weighted", N_bins_phi, phi_min, phi_max);
    auto hLambdaCounter_phiRingAngles_Weighted = ROOT::TThreadedObject<TH1D>("hLambdaCounter_phiRingAngles_Weighted", "hLambdaCounter_phiRingAngles_Weighted", N_bins_phi, phi_min_Ring, phi_max_Ring);
    auto hLambdaCounter_DeltaphiJRing_Weighted = ROOT::TThreadedObject<TH1D>("hLambdaCounter_DeltaphiJRing_Weighted", "hLambdaCounter_DeltaphiJRing_Weighted", N_bins_phi, phi_min_Ring, phi_max_Ring);
    auto hProtonCounter_phiRing_Weighted = ROOT::TThreadedObject<TH1D>("hProtonCounter_phiRing_Weighted", "hProtonCounter_phiRing_Weighted", N_bins_phi, phi_min_Ring, phi_max_Ring);
    auto hProtonStarCounter_phiRing_Weighted = ROOT::TThreadedObject<TH1D>("hProtonStarCounter_phiRing_Weighted", "hProtonStarCounter_phiRing_Weighted", N_bins_phi, phi_min_Ring, phi_max_Ring);

    auto hLambdaPolX_phi_Weighted = ROOT::TThreadedObject<TH1D>("hLambdaPolX_phi_Weighted", "hLambdaPolX_phi_Weighted", N_bins_phi, phi_min, phi_max);
    auto hLambdaPolY_phi_Weighted = ROOT::TThreadedObject<TH1D>("hLambdaPolY_phi_Weighted", "hLambdaPolY_phi_Weighted", N_bins_phi, phi_min, phi_max);
    auto hLambdaPolZ_phi_Weighted = ROOT::TThreadedObject<TH1D>("hLambdaPolZ_phi_Weighted", "hLambdaPolZ_phi_Weighted", N_bins_phi, phi_min, phi_max);

        // Now for the reconstructions in 3D:
    TH3D *hLambdaPolStarX_pT_y_phiReco_Weighted;
    TH3D *hLambdaPolStarY_pT_y_phiReco_Weighted;
    TH3D *hLambdaPolStarZ_pT_y_phiReco_Weighted;

        // The lab polarization histograms:
    TH3D *hLambdaPolX_pT_y_phiReco_Weighted = new TH3D("hLambdaPolX_pT_y_phiReco_Weighted", "hLambdaPolX_pT_y_phiReco_Weighted", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);
    TH3D *hLambdaPolY_pT_y_phiReco_Weighted = new TH3D("hLambdaPolY_pT_y_phiReco_Weighted", "hLambdaPolY_pT_y_phiReco_Weighted", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);
    TH3D *hLambdaPolZ_pT_y_phiReco_Weighted = new TH3D("hLambdaPolZ_pT_y_phiReco_Weighted", "hLambdaPolZ_pT_y_phiReco_Weighted", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max, N_bins_phi, phi_min, phi_max);

        // Some 2D lab polarization histograms to compare easily:
    TH2D *hLambdaPolX_pT_yReco_Weighted = new TH2D("hLambdaPolX_pT_yReco_Weighted", "hLambdaPolX_pT_yReco_Weighted", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max);
    TH2D *hLambdaPolY_pT_yReco_Weighted = new TH2D("hLambdaPolY_pT_yReco_Weighted", "hLambdaPolY_pT_yReco_Weighted", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max);
    TH2D *hLambdaPolZ_pT_yReco_Weighted = new TH2D("hLambdaPolZ_pT_yReco_Weighted", "hLambdaPolZ_pT_yReco_Weighted", N_bins_pT, pT_min, pT_max, N_bins_rap, -rap_max, rap_max);

            // Other 2D histograms to compare with Vitor's paper:
    TH2D *hLambdaPolX_phi_pTReco_Weighted = new TH2D("hLambdaPolX_phi_pTReco_Weighted", "hLambdaPolX_phi_pTReco_Weighted", N_bins_phi, phi_min, phi_max, N_bins_pT, pT_min, pT_max);
    TH2D *hLambdaPolY_phi_pTReco_Weighted = new TH2D("hLambdaPolY_phi_pTReco_Weighted", "hLambdaPolY_phi_pTReco_Weighted", N_bins_phi, phi_min, phi_max, N_bins_pT, pT_min, pT_max);
    TH2D *hLambdaPolZ_phi_pTReco_Weighted = new TH2D("hLambdaPolZ_phi_pTReco_Weighted", "hLambdaPolZ_phi_pTReco_Weighted", N_bins_phi, phi_min, phi_max, N_bins_pT, pT_min, pT_max);

        // Finally, 1D lab polarization histograms to compare with the plots in the paper:
    TH1D *hLambdaPolX_phiReco_Weighted = new TH1D("hLambdaPolX_phiReco_Weighted", "hLambdaPolX_phiReco_Weighted", N_bins_phi, phi_min, phi_max);
    TH1D *hLambdaPolY_phiReco_Weighted = new TH1D("hLambdaPolY_phiReco_Weighted", "hLambdaPolY_phiReco_Weighted", N_bins_phi, phi_min, phi_max);
    TH1D *hLambdaPolZ_phiReco_Weighted = new TH1D("hLambdaPolZ_phiReco_Weighted", "hLambdaPolZ_phiReco_Weighted", N_bins_phi, phi_min, phi_max);

        // Testing some ring observable calculations:
    auto hRingObservable_TrueValue = ROOT::TThreadedObject<TH1D>("hRingObservable_TrueValue", "hRingObservable_TrueValue", N_bins_phi, phi_min_Ring, phi_max_Ring);
    auto hRingObservable_TrueValuePtCuts = ROOT::TThreadedObject<TH1D>("hRingObservable_TrueValuePtCuts", "hRingObservable_TrueValuePtCuts", N_bins_phi, phi_min_Ring, phi_max_Ring);
    auto hRingObservable_NormPolTrueValuePtCuts = ROOT::TThreadedObject<TH1D>("hRingObservable_NormPolTrueValuePtCuts", "hRingObservable_NormPolTrueValuePtCuts", N_bins_phi, phi_min_Ring, phi_max_Ring);
    TH1D *hRingObservableReco = new TH1D("hRingObservableReco", "hRingObservableReco", N_bins_phi, phi_min_Ring, phi_max_Ring); // This doesn't need threading!
    TH1D *hRingObservableRecoPtCuts = new TH1D("hRingObservableRecoPtCuts", "hRingObservableRecoPtCuts", N_bins_phi, phi_min_Ring, phi_max_Ring);
    TH1D *hRingObservableNormPolRecoPtCuts = new TH1D("hRingObservableNormPolRecoPtCuts", "hRingObservableNormPolRecoPtCuts", N_bins_phi, phi_min_Ring, phi_max_Ring);
    auto hRingObservable_proxy_from_daughter = ROOT::TThreadedObject<TH1D>("hRingObservable_proxy_from_daughter", "hRingObservable_proxy_from_daughter", N_bins_phi, phi_min_Ring, phi_max_Ring);
    auto hRingObservable_proxy_from_daughter_eq_def = ROOT::TThreadedObject<TH1D>("hRingObservable_proxy_from_daughter_eq_def", "hRingObservable_proxy_from_daughter_eq_def", N_bins_phi, phi_min_Ring, phi_max_Ring);
    auto hRingObservable_proxy_from_daughter_eq_defPtCuts = ROOT::TThreadedObject<TH1D>("hRingObservable_proxy_from_daughter_eq_defPtCuts", "hRingObservable_proxy_from_daughter_eq_defPtCuts", N_bins_phi, phi_min_Ring, phi_max_Ring);
    auto hRingObservable_proxy_from_daughter_eq_def2PtCuts = ROOT::TThreadedObject<TH1D>("hRingObservable_proxy_from_daughter_eq_def2PtCuts", "hRingObservable_proxy_from_daughter_eq_def2PtCuts", N_bins_phi, phi_min_Ring, phi_max_Ring);
    auto hRingObservable_proxy_from_daughter_star_eq_def = ROOT::TThreadedObject<TH1D>("hRingObservable_proxy_from_daughter_star_eq_def", "hRingObservable_proxy_from_daughter_star_eq_def", N_bins_phi, phi_min_Ring, phi_max_Ring);
    auto hRingObservable_proxy_from_daughter_star_eq_defPtCuts = ROOT::TThreadedObject<TH1D>("hRingObservable_proxy_from_daughter_star_eq_defPtCuts", "hRingObservable_proxy_from_daughter_star_eq_defPtCuts", N_bins_phi, phi_min_Ring, phi_max_Ring);

        // Integrated ring observable histograms:
    auto hRingObservable_TrueValuePtCuts_integrated = ROOT::TThreadedObject<TH1D>("hRingObservable_TrueValuePtCuts_integrated", "hRingObservable_TrueValuePtCuts_integrated", 1, phi_min_Ring, phi_max_Ring);
    auto hRingObservable_NormPolTrueValuePtCuts_integrated = ROOT::TThreadedObject<TH1D>("hRingObservable_NormPolTrueValuePtCuts_integrated", "hRingObservable_NormPolTrueValuePtCuts_integrated", 1, phi_min_Ring, phi_max_Ring);
    TH1D *hRingObservableRecoPtCuts_integrated = new TH1D("hRingObservableRecoPtCuts_integrated", "hRingObservableRecoPtCuts_integrated", 1, phi_min_Ring, phi_max_Ring);
    TH1D *hRingObservableNormPolRecoPtCuts_integrated = new TH1D("hRingObservableNormPolRecoPtCuts_integrated", "hRingObservableNormPolRecoPtCuts_integrated", 1, phi_min_Ring, phi_max_Ring);
    auto hRingObservable_proxy_from_daughter_eq_def2PtCuts_integrated = ROOT::TThreadedObject<TH1D>("hRingObservable_proxy_from_daughter_eq_def2PtCuts_integrated", "hRingObservable_proxy_from_daughter_eq_def2PtCuts_integrated", 1, phi_min_Ring, phi_max_Ring);
    auto hRingObservable_proxy_from_daughter_star_eq_defPtCuts_integrated = ROOT::TThreadedObject<TH1D>("hRingObservable_proxy_from_daughter_star_eq_defPtCuts_integrated", "hRingObservable_proxy_from_daughter_star_eq_defPtCuts_integrated", 1, phi_min_Ring, phi_max_Ring);

    auto hDebugCounter_phi_star_sampler = ROOT::TThreadedObject<TH1D>("hDebugCounter_phi_star_sampler", "hDebugCounter_phi_star_sampler", N_bins_phi, phi_min_Ring, phi_max_Ring);

    // The following block of code was moved inside the resamplings loop to take advantage of parallelization capabilities:
    // ////////////////////
    // //// 1.1 - Defining randomized samplers for the loop -- It would make no sense to define them inside the loop, as they will not change:
    // ////////////////////
    // std::mt19937 rng {std::random_device{}()};  // Creates a random_device object (temporary) with the uniform initialization on newer C++ (the "{}"),
    //                                             // then calls it with () getting a random seed, and then that seed is passed to the random engine.
    //     // Uniform distribution for angle xi_star between 0 and pi:
    // // std::uniform_real_distribution<double> dist_x(0.0, PI); // Deprecated!
    // std::uniform_real_distribution<double> dist_unit(0.0, 1.0);
    // std::uniform_real_distribution<double> dist_azimuth(0.0, 2*PI); // Can still sample from 0 to 2pi instead of [phi_min, phi_max) because this is just a liberality within the histogram construction

    ////////////////////
    //// 2 - Getter function calls:
    ////////////////////
    std::cout << "\nFetching data..." << std::endl;
    DoubleMatrix y_matrix, phi_matrix, px_matrix, py_matrix, pz_matrix, pT_matrix, E_matrix, mult_matrix, St_matrix, Sx_matrix, Sy_matrix, Sz_matrix (N_events, std::vector<double>(N_bins_total));
    get_lambda(y_matrix, phi_matrix, px_matrix, py_matrix, pz_matrix, pT_matrix, E_matrix, mult_matrix, St_matrix, Sx_matrix, Sy_matrix, Sz_matrix, N_events, with_bullet);
    std::vector<double> n_event, n_random, phi_random, momentum_x, momentum_y;
    // if (with_bullet){
        get_jet(n_event, n_random, phi_random, momentum_x, momentum_y);

        // Transform each event's phi coordinates (convert phi_matrix into a (\phi - \phi_J) coordinate matrix):
            // A transform that makes it so \hat{t} always points in (1, 0, 0). This makes the RingObservable calculation work across events, and the average polarization is defined
            // in a coordinate system that is shared across all events properly.
            // todo: fix this rotation! Without it, the Ring Observable calculated through polarization is meaningless!
            
            ///////////////////////////////////////////////////////////////////////////////////
            // This should not cause any problems with the averaging of the no-bullet events, 
            // and is even needed for the ring observable to average to zero in an EbE basis,
            // because we need to make the v_2 expansion direction fluctuate!
            // Essentially, this is what guarantees we have (\phi - \phi_J) for the Ring
            // observable plots instead of just a meaningless \phi.
            ///////////////////////////////////////////////////////////////////////////////////
        for (size_t evt = 0; evt < phi_matrix.size(); evt++){ // Event-level loop
            double phi_jet = phi_random[evt];
            phi_random[evt] = 0.; // Rotates the jet to zero for the formulae that depend on this value
            for (size_t i = 0; i < phi_matrix[evt].size(); i++){ // Particle-level loop
                double shifted = phi_matrix[evt][i] - phi_jet;
                
                    // Then we should wrap the angles into the [phi_min, phi_max) interval, because that is the way in which the counting histograms are defined.
                    // This avoids passing data into overflow or underflow bins, and properly converts it to the same range that was used for the iSS histograms.
                phi_matrix[evt][i] = wrapToInterval(shifted, phi_min, phi_max);

                // Rotating the Px and Py matrices in a similar fashion to what will be done to the Sx and Sy matrices:
                double temp_px = px_matrix[evt][i];
                px_matrix[evt][i] = std::cos(-phi_jet) * px_matrix[evt][i] - std::sin(-phi_jet) * py_matrix[evt][i];
                py_matrix[evt][i] = std::sin(-phi_jet) * temp_px + std::cos(-phi_jet) * py_matrix[evt][i];

                // Also, rotating the Sx and Sy spatial components to match the new angle definition (Sz should not be affect by rotations curling the Z axis):
                // (this has to be an active rotation because we are rotating the vectors into a new coordinate system, not the coordinate vectors but the physical vectors themselves!)
                double temp_Sx = Sx_matrix[evt][i]; // Should use the same Sx value for the Sy_matrix rotation too, not the value from the line below after rotation!
                Sx_matrix[evt][i] = std::cos(-phi_jet) * Sx_matrix[evt][i] - std::sin(-phi_jet) * Sy_matrix[evt][i];
                Sy_matrix[evt][i] = std::sin(-phi_jet) * temp_Sx + std::cos(-phi_jet) * Sy_matrix[evt][i];
            }
        }
    // }
    std::cout << "Done fetching data!" << std::endl;
    hLambdaCounter->Fill(0., N_resamples * N_events * mult_matrix[0].size()); // In the way this was declared, it does not count particle number, but the number of total bins that we went through. Declared as a simple TH1D.
        // Need to fill this histogram above with 0. because just a 0 can get interpreted as a char, so the constructor gets ambiguous.

        // Calculating fetcher elapsed time:
    auto data_fetch_end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_fetch = data_fetch_end_time - start_time;
    std::cout << "Spent " << duration_fetch.count() << " seconds fetching data." << std::endl;

    // A loop over all entries of the DoubleMatrix objects (multiplicity, defining polarization and the sorts) to clean and speed things up in the loop:
    DoubleMatrix PolX_matrix = Sx_matrix, PolY_matrix = Sy_matrix, PolZ_matrix = Sz_matrix, PolT_matrix = St_matrix; // Declaring and initializing them in the same shape as the spin matrices

        // Defining vectors that a global bin number equivalent indicator to which each DoubleMatrix contributes to
        // (Fill() would have to search where to place each entry among the 30.000 available bins, and that is a process that
        // can be done only once for all events and all resamples!)
            // Also useful for the buffering vectors for Filling the histograms on a per-event basis instead of per-particle, to reduce memory accesses.
            // Declared as vectors instead of C arrays so that I can have a thread-local copy of them later on to increase access speeds.
    int N_bins_total = N_bins_phi*N_bins_pT*N_bins_rap;
    std::vector<int> compact_bin_idx_pT_y_phi(N_events * N_bins_total, 0); // Rewrote these as 1D vectors with linear indexing!
    std::vector<int> compact_bin_idx_phi_pT(N_events * N_bins_total, 0);
    std::vector<int> compact_bin_idx_pT_y(N_events * N_bins_total, 0);
    std::vector<int> compact_bin_idx_phi(N_events * N_bins_total, 0);
    std::vector<int> compact_bin_idx_phiRing(N_events * N_bins_total, 0);
    // int global_bin_idx_DeltaphiJRing[N_bins_total];

    for (int ev_idx = 0; ev_idx < N_events; ev_idx++){
        int ev_idx_times_N_bins_total = ev_idx * N_bins_total; // This only needs to be calculated once per event loop!
        for (int particle_idx = 0; particle_idx < mult_matrix[ev_idx].size(); particle_idx++){
            // First, converting the "multiplicity" matrix from 1/pT dN/dpTdydphi into actual multiplicity values, 
            // based on "h_Sz_delta_phi_pT_sum->Fill(delta_phi[ip], pT[ip], Sz[ip]*mult[ip]*pT[ip]*dy);" from Vitor's Pol_Analysis_One_Event_hist.C, line 263
            double bin_pT = std::sqrt(px_matrix[ev_idx][particle_idx]*px_matrix[ev_idx][particle_idx] + py_matrix[ev_idx][particle_idx]*py_matrix[ev_idx][particle_idx]);
            // std::cout << "Before: " << mult_matrix[ev_idx][particle_idx] << "(pT, dpT, drap, dphi) = (" << bin_pT << ", " << dpT << ", " << drap << ", " << dphi << ")" << std::endl;
            mult_matrix[ev_idx][particle_idx] *= bin_pT * dpT * drap * dphi;
            // std::cout << "After: " << mult_matrix[ev_idx][particle_idx] << std::endl;

            // Now attributing the values of the polarization matrices (twice the spin)
            PolX_matrix[ev_idx][particle_idx] = Sx_matrix[ev_idx][particle_idx] * 2;
            PolY_matrix[ev_idx][particle_idx] = Sy_matrix[ev_idx][particle_idx] * 2;
            PolZ_matrix[ev_idx][particle_idx] = Sz_matrix[ev_idx][particle_idx] * 2;
            PolT_matrix[ev_idx][particle_idx] = St_matrix[ev_idx][particle_idx] * 2;
   
            // Getting the global index of the current [particle_idx] for each kind of histogram:
                // Do notice that even though the phi, y and pT matrices are the same for all events
                // (the only change is in the multiplicity, and S^\mu matrices), when you calculate
                // the coordinate transform to DeltaPhiJ, each event gets a different phi binning.
            int global_bin_idx_pT_y_phi = hLambdaCounter_pT_y_phi_Weighted.Get()->FindBin(pT_matrix[ev_idx][particle_idx], y_matrix[ev_idx][particle_idx], phi_matrix[ev_idx][particle_idx]);
            int global_bin_idx_phi_pT = hLambdaCounter_phi_pT_Weighted.Get()->FindBin(phi_matrix[ev_idx][particle_idx], pT_matrix[ev_idx][particle_idx]);
            int global_bin_idx_pT_y = hLambdaCounter_pT_y_Weighted.Get()->FindBin(pT_matrix[ev_idx][particle_idx], y_matrix[ev_idx][particle_idx]);
            int global_bin_idx_phi = hLambdaCounter_phi_Weighted.Get()->FindBin(phi_matrix[ev_idx][particle_idx]);
            int global_bin_idx_phiRing = hLambdaCounter_phiRingAngles_Weighted.Get()->FindBin(wrapToInterval(phi_matrix[ev_idx][particle_idx], phi_min_Ring, phi_max_Ring)); // Be cautious with this! The global bin number will change in this coordinate system here!
            // global_bin_idx_DeltaphiJRing[particle_idx] = hLambdaCounter_DeltaphiJRing_Weighted.Get()->FindBin(wrapToInterval(phi_matrix[ev_idx][particle_idx], phi_min_Ring, phi_max_Ring));

            // Because ROOT's GlobalIndex includes overflow/underflow bins between 1 to N_bins_total (even if those bins are empty in this case!), I need to map these indices
            // into the [0, N_bins_total) interval in order to skip overflow/underflow numberings and have an index that can be used for the buffer vector:
                // (I want a linear index from 0 to N_bins_total - 1)
                // Mapping each index into a compact index -- Also notice you have to dereference the shared pointer for the function to work:
            compact_bin_idx_pT_y_phi[ev_idx_times_N_bins_total + particle_idx] = mapGlobalToCompact(*(hLambdaCounter_pT_y_phi_Weighted.Get()), global_bin_idx_pT_y_phi);
            compact_bin_idx_phi_pT[ev_idx_times_N_bins_total + particle_idx] = mapGlobalToCompact(*(hLambdaCounter_phi_pT_Weighted.Get()), global_bin_idx_phi_pT);
            compact_bin_idx_pT_y[ev_idx_times_N_bins_total + particle_idx] = mapGlobalToCompact(*(hLambdaCounter_pT_y_Weighted.Get()), global_bin_idx_pT_y);
            compact_bin_idx_phi[ev_idx_times_N_bins_total + particle_idx] = mapGlobalToCompact(*(hLambdaCounter_phi_Weighted.Get()), global_bin_idx_phi);
            compact_bin_idx_phiRing[ev_idx_times_N_bins_total + particle_idx] = mapGlobalToCompact(*(hLambdaCounter_phiRingAngles_Weighted.Get()), global_bin_idx_phiRing);

            // // DEBUG! Checking if the compact index is within the expected margins (should go from 0 to N_bins-1)
            // if(compact_bin_idx_pT_y_phi[ev_idx_times_N_bins_total + particle_idx] < 0 || compact_bin_idx_pT_y_phi[ev_idx_times_N_bins_total + particle_idx] >= N_bins_total){
            //     std::cout << "global_bin_idx_pT_y_phi: " << global_bin_idx_pT_y_phi << std::endl;
            //     std::cout << "compact_bin_idx_pT_y_phi[ev_idx_times_N_bins_total + particle_idx]: " << compact_bin_idx_pT_y_phi[ev_idx_times_N_bins_total + particle_idx] << std::endl;
            // }
            // if(compact_bin_idx_phi_pT[ev_idx_times_N_bins_total + particle_idx] < 0 || compact_bin_idx_phi_pT[ev_idx_times_N_bins_total + particle_idx] >= (N_bins_phi*N_bins_pT)){
            //     std::cout << "global_bin_idx_phi_pT: " << global_bin_idx_phi_pT << std::endl;
            //     std::cout << "compact_bin_idx_phi_pT[ev_idx_times_N_bins_total + particle_idx]: " << compact_bin_idx_phi_pT[ev_idx_times_N_bins_total + particle_idx] << std::endl;
            // }
        }
    }

    // Testing to see if the compact indices are indeed restricted from 0 to N_bins_total:
    // for (const auto& num : compact_bin_idx_pT_y_phi){std::cout << num << " ";} std::cout << std::endl;

        // Calculating pre-processer elapsed time:
    auto pre_process_end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> pre_process_duration = pre_process_end_time - start_time - duration_fetch;
    std::cout << "Spent " << pre_process_duration.count() << " seconds pre-processing data for the multi-threaded loop." << std::endl;

    // // DEBUG!
    // // Checking to see if any invalid global indices were included in the global_bin arrays:
    // // Generic checker for C arrays
    // /////////////////////////////////////////////
    // // Example usage (after arrays are filled)
    // CheckBins("pT_y_phi", global_bin_idx_pT_y_phi, N_bins_total, hLambdaCounter_pT_y_phi_Weighted.Get().get());
    // CheckBins("phi_pT",   global_bin_idx_phi_pT,   N_bins_total, hLambdaCounter_phi_pT_Weighted.Get().get());
    // CheckBins("pT_y",     global_bin_idx_pT_y,     N_bins_total, hLambdaCounter_pT_y_Weighted.Get().get());
    // CheckBins("phi",      global_bin_idx_phi,      N_bins_total, hLambdaCounter_phi_Weighted.Get().get());
    // CheckBins("phiRing",  global_bin_idx_phiRing,  N_bins_total, hLambdaCounter_phiRingAngles_Weighted.Get().get());
    // /////////////////////////////////////////////

    ////////////////////
    //// 3 - Decaying Lambdas and filling histograms:
    ////////////////////
    std::cout << "\nDecaying Lambdas and filling histograms" << std::endl;
        // Looping on all events and all bins(/particles with multiplicity) in each event:
    // int N_resamples = 150; // Goes through each particle N_resamples # of times. This is an attempt to see if the polarization estimate values become more stable!
                             // 55000 is 500 resamplings for each of the 110 threads, which should take about 450 minutes

    // This resampling loop can be parallelized!
    // const int N_threads = 1;
    omp_set_num_threads(N_threads);

        // Now actually declaring the parallelization loop:
    #pragma omp parallel for
    for (int resample_idx = 0; resample_idx < N_resamples; resample_idx++){
        #pragma omp critical
        {
            std::cout << "\tNow on resample " << std::to_string(resample_idx + 1) << " of " + std::to_string(N_resamples) << std::endl;
            // todo: fix this to print only once per resampling batch!
        }

        // Declaring thread-local variables:
            // Each thread should get its own RNG device, seeded differently, to avoid resamplings that share the same random numbers!
        std::random_device rd; // Cleaner than the earlier method
        std::mt19937 rng(rd() + omp_get_thread_num()); // Creates a random_device object then calls it with () getting a random seed (which is summed
                                                    // to a number related to the current thread, to ensure each thread has a different seed), and
                                                    // then that seed is passed to the random engine.
            // Uniform distribution for sampling cos(xi_star) using the inverse CDF method:
        std::uniform_real_distribution<double> dist_unit(0.0, 1.0);
            // Azimuth distribution:
        std::uniform_real_distribution<double> dist_azimuth(0.0, 2*PI); // Can still sample from 0 to 2pi instead of [phi_min, phi_max) because this is just a liberality within the histogram construction

        // Removed this block because it gave no real performance benefit, yet spent a lot of memory!
        // // Declaring thread-local copies of the vectors and matrices to avoid memory speed bottlenecks:
        //     // (I believe the bottleneck in parallelization was that each thread tried reading the same address in memory,
        //     // where the data matrices where stored, and each thread was trying to read a different part of that vector at
        //     // different times. This could be causing the access to become slow, as seen in the HUGE wall-time from mem 
        //     // access via perf stat command.)
        // auto px_matrix_thrd_copy = px_matrix; // Am not declaring these guys with thread_local keywords because that causes a bus error
        // auto py_matrix_thrd_copy = py_matrix;
        // auto pz_matrix_thrd_copy = pz_matrix;
        // auto E_matrix_thrd_copy = E_matrix;
        //
        // auto y_matrix_thrd_copy = y_matrix;
        // auto pT_matrix_thrd_copy = pT_matrix;
        // auto phi_matrix_thrd_copy = phi_matrix;
        //
        // auto mult_matrix_thrd_copy = mult_matrix;
        // auto PolX_matrix_thrd_copy = PolX_matrix;
        // auto PolY_matrix_thrd_copy = PolY_matrix;
        // auto PolZ_matrix_thrd_copy = PolZ_matrix;
        // auto PolT_matrix_thrd_copy = PolT_matrix;
        //
        // auto phi_random_thrd_copy = phi_random;
        //
        // auto compact_bin_idx_pT_y_phi = compact_bin_idx_pT_y_phi_non_Local; // Don't expect that much speed gain with this, so removed to avoid excessive memory usage...
        // auto compact_bin_idx_phi_pT = compact_bin_idx_phi_pT_non_Local;
        // auto compact_bin_idx_pT_y = compact_bin_idx_pT_y_non_Local;
        // auto compact_bin_idx_phi = compact_bin_idx_phi_non_Local;
        // auto compact_bin_idx_phiRing = compact_bin_idx_phiRing_non_Local;

        // Declaring buffer arrays:
            // It is faster to fill the histogram contiguously and only once for each event,instead of filling for each DoubleMatrix entry,
            // because of the other updates performed by TH3D when filling that are more complex than simply adding a number to the counters.
            // All buffers are initialized with 0 and are dumped into the related histograms after all events have been processed.
            // These are indexed with GlobalBinIndexes, which should properly go from 1 to N_bins_total.
            // Using C++ vectors instead of C arrays to use heap memory instead of stack (friendlier for large vectors inside )
        // The buffers correspond 1:1 with the histograms, so you don't need smaller histograms' buffers to go up to N_bins_total: they can have less entries!
        std::vector<double> hLambdaCounter_pT_y_phi_Weighted_buffer(N_bins_total, 0.0);
        std::vector<double> hLambdaAvgDotX_pT_y_phi_Weighted_buffer(N_bins_total, 0.0);
        std::vector<double> hLambdaAvgDotY_pT_y_phi_Weighted_buffer(N_bins_total, 0.0);
        std::vector<double> hLambdaAvgDotZ_pT_y_phi_Weighted_buffer(N_bins_total, 0.0);
        std::vector<double> hLambdaCounter_pT_y_Weighted_buffer(N_bins_pT*N_bins_rap, 0.0); 
        std::vector<double> hLambdaPolX_phi_pT_Weighted_buffer(N_bins_phi*N_bins_pT, 0.0);
        std::vector<double> hLambdaPolY_phi_pT_Weighted_buffer(N_bins_phi*N_bins_pT, 0.0);
        std::vector<double> hLambdaPolZ_phi_pT_Weighted_buffer(N_bins_phi*N_bins_pT, 0.0);
        std::vector<double> hLambdaCounter_phi_pT_Weighted_buffer(N_bins_phi*N_bins_pT, 0.0);
        std::vector<double> hLambdaCounter_phi_Weighted_buffer(N_bins_phi, 0.0);
        std::vector<double> hLambdaCounter_phiRingAngles_Weighted_buffer(N_bins_phi, 0.0);
        // double hLambdaCounter_DeltaphiJRing_Weighted_buffer[N_bins_total]{};
        // double hProtonCounter_phiRing_Weighted_buffer[N_bins_total]{}; // These three benefit from not having a buffer because their data isn't pre-binned!
        // double hProtonStarCounter_phiRing_Weighted_buffer[N_bins_total]{};
        // double hDebugCounter_phi_star_sampler_buffer[N_bins_total]{};
        std::vector<double> hLambdaPolX_phi_Weighted_buffer(N_bins_phi, 0.0);
        std::vector<double> hLambdaPolY_phi_Weighted_buffer(N_bins_phi, 0.0);
        std::vector<double> hLambdaPolZ_phi_Weighted_buffer(N_bins_phi, 0.0);
        std::vector<double> hLambdaPolX_pT_y_phi_Weighted_buffer(N_bins_total, 0.0);
        std::vector<double> hLambdaPolY_pT_y_phi_Weighted_buffer(N_bins_total, 0.0);
        std::vector<double> hLambdaPolZ_pT_y_phi_Weighted_buffer(N_bins_total, 0.0);
        std::vector<double> hRingObservable_proxy_from_daughter_buffer(N_bins_phi, 0.0);
        std::vector<double> hRingObservable_proxy_from_daughter_eq_def_buffer(N_bins_phi, 0.0);
        std::vector<double> hRingObservable_proxy_from_daughter_star_eq_def_buffer(N_bins_phi, 0.0);
        std::vector<double> hRingObservable_TrueValue_buffer(N_bins_phi, 0.0);
        std::vector<double> hRingObservable_TrueValuePtCuts_buffer(N_bins_phi, 0.0);
        std::vector<double> hRingObservable_NormPolTrueValuePtCuts_buffer(N_bins_phi, 0.0);
        std::vector<double> hRingObservable_TrueValuePtCuts_integrated_buffer(N_bins_phi, 0.0);
        std::vector<double> hRingObservable_NormPolTrueValuePtCuts_integrated_buffer(N_bins_phi, 0.0);
        std::vector<double> hRingObservable_proxy_from_daughter_eq_defPtCuts_buffer(N_bins_phi, 0.0);
        std::vector<double> hRingObservable_proxy_from_daughter_eq_def2PtCuts_buffer(N_bins_phi, 0.0);
        std::vector<double> hRingObservable_proxy_from_daughter_star_eq_defPtCuts_buffer(N_bins_phi, 0.0);
        std::vector<double> hRingObservable_proxy_from_daughter_eq_def2PtCuts_integrated_buffer(N_bins_phi, 0.0);
        std::vector<double> hRingObservable_proxy_from_daughter_star_eq_defPtCuts_integrated_buffer(N_bins_phi, 0.0);
        // std::cout << "DEBUG! Thread locals declared" << std::endl;

        // DEBUG!
        // std::cout << y_matrix_thrd_copy[100].size() << " " << N_bins_total << std::endl; // Yes, these two values did match!

        for (int ev_idx = 0; ev_idx < N_events; ev_idx++){
            int ev_idx_times_N_bins_total = ev_idx * N_bins_total; // This only needs to be calculated once per event loop!
            // if (ev_idx % int (N_events * 0.1) == 0){ // Keeping track of the proccess for every 10% of the events
            //     std::cout << "\t\tNow on event " << ev_idx << " of " << N_events  << " (" << ev_idx * 1./N_events * 100 << "%)" << std::endl; // The 1. comes to represent the percentage as a double
            // }
            for (int particle_idx = 0; particle_idx < y_matrix[ev_idx].size(); particle_idx++){ // Not truly a particle loop: the data is pre-binned, so this is a loop on bins of each event, which can be treated as "mean particles" of each (pT, phi, y) bin/interval
                // Getiting the global index for each type of buffer for this particle bin:
                    // Notice you don't have to subtract 1 from the index because mapGlobalToCompact() already does that!
                int compact_bin_idx_pT_y_phi_value = compact_bin_idx_pT_y_phi[ev_idx_times_N_bins_total + particle_idx];
                int compact_bin_idx_phi_pT_value = compact_bin_idx_phi_pT[ev_idx_times_N_bins_total + particle_idx];
                int compact_bin_idx_pT_y_value = compact_bin_idx_pT_y[ev_idx_times_N_bins_total + particle_idx];
                int compact_bin_idx_phi_value = compact_bin_idx_phi[ev_idx_times_N_bins_total + particle_idx];
                int compact_bin_idx_phiRing_value = compact_bin_idx_phiRing[ev_idx_times_N_bins_total + particle_idx];

                if ((compact_bin_idx_pT_y_phi_value < 0 || compact_bin_idx_pT_y_phi_value >= N_bins_total) ||
                    (compact_bin_idx_phi_pT_value < 0 || compact_bin_idx_phi_pT_value >= N_bins_phi*N_bins_pT) ||
                    (compact_bin_idx_pT_y_value < 0 || compact_bin_idx_pT_y_value >= N_bins_pT*N_bins_rap) ||
                    (compact_bin_idx_phi_value < 0 || compact_bin_idx_phi_value >= N_bins_phi) ||
                    (compact_bin_idx_phiRing_value < 0 || compact_bin_idx_phiRing_value >= N_bins_phi)){
                    
                    std::cout << "Testing compact bin indices:" << std::endl;
                    std::cout << "compact_bin_idx_pT_y_phi_value: " << compact_bin_idx_pT_y_phi_value << " max should be " << N_bins_total << std::endl;
                    std::cout << "compact_bin_idx_phi_pT_value: " << compact_bin_idx_phi_pT_value << " max should be " << N_bins_phi*N_bins_pT << std::endl;
                    std::cout << "compact_bin_idx_pT_y_value: " << compact_bin_idx_pT_y_value << " max should be " << N_bins_pT*N_bins_rap << std::endl;
                    std::cout << "compact_bin_idx_phi_value: " << compact_bin_idx_phi_value << " max should be " << N_bins_phi << std::endl;
                    std::cout << "compact_bin_idx_phiRing_value: " << compact_bin_idx_phiRing_value << " max should be " << N_bins_phi << std::endl;
                }

                double current_particle_multiplicity = mult_matrix[ev_idx][particle_idx];
                double true_PolX = PolX_matrix[ev_idx][particle_idx];
                double true_PolY = PolY_matrix[ev_idx][particle_idx];
                double true_PolZ = PolZ_matrix[ev_idx][particle_idx];
                double true_PolT = PolT_matrix[ev_idx][particle_idx];

                // 1 - Fetching particle information:
                TLorentzVector Lambda_4_momentum_lab(px_matrix[ev_idx][particle_idx], py_matrix[ev_idx][particle_idx], pz_matrix[ev_idx][particle_idx], E_matrix[ev_idx][particle_idx]);
                // TVector3 P_Lambda_lab(Sx_matrix[ev_idx][particle_idx], Sy_matrix[ev_idx][particle_idx], Sz_matrix[ev_idx][particle_idx]);
                    // Polarization is given by P = S/<S> in the case of Lambdas with spin-1/2, and as <S> = 1/2, you just need to multiply everything by 2.
                TLorentzVector P_Lambda_lab_4vec(true_PolX, true_PolY, true_PolZ, true_PolT); // You need the temporal component too!

                ///////////////////////////////////////////
                // // Making sure that the S^mu * P_mu product is zero.
                // // Compute scalar product S.P
                // double dotProduct = P_Lambda_lab_4vec * Lambda_4_momentum_lab;
                //
                // // Report
                // std::cout << "S^mu * P_mu = " << dotProduct << std::endl;
                // It was close enough to zero! (within float precision, which seems to be the format the data was previously stored in)
                ///////////////////////////////////////////

                // 2 - Calculating the polarization in the Lambda rest frame:
                TVector3 P_Lambda_star = boost_polarization_to_rest_frame(Lambda_4_momentum_lab, P_Lambda_lab_4vec);
                double P_Lambda_star_mag = P_Lambda_star.Mag();

                ///////////////////////////////////////////
                // // Verifying if the polarization magnitude is indeed close to one, or if I am dealing with excess polarizations already
                // std::cout << "Current bin's P_lambda_star magnitude is: " << P_Lambda_star_mag << std::endl;
                ///////////////////////////////////////////

                // 3 - Sampling decay angles:
                // // Older code that used rejection sampling:
                //     // Calculating the P_max value for the rejection sampling:
                //     // Compute P_max (the maximum possible value of P(x))
                //     // This occurs at x = 0 or x = pi depending on sign of cos(x) term
                // double P_max = (1.0 / PI) * (1.0 + std::abs(alpha_H * P_Lambda_star_mag));

                //     // Defining the uniform distribution dist_y for the sampling:
                //     // Uniform distribution for y between 0 and P_max
                // std::uniform_real_distribution<double> dist_y(0.0, P_max);

                // auto [xi_star, phi_star] = sample_P_angle_proton(P_Lambda_star_mag, rng, dist_x, dist_y, dist_azimuth); // Uses the new unpacking of C++17

                // Newer code, that samples directly from the cos_xi distribution:
                auto [xi_star, phi_star] = sample_P_angle_proton_from_cos_xi(P_Lambda_star_mag, rng, dist_unit, dist_azimuth); // Optimized version that does not require

                //////////////////////////////////////////////
                /// DEBUG ONLY! Testing to see if xi_star = 0
                /// (proton always decays in the direction of
                /// the polarization) gives us a working result,
                /// i.e., recovers the true Ring Observable
                /// as it should from both the polarization
                /// and the proton momentum proxy for the ring
                /// observable
                // force_DecayDist_mode == 0 for the usual, physical, decay.
                if (force_DecayDist_mode == 1){ // 1 for decays along Pol.
                    xi_star = 0; // Should be so that the polarization vector and the proton trimomentum point in the same direction
                }
                else if (force_DecayDist_mode == 2){ // 2 for decays perpendicular to Pol
                    xi_star = PI/2; // Should properly give us a zero ring observable signal!
                }
                //////////////////////////////////////////////


                // 4 - Generating the proton 4-momentum from the decay:
                    // Actually, I just need the angles at which the proton would decay, not the whole 4-momentum of the decay, but whatever, let's keep it!
                    // In other words, I could've stopped at the sample_P_angle_proton function, and then just rotate those angles to the XYZ axes of the lab frame.
                    // This could turn out to be useful later on, if I intend on doing some background checks for the Lambda reconstructions or something like that.
                    // You won't even need the proton's momentum for the ring observable's reconstruction! But whatever...
                auto [proton_4_momentum, proton_4_momentum_star] = Lambda_decay(Lambda_4_momentum_lab, P_Lambda_star, xi_star, phi_star); // Actually don't use the proton_4_momentum variable, but whatever

                // 5 - Extracting useful variables from the 4-momentum of the Lambda:              
                    // Don't need to calculate these numbers! They are already provided in the code! (But I checked and the results match)
                double lambda_y = y_matrix[ev_idx][particle_idx];
                double lambda_pT = pT_matrix[ev_idx][particle_idx];
                double lambda_phi = phi_matrix[ev_idx][particle_idx];

                // 6 - Summing to average the polarization vector -- Global polarization:
                    // To do the average dot product between the proton momentum direction and \hat{n} when there are (pT, y) bins involved,
                    // it would probably be easier to create a TH2D that collects the sums, and then normalize each bin by the number of times
                    // I added values into it. The bins will not have counts, but the actual value of the sum. If I could just create a matrix
                    // that has those summation values and the appropriate bin limits, then that would be equivalent.
                    // A good idea would be to have a TH2D that has the sums for the averages, and another that has the number of particles in
                    // that specific bin
                // TVector3 proton_unit_vector = (proton_4_momentum.Vect()).Unit();
                    // You actually need the proton vector in the Lambda frame of reference, and you will get a polarization vector on that frame, which will need to be boosted later.
                TVector3 proton_star_unit_vector = (proton_4_momentum_star.Vect()).Unit();
                // double X_dot = proton_star_unit_vector.Dot(x_hat);
                // double Y_dot = proton_star_unit_vector.Dot(y_hat);
                // double Z_dot = proton_star_unit_vector.Dot(z_hat);
                    // Instead of actually calculating dot products, you could just take proton_star_unit_vector.X(), proton_star_unit_vector.Y() and proton_star_unit_vector.Z() !!!
                double X_dot = proton_star_unit_vector.X();
                double Y_dot = proton_star_unit_vector.Y();
                double Z_dot = proton_star_unit_vector.Z();
                
                ///////////////////////
                // std::cout << "The (x_dot, y_dot, z_dot) values for the p_proton * hat vectors are (" << X_dot << ", " << Y_dot << ", " << Z_dot << ")" << std::endl; // Just checking the magnitude of the projections
                ///////////////////////

                // 7 - Summing to calculate polarization on subsets of the available Lambdas -- Bins of (pT, y):
                    // Filling the weighted values:
                hLambdaCounter_pT_y_phi_Weighted_buffer[compact_bin_idx_pT_y_phi_value] += current_particle_multiplicity;
                hLambdaAvgDotX_pT_y_phi_Weighted_buffer[compact_bin_idx_pT_y_phi_value] += X_dot * current_particle_multiplicity; // Weighted averages
                hLambdaAvgDotY_pT_y_phi_Weighted_buffer[compact_bin_idx_pT_y_phi_value] += Y_dot * current_particle_multiplicity;
                hLambdaAvgDotZ_pT_y_phi_Weighted_buffer[compact_bin_idx_pT_y_phi_value] += Z_dot * current_particle_multiplicity;
                
                hLambdaCounter_pT_y_Weighted_buffer[compact_bin_idx_pT_y_value] += current_particle_multiplicity;

                // 9 - Filling a peace of mind plot (plot that has the true values of the MC):
                hLambdaPolX_phi_pT_Weighted_buffer[compact_bin_idx_phi_pT_value] += true_PolX * current_particle_multiplicity;
                hLambdaPolY_phi_pT_Weighted_buffer[compact_bin_idx_phi_pT_value] += true_PolY * current_particle_multiplicity;
                hLambdaPolZ_phi_pT_Weighted_buffer[compact_bin_idx_phi_pT_value] += true_PolZ * current_particle_multiplicity;
                hLambdaCounter_phi_pT_Weighted_buffer[compact_bin_idx_phi_pT_value] += current_particle_multiplicity;

                // Another peace of mind, just in phi:
                hLambdaCounter_phi_Weighted_buffer[compact_bin_idx_phi_value] += current_particle_multiplicity;
                hLambdaCounter_phiRingAngles_Weighted_buffer[compact_bin_idx_phiRing_value] += current_particle_multiplicity; // No the prettiest way to build a histogram that is essentially equivalent to the other counter yet going from -pi to pi, but works
                hLambdaCounter_DeltaphiJRing_Weighted.Get()->Fill(wrapToInterval(lambda_phi - phi_random[ev_idx], phi_min_Ring, phi_max_Ring), current_particle_multiplicity); // You actually need to have this counter defined as lambda_phi - phi_random[ev_idx] because the ring observable has a coordinate shift!
                    // These three histograms here benefit from using the Fill() in place, because the global index can not be predetermined!
                hProtonCounter_phiRing_Weighted.Get()->Fill(wrapToInterval(proton_4_momentum.Phi(), phi_min_Ring, phi_max_Ring), current_particle_multiplicity); // A counter for the proton's angular distribution too, for the other ring estimator!
                hProtonStarCounter_phiRing_Weighted.Get()->Fill(wrapToInterval(proton_4_momentum_star.Phi(), phi_min_Ring, phi_max_Ring), current_particle_multiplicity); // For the other proxy of the ring observable, which 
                                                                                                                                                                    // should be normalized by the number of protons
                hDebugCounter_phi_star_sampler.Get()->Fill(wrapToInterval(phi_star, phi_min_Ring, phi_max_Ring)); // A debug histogram, just to know if the sampling is being truly random in phi_star, before coordinate rotation

                hLambdaPolX_phi_Weighted_buffer[compact_bin_idx_phi_value] += true_PolX * current_particle_multiplicity;
                hLambdaPolY_phi_Weighted_buffer[compact_bin_idx_phi_value] += true_PolY * current_particle_multiplicity;
                hLambdaPolZ_phi_Weighted_buffer[compact_bin_idx_phi_value] += true_PolZ * current_particle_multiplicity;

                // A third plot, this time in all three dimensions, which will later be reduced to 2D to test if the 3D-->2D conversion was implemented correctly:
                hLambdaPolX_pT_y_phi_Weighted_buffer[compact_bin_idx_pT_y_phi_value] += true_PolX * current_particle_multiplicity;
                hLambdaPolY_pT_y_phi_Weighted_buffer[compact_bin_idx_pT_y_phi_value] += true_PolY * current_particle_multiplicity;
                hLambdaPolZ_pT_y_phi_Weighted_buffer[compact_bin_idx_pT_y_phi_value] += true_PolZ * current_particle_multiplicity;

                // Calculating the Ring Observable proxy with the daughter particle's momentum instead of the polarization:
                    // See the other RP_temp variable's surrounding lines for more information
                // While the rotation formula is still not solved:
                double x_trigger = std::cos(phi_random[ev_idx]); // Defined a coordinate system such that double x_trigger = 1 is always true, but kept it general where possible, so that code is flexible to changes
                double y_trigger = std::sin(phi_random[ev_idx]);
                const double z_trigger = 0;

                double cross_x = y_trigger*pz_matrix[ev_idx][particle_idx] - z_trigger*py_matrix[ev_idx][particle_idx];
                double cross_y = z_trigger*px_matrix[ev_idx][particle_idx] - x_trigger*pz_matrix[ev_idx][particle_idx];
                double cross_z = x_trigger*py_matrix[ev_idx][particle_idx] - y_trigger*px_matrix[ev_idx][particle_idx];

                double cross_product_norm = std::sqrt(cross_x*cross_x + cross_y*cross_y + cross_z*cross_z);
                    
                double RP_temp = ((proton_4_momentum.X()*cross_x + proton_4_momentum.Y()*cross_y + proton_4_momentum.Z()*cross_z)/cross_product_norm);

                    // Another good implementation of this, that should be equivalent to the true Ring Observable value, is:
                double proton_4_momentum_star_norm = std::sqrt(proton_4_momentum_star.X()*proton_4_momentum_star.X() + proton_4_momentum_star.Y()*proton_4_momentum_star.Y()
                                                    + proton_4_momentum_star.Z()*proton_4_momentum_star.Z());
                double proton_4_momentum_norm = std::sqrt(proton_4_momentum.X()*proton_4_momentum.X() + proton_4_momentum.Y()*proton_4_momentum.Y()
                                                    + proton_4_momentum.Z()*proton_4_momentum.Z());    
                double RP_temp_eqv_def = 3./(alpha_H * proton_4_momentum_star_norm) * 
                                                ((proton_4_momentum.X()*cross_x + proton_4_momentum.Y()*cross_y + proton_4_momentum.Z()*cross_z)/cross_product_norm);
                double RP_temp_eqv_def2 = 3./(alpha_H * proton_4_momentum_norm) * 
                                                ((proton_4_momentum.X()*cross_x + proton_4_momentum.Y()*cross_y + proton_4_momentum.Z()*cross_z)/cross_product_norm);
                    // Now for the observables that take into consideration the proton's momentum in the Lambda rest frame -- Even more closely related to the polarization vector:
                double RP_temp_proton_star_eqv_def = 3./(alpha_H * proton_4_momentum_star_norm) * 
                                                ((proton_4_momentum_star.X()*cross_x + proton_4_momentum_star.Y()*cross_y + proton_4_momentum_star.Z()*cross_z)/cross_product_norm);

                // double delta_phi_J = lambda_phi - phi_random_thrd_copy[ev_idx]; // Also a general formula: in case phi_random was rotated to 0, this does nothing, but it also works in case it was not rotated.
                //     // Making sure that the values are within 0 to 2*PI appropriately:
                // delta_phi_J = wrapToInterval(delta_phi_J, phi_min_Ring, phi_max_Ring);

                hRingObservable_proxy_from_daughter_buffer[compact_bin_idx_phiRing_value] += RP_temp * current_particle_multiplicity;
                hRingObservable_proxy_from_daughter_eq_def_buffer[compact_bin_idx_phiRing_value] += RP_temp_eqv_def * current_particle_multiplicity;
                hRingObservable_proxy_from_daughter_star_eq_def_buffer[compact_bin_idx_phiRing_value] += RP_temp_proton_star_eqv_def * current_particle_multiplicity;

                // Calculating the true value:
                double RP_temp_true = ((true_PolX*cross_x + true_PolY*cross_y + true_PolZ*cross_z)/cross_product_norm);
                hRingObservable_TrueValue_buffer[compact_bin_idx_phiRing_value] += RP_temp_true * current_particle_multiplicity;

                    // An observable that carries only the direction of the polarization, not the magnitude:
                double Pol_norm = std::sqrt(true_PolX*true_PolX + true_PolY*true_PolY + true_PolZ*true_PolZ);
                double RP_temp_true_Norm = ((true_PolX*cross_x + true_PolY*cross_y + true_PolZ*cross_z)/cross_product_norm)/Pol_norm;

                if (lambda_pT > 0.5 && lambda_pT < 1.5){
                    hRingObservable_TrueValuePtCuts_buffer[compact_bin_idx_phiRing_value] += RP_temp_true * current_particle_multiplicity;
                    hRingObservable_NormPolTrueValuePtCuts_buffer[compact_bin_idx_phiRing_value] += RP_temp_true_Norm * current_particle_multiplicity;

                    hRingObservable_TrueValuePtCuts_integrated_buffer[compact_bin_idx_phiRing_value] += RP_temp_true * current_particle_multiplicity;
                    hRingObservable_NormPolTrueValuePtCuts_integrated_buffer[compact_bin_idx_phiRing_value] += RP_temp_true_Norm * current_particle_multiplicity;

                    hRingObservable_proxy_from_daughter_eq_defPtCuts_buffer[compact_bin_idx_phiRing_value] += RP_temp_eqv_def * current_particle_multiplicity;
                    hRingObservable_proxy_from_daughter_eq_def2PtCuts_buffer[compact_bin_idx_phiRing_value] += RP_temp_eqv_def2 * current_particle_multiplicity;
                    hRingObservable_proxy_from_daughter_star_eq_defPtCuts_buffer[compact_bin_idx_phiRing_value] += RP_temp_proton_star_eqv_def * current_particle_multiplicity;

                    hRingObservable_proxy_from_daughter_eq_def2PtCuts_integrated_buffer[compact_bin_idx_phiRing_value] += RP_temp_eqv_def2 * current_particle_multiplicity;
                    hRingObservable_proxy_from_daughter_star_eq_defPtCuts_integrated_buffer[compact_bin_idx_phiRing_value] += RP_temp_proton_star_eqv_def * current_particle_multiplicity;
                }

                // // Checking/debugging the angles after the coordinate rotation:
                // std::cout << "#Ev: " << ev_idx << ", phi_random[ev_idx]: " << phi_random[ev_idx] << ", std::cos(phi_random[ev_idx]): " << std::cos(phi_random[ev_idx]) << std::endl;
                // // Checking /debugging the momenta, cross products and polarizations after the rotation of phi:
                // if (ev_idx == 5 && particle_idx == 4){
                //     std::cout << "(px, py, pz) = (" << px_matrix[ev_idx][particle_idx] << ", " << py_matrix[ev_idx][particle_idx] << ", " << pz_matrix[ev_idx][particle_idx] << ")" << std::endl;
                //     std::cout << "(cross_x, cross_y, cross_z) = (" << cross_x << ", " << cross_y << ", " << cross_z << ")" << std::endl;
                //     std::cout << "(true_PolX, true_PolY, true_PolZ) = (" << true_PolX << ", " << true_PolY << ", " << true_PolZ << ")" << std::endl;
                //     std::cout << "RP_temp_true: " << RP_temp_true << std::endl;
                // }
            } // End of particle loop
            // std::cout << "Event end " << ev_idx << std::endl;
        } // End of event loop
        // Emptying buffers after all events where processed -- No true gain in doing this histogram filling on a per-event basis...
            // The BufferEmpty anf BufferFill methods don't work for TH3D's because they are protected methods.
            // Thus, created a straightforward method with C++ vectors.
            // Notice you have to use AddBinContent() instead of SetBinContent() because the projection histograms have less dimensions than N_bins_total!
        std::cout<<"DEBUG! Flushing buffers for resample " << resample_idx << std::endl;

        // Doing the flush in a sequential manner, so that there is no chance of multiple access to same object, and guaranteeing thread safety
            // (Not the best solution nor the fastest, but surely a solution)
        #pragma omp critical
        {
            // Will divide the flushing in 4 categories, one for each buffer size:
                // The buffers could all have the same N_bins_total size and work in the same loop, but that would be a wasteful use of RAM
                // 1 - Flushing 3D buffers, which have GlobalBinIndex going from 1 to N_bins_total:
                    // Fetching all histograms only ONCE per loop -- This avoids lazy recreations of pointers by ROOT
                    // Notice you have to declare the type explicitly! The compiler can't handle this by its own.
                    // Also, you have to first get the proxy (the std::shared_ptr<TH3D> that would get out of Get()), then do another .get() to get the raw pointer!
            std::shared_ptr<TH3D> hLambdaCounter_pT_y_phi_Weighted_pointer = hLambdaCounter_pT_y_phi_Weighted.Get();
            std::shared_ptr<TH3D> hLambdaAvgDotX_pT_y_phi_Weighted_pointer = hLambdaAvgDotX_pT_y_phi_Weighted.Get();
            std::shared_ptr<TH3D> hLambdaAvgDotY_pT_y_phi_Weighted_pointer = hLambdaAvgDotY_pT_y_phi_Weighted.Get();
            std::shared_ptr<TH3D> hLambdaAvgDotZ_pT_y_phi_Weighted_pointer = hLambdaAvgDotZ_pT_y_phi_Weighted.Get();
            std::shared_ptr<TH3D> hLambdaPolX_pT_y_phi_Weighted_pointer = hLambdaPolX_pT_y_phi_Weighted.Get();
            std::shared_ptr<TH3D> hLambdaPolY_pT_y_phi_Weighted_pointer = hLambdaPolY_pT_y_phi_Weighted.Get();
            std::shared_ptr<TH3D> hLambdaPolZ_pT_y_phi_Weighted_pointer = hLambdaPolZ_pT_y_phi_Weighted.Get();

            for (int CompactBinIndex_pT_y_phi = 0; CompactBinIndex_pT_y_phi < N_bins_total; CompactBinIndex_pT_y_phi++){
                // No longer using the format below -- Now the index of the buffer itself corresponds to a global bin number, so you can just use that index instead!
                    // In other words, the whole GlobalBin matching procedure was done when filling the buffers!
                // hLambdaCounter_pT_y_phi_Weighted.Get()->AddBinContent(global_bin_idx_pT_y_phi[pre_binned_idx], hLambdaCounter_pT_y_phi_Weighted_buffer[pre_binned_idx]);

                // Converting the CompactBinIndex back to a global index:
                int GlobalBinIndex_pT_y_phi = mapCompactToGlobal(*(hLambdaCounter_pT_y_phi_Weighted_pointer), CompactBinIndex_pT_y_phi);
                hLambdaCounter_pT_y_phi_Weighted_pointer->AddBinContent(GlobalBinIndex_pT_y_phi, hLambdaCounter_pT_y_phi_Weighted_buffer[CompactBinIndex_pT_y_phi]);
                hLambdaAvgDotX_pT_y_phi_Weighted_pointer->AddBinContent(GlobalBinIndex_pT_y_phi, hLambdaAvgDotX_pT_y_phi_Weighted_buffer[CompactBinIndex_pT_y_phi]);
                hLambdaAvgDotY_pT_y_phi_Weighted_pointer->AddBinContent(GlobalBinIndex_pT_y_phi, hLambdaAvgDotY_pT_y_phi_Weighted_buffer[CompactBinIndex_pT_y_phi]);
                hLambdaAvgDotZ_pT_y_phi_Weighted_pointer->AddBinContent(GlobalBinIndex_pT_y_phi, hLambdaAvgDotZ_pT_y_phi_Weighted_buffer[CompactBinIndex_pT_y_phi]);

                hLambdaPolX_pT_y_phi_Weighted_pointer->AddBinContent(GlobalBinIndex_pT_y_phi, hLambdaPolX_pT_y_phi_Weighted_buffer[CompactBinIndex_pT_y_phi]);
                hLambdaPolY_pT_y_phi_Weighted_pointer->AddBinContent(GlobalBinIndex_pT_y_phi, hLambdaPolY_pT_y_phi_Weighted_buffer[CompactBinIndex_pT_y_phi]);
                hLambdaPolZ_pT_y_phi_Weighted_pointer->AddBinContent(GlobalBinIndex_pT_y_phi, hLambdaPolZ_pT_y_phi_Weighted_buffer[CompactBinIndex_pT_y_phi]);
            }

            std::cout<<"DEBUG! Flushing buffers 2 for resample " << resample_idx << std::endl;

                // 2 - Flushing 2D buffers for histograms with pT and y binning:
            std::shared_ptr<TH2D> hLambdaCounter_pT_y_Weighted_pointer = hLambdaCounter_pT_y_Weighted.Get();
            for (int CompactBinIndex_pT_y = 0; CompactBinIndex_pT_y < N_bins_pT*N_bins_rap; CompactBinIndex_pT_y++){
                // Converting the CompactBinIndex back to a global index:
                int GlobalBinIndex_pT_y = mapCompactToGlobal(*(hLambdaCounter_pT_y_Weighted_pointer), CompactBinIndex_pT_y);
                hLambdaCounter_pT_y_Weighted_pointer->AddBinContent(GlobalBinIndex_pT_y, hLambdaCounter_pT_y_Weighted_buffer[CompactBinIndex_pT_y]);
            }

            std::cout<<"DEBUG! Flushing buffers 3 for resample " << resample_idx << std::endl;

                // 3 - Flushing 2D buffers for histograms with phi and pT binning:
            std::shared_ptr<TH2D> hLambdaPolX_phi_pT_Weighted_pointer = hLambdaPolX_phi_pT_Weighted.Get();
            std::shared_ptr<TH2D> hLambdaPolY_phi_pT_Weighted_pointer = hLambdaPolY_phi_pT_Weighted.Get();
            std::shared_ptr<TH2D> hLambdaPolZ_phi_pT_Weighted_pointer = hLambdaPolZ_phi_pT_Weighted.Get();
            std::shared_ptr<TH2D> hLambdaCounter_phi_pT_Weighted_pointer = hLambdaCounter_phi_pT_Weighted.Get();

            for (int CompactBinIndex_phi_pT = 0; CompactBinIndex_phi_pT < N_bins_phi*N_bins_pT; CompactBinIndex_phi_pT++){
                int GlobalBinIndex_phi_pT = mapCompactToGlobal(*(hLambdaPolX_phi_pT_Weighted_pointer), CompactBinIndex_phi_pT);
                hLambdaPolX_phi_pT_Weighted_pointer->AddBinContent(GlobalBinIndex_phi_pT, hLambdaPolX_phi_pT_Weighted_buffer[CompactBinIndex_phi_pT]);
                hLambdaPolY_phi_pT_Weighted_pointer->AddBinContent(GlobalBinIndex_phi_pT, hLambdaPolY_phi_pT_Weighted_buffer[CompactBinIndex_phi_pT]);
                hLambdaPolZ_phi_pT_Weighted_pointer->AddBinContent(GlobalBinIndex_phi_pT, hLambdaPolZ_phi_pT_Weighted_buffer[CompactBinIndex_phi_pT]);
                hLambdaCounter_phi_pT_Weighted_pointer->AddBinContent(GlobalBinIndex_phi_pT, hLambdaCounter_phi_pT_Weighted_buffer[CompactBinIndex_phi_pT]);
            }

            std::cout<<"DEBUG! Flushing buffers 4 for resample " << resample_idx << std::endl;
                // 4 - Flushing 1D buffers of histograms in phi:
            std::shared_ptr<TH1D> hLambdaCounter_phi_Weighted_pointer = hLambdaCounter_phi_Weighted.Get();
            std::shared_ptr<TH1D> hLambdaCounter_phiRingAngles_Weighted_pointer = hLambdaCounter_phiRingAngles_Weighted.Get();
            std::shared_ptr<TH1D> hLambdaPolX_phi_Weighted_pointer = hLambdaPolX_phi_Weighted.Get();
            std::shared_ptr<TH1D> hLambdaPolY_phi_Weighted_pointer = hLambdaPolY_phi_Weighted.Get();
            std::shared_ptr<TH1D> hLambdaPolZ_phi_Weighted_pointer = hLambdaPolZ_phi_Weighted.Get();
            std::shared_ptr<TH1D> hRingObservable_proxy_from_daughter_pointer = hRingObservable_proxy_from_daughter.Get();
            std::shared_ptr<TH1D> hRingObservable_proxy_from_daughter_eq_def_pointer = hRingObservable_proxy_from_daughter_eq_def.Get();
            std::shared_ptr<TH1D> hRingObservable_proxy_from_daughter_star_eq_def_pointer = hRingObservable_proxy_from_daughter_star_eq_def.Get();
            std::shared_ptr<TH1D> hRingObservable_TrueValue_pointer = hRingObservable_TrueValue.Get();
            std::shared_ptr<TH1D> hRingObservable_TrueValuePtCuts_pointer = hRingObservable_TrueValuePtCuts.Get();
            std::shared_ptr<TH1D> hRingObservable_NormPolTrueValuePtCuts_pointer = hRingObservable_NormPolTrueValuePtCuts.Get();
            std::shared_ptr<TH1D> hRingObservable_TrueValuePtCuts_integrated_pointer = hRingObservable_TrueValuePtCuts_integrated.Get();
            std::shared_ptr<TH1D> hRingObservable_NormPolTrueValuePtCuts_integrated_pointer = hRingObservable_NormPolTrueValuePtCuts_integrated.Get();
            std::shared_ptr<TH1D> hRingObservable_proxy_from_daughter_eq_defPtCuts_pointer = hRingObservable_proxy_from_daughter_eq_defPtCuts.Get();
            std::shared_ptr<TH1D> hRingObservable_proxy_from_daughter_eq_def2PtCuts_pointer = hRingObservable_proxy_from_daughter_eq_def2PtCuts.Get();
            std::shared_ptr<TH1D> hRingObservable_proxy_from_daughter_star_eq_defPtCuts_pointer = hRingObservable_proxy_from_daughter_star_eq_defPtCuts.Get();
            std::shared_ptr<TH1D> hRingObservable_proxy_from_daughter_eq_def2PtCuts_integrated_pointer = hRingObservable_proxy_from_daughter_eq_def2PtCuts_integrated.Get();
            std::shared_ptr<TH1D> hRingObservable_proxy_from_daughter_star_eq_defPtCuts_integrated_pointer = hRingObservable_proxy_from_daughter_star_eq_defPtCuts_integrated.Get();

            for (int CompactBinIndex_phi = 0; CompactBinIndex_phi < N_bins_phi; CompactBinIndex_phi++){
                // Converting the CompactBinIndex'es back to global indexes:
                int GlobalBinIndex_phi = mapCompactToGlobal(*(hLambdaCounter_phi_Weighted_pointer), CompactBinIndex_phi);
                int GlobalBinIndex_phiRingAngles = mapCompactToGlobal(*(hLambdaCounter_phiRingAngles_Weighted_pointer), CompactBinIndex_phi);

                hLambdaCounter_phi_Weighted_pointer->AddBinContent(GlobalBinIndex_phi, hLambdaCounter_phi_Weighted_buffer[CompactBinIndex_phi]);
                hLambdaCounter_phiRingAngles_Weighted_pointer->AddBinContent(GlobalBinIndex_phiRingAngles, hLambdaCounter_phiRingAngles_Weighted_buffer[CompactBinIndex_phi]);
                // hLambdaCounter_DeltaphiJRing_Weighted.Get()->AddBinContent(global_bin_idx_phiRing[pre_binned_idx], hLambdaCounter_DeltaphiJRing_Weighted_buffer[pre_binned_idx]);
                
                hLambdaPolX_phi_Weighted_pointer->AddBinContent(GlobalBinIndex_phi, hLambdaPolX_phi_Weighted_buffer[CompactBinIndex_phi]);
                hLambdaPolY_phi_Weighted_pointer->AddBinContent(GlobalBinIndex_phi, hLambdaPolY_phi_Weighted_buffer[CompactBinIndex_phi]);
                hLambdaPolZ_phi_Weighted_pointer->AddBinContent(GlobalBinIndex_phi, hLambdaPolZ_phi_Weighted_buffer[CompactBinIndex_phi]);

                hRingObservable_proxy_from_daughter_pointer->AddBinContent(GlobalBinIndex_phiRingAngles, hRingObservable_proxy_from_daughter_buffer[CompactBinIndex_phi]);
                hRingObservable_proxy_from_daughter_eq_def_pointer->AddBinContent(GlobalBinIndex_phiRingAngles, hRingObservable_proxy_from_daughter_eq_def_buffer[CompactBinIndex_phi]);
                hRingObservable_proxy_from_daughter_star_eq_def_pointer->AddBinContent(GlobalBinIndex_phiRingAngles, hRingObservable_proxy_from_daughter_star_eq_def_buffer[CompactBinIndex_phi]);
                hRingObservable_TrueValue_pointer->AddBinContent(GlobalBinIndex_phiRingAngles, hRingObservable_TrueValue_buffer[CompactBinIndex_phi]);
                hRingObservable_TrueValuePtCuts_pointer->AddBinContent(GlobalBinIndex_phiRingAngles, hRingObservable_TrueValuePtCuts_buffer[CompactBinIndex_phi]);
                hRingObservable_NormPolTrueValuePtCuts_pointer->AddBinContent(GlobalBinIndex_phiRingAngles, hRingObservable_NormPolTrueValuePtCuts_buffer[CompactBinIndex_phi]);
                hRingObservable_TrueValuePtCuts_integrated_pointer->AddBinContent(GlobalBinIndex_phiRingAngles, hRingObservable_TrueValuePtCuts_integrated_buffer[CompactBinIndex_phi]);
                hRingObservable_NormPolTrueValuePtCuts_integrated_pointer->AddBinContent(GlobalBinIndex_phiRingAngles, hRingObservable_NormPolTrueValuePtCuts_integrated_buffer[CompactBinIndex_phi]);
                hRingObservable_proxy_from_daughter_eq_defPtCuts_pointer->AddBinContent(GlobalBinIndex_phiRingAngles, hRingObservable_proxy_from_daughter_eq_defPtCuts_buffer[CompactBinIndex_phi]);
                hRingObservable_proxy_from_daughter_eq_def2PtCuts_pointer->AddBinContent(GlobalBinIndex_phiRingAngles, hRingObservable_proxy_from_daughter_eq_def2PtCuts_buffer[CompactBinIndex_phi]);
                hRingObservable_proxy_from_daughter_star_eq_defPtCuts_pointer->AddBinContent(GlobalBinIndex_phiRingAngles, hRingObservable_proxy_from_daughter_star_eq_defPtCuts_buffer[CompactBinIndex_phi]);
                hRingObservable_proxy_from_daughter_eq_def2PtCuts_integrated_pointer->AddBinContent(GlobalBinIndex_phiRingAngles, hRingObservable_proxy_from_daughter_eq_def2PtCuts_integrated_buffer[CompactBinIndex_phi]);
                hRingObservable_proxy_from_daughter_star_eq_defPtCuts_integrated_pointer->AddBinContent(GlobalBinIndex_phiRingAngles, hRingObservable_proxy_from_daughter_star_eq_defPtCuts_integrated_buffer[CompactBinIndex_phi]);
            }
                // These three did not use buffers:
            // hProtonCounter_phiRing_Weighted.Get()->Flush();
            // hProtonStarCounter_phiRing_Weighted.Get()->Flush();
            // hDebugCounter_phi_star_sampler.Get()->Flush();
        }
    } // End of pragma parallel resampling loop
    std::cout << "\tDone resampling!" << std::endl;

    ////////////////////////////////////////////////////////////////////
    /// Merge histograms from all threads into one before proceeding:
    ////////////////////////////////////////////////////////////////////
    std::shared_ptr<TH3D> hLambdaCounter_pT_y_phi_Weighted_merged = hLambdaCounter_pT_y_phi_Weighted.Merge();
    std::shared_ptr<TH3D> hLambdaAvgDotX_pT_y_phi_Weighted_merged = hLambdaAvgDotX_pT_y_phi_Weighted.Merge();
    std::shared_ptr<TH3D> hLambdaAvgDotY_pT_y_phi_Weighted_merged = hLambdaAvgDotY_pT_y_phi_Weighted.Merge();
    std::shared_ptr<TH3D> hLambdaAvgDotZ_pT_y_phi_Weighted_merged = hLambdaAvgDotZ_pT_y_phi_Weighted.Merge();
    std::shared_ptr<TH2D> hLambdaCounter_pT_y_Weighted_merged = hLambdaCounter_pT_y_Weighted.Merge();

    std::shared_ptr<TH1D> hLambdaCounter_phi_Weighted_merged = hLambdaCounter_phi_Weighted.Merge();
    std::shared_ptr<TH1D> hLambdaCounter_phiRingAngles_Weighted_merged = hLambdaCounter_phiRingAngles_Weighted.Merge();
    std::shared_ptr<TH1D> hLambdaCounter_DeltaphiJRing_Weighted_merged = hLambdaCounter_DeltaphiJRing_Weighted.Merge();
    std::shared_ptr<TH1D> hProtonCounter_phiRing_Weighted_merged = hProtonCounter_phiRing_Weighted.Merge();
    std::shared_ptr<TH1D> hProtonStarCounter_phiRing_Weighted_merged = hProtonStarCounter_phiRing_Weighted.Merge();

    std::shared_ptr<TH2D> hLambdaPolX_phi_pT_Weighted_merged = hLambdaPolX_phi_pT_Weighted.Merge();
    std::shared_ptr<TH2D> hLambdaPolY_phi_pT_Weighted_merged = hLambdaPolY_phi_pT_Weighted.Merge();
    std::shared_ptr<TH2D> hLambdaPolZ_phi_pT_Weighted_merged = hLambdaPolZ_phi_pT_Weighted.Merge();
    std::shared_ptr<TH2D> hLambdaCounter_phi_pT_Weighted_merged = hLambdaCounter_phi_pT_Weighted.Merge();

    std::shared_ptr<TH1D> hDebugCounter_phi_star_sampler_merged = hDebugCounter_phi_star_sampler.Merge();
    std::shared_ptr<TH1D> hLambdaPolX_phi_Weighted_merged = hLambdaPolX_phi_Weighted.Merge();
    std::shared_ptr<TH1D> hLambdaPolY_phi_Weighted_merged = hLambdaPolY_phi_Weighted.Merge();
    std::shared_ptr<TH1D> hLambdaPolZ_phi_Weighted_merged = hLambdaPolZ_phi_Weighted.Merge();

    std::shared_ptr<TH3D> hLambdaPolX_pT_y_phi_Weighted_merged = hLambdaPolX_pT_y_phi_Weighted.Merge();
    std::shared_ptr<TH3D> hLambdaPolY_pT_y_phi_Weighted_merged = hLambdaPolY_pT_y_phi_Weighted.Merge();
    std::shared_ptr<TH3D> hLambdaPolZ_pT_y_phi_Weighted_merged = hLambdaPolZ_pT_y_phi_Weighted.Merge();

    std::shared_ptr<TH1D> hRingObservable_proxy_from_daughter_merged = hRingObservable_proxy_from_daughter.Merge();
    std::shared_ptr<TH1D> hRingObservable_proxy_from_daughter_eq_def_merged = hRingObservable_proxy_from_daughter_eq_def.Merge();
    std::shared_ptr<TH1D> hRingObservable_proxy_from_daughter_star_eq_def_merged = hRingObservable_proxy_from_daughter_star_eq_def.Merge();

    std::shared_ptr<TH1D> hRingObservable_TrueValue_merged = hRingObservable_TrueValue.Merge();

    std::shared_ptr<TH1D> hRingObservable_TrueValuePtCuts_merged = hRingObservable_TrueValuePtCuts.Merge();
    std::shared_ptr<TH1D> hRingObservable_NormPolTrueValuePtCuts_merged = hRingObservable_NormPolTrueValuePtCuts.Merge();

    std::shared_ptr<TH1D> hRingObservable_TrueValuePtCuts_integrated_merged = hRingObservable_TrueValuePtCuts_integrated.Merge();
    std::shared_ptr<TH1D> hRingObservable_NormPolTrueValuePtCuts_integrated_merged = hRingObservable_NormPolTrueValuePtCuts_integrated.Merge();

    std::shared_ptr<TH1D> hRingObservable_proxy_from_daughter_eq_defPtCuts_merged = hRingObservable_proxy_from_daughter_eq_defPtCuts.Merge();
    std::shared_ptr<TH1D> hRingObservable_proxy_from_daughter_eq_def2PtCuts_merged = hRingObservable_proxy_from_daughter_eq_def2PtCuts.Merge();
    std::shared_ptr<TH1D> hRingObservable_proxy_from_daughter_star_eq_defPtCuts_merged = hRingObservable_proxy_from_daughter_star_eq_defPtCuts.Merge();

    std::shared_ptr<TH1D> hRingObservable_proxy_from_daughter_eq_def2PtCuts_integrated_merged = hRingObservable_proxy_from_daughter_eq_def2PtCuts_integrated.Merge();
    std::shared_ptr<TH1D> hRingObservable_proxy_from_daughter_star_eq_defPtCuts_integrated_merged = hRingObservable_proxy_from_daughter_star_eq_defPtCuts_integrated.Merge();


    // Finishing the average for the 3 global dot products -- In 3D:
        // (For weighted averages)
    hLambdaAvgDotX_pT_y_phi_Weighted_merged->Divide(hLambdaCounter_pT_y_phi_Weighted_merged.get()); // Do notice that in the new format, without .get(), we need to pass the address of the TH3D so that ROOT has the expected pointer
    hLambdaAvgDotY_pT_y_phi_Weighted_merged->Divide(hLambdaCounter_pT_y_phi_Weighted_merged.get());
    hLambdaAvgDotZ_pT_y_phi_Weighted_merged->Divide(hLambdaCounter_pT_y_phi_Weighted_merged.get());

    // hLambdaPolStarX_pT_y_phiReco_Weighted = (TH3D*) hLambdaAvgDotX_pT_y_phi_Weighted_merged->Clone("hLambdaPolStarX_pT_y_phiReco_Weighted");
    // hLambdaPolStarY_pT_y_phiReco_Weighted = (TH3D*) hLambdaAvgDotY_pT_y_phi_Weighted_merged->Clone("hLambdaPolStarY_pT_y_phiReco_Weighted");
    // hLambdaPolStarZ_pT_y_phiReco_Weighted = (TH3D*) hLambdaAvgDotZ_pT_y_phi_Weighted_merged->Clone("hLambdaPolStarZ_pT_y_phiReco_Weighted");
    // Cloning everything by hand:
    hLambdaPolStarX_pT_y_phiReco_Weighted = CopyTH3D(hLambdaAvgDotX_pT_y_phi_Weighted_merged.get());
    hLambdaPolStarY_pT_y_phiReco_Weighted = CopyTH3D(hLambdaAvgDotY_pT_y_phi_Weighted_merged.get());
    hLambdaPolStarZ_pT_y_phiReco_Weighted = CopyTH3D(hLambdaAvgDotZ_pT_y_phi_Weighted_merged.get());

    hLambdaPolStarX_pT_y_phiReco_Weighted->SetTitle("hLambdaPolStarX_pT_y_phiReco_Weighted");
    hLambdaPolStarY_pT_y_phiReco_Weighted->SetTitle("hLambdaPolStarY_pT_y_phiReco_Weighted");
    hLambdaPolStarZ_pT_y_phiReco_Weighted->SetTitle("hLambdaPolStarZ_pT_y_phiReco_Weighted");
    
        // Applying the transform that turns them into polarizations -- In the lambda rest frame!
    hLambdaPolStarX_pT_y_phiReco_Weighted->Scale(3.0 / alpha_H);
    hLambdaPolStarY_pT_y_phiReco_Weighted->Scale(3.0 / alpha_H);
    hLambdaPolStarZ_pT_y_phiReco_Weighted->Scale(3.0 / alpha_H);
    ///////////////////////////////////////////////////////////////

        // Now boosting each of these into the lab rest frame by calculating a mean 4-vector of all Lambdas in each (pT, y) bin:
    for (int pT_idx = 1; pT_idx <= N_bins_pT; pT_idx++){ // bin index starts at 1 (not 0)
        for (int y_idx = 1; y_idx <= N_bins_rap; y_idx++){
            for (int phi_idx = 1; phi_idx <= N_bins_phi; phi_idx++){
                // Getting all the needed values for all three spatial components:
                double pT_bin_center = hLambdaCounter_pT_y_phi_Weighted_merged->GetXaxis()->GetBinCenter(pT_idx);
                double y_bin_center = hLambdaCounter_pT_y_phi_Weighted_merged->GetYaxis()->GetBinCenter(y_idx);
                double phi_bin_center = hLambdaCounter_pT_y_phi_Weighted_merged->GetZaxis()->GetBinCenter(phi_idx);

                // std::cout << "pT bin " << pT_idx << ", center (" << pT_bin_center << "), y: " << y_bin_center << " phi: " << phi_bin_center << std::endl;

                    // For multiplicity-weighted averages of the polarization:
                double PolX_star_pTy_phi_Weighted = hLambdaPolStarX_pT_y_phiReco_Weighted->GetBinContent(pT_idx, y_idx, phi_idx);
                double PolY_star_pTy_phi_Weighted = hLambdaPolStarY_pT_y_phiReco_Weighted->GetBinContent(pT_idx, y_idx, phi_idx);
                double PolZ_star_pTy_phi_Weighted = hLambdaPolStarZ_pT_y_phiReco_Weighted->GetBinContent(pT_idx, y_idx, phi_idx);

                // Now calculating a mean 4-momentum of the Lambda in each bin:
                double mean_px = pT_bin_center * std::cos(phi_bin_center);
                double mean_py = pT_bin_center * std::sin(phi_bin_center);
                // double mean_mT = std::sqrt(pT_bin_center * pT_bin_center + mL2);
                // double mean_pz = mean_mT * std::sinh(y_bin_center);
                // double mean_E  = mean_mT * std::cosh(y_bin_center);
                    // Actually using a different definition of pz and E, as the provided variables from iSS seem to be (pT, y, phi),
                    // not (px, py, pz) which are converted into those three variables. Thus, as y is not derived from tanh^-1 (pz/E),
                    // and is an independently generated quantity, it is better to use the following:
                
                    // Actually, forget these definitions: I will use what Vitor used:
                // double mean_pz = pT_bin_center * std::sinh(y_bin_center);
                // double mean_E  = std::sqrt(pT_bin_center*pT_bin_center * std::cosh(y_bin_center)*std::cosh(y_bin_center) + mL2);
                double mean_mT = std::sqrt(mL2 + pT_bin_center * pT_bin_center);
                double mean_pz = mean_mT * std::sinh(y_bin_center);
                double mean_E = std::sqrt(mL2 + mean_px * mean_px + mean_py * mean_py + mean_pz * mean_pz);

                TLorentzVector mean_Lambda_4vec_lab(mean_px, mean_py, mean_pz, mean_E);
                    // Building boost vector for the polarization in this bin:
                TVector3 beta = mean_Lambda_4vec_lab.BoostVector(); // The boost that takes from the Lambda rest frame into the lab frame.

                // Building the mean polarization 4-momentum, with the fact that S^0 = 0 in the particle's own rest frame (no temporal component of spin in the rest frame!):
                TLorentzVector P_Lambda_star_mean_4vec_Weighted(PolX_star_pTy_phi_Weighted, PolY_star_pTy_phi_Weighted, PolZ_star_pTy_phi_Weighted, 0);
                TLorentzVector P_Lambda_lab_mean_4vec_Weighted = P_Lambda_star_mean_4vec_Weighted;
                P_Lambda_lab_mean_4vec_Weighted.Boost(beta);

                ///////////////////////////////////////////
                // // Making sure that the S^mu * P_mu product is zero (as it should, because S^0 = 0 and the trimomentum of P_mu is zero in the Lambda rest frame).
                // // Compute scalar product S.P
                // double dotProductReco = P_Lambda_star_mean_4vec_Weighted * mean_Lambda_4vec_lab;
                //
                // // Report
                // std::cout << "Reconstructed pol: 2*S^mu * P_mu = " << dotProductReco << std::endl;
                // It is also zero within machine precision!
                // ///////////////////////////////////////////

                /////////////////////////////////////////////////////
                // Checking if the boost is working and the polarization is properly being converted to the lab frame:
                // std::cout << "Before (PolStar): " << "(PolX, PolY, PolZ, PolT) = (" << P_Lambda_star_mean_4vec_Weighted.X() << ", " << P_Lambda_star_mean_4vec_Weighted.Y() << ", " << P_Lambda_star_mean_4vec_Weighted.Z() << ", " << P_Lambda_star_mean_4vec_Weighted.T() << ")" << std::endl;
                // std::cout << "After (PolLab):" << "(PolX, PolY, PolZ, PolT) = (" << P_Lambda_lab_mean_4vec_Weighted.X() << ", " << P_Lambda_lab_mean_4vec_Weighted.Y() << ", " << P_Lambda_lab_mean_4vec_Weighted.Z() << ", " << P_Lambda_lab_mean_4vec_Weighted.T() << ")" << std::endl;
                /////////////////////////////////////////////////////

                // For the weighted variation of the polarization:
                hLambdaPolX_pT_y_phiReco_Weighted->Fill(pT_bin_center, y_bin_center, phi_bin_center, P_Lambda_lab_mean_4vec_Weighted.X());
                hLambdaPolY_pT_y_phiReco_Weighted->Fill(pT_bin_center, y_bin_center, phi_bin_center, P_Lambda_lab_mean_4vec_Weighted.Y());
                hLambdaPolZ_pT_y_phiReco_Weighted->Fill(pT_bin_center, y_bin_center, phi_bin_center, P_Lambda_lab_mean_4vec_Weighted.Z());

                // Also doing a projection in 2D over pT and y:
                    // Correction -- When reducing dimensionality, you need to ponder the polarization by the multiplicity of the bin it belongs to!
                    // (summing over all bins as if they were never there requires you to sum over particles and their appropriate weights):
                double current_pTyphi_bin_multiplicity = hLambdaCounter_pT_y_phi_Weighted_merged->GetBinContent(pT_idx, y_idx, phi_idx);
                hLambdaPolX_pT_yReco_Weighted->Fill(pT_bin_center, y_bin_center, P_Lambda_lab_mean_4vec_Weighted.X() * current_pTyphi_bin_multiplicity); // Don't need to average over N_phi bins or anything like that!
                hLambdaPolY_pT_yReco_Weighted->Fill(pT_bin_center, y_bin_center, P_Lambda_lab_mean_4vec_Weighted.Y() * current_pTyphi_bin_multiplicity);
                hLambdaPolZ_pT_yReco_Weighted->Fill(pT_bin_center, y_bin_center, P_Lambda_lab_mean_4vec_Weighted.Z() * current_pTyphi_bin_multiplicity);
                
                    // A second 2D projection, in the same way it was plotted on Vitor's paper:
                hLambdaPolX_phi_pTReco_Weighted->Fill(phi_bin_center, pT_bin_center, P_Lambda_lab_mean_4vec_Weighted.X() * current_pTyphi_bin_multiplicity); // Also with the multiplicity correction
                hLambdaPolY_phi_pTReco_Weighted->Fill(phi_bin_center, pT_bin_center, P_Lambda_lab_mean_4vec_Weighted.Y() * current_pTyphi_bin_multiplicity);
                hLambdaPolZ_phi_pTReco_Weighted->Fill(phi_bin_center, pT_bin_center, P_Lambda_lab_mean_4vec_Weighted.Z() * current_pTyphi_bin_multiplicity);

                    // Finally, a 1D projection for the distribution in Phi only:
                hLambdaPolX_phiReco_Weighted->Fill(phi_bin_center, P_Lambda_lab_mean_4vec_Weighted.X() * current_pTyphi_bin_multiplicity); // Also with the multiplicity correction
                hLambdaPolY_phiReco_Weighted->Fill(phi_bin_center, P_Lambda_lab_mean_4vec_Weighted.Y() * current_pTyphi_bin_multiplicity);
                hLambdaPolZ_phiReco_Weighted->Fill(phi_bin_center, P_Lambda_lab_mean_4vec_Weighted.Z() * current_pTyphi_bin_multiplicity);

                    // Experimenting with 2D projections for the true values to see if this process would work correctly:
                // hLambdaPolX_phi_pT_Weighted_ProjTEST->Fill(phi_bin_center, pT_bin_center, hLambdaPolX_pT_y_phi_Weighted_merged.GetBinContent(pT_idx, y_idx, phi_idx) * current_pTyphi_bin_multiplicity);
                // The histograms where not yet divided by the number of counts in each bin, so you don't need to multiply by the current bins' multiplicity again!
                hLambdaPolX_phi_pT_Weighted_ProjTEST->Fill(phi_bin_center, pT_bin_center, hLambdaPolX_pT_y_phi_Weighted_merged->GetBinContent(pT_idx, y_idx, phi_idx));
                hLambdaPolY_phi_pT_Weighted_ProjTEST->Fill(phi_bin_center, pT_bin_center, hLambdaPolY_pT_y_phi_Weighted_merged->GetBinContent(pT_idx, y_idx, phi_idx));
                hLambdaPolZ_phi_pT_Weighted_ProjTEST->Fill(phi_bin_center, pT_bin_center, hLambdaPolZ_pT_y_phi_Weighted_merged->GetBinContent(pT_idx, y_idx, phi_idx));

                // Calculating the ring observable:
                    // Cross product (trigger X p)
                    const double x_trigger = 1.; // Defined coordinate system so that this is always true, in a way that I can average over many events!
                    const double y_trigger = 0.;
                    const double z_trigger = 0.;
                    double cross_x = y_trigger*mean_pz - z_trigger*mean_py;
                    double cross_y = z_trigger*mean_px - x_trigger*mean_pz;
                    double cross_z = x_trigger*mean_py - y_trigger*mean_px;

                    double cross_product_norm = std::sqrt(cross_x*cross_x + cross_y*cross_y + cross_z*cross_z);
		    			
                    double RP_temp = ((P_Lambda_lab_mean_4vec_Weighted.X()*cross_x + P_Lambda_lab_mean_4vec_Weighted.Y()*cross_y + P_Lambda_lab_mean_4vec_Weighted.Z()*cross_z)/
            					      cross_product_norm);

                    double delta_phi_J = phi_bin_center; // In this coordinate definition after the phi_matrix rotation, the jet is always facing the x axis and phi is already a delta_phi_J
                    delta_phi_J = wrapToInterval(delta_phi_J, phi_min_Ring, phi_max_Ring);

                    hRingObservableReco->Fill(delta_phi_J, RP_temp * current_pTyphi_bin_multiplicity);

                        // An observable that carries only the direction of the polarization, not the magnitude:
                    double Pol_norm = std::sqrt(P_Lambda_lab_mean_4vec_Weighted.X()*P_Lambda_lab_mean_4vec_Weighted.X()
                                                + P_Lambda_lab_mean_4vec_Weighted.Y()*P_Lambda_lab_mean_4vec_Weighted.Y()
                                                + P_Lambda_lab_mean_4vec_Weighted.Z()*P_Lambda_lab_mean_4vec_Weighted.Z());
                    double RP_temp_Norm = ((P_Lambda_lab_mean_4vec_Weighted.X()*cross_x + P_Lambda_lab_mean_4vec_Weighted.Y()*cross_y 
                                            + P_Lambda_lab_mean_4vec_Weighted.Z()*cross_z)/cross_product_norm)/Pol_norm;

                    if (pT_bin_center > 0.5 && pT_bin_center < 1.5){
                        hRingObservableRecoPtCuts->Fill(delta_phi_J, RP_temp * current_pTyphi_bin_multiplicity);
                        hRingObservableNormPolRecoPtCuts->Fill(delta_phi_J, RP_temp_Norm * current_pTyphi_bin_multiplicity);

                        hRingObservableRecoPtCuts_integrated->Fill(delta_phi_J, RP_temp * current_pTyphi_bin_multiplicity);
                        hRingObservableNormPolRecoPtCuts_integrated->Fill(delta_phi_J, RP_temp_Norm * current_pTyphi_bin_multiplicity);
                    }
            }
        }
    }

    // Averaging the true values too -- This should already mask the arrays too:
        // Peace of mind plots averaging:
    hLambdaPolX_phi_pT_Weighted_merged->Divide(hLambdaCounter_phi_pT_Weighted_merged.get());
    hLambdaPolY_phi_pT_Weighted_merged->Divide(hLambdaCounter_phi_pT_Weighted_merged.get());
    hLambdaPolZ_phi_pT_Weighted_merged->Divide(hLambdaCounter_phi_pT_Weighted_merged.get());

        // For the second type:
    hLambdaPolX_pT_y_phi_Weighted_merged->Divide(hLambdaCounter_pT_y_phi_Weighted_merged.get());
    hLambdaPolY_pT_y_phi_Weighted_merged->Divide(hLambdaCounter_pT_y_phi_Weighted_merged.get());
    hLambdaPolZ_pT_y_phi_Weighted_merged->Divide(hLambdaCounter_pT_y_phi_Weighted_merged.get());

            // For its projection test (look at the procedure employed in hLambdaPolX_phi_pTReco_Weighted):
            // Compare these with hLambdaPolX_phi_pT_Weighted_merged
    hLambdaPolX_phi_pT_Weighted_ProjTEST->Divide(hLambdaCounter_phi_pT_Weighted_merged.get());
    hLambdaPolY_phi_pT_Weighted_ProjTEST->Divide(hLambdaCounter_phi_pT_Weighted_merged.get());
    hLambdaPolZ_phi_pT_Weighted_ProjTEST->Divide(hLambdaCounter_phi_pT_Weighted_merged.get());

    // Dividing by the total number of particles to properly do a weighted average using the current_pTyphi_bin_multiplicity values:
    hLambdaPolX_pT_yReco_Weighted->Divide(hLambdaCounter_pT_y_Weighted_merged.get());
    hLambdaPolY_pT_yReco_Weighted->Divide(hLambdaCounter_pT_y_Weighted_merged.get());
    hLambdaPolZ_pT_yReco_Weighted->Divide(hLambdaCounter_pT_y_Weighted_merged.get());

    hLambdaPolX_phi_pTReco_Weighted->Divide(hLambdaCounter_phi_pT_Weighted_merged.get());
    hLambdaPolY_phi_pTReco_Weighted->Divide(hLambdaCounter_phi_pT_Weighted_merged.get());
    hLambdaPolZ_phi_pTReco_Weighted->Divide(hLambdaCounter_phi_pT_Weighted_merged.get());

    hLambdaPolX_phiReco_Weighted->Divide(hLambdaCounter_phi_Weighted_merged.get());
    hLambdaPolY_phiReco_Weighted->Divide(hLambdaCounter_phi_Weighted_merged.get());
    hLambdaPolZ_phiReco_Weighted->Divide(hLambdaCounter_phi_Weighted_merged.get());

        // The same for the expected true values of the phi-only polarization histogram (to compare with figure 3 of the paper):
    hLambdaPolX_phi_Weighted_merged->Divide(hLambdaCounter_phi_Weighted_merged.get());
    hLambdaPolY_phi_Weighted_merged->Divide(hLambdaCounter_phi_Weighted_merged.get());
    hLambdaPolZ_phi_Weighted_merged->Divide(hLambdaCounter_phi_Weighted_merged.get());

    // Normalizing all ring observables by the particle weights to complete the averaging:
    hRingObservable_TrueValue_merged->Divide(hLambdaCounter_DeltaphiJRing_Weighted_merged.get()); // Notice I'm using hLambdaCounter_DeltaphiJRing_Weighted_merged, not hLambdaCounter_phiRingAngles_Weighted_merged, to have the proper counts in DeltaPhiJ coordinates!
    hRingObservable_TrueValuePtCuts_merged->Divide(hLambdaCounter_DeltaphiJRing_Weighted_merged.get());
    hRingObservable_NormPolTrueValuePtCuts_merged->Divide(hLambdaCounter_DeltaphiJRing_Weighted_merged.get());
    hRingObservable_TrueValuePtCuts_integrated_merged->Scale(1./hLambdaCounter_DeltaphiJRing_Weighted_merged->Integral());
    hRingObservable_NormPolTrueValuePtCuts_integrated_merged->Scale(1./hLambdaCounter_DeltaphiJRing_Weighted_merged->Integral());

    hRingObservableReco->Divide(hLambdaCounter_DeltaphiJRing_Weighted_merged.get()); // These three also need to be passed as references
    hRingObservableRecoPtCuts->Divide(hLambdaCounter_DeltaphiJRing_Weighted_merged.get());
    hRingObservableNormPolRecoPtCuts->Divide(hLambdaCounter_DeltaphiJRing_Weighted_merged.get());

    hRingObservableRecoPtCuts_integrated->Scale(1./hLambdaCounter_DeltaphiJRing_Weighted_merged->Integral());
    hRingObservableNormPolRecoPtCuts_integrated->Scale(1./hLambdaCounter_DeltaphiJRing_Weighted_merged->Integral());
        // Should these following proxies be normalized by the number of protons in each phi angle, as the weighted averages somewhat have that in mind?
        // Probably not, because I fill the histograms with the Lambda's phi angle, but it would be interesting to have a version filled with the proton's
        // angle and pondering the ring observable with the proton multiplicity on each bin.
        // todo: test an averaging weight based on the proton weight on each of its own phi bins maybe? Something like an alternative ring observable.
    hRingObservable_proxy_from_daughter_merged->Divide(hLambdaCounter_DeltaphiJRing_Weighted_merged.get());
    hRingObservable_proxy_from_daughter_eq_def_merged->Divide(hLambdaCounter_DeltaphiJRing_Weighted_merged.get());
    hRingObservable_proxy_from_daughter_eq_defPtCuts_merged->Divide(hLambdaCounter_DeltaphiJRing_Weighted_merged.get());
    hRingObservable_proxy_from_daughter_eq_def2PtCuts_merged->Divide(hLambdaCounter_DeltaphiJRing_Weighted_merged.get());

    hRingObservable_proxy_from_daughter_eq_def2PtCuts_integrated_merged->Scale(1./hLambdaCounter_DeltaphiJRing_Weighted_merged->Integral());
    hRingObservable_proxy_from_daughter_star_eq_defPtCuts_integrated_merged->Scale(1./hLambdaCounter_DeltaphiJRing_Weighted_merged->Integral());

        // Another attempt (with proton momentum in the rest frame):
    hRingObservable_proxy_from_daughter_star_eq_def_merged->Divide(hLambdaCounter_DeltaphiJRing_Weighted_merged.get());
    hRingObservable_proxy_from_daughter_star_eq_defPtCuts_merged->Divide(hLambdaCounter_DeltaphiJRing_Weighted_merged.get());

    ////////////////////
    //// 4 - Exporting into a .root file for later review
    ////////////////////
    std::cout << "\nExporting results" << std::endl;
    // std::string filename = "/home/cicero/results/hydro_vorticity/LambdaPol_lumpy_events_40_50_" + std::to_string(N_events) + "ev.root";
    // std::string filename = "/home/cicero/results/hydro_vorticity/LambdaPol2_lumpy_events_40_50_" + std::to_string(N_events) + "ev.root"; // Having many different file names for comparison
    std::string filename = "/home/cicero/results/hydro_vorticity/LambdaPolResampled" + std::to_string(N_resamples) + "times_lumpy_events_40_50_" 
                           + ((with_bullet) ? "with_bullet_" : "no_bullet_") + std::to_string(N_events) + "ev"; // Used a ternary operator to symplify syntax, but beware of the proper encapsulation!
    if (force_DecayDist_mode == 1){ // 0 for the usual, physical, decay. 1 for decays along Pol. 2 for decays perpendicular to Pol
        filename += "_ForcedDecayAlongP";
    }
    else if (force_DecayDist_mode == 2){
        filename += "_ForcedDecayPerpP";
    }
    filename+= ".root";

    TFile f(filename.c_str(), "RECREATE");

        // Saving some particle counters
    hLambdaCounter->Write(); // An unweighted counter, just to know how many bins we have went through
    hLambdaCounter_pT_y_phi_Weighted_merged->Write();
    hLambdaCounter_phi_pT_Weighted_merged->Write();
    hLambdaCounter_pT_y_Weighted_merged->Write();

        // Some more counters, for the ring observables part:
    hLambdaCounter_phi_Weighted_merged->Write();
    hLambdaCounter_phiRingAngles_Weighted_merged->Write();
    hLambdaCounter_DeltaphiJRing_Weighted_merged->Write();
    hProtonCounter_phiRing_Weighted_merged->Write();
    hProtonStarCounter_phiRing_Weighted_merged->Write();
    hDebugCounter_phi_star_sampler_merged->Write();

        // Peace-of-mind plots -- True values!
    hLambdaPolX_phi_pT_Weighted_merged->Write();
    hLambdaPolY_phi_pT_Weighted_merged->Write();
    hLambdaPolZ_phi_pT_Weighted_merged->Write();

        // Some 3D versions of those, for comparison:
    hLambdaPolX_pT_y_phi_Weighted_merged->Write();
    hLambdaPolY_pT_y_phi_Weighted_merged->Write();
    hLambdaPolZ_pT_y_phi_Weighted_merged->Write();

        // Their 2D projection tests:
        // To compare with "hLambdaPolX_phi_pT_Weighted_merged" histograms
    hLambdaPolX_phi_pT_Weighted_ProjTEST->Write();
    hLambdaPolY_phi_pT_Weighted_ProjTEST->Write();
    hLambdaPolZ_phi_pT_Weighted_ProjTEST->Write();

    // For the phi-only averages and true values:
    hLambdaPolX_phi_Weighted_merged->Write();
    hLambdaPolY_phi_Weighted_merged->Write();
    hLambdaPolZ_phi_Weighted_merged->Write();

    // Weighted averages -- Reconstructed values:
    hLambdaAvgDotX_pT_y_phi_Weighted_merged->Write();
    hLambdaAvgDotY_pT_y_phi_Weighted_merged->Write();
    hLambdaAvgDotZ_pT_y_phi_Weighted_merged->Write();

    hLambdaPolStarX_pT_y_phiReco_Weighted->Write();
    hLambdaPolStarY_pT_y_phiReco_Weighted->Write();
    hLambdaPolStarZ_pT_y_phiReco_Weighted->Write();

    hLambdaPolX_pT_y_phiReco_Weighted->Write();
    hLambdaPolY_pT_y_phiReco_Weighted->Write();
    hLambdaPolZ_pT_y_phiReco_Weighted->Write();

    hLambdaPolX_pT_yReco_Weighted->Write();
    hLambdaPolY_pT_yReco_Weighted->Write();
    hLambdaPolZ_pT_yReco_Weighted->Write();

    hLambdaPolX_phi_pTReco_Weighted->Write();
    hLambdaPolY_phi_pTReco_Weighted->Write();
    hLambdaPolZ_phi_pTReco_Weighted->Write();

    hLambdaPolX_phiReco_Weighted->Write();
    hLambdaPolY_phiReco_Weighted->Write();
    hLambdaPolZ_phiReco_Weighted->Write();

    // For the ring observables:
    hRingObservable_TrueValue_merged->Write();
    hRingObservable_TrueValuePtCuts_merged->Write();
    hRingObservable_NormPolTrueValuePtCuts_merged->Write();
    
    hRingObservableReco->Write();
    hRingObservableRecoPtCuts->Write();
    hRingObservableNormPolRecoPtCuts->Write();

    hRingObservable_proxy_from_daughter_merged->Write();
    hRingObservable_proxy_from_daughter_eq_def_merged->Write();
    hRingObservable_proxy_from_daughter_eq_defPtCuts_merged->Write();
    hRingObservable_proxy_from_daughter_eq_def2PtCuts_merged->Write();

    hRingObservable_proxy_from_daughter_star_eq_def_merged->Write();
    hRingObservable_proxy_from_daughter_star_eq_defPtCuts_merged->Write();

        // Integrated observables:
    hRingObservable_TrueValuePtCuts_integrated_merged->Write();
    hRingObservable_NormPolTrueValuePtCuts_integrated_merged->Write();
    hRingObservableRecoPtCuts_integrated->Write();
    hRingObservableNormPolRecoPtCuts_integrated->Write();
    hRingObservable_proxy_from_daughter_eq_def2PtCuts_integrated_merged->Write();
    hRingObservable_proxy_from_daughter_star_eq_defPtCuts_integrated_merged->Write();

    // Deleting the Clone() histograms, which may still be kept in memory:
    // delete hLambdaPolStarX_pT_y_phiReco_Weighted;
    // delete hLambdaPolStarY_pT_y_phiReco_Weighted;
    // delete hLambdaPolStarZ_pT_y_phiReco_Weighted;

    f.Close();
    std::cout << "Done exporting!" << std::endl;

    auto end_time = std::chrono::high_resolution_clock::now();
        // Calculate elapsed time in seconds:
    std::chrono::duration<double> total_elapsed = end_time - start_time;
    std::chrono::duration<double> resampling_elapsed = end_time - start_time - duration_fetch - pre_process_duration; // Removes the time spent on fetching data

    std::cout << "\n////// TIME STATISTICS //////" << std::endl;
    std::cout << "The program ran for: " << total_elapsed.count() << " seconds." << std::endl;
    std::cout << "Of those, " << duration_fetch.count() << " seconds were used to fetch data, " << pre_process_duration.count() << " seconds to pre-process for the parallel loop, and " << resampling_elapsed.count() << " seconds to sample distributions." << std::endl;
    std::cout << "The code ran with " << N_threads << " parallel threads, through " << N_events << " events, which were resampled " << N_resamples << " times." << std::endl;
    std::cout << "The average time for each sample was " << resampling_elapsed.count() * 1./N_resamples << " seconds per sample, or " << resampling_elapsed.count() * 1./N_threads << " seconds of total time per thread." << std::endl;
    // std::cout << "Time per sample, per thread, is: " << resampling_elapsed.count() * 1./(N_resamples * N_threads) << " seconds." << std::endl;
    std::cout << "/////////////////////////////" << std::endl;

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

// Older, inefficient version:
// double wrapToInterval(double phi, double phi_min, double phi_max){
//     double range = phi_max - phi_min;
//     while (phi < phi_min) phi += range;
//     while (phi >= phi_max) phi -= range;
//     return phi;
// }

/// Checks if there are any points of the DoubleMatrices with pT, y or phi that would fall into the underflow or overflow bins of the histograms
/// (A sanity check for the buffer implementation)
void CheckBins(const char* name, const int* bins, int N_bins_total, TH1* hist){ // Accepts TH1 because TH2D and TH3D internally run on TH1
    int nbins_total = hist->GetNcells(); // includes under/overflow

    for (int i = 0; i < N_bins_total; i++) {
        int bin = bins[i];
        if (bin < 1 || bin > nbins_total - 2) { // 0 = underflow, nbins_total-1 = overflow
            std::cout << "[!] Out-of-range bin in " << name 
                    << " at index " << i 
                    << " → globalBin = " << bin 
                    << " (valid range: 1.." << nbins_total-2 << ")"
                    << std::endl;
        }
    }
    std::cout << "Check finished for " << name << std::endl;
}

/// Map ROOT global bin index to a compact index (0 … N_bins_total-1).
/// Works for TH1D, TH2D, TH3D.
/// Throws/asserts if the bin is under/overflow.
inline int mapGlobalToCompact(const TH1& h, int globalBin){ // Uses const TH1& h instead of const TH1* h because the TThreadObject returns std::shared_ptr<TH1D>, not TH1D* values
    int ix, iy, iz;
    h.GetBinXYZ(globalBin, ix, iy, iz);

    int nx = h.GetNbinsX();
    int ny = h.GetNbinsY();
    int nz = h.GetNbinsZ();

    // Check if it's inside the physical bin ranges -- If not, will crash the program
    if (ix < 1 || ix > nx){
        std::cout << "Failed on bin ix = " << ix << " iy = " << iy << " iz = " << iz << std::endl;
        throw std::out_of_range("mapGlobalToCompact: X bin under/overflow");
    }
    if (h.GetDimension() >= 2 && (iy < 1 || iy > ny)){
        std::cout << "Failed on bin ix = " << ix << " iy = " << iy << " iz = " << iz << std::endl;
        throw std::out_of_range("mapGlobalToCompact: Y bin under/overflow");
    }
    if (h.GetDimension() == 3 && (iz < 1 || iz > nz)){
        std::cout << "Failed on bin ix = " << ix << " iy = " << iy << " iz = " << iz << std::endl;
        throw std::out_of_range("mapGlobalToCompact: Z bin under/overflow");
    }

    // Convert to 0-based
    int cx = ix - 1;
    int cy = (h.GetDimension() >= 2) ? iy - 1 : 0;
    int cz = (h.GetDimension() == 3) ? iz - 1 : 0;

    // Flatten to linear index
    return cx + nx * (cy + ny * cz);
}

/// The inverse transform from mapGlobalToCompact
inline int mapCompactToGlobal(const TH1& h, int compactIndex){
    int nx = h.GetNbinsX();
    int ny = h.GetNbinsY();
    int nz = h.GetNbinsZ();

    int cx, cy, cz;
    if (h.GetDimension() == 3){
        cz = compactIndex / (nx * ny);
        int rem = compactIndex % (nx * ny);
        cy = rem / nx;
        cx = rem % nx;
    }
    else if (h.GetDimension() == 2){
        cz = 0; // dummy dimension
        cy = compactIndex / nx;
        cx = compactIndex % nx;
    }
    else if (h.GetDimension() == 1){
        cz = 0; cy = 0;
        cx = compactIndex;
    }
    else{
        throw std::logic_error("mapCompactToGlobal: unsupported histogram dimension");
    }

    // Convert back to ROOT bin numbers (1-based)
    int ix = cx + 1;
    int iy = (h.GetDimension() >= 2) ? cy + 1 : 0;
    int iz = (h.GetDimension() == 3) ? cz + 1 : 0;

    // Optional safety: assert range
    // assert(compactIndex >= 0 && compactIndex < nx * std::max(1, ny) * std::max(1, nz));

    // ROOT’s global bin
    return h.GetBin(ix, iy, iz);  // works for TH1D, TH2D, TH3D
}

// Copy contents of a TH3D into a new TH3D -- Better than the Clone() function...
TH3D* CopyTH3D(const TH3D* h){
    if (!h) return nullptr;

    // Create a new histogram with the same binning and axis titles
    TH3D* hcopy = new TH3D(
        Form("%s_copy", h->GetName()),
        h->GetTitle(),
        h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(),
        h->GetNbinsY(), h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax(),
        h->GetNbinsZ(), h->GetZaxis()->GetXmin(), h->GetZaxis()->GetXmax()
    );

    // Copy bin contents and errors
    for (int ix = 0; ix <= h->GetNbinsX() + 1; ++ix) {
        for (int iy = 0; iy <= h->GetNbinsY() + 1; ++iy) {
            for (int iz = 0; iz <= h->GetNbinsZ() + 1; ++iz) {
                int bin = h->GetBin(ix, iy, iz);
                hcopy->SetBinContent(bin, h->GetBinContent(bin));
                hcopy->SetBinError(bin, h->GetBinError(bin));
            }
        }
    }

    return hcopy;
}


///////////////////////////////////
/// Defining sampling functions ///
///////////////////////////////////
// Declared the function as an inline-able object, so that sampling can take as short a time as possible:
inline double P_angle_proton(double xi_star, double P_Lambda_star){
    return (1.0 / PI) * (1.0 + alpha_H * P_Lambda_star * std::cos(xi_star));
}

// Not quite necessary as a function (because I used its CDF, not the PDF by itself), but kept as a reference for sample_P_angle_proton_from_cos_xi
inline double P_angle_proton_as_cos_func(double cos_xi_star, double P_Lambda_star){
    return (1.0 / 2.0) * (1.0 + alpha_H * P_Lambda_star * cos_xi_star);
}

// Now the actual sampler function -- The goal is to sample an angle given a P(angle) function.
    // Will use a rejection sampling for quick samples!
    // Receives the Lambda polarization magnitude, already in the rest frame, and returns the two angles needed 
    // Usage example:
    // auto [xi_star, phi_star] = sample_P_angle_proton_for_each_PLambda(...);
std::pair<double, double> sample_P_angle_proton(double P_Lambda_star_mag, std::mt19937 &rng, std::uniform_real_distribution<double> dist_x,
                                                std::uniform_real_distribution<double> dist_y, std::uniform_real_distribution<double> &dist_azimuth){
    // First sampling the azimuthal angle, which is really simple to do:
    double phi_star = dist_azimuth(rng);

    // Rejection sampling loop - Generates one sample for the current P_Lambda_star_mag:
    double xi_star = 0;
    while (true){
        // Propose an x uniformly in [0, pi]
        double xi_star_candidate = dist_x(rng);

        // Propose a y uniformly in [0, P_max]
        double y = dist_y(rng);  // Draws from the maximum value distribution independently.
                                    // This makes it so that we always follow the appropriate limit and have correctly weighted samples of x(!)

        // Compute the actual probability density at x
        double px = P_angle_proton(xi_star_candidate, P_Lambda_star_mag);

        // Accept x if y < P(x). This ensures samples are distributed according to P(x): larger P(x) values will be accepted more easily!
        if (y < px){
            xi_star = xi_star_candidate;  // Store accepted xi_star_candidate
            break;                        // Done with this P_Lambda_star_mag. Breaks the loop
        }
        // Otherwise, reject and try again
    }
    return std::make_pair(xi_star, phi_star);
}

// A simpler implementation of sample_P_angle_proton, which samples a cos(xi_star) instead of xi_star directly.
// As cosine is invertible in the (0, PI) interval we are analyzing, then we can just invert cos(xi_star) into xi_star later,
// and the sampling procedure is far simpler for this case where there is no rejection sampling!
// SEE LOG 471, PAGE 4 !!!
std::pair<double, double> sample_P_angle_proton_from_cos_xi(double P_Lambda_star_mag, std::mt19937 &rng, std::uniform_real_distribution<double> &dist_unit,
                                                               std::uniform_real_distribution<double> &dist_azimuth){ // Should pass rng as a reference! If not, it won't update and you will always sample the same values!
    // First sampling the azimuthal angle, which is really simple to do:
    double phi_star = dist_azimuth(rng);
    double u = dist_unit(rng); // The uniform unit sample
    double A = P_Lambda_star_mag * alpha_H;
    // double A = alpha_H; // Testing to see if Mike's distribution was correct (doesn't seem like it!)

    // Calculate the sample using the inverse CDF of P_angle_proton_as_cos_func:
    double cos_xi_star;
    if (P_Lambda_star_mag < 1e-9){ // If polarization is close to zero, should use the other formula for the resampling
        // Special case for A ~ 0
        // F_inv(u) = 2u - 1
        cos_xi_star = 2.0 * u - 1.0;
    }
    else{
        // General case using the quadratic formula solution
        // F_inv(u) = (-1 + sqrt((1-A)^2 + 4Au)) / A // Took the positive root to define xi_star in [0, pi]
        cos_xi_star = (-1.0 + std::sqrt((1.0 - A) * (1.0 - A) + 4.0 * A * u))/A;
    }

    // Finally, converting cos_xi_star into xi_star, which is the variable we will use:
    double xi_star = std::acos(cos_xi_star);

    return std::make_pair(xi_star, phi_star);
}

///////////////////////////////
/// Particle decay function ///
///////////////////////////////
// Takes variables in the lab frame, converts them into 
// 1) Samples a decay angle wrlt the polarization vector (the polarization vector wouldn't be known in a real physics case, but here we do know it from the start!)
// 2) Samples an azimutal angle around that polarization (anything from 0 to 2pi, uniformly distributed)
// 3) Performs the decay with that given angle, preserving total energy and momenta of the original Lambda (check day 440)
// 4) Returns the 4-momentum of the proton (from which we can calculate phi and y for histogramming and averaging).
    // Do note that we will keep the information on which Lambda was the mother for that, which would involve background in a realistic scenario!
// After those steps, all that is left is to calculate the original polarization vector for each cell in our (y, phi, pT) histogram (or maybe y, phi, px, py, pz).
// One other interesting thing to test is to see if a (y, phi, dTheta, pT) or (y, phi, dTheta, px, py, pz) or (px, py, pz, dTheta) histogram would be enough to average 
    // around many jets with discretized angular separation dTheta, and calculate the full polarization vector.
// Initially, I will only calculate polarization in the 3 directions of the lab frame, to compare with what Vitor's code had as an output (Maybe by comparing those 4-
// dimensional histograms' slices), but later on it could be interesting to calculate

std::pair<TLorentzVector, TLorentzVector> Lambda_decay(TLorentzVector Lambda_4_momentum_lab, TVector3 P_Lambda_star, double xi_star, double phi_star){
    
    // Rewrote the code to accept a P_Lambda vector already in the rest frame of the Lambda, so no need of the following block:
        // This makes it so that we only have to calculate the Lorentz boost of this quantity once for both Lambda_decay and P_angle_proton, instead of one for each.
        // Will use this section's logic for the actual Lambda polarization boost from the lab frame to the Lambda rest frame
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // ////////////////////
    // //// 1 - Boosting the polarization to the lambda rest frame:
    // ////////////////////
    //
    // TVector3 beta_inverse = -Lambda_4_momentum_lab.BoostVector();
    // TLorentzVector P_Lambda_lab_4vec(P_Lambda_lab.X(), P_Lambda_lab.Y(), P_Lambda_lab.Z(), 0); // TLorentzVectors have boost properties!
    // TLorentzVector P_Lambda_star_4vec = P_Lambda_lab_4vec;
    // P_Lambda_star_4vec.Boost(beta_inverse); // This boosts the P_Lambda_star_4vec into the rest frame.
    //
    // // Getting just the 3 needed components of this vector, in a shorter form for the function.
    //     // Also, only the spatial part is what matters in this case
    // TVector3 P_Lambda_star = P_Lambda_star_4vec.Vect();
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////
    //// 2 - Now decaying the particle given the xi_star and phi_star decay angles:
    ////////////////////
    // double daughter_momentum = std::sqrt((mass_Lambda*mass_Lambda - (mass_proton + mass_pi)*(mass_proton + mass_pi))*
    //                               (mass_Lambda*mass_Lambda - (mass_proton - mass_pi)*(mass_proton - mass_pi)))
    //                               /(2*mass_Lambda); // This momentum is the same for the proton and the pion, only with opposite directions
        // Rewriting to use pre-calculated mass^2 values (and correcting a few parenthesis mistakes):
    double daughter_momentum = std::sqrt((mL2 - mp_plus_mpi2) *
                                         (mL2 - mp_minus_mpi2))
                               / (2.0 * mass_Lambda);
    double E_p_star = std::sqrt(mass_proton * mass_proton + daughter_momentum * daughter_momentum); // The proton energy, using the shared daughter momentum

    ////////////////////
    //// 3 - Now working with the angles and creating the proton_4_momentum:
    ////////////////////
        // xi_star is the polar angle wrlt the polarization
        // phi_star is the azimuth angle for the proton. Just a uniform sample between 0 and 2pi
        // The goal is to use both of these to determine the 4-momentum of the proton daughter

        // Method 1 - Defining the vectors in a local frame, where the polarization points in the Z axis, and then using a rotation matrix:
    // // First, calculating the momenta in a frame where the polarization is in the z axis (easier to define)
    // double px_star = daughter_momentum * std::sin(xi_star) * std::cos(phi_star);
    // double py_star = daughter_momentum * std::sin(xi_star) * std::sin(phi_star);
    // double pz_star = daughter_momentum * std::cos(xi_star); // This means the z axis points in the polarization direction! Try xi_star = 0 and xi_star = pi to convince yourself.
    // // TLorentzVector proton_4_momentum(px_star, py_star, pz_star, E_p_star);

    // ////////////////////
    // //// 4 - Now rotating this frame into the actual xyz coordinates where the polarization points
    // ////////////////////
    //
    //
    // // The idea is that we built the proton's XYZ momenta in a coordinate system where the Z axis is the direction of the polarization.
    // // Now we will take that polarization vector and use it as reference to rotate our proton momentum into the usual XYZ definition, with Z in the beam direction.
    //     // Do notice that this here might be a problem if the polarization vector is null, but just assume it isn't!
    //     // Also, do notice that you couldn't do this with a momentum_Lambda_star vector and its angle with the P* vector, because the Lambda momentum in its own rest frame is zero!
    // TVector3 P_Lambda_star_unit = P_Lambda_star.Unit(); // The direction of the polarization, which we defined to be in the z axis for the (px,py,pz) vector above.

    //     // Do notice that this code only matches the local z axis (the direction of the polarization) with the laboratory's Z axis (the direction of the beam)!
    //     // It does not match the actual x and y axes of the laboratory! 
    //     // But that is not so much of a problem: the direction of the xy component is totally randomized (it comes from the sampling of a random phi_star
    //     // angle between 0 and 2pi) and does not carry any important physical meaning in this case! So there is no problem with how exactly you define this
    //     // proton's x and y axes components: they would be randomized and you just added a systematic shift to what is already a random uniform distribution.
    //     // This makes no difference and you can interpret it as a slightly more complex way of picking the random xy component from phi_star.
    //     // Just be careful not to say that the final proton's azimuthal angle is the same as this phi_star angle in a well defined rotation!
    // // Building a rotation matrix around that will take this local z axis definition into the direction of the polarization in the usual XYZ coordinates:
    // TVector3 z_axis(0, 0, 1);
    // TVector3 rotation_axis = z_axis.Cross(P_Lambda_star_unit);
    // // double rotation_angle = std::acos(z_axis.Dot(P_Lambda_star_unit));
    //     // I will actually use the cosine of this dot product above:
    // double cos_angle = z_axis.Dot(P_Lambda_star_unit);
    // double rotation_angle = std::acos(cos_angle);

    // TRotation rotation_matrix;

    // if (rotation_axis.Mag() > 1e-6){
    //     rotation_axis = rotation_axis.Unit();
    //     rotation_matrix.Rotate(rotation_angle, rotation_axis);
    // }
    // else{
    //     // A particular case, where the angle is way too small, i.e., the local z axis is already close to the laboratory Z axis.
    //     // Already aligned (or opposite) - either no rotation or do a 180 degrees rotation:
    //     if (cos_angle < 0){
    //         // 180 degree rotation about X (or any axis perpendicular to z)
    //         rotation_matrix.Rotate(TMath::Pi(), TVector3(1, 0, 0));
    //     } 
    //     else{
    //         // No rotation needed
    //         rotation_matrix = TRotation(); // identity
    //     }
    // }
    //
    // TVector3 p_star_vec(px_star, py_star, pz_star);
    // TVector3 p_star_rotated = rotation_matrix * p_star_vec; // Rotating the local-frame into the usual lab-frame with the Z axis pointing in the beam direction.
    // // The other way around would be a transpose of this rotation matrix.


        // Method 2 - Defining the vector straight in the laboratory's frame of reference:
    TVector3 P_Lambda_star_unit = P_Lambda_star.Unit(); // The reference for our sampling, i.e., the polarization vector in the rest frame, and in the laboratory's coordinate system.

    // Pick arbitrary vector not parallel to P_Lambda_star_unit:
    TVector3 T;
    if (std::abs(P_Lambda_star_unit.Z()) < 0.99){
        T = TVector3(0, 0, 1);
    }
    else{
        T = TVector3(0, 1, 0);
    }

    // Compute orthonormal basis
    TVector3 u_hat = (T.Cross(P_Lambda_star_unit)).Unit();
    TVector3 v_hat = (P_Lambda_star_unit.Cross(u_hat)).Unit();

    // Build p_star -- The proton's momentum vector in the usual XYZ axes of the laboratory:
        // See Le Bellac's Fig. 10.1 for the idea of this rotation: we align the local Z axis (defined as the direction where P_Lambda_star_unit points) and
        // align that with the usual Z direction of the laboratory (so we take fig. 10.1's z axis and shift it into \hat{n} = Z_lab)
        // Even though u_hat and v_hat are defined in an arbitrary way for each particle's decay (and in a biased way, not a truly random way), 
        // the phi_star sampling renders that bias harmless (it acts as just another layer of multiplicative constants on the random phi_star value!)
        // See the document "Coordinate System Transformation for Lambda Polarization Studies.pdf" for more information on this discussion.
    TVector3 p_star = daughter_momentum * (std::cos(xi_star) * P_Lambda_star_unit + std::sin(xi_star) * (std::cos(phi_star) * u_hat + std::sin(phi_star) * v_hat));

    ////////////////////
    //// 5 - Building the proton 4-momentum and boosting it into the lab frame:
    ////////////////////
        // For method 1:
    // TLorentzVector proton_4_momentum(p_star_rotated, E_p_star); // Short for (p_star_rotated.X(), p_star_rotated.Y(), p_star_rotated.Z(), E_p_star);
        // Modifying the vector in place with the Boost() method:
    // proton_4_momentum.Boost(-beta_inverse); // -beta_inverse is just beta, so this properly translates to the lab frame!

        // For method 2:
    TLorentzVector proton_4_momentum_star(p_star, E_p_star);
    TLorentzVector proton_4_momentum(p_star, E_p_star); // The same definition, but this one will be boosted.
    TVector3 beta = Lambda_4_momentum_lab.BoostVector(); // The boost that takes from the Lambda rest frame into the lab frame.
    proton_4_momentum.Boost(beta);
    
    return std::make_pair(proton_4_momentum, proton_4_momentum_star);
}

TVector3 boost_polarization_to_rest_frame(TLorentzVector Lambda_4_momentum_lab, TLorentzVector P_Lambda_lab_4vec){
    TVector3 beta_inverse = -Lambda_4_momentum_lab.BoostVector(); // TLorentzVectors have boost properties!
        // The definition below is wrong! I need the St components for the polarization vector too!
        // Reworked the code to include the St component now.
    // TLorentzVector P_Lambda_lab_4vec(P_Lambda_lab.X(), P_Lambda_lab.Y(), P_Lambda_lab.Z(), 0);
    TLorentzVector P_Lambda_star_4vec = P_Lambda_lab_4vec;
    P_Lambda_star_4vec.Boost(beta_inverse); // This boosts the P_Lambda_star_4vec into the rest frame.

    // Getting just the 3 needed components of this vector, in a shorter form for the function.
        // Also, only the spatial part is what matters in this case, because we are no longer doing any Lorentz transforms, so 'am converting it to a TVector3
    TVector3 P_Lambda_star = P_Lambda_star_4vec.Vect();

    ////////////////////////////////////////////////////
    // std::cout << P_Lambda_star_4vec.T() << std::endl; // This should be zero, as it is the polarization's time component in the rest frame!
    // Also zero within machine precision (within the float precision of the provided .dat data, actually)
    ////////////////////////////////////////////////////
    return P_Lambda_star;
}




/////////////////////////////////
/// Defining getter functions ///
/////////////////////////////////
// Receives 2D vectors (matrices) that have the event index in the first dimension and the particle index in the second.
// Retrieves quantities in the lab frame!
// If you'd like to read the data Vitor processed, just look at files like these:
// TString output_ebe = "/storage1/vribeiro/work_station/analysis/Polarization/data_ebe/60_70_bullet/event"+std::to_string(i)+".root";
void get_lambda(DoubleMatrix &y_matrix, DoubleMatrix &phi_matrix, DoubleMatrix &px_matrix, DoubleMatrix &py_matrix,
                DoubleMatrix &pz_matrix, DoubleMatrix &pT_matrix, DoubleMatrix &E_matrix, DoubleMatrix& mult_matrix,
                DoubleMatrix &St_matrix, DoubleMatrix &Sx_matrix, DoubleMatrix &Sy_matrix, DoubleMatrix &Sz_matrix,
                int N_events, bool with_bullet){
    std::cout << "\tStarting event loop for Lambda getter..." << std::endl;
    for(int i = 0; i < N_events; ++i){
        if (i % int (N_events * 0.1) == 0){ // Keeping track of the proccess for every 10% of the events
            std::cout << "\tNow on event " << i << " of " << N_events  << " (" << i * 1./N_events * 100 << "%)" << std::endl; // The 1. comes to represent the percentage as a double
        }

        std::vector<double> y, phi, pT, mult, St, Sx, Sy, Sz; // Declaring variables that will store the Lambda information for the current event
        double const mass = 1.115683; // The Lambda hyperon mass
		
		TString input_file_path;
        if (with_bullet){
            input_file_path = "/storage1/vribeiro/lumpy_events_40_50/bullet/";
        }
        else{
            input_file_path = "/storage1/vribeiro/lumpy_events_40_50/no_bullet/";
        }
        input_file_path += std::to_string(i)  + "/results/Smu_dpTdphidy_Thermal_rapidity_3122.dat";

		// cout << " " << endl;
		// cout << "Reading file  "+ std::to_string(i)  +  "..." << endl;

        // Load data in a TString -- Kept the TString buffering method (these files have about 5 MB!)
        std::ifstream t(input_file_path.Data());
        std::stringstream buffer;
        buffer << t.rdbuf();
        TString iSS_data = TString(buffer.str());

        //Break into several lines
        TObjArray* iSS_lines = iSS_data.Tokenize("\n");
        iSS_lines->SetOwner(kTRUE);
        for (TObject* line_obj : *iSS_lines){
            TString line = ((TObjString*) line_obj)->GetString();
            //If not header, break the lines in columns and stores in data members
            if ( !TString(line(0,1).Data()).Contains("#") ){
                TObjArray* columns = line.Tokenize("  ");
                y.push_back   (std::stod(((TObjString*) (*columns)[0])->GetString().Data()));
                pT.push_back  (std::stod(((TObjString*) (*columns)[1])->GetString().Data()));
                phi.push_back (std::stod(((TObjString*) (*columns)[2])->GetString().Data())); // phi in [0, 2*pi]
                mult.push_back(std::stod(((TObjString*) (*columns)[3])->GetString().Data())); // dN/(pTdpTdphidy)[GeV^-2] distribution, i.e., where this particle would be in a TH3D histogram.
                St.push_back  (std::stod(((TObjString*) (*columns)[4])->GetString().Data()));
                Sx.push_back  (std::stod(((TObjString*) (*columns)[5])->GetString().Data()));
                Sy.push_back  (std::stod(((TObjString*) (*columns)[6])->GetString().Data()));
                Sz.push_back  (std::stod(((TObjString*) (*columns)[7])->GetString().Data()));
                // delete columns;
            }
        }
        delete iSS_lines;

        // Processing the data into the 4-vector of the particle:
        std::vector<double> px, py, pz, E;
        int npoints = y.size();
		  
		for(int ip=0; ip<npoints;++ip){
			px.push_back(pT[ip]*std::cos(phi[ip]));
            py.push_back(pT[ip]*std::sin(phi[ip]));
            double mT = std::sqrt(mass*mass + pT[ip]*pT[ip]);
            pz.push_back(mT*std::sinh(y[ip])); 
            E.push_back(std::sqrt(mass*mass + px[ip]*px[ip] + py[ip]*py[ip] + pz[ip]*pz[ip]));   
        }

        // Finally, passing all these values into a larger array that contains all events:
        y_matrix.push_back(y);
        phi_matrix.push_back(phi);
        mult_matrix.push_back(mult);
        
        px_matrix.push_back(px);
        py_matrix.push_back(py);
        pz_matrix.push_back(pz);
        pT_matrix.push_back(pT);
        E_matrix.push_back(E);

        St_matrix.push_back(St);
        Sx_matrix.push_back(Sx);
        Sy_matrix.push_back(Sy);
        Sz_matrix.push_back(Sz);
    }
    std::cout << "\tDone getting Lambda information!" << std::endl;
}

// This function retrieves the same set of data that Vitor used in Pol_Analysis_Random_hist_ebe.C
// It modifies the provided vectors and stores the jet information.
void get_jet(std::vector<double> &n_event, std::vector<double> &n_random, std::vector<double> &phi_random,
             std::vector<double> &momentum_x, std::vector<double> &momentum_y){
    std::cout << "\n\tStarting event loop for Jet getter (even if the w/ jets flag is off, so that v_2 is averaged out!)..." << std::endl;
    
    // Getting a file that has jets in hydro set in a random direction for the jet momentum deposition:
    const char* file_random_jet_path = "/storage1/vribeiro/lumpy_events_40_50/random_parameters.dat";

        // Skipping the whole buffering process Vitor did -- The files in only about ~10 kB, so it can be read at once:
    std::ifstream file_random_jet(file_random_jet_path);
    std::string content((std::istreambuf_iterator<char>(file_random_jet)), std::istreambuf_iterator<char>());
    TString tcontent(content.c_str()); // Converting to a TString in order to use the Tokenize() method

    //Break into several lines
    TObjArray* random_lines = tcontent.Tokenize("\n"); // tcontent is the same as the "random_data" variable from Vitor's code.
    random_lines->SetOwner(kTRUE); // If you don't do this, then ROOT will not delete the contents of TObjArray! Could introduce leaks! Btw, kTRUE is defined in TString.h
    for (TObject* line_obj : *random_lines){ // Iterates over random_lines objects with line_obj being the variable that will assume the random_lines values. This is a range-based loop.
        TString line = ((TObjString*) line_obj)->GetString(); // Getting the string in the TObject
        // If not a header or comment, break the lines in columns and stores in data members
        if ( !TString(line(0,1).Data()).Contains("#") ){
            TObjArray* columns = line.Tokenize("  ");

            n_event.push_back   (std::stod(((TObjString*) (*columns)[0])->GetString().Data()));
            n_random.push_back  (std::stod(((TObjString*) (*columns)[1])->GetString().Data()));
            phi_random.push_back (std::stod(((TObjString*) (*columns)[2])->GetString().Data())); //phi in [0, 2*pi]
            momentum_x.push_back(std::stod(((TObjString*) (*columns)[3])->GetString().Data()));
            momentum_y.push_back  (std::stod(((TObjString*) (*columns)[4])->GetString().Data()));
            // delete columns;
        }
    }
    delete random_lines;

    std::cout << "\tDone getting jet information!" << std::endl;
}

// A generalization of previous code for passing multiple particles to be sampled at once.
    // This one is not useless, just a more general version that is kinda overkill for the way I designed the later code.
//     // Will pass P_Lambda_star_array, an array with many different Lambda polarization modules, already in the rest frame, and retrieve a samples vector from it.
//     // Usage example:
//     // std::vector<double> P_Lambda_star_array = {0.2, 0.5, -0.3, 0.0, 1.0};  // One P_Lambda per particle
//     // auto [xi_star_vec, phi_vec] = sample_P_angle_proton_for_each_PLambda(P_Lambda_star_array);
// std::pair<std::vector<double>, std::vector<double>> sample_P_angle_proton_array_modification(const std::vector<double> &P_Lambda_star_array){
//     std::vector<double> samples;
//     std::vector<double> samples_azimuth;
//     samples.reserve(P_Lambda_star_array.size());  // Reserve space for performance
//     samples_azimuth.reserve(P_Lambda_star_array.size()); // The azimutal angle from 0 to 2 pi
// 
//     std::mt19937 rng {std::random_device{}()};  // Creates a random_device object (temporary) with the uniform initialization on newer C++ (the "{}"),
//                                                 // then calls it with () getting a random seed, and then that seed is passed to the random engine.
// 
//     // Uniform distribution for angle xi_star between 0 and pi
//     std::uniform_real_distribution<double> dist_x(0.0, PI);
//     std::uniform_real_distribution<double> dist_azimuth(0.0, 2*PI);
//
//     // Loop over all provided P_Lambda values
//     for (double P_Lambda : P_Lambda_star_array){
//         // Compute P_max (the maximum possible value of P(x))
//         // This occurs at x = 0 or x = pi depending on sign of cos(x) term
//         double P_max = (1.0 / PI) * (1.0 + std::abs(alpha_H * P_Lambda));
// 
//         // Uniform distribution for y between 0 and P_max
//         // This is used for the rejection sampling step
//         std::uniform_real_distribution<double> dist_y(0.0, P_max);
// 
//         // Rejection sampling loop - generate one sample for this P_Lambda
//         while (true){
//             // Propose an x uniformly in [0, pi]
//             double x = dist_x(rng);
// 
//             // Propose a y uniformly in [0, P_max]
//             double y = dist_y(rng);  // Draws from the maximum value distribution independently.
//                                      // This makes it so that we always follow the appropriate limit and have correctly weighted samples of x(!)
// 
//             // Compute the actual probability density at x
//             double px = P_angle_proton(x, P_Lambda);
// 
//             // Accept x if y < P(x). This ensures samples are distributed according to P(x): larger P(x) values will be accepted more easily!
//             if (y < px){
//                 samples.push_back(x);  // Store accepted x
//                 break;                // Done with this P_Lambda, move to next
//             }
//             // Otherwise, reject and try again
//         }
// 
//         // Now sampling the azimuthal angle:
//         double phi = dist_azimuth(rng);
//         samples_azimuth.push_back(phi);
//     }
//     return std::make_pair(samples, samples_azimuth);
// }



// Useless, older code that actually made no sense to use:
// Another version, that could sample the whole distribution from the same P_Lambda value (not quite what I want: one sample per each P_Lambda value)
// std::vector<double> sample_P_angle_proton_multiple_times(size_t N_part, double P_Lambda){
//     std::vector<double> samples;
//     samples.reserve(N_part); // Reserve space for performance

//     std::mt19937 rng {std::random_device{}()}; // Creates a random_device object (temporary) with the uniform initialization on newer C++ (the "{}"),
//                                                // then calls it with () getting a random seed, and then that seed is passed to the random engine.

//     // Uniform distribution for angle xi_star between 0 and pi
//     std::uniform_real_distribution<double> dist_x(0.0, PI);

//     // We compute P_max (the maximum possible value of P(x))
//     // This occurs at x = 0 or x = pi depending on sign of cos(x) term
//     double P_max = (1.0 / PI) * (1.0 + std::abs(alpha_H * P_Lambda));

//     // Uniform distribution for y between 0 and P_max
//     // This is used for the rejection sampling step
//     std::uniform_real_distribution<double> dist_y(0.0, P_max);

//     size_t N_part = samples.size();  // Number of particles (or samples) we want to generate
//     size_t i = 0;                    // Counter for accepted samples

//     // Rejection sampling loop
//     while (i < N_part){
//         // Propose an x uniformly in [0, pi]
//         double x = dist_x(rng);

//         // Propose a y uniformly in [0, P_max]
//         double y = dist_y(rng); // Draws from the maximum value distribution independently.
//                                 // This makes it so that we always follow the appropriate limit and have correctly weighted samples of x(!)

//         // Compute the actual probability density at x
//         double px = P_angle_proton(x, alpha_H, P_Lambda);

//         // Accept x if y < P(x). This ensures samples are distributed according to P(x): larger P(x) values will be accepted more easily!
//         if (y < px){
//             // samples[i++] = x;  // Store accepted x and increment sample counter
//             samples.push_back(x); // Store accepted x in the new non-predefined vector
//             i++;
//         }
//         // Otherwise, reject and try again
//     }
// }