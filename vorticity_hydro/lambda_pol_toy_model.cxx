#include "TROOT.h"
#include "TFile.h"
// #include "TTree.h"

#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRotation.h"

#include <vector>
#include <array>
#include <cmath>
#include <random>
// #include "math.h" // Same as cmath!
#include "TString.h"
#include "TObjArray.h"
#include "TMath.h"

#include <dirent.h>
#include <iostream>
#include <string>
#include <fstream>

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
                DoubleMatrix &pz_matrix, DoubleMatrix &E_matrix, DoubleMatrix& mult_matrix, DoubleMatrix &St_matrix,
                DoubleMatrix &Sx_matrix, DoubleMatrix &Sy_matrix, DoubleMatrix &Sz_matrix, int N_events);
void get_jet(std::vector<double> &n_event, std::vector<double> &n_random, std::vector<double> &phi_random,
             std::vector<double> &momentum_x, std::vector<double> &momentum_y);

    // Physics part functions
std::pair<double, double> sample_P_angle_proton(double P_Lambda_star_mag, std::mt19937 rng, std::uniform_real_distribution<double> dist_x,
                                                std::uniform_real_distribution<double> dist_y, std::uniform_real_distribution<double> dist_azimuth);
std::pair<TLorentzVector, TLorentzVector> Lambda_decay(TLorentzVector Lambda_4_momentum_lab, TVector3 P_Lambda_star, double xi_star, double phi_star);
TVector3 boost_polarization_to_rest_frame(TLorentzVector Lambda_4_momentum_lab, TLorentzVector P_Lambda_lab_4vec);

// This code assumes we are in Jarvis4! Also, the getter functions are based on Pol_Analysis_Random_hist_ebe.C from the HadrEx_Ph repository
int lambda_pol_toy_model(){
    ////////////////////
    //// 1 - Initializing variables and histograms:
    ////////////////////
    std::cout << "Initializing variables" << std::endl;
    int N_events = 250; // There are only 250 events in the desired folder. Possibly oversampled to get statistics.

    TH1D *hLambdaCounter = new TH1D ("hLambdaCounter", "", 1, -1, 1);
    
    // Histograms to compare the true magnitude of the Lambda polarizations and the reconstructed magnitude:
    int N_bins_pol = 200;
    TH1D *hLambdaPolMag = new TH1D("hLambdaPolMag", "hLambdaPolMag", N_bins_pol, 0, 1); // Should only go from -1 to 1, but I want to check it!
    TH1D *hLambdaPolStarMagReco = new TH1D("hLambdaPolStarMagReco", "hLambdaPolStarMagReco", N_bins_pol, 0, 1);

        // Checking component by component to see if we can properly reconstruct the average of polarization in each direction
        // (in this first reconstruction method, I am not aiming to reconstruct the polarization around a jet, but to reconstruct
        // the "global" polarization components in the XYZ axes. Thus, the Reco value will be a single entry, while the true values
        // will actually give me distributions of polarizations that I ended up averaging on)
    TH1D *hLambdaPolX = new TH1D("hLambdaPolX", "hLambdaPolX", N_bins_pol, -0.5, 0.5); // Polarization is normalized to 1, and usually does not go over 0.2 !
    TH1D *hLambdaPolStarXReco = new TH1D("hLambdaPolStarXReco", "hLambdaPolStarXReco", N_bins_pol, -0.5, 0.5);

    TH1D *hLambdaPolY = new TH1D("hLambdaPolY", "hLambdaPolY", N_bins_pol, -0.5, 0.5);
    TH1D *hLambdaPolStarYReco = new TH1D("hLambdaPolStarYReco", "hLambdaPolStarYReco", N_bins_pol, -0.5, 0.5);

    TH1D *hLambdaPolZ = new TH1D("hLambdaPolZ", "hLambdaPolZ", N_bins_pol, -0.5, 0.5);
    TH1D *hLambdaPolStarZReco = new TH1D("hLambdaPolStarZReco", "hLambdaPolStarZReco", N_bins_pol, -0.5, 0.5);

        // Now declaring 2D histograms that will have the average polarization in the Z axis, and (pT, y) in the plane:
    int N_bins_pT = 50;
    int N_bins_rap = 50;
    int N_bins_phi = 50;

        // These TH2Ds here will be defined as clones later, so it is clearer to the reader to show that they will not need full definitions:
        // (see below)
    // TH2D *hLambdaPolX_pT_yReco = new TH2D("hLambdaPolX_pT_yReco", "hLambdaPolX_pT_yReco", N_bins_pT, 0, 5, N_bins_rap, -3, 3);
    // TH2D *hLambdaPolY_pT_yReco = new TH2D("hLambdaPolY_pT_yReco", "hLambdaPolX_pT_yReco", N_bins_pT, 0, 5, N_bins_rap, -3, 3);
    // TH2D *hLambdaPolZ_pT_yReco = new TH2D("hLambdaPolZ_pT_yReco", "hLambdaPolX_pT_yReco", N_bins_pT, 0, 5, N_bins_rap, -3, 3);

        // This simulation only has |y|<0.5 and pT between 0 and 3 GeV/c.
    TH2D *hLambdaPolX_pT_y = new TH2D("hLambdaPolX_pT_y", "hLambdaPolX_pT_y", N_bins_pT, 0, 3, N_bins_rap, -0.5, 0.5);
    TH2D *hLambdaPolStarX_pT_yReco;

    TH2D *hLambdaPolY_pT_y = new TH2D("hLambdaPolY_pT_y", "hLambdaPolX_pT_y", N_bins_pT, 0, 3, N_bins_rap, -0.5, 0.5);
    TH2D *hLambdaPolStarY_pT_yReco;

    TH2D *hLambdaPolZ_pT_y = new TH2D("hLambdaPolZ_pT_y", "hLambdaPolX_pT_y", N_bins_pT, 0, 3, N_bins_rap, -0.5, 0.5);
    TH2D *hLambdaPolStarZ_pT_yReco;

            // Declaring a counter histogram and some summation-storing histograms to define polarization in each (pT, y) bin:
    TH2D *hLambdaCounter_pT_y = new TH2D("hLambdaCounter_pT_y", "hLambdaCounter_pT_y", N_bins_pT, 0, 3, N_bins_rap, -0.5, 0.5); // Counts the amount of Lambdas in each (pT, y) bin

    TH2D *hLambdaAvgDotX_pT_y = new TH2D("hLambdaAvgDotX_pT_y", "hLambdaAvgDotX_pT_y", N_bins_pT, 0, 3, N_bins_rap, -0.5, 0.5); // Receives the same sum as average_dotX, but for each (pT, y) bin
    TH2D *hLambdaAvgDotY_pT_y = new TH2D("hLambdaAvgDotY_pT_y", "hLambdaAvgDotY_pT_y", N_bins_pT, 0, 3, N_bins_rap, -0.5, 0.5);
    TH2D *hLambdaAvgDotZ_pT_y = new TH2D("hLambdaAvgDotZ_pT_y", "hLambdaAvgDotZ_pT_y", N_bins_pT, 0, 3, N_bins_rap, -0.5, 0.5);

        // Declaring some histograms to calculate the polarization in the lab frame -- You will need the phi information to build the full 4-vector of the "mean Lambda" in that bin:
        // In a later version, I can just eliminate the 2D histogram and work with this one instead, even for the P* calculations only. Just sum on all phi values.
    TH3D *hLambdaCounter_pT_y_phi = new TH3D("hLambdaCounter_pT_y_phi", "hLambdaCounter_pT_y_phi", N_bins_pT, 0, 3, N_bins_rap, -0.5, 0.5, N_bins_phi, 0, 2*PI); // Counts the amount of Lambdas in each (pT, y, phi) bin
    TH3D *hLambdaAvgDotX_pT_y_phi = new TH3D("hLambdaAvgDotX_pT_y_phi", "hLambdaAvgDotX_pT_y_phi", N_bins_pT, 0, 3, N_bins_rap, -0.5, 0.5, N_bins_phi, 0, 2*PI);
    TH3D *hLambdaAvgDotY_pT_y_phi = new TH3D("hLambdaAvgDotY_pT_y_phi", "hLambdaAvgDotY_pT_y_phi", N_bins_pT, 0, 3, N_bins_rap, -0.5, 0.5, N_bins_phi, 0, 2*PI);
    TH3D *hLambdaAvgDotZ_pT_y_phi = new TH3D("hLambdaAvgDotZ_pT_y_phi", "hLambdaAvgDotZ_pT_y_phi", N_bins_pT, 0, 3, N_bins_rap, -0.5, 0.5, N_bins_phi, 0, 2*PI);

    TH3D *hLambdaPolStarX_pT_y_phiReco;
    TH3D *hLambdaPolStarY_pT_y_phiReco;
    TH3D *hLambdaPolStarZ_pT_y_phiReco;

        // The lab polarization histograms:
    TH3D *hLambdaPolX_pT_y_phiReco = new TH3D("hLambdaPolX_pT_y_phiReco", "hLambdaPolX_pT_y_phiReco", N_bins_pT, 0, 3, N_bins_rap, -0.5, 0.5, N_bins_phi, 0, 2*PI);
    TH3D *hLambdaPolY_pT_y_phiReco = new TH3D("hLambdaPolY_pT_y_phiReco", "hLambdaPolY_pT_y_phiReco", N_bins_pT, 0, 3, N_bins_rap, -0.5, 0.5, N_bins_phi, 0, 2*PI);
    TH3D *hLambdaPolZ_pT_y_phiReco = new TH3D("hLambdaPolZ_pT_y_phiReco", "hLambdaPolZ_pT_y_phiReco", N_bins_pT, 0, 3, N_bins_rap, -0.5, 0.5, N_bins_phi, 0, 2*PI);

        // Some 2D lab polarization histograms to compare easily with hLambdaPolX_pT_y:
    TH2D *hLambdaPolX_pT_yReco = new TH2D("hLambdaPolX_pT_yReco", "hLambdaPolX_pT_yReco", N_bins_pT, 0, 3, N_bins_rap, -0.5, 0.5);
    TH2D *hLambdaPolY_pT_yReco = new TH2D("hLambdaPolY_pT_yReco", "hLambdaPolY_pT_yReco", N_bins_pT, 0, 3, N_bins_rap, -0.5, 0.5);
    TH2D *hLambdaPolZ_pT_yReco = new TH2D("hLambdaPolZ_pT_yReco", "hLambdaPolZ_pT_yReco", N_bins_pT, 0, 3, N_bins_rap, -0.5, 0.5);

            // Other 2D histograms to compare with Vitor's paper:
    TH2D *hLambdaPolX_phi_pTReco = new TH2D("hLambdaPolX_phi_pTReco", "hLambdaPolX_phi_pTReco", N_bins_phi, 0, 2*PI, N_bins_pT, 0, 3);
    TH2D *hLambdaPolY_phi_pTReco = new TH2D("hLambdaPolY_phi_pTReco", "hLambdaPolY_phi_pTReco", N_bins_phi, 0, 2*PI, N_bins_pT, 0, 3);
    TH2D *hLambdaPolZ_phi_pTReco = new TH2D("hLambdaPolZ_phi_pTReco", "hLambdaPolZ_phi_pTReco", N_bins_phi, 0, 2*PI, N_bins_pT, 0, 3);

        // Finally, 1D lab polarization histograms to compare with the plots in the paper:
    TH1D *hLambdaPolX_phiReco = new TH1D("hLambdaPolX_phiReco", "hLambdaPolX_phiReco", N_bins_phi, 0, 2*PI);
    TH1D *hLambdaPolY_phiReco = new TH1D("hLambdaPolY_phiReco", "hLambdaPolY_phiReco", N_bins_phi, 0, 2*PI);
    TH1D *hLambdaPolZ_phiReco = new TH1D("hLambdaPolZ_phiReco", "hLambdaPolZ_phiReco", N_bins_phi, 0, 2*PI);

    ////////////////////
    //// 1.1 - Defining randomized samplers for the loop -- It would make no sense to define them inside the loop, as they will not change:
    ////////////////////
    std::mt19937 rng {std::random_device{}()};  // Creates a random_device object (temporary) with the uniform initialization on newer C++ (the "{}"),
                                                // then calls it with () getting a random seed, and then that seed is passed to the random engine.
        // Uniform distribution for angle xi_star between 0 and pi:
    std::uniform_real_distribution<double> dist_x(0.0, PI);
    std::uniform_real_distribution<double> dist_azimuth(0.0, 2*PI);

    ////////////////////
    //// 2 - Getter function calls:
    ////////////////////
    std::cout << "\nFetching data..." << std::endl;
    DoubleMatrix y_matrix, phi_matrix, px_matrix, py_matrix, pz_matrix, E_matrix, mult_matrix, St_matrix, Sx_matrix, Sy_matrix, Sz_matrix;
    get_lambda(y_matrix, phi_matrix, px_matrix, py_matrix, pz_matrix, E_matrix, mult_matrix, St_matrix, Sx_matrix, Sy_matrix, Sz_matrix, N_events);
    std::cout << "\tDone!" << std::endl;

    ////////////////////
    //// 3 - Decaying Lambdas and filling histograms:
    ////////////////////
    std::cout << "\nDecaying Lambdas and filling histograms" << std::endl;
        // Looping on all events and all particles of each event:
    TVector3 x_hat(1, 0, 0);
    TVector3 y_hat(0, 1, 0);
    TVector3 z_hat(0, 0, 1);

    double average_dotX = 0;
    double average_dotY = 0;
    double average_dotZ = 0;
    int N_resamples = 1; // Goes through each particle N_resamples # of times. This is an attempt to see if the polarization estimate values become more stable!
    for (int resample_idx = 0; resample_idx < N_resamples; resample_idx++){
        std::cout << "Now on resample " << std::to_string(resample_idx + 1) << " of " + std::to_string(N_resamples) << std::endl;
        for (int ev_idx = 0; ev_idx < N_events; ev_idx++){
            for (int particle_idx = 0; particle_idx < y_matrix[ev_idx].size(); particle_idx++){
                hLambdaCounter->Fill(0);

                // 1 - Fetching particle information:
                TLorentzVector Lambda_4_momentum_lab(px_matrix[ev_idx][particle_idx], py_matrix[ev_idx][particle_idx], pz_matrix[ev_idx][particle_idx], E_matrix[ev_idx][particle_idx]);
                // TVector3 P_Lambda_lab(Sx_matrix[ev_idx][particle_idx], Sy_matrix[ev_idx][particle_idx], Sz_matrix[ev_idx][particle_idx]);
                TLorentzVector P_Lambda_lab_4vec(Sx_matrix[ev_idx][particle_idx], Sy_matrix[ev_idx][particle_idx], Sz_matrix[ev_idx][particle_idx], St_matrix[ev_idx][particle_idx]); // You need the temporal component too!
                hLambdaPolMag->Fill((P_Lambda_lab_4vec.Vect()).Mag()); // I just want the magnitude of the spatial part for this plot, so just grab the Vect() magnitude

                    // Filling the known polarizations:
                hLambdaPolX->Fill(Sx_matrix[ev_idx][particle_idx]);
                hLambdaPolY->Fill(Sy_matrix[ev_idx][particle_idx]);
                hLambdaPolZ->Fill(Sz_matrix[ev_idx][particle_idx]);

                // 2 - Calculating the polarization in the Lambda rest frame:
                TVector3 P_Lambda_star = boost_polarization_to_rest_frame(Lambda_4_momentum_lab, P_Lambda_lab_4vec);
                double P_Lambda_star_mag = P_Lambda_star.Mag();

                // 3 - Sampling decay angles:
                    // Calculating the P_max value for the rejection sampling:
                    // Compute P_max (the maximum possible value of P(x))
                    // This occurs at x = 0 or x = pi depending on sign of cos(x) term
                double P_max = (1.0 / PI) * (1.0 + std::abs(alpha_H * P_Lambda_star_mag));

                    // Defining the uniform distribution dist_y for the sampling:
                    // Uniform distribution for y between 0 and P_max
                std::uniform_real_distribution<double> dist_y(0.0, P_max);

                auto [xi_star, phi_star] = sample_P_angle_proton(P_Lambda_star_mag, rng, dist_x, dist_y, dist_azimuth); // Uses the new unpacking of C++17
                
                // 4 - Generating the proton 4-momentum from the decay:
                    // Actually, I just need the angles at which the proton would decay, not the whole 4-momentum of the decay, but whatever, let's keep it!
                    // In other words, I could've stopped at the sample_P_angle_proton function, and then just rotate those angles to the XYZ axes of the lab frame.
                    // This could turn out to be useful later on, if I intend on doing some background checks for the Lambda reconstructions or something like that.
                    // You won't even need the proton's momentum for the ring observable's reconstruction! But whatever...
                auto [proton_4_momentum, proton_4_momentum_star] = Lambda_decay(Lambda_4_momentum_lab, P_Lambda_star, xi_star, phi_star);

                // 5 - Extracting useful variables from the 4-momentum of the Lambda:
                    // Earlier getters from the proton 4-momenta (not quite what I need right now)
                // double proton_y = proton_4_momentum.Rapidity();
                // double proton_phi = proton_4_momentum.Phi();
                
                    // Don't need to calculate these numbers! They are already provided in the code! (But I checked and the results match)
                // double lambda_y = Lambda_4_momentum_lab.Rapidity();
                // double lambda_pT = Lambda_4_momentum_lab.Pt();
                double lambda_y = y_matrix[ev_idx][particle_idx];
                double lambda_pT = std::sqrt(px_matrix[ev_idx][particle_idx]*px_matrix[ev_idx][particle_idx] + py_matrix[ev_idx][particle_idx]*py_matrix[ev_idx][particle_idx]);
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
                double X_dot = proton_star_unit_vector.Dot(x_hat);
                double Y_dot = proton_star_unit_vector.Dot(y_hat);
                double Z_dot = proton_star_unit_vector.Dot(z_hat);

                average_dotX += X_dot;
                average_dotY += Y_dot;
                average_dotZ += Z_dot;

                // 7 - Summing to calculate polarization on subsets of the available Lambdas -- Bins of (pT, y):
                hLambdaCounter_pT_y->Fill(lambda_pT, lambda_y);

                hLambdaAvgDotX_pT_y->Fill(lambda_pT, lambda_y, X_dot);
                hLambdaAvgDotY_pT_y->Fill(lambda_pT, lambda_y, Y_dot);
                hLambdaAvgDotZ_pT_y->Fill(lambda_pT, lambda_y, Z_dot);

                hLambdaCounter_pT_y_phi->Fill(lambda_pT, lambda_y, lambda_phi);
                hLambdaAvgDotX_pT_y_phi->Fill(lambda_pT, lambda_y, lambda_phi, X_dot);
                hLambdaAvgDotY_pT_y_phi->Fill(lambda_pT, lambda_y, lambda_phi, Y_dot);
                hLambdaAvgDotZ_pT_y_phi->Fill(lambda_pT, lambda_y, lambda_phi, Z_dot);

                // 8 - Filling the true values of the TH2D polarization histogram:
                hLambdaPolX_pT_y->Fill(lambda_pT, lambda_y, Sx_matrix[ev_idx][particle_idx]);
                hLambdaPolY_pT_y->Fill(lambda_pT, lambda_y, Sy_matrix[ev_idx][particle_idx]);
                hLambdaPolZ_pT_y->Fill(lambda_pT, lambda_y, Sz_matrix[ev_idx][particle_idx]);
            }
        }
    }
    std::cout << "\tDone!" << std::endl;

    // Finishing the average for the 3 global dot products:
    average_dotX /= hLambdaCounter->GetBinContent(1);
    average_dotY /= hLambdaCounter->GetBinContent(1);
    average_dotZ /= hLambdaCounter->GetBinContent(1);

    // Calculating the mean polarization in each of those three global axes -- IN THE LAMBDA REST FRAME!
    double PolX_star = 3./alpha_H * average_dotX;
    double PolY_star = 3./alpha_H * average_dotY;
    double PolZ_star = 3./alpha_H * average_dotZ;

    hLambdaPolStarXReco->Fill(PolX_star);
    hLambdaPolStarYReco->Fill(PolY_star);
    hLambdaPolStarZReco->Fill(PolZ_star);

    hLambdaPolStarMagReco->Fill(std::sqrt(PolX_star*PolX_star + PolY_star*PolY_star + PolZ_star*PolZ_star));

    // Averaging the TH2D histograms that collect the mean dot product, and then converting them into polarizations:
    hLambdaAvgDotX_pT_y->Divide(hLambdaCounter_pT_y);
    hLambdaAvgDotY_pT_y->Divide(hLambdaCounter_pT_y);
    hLambdaAvgDotZ_pT_y->Divide(hLambdaCounter_pT_y);

        // Cloning the histograms, to keep an average dot product histogram and a polarization histogram as separate entities:
    hLambdaPolStarX_pT_yReco = (TH2D*) hLambdaAvgDotX_pT_y->Clone("hLambdaPolStarX_pT_yReco");
    hLambdaPolStarY_pT_yReco = (TH2D*) hLambdaAvgDotY_pT_y->Clone("hLambdaPolStarY_pT_yReco");
    hLambdaPolStarZ_pT_yReco = (TH2D*) hLambdaAvgDotZ_pT_y->Clone("hLambdaPolStarZ_pT_yReco");

    hLambdaPolStarX_pT_yReco->SetTitle("hLambdaPolStarX_pT_yReco");
    hLambdaPolStarY_pT_yReco->SetTitle("hLambdaPolStarY_pT_yReco");
    hLambdaPolStarZ_pT_yReco->SetTitle("hLambdaPolStarZ_pT_yReco");
    
        // Applying the transform that turns them into polarizations -- In the lambda rest frame!
    hLambdaPolStarX_pT_yReco->Scale(3.0 / alpha_H);
    hLambdaPolStarY_pT_yReco->Scale(3.0 / alpha_H);
    hLambdaPolStarZ_pT_yReco->Scale(3.0 / alpha_H);

            // The same for the TH3D histograms:
    ///////////////////////////////////////////////////////////////
    hLambdaAvgDotX_pT_y_phi->Divide(hLambdaCounter_pT_y_phi);
    hLambdaAvgDotY_pT_y_phi->Divide(hLambdaCounter_pT_y_phi);
    hLambdaAvgDotZ_pT_y_phi->Divide(hLambdaCounter_pT_y_phi);

    hLambdaPolStarX_pT_y_phiReco = (TH3D*) hLambdaAvgDotX_pT_y_phi->Clone("hLambdaPolStarX_pT_y_phiReco");
    hLambdaPolStarY_pT_y_phiReco = (TH3D*) hLambdaAvgDotY_pT_y_phi->Clone("hLambdaPolStarY_pT_y_phiReco");
    hLambdaPolStarZ_pT_y_phiReco = (TH3D*) hLambdaAvgDotZ_pT_y_phi->Clone("hLambdaPolStarZ_pT_y_phiReco");

    hLambdaPolStarX_pT_y_phiReco->SetTitle("hLambdaPolStarX_pT_y_phiReco");
    hLambdaPolStarY_pT_y_phiReco->SetTitle("hLambdaPolStarY_pT_y_phiReco");
    hLambdaPolStarZ_pT_y_phiReco->SetTitle("hLambdaPolStarZ_pT_y_phiReco");
    
        // Applying the transform that turns them into polarizations -- In the lambda rest frame!
    hLambdaPolStarX_pT_y_phiReco->Scale(3.0 / alpha_H);
    hLambdaPolStarY_pT_y_phiReco->Scale(3.0 / alpha_H);
    hLambdaPolStarZ_pT_y_phiReco->Scale(3.0 / alpha_H);
    ///////////////////////////////////////////////////////////////

        // Now boosting each of these into the lab rest frame by calculating a mean 4-vector of all Lambdas in each (pT, y) bin:
    for (int pT_idx = 1; pT_idx <= N_bins_pT; pT_idx++){ // bin index starts at 1 (not 0)
        for (int y_idx = 1; y_idx <= N_bins_rap; y_idx++){
            for (int phi_idx = 1; phi_idx <= N_bins_phi; phi_idx++){
                // Getting all the needed values for all three components at once:
                double PolX_star_pTy_phi = hLambdaPolStarX_pT_y_phiReco->GetBinContent(pT_idx, y_idx, phi_idx);
                double pT_bin_center = hLambdaPolStarX_pT_y_phiReco->GetXaxis()->GetBinCenter(pT_idx);
                double y_bin_center = hLambdaPolStarX_pT_y_phiReco->GetYaxis()->GetBinCenter(y_idx);
                double phi_bin_center = hLambdaPolStarX_pT_y_phiReco->GetZaxis()->GetBinCenter(phi_idx);

                    // Getting the two other components, as the (pT, y) bin centers will be the same in all three cases:
                double PolY_star_pTy_phi = hLambdaPolStarY_pT_y_phiReco->GetBinContent(pT_idx, y_idx, phi_idx);
                double PolZ_star_pTy_phi = hLambdaPolStarZ_pT_y_phiReco->GetBinContent(pT_idx, y_idx, phi_idx);

                // Now calculating a mean 4-momentum of the Lambda in each bin:
                double mean_px = pT_bin_center * std::cos(phi_bin_center);
                double mean_py = pT_bin_center * std::sin(phi_bin_center);
                double mean_mT = std::sqrt(pT_bin_center * pT_bin_center + mL2);
                double mean_pz = mean_mT * std::sinh(y_bin_center);
                double mean_E  = mean_mT * std::cosh(y_bin_center);

                TLorentzVector mean_Lambda_4vec_lab(mean_px, mean_py, mean_pz, mean_E);

                // Building the mean polarization 4-momentum, assuming that S^0 = 0 from orthogonality:
                TLorentzVector P_Lambda_star_mean_4vec(PolX_star_pTy_phi, PolY_star_pTy_phi, PolZ_star_pTy_phi, 0);

                // Boosting the polarization in this bin:
                TVector3 beta = mean_Lambda_4vec_lab.BoostVector(); // The boost that takes from the Lambda rest frame into the lab frame.
                TLorentzVector P_Lambda_lab_mean_4vec = P_Lambda_star_mean_4vec;
                P_Lambda_lab_mean_4vec.Boost(beta);
                TVector3 P_Lambda_lab_mean = P_Lambda_lab_mean_4vec.Vect();

                hLambdaPolX_pT_y_phiReco->Fill(pT_bin_center, y_bin_center, phi_bin_center, P_Lambda_lab_mean.X());
                hLambdaPolY_pT_y_phiReco->Fill(pT_bin_center, y_bin_center, phi_bin_center, P_Lambda_lab_mean.Y());
                hLambdaPolZ_pT_y_phiReco->Fill(pT_bin_center, y_bin_center, phi_bin_center, P_Lambda_lab_mean.Z());

                // Also doing a projection in 2D, averaging over phi bins:
                    // You also need to average over the number of bins that you are summing! If polarization = 1 in one, and 1 in another bin, then
                    // the sum would be 2, and polarization can't ever go up to that value! You need to properly average it by the number of N_phi bins.
                hLambdaPolX_pT_yReco->Fill(pT_bin_center, y_bin_center, P_Lambda_lab_mean.X()/N_bins_phi);
                hLambdaPolY_pT_yReco->Fill(pT_bin_center, y_bin_center, P_Lambda_lab_mean.Y()/N_bins_phi);
                hLambdaPolZ_pT_yReco->Fill(pT_bin_center, y_bin_center, P_Lambda_lab_mean.Z()/N_bins_phi);

                    // A second 2D projection, in the same way it was plotted on Vitor's paper:
                hLambdaPolX_phi_pTReco->Fill(phi_bin_center, pT_bin_center, P_Lambda_lab_mean.X()/N_bins_rap);
                hLambdaPolY_phi_pTReco->Fill(phi_bin_center, pT_bin_center, P_Lambda_lab_mean.Y()/N_bins_rap);
                hLambdaPolZ_phi_pTReco->Fill(phi_bin_center, pT_bin_center, P_Lambda_lab_mean.Z()/N_bins_rap);

                // Finally, a 1D projection for the distribution in Phi only:
                hLambdaPolX_phiReco->Fill(phi_bin_center, P_Lambda_lab_mean.X()/(N_bins_pT*N_bins_rap));
                hLambdaPolY_phiReco->Fill(phi_bin_center, P_Lambda_lab_mean.Y()/(N_bins_pT*N_bins_rap));
                hLambdaPolZ_phiReco->Fill(phi_bin_center, P_Lambda_lab_mean.Z()/(N_bins_pT*N_bins_rap));
            }
        }
    }

    // Averaging the true values too -- This should already mask the arrays too:
    hLambdaPolX_pT_y->Divide(hLambdaCounter_pT_y);
    hLambdaPolY_pT_y->Divide(hLambdaCounter_pT_y);
    hLambdaPolZ_pT_y->Divide(hLambdaCounter_pT_y);

    ////////////////////
    //// 4 - Exporting into a .root file for later review
    ////////////////////
    std::cout << "\nExporting results" << std::endl;
    // std::string filename = "/home/cicero/results/hydro_vorticity/LambdaPol_lumpy_events_60_70_" + std::to_string(N_events) + "ev.root";
    // std::string filename = "/home/cicero/results/hydro_vorticity/LambdaPol2_lumpy_events_60_70_" + std::to_string(N_events) + "ev.root"; // Having many different file names for comparison
    std::string filename = "/home/cicero/results/hydro_vorticity/LambdaPolResampled" + std::to_string(N_resamples) + "times_lumpy_events_60_70_" + std::to_string(N_events) + "ev.root";
    TFile f(filename.c_str(), "RECREATE");

    hLambdaCounter->Write();
    hLambdaPolMag->Write();
    hLambdaPolStarMagReco->Write();

    hLambdaPolX->Write();
    hLambdaPolStarXReco->Write();
    hLambdaPolY->Write();
    hLambdaPolStarYReco->Write();
    hLambdaPolZ->Write();
    hLambdaPolStarZReco->Write();
    
    hLambdaPolX_pT_y->Write();
    hLambdaPolStarX_pT_yReco->Write();
    hLambdaPolY_pT_y->Write();
    hLambdaPolStarY_pT_yReco->Write();
    hLambdaPolZ_pT_y->Write();
    hLambdaPolStarZ_pT_yReco->Write();

    hLambdaCounter_pT_y->Write();
    hLambdaAvgDotX_pT_y->Write();
    hLambdaAvgDotY_pT_y->Write();
    hLambdaAvgDotZ_pT_y->Write();

    hLambdaCounter_pT_y_phi->Write();
    hLambdaAvgDotX_pT_y_phi->Write();
    hLambdaAvgDotY_pT_y_phi->Write();
    hLambdaAvgDotZ_pT_y_phi->Write();

    hLambdaPolStarX_pT_y_phiReco->Write();
    hLambdaPolStarY_pT_y_phiReco->Write();
    hLambdaPolStarZ_pT_y_phiReco->Write();

    hLambdaPolX_pT_y_phiReco->Write();
    hLambdaPolY_pT_y_phiReco->Write();
    hLambdaPolZ_pT_y_phiReco->Write();

    hLambdaPolX_pT_yReco->Write();
    hLambdaPolY_pT_yReco->Write();
    hLambdaPolZ_pT_yReco->Write();

    hLambdaPolX_phi_pTReco->Write();
    hLambdaPolY_phi_pTReco->Write();
    hLambdaPolZ_phi_pTReco->Write();

    hLambdaPolX_phiReco->Write();
    hLambdaPolY_phiReco->Write();
    hLambdaPolZ_phiReco->Write();

    f.Close();
    std::cout << "\tDone!" << std::endl;

    std::cout << "\n\nCode execution finished with code 0. Thank you!" << std::endl;
    return 0;
}

///////////////////////////////////
/// Defining sampling functions ///
///////////////////////////////////
// Declared the function as an inline-able object, so that sampling can take as short a time as possible:
inline double P_angle_proton(double xi_star, double P_Lambda_star){
    return (1.0 / PI) * (1.0 + alpha_H * P_Lambda_star * std::cos(xi_star));
}

// Now the actual sampler function -- The goal is to sample an angle given a P(angle) function.
    // Will use a rejection sampling for quick samples!
    // Receives the Lambda polarization magnitude, already in the rest frame, and returns the two angles needed 
    // Usage example:
    // auto [xi_star, phi_star] = sample_P_angle_proton_for_each_PLambda(...);
std::pair<double, double> sample_P_angle_proton(double P_Lambda_star_mag, std::mt19937 rng, std::uniform_real_distribution<double> dist_x,
                                                std::uniform_real_distribution<double> dist_y, std::uniform_real_distribution<double> dist_azimuth){
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
        // Rewriting to use pre-calculated mass^2 values:
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
    TVector3 P_Lambda_star_unit = P_Lambda_star.Unit(); // The reference for our sampling

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
                DoubleMatrix &pz_matrix, DoubleMatrix &E_matrix, DoubleMatrix& mult_matrix, DoubleMatrix &St_matrix,
                DoubleMatrix &Sx_matrix, DoubleMatrix &Sy_matrix, DoubleMatrix &Sz_matrix, int N_events){
    std::cout << "\nStarting event loop for Lambda getter..." << std::endl;
    for(int i = 0; i < N_events; ++i){
        if (i % int (N_events * 0.1) == 0){ // Keeping track of the proccess for every 10% of the events
            std::cout << "Now on event " << i << " of " << N_events  << " (" << i * 1./N_events * 100 << "%)" << std::endl; // The 1. comes to represent the percentage as a double
        }

        std::vector<double> y, phi, pT, mult, St, Sx, Sy, Sz; // Declaring variables that will store the Lambda information for the current event
        double const mass = 1.115683; // The Lambda hyperon mass
		
		TString input_file_path = "/storage1/vribeiro/lumpy_events_60_70/bullet/" + std::to_string(i)  + "/results/Smu_dpTdphidy_Thermal_rapidity_3122.dat";

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
        E_matrix.push_back(E);

        St_matrix.push_back(St);
        Sx_matrix.push_back(Sx);
        Sy_matrix.push_back(Sy);
        Sz_matrix.push_back(Sz);
    }
    // std::cout << "\nCompleted!" << std::endl;
}

// This function retrieves the same set of data that Vitor used in Pol_Analysis_Random_hist_ebe.C
// It modifies the provided vectors and stores the jet information.
void get_jet(std::vector<double> &n_event, std::vector<double> &n_random, std::vector<double> &phi_random,
             std::vector<double> &momentum_x, std::vector<double> &momentum_y){
    
    // Getting a file that has jets in hydro set in a random direction for the jet momentum deposition:
    const char* file_random_jet_path = "/storage1/vribeiro/lumpy_events_60_70/random_parameters.dat";

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