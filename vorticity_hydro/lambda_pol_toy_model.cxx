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

#include <vector>
#include <array>
#include <cmath>
#include <random>
// #include "math.h" // Same as cmath!
#include "TString.h"
#include "TObjArray.h"

#include <dirent.h>
#include <iostream>
#include <string>
#include <fstream>

typedef std::vector<std::vector<double>> DoubleMatrix; // An alias

// Some constants for the code:
const double PI = 3.14159265358979323846;
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
void get_lambda(DoubleMatrix &y_matrix, DoubleMatrix &phi_matrix, DoubleMatrix &px_matrix, DoubleMatrix &py_matrix,
                DoubleMatrix &pz_matrix, DoubleMatrix &E_matrix, DoubleMatrix &St_matrix, DoubleMatrix &Sx_matrix,
                DoubleMatrix &Sy_matrix, DoubleMatrix &Sz_matrix, int N_events);
void get_jet(std::vector<double> &n_event, std::vector<double> &n_random, std::vector<double> &phi_random,
             std::vector<double> &momentum_x, std::vector<double> &momentum_y);

// This code assumes we are in Jarvis4! Also, the getter functions are based on Pol_Analysis_Random_hist_ebe.C from the HadrEx_Ph repository
int lambda_pol_toy_model(){
    int N_events = 250; // There are only 250 events in the desired folder. Possibly oversampled to get statistics.

    DoubleMatrix y_matrix, phi_matrix, px_matrix, py_matrix, pz_matrix, E_matrix, mult_matri, St_matrix, Sx_matrix, Sy_matrix, Sz_matrix;
    get_lambda(y_matrix, phi_matrix, px_matrix, py_matrix, pz_matrix, E_matrix, mult_matrix, St_matrix, Sx_matrix, Sy_matrix, Sz_matrix);

    return 0
}

///////////////////////////////////
/// Defining sampling functions ///
///////////////////////////////////
// Declared the function as an inline-able object, so that sampling can take as short a time as possible:
inline double P_angle_proton(double xi_star, double P_Lambda){
    return (1.0 / PI) * (1.0 + alpha_H * P_Lambda * std::cos(xi_star));
}

// Now the actual sampler function -- The goal is to sample an angle given a P(angle) function.
    // Will use a rejection sampling for quick samples!
    // Will pass P_Lambda_vec, a vector with many different Lambda polarizations, and retrieve a samples vector from it.
    // Usage example:
    // std::vector<double> P_Lambda_vec = {0.2, 0.5, -0.3, 0.0, 1.0};  // One P_Lambda per particle
    // auto [xi_star_vec, phi_vec] = sample_P_angle_proton_for_each_PLambda(P_Lambda_vec);
std::pair<std::vector<double>, std::vector<double>> sample_P_angle_proton(const std::vector<double> &P_Lambda_vec){
    std::vector<double> samples;
    samples.reserve(P_Lambda_vec.size());  // Reserve space for performance
    samples_azimuth.reserve(P_Lambda_vec.size()); // The azimutal angle from 0 to 2 pi

    std::mt19937 rng {std::random_device{}()};  // Creates a random_device object (temporary) with the uniform initialization on newer C++ (the "{}"),
                                                // then calls it with () getting a random seed, and then that seed is passed to the random engine.

    // Uniform distribution for angle xi_star between 0 and pi
    std::uniform_real_distribution<double> dist_x(0.0, PI);
    std::uniform_real_distribution<double> dist_azimuth(0.0, 2*PI);

    // Loop over all provided P_Lambda values
    for (double P_Lambda : P_Lambda_vec){
        // Compute P_max (the maximum possible value of P(x))
        // This occurs at x = 0 or x = pi depending on sign of cos(x) term
        double P_max = (1.0 / PI) * (1.0 + std::abs(alpha_H * P_Lambda));

        // Uniform distribution for y between 0 and P_max
        // This is used for the rejection sampling step
        std::uniform_real_distribution<double> dist_y(0.0, P_max);

        // Rejection sampling loop - generate one sample for this P_Lambda
        while (true){
            // Propose an x uniformly in [0, pi]
            double x = dist_x(rng);

            // Propose a y uniformly in [0, P_max]
            double y = dist_y(rng);  // Draws from the maximum value distribution independently.
                                     // This makes it so that we always follow the appropriate limit and have correctly weighted samples of x(!)

            // Compute the actual probability density at x
            double px = P_angle_proton(x, P_Lambda);

            // Accept x if y < P(x). This ensures samples are distributed according to P(x): larger P(x) values will be accepted more easily!
            if (y < px){
                samples.push_back(x);  // Store accepted x
                break;                // Done with this P_Lambda, move to next
            }
            // Otherwise, reject and try again
        }

        // Now sampling the azimuthal angle:
        double phi = dist_azimuth(rng);
        samples_azimuth.push_back(phi);
    }
    return std::make_pair(samples, samples_azimuth);
}


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

TLorentzVector Lambda_decay(TLorentzVector Lambda_4_momentum_lab, TVector3 P_Lambda_lab, double xi_star, double phi_star){
    ////////////////////
    //// 1 - Boosting the polarization to the lambda rest frame:
    ////////////////////
    TVector3 beta_inverse = -Lambda_4_momentum_lab.BoostVector();
    TLorentzVector P_Lambda_lab_4vec(P_Lambda_lab.X(), P_Lambda_lab.Y(), P_Lambda_lab.Z(), 0); // TLorentzVectors have boost properties!

    TLorentzVector P_Lambda_star_4vec = P_Lambda_lab_4vec;
    P_Lambda_star_4vec.Boost(beta_inverse); // This boosts the P_Lambda_star_4vec into the rest frame.

    // Getting just the 3 needed components of this vector, in a shorter form for the function.
        // Also, only the spatial part is what matters in this case
    TVector3 P_Lambda_star = P_Lambda_star_4vec.Vect();
    
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

    // First, calculating the momenta in a frame where the polarization is in the z axis (easier to define)
    double px_star = daughter_momentum * std::sin(xi_star) * std::cos(phi_star);
    double py_star = daughter_momentum * std::sin(xi_star) * std::sin(phi_star);
    double pz_star = daughter_momentum * std::cos(xi_star); // This means the z axis points in the polarization direction! Try xi_star = 0 and xi_star = pi to convince yourself.
    // TLorentzVector proton_4_momentum(px_star, py_star, pz_star, E_p_star);

    ////////////////////
    //// 4 - Now rotating this frame into the actual xyz coordinates where the polarization points
    ////////////////////
    P_Lambda_star_unit = P_Lambda_star.Unit(); // The direction of the polarization, which we defined to be in the z axis for the (px,py,pz) vector above.

    // Building a rotation matrix around that will take this local z axis definition into the direction of the polarization in the usual XYZ coordinates:
    TVector3 z_axis(0, 0, 1);
    TVector3 rotation_axis = z_axis.Cross(P_Lambda_star_unit);
    double rotation_angle = std::acos(z_axis.Dot(P_Lambda_star_unit));
    TMatrix3 rotation_matrix;
    if (rotation_axis.Mag() > 1e-6){
        rotation_axis = rotation_axis.Unit();
        rotation_matrix.SetRotation(rotation_axis, rotation_angle);
    }
    else{
        // Already aligned (or opposite) - either no rotation or 180 deg
        if (z_axis.Dot(P_Lambda_star_unit) < 0) {
            // 180 deg rotation about X (could be any axis perpendicular to z)
            rotation_matrix.SetRotationX(PI);
        } else {
            // No rotation
            rotation_matrix.UnitMatrix();
        }
    }

    TVector3 p_star_vec(px_star, py_star, pz_star);
    TVector3 p_star_rotated = rotation_matrix * p_star_vec;

    ////////////////////
    //// 5 - Building the proton 4-momentum and boosting it into the lab frame:
    ////////////////////
    TLorentzVector proton_4_momentum(p_star_rotated, E_p_star); // Short for (p_star_rotated.X(), p_star_rotated.Y(), p_star_rotated.Z(), E_p_star);
        // Modifying the vector in place with the Boost() method:
    proton_4_momentum.Boost(-beta_inverse); // -beta_inverse is just beta, so this properly translates to the lab frame!
    
    return proton_4_momentum;
}






/////////////////////////////////
/// Defining getter functions ///
/////////////////////////////////
// Receives 2D vectors (matrices) that have the event index in the first dimension and the particle index in the second.
// Retrieves quantities in the lab frame!
// If you'd like to read the data Vitor processed, just look at files like these:
// TString output_ebe = "/storage1/vribeiro/work_station/analysis/Polarization/data_ebe/60_70_bullet/event"+std::to_string(i)+".root";
void get_lambda(DoubleMatrix &y_matrix, DoubleMatrix &phi_matrix, DoubleMatrix &px_matrix, DoubleMatrix &py_matrix,
                DoubleMatrix &pz_matrix, DoubleMatrix &E_matrix, DoubleMatrix &St_matrix, DoubleMatrix &Sx_matrix,
                DoubleMatrix &Sy_matrix, DoubleMatrix &Sz_matrix, int N_events){
    std::cout << "\nStarting event loop..." << std::endl;
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
    std::cout << "\nCompleted!" << std::endl;
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
    TObjArray* random_lines = random_data.Tokenize("\n");
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