#include "TROOT.h"
#include "TFile.h"
// #include "TTree.h"

#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TProfile.h"

#include <vector>
#include <array>
#include <cmath>
#include <random>
#include "math.h"
#include "TString.h"
#include "TObjArray.h"

#include <dirent.h>
#include <iostream>
#include <string>
#include <fstream>

typedef std::vector<std::vector<double>> DoubleMatrix; // An alias

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
// Constant from Phys. Rev. D 110, 030001 - Published 1 August, 2024, i.e.: https://journals.aps.org/prd/pdf/10.1103/PhysRevD.110.030001
    // Declared both the constant and the function as inline-able objects, so that sampling can take as short a time as possible.
constexpr double PI = 3.14159265358979323846;
inline double P_angle_proton(double xi_star, double alpha_H, double P_Lambda) {
    return (1.0 / PI) * (1.0 + alpha_H * P_Lambda * std::cos(xi_star));
}

// Now the actual sampler function -- The goal is to sample an angle given a P(angle) function.
    // Will use a rejection sampling for quick samples!
    // Will store the samples inside the samples vector.
void sample_P_angle_proton(std::vector<double> &samples, double alpha_H, double P_Lambda) {
    std::mt19937 rng {std::random_device{}()}; // Creates a random_device object (temporary) with the uniform initialization on newer C++ (the "{}"),
                                               // then calls it with () getting a random seed, and then that seed is passed to the random engine.

    // Uniform distribution for angle xi_star between 0 and pi
    std::uniform_real_distribution<double> dist_x(0.0, PI);

    // We compute P_max (the maximum possible value of P(x))
    // This occurs at x = 0 or x = pi depending on sign of cos(x) term
    double P_max = (1.0 / PI) * (1.0 + std::abs(alpha_H * P_Lambda));

    // Uniform distribution for y between 0 and P_max
    // This is used for the rejection sampling step
    std::uniform_real_distribution<double> dist_y(0.0, P_max);

    size_t N_part = samples.size();  // Number of particles (or samples) we want to generate
    size_t i = 0;                    // Counter for accepted samples

    // Rejection sampling loop
    while (i < N_part) {
        // Propose an x uniformly in [0, pi]
        double x = dist_x(rng);

        // Propose a y uniformly in [0, P_max]
        double y = dist_y(rng); // Draws from the maximum value distribution independently.
                                // This makes it so that we always follow the appropriate limit and have correctly weighted samples of x(!)

        // Compute the actual probability density at x
        double px = P_angle_proton(x, alpha_H, P_Lambda);

        // Accept x if y < P(x). This ensures samples are distributed according to P(x): larger P(x) values will be accepted more easily!
        if (y < px) {
            samples[i++] = x;  // Store accepted x and increment sample counter
        }
        // Otherwise, reject and try again
    }
}

/////////////////////////////////
/// Defining getter functions ///
/////////////////////////////////
// Receives 2D vectors (matrices) that have the event index in the first dimension and the particle index in the second.
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