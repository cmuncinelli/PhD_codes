#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"

#include "TH2.h"
#include "TH2D.h"
#include "TMath.h"

#include <dirent.h>
#include <filesystem> // Useful for manipulation of the input_card_path

#include <iostream>

void multiplicity_to_centrality(TH1D *multiplicity_hist, TH1D *centrality_hist);
void doCentrality(std::string output_folder, std::string input_card_name);

int main(int argc, char *argv[]){
    if (argc != 3){
		std::cerr << "Usage:" << argv[0] << " input_card_name N_events" << std::endl;
        return 1; // Indicate an error.
    }

	const char *input_card_name = (const char*) argv[1];
    double N_ev_receiver = std::atof(argv[2]); // This receives input like 1e9 and converts it into a proper double.
    // uint N_events = static_cast<int>(N_ev_receiver); // Should use long long for real high statistics (even uint would overflow at 1e10 or so!), but this is enough for now
	long long N_events = static_cast<long long>(N_ev_receiver);  // Use long long because int fails when N_ev ~ 1e10 !

    std::cout << "\nCorrecting overflowed centrality estimation for input card " << input_card_name << " and running for " << N_events << " events" << std::endl;

    // Defining the path related variables    
    std::string input_folder = "/storage2/cicero/pythia_data/";
    // std::string input_folder = "/storage3/cicero/pythia_data/"; // Quick modification for jarvis5 (storage3 was freer here)
    input_folder += input_card_name;
    input_folder += "_" + std::to_string(N_events) + "ev/"; // Specific folder for this input card and its number of events!

    std::string output_folder = input_folder;

    // Done in the doCentrality part:
    // std::string centrality_conversion_file_path = input_folder + "Ntracks_to_centrality_conversion.root";

    // Now re-running the centrality conversion:
        // Processing centrality for all files:
	doCentrality(output_folder, input_card_name);
    std::cout << "Done! Thank you for waiting!" << std::endl;

    return 1;
}



// Adapted from PythiaGenMin, but instead takes the pre-summed multiplicity histograms and just iterates on them to calculate centrality:
void doCentrality(std::string output_folder, std::string input_card_name){
	std::cout << "\nProcessing centrality" << std::endl;
	int multiplicity_nbins = 50000; // The number of bins is also the number of entries, which gives one bin per possible multiplicity!

    // Getting the input data:    
    std::string output_centrality_filename =(std::string) output_folder;
    std::string input_centrality_filename = (std::string) output_folder;

    input_centrality_filename += "/Ntracks_to_centrality_conversion-Overflows.root";
	output_centrality_filename += "/Ntracks_to_centrality_conversion.root";

	std::cout << "Reading input file in path " << input_centrality_filename << std::endl;
    TFile *data_file = TFile::Open(input_centrality_filename.c_str(), "READ");
    TFile output_file(output_centrality_filename.c_str(), "RECREATE");

	// Reading all N_cores worker data and summing the hNtracks histograms:
		// First declaring histograms to receive the sums:
    std::cout << "Collecting histograms" << std::endl;
	TH1D *hNtracks_sum = dynamic_cast<TH1D*>(data_file->Get("hNtracks_sum"));
	TH1D *hNtracks_charged_forward_sum = dynamic_cast<TH1D*>(data_file->Get("hNtracks_charged_forward_sum"));
	TH1D *hNtracks_primary_sum = dynamic_cast<TH1D*>(data_file->Get("hNtracks_primary_sum"));
	TH1D *hNtracks_final_sum = dynamic_cast<TH1D*>(data_file->Get("hNtracks_final_sum"));
	TH1D *hNtracks_final_forward_sum = dynamic_cast<TH1D*>(data_file->Get("hNtracks_final_forward_sum"));
	TH1D *hNtracks_final_forward_INEL_l0_sum = dynamic_cast<TH1D*>(data_file->Get("hNtracks_final_forward_INEL_l0_sum"));
	TH1D *hNtracks_final_charged_center_sum = dynamic_cast<TH1D*>(data_file->Get("hNtracks_final_charged_center_sum"));
	TH1D *hNtracks_final_charged_sum = dynamic_cast<TH1D*>(data_file->Get("hNtracks_final_charged_sum"));
	TH1D *hNtracks_final_charged_INEL_l0_sum = dynamic_cast<TH1D*>(data_file->Get("hNtracks_final_charged_INEL_l0_sum"));
	TH1D *hNtracks_final_charged_forward_sum = dynamic_cast<TH1D*>(data_file->Get("hNtracks_final_charged_forward_sum"));
	TH1D *hNtracks_final_charged_forward_INEL_l0_sum = dynamic_cast<TH1D*>(data_file->Get("hNtracks_final_charged_forward_INEL_l0_sum"));

    // data_file->Close(); // This would make the histograms be deleted, as they are not in the gDirectory of the output_file! Don't call this!
    std::cout << "Starting centrality definition process" << std::endl;
	output_file.cd();

	// Processing these multiplicity histograms into centrality histograms:
		// First, declaring the new histograms using the information of maximum and minimum multiplicity of the TH1D's:
		// The number of bins is exactly the maximum number of multiplicity + 1 (the upper limit is non-inclusive!).
		// Index 0 will be multiplicity [0, 1), index 1 will be [1, 2) and so on until the last index.
		// If the maximum value (the last one) was 10, then we need 11 bins to get the [10, 11) bin as a regular bin, not as the overflow bin!
	// TH1D *hNtracks_to_centrality = new TH1D("hNtracks_to_centrality", "hNtracks_to_centrality", hNtracks->GetMaximum() + 1, 0, hNtracks->GetMaximum() + 1);
		// Updated to use FindLastBinAbove() instead of maximum value! I don't want the highest value in the TH1D, I want the highest bin with non-zero value!
	TH1D *hNtracks_to_centrality = new TH1D("hNtracks_to_centrality", "hNtracks_to_centrality", hNtracks_sum->FindLastBinAbove(0) + 1, 0, hNtracks_sum->FindLastBinAbove(0) + 1);
	TH1D *hNtracks_to_centrality_charged_forward = new TH1D("hNtracks_to_centrality_charged_forward", "hNtracks_to_centrality_charged_forward", hNtracks_charged_forward_sum->FindLastBinAbove(0) + 1, 0, hNtracks_charged_forward_sum->FindLastBinAbove(0) + 1);
	TH1D *hNtracks_to_centrality_primary = new TH1D("hNtracks_to_centrality_primary", "hNtracks_to_centrality_primary", hNtracks_primary_sum->FindLastBinAbove(0) + 1, 0, hNtracks_primary_sum->FindLastBinAbove(0) + 1);
	TH1D *hNtracks_to_centrality_final = new TH1D("hNtracks_to_centrality_final", "hNtracks_to_centrality_final", hNtracks_final_sum->FindLastBinAbove(0) + 1, 0, hNtracks_final_sum->FindLastBinAbove(0) + 1);
	TH1D *hNtracks_to_centrality_final_forward = new TH1D("hNtracks_to_centrality_final_forward", "hNtracks_to_centrality_final_forward", hNtracks_final_forward_sum->FindLastBinAbove(0) + 1, 0, hNtracks_final_forward_sum->FindLastBinAbove(0) + 1);
	TH1D *hNtracks_to_centrality_final_forward_INEL_l0 = new TH1D("hNtracks_to_centrality_final_forward_INEL_l0", "hNtracks_to_centrality_final_forward_INEL_l0", hNtracks_final_forward_INEL_l0_sum->FindLastBinAbove(0) + 1, 0, hNtracks_final_forward_INEL_l0_sum->FindLastBinAbove(0) + 1);
	TH1D *hNtracks_to_centrality_final_charged_center = new TH1D("hNtracks_to_centrality_final_charged_center", "hNtracks_to_centrality_final_charged_center", hNtracks_final_charged_center_sum->FindLastBinAbove(0) + 1, 0, hNtracks_final_charged_center_sum->FindLastBinAbove(0) + 1);
	TH1D *hNtracks_to_centrality_final_charged = new TH1D("hNtracks_to_centrality_final_charged", "hNtracks_to_centrality_final_charged", hNtracks_final_charged_sum->FindLastBinAbove(0) + 1, 0, hNtracks_final_charged_sum->FindLastBinAbove(0) + 1);
	TH1D *hNtracks_to_centrality_final_charged_INEL_l0 = new TH1D("hNtracks_to_centrality_final_charged_INEL_l0", "hNtracks_to_centrality_final_charged_INEL_l0", hNtracks_final_charged_INEL_l0_sum->FindLastBinAbove(0) + 1, 0, hNtracks_final_charged_INEL_l0_sum->FindLastBinAbove(0) + 1);
	TH1D *hNtracks_to_centrality_final_charged_forward = new TH1D("hNtracks_to_centrality_final_charged_forward", "hNtracks_to_centrality_final_charged_forward", hNtracks_final_charged_forward_sum->FindLastBinAbove(0) + 1, 0, hNtracks_final_charged_forward_sum->FindLastBinAbove(0) + 1);
	TH1D *hNtracks_to_centrality_final_charged_forward_INEL_l0 = new TH1D("hNtracks_to_centrality_final_charged_forward_INEL_l0", "hNtracks_to_centrality_final_charged_forward_INEL_l0", hNtracks_final_charged_forward_INEL_l0_sum->FindLastBinAbove(0) + 1, 0, hNtracks_final_charged_forward_INEL_l0_sum->FindLastBinAbove(0) + 1);

    std::cout << "Calling multiplicity_to_centrality" << std::endl;
	multiplicity_to_centrality(hNtracks_sum, hNtracks_to_centrality);
	multiplicity_to_centrality(hNtracks_charged_forward_sum, hNtracks_to_centrality_charged_forward);
	multiplicity_to_centrality(hNtracks_primary_sum, hNtracks_to_centrality_primary);
	multiplicity_to_centrality(hNtracks_final_sum, hNtracks_to_centrality_final);
	multiplicity_to_centrality(hNtracks_final_forward_sum, hNtracks_to_centrality_final_forward);
	multiplicity_to_centrality(hNtracks_final_forward_INEL_l0_sum, hNtracks_to_centrality_final_forward_INEL_l0);
	multiplicity_to_centrality(hNtracks_final_charged_center_sum, hNtracks_to_centrality_final_charged_center);
	multiplicity_to_centrality(hNtracks_final_charged_sum, hNtracks_to_centrality_final_charged);
	multiplicity_to_centrality(hNtracks_final_charged_INEL_l0_sum, hNtracks_to_centrality_final_charged_INEL_l0);
	multiplicity_to_centrality(hNtracks_final_charged_forward_sum, hNtracks_to_centrality_final_charged_forward);
	multiplicity_to_centrality(hNtracks_final_charged_forward_INEL_l0_sum, hNtracks_to_centrality_final_charged_forward_INEL_l0);

	std::cout << "Writing to output file" << std::endl;
    hNtracks_sum->Write();
	hNtracks_charged_forward_sum->Write();
	hNtracks_primary_sum->Write();
	hNtracks_final_sum->Write();
	hNtracks_final_forward_sum->Write();
	hNtracks_final_forward_INEL_l0_sum->Write();
	hNtracks_final_charged_center_sum->Write();
	hNtracks_final_charged_sum->Write();
	hNtracks_final_charged_INEL_l0_sum->Write();
	hNtracks_final_charged_forward_sum->Write();
	hNtracks_final_charged_forward_INEL_l0_sum->Write();

	hNtracks_to_centrality->Write();
	hNtracks_to_centrality_charged_forward->Write();
	hNtracks_to_centrality_primary->Write();
	hNtracks_to_centrality_final->Write();
	hNtracks_to_centrality_final_forward->Write();
	hNtracks_to_centrality_final_forward_INEL_l0->Write();
	hNtracks_to_centrality_final_charged_center->Write();
	hNtracks_to_centrality_final_charged->Write();
	hNtracks_to_centrality_final_charged_INEL_l0->Write();
	hNtracks_to_centrality_final_charged_forward->Write();
	hNtracks_to_centrality_final_charged_forward_INEL_l0->Write();

	return;
}



// A corrected version, that should deal with the high number of events using long long
// Converting the event multiplicity to centrality in the most precise way possible -- Corresponding each integer number of particles to a limit in centrality.
void multiplicity_to_centrality(TH1D *multiplicity_hist, TH1D *centrality_hist){
	long long Nev = multiplicity_hist->GetEntries();
	int Nbins = centrality_hist->GetNbinsX(); // I want to loop only on the [0 particles, maximum_number_of_particles] region, not in the [0, 10.000] particles region that the hNtracks histograms need to have to avoid missing any data.

	long long number_of_events_up_to_current_bin = 0; // The number of events that came up to this bin, including this bin.
	for (int bin = 1; bin <= Nbins; bin++){ // 1 <= bin <= Nbins because the underflow bin is 0 and the overflow bin is Nbins + 1.
		// Determining the height of the current bin, i.e., the lower centrality limit that this bin has.
		number_of_events_up_to_current_bin += multiplicity_hist->GetBinContent(bin); // Sums the contents of the current bin to see how much of the total number of events we have gone through
		double centrality_inferior_limit = (1 - number_of_events_up_to_current_bin * 1./Nev) * 100; // The 1. in the division is to convert the rightmost term to a double.
		
		centrality_hist->SetBinContent(bin, centrality_inferior_limit);
	}

	// Setting the "Entries" value manually for the statistics box to show the number of events involved in the calculation.
	centrality_hist->SetEntries(Nev);
	centrality_hist->GetXaxis()->SetTitle("Number of particles of current estimator");
	centrality_hist->GetYaxis()->SetTitle("Lower cent_lim for current bin of Npart");
}