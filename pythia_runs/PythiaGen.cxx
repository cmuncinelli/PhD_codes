// Author: Gianni S. S. Liveraro, adapted by Cicero D. Muncinelli
// Date: 10/09/24, adapted on 12/09/24 up to March/2025
// #include "/home/daviddc/pythia8/include/Pythia8/Pythia.h"
#include "/home/cicero/pythia8_3/include/Pythia8/Pythia.h" // Updated from the #include "/home/daviddc/pythia8/include/Pythia8/Pythia.h" installation from prof. David
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH1I.h"

#include <dirent.h>

#include <iostream>

// void histogram_copy(TH1D *origin_hist, TH1D *destiny_hist, int N_bins);
void multiplicity_to_centrality(TH1I *multiplicity_hist, TH1D *centrality_hist);

using namespace Pythia8;

int main(int argc, char *argv[]){
    // (const char *output_folder, int hard_bool, int N_ev, double target_rapidity)

    if (argc != 6) {
        // std::cerr << "Usage:" << argv[0] << " output_folder hard_bool N_ev target_rapidity s_nn pseudorap_int input_card_path" << std::endl;
		// std::cerr << "Usage:" << argv[0] << " output_folder card_input_folder simulation_idx s_nn target_rapidity pseudorap_int N_ev" << std::endl;
		std::cerr << "Usage:" << argv[0] << " output_folder card_input_folder simulation_idx s_nn N_ev" << std::endl;
        return 1; // Indicate an error. (This warning code block is chat-gpt made, by the way!)
    }

    const char *output_folder = (const char*) argv[1]; // Now this is the folder with all the runs in subfolders
	const char *card_input_folder = (const char*) argv[2];
	int simulation_idx = atoi(argv[3]);
	int s_nn = atoi(argv[4]); // Can receive this as int, as the units will be in GeV, not TeV
    // double target_rapidity = std::stod((std::string) ((const char*) argv[5])); // Conversion to double is trickier...
	// int pseudorap_int = atoi(argv[6]);
	// bool experimental_intervals_bool = (bool) atoi(argv[6]); // Converts 0 or 1 to a bool
	int N_ev = atoi(argv[5]);


	std::string folder_path = (std::string) card_input_folder + "/"; // Had to include this extra "/" in order for it to work.
    DIR *dir;
    struct dirent *entry;


	std::string input_card_path;
	if ((dir = opendir(folder_path.c_str())) != nullptr){
        std::cout << "Starting the loop on the input cards folder: " << folder_path << std::endl;

        while ((entry = readdir(dir)) != nullptr){
            std::string file_name = entry->d_name;
            if (string(entry->d_name) == "." || string(entry->d_name) == ".."){
                continue; // You have to skip these weird filenames that dirent will give you...
            }
            std::string file_path = folder_path + file_name; // Converting it into a string

			size_t output_index = file_path.find(std::string("pythia_card") + std::to_string(simulation_idx) + "_"); // The cards can't have a preceding 0 on numbers such as 06. It will not be read appropriatelly.
			if (output_index != string::npos){ // aka if the current file's simulation_idx is the one being searched on the loop
				input_card_path = file_path;
			}
		}
	}


	

  	Pythia pythia;
	// pythia.readString("Beams:idA = 11"); // This is for an electron-positron collision
	// pythia.readString("Beams:idB = -11");
	// pythia.readString("Beams:eCM = 209.");
	// pythia.readString("Next:numberShowEvent = 1");

	pythia.readString("Beams:idA = 2212");
	// pythia.readString("Beams:idB = -2212");
	pythia.readString("Beams:idB = 2212");
	std::string energy_string = "Beams:eCM = " + std::to_string(s_nn) + ".";
	// pythia.readString("Beams:eCM = 7000."); // Running at s_nn GeV
	pythia.readString(energy_string.c_str()); // Running at s_nn GeV
	
	// Some settings related to output in init(), next() and stat() that came from the Parnassus settings
	pythia.readString("Init:showChangedSettings = on"); // list changed settings
	pythia.readString("Init:showChangedParticleData = off"); // list changed particle data
	pythia.readString("Next:numberCount = 0"); // print message every n events
	pythia.readString("Next:numberShowInfo = 0"); // print event information n times
	pythia.readString("Next:numberShowProcess = 1"); // print process record n times
	pythia.readString("Next:numberShowEvent = 1"); // print event record n times

	// 2 - Seeds:
	pythia.readString("Random:setSeed = on"); 
	pythia.readString("Random:seed = 0"); // For random numbers on each run
	// pythia.readString("Random:seed = 42"); // Kept for testing!

	// 3 - Tunes settings:
	// pythia.readString("Tune:preferLHAPDF = 2"); // Disabled to use PYTHIA 8.313. In David's code, this used to work... maybe it is internal in PYTHIA 8.313 now?
	// Welp, that should only be necessary if LHAPDF is kept up-to-date in David's files, which it surely wasn't!
	pythia.readString("Tune:pp = 14");

	// 4 - Event settings:
	std::cout << "Reading from input file " << input_card_path << std::endl;
	pythia.readFile(input_card_path);

	// pythia.readString("");

	pythia.init();


	const Int_t kMaxTrack = 7000;
	Int_t ntrack;
	Int_t ntrack_final;
	Int_t ID[kMaxTrack];
	Int_t mother1_ID[kMaxTrack];
	Int_t mother2_ID[kMaxTrack];
	Bool_t IsFinal[kMaxTrack];
	Bool_t IsPrimary[kMaxTrack];
	Bool_t IsCarbonCopy[kMaxTrack];
	Int_t charged_in_back_forward_eta; // The multiplicity indicator for 1807.11321 [nucl-ex]'s data
	Int_t ntracks_in_back_forward_eta; // Another multiplicity indicator, trying to mimic the indirect charged particle production due to neutral's interactions with V0. Something finer than just Nch, and will also receive only final particles
	Int_t Ntracks_charged_forward;
	Int_t Ntracks_primary;
	Int_t Ntracks_final_charged_center;
	Int_t Ntracks_final_charged;

	Float_t px[kMaxTrack];
	Float_t py[kMaxTrack];
	Float_t pz[kMaxTrack];
	Float_t pt[kMaxTrack];
	Float_t m[kMaxTrack];
	Float_t e[kMaxTrack];

	Bool_t charged[kMaxTrack];
	Bool_t charged_in_central_eta; // A single boolean that is reset at the start of each new event and marks that event as INEL>0 (or not).

	Float_t Eta[kMaxTrack];
	Float_t y[kMaxTrack];
	Float_t Phi[kMaxTrack];

	TH1D *hEventCounter = new TH1D ( "hEventCounter", "", 1, -1, 1);
	TH1D *hINELEventCounter = new TH1D ( "hINELEventCounter", "", 1, -1, 1);
	TH1D *hEventCounterCharged = new TH1D ( "hEventCounterCharged", "", 1, -1, 1);
	TH1D *hEventCounterPion = new TH1D ( "hEventCounterPion", "", 1, -1, 1);
	TH1D *hEventCounterProton = new TH1D ( "hEventCounterProton", "", 1, -1, 1);
	TH1D *hEventCounterKaon = new TH1D ( "hEventCounterKaon", "", 1, -1, 1);

	std::string filename = (std::string) output_folder;
	filename += "/Pythia_card" + std::to_string(simulation_idx) + "_";
	filename += std::to_string(N_ev) + "events_";
	// if(pseudorap_int == 1){
	// 	filename += "pseudorap_selected_";
	// }
	filename += std::to_string(s_nn) + "GeV.root";
	// These files DO NOT have anything related to the experimental_intervals_bool in their names as they only became more general with the new trees, with no damage to the previous structures of data

	TFile f(filename.c_str(), "RECREATE");
	TTree *t3 = new TTree("t3","Reconst ntuple");
	
	// TTree *event_pT_tree = new TTree("event_pT_tree","event_pT_tree"); // Careful! The name in the final file is the one inside the parenthesis!
	int N_bins = 800;
	double upper_limit = 50;
	// TH1D *pointer_to_hist_for_filling = new TH1D("event_hist", "event_hist", N_bins, 0, upper_limit); // Before it was set as 700 bins from 0 to 7
	// event_pT_tree->Branch("charged_branch", "TH1D", pointer_to_hist_for_filling, 32000, 0);
	// event_pT_tree->Branch("pion_branch", "TH1D", pointer_to_hist_for_filling, 32000, 0);
	// event_pT_tree->Branch("proton_branch", "TH1D", pointer_to_hist_for_filling, 32000, 0);
	// event_pT_tree->Branch("kaon_branch", "TH1D", pointer_to_hist_for_filling, 32000, 0);

	// TBranch *charged_branch = event_pT_tree->GetBranch("charged_branch");
	// TBranch *pion_branch = event_pT_tree->GetBranch("pion_branch");
	// TBranch *proton_branch = event_pT_tree->GetBranch("proton_branch");
	// TBranch *kaon_branch = event_pT_tree->GetBranch("kaon_branch");

	t3->Branch("ntrack",&ntrack,"ntrack/I");
	t3->Branch("ntrack_final",&ntrack_final,"ntrack_final/I");
	t3->Branch("IsFinal",&IsFinal,"IsFinal[ntrack]/O"); // Actually don't need to pass this array by reference, as it is a C++ array and it would be passed by ref by default, but whatever...
	t3->Branch("IsPrimary",&IsPrimary,"IsPrimary[ntrack]/O"); 
	// Do notice that the [ntrack] part is making it so that only the values from 0 up to ntrack - 1 will be stored in your TTree, so you don't need to reset IsPrimary for each new
	// event in the loop!
	t3->Branch("IsCarbonCopy",&IsCarbonCopy,"IsCarbonCopy[ntrack]/O"); 
	t3->Branch("charged_in_back_forward_eta",&charged_in_back_forward_eta,"charged_in_back_forward_eta/I");
	t3->Branch("ntracks_in_back_forward_eta",&ntracks_in_back_forward_eta,"ntracks_in_back_forward_eta/I");
	t3->Branch("charged_in_central_eta",&charged_in_central_eta,"charged_in_central_eta/O");
	t3->Branch("Ntracks_charged_forward",&Ntracks_charged_forward,"Ntracks_charged_forward/I");
	t3->Branch("Ntracks_primary",&Ntracks_primary,"Ntracks_primary/I");
	t3->Branch("Ntracks_final_charged_center",&Ntracks_final_charged_center,"Ntracks_final_charged_center/I");
	t3->Branch("Ntracks_final_charged",&Ntracks_final_charged,"Ntracks_final_charged/I");

	t3->Branch("ID",&ID,"ID[ntrack]/I"); // A PID identifier
	t3->Branch("Mother1_ID",&mother1_ID,"Mother1_ID[ntrack]/I");
	t3->Branch("Mother2_ID",&mother2_ID,"Mother2_ID[ntrack]/I");
	t3->Branch("px",px,"px[ntrack]/F");
	t3->Branch("py",py,"py[ntrack]/F");
	t3->Branch("pz",pz,"pz[ntrack]/F");
	t3->Branch("pt",pt,"pt[ntrack]/F");
	t3->Branch("m",m,"m[ntrack]/F");
	t3->Branch("e",e,"e[ntrack]/F");
	t3->Branch("charged",charged,"charged[ntrack]/O");

	t3->Branch("Eta",Eta,"Eta[ntrack]/F");
	t3->Branch("y",y,"y[ntrack]/F");
	t3->Branch("Phi",Phi,"Phi[ntrack]/F");

	// Saving a global TH1D for each particle, just to see the information:
	TH1D *total_hist_charged_final = new TH1D("total_hist_charged_final", "total_hist_charged_final", N_bins, 0, upper_limit);
	TH1D *total_hist_pion_final = new TH1D("total_hist_pion_final", "total_hist_pion_final", N_bins, 0, upper_limit);
	TH1D *total_hist_proton_final = new TH1D("total_hist_proton_final", "total_hist_proton_final", N_bins, 0, upper_limit);
	TH1D *total_hist_kaon_final = new TH1D("total_hist_kaon_final", "total_hist_kaon_final", N_bins, 0, upper_limit);

	// Saving multiplicity information in a TH1I (which consumes at most 40 kB with this many bins!) -- Faster centrality/multiplicity class fetching:
	int multiplicity_nbins = 10000; // The number of bins is also the number of entries, which gives one bin per possible multiplicity!
		// Could make this stop at 7000 bins, as the kMaxTrack is 7000 for my pythia simulations, but this format is more general and the disk size is ridiculously small...
		// By declaring in this way, GetBin(0) will get the bin of events with 0 particles in multiplicity. The value 10000 is not stored, and it will go to the overflow bin.
		// In other words, GetBin(N) will get you the exact bin with N particles! A convenient indexing, where [0] is the index of 0 particles and [1] for 1 particle!
	TH1I *hNtracks = new TH1I("hNtracks", "hNtracks", multiplicity_nbins, 0, multiplicity_nbins);
	TH1I *hNtracks_charged_forward = new TH1I("hNtracks_charged_forward", "hNtracks_charged_forward", multiplicity_nbins, 0, multiplicity_nbins);
	TH1I *hNtracks_primary = new TH1I("hNtracks_primary", "hNtracks_primary", multiplicity_nbins, 0, multiplicity_nbins);
	TH1I *hNtracks_final = new TH1I("hNtracks_final", "hNtracks_final", multiplicity_nbins, 0, multiplicity_nbins);
	TH1I *hNtracks_final_forward = new TH1I("hNtracks_final_forward", "hNtracks_final_forward", multiplicity_nbins, 0, multiplicity_nbins);
	TH1I *hNtracks_final_forward_INEL_l0 = new TH1I("hNtracks_final_forward_INEL_l0", "hNtracks_final_forward_INEL_l0", multiplicity_nbins, 0, multiplicity_nbins);
	TH1I *hNtracks_final_charged_center = new TH1I("hNtracks_final_charged_center", "hNtracks_final_charged_center", multiplicity_nbins, 0, multiplicity_nbins);
	TH1I *hNtracks_final_charged = new TH1I("hNtracks_final_charged", "hNtracks_final_charged", multiplicity_nbins, 0, multiplicity_nbins);
	TH1I *hNtracks_final_charged_INEL_l0 = new TH1I("hNtracks_final_charged_INEL_l0", "hNtracks_final_charged_INEL_l0", multiplicity_nbins, 0, multiplicity_nbins);
	TH1I *hNtracks_final_charged_forward = new TH1I("hNtracks_final_charged_forward", "hNtracks_final_charged_forward", multiplicity_nbins, 0, multiplicity_nbins);
	TH1I *hNtracks_final_charged_forward_INEL_l0 = new TH1I("hNtracks_final_charged_forward_INEL_l0", "hNtracks_final_charged_forward_INEL_l0", multiplicity_nbins, 0, multiplicity_nbins);
	// Caution! All the non INEL>0 (INEL_l0 or "larger than 0" in the above notation) may include multiplicity = 0 events, which will be assigned a position in the centrality determination
	// as the most peripheral events!
	// This makes those estimators not so trustable, as we wouldn't usually consider 0-particle events in the centrality estimator in V0 (Or so I think. It might be good to check this!)

		// Now storing information to recover the number of events in each centrality class:
		// You don't actually need this information for the root_to_csv_custom_PYTHIA.cxx to work, and it won't make it much faster,
		// but it is nice to be able to visualize it!
	TH1I *hINELev_Ntracks_final_charged_forward = new TH1I ("hINELev_Ntracks_final_charged_forward", "", multiplicity_nbins, 0, multiplicity_nbins);
	TH1I *hINELev_Ntracks_final_charged_forward_Charged = new TH1I ("hINELev_Ntracks_final_charged_forward_Charged", "", multiplicity_nbins, 0, multiplicity_nbins);
	TH1I *hINELev_Ntracks_final_charged_forward_Pion = new TH1I ("hINELev_Ntracks_final_charged_forward_Pion", "", multiplicity_nbins, 0, multiplicity_nbins);
	TH1I *hINELev_Ntracks_final_charged_forward_Proton = new TH1I ("hINELev_Ntracks_final_charged_forward_Proton", "", multiplicity_nbins, 0, multiplicity_nbins);
	TH1I *hINELev_Ntracks_final_charged_forward_Kaon = new TH1I ("hINELev_Ntracks_final_charged_forward_Kaon", "", multiplicity_nbins, 0, multiplicity_nbins);

	// TH1D *event_hist_charged = new TH1D("event_hist_charged", "event_hist_charged", N_bins, 0, upper_limit);
	// TH1D *event_hist_pion = new TH1D("event_hist_pion", "event_hist_pion", N_bins, 0, upper_limit);
	// TH1D *event_hist_proton = new TH1D("event_hist_proton", "event_hist_proton", N_bins, 0, upper_limit);
	// TH1D *event_hist_kaon = new TH1D("event_hist_kaon", "event_hist_kaon", N_bins, 0, upper_limit);

	// pointer_to_hist_for_filling->Sumw2(); // Starting an error structure for this histogram
	// event_hist_charged->Sumw2();
	// event_hist_pion->Sumw2();
	// event_hist_proton->Sumw2();
	// event_hist_kaon->Sumw2();

	// Event loop
	for (int iEvent=0; iEvent<N_ev; ++iEvent){
		pythia.next();
		ntrack = pythia.event.size(); // This could be a bit problematic to estimate multiplicity: this includes non-final particles!
		ntrack_final = 0; // Resetting for the next event

		hEventCounter->Fill(0);
		if (iEvent % int(0.05 * N_ev) == 0){
			std::cout << "Now on event " << iEvent << " of " << N_ev << " or " << iEvent * 1./N_ev * 100 << " %" << std::endl;
		}

		charged_in_central_eta = false; // Resetting the bool for the next event
		charged_in_back_forward_eta = 0; // Actually ""charged_in_back_forward_eta_FINAL" cause it only gets final particles
		ntracks_in_back_forward_eta = 0;

		// Some variables that will not be stored, but will be used inside the loop for counting particles -- These serve only for studying and are not necessarily useful for any analysis:
		Ntracks_charged_forward = 0;
		Ntracks_primary = 0;
		Ntracks_final_charged = 0;
		Ntracks_final_charged_center = 0;

		// int Ntracks = 0; // Already accounted for and stored inside the TTree (ntrack variable)
		// int Ntracks_final = 0; // Already accounted for and stored inside the TTree (ntrack_final variable)
		// int Ntracks_final_forward = 0; // Already accounted for and stored inside the TTree (ntracks_in_back_forward_eta variable)
		// int Ntracks_final_charged_forward = 0; // Already accounted for and stored inside the TTree (charged_in_back_forward_eta variable)

		int n_charged_event = 0;
		int n_pion_event = 0;
		int n_proton_event = 0;
		int n_kaon_event = 0;
		
		// Track loop
		for (int i = 0; i < pythia.event.size(); ++i){
			bool isfinal = pythia.event[i].isFinal();
				// Selects only final particles, not intermediate ones:
			// if (!isfinal) continue;
			// ntrack_final += 1;
			// Am not skipping the final particles anymore! The mothers can be useful for studying strangeness!

			// Checking if this particle has a carbon-copy daughter. If it does, then it isn't final, and it didn't decay: it just scattered. This is not what I want to look at!
				// From PYTHIA: "daughter1 = daughter2 > 0: the particle has a "carbon copy" as its sole daughter, but with changed momentum as a "recoil" effect".
			Int_t daughter1 = pythia.event[i].daughter1();
			Int_t daughter2 = pythia.event[i].daughter2();
			
			// if ((daughter1 == daughter2) && (daughter1 > 0)){continue;} // This means this particle should not be considered!
			// Be careful though. When you do this kind of particle skipping, you need to keep track of the PYTHIA particle indexes in your exported TTree.
			// If you don't, then you have to be extra careful: the PYTHIA particle index will not be the same as the index you use in your TTree for the particles in each event!
				// If you do that kind of skipping, you will not stop ROOT from storing this particle's information in the Tree. You will only corrupt its information!
				// ROOT will still save all particles with indexes from 0 to ntrack - 1! The thing is you are not updating this particle information by doing so!
				// This could be artificially increasing the number of final state particles in your final spectra!
			// Final particles will never enter this category, but non-final particles will. To actually avoid saving this kind of particle, you would need to have two indexes:
			// the first index will actually run through the event, and the second index will run through the vectors that will be exported to the TTree:
			
				// TODO: Implement the carbon-copy selection variable for studying non-final states!
			if ((daughter1 == daughter2) && (daughter1 > 0)){
				IsCarbonCopy[i] = true; // This will be a useful bool to consider, where you will skip the particle if it is simply a carbon-copy: it may not be a final state,
				// but a carbon-copy isn't even a physical particle, and it shouldn't be considered if you want to count the number of intermediate particles in your collision.
				// It is best not to do an if ((daughter1 == daughter2) && (daughter1 > 0)){continue} if you are going to save these particles in a TTree, because the value on
				// the i'th index of the vector that will be flushed into the TTree will still keep the previous iteration's particle data!
				
					// Now doing the trick to avoid saving these particles again:
				// second_idx -= 1; // This will start as ntrack = pythia.event.size(), but will decrease with the iterations. Then for accessing the vectors you could use something like pt[pythia.event.size() - second_idx] = pT;
				// // You will increase second_idx each time you restart the particle loop, thus the need to reduce it by one if you aren't actually going to consider that particle in your final state
				// ntrack -= 1; // You have to subtract this in order to save one less particle in the TTree. You could possibly use this as part of your second index
			}
			

			// if (isfinal){ntrack_final += 1;}

			Float_t rapidity = pythia.event[i].y();
			Float_t eta_rap = pythia.event[i].eta();
			Bool_t isCharged = pythia.event[i].isCharged();

			// Starting the selections according to ALICE cuts in 1807.11321 [nucl-ex]
				// 1) The original conditions are:
			if (isfinal){
				ntrack_final += 1; // Checking if final

				if (isCharged){
					Ntracks_final_charged += 1; // Checking final charged in all eta
				}

				if ((eta_rap > -3.7 && eta_rap < -1.7) || (eta_rap > 2.8 && eta_rap < 5.1)){
					ntracks_in_back_forward_eta += 1; // Checking final in backward/forward eta
					if (isCharged){ // Selecting only the charged particles, for an Nch metric of multiplicity
						charged_in_back_forward_eta += 1; // Selects all charged particles, regardless of being from secondary decays, as V0 would get them too
						// Checking final charged in backward/forward eta
					}
				}
				else if (isCharged && (eta_rap > -1.0 && eta_rap < 1.0)){
					Ntracks_final_charged_center += 1; // Checking final charged in central eta
					charged_in_central_eta = true;
				}
			}
			else{
				if (isCharged && ((eta_rap > -3.7 && eta_rap < -1.7) || (eta_rap > 2.8 && eta_rap < 5.1))){
					Ntracks_charged_forward += 1; // Checking charged and forward, but non-final
				}
			}
			// else if (isfinal && isCharged && (eta_rap > -1.0 && eta_rap < 1.0)){ // Used an "else" because they are complementary intervals in rapidity
			// // if (isCharged && (eta_rap > -1.0 && eta_rap < 1.0)){ // Replaced the "else if" with a single "if" to be able to use superposing intervals for my quick quality checks
			// 	charged_in_central_eta = true;
			// }

			// For the new formats, this kind of cut is done in the root_to_csv_custom_PYTHIA scripts, not in the generator-level scripts!
			// if(pseudorap_int == 1){
			// 	if (eta_rap < (-1) * target_rapidity || eta_rap > target_rapidity) continue;
			// }
			// else{ // Notice that this comes AFTER the multiplicity tracking in ntracks_in_back_forward_eta and charged_in_back_forward_eta, so it should NOT influence multiplicity calculations
			// 	if (rapidity < (-1) * target_rapidity || rapidity > target_rapidity) continue; // Should not use particles out of the |y| < target_rapidity range
			// }

			int particle_PID = pythia.event[i].id();
			ID[i] = particle_PID;
			int motherIdx1 = pythia.event[i].mother1(); // Doing this, I can access it only once and don't need to re-read this information for pions!
			int motherIdx2 = pythia.event[i].mother2();
			mother1_ID[i] = motherIdx1;
			mother2_ID[i] = motherIdx2;

			IsFinal[i]   = isfinal;
			px[i]   = pythia.event[i].px();
			py[i]   = pythia.event[i].py();
			pz[i]   = pythia.event[i].pz();
			
			Float_t pT = pythia.event[i].pT();
			pt[i]   = pT;
			m[i]   = pythia.event[i].m();
			e[i]   = pythia.event[i].e();
			
			Phi[i]   = pythia.event[i].phi();
			Eta[i]   = eta_rap;
			y[i]   = rapidity;

			charged[i] = isCharged;

			// Checking if the particle would be considered primary or not
				// The particle will be considered primary by default. If it has a mother with a long lifetime, then it will be considered secondary.
				// This allows for some contamination from possible secondaries that have long decay showers, though. My hypothesis is that this effect is too small to be significant.
				// Possible secondaries from 2 -> N processes and two partons becoming a single particle will be considered as primaries.
				// A little bit of a MC cheat is that particles whose lifetime is smaller than 10 mm/c but came straight from the collision will be considered primary, 
				// and in experimental data that would not be true.
			IsPrimary[i] = true;
			// First, checking if it came from an actual decay:
				// (I want to count all particles that didn't fall in both the loop's conditions with the Ntracks_primary, not the ones that checked one but didn't check the other if)
			if (motherIdx1 > 0 && motherIdx2 == 0){
				if (pythia.event[motherIdx1].tau0() >= 10.){ // Tau0 in mm/c. Reference is slide 112 of https://indico.cern.ch/event/666222/contributions/2768780/attachments/1551303/2437229/DPG_AnalysisTutorial_20171102.pdf
					IsPrimary[i] = false;
				}
				else{
					Ntracks_primary += 1;
				}
			}
			else{
				Ntracks_primary += 1; // This does not overcount the above check: they never happen simultaneuosly for the same particle!
			}

			if (isfinal && isCharged){
				// event_hist_charged->Fill((double) pT);
				total_hist_charged_final->Fill((double) pT);
				n_charged_event += 1;
			}

			if (isfinal && (particle_PID == 211 or particle_PID == (-1)*211)){ // Could just use idAbs() == 211 too!
				// Excluding all pions that came from weak decays, in a workaround to not having the ALIPhyisics' isSecondaryFromWeakDecay() in pure PYTHIA.
					// Specifically, excluding pions from mothers with Tau0 > 10 mm/c (could have used particle.tau0() to do a general check too!),
					// that are surely secondary!
					// Did not use the ParticleDecays:limitTau0 with tau0Max = 10 mm/c because it would influence the overall event multiplicity in a different
					// way from the experiment.
				// int motherIdx = pythia.event[i].mother1();
				if (motherIdx1 > 0 && motherIdx2 == 0){ // Has a mother. We will now check the mother particle's PID.
					// int motherID = pythia.event[motherIdx].idAbs(); // Using the absolute mother index to select the mother.
					int motherID = pythia.event[motherIdx1].id(); // Using the non-absolute index
					// Exclude pions from weak decays
					if (motherID == 310 ||  // K_S^0 -- Am considering its lifetime long enough to be discarded by the DCA cuts.
					motherID == 130 ||  // K_L^0
					motherID == 3122 || // Lambda
					motherID == 3222 || // Sigma+
					// motherID == 3212 || // Sigma0 // Maybe don't include this cut! Sigma0's products will, most of the time, be considered primary particles!
					motherID == 3112 || // Sigma-
					motherID == 3312 || // Xi-
					motherID == 3334 || // Omega-
					motherID == -310 ||  // Anti-K_S^0
					motherID == -130 ||  // Anti-K_L^0
					motherID == -3122 || // Anti-Lambda
					motherID == -3222 || // Anti-Sigma+
					// motherID == -3212 || // Anti-Sigma0 // Maybe don't include this cut! Sigma0's products will, most of the time, be considered primary particles!
					motherID == -3112 || // Anti-Sigma-
					motherID == -3312 || // Anti-Xi+
					motherID == -3334){  // Anti-Omega+
						IsPrimary[i] = false;
						continue;  // Skip this pion and go to the next particle in the event.
						// Always beware of such uses of "continue"! In this case, you would just not save this particle to the total_hist_pion_final event and not
						// increase the n_pion_event variable, which is ok because those variables are legacy!
						// It would be more accurate to use something like a IsPrimary[i] = false; here too, as this knowledge is deeper than the motherIdx1 > 0 && motherIdx2 == 0 case!
							// Wait, is this really necessary with the above if check? That would already include this redundant check here, but whatever...
							// TODO: Check if this could be removed altogether, speeding up the code! A simple check into the already defined IsPrimary would be enough to do the "continue" part.
					}
				}
				// event_hist_pion->Fill((double) pT);
				total_hist_pion_final->Fill((double) pT);
				n_pion_event += 1;
			}
			else if (isfinal && (particle_PID == 2212 or particle_PID == (-1)*2212)){ // Proton
				// event_hist_proton->Fill((double) pT);
				total_hist_proton_final->Fill((double) pT);
				n_proton_event += 1;
			}
			else if (isfinal && (particle_PID == 321 or particle_PID == (-1)*321)){ // Kaon
				// event_hist_kaon->Fill((double) pT);
				total_hist_kaon_final->Fill((double) pT);
				n_kaon_event += 1;
			}
		}
		// Filling the event counters for INEL events for each multiplicity class:
			// (This is actually just to sate my curiosity. It has no practical use, as the number of events is the same for all particle species in the ALICE standard normalization of pT spectra)
		if (charged_in_central_eta){
			hINELEventCounter->Fill(0);
			hINELev_Ntracks_final_charged_forward->Fill(charged_in_back_forward_eta);

			// if (n_charged_event != 0){hINELev_Ntracks_final_charged_forward_Charged->Fill(charged_in_back_forward_eta);}
				// You don't need to check if n_charged_event != 0. That comes imediatelly from the fact that charged_in_central_eta is true.
			hINELev_Ntracks_final_charged_forward_Charged->Fill(charged_in_back_forward_eta);
			if (n_pion_event != 0){hINELev_Ntracks_final_charged_forward_Pion->Fill(charged_in_back_forward_eta);}
			if (n_proton_event != 0){hINELev_Ntracks_final_charged_forward_Proton->Fill(charged_in_back_forward_eta);}
			if (n_kaon_event != 0){hINELev_Ntracks_final_charged_forward_Kaon->Fill(charged_in_back_forward_eta);}
		}

			// First the non-final inclusive histograms:
		hNtracks->Fill(ntrack);
		hNtracks_charged_forward->Fill(Ntracks_charged_forward);
		hNtracks_primary->Fill(Ntracks_primary);
			// Then the histograms that only have isFinal == true particles:
		hNtracks_final->Fill(ntrack_final);
		hNtracks_final_forward->Fill(ntracks_in_back_forward_eta);
		if(charged_in_central_eta){hNtracks_final_forward_INEL_l0->Fill(ntracks_in_back_forward_eta);}
		hNtracks_final_charged_center->Fill(Ntracks_final_charged_center);
		hNtracks_final_charged->Fill(Ntracks_final_charged);
		if(charged_in_central_eta){hNtracks_final_charged_INEL_l0->Fill(Ntracks_final_charged);} // This histogram should only be filled if the event is INEL>0, so it is more restrictive than the histogram above!
		hNtracks_final_charged_forward->Fill(charged_in_back_forward_eta);
		if(charged_in_central_eta){hNtracks_final_charged_forward_INEL_l0->Fill(charged_in_back_forward_eta);} // This histogram should only be filled if the event is INEL>0, so it is more restrictive than the histogram above!

		t3->Fill();

		// pointer_to_hist_for_filling->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        // pointer_to_hist_for_filling->GetYaxis()->SetTitle("d^{2}N/(dp_{T} dy) (GeV/c)^{-1}");
		// std::string title = (std::string) "Event " + std::to_string(iEvent);
		// pointer_to_hist_for_filling->SetTitle(title.c_str());

		// // Now filling the histograms for each particle type:
		if (n_charged_event != 0){hEventCounterCharged->Fill(0);}
		if (n_pion_event != 0){hEventCounterPion->Fill(0);}
		if (n_proton_event != 0){hEventCounterProton->Fill(0);}
		if (n_kaon_event != 0){hEventCounterKaon->Fill(0);}

		// if (event_hist_charged->GetEntries() != 0){
		// 	hEventCounterCharged->Fill(0);

		// 	// event_hist_charged->Scale(1./(event_hist_charged->GetBinWidth(1))); // This is done in a later part of the code
		// 	histogram_copy(event_hist_charged, pointer_to_hist_for_filling, N_bins);
		// }
		// // Will fill the tree even if the histograms are totally empty. This is done in order to properly relate multiplicity to the histogram on later codes, but could be improved...
		// 	// The hEventCounterPion-like histograms will report the actual number of events that had a non-zero count of the desired particle
		// charged_branch->Fill(); // Have to fill this before entering the next loop, where the proton data will be inserted in pointer_to_hist_for_filling

		// if (event_hist_pion->GetEntries() != 0){
		// 	hEventCounterPion->Fill(0);

		// 	// event_hist_pion->Scale(1./(event_hist_pion->GetBinWidth(1))); // This is done in a later part of the code
		// 	histogram_copy(event_hist_pion, pointer_to_hist_for_filling, N_bins);
		// }
		// // Will fill the tree even if the histograms are totally empty. This is done in order to properly relate multiplicity to the histogram on later codes, but could be improved...
		// 	// The hEventCounterPion-like histograms will report the actual number of events that had a non-zero count of the desired particle
		// pion_branch->Fill(); // Have to fill this before entering the next loop, where the proton data will be inserted in pointer_to_hist_for_filling
		
		// if (event_hist_proton->GetEntries() != 0){
		// 	hEventCounterProton->Fill(0);
		// 	// event_hist_proton->Scale(1./(event_hist_proton->GetBinWidth(1)));
		// 	histogram_copy(event_hist_proton, pointer_to_hist_for_filling, N_bins);
		// }
		// proton_branch->Fill();

		// if (event_hist_kaon->GetEntries() != 0){
		// 	hEventCounterKaon->Fill(0);
		// 	// event_hist_kaon->Scale(1./(event_hist_kaon->GetBinWidth(1)));
		// 	histogram_copy(event_hist_kaon, pointer_to_hist_for_filling, N_bins);
		// }
		// kaon_branch->Fill();

		// pointer_to_hist_for_filling->Reset(); // Actually doesn't need this reset, but kept it anyways
		// event_hist_charged->Reset();
		// event_hist_pion->Reset();
		// event_hist_proton->Reset();
		// event_hist_kaon->Reset();
	}

	pythia.stat();

	t3->Print();
	f.cd();

	std::cout << "Writing information into the .root file" << std::endl;
	t3->Write();
	// event_pT_tree->Write();

	hEventCounter->Write();
	hINELEventCounter->Write();
	hEventCounterCharged->Write();
	hEventCounterPion->Write();
	hEventCounterProton->Write();
	hEventCounterKaon->Write();

	total_hist_charged_final->Sumw2();
	total_hist_pion_final->Sumw2();
	total_hist_proton_final->Sumw2();
	total_hist_kaon_final->Sumw2();

	total_hist_charged_final->Write();
	total_hist_pion_final->Write();
	total_hist_proton_final->Write();
	total_hist_kaon_final->Write();

	// Saving the multiplicity histograms that will speed up the centrality selection process!
		// First the non-final inclusive histograms:
	hNtracks->Write();
	hNtracks_charged_forward->Write();
	hNtracks_primary->Write();
		// Then the histograms that only have isFinal == true particles:
	hNtracks_final->Write();
	hNtracks_final_forward->Write();
	hNtracks_final_forward_INEL_l0->Write();
	hNtracks_final_charged_center->Write();
	hNtracks_final_charged->Write();
	hNtracks_final_charged_INEL_l0->Write();
	hNtracks_final_charged_forward->Write();
	hNtracks_final_charged_forward_INEL_l0->Write();

	// Processing these multiplicity histograms into centrality histograms:
	std::cout << "\nProcessing centrality" << std::endl;
		// First, declaring the new histograms using the information of maximum and minimum multiplicity of the TH1I's:
		// The number of bins is exactly the maximum number of multiplicity + 1 (the upper limit is non-inclusive!).
		// Index 0 will be multiplicity [0, 1), index 1 will be [1, 2) and so on until the last index.
		// If the maximum value (the last one) was 10, then we need 11 bins to get the [10, 11) bin as a regular bin, not as the overflow bin!
	// TH1D *hNtracks_to_centrality = new TH1D("hNtracks_to_centrality", "hNtracks_to_centrality", hNtracks->GetMaximum() + 1, 0, hNtracks->GetMaximum() + 1);
		// Updated to use FindLastBinAbove() instead of maximum value! I don't want the highest value in the TH1I, I want the highest bin with non-zero value!
	TH1D *hNtracks_to_centrality = new TH1D("hNtracks_to_centrality", "hNtracks_to_centrality", hNtracks->FindLastBinAbove(0) + 1, 0, hNtracks->FindLastBinAbove(0) + 1);
	TH1D *hNtracks_to_centrality_charged_forward = new TH1D("hNtracks_to_centrality_charged_forward", "hNtracks_to_centrality_charged_forward", hNtracks_charged_forward->FindLastBinAbove(0) + 1, 0, hNtracks_charged_forward->FindLastBinAbove(0) + 1);
	TH1D *hNtracks_to_centrality_primary = new TH1D("hNtracks_to_centrality_primary", "hNtracks_to_centrality_primary", hNtracks_primary->FindLastBinAbove(0) + 1, 0, hNtracks_primary->FindLastBinAbove(0) + 1);
	TH1D *hNtracks_to_centrality_final = new TH1D("hNtracks_to_centrality_final", "hNtracks_to_centrality_final", hNtracks_final->FindLastBinAbove(0) + 1, 0, hNtracks_final->FindLastBinAbove(0) + 1);
	TH1D *hNtracks_to_centrality_final_forward = new TH1D("hNtracks_to_centrality_final_forward", "hNtracks_to_centrality_final_forward", hNtracks_final_forward->FindLastBinAbove(0) + 1, 0, hNtracks_final_forward->FindLastBinAbove(0) + 1);
	TH1D *hNtracks_to_centrality_final_forward_INEL_l0 = new TH1D("hNtracks_to_centrality_final_forward_INEL_l0", "hNtracks_to_centrality_final_forward_INEL_l0", hNtracks_final_forward_INEL_l0->FindLastBinAbove(0) + 1, 0, hNtracks_final_forward_INEL_l0->FindLastBinAbove(0) + 1);
	TH1D *hNtracks_to_centrality_final_charged_center = new TH1D("hNtracks_to_centrality_final_charged_center", "hNtracks_to_centrality_final_charged_center", hNtracks_final_charged_center->FindLastBinAbove(0) + 1, 0, hNtracks_final_charged_center->FindLastBinAbove(0) + 1);
	TH1D *hNtracks_to_centrality_final_charged = new TH1D("hNtracks_to_centrality_final_charged", "hNtracks_to_centrality_final_charged", hNtracks_final_charged->FindLastBinAbove(0) + 1, 0, hNtracks_final_charged->FindLastBinAbove(0) + 1);
	TH1D *hNtracks_to_centrality_final_charged_INEL_l0 = new TH1D("hNtracks_to_centrality_final_charged_INEL_l0", "hNtracks_to_centrality_final_charged_INEL_l0", hNtracks_final_charged_INEL_l0->FindLastBinAbove(0) + 1, 0, hNtracks_final_charged_INEL_l0->FindLastBinAbove(0) + 1);
	TH1D *hNtracks_to_centrality_final_charged_forward = new TH1D("hNtracks_to_centrality_final_charged_forward", "hNtracks_to_centrality_final_charged_forward", hNtracks_final_charged_forward->FindLastBinAbove(0) + 1, 0, hNtracks_final_charged_forward->FindLastBinAbove(0) + 1);
	TH1D *hNtracks_to_centrality_final_charged_forward_INEL_l0 = new TH1D("hNtracks_to_centrality_final_charged_forward_INEL_l0", "hNtracks_to_centrality_final_charged_forward_INEL_l0", hNtracks_final_charged_forward_INEL_l0->FindLastBinAbove(0) + 1, 0, hNtracks_final_charged_forward_INEL_l0->FindLastBinAbove(0) + 1);

	multiplicity_to_centrality(hNtracks, hNtracks_to_centrality);
	multiplicity_to_centrality(hNtracks_charged_forward, hNtracks_to_centrality_charged_forward);
	multiplicity_to_centrality(hNtracks_primary, hNtracks_to_centrality_primary);
	multiplicity_to_centrality(hNtracks_final, hNtracks_to_centrality_final);
	multiplicity_to_centrality(hNtracks_final_forward, hNtracks_to_centrality_final_forward);
	multiplicity_to_centrality(hNtracks_final_forward_INEL_l0, hNtracks_to_centrality_final_forward_INEL_l0);
	multiplicity_to_centrality(hNtracks_final_charged_center, hNtracks_to_centrality_final_charged_center);
	multiplicity_to_centrality(hNtracks_final_charged, hNtracks_to_centrality_final_charged);
	multiplicity_to_centrality(hNtracks_final_charged_INEL_l0, hNtracks_to_centrality_final_charged_INEL_l0);
	multiplicity_to_centrality(hNtracks_final_charged_forward, hNtracks_to_centrality_final_charged_forward);
	multiplicity_to_centrality(hNtracks_final_charged_forward_INEL_l0, hNtracks_to_centrality_final_charged_forward_INEL_l0);

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

	hINELev_Ntracks_final_charged_forward->Write();
	hINELev_Ntracks_final_charged_forward_Charged->Write();
	hINELev_Ntracks_final_charged_forward_Pion->Write();
	hINELev_Ntracks_final_charged_forward_Proton->Write();
	hINELev_Ntracks_final_charged_forward_Kaon->Write();

	return 0;
}


// void histogram_copy(TH1D *origin_hist, TH1D *destiny_hist, int N_bins){
//     for(int i = 1; i <= N_bins; i++){ // Do notice the first bin, 0, is the underflow bin, and the n_bins bin is the last before overflow!
//         double bin_value = origin_hist->GetBinContent(i);
//         double bin_error = origin_hist->GetBinError(i);

// 		// if (bin_value != 0){ // Did this because the SetBinContent would make the number of entries equal to the number of bins it had Set the contents of.
// 		// 	for (int j = 0; j < bin_value; j++){
// 		// 		destiny_hist->AddBinContent(i); // With hope, this will avoid the previous error in the number of entries of each histogram
// 		// 		// Could probably use the TH1::AddBinContent (Int_t bin, Double_t w) version to avoid a loop, but I remind it having some problems on the ROOT version of jarvis1
// 		// 	}
// 		// }
// 		// On this weird alternative version, had to reset the histogram with pointer_to_hist_for_filling->Reset(); after each fill!

// 		destiny_hist->SetBinContent(i, bin_value);
//         destiny_hist->SetBinError(i, bin_error);
// 		destiny_hist->ResetStats(); // Hopefully this will solve the STILL WRONG counts on the final histogram... YES! This fixes it!
//     }
// }

// Converting the event multiplicity to centrality in the most precise way possible -- Corresponding each integer number of particles to a limit in centrality.
void multiplicity_to_centrality(TH1I *multiplicity_hist, TH1D *centrality_hist){
	int Nev = multiplicity_hist->GetEntries();
	int Nbins = centrality_hist->GetNbinsX(); // I want to loop only on the [0 particles, maximum_number_of_particles] region, not in the [0, 10.000] particles region that the hNtracks histograms need to have to avoid missing any data.

	int number_of_events_up_to_current_bin = 0; // The number of events that came up to this bin, including this bin.
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
