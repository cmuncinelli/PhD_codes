// Author: Gianni S. S. Liveraro, adapted by Cicero D. Muncinelli
#include "Pythia8/Pythia.h" // Not just Pythia.h. You have to take an extra step after the $PYTHIA8 folder.
#include "Pythia8/HeavyIons.h" // Needed for Angantyr
#include "Pythia8Plugins/FastJet3.h"

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"
// #include "TH1I.h" // Not using this class anymore!
#include "TH2.h"
#include "TH2D.h"
#include "TMath.h"

#include <dirent.h>
#include <random> // For the seed generation for each worker
#include <filesystem> // Useful for manipulation of the input_card_path
#include <omp.h> // My own way of manually parallelizing Pythia -- More flexible than PythiaParallel's per-event call to my RunWorker function!
#include <TROOT.h> // Necessary for the EnableThreadSafety() part of ROOT

#include <iostream>
#include <map> // For the PID tracking of the jet constituents

const double PI = TMath::Pi();

// void histogram_copy(TH1D *origin_hist, TH1D *destiny_hist, int N_bins);
void multiplicity_to_centrality(TH1D *multiplicity_hist, TH1D *centrality_hist);
inline double wrapToInterval(double phi, double phi_min, double phi_max);
void FillHistogramFromMap(TH1* hist, const std::map<int,int>& map);

using namespace Pythia8;

void RunWorker(int WorkerId, int N_ev, const std::string output_folder, const std::string input_card_path, const std::string input_card_name){
  	Pythia pythia;
	
	// // Some settings related to output in init(), next() and stat() that came from the Parnassus settings
	// pythia.readString("Init:showChangedSettings = on"); // list changed settings
	// pythia.readString("Init:showChangedParticleData = off"); // list changed particle data
	// pythia.readString("Next:numberCount = 0"); // print message every n events
	// pythia.readString("Next:numberShowInfo = 0"); // print event information n times
	// pythia.readString("Next:numberShowProcess = 1"); // print process record n times
	// pythia.readString("Next:numberShowEvent = 1"); // print event record n times
	
		// Actually, making Pythia quiet (could possibly re-route this into a .log file, but for testing this is OK):
	pythia.readString("Print:quiet = on"); // suppress banner and setup info
	pythia.readString("Init:showProcesses = off"); // don't list processes
	pythia.readString("Init:showChangedSettings = off");
	pythia.readString("Init:showMultipartonInteractions = off");
	pythia.readString("Next:numberCount = 0"); // don't print event counter
	pythia.readString("Next:numberShowLHA = 0");
	pythia.readString("Next:numberShowInfo = 0");
	pythia.readString("Next:numberShowProcess = 0");
	pythia.readString("Next:numberShowEvent = 0");

	// Event settings:
	pythia.readFile(input_card_path);

	// Seeding each worker core differently, based on a clock+WorkerId mixing:
		// If you just give seed = 0, it might get the same seed for all workers, based only on the starting clock for all of them!
	std::random_device rd;
	std::mt19937 rng(rd() + static_cast<unsigned int>(WorkerId));
    std::uniform_int_distribution<int> dist(1, 900000000); // Pythia expects values between 1 and 900.000.000 for its seed. Zero would be clock seeding, and that can be troublesome for multi-threading
    int seed = dist(rng);
	pythia.readString("Random:seed = " + std::to_string(seed));

	// pythia.readString("Random:seed = 42"); // Kept for testing!

	pythia.init();

	// Setting up FastJet finder with anti-kT method:
	// (for more information, see logs 503 and 504)
	double R = 0.4;
	double jet_min_pT = 8.0; // Changed from 2.0 to 8.0: we were having 10% of events having 20 or more jets! The 8.0 GeV/c comes from Youpeng's PAG presentation
	double ALICE_charged_jet_acceptance = 0.9; // Look at the Detector_Acceptances_and_FastJet.pdf document: TPS, ITS, TRD and TOF cover <0.9, and EMCal does not cover full azimuth!
	double jet_max_eta = ALICE_charged_jet_acceptance - R; // The jet must be entirely in the region where ALICE can actually see it!
														   // Notice that you don't need to do an extra |y|<0.5 cut if the jet pT is sufficiently high, because y ~ \eta in that limit.
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, R);
  	std::vector <fastjet::PseudoJet> FastJetInputs;

		// Creating the file to associate all current histograms and trees to the current gDirectory
		// (i.e., the current folder) instead of the global ROOT gDirectory.
	std::string filename = (std::string) output_folder;
	filename += "/output_core" + to_string(WorkerId) + ".root"; // A simple naming scheme. Left the whole event naming and info for the output_folder name.

	TFile f(filename.c_str(), "RECREATE");
	TTree *t3 = new TTree("t3","Reconst ntuple");

	// const Int_t kMaxTrack = 7000; // This could be a problem for very large multiplicity pp events and PbPb events!
	// (Memory corruption and the sort, if the vectors have more than 7000 particles, which may be final or intermediary states too!)
	// Thus, changed to a C++ approach instead of the usual "Int_t ID[kMaxTrack];" C-style arrays for the branches.
	Int_t ntrack;
	Int_t ntrack_final;
	Int_t charged_in_back_forward_eta; // The multiplicity indicator for 1807.11321 [nucl-ex]'s data
	Int_t ntracks_in_back_forward_eta; // Another multiplicity indicator, trying to mimic the indirect charged particle production due to neutral's interactions with V0. Something finer than just Nch, and will also receive only final particles
	Int_t Ntracks_charged_forward;
	Int_t Ntracks_primary;
	Int_t Ntracks_final_charged_center;
	Int_t Ntracks_final_charged;
	
	std::vector<Float_t> ProtonOrPBar_pt;
	std::vector<Float_t> ProtonOrPBar_y;
	std::vector<Float_t> ProtonOrPBar_Phi;
	std::vector<Float_t> ProtonOrPBar_m;

	std::vector<Float_t> LambdaOrLBar_pt;
	std::vector<Float_t> LambdaOrLBar_y;
	std::vector<Float_t> LambdaOrLBar_Phi;
	std::vector<Float_t> LambdaOrLBar_m;

	std::vector<Bool_t> isProton_not_PBar;
	std::vector<Bool_t> isLambda_notLBar;

	std::vector<Bool_t> isLambdaOrLBarExperimentalPrimary;

	// Don't need to store any of these for the Lambda polarization analysis! They were just occupying space in disk!
	// std::vector<Int_t> ID;
	// std::vector<Int_t> mother1_ID;
	// std::vector<Int_t> mother2_ID;
	// std::vector<Bool_t> IsFinal;
	// std::vector<Bool_t> IsPrimary;
	// std::vector<Bool_t> IsCarbonCopy;
	// // std::vector<Float_t> px; // Already saving pT,y,phi,mass
	// // std::vector<Float_t> py;
	// // std::vector<Float_t> pz;
	// std::vector<Float_t> pt;
	// std::vector<Float_t> m;
	// // std::vector<Float_t> e; // No need for charge

	// // std::vector<Bool_t> charged; // Don't actually need to save this! This is purely a PDG-lookup information, so no need to store it when the
	// 								// PID is already enough to know it! The isCharged() will be useful inside the loop, but no need to store it!
	Bool_t charged_in_central_eta; // A single boolean that is reset at the start of each new event and marks that event as INEL>0 (or not).

	// // std::vector<Float_t> Eta;
	// std::vector<Float_t> y;
	// std::vector<Float_t> Phi;

	//////////////////////////////////////////////////////////////////////
	// Including FastJet-estimated jet kinematic variables:
	// (CAUTION! Do notice that I expect to save only one jet per event!)
	Float_t pt_jet;
	Float_t m_jet;
	Float_t y_jet;
	Float_t Phi_jet;
	//////////////////////////////////////////////////////////////////////

	TH1D *hEventCounter = new TH1D ("hEventCounter", "Total number of simulated events", 1, -1, 1);
	TH1D *hEventCounterWithLambdaOrBar = new TH1D ("hEventCounterWithLambdaOrBar", "hEventCounterWithLambdaOrBar", 1, -1, 1);
	TH1D *hEventCounterWithJets_pTleqJetMinPt = new TH1D ("hEventCounterWithJets_pTleqJetMinPt", "hEventCounter(jets w/ pT>8)", 1, -1, 1); // Jets in the interval pT > jet_min_pT, all Eta
	TH1D *hEventCounterWithJets_pTleqJetMinPt_EtaCuts = new TH1D ("hEventCounterWithJets_pTleqJetMinPt_EtaCuts", "hEventCounter(jets w/ pT>8,|eta|<0.5)", 1, -1, 1); // Jets in the interval pT > jet_min_pT and |eta| < 0.9 - R
	TH1D *hEventCounterWithLambdaOrBarWithJets_pTleqJetMinPt_EtaCuts = new TH1D ("hEventCounterWithLambdaOrBarWithJets_pTleqJetMinPt_EtaCuts", "hEventCounter(jets w/ pT>8,|eta|<0.5),hasLambdaorLambdaBar", 1, -1, 1);
	TH1D *hEventCounterPassedYExtraCheck = new TH1D ("hEventCounterPassedYExtraCheck", "hEventCounter(jets w/ pT>8,|eta|<0.5,|y|<0.5),hasLambda in pTandYcuts", 1, -1, 1); // Contains Lambda in the desired pT and |y| range, contains a jet in the desired range
	TH1D *hINELEventCounter = new TH1D ("hINELEventCounter", "", 1, -1, 1);
	TH1D *hEventCounterCharged = new TH1D ("hEventCounterCharged", "", 1, -1, 1);
	TH1D *hEventCounterPion = new TH1D ("hEventCounterPion", "", 1, -1, 1);
	TH1D *hEventCounterProton = new TH1D ("hEventCounterProton", "", 1, -1, 1);
	TH1D *hEventCounterKaon = new TH1D ("hEventCounterKaon", "", 1, -1, 1);
	
	TH1D *hJetCounterPerEvent = new TH1D ("hJetCounterPerEvent", "Number of jets found in each event", 40, 0, 40);
	TH1D *hDiscardedHighPtJetCounterPerEvent = new TH1D ("hDiscardedHighPtJetCounterPerEvent", "Number of jets discarded until |eta|<jet_max_eta jet was found", 40, 0, 40);
	TH1D *hNParticlesPerJet = new TH1D ("hNParticlesPerJet", "N_particles in each jet(pT>8)", 50, 0, 50);
	TH1D *hNParticlesLeadingJet_EtaCuts = new TH1D ("hNParticlesLeadingJet_EtaCuts", "N_particles in each |eta|<0.5 jet(pT>8)", 50, 0, 50);
	TH1D *hNParticlesLeadingJet_EtaCutsWithLambda = new TH1D ("hNParticlesLeadingJet_EtaCutsWithLambda", "N_particles in each |eta|<0.5,w/Lambda,jet(pT>8)", 50, 0, 50);
	TH1D *hNParticlesLeadingJet_PassedYExtraCheck = new TH1D ("hNParticlesLeadingJet_PassedYExtraCheck", "N_particles in |y|<0.5,w/Lambda,jet(pT>8)", 50, 0, 50);
		// Some maps to know the PID of the particles on each type of jet:
	std::map<int,int> mapCountsOfPIDAllJetConstituents;
	std::map<int,int> mapCountsOfPIDLeadingJetConstituents_EtaCuts;
	std::map<int,int> mapCountsOfPIDLeadingJetConstituents_EtaCutsWithLambda;
	std::map<int,int> mapCountsOfPIDLeadingJetConstituents_PassedYExtraCheck;
	
	// The histograms themselves will be declared AFTER we know what PIDs are actually contained in the events
	// hPIDOfAllJetConstituents
	// hPIDLeadingJetConstituents_EtaCuts
	// hPIDLeadingJetConstituents_EtaCutsWithLambda
	// hPIDLeadingJetConstituents_PassedYExtraCheck


	int N_bins_pT = 800;
	double upper_limit_pT = 50;

	t3->Branch("ntrack",&ntrack,"ntrack/I");
	t3->Branch("ntrack_final",&ntrack_final,"ntrack_final/I");

	t3->Branch("charged_in_back_forward_eta",&charged_in_back_forward_eta,"charged_in_back_forward_eta/I");
	t3->Branch("ntracks_in_back_forward_eta",&ntracks_in_back_forward_eta,"ntracks_in_back_forward_eta/I");
	t3->Branch("charged_in_central_eta",&charged_in_central_eta,"charged_in_central_eta/O");
	t3->Branch("Ntracks_charged_forward",&Ntracks_charged_forward,"Ntracks_charged_forward/I");
	t3->Branch("Ntracks_primary",&Ntracks_primary,"Ntracks_primary/I");
	t3->Branch("Ntracks_final_charged_center",&Ntracks_final_charged_center,"Ntracks_final_charged_center/I");
	t3->Branch("Ntracks_final_charged",&Ntracks_final_charged,"Ntracks_final_charged/I");

		// Including FastJet-estimated jet kinematic variables:
	t3->Branch("pt_jet", &pt_jet, "pt_jet/F");
	t3->Branch("m_jet", &m_jet, "m_jet/F");
	t3->Branch("y_jet", &y_jet, "y_jet/F");
	t3->Branch("Phi_jet", &Phi_jet, "Phi_jet/F");
	
	// Only the bare minimum information is stored now:
	// -- Lambda, LambdaBar, p and pBar information, and only if the P or PBar came from a Lambda or LambdaBar!
	// (already using the new C++ style of branch declaration, as described below in the LEGACY CODE section)
	t3->Branch("ProtonOrPBar_pt", &ProtonOrPBar_pt);
	t3->Branch("ProtonOrPBar_y", &ProtonOrPBar_y);
	t3->Branch("ProtonOrPBar_Phi", &ProtonOrPBar_Phi);
	t3->Branch("ProtonOrPBar_m", &ProtonOrPBar_m);

	t3->Branch("LambdaOrLBar_pt", &LambdaOrLBar_pt);
	t3->Branch("LambdaOrLBar_y", &LambdaOrLBar_y);
	t3->Branch("LambdaOrLBar_Phi", &LambdaOrLBar_Phi);
	t3->Branch("LambdaOrLBar_m", &LambdaOrLBar_m);

	t3->Branch("isProton_not_PBar", &isProton_not_PBar);
	t3->Branch("isLambda_notLBar", &isLambda_notLBar);

	t3->Branch("isLambdaOrLBarExperimentalPrimary", &isLambdaOrLBarExperimentalPrimary);

	// LEGACY CODE! Don't need to store all of this information for the processing of Lambda polarization, so I can reduce the size of the datasets!
	// // t3->Branch("IsFinal",IsFinal,"IsFinal[ntrack]/O"); // Actually don't need to pass this array by reference, as it is a C++ array and it would be passed by ref by default!
	// 	// In the new declaration style, with C++ arrays, you don't need to declare types nor lengths! ROOT manages all that bookkeeping for you!
	// t3->Branch("IsFinal",&IsFinal);
	// t3->Branch("IsPrimary",&IsPrimary); 
	// // Do notice that the [ntrack] part is making it so that only the values from 0 up to ntrack - 1 will be stored in your TTree, so you don't need to reset IsPrimary for each new
	// // event in the loop!
	// t3->Branch("IsCarbonCopy",&IsCarbonCopy); 

	// t3->Branch("ID",&ID); // A PID identifier
	// t3->Branch("Mother1_ID",&mother1_ID);
	// t3->Branch("Mother2_ID",&mother2_ID);
	// // t3->Branch("px",&px);
	// // t3->Branch("py",&py);
	// // t3->Branch("pz",&pz);
	// t3->Branch("pt",&pt);
	// t3->Branch("m",&m);
	// // t3->Branch("e",&e);
	// // t3->Branch("charged",&charged);

	// // t3->Branch("Eta",&Eta);
	// t3->Branch("y",&y);
	// t3->Branch("Phi",&Phi);

	// Saving a global TH1D for each particle, just to see the information on their pT:
	TH1D *pT_hist_charged_final = new TH1D("pT_hist_charged_final", "pT_hist_charged_final", N_bins_pT, 0, upper_limit_pT);
	TH1D *pT_hist_pion_final = new TH1D("pT_hist_pion_final", "pT_hist_pion_final", N_bins_pT, 0, upper_limit_pT);
	TH1D *pT_hist_proton_final = new TH1D("pT_hist_proton_final", "pT_hist_proton_final", N_bins_pT, 0, upper_limit_pT);
	TH1D *pT_hist_kaon_final = new TH1D("pT_hist_kaon_final", "pT_hist_kaon_final", N_bins_pT, 0, upper_limit_pT);

	// Saving multiplicity information in a TH1I (which consumes at most 40 kB with this many bins!) -- Faster centrality/multiplicity class fetching:
	// Had to change to TH1D because was experiencing some overflows and rounding errors when trying to add a count when the bin already had many entries!
	int multiplicity_nbins = 50000; // The number of bins is also the number of entries, which gives one bin per possible multiplicity!
		// Could make this stop at 7000 bins, as the kMaxTrack is 7000 for my pythia simulations, but this format is more general and the disk size is ridiculously small...
		// By declaring in this way, GetBin(0) will get the bin of events with 0 particles in multiplicity. The value 10000 is not stored, and it will go to the overflow bin.
		// In other words, GetBin(N) will get you the exact bin with N particles! A convenient indexing, where [0] is the index of 0 particles and [1] for 1 particle!
	TH1D *hNtracks = new TH1D("hNtracks", "hNtracks", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hNtracks_charged_forward = new TH1D("hNtracks_charged_forward", "hNtracks_charged_forward", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hNtracks_primary = new TH1D("hNtracks_primary", "hNtracks_primary", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hNtracks_final = new TH1D("hNtracks_final", "hNtracks_final", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hNtracks_final_forward = new TH1D("hNtracks_final_forward", "hNtracks_final_forward", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hNtracks_final_forward_INEL_l0 = new TH1D("hNtracks_final_forward_INEL_l0", "hNtracks_final_forward_INEL_l0", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hNtracks_final_charged_center = new TH1D("hNtracks_final_charged_center", "hNtracks_final_charged_center", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hNtracks_final_charged = new TH1D("hNtracks_final_charged", "hNtracks_final_charged", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hNtracks_final_charged_INEL_l0 = new TH1D("hNtracks_final_charged_INEL_l0", "hNtracks_final_charged_INEL_l0", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hNtracks_final_charged_forward = new TH1D("hNtracks_final_charged_forward", "hNtracks_final_charged_forward", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hNtracks_final_charged_forward_INEL_l0 = new TH1D("hNtracks_final_charged_forward_INEL_l0", "hNtracks_final_charged_forward_INEL_l0", multiplicity_nbins, 0, multiplicity_nbins);
	// Caution! All the non INEL>0 (INEL_l0 or "larger than 0" in the above notation) may include multiplicity = 0 events, which will be assigned a position in the centrality determination
	// as the most peripheral events!
	// This makes those estimators not so trustable, as we wouldn't usually consider 0-particle events in the centrality estimator in V0 (Or so I think. It might be good to check this!)

		// Now storing information to recover the number of events in each centrality class:
		// You don't actually need this information for the root_to_csv_custom_PYTHIA.cxx to work, and it won't make it much faster,
		// but it is nice to be able to visualize it!
	TH1D *hINELev_Ntracks_final_charged_forward = new TH1D ("hINELev_Ntracks_final_charged_forward", "", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hINELev_Ntracks_final_charged_forward_Charged = new TH1D ("hINELev_Ntracks_final_charged_forward_Charged", "", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hINELev_Ntracks_final_charged_forward_Pion = new TH1D ("hINELev_Ntracks_final_charged_forward_Pion", "", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hINELev_Ntracks_final_charged_forward_Proton = new TH1D ("hINELev_Ntracks_final_charged_forward_Proton", "", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hINELev_Ntracks_final_charged_forward_Kaon = new TH1D ("hINELev_Ntracks_final_charged_forward_Kaon", "", multiplicity_nbins, 0, multiplicity_nbins);

		// Some values for the Jet kinematic analysis histograms:
	int N_bins_phi = 100;
    double dphi = 2*PI * 1./N_bins_phi; // Bin width
    double phi_min = -dphi/2.; // The lower edge is defined in such a way that phi = 0 is the bin center.
    double phi_max = (2*PI - dphi) + dphi/2.;
    double phi_min_Ring = phi_min - PI;
    double phi_max_Ring = phi_max - PI;

		// All jets with minimum pT threshold:
		// (This might include jets with a larger pT than the one jet we selected as the leading jet!
		// This will help to identify if we are discarding a large number of jets before we get to the
		// ideal leading jet! In other words, this is a counter that includes the selected leading jet,
		// but also includes all other jets with pT > pT_selected_leading_jet that did not pass the
		// maximum jet_max_eta selection.)
	auto hLeadJetY_DiscardedUntilLead_pTleqJetMinPt = new TH1D("hLeadJetY_DiscardedUntilLead_pTleqJetMinPt", "Jets w/ min_pT threshold, and pT>=selected_lead_jet", 100, -2.5, 2.5);
    auto hLeadJetPt_DiscardedUntilLead_pTleqJetMinPt = new TH1D("hLeadJetPt_DiscardedUntilLead_pTleqJetMinPt", "Jets w/ min_pT threshold, and pT>=selected_lead_jet", 100, 0, 50);
    auto hLeadJetPtY_DiscardedUntilLead_pTleqJetMinPt = new TH2D("hLeadJetPtY_DiscardedUntilLead_pTleqJetMinPt", "Jets w/ min_pT threshold, and pT>=selected_lead_jet", 40, 0, 50, 40, -2.5, 2.5);
    auto hLeadJetPhi_DiscardedUntilLead_pTleqJetMinPt = new TH1D("hLeadJetPhi_DiscardedUntilLead_pTleqJetMinPt", "Jets w/ min_pT threshold, and pT>=selected_lead_jet", 100, -PI, PI); // Pythia gives us Phi from -PI to PI! (thus also got the phi_std from FastJet)

    auto hLeadJetY_pTleqJetMinPt_eta_cut = new TH1D("hLeadJetY_pTleqJetMinPt_eta_cut", "(minpT and maxEta OK)", 100, -2.5, 2.5);
    auto hLeadJetPt_pTleqJetMinPt_eta_cut = new TH1D("hLeadJetPt_pTleqJetMinPt_eta_cut", "(minpT and maxEta OK)", 100, 0, 50);
    auto hLeadJetPtY_pTleqJetMinPt_eta_cut = new TH2D("hLeadJetPtY_pTleqJetMinPt_eta_cut", "(minpT and maxEta OK)", 40, 0, 50, 40, -2.5, 2.5);
    auto hLeadJetPhi_pTleqJetMinPt_eta_cut = new TH1D("hLeadJetPhi_pTleqJetMinPt_eta_cut", "(minpT and maxEta OK)", 100, -PI, PI);

	auto hfirstDiscardedJetY_minus_leadingJet_MinPtCut = new TH1D("hfirstDiscardedJetY_minus_leadingJet_MinPtCut", "(minpT_OK,yFirstDiscardedJet-ySelectedJet)", 100, -2.5, 2.5);
	auto hfirstDiscardedJetPt_MinPtCut = new TH1D("hfirstDiscardedJetPt_MinPtCut", "(minpT_OK,FirstDiscardedJet)", 100, 0, 50);
	auto hfirstDiscardedJetPtY_minus_leadingJet_MinPtCut = new TH2D("hfirstDiscardedJetPtY_minus_leadingJet_MinPtCut", "(minpT_OK,yFirstDiscardedJet-ySelectedJet)", 40, 0, 50, 40, -2.5, 2.5);
	auto hfirstDiscardedJetPhi_minus_leadingJet_MinPtCut = new TH1D("hfirstDiscardedJetPhi_minus_leadingJet_MinPtCut", "(minpT_OK,PhiFirstDiscardedJet-PhiSelectedJet)", 100, -PI, PI);

	auto hJetProxyY_pTleqJetMinPt_eta_cut_WithLambdaOrBar = new TH1D("hJetProxyY_pTleqJetMinPt_eta_cut_WithLambdaOrBar", "(minpT, maxEta OK, contains LambdaOrLBar)", 100, -2.5, 2.5);
    auto hJetProxyPt_pTleqJetMinPt_eta_cut_WithLambdaOrBar = new TH1D("hJetProxyPt_pTleqJetMinPt_eta_cut_WithLambdaOrBar", "(minpT, maxEta OK, contains LambdaOrLBar)", 100, 0, 50);
    auto hJetProxyPtY_pTleqJetMinPt_eta_cut_WithLambdaOrBar = new TH2D("hJetProxyPtY_pTleqJetMinPt_eta_cut_WithLambdaOrBar", "(minpT, maxEta OK, contains LambdaOrLBar)", 40, 0, 50, 40, -2.5, 2.5);
    auto hJetProxyPhi_pTleqJetMinPt_eta_cut_WithLambdaOrBar = new TH1D("hJetProxyPhi_pTleqJetMinPt_eta_cut_WithLambdaOrBar", "(minpT, maxEta OK, contains LambdaOrLBar)", 100, -PI, PI);

    auto hJetProxyY_YExtraCheck = new TH1D("hJetProxyY_YExtraCheck", "(minpT and maxY OK, contains LambdaOrLBar)", 100, -2.5, 2.5);
    auto hJetProxyPt_YExtraCheck = new TH1D("hJetProxyPt_YExtraCheck", "(minpT and maxY OK, contains LambdaOrLBar)", 100, 0, 50);
    auto hJetProxyPtY_YExtraCheck = new TH2D("hJetProxyPtY_YExtraCheck", "(minpT and maxY OK, contains LambdaOrLBar)", 40, 0, 50, 40, -2.5, 2.5);
    auto hJetProxyPhi_YExtraCheck = new TH1D("hJetProxyPhi_YExtraCheck", "(minpT and maxY OK, contains LambdaOrLBar)", 100, -PI, PI);

	// Event loop
	for (int iEvent=0; iEvent<N_ev; ++iEvent){
		// pythia.next();
		if (!pythia.next()){ // Skip failed events --> Angantyr can generate events that have no particles at all!
			iEvent--; // Before skipping, will reset the counter increase of this non-generated event: we need EXACTLY N_ev to be generated
			continue;
		}

		hEventCounter->Fill(0);
		int five_percent_step = 0.05 * N_ev; // Defined this outside the check to avoid some divisions by zero that could happen when using low N_ev per worker.
		if ((WorkerId < 10) && five_percent_step > 0 ? (iEvent % five_percent_step == 0) : false){ // Printing for just 10 workers is already more than enough!
			std::cout << "[Worker " << WorkerId << "] Now on event " << iEvent << " of " << N_ev << " (" << (100.0 * iEvent / N_ev) << "%)" << std::endl;
		}

		ntrack = pythia.event.size(); // This could be a bit problematic to estimate multiplicity: this includes non-final particles!
		ntrack_final = 0; // Resetting for the next event

		charged_in_central_eta = false; // Resetting the bool for the next event
		charged_in_back_forward_eta = 0; // Actually ""charged_in_back_forward_eta_FINAL" cause it only gets final particles
		ntracks_in_back_forward_eta = 0;

		// Some variables that will not be stored, but will be used inside the loop for counting particles -- These serve only for studying and are not necessarily useful for any analysis:
		Ntracks_charged_forward = 0;
		Ntracks_primary = 0;
		Ntracks_final_charged = 0;
		Ntracks_final_charged_center = 0;

		int n_charged_event = 0;
		int n_pion_event = 0;
		int n_proton_event = 0;
		int n_kaon_event = 0;

		// Information for the jet loop and needed initializations for FastJet:
		pt_jet = 0; // Reinitializing these values as well, even though this shouldn't be necessary
		m_jet = 0;
		y_jet = 0;
		Phi_jet = 0;
		Bool_t contains_lambda_or_lambdabar = false; // In this new code, will exclude events that have neither a Lambda nor a Lambda bar

		// Track loop:
			// Clearing all C++ arrays before the loop starts, so that you safely store only the information about the current event:
		ProtonOrPBar_pt.clear();
		ProtonOrPBar_y.clear();
		ProtonOrPBar_Phi.clear();
		ProtonOrPBar_m.clear();
		LambdaOrLBar_pt.clear();
		LambdaOrLBar_y.clear();
		LambdaOrLBar_Phi.clear();
		LambdaOrLBar_m.clear();
		isProton_not_PBar.clear();
		isLambda_notLBar.clear();
		isLambdaOrLBarExperimentalPrimary.clear();

			// No longer useful arrays:
		// ID.clear();
		// mother1_ID.clear();
		// mother2_ID.clear();
		// IsFinal.clear();
		// IsPrimary.clear();
		// IsCarbonCopy.clear();
		// // px.clear();
		// // py.clear();
		// // pz.clear();
		// pt.clear();
		// m.clear();
		// // e.clear();
		// // charged.clear();
		// // Eta.clear();
		// y.clear();
		// Phi.clear();

			// Clearing the list of FastJet candidates for this loop:
		FastJetInputs.clear();
		for (int i = 0; i < pythia.event.size(); ++i){
			bool isfinal = pythia.event[i].isFinal();
			bool isvisible = pythia.event[i].isVisible();
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
			
			Bool_t IsCarbonCopy_value = false; // This should be set as false by default for each particle! 
											   // If you don't declare this for each particle, the vector 
											   // won't be updated correctly and eventually will become full of trues!
			if ((daughter1 == daughter2) && (daughter1 > 0)){
				IsCarbonCopy_value = true; // This will be a useful bool to consider, where you will skip the particle if it is simply a carbon-copy: it may not be a final state,
				// but a carbon-copy isn't even a physical particle, and it shouldn't be considered if you want to count the number of intermediate particles in your collision.
				// It is best not to do an if ((daughter1 == daughter2) && (daughter1 > 0)){continue} if you are going to save these particles in a TTree, because the value on
				// the i'th index of the vector that will be flushed into the TTree will still keep the previous iteration's particle data!
			}
			// IsCarbonCopy.push_back(IsCarbonCopy_value); // Moved all push_backs to the end of the loop

			// if (isfinal){ntrack_final += 1;}

			Float_t pT = pythia.event[i].pT();
			Float_t rapidity = pythia.event[i].y();
			Float_t Phi = pythia.event[i].phi();
			Float_t mass = pythia.event[i].m();

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

			int particle_PID = pythia.event[i].id();

				// Updating the "has Lambda or Lambda Bar" flag:
				// (do notice that this does not depend on the particle being a carbon copy nor being final!
				// If a carbon copy was a Lambda, it still counts, and no Lambda should ever be final in this
				// simulation where I allow it to decay!)
			if (particle_PID == 3122 || particle_PID == -3122){contains_lambda_or_lambdabar = true;}

			int motherIdx1 = pythia.event[i].mother1(); // Doing this, I can access it only once and don't need to re-read this information for pions!
			int motherIdx2 = pythia.event[i].mother2();

			////////////////////////////////////////////////////////////////////////////////////////
			/// CAUTION! THIS PRIMARY DEFINITION IS ACTUALLY A LITTLE BIT WRONG!!!
			/// Check the Lambda isExperimentalPrimary definition below!!!
			// Checking if the particle would be considered primary or not
				// The particle will be considered primary by default. If it has a mother with a long lifetime, then it will be considered secondary.
				// This allows for some contamination from possible secondaries that have long decay showers, though. My hypothesis is that this effect is too small to be significant.
				// Possible secondaries from 2 -> N processes and two partons becoming a single particle will be considered as primaries.
				// A little bit of a MC cheat is that particles whose lifetime is smaller than 10 mm/c but came straight from the collision will be considered primary, 
				// and in experimental data that would not be true.
			Bool_t IsPrimary_value = true;
			// First, checking if it came from an actual decay:
				// (I want to count all particles that didn't fall in both the loop's conditions with the Ntracks_primary, not the ones that checked one but didn't check the other if)
			if (motherIdx1 > 0 && motherIdx2 == 0){
				if (pythia.event[motherIdx1].tau0() >= 10.){ // Tau0 in mm/c. Reference is slide 112 of https://indico.cern.ch/event/666222/contributions/2768780/attachments/1551303/2437229/DPG_AnalysisTutorial_20171102.pdf
					IsPrimary_value = false;
				}
				else{
					Ntracks_primary += 1;
				}
			}
			else{
				Ntracks_primary += 1; // This does not overcount the above check: they never happen simultaneuosly for the same particle!
			}
			////////////////////////////////////////////////////////////////////////////////////////

			if (isfinal && isCharged){
				// event_hist_charged->Fill((double) pT);
				pT_hist_charged_final->Fill((double) pT);
				n_charged_event += 1;
			}

			if (isfinal && (particle_PID == 211 or particle_PID == (-1)*211)){ // Could just use idAbs() == 211 too!
				pT_hist_pion_final->Fill((double) pT);
				n_pion_event += 1;
			}
			else if (isfinal && (particle_PID == 2212 or particle_PID == (-1)*2212)){ // Proton
				// event_hist_proton->Fill((double) pT);
				pT_hist_proton_final->Fill((double) pT);
				n_proton_event += 1;
			}
			else if (isfinal && (particle_PID == 321 or particle_PID == (-1)*321)){ // Kaon
				// event_hist_kaon->Fill((double) pT);
				pT_hist_kaon_final->Fill((double) pT);
				n_kaon_event += 1;
			}

			// Including all USEFUL particles into the arrays for storage:
			// First, should look only for protons and anti-protons:
			if (isfinal){
				if (particle_PID == 2212 || particle_PID == -2212){
					// After that, should check if the mother of this proton or anti-proton is a Lambda or LambdaBar
					Particle mother_particle_object = pythia.event[motherIdx1]; // Already using the Pythia8 namespace
					int mother_PID = mother_particle_object.id();
					if (mother_PID == 3122 || mother_PID == -3122){
						// Should not worry with this Lambda being a carbon-copy or not: the useful kinematic information comes from the actual
						// mother of the proton, not the previous stages of when that Lambda was generated, as far as I know. If we had the QGP,
						// what would be measured is not the "first" Lambda, but the final after-interactions Lambda's daughter. In iSS, the
						// Lambda would be the result of all the quark interactions inside the QGP, so it would also come after these interactions
						// that are responsible for creating carbon-copies!
						
						Bool_t isProton_not_PBar_value = false; // A bool to know if this particle is a proton or a PBar (they will be stored in the same branch)
						if (particle_PID == 2212){
							isProton_not_PBar_value = true;
						}

						Bool_t isLambda_notLBar_value = false; // A bool to know if this particle is a Lambda or an LBar (they will be stored in the same branch)
						if (mother_PID == 3122){
							isLambda_notLBar_value = true;
						}

						//////////////////////////////////////////////////////////
						// Checking if the Lambda is primary:
						Bool_t isExperimentalPrimary = true;
						while ((mother_particle_object.mother1() == mother_particle_object.mother2()) && (mother_particle_object.mother1() > 0)){
							// Should keep changing the mother particle definition if the mother is a carbon copy!
							// We want to see if the physical mother of the Lambda has a specific tau0 lifetime, not the lifetime of the carbon-copy mother!
							mother_particle_object = pythia.event[mother_particle_object.mother1()];
						}

							// Getting the mother of the Lambda:
						Particle lambda_mother_object = pythia.event[mother_particle_object.mother1()];
						// Case 1: Check lifetime of the mother
						// (long-lived weak decays should be considered non-primary)
						// (Definition of experimental primary in slide 112 of https://indico.cern.ch/event/666222/contributions/2768780/attachments/1551303/2437229/DPG_AnalysisTutorial_20171102.pdf)
						// If the mother of this lambda is the beam itself, there is no 
						if (lambda_mother_object.tau0() >= 10.0){ // tau0 in mm/c
							isExperimentalPrimary = false;
						}

						// Case 2: Direct hadronization products (fragmentation or R-hadron formation)
						// (If the Lambda's mother comes from any of these statuses, then it came from the beam's strings in some sense, so should be set as primary again)
						if ((abs(lambda_mother_object.status()) >= 81 && abs(lambda_mother_object.status()) <= 86) || (abs(lambda_mother_object.status()) >= 101 && abs(lambda_mother_object.status()) <= 106)){
							isExperimentalPrimary = true;
						}
						//////////////////////////////////////////////////////////

						// Saving the useful variables:
						ProtonOrPBar_pt.push_back(pT);
						ProtonOrPBar_y.push_back(rapidity);
						ProtonOrPBar_Phi.push_back(Phi);
						ProtonOrPBar_m.push_back(mass);

						LambdaOrLBar_pt.push_back(pythia.event[motherIdx1].pT());
						LambdaOrLBar_y.push_back(pythia.event[motherIdx1].y());
						LambdaOrLBar_Phi.push_back(pythia.event[motherIdx1].phi());
						LambdaOrLBar_m.push_back(pythia.event[motherIdx1].m());

						isProton_not_PBar.push_back(isProton_not_PBar_value);
						isLambda_notLBar.push_back(isLambda_notLBar_value);

						isLambdaOrLBarExperimentalPrimary.push_back(isExperimentalPrimary);
					}
				}
			}

			// Do not need to store these variables anymore -- Could make some std::map objects to at least know the PIDs that are produced, but that is not necessary right now!
			// // ID[i] = particle_PID;
			// ID.push_back(particle_PID); // New C++ style vector storing
			
			// IsCarbonCopy.push_back(IsCarbonCopy_value);

			// mother1_ID.push_back(motherIdx1); // Not quite the PID of the mother. It is the Idx of the mother in Pythia's list of particles
			// mother2_ID.push_back(motherIdx2);

			// IsFinal.push_back(isfinal);
			// // px.push_back(pythia.event[i].px());
			// // py.push_back(pythia.event[i].py());
			// // pz.push_back(pythia.event[i].pz());

			// pt.push_back(pT);
			// m.push_back(mass); // The actual simulated mass of the particle in the dynamic, so may be somewhat different from the PDG value. Thus the need to store it separately.
			// // e.push_back(pythia.event[i].e());
			
			// Phi.push_back(Phi);
			// // Eta.push_back(eta_rap);
			// y.push_back(rapidity);

			// // charged.push_back(isCharged); // Don't actually need to save this! This is purely a PDG-lookup information, so no need to store it when the PID is already enough to know it!

			// IsPrimary.push_back(IsPrimary_value);



			// Finally, introducing this particle to the jet candidate list if a certain number of criteria are obeyed:
			if (isfinal && isCharged && !IsCarbonCopy_value && isvisible){ // Particle should be physical and a final state only! Should also be visible and CHARGED only. Notice that 
				// Create a PseudoJet from the complete Pythia particle:
      			fastjet::PseudoJet particleTemp = pythia.event[i];
				// Store the original Pythia index in "user_index" to later on get the PID of the particle
  				particleTemp.set_user_index(i);

				// Store the PseudoJet:
				FastJetInputs.push_back(particleTemp);
			}
		}
		///////////////////////////////////////////////////////////
		// First, filling the centrality estimators to have an unbiased estimation
		// (If we defined centrality from jet-only events, we would have higher centrality events
		// by construction, wouldn't we?)
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

		// Now filling the histograms for each particle type:
		if (n_charged_event != 0){hEventCounterCharged->Fill(0);}
		if (n_pion_event != 0){hEventCounterPion->Fill(0);}
		if (n_proton_event != 0){hEventCounterProton->Fill(0);}
		if (n_kaon_event != 0){hEventCounterKaon->Fill(0);}
		///////////////////////////////////////////////////////////


		// Now applying the jet cuts:
		if(contains_lambda_or_lambdabar){hEventCounterWithLambdaOrBar->Fill(0);}

		// Checking the Jet information to see if this event will be excluded or not:
			// Run Fastjet algorithm and sort jets in pT order.
		vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
		fastjet::ClusterSequence clustSeq(FastJetInputs, jetDef); // Clusters the jet candidates
		inclusiveJets = clustSeq.inclusive_jets(jet_min_pT); // Gets all the jets with pT above this minimum threshold
		sortedJets = fastjet::sorted_by_pt(inclusiveJets); // Sort from highest to lowest pT

		hJetCounterPerEvent->Fill(sortedJets.size()); // This may be zero without a problem! I just need to know how many of the total events have a jet
		if(sortedJets.size() != 0){hEventCounterWithJets_pTleqJetMinPt->Fill(0);}

			// Now applying the eta cut for the jets (if they don't pass |eta|<jet_max_eta, then the jet is not selected)
		fastjet::PseudoJet leadingJet;
		fastjet::PseudoJet firstDiscardedJet;
		Bool_t foundJet = false; // A bool to verify if there is any jet with pT > jet_min_pT AND |eta| < 1.0 (or jet_max_eta, to be more general)
		int discarded_highpT_jets = 0;
		for (const auto& jet : sortedJets){
			// A pre-selection that is necessary to pick leading Jets that actually point in the midrapidity region, 
			// properly using events even if we have to pick a leading jet that does not have the largest pT of the 
			// event! In other words, we can discard the largest pT jet of the event in some cases, to pick only the
			// most transverse jet.
			hLeadJetY_DiscardedUntilLead_pTleqJetMinPt->Fill(jet.rap()); // Will fill for all the jets, ordered by pT, even if they are not the selected leading jet of this event!
			hLeadJetPt_DiscardedUntilLead_pTleqJetMinPt->Fill(jet.pt());
			hLeadJetPtY_DiscardedUntilLead_pTleqJetMinPt->Fill(jet.pt(), jet.rap());
			hLeadJetPhi_DiscardedUntilLead_pTleqJetMinPt->Fill(jet.phi_std());
			
			// Will update the leadingJet variable only once, at the first viable jet.
			// Since it is all sorted by pT, selects the first jet with |\eta|<\eta_max,
			// searching from the highest-pT jet of the event to the lowest one.
			// (Rewrote this loop to loop over all of the jets of the event, instead
			// of breaking at the first valid leadingJet)
			if (!foundJet){
				if (std::abs(jet.eta()) < jet_max_eta){
					leadingJet = jet;
					foundJet = true;
					hDiscardedHighPtJetCounterPerEvent->Fill(discarded_highpT_jets); // Fills the histogram with the counter of jets discarded until the desired one was found
					// The number of times this was filled can also show us how many events failed to have any jet in the |eta|<jet_max_eta region, if compared with the
					// hEventCounterWithJets_pTleqJetMinPt_EtaCuts counter
				}
				else{
					// Should only increase if the jet was not yet found.
					// It wouldn't make much difference if this counter was increased outside of the "else" statement, but this way it avoids an extra useless sum
					discarded_highpT_jets += 1; // A counter to know how many jets were thrown away until the desired eta-range jet was found
					
					// Getting the information of the first discarded jet (useful to know if it was back-to-back with the selected leadingJet!):
					if (discarded_highpT_jets == 1){firstDiscardedJet = jet;}
				}
			}

			std::vector<fastjet::PseudoJet> jet_constituents = jet.constituents();
			hNParticlesPerJet->Fill(jet_constituents.size());

			// Getting the PID of each constituent of this jet:
			for (const auto& jet_constituent : jet_constituents){
				int particle_PythiaIdx = jet_constituent.user_index(); // Using what we defined in particleTemp
				int particle_PID = pythia.event[particle_PythiaIdx].id();
				mapCountsOfPIDAllJetConstituents[particle_PID]++; // Creates a key with the PID number, and then adds 1 to that key's value
			}
		}
		if (!foundJet){ // Should restart the whole loop, as we are not saving events that don't have the minimum jets necessary for my analysis!
			iEvent--;
			continue;
		}
		hEventCounterWithJets_pTleqJetMinPt_EtaCuts->Fill(0); // Notice you don't need an "else" statement here

		std::vector<fastjet::PseudoJet> leadingJet_constituents = leadingJet.constituents();
		hNParticlesLeadingJet_EtaCuts->Fill(leadingJet_constituents.size());

		// Creating a leadingJet map for the current event -- No need to do the jet_constituent loop for each more restrictive cut!
		// In other words, loops through the jets only once, then just reads from the map for this event's leading jet
		std::map<int,int> CurrentEvent_LeadingJet_map;
		for (const auto& jet_constituent : leadingJet_constituents){
			// Getting the PID of each constituent of this jet:
			int particle_PythiaIdx = jet_constituent.user_index(); // Using what we defined in particleTemp
			int particle_PID = pythia.event[particle_PythiaIdx].id();
			CurrentEvent_LeadingJet_map[particle_PID]++; // Creates a key with the PID number, and then adds 1 to that key's value

			// Can fill this particular Eta cut's map right now:
			mapCountsOfPIDLeadingJetConstituents_EtaCuts[particle_PID]++;
		}

		// Setting the pT, m, y, phi values for the jet before filling the tree:		
		pt_jet = leadingJet.pt();
		m_jet = leadingJet.m();
		y_jet = leadingJet.rap();
		Phi_jet = leadingJet.phi_std(); // Used phi_std instead of phi(), to stay in agreement with Pythia's phi angles from -PI to PI

		hLeadJetY_pTleqJetMinPt_eta_cut->Fill(y_jet);
		hLeadJetPt_pTleqJetMinPt_eta_cut->Fill(pt_jet);
		hLeadJetPtY_pTleqJetMinPt_eta_cut->Fill(pt_jet, y_jet);
		hLeadJetPhi_pTleqJetMinPt_eta_cut->Fill(Phi_jet);

		// Saving the information of the first discarded jet's relative position to the leading jet:
		Float_t pt_firstDiscardedJet = firstDiscardedJet.pt();
		Float_t m_firstDiscardedJet = firstDiscardedJet.m();
		Float_t y_firstDiscardedJet = firstDiscardedJet.rap();
		Float_t Phi_firstDiscardedJet = firstDiscardedJet.phi_std();

		hfirstDiscardedJetY_minus_leadingJet_MinPtCut->Fill(y_firstDiscardedJet - y_jet);
		hfirstDiscardedJetPt_MinPtCut->Fill(pt_jet);
		hfirstDiscardedJetPtY_minus_leadingJet_MinPtCut->Fill(pt_jet, y_jet - y_jet);
		hfirstDiscardedJetPhi_minus_leadingJet_MinPtCut->Fill(wrapToInterval(Phi_firstDiscardedJet - Phi_jet, -PI, PI));

		if(!contains_lambda_or_lambdabar){ // I am only going to save events that have jets in the interval which is useful to my Ring analysis!
			iEvent--;
			continue;
		}
		hEventCounterWithLambdaOrBarWithJets_pTleqJetMinPt_EtaCuts->Fill(0);
		hNParticlesLeadingJet_EtaCutsWithLambda->Fill(leadingJet_constituents.size());
			// Filling the map for this "EtaCutsWithLambda" cut:
		for (const auto &key_and_value_pair : CurrentEvent_LeadingJet_map){
			int PID = key_and_value_pair.first; // the key
			int counts = key_and_value_pair.second; // the count
			mapCountsOfPIDLeadingJetConstituents_EtaCutsWithLambda[PID] += counts;
		}

		hJetProxyY_pTleqJetMinPt_eta_cut_WithLambdaOrBar->Fill(y_jet);
		hJetProxyPt_pTleqJetMinPt_eta_cut_WithLambdaOrBar->Fill(pt_jet);
		hJetProxyPtY_pTleqJetMinPt_eta_cut_WithLambdaOrBar->Fill(pt_jet, y_jet);
		hJetProxyPhi_pTleqJetMinPt_eta_cut_WithLambdaOrBar->Fill(Phi_jet);

		// The check below is just paranoid -- For a high enough pT, this should not be necessary!
		// (just removed the iEvent-- and the skipping part of it, while preserving the histogramming)
		if (std::fabs(leadingJet.rap()) < 0.5){ // If the event contains a lambda, we will do a final, less restrictive check: does this event's leading jet stay in the |y|<0.5 region?
			hEventCounterPassedYExtraCheck->Fill(0);
		}
		// else{
		// 	iEvent--;
		// 	continue;
		// }
		hNParticlesLeadingJet_PassedYExtraCheck->Fill(leadingJet_constituents.size());
			// Filling the map for this "EtaCutsWithLambda" cut:
		for (const auto &key_and_value_pair : CurrentEvent_LeadingJet_map){
			int PID = key_and_value_pair.first; // the key
			int counts = key_and_value_pair.second; // the count
			mapCountsOfPIDLeadingJetConstituents_PassedYExtraCheck[PID] += counts;
		}

		hJetProxyY_YExtraCheck->Fill(y_jet);
		hJetProxyPt_YExtraCheck->Fill(pt_jet);
		hJetProxyPtY_YExtraCheck->Fill(pt_jet, y_jet);
		hJetProxyPhi_YExtraCheck->Fill(Phi_jet);

		t3->Fill(); // Filling the tree first thing after the loop!
	}

	// Creating histograms for the PIDs contained in each specific jet subset:
	auto hPIDOfAllJetConstituents = new TH1D("hPIDOfAllJetConstituents", "PIDs of all jets' constituents", mapCountsOfPIDAllJetConstituents.size(), 0, mapCountsOfPIDAllJetConstituents.size());
	auto hPIDLeadingJetConstituents_EtaCuts = new TH1D("hPIDLeadingJetConstituents_EtaCuts", "PIDs of eta_cut_OK jets' constituents", mapCountsOfPIDLeadingJetConstituents_EtaCuts.size(), 0, mapCountsOfPIDLeadingJetConstituents_EtaCuts.size());
	auto hPIDLeadingJetConstituents_EtaCutsWithLambda = new TH1D("hPIDLeadingJetConstituents_EtaCutsWithLambda", "PIDs of eta_cut_OK+w/LambdaOrLBar jets' constituents", mapCountsOfPIDLeadingJetConstituents_EtaCutsWithLambda.size(), 0, mapCountsOfPIDLeadingJetConstituents_EtaCutsWithLambda.size());
	auto hPIDLeadingJetConstituents_PassedYExtraCheck = new TH1D("hPIDLeadingJetConstituents_PassedYExtraCheck", "PIDs of Useful(w/LambdaOrLBar,Ycut_OK) jets' constituents", mapCountsOfPIDLeadingJetConstituents_PassedYExtraCheck.size(), 0, mapCountsOfPIDLeadingJetConstituents_PassedYExtraCheck.size());
		// Filling these histograms:
	FillHistogramFromMap(hPIDOfAllJetConstituents, mapCountsOfPIDAllJetConstituents);
	FillHistogramFromMap(hPIDLeadingJetConstituents_EtaCuts, mapCountsOfPIDLeadingJetConstituents_EtaCuts);
	FillHistogramFromMap(hPIDLeadingJetConstituents_EtaCutsWithLambda, mapCountsOfPIDLeadingJetConstituents_EtaCutsWithLambda);
	FillHistogramFromMap(hPIDLeadingJetConstituents_PassedYExtraCheck, mapCountsOfPIDLeadingJetConstituents_PassedYExtraCheck);

	// pythia.stat();
	// t3->Print();

	f.cd();

	std::cout << "Writing information into the .root file" << std::endl;
	t3->Write();

	hEventCounter->Write();
	hEventCounterWithLambdaOrBar->Write();
	hEventCounterWithJets_pTleqJetMinPt->Write();
	hEventCounterWithJets_pTleqJetMinPt_EtaCuts->Write();
	hEventCounterWithLambdaOrBarWithJets_pTleqJetMinPt_EtaCuts->Write();
	hEventCounterPassedYExtraCheck->Write();
	hINELEventCounter->Write();
	hEventCounterCharged->Write();
	hEventCounterPion->Write();
	hEventCounterProton->Write();
	hEventCounterKaon->Write();
	
	hJetCounterPerEvent->Write();
	hDiscardedHighPtJetCounterPerEvent->Write();
	hNParticlesPerJet->Write();
	hNParticlesLeadingJet_EtaCuts->Write();
	hNParticlesLeadingJet_EtaCutsWithLambda->Write();
	hNParticlesLeadingJet_PassedYExtraCheck->Write();

	hPIDOfAllJetConstituents->Write();
	hPIDLeadingJetConstituents_EtaCuts->Write();
	hPIDLeadingJetConstituents_EtaCutsWithLambda->Write();
	hPIDLeadingJetConstituents_PassedYExtraCheck->Write();

	pT_hist_charged_final->Sumw2();
	pT_hist_pion_final->Sumw2();
	pT_hist_proton_final->Sumw2();
	pT_hist_kaon_final->Sumw2();

	pT_hist_charged_final->Write();
	pT_hist_pion_final->Write();
	pT_hist_proton_final->Write();
	pT_hist_kaon_final->Write();

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

	hINELev_Ntracks_final_charged_forward->Write();
	hINELev_Ntracks_final_charged_forward_Charged->Write();
	hINELev_Ntracks_final_charged_forward_Pion->Write();
	hINELev_Ntracks_final_charged_forward_Proton->Write();
	hINELev_Ntracks_final_charged_forward_Kaon->Write();

	hLeadJetY_DiscardedUntilLead_pTleqJetMinPt->Write();
	hLeadJetPt_DiscardedUntilLead_pTleqJetMinPt->Write();
	hLeadJetPtY_DiscardedUntilLead_pTleqJetMinPt->Write();
	hLeadJetPhi_DiscardedUntilLead_pTleqJetMinPt->Write();

	hLeadJetY_pTleqJetMinPt_eta_cut->Write();
	hLeadJetPt_pTleqJetMinPt_eta_cut->Write();
	hLeadJetPtY_pTleqJetMinPt_eta_cut->Write();
	hLeadJetPhi_pTleqJetMinPt_eta_cut->Write();

	hfirstDiscardedJetY_minus_leadingJet_MinPtCut->Write();
	hfirstDiscardedJetPt_MinPtCut->Write();
	hfirstDiscardedJetPtY_minus_leadingJet_MinPtCut->Write();
	hfirstDiscardedJetPhi_minus_leadingJet_MinPtCut->Write();

	hJetProxyY_pTleqJetMinPt_eta_cut_WithLambdaOrBar->Write();
	hJetProxyPt_pTleqJetMinPt_eta_cut_WithLambdaOrBar->Write();
	hJetProxyPtY_pTleqJetMinPt_eta_cut_WithLambdaOrBar->Write();
	hJetProxyPhi_pTleqJetMinPt_eta_cut_WithLambdaOrBar->Write();

	hJetProxyY_YExtraCheck->Write();
	hJetProxyPt_YExtraCheck->Write();
	hJetProxyPtY_YExtraCheck->Write();
	hJetProxyPhi_YExtraCheck->Write();

	f.Close();

	// Deleting all objects instatiated in heap:
		// Apparently ROOT already deletes them -- Deleting again made the program crash!
	// delete hEventCounter;
	// delete hINELEventCounter;
	// delete hEventCounterCharged;
	// delete hEventCounterPion;
	// delete hEventCounterProton;
	// delete hEventCounterKaon;

	// delete t3;

	// delete pT_hist_charged_final;
	// delete pT_hist_pion_final;
	// delete pT_hist_proton_final;
	// delete pT_hist_kaon_final;

	// delete hNtracks;
	// delete hNtracks_charged_forward;
	// delete hNtracks_primary;
	// delete hNtracks_final;
	// delete hNtracks_final_forward;
	// delete hNtracks_final_forward_INEL_l0;
	// delete hNtracks_final_charged_center;
	// delete hNtracks_final_charged;
	// delete hNtracks_final_charged_INEL_l0;
	// delete hNtracks_final_charged_forward;
	// delete hNtracks_final_charged_forward_INEL_l0;

	// delete hINELev_Ntracks_final_charged_forward;
	// delete hINELev_Ntracks_final_charged_forward_Charged;
	// delete hINELev_Ntracks_final_charged_forward_Pion;
	// delete hINELev_Ntracks_final_charged_forward_Proton;
	// delete hINELev_Ntracks_final_charged_forward_Kaon;

	return;
}

////////////////////////////////
/// Physics helper functions ///
////////////////////////////////
void doCentrality(std::string output_folder, int N_cores, std::string input_card_name){
	std::cout << "\nProcessing centrality" << std::endl;
	int multiplicity_nbins = 50000; // The number of bins is also the number of entries, which gives one bin per possible multiplicity!

	std::string output_centrality_filename =(std::string) output_folder;
	output_centrality_filename += "/Ntracks_to_centrality_conversion.root";

	TFile output_file(output_centrality_filename.c_str(), "RECREATE");

	// Reading all N_cores worker data and summing the hNtracks histograms:
		// First declaring histograms to receive the sums:
	TH1D *hNtracks_sum = new TH1D("hNtracks_sum", "hNtracks_sum", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hNtracks_charged_forward_sum = new TH1D("hNtracks_charged_forward_sum", "hNtracks_charged_forward_sum", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hNtracks_primary_sum = new TH1D("hNtracks_primary_sum", "hNtracks_primary_sum", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hNtracks_final_sum = new TH1D("hNtracks_final_sum", "hNtracks_final_sum", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hNtracks_final_forward_sum = new TH1D("hNtracks_final_forward_sum", "hNtracks_final_forward_sum", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hNtracks_final_forward_INEL_l0_sum = new TH1D("hNtracks_final_forward_INEL_l0_sum", "hNtracks_final_forward_INEL_l0_sum", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hNtracks_final_charged_center_sum = new TH1D("hNtracks_final_charged_center_sum", "hNtracks_final_charged_center_sum", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hNtracks_final_charged_sum = new TH1D("hNtracks_final_charged_sum", "hNtracks_final_charged_sum", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hNtracks_final_charged_INEL_l0_sum = new TH1D("hNtracks_final_charged_INEL_l0_sum", "hNtracks_final_charged_INEL_l0_sum", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hNtracks_final_charged_forward_sum = new TH1D("hNtracks_final_charged_forward_sum", "hNtracks_final_charged_forward_sum", multiplicity_nbins, 0, multiplicity_nbins);
	TH1D *hNtracks_final_charged_forward_INEL_l0_sum = new TH1D("hNtracks_final_charged_forward_INEL_l0_sum", "hNtracks_final_charged_forward_INEL_l0_sum", multiplicity_nbins, 0, multiplicity_nbins);

    for (int WorkerId = 0; WorkerId < N_cores; ++WorkerId){
        std::string data_filename = (std::string) output_folder;
		data_filename += "/output_core" + to_string(WorkerId) + ".root";

		TFile *data_file = TFile::Open(data_filename.c_str(), "READ");

			// Fetching from current file:
		TH1D *hNtracks_local = dynamic_cast<TH1D*>(data_file->Get("hNtracks"));
		TH1D *hNtracks_charged_forward_local = dynamic_cast<TH1D*>(data_file->Get("hNtracks_charged_forward"));
		TH1D *hNtracks_primary_local = dynamic_cast<TH1D*>(data_file->Get("hNtracks_primary"));
		TH1D *hNtracks_final_local = dynamic_cast<TH1D*>(data_file->Get("hNtracks_final"));
		TH1D *hNtracks_final_forward_local = dynamic_cast<TH1D*>(data_file->Get("hNtracks_final_forward"));
		TH1D *hNtracks_final_forward_INEL_l0_local = dynamic_cast<TH1D*>(data_file->Get("hNtracks_final_forward_INEL_l0"));
		TH1D *hNtracks_final_charged_center_local = dynamic_cast<TH1D*>(data_file->Get("hNtracks_final_charged_center"));
		TH1D *hNtracks_final_charged_local = dynamic_cast<TH1D*>(data_file->Get("hNtracks_final_charged"));
		TH1D *hNtracks_final_charged_INEL_l0_local = dynamic_cast<TH1D*>(data_file->Get("hNtracks_final_charged_INEL_l0"));
		TH1D *hNtracks_final_charged_forward_local = dynamic_cast<TH1D*>(data_file->Get("hNtracks_final_charged_forward"));
		TH1D *hNtracks_final_charged_forward_INEL_l0_local = dynamic_cast<TH1D*>(data_file->Get("hNtracks_final_charged_forward_INEL_l0"));

			// Adding all histograms:
		hNtracks_sum->Add(hNtracks_local);
		hNtracks_charged_forward_sum->Add(hNtracks_charged_forward_local);
		hNtracks_primary_sum->Add(hNtracks_primary_local);
		hNtracks_final_sum->Add(hNtracks_final_local);
		hNtracks_final_forward_sum->Add(hNtracks_final_forward_local);
		hNtracks_final_forward_INEL_l0_sum->Add(hNtracks_final_forward_INEL_l0_local);
		hNtracks_final_charged_center_sum->Add(hNtracks_final_charged_center_local);
		hNtracks_final_charged_sum->Add(hNtracks_final_charged_local);
		hNtracks_final_charged_INEL_l0_sum->Add(hNtracks_final_charged_INEL_l0_local);
		hNtracks_final_charged_forward_sum->Add(hNtracks_final_charged_forward_local);
		hNtracks_final_charged_forward_INEL_l0_sum->Add(hNtracks_final_charged_forward_INEL_l0_local);

        data_file->Close();
    }

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


////////////////////////
/// Helper functions ///
////////////////////////
void FillHistogramFromMap(TH1* hist, const std::map<int,int>& map){ // A helper function to make the whole process of filling maps cleaner
    int bin = 1;
    for (const auto &kv : map){
		// Fixing the label to the actual PID number:
        hist->GetXaxis()->SetBinLabel(bin, std::to_string(kv.first).c_str());
		// Setting the content:
        hist->SetBinContent(bin, kv.second);
        bin++;
    }
	hist->ResetStats(); // Have to do this to make the histogram actually have the number of particles as the NEntries, not the number of bins/fills we did in this function!
}

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


/////////////////////////////////////////
/////////////// MAIN CALL ///////////////
/////////////////////////////////////////
int main(int argc, char *argv[]){
    if (argc != 5){
		std::cerr << "Usage:" << argv[0] << " output_folder input_card_path N_ev N_cores" << std::endl;
        return 1; // Indicate an error.
    }

    const char *output_folder = (const char*) argv[1];
	const char *input_card_path = (const char*) argv[2];
    double N_ev_receiver = atof(argv[3]); // This receives input like 1e9 and converts it into a proper double.
	uint N_ev = static_cast<int>(N_ev_receiver); // Should use long long for real high statistics (even uint would overflow at 1e10 or so!), but this is enough for now
	uint N_cores = std::atoi(argv[4]);

	std::cout << "\n\nNow running " << argv[0] << " " << argv[1] << " " << argv[2] << " " << N_ev << " " << argv[4] << std::endl;

	namespace fs = std::filesystem;
	fs::path p(input_card_path);
	std::string input_card_name = p.stem().string(); // stem = filename without extension

	// Creating a folder to contain all .root files generated by the workers, based on the total number of events to generate:
	std::string output_folder_name = output_folder;
	output_folder_name += "/" + input_card_name + "_" + std::to_string(N_ev) + "ev";
	fs::path output_folder_for_current_card = output_folder_name;
	
		// Create folder (including parents if needed)
	std::error_code ec;
	if (!fs::create_directories(output_folder_for_current_card, ec) && ec){
		std::cerr << "Warning: could not create directory " << output_folder_for_current_card << " (" << ec.message() << ")\n";
	}

	ROOT::EnableThreadSafety(); // Makes ROOT thread-aware -- This solves the parallelization problems and initializes ROOT, solving all previous problems!
	omp_set_num_threads(N_cores); // Simple and short parallelization!
    #pragma omp parallel for
	for (int WorkerId = 0; WorkerId < N_cores; WorkerId++){
		// Calculating the number of events per worker, in a way that will give me exactly N_ev for whichever number of workers I use:
		int base = N_ev / N_cores;
    	int remainder = N_ev % N_cores;
		int N_ev_current_worker = base + (WorkerId < remainder ? 1 : 0);

		// std::cout << base << std::endl;
		// std::cout << remainder << std::endl;
		// std::cout << N_ev_current_worker << std::endl;

		RunWorker(WorkerId, N_ev_current_worker, output_folder_for_current_card, input_card_path, input_card_name);
	}

	std::cout << "\nPythia parallel loop ended!" << std::endl;

	// Processing centrality for all files:
	doCentrality(output_folder_for_current_card, N_cores, input_card_name);
	
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
void multiplicity_to_centrality(TH1D *multiplicity_hist, TH1D *centrality_hist){
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
