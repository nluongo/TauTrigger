using namespace std;
#include <iostream>

// This script requires functions found in the TauEtFuncs.cpp script

void TauEtStudy() {
    // May need to specify an absolute path if not running ROOT from command line
    // Open file containing background, hardcode the filename for now
    TFile fback("output_MB80.root");
	// Open file containing signal
	TFile fsig("output_Z80.root");
    
    // Get a pointer to the background TTree (name must match)
    TTree* tback = (TTree*)fback.Get("mytree");
	// Get a pointer to the signal TTree
	TTree* tsig = (TTree*)fsig.Get("mytree");
    
    // Arrays, dimensions must match size in ntuple
	// Define background variables to read in
	Float_t L0_Back_Et[3][3];
	Float_t L1_Back_Et[13][3];
	Float_t L2_Back_Et[13][3];
	Float_t L3_Back_Et[3][3];
	Float_t Had_Back_Et[3][3];
	// Background object
	TLorentzVector* mcseed = new TLorentzVector();

	// Define signal variables to read in
	Float_t L0_Sig_Et[3][3];
	Float_t L1_Sig_Et[13][3];
	Float_t L2_Sig_Et[13][3];
	Float_t L3_Sig_Et[3][3];
	Float_t Had_Sig_Et[3][3];
	// Signal object
	TLorentzVector* mctau = new TLorentzVector();
	Int_t numofpi0;	// The number of Pi0 particles involved in an event
	Int_t zeropi0num = 0;	// The number of events with no Pi0 particles involved
	Int_t onepi0num = 0;	// The number of events with one Pi0 particles involved
	Int_t twopi0num = 0;	// The number of events with two Pi0 particles involved
	Int_t threepi0num = 0;	// The number of events with three Pi0 particles involved
    
    // Make assignemnt between tree and declared variables for background
	tback->SetBranchAddress("L0CellEt[3][3]", &L0_Back_Et[0][0]);  // Assigns to address of first array element
	tback->SetBranchAddress("L1CellEt[13][3]", &L1_Back_Et[0][0]);  // Assigns to address of first array element
    tback->SetBranchAddress("L2CellEt[13][3]", &L2_Back_Et[0][0]);  // Assigns to address of first array element
	tback->SetBranchAddress("L3CellEt[3][3]", &L3_Back_Et[0][0]);  // Assigns to address of first array element
	tback->SetBranchAddress("HadCellEt[3][3]", &Had_Back_Et[0][0]);  // Assigns to address of first array element
	tback->SetBranchAddress("seed", &mcseed); // Assigns to address of first array element

	// Make assignment between tree and declared variables for signal
	tsig->SetBranchAddress("L0CellEt[3][3]", &L0_Sig_Et[0][0]);  // Assigns to address of first array element
	tsig->SetBranchAddress("L1CellEt[13][3]", &L1_Sig_Et[0][0]);  // Assigns to address of first array element
	tsig->SetBranchAddress("L2CellEt[13][3]", &L2_Sig_Et[0][0]);  // Assigns to address of first array element
	tsig->SetBranchAddress("L3CellEt[3][3]", &L3_Sig_Et[0][0]);  // Assigns to address of first array element
	tsig->SetBranchAddress("HadCellEt[3][3]", &Had_Sig_Et[0][0]);  // Assigns to address of first array element
	tsig->SetBranchAddress("mc_visibleTau", &mctau);
	tsig->SetBranchAddress("mc_pi0", &numofpi0);

	// Signal variables
	Int_t sigentries = (Int_t)tsig->GetEntries();
	Float_t L0_AverageCell_Sig_Et[3][3] = { 0 };
	Float_t L0_AverageFlipCell_Sig_Et[3][3] = { 0 };
	Float_t L1_AverageCell_Sig_Et[13][3] = { 0 };
	Float_t L1_AverageFlipCell_Sig_Et[13][3] = { 0 };
	Float_t L2_AverageCell_Sig_Et[13][3] = { 0 };
	Float_t L2_AverageFlipCell_Sig_Et[13][3] = { 0 };
	Float_t L3_AverageCell_Sig_Et[3][3] = { 0 };
	Float_t L3_AverageFlipCell_Sig_Et[3][3] = { 0 };
	Float_t Had_AverageCell_Sig_Et[3][3] = { 0 };
	Float_t Had_AverageFlipCell_Sig_Et[3][3] = { 0 };
	Float_t TrueEt[sigentries];
	Float_t Reco1_Sig_Et[sigentries]; // Reconstructed signal energy using the following substructures (1x1, 5x3, 5x3, 3x3, 3x3)
	Float_t Reco2_Sig_Et[sigentries]; // Reconstructed signal energy using the following substructures (3x3, 7x3, 5x3, 3x3, 3x3)
	Float_t Reco3_Sig_Et[sigentries]; // Reconstructed signal energy using the following substructures (1x2, 5x2, 5x2, 1x2, 3x3). Done post-pre-processing
	Float_t Reco3A_Sig_Et[sigentries]; // Reconstructed signal energy using the following substructures (1x3, 5x3, 5x3, 1x3, 3x3). Done post-pre-processing
	Float_t Reco3B_Sig_Et[sigentries]; // Reconstructed signal energy using the following substructures (1x3, 5x3, 5x3, 1x3, 3x3). Done without pre-processing
	Float_t Reco4_Sig_Et[sigentries]; // Reconstructed signal energy using the following substructures (1x2, 5x2, 5x2, 1x2, None). Done post-pre-processing
	Float_t Reco5_Sig_Et[sigentries]; // Reconstructed signal energy using the following substructures (1x2, 5x2, 5x2, None, None). Done post-pre-processing
	Float_t Reco6_Sig_Et[sigentries]; // Reconstructed signal energy using the following substructures (None, 5x2, 5x2, None, None). Done post-pre-processing
	Float_t RecoAll_Sig_Et[sigentries];
	Float_t L1_AverageLayer_0PiZero_Et = 0;
	Float_t L2_AverageLayer_0PiZero_Et = 0;
	Float_t L1_AverageLayer_1PiZero_Et = 0;
	Float_t L2_AverageLayer_1PiZero_Et = 0;
	Float_t L1_AverageLayer_2PiZero_Et = 0;
	Float_t L2_AverageLayer_2PiZero_Et = 0;
	Float_t L1_AverageLayer_3PiZero_Et = 0;
	Float_t L2_AverageLayer_3PiZero_Et = 0;

	// Create string variables holding reconstructed energy definitions
	string Reco1_DefString = "(1x1, 5x3, 5x3, 3x3, 3x3)";
	string Reco2_DefString = "(3x3, 7x3, 5x3, 3x3, 3x3)";
	string Reco3_DefString = "(1x2, 5x2, 5x2, 1x2, 3x3)";
	string Reco3A_DefString = "(1x3, 5x3, 5x3, 1x3, 3x3)";
	string Reco3B_DefString = "(1x3, 5x3, 5x3, 1x3, 3x3) (No Pre-processing)";
	string Reco4_DefString = "(1x2, 5x2, 5x2, 1x2, None)";
	string Reco5_DefString = "(1x2, 5x2, 5x2, None, None)";
	string Reco6_DefString = "(None, 5x2, 5x2, None, None)";
	string RecoAll_DefString = "(3x3, 13x3, 13x3, 3x3, 3x3)";

	// Background variables
    Int_t backentries = (Int_t)tback->GetEntries();
	Float_t L0_AverageCell_Back_Et[3][3] = { 0 };
	Float_t L0_AverageFlipCell_Back_Et[3][3] = { 0 };
	Float_t L1_AverageCell_Back_Et[13][3] = { 0 };
	Float_t L1_AverageFlipCell_Back_Et[13][3] = { 0 };
	Float_t L2_AverageCell_Back_Et[13][3] = {0};
	Float_t L2_AverageFlipCell_Back_Et[13][3] = { 0 };
	Float_t L3_AverageCell_Back_Et[3][3] = { 0 };
	Float_t L3_AverageFlipCell_Back_Et[3][3] = { 0 };
	Float_t Had_AverageCell_Back_Et[3][3] = { 0 };
	Float_t Had_AverageFlipCell_Back_Et[3][3] = { 0 };
	Float_t Reco1_Back_Et[backentries]; // Reconstructed background energy using the following substructures (1x1, 5x3, 5x3, 3x3, 3x3)
	Float_t Reco2_Back_Et[backentries]; // Reconstructed background energy using the following substructures (3x3, 7x3, 5x3, 3x3, 3x3)
	Float_t Reco3_Back_Et[backentries]; // Reconstructed background energy using the following substructures (1x2, 5x2, 5x2, 1x2, 3x3). Done post-pre-processing
	Float_t Reco3A_Back_Et[backentries]; // Reconstructed background energy using the following substructures (1x3, 5x3, 5x3, 1x3, 3x3). Done post-pre-processing
	Float_t Reco3B_Back_Et[backentries]; // Reconstructed background energy using the following substructures (1x3, 5x3, 5x3, 1x3, 3x3). Done without pre-processing
	Float_t Reco4_Back_Et[backentries]; // Reconstructed background energy using the following substructures (1x2, 5x2, 5x2, 1x2, None). Done post-pre-processing
	Float_t Reco5_Back_Et[backentries]; // Reconstructed background energy using the following substructures (1x2, 5x2, 5x2, None, None). Done post-pre-processing
	Float_t Reco6_Back_Et[backentries]; // Reconstructed background energy using the following substructures (None, 5x2, 5x2, None, None). Done post-pre-processing
	Float_t RecoAll_Back_Et[backentries];

	// Number of events above Et threshold
	const Int_t netcuts = 150;
	Int_t N_Reco1_Back_EtCut[netcuts] = { 0 };
	Int_t N_Reco1_Sig_EtCut[netcuts] = { 0 };
	Int_t N_Reco2_Back_EtCut[netcuts] = { 0 };
	Int_t N_Reco2_Sig_EtCut[netcuts] = { 0 };
	Int_t N_Reco3_Back_EtCut[netcuts] = { 0 };
	Int_t N_Reco3_Sig_EtCut[netcuts] = { 0 };
	Int_t N_Reco3A_Back_EtCut[netcuts] = { 0 };
	Int_t N_Reco3A_Sig_EtCut[netcuts] = { 0 };
	Int_t N_Reco3B_Back_EtCut[netcuts] = { 0 };
	Int_t N_Reco3B_Sig_EtCut[netcuts] = { 0 };
	Int_t N_Reco4_Back_EtCut[netcuts] = { 0 };
	Int_t N_Reco4_Sig_EtCut[netcuts] = { 0 };
	Int_t N_Reco5_Back_EtCut[netcuts] = { 0 };
	Int_t N_Reco5_Sig_EtCut[netcuts] = { 0 };
	Int_t N_Reco6_Back_EtCut[netcuts] = { 0 };
	Int_t N_Reco6_Sig_EtCut[netcuts] = { 0 };
	Int_t N_RecoAll_Back_EtCut[netcuts] = { 0 };
	Int_t N_RecoAll_Sig_EtCut[netcuts] = { 0 };

	// Loop over all signal entries
	cout << "sigentries: " << sigentries << endl;
	for (Int_t i = 0; i < sigentries; i++) {
		tsig->GetEntry(i);

		Reco1_Sig_Et[i] = 0;
		Reco2_Sig_Et[i] = 0;
		Reco3_Sig_Et[i] = 0;
		Reco3A_Sig_Et[i] = 0;
		Reco3B_Sig_Et[i] = 0;
		Reco4_Sig_Et[i] = 0;
		Reco5_Sig_Et[i] = 0;
		Reco6_Sig_Et[i] = 0;
		RecoAll_Sig_Et[i] = 0;

		//TrueEt[i] = mctau->Et() / 1000.;
		//if (TrueEt[i] < 20.) continue;

		// Set flag for preprocessing if necessary by considering sum of off-center phi cells for all layers
		Int_t SigFlip = 0;
		SigFlip = PhiFlip_Et(L0_Sig_Et, L1_Sig_Et, L2_Sig_Et, L3_Sig_Et, Had_Sig_Et);

		// Loop over all cells to calculate average and reconstructed energies
		// Performance gains can be realized by doing the conversion to Gev and division by backentries outside the i loop
		// j = eta, k = phi
		for (Int_t j = 0; j < 13; j++) {
			for (Int_t k = 0; k < 3; k++) {
				Float_t L0_Reco_Et_Holder;
				Float_t L1_Reco_Et_Holder;
				Float_t L2_Reco_Et_Holder;
				Float_t L3_Reco_Et_Holder;
				Float_t Had_Reco_Et_Holder;

				PopulateRecoEtHolders(j, k, SigFlip, L0_Sig_Et, L0_Reco_Et_Holder, L1_Sig_Et, L1_Reco_Et_Holder, L2_Sig_Et, L2_Reco_Et_Holder, 
					L3_Sig_Et, L3_Reco_Et_Holder, Had_Sig_Et, Had_Reco_Et_Holder);

				// Calculate per-cell averages for L0, L3, and Had layers i.e 3x3
				// Calculate per-cell average energy
				Add3x3AverageNoFlipCellContribution(j, k, L0_Sig_Et, L0_AverageCell_Sig_Et, L3_Sig_Et, L3_AverageCell_Sig_Et, Had_Sig_Et,
					Had_AverageCell_Sig_Et);
				// Calculate per-cell average energy when pre-processed to flip phi, loading peripheral energy to one side
				Add3x3AverageFlipCellContribution(j, k, SigFlip, L0_Sig_Et, L0_AverageFlipCell_Sig_Et, L3_Sig_Et, L3_AverageFlipCell_Sig_Et,
					Had_Sig_Et, Had_AverageFlipCell_Sig_Et);

				// Calculate per-cell averages for L1 and L2 i.e 13x3
				Add13x3AverageNoFlipCellContribution(j, k, L1_Sig_Et, L1_AverageCell_Sig_Et, L2_Sig_Et, L2_AverageCell_Sig_Et);
				// Calculate per-cell average energy when pre-processed to flip phi, loading peripheral energy to one side
				Add13x3AverageFlipCellContribution(j, k, SigFlip, L1_Sig_Et, L1_AverageFlipCell_Sig_Et, L2_Sig_Et, L2_AverageFlipCell_Sig_Et);

				// Calculate reconstructed energy
				// Reconstruction definition 1
				AddReco1Contribution(i, j, k, Reco1_Sig_Et, L0_Reco_Et_Holder, L1_Reco_Et_Holder, L2_Reco_Et_Holder,
					L3_Reco_Et_Holder, Had_Reco_Et_Holder);
				// Reconstruction definition 2
				AddReco2Contribution(i, j, k, Reco2_Sig_Et, L0_Reco_Et_Holder, L1_Reco_Et_Holder, L2_Reco_Et_Holder,
					L3_Reco_Et_Holder, Had_Reco_Et_Holder);
				// Reconstruction definition 3
				AddReco3Contribution(i, j, k, Reco3_Sig_Et, L0_Reco_Et_Holder, L1_Reco_Et_Holder, L2_Reco_Et_Holder,
					L3_Reco_Et_Holder, Had_Reco_Et_Holder);
				// Reconstruction definition 3A
				AddReco3AContribution(i, j, k, Reco3A_Sig_Et, L0_Reco_Et_Holder, L1_Reco_Et_Holder, L2_Reco_Et_Holder,
					L3_Reco_Et_Holder, Had_Reco_Et_Holder);
				// Reconstruction definition 3B, done manually as the holder variables are pre-processed
				if (j < 3) {
					Reco3B_Sig_Et[i] += Had_Sig_Et[j][k];
				}
				if (j == 1) {
					Reco3B_Sig_Et[i] += L0_Sig_Et[j][k];
					Reco3B_Sig_Et[i] += L3_Sig_Et[j][k];
				}
				if (j > 3 && j < 9) {
					Reco3B_Sig_Et[i] += L1_Sig_Et[j][k];
					Reco3B_Sig_Et[i] += L2_Sig_Et[j][k];
				}
				// Reconstruction definition 4
				AddReco4Contribution(i, j, k, Reco4_Sig_Et, L0_Reco_Et_Holder, L1_Reco_Et_Holder, L2_Reco_Et_Holder,
					L3_Reco_Et_Holder, Had_Reco_Et_Holder);
				// Reconstruction definition 5
				AddReco5Contribution(i, j, k, Reco5_Sig_Et, L0_Reco_Et_Holder, L1_Reco_Et_Holder, L2_Reco_Et_Holder,
					L3_Reco_Et_Holder, Had_Reco_Et_Holder);
				// Reconstruction definition 6
				AddReco6Contribution(i, j, k, Reco6_Sig_Et, L0_Reco_Et_Holder, L1_Reco_Et_Holder, L2_Reco_Et_Holder,
					L3_Reco_Et_Holder, Had_Reco_Et_Holder);
				// Reconstruction all definition
				AddRecoAllContribution(i, j, k, RecoAll_Sig_Et, L0_Reco_Et_Holder, L1_Reco_Et_Holder, L2_Reco_Et_Holder,
					L3_Reco_Et_Holder, Had_Reco_Et_Holder);

				if (numofpi0 == 0) {
					L1_AverageLayer_0PiZero_Et += L1_Sig_Et[j][k];
					L2_AverageLayer_0PiZero_Et += L2_Sig_Et[j][k];
				}
				if (numofpi0 == 1) {
					L1_AverageLayer_1PiZero_Et += L1_Sig_Et[j][k];
					L2_AverageLayer_1PiZero_Et += L2_Sig_Et[j][k];
				}
				if (numofpi0 == 2) {
					L1_AverageLayer_2PiZero_Et += L1_Sig_Et[j][k];
					L2_AverageLayer_2PiZero_Et += L2_Sig_Et[j][k];
				}
				if (numofpi0 == 3) {
					L1_AverageLayer_3PiZero_Et += L1_Sig_Et[j][k];
					L2_AverageLayer_3PiZero_Et += L2_Sig_Et[j][k];
				}
			}
		}

		// Convert MeV to GeV
		Reco1_Sig_Et[i] /= 1000.;
		Reco2_Sig_Et[i] /= 1000.;
		Reco3_Sig_Et[i] /= 1000.;
		Reco3A_Sig_Et[i] /= 1000.;
		Reco3B_Sig_Et[i] /= 1000.;
		Reco4_Sig_Et[i] /= 1000.;
		Reco5_Sig_Et[i] /= 1000.;
		Reco6_Sig_Et[i] /= 1000.;
		RecoAll_Sig_Et[i] /= 1000.;

		cout << "Reco3_Sig_Et: " << Reco3_Sig_Et[i] << endl;

		// For each event, increment appropriate Et counters
		IncrementEtCounters(i, netcuts, Reco1_Sig_Et, N_Reco1_Sig_EtCut);
		IncrementEtCounters(i, netcuts, Reco2_Sig_Et, N_Reco2_Sig_EtCut);
		IncrementEtCounters(i, netcuts, Reco3_Sig_Et, N_Reco3_Sig_EtCut);
		IncrementEtCounters(i, netcuts, Reco3A_Sig_Et, N_Reco3A_Sig_EtCut);
		IncrementEtCounters(i, netcuts, Reco3B_Sig_Et, N_Reco3B_Sig_EtCut);
		IncrementEtCounters(i, netcuts, Reco4_Sig_Et, N_Reco4_Sig_EtCut);
		IncrementEtCounters(i, netcuts, Reco5_Sig_Et, N_Reco5_Sig_EtCut);
		IncrementEtCounters(i, netcuts, Reco6_Sig_Et, N_Reco6_Sig_EtCut);
		IncrementEtCounters(i, netcuts, RecoAll_Sig_Et, N_RecoAll_Sig_EtCut);

		cout << "N_Reco1_Sig_EtCut[0]: " << N_Reco1_Sig_EtCut[0] << endl;

		if (numofpi0 == 0) {
			zeropi0num += 1;
		}
		if (numofpi0 == 1) {
			onepi0num += 1;
		}
		if (numofpi0 == 2) {
			twopi0num += 1;
		}
		if (numofpi0 == 3) {
			threepi0num += 1;
		}
	}

	// Convert MeV to GeV
	L1_AverageLayer_0PiZero_Et /= 1000.;
	L2_AverageLayer_0PiZero_Et /= 1000.;
	L1_AverageLayer_1PiZero_Et /= 1000.;
	L2_AverageLayer_1PiZero_Et /= 1000.;
	L1_AverageLayer_2PiZero_Et /= 1000.;
	L2_AverageLayer_2PiZero_Et /= 1000.;
	L1_AverageLayer_3PiZero_Et /= 1000.;
	L2_AverageLayer_3PiZero_Et /= 1000.;

	// Now loop over all background entries
    for (Int_t i=0; i<backentries; i++) {
        tback->GetEntry(i);  // This fills declared variables from ntuple

		// cout << "Had Back: " << Had_Back_Et[0][1] << endl;

		Reco1_Back_Et[i] = 0;
		Reco2_Back_Et[i] = 0;
		Reco3_Back_Et[i] = 0;
		Reco3A_Back_Et[i] = 0;
		Reco3B_Back_Et[i] = 0;
		Reco4_Back_Et[i] = 0;
		Reco5_Back_Et[i] = 0;
		Reco6_Back_Et[i] = 0;
		RecoAll_Back_Et[i] = 0;

		// Set flag for preprocessing if necessary by considering sum of off-center phi cells for all layers
		Int_t BackFlip = 0;
		BackFlip = PhiFlip_Et(L0_Back_Et, L1_Back_Et, L2_Back_Et, L3_Back_Et, Had_Back_Et);

		// Loop over all cells to calculate average  and reconstructed energies
		// j = eta, k = phi
		for (Int_t j = 0; j < 13; j++) {
			for (Int_t k = 0; k < 3; k++) {
				Float_t L0_Reco_Et_Holder;
				Float_t L1_Reco_Et_Holder;
				Float_t L2_Reco_Et_Holder;
				Float_t L3_Reco_Et_Holder;
				Float_t Had_Reco_Et_Holder;

				PopulateRecoEtHolders(j, k, BackFlip, L0_Back_Et, L0_Reco_Et_Holder, L1_Back_Et, L1_Reco_Et_Holder, L2_Back_Et, L2_Reco_Et_Holder,
					L3_Back_Et, L3_Reco_Et_Holder, Had_Back_Et, Had_Reco_Et_Holder);

				// Calculate for L0, L3, and Had layers i.e 3x3
				// Calculate per-cell average energy
				Add3x3AverageNoFlipCellContribution(j, k, L0_Back_Et, L0_AverageCell_Back_Et, L3_Back_Et, L3_AverageCell_Back_Et, Had_Back_Et,
					Had_AverageCell_Back_Et);
				// Calculate per-cell average energy when pre-processed to flip phi, loading peripheral energy to one side
				Add3x3AverageFlipCellContribution(j, k, BackFlip, L0_Back_Et, L0_AverageFlipCell_Back_Et, L3_Back_Et, L3_AverageFlipCell_Back_Et,
					Had_Back_Et, Had_AverageFlipCell_Back_Et);

				// Calculate per-cell averages for L1 and L2 i.e 13x3
				Add13x3AverageNoFlipCellContribution(j, k, L1_Back_Et, L1_AverageCell_Back_Et, L2_Back_Et, L2_AverageCell_Back_Et);
				// Calculate per-cell average energy when pre-processed to flip phi, loading peripheral energy to one side
				Add13x3AverageFlipCellContribution(j, k, BackFlip, L1_Back_Et, L1_AverageFlipCell_Back_Et, L2_Back_Et, L2_AverageFlipCell_Back_Et);

				// Calculate reconstructed energy
				// Reconstruction definition 1
				AddReco1Contribution(i, j, k, Reco1_Back_Et, L0_Reco_Et_Holder, L1_Reco_Et_Holder, L2_Reco_Et_Holder,
					L3_Reco_Et_Holder, Had_Reco_Et_Holder);
				// Reconstruction definition 2
				AddReco2Contribution(i, j, k, Reco2_Back_Et, L0_Reco_Et_Holder, L1_Reco_Et_Holder, L2_Reco_Et_Holder,
					L3_Reco_Et_Holder, Had_Reco_Et_Holder);
				// Reconstruction definition 3
				AddReco3Contribution(i, j, k, Reco3_Back_Et, L0_Reco_Et_Holder, L1_Reco_Et_Holder, L2_Reco_Et_Holder,
					L3_Reco_Et_Holder, Had_Reco_Et_Holder);
				// Reconstruction definition 3A
				AddReco3AContribution(i, j, k, Reco3A_Back_Et, L0_Reco_Et_Holder, L1_Reco_Et_Holder, L2_Reco_Et_Holder,
					L3_Reco_Et_Holder, Had_Reco_Et_Holder);
				// Reconstruction definition 3B, done manually as the holder variables are pre-processed
				if (j < 3) {
					Reco3B_Back_Et[i] += Had_Back_Et[j][k];
				}
				if (j == 1) {
					Reco3B_Back_Et[i] += L0_Back_Et[j][k];
					Reco3B_Back_Et[i] += L3_Back_Et[j][k];
				}
				if (j > 3 && j < 9) {
					Reco3B_Back_Et[i] += L1_Back_Et[j][k];
					Reco3B_Back_Et[i] += L2_Back_Et[j][k];
				}
				// Reconstruction definition 4
				AddReco4Contribution(i, j, k, Reco4_Back_Et, L0_Reco_Et_Holder, L1_Reco_Et_Holder, L2_Reco_Et_Holder,
					L3_Reco_Et_Holder, Had_Reco_Et_Holder);
				// Reconstruction definition 5
				AddReco5Contribution(i, j, k, Reco5_Back_Et, L0_Reco_Et_Holder, L1_Reco_Et_Holder, L2_Reco_Et_Holder,
					L3_Reco_Et_Holder, Had_Reco_Et_Holder);
				// Reconstruction definition 6
				AddReco6Contribution(i, j, k, Reco6_Back_Et, L0_Reco_Et_Holder, L1_Reco_Et_Holder, L2_Reco_Et_Holder,
					L3_Reco_Et_Holder, Had_Reco_Et_Holder);
				// Reconstruction all definition
				AddRecoAllContribution(i, j, k, RecoAll_Back_Et, L0_Reco_Et_Holder, L1_Reco_Et_Holder, L2_Reco_Et_Holder,
					L3_Reco_Et_Holder, Had_Reco_Et_Holder);
			}
		}

		// Convert reconstructed energy from MeV to GeV
		Reco1_Back_Et[i] /= 1000.;
		Reco2_Back_Et[i] /= 1000.;
		Reco3_Back_Et[i] /= 1000.;
		Reco3A_Back_Et[i] /= 1000.;
		Reco3B_Back_Et[i] /= 1000.;
		Reco4_Back_Et[i] /= 1000.;
		Reco5_Back_Et[i] /= 1000.;
		Reco6_Back_Et[i] /= 1000.;
		RecoAll_Back_Et[i] /= 1000.;

		// For each event, increment appropriate Et counters
		IncrementEtCounters(i, netcuts, Reco1_Back_Et, N_Reco1_Back_EtCut);
		IncrementEtCounters(i, netcuts, Reco2_Back_Et, N_Reco2_Back_EtCut);
		IncrementEtCounters(i, netcuts, Reco3_Back_Et, N_Reco3_Back_EtCut);
		IncrementEtCounters(i, netcuts, Reco3A_Back_Et, N_Reco3A_Back_EtCut);
		IncrementEtCounters(i, netcuts, Reco3B_Back_Et, N_Reco3B_Back_EtCut);
		IncrementEtCounters(i, netcuts, Reco4_Back_Et, N_Reco4_Back_EtCut);
		IncrementEtCounters(i, netcuts, Reco5_Back_Et, N_Reco5_Back_EtCut);
		IncrementEtCounters(i, netcuts, Reco6_Back_Et, N_Reco6_Back_EtCut);
		IncrementEtCounters(i, netcuts, RecoAll_Back_Et, N_RecoAll_Back_EtCut);
    }

	// Average signal and background cell energies and convert from MeV to GeV
	for (Int_t j = 0; j < 13; j++) {
		for (Int_t k = 0; k < 3; k++) {
			// Calculate for L0, L3, and Had layers i.e 3x3
			Average3x3CellEts(j, k, backentries, L0_AverageCell_Back_Et, L0_AverageFlipCell_Back_Et, L3_AverageCell_Back_Et,
				L3_AverageFlipCell_Back_Et, Had_AverageCell_Back_Et, Had_AverageFlipCell_Back_Et);
			Average3x3CellEts(j, k, sigentries, L0_AverageCell_Sig_Et, L0_AverageFlipCell_Sig_Et, L3_AverageCell_Sig_Et,
				L3_AverageFlipCell_Sig_Et, Had_AverageCell_Sig_Et, Had_AverageFlipCell_Sig_Et);
			Average13x3CellEts(j, k, backentries, L1_AverageCell_Back_Et, L1_AverageFlipCell_Back_Et, L2_AverageCell_Back_Et,
				L2_AverageFlipCell_Back_Et);
			Average13x3CellEts(j, k, sigentries, L1_AverageCell_Sig_Et, L1_AverageFlipCell_Sig_Et, L2_AverageCell_Sig_Et,
				L2_AverageFlipCell_Sig_Et);
		}
	}

	// Test various windows to determine smallest that reasonably captures energy
	// Construct sub-energies for 3x3 layers e.g. 0, 3, and Had
	Float_t L0_Back_AverageLayerEt = 0, L0_Back_Average1x1Et = 0, L0_Back_Average1x2Et = 0, L0_Back_Average1x3Et = 0;
	Float_t L3_Back_AverageLayerEt = 0, L3_Back_Average1x1Et = 0, L3_Back_Average1x2Et = 0, L3_Back_Average1x3Et = 0;
	Float_t Had_Back_AverageLayerEt = 0, Had_Back_Average1x1Et = 0, Had_Back_Average1x2Et = 0, Had_Back_Average1x3Et = 0;
	Float_t L0_Sig_AverageLayerEt = 0, L0_Sig_Average1x1Et = 0, L0_Sig_Average1x2Et = 0, L0_Sig_Average1x3Et = 0;
	Float_t L3_Sig_AverageLayerEt = 0, L3_Sig_Average1x1Et = 0, L3_Sig_Average1x2Et = 0, L3_Sig_Average1x3Et = 0;
	Float_t Had_Sig_AverageLayerEt = 0, Had_Sig_Average1x1Et = 0, Had_Sig_Average1x2Et = 0, Had_Sig_Average1x3Et = 0;
	// j = eta, k = phi
	for (Int_t j = 0; j < 3; j++) {
		for (Int_t k = 0; k < 3; k++) {
			L0_Back_AverageLayerEt += L0_AverageFlipCell_Back_Et[j][k];
			L3_Back_AverageLayerEt += L3_AverageFlipCell_Back_Et[j][k];
			Had_Back_AverageLayerEt += Had_AverageFlipCell_Back_Et[j][k];
			L0_Sig_AverageLayerEt += L0_AverageFlipCell_Sig_Et[j][k];
			L3_Sig_AverageLayerEt += L3_AverageFlipCell_Sig_Et[j][k];
			Had_Sig_AverageLayerEt += Had_AverageFlipCell_Sig_Et[j][k];
			if (j == 1) {
				L0_Back_Average1x3Et += L0_AverageFlipCell_Back_Et[j][k];
				L3_Back_Average1x3Et += L3_AverageFlipCell_Back_Et[j][k];
				Had_Back_Average1x3Et += Had_AverageFlipCell_Back_Et[j][k];
				L0_Sig_Average1x3Et += L0_AverageFlipCell_Sig_Et[j][k];
				L3_Sig_Average1x3Et += L3_AverageFlipCell_Sig_Et[j][k];
				Had_Sig_Average1x3Et += Had_AverageFlipCell_Sig_Et[j][k];
			}
			if (j == 1 && k == 1) {
				L0_Back_Average1x1Et += L0_AverageFlipCell_Back_Et[j][k];
				L3_Back_Average1x1Et += L3_AverageFlipCell_Back_Et[j][k];
				Had_Back_Average1x1Et += Had_AverageFlipCell_Back_Et[j][k];
				L0_Sig_Average1x1Et += L0_AverageFlipCell_Sig_Et[j][k];
				L3_Sig_Average1x1Et += L3_AverageFlipCell_Sig_Et[j][k];
				Had_Sig_Average1x1Et += Had_AverageFlipCell_Sig_Et[j][k];
			}
			if (j == 1 && k < 2) {
				L0_Back_Average1x2Et += L0_AverageFlipCell_Back_Et[j][k];
				L3_Back_Average1x2Et += L3_AverageFlipCell_Back_Et[j][k];
				Had_Back_Average1x2Et += Had_AverageFlipCell_Back_Et[j][k];
				L0_Sig_Average1x2Et += L0_AverageFlipCell_Sig_Et[j][k];
				L3_Sig_Average1x2Et += L3_AverageFlipCell_Sig_Et[j][k];
				Had_Sig_Average1x2Et += Had_AverageFlipCell_Sig_Et[j][k];
			}
		}
	}

	// Construct sub-energies for 13x3 layers e.g. 1 and 2
	Float_t L1_Back_AverageLayerEt = 0, L1_Back_Average3x3Et = 0, L1_Back_Average5x3Et = 0, L1_Back_Average7x3Et = 0, L1_Back_Average9x3Et = 0, L1_Back_Average11x3Et = 0;
	Float_t L2_Back_AverageLayerEt = 0, L2_Back_Average3x3Et = 0, L2_Back_Average5x3Et = 0, L2_Back_Average7x3Et = 0, L2_Back_Average9x3Et = 0, L2_Back_Average11x3Et = 0;
	Float_t L1_Back_MainPeakEt = 0, L1_Back_OffPeakEt = 0;
	Float_t L2_Back_MainPeakEt = 0, L2_Back_OffPeakEt = 0;
	Float_t L1_Sig_AverageLayerEt = 0, L1_Sig_Average3x3Et = 0, L1_Sig_Average5x3Et = 0, L1_Sig_Average7x3Et = 0, L1_Sig_Average9x3Et = 0, L1_Sig_Average11x3Et = 0;
	Float_t L2_Sig_AverageLayerEt = 0, L2_Sig_Average3x3Et = 0, L2_Sig_Average5x3Et = 0, L2_Sig_Average7x3Et = 0, L2_Sig_Average9x3Et = 0, L2_Sig_Average11x3Et = 0;
	Float_t L1_Sig_MainPeakEt = 0, L1_Sig_OffPeakEt = 0;
	Float_t L2_Sig_MainPeakEt = 0, L2_Sig_OffPeakEt = 0;
	for (Int_t j = 0; j < 13; j++) {
		for (Int_t k = 0; k < 3; k++) {
			L1_Back_AverageLayerEt += L1_AverageCell_Back_Et[j][k];
			L2_Back_AverageLayerEt += L2_AverageCell_Back_Et[j][k];
			L1_Sig_AverageLayerEt += L1_AverageCell_Sig_Et[j][k];
			L2_Sig_AverageLayerEt += L2_AverageCell_Sig_Et[j][k];
			if (j > 0 && j < 12) {
				L1_Back_Average11x3Et += L1_AverageCell_Back_Et[j][k];
				L2_Back_Average11x3Et += L2_AverageCell_Back_Et[j][k];
				L1_Sig_Average11x3Et += L1_AverageCell_Sig_Et[j][k];
				L2_Sig_Average11x3Et += L2_AverageCell_Sig_Et[j][k];
			}
			if (j > 1 && j < 11) {
				L1_Back_Average9x3Et += L1_AverageCell_Back_Et[j][k];
				L2_Back_Average9x3Et += L2_AverageCell_Back_Et[j][k];
				L1_Sig_Average9x3Et += L1_AverageCell_Sig_Et[j][k];
				L2_Sig_Average9x3Et += L2_AverageCell_Sig_Et[j][k];
			}
			if (j > 2 && j < 10) {
				L1_Back_Average7x3Et += L1_AverageCell_Back_Et[j][k];
				L2_Back_Average7x3Et += L2_AverageCell_Back_Et[j][k];
				L1_Sig_Average7x3Et += L1_AverageCell_Sig_Et[j][k];
				L2_Sig_Average7x3Et += L2_AverageCell_Sig_Et[j][k];
			}
			if (j > 3 && j < 9) {
				L1_Back_Average5x3Et += L1_AverageCell_Back_Et[j][k];
				L2_Back_Average5x3Et += L2_AverageCell_Back_Et[j][k];
				L1_Sig_Average5x3Et += L1_AverageCell_Sig_Et[j][k];
				L2_Sig_Average5x3Et += L2_AverageCell_Sig_Et[j][k];
			}
			if (j > 4 && j < 8) {
				L1_Back_Average3x3Et += L1_AverageCell_Back_Et[j][k];
				L2_Back_Average3x3Et += L2_AverageCell_Back_Et[j][k];
				L1_Sig_Average3x3Et += L1_AverageCell_Sig_Et[j][k];
				L2_Sig_Average3x3Et += L2_AverageCell_Sig_Et[j][k];
			}
			if (k == 1 && j > 3 && j < 9) {
				L1_Back_MainPeakEt += L1_AverageFlipCell_Back_Et[j][k];
				L2_Back_MainPeakEt += L2_AverageFlipCell_Back_Et[j][k];
				L1_Sig_MainPeakEt += L1_AverageFlipCell_Sig_Et[j][k];
				L2_Sig_MainPeakEt += L2_AverageFlipCell_Sig_Et[j][k];
			}
			if (k == 0 && j > 4 && j < 8) {
				L1_Back_OffPeakEt += L1_AverageFlipCell_Back_Et[j][k];
				L2_Back_OffPeakEt += L2_AverageFlipCell_Back_Et[j][k];
				L1_Sig_OffPeakEt += L1_AverageFlipCell_Sig_Et[j][k];
				L2_Sig_OffPeakEt += L2_AverageFlipCell_Sig_Et[j][k];
			}
		}
	}

	ofstream textfile ("TauEtStudy.txt");

	textfile << "Reco 1 = " + Reco1_DefString << endl;
	textfile << "Reco 2 = " + Reco2_DefString << endl;
	textfile << "Reco 3 = " + Reco3_DefString << endl;
	textfile << "Reco 3A = " + Reco3A_DefString << endl;
	textfile << "Reco 3B = " + Reco3B_DefString << endl;
	textfile << "Reco 4 = " + Reco4_DefString << endl;
	textfile << "Reco 5 = " + Reco5_DefString << endl;
	textfile << "Reco 6 = " + Reco6_DefString << endl;
	textfile << "Reco All = " + RecoAll_DefString << endl << endl;

	textfile << "Background Reconstruction Energy Subset Candidates" << endl << endl;

	OutputText3x3EtSubsets(textfile, "L0", L0_Back_AverageLayerEt, L0_Back_Average1x1Et, L0_Back_Average1x2Et, L0_Back_Average1x3Et);
	OutputText13x3EtSubsets(textfile, "L1", L1_Back_AverageLayerEt, L1_Back_Average3x3Et, L1_Back_Average5x3Et, L1_Back_Average7x3Et,
		L1_Back_Average9x3Et, L1_Back_Average11x3Et);
	OutputText13x3EtSubsets(textfile, "L2", L2_Back_AverageLayerEt, L2_Back_Average3x3Et, L2_Back_Average5x3Et, L2_Back_Average7x3Et,
		L2_Back_Average9x3Et, L2_Back_Average11x3Et);
	OutputText3x3EtSubsets(textfile, "L3", L3_Back_AverageLayerEt, L3_Back_Average1x1Et, L3_Back_Average1x2Et, L3_Back_Average1x3Et);
	OutputText3x3EtSubsets(textfile, "Had", Had_Back_AverageLayerEt, Had_Back_Average1x1Et, Had_Back_Average1x2Et, Had_Back_Average1x3Et);

	textfile << "Signal Reconstruction Energy Subset Candidates" << endl << endl;

	OutputText3x3EtSubsets(textfile, "L0", L0_Sig_AverageLayerEt, L0_Sig_Average1x1Et, L0_Sig_Average1x2Et, L0_Sig_Average1x3Et);
	OutputText13x3EtSubsets(textfile, "L1", L1_Sig_AverageLayerEt, L1_Sig_Average3x3Et, L1_Sig_Average5x3Et, L1_Sig_Average7x3Et,
		L1_Sig_Average9x3Et, L1_Sig_Average11x3Et);
	OutputText13x3EtSubsets(textfile, "L2", L2_Sig_AverageLayerEt, L2_Sig_Average3x3Et, L2_Sig_Average5x3Et, L2_Sig_Average7x3Et,
		L2_Sig_Average9x3Et, L2_Sig_Average11x3Et);
	OutputText3x3EtSubsets(textfile, "L3", L3_Sig_AverageLayerEt, L3_Sig_Average1x1Et, L3_Sig_Average1x2Et, L3_Sig_Average1x3Et);
	OutputText3x3EtSubsets(textfile, "Had", Had_Sig_AverageLayerEt, Had_Sig_Average1x1Et, Had_Sig_Average1x2Et, Had_Sig_Average1x3Et);

	textfile << endl << endl;

	// Break 2-D array into 3 1-D arrays to feed into Graph2D
	Int_t ncolumns_13x3 = 39;
	Int_t ncolumns_3x3 = 9;
	Float_t Eta_13x3[ncolumns_13x3];
	Float_t Phi_13x3[ncolumns_13x3];
	Float_t Eta_3x3[ncolumns_3x3];
	Float_t Phi_3x3[ncolumns_3x3];
	Float_t L0_AverageCell_Back_Et_Z[ncolumns_3x3];
	Float_t L0_AverageFlipCell_Back_Et_Z[ncolumns_3x3];
	Float_t L1_AverageCell_Back_Et_Z[ncolumns_13x3];
	Float_t L1_AverageFlipCell_Back_Et_Z[ncolumns_13x3];
	Float_t L2_AverageCell_Back_Et_Z[ncolumns_13x3];
	Float_t L2_AverageFlipCell_Back_Et_Z[ncolumns_13x3];
	Float_t L3_AverageCell_Back_Et_Z[ncolumns_3x3];
	Float_t L3_AverageFlipCell_Back_Et_Z[ncolumns_3x3];
	Float_t Had_AverageCell_Back_Et_Z[ncolumns_3x3];
	Float_t Had_AverageFlipCell_Back_Et_Z[ncolumns_3x3];
	Float_t L0_AverageCell_Sig_Et_Z[ncolumns_3x3];
	Float_t L0_AverageFlipCell_Sig_Et_Z[ncolumns_3x3];
	Float_t L1_AverageCell_Sig_Et_Z[ncolumns_13x3];
	Float_t L1_AverageFlipCell_Sig_Et_Z[ncolumns_13x3];
	Float_t L2_AverageCell_Sig_Et_Z[ncolumns_13x3];
	Float_t L2_AverageFlipCell_Sig_Et_Z[ncolumns_13x3];
	Float_t L3_AverageCell_Sig_Et_Z[ncolumns_3x3];
	Float_t L3_AverageFlipCell_Sig_Et_Z[ncolumns_3x3];
	Float_t Had_AverageCell_Sig_Et_Z[ncolumns_3x3];
	Float_t Had_AverageFlipCell_Sig_Et_Z[ncolumns_3x3];
	for (Int_t i = 0; i < ncolumns_13x3; i++) {
		Int_t eta = float(i % 13);
		Int_t phi = float(int(i / 13));
		Eta_13x3[i] = eta;
		Phi_13x3[i] = phi;
		L1_AverageCell_Back_Et_Z[i] = L1_AverageCell_Back_Et[eta][phi];
		L1_AverageFlipCell_Back_Et_Z[i] = L1_AverageFlipCell_Back_Et[eta][phi];
		L2_AverageCell_Back_Et_Z[i] = L2_AverageCell_Back_Et[eta][phi];
		L2_AverageFlipCell_Back_Et_Z[i] = L2_AverageFlipCell_Back_Et[eta][phi];
		L1_AverageCell_Sig_Et_Z[i] = L1_AverageCell_Sig_Et[eta][phi];
		L1_AverageFlipCell_Sig_Et_Z[i] = L1_AverageFlipCell_Sig_Et[eta][phi];
		L2_AverageCell_Sig_Et_Z[i] = L2_AverageCell_Sig_Et[eta][phi];
		L2_AverageFlipCell_Sig_Et_Z[i] = L2_AverageFlipCell_Sig_Et[eta][phi];
	}

	for (Int_t i = 0; i < ncolumns_3x3; i++) {
		Int_t eta = float(i % 3);
		Int_t phi = float(int(i / 3));
		Eta_3x3[i] = eta;
		Phi_3x3[i] = phi;
		L0_AverageCell_Back_Et_Z[i] = L0_AverageCell_Back_Et[eta][phi];
		L0_AverageFlipCell_Back_Et_Z[i] = L0_AverageFlipCell_Back_Et[eta][phi];
		L3_AverageCell_Back_Et_Z[i] = L3_AverageCell_Back_Et[eta][phi];
		L3_AverageFlipCell_Back_Et_Z[i] = L3_AverageFlipCell_Back_Et[eta][phi];
		Had_AverageCell_Back_Et_Z[i] = Had_AverageCell_Back_Et[eta][phi];
		Had_AverageFlipCell_Back_Et_Z[i] = Had_AverageFlipCell_Back_Et[eta][phi];
		L0_AverageCell_Sig_Et_Z[i] = L0_AverageCell_Sig_Et[eta][phi];
		L0_AverageFlipCell_Sig_Et_Z[i] = L0_AverageFlipCell_Sig_Et[eta][phi];
		L3_AverageCell_Sig_Et_Z[i] = L3_AverageCell_Sig_Et[eta][phi];
		L3_AverageFlipCell_Sig_Et_Z[i] = L3_AverageFlipCell_Sig_Et[eta][phi];
		Had_AverageCell_Sig_Et_Z[i] = Had_AverageCell_Sig_Et[eta][phi];
		Had_AverageFlipCell_Sig_Et_Z[i] = Had_AverageFlipCell_Sig_Et[eta][phi];
	}

	// Create 1-D array of middle row energies
	const Int_t ncellsinrow = 13;
	Float_t L2_AverageCell_Back_Et_MiddleRow[ncellsinrow];
	Float_t L2_AverageCell_Sig_Et_MiddleRow[ncellsinrow];
	for (Int_t i = 0; i < ncellsinrow; i++) {
		L2_AverageCell_Back_Et_MiddleRow[i] = L2_AverageCell_Back_Et[i][1];
		L2_AverageCell_Sig_Et_MiddleRow[i] = L2_AverageCell_Sig_Et[i][1];
	}

	// Create array of Et cuts and signal/background ratios
	Float_t NetCutsX[netcuts];
	Float_t Reco1_SigBackEtRatio[netcuts];
	Float_t Reco1_SigEtCutEfficiency[netcuts];
	Float_t Reco1_BackEtCutEfficiency[netcuts];
	Float_t Reco2_SigBackEtRatio[netcuts];
	Float_t Reco2_SigEtCutEfficiency[netcuts];
	Float_t Reco2_BackEtCutEfficiency[netcuts];
	Float_t Reco3_SigBackEtRatio[netcuts];
	Float_t Reco3_SigEtCutEfficiency[netcuts];
	Float_t Reco3_BackEtCutEfficiency[netcuts];
	Float_t Reco3A_SigBackEtRatio[netcuts];
	Float_t Reco3A_SigEtCutEfficiency[netcuts];
	Float_t Reco3A_BackEtCutEfficiency[netcuts];
	Float_t Reco3B_SigBackEtRatio[netcuts];
	Float_t Reco3B_SigEtCutEfficiency[netcuts];
	Float_t Reco3B_BackEtCutEfficiency[netcuts];
	Float_t Reco4_SigBackEtRatio[netcuts];
	Float_t Reco4_SigEtCutEfficiency[netcuts];
	Float_t Reco4_BackEtCutEfficiency[netcuts];
	Float_t Reco5_SigBackEtRatio[netcuts];
	Float_t Reco5_SigEtCutEfficiency[netcuts];
	Float_t Reco5_BackEtCutEfficiency[netcuts];
	Float_t Reco6_SigBackEtRatio[netcuts];
	Float_t Reco6_SigEtCutEfficiency[netcuts];
	Float_t Reco6_BackEtCutEfficiency[netcuts];
	Float_t RecoAll_SigBackEtRatio[netcuts];
	Float_t RecoAll_SigEtCutEfficiency[netcuts];
	Float_t RecoAll_BackEtCutEfficiency[netcuts];
	for (Int_t i = 0; i < netcuts; i++) {
		NetCutsX[i] = i;
		Reco1_SigBackEtRatio[i] = RecoEtRatio(N_Reco1_Sig_EtCut[i], N_Reco1_Back_EtCut[i]);
		Reco2_SigBackEtRatio[i] = RecoEtRatio(N_Reco2_Sig_EtCut[i], N_Reco2_Back_EtCut[i]);
		Reco3_SigBackEtRatio[i] = RecoEtRatio(N_Reco3_Sig_EtCut[i], N_Reco3_Back_EtCut[i]);
		Reco4_SigBackEtRatio[i] = RecoEtRatio(N_Reco4_Sig_EtCut[i], N_Reco4_Back_EtCut[i]);
		Reco5_SigBackEtRatio[i] = RecoEtRatio(N_Reco5_Sig_EtCut[i], N_Reco5_Back_EtCut[i]);
		Reco6_SigBackEtRatio[i] = RecoEtRatio(N_Reco6_Sig_EtCut[i], N_Reco6_Back_EtCut[i]);
		RecoAll_SigBackEtRatio[i] = RecoEtRatio(N_RecoAll_Sig_EtCut[i], N_RecoAll_Back_EtCut[i]);

		RecoCutEfficiencies(Reco1_SigEtCutEfficiency[i], Reco1_BackEtCutEfficiency[i], N_Reco1_Sig_EtCut[i], N_Reco1_Sig_EtCut[0],
			N_Reco1_Back_EtCut[i], N_Reco1_Back_EtCut[0]);
		RecoCutEfficiencies(Reco2_SigEtCutEfficiency[i], Reco2_BackEtCutEfficiency[i], N_Reco2_Sig_EtCut[i], N_Reco2_Sig_EtCut[0],
			N_Reco2_Back_EtCut[i], N_Reco2_Back_EtCut[0]);
		RecoCutEfficiencies(Reco3_SigEtCutEfficiency[i], Reco3_BackEtCutEfficiency[i], N_Reco3_Sig_EtCut[i], N_Reco3_Sig_EtCut[0],
			N_Reco3_Back_EtCut[i], N_Reco3_Back_EtCut[0]);
		RecoCutEfficiencies(Reco3A_SigEtCutEfficiency[i], Reco3A_BackEtCutEfficiency[i], N_Reco3A_Sig_EtCut[i], N_Reco3A_Sig_EtCut[0],
			N_Reco3A_Back_EtCut[i], N_Reco3A_Back_EtCut[0]);
		RecoCutEfficiencies(Reco3B_SigEtCutEfficiency[i], Reco3B_BackEtCutEfficiency[i], N_Reco3B_Sig_EtCut[i], N_Reco3B_Sig_EtCut[0],
			N_Reco3B_Back_EtCut[i], N_Reco3B_Back_EtCut[0]);
		RecoCutEfficiencies(Reco4_SigEtCutEfficiency[i], Reco4_BackEtCutEfficiency[i], N_Reco4_Sig_EtCut[i], N_Reco4_Sig_EtCut[0],
			N_Reco4_Back_EtCut[i], N_Reco4_Back_EtCut[0]);
		RecoCutEfficiencies(Reco5_SigEtCutEfficiency[i], Reco5_BackEtCutEfficiency[i], N_Reco5_Sig_EtCut[i], N_Reco5_Sig_EtCut[0],
			N_Reco5_Back_EtCut[i], N_Reco5_Back_EtCut[0]);
		RecoCutEfficiencies(Reco6_SigEtCutEfficiency[i], Reco6_BackEtCutEfficiency[i], N_Reco6_Sig_EtCut[i], N_Reco6_Sig_EtCut[0],
			N_Reco6_Back_EtCut[i], N_Reco6_Back_EtCut[0]);
		RecoCutEfficiencies(RecoAll_SigEtCutEfficiency[i], RecoAll_BackEtCutEfficiency[i], N_RecoAll_Sig_EtCut[i], N_RecoAll_Sig_EtCut[0],
			N_RecoAll_Back_EtCut[i], N_RecoAll_Back_EtCut[0]);
		}

	// Comparison of peak energy to off-peak energy here
	textfile << "L2 Phi Peripheral Energy - Flipped and Unflipped" << endl << endl;
	textfile << "Background" << endl;
	textfile << "Unflipped [6][0]:" << endl;
	textfile << L2_AverageCell_Back_Et[6][0] << endl;
	textfile << "Unflipped [6][2]:" << endl;
	textfile << L2_AverageCell_Back_Et[6][2] << endl;
	textfile << "Flipped [6][0]:" << endl;
	textfile << L2_AverageFlipCell_Back_Et[6][0] << endl;
	textfile << "Flipped [6][2]:" << endl;
	textfile << L2_AverageFlipCell_Back_Et[6][2] << endl << endl;
	textfile << "Signal" << endl;
	textfile << "Unflipped [6][0]:" << endl;
	textfile << L2_AverageCell_Sig_Et[6][0] << endl;
	textfile << "Unflipped [6][2]:" << endl;
	textfile << L2_AverageCell_Sig_Et[6][2] << endl;
	textfile << "Flipped [6][0]:" << endl;
	textfile << L2_AverageFlipCell_Sig_Et[6][0] << endl;
	textfile << "Flipped [6][2]:" << endl;
	textfile << L2_AverageFlipCell_Sig_Et[6][2] << endl << endl << endl;


	textfile << "Off Peak (3 Cells) / Main Peak (5 Cells) Et (with Pre-processing)" << endl << endl;
	textfile << "Background" << endl;
	textfile << "L1 Off Peak Et: " << L1_Back_OffPeakEt << endl;
	textfile << "L1 Main Peak Et: " << L1_Back_MainPeakEt << endl;
	textfile << "L1 Off/Main Ratio: " << L1_Back_OffPeakEt / L1_Back_MainPeakEt << endl << endl;
	textfile << "L2 Off Peak Et: " << L2_Back_OffPeakEt << endl;
	textfile << "L2 Main Peak Et: " << L2_Back_MainPeakEt << endl;
	textfile << "L2 Off/Main Ratio: " << L2_Back_OffPeakEt / L2_Back_MainPeakEt << endl << endl;
	textfile << "Signal" << endl;
	textfile << "L1 Off Peak Et: " << L1_Sig_OffPeakEt << endl;
	textfile << "L1 Main Peak Et: " << L1_Sig_MainPeakEt << endl;
	textfile << "L1 Off/Main Ratio: " << L1_Sig_OffPeakEt / L1_Sig_MainPeakEt << endl << endl;
	textfile << "L2 Off Peak Et: " << L2_Sig_OffPeakEt << endl;
	textfile << "L2 Main Peak Et: " << L2_Sig_MainPeakEt << endl;
	textfile << "L2 Off/Main Ratio: " << L2_Sig_OffPeakEt / L2_Sig_MainPeakEt << endl << endl << endl;

	// Output number of events after various Et cuts
	OutputEventsAfterCutsInfo(textfile, "1", N_Reco1_Back_EtCut, N_Reco1_Sig_EtCut);
	OutputEventsAfterCutsInfo(textfile, "2", N_Reco2_Back_EtCut, N_Reco2_Sig_EtCut);
	OutputEventsAfterCutsInfo(textfile, "3", N_Reco3_Back_EtCut, N_Reco3_Sig_EtCut);
	OutputEventsAfterCutsInfo(textfile, "3A", N_Reco3A_Back_EtCut, N_Reco3A_Sig_EtCut);
	OutputEventsAfterCutsInfo(textfile, "3B", N_Reco3B_Back_EtCut, N_Reco3B_Sig_EtCut);
	OutputEventsAfterCutsInfo(textfile, "4", N_Reco4_Back_EtCut, N_Reco4_Sig_EtCut);
	OutputEventsAfterCutsInfo(textfile, "5", N_Reco5_Back_EtCut, N_Reco5_Sig_EtCut);
	OutputEventsAfterCutsInfo(textfile, "6", N_Reco6_Back_EtCut, N_Reco6_Sig_EtCut);
	OutputEventsAfterCutsInfo(textfile, "All", N_RecoAll_Back_EtCut, N_RecoAll_Sig_EtCut);

	// Output info on cuts that result in 90% signal efficiency
	textfile << "90% Signal Efficiency Cuts" << endl << endl;
	Output90PercentSignalInfo(textfile, "Reco 1 " + Reco1_DefString, netcuts, Reco1_SigEtCutEfficiency, Reco1_BackEtCutEfficiency);
	Output90PercentSignalInfo(textfile, "Reco 2 " + Reco2_DefString, netcuts, Reco2_SigEtCutEfficiency, Reco2_BackEtCutEfficiency);
	Output90PercentSignalInfo(textfile, "Reco 3 " + Reco3_DefString, netcuts, Reco3_SigEtCutEfficiency, Reco3_BackEtCutEfficiency);
	Output90PercentSignalInfo(textfile, "Reco 3A " + Reco3A_DefString, netcuts, Reco3A_SigEtCutEfficiency, Reco3A_BackEtCutEfficiency);
	Output90PercentSignalInfo(textfile, "Reco 3B " + Reco3B_DefString, netcuts, Reco3B_SigEtCutEfficiency, Reco3B_BackEtCutEfficiency);
	Output90PercentSignalInfo(textfile, "Reco 4 " + Reco4_DefString, netcuts, Reco4_SigEtCutEfficiency, Reco4_BackEtCutEfficiency);
	Output90PercentSignalInfo(textfile, "Reco 5 " + Reco5_DefString, netcuts, Reco5_SigEtCutEfficiency, Reco5_BackEtCutEfficiency);
	Output90PercentSignalInfo(textfile, "Reco 6 " + Reco6_DefString, netcuts, Reco6_SigEtCutEfficiency, Reco6_BackEtCutEfficiency);
	Output90PercentSignalInfo(textfile, "Reco All " + RecoAll_DefString, netcuts, RecoAll_SigEtCutEfficiency, RecoAll_BackEtCutEfficiency);
	textfile << endl << endl;

	textfile << "L1/L2 Ratio by Pi0 Number" << endl;
	textfile << "0 Pi0s" << endl;
	textfile << "Total Events: " << zeropi0num << endl;
	textfile << "Ratio: " << L1_AverageLayer_0PiZero_Et / L2_AverageLayer_0PiZero_Et << endl;
	textfile << "1 Pi0" << endl;
	textfile << "Total Events: " << onepi0num << endl;
	textfile << "Ratio: " << L1_AverageLayer_1PiZero_Et / L2_AverageLayer_1PiZero_Et << endl;
	textfile << "2 Pi0s" << endl;
	textfile << "Total Events: " << twopi0num << endl;
	textfile << "Ratio: " << L1_AverageLayer_2PiZero_Et / L2_AverageLayer_2PiZero_Et << endl;
	textfile << "3 Pi0s" << endl;
	textfile << "Total Events: " << threepi0num << endl;
	textfile << "Ratio: " << L1_AverageLayer_3PiZero_Et / L2_AverageLayer_3PiZero_Et << endl;

	textfile.close();

    // Now draw the plots.  Put each on a separate page.

	// 1-D plot of middle row
	// Background
	TGraph *grb0 = new TGraph(ncellsinrow, Eta_13x3, L2_AverageCell_Back_Et_MiddleRow);
	grb0->SetTitle("Background Average L2 Et (Middle Row)");
	grb0->GetXaxis()->SetTitle("Eta");
	grb0->GetYaxis()->SetTitle("Et (GeV)");
	// Signal
	TGraph *grs0 = new TGraph(ncellsinrow, Eta_13x3, L2_AverageCell_Sig_Et_MiddleRow);
	grs0->SetTitle("Signal Average L2 Et (Middle Row)");
	grs0->GetXaxis()->SetTitle("Eta");
	grs0->GetYaxis()->SetTitle("Et (GeV)");

	// 2-D plots of all cells
	// Background
	TGraph2D *grb1 = new TGraph2D(ncolumns_3x3, Eta_3x3, Phi_3x3, L0_AverageCell_Back_Et_Z);
	grb1->SetName("ave_l0_et");
	grb1->SetTitle("Background Average L0 Et");
	grb1->SetMinimum(0);
	grb1->SetMaximum(.3);
	grb1->GetXaxis()->SetTitle("Eta");
	grb1->GetYaxis()->SetTitle("Phi");
	grb1->GetZaxis()->SetTitle("Et (GeV)");
	TGraph2D *grb2 = new TGraph2D(ncolumns_3x3, Eta_3x3, Phi_3x3, L0_AverageFlipCell_Back_Et_Z);
	grb2->SetName("ave_l0_flip_et");
	grb2->SetTitle("Background Average L0 Flipped Et");
	grb2->SetMinimum(0);
	grb2->SetMaximum(.3);
	grb2->GetXaxis()->SetTitle("Eta");
	grb2->GetYaxis()->SetTitle("Phi");
	grb2->GetZaxis()->SetTitle("Et (GeV)");

	TGraph2D *grb3 = new TGraph2D(ncolumns_13x3, Eta_13x3, Phi_13x3, L1_AverageCell_Back_Et_Z);
	grb3->SetName("ave_l1_et");
	grb3->SetTitle("Background Average L1 Et");
	grb3->SetMinimum(0);
	grb3->SetMaximum(0.6);
	grb3->GetXaxis()->SetTitle("Eta");
	grb3->GetYaxis()->SetTitle("Phi");
	grb3->GetZaxis()->SetTitle("Et (GeV)");
	TGraph2D *grb4 = new TGraph2D(ncolumns_13x3, Eta_13x3, Phi_13x3, L1_AverageFlipCell_Back_Et_Z);
	grb4->SetName("ave_l1_flip_et");
	grb4->SetTitle("Background Average L1 Flipped Et");
	grb4->SetMinimum(0);
	grb4->SetMaximum(0.6);
	grb4->GetXaxis()->SetTitle("Eta");
	grb4->GetYaxis()->SetTitle("Phi");
	grb4->GetZaxis()->SetTitle("Et (GeV)");

	TGraph2D *grb5 = new TGraph2D(ncolumns_13x3, Eta_13x3, Phi_13x3, L2_AverageCell_Back_Et_Z);
	grb5->SetName("ave_l2_et");
	grb5->SetTitle("Background Average L2 Et");
	grb5->SetMinimum(0);
	grb5->SetMaximum(1.6);
	grb5->GetXaxis()->SetTitle("Eta");
	grb5->GetYaxis()->SetTitle("Phi");
	grb5->GetZaxis()->SetTitle("Et (GeV)");
	TGraph2D *grb6 = new TGraph2D(ncolumns_13x3, Eta_13x3, Phi_13x3, L2_AverageFlipCell_Back_Et_Z);
	grb6->SetName("ave_l2_flip_et");
	grb6->SetTitle("Background Average L2 Flipped Et");
	grb6->SetMinimum(0);
	grb6->SetMaximum(1.6);
	grb6->GetXaxis()->SetTitle("Eta");
	grb6->GetYaxis()->SetTitle("Phi");
	grb6->GetZaxis()->SetTitle("Et (GeV)");

	TGraph2D *grb7 = new TGraph2D(ncolumns_3x3, Eta_3x3, Phi_3x3, L3_AverageCell_Back_Et_Z);
	grb7->SetName("ave_3_et");
	grb7->SetTitle("Background Average L3 Et");
	grb7->SetMinimum(0);
	grb7->SetMaximum(0.2);
	grb7->GetXaxis()->SetTitle("Eta");
	grb7->GetYaxis()->SetTitle("Phi");
	grb7->GetZaxis()->SetTitle("Et (GeV)");
	TGraph2D *grb8 = new TGraph2D(ncolumns_3x3, Eta_3x3, Phi_3x3, L3_AverageFlipCell_Back_Et_Z);
	grb8->SetName("ave_l3_flip_et");
	grb8->SetTitle("Background Average L3 Flipped Et");
	grb8->SetMinimum(0);
	grb8->SetMaximum(0.2);
	grb8->GetXaxis()->SetTitle("Eta");
	grb8->GetYaxis()->SetTitle("Phi");
	grb8->GetZaxis()->SetTitle("Et (GeV)");

	TGraph2D *grb9 = new TGraph2D(ncolumns_3x3, Eta_3x3, Phi_3x3, Had_AverageCell_Back_Et_Z);
	grb9->SetName("ave_had_et");
	grb9->SetTitle("Background Average Had Et");
	grb9->SetMinimum(0);
	grb9->SetMaximum(.2);
	grb9->GetXaxis()->SetTitle("Eta");
	grb9->GetYaxis()->SetTitle("Phi");
	grb9->GetZaxis()->SetTitle("Et (GeV)");
	TGraph2D *grb10 = new TGraph2D(ncolumns_3x3, Eta_3x3, Phi_3x3, Had_AverageFlipCell_Back_Et_Z);
	grb10->SetName("ave_had_flip_et");
	grb10->SetTitle("Background Average Had Flipped Et");
	grb10->SetMinimum(0);
	grb10->SetMaximum(.2);
	grb10->GetXaxis()->SetTitle("Eta");
	grb10->GetYaxis()->SetTitle("Phi");
	grb10->GetZaxis()->SetTitle("Et (GeV)");

	// Signal
	TGraph2D *grs1 = new TGraph2D(ncolumns_3x3, Eta_3x3, Phi_3x3, L0_AverageCell_Sig_Et_Z);
	grs1->SetName("ave_l0_et");
	grs1->SetTitle("Signal Average L0 Et");
	grs1->SetMinimum(0);
	grs1->SetMaximum(2);
	grs1->GetXaxis()->SetTitle("Eta");
	grs1->GetYaxis()->SetTitle("Phi");
	grs1->GetZaxis()->SetTitle("Et (GeV)");
	TGraph2D *grs2 = new TGraph2D(ncolumns_3x3, Eta_3x3, Phi_3x3, L0_AverageFlipCell_Sig_Et_Z);
	grs2->SetName("ave_l0_flip_et");
	grs2->SetTitle("Signal Average L0 Flipped Et");
	grs2->SetMinimum(0);
	grs2->SetMaximum(2);
	grs2->GetXaxis()->SetTitle("Eta");
	grs2->GetYaxis()->SetTitle("Phi");
	grs2->GetZaxis()->SetTitle("Et (GeV)");

	TGraph2D *grs3 = new TGraph2D(ncolumns_13x3, Eta_13x3, Phi_13x3, L1_AverageCell_Sig_Et_Z);
	grs3->SetName("ave_l1_et");
	grs3->SetTitle("Signal Average L1 Et");
	grs3->SetMinimum(0);
	grs3->SetMaximum();
	grs3->GetXaxis()->SetTitle("Eta");
	grs3->GetYaxis()->SetTitle("Phi");
	grs3->GetZaxis()->SetTitle("Et (GeV)");
	TGraph2D *grs4 = new TGraph2D(ncolumns_13x3, Eta_13x3, Phi_13x3, L1_AverageFlipCell_Sig_Et_Z);
	grs4->SetName("ave_l1_flip_et");
	grs4->SetTitle("Signal Average L1 Flipped Et");
	grs4->SetMinimum(0);
	grs4->SetMaximum(3);
	grs4->GetXaxis()->SetTitle("Eta");
	grs4->GetYaxis()->SetTitle("Phi");
	grs4->GetZaxis()->SetTitle("Et (GeV)");

	TGraph2D *grs5 = new TGraph2D(ncolumns_13x3, Eta_13x3, Phi_13x3, L2_AverageCell_Sig_Et_Z);
	grs5->SetName("ave_l2_et");
	grs5->SetTitle("Signal Average L2 Et");
	grs5->SetMinimum(0);
	grs5->SetMaximum(8);
	grs5->GetXaxis()->SetTitle("Eta");
	grs5->GetYaxis()->SetTitle("Phi");
	grs5->GetZaxis()->SetTitle("Et (GeV)");
	TGraph2D *grs6 = new TGraph2D(ncolumns_13x3, Eta_13x3, Phi_13x3, L2_AverageFlipCell_Sig_Et_Z);
	grs6->SetName("ave_l2_flip_et");
	grs6->SetTitle("Signal Average L2 Flipped Et");
	grs6->SetMinimum(0);
	grs6->SetMaximum(8);
	grs6->GetXaxis()->SetTitle("Eta");
	grs6->GetYaxis()->SetTitle("Phi");
	grs6->GetZaxis()->SetTitle("Et (GeV)");

	TGraph2D *grs7 = new TGraph2D(ncolumns_3x3, Eta_3x3, Phi_3x3, L3_AverageCell_Sig_Et_Z);
	grs7->SetName("ave_3_et");
	grs7->SetTitle("Signal Average L3 Et");
	grs7->SetMinimum(0);
	grs7->SetMaximum(1);
	grs7->GetXaxis()->SetTitle("Eta");
	grs7->GetYaxis()->SetTitle("Phi");
	grs7->GetZaxis()->SetTitle("Et (GeV)");
	TGraph2D *grs8 = new TGraph2D(ncolumns_3x3, Eta_3x3, Phi_3x3, L3_AverageFlipCell_Sig_Et_Z);
	grs8->SetName("ave_l3_flip_et");
	grs8->SetTitle("Signal Average L3 Flipped Et");
	grs8->SetMinimum(0);
	grs8->SetMaximum(1);
	grs8->GetXaxis()->SetTitle("Eta");
	grs8->GetYaxis()->SetTitle("Phi");
	grs8->GetZaxis()->SetTitle("Et (GeV)");

	TGraph2D *grs9 = new TGraph2D(ncolumns_3x3, Eta_3x3, Phi_3x3, Had_AverageCell_Sig_Et_Z);
	grs9->SetName("ave_had_et");
	grs9->SetTitle("Signal Average Had Et");
	//grs9->SetMinimum(0);
	//grs9->SetMaximum(.2);
	grs9->GetXaxis()->SetTitle("Eta");
	grs9->GetYaxis()->SetTitle("Phi");
	grs9->GetZaxis()->SetTitle("Et (GeV)");
	TGraph2D *grs10 = new TGraph2D(ncolumns_3x3, Eta_3x3, Phi_3x3, Had_AverageFlipCell_Sig_Et_Z);
	grs10->SetName("ave_had_flip_et");
	grs10->SetTitle("Signal Average Had Flipped Et");
	//grs10->SetMinimum(0);
	//grs10->SetMaximum(.2);
	grs10->GetXaxis()->SetTitle("Eta");
	grs10->GetYaxis()->SetTitle("Phi");
	grs10->GetZaxis()->SetTitle("Et (GeV)");

	// Plot of signal/background ratio
	TGraph *gr11 = new TGraph(netcuts, NetCutsX, Reco3_SigBackEtRatio);
	gr11->SetName("sig_back_ratio");
	gr11->SetTitle("Signal/Background Ratio vs Et Cut");
	gr11->GetXaxis()->SetTitle("Et (GeV)");
	gr11->GetYaxis()->SetTitle("Signal/Background Ratio");

	// Reco 1 ROC curve
	TGraph *groc1 = new TGraph(netcuts, Reco1_BackEtCutEfficiency, Reco1_SigEtCutEfficiency);
	groc1->SetName("reco_1_roc");
	groc1->SetTitle("Reco 1 ROC");
	groc1->GetXaxis()->SetTitle("Background Efficiency");
	groc1->GetXaxis()->SetLimits(0, 0.2);
	groc1->GetYaxis()->SetTitle("Signal Efficiency");

	// Reco 2 ROC curve
	TGraph *groc2 = new TGraph(netcuts, Reco2_BackEtCutEfficiency, Reco2_SigEtCutEfficiency);
	groc2->SetName("reco_2_roc");
	groc2->SetTitle("Reco 2 ROC");
	groc2->GetXaxis()->SetTitle("Background Efficiency");
	//groc2->GetXaxis()->SetLimits(0, 0.2);
	groc2->GetYaxis()->SetTitle("Signal Efficiency");

	// Reco 3 ROC curve
	TGraph *groc3 = new TGraph(netcuts, Reco3_BackEtCutEfficiency, Reco3_SigEtCutEfficiency);
	groc3->SetName("reco_3_roc");
	groc3->SetTitle("Reco 3 ROC");
	groc3->GetXaxis()->SetTitle("Background Efficiency");
	//groc3->GetXaxis()->SetLimits(0, 0.2);
	groc3->GetYaxis()->SetTitle("Signal Efficiency");

	// Reco 3A ROC curve
	TGraph *groc3a = new TGraph(netcuts, Reco3A_BackEtCutEfficiency, Reco3A_SigEtCutEfficiency);
	groc3a->SetName("reco_3A_roc");
	groc3a->SetTitle("Reco 3A ROC");
	groc3a->GetXaxis()->SetTitle("Background Efficiency");
	//groc3a->GetXaxis()->SetLimits(0, 0.2);
	groc3a->GetYaxis()->SetTitle("Signal Efficiency");

	// Reco 3B ROC curve
	TGraph *groc3b = new TGraph(netcuts, Reco3B_BackEtCutEfficiency, Reco3B_SigEtCutEfficiency);
	groc3b->SetName("reco_3B_roc");
	groc3b->SetTitle("Reco 3B ROC");
	groc3b->GetXaxis()->SetTitle("Background Efficiency");
	//groc3b->GetXaxis()->SetLimits(0, 0.2);
	groc3b->GetYaxis()->SetTitle("Signal Efficiency");

	// Reco 4 ROC curve
	TGraph *groc4 = new TGraph(netcuts, Reco4_BackEtCutEfficiency, Reco4_SigEtCutEfficiency);
	groc4->SetName("reco_4_roc");
	groc4->SetTitle("Reco 4 ROC");
	groc4->GetXaxis()->SetTitle("Background Efficiency");
	//groc4->GetXaxis()->SetLimits(0, 0.2);
	groc4->GetYaxis()->SetTitle("Signal Efficiency");

	// Reco 5 ROC curve
	TGraph *groc5 = new TGraph(netcuts, Reco5_BackEtCutEfficiency, Reco5_SigEtCutEfficiency);
	groc5->SetName("reco_5_roc");
	groc5->SetTitle("Reco 5 ROC");
	groc5->GetXaxis()->SetTitle("Background Efficiency");
	groc5->GetXaxis()->SetLimits(0, 0.2);
	groc5->GetYaxis()->SetTitle("Signal Efficiency");

	// Reco 6 ROC curve
	TGraph *groc6 = new TGraph(netcuts, Reco6_BackEtCutEfficiency, Reco6_SigEtCutEfficiency);
	groc6->SetName("reco_6_roc");
	groc6->SetTitle("Reco 6 ROC");
	groc6->GetXaxis()->SetTitle("Background Efficiency");
	//groc6->GetXaxis()->SetLimits(0, 0.2);
	groc6->GetYaxis()->SetTitle("Signal Efficiency");

	// Reco All ROC curve
	TGraph *grocall = new TGraph(netcuts, RecoAll_BackEtCutEfficiency, RecoAll_SigEtCutEfficiency);
	grocall->SetName("reco_all_roc");
	grocall->SetTitle("Reco All ROC");
	grocall->GetXaxis()->SetTitle("Background Efficiency");
	//grocall->GetXaxis()->SetLimits(0, 0.2);
	grocall->GetYaxis()->SetTitle("Signal Efficiency");
	
	// Legend for ROC curve graph
	TLegend *rocleg = new TLegend(0.6, 0.1, 0.9, 0.4);
	rocleg->SetHeader("Reco Definition");

	// Histogram of reconstructed signal and background energies
	TH1F* h1reco1 = new TH1F("reco1sig", "Reco 1 Signal; Et (GeV); Entries", 100, 0., 150.);
	TH1F* h2reco1 = new TH1F("reco1back", "Reco 1 Background; Et (GeV); Entries", 100, 0., 150.);
	TH1F* h3reco1 = new TH1F("reco1 - true", "Reco 1 Signal - True Et; Et (Gev); Entries", 100, -100., 50.);
	TH1F* h1reco2 = new TH1F("reco2sig", "Reco 2 Signal; Et (GeV); Entries", 100, 0., 150.);
	TH1F* h2reco2 = new TH1F("reco2back", "Reco 2 Background; Et (GeV); Entries", 100, 0., 150.);
	TH1F* h3reco2 = new TH1F("reco2 - true", "Reco 2 Signal - True Et; Et (Gev); Entries", 100, -100., 50.);
	TH1F* h1reco3 = new TH1F("reco3sig", "Reco 3 Signal; Et (GeV); Entries", 100, 0., 150.);
	TH1F* h2reco3 = new TH1F("reco3back", "Reco 3 Background; Et (GeV); Entries", 100, 0., 150.);
	TH1F* h3reco3 = new TH1F("reco3 - true", "Reco 3 Signal - True Et; Et (Gev); Entries", 100, -100., 50.);
	TH1F* h1reco3a = new TH1F("reco3asig", "Reco 3 Signal; Et (GeV); Entries", 100, 0., 150.);
	TH1F* h2reco3a = new TH1F("reco3aback", "Reco 3 Background; Et (GeV); Entries", 100, 0., 150.);
	TH1F* h3reco3a = new TH1F("reco3a - true", "Reco 3 Signal - True Et; Et (Gev); Entries", 100, -100., 50.);
	TH1F* h1reco3b = new TH1F("reco3bsig", "Reco 3 Signal; Et (GeV); Entries", 100, 0., 150.);
	TH1F* h2reco3b = new TH1F("reco3bback", "Reco 3 Background; Et (GeV); Entries", 100, 0., 150.);
	TH1F* h3reco3b = new TH1F("reco3b - true", "Reco 3 Signal - True Et; Et (Gev); Entries", 100, -100., 50.);
	TH1F* h1reco4 = new TH1F("reco4sig", "Reco 4 Signal; Et (GeV); Entries", 100, 0., 150.);
	TH1F* h2reco4 = new TH1F("reco4back", "Reco 4 Background; Et (GeV); Entries", 100, 0., 150.);
	TH1F* h3reco4 = new TH1F("reco4 - true", "Reco 4 Signal - True Et; Et (Gev); Entries", 100, -100., 50.);
	TH1F* h1reco5 = new TH1F("reco5sig", "Reco 5 Signal; Et (GeV); Entries", 100, 0., 150.);
	TH1F* h2reco5 = new TH1F("reco5back", "Reco 5 Background; Et (GeV); Entries", 100, 0., 150.);
	TH1F* h3reco5 = new TH1F("reco5 - true", "Reco 5 Signal - True Et; Et (Gev); Entries", 100, -100., 50.);
	TH1F* h1reco6 = new TH1F("reco6sig", "Reco 6 Signal; Et (GeV); Entries", 100, 0., 150.);
	TH1F* h2reco6 = new TH1F("reco6back", "Reco 6 Background; Et (GeV); Entries", 100, 0., 150.);
	TH1F* h3reco6 = new TH1F("reco6 - true", "Reco 6 Signal - True Et; Et (Gev); Entries", 100, -100., 50.);
	TH1F* h1recoall = new TH1F("recoallsig", "Reco All Signal; Et (GeV); Entries", 100, 0., 150.);
	TH1F* h2recoall = new TH1F("recoallback", "Reco All Background; Et (GeV); Entries", 100, 0., 150.);
	TH1F* h3recoall = new TH1F("recoall - true", "Reco All Signal - True Et; Et (Gev); Entries", 100, -100., 50.);
	TH1F* htrue = new TH1F("true", "True Et; Et (GeV); Entries", 100, 0., 100.);

	// Fill histograms
	for (Int_t i = 0; i < sigentries; i++) {
		//if (TrueEt[i] > 20.) {
		htrue->Fill(TrueEt[i]);
		//}
	}

	FillRecoHistograms(h1reco1, h2reco1, h3reco1, sigentries, backentries, Reco1_Sig_Et, Reco1_Back_Et, TrueEt);
	FillRecoHistograms(h1reco2, h2reco2, h3reco2, sigentries, backentries, Reco2_Sig_Et, Reco2_Back_Et, TrueEt);
	FillRecoHistograms(h1reco3, h2reco3, h3reco3, sigentries, backentries, Reco3_Sig_Et, Reco3_Back_Et, TrueEt);
	FillRecoHistograms(h1reco3a, h2reco3a, h3reco3a, sigentries, backentries, Reco3A_Sig_Et, Reco3A_Back_Et, TrueEt);
	FillRecoHistograms(h1reco3b, h2reco3b, h3reco3b, sigentries, backentries, Reco3B_Sig_Et, Reco3B_Back_Et, TrueEt);
	FillRecoHistograms(h1reco4, h2reco4, h3reco4, sigentries, backentries, Reco4_Sig_Et, Reco4_Back_Et, TrueEt);
	FillRecoHistograms(h1reco5, h2reco5, h3reco5, sigentries, backentries, Reco5_Sig_Et, Reco5_Back_Et, TrueEt);
	FillRecoHistograms(h1reco6, h2reco6, h3reco6, sigentries, backentries, Reco6_Sig_Et, Reco6_Back_Et, TrueEt);
	FillRecoHistograms(h1recoall, h2recoall, h3recoall, sigentries, backentries, RecoAll_Sig_Et, RecoAll_Back_Et, TrueEt);

    // Define the Canvas
    TCanvas *c1 = new TCanvas("c1", "Graph Draw Options", 200, 10, 600, 400);

	// And write out to a file
	const char *pdfname = "TauEtStudy.pdf";
	string pdfstring = pdfname;

	

	grb0->Draw();
	c1->Print((pdfstring+"(").c_str());

	grb1->Draw("P0");
	c1->Print(pdfname);
	grb2->Draw("P0");
	c1->Print(pdfname);
	grb3->Draw("P0");
	c1->Print(pdfname);
	grb4->Draw("P0");
	c1->Print(pdfname);
	grb5->Draw("P0");
	c1->Print(pdfname);
	grb6->Draw("P0");
	c1->Print(pdfname);
	grb7->Draw("P0");
	c1->Print(pdfname);
	grb8->Draw("P0");
	c1->Print(pdfname);
	grb9->Draw("P0");
	c1->Print(pdfname);
	grb10->Draw("P0");
	c1->Print(pdfname);
	// Signal
	grs0->Draw();
	c1->Print(pdfname);
	grs1->Draw("P0");
	c1->Print(pdfname);
	grs2->Draw("P0");
	c1->Print(pdfname);
	grs3->Draw("P0");
	c1->Print(pdfname);
	grs4->Draw("P0");
	c1->Print(pdfname);
	grs5->Draw("P0");
	c1->Print(pdfname);
	grs6->Draw("P0");
	c1->Print(pdfname);
	grs7->Draw("P0");
	c1->Print(pdfname);
	grs8->Draw("P0");
	c1->Print(pdfname);
	grs9->Draw("P0");
	c1->Print(pdfname);
	grs10->Draw("P0");
	c1->Print(pdfname);
	gr11->Draw();
	c1->Print(pdfname);

	groc1->Draw("ALP");
	c1->Print(pdfname);
	groc2->Draw("ALP");
	c1->Print(pdfname);
	groc3->Draw("ALP");
	c1->Print(pdfname);
	groc4->Draw("ALP");
	c1->Print(pdfname);
	groc5->Draw("ALP");
	c1->Print(pdfname);
	groc6->Draw("ALP");
	c1->Print(pdfname);
	grocall->Draw("ALP");
	c1->Print(pdfname);

	groc1->Draw("ALP");
	groc1->SetTitle("Reco ROCs");
	groc2->Draw("same");
	groc2->SetLineColor(kRed);
	groc3->Draw("same");
	groc3->SetLineColor(kBlue);
	groc4->Draw("same");
	groc4->SetLineColor(kGreen);
	groc5->Draw("same");
	groc5->SetLineColor(kYellow);
	groc6->Draw("same");
	groc6->SetLineColor(kOrange);
	grocall->Draw("same");
	grocall->SetLineColor(kPink);
	rocleg->AddEntry(groc1, ("Reco 1 "+Reco1_DefString).c_str(), "l");
	rocleg->AddEntry(groc2, ("Reco 2 "+Reco2_DefString).c_str(), "l");
	rocleg->AddEntry(groc3, ("Reco 3 "+Reco3_DefString).c_str(), "l");
	rocleg->AddEntry(groc4, ("Reco 4 "+Reco4_DefString).c_str(), "l");
	rocleg->AddEntry(groc5, ("Reco 5 "+Reco5_DefString).c_str(), "l");
	rocleg->AddEntry(groc6, ("Reco 6 "+Reco6_DefString).c_str(), "l");
	rocleg->AddEntry(grocall, ("Reco All "+RecoAll_DefString).c_str(), "l");
	rocleg->Draw("same");
	c1->Print(pdfname);
	
	groc3->Draw("ALP");
	groc3a->Draw("same");
	groc3a->SetLineColor(kRed);
	groc3b->Draw("same");
	groc3b->SetLineColor(kBlue);
	rocleg->Clear();
	rocleg->AddEntry(groc3, ("Reco 3 " + Reco3_DefString).c_str(), "l");
	rocleg->AddEntry(groc3a, ("Reco 3A " + Reco3A_DefString).c_str(), "l");
	rocleg->AddEntry(groc3b, ("Reco 3B " + Reco3B_DefString).c_str(), "l");
	rocleg->Draw("same");
	c1->Print(pdfname);

	h1reco3->Draw();
	h1reco3->SetMaximum(240);
	h1reco3->SetLineColor(kRed);
	c1->Print(pdfname);

	h2reco3->Draw();
	h2reco3->SetLineColor(kBlue);
	c1->Print(pdfname);

	h1reco1->Draw();
	h1reco1->SetLineColor(kRed);
	h1reco1->SetTitle("Reco 1 Signal + Background");
	h2reco1->Draw("same");
	h2reco1->SetLineColor(kBlue);
	c1->Print(pdfname);

	h1reco2->Draw();
	h1reco2->SetLineColor(kRed);
	h1reco2->SetTitle("Reco 2 Signal + Background");
	h2reco2->Draw("same");
	h2reco2->SetLineColor(kBlue);
	c1->Print(pdfname);

	h1reco3->Draw();
	h1reco3->SetLineColor(kRed);
	h1reco3->SetTitle("Reco 3 Signal + Background");
	h2reco3->Draw("same");
	h2reco3->SetLineColor(kBlue);
	c1->Print(pdfname);

	h1reco4->Draw();
	h1reco4->SetLineColor(kRed);
	h1reco4->SetTitle("Reco 4 Signal + Background");
	h2reco4->Draw("same");
	h2reco4->SetLineColor(kBlue);
	c1->Print(pdfname);

	h1reco5->Draw();
	h1reco5->SetLineColor(kRed);
	h1reco5->SetTitle("Reco 5 Signal + Background");
	h2reco5->Draw("same");
	h2reco5->SetLineColor(kBlue);
	c1->Print(pdfname);

	h1reco6->Draw();
	h1reco6->SetLineColor(kRed);
	h1reco6->SetTitle("Reco 6 Signal + Background");
	h2reco6->Draw("same");
	h2reco6->SetLineColor(kBlue);
	c1->Print(pdfname);

	h1recoall->Draw();
	h1recoall->SetLineColor(kRed);
	h1recoall->SetTitle("Reco All Signal + Background");
	h2recoall->Draw("same");
	h2recoall->SetLineColor(kBlue);
	c1->Print(pdfname);

	h3reco1->Draw();
	c1->Print(pdfname);
	h3reco2->Draw();
	c1->Print(pdfname);
	h3reco3->Draw();
	c1->Print(pdfname);
	h3reco4->Draw();
	c1->Print(pdfname);
	h3reco5->Draw();
	c1->Print(pdfname);
	h3reco6->Draw();
	c1->Print(pdfname);
	h3recoall->Draw();
	c1->Print(pdfname);

	htrue->Draw();
	c1->Print((pdfstring + ")").c_str());
}