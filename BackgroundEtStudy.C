#include <iostream>
using namespace std;

void BackgroundEtStudy() {
    // Open the file, hardcode the filename for now
    // May need to specify an absolute path if not running ROOT from command line
    TFile f("output_MB80.root");
    
    // Get a pointer to the TTree (name must match)
    TTree* t = (TTree*)f.Get("mytree");
    
    // Define some variables to read in
    // Arrays, dimensions must match size in ntuple
	Float_t L0Et[3][3];
	Float_t L1Et[13][3];
	Float_t L2Et[13][3];
	Float_t L3Et[3][3];
	Float_t HadEt[3][3];
    // Objects
    TLorentzVector* mctau = new TLorentzVector();
    
    // Make assignemnt between tree and declared variables
	t->SetBranchAddress("L0CellEt[3][3]", &L0Et[0][0]);  // Assigns to address of first array element
	t->SetBranchAddress("L1CellEt[13][3]", &L1Et[0][0]);  // Assigns to address of first array element
    t->SetBranchAddress("L2CellEt[13][3]", &L2Et[0][0]);  // Assigns to address of first array element
	t->SetBranchAddress("L3CellEt[3][3]", &L3Et[0][0]);  // Assigns to address of first array element
	t->SetBranchAddress("HadCellEt[3][3]", &HadEt[0][0]);  // Assigns to address of first array element

    // Now loop over all entries
    Int_t nentries = (Int_t)t->GetEntries();
	Float_t AverageL0Et[3][3] = { 0 };
	Float_t AveragePhiFlippedL0Et[3][3] = { 0 };
	Float_t AverageL1Et[13][3] = { 0 };
	Float_t AveragePhiFlippedL1Et[13][3] = { 0 };
	Float_t AverageL2Et[13][3] = {0};
	Float_t AveragePhiFlippedL2Et[13][3] = { 0 };
	Float_t AverageL3Et[3][3] = { 0 };
	Float_t AveragePhiFlippedL3Et[3][3] = { 0 };
	Float_t AverageHadEt[3][3] = { 0 };
	Float_t AveragePhiFlippedHadEt[3][3] = { 0 };
	Float_t Reco_1_Et[nentries];	// Reconstructed energy using the following substructures (1x1, 5x3, 5x3, 3x3, 3x3)
	Float_t Reco_2_Et[nentries];  // Reconstructed energy using the following substructures (3x3, 7x3, 5x3, 3x3, 3x3)
	Float_t Reco_All_Et[nentries]; // Reconstructed energy using all cells

	Int_t L0flip = 0, L1flip = 0, L2flip = 0, L3flip = 0, Hadflip = 0;

    for (Int_t i=0; i<nentries; i++) {
        t->GetEntry(i);  // This fills declared variables from ntuple
		
		Reco_1_Et[i] = 0;
		Reco_2_Et[i] = 0;
		Reco_All_Et[i] = 0;

		// Set flag for preprocessing if necessary
		L0flip = 0;
		if (L0Et[1][2] > L0Et[1][0]) {
			L0flip = 1;
		}

		L1flip = 0;
		if (L1Et[6][2] > L1Et[6][0]) {
			L1flip = 1;
		}

		L2flip = 0;
		if (L2Et[6][2] > L2Et[6][0]) {
			L2flip = 1;
		}

		L3flip = 0;
		if (L3Et[1][2] > L3Et[1][0]) {
			L3flip = 1;
		}

		Hadflip = 0;
		if (HadEt[1][2] > HadEt[1][0]) {
			Hadflip = 1;
		}
		
		// if (i % 1000 == 0) {
		// 	cout << "True Et: " << TrueEt[i] << endl;
		// 	cout << "1st Reconstructed Et (1x1, 5x3, 5x3, 3x3, 3x3): " << endl;
		// 	cout << Reco_1_Et[i] << endl;
		// 	cout << "2nd Reconstructed Et (3x3, 7x3, 5x3, 3x3, 3x3): " << endl;
		// 	cout << Reco_2_Et[i] << endl;
		// 	cout << "All Reconstructed Et: " << endl;
		// 	cout << Reco_All_Et[i] << endl << endl;
		// }

		// Loop over necessary indices for L0, L3, and Had i.e. 3x3
		for (Int_t j = 0; j < 3; j++) {
			for (Int_t k = 0; k < 3; k++) {
				// Calculate per-cell average energy
				AverageL0Et[j][k] += ((L0Et[j][k] / 1000.) / nentries);
				AverageL3Et[j][k] += ((L3Et[j][k] / 1000.) / nentries);
				AverageHadEt[j][k] += ((HadEt[j][k] / 1000.) / nentries);

				// Calculate reconstructed total energies
				if (j == 1 && k == 1)
					Reco_1_Et[i] += L0Et[j][k];
				Reco_1_Et[i] += L3Et[j][k];
				Reco_1_Et[i] += HadEt[j][k];
				Reco_2_Et[i] += L0Et[j][k];
				Reco_2_Et[i] += L3Et[j][k];
				Reco_2_Et[i] += HadEt[j][k];
				Reco_All_Et[i] += L0Et[j][k];
				Reco_All_Et[i] += L3Et[j][k];
				Reco_All_Et[i] += HadEt[j][k];

				// Calculate per-cell average energy when pre-processed to flip phi, loading peripheral energy to one side
				if (L0flip == 0) {
					AveragePhiFlippedL0Et[j][k] += ((L0Et[j][k] / 1000.) / nentries);
				}
				else {
					AveragePhiFlippedL0Et[j][2 - k] += ((L0Et[j][k] / 1000.) / nentries);
				}

				if (L3flip == 0) {
					AveragePhiFlippedL3Et[j][k] += ((L3Et[j][k] / 1000.) / nentries);
				}
				else {
					AveragePhiFlippedL3Et[j][2 - k] += ((L3Et[j][k] / 1000.) / nentries);
				}

				if (Hadflip == 0) {
					AveragePhiFlippedHadEt[j][k] += ((HadEt[j][k] / 1000.) / nentries);
				}
				else {
					AveragePhiFlippedHadEt[j][2 - k] += ((HadEt[j][k] / 1000.) / nentries);
				}
			}
		}

		//if (i % 1000 == 0) {
		//	cout << "True Et: " << TrueEt[i] << endl;
		//	cout << "1st Reconstructed Et (1x1, 5x3, 5x3, 3x3, 3x3): " << endl;
		//	cout << Reco_1_Et[i] << endl;
		//	cout << "2nd Reconstructed Et (3x3, 7x3, 5x3, 3x3, 3x3): " << endl;
		//	cout << Reco_2_Et[i] << endl;
		//	cout << "All Reconstructed Et: " << endl;
		//	cout << Reco_All_Et[i] << endl << endl;
		//}

		// Loop over necessary indices for L1 and L2 i.e. 13x3
		for (Int_t j = 0; j < 13; j++) {
			for (Int_t k = 0; k < 3; k++) {
				// Calculate per-cell average energy
				AverageL1Et[j][k] += ((L1Et[j][k] / 1000.) / nentries);
				AverageL2Et[j][k] += ((L2Et[j][k] / 1000.) / nentries);

				// Calculate reconstructed total energies
				if (j > 2 && j < 10) {
					Reco_2_Et[i] += L1Et[j][k];
					if (j > 3 && j < 9) {
						Reco_1_Et[i] += L1Et[j][k];
						Reco_1_Et[i] += L2Et[j][k];
						Reco_2_Et[i] += L2Et[j][k];
					}
				}
				Reco_All_Et[i] += L1Et[j][k];
				Reco_All_Et[i] += L2Et[j][k];

				// Calculate per-cell average energy when pre-processed to flip phi, loading peripheral energy to one side
				if (L1flip == 0) {
					AveragePhiFlippedL1Et[j][k] += ((L1Et[j][k] / 1000.) / nentries);
				}
				else {
					AveragePhiFlippedL1Et[j][2 - k] += ((L1Et[j][k] / 1000.) / nentries);
				}

				if (L2flip == 0) {
					AveragePhiFlippedL2Et[j][k] += ((L2Et[j][k] / 1000.) / nentries);
				}
				else {
					AveragePhiFlippedL2Et[j][2 - k] += ((L2Et[j][k] / 1000.) / nentries);
				}
			}
		}

		Reco_1_Et[i] /= 1000.;
		Reco_2_Et[i] /= 1000.;
		Reco_All_Et[i] /= 1000.;

		//if (i % 1000 == 0) {
		//		cout << "True Et: " << TrueEt[i] << endl;
		//		cout << "1st Reconstructed Et (1x1, 5x3, 5x3, 3x3, 3x3): " << endl;
		//		cout << Reco_1_Et[i] << endl;
		//		cout << "2nd Reconstructed Et (3x3, 7x3, 5x3, 3x3, 3x3): " << endl;
		//		cout << Reco_2_Et[i] << endl;
		//		cout << "All Reconstructed Et: " << endl;
		//		cout << Reco_All_Et[i] << endl << endl;
		//}
    }

	// cout << "Final value:" << endl;
    // cout << AverageL2Et[4][1] << endl;

	// Test various windows to determine smallest that reasonably captures energy
	// Construct sub-energies for 3x3 layers e.g. 0, 3, and Had
	Float_t L0_TotalEnergy = 0, L0_1x1Energy = 0, L0_3x1Energy = 0;
	Float_t L3_TotalEnergy = 0, L3_1x1Energy = 0, L3_3x1Energy = 0;
	Float_t Had_TotalEnergy = 0, Had_1x1Energy = 0, Had_3x1Energy = 0;
	for (Int_t j = 0; j < 3; j++) {
		for (Int_t k = 0; k < 3; k++) {
			L0_TotalEnergy += AverageL0Et[j][k];
			L3_TotalEnergy += AverageL3Et[j][k];
			Had_TotalEnergy += AverageHadEt[j][k];
			if (j == 1) {
				L0_3x1Energy += AverageL0Et[j][k];
				L3_3x1Energy += AverageL3Et[j][k];
				Had_3x1Energy += AverageHadEt[j][k];
			}
			if (j == 1 && k == 1) {
				L0_1x1Energy += AverageL0Et[j][k];
				L3_1x1Energy += AverageL3Et[j][k];
				Had_1x1Energy += AverageHadEt[j][k];
			}
		}
	}

	// Construct sub-energies for 13x3 layers e.g. 1 and 2
	Float_t L1_TotalEnergy = 0, L1_3x3Energy = 0, L1_5x3Energy = 0, L1_7x3Energy = 0, L1_9x3Energy = 0, L1_11x3Energy = 0;
	Float_t L2_TotalEnergy = 0, L2_3x3Energy = 0, L2_5x3Energy = 0, L2_7x3Energy = 0, L2_9x3Energy = 0, L2_11x3Energy = 0;
	for (Int_t j = 0; j < 13; j++) {
		for (Int_t k = 0; k < 3; k++) {
			L1_TotalEnergy += AverageL1Et[j][k];
			L2_TotalEnergy += AverageL2Et[j][k];
			if (j > 0 && j < 12) {
				L1_11x3Energy += AverageL1Et[j][k];
				L2_11x3Energy += AverageL2Et[j][k];
			}
			if (j > 1 && j < 11) {
				L1_9x3Energy += AverageL1Et[j][k];
				L2_9x3Energy += AverageL2Et[j][k];
			}
			if (j > 2 && j < 10) {
				L1_7x3Energy += AverageL1Et[j][k];
				L2_7x3Energy += AverageL2Et[j][k];
			}
			if (j > 3 && j < 9) {
				L1_5x3Energy += AverageL1Et[j][k];
				L2_5x3Energy += AverageL2Et[j][k];
			}
			if (j > 4 && j < 8) {
				L1_3x3Energy += AverageL1Et[j][k];
				L2_3x3Energy += AverageL2Et[j][k];
			}
		}
	}
	
	ofstream textfile;
	textfile.open("BackgroundEtStudy.txt");

	textfile << "Reconstruction Energy Subset Candidates" << endl << endl;
	textfile << "L0" << endl;
	textfile << "Total energy: " << L0_TotalEnergy << endl;
	textfile << "1x1 energy: " << L0_1x1Energy << "	Ratio: " << L0_1x1Energy / L0_TotalEnergy << endl;
	textfile << "3x1 energy: " << L0_3x1Energy << "	Ratio: " << L0_3x1Energy / L0_TotalEnergy << endl << endl;

	textfile << "L1" << endl;
	textfile << "Total energy: " << L1_TotalEnergy << endl;
	textfile << "3x3 energy: " << L1_3x3Energy << "		Ratio: " << L1_3x3Energy / L1_TotalEnergy << endl;
	textfile << "5x3 energy: " << L1_5x3Energy << "		Ratio: " << L1_5x3Energy / L1_TotalEnergy << endl;
	textfile << "7x3 energy: " << L1_7x3Energy << "		Ratio: " << L1_7x3Energy / L1_TotalEnergy << endl;
	textfile << "9x3 energy: " << L1_9x3Energy << "		Ratio: " << L1_9x3Energy / L1_TotalEnergy << endl;
	textfile << "11x3 energy: " << L1_11x3Energy << "		Ratio: " << L1_11x3Energy / L1_TotalEnergy << endl << endl;

	textfile << "L2" << endl;
	textfile << "Total energy: " << L2_TotalEnergy << endl;
	textfile << "3x3 energy: " << L2_3x3Energy << "		Ratio: " << L2_3x3Energy / L2_TotalEnergy << endl;
	textfile << "5x3 energy: " << L2_5x3Energy << "		Ratio: " << L2_5x3Energy / L2_TotalEnergy << endl;
	textfile << "7x3 energy: " << L2_7x3Energy << "		Ratio: " << L2_7x3Energy / L2_TotalEnergy << endl;
	textfile << "9x3 energy: " << L2_9x3Energy << "		Ratio: " << L2_9x3Energy / L2_TotalEnergy << endl;
	textfile << "11x3 energy: " << L2_11x3Energy << "		Ratio: " << L2_11x3Energy / L2_TotalEnergy << endl << endl;

	textfile << "L3" << endl;
	textfile << "Total energy: " << L3_TotalEnergy << endl;
	textfile << "1x1 energy: " << L3_1x1Energy << "	Ratio: " << L3_1x1Energy / L3_TotalEnergy << endl;
	textfile << "3x1 energy: " << L3_3x1Energy << "	Ratio: " << L3_3x1Energy / L3_TotalEnergy << endl << endl;

	textfile << "Had" << endl;
	textfile << "Total energy: " << Had_TotalEnergy << endl;
	textfile << "1x1 energy: " << Had_1x1Energy << "	Ratio: " << Had_1x1Energy / Had_TotalEnergy << endl;
	textfile << "3x1 energy: " << Had_3x1Energy << "	Ratio: " << Had_3x1Energy / Had_TotalEnergy << endl << endl;

	textfile << "Reco 1 = (1x1, 5x3, 5x3, 3x3, 3x3)" << endl;
	textfile << "Reco 2 = (3x3, 7x3, 5x3, 3x3, 3x3)" << endl;
	textfile << "Reco All = (3x3, 13x3, 13x3, 3x3, 3x3)" << endl << endl << endl;

	// Break 2-D array into 3 1-D array to feed into Graph2D
	Int_t ncolumns_13x3 = 39;
	Int_t ncolumns_3x3 = 9;
	Float_t Eta_13x3[ncolumns_13x3];
	Float_t Phi_13x3[ncolumns_13x3];
	Float_t Eta_3x3[ncolumns_3x3];
	Float_t Phi_3x3[ncolumns_3x3];
	Float_t AverageL0Et_Z[ncolumns_3x3];
	Float_t AveragePhiFlippedL0Et_Z[ncolumns_3x3];
	Float_t AverageL1Et_Z[ncolumns_13x3];
	Float_t AveragePhiFlippedL1Et_Z[ncolumns_13x3];
	Float_t AverageL2Et_Z[ncolumns_13x3];
	Float_t AveragePhiFlippedL2Et_Z[ncolumns_13x3];
	Float_t AverageL3Et_Z[ncolumns_3x3];
	Float_t AveragePhiFlippedL3Et_Z[ncolumns_3x3];
	Float_t AverageHadEt_Z[ncolumns_3x3];
	Float_t AveragePhiFlippedHadEt_Z[ncolumns_3x3];
	Int_t eta;
	Int_t phi;
	for (Int_t i = 0; i<ncolumns_13x3; i++) {
		eta = float(i % 13);
		phi = float(int(i / 13));
		Eta_13x3[i] = eta;
		Phi_13x3[i] = phi;
		AverageL1Et_Z[i] = AverageL1Et[eta][phi];
		AveragePhiFlippedL1Et_Z[i] = AveragePhiFlippedL1Et[eta][phi];
		AverageL2Et_Z[i] = AverageL2Et[eta][phi];
		AveragePhiFlippedL2Et_Z[i] = AveragePhiFlippedL2Et[eta][phi];
	}

	for (Int_t i = 0; i < ncolumns_3x3; i++) {
		eta = float(i % 3);
		phi = float(int(i / 3));
		Eta_3x3[i] = eta;
		Phi_3x3[i] = phi;
		AverageL0Et_Z[i] = AverageL0Et[eta][phi];
		AveragePhiFlippedL0Et_Z[i] = AveragePhiFlippedL0Et[eta][phi];
		AverageL3Et_Z[i] = AverageL3Et[eta][phi];
		AveragePhiFlippedL3Et_Z[i] = AveragePhiFlippedL3Et[eta][phi];
		AverageHadEt_Z[i] = AverageHadEt[eta][phi];
		AveragePhiFlippedHadEt_Z[i] = AveragePhiFlippedHadEt[eta][phi];
	}

	// Create 1-D array of middle row energies
	const Int_t ncellsinrow = 13;
	Float_t AverageL2Et_MiddleRow[ncellsinrow];
	for (Int_t i = 0; i < ncellsinrow; i++) {
		AverageL2Et_MiddleRow[i] = AverageL2Et[i][1];
	}

	textfile << "Phi Peripheral Energy - Flipped and Unflipped" << endl << endl;
	textfile << "Unflipped [6][0]:" << endl;
	textfile << AverageL2Et[6][0] << endl;
	textfile << "Unflipped [6][2]:" << endl;
	textfile << AverageL2Et[6][2] << endl;
	textfile << "Flipped [6][0]:" << endl;
	textfile << AveragePhiFlippedL2Et[6][0] << endl;
	textfile << "Flipped [6][2]:" << endl;
	textfile << AveragePhiFlippedL2Et[6][2] << endl;
	textfile.close();

    // Now draw the plots.  Put each on a separate page.

	// 1-D plot of middle row
	TGraph *gr0 = new TGraph(ncellsinrow, Eta_13x3, AverageL2Et_MiddleRow);
	gr0->SetTitle("Average L2 Et (Middle Row)");
	gr0->GetXaxis()->SetTitle("Eta");
	gr0->GetYaxis()->SetTitle("Et (GeV)");

	// 2-D plots of all cells
	TGraph2D *gr1 = new TGraph2D(ncolumns_3x3, Eta_3x3, Phi_3x3, AverageL0Et_Z);
	gr1->SetName("ave_l0_et");
	gr1->SetTitle("Average L0 Et");
	gr1->GetXaxis()->SetTitle("Eta");
	gr1->GetYaxis()->SetTitle("Phi");
	gr1->GetZaxis()->SetTitle("Et (GeV)");
	gr1->SetMinimum(0);
	gr1->SetMaximum(.3);
	TGraph2D *gr2 = new TGraph2D(ncolumns_3x3, Eta_3x3, Phi_3x3, AveragePhiFlippedL0Et_Z);
	gr2->SetName("ave_l0_flip_et");
	gr2->SetTitle("Average L0 Flipped Et");
	gr2->GetXaxis()->SetTitle("Eta");
	gr2->GetYaxis()->SetTitle("Phi");
	gr2->GetZaxis()->SetTitle("Et (GeV)");
	gr2->SetMinimum(0);
	gr2->SetMaximum(.3);

	TGraph2D *gr3 = new TGraph2D(ncolumns_13x3, Eta_13x3, Phi_13x3, AverageL1Et_Z);
	gr3->SetName("ave_l1_et");
	gr3->SetTitle("Average L1 Et");
	gr3->GetXaxis()->SetTitle("Eta");
	gr3->GetYaxis()->SetTitle("Phi");
	gr3->GetZaxis()->SetTitle("Et (GeV)");
	TGraph2D *gr4 = new TGraph2D(ncolumns_13x3, Eta_13x3, Phi_13x3, AveragePhiFlippedL1Et_Z);
	gr4->SetName("ave_l1_flip_et");
	gr4->SetTitle("Average L1 Flipped Et");
	gr4->GetXaxis()->SetTitle("Eta");
	gr4->GetYaxis()->SetTitle("Phi");
	gr4->GetZaxis()->SetTitle("Et (GeV)");

	TGraph2D *gr5 = new TGraph2D(ncolumns_13x3, Eta_13x3, Phi_13x3, AverageL2Et_Z);
	gr5->SetName("ave_l2_et");
	gr5->SetTitle("Average L2 Et");
	gr5->GetXaxis()->SetTitle("Eta");
	gr5->GetYaxis()->SetTitle("Phi");
	gr5->GetZaxis()->SetTitle("Et (GeV)");
	TGraph2D *gr6 = new TGraph2D(ncolumns_13x3, Eta_13x3, Phi_13x3, AveragePhiFlippedL2Et_Z);
	gr6->SetName("ave_l2_flip_et");
	gr6->SetTitle("Average L2 Flipped Et");
	gr6->GetXaxis()->SetTitle("Eta");
	gr6->GetYaxis()->SetTitle("Phi");
	gr6->GetZaxis()->SetTitle("Et (GeV)");

	TGraph2D *gr7 = new TGraph2D(ncolumns_3x3, Eta_3x3, Phi_3x3, AverageL3Et_Z);
	gr7->SetName("ave_3_et");
	gr7->SetTitle("Average L3 Et");
	gr7->GetXaxis()->SetTitle("Eta");
	gr7->GetYaxis()->SetTitle("Phi");
	gr7->GetZaxis()->SetTitle("Et (GeV)");
	gr7->SetMinimum(0);
	gr7->SetMaximum(0.2);
	TGraph2D *gr8 = new TGraph2D(ncolumns_3x3, Eta_3x3, Phi_3x3, AveragePhiFlippedL3Et_Z);
	gr8->SetName("ave_l3_flip_et");
	gr8->SetTitle("Average L3 Flipped Et");
	gr8->GetXaxis()->SetTitle("Eta");
	gr8->GetYaxis()->SetTitle("Phi");
	gr8->GetZaxis()->SetTitle("Et (GeV)");
	gr8->SetMinimum(0);
	gr8->SetMaximum(0.2);

	TGraph2D *gr9 = new TGraph2D(ncolumns_3x3, Eta_3x3, Phi_3x3, AverageHadEt_Z);
	gr9->SetName("ave_had_et");
	gr9->SetTitle("Average Had Et");
	gr9->GetXaxis()->SetTitle("Eta");
	gr9->GetYaxis()->SetTitle("Phi");
	gr9->GetZaxis()->SetTitle("Et (GeV)");
	gr9->SetMinimum(0);
	gr9->SetMaximum(.2);
	TGraph2D *gr10 = new TGraph2D(ncolumns_3x3, Eta_3x3, Phi_3x3, AveragePhiFlippedHadEt_Z);
	gr10->SetName("ave_had_flip_et");
	gr10->SetTitle("Average Had Flipped Et");
	gr10->GetXaxis()->SetTitle("Eta");
	gr10->GetYaxis()->SetTitle("Phi");
	gr10->GetZaxis()->SetTitle("Et (GeV)");
	gr10->SetMinimum(0);
	gr10->SetMaximum(.2);

    // Define the Canvas
    TCanvas *c1 = new TCanvas("c1", "Graph Draw Options", 200, 10, 600, 400);

	// And write out to a file
	gr0->Draw();
	c1->Print("BackgroundEtStudy.pdf(");

	gr1->Draw("P0");
	c1->Print("BackgroundEtStudy.pdf");

	gr2->Draw("P0");
	c1->Print("BackgroundEtStudy.pdf");

	gr3->Draw("P0");
	c1->Print("BackgroundEtStudy.pdf");

	gr4->Draw("P0");
	c1->Print("BackgroundEtStudy.pdf");

	gr5->Draw("P0");
	c1->Print("BackgroundEtStudy.pdf");

	gr6->Draw("P0");
	c1->Print("BackgroundEtStudy.pdf");

	gr7->Draw("P0");
	c1->Print("BackgroundEtStudy.pdf");

	gr8->Draw("P0");
	c1->Print("BackgroundEtStudy.pdf");

	gr9->Draw("P0");
	c1->Print("BackgroundEtStudy.pdf");

	gr10->Draw("P0");
	c1->Print("BackgroundEtStudy.pdf)");

}