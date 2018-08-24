using namespace std;
#include <iostream>

// This script runs over all given reconstruction energy definitions and outputs ROC curves and 90% efficiency stats for each definition
void TauAllRecos() {
	// Open file containing background, hardcode the filename for now
	TFile fback("output_MB80.root");
	// Open file containing signal
	TFile fsig("output_Z80.root");

	// Get a pointer to the background TTree (name must match)
	TTree* tback = (TTree*)fback.Get("mytree");
	// Get a pointer to the signal TTree
	TTree* tsig = (TTree*)fsig.Get("mytree");

	// Define the Canvas
	TCanvas *c1 = new TCanvas("c1", "Graph Draw Options", 200, 10, 600, 400);

	// And write out to a file
	ofstream textfile("TauAllRecos.txt");

	Int_t Reco1Def[10] = { 1, 1, 5, 3, 5, 3, 3, 3, 3, 3 };
	GetRocAndEfficiencyFromRecoDef(tsig, tback, c1, Reco1Def, textfile, "1", 0);

	Int_t Reco2Def[10] = { 3, 3, 7, 3, 5, 3, 3, 3, 3, 3 };
	GetRocAndEfficiencyFromRecoDef(tsig, tback, c1, Reco2Def, textfile, "2", 2);

	textfile.close();
}