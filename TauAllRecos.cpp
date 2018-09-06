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

	Int_t RecoDef[10];

	// i = L0 eta, j = L0 phi, k = L1 eta, l = L1 phi, m = L2 eta, n = L2 phi, p = L3 eta, q = L3 phi, r = Had eta, s = Had phi
	for (Int_t i = 1; i < 4; i=i+2) {
		for (Int_t j = 1; j < 4; j++) {
			for (Int_t k = 1; k < 14; k=k+2) {
				for (Int_t l = 1; l < 4; l++) {
					for (Int_t m = 1; m < 14; m=m+2) {
						for (Int_t n = 1; n < 4; n++) {
							for (Int_t p = 1; p < 4; p=p+2) {
								for (Int_t q = 1; q < 4; q++) {
									for (Int_t r = 1; r < 4; r=r+2) {
										for (Int_t s = 1; s < 4; s++) {
											if (i == 1 && j == 1 && k == 1 && l == 1 && m == 1 && n == 1 && p == 1 && q == 1 && r == 1 && s == 1) {
												SetRecoDef(RecoDef, i, j, k, l, m, n, p, q, r, s);
												GetRocAndEfficiencyFromRecoDef(tsig, tback, c1, RecoDef, textfile, "1", 0);
											}
											else if (i == 3 && j == 3 && k == 13 && l == 3 && m == 13 && n == 3 && p == 3 && q == 3 && r == 3 && s == 3) {
												SetRecoDef(RecoDef, i, j, k, l, m, n, p, q, r, s);
												GetRocAndEfficiencyFromRecoDef(tsig, tback, c1, RecoDef, textfile, "1", 2);
											}
											else {
												SetRecoDef(RecoDef, i, j, k, l, m, n, p, q, r, s);
												GetRocAndEfficiencyFromRecoDef(tsig, tback, c1, RecoDef, textfile, "1", 1);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	textfile.close();
}