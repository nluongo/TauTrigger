void TauFCore() {
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
	ofstream textfile("TauFCore.txt");

	Int_t CoreDef[2];
	Int_t IsolationDef[2];

	// j = isolation eta, k = isolation phi, m = core eta, n = core phi
	for (Int_t j = 1; j < 14; j = j + 2) {
		for (Int_t k = 1; k < 4; k++) {
			for (Int_t m = 1; m <= j; m = m + 2) {
				for (Int_t n = 1; n <= k; n++) {
					if (j == 1 && k == 1 && m == 1 && n == 1) {
						SetCoreAndIsolationDefs(m, n, j, k, CoreDef, IsolationDef);
						OutputFCoreTextAndHisto(tsig, tback, c1, textfile, CoreDef, IsolationDef, 0);
					}
					else if (j == 13 && k == 3 && m == 13 && n == 3) {
						SetCoreAndIsolationDefs(m, n, j, k, CoreDef, IsolationDef);
						OutputFCoreTextAndHisto(tsig, tback, c1, textfile, CoreDef, IsolationDef, 2);
					}
					else {
						SetCoreAndIsolationDefs(m, n, j, k, CoreDef, IsolationDef);
						OutputFCoreTextAndHisto(tsig, tback, c1, textfile, CoreDef, IsolationDef, 1);
					}
				}
			}
		}
	}
	textfile.close();
}