// Make signal Et histograms before and after a 95% FCore (5x3, 13x3) cut

void TauFCoreEtHisto() {
	// Open file containing signal
	TFile fsig("output_Z80.root");
	// Open file containing background
	TFile fback("output_MB80.root");

	// Get a pointer to the signal TTree
	TTree* tsig = (TTree*)fsig.Get("mytree");
	// Get a pointer to the background TTree
	TTree* tback = (TTree*)fback.Get("mytree");

	// Define the Canvas
	TCanvas *c1 = new TCanvas("c1", "Graph Draw Options", 200, 10, 600, 400);

	// Define histograms
	TH1F* fcorehisto = new TH1F("fcorehisto", "Signal FCore; FCore; Entries", 100, 0., 1.);

	// Arrays containing definitions of FCore core and isolation regions
	Int_t CoreDef[2];
	Int_t IsolationDef[2];

	// Define signal variables to read in
	Float_t L0_Sig_Et[13][3];
	Float_t L1_Sig_Et[13][3];
	Float_t L2_Sig_Et[13][3];
	Float_t L3_Sig_Et[13][3];
	Float_t Had_Sig_Et[13][3];
	TLorentzVector* mctau = new TLorentzVector();

	// Make assignment between tree and declared variables for signal
	tsig->SetBranchAddress("L0CellEt[3][3]", &L0_Sig_Et[0][0]);  // Assigns to address of first array element
	tsig->SetBranchAddress("L1CellEt[13][3]", &L1_Sig_Et[0][0]);  // Assigns to address of first array element
	tsig->SetBranchAddress("L2CellEt[13][3]", &L2_Sig_Et[0][0]);  // Assigns to address of first array element
	tsig->SetBranchAddress("L3CellEt[3][3]", &L3_Sig_Et[0][0]);  // Assigns to address of first array element
	tsig->SetBranchAddress("HadCellEt[3][3]", &Had_Sig_Et[0][0]);  // Assigns to address of first array element
	tsig->SetBranchAddress("mc_visibleTau", &mctau);


	Int_t sigentries = tsig->GetEntries();
	for (Int_t i = 0; i < sigentries; i++) {
		tsig->GetEntry(i);

		if (mctau->Et() < 20000.) continue;

		Float_t Et = 0;
		for (Int_t j = 0; j < 13; j++) {
			for (Int_t k = 0; k < 3; k++) {
				Et += L1_Sig_Et[j][k];
				Et += L2_Sig_Et[j][k];
			}
		}

		SetCoreAndIsolationDefs(5, 2, 13, 3, CoreDef, IsolationDef);
		Float_t FCore = CalculateFCore(L1_Sig_Et, L2_Sig_Et, CoreDef, IsolationDef);
		fcorehisto->Fill(FCore);
	}

	TH1F* backprehisto = new TH1F("backprehisto", "Background Et Before FCore Cut; Visible Et; Entries", 100, 0., 100.);
	TH1F* backposthisto = new TH1F("backposthisto", "Background Et After FCore Cut; Visible Et; Entries", 100, 0., 100.);
	TH1F* backratiohisto = new TH1F("backratiohisto", "Background Efficiency Due to FCore Cut; Et; Ratio", 100, 0., 100.);

	// Define signal variables to read in
	Float_t L0_Back_Et[13][3];
	Float_t L1_Back_Et[13][3];
	Float_t L2_Back_Et[13][3];
	Float_t L3_Back_Et[13][3];
	Float_t Had_Back_Et[13][3];
	 
	// Make assignment between tree and declared variables for background
	tback->SetBranchAddress("L0CellEt[3][3]", &L0_Back_Et[0][0]);		// Assigns to address of first array element
	tback->SetBranchAddress("L1CellEt[13][3]", &L1_Back_Et[0][0]);		// Assigns to address of first array element
	tback->SetBranchAddress("L2CellEt[13][3]", &L2_Back_Et[0][0]);		// Assigns to address of first array element
	tback->SetBranchAddress("L3CellEt[3][3]", &L3_Back_Et[0][0]);		// Assigns to address of first array element
	tback->SetBranchAddress("HadCellEt[3][3]", &Had_Back_Et[0][0]);	// Assigns to address of first array element

	// 95% of signal entries are placed in an FCore bin above this one
	Int_t Sig95PercentFCoreBin = Find95PercentBin(fcorehisto);

	// Loop over background entries and fill histograms before and after 95% FCore cut
	Int_t backentries = tback->GetEntries();
	for (Int_t i = 0; i < backentries; i++) {
		tback->GetEntry(i); 

		Float_t Et = 0;

		// Set flag for preprocessing if necessary by considering sum of off-center phi cells for all layers
		Int_t BackFlip = 0;
		BackFlip = PhiFlip_Et(L0_Back_Et, L1_Back_Et, L2_Back_Et, L3_Back_Et, Had_Back_Et);

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

				// Use baseline reco (3) energy definition to calculate observed energy
				AddReco3Contribution(j, k, Et, L0_Reco_Et_Holder, L1_Reco_Et_Holder, L2_Reco_Et_Holder,
					L3_Reco_Et_Holder, Had_Reco_Et_Holder);
			}
		}


		SetCoreAndIsolationDefs(5, 2, 13, 3, CoreDef, IsolationDef);
		Float_t FCore = CalculateFCore(L1_Back_Et, L2_Back_Et, CoreDef, IsolationDef);
		Int_t FCoreBin = fcorehisto->GetXaxis()->FindBin(FCore);

		backprehisto->Fill(Et / 1000.);

		if (FCoreBin > Sig95PercentFCoreBin) {
			backposthisto->Fill(Et/1000.);
		}
	}

	TH1F* sigprehisto = new TH1F("sigprehisto", "Signal Et Before FCore Cut; Visible Et; Entries", 100, 0., 100.);
	TH1F* sigposthisto = new TH1F("sigposthisto", "Signal Et After FCore Cut; Visible Et; Entries", 100, 0., 100.);
	TH1F* sigratiohisto = new TH1F("sigratiohisto", "Signal Efficiency Due to FCore Cut; Visible Et; Ratio", 100, 0., 100.);

	// Loop over signal entries and fill histograms before and after 95% FCore cut
	for (Int_t i = 0; i < sigentries; i++) {
		tsig->GetEntry(i);

		Float_t Et = 0;

		// Set flag for preprocessing if necessary by considering sum of off-center phi cells for all layers
		Int_t SigFlip = 0;
		SigFlip = PhiFlip_Et(L0_Sig_Et, L1_Sig_Et, L2_Sig_Et, L3_Sig_Et, Had_Sig_Et);

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

				// Use baseline reco (3) energy definition to calculate observed energy
				AddReco3Contribution(j, k, Et, L0_Reco_Et_Holder, L1_Reco_Et_Holder, L2_Reco_Et_Holder,
					L3_Reco_Et_Holder, Had_Reco_Et_Holder);
			}
		}

		SetCoreAndIsolationDefs(5, 2, 13, 3, CoreDef, IsolationDef);
		Float_t FCore = CalculateFCore(L1_Sig_Et, L2_Sig_Et, CoreDef, IsolationDef);
		Int_t FCoreBin = fcorehisto->GetXaxis()->FindBin(FCore);

		sigprehisto->Fill(Et / 1000.);

		if (FCoreBin > Sig95PercentFCoreBin) {
			sigposthisto->Fill(Et / 1000.);
		}
	}

	// Collision rate of detector in Hz
	const Float_t CollisionRate = 30000000;
	// Number of events in the given sample
	const Float_t EventsInSample = 5500;
	// Number of time that this sample represents
	Float_t SampleTime = EventsInSample / CollisionRate;
	
	Float_t VisibleEt[100] = { 0 };
	Float_t SingleTauRate[100] = { 0 };
	Float_t SingleTauRateWithFCoreCut[100] = { 0 };
	Float_t DiTauRate[100] = { 0 };
	Float_t DiTauRateWithFCoreCut[100] = { 0 };
	Int_t RunningBackTotal = backprehisto->GetEntries();
	Int_t RunningBackTotalWithFCoreCut = backposthisto->GetEntries();

	for (Int_t i = 0; i < 101; i++) {
		RunningBackTotal -= backprehisto->GetBinContent(i);
		RunningBackTotalWithFCoreCut -= backposthisto->GetBinContent(i);
		VisibleEt[i] = i;
		SingleTauRate[i] = RunningBackTotal / SampleTime;
		SingleTauRateWithFCoreCut[i] = RunningBackTotalWithFCoreCut / SampleTime;
		DiTauRate[i] = (SingleTauRate[i] * SingleTauRate[i]) / CollisionRate;
		DiTauRateWithFCoreCut[i] = (SingleTauRateWithFCoreCut[i] * SingleTauRateWithFCoreCut[i]) / CollisionRate;
		cout << "RunningBackTotal: " << RunningBackTotal << endl;
		cout << "RunningBackTotalWithFCoreCut: " << RunningBackTotalWithFCoreCut << endl;
		cout << "VisibleEt: " << VisibleEt[i] << endl;
		cout << "SingleTauRate: " << SingleTauRate[i] << endl;
		cout << "SingleTauRateWithFCoreCut: " << SingleTauRateWithFCoreCut[i] << endl;
		cout << "DiTauRate: " << DiTauRate[i] << endl;
		cout << "DiTauRateWithFCoreCut: " << DiTauRateWithFCoreCut[i] << endl;
	}

	TGraph *gr0 = new TGraph(100, VisibleEt, SingleTauRate);
	gr0->SetTitle("Single Tau Trigger Rate");
	gr0->GetXaxis()->SetTitle("Visible Et (GeV)");
	gr0->GetYaxis()->SetTitle("Single Tau Rate (Hz)");

	TGraph *gr1 = new TGraph(100, VisibleEt, SingleTauRateWithFCoreCut);
	gr1->SetTitle("Single Tau Trigger Rate");
	gr1->GetXaxis()->SetTitle("Visible Et (GeV)");
	gr1->GetYaxis()->SetTitle("Single Tau Rate (Hz)");

	TGraph *gr2 = new TGraph(100, VisibleEt, DiTauRate);
	gr2->SetTitle("Di-Tau Trigger Rate");
	gr2->GetXaxis()->SetTitle("Visible Et (GeV)");
	gr2->GetYaxis()->SetTitle("Di-Tau Rate (Hz)");

	TGraph *gr3 = new TGraph(100, VisibleEt, DiTauRateWithFCoreCut);
	gr3->SetTitle("Di-Tau Trigger Rate");
	gr3->GetXaxis()->SetTitle("Visible Et (GeV)");
	gr3->GetYaxis()->SetTitle("Di-Tau Rate (Hz)");

	//Legend for FCore histogram
	TLegend *fcoreleg = new TLegend(0.7, 0.1, 0.9, 0.2);

	//Legend for trigger rate
	TLegend *grleg = new TLegend(0.6, 0.7, 0.9, 0.9);

	fcorehisto->Draw("histo");
	TLine *line = new TLine(Sig95PercentFCoreBin * 0.01, 0, Sig95PercentFCoreBin * 0.01, fcorehisto->GetMaximum());
	fcoreleg->AddEntry(line, "95% FCore (True Et > 20 GeV)", "l");
	fcoreleg->Draw();
	line->Draw();
	c1->Print("TauFCoreEtHisto.pdf(");

	backprehisto->Draw("histo");
	backprehisto->SetMaximum(13000);
	c1->Print("TauFCoreEtHisto.pdf");
	backposthisto->Draw("histo");
	backposthisto->SetMaximum(13000);
	c1->Print("TauFCoreEtHisto.pdf");
	backratiohisto->Divide(backposthisto, backprehisto);
	backratiohisto->SetMarkerStyle(5);
	backratiohisto->Draw("p");
	c1->Print("TauFCoreEtHisto.pdf");

	sigprehisto->Draw("histo");
	sigprehisto->SetMaximum(170);
	c1->Print("TauFCoreEtHisto.pdf");
	sigposthisto->Draw("histo");
	sigposthisto->SetMaximum(170);
	c1->Print("TauFCoreEtHisto.pdf");
	sigratiohisto->Divide(sigposthisto, sigprehisto);
	sigratiohisto->SetMarkerStyle(5);
	sigratiohisto->Draw("p");
	sigratiohisto->SetStats(FALSE);
	c1->Print("TauFCoreEtHisto.pdf");
	
	gr1->Draw("AP");
	gr1->GetYaxis()->SetRangeUser(50000, 50000000);
	gr1->GetXaxis()->SetRangeUser(10, 50);
	gr1->SetMarkerStyle(3);
	gr1->SetMarkerColor(kGreen);
	gr0->Draw("P");
	gr0->SetMarkerStyle(2);
	gr0->GetXaxis()->SetRangeUser(10, 50);
	grleg->AddEntry(gr0, "No FCore Cut", "p");
	grleg->AddEntry(gr1, "With 95% FCore Cut", "p");
	grleg->Draw("same");
	c1->SetLogy();
	c1->Print("TauFCoreEtHisto.pdf");

	gr3->Draw("AP");
	gr3->GetYaxis()->SetRangeUser(100, 50000000);
	gr3->GetXaxis()->SetRangeUser(10, 50);
	gr3->SetMarkerStyle(3);
	gr3->SetMarkerColor(kGreen);
	gr2->Draw("P");
	gr2->SetMarkerStyle(2);
	gr2->GetXaxis()->SetRangeUser(10, 50);
	TLine *line2 = new TLine(10, 200000, 50, 200000);
	line2->Draw();
	grleg->Clear();
	grleg->AddEntry(gr2, "No FCore Cut", "p");
	grleg->AddEntry(gr3, "With 95% FCore Cut", "p");
	grleg->AddEntry(line2, "200 kHz Trigger Rate", "l");
	grleg->Draw("same");
	c1->SetLogy();
	c1->Print("TauFCoreEtHisto.pdf)");
}