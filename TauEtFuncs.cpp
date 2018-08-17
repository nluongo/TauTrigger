#include <iostream>
using namespace std;

// Function receives each layer's Et array for a given entry and returns whether it should be flipped in phi. An entry is flipped if Et is weighted
// "above" (index 2) the center cell in phi as opposed to "below" (index 0) the center cell. This ensures we always orient the extra Et toward 0 in phi. 
int PhiFlip_Et(Float_t L0_Et_Array[3][3], Float_t L1_Et_Array[13][3], Float_t L2_Et_Array[13][3], Float_t L3_Et_Array[3][3], Float_t Had_Et_Array[3][3]) {
	if (L0_Et_Array[1][2] + L1_Et_Array[6][2] + L2_Et_Array[6][2] + L3_Et_Array[1][2] + Had_Et_Array[1][2] >
		L0_Et_Array[1][0] + L1_Et_Array[6][0] + L2_Et_Array[6][0] + L3_Et_Array[1][0] + Had_Et_Array[1][0]) {
		return 1;
	}
	else {
		return 0;
	}
}

// Function populate the Holder variables depending on whether the given entry has been flagged to be flipped in phi.
void PopulateRecoEtHolders(Int_t eta, Int_t phi, Int_t FlipEt, Float_t L0_Et[3][3], Float_t& L0_Holder, Float_t L1_Et[13][3], Float_t& L1_Holder, 
	Float_t L2_Et[13][3], Float_t& L2_Holder, Float_t L3_Et[3][3], Float_t& L3_Holder, Float_t Had_Et[3][3], Float_t& Had_Holder) {
	if (FlipEt == 0) {
		L0_Holder = L0_Et[eta][phi];
		L1_Holder = L1_Et[eta][phi];
		L2_Holder = L2_Et[eta][phi];
		L3_Holder = L3_Et[eta][phi];
		Had_Holder = Had_Et[eta][phi];
	}
	else {
		L0_Holder = L0_Et[eta][2 - phi];
		L1_Holder = L1_Et[eta][2 - phi];
		L2_Holder = L2_Et[eta][2 - phi];
		L3_Holder = L3_Et[eta][2 - phi];
		Had_Holder = Had_Et[eta][2 - phi];
	}
}

// Add contribution from 3x3 layers of current entry to overall average energy per cell, no preprocessing
void Add3x3AverageNoFlipCellContribution(Int_t eta, Int_t phi, Float_t L0_Et[3][3], Float_t L0_AverageCell_Et[3][3], Float_t L3_Et[3][3], 
	Float_t L3_AverageCell_Et[3][3], Float_t Had_Et[3][3], Float_t Had_AverageCell_Et[3][3]) {
	if (eta < 3 && phi < 3) {
		L0_AverageCell_Et[eta][phi] += L0_Et[eta][phi];
		L3_AverageCell_Et[eta][phi] += L3_Et[eta][phi];
		Had_AverageCell_Et[eta][phi] += Had_Et[eta][phi];
	}
}

// Add contribution from 3x3 layers of current entry to overall average energy per cell, with preprocessing
void Add3x3AverageFlipCellContribution(Int_t eta, Int_t phi, Int_t FlipEt, Float_t L0_Et[3][3], Float_t L0_AverageFlipCell_Et[3][3], Float_t L3_Et[3][3], 
	Float_t L3_AverageFlipCell_Et[3][3], Float_t Had_Et[3][3], Float_t Had_AverageFlipCell_Et[3][3]) {
	if (eta < 3 && phi < 3) {
		if (FlipEt == 0) {
			L0_AverageFlipCell_Et[eta][phi] += L0_Et[eta][phi];
			L3_AverageFlipCell_Et[eta][phi] += L3_Et[eta][phi];
			Had_AverageFlipCell_Et[eta][phi] += Had_Et[eta][phi];
		}
		else {
			L0_AverageFlipCell_Et[eta][phi] += L0_Et[eta][2 - phi];
			L3_AverageFlipCell_Et[eta][phi] += L3_Et[eta][2 - phi];
			Had_AverageFlipCell_Et[eta][phi] += Had_Et[eta][2 - phi];
		}
	}
}

// Add contribution from 13x3 layers of current entry to overall average energy per cell, no preprocessing
void Add13x3AverageNoFlipCellContribution(Int_t eta, Int_t phi, Float_t L1_Et[13][3], Float_t L1_AverageCell_Et[13][3], Float_t L2_Et[13][3],
	Float_t L2_AverageCell_Et[13][3]) {
	L1_AverageCell_Et[eta][phi] += L1_Et[eta][phi];
	L2_AverageCell_Et[eta][phi] += L2_Et[eta][phi];
}

// Add contribution from 13x3 layers of current entry to overall average energy per cell, with preprocessing
void Add13x3AverageFlipCellContribution(Int_t eta, Int_t phi, Int_t FlipEt, Float_t L1_Et[13][3], Float_t L1_AverageFlipCell_Et[13][3], 
	Float_t L2_Et[13][3], Float_t L2_AverageFlipCell_Et[13][3]) {
	if (FlipEt == 0) {
		L1_AverageFlipCell_Et[eta][phi] += L1_Et[eta][phi];
		L2_AverageFlipCell_Et[eta][phi] += L2_Et[eta][phi];
	}
	else {
		L1_AverageFlipCell_Et[eta][phi] += L1_Et[eta][2 - phi];
		L2_AverageFlipCell_Et[eta][phi] += L2_Et[eta][2 - phi];
	}
}

// Add energies of current cells to the reconstruction energy with definition 1 (1x1, 5x3, 5x3, 3x3, 3x3)
void AddReco1Contribution(Int_t entry, Int_t eta, Int_t phi, Int_t numberofentries, Float_t Reco1_Et[], Float_t L0_Reco_Et_Holder,
	Float_t L1_Reco_Et_Holder, Float_t L2_Reco_Et_Holder, Float_t L3_Reco_Et_Holder, Float_t Had_Reco_Et_Holder) {
	if (eta < 3 && phi < 3) {
		Reco1_Et[entry] += L3_Reco_Et_Holder;
		Reco1_Et[entry] += Had_Reco_Et_Holder;
	}
	if (eta == 1 && phi == 1) {
		Reco1_Et[entry] += L0_Reco_Et_Holder;
	}
	if (eta > 3 && eta < 9) {
		Reco1_Et[entry] += L1_Reco_Et_Holder;
		Reco1_Et[entry] += L2_Reco_Et_Holder;
	}
}

// Add energies of current cells to the reconstruction energy with definition 2 (3x3, 7x3, 5x3, 3x3, 3x3)
void AddReco2Contribution(Int_t entry, Int_t eta, Int_t phi, Int_t numberofentries, Float_t Reco2_Et[], Float_t L0_Reco_Et_Holder,
	Float_t L1_Reco_Et_Holder, Float_t L2_Reco_Et_Holder, Float_t L3_Reco_Et_Holder, Float_t Had_Reco_Et_Holder) {
	if (eta < 3) {
		Reco2_Et[entry] += L0_Reco_Et_Holder;
		Reco2_Et[entry] += L3_Reco_Et_Holder;
		Reco2_Et[entry] += Had_Reco_Et_Holder;
	}
	if (eta > 2 && eta < 10) {
		Reco2_Et[entry] += L1_Reco_Et_Holder;
	}
	if (eta > 3 && eta < 9) {
		Reco2_Et[entry] += L2_Reco_Et_Holder;
	}
}

// Add energies of current cells to the reconstruction energy with definition 3 (1x2, 5x2, 5x2, 1x2, 3x3)
void AddReco3Contribution(Int_t entry, Int_t eta, Int_t phi, Int_t numberofentries, Float_t Reco3_Et[], Float_t L0_Reco_Et_Holder, 
	Float_t L1_Reco_Et_Holder, Float_t L2_Reco_Et_Holder, Float_t L3_Reco_Et_Holder, Float_t Had_Reco_Et_Holder) {
	if (eta < 3 && phi < 3) {
		Reco3_Et[entry] += Had_Reco_Et_Holder;
	}
	if (eta == 1 && phi < 2) {
		Reco3_Et[entry] += L0_Reco_Et_Holder;
		Reco3_Et[entry] += L3_Reco_Et_Holder;
	}
	if (eta > 3 && eta < 9 && phi < 2) {
		Reco3_Et[entry] += L1_Reco_Et_Holder;
		Reco3_Et[entry] += L2_Reco_Et_Holder;
	}
}

// Add energies of current cells to the reconstruction energy with definition 3 (1x3, 5x3, 5x3, 1x3, 3x3)
void AddReco3AContribution(Int_t entry, Int_t eta, Int_t phi, Int_t numberofentries, Float_t Reco3A_Et[], Float_t L0_Reco_Et_Holder,
	Float_t L1_Reco_Et_Holder, Float_t L2_Reco_Et_Holder, Float_t L3_Reco_Et_Holder, Float_t Had_Reco_Et_Holder) {
	if (eta < 3) {
		Reco3A_Et[entry] += Had_Reco_Et_Holder;
	}
	if (eta == 1) {
		Reco3A_Et[entry] += L0_Reco_Et_Holder;
		Reco3A_Et[entry] += L3_Reco_Et_Holder;
	}
	if (eta > 3 && eta < 9) {
		Reco3A_Et[entry] += L1_Reco_Et_Holder;
		Reco3A_Et[entry] += L2_Reco_Et_Holder;
	}
}

// Add energies of current cells to the reconstruction energy with definition 4 (1x2, 5x2, 5x2, 1x2, None)
void AddReco4Contribution(Int_t entry, Int_t eta, Int_t phi, Int_t numberofentries, Float_t Reco4_Et[], Float_t L0_Reco_Et_Holder,
	Float_t L1_Reco_Et_Holder, Float_t L2_Reco_Et_Holder, Float_t L3_Reco_Et_Holder, Float_t Had_Reco_Et_Holder) {
	if (eta == 1 && phi < 2) {
		Reco4_Et[entry] += L0_Reco_Et_Holder;
		Reco4_Et[entry] += L3_Reco_Et_Holder;
	}
	if (eta > 3 && eta < 9 && phi < 2) {
		Reco4_Et[entry] += L1_Reco_Et_Holder;
		Reco4_Et[entry] += L2_Reco_Et_Holder;
	}
}

// Add energies of current cells to the reconstruction energy with definition 5 (1x2, 5x2, 5x2, None, None)
void AddReco5Contribution(Int_t entry, Int_t eta, Int_t phi, Int_t numberofentries, Float_t Reco5_Et[], Float_t L0_Reco_Et_Holder,
	Float_t L1_Reco_Et_Holder, Float_t L2_Reco_Et_Holder, Float_t L3_Reco_Et_Holder, Float_t Had_Reco_Et_Holder) {
	if (eta == 1 && phi < 2) {
		Reco5_Et[entry] += L0_Reco_Et_Holder;
	}
	if (eta > 3 && eta < 9 && phi < 2) {
		Reco5_Et[entry] += L1_Reco_Et_Holder;
		Reco5_Et[entry] += L2_Reco_Et_Holder;
	}
}

// Add energies of current cells to the reconstruction energy with definition 6 (None, 5x2, 5x2, None, None)
void AddReco6Contribution(Int_t entry, Int_t eta, Int_t phi, Int_t numberofentries, Float_t Reco6_Et[], Float_t L0_Reco_Et_Holder,
	Float_t L1_Reco_Et_Holder, Float_t L2_Reco_Et_Holder, Float_t L3_Reco_Et_Holder, Float_t Had_Reco_Et_Holder) {
	if (eta > 3 && eta < 9 && phi < 2) {
		Reco6_Et[entry] += L1_Reco_Et_Holder;
		Reco6_Et[entry] += L2_Reco_Et_Holder;
	}
}

// Add energies of current cells to the reconstruction energy with all avaiable cells (3x3, 13x3, 13x3, 3x3, 3x3)
void AddRecoAllContribution(Int_t entry, Int_t eta, Int_t phi, Int_t numberofentries, Float_t RecoAll_Et[], Float_t L0_Reco_Et_Holder,
	Float_t L1_Reco_Et_Holder, Float_t L2_Reco_Et_Holder, Float_t L3_Reco_Et_Holder, Float_t Had_Reco_Et_Holder) {
	if (eta < 3) {
		RecoAll_Et[entry] += Had_Reco_Et_Holder;
		RecoAll_Et[entry] += L0_Reco_Et_Holder;
		RecoAll_Et[entry] += L3_Reco_Et_Holder;
	}
	RecoAll_Et[entry] += L1_Reco_Et_Holder;
	RecoAll_Et[entry] += L2_Reco_Et_Holder;
}

// For each event, increment appropriate Et counters
void IncrementEtCounters(Int_t entry, Int_t netcuts, Float_t Reco_Et[], Int_t N_EtCut[]) {
	for (Int_t l = 0; l < netcuts; l++) {
		if (Reco_Et[entry] >= l)
			N_EtCut[l] += 1;
	}
}

// Average total Et cells for 3x3 layers, also convert MeV to Gev
void Average3x3CellEts(Int_t eta, Int_t phi, Int_t entries, Float_t L0_AverageCell_Et[3][3], Float_t L0_AverageFlipCell_Et[3][3], 
	Float_t L3_AverageCell_Et[3][3], Float_t L3_AverageFlipCell_Et[3][3], Float_t Had_AverageCell_Et[3][3], Float_t Had_AverageFlipCell_Et[3][3]) {
	if (eta < 3 && phi < 3) {
		L0_AverageCell_Et[eta][phi] = ((L0_AverageCell_Et[eta][phi] / 1000.) / entries);
		L3_AverageCell_Et[eta][phi] = ((L3_AverageCell_Et[eta][phi] / 1000.) / entries);
		Had_AverageCell_Et[eta][phi] = ((Had_AverageCell_Et[eta][phi] / 1000.) / entries);
		L0_AverageFlipCell_Et[eta][phi] = ((L0_AverageFlipCell_Et[eta][phi] / 1000.) / entries);
		L3_AverageFlipCell_Et[eta][phi] = ((L3_AverageFlipCell_Et[eta][phi] / 1000.) / entries);
		Had_AverageFlipCell_Et[eta][phi] = ((Had_AverageFlipCell_Et[eta][phi] / 1000.) / entries);
	}
}

// Average total Et cells for 3x3 layers, also convert MeV to GeV
void Average13x3CellEts(Int_t eta, Int_t phi, Int_t entries, Float_t L1_AverageCell_Et[13][3], Float_t L1_AverageFlipCell_Et[13][3], 
	Float_t L2_AverageCell_Et[13][3], Float_t L2_AverageFlipCell_Et[13][3]) {
	L1_AverageCell_Et[eta][phi] = ((L1_AverageCell_Et[eta][phi] / 1000.) / entries);
	L2_AverageCell_Et[eta][phi] = ((L2_AverageCell_Et[eta][phi] / 1000.) / entries);
	L1_AverageFlipCell_Et[eta][phi] = ((L1_AverageFlipCell_Et[eta][phi] / 1000.) / entries);
	L2_AverageFlipCell_Et[eta][phi] = ((L2_AverageFlipCell_Et[eta][phi] / 1000.) / entries);
}

float RecoEtRatio(Int_t N_Reco_Sig_EtCut, Int_t N_Reco_Back_EtCut) {
	if (N_Reco_Back_EtCut > 0) {
		return float(N_Reco_Sig_EtCut) / float(N_Reco_Back_EtCut);
	}
	else {
		return float(N_Reco_Sig_EtCut) / 1;
	}
}

float CutEfficiency(Int_t N_Entry_EtCut, Int_t N_Total_EtCut) {
	return float(N_Entry_EtCut) / float(N_Total_EtCut);
}

void RecoCutEfficiencies(Float_t& SigCutEfficiency, Float_t& BackCutEfficiency, Int_t N_Sig_Entry_EtCut, Int_t N_Sig_Total_EtCut,
	Int_t N_Back_Entry_EtCut, Int_t N_Back_Total_EtCut) {
	SigCutEfficiency = CutEfficiency(N_Sig_Entry_EtCut, N_Sig_Total_EtCut);
	BackCutEfficiency = CutEfficiency(N_Back_Entry_EtCut, N_Back_Total_EtCut);
}

void OutputRecoEventsAfterCuts(ofstream& writefile, string RecoName, Int_t N_Reco_EtCut[]) {
	writefile << RecoName << " Events After Et Cuts" << endl;
	for (Int_t i = 0; i < 10; i = i + 2) {
		writefile << "Et > " << i << " GeV: " << N_Reco_EtCut[i] << endl;
	}
	for (Int_t i = 10; i < 150; i = i + 10) {
		writefile << "Et > " << i << " GeV: " << N_Reco_EtCut[i] << endl;
	}
	writefile << endl;
}

void OutputSigBackRatioAfterCuts(ofstream& writefile, string RecoName, Int_t N_Reco_Sig_EtCut[], Int_t N_Reco_Back_EtCut[]) {
	writefile << RecoName << " Signal/Background Ratio After Et Cuts" << endl;
	for (Int_t i = 0; i < 10; i = i + 2) {
		writefile << "Et > " << i << " GeV: " << float(N_Reco_Sig_EtCut[i]) / float(N_Reco_Back_EtCut[i]) << endl;
	}
	for (Int_t i = 10; i < 150; i = i + 10) {
		writefile << "Et > " << i << " GeV: " << float(N_Reco_Sig_EtCut[i]) / float(N_Reco_Back_EtCut[i]) << endl;
	}
	writefile << endl;
}

void OutputEventsAfterCutsInfo(ofstream& writefile, string RecoName, Int_t N_Reco_Back_EtCut[], Int_t N_Reco_Sig_EtCut[]) {
	OutputRecoEventsAfterCuts(writefile, "Reco " + RecoName + " Background", N_Reco_Back_EtCut);
	OutputRecoEventsAfterCuts(writefile, "Reco " + RecoName + " Signal", N_Reco_Sig_EtCut);
	OutputSigBackRatioAfterCuts(writefile, "Reco " + RecoName, N_Reco_Sig_EtCut, N_Reco_Back_EtCut);
}

void FillRecoHistograms(TH1* sighist, TH1* backhist, TH1* sigminustruehist, Int_t sigentries, Int_t backentries, 
	Float_t Reco_Sig_Et[], Float_t Reco_Back_Et[], Float_t TrueEt[]) {
	for (Int_t i = 0; i < sigentries; i++) {
		sighist->Fill(Reco_Sig_Et[i]);
		sigminustruehist->Fill(Reco_Sig_Et[i] - TrueEt[i]);
	}
	for (Int_t i = 0; i < backentries; i++) {
		backhist->Fill(Reco_Back_Et[i]);
	}
}

void OutputText3x3EtSubsets(ofstream& writefile, string LayerName, Float_t AverageLayerEt, Float_t Average1x1Et, Float_t Average1x2Et,
	Float_t Average1x3Et) {
	writefile << LayerName << endl;
	writefile << "Total energy: " << AverageLayerEt << endl;
	writefile << "1x1 energy: " << Average1x1Et << "	Ratio: " << Average1x1Et / AverageLayerEt << endl;
	writefile << "1x2 energy: " << Average1x2Et << "	Ratio: " << Average1x2Et / AverageLayerEt << endl;
	writefile << "1x3 energy: " << Average1x3Et << "	Ratio: " << Average1x3Et / AverageLayerEt << endl << endl;
}

void OutputText13x3EtSubsets(ofstream& writefile, string LayerName, Float_t AverageLayerEt, Float_t Average3x3Et, Float_t Average5x3Et,
	Float_t Average7x3Et, Float_t Average9x3Et, Float_t Average11x3Et) {
	writefile << LayerName << endl;
	writefile << "Total energy: " << AverageLayerEt << endl;
	writefile << "3x3 energy: " << Average3x3Et << "		Ratio: " << Average3x3Et / AverageLayerEt << endl;
	writefile << "5x3 energy: " << Average5x3Et << "		Ratio: " << Average5x3Et / AverageLayerEt << endl;
	writefile << "7x3 energy: " << Average7x3Et << "		Ratio: " << Average7x3Et / AverageLayerEt << endl;
	writefile << "9x3 energy: " << Average9x3Et << "		Ratio: " << Average9x3Et / AverageLayerEt << endl;
	writefile << "11x3 energy: " << Average11x3Et << "		Ratio: " << Average11x3Et / AverageLayerEt << endl << endl;
}

void Output90PercentSignalInfo(ofstream& writefile, string RecoName, Int_t netcuts, Float_t Sig_EtCutEfficiency[], Float_t Back_EtCutEfficiency[]) {
	for (Int_t i = 0; i < netcuts; i++) {
		if (Sig_EtCutEfficiency[netcuts-i-1] > 0.9) {
			writefile << RecoName << endl;
			writefile << "Et Cut (GeV): " << netcuts - i - 1 << endl;
			writefile << "Signal Efficiency: " << Sig_EtCutEfficiency[netcuts - i - 1] << endl;
			writefile << "Background Efficiency: " << Back_EtCutEfficiency[netcuts - i - 1] << endl << endl;
			break;
		}
	}
}