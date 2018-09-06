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

// Given the eta definition, return bounds between which valid cells must fall
void GetEtaBounds(Int_t& LowerBound, Int_t& UpperBound, Int_t EtaDef, Int_t LayerNumber) {
	if (EtaDef == 0) {
		LowerBound = -1;
		UpperBound = -1;
		return;
	}
	Int_t EtaOffsetFromCenter = (EtaDef - 1) / 2;
	if (LayerNumber == 0 || LayerNumber == 3 || LayerNumber == 4) {
		LowerBound = 1 - EtaOffsetFromCenter;
		UpperBound = 1 + EtaOffsetFromCenter;
	}
	else if (LayerNumber == 1 || LayerNumber == 2) {
		LowerBound = 6 - EtaOffsetFromCenter;
		UpperBound = 6 + EtaOffsetFromCenter;
	}
}

// Given the phi definition, return bounds between which valid cells must fall
void GetPhiBounds(Int_t& LowerBound, Int_t& UpperBound, Int_t PhiDef) {
	if (PhiDef == 1) {
		LowerBound = 1;
		UpperBound = 1;
	}
	else if (PhiDef == 2) {
		LowerBound = 0;
		UpperBound = 1;
	}
	else if (PhiDef == 3) {
		LowerBound = 0;
		UpperBound = 2;
	}
}

// Given the relevant layer and eta/phi definitions for that layer, calculate eta/phi bounds and add all cells within those bounds to the overall Et
void AddLayerContribution(Float_t& Et, Int_t eta, Int_t phi, Int_t EtaDef, Int_t PhiDef, Float_t Et_Holder, Int_t LayerNum) {
	Int_t Eta_LowerBound = 0, Eta_UpperBound = 0, Phi_LowerBound = 0, Phi_UpperBound = 0;
	GetEtaBounds(Eta_LowerBound, Eta_UpperBound, EtaDef, LayerNum);
	GetPhiBounds(Phi_LowerBound, Phi_UpperBound, PhiDef);
	if (eta >= Eta_LowerBound && eta <= Eta_UpperBound && phi >= Phi_LowerBound && phi <= Phi_UpperBound) {
		Et += Et_Holder;
	}
}

// Add energies of current cells to the reconstruction energy with a given definition
// RecoDefs values must be as follows: L0 eta, L0 phi, L1 eta, L1 phi, L2 eta, L2 phi, L3 eta, L3 phi, Had eta, Had phi
void AddRecoContribution(Int_t eta, Int_t phi, Float_t& Reco_Et, Int_t RecoDef[10], Float_t L0_Reco_Et_Holder, Float_t L1_Reco_Et_Holder,
	Float_t L2_Reco_Et_Holder, Float_t L3_Reco_Et_Holder, Float_t Had_Reco_Et_Holder) {
	// Add L0 piece
	AddLayerContribution(Reco_Et, eta, phi, RecoDef[0], RecoDef[1], L0_Reco_Et_Holder, 0);
	// Add L1 piece
	AddLayerContribution(Reco_Et, eta, phi, RecoDef[2], RecoDef[3], L1_Reco_Et_Holder, 1);
	// Add L2 piece
	AddLayerContribution(Reco_Et, eta, phi, RecoDef[4], RecoDef[5], L2_Reco_Et_Holder, 2);
	// Add L3 piece
	AddLayerContribution(Reco_Et, eta, phi, RecoDef[6], RecoDef[7], L3_Reco_Et_Holder, 3);
	// Add Had piece
	AddLayerContribution(Reco_Et, eta, phi, RecoDef[8], RecoDef[9], Had_Reco_Et_Holder, 4);
}

// Add energies of current cells to the reconstruction energy with definition 1 (1x1, 5x3, 5x3, 3x3, 3x3)
void AddReco1Contribution(Int_t eta, Int_t phi, Float_t& Reco1_Et, Float_t L0_Reco_Et_Holder,
	Float_t L1_Reco_Et_Holder, Float_t L2_Reco_Et_Holder, Float_t L3_Reco_Et_Holder, Float_t Had_Reco_Et_Holder) {
	if (eta < 3 && phi < 3) {
		Reco1_Et += L3_Reco_Et_Holder;
		Reco1_Et += Had_Reco_Et_Holder;
	}
	if (eta == 1 && phi == 1) {
		Reco1_Et += L0_Reco_Et_Holder;
	}
	if (eta > 3 && eta < 9) {
		Reco1_Et += L1_Reco_Et_Holder;
		Reco1_Et += L2_Reco_Et_Holder;
	}
}

// Add energies of current cells to the reconstruction energy with definition 2 (3x3, 7x3, 5x3, 3x3, 3x3)
void AddReco2Contribution(Int_t eta, Int_t phi, Float_t& Reco2_Et, Float_t L0_Reco_Et_Holder,
	Float_t L1_Reco_Et_Holder, Float_t L2_Reco_Et_Holder, Float_t L3_Reco_Et_Holder, Float_t Had_Reco_Et_Holder) {
	if (eta < 3) {
		Reco2_Et += L0_Reco_Et_Holder;
		Reco2_Et += L3_Reco_Et_Holder;
		Reco2_Et += Had_Reco_Et_Holder;
	}
	if (eta > 2 && eta < 10) {
		Reco2_Et += L1_Reco_Et_Holder;
	}
	if (eta > 3 && eta < 9) {
		Reco2_Et += L2_Reco_Et_Holder;
	}
}

// Add energies of current cells to the reconstruction energy with definition 3 (1x2, 5x2, 5x2, 1x2, 3x3)
void AddReco3Contribution(Int_t eta, Int_t phi, Float_t& Reco3_Et, Float_t L0_Reco_Et_Holder, 
	Float_t L1_Reco_Et_Holder, Float_t L2_Reco_Et_Holder, Float_t L3_Reco_Et_Holder, Float_t Had_Reco_Et_Holder) {
	if (eta < 3 && phi < 3) {
		Reco3_Et += Had_Reco_Et_Holder;
	}
	if (eta == 1 && phi < 2) {
		Reco3_Et += L0_Reco_Et_Holder;
		Reco3_Et += L3_Reco_Et_Holder;
	}
	if (eta > 3 && eta < 9 && phi < 2) {
		Reco3_Et += L1_Reco_Et_Holder;
		Reco3_Et += L2_Reco_Et_Holder;
	}
}

// Add energies of current cells to the reconstruction energy with definition 3 (1x3, 5x3, 5x3, 1x3, 3x3)
void AddReco3AContribution(Int_t eta, Int_t phi, Float_t& Reco3A_Et, Float_t L0_Reco_Et_Holder,
	Float_t L1_Reco_Et_Holder, Float_t L2_Reco_Et_Holder, Float_t L3_Reco_Et_Holder, Float_t Had_Reco_Et_Holder) {
	if (eta < 3) {
		Reco3A_Et += Had_Reco_Et_Holder;
	}
	if (eta == 1) {
		Reco3A_Et += L0_Reco_Et_Holder;
		Reco3A_Et += L3_Reco_Et_Holder;
	}
	if (eta > 3 && eta < 9) {
		Reco3A_Et += L1_Reco_Et_Holder;
		Reco3A_Et += L2_Reco_Et_Holder;
	}
}

// Add energies of current cells to the reconstruction energy with definition 4 (1x2, 5x2, 5x2, 1x2, None)
void AddReco4Contribution(Int_t eta, Int_t phi, Float_t& Reco4_Et, Float_t L0_Reco_Et_Holder,
	Float_t L1_Reco_Et_Holder, Float_t L2_Reco_Et_Holder, Float_t L3_Reco_Et_Holder, Float_t Had_Reco_Et_Holder) {
	if (eta == 1 && phi < 2) {
		Reco4_Et += L0_Reco_Et_Holder;
		Reco4_Et += L3_Reco_Et_Holder;
	}
	if (eta > 3 && eta < 9 && phi < 2) {
		Reco4_Et += L1_Reco_Et_Holder;
		Reco4_Et += L2_Reco_Et_Holder;
	}
}

// Add energies of current cells to the reconstruction energy with definition 5 (1x2, 5x2, 5x2, None, None)
void AddReco5Contribution(Int_t eta, Int_t phi, Float_t& Reco5_Et, Float_t L0_Reco_Et_Holder,
	Float_t L1_Reco_Et_Holder, Float_t L2_Reco_Et_Holder, Float_t L3_Reco_Et_Holder, Float_t Had_Reco_Et_Holder) {
	if (eta == 1 && phi < 2) {
		Reco5_Et += L0_Reco_Et_Holder;
	}
	if (eta > 3 && eta < 9 && phi < 2) {
		Reco5_Et += L1_Reco_Et_Holder;
		Reco5_Et += L2_Reco_Et_Holder;
	}
}

// Add energies of current cells to the reconstruction energy with definition 6 (None, 5x2, 5x2, None, None)
void AddReco6Contribution(Int_t eta, Int_t phi, Float_t& Reco6_Et, Float_t L0_Reco_Et_Holder,
	Float_t L1_Reco_Et_Holder, Float_t L2_Reco_Et_Holder, Float_t L3_Reco_Et_Holder, Float_t Had_Reco_Et_Holder) {
	if (eta > 3 && eta < 9 && phi < 2) {
		Reco6_Et += L1_Reco_Et_Holder;
		Reco6_Et += L2_Reco_Et_Holder;
	}
}

// Add energies of current cells to the reconstruction energy with all avaiable cells (3x3, 13x3, 13x3, 3x3, 3x3)
void AddRecoAllContribution(Int_t eta, Int_t phi, Float_t& RecoAll_Et, Float_t L0_Reco_Et_Holder,
	Float_t L1_Reco_Et_Holder, Float_t L2_Reco_Et_Holder, Float_t L3_Reco_Et_Holder, Float_t Had_Reco_Et_Holder) {
	if (eta < 3) {
		RecoAll_Et += Had_Reco_Et_Holder;
		RecoAll_Et += L0_Reco_Et_Holder;
		RecoAll_Et += L3_Reco_Et_Holder;
	}
	RecoAll_Et += L1_Reco_Et_Holder;
	RecoAll_Et += L2_Reco_Et_Holder;
}

// For each event, increment appropriate Et counters
void IncrementEtCounters(Int_t entry, Int_t netcuts, Float_t Reco_Et[], Int_t N_EtCut[]) {
	for (Int_t l = 0; l < netcuts; l++) {
		if (Reco_Et[entry] > l)
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

// Avoid dividing by zero by setting denominator to 1 if would be 0
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

// Output information on the highest Et cut that retains at least a 90% signal efficiency
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

// Check if the provided eta dimension is valid for a 3x3 layer
int Valid3x3EtaDim(Int_t RecoDefDim) {
	if (RecoDefDim == 0 || RecoDefDim == 1 || RecoDefDim == 3) return 0;
	else return 1;
}

// Check if the provided eta dimension is valid for a 13x3 layer
int Valid13x3EtaDim(Int_t RecoDefDim) {
	if (RecoDefDim == 0 || RecoDefDim == 1 || RecoDefDim == 3 || RecoDefDim == 5 || RecoDefDim == 7 || RecoDefDim == 9 || RecoDefDim == 11 || RecoDefDim == 13) return 0;
	else return 1;
}

// Check if the provided phi dimension is valid, does not depend on layer
int ValidPhiDim(Int_t RecoDefDim) {
	if (RecoDefDim < 0 || RecoDefDim > 3) return 1;
	else return 0;
}

void SetRecoDef(Int_t RecoDef[10], Int_t i, Int_t j, Int_t k, Int_t l, Int_t m, Int_t n, Int_t p, Int_t q, Int_t r, Int_t s) {
	RecoDef[0] = i;
	RecoDef[1] = j;
	RecoDef[2] = k;
	RecoDef[3] = l;
	RecoDef[4] = m;
	RecoDef[5] = n;
	RecoDef[6] = p;
	RecoDef[7] = q;
	RecoDef[8] = r;
	RecoDef[9] = s;
}

// Given a reconstructed energy definition return a ROC curve and 90% efficiency info
// RecoCellDefs values must be as follows: L0 eta, L0 phi, L1 eta, L1 phi, L2 eta, L2 phi, L3 eta, L3 phi, Had eta, Had phi
// IsFirstOrLast == 0 marks the first, == 1 the middle, == 2 the end
void GetRocAndEfficiencyFromRecoDef(TTree* tsig, TTree* tback, TCanvas* c1, Int_t RecoCellDef[10], ofstream& textfile, string RecoName, Int_t IsFirstOrLast) {
	// Do some error handling to make sure the definition makes sense
	if (Valid3x3EtaDim(RecoCellDef[0]) + ValidPhiDim(RecoCellDef[1]) + Valid13x3EtaDim(RecoCellDef[2]) + ValidPhiDim(RecoCellDef[3]) + Valid13x3EtaDim(RecoCellDef[4]) +
		ValidPhiDim(RecoCellDef[5]) + Valid3x3EtaDim(RecoCellDef[6]) + ValidPhiDim(RecoCellDef[7]) + Valid3x3EtaDim(RecoCellDef[8]) + ValidPhiDim(RecoCellDef[9]) > 0) {
		cout << "Invalid definition for reconstructed energy" << endl;
		return;
	}

	// Define background variables to read in
	Float_t L0_Back_Et[3][3];
	Float_t L1_Back_Et[13][3];
	Float_t L2_Back_Et[13][3];
	Float_t L3_Back_Et[3][3];
	Float_t Had_Back_Et[3][3];

	// Define signal variables to read in
	Float_t L0_Sig_Et[3][3];
	Float_t L1_Sig_Et[13][3];
	Float_t L2_Sig_Et[13][3];
	Float_t L3_Sig_Et[3][3];
	Float_t Had_Sig_Et[3][3];

	// Make assignemnt between tree and declared variables for background
	tback->SetBranchAddress("L0CellEt[3][3]", &L0_Back_Et[0][0]);  // Assigns to address of first array element
	tback->SetBranchAddress("L1CellEt[13][3]", &L1_Back_Et[0][0]);  // Assigns to address of first array element
	tback->SetBranchAddress("L2CellEt[13][3]", &L2_Back_Et[0][0]);  // Assigns to address of first array element
	tback->SetBranchAddress("L3CellEt[3][3]", &L3_Back_Et[0][0]);  // Assigns to address of first array element
	tback->SetBranchAddress("HadCellEt[3][3]", &Had_Back_Et[0][0]);  // Assigns to address of first array element

	// Make assignment between tree and declared variables for signal
	tsig->SetBranchAddress("L0CellEt[3][3]", &L0_Sig_Et[0][0]);  // Assigns to address of first array element
	tsig->SetBranchAddress("L1CellEt[13][3]", &L1_Sig_Et[0][0]);  // Assigns to address of first array element
	tsig->SetBranchAddress("L2CellEt[13][3]", &L2_Sig_Et[0][0]);  // Assigns to address of first array element
	tsig->SetBranchAddress("L3CellEt[3][3]", &L3_Sig_Et[0][0]);  // Assigns to address of first array element
	tsig->SetBranchAddress("HadCellEt[3][3]", &Had_Sig_Et[0][0]);  // Assigns to address of first array element
	 
	// Create string variables holding reconstructed energy definitions
	string Reco_DefString = "("+to_string(RecoCellDef[0])+"x"+to_string(RecoCellDef[1])+") , ("+to_string(RecoCellDef[2])+"x"+to_string(RecoCellDef[3])+") , ("+
		to_string(RecoCellDef[4])+"x"+to_string(RecoCellDef[5])+") , ("+to_string(RecoCellDef[6])+"x"+to_string(RecoCellDef[7])+") , ("+to_string(RecoCellDef[8])+
		"x"+to_string(RecoCellDef[9])+")";

	// Signal variables
	Int_t sigentries = (Int_t)tsig->GetEntries();
	Float_t Reco_Sig_Et[sigentries];

	// Background variables
	Int_t backentries = (Int_t)tback->GetEntries();
	Float_t Reco_Back_Et[backentries];

	// Number of events above Et threshold
	const Int_t netcuts = 150;
	Int_t N_Reco_Back_EtCut[netcuts] = { 0 };
	Int_t N_Reco_Sig_EtCut[netcuts] = { 0 };


	// Loop over all signal entries
	for (Int_t i = 0; i < sigentries; i++) {
		tsig->GetEntry(i);

		Reco_Sig_Et[i] = 0;

		//TrueEt[i] = mctau->Et() / 1000.;
		//if (TrueEt[i] < 20.) continue;

		// Set flag for preprocessing if necessary by considering sum of off-center phi cells for all layers
		Int_t SigFlip = 0;
		SigFlip = PhiFlip_Et(L0_Sig_Et, L1_Sig_Et, L2_Sig_Et, L3_Sig_Et, Had_Sig_Et);

		// Loop over all cells to calculate average and reconstructed energies
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

				// Calculate reconstructed energy
				AddRecoContribution(j, k, Reco_Sig_Et[i], RecoCellDef, L0_Reco_Et_Holder, L1_Reco_Et_Holder, L2_Reco_Et_Holder,
					L3_Reco_Et_Holder, Had_Reco_Et_Holder);
			}
		}

		// Convert MeV to GeV
		Reco_Sig_Et[i] /= 1000.;

		// For each event, increment appropriate Et counters
		IncrementEtCounters(i, netcuts, Reco_Sig_Et, N_Reco_Sig_EtCut);
	}

	// Now loop over all background entries
	for (Int_t i = 0; i<backentries; i++) {
		tback->GetEntry(i);  // This fills declared variables from ntuple

		Reco_Back_Et[i] = 0;

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

				// Calculate reconstructed energy
				AddRecoContribution(j, k, Reco_Back_Et[i], RecoCellDef, L0_Reco_Et_Holder, L1_Reco_Et_Holder, L2_Reco_Et_Holder,
					L3_Reco_Et_Holder, Had_Reco_Et_Holder);
			}
		}

		// Convert reconstructed energy from MeV to GeV
		Reco_Back_Et[i] /= 1000.;

		// For each event, increment appropriate Et counters
		IncrementEtCounters(i, netcuts, Reco_Back_Et, N_Reco_Back_EtCut);
	}
	
	Float_t Reco_SigEtCutEfficiency[netcuts];
	Float_t Reco_BackEtCutEfficiency[netcuts];

	for (Int_t i = 0; i < netcuts; i++) {
		RecoCutEfficiencies(Reco_SigEtCutEfficiency[i], Reco_BackEtCutEfficiency[i], N_Reco_Sig_EtCut[i], N_Reco_Sig_EtCut[0],
			N_Reco_Back_EtCut[i], N_Reco_Back_EtCut[0]);
	}

	Output90PercentSignalInfo(textfile, "Reco " + Reco_DefString, netcuts, Reco_SigEtCutEfficiency, Reco_BackEtCutEfficiency);

	TGraph *groc = new TGraph(netcuts, Reco_BackEtCutEfficiency, Reco_SigEtCutEfficiency);
	groc->SetName(("reco_"+RecoName+"_roc").c_str());
	groc->SetTitle(("Reco "+RecoName+" ROC").c_str());
	groc->GetXaxis()->SetTitle("Background Efficiency");
	groc->GetXaxis()->SetLimits(0, 0.2);
	groc->GetYaxis()->SetTitle("Signal Efficiency");

	groc->Draw("ALP");
	if (IsFirstOrLast == 0) {
		c1->Print("TauAllRecos.pdf(");
	}
	else if (IsFirstOrLast == 1) {
		c1->Print("TauAllRecos.pdf");
	}
	else if (IsFirstOrLast == 2) {
		c1->Print("TauAllRecos.pdf)");
	}
}

void CombineL1AndL2Et(Float_t CombinedEt[13][3], Float_t L1Et[13][3], Float_t L2Et[13][3]) {
	// j = eta, k = phi
	for (Int_t j = 0; j < 13; j++) {
		for (Int_t k = 0; k < 3; k++) {
			CombinedEt[j][k] = L1Et[j][k] + L2Et[j][k];
		}
	}
}

// Calculates FCore for a given entry
float CalculateFCore(Float_t L1_Et[13][3], Float_t L2_Et[13][3], Int_t CoreDef[2], Int_t IsolationDef[2]) {
	Float_t Core_Et = 0;
	Float_t Isolation_Et = 0;
	Float_t Combined_Et[13][3] = { 0 };

	// Add elements of two layers cell-wise to get composite "layer"
	CombineL1AndL2Et(Combined_Et, L1_Et, L2_Et);

	// Set flag for preprocessing if necessary by considering sum of off-center phi cells for all layers
	Int_t Flip = 0;
	if (Combined_Et[6][2] > Combined_Et[6][0]) {
		Flip = 1;
	}

	// j = eta, k = phi
	for (Int_t j = 0; j < 13; j++) {
		for (Int_t k = 0; k < 3; k++) {
			Float_t Et_Holder;

			if (Flip == 0) {
				Et_Holder = Combined_Et[j][k];
			}
			else {
				Et_Holder = Combined_Et[j][2 - k];
			}

			// Passing 1 to LayerNum to get 13x3 behavior
			AddLayerContribution(Core_Et, j, k, CoreDef[0], CoreDef[1], Et_Holder, 1);
			AddLayerContribution(Isolation_Et, j, k, IsolationDef[0], IsolationDef[1], Et_Holder, 1);
		}
	}

	return (Core_Et / Isolation_Et);
}

void CalculateFCoreAndFillHisto(TTree* tree, Int_t numofentries, Float_t L1_Et[13][3], Float_t L2_Et[13][3], Int_t CoreDef[2], Int_t IsolationDef[2], TH1F* histo) {
	for (Int_t i = 0; i < numofentries; i++) {
		tree->GetEntry(i);

		Float_t FCore = CalculateFCore(L1_Et, L2_Et, CoreDef, IsolationDef);
		histo->Fill(FCore);
	}
}

// Finds the bin of the given histogram above which 95% of the signal resides
int Find95PercentBin(TH1F* histo) {
	Int_t numofbins = histo->GetNbinsX();
	Int_t numofentries = histo->GetEntries();
	Int_t CumulativeEvents = 0;
	Int_t Bin95Percent;
	for (Int_t i = 0; i <= numofbins; i++) {
		CumulativeEvents += histo->GetBinContent(i);
		if (float(CumulativeEvents) / float(numofentries) > 0.05) {
			return i;
		}
	}
	cout << "No 95% bin found" << endl;
	return 0;
}

// Return percentage of events in the histogram above the given bin
float GetEfficiencyAboveBin(Int_t bin, TH1F* histo) {
	Int_t numofentries = histo->GetEntries();
	Int_t CumulativeEvents = 0;
	for (Int_t i = 0; i <= bin; i++) {
		CumulativeEvents += histo->GetBinContent(i);
	}
	return 1 - (float(CumulativeEvents) / float(numofentries));
}

void EvaluateFCoreForDef(TTree* tsig, TTree* tback, TCanvas* c1, TH1F* sighisto, TH1F* backhisto, Int_t CoreDef[2], Int_t IsolationDef[2], Float_t& Bin95PercentSignal, ofstream& textfile) {
	// Do some error handling to make sure the definition makes sense
	if (Valid13x3EtaDim(CoreDef[0]) + ValidPhiDim(CoreDef[1]) + Valid13x3EtaDim(IsolationDef[0]) + ValidPhiDim(IsolationDef[1]) > 0) {
		cout << "Invalid definition for core or isolation energy" << endl;
		return;
	}

	// Define background variables to read in
	Float_t L1_Back_Et[13][3];
	Float_t L2_Back_Et[13][3];

	// Define signal variables to read in
	Float_t L1_Sig_Et[13][3];
	Float_t L2_Sig_Et[13][3];

	// Make assignemnt between tree and declared variables for background
	tback->SetBranchAddress("L1CellEt[13][3]", &L1_Back_Et[0][0]);  // Assigns to address of first array element
	tback->SetBranchAddress("L2CellEt[13][3]", &L2_Back_Et[0][0]);  // Assigns to address of first array element

	// Make assignment between tree and declared variables for signal
	tsig->SetBranchAddress("L1CellEt[13][3]", &L1_Sig_Et[0][0]);  // Assigns to address of first array element
	tsig->SetBranchAddress("L2CellEt[13][3]", &L2_Sig_Et[0][0]);  // Assigns to address of first array element

	// Signal variables
	Int_t sigentries = (Int_t)tsig->GetEntries();
	Float_t sigentriesfloat = sigentries;

	// Background variables
	Int_t backentries = (Int_t)tback->GetEntries();
	Float_t backentriesfloat = backentries;

	CalculateFCoreAndFillHisto(tsig, sigentries, L1_Sig_Et, L2_Sig_Et, CoreDef, IsolationDef, sighisto);
	CalculateFCoreAndFillHisto(tback, backentries, L1_Back_Et, L2_Back_Et, CoreDef, IsolationDef, backhisto);

	Bin95PercentSignal = Find95PercentBin(sighisto) * 0.01;

	textfile << "Signal Efficiency: " << GetEfficiencyAboveBin(Bin95PercentSignal, sighisto) << endl;
	textfile << "Background Efficiency: " << GetEfficiencyAboveBin(Bin95PercentSignal, backhisto) << endl << endl;


	sighisto->Scale(1 / sigentriesfloat);
	backhisto->Scale(1 / backentriesfloat);
}

void SetCoreAndIsolationDefs(Int_t CoreEta, Int_t CorePhi, Int_t IsolationEta, Int_t IsolationPhi, Int_t CoreDef[2], Int_t IsolationDef[2]) {
	CoreDef[0] = CoreEta;
	CoreDef[1] = CorePhi;
	IsolationDef[0] = IsolationEta;
	IsolationDef[1] = IsolationPhi;
}

void OutputFCoreTextAndHisto(TTree* tsig, TTree* tback, TCanvas* c1, ofstream& textfile, Int_t CoreDef[2], Int_t IsolationDef[2], Int_t IsFirstOrLast) {
	textfile << "Core definition: " << CoreDef[0] << "x" << CoreDef[1] << "	Isolation definition: " << IsolationDef[0] << "x" << IsolationDef[1] << endl;
	string histotitle = "FCore (" + to_string(CoreDef[0]) + "x" + to_string(CoreDef[1]) + ") (" + to_string(IsolationDef[0]) + "x" + to_string(IsolationDef[1]) + ")";

	// Define histograms
	TH1F* sighisto = new TH1F("sighisto", "FCore; FCore; Entries", 100, 0., 1.);
	TH1F* backhisto = new TH1F("backhisto", "FCore; FCore; Entries", 100, 0., 1.);

	Float_t Bin95PercentSignal = 0;

	EvaluateFCoreForDef(tsig, tback, c1, sighisto, backhisto, CoreDef, IsolationDef, Bin95PercentSignal, textfile);

	sighisto->Draw("hist");
	sighisto->SetLineColor(kRed);
	sighisto->SetTitle(histotitle.c_str());
	backhisto->Draw("same hist");
	backhisto->SetLineColor(kBlue);
	TLine *line = new TLine(Bin95PercentSignal, 0, Bin95PercentSignal, backhisto->GetMaximum());
	line->Draw();

	if (IsFirstOrLast == 0) {
		c1->Print("TauFCore.pdf(");
	}
	else if (IsFirstOrLast == 1) {
		c1->Print("TauFCore.pdf");
	}
	else {
		c1->Print("TauFCore.pdf)");
	}
}