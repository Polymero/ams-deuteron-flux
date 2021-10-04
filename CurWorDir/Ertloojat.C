// Written by Sebastiaan Venendaal (University of Groningen, the Netherlands)
// C++ class for histogram analysis of AMS-02 data, used for p/d flux analysis
// Created 			02-06-21
// Last modified 	06-06-21

// Include header files
#include <iostream>
#include <string>
#include <algorithm>
#include "Header Files/Ntp.h"
#include "TChain.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TString.h"

//-----------------------------------------------------------------------------------
// CLASS DEFINITION
//-----------------------------------------------------------------------------------

// Class
class Uspo {
	public: // Access Specifier

	//-------------------------------------------------------------------------------
	// ATTRIBUTES
	//-------------------------------------------------------------------------------
	// Rigidity bins
	const int Bin_num = 32;
	double Bin_edges[32+1] = {1.00,1.16,1.33,1.51,1.71,1.92,2.15,2.40,2.67,2.97,
							  3.29,3.64,4.02,4.43,4.88,5.37,5.90,6.47,7.09,7.76,
							  8.48,9.26,10.1,11.0,12.0,13.0,14.1,15.3,16.6,18.0,
							  19.5,21.1,22.8};
	double Bin_err[32];
	double Bin_mid[32];

	// Rigidity cut-off level
	double Rig_Cut_Level = 1.2;

	//-------------------------------------------------------------------------------
	// CONSTRUCTORS
	//-------------------------------------------------------------------------------
	Uspo() { // Default constructor

		// gStyle
		gStyle->SetOptTitle(0);
		gStyle->SetOptStat(0);
		gStyle->SetOptLogx(1);

		// Bin properties
		for (int i=0; i<Bin_num; i++) {
			Bin_err[i] = (Bin_edges[i+1] - Bin_edges[i]) / 2;
			Bin_mid[i] = (Bin_edges[i+1] + Bin_edges[i]) / 2;
		}

		// Retrieve histogram ROOT file $ Hischaajat.C
		cout << "Trying to load Histogram file..." << flush
		TFile *Hist_file = new TFile("../test.root");
		cout << "   ...File loaded!" << endl;

		cout << "\nClass succesfully constructed!\n" << endl;

	};

	//-------------------------------------------------------------------------------
	// METHODS
	//-------------------------------------------------------------------------------
	void RunAnalysis();

};

//-----------------------------------------------------------------------------------
// SUPPORT FUNCTIONS
//-----------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------
// METHOD FUNCTIONS
//-----------------------------------------------------------------------------------
// Saves TGraph graphs
void Uspo::RunAnalysis() {

	cout << "Starting RunAnalysis()..." << endl;

	//-------------------------------------------------------------------------------
	cout << "Creating Events graph... (1/?)" << endl;

	// Get histograms
	H_Events_raw  = (TH1F*)Hist_file->Get("Events_raw");
	H_Events_pcut = (TH1F*)Hist_file->Get("Events_pcut");
	H_Events_dcut = (TH1F*)Hist_file->Get("Events_dcut");

	// Set arrays
	double Events_raw[Bin_num]; double Events_raw_err[Bin_num];
	double Events_pcut[Bin_num]; double Events_pcut_err[Bin_num]; 
	double Events_dcut[Bin_num]; double Events_dcut_err[Bin_num];
	for (int i=0; i<Bin_num; i++) {
		Events_raw[i]  = H_Events_raw->GetBinContent(i+1);
		Events_pcut[i] = H_Events_pcut->GetBinContent(i+1);
		Events_dcut[i] = H_Events_dcut->GetBinContent(i+1);
		Events_raw_err[i]  = TMath::Sqrt(Events_raw[i]);
		Events_pcut_err[i] = TMath::Sqrt(Events_pcut[i]);
		Events_dcut_err[i] = TMath::Sqrt(Events_dcut[i]);
	}

	// TGraphs
	TGraphErrors *G_Events_raw  = new TGraphErrors(Bin_num, Bin_mid, Events_raw, Bin_err, Events_raw_err);
	TGraphErrors *G_Events_pcut = new TGraphErrors(Bin_num, Bin_mid, Events_pcut, Bin_err, Events_pcut_err);
	TGraphErrors *G_Events_dcut = new TGraphErrors(Bin_num, Bin_mid, Events_dcut, Bin_err, Events_dcut_err);

	// Canvas
	TCanvas *C_Events = new TCanvas("C_Events", "Events per Rigidity Bin");
	G_Events_raw->Draw("AP");
	G_Events_pcut->Draw("P");
	G_Events_dcut->Draw("P");
	// Styling
	G_Events_raw->SetMarkerStyle(20);
	G_Events_raw->SetMarkerSize(1);
	G_Events_raw->SetMarkerColor(kBlack);
	G_Events_pcut->SetMarkerStyle(20);
	G_Events_pcut->SetMarkerSize(1);
	G_Events_pcut->SetMarkerColor(kBlue);
	G_Events_dcut->SetMarkerStyle(20);
	G_Events_dcut->SetMarkerSize(1);
	G_Events_dcut->SetMarkerColor(kGreen);
	// Axes
	C_Events->SetLogy();
	G_Events_raw->SetMinimum(1);
	G_Events_raw->GetXaxis()->SetTitle("R [GV]");
	G_Events_raw->GetYaxis()->SetTitle("Events");
	// Print
	C_Events->Draw();
	C_Events->Print((outdir + "Events.png").c_str());



	//-------------------------------------------------------------------------------
	cout << "Creating ExposureTime graph... (2/?)" << endl;

	// Get histograms
	H_ExposureTime = (TH1F*)Hist_file->Get("ExposureTime");

	// Set arrays
	double ExposureTime[Bin_num];
	for (int i=0; i<Bin_num; i++) {
		ExposureTime[i] = H_ExposureTime->GetBinContent(i+1);
	}

	// TGraphs
	TGraphErrors *G_ExposureTime = new TGraphErrors(Bin_num, Bin_mid, ExposureTime, Bin_err, 0);

	// Canvas
	TCanvas *C_ExposureTime = new TCanvas("C_ExposureTime", "Exposure Time per Rigidity Bin");
	G_ExposureTime->Draw("AP");
	// Styling
	G_ExposureTime->SetMarkerSize(20);
	G_ExposureTime->SetMarkerSize(1);
	G_ExposureTime->SetMarkerColor(kRed);
	// Axes
	C_ExposureTime->SetLogy();
	G_ExposureTime->GetXaxis()->SetTitle("R [GV]");
	G_ExposureTime->GetYaxis()->SetTitle("Exposure Time [s]");
	// Print
	C_ExposureTime->Draw();
	C_ExposureTime->Print((outdir + "ExposureTime.png").c_str());



	//-------------------------------------------------------------------------------
	cout << "Creating Acceptance graphs... (3/?)" << endl;

	// Get histograms
	H_MCp_gen = (TH1F*)Hist_file->Get("MCp_gen");
	H_MCp_det = (TH1F*)Hist_file->Get("MCp_det");
	H_MCp_sel = (TH1F*)Hist_file->Get("MCp_sel");
	H_MCd_gen = (TH1F*)Hist_file->Get("MCd_gen");
	H_MCd_det = (TH1F*)Hist_file->Get("MCd_det");
	H_MCd_sel = (TH1F*)Hist_file->Get("MCd_sel");

	// Set arrays
	double pAcceptance_det[Bin_num];
	double pAcceptance_sel[Bin_num];
	double dAcceptance_det[Bin_num];
	double dAcceptance_sel[Bin_num];
	for (int i=0; i<Bin_num; i++) {
		pAcceptance_det[i] = TMath::Pi()*3.9*3.9 * H_MCp_det->GetBinContent(i+1) / H_MCp_gen->GetBinContent(i+1);
		pAcceptance_sel[i] = TMath::Pi()*3.9*3.9 * H_MCp_sel->GetBinContent(i+1) / H_MCp_gen->GetBinContent(i+1);
		dAcceptance_det[i] = TMath::Pi()*3.9*3.9 * H_MCd_det->GetBinContent(i+1) / H_MCd_gen->GetBinContent(i+1);
		dAcceptance_sel[i] = TMath::Pi()*3.9*3.9 * H_MCd_sel->GetBinContent(i+1) / H_MCd_gen->GetBinContent(i+1);
	}

	// TGraphs
	TGraphErrors *G_pAcceptance_det = new TGraphErrors(Bin_num, Bin_mid, pAcceptance_det, Bin_err, 0);
	TGraphErrors *G_pAcceptance_sel = new TGraphErrors(Bin_num, Bin_mid, pAcceptance_sel, Bin_err, 0);
	TGraphErrors *G_dAcceptance_det = new TGraphErrors(Bin_num, Bin_mid, dAcceptance_det, Bin_err, 0);
	TGraphErrors *G_dAcceptance_sel = new TGraphErrors(Bin_num, Bin_mid, dAcceptance_sel, Bin_err, 0);

	// Canvas 1
	TCanvas *C_Acceptance_det = new TCanvas("C_Acceptance_det", "Uncut Acceptance per Rigidity Bin");
	G_pAcceptance_det->Draw("AP");
	G_dAcceptance_det->Draw("P");
	// Styling
	G_pAcceptance_det->SetMarkerStyle(20);
	G_pAcceptance_det->SetMarkerSize(1);
	G_pAcceptance_det->SetMarkerColor(kBlue);
	G_dAcceptance_det->SetMarkerStyle(20);
	G_dAcceptance_det->SetMarkerSize(1);
	G_dAcceptance_det->SetMarkerColor(kGreen);
	// Axes
	G_pAcceptance_det->SetMinimum(0);
	G_pAcceptance_det->GetXaxis()->SetTitle("R [GV]");
	G_pAcceptance_det->GetYaxis()->SetTitle("Acceptance [m^2 sr]");
	// Print
	C_Acceptance_det->Draw();
	C_Acceptance_det->Print((outdir + "AcceptanceUncut.png").c_str());

	// Canvas 2
	TCanvas *C_Acceptance_sel = new TCanvas("C_Acceptance_sel", "Cut Acceptance per Rigidity Bin");
	G_pAcceptance_sel->Draw("AP");
	G_dAcceptance_sel->Draw("P");
	// Styling
	G_pAcceptance_sel->SetMarkerStyle(20);
	G_pAcceptance_sel->SetMarkerSize(1);
	G_pAcceptance_sel->SetMarkerColor(kBlue);
	G_dAcceptance_sel->SetMarkerStyle(20);
	G_dAcceptance_sel->SetMarkerSize(1);
	G_dAcceptance_sel->SetMarkerColor(kGreen);
	// Axes
	G_pAcceptance_sel->SetMinimum(0);
	G_pAcceptance_sel->GetXaxis()->SetTitle("R [GV]");
	G_pAcceptance_sel->GetYaxis()->SetTitle("Acceptance [m^2 sr]");
	// Print
	C_Acceptance_sel->Draw();
	C_Acceptance_sel->Print((outdir + "AcceptanceCut.png").c_str());


	//-------------------------------------------------------------------------------
	cout << "Creating Trigger Efficiency graphs... (4/?)" << endl;

	// Get histograms
	H_P_phT = (TH1F*)Hist_file->Get("P_phT");
	H_P_unT = (TH1F*)Hist_file->Get("P_unT");
	H_D_phT = (TH1F*)Hist_file->Get("D_phT");
	H_D_unT = (TH1F*)Hist_file->Get("D_unT");
	H_MCp_phT = (TH1F*)Hist_file->Get("MCp_phT");
	H_MCp_unT = (TH1F*)Hist_file->Get("MCp_unT");
	H_MCd_phT = (TH1F*)Hist_file->Get("MCd_phT");
	H_MCd_unT = (TH1F*)Hist_file->Get("MCd_unT");

	// Set arrays
 	double p_TE[Bin_num]; double p_TE_err[Bin_num];
 	double d_TE[Bin_num]; double d_TE_err[Bin_num];
 	double MCp_TE[Bin_num]; double MCp_TE_err[Bin_num];
 	double MCd_TE[Bin_num]; double MCd_TE_err[Bin_num];
 	double p_TER[Bin_num]; double p_TER_err[Bin_num];
 	double d_TER[Bin_num]; double d_TER_err[Bin_num];
 	for (int i=0; i<Bin_num; i++) {
 		p_TE[i] = H_P_phT->GetBinContent(i+1) / (H_P_phT->GetBinContent(i+1) + H_P_unT->GetBinContent(i+1));
 		d_TE[i] = H_D_phT->GetBinContent(i+1) / (H_D_phT->GetBinContent(i+1) + H_D_unT->GetBinContent(i+1));
 		MCp_TE[i] = H_MCp_phT->GetBinContent(i+1) / (H_MCp_phT->GetBinContent(i+1) + H_MCp_unT->GetBinContent(i+1));
 		MCd_TE[i] = H_MCd_phT->GetBinContent(i+1) / (H_MCd_phT->GetBinContent(i+1) + H_MCd_unT->GetBinContent(i+1));
 		p_TER[i] = 
 	}

	// Tgraphs


	// Canvas
 	// Styling
 	// Axes
 	// Print



}