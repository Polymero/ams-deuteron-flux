// Written by Sebastiaan Venendaal (University of Groningen, the Netherlands)
// C++ class for histogram analysis of AMS-02 data, used for p/d flux analysis
// Created 			02-06-21
// Last modified 	02-06-21

// Include header files
#include <iostream>
#include <string>
#include <algorithm>
#include "Header Files/Ntp.h"
//#include "Header Files/Simple.h"
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
	cout << "Creating event graphs... (1/?)" << endl;

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


	cout << "\n All done! :) \n"

}