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
//
void Uspo::RunAnalysis() {

	cout << "Starting RunAnalysis()..." << endl;

	//-------------------------------------------------------------------------------

	cout << "\n All done! :) \n"

}