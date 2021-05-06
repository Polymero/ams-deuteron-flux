// Written by Sebastiaan Venendaal (University of Groningen, the Netherlands)
// C++ class for histogram generation of AMS-02 data, used for p/d flux analysis
// Created 			03-03-21
// Last modified 	03-03-21

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
class Mirja {
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

    // Histograms
    // RigBinner()
    // Exposure()
    // Acceptance()
    // SelEff()
    // TrigEff()

    // Data objects
    // Chains
    TChain *Comp_chain		= new TChain("Compact");
    TChain *RTII_chain		= new TChain("RTI");
    TChain *MCps_chain		= new TChain("Compact");
    TChain *MCpi_chain		= new TChain("File");
    TChain *MCds_chain		= new TChain("Compact");
    TChain *MCdi_chain		= new TChain("File");
    // Classes
    NtpCompact *Comp 		= new NtpCompact();
    NtpSheader *SHeader 	= new NtpSHeader();
    RTIInfo *RTIInfo 		= new class RTIInfo();
    NtpCompact *MCp_comp	= new NtpCompact();
    FileMCInfo *MCp_info	= new class FileMCInfo();
    NtpCompact *MCd_comp	= new NtpCompact();
    FileMCInfo *MCd_info	= new class FileMCInfo();

    //-------------------------------------------------------------------------------
	// CONSTRUCTORS
	//-------------------------------------------------------------------------------
	Mirja(int fileloc = 1) { // Default constructor

		// gStyle
		gStyle->SetOptTitle(0);
		gStyle->SetOptStat(0);
		gStyle->SetOptLogx(1);

		// Bin properties
	    for (int i=0; i<Bin_num; i++) {
	    	Bin_err[i] = (Bin_edges[i+1] - Bin_edges[i]) / 2;
	    	Bin_mid[i] = (Bin_edges[i+1] + Bin_edges[i]) / 2;
	    }

	    // Read the trees
	    if (fileloc == 1) { // Local files @ Kapteyn
	    	Comp_chain->Add("/net/dataserver3/data/users/bueno/data/iss/AO/ISS.B1130/pass7/*.root");
	    	RTII_chain->Add("/net/dataserver3/data/users/bueno/data/iss/AO/ISS.B1130/pass7/*.root");
	    	MCps_chain->Add("/net/dataserver3/data/users/bueno/data/mc/Pr.B1200/*.root");
	    	MCpi_chain->Add("/net/dataserver3/data/users/bueno/data/mc/Pr.B1200/*.root");
	    	MCds_chain->Add("/net/dataserver3/data/users/bueno/data/mc/D.B1220/*.root");
	    	MCdi_chain->Add("/net/dataserver3/data/users/bueno/data/mc/D.B1220/*.root");
	    }
	    
	    // Set branch addresses
	    Comp_chain->SetBranchAddress("Compact", &Comp);
	    Comp_chain->SetBranchAddress("SHeader", &SHeader);
	    RTII_chain->SetBranchAddress("RTIInfo", &RTIInfo);
	    MCps_chain->SetBranchAddress("Compact", &MCp_comp);
	    MCpi_chain->SetBranchAddress("FileMCInfo", MCp_info);
	    MCds_chain->SetBranchAddress("Compact", &MCd_comp);
	    MCdi_chain->SetBranchAddress("FileMCInfo", MCd_info);

	    cout << "\nClass succesfully constructed!\n" << endl;

	};

	//-------------------------------------------------------------------------------
	// METHODS
	//-------------------------------------------------------------------------------
	void RunAnalysis(const char* runbit = "1111111");

};

//-----------------------------------------------------------------------------------
// SUPPORT FUNCTIONS
//-----------------------------------------------------------------------------------




//-----------------------------------------------------------------------------------
// METHOD FUNCTIONS
//-----------------------------------------------------------------------------------
// Returns TH1D histograms
void Mirja::RunAnalysis(const char* runbit = "1111111") {

	cout << "Starting RunAnalysis()..." << endl;

	// Create RTI map
	map<int, std::pair<float, float>> rtimap = map<int, std::pair<float, float>>();

	//-------------------------------------------------------------------------------
	cout << "Looping over RTIInfo Data... (1/6)" << endl;
	cout << "Number of RTIInfo Entries: " << RTII_chain->GetEntries() << endl;
	// Looping over RTII files
	for (int i=0; i<RTII_chain->GetEntries(); i++) {

		// Get entry
		RTII_chain->GetEntry(i);

		// Fill RTI map
		rtimap.insert({RTIInfo->utime, std::pair<float, float>(RTIInfo->lf, RTIInfo->cf[0][3][1])});

		// Loop over rigidity bins
		for (int j=0; j<Bin_num; j++) {

			// Exposure() -- Get total livetime as function of rigidity
			// If the centre of the bin is above the geo-magnetic cut-off, add the livetime
			ExposureTime->SetBinContent(j+1, ExposureTime->GetBinContent(j+1) + RTIInfo->lf);
		
		}

	}

	//-------------------------------------------------------------------------------
	cout << "Looping over Compact Data... (2/6)" << endl;
	cout << "Number of Compact Entries: " << Comp_chain->GetEntries() << endl;
	// Loop over Data Compact entries
	for (int i=0; i<Comp_chain->GetEntries(); i++) {

		// Get entry
		Comp_chain->GetEntry(i);

		// Boolean cuts
		// Geomagnetic cut-off
		bool GeoC = Comp->trk_rig[0] > 1.2 * rtimap[SHeader->utime].second;
		// Proton-like events
		bool Prig = (Comp->trk_rig[0] > Bin_edges[0])&&(Comp->trk_rig[0] <= Bin_edges[Bin_num]);
		bool Ptrg = ((Comp->sublvl1 & 0x3E) != 0)&&((Comp->trigpatt & 0x2) != 0);
		bool Ppar = Comp->status % 10 == 1;
		bool Pbet = Comp->tof_beta > 0.3;
		bool Pchi = (Comp->trk_chisqn[0][0] < 10)&&(Comp->trk_chisqn[0][1] < 10)&&(Comp->trk_chisqn[0][0] > 0)&&(Comp->trk_chisqn[0][1] > 0);
		bool Pinn = (Comp->trk_q_inn > 0.80)&&(Comp->trk_q_inn < 1.30);
		// Deuteron restrictions (Aerogel only)
		bool Drig = (Comp->trk_rig[0] > Bin_edges[10])&&(Comp->trk_rig[0] <= Bin_edges[Bin_num]);
		bool Dagl = Comp->rich_select == 2;
    	bool Dcon = std::abs(Comp->tof_beta - Comp->rich_beta)/Comp->tof_beta < 0.05;
    	bool Dlay = Comp->trk_q_lay[0] <= 1.7;
    	// Get boolbit
    	int boolbit = Dlay + (Dcon << 1) + (Dagl << 2) + (Drig << 3) + 
    	              (Pinn << 4) + (Pchi << 5) + (Pbet << 6) + (Ppar << 7) + 
    	              (Ptrg << 8) + (Prig << 9) + (GeoC << 10);

		// RigBinner() -- Bin events as function of rigidity
		Events_raw->Fill(Comp->trk_rig[0]);
		// Protons
		if ((boolbit & 2032) == 2032) {
			Events_pcut->Fill(Comp->trk_rig[0]);
		}
		// Deuterons
		if ((boolbit & 2047) == 2047) { 
			Events_dcut->Fill(Comp->trk_rig[0]);
		}

		// TrigEff(): Data -- Trigger efficiency as function of rigidity
		// Trigger booleans
    	bool HasPhysTrig = ((Comp->sublvl1&0x3E)!=0)&&((Comp->trigpatt&0x2)!=0);
    	bool HasUnphTrig = ((Comp->sublvl1&0x3E)==0)&&((Comp->trigpatt&0x2)!=0);
    	// Fill histograms
    	// Protons
    	if ((boolbit & 1776) == 1776) {
    		if (HasPhysTrig) {
    			P_phT->Fill(Comp->trk_rig[0]);
    		}
    		if (HasUnphTrig) {
    			P_unT->Fill(Comp->trk_rig[0]);
    		}
    	}
    	// Deuterons
    	if ((boolbit & 1791) == 1791) {
    		if (HasPhysTrig) {
    			D_phT->Fill(Comp->trk_rig[0]);
    		}
    		if (HasUnphTrig) {
    			D_unT->Fill(Comp->trk_rig[0]);
    		}
    	}

		// SelEff(): Data -- Selection efficiency of applied cuts as function of rigidity
    	// Additional TOF charge cuts (to replace TRK charge cuts)
    	bool tof_q_p = (Comp->tof_q_lay[0] > 0.8)&&(Comp-tof_q_lay[0] < 1.5)
    	bool tof_q_d = (Comp->tof_d_lay[1] < 1.9)
		// Fill histograms
		// Protons
		// TRK Base
		if ((boolbit & 1856) == 1856)&&(tof_q_p) {
			P_TRK_Base_Data->Fill(Comp->trk_rig[0]);
		}
		// TOF Base
		if ((boolbit & 1968) == 1968) {
			P_TOF_Base_Data->Fill(Comp->trk_rig[0]);
		}
		// Ppar Cut (TRK Base)
		if ((boolbit & 1984) == 1984)&&(tof_q_p) {
			P_Ppar_Cut_Data->Fill(Comp->trk_rig[0]);
		}
		// Pbet Cut (TOF Base)
		if ((boolbit & 2032) == 2032) {
			P_Pbet_Cut_Data->Fill(Comp->trk_rig[0]);
		}
		// Pchi Cut (TRK Base)
		if ((boolbit & 1888) == 1888)&&(tof_q_p) {
			P_Pchi_Cut_Data->Fill(Comp->trk_rig[0]);
		}
		// Pinn Cut (TRK Base)
		if ((boolbit & 1872) == 1872) {
			P_Pinn_Cut_Data->Fill(Comp->trk_rig[0]);
		}
		// Deuterons
		// TRK Base
		if ((boolbit & 1870) == 1870)&&(tof_q_p)&&(tof_q_d) {
			D_TRK_Base_Data->Fill(Comp->trk_rig[0]);
		}
		// TOF Base
		if ((boolbit & 1981) == 1981) {
			D_TOF_Base_Data->Fill(Comp->trk_rig[0]);
		}
		// RCH Base
		if ((boolbit & 2041) == 2041) {
			D_RCH_Base_Data->Fill(Comp->trk_rig[0]);
		}
		// Ppar Cut
		if ((boolbit & 1998) == 1998)&&(tof_q_p)&&(tof_q_d) {
			D_Ppar_Cut_Data->Fill(Comp->trk_rig[0]);
		}
		// Pbet Cut
		if ((boolbit & 2045) == 2045) {
			D_Pbet_Cut_Data->Fill(Comp->trk_rig[0]);
		}
		// Pchi Cut
		if ((boolbit & 1902) == 1902)&&(tof_q_p)&&(tof_q_d) {
			D_Pchi_Cut_Data->Fill(Comp->trk_rig[0]);
		}
		// Pinn Cut
		if ((boolbit & 1886) == 1886)&&(tof_q_d) {
			D_Pinn_Cut_Data->Fill(Comp->trk_rig[0]);
		}
		// Dagl Cut
		if ((boolbit & 2045) == 2045) {
			D_Dagl_Cut_Data->Fill(Comp->trk_rig[0]);
		}
		// Dcon Cut
		if ((boolbit & 2043) == 2043) {
			D_Dcon_Cut_Data->Fill(Comp->trk_rig[0]);
		}
		// Dlay Cut
		if ((boolbit & 1871) == 1871)&&(tof_q_p) {
			D_Dlay_Cut_Data->Fill(Comp->trk_rig[0]);
		}

	}

	//-------------------------------------------------------------------------------
	cout << "Looping over Proton MC Compact Data... (3/6)" << endl;
	cout << "Number of MC Proton Compact Entries: " << MCps_chain->GetEntries() << endl;
	// Loop over MC Proton Compact entries
	for (int i=0; i<MCps_chain->GetEntries(); i++) {

		// Get entry
		MCps_chain->GetEntry(i);

		// Boolean cuts
		// Proton-like events
		bool Prig = (MCp_comp->trk_rig[0] > Bin_edges[0])&&(MCp_comp->trk_rig[0] <= Bin_edges[Bin_num]);
		bool Ptrg = ((MCp_comp->sublvl1 & 0x3E) != 0)&&((MCp_comp->trigpatt & 0x2) != 0);
		bool Ppar = MCp_comp->status % 10 == 1;
		bool Pbet = MCp_comp->tof_beta > 0.3;
		bool Pchi = (MCp_comp->trk_chisqn[0][0] < 10)&&(MCp_comp->trk_chisqn[0][1] < 10)&&(MCp_comp->trk_chisqn[0][0] > 0)&&(MCp_comp->trk_chisqn[0][1] > 0);
		bool Pinn = (MCp_comp->trk_q_inn > 0.80)&&(MCp_comp->trk_q_inn < 1.30);
		// Deuteron restrictions (Aerogel only)
		bool Drig = (MCp_comp->trk_rig[0] > Bin_edges[10])&&(MCp_comp->trk_rig[0] <= Bin_edges[Bin_num]);
		bool Dagl = MCp_comp->rich_select == 2;
    	bool Dcon = std::abs(MCp_comp->tof_beta - MCp_comp->rich_beta)/MCp_comp->tof_beta < 0.05;
    	bool Dlay = MCp_comp->trk_q_lay[0] <= 1.7;
    	// Get boolbit
    	int boolbit = Dlay + (Dcon << 1) + (Dagl << 2) + (Drig << 3) + 
    	              (Pinn << 4) + (Pchi << 5) + (Pbet << 6) + (Ppar << 7) + 
    	              (Ptrg << 8) + (Prig << 9);

		// Acceptance() -- Geometric Aperature Acceptance as function of rigidity
    	// Fill histograms
		MCp_det->Fill(MCp_comp->trk_rig[0]);
		if ((boolbit & 1008) == 1008) {
			MCp_sel->Fill(MCp_comp->trk_rig[0]);
		}

		// TrigEff(): MC -- Trigger efficiency as function of rigidity
		// Trigger booleans
    	bool HasPhysTrig_mc = ((MCp_comp->sublvl1&0x3E)!=0)&&((MCp_comp->trigpatt&0x2)!=0);
    	bool HasUnphTrig_mc = ((MCp_comp->sublvl1&0x3E)==0)&&((MCp_comp->trigpatt&0x2)!=0);
    	// Fill histograms
    	if ((boolbit & 752) == 752) {
    		if (HasPhysTrig_mc) {
    			MCp_phT->Fill(MCp_comp->trk_rig[0]);
    		}
    		if (HasUnphTrig_mc) {
    			MCp_unT->Fill(MCp_comp->trk_rig[0]);
    		}
    	}

		// SelEff(): MC -- Selection efficiency of applied cuts as function of rigidity

	}

	//-------------------------------------------------------------------------------
	cout << "Looping over Deuteron MC Compact Data... (4/6)" << endl;
	cout << "Number of MC Deuteron Compact Entries: " << MCds_chain->GetEntries() << endl;
	// Loop over MC Deuteron Compact entries
	for (int i=0; i<MCds_chain->GetEntries(); i++) {

		// Get entry
		MCds_chain->GetEntry(i);

		// Boolean cuts
		// Proton-like events
		bool Prig = (MCd_comp->trk_rig[0] > Bin_edges[0])&&(MCd_comp->trk_rig[0] <= Bin_edges[Bin_num]);
		bool Ptrg = ((MCd_comp->sublvl1 & 0x3E) != 0)&&((MCd_comp->trigpatt & 0x2) != 0);
		bool Ppar = MCd_comp->status % 10 == 1;
		bool Pbet = MCd_comp->tof_beta > 0.3;
		bool Pchi = (MCd_comp->trk_chisqn[0][0] < 10)&&(MCd_comp->trk_chisqn[0][1] < 10)&&(MCd_comp->trk_chisqn[0][0] > 0)&&(MCp_comp->trk_chisqn[0][1] > 0);
		bool Pinn = (MCd_comp->trk_q_inn > 0.80)&&(MCd_comp->trk_q_inn < 1.30);
		// Deuteron restrictions (Aerogel only)
		bool Drig = (MCd_comp->trk_rig[0] > Bin_edges[10])&&(MCd_comp->trk_rig[0] <= Bin_edges[Bin_num]);
		bool Dagl = MCd_comp->rich_select == 2;
    	bool Dcon = std::abs(MCd_comp->tof_beta - MCd_comp->rich_beta)/MCd_comp->tof_beta < 0.05;
    	bool Dlay = MCd_comp->trk_q_lay[0] <= 1.7;
    	// Get boolbit
    	int boolbit = Dlay + (Dcon << 1) + (Dagl << 2) + (Drig << 3) + 
    	              (Pinn << 4) + (Pchi << 5) + (Pbet << 6) + (Ppar << 7) + 
    	              (Ptrg << 8) + (Prig << 9);

		// Acceptance() -- Geometric Aperature Acceptance as function of rigidity
    	// Fill histograms
		MCd_det->Fill(MCd_comp->trk_rig[0]);
		if ((boolbit & 1023) == 1023) {
			MCd_sel->Fill(MCd_comp->trk_rig[0]);
		}

	    // TrigEff(): MC -- Trigger efficiency as function of rigidity
		// Trigger booleans
    	bool HasPhysTrig_mc = ((MCd_comp->sublvl1&0x3E)!=0)&&((MCd_comp->trigpatt&0x2)!=0);
    	bool HasUnphTrig_mc = ((MCd_comp->sublvl1&0x3E)==0)&&((MCd_comp->trigpatt&0x2)!=0);
    	// Fill histograms
    	if ((boolbit & 767) == 767) {
    		if (HasPhysTrig_mc) {
    			MCd_phT->Fill(MCd_comp->trk_rig[0]);
    		}
    		if (HasUnphTrig_mc) {
    			MCd_unT->Fill(MCd_comp->trk_rig[0]);
    		}
    	}

    	// SelEff(): MC -- Selection efficiency of applied cuts as function of rigidity

	}

	//-------------------------------------------------------------------------------
	cout << "Looping over Proton MC File Info Data... (5/6)" << endl;
	cout << "Number of Proton FileMCInfo Entries: " << MCpi_chain->GetEntries() << endl;
	// Loop over MC Proton FileMCInfo entries
	for (int i=0; i<MCpi_chain->GetEntries(); i++) {

		// Get entry
		MCpi_chain->GetEntry(i);

		// Get MC generation parameters
		double ngen = MCp_info->ngen_datacard;
		double rmin = MCp_comp->momentum[0];
		double rmax = MCp_comp->momentum[1];

		// 1/R generation spectrum
		TF1 *genFlux = new TF1("genFlux", "[0]/(x)", rmin, rmax);
		// Normalisation
		genFlux->SetParameter(0, log(rmax / rmin));

		// Loop over rigidity bins
	    for (int j=0; j<Bin_num; j++) {
	      // Fraction of spectrum in bin
	      double frac = genFlux->Integral(Bin_edges[j], Bin_edges[j+1]) / genFlux->Integral(rmin, rmax);
	      // Number of events in fraction
	      double evnum = frac * ngen;
	      // Set bin to current bin count + number of events
	      MCp_gen->SetBinContent(j+1, MCp_gen->GetBinContent(j+1) + evnum);
	    }

	}

	//-------------------------------------------------------------------------------
	cout << "Looping over Deuteron MC File Info Data... (6/6)" << endl;
	cout << "Number of Deuteron FileMCInfo Entries: " << MCdi_chain->GetEntries() << endl;
	// Loop over MC Deuteron FileMCInfo entries
	for (int i=0; i<MCdi_chain->GetEntries(); i++) {

		// Get entry
		MCdi_chain->GetEntry(i);

		// Get MC generation parameters
		double ngen = MCd_info->ngen_datacard;
		double rmin = MCd_comp->momentum[0];
		double rmax = MCd_comp->momentum[1];

		// 1/R generation spectrum
		TF1 *genFlux = new TF1("genFlux", "[0]/(x)", rmin, rmax);
		// Normalisation
		genFlux->SetParameter(0, log(rmax / rmin));

		// Loop over rigidity bins
	    for (int j=0; j<Bin_num; j++) {
	      // Fraction of spectrum in bin
	      double frac = genFlux->Integral(Bin_edges[j], Bin_edges[j+1]) / genFlux->Integral(rmin, rmax);
	      // Number of events in fraction
	      double evnum = frac * ngen;
	      // Set bin to current bin count + number of events
	      MCd_gen->SetBinContent(j+1, MCd_gen->GetBinContent(j+1) + evnum);
	    }
		
	}

	cout << "\n All done! :) \n"

}