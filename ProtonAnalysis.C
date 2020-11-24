// C++ class for full PROTON analysis
// Created        23-11-20
// Last Edited    23-11-20

// Include header file(s)
#include <iostream>
#include <string>
#include "Header Files/Ntp.h"
#include "TChain.h"
#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"

//------------------------------------------------------------------------------
// CLASS DEFINITION
//------------------------------------------------------------------------------

// Class
class Anaaqra {
  public: // Acces Specifier

    //--------------------------------------------------------------------------
    // ATTRIBUTES
    //--------------------------------------------------------------------------
    // File paths
    string run_files = "../Simp.root";
    string pmc_dir = "../MC Protons";
    // Rigidity bins
    int bin_num = 32;
    double bin_edges[33] = {1.00,1.16,1.33,1.51,1.71,1.92,2.15,2.40,2.67,2.97,
      3.29,3.64,4.02,4.43,4.88,5.37,5.90,6.47,7.09,7.76,8.48,9.26,10.1,11.0,
      12.0,13.0,14.1,15.3,16.6,18.0,19.5,21.1,22.8};

    //--------------------------------------------------------------------------
    // CONSTRUCTORS
    //--------------------------------------------------------------------------
    Anaaqra() { // Default constructor
      cout << "Class succesfully constructed!" << endl;
    }

    //--------------------------------------------------------------------------
    // LIST OF METHODS
    //--------------------------------------------------------------------------
    TH1F Acceptance(); // Returns acceptance as function of rigidity

};

//------------------------------------------------------------------------------
// METHOD FUNCTIONS
//------------------------------------------------------------------------------
// Returns TH1F of Acceptance (for flux calculation)
TH1F Anaaqra::Acceptance() {

  // Define Ngen histogram
  TH1F* NgenHist = new TH1F("NgenHist", "Generated MC Events per Rigidity Bin", 32, bin_edges);
  double total_gen = 0;
  // Create canvas
  TCanvas* cAcc = new TCanvas("cAcc", "Generated MC Events per Rigitidy Bin");

  // Loop over MC root files to fill Ngen histogram
  int start_num = 1209496744;
  for (int i=0; i<13; i++) { // !: change 13 to adaptive int

    // Import tree
    TChain mc_chain("Compact");
    TChain fi_chain("File");
    mc_chain.Add(Form("../MC Protons/%d.root", start_num + i));
    fi_chain.Add(Form("../MC Protons/%d.root", start_num + i));
    // Create empty class objects
    NtpCompact *MCompact = new class NtpCompact();
    FileMCInfo *FMCI = new class FileMCInfo();
    // Set branch address
    mc_chain.SetBranchAddress("Compact", &MCompact);
    fi_chain.SetBranchAddress("FileMCInfo", &FMCI);

    // MC parameters
    double ngen = 0;
    double rig_min = 0;
    double rig_max = 0;
    // Get FileMCInfo entry for MC information
    fi_chain.GetEntry(0);
    ngen = FMCI->ngen_datacard;
    rig_min = FMCI->momentum[0];
    rig_max = FMCI->momentum[1];
    // Add ngen to total_gen for visual check
    total_gen += ngen;

    // Generate 1/R spectrum
    TF1 *genFlux = new TF1("genFlux", "[0]/(x)", rig_min, rig_max);
    // Normalisation
    genFlux->SetParameter(0, log(rig_max / rig_min));

    // Loop over bins
    double total_frac = 0;
    for (int j=0; j<bin_num; j++) {
      // Fraction of spectrum in bin
      double frac = genFlux->Integral(bin_edges[j], bin_edges[j])/genFlux->Integral(rig_min, rig_max);
      total_frac += frac;
      // Number of events in fraction
      double nev = frac * ngen;
      // Set bin to current count + nev
      NgenHist->SetBinContent(j, NgenHist->GetBinContent(j) + nev);

    }

    // Print total fractions for visual check
    cout << "Total fraction in R=(" << bin_edges[0] << ", " << bin_edges[bin_num+1]
         << ") of " << start_num+i << ".root: " << total_frac << endl;

  }

  // Sum over bins for visual check
  double total_bin = 0;
  for (int i=0; i<bin_num; i++) {
    total_bin += NgenHist->GetBinContent(i);
  }
  // Print sum for visual check
  cout << "Total generated MC events: " << total_gen << endl;
  cout << "Total MC events in R=" << bin_edges[0] << ", " << bin_edges[bin_num+1]
       << "): " << total_bin << endl;

  return *NgenHist;
}
