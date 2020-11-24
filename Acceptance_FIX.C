// C++ macro for acceptance calculation. PNG output in ./Acceptance
// Created        22-11-20
// Last Edited    23-11-20

// Include header file(s)
#include <iostream>
#include "Header Files/Ntp.h"

void Acceptance() {
  // Define Rigidity Bins
  double bin_left[32] = {1.00,1.16,1.33,1.51,1.71,1.92,2.15,2.40,2.67,2.97,3.29,3.64,4.02,4.43,4.88,5.37,5.90,6.47,7.09,7.76,8.48,9.26,10.1,11.0,12.0,13.0,14.1,15.3,16.6,18.0,19.5,21.1};
  double bin_right[32] = {1.16,1.33,1.51,1.71,1.92,2.15,2.40,2.67,2.97,3.29,3.64,4.02,4.43,4.88,5.37,5.90,6.47,7.09,7.76,8.48,9.26,10.1,11.0,12.0,13.0,14.1,15.3,16.6,18.0,19.5,21.1,22.8};
  double bin_edges[33] = {1.00,1.16,1.33,1.51,1.71,1.92,2.15,2.40,2.67,2.97,3.29,3.64,4.02,4.43,4.88,5.37,5.90,6.47,7.09,7.76,8.48,9.26,10.1,11.0,12.0,13.0,14.1,15.3,16.6,18.0,19.5,21.1,22.8};

  // Define Ngen Histogram
  TH1F* NgenHist = new TH1F("NgenHist", "Generated Events per Rigidity Bin", 32, bin_edges);
  double total_gen = 0;

  // Create canvas
  TCanvas* c1 = new TCanvas("Ngen", "Generated Events per Rigidity Bin");



  // Loop over MC root files to fill Ngen histogram
  int start_num = 1209496744;
  for (Int_t i=0; i<13; i++) {

    // Import tree
    TChain mc_chain("Compact");
    TChain fi_chain("File");
    mc_chain.Add(Form("../MC Protons/%d.root", start_num + i));
    fi_chain.Add(Form("../MC Protons/%d.root", start_num + i));

    // Create empty class objects
    NtpCompact *Compact = new class NtpCompact();
    FileMCInfo *FMCI = new class FileMCInfo();
    // Set Branch address
    mc_chain.SetBranchAddress("Compact", &Compact);
    fi_chain.SetBranchAddress("FileMCInfo", &FMCI);

    // Initialise parameters
    double ngen = 0;
    double rig_min = 0;
    double rig_max = 0;
    // Get FileMCInfo entry for information about MC generation
    fi_chain.GetEntry(0);
    ngen = FMCI->ngen_datacard;
    rig_min = FMCI->momentum[0];
    rig_max = FMCI->momentum[1];

    // Add ngen to total_gen for check
    total_gen += ngen;

    // Generated 1/R Spectrum
    TF1 *genFlux = new TF1("genFlux", "[0]/(x)", rig_min, rig_max);
    // Normalisation
    genFlux->SetParameter(0, log(rig_max / rig_min));

    // Loop over bins
    double check = 0;
    for (Int_t j=0; j<32; j++) {
      // Fraction of spectrum in bin
      double frac = genFlux->Integral(bin_left[j], bin_right[j])/genFlux->Integral(rig_min, rig_max);
      check += frac;
      // Number of events in fraction
      double Nev = frac * ngen;
      // Set bin to current count + Nev
      NgenHist->SetBinContent(j, NgenHist->GetBinContent(j) + Nev);

    }
    cout << "Total fraction in R=(1.0, 22.8) of " << start_num+i << ".root: " << check << endl;
  }

  // Sum over bins for check
  double total_bin = 0;
  for (Int_t i=0; i<32; i++) {
    total_bin += NgenHist->GetBinContent(i);
  }

  // Draw Histogram
  NgenHist->Draw();
  c1->SetLogx();
  c1->SetLogy();

  // Draw canvas
  c1->Draw();
  // Print canvas
  c1->Print("./Acceptance/Acceptance Histograms.png");

  cout << "Total Generated Events: " << total_gen << endl;
  cout << "Total Events in R=(1.0, 22.8): " << total_bin << endl;

}
