// C(++) macro for the trigger efficiency
// Created 28-11-20
// Last edited 28-11-20

// Include header file(s)
#include <iostream>
#include "Header Files/Ntp.h"

void TriggerEfficiency() {

  // Define Rigidity Bins
  double bin_left[32] = {1.00,1.16,1.33,1.51,1.71,1.92,2.15,2.40,2.67,2.97,3.29,3.64,4.02,4.43,4.88,5.37,5.90,6.47,7.09,7.76,8.48,9.26,10.1,11.0,12.0,13.0,14.1,15.3,16.6,18.0,19.5,21.1};
  double bin_right[32] = {1.16,1.33,1.51,1.71,1.92,2.15,2.40,2.67,2.97,3.29,3.64,4.02,4.43,4.88,5.37,5.90,6.47,7.09,7.76,8.48,9.26,10.1,11.0,12.0,13.0,14.1,15.3,16.6,18.0,19.5,21.1,22.8};
  double bin_middle[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double ErrRig[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  for (Int_t i=0; i<32; i++) {
    bin_middle[i] = (bin_right[i] + bin_left[i]) / 2;
    ErrRig[i] = (bin_right[i] - bin_left[i]) / 2;
  }
  double bin_edges[33] = {1.00,1.16,1.33,1.51,1.71,1.92,2.15,2.40,2.67,2.97,3.29,3.64,4.02,4.43,4.88,5.37,5.90,6.47,7.09,7.76,8.48,9.26,10.1,11.0,12.0,13.0,14.1,15.3,16.6,18.0,19.5,21.1,22.8};

  // Import tree
  TChain mc_chain("Compact");
  mc_chain.Add("../MC Protons/*.root");
  // Create empty class objects
  NtpCompact *Compact = new class NtpCompact();
  mc_chain.SetBranchAddress("Compact", &Compact);

  // Define histograms
  TH1F* RigHist = new TH1F("RigHist", "Entries per Rigidity Bin", 32, bin_edges);
  TH1F* PhysHist = new TH1F("PhysHist", "Physical Triggers per Rigidity Bin", 32, bin_edges);
  TH1F* UnphHist = new TH1F("UnphHist", "Unphysical Triggers per Rigidity Bin", 32, bin_edges);
  TH1F* TriggEffHist = new TH1F("TriggEffHist", "Trigger Efficiency per Rigidity Bin", 32, bin_edges);

  // Loop to fill histogram 1 and 2
  for (Int_t i = 0; i < mc_chain.GetEntries(); i++) {

    mc_chain.GetEntry(i);

    bool HasPhysTrig = ((Compact->sublvl1&0x3E)!=0)&&((Compact->trigpatt&0x2)!=0);
    bool HasUnphTrig = ((Compact->sublvl1&0x3E)==0)&&((Compact->trigpatt&0x2)!=0);
    if (HasPhysTrig) {
      PhysHist->Fill(Compact->trk_rig[0]);
    }
    if (HasUnphTrig) {
      UnphHist->Fill(Compact->trk_rig[0]);
    }
    RigHist->Fill(Compact->trk_rig[0]);

  }

  // Loop to fill histogram 3
  for (Int_t i = 0; i < 32; i++) {

    cout << RigHist->GetBinContent(i+1) << "     "
    << PhysHist->GetBinContent(i+1) << "     "
    << UnphHist->GetBinContent(i+1) << "     "
    << PhysHist->GetBinContent(i+1) + UnphHist->GetBinContent(i+1) << "     "
    << PhysHist->GetBinContent(i+1) / RigHist->GetBinContent(i+1) << endl;

    TriggEffHist->SetBinContent(i+1, PhysHist->GetBinContent(i+1) / RigHist->GetBinContent(i+1));

  }

  // Create canvas 1
  TCanvas* c1 = new TCanvas("Phys", "PHysical Triggers per Rigidity Bin");
  PhysHist->Draw();
  PhysHist->SetStats(0);
  c1->SetLogx();
  c1->Draw();

  // Create canvas 2
  TCanvas* c2 = new TCanvas("Unph", "Unphysical Triggers per Rigidity Bin");
  UnphHist->Draw();
  UnphHist->SetStats(0);
  c2->SetLogx();
  c2->Draw();

  // Create canvas 3
  TCanvas* c3 = new TCanvas("TrigEff", "Trigger Efficiency per Rigidity Bin");
  TriggEffHist->Draw();
  TriggEffHist->SetStats(0);
  c3->SetLogx();
  c3->Draw();

  // Convert TH1F to TGraph
  double TrigEff[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  for (Int_t i=0; i<32; i++) {
    TrigEff[i] = TriggEffHist->GetBinContent(i+1);
  }
  // Create TGraph
  TGraphErrors* TEGraph = new TGraphErrors(32, bin_middle, TrigEff, ErrRig, 0);
  TCanvas* c4 = new TCanvas("TrigEff2", "Trigger Efficiency TGraph");
  // Styling
  TEGraph->SetMarkerStyle(20);
  TEGraph->SetMarkerSize(1);
  TEGraph->SetMarkerColor(kBlue);
  // Axes
  TEGraph->GetXaxis()->SetTitle("R [GV]");
  TEGraph->GetYaxis()->SetTitle("Trigger Efficiency");
  TEGraph->SetMaximum(1);
  TEGraph->SetMinimum(0);
  c4->SetLogx();

  // Draw
  TEGraph->Draw("AP");

}
