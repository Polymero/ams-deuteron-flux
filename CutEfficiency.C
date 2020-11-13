// C++ macro for cut efficiency calculation
// Created        12-11-20
// Last Edited    12-11-20

// Include header file(s)
#include <iostream>
#include "Ntp.h"
#include "Simple.h"

void CutEfficiency() {

  // No canvas titles
  gStyle->SetOptTitle(0);

  // Import trees from Simple.root
  TChain simp_chain("Simp");
  simp_chain.Add("../Simp.root");

  // Cuts
  TCut* Crig = new TCut("trk_rig > 0 && trk_rig <= 22");
  TCut* Cpar = new TCut("status & 10 == 1");
  TCut* Ccon = new TCut("abs(tof_beta-rich_beta)/tof_beta < 0.05");
  TCut* Cbet = new TCut("tof_beta > 0");
  TCut* Cchi = new TCut("trk_chisqn[0] < 10 && trk_chisqn[1] < 10 && trk_chisqn[0] > 0 && trk_chisqn[1] > 0");
  TCut* Cinn = new TCut("trk_q_inn > 0.80 && trk_q_inn < 1.30");
  TCut* Clay = new TCut("trk_q_lay[4]>=0 && trk_q_lay[1]>=0 && trk_q_lay[2]>=0 && trk_q_lay[3]>=0 && trk_q_lay[5]>=0 && trk_q_lay[6]>=0 && trk_q_lay[7]>=0 && trk_q_lay[8]>=0 && trk_q_lay[0]>=0");
  TCut* Cgeo = new TCut("trk_rig > 1.2*cf");

  // Cut Map
  map <class TCut, string> CutMap = {

    {*Crig, "General"},
    {*Cpar, "Combi"},
    {*Ccon, "TOF"},
    {*Cbet, "TOF"},
    {*Cchi, "Tracker"},
    {*Cinn, "Tracker"},
    {*Clay, "Tracker"},
    {*Cgeo, "Tracker"}

  };

  Int_t i = 0;
  // Loop over cuts
  for (auto entry : CutMap) {

    // Get current cut
    TCut current_cut = entry.first;
    TCut base_cut = TCut("");

    // Loop over all cuts again
    for (auto i : CutMap) {

      // Add other cuts of other instuments
      if (i.second != entry.second) {
        base_cut += i.first;
      }

    }

    // Create canvas
    TCanvas* c1 = new TCanvas(Form("c%d", i), "Ratio Histogram");
    i++;

    // Draw histograms
    simp_chain.Draw(Form("trk_rig >> hist1%d(100, 0, 22)", i), base_cut, "");
    simp_chain.Draw(Form("trk_rig >> hist2%d(100, 0, 22)", i), base_cut && current_cut, "");

    // Get TH1F objects
    TH1F* hist1 = new TH1F();
    TH1F* hist2 = new TH1F();
    gROOT->GetObject(Form("hist1%d", i), hist1);
    gROOT->GetObject(Form("hist2%d", i), hist2);

    TH1F* hist3 = new TH1F();
    hist2->Divide(hist1);
    hist2->Draw();

    // Draw Canvas
    c1->Draw();


  }

}
