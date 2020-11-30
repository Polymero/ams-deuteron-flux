// C(++) macro for cut efficiency histograms
// Created        28-11-20
// Last edited    28-11-20

// Include header file(s)
#include <iostream>
#include "Header Files/Ntp.h"
#include "Header Files/Simple.h"

void CutEffHist() {

  double bin_edges[33] = {1.00,1.16,1.33,1.51,1.71,1.92,2.15,2.40,2.67,2.97,3.29,3.64,4.02,4.43,4.88,5.37,5.90,6.47,7.09,7.76,8.48,9.26,10.1,11.0,12.0,13.0,14.1,15.3,16.6,18.0,19.5,21.1,22.8};

  // Import trees
  TChain simp_chain("Simp");
  simp_chain.Add("../Simp.root");
  // Create empty class objects
  Miiqtool *Tool = new Miiqtool();
  simp_chain.SetBranchAddress("Simp", &Tool);

  // Histograms
  // Cuts
  //TH1F* Hrig = new TH1F("Hrig", "Hrig", 32, bin_edges);
  TH1F* Hpar = new TH1F("Hpar", "Hpar", 32, bin_edges);
  TH1F* Hcon = new TH1F("Hcon", "Hcon", 32, bin_edges);
  TH1F* Hbet = new TH1F("Hbet", "Hbet", 32, bin_edges);
  TH1F* Hchi = new TH1F("Hchi", "Hchi", 32, bin_edges);
  TH1F* Hinn = new TH1F("Hinn", "Hinn", 32, bin_edges);
  TH1F* Hlay = new TH1F("Hlay", "Hlay", 32, bin_edges);
  TH1F* Hgeo = new TH1F("Hgeo", "Hgeo", 32, bin_edges);
  // Total efficiency
  TH1F* Htot = new TH1F("Htot", "Htot", 32, bin_edges);
  // Instrument bases
  TH1F* Btof = new TH1F("Btof", "Btof", 32, bin_edges);
  TH1F* Btrk = new TH1F("Btrk", "Btrk", 32, bin_edges);

  // Loop over all entries of the chain
  for (Int_t i = 0; i < simp_chain->GetEntries(); i++) {

    simp_chain.GetEntry(i);

    // Boolean cuts (instead of TCut)
    bool Crig = (Tool->trk_rig > 0)&&(Tool->trk_rig <= 22);
    bool Cpar = Tool->status % 10 == 1;
    bool Ccon = abs(Tool->tof_beta - Tool->rich_beta)/Tool->tof_beta < 0.05;
    bool Cbet = Tool->tof_beta > 0;
    bool Cchi = (Tool->trk_chisqn[0] < 10)&&(Tool->trk_chisqn[1] < 10)&&(Tool->trk_chisqn[0] > 0)&&(Tool->trk_chisqn[1] > 0);
    bool Cinn = (Tool->trk_q_inn > 0.80)&&(Tool->trk_q_inn < 1.30);
    bool Clay = (Tool->trk_q_lay[0] >= 0)&&(Tool->trk_q_lay[1] >= 0)&&(Tool->trk_q_lay[2] >= 0)&&(Tool->trk_q_lay[3] >= 0)&&(Tool->trk_q_lay[4] >= 0)&&(Tool->trk_q_lay[5] >= 0)&&(Tool->trk_q_lay[6] >= 0)&&(Tool->trk_q_lay[7] >= 0)&&(Tool->trk_q_lay[8] >= 0);
    bool Cgeo = Tool->trk_rig > 1.2 * Tool->cf;
    //bool Ctrg = ((Tool->sublvl1&0x3E)!=0)&&((Tool->trigpatt&0x2)!=0);

    if () {
      Hrig->Fill(Tool->trk_rig);
    }

  }

  cout << endl;

}
