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
  //TH1F* Htot = new TH1F("Htot", "Htot", 32, bin_edges);
  // Instrument bases
  TH1F* Btof = new TH1F("Btof", "Btof", 32, bin_edges);
  TH1F* Btrk = new TH1F("Btrk", "Btrk", 32, bin_edges);

  // Loop over all entries of the chain
  for (Int_t i = 0; i < simp_chain.GetEntries(); i++) {

    simp_chain.GetEntry(i);

    // Boolean cuts (instead of TCut)
    bool Crig = (Tool->trk_rig > 0)&&(Tool->trk_rig <= 22);
    bool Ctrg = ((Tool->sublvl1&0x3E)!=0)&&((Tool->trigpatt&0x2)!=0);
    bool Cpar = Tool->status % 10 == 1;
    bool Ccon = std::abs(Tool->tof_beta - Tool->rich_beta)/Tool->tof_beta < 0.05;
    bool Cbet = Tool->tof_beta > 0;
    bool Cchi = (Tool->trk_chisqn[0] < 10)&&(Tool->trk_chisqn[1] < 10)&&(Tool->trk_chisqn[0] > 0)&&(Tool->trk_chisqn[1] > 0);
    bool Cinn = (Tool->trk_q_inn > 0.80)&&(Tool->trk_q_inn < 1.30);
    bool Clay = (Tool->trk_q_lay[0] >= 0)&&(Tool->trk_q_lay[1] >= 0)&&(Tool->trk_q_lay[2] >= 0)&&(Tool->trk_q_lay[3] >= 0)&&(Tool->trk_q_lay[4] >= 0)&&(Tool->trk_q_lay[5] >= 0)&&(Tool->trk_q_lay[6] >= 0)&&(Tool->trk_q_lay[7] >= 0)&&(Tool->trk_q_lay[8] >= 0);
    bool Cgeo = Tool->trk_rig > 1.2 * Tool->cf;

    bool Ctof = Cpar && Ccon && Cbet;
    bool Ctrk = Cchi && Cinn && Clay && Cgeo;

    if (Ctrk && Ctrg && Crig) {Btof->Fill(Tool->trk_rig);}
    if (Ctof && Ctrg && Crig) {Btrk->Fill(Tool->trk_rig);}

    if (Cpar && Ctrk && Ctrg && Crig) {Hpar->Fill(Tool->trk_rig);}
    if (Ccon && Ctrk && Ctrg && Crig) {Hcon->Fill(Tool->trk_rig);}
    if (Cbet && Ctrk && Ctrg && Crig) {Hbet->Fill(Tool->trk_rig);}
    if (Cchi && Ctof && Ctrg && Crig) {Hchi->Fill(Tool->trk_rig);}
    if (Cinn && Ctof && Ctrg && Crig) {Hinn->Fill(Tool->trk_rig);}
    if (Clay && Ctof && Ctrg && Crig) {Hlay->Fill(Tool->trk_rig);}
    if (Cgeo && Ctof && Ctrg && Crig) {Hgeo->Fill(Tool->trk_rig);}

  }

  Hpar->Divide(Btof);
  Hcon->Divide(Btof);
  Hbet->Divide(Btof);
  Hchi->Divide(Btrk);
  Hinn->Divide(Btrk);
  Hlay->Divide(Btrk);
  Hgeo->Divide(Btrk);

  TH1F* TotCutEff = (TH1F*) Hpar->Clone();
  TotCutEff->Multiply(Hcon);
  TotCutEff->Multiply(Hbet);
  TotCutEff->Multiply(Hchi);
  TotCutEff->Multiply(Hinn);
  TotCutEff->Multiply(Hlay);
  TotCutEff->Multiply(Hgeo);

  TCanvas* c1 = new TCanvas("c1", "Cut Efficiencies");
  c1->SetLogx();

  Hpar->Draw(); Hpar->SetAxisRange(0,1,"Y"); c1->Print("./CutEfficiency/Single Particle Cut.png");
  Hcon->Draw(); Hcon->SetAxisRange(0,1,"Y"); c1->Print("./CutEfficiency/Beta Consistency Cut.png");
  Hbet->Draw(); Hbet->SetAxisRange(0,1,"Y"); c1->Print("./CutEfficiency/TOF Beta Cut.png");
  Hchi->Draw(); Hchi->SetAxisRange(0,1,"Y"); c1->Print("./CutEfficiency/Well Consructed Track Cut.png");
  Hinn->Draw(); Hinn->SetAxisRange(0,1,"Y"); c1->Print("./CutEfficiency/Inner Tracker Charge Cut.png");
  Hlay->Draw(); Hlay->SetAxisRange(0,1,"Y"); c1->Print("./CutEfficiency/Tracker Layer Charge Cut.png");
  Hgeo->Draw(); Hgeo->SetAxisRange(0,1,"Y"); c1->Print("./CutEfficiency/Geo-magnetic Cut-off Cut.png");

  TotCutEff->Draw();
  TotCutEff->SetAxisRange(0,.2,"Y");

  c1->Print("Total Cut Efficiency.png");

}
