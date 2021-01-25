// C++ class for full PROTON analysis
// Created        23-11-20
// Last Edited    24-01-21

// Include header file(s)
#include <iostream>
#include <string>
#include <algorithm>
#include "Header Files/Ntp.h"
#include "Header Files/Simple.h"
#include "TChain.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TString.h"

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// CLASS DEFINITION
//------------------------------------------------------------------------------------------------------------------------------------------------------------

// Class
class Anaaqra {
  public: // Acces Specifier

    //----------------------------------------------------------------------------------------------------------------------------------------------------
    // ATTRIBUTES
    //----------------------------------------------------------------------------------------------------------------------------------------------------
    // Analysis type
    int atype = 0;
    // Rigidity bins
    const int Bin_num = 32;
    double Bin_edges[32+1] = {1.00,1.16,1.33,1.51,1.71,1.92,2.15,2.40,2.67,2.97,
                              3.29,3.64,4.02,4.43,4.88,5.37,5.90,6.47,7.09,7.76,
                              8.48,9.26,10.1,11.0,12.0,13.0,14.1,15.3,16.6,18.0,
                              19.5,21.1,22.8};
    double Bin_err[32];
    double Bin_mid[32];
    // MC Files
    int Start_num = 1209496744;

    // Histograms
    // RigBinner()
    TH1F *Events_raw      = new TH1F("Events_raw", "Raw AMS-02 Events", 32, Bin_edges);
    TH1F *Events_cut      = new TH1F("Events_cut", "Selected AMS-02 Events", 32, Bin_edges);
    // Exposure()
    TH1F *ExposureTime    = new TH1F("ExposureTime", "Exposure Time per Rigidity Bin", 32, Bin_edges);
    // ProtonRate()
    TH1F *RateHist        = new TH1F("RateHist", "Proton Rate per Rigidity Bin", 32, Bin_edges);
    // Acceptance()
    TH1F *MC_generated    = new TH1F("MC_generated", "Generated MC Events per Rigidity Bin", 32, Bin_edges);
    TH1F *MC_detected     = new TH1F("MC_detected", "Detected MC Events per Rigidity Bin", 32, Bin_edges);
    TH1F *AcceptHist      = new TH1F("AcceptHist", "Acceptance per Rigidity Bin", 32, Bin_edges);
    // CutEff()
    TH1F *CutEff_data     = new TH1F("CutEff_data", "Selection Efficiency of Data per Rigidity Bin", 32, Bin_edges);
    TH1F *CutEff_MC       = new TH1F("CutEff_MC", "Selection Efficiency of MC per Rigidity Bin", 32, Bin_edges);
    // TrigEff()
    TH1F *PhysHist_mc     = new TH1F("PhysHist_mc", "PhysHist_mc", 32, Bin_edges);
    TH1F *UnphHist_mc     = new TH1F("UnphHist_mc", "UnphHist_mc", 32, Bin_edges);
    TH1F *PhysHist_data   = new TH1F("PhysHist_data", "PhysHist_data", 32, Bin_edges);
    TH1F *UnphHist_data   = new TH1F("UnphHist_data", "UnphHist_data", 32, Bin_edges);
    TH1F *TrigEff_data    = new TH1F("TrigEff_data", "Trigger Efficiency of Data per Rigidity Bin", 32, Bin_edges);
    TH1F *TrigEff_MC      = new TH1F("TrigEff_MC", "Trigger Efficiency of MC per Rigidity Bin", 32, Bin_edges);
    // ProtonFlux()
    TH1F *FluxHist        = new TH1F("FluxHist", "Flux per Rigidity Bin", 32, Bin_edges);

    // Data objects
    TChain *Simp_chain    = new TChain("Simp");
    TChain *RTII_chain    = new TChain("RTIInfo");
    TChain *MC_chain      = new TChain("Compact");
    Miiqtool *Tool        = new Miiqtool();
    MiiqRTI *Woi          = new MiiqRTI();
    NtpCompact *MC_comp   = new NtpCompact();

    // TGraph arrays
    // RigBinner()
    double Ev_raw[32]; double Ev_raw_err[32];
    double Ev_cut[32]; double Ev_cut_err[32];
    // Exposure()
    double ExpTime[32];
    // Acceptance()
    double Accept[32];
    // CutEff()
    double CE_mc[32]; double CE_mc_err[32];
    double CE_data[32]; double CE_data_err[32];
    // TrigEff()
    double TE_mc[32]; double TE_mc_err[32];
    double TE_data[32]; double TE_data_err[32];
    // Rate()
    double PRate[32]; double PRate_err[32];
    // Flux()
    double PFlux[32]; double PFlux_err[32];
    double PSFlux[32]; double PSFlux_err[32];
    double PSFlux_SSDC_Ratio[32];

    //----------------------------------------------------------------------------------------------------------------------------------------------------
    // CONSTRUCTORS
    //----------------------------------------------------------------------------------------------------------------------------------------------------
    Anaaqra(string type) { // Default constructor

      // Analysis type
      if (type == "proton" || type == "p" || type == "1H"){
        atype = 1;
      } else if (type == "deuteron" || type == "d" || type == "2H"){
        atype = 2;
      } else {
        cout << "Analysis type not recognised!\n" << endl;
      }

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
      Simp_chain->Add("../Simp.root");
      RTII_chain->Add("../Simp.root");
      MC_chain->Add("../MC Protons/*.root");
      // Set branch addresses
      Simp_chain->SetBranchAddress("Simp", &Tool);
      RTII_chain->SetBranchAddress("RTI", &Woi);
      MC_chain->SetBranchAddress("Compact", &MC_comp);

      cout << "\nAnalysis type: " << type << endl;
      cout << "Class succesfully constructed!\n" << endl;

    };

    //----------------------------------------------------------------------------------------------------------------------------------------------------
    // LIST OF METHODS
    //----------------------------------------------------------------------------------------------------------------------------------------------------
    // Support functions
    bool EventSelectorCompact(NtpCompact* comp, const char* cutbit);
    bool EventSelectorSimple(Miiqtool* tool, const char* cutbit);
    // Singular (independent)
    void ParameterAnalysis(const char* cutbit = "1111111_111"); // Returns raw and cut histograms of various parameters
    void RigBinner();                           // Returns number of selected events as function of rigidity
    void Exposure();                            // Returns exposure time (livetime) as function of rigidity
    void Acceptance(bool apply_cuts = 0);       // Returns (geometric) acceptance as function of rigidity
    void CutEff(bool plot_all = 0);             // Returns the cut (selection) efficiency as function of rigidity
    void TrigEff(int delta = 100);              // Returns the trigger efficiency as function of rigidity
    void AerogelSlice();
    // Plural (dependent)
    void Rate();                                // Returns the proton rate as function of rigidity
    void Flux(bool comp = 0);                   // Returns the proton flux as function of rigidity

};

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// SUPPORT FUNCTIONS
//------------------------------------------------------------------------------------------------------------------------------------------------------------
// Returns bool if event passed specified selection of cuts
bool Anaaqra::EventSelectorCompact(NtpCompact* comp, const char* cutbit) {

  bool pass = 1;

  // Boolean cuts (instead of TCut)
  // Protons (base cuts)
  bool Crig = (comp->trk_rig[0] > Bin_edges[0])&&(comp->trk_rig[0] <= Bin_edges[Bin_num]);
  bool Ctrg = ((comp->sublvl1&0x3E)!=0)&&((comp->trigpatt&0x2)!=0);
  bool Cpar = comp->status % 10 == 1;
  bool Cbet = comp->tof_beta > 0.3;
  bool Cchi = (comp->trk_chisqn[0][0] < 10)&&(comp->trk_chisqn[0][1] < 10)&&(comp->trk_chisqn[0][0] > 0)&&(comp->trk_chisqn[0][1] > 0);
  bool Cinn = (comp->trk_q_inn > 0.80)&&(comp->trk_q_inn < 1.30);
  // Adjust return bool according to cutbit
  if (cutbit[0] == '1') {pass &= Crig;}
  if (cutbit[1] == '1') {pass &= Ctrg;}
  if (cutbit[2] == '1') {pass &= Cpar;}
  if (cutbit[3] == '1') {pass &= Cbet;}
  if (cutbit[4] == '1') {pass &= Cchi;}
  if (cutbit[5] == '1') {pass &= Cinn;}
  // Deuterons (additional cuts)
  if (atype == 2){
    bool Csel = comp->rich_select == 2;
    bool Ccon = std::abs(comp->tof_beta - comp->rich_beta)/comp->tof_beta < 0.05;
    bool Clay = comp->trk_q_lay[0] <= 1.7;
    // Adjust return bool according to cutbit
    if (cutbit[8] == '1') {pass &= Csel;}
    if (cutbit[9] == '1') {pass &= Ccon;}
    if (cutbit[10] == '1') {pass &= Clay;}
  }

  // Return
  return pass;

}



// Returns bool if event passed specified selection of cuts
bool Anaaqra::EventSelectorSimple(Miiqtool* tool, const char* cutbit) {

  bool pass = 1;

  // Boolean cuts (instead of TCut)
  // Protons (base cuts)
  bool Crig = (tool->trk_rig > Bin_edges[0])&&(tool->trk_rig <= Bin_edges[Bin_num]);
  bool Ctrg = ((tool->sublvl1&0x3E)!=0)&&((tool->trigpatt&0x2)!=0);
  bool Cpar = tool->status % 10 == 1;
  bool Cbet = tool->tof_beta > 0.3;
  bool Cchi = (tool->trk_chisqn[0] < 10)&&(tool->trk_chisqn[1] < 10)&&(tool->trk_chisqn[0] > 0)&&(tool->trk_chisqn[1] > 0);
  bool Cinn = (tool->trk_q_inn > 0.80)&&(tool->trk_q_inn < 1.30);
  bool Cgeo = tool->trk_rig > 1.2 * tool->cf;
  // Adjust return bool according to cutbit
  if (cutbit[0] == '1') {pass &= Crig;}
  if (cutbit[1] == '1') {pass &= Ctrg;}
  if (cutbit[2] == '1') {pass &= Cpar;}
  if (cutbit[3] == '1') {pass &= Cbet;}
  if (cutbit[4] == '1') {pass &= Cchi;}
  if (cutbit[5] == '1') {pass &= Cinn;}
  if (cutbit[6] == '1') {pass &= Cgeo;}
  // Deuterons (additional cuts)
  if (atype == 2){
    bool Csel = tool->rich_select == 2;
    bool Ccon = std::abs(tool->tof_beta - tool->rich_beta)/tool->tof_beta < 0.05;
    bool Clay = tool->trk_q_lay[0] <= 1.7;
    // Adjust return bool according to cutbit
    if (cutbit[8] == '1') {pass &= Csel;}
    if (cutbit[9] == '1') {pass &= Ccon;}
    if (cutbit[10] == '1') {pass &= Clay;}
  }

  // Return
  return pass;

}



//------------------------------------------------------------------------------------------------------------------------------------------------------------
// METHOD FUNCTIONS
//------------------------------------------------------------------------------------------------------------------------------------------------------------
// Returns TH1F of various parameters
void Anaaqra::ParameterAnalysis(const char* cutbit = "1111111_111") {

  cout << "Running ParameterAnalysis()..." << endl;

  // Create canvasses
  TCanvas *Beta_TOF_MC     = new TCanvas("Beta_TOF_MC", "TOF Beta MC");
  TCanvas *Beta_RICH_MC    = new TCanvas("Beta_RICH_MC", "RICH Beta MC");
  TCanvas *Rig_MC          = new TCanvas("Rig_MC", "Tracker Rigidity MC");
  TCanvas *Chisq_MC        = new TCanvas("ChiSq_MC", "Chi-squared Bending plane MC");
  TCanvas *Q_inn_MC        = new TCanvas("Q_inn_MC", "Inner Tracker Charge MC");
  TCanvas *Q_lay_MC        = new TCanvas("Q_lay_MC", "First Tracker Layer Charge MC");
  TCanvas *Mass_RICH_MC    = new TCanvas("Mass_RICH_MC", "RICH Mass MC");
  TCanvas *Beta_TOF_data   = new TCanvas("Beta_TOF_data", "TOF Beta Data");
  TCanvas *Beta_RICH_data  = new TCanvas("Beta_RICH_data", "RICH Beta Data");
  TCanvas *Rig_data        = new TCanvas("Rig_data", "Tracker Rigidity Data");
  TCanvas *Chisq_data      = new TCanvas("ChiSq_data", "Chi-squared Bending plane Data");
  TCanvas *Q_inn_data      = new TCanvas("Q_inn_data", "Inner Tracker Charge Data");
  TCanvas *Q_lay_data      = new TCanvas("Q_lay_data", "First Tracker Layer Charge Data");
  TCanvas *Mass_RICH_data  = new TCanvas("Mass_RICH_data", "RICH Mass Data");
  TCanvas *UTime_data      = new TCanvas("UTime_data", "UTime Data");

  // Temporary histograms
  TH1F *tofbeta_raw = new TH1F("t1r", "", 100, 0.01, 2.0);
  TH1F *tofbeta_cut = new TH1F("t1c", "", 100, 0.01, 2.0);
  TH1F *richbeta_raw = new TH1F("t2r", "", 100, 0.74, 1.05);
  TH1F *richbeta_cut = new TH1F("t2c", "", 100, 0.74, 1.05);
  TH1F *rig_raw = new TH1F("t3r", "", 100, 0, 22);
  TH1F *rig_cut = new TH1F("t3c", "", 100, 0, 22);
  TH1F *chisq_raw = new TH1F("t4r", "", 100, -0.5, 10.5);
  TH1F *chisq_cut = new TH1F("t4c", "", 100, -0.5, 10.5);
  TH1F *qinn_raw = new TH1F("t5r", "", 100, 0.3, 3.0);
  TH1F *qinn_cut = new TH1F("t5c", "", 100, 0.3, 3.0);
  TH1F *qlay_raw = new TH1F("t6r", "", 100, 0.3, 3.0);
  TH1F *qlay_cut = new TH1F("t6c", "", 100, 0.3, 3.0);
  TH1F *mass_raw = new TH1F("t7r", "", 100, 0, 2.5);
  TH1F *mass_cut = new TH1F("t7c", "", 100, 0, 2.5);
  TH1F *utime_raw = new TH1F("t8r", "", 100, 1.30580e9, 1.30660e9);
  TH1F *utime_cut = new TH1F("t8c", "", 100, 1.30580e9, 1.30660e9);

  // Loop over data entries
  for (int i=0; i<Simp_chain->GetEntries(); i++) {

    // Get entry
    Simp_chain->GetEntry(i);

    // Fill raw histograms
    tofbeta_raw->Fill(Tool->tof_beta);
    richbeta_raw->Fill(Tool->rich_beta);
    rig_raw->Fill(Tool->trk_rig);
    chisq_raw->Fill(Tool->trk_chisqn[1]);
    qinn_raw->Fill(Tool->trk_q_inn);
    qlay_raw->Fill(Tool->trk_q_lay[0]);
    mass_raw->Fill(Tool->trk_q_inn * Tool->trk_rig * TMath::Sqrt(1 / Tool->rich_beta / Tool->rich_beta - 1));
    utime_raw->Fill(Tool->utime);

    // Fill cut histograms if selected
    if (EventSelectorSimple(Tool, cutbit)) {
      tofbeta_cut->Fill(Tool->tof_beta);
      richbeta_cut->Fill(Tool->rich_beta);
      rig_cut->Fill(Tool->trk_rig);
      chisq_cut->Fill(Tool->trk_chisqn[1]);
      qinn_cut->Fill(Tool->trk_q_inn);
      qlay_cut->Fill(Tool->trk_q_lay[0]);
      mass_cut->Fill(Tool->trk_q_inn * Tool->trk_rig * TMath::Sqrt(1 / Tool->rich_beta / Tool->rich_beta - 1));
      utime_cut->Fill(Tool->utime);
    }

  }

  // Create canvasses
  // TOF Beta Data
  Beta_TOF_data->cd();
  tofbeta_raw->Draw(); tofbeta_cut->Draw("same");
  tofbeta_raw->GetXaxis()->SetTitle("#beta (v/c)"); tofbeta_raw->GetYaxis()->SetTitle("Events");
  Beta_TOF_data->SetLogx(0); Beta_TOF_data->SetLogy(1); tofbeta_raw->SetMinimum(1);
  tofbeta_raw->SetLineColor(kRed); tofbeta_cut->SetLineColor(kBlue);
  Beta_TOF_data->Draw(); Beta_TOF_data->Print("./ProtonAnalysis/PA/TOF Beta Data.png");
  // RICH Beta Data
  Beta_RICH_data->cd();
  richbeta_raw->Draw(); richbeta_cut->Draw("same");
  richbeta_raw->GetXaxis()->SetTitle("#beta (v/c)"); richbeta_raw->GetYaxis()->SetTitle("Events");
  Beta_RICH_data->SetLogx(0); Beta_RICH_data->SetLogy(1); richbeta_raw->SetMinimum(1);
  richbeta_raw->SetLineColor(kRed); richbeta_cut->SetLineColor(kBlue);
  Beta_RICH_data->Draw(); Beta_RICH_data->Print("./ProtonAnalysis/PA/RICH Beta Data.png");
  // Rigidity Data
  Rig_data->cd();
  rig_raw->Draw(); rig_cut->Draw("same");
  rig_raw->GetXaxis()->SetTitle("R [GV]"); rig_raw->GetYaxis()->SetTitle("Events");
  Rig_data->SetLogx(0); Rig_data->SetLogy(1); rig_raw->SetMinimum(1);
  rig_raw->SetLineColor(kRed); rig_cut->SetLineColor(kBlue);
  Rig_data->Draw(); Rig_data->Print("./ProtonAnalysis/PA/Rigidity Data.png");
  // Chi-squared Bending Plane Data
  Chisq_data->cd();
  chisq_raw->Draw(); chisq_cut->Draw("same");
  chisq_raw->GetXaxis()->SetTitle("Chi-squared"); chisq_raw->GetYaxis()->SetTitle("Events");
  Chisq_data->SetLogx(0); Chisq_data->SetLogy(1); chisq_raw->SetMinimum(1);
  chisq_raw->SetLineColor(kRed); chisq_cut->SetLineColor(kBlue);
  Chisq_data->Draw(); Chisq_data->Print("./ProtonAnalysis/PA/Chi-squared Bending Plane Data.png");
  // Inner Tracker Charge Data
  Q_inn_data->cd();
  qinn_raw->Draw(); qinn_cut->Draw("same");
  qinn_raw->GetXaxis()->SetTitle("Z"); qinn_raw->GetYaxis()->SetTitle("Events");
  Q_inn_data->SetLogx(0); Q_inn_data->SetLogy(1); qinn_raw->SetMinimum(1);
  qinn_raw->SetLineColor(kRed); qinn_cut->SetLineColor(kBlue);
  Q_inn_data->Draw(); Q_inn_data->Print("./ProtonAnalysis/PA/Inner Tracker Charge Data.png");
  // First Tracker Layer Charge Data
  Q_lay_data->cd();
  qlay_raw->Draw(); qlay_cut->Draw("same");
  qlay_raw->GetXaxis()->SetTitle("Z"); qlay_raw->GetYaxis()->SetTitle("Events");
  Q_lay_data->SetLogx(0); Q_lay_data->SetLogy(1); qlay_raw->SetMinimum(1);
  qlay_raw->SetLineColor(kRed); qlay_cut->SetLineColor(kBlue);
  Q_lay_data->Draw(); Q_lay_data->Print("./ProtonAnalysis/PA/First Tracker Layer Charge Data.png");
  // RICH Mass Data
  Mass_RICH_data->cd();
  mass_raw->Draw(); mass_cut->Draw("same");
  mass_raw->GetXaxis()->SetTitle("m [GeV/c#^2]"); mass_raw->GetYaxis()->SetTitle("Events");
  Mass_RICH_data->SetLogx(0); Mass_RICH_data->SetLogy(1); mass_raw->SetMinimum(1);
  mass_raw->SetLineColor(kRed); mass_cut->SetLineColor(kBlue);
  Mass_RICH_data->Draw(); Mass_RICH_data->Print("./ProtonAnalysis/PA/RICH Mass Data.png");
  // UTime Data
  UTime_data->cd();
  utime_raw->Draw(); utime_cut->Draw("same");
  utime_raw->GetXaxis()->SetTitle("utime"); utime_raw->GetYaxis()->SetTitle("Events");
  UTime_data->SetLogx(0); UTime_data->SetLogy(1); utime_raw->SetMinimum(1);
  utime_raw->SetLineColor(kRed); utime_cut->SetLineColor(kBlue);
  utime_raw->GetXaxis()->SetTimeDisplay(1); utime_raw->GetXaxis()->SetTimeFormat("%m/%d");
  //utime_raw->GetXaxis()->SetTimeOffset(-788918400);
  UTime_data->Draw(); UTime_data->Print("./ProtonAnalysis/PA/UTime Data.png");

  // Clear temporary histograms
  tofbeta_raw->Reset("ICESM"); tofbeta_cut->Reset("ICESM");
  richbeta_raw->Reset("ICESM"); richbeta_cut->Reset("ICESM");
  rig_raw->Reset("ICESM"); rig_cut->Reset("ICESM");
  chisq_raw->Reset("ICESM"); chisq_cut->Reset("ICESM");
  qinn_raw->Reset("ICESM"); qinn_cut->Reset("ICESM");
  qlay_raw->Reset("ICESM"); qlay_cut->Reset("ICESM");
  mass_raw->Reset("ICESM"); mass_cut->Reset("ICESM");

  // Loop over MC entries
  for (int i=0; i<MC_chain->GetEntries(); i++) {

    // Get entry
    MC_chain->GetEntry(i);

    // Fill raw histograms
    tofbeta_raw->Fill(MC_comp->tof_beta);
    richbeta_raw->Fill(MC_comp->rich_beta);
    rig_raw->Fill(MC_comp->trk_rig[0]);
    chisq_raw->Fill(MC_comp->trk_chisqn[0][1]);
    qinn_raw->Fill(MC_comp->trk_q_inn);
    qlay_raw->Fill(MC_comp->trk_q_lay[0]);
    mass_raw->Fill(MC_comp->trk_q_inn * MC_comp->trk_rig[0] * TMath::Sqrt(1 / MC_comp->rich_beta / MC_comp->rich_beta - 1));

    // Fill cut histograms if selected
    if (EventSelectorCompact(MC_comp, cutbit)) {
      tofbeta_cut->Fill(MC_comp->tof_beta);
      richbeta_cut->Fill(MC_comp->rich_beta);
      rig_cut->Fill(MC_comp->trk_rig[0]);
      chisq_cut->Fill(MC_comp->trk_chisqn[0][1]);
      qinn_cut->Fill(MC_comp->trk_q_inn);
      qlay_cut->Fill(MC_comp->trk_q_lay[0]);
      mass_cut->Fill(MC_comp->trk_q_inn * MC_comp->trk_rig[0] * TMath::Sqrt(1 / MC_comp->rich_beta / MC_comp->rich_beta - 1));

    }

  }

  // Create canvasses
  // TOF Beta MC
  Beta_TOF_MC->cd();
  tofbeta_raw->Draw(); tofbeta_cut->Draw("same");
  tofbeta_raw->GetXaxis()->SetTitle("#beta (v/c)"); tofbeta_raw->GetYaxis()->SetTitle("Events");
  Beta_TOF_MC->SetLogx(0); Beta_TOF_MC->SetLogy(1); tofbeta_raw->SetMinimum(1);
  tofbeta_raw->SetLineColor(kRed); tofbeta_cut->SetLineColor(kBlue);
  Beta_TOF_MC->Draw(); Beta_TOF_MC->Print("./ProtonAnalysis/PA/TOF Beta MC.png");
  // RICH Beta MC
  Beta_RICH_MC->cd();
  richbeta_raw->Draw(); richbeta_cut->Draw("same");
  richbeta_raw->GetXaxis()->SetTitle("#beta (v/c)"); richbeta_raw->GetYaxis()->SetTitle("Events");
  Beta_RICH_MC->SetLogx(0); Beta_RICH_MC->SetLogy(1); richbeta_raw->SetMinimum(1);
  richbeta_raw->SetLineColor(kRed); richbeta_cut->SetLineColor(kBlue);
  Beta_RICH_MC->Draw(); Beta_RICH_MC->Print("./ProtonAnalysis/PA/RICH Beta MC.png");
  // Rigidity MC
  Rig_MC->cd();
  rig_raw->Draw(); rig_cut->Draw("same");
  rig_raw->GetXaxis()->SetTitle("R [GV]"); rig_raw->GetYaxis()->SetTitle("Events");
  Rig_MC->SetLogx(0); Rig_MC->SetLogy(1); rig_raw->SetMinimum(1);
  rig_raw->SetLineColor(kRed); rig_cut->SetLineColor(kBlue);
  Rig_MC->Draw(); Rig_MC->Print("./ProtonAnalysis/PA/Rigidity MC.png");
  // Chi-squared Bending Plane MC
  Chisq_MC->cd();
  chisq_raw->Draw(); chisq_cut->Draw("same");
  chisq_raw->GetXaxis()->SetTitle("Chi-squared"); chisq_raw->GetYaxis()->SetTitle("Events");
  Chisq_MC->SetLogx(0); Chisq_MC->SetLogy(1); chisq_raw->SetMinimum(1);
  chisq_raw->SetLineColor(kRed); chisq_cut->SetLineColor(kBlue);
  Chisq_MC->Draw(); Chisq_MC->Print("./ProtonAnalysis/PA/Chi-squared Bending Plane MC.png");
  // Inner Tracker Charge MC
  Q_inn_MC->cd();
  qinn_raw->Draw(); qinn_cut->Draw("same");
  qinn_raw->GetXaxis()->SetTitle("Z"); qinn_raw->GetYaxis()->SetTitle("Events");
  Q_inn_MC->SetLogx(0); Q_inn_MC->SetLogy(1); qinn_raw->SetMinimum(1);
  qinn_raw->SetLineColor(kRed); qinn_cut->SetLineColor(kBlue);
  Q_inn_MC->Draw(); Q_inn_MC->Print("./ProtonAnalysis/PA/Inner Tracker Charge MC.png");
  // First Tracker Layer Charge MC
  Q_lay_MC->cd();
  qlay_raw->Draw(); qlay_cut->Draw("same");
  qlay_raw->GetXaxis()->SetTitle("Z"); qlay_raw->GetYaxis()->SetTitle("Events");
  Q_lay_MC->SetLogx(0); Q_lay_MC->SetLogy(1); qlay_raw->SetMinimum(1);
  qlay_raw->SetLineColor(kRed); qlay_cut->SetLineColor(kBlue);
  Q_lay_MC->Draw(); Q_lay_MC->Print("./ProtonAnalysis/PA/First Tracker Layer Charge MC.png");
  // RICH Mass MC
  Mass_RICH_MC->cd();
  mass_raw->Draw(); mass_cut->Draw("same");
  mass_raw->GetXaxis()->SetTitle("m [GeV/c#^2]"); mass_cut->GetYaxis()->SetTitle("Events");
  Mass_RICH_MC->SetLogx(0); Mass_RICH_MC->SetLogy(1); mass_raw->SetMinimum(1);
  mass_raw->SetLineColor(kRed); mass_cut->SetLineColor(kBlue);
  Mass_RICH_MC->Draw(); Mass_RICH_MC->Print("./ProtonAnalysis/PA/RICH Mass MC.png");

  cout << "ParameterAnalysis() has finished!\n" << endl;

}



//------------------------------------------------------------------------------------------------------------------------------------------------------------
// Returns TH1F of the selected events per rigidity bin
void Anaaqra::RigBinner() {

  cout << "Running RigBinner()..." << endl;

  // Empty the event histograms
  Events_raw->Reset("ICESM");
  Events_cut->Reset("ICESM");

  // Loop over all Compact entries
  for (int i=0; i<Simp_chain->GetEntries(); i++) {

    // Get entry
    Simp_chain->GetEntry(i);

    // Fill histograms
    Events_raw->Fill(Tool->trk_rig);
    if (EventSelectorSimple(Tool, "1111111_111")) {
      Events_cut->Fill(Tool->trk_rig);
    }

  }

  // TGraph
  for (int i=0; i<Bin_num; i++) {
    Ev_raw[i] = Events_raw->GetBinContent(i+1);
    Ev_cut[i] = Events_cut->GetBinContent(i+1);
    Ev_raw_err[i] = TMath::Sqrt(Events_raw->GetBinContent(i+1));
    Ev_cut_err[i] = TMath::Sqrt(Events_cut->GetBinContent(i+1));
  }
  TGraphErrors* ev_raw_graph = new TGraphErrors(32, Bin_mid, Ev_raw, Bin_err, Ev_raw_err);
  TGraphErrors* ev_cut_graph = new TGraphErrors(32, Bin_mid, Ev_cut, Bin_err, Ev_cut_err);
  // Canvas
  TCanvas* c_Events = new TCanvas("c_Events", "Events per Rigitidy Bin");
  ev_raw_graph->Draw("AP");
  ev_cut_graph->Draw("P");
  // Styling
  ev_raw_graph->SetMarkerStyle(20);
  ev_raw_graph->SetMarkerSize(1);
  ev_raw_graph->SetMarkerColor(kRed);
  ev_cut_graph->SetMarkerStyle(20);
  ev_cut_graph->SetMarkerSize(1);
  ev_cut_graph->SetMarkerColor(kBlue);
  // Axes
  c_Events->SetLogy();
  ev_raw_graph->SetMinimum(1);
  ev_raw_graph->GetXaxis()->SetTitle("R [GV]");
  ev_raw_graph->GetYaxis()->SetTitle("Events");
  // Print
  c_Events->Draw();
  c_Events->Print("./ProtonAnalysis/(Selected) Events.png");

  cout << "RigBinner() has finished!\n" << endl;

}



//------------------------------------------------------------------------------------------------------------------------------------------------------------
// Returns TH1F of the exposure time
void Anaaqra::Exposure() {

  cout << "Running Exposure()..." << endl;

  // Loop over RTI entries
  for (int i=0; i<RTII_chain->GetEntries(); i++) {

    // Get entry
    RTII_chain->GetEntry(i);

    //Loop over rigidity bins
    for (int j=0; j<Bin_num; j++) {

      // If the centre of the bin is above the geo-magnetic cut-off, add the livetime
      if (Bin_mid[j] > 1.2 * Woi->cf) {
        ExposureTime->SetBinContent(j+1, ExposureTime->GetBinContent(j+1) + Woi->lf);
      }

    }

  }

  // TGraph
  for (int i=0; i<Bin_num; i++) {
    ExpTime[i] = ExposureTime->GetBinContent(i+1);
  }
  TGraphErrors* exptime_graph = new TGraphErrors(32, Bin_mid, ExpTime, Bin_err, 0);
  // Canvas
  TCanvas* c_Exposure = new TCanvas("c_Exposure", "Exposure Time per Rigitidy Bin");
  exptime_graph->Draw("AP");
  // Styling
  exptime_graph->SetMarkerStyle(20);
  exptime_graph->SetMarkerSize(1);
  exptime_graph->SetMarkerColor(kRed);
  // Axes
  c_Exposure->SetLogy();
  exptime_graph->GetXaxis()->SetTitle("R [GV]");
  exptime_graph->GetYaxis()->SetTitle("Exposure Time [s]");
  // Print
  c_Exposure->Draw();
  c_Exposure->Print("./ProtonAnalysis/Exposure Time.png");

    cout << "Exposure() has finished!\n" << endl;

}



//------------------------------------------------------------------------------------------------------------------------------------------------------------
// Returns TH1F of the geometric Acceptance
void Anaaqra::Acceptance(bool apply_cuts = 0) {

  cout << "Running Acceptance()..." << endl;

  // Clear histograms
  MC_detected->Reset("ICESM");
  MC_generated->Reset("ICESM");

  // Loop over MC root files individually
  for (int i=0; i<13; i++) {

    // Import tree
    TChain mc_chain("Compact");
    TChain fi_chain("File");
    mc_chain.Add(Form("../MC Protons/%d.root", Start_num + i));
    fi_chain.Add(Form("../MC Protons/%d.root", Start_num + i));
    // Create empty class objects
    NtpCompact *comp = new class NtpCompact();
    FileMCInfo *fmci = new class FileMCInfo();
    // Set branch addresses
    mc_chain.SetBranchAddress("Compact", &comp);
    fi_chain.SetBranchAddress("FileMCInfo", &fmci);

    // Temporary parameters
    double ngen = 0;
    double rig_min = 0;
    double rig_max = 0;
    // Get FileMCInfo entry for information about MC gegneration
    fi_chain.GetEntry(0);
    // Get parameters
    ngen = fmci->ngen_datacard;
    rig_min = fmci->momentum[0];
    rig_max = fmci->momentum[1];

    // 1/R generation spectum
    TF1 *genFlux = new TF1("genFlux", "[0]/(x)", rig_min, rig_max);
    // Normalisation
    genFlux->SetParameter(0, log(rig_max / rig_min));

    // Loop over rigidity bins
    for (int j=0; j<Bin_num; j++) {

      // Fraction of spectrum in bin
      double frac = genFlux->Integral(Bin_edges[j], Bin_edges[j+1]) / genFlux->Integral(rig_min, rig_max);
      // Number of events in fraction
      double nev = frac * ngen;

      // Set bin to current bin count + nev
      MC_generated->SetBinContent(j+1, MC_generated->GetBinContent(j+1) + nev);

    }

  }

  // Loop over MC Compact entries
  for (int i=0; i<MC_chain->GetEntries(); i++) {

    // Get entry
    MC_chain->GetEntry(i);

    // Apply cuts
    if (apply_cuts) {
      if (EventSelectorCompact(MC_comp, "1111111_111")) {
        MC_detected->Fill(MC_comp->trk_rig[0]);
      }
    } else {
      MC_detected->Fill(MC_comp->trk_rig[0]);
    }

  }

  // Fill Acceptance histogram
  // Loop over rigidity bins
  for (int i=0; i<Bin_num; i++) {
    AcceptHist->SetBinContent(i+1, MC_detected->GetBinContent(i+1) / MC_generated->GetBinContent(i+1));
  }
  // Scale histogram by generation volume
  double gen_vol = TMath::Pi() * 3.9 * 3.9;
  AcceptHist->Scale(gen_vol);

  // TGraph
  for (int i=0; i<Bin_num; i++) {
    Accept[i] = AcceptHist->GetBinContent(i+1);
  }
  TGraphErrors* accept_graph = new TGraphErrors(32, Bin_mid, Accept, Bin_err, 0);
  // Canvas
  TCanvas* c_Accept = new TCanvas("c_Accept", "Acceptance per Rigitidy Bin");
  accept_graph->Draw("AP");
  // Styling
  accept_graph->SetMarkerStyle(20);
  accept_graph->SetMarkerSize(1);
  accept_graph->SetMarkerColor(kRed);
  // Axes
  accept_graph->SetMinimum(0); accept_graph->SetMaximum(0.4);
  accept_graph->GetXaxis()->SetTitle("R [GV]");
  accept_graph->GetYaxis()->SetTitle("Acceptance [m^2 sr]");
  // Print
  c_Accept->Draw();
  c_Accept->Print("./ProtonAnalysis/Acceptance.png");

  cout << "Acceptance() has finished!\n" << endl;

}



//------------------------------------------------------------------------------------------------------------------------------------------------------------
// Returns TH1F of the selection efficiency
void Anaaqra::CutEff(bool plot_all = 0) {

  cout << "Running CutEff()..." << endl;

  // Temporary histograms
  TH1F *Cpar_MC = new TH1F("Cpar_MC", "Single Particle Cut", 32, Bin_edges);
  TH1F *Cbet_MC = new TH1F("Cbet_MC", "Downward Particle Cut", 32, Bin_edges);
  TH1F *Cchi_MC = new TH1F("Cchi_MC", "Well-constructed Track Cut", 32, Bin_edges);
  TH1F *Cinn_MC = new TH1F("Cinn_MC", "Inner Tracker Charge Cut", 32, Bin_edges);
  TH1F *Btof_MC = new TH1F("Btof_MC", "TOF Base Cut", 32, Bin_edges);
  TH1F *Btrk_MC = new TH1F("Btrk_MC", "Tracker Base Cut", 32, Bin_edges);

  // Loop over MC entries
  for (int i=0; i<MC_chain->GetEntries(); i++) {

    // Get entry
    MC_chain->GetEntry(i);

    bool tof_q = (MC_comp->tof_q_lay[0] > 0.8)&&(MC_comp->tof_q_lay[0] < 1.5);

    // Fill base histograms
    if (EventSelectorCompact(MC_comp, "110011x_111")) {Btof_MC->Fill(MC_comp->trk_rig[0]);}
    if (EventSelectorCompact(MC_comp, "111100x_111") && tof_q) {Btrk_MC->Fill(MC_comp->trk_rig[0]);}

    // Fill cut histograms
    if (EventSelectorCompact(MC_comp, "111111x_111")) {Cpar_MC->Fill(MC_comp->trk_rig[0]);}
    if (EventSelectorCompact(MC_comp, "110111x_111")) {Cbet_MC->Fill(MC_comp->trk_rig[0]);}
    if (EventSelectorCompact(MC_comp, "111100x_111") && tof_q) {Cchi_MC->Fill(MC_comp->trk_rig[0]);}
    if (EventSelectorCompact(MC_comp, "111010x_111") && tof_q) {Cinn_MC->Fill(MC_comp->trk_rig[0]);}

  }

  // Cut Efficiency MC Error
  for (int i=0; i<Bin_num; i++) {
    CE_mc_err[i] = TMath::Sqrt((1/Cpar_MC->GetBinContent(i+1) + 1/Btof_MC->GetBinContent(i+1))
                   + (1/Cbet_MC->GetBinContent(i+1) + 1/Btof_MC->GetBinContent(i+1))
                   + (1/Cchi_MC->GetBinContent(i+1) + 1/Btrk_MC->GetBinContent(i+1))
                   + (1/Cinn_MC->GetBinContent(i+1) + 1/Btrk_MC->GetBinContent(i+1)));
  }

  // Divide by corresponding instrument base
  Cpar_MC->Divide(Btof_MC);
  Cbet_MC->Divide(Btof_MC);
  Cchi_MC->Divide(Btrk_MC);
  Cinn_MC->Divide(Btrk_MC);

  // Loop over rigidity bins
  for (int i=0; i<Bin_num; i++) {
    CutEff_MC->SetBinContent(i+1, Cpar_MC->GetBinContent(i+1)
                                * Cbet_MC->GetBinContent(i+1)
                                * Cchi_MC->GetBinContent(i+1)
                                * Cinn_MC->GetBinContent(i+1));
  }

  // Temporary histograms
  TH1F *Cpar_data = new TH1F("Cpar_data", "Single Particle Cut", 32, Bin_edges);
  TH1F *Cbet_data = new TH1F("Cbet_data", "Downward Particle Cut", 32, Bin_edges);
  TH1F *Cchi_data = new TH1F("Cchi_data", "Well-constructed Track Cut", 32, Bin_edges);
  TH1F *Cinn_data = new TH1F("Cinn_data", "Inner Tracker Charge Cut", 32, Bin_edges);
  TH1F *Btof_data = new TH1F("Btof_data", "TOF Base Cut", 32, Bin_edges);
  TH1F *Btrk_data = new TH1F("Btrk_data", "Tracker Base Cut", 32, Bin_edges);

  // Loop over data entries
  for (int i=0; i<Simp_chain->GetEntries(); i++) {

    // Get entry
    Simp_chain->GetEntry(i);

    bool tof_q = (Tool->tof_q_lay[0] > 0.8)&&(Tool->tof_q_lay[0] < 1.5);

    // Fill base histograms
    if (EventSelectorSimple(Tool, "1100111_111")) {Btof_data->Fill(Tool->trk_rig);}
    if (EventSelectorSimple(Tool, "1111001_111") && tof_q) {Btrk_data->Fill(Tool->trk_rig);}

    // Fill cut histograms
    if (EventSelectorSimple(Tool, "1110111_111")) {Cpar_data->Fill(Tool->trk_rig);}
    if (EventSelectorSimple(Tool, "1101111_111")) {Cbet_data->Fill(Tool->trk_rig);}
    if (EventSelectorSimple(Tool, "1111101_111") && tof_q) {Cchi_data->Fill(Tool->trk_rig);}
    if (EventSelectorSimple(Tool, "1111011_111") && tof_q) {Cinn_data->Fill(Tool->trk_rig);}

  }

  // Cut Efficiency Data Error
  for (int i=0; i<Bin_num; i++) {
    CE_data_err[i] = TMath::Sqrt((1/Cpar_data->GetBinContent(i+1) + 1/Btof_data->GetBinContent(i+1))
                     + (1/Cbet_data->GetBinContent(i+1) + 1/Btof_data->GetBinContent(i+1))
                     + (1/Cchi_data->GetBinContent(i+1) + 1/Btrk_data->GetBinContent(i+1))
                     + (1/Cinn_data->GetBinContent(i+1) + 1/Btrk_data->GetBinContent(i+1)));
  }

  // Divide by corresponding instrument base
  Cpar_data->Divide(Btof_data);
  Cbet_data->Divide(Btof_data);
  Cchi_data->Divide(Btrk_data);
  Cinn_data->Divide(Btrk_data);

  // Loop over rigidity bins
  for (int i=0; i<Bin_num; i++) {
    CutEff_data->SetBinContent(i+1, Cpar_data->GetBinContent(i+1)
                                  * Cbet_data->GetBinContent(i+1)
                                  * Cchi_data->GetBinContent(i+1)
                                  * Cinn_data->GetBinContent(i+1));
  }

  if (plot_all) {

    // Canvasses
    // Single particle cut
    TCanvas *c_spc = new TCanvas("c_spc", "Single Particle Cut Efficiency");
    Cpar_MC->Draw(); Cpar_data->Draw("same");
    Cpar_MC->SetLineColor(kRed); Cpar_data->SetLineColor(kBlue);
    Cpar_MC->SetLineWidth(2); Cpar_data->SetLineWidth(2);
    Cpar_MC->SetMinimum(0); Cpar_MC->SetMaximum(1.05);
    Cpar_MC->GetXaxis()->SetTitle("R [GV]"); Cpar_MC->GetYaxis()->SetTitle("Cut Efficiency");
    c_spc->Draw(); c_spc->Print("./ProtonAnalysis/CE/Single Particle Cut Efficiency.png");
    // Downward going track cut
    TCanvas *c_dgc = new TCanvas("c_dgc", "Downward Particle Cut Efficiency");
    Cbet_MC->Draw(); Cbet_data->Draw("same");
    Cbet_MC->SetLineColor(kRed); Cbet_data->SetLineColor(kBlue);
    Cbet_MC->SetLineWidth(2); Cbet_data->SetLineWidth(2);
    Cbet_MC->SetMinimum(0); Cbet_MC->SetMaximum(1.05);
    Cbet_MC->GetXaxis()->SetTitle("R [GV]"); Cbet_MC->GetYaxis()->SetTitle("Cut Efficiency");
    c_dgc->Draw(); c_dgc->Print("./ProtonAnalysis/CE/Downward Particle Cut Efficiency.png");
    // Well constructed track cut
    TCanvas *c_wtc = new TCanvas("c_wtc", "Well-constructed Track Cut Efficiency");
    Cchi_MC->Draw(); Cchi_data->Draw("same");
    Cchi_MC->SetLineColor(kRed); Cchi_data->SetLineColor(kBlue);
    Cchi_MC->SetLineWidth(2); Cchi_data->SetLineWidth(2);
    Cchi_MC->SetMinimum(0); Cchi_MC->SetMaximum(1.05);
    Cchi_MC->GetXaxis()->SetTitle("R [GV]"); Cchi_MC->GetYaxis()->SetTitle("Cut Efficiency");
    c_wtc->Draw(); c_wtc->Print("./ProtonAnalysis/CE/Well-constructed Track Cut Efficiency.png");
    // Inner tracker charge cut
    TCanvas *c_itc = new TCanvas("c_itc", "Inner Tracker Charge Cut Efficiency");
    Cinn_MC->Draw(); Cinn_data->Draw("same");
    Cinn_MC->SetLineColor(kRed); Cinn_data->SetLineColor(kBlue);
    Cinn_MC->SetLineWidth(2); Cinn_data->SetLineWidth(2);
    Cinn_MC->SetMinimum(0); Cinn_MC->SetMaximum(1.05);
    Cinn_MC->GetXaxis()->SetTitle("R [GV]"); Cinn_MC->GetYaxis()->SetTitle("Cut Efficiency");
    c_itc->Draw(); c_itc->Print("./ProtonAnalysis/CE/Inner Tracker Charge Cut Efficiency.png");

  }

  // TGraph
  for (int i=0; i<Bin_num; i++) {
    CE_mc[i] = CutEff_MC->GetBinContent(i+1);
    CE_data[i] = CutEff_data->GetBinContent(i+1);
    CE_mc_err[i] *= CE_mc[i];
    CE_data_err[i] *= CE_data[i];
  }
  TGraphErrors* ce_mc_graph = new TGraphErrors(32, Bin_mid, CE_mc, Bin_err, CE_mc_err);
  TGraphErrors* ce_data_graph = new TGraphErrors(32, Bin_mid, CE_data, Bin_err, CE_data_err);
  // Canvas
  TCanvas* c_CutEff = new TCanvas("c_CutEff", "Selection Efficiency per Rigitidy Bin");
  ce_mc_graph->Draw("AP");
  ce_data_graph->Draw("P");
  // Styling
  ce_mc_graph->SetMarkerStyle(20);
  ce_mc_graph->SetMarkerSize(1);
  ce_mc_graph->SetMarkerColor(kRed);
  ce_data_graph->SetMarkerStyle(20);
  ce_data_graph->SetMarkerSize(1);
  ce_data_graph->SetMarkerColor(kBlue);
  // Axes
  ce_mc_graph->SetMaximum(1);
  ce_mc_graph->SetMinimum(0);
  ce_mc_graph->GetXaxis()->SetTitle("R [GV]");
  ce_mc_graph->GetYaxis()->SetTitle("Selection Efficiency");
  // Print
  c_CutEff->Draw();
  c_CutEff->Print("./ProtonAnalysis/Selection Efficiency.png");

  cout << "CutEff() has finished!\n" << endl;

}



//------------------------------------------------------------------------------------------------------------------------------------------------------------
// Returns TH1F of the trigger efficiency
void Anaaqra::TrigEff(int delta = 100) {

  cout << "Running TrigEff()..." << endl;

  // Loop over MC entries
  for (int i=0; i<MC_chain->GetEntries(); i++) {

    // Get entry
    MC_chain->GetEntry(i);

    // Check for physical trigger
    bool HasPhysTrig_mc = ((MC_comp->sublvl1&0x3E)!=0)&&((MC_comp->trigpatt&0x2)!=0);
    bool HasUnphTrig_mc = ((MC_comp->sublvl1&0x3E)==0)&&((MC_comp->trigpatt&0x2)!=0);

    // Fill histograms
    if (EventSelectorCompact(MC_comp, "101111x_111")) {
      if (HasPhysTrig_mc) {
        PhysHist_mc->Fill(MC_comp->trk_rig[0]);
      }
      if (HasUnphTrig_mc) {
        UnphHist_mc->Fill(MC_comp->trk_rig[0]);
      }
    }

  }

  // Loop over data entries
  for (int i=0; i<Simp_chain->GetEntries(); i++) {

    // Get entry
    Simp_chain->GetEntry(i);

    // Check for physical trigger
    bool HasPhysTrig_data = ((Tool->sublvl1&0x3E)!=0)&&((Tool->trigpatt&0x2)!=0);
    bool HasUnphTrig_data = ((Tool->sublvl1&0x3E)==0)&&((Tool->trigpatt&0x2)!=0);

    // Fill histograms
    if (EventSelectorSimple(Tool, "1011111_111")) {
      if (HasPhysTrig_data) {
        PhysHist_data->Fill(Tool->trk_rig);
      }
      if (HasUnphTrig_data) {
        UnphHist_data->Fill(Tool->trk_rig);
      }
    }

  }

  // Loop over rigidity bins
  for (int i=0; i<Bin_num; i++) {

    // Fill histograms
    if (PhysHist_mc->GetBinContent(i+1) == 0) {
      TrigEff_MC->SetBinContent(i+1, 0);
    } else {
      TrigEff_MC->SetBinContent(i+1, PhysHist_mc->GetBinContent(i+1) / (PhysHist_mc->GetBinContent(i+1) + UnphHist_mc->GetBinContent(i+1)));
    }
    if (PhysHist_data->GetBinContent(i+1) == 0) {
      TrigEff_data->SetBinContent(i+1, 0);
    } else {
      TrigEff_data->SetBinContent(i+1, PhysHist_data->GetBinContent(i+1) / (PhysHist_data->GetBinContent(i+1) + delta * UnphHist_data->GetBinContent(i+1)));
    }

  }

  // TGraph
  for (int i=0; i<Bin_num; i++) {
    TE_mc[i] = TrigEff_MC->GetBinContent(i+1);
    TE_data[i] = TrigEff_data->GetBinContent(i+1);
    TE_mc_err[i] = delta * TMath::Sqrt(PhysHist_mc->GetBinContent(i+1)*pow(UnphHist_mc->GetBinContent(i+1), 2)
                   + UnphHist_mc->GetBinContent(i+1)*pow(PhysHist_mc->GetBinContent(i+1), 2))
                   / pow((PhysHist_mc->GetBinContent(i+1) + delta * UnphHist_mc->GetBinContent(i+1)), 2);
    TE_data_err[i] = delta * TMath::Sqrt(PhysHist_data->GetBinContent(i+1)*pow(UnphHist_data->GetBinContent(i+1), 2)
                     + UnphHist_data->GetBinContent(i+1)*pow(PhysHist_data->GetBinContent(i+1), 2))
                     / pow((PhysHist_data->GetBinContent(i+1) + delta * UnphHist_data->GetBinContent(i+1)), 2);
  }
  TGraphErrors* te_mc_graph = new TGraphErrors(32, Bin_mid, TE_mc, Bin_err, TE_mc_err);
  TGraphErrors* te_data_graph = new TGraphErrors(32, Bin_mid, TE_data, Bin_err, TE_data_err);
  // Canvas
  TCanvas* c_TrigEff = new TCanvas("c_TrigEff", "Trigger Efficiency per Rigitidy Bin");
  te_mc_graph->Draw("AP");
  te_data_graph->Draw("P");
  // Styling
  te_mc_graph->SetMarkerStyle(20);
  te_mc_graph->SetMarkerSize(1);
  te_mc_graph->SetMarkerColor(kRed);
  te_data_graph->SetMarkerStyle(20);
  te_data_graph->SetMarkerSize(1);
  te_data_graph->SetMarkerColor(kBlue);
  // Axes
  //c_TrigEff->SetLogy();
  te_mc_graph->SetMaximum(1.0);
  te_mc_graph->SetMinimum(0.0);
  te_mc_graph->GetXaxis()->SetTitle("R [GV]");
  te_mc_graph->GetYaxis()->SetTitle("Trigger Efficiency");
  // Print
  c_TrigEff->Draw();
  c_TrigEff->Print("./ProtonAnalysis/Trigger Efficiency.png");

  cout << "TrigEff() has finished!\n" << endl;

}



//------------------------------------------------------------------------------------------------------------------------------------------------------------
// Returns TH1F of the proton rate
void Anaaqra::Rate() {

  // Fill empty necessary histograms
  if (Events_cut->GetEntries() == 0) {RigBinner();}
  if (ExposureTime->GetEntries() == 0) {Exposure();}

  // Loop over rigidity bins
  for (int i=0; i<Bin_num; i++) {
    if (ExposureTime->GetBinContent(i+1) == 0) {
      RateHist->SetBinContent(i+1, 0);
      PRate_err[i] = 0;
    } else {
      RateHist->SetBinContent(i+1, Events_cut->GetBinContent(i+1) / ExposureTime->GetBinContent(i+1) / (2 * Bin_err[i]));
      PRate_err[i] = RateHist->GetBinContent(i+1) / TMath::Sqrt(Events_cut->GetBinContent(i+1));
    }
  }

  // TGraph
  for (int i=0; i<Bin_num; i++) {PRate[i] = RateHist->GetBinContent(i+1);}
  TGraphErrors* p_rate_graph = new TGraphErrors(32, Bin_mid, PRate, Bin_err, PRate_err);
  // Canvas
  TCanvas* c_Rate = new TCanvas("c_Rate", "Proton Rate per Rigitidy Bin");
  p_rate_graph->Draw("AP");
  // Styling
  p_rate_graph->SetMarkerStyle(20);
  p_rate_graph->SetMarkerSize(1);
  p_rate_graph->SetMarkerColor(kRed);
  // Axes
  c_Rate->SetLogy();
  p_rate_graph->GetXaxis()->SetTitle("R [GV]");
  p_rate_graph->GetYaxis()->SetTitle("Rate [s^-1]");
  // Print
  c_Rate->Draw();
  c_Rate->Print("./ProtonAnalysis/Proton Rate.png");

  cout << "Rate() has finished!\n" << endl;

}



//------------------------------------------------------------------------------------------------------------------------------------------------------------
// Returns TH1F of the proton flux
void Anaaqra::Flux(bool comp = 0) {

  // Fill empty necessary histograms
  if (Events_cut->GetEntries() == 0) {RigBinner();}
  if (ExposureTime->GetEntries() == 0) {Exposure();}
  if (AcceptHist->GetEntries() == 0) {Acceptance();}
  if ((CutEff_MC->GetEntries() == 0) || (CutEff_data->GetEntries() == 0)) {CutEff();}
  if ((TrigEff_MC->GetEntries() == 0) || (TrigEff_data->GetEntries() == 0)) {TrigEff();}

  // Loop over rigidity bins
  for (int i=0; i<Bin_num; i++) {
    if (ExposureTime->GetBinContent(i+1) == 0) {
      FluxHist->SetBinContent(i+1, 0);
      PFlux_err[i] = 0;
    } else {
      FluxHist->SetBinContent(i+1, Events_cut->GetBinContent(i+1) / ExposureTime->GetBinContent(i+1)
                                   / 2 / Bin_err[i] / AcceptHist->GetBinContent(i+1)
                                   / CutEff_data->GetBinContent(i+1) * CutEff_MC->GetBinContent(i+1)
                                   / TrigEff_data->GetBinContent(i+1) * TrigEff_MC->GetBinContent(i+1)
                             );
      PFlux_err[i] = FluxHist->GetBinContent(i+1) * TMath::Sqrt(1/Events_cut->GetBinContent(i+1)
                     + pow(TE_mc_err[i]/TE_mc[i], 2) + pow(TE_data_err[i]/TE_data[i], 2)
                     + pow(CE_mc_err[i]/CE_mc[i], 2) + pow(CE_data_err[i]/CE_data[i], 2));
    }
  }

  // TGraph 1/2
  for (int i=0; i<Bin_num; i++) {
    PFlux[i] = FluxHist->GetBinContent(i+1);
  }
  TGraphErrors* p_flux_graph = new TGraphErrors(32, Bin_mid, PFlux, Bin_err, PFlux_err);
  // Canvas
  TCanvas* c_Flux = new TCanvas("c_Flux", "Proton Flux per Rigitidy Bin");
  p_flux_graph->Draw("AP");
  // Styling
  p_flux_graph->SetMarkerStyle(20);
  p_flux_graph->SetMarkerSize(1);
  p_flux_graph->SetMarkerColor(kRed);
  // Axes
  c_Flux->SetLogy();
  p_flux_graph->GetXaxis()->SetTitle("R [GV]");
  p_flux_graph->GetYaxis()->SetTitle("Flux [m^-2 sr^-1 s^-1 GV^-1]");
  // Print
  c_Flux->Draw();
  c_Flux->Print("./ProtonAnalysis/Proton Flux.png");

  // Comparison to SSDC data
  TFile *ssdc = new TFile("../SSDC_PubPFlux/ssdc_canvas.root","open");
  TGraphAsymmErrors *pubpflux = (TGraphAsymmErrors *)ssdc->Get("graph1");

  // TGraph 2/2
  for (int i=0; i<Bin_num; i++) {
    PSFlux[i] = FluxHist->GetBinContent(i+1) * pow(Bin_mid[i], 2.7);
    PSFlux_err[i] = PFlux_err[i] * pow(Bin_mid[i], 2.7);

    double SSDC_PFx; double SSDC_PFy;
    pubpflux->GetPoint(i, SSDC_PFx, SSDC_PFy);
    PSFlux_SSDC_Ratio[i] = PSFlux[i] / SSDC_PFy;
  }
  TGraphErrors* p_sflux_graph = new TGraphErrors(32, Bin_mid, PSFlux, Bin_err, PSFlux_err);
  TGraph* p_sflux_ratio = new TGraph(32, Bin_mid, PSFlux_SSDC_Ratio);
  // Canvas
  TCanvas* c_SFlux = new TCanvas("c_SFlux", "Scaled Proton Flux per Rigitidy Bin");
  c_SFlux->Divide(1,2);
  // First split
  c_SFlux->cd(1);
  p_sflux_graph->Draw("AP");
  pubpflux->Draw("P");
  // Styling
  p_sflux_graph->SetMarkerStyle(20);
  p_sflux_graph->SetMarkerSize(1);
  p_sflux_graph->SetMarkerColor(kRed);
  // Axes
  c_SFlux->SetLogy();
  p_sflux_graph->GetXaxis()->SetTitle("R [GV]");
  p_sflux_graph->GetYaxis()->SetTitle("Flux R^2.7 [m^-2 sr^-1 s^-1 GV^1.7]");
  p_sflux_graph->SetMaximum(11000);
  // Second split
  c_SFlux->cd(2);
  p_sflux_ratio->Draw("AP");
  // Styling
  p_sflux_ratio->SetMarkerStyle(20);
  p_sflux_ratio->SetMarkerSize(1);
  p_sflux_ratio->SetMarkerColor(kBlack);
  // Axes
  c_SFlux->SetLogy(0);
  p_sflux_ratio->GetXaxis()->SetTitle("R [GV]");
  p_sflux_ratio->GetYaxis()->SetTitle("Flux Ratio");
  p_sflux_ratio->SetMaximum(1); p_sflux_ratio->SetMinimum(0);

  // Print
  c_SFlux->Draw();
  c_SFlux->Print("./ProtonAnalysis/Scaled Proton Flux.png");

  cout << "Flux() has finished!\n" << endl;

}



//------------------------------------------------------------------------------------------------------------------------------------------------------------
// Aerogel beta mass test
void Anaaqra::AerogelSlice(){

  cout << "Running AerogelSlice()..." << endl;

  // Histograms
  TH2D *aero_beta_mass = new TH2D("abm", "", 100, 0.94, 1.02, 100, 0, 2.5);

  // Loop over data entries
  for (int i=0; i<Simp_chain->GetEntries(); i++) {

    // Get Entry
    Simp_chain->GetEntry(i);

    // Fill histograms
    if (EventSelectorSimple(Tool, "1111111_111")) {
      double beta = Tool->rich_beta;
      double mass = Tool->trk_q_inn * Tool->trk_rig * TMath::Sqrt(1 / Tool->rich_beta / Tool->rich_beta - 1);
      aero_beta_mass->Fill(beta, mass);
      }
    }

  // Create canvasses
  TCanvas *ABMcanvas = new TCanvas("ABMcanvas", "Aerogel Beta Mass");
  aero_beta_mass->Draw("COLZ");
  ABMcanvas->SetLogx(0); ABMcanvas->SetLogy(0); ABMcanvas->SetLogz(1);
  ABMcanvas->Draw(); ABMcanvas->Print("./ProtonAnalysis/Aerogel Beta Mass.png");

}
