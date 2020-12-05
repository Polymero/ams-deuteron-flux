// C++ class for full PROTON analysis
// Created        23-11-20
// Last Edited    23-11-20

// Include header file(s)
#include <iostream>
#include <string>
#include <algorithm>
#include "Header Files/Ntp.h"
#include "Header Files/Simple.h"
#include "TChain.h"
#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TString.h"

//------------------------------------------------------------------------------
// CLASS DEFINITION
//------------------------------------------------------------------------------

// Class
class Anaaqra {
  public: // Acces Specifier

    //--------------------------------------------------------------------------
    // ATTRIBUTES
    //--------------------------------------------------------------------------
    // Rigidity bins
    const int Bin_num = 32;
    double Bin_edges[32+1] = {1.00,1.16,1.33,1.51,1.71,1.92,2.15,2.40,2.67,2.97,
                            3.29,3.64,4.02,4.43,4.88,5.37,5.90,6.47,7.09,7.76,
                            8.48,9.26,10.1,11.0,12.0,13.0,14.1,15.3,16.6,18.0,
                            19.5,21.1,22.8};
    double Bin_err[32];
    double Bin_mid[32];
    // Histograms
    TH1F *Events_raw    = new TH1F("Events_raw", "Raw AMS-02 Events", 32, Bin_edges);
    TH1F *Events_cut    = new TH1F("Events_cut", "Selected AMS-02 Events", 32, Bin_edges);
    TH1F *ExposureTime  = new TH1F("ExposureTime", "Exposure Time per Rigidity Bin", 32, Bin_edges);
    TH1F *RateHist      = new TH1F("RateHist", "Proton Rate per Rigidity Bin", 32, Bin_edges);
    TH1F *MC_generated  = new TH1F("MC_generated", "Generated MC Events per Rigidity Bin", 32, Bin_edges);
    TH1F *MC_detected   = new TH1F("MC_detected", "Detected MC Events per Rigidity Bin", 32, Bin_edges);
    TH1F *AcceptHist    = new TH1F("AcceptHist", "Acceptance per Rigidity Bin", 32, Bin_edges);
    TH1F *CutEff_data   = new TH1F("CutEff_data", "Selection Efficiency of Data per Rigidity Bin", 32, Bin_edges);
    TH1F *CutEff_MC     = new TH1F("CutEff_MC", "Selection Efficiency of MC per Rigidity Bin", 32, Bin_edges);
    TH1F *TrigEff_data  = new TH1F("TrigEff_data", "Trigger Efficiency of Data per Rigidity Bin", 32, Bin_edges);
    TH1F *TrigEff_MC    = new TH1F("TrigEff_MC", "Trigger Efficiency of MC per Rigidity Bin", 32, Bin_edges);
    // Data objects
    TChain *Simp_chain  = new TChain("Simp");
    TChain *RTII_chain  = new TChain("RTIInfo");
    TChain *MC_chain    = new TChain("Compact");
    Miiqtool *Tool      = new Miiqtool();
    MiiqRTI *Woi        = new MiiqRTI();
    NtpCompact *MC_comp = new NtpCompact();

    //--------------------------------------------------------------------------
    // CONSTRUCTORS
    //--------------------------------------------------------------------------
    Anaaqra() { // Default constructor

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

      cout << "Class succesfully constructed!" << endl;

    };

    //--------------------------------------------------------------------------
    // LIST OF METHODS
    //--------------------------------------------------------------------------
    // Singular (independent)
    TH1F RigBinner();                           // Returns number of selected events as function of rigidity
    TH1F Exposure();                            // Returns exposure time (livetime) as function of rigidity
    TH1F Acceptance(bool apply_cuts = 0);       // Returns (geometric) acceptance as function of rigidity
    TH1F CutEff();                              // Returns the cut (selection) efficiency as function of rigidity
    const char* TrigEff();                             // Returns the trigger efficiency as function of rigidity
    // Plural (dependent)
    TH1F ProtonRate();                          // Returns the proton rate as function of rigidity
    TH1F ProtonFlux();                          // Returns the proton flux as function of rigidity
    // Debug
    int Debug();

};

//------------------------------------------------------------------------------
// ASSIST FUNCTIONS
//------------------------------------------------------------------------------
// Returns bool if event passed specified selection of cuts
bool EventSelectorCompact(NtpCompact* comp, const char* cutbit) {

  bool pass = 1;

  // Boolean cuts (instead of TCut)
  bool Crig = (comp->trk_rig[0] > 0)&&(comp->trk_rig[0] <= 22);
  bool Ctrg = ((comp->sublvl1&0x3E)!=0)&&((comp->trigpatt&0x2)!=0);
  bool Cpar = comp->status % 10 == 1;
  bool Ccon = std::abs(comp->tof_beta - comp->rich_beta)/comp->tof_beta < 0.05;
  bool Cbet = comp->tof_beta > 0;
  bool Cchi = (comp->trk_chisqn[0][0] < 10)&&(comp->trk_chisqn[0][1] < 10)&&(comp->trk_chisqn[0][0] > 0)&&(comp->trk_chisqn[0][1] > 0);
  bool Cinn = (comp->trk_q_inn > 0.80)&&(comp->trk_q_inn < 1.30);
  bool Clay = (comp->trk_q_lay[0] >= 0)&&(comp->trk_q_lay[1] >= 0)&&(comp->trk_q_lay[2] >= 0)&&(comp->trk_q_lay[3] >= 0)&&(comp->trk_q_lay[4] >= 0)&&(comp->trk_q_lay[5] >= 0)&&(comp->trk_q_lay[6] >= 0)&&(comp->trk_q_lay[7] >= 0)&&(comp->trk_q_lay[8] >= 0);

  // Adjust return bool according to cutbit
  if (cutbit[0] == '1') {pass &= Crig;}
  if (cutbit[1] == '1') {pass &= Ctrg;}
  if (cutbit[2] == '1') {pass &= Cpar;}
  if (cutbit[3] == '1') {pass &= Ccon;}
  if (cutbit[4] == '1') {pass &= Cbet;}
  if (cutbit[5] == '1') {pass &= Cchi;}
  if (cutbit[6] == '1') {pass &= Cinn;}
  if (cutbit[7] == '1') {pass &= Clay;}

  // Return
  return pass;

}



//------------------------------------------------------------------------------
// METHOD FUNCTIONS
//------------------------------------------------------------------------------
// Returns Bin_num as debug
int Anaaqra::Debug() {
  return Bin_num;
}



// Returns TH1F of the selected events per rigidity bin
TH1F Anaaqra::RigBinner() {

  // Empty the event histograms
  Events_raw->Reset("ICESM");
  Events_cut->Reset("ICESM");

  // Loop over all Compact entries
  for (int i=0; i<Simp_chain->GetEntries(); i++) {

    // Get entry
    Simp_chain->GetEntry(i);

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

    // Fill histograms
    Events_raw->Fill(Tool->trk_rig);
    if (Crig && Ctrg && Cpar && Ccon && Cbet && Cchi && Cinn && Clay && Cgeo) {
      Events_cut->Fill(Tool->trk_rig);
    }

  }

  // Canvas
  TCanvas* c_Events = new TCanvas("c_Events", "Events per Rigitidy Bin");
  Events_raw->Draw();
  Events_cut->Draw("Same");
  // Styling
  Events_raw->SetLineColor(kRed);
  Events_raw->SetLineWidth(2);
  Events_cut->SetLineColor(kBlue);
  Events_cut->SetLineWidth(2);
  // Axes
  c_Events->SetLogy();
  Events_raw->SetMinimum(1);
  Events_raw->GetXaxis()->SetTitle("R [GV]");
  Events_raw->GetYaxis()->SetTitle("Events");
  // Print
  c_Events->Draw();
  c_Events->Print("./ProtonAnalysis/Events Histogram.png");

  // Return
  return *Events_cut;

}



// Returns TH1F of the exposure time
TH1F Anaaqra::Exposure() {

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

  // Canvas
  TCanvas* c_Exposure = new TCanvas("c_Exposure", "Exposure Time per Rigitidy Bin");
  ExposureTime->Draw();
  // Styling
  ExposureTime->SetLineColor(kGreen);
  ExposureTime->SetLineWidth(2);
  // Axes
  c_Exposure->SetLogy();
  ExposureTime->SetMinimum(1);
  ExposureTime->GetXaxis()->SetTitle("R [GV]");
  ExposureTime->GetYaxis()->SetTitle("Exposure Time [s]");
  // Print
  c_Exposure->Draw();
  c_Exposure->Print("./ProtonAnalysis/Exposure Time Histogram.png");

  // Return
  return *ExposureTime;

}



// Returns TH1F of the geometric Acceptance
TH1F Anaaqra::Acceptance(bool apply_cuts = 0) {

  // Loop over MC root files individually (!)
  int start_num = 1209496744;
  for (int i=0; i<13; i++) { // ??? Make adaptable

    // Import tree
    TChain mc_chain("Compact");
    TChain fi_chain("File");
    mc_chain.Add(Form("../MC Protons/%d.root", start_num + i));
    fi_chain.Add(Form("../MC Protons/%d.root", start_num + i));
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
      if (EventSelectorCompact(MC_comp, "111111110")) {
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

  // Canvas
  TCanvas* c_Acceptance = new TCanvas("c_Acceptance", "Acceptance per Rigitidy Bin");
  AcceptHist->Draw("hist");
  // Styling
  AcceptHist->SetLineColor(kOrange);
  AcceptHist->SetLineWidth(2);
  // Axes
  c_Acceptance->SetLogy();
  AcceptHist->GetXaxis()->SetTitle("R [GV]");
  AcceptHist->GetYaxis()->SetTitle("Acceptance [m^2 sr]");
  // Print
  c_Acceptance->Draw();
  c_Acceptance->Print("./ProtonAnalysis/Acceptance Histogram.png");

  // Return
  return *AcceptHist;

}



// Returns TH1F of the selection efficiency
TH1F Anaaqra::CutEff() {
  TH1F* empty = new TH1F();
  return *empty;
}



// Returns TH1F of the trigger efficiency
const char* Anaaqra::TrigEff() {

  // Temporary histograms
  TH1F* physHist_mc = new TH1F("physHist_mc", "physHist_mc", 32, Bin_edges);
  TH1F* unphHist_mc = new TH1F("unphHist_mc", "unphHist_mc", 32, Bin_edges);
  TH1F* physHist_data = new TH1F("physHist_data", "physHist_data", 32, Bin_edges);
  TH1F* unphHist_data = new TH1F("unphHist_data", "unphHist_data", 32, Bin_edges);

  // Loop over MC entries
  for (int i=0; i<MC_chain->GetEntries(); i++) {

    // Get entry
    MC_chain->GetEntry(i);

    // Check for physical trigger
    bool HasPhysTrig = ((MC_comp->sublvl1&0x3E)!=0)&&((MC_comp->trigpatt&0x2)!=0);
    bool HasUnphTrig = ((MC_comp->sublvl1&0x3E)==0)&&((MC_comp->trigpatt&0x2)!=0);

    // Fill histograms
    if (HasPhysTrig) {
      physHist_mc->Fill(MC_comp->trk_rig[0]);
    }
    if (HasUnphTrig) {
      unphHist_mc->Fill(MC_comp->trk_rig[0]);
    }

  }

  // Loop over data entries
  for (int i=0; i<Simp_chain->GetEntries(); i++) {

    // Get entry
    Simp_chain->GetEntry(i);

    // Check for physical trigger
    bool HasPhysTrig = ((Tool->sublvl1&0x3E)!=0)&&((Tool->trigpatt&0x2)!=0);
    bool HasUnphTrig = ((Tool->sublvl1&0x3E)==0)&&((Tool->trigpatt&0x2)!=0);

    // Fill histograms
    if (HasPhysTrig) {
      physHist_data->Fill(Tool->trk_rig);
    }
    if (HasUnphTrig) {
      unphHist_data->Fill(Tool->trk_rig);
    }

  }

  // Loop over rigidity bins
  for (int i=0; i<Bin_num; i++) {

    TrigEff_MC->SetBinContent(i+1, physHist_mc->GetBinContent(i+1) / (physHist_mc->GetBinContent(i+1) + unphHist_mc->GetBinContent(i+1)));
    TrigEff_data->SetBinContent(i+1, physHist_data->GetBinContent(i+1) / (physHist_data->GetBinContent(i+1) + 100 * unphHist_data->GetBinContent(i+1)));

  }

  return "Assigned to the reference parameters.";

}



// Returns TH1F of the proton rate
TH1F Anaaqra::ProtonRate() {

  // Fill empty necessary histograms
  if (Events_cut->GetEntries() == 0) {
    TH1F Events_cut = RigBinner();
  }
  if (ExposureTime->GetEntries() == 0) {
    TH1F ExposureTime = Exposure();
  }

  // Loop over rigidity bins
  for (int i=0; i<Bin_num; i++) {
    if (ExposureTime->GetBinContent(i+1) == 0) {
      RateHist->SetBinContent(i+1, 0);
    } else {
      RateHist->SetBinContent(i+1, Events_cut->GetBinContent(i+1) / ExposureTime->GetBinContent(i+1) / (2 * Bin_err[i]));
    }
  }

  // Canvas
  TCanvas* c_Rate = new TCanvas("c_Rate", "Proton Rate per Rigitidy Bin");
  RateHist->Draw();
  // Styling
  RateHist->SetLineColor(kMagenta);
  RateHist->SetLineWidth(2);
  // Axes
  c_Rate->SetLogy();
  RateHist->GetXaxis()->SetTitle("R [GV]");
  RateHist->GetYaxis()->SetTitle("Rate [s^-1]");
  // Print
  c_Rate->Draw();
  c_Rate->Print("./ProtonAnalysis/Proton Rate Histogram.png");

  // Return
  return *RateHist;

}



// Returns TH1F of the proton flux
TH1F Anaaqra::ProtonFlux() {

  // Fill empty necessary histograms
  if (Events_cut->GetEntries() == 0) {
    TH1F Events_cut = RigBinner();
  }
  if (ExposureTime->GetEntries() == 0) {
    TH1F ExposureTime = Exposure();
  }
  if (AcceptHist->GetEntries() == 0) {
    TH1F AcceptHist = Acceptance();
  }

  TH1F* empty = new TH1F();
  return *empty;
}
