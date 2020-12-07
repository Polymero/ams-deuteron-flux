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
    TH1F *PhysHist_mc   = new TH1F("PhysHist_mc", "PhysHist_mc", 32, Bin_edges);
    TH1F *UnphHist_mc   = new TH1F("UnphHist_mc", "UnphHist_mc", 32, Bin_edges);
    TH1F *PhysHist_data = new TH1F("PhysHist_data", "PhysHist_data", 32, Bin_edges);
    TH1F *UnphHist_data = new TH1F("UnphHist_data", "UnphHist_data", 32, Bin_edges);
    TH1F *TrigEff_data  = new TH1F("TrigEff_data", "Trigger Efficiency of Data per Rigidity Bin", 32, Bin_edges);
    TH1F *TrigEff_MC    = new TH1F("TrigEff_MC", "Trigger Efficiency of MC per Rigidity Bin", 32, Bin_edges);
    TH1F *FluxHist      = new TH1F("FluxHist", "Flux per Rigidity Bin", 32, Bin_edges);
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
    // Support functions
    bool EventSelectorCompact(NtpCompact* comp, const char* cutbit);
    bool EventSelectorSimple(Miiqtool* tool, const char* cutbit);
    // Singular (independent)
    void RigBinner();                           // Returns number of selected events as function of rigidity
    void Exposure();                            // Returns exposure time (livetime) as function of rigidity
    void Acceptance(bool apply_cuts = 0);       // Returns (geometric) acceptance as function of rigidity
    void CutEff();                              // Returns the cut (selection) efficiency as function of rigidity
    void TrigEff();                             // Returns the trigger efficiency as function of rigidity
    // Plural (dependent)
    void ProtonRate();                          // Returns the proton rate as function of rigidity
    void ProtonFlux();                          // Returns the proton flux as function of rigidity

};

//------------------------------------------------------------------------------
// ASSIST FUNCTIONS
//------------------------------------------------------------------------------
// Returns bool if event passed specified selection of cuts
bool Anaaqra::EventSelectorCompact(NtpCompact* comp, const char* cutbit) {

  bool pass = 1;

  // Boolean cuts (instead of TCut)
  bool Crig = (comp->trk_rig[0] > Bin_edges[0])&&(comp->trk_rig[0] <= Bin_edges[Bin_num]);
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



// Returns bool if event passed specified selection of cuts
bool Anaaqra::EventSelectorSimple(Miiqtool* tool, const char* cutbit) {

  bool pass = 1;

  // Boolean cuts (instead of TCut)
  bool Crig = (tool->trk_rig > Bin_edges[0])&&(tool->trk_rig <= Bin_edges[Bin_num]);
  bool Ctrg = ((tool->sublvl1&0x3E)!=0)&&((tool->trigpatt&0x2)!=0);
  bool Cpar = tool->status % 10 == 1;
  bool Ccon = std::abs(tool->tof_beta - tool->rich_beta)/tool->tof_beta < 0.05;
  bool Cbet = tool->tof_beta > 0;
  bool Cchi = (tool->trk_chisqn[0] < 10)&&(tool->trk_chisqn[1] < 10)&&(tool->trk_chisqn[0] > 0)&&(tool->trk_chisqn[1] > 0);
  bool Cinn = (tool->trk_q_inn > 0.80)&&(tool->trk_q_inn < 1.30);
  bool Clay = (tool->trk_q_lay[0] >= 0)&&(tool->trk_q_lay[1] >= 0)&&(tool->trk_q_lay[2] >= 0)&&(tool->trk_q_lay[3] >= 0)&&(tool->trk_q_lay[4] >= 0)&&(tool->trk_q_lay[5] >= 0)&&(tool->trk_q_lay[6] >= 0)&&(tool->trk_q_lay[7] >= 0)&&(tool->trk_q_lay[8] >= 0);
  bool Cgeo = tool->trk_rig > 1.2 * tool->cf;

  // Adjust return bool according to cutbit
  if (cutbit[0] == '1') {pass &= Crig;}
  if (cutbit[1] == '1') {pass &= Ctrg;}
  if (cutbit[2] == '1') {pass &= Cpar;}
  if (cutbit[3] == '1') {pass &= Ccon;}
  if (cutbit[4] == '1') {pass &= Cbet;}
  if (cutbit[5] == '1') {pass &= Cchi;}
  if (cutbit[6] == '1') {pass &= Cinn;}
  if (cutbit[7] == '1') {pass &= Clay;}
  if (cutbit[8] == '1') {pass &= Cgeo;}

  // Return
  return pass;

}



//------------------------------------------------------------------------------
// METHOD FUNCTIONS
// Returns TH1F of the selected events per rigidity bin
void Anaaqra::RigBinner() {

  // Empty the event histograms
  Events_raw->Reset("ICESM");
  Events_cut->Reset("ICESM");

  // Loop over all Compact entries
  for (int i=0; i<Simp_chain->GetEntries(); i++) {

    // Get entry
    Simp_chain->GetEntry(i);

    // Fill histograms
    Events_raw->Fill(Tool->trk_rig);
    if (EventSelectorSimple(Tool, "111111111")) {
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

  cout << "RigBinner() has finished!" << endl;

}



// Returns TH1F of the exposure time
void Anaaqra::Exposure() {

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

    cout << "Exposure() has finished!" << endl;

}



// Returns TH1F of the geometric Acceptance
void Anaaqra::Acceptance(bool apply_cuts = 0) {

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

  cout << "Acceptance() has finished!" << endl;

}



// Returns TH1F of the selection efficiency
void Anaaqra::CutEff() {

  // Temporary histograms
  TH1F *Cpar_MC = new TH1F("Cpar_MC", "Single Particle Cut", 32, Bin_edges);
  TH1F *Ccon_MC = new TH1F("Ccon_MC", "Consistent Beta Cut", 32, Bin_edges);
  TH1F *Cbet_MC = new TH1F("Cbet_MC", "Single Particle Cut", 32, Bin_edges);
  TH1F *Cchi_MC = new TH1F("Cchi_MC", "Single Particle Cut", 32, Bin_edges);
  TH1F *Cinn_MC = new TH1F("Cinn_MC", "Single Particle Cut", 32, Bin_edges);
  TH1F *Clay_MC = new TH1F("Clay_MC", "Single Particle Cut", 32, Bin_edges);
  TH1F *Btof_MC = new TH1F("Btof_MC", "TOF Base Cut", 32, Bin_edges);
  TH1F *Btrk_MC = new TH1F("Btrk_MC", "Tracker Base Cut", 32, Bin_edges);


  // Loop over MC entries
  for (int i=0; i<MC_chain->GetEntries(); i++) {

    // Get entry
    MC_chain->GetEntry(i);

    // Fill base histograms
    if (EventSelectorCompact(MC_comp, "10000111")) {Btof_MC->Fill(MC_comp->trk_rig[0]);}
    if (EventSelectorCompact(MC_comp, "10111000")) {Btrk_MC->Fill(MC_comp->trk_rig[0]);}

    // Fill cut histograms
    if (EventSelectorCompact(MC_comp, "10100111")) {Cpar_MC->Fill(MC_comp->trk_rig[0]);}
    if (EventSelectorCompact(MC_comp, "10010111")) {Ccon_MC->Fill(MC_comp->trk_rig[0]);}
    if (EventSelectorCompact(MC_comp, "10001111")) {Cbet_MC->Fill(MC_comp->trk_rig[0]);}
    if (EventSelectorCompact(MC_comp, "10111100")) {Cchi_MC->Fill(MC_comp->trk_rig[0]);}
    if (EventSelectorCompact(MC_comp, "10111010")) {Cinn_MC->Fill(MC_comp->trk_rig[0]);}
    if (EventSelectorCompact(MC_comp, "10111001")) {Clay_MC->Fill(MC_comp->trk_rig[0]);}

  }

  // Divide by corresponding instrument base
  Cpar_MC->Divide(Btof_MC);
  Ccon_MC->Divide(Btof_MC);
  Cbet_MC->Divide(Btof_MC);
  Cchi_MC->Divide(Btrk_MC);
  Cinn_MC->Divide(Btrk_MC);
  Clay_MC->Divide(Btrk_MC);

  // Loop over rigidity bins
  for (int i=0; i<Bin_num; i++) {
    CutEff_MC->SetBinContent(i+1, Cpar_MC->GetBinContent(i+1) * Ccon_MC->GetBinContent(i+1)
                                * Cbet_MC->GetBinContent(i+1) * Cchi_MC->GetBinContent(i+1)
                                * Cinn_MC->GetBinContent(i+1) * Clay_MC->GetBinContent(i+1));
  }

  // Temporary histograms
  TH1F *Cpar_data = new TH1F("Cpar_data", "Single Particle Cut", 32, Bin_edges);
  TH1F *Ccon_data = new TH1F("Ccon_data", "Consistent Beta Cut", 32, Bin_edges);
  TH1F *Cbet_data = new TH1F("Cbet_data", "Single Particle Cut", 32, Bin_edges);
  TH1F *Cchi_data = new TH1F("Cchi_data", "Single Particle Cut", 32, Bin_edges);
  TH1F *Cinn_data = new TH1F("Cinn_data", "Single Particle Cut", 32, Bin_edges);
  TH1F *Clay_data = new TH1F("Clay_data", "Single Particle Cut", 32, Bin_edges);
  TH1F *Cgeo_data = new TH1F("Cgeo_data", "Geo-magnetic Cut", 32, Bin_edges);
  TH1F *Btof_data = new TH1F("Btof_data", "TOF Base Cut", 32, Bin_edges);
  TH1F *Btrk_data = new TH1F("Btrk_data", "Tracker Base Cut", 32, Bin_edges);

  // Loop over data entries
  for (int i=0; i<Simp_chain->GetEntries(); i++) {

    // Get entry
    Simp_chain->GetEntry(i);

    // Fill base histograms
    if (EventSelectorSimple(Tool, "100001111")) {Btof_data->Fill(Tool->trk_rig);}
    if (EventSelectorSimple(Tool, "101110000")) {Btrk_data->Fill(Tool->trk_rig);}

    // Fill cut histograms
    if (EventSelectorSimple(Tool, "101001111")) {Cpar_data->Fill(Tool->trk_rig);}
    if (EventSelectorSimple(Tool, "100101111")) {Ccon_data->Fill(Tool->trk_rig);}
    if (EventSelectorSimple(Tool, "100011111")) {Cbet_data->Fill(Tool->trk_rig);}
    if (EventSelectorSimple(Tool, "101111000")) {Cchi_data->Fill(Tool->trk_rig);}
    if (EventSelectorSimple(Tool, "101110100")) {Cinn_data->Fill(Tool->trk_rig);}
    if (EventSelectorSimple(Tool, "101110010")) {Clay_data->Fill(Tool->trk_rig);}
    if (EventSelectorSimple(Tool, "101110001")) {Cgeo_data->Fill(Tool->trk_rig);}

  }

  // Divide by corresponding instrument base
  Cpar_data->Divide(Btof_data);
  Ccon_data->Divide(Btof_data);
  Cbet_data->Divide(Btof_data);
  Cchi_data->Divide(Btrk_data);
  Cinn_data->Divide(Btrk_data);
  Clay_data->Divide(Btrk_data);
  Cgeo_data->Divide(Btrk_data);

  // Loop over rigidity bins
  for (int i=0; i<Bin_num; i++) {
    CutEff_data->SetBinContent(i+1, Cpar_data->GetBinContent(i+1) * Ccon_data->GetBinContent(i+1)
                                  * Cbet_data->GetBinContent(i+1) * Cchi_data->GetBinContent(i+1)
                                  * Cinn_data->GetBinContent(i+1) * Clay_data->GetBinContent(i+1)
                                  * Cgeo_data->GetBinContent(i+1));
  }

  cout << "CutEff() has finished!" << endl;

}



// Returns TH1F of the trigger efficiency
void Anaaqra::TrigEff() {

  // Loop over MC entries
  for (int i=0; i<MC_chain->GetEntries(); i++) {

    // Get entry
    MC_chain->GetEntry(i);

    // Check for physical trigger
    bool HasPhysTrig = ((MC_comp->sublvl1&0x3E)!=0)&&((MC_comp->trigpatt&0x2)!=0);
    bool HasUnphTrig = ((MC_comp->sublvl1&0x3E)==0)&&((MC_comp->trigpatt&0x2)!=0);

    // Fill histograms
    if (EventSelectorCompact(MC_comp, "10111111")) {
      if (HasPhysTrig) {
        PhysHist_mc->Fill(MC_comp->trk_rig[0]);
      }
      if (HasUnphTrig) {
        UnphHist_mc->Fill(MC_comp->trk_rig[0]);
      }
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
    if (EventSelectorSimple(Tool, "101111111")) {
      if (HasPhysTrig) {
        PhysHist_data->Fill(Tool->trk_rig);
      }
      if (HasUnphTrig) {
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
      TrigEff_data->SetBinContent(i+1, PhysHist_data->GetBinContent(i+1) / (PhysHist_data->GetBinContent(i+1) + 100 * UnphHist_data->GetBinContent(i+1)));
    }

  }

  cout << "TrigEff() has finished!" << endl;

}



// Returns TH1F of the proton rate
void Anaaqra::ProtonRate() {

  // Fill empty necessary histograms
  if (Events_cut->GetEntries() == 0) {
    cout << "Running RigBinner()..." << endl;
    RigBinner();
  }
  if (ExposureTime->GetEntries() == 0) {
    cout << "Running Exposure()..." << endl;
    Exposure();
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

  cout << "ProtonRate() has finished!" << endl;

}



// Returns TH1F of the proton flux
void Anaaqra::ProtonFlux() {

  // Fill empty necessary histograms
  if (Events_cut->GetEntries() == 0) {
    cout << "Running RigBinner()..." << endl;
    RigBinner();
  }
  if (ExposureTime->GetEntries() == 0) {
    cout << "Running Exposure()..." << endl;
    Exposure();
  }
  if (AcceptHist->GetEntries() == 0) {
    cout << "Running Acceptance()..." << endl;
    Acceptance();
  }
  if ((CutEff_MC->GetEntries() == 0) || (CutEff_data->GetEntries() == 0)) {
    cout << "Running CutEff()..." << endl;
    CutEff();
  }
  if ((TrigEff_MC->GetEntries() == 0) || (TrigEff_data->GetEntries() == 0)) {
    cout << "Running TrigEff()..." << endl;
    TrigEff();
  }

  // Loop over rigidity bins
  for (int i=0; i<Bin_num; i++) {
    if (ExposureTime->GetBinContent(i+1) == 0) {
      FluxHist->SetBinContent(i+1, 0);
    } else {
      FluxHist->SetBinContent(i+1, Events_cut->GetBinContent(i+1) / ExposureTime->GetBinContent(i+1)
                                   / 2 / Bin_err[i] / AcceptHist->GetBinContent(i+1)
                                   / CutEff_data->GetBinContent(i+1) * CutEff_MC->GetBinContent(i+1)
                                   / TrigEff_data->GetBinContent(i+1) * TrigEff_MC->GetBinContent(i+1)
                                   * pow(Bin_mid[i], 2.7));
    }
  }

  cout << "ProtonFlux() has finished!" << endl;

}
