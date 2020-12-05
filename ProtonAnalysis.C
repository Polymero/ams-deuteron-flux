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
    int Bin_num = 32;
    double Bin_edges[33] = {1.00,1.16,1.33,1.51,1.71,1.92,2.15,2.40,2.67,2.97,
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
    // Data objects
    TChain *Simp_chain = new TChain("Simp");
    TChain *RTII_chain = new TChain("RTIInfo");
    Miiqtool *Tool = new Miiqtool();
    MiiqRTI *Woi = new MiiqRTI();

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

      // Set branch addresses
      Simp_chain->SetBranchAddress("Simp", &Tool);
      RTII_chain->SetBranchAddress("RTI", &Woi);

      cout << "Class succesfully constructed!" << endl;

    };

    //--------------------------------------------------------------------------
    // LIST OF METHODS
    //--------------------------------------------------------------------------
    // Singular (independent)
    TH1F RigBinner();   // Returns number of selected events as function of rigidity
    TH1F Exposure();    // Returns exposure time (livetime) as function of rigidity
    TH1F Acceptance();  // Returns (geometric) acceptance as function of rigidity
    TH1F CutEff();      // Returns the cut (selection) efficiency as function of rigidity
    TH1F TrigEff();     // Returns the trigger efficiency as function of rigidity
    // Plural (dependent)
    TH1F ProtonRate();  // Returns the proton rate as function of rigidity
    TH1F ProtonFlux();  // Returns the proton flux as function of rigidity
    // Debug
    int Debug();

};

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
TH1F Anaaqra::Acceptance() {
  TH1F* empty = new TH1F();
  return *empty;
}

// Returns TH1F of the selection efficiency
TH1F Anaaqra::CutEff() {
  TH1F* empty = new TH1F();
  return *empty;
}

// Returns TH1F of the trigger efficiency
TH1F Anaaqra::TrigEff() {
  TH1F* empty = new TH1F();
  return *empty;
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
  TH1F* empty = new TH1F();
  return *empty;
}
