// C++ macro for rate calculation of cosmic protons
// Created        09-11-20
// Last Edited    14-11-20

// Include header file(s)
#include <iostream>
#include "Header Files/Ntp.h"

void ProtonRate(bool logx = 1, bool logy = 1, bool vcheck = 0) {

  // No canvas titles
  gStyle->SetOptTitle(0);

  // Import trees from Simple.root
  TChain simp_chain("Simp");
  simp_chain.Add("../Simp.root");
  TChain rtii_chain("RTIInfo");
  rtii_chain.Add("../Simp.root");

  Float_t RTIcf;
  Float_t RTIlf;
  rtii_chain.SetBranchAddress("cf", &RTIcf);
  rtii_chain.SetBranchAddress("lf", &RTIlf);

  // Creation of the bins
  Float_t bin_left[32] = {1.00,1.16,1.33,1.51,1.71,1.92,2.15,2.40,2.67,2.97,3.29,3.64,4.02,4.43,4.88,5.37,5.90,6.47,7.09,7.76,8.48,9.26,10.1,11.0,12.0,13.0,14.1,15.3,16.6,18.0,19.5,21.1};
  Float_t bin_right[32] = {1.16,1.33,1.51,1.71,1.92,2.15,2.40,2.67,2.97,3.29,3.64,4.02,4.43,4.88,5.37,5.90,6.47,7.09,7.76,8.48,9.26,10.1,11.0,12.0,13.0,14.1,15.3,16.6,18.0,19.5,21.1,22.8};
  Float_t bin_middle[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Float_t DeltaR[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  for (Int_t i = 0; i < 32; i++) {
    bin_middle[i] = (bin_left[i] + bin_right[i]) / 2;
    DeltaR[i] = bin_right[i] - bin_left[i];
  }

  // Number of events per bin
  Float_t Num[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  // Cuts
  TCut* gen_beta = new TCut("tof_beta > 0");
  TCut* recons = new TCut("trk_chisqn[0] < 10 && trk_chisqn[1] < 10");
  TCut* q_inner = new TCut("trk_q_inn > 0.7 && trk_q_inn < 1.3");
  TCut* q_layer = new TCut("trk_q_lay[0] >= 0 && trk_q_lay[1] >= 0 && trk_q_lay[2] >= 0 && trk_q_lay[3] >= 0 && trk_q_lay[4] >= 0 && trk_q_lay[5] >= 0 && trk_q_lay[6] >= 0 && trk_q_lay[7] >= 0 && trk_q_lay[8] >= 0");
  TCut* rig_geo = new TCut("trk_rig > 1.2 * cf");
  // All-in-one cut
  TCut* master_cut = new TCut(*gen_beta && *recons && *q_inner && *q_layer && *rig_geo);


  // Create canvas
  TCanvas* c1 = new TCanvas("c1", "Rigidity Bin Histogram");

  // Loop over all rigidity bins
  for (Int_t bin = 0; bin < 32; bin++) {

    // Create rigidity slice
    TCut* rigslice = new TCut(Form("trk_rig > %f && trk_rig <= %f", bin_left[bin], bin_right[bin]));
    // Draw rigidity slice histogram
    if (bin == 0) {
      simp_chain.Draw(Form("trk_rig >> hist%d(100, 0, 22)", bin), *master_cut && *rigslice, "");
    } else {
      simp_chain.Draw(Form("trk_rig >> hist%d(100, 0, 22)", bin), *master_cut && *rigslice, "same");
    }

    // Get TH1F object for styling
    TH1F* hist = new TH1F();
    gROOT->GetObject(Form("hist%d", bin), hist);
    // Styling
    // hist->SetLineWidth(3);
    // hist->SetStats(0);
    // Axes
    // hist->SetAxisRange(0, 2000, "Y");
    // hist->GetXaxis()->SetTitle("R [GV]");
    // hist->GetXaxis()->SetTitle("Counts");
    // c1->SetLogx(logx);
    // c1->SetLogy(logy);

    // Number of events
    Num[bin] = hist->GetEntries();

  }

  // Draw canvas
  // c1->Draw();
  // Print canvas
  // c1->Print("./ProtonRate/Rigidity Binning.png");

  // Exposure time as function of rigidity
  Float_t T[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  // Loop over all RTI entries
  for (Int_t i = 0; i < rtii_chain.GetEntries(); i++) {
    rtii_chain.GetEntry(i);
    // Loop over rigidity bins
    for (Int_t j = 0; j < 32; j++) {
      // If the center of the rigidity bin is above the geo-magnetic cut-off, add the livetime
      if (bin_middle[j] > 1.2 * RTIcf) {
        T[j] += RTIlf;
      }
    }
  }

  // Rate calculation
  Float_t Rate[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  // Loop over rigidity bins
  for (Int_t i = 0; i < 32; i++) {
    // Calculate rate
    Float_t rate = 0;
    rate = Num[i] / T[i] / DeltaR[i];
    if (isnan(rate)) {
      Rate[i] = 0;
    } else if (isinf(rate)) {
      Rate[i] = 0;
    } else {
      Rate[i] = rate;
    }

    // Print elements (for visual check)
    if (vcheck == 1) {
      if (i == 0) {
        cout << "R_middle" << "     " << "Bin Counts" << "     " << "Exposure Time" << "     " << "Bin Width" << "      " << "Rate" << endl;
      }
      cout << bin_middle[i] << "     " << Num[i] << "     " << T[i] << "     " << DeltaR[i] << "     " << Rate[i] << endl;
    }
  }

  // Error in Rate (purely from N)
  Float_t ErrRate[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  for (Int_t i = 0; i < 32; i++) {
    Float_t errrate = 0;
    errrate = TMath::Sqrt(Num[i]) / T[i] / DeltaR[i];
    if (isnan(errrate)) {
      ErrRate[i] = 0;
    } else if (isinf(errrate)) {
      ErrRate[i] = 0;
    } else {
      ErrRate[i] = errrate;
    }
  }
  // Error in Rigidity
  Float_t ErrRig[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  for (Int_t i = 0; i < 32; i++) {
    ErrRig[i] = DeltaR[i] / 2;
  }

  // Rate Graph
  TGraphErrors* Rgraph = new TGraphErrors(32, bin_middle, Rate, ErrRig, ErrRate);
  TCanvas* c2 = new TCanvas("c2", "Rate Graph");
  // Styling
  Rgraph->SetMarkerStyle(20);
  Rgraph->SetMarkerSize(1);
  Rgraph->SetMarkerColor(kRed);
  // Axes
  Rgraph->GetXaxis()->SetTitle("R [GV]");
  Rgraph->GetYaxis()->SetTitle("Rate [s^-1]");
  c2->SetLogx(logx);
  c2->SetLogy(logy);
  // Draw
  Rgraph->Draw("AP");
  // Print canvas
  c2->Print("./ProtonRate/Proton Rate.png");

  // Counts Graph
  TGraphErrors* Ngraph = new TGraphErrors(32, bin_middle, Num, ErrRig, 0);
  TCanvas* c3 = new TCanvas("c3", "Count Graph");
  // Styling
  Ngraph->SetMarkerStyle(20);
  Ngraph->SetMarkerSize(1);
  Ngraph->SetMarkerColor(kGreen);
  // Axes
  Ngraph->GetXaxis()->SetTitle("R [GV]");
  Ngraph->GetYaxis()->SetTitle("Bin Counts [GV]");
  c3->SetLogx(logx);
  c3->SetLogy(logy);
  // Draw
  Ngraph->Draw("AP");
  // Print canvas
  c3->Print("./ProtonRate/Bin Counts.png");

  // Livetime Graph
  TGraphErrors* Tgraph = new TGraphErrors(32, bin_middle, T, ErrRig, 0);
  TCanvas* c4 = new TCanvas("c4", "Exposure Time Graph");
  // Styling
  Tgraph->SetMarkerStyle(20);
  Tgraph->SetMarkerSize(1);
  Tgraph->SetMarkerColor(kBlue);
  // Axes
  Tgraph->GetXaxis()->SetTitle("R [GV]");
  Tgraph->GetYaxis()->SetTitle("Exposure Time [s]");
  c4->SetLogx(logx);
  c4->SetLogy(logy);
  // Draw
  Tgraph->Draw("AP");
  // Print canvas
  c4->Print("./ProtonRate/Exposure Time.png");

}
