// C++ macro's edited for remote ssh to Kapteyn
// Created        23-10-20
// Last Edited    11-11-20

// Include header file(s)
#include <iostream>
#include "Ntp.h"
#include "Simple.h"

// Create several important
void RigBinner() {

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
  Float_t bin_left[11] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20};
  Float_t bin_right[11] = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22};
  Float_t bin_middle[11];
  Float_t DeltaR[11];
  for (Int_t i = 0; i < 11; i++) {
    bin_middle[i] = (bin_left[i] + bin_right[i]) / 2;
    DeltaR[i] = bin_right[i] - bin_left[i];
  }

  // Number of events per bin
  // Cuts
  TCut* gen_beta = new TCut("tof_beta > 0");
  TCut* recons = new TCut("trk_chisqn[0] < 10 && trk_chisqn[1] < 10");
  TCut* q_inner = new TCut("trk_q_inn > 0.7 && trk_q_inn < 1.3");
  TCut* q_layer = new TCut("trk_q_lay[0] >= 0 && trk_q_lay[1] >= 0 && trk_q_lay[2] >= 0 && trk_q_lay[3] >= 0 && trk_q_lay[4] >= 0 && trk_q_lay[5] >= 0 && trk_q_lay[6] >= 0 && trk_q_lay[7] >= 0 && trk_q_lay[8] >= 0");
  TCut* rig_geo = new TCut("trk_rig > 1.2 * cf");

  TCut* master_cut = new TCut(*gen_beta && *recons && *q_inner && *q_layer && *rig_geo);

  // Slices
  TCut* slice_1 = new TCut("trk_rig > 0 && trk_rig <= 2");
  TCut* slice_2 = new TCut("trk_rig > 2 && trk_rig <= 4");
  TCut* slice_3 = new TCut("trk_rig > 4 && trk_rig <= 6");
  TCut* slice_4 = new TCut("trk_rig > 6 && trk_rig <= 8");
  TCut* slice_5 = new TCut("trk_rig > 8 && trk_rig <= 10");
  TCut* slice_6 = new TCut("trk_rig > 10 && trk_rig <= 12");
  TCut* slice_7 = new TCut("trk_rig > 12 && trk_rig <= 14");
  TCut* slice_8 = new TCut("trk_rig > 14 && trk_rig <= 16");
  TCut* slice_9 = new TCut("trk_rig > 16 && trk_rig <= 18");
  TCut* slice_10 = new TCut("trk_rig > 18 && trk_rig <= 20");
  TCut* slice_11 = new TCut("trk_rig > 20 && trk_rig <= 22");

  // Draw histograms
  TCanvas* c1 = new TCanvas();
  simp_chain.Draw("trk_rig >> hist_1(100, 0, 22)", *master_cut && *slice_1, "");
  simp_chain.Draw("trk_rig >> hist_2(100)", *master_cut && *slice_2, "same");
  simp_chain.Draw("trk_rig >> hist_3(100)", *master_cut && *slice_3, "same");
  simp_chain.Draw("trk_rig >> hist_4(100)", *master_cut && *slice_4, "same");
  simp_chain.Draw("trk_rig >> hist_5(100)", *master_cut && *slice_5, "same");
  simp_chain.Draw("trk_rig >> hist_6(100)", *master_cut && *slice_6, "same");
  simp_chain.Draw("trk_rig >> hist_7(100)", *master_cut && *slice_7, "same");
  simp_chain.Draw("trk_rig >> hist_8(100)", *master_cut && *slice_8, "same");
  simp_chain.Draw("trk_rig >> hist_9(100)", *master_cut && *slice_9, "same");
  simp_chain.Draw("trk_rig >> hist_10(100)", *master_cut && *slice_10, "same");
  simp_chain.Draw("trk_rig >> hist_11(100)", *master_cut && *slice_11, "same");

  TH1F* hist_1 = new TH1F(); gROOT->GetObject("hist_1", hist_1);
  TH1F* hist_2 = new TH1F(); gROOT->GetObject("hist_2", hist_2);
  TH1F* hist_3 = new TH1F(); gROOT->GetObject("hist_3", hist_3);
  TH1F* hist_4 = new TH1F(); gROOT->GetObject("hist_4", hist_4);
  TH1F* hist_5 = new TH1F(); gROOT->GetObject("hist_5", hist_5);
  TH1F* hist_6 = new TH1F(); gROOT->GetObject("hist_6", hist_6);
  TH1F* hist_7 = new TH1F(); gROOT->GetObject("hist_7", hist_7);
  TH1F* hist_8 = new TH1F(); gROOT->GetObject("hist_8", hist_8);
  TH1F* hist_9 = new TH1F(); gROOT->GetObject("hist_9", hist_9);
  TH1F* hist_10 = new TH1F(); gROOT->GetObject("hist_10", hist_10);
  TH1F* hist_11 = new TH1F(); gROOT->GetObject("hist_11", hist_11);

  hist_1->SetAxisRange(0, 2000, "Y");
  hist_1->SetStats(0);

  hist_1->SetLineColor(kBlack); hist_1->SetLineWidth(3);
  hist_2->SetLineColor(kGray); hist_2->SetLineWidth(3);
  hist_3->SetLineColor(kRed); hist_3->SetLineWidth(3);
  hist_4->SetLineColor(kOrange); hist_4->SetLineWidth(3);
  hist_5->SetLineColor(kGreen); hist_5->SetLineWidth(3);
  hist_6->SetLineColor(kCyan); hist_6->SetLineWidth(3);
  hist_7->SetLineColor(kBlue); hist_7->SetLineWidth(3);
  hist_8->SetLineColor(kViolet); hist_8->SetLineWidth(3);
  hist_9->SetLineColor(kMagenta); hist_9->SetLineWidth(3);
  hist_10->SetLineColor(kGray); hist_10->SetLineWidth(3);
  hist_11->SetLineColor(kBlack); hist_11->SetLineWidth(3);

  c1->Draw();

  // GetEntries count of histograms
  Int_t N[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  N[0] = hist_1->GetEntries();
  N[1] = hist_2->GetEntries();
  N[2] = hist_3->GetEntries();
  N[3] = hist_4->GetEntries();
  N[4] = hist_5->GetEntries();
  N[5] = hist_6->GetEntries();
  N[6] = hist_7->GetEntries();
  N[7] = hist_8->GetEntries();
  N[8] = hist_9->GetEntries();
  N[9] = hist_10->GetEntries();
  N[10] = hist_11->GetEntries();

  // Exposure time as function of rigidity
  Float_t T[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  for (Int_t i = 0; i < rtii_chain.GetEntries(); i++) {
    rtii_chain.GetEntry(i);
    for (Int_t j = 0; j < 11; j++) {
      if (bin_middle[j] > 1.2 * RTIcf) {
        T[j] += RTIlf;
      }
    }
  }

  // Proton rate
  Float_t Rate[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  for (Int_t i = 0; i < 11; i++) {
    Rate[i] = N[i] / T[i] / DeltaR[i];

    cout << N[i] << "     " << T[i] << "     " << DeltaR[i] << "     " << Rate[i] << endl;

  }

}
