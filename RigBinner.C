// C++ macro's edited for remote ssh to Kapteyn
// Created        23-10-20
// Last Edited    24-10-20

// Include header file(s)
#include <iostream>
#include "Ntp.h"

// Create several important
void RigBinner() {

  // Import trees from Simple.root
  TChain simp_chain("Simp");
  simp_chain.Add("../Simp.root");
  TChain rtii_chain("RTI");
  rtii_chain.Add("../Simp.root");

  //simp_chain.SetBranchAddress("Simple", &Simple);
  //rtii_chain.SetBranchAddress("RTIInfo", &RTIInfo);

  // Creation of the bins
  Float_t bin_left[11] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20};
  Float_t bin_right[11] = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22};
  Float_t bin_middle[11];
  for (Int_t i = 0; i < 11; i++) {
    bin_middle[i] = (bin_left[i] + bin_right[i]) / 2;
  }

  // Delta R, i.e. rigidity bin width
  Float_t DeltaR = bin_right - bin_left;

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

  // GetEntries count of histograms


  // Exposure time as function of rigidity

}
