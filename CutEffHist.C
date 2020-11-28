// C(++) macro for cut efficiency histograms
// Created        28-11-20
// Last edited    28-11-20

// Include header file(s)
#include <iostream>
#include "Header Files/Ntp.h"
#include "Header Files/Simple.h"

void CutEffHist() {

  // Import trees
  TChain simp_chain("Simp");
  simp_chain.Add("../Simp.root");
  // Create empty class objects
  Miiqtool *Tool = new class Miiqtool();
  // Set branch address
  simp_chain.SetBranchAddress("trk_rig", &Tool);



  // Loop over all entries of the chain
  for (Int_t i = 0; i < simp_chain.GetEntries(); i++) {

    simp_chain.GetEntry(i);

    // Boolean cuts (instead of TCut)
    bool Crig = (Tool->trk_rig > 0)&&(Tool->trk_rig <= 22);

  }

}
