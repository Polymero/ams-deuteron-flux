// Include header file(s)
#include <iostream>
#include <string>
#include "Header Files/Ntp.h"
#include "Header Files/Simple.h"

// CREATE SIMPLIFIED TREE FROM ROOT FILES
void Bruh() {

  // Create new objects
  TFile f("../Simp.root", "recreate");
  TTree *T = new TTree("Simp", "Simplified Compact Tree");
  // Reading objects
  NtpCompact *Compact = new NtpCompact();
  NtpSHeader *SHeader = new NtpSHeader();
  RTIInfo *RTIInfo = new class RTIInfo();
  // Writing objects
  Miiqtool *Tool = new class Miiqtool();

  //----------------------------------------------------------------------------
  // DATA
  //----------------------------------------------------------------------------
  // Create chains with pass7 files
  TChain comp_chain("Compact");
  TChain rtii_chain("RTI");
  // Add files to TChain objects
  comp_chain.Add("/net/dataserver3/data/users/bueno/data/iss/AO/ISS.B1130/pass7/*.root");
  rtii_chain.Add("/net/dataserver3/data/users/bueno/data/iss/AO/ISS.B1130/pass7/*.root");

  // Set branch addresses
  comp_chain.SetBranchAddress("Compact", &Compact);
  comp_chain.SetBranchAddress("SHeader", &SHeader);
  rtii_chain.SetBranchAddress("RTIInfo", &RTIInfo);

  //----------------------------------------------------------------------------
  // PARAMETERS
  //----------------------------------------------------------------------------
  Int_t status;                 ///< nParticle()+nAntiCluster()*10+nBetaH()*100+nTrTrack()*1000+nTrRecHit()*10000+nTrdCluster()*1000000+nTofClusterH()*100000000
  Int_t event;                  ///< Event
  Int_t utime;                  ///< JMDC unix time [s]
  // TRACKER
  Float_t trk_q_inn;            ///< Inner Tracker Charge
  Float_t trk_q_lay[9];         ///< Tracker Layer charge (inverted sign for bad status)
  Float_t trk_rig;              ///< Choutko fit rigidities (FS) [GV]
  Float_t trk_chisqn[2];        ///< Choutko fit normalized Chi2 (FS) (X, Y)
  // TOF
  Float_t tof_beta;             ///< TOF Beta
  // RICH
  Float_t rich_beta;            ///< RICH beta best estimator
  // RTIInfo
  Int_t utime_rti;              ///< JMDC unix time [s]
  Float_t lf;                   ///< Livetime [0,1]
  Float_t cf;                   ///< Max geomagnetic cutoff in the field of view (Stoermer|40 degrees|+) [GV

  // Initialise Single Branch
  T->Branch("Simp", &Tool, "status/I:event/I:utime/I:trk_q_inn/F:trk_q_lay[9]/F:trk_rig/F:trk_chisqn[2]/F:tof_beta/F:rich_beta/F:lf/F:cf/F");

  // Get number of entries
  Int_t nentries = comp_chain.GetEntries();
  // Print nentries (for visual check)
  cout << "Number of Entries: " << nentries << endl;


  // Loop over entries
  for(Int_t i=0; i < nentries; i++){

    comp_chain.GetEntry(i);

    // Compact parameters
    status = Compact->status;
    trk_q_inn = Compact->trk_q_inn;
    for(Int_t k=0; k < 9; k++){
      trk_q_lay[k] = Compact->trk_q_lay[k];
    }
    trk_rig = Compact->trk_rig[0];
    trk_chisqn[0] = Compact->trk_chisqn[0][0];
    trk_chisqn[1] = Compact->trk_chisqn[0][1];
    tof_beta = Compact->tof_beta;
    rich_beta = Compact->rich_beta;

    // SHeader parameters
    event = SHeader->event;
    utime = SHeader->utime;

    // RTIInfo parameters
    lf = 0;
    cf = 0;

    // Fill tree
    T->Fill();
  }

  // Write tree to file
  T->Write();
  // Print tree (for visual check)
  T->Print();
  T->Scan("tof_beta:rich_beta:utime:lf:cf", "", "colsize=15 precision=10", 24, 0);



}
