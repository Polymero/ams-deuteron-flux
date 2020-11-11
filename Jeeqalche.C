// C++ macro's edited for remote ssh to Kapteyn
// Created        23-10-20
// Last Edited    24-10-20

// Include header file(s)
#include <iostream>
#include <string>
#include "Ntp.h"

// CREATE SIMPLIFIED TREE FROM ROOT FILES
void Miiqtoolat(string rootfiles = "local") {

  // Create new objects
  TFile f("../Simp.root", "recreate");
  TTree *T = new TTree("Simp", "Simplified Compact Tree");
  TTree *C = new TTree("RTIInfo", "RTIInfo");
  NtpCompact *Compact = new NtpCompact();
  NtpSHeader *SHeader = new NtpSHeader();
  RTIInfo *RTIInfo = new class RTIInfo();

  // Create chains with pass7 files
  TChain comp_chain("Compact");
  TChain rtii_chain("RTI");
  // Print rootfiles (for visual check)
  cout << "File location: " << rootfiles << endl;
  // Add files to TChain objects
  try {
    if(rootfiles == "kapteyn") {
      comp_chain.Add("/net/dataserver3/data/users/bueno/data/iss/AO/ISS.B1130/pass7/*.root");
      rtii_chain.Add("/net/dataserver3/data/users/bueno/data/iss/AO/ISS.B1130/pass7/*.root");
    } else if(rootfiles == "local") {
      comp_chain.Add("../Runs/*.root");
      rtii_chain.Add("../Runs/*.root");
    } else {
      throw 001;
    }
  }
  catch(Int_t e) {
    cout << "Error Nr. " << e << ": File location not recognised." << endl;
    return 0;
  }
  // Set branch address to compact
  comp_chain.SetBranchAddress("Compact", &Compact);
  comp_chain.SetBranchAddress("SHeader", &SHeader);
  rtii_chain.SetBranchAddress("RTIInfo", &RTIInfo);

  // Initialise parameters
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
  Float_t cf;                   ///< Max geomagnetic cutoff in the field of view (Stoermer|40 degrees|+) [GV]

  // Initialise branches
  T->Branch("status",     &status,      "status/I");
  T->Branch("event",      &event,       "event/I");
  T->Branch("utime",      &utime,       "utime/I");
  T->Branch("trk_q_inn",  &trk_q_inn,   "trk_q_inn/F");
  T->Branch("trk_q_lay",  &trk_q_lay,   "trk_q_lay[9]/F");
  T->Branch("trk_rig",    &trk_rig,     "trk_rig/F");
  T->Branch("trk_chisqn", &trk_chisqn,  "trk_chisqn[2]/F");
  T->Branch("tof_beta",   &tof_beta,    "tof_beta/F");
  T->Branch("rich_beta",  &rich_beta,   "rich_beta/F");
  T->Branch("utime_rti",  &utime_rti,   "utime_rti/I");
  T->Branch("lf",         &lf,          "lf/F");
  T->Branch("cf",         &cf,          "cf/F");

  // Get number of entries
  Int_t nentries = comp_chain.GetEntries();
  // Print nentries (for visual check)
  cout << "Number of Entries: " << nentries << endl;

  // Loop over entries
  Int_t j=0;
  for(Int_t i=0; i < nentries; i++){

    // Get entry
    comp_chain.GetEntry(i);
    // Compact parameters
    status = Compact->status;
    trk_q_inn = Compact->trk_q_inn;
    for(Int_t j=0; j < 9; j++){
      trk_q_lay[j] = Compact->trk_q_lay[j];
    }
    trk_rig = Compact->trk_rig[0];
    trk_chisqn[0] = Compact->trk_chisqn[0][0];
    trk_chisqn[1] = Compact->trk_chisqn[0][1];
    tof_beta = Compact->tof_beta;
    rich_beta = Compact->rich_beta;

    // SHeader parameters
    event = SHeader->event;
    utime = SHeader->utime;

    // Get entry
    rtii_chain.GetEntry(j);
    utime_rti = RTIInfo->utime;
    // Math Compact utime with RTI utime
    while(utime != utime_rti){
      j++;
      rtii_chain.GetEntry(j);
      utime_rti = RTIInfo->utime;
    }
    // RTIInfo parameters
    lf = RTIInfo->lf;
    cf = RTIInfo->cf[0][3][1];

    // Fill tree
    T->Fill();
  }

  // Write tree to file
  T->Write();
  // Print tree (for visual check)
  T->Print();
  T->Scan("tof_beta:utime:utime_rti:lf:cf", "", "colsize=15 precision=10", 24, 40);


  // Get all RTIInfo values too into a seperate tree
  C->Branch("utime_rti",  &utime_rti,   "utime_rti/I");
  C->Branch("lf",         &lf,          "lf/F");
  C->Branch("cf",         &cf,          "cf/F");

  for (Int_t i = 0; i < rtii_chain.GetEntries(); i++) {
    rtii_chain.GetEntry(i);
    utime_rti = RTIInfo->utime;
    lf = RTIInfo->lf;
    cf=RTIInfo->cf[0][3][1];

    // Fill tree
    C->Fill();
  }

  // Write tree to file
  C->Write();
  // Print tree (for visual check)
  C->Print();
}
