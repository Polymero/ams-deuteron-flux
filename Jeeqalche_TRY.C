// C++ macro's edited for remote ssh to Kapteyn
// Created        23-10-20
// Last Edited    10-11-20

// Include header file(s)
#include <iostream>
#include <string>
#include "Header Files/Ntp.h"
#include "Header Files/Simple.h"

// CREATE SIMPLIFIED TREE FROM ROOT FILES
void Miiqtoolat(string rootfiles = "local") {

  // Create new objects
  TFile f("../Simp.root", "recreate");
  TTree *T = new TTree("Simp", "Simplified Compact Tree");
  TTree *C = new TTree("RTIInfo", "RTIInfo");
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
  Float_t cf;                   ///< Max geomagnetic cutoff in the field of view (Stoermer|40 degrees|+) [GV]

  // Initialise branches
  C->Branch("utime_rti",  &utime_rti,   "utime_rti/I");
  C->Branch("lf",         &lf,          "lf/F");
  C->Branch("cf",         &cf,          "cf/F");
  T->Branch("Simp", &Tool);

  // Get number of entries
  Int_t nentries = comp_chain.GetEntries();
  // Print nentries (for visual check)
  cout << "Number of Entries: " << nentries << endl;

  // RTI map
  map<int, std::pair<float,float>> rtimap = map<int, std::pair<float,float>>();

  // Loop over RTIInfo entries
  for (Int_t i = 0; i < rtii_chain.GetEntries(); i++) {

    rtii_chain.GetEntry(i);

    // Get parameters
    utime_rti = RTIInfo->utime;
    lf = RTIInfo->lf;
    cf = RTIInfo->cf[0][3][1];

    // Fill RTI map
    rtimap.insert({utime_rti, std::pair<float,float>(lf, cf)});

    // Fill tree
    C->Fill();
  }

  // Loop over Compact entries
  for(Int_t i=0; i < nentries; i++){

    comp_chain.GetEntry(i);

    // Compact parameters
    Tool->status = Compact->status;
    Tool->trk_q_inn = Compact->trk_q_inn;
    for(Int_t k=0; k < 9; k++){
      Tool->trk_q_lay[k] = Compact->trk_q_lay[k];
    }
    Tool->trk_rig = Compact->trk_rig[0];
    Tool->trk_chisqn[0] = Compact->trk_chisqn[0][0];
    Tool->trk_chisqn[1] = Compact->trk_chisqn[0][1];
    Tool->tof_beta = Compact->tof_beta;
    Tool->rich_beta = Compact->rich_beta;

    // SHeader parameters
    Tool->event = SHeader->event;
    Tool->utime = SHeader->utime;

    // RTIInfo parameters
    Tool->lf = rtimap[SHeader->utime].first;
    Tool->cf = rtimap[SHeader->utime].second;

    // Fill tree
    T->Fill();
  }

  // Write tree to file
  C->Write();
  // Print tree (for visual check)
  C->Print();
  C->Scan("utime_rti:lf:cf", "", "colsize=15 precision=10", 24, 0);

  // Write tree to file
  T->Write();
  // Print tree (for visual check)
  T->Print();
  T->Scan("tof_beta:rich_beta:utime:lf:cf", "", "colsize=15 precision=10", 24, 0);



}
