// C++ macro's edited for remote ssh to Kapteyn
// Created        23-10-20
// Last Edited    30-11-20

// Include header file(s)
#include <iostream>
#include <string>
#include "Header Files/Ntp.h"
#include "Header Files/Simple.h"

// CREATE SIMPLIFIED TREE FROM ROOT FILES
void Toimeejat(string rootfiles = "local") {

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
  MiiqRTI *Mrti = new class MiiqRTI();

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
    } else if(rootfiles == "ssh") {
      comp_chain.Add("venendaal@kapteyn.astro.rug.nl:/net/dataserver3/data/users/bueno/data/iss/AO/ISS.B1130/pass7/*.root");
      rtii_chain.Add("venendaal@kapteyn.astro.rug.nl:/net/dataserver3/data/users/bueno/data/iss/AO/ISS.B1130/pass7/*.root");
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
  // BRANCHES
  //----------------------------------------------------------------------------
  // Initialise branches
  C->Branch("RTI", &Mrti);
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

    // RTI parameters
    Mrti->utime_rti = RTIInfo->utime;
    Mrti->lf = RTIInfo->lf;
    Mrti->cf = RTIInfo->cf[0][3][1];

    // Fill RTI map
    rtimap.insert({Mrti->utime_rti, std::pair<float,float>(Mrti->lf, Mrti->cf)});

    // Fill tree
    C->Fill();

  }

  // Loop over Compact entries
  for(Int_t i=0; i < nentries; i++){

    comp_chain.GetEntry(i);

    // Compact parameters
    Tool->status          = Compact->status;
    Tool->sublvl1         = Compact->sublvl1;
    Tool->trigpatt        = Compact->trigpatt;
    Tool->trk_q_inn       = Compact->trk_q_inn;
    for(Int_t k=0; k < 9; k++){
      Tool->trk_q_lay[k]  = Compact->trk_q_lay[k];
    }
    Tool->trk_rig         = Compact->trk_rig[0];
    Tool->trk_chisqn[0]   = Compact->trk_chisqn[0][0];
    Tool->trk_chisqn[1]   = Compact->trk_chisqn[0][1];
    Tool->tof_beta        = Compact->tof_beta;
    for(Int_t k=0; k < 4; k++){
      Tool->tof_q_lay[k] = Compact->tof_q_lay[k];
    }
    Tool->rich_beta       = Compact->rich_beta;

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
  // Print and scan tree (for visual check)
  C->Print();
  C->Scan("utime_rti:lf:cf", "", "colsize=15 precision=10", 24, 0);

  // Write tree to file
  T->Write();
  // Print and scan tree (for visual check)
  T->Print();
  T->Scan("tof_beta:rich_beta:utime:lf:cf", "", "colsize=15 precision=10", 24, 0);

}
