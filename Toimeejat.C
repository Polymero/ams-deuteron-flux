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

  cout << "Running Toimeejat()..." << endl;

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
    Tool->rich_select     = Compact->rich_select;

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

  cout << "Toimeejat() has fininshed!\n" << endl;

}



// Returns bool if event passed specified selection of cuts
bool EventSelectorCompact(NtpCompact* comp, const char* cutbit, int atype) {

  bool pass = 1;
  double bin_edges[32+1] = {1.00,1.16,1.33,1.51,1.71,1.92,2.15,2.40,2.67,2.97,
                            3.29,3.64,4.02,4.43,4.88,5.37,5.90,6.47,7.09,7.76,
                            8.48,9.26,10.1,11.0,12.0,13.0,14.1,15.3,16.6,18.0,
                            19.5,21.1,22.8};

  // Boolean cuts (instead of TCut)
  // Protons (base cuts)
  bool Crig = (comp->trk_rig[0] > bin_edges[0])&&(comp->trk_rig[0] <= bin_edges[32]);
  bool Ctrg = ((comp->sublvl1&0x3E)!=0)&&((comp->trigpatt&0x2)!=0);
  bool Cpar = comp->status % 10 == 1;
  bool Cbet = comp->tof_beta > 0.3;
  bool Cchi = (comp->trk_chisqn[0][0] < 10)&&(comp->trk_chisqn[0][1] < 10)&&(comp->trk_chisqn[0][0] > 0)&&(comp->trk_chisqn[0][1] > 0);
  bool Cinn = (comp->trk_q_inn > 0.80)&&(comp->trk_q_inn < 1.30);
  // Adjust return bool according to cutbit
  if (cutbit[0] == '1') {pass &= Crig;}
  if (cutbit[1] == '1') {pass &= Ctrg;}
  if (cutbit[2] == '1') {pass &= Cpar;}
  if (cutbit[3] == '1') {pass &= Cbet;}
  if (cutbit[4] == '1') {pass &= Cchi;}
  if (cutbit[5] == '1') {pass &= Cinn;}
  // Deuterons (additional cuts)
  if (atype == 2){
    bool Cagl = comp->rich_select == 2;
    bool Cnaf = comp->rich_select == 1;
    bool Ccon = std::abs(comp->tof_beta - comp->rich_beta)/comp->tof_beta < 0.05;
    bool Clay = comp->trk_q_lay[0] <= 1.7;
    // Adjust return bool according to cutbit
    if (cutbit[8] == '1') {pass &= Cagl;}
    if (cutbit[8] == '2') {pass &= Cnaf;}
    if (cutbit[9] == '1') {pass &= Ccon;}
    if (cutbit[10] == '1') {pass &= Clay;}
  }

  // Return
  return pass;

}



void MCmeejat(string rootfiles = "kapteyn") {

  cout << "Running MCmeejat()..." << endl;

  // Create new objects
  TFile f("../MCHists.root", "recreate");
  // TTree *P = new TTree("Protons", "Proton MC Histograms");
  // TTree *D = new TTree("Deuterons", "Deuteron MC Histograms");
  // // Initialise branches
  // P->Branch("P_detected", &p_det);
  // P->Branch("P_generated", &p_gen);
  // D->Branch("D_detected", &d_det);
  // D->Branch("D_generated", %d_det);
  // Parameters
  int p_num = 52; //minus last one
  int d_num = 140;
  int p_start = 1209496744;
  int d_start = 74496930;
  double bin_edges[32+1] = {1.00,1.16,1.33,1.51,1.71,1.92,2.15,2.40,2.67,2.97,
                            3.29,3.64,4.02,4.43,4.88,5.37,5.90,6.47,7.09,7.76,
                            8.48,9.26,10.1,11.0,12.0,13.0,14.1,15.3,16.6,18.0,
                            19.5,21.1,22.8};
  // Create histograms
  // Acceptance()
  TH1D *p_det = new TH1D("p_det", "Detected MC Proton Events", 32, bin_edges);
  TH1D *p_cut = new TH1D("p_cut", "Selected MC Proton Events", 32, bin_edges);
  TH1D *p_gen = new TH1D("p_gen", "Generated MC Proton Events", 32, bin_edges);
  TH1D *d_det = new TH1D("d_det", "Detected MC Deuteron Events", 32, bin_edges);
  TH1D *d_cut = new TH1D("d_cut", "Selected MC Deuteron Events", 32, bin_edges);
  TH1D *d_gen = new TH1D("d_gen", "Generated MC Deuteron Events", 32, bin_edges);
  // TrigEff()
  TH1D *p_te_phys = new TH1D("p_te_phys", "Physical Triggers MC Proton Events", 32, bin_edges);
  TH1D *p_te_unph = new TH1D("p_te_unph", "Bias Triggers MC Proton Events", 32, bin_edges);
  TH1D *d_te_phys = new TH1D("d_te_phys", "Physical Triggers MC Deuteron Events", 32, bin_edges);
  TH1D *d_te_unph = new TH1D("d_te_unph", "Bias Triggers MC Deuteron Events", 32, bin_edges);
  // CutEff()
  TH1D *p_ce = new TH1D("p_ce", "Selection Efficiency MC Proton Events", 32, bin_edges);
  TH1D *p_ce_err = new TH1D("p_ce_err", "Selection Efficiency Error MC Proton Events", 32, bin_edges);
  TH1D *d_ce = new TH1D("d_ce", "Selection Efficiency MC Deuteron Events", 32, bin_edges);
  TH1D *d_ce_err = new TH1D("d_ce_err", "Selection Efficiency Error MC Deuteron Events", 32, bin_edges);
  // Individual histograms
  TH1F *p_Cpar_MC = new TH1F("p_Cpar_MC", "Proton Single Particle Cut", 32, bin_edges);
  TH1F *p_Cbet_MC = new TH1F("p_Cbet_MC", "Proton Downward Particle Cut", 32, bin_edges);
  TH1F *p_Cchi_MC = new TH1F("p_Cchi_MC", "Proton Well-constructed Track Cut", 32, bin_edges);
  TH1F *p_Cinn_MC = new TH1F("p_Cinn_MC", "Proton Inner Tracker Charge Cut", 32, bin_edges);
  TH1F *p_Btof_MC = new TH1F("p_Btof_MC", "Proton TOF Base Cut", 32, bin_edges);
  TH1F *p_Btrk_MC = new TH1F("p_Btrk_MC", "Proton Tracker Base Cut", 32, bin_edges);
  TH1F *d_Cpar_MC = new TH1F("d_Cpar_MC", "Deuteron Single Particle Cut", 32, bin_edges);
  TH1F *d_Cbet_MC = new TH1F("d_Cbet_MC", "Deuteron Downward Particle Cut", 32, bin_edges);
  TH1F *d_Cchi_MC = new TH1F("d_Cchi_MC", "Deuteron Well-constructed Track Cut", 32, bin_edges);
  TH1F *d_Cinn_MC = new TH1F("d_Cinn_MC", "Deuteron Inner Tracker Charge Cut", 32, bin_edges);
  TH1F *d_Ccon_MC = new TH1F("d_Ccon_MC", "Deuteron Beta Consistency Cut", 32, bin_edges);
  TH1F *d_Clay_MC = new TH1F("d_Clay_MC", "Deuteron First Tracker Layer Charge Cut", 32, bin_edges);
  TH1F *d_Cagl_MC = new TH1F("d_Cagl_MC", "Deuteron Aerogel RICH Select Cut", 32, bin_edges);
  //TH1F *d_Cnaf_MC = new TH1F("d_Cnaf_MC", "Deuteron NaF RICH Select Cut", 32, bin_edges);
  TH1F *d_Btof_MC = new TH1F("d_Btof_MC", "Deuteron TOF Base Cut", 32, bin_edges);
  TH1F *d_Btrk_MC = new TH1F("d_Btrk_MC", "Deuteron Tracker Base Cut", 32, bin_edges);
  TH1F *d_Brch_MC = new TH1F("d_Brch_MC", "Deuteron RICH Base Cut", 32, bin_edges);

  // loop over MC files individually (for gen hists)
  // Protons
  for (int i=0; i<p_num; i++) {

    // Import tree
    TChain p_fi("File");
    p_fi.Add(Form("/net/dataserver3/data/users/bueno/data/mc/Pr.B1200/%d.root", p_start + i));
    // Create empty class objects
    FileMCInfo *p_fmci = new class FileMCInfo();
    // Set branch addresses
    p_fi.SetBranchAddress("FileMCInfo", &p_fmci);

    // Temporary parameters
    double ngen = 0;
    double rig_min = 0;
    double rig_max = 0;

    // Get FileMCInfo entry for info on MC generation
    p_fi.GetEntry(0);
    // Get file parameters
    ngen = p_fmci->ngen_datacard;
    rig_min = p_fmci->momentum[0];
    rig_max = p_fmci->momentum[1];

    // 1/R generation spectrum
    TF1 *genFlux = new TF1("genFlux", "[0]/(x)", rig_min, rig_max);
    // Normalisation
    genFlux->SetParameter(0, log(rig_max / rig_min));

    // Loop over rigidity bins
    for (int j=0; j<32; j++) {
      // Fraction of spectrum in bin
      double frac = genFlux->Integral(bin_edges[j], bin_edges[j+1]) / genFlux->Integral(rig_min, rig_max);
      // Number of events in fraction
      double nev = frac * ngen;
      // Set bin to current bin count + nev
      p_gen->SetBinContent(j+1, p_gen->GetBinContent(j+1) + nev);
    }

  }

  // Deuterons
  for (int i=0; i<d_num; i++) {

    // Get chain
    TChain d_fi("File");
    d_fi.Add(Form("/net/dataserver3/data/users/bueno/data/mc/D.B1220/%d.root", d_start + i));
    // Create empty class objects
    FileMCInfo *d_fmci = new class FileMCInfo();
    // Set branch addresses
    d_fi.SetBranchAddress("FileMCInfo", &d_fmci);

    // Temporary parameters
    double ngen = 0;
    double rig_min = 0;
    double rig_max = 0;

    // Get FileMCInfo entry for info on MC generation
    d_fi.GetEntry(0);
    // Get file parameters
    ngen = d_fmci->ngen_datacard;
    rig_min = d_fmci->momentum[0];
    rig_max = d_fmci->momentum[1];

    // 1/R generation spectrum
    TF1 *genFlux = new TF1("genFlux", "[0]/(x)", rig_min, rig_max);
    // Normalisation
    genFlux->SetParameter(0, log(rig_max / rig_min));

    // Loop over rigidity bins
    for (int j=0; j<32; j++) {
      // Fraction of spectrum in bin
      double frac = genFlux->Integral(bin_edges[j], bin_edges[j+1]) / genFlux->Integral(rig_min, rig_max);
      // Number of events in fraction
      double nev = frac * ngen;
      // Set bin to current bin count + nev
      d_gen->SetBinContent(j+1, d_gen->GetBinContent(j+1) + nev);
    }

  }

  // Loop over MC Compact entries
  // Protons
  // Get chain
  TChain p_mc("Compact");
  p_mc.Add("/net/dataserver3/data/users/bueno/data/mc/Pr.B1200/*.root");
  // Create empty class objects
  NtpCompact *p_comp = new NtpCompact();
  // Set Branch address
  p_mc.SetBranchAddress("Compact", &p_comp);
  for (int i=0; i<p_mc.GetEntries(); i++) {

    // Get entry
    p_mc.GetEntry(i);

    // Fill Acceptance() cuts
    if (EventSelectorCompact(p_comp, "111111x_111", 1)) {
      p_cut->Fill(p_comp->trk_rig[0]);
    }
    p_det->Fill(p_comp->trk_rig[0]);

    // Trigger booleans
    bool HasPhysTrig_mc = ((p_comp->sublvl1&0x3E)!=0)&&((p_comp->trigpatt&0x2)!=0);
    bool HasUnphTrig_mc = ((p_comp->sublvl1&0x3E)==0)&&((p_comp->trigpatt&0x2)!=0);

    // Fill TrigEff() histograms
    if (EventSelectorCompact(p_comp, "101111x_111", 1)) {
      if (HasPhysTrig_mc) {
        p_te_phys->Fill(p_comp->trk_rig[0]);
      }
      if (HasUnphTrig_mc) {
        p_te_unph->Fill(p_comp->trk_rig[0]);
      }
    }

    // Fill CutEff() histograms
    // TOF charge cut
    bool tof_q = (p_comp->tof_q_lay[0] > 0.8)&&(p_comp->tof_q_lay[0] < 1.5);
    // Fill base histograms
    if (EventSelectorCompact(p_comp, "110011x_111", 1)) {p_Btof_MC->Fill(p_comp->trk_rig[0]);}
    if (EventSelectorCompact(p_comp, "111100x_111", 1) && tof_q) {p_Btrk_MC->Fill(p_comp->trk_rig[0]);}
    // Fill cut histograms
    if (EventSelectorCompact(p_comp, "111111x_111", 1)) {p_Cpar_MC->Fill(p_comp->trk_rig[0]);}
    if (EventSelectorCompact(p_comp, "110111x_111", 1)) {p_Cbet_MC->Fill(p_comp->trk_rig[0]);}
    if (EventSelectorCompact(p_comp, "111100x_111", 1) && tof_q) {p_Cchi_MC->Fill(p_comp->trk_rig[0]);}
    if (EventSelectorCompact(p_comp, "111010x_111", 1) && tof_q) {p_Cinn_MC->Fill(p_comp->trk_rig[0]);}

  }

  // Deuterons
  // Get chain
  TChain d_mc("Compact");
  d_mc.Add("/net/dataserver3/data/users/bueno/data/mc/D.B1220/*.root");
  // Create empty class objects
  NtpCompact *d_comp = new NtpCompact();
  // Set Branch address
  d_mc.SetBranchAddress("Compact", &d_comp);
  for (int i=0; i<d_mc.GetEntries(); i++) {

    // Get entry
    d_mc.GetEntry(i);

    // Fill Acceptance() cuts
    if (EventSelectorCompact(d_comp, "111111x_111", 1)) {
      d_cut->Fill(d_comp->trk_rig[0]);
    }
    d_det->Fill(d_comp->trk_rig[0]);

    // Trigger booleans
    bool HasPhysTrig_mc = ((d_comp->sublvl1&0x3E)!=0)&&((d_comp->trigpatt&0x2)!=0);
    bool HasUnphTrig_mc = ((d_comp->sublvl1&0x3E)==0)&&((d_comp->trigpatt&0x2)!=0);

    // Fill TrigEff() histograms
    if (EventSelectorCompact(d_comp, "101111x_111", 2)) {
      if (HasPhysTrig_mc) {
        d_te_phys->Fill(d_comp->trk_rig[0]);
      }
      if (HasUnphTrig_mc) {
        d_te_unph->Fill(d_comp->trk_rig[0]);
      }
    }

    // Fill CutEff() histograms
    // TOF charge cut
    bool tof_q = (d_comp->tof_q_lay[0] > 0.8)&&(d_comp->tof_q_lay[0] < 1.5);
    // Fill base histograms
    if (EventSelectorCompact(d_comp, "110011x_111", 2)) {d_Btof_MC->Fill(d_comp->trk_rig[0]);}
    if (EventSelectorCompact(d_comp, "111100x_111", 2) && tof_q) {d_Btrk_MC->Fill(d_comp->trk_rig[0]);}
    if (EventSelectorCompact(d_comp, "111111x_001", 2)) {d_Brch_MC->Fill(d_comp->trk_rig[0]);}
    // Fill cut histograms
    if (EventSelectorCompact(d_comp, "111111x_111", 2)) {d_Cpar_MC->Fill(d_comp->trk_rig[0]);}
    if (EventSelectorCompact(d_comp, "110111x_111", 2)) {d_Cbet_MC->Fill(d_comp->trk_rig[0]);}
    if (EventSelectorCompact(d_comp, "111100x_110", 2) && tof_q) {d_Cchi_MC->Fill(d_comp->trk_rig[0]);}
    if (EventSelectorCompact(d_comp, "111010x_110", 2) && tof_q) {d_Cinn_MC->Fill(d_comp->trk_rig[0]);}
    if (EventSelectorCompact(d_comp, "111000x_111", 2) && tof_q) {d_Clay_MC->Fill(d_comp->trk_rig[0]);}
    if (EventSelectorCompact(d_comp, "111111x_011", 2)) {d_Ccon_MC->Fill(d_comp->trk_rig[0]);}
    if (EventSelectorCompact(d_comp, "111111x_101", 2)) {d_Cagl_MC->Fill(d_comp->trk_rig[0]);}
    //if (EventSelectorCompact(d_comp, "111111x_201", 2)) {d_Cnaf_MC->Fill(d_comp->trk_rig[0]);}

    }

    // PROTONS
    // Cut Efficiency MC Error
    for (int i=0; i<32; i++) {
      p_ce_err->SetBinContent(i+1, TMath::Sqrt((1/p_Cpar_MC->GetBinContent(i+1) + 1/p_Btof_MC->GetBinContent(i+1))
                                   + (1/p_Cbet_MC->GetBinContent(i+1) + 1/p_Btof_MC->GetBinContent(i+1))
                                   + (1/p_Cchi_MC->GetBinContent(i+1) + 1/p_Btrk_MC->GetBinContent(i+1))
                                   + (1/p_Cinn_MC->GetBinContent(i+1) + 1/p_Btrk_MC->GetBinContent(i+1))));
    }
    // Divide by corresponding instrument base
    p_Cpar_MC->Divide(p_Btof_MC);
    p_Cbet_MC->Divide(p_Btof_MC);
    p_Cchi_MC->Divide(p_Btrk_MC);
    p_Cinn_MC->Divide(p_Btrk_MC);
    // Cut Efficiency MC
    for (int i=0; i<32; i++) {
      p_ce->SetBinContent(i+1, p_Cpar_MC->GetBinContent(i+1)
                               * p_Cbet_MC->GetBinContent(i+1)
                               * p_Cchi_MC->GetBinContent(i+1)
                               * p_Cinn_MC->GetBinContent(i+1));
    }

    // DEUTERONS
    // Cut Efficiency MC Error
    for (int i=0; i<32; i++) {
      d_ce_err->SetBinContent(i+1, TMath::Sqrt((1/d_Cpar_MC->GetBinContent(i+1) + 1/d_Btof_MC->GetBinContent(i+1))
                                   + (1/d_Cbet_MC->GetBinContent(i+1) + 1/d_Btof_MC->GetBinContent(i+1))
                                   + (1/d_Cchi_MC->GetBinContent(i+1) + 1/d_Btrk_MC->GetBinContent(i+1))
                                   + (1/d_Cinn_MC->GetBinContent(i+1) + 1/d_Btrk_MC->GetBinContent(i+1))
                                   + (1/d_Clay_MC->GetBinContent(i+1) + 1/d_Btrk_MC->GetBinContent(i+1))
                                   + (1/d_Ccon_MC->GetBinContent(i+1) + 1/d_Brch_MC->GetBinContent(i+1))
                                   + (1/d_Cagl_MC->GetBinContent(i+1) + 1/d_Brch_MC->GetBinContent(i+1))));
    }
    // Divide by corresponding instrument base
    d_Cpar_MC->Divide(d_Btof_MC);
    d_Cbet_MC->Divide(d_Btof_MC);
    d_Cchi_MC->Divide(d_Btrk_MC);
    d_Cinn_MC->Divide(d_Btrk_MC);
    d_Clay_MC->Divide(d_Btrk_MC);
    d_Ccon_MC->Divide(d_Brch_MC);
    d_Cagl_MC->Divide(d_Brch_MC);
    //d_Cnaf_MC->Divide(d_Brch_MC);
    // Cut Efficiency MC
    for (int i=0; i<32; i++) {
      d_ce->SetBinContent(i+1, d_Cpar_MC->GetBinContent(i+1)
                               * d_Cbet_MC->GetBinContent(i+1)
                               * d_Cchi_MC->GetBinContent(i+1)
                               * d_Cinn_MC->GetBinContent(i+1)
                               * d_Clay_MC->GetBinContent(i+1)
                               * d_Ccon_MC->GetBinContent(i+1)
                               * d_Cagl_MC->GetBinContent(i+1));
    }

  // Fill file  ( NECESSARY ??? )
  // Acceptance()
  // p_det->Write();
  // p_cut->Write();
  // p_gen->Write();
  // d_det->Write();
  // d_cut->Write();
  // d_gen->Write();
  // // TrigEff()
  // p_te_phys->Write();
  // p_te_unph->Write();
  // d_te_phys->Write();
  // d_te_unph->Write();
  // // CutEff()
  // p_ce->Write();
  // p_ce_err->Write();
  // d_ce->Write();
  // d_ce_err->Write();

  // Write file
  f.Write();
  f.Close();

  cout << "MCmeejat() has finished!\n" << endl;

}
