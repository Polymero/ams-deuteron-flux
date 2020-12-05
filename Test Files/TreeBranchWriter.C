#include <iostream>
#include "Header Files/Ntp.h"
#include "Header Files/Simple.h"

void Rip() {

  TFile f("../Try.root", "recreate");

  TTree *T = new TTree("Simp", "Simplified Compact Tree");

  // PARAMETERS
  int event;
  float tof_beta;
  float rich_beta;

  T->Branch("event", &event, "event/I");
  T->Branch("tof_beta", &tof_beta, "tof_beta/F");
  T->Branch("rich_beta", &rich_beta, "rich_beta/F");

  for (int i =0; i < 20; i++) {

    event = i;
    tof_beta = 473.54734 * i * i - 27 * i - 3719;
    rich_beta = 34.10;

    T->Fill();

  }

  T->Write();

  T->Print();

  T->Scan("event:tof_beta:rich_beta", "", "colsize=15 precision=10", 20, 0);

}



void Rip2() {

  TFile f("../Try.root", "recreate");

  TTree *T = new TTree("Simp", "Simplified Compact Tree");

  Miiqtool *Tool = new Miiqtool();

  T->Branch("Tool", &Tool);

  for (int i = 0; i < 20; i++) {

    Tool->event = i;
    Tool->tof_beta = 473.54734 * i * i - 27 * i - 3719;
    Tool->rich_beta = 34.10;

    cout << Tool->event << endl;

    T->Fill();

  }

  T->Write();

  T->Print();

  T->Scan("Tool.event:tof_beta:rich_beta", "", "colsize=15 precision=10", 20, 0);

}
