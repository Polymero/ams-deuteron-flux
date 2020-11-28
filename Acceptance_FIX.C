// C++ macro for acceptance calculation. PNG output in ./Acceptance
// Created        22-11-20
// Last Edited    23-11-20

// Include header file(s)
#include <iostream>
#include "Header Files/Ntp.h"

void Acceptance() {

  // Define Rigidity Bins
  double bin_left[32] = {1.00,1.16,1.33,1.51,1.71,1.92,2.15,2.40,2.67,2.97,3.29,3.64,4.02,4.43,4.88,5.37,5.90,6.47,7.09,7.76,8.48,9.26,10.1,11.0,12.0,13.0,14.1,15.3,16.6,18.0,19.5,21.1};
  double bin_right[32] = {1.16,1.33,1.51,1.71,1.92,2.15,2.40,2.67,2.97,3.29,3.64,4.02,4.43,4.88,5.37,5.90,6.47,7.09,7.76,8.48,9.26,10.1,11.0,12.0,13.0,14.1,15.3,16.6,18.0,19.5,21.1,22.8};
  double bin_middle[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double ErrRig[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  for (Int_t i=0; i<32; i++) {
    bin_middle[i] = (bin_right[i] + bin_left[i]) / 2;
    ErrRig[i] = (bin_right[i] - bin_left[i]) / 2;
  }
  double bin_edges[33] = {1.00,1.16,1.33,1.51,1.71,1.92,2.15,2.40,2.67,2.97,3.29,3.64,4.02,4.43,4.88,5.37,5.90,6.47,7.09,7.76,8.48,9.26,10.1,11.0,12.0,13.0,14.1,15.3,16.6,18.0,19.5,21.1,22.8};

  // Define histogram
  TH1F* NgenHist = new TH1F("NgenHist", "Generated Events per Rigidity Bin", 32, bin_edges);
  double total_gen = 0;
  // Create canvas
  TCanvas* c1 = new TCanvas("Ngen", "Generated Events per Rigidity Bin");

  // Loop over MC root files to fill Ngen histogram
  int start_num = 1209496744;
  for (Int_t i=0; i<13; i++) {

    // Import tree
    TChain mc_chain("Compact");
    TChain fi_chain("File");
    mc_chain.Add(Form("../MC Protons/%d.root", start_num + i));
    fi_chain.Add(Form("../MC Protons/%d.root", start_num + i));

    // Create empty class objects
    NtpCompact *Compact = new class NtpCompact();
    FileMCInfo *FMCI = new class FileMCInfo();
    // Set Branch address
    mc_chain.SetBranchAddress("Compact", &Compact);
    fi_chain.SetBranchAddress("FileMCInfo", &FMCI);

    // Initialise parameters
    double ngen = 0;
    double rig_min = 0;
    double rig_max = 0;
    // Get FileMCInfo entry for information about MC generation
    fi_chain.GetEntry(0);
    ngen = FMCI->ngen_datacard;
    rig_min = FMCI->momentum[0];
    rig_max = FMCI->momentum[1];

    // Add ngen to total_gen for check
    total_gen += ngen;

    // Generated 1/R Spectrum
    TF1 *genFlux = new TF1("genFlux", "[0]/(x)", rig_min, rig_max);
    // Normalisation
    genFlux->SetParameter(0, log(rig_max / rig_min));

    // Loop over bins
    double check = 0;
    for (Int_t j=0; j<32; j++) {
      // Fraction of spectrum in bin
      double frac = genFlux->Integral(bin_left[j], bin_right[j])/genFlux->Integral(rig_min, rig_max);
      check += frac;
      // Number of events in fraction
      double Nev = frac * ngen;
      // Set bin to current count + Nev
      NgenHist->SetBinContent(j+1, NgenHist->GetBinContent(j+1) + Nev);

    }
    cout << "Total fraction in R=(1.0, 22.8) of " << start_num+i << ".root: " << check << endl;
  }

  // Sum over bins for check
  double total_bin = 0;
  for (Int_t i=0; i<32; i++) {
    total_bin += NgenHist->GetBinContent(i);
  }

  // Draw histogram
  NgenHist->Draw();
  NgenHist->SetStats(0);
  c1->SetLogx();
  c1->SetLogy();
  // Draw canvas
  c1->Draw();
  // Print canvas
  c1->Print("./Acceptance/Acceptance Histograms.png");

  cout << "Total Generated Events: " << total_gen << endl;
  cout << "Total Events in R=(1.0, 22.8): " << total_bin << endl;



  // Define histogram
  TH1F* NdetHist = new TH1F("NdetHist", "Detected Events per Rigidity Bin", 32, bin_edges);
  // Create canvas
  TCanvas* c2 = new TCanvas("Ndet", "Detected Events per Rigidity Bin");

  // Import tree
  TChain ev_chain("Compact");
  ev_chain.Add("../MC Protons/*.root");
  // Create empty class objects
  NtpCompact *EvComp = new class NtpCompact();
  // Set Branch address
  ev_chain.SetBranchAddress("Compact", &EvComp);

  // Fill histogram
  for (Int_t i=0; i<ev_chain.GetEntries(); i++) {
    ev_chain.GetEntry(i);
    NdetHist->Fill(EvComp->trk_rig[0]);
  }

  // Draw histogram
  NdetHist->Draw();
  NdetHist->SetStats(0);
  c2->SetLogx();
  c2->SetLogy();
  // Draw canvas
  c2->Draw();



  // Define histogram
  TH1F* AccHist = (TH1F*) NdetHist->Clone();
  AccHist->SetName("AccHist");
  AccHist->SetTitle("Acceptance per Rigidity Bin");
  AccHist->Divide(NgenHist);
  // Create canvas
  TCanvas* c3 = new TCanvas("Acc", "Acceptance per Rigidity Bin");

  // Scale histogram
  double gen_vol = TMath::Pi() * 3.9 * 3.9;
  AccHist->Scale(gen_vol);
  cout << "MC Generation Volume: " << gen_vol << " m2 sr" << endl;

  // Draw histogram
  AccHist->Draw("hist");
  AccHist->SetStats(0);
  c3->SetLogx();
  // Draw canvas
  c3->Draw();



  // Convert TH1F to TGraph
  cout << "Acceptance: ";
  double Acc[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  for (Int_t i=0; i<32; i++) {
    Acc[i] = AccHist->GetBinContent(i+1);
    cout << Acc[i] << " ";
  }
  cout << endl;
  // Create TGraph
  TGraphErrors* Agraph = new TGraphErrors(32, bin_middle, Acc, ErrRig, 0);
  TCanvas* c4 = new TCanvas("Acc2", "Acceptance TGraph");
  // Styling
  Agraph->SetMarkerStyle(20);
  Agraph->SetMarkerSize(1);
  Agraph->SetMarkerColor(kBlue);
  // Axes
  Agraph->GetXaxis()->SetTitle("R [GV]");
  Agraph->GetYaxis()->SetTitle("Acceptance [m^2 sr]");
  c4->SetLogx();
  // Draw
  Agraph->Draw("AP");

  // Spline
  // Input nodees
  std::vector<double> nodes;
  for (int k = 0; k < 32; k = k+4) {
    nodes.emplace_back(bin_middle[k]);
  }
  nodes.emplace_back(22.08);
  double x_min = nodes[0];
  double x_max = nodes[nodes.size()-1];
  // Define TF1 using a lambda function
  // Coordinates of each node are fit parameters, incl derivs at beg and end -> nodes.size()*2+2
  TF1 *f_spline = new TF1("f_spline", [&](Double_t *x, Double_t *par){
    Double_t xx = x[0];
    // Vectors for node coordinates
    std::vector<double> xn;
    std::vector<double> yn;
    for (int k = 0; k < nodes.size(); k++){
      xn.emplace_back(par[k]);
      yn.emplace_back(par[k+nodes.size()]);
    }
    // First deriv at beg and end
    Double_t b1 = par[nodes.size()*2];
    Double_t e1 = par[nodes.size()*2+1];
    // Define spline with nodes and derivs
    TSpline3 sp3("sp3", &xn[0], &yn[0], nodes.size(), "b1e1", b1, e1);
    return
      sp3.Eval(xx);
  }, x_min, x_max, nodes.size()*2+2);
  for (int k = 0; k < nodes.size(); k++){
    f_spline->FixParameter(k, nodes[k]);
    f_spline->SetParameter(k+nodes.size(), AccHist->GetBinContent(AccHist->FindBin(nodes[k])));
  }
  f_spline->FixParameter(nodes.size()*2, 0.);
  f_spline->FixParameter(nodes.size()*2+1, 0.);
  // Draw Spline
  f_spline->SetLineColor(kBlack);
  f_spline->Draw("same");

  // Print canvas
  c4->Print("./Acceptance/Acceptance TGraph.png");


}
