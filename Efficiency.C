// C++ macro for cut efficiency calculation. PNG output in ./Efficiency
// Created        12-11-20
// Last Edited    14-11-20

// Include header file(s)
#include <iostream>
#include "Header Files/Ntp.h"

void Efficiency(bool logx = 0, bool vcheck = 0) {

  // No canvas titles
  gStyle->SetOptTitle(0);

  // Import trees from Simple.root
  TChain simp_chain("Simp");
  simp_chain.Add("../Simp.root");

  // Cuts
  TCut* Crig = new TCut("trk_rig > 0 && trk_rig <= 22");
  TCut* Cpar = new TCut("status % 10 == 1");
  TCut* Ccon = new TCut("abs(tof_beta-rich_beta)/tof_beta < 0.05");
  TCut* Cbet = new TCut("tof_beta > 0");
  TCut* Cchi = new TCut("trk_chisqn[0] < 10 && trk_chisqn[1] < 10 && trk_chisqn[0] > 0 && trk_chisqn[1] > 0");
  TCut* Cinn = new TCut("trk_q_inn > 0.80 && trk_q_inn < 1.30");
  TCut* Clay = new TCut("trk_q_lay[4]>=0 && trk_q_lay[1]>=0 && trk_q_lay[2]>=0 && trk_q_lay[3]>=0 && trk_q_lay[5]>=0 && trk_q_lay[6]>=0 && trk_q_lay[7]>=0 && trk_q_lay[8]>=0 && trk_q_lay[0]>=0");
  TCut* Cgeo = new TCut("trk_rig > 1.2*cf");

  // Cut Map
  map <class TCut, string> CutMap = {
    {*Cpar, "Combi"},
    {*Ccon, "TOF"},
    {*Cbet, "TOF"},
    {*Cchi, "Tracker"},
    {*Cinn, "Tracker"},
    {*Clay, "Tracker"},
    {*Cgeo, "Tracker"}
  };

  Int_t i = 0;
  // Loop over cuts
  for (auto entry : CutMap) {
    // Get current cut
    TCut current_cut = entry.first;
    TCut base_cut = *Crig;

    // Loop over all cuts again
    for (auto i : CutMap) {
      // Add other cuts of other instuments
      if (i.second != entry.second) {
        base_cut += i.first;
      }

    }

    // Create canvas
    TCanvas* c1 = new TCanvas(Form("c%d", i), "Ratio Histogram");

    // Print Cuts (for visual check)
    if (vcheck == 1) {
      cout << base_cut << endl;
      cout << " " << endl;
      cout << current_cut << endl;
      cout << " " << endl;
      cout << " " << endl;
    }

    // Draw histograms
    simp_chain.Draw(Form("trk_rig >> hist1%d(50, 0, 22)", i), base_cut, "");
    simp_chain.Draw(Form("trk_rig >> hist2%d(50, 0, 22)", i), base_cut && current_cut, "same");

    // Get TH1F objects
    TH1F* hist1 = new TH1F();
    TH1F* hist2 = new TH1F();
    gROOT->GetObject(Form("hist1%d", i), hist1);
    gROOT->GetObject(Form("hist2%d", i), hist2);

    // Draw made histograms (uncomment for bugfixing)
    // hist1->Draw();
    // hist1->SetLineColor(kRed);
    // hist2->Draw("same");
    // hist2->SetLineColor(kBlue);

    // Get histrogram ratio
    TH1F* hist3 = (TH1F*) hist2->Clone();
    hist3->SetName(Form("hist3%d", i));
    hist3->Divide(hist1);
    // Draw
    hist3->Draw("");
    // Styling
    hist3->SetStats(0);
    hist3->SetLineColor(kRed);
    hist3->SetLineWidth(4);
    // Axes
    hist3->GetXaxis()->SetTitle("R [GV]");
    hist3->GetYaxis()->SetTitle("Cut Efficiency");
    hist3->SetAxisRange(0,1,"Y");
    c1->SetLogx(logx);



    // SPLICE
    // Input nodees
    double bin_left[32] = {1.00,1.16,1.33,1.51,1.71,1.92,2.15,2.40,2.67,2.97,3.29,3.64,4.02,4.43,4.88,5.37,5.90,6.47,7.09,7.76,8.48,9.26,10.1,11.0,12.0,13.0,14.1,15.3,16.6,18.0,19.5,21.1};
    double bin_right[32] = {1.16,1.33,1.51,1.71,1.92,2.15,2.40,2.67,2.97,3.29,3.64,4.02,4.43,4.88,5.37,5.90,6.47,7.09,7.76,8.48,9.26,10.1,11.0,12.0,13.0,14.1,15.3,16.6,18.0,19.5,21.1,22.8};
    std::vector<double> nodes;
    for (int j = 0; j < 32; j++) {
      nodes.emplace_back((bin_right[j] + bin_left[j]) / 2);
    }
    double x_min = nodes[0];
    double x_max = nodes[nodes.size()-1];
    // Define TF1 using a lambda function
    // Coordinates of each node are fit parameters, incl derivs at beg and end -> nodes.size()*2+2
    TF1 *f_spline = new TF1(Form("func%d", i), [&](Double_t *x, Double_t *par){
      Double_t xx = x[0];
      // Vectors for node coordinates
      std::vector<double> xn;
      std::vector<double> yn;
      for (int j = 0; j < nodes.size(); j++){
        xn.emplace_back(par[j]);
        yn.emplace_back(par[j+nodes.size()]);
      }
      // First deriv at beg and end
      Double_t b1 = par[nodes.size()*2];
      Double_t e1 = par[nodes.size()*2+1];
      // Define spline with nodes and derivs
      TSpline3 sp3("sp3", &xn[0], &yn[0], nodes.size(), "b1e1", b1, e1);
      return
        sp3.Eval(xx);
    }, x_min, x_max, nodes.size()*2+2);
    for (int j = 0; j < nodes.size(); j++){
      f_spline->FixParameter(j, nodes[j]);
      f_spline->SetParameter(j+nodes.size(), hist3->GetBinContent(hist3->FindBin(nodes[j])));
    }
    f_spline->FixParameter(nodes.size()*2, 0.);
    f_spline->FixParameter(nodes.size()*2+1, 0.);
    // Draw Spline
    f_spline->SetLineColor(kBlack);
    f_spline->Draw("same");



    // Draw Canvas
    c1->Draw();

    // Get Filenames
    const char* print_path;
    if (entry.first == *Cpar) {
      print_path = "./Efficiency/Single Particle Cut.png";
    } else if (entry.first == *Ccon) {
      print_path = "./Efficiency/Beta Consistency Cut.png";
    } else if (entry.first == *Cbet) {
      print_path = "./Efficiency/TOF Beta Cut.png";
    } else if (entry.first == *Cchi) {
      print_path = "./Efficiency/Well Constructed Track Cut.png";
    } else if (entry.first == *Cinn) {
      print_path = "./Efficiency/Inner Tracker Charge Cut.png";
    } else if (entry.first == *Clay) {
      print_path = "./Efficiency/Tracker Layer Charge Cut.png";
    } else if (entry.first == *Cgeo) {
      print_path = "./Efficiency/Geo-magnetic Cut-off Cut.png";
    }

    // Print Canvas
    c1->Print(print_path);
    i++;

  }

}
