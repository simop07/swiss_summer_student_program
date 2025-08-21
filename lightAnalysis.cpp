// To compile in SHELL:
// "g++ lightAnalysis.cpp `root-config --cflags --libs`"

#include <array>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"
#include "data.hpp"

void setFitStyle() {
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(10);
  gStyle->SetOptFit(0);  // It was 1111
  gStyle->SetPalette(57);
  gStyle->SetOptTitle(1);
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.9);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.2);
  gStyle->SetTitleX(0.5);
  gStyle->SetTitleY(0.98);
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleXOffset(1.2f);
  gStyle->SetTitleYOffset(1.1f);
  gStyle->SetLineScalePS(1);
  // gStyle->SetPadTopMargin(-9.);
  // gStyle->SetPadRightMargin(-9.);
  // gStyle->SetPadBottomMargin(-9.);
  // gStyle->SetPadLeftMargin(-9.);
  // gStyle->SetTitleW(0.5f);
}

Double_t beerLambert(Double_t *x, Double_t *par) {
  // par[0] = A
  // par[1] = #lambda_t
  // par[2] = B

  Double_t xVal = x[0];
  Double_t fitVal = par[0] * TMath::Exp(-xVal / par[1]) + par[2];
  return fitVal;
}

// Alpha function (peak of gaussian specular distribution)
Double_t alphaFunc(Double_t x, Double_t aDeg) {
  Double_t aRad = aDeg * TMath::Pi() / 180.0;
  Double_t xRad = x * TMath::Pi() / 180.0;

  Double_t val = (TMath::Sin(aRad) * TMath::Sin(xRad)) +
                 (TMath::Cos(aRad) * TMath::Cos(xRad));

  return TMath::ACos(val);
}

// F(x) function for reflectance
Double_t FFunc(Double_t *x, Double_t *par) {
  // par[0] = R_1
  // par[1] = a (degrees)
  // par[2] = #beta
  // par[3] = R_2

  Double_t xVal = x[0];

  Double_t alpha = alphaFunc(xVal, par[1]);

  Double_t expo = par[0] * TMath::Exp(-TMath::Power(alpha, 2) /
                                      (2.0 * TMath::Power(par[2], 2)));

  Double_t cosine = par[3] * TMath::Cos(xVal * TMath::Pi() / 180.0);

  Double_t fitVal = (expo + cosine) / 125.6324;

  return fitVal;
}

Point lightAnalysis(std::string filePath = "./rootFiles/45Degrees3Layer/wA",
                    std::string fileLightAnalysisName =
                        "./rootFiles/lightAnalysis45Deg3Layer.root") {
  // Define useful variables
  int const nFiles{4};
  int const nRegions{4};
  TFile *files[nFiles];
  TTree *trees[nFiles];
  TCanvas *canvases[nFiles];
  PhotonData photondata[nFiles];
  TMultiGraph *mg[nFiles];
  std::string namesF[nFiles] = {"Incident_T", "Incident_R", "Transmitted",
                                "Reflected"};
  std::string namesR[nRegions] = {"PreTrig", "Trig", "PostTrig1", "PostTrig2"};
  std::string fileType{".root"};

  // Creating ROOT File
  TFile *fileLightAnalysis =
      new TFile(fileLightAnalysisName.c_str(), "RECREATE");

  setFitStyle();

  // Loop on files
  for (int i = 0; i < nFiles; ++i) {
    // Define files
    std::string fileName =
        filePath + Form("%s", std::to_string(i).c_str()) + fileType;
    files[i] = new TFile(fileName.c_str(), "READ");

    // Control whether files is open or was not opened
    if (files[i]->IsZombie() || !files[i]) {
      std::cerr << "Impossible to open file " << fileName << std::endl;
    }

    // Define canvases
    canvases[i] = new TCanvas(Form("c%s", namesF[i].c_str()),
                              Form("Pulses %s", namesF[i].c_str()), 1500, 700);

    // Define and draw graphs
    mg[i] = (TMultiGraph *)files[i]->Get("Regions of pulses");
    canvases[i]->cd();
    mg[i]->Draw("ALP");
    mg[i]->SetTitle("Pulses");
    mg[i]->SetName("Regions of pulses");
    mg[i]->GetXaxis()->SetTitle("Time since \"trigger\" [#mus]");
    mg[i]->GetYaxis()->SetTitle("Voltage [mV]");

    // Write canvases on file
    fileLightAnalysis->cd();
    canvases[i]->Write();

    // Define trees
    trees[i] = (TTree *)files[i]->Get("variablesRegion");

    // Collect variables from trees

    // PE counters
    trees[i]->SetBranchAddress("numPreTrigPE",
                               &photondata[i].preTrigger.PECounter);
    trees[i]->SetBranchAddress("numTrigPE", &photondata[i].inTrigger.PECounter);
    trees[i]->SetBranchAddress("numPostTrigPE1",
                               &photondata[i].postTrigger1.PECounter);
    trees[i]->SetBranchAddress("numPostTrigPE2",
                               &photondata[i].postTrigger2.PECounter);

    // PE per pulse
    trees[i]->SetBranchAddress("numPreTrigPEPuls",
                               &photondata[i].preTrigger.PEPulses);
    trees[i]->SetBranchAddress("numTrigPEPuls",
                               &photondata[i].inTrigger.PEPulses);
    trees[i]->SetBranchAddress("numPostTrigPE1Puls",
                               &photondata[i].postTrigger1.PEPulses);
    trees[i]->SetBranchAddress("numPostTrigPE2Puls",
                               &photondata[i].postTrigger2.PEPulses);

    // Time extension per region
    trees[i]->SetBranchAddress("deltaPreTrig",
                               &photondata[i].preTrigger.deltaT);
    trees[i]->SetBranchAddress("deltaTrig", &photondata[i].inTrigger.deltaT);
    trees[i]->SetBranchAddress("deltaPostTrig1",
                               &photondata[i].postTrigger1.deltaT);
    trees[i]->SetBranchAddress("deltaPostTrig2",
                               &photondata[i].postTrigger2.deltaT);

    // Read tree data into variables
    auto const nEntries = trees[i]->GetEntries();
    for (Long64_t j = 0; j < nEntries; ++j) {
      trees[i]->GetEntry(j);
    }
  }

  // Round printing to 10 decimal place
  std::cout << std::fixed << std::setprecision(10);

  // Define reflectance and transmittance
  std::vector<double> incT{};
  std::vector<double> incR{};
  std::vector<double> transm{};
  std::vector<double> refl{};
  incT.reserve(nRegions);
  incR.reserve(nRegions);
  transm.reserve(nRegions);
  refl.reserve(nRegions);

  // Print region properties
  int counter{};
  double rate{};
  double rateCorr{};
  double rateCorrRefl{};

  // Loop on files
  for (PhotonData const &pd : photondata) {
    std::cout << Form("\n *** %s", namesF[counter].c_str())
              << " photon data ***\n";

    // Loop on regions
    for (int i{}; i < nRegions; ++i) {
      switch (i) {
        case 0:
          rate = (pd.preTrigger.PECounter) / (pd.preTrigger.deltaT);
          std::cout << "\nPre trigger region" << '\n';
          std::cout << " Rate          = " << rate << " PE/ns\n";
          std::cout << " PE counter    = " << pd.preTrigger.PECounter
                    << " PE\n";
          std::cout << " PE per pulse  = " << pd.preTrigger.PEPulses
                    << " PE/pulse\n";
          std::cout << " Delta time    = " << pd.preTrigger.deltaT << " ns\n";
          switch (counter) {
            case 0:
              incT.push_back(rate);
              break;

            case 1:
              incR.push_back(rate);
              break;

            case 2:
              transm.push_back(rate);
              break;

            case 3:
              rateCorrRefl = rate - incR[i];
              std::cout << " Rate corr ref = " << rateCorrRefl << " PE/ns\n";
              refl.push_back(rateCorrRefl);
              break;
          }
          break;

        case 1:
          std::cout << "\nTrigger region" << '\n';
          std::cout << " PE counter    = " << pd.inTrigger.PECounter << " PE\n";
          rateCorr = (pd.inTrigger.PECounter - (rate * pd.inTrigger.deltaT)) /
                     pd.inTrigger.deltaT;
          std::cout << " Rate          = "
                    << (pd.inTrigger.PECounter) / (pd.inTrigger.deltaT)
                    << " PE/ns\n";
          std::cout << " Rate correct  = " << rateCorr << " PE/ns\n";
          std::cout << " PE per pulse  = " << pd.inTrigger.PEPulses
                    << " PE/pulse\n";
          std::cout << " Delta time    = " << pd.inTrigger.deltaT << " ns\n";
          switch (counter) {
            case 0:
              incT.push_back(rateCorr);
              break;

            case 1:
              incR.push_back(rateCorr);
              break;

            case 2:
              transm.push_back(rateCorr);
              break;

            case 3:
              rateCorrRefl = rateCorr - incR[i];
              std::cout << " Rate corr ref = " << rateCorrRefl << " PE/ns\n";
              refl.push_back(rateCorrRefl);
              break;
          }
          break;

        case 2:
          std::cout << "\nPost trigger region 1" << '\n';
          std::cout << " PE counter    = " << pd.postTrigger1.PECounter
                    << " PE\n";
          rateCorr =
              (pd.postTrigger1.PECounter - (rate * pd.postTrigger1.deltaT)) /
              pd.postTrigger1.deltaT;
          std::cout << " Rate          = "
                    << (pd.postTrigger1.PECounter) / (pd.postTrigger1.deltaT)
                    << " PE/ns\n";
          std::cout << " Rate correct  = " << rateCorr << " PE/ns\n";
          std::cout << " PE per pulse  = " << pd.postTrigger1.PEPulses
                    << " PE/pulse\n";
          std::cout << " Delta time    = " << pd.postTrigger1.deltaT << " ns\n";
          switch (counter) {
            case 0:
              incT.push_back(rateCorr);
              break;

            case 1:
              incR.push_back(rateCorr);
              break;

            case 2:
              transm.push_back(rateCorr);
              break;

            case 3:
              rateCorrRefl = rateCorr - incR[i];
              std::cout << " Rate corr ref = " << rateCorrRefl << " PE/ns\n";
              refl.push_back(rateCorrRefl);
              break;
          }
          break;

        case 3:
          std::cout << "\nPost trigger region 2" << '\n';
          std::cout << " PE counter    = " << pd.postTrigger2.PECounter
                    << " PE\n";
          rateCorr =
              (pd.postTrigger2.PECounter - (rate * pd.postTrigger2.deltaT)) /
              pd.postTrigger2.deltaT;
          std::cout << " Rate          = "
                    << (pd.postTrigger2.PECounter) / (pd.postTrigger2.deltaT)
                    << " PE/ns\n";
          std::cout << " Rate correct  = " << rateCorr << " PE/ns\n";
          std::cout << " PE per pulse  = " << pd.postTrigger2.PEPulses
                    << " PE/pulse\n";
          std::cout << " Delta time    = " << pd.postTrigger2.deltaT << " ns\n";
          switch (counter) {
            case 0:
              incT.push_back(rateCorr);
              break;

            case 1:
              incR.push_back(rateCorr);
              break;

            case 2:
              transm.push_back(rateCorr);
              break;

            case 3:
              rateCorrRefl = rateCorr - incR[i];
              std::cout << " Rate corr ref = " << rateCorrRefl << " PE/ns\n";
              refl.push_back(rateCorrRefl);
              break;
          }
          break;
      }
    }
    ++counter;
  }

  // Print photon information in each region
  std::string titles[7] = {"Region",         "Inc_T [PE/ns]", "Inc_R [PE/ns]",
                           "Transm [PE/ns]", "Refl [PE/ns]",  "Prob_T",
                           "Prob_R"};
  std::cout << "\n\n" << std::left << std::fixed << std::setprecision(3);
  for (auto const &str : titles) {
    std::cout << std::setw(20) << str;
  }
  std::cout << '\n';
  for (int i{}; i < nRegions; ++i) {
    double probT = transm[i] / incT[i];
    double probR = refl[i] / incT[i];
    std::cout << std::fixed << std::setprecision(3) << std::left;
    std::cout << std::setw(20) << namesR[i] << std::setw(20) << incT[i]
              << std::setw(20) << incR[i] << std::setw(20) << transm[i]
              << std::setw(20) << refl[i] << std::setw(20) << probT
              << std::setw(20) << probR << '\n';
  }
  std::cout << "\n\n";

  // Create point (Prob_T, Prob_R)
  Point p{transm[1] / incT[1], refl[1] / incT[1]};

  return p;
}

// Create plots using points (Prob_T, Prob_R) for different configurations
void reflTransm() {
  // Create file to save canvases
  TFile *reflTransmAnalyisis =
      new TFile("./rootFiles/reflTransmAnalyisis45Degrees.root", "RECREATE");

  // 45 DEGREES CONFIGURATION

  // Collect all points from available ROOT files
  Point p45Deg0p20mm =
      lightAnalysis("./rootFiles/45Degrees0.20mm/wA",
                    "./rootFiles/lightAnalysis45Deg0p20mm.root");
  Point p45Deg0p80mm =
      lightAnalysis("./rootFiles/45Degrees0.80mm/wA",
                    "./rootFiles/lightAnalysis45Deg0p80mm.root");
  Point p45Deg1p55mm =
      lightAnalysis("./rootFiles/45Degrees1.55mm/wA",
                    "./rootFiles/lightAnalysis45Deg1p55mm.root");
  Point p45Deg2p05mm =
      lightAnalysis("./rootFiles/45Degrees2.05mm/wA",
                    "./rootFiles/lightAnalysis45Deg2p05mm.root");
  Point p45Deg3p10mm =
      lightAnalysis("./rootFiles/45Degrees3.10mm/wA",
                    "./rootFiles/lightAnalysis45Deg3p10mm.root");
  Point p45Deg3p60mm =
      lightAnalysis("./rootFiles/45Degrees3.60mm/wA",
                    "./rootFiles/lightAnalysis45Deg3p60mm.root");
  Point p45Deg4p10mm =
      lightAnalysis("./rootFiles/45Degrees4.10mm/wA",
                    "./rootFiles/lightAnalysis45Deg4p10mm.root");
  Point p45Deg5p15mm =
      lightAnalysis("./rootFiles/45Degrees5.15mm/wA",
                    "./rootFiles/lightAnalysis45Deg5p15mm.root");

  // Define useful variables
  std::array<std::string, 2> names{"Transmittance", "Reflectance"};

  // Tune the following chunk if varying the nnumber of measurements
  std::array<double, 8> thicknesses{0.20, 0.80, 1.55, 2.05,
                                    3.10, 3.60, 4.10, 5.15};
  std::array<double, 8> probT{p45Deg0p20mm.x, p45Deg0p80mm.x, p45Deg1p55mm.x,
                              p45Deg2p05mm.x, p45Deg3p10mm.x, p45Deg3p60mm.x,
                              p45Deg4p10mm.x, p45Deg5p15mm.x};
  std::array<double, 8> probR{p45Deg0p20mm.y, p45Deg0p80mm.y, p45Deg1p55mm.y,
                              p45Deg2p05mm.y, p45Deg3p10mm.y, p45Deg3p60mm.y,
                              p45Deg4p10mm.y, p45Deg5p15mm.y};
  std::array<int, 8> colours{kBlue,   kRed,      kOrange + 2, kGreen + 2,
                             kViolet, kCyan + 1, kMagenta,    kBlack};
  std::array<int, 8> markers{20, 21, 22, 23, 24, 25, 26, 27};

  // Print table with (Prob_T, Prob_R) and thicknesses
  std::cout << "\n\n*** Results in trigger region ***\n\n";
  std::cout << std::setw(20) << std::left << "Thickness [mm]" << std::setw(20)
            << "Prob_T" << std::setw(20) << "Prob_R" << '\n';
  for (size_t i = 0; i < thicknesses.size(); ++i) {
    std::cout << std::fixed << std::setprecision(3);
    std::cout << std::setw(20) << std::left << thicknesses[i] << std::setw(20)
              << probT[i] << std::setw(20) << probR[i] << '\n';
  }

  // Create function for angular correction
  TF1 *fFFunction45 = new TF1("fFFunction45", FFunc, 0., 180., 4);

  // Fix parameters value for 45 degrees
  fFFunction45->FixParameter(0, 90.0);     // R_1
  fFFunction45->FixParameter(1, 45.0);     // a (degrees)
  fFFunction45->FixParameter(2, -0.1469);  // #beta
  fFFunction45->FixParameter(3, 50.18);    // R_2

  // Compute ratio between the PMT area and the total area [-90,90]
  double angle1{40.};
  double angle2{50.};
  auto geomCorr = fFFunction45->Integral(angle1, angle2) /
                  fFFunction45->Integral(-90., 90.);

  // Print correction factor for 45 degrees
  std::cout << std::fixed << std::setprecision(6)
            << "\nCorrection factor for 45° (" << angle1 << "°–" << angle2
            << "°) = " << geomCorr << "\n\n";

  // Create graphs for (Prob_T, Prob_R) as function of thickness
  TMultiGraph *mg45Deg = new TMultiGraph();
  TMultiGraph *mgReflVsTransm = new TMultiGraph();
  std::array<TGraph *, 2> g45Deg{
      new TGraph(thicknesses.size(), thicknesses.data(), probT.data()),
      new TGraph(thicknesses.size(), thicknesses.data(), probR.data())};

  // Draw all pulses on multigraph object
  for (size_t i = 0; i < names.size(); ++i) {
    g45Deg[i]->SetMarkerStyle(markers[i]);
    g45Deg[i]->SetMarkerColor(colours[i]);
    g45Deg[i]->SetLineColor(colours[i]);
    g45Deg[i]->SetTitle(names[i].c_str());
    mg45Deg->Add(g45Deg[i], "LP");
  }

  // Create canvases to display plots
  TCanvas *c45Deg = new TCanvas("c45Deg", "45 degrees", 1500, 700);
  TCanvas *cReflVsTransm =
      new TCanvas("cReflVsTransm", "Refl vs transm", 1500, 700);

  // Fit transmittance
  c45Deg->cd();
  TF1 *beerLambert45Transm =
      new TF1("beerLambert45Transm", beerLambert, 0.1, 4.2, 3);
  beerLambert45Transm->SetLineColor(kRed);
  beerLambert45Transm->SetLineWidth(4);
  beerLambert45Transm->SetLineStyle(2);
  beerLambert45Transm->SetTitle("Transmittance fit");
  beerLambert45Transm->SetParNames("A", "#lambda_t", "B");
  beerLambert45Transm->SetParameter(0, 0.5);     // A
  beerLambert45Transm->SetParameter(1, 350e-3);  // #lambda_t
  beerLambert45Transm->SetParameter(2, 0.0);     // B
  g45Deg[0]->Fit(beerLambert45Transm, "R");

  // Fill canvas and draw multigraph
  mg45Deg->Draw("ALP");
  mg45Deg->SetTitle("45 degrees probability");
  mg45Deg->SetName("mg45Deg");
  mg45Deg->GetXaxis()->SetTitle("Thickness [mm]");
  mg45Deg->GetYaxis()->SetTitle("Probability");
  c45Deg->BuildLegend(.70, .7, .9, .9, "Legend");

  // Create transmittance vs reflectance graph
  for (size_t i = 0; i < thicknesses.size(); ++i) {
    auto *g = new TGraph();
    g->AddPoint(probR[i], probT[i]);
    g->SetMarkerStyle(markers[i]);
    g->SetMarkerSize(2);
    g->SetMarkerColor(colours[i]);
    g->SetLineColor(colours[i]);
    g->SetTitle(Form("%.2f", thicknesses[i]));
    mgReflVsTransm->Add(g, "P");
  }

  // Draw transmittance vs reflectance graph
  cReflVsTransm->cd();
  mgReflVsTransm->Draw("ALP");
  mgReflVsTransm->SetTitle("Transmittance vs Reflectance");
  mgReflVsTransm->SetName("mgReflVsTransm");
  mgReflVsTransm->GetXaxis()->SetTitle("Reflectance prob");
  mgReflVsTransm->GetYaxis()->SetTitle("Transm prob");
  cReflVsTransm->BuildLegend(.70, .7, .9, .9, "Thickness");

  // Write everything on file
  reflTransmAnalyisis->cd();
  c45Deg->Write();
  cReflVsTransm->Write();
  reflTransmAnalyisis->Close();
}

int main() {
  reflTransm();

  return EXIT_SUCCESS;
}