#include <array>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
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

Point lightAnalysis(
    std::string filePath = "./rootFiles/45Degrees3Layer/wA3Layer",
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
              rateCorrRefl = (((rate * pd.preTrigger.deltaT) -
                               incR[i] * pd.preTrigger.deltaT) /
                              pd.preTrigger.deltaT);
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
              rateCorrRefl = (((rateCorr * pd.inTrigger.deltaT) -
                               incR[i] * pd.inTrigger.deltaT) /
                              pd.inTrigger.deltaT);
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
              rateCorrRefl = (((rateCorr * pd.postTrigger1.deltaT) -
                               incR[i] * pd.postTrigger1.deltaT) /
                              pd.postTrigger1.deltaT);
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
              rateCorrRefl = (((rateCorr * pd.postTrigger2.deltaT) -
                               incR[i] * pd.postTrigger2.deltaT) /
                              pd.postTrigger2.deltaT);
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

void reflTransm() {
  // Create points (Prob_T, Prob_R) for different configurations

  // 45 DEGREES

  Point p45Deg1Layer =
      lightAnalysis("./rootFiles/45Degrees1Layer/wA1Layer",
                    "./rootFiles/lightAnalysis45Deg1Layer.root");
  Point p45Deg2Layer =
      lightAnalysis("./rootFiles/45Degrees2Layer/wA2Layer",
                    "./rootFiles/lightAnalysis45Deg2Layer.root");
  Point p45Deg3Layer =
      lightAnalysis("./rootFiles/45Degrees3Layer/wA3Layer",
                    "./rootFiles/lightAnalysis45Deg3Layer.root");

  // Create graphs for (Prob_T, Prob_R) as function of thickness
  std::array<double, 3> thicknesses{1.55, 3.1, 5.14};
  std::array<double, 3> probT{p45Deg1Layer.x, p45Deg2Layer.x, p45Deg3Layer.x};
  std::array<double, 3> probR{p45Deg1Layer.y, p45Deg2Layer.y, p45Deg3Layer.y};

  TMultiGraph *mg45Deg = new TMultiGraph();
  std::array<TGraph *, 2> g45Deg{
      new TGraph(thicknesses.size(), thicknesses.data(), probT.data()),
      new TGraph(thicknesses.size(), thicknesses.data(), probR.data())};

  // Set colours for transmittance
  g45Deg[0]->SetMarkerStyle(20);
  g45Deg[0]->SetMarkerColor(kBlue);
  g45Deg[0]->SetLineColor(kBlue);
  g45Deg[0]->SetTitle("Prob transmittance");

  // Set colours for reflectance
  g45Deg[1]->SetMarkerStyle(21);
  g45Deg[1]->SetMarkerColor(kRed);
  g45Deg[1]->SetLineColor(kRed);
  g45Deg[1]->SetTitle("Prob reflectance");

  // Draw all pulses on multigraph object
  for (auto const &g : g45Deg) {
    mg45Deg->Add(g, "LP");
  }

  // Create canvases to display probabilities
  TCanvas *c45Deg = new TCanvas("c45Deg", "45 degrees", 1500, 700);

  // Fill canvas and draw multigraph
  c45Deg->cd();
  mg45Deg->Draw("ALP");
  mg45Deg->SetTitle("45 degrees probability");
  mg45Deg->SetName("45Deg");
  mg45Deg->GetXaxis()->SetTitle("Thickness [mm]");
  mg45Deg->GetYaxis()->SetTitle("Probability");
  c45Deg->BuildLegend(.70, .7, .9, .9, "Legend");
}

int main() {
  Point p = lightAnalysis();

  return EXIT_SUCCESS;
}