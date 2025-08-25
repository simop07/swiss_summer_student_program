// To compile in SHELL:
// "g++ lightAnalysis.cpp `root-config --cflags --libs`"

#include <algorithm>
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
#include "TLatex.h"
#include "TLegend.h"
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
  // par[1] = Angle (degrees)
  // par[2] = #sigma
  // par[3] = R_2

  Double_t xVal = x[0];

  Double_t alpha = alphaFunc(xVal, par[1]);

  Double_t expo = par[0] * TMath::Exp(-TMath::Power(alpha, 2) /
                                      (2.0 * TMath::Power(par[2], 2)));

  Double_t cosine = par[3] * TMath::Cos(xVal * TMath::Pi() / 180.0);

  Double_t fitVal = (expo + cosine);

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
  // Define angle vector
  std::vector<int> angles{30, 45, 60};

  // Define useful vector for thicknesses and colours
  std::array<double, 8> thicknesses{0.20, 0.80, 1.55, 2.05,
                                    3.10, 3.60, 4.10, 5.15};
  std::array<int, 8> colours{kBlue,   kRed,      kOrange + 2, kGreen + 2,
                             kViolet, kCyan + 1, kMagenta,    kBlack};
  std::array<int, 8> markers{20, 21, 22, 23, 24, 25, 26, 27};

  // Store probabilities for each angle here
  std::vector<double> prob30T;
  std::vector<double> prob30R;
  std::vector<double> prob45T;
  std::vector<double> prob45R;
  std::vector<double> prob60T;
  std::vector<double> prob60R;

  // Store fit parameters for transmittance
  double a30{};
  double lambda30{};
  double b30{};
  double a45{};
  double lambda45{};
  double b45{};
  double a60{};
  double lambda60{};
  double b60{};

  // Loop on angles
  for (auto angle : angles) {
    // Create file to save canvases
    TFile *reflTransmAnalysis =
        new TFile(Form("./rootFiles/reflTransmAnalyisis%dDegrees.root", angle),
                  "RECREATE");

    // Collect points for fixed angle
    std::vector<Point> points;
    points.reserve(thicknesses.size());
    for (auto thick : thicknesses) {
      points.push_back(lightAnalysis(
          Form("./rootFiles/%dDegrees%.2fmm/wA", angle, thick),
          Form("./rootFiles/%dDegrees/lightAnalysis%dDeg%.2fmm.root", angle,
               angle, thick)));
    }

    // Create refl and transm probability vectors
    std::vector<double> probT;
    std::vector<double> probR;
    for (auto &p : points) {
      probT.push_back(p.x);
      probR.push_back(p.y);
    }

    // Create function for angular correction
    TF1 *fFFunction = new TF1(Form("fFFunction%d", angle), FFunc, -90., 90., 4);

    // Fix angles parameters
    switch (angle) {
      case 30:
        fFFunction->FixParameter(0, 40.59);          // R_1
        fFFunction->FixParameter(1, (double)angle);  // Angle (degrees)
        fFFunction->FixParameter(2, 0.1768);         // #sigma
        fFFunction->FixParameter(3, 48.99);          // R_2
        break;

      case 45:
        fFFunction->FixParameter(0, 90.0);           // R_1
        fFFunction->FixParameter(1, (double)angle);  // Angle (degrees)
        fFFunction->FixParameter(2, -0.1469);        // #sigma
        fFFunction->FixParameter(3, 50.18);          // R_2
        break;

      case 60:
        fFFunction->FixParameter(0, 240);            // R_1
        fFFunction->FixParameter(1, (double)angle);  // Angle (degrees)
        fFFunction->FixParameter(2, 0.1314);         // #sigma
        fFFunction->FixParameter(3, 66.51);          // R_2
        break;
    }

    // Compute ratio between the PMT area and the total area [-90,90]
    double angle1{angle - 25.};
    double angle2{angle + 25.};
    auto geomCorr =
        fFFunction->Integral(angle1, angle2) / fFFunction->Integral(-90., 90.);

    // Print correction factor for 45 degrees
    std::cout << std::fixed << std::setprecision(6)
              << "\nCorrection factor for " << angle << "° (" << angle1 << "°–"
              << angle2 << "°) = " << geomCorr << "\n\n";

    // Apply correction factor to Reflectance
    std::transform(probR.begin(), probR.end(), probR.begin(),
                   [geomCorr](double r) { return r / geomCorr; });

    // Save probabilities into vectors
    switch (angle) {
      case 30:
        std::copy(probT.begin(), probT.end(), std::back_inserter(prob30T));
        std::copy(probR.begin(), probR.end(), std::back_inserter(prob30R));
        break;

      case 45:
        std::copy(probT.begin(), probT.end(), std::back_inserter(prob45T));
        std::copy(probR.begin(), probR.end(), std::back_inserter(prob45R));
        break;

      case 60:
        std::copy(probT.begin(), probT.end(), std::back_inserter(prob60T));
        std::copy(probR.begin(), probR.end(), std::back_inserter(prob60R));
        break;
    }

    // Print table with (Prob_T, Prob_R) and thicknesses
    std::cout << "\n\n*** Results in trigger region ***\n\n";
    std::cout << std::setw(20) << std::left << "Thickness [mm]" << std::setw(20)
              << "Prob_T" << std::setw(20) << "Prob_R" << '\n';
    for (size_t i = 0; i < thicknesses.size(); ++i) {
      std::cout << std::fixed << std::setprecision(3);
      std::cout << std::setw(20) << std::left << thicknesses[i] << std::setw(20)
                << probT[i] << std::setw(20) << probR[i] << '\n';
    }

    // Multigraphs
    TMultiGraph *mg = new TMultiGraph();
    TMultiGraph *mgReflVsTransm = new TMultiGraph();
    std::array<TGraph *, 2> g{
        new TGraph(thicknesses.size(), thicknesses.data(), probT.data()),
        new TGraph(thicknesses.size(), thicknesses.data(), probR.data())};

    std::array<std::string, 2> names{"Transmittance", "Reflectance"};
    for (size_t i = 0; i < names.size(); ++i) {
      g[i]->SetMarkerStyle(markers[i]);
      g[i]->SetMarkerColor(colours[i]);
      g[i]->SetLineColor(colours[i]);
      g[i]->SetTitle(names[i].c_str());
      mg->Add(g[i], "LP");
    }

    // Canvases
    TCanvas *cDeg = new TCanvas(Form("c%dDeg", angle),
                                Form("%d degrees", angle), 1500, 700);
    TCanvas *cReflVsTransm = new TCanvas(Form("cReflVsTransm%d", angle),
                                         "Refl vs transm", 1500, 700);

    // Fit transmittance
    cDeg->cd();
    TF1 *beerLambertTransm =
        new TF1(Form("beerLambert%dTransm", angle), beerLambert, 0.1, 4.2, 3);
    beerLambertTransm->SetLineColor(kRed);
    beerLambertTransm->SetLineWidth(4);
    beerLambertTransm->SetLineStyle(2);
    beerLambertTransm->SetTitle("Transmittance fit");
    beerLambertTransm->SetParNames("A", "#lambda_t", "B");
    beerLambertTransm->SetParameter(0, 0.5);
    beerLambertTransm->SetParameter(1, 350e-3);
    beerLambertTransm->SetParameter(2, 0.0);
    g[0]->Fit(beerLambertTransm, "R");

    // Save transmittance fit parameters
    switch (angle) {
      case 30:
        a30 = beerLambertTransm->GetParameter(0);
        lambda30 = beerLambertTransm->GetParameter(1);
        b30 = beerLambertTransm->GetParameter(2);
        break;

      case 45:
        a45 = beerLambertTransm->GetParameter(0);
        lambda45 = beerLambertTransm->GetParameter(1);
        b45 = beerLambertTransm->GetParameter(2);
        break;

      case 60:
        a60 = beerLambertTransm->GetParameter(0);
        lambda60 = beerLambertTransm->GetParameter(1);
        b60 = beerLambertTransm->GetParameter(2);
        break;
    }

    mg->Draw("ALP");
    mg->SetTitle(Form("%d degrees probability", angle));
    mg->SetName(Form("mg%dDeg", angle));
    mg->GetXaxis()->SetTitle("Thickness [mm]");
    mg->GetYaxis()->SetTitle("Probability");
    cDeg->BuildLegend(.70, .7, .9, .9, "Legend");

    // Transmittance vs Reflectance
    for (size_t i = 0; i < thicknesses.size(); ++i) {
      auto *gr = new TGraph();
      gr->AddPoint(probR[i], probT[i]);
      gr->SetMarkerStyle(markers[i]);
      gr->SetMarkerSize(2);
      gr->SetMarkerColor(colours[i]);
      gr->SetLineColor(colours[i]);
      gr->SetTitle(Form("%.2f", thicknesses[i]));
      mgReflVsTransm->Add(gr, "P");
    }

    cReflVsTransm->cd();
    mgReflVsTransm->Draw("ALP");
    mgReflVsTransm->SetTitle("Transmittance vs Reflectance");
    mgReflVsTransm->SetName("mgReflVsTransm");
    mgReflVsTransm->GetXaxis()->SetTitle("Reflectance prob");
    mgReflVsTransm->GetYaxis()->SetTitle("Transm prob");
    cReflVsTransm->BuildLegend(.70, .7, .9, .9, "Thickness");

    // Save to file
    reflTransmAnalysis->cd();
    cDeg->Write();
    cReflVsTransm->Write();
    reflTransmAnalysis->Close();
  }

  // Cast int into doubles
  std::vector<double> angleVals(angles.begin(), angles.end());

  // Draw a canvas for probabilities as function of angle fixing thicknesses
  TCanvas *cAngle = new TCanvas("cAngle", "Probabilities vs Angle", 1500, 700);
  cAngle->Divide(4, 2);

  // Loop over thicknesses
  for (size_t i = 0; i < thicknesses.size(); ++i) {
    cAngle->cd(i + 1);

    // Extract probabilities for this thickness across different angles
    double probT[3] = {prob30T[i], prob45T[i], prob60T[i]};
    double probR[3] = {prob30R[i], prob45R[i], prob60R[i]};

    // Create graphs
    TGraph *gT = new TGraph(3, angleVals.data(), probT);
    TGraph *gR = new TGraph(3, angleVals.data(), probR);
    gT->SetMarkerStyle(20);
    gT->SetMarkerColor(kBlue);
    gT->SetLineColor(kBlue);
    gT->SetLineWidth(2);
    gT->SetTitle(Form("Thickness %.2f mm", thicknesses[i]));
    gR->SetMarkerStyle(21);
    gR->SetMarkerColor(kRed);
    gR->SetLineColor(kRed);
    gR->SetLineWidth(2);

    // Draw graphs
    cAngle->Update();
    cAngle->cd();
    gT->Draw("ALP");
    gT->GetXaxis()->SetTitle("Angle [^{#circ}]");
    gT->GetYaxis()->SetTitle("Probability");

    gR->Draw("LP SAME");

    // Legend
    TLegend *legend = new TLegend(0.65, 0.75, 0.9, 0.9);
    legend->AddEntry(gT, "Transmittance", "L P");
    legend->AddEntry(gR, "Reflectance", "L P");
    legend->Draw();
  }

  // Print table with (Prob_T, Prob_R) and thicknesses for 30°
  std::cout << "\n\n*** Results in trigger region for 30° ***\n";
  std::cout << std::setw(20) << std::left << "Thickness [mm]" << std::setw(20)
            << "Prob_T" << std::setw(20) << "Prob_R" << '\n';
  for (size_t i = 0; i < thicknesses.size(); ++i) {
    std::cout << std::fixed << std::setprecision(3);
    std::cout << std::setw(20) << std::left << thicknesses[i] << std::setw(20)
              << prob30T[i] << std::setw(20) << prob30R[i] << '\n';
  }

  // Print table with (Prob_T, Prob_R) and thicknesses for 45°
  std::cout << "\n\n*** Results in trigger region for 45° ***\n";
  std::cout << std::setw(20) << std::left << "Thickness [mm]" << std::setw(20)
            << "Prob_T" << std::setw(20) << "Prob_R" << '\n';
  for (size_t i = 0; i < thicknesses.size(); ++i) {
    std::cout << std::fixed << std::setprecision(3);
    std::cout << std::setw(20) << std::left << thicknesses[i] << std::setw(20)
              << prob45T[i] << std::setw(20) << prob45R[i] << '\n';
  }

  // Print table with (Prob_T, Prob_R) and thicknesses for 60°
  std::cout << "\n\n*** Results in trigger region for 60° ***\n";
  std::cout << std::setw(20) << std::left << "Thickness [mm]" << std::setw(20)
            << "Prob_T" << std::setw(20) << "Prob_R" << '\n';
  for (size_t i = 0; i < thicknesses.size(); ++i) {
    std::cout << std::fixed << std::setprecision(3);
    std::cout << std::setw(20) << std::left << thicknesses[i] << std::setw(20)
              << prob60T[i] << std::setw(20) << prob60R[i] << '\n';
  }

  // Transmittance fit results
  std::cout << "\n\n** Fit results for 30° **\n";
  std::cout << "A           = " << a30 << '\n';
  std::cout << "lambda_t    = " << lambda30 << '\n';
  std::cout << "B           = " << b30 << '\n';
  std::cout << "\n\n** Fit results for 45° **\n";
  std::cout << "A           = " << a45 << '\n';
  std::cout << "lambda_t    = " << lambda45 << '\n';
  std::cout << "B           = " << b45 << '\n';
  std::cout << "\n\n** Fit results for 60° **\n";
  std::cout << "A           = " << a60 << '\n';
  std::cout << "lambda_t    = " << lambda60 << '\n';
  std::cout << "B           = " << b60 << '\n';

  // Creating ROOT File
  TFile *fileAngleAnalysis =
      new TFile("./rootFiles/fileAngleAnalysis.root", "RECREATE");

  // Save canvas
  fileAngleAnalysis->cd();
  cAngle->Write();
  fileAngleAnalysis->Close();
}

int main() {
  reflTransm();

  return EXIT_SUCCESS;
}