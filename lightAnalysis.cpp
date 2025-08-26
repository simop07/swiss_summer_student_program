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
#include "TH1F.h"
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
  // par[4] = Norm factor

  Double_t xVal = x[0];

  Double_t alpha = alphaFunc(xVal, par[1]);

  Double_t expo = par[0] * TMath::Exp(-TMath::Power(alpha, 2) /
                                      (2.0 * TMath::Power(par[2], 2)));

  Double_t cosine = par[3] * TMath::Cos(xVal * TMath::Pi() / 180.0);

  Double_t fitVal = (expo + cosine) / par[4];

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
    canvases[i]->Close();

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

  // Colour array
  int colours[8] = {kRed - 10, kRed - 9, kRed - 7, kRed - 4,
                    kRed,      kRed + 1, kRed + 2, kRed + 3};

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

    // Define angles for shaded region
    double angle1 = angle - 25.;
    double angle2 = angle + 25.;

    // Create canvas
    TCanvas *cFunc = new TCanvas("cFunc", "Angular Correction", 800, 600);

    // Create a frame for axes, titles, and range
    TH1F *frame =
        new TH1F("frame", Form("Correction function for %d^{#circ}", angle),
                 100, -90, 90);
    frame->SetXTitle("Angle [deg]");
    frame->SetYTitle("Arbitrary units");
    frame->SetTitle(Form("Correction function for %d^{#circ}", angle));
    frame->SetMinimum(0);
    frame->SetMaximum(1);

    // Create function for angular correction
    TF1 *fFFunction = new TF1(Form("fFFunction%d", angle), FFunc, -90., 90., 5);

    // Fix parameters depending on angle
    switch (angle) {
      case 30:
        fFFunction->FixParameter(0, 40.59);
        fFFunction->FixParameter(1, (double)angle);
        fFFunction->FixParameter(2, 0.1768);
        fFFunction->FixParameter(3, 48.99);
        fFFunction->FixParameter(4, 83.24086);
        break;
      case 45:
        fFFunction->FixParameter(0, 90.0);
        fFFunction->FixParameter(1, (double)angle);
        fFFunction->FixParameter(2, -0.1469);
        fFFunction->FixParameter(3, 50.18);
        fFFunction->FixParameter(4, 125.6324);
        break;
      case 60:
        fFFunction->FixParameter(0, 240);
        fFFunction->FixParameter(1, (double)angle);
        fFFunction->FixParameter(2, 0.1314);
        fFFunction->FixParameter(3, 66.51);
        fFFunction->FixParameter(4, 273.37408);
        break;
    }

    // Draw main function
    fFFunction->SetLineColor(kBlue);
    fFFunction->SetLineWidth(2);
    cFunc->cd();

    // Create shaded area for angle1-angle2
    TF1 *f1Range = (TF1 *)fFFunction->Clone();
    f1Range->SetRange(angle1, angle2);
    f1Range->SetFillColor(kRed);
    f1Range->SetFillStyle(1001);
    cFunc->cd();
    frame->Draw();
    frame->Draw("SAME AXIS");
    fFFunction->Draw("SAME");
    f1Range->Draw("SAME FC");

    // Add legend
    TLegend *leg = new TLegend(.70, .7, .9, .9, "Legend");
    leg->AddEntry(f1Range, "Shaded region", "F");
    leg->Draw();

    // Add text labels for angle1 and angle2
    TLatex *text1 = new TLatex(angle1, fFFunction->Eval(angle1) + 0.05,
                               Form("#color[1]{%.0f^{#circ}}", angle1));
    text1->SetTextAlign(22);  // Center alignment
    text1->Draw();

    TLatex *text2 = new TLatex(angle2, fFFunction->Eval(angle2) + 0.05,
                               Form("#color[1]{%.0f^{#circ}}", angle2));
    text2->SetTextAlign(22);
    text2->Draw();

    // Update canvas
    cFunc->Update();

    // Compute ratio between the PMT area and the total area [-90,90]
    double geomCorr =
        fFFunction->Integral(angle1, angle2) / fFFunction->Integral(-90., 90.);
    std::cout << std::fixed << std::setprecision(6)
              << "\nCorrection factor for " << angle << "° (" << angle1 << "°-"
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
      g[i]->SetMarkerStyle(20);
      g[i]->SetTitle(names[i].c_str());
      mg->Add(g[i], "LP");
    }

    // Insert colours
    g[0]->SetMarkerColor(kBlue);
    g[0]->SetLineColor(kBlue);
    g[1]->SetMarkerColor(kRed);
    g[1]->SetLineColor(kRed);

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

    setFitStyle();

    // Transmittance vs Reflectance
    for (size_t i = 0; i < thicknesses.size(); ++i) {
      auto *gr = new TGraph();
      gr->AddPoint(probR[i], probT[i]);
      gr->SetMarkerStyle(20);
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
    cFunc->Write();
    cFunc->Close();
    cDeg->Write();
    cDeg->Close();
    cReflVsTransm->Write();
    cReflVsTransm->Close();
    reflTransmAnalysis->Close();
  }

  // Cast int into doubles
  std::vector<double> angleVals(angles.begin(), angles.end());

  // Draw a canvas for probabilities as function of angle fixing thicknesses
  TCanvas *cAngle = new TCanvas("cAngle", "Probabilities vs angle", 1500, 700);
  cAngle->Divide(4, 2);

  // Create vectors for graphs
  std::vector<TGraph *> gTs;
  std::vector<TGraph *> gRs;

  // Loop over thicknesses
  for (size_t i = 0; i < thicknesses.size(); ++i) {
    // Extract probabilities for this thickness across different angles
    double probT[3] = {prob30T[i], prob45T[i], prob60T[i]};
    double probR[3] = {prob30R[i], prob45R[i], prob60R[i]};

    // Push graphs into vectors
    gTs.push_back(new TGraph(3, angleVals.data(), probT));
    gRs.push_back(new TGraph(3, angleVals.data(), probR));
  }

  // Draw graphs in vectors
  for (size_t i = 0; i < gTs.size(); ++i) {
    cAngle->cd(i + 1);
    gTs[i]->GetYaxis()->SetRangeUser(0., 1.);
    gTs[i]->GetXaxis()->SetTitle("Angle [deg]");
    gTs[i]->GetYaxis()->SetTitle("Probability");
    gTs[i]->GetYaxis()->SetTitleOffset(1.2);
    gTs[i]->SetTitle(Form("Thickness %.2f mm", thicknesses[i]));
    gTs[i]->SetMarkerStyle(20);
    gTs[i]->SetMarkerColor(kBlue);
    gTs[i]->SetLineColor(kBlue);
    gTs[i]->SetLineWidth(2);
    gTs[i]->Draw("ALP");

    // Draw the reflected graph
    gRs[i]->SetMarkerStyle(21);
    gRs[i]->SetMarkerColor(kRed);
    gRs[i]->SetLineColor(kRed);
    gRs[i]->SetLineWidth(2);
    gRs[i]->Draw("LP SAME");

    gPad->Update();

    // Legend
    TLegend *legend = new TLegend(.70, .7, .9, .9, "Legend");
    legend->AddEntry(gTs[i], "Transmittance", "L P");
    legend->AddEntry(gRs[i], "Reflectance", "L P");
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
  cAngle->Close();
  fileAngleAnalysis->Close();
}

int main() {
  reflTransm();

  return EXIT_SUCCESS;
}