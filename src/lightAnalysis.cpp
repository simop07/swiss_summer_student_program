// To compile in SHELL:
// "g++ lightAnalysis.cpp `root-config --cflags --libs`"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
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
  gStyle->SetTitleXOffset(0.9f);
  gStyle->SetTitleYOffset(0.7f);
  gStyle->SetLineScalePS(1);
  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleYSize(0.05);
  // gStyle->SetPadTopMargin(0.02);
  // gStyle->SetPadBottomMargin(0.6);  // More room for x-axis labels
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

// C(x) function for PMT geometry
Double_t CFunc(Double_t *x, Double_t *par) {
  // par[0] = Angle (degrees)

  Double_t xVal = x[0];

  Double_t fitVal =
      TMath::Sqrt(1 - ((xVal - par[0]) / (25.)) * ((xVal - par[0]) / (25.)));

  return fitVal;
}

// T(x) = C(x) * F(x)
Double_t TFunc(Double_t *x, Double_t *par) {
  // Parameters:
  // par[0] = R_1
  // par[1] = Angle (degrees)
  // par[2] = sigma
  // par[3] = R_2
  // par[4] = Norm factor

  Double_t fVal = FFunc(x, par);      // Uses par[0]-par[4]
  Double_t cVal = CFunc(x, &par[1]);  // Pass only the center angle for C

  return fVal * cVal;
}

Double_t Integrand(Double_t *xx, Double_t *par) {
  // par[0] = R_1
  // par[1] = Angle (degrees)
  // par[2] = sigma
  // par[3] = R_2
  // par[4] = Norm factor

  Double_t xVal = xx[0];

  Double_t fVal = FFunc(xx, par);      // Uses par[0]-par[4]
  Double_t cVal = CFunc(xx, &par[1]);  // Pass only the center angle for C

  Double_t dx = xVal - par[1];
  Double_t dx2 = dx * dx;

  // Final integrand
  Double_t integrand = fVal * dx2 / (25. * 25. * 25. * cVal);
  return integrand;
}

Double_t Integrand2(Double_t *xx, Double_t *par) {
  // par[0] = R_1
  // par[1] = Angle (degrees)
  // par[2] = sigma
  // par[3] = R_2
  // par[4] = Norm factor

  Double_t xVal = xx[0];

  // F(x)
  Double_t fVal = FFunc(xx, par);

  // C(x)
  Double_t cVal = CFunc(xx, &par[1]);  // Only the center angle

  Double_t dx = xVal - par[1];

  // Derivative term dC/dtheta = (x - theta) / (L^2 * sqrt(1 - ((x-theta)/L)^2))
  Double_t derivC = dx / (25. * 25. * cVal);

  // Final integrand
  Double_t integrand = fVal * derivC;

  return integrand;
}

// dF/dR1 contribution
Double_t IntegrandR1(Double_t *xx, Double_t *par) {
  // par[0] = R1
  // par[1] = angle (deg)
  // par[2] = sigma
  // par[3] = R2
  // par[4] = Norm
  Double_t f = TMath::Exp(-TMath::Power(alphaFunc(xx[0], par[1]), 2) /
                          (2.0 * TMath::Power(par[2], 2)));
  Double_t c = CFunc(xx, &par[1]);
  return (f / par[4]) * c;  // C(x) * dF/dR1
}

// dF/dR2 contribution
Double_t IntegrandR2(Double_t *xx, Double_t *par) {
  Double_t c = CFunc(xx, &par[1]);
  Double_t f = TMath::Cos(xx[0] * TMath::Pi() / 180.0);
  return (f / par[4]) * c;  // C(x) * dF/dR2
}

// dF/dsigma contribution
Double_t IntegrandSigma(Double_t *xx, Double_t *par) {
  Double_t alpha = alphaFunc(xx[0], par[1]);
  Double_t expo =
      TMath::Exp(-TMath::Power(alpha, 2) / (2.0 * TMath::Power(par[2], 2)));
  Double_t term = par[0] * expo * (alpha * alpha) /
                  (TMath::Power(par[2], 3) * par[4]);  // dF/dsigma
  Double_t c = CFunc(xx, &par[1]);
  return term * c;  // C(x) * dF/dsigma
}

// Print information
void printResults(std::string const &angleLabel,
                  std::array<double, 8> const &thicknesses,
                  std::vector<double> const &probT,
                  std::vector<double> const &probTErr,
                  std::vector<double> const &probR,
                  std::vector<double> const &probRErr,
                  std::vector<double> const &probRErrStat,
                  std::vector<double> const &probRErrSys) {
  std::cout << "\n\n*** Results in trigger region for " << angleLabel
            << " ***\n";
  std::cout << std::setw(15) << std::left << "Thickness[mm]" << std::setw(15)
            << "Prob_T" << std::setw(15) << "ErrT" << std::setw(15) << "Prob_R"
            << std::setw(15) << "ErrR" << std::setw(15) << "ErrR(stat)"
            << std::setw(15) << "ErrR(sys)" << std::setw(15) << "FracR(stat)"
            << std::setw(15) << "FracR(sys)" << '\n';
  for (size_t i = 0; i < thicknesses.size(); ++i) {
    std::cout << std::fixed << std::setprecision(3) << std::setw(15)
              << thicknesses[i] << std::setw(15) << probT[i] << std::setw(15)
              << probTErr[i] << std::setw(15) << probR[i] << std::setw(15)
              << probRErr[i] << std::setw(15) << probRErrStat[i]
              << std::setw(15) << probRErrSys[i] << std::setw(15)
              << probRErrStat[i] / probR[i] << std::setw(15)
              << probRErrSys[i] / probR[i] << '\n';
  }
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
    /*   mg[i] = (TMultiGraph *)files[i]->Get("Regions of pulses");
      canvases[i]->cd();
      mg[i]->Draw("ALPE");
      mg[i]->SetTitle("Pulses");
      mg[i]->SetName("Regions of pulses");
      mg[i]->GetXaxis()->SetTitle("Time since \"trigger\" [#mus]");
      mg[i]->GetYaxis()->SetTitle("Voltage [mV]");

      // Write canvases on file
      fileLightAnalysis->cd();
      canvases[i]->Write();
      canvases[i]->Close(); */

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

  // Define errors for reflectance and transmittance
  std::vector<double> incTErrors{};
  std::vector<double> incRErrors{};
  std::vector<double> transmErrors{};
  std::vector<double> reflErrors{};
  incTErrors.reserve(nRegions);
  incRErrors.reserve(nRegions);
  transmErrors.reserve(nRegions);
  reflErrors.reserve(nRegions);

  // Print region properties
  int counter{};
  double rate{};
  double rateCorr{};
  double rateCorrRefl{};
  double rateError{};
  double rateCorrError{};
  double rateCorrReflError{};

  // Loop on files
  for (PhotonData const &pd : photondata) {
    std::cout << Form("\n *** %s", namesF[counter].c_str())
              << " photon data ***\n";

    // Loop on regions
    for (int i{}; i < nRegions; ++i) {
      switch (i) {
        case 0:  // Pre-trigger
          rate = pd.preTrigger.PECounter / (pd.preTrigger.deltaT * 60000.0);
          rateError = rate * std::sqrt((0.02 / 0.9) * (0.02 / 0.9) +
                                       (0.03 / 0.9) * (0.03 / 0.9));

          switch (counter) {
            case 0:
              incT.push_back(rate);
              incTErrors.push_back(rateError);
              break;
            case 1:
              incR.push_back(rate);
              incRErrors.push_back(rateError);
              break;
            case 2:
              transm.push_back(rate);
              transmErrors.push_back(rateError);
              break;
            case 3:
              rateCorrRefl = rate - incR[i];
              rateCorrReflError = std::sqrt(rateError * rateError +
                                            incRErrors[i] * incRErrors[i]);
              refl.push_back(rateCorrRefl);
              reflErrors.push_back(rateCorrReflError);
              break;
          }
          break;

        case 1:  // Trigger
          rateCorr = (pd.inTrigger.PECounter -
                      (rate * pd.inTrigger.deltaT * 60000.0)) /
                     (pd.inTrigger.deltaT * 60000.0);
          rateCorrError = std::sqrt(
              rateError * rateError +
              (std::sqrt((0.02 / 0.9) * (0.02 / 0.9) +
                         (0.03 / 0.9) * (0.03 / 0.9)) *
               pd.inTrigger.PECounter / (60000. * pd.inTrigger.deltaT)) *
                  (std::sqrt((0.02 / 0.9) * (0.02 / 0.9) +
                             (0.03 / 0.9) * (0.03 / 0.9)) *
                   pd.inTrigger.PECounter / (60000. * pd.inTrigger.deltaT)));
          switch (counter) {
            case 0:
              incT.push_back(rateCorr);
              incTErrors.push_back(rateCorrError);
              break;
            case 1:
              incR.push_back(rateCorr);
              incRErrors.push_back(rateCorrError);
              break;
            case 2:
              transm.push_back(rateCorr);
              transmErrors.push_back(rateCorrError);
              break;
            case 3:
              rateCorrRefl = rateCorr - incR[i];
              rateCorrReflError = std::sqrt(rateCorrError * rateCorrError +
                                            incRErrors[i] * incRErrors[i]);
              refl.push_back(rateCorrRefl);
              reflErrors.push_back(rateCorrReflError);
              break;
          }
          break;

        case 2:  // Post trigger 1
          rateCorr = (pd.postTrigger1.PECounter -
                      (rate * pd.postTrigger1.deltaT * 60000.0)) /
                     (pd.postTrigger1.deltaT * 60000.0);
          rateCorrError = std::sqrt(
              rateError * rateError +
              (std::sqrt((0.02 / 0.9) * (0.02 / 0.9) +
                         (0.03 / 0.9) * (0.03 / 0.9)) *
               pd.postTrigger1.PECounter / (60000. * pd.postTrigger1.deltaT)) *
                  (std::sqrt((0.02 / 0.9) * (0.02 / 0.9) +
                             (0.03 / 0.9) * (0.03 / 0.9)) *
                   pd.postTrigger1.PECounter /
                   (60000. * pd.postTrigger1.deltaT)));
          switch (counter) {
            case 0:
              incT.push_back(rateCorr);
              incTErrors.push_back(rateCorrError);
              break;
            case 1:
              incR.push_back(rateCorr);
              incRErrors.push_back(rateCorrError);
              break;
            case 2:
              transm.push_back(rateCorr);
              transmErrors.push_back(rateCorrError);
              break;
            case 3:
              rateCorrRefl = rateCorr - incR[i];
              rateCorrReflError = std::sqrt(rateCorrError * rateCorrError +
                                            incRErrors[i] * incRErrors[i]);
              refl.push_back(rateCorrRefl);
              reflErrors.push_back(rateCorrReflError);
              break;
          }
          break;

        case 3:  // Post trigger 2
          rateCorr = (pd.postTrigger2.PECounter -
                      (rate * pd.postTrigger2.deltaT * 60000.0)) /
                     (pd.postTrigger2.deltaT * 60000.0);
          rateCorrError = std::sqrt(
              rateError * rateError +
              (std::sqrt((0.02 / 0.9) * (0.02 / 0.9) +
                         (0.03 / 0.9) * (0.03 / 0.9)) *
               pd.postTrigger2.PECounter / (60000. * pd.postTrigger2.deltaT)) *
                  (std::sqrt((0.02 / 0.9) * (0.02 / 0.9) +
                             (0.03 / 0.9) * (0.03 / 0.9)) *
                   pd.postTrigger2.PECounter /
                   (60000. * pd.postTrigger2.deltaT)));
          switch (counter) {
            case 0:
              incT.push_back(rateCorr);
              incTErrors.push_back(rateCorrError);
              break;
            case 1:
              incR.push_back(rateCorr);
              incRErrors.push_back(rateCorrError);
              break;
            case 2:
              transm.push_back(rateCorr);
              transmErrors.push_back(rateCorrError);
              break;
            case 3:
              rateCorrRefl = rateCorr - incR[i];
              rateCorrReflError = std::sqrt(rateCorrError * rateCorrError +
                                            incRErrors[i] * incRErrors[i]);
              refl.push_back(rateCorrRefl);
              reflErrors.push_back(rateCorrReflError);
              break;
          }
          break;
      }
    }
    ++counter;
  }

  // Print photon information in each region
  /* std::string titles[11] = {"Region",        "Inc_T [PE/ns]", "err(Inc_T)",
                            "Inc_R [PE/ns]", "err(Inc_R)",    "Transm [PE/ns]",
                            "err(Transm)",   "Refl [PE/ns]",  "err(Refl)",
                            "Prob_T",        "Prob_R"};
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
              << std::setw(20) << incTErrors[i] << std::setw(20) << incR[i]
              << std::setw(20) << incRErrors[i] << std::setw(20) << transm[i]
              << std::setw(20) << transmErrors[i] << std::setw(20) << refl[i]
              << std::setw(20) << reflErrors[i] << std::setw(20) << probT
              << std::setw(20) << probR << '\n';
  }
  std::cout << "\n\n"; */

  // Create point (Prob_T, Prob_R) with errors
  Point p;
  p.x = transm[1] / incT[1];
  p.y = refl[1] / incT[1];
  p.xErr =
      p.x *
      std::sqrt((transmErrors[1] / transm[1]) * (transmErrors[1] / transm[1]) +
                (incTErrors[1] / incT[1]) * (incTErrors[1] / incT[1]));

  p.yErr =
      p.y * std::sqrt((reflErrors[1] / refl[1]) * (reflErrors[1] / refl[1]) +
                      (incTErrors[1] / incT[1]) * (incTErrors[1] / incT[1]));

  return p;
}

// Create plots using points (Prob_T, Prob_R) for different configurations
void reflTransm() {
  // Define angle vector
  std::vector<int> angles{30, 45, 60, 61, 62};

  // Define useful vector for thicknesses and colours
  std::array<double, 8> thicknesses{0.20, 0.80, 1.55, 2.05,
                                    3.10, 3.60, 4.10, 5.15};

  // Errors on thicknesses
  std::array<double, 8> thicknessesErr{0.05, 0.05, 0.05, 0.05,
                                       0.05, 0.05, 0.05, 0.05};

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
  std::vector<double> prob61T;
  std::vector<double> prob61R;
  std::vector<double> prob62T;
  std::vector<double> prob62R;
  std::vector<double> prob30TErr;
  std::vector<double> prob30RErr;
  std::vector<double> prob30RErrStat;
  std::vector<double> prob30RErrSys;
  std::vector<double> prob45TErr;
  std::vector<double> prob45RErr;
  std::vector<double> prob45RErrStat;
  std::vector<double> prob45RErrSys;
  std::vector<double> prob60TErr;
  std::vector<double> prob60RErr;
  std::vector<double> prob60RErrStat;
  std::vector<double> prob60RErrSys;
  std::vector<double> prob61TErr;
  std::vector<double> prob61RErr;
  std::vector<double> prob61RErrStat;
  std::vector<double> prob61RErrSys;
  std::vector<double> prob62TErr;
  std::vector<double> prob62RErr;
  std::vector<double> prob62RErrStat;
  std::vector<double> prob62RErrSys;

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
  double a61{};
  double lambda61{};
  double b61{};
  double a62{};
  double lambda62{};
  double b62{};
  double a30Error{};
  double lambda30Error{};
  double b30Error{};
  double a45Error{};
  double lambda45Error{};
  double b45Error{};
  double a60Error{};
  double lambda60Error{};
  double b60Error{};
  double a61Error{};
  double lambda61Error{};
  double b61Error{};
  double a62Error{};
  double lambda62Error{};
  double b62Error{};
  std::vector<double> geoFactors{};
  std::vector<double> geoFactorsNew{};
  std::vector<double> tPlusRfactors{};
  std::vector<double> tPlusRfactorsErr{};
  std::vector<double> pValuesBeerLambert{};
  std::vector<double> chiReducedBeerLambert{};
  std::vector<double> chiReducedConstantFit{};
  std::vector<double> pValuesConstantFit{};

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
    std::vector<double> probT, probTErr;
    std::vector<double> probR, probRErr;
    std::vector<double> probRErrStat, probRErrSys;
    for (auto &p : points) {
      probT.push_back(p.x);
      probTErr.push_back(p.xErr);
      probR.push_back(p.y);
      probRErr.push_back(p.yErr);
    }

    // Define angles for shaded region
    double angle1{};
    double angle2{};
    if (angle == 61 || angle == 62) {
      angle1 = 60. - 25.;
      angle2 = 60. + 25.;
    } else {
      angle1 = angle - 25.;
      angle2 = angle + 25.;
    }

    // Create canvas
    TCanvas *cFunc = new TCanvas("cFunc", "Angular Correction", 800, 600);

    // Create a frame for axes, titles, and range
    TH1F *frame =
        new TH1F("frame", Form("Correction function for %d^{#circ}", angle),
                 100, -90, 90);
    frame->SetXTitle("#theta_{refl} [deg]");
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

      case 61:
        fFFunction->FixParameter(0, 240);
        fFFunction->FixParameter(1, 60.);
        fFFunction->FixParameter(2, 0.1314);
        fFFunction->FixParameter(3, 66.51);
        fFFunction->FixParameter(4, 273.37408);
        break;

      case 62:
        fFFunction->FixParameter(0, 240);
        fFFunction->FixParameter(1, 60.);
        fFFunction->FixParameter(2, 0.1314);
        fFFunction->FixParameter(3, 66.51);
        fFFunction->FixParameter(4, 273.37408);
        break;
    }

    // Create T(x) = C(x)*F(x) with 5 params from F
    TF1 *fTFunction = new TF1(Form("fTFunction%d", angle), TFunc, -90., 90., 5);
    for (int i = 0; i < 5; i++) {
      fTFunction->FixParameter(i, fFFunction->GetParameter(i));
    }

    // Draw F function
    fFFunction->SetLineColor(kBlue + 2);
    fFFunction->SetLineWidth(2);

    // Draw T function
    fTFunction->SetLineColor(kRed + 2);
    fTFunction->SetLineWidth(2);

    // Create shaded region for F(x)
    TF1 *fFRange = (TF1 *)fFFunction->Clone();
    fFRange->SetRange(angle1, angle2);
    fFRange->SetFillColor(kBlue + 2);
    fFRange->SetFillStyle(1001);

    // Create shaded region for T(x)
    TF1 *fTRange = (TF1 *)fTFunction->Clone();
    fTRange->SetRange(angle1, angle2);
    fTRange->SetFillColor(kRed + 2);
    fTRange->SetFillStyle(1001);

    // Draw everything
    cFunc->cd();
    frame->Draw();
    frame->Draw("SAME AXIS");
    fFFunction->Draw("SAME");
    fTFunction->Draw("SAME");
    fFRange->Draw("SAME FC");
    fTRange->Draw("SAME FC");

    // Add legend
    TLegend *leg = new TLegend(.70, .7, .9, .9);
    leg->AddEntry(fFFunction, "F(#theta_{refl})", "L");
    leg->AddEntry(fTFunction, "C(#theta_{refl})#timesF(#theta_{refl})", "L");
    leg->AddEntry(fFRange, "Without C(#theta_{refl})", "F");
    leg->AddEntry(fTRange, "With C(#theta_{refl})", "F");
    leg->Draw();

    // Add text labels for angle1 and angle2
    TLatex *text1 = new TLatex(angle1, fTFunction->Eval(angle1) + 0.05,
                               Form("#color[1]{%.0f^{#circ}}", angle1));
    text1->SetTextAlign(22);  // Center alignment
    text1->Draw();

    TLatex *text2 = new TLatex(angle2, fTFunction->Eval(angle2) + 0.05,
                               Form("#color[1]{%.0f^{#circ}}", angle2));
    text2->SetTextAlign(22);
    text2->Draw();

    // Update canvas
    cFunc->Update();

    // Compute denominators
    double denomF = fFFunction->Integral(-90., 90.);

    // Compute numerators
    double numF = fFFunction->Integral(angle1, angle2);
    double numT = fTFunction->Integral(angle1, angle2);

    // Ratios
    double corrF = numF / denomF;  // Old factor
    double corrT = numT / denomF;  // New factor

    // Define errors for the parameters
    double deltaA{1.};
    double deltaTheta{1.};
    double deltaR1{};
    double deltaR2{};
    double deltaSigma{};
    if (angle == 30) {
      deltaR1 = 2.711;
      deltaR2 = 2.067;
      deltaSigma = 0.01768;
    } else if (angle == 45) {
      deltaR1 = 0.3002;
      deltaR2 = 1.483;
      deltaSigma = 0.006123;
    } else if (angle == 60) {
      deltaR1 = 0.2845;
      deltaR2 = 3.84;
      deltaSigma = 0.003067;
    } else if (angle == 61) {
      deltaR1 = 0.2845;
      deltaR2 = 3.84;
      deltaSigma = 0.003067;
    } else if (angle == 62) {
      deltaR1 = 0.2845;
      deltaR2 = 3.84;
      deltaSigma = 0.003067;
    }

    // Build integrand functions
    TF1 *fIntegrand =
        new TF1(Form("fIntegrand%d", angle), Integrand, -90., 90., 5);
    TF1 *fIntegrand2 =
        new TF1(Form("fIntegrand2%d", angle), Integrand2, -90., 90., 5);
    TF1 *fIntegrandR1 =
        new TF1(Form("fIntegrandR1%d", angle), IntegrandR1, -90., 90., 5);
    TF1 *fIntegrandR2 =
        new TF1(Form("fIntegrandR2%d", angle), IntegrandR2, -90., 90., 5);
    TF1 *fIntegrandSigma =
        new TF1(Form("fIntegrandSigma%d", angle), IntegrandSigma, -90., 90., 5);

    // Insert parameters
    for (int i = 0; i < 5; i++) {
      fIntegrand->FixParameter(i, fFFunction->GetParameter(i));
      fIntegrand2->FixParameter(i, fFFunction->GetParameter(i));
      fIntegrandR1->FixParameter(i, fFFunction->GetParameter(i));
      fIntegrandR2->FixParameter(i, fFFunction->GetParameter(i));
      fIntegrandSigma->FixParameter(i, fFFunction->GetParameter(i));
    }

    // Store / print results
    geoFactors.push_back(corrF);
    geoFactorsNew.push_back(corrT);
    /*   std::cout << std::fixed << std::setprecision(6)
                << "\nCorrection factor F for " << angle << "° (" << angle1
                << "°-" << angle2 << "°) = " << corrF << "\n";
      std::cout << "Correction factor T for " << angle << "° (" << angle1 <<
      "°-"
                << angle2 << "°) = " << corrT << "\n\n"; */

    // Apply correction factor to reflectance
    for (size_t i = 0; i < probR.size(); ++i) {
      double R = probR[i];
      double sigmaR = probRErr[i];
      double sigmaCorrT =
          (1.0 / denomF) *
          (std::abs(fIntegrand->Integral(angle1, angle2)) * deltaA  // δa
           + std::abs(fIntegrand2->Integral(angle1, angle2, 1e-11)) *
                 deltaTheta                                              // δθ
           + std::abs(fIntegrandR1->Integral(angle1, angle2)) * deltaR1  // δR1
           + std::abs(fIntegrandR2->Integral(angle1, angle2)) * deltaR2  // δR2
           + std::abs(fIntegrandSigma->Integral(angle1, angle2)) *
                 deltaSigma  // δσ
          );
      probR[i] = R / corrT;

      // Propagate stat and sys errors separately, then total
      double stat = sigmaR / corrT;
      double sys = (R * sigmaCorrT) / (corrT * corrT);

      probRErr[i] = sqrt(stat * stat + sys * sys);
      probRErrStat.push_back(stat);
      probRErrSys.push_back(sys);
    }

    // Save probabilities into vectors
    switch (angle) {
      case 30:
        std::copy(probT.begin(), probT.end(), std::back_inserter(prob30T));
        std::copy(probR.begin(), probR.end(), std::back_inserter(prob30R));
        std::copy(probTErr.begin(), probTErr.end(),
                  std::back_inserter(prob30TErr));
        std::copy(probRErr.begin(), probRErr.end(),
                  std::back_inserter(prob30RErr));
        std::copy(probRErrStat.begin(), probRErrStat.end(),
                  std::back_inserter(prob30RErrStat));
        std::copy(probRErrSys.begin(), probRErrSys.end(),
                  std::back_inserter(prob30RErrSys));
        break;

      case 45:
        std::copy(probT.begin(), probT.end(), std::back_inserter(prob45T));
        std::copy(probR.begin(), probR.end(), std::back_inserter(prob45R));
        std::copy(probTErr.begin(), probTErr.end(),
                  std::back_inserter(prob45TErr));
        std::copy(probRErr.begin(), probRErr.end(),
                  std::back_inserter(prob45RErr));
        std::copy(probRErrStat.begin(), probRErrStat.end(),
                  std::back_inserter(prob45RErrStat));
        std::copy(probRErrSys.begin(), probRErrSys.end(),
                  std::back_inserter(prob45RErrSys));
        break;

      case 60:
        std::copy(probT.begin(), probT.end(), std::back_inserter(prob60T));
        std::copy(probR.begin(), probR.end(), std::back_inserter(prob60R));
        std::copy(probTErr.begin(), probTErr.end(),
                  std::back_inserter(prob60TErr));
        std::copy(probRErr.begin(), probRErr.end(),
                  std::back_inserter(prob60RErr));
        std::copy(probRErrStat.begin(), probRErrStat.end(),
                  std::back_inserter(prob60RErrStat));
        std::copy(probRErrSys.begin(), probRErrSys.end(),
                  std::back_inserter(prob60RErrSys));
        break;

      case 61:
        std::copy(probT.begin(), probT.end(), std::back_inserter(prob61T));
        std::copy(probR.begin(), probR.end(), std::back_inserter(prob61R));
        std::copy(probTErr.begin(), probTErr.end(),
                  std::back_inserter(prob61TErr));
        std::copy(probRErr.begin(), probRErr.end(),
                  std::back_inserter(prob61RErr));
        std::copy(probRErrStat.begin(), probRErrStat.end(),
                  std::back_inserter(prob61RErrStat));
        std::copy(probRErrSys.begin(), probRErrSys.end(),
                  std::back_inserter(prob61RErrSys));
        break;

      case 62:
        std::copy(probT.begin(), probT.end(), std::back_inserter(prob62T));
        std::copy(probR.begin(), probR.end(), std::back_inserter(prob62R));
        std::copy(probTErr.begin(), probTErr.end(),
                  std::back_inserter(prob62TErr));
        std::copy(probRErr.begin(), probRErr.end(),
                  std::back_inserter(prob62RErr));
        std::copy(probRErrStat.begin(), probRErrStat.end(),
                  std::back_inserter(prob62RErrStat));
        std::copy(probRErrSys.begin(), probRErrSys.end(),
                  std::back_inserter(prob62RErrSys));
        break;
    }

    // Print table with (Prob_T, Prob_R) and thicknesses
    /* std::cout << "\n\n*** Results in trigger region ***\n\n";
    std::cout << std::setw(20) << std::left << "Thickness [mm]" << std::setw(20)
              << "Prob_T" << std::setw(20) << "Prob_R" << '\n';
    for (size_t i = 0; i < thicknesses.size(); ++i) {
      std::cout << std::fixed << std::setprecision(3);
      std::cout << std::setw(20) << thicknesses[i] << std::setw(20) << probT[i]
                << " ± " << probTErr[i] << std::setw(20) << probR[i] << " ± "
                << probRErr[i] << '\n';
    }
 */
    // Multigraphs
    TMultiGraph *mg = new TMultiGraph();
    TMultiGraph *mgReflVsTransm = new TMultiGraph();
    std::array<TGraphErrors *, 2> g{
        new TGraphErrors(thicknesses.size(), thicknesses.data(), probT.data(),
                         thicknessesErr.data(), probTErr.data()),
        new TGraphErrors(thicknesses.size(), thicknesses.data(), probR.data(),
                         thicknessesErr.data(), probRErr.data())};
    std::array<std::string, 2> names{"Transmittance", "Reflectance"};
    for (size_t i = 0; i < names.size(); ++i) {
      g[i]->SetMarkerStyle(20);
      g[i]->SetLineWidth(2);
      g[i]->SetTitle(names[i].c_str());
      mg->Add(g[i], "LP");
    }

    // Insert colours
    g[0]->SetMarkerColor(kBlue);
    g[0]->SetLineColor(kBlue);
    g[1]->SetMarkerColor(kRed);
    g[1]->SetLineColor(kRed);

    // Draw canvas
    TCanvas *cReflVsTransm = new TCanvas(Form("cReflVsTransm%d", angle),
                                         "Refl vs transm", 1500, 700);

    // Fit transmittance
    TF1 *beerLambertTransm =
        new TF1(Form("beerLambert%dTransm", angle), beerLambert, 0.1, 5.3, 3);
    beerLambertTransm->SetLineColor(kRed);
    beerLambertTransm->SetLineWidth(4);
    beerLambertTransm->SetLineStyle(2);
    beerLambertTransm->SetTitle("Transmittance fit");
    beerLambertTransm->SetParNames("A", "#lambda_t", "B");
    beerLambertTransm->SetParameter(0, 1.);
    beerLambertTransm->SetParameter(1, 0.7);
    beerLambertTransm->SetParameter(2, 0.0);
    g[0]->Fit(beerLambertTransm, "R");

    // Save transmittance fit parameters
    switch (angle) {
      case 30:
        a30 = beerLambertTransm->GetParameter(0);
        lambda30 = beerLambertTransm->GetParameter(1);
        b30 = beerLambertTransm->GetParameter(2);
        a30Error = beerLambertTransm->GetParError(0);
        lambda30Error = beerLambertTransm->GetParError(1);
        b30Error = beerLambertTransm->GetParError(2);
        break;

      case 45:
        a45 = beerLambertTransm->GetParameter(0);
        lambda45 = beerLambertTransm->GetParameter(1);
        b45 = beerLambertTransm->GetParameter(2);
        a45Error = beerLambertTransm->GetParError(0);
        lambda45Error = beerLambertTransm->GetParError(1);
        b45Error = beerLambertTransm->GetParError(2);
        break;

      case 60:
        a60 = beerLambertTransm->GetParameter(0);
        lambda60 = beerLambertTransm->GetParameter(1);
        b60 = beerLambertTransm->GetParameter(2);
        a60Error = beerLambertTransm->GetParError(0);
        lambda60Error = beerLambertTransm->GetParError(1);
        b60Error = beerLambertTransm->GetParError(2);
        break;

      case 61:
        a61 = beerLambertTransm->GetParameter(0);
        lambda61 = beerLambertTransm->GetParameter(1);
        b61 = beerLambertTransm->GetParameter(2);
        a61Error = beerLambertTransm->GetParError(0);
        lambda61Error = beerLambertTransm->GetParError(1);
        b61Error = beerLambertTransm->GetParError(2);
        break;

      case 62:
        a62 = beerLambertTransm->GetParameter(0);
        lambda62 = beerLambertTransm->GetParameter(1);
        b62 = beerLambertTransm->GetParameter(2);
        a62Error = beerLambertTransm->GetParError(0);
        lambda62Error = beerLambertTransm->GetParError(1);
        b62Error = beerLambertTransm->GetParError(2);
        break;
    }

    // Compute chi2 and p-value
    double chi2BL = beerLambertTransm->GetChisquare();
    int ndfBL = beerLambertTransm->GetNDF();
    chiReducedBeerLambert.push_back(chi2BL / ndfBL);
    pValuesBeerLambert.push_back(beerLambertTransm->GetProb());

    // Draw canvas
    TCanvas *cDeg = new TCanvas(Form("c%dDeg", angle),
                                Form("%d degrees", angle), 1500, 900);

    // Top pad
    TPad *pad1 = new TPad("pad1", "Main", 0.0, 0.35, 1.0, 1.0);  // Make top 65%
    pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();

    // Draw main multigraph here
    mg->Draw("ALPE");
    mg->SetTitle(Form("%d^{#circ} probability", angle));
    mg->GetXaxis()->SetTitle("Thickness [mm]");
    mg->GetXaxis()->SetLabelSize(0);  // Hide labels on top graph
    mg->GetXaxis()->SetTitleSize(0);  // Hide title on top graph
    mg->GetYaxis()->SetTitle("Probability");
    mg->GetYaxis()->SetTitleOffset(0.8);
    mg->GetYaxis()->SetLabelSize(0.035);  // Smaller y labels
    mg->GetYaxis()->SetTitleSize(0.06);   // Smaller y title

    TLegend *legend2 = new TLegend(.70, .7, .9, .9);
    legend2->SetTextSize(0.035);  // Make legend text smaller
    legend2->AddEntry(g[0], "Transmittance", "LEP");
    legend2->AddEntry(g[1], "Reflectance", "LEP");
    legend2->AddEntry(beerLambertTransm, "Exponential fit", "L");
    legend2->Draw();

    cDeg->cd();

    // Bottom pad
    TPad *pad2 = new TPad("pad2", "Sum", 0.0, 0.0, 1.0, 0.35);
    // bottom 35%
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.3);  // More room for x-axis labels
    pad2->Draw();
    pad2->cd();

    // Build the sum graph
    std::vector<double> probSum(thicknesses.size());
    std::vector<double> probSumErr(thicknesses.size());
    for (size_t i = 0; i < thicknesses.size(); ++i) {
      probSum[i] = probT[i] + probR[i];
      probSumErr[i] =
          std::sqrt(probTErr[i] * probTErr[i] + probRErr[i] * probRErr[i]);
    }
    TGraphErrors *gSum =
        new TGraphErrors(thicknesses.size(), thicknesses.data(), probSum.data(),
                         thicknessesErr.data(), probSumErr.data());
    gSum->SetLineColor(kGreen + 2);
    gSum->SetMarkerColor(kGreen + 2);
    gSum->SetMarkerStyle(21);
    gSum->SetLineWidth(2);

    double xmin = 0.0, xmax = 5.3;
    mg->GetXaxis()->SetLimits(xmin, xmax);    // For multigraph
    gSum->GetXaxis()->SetLimits(xmin, xmax);  // For bottom graph
    gSum->Draw("ALPE");
    gSum->SetTitle("");
    gSum->GetXaxis()->SetTitle("Thickness [mm]");
    gSum->GetXaxis()->SetLabelSize(0.07);  // Slightly bigger (bottom axis only)
    gSum->GetXaxis()->SetTitleSize(0.1);
    gSum->GetYaxis()->SetTitle("T+R");
    gSum->GetYaxis()->SetTitleOffset(0.3);
    gSum->GetYaxis()->SetLabelSize(0.07);
    gSum->GetYaxis()->SetTitleSize(0.1);

    // Define constant fit function (just p0)
    TF1 *fConst = new TF1("fConst", "[0]", 0.1, 5.2);
    fConst->SetLineColor(kRed);
    fConst->SetLineStyle(2);
    fConst->SetLineWidth(4);
    fConst->SetParameter(0, 1.0);  // initial guess

    // Draw and fit
    gSum->Draw("ALPE");
    gSum->Fit(fConst, "R");  // "R" restricts to function range

    // Add legend for bottom pad
    TLegend *legendBottom = new TLegend(.70, .7, .9, .9);
    legendBottom->SetTextSize(0.06);
    legendBottom->AddEntry(gSum, "T + R", "LEP");
    legendBottom->AddEntry(fConst, "Constant fit", "L");
    legendBottom->Draw();

    // Optionally print the result
    double p0 = fConst->GetParameter(0);
    double p0err = fConst->GetParError(0);
    // std::cout << "Fit result: p0 = " << p0 << " ± " << p0err << std::endl;
    tPlusRfactors.push_back(p0);
    tPlusRfactorsErr.push_back(p0err);

    // Compute chi2 and p-value for constant fit
    double chi2Const = fConst->GetChisquare();
    int ndfConst = fConst->GetNDF();
    chiReducedConstantFit.push_back(chi2Const / ndfConst);
    pValuesConstantFit.push_back(fConst->GetProb());

    setFitStyle();

    // Transmittance vs Reflectance
    for (size_t i = 0; i < thicknesses.size(); ++i) {
      auto *gr = new TGraphErrors();
      gr->SetPoint(0, probR[i], probT[i]);
      gr->SetPointError(0, probRErr[i], probTErr[i]);
      gr->SetMarkerStyle(20);
      gr->SetMarkerSize(2);
      gr->SetMarkerColor(colours[i]);
      gr->SetLineColor(colours[i]);
      gr->SetTitle(Form("%.2f mm", thicknesses[i]));
      mgReflVsTransm->Add(gr, "P");
    }

    cReflVsTransm->cd();
    mgReflVsTransm->Draw("ALPE");
    mgReflVsTransm->SetTitle("Transmittance vs Reflectance");
    mgReflVsTransm->SetName("mgReflVsTransm");
    mgReflVsTransm->GetXaxis()->SetTitle("Reflectance prob");
    mgReflVsTransm->GetYaxis()->SetTitle("Transm prob");
    mgReflVsTransm->GetYaxis()->SetTitleOffset(1.2);

    // Build legend
    TLegend *leg2 = new TLegend(0.70, 0.70, 0.90, 0.90, "Thickness");
    TIter next(mgReflVsTransm->GetListOfGraphs());
    TGraphErrors *graph = nullptr;
    while ((graph = (TGraphErrors *)next())) {
      leg2->AddEntry(graph, graph->GetTitle(), "EP");
    }
    leg2->Draw();

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
  std::vector<double> angleVals(angles.begin(), angles.begin() + 3);

  // Draw a canvas for probabilities as function of angle fixing thicknesses
  TCanvas *cAngle = new TCanvas("cAngle", "Probabilities vs angle", 1500, 700);
  cAngle->Divide(4, 2);

  // Create vectors for graphs
  std::vector<TGraphErrors *> gTs;
  std::vector<TGraphErrors *> gRs;

  // Loop over thicknesses
  for (size_t i = 0; i < thicknesses.size(); ++i) {
    // Extract probabilities for this thickness across different angles
    double newProbT[3] = {prob30T[i], prob45T[i], prob61T[i]};
    double newProbR[3] = {prob30R[i], prob45R[i], prob61R[i]};
    double newProbTErr[3] = {prob30TErr[i], prob45TErr[i], prob61TErr[i]};
    double newProbRErr[3] = {prob30RErr[i], prob45RErr[i], prob61RErr[i]};

    // Push graphs into vectors
    gTs.push_back(
        new TGraphErrors(3, angleVals.data(), newProbT, nullptr, newProbTErr));
    gRs.push_back(
        new TGraphErrors(3, angleVals.data(), newProbR, nullptr, newProbRErr));
  }

  // Draw graphs in vectors
  for (size_t i = 0; i < gTs.size(); ++i) {
    cAngle->cd(i + 1);
    gTs[i]->GetYaxis()->SetRangeUser(0., 1.6);
    gTs[i]->GetXaxis()->SetTitle("Angle [deg]");
    // gTs[i]->GetYaxis()->SetTitle("Probability");
    gTs[i]->GetYaxis()->SetTitleOffset(1);
    gTs[i]->SetTitle(Form("%.2f mm", thicknesses[i]));
    gTs[i]->SetMarkerStyle(20);
    gTs[i]->SetMarkerColor(kBlue);
    gTs[i]->SetLineColor(kBlue);
    gTs[i]->SetLineWidth(2);
    gTs[i]->Draw("ALPE");

    // Draw the reflected graph
    gRs[i]->SetMarkerStyle(21);
    gRs[i]->SetMarkerColor(kRed);
    gRs[i]->SetLineColor(kRed);
    gRs[i]->SetLineWidth(2);
    gRs[i]->Draw("LPE SAME");

    gPad->Update();

    // Legend
    TLegend *legend = new TLegend(.70, .7, .9, .9);
    legend->AddEntry(gTs[i], "Transmittance", "L E P");
    legend->AddEntry(gRs[i], "Reflectance", "L E P");
    legend->Draw();
  }

  // Print info for each angle
  printResults("30°", thicknesses, prob30T, prob30TErr, prob30R, prob30RErr,
               prob30RErrStat, prob30RErrSys);
  printResults("45°", thicknesses, prob45T, prob45TErr, prob45R, prob45RErr,
               prob45RErrStat, prob45RErrSys);
  printResults("60°", thicknesses, prob60T, prob60TErr, prob60R, prob60RErr,
               prob60RErrStat, prob60RErrSys);
  printResults("61°", thicknesses, prob61T, prob61TErr, prob61R, prob61RErr,
               prob61RErrStat, prob61RErrSys);
  printResults("62°", thicknesses, prob62T, prob62TErr, prob62R, prob62RErr,
               prob62RErrStat, prob62RErrSys);

  // Transmittance fit results
  std::cout << "\n\n**Fit results for 30°**\n";
  std::cout << "- A           = " << a30 << " ± " << a30Error << '\n';
  std::cout << "- lambda_t    = " << lambda30 << " ± " << lambda30Error
            << " mm\n";
  std::cout << "- B           = " << b30 << " ± " << b30Error << '\n';

  std::cout << "\n**Fit results for 45°**\n";
  std::cout << "- A           = " << a45 << " ± " << a45Error << '\n';
  std::cout << "- lambda_t    = " << lambda45 << " ± " << lambda45Error
            << " mm\n";
  std::cout << "- B           = " << b45 << " ± " << b45Error << '\n';

  std::cout << "\n**Fit results for 60°**\n";
  std::cout << "- A           = " << a60 << " ± " << a60Error << '\n';
  std::cout << "- lambda_t    = " << lambda60 << " ± " << lambda60Error
            << " mm\n";
  std::cout << "- B           = " << b60 << " ± " << b60Error << '\n';

  std::cout << "\n**Fit results for 61°**\n";
  std::cout << "- A           = " << a61 << " ± " << a61Error << '\n';
  std::cout << "- lambda_t    = " << lambda61 << " ± " << lambda61Error
            << " mm\n";
  std::cout << "- B           = " << b61 << " ± " << b61Error << '\n';

  std::cout << "\n**Fit results for 62°**\n";
  std::cout << "- A           = " << a62 << " ± " << a60Error << '\n';
  std::cout << "- lambda_t    = " << lambda62 << " ± " << lambda60Error
            << " mm\n";
  std::cout << "- B           = " << b62 << " ± " << b60Error << '\n';

  // Print chi squared and p values
  std::cout << std::fixed << std::setprecision(20);
  std::cout << "\n=== Beer-Lambert Fit ===\n";
  for (size_t i = 0; i < chiReducedBeerLambert.size(); ++i) {
    std::cout << "Angle " << angles[i] << "° -> "
              << "Chi2/NDF = " << chiReducedBeerLambert[i]
              << ", p-value = " << pValuesBeerLambert[i] << "\n";
  }

  std::cout << "\n=== Constant Fit ===\n";
  for (size_t i = 0; i < chiReducedConstantFit.size(); ++i) {
    std::cout << "Angle " << angles[i] << "° -> "
              << "Chi2/NDF = " << chiReducedConstantFit[i]
              << ", p-value = " << pValuesConstantFit[i] << "\n";
  }

  // Print correction factors
  std::cout << "\n- Correction factor for 30° = " << geoFactors[0] << '\n';
  std::cout << "- Correction factor for 45° = " << geoFactors[1] << '\n';
  std::cout << "- Correction factor for 60° = " << geoFactors[2] << '\n';
  std::cout << "- Correction factor for 61° = " << geoFactors[3] << '\n';
  std::cout << "- Correction factor for 62° = " << geoFactors[4] << '\n';

  // Print new correction factors
  std::cout << "\n- Correction factor with geo correction for 30° = "
            << geoFactorsNew[0] << '\n';
  std::cout << "- Correction factor with geo correction for 45° = "
            << geoFactorsNew[1] << '\n';
  std::cout << "- Correction factor with geo correction for 60° = "
            << geoFactorsNew[2] << '\n';
  std::cout << "- Correction factor with geo correction for 61° = "
            << geoFactorsNew[3] << '\n';
  std::cout << "- Correction factor with geo correction for 62° = "
            << geoFactorsNew[4] << '\n';

  // Print transmittance + reflectance factors
  std::cout << "\n- Fit result 30°: p0 = " << tPlusRfactors[0] << " ± "
            << tPlusRfactorsErr[0] << std::endl;
  std::cout << "- Fit result 45°: p0 = " << tPlusRfactors[1] << " ± "
            << tPlusRfactorsErr[1] << std::endl;
  std::cout << "- Fit result 60°: p0 = " << tPlusRfactors[2] << " ± "
            << tPlusRfactorsErr[2] << std::endl;
  std::cout << "- Fit result 61°: p0 = " << tPlusRfactors[3] << " ± "
            << tPlusRfactorsErr[3] << std::endl;
  std::cout << "- Fit result 62°: p0 = " << tPlusRfactors[4] << " ± "
            << tPlusRfactorsErr[4] << std::endl;

  // Creating ROOT File
  TFile *fileAngleAnalysis =
      new TFile("./rootFiles/fileAngleAnalysis.root", "RECREATE");

  // Save canvas
  fileAngleAnalysis->cd();
  cAngle->Write();
  cAngle->Close();
  fileAngleAnalysis->Close();
}

void reflTransmNew() {
  std::string basePath = "./rootFiles/45DegreesNew";

  // Metal and scraped configs
  std::vector<std::string> metalConfigs{"metal_1", "metal_2", "metal_3",
                                        "metal_4"};
  std::vector<std::string> scrapedConfigs{"scraped_", "scraped__1",
                                          "scraped__2", "scraped__3"};

  // Angle for correction
  int angle = 45;
  double angle1 = angle - 25.;
  double angle2 = angle + 25.;

  // F and T functions for 45°
  TF1 *fFFunction = new TF1("fFFunction45", FFunc, -90., 90., 5);
  fFFunction->FixParameter(0, 90.0);
  fFFunction->FixParameter(1, 45.);
  fFFunction->FixParameter(2, -0.1469);
  fFFunction->FixParameter(3, 50.18);
  fFFunction->FixParameter(4, 125.6324);

  TF1 *fTFunction = new TF1("fTFunction45", TFunc, -90., 90., 5);
  for (int i = 0; i < 5; i++)
    fTFunction->FixParameter(i, fFFunction->GetParameter(i));

  double denomF = fFFunction->Integral(-90., 90.);
  double numT = fTFunction->Integral(angle1, angle2);
  double corrFactor45 = numT / denomF;

  // Integrand functions for reflectance correction
  double deltaA = 1., deltaTheta = 1., deltaR1 = 0.3002, deltaR2 = 1.483,
         deltaSigma = 0.006123;
  TF1 *fIntegrand = new TF1("fIntegrand45", Integrand, -90., 90., 5);
  TF1 *fIntegrand2 = new TF1("fIntegrand2_45", Integrand2, -90., 90., 5);
  TF1 *fIntegrandR1 = new TF1("fIntegrandR1_45", IntegrandR1, -90., 90., 5);
  TF1 *fIntegrandR2 = new TF1("fIntegrandR2_45", IntegrandR2, -90., 90., 5);
  TF1 *fIntegrandSigma =
      new TF1("fIntegrandSigma_45", IntegrandSigma, -90., 90., 5);
  for (int i = 0; i < 5; i++) {
    fIntegrand->FixParameter(i, fFFunction->GetParameter(i));
    fIntegrand2->FixParameter(i, fFFunction->GetParameter(i));
    fIntegrandR1->FixParameter(i, fFFunction->GetParameter(i));
    fIntegrandR2->FixParameter(i, fFFunction->GetParameter(i));
    fIntegrandSigma->FixParameter(i, fFFunction->GetParameter(i));
  }

  // Ccollect probT and probR (with correction)
  auto collectProbs = [&](std::vector<std::string> configs,
                          std::vector<double> &probT,
                          std::vector<double> &probR,
                          std::vector<double> &probTerr,
                          std::vector<double> &probRerr) {
    for (auto &conf : configs) {
      std::string folderPath = basePath + "/" + conf + "/w";
      Point p = lightAnalysis(folderPath, "");

      probT.push_back(p.x);
      probTerr.push_back(p.xErr);

      double R = p.y;
      double sigmaR = p.yErr;

      // Apply 45° correction as before
      double sigmaCorrT =
          (1.0 / denomF) *
          (std::abs(fIntegrand->Integral(angle1, angle2)) * deltaA +
           std::abs(fIntegrand2->Integral(angle1, angle2, 1e-11)) * deltaTheta +
           std::abs(fIntegrandR1->Integral(angle1, angle2)) * deltaR1 +
           std::abs(fIntegrandR2->Integral(angle1, angle2)) * deltaR2 +
           std::abs(fIntegrandSigma->Integral(angle1, angle2)) * deltaSigma);

      probR.push_back(R / corrFactor45);
      probRerr.push_back(
          sqrt(sigmaR * sigmaR / (corrFactor45 * corrFactor45) +
               sigmaCorrT * sigmaCorrT / (corrFactor45 * corrFactor45)));

      // Print table header
      std::cout << std::setw(12) << "Config" << std::setw(15) << "Prob_T"
                << std::setw(15) << "Err_T" << std::setw(15) << "Prob_R"
                << std::setw(15) << "Err_R" << std::endl;
      std::cout << std::string(72, '-') << std::endl;

      // Print values
      for (size_t i = 0; i < configs.size(); ++i) {
        std::cout << std::setw(12) << configs[i] << std::setw(15) << std::fixed
                  << std::setprecision(6) << probT[i] << std::setw(15)
                  << probTerr[i] << std::setw(15) << probR[i] << std::setw(15)
                  << probRerr[i] << std::endl;
      }
    }
  };

  // --- Collect for metal ---
  std::vector<double> probT_metal, probR_metal, probTerr_metal, probRerr_metal;
  collectProbs(metalConfigs, probT_metal, probR_metal, probTerr_metal,
               probRerr_metal);

  // --- Collect for scraped ---
  std::vector<double> probT_scraped, probR_scraped, probTerr_scraped,
      probRerr_scraped;
  collectProbs(scrapedConfigs, probT_scraped, probR_scraped, probTerr_scraped,
               probRerr_scraped);

  // --- Make multigraphs ---
  TFile *outFile = new TFile("./rootFiles/reflTransm_45Deg.root", "RECREATE");

  auto makeGraph = [&](std::vector<double> &probT, std::vector<double> &probR,
                       std::string name) {
    TMultiGraph *mg = new TMultiGraph();
    int n = probT.size();
    std::vector<double> x(n), ex(n, 0.);
    for (int i = 0; i < n; i++) x[i] = i + 1;
    TGraphErrors *gT = new TGraphErrors(n, x.data(), probT.data(), ex.data(),
                                        probTerr_metal.data());
    TGraphErrors *gR = new TGraphErrors(n, x.data(), probR.data(), ex.data(),
                                        probRerr_metal.data());
    gT->SetMarkerStyle(20);
    gT->SetMarkerColor(kBlue);
    gT->SetLineColor(kBlue);
    gT->SetTitle("Transmittance");
    gR->SetMarkerStyle(21);
    gR->SetMarkerColor(kRed);
    gR->SetLineColor(kRed);
    gR->SetTitle("Reflectance");
    mg->Add(gT, "LP");
    mg->Add(gR, "LP");
    mg->SetTitle(name.c_str());
    mg->Draw("A");
    mg->Write(name.c_str());
  };

  makeGraph(probT_metal, probR_metal, "MetalConfigs");
  makeGraph(probT_scraped, probR_scraped, "ScrapedConfigs");

  outFile->Close();
}

int main() {
  reflTransm();

  return EXIT_SUCCESS;
}