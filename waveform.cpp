// To compile in SHELL:
// "g++ waveformAnalysisPos.cpp waveform.cpp `root-config --cflags --libs`" Best
// data to show fit is DataF_CH0@DT5730S_59483_run_new_1300_2-3.5 with Refl PMT
// conversion factor inside data/miscellaneous

#include <time.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <map>
#include <numeric>
#include <sstream>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMultiGraph.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"
#include "waveformAnalysisPos.hpp"

// Define global constants
constexpr int nMinAnalysedRows{1};      // Minimum rows EXCLUDED
constexpr int nMaxAnalysedRows{30000};  // Maximum rows INCLUDED

// Asymmetric gaussian functions

Double_t asymGaussians(Double_t *x, Double_t *par) {
  // par[0] = N
  // par[1] = #mu
  // par[2] = #sigma_{1}
  // par[3] = #sigma_{2}

  Double_t xVal = x[0];
  Double_t fitVal =
      par[0] *
      (((xVal < par[1]) * (TMath::Exp(-0.5 * ((xVal - par[1]) / par[2]) *
                                      ((xVal - par[1]) / par[2])))) +
       ((xVal >= par[1]) * (TMath::Exp(-0.5 * ((xVal - par[1]) / par[3]) *
                                       ((xVal - par[1]) / par[3])))));
  return fitVal;
}

Double_t asym2Gaussians(Double_t *x, Double_t *par) {
  // par[0] = N^{1}
  // par[1] = #mu^{1}
  // par[2] = #sigma^{1}_{1}
  // par[3] = #sigma^{1}_{2}
  // par[4] = N^{2}
  // par[5] = #mu^{2}
  // par[6] = #sigma^{2}_{1}
  // par[7] = #sigma^{2}_{2}
  // par[8] = Background

  return asymGaussians(x, &par[0]) + asymGaussians(x, &par[4]) + par[8];
}

Double_t asym2GaussiansConstrained(Double_t *x, Double_t *par) {
  // par[0] = N^{1}
  // par[1] = #mu^{1}
  // par[2] = #sigma^{1}_{1}
  // par[3] = #sigma^{1}_{2}
  // par[4] = N^{2}
  // par[5] = #sigma^{2}_{1}
  // par[6] = #sigma^{2}_{2}
  // par[7] = Background

  return asymGaussians(x, new Double_t[4]{par[0], par[1], par[2], par[3]}) +
         asymGaussians(x, new Double_t[4]{par[4], 2 * par[1], par[5], par[6]}) +
         par[7];
}

Double_t asym2GaussiansConstrainedSameSigma(Double_t *x, Double_t *par) {
  // par[0] = N^{1}
  // par[1] = #mu^{1}
  // par[2] = #sigma^{1}_{1}
  // par[3] = #sigma^{1}_{2}
  // par[4] = N^{2}
  // par[5] = Background

  return asymGaussians(x, new Double_t[4]{par[0], par[1], par[2], par[3]}) +
         asymGaussians(x, new Double_t[4]{par[4], 2 * par[1], par[2], par[3]}) +
         par[5];
}

// Group functions for fitting PE histo
void fitPEHistoNoExp(TH1F *hPhotoElectrons) {
  // Import user defined function asymmetric gaussians for 1 PE
  TF1 *fAsymmetric1PE = new TF1("fAsymmetric1PE", asymGaussians, 0.6, 1.9, 4);
  fAsymmetric1PE->SetLineColor(kRed);
  fAsymmetric1PE->SetLineWidth(4);
  fAsymmetric1PE->SetLineStyle(2);
  fAsymmetric1PE->SetParNames("N^{1}", "#mu^{1}", "#sigma^{1}_{1}",
                              "#sigma^{1}_{2}");
  fAsymmetric1PE->SetParameter(0, 0.3);  // N^{1}
  fAsymmetric1PE->SetParameter(1, 1.);   // #mu^{1}
  fAsymmetric1PE->SetParameter(2, 0.5);  // #sigma^{1}_{1}
  fAsymmetric1PE->SetParameter(3, 0.5);  // #sigma^{1}_{2}

  // Import user defined function asymmetric gaussians for 2 PE
  TF1 *fAsymmetric2PE = new TF1("fAsymmetric2PE", asymGaussians, 1.9, 2.5, 4);
  fAsymmetric2PE->SetLineColor(kOrange + 2);
  fAsymmetric2PE->SetLineWidth(4);
  fAsymmetric2PE->SetLineStyle(2);
  fAsymmetric2PE->SetParNames("N^{2}", "#mu^{2}", "#sigma^{2}_{1}",
                              "#sigma^{2}_{2}");
  fAsymmetric2PE->SetParameter(0, 0.1);  // N^{2}
  fAsymmetric2PE->SetParameter(1, 2.);   // #mu^{2}
  fAsymmetric2PE->SetParameter(2, 1);    // #sigma^{2}_{1}
  fAsymmetric2PE->SetParameter(3, 0.5);  // #sigma^{2}_{2}

  // Define total function as sum of 1 PE + 2 PE
  TF1 *fTotal = new TF1("fTotal", asym2Gaussians, 0.5, 3., 9);
  fTotal->SetLineColor(kGreen + 2);
  fTotal->SetLineWidth(4);
  fTotal->SetLineStyle(2);

  // Define parameter array for total functions
  Double_t par1[9];

  // Normalise  hPhotoElectrons->Scale(1.0 / hPhotoElectrons->GetMaximum());
  hPhotoElectrons->Scale(1.0 / hPhotoElectrons->GetMaximum());

  // Fit PE graphs

  // 1PE
  hPhotoElectrons->Fit(fAsymmetric1PE, "RN");
  fAsymmetric1PE->GetParameters(&par1[0]);
  // Print pvalue and reduced chi squared
  std::cout << "\n\n**** FIT RESULT 1 PE peak ****\n\nP value       "
               "      = "
            << fAsymmetric1PE->GetProb() << '\n';
  std::cout << "Reduced chi squared = "
            << fAsymmetric1PE->GetChisquare() / fAsymmetric1PE->GetNDF()
            << "\n\n";

  // 2PE
  hPhotoElectrons->Fit(fAsymmetric2PE, "RN");
  fAsymmetric2PE->GetParameters(&par1[4]);
  // Print pvalue and reduced chi squared
  std::cout << "\n\n**** FIT RESULT 2 PE peak ****\n\nP value       "
               "      = "
            << fAsymmetric2PE->GetProb() << '\n';
  std::cout << "Reduced chi squared = "
            << fAsymmetric2PE->GetChisquare() / fAsymmetric2PE->GetNDF()
            << "\n\n";

  // Total function fit NOT CONSTRAINED

  par1[8] = 0.005;
  fTotal->SetParameters(par1);
  fTotal->SetParNames("N^{1}", "#mu^{1}", "#sigma^{1}_{1}", "#sigma^{1}_{2}",
                      "N^{2}", "#mu^{2}", "#sigma^{2}_{1}", "#sigma^{2}_{2}",
                      "Background");
  TFitResultPtr fitResult = hPhotoElectrons->Fit(fTotal, "S R+ N");

  // Get results
  std::cout
      << "\n\n**** FIT RESULT TOTAL NOT CONSTRAINED #mu ****\n\nP value       "
         "      = "
      << fTotal->GetProb() << '\n';
  std::cout << "Reduced chi squared = "
            << fTotal->GetChisquare() / fTotal->GetNDF() << "\n\n";
  TMatrixD covMatrix = fitResult->GetCorrelationMatrix();
  TMatrixD corMatrix = fitResult->GetCovarianceMatrix();
  std::cout << "\n*** Print covariance matrix ***\n" << std::endl;
  covMatrix.Print();
  std::cout << "\n*** Print correlation matrix ***\n" << std::endl;
  corMatrix.Print();

  // Total function fit CONSTRAINED

  // Define total function as sum of 1 PE + 2 PE
  TF1 *fTotalConst =
      new TF1("fTotalConst", asym2GaussiansConstrained, 0.5, 3., 8);
  fTotalConst->SetLineColor(kOrange + 2);
  fTotalConst->SetLineWidth(4);
  fTotalConst->SetLineStyle(2);

  fTotalConst->SetParameters(par1[0], par1[1], par1[2], par1[3], par1[4],
                             par1[6], par1[7], par1[8]);
  fTotalConst->SetParNames("N^{1}", "#mu^{1}", "#sigma^{1}_{1}",
                           "#sigma^{1}_{2}", "N^{2}", "#sigma^{2}_{1}",
                           "#sigma^{2}_{2}", "Background");
  TFitResultPtr fitResultConst = hPhotoElectrons->Fit(fTotalConst, "S R+ N");

  // Get results
  std::cout
      << "\n\n**** FIT RESULT TOTAL CONSTRAINED #mu ****\n\nP value       "
         "      = "
      << fTotalConst->GetProb() << '\n';
  std::cout << "Reduced chi squared = "
            << fTotalConst->GetChisquare() / fTotalConst->GetNDF() << "\n\n";
  TMatrixD covMatrixConst = fitResultConst->GetCorrelationMatrix();
  TMatrixD corMatrixConst = fitResultConst->GetCovarianceMatrix();
  std::cout << "\n*** Print covariance matrix ***\n" << std::endl;
  covMatrixConst.Print();
  std::cout << "\n*** Print correlation matrix ***\n" << std::endl;
  corMatrixConst.Print();

  // Total function fit CONSTRAINED SAME SIGMAS

  // Define total function as sum of 1 PE + 2 PE
  TF1 *fTotalConstSameSigmas = new TF1(
      "fTotalConstSameSigmas", asym2GaussiansConstrainedSameSigma, 0.5, 3., 6);
  fTotalConstSameSigmas->SetLineColor(kRed);
  fTotalConstSameSigmas->SetLineWidth(4);
  fTotalConstSameSigmas->SetLineStyle(2);

  fTotalConstSameSigmas->SetParameters(par1[0], par1[1], par1[2], par1[3],
                                       par1[4], par1[8]);
  fTotalConstSameSigmas->SetParNames("N^{1}", "#mu^{1}", "#sigma^{1}_{1}",
                                     "#sigma^{1}_{2}", "N^{2}", "Background");
  TFitResultPtr fitResultConstSameSigma =
      hPhotoElectrons->Fit(fTotalConstSameSigmas, "S R+");

  // Get results
  std::cout << "\n\n**** FIT RESULT TOTAL CONSTRAINED SAME SIGMA #mu ****\n\nP "
               "value       "
               "      = "
            << fTotalConstSameSigmas->GetProb() << '\n';
  std::cout << "Reduced chi squared = "
            << fTotalConstSameSigmas->GetChisquare() /
                   fTotalConstSameSigmas->GetNDF()
            << "\n\n";
  TMatrixD covMatrixConstSameSigma =
      fitResultConstSameSigma->GetCorrelationMatrix();
  TMatrixD corMatrixConstSameSigma =
      fitResultConstSameSigma->GetCovarianceMatrix();
  std::cout << "\n*** Print covariance matrix ***\n" << std::endl;
  covMatrixConstSameSigma.Print();
  std::cout << "\n*** Print correlation matrix ***\n" << std::endl;
  corMatrixConstSameSigma.Print();
}

// Asymmetric gaussian functions

Double_t asym2GaussiansExpo(Double_t *x, Double_t *par) {
  // par[0] = N^{1}
  // par[1] = #mu^{1}
  // par[2] = #sigma^{1}_{1}
  // par[3] = #sigma^{1}_{2}
  // par[4] = N^{2}
  // par[5] = #mu^{2}
  // par[6] = #sigma^{2}_{1}
  // par[7] = #sigma^{2}_{2}
  // par[8] = Constant
  // par[9] = Slope
  // par[10] = Background

  return asymGaussians(x, &par[0]) + asymGaussians(x, &par[4]) +
         TMath::Exp(par[8] + par[9] * x[0]) + par[10];
}

Double_t asym2GaussiansExpoConstrained(Double_t *x, Double_t *par) {
  // par[0] = N^{1}
  // par[1] = #mu^{1}
  // par[2] = #sigma^{1}_{1}
  // par[3] = #sigma^{1}_{2}
  // par[4] = N^{2}
  // par[5] = #sigma^{2}_{1}
  // par[6] = #sigma^{2}_{2}
  // par[7] = Constant
  // par[8] = Slope
  // par[9] = Background

  return asymGaussians(x, new Double_t[4]{par[0], par[1], par[2], par[3]}) +
         asymGaussians(x, new Double_t[4]{par[4], 2 * par[1], par[5], par[6]}) +
         TMath::Exp(par[7] + par[8] * x[0]) + par[9];
}

Double_t asym2GaussiansExpoConstrainedSameSigma(Double_t *x, Double_t *par) {
  // par[0] = N^{1}
  // par[1] = #mu^{1}
  // par[2] = #sigma^{1}_{1}
  // par[3] = #sigma^{1}_{2}
  // par[4] = N^{2}
  // par[5] = Constant
  // par[6] = Slope
  // par[7] = Background

  return asymGaussians(x, new Double_t[4]{par[0], par[1], par[2], par[3]}) +
         asymGaussians(x, new Double_t[4]{par[4], 2 * par[1], par[2], par[3]}) +
         TMath::Exp(par[5] + par[6] * x[0]) + par[7];
}

// Group functions for fitting PE histo

void fitPEHistoExp(TH1F *hPhotoElectrons) {
  // Import user defined function asymmetric gaussians for 1 PE
  TF1 *fAsymmetric1PE = new TF1("fAsymmetric1PE", asymGaussians, 0.5, 1.9, 4);
  fAsymmetric1PE->SetLineColor(kRed);
  fAsymmetric1PE->SetLineWidth(4);
  fAsymmetric1PE->SetLineStyle(2);
  fAsymmetric1PE->SetParNames("N^{1}", "#mu^{1}", "#sigma^{1}_{1}",
                              "#sigma^{1}_{2}");
  fAsymmetric1PE->SetParameter(0, 0.3);  // N^{1}
  fAsymmetric1PE->SetParameter(1, 1.);   // #mu^{1}
  fAsymmetric1PE->SetParameter(2, 0.5);  // #sigma^{1}_{1}
  fAsymmetric1PE->SetParameter(3, 0.5);  // #sigma^{1}_{2}

  // Import user defined function asymmetric gaussians for 2 PE
  TF1 *fAsymmetric2PE = new TF1("fAsymmetric2PE", asymGaussians, 1.9, 2.6, 4);
  fAsymmetric2PE->SetLineColor(kOrange + 2);
  fAsymmetric2PE->SetLineWidth(4);
  fAsymmetric2PE->SetLineStyle(2);
  fAsymmetric2PE->SetParNames("N^{2}", "#mu^{2}", "#sigma^{2}_{1}",
                              "#sigma^{2}_{2}");
  fAsymmetric2PE->SetParameter(0, 0.1);  // N^{2}
  fAsymmetric2PE->SetParameter(1, 2.2);  // #mu^{2}
  fAsymmetric2PE->SetParameter(2, 1);    // #sigma^{2}_{1}
  fAsymmetric2PE->SetParameter(3, 0.5);  // #sigma^{2}_{2}

  // Define expo function for noise fit
  TF1 *fExpo = new TF1("fExpo", "expo", 0.5, 3.);
  fExpo->SetLineColor(kGreen + 2);
  fExpo->SetLineWidth(4);
  fExpo->SetLineStyle(2);
  fExpo->SetParameter(0, -0.4);  // Constant
  fExpo->SetParameter(1, -1.2);  // Slope

  // Define total function as sum of 1 PE + 2 PE
  TF1 *fTotal = new TF1("fTotal", asym2GaussiansExpo, 0.5, 3., 11);
  fTotal->SetLineColor(kGreen + 2);
  fTotal->SetLineWidth(4);
  fTotal->SetLineStyle(2);

  // Define parameter array for total functions
  Double_t par1[11];

  // Normalise  hPhotoElectrons->Scale(1.0 / hPhotoElectrons->GetMaximum());
  hPhotoElectrons->Scale(1.0 / hPhotoElectrons->GetMaximum());

  // Fit PE graphs

  // 1PE
  hPhotoElectrons->Fit(fAsymmetric1PE, "RN");
  fAsymmetric1PE->GetParameters(&par1[0]);
  // Print pvalue and reduced chi squared
  std::cout << "\n\n**** FIT RESULT 1 PE peak ****\n\nP value       "
               "      = "
            << fAsymmetric1PE->GetProb() << '\n';
  std::cout << "Reduced chi squared = "
            << fAsymmetric1PE->GetChisquare() / fAsymmetric1PE->GetNDF()
            << "\n\n";

  // 2PE
  hPhotoElectrons->Fit(fAsymmetric2PE, "RN");
  fAsymmetric2PE->GetParameters(&par1[4]);
  // Print pvalue and reduced chi squared
  std::cout << "\n\n**** FIT RESULT 2 PE peak ****\n\nP value       "
               "      = "
            << fAsymmetric2PE->GetProb() << '\n';
  std::cout << "Reduced chi squared = "
            << fAsymmetric2PE->GetChisquare() / fAsymmetric2PE->GetNDF()
            << "\n\n";

  // Exponential noise
  hPhotoElectrons->Fit(fExpo, "N", "", 0., 0.2);
  fExpo->GetParameters(&par1[8]);
  // Print pvalue and reduced chi squared
  std::cout << "\n\n**** FIT RESULT EXPO NOISE ****\n\nP value       "
               "      = "
            << fExpo->GetProb() << '\n';
  std::cout << "Reduced chi squared = "
            << fExpo->GetChisquare() / fExpo->GetNDF() << "\n\n";

  // Total function fit NOT CONSTRAINED

  par1[10] = 0.001;
  fTotal->SetParameters(par1);
  fTotal->SetParNames("N^{1}", "#mu^{1}", "#sigma^{1}_{1}", "#sigma^{1}_{2}",
                      "N^{2}", "#mu^{2}", "#sigma^{2}_{1}", "#sigma^{2}_{2}",
                      "Constant", "Slope");
  fTotal->SetParName(10, "Background");
  TFitResultPtr fitResult = hPhotoElectrons->Fit(fTotal, "S R+");

  // Get results
  std::cout
      << "\n\n**** FIT RESULT TOTAL NOT CONSTRAINED #mu ****\n\nP value       "
         "      = "
      << fTotal->GetProb() << '\n';
  std::cout << "Reduced chi squared = "
            << fTotal->GetChisquare() / fTotal->GetNDF() << "\n\n";
  TMatrixD covMatrix = fitResult->GetCorrelationMatrix();
  TMatrixD corMatrix = fitResult->GetCovarianceMatrix();
  std::cout << "\n*** Print covariance matrix ***\n" << std::endl;
  covMatrix.Print();
  std::cout << "\n*** Print correlation matrix ***\n" << std::endl;
  corMatrix.Print();

  // Total function fit CONSTRAINED

  // Define total function as sum of 1 PE + 2 PE
  TF1 *fTotalConst =
      new TF1("fTotalConst", asym2GaussiansExpoConstrained, 0.5, 3., 10);
  fTotalConst->SetLineColor(kOrange + 2);
  fTotalConst->SetLineWidth(4);
  fTotalConst->SetLineStyle(2);

  fTotalConst->SetParameters(par1[0], par1[1], par1[2], par1[3], par1[4],
                             par1[6], par1[7], par1[8], par1[9], par1[10]);
  fTotalConst->SetParNames("N^{1}", "#mu^{1}", "#sigma^{1}_{1}",
                           "#sigma^{1}_{2}", "N^{2}", "#sigma^{2}_{1}",
                           "#sigma^{2}_{2}", "Constant", "Slope", "Background");
  TFitResultPtr fitResultConst = hPhotoElectrons->Fit(fTotalConst, "S R+");

  // Get results
  std::cout
      << "\n\n**** FIT RESULT TOTAL CONSTRAINED #mu ****\n\nP value       "
         "      = "
      << fTotalConst->GetProb() << '\n';
  std::cout << "Reduced chi squared = "
            << fTotalConst->GetChisquare() / fTotalConst->GetNDF() << "\n\n";
  TMatrixD covMatrixConst = fitResultConst->GetCorrelationMatrix();
  TMatrixD corMatrixConst = fitResultConst->GetCovarianceMatrix();
  std::cout << "\n*** Print covariance matrix ***\n" << std::endl;
  covMatrixConst.Print();
  std::cout << "\n*** Print correlation matrix ***\n" << std::endl;
  corMatrixConst.Print();

  // Total function fit CONSTRAINED SAME SIGMAS

  // Define total function as sum of 1 PE + 2 PE
  TF1 *fTotalConstSameSigmas =
      new TF1("fTotalConstSameSigmas", asym2GaussiansExpoConstrainedSameSigma,
              0.5, 3., 8);
  fTotalConstSameSigmas->SetLineColor(kRed);
  fTotalConstSameSigmas->SetLineWidth(4);
  fTotalConstSameSigmas->SetLineStyle(2);

  fTotalConstSameSigmas->SetParameters(par1[0], par1[1], par1[2], par1[3],
                                       par1[4], par1[8], par1[9], par1[10]);
  fTotalConstSameSigmas->SetParNames("N^{1}", "#mu^{1}", "#sigma^{1}_{1}",
                                     "#sigma^{1}_{2}", "N^{2}", "Constant",
                                     "Slope", "Background");
  TFitResultPtr fitResultConstSameSigma =
      hPhotoElectrons->Fit(fTotalConstSameSigmas, "S R+");

  // Get results
  std::cout << "\n\n**** FIT RESULT TOTAL CONSTRAINED SAME SIGMA #mu ****\n\nP "
               "value       "
               "      = "
            << fTotalConstSameSigmas->GetProb() << '\n';
  std::cout << "Reduced chi squared = "
            << fTotalConstSameSigmas->GetChisquare() /
                   fTotalConstSameSigmas->GetNDF()
            << "\n\n";
  TMatrixD covMatrixConstSameSigma =
      fitResultConstSameSigma->GetCorrelationMatrix();
  TMatrixD corMatrixConstSameSigma =
      fitResultConstSameSigma->GetCovarianceMatrix();
  std::cout << "\n*** Print covariance matrix ***\n" << std::endl;
  covMatrixConstSameSigma.Print();
  std::cout << "\n*** Print correlation matrix ***\n" << std::endl;
  corMatrixConstSameSigma.Print();
}

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

// This function analyses the waveform by building areaVStime, Noise, PE
// counts and pulseWidth histos
void waveformAnalysis(
    AreaConvFactor areaConv = Transm,
    std::string infileName =
        "./data/45Degrees/0.20mm/"
        "DataF_CH0@DT5730S_59483_run_45_0.2_TRANSM_REFL.txt",
    std::string rootFileName = "./rootFiles/45Degrees0.20mm/wA2.root") {
  // To avoid reloading manually if .so is present
  R__LOAD_LIBRARY(waveformAnalysisPos_cpp.so);

  // Area conversion factor (current assumption is 1 PE = 11000 ADC*ns)
  // For Transm PMT = 4000.
  // For Refl PMT = 12500.
  auto const areaConvFactor = static_cast<double>(areaConv);

  // Variables used later
  double const samplePeriod = 2.0;  // In [ns]
  std::ifstream infile(infileName.c_str());
  std::string line;
  std::vector<double> colours{1, 3, 4, 5, 6, 7, 8, 9};  // Colour vector
  TMultiGraph *mg = new TMultiGraph();
  TMultiGraph *mgSuperimposed = new TMultiGraph();
  std::vector<TGraph *> graphs{};
  std::vector<TGraph *> graphsSuperimposed{};
  std::map<int, double> map{};
  int row{0};
  int rowSum{0};

  // Select random generator seed for colours based on current time
  srand(time(NULL));

  // Creating TFile
  TFile *file1 = new TFile(rootFileName.c_str(), "RECREATE");

  // Define histograms
  TH2F *hAreaVsTime = new TH2F("hAreaVsTime",
                               "Pulse area vs relative time; Relative time "
                               "peak [ns]; Area [ADC #times ns]",
                               40, 100., 460., 100, 0., 100000.);
  TH1F *hNoise = new TH1F("hNoise", "Noise distribution; ADC counts; Counts",
                          30, 2755, 3000);
  TH1F *hPhotoElectrons = new TH1F(
      "hPE", "Pulse area distribution; Area [PE]; Normalized counts", 90, 0, 6);
  TH1F *hWidth =
      new TH1F("hWidth", "Width distribution; Width [ns]; Counts", 20, 2, 50);
  TH1F *hPETrigger = new TH1F(
      "hPETrigger", "Trigger region; Area [PE]; Normalized counts", 150, 0, 6);
  TH1F *hPEPreTrigger =
      new TH1F("hPEPreTrigger",
               "Pre trigger region; Area [PE]; Normalized counts", 15, 0, 2);
  TH1F *hPEPostTrigger1 = new TH1F(
      "hPEPostTrigger1", "Post trigger region 1; Area [PE]; Normalized counts",
      15, 0., 1.6);
  TH1F *hPEPostTrigger2 = new TH1F(
      "hPEPostTrigger2", "Post trigger region 2; Area [PE]; Normalized counts",
      15, 0, 2.5);

  // Plotting parameters
  int const nPulseParam{14};
  Parameter pulsePar[nPulseParam] = {
      {"Width [ns]", 2., 50.},
      {"Area [ADC #times ns]", 0., 30000.},
      {"Area [PE]", 0, 6},
      {"Relative peak time [ns]", 100., 300.},
      {"Noise RMS [ADC]", 8034, 8040},
      {"Peak amplitude [ADC]", 7900, 15500},
      {"Height over width [ADC/ns]", 0., 2000},
      {"Peak fraction position", 0., 1.},
      {"Area over full time width [ADC]", 0., 2000.},
      {"Rise time [ns]", 0., 20.},
      {"FWHM [ns]", 0., 20.},
      {"90% area time [ns]", 0., 35.},
      {"Neg area/overall area", 0., 0.05},
      {"Neg counts/overall counts", 0., 1.}};
  int const nBins{150};
  TH1F *hPulsePar[nPulseParam];
  TH2F *h2PulsePar[nPulseParam][nPulseParam];
  TGraph *gPulsePar[nPulseParam][nPulseParam];

  // Initialisation loop for plotting parameters
  // par_i labels the rows and par_j labels the columns
  for (int par_i = 0; par_i < nPulseParam; ++par_i) {
    hPulsePar[par_i] =
        new TH1F(Form("h1PulsePar_%d", par_i),
                 Form("%s; %s; Counts", pulsePar[par_i].label.c_str(),
                      pulsePar[par_i].label.c_str()),
                 nBins, pulsePar[par_i].min, pulsePar[par_i].max);
    for (int par_j = 0; par_j < nPulseParam; ++par_j) {
      if (par_i > par_j) {
        gPulsePar[par_i][par_j] = new TGraph();
        gPulsePar[par_i][par_j]->SetTitle(
            Form("%s vs %s;%s;%s", pulsePar[par_i].label.c_str(),
                 pulsePar[par_j].label.c_str(), pulsePar[par_j].label.c_str(),
                 pulsePar[par_i].label.c_str()));

      } else if (par_i < par_j) {
        h2PulsePar[par_i][par_j] = new TH2F(
            Form("h2PulsePar_%d_%d", par_i, par_j),
            Form("%s vs %s;%s;%s", pulsePar[par_i].label.c_str(),
                 pulsePar[par_j].label.c_str(), pulsePar[par_j].label.c_str(),
                 pulsePar[par_i].label.c_str()),
            nBins, pulsePar[par_j].min, pulsePar[par_j].max, nBins,
            pulsePar[par_i].min, pulsePar[par_i].max);
      }
    }
  }

  // Loop over rows (waveforms) to extract pulse sum graph
  while (std::getline(infile, line)) {
    // Control over analysed rows
    if (rowSum < nMinAnalysedRows) {
      ++rowSum;
      continue;
    }
    if (rowSum >= nMaxAnalysedRows) {
      break;
    }

    // Defining loop variables
    std::stringstream ss(line);
    std::string item;
    std::vector<double> samples;
    double timestamp = 0.;
    int column = 1;

    // Loop over columns
    while (std::getline(ss, item, '\t')) {
      if (item.empty()) continue;
      if (column == 3)
        timestamp = std::stod(item);
      else if (column >= 7)
        samples.push_back(std::stod(item));
      ++column;
    }

    // Creating WaveformAnalysis object
    WaveformAnalysisPos wf(samples, timestamp, samplePeriod);
    const auto &pulses = wf.getPulses();

    for (size_t i = 0; i < pulses.size(); ++i) {
      const auto &p = pulses[i];

      // Insert selections on pulses
      if (p.peakValue > 15900.) {
        continue;
      }

      // Generate sum of pulses using a STL map. Here the keys of the map
      // correspond to the time values of pulses, while values correpond to the
      // voltage/ADC counts. Contrary to std::vector<T>, std::map<T1,T2>'s
      // operator[] is more secure: when you have an empty map, with no
      // specified size, and you access map[X], the map safely creates the
      // key-value pair (X,0.0) - in this case I put 0.0 because it is the
      // default initialiser for double values
      for (int j = 0; j < static_cast<int>(p.times.size()); ++j) {
        map[p.times[j]] += p.values[j];
      }
    }
    ++rowSum;
  }

  // Create canvas for summing pulses
  TCanvas *cPulseSum = new TCanvas("cPulseSum", "Pulse sum", 1500, 700);
  std::vector<double> xValues{};
  std::vector<double> yValues{};
  xValues.reserve(map.size());
  yValues.reserve(map.size());
  for (auto const &[key, value] : map) {
    xValues.push_back(key);
    yValues.push_back(value);
  }

  // Create graph for summing pulses
  TGraph *gPulseSum = new TGraph(map.size(), xValues.data(), yValues.data());
  gPulseSum->SetTitle("Pulse sum; Time since \"trigger\" [ns]; ADC Counts");
  gPulseSum->SetLineColor(kBlue);
  gPulseSum->SetLineWidth(1);
  gPulseSum->SetMarkerColor(kBlack);
  gPulseSum->SetMarkerStyle(20);
  gPulseSum->SetMarkerSize(1);

  // gPulseSum relevant paramters
  auto const maxId{TMath::LocMax(gPulseSum->GetN(), gPulseSum->GetY())};
  auto const sigmaId{gPulseSum->GetRMS()};
  double const startPoint{gPulseSum->GetX()[maxId] - 0.5 * gPulseSum->GetRMS()};
  double const endPoint{gPulseSum->GetX()[maxId] + 0.5 * gPulseSum->GetRMS()};

  // Create fit function for summed pulses
  TF1 *fGaus = new TF1("fGaus", "gaus", startPoint, endPoint);
  fGaus->SetLineColor(kRed);
  fGaus->SetLineWidth(4);
  fGaus->SetLineStyle(2);
  fGaus->SetParameter(0, gPulseSum->GetY()[maxId]);  // Amplitude
  fGaus->SetParameter(1, gPulseSum->GetX()[maxId]);  // Mean
  fGaus->SetParameter(2, sigmaId);                   // Sigma
  gPulseSum->Fit(fGaus, "R");

  // Save gPulseSum fit relevant parameters
  double const meanSum{fGaus->GetParameter(1)};
  double const sigmaSum{fGaus->GetParameter(2)};

  // Variable for analysis in trigger region in 1 single file
  int pulseCounter{};
  int pulseCounterTriggerRegion{};
  double totTrigArea{};
  double numTrigPE{};

  // Below gaussian fit on "Pulse sum" graph is used (0.5 sigmas)
  double const triggerStart{meanSum - 0.5 * sigmaSum};
  double const triggerEnd{meanSum + 0.5 * sigmaSum};

  // Variable for analysis in pre-trigger region in 1 single file
  int pulseCounterPreTriggerRegion{};
  double totPreTrigArea{};
  double numPreTrigPE{};
  double const preTriggerStart{meanSum - 4.5 * sigmaSum};
  double const preTriggerEnd{meanSum - 4.5 * sigmaSum + 20};

  // Variable for analysis in post-trigger region 1 in 1 single file
  int pulseCounterPostTriggerRegion1{};
  double totPostTrigArea1{};
  double numPostTrigPE1{};
  double const postTriggerStart1{meanSum + 3 * sigmaSum};
  double const postTriggerEnd1{meanSum + 3 * sigmaSum + 20};

  // Variable for analysis in post-trigger region 2 in 1 single file
  int pulseCounterPostTriggerRegion2{};
  double totPostTrigArea2{};
  double numPostTrigPE2{};
  double const postTriggerStart2{meanSum + 3 * sigmaSum + 40};
  double const postTriggerEnd2{meanSum + 3 * sigmaSum + 60};

  // Prepare to read again the file
  infile.clear();               // Clear the failed state of the stream
  infile.seekg(0, infile.beg);  // Seek to the first character in the file

  // Loop over rows (waveforms)
  while (std::getline(infile, line)) {
    // Control over analysed rows
    if (row < nMinAnalysedRows) {
      ++row;
      continue;
    }
    if (row >= nMaxAnalysedRows) {
      break;
    }

    // Defining loop variables
    std::stringstream ss(line);
    std::string item;
    std::vector<double> samples;
    double timestamp = 0.;
    int column = 1;

    // Loop over columns
    while (std::getline(ss, item, '\t')) {
      if (item.empty()) {
        continue;
      }
      if (column == 3)
        timestamp = std::stod(item);
      else if (column >= 7)
        samples.push_back(std::stod(item));
      ++column;
    }

    // Creating WaveformAnalysis object
    WaveformAnalysisPos wf(samples, timestamp, samplePeriod);

    // Print waveform properties
    std::cout << std::fixed
              << std::setprecision(2);  // Round to 2 decimal place
    std::cout << "\n********** Waveform n. " << row + 1 << " **********\n";
    std::cout << "Timestamp        = " << wf.getTimeStamp() << " ns\n";
    std::cout << "Baseline         = " << wf.getBaseline() << " ADC counts\n";
    std::cout << "Sample period    = " << wf.getSamplePeriod() << " ns\n";

    // Get pulse vector from each single waveform
    const auto &pulses = wf.getPulses();
    std::cout << "Number of Pulses = " << pulses.size() << '\n';

    // Fill noise information
    int const nBaselineSamples{50};
    std::for_each(samples.begin(), samples.begin() + nBaselineSamples,
                  [&](double sample) { hNoise->Fill(sample); });

    // Print pulse properties
    for (size_t i = 0; i < pulses.size(); ++i) {
      const auto &p = pulses[i];

      // Insert selections on pulses
      if (p.peakValue > 15900.) {
        continue;
      }

      // Count total selected pulses
      ++pulseCounter;

      // Find area in trigger region
      if ((p.peakTime - wf.getTimeStamp()) >= triggerStart &&
          (p.peakTime - wf.getTimeStamp()) <= triggerEnd) {
        totTrigArea += p.area;
        ++pulseCounterTriggerRegion;
        hPETrigger->Fill(p.area / areaConvFactor);
      }

      // Find area in pre-trigger region
      if ((p.peakTime - wf.getTimeStamp()) >= preTriggerStart &&
          (p.peakTime - wf.getTimeStamp()) <= preTriggerEnd) {
        totPreTrigArea += p.area;
        ++pulseCounterPreTriggerRegion;
        hPEPreTrigger->Fill(p.area / areaConvFactor);
      }

      // Find area in post-trigger region 1
      if ((p.peakTime - wf.getTimeStamp()) >= postTriggerStart1 &&
          (p.peakTime - wf.getTimeStamp()) <= postTriggerEnd1) {
        totPostTrigArea1 += p.area;
        ++pulseCounterPostTriggerRegion1;
        hPEPostTrigger1->Fill(p.area / areaConvFactor);
      }

      // Find area in post-trigger region 2
      if ((p.peakTime - wf.getTimeStamp()) >= postTriggerStart2 &&
          (p.peakTime - wf.getTimeStamp()) <= postTriggerEnd2) {
        totPostTrigArea2 += p.area;
        ++pulseCounterPostTriggerRegion2;
        hPEPostTrigger2->Fill(p.area / areaConvFactor);
      }

      // Params of interest
      double heightOverWidth{p.peakValue / (p.endTime - p.startTime)};
      double peakFractionPos{(p.peakTime - p.startTime) /
                             (p.endTime - p.startTime)};
      double areaOverFullTime{p.area / ((p.endTime - p.startTime))};

      std::cout << "\n  *** Pulse n. " << i + 1 << " ***\n\n";
      std::cout << "  Overall start time           = " << p.startTime
                << " ns\n";
      std::cout << "  Overall end time             = " << p.endTime << " ns\n";
      std::cout << "  Overall peak time            = " << p.peakTime << " ns\n";
      std::cout << "  Relative start time          = "
                << p.startTime - wf.getTimeStamp() << " ns\n";
      std::cout << "  Relative end time            = "
                << p.endTime - wf.getTimeStamp() << " ns\n";
      std::cout << "  Relative peak time           = "
                << p.peakTime - wf.getTimeStamp() << " ns\n";
      std::cout << "  Peak time since startPulse   = "
                << p.peakTime - wf.getTimeStamp() - p.times[0] << " ns\n";
      std::cout << "  Peak value                   = " << p.peakValue
                << " ADC\n";
      std::cout << "  Width                        = "
                << p.endTime - p.startTime << " ns\n";
      std::cout << "  Rise time                    = " << p.riseTime << " ns\n";
      std::cout << "  FWHM                         = " << p.FWHMTime << " ns\n";
      std::cout << "  90% area time                = " << p.areaFractionTime
                << " ns\n";
      std::cout << "  Height over width            = " << heightOverWidth
                << " ADC/ns\n";
      std::cout << "  Peak fraction pos.           = " << peakFractionPos
                << '\n';
      std::cout << "  Area / full width            = " << areaOverFullTime
                << " ADC\n";
      std::cout << "  Area                         = " << p.area << " ADC*ns\n";
      std::cout << "  Area in PE                   = "
                << p.area / areaConvFactor << " PE\n";
      std::cout << "  Negative/overall area frac   = " << p.negFracArea
                << " \n";
      std::cout << "  Negative/overall counts      = " << p.negFrac << " \n";

      // Generate a random number between 0 and 7 (used for colour indices)
      int randIndex = rand() % 8;

      // Create vector to superimpose pulses
      std::vector<double> superimposedTimes = p.times;
      double shift{superimposedTimes[0]};
      for (int timeId{}; timeId < static_cast<int>(superimposedTimes.size());
           ++timeId) {
        superimposedTimes[timeId] -= shift;
      }

      // Plot each pulse using a graph object
      TGraph *g = new TGraph(p.times.size(), p.times.data(), p.values.data());
      if ((p.peakTime - wf.getTimeStamp()) >= triggerStart &&
          (p.peakTime - wf.getTimeStamp()) <= triggerEnd) {
        g->SetLineColor(kRed);
      } else if ((p.peakTime - wf.getTimeStamp()) >= preTriggerStart &&
                 (p.peakTime - wf.getTimeStamp()) <= preTriggerEnd) {
        g->SetLineColor(kRed);
      } else if (((p.peakTime - wf.getTimeStamp()) >= postTriggerStart1 &&
                  (p.peakTime - wf.getTimeStamp()) <= postTriggerEnd1) ||
                 ((p.peakTime - wf.getTimeStamp()) >= postTriggerStart2 &&
                  (p.peakTime - wf.getTimeStamp()) <= postTriggerEnd2)) {
        g->SetLineColor(kRed);
      } else {
        g->SetLineColor(colours[randIndex]);
      }
      g->SetLineWidth(1);
      g->SetMarkerColor(kBlack);
      g->SetMarkerStyle(20);
      g->SetMarkerSize(1);
      g->SetTitle(Form("Pulse %d; Time [ns]; ADC counts", pulseCounter));
      graphs.push_back(g);

      // Superimpose pulses from riseTime
      TGraph *gSuperimposed = new TGraph(
          superimposedTimes.size(), superimposedTimes.data(), p.values.data());
      if ((p.peakTime - wf.getTimeStamp()) >= triggerStart &&
          (p.peakTime - wf.getTimeStamp()) <= triggerEnd) {
        gSuperimposed->SetLineColor(kRed);
      } else if ((p.peakTime - wf.getTimeStamp()) >= preTriggerStart &&
                 (p.peakTime - wf.getTimeStamp()) <= preTriggerEnd) {
        gSuperimposed->SetLineColor(kRed);
      } else if (((p.peakTime - wf.getTimeStamp()) >= postTriggerStart1 &&
                  (p.peakTime - wf.getTimeStamp()) <= postTriggerEnd1) ||
                 ((p.peakTime - wf.getTimeStamp()) >= postTriggerStart2 &&
                  (p.peakTime - wf.getTimeStamp()) <= postTriggerEnd2)) {
        gSuperimposed->SetLineColor(kRed);
      } else {
        gSuperimposed->SetLineColor(colours[randIndex]);
      }
      gSuperimposed->SetLineWidth(1);
      gSuperimposed->SetMarkerColor(kBlack);
      gSuperimposed->SetMarkerStyle(20);
      gSuperimposed->SetMarkerSize(1);
      gSuperimposed->SetTitle(
          Form("Pulse %d; Time [ns]; ADC counts", pulseCounter));
      graphsSuperimposed.push_back(gSuperimposed);

      // Fill pulse information
      hAreaVsTime->Fill(p.peakTime - wf.getTimeStamp(), p.area);
      hWidth->Fill(p.endTime - p.startTime);

      // Convert area into PE
      double areaInPE = p.area / areaConvFactor;
      hPhotoElectrons->Fill(areaInPE);

      // Plot all parameters against each other

      // Gather params of interest from each pulse and noise from waveform
      double parValues[nPulseParam] = {p.endTime - p.startTime,
                                       p.area,
                                       p.area / areaConvFactor,
                                       p.peakTime - wf.getTimeStamp(),
                                       wf.getBaseline(),
                                       p.peakValue,
                                       heightOverWidth,
                                       peakFractionPos,
                                       areaOverFullTime,
                                       p.riseTime,
                                       p.FWHMTime,
                                       p.areaFractionTime,
                                       p.negFracArea,
                                       p.negFrac};
      for (int par_i = 0; par_i < nPulseParam; ++par_i) {
        hPulsePar[par_i]->Fill(parValues[par_i]);

        for (int par_j = 0; par_j < nPulseParam; ++par_j) {
          if (par_i > par_j) {
            {
              gPulsePar[par_i][par_j]->AddPoint(parValues[par_j],
                                                parValues[par_i]);
            }
          } else if (par_i < par_j) {
            h2PulsePar[par_i][par_j]->Fill(parValues[par_j], parValues[par_i]);
          }
        }
      }
    }
    ++row;
  }

  // Print information for each region of interest rounding to 5 decimal place
  std::cout << std::fixed << std::setprecision(5);

  // Correctly define the total number of PE per region
  numTrigPE = totTrigArea / (areaConvFactor);
  numPreTrigPE = totPreTrigArea / (areaConvFactor);
  numPostTrigPE1 = totPostTrigArea1 / (areaConvFactor);
  numPostTrigPE2 = totPostTrigArea2 / (areaConvFactor);

  // Correctly define the average number of PE per pulse
  double numTrigPEPuls = numTrigPE / (pulseCounterTriggerRegion);
  double numPreTrigPEPuls = numPreTrigPE / (pulseCounterPreTriggerRegion);
  double numPostTrigPE1Puls = numPostTrigPE1 / (pulseCounterPostTriggerRegion1);
  double numPostTrigPE2Puls = numPostTrigPE2 / (pulseCounterPostTriggerRegion2);

  // Define time intervals
  auto deltaTrig = (triggerEnd - triggerStart);
  auto deltaPreTrig = (preTriggerEnd - preTriggerStart);
  auto deltaPostTrig1 = (postTriggerEnd1 - postTriggerStart1);
  auto deltaPostTrig2 = (postTriggerEnd2 - postTriggerStart2);

  // Define rates per region
  double rateTrig = numTrigPE / deltaTrig;
  double ratePreTrig = numPreTrigPE / deltaPreTrig;
  double ratePostTrig1 = numPostTrigPE1 / deltaPostTrig1;
  double ratePostTrig2 = numPostTrigPE2 / deltaPostTrig2;

  // Saving region variables into a TTree
  TTree *tree = new TTree("variablesRegion", "Variables for each region");
  tree->Branch("numPreTrigPE", &numPreTrigPE);
  tree->Branch("numTrigPE", &numTrigPE);
  tree->Branch("numPostTrigPE1", &numPostTrigPE1);
  tree->Branch("numPostTrigPE2", &numPostTrigPE2);
  tree->Branch("numPreTrigPEPuls", &numPreTrigPEPuls);
  tree->Branch("numTrigPEPuls", &numTrigPEPuls);
  tree->Branch("numPostTrigPE1Puls", &numPostTrigPE1Puls);
  tree->Branch("numPostTrigPE2Puls", &numPostTrigPE2Puls);
  tree->Branch("deltaPreTrig", &deltaPreTrig);
  tree->Branch("deltaTrig", &deltaTrig);
  tree->Branch("deltaPostTrig1", &deltaPostTrig1);
  tree->Branch("deltaPostTrig2", &deltaPostTrig2);
  tree->Fill();

  // Printing region information

  std::cout << "\n\n *** INFORMATION ON TOTAL AREA AND NUMBER OF "
               "PHOTOELECTRONS IN PRE-TRIGGER REGION "
            << '[' << preTriggerStart << ',' << preTriggerEnd << "] ns ***\n";
  std::cout << " Pulses region / total  = " << pulseCounterPreTriggerRegion
            << " / " << pulseCounter << " = "
            << (static_cast<double>(pulseCounterPreTriggerRegion) /
                static_cast<double>(pulseCounter))
            << '\n';
  std::cout << " Total area             = " << totPreTrigArea << " [ADC*ns]\n";
  std::cout << " Total number of PE     = " << numPreTrigPE << " PE\n";
  std::cout << " Number of PE per pulse = " << numPreTrigPEPuls
            << " PE/pulse\n";
  std::cout << " Rate                   = " << ratePreTrig << " PE/ns\n";

  std::cout << "\n\n *** INFORMATION ON TOTAL AREA AND NUMBER OF "
               "PHOTOELECTRONS IN TRIGGER REGION "
            << '[' << triggerStart << ',' << triggerEnd << "] ns ***\n";
  std::cout << " Pulses region / total  = " << pulseCounterTriggerRegion
            << " / " << pulseCounter << " = "
            << (static_cast<double>(pulseCounterTriggerRegion) /
                static_cast<double>(pulseCounter))
            << '\n';
  std::cout << " Total area             = " << totTrigArea << " [ADC*ns]\n";
  std::cout << " Total number of PE     = " << numTrigPE << " PE\n";
  std::cout << " Number of PE per pulse = " << numTrigPEPuls << " PE/pulse\n";
  std::cout << " Rate                   = " << rateTrig << " PE/ns\n";

  std::cout << "\n\n *** INFORMATION ON TOTAL AREA AND NUMBER OF "
               "PHOTOELECTRONS IN POST-TRIGGER REGION 1 "
            << '[' << postTriggerStart1 << ',' << postTriggerEnd1
            << "] ns ***\n";
  std::cout << " Pulses region / total  = " << pulseCounterPostTriggerRegion1
            << " / " << pulseCounter << " = "
            << (static_cast<double>(pulseCounterPostTriggerRegion1) /
                static_cast<double>(pulseCounter))
            << '\n';
  std::cout << " Total area             = " << totPostTrigArea1
            << " [ADC*ns]\n";
  std::cout << " Total number of PE     = " << numPostTrigPE1 << " PE\n";
  std::cout << " Number of PE per pulse = " << numPostTrigPE1Puls
            << " PE/pulse\n";
  std::cout << " Rate                   = " << ratePostTrig1 << " PE/ns\n";

  std::cout << "\n\n *** INFORMATION ON TOTAL AREA AND NUMBER OF "
               "PHOTOELECTRONS IN POST-TRIGGER REGION 2 "
            << '[' << postTriggerStart2 << ',' << postTriggerEnd2
            << "] ns ***\n";
  std::cout << " Pulses region / total  = " << pulseCounterPostTriggerRegion2
            << " / " << pulseCounter << " = "
            << (static_cast<double>(pulseCounterPostTriggerRegion2) /
                static_cast<double>(pulseCounter))
            << '\n';
  std::cout << " Total area             = " << totPostTrigArea2
            << " [ADC*ns]\n";
  std::cout << " Total number of PE     = " << numPostTrigPE2 << " PE\n";
  std::cout << " Number of PE per pulse = " << numPostTrigPE2Puls
            << " PE/pulse\n";
  std::cout << " Rate                   = " << ratePostTrig2 << " PE/ns\n";

  // FITTING FUNCTIONS SPACE
  std::cout << "\n\n\n************************************\n";
  std::cout << "****** FITTING FUNCTION SPACE ******\n";
  std::cout << "************************************\n\n\n";

  // Fit noise with gaussian function
  hNoise->Fit("gaus");

  setFitStyle();

  // Fit PE histograms with several functions
  gPad->Update();
  fitPEHistoNoExp(hPhotoElectrons);

  // Draw summed pulses
  gPad->Update();
  cPulseSum->cd();
  gPad->Update();
  gPulseSum->Draw("ALP");

  // Print pulse Sum fit info
  std::cout << "\n Print pulse sum fit info:\n";
  std::cout << "  Amplitude         = " << fGaus->GetParameter(0) << " +/- "
            << fGaus->GetParError(0) << '\n';
  std::cout << "  Mean              = " << meanSum << " +/- "
            << fGaus->GetParError(1) << '\n';
  std::cout << "  Sigma             = " << sigmaSum << " +/- "
            << fGaus->GetParError(2) << '\n';
  std::cout << "  P-value           = " << fGaus->GetProb() << '\n';
  std::cout << "  Red chi-squared   = "
            << (fGaus->GetChisquare()) / (fGaus->GetNDF()) << "\n\n";

  // Draw all histograms on canvas
  TCanvas *c1 = new TCanvas("c1", "Pulse analysis", 1300, 700);
  c1->Divide(2, 2);

  c1->cd(1);
  gPad->SetLogz();
  // gPad->SetLeftMargin(0.99f);
  // gPad->SetRightMargin(0.99f);
  // gPad->SetBottomMargin(0.13);
  gPad->Update();
  hAreaVsTime->DrawCopy("COLZ");

  c1->cd(2);
  // hNoise->GetXaxis()->SetRangeUser(2750, 3150.);
  hNoise->SetLineWidth(1);
  hNoise->DrawCopy();

  c1->cd(3);
  gPad->SetLogy();
  gPad->Update();
  hPhotoElectrons->SetLineWidth(1);
  hPhotoElectrons->DrawCopy();

  c1->cd(4);
  hWidth->SetLineWidth(1);
  hWidth->DrawCopy();

  // Draw parameters on canvas
  TCanvas *c3 = new TCanvas("c3", "Parameter plots", 1300, 700);
  c3->Divide(nPulseParam, nPulseParam);
  for (int par_i = 0; par_i < nPulseParam; ++par_i) {
    for (int par_j = 0; par_j < nPulseParam; ++par_j) {
      int canvasIndex{par_i * nPulseParam + par_j + 1};
      c3->cd(canvasIndex);

      if (par_i > par_j) {
        gPulsePar[par_i][par_j]->SetMarkerStyle(20);
        gPulsePar[par_i][par_j]->SetMarkerSize(0.1f);
        gPulsePar[par_i][par_j]->SetMarkerColor(kRed);
        gPulsePar[par_i][par_j]->Draw("AP");
      } else if (par_i < par_j) {
        gPad->SetLogz();
        gPad->Update();
        h2PulsePar[par_i][par_j]->DrawCopy("COLZ");
      } else if (par_i == par_j) {
        hPulsePar[par_i]->DrawCopy("");
        file1->cd();
        hPulsePar[par_i]->Write();
      }
    }
  }

  // Create canvas to display all pulses of one file
  TCanvas *cPulses = new TCanvas("cPulses", "Pulses", 1500, 700);

  // Draw all pulses on multigraph object
  for (size_t i = 0; i < graphs.size(); ++i) {
    mg->Add(graphs[i]);
  }
  cPulses->cd();
  mg->Draw("ALP");
  mg->SetTitle("Pulses");
  mg->SetName("Regions of pulses");
  mg->GetXaxis()->SetTitle("Time since \"trigger\" [ns]");
  mg->GetYaxis()->SetTitle("ADC Counts");

  // Create canvas to superimpose all pulses of one file
  TCanvas *cPulsesSuperimp =
      new TCanvas("cPulsesSuperimp", "Superimposed pulses", 1500, 700);

  // Draw all pulses on multigraph object
  for (size_t i = 0; i < graphsSuperimposed.size(); ++i) {
    mgSuperimposed->Add(graphsSuperimposed[i]);
  }
  cPulsesSuperimp->cd();
  mgSuperimposed->Draw("ALP");
  mgSuperimposed->SetTitle("Superimposed pulses");
  mgSuperimposed->GetXaxis()->SetTitle("Time since startPulse [ns]");
  mgSuperimposed->GetYaxis()->SetTitle("ADC Counts");

  // Create canvas for areas in PE of region of interest
  TCanvas *cPEArea = new TCanvas("cPEArea", "Pulse distribution", 1500, 700);
  cPEArea->Divide(2, 2);

  // Normalise  hPhotoElectrons->Scale(1.0 / hPhotoElectrons->GetMaximum());
  hPETrigger->Scale(1.0 / hPETrigger->GetMaximum());
  hPEPreTrigger->Scale(1.0 / hPEPreTrigger->GetMaximum());
  hPEPostTrigger1->Scale(1.0 / hPEPostTrigger1->GetMaximum());
  hPEPostTrigger2->Scale(1.0 / hPEPostTrigger2->GetMaximum());

  // Draw areas in PE trigger region
  cPEArea->cd(1);
  gPad->SetLogy();
  gPad->Update();
  hPETrigger->SetLineWidth(1);
  hPETrigger->DrawCopy();

  // Draw areas in PE pre trigger region
  cPEArea->cd(2);
  gPad->SetLogy();
  gPad->Update();
  hPEPreTrigger->SetLineWidth(1);
  hPEPreTrigger->DrawCopy();

  // Draw areas in PE post trigger region 1
  cPEArea->cd(3);
  gPad->SetLogy();
  gPad->Update();
  hPEPostTrigger1->SetLineWidth(1);
  hPEPostTrigger1->DrawCopy();

  // Draw areas in PE post trigger region 2
  cPEArea->cd(4);
  gPad->SetLogy();
  gPad->Update();
  hPEPostTrigger2->SetLineWidth(1);
  hPEPostTrigger2->DrawCopy();

  // Save canvases
  c1->SaveAs("./plots/dig/pulse_analysis_results.pdf");
  c3->SaveAs("./plots/dig/params_analysis.pdf");
  cPulses->SaveAs("./plots/dig/pulses.pdf");
  cPulsesSuperimp->SaveAs("./plots/dig/pulsesSuperimposed.pdf");
  cPulseSum->SaveAs("./plots/dig/cPulseSum.pdf");
  cPEArea->SaveAs("./plots/dig/cPEArea.pdf");

  // Write objects on file
  file1->cd();
  c1->Write();
  c3->Write();
  mg->Write();
  cPulses->Write();
  cPulsesSuperimp->Write();
  cPEArea->Write();
  tree->Write();
  // tree->Print();
  file1->Close();
}

// Plot waveform amplitudes as function of time
void waveformTotal() {
  // To avoid reloading manually if .so is present
  R__LOAD_LIBRARY(waveformAnalysisPos_cpp.so);

  // Creating files and canvases
  TFile *file2 = new TFile("./rootFiles/waveform.root", "RECREATE");
  TCanvas *c2 = new TCanvas("c2", "Waveform analysis", 1500, 700);

  const double samplePeriod = 2.0;  // In [ns]
  std::ifstream infile(
      "./data/45Degrees/CH0_3PTFE-LED_45_1.3_2-3.5_70_INC_TRANSM.txt");
  std::string line;

  int row = 0;
  std::vector<TGraph *> graphs;

  // Select random generator seed for colours based on current time
  srand(time(NULL));

  // Loop on rows
  while (std::getline(infile, line)) {
    // Control over analysed rows
    if (row < nMinAnalysedRows) {
      ++row;
      continue;
    }
    if (row >= nMaxAnalysedRows) {
      break;
    }

    // Defining loop variables
    std::stringstream ss(line);
    std::string item;
    int column{1};
    double timestamp{};
    int sampleIndex{};

    std::vector<double> xValues;                          // Relative time
    std::vector<double> yValues;                          // ADC counts
    std::vector<double> colours{1, 3, 4, 5, 6, 7, 8, 9};  // Colour vector

    // Loop on columns
    while (std::getline(ss, item, '\t')) {
      if (item.empty()) {
        continue;
      }

      if (column >= 7) {
        yValues.push_back(std::stod(item));
        xValues.push_back(sampleIndex * samplePeriod);
        ++sampleIndex;
      }
      ++column;
    }

    // Generate a random number between 0 and 7 (used for colour indices)
    int randIndex = rand() % 8;

    // Plot each waveform using a graph object
    TGraph *g = new TGraph(xValues.size(), xValues.data(), yValues.data());
    g->SetLineColor(colours[randIndex]);
    g->SetLineWidth(1);
    g->SetMarkerColor(kBlack);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1);
    g->GetYaxis()->SetRangeUser(2000, 16000.);
    g->SetTitle(Form("Waveform %d; Time [ns]; ADC counts",
                     row + 1));  // Inserting placeholder
    graphs.push_back(g);
    ++row;
  }

  setFitStyle();

  // Draw all graphs
  c2->cd();
  for (size_t i = 0; i < graphs.size(); ++i) {
    if (i == 0) {
      graphs[i]->Draw("ALP");
    } else {
      graphs[i]->Draw("L P SAME");
    }
  }

  // Save canvas
  // c2->BuildLegend(.70, .7, .9, .9, "Legend");
  file2->cd();
  c2->Write();
  file2->Close();

  // Print canvas
  c2->SaveAs("./plots/dig/waveform_plot.png");
}

void rateAnalysis() {
  // 45 DEGREES CONFIGURATION

  // 0.20mm thickness of PTFE
  waveformAnalysis(Transm,
                   "./data/45Degrees/0.20mm/"
                   "DataF_CH0@DT5730S_59483_run_45_0.2_TRANSM_REFL.txt",
                   "./rootFiles/45Degrees0.20mm/wA2.root");

  waveformAnalysis(Refl,
                   "./data/45Degrees/0.20mm/"
                   "DataF_CH1@DT5730S_59483_run_45_0.2_TRANSM_REFL.txt",
                   "./rootFiles/45Degrees0.20mm/wA3.root");

  // 0.80mm thickness of PTFE
  waveformAnalysis(Transm,
                   "./data/45Degrees/0.80mm/"
                   "DataF_CH0@DT5730S_59483_run_45_0.8_TRANSM_REFL.txt",
                   "./rootFiles/45Degrees0.80mm/wA2.root");

  waveformAnalysis(Refl,
                   "./data/45Degrees/0.80mm/"
                   "DataF_CH1@DT5730S_59483_run_45_0.8_TRANSM_REFL.txt",
                   "./rootFiles/45Degrees0.80mm/wA3.root");

  // 1.55mm thickness of PTFE
  waveformAnalysis(Transm,
                   "./data/45Degrees/1.55mm/"
                   "DataF_CH0@DT5730S_59483_run_45_1.55_TRANSM_REFL.txt",
                   "./rootFiles/45Degrees1.55mm/wA2.root");

  waveformAnalysis(Refl,
                   "./data/45Degrees/1.55mm/"
                   "DataF_CH1@DT5730S_59483_run_45_1.55_TRANSM_REFL.txt",
                   "./rootFiles/45Degrees1.55mm/wA3.root");

  // 2.05mm thickness of PTFE
  waveformAnalysis(Transm,
                   "./data/45Degrees/2.05mm/"
                   "DataF_CH0@DT5730S_59483_run_45_2.05_TRANSM_REFL.txt",
                   "./rootFiles/45Degrees2.05mm/wA2.root");

  waveformAnalysis(Refl,
                   "./data/45Degrees/2.05mm/"
                   "DataF_CH1@DT5730S_59483_run_45_2.05_TRANSM_REFL.txt",
                   "./rootFiles/45Degrees2.05mm/wA3.root");
  // 3.10mm thickness of PTFE
  waveformAnalysis(Transm,
                   "./data/45Degrees/3.10mm/"
                   "DataF_CH0@DT5730S_59483_run_45_3.10_TRANSM_REFL.txt",
                   "./rootFiles/" + a + "Degrees5.15mm/"};

    // Loop over folders
    for (auto &f : folders) {
      std::string incTransmFile = f + "wA0.root";
      std::string incReflFile = f + "wA1.root";

      // Copy incidence transmittance
      if (gSystem->CopyFile(incTransm.c_str(), incTransmFile.c_str())) {
        std::cout << "Copied transmittance to " << incTransmFile << '\n';
      } else {
        std::cout << "Failed to copy Transm to " << incTransmFile << '\n';
      }

      // Copy incidence reflectance
      if (gSystem->CopyFile(incRefl.c_str(), incReflFile.c_str())) {
        std::cout << "Copied Refl to " << incReflFile << '\n';
      } else {
        std::cout << "Failed to copy Refl to " << incReflFile << '\n';
      }
    }
  }
}

// Plot working principle of pulse finder
void plotWaveform(
    AreaConvFactor areaConv = Transm,
    std::string infileName =
        "./data/miscellaneous/DataF_CH0@DT5730S_59483_run_new_1300_2-3.5.txt",
    std::string rootFileName = "./rootFiles/miscellaneous/plotWF.root") {
  // To avoid reloading manually if .so is present
  R__LOAD_LIBRARY(waveformAnalysisPos_cpp.so);

  // Area conversion factor (current assumption is 1 PE = 11000 ADC*ns)
  // For Transm PMT = 4000.
  // For Refl PMT = 12500.
  auto const areaConvFactor = static_cast<double>(areaConv);

  // Variables used later
  double const samplePeriod = 2.0;  // In [ns]
  std::ifstream infile(infileName.c_str());
  std::string line;
  std::vector<double> colours{1, 3, 4, 5, 6, 7, 8, 9};  // Colour vector
  std::vector<TGraph *> graphs{};
  std::vector<TF1 *> thresholds{};
  std::vector<TBox *> boxes{};
  std::vector<TLine *> lines{};
  TMultiGraph *mg = new TMultiGraph();
  int row{0};

  // Creating TFile
  TFile *file2 = new TFile(rootFileName.c_str(), "RECREATE");

  // Create canvas to store waveforms
  TCanvas *c2 = new TCanvas("c2", "Plot waveform", 1500, 700);

  // Select random generator seed for colours based on current time
  srand(time(NULL));

  // Loop over rows (waveforms) to extract pulse sum graph
  while (std::getline(infile, line)) {
    // Control over analysed rows
    if (row < nMinAnalysedRows) {
      ++row;
      continue;
    }
    if (row >= nMaxAnalysedRows) {
      break;
    }

    // Defining loop variables
    std::stringstream ss(line);
    std::string item;
    std::vector<double> samples;
    std::vector<double> xValues;
    double timestamp = 0.;
    int sampleIndex{};
    int column = 1;

    setFitStyle();

    // Loop over columns
    while (std::getline(ss, item, '\t')) {
      if (item.empty()) continue;
      if (column == 3) {
        timestamp = std::stod(item);
      } else if (column >= 7) {
        samples.push_back(std::stod(item));
        xValues.push_back(sampleIndex * samplePeriod);
        ++sampleIndex;
      }
      ++column;
    }

    // Creating WaveformAnalysis object
    WaveformAnalysisPos wf(samples, timestamp, samplePeriod);

    // Generate a random number between 0 and 7 (used for colour indices)
    int randIndex = rand() % 8;

    // Print waveform properties
    std::cout << std::fixed
              << std::setprecision(2);  // Round to 2 decimal place
    std::cout << "\n********** Waveform n. " << row + 1 << " **********\n";
    std::cout << "Timestamp        = " << wf.getTimeStamp() << " ns\n";
    std::cout << "Baseline         = " << wf.getBaseline() << " ADC counts\n";
    std::cout << "Sample period    = " << wf.getSamplePeriod() << " ns\n";

    // Get pulse vector from each single waveform
    const auto &pulses = wf.getPulses();
    std::cout << "Number of Pulses without selection = " << pulses.size()
              << '\n';

  // Draw threshold for given waveform
    TF1 *fThreshold = new TF1(Form("fThreshold_%d", row), "[0]", 0., 450.);
    fThreshold->SetLineWidth(4);
    fThreshold->SetLineStyle(2);
    fThreshold->SetLineColor(colours[randIndex]);
    fThreshold->SetParNames("Const");
    fThreshold->FixParameter(0, wf.getThreshold());  // Const
    thresholds.push_back(fThreshold);

    // Plot each waveform using a graph object
    TGraph *g = new TGraph(xValues.size(), xValues.data(), samples.data());
    g->SetLineColor(colours[randIndex]);
    g->SetLineWidth(1);
    g->SetMarkerColor(kBlack);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1);
    graphs.push_back(g);

    // Loop on pulses
    for (size_t i = 0; i < pulses.size(); ++i) {
      const auto &p = pulses[i];

      // Insert selections on pulses
      if (p.peakValue > 15900.) {
        continue;
      }

      // Define box coordinates for pulse
      double x1 = p.startTime - wf.getTimeStamp();
      double x2 = p.endTime - wf.getTimeStamp();
      double y1 = p.values[0];
      double y2 = p.peakValue;

      TBox *box = new TBox(x1, y1, x2, y2);
      box->SetFillColor(kRed - 9);
      box->SetFillStyle(3354);  // Hatched style
      box->SetLineColor(kRed);
      boxes.push_back(box);

      // Draw boundary lines
      TLine *l1 = new TLine(x1, y1, x1, y2);
      TLine *l2 = new TLine(x2, y1, x2, y2);
      l1->SetLineColor(kRed);
      l2->SetLineColor(kRed);
      lines.push_back(l1);
      lines.push_back(l2);

      // Params of interest
      double heightOverWidth{p.peakValue / (p.endTime - p.startTime)};
      double peakFractionPos{(p.peakTime - p.startTime) /
                             (p.endTime - p.startTime)};
      double areaOverFullTime{p.area / ((p.endTime - p.startTime))};

      std::cout << "\n  *** Pulse n. " << i + 1 << " ***\n\n";
      std::cout << "  Overall start time           = " << p.startTime
                << " ns\n";
      std::cout << "  Overall end time             = " << p.endTime << " ns\n";
      std::cout << "  Overall peak time            = " << p.peakTime << " ns\n";
      std::cout << "  Relative start time          = "
                << p.startTime - wf.getTimeStamp() << " ns\n";
      std::cout << "  Relative end time            = "
                << p.endTime - wf.getTimeStamp() << " ns\n";
      std::cout << "  Relative peak time           = "
                << p.peakTime - wf.getTimeStamp() << " ns\n";
      std::cout << "  Peak time since startPulse   = "
                << p.peakTime - wf.getTimeStamp() - p.times[0] << " ns\n";
      std::cout << "  Peak value                   = " << p.peakValue
                << " ADC\n";
      std::cout << "  Width                        = "
                << p.endTime - p.startTime << " ns\n";
      std::cout << "  Rise time                    = " << p.riseTime << " ns\n";
      std::cout << "  FWHM                         = " << p.FWHMTime << " ns\n";
      std::cout << "  90% area time                = " << p.areaFractionTime
                << " ns\n";
      std::cout << "  Height over width            = " << heightOverWidth
                << " ADC/ns\n";
      std::cout << "  Peak fraction pos.           = " << peakFractionPos
                << '\n';
      std::cout << "  Area / full width            = " << areaOverFullTime
                << " ADC\n";
      std::cout << "  Area                         = " << p.area << " ADC*ns\n";
      std::cout << "  Area in PE                   = "
                << p.area / areaConvFactor << " PE\n";
      std::cout << "  Negative/overall area frac   = " << p.negFracArea
                << " \n";
      std::cout << "  Negative/overall counts      = " << p.negFrac << " \n";
    }
    ++row;
  }

  // Draw all graphs
  for (size_t i = 0; i < graphs.size(); ++i) {
    mg->Add(graphs[i]);
  }
  mg->SetTitle("Pulse finder");
  mg->SetName("mg");
  mg->GetXaxis()->SetTitle("Time since \"trigger\" [ns]");
  mg->GetYaxis()->SetTitle("ADC Counts");
  c2->cd();
  mg->Draw("ALP");
  for (auto &t : thresholds) {
    t->Draw("SAME");
}
for (auto &l : lines) {
    l->Draw("SAME");
  }
  for (auto &b : boxes) {
    b->Draw("SAME");
  }

  // Save canvas
  c2->Update();
  file2->cd();
  c2->Write();
  file2->Close();
}

int main() {
  rateAnalysis();

  return EXIT_SUCCESS;
}