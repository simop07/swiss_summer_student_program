// To compile in SHELL: "analysis.cpp `root-config --cflags --libs`"
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
#include "TMatrixD.h"
#include "TMultiGraph.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "waveformAnalysisPos.hpp"

// Define global constants
constexpr int nMinAnalysedRows{0};  // Minimum index of analysed rows EXCLUDED
constexpr int nMaxAnalysedRows{9961};  // Maximum rows INCLUDED (9961)

// Asymmetric gaussian function
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

void fitPEHisto(TH1F *hPhotoElectrons) {
  // Import user defined function asymmetric gaussians for 1 PE
  TF1 *fAsymmetric1PE = new TF1("fAsymmetric1PE", asymGaussians, 0.5, 1.8, 4);
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
  TF1 *fAsymmetric2PE = new TF1("fAsymmetric2PE", asymGaussians, 1.8, 2.6, 4);
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
  TF1 *fExpo = new TF1("fExpo", "expo", 0., 4.2);
  fExpo->SetLineColor(kGreen + 2);
  fExpo->SetLineWidth(4);
  fExpo->SetLineStyle(2);
  fExpo->SetParameter(0, -0.4);  // Constant
  fExpo->SetParameter(1, -1.2);  // Slope

  // Define total function as sum of 1 PE + 2 PE
  TF1 *fTotal = new TF1("fTotal", asym2GaussiansExpo, 0., 4.2, 11);
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
      new TF1("fTotalConst", asym2GaussiansExpoConstrained, 0., 4.2, 10);
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
              0., 4.2, 8);
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
  gStyle->SetOptFit(1111);
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
void waveformAnalysis() {
  // To avoid reloading manually if .so is present
  R__LOAD_LIBRARY(waveformAnalysisPos_cpp.so);

  // Area conversion factor (current assumption is 1 PE = 11000 ADC*ns)
  double const areaConvFactor{11000.};

  // Variables used later
  double const samplePeriod = 2.0;  // In [ns]
  std::ifstream infile(
      "./data/DataR_CH0@DT5730S_59483_250321_led_on_no_cover_3_2.txt");
  std::string line;
  std::vector<double> colours{1, 3, 4, 5, 6, 7, 8, 9};  // Colour vector
  TMultiGraph *mg = new TMultiGraph();
  TMultiGraph *mgSuperimposed = new TMultiGraph();
  std::vector<TGraph *> graphs{};
  std::vector<TGraph *> graphsSuperimposed{};
  std::map<int, double> map{};
  int row = 0;

  // Select random generator seed for colours based on current time
  srand(time(NULL));

  // Creating TFile
  TFile *file1 = new TFile("./rootFiles/waveformAnalysis.root", "RECREATE");

  // Define histograms
  TH2F *hAreaVsTime = new TH2F("hAreaVsTime",
                               "Pulse area vs relative time; Relative time "
                               "peak [ns]; Area [ADC #times ns]",
                               40, 100., 300., 100, 0., 30000.);
  TH1F *hNoise = new TH1F("hNoise", "Noise distribution; ADC counts; Counts",
                          30, 8020, 8050);
  TH1F *hPhotoElectrons =
      new TH1F("hPE", "Pulse area distribution; Area [PE]; Normalized counts",
               150, 0, 6);
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

  // Variable for analysis in trigger region in 1 single file
  int pulseCounter{};
  int pulseCounterTriggerRegion{};
  double totTrigArea{};
  double numTrigPE{};

  // Below gaussian fit on "Pulse sum" graph is used (2 sigmas)
  double const triggerStart{179.12236 - 2 * 15.27892};
  double const triggerEnd{179.12236 + 2 * 15.27892};

  // Variable for analysis in pre-trigger region in 1 single file
  int pulseCounterPreTriggerRegion{};
  double totPreTrigArea{};
  double numPreTrigPE{};
  double const preTriggerStart{100.};
  double const preTriggerEnd{130.};

  // Variable for analysis in post-trigger region 1 in 1 single file
  int pulseCounterPostTriggerRegion1{};
  double totPostTrigArea1{};
  double numPostTrigPE1{};
  double const postTriggerStart1{219.};
  double const postTriggerEnd1{239.};

  // Variable for analysis in post-trigger region 2 in 1 single file
  int pulseCounterPostTriggerRegion2{};
  double totPostTrigArea2{};
  double numPostTrigPE2{};
  double const postTriggerStart2{264.};
  double const postTriggerEnd2{284.};

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
      if (item.empty()) continue;
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
      // if (((p.area / areaConvFactor) < (1.01477 - 1 * 0.59664)) /* ||
      //     ((p.endTime - p.startTime) < 40.) */) {
      //   continue;
      // }

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

      std::cout << "  *** Pulse n. " << i + 1 << " ***\n";
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

      // Generate sum of pulses using a STL map. Here the keys of the map
      // correspond to the time values of pulses, while values correpond to the
      // voltage/ADC counts. Contrary to std::vector<T>, std::map<T1,T2>'s
      // operator[] is more secure: when you have an empty map, with no
      // specified size, and you access map[X], the map safely creates the
      // key-value pair (X,0.0) - in this case I put 0.0 because it is the
      // default initialiser for double values
      for (int j = 0; j < p.times.size(); ++j) {
        map[p.times[j]] += p.values[j];
      }

      // Create vector to superimpose pulses
      std::vector<double> superimposedTimes = p.times;
      double shift{superimposedTimes[0]};
      for (int timeId{}; timeId < superimposedTimes.size(); ++timeId) {
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

  // Define rates per region
  double rateTrig = numTrigPE / (triggerEnd - triggerStart);
  double ratePreTrig = numPreTrigPE / (preTriggerEnd - preTriggerStart);
  double ratePostTrig1 = numPostTrigPE1 / (postTriggerEnd1 - postTriggerStart1);
  double ratePostTrig2 = numPostTrigPE2 / (postTriggerEnd2 - postTriggerStart2);

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

  // Fit PE histograms with several functions
  fitPEHisto(hPhotoElectrons);

  // Draw all histograms on canvas
  TCanvas *c1 = new TCanvas("c1", "Pulse analysis", 1300, 700);
  c1->Divide(2, 2);

  setFitStyle();

  c1->cd(1);
  gPad->SetLogz();
  // gPad->SetLeftMargin(0.99f);
  // gPad->SetRightMargin(0.99f);
  // gPad->SetBottomMargin(0.13);
  gPad->Update();
  hAreaVsTime->DrawCopy("COLZ");

  c1->cd(2);
  // hNoise->GetXaxis()->SetRangeUser(8010, 8050.);
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

  // Create fit function for summed pulses
  TF1 *fGaus = new TF1("fGaus", "gaus", 120., 220.);
  fGaus->SetLineColor(kRed);
  fGaus->SetLineWidth(4);
  fGaus->SetLineStyle(2);
  fGaus->SetParameter(0, 2e7);   // Amplitude
  fGaus->SetParameter(1, 180.);  // Mean
  fGaus->SetParameter(2, 20.);   // Sigma
  gPulseSum->Fit(fGaus, "M R");

  // Draw summed pulses
  cPulseSum->cd();
  gPulseSum->Draw("ALP");

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
  cPulses->Write();
  cPulsesSuperimp->Write();
  cPEArea->Write();
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
      "./data/DataR_CH0@DT5730S_59483_250321_led_on_no_cover_3_2.txt");
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
  c2->BuildLegend(.70, .7, .9, .9, "Legend");
  file2->cd();
  c2->Write();
  file2->Close();

  // Print canvas
  c2->SaveAs("./plots/dig/waveform_plot.png");
}

int main() {
  waveformAnalysis();
  waveformTotal();

  EXIT_SUCCESS;
}

// Unused code

// Import user defined function symmetric gaussian for 1 PE
// TF1 *fSymmetric1PE = new TF1("fSymmetric1PE", "gaus", 0.3, 1.8);
// fSymmetric1PE->SetLineColor(kGreen);
// fSymmetric1PE->SetLineWidth(4);
// fSymmetric1PE->SetLineStyle(2);
// fSymmetric1PE->SetParNames("N^{1}", "#mu^{1}", "#sigma^{1}");
// fSymmetric1PE->SetParameter(0, 0.3);  // Constant
// fSymmetric1PE->SetParameter(1, 1.);   // #mu
// fSymmetric1PE->SetParameter(2, 0.5);  // #sigma

// Import user defined function symmetric gaussian for 2 PE
// TF1 *fSymmetric2PE = new TF1("fSymmetric2PE", "gaus", 1.8, 2.8);
// fSymmetric2PE->SetLineColor(kYellow);
// fSymmetric2PE->SetLineWidth(4);
// fSymmetric2PE->SetLineStyle(2);
// fSymmetric2PE->SetParNames("N^{2}", "#mu^{2}", "#sigma^{2}");
// fSymmetric2PE->SetParameter(0, 0.1);  // Constant
// fSymmetric2PE->SetParameter(1, 2.2);  // #mu
// fSymmetric2PE->SetParameter(2, 0.5);  // #sigma

// Import user defined function symmetric gaussian for 3 PE
// TF1 *fSymmetric3PE = new TF1("fSymmetric3PE", "gaus", 2.8, 3.55);
// fSymmetric3PE->SetLineColor(kRed);
// fSymmetric3PE->SetLineWidth(4);
// fSymmetric3PE->SetLineStyle(2);
// fSymmetric3PE->SetParNames("N^{3}", "#mu^{3}", "#sigma^{3}");
// fSymmetric3PE->SetParameter(0, 0.01);  // Constant
// fSymmetric3PE->SetParameter(1, 3.1);   // #mu
// fSymmetric3PE->SetParameter(2, 2.0);   // #sigma

// Import user defined function asymmetric gaussians for 3 PE
// TF1 *fAsymmetric3PE = new TF1("fAsymmetric3PE", asymGaussians, 2.8, 3.55,
// 4); fAsymmetric3PE->SetLineColor(kViolet); fAsymmetric3PE->SetLineWidth(4);
// fAsymmetric3PE->SetLineStyle(2);
// fAsymmetric3PE->SetParNames("N^{3}", "#mu^{3}", "#sigma^{3}_{1}",
//                             "#sigma^{3}_{2}");
// fAsymmetric3PE->SetParameter(0, 0.01);  // Constant
// fAsymmetric3PE->SetParameter(1, 3.1);   // #mu
// fAsymmetric3PE->SetParameter(2, 2.0);   // #sigma_{1}
// fAsymmetric3PE->SetParameter(3, 2.0);   // #sigma_{2}