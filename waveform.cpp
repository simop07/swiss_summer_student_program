// To compile in SHELL: "analysis.cpp `root-config --cflags --libs`"
#include <time.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
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
#include "waveformAnalysis.hpp"

// Define global constants
constexpr int nMinAnalysedRows{0};  // minimum index of analysed rows excluded
constexpr int nMaxAnalysedRows{9961};  // maximum rows included (9961)

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

// Asymmetric gaussian function
Double_t asym2GaussiansExpo(Double_t *x, Double_t *par) {
  return asymGaussians(x, &par[0]) + asymGaussians(x, &par[4]) +
         TMath::Exp(par[8] + par[9] * x[0]) + par[10];
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
  fAsymmetric1PE->SetParameter(0, 0.3);  // Constant
  fAsymmetric1PE->SetParameter(1, 1.);   // #mu
  fAsymmetric1PE->SetParameter(2, 0.5);  // #sigma_{1}
  fAsymmetric1PE->SetParameter(3, 0.5);  // #sigma_{2}

  // Import user defined function asymmetric gaussians for 2 PE
  TF1 *fAsymmetric2PE = new TF1("fAsymmetric2PE", asymGaussians, 1.8, 2.6, 4);
  fAsymmetric2PE->SetLineColor(kOrange + 2);
  fAsymmetric2PE->SetLineWidth(4);
  fAsymmetric2PE->SetLineStyle(2);
  fAsymmetric2PE->SetParNames("N^{2}", "#mu^{2}", "#sigma^{2}_{1}",
                              "#sigma^{2}_{2}");
  fAsymmetric2PE->SetParameter(0, 0.1);  // Constant
  fAsymmetric2PE->SetParameter(1, 2.2);  // #mu
  fAsymmetric2PE->SetParameter(2, 1);    // #sigma_{1}
  fAsymmetric2PE->SetParameter(3, 0.5);  // #sigma_{2}

  // Define expo function for noise fit
  TF1 *fExpo = new TF1("fExpo", "expo", 0., 4.2);
  fExpo->SetLineColor(kGreen + 2);
  fExpo->SetLineWidth(4);
  fExpo->SetLineStyle(2);
  fExpo->SetParameter(0, -0.4);  // Constant
  fExpo->SetParameter(1, -1.2);  // Slope

  // Define total function as sum of 1 PE + 2 PE
  TF1 *fTotal = new TF1("fTotal", asym2GaussiansExpo, 0., 4.2, 11);
  fTotal->SetLineColor(kBlack);
  fTotal->SetLineWidth(4);
  fTotal->SetLineStyle(2);

  // Define parameter array for total function
  double par[11];

  // Normalise  hPhotoElectrons->Scale(1.0 / hPhotoElectrons->GetMaximum());
  hPhotoElectrons->Scale(1.0 / hPhotoElectrons->GetMaximum());

  // Fit PE graphs

  // 1PE
  hPhotoElectrons->Fit(fAsymmetric1PE, "R");
  fAsymmetric1PE->GetParameters(&par[0]);
  // Print pvalue and reduced chi squared
  std::cout << "\n\n**** FIT RESULT 1 PE peak ****\n\nP value       "
               "      = "
            << fAsymmetric1PE->GetProb() << "\n";
  std::cout << "Reduced chi squared = "
            << fAsymmetric1PE->GetChisquare() / fAsymmetric1PE->GetNDF()
            << "\n\n";

  // 2PE
  hPhotoElectrons->Fit(fAsymmetric2PE, "R+");
  fAsymmetric2PE->GetParameters(&par[4]);
  // Print pvalue and reduced chi squared
  std::cout << "\n\n**** FIT RESULT 2 PE peak ****\n\nP value       "
               "      = "
            << fAsymmetric2PE->GetProb() << "\n";
  std::cout << "Reduced chi squared = "
            << fAsymmetric2PE->GetChisquare() / fAsymmetric2PE->GetNDF()
            << "\n\n";

  // Exponential noise
  hPhotoElectrons->Fit(fExpo, "+", "", 0., 0.2);
  fExpo->GetParameters(&par[8]);
  // Print pvalue and reduced chi squared
  std::cout << "\n\n**** FIT RESULT EXPO NOISE ****\n\nP value       "
               "      = "
            << fExpo->GetProb() << "\n";
  std::cout << "Reduced chi squared = "
            << fExpo->GetChisquare() / fExpo->GetNDF() << "\n\n";

  // Total function fit
  par[10] = 0.001;
  fTotal->SetParameters(par);
  fTotal->SetParNames("N^{1}", "#mu^{1}", "#sigma^{1}_{1}", "#sigma^{1}_{2}",
                      "N^{2}", "#mu^{2}", "#sigma^{2}_{1}", "#sigma^{2}_{2}",
                      "Constant", "Slope");
  fTotal->SetParName(10, "Background");
  TFitResultPtr fitResult = hPhotoElectrons->Fit(fTotal, "S R+");

  // Get results
  std::cout << "\n\n**** FIT RESULT TOTAL INDPENDENT #mu ****\n\nP value       "
               "      = "
            << fTotal->GetProb() << "\n";
  std::cout << "Reduced chi squared = "
            << fTotal->GetChisquare() / fTotal->GetNDF() << "\n\n";
  TMatrixD covMatrix = fitResult->GetCorrelationMatrix();
  TMatrixD corMatrix = fitResult->GetCovarianceMatrix();
  std::cout << "\n*** Print covariance matrix ***\n" << std::endl;
  covMatrix.Print();
  std::cout << "\n*** Print correlation matrix ***\n" << std::endl;
  corMatrix.Print();
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
  // gStyle->SetPadTopMargin(-9.);
  // gStyle->SetPadRightMargin(-9.);
  // gStyle->SetPadBottomMargin(-9.);
  // gStyle->SetPadLeftMargin(-9.);
  // gStyle->SetTitleW(0.5f);
}

// This function analyses the waveform by building areaVStime, Noise, PE counts
// and pulseWidth histos
void waveformAnalysis() {
  // To avoid reloading manually if .so is present
  R__LOAD_LIBRARY(waveformAnalysis_cpp.so);

  double const samplePeriod = 2.0;  // In [ns]
  std::ifstream infile(
      "DataR_CH0@DT5730S_59483_250321_led_on_no_cover_3_2.txt");
  std::string line;
  int row = 0;

  // Creating TFile
  TFile *file1 = new TFile("waveformAnalysis.root", "RECREATE");

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
      new TH1F("hWidth", "Width distribution; Width [ns]; Counts", 40, 2, 50);

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
    WaveformAnalysis wf(samples, timestamp, samplePeriod);

    // Print waveform properties
    std::cout << std::fixed
              << std::setprecision(1);  // Round to 1 decimal place
    std::cout << "\n********** Waveform n. " << row + 1 << " **********\n";
    std::cout << "Timestamp        = " << wf.getTimeStamp() << " ns\n";
    std::cout << "Baseline         = " << wf.getBaseline() << " ADC counts\n";
    std::cout << "Sample period    = " << wf.getSamplePeriod() << " ns\n";

    // Get pulse vector from each single waveform
    const auto &pulses = wf.getPulses();
    std::cout << "Number of Pulses = " << pulses.size() << "\n";

    // Print pulse properties
    for (size_t i = 0; i < pulses.size(); ++i) {
      const auto &p = pulses[i];
      std::cout << std::fixed
                << std::setprecision(1);  // Round to 1 decimal place
      std::cout << "  *** Pulse n. " << i + 1 << " ***\n";
      std::cout << "  Overall start time  = " << p.startTime << " ns\n";
      std::cout << "  Overall end time    = " << p.endTime << " ns\n";
      std::cout << "  Overall peak time   = " << p.peakTime << " ns\n";
      std::cout << "  Relative start time = " << p.startTime - wf.getTimeStamp()
                << " ns\n";
      std::cout << "  Relative end time   = " << p.endTime - wf.getTimeStamp()
                << " ns\n";
      std::cout << "  Relative peak time  = " << p.peakTime - wf.getTimeStamp()
                << " ns\n";
      std::cout << "  Peak value          = " << p.peakValue << " ADC\n";
      std::cout << "  Width               = " << p.endTime - p.startTime
                << " ns\n";
      std::cout << "  Area                = " << p.area << " ADC*ns\n";
      std::cout << "  Area in PE          = " << p.area / 11000. << " PE\n";

      // Fill pulse information
      hAreaVsTime->Fill(p.peakTime - wf.getTimeStamp(), p.area);
      hWidth->Fill(p.endTime - p.startTime);

      // Convert area into PE (current assumption is 1 PE = 11000 ADC*ns)
      double areaInPE = p.area / 11000.;
      hPhotoElectrons->Fill(areaInPE);
    }

    // Fill noise information
    int const nBaselineSamples{50};
    std::for_each(samples.begin(), samples.begin() + nBaselineSamples,
                  [&](double sample) { hNoise->Fill(sample); });
    ++row;
  }

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

  c1->SaveAs("pulse_analysis_results.pdf");
  file1->cd();
  c1->Write();
  file1->Close();
}

// Plot waveform amplitudes as function of time
void waveformTotal() {
  // To avoid reloading manually if .so is present
  R__LOAD_LIBRARY(waveformAnalysis_cpp.so);

  // Creating files and canvases
  TFile *file2 = new TFile("waveform.root", "RECREATE");
  TCanvas *c2 = new TCanvas("c2", "Waveform analysis", 1500, 700);

  const double samplePeriod = 2.0;  // In [ns]
  std::ifstream infile(
      "DataR_CH0@DT5730S_59483_250321_led_on_no_cover_3_2.txt");
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

    std::vector<double> xValues;                             // Relative time
    std::vector<double> yValues;                             // ADC counts
    std::vector<double> colours{1, 2, 3, 4, 5, 6, 7, 8, 9};  // Colour vector

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

    // Generate a random number between 0 and 8 (used for colour indices)
    int randIndex = rand() % 9;

    // Plot each waveform using a graph object
    TGraph *g = new TGraph(xValues.size(), xValues.data(), yValues.data());
    g->SetLineColor(colours[randIndex]);
    g->SetLineWidth(1);
    g->SetMarkerColor(kBlack);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1);
    g->GetXaxis()->SetRangeUser(4., 300.);
    g->GetYaxis()->SetRangeUser(7500., 16000.);
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
  c2->SaveAs("waveform_plot.png");
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
// TF1 *fAsymmetric3PE = new TF1("fAsymmetric3PE", asymGaussians, 2.8, 3.55, 4);
// fAsymmetric3PE->SetLineColor(kViolet);
// fAsymmetric3PE->SetLineWidth(4);
// fAsymmetric3PE->SetLineStyle(2);
// fAsymmetric3PE->SetParNames("N^{3}", "#mu^{3}", "#sigma^{3}_{1}",
//                             "#sigma^{3}_{2}");
// fAsymmetric3PE->SetParameter(0, 0.01);  // Constant
// fAsymmetric3PE->SetParameter(1, 3.1);   // #mu
// fAsymmetric3PE->SetParameter(2, 2.0);   // #sigma_{1}
// fAsymmetric3PE->SetParameter(3, 2.0);   // #sigma_{2}