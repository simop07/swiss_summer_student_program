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
#include "waveformAnalysisNeg.hpp"

// Define global constants
constexpr int nMinAnalysedRows{0};  // Minimum index of analysed rows EXCLUDED
constexpr int nMaxAnalysedRows{1021};  // Maximum rows INCLUDED (1021)

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
void fitPEHisto(TH1F *hPhotoElectrons) {
  // Import user defined function asymmetric gaussians for 1 PE
  TF1 *fAsymmetric1PE = new TF1("fAsymmetric1PE", asymGaussians, 0.5, 2., 4);
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
  TF1 *fAsymmetric2PE = new TF1("fAsymmetric2PE", asymGaussians, 2., 2.6, 4);
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
  TF1 *fTotal = new TF1("fTotal", asym2Gaussians, 0., 4.2, 9);
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
            << fAsymmetric1PE->GetProb() << "\n";
  std::cout << "Reduced chi squared = "
            << fAsymmetric1PE->GetChisquare() / fAsymmetric1PE->GetNDF()
            << "\n\n";

  // 2PE
  hPhotoElectrons->Fit(fAsymmetric2PE, "RN");
  fAsymmetric2PE->GetParameters(&par1[4]);
  // Print pvalue and reduced chi squared
  std::cout << "\n\n**** FIT RESULT 2 PE peak ****\n\nP value       "
               "      = "
            << fAsymmetric2PE->GetProb() << "\n";
  std::cout << "Reduced chi squared = "
            << fAsymmetric2PE->GetChisquare() / fAsymmetric2PE->GetNDF()
            << "\n\n";

  // Total function fit NOT CONSTRAINED

  par1[8] = 0.005;
  fTotal->SetParameters(par1);
  fTotal->SetParNames("N^{1}", "#mu^{1}", "#sigma^{1}_{1}", "#sigma^{1}_{2}",
                      "N^{2}", "#mu^{2}", "#sigma^{2}_{1}", "#sigma^{2}_{2}",
                      "Background");
  TFitResultPtr fitResult = hPhotoElectrons->Fit(fTotal, "S R+");

  // Get results
  std::cout
      << "\n\n**** FIT RESULT TOTAL NOT CONSTRAINED #mu ****\n\nP value       "
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

  // Total function fit CONSTRAINED

  // Define total function as sum of 1 PE + 2 PE
  TF1 *fTotalConst =
      new TF1("fTotalConst", asym2GaussiansConstrained, 0., 4.2, 8);
  fTotalConst->SetLineColor(kOrange + 2);
  fTotalConst->SetLineWidth(4);
  fTotalConst->SetLineStyle(2);

  fTotalConst->SetParameters(par1[0], par1[1], par1[2], par1[3], par1[4],
                             par1[6], par1[7], par1[8]);
  fTotalConst->SetParNames("N^{1}", "#mu^{1}", "#sigma^{1}_{1}",
                           "#sigma^{1}_{2}", "N^{2}", "#sigma^{2}_{1}",
                           "#sigma^{2}_{2}", "Background");
  TFitResultPtr fitResultConst = hPhotoElectrons->Fit(fTotalConst, "S R+");

  // Get results
  std::cout
      << "\n\n**** FIT RESULT TOTAL CONSTRAINED #mu ****\n\nP value       "
         "      = "
      << fTotalConst->GetProb() << "\n";
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
      "fTotalConstSameSigmas", asym2GaussiansConstrainedSameSigma, 0., 4.2, 6);
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
            << fTotalConstSameSigmas->GetProb() << "\n";
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

// This function analyses the waveform by building areaVStime, Noise, PE counts
// and pulseWidth histos
void waveformAnalysis() {
  // To avoid reloading manually if .so is present
  R__LOAD_LIBRARY(waveformAnalysisNeg_cpp.so);

  double const samplePeriod = 0.005e-1;  // In [\mus]
  std::ifstream infile("4Layers.txt");
  std::string line;
  int row = 0;

  // Creating TFile
  TFile *file1 = new TFile("waveformAnalysis.root", "RECREATE");

  // Define histograms
  TH2F *hAreaVsTime = new TH2F("hAreaVsTime",
                               "Pulse area vs relative time; Relative time "
                               "peak [#mus]; Area [mV #times ns]",
                               60, 0.4, 2.1, 60, 0., 100.);
  TH1F *hNoise = new TH1F("hNoise", "Noise distribution; Voltage [mV]; Counts",
                          17, -0.6, 0.6);
  TH1F *hPhotoElectrons =
      new TH1F("hPE", "Pulse area distribution; Area [PE]; Normalized counts",
               80, 0, 3.5);
  TH1F *hWidth =
      new TH1F("hWidth", "Width distribution; Width [ns]; Counts", 40, 1., 30.);

  // Plotting parameters
  int const nPulseParam{5};
  Parameter pulsePar[nPulseParam] = {{"Width [ns]", 1., 30.},
                                     {"Area [mV #times ns]", 0., 100.},
                                     {"Area [PE]", 0, 3.5},
                                     {"Relative peak time [#mus]", 0.4, 2.1},
                                     {"Noise RMS [mV]", -0.6, 0.6}};
  int const nBins{200};
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
        timestamp = std::stod(item) * 1e6;
      else if (column >= 7)
        samples.push_back(std::stod(item) * 1e3);
      ++column;
    }

    // Creating WaveformAnalysis object
    WaveformAnalysisNeg wf(samples, timestamp, samplePeriod);

    // Print waveform properties
    std::cout << std::fixed
              << std::setprecision(10);  // Round to 10 decimal place
    std::cout << "\n********** Waveform n. " << row + 1 << " **********\n";
    std::cout << "Timestamp        = " << wf.getTimeStamp() << " #mus\n";
    std::cout << "Baseline         = " << wf.getBaseline() << " mV\n";
    std::cout << "Sample period    = " << wf.getSamplePeriod() << " #mus\n";

    // Get pulse vector from each single waveform
    const auto &pulses = wf.getPulses();
    std::cout << "Number of Pulses = " << pulses.size() << "\n";

    // Fill noise information
    int const nBaselineSamples{1100};
    std::for_each(samples.begin(), samples.begin() + nBaselineSamples,
                  [&](double sample) { hNoise->Fill(sample); });

    // Print pulse properties
    for (size_t i = 0; i < pulses.size(); ++i) {
      const auto &p = pulses[i];
      std::cout << "  *** Pulse n. " << i + 1 << " ***\n";
      std::cout << "  Overall start time  = " << p.startTime << " #mus\n";
      std::cout << "  Overall end time    = " << p.endTime << " #mus\n";
      std::cout << "  Overall peak time   = " << p.peakTime << " #mus\n";
      std::cout << "  Relative start time = " << p.startTime - wf.getTimeStamp()
                << " #mus\n";
      std::cout << "  Relative end time   = " << p.endTime - wf.getTimeStamp()
                << " #mus\n";
      std::cout << "  Relative peak time  = " << p.peakTime - wf.getTimeStamp()
                << " #mus\n";
      std::cout << "  Peak value          = " << p.peakValue << " mV\n";
      std::cout << "  Width               = "
                << (p.endTime - p.startTime) * 1000 << " ns\n";
      std::cout << "  Area                = " << p.area * 1000 << " mV*ns\n";
      std::cout << "  Area in PE          = " << p.area * 1000 / 24. << " PE\n";

      // Fill pulse information
      hAreaVsTime->Fill(p.peakTime - wf.getTimeStamp(), p.area * 1000);
      hWidth->Fill((p.endTime - p.startTime) * 1000);  // In [ns]

      // Convert area into PE (current assumption is 1 PE = 11000 ADC*ns)
      double areaInPE = p.area * 1000 / 24.;
      hPhotoElectrons->Fill(areaInPE);

      // Plot all parameters against each other

      // Gather params of interest from each pulse and noise from waveform
      double parValues[nPulseParam] = {
          (p.endTime - p.startTime) * 1000, p.area * 1000, p.area * 1000 / 24.,
          p.peakTime - wf.getTimeStamp(), wf.getBaseline()};
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

  // Round fit printing to 10 decimal place
  std::cout << std::fixed << std::setprecision(10);

  // FITTING FUNCTIONS SPACE
  std::cout << "\n\n\n************************************\n";
  std::cout << "****** FIITING FUNCTION SPACE ******\n";
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
        gPulsePar[par_i][par_j]->SetMarkerSize(0.2f);
        gPulsePar[par_i][par_j]->SetMarkerColor(kRed);
        gPulsePar[par_i][par_j]->Draw("AP");
      } else if (par_i < par_j) {
        gPad->SetLogz();
        gPad->Update();
        h2PulsePar[par_i][par_j]->DrawCopy("COLZ");
      } else if (par_i == par_j) {
        hPulsePar[par_i]->DrawCopy("");
      }
    }
  }

  c1->SaveAs("pulse_analysis_results.pdf");
  c3->SaveAs("params_analysis.pdf");

  file1->cd();
  c1->Write();
  c3->Write();
  file1->Close();
}

// Plot waveform amplitudes as function of time
void waveformTotal() {
  // To avoid reloading manually if .so is present
  R__LOAD_LIBRARY(waveformAnalysisNeg_cpp.so);

  // Creating files and canvases
  TFile *file2 = new TFile("waveform.root", "RECREATE");
  TCanvas *c2 = new TCanvas("c2", "Waveform analysis", 1500, 700);

  const double samplePeriod = 0.005e-1;  // In [\mus]
  std::ifstream infile("4Layers.txt");
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
        yValues.push_back(std::stod(item) * 1e3);
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
    // g->GetXaxis()->SetRangeUser(-0.5e-6, 2.1e-6);
    // g->GetYaxis()->SetRangeUser(7500., 16000.);
    g->SetTitle(Form("Waveform %d; Time [#mus]; Voltage [mV]",
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