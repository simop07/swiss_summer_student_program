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
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "waveformAnalysis.hpp"

// Define global constants
constexpr int nMinAnalysedRows{0};     // minimum index of analysed rows (0)
constexpr int nMaxAnalysedRows{9961};  // maximum rows (9961)

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

  double const samplePeriod_ns = 2.0;
  std::ifstream infile(
      "DataR_CH0@DT5730S_59483_250321_led_on_no_cover_3_2.txt");
  std::string line;
  int row = 0;

  // Creating TFile
  TFile *file = new TFile("waveformAnalysis.root", "RECREATE");

  // Histograms (mentioned above)
  TH2F *hAreaVsTime = new TH2F("hAreaVsTime",
                               "Pulse area vs time since start; Time since "
                               "start [ns]; Area [ADC #times ns]",
                               33, 100., 300., 100, 0., 30000.);
  TH1F *hNoise = new TH1F("hNoise", "Noise distribution; ADC counts; Counts",
                          30, 8020, 8050);
  TH1F *hPhotoElectrons =
      new TH1F("hPE", "Pulse area distribution; Area [PE]; Normalized counts",
               1000, 0, 6);
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
    WaveformAnalysis wf(samples, timestamp, samplePeriod_ns);

    // Print waveform properties
    std::cout << std::fixed << std::setprecision(1);  // Use more precision
    std::cout << "\n========== Waveform #" << row + 1 << " ==========\n";
    std::cout << "Timestamp       = " << wf.getTimestamp() << " ns\n";
    std::cout << "Baseline        = " << wf.getBaseline() << " ADC counts\n";
    std::cout << "Sample Period   = " << wf.getSamplePeriod() << " ns\n";

    // Get pulse vector from each single waveform
    const auto &pulses = wf.getPulses();
    std::cout << "Number of Pulses= " << pulses.size() << "\n";

    // Print pulse properties
    for (size_t i = 0; i < pulses.size(); ++i) {
      const auto &p = pulses[i];
      std::cout << "  --- Pulse #" << i + 1 << " ---\n";
      std::cout << "  Overall Start time    = " << p.startTime << " ns\n";
      std::cout << "  Overall End time      = " << p.endTime << " ns\n";
      std::cout << "  Overall peak time     = " << p.peakTime << " ns\n";
      std::cout << "  Start time since start= "
                << p.startTime - wf.getTimestamp() << " ns\n";
      std::cout << "  End time since start  = " << p.endTime - wf.getTimestamp()
                << " ns\n";
      std::cout << "  Peak time since start = "
                << p.peakTime - wf.getTimestamp() << " ns\n";
      std::cout << "  Peak value            = " << p.peakValue << " ADC\n";
      std::cout << "  Width                 = " << p.endTime - p.startTime
                << " ns\n";
      std::cout << "  Area                  = " << p.area << " ADC*ns\n";
      std::cout << "  Area in PE            = " << p.area / 10400. << " PE\n";

      // Fill pulse info
      hAreaVsTime->Fill(p.peakTime - wf.getTimestamp(), p.area);
      hWidth->Fill(p.endTime - p.startTime);

      // Convert area to PE (if you know 1 PE ~ X ADC*ns)
      double areaInPE = p.area / 10400.;  // Example: 1 PE = 10400 ADC*ns
      hPhotoElectrons->Fill(areaInPE);
    }

    // Fill noise info
    int const nBaselineSamples{50};
    std::for_each(samples.begin(), samples.begin() + nBaselineSamples,
                  [&](double sample) { hNoise->Fill(sample); });
    ++row;
  }

  // Fit noise with gaussian function
  hNoise->Fit("gaus");

  // Draw all histograms on canvas
  TCanvas *c1 = new TCanvas("c1", "Pulse Analysis", 1300, 700);
  c1->Divide(2, 2);

  c1->cd(1);
  gPad->SetLogz();
  // gPad->SetLeftMargin(0.99f);
  // gPad->SetRightMargin(0.99f);
  // gPad->SetBottomMargin(0.13);
  gPad->Update();
  hAreaVsTime->DrawCopy("COLZ");

  c1->cd(2);
  // hNoise->GetXaxis()->SetRangeUser(8010, 8050.);
  hNoise->DrawCopy();

  c1->cd(3);
  gPad->SetLogy();
  gPad->Update();
  hPhotoElectrons->Scale(1.0 / hPhotoElectrons->GetMaximum());  // Normalise
  hPhotoElectrons->DrawCopy();

  c1->cd(4);
  hWidth->DrawCopy();

  c1->SaveAs("pulse_analysis_results.pdf");
  file->cd();
  c1->Write();
  file->Close();
}

// Plot waveform amplitudes as function of time
void waveformTotal() {
  // To avoid reloading manually if .so is present
  R__LOAD_LIBRARY(waveformAnalysis_cpp.so);

  // Creating files and canvases
  TFile *file = new TFile("waveform.root", "RECREATE");
  TCanvas *c2 = new TCanvas("c2", "Waveform", 1500, 700);

  const double samplePeriod_ns = 2.0;
  std::ifstream infile(
      "DataR_CH0@DT5730S_59483_250321_led_on_no_cover_3_2.txt");
  std::string line;

  int row = 0;
  std::vector<TGraph *> graphs;

  // Select random generator for colours
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
    int column = 1;
    double timestamp = 0.;
    int sample_index = 0;

    std::vector<double> x_vals;                              // Time data
    std::vector<double> y_vals;                              // Amplitude data
    std::vector<double> colours{1, 2, 3, 4, 5, 6, 7, 8, 9};  // Colour vector

    // Loop on columns
    while (std::getline(ss, item, '\t')) {
      if (item.empty()) continue;

      if (column == 3) {
        timestamp = std::stod(item);
      } else if (column >= 7) {
        double sample_value = std::stod(item);
        double time_ns = sample_index * samplePeriod_ns;
        x_vals.push_back(time_ns);
        y_vals.push_back(sample_value);
        ++sample_index;
      }
      ++column;
    }

    int RandIndex = rand() % 9;  // generates a random number between 0 and 8

    // Plot each waveform with a graph object
    TGraph *g = new TGraph(x_vals.size(), x_vals.data(), y_vals.data());
    g->SetLineColor(colours[RandIndex]);
    g->SetLineWidth(2);
    g->SetMarkerColor(kBlack);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1);
    g->GetXaxis()->SetRangeUser(4., 300.);
    g->GetYaxis()->SetRangeUser(7500., 16000.);
    g->SetTitle(Form("Waveform %d; Time [ns]; ADC Counts",
                     row + 1));  // Inserting placeholder
    graphs.push_back(g);
    ++row;
  }

  setFitStyle();

  // Draw all graphs
  c2->cd();
  for (size_t i = 0; i < graphs.size(); i++) {
    if (i == 0) {
      graphs[i]->Draw("ALP");
    } else {
      graphs[i]->Draw("L P SAME");
    }
  }

  // Save canvas
  c2->BuildLegend(.70, .7, .9, .9, "Legend");
  file->cd();
  c2->Write();
  file->Close();

  // Print canvas
  c2->SaveAs("waveform_plot.png");
}

int main() {
  waveformAnalysis();
  waveformTotal();

  EXIT_SUCCESS;
}