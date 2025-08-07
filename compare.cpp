#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "waveformAnalysisPos.hpp"

void compareDig() {
  // To avoid reloading manually if .so is present
  R__LOAD_LIBRARY(waveformAnalysisPos_cpp.so);

  // Loading ROOT File
  TFile *file1 = new TFile("./rootFiles/waveformAnalysis.root", "READ");
  TFile *file2 = new TFile("./rootFiles/waveformAnalysis2.root", "READ");

  // Creating ROOT File
  TFile *fileCompare = new TFile("./rootFiles/compareDig.root", "RECREATE");

  // Reading histos from canvas
  int const nPulseParam{14};
  TH1F *hPulsePar1[nPulseParam];
  TH1F *hPulsePar2[nPulseParam];
  for (int par_i = 0; par_i < nPulseParam; ++par_i) {
    hPulsePar1[par_i] = (TH1F *)file1->Get(Form("h1PulsePar_%d", par_i));
    hPulsePar2[par_i] = (TH1F *)file2->Get(Form("h1PulsePar_%d", par_i));
  }

  // Plotting these instograms histograms in ROOT File
  TCanvas *c4 = new TCanvas("c4", "Compare params", 1300, 700);
  c4->Divide(7, 2);
  for (int par_i = 0; par_i < nPulseParam; ++par_i) {
    c4->cd(par_i + 1);
    hPulsePar1[par_i]->DrawCopy("");
    hPulsePar2[par_i]->SetLineColor(kRed);
    hPulsePar2[par_i]->DrawCopy("SAME");
    gPad->Update();
  }
  // c4->BuildLegend(.70, .7, .9, .9, "Legend");

  c4->SaveAs("./plots/params_comparison.pdf");

  fileCompare->cd();
  c4->Write();
  fileCompare->Close();
  file1->Close();
}

void compareOsc() {
  // To avoid reloading manually if .so is present
  R__LOAD_LIBRARY(waveformAnalysisNeg_cpp.so);

  // Loading ROOT File
  TFile *file1 = new TFile("./rootFiles/waveformAnalysisOsc.root", "READ");
  TFile *file2 = new TFile("./rootFiles/waveformAnalysisOsc2.root", "READ");

  // Creating ROOT File
  TFile *fileCompare = new TFile("./rootFiles/compareOsc.root", "RECREATE");

  // Reading histos from canvas
  int const nPulseParam{14};
  TH1F *hPulsePar1[nPulseParam];
  TH1F *hPulsePar2[nPulseParam];
  for (int par_i = 0; par_i < nPulseParam; ++par_i) {
    hPulsePar1[par_i] = (TH1F *)file1->Get(Form("h1PulsePar_%d", par_i));
    hPulsePar2[par_i] = (TH1F *)file2->Get(Form("h1PulsePar_%d", par_i));
  }

  // Plotting these instograms histograms in ROOT File
  TCanvas *c4 = new TCanvas("c4", "Compare params", 1300, 700);
  c4->Divide(7, 2);
  for (int par_i = 0; par_i < nPulseParam; ++par_i) {
    c4->cd(par_i + 1);
    hPulsePar2[par_i]->SetLineColor(kRed);
    hPulsePar2[par_i]->DrawCopy("");
    hPulsePar1[par_i]->DrawCopy("SAME");
    gPad->Update();
  }
  // c4->BuildLegend(.70, .7, .9, .9, "Legend");

  c4->SaveAs("./plots/params_comparisonOsc.pdf");

  fileCompare->cd();
  c4->Write();
  fileCompare->Close();
  file1->Close();
}