#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include "TStyle.h"

void setFitStyle() {
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);  // It was 10
  gStyle->SetOptFit(0);   // It was 1111
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
  gStyle->SetTitleYOffset(1.0f);
  gStyle->SetLineScalePS(1);
  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleYSize(0.05);
  // gStyle->SetPadTopMargin(-9.);
  // gStyle->SetPadRightMargin(-9.);
  // gStyle->SetPadBottomMargin(-9.);
  // gStyle->SetPadLeftMargin(-9.);
  // gStyle->SetTitleW(0.5f);
}

void compareDig() {
  // Loading ROOT File
  TFile *file1 = new TFile("./rootFiles/miscellaneous/wfAnalysis.root", "READ");
  TFile *file2 = new TFile("./rootFiles/miscellaneous/wfAnalysis.root", "READ");

  // Creating ROOT File
  TFile *fileCompare =
      new TFile("./rootFiles/miscellaneous/compareDig.root", "RECREATE");

  // Reading histos from canvas
  int const nPulseParam{14};
  TH1F *hPulsePar1[nPulseParam];
  TH1F *hPulsePar2[nPulseParam];
  for (int par_i = 0; par_i < nPulseParam; ++par_i) {
    hPulsePar1[par_i] = (TH1F *)file1->Get(Form("h1PulsePar_%d", par_i));
    hPulsePar2[par_i] = (TH1F *)file2->Get(Form("h1PulsePar_%d", par_i));
  }

  setFitStyle();

  // Plotting these instograms histograms in ROOT File
  TCanvas *c4 = new TCanvas("c4", "Compare params", 1300, 700);
  c4->Divide(6, 3);
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

  setFitStyle();

  // Plotting these instograms histograms in ROOT File
  TCanvas *c4 = new TCanvas("c4", "Compare params", 1300, 700);
  c4->Divide(6, 3);
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

void compareDigOsc() {
  // To avoid reloading manually if .so is present

  // Loading ROOT File
  TFile *file1 = new TFile("./rootFiles/waveformAnalysis.root", "READ");
  TFile *file2 = new TFile("./rootFiles/waveformAnalysisOsc.root", "READ");

  // Creating ROOT File
  TFile *fileCompare = new TFile("./rootFiles/compareDigOsc.root", "RECREATE");

  // Reading histos from canvas
  int const nPulseParam{14};
  TH1F *hPulsePar1[nPulseParam];
  TH1F *hPulsePar2[nPulseParam];
  for (int par_i = 0; par_i < nPulseParam; ++par_i) {
    hPulsePar1[par_i] = (TH1F *)file1->Get(Form("h1PulsePar_%d", par_i));
    hPulsePar2[par_i] = (TH1F *)file2->Get(Form("h1PulsePar_%d", par_i));
  }

  setFitStyle();

  // Plotting these instograms histograms in ROOT File
  TCanvas *c[nPulseParam];
  for (int i = 0; i < nPulseParam; ++i) {
    c[i] = new TCanvas(Form("c%d", i + 1), "Compare params", 1300, 700);
    c[i]->Divide(2);
    c[i]->cd(1);
    hPulsePar1[i]->DrawCopy("");
    c[i]->cd(2);
    hPulsePar2[i]->SetLineColor(kRed);
    hPulsePar2[i]->DrawCopy("");
    fileCompare->cd();
    c[i]->Write();
  }

  // Superimpose graphs
  TCanvas *cSuperimp =
      new TCanvas("cSuperimp", "Superimpose params", 1300, 700);
  cSuperimp->Divide(6, 3);
  for (int par_i = 0; par_i < nPulseParam; ++par_i) {
    cSuperimp->cd(par_i + 1);
    hPulsePar2[par_i]->SetLineColor(kRed);
    hPulsePar2[par_i]->DrawCopy("");
    hPulsePar1[par_i]->DrawCopy("SAME");
    gPad->Update();
  }

  fileCompare->cd();
  cSuperimp->Write();
  fileCompare->Close();
  file1->Close();
}