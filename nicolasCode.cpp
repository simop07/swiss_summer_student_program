#include <iostream>
#include <string>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGraph.h"
#include "TGraph2D.h"
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

using namespace std;

struct Parameter {
  string label;
  double min;
  double max;
};

int main() {
  // Ranges
  double min_TBA = -1;
  double max_TBA = 0;
  double min_S1area = 1.e4;
  double max_S1area = 10.e4;
  double min_RR = 0;
  double max_RR = pow(147.6 / 2., 2);
  double min_XY = -75.;
  double max_XY = 75.;
  double min_dt = 0;
  double max_dt = 1.4e3;
  double PI = 3.1415927;
  double min_phi = -PI;
  double max_phi = PI;
  double min_height = 0;
  double max_height = 100.;
  double min_ls2area = 4;
  double max_ls2area = 8;
  double min_chiSq = 0;
  double max_chiSq = 18;
  int n_bins = 50;

  const int n_deriv_par = 22;
  Parameter deriv_par[n_deriv_par];
  deriv_par[0] = {"s1 area [phd]", min_S1area, max_S1area};
  deriv_par[1] = {"s1 height [phd]", min_height, max_height};
  deriv_par[2] = {"s1 TBA", min_TBA, max_TBA};
  deriv_par[3] = {"s1 T centr X [cm]", min_XY, max_XY};
  deriv_par[4] = {"s1 T centr Y [cm]", min_XY, max_XY};
  deriv_par[5] = {"s1 B centr X [cm]", min_XY, max_XY};
  deriv_par[6] = {"s1 B centr Y [cm]", min_XY, max_XY};
  deriv_par[7] = {"s1 T+B centr X [cm]", 2 * min_XY, 2 * max_XY};
  deriv_par[8] = {"s1 T+B centr Y [cm]", 2 * min_XY, 2 * max_XY};
  deriv_par[9] = {"R of T+B centr [cm^2]", min_RR, 10 * max_XY};
  deriv_par[10] = {"R of T+ Ro B centr [cm^2]", min_RR, 10 * max_XY};
  deriv_par[11] = {"Cent T [cm^2]", min_RR, max_RR};
  deriv_par[12] = {"Cent B [cm^2]", min_RR, max_RR};
  deriv_par[13] = {"Cent T-B [cm^2]", min_RR, max_RR};
  deriv_par[14] = {"drift [us]", min_dt, max_dt};
  deriv_par[15] = {"s2 merc x [cm]", min_XY, max_XY};
  deriv_par[16] = {"s2 merc y [cm]", min_XY, max_XY};
  deriv_par[17] = {"RR [cm^2]", min_RR, max_RR};
  deriv_par[18] = {"phi centr [rad]", min_phi, max_phi};
  deriv_par[19] = {"phi [rad]", min_phi, max_phi};
  deriv_par[20] = {"phi diff [rad]", min_phi, max_phi};
  deriv_par[21] = {"Area fraction in max PMT", 0, 1};

  // initialise plots
  TH1F *h_deriv_par[n_deriv_par];
  TGraph *g2_deriv_par[n_deriv_par][n_deriv_par];
  TGraph2D *g_reco_RR = new TGraph2D();
  TGraph *g_reco_dt = new TGraph();
  TH2F *h2_deriv_par[n_deriv_par][n_deriv_par];
  for (int par_ii = 0; par_ii < n_deriv_par; par_ii++) {
    h_deriv_par[par_ii] =
        new TH1F("", "", n_bins, deriv_par[par_ii].min, deriv_par[par_ii].max);
    for (int par_jj = 0; par_jj < n_deriv_par; par_jj++) {
      if (par_ii < par_jj) g2_deriv_par[par_ii][par_jj] = new TGraph();
      if (par_ii > par_jj)
        h2_deriv_par[par_ii][par_jj] = new TH2F(
            "", "", n_bins, deriv_par[par_ii].min, deriv_par[par_ii].max,
            n_bins, deriv_par[par_jj].min, deriv_par[par_jj].max);
    }
  }

  // for (int evt_ii = 0; evt_ii < ceil(total_evt_entries * 4. / 10.); evt_ii++)
  // {
  double deriv_par_values[n_deriv_par] = {};
  //       s1_area_phd->at(0),
  //       s1_height_phd->at(0),
  //       s1_TBA->at(0),
  //       s1_TcentrX_cm->at(0),
  //       s1_TcentrY_cm->at(0),
  //       s1_BcentrX_cm->at(0),
  //       s1_BcentrY_cm->at(0),
  //       s1_TcentrX_cm->at(0) + s1_BcentrX_cm->at(0),
  //       s1_TcentrY_cm->at(0) + s1_BcentrY_cm->at(0),
  //       pow(s1_TcentrY_cm->at(0) + s1_BcentrY_cm->at(0), 2) +
  //           pow(s1_TcentrX_cm->at(0) + s1_BcentrX_cm->at(0), 2),
  //       CentT + CentB,
  //       CentT,
  //       CentB,
  //       CentT - CentB,
  //       drift_T_us,
  //       s2_merc_x->at(0),
  //       s2_merc_y->at(0),
  //       RR,
  //       phi_s1,
  //       phi_s2,
  //       phi_s2 - phi_s1,
  //       1
  //       // s1_maxChArea_phd->at(0)/s1_area_phd->at(0)
  //   };

  for (int par_ii = 0; par_ii < n_deriv_par; par_ii++) {
    h_deriv_par[par_ii]->Fill(deriv_par_values[par_ii]);
    for (int par_jj = 0; par_jj < n_deriv_par; par_jj++) {
      if (par_ii < par_jj)
        g2_deriv_par[par_ii][par_jj]->AddPoint(deriv_par_values[par_ii],
                                               deriv_par_values[par_jj]);
      if (par_ii > par_jj)
        h2_deriv_par[par_ii][par_jj]->Fill(deriv_par_values[par_ii],
                                           deriv_par_values[par_jj]);
    }
  }
}