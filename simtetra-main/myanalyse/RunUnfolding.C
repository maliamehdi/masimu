// RunUnfolding.C
//
// Driver pour lancer TROIS unfoldings en parallèle :
//   1) Méthode linéaire itérative (IterativeUnfoldLinear)
//   2) Méthode bayésienne itérative (IterativeUnfoldBayes)
//   3) Méthode directe (type "spectral data", DirectSpectralUnfold)
//
// + Comparaison :
//   - g_mesuré (tronqué / "histo fille") vs g_replié_lin vs g_replié_bayes vs g_replié_direct
//   - f_unfold_lin vs f_unfold_bayes vs f_unfold_direct
//   - χ²/point vs itération (lin vs bayes)
//
// Prérequis :
//   .L IterativeUnfoldLinear.C+
//   .L IterativeUnfoldBayes.C+
//   .L DirectSpectralUnfold.C+
//   .L RunUnfolding.C+
//
// Usage typique :
//   RunUnfolding("Response_PARIS50.root", "hResp",
//                "resolution0508_bckgndsub_promptsumEvents_PARIS_corrected.root",
//                "sum_Res11keVGamma_PARIS50",
//                "UnfoldCompare_PARIS50.root");
//
// ----------------------------------------------------------------------

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TSystem.h>
#include <TAxis.h>
#include <TArrayD.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <cmath>

// ----------------------------------------------------------------------
// Déclarations des fonctions d'unfolding (implémentées ailleurs)
// ----------------------------------------------------------------------

// Méthode linéaire itérative
TH1D* IterativeUnfoldLinear(const TH1*  hMeas,
                            const TH2D* hResp,
                            int   maxIter,
                            bool  normalizeResponse,
                            bool  enforcePositivity,
                            double relChi2Tol,
                            int   minIter,
                            TH1D** hRefoldOut,
                            TGraph** gChi2PerPointOut);

// Méthode bayésienne itérative
TH1D* IterativeUnfoldBayes(const TH1*  hMeas,
                           const TH2D* hResp,
                           int   maxIter,
                           double tolRelChi2,
                           bool  enforcePos,
                           bool  verbose,
                           TGraph** gChi2,
                           TGraph** gChi2NDF);

// Méthode directe (spectral data) – implémentée dans DirectSpectralUnfold.C
TH1D* DirectSpectralUnfold(const TH1*  hMeas,
                           const TH2D* hResp,
                           const char* nameOut,
                           bool        enforcePos,
                           bool        verbose,
                           TH1D**      hResidualOut,
                           TH1D**      hRefoldOut);

// ----------------------------------------------------------------------
// Helper : reconstruire le spectre replié g’(Emeas) = Σ_j R_ij f_j
// à partir de hUnfold (axe Y de hResp) et de la matrice de réponse hResp.
// On normalise colonne par colonne : PDF(Emeas | Etrue_j).
// (Utilisé pour la méthode bayésienne, quand on n’a pas déjà le g’ en sortie.)
// ----------------------------------------------------------------------
static TH1D* BuildRefoldedFromUnfold(const TH1D* hUnfold,
                                     const TH2D* hResp,
                                     const char* name = "hRefold_fromUnfold")
{
  if (!hUnfold || !hResp) return nullptr;

  const TH2D* hR2D = hResp;
  int nbinsX = hR2D->GetNbinsX(); // Emeas
  int nbinsY = hR2D->GetNbinsY(); // Etrue

  if (hUnfold->GetNbinsX() != nbinsY) {
    std::cerr << "[BuildRefoldedFromUnfold] WARNING: hUnfold nbins ("
              << hUnfold->GetNbinsX() << ") != hResp Y nbins ("
              << nbinsY << "). Results might be inconsistent." << std::endl;
  }

  // Construire TH1D sur l'axe X (Emeas) de la réponse
  const TAxis* axX = hR2D->GetXaxis();
  TH1D* hRef = nullptr;
  const TArrayD* xb = axX->GetXbins();
  if (xb && xb->GetSize() == nbinsX+1) {
    hRef = new TH1D(name, "Refolded spectrum;E_{meas};Counts",
                    nbinsX, xb->GetArray());
  } else {
    hRef = new TH1D(name, "Refolded spectrum;E_{meas};Counts",
                    nbinsX, axX->GetXmin(), axX->GetXmax());
  }
  hRef->SetDirectory(nullptr);

  // Construire R_ij normalisé colonne par colonne (PDF(Emeas | Etrue_j))
  std::vector< std::vector<double> > R(nbinsX, std::vector<double>(nbinsY, 0.0));
  std::vector<double> colSum(nbinsY, 0.0);

  for (int ix = 1; ix <= nbinsX; ++ix) {
    for (int iy = 1; iy <= nbinsY; ++iy) {
      double v = hR2D->GetBinContent(ix, iy);
      if (v < 0) v = 0.0;
      R[ix-1][iy-1] = v;
      colSum[iy-1] += v;
    }
  }

  for (int j = 0; j < nbinsY; ++j) {
    double s = colSum[j];
    if (s > 0) {
      for (int i = 0; i < nbinsX; ++i) R[i][j] /= s;
    } else {
      for (int i = 0; i < nbinsX; ++i) R[i][j] = 0.0;
    }
  }

  // Lire f_j
  std::vector<double> f(nbinsY, 0.0);
  for (int j = 0; j < nbinsY; ++j) {
    f[j] = hUnfold->GetBinContent(j+1);
    if (f[j] < 0) f[j] = 0.0;
  }

  // g'_i = Σ_j R_ij f_j
  for (int i = 0; i < nbinsX; ++i) {
    double ti = 0.0;
    for (int j = 0; j < nbinsY; ++j) {
      ti += R[i][j] * f[j];
    }
    hRef->SetBinContent(i+1, ti);
  }

  return hRef;
}

// ----------------------------------------------------------------------
// Helper : construire un "histo fille" hMeas_used qui contient
// les N premiers bins de hMeasOrig, avec le binning X de la matrice hResp.
// N = nbins X de hResp.
// ----------------------------------------------------------------------
static TH1D* MakeMeasuredChildFirstNBins(const TH1*  hMeasOrig,
                                         const TH2D* hResp,
                                         const char* name = "hMeas_used")
{
  if (!hMeasOrig || !hResp) return nullptr;

  const TAxis* axResp = hResp->GetXaxis();
  int    Nresp  = axResp->GetNbins();
  const TArrayD* edges = axResp->GetXbins();

  TH1D* h = nullptr;

  // Cas binning variable pour la réponse
  if (edges && edges->GetSize() == Nresp+1) {
    h = new TH1D(name,
                 "Measured (first N bins, response X-binning);E_{meas};Counts",
                 Nresp, edges->GetArray());
  } else {
    // Cas binning uniforme
    h = new TH1D(name,
                 "Measured (first N bins, response X-binning);E_{meas};Counts",
                 Nresp, axResp->GetXmin(), axResp->GetXmax());
  }
  h->SetDirectory(nullptr);

  int Nmother = hMeasOrig->GetNbinsX();
  int Ncopy   = std::min(Nresp, Nmother);

  std::cout << "[INFO] Building measured child histo: Nresp = " << Nresp
            << ", Nmother = " << Nmother
            << ", copying first Ncopy = " << Ncopy << " bins." << std::endl;

  for (int ib = 1; ib <= Ncopy; ++ib) {
    h->SetBinContent(ib, hMeasOrig->GetBinContent(ib));
    h->SetBinError  (ib, hMeasOrig->GetBinError  (ib));
  }
  // Les bins au-delà de Ncopy (si Nresp > Nmother) restent à 0.

  return h;
}

// ----------------------------------------------------------------------
// Driver principal
// ----------------------------------------------------------------------
void RunUnfolding(const char* responseFile      = "Response_PARIS50.root",
                  const char* hRespName         = "hResp",
                  const char* dataFile          = "resolution0508_bckgndsub_promptsumEvents_PARIS_corrected.root",
                  const char* hMeasName         = "sum_Res11keVGamma_PARIS50",
                  const char* outFile           = "UnfoldCompare_PARIS50.root",
                  // paramètres linéaire
                  int   maxIterLin              = 100,
                  int   minIterLin              = 3,
                  double relChi2TolLin          = 1e-3,
                  bool  normalizeResponseLin    = true,
                  bool  enforcePosLin           = true,
                  // paramètres bayésien
                  int   maxIterBay              = 100,
                  double tolRelChi2Bay          = 1e-4,
                  bool  enforcePosBay           = true,
                  bool  verboseBay              = true)
{
  // ------------------ 1) Ouvrir les fichiers d'entrée ------------------
  TFile* fResp = TFile::Open(responseFile, "READ");
  if (!fResp || fResp->IsZombie()) {
    std::cerr << "[ERROR] Cannot open response file: " << responseFile << std::endl;
    return;
  }

  TH2D* hResp = dynamic_cast<TH2D*>(fResp->Get(hRespName));
  if (!hResp) {
    std::cerr << "[ERROR] TH2D '" << hRespName
              << "' not found in " << responseFile << std::endl;
    fResp->Close();
    return;
  }

  TFile* fData = TFile::Open(dataFile, "READ");
  if (!fData || fData->IsZombie()) {
    std::cerr << "[ERROR] Cannot open data file: " << dataFile << std::endl;
    fResp->Close();
    return;
  }

  TH1* hMeasOrig = dynamic_cast<TH1*>(fData->Get(hMeasName));
  if (!hMeasOrig) {
    std::cerr << "[ERROR] TH1 '" << hMeasName
              << "' not found in " << dataFile << std::endl;
    fResp->Close();
    fData->Close();
    return;
  }

  std::cout << "[INFO] Loaded response '" << hRespName << "' from " << responseFile << std::endl;
  std::cout << "[INFO] Loaded measured (original) '" << hMeasName << "' from " << dataFile << std::endl;

  // Histo complet original (pour sortir)
  TH1* hMeasFull = (TH1*)hMeasOrig->Clone("hMeas_full");
  hMeasFull->SetDirectory(nullptr);

  // ------------------ 2) Construire l’histo "fille" hMeas ------------------
  TH1D* hMeas = MakeMeasuredChildFirstNBins(hMeasOrig, hResp, "hMeas_used");
  if (!hMeas) {
    std::cerr << "[ERROR] Failed to build child measured histogram." << std::endl;
    fResp->Close();
    fData->Close();
    return;
  }

  std::cout << "[INFO] Final hMeas : nbins = " << hMeas->GetNbinsX()
            << ", Xmin = " << hMeas->GetXaxis()->GetXmin()
            << ", Xmax = " << hMeas->GetXaxis()->GetXmax()
            << std::endl;

  // ------------------ 3) Unfolding linéaire ------------------
  TH1D*  hRefoldLin = nullptr;
  TGraph* gChi2Lin  = nullptr;

  std::cout << "[INFO] Running linear iterative unfolding..." << std::endl;
  TH1D* hUnfoldLin = IterativeUnfoldLinear(hMeas, hResp,
                                           maxIterLin,
                                           normalizeResponseLin,
                                           enforcePosLin,
                                           relChi2TolLin,
                                           minIterLin,
                                           &hRefoldLin,
                                           &gChi2Lin);

  if (!hUnfoldLin) {
    std::cerr << "[ERROR] Linear unfolding failed (hUnfoldLin is null)." << std::endl;
  }

  // ------------------ 4) Unfolding bayésien ------------------
  TGraph* gChi2Bay    = nullptr;
  TGraph* gChi2NDFBay = nullptr;

  std::cout << "[INFO] Running Bayesian iterative unfolding..." << std::endl;
  TH1D* hUnfoldBay = IterativeUnfoldBayes(hMeas, hResp,
                                          maxIterBay,
                                          tolRelChi2Bay,
                                          enforcePosBay,
                                          verboseBay,
                                          &gChi2Bay,
                                          &gChi2NDFBay);

  if (!hUnfoldBay) {
    std::cerr << "[ERROR] Bayesian unfolding failed (hUnfoldBay is null)." << std::endl;
  }

  // Repliage bayésien (R * f_bayes, colonnes normalisées)
  TH1D* hRefoldBay = nullptr;
  if (hUnfoldBay) {
    hRefoldBay = BuildRefoldedFromUnfold(hUnfoldBay, hResp, "hRefold_bayes");
  }

  // Version "par fission" pour Bayes (scaling a posteriori)
  TH1D* hUnfoldBayScaled = nullptr;
  if (hUnfoldBay) {
    hUnfoldBayScaled = (TH1D*)hUnfoldBay->Clone("Scaled_hUnfold_Bayes_perFission");
    hUnfoldBayScaled->SetDirectory(nullptr);
    // À adapter si tu changes le nombre de fissions
    hUnfoldBayScaled->Scale(1.0 / 2.882315e+09);
  }

  // ------------------ 5) Unfolding direct (spectral data) ------------------
  TH1D* hResDirect    = nullptr;
  TH1D* hRefoldDirect = nullptr;

  std::cout << "[INFO] Running DIRECT spectral unfolding..." << std::endl;
  TH1D* hUnfoldDirect = DirectSpectralUnfold(hMeas, hResp,
                                             "hUnfold_Direct",
                                             true,   // enforcePos
                                             true,   // verbose
                                             &hResDirect,
                                             &hRefoldDirect);

  if (!hUnfoldDirect) {
    std::cerr << "[ERROR] Direct unfolding failed (hUnfoldDirect is null)." << std::endl;
  }

  // ------------------ 6) Dessins de comparaison ------------------

  // 6.1 Mesuré (fille) vs repliés
  TCanvas* c1 = new TCanvas("cUnfold_MeasVsRefold",
                            "Measured(child) vs Refolded (Linear / Bayesian / Direct)",
                            900, 700);
  c1->SetGrid();

  TH1* hMeasClone = (TH1*)hMeas->Clone("hMeas_used_clone");
  hMeasClone->SetDirectory(nullptr);
  hMeasClone->SetLineColor(kBlack);
  hMeasClone->SetLineWidth(2);

  hMeasClone->SetTitle("Measured (child N bins) vs Refolded;E_{meas};Counts");
  hMeasClone->Draw("HIST");

  if (hRefoldLin) {
    hRefoldLin->SetLineColor(kRed+1);
    hRefoldLin->SetLineWidth(2);
    hRefoldLin->Draw("HIST SAME");
  }
  if (hRefoldBay) {
    hRefoldBay->SetLineColor(kBlue+1);
    hRefoldBay->SetLineWidth(2);
    hRefoldBay->SetLineStyle(2);
    hRefoldBay->Draw("HIST SAME");
  }
  if (hRefoldDirect) {
    hRefoldDirect->SetLineColor(kGreen+2);
    hRefoldDirect->SetLineWidth(2);
    hRefoldDirect->SetLineStyle(7);
    hRefoldDirect->Draw("HIST SAME");
  }

  TLegend* leg1 = new TLegend(0.55,0.68,0.88,0.88);
  leg1->AddEntry(hMeasClone,    "Measured (child, N first bins)", "l");
  if (hRefoldLin)    leg1->AddEntry(hRefoldLin,    "Refolded (Linear)",  "l");
  if (hRefoldBay)    leg1->AddEntry(hRefoldBay,    "Refolded (Bayesian)","l");
  if (hRefoldDirect) leg1->AddEntry(hRefoldDirect, "Refolded (Direct)",  "l");
  leg1->Draw();

  // 6.2 Spectres unfolded linéaire vs bayésien vs direct
  TCanvas* c2 = new TCanvas("cUnfold_Spectra",
                            "Unfolded spectra (Linear / Bayesian / Direct)",
                            900, 700);
  c2->SetGrid();

  bool first = true;
  if (hUnfoldLin) {
    hUnfoldLin->SetLineColor(kRed+1);
    hUnfoldLin->SetLineWidth(2);
    hUnfoldLin->SetTitle("Unfolded spectra;E_{true};Counts");
    hUnfoldLin->Draw("HIST");
    first = false;
  }
  if (hUnfoldBay) {
    hUnfoldBay->SetLineColor(kBlue+1);
    hUnfoldBay->SetLineWidth(2);
    hUnfoldBay->SetLineStyle(2);
    if (first) {
      hUnfoldBay->SetTitle("Unfolded spectra;E_{true};Counts");
      hUnfoldBay->Draw("HIST");
      first = false;
    } else {
      hUnfoldBay->Draw("HIST SAME");
    }
  }
  if (hUnfoldDirect) {
    hUnfoldDirect->SetLineColor(kGreen+2);
    hUnfoldDirect->SetLineWidth(2);
    hUnfoldDirect->SetLineStyle(7);
    if (first) {
      hUnfoldDirect->SetTitle("Unfolded spectra;E_{true};Counts");
      hUnfoldDirect->Draw("HIST");
      first = false;
    } else {
      hUnfoldDirect->Draw("HIST SAME");
    }
  }

  TLegend* leg2 = new TLegend(0.55,0.68,0.88,0.88);
  if (hUnfoldLin)    leg2->AddEntry(hUnfoldLin,    "Linear unfolding",   "l");
  if (hUnfoldBay)    leg2->AddEntry(hUnfoldBay,    "Bayesian unfolding", "l");
  if (hUnfoldDirect) leg2->AddEntry(hUnfoldDirect, "Direct unfolding",   "l");
  leg2->Draw();

  // 6.3 χ²/point vs itération (uniquement lin & bayes, direct n'est pas itératif)
  TCanvas* c3 = new TCanvas("cUnfold_Chi2",
                            "Chi2/point vs iteration (Linear vs Bayesian)",
                            900, 700);
  c3->SetGrid();

  bool firstDraw = true;
  if (gChi2Lin) {
    gChi2Lin->SetName("gChi2PerPoint_Linear");
    gChi2Lin->SetTitle("Chi2/point vs iteration;Iteration;Chi2/point");
    gChi2Lin->SetMarkerStyle(20);
    gChi2Lin->SetMarkerColor(kRed+1);
    gChi2Lin->SetLineColor(kRed+1);
    gChi2Lin->Draw("ALP");
    firstDraw = false;
  }

  if (gChi2NDFBay) {
    gChi2NDFBay->SetName("gChi2PerPoint_Bayes");
    gChi2NDFBay->SetMarkerStyle(21);
    gChi2NDFBay->SetMarkerColor(kBlue+1);
    gChi2NDFBay->SetLineColor(kBlue+1);
    if (firstDraw) {
      gChi2NDFBay->SetTitle("Chi2/point vs iteration;Iteration;Chi2/point");
      gChi2NDFBay->Draw("ALP");
      firstDraw = false;
    } else {
      gChi2NDFBay->Draw("LP SAME");
    }
  }

  TLegend* leg3 = new TLegend(0.55,0.75,0.88,0.88);
  if (gChi2Lin)     leg3->AddEntry(gChi2Lin,     "Linear (chi2/point)",  "lp");
  if (gChi2NDFBay)  leg3->AddEntry(gChi2NDFBay,  "Bayesian (chi2/point)","lp");
  leg3->Draw();

  // ------------------ 7) Sauvegarde dans un fichier ROOT ------------------
  TFile* fOut = TFile::Open(outFile, "RECREATE");
  if (!fOut || fOut->IsZombie()) {
    std::cerr << "[WARN] Cannot create output file: " << outFile << std::endl;
  } else {
    fOut->cd();
    // Spectres mesurés
    hMeasFull->Write("hMeas_full");
    hMeasClone->Write("hMeas_used");
    // Unfolded
    if (hUnfoldLin)      hUnfoldLin->Write("hUnfold_Linear");
    if (hUnfoldBay)      hUnfoldBay->Write("hUnfold_Bayes");
    if (hUnfoldDirect)   hUnfoldDirect->Write("hUnfold_Direct");
    if (hUnfoldBayScaled) hUnfoldBayScaled->Write("Scaled_hUnfold_Bayes_perFission");
    // Refolded
    if (hRefoldLin)      hRefoldLin->Write("hRefold_Linear");
    if (hRefoldBay)      hRefoldBay->Write("hRefold_Bayes");
    if (hRefoldDirect)   hRefoldDirect->Write("hRefold_Direct");
    // Résiduel direct
    if (hResDirect)      hResDirect->Write("hResidual_Direct");
    // Graphes chi2
    if (gChi2Lin)     gChi2Lin->Write("gChi2PerPoint_Linear");
    if (gChi2Bay)     gChi2Bay->Write("gChi2_Bayes_total");
    if (gChi2NDFBay)  gChi2NDFBay->Write("gChi2PerPoint_Bayes");
    fOut->Close();
    std::cout << "[INFO] Results written to " << outFile << std::endl;
  }

  // ------------------ 8) Sauvegarde des PNG ------------------
  TString outBase = outFile;
  outBase.ReplaceAll(".root","");

  c1->SaveAs(outBase + "_Meas_vs_Refold.png");
  c2->SaveAs(outBase + "_UnfoldedCompare.png");
  c3->SaveAs(outBase + "_Chi2Compare.png");

  std::cout << "[INFO] PNGs saved with base: " << outBase << "_*.png" << std::endl;

  // ------------------ 9) Nettoyage des fichiers d'entrée ------------------
  fResp->Close();
  fData->Close();
}