// RunUnfolding.C
//
// Driver pour lancer DEUX unfoldings en parallèle :
//   1) Méthode linéaire itérative (IterativeUnfoldLinear)
//   2) Méthode bayésienne itérative (IterativeUnfoldBayes)
//
// + Comparaison :
//   - g_mesuré (rebinné) vs g_replié_lin vs g_replié_bayes
//   - f_unfold_lin vs f_unfold_bayes
//   - χ²/point vs itération (lin vs bayes)
//
// Prérequis :
//   .L IterativeUnfoldLinear.C+
//   .L IterativeUnfoldBayes.C+
//   .L RunUnfolding.C+
//
// Usage typique :
//   RunUnfolding("Response_PARIS70.root", "hResp",
//                "resolution_bckgndsub_promptsumEvents_PARIS_corrected.root",
//                "sum_Res11keVGamma_PARIS70",
//                "UnfoldCompare_PARIS70.root");
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
#include <TString.h>
#include <iostream>
#include <vector>
#include <cmath>

// ----------------------------------------------------------------------
// Déclarations des fonctions d'unfolding (implémentées ailleurs)
// ⚠️ PAS d'arguments par défaut ici, juste les signatures exactes.
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

// ----------------------------------------------------------------------
// Helper : reconstruire le spectre replié g’(Emeas) = Σ_j R_ij f_j
// à partir de hUnfold (axe Y de hResp) et de la matrice de réponse hResp.
// On renormalise les colonnes comme dans la méthode bayésienne.
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

  // Construire R_ij normalisé colonne par colonne
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
// Driver principal
// ----------------------------------------------------------------------
void RunUnfolding(const char* responseFile      = "Response_PARIS70.root",
                  const char* hRespName         = "hResp",
                  const char* dataFile          = "resolution_bckgndsub_promptsumEvents_PARIS_corrected.root",
                  const char* hMeasName         = "sum_Res11keVGamma_PARIS70",
                  const char* outFile           = "UnfoldCompare_PARIS70.root",
                  // paramètres linéaire
                  int   maxIterLin              = 50,
                  int   minIterLin              = 3,
                  double relChi2TolLin          = 1e-3,
                  bool  normalizeResponseLin    = true,
                  bool  enforcePosLin           = true,
                  // paramètres bayésien
                  int   maxIterBay              = 30,
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

  // ------------------ 2) Rebin + coupure du spectre mesuré ------------------
  //
  // On projette hMeasOrig sur le binning X de la matrice de réponse (E_meas).
  // -> même nbins, mêmes bords, et coupure automatique à [Xmin, Xmax] de hResp.

  TAxis* axResp = hResp->GetXaxis();
  int    nMeasBins = axResp->GetNbins();
  double xMinResp  = axResp->GetXmin();
  double xMaxResp  = axResp->GetXmax();

  const TArrayD* xb = axResp->GetXbins();
  TH1D* hMeas = nullptr;

  if (xb && xb->GetSize() == nMeasBins+1) {
    hMeas = new TH1D("hMeas_used",
                     "Measured (rebinned to response X);E_{meas};Counts",
                     nMeasBins, xb->GetArray());
  } else {
    hMeas = new TH1D("hMeas_used",
                     "Measured (rebinned to response X);E_{meas};Counts",
                     nMeasBins, xMinResp, xMaxResp);
  }
  hMeas->SetDirectory(nullptr);

  // copie complète pour info/sortie
  TH1* hMeasFull = (TH1*)hMeasOrig->Clone("hMeas_full");
  hMeasFull->SetDirectory(nullptr);

  // Projection de hMeasOrig sur le binning X de la réponse
  for (int ib = 1; ib <= hMeasOrig->GetNbinsX(); ++ib) {
    double c = hMeasOrig->GetBinContent(ib);
    double x = hMeasOrig->GetBinCenter(ib);
    if (c <= 0.0) continue;
    if (x < xMinResp || x >= xMaxResp) continue; // coupe explicite
    hMeas->Fill(x, c);
  }

  std::cout << "[INFO] Measured spectrum rebinned:"
            << " nbins = " << hMeas->GetNbinsX()
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

  // Repliage bayésien (matrix * f_bayes)
  TH1D* hRefoldBay = nullptr;
  if (hUnfoldBay) {
    hRefoldBay = BuildRefoldedFromUnfold(hUnfoldBay, hResp, "hRefold_bayes");
  }

  // ------------------ 5) Dessins de comparaison ------------------

  // 5.1 Mesuré (rebinné) vs replié linéaire et bayésien
  TCanvas* c1 = new TCanvas("cUnfold_MeasVsRefold",
                            "Measured vs Refolded (Linear vs Bayesian)",
                            900, 700);
  c1->SetGrid();

  TH1* hMeasClone = (TH1*)hMeas->Clone("hMeas_used_clone");
  hMeasClone->SetDirectory(nullptr);
  hMeasClone->SetLineColor(kBlack);
  hMeasClone->SetLineWidth(2);

  hMeasClone->SetTitle("Measured vs Refolded;E_{meas};Counts");
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

  TLegend* leg1 = new TLegend(0.60,0.70,0.88,0.88);
  leg1->AddEntry(hMeasClone, "Measured (rebinned)", "l");
  if (hRefoldLin) leg1->AddEntry(hRefoldLin, "Refolded (Linear)", "l");
  if (hRefoldBay) leg1->AddEntry(hRefoldBay, "Refolded (Bayesian)", "l");
  leg1->Draw();

  // 5.2 Spectres unfolded linéaire vs bayésien
  TCanvas* c2 = new TCanvas("cUnfold_Spectra",
                            "Unfolded spectra (Linear vs Bayesian)",
                            900, 700);
  c2->SetGrid();

  if (hUnfoldLin) {
    hUnfoldLin->SetLineColor(kRed+1);
    hUnfoldLin->SetLineWidth(2);
    hUnfoldLin->SetTitle("Unfolded spectra;E_{true};Counts");
    hUnfoldLin->Draw("HIST");
  }
  if (hUnfoldBay) {
    hUnfoldBay->SetLineColor(kBlue+1);
    hUnfoldBay->SetLineWidth(2);
    hUnfoldBay->SetLineStyle(2);
    if (hUnfoldLin) hUnfoldBay->Draw("HIST SAME");
    else {
      hUnfoldBay->SetTitle("Unfolded spectrum (Bayesian);E_{true};Counts");
      hUnfoldBay->Draw("HIST");
    }
  }

  TLegend* leg2 = new TLegend(0.60,0.75,0.88,0.88);
  if (hUnfoldLin) leg2->AddEntry(hUnfoldLin, "Linear unfolding", "l");
  if (hUnfoldBay) leg2->AddEntry(hUnfoldBay, "Bayesian unfolding", "l");
  leg2->Draw();

  // 5.3 χ²/point vs itération
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
    if (firstDraw) {
      gChi2Lin->Draw("ALP");
      firstDraw = false;
    } else {
      gChi2Lin->Draw("LP SAME");
    }
  }

  // Pour le bayésien on trace χ²/point (gChi2NDFBay)
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

  TLegend* leg3 = new TLegend(0.60,0.75,0.88,0.88);
  if (gChi2Lin)     leg3->AddEntry(gChi2Lin,     "Linear (chi2/point)", "lp");
  if (gChi2NDFBay)  leg3->AddEntry(gChi2NDFBay,  "Bayesian (chi2/point)", "lp");
  leg3->Draw();

  // ------------------ 6) Sauvegarde dans un fichier ROOT ------------------
  TFile* fOut = TFile::Open(outFile, "RECREATE");
  if (!fOut || fOut->IsZombie()) {
    std::cerr << "[WARN] Cannot create output file: " << outFile << std::endl;
  } else {
    fOut->cd();
    // Spectres mesurés
    hMeasFull->Write("hMeas_full");
    hMeasClone->Write("hMeas_used");
    // Unfolded
    if (hUnfoldLin) hUnfoldLin->Write("hUnfold_Linear");
    if (hUnfoldBay) hUnfoldBay->Write("hUnfold_Bayes");
    // Refolded
    if (hRefoldLin) hRefoldLin->Write("hRefold_Linear");
    if (hRefoldBay) hRefoldBay->Write("hRefold_Bayes");
    // Graphes chi2
    if (gChi2Lin)     gChi2Lin->Write("gChi2PerPoint_Linear");
    if (gChi2Bay)     gChi2Bay->Write("gChi2_Bayes_total");
    if (gChi2NDFBay)  gChi2NDFBay->Write("gChi2PerPoint_Bayes");
    fOut->Close();
    std::cout << "[INFO] Results written to " << outFile << std::endl;
  }

  // ------------------ 7) Sauvegarde des PNG ------------------
  TString outBase = outFile;
  outBase.ReplaceAll(".root","");

  c1->SaveAs(outBase + "_Meas_vs_Refold.png");
  c2->SaveAs(outBase + "_UnfoldedCompare.png");
  c3->SaveAs(outBase + "_Chi2Compare.png");

  std::cout << "[INFO] PNGs saved with base: " << outBase << "_*.png" << std::endl;

  // ------------------ 8) Nettoyage des fichiers d'entrée ------------------
  fResp->Close();
  fData->Close();
}