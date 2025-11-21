// IterativeUnfoldLinear.C
//
// Unfolding itératif linéaire robuste à partir d'un TH1 mesuré et d'un TH2 de réponse.
//
// Hypothèses :
//   - hResp : TH2D avec X = Emeas, Y = Etrue
//   - hMeas : TH1 avec même binning que l'axe X de hResp
//
// Usage typique :
//   .L IterativeUnfoldLinear.C+
//   TFile fR("Response_PARIS50.root");
//   auto hResp = (TH2D*)fR.Get("hResp");
//   TFile fD("data.root");
//   auto hMeas = (TH1*)fD.Get("hMeas");
//
//   TH1D* hRefold = nullptr;
//   TGraph* gChi2  = nullptr;
//
//   TH1D* hUnfold = IterativeUnfoldLinear(hMeas, hResp,
//                                         50,    // maxIter
//                                         true,  // normalizeResponse (colonne par colonne en Etrue)
//                                         true,  // enforcePositivity
//                                         1e-3,  // relChi2Tol (sur chi2/point)
//                                         3,     // minIter
//                                         &hRefold,
//                                         &gChi2);
//
//   // Exemple pour tracer chi2/point vs itération :
//   TCanvas* cChi2 = new TCanvas("cChi2","Chi2 per point vs iteration",800,600);
//   gChi2->SetMarkerStyle(20);
//   gChi2->Draw("ALP");

#include <TH1.h>
#include <TH2.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TMath.h>
#include <TGraph.h>
#include <iostream>

// -----------------------------------------------------------------------------
// Normalisation colonne par colonne de la matrice de réponse (copie du TH2)
//
// Convention :
//   - Axe X = Emeas (i)
//   - Axe Y = Etrue (j)
//   - R_ij = P(Emeas bin i | Etrue bin j)
//
// On normalise donc chaque colonne (Etrue=j) pour que :
//   sum_i R_ij = 1
// -----------------------------------------------------------------------------
TH2D* NormalizeResponseColumns(const TH2D* hResp, bool verbose = true)
{
  if (!hResp) return nullptr;

  auto* hNorm = (TH2D*)hResp->Clone(Form("%s_norm", hResp->GetName()));
  hNorm->SetDirectory(nullptr);

  int nx = hNorm->GetNbinsX(); // Emeas
  int ny = hNorm->GetNbinsY(); // Etrue

  for (int jy = 1; jy <= ny; ++jy) {
    double colSum = 0.0;
    for (int ix = 1; ix <= nx; ++ix) {
      colSum += hNorm->GetBinContent(ix, jy);
    }
    if (colSum > 0.0) {
      for (int ix = 1; ix <= nx; ++ix) {
        double v = hNorm->GetBinContent(ix, jy);
        hNorm->SetBinContent(ix, jy, v / colSum);
      }
    } else {
      if (verbose) {
        std::cerr << "[WARN] NormalizeResponseColumns: column j=" << jy
                  << " (Etrue bin) has zero sum; left unchanged.\n";
      }
    }
  }

  return hNorm;
}

// -----------------------------------------------------------------------------
// Unfolding itératif linéaire (Laszlo-like) avec arrêt automatique.
//
// Entrées :
//   - hMeas  : TH1 mesuré g(Emeas)
//   - hResp  : TH2 de réponse R(Emeas,Etrue) (X=Emeas, Y=Etrue)
//   - maxIter          : nb max d'itérations
//   - normalizeResponse: si true -> normalise R colonne par colonne (Etrue)
//   - enforcePositivity: si true -> tronque f_j < 0 à 0
//   - relChi2Tol       : tolérance relative sur chi2/point pour l'arrêt
//   - minIter          : nb minimal d'itérations avant de pouvoir s'arrêter
//   - hRefoldOut       : si non nul, on remplit avec g' = R f (dernière iter)
//   - gChi2PerPointOut : si non nul, on remplit un TGraph chi2/point vs iter
//
// Sortie :
//   - TH1D* hUnfold : spectre déplié f(Etrue)
// -----------------------------------------------------------------------------
TH1D* IterativeUnfoldLinear(const TH1* hMeas,
                            const TH2D* hResp,
                            int maxIter = 50,
                            bool normalizeResponse = true,
                            bool enforcePositivity = true,
                            double relChi2Tol = 1e-3,
                            int minIter = 3,
                            TH1D** hRefoldOut = nullptr,
                            TGraph** gChi2PerPointOut = nullptr)
{
  if (!hMeas || !hResp) {
    std::cerr << "[ERROR] IterativeUnfoldLinear: null input histogram.\n";
    return nullptr;
  }

  int nMeasBins = hResp->GetNbinsX(); // X = Emeas
  int nTrueBins = hResp->GetNbinsY(); // Y = Etrue

  if (hMeas->GetNbinsX() != nMeasBins) {
    std::cerr << "[ERROR] IterativeUnfoldLinear: hMeas and hResp X-axis mismatch.\n";
    return nullptr;
  }

  if (maxIter <= 0) {
    std::cerr << "[ERROR] IterativeUnfoldLinear: maxIter must be > 0.\n";
    return nullptr;
  }

  if (minIter < 1) minIter = 1;
  if (minIter > maxIter) minIter = maxIter;

  // -- Clone & normaliser la réponse si demandé
  TH2D* hRused = nullptr;
  if (normalizeResponse) {
    hRused = NormalizeResponseColumns(hResp, true);
  } else {
    hRused = (TH2D*)hResp->Clone(Form("%s_used", hResp->GetName()));
    hRused->SetDirectory(nullptr);
  }

  // -- Construire les matrices R (nMeas x nTrue) et R^T (nTrue x nMeas)
  TMatrixD R(nMeasBins, nTrueBins);
  TMatrixD Rt(nTrueBins, nMeasBins);

  for (int ix = 1; ix <= nMeasBins; ++ix) {
    for (int jy = 1; jy <= nTrueBins; ++jy) {
      double v = hRused->GetBinContent(ix, jy);
      int i = ix - 1;
      int j = jy - 1;
      R(i, j)  = v;
      Rt(j, i) = v;
    }
  }

  // -- Construire le vecteur g (mesuré)
  TVectorD g(nMeasBins);
  for (int ix = 1; ix <= nMeasBins; ++ix) {
    g[ix - 1] = hMeas->GetBinContent(ix);
  }

  // -- Initialisation f^(0) : distribution plate avec même intégrale totale
  TVectorD f(nTrueBins);
  double gSum = g.Sum();
  double avg = (nTrueBins > 0) ? gSum / (double)nTrueBins : 0.0;
  for (int j = 0; j < nTrueBins; ++j) {
    f[j] = avg;
  }

  // -- Calcul d'un facteur Ncol pour stabiliser la mise à jour
  double Ncol = 0.0;
  for (int j = 0; j < nTrueBins; ++j) {
    double colSum = 0.0;
    for (int i = 0; i < nMeasBins; ++i) {
      colSum += TMath::Abs(R(i, j));
    }
    if (colSum > Ncol) Ncol = colSum;
  }
  if (Ncol <= 0.0) {
    std::cerr << "[ERROR] IterativeUnfoldLinear: Ncol <= 0, empty response?\n";
    delete hRused;
    return nullptr;
  }

  std::cout << "[INFO] IterativeUnfoldLinear: Ncol = " << Ncol << std::endl;

  // -- Vecteurs de travail
  TVectorD Rf(nMeasBins);
  TVectorD delta(nMeasBins);
  TVectorD corr(nTrueBins);

  double prevChi2PerPoint = -1.0;
  int    usedIter         = maxIter;

  // -- Graph pour chi2/point vs iteration (si demandé)
  TGraph* gChi2PerPoint = nullptr;
  if (gChi2PerPointOut) {
    gChi2PerPoint = new TGraph();
    gChi2PerPoint->SetName("gChi2PerPoint");
    gChi2PerPoint->SetTitle("Chi2 per point vs iteration;Iteration;Chi2 / point");
  }

  // ====================== BOUCLE D'ITÉRATIONS ==========================
  for (int k = 0; k < maxIter; ++k) {

    // 1) Rf = R * f (prédiction mesurée pour f^(k))
    Rf = R * f;

    // 2) delta = g - Rf
    for (int i = 0; i < nMeasBins; ++i) {
      delta[i] = g[i] - Rf[i];
    }

    // 3) corr = (1/Ncol) * R^T * delta
    corr = Rt * delta;
    corr *= (1.0 / Ncol);

    // 4) f^{(k+1)} = f^{(k)} + corr
    for (int j = 0; j < nTrueBins; ++j) {
      f[j] += corr[j];
      if (enforcePositivity && f[j] < 0.0) f[j] = 0.0;
    }

    // 5) Recalcul Rf avec la nouvelle f pour évaluer chi2
    Rf = R * f;

    // 6) Calcul du χ², χ²/ndf et χ²/point
    double chi2 = 0.0;
    int    nPoints = 0;

    for (int ix = 1; ix <= nMeasBins; ++ix) {
      int i = ix - 1;

      double gi   = g[i];
      double gpi  = Rf[i];  // g' = R f

      // Erreur sur gi : si le TH1 a des erreurs, on les prend,
      // sinon on prend sigma^2 = gi (Poisson) ou 1 si gi=0.
      double err    = hMeas->GetBinError(ix);
      double sigma2 = err*err;
      if (sigma2 <= 0.0) {
        sigma2 = (gi > 0.0) ? gi : 1.0;
      }

      if (sigma2 > 0.0) {
        double d = gi - gpi;
        chi2 += (d*d) / sigma2;
        ++nPoints;
      }
    }

    if (nPoints <= 0) {
      std::cerr << "[WARN] IterativeUnfoldLinear: no valid points for chi2 at iter "
                << (k+1) << ". Aborting.\n";
      usedIter = k+1;
      break;
    }

    int    ndf          = (nPoints > 1) ? (nPoints - 1) : 1; // ndf ~ Npoints-1
    double chi2red      = chi2 / (double)ndf;
    double chi2PerPoint = chi2 / (double)nPoints;

    std::cout << "[INFO] Iter " << (k+1)
              << " / " << maxIter
              << " : chi2 = " << chi2
              << ", nPoints = " << nPoints
              << ", chi2/ndf = " << chi2red
              << ", chi2/point = " << chi2PerPoint
              << std::endl;

    // Ajouter au graphe χ²/point vs itération
    if (gChi2PerPoint) {
      int ip = gChi2PerPoint->GetN();
      gChi2PerPoint->SetPoint(ip, (double)(k+1), chi2PerPoint);
    }

    // 7) Critère d'arrêt : après minIter, regarder la variation relative de chi2/point
    if (k+1 >= minIter && prevChi2PerPoint > 0.0) {
      double relDiff = TMath::Abs(chi2PerPoint - prevChi2PerPoint) / prevChi2PerPoint;
      if (relDiff < relChi2Tol) {
        usedIter = k+1;
        std::cout << "[INFO] Early stopping at iteration " << usedIter
                  << " (relative change in chi2/point = " << relDiff
                  << " < " << relChi2Tol << ")\n";
        break;
      }
    }

    prevChi2PerPoint = chi2PerPoint;
  }

  std::cout << "[INFO] IterativeUnfoldLinear: used iterations = " << usedIter
            << " (maxIter = " << maxIter << ")\n";

  // ====================== CONSTRUCTION DES SORTIES =====================

  // -- Spectre déplié hUnfold (axe Etrue)
  TH1D* hUnfold = nullptr;
  {
    const TArrayD* ybins = hResp->GetYaxis()->GetXbins();
    if (ybins && ybins->GetSize() > 0) {
      hUnfold = new TH1D("hUnfold", "Unfolded spectrum;E_{true} [keV];Counts",
                         nTrueBins, ybins->GetArray());
    } else {
      double yMin = hResp->GetYaxis()->GetXmin();
      double yMax = hResp->GetYaxis()->GetXmax();
      hUnfold = new TH1D("hUnfold", "Unfolded spectrum;E_{true} [keV];Counts",
                         nTrueBins, yMin, yMax);
    }
    hUnfold->SetDirectory(nullptr);
    for (int jy = 1; jy <= nTrueBins; ++jy) {
      int j = jy - 1;
      hUnfold->SetBinContent(jy, f[j]);
    }
  }

  // -- Spectre replié g' = R f (dernière itération)
  if (hRefoldOut) {
    TH1D* hRef = (TH1D*)hMeas->Clone("hRefold");
    hRef->Reset();
    hRef->SetDirectory(nullptr);

    // Rf contient déjà R f de la dernière itération
    for (int ix = 1; ix <= nMeasBins; ++ix) {
      int i = ix - 1;
      hRef->SetBinContent(ix, Rf[i]);
    }
    *hRefoldOut = hRef;
  }

  delete hRused;

  // -- Sortie du graphe (ou nettoyage)
  if (gChi2PerPointOut) {
    *gChi2PerPointOut = gChi2PerPoint;   // on le renvoie à l'appelant
  } else if (gChi2PerPoint) {
    delete gChi2PerPoint;                // pas demandé -> on nettoie
  }

  return hUnfold;
}