// IterativeUnfoldBayes.C
// ----------------------
// Implémentation ROOT d'un unfolding bayésien itératif (D'Agostini / Qi).
// Convention :
//   - hResp : TH2D avec E_meas sur X, E_true sur Y
//   - hMeas : TH1 (mesuré) avec le même binning X que hResp
//
// Signature :
//   TH1D* IterativeUnfoldBayes(const TH1*  hMeas,
//                              const TH2D* hResp,
//                              int   maxIter      = 30,
//                              double tolRelChi2  = 1e-4,
//                              bool  enforcePos   = true,
//                              bool  verbose      = true,
//                              TGraph** gChi2     = nullptr,
//                              TGraph** gChi2NDF  = nullptr);

#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <cmath>

TH1D* IterativeUnfoldBayes(const TH1*  hMeas,
                           const TH2D* hResp,
                           int   maxIter     = 30,
                           double tolRelChi2 = 1e-4,
                           bool  enforcePos  = true,
                           bool  verbose     = true,
                           TGraph** gChi2    = nullptr,
                           TGraph** gChi2NDF = nullptr)
{
  if (!hMeas || !hResp) {
    std::cerr << "[ERROR] IterativeUnfoldBayes: null pointer(s)." << std::endl;
    return nullptr;
  }

  const TH1* hMeas1D = dynamic_cast<const TH1*>(hMeas);
  const TH2D* hR2D   = dynamic_cast<const TH2D*>(hResp);
  if (!hMeas1D || !hR2D) {
    std::cerr << "[ERROR] IterativeUnfoldBayes: wrong histogram types." << std::endl;
    return nullptr;
  }

  // Dimensions : X = Emeas, Y = Etrue
  int nbinsX = hR2D->GetNbinsX(); // measured
  int nbinsY = hR2D->GetNbinsY(); // true

  if (hMeas1D->GetNbinsX() != nbinsX) {
    std::cerr << "[ERROR] IterativeUnfoldBayes: hMeas X-axis bins ("
              << hMeas1D->GetNbinsX() << ") != hResp X-axis bins ("
              << nbinsX << ")." << std::endl;
    return nullptr;
  }

  // Vérifier compatibilité en énergie (borne min/max)
  double xMinMeas = hMeas1D->GetXaxis()->GetXmin();
  double xMaxMeas = hMeas1D->GetXaxis()->GetXmax();
  double xMinResp = hR2D   ->GetXaxis()->GetXmin();
  double xMaxResp = hR2D   ->GetXaxis()->GetXmax();

  if (std::fabs(xMinMeas - xMinResp) > 1e-6 || std::fabs(xMaxMeas - xMaxResp) > 1e-6) {
    std::cerr << "[WARN] IterativeUnfoldBayes: hMeas and hResp X-axis ranges differ: "
              << "[" << xMinMeas << "," << xMaxMeas << "] vs. ["
              << xMinResp << "," << xMaxResp << "]." << std::endl;
    std::cerr << "       -> Les χ² seront calculés seulement sur l'overlap." << std::endl;
  }

  // ----------------------------
  // 1) Extraire g_i et construire R_ij
  // ----------------------------
  std::vector<double> g(nbinsX, 0.0);
  double sumG = 0.0;
  for (int ix = 1; ix <= nbinsX; ++ix) {
    double val = hMeas1D->GetBinContent(ix);
    if (val < 0) val = 0; // sécurité
    g[ix-1] = val;
    sumG   += val;
  }

  // R_ij brut (à partir de hResp)
  std::vector< std::vector<double> > R(nbinsX, std::vector<double>(nbinsY, 0.0));
  for (int ix = 1; ix <= nbinsX; ++ix) {
    for (int iy = 1; iy <= nbinsY; ++iy) {
      R[ix-1][iy-1] = hR2D->GetBinContent(ix, iy);
      if (R[ix-1][iy-1] < 0) R[ix-1][iy-1] = 0.0;
    }
  }

  // Normaliser chaque colonne j (Etrue) pour avoir sum_i R_ij = 1 (PDF(Emeas|Etrue))
  std::vector<double> eff(nbinsY, 0.0);
  for (int j = 0; j < nbinsY; ++j) {
    double colSum = 0.0;
    for (int i = 0; i < nbinsX; ++i) colSum += R[i][j];
    if (colSum > 0) {
      for (int i = 0; i < nbinsX; ++i) R[i][j] /= colSum;
      eff[j] = 1.0; // efficacité "interne" (tout ce qui est dans le TH2)
    } else {
      eff[j] = 0.0;
    }
  }

  // ----------------------------
  // 2) Initialisation de f_j^(0) (prior)
  // ----------------------------
  // Choix simple : uniforme et normalisé à sumG.
  std::vector<double> f(nbinsY, 0.0);
  double initVal = (nbinsY > 0) ? (sumG / nbinsY) : 0.0;
  for (int j = 0; j < nbinsY; ++j) f[j] = initVal;

  // Pour stocker la précédente itération si on veut
  std::vector<double> f_prev(nbinsY, 0.0);

  // ----------------------------
  // 3) Préparation des graphes χ²
  // ----------------------------
  TGraph* grChi2    = nullptr;
  TGraph* grChi2NDF = nullptr;
  if (gChi2)    grChi2    = new TGraph();
  if (gChi2NDF) grChi2NDF = new TGraph();

  auto computeChi2 = [&](const std::vector<double>& fcur, double& chi2, double& chi2NDF) {
    chi2    = 0.0;
    chi2NDF = 0.0;
    int nUsed = 0;
    // t_i = sum_j R_ij * f_j
    for (int i = 0; i < nbinsX; ++i) {
      double ti = 0.0;
      for (int j = 0; j < nbinsY; ++j) {
        if (eff[j] <= 0) continue;
        ti += R[i][j] * fcur[j];
      }
      double gi = g[i];
      // on ne compte que les bins avec un peu de stats (ou prédiction non nulle)
      if (gi > 0 || ti > 0) {
        double denom = (gi > 0 ? gi : 1.0);
        double diff  = (gi - ti);
        chi2 += (diff*diff)/denom;
        ++nUsed;
      }
    }
    chi2NDF = (nUsed > 0 ? chi2 / nUsed : 0.0);
  };

  // ----------------------------
  // 4) Boucle d'itération bayésienne
  // ----------------------------
  double prevChi2NDF = -1.0;
  int    nIterDone   = 0;

  for (int k = 0; k < maxIter; ++k) {
    f_prev = f; // sauvegarde

    // 4.1) Calcul des dénominateurs D_i = sum_j R_ij f_j^(k)
    std::vector<double> D(nbinsX, 0.0);
    for (int i = 0; i < nbinsX; ++i) {
      double sum = 0.0;
      for (int j = 0; j < nbinsY; ++j) {
        if (eff[j] <= 0) continue;
        sum += R[i][j] * f_prev[j];
      }
      D[i] = sum;
    }

    // 4.2) Mise à jour de f_j^(k+1)
    std::vector<double> f_new(nbinsY, 0.0);
    for (int j = 0; j < nbinsY; ++j) {
      if (eff[j] <= 0) {
        f_new[j] = 0.0;
        continue;
      }
      double sumOverMeas = 0.0;
      for (int i = 0; i < nbinsX; ++i) {
        if (D[i] <= 0) continue;  // si prédiction nulle, on ne met pas à jour avec ce bin
        double w_ij = R[i][j] * f_prev[j] / D[i]; // P(true j | meas i)
        sumOverMeas += w_ij * g[i];
      }
      // facteur 1/eff[j] : ici eff[j] ~ 1 (car on a renormalisé), mais on garde la forme générale
      f_new[j] = (eff[j] > 0 ? sumOverMeas / eff[j] : 0.0);
    }

    // 4.3) Enforcement positivité s'il y a des soucis num
    if (enforcePos) {
      for (int j = 0; j < nbinsY; ++j) {
        if (f_new[j] < 0) f_new[j] = 0.0;
      }
    }

    // 4.4) Re-normaliser sur la somme totale mesurée (facultatif mais pratique)
    double sumF = 0.0;
    for (int j = 0; j < nbinsY; ++j) sumF += f_new[j];
    if (sumF > 0 && sumG > 0) {
      double factor = sumG / sumF;
      for (int j = 0; j < nbinsY; ++j) f_new[j] *= factor;
    }

    // Remplacer f par f_new
    f = f_new;
    nIterDone = k+1;

    // 4.5) Calcul du χ² pour cette itération
    double chi2 = 0.0, chi2NDF = 0.0;
    computeChi2(f, chi2, chi2NDF);

    if (verbose) {
      std::cout << "[Bayes] Iter " << (k+1)
                << " : chi2 = " << chi2
                << " , chi2/Npoints = " << chi2NDF << std::endl;
    }

    // Stocker dans les graphes si demandés
    if (grChi2)    grChi2   ->SetPoint(grChi2   ->GetN(), k+1, chi2);
    if (grChi2NDF) grChi2NDF->SetPoint(grChi2NDF->GetN(), k+1, chi2NDF);

    // 4.6) Critère d'arrêt sur chi2NDF
    if (prevChi2NDF > 0 && tolRelChi2 > 0) {
      double relDiff = std::fabs(chi2NDF - prevChi2NDF) /
                       ( (chi2NDF + prevChi2NDF)/2.0 );
      if (relDiff < tolRelChi2) {
        if (verbose) {
          std::cout << "[Bayes] Stopping after " << (k+1)
                    << " iterations (rel Δchi2/Npts = " << relDiff
                    << " < " << tolRelChi2 << ")" << std::endl;
        }
        break;
      }
    }
    prevChi2NDF = chi2NDF;
  }

  // ----------------------------
  // 5) Construire l'histo TH1D de sortie (axe Y du TH2)
  // ----------------------------
  TH1D* hUnfold = new TH1D("hUnfold_Bayes",
                           "Unfolded spectrum (Bayesian);E_{true};Counts",
                           nbinsY,
                           hR2D->GetYaxis()->GetXbins()->GetArray());
  // Si Y n'a pas de binning variable :
  if (hR2D->GetYaxis()->GetXbins()->GetSize() == 0) {
    hUnfold->SetBins(nbinsY,
                     hR2D->GetYaxis()->GetXmin(),
                     hR2D->GetYaxis()->GetXmax());
  }

  for (int j = 0; j < nbinsY; ++j) {
    hUnfold->SetBinContent(j+1, f[j]);
  }

  if (verbose) {
    std::cout << "[Bayes] Done. Iterations performed: " << nIterDone << std::endl;
  }

  // Sortie des graphes si besoin
  if (gChi2)    *gChi2    = grChi2;
  if (gChi2NDF) *gChi2NDF = grChi2NDF;

  return hUnfold;
}