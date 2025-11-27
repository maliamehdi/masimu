// DirectSpectralUnfold.C
//
// Implémentation "spectral data" à la Billnert/Oberstedt :
// - travail sur le résiduel d(i), initialement = g_mes(i)
// - on descend en énergie vraie par paquets de 3 bins (j-2, j-1, j)
// - pour chaque triplet, on construit un système 3x3 local
//     A * f = b
//   et on soustrait la contribution alpha_j au résiduel.
//
// Hypothèses :
// - hResp : TH2D avec axe X = E_meas, axe Y = E_true
// - hMeas : TH1 (E_meas) avec même nbinsX que hResp
//
// Signature utilisée par RunUnfolding.C :
//   TH1D* DirectSpectralUnfold(const TH1*  hMeas,
//                              const TH2D* hResp,
//                              const char* nameOut,
//                              bool        enforcePos,
//                              bool        verbose,
//                              TH1D**      hResidualOut,
//                              TH1D**      hRefoldOut);

#include <TH1D.h>
#include <TH2D.h>
#include <TAxis.h>
#include <TArrayD.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <cmath>

// ----------------------------------------------------------------------
// Petit solveur 3x3 par élimination de Gauss avec pivot partiel
// ----------------------------------------------------------------------
static bool Solve3x3(const double A_[3][3], const double b_[3], double x[3])
{
  double A[3][3];
  double b[3];
  for (int i = 0; i < 3; ++i) {
    b[i] = b_[i];
    for (int j = 0; j < 3; ++j) A[i][j] = A_[i][j];
  }

  const double eps = 1e-12;

  // Pivot sur la 1ère colonne
  int piv = 0;
  if (std::fabs(A[1][0]) > std::fabs(A[piv][0])) piv = 1;
  if (std::fabs(A[2][0]) > std::fabs(A[piv][0])) piv = 2;
  if (piv != 0) {
    for (int j = 0; j < 3; ++j) std::swap(A[0][j], A[piv][j]);
    std::swap(b[0], b[piv]);
  }
  if (std::fabs(A[0][0]) < eps) return false;

  // Eliminer ligne 1,2 sur colonne 0
  for (int i = 1; i < 3; ++i) {
    double f = A[i][0] / A[0][0];
    for (int j = 0; j < 3; ++j) A[i][j] -= f * A[0][j];
    b[i] -= f * b[0];
  }

  // Pivot sur la 2ème colonne
  piv = 1;
  if (std::fabs(A[2][1]) > std::fabs(A[piv][1])) piv = 2;
  if (piv != 1) {
    for (int j = 0; j < 3; ++j) std::swap(A[1][j], A[piv][j]);
    std::swap(b[1], b[piv]);
  }
  if (std::fabs(A[1][1]) < eps) return false;

  // Eliminer ligne 2 sur colonne 1
  {
    double f = A[2][1] / A[1][1];
    for (int j = 0; j < 3; ++j) A[2][j] -= f * A[1][j];
    b[2] -= f * b[1];
  }

  if (std::fabs(A[2][2]) < eps) return false;

  // Remontée
  x[2] = b[2] / A[2][2];
  x[1] = (b[1] - A[1][2] * x[2]) / A[1][1];
  x[0] = (b[0] - A[0][1] * x[1] - A[0][2] * x[2]) / A[0][0];

  return true;
}

// ----------------------------------------------------------------------
// DirectSpectralUnfold : version "3 bins par 3 bins"
// ----------------------------------------------------------------------
TH1D* DirectSpectralUnfold(const TH1*  hMeas,
                           const TH2D* hResp,
                           const char* nameOut,
                           bool        enforcePos,
                           bool        verbose,
                           TH1D**      hResidualOut,
                           TH1D**      hRefoldOut)
{
  if (!hMeas || !hResp) {
    std::cerr << "[DirectSpectralUnfold] ERROR: null input histogram(s)." << std::endl;
    return nullptr;
  }

  // Vérifier binning X
  const TAxis* axMeas  = hMeas->GetXaxis();
  const TAxis* axRespX = hResp->GetXaxis();

  int nbinsMeas  = axMeas->GetNbins();
  int nbinsRespX = axRespX->GetNbins();
  if (nbinsMeas != nbinsRespX) {
    std::cerr << "[DirectSpectralUnfold] ERROR: hMeas nbinsX = " << nbinsMeas
              << " != hResp nbinsX = " << nbinsRespX << std::endl;
    return nullptr;
  }

  // Axe Y (E_true) pour le spectre déplié
  const TAxis* axRespY = hResp->GetYaxis();
  int nbinsTrue = axRespY->GetNbins();
  const TArrayD* yb = axRespY->GetXbins();
  TH1D* hTrue = nullptr;
  if (yb && yb->GetSize() == nbinsTrue+1) {
    hTrue = new TH1D(nameOut,
                     "Direct unfolded spectrum;E_{true};Arbitrary units",
                     nbinsTrue, yb->GetArray());
  } else {
    hTrue = new TH1D(nameOut,
                     "Direct unfolded spectrum;E_{true};Arbitrary units",
                     nbinsTrue, axRespY->GetXmin(), axRespY->GetXmax());
  }
  hTrue->SetDirectory(nullptr);

  // Résiduel initial = spectre mesuré
  TH1D* hRes = (TH1D*)hMeas->Clone("hResidual_direct");
  hRes->SetDirectory(nullptr);

  if (verbose) {
    std::cout << "[DirectSpectralUnfold] Start unfolding with "
              << nbinsTrue << " true bins and "
              << nbinsMeas << " measured bins." << std::endl;
  }

  // On stocke f_j dans un tableau temporaire (au cas où on veuille l'exploiter)
  std::vector<double> fTrue(nbinsTrue, 0.0);

  // Boucle sur les true bins par groupes de 3 : j0 = j-2, j1 = j-1, j2 = j
  // On commence à nbinsTrue, on descend par pas de 3.
  for (int j2 = nbinsTrue; j2 >= 1; j2 -= 3) {
    int j1 = j2 - 1;
    int j0 = j2 - 2;

    // Gérer les bords : si j2<3, on peut réduire à 2 bins ou 1 bin
    int nUnknowns = 0;
    int idx[3]    = {-1, -1, -1}; // indices des bins vrais utilisés

    if (j0 >= 1) {
      nUnknowns = 3;
      idx[0] = j0; idx[1] = j1; idx[2] = j2;
    } else if (j1 >= 1) {
      nUnknowns = 2;
      idx[0] = j1; idx[1] = j2;
    } else {
      nUnknowns = 1;
      idx[0] = j2;
    }

    if (verbose) {
      std::cout << "[DirectSpectralUnfold] Solving group: ";
      for (int u = 0; u < nUnknowns; ++u)
        std::cout << idx[u] << " ";
      std::cout << std::endl;
    }

    // Construire la petite matrice A (3x3 max) et le vecteur b (3),
    // en choisissant 3 bins mesurés "autour de la diagonale".
    // Stratégie simple : i0 = idx[0], i1 = idx[0]+1, i2 = idx[0]+2
    // en clampant dans [1, nbinsMeas].
    double A[3][3] = {{0.0}};
    double b[3]    = {0.0, 0.0, 0.0};
    int    nEq     = std::min(3, nbinsMeas); // on ne dépassera pas 3 équations

    for (int r = 0; r < nEq; ++r) {
      int i = idx[0] + r;
      if (i < 1) i = 1;
      if (i > nbinsMeas) i = nbinsMeas;

      // b_r = résiduel à i
      b[r] = hRes->GetBinContent(i);

      // A[r][u] = R_{i, idx[u]}
      for (int u = 0; u < nUnknowns; ++u) {
        A[r][u] = hResp->GetBinContent(i, idx[u]);
      }
    }

    // Si moins de 3 eq utiles (ex. nbinsMeas < 3), on s'arrête proprement
    if (nEq < nUnknowns) {
      if (verbose) {
        std::cout << "[DirectSpectralUnfold] Not enough equations for group, skipping."
                  << std::endl;
      }
      continue;
    }

    // On complète les colonnes restantes de A par 0 pour ne pas polluer Solve3x3
    for (int r = 0; r < 3; ++r) {
      for (int u = nUnknowns; u < 3; ++u) A[r][u] = 0.0;
    }

    double x[3] = {0.0, 0.0, 0.0};
    bool ok = Solve3x3(A, b, x);
    if (!ok) {
      if (verbose) {
        std::cout << "[DirectSpectralUnfold] 3x3 solve failed for group ";
        for (int u = 0; u < nUnknowns; ++u) std::cout << idx[u] << " ";
        std::cout << " -> skipping." << std::endl;
      }
      continue;
    }

    // Enforce positivité + remplir fTrue, hTrue
    for (int u = 0; u < nUnknowns; ++u) {
      double val = x[u];
      if (enforcePos && val < 0.0) val = 0.0;
      int j = idx[u];
      fTrue[j-1] = val;  // stockage interne
      hTrue->SetBinContent(j, val);
    }

    // Mise à jour du résiduel sur TOUTES les bins mesurées :
    // d_i <- d_i - sum_u fTrue(idx[u]) * R_{i, idx[u]}
    for (int i = 1; i <= nbinsMeas; ++i) {
      double old = hRes->GetBinContent(i);
      double sub = 0.0;
      for (int u = 0; u < nUnknowns; ++u) {
        int j = idx[u];
        double Rij = hResp->GetBinContent(i, j);
        if (Rij > 0.0) sub += fTrue[j-1] * Rij;
      }
      double newv = old - sub;
      if (enforcePos && newv < 0.0) newv = 0.0;
      hRes->SetBinContent(i, newv);
    }

    if (verbose) {
      double Ecenter2 = axRespY->GetBinCenter(j2);
      std::cout << "[DirectSpectralUnfold] Group around E_true~"
                << Ecenter2 << " MeV solved, unknowns = ";
      for (int u = 0; u < nUnknowns; ++u)
        std::cout << "f(" << idx[u] << ")=" << fTrue[idx[u]-1] << "  ";
      std::cout << std::endl;
    }
  }

  // Construire un "replié direct" g'_i = sum_j R_ij fTrue_j (sans normalisation)
  TH1D* hRefold = nullptr;
  {
    const TArrayD* xb = axMeas->GetXbins();
    if (xb && xb->GetSize()==nbinsMeas+1) {
      hRefold = new TH1D("hRefold_Direct",
                         "Refolded (Direct);E_{meas};Counts",
                         nbinsMeas, xb->GetArray());
    } else {
      hRefold = new TH1D("hRefold_Direct",
                         "Refolded (Direct);E_{meas};Counts",
                         nbinsMeas, axMeas->GetXmin(), axMeas->GetXmax());
    }
    hRefold->SetDirectory(nullptr);

    for (int i = 1; i <= nbinsMeas; ++i) {
      double sum = 0.0;
      for (int j = 1; j <= nbinsTrue; ++j) {
        double Rij = hResp->GetBinContent(i, j);
        if (Rij <= 0.0) continue;
        sum += Rij * fTrue[j-1];
      }
      hRefold->SetBinContent(i, sum);
    }
  }

  if (hResidualOut) *hResidualOut = hRes;
  if (hRefoldOut)   *hRefoldOut   = hRefold;

  if (verbose) {
    std::cout << "[DirectSpectralUnfold] Done." << std::endl;
  }

  return hTrue;
}