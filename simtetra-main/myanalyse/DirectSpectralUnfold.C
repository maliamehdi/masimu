// DirectSpectralUnfold.C  (VERSION "NO-NORM" + efficacité conservée)
//
// Convention :
//  - hResp: TH2D avec X=Emeas (i), Y=Etrue (j)
//  - On NE normalise PAS les colonnes.
//  - On construit A(i,j) = hResp(i,j) / NgenPerBinTrue  (ici NgenPerBinTrue = 1e8)
//    => sum_i A(i,j) = eff(j) = efficacité (probabilité de détecter quelque chose dans la matrice)
//  - L'unfold direct et le refold utilisent A partout.
//
// Remarque : cet algorithme "direct" est heuristique (stripping). Garder l'efficacité
// est cohérent via A, mais ce n'est pas une inversion générale comme Bayes.

#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"
#include "TArrayD.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include <algorithm>

TH1D* DirectSpectralUnfold(const TH1*  hMeas,
                           const TH2D* hResp,
                           const char* nameOut,
                           bool        enforcePos,
                           bool        verbose,
                           TH1D**      hResidualOut,
                           TH1D**      hRefoldOut)
{
  if (!hMeas || !hResp) {
    std::cerr << "[DirectSpectralUnfold] ERROR: null input histogram(s).\n";
    return nullptr;
  }

  const TAxis* axMeas  = hMeas->GetXaxis();
  const TAxis* axRespX = hResp->GetXaxis();
  const TAxis* axRespY = hResp->GetYaxis();

  const int nbinsMeas  = axMeas->GetNbins();
  const int nbinsRespX = axRespX->GetNbins();
  const int nbinsTrue  = axRespY->GetNbins();

  if (nbinsMeas != nbinsRespX) {
    std::cerr << "[DirectSpectralUnfold] ERROR: hMeas nbinsX=" << nbinsMeas
              << " != hResp nbinsX=" << nbinsRespX << "\n";
    return nullptr;
  }

  const int N = std::min(nbinsRespX, nbinsTrue);

  // ---- Construire hTrue (axe Y) ----
  TH1D* hTrue = nullptr;
  if (const TArrayD* yb = axRespY->GetXbins(); yb && yb->GetSize() == nbinsTrue+1) {
    hTrue = new TH1D(nameOut, "Direct unfolded spectrum;E_{true};Counts",
                     nbinsTrue, yb->GetArray());
  } else {
    hTrue = new TH1D(nameOut, "Direct unfolded spectrum;E_{true};Counts",
                     nbinsTrue, axRespY->GetXmin(), axRespY->GetXmax());
  }
  hTrue->SetDirectory(nullptr);
  hTrue->Reset();

  // ---- Copier le mesuré dans source[1..N] ----
  std::vector<double> source(N+1, 0.0);
  for (int i = 1; i <= N; ++i) source[i] = std::max(0.0, hMeas->GetBinContent(i));

  // ---- Construire matrice A non normalisée + efficacité eff[j] ----
  // A[i][j] = P(Emeas bin i | Etrue bin j) *avec efficacité*
  const double NgenPerBinTrue = 1e07; // <-- ton cas : même stat par bin vrai

  std::vector< std::vector<double> > A(N+1, std::vector<double>(N+1, 0.0));
  std::vector<double> eff(N+1, 0.0);

  for (int j = 1; j <= N; ++j) {
    double s = 0.0;
    for (int i = 1; i <= N; ++i) {
      double v = hResp->GetBinContent(i, j) / NgenPerBinTrue;
      if (v < 0) v = 0.0;
      A[i][j] = v;
      s += v;
    }
    eff[j] = s; // sum_i A(i,j) = efficacité du bin vrai j
  }

  if (verbose) {
    int jtest = std::min(N, 5);
    std::cout << "[DirectSpectralUnfold] Using NON-normalized response A = hResp/1e8.\n";
    std::cout << "  Example: eff[" << jtest << "]=" << eff[jtest]
              << " ; A(" << jtest << "," << jtest << ")=" << A[jtest][jtest] << "\n";
  }

  // ---- Résiduel hRes ----
  TH1D* hRes = nullptr;
  {
    const TArrayD* xb = axMeas->GetXbins();
    if (xb && xb->GetSize() == nbinsMeas+1) {
      hRes = new TH1D("hResidual_Direct",
                      "Residual after direct unfolding;E_{meas};Counts",
                      nbinsMeas, xb->GetArray());
    } else {
      hRes = new TH1D("hResidual_Direct",
                      "Residual after direct unfolding;E_{meas};Counts",
                      nbinsMeas, axMeas->GetXmin(), axMeas->GetXmax());
    }
    hRes->SetDirectory(nullptr);
    hRes->Reset();
    for (int i = 1; i <= N; ++i) hRes->SetBinContent(i, source[i]);
  }

  // ---- Boucle direct unfolding (haut -> bas) ----
  const double epsCut = 1e-12;

  for (int i = N; i >= 1; --i) {
    double Si = source[i];
    if (Si <= 0.0) continue;

    double Adiag = A[i][i];
    if (Adiag <= 0.0) {
      if (verbose) std::cout << "[DirectSpectralUnfold] i=" << i << " Adiag<=0 skip\n";
      continue;
    }

    // Estimation du "vrai" dans le bin i (avec matrice absolue A)
    double norm = Si / Adiag;

    // Option "efficacité" : si tu veux estimer l'émis, tu peux corriger par eff[i]
    // -> à activer uniquement si tu assumes que A inclut toute la perte d'efficacité
    if (eff[i] > 0.0) norm /= eff[i];

    if (enforcePos && norm < 0.0) norm = 0.0;
    hTrue->SetBinContent(i, norm);

    if (verbose) {
      std::cout << "[DirectSpectralUnfold] i=" << i
                << " Etrue~" << axRespY->GetBinCenter(i)
                << " norm=" << norm
                << " (eff=" << eff[i] << ")\n";
    }

    // Soustraction : source[j] -= norm * A[j][i]
    for (int j = 1; j <= i; ++j) {
      double newv = source[j] - norm * A[j][i];
      if (enforcePos && newv < epsCut) newv = 0.0;
      source[j] = newv;
      hRes->SetBinContent(j, newv);
    }
  }

  // ---- Refold cohérent : g'_i = sum_j A(i,j) * f_j ----
  TH1D* hRefold = nullptr;
  {
    const TArrayD* xb = axMeas->GetXbins();
    if (xb && xb->GetSize() == nbinsMeas+1) {
      hRefold = new TH1D("hRefold_Direct",
                         "Refolded (Direct, no-norm, eff kept);E_{meas};Counts",
                         nbinsMeas, xb->GetArray());
    } else {
      hRefold = new TH1D("hRefold_Direct",
                         "Refolded (Direct, no-norm, eff kept);E_{meas};Counts",
                         nbinsMeas, axMeas->GetXmin(), axMeas->GetXmax());
    }
    hRefold->SetDirectory(nullptr);
    hRefold->Reset();

    for (int i = 1; i <= N; ++i) {
      double sum = 0.0;
      for (int j = 1; j <= N; ++j) {
        double ftrue = hTrue->GetBinContent(j);
        if (ftrue <= 0.0) continue;
        sum += A[i][j] * ftrue;
      }
      hRefold->SetBinContent(i, sum);
    }
  }

  if (hResidualOut) *hResidualOut = hRes;
  if (hRefoldOut)   *hRefoldOut  = hRefold;

  return hTrue;
}