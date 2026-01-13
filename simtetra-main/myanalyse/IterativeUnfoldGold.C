// IterativeUnfoldGold.C
// ----------------------
// Gold unfolding + option prior (hPrior sur axe TRUE)

#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

TH1D* IterativeUnfoldGold(const TH1*  hMeas,
                          const TH2D* hResp,
                          const TH1*  hPrior = nullptr,   // <-- NEW (optional)
                          double ngenPerEtrue = 1e7,
                          int   maxIter       = 40,
                          double tolRelChi2   = 1e-6,
                          bool  enforcePos    = true,
                          bool  verbose       = true,
                          TGraph** gChi2      = nullptr,
                          TGraph** gChi2NDF   = nullptr)
{
  if (!hMeas || !hResp) {
    std::cerr << "[ERROR] IterativeUnfoldGold: null pointer(s).\n";
    return nullptr;
  }
  if (ngenPerEtrue <= 0) {
    std::cerr << "[ERROR] IterativeUnfoldGold: ngenPerEtrue must be > 0.\n";
    return nullptr;
  }

  const TH1*  hMeas1D = dynamic_cast<const TH1*>(hMeas);
  const TH2D* hR2D    = dynamic_cast<const TH2D*>(hResp);
  if (!hMeas1D || !hR2D) {
    std::cerr << "[ERROR] IterativeUnfoldGold: wrong histogram types.\n";
    return nullptr;
  }

  const int nbinsX = hR2D->GetNbinsX(); // meas
  const int nbinsY = hR2D->GetNbinsY(); // true

  if (hMeas1D->GetNbinsX() != nbinsX) {
    std::cerr << "[ERROR] IterativeUnfoldGold: hMeas nbinsX != hResp nbinsX.\n";
    return nullptr;
  }

  // --- optional prior sanity (must live on TRUE axis) ---
  const TH1* hPrior1D = nullptr;
  if (hPrior) {
    hPrior1D = dynamic_cast<const TH1*>(hPrior);
    if (!hPrior1D) {
      std::cerr << "[ERROR] IterativeUnfoldGold: hPrior provided but not TH1.\n";
      return nullptr;
    }
    if (hPrior1D->GetNbinsX() != nbinsY) {
      std::cerr << "[ERROR] IterativeUnfoldGold: hPrior nbinsX ("
                << hPrior1D->GetNbinsX() << ") != nbinsY ("
                << nbinsY << ").\n";
      return nullptr;
    }
  }

  // r_i (measured)
  std::vector<double> r(nbinsX, 0.0);
  double sumR = 0.0;
  for (int i=1; i<=nbinsX; ++i) {
    double v = hMeas1D->GetBinContent(i);
    if (v < 0) v = 0.0;
    r[i-1] = v;
    sumR  += v;
  }

  // R_ij = hResp(i,j)/ngenPerEtrue
  std::vector< std::vector<double> > R(nbinsX, std::vector<double>(nbinsY, 0.0));
  for (int i=1; i<=nbinsX; ++i) {
    for (int j=1; j<=nbinsY; ++j) {
      double v = hR2D->GetBinContent(i,j) / ngenPerEtrue;
      if (v < 0) v = 0.0;
      R[i-1][j-1] = v;
    }
  }

  // eff_j = Σ_i R_ij
  std::vector<double> eff(nbinsY, 0.0);
  for (int j=0; j<nbinsY; ++j) {
    double s = 0.0;
    for (int i=0; i<nbinsX; ++i) s += R[i][j];
    eff[j] = s;
  }

  // =========================================================
  // Initialization u^(0)
  //   - if hPrior given: u = prior scaled to sumR
  //   - else: keep your original backprojection init
  // =========================================================
  std::vector<double> u(nbinsY, 0.0);

  if (hPrior1D) {
    double sumP = 0.0;
    for (int j=1; j<=nbinsY; ++j) {
      double v = hPrior1D->GetBinContent(j);
      if (v < 0) v = 0.0;
      u[j-1] = v;
      sumP  += v;
    }
    if (sumP > 0.0 && sumR > 0.0) {
      double scale = sumR / sumP;
      for (int j=0; j<nbinsY; ++j) {
        u[j] *= scale;
        if (enforcePos && u[j] < 0) u[j] = 0.0;
      }
    } else {
      // fallback: uniform tiny prior if prior is empty
      double init = (nbinsY > 0 ? sumR / nbinsY : 0.0);
      for (int j=0; j<nbinsY; ++j) u[j] = init;
    }
    if (verbose) std::cout << "[Gold] Using external prior (scaled to sum(meas)).\n";
  } else {
    // original backprojection init: u_j = Σ_i R_ij r_i / eff_j
    for (int j=0; j<nbinsY; ++j) {
      if (eff[j] <= 0) { u[j] = 0.0; continue; }
      double s = 0.0;
      for (int i=0; i<nbinsX; ++i) s += R[i][j] * r[i];
      u[j] = s / eff[j];
      if (enforcePos && u[j] < 0) u[j] = 0.0;
    }
  }

  // Graphes chi2
  TGraph* grChi2    = (gChi2    ? new TGraph() : nullptr);
  TGraph* grChi2NDF = (gChi2NDF ? new TGraph() : nullptr);

  auto computeChi2 = [&](const std::vector<double>& ucur, double& chi2, double& chi2NDF) {
    chi2 = 0.0;
    int nUsed = 0;
    for (int i=0; i<nbinsX; ++i) {
      double fi = 0.0;
      for (int j=0; j<nbinsY; ++j) fi += R[i][j] * ucur[j];
      double ri = r[i];
      if (ri > 0 || fi > 0) {
        double denom = (ri > 0 ? ri : 1.0);
        double diff  = (ri - fi);
        chi2 += (diff*diff)/denom;
        ++nUsed;
      }
    }
    chi2NDF = (nUsed > 0 ? chi2 / nUsed : 0.0);
  };

  double prevChi2NDF = -1.0;
  int nIterDone = 0;

  for (int it=0; it<maxIter; ++it) {
    // fold: f_i = Σ_j R_ij u_j
    std::vector<double> f(nbinsX, 0.0);
    for (int i=0; i<nbinsX; ++i) {
      double s = 0.0;
      for (int j=0; j<nbinsY; ++j) s += R[i][j] * u[j];
      f[i] = s;
    }

    // ratio_i = r_i / f_i (si f_i==0 -> ratio=1)
    std::vector<double> ratio(nbinsX, 1.0);
    for (int i=0; i<nbinsX; ++i) {
      if (f[i] > 0.0) ratio[i] = r[i] / f[i];
      else            ratio[i] = 1.0;
    }

    // correction sur u_j : corr_j = (Σ_i R_ij * ratio_i) / eff_j
    std::vector<double> u_new = u;
    for (int j=0; j<nbinsY; ++j) {
      if (eff[j] <= 0.0) { u_new[j] = 0.0; continue; }
      double s = 0.0;
      for (int i=0; i<nbinsX; ++i) s += R[i][j] * ratio[i];
      double corr = s / eff[j];
      u_new[j] = u[j] * corr;
      if (enforcePos && u_new[j] < 0) u_new[j] = 0.0;
    }

    u = u_new;
    nIterDone = it+1;

    double chi2=0.0, chi2NDFv=0.0;
    computeChi2(u, chi2, chi2NDFv);

    if (verbose) {
      std::cout << "[Gold] Iter " << (it+1)
                << " : chi2=" << chi2
                << " ; chi2/Npts=" << chi2NDFv << "\n";
    }
    if (grChi2)    grChi2   ->SetPoint(grChi2   ->GetN(), it+1, chi2);
    if (grChi2NDF) grChi2NDF->SetPoint(grChi2NDF->GetN(), it+1, chi2NDFv);

    if (prevChi2NDF > 0 && tolRelChi2 > 0) {
      double relDiff = std::fabs(chi2NDFv - prevChi2NDF) / ((chi2NDFv + prevChi2NDF)/2.0);
      if (relDiff < tolRelChi2) {
        if (verbose) {
          std::cout << "[Gold] Stop at iter " << (it+1)
                    << " (rel Δchi2/Npts=" << relDiff << " < " << tolRelChi2 << ")\n";
        }
        break;
      }
    }
    prevChi2NDF = chi2NDFv;
  }

  // output histo (TRUE axis)
  TH1D* hUnfold = new TH1D("hUnfold_Gold",
                           "Unfolded spectrum (Gold);E_{true};Counts",
                           nbinsY,
                           hR2D->GetYaxis()->GetXbins()->GetArray());

  if (hR2D->GetYaxis()->GetXbins()->GetSize() == 0) {
    hUnfold->SetBins(nbinsY, hR2D->GetYaxis()->GetXmin(), hR2D->GetYaxis()->GetXmax());
  }
  hUnfold->SetDirectory(nullptr);

  for (int j=0; j<nbinsY; ++j) hUnfold->SetBinContent(j+1, u[j]);

  if (verbose) std::cout << "[Gold] Done. Iterations performed: " << nIterDone << "\n";

  if (gChi2)    *gChi2    = grChi2;
  if (gChi2NDF) *gChi2NDF = grChi2NDF;

  return hUnfold;
}