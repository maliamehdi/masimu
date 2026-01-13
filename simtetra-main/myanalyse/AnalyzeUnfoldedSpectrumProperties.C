#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TParameter.h"
#include "TStyle.h"
#include "TLine.h"
#include <iostream>
#include <vector>
#include <string>

// ------------------------------------------------------------
// Helper: build per-bin uncertainty = sqrt(stat^2 + syst^2)
// stat assumed ~ sqrt(N) (Poisson-like) even if N is float
// syst is fracSyst * N
// ------------------------------------------------------------
static TH1D* BuildUncertaintyFromCounts(const TH1* hCounts,
                                       double fracSyst,
                                       const char* nameOut)
{
  if (!hCounts) return nullptr;
  TH1D* hSig = (TH1D*)hCounts->Clone(nameOut);
  hSig->SetDirectory(nullptr);
  hSig->Reset();

  for (int i=0; i<=hCounts->GetNbinsX()+1; ++i) { // include under/over
    double N = hCounts->GetBinContent(i);
    if (N < 0) N = 0; // safety
    double stat = TMath::Sqrt(N);
    double syst = fracSyst * N;
    double sig  = TMath::Sqrt(stat*stat + syst*syst);
    hSig->SetBinContent(i, sig);
  }
  return hSig;
}

// ------------------------------------------------------------
// Correct propagation for y = N/F
// sigma_y^2 = (sigmaN/F)^2 + (N*sigmaF/F^2)^2
// ------------------------------------------------------------
static void NormalizePerFission(TH1D* h, TH1D* hSigma,
                               double nfissions = 1.926593664e+09,
                               double sigmaFissions=0)
{
  if (!h || !hSigma) return;
  if (nfissions <= 0) return;

  for (int i=0; i<=h->GetNbinsX()+1; ++i) {
    double N     = h->GetBinContent(i);
    double sigN  = hSigma->GetBinContent(i);

    double y = N / nfissions;

    double sigY = 0.0;
    if (nfissions > 0.0) {
      double term1 = (sigN / nfissions);
      double term2 = (N * sigmaFissions / (nfissions*nfissions));
      sigY = TMath::Sqrt(term1*term1 + term2*term2);
    }

    h->SetBinContent(i, y);
    h->SetBinError(i, sigY);
    hSigma->SetBinContent(i, sigY);
  }
}

// ------------------------------------------------------------
// Divide content AND errors by bin width (MeV)
// ------------------------------------------------------------
static void NormalizePerMeV(TH1D* h, TH1D* hSigma)
{
  if (!h || !hSigma) return;
  for (int i=0; i<=h->GetNbinsX()+1; ++i) {
    double bw_keV = h->GetBinWidth(i);
    if (bw_keV <= 0) continue;
    double bw_MeV = bw_keV / 1000.0;

    double y    = h->GetBinContent(i) / bw_MeV;
    double sigY = hSigma->GetBinContent(i) / bw_MeV;

    h->SetBinContent(i, y);
    h->SetBinError(i, sigY);
    hSigma->SetBinContent(i, sigY);
  }
}

// ------------------------------------------------------------
// Compute Ngamma/fission and its error from histogram errors
// (simple quadrature of bin errors, assuming independence)
// For an histogram whose bin contents are "per bin" already.
// ------------------------------------------------------------
static void IntegralWithSigma(const TH1D* h, double& I, double& sigI,
                              int b1=1, int b2=-1)
{
  if (!h) { I=0; sigI=0; return; }
  if (b2 < 0) b2 = h->GetNbinsX();

  I = 0.0;
  double v = 0.0;
  for (int i=b1; i<=b2; ++i) {
    I += h->GetBinContent(i);
    double e = h->GetBinError(i);
    v += e*e;
  }
  sigI = TMath::Sqrt(v);
}

// ------------------------------------------------------------
// NEW: Integral of a density histogram (per MeV)
// I = sum_i (dN/dE)_i * dE_i(MeV)
// sigma^2 = sum_i (sigma_i * dE_i)^2
// ------------------------------------------------------------
static void IntegralWithSigmaWidthMeV(const TH1D* h, double& I, double& sigI,
                                      int b1=1, int b2=-1)
{
  if (!h) { I=0; sigI=0; return; }
  if (b2 < 0) b2 = h->GetNbinsX();

  I = 0.0;
  double v = 0.0;
  for (int i=b1; i<=b2; ++i) {
    double bw_keV = h->GetBinWidth(i);
    if (bw_keV <= 0) continue;
    double dE_MeV = bw_keV / 1000.0;

    double y  = h->GetBinContent(i); // #/(fission.MeV)
    double ey = h->GetBinError(i);

    I += y * dE_MeV;
    v += (ey * dE_MeV) * (ey * dE_MeV);
  }
  sigI = TMath::Sqrt(v);
}

// ------------------------------------------------------------
// Compute mean energy and its uncertainty with proper propagation
// using per-bin errors on N_i (already /fission here).
// A = sum N_i E_i ; B = sum N_i ; Ebar = A/B
// sigmaA^2 = sum (E_i^2 sigmaNi^2)
// sigmaB^2 = sum (sigmaNi^2)
// sigmaEbar^2 ~ (sigmaA/B)^2 + (A*sigmaB/B^2)^2
// ------------------------------------------------------------
static void MeanEnergyWithSigma(const TH1D* h, double& Ebar_keV, double& sigEbar_keV,
                                int b1=1, int b2=-1)
{
  if (!h) { Ebar_keV=0; sigEbar_keV=0; return; }
  if (b2 < 0) b2 = h->GetNbinsX();

  double A = 0.0;
  double B = 0.0;
  double varA = 0.0;
  double varB = 0.0;

  for (int i=b1; i<=b2; ++i) {
    double Ni  = h->GetBinContent(i);
    double sig = h->GetBinError(i);
    double Ei  = h->GetBinCenter(i);

    A += Ni * Ei;
    B += Ni;

    varA += (Ei*Ei) * (sig*sig);
    varB += (sig*sig);
  }

  if (B <= 0.0) { Ebar_keV=0; sigEbar_keV=0; return; }

  Ebar_keV = A / B;

  double sigA = TMath::Sqrt(varA);
  double sigB = TMath::Sqrt(varB);

  double term1 = (sigA / B);
  double term2 = (A * sigB / (B*B));
  sigEbar_keV = TMath::Sqrt(term1*term1 + term2*term2);
}

// ------------------------------------------------------------
// NEW: ratio bin-by-bin: hNum/hDen with safe handling for zero denom
// Also returns mean ratio over [b1,b2] ignoring denom<=0 or num<0.
// By default uses simple arithmetic mean over valid bins.
// ------------------------------------------------------------
static TH1D* MakeSafeRatio(const TH1D* hNum, const TH1D* hDen,
                           const char* nameOut,
                           double& meanRatio,
                           int b1=1, int b2=-1)
{
  meanRatio = 0.0;
  if (!hNum || !hDen) return nullptr;
  if (hNum->GetNbinsX() != hDen->GetNbinsX()) return nullptr;

  TH1D* hR = (TH1D*)hNum->Clone(nameOut);
  hR->SetDirectory(nullptr);
  hR->Reset();
  hR->Sumw2(false);

  if (b2 < 0) b2 = hNum->GetNbinsX();
  b1 = std::max(1, b1);
  b2 = std::min(hNum->GetNbinsX(), b2);

  int nValid = 0;
  double sum = 0.0;

  for (int i=0; i<=hNum->GetNbinsX()+1; ++i) {
    double num = hNum->GetBinContent(i);
    double den = hDen->GetBinContent(i);

    if (den > 0.0) {
      double r = num / den;
      hR->SetBinContent(i, r);

      // mean only over requested physical range
      if (i >= b1 && i <= b2) {
        sum += r;
        nValid++;
      }
    } else {
      hR->SetBinContent(i, 0.0);
    }
  }

  meanRatio = (nValid > 0) ? (sum / nValid) : 0.0;
  return hR;
}

// ============================================================
// MAIN MACRO
// ============================================================
void AnalyzeUnfoldedSpectrumProperties(const char* file1="new_UnfoldCompare_PARIS235.root",
                                      const char* histName="hUnfold_Direct",
                                      const char* outFile = "new_PARIS235_Direct_SpectrumProperties.root",
                                      const char* file2 = nullptr,
                                      int   zeroBelowBin = 1,
                                      double fracSystGamma = 0.025,
                                      double nfissions = 1.926593664e9,
                                      double fracSystFission = 0.019,
                                      int bIntMin = 1,
                                      int bIntMax = -1,
                                      bool doPerMeV = true,
                                      // ---- NEW (ratios) ----
                                      const char* dataFileForRatios = "new_UnfoldCompare_PARIS235.root",
                                      const char* hMeasUsedName = "hMeas_used",
                                      const char* hRefoldName   = "hRefold_Direct")
{
  gStyle->SetOptStat(0);

  // ------------------- Load input(s) -------------------
  TFile* f1 = TFile::Open(file1, "READ");
  if (!f1 || f1->IsZombie()) { std::cerr << "[ERROR] Cannot open " << file1 << "\n"; return; }

  TH1* h1_in = dynamic_cast<TH1*>(f1->Get(histName));
  if (!h1_in) { std::cerr << "[ERROR] Cannot find hist '" << histName << "' in " << file1 << "\n"; f1->Close(); return; }

  TH1D* h1 = (TH1D*)h1_in->Clone("init_spectrum_1");
  h1->SetDirectory(nullptr);

  TH1D* h2 = nullptr;
  if (file2) {
    TFile* f2 = TFile::Open(file2, "READ");
    if (!f2 || f2->IsZombie()) { std::cerr << "[ERROR] Cannot open " << file2 << "\n"; f1->Close(); return; }

    TH1* h2_in = dynamic_cast<TH1*>(f2->Get(histName));
    if (!h2_in) { std::cerr << "[ERROR] Cannot find hist '" << histName << "' in " << file2 << "\n"; f2->Close(); f1->Close(); return; }

    h2 = (TH1D*)h2_in->Clone("init_spectrum_2");
    h2->SetDirectory(nullptr);
    f2->Close();
  }

  // ============================================================
  // NEW: Load measured used histogram and refolded histogram
  // and build bin-by-bin ratios
  // ============================================================
  TH1D* hMeasUsed = nullptr;
  TH1D* hRefold   = nullptr;

  TH1D* hRatio_RefoldOverMeas = nullptr;
  TH1D* hRatio_MeasOverRefold = nullptr;
  double mean_R_over_M = 0.0;
  double mean_M_over_R = 0.0;

  TFile* fRat = nullptr;
  if (dataFileForRatios && dataFileForRatios[0] != '\0') {
    fRat = TFile::Open(dataFileForRatios, "READ");
    if (!fRat || fRat->IsZombie()) {
      std::cerr << "[WARN] Cannot open dataFileForRatios=" << dataFileForRatios
                << " -> ratios disabled.\n";
      fRat = nullptr;
    } else {
      TH1* hM_in = dynamic_cast<TH1*>(fRat->Get(hMeasUsedName));
      TH1* hR_in = dynamic_cast<TH1*>(fRat->Get(hRefoldName));

      if (!hM_in || !hR_in) {
        std::cerr << "[WARN] Cannot find hMeasUsedName='" << hMeasUsedName
                  << "' or hRefoldName='" << hRefoldName
                  << "' in " << dataFileForRatios << " -> ratios disabled.\n";
      } else {
        hMeasUsed = dynamic_cast<TH1D*>(hM_in->Clone("hMeas_used_clone_forRatio"));
        hRefold   = dynamic_cast<TH1D*>(hR_in->Clone("hRefold_clone_forRatio"));
        if (hMeasUsed) hMeasUsed->SetDirectory(nullptr);
        if (hRefold)   hRefold->SetDirectory(nullptr);

        if (!hMeasUsed || !hRefold) {
          std::cerr << "[WARN] hMeasUsed or hRefold is not TH1D after clone -> ratios disabled.\n";
        } else if (hMeasUsed->GetNbinsX() != hRefold->GetNbinsX()) {
          std::cerr << "[WARN] Ratio hist nbins mismatch: meas=" << hMeasUsed->GetNbinsX()
                    << " refold=" << hRefold->GetNbinsX() << " -> ratios disabled.\n";
        } else {
          hRatio_RefoldOverMeas = MakeSafeRatio(hRefold,   hMeasUsed,
                                                "hRatio_RefoldOverMeas",
                                                mean_R_over_M,
                                                bIntMin, bIntMax);

          hRatio_MeasOverRefold = MakeSafeRatio(hMeasUsed, hRefold,
                                                "hRatio_MeasOverRefold",
                                                mean_M_over_R,
                                                bIntMin, bIntMax);

          std::cout << "[INFO] Ratio means over bins [" << bIntMin << "," << (bIntMax<0? hMeasUsed->GetNbinsX():bIntMax) << "] :\n"
                    << "       Refold/Meas = " << mean_R_over_M << "\n"
                    << "       Meas/Refold = " << mean_M_over_R << "\n";
        }
      }
    }
  }
  // ------------------- Canvases check -------------------
  TCanvas* c1 = new TCanvas("c1", "Initial spectrum 1", 700, 600);
  c1->SetLogy(false);
  h1->Draw("HIST");

  TCanvas* c2 = new TCanvas("c2", "Initial spectrum 2 (optional)", 700, 600);
  if (h2) h2->Draw("HIST");

  // ------------------- Sum (or copy) -------------------
  TH1D* summed = (TH1D*)h1->Clone("summed_spectrum");
  summed->SetDirectory(nullptr);
  if (h2) summed->Add(h2);

  for (int i=0; i<zeroBelowBin; ++i) {
    summed->SetBinContent(i, 0.0);
  }

  // Uncertainties on counts (stat + syst)
  TH1D* sigCounts = BuildUncertaintyFromCounts(summed, fracSystGamma, "uncertainty_spectrum_counts");

  TCanvas* c3 = new TCanvas("c3", "Summed spectrum (counts)", 700, 600);
  summed->Draw("HIST");

  // ------------------- Normalize per fission -------------------
  double statF = TMath::Sqrt(nfissions);
  double systF = fracSystFission * nfissions;
  double sigF  = TMath::Sqrt(statF*statF + systF*systF);

  TH1D* final = (TH1D*)summed->Clone("final_spectrum_per_fission");
  final->SetDirectory(nullptr);

  TH1D* sigFinal = (TH1D*)sigCounts->Clone("uncertainty_final_spectrum_per_fission");
  sigFinal->SetDirectory(nullptr);

  NormalizePerFission(final, sigFinal, nfissions*2, sigF);

  // ------------------- Properties BEFORE /MeV (per-bin) -------------------
  double ng_bins=0, sigNg_bins=0;
  IntegralWithSigma(final, ng_bins, sigNg_bins, bIntMin, bIntMax);

  double Ebar=0, sigEbar=0;
  MeanEnergyWithSigma(final, Ebar, sigEbar, bIntMin, bIntMax);

  // ------------------- Convert to /MeV and compute integral AFTER /MeV -------------------
  double ng_fromDensity = ng_bins;
  double sigNg_fromDensity = sigNg_bins;

  if (doPerMeV) {
    NormalizePerMeV(final, sigFinal);
    final->GetYaxis()->SetTitle("#gamma / (fission.MeV)");

    // This is the integral you want after /MeV
    IntegralWithSigmaWidthMeV(final, ng_fromDensity, sigNg_fromDensity, bIntMin, bIntMax);
  } else {
    final->GetYaxis()->SetTitle("#gamma / fission");
  }

  final->GetXaxis()->SetTitle("Energy (keV)");

  double Etot_MeV = ng_fromDensity * (Ebar/1000.0);
  double sigEtot_MeV = TMath::Sqrt( TMath::Power((Ebar/1000.0)*sigNg_fromDensity,2)
                                 + TMath::Power(ng_fromDensity*(sigEbar/1000.0),2) );

  std::cout << "========================================\n";
  std::cout << "Histogram used: " << histName << "\n";
  std::cout << "nfissions = " << nfissions << " +/- " << sigF << "\n";
  std::cout << "Ngamma/fission (sum of bins, before /MeV) = " << ng_bins << " +/- " << sigNg_bins << "\n";
  if (doPerMeV) {
    std::cout << "Ngamma/fission (AFTER /MeV, using bin widths) = "
              << ng_fromDensity << " +/- " << sigNg_fromDensity << "\n";
  }
  std::cout << "Mean energy per gamma = " << Ebar << " +/- " << sigEbar << " keV\n";
  std::cout << "Total energy released = " << Etot_MeV << " +/- " << sigEtot_MeV << " MeV\n";
  std::cout << "========================================\n";

  // ------------------- Plot final spectrum -------------------
  TCanvas* c4 = new TCanvas("c4", "Final spectrum", 700, 600);
  c4->SetLogy(true);
  final->Draw("PE1");
    // ============================================================
  // NEW: Plot ratios (Refold/Meas and Meas/Refold)
  // ============================================================
  TCanvas* cR1 = nullptr;
  TCanvas* cR2 = nullptr;

  if (hRatio_RefoldOverMeas && hRatio_MeasOverRefold) {
    // 1) Refold / Meas
    cR1 = new TCanvas("cRatio_RefoldOverMeas", "Ratio: Refold / Meas (bin-by-bin)", 800, 600);
    cR1->SetGrid();
    hRatio_RefoldOverMeas->SetTitle("Bin-by-bin ratio;Energy bin;Refold / Meas");
    hRatio_RefoldOverMeas->SetLineWidth(2);
    hRatio_RefoldOverMeas->Draw("HIST");

    TLegend* legR1 = new TLegend(0.50, 0.78, 0.88, 0.88);
    legR1->SetBorderSize(0);
    legR1->SetFillStyle(0);
    legR1->AddEntry(hRatio_RefoldOverMeas,
                    Form("Refold / Meas ; mean=%.4g", mean_R_over_M),
                    "l");
    legR1->Draw();

    TLine* l1 = new TLine(hRatio_RefoldOverMeas->GetXaxis()->GetXmin(), 1.0,
                          hRatio_RefoldOverMeas->GetXaxis()->GetXmax(), 1.0);
    l1->SetLineStyle(2);
    l1->Draw("SAME");

    // 2) Meas / Refold
    cR2 = new TCanvas("cRatio_MeasOverRefold", "Ratio: Meas / Refold (bin-by-bin)", 800, 600);
    cR2->SetGrid();
    hRatio_MeasOverRefold->SetTitle("Bin-by-bin ratio;Energy bin;Meas / Refold");
    hRatio_MeasOverRefold->SetLineWidth(2);
    hRatio_MeasOverRefold->Draw("HIST");

    TLegend* legR2 = new TLegend(0.50, 0.78, 0.88, 0.88);
    legR2->SetBorderSize(0);
    legR2->SetFillStyle(0);
    legR2->AddEntry(hRatio_MeasOverRefold,
                    Form("Meas / Refold ; mean=%.4g", mean_M_over_R),
                    "l");
    legR2->Draw();

    TLine* l2 = new TLine(hRatio_MeasOverRefold->GetXaxis()->GetXmin(), 1.0,
                          hRatio_MeasOverRefold->GetXaxis()->GetXmax(), 1.0);
    l2->SetLineStyle(2);
    l2->Draw("SAME");
  }
  // ------------------- Save outputs -------------------
  TFile* fout = TFile::Open(outFile, "RECREATE");
  if (fout && !fout->IsZombie()) {
    fout->cd();
    h1->Write("init_spectrum_1");
    if (h2) h2->Write("init_spectrum_2");

    summed->Write("summed_spectrum_counts");
    sigCounts->Write("uncertainty_summed_counts");

    final->Write("final_spectrum_per_fission_per_MeV");
    sigFinal->Write("uncertainty_final_spectrum");

    // store both integrals
    TParameter<double> pNgBins   ("Ngamma_per_fission_fromBinSum",        ng_bins);
    TParameter<double> pSigNgBins("Sigma_Ngamma_per_fission_fromBinSum",  sigNg_bins);

    TParameter<double> pNgDens   ("Ngamma_per_fission_fromDensity",       ng_fromDensity);
    TParameter<double> pSigNgDens("Sigma_Ngamma_per_fission_fromDensity", sigNg_fromDensity);

    TParameter<double> pEbar ("MeanEnergy_keV",            Ebar);
    TParameter<double> pSigE ("Sigma_MeanEnergy_keV",      sigEbar);
    TParameter<double> pEtot ("TotalEnergy_MeV",           Etot_MeV);
    TParameter<double> pSigEt("Sigma_TotalEnergy_MeV",     sigEtot_MeV);

    pNgBins.Write();
    pSigNgBins.Write();
    pNgDens.Write();
    pSigNgDens.Write();
    pEbar.Write();
    pSigE.Write();
    pEtot.Write();
    pSigEt.Write();
        // --- NEW: write ratios + canvases (if available) ---
    if (hMeasUsed) hMeasUsed->Write("hMeas_used_forRatio");
    if (hRefold)   hRefold->Write("hRefold_forRatio");

    if (hRatio_RefoldOverMeas) hRatio_RefoldOverMeas->Write();
    if (hRatio_MeasOverRefold) hRatio_MeasOverRefold->Write();

    TParameter<double> pMeanRoverM("MeanRatio_RefoldOverMeas", mean_R_over_M);
    TParameter<double> pMeanMoverR("MeanRatio_MeasOverRefold", mean_M_over_R);
    pMeanRoverM.Write();
    pMeanMoverR.Write();

    if (cR1) cR1->Write();
    if (cR2) cR2->Write();

    fout->Close();
    std::cout << "[INFO] Wrote outputs to: " << outFile << "\n";
  }
  if (fRat) fRat->Close();

  f1->Close();
}