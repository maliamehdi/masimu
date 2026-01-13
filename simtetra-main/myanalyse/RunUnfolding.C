// RunUnfolding.C  (ADD: Gold + RAW refold with efficiency + optional prior + Emax cut on measured)
// --------------------------------------------------------

#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TArrayD.h>
#include <TString.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

// --- prototypes ---
// Bayes: ngenPerEtrue + optional prior (TRUE axis) at the end
TH1D* IterativeUnfoldBayes(const TH1*  hMeas,
                           const TH2D* hResp,
                           double ngenPerEtrue,
                           int   maxIter,
                           double tolRelChi2,
                           bool  enforcePos,
                           bool  verbose,
                           TGraph** gChi2,
                           TGraph** gChi2NDF,
                           const TH1*  hPrior /*= nullptr*/);

// Direct
TH1D* DirectSpectralUnfold(const TH1*  hMeas,
                           const TH2D* hResp,
                           const char* nameOut,
                           bool        enforcePos,
                           bool        verbose,
                           TH1D**      hResidualOut,
                           TH1D**      hRefoldOut);

// Gold: ngenPerEtrue + optional prior (TRUE axis) at the end
TH1D* IterativeUnfoldGold(const TH1*  hMeas,
                          const TH2D* hResp,
                          const TH1*  hPrior /*= nullptr*/,
                          double ngenPerEtrue,
                          int   maxIter,
                          double tolRelChi2,
                          bool  enforcePos,
                          bool  verbose,
                          TGraph** gChi2,
                          TGraph** gChi2NDF);

// ----------------- small helpers -----------------
static double SafeIntegral(const TH1* h, int b1=1, int b2=-1) {
  if (!h) return 0.0;
  if (b2 < 0) b2 = h->GetNbinsX();
  double s = 0.0;
  for (int i=b1;i<=b2;++i) s += h->GetBinContent(i);
  return s;
}

// Apply an Emax cut on BIN CENTER: keep only bins with center < Emax_keV (content+error kept), others set to 0.
static void ApplyBinCenterCut(TH1* h, double Emax_keV) {
  if (!h) return;
  if (!(Emax_keV > 0.0)) return; // if Emax<=0 => do nothing
  const int nb = h->GetNbinsX();
  for (int ib=1; ib<=nb; ++ib) {
    const double c = h->GetBinCenter(ib);
    if (c >= Emax_keV) {
      h->SetBinContent(ib, 0.0);
      h->SetBinError  (ib, 0.0);
    }
  }
}

// Count how many bins on an axis satisfy center < Emax_keV.
// If Emax_keV <= 0 => keep all.
static int NbinsWithCenterBelow(const TAxis* ax, double Emax_keV) {
  if (!ax) return 0;
  const int nb = ax->GetNbins();
  if (!(Emax_keV > 0.0)) return nb;
  int n = 0;
  for (int j=1; j<=nb; ++j) {
    if (ax->GetBinCenter(j) < Emax_keV) ++n;
    else break; // assumes increasing centers
  }
  return n;
}

// Compare ONLY bins with axis bin center < Emax_keV, by checking low/up edges bin-by-bin.
// This is robust even if one histogram stores variable edges and the other is uniform.
static bool SameBinningAsAxisUpToEmaxCenter_Edges(const TH1* h,
                                                  const TAxis* ax,
                                                  double Emax_keV,
                                                  double tol=1e-6)
{
  if (!h || !ax) return false;
  const int nbAx = ax->GetNbins();
  const int nbH  = h->GetNbinsX();
  if (nbH != nbAx) return false; // IMPORTANT: we require same nbins overall; we just compare only up-to-Emax

  const TAxis* ah = h->GetXaxis();

  const int Nkeep = NbinsWithCenterBelow(ax, Emax_keV);
  if (Nkeep <= 0) return false;

  for (int ib = 1; ib <= Nkeep; ++ib) {
    const double lowH = ah->GetBinLowEdge(ib);
    const double upH  = ah->GetBinUpEdge(ib);
    const double lowA = ax->GetBinLowEdge(ib);
    const double upA  = ax->GetBinUpEdge(ib);

    if (std::fabs(lowH - lowA) > tol) return false;
    if (std::fabs(upH  - upA) > tol) return false;
  }
  return true;
}

// Make a cleaned prior clone (NO rebin), but binning check is ONLY up-to-Emax on TRUE bin centers.
// Output is ALWAYS a FULL-size TRUE-axis histogram; we copy only bins with center < Emax_keV, rest = 0.
static TH1D* MakePriorCloneNoRebin_UpToEmaxCenterOnly(const TH1* hPriorIn,
                                                      const TH2D* hResp,
                                                      double Emax_keV,
                                                      const char* name="hPrior_trueAxis")
{
  if (!hPriorIn || !hResp) return nullptr;

  const TAxis* ay = hResp->GetYaxis();
  const int nbY   = ay->GetNbins();

  // If Emax<=0 => "UpTo" means "all bins"
  const double Ecut = (Emax_keV > 0.0 ? Emax_keV : -1.0);
  const int Nkeep = NbinsWithCenterBelow(ay, Ecut);
  if (Nkeep <= 0) return nullptr;

  // Require same nbins overall AND matching edges for bins 1..Nkeep
  if (!SameBinningAsAxisUpToEmaxCenter_Edges(hPriorIn, ay, Ecut)) return nullptr;

  // Create FULL TRUE-axis output histogram (same as response Y axis)
  TH1D* h = nullptr;
  const TArrayD* yb = ay->GetXbins();
  if (yb && yb->GetSize() == nbY+1) h = new TH1D(name, name, nbY, yb->GetArray());
  else                             h = new TH1D(name, name, nbY, ay->GetXmin(), ay->GetXmax());
  h->SetDirectory(nullptr);
  h->Reset();

  // Copy only bins 1..Nkeep; rest stays 0
  for (int j=1; j<=Nkeep; ++j) {
    double v = hPriorIn->GetBinContent(j);
    if (v < 0) v = 0.0;
    h->SetBinContent(j, v);
    // optional: copy errors if useful
    // h->SetBinError(j, hPriorIn->GetBinError(j));
  }

  // safety: no negatives anywhere
  for (int j=0; j<=h->GetNbinsX()+1; ++j) {
    if (h->GetBinContent(j) < 0) h->SetBinContent(j, 0.0);
  }
  return h;
}

static void DebugCompareBinningUpToEmax(const TH1* hPriorIn, const TH2D* hResp, double Emax_keV)
{
  if (!hPriorIn || !hResp) return;

  const TAxis* ay = hResp->GetYaxis();
  const int nbY   = ay->GetNbins();

  // Nkeep = nbins TRUE dont le centre < Emax
  int Nkeep = nbY;
  if (Emax_keV > 0.0) {
    Nkeep = 0;
    for (int j=1; j<=nbY; ++j) {
      if (ay->GetBinCenter(j) < Emax_keV) ++Nkeep;
      else break;
    }
  }

  std::cout << "----- DEBUG BINNING (up-to-Emax-center) -----\n";
  std::cout << "Emax_keV=" << Emax_keV << " ; Nkeep(TRUE)=" << Nkeep << "\n";
  std::cout << "prior nbins=" << hPriorIn->GetNbinsX() << "\n";
  std::cout << "resp  Y nbins=" << nbY << "\n";

  const TArrayD* yb = ay->GetXbins();
  const TArrayD* pb = hPriorIn->GetXaxis()->GetXbins();

  std::cout << "resp Y has variable bins? " << (yb && yb->GetSize()>0 ? "YES":"NO") << "\n";
  std::cout << "prior has variable bins?  " << (pb && pb->GetSize()>0 ? "YES":"NO") << "\n";

  // Compare edges in the safest way: use GetBinLowEdge/UpEdge (works for uniform & variable)
  const int Nedge = Nkeep + 1; // edges 0..Nkeep
  const double tol = 1e-9;

  for (int i=1; i<=Nedge; ++i) {
    // Edge i-1 in [lowEdge(1), upEdge(1)=edge1, ...]
    double eResp = (i==1) ? ay->GetBinLowEdge(1) : ay->GetBinUpEdge(i-1);
    double ePrior;
    // if prior has >= Nkeep bins, compare same index; otherwise stop
    if (hPriorIn->GetNbinsX() >= Nkeep) {
      ePrior = (i==1) ? hPriorIn->GetXaxis()->GetBinLowEdge(1) : hPriorIn->GetXaxis()->GetBinUpEdge(i-1);
    } else {
      std::cout << "[DEBUG] prior is truncated (nbins < Nkeep). nbinsPrior=" << hPriorIn->GetNbinsX()
                << " < Nkeep=" << Nkeep << "\n";
      break;
    }

    if (std::fabs(eResp - ePrior) > tol) {
      std::cout << "[MISMATCH] edge index i=" << i
                << "  respY=" << eResp
                << "  prior =" << ePrior
                << "  diff=" << (eResp - ePrior) << "\n";
      // also print the bin center for context (bin i-1)
      if (i>=2) {
        int b = i-1;
        std::cout << "          bin=" << b
                  << " respCenter=" << ay->GetBinCenter(b)
                  << " priorCenter=" << hPriorIn->GetXaxis()->GetBinCenter(b) << "\n";
      }
      std::cout << "--------------------------------------------\n";
      return;
    }
  }

  std::cout << "[DEBUG] OK: edges match up-to-Emax-center (within tol=" << tol << ")\n";
  std::cout << "--------------------------------------------\n";
}

// RAW refold: g'_i = Σ_j (hResp(i,j)/ngenPerEtrue) * u_j
static TH1D* BuildRefoldedFromUnfoldRaw(const TH1D* hUnfold,
                                        const TH2D* hResp,
                                        double ngenPerEtrue,
                                        const char* name = "hRefold_fromUnfold")
{
  if (!hUnfold || !hResp) return nullptr;
  if (ngenPerEtrue <= 0) return nullptr;

  int nbinsX = hResp->GetNbinsX();
  int nbinsY = hResp->GetNbinsY();

  const TAxis* axX = hResp->GetXaxis();
  TH1D* hRef = nullptr;
  const TArrayD* xb = axX->GetXbins();
  if (xb && xb->GetSize() == nbinsX+1) {
    hRef = new TH1D(name, "Refolded spectrum (RAW);E_{meas};Counts", nbinsX, xb->GetArray());
  } else {
    hRef = new TH1D(name, "Refolded spectrum (RAW);E_{meas};Counts", nbinsX, axX->GetXmin(), axX->GetXmax());
  }
  hRef->SetDirectory(nullptr);
  hRef->Reset();

  for (int i=1; i<=nbinsX; ++i) {
    double ti = 0.0;
    for (int j=1; j<=nbinsY; ++j) {
      double u = hUnfold->GetBinContent(j);
      if (u <= 0) continue;
      double Rij = hResp->GetBinContent(i,j) / ngenPerEtrue;
      if (Rij < 0) Rij = 0.0;
      ti += Rij * u;
    }
    hRef->SetBinContent(i, ti);
  }
  return hRef;
}

static TH1D* MakeMeasuredChildFirstNBins(const TH1*  hMeasOrig,
                                         const TH2D* hResp,
                                         const char* name = "hMeas_used")
{
  if (!hMeasOrig || !hResp) return nullptr;

  const TAxis* axResp = hResp->GetXaxis();
  int    Nresp  = axResp->GetNbins();
  const TArrayD* edges = axResp->GetXbins();

  TH1D* h = nullptr;
  if (edges && edges->GetSize() == Nresp+1) {
    h = new TH1D(name,
                 "Measured (first N bins, response X-binning);E_{meas};Counts",
                 Nresp, edges->GetArray());
  } else {
    h = new TH1D(name,
                 "Measured (first N bins, response X-binning);E_{meas};Counts",
                 Nresp, axResp->GetXmin(), axResp->GetXmax());
  }
  h->SetDirectory(nullptr);
  h->Reset();

  int Nmother = hMeasOrig->GetNbinsX();
  int Ncopy   = std::min(Nresp, Nmother);

  for (int ib=1; ib<=Ncopy; ++ib) {
    h->SetBinContent(ib, hMeasOrig->GetBinContent(ib));
    h->SetBinError  (ib, hMeasOrig->GetBinError  (ib));
  }
  return h;
}

// ----------------- main driver -----------------
void RunUnfolding(const char* responseFile      = "new_Response_PARIS235.root",
                  const char* hRespName         = "hResp",
                  const char* dataFile          = "run6CheckEnergyspectra_all.root",//"Res0508_promptsumEvents_PARIS_january.root",
                  const char* hMeasName         = "resbinspectrumPARIS235",//"sum_Res11keVGamma_PARIS235",
                  const char* outFile           = "new_testunfold_Co60_PARIS235.root",//"new_UnfoldCompare_PARIS235.root",
                  // --- optional prior inputs ---
                  const char* priorFile         = "run6CheckEnergyspectra_all.root",//"Res0508_promptsumEvents_PARIS_january.root",
                  const char* hPriorName        = "resbinspectrumPARIS235",//"sum_Res11keVGamma_PARIS235",
                  // measured-energy cut on BIN CENTER (keV). If <=0, no cut. -1 for no cut.
                  double EmaxMeas_keV           = 15000.0,
                  // scaling of response
                  double ngenPerEtrue           = 1e7,
                  // Bayes
                  int   maxIterBay              = 50,
                  double tolRelChi2Bay          = 1e-9,
                  bool  enforcePosBay           = true,
                  bool  verboseBay              = true,
                  // Gold
                  int   maxIterGold             = 50,
                  double tolRelChi2Gold         = 1e-6,
                  bool  enforcePosGold          = true,
                  bool  verboseGold             = true,
                  // Direct
                  bool  enforcePosDirect        = true,
                  bool  verboseDirect           = true,
                  // Physically-motivated scaling
                  double nFission = 729729640)//1926593664)
{
  TFile* fResp = TFile::Open(responseFile, "READ");
  if (!fResp || fResp->IsZombie()) { std::cerr << "[ERROR] Cannot open " << responseFile << "\n"; return; }

  TH2D* hResp = dynamic_cast<TH2D*>(fResp->Get(hRespName));
  if (!hResp) { std::cerr << "[ERROR] TH2D '" << hRespName << "' not found\n"; fResp->Close(); return; }

  // optional: quick efficiency sanity
  {
    int jtest = 50;
    double col = 0;
    for (int i=1;i<=hResp->GetNbinsX();++i) col += hResp->GetBinContent(i,jtest);
    std::cout << "colSum raw=" << col << " ; eff=" << col/ngenPerEtrue << std::endl;
  }

  TFile* fData = TFile::Open(dataFile, "READ");
  if (!fData || fData->IsZombie()) { std::cerr << "[ERROR] Cannot open " << dataFile << "\n"; fResp->Close(); return; }

  TH1* hMeasOrig = dynamic_cast<TH1*>(fData->Get(hMeasName));
  if (!hMeasOrig) { std::cerr << "[ERROR] TH1 '" << hMeasName << "' not found\n"; fResp->Close(); fData->Close(); return; }

  TH1* hMeasFull = (TH1*)hMeasOrig->Clone("hMeas_full");
  hMeasFull->SetDirectory(nullptr);

  TH1D* hMeas = MakeMeasuredChildFirstNBins(hMeasOrig, hResp, "hMeas_used");
  if (!hMeas) { std::cerr << "[ERROR] Failed to build hMeas_used\n"; fResp->Close(); fData->Close(); return; }

  // apply measured-energy cut on BIN CENTER
  if (EmaxMeas_keV > 0.0) {
    ApplyBinCenterCut(hMeas, EmaxMeas_keV);
    std::cout << "[INFO] Applied measured cut: keep bins with center < " << EmaxMeas_keV
              << " keV (bins above set to 0, no rebin)\n";
  }

  // =========================================================
  // Load optional prior histogram (NO REBIN), compare ONLY UpToEmaxCenter
  // =========================================================
  TH1*   hPriorIn = nullptr;
  TH1D*  hPriorY  = nullptr;
  TFile* fPrior   = nullptr;

  if (priorFile && priorFile[0] != '\0' && hPriorName && hPriorName[0] != '\0') {
    fPrior = TFile::Open(priorFile, "READ");
    if (!fPrior || fPrior->IsZombie()) {
      std::cerr << "[WARN] Cannot open priorFile=" << priorFile << " -> no prior used.\n";
      fPrior = nullptr;
    } else {
      hPriorIn = dynamic_cast<TH1*>(fPrior->Get(hPriorName));
      if (!hPriorIn) {
        std::cerr << "[WARN] Cannot find hPriorName='" << hPriorName << "' in " << priorFile
                  << " -> no prior used.\n";
      } else {
        // IMPORTANT: only up-to-EmaxCenter check (no full check), no rebin
        hPriorY = MakePriorCloneNoRebin_UpToEmaxCenterOnly(hPriorIn, hResp, EmaxMeas_keV, "hPrior_trueAxis");
        if (!hPriorY) {
          std::cerr << "[WARN] Prior found but binning does NOT match TRUE axis "
                    << "(ONLY up-to-Emax on TRUE bin centers). No rebin requested -> prior disabled.\n";
        } else {
          const int Nkeep = NbinsWithCenterBelow(hResp->GetYaxis(), EmaxMeas_keV);
          std::cout << "[INFO] Prior loaded (NO rebin, UpToEmaxCenter ONLY): " << priorFile << " : " << hPriorName
                    << " -> using as 'hPrior_trueAxis' (copied bins = " << Nkeep
                    << ", remaining TRUE bins set to 0)\n";
        }
        DebugCompareBinningUpToEmax(hPriorIn, hResp, EmaxMeas_keV);
      }
    }
  }

  // -------- BAYES --------
  TGraph* gChi2Bay    = nullptr;
  TGraph* gChi2NDFBay = nullptr;

  std::cout << "[INFO] Running Bayesian unfolding...\n";
  TH1D* hUnfoldBay = IterativeUnfoldBayes(hMeas, hResp, ngenPerEtrue,
                                         maxIterBay, tolRelChi2Bay,
                                         enforcePosBay, verboseBay,
                                         &gChi2Bay, &gChi2NDFBay,
                                         hPriorY);

  TH1D* hRefoldBay = nullptr;
  if (hUnfoldBay) hRefoldBay = BuildRefoldedFromUnfoldRaw(hUnfoldBay, hResp, ngenPerEtrue, "hRefold_Bayes");

  // -------- GOLD --------
  TGraph* gChi2Gold    = nullptr;
  TGraph* gChi2NDFGold = nullptr;

  std::cout << "[INFO] Running Gold unfolding...\n";
  TH1D* hUnfoldGold = IterativeUnfoldGold(hMeas, hResp, hPriorY, ngenPerEtrue,
                                         maxIterGold, tolRelChi2Gold,
                                         enforcePosGold, verboseGold,
                                         &gChi2Gold, &gChi2NDFGold);

  TH1D* hRefoldGold = nullptr;
  if (hUnfoldGold) hRefoldGold = BuildRefoldedFromUnfoldRaw(hUnfoldGold, hResp, ngenPerEtrue, "hRefold_Gold");

  // -------- DIRECT --------
  TH1D* hResDirect    = nullptr;
  TH1D* hRefoldDirect = nullptr;

  std::cout << "[INFO] Running Direct spectral unfolding...\n";
  TH1D* hUnfoldDirect = DirectSpectralUnfold(hMeas, hResp, "hUnfold_Direct",
                                            enforcePosDirect, verboseDirect,
                                            &hResDirect, &hRefoldDirect);

  // =====================================================================
  // Per-fission outputs
  // =====================================================================
  auto MakePerFission = [&](TH1D* hin, const char* newname) -> TH1D* {
    if (!hin || nFission<=0) return nullptr;
    TH1D* h = (TH1D*)hin->Clone(newname);
    h->SetDirectory(nullptr);
    h->Scale(1.0/nFission);
    return h;
  };

  TH1D* hUnfoldBay_perFission    = MakePerFission(hUnfoldBay,    "hUnfold_Bayes_perFission");
  TH1D* hUnfoldGold_perFission   = MakePerFission(hUnfoldGold,   "hUnfold_Gold_perFission");
  TH1D* hUnfoldDirect_perFission = MakePerFission(hUnfoldDirect, "hUnfold_Direct_perFission");

  TH1D* hRefoldBay_perFission    = MakePerFission(hRefoldBay,    "hRefold_Bayes_perFission");
  TH1D* hRefoldGold_perFission   = MakePerFission(hRefoldGold,   "hRefold_Gold_perFission");
  TH1D* hRefoldDirect_perFission = MakePerFission(hRefoldDirect, "hRefold_Direct_perFission");

  // -------- comparison canvas: Meas vs refolds --------
  TCanvas* c1 = new TCanvas("cUnfold_MeasVsRefold",
                            "Measured(child) vs Refolded (Bayes/Gold/Direct)", 900, 700);
  c1->SetGrid();

  TH1D* hMeasClone = (TH1D*)hMeas->Clone("hMeas_used_clone");
  hMeasClone->SetDirectory(nullptr);
  hMeasClone->SetLineColor(kBlack);
  hMeasClone->SetLineWidth(2);
  hMeasClone->SetTitle("Measured vs Refolded;E_{meas};Counts");
  hMeasClone->Draw("HIST");

  if (hRefoldBay)    { hRefoldBay->SetLineColor(kBlue+1);  hRefoldBay->SetLineStyle(2);  hRefoldBay->SetLineWidth(2); hRefoldBay->Draw("HIST SAME"); }
  if (hRefoldGold)   { hRefoldGold->SetLineColor(kRed+1);  hRefoldGold->SetLineStyle(9); hRefoldGold->SetLineWidth(2); hRefoldGold->Draw("HIST SAME"); }
  if (hRefoldDirect) { hRefoldDirect->SetLineColor(kGreen+2); hRefoldDirect->SetLineStyle(7); hRefoldDirect->SetLineWidth(2); hRefoldDirect->Draw("HIST SAME"); }

  TLegend* leg1 = new TLegend(0.52,0.70,0.88,0.88);
  leg1->AddEntry(hMeasClone, "Measured (child, cut on bin center)", "l");
  if (hRefoldBay)    leg1->AddEntry(hRefoldBay,  "Refold (Bayes, RAW)", "l");
  if (hRefoldGold)   leg1->AddEntry(hRefoldGold, "Refold (Gold, RAW)", "l");
  if (hRefoldDirect) leg1->AddEntry(hRefoldDirect,"Refold (Direct)", "l");
  leg1->Draw();

  // -------- unfolded spectra canvas --------
  TCanvas* c2 = new TCanvas("cUnfold_Spectra", "Unfolded spectra (Bayes/Gold/Direct)", 900, 700);
  c2->SetGrid();

  bool first = true;
  if (hUnfoldBay) {
    hUnfoldBay->SetLineColor(kBlue+1); hUnfoldBay->SetLineStyle(2); hUnfoldBay->SetLineWidth(2);
    hUnfoldBay->SetTitle("Unfolded spectra;E_{true};Counts");
    hUnfoldBay->Draw("HIST");
    first = false;
  }
  if (hUnfoldGold) {
    hUnfoldGold->SetLineColor(kRed+1); hUnfoldGold->SetLineStyle(9); hUnfoldGold->SetLineWidth(2);
    if (first) hUnfoldGold->Draw("HIST"); else hUnfoldGold->Draw("HIST SAME");
    first = false;
  }
  if (hUnfoldDirect) {
    hUnfoldDirect->SetLineColor(kGreen+2); hUnfoldDirect->SetLineStyle(7); hUnfoldDirect->SetLineWidth(2);
    if (first) hUnfoldDirect->Draw("HIST"); else hUnfoldDirect->Draw("HIST SAME");
  }

  TLegend* leg2 = new TLegend(0.55,0.70,0.88,0.88);
  if (hUnfoldBay)    leg2->AddEntry(hUnfoldBay,  "Bayes", "l");
  if (hUnfoldGold)   leg2->AddEntry(hUnfoldGold, "Gold", "l");
  if (hUnfoldDirect) leg2->AddEntry(hUnfoldDirect,"Direct", "l");
  leg2->Draw();

  // -------- save outputs --------
  TFile* fOut = TFile::Open(outFile, "RECREATE");
  if (fOut && !fOut->IsZombie()) {
    fOut->cd();
    hMeasFull->Write("hMeas_full");
    hMeasClone->Write("hMeas_used");

    if (hPriorY) hPriorY->Write("hPrior_trueAxis");

    if (hUnfoldBay)    hUnfoldBay->Write("hUnfold_Bayes");
    if (hUnfoldGold)   hUnfoldGold->Write("hUnfold_Gold");
    if (hUnfoldDirect) hUnfoldDirect->Write("hUnfold_Direct");

    if (hRefoldBay)    hRefoldBay->Write("hRefold_Bayes");
    if (hRefoldGold)   hRefoldGold->Write("hRefold_Gold");
    if (hRefoldDirect) hRefoldDirect->Write("hRefold_Direct");

    if (hUnfoldBay_perFission)    hUnfoldBay_perFission->Write();
    if (hUnfoldGold_perFission)   hUnfoldGold_perFission->Write();
    if (hUnfoldDirect_perFission) hUnfoldDirect_perFission->Write();

    if (hRefoldBay_perFission)    hRefoldBay_perFission->Write();
    if (hRefoldGold_perFission)   hRefoldGold_perFission->Write();
    if (hRefoldDirect_perFission) hRefoldDirect_perFission->Write();

    if (hResDirect) hResDirect->Write("hResidual_Direct");

    if (gChi2Bay)       gChi2Bay->Write("gChi2_Bayes_total");
    if (gChi2NDFBay)    gChi2NDFBay->Write("gChi2PerPoint_Bayes");
    if (gChi2Gold)      gChi2Gold->Write("gChi2_Gold_total");
    if (gChi2NDFGold)   gChi2NDFGold->Write("gChi2PerPoint_Gold");

    fOut->Close();
    std::cout << "[INFO] Results written to " << outFile << "\n";
  }

  if (fPrior) fPrior->Close();
  fResp->Close();
  fData->Close();
}