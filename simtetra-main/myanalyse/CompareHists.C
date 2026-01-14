// CompareHists.C
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TError.h"
#include "TString.h"
#include <cmath>
#include <memory>
#include <iostream>
#include <vector>

// === Helpers =============================================================

double BinWidth(const TH1* h){
  const auto *ax=h->GetXaxis();
  return (ax->GetXmax()-ax->GetXmin())/ax->GetNbins();
}

bool SameBinWidth(const TH1* a,const TH1* b,double eps=1e-9){
  return std::fabs(BinWidth(a)-BinWidth(b))<eps;
}

bool SameBinning(const TH1* a,const TH1* b,double eps=1e-9){
  const auto *xa=a->GetXaxis(), *xb=b->GetXaxis();
  return xa->GetNbins()==xb->GetNbins()
      && std::fabs(xa->GetXmin()-xb->GetXmin())<eps
      && std::fabs(xa->GetXmax()-xb->GetXmax())<eps;
}

double RangeLen(const TH1* h){
  const auto* x=h->GetXaxis(); return x->GetXmax()-x->GetXmin();
}

// clone + normalise à l’intégrale (aire=1, width-weighted)
TH1* CloneNormIntegral(const TH1* h,const char* newname){
  TH1* c=(TH1*)h->Clone(newname); c->SetDirectory(nullptr);
  const double integ=h->Integral("width"); if(integ!=0) c->Scale(1.0/integ);
  return c;
}

// clone + scale par facteur utilisateur
TH1* CloneScale(const TH1* h,const char* newname,double f){
  TH1* c=(TH1*)h->Clone(newname); c->SetDirectory(nullptr); c->Scale(f);
  return c;
}

// extrait un sous-histo sur range [xmin,xmax] en conservant le binning (uniforme)
TH1* ExtractSubHist(const TH1* h,double xmin,double xmax){
  const auto *ax=h->GetXaxis();
  const int n=ax->GetNbins();
  const double bw=BinWidth(h);
  const double base=ax->GetXmin();
  const int iStart=std::max(1,(int)std::llround((xmin-base)/bw)+1);
  const int iEnd  =std::min(n,(int)std::llround((xmax-base)/bw));
  if(iEnd<iStart) return nullptr;

  const int nSub=iEnd-iStart+1;
  const double xlow=base+(iStart-1)*bw;
  const double xup =base+iEnd*bw;

  TH1* sub=(TH1*)h->Clone(TString(h->GetName())+"_sub");
  sub->SetDirectory(nullptr); sub->Reset();
  sub->SetBins(nSub,xlow,xup);

  for(int i=0;i<nSub;++i){
    const int isrc=iStart+i;
    sub->SetBinContent(i+1, h->GetBinContent(isrc));
    sub->SetBinError  (i+1, h->GetBinError(isrc));
  }
  return sub;
}

// Rebin uniforme des deux hists vers un binning commun (coarsening vers le plus large binwidth)
std::pair<TH1*,TH1*> RebinToCommon(const TH1* h1,const TH1* h2){
  const double bw1=BinWidth(h1), bw2=BinWidth(h2);
  const double bw = std::max(bw1,bw2); // on prend le plus large pour éviter les interpolations
  const double xmin=std::max(h1->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmin());
  const double xmax=std::min(h1->GetXaxis()->GetXmax(), h2->GetXaxis()->GetXmax());
  const int nbins=(int)std::floor((xmax-xmin)/bw);
  if(nbins<=0) return {nullptr,nullptr};

  std::vector<double> edges(nbins+1);
  for(int i=0;i<=nbins;++i) edges[i]=xmin+i*bw;

  // Rebin retourne un nouvel histogramme; on détache du fichier.
  TH1* r1=(TH1*)h1->Rebin(nbins,(TString)h1->GetName()+"_rb", edges.data());
  TH1* r2=(TH1*)h2->Rebin(nbins,(TString)h2->GetName()+"_rb", edges.data());
  if(r1) r1->SetDirectory(nullptr);
  if(r2) r2->SetDirectory(nullptr);
  return {r1,r2};
}

// Dessin overlay + ratio si possible (h1 et h2 doivent avoir exactement le même binning)
void DrawOverlayWithOptionalRatio(TH1* h1, TH1* h2,
                                  const char* ctitle,
                                  const char* l1, const char* l2){
  const bool canRatio = SameBinning(h1,h2);
  if(!canRatio){
    // canvas simple sans ratio
    TCanvas* c=new TCanvas(TString("c_")+ctitle, ctitle, 900, 700);
    h1->Draw("hist");
    h2->Draw("hist same");
    TLegend* leg=new TLegend(0.60,0.74,0.88,0.88);
    leg->AddEntry(h1,l1,"l");
    leg->AddEntry(h2,l2,"l");
    leg->Draw();
    return;
  }

  // Canvas avec 2 pads : haut (spectres) / bas (ratio)
  TCanvas* c=new TCanvas(TString("c_")+ctitle, ctitle, 900, 800);
  c->Divide(1,2);
  TPad* p1=(TPad*)c->cd(1); p1->SetPad(0,0.30,1,1.0); p1->SetBottomMargin(0.05);
  TPad* p2=(TPad*)c->cd(2); p2->SetPad(0,0.00,1,0.30); p2->SetTopMargin(0.08); p2->SetBottomMargin(0.30);
  p1->cd();

  // Overlay
  h1->Draw("hist");
  h2->Draw("hist same");
  TLegend* leg=new TLegend(0.60,0.74,0.88,0.88);
  leg->AddEntry(h1,l1,"l");
  leg->AddEntry(h2,l2,"l");
  leg->Draw();

  // Ratio h1/h2
  p2->cd();
  TH1* r=(TH1*)h1->Clone(TString(h1->GetName())+"_ratio");
  r->SetDirectory(nullptr);
  // éviter divisions par zéro
  for(int i=1;i<=r->GetNbinsX();++i){
    const double d = h2->GetBinContent(i);
    const double n = h1->GetBinContent(i);
    if(d!=0) r->SetBinContent(i, n/d);
    else     r->SetBinContent(i, 0.);
    // Erreur approximative (propagation simple si stats indépendantes)
    double en = h1->GetBinError(i), ed = h2->GetBinError(i);
    if(d!=0) r->SetBinError(i, std::sqrt( (en*en)/(d*d) + (n*n*ed*ed)/(d*d*d*d) ));
    else     r->SetBinError(i, 0.);
  }
  r->SetTitle("Ratio h1/h2");
  r->GetYaxis()->SetTitle("h1 / h2");
  r->GetYaxis()->SetNdivisions(505);
  r->GetYaxis()->SetTitleSize(0.10);
  r->GetYaxis()->SetLabelSize(0.09);
  r->GetXaxis()->SetTitleSize(0.12);
  r->GetXaxis()->SetLabelSize(0.11);
  r->SetMinimum(0.0);
  r->SetMaximum(2.0);
  p2->SetGridy(true);
  r->Draw("E1");
}

// === Main ================================================================

void CompareHists(const char* file1,const char* hist1,
                  const char* file2,const char* hist2,
                  double factor1,double factor2)
{
  gStyle->SetOptStat(0);

  std::unique_ptr<TFile> f1(TFile::Open(file1));
  std::unique_ptr<TFile> f2(TFile::Open(file2));
  if(!f1||f1->IsZombie()){Error("CompareHists","Cannot open %s",file1);return;}
  if(!f2||f2->IsZombie()){Error("CompareHists","Cannot open %s",file2);return;}
  TH1* h1=dynamic_cast<TH1*>(f1->Get(hist1));
  TH1* h2=dynamic_cast<TH1*>(f2->Get(hist2));
  if(!h1||!h2){Error("CompareHists","Histogram(s) not found");return;}

  // ---- A) Normalisation par intégrale + harmonisation (sous-range + éventuel rebin) ----
  TH1 *hA1=nullptr, *hA2=nullptr;  // versions prêtes pour comparaison A

  if(SameBinning(h1,h2)){
    hA1 = CloneNormIntegral(h1,"h1_norm");
    hA2 = CloneNormIntegral(h2,"h2_norm");
  } else {
    // sous-range commun = plus petit range
    const double xmin = std::max(h1->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmin());
    const double xmax = std::min(h1->GetXaxis()->GetXmax(), h2->GetXaxis()->GetXmax());
    TH1* h1sub = ExtractSubHist(h1, xmin, xmax);
    TH1* h2sub = ExtractSubHist(h2, xmin, xmax);

    if(h1sub && h2sub){
      if(SameBinWidth(h1sub,h2sub)){
        hA1 = CloneNormIntegral(h1sub,"h1_norm_sub");
        hA2 = CloneNormIntegral(h2sub,"h2_norm_sub");
      } else {
        // Rebin des sous-histos vers un binwidth commun (coarsening)
        auto rb = RebinToCommon(h1sub,h2sub);
        if(rb.first && rb.second){
          hA1 = CloneNormIntegral(rb.first,"h1_norm_rb");
          hA2 = CloneNormIntegral(rb.second,"h2_norm_rb");
        } else {
          Warning("CompareHists","Unable to harmonize even after sub-range; drawing non-harmonised.");
          hA1 = CloneNormIntegral(h1,"h1_norm_fallback");
          hA2 = CloneNormIntegral(h2,"h2_norm_fallback");
        }
        if(rb.first)  delete rb.first;
        if(rb.second) delete rb.second;
      }
    } else {
      Warning("CompareHists","Sub-range extraction failed; drawing non-harmonised.");
      hA1 = CloneNormIntegral(h1,"h1_norm_fallback2");
      hA2 = CloneNormIntegral(h2,"h2_norm_fallback2");
    }
    if(h1sub) delete h1sub;
    if(h2sub) delete h2sub;
  }

  // Style
  hA1->SetLineColor(kBlue+1); hA1->SetLineWidth(2);
  hA2->SetLineColor(kRed+1);  hA2->SetLineWidth(2);

  // Dessin A : overlay + ratio si possible
  DrawOverlayWithOptionalRatio(hA1, hA2,
      "A) Integrale (harmonisation range/binwidth si besoin)",
      TString::Format("%s::%s", gSystem->BaseName(file1), hist1),
      TString::Format("%s::%s", gSystem->BaseName(file2), hist2));

  // ---- B) Normalisation par facteurs utilisateur (binning libre) -----------------------
  TH1* hB1 = CloneScale(h1,"h1_scaled",factor1);
  TH1* hB2 = CloneScale(h2,"h2_scaled",factor2);
  hB1->SetLineColor(kBlue+2); hB1->SetLineWidth(2);
  hB2->SetLineColor(kRed+2);  hB2->SetLineWidth(2);

  // Ici on tente aussi un ratio si (par chance) binning identique
  DrawOverlayWithOptionalRatio(hB1, hB2,
      "B) Facteurs utilisateur (binning libre)",
      TString::Format("%s (x %.6g)", hist1, factor1),
      TString::Format("%s (x %.6g)", hist2, factor2));

  // Info console
  std::cout << "\n=== INFO ===\n";
  std::cout << "Same binning (raw)? " << (SameBinning(h1,h2)?"yes":"no") << "\n";
  std::cout << "Raw integrals (width): H1=" << h1->Integral("width")
            << "  H2=" << h2->Integral("width") << "\n";
}