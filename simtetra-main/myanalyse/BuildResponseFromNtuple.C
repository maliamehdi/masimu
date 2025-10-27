// BuildResponseFromNtuple.C
//
// Construit les matrices de réponse R_{j i} = P(Emeas_j | Etrue_i)
// depuis un ntuple "resp" (colonnes: eventID, parisIndex, Etrue_keV, Emeas_keV, ...).
//
// Axes des 2D produits ici :
//   X = Emeas [keV], Y = Etrue [keV]
//   => Normalisation "par ligne Y" pour obtenir P(Emeas | Etrue).
//
// Usage (dans ROOT) :
//   .L BuildResponseFromNtuple.C+
//   BuildResponseFromNtuple("../../myanalyse/output_137Cs_run0_smeared.root",
//                           "resp", "RespMatrices.root", "resp_txt",
//                           /*emin=*/0, /*emax=*/8000, /*nbin=*/200,
//                           /*ids=*/{50,70,90,110,130,235,262,278,305});
//
// Le paramètre 'ids' définit l’ordre d’affichage / mapping lisible ; on suppose
// que ton "parisIndex" interne est 0..8 correspondant à ce vecteur 'ids'.
//

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TString.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <map>
#include <memory>

namespace {

void NormalizeRows_Y_as_Conditional(TH2D* hCounts, TH2D* hProb) {
  // hProb = hCounts normalisée par ligne Y (i.e. pour chaque bin Y fixé,
  // on divise les X par la somme de la ligne).
  const int nx = hCounts->GetNbinsX();
  const int ny = hCounts->GetNbinsY();
  for (int iy=1; iy<=ny; ++iy) {
    double rowsum = 0.0;
    for (int ix=1; ix<=nx; ++ix) rowsum += hCounts->GetBinContent(ix, iy);
    for (int ix=1; ix<=nx; ++ix) {
      double v = hCounts->GetBinContent(ix, iy);
      hProb->SetBinContent(ix, iy, (rowsum>0.0 ? v/rowsum : 0.0));
    }
  }
}

void SaveMatrixAsTXT(const TH2D* hProb, const TString& pathTxt) {
  // Écrit une matrice texte : une ligne par bin Y (Etrue), colonnes = bins X (Emeas)
  // En-tête : nx ny  xlow xhigh  ylow yhigh
  std::ofstream out(pathTxt.Data());
  if (!out) {
    std::cerr << "!!! Impossible d'ouvrir en écriture : " << pathTxt << "\n";
    return;
  }
  const int nx = hProb->GetNbinsX();
  const int ny = hProb->GetNbinsY();
  out.setf(std::ios::scientific);
  out << "# nx ny xlow xhigh ylow yhigh\n";
  out << nx << " " << ny << " "
      << hProb->GetXaxis()->GetXmin() << " " << hProb->GetXaxis()->GetXmax() << " "
      << hProb->GetYaxis()->GetXmin() << " " << hProb->GetYaxis()->GetXmax() << "\n";
  out << "# Rows: Y=Etrue bins (low->high), Cols: X=Emeas bins (low->high)\n";

  for (int iy=1; iy<=ny; ++iy) {
    for (int ix=1; ix<=nx; ++ix) {
      double v = hProb->GetBinContent(ix, iy);
      out << v;
      if (ix<nx) out << " ";
    }
    out << "\n";
  }
}

void BuildProfiles(const TH2D* hProb, TH1D*& hMean, TH1D*& hSigma) {
  // Calcule pour chaque Y (Etrue bin) la moyenne et l'écart-type de Emeas
  // sous la distribution P(Emeas|Etrue)
  const int nx = hProb->GetNbinsX();
  const int ny = hProb->GetNbinsY();
  const double ylo = hProb->GetYaxis()->GetXmin();
  const double yhi = hProb->GetYaxis()->GetXmax();

  hMean  = new TH1D(TString(hProb->GetName())+"_MeanEmeas_vs_Etrue",
                    "Mean(Emeas)|Etrue;Etrue [keV];Mean Emeas [keV]",
                    ny, ylo, yhi);
  hSigma = new TH1D(TString(hProb->GetName())+"_SigmaEmeas_vs_Etrue",
                    "Sigma(Emeas)|Etrue;Etrue [keV];Sigma Emeas [keV]",
                    ny, ylo, yhi);

  for (int iy=1; iy<=ny; ++iy) {
    double sumP  = 0.0;
    double sumX  = 0.0;
    double sumX2 = 0.0;
    for (int ix=1; ix<=nx; ++ix) {
      const double p  = hProb->GetBinContent(ix, iy);
      const double xc = hProb->GetXaxis()->GetBinCenter(ix);
      sumP  += p;
      sumX  += p * xc;
      sumX2 += p * xc * xc;
    }
    double mean = 0.0, sig = 0.0;
    if (sumP>0.0) {
      mean = sumX / sumP;
      double var = sumX2/sumP - mean*mean;
      sig = (var>0.0 ? std::sqrt(var) : 0.0);
    }
    hMean ->SetBinContent(iy, mean);
    hSigma->SetBinContent(iy, sig);
  }
}

} // namespace


void BuildResponseFromNtuple(const char* inFile,
                             const char* treeName   = "resp",
                             const char* outRoot    = "RespMatrices.root",
                             const char* outTxtDir  = "resp_txt",
                             double emin_keV        = 0.0,
                             double emax_keV        = 8000.0,
                             int    nbin            = 200,
                             std::vector<int> parisIds = {50,70,90,110,130,235,262,278,305})
{
  // --- ouverture & lecture ntuple ---
  std::unique_ptr<TFile> fin(TFile::Open(inFile,"READ"));
  if (!fin || fin->IsZombie()) { std::cerr << "!!! Impossible d'ouvrir " << inFile << "\n"; return; }

  TTree* T = dynamic_cast<TTree*>(fin->Get(treeName));
  if (!T) { std::cerr << "!!! TTree '" << treeName << "' introuvable\n"; return; }

  // Adresses des branches
  Int_t    parisIndex = -1;
  Double_t Etrue_keV=0, Emeas_keV=0;
  // Les noms doivent matcher exactement tes ntuples :
  if (T->SetBranchAddress("parisIndex", &parisIndex)!=0 ||
      T->SetBranchAddress("Etrue_keV",  &Etrue_keV )!=0 ||
      T->SetBranchAddress("Emeas_keV",  &Emeas_keV )!=0) {
    std::cerr << "!!! Branches manquantes (parisIndex, Etrue_keV, Emeas_keV)\n";
    return;
  }

  // Vérification mapping index->label lisible
  const int nDet = static_cast<int>(parisIds.size());
  std::cout << "Detected " << nDet << " PARIS ids: ";
  for (auto id : parisIds) std::cout << id << " ";
  std::cout << "\n";

  // Création dossier txt
  gSystem->mkdir(outTxtDir, kTRUE);

  // Fichier ROOT de sortie
  std::unique_ptr<TFile> fout(TFile::Open(outRoot,"RECREATE"));
  if (!fout || fout->IsZombie()) { std::cerr << "!!! Impossible d'ouvrir " << outRoot << "\n"; return; }

  // Prépare un conteneur de TH2 par index (0..nDet-1)
  std::vector<TH2D*> hCounts(nDet, nullptr);
  std::vector<TH2D*> hProb(nDet,   nullptr);

  // Noms & titres
  for (int i=0; i<nDet; ++i) {
    const int pid = parisIds[i];
    TString nameC = Form("Rcounts_PARIS%d", pid);
    TString nameP = Form("Rprob_PARIS%d",   pid);
    TString titC  = Form("Counts: Emeas vs Etrue (PARIS%d);E_{meas} [keV];E_{true} [keV]", pid);
    TString titP  = Form("P(Emeas|Etrue) (PARIS%d);E_{meas} [keV];E_{true} [keV]",         pid);

    hCounts[i] = new TH2D(nameC, titC, nbin, emin_keV, emax_keV, nbin, emin_keV, emax_keV);
    hProb[i]   = new TH2D(nameP, titP, nbin, emin_keV, emax_keV, nbin, emin_keV, emax_keV);
    hCounts[i]->SetOption("COLZ");
    hProb[i]  ->SetOption("COLZ");
  }

  // --- Remplissage des 2D bruts (X=Emeas, Y=Etrue) ---
  const Long64_t nentries = T->GetEntries();
  std::cout << "→ Lecture " << nentries << " events ...\n";

  for (Long64_t i=0; i<nentries; ++i) {
    T->GetEntry(i);
    if (parisIndex < 0 || parisIndex >= nDet) continue;
    if (Etrue_keV < emin_keV || Etrue_keV >= emax_keV) continue;
    if (Emeas_keV < emin_keV || Emeas_keV >= emax_keV) continue;

    hCounts[parisIndex]->Fill(Emeas_keV, Etrue_keV, 1.0);
  }

  // --- Normalisation par ligne Y pour obtenir P(Emeas|Etrue) ---
  for (int i=0; i<nDet; ++i) {
    NormalizeRows_Y_as_Conditional(hCounts[i], hProb[i]);
  }

  // --- Sauvegarde des profils (moyenne & sigma) ---
  for (int i=0; i<nDet; ++i) {
    TH1D *hMean=nullptr, *hSigma=nullptr;
    BuildProfiles(hProb[i], hMean, hSigma);

    fout->mkdir(Form("PARIS%d", parisIds[i]));
    fout->cd(Form("PARIS%d", parisIds[i]));
    hCounts[i]->Write();   // brut
    hProb[i]->Write();     // probas
    hMean->Write();
    hSigma->Write();
    fout->cd();

    // TXT pour la probabilité normalisée
    TString pathTxt = Form("%s/Rprob_PARIS%d.txt", outTxtDir, parisIds[i]);
    SaveMatrixAsTXT(hProb[i], pathTxt);
  }

  // On écrit aussi un petit TDirectory avec les métadonnées
  fout->cd();
  TDirectory* meta = fout->mkdir("meta");
  meta->cd();
  TH1D hInfo("info","info",1,0,1);
  hInfo.SetBinContent(1, nDet);
  hInfo.GetXaxis()->SetTitle("dummy");
  hInfo.Write("nDet");
  fout->Write();
  fout->Close();

  std::cout << "✔ Terminé. Matrices ROOT : " << outRoot
            << "  |  TXT par détecteur dans : " << outTxtDir << "\n";
}