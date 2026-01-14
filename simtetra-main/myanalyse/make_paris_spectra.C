// make_paris_spectra.C
// Produit un histogramme d'énergie Ce (smearé) pour chaque détecteur PARIS,
// à partir du ntuple "ParisEdep" (sortie de ton code Geant4).

#include <map>
#include <string>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TString.h"

// (Optionnel) Si tu veux des labels explicites plutôt que PARIS0..8 :
static std::map<int,std::string> gCopyToLabel = {
  // {0, "PARIS50"}, {1, "PARIS70"}, {2, "PARIS90"}, ...
};

void make_paris_spectra(const char* infile = "eu152_104_10_7.root",
                        int nbins = 30000,
                        double emax_keV = 30000)
{
  // --- Ouvrir le fichier ---
  TFile* fin = TFile::Open(infile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "❌ Impossible d'ouvrir " << infile << std::endl;
    return;
  }

  TTree* t = dynamic_cast<TTree*>(fin->Get("ParisEdep"));
  if (!t) {
    std::cerr << "❌ TTree 'ParisEdep' introuvable dans " << infile << std::endl;
    fin->Close();
    return;
  }

  // --- Branches ---
  Int_t    eventID = 0;
  Int_t    copy = -1;
  Double_t eCe_keV = 0.;
  Double_t eNaI_keV = 0.;

  t->SetBranchAddress("eventID", &eventID);
  t->SetBranchAddress("copy",    &copy);
  t->SetBranchAddress("eCe_keV", &eCe_keV);
  t->SetBranchAddress("eNaI_keV",&eNaI_keV);

  // --- Création des histogrammes par détecteur ---
  std::map<int, TH1D*> hCeByCopy;

  const Long64_t nent = t->GetEntries();
  std::cout << "→ Lecture de " << nent << " entrées..." << std::endl;

  for (Long64_t i = 0; i < nent; ++i) {
    t->GetEntry(i);
    if (eCe_keV <= 0) continue; // ignorer les zéros

    if (!hCeByCopy.count(copy)) {
      TString lbl = gCopyToLabel.count(copy)
                      ? gCopyToLabel[copy].c_str()
                      : Form("PARIS%d", copy);

      hCeByCopy[copy] = new TH1D(Form("hCe_%s", lbl.Data()),
                                 Form("%s - Edep Ce (smeared);E_{dep}^{Ce} [keV];Counts", lbl.Data()),
                                 nbins, 0., emax_keV);
      hCeByCopy[copy]->Sumw2();
    }

    hCeByCopy[copy]->Fill(eCe_keV);
  }

  // --- Écriture dans un nouveau fichier ROOT (pas de dossier) ---
  TString outRoot = infile;
  outRoot.ReplaceAll(".root", "_PARIS_spectra.root");
  TFile* fout = TFile::Open(outRoot, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "❌ Impossible de créer " << outRoot << std::endl;
    fin->Close();
    return;
  }

  for (auto& kv : hCeByCopy) {
    kv.second->Write();  // directement à la racine du fichier
  }

  fout->Write();
  fout->Close();
  fin->Close();

  std::cout << "✅ Spectres écrits dans : " << outRoot << std::endl;
}