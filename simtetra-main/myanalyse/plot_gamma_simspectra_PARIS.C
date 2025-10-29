// plot_paris.C
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TSystem.h"

// Si tu connais un mapping copy->label, renseigne-le ici (facultatif).
// Laisse vide pour utiliser "PARIS<copy>" par défaut.
static std::map<int,std::string> gCopyToLabel = {
  // {0, "PARIS50"},
  // {1, "PARIS70"},
  // ...
};

void plot_paris(const char* infile = "output_152Eu_smeared0.root",
                int nbins = 3000,       // binning par défaut
                double emax_keV = 30000 // xmax en keV
                )
{
  // --- Ouvrir le fichier et récupérer le ntuple ---
  TFile* f = TFile::Open(infile, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Impossible d'ouvrir: " << infile << std::endl;
    return;
  }

  TTree* t = dynamic_cast<TTree*>(f->Get("ParisEdep"));
  if (!t) {
    std::cerr << "TTree 'ParisEdep' introuvable dans " << infile << std::endl;
    return;
  }

  // Branches attendues: eventID (I), copy (I), eCe_keV (D), eNaI_keV (D)
  Int_t   eventID = 0;
  Int_t   copy    = -1;
  Double_t eCe_keV = 0.;
  Double_t eNaI_keV = 0.;

  t->SetBranchAddress("eventID", &eventID);
  t->SetBranchAddress("copy",    &copy);
  t->SetBranchAddress("eCe_keV", &eCe_keV);
  t->SetBranchAddress("eNaI_keV",&eNaI_keV);

  // --- Histos par copie ---
  std::map<int, TH1D*> hCeByCopy;
  std::map<int, TH1D*> hNaIByCopy;

  // --- Boucle d'entrées ---
  const Long64_t nent = t->GetEntries();
  for (Long64_t i=0; i<nent; ++i) {
    t->GetEntry(i);
    // Ignore les events à 0 
    if (eCe_keV <= 0 && eNaI_keV <= 0) continue;

    // Création à la volée si besoin
    if (!hCeByCopy.count(copy)) {
      TString lbl = gCopyToLabel.count(copy) ? gCopyToLabel[copy].c_str()
                                             : Form("PARIS%d", copy);
      hCeByCopy[copy]  = new TH1D(Form("hCe_%s", lbl.Data()),
                                  Form("%s - E_{dep}^{Ce} (keV);E_{dep}^{Ce} [keV];Counts", lbl.Data()),
                                  nbins, 0., emax_keV);
      hNaIByCopy[copy] = new TH1D(Form("hNaI_%s", lbl.Data()),
                                  Form("%s - E_{dep}^{NaI} (keV);E_{dep}^{NaI} [keV];Counts", lbl.Data()),
                                  nbins, 0., emax_keV);
    }

    // Remplissage
    if (eCe_keV  > 0) hCeByCopy[copy]->Fill(eCe_keV);
    if (eNaI_keV > 0) hNaIByCopy[copy]->Fill(eNaI_keV);
  }

  // --- Fichier de sortie (histos + images) ---
  TString outRoot = infile;
  outRoot.ReplaceAll(".root", "_hists.root");
  TFile* fout = TFile::Open(outRoot, "RECREATE");

  // Dossiers de rangement
  TDirectory* dCe  = fout->mkdir("CeBr3");
  TDirectory* dNaI = fout->mkdir("NaI");

  // Dessiner & sauvegarder PNG
  gSystem->mkdir("plots", kTRUE);

  // Canvases multi-plots facultatifs
  TCanvas c1("cCe","Edep Ce", 900, 700);
  TCanvas c2("cNaI","Edep NaI", 900, 700);

  // Dessin par copie (un PNG par histo)
  for (const auto& kv : hCeByCopy) {
    int cpy = kv.first;
    TString lbl = gCopyToLabel.count(cpy) ? gCopyToLabel[cpy].c_str()
                                          : Form("PARIS%d", cpy);

    // Ce
    c1.cd();
    kv.second->SetLineWidth(2);
    kv.second->Draw("HIST");
    c1.SaveAs(Form("plots/EdepCe_%s.png", lbl.Data()));
    dCe->cd();
    kv.second->Write();

    // NaI
    c2.cd();
    TH1D* hN = hNaIByCopy.at(cpy);
    hN->SetLineWidth(2);
    hN->Draw("HIST");
    c2.SaveAs(Form("plots/EdepNaI_%s.png", lbl.Data()));
    dNaI->cd();
    hN->Write();
  }

  fout->Write();
  fout->Close();
  f->Close();

  std::cout << "✓ Terminé. Histos écrits dans: " << outRoot << std::endl
            << "  PNGs dans le dossier ./plots/" << std::endl;
}