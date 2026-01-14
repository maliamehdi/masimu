// plot_paris.C
//
// Macro ROOT (simple & robuste) utilisant TES 2 arbres :
//   - TTree "resp"       : trace Emeas_keV par PARIS (parisIndex 0..8)
//   - TTree "paris_time" : trace les temps de vol (tFirstCe_ns et/ou tFirstNaI_ns) par PARIS (parisIdx 0..8)
//
// Sorties :
//   - Fichier ROOT : <infile>_parisResp_time_hists.root
//   - PNGs : ./plots/Emeas/  et ./plots/TOF_Ce/  (et ./plots/TOF_NaI/ si dispo)
//
// Usage :
//   root -l -q 'plot_paris.C("output.root")'
//   root -l -q 'plot_paris.C("output.root", 6000, 30000, 2000, -200, 1800)'
//
// Notes :
//   - "TOF" ici = tFirstCe_ns (et/ou tFirstNaI_ns) tel que enregistré dans paris_time.
//   - Si tu veux un vrai TOF relatif à une référence (cathode, source, etc.), il faut stocker cette référence
//     ou calculer Δt dans l’EventAction et l’écrire dans le ntuple.

#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TSystem.h"

// ---------------------------
// Mapping parisIndex -> label
// ---------------------------
// Par défaut, j’utilise ton set d’angles (9 détecteurs):
// {50, 70, 90, 110, 130, 235, 262, 278, 305} deg
static std::map<int, std::string> gParisIdxToLabel = {
  {0, "PARIS50"},
  {1, "PARIS70"},
  {2, "PARIS90"},
  {3, "PARIS110"},
  {4, "PARIS130"},
  {5, "PARIS235"},
  {6, "PARIS262"},
  {7, "PARIS278"},
  {8, "PARIS305"}
};

static TString LabelForIdx(int idx) {
  auto it = gParisIdxToLabel.find(idx);
  if (it != gParisIdxToLabel.end()) return it->second.c_str();
  return Form("PARIS%d", idx);
}

static bool HasBranch(TTree* t, const char* br) {
  return (t && t->GetBranch(br) != nullptr);
}

void plot_paris(const char* infile = "output.root",
                int    nbinsE      = 30000,
                double emax_keV    = 30000.0,
                int    nbinsTOF    = 2000,
                double tofmin_ns   = -200.0,
                double tofmax_ns   = 1800.0)
{
  // ------------------------
  // Ouvrir fichier
  // ------------------------
  TFile* f = TFile::Open(infile, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Impossible d'ouvrir: " << infile << std::endl;
    return;
  }

  // ------------------------
  // Récupérer les TTrees
  // ------------------------
  TTree* tResp = dynamic_cast<TTree*>(f->Get("resp"));
  if (!tResp) {
    std::cerr << "TTree 'resp' introuvable dans " << infile << std::endl;
    f->Close();
    return;
  }

  TTree* tTime = dynamic_cast<TTree*>(f->Get("paris_time"));
  if (!tTime) {
    std::cerr << "TTree 'paris_time' introuvable dans " << infile << std::endl;
    f->Close();
    return;
  }

  std::cout << "TTrees utilisés: resp=" << tResp->GetName()
            << "  paris_time=" << tTime->GetName() << std::endl;

  // ============================================================
  // 1) RESP : Emeas_keV par PARIS
  //    Branches attendues:
  //      eventID (I), parisIndex (I), Emeas_keV (D), ...
  // ============================================================
  if (!HasBranch(tResp, "parisIndex") || !HasBranch(tResp, "Emeas_keV")) {
    std::cerr << "[resp] Branches manquantes. Attendu: parisIndex (I) et Emeas_keV (D)\n";
    std::cerr << "Branches disponibles dans 'resp':\n";
    tResp->Print();
    f->Close();
    return;
  }

  Int_t   resp_eventID   = 0;
  Int_t   resp_parisIdx  = -1;
  Double_t resp_Emeas_keV = 0.0;

  tResp->SetBranchStatus("*", 0);
  if (HasBranch(tResp, "eventID"))    tResp->SetBranchStatus("eventID", 1);
  tResp->SetBranchStatus("parisIndex", 1);
  tResp->SetBranchStatus("Emeas_keV", 1);

  if (HasBranch(tResp, "eventID")) tResp->SetBranchAddress("eventID", &resp_eventID);
  tResp->SetBranchAddress("parisIndex", &resp_parisIdx);
  tResp->SetBranchAddress("Emeas_keV", &resp_Emeas_keV);

  std::map<int, TH1D*> hEmeasByIdx;

  const Long64_t nResp = tResp->GetEntries();
  for (Long64_t i = 0; i < nResp; ++i) {
    tResp->GetEntry(i);

    if (resp_parisIdx < 0) continue;
    if (!(resp_Emeas_keV > 0.0)) continue;

    const TString lbl = LabelForIdx(resp_parisIdx);

    if (!hEmeasByIdx.count(resp_parisIdx)) {
      hEmeasByIdx[resp_parisIdx] =
        new TH1D(Form("hEmeas_%s", lbl.Data()),
                 Form("%s - E_{meas} (keV);E_{meas} [keV];Counts", lbl.Data()),
                 nbinsE, 0.0, emax_keV);
      hEmeasByIdx[resp_parisIdx]->SetLineWidth(2);
    }

    hEmeasByIdx[resp_parisIdx]->Fill(resp_Emeas_keV);
  }

  // ============================================================
  // 2) PARIS_TIME : TOF par PARIS
  //    Branches attendues:
  //      eventID (I), parisIdx (I), tFirstCe_ns (D), tFirstNaI_ns (D)
  //    + (optionnel) Ece_keV / Enai_keV pour filtrer
  // ============================================================
  if (!HasBranch(tTime, "parisIdx") || !HasBranch(tTime, "tFirstCe_ns")) {
    std::cerr << "[paris_time] Branches manquantes. Attendu: parisIdx (I) et tFirstCe_ns (D)\n";
    std::cerr << "Branches disponibles dans 'paris_time':\n";
    tTime->Print();
    f->Close();
    return;
  }

  const bool hasEce   = HasBranch(tTime, "Ece_keV");
  const bool hasEnai  = HasBranch(tTime, "Enai_keV");
  const bool hasTNaI  = HasBranch(tTime, "tFirstNaI_ns");

  Int_t   time_eventID  = 0;
  Int_t   time_parisIdx = -1;
  Double_t time_Ece_keV = 0.0;
  Double_t time_Enai_keV = 0.0;
  Double_t time_tCe_ns  = 0.0;
  Double_t time_tNaI_ns = 0.0;

  tTime->SetBranchStatus("*", 0);
  if (HasBranch(tTime, "eventID")) tTime->SetBranchStatus("eventID", 1);
  tTime->SetBranchStatus("parisIdx", 1);
  if (hasEce)  tTime->SetBranchStatus("Ece_keV", 1);
  if (hasEnai) tTime->SetBranchStatus("Enai_keV", 1);
  tTime->SetBranchStatus("tFirstCe_ns", 1);
  if (hasTNaI) tTime->SetBranchStatus("tFirstNaI_ns", 1);

  if (HasBranch(tTime, "eventID")) tTime->SetBranchAddress("eventID", &time_eventID);
  tTime->SetBranchAddress("parisIdx", &time_parisIdx);
  if (hasEce)  tTime->SetBranchAddress("Ece_keV", &time_Ece_keV);
  if (hasEnai) tTime->SetBranchAddress("Enai_keV", &time_Enai_keV);
  tTime->SetBranchAddress("tFirstCe_ns", &time_tCe_ns);
  if (hasTNaI) tTime->SetBranchAddress("tFirstNaI_ns", &time_tNaI_ns);

  std::map<int, TH1D*> hTOF_Ce_ByIdx;
  std::map<int, TH1D*> hTOF_NaI_ByIdx;

  const Long64_t nTime = tTime->GetEntries();
  for (Long64_t i = 0; i < nTime; ++i) {
    tTime->GetEntry(i);

    if (time_parisIdx < 0) continue;

    const TString lbl = LabelForIdx(time_parisIdx);

    // Création à la volée
    if (!hTOF_Ce_ByIdx.count(time_parisIdx)) {
      hTOF_Ce_ByIdx[time_parisIdx] =
        new TH1D(Form("hTOF_Ce_%s", lbl.Data()),
                 Form("%s - tFirstCe (ns);tFirstCe [ns];Counts", lbl.Data()),
                 nbinsTOF, tofmin_ns, tofmax_ns);
      hTOF_Ce_ByIdx[time_parisIdx]->SetLineWidth(2);

      if (hasTNaI) {
        hTOF_NaI_ByIdx[time_parisIdx] =
          new TH1D(Form("hTOF_NaI_%s", lbl.Data()),
                   Form("%s - tFirstNaI (ns);tFirstNaI [ns];Counts", lbl.Data()),
                   nbinsTOF, tofmin_ns, tofmax_ns);
        hTOF_NaI_ByIdx[time_parisIdx]->SetLineWidth(2);
      }
    }

    // Filtrage simple (optionnel) : ne remplir que si énergie déposée > 0
    // (ça évite des temps "vides" si tu stockes des lignes sans dépôt)
    const bool passCe  = (!hasEce)  ? true : (time_Ece_keV  > 0.0);
    const bool passNaI = (!hasEnai) ? true : (time_Enai_keV > 0.0);

    if (passCe)  hTOF_Ce_ByIdx[time_parisIdx]->Fill(time_tCe_ns);
    if (hasTNaI && passNaI) hTOF_NaI_ByIdx[time_parisIdx]->Fill(time_tNaI_ns);
  }

  // ------------------------
  // Sorties : ROOT + PNG
  // ------------------------
  TString outRoot = infile;
  outRoot.ReplaceAll(".root", "_parisResp_time_hists.root");
  TFile* fout = TFile::Open(outRoot, "RECREATE");

  TDirectory* dE   = fout->mkdir("resp_Emeas");
  TDirectory* dTCe = fout->mkdir("paris_time_tFirstCe");
  TDirectory* dTNaI = hasTNaI ? fout->mkdir("paris_time_tFirstNaI") : nullptr;

  gSystem->mkdir("plots", kTRUE);
  gSystem->mkdir("plots/Emeas", kTRUE);
  gSystem->mkdir("plots/TOF_Ce", kTRUE);
  if (hasTNaI) gSystem->mkdir("plots/TOF_NaI", kTRUE);

  TCanvas cE("cE", "Emeas", 900, 700);
  TCanvas cT("cT", "TOF",   900, 700);

  // Liste triée des indices rencontrés (union des maps)
  std::vector<int> idxs;
  idxs.reserve(32);
  for (const auto& kv : hEmeasByIdx) idxs.push_back(kv.first);
  for (const auto& kv : hTOF_Ce_ByIdx) idxs.push_back(kv.first);
  std::sort(idxs.begin(), idxs.end());
  idxs.erase(std::unique(idxs.begin(), idxs.end()), idxs.end());

  for (int idx : idxs) {
    const TString lbl = LabelForIdx(idx);

    // Emeas (si présent)
    if (hEmeasByIdx.count(idx)) {
      cE.cd();
      hEmeasByIdx[idx]->Draw("HIST");
      cE.SaveAs(Form("plots/Emeas/Emeas_%s.png", lbl.Data()));
      dE->cd();
      hEmeasByIdx[idx]->Write();
    }

    // TOF Ce (si présent)
    if (hTOF_Ce_ByIdx.count(idx)) {
      cT.cd();
      hTOF_Ce_ByIdx[idx]->Draw("HIST");
      cT.SaveAs(Form("plots/TOF_Ce/TOF_Ce_%s.png", lbl.Data()));
      dTCe->cd();
      hTOF_Ce_ByIdx[idx]->Write();
    }

    // TOF NaI (optionnel)
    if (hasTNaI && hTOF_NaI_ByIdx.count(idx)) {
      cT.cd();
      hTOF_NaI_ByIdx[idx]->Draw("HIST");
      cT.SaveAs(Form("plots/TOF_NaI/TOF_NaI_%s.png", lbl.Data()));
      dTNaI->cd();
      hTOF_NaI_ByIdx[idx]->Write();
    }
  }

  fout->Write();
  fout->Close();
  f->Close();

  std::cout << "✓ Terminé.\n"
            << "  - Input : " << infile << "\n"
            << "  - Output ROOT : " << outRoot << "\n"
            << "  - PNGs : ./plots/Emeas/  ./plots/TOF_Ce/"
            << (hasTNaI ? "  ./plots/TOF_NaI/" : "")
            << "\n" << std::endl;
}