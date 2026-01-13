// WriteRW05Inputs.C
//
// Macro ROOT pour exporter un spectre mesuré et une matrice de réponse
// au format texte simple compatible RW05 (.spe et .rsp).
//
// Usage dans ROOT :
//   .L WriteRW05Inputs.C+
//
//   WriteRW05Inputs("Response_PARIS50.root", "hResp",
//                   "resolution0508_bckgndsub_promptsumEvents_PARIS_corrected.root",
//                   "sum_Res11keVGamma_PARIS50",
//                   "paris50_rw05");   // donnera paris50_rw05.spe et paris50_rw05.rsp
//
// Ensuite :
//   - mettre ces fichiers dans le répertoire où tu lances rw05
//   - lancer rw05 et lui donner ces noms de fichiers aux prompts.

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TArrayD.h>
#include <TString.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

// ----------------------------------------------------------------------
// Helper : construit un "histo fille" qui contient les N premiers bins
// du spectre mesuré, avec N = nbinsX de la matrice de réponse.
// ----------------------------------------------------------------------
static TH1D* MakeMeasuredChildForResponse(const TH1*  hMeasOrig,
                                          const TH2*  hResp,
                                          const char* name = "hMeas_rw05")
{
  if (!hMeasOrig || !hResp) return nullptr;

  const TAxis* axRespX = hResp->GetXaxis();
  int NrespX = axRespX->GetNbins();
  const TArrayD* xedges = axRespX->GetXbins();

  TH1D* h = nullptr;
  if (xedges && xedges->GetSize() == NrespX+1) {
    h = new TH1D(name, "Measured spectrum (RW05 export);E_{meas};Counts",
                 NrespX, xedges->GetArray());
  } else {
    h = new TH1D(name, "Measured spectrum (RW05 export);E_{meas};Counts",
                 NrespX, axRespX->GetXmin(), axRespX->GetXmax());
  }
  h->SetDirectory(nullptr);

  int Nmother = hMeasOrig->GetNbinsX();
  int Ncopy   = std::min(NrespX, Nmother);

  std::cout << "[RW05] Build child measured histogram: NrespX = " << NrespX
            << ", Nmother = " << Nmother
            << ", copying first Ncopy = " << Ncopy << " bins." << std::endl;

  for (int ib = 1; ib <= Ncopy; ++ib) {
    h->SetBinContent(ib, hMeasOrig->GetBinContent(ib));
    h->SetBinError  (ib, hMeasOrig->GetBinError  (ib));
  }
  // bins au-delà de Ncopy restent à 0

  return h;
}

// ----------------------------------------------------------------------
// Fonction principale : écrit outPrefix.spe et outPrefix.rsp
// ----------------------------------------------------------------------
void WriteRW05Inputs(const char* responseFile  = "Response_PARIS50.root",
                     const char* hRespName     = "hResp",
                     const char* dataFile      = "resolution0508_bckgndsub_promptsumEvents_PARIS_corrected.root",
                     const char* hMeasName     = "sum_Res11keVGamma_PARIS50",
                     const char* outPrefix     = "rw05_input")
{
  // ---------- 1) Charger matrice de réponse ----------
  TFile* fResp = TFile::Open(responseFile, "READ");
  if (!fResp || fResp->IsZombie()) {
    std::cerr << "[RW05] ERROR: cannot open response file " << responseFile << std::endl;
    return;
  }

  TH2* hRespBase = dynamic_cast<TH2*>(fResp->Get(hRespName));
  if (!hRespBase) {
    std::cerr << "[RW05] ERROR: TH2 '" << hRespName
              << "' not found in " << responseFile << std::endl;
    fResp->Close();
    return;
  }

  // Forcer TH2D pour la sécurité
  TH2D* hResp = dynamic_cast<TH2D*>(hRespBase);
  if (!hResp) {
    // si ce n'est pas un TH2D, on clone en TH2D
    hResp = (TH2D*)hRespBase->Clone("hResp_as_TH2D");
  }

  int NX = hResp->GetNbinsX(); // measured bins
  int NY = hResp->GetNbinsY(); // true bins

  std::cout << "[RW05] Response matrix: NX (meas) = " << NX
            << ", NY (true) = " << NY << std::endl;

  // ---------- 2) Charger spectre mesuré ----------
  TFile* fData = TFile::Open(dataFile, "READ");
  if (!fData || fData->IsZombie()) {
    std::cerr << "[RW05] ERROR: cannot open data file " << dataFile << std::endl;
    fResp->Close();
    return;
  }

  TH1* hMeasOrig = dynamic_cast<TH1*>(fData->Get(hMeasName));
  if (!hMeasOrig) {
    std::cerr << "[RW05] ERROR: TH1 '" << hMeasName
              << "' not found in " << dataFile << std::endl;
    fResp->Close();
    fData->Close();
    return;
  }

  std::cout << "[RW05] Loaded measured histogram '" << hMeasName
            << "' with " << hMeasOrig->GetNbinsX() << " bins." << std::endl;

  // Construire un spectre mesuré compatible avec NX (comme dans RunUnfolding)
  TH1D* hMeas = MakeMeasuredChildForResponse(hMeasOrig, hResp, "hMeas_rw05");
  if (!hMeas) {
    std::cerr << "[RW05] ERROR: failed to build child measured histogram.\n";
    fResp->Close();
    fData->Close();
    return;
  }

  // ---------- 3) Ecriture du fichier .spe ----------
  // Format texte simple :
  //   ligne 1 : NX       (nombre de canaux)
  //   lignes suivantes : contenus par bin (0..NX-1), un entier (ou double) par ligne
  // RW05 lit ça en free-format, donc les retours à la ligne sont libres.

  TString speName = TString(outPrefix) + ".spe";
  std::ofstream ofsSpe(speName.Data());
  if (!ofsSpe) {
    std::cerr << "[RW05] ERROR: cannot open " << speName << " for writing.\n";
  } else {
    ofsSpe << "# RW05 spectrum file generated from ROOT\n";
    ofsSpe << NX << "\n";  // nombre de canaux (sans under/overflow)

    for (int i = 1; i <= NX; ++i) {
      double c = hMeas->GetBinContent(i);
      // RW05 préfère des entiers — on arrondit :
      long long ic = (long long) std::llround(c);
      ofsSpe << ic << "\n";
    }
    ofsSpe.close();
    std::cout << "[RW05] Wrote spectrum to " << speName << std::endl;
  }

  // ---------- 4) Ecriture du fichier .rsp ----------
  // Format texte simple :
  //   ligne 1 : NX NY
  //   ensuite NY lignes, une par "bins vrai" j (Y),
  //   chaque ligne contient NX valeurs : R_ij pour i=1..NX (measured)
  //
  // Ici on prend directement les bin contents de hResp (non normalisés),
  // RW05 gère lui-même ses normalisations internes.

  TString rspName = TString(outPrefix) + ".rsp";
  std::ofstream ofsRsp(rspName.Data());
  if (!ofsRsp) {
    std::cerr << "[RW05] ERROR: cannot open " << rspName << " for writing.\n";
  } else {
    ofsRsp << "# RW05 response matrix generated from ROOT\n";
    ofsRsp << NX << " " << NY << "\n";  // dimensions

    // Boucle sur les bins "true" (axe Y) :
    for (int j = 1; j <= NY; ++j) {
      // pour chaque énergie vraie j, on écrit la ligne de réponse en E_mesurée i
      for (int i = 1; i <= NX; ++i) {
        double v = hResp->GetBinContent(i, j);
        if (v < 0) v = 0.0;
        ofsRsp << v;
        if (i < NX) ofsRsp << " ";
      }
      ofsRsp << "\n";
    }
    ofsRsp.close();
    std::cout << "[RW05] Wrote response matrix to " << rspName << std::endl;
  }

  fResp->Close();
  fData->Close();

  std::cout << "[RW05] Done. Now you can run rw05 with:\n"
            << "  - spectrum file : " << speName << "\n"
            << "  - matrix file   : " << rspName << "\n"
            << "en te basant sur l'exemple 'demo' du dépôt.\n";
}