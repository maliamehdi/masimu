// CompareUnfoldWithReference.C
//
// Macro ROOT pour comparer :
//   - hUnfold_Bayes_perFission       (depuis UnfoldCompare_*.root)
//   - hUnfold_Direct                 (rescalé "per fission" ici)
//   - spectre de référence importé depuis un fichier texte/CSV
//     contenant 3 colonnes : Energy(MeV), Value, Uncertainty.
//
// Usage typique dans ROOT :
//   .L CompareUnfoldWithReference.C+
//   CompareUnfoldWithReference("UnfoldCompare_PARIS50.root",
//                              "hUnfold_Bayes_perFission",
//                              "hUnfold_Direct",
//                              2.882315e9,
//                              "Cf252_LaCl3_ref.csv",
//                              true,   // refInMeV
//                              "Cf252_Compare_PARIS50");
//
// ----------------------------------------------------------------------

#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

void CompareUnfoldWithReference(const char* unfoldFile   = "new_PARIS235_Bayes_SpectrumProperties.root",
                                const char* hBayesName   = "final_spectrum_per_fission_per_MeV",
                                const char* hDirectName  = "final_spectrum_per_fission_per_MeV",
                                double nFissionsDirect   = 1.926593664e+09,
                                const char* refFile      = "252Cf_PFG_ref.txt",
                                bool refEnergyInMeV      = true,
                                const char* outPrefix    = "Cf252_Compare_new_PARIS235")
{
  // ------------------ 1) Ouvrir le fichier ROOT d’unfolding ------------------
  TFile* fUnfold = TFile::Open(unfoldFile, "READ");
  if (!fUnfold || fUnfold->IsZombie()) {
    std::cerr << "[ERROR] Cannot open unfolding file: " << unfoldFile << std::endl;
    return;
  }

  // Histogramme Bayes (déjà "per fission" dans RunUnfolding)
  TH1D* hBayes = dynamic_cast<TH1D*>(fUnfold->Get(hBayesName));
  if (!hBayes) {
    std::cerr << "[WARN] '" << hBayesName
              << "' not found in " << unfoldFile
              << ". Trying fallback 'hUnfold_Bayes'..." << std::endl;
    hBayes = dynamic_cast<TH1D*>(fUnfold->Get("hUnfold_Bayes"));
  }
  if (!hBayes) {
    std::cerr << "[ERROR] No Bayesian unfolded histogram found." << std::endl;
    fUnfold->Close();
    return;
  }

  // Histogramme Direct (sera rescalé per fission ici)
  TH1D* hDirect = dynamic_cast<TH1D*>(fUnfold->Get(hDirectName));
  if (!hDirect) {
    std::cerr << "[ERROR] '" << hDirectName << "' not found in " << unfoldFile << std::endl;
    fUnfold->Close();
    return;
  }

  // Clones détachés du fichier
  TH1D* hBayesClone = (TH1D*)hBayes->Clone("hUnfold_Bayes_perFission_clone");
  hBayesClone->SetDirectory(nullptr);

  TH1D* hDirectPerFission = (TH1D*)hDirect->Clone("hUnfold_Direct_perFission");
  hDirectPerFission->SetDirectory(nullptr);
  if (nFissionsDirect > 0.0) {
    hDirectPerFission->Scale(1.0 / nFissionsDirect);
  } else {
    std::cerr << "[WARN] nFissionsDirect <= 0, no scaling applied to Direct spectrum." << std::endl;
  }

  std::cout << "[INFO] Loaded unfolded spectra from " << unfoldFile << std::endl;

  // ------------------ 2) Lire le fichier de référence (CSV / texte) ------------------
  std::ifstream fin(refFile);
  if (!fin.is_open()) {
    std::cerr << "[ERROR] Cannot open reference file: " << refFile << std::endl;
    fUnfold->Close();
    return;
  }

  std::vector<double> vx, vy, vex, vey;

  std::string line;
  bool firstLine = true;

  while (std::getline(fin, line)) {
    if (line.empty()) continue;

    // On saute la première ligne si c'est un header
    if (firstLine) {
      // Si tu veux forcer la lecture de la première ligne comme data,
      // commente la ligne suivante
      firstLine = false;
      // On peut inspecter la ligne pour voir si elle contient des lettres
      bool hasAlpha = false;
      for (char c : line) {
        if (std::isalpha((unsigned char)c)) { hasAlpha = true; break; }
      }
      if (hasAlpha) continue; // header textuel -> on le skip
    }

    std::istringstream iss(line);
    double E, val, err;
    if (!(iss >> E >> val >> err)) {
      // Si le séparateur est ';' ou ',', on peut faire un split manuel
      // mais ici on reste simple : lignes mal formées sont ignorées
      continue;
    }

    // Conversion MeV -> keV si nécessaire
    double EkeV = refEnergyInMeV ? (E * 1000.0) : E;

    vx.push_back(EkeV);
    vy.push_back(val);
    vex.push_back(0.0);   // pas d'incertitude en X
    vey.push_back(err);
  }
  fin.close();

  if (vx.empty()) {
    std::cerr << "[ERROR] No valid reference points read from " << refFile << std::endl;
    fUnfold->Close();
    return;
  }

  TGraphErrors* gRef = new TGraphErrors((int)vx.size(),
                                        vx.data(), vy.data(),
                                        vex.data(), vey.data());
  gRef->SetName("gRef_SEB347");
  gRef->SetTitle("Reference  SEB347 LaCl3 spectrum");

  gRef->SetMarkerStyle(20);
  gRef->SetMarkerSize(1.1);
  gRef->SetMarkerColor(kBlack);
  gRef->SetLineColor(kBlack);

  std::cout << "[INFO] Loaded " << vx.size()
            << " reference points from " << refFile << std::endl;

  // ------------------ 3) Ajuster styles des histos ------------------
  hBayesClone->SetLineColor(kBlue+1);
  hBayesClone->SetLineWidth(2);

  hDirectPerFission->SetLineColor(kRed+1);
  hDirectPerFission->SetLineWidth(2);
  hDirectPerFission->SetLineStyle(2);

  // On choisit une échelle Y adaptée (option logY pratique)
  gStyle->SetOptStat(0);

  TCanvas* c = new TCanvas("cCompare_Unfold_Ref",
                           "Unfolded vs Reference",
                           1000, 800);
  c->SetGrid();
  c->SetLogy(); // souvent utile pour les spectres en énergie

  // On commence par Bayes pour définir le cadre
  hBayesClone->SetTitle("Cf-252 reference PFG  spectrum;E_{true} [keV];(arb. units or per fission)");
  hBayesClone->Draw("HIST");

  // Pour être sûr que le graphe soit visible
  double ymin = hBayesClone->GetMinimum(1e-30); // ignore bins vides
  double ymax = hBayesClone->GetMaximum();
  if (hDirectPerFission->GetMaximum() > ymax) ymax = hDirectPerFission->GetMaximum();
  // On peut aussi regarder la valeur min positive
  if (ymin <= 0) ymin = ymax*1e-6;
  hBayesClone->SetMinimum(ymin*0.5);
  hBayesClone->SetMaximum(ymax*5.0);

  // On dessine la méthode directe
  hDirectPerFission->Draw("HIST SAME");

  // Puis les points expérimentaux
  gRef->Draw("P SAME");

  // Légende
  TLegend* leg = new TLegend(0.55, 0.68, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hBayesClone,       "Bayesian unfolded (per fission)", "l");
  leg->AddEntry(hDirectPerFission, "Direct unfolded (per fission)",   "l");
  leg->AddEntry(gRef,              "Reference data (LaCl_{3}:Ce)",    "pe");
  leg->Draw();

  // ------------------ 4) Sauvegarde ------------------
  TString base = outPrefix;
  c->SaveAs(base + ".png");
  c->SaveAs(base + ".pdf");

  std::cout << "[INFO] Comparison plot saved as: "
            << base << ".png/.pdf" << std::endl;

  // Optionnel : sauver aussi dans un ROOT
  TFile* fOut = TFile::Open((base + ".root").Data(), "RECREATE");
  if (fOut && !fOut->IsZombie()) {
    fOut->cd();
    hBayesClone->Write("hUnfold_Bayes_perFission");
    hDirectPerFission->Write("hUnfold_Direct_perFission");
    gRef->Write("gRef_LaCl3");
    fOut->Close();
    std::cout << "[INFO] ROOT file written: " << base << ".root" << std::endl;
  }

  fUnfold->Close();
}