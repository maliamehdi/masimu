// PlotResponseMatrix.C
// Macro ROOT pour tracer la matrice de réponse Etrue vs Emeas
// Source: ntuple "resp" avec colonnes: eventID, parisIndex, Etrue_keV, Emeas_keV, EdepCe_keV, EdepNaI_keV
//
// Axes (version actuelle) :
//   X = Emeas_keV (énergie mesurée)
//   Y = Etrue_keV (énergie vraie)
//
// Usage typique :
// root -l PlotResponseMatrix.C'+
//   PlotResponseMatrix("out_PARIS50_100000.root",
//                      -1, 400, 0, 15000,
//                      400, 0, 15000,
//                      false, "PARIS50_resp_resbin",
//                      true,  true,  false,
//                      true, "PARIS50");'
//
// Paramètres :
//   fname        : fichier ROOT contenant l'ntuple resp
//   parisIndex   : -1 = tous les détecteurs, sinon filtre sur un index PARIS
//   nbinsTrue    : nbins axe énergie vraie (Etrue)  [Y]
//   trueMin/Max  : bornes Etrue (keV)
//   nbinsMeas    : nbins axe énergie mesurée (Emeas) [X]
//   measMin/Max  : bornes Emeas (keV)
//   normRows     : si true, normalise chaque ligne en Y pour somme=1 (PDF(Etrue|Emeas))
//   savePrefix   : préfixe pour les fichiers de sortie
//   useDraw      : true => TTree::Draw, false => boucle explicite
//   logz         : échelle log en Z
//   smooth       : lissage léger
//   resolutionbin: si true, binning variable sur X et Y basé sur (resA,respower)
//   detName      : nom du détecteur PARIS ("PARIS50", etc.) pour récupérer (resA,respower)

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TVectorD.h>
#include <iostream>
#include <sstream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>

struct ResParams {
  double A;
  double power;
};

// Paramètres de résolution PARIS (A, power)
static const std::map<std::string, ResParams> parisRes = { // run du 05/08/2024
  {"PARIS50",  {1.12145,  -0.441244}},
  {"PARIS70",  {1.80973,  -0.550685}},
  {"PARIS90",  {1.94868,  -0.564616}},
  {"PARIS110", {2.11922,  -0.582147}},
  {"PARIS130", {0.794233, -0.377311}},
  {"PARIS235", {1.30727,  -0.477402}},
  {"PARIS262", {1.76345,  -0.542769}},
  {"PARIS278", {1.98579,  -0.559095}},
  {"PARIS305", {1.9886,   -0.574021}}
};

void PlotResponseMatrix(const char* fname = "PARIS50/output_PARIS50_12million.root",
                        int  parisIndex   = 1,
                        int  nbinsTrue    = 155,
                        double trueMin    = 0.0,
                        double trueMax    = 15000.0,
                        int  nbinsMeas    = 155,
                        double measMin    = 0.0,
                        double measMax    = 15000.0,
                        bool normRows     = false,
                        const char* savePrefix = "plot_Response12M",
                        bool useDraw      = true,
                        bool logz         = false,
                        bool smooth       = false,
                        bool resolutionbin = true,
                        const char* detName = "PARIS50")
{
  // --- Ouverture du fichier et de l'arbre ---
  TFile* f = TFile::Open(fname, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "[ERROR] Impossible d'ouvrir le fichier " << fname << std::endl;
    return;
  }

  TTree* t = dynamic_cast<TTree*>(f->Get("resp"));
  if (!t) {
    std::cerr << "[ERROR] Ntuple 'resp' introuvable dans " << fname << std::endl;
    f->Close();
    return;
  }

  std::cout << "[INFO] Opened file " << fname << std::endl;
  std::cout << "[INFO] Entries in resp = " << t->GetEntries() << std::endl;

  // --- Construction du binning variable (si demandé) ---
  std::vector<double> trueEdges;
  std::vector<double> measEdges;

  if (resolutionbin) {
    std::string key = detName ? std::string(detName) : std::string();
    auto it = parisRes.find(key);
    if (it == parisRes.end()) {
      std::cerr << "[WARN] resolutionbin=true mais detName='" << key
                << "' introuvable dans la map, binning uniforme conserve." << std::endl;
    } else {
      double resA = it->second.A;
      double respower = it->second.power;
      std::cout << "[INFO] Using resolution binning for " << key
                << " : resA=" << resA << " respower=" << respower << std::endl;

      if (resA != 0.0 && resA < 100.0 && respower != 0.0 && respower < 1.0) {
        auto buildEdges = [&](int nbins, double Emin) {
          std::vector<double> edges(nbins+1);
          edges[0] = Emin;
          edges[1] = Emin + 11.0;  // bin initial ~11 keV
          for (int i = 2; i < nbins+1; ++i) {
            double prev = edges[i-1];
            double step = resA * TMath::Power(prev, respower) * prev;
            edges[i] = prev + step;
          }
          return edges;
        };

        trueEdges = buildEdges(nbinsTrue, trueMin);
        measEdges = buildEdges(nbinsMeas, measMin);

        std::cout << "[DEBUG] First trueEdges (keV): ";
        for (int i = 0; i < std::min(10, nbinsTrue+1); ++i)
          std::cout << trueEdges[i] << " ";
        std::cout << "\n[INFO] Last true edge = " << trueEdges.back() << " keV" << std::endl;

        std::cout << "[DEBUG] First measEdges (keV): ";
        for (int i = 0; i < std::min(10, nbinsMeas+1); ++i)
          std::cout << measEdges[i] << " ";
        std::cout << "\n[INFO] Last meas edge = " << measEdges.back() << " keV" << std::endl;

      } else {
        std::cerr << "[WARN] resA/respower hors plage raisonnable, binning uniforme conserve." << std::endl;
        trueEdges.clear();
        measEdges.clear();
      }
    }
  }

  // --- Construction de l'histogramme ---
  std::ostringstream hname; hname << "hResp_Etrue_vs_Emeas";
  TH2D* h = nullptr;

  if (resolutionbin && !trueEdges.empty() && !measEdges.empty()) {
    // Binning variable sur les deux axes
    // X = Emeas, Y = Etrue
    h = new TH2D(hname.str().c_str(),
                 "Response matrix;E_{mes} [keV];E_{true} [keV]",
                 nbinsMeas, measEdges.data(),
                 nbinsTrue, trueEdges.data());
  } else {
    // Binning uniforme (comportement original)
    h = new TH2D(hname.str().c_str(),
                 "Response matrix;E_{mes} [keV];E_{true} [keV]",
                 nbinsMeas, measMin, measMax,   // X
                 nbinsTrue, trueMin, trueMax);  // Y
  }

  // --- Sélection éventuelle sur parisIndex ---
  std::string sel;
  if (parisIndex >= 0) {
    std::ostringstream ss; ss << "parisIndex==" << parisIndex;
    sel = ss.str();
  }

  // --- Remplissage ---
  if (useDraw) {
    // Attention : syntaxe Y:X >> hist
    // Ici, Y = Etrue, X = Emeas
    std::ostringstream drawCmd;
    drawCmd << "Etrue_keV:Emeas_keV>>" << h->GetName();
    t->Draw(drawCmd.str().c_str(), sel.empty() ? nullptr : sel.c_str(), "goff");
    std::cout << "[INFO] Filled via TTree::Draw with command: "
              << drawCmd.str() << std::endl;
  } else {
    // Boucle explicite
    double Etrue=0, Emeas=0; int pIdx=0;
    t->SetBranchAddress("Etrue_keV", &Etrue);
    t->SetBranchAddress("Emeas_keV", &Emeas);
    t->SetBranchAddress("parisIndex", &pIdx);
    Long64_t n = t->GetEntries();
    for (Long64_t i=0;i<n;++i) {
      t->GetEntry(i);
      if (parisIndex >= 0 && parisIndex != pIdx) continue;
      // X = Emeas, Y = Etrue
      h->Fill(Emeas, Etrue);
    }
    std::cout << "[INFO] Filled via explicit loop, n entries = " << n << std::endl;
  }

  std::cout << "[INFO] Integral(h) after filling = " << h->Integral() << std::endl;

  // Clone brut avant normalisation
  TH2D* hRaw = (TH2D*)h->Clone("ResponseMatrix_raw");
  hRaw->SetDirectory(nullptr);

  // Maintenant qu'il est rempli, on peut détacher h aussi
  h->SetDirectory(nullptr);

  // --- Normalisation par "ligne" (ici : à Y fixe, on normalise sur X) ---
  // NB : comme on a inversé les axes, il faut décider :
  //   - si tu veux PDF(Emeas | Etrue)  -> normalise sur X à Y fixé
  // Ici on garde cette logique : chaque Y (Etrue) => somme sur tous les X.
  if (normRows) {
    for (int iy=1; iy<=h->GetNbinsY(); ++iy) {
      double rowSum = 0.0;
      for (int ix=1; ix<=h->GetNbinsX(); ++ix) rowSum += h->GetBinContent(ix, iy);
      if (rowSum > 0) {
        for (int ix=1; ix<=h->GetNbinsX(); ++ix) {
          double v = h->GetBinContent(ix, iy);
          h->SetBinContent(ix, iy, v/rowSum);
        }
      }
    }
    h->SetName("ResponseMatrix"); // nom propre pour la version normalisée
    h->SetTitle("Response matrix (row-normalized);E_{mes} [keV];E_{true} [keV]");
  } else {
    h->SetName("ResponseMatrix"); // version brute mais avec nom explicite
  }

  if (smooth) {
    h->Smooth(1);
  }

  // --- Style & affichage ---
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);

  TCanvas* c = new TCanvas("cResp", "Response Matrix", 1000, 800);
  if (logz) c->SetLogz();
  h->Draw("COLZ");

  // Info pave
  TPaveText* pt = new TPaveText(0.13,0.82,0.42,0.93, "NDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  pt->AddText(Form("File: %s", fname));
  if (parisIndex >= 0) pt->AddText(Form("parisIndex = %d", parisIndex));
  if (detName && detName[0] != '\0') pt->AddText(Form("Det: %s", detName));
  pt->AddText(Form("Entries (filled): %.0f", h->GetEntries()));
  if (normRows) pt->AddText("Row norm: yes");
  if (resolutionbin && !trueEdges.empty() && !measEdges.empty())
    pt->AddText("Binning: resolution (X & Y)");
  pt->Draw();

  // Fichiers image
  std::string base = savePrefix;
  c->SaveAs((base+".png").c_str());
  c->SaveAs((base+".pdf").c_str());

  // --- Sauvegarde ROOT ---
  {
    std::string rootOut = base + ".root";
    TFile fout(rootOut.c_str(), "RECREATE");

    // matrice brute
    if (hRaw) hRaw->Write();        // nom = "ResponseMatrix_raw"
    // matrice éventuellement normalisée
    if (h)    h->Write();           // nom = "ResponseMatrix"

    // Sauvegarde des edges si binning résolution
    if (resolutionbin && !trueEdges.empty() && !measEdges.empty()) {
      TVectorD vTrue(trueEdges.size());
      TVectorD vMeas(measEdges.size());
      for (size_t i=0;i<trueEdges.size(); ++i) vTrue[i]=trueEdges[i];
      for (size_t i=0;i<measEdges.size();  ++i) vMeas[i]=measEdges[i];
      vTrue.Write("TrueBinEdges");
      vMeas.Write("MeasBinEdges");
    }

    fout.Close();
    std::cout << "[INFO] ROOT output saved to: " << rootOut << std::endl;
  }

  std::cout << "[INFO] Saved: " << base << ".png/.pdf/.root" << std::endl;

  // Nettoyage fichier d'entrée
  f->Close();
}