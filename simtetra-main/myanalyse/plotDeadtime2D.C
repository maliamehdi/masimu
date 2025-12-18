// plotDeadtime2D.C
// Construit un TH2D : X = multiplicité nominale (multN), Y = multiplicité mesurée (sum rings deadtime)
// en lisant les *_hists.root générés par plotTetraFromList().

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TString.h>
#include <TStyle.h>
#include <TLatex.h>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cctype>

static std::string trim_copy(std::string s) {
  auto notspace = [](int ch){ return !std::isspace(ch); };
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), notspace));
  s.erase(std::find_if(s.rbegin(), s.rend(), notspace).base(), s.end());
  return s;
}

// Essaie d'extraire la multiplicité nominale depuis le chemin
// ex: ".../mult12/..." -> 12, ou "...fixed_mult_07..." -> 7
static int extractNominalMult(const std::string& path) {
  // 1) pattern /multN/
  {
    const std::string key = "/mult";
    auto pos = path.find(key);
    if (pos != std::string::npos) {
      pos += key.size();
      int val = 0;
      bool ok = false;
      while (pos < path.size() && std::isdigit(path[pos])) {
        ok = true;
        val = 10*val + (path[pos]-'0');
        ++pos;
      }
      if (ok) return val;
    }
  }
  // 2) pattern fixed_mult_XX
  {
    const std::string key = "fixed_mult_";
    auto pos = path.find(key);
    if (pos != std::string::npos) {
      pos += key.size();
      int val = 0;
      bool ok = false;
      while (pos < path.size() && std::isdigit(path[pos])) {
        ok = true;
        val = 10*val + (path[pos]-'0');
        ++pos;
      }
      if (ok) return val;
    }
  }
  return -1;
}

void plotDeadtime2D(const char* listfile = "filelist.dat")
{
  gStyle->SetOptStat(0);

  // --- Lire la liste ---
  std::ifstream fin(listfile);
  if (!fin.is_open()) {
    std::cerr << "Cannot open list file: " << listfile << std::endl;
    return;
  }

  struct Item { std::string infile; double tGate=150000; double budget=14; int multNom=-1; std::string hfile; };
  std::vector<Item> items;

  std::string line;
  int lineNo = 0;
  while (std::getline(fin, line)) {
    ++lineNo;
    line = trim_copy(line);
    if (line.empty()) continue;
    if (line[0]=='#') continue;
    if (line.size()>1 && line[0]=='/' && line[1]=='/') continue;

    std::istringstream iss(line);
    Item it;
    if (!(iss >> it.infile)) continue;
    if (!(iss >> it.tGate))   it.tGate = 150000.0;
    if (!(iss >> it.budget))  it.budget = 14.0;

    it.multNom = extractNominalMult(it.infile);
    if (it.multNom < 0) {
      std::cerr << "WARN line " << lineNo << ": cannot infer nominal multiplicity from path: "
                << it.infile << "\n"
                << " -> will skip this file.\n";
      continue;
    }

    // Reconstruire le nom du *_hists.root produit par plotTetra()
    // plotTetra écrit: outBase = "<dir>/plot_<base>" puis rootOut = "<outBase>_hists.root"
    TString tin(it.infile.c_str());
    TString dir  = gSystem->DirName(tin);
    TString base = gSystem->BaseName(tin);
    if (base.EndsWith(".root")) base.ReplaceAll(".root", "");
    TString outBase = Form("%s/plot_%s", dir.Data(), base.Data());
    it.hfile = std::string(Form("%s_hists.root", outBase.Data()));

    items.push_back(it);
  }

  if (items.empty()) {
    std::cerr << "No valid entries found in " << listfile << std::endl;
    return;
  }

  // Trier par multiplicité nominale
  std::sort(items.begin(), items.end(), [](const Item& a, const Item& b){
    return a.multNom < b.multNom;
  });

  // --- Déterminer Y max (multiplicité mesurée) en scannant les histos ---
  int yMax = 0;
  for (const auto& it : items) {
    std::unique_ptr<TFile> f(TFile::Open(it.hfile.c_str(), "READ"));
    if (!f || f->IsZombie()) {
      std::cerr << "WARN cannot open: " << it.hfile << " (did you run plotTetraFromList ?) \n";
      continue;
    }
    auto* h = dynamic_cast<TH1*>(f->Get("hSumRings_deadtime"));
    if (!h) {
      std::cerr << "WARN missing hSumRings_deadtime in " << it.hfile << "\n";
      continue;
    }
    yMax = std::max(yMax, (int)std::ceil(h->GetXaxis()->GetXmax()));
  }
  if (yMax <= 0) yMax = 20;

  // --- Construire TH2D : X = nominal multiplicity, Y = measured(deadtime) ---
  int xMin = items.front().multNom;
  int xMax = items.back().multNom;

  auto* h2 = new TH2D("h2_deadtimeResponse",
                      "Deadtime response: sum of rings;Nominal multiplicity;Measured multiplicity (sum rings, deadtime);Events",
                      (xMax - xMin + 1), xMin - 0.5, xMax + 0.5,
                      yMax + 1, -0.5, yMax + 0.5);

  // --- Remplir h2 depuis chaque hSumRings_deadtime ---
  for (const auto& it : items) {
    std::unique_ptr<TFile> f(TFile::Open(it.hfile.c_str(), "READ"));
    if (!f || f->IsZombie()) continue;
    auto* h = dynamic_cast<TH1*>(f->Get("hSumRings_deadtime"));
    if (!h) continue;

    const int nb = h->GetNbinsX();
    for (int b = 1; b <= nb; ++b) {
      double y = h->GetXaxis()->GetBinCenter(b);
      double w = h->GetBinContent(b);
      if (w <= 0) continue;
      h2->Fill((double)it.multNom, y, w);
    }
  }

  // --- Plot ---
  TCanvas* c = new TCanvas("c_deadtime2D", "Deadtime multiplicity 2D", 1100, 750);
  c->SetRightMargin(0.14);
  c->SetLeftMargin(0.12);
  c->SetBottomMargin(0.12);
  h2->Draw("COLZ");

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.035);
  lat.DrawLatex(0.15, 0.94, Form("Built from %s", listfile));

  // Sauvegardes
  TString outPng = TString::Format("deadtime2D_%s.png", gSystem->BaseName(listfile));
  outPng.ReplaceAll(".txt","");
  c->SaveAs(outPng);

  TString outRoot = TString::Format("deadtime2D_%s.root", gSystem->BaseName(listfile));
  outRoot.ReplaceAll(".txt","");
  TFile fout(outRoot, "RECREATE");
  h2->Write();
  fout.Close();

  std::cout << "Saved: " << outPng << "\n";
  std::cout << "Saved: " << outRoot << " (contains h2_deadtimeResponse)\n";
}