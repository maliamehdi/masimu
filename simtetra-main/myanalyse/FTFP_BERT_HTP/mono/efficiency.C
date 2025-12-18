// efficiency.C - Version corrigée
// ROOT macro: parse summary .txt files and plot efficiency vs neutron energy
// Usage (from myanalyse):
//   root -l
//   .L efficiency.C+
//   plot_efficiency("FTFP_BERT_HTP/mono", "FTFP_BERT_HTP/mono/mono_eff_scan")

#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TLegend.h>

#include <fstream>
#include <iostream>
#include <regex>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>

// ----------------------------------------------------------------------
// Structure pour stocker un point (un fichier summary)
// ----------------------------------------------------------------------
struct Point {
  double E_MeV       = -1.0;   // énergie neutron (MeV)
  double eff_percent = -1.0;   // efficacité (%) 
  double eff_unc_percent = 0.0; // incertitude sur l’efficacité (%)
  double meanMult    = -1.0;   // <multiplicité mesurée>
  long long totalEvents = 0;
  long long eventsAllDet = 0;
  std::string file;            // nom du fichier summary
};

// ----------------------------------------------------------------------
// Helper: parse un token d'énergie comme "500keV" ou "10" -> MeV
// ----------------------------------------------------------------------
static double parseEnergyToken(const std::string &tok) {
  std::smatch mk;
  std::regex kevrg("^([0-9]+)keV$");
  std::regex mevrg("^([0-9]+)$");

  if (std::regex_match(tok, mk, kevrg)) {
    return std::stod(mk[1]) / 1000.0; // keV -> MeV
  }
  if (std::regex_match(tok, mk, mevrg)) {
    return std::stod(mk[1]);          // déjà en MeV
  }
  return -1.0;
}

// ----------------------------------------------------------------------
// Parse l'énergie à partir du nom du fichier summary
//   Accepte par ex. :
//   - output_mono_n_500keV_Tetra_summary.txt
//   - output_mono_n_500keV_summary.txt
//   - plot_output_mono_n_500keV_Tetra_summary.txt
//   - plot_output_mono_n_500keV_summary.txt
// ----------------------------------------------------------------------
static double parseEnergyFromName(const std::string &name) {
  std::regex rg("^(?:plot_)?output_mono_n_([^_]+)(?:_Tetra)?_summary\\.txt$");
  std::smatch m;
  if (!std::regex_match(name, m, rg)) return -1.0;

  std::string tok = m[1]; // ex: "500keV" ou "10"
  return parseEnergyToken(tok);
}

// ----------------------------------------------------------------------
// Parse un fichier *_summary.txt et remplit un Point
//   Exemple de format :
//
//   Input: FTFP_BERT_HTP/mono/output_mono_n_10.root
//   Total events: 100000
//   ...
//   Events (all detections <= gate): 8408 (eff = 8.408000 +/- 0.087756 %)
//   ...
//   Mean measured multiplicity in gate: 0.084080
// ----------------------------------------------------------------------
static Point parseSummaryFile(const std::string &path) {
  Point p;
  p.file = path;

  // Récupérer juste le nom de fichier (sans le path)
  size_t pos = path.find_last_of('/');
  std::string fname = (pos == std::string::npos) ? path : path.substr(pos+1);
  // 1ère tentative : énergie à partir du nom du summary
  p.E_MeV = parseEnergyFromName(fname);

  std::ifstream in(path);
  if (!in.is_open()) {
    std::cerr << "Cannot open " << path << "\n";
    return p;
  }

  std::string line;

  // regex pour les différentes lignes

  // Input: FTFP_BERT_HTP/mono/output_mono_n_10.root
  std::regex re_input(
    "^Input:\\s+.*output_mono_n_([^./]+)\\.root\\s*$"
  );

  // Mean measured multiplicity
  std::regex re_mean("Mean measured multiplicity.*:\\s*([0-9eE+\\.-]+)");

  // Total events: N
  std::regex re_tot("^Total events:\\s*([0-9]+)");

  // Events (all detections <= gate): N (eff = X +/- dX %)
  std::regex re_all_with_err(
      "^Events \\(all detections .*\\):\\s*([0-9]+) "
      "\\(eff = ([0-9eE+\\.-]+) \\+/- ([0-9eE+\\.-]+) %\\)");

  // Events (all detections <= gate): N (eff = X %)
  std::regex re_all_no_err(
      "^Events \\(all detections .*\\):\\s*([0-9]+) "
      "\\(eff = ([0-9eE+\\.-]+) %\\)");

  std::smatch m;

  double eff_from_line  = -1.0;
  double err_from_line  = -1.0;

  while (std::getline(in, line)) {

    // Si l’énergie n’a pas encore été trouvée,
    // tente de la lire dans la ligne Input:
    if (p.E_MeV <= 0.0 && std::regex_search(line, m, re_input)) {
      std::string tok = m[1]; // "10", "500keV", etc.
      p.E_MeV = parseEnergyToken(tok);
      continue;
    }

    // Mean multiplicity
    if (std::regex_search(line, m, re_mean)) {
      p.meanMult = std::stod(m[1]);
      continue;
    }

    // Total events
    if (std::regex_search(line, m, re_tot)) {
      p.totalEvents = std::stoll(m[1]);
      continue;
    }

    // Events (all detections ...) avec erreur
    if (std::regex_search(line, m, re_all_with_err)) {
      p.eventsAllDet = std::stoll(m[1]);
      eff_from_line  = std::stod(m[2]); // efficacité %
      err_from_line  = std::stod(m[3]); // erreur absolue en points de %
      continue;
    }

    // Events (all detections ...) sans erreur
    if (std::regex_search(line, m, re_all_no_err)) {
      p.eventsAllDet = std::stoll(m[1]);
      eff_from_line  = std::stod(m[2]);
      continue;
    }
  }

  // Déterminer l’efficacité
  if (eff_from_line >= 0.0) {
    p.eff_percent = eff_from_line;
  } else if (p.totalEvents > 0) {
    // fallback si la ligne n'était pas trouvée : ratio simple
    p.eff_percent = 100.0 * double(p.eventsAllDet) / double(p.totalEvents);
  } else {
    p.eff_percent = -1.0;
  }

  // Déterminer l’incertitude
  if (err_from_line >= 0.0) {
    // si l’erreur est écrite dans le summary, on la prend telle quelle
    p.eff_unc_percent = err_from_line;
  } else {
    // sinon, on recalcule une erreur statistique
    if (p.meanMult >= 0.0 && p.totalEvents > 0) {
      p.eff_unc_percent = 100.0 * std::sqrt(p.meanMult / double(p.totalEvents));
    }
    else if (p.eventsAllDet > 0 && p.totalEvents > 0) {
      double frac = double(p.eventsAllDet) / double(p.totalEvents);
      p.eff_unc_percent = 100.0 * std::sqrt(frac * (1.0 - frac) / double(p.totalEvents));
    } else {
      p.eff_unc_percent = 0.0;
    }
  }

  return p;
}

// ----------------------------------------------------------------------
// Fonction principale de tracé
//   dir      : répertoire contenant les *_summary.txt
//   outprefix: préfixe de sortie pour les fichiers PNG / ROOT / TXT
// ----------------------------------------------------------------------
void plot_efficiency(const char* dir = ".",
                     const char* outprefix = "./mono_eff_scan")
{
  // --------- Liste les fichiers dans le répertoire ---------
  void *dp = gSystem->OpenDirectory(dir);
  if (!dp) {
    std::cerr << "Cannot open dir: " << dir << "\n";
    return;
  }

  const char* ent;
  std::vector<std::string> files;

  while ((ent = gSystem->GetDirEntry(dp))) {
    std::string s(ent);

    // On accepte :
    //  - plot_output_mono_n_*_summary.txt
    //  - output_mono_n_*_summary.txt (au cas où)
    bool isSummary = (s.size() > 12 && s.rfind("_summary.txt") == s.size() - 12);
    bool isMonoPrefix =
        (s.find("plot_output_mono_n_") == 0) ||
        (s.find("output_mono_n_") == 0);

    if (isSummary && isMonoPrefix) {
      files.push_back(std::string(dir) + "/" + s);
    }
  }
  gSystem->FreeDirectory(dp);

  if (files.empty()) {
    std::cerr << "No summary files found in " << dir << "\n";
    return;
  }

  // --------- Parse tous les fichiers ---------
  std::vector<Point> pts;
  pts.reserve(files.size());

  for (auto &f : files) {
    Point p = parseSummaryFile(f);
    if (p.E_MeV > 0 && p.eff_percent >= 0.0) {
      pts.push_back(p);
    } else if (p.eff_percent >= 0.0 && p.E_MeV <= 0) {
      std::cerr << "Skipping (no energy parsed): " << f << "\n";
    } else {
      std::cerr << "Skipping (no valid eff): " << f << "\n";
    }
  }

  if (pts.empty()) {
    std::cerr << "No valid points parsed.\n";
    return;
  }

  // --------- Trie par énergie croissante ---------
  std::sort(pts.begin(), pts.end(),
            [](const Point &a, const Point &b) {
              return a.E_MeV < b.E_MeV;
            });

  // --------- Crée le TGraphErrors ---------
  int n = (int)pts.size();
  TGraphErrors *g = new TGraphErrors(n);

  for (int i = 0; i < n; ++i) {
    g->SetPoint(i, pts[i].E_MeV, pts[i].eff_percent);
    g->SetPointError(i, 0.0, pts[i].eff_unc_percent); // barres d'erreur en Y
  }

  gStyle->SetOptStat(0);

  TCanvas *c = new TCanvas("cEff","Efficiency vs neutron energy",900,600);
  c->SetGridx();
  c->SetGridy();
  c->SetLogx(); // axe énergie en log

  g->SetTitle("Detection efficiency;Neutron energy (MeV);Efficiency (%)");
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.1);
  g->SetLineWidth(2);

  // Axe X : plage un peu élargie
  double xmin = pts.front().E_MeV * 0.8;
  if (xmin <= 0) xmin = pts.front().E_MeV * 0.5;
  if (xmin <= 0) xmin = 1e-4;
  double xmax = pts.back().E_MeV * 1.2;

  g->GetXaxis()->SetMoreLogLabels();
  g->GetXaxis()->SetNoExponent();
  g->GetXaxis()->SetRangeUser(xmin, xmax);

  g->Draw("APLE");

  TLegend *leg = new TLegend(0.15,0.75,0.45,0.88);
  leg->AddEntry(g,"Detection efficiency","p");
  leg->AddEntry((TObject*)nullptr,Form("N points = %d", n),"");
  leg->Draw();

  // --------- Sauvegardes ---------
  // PNG
  std::string png = std::string(outprefix) + "_eff_vs_energy.png";
  c->SaveAs(png.c_str());

  // ROOT
  std::string rootf = std::string(outprefix) + "_eff_vs_energy.root";
  TFile fout(rootf.c_str(), "RECREATE");
  g->Write("eff_vs_energy");

  // Tableau dans un TTree
  TTree t("EffTable","Efficiency points");

  double e, eff, eff_unc;
  long long tot, det;
  std::string fname;

  t.Branch("E_MeV",            &e);
  t.Branch("eff_percent",      &eff);
  t.Branch("eff_unc_percent",  &eff_unc);
  t.Branch("totalEvents",      &tot);
  t.Branch("eventsAllDet",     &det);
  t.Branch("filename",         &fname);

  for (auto &p : pts) {
    e       = p.E_MeV;
    eff     = p.eff_percent;
    eff_unc = p.eff_unc_percent;
    tot     = p.totalEvents;
    det     = p.eventsAllDet;
    fname   = p.file;
    t.Fill();
  }
  t.Write();
  fout.Close();

  // TXT
  std::string txt = std::string(outprefix) + "_eff_vs_energy.txt";
  std::ofstream ofs(txt);
  ofs << "#E_MeV\teff_percent\teff_unc_percent\tTotalEvents\tEventsAllDet\tFile\n";
  for (auto &p : pts) {
    ofs << p.E_MeV << "\t"
        << p.eff_percent << "\t"
        << p.eff_unc_percent << "\t"
        << p.totalEvents << "\t"
        << p.eventsAllDet << "\t"
        << p.file << "\n";
  }
  ofs.close();

  std::cout << "Wrote:\n"
            << "  " << png << "\n"
            << "  " << rootf << "\n"
            << "  " << txt << "\n";
}