#include <TGraph.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TPaveText.h>
#include <fstream>
#include <iostream>
#include <regex>
#include <vector>
#include <algorithm>
#include <cctype>

struct EffPoint {
  double E_MeV;
  double eff_percent;
  long long totalEvents;
  long long eventsAllDet;
  std::string filename;
};

double ParseEnergyMeV(const std::string& fname) {
  // fname pattern: output_mono_n_<token>_Tetra_summary.txt
  // token examples: 100keV, 500keV, 6, 14, 3, etc.
  std::regex rg("output_mono_n_([^_]+)_Tetra_summary\\.txt");
  std::smatch m;
  if (!std::regex_search(fname, m, rg)) return -1.0;
  std::string token = m[1];
  // keV case
  std::regex kevrg("([0-9]+)keV");
  std::smatch mk;
  if (std::regex_match(token, mk, kevrg)) {
    double kev = std::stod(mk[1]);
    return kev/1000.0;
  }
  // MeV integer
  std::regex mevrg("([0-9]+)");
  std::smatch mm;
  if (std::regex_match(token, mm, mevrg)) {
    return std::stod(mm[1]); // already MeV
  }
  return -1.0;
}

EffPoint ParseFile(const std::string& path) {
  EffPoint p;
  p.filename = path;
  p.E_MeV = ParseEnergyMeV(path);
  p.eff_percent = -1.0;
  p.totalEvents = 0;
  p.eventsAllDet = 0;

  std::ifstream in(path);
  if (!in) {
    std::cerr << "Cannot open " << path << std::endl;
    return p;
  }
  std::string line;
  while (std::getline(in, line)) {
    // Total events
    {
      std::regex rg("^Total events:\\s*([0-9]+)");
      std::smatch m;
      if (std::regex_search(line, m, rg)) {
        p.totalEvents = std::stoll(m[1]);
        continue;
      }
    }
    // Events (all detections <= gate): N (eff = X %)
    {
      std::regex rg("^Events \\(all detections .*\\):\\s*([0-9]+) \\(eff = ([0-9eE+\\.-]+) %\\)");
      std::smatch m;
      if (std::regex_search(line, m, rg)) {
        p.eventsAllDet = std::stoll(m[1]);
        p.eff_percent = std::stod(m[2]);
        continue;
      }
    }
  }
  // Fallback compute if eff not parsed
  if (p.eff_percent < 0.0 && p.totalEvents > 0)
    p.eff_percent = 100.0 * double(p.eventsAllDet) / double(p.totalEvents);
  return p;
}

void plot_mono_efficiency(const char* dir = ".", const char* outprefix = "mono_eff") {
  // Collect summary files
  std::vector<std::string> files;
  void* dirp = gSystem->OpenDirectory(dir);
  if (!dirp) {
    std::cerr << "Cannot open directory: " << dir << std::endl;
    return;
  }
  const char* fname;
  while ((fname = gSystem->GetDirEntry(dirp))) {
    std::string s = fname;
    if (s.find("output_mono_n_") == 0 && s.find("_Tetra_summary.txt") != std::string::npos) {
      files.push_back(std::string(dir) + "/" + s);
    }
  }
  gSystem->FreeDirectory(dirp);
  if (files.empty()) {
    std::cerr << "No summary files found in " << dir << std::endl;
    return;
  }

  std::vector<EffPoint> points;
  points.reserve(files.size());
  for (auto& f : files) {
    auto p = ParseFile(f);
    if (p.E_MeV > 0) points.push_back(p);
  }
  if (points.empty()) {
    std::cerr << "No valid energy points parsed." << std::endl;
    return;
  }

  std::sort(points.begin(), points.end(),
            [](const EffPoint& a, const EffPoint& b){ return a.E_MeV < b.E_MeV; });

  // Prepare graph
  auto* g = new TGraph(points.size());
  for (size_t i=0;i<points.size();++i) {
    g->SetPoint(i, points[i].E_MeV, points[i].eff_percent);
  }

  gStyle->SetOptStat(0);
  auto* c = new TCanvas("cEff","Efficiency vs neutron energy",900,600);
  c->SetGrid();
c->SetLogx(); // semi-log scale on X
  g->SetTitle("Detection efficiency (all detections within gate);Neutron energy (MeV);Efficiency (%)");
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.1);
  g->SetLineWidth(2);
  g->Draw("APL");

  // Annotate low-energy points if desired
  auto* leg = new TLegend(0.12,0.72,0.42,0.88);
  leg->AddEntry(g,"Efficiency (all detections)","");
  leg->AddEntry((TObject*)nullptr,Form("Points: %zu", points.size()),"");
  leg->Draw();

  // Save PNG
  std::string png = std::string(outprefix) + "_eff_vs_energy.png";
  c->SaveAs(png.c_str());

  // Write ROOT file
  std::string rootf = std::string(outprefix) + "_eff_vs_energy.root";
  TFile fout(rootf.c_str(),"RECREATE");
  g->Write("eff_vs_energy");
  // Store table
  TTree t("EffTable","Efficiency per energy point");
  double te_E, te_eff;
  long long te_tot, te_det;
  std::string te_file;
  std::string* pFile = &te_file;
  t.Branch("E_MeV",&te_E);
  t.Branch("eff_percent",&te_eff);
  t.Branch("totalEvents",&te_tot);
  t.Branch("eventsAllDet",&te_det);
  t.Branch("filename",&pFile);
  for (auto& p : points) {
    te_E   = p.E_MeV;
    te_eff = p.eff_percent;
    te_tot = p.totalEvents;
    te_det = p.eventsAllDet;
    te_file= p.filename;
    t.Fill();
  }
  t.Write();
  fout.Close();

  // Text summary
  std::string txt = std::string(outprefix) + "_eff_vs_energy.txt";
  std::ofstream ofs(txt);
  ofs << "# Energy(MeV)  Efficiency(%)  TotalEvents  EventsAllDet  File\n";
  for (auto& p : points) {
    ofs << p.E_MeV << "  "
        << p.eff_percent << "  "
        << p.totalEvents << "  "
        << p.eventsAllDet << "  "
        << p.filename << "\n";
  }
  ofs.close();

  std::cout << "DONE: " << png << "\n"
            << "      " << rootf << "\n"
            << "      " << txt << std::endl;
}