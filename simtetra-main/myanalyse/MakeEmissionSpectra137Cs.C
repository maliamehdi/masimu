// MakeEmissionSpectra137Cs.C
// Crée un spectre TH1D pour chaque détecteur PARIS avec le binning "résolution"
// (même règle que PlotResponseMatrix.C), MAIS avec un nombre de bins choisi
// automatiquement pour que le DERNIER BIN CENTER soit STRICTEMENT < Emax_keV (15000 keV par défaut).
// Puis remplit N fois l'énergie Emono_keV (662 keV), et sauve tout dans un fichier ROOT.
//
// Usage :
//   root -l -q 'MakeEmissionSpectra137Cs.C("Emit137Cs_AllPARIS.root", 280146585, 662.0, 15000.0, 0.0)'

#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <TVectorD.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>

struct ResParams { double A; double power; };

static const std::map<std::string, ResParams> parisRes = {
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

// Binning "résolution" :
// edges[0]=Emin, edges[1]=Emin+11, puis next = prev + A * prev^power * prev
// MAIS on s'arrête dès que le centre du bin (prev,next) serait >= Emax_keV.
static std::vector<double> BuildEdgesResolutionUpToLastCenter(double Emin, double Emax_keV,
                                                              double resA, double respower)
{
  std::vector<double> edges;
  edges.reserve(512);

  edges.push_back(Emin);
  edges.push_back(Emin + 11.0);

  // Si on n'a même pas un bin valide avec center<Emax, on renvoie 0 bin => on skip ensuite.
  {
    double prev = edges[0];
    double next = edges[1];
    double center = 0.5 * (prev + next);
    if (center >= Emax_keV) {
      edges.resize(1); // => nbins = 0
      return edges;
    }
  }

  while (true) {
    double prev = edges.back();
    double step = resA * TMath::Power(prev, respower) * prev;
    double next = prev + step;

    double center = 0.5 * (prev + next);
    if (center >= Emax_keV) break; // <-- critère demandé (sur le BIN CENTER)

    edges.push_back(next);

    if (edges.size() > 20000) {
      std::cerr << "[WARN] Too many bins, stopping at " << edges.size() << "\n";
      break;
    }
  }

  return edges; // taille = nbins+1
}

void MakeEmissionSpectra137Cs(const char* outFile   = "Emit137Cs_AllPARIS.root",
                              long long   Nfill     = 280146585LL,
                              double      Emono_keV = 662.0,
                              double      Emax_keV  = 15000.0,
                              double      Emin      = 0.0)
{
  if (Nfill <= 0) { std::cerr << "[ERROR] Nfill must be > 0\n"; return; }
  if (Emax_keV <= Emin) { std::cerr << "[ERROR] Emax_keV must be > Emin\n"; return; }

  TFile* fOut = TFile::Open(outFile, "RECREATE");
  if (!fOut || fOut->IsZombie()) {
    std::cerr << "[ERROR] Cannot create output file: " << outFile << "\n";
    return;
  }

  std::cout << "[INFO] Writing spectra to " << outFile
            << " | Nfill=" << Nfill
            << " | Emono_keV=" << Emono_keV
            << " | condition: last BIN CENTER < " << Emax_keV << " keV\n";

  for (const auto& kv : parisRes) {
    const std::string& det = kv.first;
    const double resA      = kv.second.A;
    const double respower  = kv.second.power;

    std::vector<double> edges = BuildEdgesResolutionUpToLastCenter(Emin, Emax_keV, resA, respower);
    const int nbins = (int)edges.size() - 1;

    if (nbins < 1) {
      std::cerr << "[WARN] " << det << " : nbins<1 (no bin center < Emax), skipping.\n";
      continue;
    }

    std::string hname  = "hMono662_" + det;
    std::string htitle = "Mono 662 keV spectrum (" + det + ");Energy [keV];Counts";

    TH1D* h = new TH1D(hname.c_str(), htitle.c_str(), nbins, edges.data());
    h->SetDirectory(nullptr);

    int b = h->FindBin(Emono_keV);
    if (b < 1 || b > nbins) {
      std::cerr << "[WARN] " << det << " : Emono_keV=" << Emono_keV
                << " outside range [" << edges.front() << ", " << edges.back() << "] keV. Skipping.\n";
      delete h;
      continue;
    }

    h->SetBinContent(b, (double)Nfill);
    h->SetBinError  (b, TMath::Sqrt((double)Nfill));
    h->SetEntries((double)Nfill);   // <-- IMPORTANT: Entries = Nfill (comme si Fill N fois)

    // debug: store edges
    TVectorD vEdges(edges.size());
    for (size_t i=0; i<edges.size(); ++i) vEdges[i] = edges[i];

    // check last bin center
    double lastCenter = 0.5 * (edges[nbins-1] + edges[nbins]);

    fOut->cd();
    h->Write(hname.c_str());
    vEdges.Write(("Edges_" + det).c_str());

    std::cout << "[OK] " << det
              << " : nbins=" << nbins
              << " ; lastEdge=" << edges.back()
              << " ; lastCenter=" << lastCenter << " (<" << Emax_keV << ")"
              << " ; filled bin=" << b
              << " (center=" << h->GetBinCenter(b) << " keV, width=" << h->GetBinWidth(b) << " keV)\n";

    delete h;
  }

  fOut->Close();
  std::cout << "[INFO] Done.\n";
}