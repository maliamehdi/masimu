// ExtractBinCenters.C
// ------------------------------------------------------------
// Extrait les centres de bins d'un histogramme TH1 et écrit une macro Geant4:
//   /control/alias Elist { c1 c2 c3 ... }
//
// - Si Emax_keV < 0 : garde TOUS les bin centers
// - Si Emax_keV >= 0 : garde uniquement bin center < Emax_keV
// Imprime le nombre total de bins et le nombre conservé.
//
// Usage :
//   root -l -q 'ExtractBinCenters.C("input.root","resbinspectrumPARIS235",-1,"")'
//   root -l -q 'ExtractBinCenters.C("input.root","resbinspectrumPARIS235",15000.0,"energies_list_PARIS235.mac")'

#include <TFile.h>
#include <TH1.h>
#include <TAxis.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>

static std::string GuessParisTagFromName(const std::string& s)
{
  // cherche "PARIS" puis lit les chiffres qui suivent
  const std::string key = "PARIS";
  size_t p = s.find(key);
  if (p == std::string::npos) return "UNKNOWN";
  p += key.size();
  std::string digits;
  while (p < s.size() && std::isdigit((unsigned char)s[p])) {
    digits.push_back(s[p]);
    ++p;
  }
  if (digits.empty()) return "UNKNOWN";
  return "PARIS" + digits;
}

void ExtractBinCenters(const char* inFile     = "137Cs_0305_CheckEnergyspectra_all.root",
                       const char* histName   = "resbinspectrumPARIS90",
                       double      Emax_keV   = 15000.0,   // <0 => pas de coupe énergie
                       const char* outMacFile = "")        // "" => auto energies_list_PARIS###.mac
{
  // ---------------- Open input ----------------
  TFile* fIn = TFile::Open(inFile, "READ");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "[ERROR] Cannot open input file: " << inFile << "\n";
    return;
  }

  TH1* h = dynamic_cast<TH1*>(fIn->Get(histName));
  if (!h) {
    std::cerr << "[ERROR] Histogram '" << histName << "' not found in " << inFile << "\n";
    fIn->Close();
    return;
  }

  const int nbinsTot = h->GetNbinsX();
  const TAxis* ax = h->GetXaxis();

  // ---------------- Select bins ----------------
  std::vector<double> centers;
  centers.reserve(nbinsTot);

  for (int i = 1; i <= nbinsTot; ++i) {
    double center = ax->GetBinCenter(i);
    if (Emax_keV < 0.0 || center < Emax_keV) {
      centers.push_back(center);
    }
  }

  const int nbinsKept = (int)centers.size();

  // ---------------- Print summary ----------------
  std::cout << "=====================================\n";
  std::cout << "[INFO] Histogram        : " << histName << "\n";
  std::cout << "[INFO] Total nbins      : " << nbinsTot << "\n";
  if (Emax_keV < 0.0)
    std::cout << "[INFO] Energy cut       : none (all bin centers kept)\n";
  else
    std::cout << "[INFO] Energy cut       : bin center < " << Emax_keV << " keV\n";
  std::cout << "[INFO] Bins kept        : " << nbinsKept << "\n";
  if (nbinsKept > 0) std::cout << "[INFO] Last kept center : " << centers.back() << " keV\n";
  std::cout << "=====================================\n";

  if (nbinsKept == 0) {
    std::cerr << "[WARN] No bin satisfies the condition. Nothing written.\n";
    fIn->Close();
    return;
  }

  // ---------------- Output file name ----------------
  std::string outName;
  if (outMacFile && outMacFile[0] != '\0') {
    outName = outMacFile;
  } else {
    std::string tag = GuessParisTagFromName(histName ? std::string(histName) : std::string());
    outName = "../simu/build/new_energies_list_" + tag + ".mac";
  }

  // ---------------- Write .mac ----------------
  std::ofstream out(outName.c_str());
  if (!out) {
    std::cerr << "[ERROR] Cannot create output file: " << outName << "\n";
    fIn->Close();
    return;
  }

  out << "/control/alias Elist { ";

  // format proche de ton exemple
  out << std::setprecision(6) << std::noshowpoint;

  for (int i = 0; i < nbinsKept; ++i) {
    // éviter la notation scientifique si possible
    out << std::fixed << std::setprecision(4) << centers[i];
    if (i+1 < nbinsKept) out << " ";
  }

  out << " }\n";
  out.close();
  fIn->Close();

  std::cout << "[INFO] Wrote Geant4 alias list to: " << outName << "\n";
}