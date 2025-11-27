// MakeResponseForUnfolding.C
//
// Construit un TH2D de réponse (Emeas vs Etrue) à partir du TTree "resp"
// et le sauve dans un fichier ROOT "prêt pour unfolding".
//
// Hypothèses sur le TTree "resp" :
//   branches : Etrue_keV, Emeas_keV, parisIndex (optionnel)
//
// Usage typique :
//   .L MakeResponseForUnfolding.C+
//   MakeResponseForUnfolding("out_PARIS50_100000.root",
//                            "Response_PARIS50.root",
//                            400, 0, 15000,  // nbinsTrue, trueMin, trueMax
//                            400, 0, 15000,  // nbinsMeas, measMin, measMax
//                            true, "PARIS50", // resolutionbin, detName
//                            -1);             // parisIndex = -1 => tous

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TVectorD.h>
#include <TMath.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>

struct ResParams {
  double A;
  double power;
};

// Paramètres de résolution PARIS (A, power)
static const std::map<std::string, ResParams> parisRes_MR = { // MR = MakeResponse
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

void MakeResponseForUnfolding(const char* inFile  = "PARIS50/output_PARIS50_12million.root",
                              const char* outFile = "Response_PARIS50.root",
                              int nbinsTrue       = 155,
                              double trueMin      = 0.0,
                              double trueMax      = 15000.0,
                              int nbinsMeas       = 155,
                              double measMin      = 0.0,
                              double measMax      = 15000.0,
                              bool resolutionbin  = true,
                              const char* detName = "PARIS50",
                              int parisIndex      = 0)
{
  // Ouvrir fichier d'entrée
  TFile* fIn = TFile::Open(inFile, "READ");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "[ERROR] Cannot open input file: " << inFile << std::endl;
    return;
  }

  TTree* t = dynamic_cast<TTree*>(fIn->Get("resp"));
  if (!t) {
    std::cerr << "[ERROR] TTree 'resp' not found in " << inFile << std::endl;
    fIn->Close();
    return;
  }

  std::cout << "[INFO] Opened " << inFile
            << ", resp entries = " << t->GetEntries() << std::endl;

  // --- Construire le binning (uniforme ou résolution) ---
  std::vector<double> trueEdges;
  std::vector<double> measEdges;

  if (resolutionbin) {
    std::string key = detName ? std::string(detName) : std::string();
    auto it = parisRes_MR.find(key);
    if (it == parisRes_MR.end()) {
      std::cerr << "[WARN] resolutionbin=true but detName='" << key
                << "' not found in map. Falling back to uniform binning."
                << std::endl;
      resolutionbin = false;
    } else {
      double resA = it->second.A;
      double respower = it->second.power;
      std::cout << "[INFO] Using resolution binning for " << key
                << " : resA=" << resA << " respower=" << respower << std::endl;

      if (!(resA != 0.0 && resA < 100.0 && respower != 0.0 && respower < 1.0)) {
        std::cerr << "[WARN] resA/respower out of range, fallback to uniform."
                  << std::endl;
        resolutionbin = false;
      } else {
        auto buildEdges = [&](int nbins, double Emin) {
          std::vector<double> edges(nbins+1);
          edges[0] = Emin;
          edges[1] = Emin + 11.0;  // bin initial
          for (int i = 2; i < nbins+1; ++i) {
            double prev = edges[i-1];
            double step = resA * TMath::Power(prev, respower) * prev;
            edges[i] = prev + step;
          }
          return edges;
        };

        trueEdges = buildEdges(nbinsTrue, trueMin);
        measEdges = buildEdges(nbinsMeas, measMin);

        std::cout << "[INFO] True edges: from " << trueEdges.front()
                  << " to " << trueEdges.back() << " keV" << std::endl;
        std::cout << "[INFO] Meas edges: from " << measEdges.front()
                  << " to " << measEdges.back() << " keV" << std::endl;
      }
    }
  }

  // Construire le TH2D réponse (X=Emeas, Y=Etrue)
  TH2D* hResp = nullptr;
  if (resolutionbin && !trueEdges.empty() && !measEdges.empty()) {
    hResp = new TH2D("hResp", "Response;E_{meas} [keV];E_{true} [keV]",
                     nbinsMeas, measEdges.data(),
                     nbinsTrue, trueEdges.data());
  } else {
    hResp = new TH2D("hResp", "Response;E_{meas} [keV];E_{true} [keV]",
                     nbinsMeas, measMin, measMax,
                     nbinsTrue, trueMin, trueMax);
  }
  hResp->SetDirectory(nullptr);

  // Branches
  double Etrue=0.0, Emeas=0.0;
  int pIdx = 0;

  t->SetBranchAddress("Etrue_keV", &Etrue);
  t->SetBranchAddress("Emeas_keV", &Emeas);
  // parisIndex est optionnel : on teste la présence
  bool hasParisIndex = (t->GetBranch("parisIndex") != nullptr);
  if (hasParisIndex) {
    t->SetBranchAddress("parisIndex", &pIdx);
    std::cout << "[INFO] Branch 'parisIndex' found, will apply filter if parisIndex>=0"
              << std::endl;
  } else {
    std::cout << "[INFO] Branch 'parisIndex' NOT found, ignoring filter." << std::endl;
  }

  // Remplissage
  Long64_t n = t->GetEntries();
  for (Long64_t i=0; i<n; ++i) {
    t->GetEntry(i);
    if (hasParisIndex && parisIndex >= 0 && pIdx != parisIndex) continue;
    // Here we fill X=Emeas, Y=Etrue
    hResp->Fill(Emeas, Etrue);
  }

  std::cout << "[INFO] Filled hResp, total integral = " << hResp->Integral() << std::endl;

  // Sauvegarde dans outFile
  TFile fOut(outFile, "RECREATE");
  hResp->Write("hResp");

  if (resolutionbin && !trueEdges.empty() && !measEdges.empty()) {
    TVectorD vTrue(trueEdges.size());
    TVectorD vMeas(measEdges.size());
    for (size_t i=0;i<trueEdges.size();++i) vTrue[i]=trueEdges[i];
    for (size_t i=0;i<measEdges.size();++i) vMeas[i]=measEdges[i];
    vTrue.Write("TrueBinEdges");
    vMeas.Write("MeasBinEdges");
  }

  fOut.Close();
  fIn->Close();

  std::cout << "[INFO] Response histogram saved to " << outFile << std::endl;
}