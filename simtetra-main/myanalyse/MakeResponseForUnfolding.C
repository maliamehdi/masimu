// MakeResponseForUnfolding.C

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

static const std::map<std::string, ResParams> parisRes_MR = {
  {"PARIS50",  {1.12145,  -0.441244}}, // run du 05/08/2024 paris index 0
  {"PARIS70",  {1.80973,  -0.550685}},// paris index 1
  {"PARIS90",  {1.94868,  -0.564616}},
  {"PARIS110", {2.11922,  -0.582147}},
  {"PARIS130", {0.794233, -0.377311}},
  {"PARIS235", {1.30727,  -0.477402}},
  {"PARIS262", {1.76345,  -0.542769}},
  {"PARIS278", {1.98579,  -0.559095}},
  {"PARIS305", {1.9886,   -0.574021}}
};

// [OK] PARIS110 : nbins=218 ; lastEdge=15004.8 ; lastCenter=14946.1 (<15000) ; b1=49 (E1 center=1180.5 keV) ; b2=53 (E2 center=1346.86 keV) ; Entries=1.45946e+09
// [OK] PARIS130 : nbins=121 ; lastEdge=15124.8 ; lastCenter=14967.8 (<15000) ; b1=43 (E1 center=1191.24 keV) ; b2=45 (E2 center=1324.17 keV) ; Entries=1.45946e+09
// [OK] PARIS235 : nbins=156 ; lastEdge=15046.6 ; lastCenter=14947.7 (<15000) ; b1=44 (E1 center=1149.32 keV) ; b2=47 (E2 center=1308.86 keV) ; Entries=1.45946e+09
// [OK] PARIS262 : nbins=192 ; lastEdge=15019.3 ; lastCenter=14948 (<15000) ; b1=47 (E1 center=1166.89 keV) ; b2=51 (E2 center=1349.62 keV) ; Entries=1.45946e+09
// [OK] PARIS278 : nbins=194 ; lastEdge=15017.5 ; lastCenter=14948.8 (<15000) ; b1=46 (E1 center=1176.24 keV) ; b2=49 (E2 center=1313.01 keV) ; Entries=1.45946e+09
// [OK] PARIS305 : nbins=218 ; lastEdge=15037 ; lastCenter=14977.3 (<15000) ; b1=50 (E1 center=1187.01 keV) ; b2=54 (E2 center=1352.78 keV) ; Entries=1.45946e+09
// [OK] PARIS50 : nbins=138 ; lastEdge=15032.6 ; lastCenter=14912.7 (<15000) ; b1=43 (E1 center=1192.16 keV) ; b2=45 (E2 center=1311.16 keV) ; Entries=1.45946e+09
// [OK] PARIS70 : nbins=199 ; lastEdge=15008 ; lastCenter=14940.2 (<15000) ; b1=48 (E1 center=1174.68 keV) ; b2=52 (E2 center=1352.33 keV) ; Entries=1.45946e+09
// [OK] PARIS90 : nbins=206 ; lastEdge=14967.1 ; lastCenter=14903.3 (<15000) ; b1=48 (E1 center=1167.01 keV) ; b2=52 (E2 center=1339.65 keV) ; Entries=1.45946e+09

void MakeResponseForUnfolding(const char* inFile  = "PARIS235/new_output_PARIS235_10million.root",//"PARIS110/output_PARIS110_10million_complete.root",//PARIS70/output_PARIS70_complete.root",//PARIS50, //PARIS90, //PARIS130, // PARIS262, //PARIS278, //PARIS305
                              const char* outFile = "new_Response_PARIS235.root",
                              int nbinsTrue       = 156,//218,//199,//138,//206, //121, //192, //194, //218
                              double trueMin      = 0.0,
                              double trueMax      = 15000.0,
                              int nbinsMeas       = 156,//207,//199,
                              double measMin      = 0.0,
                              double measMax      = 15000.0,
                              bool resolutionbin  = true,
                              const char* detName = "PARIS235",
                              int parisIndex      = 5,
                              // --- NEW: smoothing option ---
                              bool doSmooth       = false,
                              int  nSmooth        = 1)
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
          edges[1] = Emin + 11.0;
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
    hResp->Fill(Emeas, Etrue);
  }

  std::cout << "[INFO] Filled hResp, total integral = " << hResp->Integral() << std::endl;

  // --- NEW: smoothing (optional) ---
  if (doSmooth) {
    if (nSmooth < 1) nSmooth = 1;
    std::cout << "[INFO] Applying TH2::Smooth(" << nSmooth << ") on response matrix." << std::endl;
    hResp->Smooth(nSmooth);
    std::cout << "[INFO] After Smooth, integral = " << hResp->Integral() << std::endl;
  }

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