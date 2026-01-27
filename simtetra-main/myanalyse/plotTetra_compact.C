// plotTetra_compact.C
// Analyse TETRA (version compacte, sans PARIS) :
// - Multiplicité détectée dans une fenêtre temporelle t_det <= tGate_ns
// - Pas de double comptage : dédoublonnage par (eventID, trackID) dans NeutronPrimaries
// - Deadtime appliqué indépendamment sur chaque ring (sur les temps détectés dans la fenêtre)
//
// Usage ROOT :
//   .L plotTetra_compact.C+
//   plotTetraCompact("output.root");
//   plotTetraCompactFromList("filelist.txt");
//
// Notes importantes :
// - On suppose que "NeutronPrimaries" contient au moins : eventID, trackID, Eemit_MeV, ring, detectTime_ns
// - S'il existe plusieurs entrées NeutronPrimaries pour un même trackID (ex: un neutron qui donne plusieurs enregistrements),
//   on ne compte qu'UNE seule détection par neutron : la détection la plus tardive (plus grand detectTime_ns) avec ring>0.
// - La multiplicité émise est déduite aussi par trackID unique (pas par nombre de lignes).
//
// (c) patched for: éviter le double comptage, garder deadtime par ring.

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TEfficiency.h>
#include <TLegend.h>
#include <TLine.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TString.h>
#include <TLatex.h>

#include <map>
#include <vector>
#include <array>
#include <algorithm>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <memory>

// ---------- helpers ----------
static inline std::string TrimStr(const std::string& s) {
  const char* ws = " \t\r\n";
  size_t b = s.find_first_not_of(ws);
  if (b == std::string::npos) return "";
  size_t e = s.find_last_not_of(ws);
  return s.substr(b, e - b + 1);
}

static std::vector<std::string> ReadFileList(const char* listfile) {
  std::vector<std::string> files;
  std::ifstream fin(listfile);
  if (!fin.is_open()) {
    std::cerr << "Cannot open list file: " << listfile << "\n";
    return files;
  }
  std::string line;
  while (std::getline(fin, line)) {
    std::string s = TrimStr(line);
    if (s.empty()) continue;
    if (s[0] == '#') continue;
    if (s.size() > 1 && s[0] == '/' && s[1] == '/') continue;
    std::istringstream iss(s);
    std::string path;
    if (!(iss >> path)) continue;
    files.push_back(path);
  }
  return files;
}

// ---------- main ----------
void plotTetraCompact(const char* infile = "output_neutron_run01_smeared.root",
                      double tGate_ns = 150000.0,
                      double deadtime_ns = 480.0,
                      double budgetMEV = 20.0)
{
  gStyle->SetOptStat(1110);

  std::unique_ptr<TFile> f(TFile::Open(infile, "READ"));
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open file: " << infile << "\n";
    return;
  }

  auto* tEvents = dynamic_cast<TTree*>(f->Get("Events"));
  auto* tPrim   = dynamic_cast<TTree*>(f->Get("NeutronPrimaries"));

  if (!tEvents) {
    std::cerr << "Missing tree 'Events' in file\n";
    return;
  }
  if (!tPrim) {
    std::cerr << "Missing tree 'NeutronPrimaries' in file (needed for de-doubling + gate)\n";
    return;
  }

  // ----------------- Events branches (pour totaux + contrôle) -----------------
  Int_t    EventID = 0;
  Int_t    NNeutronsEmitted = 0;
  Double_t nDetected = 0.0;     // stored as double in ntuple
  Int_t    NNeutronsEscaped = 0;

  // optionnel : si présent
  Double_t lastNeutronTime_ns = -1.0;
  bool hasLast = (tEvents->GetBranch("lastNeutronTime_ns") != nullptr);

  tEvents->SetBranchAddress("EventID", &EventID);
  tEvents->SetBranchAddress("NNeutronsEmitted", &NNeutronsEmitted);
  tEvents->SetBranchAddress("nDetected", &nDetected);
  tEvents->SetBranchAddress("NNeutronsEscaped", &NNeutronsEscaped);
  if (hasLast) tEvents->SetBranchAddress("lastNeutronTime_ns", &lastNeutronTime_ns);

  const Long64_t totEvents = tEvents->GetEntries();

  // Totaux globaux (Events)
  long long totEmitted_evt  = 0;
  long long totDetected_evt = 0;
  long long totEscaped_evt  = 0;

  // ----------------- Primaries branches -----------------
  Int_t    np_eventID       = 0;
  Int_t    np_trackID       = 0;
  Double_t np_Eemit_MeV     = 0.0;
  Int_t    np_ring          = 0;
  Double_t np_detectTime_ns = -1.0;

  tPrim->SetBranchAddress("eventID",       &np_eventID);
  tPrim->SetBranchAddress("trackID",       &np_trackID);
  tPrim->SetBranchAddress("Eemit_MeV",     &np_Eemit_MeV);
  tPrim->SetBranchAddress("ring",          &np_ring);
  tPrim->SetBranchAddress("detectTime_ns", &np_detectTime_ns);

  // ----------------- Structures de dédoublonnage -----------------
  struct NeutronRec {
    bool   hasEmit = false;     // a été vu au moins une fois (ligne primaire)
    double Eemit   = 0.0;

    bool   hasDet  = false;     // au moins une détection (ring>0, time>=0)
    int    ringDet = 0;         // ring de la détection retenue
    double tDet    = -1.0;      // temps de la détection retenue (min)
  };

  // event -> track -> neutron record
  std::map<int, std::map<int, NeutronRec>> recByEvt;

  // accum max energy + sum energy per event (dédoublés par track)
  double maxEemit = 0.0;
  std::map<int,double> sumEperEvt;

  // Pass unique : remplir recByEvt
  const Long64_t nPrim = tPrim->GetEntries();
  for (Long64_t i=0; i<nPrim; ++i) {
    tPrim->GetEntry(i);

    auto &byTrack = recByEvt[np_eventID];
    auto &rec = byTrack[np_trackID];

    // emitted info (une fois)
    if (!rec.hasEmit) {
      rec.hasEmit = true;
      rec.Eemit   = np_Eemit_MeV;
      if (rec.Eemit > maxEemit) maxEemit = rec.Eemit;
      sumEperEvt[np_eventID] += rec.Eemit;
    }

    // detected info: garder uniquement la détection la plus tardive pour ce track
    if (np_ring > 0 && np_ring <= 4 && np_detectTime_ns >= 0.0) {
      if (!rec.hasDet || np_detectTime_ns > rec.tDet) {
        rec.hasDet  = true;
        rec.ringDet = np_ring;
        rec.tDet    = np_detectTime_ns;
      }
    }
  }

  if (maxEemit <= 0) maxEemit = 20.0;
  const int eBinsEmit = std::clamp((int)std::ceil(maxEemit*10.0), 100, 2000);

  // ----------------- Histos -----------------
  // Ranges multiplicité (on ne connaît pas a priori : on prend le max émis/détecté dédoublé)
  int maxMultEm = 1;
  int maxMultDetGate = 1;
  int maxSumRingsGate = 1;

  // Précalc maxima sur les maps (dédoublées)
  for (auto &kvE : recByEvt) {
    auto &m = kvE.second;

    int nEm = (int)m.size(); // track uniques = émis
    maxMultEm = std::max(maxMultEm, nEm);

    std::array<int,4> ringCounts = {0,0,0,0};
    int nDetGate = 0;
    for (auto &kvT : m) {
      const NeutronRec &r = kvT.second;
      if (r.hasDet && r.tDet >= 0.0 && r.tDet <= tGate_ns) {
        ++nDetGate;
        if (r.ringDet>=1 && r.ringDet<=4) ringCounts[r.ringDet-1]++;
      }
    }
    maxMultDetGate   = std::max(maxMultDetGate, nDetGate);
    maxSumRingsGate  = std::max(maxSumRingsGate, ringCounts[0]+ringCounts[1]+ringCounts[2]+ringCounts[3]);
  }

  // Find max detected multiplicity from Events tree (fast scan) to match binning
  int maxMultEventsTree = 1;
  for (Long64_t i=0; i<totEvents; ++i) {
    tEvents->GetEntry(i);
    int nd = (int)std::llround(nDetected);
    if (nd > maxMultEventsTree) maxMultEventsTree = nd;
  }
  int maxMultBins = std::max(maxMultDetGate, maxMultEventsTree);

  auto hMultEmitted = new TH1I("hMultEmitted",
                               "Emitted neutron multiplicity per event (dedup by trackID);N emitted;Events",
                               maxMultEm+1, -0.5, maxMultEm+0.5);

  auto hMultGate = new TH1I("hMultGate",
                            "Measured neutron multiplicity per event in gate (dedup by trackID);N detected (t_{det}<=gate);Events",
                            maxMultBins+1, -0.5, maxMultBins+0.5);

  // Direct multiplicity as recorded in Events ntuple (fast, no dedup)
  auto hMultEvents = new TH1I("hMultEvents",
                              "Detected multiplicity from Events tree;N detected;Events",
                              maxMultBins+1, -0.5, maxMultBins+0.5);

  auto hR1 = new TH1I("hR1", "Ring 1 multiplicity in gate (dedup);HitsRing1;Events",
                      maxMultBins+1, -0.5, maxMultBins+0.5);
  auto hR2 = new TH1I("hR2", "Ring 2 multiplicity in gate (dedup);HitsRing2;Events",
                      maxMultBins+1, -0.5, maxMultBins+0.5);
  auto hR3 = new TH1I("hR3", "Ring 3 multiplicity in gate (dedup);HitsRing3;Events",
                      maxMultBins+1, -0.5, maxMultBins+0.5);
  auto hR4 = new TH1I("hR4", "Ring 4 multiplicity in gate (dedup);HitsRing4;Events",
                      maxMultBins+1, -0.5, maxMultBins+0.5);

  auto hSumRings = new TH1I("hSumRings",
                            "Sum of hits over 4 rings in gate (dedup);Sum rings (t_{det}<=gate);Events",
                            maxSumRingsGate+1, -0.5, maxSumRingsGate+0.5);

  auto hSumRingsDead = new TH1I("hSumRings_deadtime",
                                Form("Sum of hits over 4 rings with deadtime %.0f ns per ring (t_{det}<=gate);Deadtime-corrected sum;Events", deadtime_ns),
                                maxSumRingsGate+1, -0.5, maxSumRingsGate+0.5);

  auto hEemitAll = new TH1D("hEemitAll",
                            "Emitted neutron energy (dedup by trackID);E_{emit} [MeV];Counts",
                            eBinsEmit, 0, maxEemit);

  auto hEemitDet = new TH1D("hEemitDet",
                            "Detected neutrons energy (dedup, ring>0);E_{emit} [MeV];Counts",
                            eBinsEmit, 0, maxEemit);

  auto hEffVsE = new TEfficiency("hEffVsE",
                                 "Detection efficiency vs emitted energy (dedup);E_{emit} [MeV];Efficiency",
                                 eBinsEmit, 0, maxEemit);

  auto hEffVsE_gateAllDetInGate = new TEfficiency(
      "hEffVsE_gateAllDetInGate",
      "Efficiency vs E_{emit} for events where all detections are in gate (dedup);E_{emit} [MeV];Efficiency",
      eBinsEmit, 0, maxEemit);

  // Sum emitted energy per event
  double maxSum = 0.0;
  for (auto &kv : sumEperEvt) if (kv.second > maxSum) maxSum = kv.second;
  double sumMaxAxis = std::max(maxSum*1.1, budgetMEV*1.5);
  auto hSumEemit = new TH1D("hSumEemit",
                            "Sum emitted energy per event (dedup);Sum E [MeV];Events",
                            120, 0, sumMaxAxis);

  // ----------------- Remplissage des histos via recByEvt -----------------
  long long nEventsAllDetBeforeGate = 0;
  long long nEventsAnyBeforeGate    = 0;

  long long overBudgetEvents = 0;

  // Pour moyennes
  double sumMultGate = 0.0;

  for (auto &kvE : recByEvt) {
    auto &m = kvE.second;

    // emitted multiplicity
    hMultEmitted->Fill((int)m.size());

    // per event energy sum
    double sumE = 0.0;
    for (auto &kvT : m) sumE += kvT.second.Eemit;
    hSumEemit->Fill(sumE);
    if (sumE > budgetMEV) ++overBudgetEvents;

    // Build per-ring times in gate for deadtime
    std::array<std::vector<double>,4> timesRing;
    std::array<int,4> ringsGate = {0,0,0,0};

    bool anyBefore = false;
    bool anyAfter  = false;

    int multGate = 0;

    for (auto &kvT : m) {
      const NeutronRec &r = kvT.second;

      // Energy per neutron
      hEemitAll->Fill(r.Eemit);

      const bool detected = r.hasDet;
      hEffVsE->Fill(detected, r.Eemit);

      if (detected) {
        hEemitDet->Fill(r.Eemit);

        if (r.tDet >= 0.0 && r.tDet <= tGate_ns) anyBefore = true;
        if (r.tDet >  tGate_ns)                 anyAfter  = true;

        if (r.tDet >= 0.0 && r.tDet <= tGate_ns) {
          ++multGate;
          if (r.ringDet>=1 && r.ringDet<=4) {
            ringsGate[r.ringDet-1]++;
            timesRing[r.ringDet-1].push_back(r.tDet);
          }
        }
      }
    }

    if (anyBefore) ++nEventsAnyBeforeGate;
    if (anyBefore && !anyAfter) ++nEventsAllDetBeforeGate;

    // Fill gated multiplicities
    hMultGate->Fill(multGate);
    sumMultGate += multGate;

    hR1->Fill(ringsGate[0]);
    hR2->Fill(ringsGate[1]);
    hR3->Fill(ringsGate[2]);
    hR4->Fill(ringsGate[3]);

    hSumRings->Fill(ringsGate[0]+ringsGate[1]+ringsGate[2]+ringsGate[3]);

    // Deadtime per ring
    int sumDeadEvt = 0;
    for (int r=0; r<4; ++r) {
      auto &v = timesRing[r];
      if (v.empty()) continue;
      std::sort(v.begin(), v.end());
      double lastCounted = -1e30;
      int countRing = 0;
      for (double t : v) {
        if (countRing==0 || (t - lastCounted) >= deadtime_ns) {
          ++countRing;
          lastCounted = t;
        }
      }
      sumDeadEvt += countRing;
    }
    hSumRingsDead->Fill(sumDeadEvt);

    // Gated efficiency (event-level condition) : remplir un point par neutron
    const bool eventAllDetInGate = (anyBefore && !anyAfter);
    for (auto &kvT : m) {
      const NeutronRec &r = kvT.second;
      hEffVsE_gateAllDetInGate->Fill(r.hasDet && eventAllDetInGate, r.Eemit);
    }
  }

  // ----------------- Totaux depuis Events tree (pour contrôle) -----------------
  for (Long64_t i=0; i<totEvents; ++i) {
    tEvents->GetEntry(i);
    totEmitted_evt  += NNeutronsEmitted;
    totDetected_evt += (long long)std::llround(nDetected);
    if (hMultEvents) hMultEvents->Fill((int)std::llround(nDetected));
    totEscaped_evt  += NNeutronsEscaped;
  }

  // ----------------- Résumés / efficacités -----------------
  double effAllDetInGate = (totEvents>0) ? (double)nEventsAllDetBeforeGate/(double)totEvents : 0.0;
  double errAllDetInGate = (totEvents>0) ? std::sqrt(effAllDetInGate*(1.0-effAllDetInGate)/(double)totEvents) : 0.0;

  double effAnyInGate = (totEvents>0) ? (double)nEventsAnyBeforeGate/(double)totEvents : 0.0;
  double errAnyInGate = (totEvents>0) ? std::sqrt(effAnyInGate*(1.0-effAnyInGate)/(double)totEvents) : 0.0;

  double meanMultGate = (totEvents>0) ? (sumMultGate/(double)totEvents) : 0.0;

  double percentEscaped = (totEmitted_evt>0) ? 100.0*(double)totEscaped_evt/(double)totEmitted_evt : 0.0;
  long long diffEmDet_evt = totEmitted_evt - totDetected_evt;

  std::cout << "\n===== plotTetraCompact SUMMARY =====\n";
  std::cout << "File: " << infile << "\n";
  std::cout << "Gate: t_det <= " << tGate_ns << " ns\n";
  std::cout << "Deadtime: " << deadtime_ns << " ns per ring (applied only to detections in gate)\n";
  std::cout << "Events (Entries in Events tree): " << totEvents << "\n";
  std::cout << "Events with ANY detection in gate : " << nEventsAnyBeforeGate
            << " (eff=" << 100.0*effAnyInGate << " +/- " << 100.0*errAnyInGate << " %)\n";
  std::cout << "Events with ALL detections in gate : " << nEventsAllDetBeforeGate
            << " (eff=" << 100.0*effAllDetInGate << " +/- " << 100.0*errAllDetInGate << " %)\n";
  std::cout << "Mean measured multiplicity in gate (dedup): " << meanMultGate << "\n";
  std::cout << "Budget: " << budgetMEV << " MeV, events over budget: " << overBudgetEvents << "\n";
  std::cout << "\n[Control from Events tree]\n";
  std::cout << "Total emitted (Events)  : " << totEmitted_evt << "\n";
  std::cout << "Total detected (Events) : " << totDetected_evt << "\n";
  std::cout << "Total escaped (Events)  : " << totEscaped_evt << " (" << percentEscaped << " %)\n";
  std::cout << "Emitted-Detected (Events): " << diffEmDet_evt << " (compare to escaped)\n";

  // ----------------- Plots -----------------
  TString tin(infile);
  TString dir  = gSystem->DirName(tin);
  TString base = gSystem->BaseName(tin);
  if (base.EndsWith(".root")) base.ReplaceAll(".root", "");
  TString outBase = Form("%s/plot_%s_compact", dir.Data(), base.Data());
  gSystem->Exec(TString::Format("mkdir -p %s", dir.Data()));

  auto c1 = new TCanvas("c1_compact", "Multiplicity (compact)", 1400, 900);
  c1->Divide(3,3);
  c1->cd(1); hMultEmitted->Draw("HIST");
  c1->cd(2);
  // Draw Events-tree multiplicity (raw) and deduped multiplicity (from NeutronPrimaries)
  hMultEvents->SetLineColor(kGray+2);
  hMultEvents->SetFillColor(kGray+2);
  hMultEvents->SetFillStyle(3004);
  hMultEvents->Draw("HIST");
  hMultGate->SetLineColor(kBlue+2);
  hMultGate->SetLineWidth(2);
  hMultGate->Draw("HIST SAME");
  {
    auto legM = new TLegend(0.60,0.65,0.90,0.88);
    legM->AddEntry(hMultEvents, "Events::nDetected (raw)", "f");
    legM->AddEntry(hMultGate, "NeutronPrimaries dedup (in gate)", "l");
    legM->SetBorderSize(0);
    legM->Draw();
  }
  c1->cd(3); hSumRings->Draw("HIST");
  c1->cd(4); hSumRingsDead->Draw("HIST");
  c1->cd(5); hR1->Draw("HIST");
  c1->cd(6); hR2->Draw("HIST");
  c1->cd(7); hR3->Draw("HIST");
  c1->cd(8); hR4->Draw("HIST");
  c1->cd(9);
  {
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.04);
    lat.DrawLatex(0.05,0.85,Form("Gate: t_{det} <= %.0f ns", tGate_ns));
    lat.DrawLatex(0.05,0.78,Form("Deadtime: %.0f ns / ring", deadtime_ns));
    lat.DrawLatex(0.05,0.71,Form("Mean mult (gate): %.4f", meanMultGate));
    lat.DrawLatex(0.05,0.64,Form("Any in gate: %.2f %%", 100.0*effAnyInGate));
    lat.DrawLatex(0.05,0.57,Form("All in gate: %.2f %%", 100.0*effAllDetInGate));
  }
  c1->SaveAs(Form("%s_multiplicities.png", outBase.Data()));

  auto cE = new TCanvas("cE_compact", "Emitted energy + efficiency (compact)", 1200, 500);
  cE->Divide(2,1);

  cE->cd(1);
  hEemitAll->SetLineColor(kBlue+1);
  hEemitDet->SetLineColor(kRed+1);
  hEemitAll->Draw("HIST");
  hEemitDet->Draw("HIST SAME");
  {
    auto leg = new TLegend(0.58,0.75,0.88,0.88);
    leg->AddEntry(hEemitAll, "All emitted (dedup)", "l");
    leg->AddEntry(hEemitDet, "Detected (dedup)", "l");
    leg->Draw();
  }

  cE->cd(2);
  cE->cd(2)->SetGrid();
  hEffVsE->SetMarkerStyle(20);
  hEffVsE->Draw("AP");
  hEffVsE_gateAllDetInGate->SetMarkerStyle(24);
  hEffVsE_gateAllDetInGate->Draw("P SAME");
  {
    auto leg = new TLegend(0.48,0.75,0.88,0.88);
    leg->AddEntry(hEffVsE, "Eff (detected/emit)", "p");
    leg->AddEntry(hEffVsE_gateAllDetInGate, "Eff (event all det in gate)", "p");
    leg->Draw();
  }

  cE->SaveAs(Form("%s_energy_eff.png", outBase.Data()));

  auto cS = new TCanvas("cS_compact", "Sum emitted energy per event (compact)", 900, 600);
  hSumEemit->SetLineColor(kGreen+2);
  hSumEemit->Draw("HIST");
  auto lineB = new TLine(budgetMEV, 0, budgetMEV, hSumEemit->GetMaximum()*1.05);
  lineB->SetLineColor(kRed+1); lineB->SetLineStyle(2); lineB->Draw();
  {
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.036);
    lat.DrawLatex(0.60,0.82,Form("Budget = %.2f MeV", budgetMEV));
    lat.DrawLatex(0.60,0.76,Form("Events over budget: %lld", overBudgetEvents));
  }
  cS->SaveAs(Form("%s_sumEemit.png", outBase.Data()));

  // ----------------- Save ROOT output -----------------
  TString rootOut = Form("%s_hists.root", outBase.Data());
  std::unique_ptr<TFile> fout(TFile::Open(rootOut, "RECREATE"));
  if (!fout || fout->IsZombie()) {
    std::cerr << "[ERROR] Cannot create output ROOT file: " << rootOut << "\n";
    return;
  }

  hMultEmitted->Write();
  hMultEvents->Write();
  hMultGate->Write();
  hR1->Write(); hR2->Write(); hR3->Write(); hR4->Write();
  hSumRings->Write();
  hSumRingsDead->Write();

  hEemitAll->Write();
  hEemitDet->Write();
  hEffVsE->Write();
  hEffVsE_gateAllDetInGate->Write();
  hSumEemit->Write();

  // Summary tree
  auto* tSummary = new TTree("Summary", "Compact summary (dedup by trackID)");
  Long64_t out_totEvents = totEvents;
  Long64_t out_totEmitted_evt = totEmitted_evt;
  Long64_t out_totDetected_evt = totDetected_evt;
  Long64_t out_totEscaped_evt = totEscaped_evt;
  Double_t out_percentEscaped = percentEscaped;

  Double_t out_tGate = tGate_ns;
  Double_t out_deadtime = deadtime_ns;
  Double_t out_meanMultGate = meanMultGate;

  Long64_t out_nAnyGate = nEventsAnyBeforeGate;
  Long64_t out_nAllGate = nEventsAllDetBeforeGate;
  Double_t out_effAny = 100.0*effAnyInGate;
  Double_t out_errAny = 100.0*errAnyInGate;
  Double_t out_effAll = 100.0*effAllDetInGate;
  Double_t out_errAll = 100.0*errAllDetInGate;

  Double_t out_budget = budgetMEV;
  Long64_t out_overBudget = overBudgetEvents;

  tSummary->Branch("totEvents", &out_totEvents, "totEvents/L");
  tSummary->Branch("totEmitted_EventsTree", &out_totEmitted_evt, "totEmitted_EventsTree/L");
  tSummary->Branch("totDetected_EventsTree", &out_totDetected_evt, "totDetected_EventsTree/L");
  tSummary->Branch("totEscaped_EventsTree", &out_totEscaped_evt, "totEscaped_EventsTree/L");
  tSummary->Branch("percentEscaped_EventsTree", &out_percentEscaped, "percentEscaped_EventsTree/D");

  tSummary->Branch("tGate_ns", &out_tGate, "tGate_ns/D");
  tSummary->Branch("deadtime_ns", &out_deadtime, "deadtime_ns/D");
  tSummary->Branch("meanMeasuredMultiplicityInGate", &out_meanMultGate, "meanMeasuredMultiplicityInGate/D");

  tSummary->Branch("eventsAnyDetectionBeforeGate", &out_nAnyGate, "eventsAnyDetectionBeforeGate/L");
  tSummary->Branch("effAnyDetectionBeforeGate_percent", &out_effAny, "effAnyDetectionBeforeGate_percent/D");
  tSummary->Branch("errEffAnyDetectionBeforeGate_percent", &out_errAny, "errEffAnyDetectionBeforeGate_percent/D");

  tSummary->Branch("eventsAllDetectionsBeforeGate", &out_nAllGate, "eventsAllDetectionsBeforeGate/L");
  tSummary->Branch("effAllDetectionsBeforeGate_percent", &out_effAll, "effAllDetectionsBeforeGate_percent/D");
  tSummary->Branch("errEffAllDetectionsBeforeGate_percent", &out_errAll, "errEffAllDetectionsBeforeGate_percent/D");

  tSummary->Branch("budgetMEV", &out_budget, "budgetMEV/D");
  tSummary->Branch("eventsOverBudget", &out_overBudget, "eventsOverBudget/L");

  tSummary->Fill();
  tSummary->Write();

  fout->Close();

  std::cout << "\nSaved:\n";
  std::cout << "  PNG : " << outBase << "_*.png\n";
  std::cout << "  ROOT: " << rootOut << "\n";
}

// ----------------------------------------------------------------------
// Wrapper : parse a file list and call plotTetraCompact for each file
// Format (txt):
//   path/to/file.root   [tGate_ns]  [deadtime_ns]  [budgetMEV]
// Ici, pour rester compact, on utilise les valeurs par défaut fournies.
// ----------------------------------------------------------------------
void plotTetraCompactFromList(const char* listfile,
                              double default_tGate_ns = 150000.0,
                              double default_deadtime_ns = 720.0,
                              double default_budgetMEV = 20.0)
{
  std::vector<std::string> files = ReadFileList(listfile);
  if (files.empty()) return;

  for (size_t i=0; i<files.size(); ++i) {
    const std::string& path = files[i];
    std::cout << "--------------------------------------------------\n";
    std::cout << "File " << (i+1) << "/" << files.size() << " : " << path << "\n";
    plotTetraCompact(path.c_str(), default_tGate_ns, default_deadtime_ns, default_budgetMEV);
  }
}
