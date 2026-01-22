// plotTetra.C
// Analyse TETRA : multiplicité détectée dans une fenêtre temporelle,
// compatibilité avec NeutronPrimaries (detectTime_ns) et deadtime par ring.

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TEfficiency.h>
#include <TLegend.h>
#include <TLine.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TString.h>
#include <TROOT.h>
#include <TLatex.h>

#include <map>
#include <vector>
#include <array>
#include <algorithm>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>    // sqrt
#include <memory>   // std::unique_ptr

// Usage ROOT :
// .L plotTetra.C+
// plotTetra("output_neutron_run01_smeared.root")
// plotTetraFromList("filelist.txt")

void plotTetra(const char* infile = "output_neutron_run01_smeared.root",
               const char* /*outprefix*/ = "../myanalyse/plots/plotTetra_analysis",
               double tGate_ns = 150000.0,
               double budgetMEV = 14.0)
{
  gStyle->SetOptStat(1110);

  // ----------------- Open file and get trees -----------------
  std::unique_ptr<TFile> f(TFile::Open(infile, "READ"));
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open file: " << infile << std::endl;
    return;
  }

  auto* tEvents = dynamic_cast<TTree*>(f->Get("Events"));
  auto* tHits   = dynamic_cast<TTree*>(f->Get("TritonHits"));
  auto* tParis  = dynamic_cast<TTree*>(f->Get("ParisEdep"));
  auto* tPrim   = dynamic_cast<TTree*>(f->Get("NeutronPrimaries"));

  if (!tEvents) {
    std::cerr << "Missing tree 'Events' in file" << std::endl;
    return;
  }
  if (!tHits) {
    std::cerr << "Missing tree 'TritonHits' in file (needed for avg detection time / neutron times)" << std::endl;
  }

  // ----------------- Branch setup for Events -----------------
  Int_t    EventID = 0;
  Int_t    NNeutronsEmitted = 0;
  Double_t nDetected = 0.0;  // stored as double in ntuple
  Double_t EdepCe_keV = 0.0;
  Double_t EdepNaI_keV = 0.0;
  Int_t    HitsRing1 = 0, HitsRing2 = 0, HitsRing3 = 0, HitsRing4 = 0;
  Double_t MeanThermTime_ns = -1.0;
  Int_t    NNeutronsEscaped = 0;
  Double_t lastNeutronTime_ns = -1.0;

  tEvents->SetBranchAddress("EventID", &EventID);
  tEvents->SetBranchAddress("NNeutronsEmitted", &NNeutronsEmitted);
  tEvents->SetBranchAddress("nDetected", &nDetected);
  tEvents->SetBranchAddress("EdepCe_keV", &EdepCe_keV);
  tEvents->SetBranchAddress("EdepNaI_keV", &EdepNaI_keV);
  tEvents->SetBranchAddress("HitsRing1", &HitsRing1);
  tEvents->SetBranchAddress("HitsRing2", &HitsRing2);
  tEvents->SetBranchAddress("HitsRing3", &HitsRing3);
  tEvents->SetBranchAddress("HitsRing4", &HitsRing4);
  tEvents->SetBranchAddress("MeanThermTime_ns", &MeanThermTime_ns);
  tEvents->SetBranchAddress("NNeutronsEscaped", &NNeutronsEscaped);
  tEvents->SetBranchAddress("lastNeutronTime_ns", &lastNeutronTime_ns);

  // ----------------- Determine histogram ranges -----------------
  auto getMaxI = [&](const char* br) { return (int)std::max(0.0, tEvents->GetMaximum(br)); };
  auto getMaxD = [&](const char* br) { return std::max(0.0, tEvents->GetMaximum(br)); };

  int maxDet = std::max(10, getMaxI("nDetected"));
  int maxR1  = std::max(10, getMaxI("HitsRing1"));
  int maxR2  = std::max(10, getMaxI("HitsRing2"));
  int maxR3  = std::max(10, getMaxI("HitsRing3"));
  int maxR4  = std::max(10, getMaxI("HitsRing4"));

  double maxEce  = getMaxD("EdepCe_keV");          if (maxEce  <= 0) maxEce  = 5000.0; // keV
  double maxLast = getMaxD("lastNeutronTime_ns");  if (maxLast <= 0) maxLast = 1.0;    // ns

  int maxSumRings = maxR1 + maxR2 + maxR3 + maxR4;

  // ----------------- Histograms -----------------
  auto hMult = new TH1I(
      "hMult",
      "Measured neutron multiplicity per event in gate;N detected (t_{det} <= gate);Events",
      maxDet+1, -0.5, maxDet+0.5);

  auto hR1 = new TH1I("hR1", "Ring 1 multiplicity in gate;HitsRing1 (t_{det} <= gate);Events",
                      maxR1+1, -0.5, maxR1+0.5);
  auto hR2 = new TH1I("hR2", "Ring 2 multiplicity in gate;HitsRing2 (t_{det} <= gate);Events",
                      maxR2+1, -0.5, maxR2+0.5);
  auto hR3 = new TH1I("hR3", "Ring 3 multiplicity in gate;HitsRing3 (t_{det} <= gate);Events",
                      maxR3+1, -0.5, maxR3+0.5);
  auto hR4 = new TH1I("hR4", "Ring 4 multiplicity in gate;HitsRing4 (t_{det} <= gate);Events",
                      maxR4+1, -0.5, maxR4+0.5);

  auto hSumRings = new TH1I(
      "hSumRings",
      "Sum of hits over 4 rings in gate;HitsRing1+2+3+4 (t_{det} <= gate);Events",
      maxSumRings+1, -0.5, maxSumRings+0.5);

  auto hSumRingsDead = new TH1I(
      "hSumRings_deadtime",
      "Sum of hits over 4 rings with 400 ns deadtime per ring (t_{det} <= gate);Deadtime-corrected sum;Events",
      maxSumRings+1, -0.5, maxSumRings+0.5);

  int eBins = std::clamp((int)std::ceil(maxEce/2.0), 100, 5000);
  auto hEce = new TH1D("hEce", "Energy in Ce;EdepCe [keV];Events",
                       eBins, 0, maxEce);

  int tBins = 200;
  auto hLast = new TH1D("hLast",
                        "Last neutron detection time per event;time [ns];Events",
                        tBins, 0, maxLast);

  // Totals
  long long totEvents  = tEvents->GetEntries();
  long long totEmitted = 0;
  long long totDetected = 0;
  long long totEscaped  = 0;

  // Event-level last detection time (from Events tree)
  std::map<int,double> lastByEvt;

  // Mean emitted multiplicity (filled if NeutronPrimaries is present)
  double meanEmittedMultiplicity = -1.0;

  // For means
  double sumDetectedD  = 0.0;  // sum of multiplicity (gated, sera rempli via NeutronPrimaries)
  double sumLastDetect = 0.0;
  int    cntLastDetect = 0;

  // Time-gated efficiency counters (calculés via NeutronPrimaries)
  long long nEventsAllDetBeforeGate = 0; // events where all detections are in gate
  long long nEventsAnyBeforeGate    = -1; // events with at least one detection in gate

  // ----------------- Loop over Events tree (totaux + hLast + ECe) -----------------
  for (Long64_t i = 0; i < totEvents; ++i) {
    tEvents->GetEntry(i);

    // Ce spectrum (ignore zeros)
    if (EdepCe_keV > 0) hEce->Fill(EdepCe_keV);

    // last detection time (ignore -1/no-detection)
    if (lastNeutronTime_ns >= 0) {
      hLast->Fill(lastNeutronTime_ns);
      sumLastDetect += lastNeutronTime_ns;
      ++cntLastDetect;
    }

    // Store last time per event (encore utilisé pour certaines choses si besoin)
    lastByEvt[EventID] = lastNeutronTime_ns;

    // Totaux "physiques" sur toute la fenêtre
    totEmitted  += NNeutronsEmitted;
    totDetected += (long long)std::llround(nDetected);
    totEscaped  += NNeutronsEscaped;
  }

  double percentEscaped = (totEmitted > 0)
                          ? 100.0 * (double)totEscaped / (double)totEmitted
                          : 0.0;
  long long diffEmDet = totEmitted - totDetected;

  // ----------------- TritonHits: moyennes de temps et hDetTime -----------------
  TH1D* hAvgDetTime = nullptr;
  TH1D* hDetTime    = nullptr;
  double meanAvgDetTime = -1.0;

  if (tHits) {
    Int_t    eID  = 0;
    Double_t t_ns = 0.0;

    tHits->SetBranchAddress("EventID", &eID);
    tHits->SetBranchAddress("time_ns", &t_ns);

    std::map<int, std::vector<double>> timesByEvtAll;

    double maxAvg  = 0.0;
    double maxTime = 0.0;

    const Long64_t nHits = tHits->GetEntries();
    for (Long64_t i = 0; i < nHits; ++i) {
      tHits->GetEntry(i);
      if (t_ns < 0) continue;

      timesByEvtAll[eID].push_back(t_ns);
      if (t_ns > maxTime) maxTime = t_ns;
    }

    // Moyenne par event
    std::vector<double> avgs;
    avgs.reserve(timesByEvtAll.size());

    for (auto& kv : timesByEvtAll) {
      const auto& v = kv.second;
      if (v.empty()) continue;

      double s = 0.0;
      for (double x : v) s += x;
      double m = s / (double)v.size();
      avgs.push_back(m);
      if (m > maxAvg) maxAvg = m;
    }

    if (maxAvg <= 0) maxAvg = 1.0;
    hAvgDetTime = new TH1D(
        "hAvgDetTime",
        "Average neutron detection time per event;time [ns];Events",
        200, 0, maxAvg);

    double sumAvg = 0.0;
    int    cntAvg = 0;
    for (double m : avgs) {
      hAvgDetTime->Fill(m);
      sumAvg += m;
      ++cntAvg;
    }
    if (cntAvg > 0) meanAvgDetTime = sumAvg / (double)cntAvg;

    // Histogramme des temps individuels (bin ~1 ns)
    if (maxTime <= 0) maxTime = 1.0;
    int nBinsTime = (int)std::ceil(maxTime) + 1;

    hDetTime = new TH1D(
        "hDetTime",
        "Neutron detection time (all TritonHits);time [ns];Counts",
        nBinsTime, 0, (double)nBinsTime);

    // 2e passe pour remplir hDetTime
    for (Long64_t i = 0; i < nHits; ++i) {
      tHits->GetEntry(i);
      if (t_ns < 0) continue;
      hDetTime->Fill(t_ns);
    }
  }

  // ----------------- Efficiences (bins) : seront remplies via NeutronPrimaries -----------------
  double effAllDetBeforeGate_frac = 0.0;
  double effAllDetBeforeGate_percent = 0.0;
  double errEffAllDetBeforeGate_percent = 0.0;

  double effAnyBeforeGate_frac   = -1.0;
  double effAnyBeforeGate_percent = -1.0;
  double errEffAnyBeforeGate_percent = -1.0;

  // ----------------- Emitted energy and multiplicity (NeutronPrimaries) -----------------
  TH1D* hEemitAll = nullptr;
  TH1D* hEemitDet = nullptr;
  TH1D* hEemitRing[4] = {nullptr, nullptr, nullptr, nullptr};
  TEfficiency* hEffVsE      = nullptr;
  TEfficiency* hEffVsE_gate = nullptr;
  TH1I* hMultEmitted = nullptr;
  TH1D* hSumEemit    = nullptr;
  int   overBudgetEvents = -1;

  if (tPrim) {
    // Branches
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

    // --- structures de travail ---
    // multiplicité émise
    std::map<int,int>    multByEvent;
    // somme d'énergie émise par évènement
    std::map<int,double> sumEperEvt;
    double maxEemit = 0.0;

    // multiplicité détectée dans la fenêtre t_detect <= tGate_ns
    std::map<int,int> multGateByEvt;

    // HitsRing dans la fenêtre, par évènement
    std::map<int, std::array<int,4>> hitsRingGateByEvt; // [0..3] -> rings 1..4

    // Temps de détection (dans la fenêtre) pour le deadtime par ring
    std::map<int, std::array<std::vector<double>,4>> timesByEvtPerRing;

    // Pour l'efficacité par évènement (all/any in gate)
    struct DetGateInfo {
      bool hasDet = false;
      bool anyBeforeGate = false;
      bool anyAfterGate  = false;
    };
    std::map<int,DetGateInfo> detGateByEvt;

    const Long64_t nPrim = tPrim->GetEntries();
    for (Long64_t i = 0; i < nPrim; ++i) {
      tPrim->GetEntry(i);

      // Multiplicité émise
      multByEvent[np_eventID]++;
      if (np_Eemit_MeV > maxEemit) maxEemit = np_Eemit_MeV;
      sumEperEvt[np_eventID] += np_Eemit_MeV;

      // Neutrons détectés (ring > 0)
      if (np_ring > 0 && np_ring <= 4 && np_detectTime_ns >= 0.0) {
        // Info pour les évènements (efficacités)
        auto &info = detGateByEvt[np_eventID];
        info.hasDet = true;
        if (np_detectTime_ns <= tGate_ns) {
          info.anyBeforeGate = true;
        } else {
          info.anyAfterGate = true;
        }

        // Gate sur t_detect pour la multiplicité et les rings :
        // on ne compte que les neutrons dont le temps de détection est <= tGate_ns
        if (np_detectTime_ns <= tGate_ns) {
          multGateByEvt[np_eventID]++;

          auto &arr = hitsRingGateByEvt[np_eventID]; // std::array<int,4> init à {0,0,0,0}
          arr[np_ring - 1]++;

          // Pour le deadtime on ne garde aussi que les temps dans la fenêtre
          timesByEvtPerRing[np_eventID][np_ring - 1].push_back(np_detectTime_ns);
        }
      }
    }

    // Echelle d'énergie pour les histogrammes Eemit
    if (maxEemit <= 0) maxEemit = 20.0;
    const int eBinsEmit = std::clamp((int)std::ceil(maxEemit*10.0), 100, 2000);

    // Create histograms en énergie
    hEemitAll = new TH1D("hEemitAll",
                         "Emitted neutron energy;E_{emit} [MeV];Counts",
                         eBinsEmit, 0, maxEemit);
    hEemitDet = new TH1D("hEemitDet",
                         "Detected neutrons (ring>0);E_{emit} [MeV];Counts",
                         eBinsEmit, 0, maxEemit);
    for (int r = 0; r < 4; ++r) {
      hEemitRing[r] = new TH1D(Form("hEemitRing%d", r+1),
                               Form("Detected neutrons in ring %d;E_{emit} [MeV];Counts", r+1),
                               eBinsEmit, 0, maxEemit);
    }

    // Sum emitted energy per event
    double maxSum = 0.0;
    for (auto &kv : sumEperEvt)
      if (kv.second > maxSum) maxSum = kv.second;
    double sumMaxAxis = std::max(maxSum*1.1, budgetMEV*1.5);
    hSumEemit = new TH1D("hSumEemit",
                         "Sum emitted energy per event;Sum E [MeV];Events",
                         120, 0, sumMaxAxis);
    overBudgetEvents = 0;
    for (auto &kv : sumEperEvt) {
      hSumEemit->Fill(kv.second);
      if (kv.second > budgetMEV) ++overBudgetEvents;
    }

    hEffVsE = new TEfficiency("hEffVsE",
                              "Detection efficiency vs emitted energy;E_{emit} [MeV];Efficiency",
                              eBinsEmit, 0, maxEemit);
    hEffVsE_gate = new TEfficiency(
        "hEffVsE_gate",
        Form("Detection efficiency vs E (all t_{det} in gate);E_{emit} [MeV];Efficiency"),
        eBinsEmit, 0, maxEemit);

    // ----------------- Deadtime 400 ns / ring => hSumRings_deadtime -----------------
    std::map<int,int> sumDeadPerEvt;
    const double deadtime_ns = 720.0;

    for (auto &kv : timesByEvtPerRing) {
      int evtID = kv.first;
      auto &perRing = kv.second;

      int sumDeadEvt = 0;

      for (int r = 0; r < 4; ++r) {
        auto &v = perRing[r];
        if (v.empty()) continue;

        std::sort(v.begin(), v.end());
        double lastCounted = -1e30;
        int countRing = 0;

        for (double t : v) {
          if (t < 0) continue;
          if (countRing == 0 || (t - lastCounted) >= deadtime_ns) {
            ++countRing;
            lastCounted = t;
          }
        }
        sumDeadEvt += countRing;
      }
      sumDeadPerEvt[evtID] = sumDeadEvt;
    }

    // ----------------- Remplissage hSumRingsDead, hSumRings, hR1..4 et hMult -----------------
    for (Long64_t iEvt = 0; iEvt < totEvents; ++iEvt) {
      tEvents->GetEntry(iEvt);

      // Deadtime
      int sumDeadEvt = 0;
      auto itD = sumDeadPerEvt.find(EventID);
      if (itD != sumDeadPerEvt.end()) sumDeadEvt = itD->second;
      hSumRingsDead->Fill(sumDeadEvt);

      // Multiplicité dans la fenêtre (t_detect <= tGate_ns)
      int multGateEvt = 0;
      auto itM = multGateByEvt.find(EventID);
      if (itM != multGateByEvt.end()) multGateEvt = itM->second;
      hMult->Fill(multGateEvt);
      sumDetectedD += multGateEvt;

      // HitsRing1..4 dans la fenêtre
      std::array<int,4> ringsGate = {0,0,0,0};
      auto itR = hitsRingGateByEvt.find(EventID);
      if (itR != hitsRingGateByEvt.end()) ringsGate = itR->second;

      hR1->Fill(ringsGate[0]);
      hR2->Fill(ringsGate[1]);
      hR3->Fill(ringsGate[2]);
      hR4->Fill(ringsGate[3]);

      int sumR_gate = ringsGate[0] + ringsGate[1] + ringsGate[2] + ringsGate[3];
      hSumRings->Fill(sumR_gate);
    }

    // ----------------- Efficiences par événement via NeutronPrimaries -----------------
    nEventsAllDetBeforeGate = 0;
    nEventsAnyBeforeGate    = 0;
    for (auto &kv : detGateByEvt) {
      const DetGateInfo &info = kv.second;
      if (!info.hasDet) continue;

      if (info.anyBeforeGate) {
        ++nEventsAnyBeforeGate;  // au moins une détection dans la fenêtre
      }
      // "all detections in gate" => aucune détection après gate, et au moins une avant
      if (!info.anyAfterGate && info.anyBeforeGate) {
        ++nEventsAllDetBeforeGate;
      }
    }

    if (totEvents > 0) {
      effAllDetBeforeGate_frac    = (double)nEventsAllDetBeforeGate / (double)totEvents;
      effAllDetBeforeGate_percent = 100.0 * effAllDetBeforeGate_frac;
      double errFrac = std::sqrt(effAllDetBeforeGate_frac * (1.0 - effAllDetBeforeGate_frac)
                                 / (double)totEvents);
      errEffAllDetBeforeGate_percent = 100.0 * errFrac;

      if (nEventsAnyBeforeGate >= 0) {
        effAnyBeforeGate_frac    = (double)nEventsAnyBeforeGate / (double)totEvents;
        effAnyBeforeGate_percent = 100.0 * effAnyBeforeGate_frac;
        double errFrac2 = std::sqrt(effAnyBeforeGate_frac * (1.0 - effAnyBeforeGate_frac)
                                    / (double)totEvents);
        errEffAnyBeforeGate_percent = 100.0 * errFrac2;
      }
    }

    // ----------------- Second pass: fill energy/efficiency histos -----------------
    for (Long64_t i = 0; i < nPrim; ++i) {
      tPrim->GetEntry(i);
      hEemitAll->Fill(np_Eemit_MeV);

      const bool detected = (np_ring > 0);
      hEffVsE->Fill(detected, np_Eemit_MeV);

      // Pour l'efficacité "gated", on compte seulement les neutrons
      // des évènements dont toutes les détections sont dans la fenêtre.
      bool eventAllDetInGate = false;
      auto itInfo = detGateByEvt.find(np_eventID);
      if (itInfo != detGateByEvt.end()) {
        const DetGateInfo &inf = itInfo->second;
        eventAllDetInGate = (inf.hasDet && inf.anyBeforeGate && !inf.anyAfterGate);
      }
      hEffVsE_gate->Fill(detected && eventAllDetInGate, np_Eemit_MeV);

      if (detected) {
        hEemitDet->Fill(np_Eemit_MeV);
        if (np_ring >= 1 && np_ring <= 4) {
          hEemitRing[np_ring-1]->Fill(np_Eemit_MeV);
        }
      }
    }

    // ----------------- Emitted multiplicity distribution -----------------
    int maxMult = 0;
    long long sumMult = 0;
    int nEvtMult = 0;
    for (auto& kv : multByEvent) {
      if (kv.second > maxMult) maxMult = kv.second;
      sumMult += kv.second;
      ++nEvtMult;
    }
    if (maxMult < 1) maxMult = 1;

    hMultEmitted = new TH1I("hMultEmitted",
                            "Emitted neutron multiplicity per event;N emitted;Events",
                            maxMult+1, -0.5, maxMult+0.5);
    for (auto& kv : multByEvent)
      hMultEmitted->Fill(kv.second);

    if (nEvtMult > 0)
      meanEmittedMultiplicity = (double)sumMult / (double)nEvtMult;


    // ----------------- Canvas 1: multiplicities, rings, sums -----------------
    auto c1 = new TCanvas("c1", "Multiplicities", 1400, 900);
    c1->Divide(3,3);
    c1->cd(1); hMult->SetLineColor(kBlue+1);        hMult->Draw("HIST");
    c1->cd(2); hR1->SetLineColor(kRed+1);           hR1->Draw("HIST");
    c1->cd(3); hR2->SetLineColor(kGreen+2);         hR2->Draw("HIST");
    c1->cd(4); hR3->SetLineColor(kOrange+7);        hR3->Draw("HIST");
    c1->cd(5); hR4->SetLineColor(kMagenta+1);       hR4->Draw("HIST");
    c1->cd(6); hSumRings->SetLineColor(kBlue+2);    hSumRings->Draw("HIST");
    c1->cd(7); hSumRingsDead->SetLineColor(kRed+2); hSumRingsDead->Draw("HIST");
    c1->cd(8); hEce->SetLineColor(kBlue+3);         hEce->Draw("HIST");
    c1->cd(9); // vide
    c1->SaveAs(TString::Format("%s_multiplicities.png", TString(infile).ReplaceAll(".root","").Data()));

    // --- Emitted-energy canvases ---
    auto cE = new TCanvas("cE", "Emitted energy spectra", 1000, 700);
    hEemitAll->SetLineColor(kBlue+1);
    hEemitDet->SetLineColor(kRed+1);
    hEemitAll->Draw("HIST");
    hEemitDet->Draw("HIST SAME");
    {
      auto legE = new TLegend(0.6,0.75,0.88,0.88);
      legE->AddEntry(hEemitAll, "All emitted", "l");
      legE->AddEntry(hEemitDet, "Detected (ring>0)", "l");
      legE->Draw();
    }
    cE->SaveAs(TString::Format("%s_emitE.png", TString(infile).ReplaceAll(".root","").Data()));

    auto cR = new TCanvas("cR", "Emitted energy per ring", 1200, 800);
    cR->Divide(2,2);
    {
      int padR=1;
      int colsR[4] = {kRed+1, kGreen+2, kOrange+7, kMagenta+1};
      for (int r=0; r<4; ++r) {
        cR->cd(padR++);
        hEemitRing[r]->SetLineColor(colsR[r]);
        hEemitRing[r]->Draw("HIST");
      }
    }
    cR->SaveAs(TString::Format("%s_emitE_rings.png", TString(infile).ReplaceAll(".root","").Data()));

    auto cEff = new TCanvas("cEff", "Efficiency vs E_{emit}", 1000, 700);
    cEff->SetGrid();
    hEffVsE->SetLineColor(kBlue+2);
    hEffVsE->SetMarkerStyle(20);
    hEffVsE->SetMarkerColor(kBlue+2);
    hEffVsE->Draw("AP");
    cEff->SaveAs(TString::Format("%s_effVsE.png", TString(infile).ReplaceAll(".root","").Data()));

    auto cEffG = new TCanvas("cEffG", "Efficiency vs E_{emit} (gated)", 1000, 700);
    cEffG->SetGrid();
    hEffVsE_gate->SetLineColor(kRed+1);
    hEffVsE_gate->SetMarkerStyle(20);
    hEffVsE_gate->SetMarkerColor(kRed+1);
    hEffVsE_gate->Draw("AP");
    cEffG->SaveAs(TString::Format("%s_effVsE_gate.png", TString(infile).ReplaceAll(".root","").Data()));

    auto cM = new TCanvas("cM", "Emitted multiplicity", 800, 600);
    hMultEmitted->SetLineColor(kBlue+1);
    hMultEmitted->Draw("HIST");
    cM->SaveAs(TString::Format("%s_multEmitted.png", TString(infile).ReplaceAll(".root","").Data()));

    auto cS = new TCanvas("cS", "Sum emitted energy per event", 900, 600);
    hSumEemit->SetLineColor(kGreen+2);
    hSumEemit->Draw("HIST");
    auto lineB = new TLine(budgetMEV, 0, budgetMEV, hSumEemit->GetMaximum()*1.05);
    lineB->SetLineColor(kRed+1); lineB->SetLineStyle(2); lineB->Draw();
    {
      TLatex lat;
      lat.SetNDC();
      lat.SetTextSize(0.036);
      lat.DrawLatex(0.60,0.82,Form("Budget = %.2f MeV", budgetMEV));
      lat.DrawLatex(0.60,0.76,Form("Events over budget: %d", overBudgetEvents));
    }
    cS->SaveAs(TString::Format("%s_sumEemit.png", TString(infile).ReplaceAll(".root","").Data()));
  }

  // ----------------- Canvas 2: times (avg + last) -----------------
  TString tin(infile);
  TString dir  = gSystem->DirName(tin);
  TString base = gSystem->BaseName(tin);
  if (base.EndsWith(".root")) base.ReplaceAll(".root", "");
  TString outBase = Form("%s/plot_%s", dir.Data(), base.Data());

  gSystem->Exec(TString::Format("mkdir -p %s", dir.Data()));

  auto c2 = new TCanvas("c2", "Times", 1200, 500);
  c2->Divide(2,1);
  c2->cd(1);
  if (hAvgDetTime) {
    hAvgDetTime->SetLineColor(kBlue+1);
    hAvgDetTime->Draw("HIST");
  } else {
    auto t = new TLatex(0.2,0.5,"No TritonHits tree to compute average detection time");
    t->SetNDC();
    t->Draw();
  }
  c2->cd(2);
  hLast->SetLineColor(kRed+1);
  hLast->Draw("HIST");
  c2->SaveAs(TString::Format("%s_times.png", outBase.Data()));

  // ----------------- Canvas 3: neutron detection times (all hits) -----------------
  if (hDetTime) {
    auto cDet = new TCanvas("cDet", "Neutron detection times", 1000, 600);
    hDetTime->SetLineColor(kBlue+1);
    hDetTime->Draw("HIST");
    cDet->SaveAs(TString::Format("%s_neutronDetTimes.png", outBase.Data()));
  }

  // ----------------- Final efficiencies (si NeutronPrimaries présent) -----------------
  if (totEvents > 0 && nEventsAllDetBeforeGate >= 0) {
    effAllDetBeforeGate_frac    = (double)nEventsAllDetBeforeGate / (double)totEvents;
    effAllDetBeforeGate_percent = 100.0 * effAllDetBeforeGate_frac;
    double errFrac = std::sqrt(effAllDetBeforeGate_frac * (1.0 - effAllDetBeforeGate_frac)
                               / (double)totEvents);
    errEffAllDetBeforeGate_percent = 100.0 * errFrac;
  }
  if (totEvents > 0 && nEventsAnyBeforeGate >= 0) {
    effAnyBeforeGate_frac    = (double)nEventsAnyBeforeGate / (double)totEvents;
    effAnyBeforeGate_percent = 100.0 * effAnyBeforeGate_frac;
    double errFrac = std::sqrt(effAnyBeforeGate_frac * (1.0 - effAnyBeforeGate_frac)
                               / (double)totEvents);
    errEffAnyBeforeGate_percent = 100.0 * errFrac;
  }

  // ----------------- Summary stdout -----------------
  std::cout << "Summary:\n"
            << "  Total events         : " << totEvents << "\n"
            << "  Total emitted        : " << totEmitted << "\n"
            << "  Total detected (nDetected field) : " << totDetected << "\n"
            << "  Total escaped        : " << totEscaped << " (" << percentEscaped << " %)\n"
            << "  Emitted-Detected     : " << diffEmDet << "\n"
            << "  Check (Em-Det vs Esc): " << diffEmDet << " vs " << totEscaped << "\n"
            << "  Gate t <= " << tGate_ns << " ns\n"
            << "    Events with all detections in gate : " << nEventsAllDetBeforeGate
            << " (eff = " << effAllDetBeforeGate_percent
            << " +/- " << errEffAllDetBeforeGate_percent << " %)\n";

  std::cout << "    Events with any detection in gate  : ";
  if (nEventsAnyBeforeGate >= 0) {
    std::cout << nEventsAnyBeforeGate;
    if (effAnyBeforeGate_percent >= 0.0) {
      std::cout << " (eff = " << effAnyBeforeGate_percent;
      if (errEffAnyBeforeGate_percent >= 0.0)
        std::cout << " +/- " << errEffAnyBeforeGate_percent;
      std::cout << " %)";
    }
    std::cout << "\n";
  } else {
    std::cout << "N/A (no NeutronPrimaries)\n";
  }

  if (meanEmittedMultiplicity >= 0.0) {
    std::cout << "  Mean emitted multiplicity: " << meanEmittedMultiplicity << "\n";
  }

  // ----------------- Text summary file -----------------
  FILE* fp = fopen(TString::Format("%s_summary.txt", outBase.Data()), "w");
  if (fp) {
    fprintf(fp, "Input: %s\n", infile);
    fprintf(fp, "Total events: %lld\n", totEvents);
    fprintf(fp, "Total emitted: %lld\n", totEmitted);
    fprintf(fp, "Total detected (nDetected field): %lld\n", totDetected);
    fprintf(fp, "Total escaped: %lld (%.3f %%)\n", totEscaped, percentEscaped);
    fprintf(fp, "Emitted - Detected: %lld\n", diffEmDet);
    fprintf(fp, "Check (Em-Det vs Esc): %lld vs %lld\n", diffEmDet, totEscaped);
    fprintf(fp, "Gate t <= %.3f ns\n", tGate_ns);
    fprintf(fp, "Events (all detections in gate): %lld (eff = %.6f +/- %.6f %%)\n",
            nEventsAllDetBeforeGate,
            effAllDetBeforeGate_percent,
            errEffAllDetBeforeGate_percent);

    if (nEventsAnyBeforeGate >= 0 && effAnyBeforeGate_percent >= 0.0) {
      fprintf(fp, "Events (any detection in gate): %lld (eff = %.6f",
              nEventsAnyBeforeGate,
              effAnyBeforeGate_percent);
      if (errEffAnyBeforeGate_percent >= 0.0)
        fprintf(fp, " +/- %.6f", errEffAnyBeforeGate_percent);
      fprintf(fp, " %%)\n");
    } else {
      fprintf(fp, "Events (any detection in gate): N/A (no NeutronPrimaries)\n");
    }

    const double meanMultGate = (totEvents > 0)
                                ? (sumDetectedD / (double)totEvents)
                                : 0.0;
    const double meanLast = (cntLastDetect > 0)
                            ? (sumLastDetect / (double)cntLastDetect)
                            : -1.0;
    fprintf(fp, "Mean measured multiplicity in gate: %.6f\n", meanMultGate);
    fprintf(fp, "Mean average detection time [ns]: %.6f\n", meanAvgDetTime);
    fprintf(fp, "Mean last detection time [ns]: %.6f\n", meanLast);
    if (meanEmittedMultiplicity >= 0.0)
      fprintf(fp, "Mean emitted multiplicity: %.6f\n", meanEmittedMultiplicity);
    fclose(fp);
  }

  std::cout << "Saved plots to: " << outBase << "_*.png\n";

  // ----------------- Save main histograms to ROOT file -----------------
  TString rootOut = TString::Format("%s_hists.root", outBase.Data());
  std::unique_ptr<TFile> fout(TFile::Open(rootOut, "RECREATE"));
  if (fout && !fout->IsZombie()) {
    hMult->Write();
    hR1->Write(); hR2->Write(); hR3->Write(); hR4->Write();
    hSumRings->Write(); hSumRingsDead->Write();
    hEce->Write(); hLast->Write();
    if (hAvgDetTime) hAvgDetTime->Write();
    if (hDetTime)    hDetTime->Write();

    // Summary tree
    Long64_t out_totEvents  = totEvents;
    Long64_t out_totEmitted = totEmitted;
    Long64_t out_totDetected = totDetected;
    Long64_t out_totEscaped  = totEscaped;
    Double_t out_percentEscaped = percentEscaped;
    Long64_t out_diffEmDet = diffEmDet;
    Double_t out_meanMult = (totEvents > 0) ? (sumDetectedD / (double)totEvents) : 0.0;
    Double_t out_meanAvgDetTime = meanAvgDetTime;
    Double_t out_meanLastDetTime = (cntLastDetect > 0) ? (sumLastDetect / (double)cntLastDetect) : -1.0;

    auto* tSummary = new TTree("Summary", "Run summary and means");
    tSummary->Branch("totEvents", &out_totEvents, "totEvents/L");
    tSummary->Branch("totEmitted", &out_totEmitted, "totEmitted/L");
    tSummary->Branch("totDetected", &out_totDetected, "totDetected/L");
    tSummary->Branch("totEscaped", &out_totEscaped, "totEscaped/L");
    tSummary->Branch("percentEscaped", &out_percentEscaped, "percentEscaped/D");
    tSummary->Branch("diffEmDet", &out_diffEmDet, "diffEmDet/L");
    tSummary->Branch("meanMeasuredMultiplicityInGate", &out_meanMult,
                     "meanMeasuredMultiplicityInGate/D");
    tSummary->Branch("meanAvgDetectionTime_ns", &out_meanAvgDetTime,
                     "meanAvgDetectionTime_ns/D");
    tSummary->Branch("meanLastDetectionTime_ns", &out_meanLastDetTime,
                     "meanLastDetectionTime_ns/D");

    Double_t out_tGate_ns = tGate_ns;
    Long64_t out_nEventsAllDetBeforeGate = nEventsAllDetBeforeGate;
    Double_t out_effAllDetBeforeGate = effAllDetBeforeGate_percent;
    Long64_t out_nEventsAnyBeforeGate = (nEventsAnyBeforeGate>=0 ? nEventsAnyBeforeGate : -1);
    Double_t out_effAnyBeforeGate = effAnyBeforeGate_percent;
    tSummary->Branch("tGate_ns", &out_tGate_ns, "tGate_ns/D");
    tSummary->Branch("eventsAllDetectionsBeforeGate", &out_nEventsAllDetBeforeGate,
                     "eventsAllDetectionsBeforeGate/L");
    tSummary->Branch("effAllDetectionsBeforeGate_percent", &out_effAllDetBeforeGate,
                     "effAllDetectionsBeforeGate_percent/D");
    tSummary->Branch("eventsAnyDetectionBeforeGate", &out_nEventsAnyBeforeGate,
                     "eventsAnyDetectionBeforeGate/L");
    tSummary->Branch("effAnyDetectionBeforeGate_percent", &out_effAnyBeforeGate,
                     "effAnyDetectionBeforeGate_percent/D");

    Double_t out_budgetMEV = budgetMEV;
    Long64_t out_eventsOverBudget = -1; // mis à jour plus bas si tPrim
    tSummary->Branch("budgetMEV", &out_budgetMEV, "budgetMEV/D");
    tSummary->Branch("eventsOverBudget", &out_eventsOverBudget, "eventsOverBudget/L");
    tSummary->Fill();
    tSummary->Write();

    // NeutronPrimaries-derived objects
    if (hEemitAll)  hEemitAll->Write();
    if (hEemitDet)  hEemitDet->Write();
    for (int r=0; r<4; ++r)
      if (hEemitRing[r]) hEemitRing[r]->Write();
    if (hEffVsE)       hEffVsE->Write();
    if (hEffVsE_gate)  hEffVsE_gate->Write();
    if (hMultEmitted)  hMultEmitted->Write();
    if (hSumEemit)     hSumEemit->Write();

    // Budget info tree
    if (hSumEemit) {
      Double_t out_budgetMEV2 = budgetMEV;
      Long64_t out_overBudgetEvents = overBudgetEvents;
      auto* tBudget = new TTree("SummaryBudget", "Budget summary");
      tBudget->Branch("budgetMEV", &out_budgetMEV2, "budgetMEV/D");
      tBudget->Branch("eventsOverBudget", &out_overBudgetEvents, "eventsOverBudget/L");
      tBudget->Fill();
      tBudget->Write();
    }

    // Extra summary for emitted multiplicity
    if (tPrim && meanEmittedMultiplicity >= 0.0) {
      Double_t out_meanEmittedMultiplicity = meanEmittedMultiplicity;
      auto* tEmit = new TTree("SummaryEmit", "Emitted multiplicity summary");
      tEmit->Branch("meanEmittedMultiplicity",
                    &out_meanEmittedMultiplicity,
                    "meanEmittedMultiplicity/D");
      tEmit->Fill();
      tEmit->Write();
    }

    // ----------------- PARIS per-copy spectra -----------------
    if (tParis) {
      Int_t    pe_eventID = 0;
      Int_t    pe_copy = -1;
      Double_t pe_eCe = 0.0;
      Double_t pe_eNaI = 0.0;

      tParis->SetBranchAddress("eventID",  &pe_eventID);
      tParis->SetBranchAddress("copy",     &pe_copy);
      tParis->SetBranchAddress("eCe_keV",  &pe_eCe);
      tParis->SetBranchAddress("eNaI_keV", &pe_eNaI);

      std::vector<int> copies;
      copies.reserve(16);
      double maxCe = 0.0, maxNaI = 0.0;
      const Long64_t nParis = tParis->GetEntries();
      for (Long64_t i = 0; i < nParis; ++i) {
        tParis->GetEntry(i);
        if (std::find(copies.begin(), copies.end(), pe_copy) == copies.end())
          copies.push_back(pe_copy);
        if (pe_eCe  > maxCe)  maxCe  = pe_eCe;
        if (pe_eNaI > maxNaI) maxNaI = pe_eNaI;
      }
      std::sort(copies.begin(), copies.end());
      if (maxCe  <= 0) maxCe  = 5000.0;
      if (maxNaI <= 0) maxNaI = 5000.0;

      int ceBins  = std::clamp((int)std::ceil(maxCe/2.0), 100, 5000);
      int naiBins = std::clamp((int)std::ceil(maxNaI/2.0), 100, 5000);

      std::map<int, TH1D*> hCeByCopy;
      std::map<int, TH1D*> hNaIByCopy;
      for (int c : copies) {
        hCeByCopy[c]  = new TH1D(Form("hCe_copy%d", c),
                                 Form("PARIS copy %d;EdepCe [keV];Counts", c),
                                 ceBins, 0, maxCe);
        hNaIByCopy[c] = new TH1D(Form("hNaI_copy%d", c),
                                 Form("PARIS copy %d;EdepNaI [keV];Counts", c),
                                 naiBins, 0, maxNaI);
      }

      for (Long64_t i = 0; i < nParis; ++i) {
        tParis->GetEntry(i);
        if (pe_copy < 0) continue;
        if (pe_eCe  > 0) hCeByCopy[pe_copy]->Fill(pe_eCe);
        if (pe_eNaI > 0) hNaIByCopy[pe_copy]->Fill(pe_eNaI);
      }

      auto* dPARIS = fout->mkdir("PARIS");
      dPARIS->cd();
      auto* dCe = dPARIS->mkdir("Ce"); dCe->cd();
      for (auto& kv : hCeByCopy) kv.second->Write();
      dPARIS->cd();

      bool anyNaI = false;
      for (auto& kv : hNaIByCopy)
        if (kv.second->GetEntries() > 0) { anyNaI = true; break; }

      if (anyNaI) {
        auto* dNaI = dPARIS->mkdir("NaI"); dNaI->cd();
        for (auto& kv : hNaIByCopy) kv.second->Write();
      }
      fout->cd();

      // Canvases PARIS (juste pour visualisation si tu le souhaites)
      int N = (int)copies.size();
      int colsC = (N >= 9 ? 3 : (N >= 6 ? 3 : std::min(3, N > 0 ? N : 1)));
      int rowsC = (N + colsC - 1) / colsC;

      auto c3 = new TCanvas("c3", "PARIS Ce per copy", 1200, 400*rowsC);
      c3->Divide(colsC, rowsC);
      {
        int pad = 1;
        int colorBase = 600; // kBlue
        for (int c : copies) {
          c3->cd(pad++);
          hCeByCopy[c]->SetLineColor((colorBase + c) % 9 + 1);
          hCeByCopy[c]->Draw("HIST");
        }
      }
      c3->SaveAs(TString::Format("%s_paris_ce.png", outBase.Data()));

      bool anyNaI2 = false;
      for (auto& kv : hNaIByCopy)
        if (kv.second->GetEntries() > 0) { anyNaI2 = true; break; }
      if (anyNaI2) {
        auto c4 = new TCanvas("c4", "PARIS NaI per copy", 1200, 400*rowsC);
        c4->Divide(colsC, rowsC);
        int pad = 1;
        int colorBase = 600;
        for (int c : copies) {
          c4->cd(pad++);
          hNaIByCopy[c]->SetLineColor((colorBase + c) % 9 + 1);
          hNaIByCopy[c]->Draw("HIST");
        }
        c4->SaveAs(TString::Format("%s_paris_nai.png", outBase.Data()));
      }
    }

    fout->Write();
    fout->Close();
  }
}

// ----------------------------------------------------------------------
// Wrapper : parse a file list and call plotTetra for each file
// Format (txt):
//   path/to/file1.root   [tGate_ns]  [budgetMEV]
//   path/to/file2.root   [tGate_ns]  [budgetMEV]
// ----------------------------------------------------------------------
void plotTetraFromList(const char* listfile,
                       double default_tGate_ns = 150000.0,
                       double default_budgetMEV = 14.0)
{
  std::ifstream fin(listfile);
  if (!fin.is_open()) {
    std::cerr << "Cannot open list file: " << listfile << std::endl;
    return;
  }

  std::cout << "Reading file list: " << listfile << std::endl;

  std::string line;
  int lineNumber = 0;

  while (std::getline(fin, line)) {
    ++lineNumber;

    if (line.empty()) continue;
    std::string trimmed = line;
    trimmed.erase(0, trimmed.find_first_not_of(" \t\r\n"));
    if (trimmed.empty()) continue;
    if (trimmed[0] == '#' ||
        (trimmed.size() > 1 && trimmed[0] == '/' && trimmed[1] == '/'))
      continue;

    std::istringstream iss(trimmed);
    std::string infile;
    double tGate_ns = default_tGate_ns;
    double budgetMEV = default_budgetMEV;

    if (!(iss >> infile)) {
      std::cerr << "Line " << lineNumber
                << ": cannot read input file path, skipping.\n";
      continue;
    }

    if (!(iss >> tGate_ns))   tGate_ns   = default_tGate_ns;
    if (!(iss >> budgetMEV))  budgetMEV  = default_budgetMEV;

    TString tin(infile.c_str());
    TString dir  = gSystem->DirName(tin);
    TString base = gSystem->BaseName(tin);
    if (base.EndsWith(".root")) base.ReplaceAll(".root", "");
    TString outprefix = Form("%s/plot_%s", dir.Data(), base.Data());

    std::cout << "--------------------------------------------------\n";
    std::cout << "Line " << lineNumber << " :\n"
              << "  infile   = " << infile << "\n"
              << "  tGate_ns = " << tGate_ns << "\n"
              << "  budget   = " << budgetMEV << " MeV\n"
              << "  outpref  = " << outprefix
              << " (automatically used inside plotTetra)\n";

    plotTetra(infile.c_str(), outprefix.Data(), tGate_ns, budgetMEV);
  }

  std::cout << "Done processing list: " << listfile << std::endl;
}

// ============================================================
// MERGE MODE (NO TChain): read many files from filelist.txt
// EventID restarts at 0 in each file -> we build a globalEventID using per-file offset.
// Histograms are filled with ALL files combined.
// ============================================================


// --- helper: trim ---
static inline std::string TrimStr2(const std::string& s) {
  const char* ws = " \t\r\n";
  size_t b = s.find_first_not_of(ws);
  if (b == std::string::npos) return "";
  size_t e = s.find_last_not_of(ws);
  return s.substr(b, e - b + 1);
}

// --- helper: read file list into vector<string> ---
static std::vector<std::string> ReadFileList2(const char* listfile) {
  std::vector<std::string> files;
  std::ifstream fin(listfile);
  if (!fin.is_open()) {
    std::cerr << "Cannot open list file: " << listfile << std::endl;
    return files;
  }

  std::string line;
  while (std::getline(fin, line)) {
    std::string s = TrimStr2(line);
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

// --- helper: safe get tree ---
static TTree* GetTreeOrNull(TFile* f, const char* name) {
  if (!f) return nullptr;
  return dynamic_cast<TTree*>(f->Get(name));
}

// ============================================================
// MAIN: merged analysis without TChain
// Usage:
// .L plotTetra.C+
// plotTetraMergedFromList_NoChain("filelist.txt", "plots/plotTetra_merged", 150000.0, 14.0);
// ============================================================
void plotTetraMergedFromList_NoChain(const char* listfile = "filelist_11.txt",
                                     const char* outprefix = "plotTetra_merged/plotTetra_09",
                                     double tGate_ns = 150e3,
                                     double budgetMEV = 100.0)
{
  gStyle->SetOptStat(1110);

  // ----------------- read list -----------------
  std::vector<std::string> files = ReadFileList2(listfile);
  if (files.empty()) {
    std::cerr << "[ERROR] No files in list: " << listfile << "\n";
    return;
  }

  // ----------------- PASS 1: find global maxima for histogram ranges -----------------
  int    maxDet = 10, maxR1 = 10, maxR2 = 10, maxR3 = 10, maxR4 = 10;
  double maxEce = 5000.0, maxLast = 1.0;

  double maxEemit = 20.0;
  double maxSumEevt = budgetMEV * 1.5;

  long long totalEventsAllFiles = 0;

  for (size_t ifile=0; ifile<files.size(); ++ifile) {
    const std::string& fname = files[ifile];
    std::unique_ptr<TFile> f(TFile::Open(fname.c_str(), "READ"));
    if (!f || f->IsZombie()) {
      std::cerr << "[WARN] Cannot open: " << fname << " (skip)\n";
      continue;
    }

    TTree* tEvents = GetTreeOrNull(f.get(), "Events");
    TTree* tPrim   = GetTreeOrNull(f.get(), "NeutronPrimaries");

    if (!tEvents) {
      std::cerr << "[WARN] Missing 'Events' in " << fname << " (skip)\n";
      continue;
    }

    totalEventsAllFiles += tEvents->GetEntries();

    auto getMaxI = [&](const char* br) { return (int)std::max(0.0, tEvents->GetMaximum(br)); };
    auto getMaxD = [&](const char* br) { return std::max(0.0, tEvents->GetMaximum(br)); };

    maxDet = std::max(maxDet, getMaxI("nDetected"));
    maxR1  = std::max(maxR1,  getMaxI("HitsRing1"));
    maxR2  = std::max(maxR2,  getMaxI("HitsRing2"));
    maxR3  = std::max(maxR3,  getMaxI("HitsRing3"));
    maxR4  = std::max(maxR4,  getMaxI("HitsRing4"));

    maxEce  = std::max(maxEce,  getMaxD("EdepCe_keV"));
    maxLast = std::max(maxLast, getMaxD("lastNeutronTime_ns"));

    // Primaries: estimate global max Eemit and max sumE per event
    if (tPrim) {
      Int_t    np_eventID = 0;
      Double_t np_Eemit_MeV = 0.0;
      tPrim->SetBranchAddress("eventID",   &np_eventID);
      tPrim->SetBranchAddress("Eemit_MeV", &np_Eemit_MeV);

      std::map<int,double> sumEperEvt_local;
      const Long64_t nP = tPrim->GetEntries();
      for (Long64_t i=0;i<nP;++i) {
        tPrim->GetEntry(i);
        if (np_Eemit_MeV > maxEemit) maxEemit = np_Eemit_MeV;
        sumEperEvt_local[np_eventID] += np_Eemit_MeV;
      }
      for (auto& kv : sumEperEvt_local) {
        if (kv.second > maxSumEevt) maxSumEevt = kv.second;
      }
    }
  }

  if (totalEventsAllFiles <= 0) {
    std::cerr << "[ERROR] No usable 'Events' entries across files.\n";
    return;
  }

  if (maxEce <= 0)  maxEce  = 5000.0;
  if (maxLast <= 0) maxLast = 1.0;
  if (maxEemit <= 0) maxEemit = 20.0;

  int maxSumRings = maxR1 + maxR2 + maxR3 + maxR4;

  // ----------------- create histograms once (global) -----------------
  auto hMult = new TH1I(
      "hMult",
      "Measured neutron multiplicity per event in gate;N detected (t_{det} <= gate);Events",
      maxDet+1, -0.5, maxDet+0.5);

  auto hR1 = new TH1I("hR1", "Ring 1 multiplicity in gate;HitsRing1 (t_{det} <= gate);Events",
                      maxR1+1, -0.5, maxR1+0.5);
  auto hR2 = new TH1I("hR2", "Ring 2 multiplicity in gate;HitsRing2 (t_{det} <= gate);Events",
                      maxR2+1, -0.5, maxR2+0.5);
  auto hR3 = new TH1I("hR3", "Ring 3 multiplicity in gate;HitsRing3 (t_{det} <= gate);Events",
                      maxR3+1, -0.5, maxR3+0.5);
  auto hR4 = new TH1I("hR4", "Ring 4 multiplicity in gate;HitsRing4 (t_{det} <= gate);Events",
                      maxR4+1, -0.5, maxR4+0.5);

  auto hSumRings = new TH1I(
      "hSumRings",
      "Sum of hits over 4 rings in gate;HitsRing1+2+3+4 (t_{det} <= gate);Events",
      maxSumRings+1, -0.5, maxSumRings+0.5);

  auto hSumRingsDead = new TH1I(
      "hSumRings_deadtime",
      "Sum of hits over 4 rings with 400 ns deadtime per ring (t_{det} <= gate);Deadtime-corrected sum;Events",
      maxSumRings+1, -0.5, maxSumRings+0.5);

  int eBins = std::clamp((int)std::ceil(maxEce/2.0), 100, 5000);
  auto hEce = new TH1D("hEce", "Energy in Ce;EdepCe [keV];Events", eBins, 0, maxEce);

  int tBins = 200;
  auto hLast = new TH1D("hLast", "Last neutron detection time per event;time [ns];Events", tBins, 0, maxLast);

  // TritonHits merged
  TH1D* hAvgDetTime = nullptr;
  TH1D* hDetTime    = nullptr;

  // Primaries merged
  const int eBinsEmit = std::clamp((int)std::ceil(maxEemit*10.0), 100, 2000);
  TH1D* hEemitAll = nullptr;
  TH1D* hEemitDet = nullptr;
  TH1D* hEemitRing[4] = {nullptr,nullptr,nullptr,nullptr};
  TEfficiency* hEffVsE      = nullptr;
  TEfficiency* hEffVsE_gate = nullptr;
  TH1I* hMultEmitted = nullptr;
  TH1D* hSumEemit = nullptr;

  hEemitAll = new TH1D("hEemitAll", "Emitted neutron energy;E_{emit} [MeV];Counts", eBinsEmit, 0, maxEemit);
  hEemitDet = new TH1D("hEemitDet", "Detected neutrons (ring>0);E_{emit} [MeV];Counts", eBinsEmit, 0, maxEemit);
  for (int r=0;r<4;++r) {
    hEemitRing[r] = new TH1D(Form("hEemitRing%d", r+1),
                             Form("Detected neutrons in ring %d;E_{emit} [MeV];Counts", r+1),
                             eBinsEmit, 0, maxEemit);
  }
  hEffVsE = new TEfficiency("hEffVsE",
                            "Detection efficiency vs emitted energy;E_{emit} [MeV];Efficiency",
                            eBinsEmit, 0, maxEemit);
  hEffVsE_gate = new TEfficiency("hEffVsE_gate",
                                 "Detection efficiency vs E (all t_{det} in gate);E_{emit} [MeV];Efficiency",
                                 eBinsEmit, 0, maxEemit);

  double sumMaxAxis = std::max(maxSumEevt*1.1, budgetMEV*1.5);
  hSumEemit = new TH1D("hSumEemit", "Sum emitted energy per event;Sum E [MeV];Events", 120, 0, sumMaxAxis);

  // ----------------- global accumulators -----------------
  long long totEvents  = 0;
  long long totEmitted = 0;
  long long totDetected = 0;
  long long totEscaped  = 0;

  // For means
  double sumDetectedGate = 0.0;
  double sumLastDetect   = 0.0;
  long long cntLastDetect = 0;

  // For TritonHits avg per event (global event id => vector times)
  std::map<long long, std::vector<double>> timesByEvtAll_global;
  double maxHitTime_global = 0.0;

  // For Primaries: need per event structures to compute gate multiplicities, deadtime, efficiencies, sumE, mult emitted
  struct DetGateInfo { bool hasDet=false; bool anyBeforeGate=false; bool anyAfterGate=false; };
  std::map<long long, DetGateInfo> detGateByEvt;

  std::map<long long, int> multEmittedByEvt;
  std::map<long long, double> sumEperEvt_global;

  std::map<long long, int> multGateByEvt;
  std::map<long long, std::array<int,4>> hitsRingGateByEvt;
  std::map<long long, std::array<std::vector<double>,4>> timesByEvtPerRing;

  // store last time per event (from Events)
  std::map<long long, double> lastByEvt;

  // ----------------- PASS 2: loop all files, build globalEventID with offsets -----------------
  long long eventOffset = 0; // increases per file to avoid collisions

  for (size_t ifile=0; ifile<files.size(); ++ifile) {
    const std::string& fname = files[ifile];
    std::unique_ptr<TFile> f(TFile::Open(fname.c_str(), "READ"));
    if (!f || f->IsZombie()) {
      std::cerr << "[WARN] Cannot open: " << fname << " (skip)\n";
      continue;
    }

    TTree* tEvents = GetTreeOrNull(f.get(), "Events");
    TTree* tHits   = GetTreeOrNull(f.get(), "TritonHits");
    TTree* tPrim   = GetTreeOrNull(f.get(), "NeutronPrimaries");

    if (!tEvents) {
      std::cerr << "[WARN] Missing 'Events' in " << fname << " (skip)\n";
      continue;
    }

    // Determine this file offset increment (use max EventID in file)
    double maxEvtIdD = tEvents->GetMaximum("EventID");
    long long maxEvtId = (maxEvtIdD > 0 ? (long long)std::llround(maxEvtIdD) : 0LL);
    long long thisFileSpan = maxEvtId + 1; // safe even if some IDs missing
    if (thisFileSpan < 1) thisFileSpan = (long long)tEvents->GetEntries(); // fallback
    if (thisFileSpan < 1) thisFileSpan = 1;

    // --- Events branches ---
    Int_t    EventID = 0;
    Int_t    NNeutronsEmitted = 0;
    Double_t nDetected = 0.0;
    Double_t EdepCe_keV = 0.0;
    Int_t    NNeutronsEscaped = 0;
    Double_t lastNeutronTime_ns = -1.0;

    tEvents->SetBranchAddress("EventID", &EventID);
    tEvents->SetBranchAddress("NNeutronsEmitted", &NNeutronsEmitted);
    tEvents->SetBranchAddress("nDetected", &nDetected);
    tEvents->SetBranchAddress("EdepCe_keV", &EdepCe_keV);
    tEvents->SetBranchAddress("NNeutronsEscaped", &NNeutronsEscaped);
    tEvents->SetBranchAddress("lastNeutronTime_ns", &lastNeutronTime_ns);

    const Long64_t nEv = tEvents->GetEntries();
    totEvents += nEv;

    for (Long64_t i=0;i<nEv;++i) {
      tEvents->GetEntry(i);

      const long long gID = eventOffset + (long long)EventID;

      if (EdepCe_keV > 0) hEce->Fill(EdepCe_keV);

      if (lastNeutronTime_ns >= 0) {
        hLast->Fill(lastNeutronTime_ns);
        sumLastDetect += lastNeutronTime_ns;
        ++cntLastDetect;
      }

      lastByEvt[gID] = lastNeutronTime_ns;

      totEmitted  += NNeutronsEmitted;
      totDetected += (long long)std::llround(nDetected);
      totEscaped  += NNeutronsEscaped;
    }

    // --- TritonHits (store per global event times) ---
    if (tHits) {
      Int_t    eID  = 0;
      Double_t t_ns = 0.0;
      tHits->SetBranchAddress("EventID", &eID);
      tHits->SetBranchAddress("time_ns", &t_ns);

      const Long64_t nH = tHits->GetEntries();
      for (Long64_t i=0;i<nH;++i) {
        tHits->GetEntry(i);
        if (t_ns < 0) continue;
        const long long gID = eventOffset + (long long)eID;
        timesByEvtAll_global[gID].push_back(t_ns);
        if (t_ns > maxHitTime_global) maxHitTime_global = t_ns;
      }
    }

    // --- NeutronPrimaries (build all per-event structures using globalEventID) ---
    if (tPrim) {
      Int_t    np_eventID       = 0;
      Int_t    np_ring          = 0;
      Double_t np_Eemit_MeV     = 0.0;
      Double_t np_detectTime_ns = -1.0;

      tPrim->SetBranchAddress("eventID",       &np_eventID);
      tPrim->SetBranchAddress("ring",          &np_ring);
      tPrim->SetBranchAddress("Eemit_MeV",     &np_Eemit_MeV);
      tPrim->SetBranchAddress("detectTime_ns", &np_detectTime_ns);

      const Long64_t nP = tPrim->GetEntries();
      for (Long64_t i=0;i<nP;++i) {
        tPrim->GetEntry(i);

        const long long gID = eventOffset + (long long)np_eventID;

        // emitted multiplicity + sum energy per event
        multEmittedByEvt[gID] += 1;
        sumEperEvt_global[gID] += np_Eemit_MeV;

        // efficiency per neutron
        hEemitAll->Fill(np_Eemit_MeV);
        const bool detected = (np_ring > 0);
        hEffVsE->Fill(detected, np_Eemit_MeV);

        if (detected) {
          hEemitDet->Fill(np_Eemit_MeV);
          if (np_ring>=1 && np_ring<=4) hEemitRing[np_ring-1]->Fill(np_Eemit_MeV);
        }

        // gate bookkeeping (only if detected and time valid)
        if (np_ring > 0 && np_ring <= 4 && np_detectTime_ns >= 0.0) {
          auto &info = detGateByEvt[gID];
          info.hasDet = true;
          if (np_detectTime_ns <= tGate_ns) info.anyBeforeGate = true;
          else                              info.anyAfterGate  = true;

          if (np_detectTime_ns <= tGate_ns) {
            multGateByEvt[gID]++;

            auto &arr = hitsRingGateByEvt[gID]; // default init {0,0,0,0}
            arr[np_ring-1]++;

            timesByEvtPerRing[gID][np_ring-1].push_back(np_detectTime_ns);
          }
        }
      }
    }

    // next file offset
    eventOffset += thisFileSpan;
  }

  // ----------------- Build TritonHits histos from collected times -----------------
  if (!timesByEvtAll_global.empty()) {
    // avg per event distribution
    double maxAvg = 0.0;
    std::vector<double> avgs;
    avgs.reserve(timesByEvtAll_global.size());

    for (auto &kv : timesByEvtAll_global) {
      auto &v = kv.second;
      if (v.empty()) continue;
      double s = 0.0;
      for (double x : v) s += x;
      double m = s / (double)v.size();
      avgs.push_back(m);
      if (m > maxAvg) maxAvg = m;
    }

    if (maxAvg <= 0) maxAvg = 1.0;
    hAvgDetTime = new TH1D("hAvgDetTime",
                           "Average neutron detection time per event;time [ns];Events",
                           200, 0, maxAvg);

    for (double m : avgs) hAvgDetTime->Fill(m);

    // all hits time spectrum
    if (maxHitTime_global <= 0) maxHitTime_global = 1.0;
    int nBinsTime = (int)std::ceil(maxHitTime_global) + 1;
    hDetTime = new TH1D("hDetTime",
                        "Neutron detection time (all TritonHits);time [ns];Counts",
                        nBinsTime, 0, (double)nBinsTime);

    for (auto &kv : timesByEvtAll_global) {
      for (double t : kv.second) hDetTime->Fill(t);
    }
  }

  // ----------------- Deadtime + fill multiplicity histos from maps -----------------
  const double deadtime_ns = 720.0;

  // deadtime corrected sum per event
  std::map<long long,int> sumDeadPerEvt;

  for (auto &kv : timesByEvtPerRing) {
    const long long gID = kv.first;
    auto &perRing = kv.second;

    int sumDeadEvt = 0;
    for (int r=0;r<4;++r) {
      auto &v = perRing[r];
      if (v.empty()) continue;
      std::sort(v.begin(), v.end());

      double lastCounted = -1e30;
      int countRing = 0;
      for (double t : v) {
        if (t < 0) continue;
        if (countRing==0 || (t - lastCounted) >= deadtime_ns) {
          ++countRing;
          lastCounted = t;
        }
      }
      sumDeadEvt += countRing;
    }
    sumDeadPerEvt[gID] = sumDeadEvt;
  }

  // Fill per-event gate multiplicities using the keys we have (detGateByEvt gives us events with detections)
  for (auto &kv : detGateByEvt) {
    const long long gID = kv.first;

    // deadtime corrected
    int sumDead = 0;
    auto itD = sumDeadPerEvt.find(gID);
    if (itD != sumDeadPerEvt.end()) sumDead = itD->second;
    hSumRingsDead->Fill(sumDead);

    // mult in gate
    int multG = 0;
    auto itM = multGateByEvt.find(gID);
    if (itM != multGateByEvt.end()) multG = itM->second;
    hMult->Fill(multG);
    sumDetectedGate += multG;

    // ring hits in gate
    std::array<int,4> rings = {0,0,0,0};
    auto itR = hitsRingGateByEvt.find(gID);
    if (itR != hitsRingGateByEvt.end()) rings = itR->second;

    hR1->Fill(rings[0]);
    hR2->Fill(rings[1]);
    hR3->Fill(rings[2]);
    hR4->Fill(rings[3]);

    hSumRings->Fill(rings[0]+rings[1]+rings[2]+rings[3]);
  }

  // ----------------- Efficiencies per event (all det in gate / any in gate) -----------------
  long long nEventsAllDetBeforeGate = 0;
  long long nEventsAnyBeforeGate    = 0;

  for (auto &kv : detGateByEvt) {
    const DetGateInfo &info = kv.second;
    if (!info.hasDet) continue;
    if (info.anyBeforeGate) ++nEventsAnyBeforeGate;
    if (!info.anyAfterGate && info.anyBeforeGate) ++nEventsAllDetBeforeGate;
  }

  // Fill gated efficiency per neutron (needs eventAllDetInGate)
  // We can approximate by re-looping on sumEperEvt_global keys is not enough (we need neutrons).
  // -> Best: do it during primaries loop, but we delayed. Here: we cannot reconstruct neutron-level info without re-reading files.
  // So: simplest & correct: re-read primaries one more time ONLY to fill hEffVsE_gate with eventAllDetInGate.

  // Re-read primaries for gated eff (third light pass)
  {
    long long eventOffset2 = 0;
    for (size_t ifile=0; ifile<files.size(); ++ifile) {
      const std::string& fname = files[ifile];
      std::unique_ptr<TFile> f(TFile::Open(fname.c_str(), "READ"));
      if (!f || f->IsZombie()) continue;

      TTree* tEvents = GetTreeOrNull(f.get(), "Events");
      TTree* tPrim   = GetTreeOrNull(f.get(), "NeutronPrimaries");
      if (!tEvents || !tPrim) {
        // still need offset increment
        if (tEvents) {
          double maxEvtIdD = tEvents->GetMaximum("EventID");
          long long maxEvtId = (maxEvtIdD > 0 ? (long long)std::llround(maxEvtIdD) : 0LL);
          long long span = maxEvtId + 1;
          if (span < 1) span = (long long)tEvents->GetEntries();
          if (span < 1) span = 1;
          eventOffset2 += span;
        }
        continue;
      }

      double maxEvtIdD = tEvents->GetMaximum("EventID");
      long long maxEvtId = (maxEvtIdD > 0 ? (long long)std::llround(maxEvtIdD) : 0LL);
      long long span = maxEvtId + 1;
      if (span < 1) span = (long long)tEvents->GetEntries();
      if (span < 1) span = 1;

      Int_t    np_eventID       = 0;
      Int_t    np_ring          = 0;
      Double_t np_Eemit_MeV     = 0.0;

      tPrim->SetBranchAddress("eventID",   &np_eventID);
      tPrim->SetBranchAddress("ring",      &np_ring);
      tPrim->SetBranchAddress("Eemit_MeV", &np_Eemit_MeV);

      const Long64_t nP = tPrim->GetEntries();
      for (Long64_t i=0;i<nP;++i) {
        tPrim->GetEntry(i);
        const long long gID = eventOffset2 + (long long)np_eventID;

        bool eventAllDetInGate = false;
        auto itInfo = detGateByEvt.find(gID);
        if (itInfo != detGateByEvt.end()) {
          const DetGateInfo &inf = itInfo->second;
          eventAllDetInGate = (inf.hasDet && inf.anyBeforeGate && !inf.anyAfterGate);
        }
        const bool detected = (np_ring > 0);
        hEffVsE_gate->Fill(detected && eventAllDetInGate, np_Eemit_MeV);
      }

      eventOffset2 += span;
    }
  }

  // ----------------- Sum emitted energy per event + over budget -----------------
  int overBudgetEvents = 0;
  for (auto &kv : sumEperEvt_global) {
    hSumEemit->Fill(kv.second);
    if (kv.second > budgetMEV) ++overBudgetEvents;
  }

  // ----------------- Emitted multiplicity distribution -----------------
  int maxMult = 1;
  long long sumMult = 0;
  long long nEvtMult = 0;
  for (auto &kv : multEmittedByEvt) {
    maxMult = std::max(maxMult, kv.second);
    sumMult += kv.second;
    ++nEvtMult;
  }
  hMultEmitted = new TH1I("hMultEmitted",
                          "Emitted neutron multiplicity per event;N emitted;Events",
                          maxMult+1, -0.5, maxMult+0.5);
  for (auto &kv : multEmittedByEvt) hMultEmitted->Fill(kv.second);

  double meanEmittedMultiplicity = (nEvtMult>0 ? (double)sumMult/(double)nEvtMult : -1.0);

  // ----------------- Output directory -----------------
  TString outp(outprefix);
  TString outDir = gSystem->DirName(outp);
  if (outDir.Length() > 0) gSystem->Exec(TString::Format("mkdir -p %s", outDir.Data()));

  // ----------------- CANVASES -----------------
  {
    auto c1 = new TCanvas("c1_merged", "Multiplicities (merged)", 1400, 900);
    c1->Divide(3,3);
    c1->cd(1); hMult->SetLineColor(kBlue+1);        hMult->Draw("HIST");
    c1->cd(2); hR1->SetLineColor(kRed+1);           hR1->Draw("HIST");
    c1->cd(3); hR2->SetLineColor(kGreen+2);         hR2->Draw("HIST");
    c1->cd(4); hR3->SetLineColor(kOrange+7);        hR3->Draw("HIST");
    c1->cd(5); hR4->SetLineColor(kMagenta+1);       hR4->Draw("HIST");
    c1->cd(6); hSumRings->SetLineColor(kBlue+2);    hSumRings->Draw("HIST");
    c1->cd(7); hSumRingsDead->SetLineColor(kRed+2); hSumRingsDead->Draw("HIST");
    c1->cd(8); hEce->SetLineColor(kBlue+3);         hEce->Draw("HIST");
    c1->cd(9); hLast->SetLineColor(kRed+1);         hLast->Draw("HIST");
    c1->SaveAs(Form("%s_multiplicities.png", outp.Data()));
  }

  {
    auto c2 = new TCanvas("c2_merged", "Times (merged)", 1200, 500);
    c2->Divide(2,1);
    c2->cd(1);
    if (hAvgDetTime) { hAvgDetTime->SetLineColor(kBlue+1); hAvgDetTime->Draw("HIST"); }
    else {
      auto t = new TLatex(0.2,0.5,"No TritonHits information (merged)");
      t->SetNDC(); t->Draw();
    }
    c2->cd(2);
    hLast->SetLineColor(kRed+1); hLast->Draw("HIST");
    c2->SaveAs(Form("%s_times.png", outp.Data()));
  }

  if (hDetTime) {
    auto cDet = new TCanvas("cDet_merged", "Neutron detection times (merged)", 1000, 600);
    hDetTime->SetLineColor(kBlue+1);
    hDetTime->Draw("HIST");
    cDet->SaveAs(Form("%s_neutronDetTimes.png", outp.Data()));
  }

  {
    auto cE = new TCanvas("cE_merged", "Emitted energy spectra (merged)", 1000, 700);
    hEemitAll->SetLineColor(kBlue+1);
    hEemitDet->SetLineColor(kRed+1);
    hEemitAll->Draw("HIST");
    hEemitDet->Draw("HIST SAME");
    auto legE = new TLegend(0.6,0.75,0.88,0.88);
    legE->AddEntry(hEemitAll, "All emitted", "l");
    legE->AddEntry(hEemitDet, "Detected (ring>0)", "l");
    legE->Draw();
    cE->SaveAs(Form("%s_emitE.png", outp.Data()));

    auto cEff = new TCanvas("cEff_merged", "Efficiency vs Eemit (merged)", 1000, 700);
    cEff->SetGrid();
    hEffVsE->SetLineColor(kBlue+2);
    hEffVsE->SetMarkerStyle(20);
    hEffVsE->SetMarkerColor(kBlue+2);
    hEffVsE->Draw("AP");
    cEff->SaveAs(Form("%s_effVsE.png", outp.Data()));

    auto cEffG = new TCanvas("cEffG_merged", "Efficiency vs Eemit gated (merged)", 1000, 700);
    cEffG->SetGrid();
    hEffVsE_gate->SetLineColor(kRed+1);
    hEffVsE_gate->SetMarkerStyle(20);
    hEffVsE_gate->SetMarkerColor(kRed+1);
    hEffVsE_gate->Draw("AP");
    cEffG->SaveAs(Form("%s_effVsE_gate.png", outp.Data()));

    auto cM = new TCanvas("cM_merged", "Emitted multiplicity (merged)", 800, 600);
    hMultEmitted->SetLineColor(kBlue+1);
    hMultEmitted->Draw("HIST");
    cM->SaveAs(Form("%s_multEmitted.png", outp.Data()));

    auto cS = new TCanvas("cS_merged", "Sum emitted energy per event (merged)", 900, 600);
    hSumEemit->SetLineColor(kGreen+2);
    hSumEemit->Draw("HIST");
    auto lineB = new TLine(budgetMEV, 0, budgetMEV, hSumEemit->GetMaximum()*1.05);
    lineB->SetLineColor(kRed+1); lineB->SetLineStyle(2); lineB->Draw();
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.036);
    lat.DrawLatex(0.60,0.82,Form("Budget = %.2f MeV", budgetMEV));
    lat.DrawLatex(0.60,0.76,Form("Events over budget: %d", overBudgetEvents));
    cS->SaveAs(Form("%s_sumEemit.png", outp.Data()));
  }

  // ----------------- Summary -----------------
  double percentEscaped = (totEmitted > 0) ? 100.0 * (double)totEscaped / (double)totEmitted : 0.0;
  long long diffEmDet = totEmitted - totDetected;

  double meanMultGate = (totEvents > 0) ? (sumDetectedGate / (double)totEvents) : 0.0;
  double meanLast = (cntLastDetect > 0) ? (sumLastDetect / (double)cntLastDetect) : -1.0;

  double effAllDetBeforeGate_frac = (totEvents>0) ? (double)nEventsAllDetBeforeGate/(double)totEvents : 0.0;
  double effAllDetBeforeGate_percent = 100.0 * effAllDetBeforeGate_frac;
  double errEffAllDetBeforeGate_percent = (totEvents>0)
    ? 100.0 * std::sqrt(effAllDetBeforeGate_frac*(1.0-effAllDetBeforeGate_frac)/(double)totEvents)
    : 0.0;

  double effAnyBeforeGate_frac = (totEvents>0) ? (double)nEventsAnyBeforeGate/(double)totEvents : 0.0;
  double effAnyBeforeGate_percent = 100.0 * effAnyBeforeGate_frac;
  double errEffAnyBeforeGate_percent = (totEvents>0)
    ? 100.0 * std::sqrt(effAnyBeforeGate_frac*(1.0-effAnyBeforeGate_frac)/(double)totEvents)
    : 0.0;

  std::cout << "\n=== MERGED (NO CHAIN) SUMMARY ===\n"
            << "Files in list: " << files.size() << "\n"
            << "Total events         : " << totEvents << "\n"
            << "Total emitted        : " << totEmitted << "\n"
            << "Total detected(field): " << totDetected << "\n"
            << "Total escaped        : " << totEscaped << " (" << percentEscaped << " %)\n"
            << "Emitted-Detected     : " << diffEmDet << "\n"
            << "Gate t <= " << tGate_ns << " ns\n"
            << "Events all det in gate : " << nEventsAllDetBeforeGate
            << " (eff=" << effAllDetBeforeGate_percent << " +/- " << errEffAllDetBeforeGate_percent << " %)\n"
            << "Events any det in gate : " << nEventsAnyBeforeGate
            << " (eff=" << effAnyBeforeGate_percent << " +/- " << errEffAnyBeforeGate_percent << " %)\n"
            << "Mean measured mult in gate: " << meanMultGate << "\n"
            << "Mean last det time [ns]   : " << meanLast << "\n"
            << "Mean emitted multiplicity : " << meanEmittedMultiplicity << "\n";

  // ----------------- Save ROOT output -----------------
  TString rootOut = Form("%s_hists.root", outp.Data());
  std::unique_ptr<TFile> fout(TFile::Open(rootOut, "RECREATE"));
  if (!fout || fout->IsZombie()) {
    std::cerr << "[ERROR] Cannot create output ROOT file: " << rootOut << "\n";
    return;
  }

  hMult->Write();
  hR1->Write(); hR2->Write(); hR3->Write(); hR4->Write();
  hSumRings->Write(); hSumRingsDead->Write();
  hEce->Write(); hLast->Write();
  if (hAvgDetTime) hAvgDetTime->Write();
  if (hDetTime)    hDetTime->Write();

  hEemitAll->Write();
  hEemitDet->Write();
  for (int r=0;r<4;++r) hEemitRing[r]->Write();
  hEffVsE->Write();
  hEffVsE_gate->Write();
  hMultEmitted->Write();
  hSumEemit->Write();

  // summary tree
  auto* tSummary = new TTree("Summary", "Merged summary (no chain)");
  Long64_t out_totEvents = totEvents, out_totEmitted = totEmitted, out_totDetected = totDetected, out_totEscaped = totEscaped;
  Double_t out_percentEscaped = percentEscaped;
  Long64_t out_diffEmDet = diffEmDet;
  Double_t out_tGate = tGate_ns;
  Double_t out_meanMultGate = meanMultGate;
  Double_t out_meanLast = meanLast;
  Double_t out_meanEmit = meanEmittedMultiplicity;

  Long64_t out_allGate = nEventsAllDetBeforeGate;
  Long64_t out_anyGate = nEventsAnyBeforeGate;
  Double_t out_effAll = effAllDetBeforeGate_percent;
  Double_t out_errAll = errEffAllDetBeforeGate_percent;
  Double_t out_effAny = effAnyBeforeGate_percent;
  Double_t out_errAny = errEffAnyBeforeGate_percent;

  Double_t out_budget = budgetMEV;
  Long64_t out_overBudget = overBudgetEvents;

  tSummary->Branch("totEvents", &out_totEvents, "totEvents/L");
  tSummary->Branch("totEmitted", &out_totEmitted, "totEmitted/L");
  tSummary->Branch("totDetected", &out_totDetected, "totDetected/L");
  tSummary->Branch("totEscaped", &out_totEscaped, "totEscaped/L");
  tSummary->Branch("percentEscaped", &out_percentEscaped, "percentEscaped/D");
  tSummary->Branch("diffEmDet", &out_diffEmDet, "diffEmDet/L");
  tSummary->Branch("tGate_ns", &out_tGate, "tGate_ns/D");
  tSummary->Branch("meanMeasuredMultiplicityInGate", &out_meanMultGate, "meanMeasuredMultiplicityInGate/D");
  tSummary->Branch("meanLastDetectionTime_ns", &out_meanLast, "meanLastDetectionTime_ns/D");
  tSummary->Branch("meanEmittedMultiplicity", &out_meanEmit, "meanEmittedMultiplicity/D");

  tSummary->Branch("eventsAllDetectionsBeforeGate", &out_allGate, "eventsAllDetectionsBeforeGate/L");
  tSummary->Branch("effAllDetectionsBeforeGate_percent", &out_effAll, "effAllDetectionsBeforeGate_percent/D");
  tSummary->Branch("errEffAllDetectionsBeforeGate_percent", &out_errAll, "errEffAllDetectionsBeforeGate_percent/D");

  tSummary->Branch("eventsAnyDetectionBeforeGate", &out_anyGate, "eventsAnyDetectionBeforeGate/L");
  tSummary->Branch("effAnyDetectionBeforeGate_percent", &out_effAny, "effAnyDetectionBeforeGate_percent/D");
  tSummary->Branch("errEffAnyDetectionBeforeGate_percent", &out_errAny, "errEffAnyDetectionBeforeGate_percent/D");

  tSummary->Branch("budgetMEV", &out_budget, "budgetMEV/D");
  tSummary->Branch("eventsOverBudget", &out_overBudget, "eventsOverBudget/L");

  tSummary->Fill();
  tSummary->Write();

  fout->Close();

  std::cout << "[INFO] Merged outputs saved:\n"
            << "  PNG:  " << outp << "_*.png\n"
            << "  ROOT: " << rootOut << "\n";
}