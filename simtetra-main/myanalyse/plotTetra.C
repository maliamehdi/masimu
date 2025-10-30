#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TString.h>
#include <TROOT.h>
#include <TLatex.h>
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>

// Usage in ROOT:
// .L plotTetra.C
// plotTetra("../myanalyse/output_cf252_neutron_primarygenerator.root")
// plotTetra("path/to/output.root", "../myanalyse/plots/plotTetra")

void plotTetra(const char* infile = "output_neutron_run0_smeared.root",
               const char* outprefix = "../myanalyse/plots/150usplotTetra_analysis",
               double tGate_ns = 150000.0)
{
  gStyle->SetOptStat(1110);

  // Open file and get trees
  std::unique_ptr<TFile> f(TFile::Open(infile, "READ"));
  if (!f || f->IsZombie()) { std::cerr << "Cannot open file: " << infile << std::endl; return; }

  auto* tEvents = dynamic_cast<TTree*>(f->Get("Events"));
  auto* tHits   = dynamic_cast<TTree*>(f->Get("TritonHits"));
  auto* tParis  = dynamic_cast<TTree*>(f->Get("ParisEdep"));
  if (!tEvents) { std::cerr << "Missing tree 'Events' in file" << std::endl; return; }
  if (!tHits)   { std::cerr << "Missing tree 'TritonHits' in file (needed for avg detection time)" << std::endl; }

  // Branch setup for Events
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

  // Determine sensible histogram ranges from data
  auto getMaxI = [&](const char* br) { return (int)std::max(0.0, tEvents->GetMaximum(br)); };
  auto getMaxD = [&](const char* br) { return std::max(0.0,   tEvents->GetMaximum(br)); };

  int maxDet = std::max(10, getMaxI("nDetected"));
  int maxR1  = std::max(10, getMaxI("HitsRing1"));
  int maxR2  = std::max(10, getMaxI("HitsRing2"));
  int maxR3  = std::max(10, getMaxI("HitsRing3"));
  int maxR4  = std::max(10, getMaxI("HitsRing4"));

  double maxEce = getMaxD("EdepCe_keV"); if (maxEce <= 0) maxEce = 5000.0; // keV
  double maxLast = getMaxD("lastNeutronTime_ns"); if (maxLast <= 0) maxLast = 1.0; // ns

  // Histograms
  auto hMult = new TH1I("hMult", "Measured neutron multiplicity per event;N detected;Events", maxDet+1, -0.5, maxDet+0.5);
  auto hR1   = new TH1I("hR1",   "Ring 1 multiplicity;HitsRing1;Events", maxR1+1, -0.5, maxR1+0.5);
  auto hR2   = new TH1I("hR2",   "Ring 2 multiplicity;HitsRing2;Events", maxR2+1, -0.5, maxR2+0.5);
  auto hR3   = new TH1I("hR3",   "Ring 3 multiplicity;HitsRing3;Events", maxR3+1, -0.5, maxR3+0.5);
  auto hR4   = new TH1I("hR4",   "Ring 4 multiplicity;HitsRing4;Events", maxR4+1, -0.5, maxR4+0.5);

  int eBins = std::clamp((int)std::ceil(maxEce/2.0), 100, 5000);
  auto hEce = new TH1D("hEce",   "Energy in Ce;EdepCe [keV];Events", eBins, 0, maxEce);

  int tBins = 200;
  auto hLast = new TH1D("hLast", "Last neutron detection time;time [ns];Events", tBins, 0, maxLast);

  // Totals
  long long totEvents = tEvents->GetEntries();
  long long totEmitted = 0, totDetected = 0, totEscaped = 0;
  // For means
  double sumDetectedD = 0.0; // sum of measured multiplicity (as double)
  double sumLastDetect = 0.0; int cntLastDetect = 0;
  // Time-gated efficiency counters
  long long nEventsAllDetBeforeGate = 0; // event counted if all detections are before the time gate (using last time)
  long long nEventsAnyBeforeGate = -1;   // will be computed from TritonHits (min time per event)

  // Loop over Events tree
  for (Long64_t i = 0; i < totEvents; ++i) {
    tEvents->GetEntry(i);

  // multiplicity (gate: fill only if detections occur before tGate_ns using lastNeutronTime_ns)
  if (nDetected == 0 || ( nDetected > 0 && lastNeutronTime_ns >= 0 && lastNeutronTime_ns <= tGate_ns)) {
    hMult->Fill((int)std::llround(nDetected));
    sumDetectedD += nDetected;
  }
  // keep overall mean multiplicity on all events (not gated)
  
    hR1->Fill(HitsRing1); hR2->Fill(HitsRing2); hR3->Fill(HitsRing3); hR4->Fill(HitsRing4);

    // Ce spectrum (ignore zeros)
    if (EdepCe_keV > 0) hEce->Fill(EdepCe_keV);

    // last detection time (ignore -1/no-detection)
  if (lastNeutronTime_ns >= 0) { hLast->Fill(lastNeutronTime_ns); sumLastDetect += lastNeutronTime_ns; ++cntLastDetect; }
    // Count events where all detections are before the gate (using last time)
    if (nDetected > 0 && lastNeutronTime_ns >= 0 && lastNeutronTime_ns <= tGate_ns) {
      ++nEventsAllDetBeforeGate;
    }

    // totals
    totEmitted  += NNeutronsEmitted;
    totDetected += (long long) std::llround(nDetected);
    totEscaped  += NNeutronsEscaped;
  }

  double percentEscaped = (totEmitted > 0) ? (100.0 * (double)totEscaped / (double)totEmitted) : 0.0;
  long long diffEmDet = totEmitted - totDetected;
  const double effAllDetBeforeGate = (totEvents > 0) ? (100.0 * (double)nEventsAllDetBeforeGate / (double)totEvents) : 0.0;
  const double effAnyBeforeGate = (tHits && totEvents > 0 && nEventsAnyBeforeGate >= 0)
                                  ? (100.0 * (double)nEventsAnyBeforeGate / (double)totEvents) : -1.0;

  std::cout << "Summary:\n"
            << "  Total events         : " << totEvents << "\n"
            << "  Total emitted        : " << totEmitted << "\n"
            << "  Total detected       : " << nEventsAllDetBeforeGate <<"\n"//totDetected << "\n"
            << "  Total escaped        : " << totEscaped << " (" << percentEscaped << " %)\n"
            << "  Emitted-Detected     : " << diffEmDet << "\n"
            << "  Check (Em-Det vs Esc): " << diffEmDet << " vs " << totEscaped << "\n"
            << "  Gate t <= " << tGate_ns << " ns\n"
            << "    Events with all detections <= gate : " << nEventsAllDetBeforeGate
            << " (eff = " << effAllDetBeforeGate << " %)\n"
            << "    Events with any detection <= gate  : "
            << ( (nEventsAnyBeforeGate>=0) ? std::to_string(nEventsAnyBeforeGate) : std::string("N/A (no TritonHits)") )
            << ( (effAnyBeforeGate>=0) ? (std::string(" (eff = ") + std::to_string(effAnyBeforeGate) + " %)") : std::string("") )
            << "\n";

  // Average detection time per event from TritonHits
  TH1D* hAvgDetTime = nullptr;
  double meanAvgDetTime = -1.0; // mean of per-event avg detection time
  if (tHits) {
    Int_t    eID = 0; Double_t t_ns = 0.0;
    tHits->SetBranchAddress("EventID", &eID);
    tHits->SetBranchAddress("time_ns", &t_ns);

    std::map<int, std::vector<double>> timesByEvt;
    double maxAvg = 0.0;

    const Long64_t nHits = tHits->GetEntries();
    for (Long64_t i = 0; i < nHits; ++i) {
      tHits->GetEntry(i);
      if (t_ns >= 0) timesByEvt[eID].push_back(t_ns);
    }
    // Compute averages and find max
    std::vector<double> avgs; avgs.reserve(timesByEvt.size());
    long long anyBefore = 0;
    for (auto& kv : timesByEvt) {
      const auto& v = kv.second;
      if (v.empty()) continue;
      double s = 0; for (double x : v) s += x;
      double m = s / (double)v.size();
      avgs.push_back(m);
      if (m > maxAvg) maxAvg = m;
      // any detection before gate -> min time criterion
      double vmin = *std::min_element(v.begin(), v.end());
      if (vmin <= tGate_ns) ++anyBefore;
    }
    nEventsAnyBeforeGate = anyBefore;
    if (maxAvg <= 0) maxAvg = 1.0;
    hAvgDetTime = new TH1D("hAvgDetTime", "Average neutron detection time per event;time [ns];Events", 200, 0, maxAvg);
    double sumAvg = 0.0; int cntAvg = 0;
    for (double m : avgs) { hAvgDetTime->Fill(m); sumAvg += m; ++cntAvg; }
    if (cntAvg > 0) meanAvgDetTime = sumAvg / (double)cntAvg;
  }

  // Draw canvases
  TString outBase(outprefix);
  gSystem->Exec(TString::Format("mkdir -p %s", gSystem->DirName(outBase)));

  // Canvas 1: multiplicities
  auto c1 = new TCanvas("c1", "Multiplicities", 1200, 900);
  c1->Divide(2,3);
  c1->cd(1); hMult->SetLineColor(kBlue+1); hMult->Draw("HIST");
  c1->cd(2); hR1->SetLineColor(kRed+1);   hR1->Draw("HIST");
  c1->cd(3); hR2->SetLineColor(kGreen+2); hR2->Draw("HIST");
  c1->cd(4); hR3->SetLineColor(kOrange+7);hR3->Draw("HIST");
  c1->cd(5); hR4->SetLineColor(kMagenta+1);hR4->Draw("HIST");
  c1->cd(6); hEce->SetLineColor(kBlue+2); hEce->Draw("HIST");
  c1->SaveAs(TString::Format("%s_multiplicities.png", outBase.Data()));

  // Canvas 2: times
  auto c2 = new TCanvas("c2", "Times", 1200, 500);
  c2->Divide(2,1);
  c2->cd(1); if (hAvgDetTime) { hAvgDetTime->SetLineColor(kBlue+1); hAvgDetTime->Draw("HIST"); }
             else { auto t = new TLatex(0.2,0.5,"No TritonHits tree to compute average detection time"); t->SetNDC(); t->Draw(); }
  c2->cd(2); hLast->SetLineColor(kRed+1); hLast->Draw("HIST");
  c2->SaveAs(TString::Format("%s_times.png", outBase.Data()));

  // Also write a small text summary
  FILE* fp = fopen(TString::Format("%s_summary.txt", outBase.Data()), "w");
  if (fp) {
    fprintf(fp, "Input: %s\n", infile);
    fprintf(fp, "Total events: %lld\n", totEvents);
    fprintf(fp, "Total emitted: %lld\n", totEmitted);
    fprintf(fp, "Total detected: %lld\n", totDetected);
    fprintf(fp, "Total escaped: %lld (%.3f %%)\n", totEscaped, percentEscaped);
    fprintf(fp, "Emitted - Detected: %lld\n", diffEmDet);
  fprintf(fp, "Check (Em-Det vs Esc): %lld vs %lld\n", diffEmDet, totEscaped);
  fprintf(fp, "Gate t <= %.3f ns\n", tGate_ns);
  fprintf(fp, "Events (all detections <= gate): %lld (eff = %.6f %%)\n", nEventsAllDetBeforeGate, effAllDetBeforeGate);
  if (nEventsAnyBeforeGate >= 0) fprintf(fp, "Events (any detection <= gate): %lld (eff = %.6f %%)\n", nEventsAnyBeforeGate, effAnyBeforeGate);
    const double meanMult = (totEvents > 0) ? (sumDetectedD / (double)totEvents) : 0.0;
    const double meanLast = (cntLastDetect > 0) ? (sumLastDetect / (double)cntLastDetect) : -1.0;
    fprintf(fp, "Mean measured multiplicity: %.6f\n", meanMult);
    fprintf(fp, "Mean average detection time [ns]: %.6f\n", meanAvgDetTime);
    fprintf(fp, "Mean last detection time [ns]: %.6f\n", meanLast);
    fclose(fp);
  }

  std::cout << "Saved plots to: " << outBase << "_multiplicities.png and _times.png\n";

  // ================= Save all histograms to a ROOT file =================
  TString rootOut = TString::Format("%s_hists.root", outBase.Data());
  std::unique_ptr<TFile> fout(TFile::Open(rootOut, "RECREATE"));
  if (fout && !fout->IsZombie()) {
    // Write main histograms
    hMult->Write(); hR1->Write(); hR2->Write(); hR3->Write(); hR4->Write();
    hEce->Write(); hLast->Write(); if (hAvgDetTime) hAvgDetTime->Write();

    // Summary tree with means and totals
    Long64_t out_totEvents = totEvents;
    Long64_t out_totEmitted = totEmitted;
    Long64_t out_totDetected = totDetected;
    Long64_t out_totEscaped = totEscaped;
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
    tSummary->Branch("meanMeasuredMultiplicity", &out_meanMult, "meanMeasuredMultiplicity/D");
    tSummary->Branch("meanAvgDetectionTime_ns", &out_meanAvgDetTime, "meanAvgDetectionTime_ns/D");
  tSummary->Branch("meanLastDetectionTime_ns", &out_meanLastDetTime, "meanLastDetectionTime_ns/D");
  // Gate outputs
  Double_t out_tGate_ns = tGate_ns;
  Long64_t out_nEventsAllDetBeforeGate = nEventsAllDetBeforeGate;
  Double_t out_effAllDetBeforeGate = effAllDetBeforeGate;
  Long64_t out_nEventsAnyBeforeGate = (nEventsAnyBeforeGate>=0 ? nEventsAnyBeforeGate : -1);
  Double_t out_effAnyBeforeGate = (effAnyBeforeGate>=0 ? effAnyBeforeGate : -1.0);
  tSummary->Branch("tGate_ns", &out_tGate_ns, "tGate_ns/D");
  tSummary->Branch("eventsAllDetectionsBeforeGate", &out_nEventsAllDetBeforeGate, "eventsAllDetectionsBeforeGate/L");
  tSummary->Branch("effAllDetectionsBeforeGate_percent", &out_effAllDetBeforeGate, "effAllDetectionsBeforeGate_percent/D");
  tSummary->Branch("eventsAnyDetectionBeforeGate", &out_nEventsAnyBeforeGate, "eventsAnyDetectionBeforeGate/L");
  tSummary->Branch("effAnyDetectionBeforeGate_percent", &out_effAnyBeforeGate, "effAnyDetectionBeforeGate_percent/D");
    tSummary->Fill();
    tSummary->Write();
  }

  // ================= PARIS per-copy spectra =================
  if (tParis) {
    Int_t    pe_eventID = 0;
    Int_t    pe_copy = -1;
    Double_t pe_eCe = 0.0;
    Double_t pe_eNaI = 0.0;

    tParis->SetBranchAddress("eventID", &pe_eventID);
    tParis->SetBranchAddress("copy",    &pe_copy);
    tParis->SetBranchAddress("eCe_keV", &pe_eCe);
    tParis->SetBranchAddress("eNaI_keV", &pe_eNaI);

    // Discover copies and ranges
    std::vector<int> copies;
    copies.reserve(16);
    double maxCe = 0.0, maxNaI = 0.0;
    const Long64_t nParis = tParis->GetEntries();
    for (Long64_t i = 0; i < nParis; ++i) {
      tParis->GetEntry(i);
      if (std::find(copies.begin(), copies.end(), pe_copy) == copies.end()) copies.push_back(pe_copy);
      if (pe_eCe  > maxCe)  maxCe  = pe_eCe;
      if (pe_eNaI > maxNaI) maxNaI = pe_eNaI;
    }
    std::sort(copies.begin(), copies.end());
    if (maxCe  <= 0) maxCe  = 5000.0;
    if (maxNaI <= 0) maxNaI = 5000.0;

    int ceBins  = std::clamp((int)std::ceil(maxCe/2.0), 100, 5000);
    int naiBins = std::clamp((int)std::ceil(maxNaI/2.0), 100, 5000);

    // Create histograms per copy
    std::map<int, TH1D*> hCeByCopy;
    std::map<int, TH1D*> hNaIByCopy;
    for (int c : copies) {
      hCeByCopy[c]  = new TH1D(Form("hCe_copy%d", c),  Form("PARIS copy %d;EdepCe [keV];Counts", c),  ceBins,  0, maxCe);
      hNaIByCopy[c] = new TH1D(Form("hNaI_copy%d", c), Form("PARIS copy %d;EdepNaI [keV];Counts", c), naiBins, 0, maxNaI);
    }

    // Fill
    for (Long64_t i = 0; i < nParis; ++i) {
      tParis->GetEntry(i);
      if (pe_copy < 0) continue;
      if (pe_eCe  > 0) hCeByCopy[pe_copy]->Fill(pe_eCe);
      if (pe_eNaI > 0) hNaIByCopy[pe_copy]->Fill(pe_eNaI);
    }

  // Draw Ce per-copy in a grid
    int N = (int)copies.size();
    int cols = (N >= 9 ? 3 : (N >= 6 ? 3 : std::min(3, N > 0 ? N : 1)));
    int rows = (N + cols - 1) / cols;
    auto c3 = new TCanvas("c3", "PARIS Ce per copy", 1200, 400*rows);
    c3->Divide(cols, rows);
    int pad = 1;
    int colorBase = 600; // kBlue
    for (int c : copies) {
      c3->cd(pad++);
      hCeByCopy[c]->SetLineColor((colorBase + c) % 9 + 1); // simple color variation
      hCeByCopy[c]->Draw("HIST");
    }
    c3->SaveAs(TString::Format("%s_paris_ce.png", outBase.Data()));
    // Save to ROOT file (if open)
    if (fout && !fout->IsZombie()) {
      auto* dPARIS = fout->mkdir("PARIS");
      dPARIS->cd();
      auto* dCe = dPARIS->mkdir("Ce"); dCe->cd();
      for (auto& kv : hCeByCopy) kv.second->Write();
      // Go back up for NaI
      dPARIS->cd();
      bool anyNaI = false; for (auto& kv : hNaIByCopy) if (kv.second->GetEntries() > 0) { anyNaI = true; break; }
      if (anyNaI) {
        auto* dNaI = dPARIS->mkdir("NaI"); dNaI->cd();
        for (auto& kv : hNaIByCopy) kv.second->Write();
      }
      fout->cd();
    }

    // Draw NaI per-copy if there is content (optional)
    bool anyNaI = false;
    for (auto& kv : hNaIByCopy) if (kv.second->GetEntries() > 0) { anyNaI = true; break; }
    if (anyNaI) {
      auto c4 = new TCanvas("c4", "PARIS NaI per copy", 1200, 400*rows);
      c4->Divide(cols, rows);
      pad = 1;
      for (int c : copies) {
        c4->cd(pad++);
        hNaIByCopy[c]->SetLineColor((colorBase + c) % 9 + 1);
        hNaIByCopy[c]->Draw("HIST");
      }
      c4->SaveAs(TString::Format("%s_paris_nai.png", outBase.Data()));
    }
  }

  if (fout && !fout->IsZombie()) { fout->Write(); fout->Close(); }
}
