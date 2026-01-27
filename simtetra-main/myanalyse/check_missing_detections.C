#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <map>
#include <set>

void check_missing_detections(const char* infile="output.root", int sample=20) {
  TFile f(infile);
  if (f.IsZombie()) { std::cerr<<"Cannot open "<<infile<<"\n"; return; }
  TTree* tp = (TTree*)f.Get("NeutronPrimaries");
  if (!tp) { std::cerr<<"No NeutronPrimaries tree\n"; return; }

  int eventID; int trackID; double detectTime;
  tp->SetBranchAddress("eventID", &eventID);
  tp->SetBranchAddress("trackID", &trackID);
  tp->SetBranchAddress("detectTime_ns", &detectTime);

  std::map<int, std::set<int>> emittedByEvent;
  std::map<int, std::set<int>> detectedByEvent;

  Long64_t N = tp->GetEntries();
  for (Long64_t i=0;i<N;++i) {
    tp->GetEntry(i);
    emittedByEvent[eventID].insert(trackID);
    if (detectTime >= 0.0) detectedByEvent[eventID].insert(trackID);
  }

  int nEvents = 0; int nProblems = 0;
  for (const auto &kv : emittedByEvent) {
    ++nEvents;
    int eid = kv.first;
    int nEm = kv.second.size();
    int nDet = detectedByEvent[eid].size();
    if (nEm != nDet) {
      ++nProblems;
      if (sample>0) {
        std::cout<<"Event "<<eid<<": emitted="<<nEm<<" detected="<<nDet<<"\n";
        std::cout<<"  emitted trackIDs: ";
        for (int t: kv.second) std::cout<<t<<" ";
        std::cout<<"\n  detected trackIDs: ";
        for (int t: detectedByEvent[eid]) std::cout<<t<<" ";
        // show missing
        std::cout<<"\n  missing: ";
        for (int t: kv.second) if (detectedByEvent[eid].count(t)==0) std::cout<<t<<" ";
        std::cout<<"\n---\n";
        --sample;
      }
    }
  }
  std::cout<<"Scanned events: "<<nEvents<<" problematic: "<<nProblems<<"\n";
}
