// -*- C++ -*-
#pragma GCC diagnostic ignored "-Wsign-compare"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"
#include "TTreeCache.h"
#include "TTreeCacheUnzip.h"
#include "TTreePerfStats.h"
#include "Math/VectorUtil.h"

#include <iostream>
#include <iomanip>
#include <map>

#include "../NanoCORE/Nano.h"
#include "../NanoCORE/Config.h"
#include "../NanoCORE/MetSelections.h"
#include "../NanoCORE/goodrun.h"
#include "../NanoCORE/dorky.h"
#include "../NanoCORE/tqdm.h"

#include "PhysicsObjects.h"
#include "HZZSelections.h"
#include "Utilities.h"

#define SUM(vec) std::accumulate((vec).begin(), (vec).end(), 0);
#define SUM_GT(vec,num) std::accumulate((vec).begin(), (vec).end(), 0, [](float x,float y){return ((y > (num)) ? x+y : x); });
#define COUNT_GT(vec,num) std::count_if((vec).begin(), (vec).end(), [](float x) { return x > (num); });
#define COUNT_LT(vec,num) std::count_if((vec).begin(), (vec).end(), [](float x) { return x < (num); });


#define H1(name,nbins,low,high) TH1F *h_##name = new TH1F(#name,#name,nbins,low,high);

// #define DEBUG

struct debugger { template<typename T> debugger& operator , (const T& v) { cerr<<v<<" "; return *this; } } dbg;
#ifdef DEBUG
  #define debug(args...) do {cerr << #args << ": "; dbg,args; cerr << endl;} while(0)
#else
  #define debug(args...)
#endif

using namespace std;
using namespace tas;
using namespace utils;

int ScanChain(TChain *ch, string sample, string outdir) {

  TFile* fout = new TFile(Form("%s/%s.root", outdir.c_str(), sample.c_str()), "RECREATE");
  // H1(met, 50, 0, -1);

  int nEventsTotal = 0;
  int nEventsChain = ch->GetEntries();
  TFile *currentFile = 0;
  TObjArray *listOfFiles = ch->GetListOfFiles();
  TIter fileIter(listOfFiles);
  tqdm bar;

  map<string,TH1*> hvec;
  // TFile* outfile = new TFile("hists.root", "RECREATE");

  // set configuration parameters
  gconf.GetConfigs(2017);
  gconf.GetSampleType(sample);

  while ( (currentFile = (TFile*)fileIter.Next()) ) {
    TFile *file = TFile::Open( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("Events");
    TString filename(currentFile->GetTitle());

    tree->SetCacheSize(128*1024*1024);
    tree->SetCacheLearnEntries(100);

    // auto psRead = new TTreePerfStats("readPerf", tree);
    nt.Init(tree);

    for( unsigned int event = 0; event < tree->GetEntriesFast(); ++event) {

      nt.GetEntry(event);
      tree->LoadTree(event);

      nEventsTotal++;
      bar.progress(nEventsTotal, nEventsChain);

      // if (event > 50000) break;

      auto tightElectrons = getElectrons();
      auto tightMuons = getMuons();
      auto photons = getPhotons();
      auto jets = getJets();

      bool isPhotonDatadriven = false;  // 
      bool isEE = (tightElectrons.size() >= 2 && !isPhotonDatadriven); //2 good electrons
      bool isMuMu = (tightMuons.size() >= 2 && !isPhotonDatadriven); //2 good muons
      bool isGamma = (photons.size() == 1 && isPhotonDatadriven); //1 good photon

      if(!isEE && !isMuMu && !isGamma) //not a good lepton pair or photon (if datadriven)
        continue;

      std::vector<Lepton> tightLeptons;

      if (isMuMu) {
        for (auto const &mu : tightMuons) {
          tightLeptons.emplace_back(mu);
        }
      } else if (isEE) {
        for (auto const &e : tightElectrons) {
          tightLeptons.emplace_back(e);
        }
      }

      TLorentzVector boson;
      if (isPhotonDatadriven) {
        boson = photons[0].p4;
        giveMassToPhoton(boson);
      } else {
        boson = tightLeptons[0].p4 + tightLeptons[1].p4;
      }

      int jetcat = geq1j;
      if (jets.size() == 0)
        jetcat = eq0j;
      else if (PassVBFcuts(jets, boson))
        jetcat = vbf;

      TLorentzVector met_p4;
      met_p4.SetPtEtaPhiM(MET_pt(), 0., MET_phi(), 0.);

      float met_sig = MET_significance();

      if (MET_pt() < 80) continue;
      if (fabs(boson.M() - MZ) > 15) continue;
      if (boson.Pt() < 55) continue;

      float deltaPhi_MET_Boson = deltaPhi(boson.Phi(), MET_phi());
      if (deltaPhi_MET_Boson < 0.5) continue;
      
      auto looseElectrons = getElectrons(ID_level::idLoose);
      auto looseMuons = getMuons(ID_level::idLoose);

      bool passLeptonVeto = true;

      if (isPhotonDatadriven)
        passLeptonVeto = (looseElectrons.empty() and looseMuons.empty());
      else if (isEE)
        passLeptonVeto = (looseElectrons.size() == 2 and looseMuons.empty());
      else if (isMuMu)
        passLeptonVeto = (looseElectrons.empty() and looseMuons.size() == 2);

      if (!passLeptonVeto) continue;

      bool passBTagVeto = true;
      for (auto const &jet : jets) {
        if (jet.bTag > gconf.WP_DeepCSV_medium) {
          passBTagVeto = false;
          break;
        }
      }

      if (!passBTagVeto) continue;

      bool passDeltaPhiJetMET = true;

      for (auto const &jet : jets) {
        if (deltaPhi(jet.p4.Phi(), MET_phi()) < 0.5) {
          passDeltaPhiJetMET = false;
          break;
        }
      }
      if (!passDeltaPhiJetMET) continue;

      if (MET_pt() < 80) continue;


      plot1d("h_met",    MET_pt()  , 1, hvec, ";E_{T}^{miss} [GeV]"  , 30,  0, 750);
      plot1d("h_metphi", MET_phi() , 1, hvec, ";#phi(E_{T}^{miss})"  , 34, -3.4, 3.4);
      plot1d("h_metsig", met_sig   , 1, hvec, ";#sigma(E_{T}^{miss}) " , 20,  0, 40);
      plot1d("h_mll",    boson.M() , 1, hvec, ";M_{ll} [GeV]"        , 20,  0, 500);

      // int njets, nbtags;
      // float ht;
      // std::tie(njets,nbtags,ht) = getJetInfo(leps);

      // float met = MET_pt();
      // bool passfilt = passesMETfilters(false);

      // debug(passfilt,nbtags,met,njets,nleps);

      // if (lep1.is_el() && !isTriggerSafeIso_v1(lep1.idx())) continue;
      // if (lep2.is_el() && !isTriggerSafeIso_v1(lep2.idx())) continue;

      // // if (hyp_class != 3 && hyp_class != 4) continue;

      // float weight = genWeight();

      // h_mll->Fill(mll);
      // h_hyp_class->Fill(hyp_class);
      // h_filt->Fill(passfilt);
      // h_nbtags->Fill(nbtags);
      // h_met->Fill(met);
      // h_njets->Fill(njets);
      // h_ht->Fill(ht);
      // h_nleps->Fill(leps.size());
      // h_weight->Fill(weight);

      // h_ptratio->Fill(getPtRatio(lep1.id(),lep1.idx()));
      // h_ptratio->Fill(getPtRatio(lep2.id(),lep2.idx()));
      // h_ptrel->Fill(getPtRel(lep1.id(),lep1.idx()));
      // h_ptrel->Fill(getPtRel(lep2.id(),lep2.idx()));


    } // Event loop

    delete file;


  } // File loop
  bar.finish();

  fout->cd();
  for (auto& h : hvec) {
    if (h.first.find("phi") != string::npos) {
      int lastbin = h.second->GetNbinsX(); // to move the overflow bin to the last bin
      h.second->SetBinContent(lastbin, h.second->Integral(lastbin, -1));
    }
    h.second->Write();
  }

  // fout->Write();
  fout->Close();
  return 0;
}
