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
#include "config/list_xsecs.icc"

#define SUM(vec) std::accumulate((vec).begin(), (vec).end(), 0);
#define SUM_GT(vec,num) std::accumulate((vec).begin(), (vec).end(), 0, [](float x,float y) { return ((y > (num)) ? x+y : x); });
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

// turn on to apply json file to data
const bool applyGoodRunList = true;

int ScanChain(TChain *ch, string sample, string outdir, int nEventsSample = -1) {

  TFile* fout = new TFile(Form("%s/%s.root", outdir.c_str(), sample.c_str()), "RECREATE");
  // H1(met, 50, 0, -1);

  unsigned int nEventsChain = ch->GetEntries();
  unsigned int nEventsTotal = 0;
  unsigned int nPassedTotal = 0;
  unsigned int nDuplicates = 0;
  TFile *currentFile = 0;
  TObjArray *listOfFiles = ch->GetListOfFiles();
  TIter fileIter(listOfFiles);
  tqdm bar;
  const bool show_progress = (sample.find("_") == string::npos);

  map<string,TH1*> hvec;
  // TFile* outfile = new TFile("hists.root", "RECREATE");

  // set configuration parameters
  gconf.GetConfigs(2017);
  gconf.GetSampleType("/"+sample);

  if (applyGoodRunList && gconf.is_data) {
    const char* json_file = "../NanoCORE/data/goodrun/Cert_271036-325175_13TeV_Combined161718_JSON_snt.txt";
    cout << ">>> Loading goodrun json file: " << json_file << endl;
    set_goodrun_file(json_file);
  }

  float scaleToLumi = 1;
  if (!gconf.is_data) {
    if (nEventsSample <= 0) {
      cerr << ">>> WARNING! The input total number of sample " << sample << " is invalid!! Setting it 1000." << endl;
      nEventsSample = 1000;
    }
    scaleToLumi = gconf.lumi * 1000 / nEventsSample;
    if (auto itxs = default_xsec_list.find(sample); itxs != default_xsec_list.end()) {
      scaleToLumi *= itxs->second;
    } else {
      bool found = false;
      for (auto& xs : default_xsec_list) {
        if (sample.find(xs.first) == 0) {
          scaleToLumi *= xs.second;
          found = true;
          break;
        }
      }
      if (!found) {
        cerr << ">>> WARNING! Cannot find the xsec in the default list for sample " << sample << "!! Exiting!" << endl;
        exit(4);
      }
    }
  }

  muoncorr = new RoccoR();
  randomGenerator = new TRandom3();

  while ( (currentFile = (TFile*)fileIter.Next()) ) {
    TFile *file = TFile::Open( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("Events");
    TString filename(currentFile->GetTitle());

    tree->SetCacheSize(128*1024*1024);
    tree->SetCacheLearnEntries(100);

    gconf.GetConfigsFromDatasetName(filename.Data());

    muoncorr->reset();
    muoncorr->init(Form("../NanoCORE/data/RocesterCorrection/RoccoR%d.txt", gconf.year));

    // auto psRead = new TTreePerfStats("readPerf", tree);
    nt.Init(tree);

    for( unsigned int evt = 0; evt < tree->GetEntriesFast(); ++evt) {

      nt.GetEntry(evt);
      tree->LoadTree(evt);

      nEventsTotal++;
      if (show_progress) bar.progress(nEventsTotal, nEventsChain);

      // if (event > 50000) break;

      if (gconf.is_data) {
        if (applyGoodRunList && !goodrun(run(), luminosityBlock())) continue;
        duplicate_removal::DorkyEventIdentifier id(run(), event(), luminosityBlock());
        if (is_duplicate(id)) {
          ++nDuplicates;
          continue;
        }
      }

      bool passMETfilt = passesMETfilters(gconf.is_data);
      if (!passMETfilt) continue;

      nPassedTotal++;

      auto [tightElectrons, looseElectrons] = getElectrons();
      auto [tightMuons, looseMuons] = getMuons();
      auto photons = getPhotons();

      // bool isPhotonDatadriven = false;  //
      // bool isEE = (tightElectrons.size() >= 2 && !isPhotonDatadriven); //2 good electrons
      // bool isMuMu = (tightMuons.size() >= 2 && !isPhotonDatadriven); //2 good muons
      // bool isGamma = (photons.size() == 1 && isPhotonDatadriven); //1 good photon

      bool isEE = (tightElectrons.size() == 2); //2 good electrons
      bool isMuMu = (tightMuons.size() == 2); //2 good muons
      bool isGamma = (photons.size() == 1); //1 good photon

      if (!isEE && !isMuMu && !isGamma) //not a good lepton pair or photon (if datadriven)
        continue;

      // if (isGamma) continue;  // don't study photon for now

      if (isGamma && !passTriggerSelections(4)) continue;
      else if (isEE && !passTriggerSelections(2)) continue;
      else if (isMuMu && !passTriggerSelections(1)) continue;

      bool passLeptonVeto = true;
      if (isGamma)
        passLeptonVeto = (looseElectrons.empty() and looseMuons.empty());
      else if (isEE)
        passLeptonVeto = (looseElectrons.size() == 2 and looseMuons.empty());
      else if (isMuMu)
        passLeptonVeto = (looseElectrons.empty() and looseMuons.size() == 2);

      if (!passLeptonVeto) continue;
      // Add track+tau veto also?? <-- to be studied


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
      if (isGamma) {
        boson = photons[0].p4;
        giveMassToPhoton(boson);
      } else {
        boson = tightLeptons[0].p4 + tightLeptons[1].p4;
      }

      // vector<Particle*> isoobjs;
      auto jets = getJets(looseMuons, looseElectrons, photons);

      bool passBTagVeto = true;
      for (auto const &jet : jets) {
        if (jet.bTag > gconf.WP_DeepCSV_loose && fabs(jet.p4.Eta()) < 2.5) {
          passBTagVeto = false;
          break;
        }
      }

      string jetcat = "geq1j";
      if (jets.size() == 0)
        jetcat = "eq0j";
      else if (PassVBFcuts(jets, boson))
        jetcat = "vbf";

      string lepcat = "gamma";
      if (isEE) lepcat = "ee";
      else if (isMuMu) lepcat = "mumu";

      TLorentzVector met_p4;
      met_p4.SetPtEtaPhiM(MET_pt(), 0., MET_phi(), 0.);

      // float met_sig = MET_significance();

      // if (MET_pt() < 125) continue;
      if (fabs(boson.M() - mZ) > 15) continue;
      if (boson.Pt() < 55) continue;

      float deltaPhi_MET_Boson = deltaPhi(boson.Phi(), MET_phi());
      if (deltaPhi_MET_Boson < 0.5) continue;


      if (!passBTagVeto) continue;

      bool passDeltaPhiJetMET = true;

      for (auto const &jet : jets) {
        if (deltaPhi(jet.p4.Phi(), MET_phi()) < 0.5) {
          passDeltaPhiJetMET = false;
          break;
        }
      }
      if (!passDeltaPhiJetMET) continue;

      if (MET_pt() < 125) continue;

      // Events passes all selections

      float mtll = getDileptonMT(boson, met_p4);

      float weight = scaleToLumi;

      auto fillhists = [&](string s) {
        plot1d("h_njets"+s,  jets.size(), weight, hvec, ";N(jets)"  , 6,  0, 6);
        plot1d("h_met"+s,    MET_pt()  , weight, hvec, ";E_{T}^{miss} [GeV]"  , 30,  0, 750);
        plot1d("h_metphi"+s, MET_phi() , weight, hvec, ";#phi(E_{T}^{miss})"  , 34, -3.4, 3.4);
        // plot1d("h_metsig"+s, met_sig   , weight, hvec, ";#sigma(E_{T}^{miss}) " , 20,  0, 40);
        // Z quantities
        if (isEE || isMuMu) {
          plot1d("h_mll"+s,    boson.M() , weight, hvec, ";M_{ll} [GeV]"        , 50,  0, 500);
          plot1d("h_ptll"+s,  boson.Pt() , weight, hvec, ";p_{T}(ll) [GeV]" , 40,  0, 800);
          plot1d("h_etall"+s, boson.Eta(), weight, hvec, ";#eta(ll) [GeV]"     , 40,  -5, 5);

          plot1d("h_lep1pt"+s,   tightLeptons[0].p4.Pt() , weight, hvec, ";p_{T}(lep1) [GeV]"  , 25,  0, 500);
          plot1d("h_lep2pt"+s,   tightLeptons[1].p4.Pt() , weight, hvec, ";p_{T}(lep2) [GeV]"  , 20,  0, 400);
          plot1d("h_lep1eta"+s,  tightLeptons[0].p4.Eta() , weight, hvec, ";#eta(lep1)"        , 36, -2.4, 2.4);
          plot1d("h_lep2eta"+s,  tightLeptons[1].p4.Eta() , weight, hvec, ";#eta(lep2)"        , 36, -2.4, 2.4);
        } else if (isGamma) {
          plot1d("h_Zmass"+s, boson.M()  , weight, hvec, ";fake M(#gamma)[GeV]"       , 20,  0, 500);
          plot1d("h_ptgamma"+s, boson.Pt()  , weight, hvec, ";p_{T}(#gamma) [GeV]"    , 40,  0, 800);
          plot1d("h_etagamma"+s, boson.Eta(), weight, hvec, ";#eta(#gamma) [GeV]"     , 40,  -5, 5);
        }

        plot1d("h_mtll"+s, mtll  , weight, hvec, ";M_{T}(ll) [GeV]"        , 60,  0, 1500);

        const vector<float> mtbin1 = {0, 75, 150, 225, 300, 375, 450, 525, 600, 725, 850, 975, 1100, 1350, 1600, 2100, 3000};
        plot1d("h_mtll_b1"+s,   mtll , weight, hvec, ";M_{T}(ll) [GeV]" , mtbin1.size()-1, mtbin1.data());
        const vector<float> mtbin3 = {0, 75, 150, 225, 300, 375, 450, 525, 600, 725, 850, 975, 1100, 1350, 1600, 2100, 3000};
        plot1d("h_mtll_b3"+s,   mtll , weight, hvec, ";M_{T}(ll) [GeV]" , mtbin3.size()-1, mtbin3.data());
      };

      fillhists("");
      fillhists("_"+jetcat);
      fillhists("_"+lepcat);
      fillhists("_"+jetcat+"_"+lepcat);

      // int njets, nbtags;
      // float ht;
      // std::tie(njets,nbtags,ht) = getJetInfo(leps);

      // float met = MET_pt();

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

    } // Event loop

    delete file;


  } // File loop
  if (show_progress) bar.finish();

  fout->cd();
  for (auto& h : hvec) {
    if (h.first.find("phi") != string::npos) {
      int lastbin = h.second->GetNbinsX(); // to move the overflow bin to the last bin
      h.second->SetBinContent(lastbin, h.second->Integral(lastbin, -1));
    }
    string dirname = "hzz2l2nu";
    for (string dirsuf : {"_eq0j_ee", "_eq0j_mumu", "_geq1j_ee", "_geq1j_mumu", "_vbf_ee", "_vbf_mumu",
            "_ee", "_mumu", "_eq0j", "_geq1j",  "_vbf" }) {
      if (h.first.find(dirsuf) != string::npos) {
        dirname += dirsuf;
        break;
      }
    }
    TDirectory* dir = (TDirectory*) fout->Get(dirname.c_str());
    if (dir == 0) dir = fout->mkdir(dirname.c_str());
    dir->cd();

    h.second->Write();
  }

  // fout->Write();
  fout->Close();

  delete fout; fout = nullptr;
  delete muoncorr; muoncorr = nullptr;

  cout << "\n---------------------------------" << endl;
  cout << nEventsTotal << " Events Processed, where " << nDuplicates << " duplicates were skipped, and ";
  cout << nPassedTotal << " Events passed all selections." << endl;
  cout << "---------------------------------\n" << endl;

  return 0;
}
