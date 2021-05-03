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
#include <bitset>
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
const bool doDuplicateRemoval = false;
const bool applySyncCuts = true;
const bool ZGestCheck = false;

int ScanChain(TChain *ch, string sample, string outdir, int nevtsamp = -1) {

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
    const char* json_file = "config/goodrun/Cert_271036-325175_13TeV_Combined161718_JSON_snt.txt";
    cout << ">>> Loading goodrun json file: " << json_file << endl;
    set_goodrun_file(json_file);
  }

  // double sum_weights(0.);
  double sum_effevts(0.);
  double sum_genwgts(0.);

  muoncorr = new RoccoR();
  randomGenerator = new TRandom3();

  bool is_GJets_HT = (sample.find("GJets_HT") == 0);

  while ( (currentFile = (TFile*)fileIter.Next()) ) {
    TFile *file = TFile::Open( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("Events");
    TString filename(currentFile->GetTitle());

    tree->SetCacheSize(128*1024*1024);
    tree->SetCacheLearnEntries(100);

    gconf.GetConfigsFromDatasetName(filename.Data());

    double scaleToLumi = 1;
    if (!gconf.is_data) {
      if (nevtsamp <= 0) {
        cerr << ">>> WARNING! The input total number of sample " << sample << " is invalid!! Setting it 1000." << endl;
        nevtsamp = 1000;
      }
      scaleToLumi = gconf.lumi * 1000 / nevtsamp;
      if (auto itxs = default_xsec_list.find(sample); itxs != default_xsec_list.end()) {
        scaleToLumi *= itxs->second;
      } else {
        bool found = false;
        for (auto& xs : default_xsec_list) {
          if (sample.find(xs.first) == 0) {
            scaleToLumi *= xs.second;
            found = true;
            cout << ">>> Loading xsec from default list for sample " << sample << " and we get " << xs.second << "pb." << endl;
            break;
          }
        }
        if (!found) {
          cerr << ">>> WARNING! Cannot find the xsec in the default list for sample " << sample << "!! Exiting!" << endl;
          cout << "             gconf.lumi= " << gconf.lumi << ", nevtsamp= " << nevtsamp << ", scaleToLumi= " << scaleToLumi << endl;
          exit(4);
        }
      }
    }

    muoncorr->reset();
    muoncorr->init(Form("config/RocesterCorrection/RoccoR%d.txt", gconf.year));

    // auto psRead = new TTreePerfStats("readPerf", tree);

    nt.Init(tree);

    for( unsigned int evt = 0; evt < tree->GetEntriesFast(); ++evt) {
      // if (evt > 500) break;

      nt.GetEntry(evt);
      tree->LoadTree(evt);

      nEventsTotal++;
      if (!gconf.is_data) {
        sum_effevts += (genWeight() > 0)? 1 : -1;
        sum_genwgts += genWeight();
      }

      if (show_progress) bar.progress(nEventsTotal, nEventsChain);

      if (gconf.is_data) {
        if (applyGoodRunList && !goodrun(run(), luminosityBlock())) continue;
        if (doDuplicateRemoval) {
          duplicate_removal::DorkyEventIdentifier id(run(), event(), luminosityBlock());
          if (is_duplicate(id)) {
            ++nDuplicates;
            continue;
          }
        }
      }

      float weight = (gconf.is_data)? 1 : (genWeight() > 0)? scaleToLumi : -1.0*scaleToLumi;

      if ((event() == 165136538 || event() == 165588499 || event() == 569563314 || event() == 175850032 || event() == 922337055 || event() == 1257305638 || event() == 532551390 || event() == 628303537 || event() == 2803911702 || event() == 3077899767 || event() == 1055641106 || event() == 932848436 || event() == 120762842 || event() == 2485139493 || event() == 2303213037))
        cout << __LINE__ << ": run()= " << run() << ", luminosityBlock()= " << luminosityBlock() << ", event()= " << event() << endl;

      // The following skim cuts are in the ntuple skiming steps
      // 2l2nu cut: nTightElection + nTightMuon >= 2, Z pt > 50, mZ > 50

      plot1d("h_preselec_steps", 0, 1, hvec, ";M_{ll} [GeV]" , 20,  0, 20);

      bool passMETfilt = passesMETfilters(gconf.is_data);
      if (!passMETfilt) continue;

      plot1d("h_preselec_steps", 1, 1, hvec, ";M_{ll} [GeV]" , 20,  0, 20);

      if ((event() == 165136538 || event() == 165588499 || event() == 569563314 || event() == 175850032 || event() == 922337055 || event() == 1257305638 || event() == 532551390 || event() == 628303537 || event() == 2803911702 || event() == 3077899767 || event() == 1055641106 || event() == 932848436 || event() == 120762842 || event() == 2485139493 || event() == 2303213037))
        cout << __LINE__ << ": run()= " << run() << ", luminosityBlock()= " << luminosityBlock() << ", event()= " << event() << endl;

      // To apply the same skim cut as the NanoAOD
      if (applySyncCuts) {
        // if (isGamma)
        //   if (!InstrMETPrunerCuts()) continue;
        // else
        // if (!ZZ2l2vPrunerCuts()) continue;
        if (!InstrMETPrunerCuts()) continue;
      }

      if ((event() == 165136538 || event() == 165588499 || event() == 569563314 || event() == 175850032 || event() == 922337055 || event() == 1257305638 || event() == 532551390 || event() == 628303537 || event() == 2803911702 || event() == 3077899767 || event() == 1055641106 || event() == 932848436 || event() == 120762842 || event() == 2485139493 || event() == 2303213037))
        cout << __LINE__ << ": run()= " << run() << ", luminosityBlock()= " << luminosityBlock() << ", event()= " << event() << endl;

      plot1d("h_preselec_steps", 2, 1, hvec, ";M_{ll} [GeV]" , 20,  0, 20);

      nPassedTotal++;

      int nFailCuts = 0;

      auto [tightElectrons, looseElectrons] = getElectrons();
      auto [tightMuons, looseMuons] = getMuons();
      auto photons = getPhotons(looseMuons, looseElectrons);
      auto jets = getJets(looseMuons, looseElectrons, photons);  // TODO: temporary move to later

      // if ((sample == "ZGJets_nunu_nlo_incl" || sample == "ZGJets_ll_nlo_incl") && pt_photon() > 135) continue;
      // if ((sample == "ZGJets_nunu_nlo_pTG_130" || sample == "ZGJets_ll_nlo_pTG_130") && pt_photon() < 135) continue;

      if ((event() == 165136538 || event() == 165588499 || event() == 569563314 || event() == 175850032 || event() == 922337055 || event() == 1257305638 || event() == 532551390 || event() == 628303537 || event() == 2803911702 || event() == 3077899767 || event() == 1055641106 || event() == 932848436 || event() == 120762842 || event() == 2485139493 || event() == 2303213037))
        cout << __LINE__ << ": run()= " << run() << ", luminosityBlock()= " << luminosityBlock() << ", event()= " << event() << endl;

      auto fillZGhists = [&](string s = "ZGcheck") {
        if (!ZGestCheck) return;
        if (gconf.is_data) return;
        if (nFailCuts > 0) return;

        plot1d("h_nphotons_"+s, photons.size(), weight, hvec, ";N(photon)"  , 6,  0, 6);
        if (photons.size() < 1) return;

        plot1d("h_njets_"+s, jets.size(), weight, hvec, ";N(jets)"  , 6,  0, 6);
        plot1d("h_nleptons_"+s, (looseMuons.size()+looseElectrons.size()), weight, hvec, ";N(leps)"  , 6,  0, 6);
        plot1d("h_nmuons_"+s, looseMuons.size(), weight, hvec, ";N(muon)"  , 6,  0, 6);
        plot1d("h_nelecs_"+s, looseElectrons.size(), weight, hvec, ";N(elec)"  , 6,  0, 6);

        int nquarks(0), nleps(0), nnus(0), nphs(0), nqhard(0), nlephard(0), nphhard(0), nZs(0), nleps_fc(0);
        int idxZ(-1), idxPh(-1), idxl1(-1), idxl2(-1);
        int idxsumll(0), idxsumll_fc(0);
        LorentzVector p4ll, p4ll_fc;
        for (size_t i = 0; i < GenPart_pdgId().size(); ++i) {
          int id = abs(GenPart_pdgId()[i]);

          if ((id >= 11 && id <= 16) && (GenPart_statusFlags()[i] & 1<<12) != 0 && (GenPart_statusFlags()[i] & 1<<8) != 0) {
            nleps_fc++;
            p4ll_fc += GenPart_p4()[i];
            idxsumll_fc += i;
          }
          if ((GenPart_statusFlags()[i] & 1<<13) == 0) continue;

          if (id == 11 || id == 13 || id == 15) nleps++;
          if (id == 22) nphs++;
          if (id >= 1 && id <= 5) nquarks++;


          if ((GenPart_statusFlags()[i] & 1<<8) == 0) continue;

          if (id >= 11 && id <= 16) {
            if (id % 2 == 1) nlephard++;
            else nnus++;
            p4ll += GenPart_p4()[i];
            idxsumll += i;
            if (idxl1 < 0) idxl1 = i;
            else idxl2 = i;
          }
          if (id >= 1 && id <= 5) nqhard++;
          if (id == 22) {
            nphhard++;
            idxPh = i;
          }
          if (id == 23) {
            nZs++;
            idxZ = i;
          }
        }
        if (nlephard + nnus < 2) return;
        if (GenPart_pt()[idxl2] > GenPart_pt()[idxl1]) std::swap(idxl1, idxl2);

        plot1d("h_nZs_"+s, nZs, weight, hvec, ";N(gen-Zs)"  , 6,  0, 6);
        plot1d("h_nleps_"+s, nleps, weight, hvec, ";N(gen-leps)"  , 6,  0, 6);
        plot1d("h_nquarks_"+s, nquarks, weight, hvec, ";N(gen-quarks)"  , 6,  0, 6);
        plot1d("h_nnus_"+s, nnus, weight, hvec, ";N(gen-#nus)"  , 6,  0, 6);
        plot1d("h_nphs_"+s, nphs, weight, hvec, ";N(gen-#gamma)"  , 6,  0, 6);

        plot1d("h_nleps_hard_"+s, nleps, weight, hvec, ";N(gen-leps hard)"  , 6,  0, 6);
        plot1d("h_nqrks_hard_"+s, nqhard, weight, hvec, ";N(gen-quarks hard)"  , 6,  0, 6);
        plot1d("h_nphs_hard_"+s, nphhard, weight, hvec, ";N(gen-#gamma hard)"  , 6,  0, 6);

        plot1d("h_gen_lep1pt_"+s, GenPart_pt()[idxl1], weight, hvec, ";pT(gen-lep1)"  , 90,  0, 450);
        plot1d("h_gen_lep2pt_"+s, GenPart_pt()[idxl2], weight, hvec, ";pT(gen-lep2)"  , 90,  0, 450);
        
        plot1d("h_gen_lep1eta_"+s, std::max(GenPart_eta()[idxl1], -4.99f), weight, hvec, ";#eta(gen-lep1)"  , 100,  -5.0f, 5.0f);
        plot1d("h_gen_lep2eta_"+s, std::max(GenPart_eta()[idxl2], -4.99f), weight, hvec, ";#eta(gen-lep2)"  , 100,  -5.0f, 5.0f);

        // plot1d("h_genmet_"+s, genmet_met(), weight, hvec, ";pT(gen-MET)"  , 90,  0, 450);

        if (nZs == 0) {
          plot1d("h_gen_mll_nonZ_"+s, p4ll.M(), weight, hvec, ";M(gen-ll)"  , 90,  0, 450);
          plot1d("h_gen_llpt_nonZ_"+s, p4ll.pt(), weight, hvec, ";pT(gen-ll)"  , 90,  0, 450);
        }

        if (nZs < 1 || nphhard < 1) return;
        // if (fabs(p4ll.M() - 91.2) > 10) return;

        float dphiZG = deltaPhi(GenPart_phi()[idxPh], GenPart_phi()[idxZ]);
        float dphillG = deltaPhi(p4ll.phi(), GenPart_phi()[idxPh]);
        plot1d("h_dphiZG_"+s, dphiZG, weight, hvec, ";#Delta#phi(Z, #gamma)"  , 64,  0, 3.2f);
        plot1d("h_dphillG_"+s, dphillG, weight, hvec, ";#Delta#phi(ll, #gamma)"  , 64,  0, 3.2f);
        // if (dphiZG < 1.0) return;

        if (nleps == 2 || nnus == 2) {
          plot1d("h_genZ_mll_"+s, p4ll.M(), weight, hvec, ";M(gen-ll)"  , 90,  0, 180);
          plot1d("h_genZ_llpt_"+s, p4ll.pt(), weight, hvec, ";pT(gen-ll)"  , 90,  0, 450);

          plot1d("h_genZ_fc_mll_"+s, p4ll_fc.M(), weight, hvec, ";M(gen-ll) 1st copy"  , 90,  0, 180);
          plot1d("h_genZ_fc_llpt_"+s, p4ll_fc.pt(), weight, hvec, ";pT(gen-ll) 1st copy"  , 90,  0, 450);
        } 

        if (static int pct284=0; fabs(p4ll.M()-GenPart_mass()[idxZ]) > 20 && s == "ZGcheck" && pct284++<20) {
          plot1d("h_dphillG_dmlarge_"+s, dphillG, weight, hvec, ";#Delta#phi(ll, #gamma)"  , 64,  0, 3.2f);

          cout << __LINE__ << ": evt= " << evt << ", p4ll.M()= " << p4ll.M() << ", p4ll_fc.M()= " << p4ll_fc.M() << ", GenPart_mass()[idxZ]= " << GenPart_mass()[idxZ] << ", idxsumll= " << idxsumll << ", idxsumll_fc= " << idxsumll_fc << ", dphillG= " << dphillG << endl;
          for (size_t i = 0; i < GenPart_pdgId().size(); ++i) {
            // if ((GenPart_statusFlags()[i] & 1<<8) == 0) continue;
            int id = abs(GenPart_pdgId()[i]);
            if ((id < 11 || id > 16) && id != 23 && id != 22 && GenPart_genPartIdxMother()[i] != idxZ) continue;
            cout << __LINE__ << ": evt= " << evt << ", i= " << i << ", GenPart_pdgId[i]= " << GenPart_pdgId()[i] << ", pt[i]= " << GenPart_pt()[i] << ", eta[i]= " << GenPart_eta()[i] << ", phi[i]= " << GenPart_phi()[i] << ", mass[i]= " << GenPart_mass()[i] << ", GenPart_genPartIdxMother[i]= " << GenPart_genPartIdxMother()[i] << ", fromHardProcess= " << hex << (GenPart_statusFlags()[i] & (1<<8 | 1<<12 | 1<<13)) << dec << endl;
          }
          for (size_t i = 0; i < Electron_pt().size(); ++i) {
            cout << __LINE__ << ": evt= " << evt << ", i= " << i << ", Electron_pdgId[i]= " << 11 << ", pt[i]= " << Electron_pt()[i] << ", eta[i]= " << Electron_eta()[i] << ", phi[i]= " << Electron_phi()[i] << ", mass[i]= " << Electron_mass()[i] << ", GenPart_genPartIdx[i]= " << Electron_genPartIdx()[i] << endl;
          }
          for (size_t i = 0; i < Photon_pt().size(); ++i) {
            cout << __LINE__ << ": evt= " << evt << ", i= " << i << ", Photon_pdgId[i]= " << 22 << ", pt[i]= " << Photon_pt()[i] << ", eta[i]= " << Photon_eta()[i] << ", phi[i]= " << Photon_phi()[i] << ", mass[i]= " << Photon_mass()[i] << ", GenPart_genPartIdx[i]= " << Photon_genPartIdx()[i] << endl;
          }
          for (auto ph : photons) {
            cout << __LINE__ << ": evt= " << evt << ", ph.p4.Pt()= " << ph.p4.Pt() << ", ph.p4.Eta()= " << ph.p4.Eta() << ", ph.p4.Phi()= " << ph.p4.Phi() << endl;
          }

        }

        plot1d("h_genph_pt_"+s, GenPart_pt()[idxPh], weight, hvec, ";pT(gen-#gamma)"  , 90,  0, 450);
        plot1d("h_genZ_pt_"+s, GenPart_pt()[idxZ], weight, hvec, ";pT(gen-Z)"  , 90,  0, 450);
        plot1d("h_genZ_mZ_"+s, GenPart_mass()[idxZ], weight, hvec, ";M(gen-Z)"  , 90,  0, 180);

      };
      fillZGhists("ZGcheck");
      
      bool isEE = (tightElectrons.size() == 2); //2 good Photon
      bool isMuMu = (tightMuons.size() == 2); //2 good muons
      bool isGamma = (photons.size() == 1 && !isEE && !isMuMu); //1 good photon
      bool isllG = (isMuMu || isEE) && (photons.size() == 1);

      if ((event() == 165136538 || event() == 165588499 || event() == 569563314 || event() == 175850032 || event() == 922337055 || event() == 1257305638 || event() == 532551390 || event() == 628303537 || event() == 2803911702 || event() == 3077899767 || event() == 1055641106 || event() == 932848436 || event() == 120762842 || event() == 2485139493 || event() == 2303213037))
        cout << __LINE__ << ": run()= " << run() << ", luminosityBlock()= " << luminosityBlock() << ", event()= " << event() << endl;

      string lepcat = "gamma";
      if (isEE) lepcat = "ee";
      else if (isMuMu) lepcat = "mumu";
      else if (isGamma) lepcat = "gamma";
      if (isllG) lepcat = "llg";

      int istep=0;
      plot1d("h_weight_"+to_string(istep), weight, 1, hvec, ";Event weight"  , 200,  0, 20);
      plot1d("h_met_step"+to_string(istep), MET_pt(), weight, hvec, ";E_{T}^{miss} [GeV]"  , 120,  0, 600);
      plot1d("h_met_unwgtd_step"+to_string(istep), MET_pt(), 1, hvec, ";E_{T}^{miss} [GeV]"  , 120,  0, 600);

      auto fill_passedsteps = [&](string s="") {
        if (nFailCuts == 0) {
          if ((event() == 165136538 || event() == 165588499 || event() == 569563314 || event() == 175850032 || event() == 922337055 || event() == 1257305638 || event() == 532551390 || event() == 628303537 || event() == 2803911702 || event() == 3077899767 || event() == 1055641106 || event() == 932848436 || event() == 120762842 || event() == 2485139493 || event() == 2303213037))
            cout << __LINE__ << ": run()= " << run() << ", luminosityBlock()= " << luminosityBlock() << ", event()= " << event() << ", istep= " << istep << ", MET_pt()= " << MET_pt() << ", s= " << s << endl;
          
          plot1d("h_met_step"+to_string(istep), MET_pt(), weight, hvec, ";E_{T}^{miss} [GeV]"  , 120,  0, 600);
          plot1d("h_passed_steps_"+lepcat, istep , weight, hvec, ";M_{ll} [GeV]" , 20,  0, 20);
          plot1d("h_passed_steps", istep++ , weight, hvec, ";M_{ll} [GeV]" , 20,  0, 20);
        }
      };
      fill_passedsteps();

      if (!isEE && !isMuMu && !isGamma) //not a good lepton pair or photon (if datadriven)
        continue;

      bool isOppositeSign = ((isEE && (tightElectrons[0].charge * tightElectrons[1].charge == -1)) ||
                             (isMuMu && (tightMuons[0].charge * tightMuons[1].charge == -1)) || isGamma);
      if (!isOppositeSign) continue;  // reject events that are not opposite sign
      fill_passedsteps();

      if (isGamma) {
        if (gconf.is_data) {
          double prescale = getPhotonTrigPrescale(photons[0].p4.Pt());
          if (prescale <= 0.0) continue;
          weight *= prescale;
        } else {
          // do nothing for now
        }
      }
      else if (isEE && !passTriggerSelections(2))
        continue;
      else if (isMuMu && !passTriggerSelections(1))
        continue;

      bool passLeptonVeto = true;
      if (isGamma)
        passLeptonVeto = (looseElectrons.empty() and looseMuons.empty());
      else if (isEE)
        passLeptonVeto = (looseElectrons.size() == 2 and looseMuons.empty());
      else if (isMuMu)
        passLeptonVeto = (looseElectrons.empty() and looseMuons.size() == 2);

      if (!passLeptonVeto) continue;
      fill_passedsteps();

      // Add track+tau veto also?? <-- to be studied
      vector<int> isotrack_idxs = getIsoTrackIndices(looseMuons, looseElectrons);
      bool passTrackVeto = (isotrack_idxs.size() > 0);
      if (!passTrackVeto) nFailCuts++;
      if (nFailCuts > 1) continue;
      fill_passedsteps();

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
      if (isGamma || isllG) {
        boson = photons[0].p4;
        giveMassToPhoton(boson);
      } else {
        boson = tightLeptons[0].p4 + tightLeptons[1].p4;
      }

      // Apply the LO to NLO scale factor for GJets_HT samples here
      if(is_GJets_HT) weight *= std::max(1., 1.716910-0.001221*boson.Pt());

      // auto jets = getJets(looseMuons, looseElectrons, photons); // should be here, commented out temporarily

      bool passBTagVeto = true;
      for (auto const &jet : jets) {
        if (jet.bTag > gconf.WP_DeepFlav_loose && fabs(jet.p4.Eta()) < 2.5) {
          passBTagVeto = false;
          break;
        }
      }
      if (!passBTagVeto) nFailCuts++;
      if (nFailCuts > 1) continue;
      fill_passedsteps();

      bool is_VBFcat = PassVBFcuts(jets, boson);

      string jetcat = "geq1j";
      if (jets.size() == 0)
        jetcat = "eq0j";
      else if (jets.size() == 1)
        jetcat = "eq1j";
      else if (jets.size() == 2)
        jetcat = "eq2j";

      TLorentzVector met_p4;
      met_p4.SetPtEtaPhiM(MET_pt(), 0., MET_phi(), 0.);
      if (isllG) met_p4 += tightLeptons[0].p4 + tightLeptons[1].p4;

      auto fill_mll_hists = [&](string s) {
        if (!isMuMu && !isEE) return;
        plot1d("h_mll_"+s, boson.M() , scaleToLumi, hvec, ";M_{ll} [GeV]" , 125,  0, 500);
        plot1d(Form("h_mll_%s_%s", s.c_str(), lepcat.c_str()), boson.M() , scaleToLumi, hvec, ";M_{ll} [GeV]" , 125,  0, 500);
        plot1d(Form("h_mll_%s_%s_%s", s.c_str(), jetcat.c_str(), lepcat.c_str()), boson.M() , scaleToLumi, hvec, ";M_{ll} [GeV]" , 125,  0, 500);
      };
      fill_mll_hists("full");

      bool passZmassWindow = (fabs(boson.M() - mZ) < 15);
      if (isllG) passZmassWindow = (fabs((tightLeptons[0].p4 + tightLeptons[1].p4).M() - mZ) < 15);
      if (!passZmassWindow) nFailCuts++;
      if (nFailCuts > 1) continue;
      fill_passedsteps();

      bool passBosonPtCut = (boson.Pt() >= 55.);
      if (!passBosonPtCut) nFailCuts++;
      if (nFailCuts > 1) continue;
      fill_passedsteps();

      // fill_mll_hists("Zpt55");

      float deltaPhi_MET_Boson = deltaPhi(boson.Phi(), met_p4.Phi());
      bool passDeltaPhiZMET = (deltaPhi_MET_Boson > 1.0);
      if (!passDeltaPhiZMET) nFailCuts++;
      if (nFailCuts > 1) continue;
      fill_passedsteps();

      fillZGhists("ZGaDphiZmet");

      TLorentzVector p4_lljets = boson;

      float min_dphijmet = 4.0;
      float ht_reco_j30 = 0;  // extra variable
      for (auto const &jet : jets) {
        ht_reco_j30 += jet.p4.Pt();
        p4_lljets += jet.p4;
        float dphi = deltaPhi(jet.p4.Phi(), met_p4.Phi());
        if (dphi < min_dphijmet) min_dphijmet = dphi;
      }
      bool passDeltaPhiJetMET = (min_dphijmet > 0.25);
      
      if (!passDeltaPhiJetMET) nFailCuts++;
      if (nFailCuts > 1) continue;
      fill_passedsteps();

      float deltaPhi_MET_lljets = deltaPhi(p4_lljets.Phi(), met_p4.Phi());
      bool passDeltaPhiMETlljets = (deltaPhi_MET_lljets > 2.5);
      if (!passDeltaPhiMETlljets) nFailCuts++;
      if (nFailCuts > 1) continue;
      fill_passedsteps();

      /// Get crude HT value from all jets
      float ht_reco_all = 0;
      for (auto jetpt : Jet_pt()) {
        ht_reco_all += jetpt;
      }

      auto fill_Zmet_hists = [&](string s) {
        if (nFailCuts == 0) {
          plot1d("h_met_"+s, met_p4.Pt(), weight, hvec, ";E_{T}^{miss} [GeV]"  , 120,  0, 600);
          plot1d("h_njets_"+s,  jets.size(), weight, hvec, ";N(jets)"  , 6,  0, 6);
          plot1d("h_ht_reco_all_"+s, ht_reco_all, weight, hvec, ";H_{T}(all jets) [GeV]"  , 80,  0, 800);
          plot1d("h_ht_reco_j30_"+s, ht_reco_j30, weight, hvec, ";H_{T}(good jets) [GeV]"  , 80,  0, 800);
          plot1d("h_boson_pt_"+s, boson.Pt(), weight, hvec, ";p_{T}(boson) [GeV]"  , 120,  0, 600);
          plot1d("h_boson_eta_"+s, boson.Eta(), weight, hvec, ";#eta(boson) "  , 100,  -5.0f, 5.0f);
          plot1d("h_dphiZmet_"+s, deltaPhi_MET_Boson , weight, hvec, ";#Delta#phi(ll, E_{T}^{miss}) ", 32,  0, 3.2);
          plot1d("h_nvtxs_"+s,  PV_npvsGood(), weight, hvec, ";N(vtx good)"  , 80,  0, 80);
        }
      };

      string metsuf = (met_p4.Pt() < 50)? "metlt50_" : (met_p4.Pt() < 125)? "metlt125_" : "metge125_";

      fill_Zmet_hists("fullMET_unwgtd_"+jetcat);
      fill_Zmet_hists(metsuf+"unwgtd_"+jetcat);

      // const bool doNvtxReweight = false;
      const bool doBosonPtReweight = false;
      if (isGamma && !gconf.is_data) {
        // if (doNvtxReweight) {
        //   int nvtx = PV_npvsGood();
        //   if (nvtx > 79) nvtx = 79;
        //   float scale = 1.0;
        //   if (jetcat == "geq1j") scale = nvtxScales_metlt125_geq1j.at(nvtx);
        //   else if (jetcat == "eq0j") scale = nvtxScales_metlt125_eq0j.at(nvtx);
        //   else if (jetcat == "vbf") scale = nvtxScales_metlt125_vbf.at(nvtx);
        //   weight *= scale;
        // }
        if (doBosonPtReweight) {
          // int icat = std::upper_bound(ptRanges.begin(), ptRanges.end(), boson.pt()) - ptRanges.begin() - 1;
          int icat = (boson.Pt() > 440)? 39 : (boson.Pt() - 50) / 10;
          float scale = 1.0;

          if (jetcat == "geq1j") scale = ZptScales_metlt125_geq1j.at(icat);
          else if (jetcat == "eq0j") scale = ZptScales_metlt125_eq0j.at(icat);
          else if (jetcat == "vbf") scale = ZptScales_metlt125_vbf.at(icat);

          weight *= scale;
        }
      }

      fill_Zmet_hists("fullMET_"+jetcat);
      fill_Zmet_hists("fullMET_"+jetcat+"_"+lepcat);
      if (jets.size() >= 2) {
        fill_Zmet_hists("fullMET_geq2j");
        fill_Zmet_hists("fullMET_geq2j_"+lepcat);
      }
      fill_Zmet_hists(metsuf+jetcat);
      fill_Zmet_hists(metsuf+jetcat+"_"+lepcat);
      if (jets.size() >= 2) {
        fill_Zmet_hists(metsuf+"_geq2j_"+lepcat);
      }
      if (is_VBFcat) {
        fill_Zmet_hists(metsuf+"_vbf_"+lepcat);
      }

      fillZGhists("ZGfullMET");
      fillZGhists("ZGfullMET_"+lepcat);

      bool passMETcut = (met_p4.Pt() > 125.);
      if (!passMETcut) nFailCuts++;
      if (nFailCuts > 1) continue;

      fillZGhists("ZGfinal");
      fillZGhists("ZGfinal_"+lepcat);

      // Now fill all the N-1 plots
      if (nFailCuts == 1) {
        auto fillNminus1 = [&](string s) {
          if (!passZmassWindow && (isEE || isMuMu)) plot1d("h_mll"+s, boson.M() , weight, hvec, ";M(ll) [GeV]" , 125,  0, 500);
          else if (!passBosonPtCut && (isEE || isMuMu)) plot1d("h_ptll"+s,  boson.Pt() , weight, hvec, ";p_{T}(ll) [GeV]" , 200,  0, 800);
          // else if (!passMET85) plot1d("h_met"+s, met_p4.Pt(), weight, hvec, ";E_{T}^{miss} [GeV]"  , 160,  0, 800);
          else if (!passMETcut) plot1d("h_met"+s, met_p4.Pt(), weight, hvec, ";E_{T}^{miss} [GeV]"  , 160,  0, 800);
          else if (!passDeltaPhiJetMET) plot1d("h_min_dphijmet"+s, min_dphijmet , weight, hvec, ";min #Delta#phi(j, E_{T}^{miss}) ", 32,  0, 3.2);
          else if (!passDeltaPhiZMET) plot1d("h_dphiZmet"+s, deltaPhi_MET_Boson , weight, hvec, ";#Delta#phi(ll, E_{T}^{miss}) ", 32,  0, 3.2);
          else if (!passDeltaPhiMETlljets) plot1d("h_dphiMETlljets"+s, deltaPhi_MET_lljets , weight, hvec, ";#Delta#phi(ll+jets, E_{T}^{miss}) ", 32,  0, 3.2);
        };
        fillNminus1("");
        fillNminus1("_"+jetcat);
        fillNminus1("_"+lepcat);
        fillNminus1("_"+jetcat+"_"+lepcat);
        if (is_VBFcat) {
          fillNminus1(metsuf+"_vbf_"+lepcat);
        }
        // if (jets.size() > 0) fillNminus1("_geq1j");

      }

      if (nFailCuts > 0) continue;
      // Events passes all sphotontions

      float mtZZ = getDileptonMT(boson, met_p4);

      auto fillhists = [&](string s) {
        plot1d("h_njets"+s,  jets.size(), weight, hvec, ";N(jets)"  , 6,  0, 6);
        plot1d("h_met"+s,    met_p4.Pt()  , weight, hvec, ";E_{T}^{miss} [GeV]"  , 30,  0, 750);
        plot1d("h_metphi"+s, met_p4.Phi() , weight, hvec, ";#phi(E_{T}^{miss})"  , 34, -3.4, 3.4);
        // plot1d("h_metsig"+s, met_sig   , weight, hvec, ";#sigma(E_{T}^{miss}) " , 20,  0, 40);
        // Z quantities
        if (isEE || isMuMu) {
          plot1d("h_mll"+s,    boson.M() , weight, hvec, ";M_{ll} [GeV]"        , 125,  0, 500);
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

        plot1d("h_mtZZ"+s, mtZZ  , weight, hvec, ";M_{T}(ll) [GeV]"        , 60,  0, 1500);
        // const vector<float> mtbin1 = {0, 75, 150, 225, 300, 375, 450, 525, 600, 725, 850, 975, 1100, 1350, 1600, 2100, 3000};
        // plot1d("h_mtll_b1"+s,   mtll , weight, hvec, ";M_{T}(ll) [GeV]" , mtbin1.size()-1, mtbin1.data());
        // const vector<float> mtbin3 = {0, 75, 150, 225, 300, 375, 450, 525, 600, 725, 850, 975, 1100, 1350, 1600, 2100, 3000};
        // plot1d("h_mtll_b3"+s,   mtll , weight, hvec, ";M_{T}(ll) [GeV]" , mtbin3.size()-1, mtbin3.data());

      };

      fillhists("");
      fillhists("_"+jetcat);
      fillhists("_"+lepcat);
      fillhists("_"+jetcat+"_"+lepcat);
      if (is_VBFcat) {
        fillhists(metsuf+"_vbf_"+lepcat);
      } else if (jets.size() >= 1) {
        fillhists(metsuf+"_geq1j_"+lepcat);
      }

      if (isGamma && met_p4.Pt() > 500 && jetcat == "eq0j") {
        cout << __LINE__ << ": met_p4.Pt()= " << met_p4.Pt() << ", boson.Pt()= " << boson.Pt() << ", eta_boson()= " << boson.Eta() << ", run()= " << run() << ", luminosityBlock()= " << luminosityBlock() << ", event()= " << event() << endl;
      }

    } // Event loop

    delete file;
  } // File loop
  if (show_progress) bar.finish();

  plot1d("h_sum_wgts", 0, nEventsTotal, hvec, ";Bin;Sum of weights" , 5, 0, 5);
  plot1d("h_sum_wgts", 1, sum_effevts, hvec, ";Bin;Sum of weights" , 5, 0, 5);
  plot1d("h_sum_wgts", 2, sum_genwgts, hvec, ";Bin;Sum of weights" , 5, 0, 5);
  // plot1d("h_sum_wgts", 3, sum_nwgt, hvec, ";Bin;Sum of weights" , 5, 0, 5);

  fout->cd();
  for (const auto& h : hvec) {
    if (h.first.find("hnum") != 0) continue;
    if (h.first.find("2d") == string::npos)
      moveOverFlowToLastBin1D(h.second);
    string hname = h.first;
    hname.erase(0, 4);
    // dummy.cd();
    TH1F* h_ratio = (TH1F*) h.second->Clone(("ratio"+hname).c_str());
    h_ratio->SetDirectory(0);
    h_ratio->Divide(h_ratio, hvec.at("hden"+hname), 1, 1, "B");

    const string dirname = "effstudy_eff";
    TDirectory* dir = (TDirectory*) fout->Get(dirname.c_str());
    if (dir == nullptr) dir = fout->mkdir(dirname.c_str());
    dir->cd();
    h_ratio->Write();

    dir = (TDirectory*) fout->Get("effstudy_num");
    if (dir == nullptr) dir = fout->mkdir("effstudy_num");
    dir->cd();
    h.second->Write();

    dir = (TDirectory*) fout->Get("effstudy_den");
    if (dir == nullptr) dir = fout->mkdir("effstudy_den");
    dir->cd();
    hvec.at("hden"+hname)->Write(("hden"+hname).c_str());
  }

  for (auto& h : hvec) {
    if (h.first.find("hnum") == 0 || h.first.find("hden") == 0) continue;
    if (h.first.find("phi") == string::npos && h.first.find("2d") == string::npos)
      moveOverFlowToLastBin1D(h.second);
    string dirname = "hzz2l2nu";
    vector<string> jetsufs = {"_eq0j", "_eq1j", "_eq2j", "_1j20", "_vbf", "_eq0j_lowdphi", };
    vector<string> lepsufs = {"_gamma", "_ee", "_mumu", "_emu", "_ll", "_llg"};
    vector<string> dirsufs;
    for (string jsuf : jetsufs) {
      for (string lsuf : lepsufs) {
        dirsufs.push_back(jsuf+lsuf);
      }
    }
    dirsufs.insert(dirsufs.end(), jetsufs.begin(), jetsufs.end());
    dirsufs.insert(dirsufs.end(), lepsufs.begin(), lepsufs.end());

    for (string dsuf : dirsufs) {
      if (TString(h.first).EndsWith(dsuf.c_str())) {
        dirname += dsuf;
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
  cout << nEventsTotal << " Events Processed, where " << int(sum_effevts) << " effective events recored, " << endl;
  cout << nDuplicates << " duplicates were skipped, and " << nPassedTotal << " Events passed all selections." << endl;
  cout << "---------------------------------\n" << endl;

  return 0;
}
