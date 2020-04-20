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

const bool applySyncCuts = false;

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
          cout << __LINE__ << ": gconf.lumi= " << gconf.lumi << ", nevtsamp= " << nevtsamp << ", scaleToLumi= " << scaleToLumi << endl;
          exit(4);
        }
      }
    }

    muoncorr->reset();
    muoncorr->init(Form("config/RocesterCorrection/RoccoR%d.txt", gconf.year));

    // auto psRead = new TTreePerfStats("readPerf", tree);

    nt.Init(tree);

    for( unsigned int evt = 0; evt < tree->GetEntriesFast(); ++evt) {

      nt.GetEntry(evt);
      tree->LoadTree(evt);

      nEventsTotal++;
      if (!gconf.is_data) {
        sum_effevts += (genWeight() > 0)? 1 : -1;
        sum_genwgts += genWeight();
      }

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

      float weight = (gconf.is_data)? 1 : (genWeight() > 0)? scaleToLumi : -1.0*scaleToLumi;

      // The following skim cuts are in the ntuple skiming steps
      // 2l2nu cut: nTightElection + nTightMuon >= 2, Z pt > 50, mZ > 50


      plot1d("h_preselec_steps", 0, 1, hvec, ";M_{ll} [GeV]" , 20,  0, 20);

      bool passMETfilt = passesMETfilters(gconf.is_data);
      if (!passMETfilt) continue;

      plot1d("h_preselec_steps", 1, 1, hvec, ";M_{ll} [GeV]" , 20,  0, 20);

      plot1d("h_preselec_steps", 2, 1, hvec, ";M_{ll} [GeV]" , 20,  0, 20);

      // To apply the same skim cut as the NanoAOD
      if (applySyncCuts) {
        // if (isGamma)
        //   if (!InstrMETPrunerCuts()) continue;
        // else
        if (!ZZ2l2vPrunerCuts()) continue;
        // if (!InstrMETPrunerCuts()) continue;
      }

      plot1d("h_preselec_steps", 3, 1, hvec, ";M_{ll} [GeV]" , 20,  0, 20);

      nPassedTotal++;

      int nFailCuts = 0;

      auto [tightElectrons, looseElectrons] = getElectrons();
      auto [tightMuons, looseMuons] = getMuons();
      auto photons = getPhotons();

      bool isEE = (tightElectrons.size() == 2); //2 good Photon
      bool isMuMu = (tightMuons.size() == 2); //2 good muons
      bool isGamma = (photons.size() == 1 && !isEE && !isMuMu); //1 good photon

      string lepcat = "gamma";
      if (isEE) lepcat = "ee";
      else if (isMuMu) lepcat = "mumu";
      else if (isGamma) lepcat = "gamma";

      int istep=0;
      plot1d("h_weight_"+to_string(istep), weight, 1, hvec, ";Event weight"  , 200,  0, 20);
      plot1d("h_met_step"+to_string(istep), MET_pt(), weight, hvec, ";E_{T}^{miss} [GeV]"  , 120,  0, 600);
      plot1d("h_met_unwgtd_step"+to_string(istep), MET_pt(), 1, hvec, ";E_{T}^{miss} [GeV]"  , 120,  0, 600);

      auto fill_passedsteps = [&]() {
        if (nFailCuts == 0) {
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

      if (gconf.is_data) {
        if (isGamma) {
          double prescale = getPhotonTrigPrescale(photons[0].p4.Pt());
          if (prescale <= 0.0) continue;
          weight *= prescale;
        }
        else if (isEE && !passTriggerSelections(2))
          continue;
        else if (isMuMu && !passTriggerSelections(1))
          continue;
      } else {

      }

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
      if (isGamma) {
        boson = photons[0].p4;
        giveMassToPhoton(boson);
      } else {
        boson = tightLeptons[0].p4 + tightLeptons[1].p4;
      }

      if (isGamma && !gconf.is_data) {
        // Apply the LO to NLO scale factor here
        weight *= (boson.Pt() >= 600)? 1 : (1.72 - 0.0012 * boson.Pt());
      }

      // vector<Particle*> isoobjs;
      auto jets = getJets(looseMuons, looseElectrons, photons);

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

      // TLorentzVector puppimet_p4;
      // puppimet_p4.SetPtEtaPhiM(PuppiMET_pt(), 0., PuppiMET_phi(), 0.);
      // float met_sig = MET_significance();

      auto fill_mll_hists = [&](string s) {
        if (!isMuMu && !isEE) return;
        plot1d("h_mll_"+s, boson.M() , scaleToLumi, hvec, ";M_{ll} [GeV]" , 125,  0, 500);
        plot1d(Form("h_mll_%s_%s", s.c_str(), lepcat.c_str()), boson.M() , scaleToLumi, hvec, ";M_{ll} [GeV]" , 125,  0, 500);
        plot1d(Form("h_mll_%s_%s_%s", s.c_str(), jetcat.c_str(), lepcat.c_str()), boson.M() , scaleToLumi, hvec, ";M_{ll} [GeV]" , 125,  0, 500);
      };
      fill_mll_hists("full");

      bool passBosonPtCut = (boson.Pt() >= 55.);
      if (!passBosonPtCut) nFailCuts++;
      if (nFailCuts > 1) continue;
      fill_passedsteps();

      // fill_mll_hists("Zpt55");

      bool passZmassWindow = (fabs(boson.M() - mZ) < 15);
      if (!passZmassWindow) nFailCuts++;
      if (nFailCuts > 1) continue;
      fill_passedsteps();

      float deltaPhi_MET_Boson = deltaPhi(boson.Phi(), MET_phi());
      bool passDeltaPhiZMET = (deltaPhi_MET_Boson > 1.0);
      if (!passDeltaPhiZMET) nFailCuts++;
      if (nFailCuts > 1) continue;
      fill_passedsteps();

      TLorentzVector p4_lljets = boson;

      float min_dphijmet = 4.0;
      for (auto const &jet : jets) {
        p4_lljets += jet.p4;
        float dphi = deltaPhi(jet.p4.Phi(), MET_phi());
        if (dphi < min_dphijmet) min_dphijmet = dphi;
      }
      bool passDeltaPhiJetMET = (min_dphijmet > 0.25);
      
      if (!passDeltaPhiJetMET) nFailCuts++;
      if (nFailCuts > 1) continue;
      fill_passedsteps();

      float deltaPhi_MET_lljets = deltaPhi(p4_lljets.Phi(), MET_phi());
      bool passDeltaPhiMETlljets = (deltaPhi_MET_lljets > 2.5);
      if (!passDeltaPhiMETlljets) nFailCuts++;
      if (nFailCuts > 1) continue;
      fill_passedsteps();

      /// Get crude HT value from all jets
      float ht_all = 0;
      for (auto jetpt : Jet_pt()) {
        ht_all += jetpt;
      }

      auto fill_Zmet_hists = [&](string s) {
        if (nFailCuts == 0) {
          plot1d("h_met_"+s, met_p4.Pt(), weight, hvec, ";E_{T}^{miss} [GeV]"  , 120,  0, 600);
          plot1d("h_njets_"+s,  jets.size(), weight, hvec, ";N(jets)"  , 6,  0, 6);
          plot1d("h_ht_all_"+s, ht_all, weight, hvec, ";H_{T}(all jets) [GeV]"  , 80,  0, 800);
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
      const bool doBosonPtReweight = true;
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

      bool passMETcut = (MET_pt() > 125.);
      if (!passMETcut) nFailCuts++;
      if (nFailCuts > 1) continue;

      // Now fill all the N-1 plots
      if (nFailCuts == 1) {
        auto fillNminus1 = [&](string s) {
          if (!passZmassWindow && (isEE || isMuMu)) plot1d("h_mll"+s, boson.M() , weight, hvec, ";M(ll) [GeV]" , 125,  0, 500);
          else if (!passBosonPtCut && (isEE || isMuMu)) plot1d("h_ptll"+s,  boson.Pt() , weight, hvec, ";p_{T}(ll) [GeV]" , 200,  0, 800);
          // else if (!passMET85) plot1d("h_met"+s, met_p4.Pt(), weight, hvec, ";E_{T}^{miss} [GeV]"  , 160,  0, 800);
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

      // if (!passTrackVeto) continue;

      float mtZZ = getDileptonMT(boson, met_p4);

      auto fillhists = [&](string s) {
        plot1d("h_njets"+s,  jets.size(), weight, hvec, ";N(jets)"  , 6,  0, 6);
        plot1d("h_met"+s,    MET_pt()  , weight, hvec, ";E_{T}^{miss} [GeV]"  , 30,  0, 750);
        plot1d("h_metphi"+s, MET_phi() , weight, hvec, ";#phi(E_{T}^{miss})"  , 34, -3.4, 3.4);
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

      // float weight = genWeight();

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
    for (string dirsuf : {
        "_eq0j_ee", "_eq0j_mumu", "_geq1j_ee", "_geq1j_mumu",
            "_eq1j_ee", "_eq1j_mumu", "_eq2j_ee", "_eq2j_mumu", "_vbf_ee", "_vbf_mumu",
            "_eq0j_gamma", "_geq1j_gamma", "_eq1j_gamma", "_eq2j_gamma", "_vbf_gamma", 
            "_geq2j_ee", "_geq2j_mumu", "_geq2j_gamma",
            "_ee", "_mumu", "_gamma", "_eq0j", "_eq1j", "_eq2j", "_geq1j", "_vbf", "_geq2j",
            }) {
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
  cout << nEventsTotal << " Events Processed, where " << int(sum_effevts) << " effective events recored, " << endl;
  cout << nDuplicates << " duplicates were skipped, and " << nPassedTotal << " Events passed all selections." << endl;
  cout << "---------------------------------\n" << endl;

  return 0;
}
