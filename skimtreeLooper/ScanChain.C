#pragma GCC diagnostic ignored "-Wsign-compare"

#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1.h"
#include "TChain.h"
#include "TSystem.h"

#include "SkimTree.h"
#include "Utilities.h"

using namespace std;
using namespace tas;

bool doBosonPtReweight = false;
bool extendEEphoton2j = false;
bool doNvtxReweight = false;

const bool produceResultTree = true;
const bool blind = true;

vector<float> nvtxscale(100, 0);

int ScanChain(TString indir, TString sample, TString tag){

  TChain *ch = new TChain("SkimTree");
  TString files_in = Form("%s/%s*.root", indir.Data(), sample.Data());
  ch->Add(files_in);
  cout << ">> Adding " << files_in << " into the chain." << endl;

  TString samplever = tag(0, min(tag.Index("_",tag.Index("_")+1),5));
  double skimver = TString(samplever).ReplaceAll("v","").ReplaceAll("_",".").Atof();
  cout << ">> Runing from skimver " << samplever  << ", i.e. " << skimver << endl;

  TString outname = sample;
  TString fileout = Form("output/%s/%s.root", tag.Data(), outname.Data());
  gSystem->Exec(Form("mkdir -p output/%s", tag.Data()));
  TFile* fout = new TFile(fileout, "RECREATE");
  cout << ">> Outputting to " << fileout << endl;

  if (tag.Contains("nvtxrwgt")) {
    doNvtxReweight = true;
    cout << ">> doNvtxReweight is set to true." << endl;
  }
  else if (tag.Contains("rwgtd")) {
    doBosonPtReweight = true;
    cout << ">> doBosonPtReweight is set to true." << endl;
  }
  extendEEphoton2j = tag.Contains("ee2j");

  bool is_data = (sample.BeginsWith("201")) || sample.BeginsWith("photon") || sample.BeginsWith("egamma");
  bool is_gjets = sample.BeginsWith("GJets");

  int year = (tag.Contains("2016"))? 2016 : (tag.Contains("2017"))? 2017 : 2018;

  map<string,TH1*> hvec;

  unsigned nEventsTotal = 0;
  unsigned nEventsChain = ch->GetEntries();

  TTree* tout = nullptr;
  float evt_weight;
  float ll_pt, ll_eta, ll_phi, ll_mass;
  float ptmiss, ptmiss_phi, mT, M_ZZ;
  float dphijmet, dphiVmet, dphilljmet;
  int lepton_cat, jet_cat, njet, nvtxs;
  int njet_tightPUjetId, njet_medPUjetId, njet_loosePUjetId;
  unsigned run, lumi, event;

  bool photon_convveto;
  bool photon_isEB;
  bool photon_id_cutbase;
  bool photon_id_mva90;
  bool photon_id_mva80;

  float photon_pt, photon_eta, photon_phi;
  float photon_r9, photon_sieie, photon_sipip;
  float photon_Emip, photon_seedtime, photon_e4oe1;

  if (produceResultTree) {
    tout = new TTree("SkimTree", "");

    tout->Branch("weight", &evt_weight);
    tout->Branch("ll_pt", &ll_pt);
    tout->Branch("ll_eta", &ll_eta);
    tout->Branch("ll_phi", &ll_phi);
    tout->Branch("ll_mass", &ll_mass);
    tout->Branch("ptmiss", &ptmiss);
    tout->Branch("ptmiss_phi", &ptmiss_phi);
    tout->Branch("mT", &mT);
    tout->Branch("lepton_cat", &lepton_cat);
    tout->Branch("jet_cat", &jet_cat);

    tout->Branch("mindphi_jet_met", &dphijmet);
    tout->Branch("dphi_boson_met", &dphiVmet);
    tout->Branch("dphi_lljets_met", &dphilljmet);
    tout->Branch("mZZ", &M_ZZ);
    tout->Branch("njet", &njet);
    tout->Branch("nvtxs", &nvtxs);
    tout->Branch("run", &run);
    tout->Branch("lumi", &lumi);
    tout->Branch("event", &event);

    tout->Branch("photon_convveto", &photon_convveto);
    tout->Branch("photon_isEB", &photon_isEB);
    tout->Branch("photon_id_cutbase", &photon_id_cutbase);
    tout->Branch("photon_id_mva80", &photon_id_mva80);

    tout->Branch("photon_pt", &photon_pt);
    tout->Branch("photon_eta", &photon_eta);
    tout->Branch("photon_phi", &photon_phi);
    tout->Branch("photon_r9", &photon_r9);
    tout->Branch("photon_sieie", &photon_sieie);
    tout->Branch("photon_sipip", &photon_sipip);
    tout->Branch("photon_Emip", &photon_Emip);
    tout->Branch("photon_seedtime", &photon_seedtime);
  }

  // Make photon-pt bins for reweight
  vector<float> ptbin0;
  for (int i = 55; i < 160; i += 5) ptbin0.emplace_back(i);
  for (int i = 160; i <= 200; i += 10) ptbin0.emplace_back(i);
  ptbin0.emplace_back(250); ptbin0.emplace_back(450);

  vector<float> ptbin1;
  for (int i = 55; i < 160; i += 5) ptbin1.emplace_back(i);
  for (int i = 160; i < 300; i += 10) ptbin1.emplace_back(i);
  for (int i = 300; i <= 350; i += 25) ptbin1.emplace_back(i);
  ptbin1.emplace_back(400); ptbin1.emplace_back(600);

  // Setup pileup re-weighting for comparing data of different years
  if (doNvtxReweight) {
    TFile f_purw("../offshellLooper/input/nvtx_reweighting_outhist_to16.root");
    TString scaletype = "17to16";
    if (year == 2018) scaletype = "18to16";
    else if (year == 2017) scaletype = "17to16";
    // else if (samplestr.find("data_2017F_09May") != string::npos) scaletype = "FToBtoE";
    TH1F* h_nvtxscale = (TH1F*) f_purw.Get("h_nvtxscale_"+scaletype);
    if (!h_nvtxscale) throw invalid_argument("No nvtx reweighting hist found for " + scaletype);
    cout << "Doing nvtx reweighting! Scaling " << scaletype << endl;
    for (int i = 1; i < h_nvtxscale->GetNbinsX(); ++i) {
      nvtxscale[i] = h_nvtxscale->GetBinContent(h_nvtxscale->FindBin(i));
    }
  }

  cout << ">> To run over nEventsChain = " << nEventsChain << endl;

  TFile *currentFile = 0;
  TObjArray *listOfFiles = ch->GetListOfFiles();
  TIter fileIter(listOfFiles);

  while ( (currentFile = (TFile*)fileIter.Next()) ) { 

    TFile *file = new TFile( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("SkimTree");
    st.Init(tree);

    TString filename(currentFile->GetTitle());

    for( unsigned int evt = 0; evt < tree->GetEntriesFast(); ++evt) {
      
      st.GetEntry(evt);
      nEventsTotal++;

      SkimTree::progress(nEventsTotal, nEventsChain);

      double weight = event_wgt();
      if (is_data) weight = 1;

      bool isllg = (is_mumu() || is_ee()) && (event_Nphotons() == 1);
      bool isgamma = is_gamma();
      // bool isgamma = event_Nphotons() == 1 && !is_mumu() && !is_ee();
      // if (extendEEphoton2j) isgamma = is_gamma() || (event_Njets() >= 2 && event_Nphotons() == 1 && !is_mumu() && !is_ee());
      if (extendEEphoton2j) isgamma = is_gamma() || (event_Njets() >= 2 && event_Nphotons() == 1 && convveto_photon() && !is_mumu() && !is_ee());

      double wgt = weight;
      auto fillPhotonHitMap = [&](string suf) {
        plot2d("h2d_photon_eta_phi"+suf, eta_photon(), phi_photon(), wgt, hvec, ";#eta(#gamma);#phi(#gamma)" , 50, -2.5f, 2.5f, 64, -3.2f, 3.2f);
        plot1d("h_photon_eta"+suf, eta_photon(), wgt, hvec, ";#eta(#gamma)" , 50, -2.5f, 2.5f);
        plot1d("h_photon_pt"+suf, pt_photon(), wgt, hvec, ";p_{T}($gamma) [GeV]" , 160, 0, 800);
        plot1d("h_met_phm"+suf, pTmiss(), wgt, hvec, ";p_{T}^{miss} [GeV]" , 160, 0, 800);
        if (skimver >= 2.1) {
          plot1d("h_photon_sieie"+suf, sieie_photon(), weight, hvec, ";#sigma_{i#etai#eta}(#gamma)" , 80,  0, 0.02);
          plot1d("h_photon_sipip"+suf, sipip_photon(), weight, hvec, ";#sigma_{i#phii#phi}(#gamma)" , 80,  0, 0.02);
        }
        if (fabs(eta_photon()) > 1.442) {
          plot1d("h_photon_endcap_phi"+suf, phi_photon(), wgt, hvec, ";#phi(#gamma)" , 64, -3.2f, 3.2f);
        }
      };

      if (skimver >= 1.11 && event_Nphotons() == 1) {
        wgt *= event_wgt_trig_photon();
        fillPhotonHitMap("_raw");
        if (convveto_photon())
          fillPhotonHitMap("_raw_convveto");
        if (id_photon_Hgg())
          fillPhotonHitMap("_raw_Hggid");
        if (convveto_photon() && id_photon_Hgg())
          fillPhotonHitMap("_raw_convAndHggId");
      }
      if (skimver == 2.02) {
        if (sample == "ZGJets_nunu_pTG_130" )
          weight *= 0.1926 / 0.2828;  // "ZGJets_nunu_nlo_pTG_130-inf"
        if (sample == "ZGJets_nunu_nlo_incl" )
          weight *= ((year == 2016)? 27.98 : 30.05) / 5.485;   // "ZGJets_ll_pTG_40-130"
        if (sample == "ZGJets_nunu_nlo_pTG_130" )
          weight *= ((year == 2016)? 0.278 : 0.2828) / 0.1472;  // "ZGJets_ll_pTG_130-inf"
        if (sample == "ZGJets_ll_pTG_40-130" )
          weight *= 5.485 / 0.1595;  // "ZGJets_ll_nlo_pTG_130-inf"
      }

      if ((sample == "ZGJets_nunu_nlo_incl" || sample == "ZGJets_ll_nlo_incl") && pt_photon() > 135) continue;
      if ((sample == "ZGJets_nunu_nlo_pTG_130" || sample == "ZGJets_ll_nlo_pTG_130") && pt_photon() < 135) continue;
      // if (sample == "ZGJets_ll_nlo_incl") weight *= 2.56;
      if (sample == "ZGJets_ll_nlo_incl") weight *= 0.8;

      if (isgamma) weight *= event_wgt_trig_photon();
      // if (isllg || isgamma) weight *= event_wgt_trig_photon();
      else if (is_ee()) weight *= event_wgt_trig_electron();
      else if (is_mumu()) weight *= event_wgt_trig_muon();
      // else if (is_emu()) weight *= ((event_wgt_trig_muon() + event_wgt_trig_electron()) < 0.001)? 1.0 : 0.0;

      string lepcat = "_lepcat";
      if (isllg) lepcat = "_llg";
      else if (isgamma) lepcat = "_gamma";
      else if (is_ee()) lepcat = "_ee";
      else if (is_mumu()) lepcat = "_mumu";
      else if (is_emu()) lepcat = "_emu";

      string possuf = (fabs(eta_boson()) < 1.4442)? "_barrel" : "_endcap";

      if (is_data && isgamma && event_Njets() == 0 && fabs(eta_boson()) > 1.4442) {
        if (fabs(phi_boson()) < 0.8 || fabs(phi_boson()) > 2.6) continue;
        else weight *= 3.14159 / 1.6;
      }

      if (!is_data && weight >= 1e6) {
        cout << __LINE__ << ": weight= " << weight << ", event= " << event << ", sample= " << sample << ", currentFile->GetTitle() = " << currentFile->GetTitle() << endl;
        continue;
      }
      
      if ((is_ee() || is_mumu()) && fabs(mass_boson() - 91.2) > 15) continue;
      if (!isllg) {
        if (dphi_boson_met() < 1.0) continue;
        if (dphi_lljets_met() < 2.5) continue;
        if (event_Njets() > 0) {
          if (mindphi_jet_met() < 0.25) continue;
        }
        if (pt_boson() < 55.) continue;
        if (fabs(mass_boson() - 91.2) > 15.) continue;
        if (isgamma && !convveto_photon()) continue;  // TODO: to confirm the effect of this cut
        if (isgamma && skimver >= 2.1 && (sieie_photon() < 0.001 || sipip_photon() < 0.001)) continue;
        if (isgamma && skimver >= 2.1 && (mipE_photon() > 4.9 || fabs(seedtime_photon()) > 1.2)) continue;
        // if (isgamma && mipE_photon() > 4.9) continue;
        if (blind && is_data && (is_ee() || is_mumu()) && pTmiss() >= 125) continue;
      }

      if (doNvtxReweight && year > 2016) {
        if (event_nvtxs_good() < 100) weight *= nvtxscale[event_nvtxs_good()];  // only scale for data
        plot1d("h_nvtxs_rwtd", event_nvtxs_good(), weight, hvec, ";Number of vertices", 100, 1, 101);
      }

      float met = pTmiss();
      float metphi = phimiss();
      M_ZZ = mZZ();
      float mtZZ = mTZZ();
      float mboson = mass_boson();
      float Vpt  = pt_boson();
      float Veta = eta_boson();
      float Vphi = phi_boson();
      dphijmet = mindphi_jet_met();
      dphiVmet = dphi_boson_met();
      dphilljmet = dphi_lljets_met();
      float dphiVjet = 4.0;

      if (isllg) {
        if (pt_photon() < 55.) continue;
        if (fabs(eta_photon()) > 1.4442) continue;
        // if (!convveto_photon()) continue;
        // And a lot of other cuts
        LorentzVector met_p4(pTmiss(), 0, phimiss(), 0);
        met_p4 += LorentzVector(pt_boson(), eta_boson(), phi_boson(), 0);
        LorentzVector boson(pt_photon(), eta_photon(), phi_photon(), mass_boson());
        // Only to look at njet == 0 now, where it matters
        dphiVmet = deltaPhi(met_p4.phi(), boson.phi());
        if (dphiVmet < 1.0) continue;
        if (fabs(mass_boson() - 91.2) > 15.) boson.SetM(91.2);

        LorentzVector lljets = boson;
        if (event_Njets() > 0) {
          dphijmet = deltaPhi(phi_jet1(), met_p4.phi());
          lljets += LorentzVector(pt_jet1(), eta_jet1(), phi_jet1(), 0);
        }
        if (event_Njets() > 1) {
          dphijmet = min(dphijmet, deltaPhi(phi_jet2(), met_p4.phi()));
          lljets += LorentzVector(pt_jet2(), eta_jet2(), phi_jet2(), 0);
        }
        dphilljmet = deltaPhi(lljets.phi(), met_p4.phi());

        // to keep dphiVmet 1.0 - 2.5 as x-check region
        if (event_Njets() > 0) {
          if (dphilljmet < 2.5) continue;
          if (dphijmet < 0.25) continue;
        }

        met = met_p4.pt();
        metphi = met_p4.phi();
        Vpt = boson.pt();
        Veta = boson.eta();
        Vphi = boson.phi();
        mboson = boson.M();
        mtZZ = getDileptonMT(boson, met_p4);
        LorentzVector Zmiss_approx_p4(pTmiss(), Veta, phimiss(), 91.2);
        M_ZZ = (Zmiss_approx_p4+boson).M();

      }

      if (event_Njets() > 0) dphiVjet = deltaPhi(phi_jet1(), Vphi);          
      if (event_Njets() > 1) dphiVjet = min(dphiVjet, deltaPhi(phi_jet2(), Vphi));

      if (skimver >= 1.11 && event_Nphotons() == 1) {
        wgt = weight;
        fillPhotonHitMap("_Vcuts");
        if (convveto_photon())
          fillPhotonHitMap("_Vcuts_convveto");
        if (id_photon_Hgg())
          fillPhotonHitMap("_Vcuts_Hggid");
        if (convveto_photon() && id_photon_Hgg())
          fillPhotonHitMap("_Vcuts_convAndHggId");
      }

      // if (!pass_lepveto()) continue;
      if ((is_ee() || is_mumu() || is_gamma())) {
        if (!pass_lepveto()) continue;
      } else if (event_Nphotons() == 1) {
      
      }
      if (!pass_trackveto()) continue;
      if (year == 2018 && has_ak4jets_inHEM1516()) continue;
      if (event_Njets_btagged() > 0) continue;

      string jetcat = "_jetcat";
      if (event_Njets() == 0) jetcat = "_eq0j";
      else if (event_Njets() == 1) jetcat = "_eq1j";
      else if (event_Njets() == 2) jetcat = "_eq2j";
      else if (event_Njets() >= 3) jetcat = "_ge3j";

      if (event_Njets() == 0 && dphiVmet < 2.5) jetcat = "_eq0j_lowdphi";

      if (isgamma) {
        // if (fabs(eta_boson()) > 1.4442) continue;  // TODO: to remove this line later
        if (is_gjets) {
          weight *= std::max(1., 1.716910-0.001221*pt_boson());
        }
        int icat = (pt_boson() > 440)? 39 : (pt_boson() - 50) / 10;
        int icat0 = std::upper_bound(ptbin0.begin(), ptbin0.end(), pt_boson()) - ptbin0.begin() - 1;
        int icat1 = std::upper_bound(ptbin1.begin(), ptbin1.end(), pt_boson()) - ptbin1.begin() - 1;

        if (doBosonPtReweight) {
          float scale = 1.0;
          if (year == 2018) {
            if (is_data) {
              if      (jetcat == "_eq1j") scale = sf_Vpt_data18_metlt80_eq1j.at(icat);
              else if (jetcat == "_eq0j") scale = sf_Vpt_data18_metlt80_eq0j.at(icat);
              else if (jetcat == "_eq2j") scale = sf_Vpt_data18_metlt80_eq2j.at(icat);
            } else {
              if      (jetcat == "_eq1j") scale = sf_Vpt_closure18_metlt80_eq1j.at(icat);
              else if (jetcat == "_eq0j") scale = sf_Vpt_closure18_metlt80_eq0j.at(icat);
              else if (jetcat == "_eq2j") scale = sf_Vpt_closure18_metlt80_eq2j.at(icat);
            }
          } else if (year == 2016) {
            bool DYclosureTest = false;
            if (is_data || !DYclosureTest) {
              cout << __LINE__ << ": icat0= " << icat0 << ", pt_boson()= " << pt_boson() << ", ptbin0.size()= " << ptbin0.size() << endl;
              if (extendEEphoton2j && jetcat == "_eq2j") {
                int ibin = (fabs(eta_boson()) >= 2.5)? 0 : int((eta_boson() + 2.5) / 0.1);
                scale = sf_Veta_data16_metlt80_ee2j.at(ibin);
                scale *= sf_Vpt_flateta_data16_metlt80_ee2j.at(icat);
              }
              else if (jetcat == "_eq0j") scale = sf_Vpt_data16_metlt125_eq0j.at(icat0);
              else if (jetcat == "_eq1j") scale = sf_Vpt_data16_metlt125_eq1j.at(icat0);
              else if (jetcat == "_eq2j") scale = sf_Vpt_data16_metlt125_eq2j.at(icat0);
            } else {
              if (extendEEphoton2j && jetcat == "_eq2j") scale = sf_Vpt_closure16_metlt125_ee2j.at(icat);
              else if (jetcat == "_eq0j") scale = sf_Vpt_closure_metlt80_eq0j.at(icat);
              else if (jetcat == "_eq1j") scale = sf_Vpt_closure_metlt80_eq1j.at(icat);
              else if (jetcat == "_eq2j") scale = sf_Vpt_closure_metlt80_eq2j.at(icat);
            }
          }
          weight *= scale;
        }
      } else {
      }

      if (skimver >= 1.11 && event_Nphotons() == 1) {
        wgt = weight;
        fillPhotonHitMap("_fullMET");
        if (convveto_photon())
          fillPhotonHitMap("_fullMET_convveto");
        if (id_photon_Hgg())
          fillPhotonHitMap("_fullMET_Hggid");
        if (convveto_photon() && id_photon_Hgg())
          fillPhotonHitMap("_fullMET_convAndHggId");
        if (met > 125) {
          fillPhotonHitMap("_final");
          if (convveto_photon())
            fillPhotonHitMap("_final_convveto");
          if (isgamma)
            fillPhotonHitMap("_final_isgamma");
        }
      }

      auto fillhists = [&](string s = "") {
        plot1d("h_met"+s, met, weight, hvec, ";E_{T}^{miss} [GeV]"  , 160,  0, 800);
        plot1d("h_metphi"+s, metphi, weight, hvec, ";#phi(E_{T}^{miss})"  , 68, -3.4, 3.4);

        plot1d("h_mZZ"+s, M_ZZ, weight, hvec, ";m(ZZ) [GeV]"  , 160,  0, 800);
        plot1d("h_mtZZ"+s, mtZZ, weight, hvec, ";m_{T}(ZZ) [GeV]"  , 160,  0, 800);
        plot1d("h_njets"+s, event_Njets(), weight, hvec, ";N(jets)"  , 6,  0, 6);
        plot1d("h_njets20"+s, event_Njets20(), weight, hvec, ";N(jets)"  , 6,  0, 6);
        plot1d("h_nvtxs"+s, event_nvtxs_good(), weight, hvec, ";N(vtxs)"  , 80,  0, 80);

        plot1d("h_boson_mass"+s, mboson , weight, hvec, ";M_{ll} [GeV]" , 125,  0, 500);
        plot1d("h_boson_mass_finebin"+s, mboson , weight, hvec, ";M_{ll} [GeV]" , 160,  75, 107);
        plot1d("h_boson_pt"+s, Vpt, weight, hvec, ";p_{T}(boson) [GeV]" , 160, 0, 800);
        plot1d("h_boson_eta"+s, Veta, weight, hvec, ";#eta(boson)" , 100,  -5.f, 5.f);
        plot1d("h_boson_phi"+s, Vphi, weight, hvec, ";#phi(boson)" , 64,  -3.2f, 3.2f);
        string possuf = (fabs(Veta) < 1.4442)? "_barrel" : "_endcap";
        plot1d("h_boson_phi"+possuf+s, phi_boson(), weight, hvec, ";#phi(boson)" , 64,  -3.2f, 3.2f);

        plot1d("h_min_dphijmet"+s, dphijmet, weight, hvec, ";min #Delta#phi(j, E_{T}^{miss}) ", 32,  0, 3.2);
        plot1d("h_dphi_boson_met"+s, dphiVmet, weight, hvec, ";#Delta#phi(boson, E_{T}^{miss}) ", 32,  0, 3.2);
        plot1d("h_dphi_lljets_met"+s, dphilljmet, weight, hvec, ";#Delta#phi(ll+j, E_{T}^{miss}) ", 32,  0, 3.2);
        plot1d("h_min_dphiVjet"+s, dphiVjet , weight, hvec, ";min #Delta#phi(j, boson)", 32,  0, 3.2);;
        if (event_Njets() >= 2) {
          plot1d("h_DjjVBF"+s, event_DjjVBF(), weight, hvec, ";D_{jj}(VBF) ", 52,  -0.02f, 1.02f);
        }

        if (event_Njets() == 0)
          plot1d("h_boson_pt_b0"+s, Vpt, weight, hvec, ";p_{T}(boson) [GeV]" , ptbin0.size()-1, ptbin0.data());
        plot1d("h_boson_pt_b1"+s, Vpt, weight, hvec, ";p_{T}(boson) [GeV]" , ptbin1.size()-1, ptbin1.data());
        
        const vector<float> mtbin1 = {0, 75, 150, 225, 300, 375, 450, 525, 600, 725, 850, 975, 1100};
        plot1d("h_mtZZ_b1"+s,  mtZZ , weight, hvec, ";M_{T}(ZZ) [GeV]" , mtbin1.size()-1, mtbin1.data());
        const vector<float> mtbin3 = {0, 75, 150, 225, 300, 375, 450, 525, 600, 725, 850, 975, 1100};
        plot1d("h_mtZZ_b3"+s,  mtZZ , weight, hvec, ";M_{T}(ZZ) [GeV]" , mtbin3.size()-1, mtbin3.data());

        plot1d("h_nphotons"+s, event_Nphotons(), weight, hvec, ";N(#gamma) ", 5,  0, 5);
        plot1d("h_ht"+s, event_HT() , weight, hvec, ";HT [GeV]" , 150,  0, 750);
        if (!is_data) {
          plot1d("h_genmet"+s, genmet_pt(), weight, hvec, ";gen-E_{T}^{miss} [GeV]"  , 160,  0, 800);
          plot1d("h_genmetphi"+s, genmet_phi(), weight, hvec, ";#phi(gen-E_{T}^{miss})"  , 64,  -3.2f, 3.2f);
        }
        if (isllg) {
          plot1d("h_ptll"+s, pt_boson(), weight, hvec, ";p_{T}^{ll} [GeV]"  , 160,  0, 800);
          plot1d("h_etall"+s, eta_boson(), weight, hvec, ";#eta(ll)"  , 64,  -3.2f, 3.2f);
          plot1d("h_mll"+s, mass_boson() , weight, hvec, ";M_{ll} [GeV]" , 125,  0, 500);
        } else {
          LorentzVector p4ll = LorentzVector(pt_l1(), eta_l1(), phi_l1(), 0) + LorentzVector(pt_l2(), eta_l2(), phi_l2(), 0);
          plot1d("h_ptll"+s, p4ll.pt(), weight, hvec, ";p_{T}^{ll} [GeV]"  , 160,  0, 800);
          plot1d("h_etall"+s, p4ll.eta(), weight, hvec, ";#eta(ll)"  , 64,  -3.2f, 3.2f);
          plot1d("h_mll"+s, p4ll.M() , weight, hvec, ";M_{ll} [GeV]" , 125,  0, 500);
        }

        if (event_Njets() > 0) {
          plot1d("h_jet1pt"+s, pt_jet1(), weight, hvec, ";p_{T}^{jet1} [GeV]" , 160,  0, 800);
          plot1d("h_jet1eta"+s, eta_jet1(), weight, hvec, ";#eta(jet1)"  , 64,  -3.2f, 3.2f);
        }
        if (event_Njets() > 1) {
          plot1d("h_jet2pt"+s, pt_jet2(), weight, hvec, ";p_{T}^{jet2} [GeV]" , 160,  0, 800);
          plot1d("h_jet2eta"+s, eta_jet2(), weight, hvec, ";#eta(jet2)"  , 64,  -3.2f, 3.2f);
        }

        if (pt_l1() > pt_l2()) {
          plot1d("h_lep1pt"+s, pt_l1(), weight, hvec, ";p_{T}^{lep1} [GeV]" , 160,  0, 800);
          plot1d("h_lep1eta"+s, eta_l1(), weight, hvec, ";#eta(lep1)"  , 64,  -3.2f, 3.2f);
          plot1d("h_lep2pt"+s, pt_l2(), weight, hvec, ";p_{T}^{lep1} [GeV]" , 160,  0, 800);
          plot1d("h_lep2eta"+s, eta_l2(), weight, hvec, ";#eta(lep1)"  , 64,  -3.2f, 3.2f);
        } else {
          plot1d("h_lep1pt"+s, pt_l2(), weight, hvec, ";p_{T}^{lep2} [GeV]" , 160,  0, 800);
          plot1d("h_lep1eta"+s, eta_l2(), weight, hvec, ";#eta(lep2)"  , 64,  -3.2f, 3.2f);
          plot1d("h_lep2pt"+s, pt_l1(), weight, hvec, ";p_{T}^{lep2} [GeV]" , 160,  0, 800);
          plot1d("h_lep2eta"+s, eta_l1(), weight, hvec, ";#eta(lep2)"  , 64,  -3.2f, 3.2f);
        }
        float minDR_Vlep(4.);
        isCloseObject(eta_l1(), phi_l1(), Veta, Vphi, minDR_Vlep, &minDR_Vlep);
        isCloseObject(eta_l2(), phi_l2(), Veta, Vphi, minDR_Vlep, &minDR_Vlep);
        plot1d("h_minDRVlep"+s, minDR_Vlep, weight, hvec, ";min#DeltaR(lep,boson)" , 64,  0.f, 3.2f);

        if (event_Nphotons() > 0 && skimver >= 2.08) {
          plot1d("h_photon_r9"+s, r9_photon(), weight, hvec, ";R9 (5x5) " , 70,  0, 1.4);
          plot1d("h_photon_pfiso"+s, pfiso_photon(), weight, hvec, ";pfIso" , 100,  0, 5);
          plot1d("h_photon_chiso"+s, chiso_photon(), weight, hvec, ";pfIso" , 50,  0, 2.5);
          plot1d("h_photon_nhiso"+s, pfiso_photon(), weight, hvec, ";pfIso" , 50,  0, 2.5);

          const vector<float> phptbin1 = {55, 82.5, 100, 135, 180, 220, 450};
          const vector<float> phptbin2 = {55, 82.5, 100, 135, 180, 190, 450};          

          int icat1 = std::upper_bound(phptbin1.begin(), phptbin1.end(), pt_photon()) - phptbin1.begin() - 1;
          int icat2 = std::upper_bound(phptbin2.begin(), phptbin2.end(), pt_photon()) - phptbin2.begin() - 1;
          plot1d(Form("h_R9_g_phptbin%d%s", icat1, s.c_str()), r9_photon(), weight, hvec, ";R9_{5x5}(#gamma)"  , 70,  0.f, 1.4f);
          plot1d(Form("h_R9_g_phpt2bin%d%s", icat2, s.c_str()), r9_photon(), weight, hvec, ";R9_{5x5}(#gamma)"  , 70,  0.f, 1.4f);
          plot1d(Form("h_pt_g_phptbin%d%s", icat1, s.c_str()), pt_photon(), weight, hvec, ";p_{T}^{#gamma} [GeV]" , 160,  0, 800);
          plot1d("hden_pt_b1_g_r9id90"+s,  pt_photon(), weight, hvec, ";p_{T}^{#gamma} [GeV]" , phptbin1.size()-1, phptbin1.data());
          plot1d("hden_pt_b2_g_r9id90"+s,  pt_photon(), weight, hvec, ";p_{T}^{#gamma} [GeV]" , phptbin2.size()-1, phptbin2.data());
          if (r9_photon() > 0.9) {
            plot1d("hnum_pt_b1_g_r9id90"+s,  pt_photon(), weight, hvec, ";p_{T}^{#gamma} [GeV]" , phptbin1.size()-1, phptbin1.data());
            plot1d("hnum_pt_b2_g_r9id90"+s,  pt_photon(), weight, hvec, ";p_{T}^{#gamma} [GeV]" , phptbin2.size()-1, phptbin2.data());
          }
        }
        if (event_Nphotons() > 0 && skimver >= 2.09) {
          plot1d("h_photon_phi"+s, phi_photon(), weight, hvec, ";#phi(#gamma)"     , 64,  -3.2, 3.2);
          plot1d("h_photon_R9"+s,  r9_photon(), weight, hvec, ";R9(#gamma)" , 80,  0.6, 1.4);
          plot1d("h_photon_sieie"+s,  sieie_photon(), weight, hvec, ";#sigma_{i#etai#eta}(#gamma)" , 80,  0, 0.02);
          plot1d("h_photon_sipip"+s,  sipip_photon(), weight, hvec, ";#sigma_{i#phii#phi}(#gamma)" , 80,  0, 0.02);
          plot1d("h_photon_seedtime"+s,  seedtime_photon(), weight, hvec, ";t_{seed}(#gamma) [ns]" , 120,  -15, 15);
          plot1d("h_photon_e4oe1"+s,  e4oe1_photon(), weight, hvec, ";E4/E1(#gamma)" , 80,  0, 2);
          plot1d("h_photon_mipE"+s,  mipE_photon(), weight, hvec, ";E_{MIP}(#gamma) [GeV]" , 60,  0, 12);
        }

        string phptsuf = "55to83";
        if (pt_boson() < 83) phptsuf = "55to83";
        else if (pt_boson() < 100) phptsuf = "83to100";
        else if (pt_boson() < 135) phptsuf = "100to135";
        else if (pt_boson() < 180) phptsuf = "135to180";
        else if (pt_boson() < 230) phptsuf = "180to230";
        else phptsuf = "230toInf";
        plot1d("h_mtZZ_phpt"+phptsuf+s, mtZZ, weight, hvec, ";m_{T}(ZZ) [GeV]"  , 160,  0, 800);

      };

      string metsuf = "_fullMET";
      fillhists(metsuf);
      fillhists(metsuf+jetcat);
      fillhists(metsuf+lepcat);
      fillhists(metsuf+jetcat+lepcat);

      if (met < 50) metsuf = "_metlt50";
      else if (met < 125) metsuf = "_met50to125";
      else metsuf = "_metge125";
      fillhists(metsuf);
      fillhists(metsuf+jetcat);
      fillhists(metsuf+lepcat);
      fillhists(metsuf+jetcat+lepcat);

      if (125 < pt_boson() && pt_boson() < 275) {
        fillhists("_Vpt125to275"+metsuf+jetcat+lepcat);
      }

      if (is_VBFcat()) 
        fillhists(metsuf+"_vbf"+lepcat);

      if (met < 125)
        fillhists("_metlt125"+jetcat+lepcat);
      if (met > 50)
        fillhists("_metge50"+jetcat+lepcat);

      if (met > 80)
        fillhists("_metge80"+jetcat+lepcat);
      else
        fillhists("_metlt80"+jetcat+lepcat);
      if (met > 80 && met < 125)
        fillhists("_met80to125"+jetcat+lepcat);

      if (isgamma) {
        if (event_Njets() == 0 && event_Njets20() == 1) {
          fillhists("_fullMET_1j20"+lepcat);
        }
      }

      if (produceResultTree) {
        evt_weight = weight;
        ll_pt = Vpt;
        ll_eta = Veta;
        ll_phi = Vphi;
        ll_mass = mboson;
        ptmiss = met;
        ptmiss_phi = metphi;
        mT = mtZZ;
        lepton_cat = is_ee()? 0 : is_mumu()? 1 : is_emu()? 2 : isgamma? 3 : -1;
        jet_cat = (is_VBFcat())? 2 : (event_Njets() == 0)? 0 : 1;

        // mindphi_jet_met = dphijmet; <-- already assigned
        // dphi_boson_met = dphiVmet; <-- already assigned
        // dphi_lljets_met = dphilljmet; <-- already assigned
        // mZZ = M_ZZ; <-- already assigned
        njet = event_Njets();
        nvtxs = event_nvtxs_good();
        if (skimver >= 2.10) {
          run = event_run();
          lumi = event_lumi();
          event = event_event();
        }

        if (isgamma) {
          // Reject photons that could come from electrons
          photon_convveto = convveto_photon();
          photon_id_cutbase = true;
          photon_isEB = true;
          // photon_isEB = ((photons_fid_mask()[i] & 1) == 1);
          // photon_id_cutbase = ((photons_id_cutBased_Fall17V2_Tight_Bits()[i] & 0x7f) == 0x7f);
          // photon_id_mva90 = photons_id_MVA_Fall17V2_pass_wp90()[i];
          // photon_id_mva80 = photons_id_MVA_Fall17V2_pass_wp80()[i];
          // photon_id_Hgg = photons_id_cutBased_HGG_Bits()[i];

          photon_pt = pt_photon();
          photon_eta = eta_photon();
          photon_phi = phi_photon();
          photon_r9 = r9_photon();
          photon_sieie = (skimver > 2.1)? sieie_photon() : 999;
          photon_sipip = (skimver > 2.1)? sipip_photon() : 999;
          photon_Emip = (skimver > 2.1)? mipE_photon() : -999;
          photon_seedtime = (skimver > 2.1)? seedtime_photon() : 0;
        } else {
          photon_convveto = false;
          photon_isEB = false;
          photon_id_cutbase = false;
          // photon_id_mva90 = false;
          // photon_id_mva80 = false;
          // photon_id_Hgg = false;
          photon_pt = -999.;
          photon_eta = -999.;
          photon_phi = -999.;
          photon_r9 = -999.;
          photon_sieie = -999.;
          photon_sipip = -999.;
          photon_Emip = -999.;
          photon_seedtime = -999.;
        }

        tout->Fill();
      }


      // // if (isgamma && met > 500 && sieie_photon() > 0.001) {
      // if (isgamma && is_data && skimver > 2.1 && met > 500 && sieie_photon() > 0.001 && event_Njets() == 0) {
      //   // cout << __LINE__ << ": met= " << met << ", pt_photon()= " << pt_photon() << ", pt_boson= " << pt_boson() << ", eta_boson()= " << eta_boson() << ", event_run= " << event_run() << ", event_lumi()= " << event_lumi() << ", event_event()= " << event_event() << ", event= " << event << ", currentFile->GetTitle()= " << currentFile->GetTitle() << endl;
      //   cout << __LINE__ << ": met= " << met << ", pt_photon()= " << pt_photon() << ", pt_boson= " << pt_boson() << ", eta_boson()= " << eta_boson() << ", event_run= " << event_run() << ", event_lumi()= " << event_lumi() << ", event_event()= " << event_event() << endl;
      //   if (!pass_trackveto())
      //     cout << __LINE__ << ": met= " << met << ", pt_photon()= " << pt_photon() << ", pt_boson= " << pt_boson() << ", eta_boson()= " << eta_boson() << ", event_run= " << event_run() << ", event_lumi()= " << event_lumi() << ", event_event()= " << event_event() << endl;
      // }
      // else if (isgamma && is_data && skimver < 2.1 && met > 500 && event_Njets() == 0) {
      //   cout << __LINE__ << ": met= " << met << ", pt_photon()= " << pt_photon() << ", pt_boson= " << pt_boson() << ", eta_boson()= " << eta_boson() << endl;
      //   if (!pass_trackveto())
      //     cout << __LINE__ << ": met= " << met << ", pt_photon()= " << pt_photon() << ", pt_boson= " << pt_boson() << ", eta_boson()= " << eta_boson() << endl;
      // }

      if (is_VBFcat())
        fillhists(metsuf+"_vbf"+lepcat);
      else if (event_Njets() >= 1)
        fillhists(metsuf+"_geq1j"+lepcat);

    }//event loop

    delete file;
  }//file loop

  for (auto& h : hvec) {
    // if (h.first.find("hnum") == 0 || h.first.find("hden") == 0) continue;
    if (h.first.find("phi") == string::npos && h.first.find("2d") == string::npos)
      moveOverFlowToLastBin1D(h.second);
    string dirname = "OffShell";
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
    if (dir == nullptr) dir = fout->mkdir(dirname.c_str());
    dir->cd();
    h.second->Write();
  }

  if (produceResultTree) {
    fout->cd();
    tout->Write();
  }
  fout->Close();

  return 0;
}

