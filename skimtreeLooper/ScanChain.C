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
const bool blind = true;

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

  if (tag.Contains("rwgtd")) {
    doBosonPtReweight = true;
    cout << ">> doBosonPtReweight is set to true." << endl;
  }
  extendEEphoton2j = tag.Contains("ee2j");

  bool is_data = (sample.BeginsWith("201"));
  bool is_gjets = sample.BeginsWith("GJets");

  int year = (tag.Contains("2016"))? 2016 : (tag.Contains("2017"))? 2017 : 2018;

  map<string,TH1*> hvec;

  unsigned nEventsTotal = 0;
  unsigned nEventsChain = ch->GetEntries();

  cout << ">> To run over nEventsChain = " << nEventsChain << endl;

  TFile *currentFile = 0;
  TObjArray *listOfFiles = ch->GetListOfFiles();
  TIter fileIter(listOfFiles);

  while ( (currentFile = (TFile*)fileIter.Next()) ) { 

    TFile *file = new TFile( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("SkimTree");
    st.Init(tree);

    TString filename(currentFile->GetTitle());

    for( unsigned int event = 0; event < tree->GetEntriesFast(); ++event) {

      st.GetEntry(event);
      nEventsTotal++;

      SkimTree::progress(nEventsTotal, nEventsChain);

      double weight = event_wgt();
      if (is_data) weight = 1;

      bool isllg = (is_mumu() || is_ee()) && (event_Nphotons() == 1);
      bool isgamma = is_gamma();
      // bool isgamma = event_Nphotons() == 1 && !is_mumu() && !is_ee();
      if (extendEEphoton2j) isgamma = is_gamma() || (event_Njets() >= 2 && event_Nphotons() == 1 && !is_mumu() && !is_ee());

      double wgt = weight;
      auto fillPhotonHitMap = [&](string suf) {
        plot2d("h2d_photon_eta_phi"+suf, eta_photon(), phi_photon(), wgt, hvec, ";#eta(#gamma);#phi(#gamma)" , 50, -2.5f, 2.5f, 64, -3.2f, 3.2f);
        plot1d("h_photon_eta"+suf, eta_photon(), wgt, hvec, ";#eta(#gamma)" , 50, -2.5f, 2.5f);
        plot1d("h_photon_pt"+suf, pt_photon(), wgt, hvec, ";p_{T}($gamma) [GeV]" , 120,  0, 600);
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
        if (mindphi_jet_met() < 0.25) continue;
        if (dphi_boson_met() < 1.0) continue;
        if (dphi_lljets_met() < 2.5) continue;
        if (pt_boson() < 55.) continue;
        if (fabs(mass_boson() - 91.2) > 15.) continue;
        // if (isgamma && !convveto_photon()) continue;  // TODO: to confirm the effect of this cut

        if (blind && is_data && (is_ee() || is_mumu()) && pTmiss() >= 125) continue;
      }

      float met = pTmiss();
      float metphi = phimiss();
      float M_ZZ = mZZ();
      float mtZZ = mTZZ();
      float mboson = mass_boson();
      float Vpt  = pt_boson();
      float Veta = eta_boson();
      float Vphi = phi_boson();
      float dphijmet = mindphi_jet_met();
      float dphiVmet = dphi_boson_met();
      float dphilljmet = dphi_lljets_met();
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

        if (dphilljmet < 2.5) continue;
        if (dphijmet < 0.25) continue;

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
      if (has_ak4jets_inHEM1516()) continue;
      if (event_Njets_btagged() > 0) continue;

      string jetcat = "_jetcat";
      if (event_Njets() == 0) jetcat = "_eq0j";
      else if (event_Njets() == 1) jetcat = "_eq1j";
      else if (event_Njets() == 2) jetcat = "_eq2j";
      else if (event_Njets() >= 3) jetcat = "_ge3j";

      if (isgamma) {
        // if (fabs(eta_boson()) > 1.4442) continue;  // TODO: to remove this line later
        if (is_gjets) {
          weight *= std::max(1., 1.716910-0.001221*pt_boson());
        }
        int icat = (pt_boson() > 440)? 39 : (pt_boson() - 50) / 10;
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
            if (is_data) {
              if      (jetcat == "_eq1j") scale = sf_Vpt_data_metlt80_eq1j.at(icat);
              else if (jetcat == "_eq0j") scale = sf_Vpt_data_metlt80_eq0j.at(icat);
              else if (jetcat == "_eq2j") scale = sf_Vpt_data_metlt80_eq2j.at(icat);
            } else {
              if      (jetcat == "_eq1j") scale = sf_Vpt_closure_metlt80_eq1j.at(icat);
              else if (jetcat == "_eq0j") scale = sf_Vpt_closure_metlt80_eq0j.at(icat);
              else if (jetcat == "_eq2j") scale = sf_Vpt_closure_metlt80_eq2j.at(icat);
            }
          }
          weight *= scale;
        }
      } else {
      }

      // if (skimver >= 1.11 && event_Nphotons() == 1) {
      //   wgt = weight;
      //   fillPhotonHitMap("_fullMET");
      //   if (convveto_photon())
      //     fillPhotonHitMap("_fullMET_convveto");
      //   if (id_photon_Hgg())
      //     fillPhotonHitMap("_fullMET_Hggid");
      //   if (convveto_photon() && id_photon_Hgg())
      //     fillPhotonHitMap("_fullMET_convAndHggId");
      // }

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
        plot1d("h_boson_pt"+s, Vpt, weight, hvec, ";p_{T}(boson) [GeV]" , 120,  0, 600);
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
        const vector<float> mtbin1 = {0, 75, 150, 225, 300, 375, 450, 525, 600, 725, 850, 975, 1100, 1350, 1600, 2100, 3000};
        plot1d("h_mtZZ_b1"+s,  mtZZ , weight, hvec, ";M_{T}(ZZ) [GeV]" , mtbin1.size()-1, mtbin1.data());
        const vector<float> mtbin3 = {0, 75, 150, 225, 300, 375, 450, 525, 600, 725, 850, 975, 1100, 1350, 1600, 2100, 3000};
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

      if (is_VBFcat())
        fillhists(metsuf+"_vbf"+lepcat);
      else if (event_Njets() >= 1)
        fillhists(metsuf+"_geq1j"+lepcat);

    }//event loop

    delete file;
  }//file loop

  for (auto& h : hvec) {
    if (h.first.find("hnum") == 0 || h.first.find("hden") == 0) continue;
    if (h.first.find("phi") == string::npos && h.first.find("2d") == string::npos)
      moveOverFlowToLastBin1D(h.second);
    string dirname = "OffShell";
    vector<string> jetsufs = {"_eq0j", "_eq1j", "_eq2j", "_1j20", "_vbf"};
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

  fout->Close();

  return 0;
}

