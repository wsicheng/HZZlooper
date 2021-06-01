#pragma GCC diagnostic ignored "-Wsign-compare"

#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1.h"
#include "TChain.h"
#include "TSystem.h"
#include "TSpline.h"

#include "TriLepTree.h"
#include "Utilities.h"

using namespace std;
using namespace tas;

bool doBosonPtReweight = false;
bool extendEEphoton2j = false;
bool tightDphiIn2j = false;
bool tightMETin2j = false;
bool doBosonEtaReweight = false;
bool doNvtxReweight = false;
bool doHllgTest = false;
bool DYclosureTest = false;
bool applyZGweightFromLLG = true;
bool applyExternalAsMZ = false;
bool applyExternalPythiaScale = false;
bool fillAdjustWeightHists = false;
bool produceResultPlots = true;
bool produceInvtdMETtree = false;

const bool makeStepPlots = true;
const bool produceResultTree = true;
const bool separateFileForTree = true;
const bool blind = false;

int ScanChain(TChain* ch, TString sample, TString tag, TString systype = "", TString specifiers = "", TString extrargs = "") {

  TString samplever = tag(0, min(tag.Index("_",tag.Index("_")+1),5));
  double skimver = TString(samplever).ReplaceAll("v","").ReplaceAll("_",".").Atof();
  cout << ">> Runing from skimver " << samplever  << ", i.e. " << skimver << endl;

  bool is_data = (sample.BeginsWith("Run201")) || sample.BeginsWith("photon") || sample.BeginsWith("EGamma");

  int year = (tag.Contains("2016"))? 2016 : (tag.Contains("2017"))? 2017 : (tag.Contains("2018"))? 2018 : -1;
  if (year < 0) cout << ">> WARNING: not able to determine year from tag: " << tag << endl;

  vector<string> systnames = {"PDFScale", "QCDScale", "AsMZ", "PDFReplica", "PU", "BtagSF", "L1Prefiring", "PUJetIdEff", "PythiaScale",
    "ElectronEff_AltMC", "ElectronEff_Syst", "ElectronEff_Stat", "MuonEff_AltMC", "MuonEff_Syst", "MuonEff_Stat"};

  vector<string> systematics = {"Nominal"};
  if (!is_data) {
    for (string isys : systnames) {
      systematics.push_back(isys+"Up");
      systematics.push_back(isys+"Dn");
    }
  }

  if (systype.BeginsWith("NominalOnly")) {
    systematics = vector<string>{"Nominal"};
    systype = "Nominal";
  }
  // Make systype to the outtree name
  if (!systype.BeginsWith("Nominal")) {
    produceResultPlots = false;
    systematics = vector<string>{systype.Data()};
    if (is_data) return 1;
  }

  TFile* fout = nullptr;
  TString outname = sample;

  if (produceResultPlots) {
    TString fileout = Form("output/%s/%s.root", tag.Data(), outname.Data());
    gSystem->Exec(Form("mkdir -p output/%s", tag.Data()));
    fout = new TFile(fileout, "RECREATE");
    cout << ">> Outputting to " << fileout << endl;
  }

  TFile* f_trigSF = new TFile(Form("/usskimLooper/data/TriggerSF/%d/trigger_efficiencies_leptons.root", year));

  map<string,TH1*> hvec;

  unsigned nEventsTotal = 0;
  unsigned nEventsChain = ch->GetEntries();

  // Result tree producer
  map<string, TFile*> fouts;
  map<string, TTree*> touts;
  float dphijmet, dphiVmet, dphilljmet, dphi3ljmet, dphil3met;
  float trigwgt(1.);

  // Variable only store for the output tree
  float evt_weight;
  // float ll_pt, ll_eta, ll_phi, ll_mass;
  float met(0.), metphi(-99.), rlmet(0.), rlmetphi(-99.);
  float mtl3(0.), m3l(0.), mllmin(9999.), mll13(0.), mll23(0.), mllclose(0.), mllsfx3(0.);
  float mZZ(0.), mWZ(0.), mtZZ(0.), mtWZ(0.), rlmtZZ(0.);
  float Vpt(0.), Veta(-99.), Vphi(-99.), Vmass(0.);
  float lep3pt(0.), lep3eta(-99.), lep3phi(-99.), lep3mass(0.);
  int llid(0), lep3id(0);

  float jet_pt[8], jet_eta[8], jet_phi[8], jet_mass[8];

  unsigned njet(0);

  bool is_vbf(false), isSF13(false), isSF23(false);

  map<string,float> evtwgt_sys;

  map<string,TH2D*> hmap_trigSF_dilep;
  map<string,TH2D*> hmap_trigEff_sel;
  map<string,TH1*> hmap_externwgt;

  if (produceResultTree) {
    for (string isys : systematics) {
      if (separateFileForTree) {
        gSystem->Exec(Form("mkdir -p outtree/%s", tag.Data()));
        TString otname = Form("outtree/%s/trileptree_%s_%s.root", tag.Data(), outname.Data(), isys.c_str());
        cout << ">> Outputting result tree to " << otname << endl;
        fouts[isys] = new TFile(otname, "RECREATE");
        fouts[isys]->cd();
      }
      bool extraBranches = (isys == "Nominal");

      TTree* tout = new TTree("FinalTree", "");

      tout->Branch("weight", &evt_weight);
      tout->Branch("mTZZ", &mtZZ);
      tout->Branch("mTWZ", &mtWZ);
      tout->Branch("mWZ", &mWZ);
      tout->Branch("mTl", &mtl3);

      tout->Branch("dilepton_id", &llid);
      tout->Branch("dilepton_pt", &Vpt);
      tout->Branch("dilepton_eta", &Veta);
      tout->Branch("dilepton_mass", &Vmass);

      tout->Branch("lepW_id", &lep3id);
      tout->Branch("lepW_pt", &lep3pt);
      tout->Branch("lepW_eta", &lep3eta);
      tout->Branch("lepW_phi", &lep3phi);
      tout->Branch("lepW_mass", &lep3mass);

      tout->Branch("pTmiss", &met);
      tout->Branch("phimiss", &metphi);
      tout->Branch("pTmiss_addl3", &rlmet);
      tout->Branch("phimiss_addl3", &rlmetphi);

      tout->Branch("n_ak4jets_pt30", &njet);

      tout->Branch("ak4jet_leading_pt", jet_pt);
      tout->Branch("ak4jet_leading_eta", jet_eta);
      tout->Branch("ak4jet_leading_phi", jet_phi);
      tout->Branch("ak4jet_leading_mass", jet_mass);

      tout->Branch("ak4jet_subleading_pt", jet_pt+1);
      tout->Branch("ak4jet_subleading_eta", jet_eta+1);
      tout->Branch("ak4jet_subleading_phi", jet_phi+1);
      tout->Branch("ak4jet_subleading_mass", jet_mass+1);

      if (extraBranches) {
        tout->Branch("m3l", &m3l);
        tout->Branch("mllmin", &mllmin);
      }

      touts[isys] = tout;
    }
  }

  // Make photon-pt bins for reweight
  vector<float> ptbin0, ptbin1;
  std::tie(ptbin0, ptbin1) = getPtBins();

  cout << ">> To run over nEventsChain = " << nEventsChain << endl;

  TFile *currentFile = 0;
  TObjArray *listOfFiles = ch->GetListOfFiles();
  TIter fileIter(listOfFiles);

  while ( (currentFile = (TFile*)fileIter.Next()) ) {
    TString filename(currentFile->GetTitle());
    // Skip other systematics other than specified
    if (systype != "" && !filename.Contains("_"+systype))
      continue;
    // Skip files that doesn't contain any one of the specifed string
    TString fullstr = specifiers;
    bool skip = false;
    while (fullstr.Contains(',')) {
      TString substr = fullstr(0, fullstr.Index(','));
      if (!filename.Contains(substr)) skip = true;
      fullstr.Remove(0, substr.Length()+1);
    }
    if (skip || (fullstr != "" && !filename.Contains(fullstr))) {
      cout << "Skipping: " << filename << endl;
      continue;
    }

    bool isNLO = (filename.Contains("amcatnlo") || filename.Contains("powheg"));
    bool isWG_nlo = (filename.Contains("WGToLNuG") && filename.Contains("amcatnlo"));
    bool isZG_nlo_incl = (filename.Contains("ZGTo2NuG_Tune") || filename.Contains("ZGTo2LG_Tune"));
    bool isZG_nlo_ptG130 = (filename.Contains("ZGTo2NuG_PtG-130") || filename.Contains("ZGTo2LG_PtG-130"));
    bool isZllG_lo = filename.Contains("ZLLGJets"); // for sync with nlo
    bool isZllG_lo_ptG130 = (isZllG_lo && filename.Contains("PtG-130")); // for sync with nlo
    bool isWZ3l_mllmin01 = (filename.Contains("WZTo3LNu_mllmin01"));
    bool useNNPDF30 = isNLO && !(isZG_nlo_incl || isZG_nlo_ptG130);
    bool applyKfactor = filename.Contains("ZZTo") || filename.Contains("WZTo");

    TFile *file = new TFile(filename);
    TTree *tree = (TTree*)file->Get("SkimTree");
    st.Init(tree);

    bool isdilep(false), isgamma(false), isllg(false), islg(false), is1el(false), is3l;
    bool fnloaded = false;
    string lepcat = "_lepcat";
    double event_wgt_sample = 1.0;
    double event_wgt_adjust = 1.0;

    for (unsigned int evt = 0; evt < tree->GetEntriesFast(); ++evt) {

      st.GetEntry(evt);
      nEventsTotal++;

      TriLepTree::progress(nEventsTotal, nEventsChain);

      double weight = event_wgt();
      if (is_data) weight = 1;
      if (weight == 0.) continue;

      if (!fnloaded) {
        if (filename.Contains("/3LEvents/")) {
          is3l = true;
          hmap_trigSF_dilep = loadScaleFactorHist("dilep", f_trigSF);
        }
        // event_wgt_sample = getExtSampleWeight(filename, year);
        event_wgt_sample = 1.0;
        fnloaded = true;
      }

      if (is3l) {
        switch (dilepton_id()) {
          case -121:
            lepcat = "_eel";
            break;
          case -169:
            lepcat = "_mumul";
            break;
          case -143:
          case 143:
            lepcat = "_emul";
            break;
          case 121:
          case 169:
            lepcat = "_sssf";
            break;
          default:
            cout << "WARNING: Invalid dilepton ID: " << dilepton_id() << "!!" << endl;
        }

        // Simplified trigger verification: any firing one is taken
        for (auto trwgt : event_wgt_triggers_SingleLepton()) {
          if (trwgt > 0) trigwgt = 1.0;
        }
        for (auto trwgt : event_wgt_triggers_Dilepton_SF()) {
          if (trwgt > 0) trigwgt = 1.0;
        }
      }

      if (trigwgt == 0.f) continue;

      double event_wgt_SFs = 1.0;
      double wgt_EWup(1.0), wgt_EWdn(1.0);
      if (!is_data) {
        // weight *= event_wgt_pileup(); // pileup weight included in event_wgt already
        event_wgt_SFs = event_wgt_SFs_electrons() * event_wgt_SFs_muons() * event_wgt_SFs_photons();
        event_wgt_SFs *= event_wgt_SFs_PUJetId() * event_wgt_SFs_btagging();
        event_wgt_SFs = std::min(event_wgt_SFs, 5.0);
        weight *= event_wgt_SFs;
        weight *= event_wgt_sample;

        event_wgt_adjust = (useNNPDF30)? event_wgt_adjustment_NNPDF30() : 1;
        if (applyKfactor) {
          weight *= KFactor_EW_NLO_qqVV_Bkg_Nominal() * KFactor_QCD_NNLO_qqVV_Bkg_Nominal();
          wgt_EWup = KFactor_EW_NLO_qqVV_Bkg_EWUp() / KFactor_EW_NLO_qqVV_Bkg_Nominal();
          wgt_EWdn = KFactor_EW_NLO_qqVV_Bkg_EWDn() / KFactor_EW_NLO_qqVV_Bkg_Nominal();
        }
      }

      // Check sample names
      // Fix the cross section for the mllmin01 samples
      if (isWZ3l_mllmin01) weight *= 0.6;

      int istep = 0;
      auto fill_passedsteps = [&](string s="") {
        if (!makeStepPlots) return;
        plot1d("h_metphi_step"+to_string(istep)+s, event_phimiss(), weight, hvec, ";#phi(E_{T}^{miss})"  , 64, -3.2, 3.2);
        plot1d("h_met_step"+to_string(istep)+s, event_pTmiss(), weight, hvec, ";E_{T}^{miss} [GeV]"  , 160, 0, 800);
        // plot1d("h_passed_steps"+lepcat, istep , weight, hvec, ";step" , 20,  0, 20);
        plot1d("h_passed_steps", istep , weight, hvec, ";step" , 20,  0, 20);
        istep++;
      };
      fill_passedsteps("_nocut");

      // btag veto, already applied on SinglePhotonEvents but not others
      if (event_n_ak4jets_pt30_btagged_loose() > 0) continue;

      // extra lepton veto
      if (isgamma && event_n_leptons_fakeableBase() > 0) continue;
      if ((is1el || islg) && event_n_leptons_fakeableBase() > 1) continue;
      if ((isllg || isdilep) && event_n_leptons_fakeableBase() > 2) continue;
      if (is3l && event_n_leptons_fakeableBase() > 0) continue;
      fill_passedsteps("_ljveto");

      met = event_pTmiss();
      metphi = event_phimiss();
      mZZ = -999;
      mWZ = -999;
      mtZZ = -999;
      mtWZ = -999;
      rlmet = -999;
      rlmetphi = -999;
      rlmtZZ = -999;
      Vpt   = dilepton_pt();
      Veta  = dilepton_eta();
      Vphi  = dilepton_phi();
      Vmass = dilepton_mass();

      mtl3 = event_mTl();
      m3l = event_m3l();
      mllmin = 9999.;
      lep3id = -1;

      int ilep1(-1), ilep2(-1), ilep3(-1);
      int itrig = -1;

      if (is3l) {
        if (leptons_id().size() != 3) continue;

        LorentzVector p4dilep(dilepton_pt(), dilepton_eta(), dilepton_phi(), dilepton_mass());

        if (dilepton_id() == -121 || dilepton_id() == -169) {
          // case that the dilepton branch get it right?
          ilep1 = dilepton_daughter_indices().at(0);
          ilep2 = dilepton_daughter_indices().at(1);
          if (leptons_pt().at(ilep1) < leptons_pt().at(ilep2)) std::swap(ilep1, ilep2);
          ilep3 = 3 - ilep1 - ilep2;
          itrig = ilep1 + ilep2 - 1;

          string s = "_cat0";
          plot1d("h_mdilep"+s, dilepton_mass(), weight, hvec, ";M(ll) [GeV]" , 160, 0, 800);
          plot1d("h_mTW"+s, event_mTl(), weight, hvec, ";M_{T}(W) [GeV]" , 160, 0, 800);
          plot1d("h_mWvis"+s, event_mWVis(), weight, hvec, ";M(W vis.) [GeV]" , 160, 0, 800);

          plot1d("h_trilep_cat", 0, weight, hvec, "; trilep category" , 8, 0, 8);

          lep3id = leptons_id().at(ilep3);

          if (dilepton_id() == -169) {
            plot1d("h_trilep_cat", 1, weight, hvec, "; trilep category" , 8, 0, 8);
            if (fabs(lep3id) == 13) {
              lepcat = "_mumumu";
              plot1d("h_trilep_cat", 3, weight, hvec, "; trilep category" , 8, 0, 8);
            } else {
              lepcat = "_mumue";
            }
          }
          if (dilepton_id() == -121) {
            plot1d("h_trilep_cat", 2, weight, hvec, "; trilep category" , 8, 0, 8);
            if (fabs(lep3id) == 13) {
              lepcat = "_eemu";
              plot1d("h_trilep_cat", 4, weight, hvec, "; trilep category" , 8, 0, 8);
            } else {
              lepcat = "_eee";
            }
          }

        } else {
          // Find the Z by ossf
          int nOSSFpair = 0;
          float mdilep = -999.;
          LorentzVector p4lep1, p4lep2;
          const vector<int> ilep1s = {1, 0, 0};
          const vector<int> ilep2s = {2, 2, 1};
          vector<int> idprod2 = {leptons_id()[1]*leptons_id()[2], leptons_id()[0]*leptons_id()[2], leptons_id()[0]*leptons_id()[1]};
          for (int i = 0; i < 3; ++i) {
            auto p4lep1_tmp = LorentzVector(leptons_pt().at(ilep1s[i]), leptons_eta().at(ilep1s[i]), leptons_phi().at(ilep1s[i]), leptons_mass().at(ilep1s[i]));
            auto p4lep2_tmp = LorentzVector(leptons_pt().at(ilep2s[i]), leptons_eta().at(ilep2s[i]), leptons_phi().at(ilep2s[i]), leptons_mass().at(ilep2s[i]));
            auto p4dilep_tmp = p4lep1_tmp + p4lep2_tmp;
            auto mdilep_tmp = p4dilep_tmp.M();

            if (mdilep_tmp < mllmin)
              mllmin = mdilep_tmp;

            if (idprod2[i] == -121 || idprod2[i] == -169) {
              if (fabs(mdilep_tmp - 91.2) < fabs(mdilep - 91.2)) {
                std::swap(p4lep1, p4lep1_tmp);
                std::swap(p4lep2, p4lep2_tmp);
                std::swap(p4dilep, p4dilep_tmp);
                std::swap(mdilep, mdilep_tmp);
                ilep3 = i;
              }
              nOSSFpair++;
            }
          }
          if (nOSSFpair == 0) {
            plot1d("h_trilep_cat", 7, weight, hvec, "; trilep category" , 8, 0, 8);
            continue;
          }
          plot1d("h_mll_cat6", mdilep, weight, hvec, ";M(ll) [GeV]" , 160, 0, 800);
          plot1d("h_trilep_cat", 6, weight, hvec, "; trilep category" , 8, 0, 8);
          continue;
        }

        lep3pt = leptons_pt().at(ilep3);
        // float mtlep3 = calculateMT(lep3pt, leptons_phi().at(ilep3), event_pTmiss(), event_phimiss());
        LorentzVector p4lep1(leptons_pt().at(ilep1), leptons_eta().at(ilep1), leptons_phi().at(ilep1), leptons_mass().at(ilep1));
        LorentzVector p4lep2(leptons_pt().at(ilep2), leptons_eta().at(ilep2), leptons_phi().at(ilep2), leptons_mass().at(ilep2));
        LorentzVector p4lep3(leptons_pt().at(ilep3), leptons_eta().at(ilep3), leptons_phi().at(ilep3), leptons_mass().at(ilep3));
        mll13 = (p4lep1 + p4lep3).M();
        mll23 = (p4lep2 + p4lep3).M();
        mllmin = std::min(dilepton_mass(), std::min(mll13, mll23));
        mllclose = (fabs(mll13 - 91.2) < fabs(mll23 - 91.2))? mll13 : mll23;
        isSF13 = abs(leptons_id().at(ilep1)) == abs(leptons_id().at(ilep3));
        isSF23 = abs(leptons_id().at(ilep2)) == abs(leptons_id().at(ilep3));
        mllsfx3 = (isSF13 && isSF23)? mllclose : (isSF13)? mll13 : (isSF23)? mll23 : -1;

        // Calculate variables with 3rd lepton added to MET
        LorentzVector met_p4(event_pTmiss(), 0, event_phimiss(), 0);
        LorentzVector rlmet_p4 = met_p4 + p4lep3;

        rlmet = rlmet_p4.pt();
        rlmetphi = rlmet_p4.phi();
        mtZZ = getDileptonMT(p4dilep, met_p4);
        mtWZ = getDileptonMT(p4dilep, rlmet_p4, 80.4);
        rlmtZZ = getDileptonMT(p4dilep, rlmet_p4);
        LorentzVector p4nunu_approx(met_p4.pt(), p4dilep.eta(), met_p4.phi(), 91.2);
        mZZ = (p4nunu_approx + p4dilep).M();
        mWZ = (rlmet_p4 + p4dilep).M();
        // assuming that the dilepton branches and the p4dilep is consistent

        LorentzVector p4lljets = p4dilep;
        dphiVmet = deltaPhi(p4dilep.phi(), metphi);
        dphil3met = deltaPhi(p4lep3.phi(), metphi);
        dphijmet = 999.;
        for (int j = 0; j < event_n_ak4jets_pt30(); ++j) {
          p4lljets += LorentzVector(ak4jets_pt()[j], ak4jets_eta()[j], ak4jets_phi()[j], ak4jets_mass()[j]);
          float dphi_met = deltaPhi(ak4jets_phi()[j], metphi);
          if (dphi_met < dphijmet) dphijmet = dphi_met;
        }
        dphilljmet = deltaPhi(p4lljets.phi(), metphi);;
        dphi3ljmet = deltaPhi((p4lljets+p4lep3).phi(), metphi);;
      }
      bool passDilepTrig_Zll = event_wgt_triggers_Dilepton_SF().at(itrig);
      bool passSinglepTrig_Zll = event_wgt_triggers_SingleLepton().at(ilep1) || event_wgt_triggers_SingleLepton().at(ilep2);
      bool passSinglepTrig_lep3 = event_wgt_triggers_SingleLepton().at(ilep3);
      bool passDilepTrig_any = false;
      bool passSinglepTrig_any = false;
      for (bool triggered : event_wgt_triggers_Dilepton_SF()) {
        passDilepTrig_any |= triggered;
      }
      for (bool triggered : event_wgt_triggers_SingleLepton()) {
        passSinglepTrig_any |= triggered;
      }
      bool passAnyTrig = passDilepTrig_any || passSinglepTrig_any;

      auto fillmasshists = [&](string s = "") {
        plot1d("h_met"+s, event_pTmiss(), weight, hvec, ";E_{T}^{miss} [GeV]"  , 160,  0, 800);
        plot1d("h_rlmet"+s, met, weight, hvec, ";E_{T}^{miss} (rl) [GeV]"  , 160,  0, 800);
        plot1d("h_mll"+s, Vmass, weight, hvec, ";M_{ll} [GeV]" , 160, 0, 800);
        plot1d("h_m3l"+s, m3l, weight, hvec, ";M_{lll} [GeV]" , 160, 0, 800);
        plot1d("h_mtl3"+s, mtl3, weight, hvec, ";M_{T} (lep3) [GeV]" , 50, 0, 250);
        plot1d("h_mllmin"+s, mllmin, weight, hvec, ";min M_{ll} [GeV]" , 100, 0, 250);
      };
      fillmasshists("_nocut");

      // Lepton pT cuts
      float lep1pt(leptons_pt().at(ilep1)), lep2pt(leptons_pt().at(ilep2));
      if (is3l && (lep1pt < 30 || lep2pt < 20 || lep3pt < 20)) continue;
      fill_passedsteps("_leppt");

      fillmasshists("_leppt");

      if (met < 20) continue;  // explicit one to be orthogonal to 4l 
      // if ((1.4*mtl3 + met) < 120) continue;
      if ((fabs(lep3id) == 13) && mtl3 < 20) continue;
      if ((fabs(lep3id) == 11) && mtl3 < 10) continue;
      if ((fabs(lep3id) == 13) && (1.6*mtl3 + met) < 120) continue;
      if ((fabs(lep3id) == 11) && (4*mtl3/3 + met) < 120) continue;
      fill_passedsteps("_mtmet");

      fillmasshists("_mtmet");

      auto fillTrigEffPlots = [&](string s, string name, bool passTrig, bool passDen=true) {
        if (!passDen) return;
        plot1d("hden_"+name+"_lepZ1pt"+s, ((lep1pt > lep2pt)? lep1pt : lep2pt)  , weight, hvec, ";p_{T}(lepZ1) [GeV]"  , 25,  0, 500);
        plot1d("hden_"+name+"_lepZ2pt"+s, ((lep1pt > lep2pt)? lep2pt : lep1pt)  , weight, hvec, ";p_{T}(lepZ2) [GeV]"  , 20,  0, 400);
        plot1d("hden_"+name+"_lepWpt"+s,   lep3pt  , weight, hvec, ";p_{T}(lepW) [GeV]"  , 20,  0, 400);
        plot1d("hden_"+name+"_ptll"+s, dilepton_pt(), weight, hvec, ";p_{T}^{ll} [GeV]"  , 160,  0, 800);
        if (passTrig) {
          plot1d("hnum_"+name+"_lepZ1pt"+s, ((lep1pt > lep2pt)? lep1pt : lep2pt)  , weight, hvec, ";p_{T}(lepZ1) [GeV]"  , 25,  0, 500);
          plot1d("hnum_"+name+"_lepZ2pt"+s, ((lep1pt > lep2pt)? lep2pt : lep1pt)  , weight, hvec, ";p_{T}(lepZ2) [GeV]"  , 20,  0, 400);
          plot1d("hnum_"+name+"_lepWpt"+s,   lep3pt  , weight, hvec, ";p_{T}(lepW) [GeV]"  , 20,  0, 400);
          plot1d("hnum_"+name+"_ptll"+s, dilepton_pt(), weight, hvec, ";p_{T}^{ll} [GeV]"  , 160,  0, 800);
        }
      };
      fillTrigEffPlots("final", "trigeff_dilepZll", passDilepTrig_Zll, passAnyTrig);
      fillTrigEffPlots("final", "trigeff_combZll", (passDilepTrig_Zll || passSinglepTrig_Zll), passAnyTrig);

      if (!passDilepTrig_Zll && !passSinglepTrig_Zll) continue;
      fill_passedsteps("_trig");
      fillmasshists("_trig");

      njet = event_n_ak4jets_pt30();
      float dphiVjet = 4.0;

      string jetcat = "jetcat";
      if (njet == 0) jetcat = "_eq0j";
      else if (njet == 1) jetcat = "_eq1j";
      else if (njet >= 2) jetcat = "_ge2j";

      // if (Vpt < 55.0) continue;  // no boson pt requirement
      if (fabs(Vmass - 91.2) > 15) continue;  // only for dilepton events
      fill_passedsteps("_VptVmass");

      if (dphiVmet < 1.0) continue;
      if (fabs(dPhi_pTleptonsjets_pTmiss()) < 2.5) continue;
      if (njet > 0 && min_abs_dPhi_pTj_pTmiss() < 0.25) continue;
      fill_passedsteps("_dphis");

      // Record jet veto at horn regions
      bool jetveto = false;
      if (year > 2016) {
        for (int j = 0; j < event_n_ak4jets_pt30(); ++j) {
          if (2.65 < fabs(ak4jets_eta()[j]) && fabs(ak4jets_eta()[j]) < 3.139) {
            float nemf_thr = (year == 2018)? 0.7 : (year == 2017)? 0.5 : 1.;
            if (ak4jets_NEMF()[j] > nemf_thr) {
              jetveto = true;
            }
          }
        }
      }

      // End of cuts (except MET)
      // --------------------------
      // Begin of scale factors

      if (njet > 0) dphiVjet = deltaPhi(ak4jets_phi()[0], Vphi);
      if (njet > 1) dphiVjet = min(dphiVjet, deltaPhi(ak4jets_phi()[1], Vphi));

      // End of scale factors
      // --------------------------
      // Begin of final plots

      const vector<double> mtbins = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 1000, 1200};

      auto fillhists = [&](string s = "", bool coreonly = false) {
        if (!produceResultPlots) return;
        // core hists
        plot1d("h_met"+s, met, weight, hvec, ";E_{T}^{miss} [GeV]"  , 160,  0, 800);
        plot1d("h_rlmet"+s, rlmet, weight, hvec, ";E_{T}^{miss} (rl) [GeV]"  , 160,  0, 800);
        plot1d("h_boson_pt"+s, Vpt, weight, hvec, ";p_{T}(boson) [GeV]" , 160, 0, 800);
        plot1d("h_njets"+s, njet, weight, hvec, ";N(jets)"  , 6,  0, 6);
        plot1d("h_mtZZ"+s, mtZZ, weight, hvec, ";m_{T}^{ZZ} [GeV]"  , 160,  0, 800);
        plot1d("h_rlmtZZ"+s, rlmtZZ, weight, hvec, ";m_{T}^{ZZ} (rl) [GeV]"  , 160,  0, 800);
        plot1d("h_mtWZ"+s, mtWZ, weight, hvec, ";m_{T}^{WZ} [GeV]"  , 160,  0, 1200);
        plot1d("h_mWZ"+s, mWZ, weight, hvec, ";m^{WZ} [GeV]"  , 160,  0, 1200);
        plot1d("h_mtWZbins"+s, mtWZ, weight, hvec, ";m_{T}^{WZ} [GeV]"  , mtbins.size()-1, mtbins.data());
        plot1d("h_mWZbins"+s, mWZ, weight, hvec, ";m^{WZ} [GeV]"  , mtbins.size()-1, mtbins.data());
        plot1d("h_mtl3"+s, event_mTl(), weight, hvec, ";M_{T} (lep3) [GeV]" , 50, 0, 250);
        plot1d("h_mllmin"+s, mllmin, weight, hvec, ";min M_{ll} [GeV]" , 100, 0, 250);
        if (coreonly) return;

        plot1d("h_metphi"+s, metphi, weight, hvec, ";#phi(E_{T}^{miss})"  , 68, -3.4, 3.4);
        plot1d("h_rlmetphi"+s, rlmetphi, weight, hvec, ";#phi(E_{T}^{miss}) (rl)"  , 68, -3.4, 3.4);
        plot1d("h_njets20"+s, event_n_ak4jets_pt20(), weight, hvec, ";N(jets)"  , 6,  0, 6);
        plot1d("h_nvtxs"+s, event_n_vtxs_good(), weight, hvec, ";N(vtxs)"  , 80,  0, 80);

        plot1d("h_boson_mass"+s, Vmass , weight, hvec, ";M_{ll} [GeV]" , 100,  40, 140);
        plot1d("h_boson_mass_finebin"+s, Vmass , weight, hvec, ";M_{ll} [GeV]" , 160,  75, 107);
        plot1d("h_boson_eta"+s, Veta, weight, hvec, ";#eta(boson)" , 100,  -5.f, 5.f);
        plot1d("h_boson_phi"+s, Vphi, weight, hvec, ";#phi(boson)" , 64,  -3.2f, 3.2f);

        plot1d("h_min_dphijmet"+s, dphijmet, weight, hvec, ";min #Delta#phi(j, E_{T}^{miss}) ", 32,  0, 3.2);
        plot1d("h_dphi_boson_met"+s, dphiVmet, weight, hvec, ";#Delta#phi(boson, E_{T}^{miss}) ", 32,  0, 3.2);
        plot1d("h_dphi_lljets_met"+s, dphilljmet, weight, hvec, ";#Delta#phi(ll+j, E_{T}^{miss}) ", 32,  0, 3.2);
        plot1d("h_dphi_lepW_met"+s, dPhi_lepW_pTmiss(), weight, hvec, ";#Delta#phi(lepW, E_{T}^{miss}) ", 32,  0, 3.2);
        plot1d("h_dphi_ZW"+s, dPhi_Z_W(), weight, hvec, ";#Delta#phi(W, Z) ", 32,  0, 3.2);
        plot1d("h_min_dphiVjet"+s, dphiVjet , weight, hvec, ";min #Delta#phi(j, boson)", 32,  0, 3.2);;

        plot2d("h2d_mtl3_met"+s, met, event_mTl(), weight, hvec, ";E_{T}^{miss} [GeV];M_{T} (lep3) [GeV]" , 100, 0, 500, 50, 0, 250);;

        if (njet >= 1) {
          plot1d("h_jet1pt"+s, ak4jets_pt()[0], weight, hvec, ";p_{T}^{jet1} [GeV]" , 160,  0, 800);
          plot1d("h_jet1eta"+s, ak4jets_eta()[0], weight, hvec, ";#eta(jet1)"  , 96,  -4.8f, 4.8f);
        }
        if (njet >= 2) {
          plot1d("h_jet2pt"+s, ak4jets_pt()[1], weight, hvec, ";p_{T}^{jet2} [GeV]" , 160,  0, 800);
          plot1d("h_jet2eta"+s, ak4jets_eta()[1], weight, hvec, ";#eta(jet2)"  , 96,  -4.8f, 4.8f);
        }

        // Additional binnings used to derive gamma to ll boson pT reweighting
        plot1d("h_boson_pt_b0"+s, Vpt, weight, hvec, ";p_{T}(boson) [GeV]" , ptbin0.size()-1, ptbin0.data());
        plot1d("h_boson_pt_b1"+s, Vpt, weight, hvec, ";p_{T}(boson) [GeV]" , ptbin1.size()-1, ptbin1.data());
        plot1d("h_boson_aeta"+s, fabs(Veta), weight, hvec, ";#eta(boson)" , 26,  0.f, 2.6f);

        plot1d("h_ht"+s, ak4jets_HT() , weight, hvec, ";HT [GeV]" , 150,  0, 750);

        // Generator quantities, MC only
        if (!is_data) {
          plot1d("h_genmet"+s, genmet_pTmiss(), weight, hvec, ";gen-E_{T}^{miss} [GeV]"  , 160,  0, 800);
          plot1d("h_genmetphi"+s, genmet_phimiss(), weight, hvec, ";#phi(gen-E_{T}^{miss})"  , 64,  -3.2f, 3.2f);
        }

        // Dilepton quantities
        if (is3l) {
          plot1d("h_nlep"+s, leptons_id().size(), weight, hvec, ";N(lep)"  , 4,  0, 4);
          plot1d("h_llid"+s, dilepton_id(), weight, hvec, ";ll ID"  , 36,  -180, 180);
          plot1d("h_ptll"+s, dilepton_pt(), weight, hvec, ";p_{T}^{ll} [GeV]"  , 160,  0, 800);
          plot1d("h_etall"+s, dilepton_eta(), weight, hvec, ";#eta(ll)"  , 64,  -3.2f, 3.2f);
          plot1d("h_mll"+s, dilepton_mass() , weight, hvec, ";M_{ll} [GeV]" , 125,  0, 500);

          plot1d("h_mll13"+s, mll13 , weight, hvec, ";M_{l1l3} [GeV]" , 125,  0, 500);
          plot1d("h_mll23"+s, mll23 , weight, hvec, ";M_{l2l3} [GeV]" , 125,  0, 500);
          plot1d("h_mllx3"+s, mll13 , weight, hvec, ";M_{lxl3} [GeV]" , 125,  0, 500);
          plot1d("h_mllx3"+s, mll23 , weight, hvec, ";M_{lxl3} [GeV]" , 125,  0, 500);
          plot1d("h_mllclose"+s, mllclose , weight, hvec, ";M_{ll} (close) [GeV]" , 125,  0, 500);
          plot1d("h_mllsfx3"+s, mllsfx3 , (mllsfx3 > 0) * weight, hvec, ";M_{ll} (SF, close) [GeV]" , 125,  0, 500);

          plot1d("h_lepZ1pt"+s, ((lep1pt > lep2pt)? lep1pt : lep2pt)  , weight, hvec, ";p_{T}(lepZ1) [GeV]"  , 25,  0, 500);
          plot1d("h_lepZ2pt"+s, ((lep1pt > lep2pt)? lep2pt : lep1pt)  , weight, hvec, ";p_{T}(lepZ2) [GeV]"  , 20,  0, 400);
          plot1d("h_lepWpt"+s,   lep3pt  , weight, hvec, ";p_{T}(lepW) [GeV]"  , 20,  0, 400);
          plot1d("h_lep1pt"+s,   leptons_pt().at(0)  , weight, hvec, ";p_{T}(lep1) [GeV]"  , 25,  0, 500);
          plot1d("h_lep2pt"+s,   leptons_pt().at(1)  , weight, hvec, ";p_{T}(lep2) [GeV]"  , 20,  0, 400);
          plot1d("h_lep3pt"+s,   leptons_pt().at(2)  , weight, hvec, ";p_{T}(lep3) [GeV]"  , 20,  0, 400);
          plot1d("h_lep1eta"+s,  leptons_eta().at(0) , weight, hvec, ";#eta(lep1)"        , 36, -2.4, 2.4);
          plot1d("h_lep2eta"+s,  leptons_eta().at(1) , weight, hvec, ";#eta(lep2)"        , 36, -2.4, 2.4);

          plot1d("h_lepZ1pt_finebin"+s, ((lep1pt > lep2pt)? lep1pt : lep2pt)  , weight, hvec, ";p_{T}(lepZ1) [GeV]"  , 40,  0, 200);
          plot1d("h_lepZ2pt_finebin"+s, ((lep1pt > lep2pt)? lep2pt : lep1pt)  , weight, hvec, ";p_{T}(lepZ2) [GeV]"  , 40,  0, 200);
          plot1d("h_lepWpt_finebin"+s,   lep3pt  , weight, hvec, ";p_{T}(lepW) [GeV]"  , 40,  0, 200);
          plot1d("h_lep1pt_finebin"+s,   leptons_pt().at(0)  , weight, hvec, ";p_{T}(lep1) [GeV]"  , 40,  0, 200);
          plot1d("h_lep2pt_finebin"+s,   leptons_pt().at(1)  , weight, hvec, ";p_{T}(lep2) [GeV]"  , 40,  0, 200);
          plot1d("h_lep3pt_finebin"+s,   leptons_pt().at(2)  , weight, hvec, ";p_{T}(lep3) [GeV]"  , 40,  0, 200);
        }

      };

      // Assign weights
      double weight_unadjust = weight;
      weight *= capVal(event_wgt_adjust);
      evtwgt_sys[systype.Data()] = weight;
      if (!is_data && systype == "Nominal") {
        evtwgt_sys["PDFScaleUp"] = weight_unadjust * capVal(event_wgt_adjust * event_wgt_adjustment_PDFScaleUp());
        evtwgt_sys["PDFScaleDn"] = weight_unadjust * capVal(event_wgt_adjust * event_wgt_adjustment_PDFScaleDn());
        evtwgt_sys["QCDScaleUp"] = weight_unadjust * capVal(event_wgt_adjust * event_wgt_adjustment_QCDScaleUp());
        evtwgt_sys["QCDScaleDn"] = weight_unadjust * capVal(event_wgt_adjust * event_wgt_adjustment_QCDScaleDn());
        evtwgt_sys["PDFReplicaUp"] = weight_unadjust * capVal(event_wgt_adjust * (useNNPDF30? event_wgt_adjustment_NNPDF30_PDFReplicaUp() : event_wgt_adjustment_PDFReplicaUp()));
        evtwgt_sys["PDFReplicaDn"] = weight_unadjust * capVal(event_wgt_adjust * (useNNPDF30? event_wgt_adjustment_NNPDF30_PDFReplicaDn() : event_wgt_adjustment_PDFReplicaDn()));
        evtwgt_sys["PUJetIdEffUp"] = weight * event_wgt_SFs_PUJetId_EffUp() / event_wgt_SFs_PUJetId();
        evtwgt_sys["PUJetIdEffDn"] = weight * event_wgt_SFs_PUJetId_EffDn() / event_wgt_SFs_PUJetId();
        evtwgt_sys["BtagSFUp"] = weight * event_wgt_SFs_btagging_EffUp() / event_wgt_SFs_btagging();
        evtwgt_sys["BtagSFDn"] = weight * event_wgt_SFs_btagging_EffDn() / event_wgt_SFs_btagging();
        evtwgt_sys["ElectronEff_AltMCUp"] = weight * event_wgt_SFs_electrons_AltMCUp() / event_wgt_SFs_electrons();
        evtwgt_sys["ElectronEff_AltMCDn"] = weight * event_wgt_SFs_electrons_AltMCDn() / event_wgt_SFs_electrons();
        evtwgt_sys["ElectronEff_SystUp"] = weight * event_wgt_SFs_electrons_SystUp() / event_wgt_SFs_electrons();
        evtwgt_sys["ElectronEff_SystDn"] = weight * event_wgt_SFs_electrons_SystDn() / event_wgt_SFs_electrons();
        evtwgt_sys["ElectronEff_StatUp"] = weight * event_wgt_SFs_electrons_StatUp() / event_wgt_SFs_electrons();
        evtwgt_sys["ElectronEff_StatDn"] = weight * event_wgt_SFs_electrons_StatDn() / event_wgt_SFs_electrons();
        evtwgt_sys["MuonEff_AltMCUp"] = weight * event_wgt_SFs_muons_AltMCUp() / event_wgt_SFs_muons();
        evtwgt_sys["MuonEff_AltMCDn"] = weight * event_wgt_SFs_muons_AltMCDn() / event_wgt_SFs_muons();
        evtwgt_sys["MuonEff_SystUp"] = weight * event_wgt_SFs_muons_SystUp() / event_wgt_SFs_muons();
        evtwgt_sys["MuonEff_SystDn"] = weight * event_wgt_SFs_muons_SystDn() / event_wgt_SFs_muons();
        evtwgt_sys["MuonEff_StatUp"] = weight * event_wgt_SFs_muons_StatUp() / event_wgt_SFs_muons();
        evtwgt_sys["MuonEff_StatDn"] = weight * event_wgt_SFs_muons_StatDn() / event_wgt_SFs_muons();
        evtwgt_sys["PUUp"] = weight * event_wgt_PUUp() / event_wgt();
        evtwgt_sys["PUDn"] = weight * event_wgt_PUDn() / event_wgt();
        evtwgt_sys["L1PrefiringUp"] = weight * event_wgt_L1PrefiringUp() / event_wgt();
        evtwgt_sys["L1PrefiringDn"] = weight * event_wgt_L1PrefiringDn() / event_wgt();

        if (applyExternalAsMZ) {
          evtwgt_sys["AsMZUp"] = weight_unadjust * capVal(event_wgt_adjust * getValFromHist1D(hmap_externwgt["AsMZUp"], genpromptparticles_sump4_pt()));
          evtwgt_sys["AsMZDn"] = weight_unadjust * capVal(event_wgt_adjust * getValFromHist1D(hmap_externwgt["AsMZDn"], genpromptparticles_sump4_pt()));
        } else {
          evtwgt_sys["AsMZUp"] = weight_unadjust * capVal(event_wgt_adjust * (useNNPDF30? event_wgt_adjustment_NNPDF30_AsMZUp() : event_wgt_adjustment_AsMZUp()));
          evtwgt_sys["AsMZDn"] = weight_unadjust * capVal(event_wgt_adjust * (useNNPDF30? event_wgt_adjustment_NNPDF30_AsMZDn() : event_wgt_adjustment_AsMZDn()));
        }
        if (applyExternalPythiaScale) {
          evtwgt_sys["PythiaScaleUp"] = weight_unadjust * capVal(event_wgt_adjust * getValFromHist1D(hmap_externwgt["PythiaScaleUp"], genpromptparticles_sump4_pt()));
          evtwgt_sys["PythiaScaleDn"] = weight_unadjust * capVal(event_wgt_adjust * getValFromHist1D(hmap_externwgt["PythiaScaleDn"], genpromptparticles_sump4_pt()));
        } else {
          evtwgt_sys["PythiaScaleUp"] = weight_unadjust * capVal(event_wgt_adjust * event_wgt_adjustment_PythiaScaleUp());
          evtwgt_sys["PythiaScaleDn"] = weight_unadjust * capVal(event_wgt_adjust * event_wgt_adjustment_PythiaScaleDn());
        }
      }
      for (auto wgt : evtwgt_sys) {
        if (isnan(wgt.second)) {
          cout << "[WARNING] >> Systematic " << wgt.first << " is seen to have NaN value!!" << endl;
          wgt.second = 0;
        }
      }

      string metsuf = "_fullMET";

      // Simplified category -- version 1
      if (lepcat == "_eee"  || lepcat == "_mumue" ) lepcat = "_lle";
      if (lepcat == "_eemu" || lepcat == "_mumumu") lepcat = "_llmu";

      double origwgt = weight;

      const bool fillFullMETplots = false;
      if (fillFullMETplots) {

        fillhists(metsuf);
        fillhists(metsuf+jetcat);
        fillhists(metsuf+lepcat);
        fillhists(metsuf+jetcat+lepcat);

        for (string isys : systematics) {
          fillhists("_"+isys+metsuf, true);
          fillhists("_"+isys+metsuf+lepcat, true);
          fillhists("_"+isys+metsuf+jetcat+lepcat, true);
        }

        // Fill EW uncertainties for the Nominal only study
        if (!is_data && systematics.size() == 1 && systematics[0] == "Nominal") {
          weight = origwgt * wgt_EWup;
          fillhists("_EWUp"+metsuf);
          fillhists("_EWUp"+metsuf+jetcat+lepcat);
          weight = origwgt * wgt_EWdn;
          fillhists("_EWDn"+metsuf);
          fillhists("_EWDn"+metsuf+jetcat+lepcat);
          weight = origwgt;
        }
      }

      metsuf = "_final";

      fillhists(metsuf);
      fillhists(metsuf+jetcat);
      fillhists(metsuf+lepcat);
      fillhists(metsuf+jetcat+lepcat);

      if (!is_data && systematics.size() == 1 && systematics[0] == "Nominal") {
        weight = origwgt * wgt_EWup;
        fillhists("_EWUp"+metsuf);
        fillhists("_EWUp"+metsuf+jetcat+lepcat);
        weight = origwgt * wgt_EWdn;
        fillhists("_EWDn"+metsuf);
        fillhists("_EWDn"+metsuf+jetcat+lepcat);
        weight = origwgt;
      }

      if (produceResultTree) {
        // Common branches
        evt_weight = weight;

        llid = dilepton_id();
        lep3eta = leptons_eta().at(ilep3);
        lep3phi = leptons_phi().at(ilep3);
        lep3mass = leptons_mass().at(ilep3);

        int j = 0;
        for (; j < std::min(njet, 2U); ++j) {
          jet_pt[j]= ak4jets_pt().at(j);
          jet_eta[j]= ak4jets_eta().at(j);
          jet_phi[j]= ak4jets_phi().at(j);
          jet_mass[j]= ak4jets_mass().at(j);
        }
        for (; j < 2; ++j) {
          jet_pt[j] = -999.;
          jet_eta[j] = -999.;
          jet_phi[j] = -999.;
          jet_mass[j] = -999.;
        }

        // Weight branches
        for (string isys : systematics) {
          fouts[isys]->cd();
          evt_weight = evtwgt_sys.at(isys);
          touts[isys]->Fill();
        }
      }

    } //event loop


    delete file;
  }//file loop

  // Produce the ratio plots from the numerator and denominator hists
  for (const auto& h : hvec) {
    if (!produceResultPlots) break;
    if (h.first.find("hnum") != 0) continue;
    if (h.first.find("2d") == string::npos)
      moveOverFlowToLastBin1D(h.second);
    string hname = h.first;
    hname.erase(0, 4);
    TH1F* h_ratio = (TH1F*) h.second->Clone(("ratio"+hname).c_str());
    h_ratio->SetDirectory(0);
    h_ratio->Divide(h_ratio, hvec.at("hden"+hname), 1, 1, "B");
    // h_ratio->Divide(hvec.at("hden"+hname));  // not eff

    const string dirname = "effstudy_eff";
    TDirectory* dir = (TDirectory*) fout->Get(dirname.c_str());
    if (dir == nullptr) dir = fout->mkdir(dirname.c_str());
    dir->cd();
    h_ratio->Write();

    dir = (TDirectory*) fout->Get("effstudy_num");
    if (dir == nullptr) dir = fout->mkdir("effstudy_num");
    dir->cd();
    h.second->Write(("hnum"+hname).c_str());

    dir = (TDirectory*) fout->Get("effstudy_den");
    if (dir == nullptr) dir = fout->mkdir("effstudy_den");
    dir->cd();
    hvec.at("hden"+hname)->Write(("hden"+hname).c_str());
  }

  // Make conditionalized 2D distributions
  for (const auto& h : hvec) {
    if (!produceResultPlots) break;
    if (h.first.find("h2d") == string::npos) continue;
    if (h.first.find("genparts_") == string::npos) continue;

    const string dirname = "OffShell";
    TDirectory* dir = (TDirectory*) fout->Get(dirname.c_str());
    if (dir == nullptr) dir = fout->mkdir(dirname.c_str());
    dir->cd();
    // h.seocnd->Write();
    if (h.first.find("Up_") != string::npos || h.first.find("Dn_") != string::npos) {
      TH2* hnew = (TH2F*) h.second->Clone(Form("%s_norm", h.second->GetName()));
      conditionalizeHistInX(hnew);
      hnew->Write();
    }
  }

  // Put hists into different folders based on their suffixes
  for (auto& h : hvec) {
    if (!produceResultPlots) break;
    if (h.first.find("hnum") == 0 || h.first.find("hden") == 0) continue;
    if (h.first.find("phi") == string::npos && h.first.find("h_") == 0)
      moveOverFlowToLastBin1D(h.second);
    string dirname = "OffShell";
    vector<string> jetsufs = {"_eq0j", "_eq1j", "_eq2j", "_ge2j", "_vbf", "_eq0j_lowdphi", };
    vector<string> lepsufs = {"_eel", "_mumul", "_lle", "_llmu", "_emul", "_sssf", };
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
    if (separateFileForTree) {
      for (string isys : systematics) {
        fouts[isys]->cd();
        touts[isys]->Write();
        fouts[isys]->Close();
      }
    } else {
      fout->cd();
      touts[0]->Write();
    }
  }

  if (produceResultPlots && fout) fout->Close();

  return 0;
}
