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

#include "SkimTree.h"
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
bool applyZGweightFromLLG = false;
bool produceResultPlots = true;

const bool makeStepPlots = false;
const bool produceFilterList = false;
const bool produceResultTree = false;
const bool separateFileForTree = true;
const bool usStyleTree = true;
const bool fastRun = false;
const bool blind = false;

vector<float> nvtxscale(100, 0);

std::ofstream ofile;

int ScanChain(TString indir, TString sample, TString tag, TString systype = "", TString specifiers = "", TString extrargs = "") {

  TChain *ch = new TChain("SkimTree");
  TString files_in = Form("%s/%s*.root", indir.Data(), sample.Data());
  ch->Add(files_in);
  cout << ">> Adding " << files_in << " into the chain." << endl;

  TString samplever = tag(0, min(tag.Index("_",tag.Index("_")+1),5));
  double skimver = TString(samplever).ReplaceAll("v","").ReplaceAll("_",".").Atof();
  cout << ">> Runing from skimver " << samplever  << ", i.e. " << skimver << endl;
  // skimver = 1.0;

  // FIXME: for synchronization debug only
  // vector<float> mlllist = {79.793434143066, 85.465690612793, 76.578109741211, 94.826354980469,};
  vector<float> mlllist{0};
  vector<float> ptlllist{0};
  vector<int> njetlist{0};
  if (false) {
    TFile* fsync = new TFile("/home/users/sicheng/working/HZZlooper/usskimLooper/example/finaltree_data_Nominal_2018.root");
    TTree* ft = (TTree*) fsync->Get("FinalTree");
    float mll(0);
    float ptll(0);
    unsigned nj30(0);
    ft->SetBranchStatus("*", 0);
    ft->SetBranchStatus("dilepton_mass", 1);
    ft->SetBranchStatus("dilepton_pt", 1);
    ft->SetBranchStatus("n_ak4jets_pt30", 1);
    ft->SetBranchAddress("dilepton_mass", &mll);
    ft->SetBranchAddress("dilepton_pt", &ptll);
    ft->SetBranchAddress("n_ak4jets_pt30", &nj30);
    for (int evt = 0; evt < ft->GetEntries(); evt += 2) {
      ft->GetEntry(evt);
      mlllist.push_back(mll);
      ptlllist.push_back(ptll);
      njetlist.push_back(nj30);
    }
  }

  bool is_data = (sample.BeginsWith("Run201")) || sample.BeginsWith("photon") || sample.BeginsWith("EGamma");
  bool is_gjets = sample.BeginsWith("GJets");

  int year = (tag.Contains("2016"))? 2016 : (tag.Contains("2017"))? 2017 : (tag.Contains("2018"))? 2018 : -1;
  if (year < 0) cout << ">> WARNING: not able to determine year from tag: " << tag << endl;

  vector<string> systematics = {"Nominal"};
  vector<string> systnames = {"PhoTriggerEff", "TransferFactor", "PDFScale", "QCDScale", "AsMZ", "PDFReplica", "PhoEff", "PU", "BtagSF", "L1Prefiring"};
  for (string isys : systnames) {
    systematics.push_back(isys+"Up");
    systematics.push_back(isys+"Dn");
  }

  if (is_data) systematics = vector<string>{"Nominal",  "TransferFactorUp", "TransferFactorDn"};
  if (is_data && sample.BeginsWith("EGamma")) {
    systematics.push_back("ElecToPhotonScaleUp");
    systematics.push_back("ElecToPhotonScaleDn");
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
    TString systname = "Nominal";
    TString fileout = Form("output/%s/%s_%s.root", tag.Data(), outname.Data(), systname.Data());
    gSystem->Exec(Form("mkdir -p output/%s", tag.Data()));
    fout = new TFile(fileout, "RECREATE");
    cout << ">> Outputting to " << fileout << endl;
  }

  if (tag.Contains("nvtxrwgt")) {
    doNvtxReweight = true;
    cout << ">> doNvtxReweight is set to true." << endl;
  }
  else if (tag.Contains("rwgtd")) {
    doBosonPtReweight = true;
    cout << ">> doBosonPtReweight is set to true." << endl;
  }
  if (tag.Contains("flateta")) {
    doBosonEtaReweight = true;
    cout << ">> doBosonEtaReweight is set to true." << endl;
  }
  extendEEphoton2j = tag.Contains("ee2j");
  tightDphiIn2j = tag.Contains("2jdphi0p5");
  tightMETin2j = tag.Contains("2jmet140");
  doHllgTest = tag.Contains("Hllgtest");

  if (tag.Contains("all2jsel")) {
    extendEEphoton2j = doBosonEtaReweight = tightMETin2j = tightDphiIn2j = true;
  }

  map<string,TH1*> hvec;

  unsigned nEventsTotal = 0;
  unsigned nEventsChain = ch->GetEntries();

  // Result tree producer
  map<string, TFile*> fouts;
  map<string, TTree*> touts;
  float evt_weight;
  float ll_pt, ll_eta, ll_phi, ll_mass;
  float ptmiss, phimiss, mT, M_ZZ;
  float dphijmet, dphiVmet, dphilljmet;
  float trigwgt(1.), lumiwgt(1.), sf_elec(1.), sf_muon(1.), sf_photon(1.), sf_btag(1.), sf_pujetid(1.), sf_zgnjet(1.);
  float DjjVBF(-1.), DjjVBFL1(-1.), DjjVBFa2(-1.), DjjVBFa3(-1.), DjjVBFL1ZGs(-1.);

  int ll_id(0), lepton_cat(-1), jet_cat(-1);
  unsigned njet(0), nvtxs(0);
  float jet_pt[32], jet_eta[32], jet_phi[32], jet_mass[32];
  unsigned run(0), lumi(0), event(0);
  unsigned njet_m60(0), nak8jet200(0), nak8jet200_m60(0), nak8jet200_m140(0);
  float ak8jet_pt(0.), ak8jet_eta(-999.), ak8jet_phi(-999.), ak8jet_mass(0.);

  float dijet_mass(0.), dijet_pt(0.), dijet_dEta(-99.), dijet_dPhi(-99.);

  bool ph_convveto;
  bool ph_passPFid;
  bool ph_isEB;
  bool is_vbf(false);

  float tf_sgtoll(1.), tf_sgtoee, tf_sgtomumu, tf_etog(1.);
  float tferr_sgtoll(0.), tferr_etog(0.), tferr_sgtoee(0.), tferr_sgtomumu(0.);
  float ph_pt, ph_eta, ph_phi;
  float ph_r9, ph_sieie, ph_sipip;
  float ph_Emip, ph_seedtime, ph_e4oe1;

  map<string,float> evtwgt_sys;

  if (produceResultTree) {
    for (string isys : systematics) {
      if (separateFileForTree) {
        gSystem->Exec(Form("mkdir -p outtree/%s", tag.Data()));
        TString otname = Form("outtree/%s/finaltree_%s_%s.root", tag.Data(), outname.Data(), isys.c_str());
        cout << ">> Outputting result tree to " << otname << endl;
        fouts[isys] = new TFile(otname, "RECREATE");
        fouts[isys]->cd();
      }
      bool extraBranches = (isys == "Nominal");

      TTree* tout = new TTree("FinalTree", "");
      if (usStyleTree) {
        tout->Branch("weight", &evt_weight);
        tout->Branch("mTZZ", &mT);

        tout->Branch("dilepton_id", &ll_id);
        tout->Branch("dilepton_pt", &ll_pt);
        tout->Branch("dilepton_eta", &ll_eta);
        tout->Branch("dilepton_mass", &ll_mass);

        tout->Branch("pTmiss", &ptmiss);
        tout->Branch("phimiss", &phimiss);

        tout->Branch("n_ak4jets_pt30", &njet);
        tout->Branch("n_ak4jets_pt30_mass60", &njet_m60);

        tout->Branch("dijet_mass", &dijet_mass);
        tout->Branch("dijet_pt", &dijet_pt);
        tout->Branch("dijet_dEta", &dijet_dEta);
        tout->Branch("dijet_dPhi", &dijet_dPhi);

        tout->Branch("ak4jet_leading_pt", jet_pt);
        tout->Branch("ak4jet_leading_eta", jet_eta);
        tout->Branch("ak4jet_leading_phi", jet_phi);
        tout->Branch("ak4jet_leading_mass", jet_mass);

        tout->Branch("ak4jet_subleading_pt", jet_pt+1);
        tout->Branch("ak4jet_subleading_eta", jet_eta+1);
        tout->Branch("ak4jet_subleading_phi", jet_phi+1);
        tout->Branch("ak4jet_subleading_mass", jet_mass+1);

        tout->Branch("n_ak8jets_pt200", &nak8jet200);
        tout->Branch("n_ak8jets_pt200_mass60to110", &nak8jet200_m60);
        tout->Branch("n_ak8jets_pt200_mass140", &nak8jet200_m140);

        tout->Branch("ak8jet_leading_pt", ak8jet_pt);
        tout->Branch("ak8jet_leading_eta", ak8jet_eta);
        tout->Branch("ak8jet_leading_mass", ak8jet_mass);

        tout->Branch("DjjVBF", &DjjVBF);
        tout->Branch("DjjVBFL1", &DjjVBFL1);
        tout->Branch("DjjVBFa2", &DjjVBFa2);
        tout->Branch("DjjVBFa3", &DjjVBFa3);
        tout->Branch("DjjVBFL1ZGs", &DjjVBFL1ZGs);

        // Extra branches for studies
        if (extraBranches) {
          tout->Branch("tf_gamma_to_ll", &tf_sgtoll);
          tout->Branch("mindphi_jet_met", &dphijmet);
          tout->Branch("dphi_boson_met", &dphiVmet);
          tout->Branch("dphi_lljets_met", &dphilljmet);
          tout->Branch("run", &run);
          tout->Branch("lumi", &lumi);
          tout->Branch("event", &event);
          tout->Branch("is_vbf", &is_vbf);
        }
      } else {
        tout->Branch("weight", &evt_weight);
        tout->Branch("ll_pt", &ll_pt);
        tout->Branch("ll_eta", &ll_eta);
        tout->Branch("ll_phi", &ll_phi);
        tout->Branch("ll_mass", &ll_mass);
        tout->Branch("ptmiss", &ptmiss);
        tout->Branch("ptmiss_phi", &phimiss);
        tout->Branch("mT", &mT);
        tout->Branch("lepton_cat", &lepton_cat);
        tout->Branch("jet_cat", &jet_cat);

        tout->Branch("mindphi_jet_met", &dphijmet);
        tout->Branch("dphi_boson_met", &dphiVmet);
        tout->Branch("dphi_lljets_met", &dphilljmet);
        tout->Branch("mZZ", &M_ZZ);
        tout->Branch("DjjVBF", &DjjVBF);
        tout->Branch("num_pv_good", &nvtxs);
        tout->Branch("weight_trigger", &trigwgt);
        tout->Branch("weight_luminosity", &lumiwgt);
        tout->Branch("weight_SF_elec", &sf_elec);
        tout->Branch("weight_SF_muon", &sf_muon);
        tout->Branch("weight_SF_photon", &sf_photon);
        tout->Branch("weight_SF_btag", &sf_btag);
        tout->Branch("weight_SF_pujetid", &sf_pujetid);

        tout->Branch("run", &run);
        tout->Branch("lumi", &lumi);
        tout->Branch("event", &event);

        tout->Branch("photon_convveto", &ph_convveto);
        tout->Branch("photon_passPFid", &ph_passPFid);
        tout->Branch("photon_isEB", &ph_isEB);

        tout->Branch("photon_pt", &ph_pt);
        tout->Branch("photon_eta", &ph_eta);
        tout->Branch("photon_phi", &ph_phi);
        tout->Branch("photon_r9", &ph_r9);
        tout->Branch("photon_sieie", &ph_sieie);
        tout->Branch("photon_sipip", &ph_sipip);
        tout->Branch("photon_Emip", &ph_Emip);
        tout->Branch("photon_seedtime", &ph_seedtime);

        tout->Branch("tf_gamma_to_ll", &tf_sgtoll);
        // tout->Branch("tf_gamma_to_ee", &tf_sgtoee);
        // tout->Branch("tf_gamma_to_mumu", &tf_sgtomumu);
        tout->Branch("jet_size", &njet);
        tout->Branch("jet_pt", jet_pt, "jet_pt[jet_size]/F");
        tout->Branch("jet_eta", jet_eta, "jet_eta[jet_size]/F");
        tout->Branch("jet_phi", jet_phi, "jet_phi[jet_size]/F");
        tout->Branch("jet_mass", jet_mass, "jet_mass[jet_size]/F");

        if (indir.Contains("/SingleLeptonEvents/")) {
          tout->Branch("tf_el_to_gamma", &tf_etog);
        }
      }

      touts[isys] = tout;
    }
  }

  if (produceFilterList && is_data) {
    gSystem->Exec(Form("mkdir -p outfile/%s", tag.Data()));
    ofile.open(Form("outfile/%s/photonFilterList_%s.txt", tag.Data(), sample.Data()));
  }

  // Make photon-pt bins for reweight
  vector<float> ptbin0, ptbin1;
  std::tie(ptbin0, ptbin1) = getPtBins();

  // For the DjjVBF calculation
  TSpline3* h_vbfcval = fetchHistCopy<TSpline3>("data/SmoothKDConstant_m4l_DjjVBF_13TeV.root", "sp_gr_varReco_Constant_Smooth");
  TSpline3* h_vbfgval_L1 = fetchHistCopy<TSpline3>("data/gConstant_VBF_L1.root", "sp_tgfinal_VBF_SM_over_tgfinal_VBF_L1");
  TSpline3* h_vbfgval_a2 = fetchHistCopy<TSpline3>("data/gConstant_VBF_g2.root", "sp_tgfinal_VBF_SM_over_tgfinal_VBF_g2");
  TSpline3* h_vbfgval_a3 = fetchHistCopy<TSpline3>("data/gConstant_VBF_g4.root", "sp_tgfinal_VBF_SM_over_tgfinal_VBF_g4");
  TSpline3* h_vbfgval_L1ZGs = fetchHistCopy<TSpline3>("data/gConstant_VBF_L1Zgs.root", "sp_tgfinal_VBF_SM_photoncut_over_tgfinal_VBF_L1Zgs");

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
    bool useNNPDF30 = isNLO && !(isZG_nlo_incl || isZG_nlo_ptG130);
    bool applyKfactor = filename.Contains("ZZTo") || filename.Contains("WZTo");


    TFile *file = new TFile(filename);
    TTree *tree = (TTree*)file->Get("SkimTree");
    st.Init(tree);

    for (unsigned int evt = 0; evt < tree->GetEntriesFast(); ++evt) {

      st.GetEntry(evt);
      nEventsTotal++;

      SkimTree::progress(nEventsTotal, nEventsChain);

      double weight = event_wgt();
      if (is_data) weight = 1;
      if (fastRun && event_pTmiss() < 125) continue;

      bool isdilep(false), isgamma(false), isllg(false), islg(false), is1el(false);

      string lepcat = "_lepcat";
      lepton_cat = -1;
      if (indir.Contains("/DileptonEvents/")) {
        // skimtype = Dilepton;
        isdilep = true;
        switch (dilepton_id()) {
          case -121:
            lepton_cat = 0;
            lepcat = "_ee";
            break;
          case -169:
            lepton_cat = 1;
            lepcat = "_mumu";
            break;
          case -143:
            lepton_cat = 2;
            lepcat = "_emu";
            break;
          default:
            lepton_cat = -1;
            cout << "WARNING: Invalid dilepton ID: " << dilepton_id() << "!!" << endl;
        }
        trigwgt = event_wgt_triggers_SingleLepton() + event_wgt_triggers_Dilepton();
        if (trigwgt == 0.f) continue;
        if (tag.Contains("emuCR") && lepton_cat != 2) continue;
        if (tag.Contains("2lSR") && (lepton_cat | 1) != 1) continue;
      }
      if (indir.Contains("/SinglePhotonEvents/")) {
        // skimtype = SinglePhoton;
        isgamma = true;
        lepton_cat = 3;
        lepcat = "_gamma";
        trigwgt = event_wgt_triggers();

        if (trigwgt == 0.f) continue;
        weight *= trigwgt;
      }
      if (indir.Contains("/SinglePhotonEvents/")) {
        // skimtype = SinglePhoton;
        islg = true;
        lepton_cat = 3;
        lepcat = "_gamma";
        trigwgt = event_wgt_triggers();

        if (trigwgt == 0.f) continue;
        weight *= trigwgt;
      }
      if (indir.Contains("/LLGEvents/")) {
        // skimtype = LLG;
        isllg = true;
        lepton_cat = 4;
        switch (dilepton_id()) {
          case -121:
            lepcat = "_eeg";
            break;
          case -169:
            lepcat = "_mumug";
            break;
          default:
            lepcat = "_mumug";
            cout << "WARNING: Invalid dilepton ID: " << dilepton_id() << "!!" << endl;
        }
        trigwgt = event_wgt_triggers_SingleLepton() + event_wgt_triggers_Dilepton();
        if (trigwgt == 0.f) continue;
        applyZGweightFromLLG = tag.Contains("zgwgtd");
      }
      if (indir.Contains("/SingleLeptonEvents/")) {
        is1el = true;
        lepton_cat = 5;
        lepcat = "_1el";

        trigwgt = event_wgt_triggers();
        weight *= trigwgt;

        // Only loop over single-e data events
        if (abs(lepton_id()) != 11) continue;
      }

      double event_wgt_SFs = 1.0;
      double event_wgt_sample = 1.0;
      if (!is_data) {
        // if (year == 2018 && event_wgt_pileup() > 100) continue;
        // weight *= event_wgt_pileup(); // pileup weight included in event_wgt for US samples
        event_wgt_SFs = event_wgt_SFs_electrons() * event_wgt_SFs_muons() * event_wgt_SFs_photons();
        event_wgt_SFs *= event_wgt_SFs_PUJetId() * event_wgt_SFs_btagging();
        event_wgt_SFs = std::min(event_wgt_SFs, 5.0);
        weight *= event_wgt_SFs;

        event_wgt_sample = getExtSampleWeight(filename, year);
        weight *= event_wgt_sample;

        if (useNNPDF30) weight *= event_wgt_adjustment_NNPDF30();
        if (applyKfactor) weight *= KFactor_EW_NLO_qqVV_Bkg_Nominal() * KFactor_QCD_NNLO_qqVV_Bkg_Nominal();

        // Fill histogram of the weights
        double scale = pow(10.0, 1 - ceil(log10(fabs(weight))));
        double wgt_range = 10 * round(weight * scale) / scale;
        plot1d("h_event_wgts", weight, 1, hvec, "; event weights" , 1000,  -wgt_range, wgt_range);
        plot1d("h_event_wgts_raw", event_wgt(), 1, hvec, "; event weights" , 1000,  -wgt_range, wgt_range);
        // plot1d("h_event_wgts_pileup", event_wgt_pileup(), 1, hvec, "; event weights" , 1000,  0, 100);
        plot1d("h_event_wgts_SFs", event_wgt_SFs, 1, hvec, "; event weights" , 100,  0, 5);
        plot1d("h_event_wgt_SFs_PUJetId",  event_wgt_SFs_PUJetId(), 1, hvec, "; event weights PUJetId" , 100,  0, 5);
        plot1d("h_event_wgt_SFs_btagging",  event_wgt_SFs_btagging(), 1, hvec, "; event weights btagging" , 100,  0, 5);
        plot1d("h_event_wgt_SFs_electrons",  event_wgt_SFs_electrons(), 1, hvec, "; event weights electrons" , 100,  0, 5);
        plot1d("h_event_wgt_SFs_muons",  event_wgt_SFs_muons(), 1, hvec, "; event weights muons" , 100,  0, 5);
        plot1d("h_event_wgt_SFs_photons",  event_wgt_SFs_photons(), 1, hvec, "; event weights photons" , 100,  0, 5);
        plot1d("h_event_wgt_adjustment",  event_wgt_adjustment_NNPDF30(), 1, hvec, "; event weights adjustment" , 100,  0, 5);
      }

      // Fill histogram of trigger weights
      if (isgamma) {
        plot1d("h_event_wgts_trigs", event_wgt_triggers(), 1, hvec, "; event weights" , 1000,  0, 1000);
      }

      double wgt = weight;
      auto fillPhotonHitMap = [&](string suf) {
        if (!produceResultPlots) return;
        if (!isgamma && !isllg) return;
        plot1d("h_lepton_cat"+suf, lepton_cat, weight, hvec, "; event cat" , 8,  -1, 7);
        plot1d("h_photon_eta"+suf, photon_eta(), wgt, hvec, ";#eta(#gamma)" , 50, -2.5f, 2.5f);
        plot1d("h_photon_pt"+suf, photon_pt(), wgt, hvec, ";p_{T}($gamma) [GeV]" , 160, 0, 800);
        plot1d("h_met_phm"+suf, event_pTmiss(), wgt, hvec, ";p_{T}^{miss} [GeV]" , 160, 0, 800);
        plot1d("h_photon_sieie"+suf, photon_full5x5_sigmaIEtaIEta(), weight, hvec, ";#sigma_{i#etai#eta}(#gamma)" , 80,  0, 0.02);
        plot1d("h_photon_sipip"+suf, photon_full5x5_sigmaIPhiIPhi(), weight, hvec, ";#sigma_{i#phii#phi}(#gamma)" , 80,  0, 0.02);
        if (!photon_isEB()) {
          plot1d("h_photon_endcap_phi"+suf, photon_phi(), wgt, hvec, ";#phi(#gamma)" , 64, -3.2f, 3.2f);
          plot1d("h_photon_endcap_phif"+suf, phiFolding(photon_phi()), wgt, hvec, ";folded-#phi(#gamma)" , 51, -1.7f, 1.7f);
        }
        if (isllg) {
          plot1d("h_llg_mass"+suf, event_mllg(), weight, hvec, ";M(ll#gamma) [GeV]" , 100, 0, 250);
        }
      };
      fillPhotonHitMap("_raw");

      auto fillJetMassHists = [&](string s = "") {
        if (!produceResultPlots) return;
        if (event_n_ak4jets_pt30() > 0) {
          plot1d("h_jet1_pt"+s, ak4jets_pt()[0], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
          plot1d("h_jet1_mass"+s, ak4jets_mass()[0], weight, hvec, ";M_{T}^{jet1} [GeV]" , 80,  0, 400);
          plot1d("h_jet1_eta"+s, ak4jets_eta()[0], weight, hvec, ";#eta(jet1)"  , 96,  -4.8f, 4.8f);
          plot1d("h_jet1_nemf"+s, ak4jets_NEMF()[0], weight, hvec, ";jet neutral EM fraction" , 50,  0, 1);
        }
        if (event_n_ak4jets_pt30() > 1) {
          plot1d("h_jet2_nemf"+s, ak4jets_NEMF()[1], weight, hvec, ";jet neutral EM fraction" , 50,  0, 1);

          float mindiffMjjW(9999), mindiffMjjZ(9999), mindiffMjjtop(9999);
          for (int i = 0; i <event_n_ak4jets_pt30(); ++i) {
            for (int j = i+1; j <event_n_ak4jets_pt30(); ++j) {
              LorentzVector p4jj = LorentzVector(ak4jets_pt()[i], ak4jets_eta()[i], ak4jets_phi()[i], ak4jets_mass()[i])
                  + LorentzVector(ak4jets_pt()[j], ak4jets_eta()[j], ak4jets_phi()[j], ak4jets_mass()[j]);
              if (fabs(p4jj.M() - 80.3) < fabs(mindiffMjjW)) mindiffMjjW = (p4jj.M() - 80.3);
              if (fabs(p4jj.M() - 91.2) < fabs(mindiffMjjZ)) mindiffMjjZ = (p4jj.M() - 91.2);
              if (fabs(p4jj.M() - 173) < fabs(mindiffMjjtop)) mindiffMjjtop = (p4jj.M() - 173);
            }
          }
          plot1d("h_jet_mass_closest_W"+s, mindiffMjjW+80.3, weight, hvec, ";M_{jj} [GeV]" , 60,  0, 1200);
          plot1d("h_jet_mass_closest_Z"+s, mindiffMjjZ+91.2, weight, hvec, ";M_{jj} [GeV]" , 40,  0, 1600);
          plot1d("h_jet_mass_closest_top"+s, mindiffMjjtop+173, weight, hvec, ";M_{jj} [GeV]" , 40,  0, 2000);
        }

        if (ak8jets_mass().size() > 0) {
          plot1d("h_ak8jet1_mass"+s, ak8jets_mass()[0], weight, hvec, ";M_{T}^{ak8jet1} [GeV]" , 80,  0, 400);
          float mindiffMjW(9999), mindiffMjZ(9999), mindiffMjtop(9999);
          for (size_t j = 0; j < ak8jets_mass().size(); ++j) {
            if (fabs(ak8jets_mass()[j] - 80.3) < fabs(mindiffMjW)) mindiffMjW = (ak8jets_mass()[j] - 80.3);
            if (fabs(ak8jets_mass()[j] - 91.2) < fabs(mindiffMjZ)) mindiffMjZ = (ak8jets_mass()[j] - 91.2);
            if (fabs(ak8jets_mass()[j] - 173) < fabs(mindiffMjtop)) mindiffMjtop = (ak8jets_mass()[j] - 173);
          }
          plot1d("h_ak8jet_mass_closest_W"+s, mindiffMjW+80.3, weight, hvec, ";M_{J} [GeV]" , 80,  0, 400);
          plot1d("h_ak8jet_mass_closest_Z"+s, mindiffMjZ+91.2, weight, hvec, ";M_{J} [GeV]" , 80,  0, 400);
          plot1d("h_ak8jet_mass_closest_top"+s, mindiffMjZ+173, weight, hvec, ";M_{J} [GeV]" , 80,  0, 400);
        }
            
        float jet_nemf_inHEHFgap = -1;
        for (int j = 0; j < event_n_ak4jets_pt30(); ++j) {
          plot1d("h_jets_nemf"+s, ak4jets_NEMF()[j], weight, hvec, ";jet neutral EM fraction" , 50,  0, 1);
          float nemf_thr = std::min(1., std::max(0., -0.03405+0.001472*ak4jets_pt()[j]+0.0001892*ak4jets_pt()[j]*ak4jets_pt()[j]));
          if (2.65 < fabs(ak4jets_eta()[j]) && fabs(ak4jets_eta()[j]) < 3.139) {
            jet_nemf_inHEHFgap = ak4jets_NEMF()[j];
            plot1d("h_jets_nemf_HEHFgap"+s, ak4jets_NEMF()[j], weight, hvec, ";jet neutral EM fraction" , 50,  0, 1);
            plot1d("h_jets_pt_HEHFgap"+s, ak4jets_pt()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
            plot1d("h_jets_empt_HEHFgap"+s, ak4jets_pt()[j] * ak4jets_NEMF()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
            if (ak4jets_NEMF()[j] > 0.6) {
              plot1d("h_jets_pt_HEHFgap_nemfgt0p6"+s, ak4jets_pt()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
              plot1d("h_jets_empt_HEHFgap_nemfgt0p6"+s, ak4jets_pt()[j] * ak4jets_NEMF()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
            } else {
              plot1d("h_jets_pt_HEHFgap_nemflt0p6"+s, ak4jets_pt()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
              plot1d("h_jets_empt_HEHFgap_nemflt0p6"+s, ak4jets_pt()[j] * ak4jets_NEMF()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
            }
            if (ak4jets_NEMF()[j] > 0.5) {
              plot1d("h_jets_pt_HEHFgap_nemfgt0p5"+s, ak4jets_pt()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
              plot1d("h_jets_empt_HEHFgap_nemfgt0p5"+s, ak4jets_pt()[j] * ak4jets_NEMF()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
            } else {
              plot1d("h_jets_pt_HEHFgap_nemflt0p5"+s, ak4jets_pt()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
              plot1d("h_jets_empt_HEHFgap_nemflt0p5"+s, ak4jets_pt()[j] * ak4jets_NEMF()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
            }
            if (ak4jets_NEMF()[j] > 0.7) {
              plot1d("h_jets_pt_HEHFgap_nemfgt0p7"+s, ak4jets_pt()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
              plot1d("h_jets_empt_HEHFgap_nemfgt0p7"+s, ak4jets_pt()[j] * ak4jets_NEMF()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
            } else {
              plot1d("h_jets_pt_HEHFgap_nemflt0p7"+s, ak4jets_pt()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
              plot1d("h_jets_empt_HEHFgap_nemflt0p7"+s, ak4jets_pt()[j] * ak4jets_NEMF()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
            }
            if (ak4jets_NEMF()[j] > nemf_thr) {
              plot1d("h_jets_pt_HEHFgap_nemfgtthr"+s, ak4jets_pt()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
              plot1d("h_jets_empt_HEHFgap_nemfgtthr"+s, ak4jets_pt()[j] * ak4jets_NEMF()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
            } else {
              plot1d("h_jets_pt_HEHFgap_nemfltthr"+s, ak4jets_pt()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
              plot1d("h_jets_empt_HEHFgap_nemfltthr"+s, ak4jets_pt()[j] * ak4jets_NEMF()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
            }
          }
          if (ak4jets_NEMF()[j] < 0.6) {
            plot1d("h_jetpt_nemflt0p6"+s, ak4jets_pt()[j], weight, hvec, ";p_{T}^{jet} [GeV]" , 80,  0, 400);
            plot1d("h_jeteta_nemflt0p6"+s, ak4jets_eta()[j], weight, hvec, ";#eta(jet)"  , 96,  -4.8f, 4.8f);
            plot1d("h_dphiVmet_nemflt0p6"+s, dPhi_pTboson_pTmiss(), weight, hvec, ";#Delta#phi(boson, E_{T}^{miss}) ", 32,  0, 3.2);
            plot1d("h_dphijmet_nemflt0p6"+s, min_abs_dPhi_pTj_pTmiss(), weight, hvec, ";min#Delta#phi(jet, E_{T}^{miss}) ", 32,  0, 3.2);
          } else {
            plot1d("h_jetpt_nemfgt0p6"+s, ak4jets_pt()[j], weight, hvec, ";p_{T}^{jet} [GeV]" , 80,  0, 400);
            plot1d("h_jeteta_nemfgt0p6"+s, ak4jets_eta()[j], weight, hvec, ";#eta(jet)"  , 96,  -4.8f, 4.8f);
            plot1d("h_dphiVmet_nemfgt0p6"+s, dPhi_pTboson_pTmiss(), weight, hvec, ";#Delta#phi(boson, E_{T}^{miss}) ", 32,  0, 3.2);
            plot1d("h_dphijmet_nemfgt0p6"+s, min_abs_dPhi_pTj_pTmiss(), weight, hvec, ";min#Delta#phi(jet, E_{T}^{miss}) ", 32,  0, 3.2);
            plot1d("h_jets_empt_nemfgt0p6"+s, ak4jets_pt()[j] * ak4jets_NEMF()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
            if (ak4jets_pt()[j] < 50) {
              plot1d("h_jets_eta_nemfgt0p6_jptlt50"+s, ak4jets_eta()[j], weight, hvec, ";#eta(jet)"  , 96,  -4.8f, 4.8f);
              plot1d("h_dphijmet_nemfgt0p6_jptlt50"+s, min_abs_dPhi_pTj_pTmiss(), weight, hvec, ";min#Delta#phi(jet, E_{T}^{miss}) ", 32,  0, 3.2);
            } else {
              plot1d("h_jets_eta_nemfgt0p6_jptgt50"+s, ak4jets_eta()[j], weight, hvec, ";#eta(jet)"  , 96,  -4.8f, 4.8f);
              plot1d("h_dphijmet_nemfgt0p6_jptgt50"+s, min_abs_dPhi_pTj_pTmiss(), weight, hvec, ";min#Delta#phi(jet, E_{T}^{miss}) ", 32,  0, 3.2);
            }
          }
          if (ak4jets_NEMF()[j] < 0.5) {
            plot1d("h_jetpt_nemflt0p5"+s, ak4jets_pt()[j], weight, hvec, ";p_{T}^{jet} [GeV]" , 80,  0, 400);
            plot1d("h_jeteta_nemflt0p5"+s, ak4jets_eta()[j], weight, hvec, ";#eta(jet)"  , 96,  -4.8f, 4.8f);
            plot1d("h_dphiVmet_nemflt0p5"+s, dPhi_pTboson_pTmiss(), weight, hvec, ";#Delta#phi(boson, E_{T}^{miss}) ", 32,  0, 3.2);
            plot1d("h_dphijmet_nemflt0p5"+s, min_abs_dPhi_pTj_pTmiss(), weight, hvec, ";min#Delta#phi(jet, E_{T}^{miss}) ", 32,  0, 3.2);
            plot1d("h_jets_empt_nemflt0p5"+s, ak4jets_pt()[j] * ak4jets_NEMF()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
          } else {
            plot1d("h_jetpt_nemfgt0p5"+s, ak4jets_pt()[j], weight, hvec, ";p_{T}^{jet} [GeV]" , 80,  0, 400);
            plot1d("h_jeteta_nemfgt0p5"+s, ak4jets_eta()[j], weight, hvec, ";#eta(jet)"  , 96,  -4.8f, 4.8f);
            plot1d("h_dphiVmet_nemfgt0p5"+s, dPhi_pTboson_pTmiss(), weight, hvec, ";#Delta#phi(boson, E_{T}^{miss}) ", 32,  0, 3.2);
            plot1d("h_dphijmet_nemfgt0p5"+s, min_abs_dPhi_pTj_pTmiss(), weight, hvec, ";min#Delta#phi(jet, E_{T}^{miss}) ", 32,  0, 3.2);
            plot1d("h_jets_empt_nemfgt0p5"+s, ak4jets_pt()[j] * ak4jets_NEMF()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
          }
          if (ak4jets_NEMF()[j] < 0.7) {
            plot1d("h_jetpt_nemflt0p7"+s, ak4jets_pt()[j], weight, hvec, ";p_{T}^{jet} [GeV]" , 80,  0, 400);
            plot1d("h_jeteta_nemflt0p7"+s, ak4jets_eta()[j], weight, hvec, ";#eta(jet)"  , 96,  -4.8f, 4.8f);
            plot1d("h_dphiVmet_nemflt0p7"+s, dPhi_pTboson_pTmiss(), weight, hvec, ";#Delta#phi(boson, E_{T}^{miss}) ", 32,  0, 3.2);
            plot1d("h_dphijmet_nemflt0p7"+s, min_abs_dPhi_pTj_pTmiss(), weight, hvec, ";min#Delta#phi(jet, E_{T}^{miss}) ", 32,  0, 3.2);
            plot1d("h_jets_empt_nemflt0p7"+s, ak4jets_pt()[j] * ak4jets_NEMF()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
          } else {
            plot1d("h_jetpt_nemfgt0p7"+s, ak4jets_pt()[j], weight, hvec, ";p_{T}^{jet} [GeV]" , 80,  0, 400);
            plot1d("h_jeteta_nemfgt0p7"+s, ak4jets_eta()[j], weight, hvec, ";#eta(jet)"  , 96,  -4.8f, 4.8f);
            plot1d("h_dphiVmet_nemfgt0p7"+s, dPhi_pTboson_pTmiss(), weight, hvec, ";#Delta#phi(boson, E_{T}^{miss}) ", 32,  0, 3.2);
            plot1d("h_dphijmet_nemfgt0p7"+s, min_abs_dPhi_pTj_pTmiss(), weight, hvec, ";min#Delta#phi(jet, E_{T}^{miss}) ", 32,  0, 3.2);
            plot1d("h_jets_empt_nemfgt0p7"+s, ak4jets_pt()[j] * ak4jets_NEMF()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
          }
          if (ak4jets_NEMF()[j] < nemf_thr) {
            plot1d("h_jetpt_nemfltthr"+s, ak4jets_pt()[j], weight, hvec, ";p_{T}^{jet} [GeV]" , 80,  0, 400);
            plot1d("h_jeteta_nemfltthr"+s, ak4jets_eta()[j], weight, hvec, ";#eta(jet)"  , 96,  -4.8f, 4.8f);
            plot1d("h_dphiVmet_nemfltthr"+s, dPhi_pTboson_pTmiss(), weight, hvec, ";#Delta#phi(boson, E_{T}^{miss}) ", 32,  0, 3.2);
            plot1d("h_dphijmet_nemfltthr"+s, min_abs_dPhi_pTj_pTmiss(), weight, hvec, ";min#Delta#phi(jet, E_{T}^{miss}) ", 32,  0, 3.2);
            plot1d("h_jets_empt_nemfltthr"+s, ak4jets_pt()[j] * ak4jets_NEMF()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
          } else {
            plot1d("h_jetpt_nemfgtthr"+s, ak4jets_pt()[j], weight, hvec, ";p_{T}^{jet} [GeV]" , 80,  0, 400);
            plot1d("h_jeteta_nemfgtthr"+s, ak4jets_eta()[j], weight, hvec, ";#eta(jet)"  , 96,  -4.8f, 4.8f);
            plot1d("h_dphiVmet_nemfgtthr"+s, dPhi_pTboson_pTmiss(), weight, hvec, ";#Delta#phi(boson, E_{T}^{miss}) ", 32,  0, 3.2);
            plot1d("h_dphijmet_nemfgtthr"+s, min_abs_dPhi_pTj_pTmiss(), weight, hvec, ";min#Delta#phi(jet, E_{T}^{miss}) ", 32,  0, 3.2);
            plot1d("h_jets_empt_nemfgtthr"+s, ak4jets_pt()[j] * ak4jets_NEMF()[j], weight, hvec, ";p_{T}^{jet1} [GeV]" , 80,  0, 400);
          }
          if (ak4jets_pt()[j] < 50) {
            plot1d("h_jets_eta_jptlt50"+s, ak4jets_eta()[j], weight, hvec, ";#eta(jet)"  , 96,  -4.8f, 4.8f);
            plot1d("h_jets_nemf_jptlt50"+s, ak4jets_NEMF()[j], weight, hvec, ";jet neutral EM fraction" , 50,  0, 1);
            plot1d("h_jets_nemfthr_jptlt50"+s, nemf_thr, weight, hvec, ";jet neutral EM fraction threshold" , 50,  0, 1);
          } else {
            plot1d("h_jets_eta_jptgt50"+s, ak4jets_eta()[j], weight, hvec, ";#eta(jet)"  , 96,  -4.8f, 4.8f);
            plot1d("h_jets_nemf_jptgt50"+s, ak4jets_NEMF()[j], weight, hvec, ";jet neutral EM fraction" , 50,  0, 1);
            plot1d("h_jets_nemfthr_jptgt50"+s, nemf_thr, weight, hvec, ";jet neutral EM fraction threshold" , 50,  0, 1);
          }
        }
        if (jet_nemf_inHEHFgap > 0) {
          plot1d("h_met_jetinHEHFgap"+s, event_pTmiss(), wgt, hvec, ";p_{T}^{miss} [GeV]" , 160, 0, 800);
          if (jet_nemf_inHEHFgap < 0.6)
            plot1d("h_met_jetinHEHFgap_nemflt0p6"+s, event_pTmiss(), wgt, hvec, ";p_{T}^{miss} [GeV]" , 160, 0, 800);
          else
            plot1d("h_met_jetinHEHFgap_nemfgt0p6"+s, event_pTmiss(), wgt, hvec, ";p_{T}^{miss} [GeV]" , 160, 0, 800);
        }
        plot1d("h_dphiVmet_jchk"+s, dPhi_pTboson_pTmiss(), weight, hvec, ";#Delta#phi(boson, E_{T}^{miss}) ", 32,  0, 3.2);
        plot1d("h_dphijmet_jchk"+s, min_abs_dPhi_pTj_pTmiss(), weight, hvec, ";min#Delta#phi(j, E_{T}^{miss}) ", 32,  0, 3.2);
      };
      // fillJetMassHists("_raw_fullMET");
      // if (event_pTmiss() > 50) fillJetMassHists("_raw_metge50");
      // if (event_pTmiss() > 125) fillJetMassHists("_raw_metge125");

      //  Chekc the sample names
      if (isZG_nlo_incl   && (isgamma || isllg) && photon_pt() > 135) continue;
      if (isZG_nlo_ptG130 && (isgamma || isllg) && photon_pt() < 135) continue;
      // if (isZllG_lo_ptG130 && isllg && photon_pt() < 135) continue;
      if (isZllG_lo && isllg && photon_pt() > 135) continue;
      if (isZG_nlo_ptG130 && isllg && year > 2016) weight *= 0.776;

      // // FIXME: debug only
      // auto checkevt = [&]() {
      //   for (int i = 1; i < mlllist.size(); ++i) {
      //     // if (fabs(dilepton_mass() - mlllist[i]) < 0.00001) return true;
      //     if (fabs(dilepton_mass() - mlllist[i]) < 0.00001 && fabs(dilepton_pt() - ptlllist[i]) < 0.0001) return i;
      //   }
      //   return 0;
      // };

      int istep = 0;
      auto fill_passedsteps = [&](string s="") {
        if (!makeStepPlots) return;
        plot1d("h_metphi_step"+to_string(istep)+s, event_phimiss(), weight, hvec, ";#phi(E_{T}^{miss})"  , 64, -3.2, 3.2);
        plot1d("h_met_step"+to_string(istep)+s, event_pTmiss(), weight, hvec, ";E_{T}^{miss} [GeV]"  , 160, 0, 800);
        plot1d("h_passed_steps"+lepcat, istep , weight, hvec, ";step" , 20,  0, 20);
        plot1d("h_passed_steps", istep , weight, hvec, ";step" , 20,  0, 20);
        istep++;
        // if (checkevt()) cout << __LINE__ << ": dilepton_mass()= " << setprecision(14) << dilepton_mass() << ", s= " << s << ", dilepton_pt()= " << dilepton_pt() << endl;
      };
      fill_passedsteps("_nocut");

      // spurious weight veto
      if (!is_data && fabs(weight) >= 1e5) {
        if (!sample.BeginsWith("QCD"))
        cout << "Suprious weight veto: sample= " << sample << ", weight= " << weight << ", event= " << evt << ", currentFile->GetTitle() = " << currentFile->GetTitle() << endl;
        continue;
      }
      
      // btag veto, already applied on SinglePhotonEvents but not others
      if (!isgamma && event_n_ak4jets_pt30_btagged_loose() != 0) continue;
      if (isgamma && event_n_leptons_fakeableBase() > 0) continue;
      if ((is1el || islg) && event_n_leptons_fakeableBase() > 1) continue;
      if ((isllg || isdilep) && event_n_leptons_fakeableBase() > 2) continue;

      // photon quality cuts
      bool veto = false;
      if (isgamma || isllg) {
        if (extendEEphoton2j) {
          if (event_n_ak4jets_pt30() < 2 && !photon_isEB()) continue;  // barrel only for njet < 2
        } else {
          if (!photon_isEB()) continue;  // barrel photon only for everything
        }
        if (!photon_is_conversionSafe()) continue;
        if (photon_full5x5_sigmaIEtaIEta() < 0.001) continue;
        if (!photon_is_PFID() || !photon_is_METSafe()) veto = true;
        if (photon_full5x5_sigmaIPhiIPhi() < 0.001) veto = true;
        if (photon_MIPTotalEnergy() > 4.9) veto = true;
        if (fabs(photon_seedTime()) > 2.0) veto = true;
        if (year == 2018 && photon_seedTime() > 1.0) veto = true;
      } else if (is1el) {
        if (extendEEphoton2j) {
          if (event_n_ak4jets_pt30() < 2 && (fabs(lepton_eta()) > 1.4442)) continue;  // barrel only for njet < 2
        } else {
          if (fabs(lepton_eta()) > 1.4442) continue;  // barrel photon only for everything
        }
        // R9 cuts to mimic the photon trigger
        if (year == 2016) {
          if (lepton_pt() < 190 && electron_full5x5_r9() < 0.9) continue;
        } else {
          if (lepton_pt() < 220 && electron_full5x5_r9() < 0.9) continue;
        }
      }
      if (!produceFilterList && veto) continue;
      fill_passedsteps("_photoncuts");

      if (doNvtxReweight && year > 2016) {
        if (event_n_vtxs_good() < 100) weight *= nvtxscale[event_n_vtxs_good()];  // only scale for data
        plot1d("h_nvtxs_rwtd", event_n_vtxs_good(), weight, hvec, ";Number of vertices", 100, 1, 101);
      }

      float met = event_pTmiss();
      float metphi = event_phimiss();
      float mZZ = event_mZZ();
      float mtZZ = event_mTZZ();
      float Vpt   = (isgamma || isllg)? photon_pt()  : (is1el)? lepton_pt()  : dilepton_pt();
      float Veta  = (isgamma || isllg)? photon_eta() : (is1el)? lepton_eta() : dilepton_eta();
      float Vphi  = (isgamma || isllg)? photon_phi() : (is1el)? lepton_phi() : dilepton_phi();
      float Vmass = (isdilep)? dilepton_mass() : 91.2;

      njet = event_n_ak4jets_pt30();
      dphijmet = min_abs_dPhi_pTj_pTmiss();
      dphiVmet = fabs(dPhi_pTboson_pTmiss());
      dphilljmet = fabs(dPhi_pTbosonjets_pTmiss());
      float dphiVjet = 4.0;

      // use while as an exitable if command
      while (isllg && doHllgTest) {
        // if (dphilljmet < 2.5) continue;
        // if (njet > 0 && dphijmet < 0.25) continue;
        // for (float Vpt_thr = 20; Vpt_thr < 55; Vpt_thr += 5) {
        for (string s : {"_nocut", "_mllZwin", "_mllgZwin"}) {
          if (s == "_mllZwin" && (fabs(dilepton_mass() - 91.2) > 15)) continue;
          if (s == "_mllgZwin" && (fabs(event_mllg() - 91.2) > 15)) continue;
          plot1d("h_met"+s, met, weight, hvec, ";E_{T}^{miss} [GeV]"  , 160,  0, 800);
          plot1d("h_mllg_Hwin"+s, event_mllg(), weight, hvec, ";M(ll#gamma) [GeV]" , 50, 100, 150);
          plot1d("h_mllg"+s, event_mllg(), weight, hvec, ";M(ll#gamma) [GeV]" , 160, 0, 800);
          plot1d("h_mll"+s, dilepton_mass(), weight, hvec, ";M(ll) [GeV]" , 160, 0, 800);
          plot1d("h_phopt"+s, Vpt, weight, hvec, ";p_{T}(#gamma) [GeV]" , 160, 0, 800);
          plot1d("h_llpt"+s, dilepton_pt(), weight, hvec, ";p_{T}(#gamma) [GeV]" , 160, 0, 800);
          plot1d("h_lep1pt"+s, leptons_pt().at(0) , weight, hvec, ";p_{T}(lep1) [GeV]", 25, 0, 500);
          plot1d("h_lep2pt"+s, leptons_pt().at(1) , weight, hvec, ";p_{T}(lep2) [GeV]", 20, 0, 400);
        }
        break;
      }

      double cval = (njet >= 2)? h_vbfcval->Eval(mZZ) : -1;
      auto getDjjVBF = [&](float num, TSpline3* hgval=nullptr, float gscale=1.) -> float {
        if (njet < 2) return -1.;
        float gsq = (hgval)? 1./(gscale * hgval->Eval(mZZ)) : 1;
        float res = num / (num +  (cval * gsq * gsq) * p_JJQCD_SIG_ghg2_1_JHUGen());
        return res;
      };
      DjjVBF = getDjjVBF(p_JJVBF_SIG_ghv1_1_JHUGen());
      DjjVBFL1 = getDjjVBF(p_JJVBF_SIG_ghv1prime2_1E4_JHUGen(), h_vbfgval_L1, 1e-4);
      DjjVBFa2 = getDjjVBF(p_JJVBF_SIG_ghv2_1_JHUGen(), h_vbfgval_a2);
      DjjVBFa2 = getDjjVBF(p_JJVBF_SIG_ghv4_1_JHUGen(), h_vbfgval_a3);
      DjjVBFL1ZGs = getDjjVBF(p_JJVBF_SIG_ghza1prime2_1E4_JHUGen(), h_vbfgval_L1ZGs, 1e-4);

      string jetcat = "jetcat";
      if (njet == 0) jetcat = "_eq0j";
      else if (njet == 1) jetcat = "_eq1j";
      else if (njet >= 2) jetcat = "_ge2j";
      // else if (njet == 2) jetcat = "_eq2j"; // alternative jet categories
      // else if (njet >= 3) jetcat = "_ge3j"; // alternative jet categories

      // Analysis selections
      if (Vpt < 55.0) continue;
      if (fabs(Vmass - 91.2) > 15) continue;  // only for dilepton events

      fill_passedsteps("_VptVmass");

      fillJetMassHists("_nodphi_fullMET"+jetcat);
      if (event_pTmiss() > 50) fillJetMassHists("_nodphi_metge50"+jetcat);
      if (event_pTmiss() > 125) fillJetMassHists("_nodphi_metge125"+jetcat);

      if (dphiVmet < 1.0) continue;
      if (dphilljmet < 2.5) continue;
      if (njet > 0 && dphijmet < 0.25) continue;
      if (tightDphiIn2j && njet >= 2 && dphijmet < 0.5) continue;

      // if (dphilljmet < 2.5 || (njet > 0 && dphijmet < 0.25)) {
      //   fillJetMassHists("_lowDphijmet_fullMET"+jetcat);
      //   if (event_pTmiss() > 50) fillJetMassHists("_lowDphijmet_metge50"+jetcat);
      //   if (event_pTmiss() > 125) fillJetMassHists("_lowDphijmet_metge125"+jetcat);
      //   continue;
      // }
      fill_passedsteps("_dphis");

      fillJetMassHists("_fullMET"+jetcat);
      if (event_pTmiss() > 50) fillJetMassHists("_metge50"+jetcat);
      if (event_pTmiss() > 125) fillJetMassHists("_metge125"+jetcat);

      // Make jet veto at horn regions
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
        // if (checkevt() && jetveto) cout << __LINE__ << ": dilepton_mass()= " << setprecision(14) << dilepton_mass() << ", dilepton_pt()= " << dilepton_pt() << endl;
        // if (jetveto) continue; // FIXME: to apply jet veto later, not now
      }

      if (!jetveto) {
        fillJetMassHists("_vetoed_fullMET"+jetcat);
        if (event_pTmiss() > 50) fillJetMassHists("_vetoed_metge50"+jetcat);
        if (event_pTmiss() > 125) fillJetMassHists("_vetoed_metge125"+jetcat);
      }

      if (isllg) {
        // make mass cut 
        if (!doHllgTest && fabs(dilepton_mass() - 91.2) > 15) continue;  // version one <-- use this
        if (doHllgTest && fabs(event_mllg() - 91.2) > 15) continue;  // version two
        // if (fabs(event_mllg - 91.2) > 15 && fabs(dilepton_mass() - 91.2) > 15) continue; // version three

        // Only MET doesn't include the dilepton
        LorentzVector met_p4(event_pTmiss(), 0, event_phimiss(), 0);
        met_p4 += LorentzVector(dilepton_pt(), dilepton_eta(), dilepton_phi(), dilepton_mass());
        met = met_p4.pt();
        metphi = met_p4.phi();
      }
      // End of cuts
      // --------------------------
      // Begin of scale factors

      if (blind && is_data && (lepton_cat | 1) == 1 && event_pTmiss() >= 125) continue;
      if (njet > 0) dphiVjet = deltaPhi(ak4jets_phi()[0], Vphi);          
      if (njet > 1) dphiVjet = min(dphiVjet, deltaPhi(ak4jets_phi()[1], Vphi));

      if (njet == 0 && dphiVmet < 2.5) jetcat = "_eq0j_lowdphi";  // shouldn't happen now

      tf_sgtoll = tf_sgtoee = tf_sgtomumu = 1;
      if (isgamma) {
        if (is_gjets && skimver < 4.02) {
          // apply the NLO/LO weight for the MC
          weight *= std::max(1., 1.716910-0.001221*photon_pt());
        }
        int icat = (photon_pt() > 440)? 39 : (photon_pt() - 50) / 10;
        if (doBosonPtReweight) {
          std::tie(tf_sgtoll, tferr_sgtoll) = getBosonPtScale(photon_pt(), photon_eta(), year, jetcat, ptbin0, ptbin1, (is_data || !DYclosureTest), extendEEphoton2j, doBosonEtaReweight);
          std::tie(tf_sgtoee, tf_sgtomumu, tferr_sgtoee, tferr_sgtomumu) = getBosonPtScale2(photon_pt(), photon_eta(), year, njet, ptbin0, ptbin1, doBosonEtaReweight);
        }
      }
      fill_passedsteps("_VptSF");

      if (produceFilterList && is_data) {
        ofile << RunNumber() << ':' << LuminosityBlock() << ':' << EventNumber() << ':' << veto << endl;
      }
      if (produceFilterList && veto) continue;

      // Get the old VBF jet category
      LorentzVector boson_p4(Vpt, Veta, Vphi, Vmass);
      vector<LorentzVector> jets;
      for (int j = 0; j < njet; ++j) {
        jets.emplace_back(ak4jets_pt()[j], ak4jets_eta()[j], ak4jets_phi()[j], ak4jets_mass()[j]);
      }
      is_vbf = (njet >= 2)? passVBFcuts(jets, boson_p4) : false;
      jet_cat = (njet >= 2)? (is_vbf? 2 : 1) : njet;

      // Photon plots again
      wgt = weight;
      fillPhotonHitMap("_fullMET");

      if (met > 125) {
        fillPhotonHitMap("_final");
        fill_passedsteps("_final");
      }

      auto fillhists = [&](string s = "", bool coreonly = false) {
        if (!produceResultPlots) return;
        // core hists
        plot1d("h_met"+s, met, weight, hvec, ";E_{T}^{miss} [GeV]"  , 160,  0, 800);
        plot1d("h_boson_pt"+s, Vpt, weight, hvec, ";p_{T}(boson) [GeV]" , 160, 0, 800);
        plot1d("h_njets"+s, njet, weight, hvec, ";N(jets)"  , 6,  0, 6);
        plot1d("h_mtZZ"+s, mtZZ, weight, hvec, ";m_{T}(ZZ) [GeV]"  , 160,  0, 800);
        if (coreonly) return;

        plot1d("h_metphi"+s, metphi, weight, hvec, ";#phi(E_{T}^{miss})"  , 68, -3.4, 3.4);
        plot1d("h_mZZ"+s, mZZ, weight, hvec, ";m(ZZ) [GeV]"  , 160,  0, 800);
        plot1d("h_diff_mZZ_mtZZ"+s, fabs(mZZ-mtZZ), weight, hvec, ";#Delta(m(ZZ),m_{T}(ZZ)) [GeV]"  , 80,  0, 160);
        plot1d("h_precdiff_mZZ_mtZZ"+s, fabs(mZZ-mtZZ)/mtZZ, weight, hvec, ";m(ZZ) [GeV]"  , 50,  0, 2);
        plot1d("h_njets20"+s, event_n_ak4jets_pt20(), weight, hvec, ";N(jets)"  , 6,  0, 6);
        plot1d("h_nvtxs"+s, event_n_vtxs_good(), weight, hvec, ";N(vtxs)"  , 80,  0, 80);

        plot1d("h_boson_mass"+s, Vmass , weight, hvec, ";M_{ll} [GeV]" , 125,  0, 500);
        plot1d("h_boson_mass_finebin"+s, Vmass , weight, hvec, ";M_{ll} [GeV]" , 160,  75, 107);
        plot1d("h_boson_eta"+s, Veta, weight, hvec, ";#eta(boson)" , 100,  -5.f, 5.f);
        plot1d("h_boson_phi"+s, Vphi, weight, hvec, ";#phi(boson)" , 64,  -3.2f, 3.2f);
        string possuf = (fabs(Veta) < 1.4442)? "_barrel" : "_endcap";
        plot1d("h_boson_phi"+possuf+s, Vphi, weight, hvec, ";#phi(boson)" , 64,  -3.2f, 3.2f);
        plot1d("h_boson_phif"+possuf+s, phiFolding(Vphi), weight, hvec, ";folded-#phi(boson)" , 51, -1.7f, 1.7f);

        plot1d("h_min_dphijmet"+s, dphijmet, weight, hvec, ";min #Delta#phi(j, E_{T}^{miss}) ", 32,  0, 3.2);
        plot1d("h_dphi_boson_met"+s, dphiVmet, weight, hvec, ";#Delta#phi(boson, E_{T}^{miss}) ", 32,  0, 3.2);
        plot1d("h_dphi_lljets_met"+s, dphilljmet, weight, hvec, ";#Delta#phi(ll+j, E_{T}^{miss}) ", 32,  0, 3.2);
        plot1d("h_min_dphiVjet"+s, dphiVjet , weight, hvec, ";min #Delta#phi(j, boson)", 32,  0, 3.2);;

        if (njet >= 1) {
          plot1d("h_jet1pt"+s, ak4jets_pt()[0], weight, hvec, ";p_{T}^{jet1} [GeV]" , 160,  0, 800);
          plot1d("h_jet1eta"+s, ak4jets_eta()[0], weight, hvec, ";#eta(jet1)"  , 96,  -4.8f, 4.8f);
        }
        if (njet >= 2) {
          plot1d("h_jet2pt"+s, ak4jets_pt()[1], weight, hvec, ";p_{T}^{jet2} [GeV]" , 160,  0, 800);
          plot1d("h_jet2eta"+s, ak4jets_eta()[1], weight, hvec, ";#eta(jet2)"  , 96,  -4.8f, 4.8f);

          // VBF categories
          plot1d("h_DjjVBF"+s, DjjVBF, weight, hvec, ";D_{jj}(VBF) ", 52,  -0.02f, 1.02f);

          LorentzVector lljets = boson_p4;
          for (auto jet : jets) lljets += jet;
          float usmetrat = met * sin(dphilljmet) / lljets.pt();
          plot1d("h_usmetrat"+s, usmetrat, weight, hvec, "; p_{T}^{miss} sin(#Delta#phi(lljets, p_{T}^{miss}))/ p_{T}^{lljets} ", 52,  0.f, 5.2f);
          if (mtZZ > 350)
            plot1d("h_usmetrat_mtge350"+s, usmetrat, weight, hvec, "; p_{T}^{miss} sin(#Delta#phi(lljets, p_{T}^{miss}))/ p_{T}^{lljets} ", 52,  0.f, 5.2f);
        }
        plot1d("h_jetcat"+s, jet_cat, weight, hvec, ";jet cat.", 3,  0, 3.f);

        // Additional binnings used to derive gamma to ll boson pT reweighting
        plot1d("h_boson_pt_b0"+s, Vpt, weight, hvec, ";p_{T}(boson) [GeV]" , ptbin0.size()-1, ptbin0.data());
        plot1d("h_boson_pt_b1"+s, Vpt, weight, hvec, ";p_{T}(boson) [GeV]" , ptbin1.size()-1, ptbin1.data());
        plot1d("h_boson_aeta"+s, fabs(Veta), weight, hvec, ";#eta(boson)" , 26,  0.f, 2.6f);

        const vector<float> mtbin1 = {0, 75, 150, 225, 300, 375, 450, 525, 600, 725, 850, 975, 1100};
        plot1d("h_mtZZ_b1"+s,  mtZZ , weight, hvec, ";M_{T}(ZZ) [GeV]" , mtbin1.size()-1, mtbin1.data());
        const vector<float> mtbin3 = {0, 75, 150, 225, 300, 375, 450, 525, 600, 725, 850, 975, 1100};
        plot1d("h_mtZZ_b3"+s,  mtZZ , weight, hvec, ";M_{T}(ZZ) [GeV]" , mtbin3.size()-1, mtbin3.data());

        plot1d("h_ht"+s, ak4jets_HT() , weight, hvec, ";HT [GeV]" , 150,  0, 750);

        // Generator quantities, MC only
        if (!is_data) {
          plot1d("h_genmet"+s, genmet_pTmiss(), weight, hvec, ";gen-E_{T}^{miss} [GeV]"  , 160,  0, 800);
          plot1d("h_genmetphi"+s, genmet_phimiss(), weight, hvec, ";#phi(gen-E_{T}^{miss})"  , 64,  -3.2f, 3.2f);
        }

        if (isllg || isgamma) {
          plot1d("h_photon_r9"+s, photon_full5x5_r9(), weight, hvec, ";R9 (5x5) " , 70,  0, 1.4);
        
          plot1d("h_photon_sieie"+s,  photon_full5x5_sigmaIEtaIEta(), weight, hvec, ";#sigma_{i#etai#eta}(#gamma)" , 80,  0, 0.02);
          plot1d("h_photon_sipip"+s,  photon_full5x5_sigmaIPhiIPhi(), weight, hvec, ";#sigma_{i#phii#phi}(#gamma)" , 80,  0, 0.02);
          plot1d("h_photon_tseed"+s,  photon_seedTime(), weight, hvec, ";t_{seed}(#gamma) [ns]" , 120,  -15, 15);
          plot1d("h_photon_Emip"+s,   photon_MIPTotalEnergy(), weight, hvec, ";E_{MIP}(#gamma) [GeV]" , 60,  0, 12);

          const vector<float> phptbin1 = {55, 82.5, 100, 135, 180, 220, 450};
          const vector<float> phptbin2 = {55, 82.5, 100, 135, 180, 190, 450};          

          float avgwgt = weight;
          int icat1 = std::upper_bound(phptbin1.begin(), phptbin1.end(), photon_pt())- phptbin1.begin() - 1;
          int icat2 = std::upper_bound(phptbin2.begin(), phptbin2.end(), photon_pt()) - phptbin2.begin() - 1;
          plot1d(Form("h_R9_photon_phptbin%d%s", icat1, s.c_str()), photon_full5x5_r9(), avgwgt, hvec, ";R9_{5x5}(#gamma)"  , 70,  0.f, 1.4f);
          plot1d(Form("h_R9_photon_phpt2bin%d%s", icat2, s.c_str()), photon_full5x5_r9(), avgwgt, hvec, ";R9_{5x5}(#gamma)"  , 70,  0.f, 1.4f);
          plot1d(Form("h_pt_photon_phptbin%d%s", icat1, s.c_str()), photon_pt(), avgwgt, hvec, ";p_{T}^{#gamma} [GeV]" , 160,  0, 800);
          plot1d("hden_pt_b1_photon_r9id90"+s,  photon_pt(), avgwgt, hvec, ";p_{T}^{#gamma} [GeV]" , phptbin1.size()-1, phptbin1.data());
          plot1d("hden_pt_b2_photon_r9id90"+s,  photon_pt(), avgwgt, hvec, ";p_{T}^{#gamma} [GeV]" , phptbin2.size()-1, phptbin2.data());
          plot1d("hden_pt_b1_photon_r9id90_5x5ver"+s,  photon_pt(), avgwgt, hvec, ";p_{T}^{#gamma} [GeV]" , phptbin1.size()-1, phptbin1.data());
          plot1d("hden_pt_b2_photon_r9id90_5x5ver"+s,  photon_pt(), avgwgt, hvec, ";p_{T}^{#gamma} [GeV]" , phptbin2.size()-1, phptbin2.data());
          if (photon_full5x5_r9() > 0.9) {
            plot1d("hnum_pt_b1_photon_r9id90"+s,  photon_pt(), avgwgt, hvec, ";p_{T}^{#gamma} [GeV]" , phptbin1.size()-1, phptbin1.data());
            plot1d("hnum_pt_b2_photon_r9id90"+s,  photon_pt(), avgwgt, hvec, ";p_{T}^{#gamma} [GeV]" , phptbin2.size()-1, phptbin2.data());
          }
        }

        // Dilepton quantities
        if (isdilep || isllg) {
          plot1d("h_nlep"+s, leptons_id().size(), weight, hvec, ";N(lep)"  , 4,  0, 4);
          plot1d("h_llid"+s, dilepton_id(), weight, hvec, ";ll ID"  , 36,  -180, 180);
          plot1d("h_ptll"+s, dilepton_pt(), weight, hvec, ";p_{T}^{ll} [GeV]"  , 160,  0, 800);
          plot1d("h_etall"+s, dilepton_eta(), weight, hvec, ";#eta(ll)"  , 64,  -3.2f, 3.2f);
          plot1d("h_mll"+s, dilepton_mass() , weight, hvec, ";M_{ll} [GeV]" , 125,  0, 500);

          float lep1pt(leptons_pt().at(0)), lep2pt(leptons_pt().at(1));
          plot1d("h_lep1pt"+s, ((lep1pt > lep2pt)? lep1pt : lep2pt)  , weight, hvec, ";p_{T}(lep1) [GeV]"  , 25,  0, 500);
          plot1d("h_lep2pt"+s, ((lep1pt > lep2pt)? lep2pt : lep1pt)  , weight, hvec, ";p_{T}(lep2) [GeV]"  , 20,  0, 400);
          plot1d("h_lep1eta"+s,  leptons_eta().at(0) , weight, hvec, ";#eta(lep1)"        , 36, -2.4, 2.4);
          plot1d("h_lep2eta"+s,  leptons_eta().at(1) , weight, hvec, ";#eta(lep2)"        , 36, -2.4, 2.4);
        }

        if (isllg) {
          plot1d("h_mllg"+s, event_mllg(), weight, hvec, ";M(ll#gamma) [GeV]" , 100, 0, 250);
          if (dilepton_id() == -121) {
            plot1d("h_el1_pt"+s, leptons_pt().at(0), weight, hvec, ";p_{T}(el1) [GeV]" , 160, 0, 800);
            plot1d("h_el2_pt"+s, leptons_pt().at(1), weight, hvec, ";p_{T}(el1) [GeV]" , 160, 0, 800);
            plot1d("h_el1_r9"+s, electrons_full5x5_r9().at(0), weight, hvec, ";R9 (5x5) " , 70,  0, 1.4);
            plot1d("h_el2_r9"+s, electrons_full5x5_r9().at(1), weight, hvec, ";R9 (5x5) " , 70,  0, 1.4);
            if (leptons_pt().at(0) > 55) {
              int ptcat1_el = std::min(std::upper_bound(ptbin1.begin(), ptbin1.end(), leptons_pt().at(0)) - ptbin1.begin() - 1, 18L);
              float scale = 1./sf_pt_phvsel_llg.at(ptcat1_el);
              plot1d("h_els_pt_b1"+s, leptons_pt().at(0), weight*scale, hvec, ";p_{T}(els) [GeV]" , ptbin1.size()-1, ptbin1.data());
              plot1d("h_els_r9"+s, electrons_full5x5_r9().at(0), weight*scale, hvec, ";R9 (5x5) " , 70,  0, 1.4);
            }
            if (leptons_pt().at(1) > 55) {
              int ptcat1_el = std::min(std::upper_bound(ptbin1.begin(), ptbin1.end(), leptons_pt().at(1)) - ptbin1.begin() - 1, 18L);
              float scale = 1./sf_pt_phvsel_llg.at(ptcat1_el);
              plot1d("h_els_pt_b1"+s, leptons_pt().at(1), weight*scale, hvec, ";p_{T}(els) [GeV]" , ptbin1.size()-1, ptbin1.data());
              plot1d("h_els_r9"+s, electrons_full5x5_r9().at(1), weight*scale, hvec, ";R9 (5x5) " , 70,  0, 1.4);
            }
          }
        }
        
        plot1d("h_EEjetveto"+s, jetveto, weight, hvec, ";N(lep)"  , 2,  0, 2);
      };

      evtwgt_sys[systype.Data()] = weight;
      if (!is_data) {
        evtwgt_sys["PDFScaleUp"] = weight * event_wgt_adjustment_PDFScaleUp();
        evtwgt_sys["PDFScaleDn"] = weight * event_wgt_adjustment_PDFScaleDn();
        evtwgt_sys["QCDScaleUp"] = weight * event_wgt_adjustment_QCDScaleUp();
        evtwgt_sys["QCDScaleDn"] = weight * event_wgt_adjustment_QCDScaleDn();
        evtwgt_sys["AsMZUp"] = weight * event_wgt_adjustment_AsMZUp();
        evtwgt_sys["AsMZDn"] = weight * event_wgt_adjustment_AsMZDn();
        evtwgt_sys["PDFReplicaUp"] = weight * event_wgt_adjustment_PDFReplicaUp();
        evtwgt_sys["PDFReplicaDn"] = weight * event_wgt_adjustment_PDFReplicaDn();
        evtwgt_sys["PhoEffUp"] = weight * event_wgt_SFs_photons_EffUp() / event_wgt_SFs_photons();
        evtwgt_sys["PhoEffDn"] = weight * event_wgt_SFs_photons_EffDn() / event_wgt_SFs_photons();
        evtwgt_sys["PUJetIdEffUp"] = weight * event_wgt_SFs_PUJetId_EffUp() / event_wgt_SFs_PUJetId();
        evtwgt_sys["PUJetIdEffDn"] = weight * event_wgt_SFs_PUJetId_EffDn() / event_wgt_SFs_PUJetId();
        evtwgt_sys["BtagSFUp"] = weight * event_wgt_SFs_btagging_EffUp() / event_wgt_SFs_btagging();
        evtwgt_sys["BtagSFDn"] = weight * event_wgt_SFs_btagging_EffDn() / event_wgt_SFs_btagging();
        evtwgt_sys["PUUp"] = weight * event_wgt_PUUp() / event_wgt();
        evtwgt_sys["PUDn"] = weight * event_wgt_PUDn() / event_wgt();
        evtwgt_sys["L1PrefiringUp"] = weight * event_wgt_L1PrefiringUp() / event_wgt();
        evtwgt_sys["L1PrefiringDn"] = weight * event_wgt_L1PrefiringDn() / event_wgt();
      }

      // Photon trigger efficiencies
      float trigeff(1.0), tefferr(0.0);
      std::tie(trigeff, tefferr) = getPhotonTrigEffs(Vpt, year);
      evtwgt_sys["PhoTriggerEffUp"] = weight * (trigeff + tefferr) / trigeff;
      evtwgt_sys["PhoTriggerEffDn"] = weight * (trigeff - tefferr) / trigeff;
      evtwgt_sys["TransferFactorUp"] = weight * (tf_sgtoll + tferr_sgtoll) / tf_sgtoll;
      evtwgt_sys["TransferFactorDn"] = weight * (tf_sgtoll - tferr_sgtoll) / tf_sgtoll;
      if (is1el) {
        evtwgt_sys["ElecToPhotonScaleUp"] = weight * (tf_etog + tferr_etog) / tf_etog;
        evtwgt_sys["ElecToPhotonScaleDn"] = weight * (tf_etog - tferr_etog) / tf_etog;
      }

      if (applyZGweightFromLLG && (isZG_nlo_incl || isZG_nlo_ptG130 || isZllG_lo) && met > 0) {
        for (string isys : systematics) {
          evtwgt_sys[isys] *= getNjetSFfromLLG(njet, year, isys);
        }
        weight = evtwgt_sys[systype.Data()];
      }

      string metsuf = "_fullMET";

      if (produceResultPlots) 
      fillhists(metsuf);
      fillhists(metsuf+jetcat);
      fillhists(metsuf+lepcat);
      fillhists(metsuf+jetcat+lepcat);

      for (string isys : systematics) {
        fillhists(metsuf+"_"+isys, true);
        fillhists(metsuf+"_"+isys+lepcat, true);
        fillhists(metsuf+"_"+isys+jetcat+lepcat, true);
      }

      if (lepcat == "_mumug" || lepcat == "_eeg") {
        fillhists(metsuf+"_llg");
        fillhists(metsuf+jetcat+"_llg");
        if (met > 50) fillhists("_metge50_llg");
        for (string isys : systematics) {
          fillhists(metsuf+"_"+isys+"_llg", true);
          fillhists(metsuf+"_"+isys+jetcat+"_llg", true);
        }
      }

      if (met > 50) {
        fillhists("_metge50"+jetcat+lepcat);
      }

      // if (met < 50) metsuf = "_metlt50";
      // else if (met < 125) metsuf = "_met50to125";
      if (met < 125) metsuf = "_metlt125";
      else if (tightMETin2j && njet >= 2 && met < 140) metsuf = "_metge125";
      else metsuf = "_final";
      fillhists(metsuf);
      fillhists(metsuf+jetcat);
      fillhists(metsuf+lepcat);
      fillhists(metsuf+jetcat+lepcat);

      for (string isys : systematics) {
        weight = evtwgt_sys.at(isys);
        fillhists(metsuf+"_"+isys, true);
        fillhists(metsuf+"_"+isys+jetcat+lepcat, true);
      }
      weight = evtwgt_sys[systype.Data()];

      if (lepcat == "_mumu" || lepcat == "_ee" || (isgamma && doBosonPtReweight)) {
        if (isgamma && doBosonPtReweight) weight *= tf_sgtoll;
        fillhists("_fullMET"+jetcat+"_ll");
        // fillhists(metsuf+"_ll");
        fillhists(metsuf+jetcat+"_ll");
        for (string isys : systematics) {
          fillhists(metsuf+"_"+isys+jetcat+"_ll", true);
        }
        if (jet_cat == 2) fillhists(metsuf+"_vbf"+"_ll");
        if (jet_cat == 1) fillhists(metsuf+"_geq1j"+"_ll");
        if (mtZZ > 350) {
          fillhists(metsuf+"_mtZZge350"+jetcat+"_ll"); 
        } else {
          fillhists(metsuf+"_mtZZlt350"+jetcat+"_ll");        
        }
      }
      // Single electron to gamma transfer factors
      else if (is1el) {
        float origwgt = weight;  // should just be 1
        std::tie(tf_etog, tferr_etog) = getElecToGammaRate(Vpt, year);

        metsuf = "_fullMET";
        weight = origwgt * tf_etog;
        fillhists(metsuf+jetcat+"_gamma");
        weight = origwgt * (tf_etog+tferr_etog);
        fillhists(metsuf+jetcat+"_etogSFUp_gamma");
        weight = origwgt * (tf_etog-tferr_etog);
        fillhists(metsuf+jetcat+"_etogSFDn_gamma");

        if (met < 125) metsuf = "_metlt125";
        else if (tightMETin2j && njet >= 2 && met < 140) metsuf = "_metge125";
        else metsuf = "_final";
        weight = origwgt * tf_etog;
        fillhists(metsuf+jetcat+"_gamma");
        weight = origwgt * (tf_etog+tferr_etog);
        fillhists(metsuf+jetcat+"_etogSFUp_gamma");
        weight = origwgt * (tf_etog-tferr_etog);
        fillhists(metsuf+jetcat+"_etogSFDn_gamma");

        if (doBosonPtReweight) {
          std::tie(tf_sgtoll, tferr_sgtoll) = getBosonPtScale(lepton_pt(), lepton_eta(), year, jetcat, ptbin0, ptbin1, (is_data || !DYclosureTest), extendEEphoton2j, doBosonEtaReweight);
          std::tie(tf_sgtoee, tf_sgtomumu, tferr_sgtoee, tferr_sgtomumu) = getBosonPtScale2(lepton_pt(), lepton_eta(), year, njet, ptbin0, ptbin1, doBosonEtaReweight);

          metsuf = "_fullMET";
          weight = origwgt * tf_etog * tf_sgtoll;
          fillhists(metsuf+jetcat+"_ll");
          weight = origwgt * (tf_etog+tferr_etog) * tf_sgtoll;
          fillhists(metsuf+jetcat+"_etogSFUp_ll");
          weight = origwgt * (tf_etog-tferr_etog) * tf_sgtoll;
          fillhists(metsuf+jetcat+"_etogSFDn_ll");

          if (met < 125) metsuf = "_metlt125";
          else if (tightMETin2j && njet >= 2 && met < 140) metsuf = "_metge125";
          else metsuf = "_final";
          weight = origwgt * tf_etog * tf_sgtoll;
          fillhists(metsuf+jetcat+"_ll");
          weight = origwgt * (tf_etog+tferr_etog) * tf_sgtoll;
          fillhists(metsuf+jetcat+"_etogSFUp_ll");
          weight = origwgt * (tf_etog-tferr_etog) * tf_sgtoll;
          fillhists(metsuf+jetcat+"_etogSFDn_ll");

        }

        weight = origwgt * tf_etog * tf_sgtoll;
      }

      float met_thr = (tightMETin2j && njet >= 2)? 140 : 125.;

      if (produceResultTree && met > met_thr) {
        // Common branches
        // evt_weight = weight; <-- already assigned
        ll_pt = Vpt;
        ll_eta = Veta;
        ll_phi = Vphi;
        ll_mass = Vmass;
        ptmiss = met;
        phimiss = metphi;

        mT = mtZZ;
        M_ZZ = mZZ;

        // njet = event_Njets(); <-- already assigned above
        // DjjVBF; <-- already assigned above
        // DjjVBFL1; <-- already assigned above
        // DjjVBFa2; <-- already assigned above
        // DjjVBFa3; <-- already assigned above
        // DjjVBFL1ZGs; <-- already assigned above

        if (usStyleTree) {
          ll_id = (lepton_cat == 3 || lepton_cat == 5)? 22 : dilepton_id();
          evt_weight = weight;

          njet_m60 = 0;
          for (int j = 0; j < std::min(njet, 2U); ++j) {
            if (ak4jets_mass().at(j) > 60) njet_m60++;
          }

          if (is_data && !isllg) {
            run = RunNumber();
            lumi = LuminosityBlock();
            event = EventNumber();
          }

          int j = 0;
          LorentzVector p4_dijet;
          for (; j < std::min(njet, 2U); ++j) {
            jet_pt[j] = ak4jets_pt().at(j);
            jet_eta[j] = ak4jets_eta().at(j);
            jet_phi[j] = ak4jets_phi().at(j);
            jet_mass[j] = ak4jets_mass().at(j);
            if (njet >= 2) p4_dijet += LorentzVector(jet_pt[j], jet_eta[j], jet_phi[j], jet_mass[j]);
          }
          for (; j < 2; ++j) {
            jet_pt[j] = -999.;
            jet_eta[j] = -999.;
            jet_phi[j] = -999.;
            jet_mass[j] = -999.;
          }

          if (njet >= 2) {
            dijet_pt = p4_dijet.pt();
            dijet_mass = p4_dijet.M();
            dijet_dEta = fabs(jet_eta[0] - jet_eta[1]);
            dijet_dPhi = deltaPhi(jet_phi[0], jet_phi[1]);
          } else {
            dijet_pt = -999.;
            dijet_mass = -999.;
            dijet_dEta = -999.;
            dijet_dPhi = -999.;
          }

          nak8jet200 = ak8jets_mass().size();
          nak8jet200_m60 = 0;
          nak8jet200_m140 = 0;
          for (int j = 0; j < nak8jet200; ++j) {
            if (ak8jets_mass()[j] >= 60 && ak8jets_mass()[j] <110) nak8jet200_m60++;
            else if (ak8jets_mass()[j] >= 140) nak8jet200_m140++;
          }
        } else {
          // lepton_cat = lepton_cat; <-- already assigned above
          // mindphi_jet_met = dphijmet; <-- already assigned above
          // dphi_boson_met = dphiVmet; <-- already assigned above
          // dphi_lljets_met = dphilljmet; <-- already assigned above
          // mZZ = M_ZZ; <-- already assigned above
          nvtxs = event_n_vtxs_good();

          lumiwgt = event_wgt();
          sf_elec = event_wgt_SFs_electrons();
          sf_muon = event_wgt_SFs_muons();
          sf_photon = event_wgt_SFs_photons();
          sf_btag = event_wgt_SFs_btagging();
          sf_pujetid = event_wgt_SFs_PUJetId();

          if (is_data) {
            run = RunNumber();
            lumi = LuminosityBlock();
            event = EventNumber();
          }

          if (isllg || isgamma) {
            ph_convveto = photon_is_conversionSafe();
            ph_passPFid = photon_is_PFID();
            ph_isEB     = photon_isEB();
            ph_pt       = photon_pt();
            ph_eta      = photon_eta();
            ph_phi      = photon_phi();
            ph_r9       = photon_full5x5_r9();
            ph_sieie    = photon_full5x5_sigmaIEtaIEta();
            ph_sipip    = photon_full5x5_sigmaIPhiIPhi();
            ph_Emip     = photon_MIPTotalEnergy();
            ph_seedtime = photon_seedTime();
            // tf_sgtoll; <-- already assigned above
          } else if (is1el) {
            ph_convveto = 1;
            ph_passPFid = 1;
            ph_isEB     = fabs(lepton_eta()) < 1.4442;
            ph_pt       = lepton_pt();
            ph_eta      = lepton_eta();
            ph_phi      = lepton_phi();
            ph_r9       = electron_full5x5_r9();
            ph_sieie    = electron_full5x5_sigmaIEtaIEta();
            ph_sipip    = electron_full5x5_sigmaIPhiIPhi();
            ph_Emip     = -1;
            ph_seedtime = electron_seedTime();
            // tf_etog; <-- already assigned above
          }

          for (int j = 0; j < std::min(njet, 32U); ++j) {
            jet_pt[j] = ak4jets_pt().at(j);
            jet_eta[j] = ak4jets_eta().at(j);
            jet_phi[j] = ak4jets_phi().at(j);
            jet_mass[j] = ak4jets_mass().at(j);
          }
        }

        for (string isys : systematics) {
          float evtwgt_ll = evtwgt_sys.at(isys);
          fouts[isys]->cd();
          if (is1el) evtwgt_ll *= tf_etog;
          if (lepton_cat == 3 || lepton_cat == 5) {
            ll_id = -121;
            evt_weight = evtwgt_ll * tf_sgtoee;
            touts[isys]->Fill();

            ll_id = -169;
            evt_weight = evtwgt_ll * tf_sgtomumu;
          }
          touts[isys]->Fill();
        }
      }

    } //event loop


    delete file;
  }//file loop

  for (const auto& h : hvec) {
    if (!produceResultPlots) continue;
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
    h.second->Write(("hnum"+hname).c_str());

    dir = (TDirectory*) fout->Get("effstudy_den");
    if (dir == nullptr) dir = fout->mkdir("effstudy_den");
    dir->cd();
    hvec.at("hden"+hname)->Write(("hden"+hname).c_str());
  }

  for (auto& h : hvec) {
    if (!produceResultPlots) continue;
    if (h.first.find("hnum") == 0 || h.first.find("hden") == 0) continue;
    if (h.first.find("phi") == string::npos && h.first.find("2d") == string::npos)
      moveOverFlowToLastBin1D(h.second);
    string dirname = "OffShell";
    vector<string> jetsufs = {"_eq0j", "_eq1j", "_eq2j", "_ge2j", "_1j20", "_vbf", "_geq1j", "_eq0j_lowdphi", };
    vector<string> lepsufs = {"_gamma", "_ee", "_mumu", "_emu", "_ll", "_llg", "_eeg", "_mumug"};
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

  if (produceFilterList) ofile.close();

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

