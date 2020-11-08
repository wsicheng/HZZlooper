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
#include "TSpline.h"

#include "SkimTree.h"
#include "Utilities.h"

using namespace std;
using namespace tas;

bool doBosonPtReweight = false;
bool extendEEphoton2j = false;
bool doBosonEtaReweight = false;
bool doNvtxReweight = false;
bool DYclosureTest = false;

const bool makeStepPlots = false;
const bool produceResultTree = false;
const bool blind = true;

int pct = 0;

vector<float> nvtxscale(100, 0);

// enum SkimType { DiLepton, SinglePhoton, LLG, unknown};
// enum EventCat { ee, mumu, emu, gamma, llg};

int ScanChain(TString indir, TString sample, TString tag, TString specifiers = "", TString extrargs = ""){

  TChain *ch = new TChain("SkimTree");
  TString files_in = Form("%s/%s*.root", indir.Data(), sample.Data());
  ch->Add(files_in);
  cout << ">> Adding " << files_in << " into the chain." << endl;

  TString samplever = tag(0, min(tag.Index("_",tag.Index("_")+1),5));
  double skimver = TString(samplever).ReplaceAll("v","").ReplaceAll("_",".").Atof();
  cout << ">> Runing from skimver " << samplever  << ", i.e. " << skimver << endl;
  // skimver = 1.0;

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
  if (tag.Contains("flateta")) {
    doBosonEtaReweight = true;
    cout << ">> doBosonEtaReweight is set to true." << endl;
  }

  extendEEphoton2j = tag.Contains("ee2j");

  bool is_data = (sample.BeginsWith("Run201")) || sample.BeginsWith("photon") || sample.BeginsWith("egamma");
  bool is_gjets = sample.BeginsWith("GJets");

  int year = (tag.Contains("2016"))? 2016 : (tag.Contains("2017"))? 2017 : (tag.Contains("2018"))? 2018 : -1;
  if (year < 0) cout << ">> WARNING: not able to determine year from tag: " << tag << endl;

  map<string,TH1*> hvec;

  unsigned nEventsTotal = 0;
  unsigned nEventsChain = ch->GetEntries();

  // Result tree producer
  TTree* tout = nullptr;
  float evt_weight;
  float ll_pt, ll_eta, ll_phi, ll_mass;
  float ptmiss, ptmiss_phi, mT, M_ZZ;
  float dphijmet, dphiVmet, dphilljmet;
  float DjjVBF;
  int lepton_cat, jet_cat, njet, nvtxs;
  float jet1_pt, jet1_eta, jet1_phi, jet2_pt, jet2_eta, jet2_phi;
  unsigned run, lumi, event;

  bool ph_convveto;
  bool ph_passPFid;
  bool ph_isEB;

  float tf_sgtoll, tf_sgtoee, tf_sgtomumu;
  float ph_pt, ph_eta, ph_phi;
  float ph_r9, ph_sieie, ph_sipip;
  float ph_Emip, ph_seedtime, ph_e4oe1;

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
    tout->Branch("DjjVBF", &DjjVBF);
    tout->Branch("njet", &njet);
    tout->Branch("nvtxs", &nvtxs);
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

    tout->Branch("tf_sgtoll", &tf_sgtoll);
    // tout->Branch("tf_sgtoee", &tf_sgtoee);
    // tout->Branch("tf_sgtomumu", &tf_sgtomumu);

    tout->Branch("jet1_pt", &jet1_pt);
    tout->Branch("jet1_eta", &jet1_eta);
    tout->Branch("jet1_phi", &jet1_phi);
    tout->Branch("jet2_pt", &jet2_pt);
    tout->Branch("jet2_eta", &jet2_eta);
    tout->Branch("jet2_phi", &jet2_phi);
  }

  // Make photon-pt bins for reweight
  vector<float> ptbin0, ptbin1;
  std::tie(ptbin0, ptbin1) = getPtBins();
  // auto [ptbin0, ptbin1] = getPtBins();

  // For the DjjVBF calculation
  TSpline3* h_vbfcval = fetchHistCopy<TSpline3>("data/SmoothKDConstant_m4l_DjjVBF_13TeV.root", "sp_gr_varReco_Constant_Smooth");

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
  pct = 0;

  TFile *currentFile = 0;
  TObjArray *listOfFiles = ch->GetListOfFiles();
  TIter fileIter(listOfFiles);

  while ( (currentFile = (TFile*)fileIter.Next()) ) { 
    TString filename(currentFile->GetTitle());
    if (specifiers != "") {
      if (!filename.Contains(specifiers)) {
        cout << "Skipping: " << filename << endl;
        continue;
      }
    }
    bool isNLO = (filename.Contains("amcatnlo") || filename.Contains("powheg"));
    bool isZG_nlo_incl   = (filename.Contains("ZGTo2NuG_Tune") || filename.Contains("ZGTo2LG_Tune") );
    bool isZG_nlo_ptG130 = (filename.Contains("ZGTo2NuG_PtG-130") || filename.Contains("ZGTo2LG_PtG-130") );

    TFile *file = new TFile(filename);
    TTree *tree = (TTree*)file->Get("SkimTree");
    st.Init(tree);

    for (unsigned int evt = 0; evt < tree->GetEntriesFast(); ++evt) {

      st.GetEntry(evt);
      nEventsTotal++;

      SkimTree::progress(nEventsTotal, nEventsChain);

      double weight = event_wgt();
      if (is_data) weight = 1;

      bool isgamma(false), isllg(false), is1el(false);

      string lepcat = "_lepcat";
      lepton_cat = -1;
      if (indir.Contains("/DileptonEvents/")) {
        // skimtype = Dilepton;
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
      }
      if (indir.Contains("/SinglePhotonEvents/")) {
        // skimtype = SinglePhoton;
        isgamma = true;
        lepton_cat = 3;
        lepcat = "_gamma";
        weight *= event_wgt_triggers();
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
      }
      // if (indir.Contains("/SingleLeptonEvents/")) {

      // }

      if (!is_data) {
        // if (year == 2018 && event_wgt_pileup() > 100) continue;
        // weight *= event_wgt_pileup(); // pileup weight included in event_wgt for US samples
        // weight *= event_wgt_SFs(); // FIXME: fix the pileup and 
        double wgt_SFs = event_wgt_SFs_electrons() *  event_wgt_SFs_muons() * event_wgt_SFs_photons();
        wgt_SFs = std::min(wgt_SFs, 5.0);
        weight *= wgt_SFs;
      }
      // FIXME: Add block for event_wgt_adjustment_NNPDF30

      // Fill histogram of the weights
      double scale = pow(10.0, 1 - ceil(log10(fabs(weight))));
      double wgt_range = 10 * round(weight * scale) / scale;
      plot1d("h_event_wgts", weight, 1, hvec, "; event weights" , 1000,  -wgt_range, wgt_range);
      plot1d("h_event_wgts_raw", event_wgt(), 1, hvec, "; event weights" , 1000,  -wgt_range, wgt_range);
      // plot1d("h_event_wgts_pileup", event_wgt_pileup(), 1, hvec, "; event weights" , 1000,  0, 100);
      plot1d("h_event_wgts_SFs", event_wgt_SFs(), 1, hvec, "; event weights" , 100,  0, 5);
      plot1d("h_event_wgt_SFs_PUJetId",  event_wgt_SFs_PUJetId(), 1, hvec, "; event weights PUJetId" , 100,  0, 5);
      plot1d("h_event_wgt_SFs_btagging",  event_wgt_SFs_btagging(), 1, hvec, "; event weights btagging" , 100,  0, 5);
      plot1d("h_event_wgt_SFs_electrons",  event_wgt_SFs_electrons(), 1, hvec, "; event weights electrons" , 100,  0, 5);
      plot1d("h_event_wgt_SFs_muons",  event_wgt_SFs_muons(), 1, hvec, "; event weights muons" , 100,  0, 5);
      plot1d("h_event_wgt_SFs_photons",  event_wgt_SFs_photons(), 1, hvec, "; event weights photons" , 100,  0, 5);

      if (isgamma)
        plot1d("h_event_wgts_trigs", event_wgt_triggers(), 1, hvec, "; event weights" , 1000,  0, 1000);

      // plot1d("h_event_wgt_SFs_PUJetId_"+to_string(event_n_ak4jets_pt30())+"jets",  event_wgt_SFs_PUJetId(), 1, hvec, "; event weights PUJetId" , 100,  0, 5);
      // plot1d("h_event_wgt_SFs_btagging_"+to_string(event_n_ak4jets_pt30())+"jets",  event_wgt_SFs_btagging(), 1, hvec, "; event weights btagging" , 100,  0, 5);
      // if (fabs(weight) > 1000 && pct++ < 20) {
      //   cout << __LINE__ << ": weight= " << weight << ", event_wgt() = " << event_wgt() << ", event_wgt_SFs()= " << event_wgt_SFs() << endl;
      // }
      // if (fabs(event_wgt_SFs()) > 10 && pct++ < 40) {
      //   cout << __LINE__ << ": weight= " << weight << ", event_wgt() = " << event_wgt() << ", event_wgt_SFs()= " << event_wgt_SFs() << ", event_wgt_SFs_PUJetId()= " << event_wgt_SFs_PUJetId() << ", event_wgt_SFs_btagging()= " << event_wgt_SFs_btagging() << ", event_wgt_SFs_electrons()= " << event_wgt_SFs_electrons() << ", event_wgt_SFs_muons()= " << event_wgt_SFs_muons() << ", event_wgt_SFs_photons()= " << event_wgt_SFs_photons() << endl;
      // }

      double wgt = weight;
      auto fillPhotonHitMap = [&](string suf) {
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

      //  Chekc the sample names
      if (isZG_nlo_incl   && (isgamma || isllg) && photon_pt() > 135) continue;
      if (isZG_nlo_ptG130 && (isgamma || isllg) && photon_pt() < 135) continue;

      int istep = 0;
      auto fill_passedsteps = [&](string s="") {
        if (!makeStepPlots) return;
        plot1d("h_metphi_step"+to_string(istep)+s, event_phimiss(), weight, hvec, ";#phi(E_{T}^{miss})"  , 64, -3.2, 3.2);
        plot1d("h_met_step"+to_string(istep)+s, event_pTmiss(), weight, hvec, ";E_{T}^{miss} [GeV]"  , 160, 0, 800);
        plot1d("h_passed_steps_"+lepcat, istep , weight, hvec, ";step" , 20,  0, 20);
        plot1d("h_passed_steps", istep , weight, hvec, ";step" , 20,  0, 20);
        istep++;
      };
      fill_passedsteps();

      // spurious weight veto
      if (!is_data && fabs(weight) >= 1e5) {
        cout << "Suprious weight veto: sample= " << sample << ", weight= " << weight << ", event= " << event << ", currentFile->GetTitle() = " << currentFile->GetTitle() << endl;
        continue;
      }
      
      // photon quality cut
      if (isgamma || isllg) {
        if (extendEEphoton2j) {
          if (event_n_ak4jets_pt30() < 2 && !photon_isEB()) continue;  // barrel only for njet < 2
          // // Cut out beam halo spot
          // if (is_data && isgamma && !photon_isEB()) {
          //   if (fabs(photon_phi()) < 0.8 || fabs(photon_phi()) > 2.6) continue;
          //   else weight *= 3.14159 / 1.6;
          // }
        } else {
          if (!photon_isEB()) continue;  // barrel photon only for everything
        }
        if (!photon_is_conversionSafe()) continue;
        if (!photon_is_PFID()) continue;
        if (photon_MIPTotalEnergy() > 4.9) continue;
        if (photon_full5x5_sigmaIEtaIEta() < 0.001 || photon_full5x5_sigmaIPhiIPhi() < 0.001) continue;
        
        if (fabs(photon_seedTime()) > 2.0) continue;
        if (year == 2018 && photon_seedTime() > 1.0) continue;
      }
      fill_passedsteps("_photoncuts");

      if (doNvtxReweight && year > 2016) {
        if (event_n_vtxs_good() < 100) weight *= nvtxscale[event_n_vtxs_good()];  // only scale for data
        plot1d("h_nvtxs_rwtd", event_n_vtxs_good(), weight, hvec, ";Number of vertices", 100, 1, 101);
      }

      float met = event_pTmiss();
      float metphi = event_phimiss();
      float mZZ = event_mZZ();
      float mtZZ = event_mTZZ();
      float Vpt  = (isgamma || isllg)? photon_pt()  : dilepton_pt();
      float Veta = (isgamma || isllg)? photon_eta() : dilepton_eta();
      float Vphi = (isgamma || isllg)? photon_phi() : dilepton_phi();
      float Vmass = (isgamma || isllg)? 91.2 : dilepton_mass();

      njet = event_n_ak4jets_pt30();
      dphijmet = min_abs_dPhi_pTj_pTmiss();
      dphiVmet = dPhi_pTboson_pTmiss();
      dphilljmet = dPhi_pTbosonjets_pTmiss();
      float dphiVjet = 4.0;

      DjjVBF = -999.;
      if (njet >= 2) {
        double constant = h_vbfcval->Eval(mZZ);
        DjjVBF = p_JJVBF_SIG_ghv1_1_JHUGen() / (p_JJVBF_SIG_ghv1_1_JHUGen() + constant*p_JJQCD_SIG_ghg2_1_JHUGen());
      }

      // Analysis selections
      if (Vpt < 55.0) continue;
      if (fabs(Vmass - 91.2) > 15) continue;

      fill_passedsteps("_VptVmass");

      if (dphiVmet < 1.0) continue;
      if (dphilljmet < 2.5) continue;
      if (njet > 0 && dphijmet < 0.25) continue;
      fill_passedsteps("_dphis");
 
      // FIXME: check the following block
      if (isllg) {
        // make mass cut 
        if (fabs(dilepton_mass() - 91.2) > 15) continue;  // version one
        // if (fabs(event_mllg - 91.2) > 15) continue;  // version two
        // if (fabs(event_mllg - 91.2) > 15 && fabs(dilepton_mass() - 91.2) > 15) continue; // version three

        // Only MET doesn't include the dilepton
        LorentzVector met_p4(event_pTmiss(), 0, event_phimiss(), 0);
        met_p4 += LorentzVector(dilepton_pt(), dilepton_eta(), dilepton_phi(), dilepton_mass());
        met = met_p4.pt();
        metphi = met_p4.phi();
      }

      if (blind && is_data && (lepton_cat | 1) == 1 && event_pTmiss() >= 125) continue;
      if (njet > 0) dphiVjet = deltaPhi(ak4jets_phi()[0], Vphi);          
      if (njet > 1) dphiVjet = min(dphiVjet, deltaPhi(ak4jets_phi()[1], Vphi));

      string jetcat = "jetcat";
      if (njet == 0) jetcat = "_eq0j";
      else if (njet == 1) jetcat = "_eq1j";
      else if (njet == 2) jetcat = "_eq2j";
      else if (njet >= 3) jetcat = "_ge3j";

      if (njet == 0 && dphiVmet < 2.5) jetcat = "_eq0j_lowdphi";  // shouldn't happen now

      tf_sgtoll = tf_sgtoee = tf_sgtomumu = 1;
      if (isgamma) {
        if (is_gjets) {
          weight *= std::max(1., 1.716910-0.001221*photon_pt());
        }
        int icat = (photon_pt() > 440)? 39 : (photon_pt() - 50) / 10;
        if (doBosonPtReweight) {
          tf_sgtoll = getBosonPtScale(photon_pt(), photon_eta(), year, jetcat, ptbin0, ptbin1, (is_data || !DYclosureTest), extendEEphoton2j, doBosonEtaReweight);
          // weight *= scale;
        }
      }
      fill_passedsteps("_VptSF");

      bool applyZGweightFromLLG = true;
      if (applyZGweightFromLLG && (isZG_nlo_incl || isZG_nlo_ptG130) && met > 125) {
        if (year == 2018) {
          if (njet == 0) weight *= 1.07;
          if (njet == 1) weight *= 0.76;
          if (njet == 2) weight *= 1.19;
        }
      }

      // if (true) {
      //   string possuf = photon_isEB()? "_barrel" : "_endcap";
      //   wgt = fabs(weight);
      //   fillPhotonHitMap("_none"+possuf);
      //   if (!photon_is_conversionSafe()) continue;
      //   fillPhotonHitMap("_convveto"+possuf);
      //   if (!photon_is_PFID() || !photon_is_METSafe()) continue;
      //   fillPhotonHitMap("_pfid"+possuf);
      //   if (photon_MIPTotalEnergy() > 4.9) continue;
      //   fillPhotonHitMap("_mip49"+possuf);
      //   if (photon_full5x5_sigmaIEtaIEta() < 0.001 || photon_full5x5_sigmaIPhiIPhi() < 0.001) continue;
      //   fillPhotonHitMap("_sieie"+possuf);
      //   if (fabs(photon_seedTime()) > 2.0) continue;
      //   fillPhotonHitMap("_seedtime"+possuf);
      //   // if (year == 2018 && photon_seedTime() > 1.0) continue;
      //   if (photon_full5x5_r9() > 0.9)
      //     fillPhotonHitMap("_r9"+possuf);

      //   if (extendEEphoton2j) {
      //     if (event_n_ak4jets_pt30() < 2 && !photon_isEB()) continue;  // barrel only for njet < 2
      //   } else {
      //     if (!photon_isEB()) continue;  // barrel photon only for everything
      //   }
      // }

      // Photon plots again
      wgt = weight;
      fillPhotonHitMap("_fullMET");
      if (met > 125)
        fillPhotonHitMap("_final");

      auto fillhists = [&](string s = "") {
        plot1d("h_met"+s, met, weight, hvec, ";E_{T}^{miss} [GeV]"  , 160,  0, 800);
        plot1d("h_metphi"+s, metphi, weight, hvec, ";#phi(E_{T}^{miss})"  , 68, -3.4, 3.4);

        plot1d("h_mZZ"+s, mZZ, weight, hvec, ";m(ZZ) [GeV]"  , 160,  0, 800);
        plot1d("h_mtZZ"+s, mtZZ, weight, hvec, ";m_{T}(ZZ) [GeV]"  , 160,  0, 800);
        plot1d("h_diff_mZZ_mtZZ"+s, fabs(mZZ-mtZZ), weight, hvec, ";#Delta(m(ZZ),m_{T}(ZZ)) [GeV]"  , 80,  0, 160);
        plot1d("h_precdiff_mZZ_mtZZ"+s, fabs(mZZ-mtZZ)/mtZZ, weight, hvec, ";m(ZZ) [GeV]"  , 50,  0, 2);
        plot1d("h_njets"+s, njet, weight, hvec, ";N(jets)"  , 6,  0, 6);
        plot1d("h_njets20"+s, event_n_ak4jets_pt20(), weight, hvec, ";N(jets)"  , 6,  0, 6);
        plot1d("h_nvtxs"+s, event_n_vtxs_good(), weight, hvec, ";N(vtxs)"  , 80,  0, 80);

        plot1d("h_boson_mass"+s, Vmass , weight, hvec, ";M_{ll} [GeV]" , 125,  0, 500);
        plot1d("h_boson_mass_finebin"+s, Vmass , weight, hvec, ";M_{ll} [GeV]" , 160,  75, 107);
        plot1d("h_boson_pt"+s, Vpt, weight, hvec, ";p_{T}(boson) [GeV]" , 160, 0, 800);
        plot1d("h_boson_eta"+s, Veta, weight, hvec, ";#eta(boson)" , 100,  -5.f, 5.f);
        plot1d("h_boson_phi"+s, Vphi, weight, hvec, ";#phi(boson)" , 64,  -3.2f, 3.2f);
        string possuf = (fabs(Veta) < 1.4442)? "_barrel" : "_endcap";
        plot1d("h_boson_phi"+possuf+s, Vphi, weight, hvec, ";#phi(boson)" , 64,  -3.2f, 3.2f);
        plot1d("h_boson_phif"+possuf+s, phiFolding(Vphi), weight, hvec, ";folded-#phi(boson)" , 51, -1.7f, 1.7f);

        plot1d("h_min_dphijmet"+s, dphijmet, weight, hvec, ";min #Delta#phi(j, E_{T}^{miss}) ", 32,  0, 3.2);
        plot1d("h_dphi_boson_met"+s, dphiVmet, weight, hvec, ";#Delta#phi(boson, E_{T}^{miss}) ", 32,  0, 3.2);
        plot1d("h_dphi_lljets_met"+s, dphilljmet, weight, hvec, ";#Delta#phi(ll+j, E_{T}^{miss}) ", 32,  0, 3.2);
        plot1d("h_min_dphiVjet"+s, dphiVjet , weight, hvec, ";min #Delta#phi(j, boson)", 32,  0, 3.2);;
        if (njet >= 2) {
          plot1d("h_DjjVBF"+s, DjjVBF, weight, hvec, ";D_{jj}(VBF) ", 52,  -0.02f, 1.02f);
        }

        plot1d("h_boson_pt_b0"+s, Vpt, weight, hvec, ";p_{T}(boson) [GeV]" , ptbin0.size()-1, ptbin0.data());
        plot1d("h_boson_pt_b1"+s, Vpt, weight, hvec, ";p_{T}(boson) [GeV]" , ptbin1.size()-1, ptbin1.data());
        plot1d("h_boson_aeta"+s, fabs(Veta), weight, hvec, ";#eta(boson)" , 26,  0.f, 2.6f);
        
        const vector<float> mtbin1 = {0, 75, 150, 225, 300, 375, 450, 525, 600, 725, 850, 975, 1100};
        plot1d("h_mtZZ_b1"+s,  mtZZ , weight, hvec, ";M_{T}(ZZ) [GeV]" , mtbin1.size()-1, mtbin1.data());
        const vector<float> mtbin3 = {0, 75, 150, 225, 300, 375, 450, 525, 600, 725, 850, 975, 1100};
        plot1d("h_mtZZ_b3"+s,  mtZZ , weight, hvec, ";M_{T}(ZZ) [GeV]" , mtbin3.size()-1, mtbin3.data());

        plot1d("h_ht"+s, ak4jets_HT() , weight, hvec, ";HT [GeV]" , 150,  0, 750);
        if (!is_data) {
          plot1d("h_genmet"+s, genmet_pTmiss(), weight, hvec, ";gen-E_{T}^{miss} [GeV]"  , 160,  0, 800);
          plot1d("h_genmetphi"+s, genmet_phimiss(), weight, hvec, ";#phi(gen-E_{T}^{miss})"  , 64,  -3.2f, 3.2f);
        }
        if (isllg || isgamma) {
          plot1d("h_photon_r9"+s, photon_full5x5_r9(), weight, hvec, ";R9 (5x5) " , 70,  0, 1.4);
          // plot1d("h_photon_3x3_r9"+s, photon_r9(), weight, hvec, ";R9 (3x3) " , 70,  0, 1.4);
        
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
        if (!isgamma) {
          plot1d("h_nlep"+s, leptons_id().size(), weight, hvec, ";N(lep)"  , 4,  0, 4);
          plot1d("h_llid"+s, dilepton_id(), weight, hvec, ";ll ID"  , 36,  -180, 180);
          plot1d("h_ptll"+s, dilepton_pt(), weight, hvec, ";p_{T}^{ll} [GeV]"  , 160,  0, 800);
          plot1d("h_etall"+s, dilepton_eta(), weight, hvec, ";#eta(ll)"  , 64,  -3.2f, 3.2f);
          plot1d("h_mll"+s, dilepton_mass() , weight, hvec, ";M_{ll} [GeV]" , 125,  0, 500);

          plot1d("h_lep1pt"+s,   leptons_pt().at(0) , weight, hvec, ";p_{T}(lep1) [GeV]"  , 25,  0, 500);
          plot1d("h_lep2pt"+s,   leptons_pt().at(1) , weight, hvec, ";p_{T}(lep2) [GeV]"  , 20,  0, 400);
          plot1d("h_lep1eta"+s,  leptons_eta().at(0) , weight, hvec, ";#eta(lep1)"        , 36, -2.4, 2.4);
          plot1d("h_lep2eta"+s,  leptons_eta().at(1) , weight, hvec, ";#eta(lep2)"        , 36, -2.4, 2.4);
        }
        if (isllg) {
          plot1d("h_mllg", event_mllg(), weight, hvec, ";M(ll#gamma) [GeV]" , 100, 0, 250);
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

        if (njet > 0) {
          plot1d("jet1pt", ak4jets_pt()[0], weight, hvec, ";p_{T}^{jet1} [GeV]" , 160,  0, 800);
          plot1d("jet1eta", ak4jets_eta()[0], weight, hvec, ";#eta(jet1)"  , 64,  -3.2f, 3.2f);
        }
        if (njet > 1) {
          plot1d("jet2pt", ak4jets_pt()[1], weight, hvec, ";p_{T}^{jet2} [GeV]" , 160,  0, 800);
          plot1d("jet2eta", ak4jets_eta()[1], weight, hvec, ";#eta(jet2)"  , 64,  -3.2f, 3.2f);

          if (ak4jets_pt()[1] > ak4jets_pt()[0])
            cout << "Sanity check: ak4jets_pt()[1] = " << ak4jets_pt()[1]  << ", ak4jets_pt()[0]= " << ak4jets_pt()[0] << endl; // shouldn't happen
        }

        string Vptsuf = "55to83";
        if      (Vpt <  83) Vptsuf = "55to83";
        else if (Vpt < 100) Vptsuf = "83to100";
        else if (Vpt < 135) Vptsuf = "100to135";
        else if (Vpt < 180) Vptsuf = "135to180";
        else if (Vpt < 230) Vptsuf = "180to230";
        else                Vptsuf = "230toInf";
        plot1d("h_mtZZ_phpt"+Vptsuf+s, mtZZ, weight, hvec, ";m_{T}(ZZ) [GeV]"  , 160,  0, 800);

      };

      string metsuf = "_fullMET";
      fillhists(metsuf);
      fillhists(metsuf+jetcat);
      fillhists(metsuf+lepcat);
      fillhists(metsuf+jetcat+lepcat);

      // To simplify the life a little bit
      if (lepcat == "_mumu" || lepcat == "_ee") {
        fillhists(metsuf+"_ll");
        fillhists(metsuf+jetcat+"_ll");
      }
      if (lepcat == "_mumug" || lepcat == "_eeg") {
        fillhists(metsuf+"_llg");
        fillhists(metsuf+jetcat+"_llg");
      }

      // if (met < 50) metsuf = "_metlt50";
      // else if (met < 125) metsuf = "_met50to125";
      if (met < 125) metsuf = "_metlt125";
      else metsuf = "_metge125";
      fillhists(metsuf);
      fillhists(metsuf+jetcat);
      fillhists(metsuf+lepcat);
      fillhists(metsuf+jetcat+lepcat);
      if (lepcat == "_mumu" || lepcat == "_ee") {
        fillhists(metsuf+"_ll");
        fillhists(metsuf+jetcat+"_ll");
      }

      if (isgamma && doBosonPtReweight) {
        weight *= tf_sgtoll;
        fillhists("_fullMET"+jetcat+"_ll");
        fillhists(metsuf+jetcat+"_ll");
      }

      // if (125 < photon_pt() && photon_pt() < 275) {
      //   fillhists("_Vpt125to275"+metsuf+jetcat+lepcat);
      // }

      // if (is_VBFcat()) 
      //   fillhists(metsuf+"_vbf"+lepcat);

      // if (met < 125)
      //   fillhists("_metlt125"+jetcat+lepcat);
      // if (met > 50)
      //   fillhists("_metge50"+jetcat+lepcat);

      // if (met > 80)
      //   fillhists("_metge80"+jetcat+lepcat);
      // else
      //   fillhists("_metlt80"+jetcat+lepcat);
      // if (met > 80 && met < 125)
      //   fillhists("_met80to125"+jetcat+lepcat);

      // if (isgamma) {
      //   if (njet == 0 && event_n_ak4jets_pt20() == 1) {
      //     fillhists("_fullMET_1j20"+lepcat);
      //   }
      // }

      if (produceResultTree && met > 125) {
        evt_weight = weight;
        ll_pt = Vpt;
        ll_eta = Veta;
        ll_phi = Vphi;
        ll_mass = Vmass;
        ptmiss = met;
        ptmiss_phi = metphi;
        mT = mtZZ;
        M_ZZ = mZZ;
        jet_cat = njet;

        // lepton_cat = lepton_cat; <-- already assigned
        // mindphi_jet_met = dphijmet; <-- already assigned
        // dphi_boson_met = dphiVmet; <-- already assigned
        // dphi_lljets_met = dphilljmet; <-- already assigned
        // mZZ = M_ZZ; <-- already assigned
        // njet = event_Njets(); <-- already assigned
        nvtxs = event_n_vtxs_good();

        // run = event_run(); // FIXME: these are needed!
        // lumi = event_lumi(); // FIXME: these are needed!
        // event = event_event(); // FIXME: these are needed!

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
        }

        jet1_pt  = (njet > 0)? ak4jets_pt()[0]  : -999;
        jet1_eta = (njet > 0)? ak4jets_eta()[0] : -999;
        jet1_phi = (njet > 0)? ak4jets_phi()[0] : -999;
        jet2_pt  = (njet > 1)? ak4jets_pt()[1]  : -999;
        jet2_eta = (njet > 1)? ak4jets_eta()[1] : -999;
        jet2_phi = (njet > 1)? ak4jets_phi()[1] : -999;

        tout->Fill();
      }

      // // if (isgamma && met > 500 && photon_sieie() > 0.001) {
      // if (isgamma && is_data && skimver > 2.1 && met > 500 && photon_sieie() > 0.001 && event_n_ak4jets_pt30() == 0) {
      //   cout << "HighMET event: met= " << met << ", photon_pt()= " << photon_pt() << ", pt_boson= " << photon_pt() << ", photon_eta()= " << photon_eta() << ", event_run= " << event_run() << ", event_lumi()= " << event_lumi() << ", event_event()= " << event_event() << endl;
      // }
      // else if (isgamma && is_data && skimver < 2.1 && met > 500 && event_n_ak4jets_pt30() == 0) {
      //   cout << "HighMET event: met= " << met << ", photon_pt()= " << photon_pt() << ", pt_boson= " << photon_pt() << ", photon_eta()= " << photon_eta() << endl;
      // }

      // if (is_VBFcat())
      //   fillhists(metsuf+"_vbf"+lepcat);
      // else if (event_n_ak4jets_pt30() >= 1)
      //   fillhists(metsuf+"_geq1j"+lepcat);

    }//event loop

    delete file;
  }//file loop

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
    h.second->Write(("hnum"+hname).c_str());

    dir = (TDirectory*) fout->Get("effstudy_den");
    if (dir == nullptr) dir = fout->mkdir("effstudy_den");
    dir->cd();
    hvec.at("hden"+hname)->Write(("hden"+hname).c_str());
  }

  for (auto& h : hvec) {
    if (h.first.find("hnum") == 0 || h.first.find("hden") == 0) continue;
    if (h.first.find("phi") == string::npos && h.first.find("2d") == string::npos)
      moveOverFlowToLastBin1D(h.second);
    string dirname = "OffShell";
    vector<string> jetsufs = {"_eq0j", "_eq1j", "_eq2j", "_1j20", "_vbf", "_eq0j_lowdphi", };
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

  if (produceResultTree) {
    fout->cd();
    tout->Write();
  }
  fout->Close();

  return 0;
}

