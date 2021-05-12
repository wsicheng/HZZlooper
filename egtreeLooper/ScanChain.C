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

#include "EGTree.h"
#include "Utilities.h"

using namespace std;
using namespace tas;

bool doBosonPtReweight = false;
const bool blind = true;

int ScanChain(TString indir, TString sample, TString tag){

  TChain *ch = new TChain("EGTree");
  TString files_in = Form("%s/%s/%s*.root", indir.Data(), tag.Data(), sample.Data());
  int a = ch->Add(files_in);
  cout << ">> Adding " << files_in << " into the chain." << endl;

  TString outname = sample;
  TString fileout = Form("output/%s_%s.root", outname.Data(), tag.Data());

  TFile* fout = new TFile(fileout, "RECREATE");
  cout << ">> Outputting to " << fileout << endl;

  if (tag.Contains("rwgtd")) {
    doBosonPtReweight = true;
    cout << ">> doBosonPtReweight is set to true." << endl;
  }

  bool is_data = (sample.BeginsWith("SingleE") || sample.BeginsWith("DoubleE") || sample.BeginsWith("EGamma"));

  int year = (tag.Contains("2016"))? 2016 : (tag.Contains("2017"))? 2017 : 2018;

  map<string,TH1*> hvec;

  unsigned nEventsTotal = 0;
  unsigned nEventsChain = ch->GetEntries();

  unsigned nEventsPositive = 0;
  unsigned nEventsNegative = 0;

  cout << ">> To run over nEventsChain = " << nEventsChain << endl;

  TFile *currentFile = 0;
  TObjArray *listOfFiles = ch->GetListOfFiles();
  TIter fileIter(listOfFiles);

  while ( (currentFile = (TFile*)fileIter.Next()) ) { 

    TString fname(currentFile->GetTitle());
    if (fname.Contains("_mc2017_") && !fname.Contains("_ext1-") ) continue;  // for only the high stat DY sample
    if (fname.Contains("Autumn18") && !fname.Contains("_ext2-") ) continue;  // for only the high stat DY sample

    TFile *file = new TFile( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("EGTree");
    st.Init(tree);

    TString filename(currentFile->GetTitle());

    for( unsigned int event = 0; event < tree->GetEntriesFast(); ++event) {

      st.GetEntry(event);
      nEventsTotal++;
      EGTree::progress(nEventsTotal, nEventsChain);

      double weight = 1;
      // double weight = event_wgt();
      if (is_data) weight = 1;
      if (weight > 0) nEventsPositive++;
      else nEventsNegative++;

      plot1d("h_num_egpair", mass_eg().size(), weight, hvec, ";N(e#gamma)" , 4,  0, 4);

      if (mass_eg().size() < 1) continue;
      if (event_pTmiss() > 125) continue;  // avoid signal region

      for (int idx = 0; idx < mass_eg().size(); ++idx) {
        // if (!isNominalTrigger()[idx]) continue;
        plot1d("h_mass_eg_raw", mass_eg()[idx], weight, hvec, ";M(e#gamma) [GeV]" , 90,  0, 180);
        if (fabs(mass_eg()[idx] - 91.2) > 15) continue;
        if (photon_full5x5_r9()[idx] < 0.9) continue;       // R9 cut for denominator
        if (!photon_is_spikeSafe()[idx]) continue;
        if (!photon_is_inTime()[idx]) continue;
        if (!photon_is_beamHaloSafe()[idx]) continue;

        // if (!photon_is_PFID()[idx]) continue;            // PFID kills photon with tracks
        // if (!photon_is_METSafe()[idx]) continue;
        // if (!photon_is_conversionSafe()[idx]) continue;
        if (photon_pt()[idx] < 55) continue;
        if (electron_pt()[idx] < 40) continue;

        const vector<float> ptbin1 = {55, 82.5, 100, 132, 181.5, 230, 450};
        const vector<float> ptbin2 = {55, 82.5, 100, 135, 200, 220, 450};

        auto fillhists = [&](int idx=0, string s="") {
          plot1d("h_pt_e"+s,  electron_pt()[idx], weight, hvec, ";p_{T}^{e} [GeV]" , 160,  0, 800);
          plot1d("h_eta_e"+s, electron_eta()[idx], weight, hvec, ";#eta(e)"  , 64,  -3.2f, 3.2f);
          plot1d("h_pt_g"+s,  photon_pt()[idx], weight, hvec, ";p_{T}^{#gamma} [GeV]" , 160,  0, 800);
          plot1d("h_eta_g"+s, photon_eta()[idx], weight, hvec, ";#eta(#gamma)"  , 64,  -3.2f, 3.2f);

          plot1d("h_pt_eg"+s,  pt_eg()[idx], weight, hvec, ";p_{T}^{e#gamma} [GeV]" , 160,  0, 800);
          plot1d("h_eta_eg"+s, eta_eg()[idx], weight, hvec, ";#eta(e#gamma)"  , 64,  -3.2f, 3.2f);
          plot1d("h_mass_eg"+s, mass_eg()[idx], weight, hvec, ";M(e#gamma) [GeV]" , 90,  0, 180);

          plot1d("h_dR_eg"+s, dR_e_g()[idx], weight, hvec, ";#eta(#gamma)"  , 64,  0.f, 3.2f);
          plot1d("h_R9_g"+s, photon_full5x5_r9()[idx], weight, hvec, ";R9_{5x5}(#gamma)"  , 70,  0.f, 1.4f);

          int icat = std::upper_bound(ptbin1.begin(), ptbin1.end(), photon_pt()[idx]) - ptbin1.begin() - 1;
          plot1d(Form("h_R9_g_phptbin%d%s", icat, s.c_str()), photon_full5x5_r9()[idx], weight, hvec, ";R9_{5x5}(#gamma)"  , 70,  0.f, 1.4f);
          plot1d(Form("h_pt_g_phptbin%d%s", icat, s.c_str()), photon_pt()[idx], weight, hvec, ";p_{T}^{#gamma} [GeV]" , 160,  0, 800);

          plot1d("hden_HLT_Photon50_R9Id90_HE10_IsoM_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
          plot1d("hden_HLT_Photon75_R9Id90_HE10_IsoM_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
          plot1d("hden_HLT_Photon90_R9Id90_HE10_IsoM_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
          plot1d("hden_HLT_Photon120_R9Id90_HE10_IsoM_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
          plot1d("hden_HLT_Photon165_R9Id90_HE10_IsoM_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);

          plot1d("hden_HLT_Photon50_R9Id90_HE10_IsoM_x_prescale_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
          plot1d("hden_HLT_Photon75_R9Id90_HE10_IsoM_x_prescale_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
          plot1d("hden_HLT_Photon90_R9Id90_HE10_IsoM_x_prescale_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
          plot1d("hden_HLT_Photon120_R9Id90_HE10_IsoM_x_prescale_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
          plot1d("hden_HLT_Photon165_R9Id90_HE10_IsoM_x_prescale_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);

          if (year == 2016) {
            plot1d("hden_HLT_Photon175_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
            plot1d("hden_HLT_Photon175_x_prescale_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
          } else {
            plot1d("hden_HLT_Photon200_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
            plot1d("hden_HLT_Photon200_x_prescale_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
          }

          if (weight_HLT_Photon50_R9Id90_HE10_IsoM()[idx] >= 1.0f) {
            plot1d("hnum_HLT_Photon50_R9Id90_HE10_IsoM_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
            plot1d("hnum_HLT_Photon50_R9Id90_HE10_IsoM_x_prescale_photon_pt"+s, photon_pt()[idx], weight*weight_HLT_Photon50_R9Id90_HE10_IsoM()[idx], hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
          }
          if (weight_HLT_Photon75_R9Id90_HE10_IsoM()[idx] >= 1.0f) {
            plot1d("hnum_HLT_Photon75_R9Id90_HE10_IsoM_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
            plot1d("hnum_HLT_Photon75_R9Id90_HE10_IsoM_x_prescale_photon_pt"+s, photon_pt()[idx], weight*weight_HLT_Photon75_R9Id90_HE10_IsoM()[idx], hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
          }
          if (weight_HLT_Photon90_R9Id90_HE10_IsoM()[idx] >= 1.0f) {
            plot1d("hnum_HLT_Photon90_R9Id90_HE10_IsoM_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
            plot1d("hnum_HLT_Photon90_R9Id90_HE10_IsoM_x_prescale_photon_pt"+s, photon_pt()[idx], weight*weight_HLT_Photon90_R9Id90_HE10_IsoM()[idx], hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
          }
          if (weight_HLT_Photon120_R9Id90_HE10_IsoM()[idx] >= 1.0f) {
            plot1d("hnum_HLT_Photon120_R9Id90_HE10_IsoM_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
            plot1d("hnum_HLT_Photon120_R9Id90_HE10_IsoM_x_prescale_photon_pt"+s, photon_pt()[idx], weight*weight_HLT_Photon120_R9Id90_HE10_IsoM()[idx], hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
          }
          if (weight_HLT_Photon165_R9Id90_HE10_IsoM()[idx] >= 1.0f) {
            plot1d("hnum_HLT_Photon165_R9Id90_HE10_IsoM_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
            plot1d("hnum_HLT_Photon165_R9Id90_HE10_IsoM_x_prescale_photon_pt"+s, photon_pt()[idx], weight*weight_HLT_Photon165_R9Id90_HE10_IsoM()[idx], hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
          }

          if (year == 2016) {
            if (weight_HLT_Photon175()[idx] >= 1.0f) {
              plot1d("hnum_HLT_Photon175_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
              plot1d("hnum_HLT_Photon175_x_prescale_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
            }
          } else {
            if (weight_HLT_Photon200()[idx] >= 1.0f) {
              plot1d("hnum_HLT_Photon200_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
              plot1d("hnum_HLT_Photon200_x_prescale_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
            } else {
              plot1d("h_HLT_Photon200_photon_failed_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 160, 0.f, 400.f);
              plot1d("h_HLT_Photon200_photon_failed_weight"+s, weight_HLT_Photon200()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", 150, 0.f, 1.5);
            }
          }
        };

        auto fillhists_b1 = [&](const vector<float>& ptbins, int idx=0, string s="_b1") {
          plot1d("hden_HLT_Photon50_R9Id90_HE10_IsoM_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
          plot1d("hden_HLT_Photon75_R9Id90_HE10_IsoM_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
          plot1d("hden_HLT_Photon90_R9Id90_HE10_IsoM_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
          plot1d("hden_HLT_Photon120_R9Id90_HE10_IsoM_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
          plot1d("hden_HLT_Photon165_R9Id90_HE10_IsoM_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());

          plot1d("hden_HLT_Photon50_R9Id90_HE10_IsoM_x_prescale_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
          plot1d("hden_HLT_Photon75_R9Id90_HE10_IsoM_x_prescale_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
          plot1d("hden_HLT_Photon90_R9Id90_HE10_IsoM_x_prescale_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
          plot1d("hden_HLT_Photon120_R9Id90_HE10_IsoM_x_prescale_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
          plot1d("hden_HLT_Photon165_R9Id90_HE10_IsoM_x_prescale_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());

          if (year == 2016) {
            plot1d("hden_HLT_Photon175_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
            plot1d("hden_HLT_Photon175_x_prescale_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
          } else {
            plot1d("hden_HLT_Photon200_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
            plot1d("hden_HLT_Photon200_x_prescale_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
          }

          if (weight_HLT_Photon50_R9Id90_HE10_IsoM()[idx] >= 1.0) {
            plot1d("hnum_HLT_Photon50_R9Id90_HE10_IsoM_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
            plot1d("hnum_HLT_Photon50_R9Id90_HE10_IsoM_x_prescale_photon_pt"+s, photon_pt()[idx], weight*weight_HLT_Photon50_R9Id90_HE10_IsoM()[idx], hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
          }
          if (weight_HLT_Photon75_R9Id90_HE10_IsoM()[idx] >= 1.0) {
            plot1d("hnum_HLT_Photon75_R9Id90_HE10_IsoM_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
            plot1d("hnum_HLT_Photon75_R9Id90_HE10_IsoM_x_prescale_photon_pt"+s, photon_pt()[idx], weight*weight_HLT_Photon75_R9Id90_HE10_IsoM()[idx], hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
          }
          if (weight_HLT_Photon90_R9Id90_HE10_IsoM()[idx] >= 1.0) {
            plot1d("hnum_HLT_Photon90_R9Id90_HE10_IsoM_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
            plot1d("hnum_HLT_Photon90_R9Id90_HE10_IsoM_x_prescale_photon_pt"+s, photon_pt()[idx], weight*weight_HLT_Photon90_R9Id90_HE10_IsoM()[idx], hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
          }
          if (weight_HLT_Photon120_R9Id90_HE10_IsoM()[idx] >= 1.0) {
            plot1d("hnum_HLT_Photon120_R9Id90_HE10_IsoM_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
            plot1d("hnum_HLT_Photon120_R9Id90_HE10_IsoM_x_prescale_photon_pt"+s, photon_pt()[idx], weight*weight_HLT_Photon120_R9Id90_HE10_IsoM()[idx], hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
          }
          if (weight_HLT_Photon165_R9Id90_HE10_IsoM()[idx] >= 1.0) {
            plot1d("hnum_HLT_Photon165_R9Id90_HE10_IsoM_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
            plot1d("hnum_HLT_Photon165_R9Id90_HE10_IsoM_x_prescale_photon_pt"+s, photon_pt()[idx], weight*weight_HLT_Photon165_R9Id90_HE10_IsoM()[idx], hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
          }

          if (year == 2016) {
            if (weight_HLT_Photon175()[idx] >= 1.0) {
              plot1d("hnum_HLT_Photon175_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
              plot1d("hnum_HLT_Photon175_x_prescale_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
            }
          } else {
            if (weight_HLT_Photon200()[idx] >= 1.0) {
              plot1d("hnum_HLT_Photon200_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
              plot1d("hnum_HLT_Photon200_x_prescale_photon_pt"+s, photon_pt()[idx], weight, hvec, "; p_{T}(#gamma) [GeV]", ptbins.size()-1, ptbins.data());
            }
          }

        };

        fillhists(idx);
        fillhists_b1(ptbin1, idx, "_b1");
        fillhists_b1(ptbin2, idx, "_b2");

        string possuf = ((photon_fid_mask()[idx] & 1) == 1)? "_barrel" : "_endcap";
        fillhists(idx, possuf);
        fillhists_b1(ptbin1, idx, "_b1"+possuf);
        fillhists_b1(ptbin2, idx, "_b2"+possuf);
      }

      // Measure electron faking photon probability
      if (mass_eg().size() == 2) {
        // if ((isNominalTrigger()[0] || isNominalTrigger()[1])) {
        if (fabs(mass_eg()[0] - 91.2) < 15 &&
            (((electron_fid_mask()[0] & 1) == 1 && electron_pt()[0] > 55) ||
             ((electron_fid_mask()[1] & 1) == 1 && electron_pt()[1] > 55) )) {
          plot1d("h_2eg_mass_eg", mass_eg()[0], weight, hvec, ";M(e#gamma) [GeV]" , 60,  60, 120);
          plot1d("h_2eg_mdiff_eg", mass_eg()[0]-mass_eg()[1], weight, hvec, ";#Delta M(e#gamma) [GeV]" , 50,  -10, 10);
          plot1d("h_2eg_e1_eta", electron_eta()[0], weight, hvec, ";#eta(#gamma)" , 48,  -2.4, 2.4);
          plot1d("h_2eg_e2_eta", electron_eta()[1], weight, hvec, ";#eta(#gamma)" , 48,  -2.4, 2.4);
          plot1d("h_2eg_pt_eg", pt_eg()[0], weight, hvec, ";M(e#gamma) [GeV]" , 80,  0, 800);
          plot1d("h_2eg_ptdiff_eg", pt_eg()[0]-pt_eg()[1], weight, hvec, ";#Delta M(e#gamma) [GeV]" , 50,  -10, 10);
          // if (fabs(mass_eg()[idx] - 91.2) > 15) continue;
          // if (full5x5_photon_r9()[idx] < 0.9) continue;
          // if (photon_pt()[idx] < 55) continue;
          // if (electron_pt()[idx] < 40) continue;
        }
      } else if (mass_eg().size() == 1) {
        string possuf = ((photon_fid_mask()[idx] & 1) == 1)? "_barrel" : "_endcap";
        bool passPhotonClean = (photon_is_spikeSafe()[0] && photon_is_inTime()[0] && photon_is_beamHaloSafe()[0]);

        if (fabs(mass_eg()[0] - 91.2) < 15) {
          plot1d("h_1eg_g_convveto"+possuf, photon_is_conversionSafe()[0], weight, hvec, ";#gamma convveto" , 2,  0, 2);
          plot1d("h_1eg_g_pt"+possuf, photon_pt()[0], weight, hvec, ";pt_{T}(#gamma)" , 80,  0, 500);
          plot1d("h_1eg_g_eta"+possuf, photon_eta()[0], weight, hvec, ";#eta(#gamma)" , 48,  -2.4, 2.4);
          plot1d("h_1eg_e_pt"+possuf, photon_is_conversionSafe()[0], weight, hvec, ";pt_{T}(e)" , 80,  0, 500);
          plot1d("h_1eg_mass_eg_allpairs"+possuf, mass_eg()[0], weight, hvec, ";M(e#gamma) [GeV]" , 60,  60, 120);
          plot1d("h_1eg_pt_eg_allpairs"+possuf, pt_eg()[0], weight, hvec, ";pt(e#gamma)" , 80,  0, 800);
          if (photon_is_conversionSafe()[0]) {
            plot1d("h_1eg_mass_eg_convveto"+possuf, mass_eg()[0], weight, hvec, ";M(e#gamma) [GeV]" , 60,  60, 120);
            plot1d("h_1eg_pt_eg_convveto"+possuf, pt_eg()[0], weight, hvec, ";pt(e#gamma) [GeV]" , 80,  0, 800);
          }        
        }
      }

    }//event loop

    delete file;
  }//file loop

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

  cout << ">> nEventsPositive= " << nEventsPositive << ", nEventsNegative= " << nEventsNegative << endl;

  fout->Close();

  return 0;
}

