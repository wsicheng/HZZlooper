#include "SkimTree.h"
SkimTree st;

void SkimTree::Init(TTree *tree) {
  tree->SetMakeClass(1);
  pt_photon_branch = 0;
  if (tree->GetBranch("pt_photon") != 0) {
    pt_photon_branch = tree->GetBranch("pt_photon");
    if (pt_photon_branch) { pt_photon_branch->SetAddress(&pt_photon_); }
  }
  chiso_photon_branch = 0;
  if (tree->GetBranch("chiso_photon") != 0) {
    chiso_photon_branch = tree->GetBranch("chiso_photon");
    if (chiso_photon_branch) { chiso_photon_branch->SetAddress(&chiso_photon_); }
  }
  id_photon_Hgg_branch = 0;
  if (tree->GetBranch("id_photon_Hgg") != 0) {
    id_photon_Hgg_branch = tree->GetBranch("id_photon_Hgg");
    if (id_photon_Hgg_branch) { id_photon_Hgg_branch->SetAddress(&id_photon_Hgg_); }
  }
  is_emu_branch = 0;
  if (tree->GetBranch("is_emu") != 0) {
    is_emu_branch = tree->GetBranch("is_emu");
    if (is_emu_branch) { is_emu_branch->SetAddress(&is_emu_); }
  }
  mipE_photon_branch = 0;
  if (tree->GetBranch("mipE_photon") != 0) {
    mipE_photon_branch = tree->GetBranch("mipE_photon");
    if (mipE_photon_branch) { mipE_photon_branch->SetAddress(&mipE_photon_); }
  }
  met_uncorr_phi_branch = 0;
  if (tree->GetBranch("met_uncorr_phi") != 0) {
    met_uncorr_phi_branch = tree->GetBranch("met_uncorr_phi");
    if (met_uncorr_phi_branch) { met_uncorr_phi_branch->SetAddress(&met_uncorr_phi_); }
  }
  mZZ_branch = 0;
  if (tree->GetBranch("mZZ") != 0) {
    mZZ_branch = tree->GetBranch("mZZ");
    if (mZZ_branch) { mZZ_branch->SetAddress(&mZZ_); }
  }
  convveto_photon_branch = 0;
  if (tree->GetBranch("convveto_photon") != 0) {
    convveto_photon_branch = tree->GetBranch("convveto_photon");
    if (convveto_photon_branch) { convveto_photon_branch->SetAddress(&convveto_photon_); }
  }
  dphi_boson_met_branch = 0;
  if (tree->GetBranch("dphi_boson_met") != 0) {
    dphi_boson_met_branch = tree->GetBranch("dphi_boson_met");
    if (dphi_boson_met_branch) { dphi_boson_met_branch->SetAddress(&dphi_boson_met_); }
  }
  dphi_lljets20_met_branch = 0;
  if (tree->GetBranch("dphi_lljets20_met") != 0) {
    dphi_lljets20_met_branch = tree->GetBranch("dphi_lljets20_met");
    if (dphi_lljets20_met_branch) { dphi_lljets20_met_branch->SetAddress(&dphi_lljets20_met_); }
  }
  event_event_branch = 0;
  if (tree->GetBranch("event_event") != 0) {
    event_event_branch = tree->GetBranch("event_event");
    if (event_event_branch) { event_event_branch->SetAddress(&event_event_); }
  }
  eta_jet2_branch = 0;
  if (tree->GetBranch("eta_jet2") != 0) {
    eta_jet2_branch = tree->GetBranch("eta_jet2");
    if (eta_jet2_branch) { eta_jet2_branch->SetAddress(&eta_jet2_); }
  }
  event_run_branch = 0;
  if (tree->GetBranch("event_run") != 0) {
    event_run_branch = tree->GetBranch("event_run");
    if (event_run_branch) { event_run_branch->SetAddress(&event_run_); }
  }
  event_Njets20_btagged_medium_branch = 0;
  if (tree->GetBranch("event_Njets20_btagged_medium") != 0) {
    event_Njets20_btagged_medium_branch = tree->GetBranch("event_Njets20_btagged_medium");
    if (event_Njets20_btagged_medium_branch) { event_Njets20_btagged_medium_branch->SetAddress(&event_Njets20_btagged_medium_); }
  }
  pt_boson_branch = 0;
  if (tree->GetBranch("pt_boson") != 0) {
    pt_boson_branch = tree->GetBranch("pt_boson");
    if (pt_boson_branch) { pt_boson_branch->SetAddress(&pt_boson_); }
  }
  id_track_branch = 0;
  if (tree->GetBranch("id_track") != 0) {
    id_track_branch = tree->GetBranch("id_track");
    if (id_track_branch) { id_track_branch->SetAddress(&id_track_); }
  }
  is_mumu_branch = 0;
  if (tree->GetBranch("is_mumu") != 0) {
    is_mumu_branch = tree->GetBranch("is_mumu");
    if (is_mumu_branch) { is_mumu_branch->SetAddress(&is_mumu_); }
  }
  event_wgt_trig_muon_branch = 0;
  if (tree->GetBranch("event_wgt_trig_muon") != 0) {
    event_wgt_trig_muon_branch = tree->GetBranch("event_wgt_trig_muon");
    if (event_wgt_trig_muon_branch) { event_wgt_trig_muon_branch->SetAddress(&event_wgt_trig_muon_); }
  }
  eta_photon_branch = 0;
  if (tree->GetBranch("eta_photon") != 0) {
    eta_photon_branch = tree->GetBranch("eta_photon");
    if (eta_photon_branch) { eta_photon_branch->SetAddress(&eta_photon_); }
  }
  is_VBFcat_branch = 0;
  if (tree->GetBranch("is_VBFcat") != 0) {
    is_VBFcat_branch = tree->GetBranch("is_VBFcat");
    if (is_VBFcat_branch) { is_VBFcat_branch->SetAddress(&is_VBFcat_); }
  }
  event_wgt_pileup_branch = 0;
  if (tree->GetBranch("event_wgt_pileup") != 0) {
    event_wgt_pileup_branch = tree->GetBranch("event_wgt_pileup");
    if (event_wgt_pileup_branch) { event_wgt_pileup_branch->SetAddress(&event_wgt_pileup_); }
  }
  event_Njets_btagged_medium_branch = 0;
  if (tree->GetBranch("event_Njets_btagged_medium") != 0) {
    event_Njets_btagged_medium_branch = tree->GetBranch("event_Njets_btagged_medium");
    if (event_Njets_btagged_medium_branch) { event_Njets_btagged_medium_branch->SetAddress(&event_Njets_btagged_medium_); }
  }
  event_Njets20_btagged_loose_branch = 0;
  if (tree->GetBranch("event_Njets20_btagged_loose") != 0) {
    event_Njets20_btagged_loose_branch = tree->GetBranch("event_Njets20_btagged_loose");
    if (event_Njets20_btagged_loose_branch) { event_Njets20_btagged_loose_branch->SetAddress(&event_Njets20_btagged_loose_); }
  }
  eta_track_branch = 0;
  if (tree->GetBranch("eta_track") != 0) {
    eta_track_branch = tree->GetBranch("eta_track");
    if (eta_track_branch) { eta_track_branch->SetAddress(&eta_track_); }
  }
  pt_track_branch = 0;
  if (tree->GetBranch("pt_track") != 0) {
    pt_track_branch = tree->GetBranch("pt_track");
    if (pt_track_branch) { pt_track_branch->SetAddress(&pt_track_); }
  }
  eta_jet1_branch = 0;
  if (tree->GetBranch("eta_jet1") != 0) {
    eta_jet1_branch = tree->GetBranch("eta_jet1");
    if (eta_jet1_branch) { eta_jet1_branch->SetAddress(&eta_jet1_); }
  }
  event_wgt_branch = 0;
  if (tree->GetBranch("event_wgt") != 0) {
    event_wgt_branch = tree->GetBranch("event_wgt");
    if (event_wgt_branch) { event_wgt_branch->SetAddress(&event_wgt_); }
  }
  pass_trackveto_branch = 0;
  if (tree->GetBranch("pass_trackveto") != 0) {
    pass_trackveto_branch = tree->GetBranch("pass_trackveto");
    if (pass_trackveto_branch) { pass_trackveto_branch->SetAddress(&pass_trackveto_); }
  }
  is_ee_branch = 0;
  if (tree->GetBranch("is_ee") != 0) {
    is_ee_branch = tree->GetBranch("is_ee");
    if (is_ee_branch) { is_ee_branch->SetAddress(&is_ee_); }
  }
  mass_boson_branch = 0;
  if (tree->GetBranch("mass_boson") != 0) {
    mass_boson_branch = tree->GetBranch("mass_boson");
    if (mass_boson_branch) { mass_boson_branch->SetAddress(&mass_boson_); }
  }
  phi_boson_branch = 0;
  if (tree->GetBranch("phi_boson") != 0) {
    phi_boson_branch = tree->GetBranch("phi_boson");
    if (phi_boson_branch) { phi_boson_branch->SetAddress(&phi_boson_); }
  }
  has_photons_inHEM1516_branch = 0;
  if (tree->GetBranch("has_photons_inHEM1516") != 0) {
    has_photons_inHEM1516_branch = tree->GetBranch("has_photons_inHEM1516");
    if (has_photons_inHEM1516_branch) { has_photons_inHEM1516_branch->SetAddress(&has_photons_inHEM1516_); }
  }
  mTZZ_branch = 0;
  if (tree->GetBranch("mTZZ") != 0) {
    mTZZ_branch = tree->GetBranch("mTZZ");
    if (mTZZ_branch) { mTZZ_branch->SetAddress(&mTZZ_); }
  }
  eta_l1_branch = 0;
  if (tree->GetBranch("eta_l1") != 0) {
    eta_l1_branch = tree->GetBranch("eta_l1");
    if (eta_l1_branch) { eta_l1_branch->SetAddress(&eta_l1_); }
  }
  has_ak4jets_inHEM1516_branch = 0;
  if (tree->GetBranch("has_ak4jets_inHEM1516") != 0) {
    has_ak4jets_inHEM1516_branch = tree->GetBranch("has_ak4jets_inHEM1516");
    if (has_ak4jets_inHEM1516_branch) { has_ak4jets_inHEM1516_branch->SetAddress(&has_ak4jets_inHEM1516_); }
  }
  id_l1_branch = 0;
  if (tree->GetBranch("id_l1") != 0) {
    id_l1_branch = tree->GetBranch("id_l1");
    if (id_l1_branch) { id_l1_branch->SetAddress(&id_l1_); }
  }
  seedtime_photon_branch = 0;
  if (tree->GetBranch("seedtime_photon") != 0) {
    seedtime_photon_branch = tree->GetBranch("seedtime_photon");
    if (seedtime_photon_branch) { seedtime_photon_branch->SetAddress(&seedtime_photon_); }
  }
  pt_jet2_branch = 0;
  if (tree->GetBranch("pt_jet2") != 0) {
    pt_jet2_branch = tree->GetBranch("pt_jet2");
    if (pt_jet2_branch) { pt_jet2_branch->SetAddress(&pt_jet2_); }
  }
  sipip_photon_branch = 0;
  if (tree->GetBranch("sipip_photon") != 0) {
    sipip_photon_branch = tree->GetBranch("sipip_photon");
    if (sipip_photon_branch) { sipip_photon_branch->SetAddress(&sipip_photon_); }
  }
  pass_lepveto_branch = 0;
  if (tree->GetBranch("pass_lepveto") != 0) {
    pass_lepveto_branch = tree->GetBranch("pass_lepveto");
    if (pass_lepveto_branch) { pass_lepveto_branch->SetAddress(&pass_lepveto_); }
  }
  phi_jet1_branch = 0;
  if (tree->GetBranch("phi_jet1") != 0) {
    phi_jet1_branch = tree->GetBranch("phi_jet1");
    if (phi_jet1_branch) { phi_jet1_branch->SetAddress(&phi_jet1_); }
  }
  sieie_photon_branch = 0;
  if (tree->GetBranch("sieie_photon") != 0) {
    sieie_photon_branch = tree->GetBranch("sieie_photon");
    if (sieie_photon_branch) { sieie_photon_branch->SetAddress(&sieie_photon_); }
  }
  is_gamma_branch = 0;
  if (tree->GetBranch("is_gamma") != 0) {
    is_gamma_branch = tree->GetBranch("is_gamma");
    if (is_gamma_branch) { is_gamma_branch->SetAddress(&is_gamma_); }
  }
  nhiso_photon_branch = 0;
  if (tree->GetBranch("nhiso_photon") != 0) {
    nhiso_photon_branch = tree->GetBranch("nhiso_photon");
    if (nhiso_photon_branch) { nhiso_photon_branch->SetAddress(&nhiso_photon_); }
  }
  isEB_photon_branch = 0;
  if (tree->GetBranch("isEB_photon") != 0) {
    isEB_photon_branch = tree->GetBranch("isEB_photon");
    if (isEB_photon_branch) { isEB_photon_branch->SetAddress(&isEB_photon_); }
  }
  has_electrons_inHEM1516_branch = 0;
  if (tree->GetBranch("has_electrons_inHEM1516") != 0) {
    has_electrons_inHEM1516_branch = tree->GetBranch("has_electrons_inHEM1516");
    if (has_electrons_inHEM1516_branch) { has_electrons_inHEM1516_branch->SetAddress(&has_electrons_inHEM1516_); }
  }
  id_l2_branch = 0;
  if (tree->GetBranch("id_l2") != 0) {
    id_l2_branch = tree->GetBranch("id_l2");
    if (id_l2_branch) { id_l2_branch->SetAddress(&id_l2_); }
  }
  event_wgt_OSSF_triggers_branch = 0;
  if (tree->GetBranch("event_wgt_OSSF_triggers") != 0) {
    event_wgt_OSSF_triggers_branch = tree->GetBranch("event_wgt_OSSF_triggers");
    if (event_wgt_OSSF_triggers_branch) { event_wgt_OSSF_triggers_branch->SetAddress(&event_wgt_OSSF_triggers_); }
  }
  eta_boson_branch = 0;
  if (tree->GetBranch("eta_boson") != 0) {
    eta_boson_branch = tree->GetBranch("eta_boson");
    if (eta_boson_branch) { eta_boson_branch->SetAddress(&eta_boson_); }
  }
  phi_photon_branch = 0;
  if (tree->GetBranch("phi_photon") != 0) {
    phi_photon_branch = tree->GetBranch("phi_photon");
    if (phi_photon_branch) { phi_photon_branch->SetAddress(&phi_photon_); }
  }
  event_DjjVBF_branch = 0;
  if (tree->GetBranch("event_DjjVBF") != 0) {
    event_DjjVBF_branch = tree->GetBranch("event_DjjVBF");
    if (event_DjjVBF_branch) { event_DjjVBF_branch->SetAddress(&event_DjjVBF_); }
  }
  event_Njets_branch = 0;
  if (tree->GetBranch("event_Njets") != 0) {
    event_Njets_branch = tree->GetBranch("event_Njets");
    if (event_Njets_branch) { event_Njets_branch->SetAddress(&event_Njets_); }
  }
  event_DjjVBF_rl_branch = 0;
  if (tree->GetBranch("event_DjjVBF_rl") != 0) {
    event_DjjVBF_rl_branch = tree->GetBranch("event_DjjVBF_rl");
    if (event_DjjVBF_rl_branch) { event_DjjVBF_rl_branch->SetAddress(&event_DjjVBF_rl_); }
  }
  pfiso_photon_branch = 0;
  if (tree->GetBranch("pfiso_photon") != 0) {
    pfiso_photon_branch = tree->GetBranch("pfiso_photon");
    if (pfiso_photon_branch) { pfiso_photon_branch->SetAddress(&pfiso_photon_); }
  }
  event_Njets20_branch = 0;
  if (tree->GetBranch("event_Njets20") != 0) {
    event_Njets20_branch = tree->GetBranch("event_Njets20");
    if (event_Njets20_branch) { event_Njets20_branch->SetAddress(&event_Njets20_); }
  }
  event_wgt_trig_electron_branch = 0;
  if (tree->GetBranch("event_wgt_trig_electron") != 0) {
    event_wgt_trig_electron_branch = tree->GetBranch("event_wgt_trig_electron");
    if (event_wgt_trig_electron_branch) { event_wgt_trig_electron_branch->SetAddress(&event_wgt_trig_electron_); }
  }
  genmet_pt_branch = 0;
  if (tree->GetBranch("genmet_pt") != 0) {
    genmet_pt_branch = tree->GetBranch("genmet_pt");
    if (genmet_pt_branch) { genmet_pt_branch->SetAddress(&genmet_pt_); }
  }
  r9_photon_branch = 0;
  if (tree->GetBranch("r9_photon") != 0) {
    r9_photon_branch = tree->GetBranch("r9_photon");
    if (r9_photon_branch) { r9_photon_branch->SetAddress(&r9_photon_); }
  }
  phi_l2_branch = 0;
  if (tree->GetBranch("phi_l2") != 0) {
    phi_l2_branch = tree->GetBranch("phi_l2");
    if (phi_l2_branch) { phi_l2_branch->SetAddress(&phi_l2_); }
  }
  event_wgt_L1prefire_branch = 0;
  if (tree->GetBranch("event_wgt_L1prefire") != 0) {
    event_wgt_L1prefire_branch = tree->GetBranch("event_wgt_L1prefire");
    if (event_wgt_L1prefire_branch) { event_wgt_L1prefire_branch->SetAddress(&event_wgt_L1prefire_); }
  }
  event_Nphotons_branch = 0;
  if (tree->GetBranch("event_Nphotons") != 0) {
    event_Nphotons_branch = tree->GetBranch("event_Nphotons");
    if (event_Nphotons_branch) { event_Nphotons_branch->SetAddress(&event_Nphotons_); }
  }
  phi_l1_branch = 0;
  if (tree->GetBranch("phi_l1") != 0) {
    phi_l1_branch = tree->GetBranch("phi_l1");
    if (phi_l1_branch) { phi_l1_branch->SetAddress(&phi_l1_); }
  }
  event_wgt_trig_photon_branch = 0;
  if (tree->GetBranch("event_wgt_trig_photon") != 0) {
    event_wgt_trig_photon_branch = tree->GetBranch("event_wgt_trig_photon");
    if (event_wgt_trig_photon_branch) { event_wgt_trig_photon_branch->SetAddress(&event_wgt_trig_photon_); }
  }
  dphi_jet20_met_branch = 0;
  if (tree->GetBranch("dphi_jet20_met") != 0) {
    dphi_jet20_met_branch = tree->GetBranch("dphi_jet20_met");
    if (dphi_jet20_met_branch) { dphi_jet20_met_branch->SetAddress(&dphi_jet20_met_); }
  }
  mindphi_jet_met_branch = 0;
  if (tree->GetBranch("mindphi_jet_met") != 0) {
    mindphi_jet_met_branch = tree->GetBranch("mindphi_jet_met");
    if (mindphi_jet_met_branch) { mindphi_jet_met_branch->SetAddress(&mindphi_jet_met_); }
  }
  pt_l2_branch = 0;
  if (tree->GetBranch("pt_l2") != 0) {
    pt_l2_branch = tree->GetBranch("pt_l2");
    if (pt_l2_branch) { pt_l2_branch->SetAddress(&pt_l2_); }
  }
  genmet_phi_branch = 0;
  if (tree->GetBranch("genmet_phi") != 0) {
    genmet_phi_branch = tree->GetBranch("genmet_phi");
    if (genmet_phi_branch) { genmet_phi_branch->SetAddress(&genmet_phi_); }
  }
  pt_l1_branch = 0;
  if (tree->GetBranch("pt_l1") != 0) {
    pt_l1_branch = tree->GetBranch("pt_l1");
    if (pt_l1_branch) { pt_l1_branch->SetAddress(&pt_l1_); }
  }
  event_nvtxs_good_branch = 0;
  if (tree->GetBranch("event_nvtxs_good") != 0) {
    event_nvtxs_good_branch = tree->GetBranch("event_nvtxs_good");
    if (event_nvtxs_good_branch) { event_nvtxs_good_branch->SetAddress(&event_nvtxs_good_); }
  }
  e4oe1_photon_branch = 0;
  if (tree->GetBranch("e4oe1_photon") != 0) {
    e4oe1_photon_branch = tree->GetBranch("e4oe1_photon");
    if (e4oe1_photon_branch) { e4oe1_photon_branch->SetAddress(&e4oe1_photon_); }
  }
  phi_track_branch = 0;
  if (tree->GetBranch("phi_track") != 0) {
    phi_track_branch = tree->GetBranch("phi_track");
    if (phi_track_branch) { phi_track_branch->SetAddress(&phi_track_); }
  }
  dphi_lljets_met_branch = 0;
  if (tree->GetBranch("dphi_lljets_met") != 0) {
    dphi_lljets_met_branch = tree->GetBranch("dphi_lljets_met");
    if (dphi_lljets_met_branch) { dphi_lljets_met_branch->SetAddress(&dphi_lljets_met_); }
  }
  pTmiss_branch = 0;
  if (tree->GetBranch("pTmiss") != 0) {
    pTmiss_branch = tree->GetBranch("pTmiss");
    if (pTmiss_branch) { pTmiss_branch->SetAddress(&pTmiss_); }
  }
  met_uncorr_pt_branch = 0;
  if (tree->GetBranch("met_uncorr_pt") != 0) {
    met_uncorr_pt_branch = tree->GetBranch("met_uncorr_pt");
    if (met_uncorr_pt_branch) { met_uncorr_pt_branch->SetAddress(&met_uncorr_pt_); }
  }
  pt_jet1_branch = 0;
  if (tree->GetBranch("pt_jet1") != 0) {
    pt_jet1_branch = tree->GetBranch("pt_jet1");
    if (pt_jet1_branch) { pt_jet1_branch->SetAddress(&pt_jet1_); }
  }
  event_lumi_branch = 0;
  if (tree->GetBranch("event_lumi") != 0) {
    event_lumi_branch = tree->GetBranch("event_lumi");
    if (event_lumi_branch) { event_lumi_branch->SetAddress(&event_lumi_); }
  }
  eta_l2_branch = 0;
  if (tree->GetBranch("eta_l2") != 0) {
    eta_l2_branch = tree->GetBranch("eta_l2");
    if (eta_l2_branch) { eta_l2_branch->SetAddress(&eta_l2_); }
  }
  phi_jet2_branch = 0;
  if (tree->GetBranch("phi_jet2") != 0) {
    phi_jet2_branch = tree->GetBranch("phi_jet2");
    if (phi_jet2_branch) { phi_jet2_branch->SetAddress(&phi_jet2_); }
  }
  event_wgt_OSDF_triggers_branch = 0;
  if (tree->GetBranch("event_wgt_OSDF_triggers") != 0) {
    event_wgt_OSDF_triggers_branch = tree->GetBranch("event_wgt_OSDF_triggers");
    if (event_wgt_OSDF_triggers_branch) { event_wgt_OSDF_triggers_branch->SetAddress(&event_wgt_OSDF_triggers_); }
  }
  phimiss_branch = 0;
  if (tree->GetBranch("phimiss") != 0) {
    phimiss_branch = tree->GetBranch("phimiss");
    if (phimiss_branch) { phimiss_branch->SetAddress(&phimiss_); }
  }
  event_HT_branch = 0;
  if (tree->GetBranch("event_HT") != 0) {
    event_HT_branch = tree->GetBranch("event_HT");
    if (event_HT_branch) { event_HT_branch->SetAddress(&event_HT_); }
  }
  event_wgt_SFs_branch = 0;
  if (tree->GetBranch("event_wgt_SFs") != 0) {
    event_wgt_SFs_branch = tree->GetBranch("event_wgt_SFs");
    if (event_wgt_SFs_branch) { event_wgt_SFs_branch->SetAddress(&event_wgt_SFs_); }
  }
  event_Njets_btagged_branch = 0;
  if (tree->GetBranch("event_Njets_btagged") != 0) {
    event_Njets_btagged_branch = tree->GetBranch("event_Njets_btagged");
    if (event_Njets_btagged_branch) { event_Njets_btagged_branch->SetAddress(&event_Njets_btagged_); }
  }
  tree->SetMakeClass(0);
}
void SkimTree::GetEntry(unsigned int idx) {
  index = idx;
  pt_photon_isLoaded = false;
  chiso_photon_isLoaded = false;
  id_photon_Hgg_isLoaded = false;
  is_emu_isLoaded = false;
  mipE_photon_isLoaded = false;
  met_uncorr_phi_isLoaded = false;
  mZZ_isLoaded = false;
  convveto_photon_isLoaded = false;
  dphi_boson_met_isLoaded = false;
  dphi_lljets20_met_isLoaded = false;
  event_event_isLoaded = false;
  eta_jet2_isLoaded = false;
  event_run_isLoaded = false;
  event_Njets20_btagged_medium_isLoaded = false;
  pt_boson_isLoaded = false;
  id_track_isLoaded = false;
  is_mumu_isLoaded = false;
  event_wgt_trig_muon_isLoaded = false;
  eta_photon_isLoaded = false;
  is_VBFcat_isLoaded = false;
  event_wgt_pileup_isLoaded = false;
  event_Njets_btagged_medium_isLoaded = false;
  event_Njets20_btagged_loose_isLoaded = false;
  eta_track_isLoaded = false;
  pt_track_isLoaded = false;
  eta_jet1_isLoaded = false;
  event_wgt_isLoaded = false;
  pass_trackveto_isLoaded = false;
  is_ee_isLoaded = false;
  mass_boson_isLoaded = false;
  phi_boson_isLoaded = false;
  has_photons_inHEM1516_isLoaded = false;
  mTZZ_isLoaded = false;
  eta_l1_isLoaded = false;
  has_ak4jets_inHEM1516_isLoaded = false;
  id_l1_isLoaded = false;
  seedtime_photon_isLoaded = false;
  pt_jet2_isLoaded = false;
  sipip_photon_isLoaded = false;
  pass_lepveto_isLoaded = false;
  phi_jet1_isLoaded = false;
  sieie_photon_isLoaded = false;
  is_gamma_isLoaded = false;
  nhiso_photon_isLoaded = false;
  isEB_photon_isLoaded = false;
  has_electrons_inHEM1516_isLoaded = false;
  id_l2_isLoaded = false;
  event_wgt_OSSF_triggers_isLoaded = false;
  eta_boson_isLoaded = false;
  phi_photon_isLoaded = false;
  event_DjjVBF_isLoaded = false;
  event_Njets_isLoaded = false;
  event_DjjVBF_rl_isLoaded = false;
  pfiso_photon_isLoaded = false;
  event_Njets20_isLoaded = false;
  event_wgt_trig_electron_isLoaded = false;
  genmet_pt_isLoaded = false;
  r9_photon_isLoaded = false;
  phi_l2_isLoaded = false;
  event_wgt_L1prefire_isLoaded = false;
  event_Nphotons_isLoaded = false;
  phi_l1_isLoaded = false;
  event_wgt_trig_photon_isLoaded = false;
  dphi_jet20_met_isLoaded = false;
  mindphi_jet_met_isLoaded = false;
  pt_l2_isLoaded = false;
  genmet_phi_isLoaded = false;
  pt_l1_isLoaded = false;
  event_nvtxs_good_isLoaded = false;
  e4oe1_photon_isLoaded = false;
  phi_track_isLoaded = false;
  dphi_lljets_met_isLoaded = false;
  pTmiss_isLoaded = false;
  met_uncorr_pt_isLoaded = false;
  pt_jet1_isLoaded = false;
  event_lumi_isLoaded = false;
  eta_l2_isLoaded = false;
  phi_jet2_isLoaded = false;
  event_wgt_OSDF_triggers_isLoaded = false;
  phimiss_isLoaded = false;
  event_HT_isLoaded = false;
  event_wgt_SFs_isLoaded = false;
  event_Njets_btagged_isLoaded = false;
}
void SkimTree::LoadAllBranches() {
  if (pt_photon_branch != 0) pt_photon();
  if (chiso_photon_branch != 0) chiso_photon();
  if (id_photon_Hgg_branch != 0) id_photon_Hgg();
  if (is_emu_branch != 0) is_emu();
  if (mipE_photon_branch != 0) mipE_photon();
  if (met_uncorr_phi_branch != 0) met_uncorr_phi();
  if (mZZ_branch != 0) mZZ();
  if (convveto_photon_branch != 0) convveto_photon();
  if (dphi_boson_met_branch != 0) dphi_boson_met();
  if (dphi_lljets20_met_branch != 0) dphi_lljets20_met();
  if (event_event_branch != 0) event_event();
  if (eta_jet2_branch != 0) eta_jet2();
  if (event_run_branch != 0) event_run();
  if (event_Njets20_btagged_medium_branch != 0) event_Njets20_btagged_medium();
  if (pt_boson_branch != 0) pt_boson();
  if (id_track_branch != 0) id_track();
  if (is_mumu_branch != 0) is_mumu();
  if (event_wgt_trig_muon_branch != 0) event_wgt_trig_muon();
  if (eta_photon_branch != 0) eta_photon();
  if (is_VBFcat_branch != 0) is_VBFcat();
  if (event_wgt_pileup_branch != 0) event_wgt_pileup();
  if (event_Njets_btagged_medium_branch != 0) event_Njets_btagged_medium();
  if (event_Njets20_btagged_loose_branch != 0) event_Njets20_btagged_loose();
  if (eta_track_branch != 0) eta_track();
  if (pt_track_branch != 0) pt_track();
  if (eta_jet1_branch != 0) eta_jet1();
  if (event_wgt_branch != 0) event_wgt();
  if (pass_trackveto_branch != 0) pass_trackveto();
  if (is_ee_branch != 0) is_ee();
  if (mass_boson_branch != 0) mass_boson();
  if (phi_boson_branch != 0) phi_boson();
  if (has_photons_inHEM1516_branch != 0) has_photons_inHEM1516();
  if (mTZZ_branch != 0) mTZZ();
  if (eta_l1_branch != 0) eta_l1();
  if (has_ak4jets_inHEM1516_branch != 0) has_ak4jets_inHEM1516();
  if (id_l1_branch != 0) id_l1();
  if (seedtime_photon_branch != 0) seedtime_photon();
  if (pt_jet2_branch != 0) pt_jet2();
  if (sipip_photon_branch != 0) sipip_photon();
  if (pass_lepveto_branch != 0) pass_lepveto();
  if (phi_jet1_branch != 0) phi_jet1();
  if (sieie_photon_branch != 0) sieie_photon();
  if (is_gamma_branch != 0) is_gamma();
  if (nhiso_photon_branch != 0) nhiso_photon();
  if (isEB_photon_branch != 0) isEB_photon();
  if (has_electrons_inHEM1516_branch != 0) has_electrons_inHEM1516();
  if (id_l2_branch != 0) id_l2();
  if (event_wgt_OSSF_triggers_branch != 0) event_wgt_OSSF_triggers();
  if (eta_boson_branch != 0) eta_boson();
  if (phi_photon_branch != 0) phi_photon();
  if (event_DjjVBF_branch != 0) event_DjjVBF();
  if (event_Njets_branch != 0) event_Njets();
  if (event_DjjVBF_rl_branch != 0) event_DjjVBF_rl();
  if (pfiso_photon_branch != 0) pfiso_photon();
  if (event_Njets20_branch != 0) event_Njets20();
  if (event_wgt_trig_electron_branch != 0) event_wgt_trig_electron();
  if (genmet_pt_branch != 0) genmet_pt();
  if (r9_photon_branch != 0) r9_photon();
  if (phi_l2_branch != 0) phi_l2();
  if (event_wgt_L1prefire_branch != 0) event_wgt_L1prefire();
  if (event_Nphotons_branch != 0) event_Nphotons();
  if (phi_l1_branch != 0) phi_l1();
  if (event_wgt_trig_photon_branch != 0) event_wgt_trig_photon();
  if (dphi_jet20_met_branch != 0) dphi_jet20_met();
  if (mindphi_jet_met_branch != 0) mindphi_jet_met();
  if (pt_l2_branch != 0) pt_l2();
  if (genmet_phi_branch != 0) genmet_phi();
  if (pt_l1_branch != 0) pt_l1();
  if (event_nvtxs_good_branch != 0) event_nvtxs_good();
  if (e4oe1_photon_branch != 0) e4oe1_photon();
  if (phi_track_branch != 0) phi_track();
  if (dphi_lljets_met_branch != 0) dphi_lljets_met();
  if (pTmiss_branch != 0) pTmiss();
  if (met_uncorr_pt_branch != 0) met_uncorr_pt();
  if (pt_jet1_branch != 0) pt_jet1();
  if (event_lumi_branch != 0) event_lumi();
  if (eta_l2_branch != 0) eta_l2();
  if (phi_jet2_branch != 0) phi_jet2();
  if (event_wgt_OSDF_triggers_branch != 0) event_wgt_OSDF_triggers();
  if (phimiss_branch != 0) phimiss();
  if (event_HT_branch != 0) event_HT();
  if (event_wgt_SFs_branch != 0) event_wgt_SFs();
  if (event_Njets_btagged_branch != 0) event_Njets_btagged();
}
const float &SkimTree::pt_photon() {
  if (not pt_photon_isLoaded) {
    if (pt_photon_branch != 0) {
      pt_photon_branch->GetEntry(index);
    } else {
      printf("branch pt_photon_branch does not exist!\n");
      exit(1);
    }
    pt_photon_isLoaded = true;
  }
  return pt_photon_;
}
const float &SkimTree::chiso_photon() {
  if (not chiso_photon_isLoaded) {
    if (chiso_photon_branch != 0) {
      chiso_photon_branch->GetEntry(index);
    } else {
      printf("branch chiso_photon_branch does not exist!\n");
      exit(1);
    }
    chiso_photon_isLoaded = true;
  }
  return chiso_photon_;
}
const bool &SkimTree::id_photon_Hgg() {
  if (not id_photon_Hgg_isLoaded) {
    if (id_photon_Hgg_branch != 0) {
      id_photon_Hgg_branch->GetEntry(index);
    } else {
      printf("branch id_photon_Hgg_branch does not exist!\n");
      exit(1);
    }
    id_photon_Hgg_isLoaded = true;
  }
  return id_photon_Hgg_;
}
const bool &SkimTree::is_emu() {
  if (not is_emu_isLoaded) {
    if (is_emu_branch != 0) {
      is_emu_branch->GetEntry(index);
    } else {
      printf("branch is_emu_branch does not exist!\n");
      exit(1);
    }
    is_emu_isLoaded = true;
  }
  return is_emu_;
}
const float &SkimTree::mipE_photon() {
  if (not mipE_photon_isLoaded) {
    if (mipE_photon_branch != 0) {
      mipE_photon_branch->GetEntry(index);
    } else {
      printf("branch mipE_photon_branch does not exist!\n");
      exit(1);
    }
    mipE_photon_isLoaded = true;
  }
  return mipE_photon_;
}
const float &SkimTree::met_uncorr_phi() {
  if (not met_uncorr_phi_isLoaded) {
    if (met_uncorr_phi_branch != 0) {
      met_uncorr_phi_branch->GetEntry(index);
    } else {
      printf("branch met_uncorr_phi_branch does not exist!\n");
      exit(1);
    }
    met_uncorr_phi_isLoaded = true;
  }
  return met_uncorr_phi_;
}
const float &SkimTree::mZZ() {
  if (not mZZ_isLoaded) {
    if (mZZ_branch != 0) {
      mZZ_branch->GetEntry(index);
    } else {
      printf("branch mZZ_branch does not exist!\n");
      exit(1);
    }
    mZZ_isLoaded = true;
  }
  return mZZ_;
}
const bool &SkimTree::convveto_photon() {
  if (not convveto_photon_isLoaded) {
    if (convveto_photon_branch != 0) {
      convveto_photon_branch->GetEntry(index);
    } else {
      printf("branch convveto_photon_branch does not exist!\n");
      exit(1);
    }
    convveto_photon_isLoaded = true;
  }
  return convveto_photon_;
}
const float &SkimTree::dphi_boson_met() {
  if (not dphi_boson_met_isLoaded) {
    if (dphi_boson_met_branch != 0) {
      dphi_boson_met_branch->GetEntry(index);
    } else {
      printf("branch dphi_boson_met_branch does not exist!\n");
      exit(1);
    }
    dphi_boson_met_isLoaded = true;
  }
  return dphi_boson_met_;
}
const float &SkimTree::dphi_lljets20_met() {
  if (not dphi_lljets20_met_isLoaded) {
    if (dphi_lljets20_met_branch != 0) {
      dphi_lljets20_met_branch->GetEntry(index);
    } else {
      printf("branch dphi_lljets20_met_branch does not exist!\n");
      exit(1);
    }
    dphi_lljets20_met_isLoaded = true;
  }
  return dphi_lljets20_met_;
}
const unsigned long long &SkimTree::event_event() {
  if (not event_event_isLoaded) {
    if (event_event_branch != 0) {
      event_event_branch->GetEntry(index);
    } else {
      printf("branch event_event_branch does not exist!\n");
      exit(1);
    }
    event_event_isLoaded = true;
  }
  return event_event_;
}
const float &SkimTree::eta_jet2() {
  if (not eta_jet2_isLoaded) {
    if (eta_jet2_branch != 0) {
      eta_jet2_branch->GetEntry(index);
    } else {
      printf("branch eta_jet2_branch does not exist!\n");
      exit(1);
    }
    eta_jet2_isLoaded = true;
  }
  return eta_jet2_;
}
const unsigned int &SkimTree::event_run() {
  if (not event_run_isLoaded) {
    if (event_run_branch != 0) {
      event_run_branch->GetEntry(index);
    } else {
      printf("branch event_run_branch does not exist!\n");
      exit(1);
    }
    event_run_isLoaded = true;
  }
  return event_run_;
}
const unsigned int &SkimTree::event_Njets20_btagged_medium() {
  if (not event_Njets20_btagged_medium_isLoaded) {
    if (event_Njets20_btagged_medium_branch != 0) {
      event_Njets20_btagged_medium_branch->GetEntry(index);
    } else {
      printf("branch event_Njets20_btagged_medium_branch does not exist!\n");
      exit(1);
    }
    event_Njets20_btagged_medium_isLoaded = true;
  }
  return event_Njets20_btagged_medium_;
}
const float &SkimTree::pt_boson() {
  if (not pt_boson_isLoaded) {
    if (pt_boson_branch != 0) {
      pt_boson_branch->GetEntry(index);
    } else {
      printf("branch pt_boson_branch does not exist!\n");
      exit(1);
    }
    pt_boson_isLoaded = true;
  }
  return pt_boson_;
}
const int &SkimTree::id_track() {
  if (not id_track_isLoaded) {
    if (id_track_branch != 0) {
      id_track_branch->GetEntry(index);
    } else {
      printf("branch id_track_branch does not exist!\n");
      exit(1);
    }
    id_track_isLoaded = true;
  }
  return id_track_;
}
const bool &SkimTree::is_mumu() {
  if (not is_mumu_isLoaded) {
    if (is_mumu_branch != 0) {
      is_mumu_branch->GetEntry(index);
    } else {
      printf("branch is_mumu_branch does not exist!\n");
      exit(1);
    }
    is_mumu_isLoaded = true;
  }
  return is_mumu_;
}
const float &SkimTree::event_wgt_trig_muon() {
  if (not event_wgt_trig_muon_isLoaded) {
    if (event_wgt_trig_muon_branch != 0) {
      event_wgt_trig_muon_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_trig_muon_branch does not exist!\n");
      exit(1);
    }
    event_wgt_trig_muon_isLoaded = true;
  }
  return event_wgt_trig_muon_;
}
const float &SkimTree::eta_photon() {
  if (not eta_photon_isLoaded) {
    if (eta_photon_branch != 0) {
      eta_photon_branch->GetEntry(index);
    } else {
      printf("branch eta_photon_branch does not exist!\n");
      exit(1);
    }
    eta_photon_isLoaded = true;
  }
  return eta_photon_;
}
const bool &SkimTree::is_VBFcat() {
  if (not is_VBFcat_isLoaded) {
    if (is_VBFcat_branch != 0) {
      is_VBFcat_branch->GetEntry(index);
    } else {
      printf("branch is_VBFcat_branch does not exist!\n");
      exit(1);
    }
    is_VBFcat_isLoaded = true;
  }
  return is_VBFcat_;
}
const float &SkimTree::event_wgt_pileup() {
  if (not event_wgt_pileup_isLoaded) {
    if (event_wgt_pileup_branch != 0) {
      event_wgt_pileup_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_pileup_branch does not exist!\n");
      exit(1);
    }
    event_wgt_pileup_isLoaded = true;
  }
  return event_wgt_pileup_;
}
const unsigned int &SkimTree::event_Njets_btagged_medium() {
  if (not event_Njets_btagged_medium_isLoaded) {
    if (event_Njets_btagged_medium_branch != 0) {
      event_Njets_btagged_medium_branch->GetEntry(index);
    } else {
      printf("branch event_Njets_btagged_medium_branch does not exist!\n");
      exit(1);
    }
    event_Njets_btagged_medium_isLoaded = true;
  }
  return event_Njets_btagged_medium_;
}
const unsigned int &SkimTree::event_Njets20_btagged_loose() {
  if (not event_Njets20_btagged_loose_isLoaded) {
    if (event_Njets20_btagged_loose_branch != 0) {
      event_Njets20_btagged_loose_branch->GetEntry(index);
    } else {
      printf("branch event_Njets20_btagged_loose_branch does not exist!\n");
      exit(1);
    }
    event_Njets20_btagged_loose_isLoaded = true;
  }
  return event_Njets20_btagged_loose_;
}
const float &SkimTree::eta_track() {
  if (not eta_track_isLoaded) {
    if (eta_track_branch != 0) {
      eta_track_branch->GetEntry(index);
    } else {
      printf("branch eta_track_branch does not exist!\n");
      exit(1);
    }
    eta_track_isLoaded = true;
  }
  return eta_track_;
}
const float &SkimTree::pt_track() {
  if (not pt_track_isLoaded) {
    if (pt_track_branch != 0) {
      pt_track_branch->GetEntry(index);
    } else {
      printf("branch pt_track_branch does not exist!\n");
      exit(1);
    }
    pt_track_isLoaded = true;
  }
  return pt_track_;
}
const float &SkimTree::eta_jet1() {
  if (not eta_jet1_isLoaded) {
    if (eta_jet1_branch != 0) {
      eta_jet1_branch->GetEntry(index);
    } else {
      printf("branch eta_jet1_branch does not exist!\n");
      exit(1);
    }
    eta_jet1_isLoaded = true;
  }
  return eta_jet1_;
}
const float &SkimTree::event_wgt() {
  if (not event_wgt_isLoaded) {
    if (event_wgt_branch != 0) {
      event_wgt_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_branch does not exist!\n");
      exit(1);
    }
    event_wgt_isLoaded = true;
  }
  return event_wgt_;
}
const bool &SkimTree::pass_trackveto() {
  if (not pass_trackveto_isLoaded) {
    if (pass_trackveto_branch != 0) {
      pass_trackveto_branch->GetEntry(index);
    } else {
      printf("branch pass_trackveto_branch does not exist!\n");
      exit(1);
    }
    pass_trackveto_isLoaded = true;
  }
  return pass_trackveto_;
}
const bool &SkimTree::is_ee() {
  if (not is_ee_isLoaded) {
    if (is_ee_branch != 0) {
      is_ee_branch->GetEntry(index);
    } else {
      printf("branch is_ee_branch does not exist!\n");
      exit(1);
    }
    is_ee_isLoaded = true;
  }
  return is_ee_;
}
const float &SkimTree::mass_boson() {
  if (not mass_boson_isLoaded) {
    if (mass_boson_branch != 0) {
      mass_boson_branch->GetEntry(index);
    } else {
      printf("branch mass_boson_branch does not exist!\n");
      exit(1);
    }
    mass_boson_isLoaded = true;
  }
  return mass_boson_;
}
const float &SkimTree::phi_boson() {
  if (not phi_boson_isLoaded) {
    if (phi_boson_branch != 0) {
      phi_boson_branch->GetEntry(index);
    } else {
      printf("branch phi_boson_branch does not exist!\n");
      exit(1);
    }
    phi_boson_isLoaded = true;
  }
  return phi_boson_;
}
const bool &SkimTree::has_photons_inHEM1516() {
  if (not has_photons_inHEM1516_isLoaded) {
    if (has_photons_inHEM1516_branch != 0) {
      has_photons_inHEM1516_branch->GetEntry(index);
    } else {
      printf("branch has_photons_inHEM1516_branch does not exist!\n");
      exit(1);
    }
    has_photons_inHEM1516_isLoaded = true;
  }
  return has_photons_inHEM1516_;
}
const float &SkimTree::mTZZ() {
  if (not mTZZ_isLoaded) {
    if (mTZZ_branch != 0) {
      mTZZ_branch->GetEntry(index);
    } else {
      printf("branch mTZZ_branch does not exist!\n");
      exit(1);
    }
    mTZZ_isLoaded = true;
  }
  return mTZZ_;
}
const float &SkimTree::eta_l1() {
  if (not eta_l1_isLoaded) {
    if (eta_l1_branch != 0) {
      eta_l1_branch->GetEntry(index);
    } else {
      printf("branch eta_l1_branch does not exist!\n");
      exit(1);
    }
    eta_l1_isLoaded = true;
  }
  return eta_l1_;
}
const bool &SkimTree::has_ak4jets_inHEM1516() {
  if (not has_ak4jets_inHEM1516_isLoaded) {
    if (has_ak4jets_inHEM1516_branch != 0) {
      has_ak4jets_inHEM1516_branch->GetEntry(index);
    } else {
      printf("branch has_ak4jets_inHEM1516_branch does not exist!\n");
      exit(1);
    }
    has_ak4jets_inHEM1516_isLoaded = true;
  }
  return has_ak4jets_inHEM1516_;
}
const int &SkimTree::id_l1() {
  if (not id_l1_isLoaded) {
    if (id_l1_branch != 0) {
      id_l1_branch->GetEntry(index);
    } else {
      printf("branch id_l1_branch does not exist!\n");
      exit(1);
    }
    id_l1_isLoaded = true;
  }
  return id_l1_;
}
const float &SkimTree::seedtime_photon() {
  if (not seedtime_photon_isLoaded) {
    if (seedtime_photon_branch != 0) {
      seedtime_photon_branch->GetEntry(index);
    } else {
      printf("branch seedtime_photon_branch does not exist!\n");
      exit(1);
    }
    seedtime_photon_isLoaded = true;
  }
  return seedtime_photon_;
}
const float &SkimTree::pt_jet2() {
  if (not pt_jet2_isLoaded) {
    if (pt_jet2_branch != 0) {
      pt_jet2_branch->GetEntry(index);
    } else {
      printf("branch pt_jet2_branch does not exist!\n");
      exit(1);
    }
    pt_jet2_isLoaded = true;
  }
  return pt_jet2_;
}
const float &SkimTree::sipip_photon() {
  if (not sipip_photon_isLoaded) {
    if (sipip_photon_branch != 0) {
      sipip_photon_branch->GetEntry(index);
    } else {
      printf("branch sipip_photon_branch does not exist!\n");
      exit(1);
    }
    sipip_photon_isLoaded = true;
  }
  return sipip_photon_;
}
const bool &SkimTree::pass_lepveto() {
  if (not pass_lepveto_isLoaded) {
    if (pass_lepveto_branch != 0) {
      pass_lepveto_branch->GetEntry(index);
    } else {
      printf("branch pass_lepveto_branch does not exist!\n");
      exit(1);
    }
    pass_lepveto_isLoaded = true;
  }
  return pass_lepveto_;
}
const float &SkimTree::phi_jet1() {
  if (not phi_jet1_isLoaded) {
    if (phi_jet1_branch != 0) {
      phi_jet1_branch->GetEntry(index);
    } else {
      printf("branch phi_jet1_branch does not exist!\n");
      exit(1);
    }
    phi_jet1_isLoaded = true;
  }
  return phi_jet1_;
}
const float &SkimTree::sieie_photon() {
  if (not sieie_photon_isLoaded) {
    if (sieie_photon_branch != 0) {
      sieie_photon_branch->GetEntry(index);
    } else {
      printf("branch sieie_photon_branch does not exist!\n");
      exit(1);
    }
    sieie_photon_isLoaded = true;
  }
  return sieie_photon_;
}
const bool &SkimTree::is_gamma() {
  if (not is_gamma_isLoaded) {
    if (is_gamma_branch != 0) {
      is_gamma_branch->GetEntry(index);
    } else {
      printf("branch is_gamma_branch does not exist!\n");
      exit(1);
    }
    is_gamma_isLoaded = true;
  }
  return is_gamma_;
}
const float &SkimTree::nhiso_photon() {
  if (not nhiso_photon_isLoaded) {
    if (nhiso_photon_branch != 0) {
      nhiso_photon_branch->GetEntry(index);
    } else {
      printf("branch nhiso_photon_branch does not exist!\n");
      exit(1);
    }
    nhiso_photon_isLoaded = true;
  }
  return nhiso_photon_;
}
const bool &SkimTree::isEB_photon() {
  if (not isEB_photon_isLoaded) {
    if (isEB_photon_branch != 0) {
      isEB_photon_branch->GetEntry(index);
    } else {
      printf("branch isEB_photon_branch does not exist!\n");
      exit(1);
    }
    isEB_photon_isLoaded = true;
  }
  return isEB_photon_;
}
const bool &SkimTree::has_electrons_inHEM1516() {
  if (not has_electrons_inHEM1516_isLoaded) {
    if (has_electrons_inHEM1516_branch != 0) {
      has_electrons_inHEM1516_branch->GetEntry(index);
    } else {
      printf("branch has_electrons_inHEM1516_branch does not exist!\n");
      exit(1);
    }
    has_electrons_inHEM1516_isLoaded = true;
  }
  return has_electrons_inHEM1516_;
}
const int &SkimTree::id_l2() {
  if (not id_l2_isLoaded) {
    if (id_l2_branch != 0) {
      id_l2_branch->GetEntry(index);
    } else {
      printf("branch id_l2_branch does not exist!\n");
      exit(1);
    }
    id_l2_isLoaded = true;
  }
  return id_l2_;
}
const float &SkimTree::event_wgt_OSSF_triggers() {
  if (not event_wgt_OSSF_triggers_isLoaded) {
    if (event_wgt_OSSF_triggers_branch != 0) {
      event_wgt_OSSF_triggers_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_OSSF_triggers_branch does not exist!\n");
      exit(1);
    }
    event_wgt_OSSF_triggers_isLoaded = true;
  }
  return event_wgt_OSSF_triggers_;
}
const float &SkimTree::eta_boson() {
  if (not eta_boson_isLoaded) {
    if (eta_boson_branch != 0) {
      eta_boson_branch->GetEntry(index);
    } else {
      printf("branch eta_boson_branch does not exist!\n");
      exit(1);
    }
    eta_boson_isLoaded = true;
  }
  return eta_boson_;
}
const float &SkimTree::phi_photon() {
  if (not phi_photon_isLoaded) {
    if (phi_photon_branch != 0) {
      phi_photon_branch->GetEntry(index);
    } else {
      printf("branch phi_photon_branch does not exist!\n");
      exit(1);
    }
    phi_photon_isLoaded = true;
  }
  return phi_photon_;
}
const float &SkimTree::event_DjjVBF() {
  if (not event_DjjVBF_isLoaded) {
    if (event_DjjVBF_branch != 0) {
      event_DjjVBF_branch->GetEntry(index);
    } else {
      printf("branch event_DjjVBF_branch does not exist!\n");
      exit(1);
    }
    event_DjjVBF_isLoaded = true;
  }
  return event_DjjVBF_;
}
const unsigned int &SkimTree::event_Njets() {
  if (not event_Njets_isLoaded) {
    if (event_Njets_branch != 0) {
      event_Njets_branch->GetEntry(index);
    } else {
      printf("branch event_Njets_branch does not exist!\n");
      exit(1);
    }
    event_Njets_isLoaded = true;
  }
  return event_Njets_;
}
const float &SkimTree::event_DjjVBF_rl() {
  if (not event_DjjVBF_rl_isLoaded) {
    if (event_DjjVBF_rl_branch != 0) {
      event_DjjVBF_rl_branch->GetEntry(index);
    } else {
      printf("branch event_DjjVBF_rl_branch does not exist!\n");
      exit(1);
    }
    event_DjjVBF_rl_isLoaded = true;
  }
  return event_DjjVBF_rl_;
}
const float &SkimTree::pfiso_photon() {
  if (not pfiso_photon_isLoaded) {
    if (pfiso_photon_branch != 0) {
      pfiso_photon_branch->GetEntry(index);
    } else {
      printf("branch pfiso_photon_branch does not exist!\n");
      exit(1);
    }
    pfiso_photon_isLoaded = true;
  }
  return pfiso_photon_;
}
const unsigned int &SkimTree::event_Njets20() {
  if (not event_Njets20_isLoaded) {
    if (event_Njets20_branch != 0) {
      event_Njets20_branch->GetEntry(index);
    } else {
      printf("branch event_Njets20_branch does not exist!\n");
      exit(1);
    }
    event_Njets20_isLoaded = true;
  }
  return event_Njets20_;
}
const float &SkimTree::event_wgt_trig_electron() {
  if (not event_wgt_trig_electron_isLoaded) {
    if (event_wgt_trig_electron_branch != 0) {
      event_wgt_trig_electron_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_trig_electron_branch does not exist!\n");
      exit(1);
    }
    event_wgt_trig_electron_isLoaded = true;
  }
  return event_wgt_trig_electron_;
}
const float &SkimTree::genmet_pt() {
  if (not genmet_pt_isLoaded) {
    if (genmet_pt_branch != 0) {
      genmet_pt_branch->GetEntry(index);
    } else {
      printf("branch genmet_pt_branch does not exist!\n");
      exit(1);
    }
    genmet_pt_isLoaded = true;
  }
  return genmet_pt_;
}
const float &SkimTree::r9_photon() {
  if (not r9_photon_isLoaded) {
    if (r9_photon_branch != 0) {
      r9_photon_branch->GetEntry(index);
    } else {
      printf("branch r9_photon_branch does not exist!\n");
      exit(1);
    }
    r9_photon_isLoaded = true;
  }
  return r9_photon_;
}
const float &SkimTree::phi_l2() {
  if (not phi_l2_isLoaded) {
    if (phi_l2_branch != 0) {
      phi_l2_branch->GetEntry(index);
    } else {
      printf("branch phi_l2_branch does not exist!\n");
      exit(1);
    }
    phi_l2_isLoaded = true;
  }
  return phi_l2_;
}
const float &SkimTree::event_wgt_L1prefire() {
  if (not event_wgt_L1prefire_isLoaded) {
    if (event_wgt_L1prefire_branch != 0) {
      event_wgt_L1prefire_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_L1prefire_branch does not exist!\n");
      exit(1);
    }
    event_wgt_L1prefire_isLoaded = true;
  }
  return event_wgt_L1prefire_;
}
const unsigned int &SkimTree::event_Nphotons() {
  if (not event_Nphotons_isLoaded) {
    if (event_Nphotons_branch != 0) {
      event_Nphotons_branch->GetEntry(index);
    } else {
      printf("branch event_Nphotons_branch does not exist!\n");
      exit(1);
    }
    event_Nphotons_isLoaded = true;
  }
  return event_Nphotons_;
}
const float &SkimTree::phi_l1() {
  if (not phi_l1_isLoaded) {
    if (phi_l1_branch != 0) {
      phi_l1_branch->GetEntry(index);
    } else {
      printf("branch phi_l1_branch does not exist!\n");
      exit(1);
    }
    phi_l1_isLoaded = true;
  }
  return phi_l1_;
}
const float &SkimTree::event_wgt_trig_photon() {
  if (not event_wgt_trig_photon_isLoaded) {
    if (event_wgt_trig_photon_branch != 0) {
      event_wgt_trig_photon_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_trig_photon_branch does not exist!\n");
      exit(1);
    }
    event_wgt_trig_photon_isLoaded = true;
  }
  return event_wgt_trig_photon_;
}
const float &SkimTree::dphi_jet20_met() {
  if (not dphi_jet20_met_isLoaded) {
    if (dphi_jet20_met_branch != 0) {
      dphi_jet20_met_branch->GetEntry(index);
    } else {
      printf("branch dphi_jet20_met_branch does not exist!\n");
      exit(1);
    }
    dphi_jet20_met_isLoaded = true;
  }
  return dphi_jet20_met_;
}
const float &SkimTree::mindphi_jet_met() {
  if (not mindphi_jet_met_isLoaded) {
    if (mindphi_jet_met_branch != 0) {
      mindphi_jet_met_branch->GetEntry(index);
    } else {
      printf("branch mindphi_jet_met_branch does not exist!\n");
      exit(1);
    }
    mindphi_jet_met_isLoaded = true;
  }
  return mindphi_jet_met_;
}
const float &SkimTree::pt_l2() {
  if (not pt_l2_isLoaded) {
    if (pt_l2_branch != 0) {
      pt_l2_branch->GetEntry(index);
    } else {
      printf("branch pt_l2_branch does not exist!\n");
      exit(1);
    }
    pt_l2_isLoaded = true;
  }
  return pt_l2_;
}
const float &SkimTree::genmet_phi() {
  if (not genmet_phi_isLoaded) {
    if (genmet_phi_branch != 0) {
      genmet_phi_branch->GetEntry(index);
    } else {
      printf("branch genmet_phi_branch does not exist!\n");
      exit(1);
    }
    genmet_phi_isLoaded = true;
  }
  return genmet_phi_;
}
const float &SkimTree::pt_l1() {
  if (not pt_l1_isLoaded) {
    if (pt_l1_branch != 0) {
      pt_l1_branch->GetEntry(index);
    } else {
      printf("branch pt_l1_branch does not exist!\n");
      exit(1);
    }
    pt_l1_isLoaded = true;
  }
  return pt_l1_;
}
const unsigned int &SkimTree::event_nvtxs_good() {
  if (not event_nvtxs_good_isLoaded) {
    if (event_nvtxs_good_branch != 0) {
      event_nvtxs_good_branch->GetEntry(index);
    } else {
      printf("branch event_nvtxs_good_branch does not exist!\n");
      exit(1);
    }
    event_nvtxs_good_isLoaded = true;
  }
  return event_nvtxs_good_;
}
const float &SkimTree::e4oe1_photon() {
  if (not e4oe1_photon_isLoaded) {
    if (e4oe1_photon_branch != 0) {
      e4oe1_photon_branch->GetEntry(index);
    } else {
      printf("branch e4oe1_photon_branch does not exist!\n");
      exit(1);
    }
    e4oe1_photon_isLoaded = true;
  }
  return e4oe1_photon_;
}
const float &SkimTree::phi_track() {
  if (not phi_track_isLoaded) {
    if (phi_track_branch != 0) {
      phi_track_branch->GetEntry(index);
    } else {
      printf("branch phi_track_branch does not exist!\n");
      exit(1);
    }
    phi_track_isLoaded = true;
  }
  return phi_track_;
}
const float &SkimTree::dphi_lljets_met() {
  if (not dphi_lljets_met_isLoaded) {
    if (dphi_lljets_met_branch != 0) {
      dphi_lljets_met_branch->GetEntry(index);
    } else {
      printf("branch dphi_lljets_met_branch does not exist!\n");
      exit(1);
    }
    dphi_lljets_met_isLoaded = true;
  }
  return dphi_lljets_met_;
}
const float &SkimTree::pTmiss() {
  if (not pTmiss_isLoaded) {
    if (pTmiss_branch != 0) {
      pTmiss_branch->GetEntry(index);
    } else {
      printf("branch pTmiss_branch does not exist!\n");
      exit(1);
    }
    pTmiss_isLoaded = true;
  }
  return pTmiss_;
}
const float &SkimTree::met_uncorr_pt() {
  if (not met_uncorr_pt_isLoaded) {
    if (met_uncorr_pt_branch != 0) {
      met_uncorr_pt_branch->GetEntry(index);
    } else {
      printf("branch met_uncorr_pt_branch does not exist!\n");
      exit(1);
    }
    met_uncorr_pt_isLoaded = true;
  }
  return met_uncorr_pt_;
}
const float &SkimTree::pt_jet1() {
  if (not pt_jet1_isLoaded) {
    if (pt_jet1_branch != 0) {
      pt_jet1_branch->GetEntry(index);
    } else {
      printf("branch pt_jet1_branch does not exist!\n");
      exit(1);
    }
    pt_jet1_isLoaded = true;
  }
  return pt_jet1_;
}
const unsigned int &SkimTree::event_lumi() {
  if (not event_lumi_isLoaded) {
    if (event_lumi_branch != 0) {
      event_lumi_branch->GetEntry(index);
    } else {
      printf("branch event_lumi_branch does not exist!\n");
      exit(1);
    }
    event_lumi_isLoaded = true;
  }
  return event_lumi_;
}
const float &SkimTree::eta_l2() {
  if (not eta_l2_isLoaded) {
    if (eta_l2_branch != 0) {
      eta_l2_branch->GetEntry(index);
    } else {
      printf("branch eta_l2_branch does not exist!\n");
      exit(1);
    }
    eta_l2_isLoaded = true;
  }
  return eta_l2_;
}
const float &SkimTree::phi_jet2() {
  if (not phi_jet2_isLoaded) {
    if (phi_jet2_branch != 0) {
      phi_jet2_branch->GetEntry(index);
    } else {
      printf("branch phi_jet2_branch does not exist!\n");
      exit(1);
    }
    phi_jet2_isLoaded = true;
  }
  return phi_jet2_;
}
const float &SkimTree::event_wgt_OSDF_triggers() {
  if (not event_wgt_OSDF_triggers_isLoaded) {
    if (event_wgt_OSDF_triggers_branch != 0) {
      event_wgt_OSDF_triggers_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_OSDF_triggers_branch does not exist!\n");
      exit(1);
    }
    event_wgt_OSDF_triggers_isLoaded = true;
  }
  return event_wgt_OSDF_triggers_;
}
const float &SkimTree::phimiss() {
  if (not phimiss_isLoaded) {
    if (phimiss_branch != 0) {
      phimiss_branch->GetEntry(index);
    } else {
      printf("branch phimiss_branch does not exist!\n");
      exit(1);
    }
    phimiss_isLoaded = true;
  }
  return phimiss_;
}
const unsigned int &SkimTree::event_HT() {
  if (not event_HT_isLoaded) {
    if (event_HT_branch != 0) {
      event_HT_branch->GetEntry(index);
    } else {
      printf("branch event_HT_branch does not exist!\n");
      exit(1);
    }
    event_HT_isLoaded = true;
  }
  return event_HT_;
}
const float &SkimTree::event_wgt_SFs() {
  if (not event_wgt_SFs_isLoaded) {
    if (event_wgt_SFs_branch != 0) {
      event_wgt_SFs_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_isLoaded = true;
  }
  return event_wgt_SFs_;
}
const unsigned int &SkimTree::event_Njets_btagged() {
  if (not event_Njets_btagged_isLoaded) {
    if (event_Njets_btagged_branch != 0) {
      event_Njets_btagged_branch->GetEntry(index);
    } else {
      printf("branch event_Njets_btagged_branch does not exist!\n");
      exit(1);
    }
    event_Njets_btagged_isLoaded = true;
  }
  return event_Njets_btagged_;
}
void SkimTree::progress( int nEventsTotal, int nEventsChain ){
  int period = 1000;
  if(nEventsTotal%1000 == 0) {
    if (isatty(1)) {
      if( ( nEventsChain - nEventsTotal ) > period ){
        float frac = (float)nEventsTotal/(nEventsChain*0.01);
        printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
               "\033[0m\033[32m <---\033[0m\015", frac);
        fflush(stdout);
      }
      else {
        printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
               "\033[0m\033[32m <---\033[0m\015", 100.);
        cout << endl;
      }
    }
  }
}
namespace tas {
  const float &pt_photon() { return st.pt_photon(); }
  const float &chiso_photon() { return st.chiso_photon(); }
  const bool &id_photon_Hgg() { return st.id_photon_Hgg(); }
  const bool &is_emu() { return st.is_emu(); }
  const float &mipE_photon() { return st.mipE_photon(); }
  const float &met_uncorr_phi() { return st.met_uncorr_phi(); }
  const float &mZZ() { return st.mZZ(); }
  const bool &convveto_photon() { return st.convveto_photon(); }
  const float &dphi_boson_met() { return st.dphi_boson_met(); }
  const float &dphi_lljets20_met() { return st.dphi_lljets20_met(); }
  const unsigned long long &event_event() { return st.event_event(); }
  const float &eta_jet2() { return st.eta_jet2(); }
  const unsigned int &event_run() { return st.event_run(); }
  const unsigned int &event_Njets20_btagged_medium() { return st.event_Njets20_btagged_medium(); }
  const float &pt_boson() { return st.pt_boson(); }
  const int &id_track() { return st.id_track(); }
  const bool &is_mumu() { return st.is_mumu(); }
  const float &event_wgt_trig_muon() { return st.event_wgt_trig_muon(); }
  const float &eta_photon() { return st.eta_photon(); }
  const bool &is_VBFcat() { return st.is_VBFcat(); }
  const float &event_wgt_pileup() { return st.event_wgt_pileup(); }
  const unsigned int &event_Njets_btagged_medium() { return st.event_Njets_btagged_medium(); }
  const unsigned int &event_Njets20_btagged_loose() { return st.event_Njets20_btagged_loose(); }
  const float &eta_track() { return st.eta_track(); }
  const float &pt_track() { return st.pt_track(); }
  const float &eta_jet1() { return st.eta_jet1(); }
  const float &event_wgt() { return st.event_wgt(); }
  const bool &pass_trackveto() { return st.pass_trackveto(); }
  const bool &is_ee() { return st.is_ee(); }
  const float &mass_boson() { return st.mass_boson(); }
  const float &phi_boson() { return st.phi_boson(); }
  const bool &has_photons_inHEM1516() { return st.has_photons_inHEM1516(); }
  const float &mTZZ() { return st.mTZZ(); }
  const float &eta_l1() { return st.eta_l1(); }
  const bool &has_ak4jets_inHEM1516() { return st.has_ak4jets_inHEM1516(); }
  const int &id_l1() { return st.id_l1(); }
  const float &seedtime_photon() { return st.seedtime_photon(); }
  const float &pt_jet2() { return st.pt_jet2(); }
  const float &sipip_photon() { return st.sipip_photon(); }
  const bool &pass_lepveto() { return st.pass_lepveto(); }
  const float &phi_jet1() { return st.phi_jet1(); }
  const float &sieie_photon() { return st.sieie_photon(); }
  const bool &is_gamma() { return st.is_gamma(); }
  const float &nhiso_photon() { return st.nhiso_photon(); }
  const bool &isEB_photon() { return st.isEB_photon(); }
  const bool &has_electrons_inHEM1516() { return st.has_electrons_inHEM1516(); }
  const int &id_l2() { return st.id_l2(); }
  const float &event_wgt_OSSF_triggers() { return st.event_wgt_OSSF_triggers(); }
  const float &eta_boson() { return st.eta_boson(); }
  const float &phi_photon() { return st.phi_photon(); }
  const float &event_DjjVBF() { return st.event_DjjVBF(); }
  const unsigned int &event_Njets() { return st.event_Njets(); }
  const float &event_DjjVBF_rl() { return st.event_DjjVBF_rl(); }
  const float &pfiso_photon() { return st.pfiso_photon(); }
  const unsigned int &event_Njets20() { return st.event_Njets20(); }
  const float &event_wgt_trig_electron() { return st.event_wgt_trig_electron(); }
  const float &genmet_pt() { return st.genmet_pt(); }
  const float &r9_photon() { return st.r9_photon(); }
  const float &phi_l2() { return st.phi_l2(); }
  const float &event_wgt_L1prefire() { return st.event_wgt_L1prefire(); }
  const unsigned int &event_Nphotons() { return st.event_Nphotons(); }
  const float &phi_l1() { return st.phi_l1(); }
  const float &event_wgt_trig_photon() { return st.event_wgt_trig_photon(); }
  const float &dphi_jet20_met() { return st.dphi_jet20_met(); }
  const float &mindphi_jet_met() { return st.mindphi_jet_met(); }
  const float &pt_l2() { return st.pt_l2(); }
  const float &genmet_phi() { return st.genmet_phi(); }
  const float &pt_l1() { return st.pt_l1(); }
  const unsigned int &event_nvtxs_good() { return st.event_nvtxs_good(); }
  const float &e4oe1_photon() { return st.e4oe1_photon(); }
  const float &phi_track() { return st.phi_track(); }
  const float &dphi_lljets_met() { return st.dphi_lljets_met(); }
  const float &pTmiss() { return st.pTmiss(); }
  const float &met_uncorr_pt() { return st.met_uncorr_pt(); }
  const float &pt_jet1() { return st.pt_jet1(); }
  const unsigned int &event_lumi() { return st.event_lumi(); }
  const float &eta_l2() { return st.eta_l2(); }
  const float &phi_jet2() { return st.phi_jet2(); }
  const float &event_wgt_OSDF_triggers() { return st.event_wgt_OSDF_triggers(); }
  const float &phimiss() { return st.phimiss(); }
  const unsigned int &event_HT() { return st.event_HT(); }
  const float &event_wgt_SFs() { return st.event_wgt_SFs(); }
  const unsigned int &event_Njets_btagged() { return st.event_Njets_btagged(); }
}
