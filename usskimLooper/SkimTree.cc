#include "SkimTree.h"
SkimTree st;

void SkimTree::Init(TTree *tree) {
  tree->SetMakeClass(1);
  EventNumber_branch = 0;
  if (tree->GetBranch("EventNumber") != 0) {
    EventNumber_branch = tree->GetBranch("EventNumber");
    if (EventNumber_branch) { EventNumber_branch->SetAddress(&EventNumber_); }
  }
  KFactor_EW_NLO_qqVV_Bkg_EWDn_branch = 0;
  if (tree->GetBranch("KFactor_EW_NLO_qqVV_Bkg_EWDn") != 0) {
    KFactor_EW_NLO_qqVV_Bkg_EWDn_branch = tree->GetBranch("KFactor_EW_NLO_qqVV_Bkg_EWDn");
    if (KFactor_EW_NLO_qqVV_Bkg_EWDn_branch) { KFactor_EW_NLO_qqVV_Bkg_EWDn_branch->SetAddress(&KFactor_EW_NLO_qqVV_Bkg_EWDn_); }
  }
  KFactor_EW_NLO_qqVV_Bkg_EWUp_branch = 0;
  if (tree->GetBranch("KFactor_EW_NLO_qqVV_Bkg_EWUp") != 0) {
    KFactor_EW_NLO_qqVV_Bkg_EWUp_branch = tree->GetBranch("KFactor_EW_NLO_qqVV_Bkg_EWUp");
    if (KFactor_EW_NLO_qqVV_Bkg_EWUp_branch) { KFactor_EW_NLO_qqVV_Bkg_EWUp_branch->SetAddress(&KFactor_EW_NLO_qqVV_Bkg_EWUp_); }
  }
  KFactor_EW_NLO_qqVV_Bkg_Nominal_branch = 0;
  if (tree->GetBranch("KFactor_EW_NLO_qqVV_Bkg_Nominal") != 0) {
    KFactor_EW_NLO_qqVV_Bkg_Nominal_branch = tree->GetBranch("KFactor_EW_NLO_qqVV_Bkg_Nominal");
    if (KFactor_EW_NLO_qqVV_Bkg_Nominal_branch) { KFactor_EW_NLO_qqVV_Bkg_Nominal_branch->SetAddress(&KFactor_EW_NLO_qqVV_Bkg_Nominal_); }
  }
  KFactor_EW_NLO_qqVV_Bkg_arg_mass_branch = 0;
  if (tree->GetBranch("KFactor_EW_NLO_qqVV_Bkg_arg_mass") != 0) {
    KFactor_EW_NLO_qqVV_Bkg_arg_mass_branch = tree->GetBranch("KFactor_EW_NLO_qqVV_Bkg_arg_mass");
    if (KFactor_EW_NLO_qqVV_Bkg_arg_mass_branch) { KFactor_EW_NLO_qqVV_Bkg_arg_mass_branch->SetAddress(&KFactor_EW_NLO_qqVV_Bkg_arg_mass_); }
  }
  KFactor_EW_NLO_qqVV_Bkg_arg_pthat_branch = 0;
  if (tree->GetBranch("KFactor_EW_NLO_qqVV_Bkg_arg_pthat") != 0) {
    KFactor_EW_NLO_qqVV_Bkg_arg_pthat_branch = tree->GetBranch("KFactor_EW_NLO_qqVV_Bkg_arg_pthat");
    if (KFactor_EW_NLO_qqVV_Bkg_arg_pthat_branch) { KFactor_EW_NLO_qqVV_Bkg_arg_pthat_branch->SetAddress(&KFactor_EW_NLO_qqVV_Bkg_arg_pthat_); }
  }
  KFactor_EW_NLO_qqVV_Bkg_arg_rho_branch = 0;
  if (tree->GetBranch("KFactor_EW_NLO_qqVV_Bkg_arg_rho") != 0) {
    KFactor_EW_NLO_qqVV_Bkg_arg_rho_branch = tree->GetBranch("KFactor_EW_NLO_qqVV_Bkg_arg_rho");
    if (KFactor_EW_NLO_qqVV_Bkg_arg_rho_branch) { KFactor_EW_NLO_qqVV_Bkg_arg_rho_branch->SetAddress(&KFactor_EW_NLO_qqVV_Bkg_arg_rho_); }
  }
  KFactor_EW_NLO_qqVV_Bkg_arg_that_branch = 0;
  if (tree->GetBranch("KFactor_EW_NLO_qqVV_Bkg_arg_that") != 0) {
    KFactor_EW_NLO_qqVV_Bkg_arg_that_branch = tree->GetBranch("KFactor_EW_NLO_qqVV_Bkg_arg_that");
    if (KFactor_EW_NLO_qqVV_Bkg_arg_that_branch) { KFactor_EW_NLO_qqVV_Bkg_arg_that_branch->SetAddress(&KFactor_EW_NLO_qqVV_Bkg_arg_that_); }
  }
  KFactor_QCD_NNLO_qqVV_Bkg_Nominal_branch = 0;
  if (tree->GetBranch("KFactor_QCD_NNLO_qqVV_Bkg_Nominal") != 0) {
    KFactor_QCD_NNLO_qqVV_Bkg_Nominal_branch = tree->GetBranch("KFactor_QCD_NNLO_qqVV_Bkg_Nominal");
    if (KFactor_QCD_NNLO_qqVV_Bkg_Nominal_branch) { KFactor_QCD_NNLO_qqVV_Bkg_Nominal_branch->SetAddress(&KFactor_QCD_NNLO_qqVV_Bkg_Nominal_); }
  }
  KFactor_QCD_NNLO_qqVV_Bkg_arg_mass_branch = 0;
  if (tree->GetBranch("KFactor_QCD_NNLO_qqVV_Bkg_arg_mass") != 0) {
    KFactor_QCD_NNLO_qqVV_Bkg_arg_mass_branch = tree->GetBranch("KFactor_QCD_NNLO_qqVV_Bkg_arg_mass");
    if (KFactor_QCD_NNLO_qqVV_Bkg_arg_mass_branch) { KFactor_QCD_NNLO_qqVV_Bkg_arg_mass_branch->SetAddress(&KFactor_QCD_NNLO_qqVV_Bkg_arg_mass_); }
  }
  KFactor_QCD_NNLO_qqVV_Bkg_arg_pthat_branch = 0;
  if (tree->GetBranch("KFactor_QCD_NNLO_qqVV_Bkg_arg_pthat") != 0) {
    KFactor_QCD_NNLO_qqVV_Bkg_arg_pthat_branch = tree->GetBranch("KFactor_QCD_NNLO_qqVV_Bkg_arg_pthat");
    if (KFactor_QCD_NNLO_qqVV_Bkg_arg_pthat_branch) { KFactor_QCD_NNLO_qqVV_Bkg_arg_pthat_branch->SetAddress(&KFactor_QCD_NNLO_qqVV_Bkg_arg_pthat_); }
  }
  LuminosityBlock_branch = 0;
  if (tree->GetBranch("LuminosityBlock") != 0) {
    LuminosityBlock_branch = tree->GetBranch("LuminosityBlock");
    if (LuminosityBlock_branch) { LuminosityBlock_branch->SetAddress(&LuminosityBlock_); }
  }
  RunNumber_branch = 0;
  if (tree->GetBranch("RunNumber") != 0) {
    RunNumber_branch = tree->GetBranch("RunNumber");
    if (RunNumber_branch) { RunNumber_branch->SetAddress(&RunNumber_); }
  }
  ak4jets_CEMF_branch = 0;
  if (tree->GetBranch("ak4jets_CEMF") != 0) {
    ak4jets_CEMF_branch = tree->GetBranch("ak4jets_CEMF");
    if (ak4jets_CEMF_branch) { ak4jets_CEMF_branch->SetAddress(&ak4jets_CEMF_); }
  }
  ak4jets_HT_branch = 0;
  if (tree->GetBranch("ak4jets_HT") != 0) {
    ak4jets_HT_branch = tree->GetBranch("ak4jets_HT");
    if (ak4jets_HT_branch) { ak4jets_HT_branch->SetAddress(&ak4jets_HT_); }
  }
  ak4jets_MHT_branch = 0;
  if (tree->GetBranch("ak4jets_MHT") != 0) {
    ak4jets_MHT_branch = tree->GetBranch("ak4jets_MHT");
    if (ak4jets_MHT_branch) { ak4jets_MHT_branch->SetAddress(&ak4jets_MHT_); }
  }
  ak4jets_NEMF_branch = 0;
  if (tree->GetBranch("ak4jets_NEMF") != 0) {
    ak4jets_NEMF_branch = tree->GetBranch("ak4jets_NEMF");
    if (ak4jets_NEMF_branch) { ak4jets_NEMF_branch->SetAddress(&ak4jets_NEMF_); }
  }
  ak4jets_btagWP_Bits_branch = 0;
  if (tree->GetBranch("ak4jets_btagWP_Bits") != 0) {
    ak4jets_btagWP_Bits_branch = tree->GetBranch("ak4jets_btagWP_Bits");
    if (ak4jets_btagWP_Bits_branch) { ak4jets_btagWP_Bits_branch->SetAddress(&ak4jets_btagWP_Bits_); }
  }
  ak4jets_eta_branch = 0;
  if (tree->GetBranch("ak4jets_eta") != 0) {
    ak4jets_eta_branch = tree->GetBranch("ak4jets_eta");
    if (ak4jets_eta_branch) { ak4jets_eta_branch->SetAddress(&ak4jets_eta_); }
  }
  ak4jets_is_genMatched_branch = 0;
  if (tree->GetBranch("ak4jets_is_genMatched") != 0) {
    ak4jets_is_genMatched_branch = tree->GetBranch("ak4jets_is_genMatched");
    if (ak4jets_is_genMatched_branch) { ak4jets_is_genMatched_branch->SetAddress(&ak4jets_is_genMatched_); }
  }
  ak4jets_is_genMatched_fullCone_branch = 0;
  if (tree->GetBranch("ak4jets_is_genMatched_fullCone") != 0) {
    ak4jets_is_genMatched_fullCone_branch = tree->GetBranch("ak4jets_is_genMatched_fullCone");
    if (ak4jets_is_genMatched_fullCone_branch) { ak4jets_is_genMatched_fullCone_branch->SetAddress(&ak4jets_is_genMatched_fullCone_); }
  }
  ak4jets_mass_branch = 0;
  if (tree->GetBranch("ak4jets_mass") != 0) {
    ak4jets_mass_branch = tree->GetBranch("ak4jets_mass");
    if (ak4jets_mass_branch) { ak4jets_mass_branch->SetAddress(&ak4jets_mass_); }
  }
  ak4jets_phi_branch = 0;
  if (tree->GetBranch("ak4jets_phi") != 0) {
    ak4jets_phi_branch = tree->GetBranch("ak4jets_phi");
    if (ak4jets_phi_branch) { ak4jets_phi_branch->SetAddress(&ak4jets_phi_); }
  }
  ak4jets_pt_branch = 0;
  if (tree->GetBranch("ak4jets_pt") != 0) {
    ak4jets_pt_branch = tree->GetBranch("ak4jets_pt");
    if (ak4jets_pt_branch) { ak4jets_pt_branch->SetAddress(&ak4jets_pt_); }
  }
  ak8jets_eta_branch = 0;
  if (tree->GetBranch("ak8jets_eta") != 0) {
    ak8jets_eta_branch = tree->GetBranch("ak8jets_eta");
    if (ak8jets_eta_branch) { ak8jets_eta_branch->SetAddress(&ak8jets_eta_); }
  }
  ak8jets_mass_branch = 0;
  if (tree->GetBranch("ak8jets_mass") != 0) {
    ak8jets_mass_branch = tree->GetBranch("ak8jets_mass");
    if (ak8jets_mass_branch) { ak8jets_mass_branch->SetAddress(&ak8jets_mass_); }
  }
  ak8jets_phi_branch = 0;
  if (tree->GetBranch("ak8jets_phi") != 0) {
    ak8jets_phi_branch = tree->GetBranch("ak8jets_phi");
    if (ak8jets_phi_branch) { ak8jets_phi_branch->SetAddress(&ak8jets_phi_); }
  }
  ak8jets_pt_branch = 0;
  if (tree->GetBranch("ak8jets_pt") != 0) {
    ak8jets_pt_branch = tree->GetBranch("ak8jets_pt");
    if (ak8jets_pt_branch) { ak8jets_pt_branch->SetAddress(&ak8jets_pt_); }
  }
  dPhi_pTboson_pTmiss_branch = 0;
  if (tree->GetBranch("dPhi_pTboson_pTmiss") != 0) {
    dPhi_pTboson_pTmiss_branch = tree->GetBranch("dPhi_pTboson_pTmiss");
    if (dPhi_pTboson_pTmiss_branch) { dPhi_pTboson_pTmiss_branch->SetAddress(&dPhi_pTboson_pTmiss_); }
  }
  dPhi_pTbosonjets_pTmiss_branch = 0;
  if (tree->GetBranch("dPhi_pTbosonjets_pTmiss") != 0) {
    dPhi_pTbosonjets_pTmiss_branch = tree->GetBranch("dPhi_pTbosonjets_pTmiss");
    if (dPhi_pTbosonjets_pTmiss_branch) { dPhi_pTbosonjets_pTmiss_branch->SetAddress(&dPhi_pTbosonjets_pTmiss_); }
  }
  dilepton_eta_branch = 0;
  if (tree->GetBranch("dilepton_eta") != 0) {
    dilepton_eta_branch = tree->GetBranch("dilepton_eta");
    if (dilepton_eta_branch) { dilepton_eta_branch->SetAddress(&dilepton_eta_); }
  }
  dilepton_id_branch = 0;
  if (tree->GetBranch("dilepton_id") != 0) {
    dilepton_id_branch = tree->GetBranch("dilepton_id");
    if (dilepton_id_branch) { dilepton_id_branch->SetAddress(&dilepton_id_); }
  }
  dilepton_mass_branch = 0;
  if (tree->GetBranch("dilepton_mass") != 0) {
    dilepton_mass_branch = tree->GetBranch("dilepton_mass");
    if (dilepton_mass_branch) { dilepton_mass_branch->SetAddress(&dilepton_mass_); }
  }
  dilepton_phi_branch = 0;
  if (tree->GetBranch("dilepton_phi") != 0) {
    dilepton_phi_branch = tree->GetBranch("dilepton_phi");
    if (dilepton_phi_branch) { dilepton_phi_branch->SetAddress(&dilepton_phi_); }
  }
  dilepton_pt_branch = 0;
  if (tree->GetBranch("dilepton_pt") != 0) {
    dilepton_pt_branch = tree->GetBranch("dilepton_pt");
    if (dilepton_pt_branch) { dilepton_pt_branch->SetAddress(&dilepton_pt_); }
  }
  electron_full5x5_r9_branch = 0;
  if (tree->GetBranch("electron_full5x5_r9") != 0) {
    electron_full5x5_r9_branch = tree->GetBranch("electron_full5x5_r9");
    if (electron_full5x5_r9_branch) { electron_full5x5_r9_branch->SetAddress(&electron_full5x5_r9_); }
  }
  electron_full5x5_sigmaIEtaIEta_branch = 0;
  if (tree->GetBranch("electron_full5x5_sigmaIEtaIEta") != 0) {
    electron_full5x5_sigmaIEtaIEta_branch = tree->GetBranch("electron_full5x5_sigmaIEtaIEta");
    if (electron_full5x5_sigmaIEtaIEta_branch) { electron_full5x5_sigmaIEtaIEta_branch->SetAddress(&electron_full5x5_sigmaIEtaIEta_); }
  }
  electron_full5x5_sigmaIPhiIPhi_branch = 0;
  if (tree->GetBranch("electron_full5x5_sigmaIPhiIPhi") != 0) {
    electron_full5x5_sigmaIPhiIPhi_branch = tree->GetBranch("electron_full5x5_sigmaIPhiIPhi");
    if (electron_full5x5_sigmaIPhiIPhi_branch) { electron_full5x5_sigmaIPhiIPhi_branch->SetAddress(&electron_full5x5_sigmaIPhiIPhi_); }
  }
  electron_seedTime_branch = 0;
  if (tree->GetBranch("electron_seedTime") != 0) {
    electron_seedTime_branch = tree->GetBranch("electron_seedTime");
    if (electron_seedTime_branch) { electron_seedTime_branch->SetAddress(&electron_seedTime_); }
  }
  electrons_full5x5_r9_branch = 0;
  if (tree->GetBranch("electrons_full5x5_r9") != 0) {
    electrons_full5x5_r9_branch = tree->GetBranch("electrons_full5x5_r9");
    if (electrons_full5x5_r9_branch) { electrons_full5x5_r9_branch->SetAddress(&electrons_full5x5_r9_); }
  }
  electrons_full5x5_sigmaIEtaIEta_branch = 0;
  if (tree->GetBranch("electrons_full5x5_sigmaIEtaIEta") != 0) {
    electrons_full5x5_sigmaIEtaIEta_branch = tree->GetBranch("electrons_full5x5_sigmaIEtaIEta");
    if (electrons_full5x5_sigmaIEtaIEta_branch) { electrons_full5x5_sigmaIEtaIEta_branch->SetAddress(&electrons_full5x5_sigmaIEtaIEta_); }
  }
  electrons_full5x5_sigmaIPhiIPhi_branch = 0;
  if (tree->GetBranch("electrons_full5x5_sigmaIPhiIPhi") != 0) {
    electrons_full5x5_sigmaIPhiIPhi_branch = tree->GetBranch("electrons_full5x5_sigmaIPhiIPhi");
    if (electrons_full5x5_sigmaIPhiIPhi_branch) { electrons_full5x5_sigmaIPhiIPhi_branch->SetAddress(&electrons_full5x5_sigmaIPhiIPhi_); }
  }
  electrons_seedTime_branch = 0;
  if (tree->GetBranch("electrons_seedTime") != 0) {
    electrons_seedTime_branch = tree->GetBranch("electrons_seedTime");
    if (electrons_seedTime_branch) { electrons_seedTime_branch->SetAddress(&electrons_seedTime_); }
  }
  event_mTZZ_branch = 0;
  if (tree->GetBranch("event_mTZZ") != 0) {
    event_mTZZ_branch = tree->GetBranch("event_mTZZ");
    if (event_mTZZ_branch) { event_mTZZ_branch->SetAddress(&event_mTZZ_); }
  }
  event_mZZ_branch = 0;
  if (tree->GetBranch("event_mZZ") != 0) {
    event_mZZ_branch = tree->GetBranch("event_mZZ");
    if (event_mZZ_branch) { event_mZZ_branch->SetAddress(&event_mZZ_); }
  }
  event_mlg_branch = 0;
  if (tree->GetBranch("event_mlg") != 0) {
    event_mlg_branch = tree->GetBranch("event_mlg");
    if (event_mlg_branch) { event_mlg_branch->SetAddress(&event_mlg_); }
  }
  event_mllg_branch = 0;
  if (tree->GetBranch("event_mllg") != 0) {
    event_mllg_branch = tree->GetBranch("event_mllg");
    if (event_mllg_branch) { event_mllg_branch->SetAddress(&event_mllg_); }
  }
  event_n_ak4jets_pt20_branch = 0;
  if (tree->GetBranch("event_n_ak4jets_pt20") != 0) {
    event_n_ak4jets_pt20_branch = tree->GetBranch("event_n_ak4jets_pt20");
    if (event_n_ak4jets_pt20_branch) { event_n_ak4jets_pt20_branch->SetAddress(&event_n_ak4jets_pt20_); }
  }
  event_n_ak4jets_pt20_btagged_loose_branch = 0;
  if (tree->GetBranch("event_n_ak4jets_pt20_btagged_loose") != 0) {
    event_n_ak4jets_pt20_btagged_loose_branch = tree->GetBranch("event_n_ak4jets_pt20_btagged_loose");
    if (event_n_ak4jets_pt20_btagged_loose_branch) { event_n_ak4jets_pt20_btagged_loose_branch->SetAddress(&event_n_ak4jets_pt20_btagged_loose_); }
  }
  event_n_ak4jets_pt20_btagged_medium_branch = 0;
  if (tree->GetBranch("event_n_ak4jets_pt20_btagged_medium") != 0) {
    event_n_ak4jets_pt20_btagged_medium_branch = tree->GetBranch("event_n_ak4jets_pt20_btagged_medium");
    if (event_n_ak4jets_pt20_btagged_medium_branch) { event_n_ak4jets_pt20_btagged_medium_branch->SetAddress(&event_n_ak4jets_pt20_btagged_medium_); }
  }
  event_n_ak4jets_pt30_branch = 0;
  if (tree->GetBranch("event_n_ak4jets_pt30") != 0) {
    event_n_ak4jets_pt30_branch = tree->GetBranch("event_n_ak4jets_pt30");
    if (event_n_ak4jets_pt30_branch) { event_n_ak4jets_pt30_branch->SetAddress(&event_n_ak4jets_pt30_); }
  }
  event_n_ak4jets_pt30_btagged_loose_branch = 0;
  if (tree->GetBranch("event_n_ak4jets_pt30_btagged_loose") != 0) {
    event_n_ak4jets_pt30_btagged_loose_branch = tree->GetBranch("event_n_ak4jets_pt30_btagged_loose");
    if (event_n_ak4jets_pt30_btagged_loose_branch) { event_n_ak4jets_pt30_btagged_loose_branch->SetAddress(&event_n_ak4jets_pt30_btagged_loose_); }
  }
  event_n_ak4jets_pt30_btagged_medium_branch = 0;
  if (tree->GetBranch("event_n_ak4jets_pt30_btagged_medium") != 0) {
    event_n_ak4jets_pt30_btagged_medium_branch = tree->GetBranch("event_n_ak4jets_pt30_btagged_medium");
    if (event_n_ak4jets_pt30_btagged_medium_branch) { event_n_ak4jets_pt30_btagged_medium_branch->SetAddress(&event_n_ak4jets_pt30_btagged_medium_); }
  }
  event_n_leptons_fakeableBase_branch = 0;
  if (tree->GetBranch("event_n_leptons_fakeableBase") != 0) {
    event_n_leptons_fakeableBase_branch = tree->GetBranch("event_n_leptons_fakeableBase");
    if (event_n_leptons_fakeableBase_branch) { event_n_leptons_fakeableBase_branch->SetAddress(&event_n_leptons_fakeableBase_); }
  }
  event_n_vtxs_good_branch = 0;
  if (tree->GetBranch("event_n_vtxs_good") != 0) {
    event_n_vtxs_good_branch = tree->GetBranch("event_n_vtxs_good");
    if (event_n_vtxs_good_branch) { event_n_vtxs_good_branch->SetAddress(&event_n_vtxs_good_); }
  }
  event_pTmiss_branch = 0;
  if (tree->GetBranch("event_pTmiss") != 0) {
    event_pTmiss_branch = tree->GetBranch("event_pTmiss");
    if (event_pTmiss_branch) { event_pTmiss_branch->SetAddress(&event_pTmiss_); }
  }
  event_pass_tightMETFilters_branch = 0;
  if (tree->GetBranch("event_pass_tightMETFilters") != 0) {
    event_pass_tightMETFilters_branch = tree->GetBranch("event_pass_tightMETFilters");
    if (event_pass_tightMETFilters_branch) { event_pass_tightMETFilters_branch->SetAddress(&event_pass_tightMETFilters_); }
  }
  event_phimiss_branch = 0;
  if (tree->GetBranch("event_phimiss") != 0) {
    event_phimiss_branch = tree->GetBranch("event_phimiss");
    if (event_phimiss_branch) { event_phimiss_branch->SetAddress(&event_phimiss_); }
  }
  event_wgt_branch = 0;
  if (tree->GetBranch("event_wgt") != 0) {
    event_wgt_branch = tree->GetBranch("event_wgt");
    if (event_wgt_branch) { event_wgt_branch->SetAddress(&event_wgt_); }
  }
  event_wgt_L1PrefiringDn_branch = 0;
  if (tree->GetBranch("event_wgt_L1PrefiringDn") != 0) {
    event_wgt_L1PrefiringDn_branch = tree->GetBranch("event_wgt_L1PrefiringDn");
    if (event_wgt_L1PrefiringDn_branch) { event_wgt_L1PrefiringDn_branch->SetAddress(&event_wgt_L1PrefiringDn_); }
  }
  event_wgt_L1PrefiringUp_branch = 0;
  if (tree->GetBranch("event_wgt_L1PrefiringUp") != 0) {
    event_wgt_L1PrefiringUp_branch = tree->GetBranch("event_wgt_L1PrefiringUp");
    if (event_wgt_L1PrefiringUp_branch) { event_wgt_L1PrefiringUp_branch->SetAddress(&event_wgt_L1PrefiringUp_); }
  }
  event_wgt_PUDn_branch = 0;
  if (tree->GetBranch("event_wgt_PUDn") != 0) {
    event_wgt_PUDn_branch = tree->GetBranch("event_wgt_PUDn");
    if (event_wgt_PUDn_branch) { event_wgt_PUDn_branch->SetAddress(&event_wgt_PUDn_); }
  }
  event_wgt_PUUp_branch = 0;
  if (tree->GetBranch("event_wgt_PUUp") != 0) {
    event_wgt_PUUp_branch = tree->GetBranch("event_wgt_PUUp");
    if (event_wgt_PUUp_branch) { event_wgt_PUUp_branch->SetAddress(&event_wgt_PUUp_); }
  }
  event_wgt_SFs_PUJetId_branch = 0;
  if (tree->GetBranch("event_wgt_SFs_PUJetId") != 0) {
    event_wgt_SFs_PUJetId_branch = tree->GetBranch("event_wgt_SFs_PUJetId");
    if (event_wgt_SFs_PUJetId_branch) { event_wgt_SFs_PUJetId_branch->SetAddress(&event_wgt_SFs_PUJetId_); }
  }
  event_wgt_SFs_PUJetId_EffDn_branch = 0;
  if (tree->GetBranch("event_wgt_SFs_PUJetId_EffDn") != 0) {
    event_wgt_SFs_PUJetId_EffDn_branch = tree->GetBranch("event_wgt_SFs_PUJetId_EffDn");
    if (event_wgt_SFs_PUJetId_EffDn_branch) { event_wgt_SFs_PUJetId_EffDn_branch->SetAddress(&event_wgt_SFs_PUJetId_EffDn_); }
  }
  event_wgt_SFs_PUJetId_EffUp_branch = 0;
  if (tree->GetBranch("event_wgt_SFs_PUJetId_EffUp") != 0) {
    event_wgt_SFs_PUJetId_EffUp_branch = tree->GetBranch("event_wgt_SFs_PUJetId_EffUp");
    if (event_wgt_SFs_PUJetId_EffUp_branch) { event_wgt_SFs_PUJetId_EffUp_branch->SetAddress(&event_wgt_SFs_PUJetId_EffUp_); }
  }
  event_wgt_SFs_btagging_branch = 0;
  if (tree->GetBranch("event_wgt_SFs_btagging") != 0) {
    event_wgt_SFs_btagging_branch = tree->GetBranch("event_wgt_SFs_btagging");
    if (event_wgt_SFs_btagging_branch) { event_wgt_SFs_btagging_branch->SetAddress(&event_wgt_SFs_btagging_); }
  }
  event_wgt_SFs_btagging_EffDn_branch = 0;
  if (tree->GetBranch("event_wgt_SFs_btagging_EffDn") != 0) {
    event_wgt_SFs_btagging_EffDn_branch = tree->GetBranch("event_wgt_SFs_btagging_EffDn");
    if (event_wgt_SFs_btagging_EffDn_branch) { event_wgt_SFs_btagging_EffDn_branch->SetAddress(&event_wgt_SFs_btagging_EffDn_); }
  }
  event_wgt_SFs_btagging_EffUp_branch = 0;
  if (tree->GetBranch("event_wgt_SFs_btagging_EffUp") != 0) {
    event_wgt_SFs_btagging_EffUp_branch = tree->GetBranch("event_wgt_SFs_btagging_EffUp");
    if (event_wgt_SFs_btagging_EffUp_branch) { event_wgt_SFs_btagging_EffUp_branch->SetAddress(&event_wgt_SFs_btagging_EffUp_); }
  }
  event_wgt_SFs_electrons_branch = 0;
  if (tree->GetBranch("event_wgt_SFs_electrons") != 0) {
    event_wgt_SFs_electrons_branch = tree->GetBranch("event_wgt_SFs_electrons");
    if (event_wgt_SFs_electrons_branch) { event_wgt_SFs_electrons_branch->SetAddress(&event_wgt_SFs_electrons_); }
  }
  event_wgt_SFs_electrons_AltMCDn_branch = 0;
  if (tree->GetBranch("event_wgt_SFs_electrons_AltMCDn") != 0) {
    event_wgt_SFs_electrons_AltMCDn_branch = tree->GetBranch("event_wgt_SFs_electrons_AltMCDn");
    if (event_wgt_SFs_electrons_AltMCDn_branch) { event_wgt_SFs_electrons_AltMCDn_branch->SetAddress(&event_wgt_SFs_electrons_AltMCDn_); }
  }
  event_wgt_SFs_electrons_AltMCUp_branch = 0;
  if (tree->GetBranch("event_wgt_SFs_electrons_AltMCUp") != 0) {
    event_wgt_SFs_electrons_AltMCUp_branch = tree->GetBranch("event_wgt_SFs_electrons_AltMCUp");
    if (event_wgt_SFs_electrons_AltMCUp_branch) { event_wgt_SFs_electrons_AltMCUp_branch->SetAddress(&event_wgt_SFs_electrons_AltMCUp_); }
  }
  event_wgt_SFs_electrons_StatDn_branch = 0;
  if (tree->GetBranch("event_wgt_SFs_electrons_StatDn") != 0) {
    event_wgt_SFs_electrons_StatDn_branch = tree->GetBranch("event_wgt_SFs_electrons_StatDn");
    if (event_wgt_SFs_electrons_StatDn_branch) { event_wgt_SFs_electrons_StatDn_branch->SetAddress(&event_wgt_SFs_electrons_StatDn_); }
  }
  event_wgt_SFs_electrons_StatUp_branch = 0;
  if (tree->GetBranch("event_wgt_SFs_electrons_StatUp") != 0) {
    event_wgt_SFs_electrons_StatUp_branch = tree->GetBranch("event_wgt_SFs_electrons_StatUp");
    if (event_wgt_SFs_electrons_StatUp_branch) { event_wgt_SFs_electrons_StatUp_branch->SetAddress(&event_wgt_SFs_electrons_StatUp_); }
  }
  event_wgt_SFs_electrons_SystDn_branch = 0;
  if (tree->GetBranch("event_wgt_SFs_electrons_SystDn") != 0) {
    event_wgt_SFs_electrons_SystDn_branch = tree->GetBranch("event_wgt_SFs_electrons_SystDn");
    if (event_wgt_SFs_electrons_SystDn_branch) { event_wgt_SFs_electrons_SystDn_branch->SetAddress(&event_wgt_SFs_electrons_SystDn_); }
  }
  event_wgt_SFs_electrons_SystUp_branch = 0;
  if (tree->GetBranch("event_wgt_SFs_electrons_SystUp") != 0) {
    event_wgt_SFs_electrons_SystUp_branch = tree->GetBranch("event_wgt_SFs_electrons_SystUp");
    if (event_wgt_SFs_electrons_SystUp_branch) { event_wgt_SFs_electrons_SystUp_branch->SetAddress(&event_wgt_SFs_electrons_SystUp_); }
  }
  event_wgt_SFs_muons_branch = 0;
  if (tree->GetBranch("event_wgt_SFs_muons") != 0) {
    event_wgt_SFs_muons_branch = tree->GetBranch("event_wgt_SFs_muons");
    if (event_wgt_SFs_muons_branch) { event_wgt_SFs_muons_branch->SetAddress(&event_wgt_SFs_muons_); }
  }
  event_wgt_SFs_muons_AltMCDn_branch = 0;
  if (tree->GetBranch("event_wgt_SFs_muons_AltMCDn") != 0) {
    event_wgt_SFs_muons_AltMCDn_branch = tree->GetBranch("event_wgt_SFs_muons_AltMCDn");
    if (event_wgt_SFs_muons_AltMCDn_branch) { event_wgt_SFs_muons_AltMCDn_branch->SetAddress(&event_wgt_SFs_muons_AltMCDn_); }
  }
  event_wgt_SFs_muons_AltMCUp_branch = 0;
  if (tree->GetBranch("event_wgt_SFs_muons_AltMCUp") != 0) {
    event_wgt_SFs_muons_AltMCUp_branch = tree->GetBranch("event_wgt_SFs_muons_AltMCUp");
    if (event_wgt_SFs_muons_AltMCUp_branch) { event_wgt_SFs_muons_AltMCUp_branch->SetAddress(&event_wgt_SFs_muons_AltMCUp_); }
  }
  event_wgt_SFs_muons_StatDn_branch = 0;
  if (tree->GetBranch("event_wgt_SFs_muons_StatDn") != 0) {
    event_wgt_SFs_muons_StatDn_branch = tree->GetBranch("event_wgt_SFs_muons_StatDn");
    if (event_wgt_SFs_muons_StatDn_branch) { event_wgt_SFs_muons_StatDn_branch->SetAddress(&event_wgt_SFs_muons_StatDn_); }
  }
  event_wgt_SFs_muons_StatUp_branch = 0;
  if (tree->GetBranch("event_wgt_SFs_muons_StatUp") != 0) {
    event_wgt_SFs_muons_StatUp_branch = tree->GetBranch("event_wgt_SFs_muons_StatUp");
    if (event_wgt_SFs_muons_StatUp_branch) { event_wgt_SFs_muons_StatUp_branch->SetAddress(&event_wgt_SFs_muons_StatUp_); }
  }
  event_wgt_SFs_muons_SystDn_branch = 0;
  if (tree->GetBranch("event_wgt_SFs_muons_SystDn") != 0) {
    event_wgt_SFs_muons_SystDn_branch = tree->GetBranch("event_wgt_SFs_muons_SystDn");
    if (event_wgt_SFs_muons_SystDn_branch) { event_wgt_SFs_muons_SystDn_branch->SetAddress(&event_wgt_SFs_muons_SystDn_); }
  }
  event_wgt_SFs_muons_SystUp_branch = 0;
  if (tree->GetBranch("event_wgt_SFs_muons_SystUp") != 0) {
    event_wgt_SFs_muons_SystUp_branch = tree->GetBranch("event_wgt_SFs_muons_SystUp");
    if (event_wgt_SFs_muons_SystUp_branch) { event_wgt_SFs_muons_SystUp_branch->SetAddress(&event_wgt_SFs_muons_SystUp_); }
  }
  event_wgt_SFs_photons_branch = 0;
  if (tree->GetBranch("event_wgt_SFs_photons") != 0) {
    event_wgt_SFs_photons_branch = tree->GetBranch("event_wgt_SFs_photons");
    if (event_wgt_SFs_photons_branch) { event_wgt_SFs_photons_branch->SetAddress(&event_wgt_SFs_photons_); }
  }
  event_wgt_SFs_photons_EffDn_branch = 0;
  if (tree->GetBranch("event_wgt_SFs_photons_EffDn") != 0) {
    event_wgt_SFs_photons_EffDn_branch = tree->GetBranch("event_wgt_SFs_photons_EffDn");
    if (event_wgt_SFs_photons_EffDn_branch) { event_wgt_SFs_photons_EffDn_branch->SetAddress(&event_wgt_SFs_photons_EffDn_); }
  }
  event_wgt_SFs_photons_EffUp_branch = 0;
  if (tree->GetBranch("event_wgt_SFs_photons_EffUp") != 0) {
    event_wgt_SFs_photons_EffUp_branch = tree->GetBranch("event_wgt_SFs_photons_EffUp");
    if (event_wgt_SFs_photons_EffUp_branch) { event_wgt_SFs_photons_EffUp_branch->SetAddress(&event_wgt_SFs_photons_EffUp_); }
  }
  event_wgt_adjustment_AsMZDn_branch = 0;
  if (tree->GetBranch("event_wgt_adjustment_AsMZDn") != 0) {
    event_wgt_adjustment_AsMZDn_branch = tree->GetBranch("event_wgt_adjustment_AsMZDn");
    if (event_wgt_adjustment_AsMZDn_branch) { event_wgt_adjustment_AsMZDn_branch->SetAddress(&event_wgt_adjustment_AsMZDn_); }
  }
  event_wgt_adjustment_AsMZUp_branch = 0;
  if (tree->GetBranch("event_wgt_adjustment_AsMZUp") != 0) {
    event_wgt_adjustment_AsMZUp_branch = tree->GetBranch("event_wgt_adjustment_AsMZUp");
    if (event_wgt_adjustment_AsMZUp_branch) { event_wgt_adjustment_AsMZUp_branch->SetAddress(&event_wgt_adjustment_AsMZUp_); }
  }
  event_wgt_adjustment_NNPDF30_branch = 0;
  if (tree->GetBranch("event_wgt_adjustment_NNPDF30") != 0) {
    event_wgt_adjustment_NNPDF30_branch = tree->GetBranch("event_wgt_adjustment_NNPDF30");
    if (event_wgt_adjustment_NNPDF30_branch) { event_wgt_adjustment_NNPDF30_branch->SetAddress(&event_wgt_adjustment_NNPDF30_); }
  }
  event_wgt_adjustment_NNPDF30_AsMZDn_branch = 0;
  if (tree->GetBranch("event_wgt_adjustment_NNPDF30_AsMZDn") != 0) {
    event_wgt_adjustment_NNPDF30_AsMZDn_branch = tree->GetBranch("event_wgt_adjustment_NNPDF30_AsMZDn");
    if (event_wgt_adjustment_NNPDF30_AsMZDn_branch) { event_wgt_adjustment_NNPDF30_AsMZDn_branch->SetAddress(&event_wgt_adjustment_NNPDF30_AsMZDn_); }
  }
  event_wgt_adjustment_NNPDF30_AsMZUp_branch = 0;
  if (tree->GetBranch("event_wgt_adjustment_NNPDF30_AsMZUp") != 0) {
    event_wgt_adjustment_NNPDF30_AsMZUp_branch = tree->GetBranch("event_wgt_adjustment_NNPDF30_AsMZUp");
    if (event_wgt_adjustment_NNPDF30_AsMZUp_branch) { event_wgt_adjustment_NNPDF30_AsMZUp_branch->SetAddress(&event_wgt_adjustment_NNPDF30_AsMZUp_); }
  }
  event_wgt_adjustment_NNPDF30_PDFReplicaDn_branch = 0;
  if (tree->GetBranch("event_wgt_adjustment_NNPDF30_PDFReplicaDn") != 0) {
    event_wgt_adjustment_NNPDF30_PDFReplicaDn_branch = tree->GetBranch("event_wgt_adjustment_NNPDF30_PDFReplicaDn");
    if (event_wgt_adjustment_NNPDF30_PDFReplicaDn_branch) { event_wgt_adjustment_NNPDF30_PDFReplicaDn_branch->SetAddress(&event_wgt_adjustment_NNPDF30_PDFReplicaDn_); }
  }
  event_wgt_adjustment_NNPDF30_PDFReplicaUp_branch = 0;
  if (tree->GetBranch("event_wgt_adjustment_NNPDF30_PDFReplicaUp") != 0) {
    event_wgt_adjustment_NNPDF30_PDFReplicaUp_branch = tree->GetBranch("event_wgt_adjustment_NNPDF30_PDFReplicaUp");
    if (event_wgt_adjustment_NNPDF30_PDFReplicaUp_branch) { event_wgt_adjustment_NNPDF30_PDFReplicaUp_branch->SetAddress(&event_wgt_adjustment_NNPDF30_PDFReplicaUp_); }
  }
  event_wgt_adjustment_PDFReplicaDn_branch = 0;
  if (tree->GetBranch("event_wgt_adjustment_PDFReplicaDn") != 0) {
    event_wgt_adjustment_PDFReplicaDn_branch = tree->GetBranch("event_wgt_adjustment_PDFReplicaDn");
    if (event_wgt_adjustment_PDFReplicaDn_branch) { event_wgt_adjustment_PDFReplicaDn_branch->SetAddress(&event_wgt_adjustment_PDFReplicaDn_); }
  }
  event_wgt_adjustment_PDFReplicaUp_branch = 0;
  if (tree->GetBranch("event_wgt_adjustment_PDFReplicaUp") != 0) {
    event_wgt_adjustment_PDFReplicaUp_branch = tree->GetBranch("event_wgt_adjustment_PDFReplicaUp");
    if (event_wgt_adjustment_PDFReplicaUp_branch) { event_wgt_adjustment_PDFReplicaUp_branch->SetAddress(&event_wgt_adjustment_PDFReplicaUp_); }
  }
  event_wgt_adjustment_PDFScaleDn_branch = 0;
  if (tree->GetBranch("event_wgt_adjustment_PDFScaleDn") != 0) {
    event_wgt_adjustment_PDFScaleDn_branch = tree->GetBranch("event_wgt_adjustment_PDFScaleDn");
    if (event_wgt_adjustment_PDFScaleDn_branch) { event_wgt_adjustment_PDFScaleDn_branch->SetAddress(&event_wgt_adjustment_PDFScaleDn_); }
  }
  event_wgt_adjustment_PDFScaleUp_branch = 0;
  if (tree->GetBranch("event_wgt_adjustment_PDFScaleUp") != 0) {
    event_wgt_adjustment_PDFScaleUp_branch = tree->GetBranch("event_wgt_adjustment_PDFScaleUp");
    if (event_wgt_adjustment_PDFScaleUp_branch) { event_wgt_adjustment_PDFScaleUp_branch->SetAddress(&event_wgt_adjustment_PDFScaleUp_); }
  }
  event_wgt_adjustment_PythiaScaleDn_branch = 0;
  if (tree->GetBranch("event_wgt_adjustment_PythiaScaleDn") != 0) {
    event_wgt_adjustment_PythiaScaleDn_branch = tree->GetBranch("event_wgt_adjustment_PythiaScaleDn");
    if (event_wgt_adjustment_PythiaScaleDn_branch) { event_wgt_adjustment_PythiaScaleDn_branch->SetAddress(&event_wgt_adjustment_PythiaScaleDn_); }
  }
  event_wgt_adjustment_PythiaScaleUp_branch = 0;
  if (tree->GetBranch("event_wgt_adjustment_PythiaScaleUp") != 0) {
    event_wgt_adjustment_PythiaScaleUp_branch = tree->GetBranch("event_wgt_adjustment_PythiaScaleUp");
    if (event_wgt_adjustment_PythiaScaleUp_branch) { event_wgt_adjustment_PythiaScaleUp_branch->SetAddress(&event_wgt_adjustment_PythiaScaleUp_); }
  }
  event_wgt_adjustment_QCDScaleDn_branch = 0;
  if (tree->GetBranch("event_wgt_adjustment_QCDScaleDn") != 0) {
    event_wgt_adjustment_QCDScaleDn_branch = tree->GetBranch("event_wgt_adjustment_QCDScaleDn");
    if (event_wgt_adjustment_QCDScaleDn_branch) { event_wgt_adjustment_QCDScaleDn_branch->SetAddress(&event_wgt_adjustment_QCDScaleDn_); }
  }
  event_wgt_adjustment_QCDScaleUp_branch = 0;
  if (tree->GetBranch("event_wgt_adjustment_QCDScaleUp") != 0) {
    event_wgt_adjustment_QCDScaleUp_branch = tree->GetBranch("event_wgt_adjustment_QCDScaleUp");
    if (event_wgt_adjustment_QCDScaleUp_branch) { event_wgt_adjustment_QCDScaleUp_branch->SetAddress(&event_wgt_adjustment_QCDScaleUp_); }
  }
  event_wgt_triggers_branch = 0;
  if (tree->GetBranch("event_wgt_triggers") != 0) {
    event_wgt_triggers_branch = tree->GetBranch("event_wgt_triggers");
    if (event_wgt_triggers_branch) { event_wgt_triggers_branch->SetAddress(&event_wgt_triggers_); }
  }
  event_wgt_triggers_Dilepton_branch = 0;
  if (tree->GetBranch("event_wgt_triggers_Dilepton") != 0) {
    event_wgt_triggers_Dilepton_branch = tree->GetBranch("event_wgt_triggers_Dilepton");
    if (event_wgt_triggers_Dilepton_branch) { event_wgt_triggers_Dilepton_branch->SetAddress(&event_wgt_triggers_Dilepton_); }
  }
  event_wgt_triggers_Dilepton_DF_Extra_branch = 0;
  if (tree->GetBranch("event_wgt_triggers_Dilepton_DF_Extra") != 0) {
    event_wgt_triggers_Dilepton_DF_Extra_branch = tree->GetBranch("event_wgt_triggers_Dilepton_DF_Extra");
    if (event_wgt_triggers_Dilepton_DF_Extra_branch) { event_wgt_triggers_Dilepton_DF_Extra_branch->SetAddress(&event_wgt_triggers_Dilepton_DF_Extra_); }
  }
  event_wgt_triggers_PFHT_Control_branch = 0;
  if (tree->GetBranch("event_wgt_triggers_PFHT_Control") != 0) {
    event_wgt_triggers_PFHT_Control_branch = tree->GetBranch("event_wgt_triggers_PFHT_Control");
    if (event_wgt_triggers_PFHT_Control_branch) { event_wgt_triggers_PFHT_Control_branch->SetAddress(&event_wgt_triggers_PFHT_Control_); }
  }
  event_wgt_triggers_PFMET_MHT_Control_branch = 0;
  if (tree->GetBranch("event_wgt_triggers_PFMET_MHT_Control") != 0) {
    event_wgt_triggers_PFMET_MHT_Control_branch = tree->GetBranch("event_wgt_triggers_PFMET_MHT_Control");
    if (event_wgt_triggers_PFMET_MHT_Control_branch) { event_wgt_triggers_PFMET_MHT_Control_branch->SetAddress(&event_wgt_triggers_PFMET_MHT_Control_); }
  }
  event_wgt_triggers_SingleLepton_branch = 0;
  if (tree->GetBranch("event_wgt_triggers_SingleLepton") != 0) {
    event_wgt_triggers_SingleLepton_branch = tree->GetBranch("event_wgt_triggers_SingleLepton");
    if (event_wgt_triggers_SingleLepton_branch) { event_wgt_triggers_SingleLepton_branch->SetAddress(&event_wgt_triggers_SingleLepton_); }
  }
  event_wgt_triggers_SinglePhoton_branch = 0;
  if (tree->GetBranch("event_wgt_triggers_SinglePhoton") != 0) {
    event_wgt_triggers_SinglePhoton_branch = tree->GetBranch("event_wgt_triggers_SinglePhoton");
    if (event_wgt_triggers_SinglePhoton_branch) { event_wgt_triggers_SinglePhoton_branch->SetAddress(&event_wgt_triggers_SinglePhoton_); }
  }
  genak4jets_eta_branch = 0;
  if (tree->GetBranch("genak4jets_eta") != 0) {
    genak4jets_eta_branch = tree->GetBranch("genak4jets_eta");
    if (genak4jets_eta_branch) { genak4jets_eta_branch->SetAddress(&genak4jets_eta_); }
  }
  genak4jets_mass_branch = 0;
  if (tree->GetBranch("genak4jets_mass") != 0) {
    genak4jets_mass_branch = tree->GetBranch("genak4jets_mass");
    if (genak4jets_mass_branch) { genak4jets_mass_branch->SetAddress(&genak4jets_mass_); }
  }
  genak4jets_phi_branch = 0;
  if (tree->GetBranch("genak4jets_phi") != 0) {
    genak4jets_phi_branch = tree->GetBranch("genak4jets_phi");
    if (genak4jets_phi_branch) { genak4jets_phi_branch->SetAddress(&genak4jets_phi_); }
  }
  genak4jets_pt_branch = 0;
  if (tree->GetBranch("genak4jets_pt") != 0) {
    genak4jets_pt_branch = tree->GetBranch("genak4jets_pt");
    if (genak4jets_pt_branch) { genak4jets_pt_branch->SetAddress(&genak4jets_pt_); }
  }
  genmet_pTmiss_branch = 0;
  if (tree->GetBranch("genmet_pTmiss") != 0) {
    genmet_pTmiss_branch = tree->GetBranch("genmet_pTmiss");
    if (genmet_pTmiss_branch) { genmet_pTmiss_branch->SetAddress(&genmet_pTmiss_); }
  }
  genmet_phimiss_branch = 0;
  if (tree->GetBranch("genmet_phimiss") != 0) {
    genmet_phimiss_branch = tree->GetBranch("genmet_phimiss");
    if (genmet_phimiss_branch) { genmet_phimiss_branch->SetAddress(&genmet_phimiss_); }
  }
  lepton_eff_branch = 0;
  if (tree->GetBranch("lepton_eff") != 0) {
    lepton_eff_branch = tree->GetBranch("lepton_eff");
    if (lepton_eff_branch) { lepton_eff_branch->SetAddress(&lepton_eff_); }
  }
  lepton_eff_StatDn_branch = 0;
  if (tree->GetBranch("lepton_eff_StatDn") != 0) {
    lepton_eff_StatDn_branch = tree->GetBranch("lepton_eff_StatDn");
    if (lepton_eff_StatDn_branch) { lepton_eff_StatDn_branch->SetAddress(&lepton_eff_StatDn_); }
  }
  lepton_eff_StatUp_branch = 0;
  if (tree->GetBranch("lepton_eff_StatUp") != 0) {
    lepton_eff_StatUp_branch = tree->GetBranch("lepton_eff_StatUp");
    if (lepton_eff_StatUp_branch) { lepton_eff_StatUp_branch->SetAddress(&lepton_eff_StatUp_); }
  }
  lepton_eff_SystDn_branch = 0;
  if (tree->GetBranch("lepton_eff_SystDn") != 0) {
    lepton_eff_SystDn_branch = tree->GetBranch("lepton_eff_SystDn");
    if (lepton_eff_SystDn_branch) { lepton_eff_SystDn_branch->SetAddress(&lepton_eff_SystDn_); }
  }
  lepton_eff_SystUp_branch = 0;
  if (tree->GetBranch("lepton_eff_SystUp") != 0) {
    lepton_eff_SystUp_branch = tree->GetBranch("lepton_eff_SystUp");
    if (lepton_eff_SystUp_branch) { lepton_eff_SystUp_branch->SetAddress(&lepton_eff_SystUp_); }
  }
  lepton_eta_branch = 0;
  if (tree->GetBranch("lepton_eta") != 0) {
    lepton_eta_branch = tree->GetBranch("lepton_eta");
    if (lepton_eta_branch) { lepton_eta_branch->SetAddress(&lepton_eta_); }
  }
  lepton_id_branch = 0;
  if (tree->GetBranch("lepton_id") != 0) {
    lepton_id_branch = tree->GetBranch("lepton_id");
    if (lepton_id_branch) { lepton_id_branch->SetAddress(&lepton_id_); }
  }
  lepton_is_genMatched_prompt_branch = 0;
  if (tree->GetBranch("lepton_is_genMatched_prompt") != 0) {
    lepton_is_genMatched_prompt_branch = tree->GetBranch("lepton_is_genMatched_prompt");
    if (lepton_is_genMatched_prompt_branch) { lepton_is_genMatched_prompt_branch->SetAddress(&lepton_is_genMatched_prompt_); }
  }
  lepton_mass_branch = 0;
  if (tree->GetBranch("lepton_mass") != 0) {
    lepton_mass_branch = tree->GetBranch("lepton_mass");
    if (lepton_mass_branch) { lepton_mass_branch->SetAddress(&lepton_mass_); }
  }
  lepton_phi_branch = 0;
  if (tree->GetBranch("lepton_phi") != 0) {
    lepton_phi_branch = tree->GetBranch("lepton_phi");
    if (lepton_phi_branch) { lepton_phi_branch->SetAddress(&lepton_phi_); }
  }
  lepton_pt_branch = 0;
  if (tree->GetBranch("lepton_pt") != 0) {
    lepton_pt_branch = tree->GetBranch("lepton_pt");
    if (lepton_pt_branch) { lepton_pt_branch->SetAddress(&lepton_pt_); }
  }
  genpromptparticles_sump4_mass_branch = 0;
  if (tree->GetBranch("genpromptparticles_sump4_mass") != 0) {
    genpromptparticles_sump4_mass_branch = tree->GetBranch("genpromptparticles_sump4_mass");
    if (genpromptparticles_sump4_mass_branch) { genpromptparticles_sump4_mass_branch->SetAddress(&genpromptparticles_sump4_mass_); }
  }
  genpromptparticles_sump4_pt_branch = 0;
  if (tree->GetBranch("genpromptparticles_sump4_pt") != 0) {
    genpromptparticles_sump4_pt_branch = tree->GetBranch("genpromptparticles_sump4_pt");
    if (genpromptparticles_sump4_pt_branch) { genpromptparticles_sump4_pt_branch->SetAddress(&genpromptparticles_sump4_pt_); }
  }
  genpromptparticles_sump4_rapidity_branch = 0;
  if (tree->GetBranch("genpromptparticles_sump4_rapidity") != 0) {
    genpromptparticles_sump4_rapidity_branch = tree->GetBranch("genpromptparticles_sump4_rapidity");
    if (genpromptparticles_sump4_rapidity_branch) { genpromptparticles_sump4_rapidity_branch->SetAddress(&genpromptparticles_sump4_rapidity_); }
  }
  leptons_eff_branch = 0;
  if (tree->GetBranch("leptons_eff") != 0) {
    leptons_eff_branch = tree->GetBranch("leptons_eff");
    if (leptons_eff_branch) { leptons_eff_branch->SetAddress(&leptons_eff_); }
  }
  leptons_eff_DF_branch = 0;
  if (tree->GetBranch("leptons_eff_DF") != 0) {
    leptons_eff_DF_branch = tree->GetBranch("leptons_eff_DF");
    if (leptons_eff_DF_branch) { leptons_eff_DF_branch->SetAddress(&leptons_eff_DF_); }
  }
  leptons_eff_DF_StatDn_branch = 0;
  if (tree->GetBranch("leptons_eff_DF_StatDn") != 0) {
    leptons_eff_DF_StatDn_branch = tree->GetBranch("leptons_eff_DF_StatDn");
    if (leptons_eff_DF_StatDn_branch) { leptons_eff_DF_StatDn_branch->SetAddress(&leptons_eff_DF_StatDn_); }
  }
  leptons_eff_DF_StatUp_branch = 0;
  if (tree->GetBranch("leptons_eff_DF_StatUp") != 0) {
    leptons_eff_DF_StatUp_branch = tree->GetBranch("leptons_eff_DF_StatUp");
    if (leptons_eff_DF_StatUp_branch) { leptons_eff_DF_StatUp_branch->SetAddress(&leptons_eff_DF_StatUp_); }
  }
  leptons_eff_DF_SystDn_branch = 0;
  if (tree->GetBranch("leptons_eff_DF_SystDn") != 0) {
    leptons_eff_DF_SystDn_branch = tree->GetBranch("leptons_eff_DF_SystDn");
    if (leptons_eff_DF_SystDn_branch) { leptons_eff_DF_SystDn_branch->SetAddress(&leptons_eff_DF_SystDn_); }
  }
  leptons_eff_DF_SystUp_branch = 0;
  if (tree->GetBranch("leptons_eff_DF_SystUp") != 0) {
    leptons_eff_DF_SystUp_branch = tree->GetBranch("leptons_eff_DF_SystUp");
    if (leptons_eff_DF_SystUp_branch) { leptons_eff_DF_SystUp_branch->SetAddress(&leptons_eff_DF_SystUp_); }
  }
  leptons_eff_StatDn_branch = 0;
  if (tree->GetBranch("leptons_eff_StatDn") != 0) {
    leptons_eff_StatDn_branch = tree->GetBranch("leptons_eff_StatDn");
    if (leptons_eff_StatDn_branch) { leptons_eff_StatDn_branch->SetAddress(&leptons_eff_StatDn_); }
  }
  leptons_eff_StatUp_branch = 0;
  if (tree->GetBranch("leptons_eff_StatUp") != 0) {
    leptons_eff_StatUp_branch = tree->GetBranch("leptons_eff_StatUp");
    if (leptons_eff_StatUp_branch) { leptons_eff_StatUp_branch->SetAddress(&leptons_eff_StatUp_); }
  }
  leptons_eff_SystDn_branch = 0;
  if (tree->GetBranch("leptons_eff_SystDn") != 0) {
    leptons_eff_SystDn_branch = tree->GetBranch("leptons_eff_SystDn");
    if (leptons_eff_SystDn_branch) { leptons_eff_SystDn_branch->SetAddress(&leptons_eff_SystDn_); }
  }
  leptons_eff_SystUp_branch = 0;
  if (tree->GetBranch("leptons_eff_SystUp") != 0) {
    leptons_eff_SystUp_branch = tree->GetBranch("leptons_eff_SystUp");
    if (leptons_eff_SystUp_branch) { leptons_eff_SystUp_branch->SetAddress(&leptons_eff_SystUp_); }
  }
  leptons_eta_branch = 0;
  if (tree->GetBranch("leptons_eta") != 0) {
    leptons_eta_branch = tree->GetBranch("leptons_eta");
    if (leptons_eta_branch) { leptons_eta_branch->SetAddress(&leptons_eta_); }
  }
  leptons_id_branch = 0;
  if (tree->GetBranch("leptons_id") != 0) {
    leptons_id_branch = tree->GetBranch("leptons_id");
    if (leptons_id_branch) { leptons_id_branch->SetAddress(&leptons_id_); }
  }
  leptons_is_TOmatched_SingleLepton_branch = 0;
  if (tree->GetBranch("leptons_is_TOmatched_SingleLepton") != 0) {
    leptons_is_TOmatched_SingleLepton_branch = tree->GetBranch("leptons_is_TOmatched_SingleLepton");
    if (leptons_is_TOmatched_SingleLepton_branch) { leptons_is_TOmatched_SingleLepton_branch->SetAddress(&leptons_is_TOmatched_SingleLepton_); }
  }
  leptons_is_genMatched_prompt_branch = 0;
  if (tree->GetBranch("leptons_is_genMatched_prompt") != 0) {
    leptons_is_genMatched_prompt_branch = tree->GetBranch("leptons_is_genMatched_prompt");
    if (leptons_is_genMatched_prompt_branch) { leptons_is_genMatched_prompt_branch->SetAddress(&leptons_is_genMatched_prompt_); }
  }
  leptons_mass_branch = 0;
  if (tree->GetBranch("leptons_mass") != 0) {
    leptons_mass_branch = tree->GetBranch("leptons_mass");
    if (leptons_mass_branch) { leptons_mass_branch->SetAddress(&leptons_mass_); }
  }
  leptons_phi_branch = 0;
  if (tree->GetBranch("leptons_phi") != 0) {
    leptons_phi_branch = tree->GetBranch("leptons_phi");
    if (leptons_phi_branch) { leptons_phi_branch->SetAddress(&leptons_phi_); }
  }
  leptons_pt_branch = 0;
  if (tree->GetBranch("leptons_pt") != 0) {
    leptons_pt_branch = tree->GetBranch("leptons_pt");
    if (leptons_pt_branch) { leptons_pt_branch->SetAddress(&leptons_pt_); }
  }
  min_abs_dPhi_pTj_pTmiss_branch = 0;
  if (tree->GetBranch("min_abs_dPhi_pTj_pTmiss") != 0) {
    min_abs_dPhi_pTj_pTmiss_branch = tree->GetBranch("min_abs_dPhi_pTj_pTmiss");
    if (min_abs_dPhi_pTj_pTmiss_branch) { min_abs_dPhi_pTj_pTmiss_branch->SetAddress(&min_abs_dPhi_pTj_pTmiss_); }
  }
  pConst_JJQCD_SIG_ghg2_1_JHUGen_branch = 0;
  if (tree->GetBranch("pConst_JJQCD_SIG_ghg2_1_JHUGen") != 0) {
    pConst_JJQCD_SIG_ghg2_1_JHUGen_branch = tree->GetBranch("pConst_JJQCD_SIG_ghg2_1_JHUGen");
    if (pConst_JJQCD_SIG_ghg2_1_JHUGen_branch) { pConst_JJQCD_SIG_ghg2_1_JHUGen_branch->SetAddress(&pConst_JJQCD_SIG_ghg2_1_JHUGen_); }
  }
  pConst_JJVBF_SIG_ghv1_1_JHUGen_branch = 0;
  if (tree->GetBranch("pConst_JJVBF_SIG_ghv1_1_JHUGen") != 0) {
    pConst_JJVBF_SIG_ghv1_1_JHUGen_branch = tree->GetBranch("pConst_JJVBF_SIG_ghv1_1_JHUGen");
    if (pConst_JJVBF_SIG_ghv1_1_JHUGen_branch) { pConst_JJVBF_SIG_ghv1_1_JHUGen_branch->SetAddress(&pConst_JJVBF_SIG_ghv1_1_JHUGen_); }
  }
  p_JJQCD_SIG_ghg2_1_JHUGen_branch = 0;
  if (tree->GetBranch("p_JJQCD_SIG_ghg2_1_JHUGen") != 0) {
    p_JJQCD_SIG_ghg2_1_JHUGen_branch = tree->GetBranch("p_JJQCD_SIG_ghg2_1_JHUGen");
    if (p_JJQCD_SIG_ghg2_1_JHUGen_branch) { p_JJQCD_SIG_ghg2_1_JHUGen_branch->SetAddress(&p_JJQCD_SIG_ghg2_1_JHUGen_); }
  }
  p_JJVBF_SIG_ghv1_1_JHUGen_branch = 0;
  if (tree->GetBranch("p_JJVBF_SIG_ghv1_1_JHUGen") != 0) {
    p_JJVBF_SIG_ghv1_1_JHUGen_branch = tree->GetBranch("p_JJVBF_SIG_ghv1_1_JHUGen");
    if (p_JJVBF_SIG_ghv1_1_JHUGen_branch) { p_JJVBF_SIG_ghv1_1_JHUGen_branch->SetAddress(&p_JJVBF_SIG_ghv1_1_JHUGen_); }
  }
  p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_branch = 0;
  if (tree->GetBranch("p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen") != 0) {
    p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_branch = tree->GetBranch("p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen");
    if (p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_branch) { p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_branch->SetAddress(&p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_); }
  }
  p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_branch = 0;
  if (tree->GetBranch("p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen") != 0) {
    p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_branch = tree->GetBranch("p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen");
    if (p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_branch) { p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_branch->SetAddress(&p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_); }
  }
  p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_branch = 0;
  if (tree->GetBranch("p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen") != 0) {
    p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_branch = tree->GetBranch("p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen");
    if (p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_branch) { p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_branch->SetAddress(&p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_); }
  }
  p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_branch = 0;
  if (tree->GetBranch("p_JJVBF_SIG_ghv1prime2_1E4_JHUGen") != 0) {
    p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_branch = tree->GetBranch("p_JJVBF_SIG_ghv1prime2_1E4_JHUGen");
    if (p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_branch) { p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_branch->SetAddress(&p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_); }
  }
  p_JJVBF_SIG_ghv2_1_JHUGen_branch = 0;
  if (tree->GetBranch("p_JJVBF_SIG_ghv2_1_JHUGen") != 0) {
    p_JJVBF_SIG_ghv2_1_JHUGen_branch = tree->GetBranch("p_JJVBF_SIG_ghv2_1_JHUGen");
    if (p_JJVBF_SIG_ghv2_1_JHUGen_branch) { p_JJVBF_SIG_ghv2_1_JHUGen_branch->SetAddress(&p_JJVBF_SIG_ghv2_1_JHUGen_); }
  }
  p_JJVBF_SIG_ghv4_1_JHUGen_branch = 0;
  if (tree->GetBranch("p_JJVBF_SIG_ghv4_1_JHUGen") != 0) {
    p_JJVBF_SIG_ghv4_1_JHUGen_branch = tree->GetBranch("p_JJVBF_SIG_ghv4_1_JHUGen");
    if (p_JJVBF_SIG_ghv4_1_JHUGen_branch) { p_JJVBF_SIG_ghv4_1_JHUGen_branch->SetAddress(&p_JJVBF_SIG_ghv4_1_JHUGen_); }
  }
  p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_branch = 0;
  if (tree->GetBranch("p_JJVBF_SIG_ghza1prime2_1E4_JHUGen") != 0) {
    p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_branch = tree->GetBranch("p_JJVBF_SIG_ghza1prime2_1E4_JHUGen");
    if (p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_branch) { p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_branch->SetAddress(&p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_); }
  }
  photon_MIPTotalEnergy_branch = 0;
  if (tree->GetBranch("photon_MIPTotalEnergy") != 0) {
    photon_MIPTotalEnergy_branch = tree->GetBranch("photon_MIPTotalEnergy");
    if (photon_MIPTotalEnergy_branch) { photon_MIPTotalEnergy_branch->SetAddress(&photon_MIPTotalEnergy_); }
  }
  photon_eta_branch = 0;
  if (tree->GetBranch("photon_eta") != 0) {
    photon_eta_branch = tree->GetBranch("photon_eta");
    if (photon_eta_branch) { photon_eta_branch->SetAddress(&photon_eta_); }
  }
  photon_full5x5_r9_branch = 0;
  if (tree->GetBranch("photon_full5x5_r9") != 0) {
    photon_full5x5_r9_branch = tree->GetBranch("photon_full5x5_r9");
    if (photon_full5x5_r9_branch) { photon_full5x5_r9_branch->SetAddress(&photon_full5x5_r9_); }
  }
  photon_full5x5_sigmaIEtaIEta_branch = 0;
  if (tree->GetBranch("photon_full5x5_sigmaIEtaIEta") != 0) {
    photon_full5x5_sigmaIEtaIEta_branch = tree->GetBranch("photon_full5x5_sigmaIEtaIEta");
    if (photon_full5x5_sigmaIEtaIEta_branch) { photon_full5x5_sigmaIEtaIEta_branch->SetAddress(&photon_full5x5_sigmaIEtaIEta_); }
  }
  photon_full5x5_sigmaIPhiIPhi_branch = 0;
  if (tree->GetBranch("photon_full5x5_sigmaIPhiIPhi") != 0) {
    photon_full5x5_sigmaIPhiIPhi_branch = tree->GetBranch("photon_full5x5_sigmaIPhiIPhi");
    if (photon_full5x5_sigmaIPhiIPhi_branch) { photon_full5x5_sigmaIPhiIPhi_branch->SetAddress(&photon_full5x5_sigmaIPhiIPhi_); }
  }
  photon_isEB_branch = 0;
  if (tree->GetBranch("photon_isEB") != 0) {
    photon_isEB_branch = tree->GetBranch("photon_isEB");
    if (photon_isEB_branch) { photon_isEB_branch->SetAddress(&photon_isEB_); }
  }
  photon_isEBEEGap_branch = 0;
  if (tree->GetBranch("photon_isEBEEGap") != 0) {
    photon_isEBEEGap_branch = tree->GetBranch("photon_isEBEEGap");
    if (photon_isEBEEGap_branch) { photon_isEBEEGap_branch->SetAddress(&photon_isEBEEGap_); }
  }
  photon_isEE_branch = 0;
  if (tree->GetBranch("photon_isEE") != 0) {
    photon_isEE_branch = tree->GetBranch("photon_isEE");
    if (photon_isEE_branch) { photon_isEE_branch->SetAddress(&photon_isEE_); }
  }
  photon_isGap_branch = 0;
  if (tree->GetBranch("photon_isGap") != 0) {
    photon_isGap_branch = tree->GetBranch("photon_isGap");
    if (photon_isGap_branch) { photon_isGap_branch->SetAddress(&photon_isGap_); }
  }
  photon_is_METSafe_branch = 0;
  if (tree->GetBranch("photon_is_METSafe") != 0) {
    photon_is_METSafe_branch = tree->GetBranch("photon_is_METSafe");
    if (photon_is_METSafe_branch) { photon_is_METSafe_branch->SetAddress(&photon_is_METSafe_); }
  }
  photon_is_PFID_branch = 0;
  if (tree->GetBranch("photon_is_PFID") != 0) {
    photon_is_PFID_branch = tree->GetBranch("photon_is_PFID");
    if (photon_is_PFID_branch) { photon_is_PFID_branch->SetAddress(&photon_is_PFID_); }
  }
  photon_is_beamHaloSafe_branch = 0;
  if (tree->GetBranch("photon_is_beamHaloSafe") != 0) {
    photon_is_beamHaloSafe_branch = tree->GetBranch("photon_is_beamHaloSafe");
    if (photon_is_beamHaloSafe_branch) { photon_is_beamHaloSafe_branch->SetAddress(&photon_is_beamHaloSafe_); }
  }
  photon_is_conversionSafe_branch = 0;
  if (tree->GetBranch("photon_is_conversionSafe") != 0) {
    photon_is_conversionSafe_branch = tree->GetBranch("photon_is_conversionSafe");
    if (photon_is_conversionSafe_branch) { photon_is_conversionSafe_branch->SetAddress(&photon_is_conversionSafe_); }
  }
  photon_is_genMatched_prompt_branch = 0;
  if (tree->GetBranch("photon_is_genMatched_prompt") != 0) {
    photon_is_genMatched_prompt_branch = tree->GetBranch("photon_is_genMatched_prompt");
    if (photon_is_genMatched_prompt_branch) { photon_is_genMatched_prompt_branch->SetAddress(&photon_is_genMatched_prompt_); }
  }
  photon_is_inTime_branch = 0;
  if (tree->GetBranch("photon_is_inTime") != 0) {
    photon_is_inTime_branch = tree->GetBranch("photon_is_inTime");
    if (photon_is_inTime_branch) { photon_is_inTime_branch->SetAddress(&photon_is_inTime_); }
  }
  photon_is_spikeSafe_branch = 0;
  if (tree->GetBranch("photon_is_spikeSafe") != 0) {
    photon_is_spikeSafe_branch = tree->GetBranch("photon_is_spikeSafe");
    if (photon_is_spikeSafe_branch) { photon_is_spikeSafe_branch->SetAddress(&photon_is_spikeSafe_); }
  }
  photon_mass_branch = 0;
  if (tree->GetBranch("photon_mass") != 0) {
    photon_mass_branch = tree->GetBranch("photon_mass");
    if (photon_mass_branch) { photon_mass_branch->SetAddress(&photon_mass_); }
  }
  photon_pass_HGGSelection_branch = 0;
  if (tree->GetBranch("photon_pass_HGGSelection") != 0) {
    photon_pass_HGGSelection_branch = tree->GetBranch("photon_pass_HGGSelection");
    if (photon_pass_HGGSelection_branch) { photon_pass_HGGSelection_branch->SetAddress(&photon_pass_HGGSelection_); }
  }
  photon_phi_branch = 0;
  if (tree->GetBranch("photon_phi") != 0) {
    photon_phi_branch = tree->GetBranch("photon_phi");
    if (photon_phi_branch) { photon_phi_branch->SetAddress(&photon_phi_); }
  }
  photon_pt_branch = 0;
  if (tree->GetBranch("photon_pt") != 0) {
    photon_pt_branch = tree->GetBranch("photon_pt");
    if (photon_pt_branch) { photon_pt_branch->SetAddress(&photon_pt_); }
  }
  photon_seedTime_branch = 0;
  if (tree->GetBranch("photon_seedTime") != 0) {
    photon_seedTime_branch = tree->GetBranch("photon_seedTime");
    if (photon_seedTime_branch) { photon_seedTime_branch->SetAddress(&photon_seedTime_); }
  }
  tree->SetMakeClass(0);
}
void SkimTree::GetEntry(unsigned int idx) {
  index = idx;
  EventNumber_isLoaded = false;
  KFactor_EW_NLO_qqVV_Bkg_EWDn_isLoaded = false;
  KFactor_EW_NLO_qqVV_Bkg_EWUp_isLoaded = false;
  KFactor_EW_NLO_qqVV_Bkg_Nominal_isLoaded = false;
  KFactor_EW_NLO_qqVV_Bkg_arg_mass_isLoaded = false;
  KFactor_EW_NLO_qqVV_Bkg_arg_pthat_isLoaded = false;
  KFactor_EW_NLO_qqVV_Bkg_arg_rho_isLoaded = false;
  KFactor_EW_NLO_qqVV_Bkg_arg_that_isLoaded = false;
  KFactor_QCD_NNLO_qqVV_Bkg_Nominal_isLoaded = false;
  KFactor_QCD_NNLO_qqVV_Bkg_arg_mass_isLoaded = false;
  KFactor_QCD_NNLO_qqVV_Bkg_arg_pthat_isLoaded = false;
  LuminosityBlock_isLoaded = false;
  RunNumber_isLoaded = false;
  ak4jets_CEMF_isLoaded = false;
  ak4jets_HT_isLoaded = false;
  ak4jets_MHT_isLoaded = false;
  ak4jets_NEMF_isLoaded = false;
  ak4jets_btagWP_Bits_isLoaded = false;
  ak4jets_eta_isLoaded = false;
  ak4jets_is_genMatched_isLoaded = false;
  ak4jets_is_genMatched_fullCone_isLoaded = false;
  ak4jets_mass_isLoaded = false;
  ak4jets_phi_isLoaded = false;
  ak4jets_pt_isLoaded = false;
  ak8jets_eta_isLoaded = false;
  ak8jets_mass_isLoaded = false;
  ak8jets_phi_isLoaded = false;
  ak8jets_pt_isLoaded = false;
  dPhi_pTboson_pTmiss_isLoaded = false;
  dPhi_pTbosonjets_pTmiss_isLoaded = false;
  dilepton_eta_isLoaded = false;
  dilepton_id_isLoaded = false;
  dilepton_mass_isLoaded = false;
  dilepton_phi_isLoaded = false;
  dilepton_pt_isLoaded = false;
  electron_full5x5_r9_isLoaded = false;
  electron_full5x5_sigmaIEtaIEta_isLoaded = false;
  electron_full5x5_sigmaIPhiIPhi_isLoaded = false;
  electron_seedTime_isLoaded = false;
  electrons_full5x5_r9_isLoaded = false;
  electrons_full5x5_sigmaIEtaIEta_isLoaded = false;
  electrons_full5x5_sigmaIPhiIPhi_isLoaded = false;
  electrons_seedTime_isLoaded = false;
  event_mTZZ_isLoaded = false;
  event_mZZ_isLoaded = false;
  event_mlg_isLoaded = false;
  event_mllg_isLoaded = false;
  event_n_ak4jets_pt20_isLoaded = false;
  event_n_ak4jets_pt20_btagged_loose_isLoaded = false;
  event_n_ak4jets_pt20_btagged_medium_isLoaded = false;
  event_n_ak4jets_pt30_isLoaded = false;
  event_n_ak4jets_pt30_btagged_loose_isLoaded = false;
  event_n_ak4jets_pt30_btagged_medium_isLoaded = false;
  event_n_leptons_fakeableBase_isLoaded = false;
  event_n_vtxs_good_isLoaded = false;
  event_pTmiss_isLoaded = false;
  event_pass_tightMETFilters_isLoaded = false;
  event_phimiss_isLoaded = false;
  event_wgt_isLoaded = false;
  event_wgt_L1PrefiringDn_isLoaded = false;
  event_wgt_L1PrefiringUp_isLoaded = false;
  event_wgt_PUDn_isLoaded = false;
  event_wgt_PUUp_isLoaded = false;
  event_wgt_SFs_PUJetId_isLoaded = false;
  event_wgt_SFs_PUJetId_EffDn_isLoaded = false;
  event_wgt_SFs_PUJetId_EffUp_isLoaded = false;
  event_wgt_SFs_btagging_isLoaded = false;
  event_wgt_SFs_btagging_EffDn_isLoaded = false;
  event_wgt_SFs_btagging_EffUp_isLoaded = false;
  event_wgt_SFs_electrons_isLoaded = false;
  event_wgt_SFs_electrons_AltMCDn_isLoaded = false;
  event_wgt_SFs_electrons_AltMCUp_isLoaded = false;
  event_wgt_SFs_electrons_StatDn_isLoaded = false;
  event_wgt_SFs_electrons_StatUp_isLoaded = false;
  event_wgt_SFs_electrons_SystDn_isLoaded = false;
  event_wgt_SFs_electrons_SystUp_isLoaded = false;
  event_wgt_SFs_muons_isLoaded = false;
  event_wgt_SFs_muons_AltMCDn_isLoaded = false;
  event_wgt_SFs_muons_AltMCUp_isLoaded = false;
  event_wgt_SFs_muons_StatDn_isLoaded = false;
  event_wgt_SFs_muons_StatUp_isLoaded = false;
  event_wgt_SFs_muons_SystDn_isLoaded = false;
  event_wgt_SFs_muons_SystUp_isLoaded = false;
  event_wgt_SFs_photons_isLoaded = false;
  event_wgt_SFs_photons_EffDn_isLoaded = false;
  event_wgt_SFs_photons_EffUp_isLoaded = false;
  event_wgt_adjustment_AsMZDn_isLoaded = false;
  event_wgt_adjustment_AsMZUp_isLoaded = false;
  event_wgt_adjustment_NNPDF30_isLoaded = false;
  event_wgt_adjustment_NNPDF30_AsMZDn_isLoaded = false;
  event_wgt_adjustment_NNPDF30_AsMZUp_isLoaded = false;
  event_wgt_adjustment_NNPDF30_PDFReplicaDn_isLoaded = false;
  event_wgt_adjustment_NNPDF30_PDFReplicaUp_isLoaded = false;
  event_wgt_adjustment_PDFReplicaDn_isLoaded = false;
  event_wgt_adjustment_PDFReplicaUp_isLoaded = false;
  event_wgt_adjustment_PDFScaleDn_isLoaded = false;
  event_wgt_adjustment_PDFScaleUp_isLoaded = false;
  event_wgt_adjustment_PythiaScaleDn_isLoaded = false;
  event_wgt_adjustment_PythiaScaleUp_isLoaded = false;
  event_wgt_adjustment_QCDScaleDn_isLoaded = false;
  event_wgt_adjustment_QCDScaleUp_isLoaded = false;
  event_wgt_triggers_isLoaded = false;
  event_wgt_triggers_Dilepton_isLoaded = false;
  event_wgt_triggers_Dilepton_DF_Extra_isLoaded = false;
  event_wgt_triggers_PFHT_Control_isLoaded = false;
  event_wgt_triggers_PFMET_MHT_Control_isLoaded = false;
  event_wgt_triggers_SingleLepton_isLoaded = false;
  event_wgt_triggers_SinglePhoton_isLoaded = false;
  genak4jets_eta_isLoaded = false;
  genak4jets_mass_isLoaded = false;
  genak4jets_phi_isLoaded = false;
  genak4jets_pt_isLoaded = false;
  genmet_pTmiss_isLoaded = false;
  genmet_phimiss_isLoaded = false;
  lepton_eff_isLoaded = false;
  lepton_eff_StatDn_isLoaded = false;
  lepton_eff_StatUp_isLoaded = false;
  lepton_eff_SystDn_isLoaded = false;
  lepton_eff_SystUp_isLoaded = false;
  lepton_eta_isLoaded = false;
  lepton_id_isLoaded = false;
  lepton_is_genMatched_prompt_isLoaded = false;
  lepton_mass_isLoaded = false;
  lepton_phi_isLoaded = false;
  lepton_pt_isLoaded = false;
  genpromptparticles_sump4_mass_isLoaded = false;
  genpromptparticles_sump4_pt_isLoaded = false;
  genpromptparticles_sump4_rapidity_isLoaded = false;
  leptons_eff_isLoaded = false;
  leptons_eff_DF_isLoaded = false;
  leptons_eff_DF_StatDn_isLoaded = false;
  leptons_eff_DF_StatUp_isLoaded = false;
  leptons_eff_DF_SystDn_isLoaded = false;
  leptons_eff_DF_SystUp_isLoaded = false;
  leptons_eff_StatDn_isLoaded = false;
  leptons_eff_StatUp_isLoaded = false;
  leptons_eff_SystDn_isLoaded = false;
  leptons_eff_SystUp_isLoaded = false;
  leptons_eta_isLoaded = false;
  leptons_id_isLoaded = false;
  leptons_is_TOmatched_SingleLepton_isLoaded = false;
  leptons_is_genMatched_prompt_isLoaded = false;
  leptons_mass_isLoaded = false;
  leptons_phi_isLoaded = false;
  leptons_pt_isLoaded = false;
  min_abs_dPhi_pTj_pTmiss_isLoaded = false;
  pConst_JJQCD_SIG_ghg2_1_JHUGen_isLoaded = false;
  pConst_JJVBF_SIG_ghv1_1_JHUGen_isLoaded = false;
  p_JJQCD_SIG_ghg2_1_JHUGen_isLoaded = false;
  p_JJVBF_SIG_ghv1_1_JHUGen_isLoaded = false;
  p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_isLoaded = false;
  p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_isLoaded = false;
  p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_isLoaded = false;
  p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_isLoaded = false;
  p_JJVBF_SIG_ghv2_1_JHUGen_isLoaded = false;
  p_JJVBF_SIG_ghv4_1_JHUGen_isLoaded = false;
  p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_isLoaded = false;
  photon_MIPTotalEnergy_isLoaded = false;
  photon_eta_isLoaded = false;
  photon_full5x5_r9_isLoaded = false;
  photon_full5x5_sigmaIEtaIEta_isLoaded = false;
  photon_full5x5_sigmaIPhiIPhi_isLoaded = false;
  photon_isEB_isLoaded = false;
  photon_isEBEEGap_isLoaded = false;
  photon_isEE_isLoaded = false;
  photon_isGap_isLoaded = false;
  photon_is_METSafe_isLoaded = false;
  photon_is_PFID_isLoaded = false;
  photon_is_beamHaloSafe_isLoaded = false;
  photon_is_conversionSafe_isLoaded = false;
  photon_is_genMatched_prompt_isLoaded = false;
  photon_is_inTime_isLoaded = false;
  photon_is_spikeSafe_isLoaded = false;
  photon_mass_isLoaded = false;
  photon_pass_HGGSelection_isLoaded = false;
  photon_phi_isLoaded = false;
  photon_pt_isLoaded = false;
  photon_seedTime_isLoaded = false;
}
void SkimTree::LoadAllBranches() {
  if (EventNumber_branch != 0) EventNumber();
  if (KFactor_EW_NLO_qqVV_Bkg_EWDn_branch != 0) KFactor_EW_NLO_qqVV_Bkg_EWDn();
  if (KFactor_EW_NLO_qqVV_Bkg_EWUp_branch != 0) KFactor_EW_NLO_qqVV_Bkg_EWUp();
  if (KFactor_EW_NLO_qqVV_Bkg_Nominal_branch != 0) KFactor_EW_NLO_qqVV_Bkg_Nominal();
  if (KFactor_EW_NLO_qqVV_Bkg_arg_mass_branch != 0) KFactor_EW_NLO_qqVV_Bkg_arg_mass();
  if (KFactor_EW_NLO_qqVV_Bkg_arg_pthat_branch != 0) KFactor_EW_NLO_qqVV_Bkg_arg_pthat();
  if (KFactor_EW_NLO_qqVV_Bkg_arg_rho_branch != 0) KFactor_EW_NLO_qqVV_Bkg_arg_rho();
  if (KFactor_EW_NLO_qqVV_Bkg_arg_that_branch != 0) KFactor_EW_NLO_qqVV_Bkg_arg_that();
  if (KFactor_QCD_NNLO_qqVV_Bkg_Nominal_branch != 0) KFactor_QCD_NNLO_qqVV_Bkg_Nominal();
  if (KFactor_QCD_NNLO_qqVV_Bkg_arg_mass_branch != 0) KFactor_QCD_NNLO_qqVV_Bkg_arg_mass();
  if (KFactor_QCD_NNLO_qqVV_Bkg_arg_pthat_branch != 0) KFactor_QCD_NNLO_qqVV_Bkg_arg_pthat();
  if (LuminosityBlock_branch != 0) LuminosityBlock();
  if (RunNumber_branch != 0) RunNumber();
  if (ak4jets_CEMF_branch != 0) ak4jets_CEMF();
  if (ak4jets_HT_branch != 0) ak4jets_HT();
  if (ak4jets_MHT_branch != 0) ak4jets_MHT();
  if (ak4jets_NEMF_branch != 0) ak4jets_NEMF();
  if (ak4jets_btagWP_Bits_branch != 0) ak4jets_btagWP_Bits();
  if (ak4jets_eta_branch != 0) ak4jets_eta();
  if (ak4jets_is_genMatched_branch != 0) ak4jets_is_genMatched();
  if (ak4jets_is_genMatched_fullCone_branch != 0) ak4jets_is_genMatched_fullCone();
  if (ak4jets_mass_branch != 0) ak4jets_mass();
  if (ak4jets_phi_branch != 0) ak4jets_phi();
  if (ak4jets_pt_branch != 0) ak4jets_pt();
  if (ak8jets_eta_branch != 0) ak8jets_eta();
  if (ak8jets_mass_branch != 0) ak8jets_mass();
  if (ak8jets_phi_branch != 0) ak8jets_phi();
  if (ak8jets_pt_branch != 0) ak8jets_pt();
  if (dPhi_pTboson_pTmiss_branch != 0) dPhi_pTboson_pTmiss();
  if (dPhi_pTbosonjets_pTmiss_branch != 0) dPhi_pTbosonjets_pTmiss();
  if (dilepton_eta_branch != 0) dilepton_eta();
  if (dilepton_id_branch != 0) dilepton_id();
  if (dilepton_mass_branch != 0) dilepton_mass();
  if (dilepton_phi_branch != 0) dilepton_phi();
  if (dilepton_pt_branch != 0) dilepton_pt();
  if (electron_full5x5_r9_branch != 0) electron_full5x5_r9();
  if (electron_full5x5_sigmaIEtaIEta_branch != 0) electron_full5x5_sigmaIEtaIEta();
  if (electron_full5x5_sigmaIPhiIPhi_branch != 0) electron_full5x5_sigmaIPhiIPhi();
  if (electron_seedTime_branch != 0) electron_seedTime();
  if (electrons_full5x5_r9_branch != 0) electrons_full5x5_r9();
  if (electrons_full5x5_sigmaIEtaIEta_branch != 0) electrons_full5x5_sigmaIEtaIEta();
  if (electrons_full5x5_sigmaIPhiIPhi_branch != 0) electrons_full5x5_sigmaIPhiIPhi();
  if (electrons_seedTime_branch != 0) electrons_seedTime();
  if (event_mTZZ_branch != 0) event_mTZZ();
  if (event_mZZ_branch != 0) event_mZZ();
  if (event_mlg_branch != 0) event_mlg();
  if (event_mllg_branch != 0) event_mllg();
  if (event_n_ak4jets_pt20_branch != 0) event_n_ak4jets_pt20();
  if (event_n_ak4jets_pt20_btagged_loose_branch != 0) event_n_ak4jets_pt20_btagged_loose();
  if (event_n_ak4jets_pt20_btagged_medium_branch != 0) event_n_ak4jets_pt20_btagged_medium();
  if (event_n_ak4jets_pt30_branch != 0) event_n_ak4jets_pt30();
  if (event_n_ak4jets_pt30_btagged_loose_branch != 0) event_n_ak4jets_pt30_btagged_loose();
  if (event_n_ak4jets_pt30_btagged_medium_branch != 0) event_n_ak4jets_pt30_btagged_medium();
  if (event_n_leptons_fakeableBase_branch != 0) event_n_leptons_fakeableBase();
  if (event_n_vtxs_good_branch != 0) event_n_vtxs_good();
  if (event_pTmiss_branch != 0) event_pTmiss();
  if (event_pass_tightMETFilters_branch != 0) event_pass_tightMETFilters();
  if (event_phimiss_branch != 0) event_phimiss();
  if (event_wgt_branch != 0) event_wgt();
  if (event_wgt_L1PrefiringDn_branch != 0) event_wgt_L1PrefiringDn();
  if (event_wgt_L1PrefiringUp_branch != 0) event_wgt_L1PrefiringUp();
  if (event_wgt_PUDn_branch != 0) event_wgt_PUDn();
  if (event_wgt_PUUp_branch != 0) event_wgt_PUUp();
  if (event_wgt_SFs_PUJetId_branch != 0) event_wgt_SFs_PUJetId();
  if (event_wgt_SFs_PUJetId_EffDn_branch != 0) event_wgt_SFs_PUJetId_EffDn();
  if (event_wgt_SFs_PUJetId_EffUp_branch != 0) event_wgt_SFs_PUJetId_EffUp();
  if (event_wgt_SFs_btagging_branch != 0) event_wgt_SFs_btagging();
  if (event_wgt_SFs_btagging_EffDn_branch != 0) event_wgt_SFs_btagging_EffDn();
  if (event_wgt_SFs_btagging_EffUp_branch != 0) event_wgt_SFs_btagging_EffUp();
  if (event_wgt_SFs_electrons_branch != 0) event_wgt_SFs_electrons();
  if (event_wgt_SFs_electrons_AltMCDn_branch != 0) event_wgt_SFs_electrons_AltMCDn();
  if (event_wgt_SFs_electrons_AltMCUp_branch != 0) event_wgt_SFs_electrons_AltMCUp();
  if (event_wgt_SFs_electrons_StatDn_branch != 0) event_wgt_SFs_electrons_StatDn();
  if (event_wgt_SFs_electrons_StatUp_branch != 0) event_wgt_SFs_electrons_StatUp();
  if (event_wgt_SFs_electrons_SystDn_branch != 0) event_wgt_SFs_electrons_SystDn();
  if (event_wgt_SFs_electrons_SystUp_branch != 0) event_wgt_SFs_electrons_SystUp();
  if (event_wgt_SFs_muons_branch != 0) event_wgt_SFs_muons();
  if (event_wgt_SFs_muons_AltMCDn_branch != 0) event_wgt_SFs_muons_AltMCDn();
  if (event_wgt_SFs_muons_AltMCUp_branch != 0) event_wgt_SFs_muons_AltMCUp();
  if (event_wgt_SFs_muons_StatDn_branch != 0) event_wgt_SFs_muons_StatDn();
  if (event_wgt_SFs_muons_StatUp_branch != 0) event_wgt_SFs_muons_StatUp();
  if (event_wgt_SFs_muons_SystDn_branch != 0) event_wgt_SFs_muons_SystDn();
  if (event_wgt_SFs_muons_SystUp_branch != 0) event_wgt_SFs_muons_SystUp();
  if (event_wgt_SFs_photons_branch != 0) event_wgt_SFs_photons();
  if (event_wgt_SFs_photons_EffDn_branch != 0) event_wgt_SFs_photons_EffDn();
  if (event_wgt_SFs_photons_EffUp_branch != 0) event_wgt_SFs_photons_EffUp();
  if (event_wgt_adjustment_AsMZDn_branch != 0) event_wgt_adjustment_AsMZDn();
  if (event_wgt_adjustment_AsMZUp_branch != 0) event_wgt_adjustment_AsMZUp();
  if (event_wgt_adjustment_NNPDF30_branch != 0) event_wgt_adjustment_NNPDF30();
  if (event_wgt_adjustment_NNPDF30_AsMZDn_branch != 0) event_wgt_adjustment_NNPDF30_AsMZDn();
  if (event_wgt_adjustment_NNPDF30_AsMZUp_branch != 0) event_wgt_adjustment_NNPDF30_AsMZUp();
  if (event_wgt_adjustment_NNPDF30_PDFReplicaDn_branch != 0) event_wgt_adjustment_NNPDF30_PDFReplicaDn();
  if (event_wgt_adjustment_NNPDF30_PDFReplicaUp_branch != 0) event_wgt_adjustment_NNPDF30_PDFReplicaUp();
  if (event_wgt_adjustment_PDFReplicaDn_branch != 0) event_wgt_adjustment_PDFReplicaDn();
  if (event_wgt_adjustment_PDFReplicaUp_branch != 0) event_wgt_adjustment_PDFReplicaUp();
  if (event_wgt_adjustment_PDFScaleDn_branch != 0) event_wgt_adjustment_PDFScaleDn();
  if (event_wgt_adjustment_PDFScaleUp_branch != 0) event_wgt_adjustment_PDFScaleUp();
  if (event_wgt_adjustment_PythiaScaleDn_branch != 0) event_wgt_adjustment_PythiaScaleDn();
  if (event_wgt_adjustment_PythiaScaleUp_branch != 0) event_wgt_adjustment_PythiaScaleUp();
  if (event_wgt_adjustment_QCDScaleDn_branch != 0) event_wgt_adjustment_QCDScaleDn();
  if (event_wgt_adjustment_QCDScaleUp_branch != 0) event_wgt_adjustment_QCDScaleUp();
  if (event_wgt_triggers_branch != 0) event_wgt_triggers();
  if (event_wgt_triggers_Dilepton_branch != 0) event_wgt_triggers_Dilepton();
  if (event_wgt_triggers_Dilepton_DF_Extra_branch != 0) event_wgt_triggers_Dilepton_DF_Extra();
  if (event_wgt_triggers_PFHT_Control_branch != 0) event_wgt_triggers_PFHT_Control();
  if (event_wgt_triggers_PFMET_MHT_Control_branch != 0) event_wgt_triggers_PFMET_MHT_Control();
  if (event_wgt_triggers_SingleLepton_branch != 0) event_wgt_triggers_SingleLepton();
  if (event_wgt_triggers_SinglePhoton_branch != 0) event_wgt_triggers_SinglePhoton();
  if (genak4jets_eta_branch != 0) genak4jets_eta();
  if (genak4jets_mass_branch != 0) genak4jets_mass();
  if (genak4jets_phi_branch != 0) genak4jets_phi();
  if (genak4jets_pt_branch != 0) genak4jets_pt();
  if (genmet_pTmiss_branch != 0) genmet_pTmiss();
  if (genmet_phimiss_branch != 0) genmet_phimiss();
  if (lepton_eff_branch != 0) lepton_eff();
  if (lepton_eff_StatDn_branch != 0) lepton_eff_StatDn();
  if (lepton_eff_StatUp_branch != 0) lepton_eff_StatUp();
  if (lepton_eff_SystDn_branch != 0) lepton_eff_SystDn();
  if (lepton_eff_SystUp_branch != 0) lepton_eff_SystUp();
  if (lepton_eta_branch != 0) lepton_eta();
  if (lepton_id_branch != 0) lepton_id();
  if (lepton_is_genMatched_prompt_branch != 0) lepton_is_genMatched_prompt();
  if (lepton_mass_branch != 0) lepton_mass();
  if (lepton_phi_branch != 0) lepton_phi();
  if (lepton_pt_branch != 0) lepton_pt();
  if (genpromptparticles_sump4_mass_branch != 0) genpromptparticles_sump4_mass();
  if (genpromptparticles_sump4_pt_branch != 0) genpromptparticles_sump4_pt();
  if (genpromptparticles_sump4_rapidity_branch != 0) genpromptparticles_sump4_rapidity();
  if (leptons_eff_branch != 0) leptons_eff();
  if (leptons_eff_DF_branch != 0) leptons_eff_DF();
  if (leptons_eff_DF_StatDn_branch != 0) leptons_eff_DF_StatDn();
  if (leptons_eff_DF_StatUp_branch != 0) leptons_eff_DF_StatUp();
  if (leptons_eff_DF_SystDn_branch != 0) leptons_eff_DF_SystDn();
  if (leptons_eff_DF_SystUp_branch != 0) leptons_eff_DF_SystUp();
  if (leptons_eff_StatDn_branch != 0) leptons_eff_StatDn();
  if (leptons_eff_StatUp_branch != 0) leptons_eff_StatUp();
  if (leptons_eff_SystDn_branch != 0) leptons_eff_SystDn();
  if (leptons_eff_SystUp_branch != 0) leptons_eff_SystUp();
  if (leptons_eta_branch != 0) leptons_eta();
  if (leptons_id_branch != 0) leptons_id();
  if (leptons_is_TOmatched_SingleLepton_branch != 0) leptons_is_TOmatched_SingleLepton();
  if (leptons_is_genMatched_prompt_branch != 0) leptons_is_genMatched_prompt();
  if (leptons_mass_branch != 0) leptons_mass();
  if (leptons_phi_branch != 0) leptons_phi();
  if (leptons_pt_branch != 0) leptons_pt();
  if (min_abs_dPhi_pTj_pTmiss_branch != 0) min_abs_dPhi_pTj_pTmiss();
  if (pConst_JJQCD_SIG_ghg2_1_JHUGen_branch != 0) pConst_JJQCD_SIG_ghg2_1_JHUGen();
  if (pConst_JJVBF_SIG_ghv1_1_JHUGen_branch != 0) pConst_JJVBF_SIG_ghv1_1_JHUGen();
  if (p_JJQCD_SIG_ghg2_1_JHUGen_branch != 0) p_JJQCD_SIG_ghg2_1_JHUGen();
  if (p_JJVBF_SIG_ghv1_1_JHUGen_branch != 0) p_JJVBF_SIG_ghv1_1_JHUGen();
  if (p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_branch != 0) p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen();
  if (p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_branch != 0) p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen();
  if (p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_branch != 0) p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen();
  if (p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_branch != 0) p_JJVBF_SIG_ghv1prime2_1E4_JHUGen();
  if (p_JJVBF_SIG_ghv2_1_JHUGen_branch != 0) p_JJVBF_SIG_ghv2_1_JHUGen();
  if (p_JJVBF_SIG_ghv4_1_JHUGen_branch != 0) p_JJVBF_SIG_ghv4_1_JHUGen();
  if (p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_branch != 0) p_JJVBF_SIG_ghza1prime2_1E4_JHUGen();
  if (photon_MIPTotalEnergy_branch != 0) photon_MIPTotalEnergy();
  if (photon_eta_branch != 0) photon_eta();
  if (photon_full5x5_r9_branch != 0) photon_full5x5_r9();
  if (photon_full5x5_sigmaIEtaIEta_branch != 0) photon_full5x5_sigmaIEtaIEta();
  if (photon_full5x5_sigmaIPhiIPhi_branch != 0) photon_full5x5_sigmaIPhiIPhi();
  if (photon_isEB_branch != 0) photon_isEB();
  if (photon_isEBEEGap_branch != 0) photon_isEBEEGap();
  if (photon_isEE_branch != 0) photon_isEE();
  if (photon_isGap_branch != 0) photon_isGap();
  if (photon_is_METSafe_branch != 0) photon_is_METSafe();
  if (photon_is_PFID_branch != 0) photon_is_PFID();
  if (photon_is_beamHaloSafe_branch != 0) photon_is_beamHaloSafe();
  if (photon_is_conversionSafe_branch != 0) photon_is_conversionSafe();
  if (photon_is_genMatched_prompt_branch != 0) photon_is_genMatched_prompt();
  if (photon_is_inTime_branch != 0) photon_is_inTime();
  if (photon_is_spikeSafe_branch != 0) photon_is_spikeSafe();
  if (photon_mass_branch != 0) photon_mass();
  if (photon_pass_HGGSelection_branch != 0) photon_pass_HGGSelection();
  if (photon_phi_branch != 0) photon_phi();
  if (photon_pt_branch != 0) photon_pt();
  if (photon_seedTime_branch != 0) photon_seedTime();
}
const unsigned long long &SkimTree::EventNumber() {
  if (not EventNumber_isLoaded) {
    if (EventNumber_branch != 0) {
      EventNumber_branch->GetEntry(index);
    } else {
      printf("branch EventNumber_branch does not exist!\n");
      exit(1);
    }
    EventNumber_isLoaded = true;
  }
  return EventNumber_;
}
const float &SkimTree::KFactor_EW_NLO_qqVV_Bkg_EWDn() {
  if (not KFactor_EW_NLO_qqVV_Bkg_EWDn_isLoaded) {
    if (KFactor_EW_NLO_qqVV_Bkg_EWDn_branch != 0) {
      KFactor_EW_NLO_qqVV_Bkg_EWDn_branch->GetEntry(index);
    } else {
      printf("branch KFactor_EW_NLO_qqVV_Bkg_EWDn_branch does not exist!\n");
      exit(1);
    }
    KFactor_EW_NLO_qqVV_Bkg_EWDn_isLoaded = true;
  }
  return KFactor_EW_NLO_qqVV_Bkg_EWDn_;
}
const float &SkimTree::KFactor_EW_NLO_qqVV_Bkg_EWUp() {
  if (not KFactor_EW_NLO_qqVV_Bkg_EWUp_isLoaded) {
    if (KFactor_EW_NLO_qqVV_Bkg_EWUp_branch != 0) {
      KFactor_EW_NLO_qqVV_Bkg_EWUp_branch->GetEntry(index);
    } else {
      printf("branch KFactor_EW_NLO_qqVV_Bkg_EWUp_branch does not exist!\n");
      exit(1);
    }
    KFactor_EW_NLO_qqVV_Bkg_EWUp_isLoaded = true;
  }
  return KFactor_EW_NLO_qqVV_Bkg_EWUp_;
}
const float &SkimTree::KFactor_EW_NLO_qqVV_Bkg_Nominal() {
  if (not KFactor_EW_NLO_qqVV_Bkg_Nominal_isLoaded) {
    if (KFactor_EW_NLO_qqVV_Bkg_Nominal_branch != 0) {
      KFactor_EW_NLO_qqVV_Bkg_Nominal_branch->GetEntry(index);
    } else {
      printf("branch KFactor_EW_NLO_qqVV_Bkg_Nominal_branch does not exist!\n");
      exit(1);
    }
    KFactor_EW_NLO_qqVV_Bkg_Nominal_isLoaded = true;
  }
  return KFactor_EW_NLO_qqVV_Bkg_Nominal_;
}
const float &SkimTree::KFactor_EW_NLO_qqVV_Bkg_arg_mass() {
  if (not KFactor_EW_NLO_qqVV_Bkg_arg_mass_isLoaded) {
    if (KFactor_EW_NLO_qqVV_Bkg_arg_mass_branch != 0) {
      KFactor_EW_NLO_qqVV_Bkg_arg_mass_branch->GetEntry(index);
    } else {
      printf("branch KFactor_EW_NLO_qqVV_Bkg_arg_mass_branch does not exist!\n");
      exit(1);
    }
    KFactor_EW_NLO_qqVV_Bkg_arg_mass_isLoaded = true;
  }
  return KFactor_EW_NLO_qqVV_Bkg_arg_mass_;
}
const float &SkimTree::KFactor_EW_NLO_qqVV_Bkg_arg_pthat() {
  if (not KFactor_EW_NLO_qqVV_Bkg_arg_pthat_isLoaded) {
    if (KFactor_EW_NLO_qqVV_Bkg_arg_pthat_branch != 0) {
      KFactor_EW_NLO_qqVV_Bkg_arg_pthat_branch->GetEntry(index);
    } else {
      printf("branch KFactor_EW_NLO_qqVV_Bkg_arg_pthat_branch does not exist!\n");
      exit(1);
    }
    KFactor_EW_NLO_qqVV_Bkg_arg_pthat_isLoaded = true;
  }
  return KFactor_EW_NLO_qqVV_Bkg_arg_pthat_;
}
const float &SkimTree::KFactor_EW_NLO_qqVV_Bkg_arg_rho() {
  if (not KFactor_EW_NLO_qqVV_Bkg_arg_rho_isLoaded) {
    if (KFactor_EW_NLO_qqVV_Bkg_arg_rho_branch != 0) {
      KFactor_EW_NLO_qqVV_Bkg_arg_rho_branch->GetEntry(index);
    } else {
      printf("branch KFactor_EW_NLO_qqVV_Bkg_arg_rho_branch does not exist!\n");
      exit(1);
    }
    KFactor_EW_NLO_qqVV_Bkg_arg_rho_isLoaded = true;
  }
  return KFactor_EW_NLO_qqVV_Bkg_arg_rho_;
}
const float &SkimTree::KFactor_EW_NLO_qqVV_Bkg_arg_that() {
  if (not KFactor_EW_NLO_qqVV_Bkg_arg_that_isLoaded) {
    if (KFactor_EW_NLO_qqVV_Bkg_arg_that_branch != 0) {
      KFactor_EW_NLO_qqVV_Bkg_arg_that_branch->GetEntry(index);
    } else {
      printf("branch KFactor_EW_NLO_qqVV_Bkg_arg_that_branch does not exist!\n");
      exit(1);
    }
    KFactor_EW_NLO_qqVV_Bkg_arg_that_isLoaded = true;
  }
  return KFactor_EW_NLO_qqVV_Bkg_arg_that_;
}
const float &SkimTree::KFactor_QCD_NNLO_qqVV_Bkg_Nominal() {
  if (not KFactor_QCD_NNLO_qqVV_Bkg_Nominal_isLoaded) {
    if (KFactor_QCD_NNLO_qqVV_Bkg_Nominal_branch != 0) {
      KFactor_QCD_NNLO_qqVV_Bkg_Nominal_branch->GetEntry(index);
    } else {
      printf("branch KFactor_QCD_NNLO_qqVV_Bkg_Nominal_branch does not exist!\n");
      exit(1);
    }
    KFactor_QCD_NNLO_qqVV_Bkg_Nominal_isLoaded = true;
  }
  return KFactor_QCD_NNLO_qqVV_Bkg_Nominal_;
}
const float &SkimTree::KFactor_QCD_NNLO_qqVV_Bkg_arg_mass() {
  if (not KFactor_QCD_NNLO_qqVV_Bkg_arg_mass_isLoaded) {
    if (KFactor_QCD_NNLO_qqVV_Bkg_arg_mass_branch != 0) {
      KFactor_QCD_NNLO_qqVV_Bkg_arg_mass_branch->GetEntry(index);
    } else {
      printf("branch KFactor_QCD_NNLO_qqVV_Bkg_arg_mass_branch does not exist!\n");
      exit(1);
    }
    KFactor_QCD_NNLO_qqVV_Bkg_arg_mass_isLoaded = true;
  }
  return KFactor_QCD_NNLO_qqVV_Bkg_arg_mass_;
}
const float &SkimTree::KFactor_QCD_NNLO_qqVV_Bkg_arg_pthat() {
  if (not KFactor_QCD_NNLO_qqVV_Bkg_arg_pthat_isLoaded) {
    if (KFactor_QCD_NNLO_qqVV_Bkg_arg_pthat_branch != 0) {
      KFactor_QCD_NNLO_qqVV_Bkg_arg_pthat_branch->GetEntry(index);
    } else {
      printf("branch KFactor_QCD_NNLO_qqVV_Bkg_arg_pthat_branch does not exist!\n");
      exit(1);
    }
    KFactor_QCD_NNLO_qqVV_Bkg_arg_pthat_isLoaded = true;
  }
  return KFactor_QCD_NNLO_qqVV_Bkg_arg_pthat_;
}
const unsigned int &SkimTree::LuminosityBlock() {
  if (not LuminosityBlock_isLoaded) {
    if (LuminosityBlock_branch != 0) {
      LuminosityBlock_branch->GetEntry(index);
    } else {
      printf("branch LuminosityBlock_branch does not exist!\n");
      exit(1);
    }
    LuminosityBlock_isLoaded = true;
  }
  return LuminosityBlock_;
}
const unsigned int &SkimTree::RunNumber() {
  if (not RunNumber_isLoaded) {
    if (RunNumber_branch != 0) {
      RunNumber_branch->GetEntry(index);
    } else {
      printf("branch RunNumber_branch does not exist!\n");
      exit(1);
    }
    RunNumber_isLoaded = true;
  }
  return RunNumber_;
}
const vector<float> &SkimTree::ak4jets_CEMF() {
  if (not ak4jets_CEMF_isLoaded) {
    if (ak4jets_CEMF_branch != 0) {
      ak4jets_CEMF_branch->GetEntry(index);
    } else {
      printf("branch ak4jets_CEMF_branch does not exist!\n");
      exit(1);
    }
    ak4jets_CEMF_isLoaded = true;
  }
  return *ak4jets_CEMF_;
}
const float &SkimTree::ak4jets_HT() {
  if (not ak4jets_HT_isLoaded) {
    if (ak4jets_HT_branch != 0) {
      ak4jets_HT_branch->GetEntry(index);
    } else {
      printf("branch ak4jets_HT_branch does not exist!\n");
      exit(1);
    }
    ak4jets_HT_isLoaded = true;
  }
  return ak4jets_HT_;
}
const float &SkimTree::ak4jets_MHT() {
  if (not ak4jets_MHT_isLoaded) {
    if (ak4jets_MHT_branch != 0) {
      ak4jets_MHT_branch->GetEntry(index);
    } else {
      printf("branch ak4jets_MHT_branch does not exist!\n");
      exit(1);
    }
    ak4jets_MHT_isLoaded = true;
  }
  return ak4jets_MHT_;
}
const vector<float> &SkimTree::ak4jets_NEMF() {
  if (not ak4jets_NEMF_isLoaded) {
    if (ak4jets_NEMF_branch != 0) {
      ak4jets_NEMF_branch->GetEntry(index);
    } else {
      printf("branch ak4jets_NEMF_branch does not exist!\n");
      exit(1);
    }
    ak4jets_NEMF_isLoaded = true;
  }
  return *ak4jets_NEMF_;
}
const vector<unsigned char> &SkimTree::ak4jets_btagWP_Bits() {
  if (not ak4jets_btagWP_Bits_isLoaded) {
    if (ak4jets_btagWP_Bits_branch != 0) {
      ak4jets_btagWP_Bits_branch->GetEntry(index);
    } else {
      printf("branch ak4jets_btagWP_Bits_branch does not exist!\n");
      exit(1);
    }
    ak4jets_btagWP_Bits_isLoaded = true;
  }
  return *ak4jets_btagWP_Bits_;
}
const vector<float> &SkimTree::ak4jets_eta() {
  if (not ak4jets_eta_isLoaded) {
    if (ak4jets_eta_branch != 0) {
      ak4jets_eta_branch->GetEntry(index);
    } else {
      printf("branch ak4jets_eta_branch does not exist!\n");
      exit(1);
    }
    ak4jets_eta_isLoaded = true;
  }
  return *ak4jets_eta_;
}
const vector<bool> &SkimTree::ak4jets_is_genMatched() {
  if (not ak4jets_is_genMatched_isLoaded) {
    if (ak4jets_is_genMatched_branch != 0) {
      ak4jets_is_genMatched_branch->GetEntry(index);
    } else {
      printf("branch ak4jets_is_genMatched_branch does not exist!\n");
      exit(1);
    }
    ak4jets_is_genMatched_isLoaded = true;
  }
  return *ak4jets_is_genMatched_;
}
const vector<bool> &SkimTree::ak4jets_is_genMatched_fullCone() {
  if (not ak4jets_is_genMatched_fullCone_isLoaded) {
    if (ak4jets_is_genMatched_fullCone_branch != 0) {
      ak4jets_is_genMatched_fullCone_branch->GetEntry(index);
    } else {
      printf("branch ak4jets_is_genMatched_fullCone_branch does not exist!\n");
      exit(1);
    }
    ak4jets_is_genMatched_fullCone_isLoaded = true;
  }
  return *ak4jets_is_genMatched_fullCone_;
}
const vector<float> &SkimTree::ak4jets_mass() {
  if (not ak4jets_mass_isLoaded) {
    if (ak4jets_mass_branch != 0) {
      ak4jets_mass_branch->GetEntry(index);
    } else {
      printf("branch ak4jets_mass_branch does not exist!\n");
      exit(1);
    }
    ak4jets_mass_isLoaded = true;
  }
  return *ak4jets_mass_;
}
const vector<float> &SkimTree::ak4jets_phi() {
  if (not ak4jets_phi_isLoaded) {
    if (ak4jets_phi_branch != 0) {
      ak4jets_phi_branch->GetEntry(index);
    } else {
      printf("branch ak4jets_phi_branch does not exist!\n");
      exit(1);
    }
    ak4jets_phi_isLoaded = true;
  }
  return *ak4jets_phi_;
}
const vector<float> &SkimTree::ak4jets_pt() {
  if (not ak4jets_pt_isLoaded) {
    if (ak4jets_pt_branch != 0) {
      ak4jets_pt_branch->GetEntry(index);
    } else {
      printf("branch ak4jets_pt_branch does not exist!\n");
      exit(1);
    }
    ak4jets_pt_isLoaded = true;
  }
  return *ak4jets_pt_;
}
const vector<float> &SkimTree::ak8jets_eta() {
  if (not ak8jets_eta_isLoaded) {
    if (ak8jets_eta_branch != 0) {
      ak8jets_eta_branch->GetEntry(index);
    } else {
      printf("branch ak8jets_eta_branch does not exist!\n");
      exit(1);
    }
    ak8jets_eta_isLoaded = true;
  }
  return *ak8jets_eta_;
}
const vector<float> &SkimTree::ak8jets_mass() {
  if (not ak8jets_mass_isLoaded) {
    if (ak8jets_mass_branch != 0) {
      ak8jets_mass_branch->GetEntry(index);
    } else {
      printf("branch ak8jets_mass_branch does not exist!\n");
      exit(1);
    }
    ak8jets_mass_isLoaded = true;
  }
  return *ak8jets_mass_;
}
const vector<float> &SkimTree::ak8jets_phi() {
  if (not ak8jets_phi_isLoaded) {
    if (ak8jets_phi_branch != 0) {
      ak8jets_phi_branch->GetEntry(index);
    } else {
      printf("branch ak8jets_phi_branch does not exist!\n");
      exit(1);
    }
    ak8jets_phi_isLoaded = true;
  }
  return *ak8jets_phi_;
}
const vector<float> &SkimTree::ak8jets_pt() {
  if (not ak8jets_pt_isLoaded) {
    if (ak8jets_pt_branch != 0) {
      ak8jets_pt_branch->GetEntry(index);
    } else {
      printf("branch ak8jets_pt_branch does not exist!\n");
      exit(1);
    }
    ak8jets_pt_isLoaded = true;
  }
  return *ak8jets_pt_;
}
const float &SkimTree::dPhi_pTboson_pTmiss() {
  if (not dPhi_pTboson_pTmiss_isLoaded) {
    if (dPhi_pTboson_pTmiss_branch != 0) {
      dPhi_pTboson_pTmiss_branch->GetEntry(index);
    } else {
      printf("branch dPhi_pTboson_pTmiss_branch does not exist!\n");
      exit(1);
    }
    dPhi_pTboson_pTmiss_isLoaded = true;
  }
  return dPhi_pTboson_pTmiss_;
}
const float &SkimTree::dPhi_pTbosonjets_pTmiss() {
  if (not dPhi_pTbosonjets_pTmiss_isLoaded) {
    if (dPhi_pTbosonjets_pTmiss_branch != 0) {
      dPhi_pTbosonjets_pTmiss_branch->GetEntry(index);
    } else {
      printf("branch dPhi_pTbosonjets_pTmiss_branch does not exist!\n");
      exit(1);
    }
    dPhi_pTbosonjets_pTmiss_isLoaded = true;
  }
  return dPhi_pTbosonjets_pTmiss_;
}
const float &SkimTree::dilepton_eta() {
  if (not dilepton_eta_isLoaded) {
    if (dilepton_eta_branch != 0) {
      dilepton_eta_branch->GetEntry(index);
    } else {
      printf("branch dilepton_eta_branch does not exist!\n");
      exit(1);
    }
    dilepton_eta_isLoaded = true;
  }
  return dilepton_eta_;
}
const int &SkimTree::dilepton_id() {
  if (not dilepton_id_isLoaded) {
    if (dilepton_id_branch != 0) {
      dilepton_id_branch->GetEntry(index);
    } else {
      printf("branch dilepton_id_branch does not exist!\n");
      exit(1);
    }
    dilepton_id_isLoaded = true;
  }
  return dilepton_id_;
}
const float &SkimTree::dilepton_mass() {
  if (not dilepton_mass_isLoaded) {
    if (dilepton_mass_branch != 0) {
      dilepton_mass_branch->GetEntry(index);
    } else {
      printf("branch dilepton_mass_branch does not exist!\n");
      exit(1);
    }
    dilepton_mass_isLoaded = true;
  }
  return dilepton_mass_;
}
const float &SkimTree::dilepton_phi() {
  if (not dilepton_phi_isLoaded) {
    if (dilepton_phi_branch != 0) {
      dilepton_phi_branch->GetEntry(index);
    } else {
      printf("branch dilepton_phi_branch does not exist!\n");
      exit(1);
    }
    dilepton_phi_isLoaded = true;
  }
  return dilepton_phi_;
}
const float &SkimTree::dilepton_pt() {
  if (not dilepton_pt_isLoaded) {
    if (dilepton_pt_branch != 0) {
      dilepton_pt_branch->GetEntry(index);
    } else {
      printf("branch dilepton_pt_branch does not exist!\n");
      exit(1);
    }
    dilepton_pt_isLoaded = true;
  }
  return dilepton_pt_;
}
const float &SkimTree::electron_full5x5_r9() {
  if (not electron_full5x5_r9_isLoaded) {
    if (electron_full5x5_r9_branch != 0) {
      electron_full5x5_r9_branch->GetEntry(index);
    } else {
      printf("branch electron_full5x5_r9_branch does not exist!\n");
      exit(1);
    }
    electron_full5x5_r9_isLoaded = true;
  }
  return electron_full5x5_r9_;
}
const float &SkimTree::electron_full5x5_sigmaIEtaIEta() {
  if (not electron_full5x5_sigmaIEtaIEta_isLoaded) {
    if (electron_full5x5_sigmaIEtaIEta_branch != 0) {
      electron_full5x5_sigmaIEtaIEta_branch->GetEntry(index);
    } else {
      printf("branch electron_full5x5_sigmaIEtaIEta_branch does not exist!\n");
      exit(1);
    }
    electron_full5x5_sigmaIEtaIEta_isLoaded = true;
  }
  return electron_full5x5_sigmaIEtaIEta_;
}
const float &SkimTree::electron_full5x5_sigmaIPhiIPhi() {
  if (not electron_full5x5_sigmaIPhiIPhi_isLoaded) {
    if (electron_full5x5_sigmaIPhiIPhi_branch != 0) {
      electron_full5x5_sigmaIPhiIPhi_branch->GetEntry(index);
    } else {
      printf("branch electron_full5x5_sigmaIPhiIPhi_branch does not exist!\n");
      exit(1);
    }
    electron_full5x5_sigmaIPhiIPhi_isLoaded = true;
  }
  return electron_full5x5_sigmaIPhiIPhi_;
}
const float &SkimTree::electron_seedTime() {
  if (not electron_seedTime_isLoaded) {
    if (electron_seedTime_branch != 0) {
      electron_seedTime_branch->GetEntry(index);
    } else {
      printf("branch electron_seedTime_branch does not exist!\n");
      exit(1);
    }
    electron_seedTime_isLoaded = true;
  }
  return electron_seedTime_;
}
const vector<float> &SkimTree::electrons_full5x5_r9() {
  if (not electrons_full5x5_r9_isLoaded) {
    if (electrons_full5x5_r9_branch != 0) {
      electrons_full5x5_r9_branch->GetEntry(index);
    } else {
      printf("branch electrons_full5x5_r9_branch does not exist!\n");
      exit(1);
    }
    electrons_full5x5_r9_isLoaded = true;
  }
  return *electrons_full5x5_r9_;
}
const vector<float> &SkimTree::electrons_full5x5_sigmaIEtaIEta() {
  if (not electrons_full5x5_sigmaIEtaIEta_isLoaded) {
    if (electrons_full5x5_sigmaIEtaIEta_branch != 0) {
      electrons_full5x5_sigmaIEtaIEta_branch->GetEntry(index);
    } else {
      printf("branch electrons_full5x5_sigmaIEtaIEta_branch does not exist!\n");
      exit(1);
    }
    electrons_full5x5_sigmaIEtaIEta_isLoaded = true;
  }
  return *electrons_full5x5_sigmaIEtaIEta_;
}
const vector<float> &SkimTree::electrons_full5x5_sigmaIPhiIPhi() {
  if (not electrons_full5x5_sigmaIPhiIPhi_isLoaded) {
    if (electrons_full5x5_sigmaIPhiIPhi_branch != 0) {
      electrons_full5x5_sigmaIPhiIPhi_branch->GetEntry(index);
    } else {
      printf("branch electrons_full5x5_sigmaIPhiIPhi_branch does not exist!\n");
      exit(1);
    }
    electrons_full5x5_sigmaIPhiIPhi_isLoaded = true;
  }
  return *electrons_full5x5_sigmaIPhiIPhi_;
}
const vector<float> &SkimTree::electrons_seedTime() {
  if (not electrons_seedTime_isLoaded) {
    if (electrons_seedTime_branch != 0) {
      electrons_seedTime_branch->GetEntry(index);
    } else {
      printf("branch electrons_seedTime_branch does not exist!\n");
      exit(1);
    }
    electrons_seedTime_isLoaded = true;
  }
  return *electrons_seedTime_;
}
const float &SkimTree::event_mTZZ() {
  if (not event_mTZZ_isLoaded) {
    if (event_mTZZ_branch != 0) {
      event_mTZZ_branch->GetEntry(index);
    } else {
      printf("branch event_mTZZ_branch does not exist!\n");
      exit(1);
    }
    event_mTZZ_isLoaded = true;
  }
  return event_mTZZ_;
}
const float &SkimTree::event_mZZ() {
  if (not event_mZZ_isLoaded) {
    if (event_mZZ_branch != 0) {
      event_mZZ_branch->GetEntry(index);
    } else {
      printf("branch event_mZZ_branch does not exist!\n");
      exit(1);
    }
    event_mZZ_isLoaded = true;
  }
  return event_mZZ_;
}
const float &SkimTree::event_mlg() {
  if (not event_mlg_isLoaded) {
    if (event_mlg_branch != 0) {
      event_mlg_branch->GetEntry(index);
    } else {
      printf("branch event_mlg_branch does not exist!\n");
      exit(1);
    }
    event_mlg_isLoaded = true;
  }
  return event_mlg_;
}
const float &SkimTree::event_mllg() {
  if (not event_mllg_isLoaded) {
    if (event_mllg_branch != 0) {
      event_mllg_branch->GetEntry(index);
    } else {
      printf("branch event_mllg_branch does not exist!\n");
      exit(1);
    }
    event_mllg_isLoaded = true;
  }
  return event_mllg_;
}
const unsigned int &SkimTree::event_n_ak4jets_pt20() {
  if (not event_n_ak4jets_pt20_isLoaded) {
    if (event_n_ak4jets_pt20_branch != 0) {
      event_n_ak4jets_pt20_branch->GetEntry(index);
    } else {
      printf("branch event_n_ak4jets_pt20_branch does not exist!\n");
      exit(1);
    }
    event_n_ak4jets_pt20_isLoaded = true;
  }
  return event_n_ak4jets_pt20_;
}
const unsigned int &SkimTree::event_n_ak4jets_pt20_btagged_loose() {
  if (not event_n_ak4jets_pt20_btagged_loose_isLoaded) {
    if (event_n_ak4jets_pt20_btagged_loose_branch != 0) {
      event_n_ak4jets_pt20_btagged_loose_branch->GetEntry(index);
    } else {
      printf("branch event_n_ak4jets_pt20_btagged_loose_branch does not exist!\n");
      exit(1);
    }
    event_n_ak4jets_pt20_btagged_loose_isLoaded = true;
  }
  return event_n_ak4jets_pt20_btagged_loose_;
}
const unsigned int &SkimTree::event_n_ak4jets_pt20_btagged_medium() {
  if (not event_n_ak4jets_pt20_btagged_medium_isLoaded) {
    if (event_n_ak4jets_pt20_btagged_medium_branch != 0) {
      event_n_ak4jets_pt20_btagged_medium_branch->GetEntry(index);
    } else {
      printf("branch event_n_ak4jets_pt20_btagged_medium_branch does not exist!\n");
      exit(1);
    }
    event_n_ak4jets_pt20_btagged_medium_isLoaded = true;
  }
  return event_n_ak4jets_pt20_btagged_medium_;
}
const unsigned int &SkimTree::event_n_ak4jets_pt30() {
  if (not event_n_ak4jets_pt30_isLoaded) {
    if (event_n_ak4jets_pt30_branch != 0) {
      event_n_ak4jets_pt30_branch->GetEntry(index);
    } else {
      printf("branch event_n_ak4jets_pt30_branch does not exist!\n");
      exit(1);
    }
    event_n_ak4jets_pt30_isLoaded = true;
  }
  return event_n_ak4jets_pt30_;
}
const unsigned int &SkimTree::event_n_ak4jets_pt30_btagged_loose() {
  if (not event_n_ak4jets_pt30_btagged_loose_isLoaded) {
    if (event_n_ak4jets_pt30_btagged_loose_branch != 0) {
      event_n_ak4jets_pt30_btagged_loose_branch->GetEntry(index);
    } else {
      printf("branch event_n_ak4jets_pt30_btagged_loose_branch does not exist!\n");
      exit(1);
    }
    event_n_ak4jets_pt30_btagged_loose_isLoaded = true;
  }
  return event_n_ak4jets_pt30_btagged_loose_;
}
const unsigned int &SkimTree::event_n_ak4jets_pt30_btagged_medium() {
  if (not event_n_ak4jets_pt30_btagged_medium_isLoaded) {
    if (event_n_ak4jets_pt30_btagged_medium_branch != 0) {
      event_n_ak4jets_pt30_btagged_medium_branch->GetEntry(index);
    } else {
      printf("branch event_n_ak4jets_pt30_btagged_medium_branch does not exist!\n");
      exit(1);
    }
    event_n_ak4jets_pt30_btagged_medium_isLoaded = true;
  }
  return event_n_ak4jets_pt30_btagged_medium_;
}
const unsigned int &SkimTree::event_n_leptons_fakeableBase() {
  if (not event_n_leptons_fakeableBase_isLoaded) {
    if (event_n_leptons_fakeableBase_branch != 0) {
      event_n_leptons_fakeableBase_branch->GetEntry(index);
    } else {
      printf("branch event_n_leptons_fakeableBase_branch does not exist!\n");
      exit(1);
    }
    event_n_leptons_fakeableBase_isLoaded = true;
  }
  return event_n_leptons_fakeableBase_;
}
const unsigned int &SkimTree::event_n_vtxs_good() {
  if (not event_n_vtxs_good_isLoaded) {
    if (event_n_vtxs_good_branch != 0) {
      event_n_vtxs_good_branch->GetEntry(index);
    } else {
      printf("branch event_n_vtxs_good_branch does not exist!\n");
      exit(1);
    }
    event_n_vtxs_good_isLoaded = true;
  }
  return event_n_vtxs_good_;
}
const float &SkimTree::event_pTmiss() {
  if (not event_pTmiss_isLoaded) {
    if (event_pTmiss_branch != 0) {
      event_pTmiss_branch->GetEntry(index);
    } else {
      printf("branch event_pTmiss_branch does not exist!\n");
      exit(1);
    }
    event_pTmiss_isLoaded = true;
  }
  return event_pTmiss_;
}
const bool &SkimTree::event_pass_tightMETFilters() {
  if (not event_pass_tightMETFilters_isLoaded) {
    if (event_pass_tightMETFilters_branch != 0) {
      event_pass_tightMETFilters_branch->GetEntry(index);
    } else {
      printf("branch event_pass_tightMETFilters_branch does not exist!\n");
      exit(1);
    }
    event_pass_tightMETFilters_isLoaded = true;
  }
  return event_pass_tightMETFilters_;
}
const float &SkimTree::event_phimiss() {
  if (not event_phimiss_isLoaded) {
    if (event_phimiss_branch != 0) {
      event_phimiss_branch->GetEntry(index);
    } else {
      printf("branch event_phimiss_branch does not exist!\n");
      exit(1);
    }
    event_phimiss_isLoaded = true;
  }
  return event_phimiss_;
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
const float &SkimTree::event_wgt_L1PrefiringDn() {
  if (not event_wgt_L1PrefiringDn_isLoaded) {
    if (event_wgt_L1PrefiringDn_branch != 0) {
      event_wgt_L1PrefiringDn_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_L1PrefiringDn_branch does not exist!\n");
      exit(1);
    }
    event_wgt_L1PrefiringDn_isLoaded = true;
  }
  return event_wgt_L1PrefiringDn_;
}
const float &SkimTree::event_wgt_L1PrefiringUp() {
  if (not event_wgt_L1PrefiringUp_isLoaded) {
    if (event_wgt_L1PrefiringUp_branch != 0) {
      event_wgt_L1PrefiringUp_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_L1PrefiringUp_branch does not exist!\n");
      exit(1);
    }
    event_wgt_L1PrefiringUp_isLoaded = true;
  }
  return event_wgt_L1PrefiringUp_;
}
const float &SkimTree::event_wgt_PUDn() {
  if (not event_wgt_PUDn_isLoaded) {
    if (event_wgt_PUDn_branch != 0) {
      event_wgt_PUDn_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_PUDn_branch does not exist!\n");
      exit(1);
    }
    event_wgt_PUDn_isLoaded = true;
  }
  return event_wgt_PUDn_;
}
const float &SkimTree::event_wgt_PUUp() {
  if (not event_wgt_PUUp_isLoaded) {
    if (event_wgt_PUUp_branch != 0) {
      event_wgt_PUUp_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_PUUp_branch does not exist!\n");
      exit(1);
    }
    event_wgt_PUUp_isLoaded = true;
  }
  return event_wgt_PUUp_;
}
const float &SkimTree::event_wgt_SFs_PUJetId() {
  if (not event_wgt_SFs_PUJetId_isLoaded) {
    if (event_wgt_SFs_PUJetId_branch != 0) {
      event_wgt_SFs_PUJetId_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_PUJetId_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_PUJetId_isLoaded = true;
  }
  return event_wgt_SFs_PUJetId_;
}
const float &SkimTree::event_wgt_SFs_PUJetId_EffDn() {
  if (not event_wgt_SFs_PUJetId_EffDn_isLoaded) {
    if (event_wgt_SFs_PUJetId_EffDn_branch != 0) {
      event_wgt_SFs_PUJetId_EffDn_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_PUJetId_EffDn_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_PUJetId_EffDn_isLoaded = true;
  }
  return event_wgt_SFs_PUJetId_EffDn_;
}
const float &SkimTree::event_wgt_SFs_PUJetId_EffUp() {
  if (not event_wgt_SFs_PUJetId_EffUp_isLoaded) {
    if (event_wgt_SFs_PUJetId_EffUp_branch != 0) {
      event_wgt_SFs_PUJetId_EffUp_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_PUJetId_EffUp_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_PUJetId_EffUp_isLoaded = true;
  }
  return event_wgt_SFs_PUJetId_EffUp_;
}
const float &SkimTree::event_wgt_SFs_btagging() {
  if (not event_wgt_SFs_btagging_isLoaded) {
    if (event_wgt_SFs_btagging_branch != 0) {
      event_wgt_SFs_btagging_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_btagging_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_btagging_isLoaded = true;
  }
  return event_wgt_SFs_btagging_;
}
const float &SkimTree::event_wgt_SFs_btagging_EffDn() {
  if (not event_wgt_SFs_btagging_EffDn_isLoaded) {
    if (event_wgt_SFs_btagging_EffDn_branch != 0) {
      event_wgt_SFs_btagging_EffDn_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_btagging_EffDn_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_btagging_EffDn_isLoaded = true;
  }
  return event_wgt_SFs_btagging_EffDn_;
}
const float &SkimTree::event_wgt_SFs_btagging_EffUp() {
  if (not event_wgt_SFs_btagging_EffUp_isLoaded) {
    if (event_wgt_SFs_btagging_EffUp_branch != 0) {
      event_wgt_SFs_btagging_EffUp_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_btagging_EffUp_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_btagging_EffUp_isLoaded = true;
  }
  return event_wgt_SFs_btagging_EffUp_;
}
const float &SkimTree::event_wgt_SFs_electrons() {
  if (not event_wgt_SFs_electrons_isLoaded) {
    if (event_wgt_SFs_electrons_branch != 0) {
      event_wgt_SFs_electrons_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_electrons_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_electrons_isLoaded = true;
  }
  return event_wgt_SFs_electrons_;
}
const float &SkimTree::event_wgt_SFs_electrons_AltMCDn() {
  if (not event_wgt_SFs_electrons_AltMCDn_isLoaded) {
    if (event_wgt_SFs_electrons_AltMCDn_branch != 0) {
      event_wgt_SFs_electrons_AltMCDn_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_electrons_AltMCDn_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_electrons_AltMCDn_isLoaded = true;
  }
  return event_wgt_SFs_electrons_AltMCDn_;
}
const float &SkimTree::event_wgt_SFs_electrons_AltMCUp() {
  if (not event_wgt_SFs_electrons_AltMCUp_isLoaded) {
    if (event_wgt_SFs_electrons_AltMCUp_branch != 0) {
      event_wgt_SFs_electrons_AltMCUp_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_electrons_AltMCUp_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_electrons_AltMCUp_isLoaded = true;
  }
  return event_wgt_SFs_electrons_AltMCUp_;
}
const float &SkimTree::event_wgt_SFs_electrons_StatDn() {
  if (not event_wgt_SFs_electrons_StatDn_isLoaded) {
    if (event_wgt_SFs_electrons_StatDn_branch != 0) {
      event_wgt_SFs_electrons_StatDn_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_electrons_StatDn_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_electrons_StatDn_isLoaded = true;
  }
  return event_wgt_SFs_electrons_StatDn_;
}
const float &SkimTree::event_wgt_SFs_electrons_StatUp() {
  if (not event_wgt_SFs_electrons_StatUp_isLoaded) {
    if (event_wgt_SFs_electrons_StatUp_branch != 0) {
      event_wgt_SFs_electrons_StatUp_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_electrons_StatUp_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_electrons_StatUp_isLoaded = true;
  }
  return event_wgt_SFs_electrons_StatUp_;
}
const float &SkimTree::event_wgt_SFs_electrons_SystDn() {
  if (not event_wgt_SFs_electrons_SystDn_isLoaded) {
    if (event_wgt_SFs_electrons_SystDn_branch != 0) {
      event_wgt_SFs_electrons_SystDn_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_electrons_SystDn_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_electrons_SystDn_isLoaded = true;
  }
  return event_wgt_SFs_electrons_SystDn_;
}
const float &SkimTree::event_wgt_SFs_electrons_SystUp() {
  if (not event_wgt_SFs_electrons_SystUp_isLoaded) {
    if (event_wgt_SFs_electrons_SystUp_branch != 0) {
      event_wgt_SFs_electrons_SystUp_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_electrons_SystUp_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_electrons_SystUp_isLoaded = true;
  }
  return event_wgt_SFs_electrons_SystUp_;
}
const float &SkimTree::event_wgt_SFs_muons() {
  if (not event_wgt_SFs_muons_isLoaded) {
    if (event_wgt_SFs_muons_branch != 0) {
      event_wgt_SFs_muons_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_muons_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_muons_isLoaded = true;
  }
  return event_wgt_SFs_muons_;
}
const float &SkimTree::event_wgt_SFs_muons_AltMCDn() {
  if (not event_wgt_SFs_muons_AltMCDn_isLoaded) {
    if (event_wgt_SFs_muons_AltMCDn_branch != 0) {
      event_wgt_SFs_muons_AltMCDn_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_muons_AltMCDn_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_muons_AltMCDn_isLoaded = true;
  }
  return event_wgt_SFs_muons_AltMCDn_;
}
const float &SkimTree::event_wgt_SFs_muons_AltMCUp() {
  if (not event_wgt_SFs_muons_AltMCUp_isLoaded) {
    if (event_wgt_SFs_muons_AltMCUp_branch != 0) {
      event_wgt_SFs_muons_AltMCUp_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_muons_AltMCUp_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_muons_AltMCUp_isLoaded = true;
  }
  return event_wgt_SFs_muons_AltMCUp_;
}
const float &SkimTree::event_wgt_SFs_muons_StatDn() {
  if (not event_wgt_SFs_muons_StatDn_isLoaded) {
    if (event_wgt_SFs_muons_StatDn_branch != 0) {
      event_wgt_SFs_muons_StatDn_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_muons_StatDn_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_muons_StatDn_isLoaded = true;
  }
  return event_wgt_SFs_muons_StatDn_;
}
const float &SkimTree::event_wgt_SFs_muons_StatUp() {
  if (not event_wgt_SFs_muons_StatUp_isLoaded) {
    if (event_wgt_SFs_muons_StatUp_branch != 0) {
      event_wgt_SFs_muons_StatUp_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_muons_StatUp_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_muons_StatUp_isLoaded = true;
  }
  return event_wgt_SFs_muons_StatUp_;
}
const float &SkimTree::event_wgt_SFs_muons_SystDn() {
  if (not event_wgt_SFs_muons_SystDn_isLoaded) {
    if (event_wgt_SFs_muons_SystDn_branch != 0) {
      event_wgt_SFs_muons_SystDn_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_muons_SystDn_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_muons_SystDn_isLoaded = true;
  }
  return event_wgt_SFs_muons_SystDn_;
}
const float &SkimTree::event_wgt_SFs_muons_SystUp() {
  if (not event_wgt_SFs_muons_SystUp_isLoaded) {
    if (event_wgt_SFs_muons_SystUp_branch != 0) {
      event_wgt_SFs_muons_SystUp_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_muons_SystUp_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_muons_SystUp_isLoaded = true;
  }
  return event_wgt_SFs_muons_SystUp_;
}
const float &SkimTree::event_wgt_SFs_photons() {
  if (not event_wgt_SFs_photons_isLoaded) {
    if (event_wgt_SFs_photons_branch != 0) {
      event_wgt_SFs_photons_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_photons_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_photons_isLoaded = true;
  }
  return event_wgt_SFs_photons_;
}
const float &SkimTree::event_wgt_SFs_photons_EffDn() {
  if (not event_wgt_SFs_photons_EffDn_isLoaded) {
    if (event_wgt_SFs_photons_EffDn_branch != 0) {
      event_wgt_SFs_photons_EffDn_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_photons_EffDn_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_photons_EffDn_isLoaded = true;
  }
  return event_wgt_SFs_photons_EffDn_;
}
const float &SkimTree::event_wgt_SFs_photons_EffUp() {
  if (not event_wgt_SFs_photons_EffUp_isLoaded) {
    if (event_wgt_SFs_photons_EffUp_branch != 0) {
      event_wgt_SFs_photons_EffUp_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_SFs_photons_EffUp_branch does not exist!\n");
      exit(1);
    }
    event_wgt_SFs_photons_EffUp_isLoaded = true;
  }
  return event_wgt_SFs_photons_EffUp_;
}
const float &SkimTree::event_wgt_adjustment_AsMZDn() {
  if (not event_wgt_adjustment_AsMZDn_isLoaded) {
    if (event_wgt_adjustment_AsMZDn_branch != 0) {
      event_wgt_adjustment_AsMZDn_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_adjustment_AsMZDn_branch does not exist!\n");
      exit(1);
    }
    event_wgt_adjustment_AsMZDn_isLoaded = true;
  }
  return event_wgt_adjustment_AsMZDn_;
}
const float &SkimTree::event_wgt_adjustment_AsMZUp() {
  if (not event_wgt_adjustment_AsMZUp_isLoaded) {
    if (event_wgt_adjustment_AsMZUp_branch != 0) {
      event_wgt_adjustment_AsMZUp_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_adjustment_AsMZUp_branch does not exist!\n");
      exit(1);
    }
    event_wgt_adjustment_AsMZUp_isLoaded = true;
  }
  return event_wgt_adjustment_AsMZUp_;
}
const float &SkimTree::event_wgt_adjustment_NNPDF30() {
  if (not event_wgt_adjustment_NNPDF30_isLoaded) {
    if (event_wgt_adjustment_NNPDF30_branch != 0) {
      event_wgt_adjustment_NNPDF30_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_adjustment_NNPDF30_branch does not exist!\n");
      exit(1);
    }
    event_wgt_adjustment_NNPDF30_isLoaded = true;
  }
  return event_wgt_adjustment_NNPDF30_;
}
const float &SkimTree::event_wgt_adjustment_NNPDF30_AsMZDn() {
  if (not event_wgt_adjustment_NNPDF30_AsMZDn_isLoaded) {
    if (event_wgt_adjustment_NNPDF30_AsMZDn_branch != 0) {
      event_wgt_adjustment_NNPDF30_AsMZDn_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_adjustment_NNPDF30_AsMZDn_branch does not exist!\n");
      exit(1);
    }
    event_wgt_adjustment_NNPDF30_AsMZDn_isLoaded = true;
  }
  return event_wgt_adjustment_NNPDF30_AsMZDn_;
}
const float &SkimTree::event_wgt_adjustment_NNPDF30_AsMZUp() {
  if (not event_wgt_adjustment_NNPDF30_AsMZUp_isLoaded) {
    if (event_wgt_adjustment_NNPDF30_AsMZUp_branch != 0) {
      event_wgt_adjustment_NNPDF30_AsMZUp_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_adjustment_NNPDF30_AsMZUp_branch does not exist!\n");
      exit(1);
    }
    event_wgt_adjustment_NNPDF30_AsMZUp_isLoaded = true;
  }
  return event_wgt_adjustment_NNPDF30_AsMZUp_;
}
const float &SkimTree::event_wgt_adjustment_NNPDF30_PDFReplicaDn() {
  if (not event_wgt_adjustment_NNPDF30_PDFReplicaDn_isLoaded) {
    if (event_wgt_adjustment_NNPDF30_PDFReplicaDn_branch != 0) {
      event_wgt_adjustment_NNPDF30_PDFReplicaDn_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_adjustment_NNPDF30_PDFReplicaDn_branch does not exist!\n");
      exit(1);
    }
    event_wgt_adjustment_NNPDF30_PDFReplicaDn_isLoaded = true;
  }
  return event_wgt_adjustment_NNPDF30_PDFReplicaDn_;
}
const float &SkimTree::event_wgt_adjustment_NNPDF30_PDFReplicaUp() {
  if (not event_wgt_adjustment_NNPDF30_PDFReplicaUp_isLoaded) {
    if (event_wgt_adjustment_NNPDF30_PDFReplicaUp_branch != 0) {
      event_wgt_adjustment_NNPDF30_PDFReplicaUp_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_adjustment_NNPDF30_PDFReplicaUp_branch does not exist!\n");
      exit(1);
    }
    event_wgt_adjustment_NNPDF30_PDFReplicaUp_isLoaded = true;
  }
  return event_wgt_adjustment_NNPDF30_PDFReplicaUp_;
}
const float &SkimTree::event_wgt_adjustment_PDFReplicaDn() {
  if (not event_wgt_adjustment_PDFReplicaDn_isLoaded) {
    if (event_wgt_adjustment_PDFReplicaDn_branch != 0) {
      event_wgt_adjustment_PDFReplicaDn_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_adjustment_PDFReplicaDn_branch does not exist!\n");
      exit(1);
    }
    event_wgt_adjustment_PDFReplicaDn_isLoaded = true;
  }
  return event_wgt_adjustment_PDFReplicaDn_;
}
const float &SkimTree::event_wgt_adjustment_PDFReplicaUp() {
  if (not event_wgt_adjustment_PDFReplicaUp_isLoaded) {
    if (event_wgt_adjustment_PDFReplicaUp_branch != 0) {
      event_wgt_adjustment_PDFReplicaUp_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_adjustment_PDFReplicaUp_branch does not exist!\n");
      exit(1);
    }
    event_wgt_adjustment_PDFReplicaUp_isLoaded = true;
  }
  return event_wgt_adjustment_PDFReplicaUp_;
}
const float &SkimTree::event_wgt_adjustment_PDFScaleDn() {
  if (not event_wgt_adjustment_PDFScaleDn_isLoaded) {
    if (event_wgt_adjustment_PDFScaleDn_branch != 0) {
      event_wgt_adjustment_PDFScaleDn_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_adjustment_PDFScaleDn_branch does not exist!\n");
      exit(1);
    }
    event_wgt_adjustment_PDFScaleDn_isLoaded = true;
  }
  return event_wgt_adjustment_PDFScaleDn_;
}
const float &SkimTree::event_wgt_adjustment_PDFScaleUp() {
  if (not event_wgt_adjustment_PDFScaleUp_isLoaded) {
    if (event_wgt_adjustment_PDFScaleUp_branch != 0) {
      event_wgt_adjustment_PDFScaleUp_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_adjustment_PDFScaleUp_branch does not exist!\n");
      exit(1);
    }
    event_wgt_adjustment_PDFScaleUp_isLoaded = true;
  }
  return event_wgt_adjustment_PDFScaleUp_;
}
const float &SkimTree::event_wgt_adjustment_PythiaScaleDn() {
  if (not event_wgt_adjustment_PythiaScaleDn_isLoaded) {
    if (event_wgt_adjustment_PythiaScaleDn_branch != 0) {
      event_wgt_adjustment_PythiaScaleDn_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_adjustment_PythiaScaleDn_branch does not exist!\n");
      exit(1);
    }
    event_wgt_adjustment_PythiaScaleDn_isLoaded = true;
  }
  return event_wgt_adjustment_PythiaScaleDn_;
}
const float &SkimTree::event_wgt_adjustment_PythiaScaleUp() {
  if (not event_wgt_adjustment_PythiaScaleUp_isLoaded) {
    if (event_wgt_adjustment_PythiaScaleUp_branch != 0) {
      event_wgt_adjustment_PythiaScaleUp_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_adjustment_PythiaScaleUp_branch does not exist!\n");
      exit(1);
    }
    event_wgt_adjustment_PythiaScaleUp_isLoaded = true;
  }
  return event_wgt_adjustment_PythiaScaleUp_;
}
const float &SkimTree::event_wgt_adjustment_QCDScaleDn() {
  if (not event_wgt_adjustment_QCDScaleDn_isLoaded) {
    if (event_wgt_adjustment_QCDScaleDn_branch != 0) {
      event_wgt_adjustment_QCDScaleDn_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_adjustment_QCDScaleDn_branch does not exist!\n");
      exit(1);
    }
    event_wgt_adjustment_QCDScaleDn_isLoaded = true;
  }
  return event_wgt_adjustment_QCDScaleDn_;
}
const float &SkimTree::event_wgt_adjustment_QCDScaleUp() {
  if (not event_wgt_adjustment_QCDScaleUp_isLoaded) {
    if (event_wgt_adjustment_QCDScaleUp_branch != 0) {
      event_wgt_adjustment_QCDScaleUp_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_adjustment_QCDScaleUp_branch does not exist!\n");
      exit(1);
    }
    event_wgt_adjustment_QCDScaleUp_isLoaded = true;
  }
  return event_wgt_adjustment_QCDScaleUp_;
}
const float &SkimTree::event_wgt_triggers() {
  if (not event_wgt_triggers_isLoaded) {
    if (event_wgt_triggers_branch != 0) {
      event_wgt_triggers_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_triggers_branch does not exist!\n");
      exit(1);
    }
    event_wgt_triggers_isLoaded = true;
  }
  return event_wgt_triggers_;
}
const float &SkimTree::event_wgt_triggers_Dilepton() {
  if (not event_wgt_triggers_Dilepton_isLoaded) {
    if (event_wgt_triggers_Dilepton_branch != 0) {
      event_wgt_triggers_Dilepton_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_triggers_Dilepton_branch does not exist!\n");
      exit(1);
    }
    event_wgt_triggers_Dilepton_isLoaded = true;
  }
  return event_wgt_triggers_Dilepton_;
}
const float &SkimTree::event_wgt_triggers_Dilepton_DF_Extra() {
  if (not event_wgt_triggers_Dilepton_DF_Extra_isLoaded) {
    if (event_wgt_triggers_Dilepton_DF_Extra_branch != 0) {
      event_wgt_triggers_Dilepton_DF_Extra_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_triggers_Dilepton_DF_Extra_branch does not exist!\n");
      exit(1);
    }
    event_wgt_triggers_Dilepton_DF_Extra_isLoaded = true;
  }
  return event_wgt_triggers_Dilepton_DF_Extra_;
}
const float &SkimTree::event_wgt_triggers_PFHT_Control() {
  if (not event_wgt_triggers_PFHT_Control_isLoaded) {
    if (event_wgt_triggers_PFHT_Control_branch != 0) {
      event_wgt_triggers_PFHT_Control_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_triggers_PFHT_Control_branch does not exist!\n");
      exit(1);
    }
    event_wgt_triggers_PFHT_Control_isLoaded = true;
  }
  return event_wgt_triggers_PFHT_Control_;
}
const float &SkimTree::event_wgt_triggers_PFMET_MHT_Control() {
  if (not event_wgt_triggers_PFMET_MHT_Control_isLoaded) {
    if (event_wgt_triggers_PFMET_MHT_Control_branch != 0) {
      event_wgt_triggers_PFMET_MHT_Control_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_triggers_PFMET_MHT_Control_branch does not exist!\n");
      exit(1);
    }
    event_wgt_triggers_PFMET_MHT_Control_isLoaded = true;
  }
  return event_wgt_triggers_PFMET_MHT_Control_;
}
const float &SkimTree::event_wgt_triggers_SingleLepton() {
  if (not event_wgt_triggers_SingleLepton_isLoaded) {
    if (event_wgt_triggers_SingleLepton_branch != 0) {
      event_wgt_triggers_SingleLepton_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_triggers_SingleLepton_branch does not exist!\n");
      exit(1);
    }
    event_wgt_triggers_SingleLepton_isLoaded = true;
  }
  return event_wgt_triggers_SingleLepton_;
}
const float &SkimTree::event_wgt_triggers_SinglePhoton() {
  if (not event_wgt_triggers_SinglePhoton_isLoaded) {
    if (event_wgt_triggers_SinglePhoton_branch != 0) {
      event_wgt_triggers_SinglePhoton_branch->GetEntry(index);
    } else {
      printf("branch event_wgt_triggers_SinglePhoton_branch does not exist!\n");
      exit(1);
    }
    event_wgt_triggers_SinglePhoton_isLoaded = true;
  }
  return event_wgt_triggers_SinglePhoton_;
}
const float &SkimTree::genmet_pTmiss() {
  if (not genmet_pTmiss_isLoaded) {
    if (genmet_pTmiss_branch != 0) {
      genmet_pTmiss_branch->GetEntry(index);
    } else {
      printf("branch genmet_pTmiss_branch does not exist!\n");
      exit(1);
    }
    genmet_pTmiss_isLoaded = true;
  }
  return genmet_pTmiss_;
}
const float &SkimTree::genmet_phimiss() {
  if (not genmet_phimiss_isLoaded) {
    if (genmet_phimiss_branch != 0) {
      genmet_phimiss_branch->GetEntry(index);
    } else {
      printf("branch genmet_phimiss_branch does not exist!\n");
      exit(1);
    }
    genmet_phimiss_isLoaded = true;
  }
  return genmet_phimiss_;
}
const float &SkimTree::lepton_eff() {
  if (not lepton_eff_isLoaded) {
    if (lepton_eff_branch != 0) {
      lepton_eff_branch->GetEntry(index);
    } else {
      printf("branch lepton_eff_branch does not exist!\n");
      exit(1);
    }
    lepton_eff_isLoaded = true;
  }
  return lepton_eff_;
}
const float &SkimTree::lepton_eff_StatDn() {
  if (not lepton_eff_StatDn_isLoaded) {
    if (lepton_eff_StatDn_branch != 0) {
      lepton_eff_StatDn_branch->GetEntry(index);
    } else {
      printf("branch lepton_eff_StatDn_branch does not exist!\n");
      exit(1);
    }
    lepton_eff_StatDn_isLoaded = true;
  }
  return lepton_eff_StatDn_;
}
const float &SkimTree::lepton_eff_StatUp() {
  if (not lepton_eff_StatUp_isLoaded) {
    if (lepton_eff_StatUp_branch != 0) {
      lepton_eff_StatUp_branch->GetEntry(index);
    } else {
      printf("branch lepton_eff_StatUp_branch does not exist!\n");
      exit(1);
    }
    lepton_eff_StatUp_isLoaded = true;
  }
  return lepton_eff_StatUp_;
}
const float &SkimTree::lepton_eff_SystDn() {
  if (not lepton_eff_SystDn_isLoaded) {
    if (lepton_eff_SystDn_branch != 0) {
      lepton_eff_SystDn_branch->GetEntry(index);
    } else {
      printf("branch lepton_eff_SystDn_branch does not exist!\n");
      exit(1);
    }
    lepton_eff_SystDn_isLoaded = true;
  }
  return lepton_eff_SystDn_;
}
const float &SkimTree::lepton_eff_SystUp() {
  if (not lepton_eff_SystUp_isLoaded) {
    if (lepton_eff_SystUp_branch != 0) {
      lepton_eff_SystUp_branch->GetEntry(index);
    } else {
      printf("branch lepton_eff_SystUp_branch does not exist!\n");
      exit(1);
    }
    lepton_eff_SystUp_isLoaded = true;
  }
  return lepton_eff_SystUp_;
}
const float &SkimTree::lepton_eta() {
  if (not lepton_eta_isLoaded) {
    if (lepton_eta_branch != 0) {
      lepton_eta_branch->GetEntry(index);
    } else {
      printf("branch lepton_eta_branch does not exist!\n");
      exit(1);
    }
    lepton_eta_isLoaded = true;
  }
  return lepton_eta_;
}
const int &SkimTree::lepton_id() {
  if (not lepton_id_isLoaded) {
    if (lepton_id_branch != 0) {
      lepton_id_branch->GetEntry(index);
    } else {
      printf("branch lepton_id_branch does not exist!\n");
      exit(1);
    }
    lepton_id_isLoaded = true;
  }
  return lepton_id_;
}
const bool &SkimTree::lepton_is_genMatched_prompt() {
  if (not lepton_is_genMatched_prompt_isLoaded) {
    if (lepton_is_genMatched_prompt_branch != 0) {
      lepton_is_genMatched_prompt_branch->GetEntry(index);
    } else {
      printf("branch lepton_is_genMatched_prompt_branch does not exist!\n");
      exit(1);
    }
    lepton_is_genMatched_prompt_isLoaded = true;
  }
  return lepton_is_genMatched_prompt_;
}
const float &SkimTree::lepton_mass() {
  if (not lepton_mass_isLoaded) {
    if (lepton_mass_branch != 0) {
      lepton_mass_branch->GetEntry(index);
    } else {
      printf("branch lepton_mass_branch does not exist!\n");
      exit(1);
    }
    lepton_mass_isLoaded = true;
  }
  return lepton_mass_;
}
const float &SkimTree::lepton_phi() {
  if (not lepton_phi_isLoaded) {
    if (lepton_phi_branch != 0) {
      lepton_phi_branch->GetEntry(index);
    } else {
      printf("branch lepton_phi_branch does not exist!\n");
      exit(1);
    }
    lepton_phi_isLoaded = true;
  }
  return lepton_phi_;
}
const float &SkimTree::lepton_pt() {
  if (not lepton_pt_isLoaded) {
    if (lepton_pt_branch != 0) {
      lepton_pt_branch->GetEntry(index);
    } else {
      printf("branch lepton_pt_branch does not exist!\n");
      exit(1);
    }
    lepton_pt_isLoaded = true;
  }
  return lepton_pt_;
}
const vector<float> &SkimTree::genak4jets_eta() {
  if (not genak4jets_eta_isLoaded) {
    if (genak4jets_eta_branch != 0) {
      genak4jets_eta_branch->GetEntry(index);
    } else {
      printf("branch genak4jets_eta_branch does not exist!\n");
      exit(1);
    }
    genak4jets_eta_isLoaded = true;
  }
  return *genak4jets_eta_;
}
const vector<float> &SkimTree::genak4jets_mass() {
  if (not genak4jets_mass_isLoaded) {
    if (genak4jets_mass_branch != 0) {
      genak4jets_mass_branch->GetEntry(index);
    } else {
      printf("branch genak4jets_mass_branch does not exist!\n");
      exit(1);
    }
    genak4jets_mass_isLoaded = true;
  }
  return *genak4jets_mass_;
}
const vector<float> &SkimTree::genak4jets_phi() {
  if (not genak4jets_phi_isLoaded) {
    if (genak4jets_phi_branch != 0) {
      genak4jets_phi_branch->GetEntry(index);
    } else {
      printf("branch genak4jets_phi_branch does not exist!\n");
      exit(1);
    }
    genak4jets_phi_isLoaded = true;
  }
  return *genak4jets_phi_;
}
const vector<float> &SkimTree::genak4jets_pt() {
  if (not genak4jets_pt_isLoaded) {
    if (genak4jets_pt_branch != 0) {
      genak4jets_pt_branch->GetEntry(index);
    } else {
      printf("branch genak4jets_pt_branch does not exist!\n");
      exit(1);
    }
    genak4jets_pt_isLoaded = true;
  }
  return *genak4jets_pt_;
}
const float &SkimTree::genpromptparticles_sump4_mass() {
  if (not genpromptparticles_sump4_mass_isLoaded) {
    if (genpromptparticles_sump4_mass_branch != 0) {
      genpromptparticles_sump4_mass_branch->GetEntry(index);
    } else {
      printf("branch genpromptparticles_sump4_mass_branch does not exist!\n");
      exit(1);
    }
    genpromptparticles_sump4_mass_isLoaded = true;
  }
  return genpromptparticles_sump4_mass_;
}
const float &SkimTree::genpromptparticles_sump4_pt() {
  if (not genpromptparticles_sump4_pt_isLoaded) {
    if (genpromptparticles_sump4_pt_branch != 0) {
      genpromptparticles_sump4_pt_branch->GetEntry(index);
    } else {
      printf("branch genpromptparticles_sump4_pt_branch does not exist!\n");
      exit(1);
    }
    genpromptparticles_sump4_pt_isLoaded = true;
  }
  return genpromptparticles_sump4_pt_;
}
const float &SkimTree::genpromptparticles_sump4_rapidity() {
  if (not genpromptparticles_sump4_rapidity_isLoaded) {
    if (genpromptparticles_sump4_rapidity_branch != 0) {
      genpromptparticles_sump4_rapidity_branch->GetEntry(index);
    } else {
      printf("branch genpromptparticles_sump4_rapidity_branch does not exist!\n");
      exit(1);
    }
    genpromptparticles_sump4_rapidity_isLoaded = true;
  }
  return genpromptparticles_sump4_rapidity_;
}
const vector<float> &SkimTree::leptons_eff() {
  if (not leptons_eff_isLoaded) {
    if (leptons_eff_branch != 0) {
      leptons_eff_branch->GetEntry(index);
    } else {
      printf("branch leptons_eff_branch does not exist!\n");
      exit(1);
    }
    leptons_eff_isLoaded = true;
  }
  return *leptons_eff_;
}
const vector<float> &SkimTree::leptons_eff_DF() {
  if (not leptons_eff_DF_isLoaded) {
    if (leptons_eff_DF_branch != 0) {
      leptons_eff_DF_branch->GetEntry(index);
    } else {
      printf("branch leptons_eff_DF_branch does not exist!\n");
      exit(1);
    }
    leptons_eff_DF_isLoaded = true;
  }
  return *leptons_eff_DF_;
}
const vector<float> &SkimTree::leptons_eff_DF_StatDn() {
  if (not leptons_eff_DF_StatDn_isLoaded) {
    if (leptons_eff_DF_StatDn_branch != 0) {
      leptons_eff_DF_StatDn_branch->GetEntry(index);
    } else {
      printf("branch leptons_eff_DF_StatDn_branch does not exist!\n");
      exit(1);
    }
    leptons_eff_DF_StatDn_isLoaded = true;
  }
  return *leptons_eff_DF_StatDn_;
}
const vector<float> &SkimTree::leptons_eff_DF_StatUp() {
  if (not leptons_eff_DF_StatUp_isLoaded) {
    if (leptons_eff_DF_StatUp_branch != 0) {
      leptons_eff_DF_StatUp_branch->GetEntry(index);
    } else {
      printf("branch leptons_eff_DF_StatUp_branch does not exist!\n");
      exit(1);
    }
    leptons_eff_DF_StatUp_isLoaded = true;
  }
  return *leptons_eff_DF_StatUp_;
}
const vector<float> &SkimTree::leptons_eff_DF_SystDn() {
  if (not leptons_eff_DF_SystDn_isLoaded) {
    if (leptons_eff_DF_SystDn_branch != 0) {
      leptons_eff_DF_SystDn_branch->GetEntry(index);
    } else {
      printf("branch leptons_eff_DF_SystDn_branch does not exist!\n");
      exit(1);
    }
    leptons_eff_DF_SystDn_isLoaded = true;
  }
  return *leptons_eff_DF_SystDn_;
}
const vector<float> &SkimTree::leptons_eff_DF_SystUp() {
  if (not leptons_eff_DF_SystUp_isLoaded) {
    if (leptons_eff_DF_SystUp_branch != 0) {
      leptons_eff_DF_SystUp_branch->GetEntry(index);
    } else {
      printf("branch leptons_eff_DF_SystUp_branch does not exist!\n");
      exit(1);
    }
    leptons_eff_DF_SystUp_isLoaded = true;
  }
  return *leptons_eff_DF_SystUp_;
}
const vector<float> &SkimTree::leptons_eff_StatDn() {
  if (not leptons_eff_StatDn_isLoaded) {
    if (leptons_eff_StatDn_branch != 0) {
      leptons_eff_StatDn_branch->GetEntry(index);
    } else {
      printf("branch leptons_eff_StatDn_branch does not exist!\n");
      exit(1);
    }
    leptons_eff_StatDn_isLoaded = true;
  }
  return *leptons_eff_StatDn_;
}
const vector<float> &SkimTree::leptons_eff_StatUp() {
  if (not leptons_eff_StatUp_isLoaded) {
    if (leptons_eff_StatUp_branch != 0) {
      leptons_eff_StatUp_branch->GetEntry(index);
    } else {
      printf("branch leptons_eff_StatUp_branch does not exist!\n");
      exit(1);
    }
    leptons_eff_StatUp_isLoaded = true;
  }
  return *leptons_eff_StatUp_;
}
const vector<float> &SkimTree::leptons_eff_SystDn() {
  if (not leptons_eff_SystDn_isLoaded) {
    if (leptons_eff_SystDn_branch != 0) {
      leptons_eff_SystDn_branch->GetEntry(index);
    } else {
      printf("branch leptons_eff_SystDn_branch does not exist!\n");
      exit(1);
    }
    leptons_eff_SystDn_isLoaded = true;
  }
  return *leptons_eff_SystDn_;
}
const vector<float> &SkimTree::leptons_eff_SystUp() {
  if (not leptons_eff_SystUp_isLoaded) {
    if (leptons_eff_SystUp_branch != 0) {
      leptons_eff_SystUp_branch->GetEntry(index);
    } else {
      printf("branch leptons_eff_SystUp_branch does not exist!\n");
      exit(1);
    }
    leptons_eff_SystUp_isLoaded = true;
  }
  return *leptons_eff_SystUp_;
}
const vector<float> &SkimTree::leptons_eta() {
  if (not leptons_eta_isLoaded) {
    if (leptons_eta_branch != 0) {
      leptons_eta_branch->GetEntry(index);
    } else {
      printf("branch leptons_eta_branch does not exist!\n");
      exit(1);
    }
    leptons_eta_isLoaded = true;
  }
  return *leptons_eta_;
}
const vector<int> &SkimTree::leptons_id() {
  if (not leptons_id_isLoaded) {
    if (leptons_id_branch != 0) {
      leptons_id_branch->GetEntry(index);
    } else {
      printf("branch leptons_id_branch does not exist!\n");
      exit(1);
    }
    leptons_id_isLoaded = true;
  }
  return *leptons_id_;
}
const vector<bool> &SkimTree::leptons_is_TOmatched_SingleLepton() {
  if (not leptons_is_TOmatched_SingleLepton_isLoaded) {
    if (leptons_is_TOmatched_SingleLepton_branch != 0) {
      leptons_is_TOmatched_SingleLepton_branch->GetEntry(index);
    } else {
      printf("branch leptons_is_TOmatched_SingleLepton_branch does not exist!\n");
      exit(1);
    }
    leptons_is_TOmatched_SingleLepton_isLoaded = true;
  }
  return *leptons_is_TOmatched_SingleLepton_;
}
const vector<bool> &SkimTree::leptons_is_genMatched_prompt() {
  if (not leptons_is_genMatched_prompt_isLoaded) {
    if (leptons_is_genMatched_prompt_branch != 0) {
      leptons_is_genMatched_prompt_branch->GetEntry(index);
    } else {
      printf("branch leptons_is_genMatched_prompt_branch does not exist!\n");
      exit(1);
    }
    leptons_is_genMatched_prompt_isLoaded = true;
  }
  return *leptons_is_genMatched_prompt_;
}
const vector<float> &SkimTree::leptons_mass() {
  if (not leptons_mass_isLoaded) {
    if (leptons_mass_branch != 0) {
      leptons_mass_branch->GetEntry(index);
    } else {
      printf("branch leptons_mass_branch does not exist!\n");
      exit(1);
    }
    leptons_mass_isLoaded = true;
  }
  return *leptons_mass_;
}
const vector<float> &SkimTree::leptons_phi() {
  if (not leptons_phi_isLoaded) {
    if (leptons_phi_branch != 0) {
      leptons_phi_branch->GetEntry(index);
    } else {
      printf("branch leptons_phi_branch does not exist!\n");
      exit(1);
    }
    leptons_phi_isLoaded = true;
  }
  return *leptons_phi_;
}
const vector<float> &SkimTree::leptons_pt() {
  if (not leptons_pt_isLoaded) {
    if (leptons_pt_branch != 0) {
      leptons_pt_branch->GetEntry(index);
    } else {
      printf("branch leptons_pt_branch does not exist!\n");
      exit(1);
    }
    leptons_pt_isLoaded = true;
  }
  return *leptons_pt_;
}
const float &SkimTree::min_abs_dPhi_pTj_pTmiss() {
  if (not min_abs_dPhi_pTj_pTmiss_isLoaded) {
    if (min_abs_dPhi_pTj_pTmiss_branch != 0) {
      min_abs_dPhi_pTj_pTmiss_branch->GetEntry(index);
    } else {
      printf("branch min_abs_dPhi_pTj_pTmiss_branch does not exist!\n");
      exit(1);
    }
    min_abs_dPhi_pTj_pTmiss_isLoaded = true;
  }
  return min_abs_dPhi_pTj_pTmiss_;
}
const float &SkimTree::pConst_JJQCD_SIG_ghg2_1_JHUGen() {
  if (not pConst_JJQCD_SIG_ghg2_1_JHUGen_isLoaded) {
    if (pConst_JJQCD_SIG_ghg2_1_JHUGen_branch != 0) {
      pConst_JJQCD_SIG_ghg2_1_JHUGen_branch->GetEntry(index);
    } else {
      printf("branch pConst_JJQCD_SIG_ghg2_1_JHUGen_branch does not exist!\n");
      exit(1);
    }
    pConst_JJQCD_SIG_ghg2_1_JHUGen_isLoaded = true;
  }
  return pConst_JJQCD_SIG_ghg2_1_JHUGen_;
}
const float &SkimTree::pConst_JJVBF_SIG_ghv1_1_JHUGen() {
  if (not pConst_JJVBF_SIG_ghv1_1_JHUGen_isLoaded) {
    if (pConst_JJVBF_SIG_ghv1_1_JHUGen_branch != 0) {
      pConst_JJVBF_SIG_ghv1_1_JHUGen_branch->GetEntry(index);
    } else {
      printf("branch pConst_JJVBF_SIG_ghv1_1_JHUGen_branch does not exist!\n");
      exit(1);
    }
    pConst_JJVBF_SIG_ghv1_1_JHUGen_isLoaded = true;
  }
  return pConst_JJVBF_SIG_ghv1_1_JHUGen_;
}
const float &SkimTree::p_JJQCD_SIG_ghg2_1_JHUGen() {
  if (not p_JJQCD_SIG_ghg2_1_JHUGen_isLoaded) {
    if (p_JJQCD_SIG_ghg2_1_JHUGen_branch != 0) {
      p_JJQCD_SIG_ghg2_1_JHUGen_branch->GetEntry(index);
    } else {
      printf("branch p_JJQCD_SIG_ghg2_1_JHUGen_branch does not exist!\n");
      exit(1);
    }
    p_JJQCD_SIG_ghg2_1_JHUGen_isLoaded = true;
  }
  return p_JJQCD_SIG_ghg2_1_JHUGen_;
}
const float &SkimTree::p_JJVBF_SIG_ghv1_1_JHUGen() {
  if (not p_JJVBF_SIG_ghv1_1_JHUGen_isLoaded) {
    if (p_JJVBF_SIG_ghv1_1_JHUGen_branch != 0) {
      p_JJVBF_SIG_ghv1_1_JHUGen_branch->GetEntry(index);
    } else {
      printf("branch p_JJVBF_SIG_ghv1_1_JHUGen_branch does not exist!\n");
      exit(1);
    }
    p_JJVBF_SIG_ghv1_1_JHUGen_isLoaded = true;
  }
  return p_JJVBF_SIG_ghv1_1_JHUGen_;
}
const float &SkimTree::p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen() {
  if (not p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_isLoaded) {
    if (p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_branch != 0) {
      p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_branch->GetEntry(index);
    } else {
      printf("branch p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_branch does not exist!\n");
      exit(1);
    }
    p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_isLoaded = true;
  }
  return p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_;
}
const float &SkimTree::p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen() {
  if (not p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_isLoaded) {
    if (p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_branch != 0) {
      p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_branch->GetEntry(index);
    } else {
      printf("branch p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_branch does not exist!\n");
      exit(1);
    }
    p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_isLoaded = true;
  }
  return p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_;
}
const float &SkimTree::p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen() {
  if (not p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_isLoaded) {
    if (p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_branch != 0) {
      p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_branch->GetEntry(index);
    } else {
      printf("branch p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_branch does not exist!\n");
      exit(1);
    }
    p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_isLoaded = true;
  }
  return p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_;
}
const float &SkimTree::p_JJVBF_SIG_ghv1prime2_1E4_JHUGen() {
  if (not p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_isLoaded) {
    if (p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_branch != 0) {
      p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_branch->GetEntry(index);
    } else {
      printf("branch p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_branch does not exist!\n");
      exit(1);
    }
    p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_isLoaded = true;
  }
  return p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_;
}
const float &SkimTree::p_JJVBF_SIG_ghv2_1_JHUGen() {
  if (not p_JJVBF_SIG_ghv2_1_JHUGen_isLoaded) {
    if (p_JJVBF_SIG_ghv2_1_JHUGen_branch != 0) {
      p_JJVBF_SIG_ghv2_1_JHUGen_branch->GetEntry(index);
    } else {
      printf("branch p_JJVBF_SIG_ghv2_1_JHUGen_branch does not exist!\n");
      exit(1);
    }
    p_JJVBF_SIG_ghv2_1_JHUGen_isLoaded = true;
  }
  return p_JJVBF_SIG_ghv2_1_JHUGen_;
}
const float &SkimTree::p_JJVBF_SIG_ghv4_1_JHUGen() {
  if (not p_JJVBF_SIG_ghv4_1_JHUGen_isLoaded) {
    if (p_JJVBF_SIG_ghv4_1_JHUGen_branch != 0) {
      p_JJVBF_SIG_ghv4_1_JHUGen_branch->GetEntry(index);
    } else {
      printf("branch p_JJVBF_SIG_ghv4_1_JHUGen_branch does not exist!\n");
      exit(1);
    }
    p_JJVBF_SIG_ghv4_1_JHUGen_isLoaded = true;
  }
  return p_JJVBF_SIG_ghv4_1_JHUGen_;
}
const float &SkimTree::p_JJVBF_SIG_ghza1prime2_1E4_JHUGen() {
  if (not p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_isLoaded) {
    if (p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_branch != 0) {
      p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_branch->GetEntry(index);
    } else {
      printf("branch p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_branch does not exist!\n");
      exit(1);
    }
    p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_isLoaded = true;
  }
  return p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_;
}
const float &SkimTree::photon_MIPTotalEnergy() {
  if (not photon_MIPTotalEnergy_isLoaded) {
    if (photon_MIPTotalEnergy_branch != 0) {
      photon_MIPTotalEnergy_branch->GetEntry(index);
    } else {
      printf("branch photon_MIPTotalEnergy_branch does not exist!\n");
      exit(1);
    }
    photon_MIPTotalEnergy_isLoaded = true;
  }
  return photon_MIPTotalEnergy_;
}
const float &SkimTree::photon_eta() {
  if (not photon_eta_isLoaded) {
    if (photon_eta_branch != 0) {
      photon_eta_branch->GetEntry(index);
    } else {
      printf("branch photon_eta_branch does not exist!\n");
      exit(1);
    }
    photon_eta_isLoaded = true;
  }
  return photon_eta_;
}
const float &SkimTree::photon_full5x5_r9() {
  if (not photon_full5x5_r9_isLoaded) {
    if (photon_full5x5_r9_branch != 0) {
      photon_full5x5_r9_branch->GetEntry(index);
    } else {
      printf("branch photon_full5x5_r9_branch does not exist!\n");
      exit(1);
    }
    photon_full5x5_r9_isLoaded = true;
  }
  return photon_full5x5_r9_;
}
const float &SkimTree::photon_full5x5_sigmaIEtaIEta() {
  if (not photon_full5x5_sigmaIEtaIEta_isLoaded) {
    if (photon_full5x5_sigmaIEtaIEta_branch != 0) {
      photon_full5x5_sigmaIEtaIEta_branch->GetEntry(index);
    } else {
      printf("branch photon_full5x5_sigmaIEtaIEta_branch does not exist!\n");
      exit(1);
    }
    photon_full5x5_sigmaIEtaIEta_isLoaded = true;
  }
  return photon_full5x5_sigmaIEtaIEta_;
}
const float &SkimTree::photon_full5x5_sigmaIPhiIPhi() {
  if (not photon_full5x5_sigmaIPhiIPhi_isLoaded) {
    if (photon_full5x5_sigmaIPhiIPhi_branch != 0) {
      photon_full5x5_sigmaIPhiIPhi_branch->GetEntry(index);
    } else {
      printf("branch photon_full5x5_sigmaIPhiIPhi_branch does not exist!\n");
      exit(1);
    }
    photon_full5x5_sigmaIPhiIPhi_isLoaded = true;
  }
  return photon_full5x5_sigmaIPhiIPhi_;
}
const bool &SkimTree::photon_isEB() {
  if (not photon_isEB_isLoaded) {
    if (photon_isEB_branch != 0) {
      photon_isEB_branch->GetEntry(index);
    } else {
      printf("branch photon_isEB_branch does not exist!\n");
      exit(1);
    }
    photon_isEB_isLoaded = true;
  }
  return photon_isEB_;
}
const bool &SkimTree::photon_isEBEEGap() {
  if (not photon_isEBEEGap_isLoaded) {
    if (photon_isEBEEGap_branch != 0) {
      photon_isEBEEGap_branch->GetEntry(index);
    } else {
      printf("branch photon_isEBEEGap_branch does not exist!\n");
      exit(1);
    }
    photon_isEBEEGap_isLoaded = true;
  }
  return photon_isEBEEGap_;
}
const bool &SkimTree::photon_isEE() {
  if (not photon_isEE_isLoaded) {
    if (photon_isEE_branch != 0) {
      photon_isEE_branch->GetEntry(index);
    } else {
      printf("branch photon_isEE_branch does not exist!\n");
      exit(1);
    }
    photon_isEE_isLoaded = true;
  }
  return photon_isEE_;
}
const bool &SkimTree::photon_isGap() {
  if (not photon_isGap_isLoaded) {
    if (photon_isGap_branch != 0) {
      photon_isGap_branch->GetEntry(index);
    } else {
      printf("branch photon_isGap_branch does not exist!\n");
      exit(1);
    }
    photon_isGap_isLoaded = true;
  }
  return photon_isGap_;
}
const bool &SkimTree::photon_is_METSafe() {
  if (not photon_is_METSafe_isLoaded) {
    if (photon_is_METSafe_branch != 0) {
      photon_is_METSafe_branch->GetEntry(index);
    } else {
      printf("branch photon_is_METSafe_branch does not exist!\n");
      exit(1);
    }
    photon_is_METSafe_isLoaded = true;
  }
  return photon_is_METSafe_;
}
const bool &SkimTree::photon_is_PFID() {
  if (not photon_is_PFID_isLoaded) {
    if (photon_is_PFID_branch != 0) {
      photon_is_PFID_branch->GetEntry(index);
    } else {
      printf("branch photon_is_PFID_branch does not exist!\n");
      exit(1);
    }
    photon_is_PFID_isLoaded = true;
  }
  return photon_is_PFID_;
}
const bool &SkimTree::photon_is_beamHaloSafe() {
  if (not photon_is_beamHaloSafe_isLoaded) {
    if (photon_is_beamHaloSafe_branch != 0) {
      photon_is_beamHaloSafe_branch->GetEntry(index);
    } else {
      printf("branch photon_is_beamHaloSafe_branch does not exist!\n");
      exit(1);
    }
    photon_is_beamHaloSafe_isLoaded = true;
  }
  return photon_is_beamHaloSafe_;
}
const bool &SkimTree::photon_is_conversionSafe() {
  if (not photon_is_conversionSafe_isLoaded) {
    if (photon_is_conversionSafe_branch != 0) {
      photon_is_conversionSafe_branch->GetEntry(index);
    } else {
      printf("branch photon_is_conversionSafe_branch does not exist!\n");
      exit(1);
    }
    photon_is_conversionSafe_isLoaded = true;
  }
  return photon_is_conversionSafe_;
}
const bool &SkimTree::photon_is_genMatched_prompt() {
  if (not photon_is_genMatched_prompt_isLoaded) {
    if (photon_is_genMatched_prompt_branch != 0) {
      photon_is_genMatched_prompt_branch->GetEntry(index);
    } else {
      printf("branch photon_is_genMatched_prompt_branch does not exist!\n");
      exit(1);
    }
    photon_is_genMatched_prompt_isLoaded = true;
  }
  return photon_is_genMatched_prompt_;
}
const bool &SkimTree::photon_is_inTime() {
  if (not photon_is_inTime_isLoaded) {
    if (photon_is_inTime_branch != 0) {
      photon_is_inTime_branch->GetEntry(index);
    } else {
      printf("branch photon_is_inTime_branch does not exist!\n");
      exit(1);
    }
    photon_is_inTime_isLoaded = true;
  }
  return photon_is_inTime_;
}
const bool &SkimTree::photon_is_spikeSafe() {
  if (not photon_is_spikeSafe_isLoaded) {
    if (photon_is_spikeSafe_branch != 0) {
      photon_is_spikeSafe_branch->GetEntry(index);
    } else {
      printf("branch photon_is_spikeSafe_branch does not exist!\n");
      exit(1);
    }
    photon_is_spikeSafe_isLoaded = true;
  }
  return photon_is_spikeSafe_;
}
const float &SkimTree::photon_mass() {
  if (not photon_mass_isLoaded) {
    if (photon_mass_branch != 0) {
      photon_mass_branch->GetEntry(index);
    } else {
      printf("branch photon_mass_branch does not exist!\n");
      exit(1);
    }
    photon_mass_isLoaded = true;
  }
  return photon_mass_;
}
const bool &SkimTree::photon_pass_HGGSelection() {
  if (not photon_pass_HGGSelection_isLoaded) {
    if (photon_pass_HGGSelection_branch != 0) {
      photon_pass_HGGSelection_branch->GetEntry(index);
    } else {
      printf("branch photon_pass_HGGSelection_branch does not exist!\n");
      exit(1);
    }
    photon_pass_HGGSelection_isLoaded = true;
  }
  return photon_pass_HGGSelection_;
}
const float &SkimTree::photon_phi() {
  if (not photon_phi_isLoaded) {
    if (photon_phi_branch != 0) {
      photon_phi_branch->GetEntry(index);
    } else {
      printf("branch photon_phi_branch does not exist!\n");
      exit(1);
    }
    photon_phi_isLoaded = true;
  }
  return photon_phi_;
}
const float &SkimTree::photon_pt() {
  if (not photon_pt_isLoaded) {
    if (photon_pt_branch != 0) {
      photon_pt_branch->GetEntry(index);
    } else {
      printf("branch photon_pt_branch does not exist!\n");
      exit(1);
    }
    photon_pt_isLoaded = true;
  }
  return photon_pt_;
}
const float &SkimTree::photon_seedTime() {
  if (not photon_seedTime_isLoaded) {
    if (photon_seedTime_branch != 0) {
      photon_seedTime_branch->GetEntry(index);
    } else {
      printf("branch photon_seedTime_branch does not exist!\n");
      exit(1);
    }
    photon_seedTime_isLoaded = true;
  }
  return photon_seedTime_;
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
  const unsigned long long &EventNumber() { return st.EventNumber(); }
  const float &KFactor_EW_NLO_qqVV_Bkg_EWDn() { return st.KFactor_EW_NLO_qqVV_Bkg_EWDn(); }
  const float &KFactor_EW_NLO_qqVV_Bkg_EWUp() { return st.KFactor_EW_NLO_qqVV_Bkg_EWUp(); }
  const float &KFactor_EW_NLO_qqVV_Bkg_Nominal() { return st.KFactor_EW_NLO_qqVV_Bkg_Nominal(); }
  const float &KFactor_EW_NLO_qqVV_Bkg_arg_mass() { return st.KFactor_EW_NLO_qqVV_Bkg_arg_mass(); }
  const float &KFactor_EW_NLO_qqVV_Bkg_arg_pthat() { return st.KFactor_EW_NLO_qqVV_Bkg_arg_pthat(); }
  const float &KFactor_EW_NLO_qqVV_Bkg_arg_rho() { return st.KFactor_EW_NLO_qqVV_Bkg_arg_rho(); }
  const float &KFactor_EW_NLO_qqVV_Bkg_arg_that() { return st.KFactor_EW_NLO_qqVV_Bkg_arg_that(); }
  const float &KFactor_QCD_NNLO_qqVV_Bkg_Nominal() { return st.KFactor_QCD_NNLO_qqVV_Bkg_Nominal(); }
  const float &KFactor_QCD_NNLO_qqVV_Bkg_arg_mass() { return st.KFactor_QCD_NNLO_qqVV_Bkg_arg_mass(); }
  const float &KFactor_QCD_NNLO_qqVV_Bkg_arg_pthat() { return st.KFactor_QCD_NNLO_qqVV_Bkg_arg_pthat(); }
  const unsigned int &LuminosityBlock() { return st.LuminosityBlock(); }
  const unsigned int &RunNumber() { return st.RunNumber(); }
  const vector<float> &ak4jets_CEMF() { return st.ak4jets_CEMF(); }
  const float &ak4jets_HT() { return st.ak4jets_HT(); }
  const float &ak4jets_MHT() { return st.ak4jets_MHT(); }
  const vector<float> &ak4jets_NEMF() { return st.ak4jets_NEMF(); }
  const vector<unsigned char> &ak4jets_btagWP_Bits() { return st.ak4jets_btagWP_Bits(); }
  const vector<float> &ak4jets_eta() { return st.ak4jets_eta(); }
  const vector<bool> &ak4jets_is_genMatched() { return st.ak4jets_is_genMatched(); }
  const vector<bool> &ak4jets_is_genMatched_fullCone() { return st.ak4jets_is_genMatched_fullCone(); }
  const vector<float> &ak4jets_mass() { return st.ak4jets_mass(); }
  const vector<float> &ak4jets_phi() { return st.ak4jets_phi(); }
  const vector<float> &ak4jets_pt() { return st.ak4jets_pt(); }
  const vector<float> &ak8jets_eta() { return st.ak8jets_eta(); }
  const vector<float> &ak8jets_mass() { return st.ak8jets_mass(); }
  const vector<float> &ak8jets_phi() { return st.ak8jets_phi(); }
  const vector<float> &ak8jets_pt() { return st.ak8jets_pt(); }
  const float &dPhi_pTboson_pTmiss() { return st.dPhi_pTboson_pTmiss(); }
  const float &dPhi_pTbosonjets_pTmiss() { return st.dPhi_pTbosonjets_pTmiss(); }
  const float &dilepton_eta() { return st.dilepton_eta(); }
  const int &dilepton_id() { return st.dilepton_id(); }
  const float &dilepton_mass() { return st.dilepton_mass(); }
  const float &dilepton_phi() { return st.dilepton_phi(); }
  const float &dilepton_pt() { return st.dilepton_pt(); }
  const float &electron_full5x5_r9() { return st.electron_full5x5_r9(); }
  const float &electron_full5x5_sigmaIEtaIEta() { return st.electron_full5x5_sigmaIEtaIEta(); }
  const float &electron_full5x5_sigmaIPhiIPhi() { return st.electron_full5x5_sigmaIPhiIPhi(); }
  const float &electron_seedTime() { return st.electron_seedTime(); }
  const vector<float> &electrons_full5x5_r9() { return st.electrons_full5x5_r9(); }
  const vector<float> &electrons_full5x5_sigmaIEtaIEta() { return st.electrons_full5x5_sigmaIEtaIEta(); }
  const vector<float> &electrons_full5x5_sigmaIPhiIPhi() { return st.electrons_full5x5_sigmaIPhiIPhi(); }
  const vector<float> &electrons_seedTime() { return st.electrons_seedTime(); }
  const float &event_mTZZ() { return st.event_mTZZ(); }
  const float &event_mZZ() { return st.event_mZZ(); }
  const float &event_mlg() { return st.event_mlg(); }
  const float &event_mllg() { return st.event_mllg(); }
  const unsigned int &event_n_ak4jets_pt20() { return st.event_n_ak4jets_pt20(); }
  const unsigned int &event_n_ak4jets_pt20_btagged_loose() { return st.event_n_ak4jets_pt20_btagged_loose(); }
  const unsigned int &event_n_ak4jets_pt20_btagged_medium() { return st.event_n_ak4jets_pt20_btagged_medium(); }
  const unsigned int &event_n_ak4jets_pt30() { return st.event_n_ak4jets_pt30(); }
  const unsigned int &event_n_ak4jets_pt30_btagged_loose() { return st.event_n_ak4jets_pt30_btagged_loose(); }
  const unsigned int &event_n_ak4jets_pt30_btagged_medium() { return st.event_n_ak4jets_pt30_btagged_medium(); }
  const unsigned int &event_n_leptons_fakeableBase() { return st.event_n_leptons_fakeableBase(); }
  const unsigned int &event_n_vtxs_good() { return st.event_n_vtxs_good(); }
  const float &event_pTmiss() { return st.event_pTmiss(); }
  const bool &event_pass_tightMETFilters() { return st.event_pass_tightMETFilters(); }
  const float &event_phimiss() { return st.event_phimiss(); }
  const float &event_wgt() { return st.event_wgt(); }
  const float &event_wgt_L1PrefiringDn() { return st.event_wgt_L1PrefiringDn(); }
  const float &event_wgt_L1PrefiringUp() { return st.event_wgt_L1PrefiringUp(); }
  const float &event_wgt_PUDn() { return st.event_wgt_PUDn(); }
  const float &event_wgt_PUUp() { return st.event_wgt_PUUp(); }
  const float &event_wgt_SFs_PUJetId() { return st.event_wgt_SFs_PUJetId(); }
  const float &event_wgt_SFs_PUJetId_EffDn() { return st.event_wgt_SFs_PUJetId_EffDn(); }
  const float &event_wgt_SFs_PUJetId_EffUp() { return st.event_wgt_SFs_PUJetId_EffUp(); }
  const float &event_wgt_SFs_btagging() { return st.event_wgt_SFs_btagging(); }
  const float &event_wgt_SFs_btagging_EffDn() { return st.event_wgt_SFs_btagging_EffDn(); }
  const float &event_wgt_SFs_btagging_EffUp() { return st.event_wgt_SFs_btagging_EffUp(); }
  const float &event_wgt_SFs_electrons() { return st.event_wgt_SFs_electrons(); }
  const float &event_wgt_SFs_electrons_AltMCDn() { return st.event_wgt_SFs_electrons_AltMCDn(); }
  const float &event_wgt_SFs_electrons_AltMCUp() { return st.event_wgt_SFs_electrons_AltMCUp(); }
  const float &event_wgt_SFs_electrons_StatDn() { return st.event_wgt_SFs_electrons_StatDn(); }
  const float &event_wgt_SFs_electrons_StatUp() { return st.event_wgt_SFs_electrons_StatUp(); }
  const float &event_wgt_SFs_electrons_SystDn() { return st.event_wgt_SFs_electrons_SystDn(); }
  const float &event_wgt_SFs_electrons_SystUp() { return st.event_wgt_SFs_electrons_SystUp(); }
  const float &event_wgt_SFs_muons() { return st.event_wgt_SFs_muons(); }
  const float &event_wgt_SFs_muons_AltMCDn() { return st.event_wgt_SFs_muons_AltMCDn(); }
  const float &event_wgt_SFs_muons_AltMCUp() { return st.event_wgt_SFs_muons_AltMCUp(); }
  const float &event_wgt_SFs_muons_StatDn() { return st.event_wgt_SFs_muons_StatDn(); }
  const float &event_wgt_SFs_muons_StatUp() { return st.event_wgt_SFs_muons_StatUp(); }
  const float &event_wgt_SFs_muons_SystDn() { return st.event_wgt_SFs_muons_SystDn(); }
  const float &event_wgt_SFs_muons_SystUp() { return st.event_wgt_SFs_muons_SystUp(); }
  const float &event_wgt_SFs_photons() { return st.event_wgt_SFs_photons(); }
  const float &event_wgt_SFs_photons_EffDn() { return st.event_wgt_SFs_photons_EffDn(); }
  const float &event_wgt_SFs_photons_EffUp() { return st.event_wgt_SFs_photons_EffUp(); }
  const float &event_wgt_adjustment_AsMZDn() { return st.event_wgt_adjustment_AsMZDn(); }
  const float &event_wgt_adjustment_AsMZUp() { return st.event_wgt_adjustment_AsMZUp(); }
  const float &event_wgt_adjustment_NNPDF30() { return st.event_wgt_adjustment_NNPDF30(); }
  const float &event_wgt_adjustment_NNPDF30_AsMZDn() { return st.event_wgt_adjustment_NNPDF30_AsMZDn(); }
  const float &event_wgt_adjustment_NNPDF30_AsMZUp() { return st.event_wgt_adjustment_NNPDF30_AsMZUp(); }
  const float &event_wgt_adjustment_NNPDF30_PDFReplicaDn() { return st.event_wgt_adjustment_NNPDF30_PDFReplicaDn(); }
  const float &event_wgt_adjustment_NNPDF30_PDFReplicaUp() { return st.event_wgt_adjustment_NNPDF30_PDFReplicaUp(); }
  const float &event_wgt_adjustment_PDFReplicaDn() { return st.event_wgt_adjustment_PDFReplicaDn(); }
  const float &event_wgt_adjustment_PDFReplicaUp() { return st.event_wgt_adjustment_PDFReplicaUp(); }
  const float &event_wgt_adjustment_PDFScaleDn() { return st.event_wgt_adjustment_PDFScaleDn(); }
  const float &event_wgt_adjustment_PDFScaleUp() { return st.event_wgt_adjustment_PDFScaleUp(); }
  const float &event_wgt_adjustment_PythiaScaleDn() { return st.event_wgt_adjustment_PythiaScaleDn(); }
  const float &event_wgt_adjustment_PythiaScaleUp() { return st.event_wgt_adjustment_PythiaScaleUp(); }
  const float &event_wgt_adjustment_QCDScaleDn() { return st.event_wgt_adjustment_QCDScaleDn(); }
  const float &event_wgt_adjustment_QCDScaleUp() { return st.event_wgt_adjustment_QCDScaleUp(); }
  const float &event_wgt_triggers() { return st.event_wgt_triggers(); }
  const float &event_wgt_triggers_Dilepton() { return st.event_wgt_triggers_Dilepton(); }
  const float &event_wgt_triggers_Dilepton_DF_Extra() { return st.event_wgt_triggers_Dilepton_DF_Extra(); }
  const float &event_wgt_triggers_PFHT_Control() { return st.event_wgt_triggers_PFHT_Control(); }
  const float &event_wgt_triggers_PFMET_MHT_Control() { return st.event_wgt_triggers_PFMET_MHT_Control(); }
  const float &event_wgt_triggers_SingleLepton() { return st.event_wgt_triggers_SingleLepton(); }
  const float &event_wgt_triggers_SinglePhoton() { return st.event_wgt_triggers_SinglePhoton(); }
  const vector<float> &genak4jets_eta() { return st.genak4jets_eta(); }
  const vector<float> &genak4jets_mass() { return st.genak4jets_mass(); }
  const vector<float> &genak4jets_phi() { return st.genak4jets_phi(); }
  const vector<float> &genak4jets_pt() { return st.genak4jets_pt(); }
  const float &genmet_pTmiss() { return st.genmet_pTmiss(); }
  const float &genmet_phimiss() { return st.genmet_phimiss(); }
  const float &lepton_eff() { return st.lepton_eff(); }
  const float &lepton_eff_StatDn() { return st.lepton_eff_StatDn(); }
  const float &lepton_eff_StatUp() { return st.lepton_eff_StatUp(); }
  const float &lepton_eff_SystDn() { return st.lepton_eff_SystDn(); }
  const float &lepton_eff_SystUp() { return st.lepton_eff_SystUp(); }
  const float &lepton_eta() { return st.lepton_eta(); }
  const int &lepton_id() { return st.lepton_id(); }
  const bool &lepton_is_genMatched_prompt() { return st.lepton_is_genMatched_prompt(); }
  const float &lepton_mass() { return st.lepton_mass(); }
  const float &lepton_phi() { return st.lepton_phi(); }
  const float &lepton_pt() { return st.lepton_pt(); }
  const float &genpromptparticles_sump4_mass() { return st.genpromptparticles_sump4_mass(); }
  const float &genpromptparticles_sump4_pt() { return st.genpromptparticles_sump4_pt(); }
  const float &genpromptparticles_sump4_rapidity() { return st.genpromptparticles_sump4_rapidity(); }
  const vector<float> &leptons_eff() { return st.leptons_eff(); }
  const vector<float> &leptons_eff_DF() { return st.leptons_eff_DF(); }
  const vector<float> &leptons_eff_DF_StatDn() { return st.leptons_eff_DF_StatDn(); }
  const vector<float> &leptons_eff_DF_StatUp() { return st.leptons_eff_DF_StatUp(); }
  const vector<float> &leptons_eff_DF_SystDn() { return st.leptons_eff_DF_SystDn(); }
  const vector<float> &leptons_eff_DF_SystUp() { return st.leptons_eff_DF_SystUp(); }
  const vector<float> &leptons_eff_StatDn() { return st.leptons_eff_StatDn(); }
  const vector<float> &leptons_eff_StatUp() { return st.leptons_eff_StatUp(); }
  const vector<float> &leptons_eff_SystDn() { return st.leptons_eff_SystDn(); }
  const vector<float> &leptons_eff_SystUp() { return st.leptons_eff_SystUp(); }
  const vector<float> &leptons_eta() { return st.leptons_eta(); }
  const vector<int> &leptons_id() { return st.leptons_id(); }
  const vector<bool> &leptons_is_TOmatched_SingleLepton() { return st.leptons_is_TOmatched_SingleLepton(); }
  const vector<bool> &leptons_is_genMatched_prompt() { return st.leptons_is_genMatched_prompt(); }
  const vector<float> &leptons_mass() { return st.leptons_mass(); }
  const vector<float> &leptons_phi() { return st.leptons_phi(); }
  const vector<float> &leptons_pt() { return st.leptons_pt(); }
  const float &min_abs_dPhi_pTj_pTmiss() { return st.min_abs_dPhi_pTj_pTmiss(); }
  const float &pConst_JJQCD_SIG_ghg2_1_JHUGen() { return st.pConst_JJQCD_SIG_ghg2_1_JHUGen(); }
  const float &pConst_JJVBF_SIG_ghv1_1_JHUGen() { return st.pConst_JJVBF_SIG_ghv1_1_JHUGen(); }
  const float &p_JJQCD_SIG_ghg2_1_JHUGen() { return st.p_JJQCD_SIG_ghg2_1_JHUGen(); }
  const float &p_JJVBF_SIG_ghv1_1_JHUGen() { return st.p_JJVBF_SIG_ghv1_1_JHUGen(); }
  const float &p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen() { return st.p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen(); }
  const float &p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen() { return st.p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen(); }
  const float &p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen() { return st.p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen(); }
  const float &p_JJVBF_SIG_ghv1prime2_1E4_JHUGen() { return st.p_JJVBF_SIG_ghv1prime2_1E4_JHUGen(); }
  const float &p_JJVBF_SIG_ghv2_1_JHUGen() { return st.p_JJVBF_SIG_ghv2_1_JHUGen(); }
  const float &p_JJVBF_SIG_ghv4_1_JHUGen() { return st.p_JJVBF_SIG_ghv4_1_JHUGen(); }
  const float &p_JJVBF_SIG_ghza1prime2_1E4_JHUGen() { return st.p_JJVBF_SIG_ghza1prime2_1E4_JHUGen(); }
  const float &photon_MIPTotalEnergy() { return st.photon_MIPTotalEnergy(); }
  const float &photon_eta() { return st.photon_eta(); }
  const float &photon_full5x5_r9() { return st.photon_full5x5_r9(); }
  const float &photon_full5x5_sigmaIEtaIEta() { return st.photon_full5x5_sigmaIEtaIEta(); }
  const float &photon_full5x5_sigmaIPhiIPhi() { return st.photon_full5x5_sigmaIPhiIPhi(); }
  const bool &photon_isEB() { return st.photon_isEB(); }
  const bool &photon_isEBEEGap() { return st.photon_isEBEEGap(); }
  const bool &photon_isEE() { return st.photon_isEE(); }
  const bool &photon_isGap() { return st.photon_isGap(); }
  const bool &photon_is_METSafe() { return st.photon_is_METSafe(); }
  const bool &photon_is_PFID() { return st.photon_is_PFID(); }
  const bool &photon_is_beamHaloSafe() { return st.photon_is_beamHaloSafe(); }
  const bool &photon_is_conversionSafe() { return st.photon_is_conversionSafe(); }
  const bool &photon_is_genMatched_prompt() { return st.photon_is_genMatched_prompt(); }
  const bool &photon_is_inTime() { return st.photon_is_inTime(); }
  const bool &photon_is_spikeSafe() { return st.photon_is_spikeSafe(); }
  const float &photon_mass() { return st.photon_mass(); }
  const bool &photon_pass_HGGSelection() { return st.photon_pass_HGGSelection(); }
  const float &photon_phi() { return st.photon_phi(); }
  const float &photon_pt() { return st.photon_pt(); }
  const float &photon_seedTime() { return st.photon_seedTime(); }
}
