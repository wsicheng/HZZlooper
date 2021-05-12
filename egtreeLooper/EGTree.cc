#include "EGTree.h"
EGTree st;

void EGTree::Init(TTree *tree) {
  tree->SetMakeClass(1);
  dR_e_g_branch = 0;
  if (tree->GetBranch("dR_e_g") != 0) {
    dR_e_g_branch = tree->GetBranch("dR_e_g");
    if (dR_e_g_branch) { dR_e_g_branch->SetAddress(&dR_e_g_); }
  }
  electron_dR_genMatch_branch = 0;
  if (tree->GetBranch("electron_dR_genMatch") != 0) {
    electron_dR_genMatch_branch = tree->GetBranch("electron_dR_genMatch");
    if (electron_dR_genMatch_branch) { electron_dR_genMatch_branch->SetAddress(&electron_dR_genMatch_); }
  }
  electron_dxy_branch = 0;
  if (tree->GetBranch("electron_dxy") != 0) {
    electron_dxy_branch = tree->GetBranch("electron_dxy");
    if (electron_dxy_branch) { electron_dxy_branch->SetAddress(&electron_dxy_); }
  }
  electron_dz_branch = 0;
  if (tree->GetBranch("electron_dz") != 0) {
    electron_dz_branch = tree->GetBranch("electron_dz");
    if (electron_dz_branch) { electron_dz_branch->SetAddress(&electron_dz_); }
  }
  electron_eta_branch = 0;
  if (tree->GetBranch("electron_eta") != 0) {
    electron_eta_branch = tree->GetBranch("electron_eta");
    if (electron_eta_branch) { electron_eta_branch->SetAddress(&electron_eta_); }
  }
  electron_etaSC_branch = 0;
  if (tree->GetBranch("electron_etaSC") != 0) {
    electron_etaSC_branch = tree->GetBranch("electron_etaSC");
    if (electron_etaSC_branch) { electron_etaSC_branch->SetAddress(&electron_etaSC_); }
  }
  electron_fid_mask_branch = 0;
  if (tree->GetBranch("electron_fid_mask") != 0) {
    electron_fid_mask_branch = tree->GetBranch("electron_fid_mask");
    if (electron_fid_mask_branch) { electron_fid_mask_branch->SetAddress(&electron_fid_mask_); }
  }
  electron_hasTightCharge_branch = 0;
  if (tree->GetBranch("electron_hasTightCharge") != 0) {
    electron_hasTightCharge_branch = tree->GetBranch("electron_hasTightCharge");
    if (electron_hasTightCharge_branch) { electron_hasTightCharge_branch->SetAddress(&electron_hasTightCharge_); }
  }
  electron_id_branch = 0;
  if (tree->GetBranch("electron_id") != 0) {
    electron_id_branch = tree->GetBranch("electron_id");
    if (electron_id_branch) { electron_id_branch->SetAddress(&electron_id_); }
  }
  electron_id_genMatch_branch = 0;
  if (tree->GetBranch("electron_id_genMatch") != 0) {
    electron_id_genMatch_branch = tree->GetBranch("electron_id_genMatch");
    if (electron_id_genMatch_branch) { electron_id_genMatch_branch->SetAddress(&electron_id_genMatch_); }
  }
  electron_is_conversionSafe_branch = 0;
  if (tree->GetBranch("electron_is_conversionSafe") != 0) {
    electron_is_conversionSafe_branch = tree->GetBranch("electron_is_conversionSafe");
    if (electron_is_conversionSafe_branch) { electron_is_conversionSafe_branch->SetAddress(&electron_is_conversionSafe_); }
  }
  electron_is_extraTight_branch = 0;
  if (tree->GetBranch("electron_is_extraTight") != 0) {
    electron_is_extraTight_branch = tree->GetBranch("electron_is_extraTight");
    if (electron_is_extraTight_branch) { electron_is_extraTight_branch->SetAddress(&electron_is_extraTight_); }
  }
  electron_is_genMatched_prompt_branch = 0;
  if (tree->GetBranch("electron_is_genMatched_prompt") != 0) {
    electron_is_genMatched_prompt_branch = tree->GetBranch("electron_is_genMatched_prompt");
    if (electron_is_genMatched_prompt_branch) { electron_is_genMatched_prompt_branch->SetAddress(&electron_is_genMatched_prompt_); }
  }
  electron_minDR_electron_branch = 0;
  if (tree->GetBranch("electron_minDR_electron") != 0) {
    electron_minDR_electron_branch = tree->GetBranch("electron_minDR_electron");
    if (electron_minDR_electron_branch) { electron_minDR_electron_branch->SetAddress(&electron_minDR_electron_); }
  }
  electron_minDR_muon_branch = 0;
  if (tree->GetBranch("electron_minDR_muon") != 0) {
    electron_minDR_muon_branch = tree->GetBranch("electron_minDR_muon");
    if (electron_minDR_muon_branch) { electron_minDR_muon_branch->SetAddress(&electron_minDR_muon_); }
  }
  electron_minDR_photon_branch = 0;
  if (tree->GetBranch("electron_minDR_photon") != 0) {
    electron_minDR_photon_branch = tree->GetBranch("electron_minDR_photon");
    if (electron_minDR_photon_branch) { electron_minDR_photon_branch->SetAddress(&electron_minDR_photon_); }
  }
  electron_phi_branch = 0;
  if (tree->GetBranch("electron_phi") != 0) {
    electron_phi_branch = tree->GetBranch("electron_phi");
    if (electron_phi_branch) { electron_phi_branch->SetAddress(&electron_phi_); }
  }
  electron_pt_branch = 0;
  if (tree->GetBranch("electron_pt") != 0) {
    electron_pt_branch = tree->GetBranch("electron_pt");
    if (electron_pt_branch) { electron_pt_branch->SetAddress(&electron_pt_); }
  }
  eta_eg_branch = 0;
  if (tree->GetBranch("eta_eg") != 0) {
    eta_eg_branch = tree->GetBranch("eta_eg");
    if (eta_eg_branch) { eta_eg_branch->SetAddress(&eta_eg_); }
  }
  event_NGenPromptParticles_branch = 0;
  if (tree->GetBranch("event_NGenPromptParticles") != 0) {
    event_NGenPromptParticles_branch = tree->GetBranch("event_NGenPromptParticles");
    if (event_NGenPromptParticles_branch) { event_NGenPromptParticles_branch->SetAddress(&event_NGenPromptParticles_); }
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
  event_n_leptons_fakeableBase_branch = 0;
  if (tree->GetBranch("event_n_leptons_fakeableBase") != 0) {
    event_n_leptons_fakeableBase_branch = tree->GetBranch("event_n_leptons_fakeableBase");
    if (event_n_leptons_fakeableBase_branch) { event_n_leptons_fakeableBase_branch->SetAddress(&event_n_leptons_fakeableBase_); }
  }
  event_nvtxs_good_branch = 0;
  if (tree->GetBranch("event_nvtxs_good") != 0) {
    event_nvtxs_good_branch = tree->GetBranch("event_nvtxs_good");
    if (event_nvtxs_good_branch) { event_nvtxs_good_branch->SetAddress(&event_nvtxs_good_); }
  }
  event_pTmiss_branch = 0;
  if (tree->GetBranch("event_pTmiss") != 0) {
    event_pTmiss_branch = tree->GetBranch("event_pTmiss");
    if (event_pTmiss_branch) { event_pTmiss_branch->SetAddress(&event_pTmiss_); }
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
  event_wgt_SFs_branch = 0;
  if (tree->GetBranch("event_wgt_SFs") != 0) {
    event_wgt_SFs_branch = tree->GetBranch("event_wgt_SFs");
    if (event_wgt_SFs_branch) { event_wgt_SFs_branch->SetAddress(&event_wgt_SFs_); }
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
  mass_eg_branch = 0;
  if (tree->GetBranch("mass_eg") != 0) {
    mass_eg_branch = tree->GetBranch("mass_eg");
    if (mass_eg_branch) { mass_eg_branch->SetAddress(&mass_eg_); }
  }
  phi_eg_branch = 0;
  if (tree->GetBranch("phi_eg") != 0) {
    phi_eg_branch = tree->GetBranch("phi_eg");
    if (phi_eg_branch) { phi_eg_branch->SetAddress(&phi_eg_); }
  }
  photon_MIPTotalEnergy_branch = 0;
  if (tree->GetBranch("photon_MIPTotalEnergy") != 0) {
    photon_MIPTotalEnergy_branch = tree->GetBranch("photon_MIPTotalEnergy");
    if (photon_MIPTotalEnergy_branch) { photon_MIPTotalEnergy_branch->SetAddress(&photon_MIPTotalEnergy_); }
  }
  photon_dR_genMatch_branch = 0;
  if (tree->GetBranch("photon_dR_genMatch") != 0) {
    photon_dR_genMatch_branch = tree->GetBranch("photon_dR_genMatch");
    if (photon_dR_genMatch_branch) { photon_dR_genMatch_branch->SetAddress(&photon_dR_genMatch_); }
  }
  photon_eta_branch = 0;
  if (tree->GetBranch("photon_eta") != 0) {
    photon_eta_branch = tree->GetBranch("photon_eta");
    if (photon_eta_branch) { photon_eta_branch->SetAddress(&photon_eta_); }
  }
  photon_etaSC_branch = 0;
  if (tree->GetBranch("photon_etaSC") != 0) {
    photon_etaSC_branch = tree->GetBranch("photon_etaSC");
    if (photon_etaSC_branch) { photon_etaSC_branch->SetAddress(&photon_etaSC_); }
  }
  photon_fid_mask_branch = 0;
  if (tree->GetBranch("photon_fid_mask") != 0) {
    photon_fid_mask_branch = tree->GetBranch("photon_fid_mask");
    if (photon_fid_mask_branch) { photon_fid_mask_branch->SetAddress(&photon_fid_mask_); }
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
  photon_id_branch = 0;
  if (tree->GetBranch("photon_id") != 0) {
    photon_id_branch = tree->GetBranch("photon_id");
    if (photon_id_branch) { photon_id_branch->SetAddress(&photon_id_); }
  }
  photon_id_genMatch_branch = 0;
  if (tree->GetBranch("photon_id_genMatch") != 0) {
    photon_id_genMatch_branch = tree->GetBranch("photon_id_genMatch");
    if (photon_id_genMatch_branch) { photon_id_genMatch_branch->SetAddress(&photon_id_genMatch_); }
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
  photon_minDR_electron_branch = 0;
  if (tree->GetBranch("photon_minDR_electron") != 0) {
    photon_minDR_electron_branch = tree->GetBranch("photon_minDR_electron");
    if (photon_minDR_electron_branch) { photon_minDR_electron_branch->SetAddress(&photon_minDR_electron_); }
  }
  photon_minDR_muon_branch = 0;
  if (tree->GetBranch("photon_minDR_muon") != 0) {
    photon_minDR_muon_branch = tree->GetBranch("photon_minDR_muon");
    if (photon_minDR_muon_branch) { photon_minDR_muon_branch->SetAddress(&photon_minDR_muon_); }
  }
  photon_minDR_photon_branch = 0;
  if (tree->GetBranch("photon_minDR_photon") != 0) {
    photon_minDR_photon_branch = tree->GetBranch("photon_minDR_photon");
    if (photon_minDR_photon_branch) { photon_minDR_photon_branch->SetAddress(&photon_minDR_photon_); }
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
  pt_eg_branch = 0;
  if (tree->GetBranch("pt_eg") != 0) {
    pt_eg_branch = tree->GetBranch("pt_eg");
    if (pt_eg_branch) { pt_eg_branch->SetAddress(&pt_eg_); }
  }
  weight_HLT_Photon120_R9Id90_HE10_IsoM_branch = 0;
  if (tree->GetBranch("weight_HLT_Photon120_R9Id90_HE10_IsoM") != 0) {
    weight_HLT_Photon120_R9Id90_HE10_IsoM_branch = tree->GetBranch("weight_HLT_Photon120_R9Id90_HE10_IsoM");
    if (weight_HLT_Photon120_R9Id90_HE10_IsoM_branch) { weight_HLT_Photon120_R9Id90_HE10_IsoM_branch->SetAddress(&weight_HLT_Photon120_R9Id90_HE10_IsoM_); }
  }
  weight_HLT_Photon165_R9Id90_HE10_IsoM_branch = 0;
  if (tree->GetBranch("weight_HLT_Photon165_R9Id90_HE10_IsoM") != 0) {
    weight_HLT_Photon165_R9Id90_HE10_IsoM_branch = tree->GetBranch("weight_HLT_Photon165_R9Id90_HE10_IsoM");
    if (weight_HLT_Photon165_R9Id90_HE10_IsoM_branch) { weight_HLT_Photon165_R9Id90_HE10_IsoM_branch->SetAddress(&weight_HLT_Photon165_R9Id90_HE10_IsoM_); }
  }
  weight_HLT_Photon200_branch = 0;
  if (tree->GetBranch("weight_HLT_Photon200") != 0) {
    weight_HLT_Photon200_branch = tree->GetBranch("weight_HLT_Photon200");
    if (weight_HLT_Photon200_branch) { weight_HLT_Photon200_branch->SetAddress(&weight_HLT_Photon200_); }
  }
  weight_HLT_Photon175_branch = 0;
  if (tree->GetBranch("weight_HLT_Photon175") != 0) {
    weight_HLT_Photon175_branch = tree->GetBranch("weight_HLT_Photon175");
    if (weight_HLT_Photon175_branch) { weight_HLT_Photon175_branch->SetAddress(&weight_HLT_Photon175_); }
  }
  weight_HLT_Photon20_HoverELoose_branch = 0;
  if (tree->GetBranch("weight_HLT_Photon20_HoverELoose") != 0) {
    weight_HLT_Photon20_HoverELoose_branch = tree->GetBranch("weight_HLT_Photon20_HoverELoose");
    if (weight_HLT_Photon20_HoverELoose_branch) { weight_HLT_Photon20_HoverELoose_branch->SetAddress(&weight_HLT_Photon20_HoverELoose_); }
  }
  weight_HLT_Photon30_HoverELoose_branch = 0;
  if (tree->GetBranch("weight_HLT_Photon30_HoverELoose") != 0) {
    weight_HLT_Photon30_HoverELoose_branch = tree->GetBranch("weight_HLT_Photon30_HoverELoose");
    if (weight_HLT_Photon30_HoverELoose_branch) { weight_HLT_Photon30_HoverELoose_branch->SetAddress(&weight_HLT_Photon30_HoverELoose_); }
  }
  weight_HLT_Photon33_branch = 0;
  if (tree->GetBranch("weight_HLT_Photon33") != 0) {
    weight_HLT_Photon33_branch = tree->GetBranch("weight_HLT_Photon33");
    if (weight_HLT_Photon33_branch) { weight_HLT_Photon33_branch->SetAddress(&weight_HLT_Photon33_); }
  }
  weight_HLT_Photon50_R9Id90_HE10_IsoM_branch = 0;
  if (tree->GetBranch("weight_HLT_Photon50_R9Id90_HE10_IsoM") != 0) {
    weight_HLT_Photon50_R9Id90_HE10_IsoM_branch = tree->GetBranch("weight_HLT_Photon50_R9Id90_HE10_IsoM");
    if (weight_HLT_Photon50_R9Id90_HE10_IsoM_branch) { weight_HLT_Photon50_R9Id90_HE10_IsoM_branch->SetAddress(&weight_HLT_Photon50_R9Id90_HE10_IsoM_); }
  }
  weight_HLT_Photon75_R9Id90_HE10_IsoM_branch = 0;
  if (tree->GetBranch("weight_HLT_Photon75_R9Id90_HE10_IsoM") != 0) {
    weight_HLT_Photon75_R9Id90_HE10_IsoM_branch = tree->GetBranch("weight_HLT_Photon75_R9Id90_HE10_IsoM");
    if (weight_HLT_Photon75_R9Id90_HE10_IsoM_branch) { weight_HLT_Photon75_R9Id90_HE10_IsoM_branch->SetAddress(&weight_HLT_Photon75_R9Id90_HE10_IsoM_); }
  }
  weight_HLT_Photon90_R9Id90_HE10_IsoM_branch = 0;
  if (tree->GetBranch("weight_HLT_Photon90_R9Id90_HE10_IsoM") != 0) {
    weight_HLT_Photon90_R9Id90_HE10_IsoM_branch = tree->GetBranch("weight_HLT_Photon90_R9Id90_HE10_IsoM");
    if (weight_HLT_Photon90_R9Id90_HE10_IsoM_branch) { weight_HLT_Photon90_R9Id90_HE10_IsoM_branch->SetAddress(&weight_HLT_Photon90_R9Id90_HE10_IsoM_); }
  }
  tree->SetMakeClass(0);
}
void EGTree::GetEntry(unsigned int idx) {
  index = idx;
  dR_e_g_isLoaded = false;
  electron_dR_genMatch_isLoaded = false;
  electron_dxy_isLoaded = false;
  electron_dz_isLoaded = false;
  electron_eta_isLoaded = false;
  electron_etaSC_isLoaded = false;
  electron_fid_mask_isLoaded = false;
  electron_hasTightCharge_isLoaded = false;
  electron_id_isLoaded = false;
  electron_id_genMatch_isLoaded = false;
  electron_is_conversionSafe_isLoaded = false;
  electron_is_extraTight_isLoaded = false;
  electron_is_genMatched_prompt_isLoaded = false;
  electron_minDR_electron_isLoaded = false;
  electron_minDR_muon_isLoaded = false;
  electron_minDR_photon_isLoaded = false;
  electron_phi_isLoaded = false;
  electron_pt_isLoaded = false;
  eta_eg_isLoaded = false;
  event_NGenPromptParticles_isLoaded = false;
  event_n_ak4jets_pt20_isLoaded = false;
  event_n_ak4jets_pt20_btagged_loose_isLoaded = false;
  event_n_ak4jets_pt30_isLoaded = false;
  event_n_ak4jets_pt30_btagged_loose_isLoaded = false;
  event_n_leptons_fakeableBase_isLoaded = false;
  event_nvtxs_good_isLoaded = false;
  event_pTmiss_isLoaded = false;
  event_phimiss_isLoaded = false;
  event_wgt_isLoaded = false;
  event_wgt_SFs_isLoaded = false;
  genmet_pTmiss_isLoaded = false;
  genmet_phimiss_isLoaded = false;
  mass_eg_isLoaded = false;
  phi_eg_isLoaded = false;
  photon_MIPTotalEnergy_isLoaded = false;
  photon_dR_genMatch_isLoaded = false;
  photon_eta_isLoaded = false;
  photon_etaSC_isLoaded = false;
  photon_fid_mask_isLoaded = false;
  photon_full5x5_r9_isLoaded = false;
  photon_full5x5_sigmaIEtaIEta_isLoaded = false;
  photon_full5x5_sigmaIPhiIPhi_isLoaded = false;
  photon_id_isLoaded = false;
  photon_id_genMatch_isLoaded = false;
  photon_is_METSafe_isLoaded = false;
  photon_is_PFID_isLoaded = false;
  photon_is_beamHaloSafe_isLoaded = false;
  photon_is_conversionSafe_isLoaded = false;
  photon_is_genMatched_prompt_isLoaded = false;
  photon_is_inTime_isLoaded = false;
  photon_is_spikeSafe_isLoaded = false;
  photon_minDR_electron_isLoaded = false;
  photon_minDR_muon_isLoaded = false;
  photon_minDR_photon_isLoaded = false;
  photon_phi_isLoaded = false;
  photon_pt_isLoaded = false;
  photon_seedTime_isLoaded = false;
  pt_eg_isLoaded = false;
  weight_HLT_Photon120_R9Id90_HE10_IsoM_isLoaded = false;
  weight_HLT_Photon165_R9Id90_HE10_IsoM_isLoaded = false;
  weight_HLT_Photon200_isLoaded = false;
  weight_HLT_Photon175_isLoaded = false;
  weight_HLT_Photon20_HoverELoose_isLoaded = false;
  weight_HLT_Photon30_HoverELoose_isLoaded = false;
  weight_HLT_Photon33_isLoaded = false;
  weight_HLT_Photon50_R9Id90_HE10_IsoM_isLoaded = false;
  weight_HLT_Photon75_R9Id90_HE10_IsoM_isLoaded = false;
  weight_HLT_Photon90_R9Id90_HE10_IsoM_isLoaded = false;
}
void EGTree::LoadAllBranches() {
  if (dR_e_g_branch != 0) dR_e_g();
  if (electron_dR_genMatch_branch != 0) electron_dR_genMatch();
  if (electron_dxy_branch != 0) electron_dxy();
  if (electron_dz_branch != 0) electron_dz();
  if (electron_eta_branch != 0) electron_eta();
  if (electron_etaSC_branch != 0) electron_etaSC();
  if (electron_fid_mask_branch != 0) electron_fid_mask();
  if (electron_hasTightCharge_branch != 0) electron_hasTightCharge();
  if (electron_id_branch != 0) electron_id();
  if (electron_id_genMatch_branch != 0) electron_id_genMatch();
  if (electron_is_conversionSafe_branch != 0) electron_is_conversionSafe();
  if (electron_is_extraTight_branch != 0) electron_is_extraTight();
  if (electron_is_genMatched_prompt_branch != 0) electron_is_genMatched_prompt();
  if (electron_minDR_electron_branch != 0) electron_minDR_electron();
  if (electron_minDR_muon_branch != 0) electron_minDR_muon();
  if (electron_minDR_photon_branch != 0) electron_minDR_photon();
  if (electron_phi_branch != 0) electron_phi();
  if (electron_pt_branch != 0) electron_pt();
  if (eta_eg_branch != 0) eta_eg();
  if (event_NGenPromptParticles_branch != 0) event_NGenPromptParticles();
  if (event_n_ak4jets_pt20_branch != 0) event_n_ak4jets_pt20();
  if (event_n_ak4jets_pt20_btagged_loose_branch != 0) event_n_ak4jets_pt20_btagged_loose();
  if (event_n_ak4jets_pt30_branch != 0) event_n_ak4jets_pt30();
  if (event_n_ak4jets_pt30_btagged_loose_branch != 0) event_n_ak4jets_pt30_btagged_loose();
  if (event_n_leptons_fakeableBase_branch != 0) event_n_leptons_fakeableBase();
  if (event_nvtxs_good_branch != 0) event_nvtxs_good();
  if (event_pTmiss_branch != 0) event_pTmiss();
  if (event_phimiss_branch != 0) event_phimiss();
  if (event_wgt_branch != 0) event_wgt();
  if (event_wgt_SFs_branch != 0) event_wgt_SFs();
  if (genmet_pTmiss_branch != 0) genmet_pTmiss();
  if (genmet_phimiss_branch != 0) genmet_phimiss();
  if (mass_eg_branch != 0) mass_eg();
  if (phi_eg_branch != 0) phi_eg();
  if (photon_MIPTotalEnergy_branch != 0) photon_MIPTotalEnergy();
  if (photon_dR_genMatch_branch != 0) photon_dR_genMatch();
  if (photon_eta_branch != 0) photon_eta();
  if (photon_etaSC_branch != 0) photon_etaSC();
  if (photon_fid_mask_branch != 0) photon_fid_mask();
  if (photon_full5x5_r9_branch != 0) photon_full5x5_r9();
  if (photon_full5x5_sigmaIEtaIEta_branch != 0) photon_full5x5_sigmaIEtaIEta();
  if (photon_full5x5_sigmaIPhiIPhi_branch != 0) photon_full5x5_sigmaIPhiIPhi();
  if (photon_id_branch != 0) photon_id();
  if (photon_id_genMatch_branch != 0) photon_id_genMatch();
  if (photon_is_METSafe_branch != 0) photon_is_METSafe();
  if (photon_is_PFID_branch != 0) photon_is_PFID();
  if (photon_is_beamHaloSafe_branch != 0) photon_is_beamHaloSafe();
  if (photon_is_conversionSafe_branch != 0) photon_is_conversionSafe();
  if (photon_is_genMatched_prompt_branch != 0) photon_is_genMatched_prompt();
  if (photon_is_inTime_branch != 0) photon_is_inTime();
  if (photon_is_spikeSafe_branch != 0) photon_is_spikeSafe();
  if (photon_minDR_electron_branch != 0) photon_minDR_electron();
  if (photon_minDR_muon_branch != 0) photon_minDR_muon();
  if (photon_minDR_photon_branch != 0) photon_minDR_photon();
  if (photon_phi_branch != 0) photon_phi();
  if (photon_pt_branch != 0) photon_pt();
  if (photon_seedTime_branch != 0) photon_seedTime();
  if (pt_eg_branch != 0) pt_eg();
  if (weight_HLT_Photon120_R9Id90_HE10_IsoM_branch != 0) weight_HLT_Photon120_R9Id90_HE10_IsoM();
  if (weight_HLT_Photon165_R9Id90_HE10_IsoM_branch != 0) weight_HLT_Photon165_R9Id90_HE10_IsoM();
  if (weight_HLT_Photon200_branch != 0) weight_HLT_Photon200();
  if (weight_HLT_Photon175_branch != 0) weight_HLT_Photon175();
  if (weight_HLT_Photon20_HoverELoose_branch != 0) weight_HLT_Photon20_HoverELoose();
  if (weight_HLT_Photon30_HoverELoose_branch != 0) weight_HLT_Photon30_HoverELoose();
  if (weight_HLT_Photon33_branch != 0) weight_HLT_Photon33();
  if (weight_HLT_Photon50_R9Id90_HE10_IsoM_branch != 0) weight_HLT_Photon50_R9Id90_HE10_IsoM();
  if (weight_HLT_Photon75_R9Id90_HE10_IsoM_branch != 0) weight_HLT_Photon75_R9Id90_HE10_IsoM();
  if (weight_HLT_Photon90_R9Id90_HE10_IsoM_branch != 0) weight_HLT_Photon90_R9Id90_HE10_IsoM();
}
const vector<float> &EGTree::dR_e_g() {
  if (not dR_e_g_isLoaded) {
    if (dR_e_g_branch != 0) {
      dR_e_g_branch->GetEntry(index);
    } else {
      printf("branch dR_e_g_branch does not exist!\n");
      exit(1);
    }
    dR_e_g_isLoaded = true;
  }
  return *dR_e_g_;
}
const vector<float> &EGTree::electron_dR_genMatch() {
  if (not electron_dR_genMatch_isLoaded) {
    if (electron_dR_genMatch_branch != 0) {
      electron_dR_genMatch_branch->GetEntry(index);
    } else {
      printf("branch electron_dR_genMatch_branch does not exist!\n");
      exit(1);
    }
    electron_dR_genMatch_isLoaded = true;
  }
  return *electron_dR_genMatch_;
}
const vector<float> &EGTree::electron_dxy() {
  if (not electron_dxy_isLoaded) {
    if (electron_dxy_branch != 0) {
      electron_dxy_branch->GetEntry(index);
    } else {
      printf("branch electron_dxy_branch does not exist!\n");
      exit(1);
    }
    electron_dxy_isLoaded = true;
  }
  return *electron_dxy_;
}
const vector<float> &EGTree::electron_dz() {
  if (not electron_dz_isLoaded) {
    if (electron_dz_branch != 0) {
      electron_dz_branch->GetEntry(index);
    } else {
      printf("branch electron_dz_branch does not exist!\n");
      exit(1);
    }
    electron_dz_isLoaded = true;
  }
  return *electron_dz_;
}
const vector<float> &EGTree::electron_eta() {
  if (not electron_eta_isLoaded) {
    if (electron_eta_branch != 0) {
      electron_eta_branch->GetEntry(index);
    } else {
      printf("branch electron_eta_branch does not exist!\n");
      exit(1);
    }
    electron_eta_isLoaded = true;
  }
  return *electron_eta_;
}
const vector<float> &EGTree::electron_etaSC() {
  if (not electron_etaSC_isLoaded) {
    if (electron_etaSC_branch != 0) {
      electron_etaSC_branch->GetEntry(index);
    } else {
      printf("branch electron_etaSC_branch does not exist!\n");
      exit(1);
    }
    electron_etaSC_isLoaded = true;
  }
  return *electron_etaSC_;
}
const vector<unsigned short> &EGTree::electron_fid_mask() {
  if (not electron_fid_mask_isLoaded) {
    if (electron_fid_mask_branch != 0) {
      electron_fid_mask_branch->GetEntry(index);
    } else {
      printf("branch electron_fid_mask_branch does not exist!\n");
      exit(1);
    }
    electron_fid_mask_isLoaded = true;
  }
  return *electron_fid_mask_;
}
const vector<bool> &EGTree::electron_hasTightCharge() {
  if (not electron_hasTightCharge_isLoaded) {
    if (electron_hasTightCharge_branch != 0) {
      electron_hasTightCharge_branch->GetEntry(index);
    } else {
      printf("branch electron_hasTightCharge_branch does not exist!\n");
      exit(1);
    }
    electron_hasTightCharge_isLoaded = true;
  }
  return *electron_hasTightCharge_;
}
const vector<int> &EGTree::electron_id() {
  if (not electron_id_isLoaded) {
    if (electron_id_branch != 0) {
      electron_id_branch->GetEntry(index);
    } else {
      printf("branch electron_id_branch does not exist!\n");
      exit(1);
    }
    electron_id_isLoaded = true;
  }
  return *electron_id_;
}
const vector<int> &EGTree::electron_id_genMatch() {
  if (not electron_id_genMatch_isLoaded) {
    if (electron_id_genMatch_branch != 0) {
      electron_id_genMatch_branch->GetEntry(index);
    } else {
      printf("branch electron_id_genMatch_branch does not exist!\n");
      exit(1);
    }
    electron_id_genMatch_isLoaded = true;
  }
  return *electron_id_genMatch_;
}
const vector<bool> &EGTree::electron_is_conversionSafe() {
  if (not electron_is_conversionSafe_isLoaded) {
    if (electron_is_conversionSafe_branch != 0) {
      electron_is_conversionSafe_branch->GetEntry(index);
    } else {
      printf("branch electron_is_conversionSafe_branch does not exist!\n");
      exit(1);
    }
    electron_is_conversionSafe_isLoaded = true;
  }
  return *electron_is_conversionSafe_;
}
const vector<bool> &EGTree::electron_is_extraTight() {
  if (not electron_is_extraTight_isLoaded) {
    if (electron_is_extraTight_branch != 0) {
      electron_is_extraTight_branch->GetEntry(index);
    } else {
      printf("branch electron_is_extraTight_branch does not exist!\n");
      exit(1);
    }
    electron_is_extraTight_isLoaded = true;
  }
  return *electron_is_extraTight_;
}
const vector<bool> &EGTree::electron_is_genMatched_prompt() {
  if (not electron_is_genMatched_prompt_isLoaded) {
    if (electron_is_genMatched_prompt_branch != 0) {
      electron_is_genMatched_prompt_branch->GetEntry(index);
    } else {
      printf("branch electron_is_genMatched_prompt_branch does not exist!\n");
      exit(1);
    }
    electron_is_genMatched_prompt_isLoaded = true;
  }
  return *electron_is_genMatched_prompt_;
}
const vector<float> &EGTree::electron_minDR_electron() {
  if (not electron_minDR_electron_isLoaded) {
    if (electron_minDR_electron_branch != 0) {
      electron_minDR_electron_branch->GetEntry(index);
    } else {
      printf("branch electron_minDR_electron_branch does not exist!\n");
      exit(1);
    }
    electron_minDR_electron_isLoaded = true;
  }
  return *electron_minDR_electron_;
}
const vector<float> &EGTree::electron_minDR_muon() {
  if (not electron_minDR_muon_isLoaded) {
    if (electron_minDR_muon_branch != 0) {
      electron_minDR_muon_branch->GetEntry(index);
    } else {
      printf("branch electron_minDR_muon_branch does not exist!\n");
      exit(1);
    }
    electron_minDR_muon_isLoaded = true;
  }
  return *electron_minDR_muon_;
}
const vector<float> &EGTree::electron_minDR_photon() {
  if (not electron_minDR_photon_isLoaded) {
    if (electron_minDR_photon_branch != 0) {
      electron_minDR_photon_branch->GetEntry(index);
    } else {
      printf("branch electron_minDR_photon_branch does not exist!\n");
      exit(1);
    }
    electron_minDR_photon_isLoaded = true;
  }
  return *electron_minDR_photon_;
}
const vector<float> &EGTree::electron_phi() {
  if (not electron_phi_isLoaded) {
    if (electron_phi_branch != 0) {
      electron_phi_branch->GetEntry(index);
    } else {
      printf("branch electron_phi_branch does not exist!\n");
      exit(1);
    }
    electron_phi_isLoaded = true;
  }
  return *electron_phi_;
}
const vector<float> &EGTree::electron_pt() {
  if (not electron_pt_isLoaded) {
    if (electron_pt_branch != 0) {
      electron_pt_branch->GetEntry(index);
    } else {
      printf("branch electron_pt_branch does not exist!\n");
      exit(1);
    }
    electron_pt_isLoaded = true;
  }
  return *electron_pt_;
}
const vector<float> &EGTree::eta_eg() {
  if (not eta_eg_isLoaded) {
    if (eta_eg_branch != 0) {
      eta_eg_branch->GetEntry(index);
    } else {
      printf("branch eta_eg_branch does not exist!\n");
      exit(1);
    }
    eta_eg_isLoaded = true;
  }
  return *eta_eg_;
}
const unsigned int &EGTree::event_NGenPromptParticles() {
  if (not event_NGenPromptParticles_isLoaded) {
    if (event_NGenPromptParticles_branch != 0) {
      event_NGenPromptParticles_branch->GetEntry(index);
    } else {
      printf("branch event_NGenPromptParticles_branch does not exist!\n");
      exit(1);
    }
    event_NGenPromptParticles_isLoaded = true;
  }
  return event_NGenPromptParticles_;
}
const unsigned int &EGTree::event_n_ak4jets_pt20() {
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
const unsigned int &EGTree::event_n_ak4jets_pt20_btagged_loose() {
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
const unsigned int &EGTree::event_n_ak4jets_pt30() {
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
const unsigned int &EGTree::event_n_ak4jets_pt30_btagged_loose() {
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
const unsigned int &EGTree::event_n_leptons_fakeableBase() {
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
const unsigned int &EGTree::event_nvtxs_good() {
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
const float &EGTree::event_pTmiss() {
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
const float &EGTree::event_phimiss() {
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
const float &EGTree::event_wgt() {
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
const float &EGTree::event_wgt_SFs() {
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
const float &EGTree::genmet_pTmiss() {
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
const float &EGTree::genmet_phimiss() {
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
const vector<float> &EGTree::mass_eg() {
  if (not mass_eg_isLoaded) {
    if (mass_eg_branch != 0) {
      mass_eg_branch->GetEntry(index);
    } else {
      printf("branch mass_eg_branch does not exist!\n");
      exit(1);
    }
    mass_eg_isLoaded = true;
  }
  return *mass_eg_;
}
const vector<float> &EGTree::phi_eg() {
  if (not phi_eg_isLoaded) {
    if (phi_eg_branch != 0) {
      phi_eg_branch->GetEntry(index);
    } else {
      printf("branch phi_eg_branch does not exist!\n");
      exit(1);
    }
    phi_eg_isLoaded = true;
  }
  return *phi_eg_;
}
const vector<float> &EGTree::photon_MIPTotalEnergy() {
  if (not photon_MIPTotalEnergy_isLoaded) {
    if (photon_MIPTotalEnergy_branch != 0) {
      photon_MIPTotalEnergy_branch->GetEntry(index);
    } else {
      printf("branch photon_MIPTotalEnergy_branch does not exist!\n");
      exit(1);
    }
    photon_MIPTotalEnergy_isLoaded = true;
  }
  return *photon_MIPTotalEnergy_;
}
const vector<float> &EGTree::photon_dR_genMatch() {
  if (not photon_dR_genMatch_isLoaded) {
    if (photon_dR_genMatch_branch != 0) {
      photon_dR_genMatch_branch->GetEntry(index);
    } else {
      printf("branch photon_dR_genMatch_branch does not exist!\n");
      exit(1);
    }
    photon_dR_genMatch_isLoaded = true;
  }
  return *photon_dR_genMatch_;
}
const vector<float> &EGTree::photon_eta() {
  if (not photon_eta_isLoaded) {
    if (photon_eta_branch != 0) {
      photon_eta_branch->GetEntry(index);
    } else {
      printf("branch photon_eta_branch does not exist!\n");
      exit(1);
    }
    photon_eta_isLoaded = true;
  }
  return *photon_eta_;
}
const vector<float> &EGTree::photon_etaSC() {
  if (not photon_etaSC_isLoaded) {
    if (photon_etaSC_branch != 0) {
      photon_etaSC_branch->GetEntry(index);
    } else {
      printf("branch photon_etaSC_branch does not exist!\n");
      exit(1);
    }
    photon_etaSC_isLoaded = true;
  }
  return *photon_etaSC_;
}
const vector<unsigned short> &EGTree::photon_fid_mask() {
  if (not photon_fid_mask_isLoaded) {
    if (photon_fid_mask_branch != 0) {
      photon_fid_mask_branch->GetEntry(index);
    } else {
      printf("branch photon_fid_mask_branch does not exist!\n");
      exit(1);
    }
    photon_fid_mask_isLoaded = true;
  }
  return *photon_fid_mask_;
}
const vector<float> &EGTree::photon_full5x5_r9() {
  if (not photon_full5x5_r9_isLoaded) {
    if (photon_full5x5_r9_branch != 0) {
      photon_full5x5_r9_branch->GetEntry(index);
    } else {
      printf("branch photon_full5x5_r9_branch does not exist!\n");
      exit(1);
    }
    photon_full5x5_r9_isLoaded = true;
  }
  return *photon_full5x5_r9_;
}
const vector<float> &EGTree::photon_full5x5_sigmaIEtaIEta() {
  if (not photon_full5x5_sigmaIEtaIEta_isLoaded) {
    if (photon_full5x5_sigmaIEtaIEta_branch != 0) {
      photon_full5x5_sigmaIEtaIEta_branch->GetEntry(index);
    } else {
      printf("branch photon_full5x5_sigmaIEtaIEta_branch does not exist!\n");
      exit(1);
    }
    photon_full5x5_sigmaIEtaIEta_isLoaded = true;
  }
  return *photon_full5x5_sigmaIEtaIEta_;
}
const vector<float> &EGTree::photon_full5x5_sigmaIPhiIPhi() {
  if (not photon_full5x5_sigmaIPhiIPhi_isLoaded) {
    if (photon_full5x5_sigmaIPhiIPhi_branch != 0) {
      photon_full5x5_sigmaIPhiIPhi_branch->GetEntry(index);
    } else {
      printf("branch photon_full5x5_sigmaIPhiIPhi_branch does not exist!\n");
      exit(1);
    }
    photon_full5x5_sigmaIPhiIPhi_isLoaded = true;
  }
  return *photon_full5x5_sigmaIPhiIPhi_;
}
const vector<int> &EGTree::photon_id() {
  if (not photon_id_isLoaded) {
    if (photon_id_branch != 0) {
      photon_id_branch->GetEntry(index);
    } else {
      printf("branch photon_id_branch does not exist!\n");
      exit(1);
    }
    photon_id_isLoaded = true;
  }
  return *photon_id_;
}
const vector<int> &EGTree::photon_id_genMatch() {
  if (not photon_id_genMatch_isLoaded) {
    if (photon_id_genMatch_branch != 0) {
      photon_id_genMatch_branch->GetEntry(index);
    } else {
      printf("branch photon_id_genMatch_branch does not exist!\n");
      exit(1);
    }
    photon_id_genMatch_isLoaded = true;
  }
  return *photon_id_genMatch_;
}
const vector<bool> &EGTree::photon_is_METSafe() {
  if (not photon_is_METSafe_isLoaded) {
    if (photon_is_METSafe_branch != 0) {
      photon_is_METSafe_branch->GetEntry(index);
    } else {
      printf("branch photon_is_METSafe_branch does not exist!\n");
      exit(1);
    }
    photon_is_METSafe_isLoaded = true;
  }
  return *photon_is_METSafe_;
}
const vector<bool> &EGTree::photon_is_PFID() {
  if (not photon_is_PFID_isLoaded) {
    if (photon_is_PFID_branch != 0) {
      photon_is_PFID_branch->GetEntry(index);
    } else {
      printf("branch photon_is_PFID_branch does not exist!\n");
      exit(1);
    }
    photon_is_PFID_isLoaded = true;
  }
  return *photon_is_PFID_;
}
const vector<bool> &EGTree::photon_is_beamHaloSafe() {
  if (not photon_is_beamHaloSafe_isLoaded) {
    if (photon_is_beamHaloSafe_branch != 0) {
      photon_is_beamHaloSafe_branch->GetEntry(index);
    } else {
      printf("branch photon_is_beamHaloSafe_branch does not exist!\n");
      exit(1);
    }
    photon_is_beamHaloSafe_isLoaded = true;
  }
  return *photon_is_beamHaloSafe_;
}
const vector<bool> &EGTree::photon_is_conversionSafe() {
  if (not photon_is_conversionSafe_isLoaded) {
    if (photon_is_conversionSafe_branch != 0) {
      photon_is_conversionSafe_branch->GetEntry(index);
    } else {
      printf("branch photon_is_conversionSafe_branch does not exist!\n");
      exit(1);
    }
    photon_is_conversionSafe_isLoaded = true;
  }
  return *photon_is_conversionSafe_;
}
const vector<bool> &EGTree::photon_is_genMatched_prompt() {
  if (not photon_is_genMatched_prompt_isLoaded) {
    if (photon_is_genMatched_prompt_branch != 0) {
      photon_is_genMatched_prompt_branch->GetEntry(index);
    } else {
      printf("branch photon_is_genMatched_prompt_branch does not exist!\n");
      exit(1);
    }
    photon_is_genMatched_prompt_isLoaded = true;
  }
  return *photon_is_genMatched_prompt_;
}
const vector<bool> &EGTree::photon_is_inTime() {
  if (not photon_is_inTime_isLoaded) {
    if (photon_is_inTime_branch != 0) {
      photon_is_inTime_branch->GetEntry(index);
    } else {
      printf("branch photon_is_inTime_branch does not exist!\n");
      exit(1);
    }
    photon_is_inTime_isLoaded = true;
  }
  return *photon_is_inTime_;
}
const vector<bool> &EGTree::photon_is_spikeSafe() {
  if (not photon_is_spikeSafe_isLoaded) {
    if (photon_is_spikeSafe_branch != 0) {
      photon_is_spikeSafe_branch->GetEntry(index);
    } else {
      printf("branch photon_is_spikeSafe_branch does not exist!\n");
      exit(1);
    }
    photon_is_spikeSafe_isLoaded = true;
  }
  return *photon_is_spikeSafe_;
}
const vector<float> &EGTree::photon_minDR_electron() {
  if (not photon_minDR_electron_isLoaded) {
    if (photon_minDR_electron_branch != 0) {
      photon_minDR_electron_branch->GetEntry(index);
    } else {
      printf("branch photon_minDR_electron_branch does not exist!\n");
      exit(1);
    }
    photon_minDR_electron_isLoaded = true;
  }
  return *photon_minDR_electron_;
}
const vector<float> &EGTree::photon_minDR_muon() {
  if (not photon_minDR_muon_isLoaded) {
    if (photon_minDR_muon_branch != 0) {
      photon_minDR_muon_branch->GetEntry(index);
    } else {
      printf("branch photon_minDR_muon_branch does not exist!\n");
      exit(1);
    }
    photon_minDR_muon_isLoaded = true;
  }
  return *photon_minDR_muon_;
}
const vector<float> &EGTree::photon_minDR_photon() {
  if (not photon_minDR_photon_isLoaded) {
    if (photon_minDR_photon_branch != 0) {
      photon_minDR_photon_branch->GetEntry(index);
    } else {
      printf("branch photon_minDR_photon_branch does not exist!\n");
      exit(1);
    }
    photon_minDR_photon_isLoaded = true;
  }
  return *photon_minDR_photon_;
}
const vector<float> &EGTree::photon_phi() {
  if (not photon_phi_isLoaded) {
    if (photon_phi_branch != 0) {
      photon_phi_branch->GetEntry(index);
    } else {
      printf("branch photon_phi_branch does not exist!\n");
      exit(1);
    }
    photon_phi_isLoaded = true;
  }
  return *photon_phi_;
}
const vector<float> &EGTree::photon_pt() {
  if (not photon_pt_isLoaded) {
    if (photon_pt_branch != 0) {
      photon_pt_branch->GetEntry(index);
    } else {
      printf("branch photon_pt_branch does not exist!\n");
      exit(1);
    }
    photon_pt_isLoaded = true;
  }
  return *photon_pt_;
}
const vector<float> &EGTree::photon_seedTime() {
  if (not photon_seedTime_isLoaded) {
    if (photon_seedTime_branch != 0) {
      photon_seedTime_branch->GetEntry(index);
    } else {
      printf("branch photon_seedTime_branch does not exist!\n");
      exit(1);
    }
    photon_seedTime_isLoaded = true;
  }
  return *photon_seedTime_;
}
const vector<float> &EGTree::pt_eg() {
  if (not pt_eg_isLoaded) {
    if (pt_eg_branch != 0) {
      pt_eg_branch->GetEntry(index);
    } else {
      printf("branch pt_eg_branch does not exist!\n");
      exit(1);
    }
    pt_eg_isLoaded = true;
  }
  return *pt_eg_;
}
const vector<float> &EGTree::weight_HLT_Photon120_R9Id90_HE10_IsoM() {
  if (not weight_HLT_Photon120_R9Id90_HE10_IsoM_isLoaded) {
    if (weight_HLT_Photon120_R9Id90_HE10_IsoM_branch != 0) {
      weight_HLT_Photon120_R9Id90_HE10_IsoM_branch->GetEntry(index);
    } else {
      printf("branch weight_HLT_Photon120_R9Id90_HE10_IsoM_branch does not exist!\n");
      exit(1);
    }
    weight_HLT_Photon120_R9Id90_HE10_IsoM_isLoaded = true;
  }
  return *weight_HLT_Photon120_R9Id90_HE10_IsoM_;
}
const vector<float> &EGTree::weight_HLT_Photon165_R9Id90_HE10_IsoM() {
  if (not weight_HLT_Photon165_R9Id90_HE10_IsoM_isLoaded) {
    if (weight_HLT_Photon165_R9Id90_HE10_IsoM_branch != 0) {
      weight_HLT_Photon165_R9Id90_HE10_IsoM_branch->GetEntry(index);
    } else {
      printf("branch weight_HLT_Photon165_R9Id90_HE10_IsoM_branch does not exist!\n");
      exit(1);
    }
    weight_HLT_Photon165_R9Id90_HE10_IsoM_isLoaded = true;
  }
  return *weight_HLT_Photon165_R9Id90_HE10_IsoM_;
}
const vector<float> &EGTree::weight_HLT_Photon200() {
  if (not weight_HLT_Photon200_isLoaded) {
    if (weight_HLT_Photon200_branch != 0) {
      weight_HLT_Photon200_branch->GetEntry(index);
    } else {
      printf("branch weight_HLT_Photon200_branch does not exist!\n");
      exit(1);
    }
    weight_HLT_Photon200_isLoaded = true;
  }
  return *weight_HLT_Photon200_;
}
const vector<float> &EGTree::weight_HLT_Photon175() {
  if (not weight_HLT_Photon175_isLoaded) {
    if (weight_HLT_Photon175_branch != 0) {
      weight_HLT_Photon175_branch->GetEntry(index);
    } else {
      printf("branch weight_HLT_Photon175_branch does not exist!\n");
      exit(1);
    }
    weight_HLT_Photon175_isLoaded = true;
  }
  return *weight_HLT_Photon175_;
}
const vector<float> &EGTree::weight_HLT_Photon20_HoverELoose() {
  if (not weight_HLT_Photon20_HoverELoose_isLoaded) {
    if (weight_HLT_Photon20_HoverELoose_branch != 0) {
      weight_HLT_Photon20_HoverELoose_branch->GetEntry(index);
    } else {
      printf("branch weight_HLT_Photon20_HoverELoose_branch does not exist!\n");
      exit(1);
    }
    weight_HLT_Photon20_HoverELoose_isLoaded = true;
  }
  return *weight_HLT_Photon20_HoverELoose_;
}
const vector<float> &EGTree::weight_HLT_Photon30_HoverELoose() {
  if (not weight_HLT_Photon30_HoverELoose_isLoaded) {
    if (weight_HLT_Photon30_HoverELoose_branch != 0) {
      weight_HLT_Photon30_HoverELoose_branch->GetEntry(index);
    } else {
      printf("branch weight_HLT_Photon30_HoverELoose_branch does not exist!\n");
      exit(1);
    }
    weight_HLT_Photon30_HoverELoose_isLoaded = true;
  }
  return *weight_HLT_Photon30_HoverELoose_;
}
const vector<float> &EGTree::weight_HLT_Photon33() {
  if (not weight_HLT_Photon33_isLoaded) {
    if (weight_HLT_Photon33_branch != 0) {
      weight_HLT_Photon33_branch->GetEntry(index);
    } else {
      printf("branch weight_HLT_Photon33_branch does not exist!\n");
      exit(1);
    }
    weight_HLT_Photon33_isLoaded = true;
  }
  return *weight_HLT_Photon33_;
}
const vector<float> &EGTree::weight_HLT_Photon50_R9Id90_HE10_IsoM() {
  if (not weight_HLT_Photon50_R9Id90_HE10_IsoM_isLoaded) {
    if (weight_HLT_Photon50_R9Id90_HE10_IsoM_branch != 0) {
      weight_HLT_Photon50_R9Id90_HE10_IsoM_branch->GetEntry(index);
    } else {
      printf("branch weight_HLT_Photon50_R9Id90_HE10_IsoM_branch does not exist!\n");
      exit(1);
    }
    weight_HLT_Photon50_R9Id90_HE10_IsoM_isLoaded = true;
  }
  return *weight_HLT_Photon50_R9Id90_HE10_IsoM_;
}
const vector<float> &EGTree::weight_HLT_Photon75_R9Id90_HE10_IsoM() {
  if (not weight_HLT_Photon75_R9Id90_HE10_IsoM_isLoaded) {
    if (weight_HLT_Photon75_R9Id90_HE10_IsoM_branch != 0) {
      weight_HLT_Photon75_R9Id90_HE10_IsoM_branch->GetEntry(index);
    } else {
      printf("branch weight_HLT_Photon75_R9Id90_HE10_IsoM_branch does not exist!\n");
      exit(1);
    }
    weight_HLT_Photon75_R9Id90_HE10_IsoM_isLoaded = true;
  }
  return *weight_HLT_Photon75_R9Id90_HE10_IsoM_;
}
const vector<float> &EGTree::weight_HLT_Photon90_R9Id90_HE10_IsoM() {
  if (not weight_HLT_Photon90_R9Id90_HE10_IsoM_isLoaded) {
    if (weight_HLT_Photon90_R9Id90_HE10_IsoM_branch != 0) {
      weight_HLT_Photon90_R9Id90_HE10_IsoM_branch->GetEntry(index);
    } else {
      printf("branch weight_HLT_Photon90_R9Id90_HE10_IsoM_branch does not exist!\n");
      exit(1);
    }
    weight_HLT_Photon90_R9Id90_HE10_IsoM_isLoaded = true;
  }
  return *weight_HLT_Photon90_R9Id90_HE10_IsoM_;
}
void EGTree::progress( int nEventsTotal, int nEventsChain ){
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
  const vector<float> &dR_e_g() { return st.dR_e_g(); }
  const vector<float> &electron_dR_genMatch() { return st.electron_dR_genMatch(); }
  const vector<float> &electron_dxy() { return st.electron_dxy(); }
  const vector<float> &electron_dz() { return st.electron_dz(); }
  const vector<float> &electron_eta() { return st.electron_eta(); }
  const vector<float> &electron_etaSC() { return st.electron_etaSC(); }
  const vector<unsigned short> &electron_fid_mask() { return st.electron_fid_mask(); }
  const vector<bool> &electron_hasTightCharge() { return st.electron_hasTightCharge(); }
  const vector<int> &electron_id() { return st.electron_id(); }
  const vector<int> &electron_id_genMatch() { return st.electron_id_genMatch(); }
  const vector<bool> &electron_is_conversionSafe() { return st.electron_is_conversionSafe(); }
  const vector<bool> &electron_is_extraTight() { return st.electron_is_extraTight(); }
  const vector<bool> &electron_is_genMatched_prompt() { return st.electron_is_genMatched_prompt(); }
  const vector<float> &electron_minDR_electron() { return st.electron_minDR_electron(); }
  const vector<float> &electron_minDR_muon() { return st.electron_minDR_muon(); }
  const vector<float> &electron_minDR_photon() { return st.electron_minDR_photon(); }
  const vector<float> &electron_phi() { return st.electron_phi(); }
  const vector<float> &electron_pt() { return st.electron_pt(); }
  const vector<float> &eta_eg() { return st.eta_eg(); }
  const unsigned int &event_NGenPromptParticles() { return st.event_NGenPromptParticles(); }
  const unsigned int &event_n_ak4jets_pt20() { return st.event_n_ak4jets_pt20(); }
  const unsigned int &event_n_ak4jets_pt20_btagged_loose() { return st.event_n_ak4jets_pt20_btagged_loose(); }
  const unsigned int &event_n_ak4jets_pt30() { return st.event_n_ak4jets_pt30(); }
  const unsigned int &event_n_ak4jets_pt30_btagged_loose() { return st.event_n_ak4jets_pt30_btagged_loose(); }
  const unsigned int &event_n_leptons_fakeableBase() { return st.event_n_leptons_fakeableBase(); }
  const unsigned int &event_nvtxs_good() { return st.event_nvtxs_good(); }
  const float &event_pTmiss() { return st.event_pTmiss(); }
  const float &event_phimiss() { return st.event_phimiss(); }
  const float &event_wgt() { return st.event_wgt(); }
  const float &event_wgt_SFs() { return st.event_wgt_SFs(); }
  const float &genmet_pTmiss() { return st.genmet_pTmiss(); }
  const float &genmet_phimiss() { return st.genmet_phimiss(); }
  const vector<float> &mass_eg() { return st.mass_eg(); }
  const vector<float> &phi_eg() { return st.phi_eg(); }
  const vector<float> &photon_MIPTotalEnergy() { return st.photon_MIPTotalEnergy(); }
  const vector<float> &photon_dR_genMatch() { return st.photon_dR_genMatch(); }
  const vector<float> &photon_eta() { return st.photon_eta(); }
  const vector<float> &photon_etaSC() { return st.photon_etaSC(); }
  const vector<unsigned short> &photon_fid_mask() { return st.photon_fid_mask(); }
  const vector<float> &photon_full5x5_r9() { return st.photon_full5x5_r9(); }
  const vector<float> &photon_full5x5_sigmaIEtaIEta() { return st.photon_full5x5_sigmaIEtaIEta(); }
  const vector<float> &photon_full5x5_sigmaIPhiIPhi() { return st.photon_full5x5_sigmaIPhiIPhi(); }
  const vector<int> &photon_id() { return st.photon_id(); }
  const vector<int> &photon_id_genMatch() { return st.photon_id_genMatch(); }
  const vector<bool> &photon_is_METSafe() { return st.photon_is_METSafe(); }
  const vector<bool> &photon_is_PFID() { return st.photon_is_PFID(); }
  const vector<bool> &photon_is_beamHaloSafe() { return st.photon_is_beamHaloSafe(); }
  const vector<bool> &photon_is_conversionSafe() { return st.photon_is_conversionSafe(); }
  const vector<bool> &photon_is_genMatched_prompt() { return st.photon_is_genMatched_prompt(); }
  const vector<bool> &photon_is_inTime() { return st.photon_is_inTime(); }
  const vector<bool> &photon_is_spikeSafe() { return st.photon_is_spikeSafe(); }
  const vector<float> &photon_minDR_electron() { return st.photon_minDR_electron(); }
  const vector<float> &photon_minDR_muon() { return st.photon_minDR_muon(); }
  const vector<float> &photon_minDR_photon() { return st.photon_minDR_photon(); }
  const vector<float> &photon_phi() { return st.photon_phi(); }
  const vector<float> &photon_pt() { return st.photon_pt(); }
  const vector<float> &photon_seedTime() { return st.photon_seedTime(); }
  const vector<float> &pt_eg() { return st.pt_eg(); }
  const vector<float> &weight_HLT_Photon120_R9Id90_HE10_IsoM() { return st.weight_HLT_Photon120_R9Id90_HE10_IsoM(); }
  const vector<float> &weight_HLT_Photon165_R9Id90_HE10_IsoM() { return st.weight_HLT_Photon165_R9Id90_HE10_IsoM(); }
  const vector<float> &weight_HLT_Photon200() { return st.weight_HLT_Photon200(); }
  const vector<float> &weight_HLT_Photon175() { return st.weight_HLT_Photon175(); }
  const vector<float> &weight_HLT_Photon20_HoverELoose() { return st.weight_HLT_Photon20_HoverELoose(); }
  const vector<float> &weight_HLT_Photon30_HoverELoose() { return st.weight_HLT_Photon30_HoverELoose(); }
  const vector<float> &weight_HLT_Photon33() { return st.weight_HLT_Photon33(); }
  const vector<float> &weight_HLT_Photon50_R9Id90_HE10_IsoM() { return st.weight_HLT_Photon50_R9Id90_HE10_IsoM(); }
  const vector<float> &weight_HLT_Photon75_R9Id90_HE10_IsoM() { return st.weight_HLT_Photon75_R9Id90_HE10_IsoM(); }
  const vector<float> &weight_HLT_Photon90_R9Id90_HE10_IsoM() { return st.weight_HLT_Photon90_R9Id90_HE10_IsoM(); }
}
