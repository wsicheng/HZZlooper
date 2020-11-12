// -*- C++ -*-
#ifndef SkimTree_H
#define SkimTree_H
#include "Math/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"
#include "Math/Point3D.h"
#include "TMath.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TBits.h"
#include <vector>
#include <unistd.h>
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> > LorentzVector;

// Generated with command: -t SkimTree -n tas -o st -c SkimTree -plv /hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/SingleLeptonEvents/SkimTrees/201108/AK4Jets_WithPUJetId_NoTightLeptonJetId_ParticleCleaned/PFMET_WithXY_WithPartMomCorr_P4Preserved/DefaultLeptons/2017E/EGamma_2017E_Nominal.root 

using namespace std;
class SkimTree {
private:
protected:
  unsigned int index;
  unsigned long long EventNumber_;
  TBranch *EventNumber_branch;
  bool EventNumber_isLoaded;
  unsigned int LuminosityBlock_;
  TBranch *LuminosityBlock_branch;
  bool LuminosityBlock_isLoaded;
  unsigned int RunNumber_;
  TBranch *RunNumber_branch;
  bool RunNumber_isLoaded;
  vector<float> *ak4jets_CEMF_;
  TBranch *ak4jets_CEMF_branch;
  bool ak4jets_CEMF_isLoaded;
  float ak4jets_HT_;
  TBranch *ak4jets_HT_branch;
  bool ak4jets_HT_isLoaded;
  float ak4jets_MHT_;
  TBranch *ak4jets_MHT_branch;
  bool ak4jets_MHT_isLoaded;
  vector<float> *ak4jets_NEMF_;
  TBranch *ak4jets_NEMF_branch;
  bool ak4jets_NEMF_isLoaded;
  vector<unsigned char> *ak4jets_btagWP_Bits_;
  TBranch *ak4jets_btagWP_Bits_branch;
  bool ak4jets_btagWP_Bits_isLoaded;
  vector<float> *ak4jets_eta_;
  TBranch *ak4jets_eta_branch;
  bool ak4jets_eta_isLoaded;
  vector<bool> *ak4jets_is_genMatched_;
  TBranch *ak4jets_is_genMatched_branch;
  bool ak4jets_is_genMatched_isLoaded;
  vector<bool> *ak4jets_is_genMatched_fullCone_;
  TBranch *ak4jets_is_genMatched_fullCone_branch;
  bool ak4jets_is_genMatched_fullCone_isLoaded;
  vector<float> *ak4jets_mass_;
  TBranch *ak4jets_mass_branch;
  bool ak4jets_mass_isLoaded;
  vector<float> *ak4jets_phi_;
  TBranch *ak4jets_phi_branch;
  bool ak4jets_phi_isLoaded;
  vector<float> *ak4jets_pt_;
  TBranch *ak4jets_pt_branch;
  bool ak4jets_pt_isLoaded;
  vector<float> *ak8jets_eta_;
  TBranch *ak8jets_eta_branch;
  bool ak8jets_eta_isLoaded;
  vector<float> *ak8jets_mass_;
  TBranch *ak8jets_mass_branch;
  bool ak8jets_mass_isLoaded;
  vector<float> *ak8jets_phi_;
  TBranch *ak8jets_phi_branch;
  bool ak8jets_phi_isLoaded;
  vector<float> *ak8jets_pt_;
  TBranch *ak8jets_pt_branch;
  bool ak8jets_pt_isLoaded;
  float dPhi_pTboson_pTmiss_;
  TBranch *dPhi_pTboson_pTmiss_branch;
  bool dPhi_pTboson_pTmiss_isLoaded;
  float dPhi_pTbosonjets_pTmiss_;
  TBranch *dPhi_pTbosonjets_pTmiss_branch;
  bool dPhi_pTbosonjets_pTmiss_isLoaded;
  float dilepton_eta_;
  TBranch *dilepton_eta_branch;
  bool dilepton_eta_isLoaded;
  int dilepton_id_;
  TBranch *dilepton_id_branch;
  bool dilepton_id_isLoaded;
  float dilepton_mass_;
  TBranch *dilepton_mass_branch;
  bool dilepton_mass_isLoaded;
  float dilepton_phi_;
  TBranch *dilepton_phi_branch;
  bool dilepton_phi_isLoaded;
  float dilepton_pt_;
  TBranch *dilepton_pt_branch;
  bool dilepton_pt_isLoaded;
  float electron_full5x5_r9_;
  TBranch *electron_full5x5_r9_branch;
  bool electron_full5x5_r9_isLoaded;
  float electron_full5x5_sigmaIEtaIEta_;
  TBranch *electron_full5x5_sigmaIEtaIEta_branch;
  bool electron_full5x5_sigmaIEtaIEta_isLoaded;
  float electron_full5x5_sigmaIPhiIPhi_;
  TBranch *electron_full5x5_sigmaIPhiIPhi_branch;
  bool electron_full5x5_sigmaIPhiIPhi_isLoaded;
  float electron_seedTime_;
  TBranch *electron_seedTime_branch;
  bool electron_seedTime_isLoaded;
  vector<float> *electrons_full5x5_r9_;
  TBranch *electrons_full5x5_r9_branch;
  bool electrons_full5x5_r9_isLoaded;
  vector<float> *electrons_full5x5_sigmaIEtaIEta_;
  TBranch *electrons_full5x5_sigmaIEtaIEta_branch;
  bool electrons_full5x5_sigmaIEtaIEta_isLoaded;
  vector<float> *electrons_full5x5_sigmaIPhiIPhi_;
  TBranch *electrons_full5x5_sigmaIPhiIPhi_branch;
  bool electrons_full5x5_sigmaIPhiIPhi_isLoaded;
  vector<float> *electrons_seedTime_;
  TBranch *electrons_seedTime_branch;
  bool electrons_seedTime_isLoaded;
  float event_mTZZ_;
  TBranch *event_mTZZ_branch;
  bool event_mTZZ_isLoaded;
  float event_mZZ_;
  TBranch *event_mZZ_branch;
  bool event_mZZ_isLoaded;
  float event_mllg_;
  TBranch *event_mllg_branch;
  bool event_mllg_isLoaded;
  unsigned int event_n_ak4jets_pt20_;
  TBranch *event_n_ak4jets_pt20_branch;
  bool event_n_ak4jets_pt20_isLoaded;
  unsigned int event_n_ak4jets_pt20_btagged_loose_;
  TBranch *event_n_ak4jets_pt20_btagged_loose_branch;
  bool event_n_ak4jets_pt20_btagged_loose_isLoaded;
  unsigned int event_n_ak4jets_pt20_btagged_medium_;
  TBranch *event_n_ak4jets_pt20_btagged_medium_branch;
  bool event_n_ak4jets_pt20_btagged_medium_isLoaded;
  unsigned int event_n_ak4jets_pt30_;
  TBranch *event_n_ak4jets_pt30_branch;
  bool event_n_ak4jets_pt30_isLoaded;
  unsigned int event_n_ak4jets_pt30_btagged_loose_;
  TBranch *event_n_ak4jets_pt30_btagged_loose_branch;
  bool event_n_ak4jets_pt30_btagged_loose_isLoaded;
  unsigned int event_n_ak4jets_pt30_btagged_medium_;
  TBranch *event_n_ak4jets_pt30_btagged_medium_branch;
  bool event_n_ak4jets_pt30_btagged_medium_isLoaded;
  unsigned int event_n_vtxs_good_;
  TBranch *event_n_vtxs_good_branch;
  bool event_n_vtxs_good_isLoaded;
  float event_pTmiss_;
  TBranch *event_pTmiss_branch;
  bool event_pTmiss_isLoaded;
  bool event_pass_tightMETFilters_;
  TBranch *event_pass_tightMETFilters_branch;
  bool event_pass_tightMETFilters_isLoaded;
  float event_phimiss_;
  TBranch *event_phimiss_branch;
  bool event_phimiss_isLoaded;
  float event_wgt_;
  TBranch *event_wgt_branch;
  bool event_wgt_isLoaded;
  float event_wgt_L1PrefiringDn_;
  TBranch *event_wgt_L1PrefiringDn_branch;
  bool event_wgt_L1PrefiringDn_isLoaded;
  float event_wgt_L1PrefiringUp_;
  TBranch *event_wgt_L1PrefiringUp_branch;
  bool event_wgt_L1PrefiringUp_isLoaded;
  float event_wgt_SFs_PUJetId_;
  TBranch *event_wgt_SFs_PUJetId_branch;
  bool event_wgt_SFs_PUJetId_isLoaded;
  float event_wgt_SFs_PUJetId_EffDn_;
  TBranch *event_wgt_SFs_PUJetId_EffDn_branch;
  bool event_wgt_SFs_PUJetId_EffDn_isLoaded;
  float event_wgt_SFs_PUJetId_EffUp_;
  TBranch *event_wgt_SFs_PUJetId_EffUp_branch;
  bool event_wgt_SFs_PUJetId_EffUp_isLoaded;
  float event_wgt_SFs_btagging_;
  TBranch *event_wgt_SFs_btagging_branch;
  bool event_wgt_SFs_btagging_isLoaded;
  float event_wgt_SFs_btagging_EffDn_;
  TBranch *event_wgt_SFs_btagging_EffDn_branch;
  bool event_wgt_SFs_btagging_EffDn_isLoaded;
  float event_wgt_SFs_btagging_EffUp_;
  TBranch *event_wgt_SFs_btagging_EffUp_branch;
  bool event_wgt_SFs_btagging_EffUp_isLoaded;
  float event_wgt_SFs_electrons_;
  TBranch *event_wgt_SFs_electrons_branch;
  bool event_wgt_SFs_electrons_isLoaded;
  float event_wgt_SFs_electrons_AltMCDn_;
  TBranch *event_wgt_SFs_electrons_AltMCDn_branch;
  bool event_wgt_SFs_electrons_AltMCDn_isLoaded;
  float event_wgt_SFs_electrons_AltMCUp_;
  TBranch *event_wgt_SFs_electrons_AltMCUp_branch;
  bool event_wgt_SFs_electrons_AltMCUp_isLoaded;
  float event_wgt_SFs_electrons_StatDn_;
  TBranch *event_wgt_SFs_electrons_StatDn_branch;
  bool event_wgt_SFs_electrons_StatDn_isLoaded;
  float event_wgt_SFs_electrons_StatUp_;
  TBranch *event_wgt_SFs_electrons_StatUp_branch;
  bool event_wgt_SFs_electrons_StatUp_isLoaded;
  float event_wgt_SFs_electrons_SystDn_;
  TBranch *event_wgt_SFs_electrons_SystDn_branch;
  bool event_wgt_SFs_electrons_SystDn_isLoaded;
  float event_wgt_SFs_electrons_SystUp_;
  TBranch *event_wgt_SFs_electrons_SystUp_branch;
  bool event_wgt_SFs_electrons_SystUp_isLoaded;
  float event_wgt_SFs_muons_;
  TBranch *event_wgt_SFs_muons_branch;
  bool event_wgt_SFs_muons_isLoaded;
  float event_wgt_SFs_muons_AltMCDn_;
  TBranch *event_wgt_SFs_muons_AltMCDn_branch;
  bool event_wgt_SFs_muons_AltMCDn_isLoaded;
  float event_wgt_SFs_muons_AltMCUp_;
  TBranch *event_wgt_SFs_muons_AltMCUp_branch;
  bool event_wgt_SFs_muons_AltMCUp_isLoaded;
  float event_wgt_SFs_muons_StatDn_;
  TBranch *event_wgt_SFs_muons_StatDn_branch;
  bool event_wgt_SFs_muons_StatDn_isLoaded;
  float event_wgt_SFs_muons_StatUp_;
  TBranch *event_wgt_SFs_muons_StatUp_branch;
  bool event_wgt_SFs_muons_StatUp_isLoaded;
  float event_wgt_SFs_muons_SystDn_;
  TBranch *event_wgt_SFs_muons_SystDn_branch;
  bool event_wgt_SFs_muons_SystDn_isLoaded;
  float event_wgt_SFs_muons_SystUp_;
  TBranch *event_wgt_SFs_muons_SystUp_branch;
  bool event_wgt_SFs_muons_SystUp_isLoaded;
  float event_wgt_SFs_photons_;
  TBranch *event_wgt_SFs_photons_branch;
  bool event_wgt_SFs_photons_isLoaded;
  float event_wgt_SFs_photons_EffDn_;
  TBranch *event_wgt_SFs_photons_EffDn_branch;
  bool event_wgt_SFs_photons_EffDn_isLoaded;
  float event_wgt_SFs_photons_EffUp_;
  TBranch *event_wgt_SFs_photons_EffUp_branch;
  bool event_wgt_SFs_photons_EffUp_isLoaded;
  float event_wgt_adjustment_NNPDF30_;
  TBranch *event_wgt_adjustment_NNPDF30_branch;
  bool event_wgt_adjustment_NNPDF30_isLoaded;
  float event_wgt_triggers_;
  TBranch *event_wgt_triggers_branch;
  bool event_wgt_triggers_isLoaded;
  float event_wgt_triggers_Dilepton_;
  TBranch *event_wgt_triggers_Dilepton_branch;
  bool event_wgt_triggers_Dilepton_isLoaded;
  float event_wgt_triggers_SingleLepton_;
  TBranch *event_wgt_triggers_SingleLepton_branch;
  bool event_wgt_triggers_SingleLepton_isLoaded;
  float event_wgt_triggers_SinglePhoton_;
  TBranch *event_wgt_triggers_SinglePhoton_branch;
  bool event_wgt_triggers_SinglePhoton_isLoaded;
  float genmet_pTmiss_;
  TBranch *genmet_pTmiss_branch;
  bool genmet_pTmiss_isLoaded;
  float genmet_phimiss_;
  TBranch *genmet_phimiss_branch;
  bool genmet_phimiss_isLoaded;
  float lepton_eta_;
  TBranch *lepton_eta_branch;
  bool lepton_eta_isLoaded;
  int lepton_id_;
  TBranch *lepton_id_branch;
  bool lepton_id_isLoaded;
  bool lepton_is_genMatched_prompt_;
  TBranch *lepton_is_genMatched_prompt_branch;
  bool lepton_is_genMatched_prompt_isLoaded;
  float lepton_mass_;
  TBranch *lepton_mass_branch;
  bool lepton_mass_isLoaded;
  float lepton_phi_;
  TBranch *lepton_phi_branch;
  bool lepton_phi_isLoaded;
  float lepton_pt_;
  TBranch *lepton_pt_branch;
  bool lepton_pt_isLoaded;
  float lepton_relIso_;
  TBranch *lepton_relIso_branch;
  bool lepton_relIso_isLoaded;
  vector<float> *leptons_eff_;
  TBranch *leptons_eff_branch;
  bool leptons_eff_isLoaded;
  vector<float> *leptons_eff_StatDn_;
  TBranch *leptons_eff_StatDn_branch;
  bool leptons_eff_StatDn_isLoaded;
  vector<float> *leptons_eff_StatUp_;
  TBranch *leptons_eff_StatUp_branch;
  bool leptons_eff_StatUp_isLoaded;
  vector<float> *leptons_eff_SystDn_;
  TBranch *leptons_eff_SystDn_branch;
  bool leptons_eff_SystDn_isLoaded;
  vector<float> *leptons_eff_SystUp_;
  TBranch *leptons_eff_SystUp_branch;
  bool leptons_eff_SystUp_isLoaded;
  vector<float> *leptons_eta_;
  TBranch *leptons_eta_branch;
  bool leptons_eta_isLoaded;
  vector<int> *leptons_id_;
  TBranch *leptons_id_branch;
  bool leptons_id_isLoaded;
  vector<bool> *leptons_is_TOmatched_SingleLepton_;
  TBranch *leptons_is_TOmatched_SingleLepton_branch;
  bool leptons_is_TOmatched_SingleLepton_isLoaded;
  vector<bool> *leptons_is_genMatched_prompt_;
  TBranch *leptons_is_genMatched_prompt_branch;
  bool leptons_is_genMatched_prompt_isLoaded;
  vector<float> *leptons_mass_;
  TBranch *leptons_mass_branch;
  bool leptons_mass_isLoaded;
  vector<float> *leptons_phi_;
  TBranch *leptons_phi_branch;
  bool leptons_phi_isLoaded;
  vector<float> *leptons_pt_;
  TBranch *leptons_pt_branch;
  bool leptons_pt_isLoaded;
  float min_abs_dPhi_pTj_pTmiss_;
  TBranch *min_abs_dPhi_pTj_pTmiss_branch;
  bool min_abs_dPhi_pTj_pTmiss_isLoaded;
  float pConst_JJQCD_SIG_ghg2_1_JHUGen_;
  TBranch *pConst_JJQCD_SIG_ghg2_1_JHUGen_branch;
  bool pConst_JJQCD_SIG_ghg2_1_JHUGen_isLoaded;
  float pConst_JJVBF_SIG_ghv1_1_JHUGen_;
  TBranch *pConst_JJVBF_SIG_ghv1_1_JHUGen_branch;
  bool pConst_JJVBF_SIG_ghv1_1_JHUGen_isLoaded;
  float p_JJQCD_SIG_ghg2_1_JHUGen_;
  TBranch *p_JJQCD_SIG_ghg2_1_JHUGen_branch;
  bool p_JJQCD_SIG_ghg2_1_JHUGen_isLoaded;
  float p_JJVBF_SIG_ghv1_1_JHUGen_;
  TBranch *p_JJVBF_SIG_ghv1_1_JHUGen_branch;
  bool p_JJVBF_SIG_ghv1_1_JHUGen_isLoaded;
  float p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_;
  TBranch *p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_branch;
  bool p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_isLoaded;
  float p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_;
  TBranch *p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_branch;
  bool p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_isLoaded;
  float p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_;
  TBranch *p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_branch;
  bool p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_isLoaded;
  float p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_;
  TBranch *p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_branch;
  bool p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_isLoaded;
  float p_JJVBF_SIG_ghv2_1_JHUGen_;
  TBranch *p_JJVBF_SIG_ghv2_1_JHUGen_branch;
  bool p_JJVBF_SIG_ghv2_1_JHUGen_isLoaded;
  float p_JJVBF_SIG_ghv4_1_JHUGen_;
  TBranch *p_JJVBF_SIG_ghv4_1_JHUGen_branch;
  bool p_JJVBF_SIG_ghv4_1_JHUGen_isLoaded;
  float p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_;
  TBranch *p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_branch;
  bool p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_isLoaded;
  float photon_MIPTotalEnergy_;
  TBranch *photon_MIPTotalEnergy_branch;
  bool photon_MIPTotalEnergy_isLoaded;
  float photon_eta_;
  TBranch *photon_eta_branch;
  bool photon_eta_isLoaded;
  float photon_full5x5_r9_;
  TBranch *photon_full5x5_r9_branch;
  bool photon_full5x5_r9_isLoaded;
  float photon_full5x5_sigmaIEtaIEta_;
  TBranch *photon_full5x5_sigmaIEtaIEta_branch;
  bool photon_full5x5_sigmaIEtaIEta_isLoaded;
  float photon_full5x5_sigmaIPhiIPhi_;
  TBranch *photon_full5x5_sigmaIPhiIPhi_branch;
  bool photon_full5x5_sigmaIPhiIPhi_isLoaded;
  bool photon_isEB_;
  TBranch *photon_isEB_branch;
  bool photon_isEB_isLoaded;
  bool photon_isEBEEGap_;
  TBranch *photon_isEBEEGap_branch;
  bool photon_isEBEEGap_isLoaded;
  bool photon_isEE_;
  TBranch *photon_isEE_branch;
  bool photon_isEE_isLoaded;
  bool photon_isGap_;
  TBranch *photon_isGap_branch;
  bool photon_isGap_isLoaded;
  bool photon_is_METSafe_;
  TBranch *photon_is_METSafe_branch;
  bool photon_is_METSafe_isLoaded;
  bool photon_is_PFID_;
  TBranch *photon_is_PFID_branch;
  bool photon_is_PFID_isLoaded;
  bool photon_is_beamHaloSafe_;
  TBranch *photon_is_beamHaloSafe_branch;
  bool photon_is_beamHaloSafe_isLoaded;
  bool photon_is_conversionSafe_;
  TBranch *photon_is_conversionSafe_branch;
  bool photon_is_conversionSafe_isLoaded;
  bool photon_is_genMatched_prompt_;
  TBranch *photon_is_genMatched_prompt_branch;
  bool photon_is_genMatched_prompt_isLoaded;
  bool photon_is_inTime_;
  TBranch *photon_is_inTime_branch;
  bool photon_is_inTime_isLoaded;
  bool photon_is_spikeSafe_;
  TBranch *photon_is_spikeSafe_branch;
  bool photon_is_spikeSafe_isLoaded;
  float photon_mass_;
  TBranch *photon_mass_branch;
  bool photon_mass_isLoaded;
  bool photon_pass_HGGSelection_;
  TBranch *photon_pass_HGGSelection_branch;
  bool photon_pass_HGGSelection_isLoaded;
  float photon_phi_;
  TBranch *photon_phi_branch;
  bool photon_phi_isLoaded;
  float photon_pt_;
  TBranch *photon_pt_branch;
  bool photon_pt_isLoaded;
  float photon_seedTime_;
  TBranch *photon_seedTime_branch;
  bool photon_seedTime_isLoaded;
public:
  void Init(TTree *tree);
  void GetEntry(unsigned int idx);
  void LoadAllBranches();
  const unsigned long long &EventNumber();
  const unsigned int &LuminosityBlock();
  const unsigned int &RunNumber();
  const vector<float> &ak4jets_CEMF();
  const float &ak4jets_HT();
  const float &ak4jets_MHT();
  const vector<float> &ak4jets_NEMF();
  const vector<unsigned char> &ak4jets_btagWP_Bits();
  const vector<float> &ak4jets_eta();
  const vector<bool> &ak4jets_is_genMatched();
  const vector<bool> &ak4jets_is_genMatched_fullCone();
  const vector<float> &ak4jets_mass();
  const vector<float> &ak4jets_phi();
  const vector<float> &ak4jets_pt();
  const vector<float> &ak8jets_eta();
  const vector<float> &ak8jets_mass();
  const vector<float> &ak8jets_phi();
  const vector<float> &ak8jets_pt();
  const float &dPhi_pTboson_pTmiss();
  const float &dPhi_pTbosonjets_pTmiss();
  const float &dilepton_eta();
  const int &dilepton_id();
  const float &dilepton_mass();
  const float &dilepton_phi();
  const float &dilepton_pt();
  const float &electron_full5x5_r9();
  const float &electron_full5x5_sigmaIEtaIEta();
  const float &electron_full5x5_sigmaIPhiIPhi();
  const float &electron_seedTime();
  const vector<float> &electrons_full5x5_r9();
  const vector<float> &electrons_full5x5_sigmaIEtaIEta();
  const vector<float> &electrons_full5x5_sigmaIPhiIPhi();
  const vector<float> &electrons_seedTime();
  const float &event_mTZZ();
  const float &event_mZZ();
  const float &event_mllg();
  const unsigned int &event_n_ak4jets_pt20();
  const unsigned int &event_n_ak4jets_pt20_btagged_loose();
  const unsigned int &event_n_ak4jets_pt20_btagged_medium();
  const unsigned int &event_n_ak4jets_pt30();
  const unsigned int &event_n_ak4jets_pt30_btagged_loose();
  const unsigned int &event_n_ak4jets_pt30_btagged_medium();
  const unsigned int &event_n_vtxs_good();
  const float &event_pTmiss();
  const bool &event_pass_tightMETFilters();
  const float &event_phimiss();
  const float &event_wgt();
  const float &event_wgt_L1PrefiringDn();
  const float &event_wgt_L1PrefiringUp();
  const float &event_wgt_SFs_PUJetId();
  const float &event_wgt_SFs_PUJetId_EffDn();
  const float &event_wgt_SFs_PUJetId_EffUp();
  const float &event_wgt_SFs_btagging();
  const float &event_wgt_SFs_btagging_EffDn();
  const float &event_wgt_SFs_btagging_EffUp();
  const float &event_wgt_SFs_electrons();
  const float &event_wgt_SFs_electrons_AltMCDn();
  const float &event_wgt_SFs_electrons_AltMCUp();
  const float &event_wgt_SFs_electrons_StatDn();
  const float &event_wgt_SFs_electrons_StatUp();
  const float &event_wgt_SFs_electrons_SystDn();
  const float &event_wgt_SFs_electrons_SystUp();
  const float &event_wgt_SFs_muons();
  const float &event_wgt_SFs_muons_AltMCDn();
  const float &event_wgt_SFs_muons_AltMCUp();
  const float &event_wgt_SFs_muons_StatDn();
  const float &event_wgt_SFs_muons_StatUp();
  const float &event_wgt_SFs_muons_SystDn();
  const float &event_wgt_SFs_muons_SystUp();
  const float &event_wgt_SFs_photons();
  const float &event_wgt_SFs_photons_EffDn();
  const float &event_wgt_SFs_photons_EffUp();
  const float &event_wgt_adjustment_NNPDF30();
  const float &event_wgt_triggers();
  const float &event_wgt_triggers_Dilepton();
  const float &event_wgt_triggers_SingleLepton();
  const float &event_wgt_triggers_SinglePhoton();
  const float &genmet_pTmiss();
  const float &genmet_phimiss();
  const float &lepton_eta();
  const int &lepton_id();
  const bool &lepton_is_genMatched_prompt();
  const float &lepton_mass();
  const float &lepton_phi();
  const float &lepton_pt();
  const float &lepton_relIso();
  const vector<float> &leptons_eff();
  const vector<float> &leptons_eff_StatDn();
  const vector<float> &leptons_eff_StatUp();
  const vector<float> &leptons_eff_SystDn();
  const vector<float> &leptons_eff_SystUp();
  const vector<float> &leptons_eta();
  const vector<int> &leptons_id();
  const vector<bool> &leptons_is_TOmatched_SingleLepton();
  const vector<bool> &leptons_is_genMatched_prompt();
  const vector<float> &leptons_mass();
  const vector<float> &leptons_phi();
  const vector<float> &leptons_pt();
  const float &min_abs_dPhi_pTj_pTmiss();
  const float &pConst_JJQCD_SIG_ghg2_1_JHUGen();
  const float &pConst_JJVBF_SIG_ghv1_1_JHUGen();
  const float &p_JJQCD_SIG_ghg2_1_JHUGen();
  const float &p_JJVBF_SIG_ghv1_1_JHUGen();
  const float &p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen();
  const float &p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen();
  const float &p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen();
  const float &p_JJVBF_SIG_ghv1prime2_1E4_JHUGen();
  const float &p_JJVBF_SIG_ghv2_1_JHUGen();
  const float &p_JJVBF_SIG_ghv4_1_JHUGen();
  const float &p_JJVBF_SIG_ghza1prime2_1E4_JHUGen();
  const float &photon_MIPTotalEnergy();
  const float &photon_eta();
  const float &photon_full5x5_r9();
  const float &photon_full5x5_sigmaIEtaIEta();
  const float &photon_full5x5_sigmaIPhiIPhi();
  const bool &photon_isEB();
  const bool &photon_isEBEEGap();
  const bool &photon_isEE();
  const bool &photon_isGap();
  const bool &photon_is_METSafe();
  const bool &photon_is_PFID();
  const bool &photon_is_beamHaloSafe();
  const bool &photon_is_conversionSafe();
  const bool &photon_is_genMatched_prompt();
  const bool &photon_is_inTime();
  const bool &photon_is_spikeSafe();
  const float &photon_mass();
  const bool &photon_pass_HGGSelection();
  const float &photon_phi();
  const float &photon_pt();
  const float &photon_seedTime();
  static void progress( int nEventsTotal, int nEventsChain );
};

#ifndef __CINT__
extern SkimTree st;
#endif

namespace tas {

  const unsigned long long &EventNumber();
  const unsigned int &LuminosityBlock();
  const unsigned int &RunNumber();
  const vector<float> &ak4jets_CEMF();
  const float &ak4jets_HT();
  const float &ak4jets_MHT();
  const vector<float> &ak4jets_NEMF();
  const vector<unsigned char> &ak4jets_btagWP_Bits();
  const vector<float> &ak4jets_eta();
  const vector<bool> &ak4jets_is_genMatched();
  const vector<bool> &ak4jets_is_genMatched_fullCone();
  const vector<float> &ak4jets_mass();
  const vector<float> &ak4jets_phi();
  const vector<float> &ak4jets_pt();
  const vector<float> &ak8jets_eta();
  const vector<float> &ak8jets_mass();
  const vector<float> &ak8jets_phi();
  const vector<float> &ak8jets_pt();
  const float &dPhi_pTboson_pTmiss();
  const float &dPhi_pTbosonjets_pTmiss();
  const float &dilepton_eta();
  const int &dilepton_id();
  const float &dilepton_mass();
  const float &dilepton_phi();
  const float &dilepton_pt();
  const float &electron_full5x5_r9();
  const float &electron_full5x5_sigmaIEtaIEta();
  const float &electron_full5x5_sigmaIPhiIPhi();
  const float &electron_seedTime();
  const vector<float> &electrons_full5x5_r9();
  const vector<float> &electrons_full5x5_sigmaIEtaIEta();
  const vector<float> &electrons_full5x5_sigmaIPhiIPhi();
  const vector<float> &electrons_seedTime();
  const float &event_mTZZ();
  const float &event_mZZ();
  const float &event_mllg();
  const unsigned int &event_n_ak4jets_pt20();
  const unsigned int &event_n_ak4jets_pt20_btagged_loose();
  const unsigned int &event_n_ak4jets_pt20_btagged_medium();
  const unsigned int &event_n_ak4jets_pt30();
  const unsigned int &event_n_ak4jets_pt30_btagged_loose();
  const unsigned int &event_n_ak4jets_pt30_btagged_medium();
  const unsigned int &event_n_vtxs_good();
  const float &event_pTmiss();
  const bool &event_pass_tightMETFilters();
  const float &event_phimiss();
  const float &event_wgt();
  const float &event_wgt_L1PrefiringDn();
  const float &event_wgt_L1PrefiringUp();
  const float &event_wgt_SFs_PUJetId();
  const float &event_wgt_SFs_PUJetId_EffDn();
  const float &event_wgt_SFs_PUJetId_EffUp();
  const float &event_wgt_SFs_btagging();
  const float &event_wgt_SFs_btagging_EffDn();
  const float &event_wgt_SFs_btagging_EffUp();
  const float &event_wgt_SFs_electrons();
  const float &event_wgt_SFs_electrons_AltMCDn();
  const float &event_wgt_SFs_electrons_AltMCUp();
  const float &event_wgt_SFs_electrons_StatDn();
  const float &event_wgt_SFs_electrons_StatUp();
  const float &event_wgt_SFs_electrons_SystDn();
  const float &event_wgt_SFs_electrons_SystUp();
  const float &event_wgt_SFs_muons();
  const float &event_wgt_SFs_muons_AltMCDn();
  const float &event_wgt_SFs_muons_AltMCUp();
  const float &event_wgt_SFs_muons_StatDn();
  const float &event_wgt_SFs_muons_StatUp();
  const float &event_wgt_SFs_muons_SystDn();
  const float &event_wgt_SFs_muons_SystUp();
  const float &event_wgt_SFs_photons();
  const float &event_wgt_SFs_photons_EffDn();
  const float &event_wgt_SFs_photons_EffUp();
  const float &event_wgt_adjustment_NNPDF30();
  const float &event_wgt_triggers();
  const float &event_wgt_triggers_Dilepton();
  const float &event_wgt_triggers_SingleLepton();
  const float &event_wgt_triggers_SinglePhoton();
  const float &genmet_pTmiss();
  const float &genmet_phimiss();
  const float &lepton_eta();
  const int &lepton_id();
  const bool &lepton_is_genMatched_prompt();
  const float &lepton_mass();
  const float &lepton_phi();
  const float &lepton_pt();
  const float &lepton_relIso();
  const vector<float> &leptons_eff();
  const vector<float> &leptons_eff_StatDn();
  const vector<float> &leptons_eff_StatUp();
  const vector<float> &leptons_eff_SystDn();
  const vector<float> &leptons_eff_SystUp();
  const vector<float> &leptons_eta();
  const vector<int> &leptons_id();
  const vector<bool> &leptons_is_TOmatched_SingleLepton();
  const vector<bool> &leptons_is_genMatched_prompt();
  const vector<float> &leptons_mass();
  const vector<float> &leptons_phi();
  const vector<float> &leptons_pt();
  const float &min_abs_dPhi_pTj_pTmiss();
  const float &pConst_JJQCD_SIG_ghg2_1_JHUGen();
  const float &pConst_JJVBF_SIG_ghv1_1_JHUGen();
  const float &p_JJQCD_SIG_ghg2_1_JHUGen();
  const float &p_JJVBF_SIG_ghv1_1_JHUGen();
  const float &p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen();
  const float &p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen();
  const float &p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen();
  const float &p_JJVBF_SIG_ghv1prime2_1E4_JHUGen();
  const float &p_JJVBF_SIG_ghv2_1_JHUGen();
  const float &p_JJVBF_SIG_ghv4_1_JHUGen();
  const float &p_JJVBF_SIG_ghza1prime2_1E4_JHUGen();
  const float &photon_MIPTotalEnergy();
  const float &photon_eta();
  const float &photon_full5x5_r9();
  const float &photon_full5x5_sigmaIEtaIEta();
  const float &photon_full5x5_sigmaIPhiIPhi();
  const bool &photon_isEB();
  const bool &photon_isEBEEGap();
  const bool &photon_isEE();
  const bool &photon_isGap();
  const bool &photon_is_METSafe();
  const bool &photon_is_PFID();
  const bool &photon_is_beamHaloSafe();
  const bool &photon_is_conversionSafe();
  const bool &photon_is_genMatched_prompt();
  const bool &photon_is_inTime();
  const bool &photon_is_spikeSafe();
  const float &photon_mass();
  const bool &photon_pass_HGGSelection();
  const float &photon_phi();
  const float &photon_pt();
  const float &photon_seedTime();
}
#endif
