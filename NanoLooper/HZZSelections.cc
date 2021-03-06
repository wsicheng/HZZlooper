#include "HZZSelections.h"
// #include "Utilities.h"
#include "config/prescales/photon_prescales_2018.hh"
// #include "config/prescales/photon_prescales_run2.hh"


inline bool isCloseObject(const float eta1, const float phi1, const float eta2, const float phi2, const float conesize, float* deltaR = nullptr) {
  const float PI = TMath::Pi();
  float deltaEta = fabs(eta1 - eta2);
  if (deltaEta > conesize) return false;
  float deltaPhi = fabs(phi1 - phi2);
  if (deltaPhi > PI) deltaPhi = 2*PI - deltaPhi;
  if (deltaPhi > conesize) return false;
  float deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
  if (deltaR2 > conesize*conesize) return false;
  if (deltaR) *deltaR = sqrt(deltaR2);

  return true;
}

using namespace std;
using namespace tas;

RoccoR* muoncorr = nullptr;
TRandom3* randomGenerator = nullptr;

bool passTriggerSelections(int trigtype) {
  if (!gconf.is_data) return true;  // not using the trigger emulatino in MC

  if (trigtype == 1) {
    // Muon triggers
    switch (gconf.year) {
      case 2016:
        return ( HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL() || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ() ||
                 HLT_IsoMu24() || HLT_IsoTkMu24());
      case 2017:
        return ( HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ() ||
                 ((run() >= 299337)? HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8() : false) ||
                 HLT_IsoMu24() || HLT_IsoMu27());
      case 2018:
      default:
        return ( HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8() ||
                 HLT_IsoMu24());
    }
  } else if (trigtype == 2) {
    // Electron triggers
    switch (gconf.year) {
      case 2016:
        return ( HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() ||
                 HLT_Ele27_WPTight_Gsf() );
      case 2017:
      case 2018:
      default:
        return ( HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() ||
                 ((run() >= 302026)? HLT_Ele32_WPTight_Gsf() : HLT_Ele35_WPTight_Gsf()) ||
                 HLT_Ele32_WPTight_Gsf_L1DoubleEG() );
    }
  } else if (trigtype == 3) {
    // emu triggers
    switch (gconf.year) {
      case 2016:
      case 2017:
      case 2018:
      default:
        return (
            HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL() || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() ||
            HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL() || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ() ||
            HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL()  || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ()  );
    }
  } else if (trigtype == 4) {
    // Photon triggers
    bool prescaled = false;
    switch (gconf.year) {
      case 2016:
        // prescaled photon triggers
        if (false) {
          prescaled = (
              HLT_Photon22_R9Id90_HE10_IsoM() ||
              HLT_Photon30_R9Id90_HE10_IsoM() ||
              HLT_Photon36_R9Id90_HE10_IsoM() ||
              HLT_Photon50_R9Id90_HE10_IsoM() ||
              HLT_Photon75_R9Id90_HE10_IsoM() ||
              HLT_Photon90_R9Id90_HE10_IsoM() ||
              HLT_Photon120_R9Id90_HE10_IsoM() ||
              HLT_Photon165_R9Id90_HE10_IsoM()
              );
        }

        // unprescaled in 2016, prescaled in 2017 & 2018
        return ( prescaled || HLT_Photon165_R9Id90_HE10_IsoM() || HLT_Photon250_NoHE() || HLT_Photon300_NoHE());
        // HLT_Photon175() || HLT_Photon200()  // unprescaled, but not use
        // HLT_Photon250_NoHE(); // in 2016 only that doesn't exist
      case 2017:
      case 2018:
      default:
        if (false) {
          prescaled = (
              HLT_Photon50_R9Id90_HE10_IsoM() ||  // eff lumi: 0.
              // HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50() ||
              HLT_Photon75_R9Id90_HE10_IsoM() ||
              // HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3() ||
              // HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3() ||
              // HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3() ||
              // HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3() ||
              HLT_Photon90_R9Id90_HE10_IsoM() ||
              HLT_Photon120_R9Id90_HE10_IsoM() ||
              HLT_Photon165_R9Id90_HE10_IsoM()
              );
        }

        return ( prescaled || HLT_Photon200() );
    }

    // Lepton triggers that doesn't exist
    // HLT_IsoTkMu24() ||
    // HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL() ||
    // HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ() ||
  }

  return false;
}

double getPhotonTrigPrescale(double objpt) {
  // precales are taken from http://homepage.iihe.ac.be/~mmahdavi/Analysis/trigsLums/trigsLumis_2016
  // Take 10% more, so that we are on the plateau (the pT in the name is the one at the middle of the turn-on curve, so at 50% efficiency).
  if (!gconf.is_data) return 1.0;

  bool useSinglePrescale = true;
  const std::map<unsigned,std::map<unsigned,int>>* runMap = nullptr;

  switch (gconf.year) {
    case 2016:
      if ( HLT_Photon250_NoHE() || HLT_Photon300_NoHE() )
        return 1.0;

      if (!useSinglePrescale) {
        if ( objpt > 181.5 ) {
          if (HLT_Photon165_R9Id90_HE10_IsoM())
            return 1.0;
          else
            return 0.0;
        } else if ( objpt > 132. ) {
          if (HLT_Photon120_R9Id90_HE10_IsoM())
            runMap = &tab_prescales_HLT_Photon120_R9Id90_HE10_IsoM;
          else
            return 0.0; // 1./(1.-0.59633)
        } else if ( objpt > 99. ) {
          if (HLT_Photon90_R9Id90_HE10_IsoM())
            runMap = &tab_prescales_HLT_Photon90_R9Id90_HE10_IsoM;
          else
            return 0.0; // 1./(1.-0.85608)
        } else if ( objpt > 82.5 ) {
          if (HLT_Photon75_R9Id90_HE10_IsoM())
            runMap = &tab_prescales_HLT_Photon75_R9Id90_HE10_IsoM;
          else
            return 0.0; // 1./(1.-0.92846)
        } else if ( objpt > 55. ) {
          if (HLT_Photon50_R9Id90_HE10_IsoM())
            runMap = &tab_prescales_HLT_Photon50_R9Id90_HE10_IsoM;
          else
            return 0.0; // 1./(1.-0.98606)
        } else {
          return 0.0;
        }
      }
      else if (useSinglePrescale) {
        if ( objpt > 190 )
          return (HLT_Photon175())? 1.0 : 0.0;
        else if ( objpt > 181.5 )
          return (HLT_Photon165_R9Id90_HE10_IsoM())? 1.0 : 0.0;
        else if ( objpt > 132. )
          return (HLT_Photon120_R9Id90_HE10_IsoM())? 2.4773 : 0.0; // 1./(1.-0.59633)
        else if ( objpt > 99. )
          return (HLT_Photon90_R9Id90_HE10_IsoM())? 6.9483 : 0.0; // 1./(1.-0.85608)
        else if ( objpt > 82.5 )
          return (HLT_Photon75_R9Id90_HE10_IsoM())? 13.978 : 0.0; // 1./(1.-0.92846)
        else if ( objpt > 55. )
          return (HLT_Photon50_R9Id90_HE10_IsoM())? 71.736 : 0.0; // 1./(1.-0.98606)
        else if ( objpt > 39.3 )
          return (HLT_Photon36_R9Id90_HE10_IsoM())? 164.47 : 0.0; // 1./(1.-0.99392)
        else if ( objpt > 33. )
          return (HLT_Photon30_R9Id90_HE10_IsoM())? 367.65 : 0.0; // 1./(1.-0.99728)
        else if ( objpt > 24.2 )
          return (HLT_Photon22_R9Id90_HE10_IsoM())? 1851.9 : 0.0; // 1./(1.-0.99946)
      }
      break;

    case 2017:
    case 2018:
      if (!useSinglePrescale) {
        if ( objpt > 230 ) {
          return (HLT_Photon200())? 1.0 : 0.0;
        } else if ( objpt > 181.5 ) {
          if (HLT_Photon165_R9Id90_HE10_IsoM())
            runMap = &tab_prescales_HLT_Photon165_R9Id90_HE10_IsoM;
          else
            return 0.0;
        } else if ( objpt > 132. ) {
          if (HLT_Photon120_R9Id90_HE10_IsoM())
            runMap = &tab_prescales_HLT_Photon120_R9Id90_HE10_IsoM;
          else
            return 0.0;
        } else if ( objpt > 99. ) {
          if (HLT_Photon90_R9Id90_HE10_IsoM())
            runMap = &tab_prescales_HLT_Photon90_R9Id90_HE10_IsoM;
          else
            return 0.0;
        } else if ( objpt > 82.5 ) {
          if (HLT_Photon75_R9Id90_HE10_IsoM())
            runMap = &tab_prescales_HLT_Photon75_R9Id90_HE10_IsoM;
          else
            return 0.0;
        } else if ( objpt > 55. ) {
          if (HLT_Photon50_R9Id90_HE10_IsoM())
            runMap = &tab_prescales_HLT_Photon50_R9Id90_HE10_IsoM;
          else
            return 0.0;
        } else {
          return 0.0;
        }
      }
    default:
      break;
  }

  if (runMap) {
    // For run number, every run shall exisit in the list
    map<unsigned,map<unsigned,int>>::const_iterator lumiMap = runMap->find(run());
    if (lumiMap == runMap->end()) {
      cout << "[getPhotonTrigPrescale] >> Cannot find run " << run() << " in the prescale table!" << endl;
      return 1.0;
    }
    // For lumi number, only the lowest lumi of each prescale is recorded
    auto ilumi = (lumiMap->second).upper_bound(luminosityBlock());
    if (ilumi != (lumiMap->second).begin()) ilumi--;
    if (luminosityBlock() < ilumi->first) {
      cout << "[getPhotonTrigPrescale] >> Cannot find luminosity block " << luminosityBlock() << " in run " << run()
           << " in the prescale table! Reverting to old method!" << endl;
      return 1.0;
    }
    int prescale = ilumi->second;

    return double(prescale);
  }

  return 0.0;
}


float getDileptonMT(const TLorentzVector& boson, const TLorentzVector& metvec) {
  double ptll = boson.Pt(), mll = boson.M(), met = metvec.Pt();
  double term1 = sqrt(ptll*ptll + mll*mll) + sqrt(met*met + mZ*mZ);
  double term2 = (boson + metvec).Pt();
  return sqrt( term1*term1 - term2*term2 );
}

void giveMassToPhoton(TLorentzVector& boson, TH1* h_weight) {
  double mass = mZ;
  if (h_weight != nullptr) {
    do {
      mass = h_weight->GetRandom();
    } while (fabs(mass-mZ) > 15);
  }
  boson.SetPtEtaPhiE(boson.Pt(), boson.Eta(), boson.Phi(), sqrt(pow(mass,2)+pow(boson.P(),2)) );
}


std::tuple<vector<Electron>, vector<Electron>> getElectrons() {

  vector<Electron> looseElectrons;
  vector<Electron> tightElectrons;

  for (unsigned i = 0; i < Electron_pt().size(); ++i) {
    double const etaSc = Electron_deltaEtaSC()[i] + Electron_eta()[i];
    double const absEtaSc = std::abs(etaSc);
    bool const passLooseId = Electron_mvaFall17V2noIso_WP90()[i];
    bool const passLooseIso = (Electron_pfRelIso03_all()[i] < 0.1);

    if (Electron_pt()[i] < k_minPt_el_loose or absEtaSc > 2.5 or not passLooseId or not passLooseIso)
      continue;

    Electron electron;
    electron.p4.SetPtEtaPhiM(Electron_pt()[i], Electron_eta()[i], Electron_phi()[i], Electron_mass()[i]);
    electron.charge = Electron_charge()[i];
    electron.etaSc = etaSc;

    // if (IsDuplicate(electron.p4, 0.1)) continue;  // for jet cleaning

    looseElectrons.emplace_back(electron);

    bool const passTightId = Electron_mvaFall17V2noIso_WP90()[i];
    bool const passTightIso = (Electron_pfRelIso03_all()[i] < 0.1);

    if (Electron_pt()[i] < k_minPt_lep_tight or not passTightId or not passTightIso)
      continue;

    if (absEtaSc > 1.4442 and absEtaSc < 1.5660)  // EB-EE gap
      continue;

    tightElectrons.emplace_back(electron);
  }

  std::sort(looseElectrons.begin(), looseElectrons.end(), PtOrdered);
  std::sort(tightElectrons.begin(), tightElectrons.end(), PtOrdered);
  return {tightElectrons, looseElectrons};
}


void ApplyRochesterCorrectionToMuon(Muon *muon, int idx) {

  // Apply the correction only in its domain of validity
  if (muon->p4.Pt() > 200. or std::abs(muon->p4.Eta()) > 2.4)
    return;

  double scale = 1.;

  if (gconf.is_data) {
    scale = muoncorr->kScaleDT(Muon_charge()[idx], Muon_pt()[idx], Muon_eta()[idx], Muon_phi()[idx]);
  } else {
    // Flavour of genParticle for MC matching to status==1 muons: 1 = prompt muon (including
    // gamma*->mu mu), 15 = muon from prompt tau, 5 = muon from b, 4 = muon from c, 3 = muon from
    // light or unknown, 0 = unmatched
    if (Muon_genPartFlav().at(idx) > 0) {
      scale = muoncorr->kSpreadMC(Muon_charge()[idx], Muon_pt()[idx], Muon_eta()[idx], Muon_phi()[idx],
                                  GenPart_pt().at(Muon_genPartIdx()[idx]));
    } else {
      randomGenerator->SetSeed(int(100000 * (Muon_phi()[idx] + 3.2)));
      scale = muoncorr->kSmearMC(Muon_charge()[idx], Muon_pt()[idx], Muon_eta()[idx], Muon_phi()[idx],
                                 Muon_nTrackerLayers()[idx], randomGenerator->Uniform());
    }
  }

  muon->p4.SetPtEtaPhiM(Muon_pt()[idx] * scale, Muon_eta()[idx], Muon_phi()[idx], Muon_mass()[idx]);  // Q: no scale for mass?
}

std::tuple<vector<Muon>, vector<Muon>> getMuons(bool applyRocCorr, float* shiftx, float* shifty) {

  vector<Muon> looseMuons;
  vector<Muon> tightMuons;

  for (unsigned i = 0; i < Muon_pt().size(); ++i) {
    // Loose ID as per https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Loose_Muon
    bool const passLooseId = Muon_isPFcand()[i] && Muon_isGlobal()[i] && Muon_isTracker()[i];
    bool const passLooseIso = (Muon_pfRelIso03_all()[i] < 0.25);

    if (std::abs(Muon_eta()[i]) > 2.4 or not passLooseId or not passLooseIso)
      continue;

    Muon muon;
    muon.uncorrP4.SetPtEtaPhiM(Muon_pt()[i], Muon_eta()[i], Muon_phi()[i], Muon_mass()[i]);
    muon.charge = Muon_charge()[i];

    if (applyRocCorr)
      ApplyRochesterCorrectionToMuon(&muon, i);
    else
      muon.p4 = muon.uncorrP4;

    if (muon.p4.Pt() < k_minPt_mu_loose)  // minPtLoose
      continue;

    looseMuons.emplace_back(muon);

    // Propagate changes in momenta of loose muons into ptmiss
    if (applyRocCorr) {
      if (shiftx) *shiftx = muon.p4.Px() - muon.uncorrP4.Px();
      if (shifty) *shifty = muon.p4.Py() - muon.uncorrP4.Py();
    }

    bool const passTightId = Muon_mediumPromptId()[i];
    bool const passTightIso = (Muon_pfRelIso03_all()[i] < 0.15);

    if (muon.p4.Pt() < k_minPt_lep_tight or not passTightId or not passTightIso) // minPtTight
      continue;

    tightMuons.emplace_back(muon);
  }

  std::sort(looseMuons.begin(), looseMuons.end(), PtOrdered);
  std::sort(tightMuons.begin(), tightMuons.end(), PtOrdered);
  return {tightMuons, looseMuons};

}

vector<int> getIsoTrackIndices(const vector<Muon>& mus, const vector<Electron>& els) {
  vector<int> trackIdxs;
  for (unsigned i = 0; i < nIsoTrack(); ++i) {
    if (!IsoTrack_isPFcand()[i]) continue;  // only consider pfcandidates
    int id = abs(IsoTrack_pdgId()[i]);
    if (id != 11 && id != 13 && id < 100) continue;
    float pt = IsoTrack_pt()[i];
    if (pt < ((id < 14)? 5 : 10)) continue;
    if (fabs(IsoTrack_eta()[i]) > ((id == 11)? 3.0 : 2.4) ) continue;
    if (fabs(IsoTrack_dz()[i]) > 0.1) continue;
    if (pt > 60) {
      if (IsoTrack_pfRelIso03_chg()[i] * pt > 6) continue;
    } else {
      if (IsoTrack_pfRelIso03_chg()[i] > 0.1) continue;
    }

    // Perform angular cleaning w.r.t. recognized leptons
    for (auto lep : mus) {
      if (isCloseObject(IsoTrack_eta()[i], IsoTrack_phi()[i], lep.p4.Eta(), lep.p4.Phi(), 0.4))
        goto end_of_loop_tracks;
    }
    for (auto lep : els) {
      if (isCloseObject(IsoTrack_eta()[i], IsoTrack_phi()[i], lep.p4.Eta(), lep.p4.Phi(), 0.4))
        goto end_of_loop_tracks;
    }

    trackIdxs.push_back(i);

 end_of_loop_tracks:;
  }

  return trackIdxs;
}

vector<Photon> getPhotons(const vector<Muon>& mus, const vector<Electron>& els) {
  vector<Photon> photons;

  for (unsigned i = 0; i < Photon_pt().size(); ++i) {
    // Tight ID
    // const bool passId = (gconf.year == 2016)? (Photon_cutBased()[i] >= 3) : (Photon_cutBasedBitmap()[i] & 0b0100);
    const bool passId = (Photon_cutBased()[i] >= 3);  // new format for 2018 data...

    if (Photon_pt()[i] < k_minPt_photon or not passId)
      continue;

    // Only consider photons in the barrel
    if (!Photon_isScEtaEB()[i])
      continue;

    Photon gamma;
    // Perform angular cleaning
    // if (IsDuplicate(gamma.p4, 0.1)) continue; // <-- to be revisit
    for (auto lep : mus) {
      if (isCloseObject(Photon_eta()[i], Photon_phi()[i], lep.p4.Eta(), lep.p4.Phi(), 0.4))
        goto end_of_loop_photons;
    }
    for (auto lep : els) {
      if (isCloseObject(Photon_eta()[i], Photon_phi()[i], lep.p4.Eta(), lep.p4.Phi(), 0.4))
        goto end_of_loop_photons;
    }

    gamma.p4.SetPtEtaPhiM(Photon_pt()[i], Photon_eta()[i], Photon_phi()[i], 0.);
    photons.emplace_back(gamma);

 end_of_loop_photons:;
  }

  // Make sure the collection is ordered in pt
  std::sort(photons.begin(), photons.end(), PtOrdered);

  return photons;
}


float getJetCorrectionFactorFromFile(int jetidx, Jet injet, bool applyJER) {

    double corrFactor = 1.;
    // const double corrPt = injet.p4.Pt();  // pt with nominal JEC

    /*
    // Evaluate JEC uncertainty
    if (syst_ == Syst::JEC) {
      jecUncProvider_->setJetEta(jet.p4.Eta());
      jecUncProvider_->setJetPt(corrPt);
      const double uncertainty = jecUncProvider_->getUncertainty(true);

      if (systDirection_ == SystDirection::Up)
        corrFactor *= (1. + uncertainty);
      else
        corrFactor *= (1. - uncertainty);
    }

    // Apply JER smearing. Corresponding correction factor is always evaluated
    // with nominal JEC applied, even if a JEC variation has been requested.
    // This aligns with how the JER smearing is usually applied in CMSSW.
    if (isSim_) {
      // Relative jet pt resolution in simulation
      double const ptResolution = jerProvider_->getResolution(
        {{JME::Binning::JetPt, corrPt}, {JME::Binning::JetEta, jet.p4.Eta()},
         {JME::Binning::Rho, *puRho_}});

      // Find data-to-simulation scale factor
      Variation jerDirection;

      if (syst_ == Syst::JER) {
        if (systDirection_ == SystDirection::Up)
          jerDirection = Variation::UP;
        else
          jerDirection = Variation::DOWN;
      } else
        jerDirection = Variation::NOMINAL;

      double const jerSF = jerSFProvider_->getScaleFactor(
        {{JME::Binning::JetEta, jet.p4.Eta()}}, jerDirection);


      // Depending on the presence of a matching generator-level jet, perform
      // deterministic or stochastic smearing [1]
      // [1] https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution?rev=71#Smearing_procedures
      GenJet const *genJet = FindGenMatch(jet, ptResolution);


      if (genJet) {
        double const jerFactor = 1. +
          (jerSF - 1.) * (corrPt - genJet->p4.Pt()) / corrPt;
        corrFactor *= jerFactor;
      } else {
        double const jerFactor = 1. + randomGenerator_.Gaus(0., ptResolution) *
          std::sqrt(std::max(std::pow(jerSF, 2) - 1., 0.));
        corrFactor *= jerFactor;
      }
    }

    TLorentzVector const originalP4 = jet.p4;
    jet.p4 *= corrFactor;

    // Propagate the change in jet momentum for the use in ptmiss. The type 1
    // correction to ptmiss has been applied using jets with originalP4.Pt()
    // above 15 GeV. After the additional corrections applied above, the set of
    // jets with pt > 15 GeV has changed. However, without the access to raw
    // momenta of jets, it is not possible to undo the contribution to the
    // type 1 correction from jets whose pt has changed from above to below
    // 15 GeV. For simplicity, use the same set of jets as in the original
    // type 1 correction (modulus the different jet ID).

    // if (originalP4.Pt() > 15.)
    //   AddMomentumShift(originalP4, jet.p4);
    */
    return corrFactor;
}

bool PassVBFcuts(const vector<Jet> &selJets, const TLorentzVector &boson) {
  if (selJets.size() < 2) return false;

  float etamax = selJets[0].p4.Eta();
  float etamin = selJets[1].p4.Eta();
  if (etamax < etamin) std::swap(etamin, etamax);

  bool centralJetVeto = false;
  if (selJets.size() > 2) {
    for (unsigned int i = 2; i < selJets.size(); i++) {
      if (selJets[i].p4.Eta() > etamin && selJets[i].p4.Eta() < etamax)
        centralJetVeto = true;
    }
  }
  bool centralBoson = (boson.Eta() > etamin && boson.Eta() < etamax);
  bool passDeltaEta = ((etamax - etamin) > 4.);
  bool passMjj = ((selJets[0].p4 + selJets[1].p4).M() > 500);

  if (!centralJetVeto && centralBoson && passDeltaEta && passMjj) return true;

  return false;
}

vector<Jet> getJets(const vector<Muon>& mus, const vector<Electron>& els, const vector<Photon>& phs, bool reapplyJEC, bool isSim) {
  vector<Jet> jets;

  for (unsigned i = 0; i < Jet_pt().size(); ++i) {
    if (not (Jet_jetId()[i] & 0b0010))  // bit1 always false in 2017 since it does not exist
      continue;

    Jet jet;
    jet.p4.SetPtEtaPhiM(Jet_pt()[i], Jet_eta()[i], Jet_phi()[i], Jet_mass()[i]);

    // Perform angular cleaning w.r.t. recognized leptons and photons
    for (auto lep : mus) {
      if (isCloseObject(jet.p4.Eta(), jet.p4.Phi(), lep.p4.Eta(), lep.p4.Phi(), 0.4))
        goto end_of_loop_jets;
    }
    for (auto lep : els) {
      if (isCloseObject(jet.p4.Eta(), jet.p4.Phi(), lep.p4.Eta(), lep.p4.Phi(), 0.4))
        goto end_of_loop_jets;
    }
    for (auto photon : phs) {
      if (isCloseObject(jet.p4.Eta(), jet.p4.Phi(), photon.p4.Eta(), photon.p4.Phi(), 0.4))
        goto end_of_loop_jets;
    }

    jet.bTag = Jet_btagDeepFlavB()[i];
    if (!gconf.is_data)
      jet.hadronFlavour = Jet_hadronFlavour()[i];
    else
      jet.hadronFlavour = 0;

    // if (reapplyJEC)
    //   jet.p4 *= getJetCorrectionFactorFromFile();

    // Kinematical cuts for jets to be stored in the collection
    if (jet.p4.Pt() < k_minPt_jet or std::abs(jet.p4.Eta()) > k_maxAbsEta_jet)
      continue;

    jets.emplace_back(jet);

 end_of_loop_jets:;
  }

  // Make sure jets are sorted in pt
  std::sort(jets.begin(), jets.end(), PtOrdered);

  return jets;
}

bool ZZ2l2vPrunerCuts() {
  double const minLepPt = 20.;
  std::vector<LorentzVector> leptonMomenta;

  for (unsigned i = 0; i < nElectron(); ++i) {
    if (Electron_pt()[i] < minLepPt or Electron_cutBased()[i] < 4  /* fails tight working point */)
      continue;

    leptonMomenta.emplace_back(Electron_pt()[i], Electron_eta()[i], Electron_phi()[i], Electron_mass()[i]);
  }

  for (unsigned i = 0; i < nMuon(); ++i) {
    if (Muon_pt()[i] < 0.9 * minLepPt or not Muon_tightId()[i])
      continue;

    leptonMomenta.emplace_back(Muon_pt()[i], Muon_eta()[i], Muon_phi()[i], Muon_mass()[i]);
  }

  if (leptonMomenta.size() < 2)
    return false;

  for (int i = 0; i < int(leptonMomenta.size()) - 1; ++i) {
    for (int j = i + 1; j < int(leptonMomenta.size()); ++j) {
      LorentzVector const p4Z = leptonMomenta[i] + leptonMomenta[j];
      if (p4Z.M() > 50 and p4Z.Pt() > 50)
        return true;
    }
  }

  return false;
}

bool InstrMETPrunerCuts() {
  // for (unsigned i = 0; i < nPhoton(); ++i) {
  for (unsigned i = 0; i < Photon_pt().size(); ++i) {
    // Tight ID
    if (Photon_pt()[i] < 50) continue;

    // const bool passId = (gconf.year == 2016)? (Photon_cutBased()[i] >= 3) : (Photon_cutBasedBitmap()[i] & 0b0100);
    const bool passId = (Photon_cutBased()[i] >= 3);

    if (!passId) continue;
    // Only consider photons in the barrel
    // if (photons_eta()[i] > 2.5) continue;
    return true;
  }

  return false;
}
