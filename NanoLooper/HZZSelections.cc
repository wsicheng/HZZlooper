#include "HZZSelections.h"
// #include "Utilities.h"

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
              HLT_Photon50_R9Id90_HE10_IsoM() ||  // eff lumi: 0.
              HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50() ||
              HLT_Photon75_R9Id90_HE10_IsoM() ||
              HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3() ||
              HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3() ||
              HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3() ||
              HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3() ||
              HLT_Photon90_R9Id90_HE10_IsoM() ||
              HLT_Photon120_R9Id90_HE10_IsoM() );
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
              HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50() ||
              HLT_Photon75_R9Id90_HE10_IsoM() ||
              HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3() ||
              HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3() ||
              HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3() ||
              HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3() ||
              HLT_Photon90_R9Id90_HE10_IsoM() ||
              HLT_Photon120_R9Id90_HE10_IsoM() );
        }

        // unprescaled in 2016, prescaled in 2017 & 2018
        return ( prescaled || HLT_Photon300_NoHE() );
        // HLT_Photon175() || HLT_Photon200()  // unprescaled, but not use
    }
    // Lepton triggers that doesn't exist
    // HLT_IsoTkMu24() ||
    // HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL() ||
    // HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ() ||
  }

  return false;
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
    bool const passLooseId = (Electron_cutBased()[i] >= 2);

    if (Electron_pt()[i] < k_minPt_el_loose or absEtaSc > 2.5 or not passLooseId)
      continue;

    Electron electron;
    electron.p4.SetPtEtaPhiM(Electron_pt()[i], Electron_eta()[i], Electron_phi()[i], Electron_mass()[i]);
    electron.charge = Electron_charge()[i];
    electron.etaSc = etaSc;

    // if (IsDuplicate(electron.p4, 0.1)) continue;  // for jet cleaning

    looseElectrons.emplace_back(electron);

    bool const passTightId = Electron_cutBased()[i] >= 4;

    if (Electron_pt()[i] < k_minPt_lep_tight or not passTightId)
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

    if (std::abs(Muon_eta()[i]) > 2.4 or not passLooseId or Muon_pfRelIso04_all()[i] > 0.25)
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

    bool const passTightId = Muon_tightId()[i];

    if (muon.p4.Pt() < k_minPt_lep_tight or not passTightId or Muon_pfRelIso04_all()[i] > 0.15) // minPtTight
      continue;

    tightMuons.emplace_back(muon);
  }

  std::sort(looseMuons.begin(), looseMuons.end(), PtOrdered);
  std::sort(tightMuons.begin(), tightMuons.end(), PtOrdered);
  return {tightMuons, looseMuons};

}

vector<Photon> getPhotons() {
  vector<Photon> photons;

  for (unsigned i = 0; i < Photon_pt().size(); ++i) {
    // Tight ID
    const bool passId = (gconf.year == 2016)? (Photon_cutBased()[i] >= 3) : (Photon_cutBasedBitmap()[i] & 0b0100);

    if (Photon_pt()[i] < k_minPt_photon or not passId)
      continue;

    // Only consider photons in the barrel
    if (!Photon_isScEtaEB()[i])
      continue;

    Photon gamma;
    gamma.p4.SetPtEtaPhiM(Photon_pt()[i], Photon_eta()[i], Photon_phi()[i], 0.);

    // Perform angular cleaning
    // if (IsDuplicate(gamma.p4, 0.1)) continue; // <-- to be revisit

    photons.emplace_back(gamma);
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
