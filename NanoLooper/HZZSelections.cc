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

void giveMassToPhoton(TLorentzVector& boson, TH1* h_weight) {
  double mass = MZ;
  if (h_weight != nullptr) {
    do {
      mass = h_weight->GetRandom();
    } while (fabs(mass-MZ) > 15);
  }

  boson.SetPtEtaPhiE(boson.Pt(), boson.Eta(), boson.Phi(), sqrt(pow(mass,2)+pow(boson.P(),2)) );
}

vector<Electron> getElectrons(int id_level) {

  vector<Electron> looseElectrons;
  vector<Electron> tightElectrons;

  for (unsigned i = 0; i < Electron_pt().size(); ++i) {
    double const etaSc = Electron_deltaEtaSC()[i] + Electron_eta()[i];
    double const absEtaSc = std::abs(etaSc);
    bool const passLooseId = (Electron_cutBased()[i] >= 2);

    if (Electron_pt()[i] < k_minPt_lep_loose or absEtaSc > 2.5 or not passLooseId)
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

  // Make sure the collections are ordered in pt
  switch (id_level) {
    case idLoose:
      std::sort(looseElectrons.begin(), looseElectrons.end(), PtOrdered);
      return tightElectrons;
    case idTight:
      std::sort(tightElectrons.begin(), tightElectrons.end(), PtOrdered);
      return tightElectrons;
  }
  return tightElectrons;
}


std::optional<GenParticle> FindGenMatchMuon(Muon const &muon, double maxDR) {

  unsigned iClosest = -1;
  float minDR = maxDR;

  for (unsigned i = 0; i < GenPart_pdgId().size(); ++i) {
    if (std::abs(GenPart_pdgId().at(i)) != 13)
      // Only consider muons
      continue;
    if (isCloseObject(muon.p4.Eta(), muon.p4.Phi(), GenPart_eta().at(i), GenPart_phi().at(i), minDR, &minDR)) {
      iClosest = i;
    }
  }

  if (iClosest != unsigned(-1)) {
    GenParticle matchedParticle{GenPart_pdgId().at(iClosest)};
    matchedParticle.p4.SetPtEtaPhiM(
        GenPart_pt().at(iClosest), GenPart_eta().at(iClosest), GenPart_phi().at(iClosest),
        0.1057);

    return matchedParticle;
  } else {
    return {};
  }
} 

void ApplyRochesterCorrectionToMuon(Muon *muon, int trackerLayers, bool isSim) {

  // Apply the correction only in its domain of validity
  if (muon->p4.Pt() > 200. or std::abs(muon->p4.Eta()) > 2.4)
    return;

  double scaleFactor = 1.;

  // if (isSim) {
  //   auto const genMatch = FindGenMatchMuon(*muon, 0.01);
  //   if (genMatch)
  //     scaleFactor = rochesterCorrection_->kScaleFromGenMC(
  //       muon->charge, muon->p4.Pt(), muon->p4.Eta(), muon->p4.Phi(),
  //       trackerLayers, genMatch->p4.Pt(), randomGenerator_.Uniform());
  //   else
  //     scaleFactor = rochesterCorrection_->kScaleAndSmearMC(
  //       muon->charge, muon->p4.Pt(), muon->p4.Eta(), muon->p4.Phi(),
  //       trackerLayers, randomGenerator_.Uniform(), randomGenerator_.Uniform());
  // } else {
  //   scaleFactor = rochesterCorrection_->kScaleDT(
  //     muon->charge, muon->p4.Pt(), muon->p4.Eta(), muon->p4.Phi());
  // }
  
  muon->p4.SetPtEtaPhiM(muon->p4.Pt() * scaleFactor, muon->p4.Eta(),
                        muon->p4.Phi(), muon->p4.M());  // Q: no scale for mass?
}

vector<Muon> getMuons(int id_level) {

  vector<Muon> looseMuons;
  vector<Muon> tightMuons;

  for (unsigned i = 0; i < Muon_pt().size(); ++i) {
    // Loose ID as per https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Loose_Muon
    bool const passLooseId = Muon_isPFcand()[i] && Muon_isGlobal()[i] && Muon_isTracker()[i];

    if (std::abs(Muon_eta()[i]) > 2.4 or not passLooseId or Muon_pfRelIso04_all()[i] > 0.25)
      continue;

    Muon muon;
    muon.p4.SetPtEtaPhiM(Muon_pt()[i], Muon_eta()[i], Muon_phi()[i], Muon_mass()[i]);
    muon.uncorrP4 = muon.p4;
    muon.charge = Muon_charge()[i];

    // ApplyRochesterCorrectionToMuon(&muon, Muon_nTrackerLayers()[i]);

    if (muon.p4.Pt() < k_minPt_lep_loose)  // minPtLoose
      continue;

    looseMuons.emplace_back(muon);

    // Propagate changes in momenta of loose muons into ptmiss
    // AddMomentumShift(muon.uncorrP4, muon.p4);

    bool const passTightId = Muon_tightId()[i];

    if (muon.p4.Pt() < k_minPt_lep_tight or not passTightId or Muon_pfRelIso04_all()[i] > 0.15) // minPtTight
      continue;

    tightMuons.emplace_back(muon);
  }


  // Make sure the collections are sorted in pt
  switch (id_level) {
    case idLoose:
      std::sort(looseMuons.begin(), looseMuons.end(), PtOrdered);
      return looseMuons;
    case idTight:
      std::sort(tightMuons.begin(), tightMuons.end(), PtOrdered);
      return tightMuons;
  }
  return tightMuons;

}

vector<Photon> getPhotons(int id_level) {
  vector<Photon> photons;

  for (unsigned i = 0; i < Photon_pt().size(); ++i) {
    // Tight ID
    bool const passId = (Photon_cutBasedBitmap()[i] & 0b0100);
    
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
    const double corrPt = injet.p4.Pt();  // pt with nominal JEC

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

vector<Jet> getJets(bool reapplyJEC, bool isSim) {
  vector<Jet> jets;

  for (unsigned i = 0; i < Jet_pt().size(); ++i) {
    if (not (Jet_jetId()[i] & 0b0010))  // bit1 always false in 2017 since it does not exist
      continue;

    Jet jet;
    jet.p4.SetPtEtaPhiM(Jet_pt()[i], Jet_eta()[i], Jet_phi()[i], Jet_mass()[i]);
    jet.bTag = Jet_btagDeepFlavB()[i];
    if (isSim)
      jet.hadronFlavour = Jet_hadronFlavour()[i];
    else
      jet.hadronFlavour = 0;

    // Perform angular cleaning
    // if (IsDuplicate(jet.p4, 0.4)) continue;

    // if (reapplyJEC)
    //   jet.p4 *= getJetCorrectionFactorFromFile();

    // Kinematical cuts for jets to be stored in the collection
    if (jet.p4.Pt() < k_minPt_jet or std::abs(jet.p4.Eta()) > k_maxAbsEta_jet)
      continue;

    jets.emplace_back(jet);
  }

  // Make sure jets are sorted in pt
  std::sort(jets.begin(), jets.end(), PtOrdered);

  return jets;
}

