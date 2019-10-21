#ifndef NZZSelections_H
#define NZZSelections_H

#include "../NanoCORE/Nano.h"
#include "../NanoCORE/Config.h"
#include "PhysicsObjects.h"

const float k_minPt_photon = 55;
const float k_minPt_lep_loose = 10;
const float k_minPt_lep_tight = 25;
const float k_minPt_jet = 30;
const float k_maxAbsEta_jet = 4.7;

const float MZ = 91.1876;

enum ID_level { idVeto, idLoose, idMedium, idTight };
enum jet_cat { eq0j, geq1j, vbf };


vector<Jet> getJets(bool reapplyJEC = false, bool isSim = true);
vector<Muon> getMuons(int id_level = idTight);
vector<Photon> getPhotons(int id_level = idTight);
vector<Electron> getElectrons(int id_level = idTight);

void giveMassToPhoton(TLorentzVector& boson, TH1* h_weight = nullptr);
std::optional<GenParticle> FindGenMatchMuon(Muon const &muon, double maxDR); 
void ApplyRochesterCorrectionToMuon(Muon *muon, int trackerLayers, bool isSim);
float getJetCorrectionFactorFromFile(int jetidx, Jet injet, bool applyJER = true);
bool PassVBFcuts(const vector<Jet> &selJets, const TLorentzVector &boson);

#endif
