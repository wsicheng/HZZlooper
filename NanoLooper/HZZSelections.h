#ifndef NZZSelections_H
#define NZZSelections_H

#include "TRandom3.h"
#include "../NanoCORE/Nano.h"
#include "../NanoCORE/Config.h"
#include "../NanoCORE/RoccoR.h"
#include "PhysicsObjects.h"

const float k_minPt_photon = 55;
const float k_minPt_el_loose = 10;
const float k_minPt_mu_loose = 3;
const float k_minPt_lep_tight = 25;
const float k_minPt_jet = 30;
const float k_maxAbsEta_jet = 4.7;

const float mZ = 91.1876;

extern RoccoR* muoncorr;
extern TRandom3* randomGenerator;

enum ID_level { idVeto, idLoose, idMedium, idTight };
enum jet_cat { eq0j, geq1j, vbf };

float getDileptonMT(const TLorentzVector& boson, const TLorentzVector& metvec);
bool passTriggerSelections(int trigtype);

vector<Jet> getJets(const vector<Muon>& mus, const vector<Electron>& els, const vector<Photon>& phs, bool reapplyJEC = false, bool isSim = true);
std::tuple<vector<Muon>, vector<Muon>> getMuons(bool applyRocCorr = true, float* shiftx = nullptr, float* shifty = nullptr);
std::tuple<vector<Electron>, vector<Electron>> getElectrons();
vector<Photon> getPhotons();

void giveMassToPhoton(TLorentzVector& boson, TH1* h_weight = nullptr);
void ApplyRochesterCorrectionToMuon(Muon *muon, int muidx);
float getJetCorrectionFactorFromFile(int jetidx, Jet injet, bool applyJER = true);
bool PassVBFcuts(const vector<Jet> &selJets, const TLorentzVector &boson);

#endif
