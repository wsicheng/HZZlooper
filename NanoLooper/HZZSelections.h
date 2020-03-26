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
double getPhotonTrigPrescale(double objpt);

vector<Jet> getJets(const vector<Muon>& mus, const vector<Electron>& els, const vector<Photon>& phs, bool reapplyJEC = false, bool isSim = true);
std::tuple<vector<Muon>, vector<Muon>> getMuons(bool applyRocCorr = true, float* shiftx = nullptr, float* shifty = nullptr);
std::tuple<vector<Electron>, vector<Electron>> getElectrons();
vector<Photon> getPhotons();
vector<int> getIsoTrackIndices(const vector<Muon>& mus, const vector<Electron>& els);

void giveMassToPhoton(TLorentzVector& boson, TH1* h_weight = nullptr);
void ApplyRochesterCorrectionToMuon(Muon *muon, int muidx);
float getJetCorrectionFactorFromFile(int jetidx, Jet injet, bool applyJER = true);
bool PassVBFcuts(const vector<Jet> &selJets, const TLorentzVector &boson);

bool ZZ2l2vPrunerCuts();

const vector<float> ZptScales_metlt125_eq0j  = {0.040,0.041,0.053,0.064,0.080,0.057,0.041,0.139,0.029,0.105,0.058,0.020,0.035,0.018,0.213,1.000,1.000,1.000,0.122,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,};
const vector<float> ZptScales_metlt125_geq1j = {0.034,0.041,0.050,0.063,0.074,0.091,0.096,0.112,0.122,0.138,0.141,0.157,0.138,0.147,0.149,0.153,0.182,0.170,0.174,0.171,0.200,0.218,0.213,0.173,0.213,0.185,0.213,0.186,0.271,0.218,0.244,0.268,0.192,0.210,0.280,0.186,0.321,0.280,0.222,0.161,0.309,0.253,0.227,0.192,0.267,0.262,0.367,0.121,0.151,0.079,0.311,0.308,0.120,0.220,0.241,};
const vector<float> ZptScales_metlt125_vbf   = {0.048,0.066,0.081,0.090,0.129,0.134,0.162,0.209,0.168,0.316,0.487,0.157,0.137,0.139,0.226,0.261,0.388,0.231,0.421,0.425,0.451,0.342,0.527,0.267,0.865,0.095,2.641,10.380,0.880,1.000,0.287,0.880,0.174,0.704,1.000,1.000,1.000,1.000,1.000,0.880,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,0.519,};

#endif
