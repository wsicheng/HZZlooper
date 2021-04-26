#ifndef HZZLOOPER_UTILITIES_H
#define HZZLOOPER_UTILITIES_H

#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"

#include <map>
#include <tuple>
#include <string>

namespace {

// Histogram manipulation
inline void moveOverFlowToLastBin1D(TH1* hist) {
  int nbin = hist->GetNbinsX();
  if (hist->GetBinLowEdge(nbin+1) < 1499.9) return;  // only when the last bin is the infinity bin
  if (hist->GetBinContent(nbin+1) > 0) {
    double err = 0;
    hist->SetBinContent(nbin, hist->IntegralAndError(nbin, -1, err));
    hist->SetBinError(nbin, err);
    hist->SetBinContent(nbin+1, 0);
    hist->SetBinError(nbin+1, 0);
  }
}

inline void moveXOverFlowToLastBin3D(TH3* hist) {
  int nbinx = hist->GetNbinsX();
  if (hist->GetXaxis()->GetBinLowEdge(nbinx+1) < 1499.9) return; // only when the last bin is infinity
  for (int ibiny = 0; ibiny < hist->GetNbinsY(); ++ibiny) {
    for (int ibinz = 0; ibinz < hist->GetNbinsZ(); ++ibinz) {
      if (hist->GetBinContent(nbinx+1, ibiny, ibinz) <= 0) continue;
      double err = 0;
      double yield = hist->IntegralAndError(nbinx, -1, ibiny, ibiny, ibinz, ibinz, err);
      hist->SetBinContent(nbinx, ibiny, ibinz, yield);
      hist->SetBinError(nbinx, ibiny, ibinz, err);
      hist->SetBinContent(nbinx+1, ibiny, ibinz, 0);
      hist->SetBinError(nbinx+1, ibiny, ibinz, 0);
    }
  }
}

// Weight each bin by the integral of the row 
inline void conditionalizeHistInY(TH2* h2d) {
  int nbinsX = h2d->GetNbinsX();
  int nbinsY = h2d->GetNbinsY();
  for (int ybin = 1; ybin <= nbinsY; ++ybin) {
    double xint = h2d->Integral(1, nbinsX, ybin, ybin);
    for (int xbin = 1; xbin <= nbinsX; ++xbin) {
      double bcorg = h2d->GetBinContent(xbin,ybin);
      h2d->SetBinContent(xbin, ybin, bcorg / xint);
    }
  }
}

// Weight each bin by the integral of the column
inline void conditionalizeHistInX(TH2* h2d) {
  int nbinsX = h2d->GetNbinsX();
  int nbinsY = h2d->GetNbinsY();
  for (int xbin = 1; xbin <= nbinsX; ++xbin) {
    double yint = h2d->Integral(xbin, xbin, 1, nbinsY);
    for (int ybin = 1; ybin <= nbinsY; ++ybin) {
      double bcorg = h2d->GetBinContent(xbin,ybin);
      h2d->SetBinContent(xbin, ybin, bcorg / yint);
    }
  }
}

// Quick helper functions
inline float deltaPhi(float phi1, float phi2) {
  float dphi = fabs(phi1 - phi2);
  if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
  return dphi;
}

inline float min2deltaPhi(float phi0, float phi1, float phi2) {
  return std::min(deltaPhi(phi0, phi1), deltaPhi(phi0, phi2));
}

inline float calculateMT(double Et1, double phi1, double Et2, double phi2) {
  return sqrt(2*Et1*Et2*(1.0 - cos(phi1-phi2)));
}

inline bool isCloseObject(const float eta1, const float phi1, const float eta2, const float phi2, const float conesize = 0.4, float* deltaR = nullptr) {
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

inline float phiFolding(float phi, float shift = 0.05) {
  phi += shift;
  while (phi > TMath::Pi()) phi -= TMath::TwoPi();
  while (phi < -TMath::Pi()) phi += TMath::TwoPi();
  phi -= TMath::Pi()/2;
  if (phi < -TMath::Pi()) phi += TMath::TwoPi();
  phi = fabs(phi) - TMath::Pi()/2;

  return phi;
}

// Use the same hist for a shared common denominator
inline void linkHist(vector<string> hnames, string hexist, std::map<std::string, TH1*> &allhistos) {
  for (string hnew : hnames) {
    if (allhistos.count(hnew)) return;
    auto iter = allhistos.find(hexist);
    if (iter == allhistos.end()) throw std::logic_error("linkHist(): Histogram "+hexist+" need to be plotted first");
    allhistos.insert( std::pair<std::string, TH1*>(hnew, iter->second) );
  }
}

// Templated function
template<class LorentzVectorType>
inline float calculateMt(const LorentzVectorType& l1p4, const LorentzVectorType& l2p4, double met, double metphi) {
  LorentzVectorType zp4 = l1p4 + l2p4;

  return calculateMT(zp4.Et(), zp4.Phi(), met, metphi);
}

template<class LorentzVectorType>
inline float getDileptonMT(const LorentzVectorType& boson, const LorentzVectorType& metvec, const double mV2 = 91.2) {
  double ptll = boson.Pt(), mll = boson.M(), met = metvec.Pt();
  double term1 = sqrt(ptll*ptll + mll*mll) + sqrt(met*met + mV2*mV2);
  double term2 = (boson + metvec).Pt();
  return sqrt( term1*term1 - term2*term2 );
}


template<class LorentzVectorType>
bool passVBFcuts(const vector<LorentzVectorType>& selJets, const LorentzVectorType &boson) {
  // The jets passed in are assumed to have already been pt sorted
  if (selJets.size() < 2) return false;

  float etamax = selJets[0].Eta();
  float etamin = selJets[1].Eta();
  if (etamax < etamin) std::swap(etamin, etamax);

  bool centralJetVeto = false;
  if (selJets.size() > 2) {
    for (unsigned int i = 2; i < selJets.size(); i++) {
      if (selJets[i].Eta() > etamin && selJets[i].Eta() < etamax)
        centralJetVeto = true;
    }
  }
  bool centralBoson = (boson.Eta() > etamin && boson.Eta() < etamax);
  bool passDeltaEta = ((etamax - etamin) > 4.);
  bool passMjj = ((selJets[0] + selJets[1]).M() > 500);

  if (!centralJetVeto && centralBoson && passDeltaEta && passMjj) return true;
  return false;
}

template<class HistType>
HistType* fetchHistCopy(std::string fname, std::string hname) {
  TFile file(fname.c_str());
  if (file.IsZombie()) {
    std::cout << "[fetchHistCopy] >> Cannot open file: " << fname << std::endl;
    return nullptr;
  }
  HistType* h = (HistType*) file.Get(hname.c_str());
  if (!h) std::cout << "[fetchHistCopy] >> Cannot find hist: " << hname << " in file: " << fname << std::endl;
  else h = (HistType*) h->Clone((hname+"_copy").c_str());

  file.Close();
  return h;
}

template<class LorentzVectorType>
bool isCloseObject(const LorentzVectorType& p1, const LorentzVectorType& p2, const float conesize = 0.4, float* deltaR = nullptr)
{
  const float PI = TMath::Pi();
  float deltaEta = fabs(p1.Eta() - p2.Eta());
  if (deltaEta > conesize) return false;
  float deltaPhi = fabs(p1.Phi() - p2.Phi());
  if (deltaPhi > PI) deltaPhi = 2*PI - deltaPhi;
  if (deltaPhi > conesize) return false;
  float deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
  if (deltaR2 > conesize*conesize) return false;
  if (deltaR) *deltaR = sqrt(deltaR2);

  return true;
}

template<typename... TArgs>
void plot1d(std::string name, double xval, double weight, std::map<std::string, TH1*> &allhistos, TArgs... args)
{
  auto iter = allhistos.find(name);
  if (iter == allhistos.end()) {
    TH1D* currentHisto= new TH1D(name.c_str(), args...);
    currentHisto->SetDirectory(0);
    currentHisto->Sumw2();
    currentHisto->Fill(xval, weight);
    allhistos.insert( std::pair<std::string, TH1*>(name, currentHisto) );
  } else {
    iter->second->Fill(xval, weight);
  }
}

template<typename... TArgs>
void plot2d(std::string name, double xval, double yval, double weight, std::map<std::string, TH1*> &allhistos, TArgs... args)
{
  auto iter = allhistos.find(name);
  if (iter == allhistos.end()) {
    TH2D* currentHisto= new TH2D(name.c_str(), args...);
    currentHisto->SetDirectory(0);
    currentHisto->Sumw2();
    currentHisto->Fill(xval, yval, weight);
    allhistos.insert( std::pair<std::string, TH1*>(name, currentHisto) );
  } else {
    ((TH2D*) iter->second)->Fill(xval, yval, weight);
  }
}

template<typename... TArgs>
void plot3d(std::string name, double xval, double yval, double zval, double weight, std::map<std::string, TH1*> &allhistos, TArgs... args)
{
  auto iter = allhistos.find(name);
  if (iter == allhistos.end()) {
    TH3D* currentHisto= new TH3D(name.c_str(), args...);
    currentHisto->SetDirectory(0);
    currentHisto->Sumw2();
    currentHisto->Fill(xval, yval, zval, weight);
    allhistos.insert( std::pair<std::string, TH1*>(name, currentHisto) );
  } else {
    ((TH3D*) iter->second)->Fill(xval, yval, zval, weight);
  }
}

// --------------------------------------------------------------
// Constants and hard coded scale factors
// --------------------------------------------------------------

std::tuple<vector<float>, vector<float>> getPtBins() {
  // vector<float> ptbin0 = {55, 60, 65, 70, 75, 80, 90, 105, 130, 160, 200, 250, 450};
  vector<float> ptbin0 = {55, 60, 65, 70, 75, 80, 90, 110, 130, 175, 450};

  vector<float> ptbin1;
  for (int i = 55;  i < 120; i += 5)  ptbin1.emplace_back(i);
  for (int i = 120; i < 150; i += 10) ptbin1.emplace_back(i);
  for (int i = 150; i < 180; i += 15) ptbin1.emplace_back(i);
  ptbin1.insert(ptbin1.end(), {180, 200, 230, 300, 600});

  return std::make_tuple(ptbin0, ptbin1);
}

std::tuple<vector<float>, vector<float>> getGenPtBins() {
  vector<float> ptbin1;
  for (int i =  0;  i <  60; i += 5)  ptbin1.emplace_back(i);
  for (int i = 60;  i < 100; i += 10) ptbin1.emplace_back(i);
  for (int i = 100; i < 140; i += 20) ptbin1.emplace_back(i);
  ptbin1.insert(ptbin1.end(), {140, 180, 220, 300, 600});

  vector<float> ptbin2{0, 5, 10, 20, 40, 60, 120, 200, 450};

  return std::make_tuple(ptbin1, ptbin2);
}

// Transfer factors from elec to photon
vector<float> tnpbin1 = { 55,    82.5,    100,    135,    175,    190,    220,   270,   450};
vector<float> tnpbin2 = { 55,    82.5,     100,     135,     220,      800}; // 800 = Inf

const map<int, vector<float>> sfval_elpt_tophCR_v4 = {
  {2016, {0.0191, 0.0158, 0.0127, 0.0093, 0.0111}},
  {2017, {0.0179, 0.0157, 0.0148, 0.0123, 0.0131}},
  {2018, {0.0190, 0.0169, 0.0152, 0.0140, 0.0125}},
};

const map<int, vector<float>> sferr_elpt_tophCR_v4 = {
  {2016, {0.0003, 0.0007, 0.0007, 0.0007, 0.0014}},
  {2017, {0.0003, 0.0006, 0.0006, 0.0007, 0.0014}},
  {2018, {0.0002, 0.0005, 0.0005, 0.0007, 0.0011}},
};

// Single Electron trigger eff, refer to https://indico.cern.ch/event/963616/contributions/4096636/ page 28,51,74 from Ulascan
// binned in pt as  vector<float> trigptbin1 = { 55,     100,     220,     inf};
const vector<float> trigeff_el_data_barrel_v1_16 = {0.896,   0.943,   0.995, };
const vector<float> trigeff_el_data_barrel_v1_17 = {0.841,   0.890,   0.995, };
const vector<float> trigeff_el_data_barrel_v1_18 = {0.867,   0.909,   0.995, };

// Single photon trigger eff, efer to: https://indico.cern.ch/event/879930/contributions/3974600/ page 10-12 from Sicheng
// binned in pt as vector<float> trigptbin0    = { 55,    82.5,    100,    135,    180,      220,      450};
const vector<float> trigeff_ph_data_barrel_v1_16 = {0.8935 , 0.9466 , 0.8981 , 0.8936 , 0.9196 , 0.9592 , };
const vector<float> trigeff_ph_data_barrel_v1_17 = {0.8938 , 0.9941 , 0.8729 , 0.9643 , 0.9584 , 0.9600 , };
const vector<float> trigeff_ph_data_barrel_v1_18 = {0.9376 , 0.9502 , 0.9565 , 0.9562 , 0.9571 , 0.9643 , };

const vector<float> tefferr_ph_data_barrel_v1_16 = {0.0212 , 0.0251 , 0.0169 , 0.0120 , 0.0048 , 0.0037 , };
const vector<float> tefferr_ph_data_barrel_v1_17 = {0.0265 , 0.0368 , 0.0222 , 0.0216 , 0.0136 , 0.0032 , };
const vector<float> tefferr_ph_data_barrel_v1_18 = {0.0294 , 0.0363 , 0.0275 , 0.0197 , 0.0197 , 0.0026 , };
 
const vector<float> trigSF_ph_barrel_v1_16 = { 0.9675, 1.0189, 0.9706, 0.9574, 0.9910, 0.9971, };
const vector<float> trigSF_ph_barrel_v1_17 = { 0.9427, 1.0460, 0.9222, 1.0188, 1.0169, 0.9971, };
const vector<float> trigSF_ph_barrel_v1_18 = { 0.9708, 0.9786, 0.9878, 0.9860, 0.9920, 0.9898, };

const vector<float> tSFerr_ph_barrel_v1_16 = { 0.0229, 0.0270, 0.0184, 0.0130, 0.0068, 0.0050, };
const vector<float> tSFerr_ph_barrel_v1_17 = { 0.0280, 0.0388, 0.0234, 0.0229, 0.0149, 0.0043, };
const vector<float> tSFerr_ph_barrel_v1_18 = { 0.0305, 0.0374, 0.0284, 0.0203, 0.0205, 0.0034, };

const vector<float> trigptbin0 = {55, 82.5, 100, 135, 180, 220, 450};
const vector<float> trigptbin1 = {55,       100,           220, 450};
const vector<float> trigptbin2 = {55, 82.5,  99, 135, 200, 220, 450};
const vector<float> trigptbin3 = {55, 82.5,  99, 132, 181.5, 230, 450};

std::pair<float,float> getPhotonTrigEffs(float Vpt, int year) {
  float scale = 1, sferr = 0;

  int trigcat0 = std::min(std::upper_bound(trigptbin0.begin(), trigptbin0.end(), Vpt) - trigptbin0.begin() - 1, 5L);
  if (year == 2016) {
    scale = trigeff_ph_data_barrel_v1_16.at(trigcat0);
    sferr = tefferr_ph_data_barrel_v1_16.at(trigcat0);
  }
  if (year == 2017) {
    scale = trigeff_ph_data_barrel_v1_17.at(trigcat0);
    sferr = tefferr_ph_data_barrel_v1_17.at(trigcat0);
  }
  if (year == 2018) {
    scale = trigeff_ph_data_barrel_v1_18.at(trigcat0);
    sferr = tefferr_ph_data_barrel_v1_18.at(trigcat0);
  }

  return std::make_pair(scale, sferr);
}

std::pair<float,float> getPhotonTrigSFs(float Vpt, int year) {
  float scale = 1, sferr = 0;

  int trigcat0 = std::min(std::upper_bound(trigptbin0.begin(), trigptbin0.end(), Vpt) - trigptbin0.begin() - 1, 5L);
  if (year == 2016) {
    scale = trigSF_ph_barrel_v1_16.at(trigcat0);
    sferr = tSFerr_ph_barrel_v1_16.at(trigcat0);
  }
  if (year == 2017) {
    scale = trigSF_ph_barrel_v1_17.at(trigcat0);
    sferr = tSFerr_ph_barrel_v1_17.at(trigcat0);
  }
  if (year == 2018) {
    scale = trigSF_ph_barrel_v1_18.at(trigcat0);
    sferr = tSFerr_ph_barrel_v1_18.at(trigcat0);
  }

  return std::make_pair(scale, sferr);
}

std::pair<float,float> getElecToGammaRate(float Vpt, int year, bool incl_trigeff = false) {
  // Apply the single electron to photon fake rate measured from tag and probe
  int icat = std::min(std::upper_bound(tnpbin2.begin(), tnpbin2.end(), Vpt) - tnpbin2.begin() - 1, 4L);
  if (icat < 0 || Vpt < 54) {
    cout << "[getElecToGammaRate] >> WARNING: Input Vpt= " << Vpt << ", icat= " << icat << endl;
    return std::make_pair(1, 0);
  }
  float scale = 1, sferr = 0;
  scale = sfval_elpt_tophCR_v4.at(year).at(icat);
  sferr = sferr_elpt_tophCR_v4.at(year).at(icat);

  if (!incl_trigeff) return std::make_pair(scale, sferr);

  // Apply the single-photon/single-electron trigger weight
  int trigcat0 = std::min(std::upper_bound(trigptbin0.begin(), trigptbin0.end(), Vpt) - trigptbin0.begin() - 1, 5L);
  int trigcat1 = std::min(std::upper_bound(trigptbin1.begin(), trigptbin1.end(), Vpt) - trigptbin1.begin() - 1, 2L);
  if (year == 2016) {
    scale *= trigeff_ph_data_barrel_v1_16.at(trigcat0) / trigeff_el_data_barrel_v1_16.at(trigcat1);
    if (incl_trigeff) sferr = sqrt(sferr*sferr + pow(tefferr_ph_data_barrel_v1_16.at(trigcat0), 2));
  }
  if (year == 2017) {
    scale *= trigeff_ph_data_barrel_v1_17.at(trigcat0) / trigeff_el_data_barrel_v1_17.at(trigcat1);
    if (incl_trigeff) sferr = sqrt(sferr*sferr + pow(tefferr_ph_data_barrel_v1_17.at(trigcat0), 2));
  }
  if (year == 2018) {
    scale *= trigeff_ph_data_barrel_v1_18.at(trigcat0) / trigeff_el_data_barrel_v1_18.at(trigcat1);
    if (incl_trigeff) sferr = sqrt(sferr*sferr + pow(tefferr_ph_data_barrel_v1_18.at(trigcat0), 2));
  }

  return std::make_pair(scale, sferr);
}

const vector<float> sfvals_phid_v1_16 = {0.939, 0.940, 0.949, 0.938, 0.965, 0.927, };
const vector<float> sfvals_phid_v1_17 = {0.934, 0.926, 0.923, 0.975, 0.897, 0.957, };
const vector<float> sfvals_phid_v1_18 = {0.940, 0.957, 0.904, 0.945, 0.970, 0.938, };
const vector<float> sfvals_phid_v1    = {0.938, 0.941, 0.929, 0.951, 0.947, 0.939, };
const vector<float> sferrs_phID_v1_16 = {0.007, 0.012, 0.012, 0.017, 0.026, 0.029, };
const vector<float> sferrs_phID_v1_17 = {0.008, 0.017, 0.016, 0.018, 0.041, 0.032, };
const vector<float> sferrs_phID_v1_18 = {0.008, 0.015, 0.018, 0.020, 0.036, 0.038, };

std::pair<float,float> getPhotonExtraIDSFs(float Vpt, int year) {
  float scale = 1, sferr = 0;

  int trigcat0 = std::min(std::upper_bound(trigptbin0.begin(), trigptbin0.end(), Vpt) - trigptbin0.begin() - 1, 5L);
  if (false) {
    if (year == 2016) {
      scale = sfvals_phid_v1_16.at(trigcat0);
      sferr = sferrs_phID_v1_16.at(trigcat0);
    }
    if (year == 2017) {
      scale = sfvals_phid_v1_17.at(trigcat0);
      sferr = sferrs_phID_v1_17.at(trigcat0);
    }
    if (year == 2018) {
      scale = sfvals_phid_v1_18.at(trigcat0);
      sferr = sferrs_phID_v1_18.at(trigcat0);
    }
  } else {
    scale = sfvals_phid_v1.at(trigcat0);
    sferr = (1.-scale) / 2.;
  }
  return std::make_pair(scale, sferr);
}

std::pair<float,float> getNjetSFfromLLG(int njet, int year, string systype) {
  const bool useDYsubtract = true;
  const bool useDiffInErr = true;
  const map<int, vector<float>> sfval_llgCR_v1 = {
    {2016, {1.060, 0.907, 0.785}},
    {2017, {1.036, 0.986, 0.831}},
    {2018, {0.954, 0.926, 0.872}},
  };
  const map<int, vector<float>> sferr_llgCR_v1 = {
    {2016, {0.0345, 0.0520, 0.0906}},
    {2017, {0.0348, 0.0542, 0.0781}},
    {2018, {0.0252, 0.0437, 0.0714}},
  };

  const map<int, vector<float>> sfval_llgCR_v2 = {
    {2016, {1.085, 0.874, 0.710}},
    {2017, {1.053, 0.994, 0.708}},
    {2018, {0.942, 0.896, 0.805}},
  };
  const map<int, vector<float>> sferr_llgCR_v2 = {
    {2016, {0.0449, 0.0792, 0.1592}},
    {2017, {0.0458, 0.0765, 0.1699}},
    {2018, {0.0348, 0.0733, 0.1423}},
  };

  float sfval = 1.0;
  float sferr = 0.0;

  if (year < 2016 || year > 2018) cout << "[getNjetSFfromLLG] Wrong year" << endl;

  if (!useDYsubtract) {
    sfval = sfval_llgCR_v1.at(year).at(std::min(njet, 2));
    sferr = sferr_llgCR_v1.at(year).at(std::min(njet, 2));
  } else {
    sfval = sfval_llgCR_v2.at(year).at(std::min(njet, 2));
    sferr = sferr_llgCR_v2.at(year).at(std::min(njet, 2));
  }
  if (useDiffInErr) {
    float diff = sfval_llgCR_v1.at(year).at(std::min(njet, 2)) - sfval_llgCR_v2.at(year).at(std::min(njet, 2));
    sferr = sqrt(sferr*sferr + diff*diff);
  }

  return std::make_pair(sfval, sferr);
}

std::map<std::string, TH2D*> loadScaleFactorHist(string type, TFile* f_sf) {
  std::map<std::string, TH2D*> hmap_sfs;
  if (type == "dilep") {
    hmap_sfs["ee:bb"]   = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_nominal_wcuts_ee_barrel_barrel_Nominal");
    hmap_sfs["ee:be"]   = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_nominal_wcuts_ee_barrel_endcap_Nominal");
    hmap_sfs["ee:eb"]   = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_nominal_wcuts_ee_endcap_barrel_Nominal");
    hmap_sfs["ee:ee"]   = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_nominal_wcuts_ee_endcap_endcap_Nominal");
    hmap_sfs["emu:bb"]  = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_nominal_wcuts_mue_barrel_barrel_Nominal");
    hmap_sfs["emu:be"]  = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_nominal_wcuts_mue_barrel_endcap_Nominal");
    hmap_sfs["emu:eb"]  = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_nominal_wcuts_mue_endcap_barrel_Nominal");
    hmap_sfs["emu:ee"]  = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_nominal_wcuts_mue_endcap_endcap_Nominal");
    hmap_sfs["mumu:bb"] = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_nominal_wcuts_mumu_barrel_barrel_Nominal");
    hmap_sfs["mumu:be"] = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_nominal_wcuts_mumu_barrel_endcap_Nominal");
    hmap_sfs["mumu:eb"] = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_nominal_wcuts_mumu_endcap_barrel_Nominal");
    hmap_sfs["mumu:ee"] = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_nominal_wcuts_mumu_endcap_endcap_Nominal");

    hmap_sfs["ee:bb:up"]   = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_up_wcuts_ee_barrel_barrel_Nominal");
    hmap_sfs["ee:be:up"]   = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_up_wcuts_ee_barrel_endcap_Nominal");
    hmap_sfs["ee:eb:up"]   = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_up_wcuts_ee_endcap_barrel_Nominal");
    hmap_sfs["ee:ee:up"]   = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_up_wcuts_ee_endcap_endcap_Nominal");
    hmap_sfs["emu:bb:up"]  = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_up_wcuts_mue_barrel_barrel_Nominal");
    hmap_sfs["emu:be:up"]  = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_up_wcuts_mue_barrel_endcap_Nominal");
    hmap_sfs["emu:eb:up"]  = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_up_wcuts_mue_endcap_barrel_Nominal");
    hmap_sfs["emu:ee:up"]  = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_up_wcuts_mue_endcap_endcap_Nominal");
    hmap_sfs["mumu:bb:up"] = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_up_wcuts_mumu_barrel_barrel_Nominal");
    hmap_sfs["mumu:be:up"] = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_up_wcuts_mumu_barrel_endcap_Nominal");
    hmap_sfs["mumu:eb:up"] = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_up_wcuts_mumu_endcap_barrel_Nominal");
    hmap_sfs["mumu:ee:up"] = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_up_wcuts_mumu_endcap_endcap_Nominal");
    hmap_sfs["ee:bb:dn"]   = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_dn_wcuts_ee_barrel_barrel_Nominal");
    hmap_sfs["ee:be:dn"]   = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_dn_wcuts_ee_barrel_endcap_Nominal");
    hmap_sfs["ee:eb:dn"]   = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_dn_wcuts_ee_endcap_barrel_Nominal");
    hmap_sfs["ee:ee:dn"]   = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_dn_wcuts_ee_endcap_endcap_Nominal");
    hmap_sfs["emu:bb:dn"]  = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_dn_wcuts_mue_barrel_barrel_Nominal");
    hmap_sfs["emu:be:dn"]  = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_dn_wcuts_mue_barrel_endcap_Nominal");
    hmap_sfs["emu:eb:dn"]  = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_dn_wcuts_mue_endcap_barrel_Nominal");
    hmap_sfs["emu:ee:dn"]  = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_dn_wcuts_mue_endcap_endcap_Nominal");
    hmap_sfs["mumu:bb:dn"] = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_dn_wcuts_mumu_barrel_barrel_Nominal");
    hmap_sfs["mumu:be:dn"] = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_dn_wcuts_mumu_barrel_endcap_Nominal");
    hmap_sfs["mumu:eb:dn"] = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_dn_wcuts_mumu_endcap_barrel_Nominal");
    hmap_sfs["mumu:ee:dn"] = (TH2D*) f_sf->Get("Dilepton_Combined/h_Combined_SF_dn_wcuts_mumu_endcap_endcap_Nominal");
  } else if (type == "sel_data_eff") {
    hmap_sfs["e:deff"]    = (TH2D*) f_sf->Get("SingleLepton_Combined/h_SingleElectron_eff_nominal_data_Nominal");
    hmap_sfs["e:deff:up"] = (TH2D*) f_sf->Get("SingleLepton_Combined/h_SingleElectron_eff_up_data_Nominal");
    hmap_sfs["e:deff:dn"] = (TH2D*) f_sf->Get("SingleLepton_Combined/h_SingleElectron_eff_dn_data_Nominal");
  }

  return hmap_sfs;
}


double getValFromHist2D(TH2D* h, double xval, double yval) {
  double xmin = h->GetXaxis()->GetBinLowEdge(1)+0.001;
  double xmax = h->GetXaxis()->GetBinLowEdge(h->GetNbinsX()+1)-0.001;
  double ymin = h->GetYaxis()->GetBinLowEdge(1)+0.001;
  double ymax = h->GetYaxis()->GetBinLowEdge(h->GetNbinsY()+1)-0.001;
  int binX = h->GetXaxis()->FindBin( std::max(xmin, std::min(xval, xmax)) ); // SCeta
  int binY = h->GetYaxis()->FindBin( std::max(ymin, std::min(yval, ymax)) );

  double scale = h->GetBinContent( binX, binY);
  return scale;
}

double getValFromHist1D(TH1* h, double val) {
  double xmin = h->GetXaxis()->GetBinLowEdge(1)+0.001;
  double xmax = h->GetXaxis()->GetBinLowEdge(h->GetNbinsX()+1)-0.001;
  int bin = h->GetXaxis()->FindBin( std::max(xmin, std::min(val, xmax)) ); // SCeta
  double scale = h->GetBinContent( bin);
  return scale;
}

std::tuple<float,float,float> getDilepTrigSF(int dilepid, float lep1pt, float lep1eta, float lep2pt, float lep2eta, std::map<std::string, TH2D*> hmap) {
  string htype;
  if (lep1pt < lep2pt) {
    std::swap(lep1pt, lep2pt);
    std::swap(lep1eta, lep2eta);
  }

  if (dilepid == -121) {
    htype = "ee:";
    if (fabs(lep1eta) < 1.479) htype += "b";
    else htype += "e";
    if (fabs(lep2eta) < 1.479) htype += "b";
    else htype += "e";
  }
  if (dilepid == -143) {
    htype = "emu:";
    if (fabs(lep1eta) < 1.2) htype += "b";
    else htype += "e";
     if (fabs(lep2eta) < 1.479) htype += "b";
    else htype += "e";
  }
  if (dilepid == -169) {
    htype = "mumu:";
    if (fabs(lep1eta) < 1.2) htype += "b";
    else htype += "e";
    if (fabs(lep2eta) < 1.2) htype += "b";
    else htype += "e";
  }
  
  float sfval = getValFromHist2D(hmap.at(htype), lep1pt, lep2pt);
  float sfval_up = getValFromHist2D(hmap.at(htype+":up"), lep1pt, lep2pt);
  float sfval_dn = getValFromHist2D(hmap.at(htype+":dn"), lep1pt, lep2pt);
  return std::make_tuple(sfval, sfval_up, sfval_dn);
}

std::tuple<float,float,float> getElecTrigEff(float pt, float eta, std::map<std::string, TH2D*> hmap) {
  string htype = "e:deff";
  float sfval = getValFromHist2D(hmap.at(htype), fabs(eta), pt);
  float sfval_up = getValFromHist2D(hmap.at(htype+":up"), fabs(eta), pt);
  float sfval_dn = getValFromHist2D(hmap.at(htype+":dn"), fabs(eta), pt);
  
  return std::make_tuple(sfval, sfval_up, sfval_dn);
}


bool loadExternalWeightHist(TString sample, string syst, int year, std::map<string,TH1*>& hmap) {
  if (year == 2018) {
    if (!(syst == "AsMZ" && sample == "WGToLNuG"))
      return false;
  }
  if (year == 2017) {
    if (syst == "AsMZ" && (sample != "WGToLNuG" && sample != "TGJets" && sample != "WZG" ))
      return false;
    else if (syst == "PythiaScale" && (sample == "TGJets" ))
      return false;
  }
  if (year == 2016) {
    if (syst == "AsMZ" && (sample != "WGToLNuG" && sample != "GJets" && sample != "ZJetsToNuNu" && sample != "QCD"))
      return false;
  }
  string bsuf = (syst == "AsMZ")? "_b1" : (syst == "PythiaScale")? "_b2" : "";

  TFile* fwgt = new TFile(Form("output/v4_08_phCR_all2jsel_rwgtd_wgtsbkup_2018/%s.root", sample.Data()));
  TH1F* ratup = (TH1F*) fwgt->Get(Form("effstudy_eff/ratio_genparts_pt%s_%sUp", bsuf.c_str(), syst.c_str()));
  TH1F* ratdn = (TH1F*) fwgt->Get(Form("effstudy_eff/ratio_genparts_pt%s_%sDn", bsuf.c_str(), syst.c_str()));
  if (!ratdn || !ratup) {
    cout << "[loadExternalWeightHist] >> Cannot find ratio plot! up: " << ratup << ", dn: " << ratdn << ", hname= " << Form("OffShell/ratio_genparts_pt_%sUp", syst.c_str()) << endl;
    return false;
  }
  hmap[syst+"Up"] = (TH1F*) ratup->Clone(Form("%s_%sUp", sample.Data(), syst.c_str())); hmap[syst+"Up"]->SetDirectory(0);
  hmap[syst+"Dn"] = (TH1F*) ratdn->Clone(Form("%s_%sDn", sample.Data(), syst.c_str())); hmap[syst+"Dn"]->SetDirectory(0);

  return true;
}

double getExtSampleWeight(TString fname, int year) {
  double sf = 1.0;
  if (year == 2018) {
    if (fname.Contains("WGToLNuG_01J_5f_PDFWeights_TuneCP5"))
      sf = 10259733782.9 / (10259733782.9 + 8994284552.6);
    else if (fname.Contains("WGToLNuG_01J_5f_TuneCP5"))
      sf = 8994284552.6 / (10259733782.9 + 8994284552.6);
    else if (fname.Contains("WGToLNuG"))
      sf = 0.;
  }
  else if (year == 2017) {
    if (fname.Contains("WGToLNuG_01J_5f_PDFWeights_TuneCP5"))
      sf = 9298794047.96 / (9298794047.96 + 9250613628.25);
    else if (fname.Contains("WGToLNuG_01J_5f_TuneCP5"))
      sf = 9250613628.25 / (9298794047.96 + 9250613628.25);
    else if (fname.Contains("WGToLNuG"))
      sf = 0.;
  }
  else if (year == 2016) {
    if (fname.Contains("WGToLNuG_01J_5f_Tune") && fname.Contains("_v3-v1"))
      sf = 8998720045.93 / (8998720045.93 + 17849780627.2);
    else if (fname.Contains("WGToLNuG_01J_5f_Tune") && fname.Contains("_v3_ext1-v1"))
      sf = 17849780627.2 / (8998720045.93 + 17849780627.2);
    else if (fname.Contains("WGToLNuG"))
      sf = 0.;
  }

  if (sf == 0.)
    // cout << "[getExtSampleWeight] >> WARNING: sample " << fname << " is not recognized!!"  << endl;
    throw std::invalid_argument("[getExtSampleWeight] >> Error: sample " + fname + " is not recognized!!");

  return sf;
}

double capVal(double val, double max=10.) {
  return (fabs(val) > max)? copysign(max, val) : val;
}


}

#endif  // HZZLOOPER_UTILITIES_H
