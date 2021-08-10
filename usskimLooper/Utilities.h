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
    // cout << "Moving the overflow for hist: " << hist->GetTitle() << " to its last bin!" << endl;
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

inline void zeroOutNegativeYields(TH1* hist) {
  int nbin = hist->GetNbinsX();
  for (int ibin = 1; ibin <= nbin; ++ibin) {
    if (hist->GetBinContent(ibin) < 0) {
      if (string(hist->GetName()).find("h_metbins_") == string::npos)  // only print out for central hist
        cout << "Reverting negative yield " << hist->GetBinContent(ibin) << " in: " << hist->GetTitle() << " bin " << ibin << " to 0!" << endl;
      hist->SetBinContent(ibin, 0);
      // hist->SetBinError(ibin, 0); // should we set the error to 0 also?
    }
  }
}

inline void conditionalizeHistInY(TH2* h2d) {
  int nbinsX = h2d->GetNbinsX();
  int nbinsY = h2d->GetNbinsY();
  for (int ybin = 1; ybin <= nbinsY; ++ybin) {
    double xint = h2d->Integral(1,nbinsX, ybin, ybin);
    for (int xbin = 1; xbin <= nbinsX; ++xbin) {
      double bcorg = h2d->GetBinContent(xbin,ybin);
      h2d->SetBinContent(xbin, ybin, bcorg / xint);
    }
  }
}

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

inline float phiFolding(float phi, float shift = 0) {
  phi += shift;
  while (phi > TMath::Pi()) phi -= TMath::TwoPi();
  while (phi < -TMath::Pi()) phi += TMath::TwoPi();
  phi -= TMath::Pi()/2;
  if (phi < -TMath::Pi()) phi += TMath::TwoPi();
  phi = fabs(phi) - TMath::Pi()/2;

  return phi;
}

// Old functions that enforce float for ranges to be consistent with xval for floating point errors
void plot1D(string name, float xval, double weight, std::map<string, TH1*> &allhistos, string title, int numbinsx, float xmin, float xmax)
{
  if (title=="") title=name;
  std::map<string, TH1*>::iterator iter= allhistos.find(name);
  if (iter == allhistos.end()) { //no histo for this yet, so make a new one
    TH1D* currentHisto= new TH1D(name.c_str(), title.c_str(), numbinsx, xmin, xmax);
    currentHisto->SetDirectory(0);
    currentHisto->Sumw2();
    currentHisto->Fill(xval, weight);
    allhistos.insert(std::pair<string, TH1*>(name, currentHisto) );
  } else {
    iter->second->Fill(xval, weight);
  }
}

void linkHist(vector<string> hnames, string hexist, std::map<std::string, TH1*> &allhistos)
{
  // Could be useful when mutiple ratio hists sharing a common denominator
  for (string hnew : hnames) {
    if (allhistos.count(hnew)) return;
    auto iter = allhistos.find(hexist);
    if (iter == allhistos.end()) throw std::logic_error("linkHist(): Histogram "+hexist+" need to be plotted first");
    allhistos.insert( std::pair<std::string, TH1*>(hnew, iter->second) );
  }
}


// Templated function
template<class LorentzVectorType>
inline float calculateMt(const LorentzVectorType& l1p4, const LorentzVectorType& l2p4, double met, double met_phi){
  LorentzVectorType zp4 = l1p4 + l2p4;
  float phi1 = zp4.Phi();
  float phi2 = met_phi;
  float Et1  = zp4.Et();
  float Et2  = met;

  return sqrt(2*Et1*Et2*(1.0 - cos(phi1-phi2)));
}

template<class LorentzVectorType>
inline float getDileptonMT(const LorentzVectorType& boson, const LorentzVectorType& metvec) {
  double ptll = boson.Pt(), mll = boson.M(), met = metvec.Pt();
  double term1 = sqrt(ptll*ptll + mll*mll) + sqrt(met*met + 91.2*91.2);
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

// Constants 
const vector<float> ptRanges = {
  55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290,
  300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 7000
};
const vector<float> ptScales_geq1j = {
  0.2624, 0.4459, 1.0055, 1.5919, 2.0895, 2.5699, 3.0485, 3.5010, 3.8843, 4.2915, 4.6683, 5.1468, 5.4683, 5.7617, 5.9639,
  6.2891, 6.4758, 6.8194, 7.1370, 7.5179, 7.8242, 7.8640, 8.1332, 8.7345, 8.5413, 8.9472, 9.5241, 9.6099, 9.5511, 9.9231,
  9.8243, 10.152, 10.586, 10.708, 9.7402, 9.3267, 10.357, 10.338, 11.061, 12.649
}; // Zpeak scales

std::tuple<vector<float>, vector<float>> getPtBins() {


  // vector<float> ptbin0 = {55, 60, 65, 70, 75, 80, 90, 105, 130, 160, 200, 250, 450};
  vector<float> ptbin0 = {55, 60, 65, 70, 75, 80, 90, 110, 130, 175, 450};

  vector<float> ptbin1;
  for (int i = 55;  i < 120; i += 5)  ptbin1.emplace_back(i);
  for (int i = 120; i < 150; i += 10) ptbin1.emplace_back(i);
  for (int i = 150; i < 180; i += 15) ptbin1.emplace_back(i);
  ptbin1.insert(ptbin1.end(), {180, 200, 230, 300, 600});
  // for (int i = 180; i < 210; i += 30) ptbin1.emplace_back(i);
  // for (int i = 220; i < 270; i += 25) ptbin1.emplace_back(i);
  // ptbin1.emplace_back(270); ptbin1.emplace_back(300); ptbin1.emplace_back(350);
  // ptbin1.emplace_back(425); ptbin1.emplace_back(600);

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
vector<float> tnpbin1                    = { 55,    82.5,    100,    135,    175,    190,    220,   270,   450};
const vector<float> sfval_elpt_tophCR_v2_16 = {0.0320, 0.0260, 0.0184, 0.0136, 0.0154, 0.0095, 0.013,  0.0142,};
const vector<float> sfval_elpt_tophCR_v2_17 = {0.0352, 0.0289, 0.0242, 0.0196, 0.0214, 0.0183, 0.0155, 0.0177,};
const vector<float> sfval_elpt_tophCR_v2_18 = {0.0334, 0.0276, 0.024,  0.0202, 0.0209, 0.021,  0.0168, 0.0153,};
const vector<float> sferr_elpt_tophCR_v2_16 = {0.0004, 0.0009, 0.0008, 0.0011, 0.0025, 0.0016, 0.0021, 0.0025};
const vector<float> sferr_elpt_tophCR_v2_17 = {0.0004, 0.0008, 0.0008, 0.0011, 0.0027, 0.0022, 0.0020, 0.0024};
const vector<float> sferr_elpt_tophCR_v2_18 = {0.0003, 0.0007, 0.0007, 0.0010, 0.0022, 0.0019, 0.0018, 0.0019};

vector<float> tnpbin2                    = { 55,    82.5,     100,     135,     220,      800}; // 800 = Inf
const vector<float> sfval_elpt_tophCR_v3_16 = {0.0301,  0.0241,  0.0166,  0.0112,  0.0111};
const vector<float> sferr_elpt_tophCR_v3_16 = {0.0004,  0.0009,  0.0008,  0.0008,  0.0014};
const vector<float> sfval_elpt_tophCR_v3_17 = {0.0318,  0.0257,  0.0216,  0.0168,  0.0131};
const vector<float> sferr_elpt_tophCR_v3_17 = {0.0004,  0.0008,  0.0008,  0.0009,  0.0013};
const vector<float> sfval_elpt_tophCR_v3_18 = {0.0306,  0.0253,  0.0218,  0.0178,  0.0128};
const vector<float> sferr_elpt_tophCR_v3_18 = {0.0003,  0.0007,  0.0006,  0.0008,  0.0011};

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

// Single photon trigger eff, efer to: https://indico.cern.ch/event/879930/contributions/3974600/h page 10-12 from Sicheng
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


/// skimsamp scales, barrel tight photon
const vector<float> ZptScales_metlt125_eq0j = {0.140,0.156,0.179,0.213,0.296,0.307,0.366,0.635,0.529,0.585,1.319,1.323,0.501,1.020,3.260,7.808,0.213,0.435,1.000,1.000,1.000,1.611,1.000,1.000,1.817,1.000,1.341,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,};
const vector<float> ZptScales_metlt125_eq1j = {0.131,0.164,0.205,0.247,0.298,0.344,0.371,0.416,0.436,0.433,0.483,0.518,0.480,0.600,0.585,0.526,0.599,0.552,0.589,0.717,0.536,0.606,0.571,0.721,0.553,0.450,0.397,0.697,0.723,0.882,0.899,0.386,0.329,0.397,0.485,0.647,1.225,0.546,1.079,0.574,1.000,1.000,};
const vector<float> ZptScales_metlt125_eq2j = {0.100,0.125,0.169,0.209,0.248,0.294,0.336,0.380,0.430,0.449,0.463,0.415,0.515,0.576,0.579,0.615,0.676,0.715,0.608,0.618,0.784,0.523,0.444,0.482,0.779,0.452,0.462,0.627,0.657,0.484,0.642,0.270,0.641,0.771,1.034,0.567,0.565,0.371,0.161,0.525,1.000,1.000,};

const vector<float> ZptScales_met50to125_eq0j = {0.144,0.152,0.160,0.204,0.322,0.425,0.472,0.317,0.384,0.338,0.802,0.844,1.141,1.000,0.589,0.756,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,};
const vector<float> ZptScales_met50to125_eq1j = {0.121,0.172,0.207,0.259,0.275,0.352,0.375,0.426,0.458,0.587,0.548,0.517,0.471,0.779,0.546,0.545,0.611,0.580,0.573,0.874,0.364,1.076,0.651,0.687,0.618,0.595,0.280,0.797,1.797,0.490,0.565,0.498,0.279,0.206,1.533,0.040,0.898,1.245,0.475,0.683,1.000,1.000,};
const vector<float> ZptScales_met50to125_eq2j = {0.105,0.128,0.172,0.208,0.250,0.295,0.376,0.394,0.442,0.463,0.466,0.439,0.534,0.542,0.442,0.621,0.818,0.578,0.743,0.451,0.540,0.725,1.004,0.507,0.325,0.258,0.849,0.515,1.825,0.647,1.237,0.524,0.842,0.840,1.420,0.492,0.327,0.415,0.947,0.105,1.000,1.000,};

vector<float> sf_Vpt_allbkg_metlt80_eq0j  = {1.845e-02,8.454e-03,2.593e-02,2.279e-02,2.647e-02,3.508e-02,1.312e-02,3.422e-02,1.967e-02,1.698e-02,1.107e-02,7.435e-03,3.808e-02,1.408e-02,1.954e-04,1.000e+00,3.002e-05,7.189e-04,3.427e-02,9.150e-01,1.272e-03,1.113e-03,2.319e-04,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,};
vector<float> sf_Vpt_allbkg_metlt80_eq1j  = {1.680e-02,2.142e-02,2.802e-02,3.373e-02,3.791e-02,4.324e-02,4.942e-02,5.232e-02,5.951e-02,6.183e-02,6.399e-02,6.656e-02,7.454e-02,6.727e-02,7.755e-02,5.231e-02,8.793e-02,8.247e-02,7.848e-02,7.623e-02,8.040e-02,7.320e-02,6.604e-02,5.754e-02,7.805e-02,1.143e-01,5.071e-02,1.091e-01,1.175e-01,1.181e-01,1.015e-01,2.629e-01,1.180e-01,1.180e-01,7.835e-03,1.000e+00,8.876e-02,1.035e-01,2.168e-01,8.095e-02,1.000e+00,1.000e+00,1.545e-01,6.165e-02,1.507e-02,1.000e+00,1.014e-03,8.514e-01,1.419e-01,8.622e-02,1.031e-01,1.766e-01,1.000e+00,2.121e-01,2.281e-01,};
vector<float> sf_Vpt_allbkg_metlt80_eq2j  = {1.254e-02,1.567e-02,2.025e-02,2.843e-02,2.985e-02,3.841e-02,4.065e-02,4.851e-02,4.710e-02,5.875e-02,6.080e-02,7.860e-02,6.877e-02,5.523e-02,6.168e-02,6.086e-02,7.655e-02,7.825e-02,6.608e-02,9.100e-02,5.355e-02,1.003e-01,6.197e-02,9.371e-02,6.870e-02,1.130e-01,8.065e-02,1.029e-01,1.273e-01,8.087e-02,1.028e-01,5.363e-02,1.856e-01,4.622e-02,2.297e-01,1.673e-01,3.960e-01,2.558e-01,5.644e-01,6.793e-02,3.463e-01,6.684e-02,1.000e+00,1.395e-01,3.518e-02,2.737e-01,1.000e+00,7.010e-02,1.370e+00,2.972e-01,1.000e+00,3.544e-01,1.000e+00,4.505e-03,1.457e-01,};
vector<float> sf_Vpt_allbkg_metlt125_eq0j = {1.841e-02,8.598e-03,2.689e-02,2.501e-02,3.153e-02,4.382e-02,1.840e-02,4.869e-02,3.829e-02,3.059e-02,2.140e-02,1.559e-02,3.746e-02,1.597e-02,2.051e-04,1.000e+00,1.061e-03,1.708e-03,3.440e-02,3.177e-01,1.136e-03,1.108e-03,1.530e-04,1.000e+00,7.327e+00,1.000e+00,1.692e-03,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,};
vector<float> sf_Vpt_allbkg_metlt125_eq1j = {1.687e-02,2.156e-02,2.821e-02,3.404e-02,3.837e-02,4.357e-02,4.993e-02,5.294e-02,6.021e-02,6.204e-02,6.444e-02,6.680e-02,7.505e-02,6.815e-02,3.657e-02,5.201e-02,8.713e-02,8.233e-02,7.564e-02,7.690e-02,8.149e-02,7.363e-02,7.017e-02,6.031e-02,7.531e-02,1.127e-01,5.225e-02,1.065e-01,1.146e-01,1.096e-01,1.154e-01,2.807e-01,1.186e-01,1.310e-01,7.794e-03,2.017e-02,7.692e-02,1.036e-01,2.227e-01,8.137e-02,1.000e+00,1.000e+00,1.552e-01,4.774e-02,1.506e-02,1.000e+00,1.013e-03,8.512e-01,1.419e-01,8.758e-02,1.033e-01,1.766e-01,1.000e+00,2.120e-01,2.478e-01,};
vector<float> sf_Vpt_allbkg_metlt125_eq2j = {1.261e-02,1.585e-02,2.050e-02,2.858e-02,3.031e-02,3.857e-02,4.073e-02,4.909e-02,4.765e-02,5.970e-02,6.077e-02,7.874e-02,6.975e-02,5.741e-02,6.262e-02,5.904e-02,7.573e-02,8.063e-02,6.595e-02,8.743e-02,5.829e-02,9.744e-02,6.785e-02,9.083e-02,8.242e-02,1.078e-01,8.594e-02,1.091e-01,1.214e-01,8.054e-02,9.992e-02,6.321e-02,1.718e-01,6.847e-02,2.143e-01,1.846e-01,3.342e-01,2.747e-01,4.825e-01,6.487e-02,3.142e-01,3.398e-02,1.000e+00,1.399e-01,3.567e-02,2.477e-01,1.000e+00,6.665e-02,1.367e+00,2.970e-01,9.967e-03,2.926e-01,1.000e+00,2.902e-03,1.826e-01,};


// Extended photon eta at 2j
vector<float> sf_Vpt_closure16_metlt125_ee2j = {1.066e-02,1.412e-02,1.842e-02,2.435e-02,2.910e-02,3.755e-02,3.712e-02,4.351e-02,4.126e-02,5.259e-02,5.128e-02,6.832e-02,5.988e-02,5.568e-02,6.330e-02,5.212e-02,7.501e-02,6.832e-02,5.499e-02,6.995e-02,5.359e-02,8.321e-02,5.665e-02,8.114e-02,7.398e-02,9.243e-02,7.856e-02,8.470e-02,1.069e-01,8.410e-02,8.371e-02,7.157e-02,1.689e-01,4.014e-02,1.934e-01,1.924e-01,2.456e-01,1.982e-01,5.248e-01,6.587e-02,2.114e-01,4.370e-02,1.000e+00,1.091e-01,4.743e-02,2.845e-01,1.000e+00,6.882e-02,9.676e-01,2.526e-01,1.205e-02,1.147e-01,1.000e+00,1.000e+00,1.266e-01,};

vector<float> sf_Veta_data16_metlt80_ee2j = {5.267e-03,7.880e-03,8.682e-03,8.716e-03,1.179e-02,9.973e-03,1.170e-02,1.197e-02,1.453e-02,2.302e-02,4.678e-02,3.925e-02,3.778e-02,3.663e-02,3.146e-02,2.620e-02,2.778e-02,2.556e-02,2.424e-02,2.415e-02,2.498e-02,2.209e-02,2.090e-02,2.204e-02,2.089e-02,2.145e-02,2.257e-02,2.237e-02,2.124e-02,2.437e-02,2.474e-02,2.378e-02,2.493e-02,2.767e-02,2.670e-02,2.871e-02,3.785e-02,3.809e-02,4.122e-02,4.560e-02,2.054e-02,1.432e-02,1.254e-02,1.097e-02,1.010e-02,9.931e-03,9.272e-03,9.774e-03,8.030e-03,5.091e-03,};
vector<float> sf_Veta_data16_metlt125_ee2j = {5.242e-03,7.867e-03,8.671e-03,8.688e-03,1.171e-02,9.974e-03,1.166e-02,1.194e-02,1.441e-02,2.282e-02,4.710e-02,3.966e-02,3.827e-02,3.681e-02,3.177e-02,2.629e-02,2.791e-02,2.571e-02,2.430e-02,2.426e-02,2.507e-02,2.229e-02,2.104e-02,2.204e-02,2.112e-02,2.149e-02,2.268e-02,2.249e-02,2.120e-02,2.450e-02,2.477e-02,2.405e-02,2.496e-02,2.766e-02,2.686e-02,2.882e-02,3.838e-02,3.817e-02,4.109e-02,4.540e-02,2.023e-02,1.439e-02,1.245e-02,1.092e-02,1.001e-02,9.880e-03,9.235e-03,9.801e-03,7.977e-03,5.064e-03,};

vector<float> sf_Vpt_flateta_data16_metlt80_ee2j = {5.099e-01,6.867e-01,9.581e-01,1.205e+00,1.421e+00,1.648e+00,1.799e+00,2.007e+00,2.170e+00,2.327e+00,2.425e+00,2.517e+00,2.480e+00,2.567e+00,2.742e+00,2.984e+00,2.646e+00,2.731e+00,2.608e+00,2.770e+00,2.272e+00,2.548e+00,3.079e+00,3.082e+00,3.471e+00,2.773e+00,2.646e+00,3.045e+00,2.288e+00,2.230e+00,2.355e+00,3.300e+00,2.978e+00,1.832e+00,3.796e+00,4.943e+00,2.504e+00,1.209e+00,2.747e+00,4.875e+00,3.708e+00,3.233e+00,3.584e+00,1.864e+00,1.168e+00,4.849e+00,2.980e+00,5.705e+00,1.000e+00,2.529e+00,4.202e+00,2.507e+00,7.759e+00,1.000e+00,1.753e+00,};
vector<float> sf_Vpt_flateta_data16_metlt125_ee2j = {5.139e-01,6.927e-01,9.671e-01,1.221e+00,1.444e+00,1.657e+00,1.805e+00,2.027e+00,2.188e+00,2.355e+00,2.500e+00,2.542e+00,2.491e+00,2.661e+00,2.800e+00,3.099e+00,2.815e+00,2.538e+00,2.607e+00,2.828e+00,2.251e+00,3.084e+00,3.160e+00,3.174e+00,3.471e+00,2.743e+00,2.705e+00,3.182e+00,2.245e+00,2.187e+00,2.421e+00,3.213e+00,3.050e+00,1.820e+00,4.040e+00,4.816e+00,2.433e+00,1.188e+00,3.634e+00,4.860e+00,3.485e+00,3.118e+00,3.511e+00,1.795e+00,1.163e+00,4.713e+00,2.918e+00,5.747e+00,1.000e+00,2.265e+00,4.523e+00,2.330e+00,6.380e+00,1.000e+00,2.004e+00,};

// With conversion veto, barrel only
vector<float> sf_Vpt_v0_data16_metlt80_eq0j  = {2.493e-02,2.846e-02,3.157e-02,3.216e-02,3.010e-02,2.516e-02,2.037e-02,1.483e-02,1.547e-02,1.185e-02,1.131e-02,4.923e-03,9.235e-03,6.851e-03,8.216e-03,6.216e-03,3.782e-03,2.288e-03,2.809e-03,3.922e-03,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,3.006e-02,0.000e+00,3.448e-02,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.429e-00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,};
vector<float> sf_Vpt_v0_data16_metlt80_eq1j  = {3.052e-02,3.824e-02,4.909e-02,6.054e-02,6.882e-02,7.855e-02,8.582e-02,9.047e-02,9.460e-02,1.021e-01,1.048e-01,1.033e-01,1.091e-01,1.204e-01,1.165e-01,1.121e-01,1.110e-01,1.062e-01,1.158e-01,1.101e-01,1.112e-01,1.117e-01,1.163e-01,1.316e-01,1.281e-01,1.203e-01,1.211e-01,1.411e-01,1.116e-01,1.029e-01,1.185e-01,1.365e-01,1.322e-01,1.261e-01,7.212e-02,8.509e-02,1.906e-01,6.258e-02,2.144e-01,1.042e-01,1.555e-01,1.001e-01,3.108e-01,9.386e-02,1.667e-01,1.251e-01,2.941e-01,1.500e-01,9.091e-02,1.113e-01,1.000e+00,2.000e-01,8.333e-02,1.000e+00,1.134e-01,};
vector<float> sf_Vpt_v0_data16_metlt80_eq2j  = {2.836e-02,3.613e-02,4.722e-02,5.822e-02,6.798e-02,7.803e-02,8.576e-02,9.444e-02,9.864e-02,1.048e-01,1.078e-01,1.116e-01,1.133e-01,1.099e-01,1.284e-01,1.335e-01,1.225e-01,1.175e-01,1.118e-01,1.245e-01,9.064e-02,1.124e-01,1.340e-01,1.353e-01,1.527e-01,1.186e-01,1.165e-01,1.347e-01,1.042e-01,9.580e-02,1.012e-01,1.458e-01,1.114e-01,7.567e-02,1.648e-01,2.112e-01,1.001e-01,1.368e-01,1.001e-01,2.222e-01,1.516e-01,1.483e-01,1.350e-01,7.703e-02,5.259e-02,1.744e-01,1.055e-01,2.003e-01,1.000e+00,1.004e-01,1.442e-01,7.176e-02,2.500e-01,1.000e+00,7.271e-02,};
vector<float> sf_Vpt_v0_data16_metlt125_eq0j = {2.499e-02,2.872e-02,3.234e-02,3.384e-02,3.264e-02,2.968e-02,2.603e-02,1.952e-02,2.170e-02,1.524e-02,1.358e-02,6.620e-03,9.483e-03,6.423e-03,8.856e-03,1.537e-03,1.643e-02,2.092e-03,2.532e-03,6.735e-03,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,2.758e-02,0.000e+00,2.857e-02,0.000e+00,0.000e+00,0.000e+00,0.000e+00,8.333e-02,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,};
vector<float> sf_Vpt_v0_data16_metlt125_eq1j = {3.056e-02,3.831e-02,4.926e-02,6.067e-02,6.901e-02,7.895e-02,8.648e-02,9.090e-02,9.505e-02,1.029e-01,1.059e-01,1.039e-01,1.105e-01,1.206e-01,1.179e-01,1.129e-01,1.100e-01,1.049e-01,1.163e-01,1.109e-01,1.117e-01,1.117e-01,1.144e-01,1.311e-01,1.268e-01,1.200e-01,1.181e-01,1.412e-01,1.158e-01,1.042e-01,1.193e-01,1.339e-01,1.282e-01,1.214e-01,8.046e-02,9.000e-02,1.820e-01,5.981e-02,2.071e-01,9.620e-02,1.523e-01,9.681e-02,3.005e-01,8.828e-02,1.622e-01,1.073e-01,2.946e-01,1.500e-01,8.348e-02,1.001e-01,1.000e+00,2.000e-01,7.718e-02,1.000e+00,1.003e-01,};
vector<float> sf_Vpt_v0_data16_metlt125_eq2j = {2.844e-02,3.628e-02,4.745e-02,5.849e-02,6.834e-02,7.835e-02,8.606e-02,9.468e-02,9.968e-02,1.054e-01,1.089e-01,1.116e-01,1.131e-01,1.104e-01,1.301e-01,1.329e-01,1.241e-01,1.165e-01,1.113e-01,1.262e-01,8.975e-02,1.183e-01,1.355e-01,1.387e-01,1.544e-01,1.154e-01,1.181e-01,1.347e-01,1.015e-01,9.104e-02,1.017e-01,1.403e-01,1.121e-01,7.327e-02,1.709e-01,1.976e-01,9.738e-02,1.328e-01,1.275e-01,2.112e-01,1.431e-01,1.431e-01,1.317e-01,7.172e-02,5.264e-02,1.543e-01,1.055e-01,2.003e-01,1.000e+00,1.004e-01,1.446e-01,6.695e-02,2.222e-01,1.000e+00,8.771e-02,};

vector<float> sf_Vpt_v1_data16_metlt125_eq0j = {2.580e-02,3.082e-02,3.382e-02,3.346e-02,3.117e-02,2.759e-02,2.216e-02,2.065e-02,1.697e-02,1.569e-02,7.314e-03,8.030e-03,5.323e-03,};
vector<float> sf_Vpt_v1_data16_metlt125_eq1j = {3.302e-02,4.339e-02,5.473e-02,6.488e-02,7.416e-02,8.262e-02,8.959e-02,9.296e-02,9.858e-02,1.025e-01,1.051e-01,1.146e-01,1.144e-01,};
vector<float> sf_Vpt_v1_data16_metlt125_eq2j = {3.134e-02,4.069e-02,5.259e-02,6.477e-02,7.251e-02,8.110e-02,9.195e-02,9.669e-02,1.018e-01,1.076e-01,1.113e-01,1.117e-01,1.235e-01,};

vector<float> sf_Vpt_data16_metlt125_eq0j = {2.580e-02,3.082e-02,3.382e-02,3.346e-02,3.117e-02,2.759e-02,2.216e-02,2.065e-02,1.697e-02,1.569e-02,7.314e-03,8.030e-03,5.323e-03,};
vector<float> sf_Vpt_data16_metlt125_eq1j = {3.302e-02,4.339e-02,5.473e-02,6.488e-02,7.416e-02,8.262e-02,8.959e-02,9.296e-02,9.858e-02,1.025e-01,1.051e-01,1.146e-01,1.144e-01,};
vector<float> sf_Vpt_data16_metlt125_eq2j = {3.134e-02,4.069e-02,5.259e-02,6.477e-02,7.251e-02,8.110e-02,9.195e-02,9.669e-02,1.018e-01,1.076e-01,1.113e-01,1.117e-01,1.235e-01,};

vector<float> sf_Vpt_b0_v4_closure16_metlt125_eq0j = {3.340e-02,3.726e-02,4.015e-02,4.457e-02,4.153e-02,4.555e-02,5.128e-02,4.919e-02,2.827e-02,3.570e-02,5.624e-02,4.794e-02,5.053e-02,1.776e-02,4.113e-02,1.462e-02,1.000e+00,4.991e-02,1.995e-02,8.085e-03,7.224e-02,1.000e+00,1.266e-01,1.038e-01,3.071e-02,1.000e+00,};
vector<float> sf_Vpt_b1_v4_closure16_metlt125_eq1j = {3.539e-02,4.184e-02,4.955e-02,5.383e-02,5.704e-02,7.347e-02,7.110e-02,7.613e-02,8.284e-02,8.247e-02,8.994e-02,9.677e-02,7.956e-02,9.150e-02,1.002e-01,9.490e-02,1.180e-01,1.042e-01,1.284e-01,1.061e-01,1.276e-01,1.308e-01,1.314e-01,9.212e-02,1.354e-01,1.347e-01,1.639e-01,1.057e-01,1.794e-01,1.812e-01,2.018e-01,1.232e-01,8.036e-02,7.885e-02,1.610e-01,1.619e-01,1.723e-01,3.611e-01,1.468e-01,};
vector<float> sf_Vpt_b1_v4_closure16_metlt125_eq2j = {3.229e-02,3.593e-02,3.989e-02,5.023e-02,6.435e-02,7.186e-02,7.058e-02,7.464e-02,9.260e-02,1.224e-01,8.868e-02,1.239e-01,8.930e-02,9.803e-02,1.130e-01,1.206e-01,1.492e-01,1.695e-01,1.107e-01,1.301e-01,1.174e-01,2.040e-01,1.852e-01,1.539e-01,1.305e-01,1.574e-01,1.177e-01,1.784e-01,1.488e-01,2.285e-01,9.234e-02,1.555e-01,1.043e-01,3.673e-01,2.273e-01,2.094e-01,1.699e-01,1.834e-01,1.586e-01,};
vector<float> sf_Vpt_b0_v4_closure17_metlt125_eq0j = {4.329e-02,4.441e-02,4.431e-02,4.771e-02,4.342e-02,6.675e-02,1.396e-01,7.701e-02,6.584e-02,1.851e-01,3.419e-01,8.314e-02,6.130e-01,1.000e+00,1.289e+00,1.000e+00,8.849e-01,1.000e+00,1.000e+00,1.000e+00,1.064e+00,9.023e-01,2.673e-01,1.000e+00,1.000e+00,3.394e+00,};
vector<float> sf_Vpt_b1_v4_closure17_metlt125_eq1j = {4.536e-02,4.766e-02,5.395e-02,5.928e-02,1.017e-01,8.207e-02,1.120e-01,1.015e-01,1.058e-01,1.112e-01,1.298e-01,1.199e-01,1.300e-01,8.973e-02,1.250e-01,1.673e-01,1.274e-01,6.253e-02,1.273e-01,1.378e-01,1.967e-01,1.664e-01,1.406e-01,9.929e-02,1.385e-01,1.645e-01,1.815e-01,3.386e-01,2.018e-01,2.497e-01,2.276e-01,1.512e-01,2.700e-01,7.895e-02,1.766e-01,3.111e-01,1.583e-01,2.277e-01,3.414e-01,};
vector<float> sf_Vpt_b1_v4_closure17_metlt125_eq2j = {4.083e-02,5.111e-02,5.772e-02,7.370e-02,1.206e-01,9.112e-02,8.211e-02,9.564e-02,1.034e-01,1.469e-01,1.762e-01,1.092e-01,1.661e-01,1.272e-01,1.594e-01,2.129e-01,1.000e+00,1.557e-01,1.492e-01,1.169e-01,1.202e-01,1.570e-01,2.860e-01,1.922e-01,6.517e-01,2.692e-01,6.144e-01,2.229e-01,2.169e-01,1.778e-01,3.861e+00,3.348e-01,3.259e-01,2.874e-01,3.657e-01,2.669e-01,2.059e-01,3.865e-01,1.794e-01,};
vector<float> sf_Vpt_b0_v4_closure18_metlt125_eq0j = {4.764e-02,5.272e-02,4.753e-02,5.511e-02,6.525e-02,7.742e-02,7.686e-02,1.065e-01,1.000e+00,9.942e-02,6.199e-02,6.051e-02,5.216e-01,1.000e+00,3.139e-01,2.308e+00,2.080e+00,2.149e-01,1.000e+00,1.000e+00,1.989e-01,1.527e+01,2.384e+02,6.834e-02,3.965e-01,2.123e-01,};
vector<float> sf_Vpt_b1_v4_closure18_metlt125_eq1j = {4.796e-02,5.990e-02,6.325e-02,7.313e-02,7.826e-02,8.567e-02,9.590e-02,9.594e-02,1.093e-01,1.297e-01,1.257e-01,1.315e-01,1.245e-01,1.316e-01,1.448e-01,1.298e-01,1.159e-01,1.730e-01,5.949e-02,1.856e-01,2.101e-01,1.580e-01,1.746e-01,2.372e-01,1.676e-01,2.449e-01,1.719e-01,1.837e-01,2.568e-01,4.536e-01,1.515e-01,7.761e-02,2.437e-01,8.897e-02,4.024e-01,1.286e-01,1.336e-01,1.884e-01,2.596e-01,};
vector<float> sf_Vpt_b1_v4_closure18_metlt125_eq2j = {4.607e-02,5.761e-02,6.833e-02,5.843e-02,9.529e-02,9.416e-02,1.103e-01,1.394e-01,1.046e-01,1.182e-01,1.125e-01,1.446e-01,1.345e-01,1.457e-01,1.858e-01,2.943e-01,1.805e-01,1.409e-01,1.474e-01,1.684e-01,1.583e-01,1.348e-01,2.298e-01,2.321e-01,3.028e-01,2.385e-01,3.246e-01,1.314e-01,2.702e-01,1.942e-01,2.515e-01,2.300e-01,1.818e-01,3.234e-01,1.735e-01,2.466e-01,2.718e-01,2.518e-01,4.276e-01,};
vector<float> sf_Vpt_b1_v4_closure16_metlt125_ee2j = {1.977e-02,2.331e-02,2.534e-02,3.189e-02,3.827e-02,4.614e-02,4.361e-02,4.968e-02,5.593e-02,7.599e-02,5.849e-02,7.718e-02,6.025e-02,6.864e-02,7.016e-02,6.952e-02,8.735e-02,1.044e-01,6.737e-02,7.907e-02,8.110e-02,1.349e-01,1.242e-01,9.754e-02,7.729e-02,1.034e-01,8.215e-02,1.222e-01,1.087e-01,1.576e-01,6.601e-02,1.344e-01,7.067e-02,2.550e-01,1.701e-01,1.425e-01,1.251e-01,1.457e-01,1.282e-01,};
vector<float> sf_Vpt_b1_v4_closure17_metlt125_ee2j = {2.734e-02,3.236e-02,3.834e-02,4.785e-02,7.736e-02,5.636e-02,5.505e-02,6.407e-02,7.106e-02,9.464e-02,1.195e-01,7.687e-02,1.162e-01,9.044e-02,1.139e-01,1.376e-01,1.000e+00,1.054e-01,9.407e-02,8.615e-02,9.484e-02,1.189e-01,2.025e-01,1.298e-01,4.699e-01,1.895e-01,4.295e-01,1.650e-01,1.685e-01,1.161e-01,2.880e+00,2.513e-01,2.433e-01,1.967e-01,2.505e-01,2.039e-01,1.720e-01,3.197e-01,1.391e-01,};
vector<float> sf_Vpt_b1_v4_closure18_metlt125_ee2j = {2.898e-02,3.802e-02,4.483e-02,3.833e-02,6.267e-02,6.245e-02,7.337e-02,9.136e-02,6.573e-02,7.919e-02,7.890e-02,9.094e-02,9.427e-02,1.010e-01,1.308e-01,1.975e-01,1.234e-01,9.305e-02,9.074e-02,1.014e-01,1.175e-01,8.904e-02,1.507e-01,1.670e-01,2.071e-01,1.606e-01,2.080e-01,9.149e-02,1.944e-01,1.332e-01,1.835e-01,1.730e-01,1.182e-01,2.306e-01,1.424e-01,1.879e-01,2.309e-01,1.950e-01,3.506e-01,};

vector<float> sf_Veta_flatpt_v4_closure16_metlt125_ee2j = {5.408e-01,5.519e-01,8.021e-01,6.709e-01,8.509e-01,7.651e-01,8.303e-01,9.830e-01,9.377e-01,1.865e+00,1.444e+00,1.690e+00,1.411e+00,1.574e+00,1.148e+00,1.176e+00,7.911e-01,9.934e-01,7.993e-01,9.082e-01,9.154e-01,7.783e-01,9.317e-01,8.342e-01,6.577e-01,8.561e-01,6.104e-01,1.086e+00,9.066e-01,1.388e+00,1.048e+00,9.040e-01,9.010e-01,1.005e+00,7.981e-01,1.296e+00,1.647e+00,1.306e+00,1.469e+00,1.409e+00,1.855e+00,9.924e-01,1.167e+00,9.898e-01,8.095e-01,7.008e-01,5.537e-01,7.354e-01,7.837e-01,5.326e-01,};
vector<float> sf_Veta_flatpt_v4_closure17_metlt125_ee2j = {8.413e-01,9.245e-01,5.722e-01,6.445e-01,5.744e-01,6.844e-01,7.393e-01,4.783e-01,6.931e-01,1.633e+00,1.561e+00,1.172e+00,1.296e+00,8.927e-01,9.621e-01,9.238e-01,6.836e-01,1.115e+00,7.414e-01,6.042e-01,6.883e-01,5.048e-01,6.408e-01,4.975e-01,6.226e-01,6.297e-01,6.860e-01,5.848e-01,3.343e-02,8.137e-01,1.065e+00,1.255e+00,6.586e-01,7.203e-01,8.750e-01,1.091e+00,1.060e+00,8.872e-01,2.235e+00,1.810e+00,3.835e+00,8.908e-01,1.143e+00,8.113e-01,7.211e-01,6.975e-01,6.851e-01,1.074e+00,8.213e-01,6.215e-01,};
vector<float> sf_Veta_flatpt_v4_closure18_metlt125_ee2j = {1.073e+00,5.090e-01,9.995e-01,8.511e-01,1.311e+00,5.222e-01,9.138e-01,1.022e+00,4.695e-01,2.936e+00,1.644e+00,1.226e+00,1.420e+00,1.226e+00,8.220e-01,9.426e-01,8.685e-01,1.276e+00,6.188e-01,8.162e-01,8.236e-01,6.559e-01,6.411e-01,7.533e-01,1.027e+00,7.094e-01,6.416e-01,8.091e-01,6.134e-01,1.080e+00,8.654e-01,8.255e-01,1.008e+00,7.856e-01,1.046e+00,1.214e+00,1.221e+00,2.172e+00,1.468e+00,1.592e+00,1.493e+00,1.429e+00,6.321e-01,1.095e+00,6.680e-01,8.524e-01,8.482e-01,9.259e-01,4.424e-01,1.529e+00,};
 
vector<float> sf_Vpt_b0_v4_data16_metlt125_eq0j = {2.381e-02,2.765e-02,3.030e-02,3.389e-02,3.443e-02,4.038e-02,3.968e-02,4.650e-02,5.029e-02,4.865e-02,5.354e-02,4.835e-02,5.514e-02,5.263e-02,5.787e-02,7.545e-02,6.092e-02,7.436e-02,5.600e-02,6.601e-02,2.979e-02,2.419e-02,4.464e-02,1.875e-02,3.715e-02,2.679e-02,};
vector<float> sf_Vpt_b1_v4_data16_metlt125_eq1j = {3.142e-02,3.825e-02,4.163e-02,4.709e-02,5.080e-02,5.595e-02,6.152e-02,6.565e-02,7.027e-02,8.022e-02,7.968e-02,8.374e-02,8.872e-02,8.778e-02,9.323e-02,9.339e-02,9.890e-02,1.003e-01,1.059e-01,9.846e-02,1.026e-01,1.024e-01,1.006e-01,1.124e-01,1.123e-01,1.084e-01,1.045e-01,9.649e-02,1.097e-01,1.080e-01,1.046e-01,1.056e-01,7.993e-02,1.079e-01,1.503e-01,1.191e-01,1.034e-01,1.023e-01,1.098e-01,};
vector<float> sf_Vpt_b1_v4_data16_metlt125_eq2j = {3.274e-02,3.962e-02,4.259e-02,5.066e-02,5.615e-02,6.142e-02,6.945e-02,7.713e-02,7.966e-02,8.499e-02,8.977e-02,8.661e-02,9.628e-02,1.014e-01,1.026e-01,1.016e-01,1.030e-01,1.142e-01,1.124e-01,1.131e-01,1.118e-01,1.089e-01,1.176e-01,1.172e-01,1.201e-01,1.261e-01,1.256e-01,1.017e-01,1.095e-01,1.304e-01,1.076e-01,1.300e-01,1.016e-01,1.256e-01,1.114e-01,1.191e-01,1.213e-01,1.215e-01,1.311e-01,};
vector<float> sf_Vpt_b1_v4_data16_metlt125_ee2j = {2.056e-02,2.417e-02,2.623e-02,3.126e-02,3.518e-02,3.858e-02,4.366e-02,4.733e-02,5.040e-02,5.327e-02,5.603e-02,5.492e-02,6.064e-02,6.255e-02,6.622e-02,6.516e-02,6.911e-02,7.328e-02,7.367e-02,7.457e-02,7.509e-02,7.351e-02,7.808e-02,7.714e-02,8.091e-02,8.611e-02,8.641e-02,7.106e-02,7.723e-02,9.327e-02,7.511e-02,9.025e-02,7.374e-02,8.503e-02,8.007e-02,8.362e-02,8.673e-02,9.188e-02,1.040e-01,};

vector<float> sf_Vpt_b0_v4_data17_metlt125_eq0j = {2.463e-02,2.778e-02,3.132e-02,4.066e-02,3.656e-02,5.203e-02,5.037e-02,5.960e-02,6.477e-02,7.409e-02,1.077e-01,8.105e-02,1.009e-01,1.159e-01,1.091e-01,1.462e-01,1.646e-01,1.244e-01,1.765e-01,2.125e-01,1.897e-01,2.456e-01,1.489e-01,1.333e-01,1.882e-01,7.843e-02,};
vector<float> sf_Vpt_b1_v4_data17_metlt125_eq1j = {3.457e-02,4.172e-02,4.809e-02,5.344e-02,6.112e-02,6.298e-02,6.903e-02,7.483e-02,7.657e-02,8.199e-02,8.892e-02,9.252e-02,1.021e-01,1.041e-01,9.858e-02,9.720e-02,1.079e-01,1.093e-01,1.075e-01,1.209e-01,1.146e-01,1.219e-01,1.182e-01,1.257e-01,1.275e-01,1.230e-01,1.122e-01,1.125e-01,1.079e-01,1.212e-01,1.270e-01,1.192e-01,1.205e-01,1.041e-01,1.269e-01,1.089e-01,1.192e-01,1.323e-01,1.117e-01,};
vector<float> sf_Vpt_b1_v4_data17_metlt125_eq2j = {3.466e-02,4.075e-02,4.866e-02,5.831e-02,5.943e-02,6.441e-02,7.384e-02,8.051e-02,8.925e-02,9.705e-02,1.003e-01,1.120e-01,1.126e-01,1.111e-01,1.190e-01,1.217e-01,1.291e-01,1.139e-01,1.295e-01,1.300e-01,1.309e-01,1.303e-01,1.383e-01,1.363e-01,1.404e-01,1.341e-01,1.447e-01,1.525e-01,1.419e-01,1.496e-01,1.390e-01,1.486e-01,1.354e-01,1.257e-01,1.687e-01,1.430e-01,1.366e-01,1.638e-01,1.393e-01,};
vector<float> sf_Vpt_b1_v4_data17_metlt125_ee2j = {2.269e-02,2.697e-02,3.168e-02,3.702e-02,4.114e-02,4.359e-02,4.790e-02,5.224e-02,5.930e-02,6.360e-02,6.634e-02,7.524e-02,7.459e-02,7.548e-02,8.015e-02,8.055e-02,8.965e-02,7.828e-02,8.777e-02,8.824e-02,8.781e-02,9.021e-02,9.299e-02,9.470e-02,9.563e-02,9.417e-02,1.004e-01,1.085e-01,1.019e-01,1.046e-01,9.946e-02,1.075e-01,9.793e-02,8.788e-02,1.209e-01,1.058e-01,1.029e-01,1.209e-01,1.136e-01,};

vector<float> sf_Vpt_b0_v4_data18_metlt125_eq0j = {2.530e-02,2.902e-02,3.279e-02,3.765e-02,4.587e-02,4.759e-02,4.829e-02,5.727e-02,6.227e-02,7.069e-02,1.007e-01,1.129e-01,1.005e-01,9.646e-02,1.627e-01,1.443e-01,1.342e-01,1.654e-01,1.505e-01,1.594e-01,9.467e-02,1.698e-01,1.318e-01,2.326e-01,8.108e-02,1.286e-01,};
vector<float> sf_Vpt_b1_v4_data18_metlt125_eq1j = {3.645e-02,4.078e-02,4.956e-02,5.268e-02,5.892e-02,6.413e-02,6.585e-02,7.310e-02,8.156e-02,8.617e-02,8.801e-02,9.361e-02,9.321e-02,1.010e-01,9.867e-02,1.059e-01,1.115e-01,1.096e-01,1.097e-01,1.136e-01,1.178e-01,1.158e-01,1.168e-01,1.213e-01,1.217e-01,1.237e-01,1.261e-01,1.072e-01,1.044e-01,1.340e-01,1.056e-01,1.078e-01,1.036e-01,1.107e-01,1.194e-01,1.017e-01,1.083e-01,8.902e-02,1.074e-01,};
vector<float> sf_Vpt_b1_v4_data18_metlt125_eq2j = {3.685e-02,4.506e-02,5.315e-02,5.733e-02,7.166e-02,7.277e-02,8.374e-02,8.873e-02,9.042e-02,9.162e-02,1.027e-01,1.081e-01,1.217e-01,1.152e-01,1.321e-01,1.231e-01,1.244e-01,1.290e-01,1.253e-01,1.272e-01,1.226e-01,1.336e-01,1.210e-01,1.303e-01,1.335e-01,1.437e-01,1.423e-01,1.470e-01,1.489e-01,1.396e-01,1.394e-01,1.615e-01,1.476e-01,1.362e-01,1.272e-01,1.399e-01,1.303e-01,1.335e-01,1.167e-01,};
vector<float> sf_Vpt_b1_v4_data18_metlt125_ee2j = {2.391e-02,2.903e-02,3.245e-02,3.632e-02,4.479e-02,4.314e-02,5.426e-02,5.528e-02,5.793e-02,5.823e-02,6.555e-02,7.031e-02,7.435e-02,7.454e-02,8.054e-02,8.061e-02,8.197e-02,8.371e-02,7.930e-02,8.360e-02,7.783e-02,9.014e-02,8.081e-02,8.521e-02,8.880e-02,9.396e-02,9.773e-02,1.008e-01,1.035e-01,9.690e-02,9.604e-02,1.132e-01,1.033e-01,9.832e-02,9.072e-02,9.934e-02,9.466e-02,1.024e-01,9.337e-02,};
    
vector<float> sf_Veta_flatpt_v4_data16_metlt125_ee2j = {4.141e-01,5.097e-01,7.067e-01,6.276e-01,7.627e-01,6.509e-01,8.585e-01,9.241e-01,1.105e+00,1.642e+00,1.689e+00,1.405e+00,1.378e+00,1.414e+00,1.074e+00,1.051e+00,1.047e+00,1.088e+00,9.869e-01,9.578e-01,1.001e+00,9.362e-01,9.309e-01,9.922e-01,8.810e-01,8.264e-01,9.600e-01,8.776e-01,8.338e-01,9.992e-01,1.014e+00,9.653e-01,1.059e+00,1.070e+00,1.081e+00,1.210e+00,1.317e+00,1.463e+00,1.409e+00,1.856e+00,1.723e+00,1.112e+00,1.018e+00,9.116e-01,6.298e-01,7.084e-01,6.963e-01,6.710e-01,6.097e-01,3.856e-01,};
vector<float> sf_Veta_flatpt_v4_data17_metlt125_ee2j = {6.647e-01,8.334e-01,8.248e-01,7.288e-01,9.871e-01,8.588e-01,9.191e-01,9.622e-01,1.076e+00,1.425e+00,1.616e+00,1.382e+00,1.502e+00,1.261e+00,1.099e+00,9.520e-01,1.023e+00,9.414e-01,9.208e-01,8.711e-01,8.639e-01,8.256e-01,8.139e-01,8.545e-01,8.415e-01,8.011e-01,8.294e-01,8.635e-01,8.364e-01,8.199e-01,8.457e-01,8.523e-01,9.073e-01,9.558e-01,1.073e+00,1.236e+00,1.379e+00,1.246e+00,1.416e+00,1.802e+00,1.484e+00,1.009e+00,1.011e+00,8.219e-01,7.215e-01,7.805e-01,8.393e-01,1.030e+00,1.099e+00,5.783e-01,};
vector<float> sf_Veta_flatpt_v4_data18_metlt125_ee2j = {4.601e-01,7.073e-01,9.607e-01,8.344e-01,9.273e-01,7.597e-01,8.256e-01,9.626e-01,1.028e+00,1.530e+00,1.674e+00,1.527e+00,1.413e+00,1.240e+00,1.125e+00,1.079e+00,9.906e-01,1.057e+00,9.478e-01,8.581e-01,9.073e-01,9.064e-01,8.605e-01,9.021e-01,8.448e-01,8.510e-01,7.670e-01,9.261e-01,9.139e-01,1.004e+00,9.063e-01,8.764e-01,1.098e+00,1.053e+00,1.053e+00,1.126e+00,1.408e+00,1.364e+00,1.162e+00,1.662e+00,1.565e+00,9.189e-01,9.698e-01,8.892e-01,7.256e-01,6.800e-01,8.515e-01,8.803e-01,6.053e-01,4.040e-01,};
    
vector<float> sf_Vpt_b0_v5_data18_metlt125_eq0j = {2.530e-02,2.902e-02,3.279e-02,3.765e-02,4.587e-02,4.786e-02,6.143e-02,1.075e-01,1.428e-01,1.614e-01,8.148e-02,1.286e-01,};
vector<float> sf_Vpt_b1_v5_data18_metlt125_eq1j = {3.645e-02,4.078e-02,4.956e-02,5.268e-02,5.892e-02,6.413e-02,6.585e-02,7.310e-02,8.156e-02,8.617e-02,8.801e-02,9.361e-02,9.321e-02,9.991e-02,1.084e-01,1.097e-01,1.164e-01,1.153e-01,1.215e-01,1.248e-01,1.105e-01,1.125e-01,1.102e-01,1.042e-01,9.167e-02,1.098e-01,};
vector<float> sf_Vpt_b1_v5_data18_metlt125_ee2j = {2.391e-02,2.903e-02,3.245e-02,3.632e-02,4.479e-02,4.314e-02,5.426e-02,5.528e-02,5.793e-02,5.823e-02,6.555e-02,7.031e-02,7.435e-02,7.729e-02,8.126e-02,8.169e-02,8.445e-02,8.282e-02,8.681e-02,9.565e-02,1.004e-01,1.033e-01,9.823e-02,9.746e-02,9.877e-02,9.556e-02,};
vector<float> sf_Vpt_b0_v5_data17_metlt125_eq0j = {2.463e-02,2.778e-02,3.132e-02,4.066e-02,3.656e-02,5.137e-02,6.389e-02,9.983e-02,1.573e-01,2.037e-01,1.389e-01,7.843e-02,};
vector<float> sf_Vpt_b1_v5_data17_metlt125_eq1j = {3.457e-02,4.172e-02,4.809e-02,5.344e-02,6.112e-02,6.298e-02,6.903e-02,7.483e-02,7.657e-02,8.199e-02,8.892e-02,9.252e-02,1.021e-01,1.016e-01,1.020e-01,1.085e-01,1.180e-01,1.211e-01,1.265e-01,1.183e-01,1.109e-01,1.253e-01,1.167e-01,1.131e-01,1.248e-01,1.164e-01,};
vector<float> sf_Vpt_b1_v5_data17_metlt125_ee2j = {2.269e-02,2.697e-02,3.168e-02,3.702e-02,4.114e-02,4.359e-02,4.790e-02,5.224e-02,5.930e-02,6.360e-02,6.634e-02,7.524e-02,7.459e-02,7.765e-02,8.456e-02,8.270e-02,8.689e-02,9.470e-02,9.511e-02,9.688e-02,1.064e-01,1.014e-01,1.005e-01,1.046e-01,1.114e-01,1.260e-01,};
vector<float> sf_Vpt_b0_v5_data16_metlt125_eq0j = {2.381e-02,2.765e-02,3.030e-02,3.389e-02,3.443e-02,4.009e-02,4.808e-02,5.275e-02,6.480e-02,3.249e-02,2.823e-02,2.679e-02,};
vector<float> sf_Vpt_b1_v5_data16_metlt125_eq1j = {3.142e-02,3.825e-02,4.163e-02,4.709e-02,5.080e-02,5.595e-02,6.152e-02,6.565e-02,7.027e-02,8.022e-02,7.968e-02,8.374e-02,8.872e-02,9.020e-02,9.585e-02,1.029e-01,1.014e-01,1.007e-01,1.124e-01,1.067e-01,1.039e-01,1.042e-01,1.067e-01,1.130e-01,1.000e-01,1.171e-01,};
vector<float> sf_Vpt_b1_v5_data16_metlt125_ee2j = {2.056e-02,2.417e-02,2.623e-02,3.126e-02,3.518e-02,3.858e-02,4.366e-02,4.733e-02,5.040e-02,5.327e-02,5.603e-02,5.492e-02,6.064e-02,6.424e-02,6.693e-02,7.346e-02,7.376e-02,7.734e-02,7.880e-02,8.624e-02,7.631e-02,8.580e-02,7.904e-02,8.484e-02,9.563e-02,1.024e-01,};

vector<float> sf_Vaeta_flatpt_v5_data16_metlt125_ee2j = {8.529e-01,9.756e-01,9.035e-01,8.834e-01,1.000e+00,9.855e-01,9.761e-01,1.073e+00,1.059e+00,1.065e+00,1.140e+00,1.364e+00,1.421e+00,1.408e+00,1.769e+00,1.682e+00,1.108e+00,9.705e-01,8.845e-01,6.401e-01,7.347e-01,6.611e-01,6.891e-01,5.568e-01,3.992e-01,};
vector<float> sf_Vaeta_flatpt_v5_data17_metlt125_ee2j = {8.203e-01,8.419e-01,8.381e-01,8.310e-01,8.411e-01,8.579e-01,8.854e-01,9.244e-01,9.886e-01,1.008e+00,1.165e+00,1.319e+00,1.367e+00,1.399e+00,1.704e+00,1.453e+00,1.042e+00,9.854e-01,8.680e-01,7.845e-01,8.785e-01,7.777e-01,9.145e-01,9.414e-01,6.193e-01,};
vector<float> sf_Vaeta_flatpt_v5_data18_metlt125_ee2j = {8.487e-01,8.302e-01,8.933e-01,9.103e-01,9.534e-01,8.811e-01,9.106e-01,1.077e+00,1.020e+00,1.065e+00,1.126e+00,1.320e+00,1.388e+00,1.323e+00,1.668e+00,1.547e+00,9.709e-01,9.660e-01,8.572e-01,7.417e-01,7.875e-01,8.432e-01,9.180e-01,6.521e-01,4.297e-01,};

vector<float> sf_Vpt_b0_v6_data18_metlt125_eq0j = {1.696e-02,1.885e-02,2.087e-02,2.405e-02,2.855e-02,3.047e-02,3.884e-02,7.013e-02,8.641e-02,1.069e-01,5.926e-02,1.143e-01,};
vector<float> sf_Vpt_b1_v6_data18_metlt125_eq1j = {2.018e-02,2.181e-02,2.590e-02,2.733e-02,2.994e-02,3.346e-02,3.449e-02,3.928e-02,4.292e-02,4.655e-02,4.724e-02,5.183e-02,4.937e-02,5.548e-02,6.256e-02,6.294e-02,6.985e-02,6.883e-02,7.439e-02,7.991e-02,6.673e-02,7.230e-02,7.066e-02,7.689e-02,5.959e-02,6.393e-02,};
vector<float> sf_Vpt_b1_v6_data18_metlt125_ee2j = {1.090e-02,1.289e-02,1.427e-02,1.595e-02,1.926e-02,1.871e-02,2.338e-02,2.348e-02,2.539e-02,2.511e-02,2.998e-02,3.286e-02,3.540e-02,3.520e-02,3.845e-02,4.093e-02,4.145e-02,4.448e-02,4.334e-02,5.047e-02,5.696e-02,5.807e-02,5.899e-02,4.954e-02,5.525e-02,5.113e-02,};
vector<float> sf_Vpt_b0_v6_data17_metlt125_eq0j = {1.667e-02,1.815e-02,2.009e-02,2.634e-02,2.299e-02,3.331e-02,4.068e-02,6.241e-02,8.989e-02,1.235e-01,8.333e-02,7.843e-02,};
vector<float> sf_Vpt_b1_v6_data17_metlt125_eq1j = {1.926e-02,2.262e-02,2.579e-02,2.811e-02,3.232e-02,3.332e-02,3.685e-02,3.940e-02,4.118e-02,4.484e-02,4.832e-02,4.932e-02,5.771e-02,5.643e-02,5.756e-02,6.251e-02,7.030e-02,7.419e-02,7.757e-02,7.237e-02,7.130e-02,8.614e-02,7.782e-02,7.436e-02,7.956e-02,7.991e-02,};
vector<float> sf_Vpt_b1_v6_data17_metlt125_ee2j = {1.068e-02,1.222e-02,1.426e-02,1.642e-02,1.811e-02,1.821e-02,2.101e-02,2.380e-02,2.631e-02,2.861e-02,3.014e-02,3.360e-02,3.484e-02,3.681e-02,4.185e-02,4.193e-02,4.399e-02,5.013e-02,4.988e-02,5.671e-02,5.681e-02,5.454e-02,5.687e-02,5.742e-02,6.202e-02,5.837e-02,};
vector<float> sf_Vpt_b0_v6_data16_metlt125_eq0j = {1.604e-02,1.806e-02,1.980e-02,2.182e-02,2.276e-02,2.552e-02,3.228e-02,3.413e-02,3.922e-02,1.685e-02,2.016e-02,1.786e-02,};
vector<float> sf_Vpt_b1_v6_data16_metlt125_eq1j = {1.818e-02,2.110e-02,2.295e-02,2.581e-02,2.783e-02,3.067e-02,3.417e-02,3.668e-02,3.866e-02,4.506e-02,4.550e-02,4.767e-02,5.113e-02,5.339e-02,5.639e-02,6.258e-02,6.288e-02,6.161e-02,7.032e-02,6.980e-02,7.218e-02,6.901e-02,6.531e-02,7.005e-02,6.825e-02,6.286e-02,};
vector<float> sf_Vpt_b1_v6_data16_metlt125_ee2j = {9.913e-03,1.166e-02,1.215e-02,1.464e-02,1.629e-02,1.729e-02,2.097e-02,2.097e-02,2.362e-02,2.555e-02,2.677e-02,2.666e-02,3.014e-02,3.259e-02,3.572e-02,3.877e-02,4.018e-02,4.195e-02,4.466e-02,4.780e-02,4.462e-02,4.557e-02,4.784e-02,4.659e-02,5.785e-02,6.024e-02,};

vector<float> sf_Vaeta_flatpt_v6_data18_metlt125_ee2j = {7.996e-01,7.911e-01,8.306e-01,8.631e-01,9.188e-01,8.150e-01,8.678e-01,1.001e+00,9.459e-01,1.032e+00,1.061e+00,1.285e+00,1.377e+00,1.314e+00,1.739e+00,1.622e+00,1.027e+00,1.049e+00,9.106e-01,8.076e-01,8.220e-01,8.914e-01,1.002e+00,7.023e-01,4.846e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v6_data17_metlt125_ee2j = {7.659e-01,8.470e-01,8.066e-01,8.087e-01,8.094e-01,8.049e-01,8.190e-01,9.220e-01,9.454e-01,9.621e-01,1.152e+00,1.291e+00,1.275e+00,1.405e+00,1.736e+00,1.555e+00,1.041e+00,1.055e+00,9.258e-01,8.113e-01,8.972e-01,8.394e-01,9.667e-01,1.000e+00,6.765e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v6_data16_metlt125_ee2j = {8.127e-01,9.357e-01,8.470e-01,8.179e-01,9.157e-01,8.599e-01,9.443e-01,1.043e+00,1.003e+00,1.035e+00,1.148e+00,1.292e+00,1.446e+00,1.404e+00,1.842e+00,1.738e+00,1.178e+00,1.026e+00,9.676e-01,6.854e-01,7.862e-01,7.172e-01,7.322e-01,5.871e-01,4.590e-01,1.000e+00,};

vector<float> sf_Vpt_b0_v7_data16_metlt125_eq0j = {1.666e-02,1.867e-02,2.045e-02,2.264e-02,2.402e-02,2.693e-02,3.487e-02,3.722e-02,3.635e-02,2.659e-02,};
vector<float> sf_Vpt_b1_v7_data16_metlt125_eq1j = {1.853e-02,2.172e-02,2.359e-02,2.659e-02,2.901e-02,3.176e-02,3.605e-02,3.861e-02,4.156e-02,4.830e-02,4.951e-02,5.312e-02,5.537e-02,5.907e-02,6.217e-02,7.517e-02,7.526e-02,7.870e-02,8.879e-02,8.705e-02,1.012e-01,1.115e-01,};
vector<float> sf_Vpt_b1_v7_data16_metlt125_ee2j = {9.911e-03,1.157e-02,1.223e-02,1.447e-02,1.644e-02,1.738e-02,2.159e-02,2.139e-02,2.384e-02,2.660e-02,2.815e-02,2.938e-02,3.294e-02,3.466e-02,3.814e-02,4.296e-02,4.576e-02,5.008e-02,5.466e-02,5.894e-02,6.155e-02,8.113e-02,};
vector<float> sf_Vpt_b0_v7_data17_metlt125_eq0j = {1.763e-02,1.926e-02,2.128e-02,2.801e-02,2.480e-02,3.609e-02,4.716e-02,6.023e-02,1.051e-01,1.481e-01,};
vector<float> sf_Vpt_b1_v7_data17_metlt125_eq1j = {1.972e-02,2.337e-02,2.678e-02,2.938e-02,3.384e-02,3.502e-02,3.902e-02,4.199e-02,4.434e-02,4.787e-02,5.164e-02,5.286e-02,6.359e-02,6.239e-02,6.441e-02,6.818e-02,8.145e-02,8.530e-02,9.429e-02,8.923e-02,1.017e-01,1.075e-01,};
vector<float> sf_Vpt_b1_v7_data17_metlt125_ee2j = {1.093e-02,1.255e-02,1.446e-02,1.670e-02,1.850e-02,1.895e-02,2.171e-02,2.504e-02,2.649e-02,3.065e-02,3.186e-02,3.415e-02,3.644e-02,3.846e-02,4.113e-02,4.505e-02,4.972e-02,5.556e-02,5.971e-02,6.382e-02,6.960e-02,8.394e-02,};
vector<float> sf_Vpt_b0_v7_data18_metlt125_eq0j = {1.804e-02,2.018e-02,2.260e-02,2.551e-02,3.053e-02,3.271e-02,4.606e-02,8.163e-02,1.082e-01,1.507e-01,};
vector<float> sf_Vpt_b1_v7_data18_metlt125_eq1j = {2.097e-02,2.256e-02,2.699e-02,2.872e-02,3.151e-02,3.578e-02,3.656e-02,4.211e-02,4.602e-02,5.055e-02,5.222e-02,5.549e-02,5.524e-02,6.039e-02,6.939e-02,7.098e-02,7.893e-02,8.067e-02,9.044e-02,9.185e-02,9.476e-02,9.639e-02,};
vector<float> sf_Vpt_b1_v7_data18_metlt125_ee2j = {1.130e-02,1.273e-02,1.446e-02,1.658e-02,1.958e-02,1.988e-02,2.414e-02,2.468e-02,2.591e-02,2.678e-02,3.188e-02,3.431e-02,3.577e-02,3.739e-02,4.064e-02,4.318e-02,4.790e-02,5.038e-02,5.288e-02,6.046e-02,6.695e-02,7.112e-02,};

vector<float> sf_Vaeta_flatpt_v7_data16_metlt125_ee2j = {7.705e-01,8.681e-01,8.228e-01,8.016e-01,8.715e-01,8.701e-01,9.485e-01,9.947e-01,9.777e-01,9.883e-01,1.135e+00,1.275e+00,1.435e+00,1.442e+00,1.892e+00,1.835e+00,1.219e+00,1.064e+00,9.941e-01,7.329e-01,8.210e-01,7.600e-01,6.942e-01,6.380e-01,4.438e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v7_data17_metlt125_ee2j = {7.327e-01,8.060e-01,7.796e-01,7.665e-01,7.707e-01,7.847e-01,8.078e-01,8.821e-01,9.329e-01,9.437e-01,1.101e+00,1.303e+00,1.249e+00,1.404e+00,1.741e+00,1.669e+00,1.054e+00,1.075e+00,9.370e-01,8.878e-01,9.725e-01,9.234e-01,1.026e+00,1.048e+00,6.945e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v7_data18_metlt125_ee2j = {7.752e-01,7.571e-01,8.013e-01,8.389e-01,8.595e-01,7.891e-01,8.600e-01,9.790e-01,9.095e-01,9.811e-01,1.001e+00,1.283e+00,1.349e+00,1.342e+00,1.653e+00,1.635e+00,1.066e+00,1.099e+00,9.339e-01,8.776e-01,8.778e-01,9.469e-01,1.024e+00,7.568e-01,5.391e-01,1.000e+00,};

vector<float> sf_Vpt_b0_v8_data16_metlt125_eq0j = {1.680e-02,1.842e-02,2.083e-02,2.223e-02,2.478e-02,2.682e-02,3.350e-02,3.598e-02,3.417e-02,3.193e-02,};
vector<float> sf_Vpt_b1_v8_data16_metlt125_eq1j = {1.875e-02,2.152e-02,2.345e-02,2.720e-02,2.931e-02,3.297e-02,3.630e-02,3.969e-02,4.296e-02,4.807e-02,4.883e-02,5.256e-02,5.464e-02,5.931e-02,6.324e-02,7.005e-02,7.422e-02,7.885e-02,8.651e-02,8.831e-02,9.560e-02,1.172e-01,};
vector<float> sf_Vpt_b1_v8_data16_metlt125_ee2j = {9.877e-03,1.156e-02,1.280e-02,1.505e-02,1.683e-02,1.780e-02,2.106e-02,2.190e-02,2.416e-02,2.747e-02,2.910e-02,3.094e-02,3.265e-02,3.600e-02,3.873e-02,4.356e-02,4.706e-02,5.201e-02,5.384e-02,5.990e-02,6.496e-02,8.142e-02,};
vector<float> sf_Vpt_b0_v8_data17_metlt125_eq0j = {1.778e-02,1.938e-02,2.177e-02,2.615e-02,2.588e-02,3.561e-02,4.583e-02,6.919e-02,1.001e-01,1.357e-01,};
vector<float> sf_Vpt_b1_v8_data17_metlt125_eq1j = {1.996e-02,2.336e-02,2.621e-02,2.981e-02,3.249e-02,3.583e-02,4.004e-02,4.226e-02,4.512e-02,4.913e-02,5.127e-02,5.467e-02,5.885e-02,6.235e-02,6.655e-02,7.298e-02,8.074e-02,8.774e-02,8.993e-02,9.349e-02,1.044e-01,1.108e-01,};
vector<float> sf_Vpt_b1_v8_data17_metlt125_ee2j = {1.190e-02,1.345e-02,1.536e-02,1.785e-02,1.997e-02,2.146e-02,2.350e-02,2.560e-02,2.882e-02,3.023e-02,3.361e-02,3.439e-02,3.624e-02,4.072e-02,4.398e-02,4.780e-02,5.220e-02,5.412e-02,6.042e-02,6.745e-02,7.089e-02,8.449e-02,};
vector<float> sf_Vpt_b0_v8_data18_metlt125_eq0j = {1.836e-02,2.063e-02,2.275e-02,2.640e-02,2.922e-02,3.370e-02,4.593e-02,7.491e-02,1.003e-01,1.425e-01,};
vector<float> sf_Vpt_b1_v8_data18_metlt125_eq1j = {2.076e-02,2.323e-02,2.678e-02,2.925e-02,3.278e-02,3.567e-02,3.851e-02,4.304e-02,4.658e-02,5.130e-02,5.250e-02,5.574e-02,5.822e-02,6.279e-02,6.829e-02,7.193e-02,7.935e-02,8.167e-02,8.956e-02,9.427e-02,9.834e-02,9.764e-02,};
vector<float> sf_Vpt_b1_v8_data18_metlt125_ee2j = {1.163e-02,1.315e-02,1.461e-02,1.739e-02,1.961e-02,2.121e-02,2.313e-02,2.501e-02,2.700e-02,2.974e-02,3.223e-02,3.468e-02,3.682e-02,3.823e-02,4.249e-02,4.470e-02,5.012e-02,5.421e-02,5.699e-02,6.337e-02,6.873e-02,7.288e-02,};
vector<float> sferr_Vpt_b0_v8_data16_metlt125_eq0j = {1.811e-04,2.589e-04,3.814e-04,5.266e-04,7.441e-04,6.339e-04,8.029e-04,1.694e-03,2.466e-03,4.969e-03,};
vector<float> sferr_Vpt_b1_v8_data16_metlt125_eq1j = {1.614e-04,2.128e-04,2.643e-04,3.504e-04,4.261e-04,4.668e-04,4.219e-04,5.127e-04,6.110e-04,6.764e-04,7.584e-04,8.976e-04,1.026e-03,8.914e-04,1.074e-03,1.263e-03,1.319e-03,1.740e-03,2.048e-03,2.154e-03,2.475e-03,4.587e-03,};
vector<float> sferr_Vpt_b1_v8_data16_metlt125_ee2j = {1.502e-04,1.903e-04,2.274e-04,2.881e-04,3.480e-04,3.504e-04,3.478e-04,3.910e-04,4.481e-04,4.977e-04,5.607e-04,6.311e-04,7.052e-04,5.970e-04,6.983e-04,8.238e-04,8.244e-04,1.071e-03,1.165e-03,1.256e-03,1.358e-03,2.322e-03,};
vector<float> sferr_Vpt_b0_v8_data17_metlt125_eq0j = {2.443e-04,3.533e-04,5.270e-04,8.634e-04,1.059e-03,1.173e-03,1.567e-03,4.902e-03,9.451e-03,2.089e-02,};
vector<float> sferr_Vpt_b1_v8_data17_metlt125_eq1j = {1.999e-04,2.736e-04,3.549e-04,4.679e-04,5.820e-04,6.157e-04,5.491e-04,6.496e-04,7.552e-04,7.794e-04,9.039e-04,1.061e-03,1.261e-03,1.081e-03,1.166e-03,1.402e-03,1.525e-03,2.030e-03,1.995e-03,2.266e-03,2.529e-03,4.127e-03,};
vector<float> sferr_Vpt_b1_v8_data17_metlt125_ee2j = {1.781e-04,2.237e-04,2.810e-04,3.603e-04,4.407e-04,4.397e-04,3.870e-04,4.536e-04,5.372e-04,5.298e-04,6.251e-04,6.902e-04,7.773e-04,6.746e-04,7.306e-04,8.493e-04,8.652e-04,1.081e-03,1.125e-03,1.279e-03,1.305e-03,2.079e-03,};
vector<float> sferr_Vpt_b0_v8_data18_metlt125_eq0j = {2.390e-04,3.606e-04,5.269e-04,8.057e-04,1.188e-03,1.388e-03,1.398e-03,4.456e-03,7.415e-03,1.671e-02,};
vector<float> sferr_Vpt_b1_v8_data18_metlt125_eq1j = {2.065e-04,2.670e-04,3.612e-04,4.529e-04,5.878e-04,7.263e-04,8.937e-04,6.386e-04,7.511e-04,7.587e-04,8.618e-04,1.008e-03,1.168e-03,1.031e-03,1.037e-03,1.159e-03,1.255e-03,1.600e-03,1.698e-03,1.948e-03,1.918e-03,2.926e-03,};
vector<float> sferr_Vpt_b1_v8_data18_metlt125_ee2j = {1.786e-04,2.250e-04,2.733e-04,3.620e-04,4.511e-04,5.258e-04,6.292e-04,4.361e-04,4.924e-04,4.984e-04,5.739e-04,6.631e-04,7.472e-04,6.061e-04,6.309e-04,6.945e-04,7.205e-04,9.205e-04,9.295e-04,1.056e-03,1.045e-03,1.520e-03,};
vector<float> sf_Vaeta_flatpt_v8_data18_metlt125_ee2j = {7.646e-01,7.572e-01,7.656e-01,7.859e-01,8.452e-01,7.867e-01,8.103e-01,8.736e-01,9.080e-01,9.511e-01,1.023e+00,1.246e+00,1.293e+00,1.344e+00,1.675e+00,1.685e+00,1.125e+00,1.112e+00,1.031e+00,9.255e-01,9.611e-01,9.965e-01,1.030e+00,8.095e-01,5.233e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v8_data17_metlt125_ee2j = {7.529e-01,7.513e-01,7.289e-01,7.716e-01,8.176e-01,7.912e-01,8.005e-01,8.798e-01,9.088e-01,9.271e-01,1.041e+00,1.277e+00,1.304e+00,1.315e+00,1.771e+00,1.717e+00,1.117e+00,1.048e+00,9.538e-01,9.004e-01,9.929e-01,9.589e-01,1.067e+00,1.063e+00,7.168e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v8_data16_metlt125_ee2j = {7.886e-01,8.250e-01,8.198e-01,8.152e-01,8.711e-01,8.627e-01,8.883e-01,9.748e-01,9.770e-01,9.585e-01,1.086e+00,1.269e+00,1.408e+00,1.450e+00,1.805e+00,1.855e+00,1.248e+00,1.071e+00,9.655e-01,8.031e-01,8.567e-01,7.774e-01,7.541e-01,6.653e-01,4.659e-01,1.000e+00,};

vector<float> sf_Vpt_b0_v8_ee_data16_metlt125_eq0j = {6.331e-03,7.052e-03,8.092e-03,8.657e-03,9.508e-03,1.046e-02,1.335e-02,1.480e-02,1.532e-02,1.597e-02,};
vector<float> sf_Vpt_b1_v8_ee_data16_metlt125_eq1j = {7.139e-03,8.269e-03,8.949e-03,1.043e-02,1.113e-02,1.270e-02,1.407e-02,1.550e-02,1.731e-02,1.862e-02,1.971e-02,2.079e-02,2.193e-02,2.378e-02,2.462e-02,2.826e-02,2.912e-02,3.197e-02,3.529e-02,3.561e-02,3.824e-02,4.501e-02,};
vector<float> sf_Vpt_b1_v8_ee_data16_metlt125_ee2j = {3.792e-03,4.514e-03,4.835e-03,5.712e-03,6.447e-03,6.866e-03,8.446e-03,8.318e-03,9.601e-03,1.094e-02,1.165e-02,1.176e-02,1.315e-02,1.442e-02,1.551e-02,1.727e-02,1.965e-02,2.043e-02,2.130e-02,2.443e-02,2.630e-02,3.238e-02,};
vector<float> sf_Vpt_b0_v8_ee_data17_metlt125_eq0j = {6.968e-03,7.570e-03,8.479e-03,1.018e-02,9.934e-03,1.432e-02,1.755e-02,2.558e-02,3.272e-02,4.524e-02,};
vector<float> sf_Vpt_b1_v8_ee_data17_metlt125_eq1j = {7.829e-03,9.171e-03,1.040e-02,1.177e-02,1.294e-02,1.411e-02,1.633e-02,1.708e-02,1.844e-02,1.964e-02,2.087e-02,2.233e-02,2.415e-02,2.565e-02,2.754e-02,3.081e-02,3.187e-02,3.747e-02,3.634e-02,3.974e-02,4.278e-02,4.898e-02,};
vector<float> sf_Vpt_b1_v8_ee_data17_metlt125_ee2j = {4.680e-03,5.232e-03,6.081e-03,7.046e-03,7.871e-03,8.468e-03,9.335e-03,1.016e-02,1.176e-02,1.211e-02,1.357e-02,1.415e-02,1.504e-02,1.640e-02,1.803e-02,1.969e-02,2.191e-02,2.253e-02,2.418e-02,2.834e-02,3.102e-02,3.651e-02,};
vector<float> sf_Vpt_b0_v8_ee_data18_metlt125_eq0j = {7.287e-03,8.227e-03,8.928e-03,1.073e-02,1.172e-02,1.329e-02,1.901e-02,2.860e-02,4.393e-02,6.405e-02,};
vector<float> sf_Vpt_b1_v8_ee_data18_metlt125_eq1j = {8.259e-03,9.368e-03,1.069e-02,1.179e-02,1.324e-02,1.447e-02,1.561e-02,1.734e-02,1.897e-02,2.136e-02,2.209e-02,2.321e-02,2.415e-02,2.658e-02,2.852e-02,3.034e-02,3.361e-02,3.508e-02,3.807e-02,3.963e-02,4.229e-02,4.555e-02,};
vector<float> sf_Vpt_b1_v8_ee_data18_metlt125_ee2j = {4.627e-03,5.309e-03,5.854e-03,7.057e-03,8.053e-03,8.833e-03,9.275e-03,1.039e-02,1.093e-02,1.211e-02,1.372e-02,1.482e-02,1.526e-02,1.591e-02,1.738e-02,1.898e-02,2.121e-02,2.343e-02,2.390e-02,2.717e-02,2.954e-02,3.302e-02,};
vector<float> sferr_Vpt_b0_v8_ee_data16_metlt125_eq0j = {8.738e-05,1.248e-04,1.831e-04,2.512e-04,3.466e-04,3.227e-04,4.501e-04,9.999e-04,1.587e-03,3.459e-03,};
vector<float> sferr_Vpt_b1_v8_ee_data16_metlt125_eq1j = {7.776e-05,1.014e-04,1.241e-04,1.626e-04,1.944e-04,2.238e-04,2.282e-04,2.764e-04,3.361e-04,3.803e-04,4.352e-04,5.063e-04,5.821e-04,5.015e-04,6.084e-04,7.593e-04,7.779e-04,1.042e-03,1.224e-03,1.334e-03,1.524e-03,2.750e-03,};
vector<float> sferr_Vpt_b1_v8_ee_data16_metlt125_ee2j = {7.854e-05,9.875e-05,1.141e-04,1.424e-04,1.709e-04,1.806e-04,2.016e-04,2.190e-04,2.575e-04,2.949e-04,3.323e-04,3.622e-04,4.164e-04,3.494e-04,4.158e-04,5.004e-04,5.134e-04,6.438e-04,7.015e-04,7.888e-04,8.486e-04,1.431e-03,};
vector<float> sferr_Vpt_b0_v8_ee_data17_metlt125_eq0j = {1.127e-04,1.614e-04,2.369e-04,3.811e-04,4.638e-04,5.572e-04,7.650e-04,2.373e-03,4.488e-03,1.097e-02,};
vector<float> sferr_Vpt_b1_v8_ee_data17_metlt125_eq1j = {9.160e-05,1.234e-04,1.599e-04,2.075e-04,2.582e-04,2.790e-04,2.811e-04,3.282e-04,3.856e-04,4.148e-04,4.845e-04,5.657e-04,6.697e-04,5.718e-04,6.586e-04,8.225e-04,8.525e-04,1.181e-03,1.203e-03,1.415e-03,1.573e-03,2.666e-03,};
vector<float> sferr_Vpt_b1_v8_ee_data17_metlt125_ee2j = {8.640e-05,1.062e-04,1.330e-04,1.678e-04,2.028e-04,2.107e-04,2.064e-04,2.394e-04,2.881e-04,2.962e-04,3.481e-04,3.880e-04,4.371e-04,3.683e-04,4.251e-04,5.062e-04,5.181e-04,6.431e-04,6.861e-04,8.026e-04,8.474e-04,1.336e-03,};
vector<float> sferr_Vpt_b0_v8_ee_data18_metlt125_eq0j = {1.068e-04,1.603e-04,2.291e-04,3.572e-04,5.170e-04,5.897e-04,6.872e-04,2.069e-03,4.084e-03,1.012e-02,};
vector<float> sferr_Vpt_b1_v8_ee_data18_metlt125_eq1j = {9.154e-05,1.187e-04,1.575e-04,1.979e-04,2.554e-04,3.154e-04,3.860e-04,3.042e-04,3.613e-04,3.912e-04,4.466e-04,5.153e-04,5.928e-04,5.264e-04,5.676e-04,6.680e-04,7.191e-04,9.221e-04,1.021e-03,1.168e-03,1.225e-03,1.951e-03,};
vector<float> sferr_Vpt_b1_v8_ee_data18_metlt125_ee2j = {8.355e-05,1.052e-04,1.256e-04,1.655e-04,2.063e-04,2.418e-04,2.784e-04,2.251e-04,2.502e-04,2.684e-04,3.154e-04,3.627e-04,3.980e-04,3.225e-04,3.578e-04,4.165e-04,4.280e-04,5.503e-04,5.687e-04,6.544e-04,6.724e-04,1.004e-03,};
vector<float> sf_Vaeta_flatpt_v8_ee_data16_metlt125_ee2j = {3.250e-01,3.282e-01,3.350e-01,3.362e-01,3.460e-01,3.474e-01,3.592e-01,4.025e-01,3.831e-01,3.641e-01,4.416e-01,4.798e-01,5.323e-01,5.428e-01,6.875e-01,6.891e-01,4.709e-01,4.060e-01,3.820e-01,3.057e-01,3.191e-01,3.035e-01,3.044e-01,2.575e-01,1.819e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v8_ee_data17_metlt125_ee2j = {3.117e-01,3.116e-01,3.029e-01,3.351e-01,3.383e-01,3.273e-01,3.297e-01,3.649e-01,3.548e-01,3.870e-01,4.181e-01,5.119e-01,5.214e-01,5.208e-01,6.954e-01,6.762e-01,4.245e-01,4.077e-01,3.760e-01,3.496e-01,3.876e-01,3.747e-01,4.227e-01,4.340e-01,2.909e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v8_ee_data18_metlt125_ee2j = {3.099e-01,3.082e-01,3.029e-01,3.263e-01,3.512e-01,3.227e-01,3.330e-01,3.512e-01,3.630e-01,3.752e-01,4.071e-01,4.860e-01,5.059e-01,5.256e-01,6.574e-01,6.626e-01,4.581e-01,4.448e-01,4.208e-01,3.885e-01,4.061e-01,4.208e-01,4.721e-01,3.641e-01,2.502e-01,1.000e+00,};

vector<float> sf_Vpt_b0_v8_mumu_data16_metlt125_eq0j = {1.047e-02,1.137e-02,1.274e-02,1.358e-02,1.527e-02,1.636e-02,2.015e-02,2.117e-02,1.885e-02,1.597e-02,};
vector<float> sf_Vpt_b1_v8_mumu_data16_metlt125_eq1j = {1.161e-02,1.325e-02,1.450e-02,1.677e-02,1.818e-02,2.027e-02,2.222e-02,2.419e-02,2.565e-02,2.944e-02,2.911e-02,3.177e-02,3.271e-02,3.552e-02,3.863e-02,4.178e-02,4.510e-02,4.688e-02,5.122e-02,5.270e-02,5.736e-02,7.217e-02,};
vector<float> sf_Vpt_b1_v8_mumu_data16_metlt125_ee2j = {6.085e-03,7.041e-03,7.964e-03,9.334e-03,1.038e-02,1.093e-02,1.261e-02,1.358e-02,1.456e-02,1.654e-02,1.745e-02,1.918e-02,1.949e-02,2.158e-02,2.323e-02,2.629e-02,2.741e-02,3.159e-02,3.254e-02,3.547e-02,3.865e-02,4.903e-02,};
vector<float> sf_Vpt_b0_v8_mumu_data17_metlt125_eq0j = {1.081e-02,1.181e-02,1.329e-02,1.597e-02,1.594e-02,2.129e-02,2.828e-02,4.360e-02,6.737e-02,9.048e-02,};
vector<float> sf_Vpt_b1_v8_mumu_data17_metlt125_eq1j = {1.213e-02,1.419e-02,1.581e-02,1.805e-02,1.956e-02,2.172e-02,2.371e-02,2.518e-02,2.668e-02,2.949e-02,3.040e-02,3.235e-02,3.470e-02,3.670e-02,3.901e-02,4.217e-02,4.888e-02,5.027e-02,5.359e-02,5.375e-02,6.158e-02,6.185e-02,};
vector<float> sf_Vpt_b1_v8_mumu_data17_metlt125_ee2j = {7.219e-03,8.216e-03,9.278e-03,1.080e-02,1.210e-02,1.299e-02,1.417e-02,1.544e-02,1.706e-02,1.812e-02,2.004e-02,2.024e-02,2.120e-02,2.432e-02,2.595e-02,2.812e-02,3.029e-02,3.159e-02,3.624e-02,3.911e-02,3.987e-02,4.797e-02,};
vector<float> sf_Vpt_b0_v8_mumu_data18_metlt125_eq0j = {1.107e-02,1.240e-02,1.383e-02,1.567e-02,1.749e-02,2.041e-02,2.692e-02,4.631e-02,5.634e-02,7.843e-02,};
vector<float> sf_Vpt_b1_v8_mumu_data18_metlt125_eq1j = {1.250e-02,1.386e-02,1.609e-02,1.746e-02,1.955e-02,2.120e-02,2.291e-02,2.570e-02,2.761e-02,2.994e-02,3.041e-02,3.254e-02,3.407e-02,3.621e-02,3.977e-02,4.160e-02,4.574e-02,4.659e-02,5.149e-02,5.464e-02,5.605e-02,5.210e-02,};
vector<float> sf_Vpt_b1_v8_mumu_data18_metlt125_ee2j = {7.006e-03,7.844e-03,8.759e-03,1.033e-02,1.156e-02,1.238e-02,1.385e-02,1.462e-02,1.607e-02,1.763e-02,1.851e-02,1.986e-02,2.156e-02,2.231e-02,2.511e-02,2.572e-02,2.891e-02,3.078e-02,3.309e-02,3.620e-02,3.919e-02,3.985e-02,};
vector<float> sferr_Vpt_b0_v8_mumu_data16_metlt125_eq0j = {1.254e-04,1.769e-04,2.569e-04,3.528e-04,4.989e-04,4.388e-04,5.774e-04,1.228e-03,1.774e-03,3.459e-03,};
vector<float> sferr_Vpt_b1_v8_mumu_data16_metlt125_eq1j = {1.107e-04,1.441e-04,1.787e-04,2.345e-04,2.854e-04,3.164e-04,3.034e-04,3.660e-04,4.307e-04,4.976e-04,5.479e-04,6.516e-04,7.394e-04,6.394e-04,7.910e-04,9.405e-04,9.898e-04,1.288e-03,1.507e-03,1.636e-03,1.883e-03,3.527e-03,};
vector<float> sferr_Vpt_b1_v8_mumu_data16_metlt125_ee2j = {1.068e-04,1.329e-04,1.602e-04,2.006e-04,2.398e-04,2.463e-04,2.541e-04,2.910e-04,3.279e-04,3.708e-04,4.160e-04,4.760e-04,5.195e-04,4.392e-04,5.197e-04,6.253e-04,6.128e-04,8.129e-04,8.805e-04,9.555e-04,1.035e-03,1.774e-03,};
vector<float> sferr_Vpt_b0_v8_mumu_data17_metlt125_eq0j = {1.600e-04,2.310e-04,3.428e-04,5.571e-04,6.896e-04,7.611e-04,1.077e-03,3.448e-03,7.146e-03,1.630e-02,};
vector<float> sferr_Vpt_b1_v8_mumu_data17_metlt125_eq1j = {1.304e-04,1.770e-04,2.270e-04,2.987e-04,3.683e-04,3.980e-04,3.669e-04,4.344e-04,5.028e-04,5.421e-04,6.219e-04,7.254e-04,8.565e-04,7.309e-04,8.172e-04,9.912e-04,1.103e-03,1.413e-03,1.487e-03,1.665e-03,1.905e-03,3.015e-03,};
vector<float> sferr_Vpt_b1_v8_mumu_data17_metlt125_ee2j = {1.193e-04,1.496e-04,1.847e-04,2.355e-04,2.868e-04,2.918e-04,2.709e-04,3.158e-04,3.689e-04,3.790e-04,4.431e-04,4.846e-04,5.419e-04,4.735e-04,5.261e-04,6.193e-04,6.234e-04,7.807e-04,8.506e-04,9.514e-04,9.647e-04,1.540e-03,};
vector<float> sferr_Vpt_b0_v8_mumu_data18_metlt125_eq0j = {1.523e-04,2.281e-04,3.350e-04,4.990e-04,7.387e-04,8.688e-04,8.993e-04,2.996e-03,4.846e-03,1.143e-02,};
vector<float> sferr_Vpt_b1_v8_mumu_data18_metlt125_eq1j = {1.307e-04,1.669e-04,2.261e-04,2.809e-04,3.629e-04,4.460e-04,5.480e-04,4.143e-04,4.849e-04,4.994e-04,5.631e-04,6.600e-04,7.642e-04,6.634e-04,7.066e-04,8.103e-04,8.707e-04,1.101e-03,1.214e-03,1.403e-03,1.420e-03,2.093e-03,};
vector<float> sferr_Vpt_b1_v8_mumu_data18_metlt125_ee2j = {1.162e-04,1.443e-04,1.750e-04,2.283e-04,2.809e-04,3.235e-04,3.948e-04,2.879e-04,3.297e-04,3.438e-04,3.853e-04,4.422e-04,5.039e-04,4.069e-04,4.477e-04,4.963e-04,5.127e-04,6.463e-04,6.802e-04,7.662e-04,7.781e-04,1.107e-03,};
vector<float> sf_Vaeta_flatpt_v8_mumu_data16_metlt125_ee2j = {4.636e-01,4.968e-01,4.848e-01,4.790e-01,5.251e-01,5.153e-01,5.291e-01,5.723e-01,5.939e-01,5.943e-01,6.447e-01,7.892e-01,8.754e-01,9.071e-01,1.118e+00,1.166e+00,7.767e-01,6.652e-01,5.836e-01,4.975e-01,5.376e-01,4.739e-01,4.497e-01,4.078e-01,2.839e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v8_mumu_data17_metlt125_ee2j = {4.413e-01,4.397e-01,4.260e-01,4.366e-01,4.793e-01,4.638e-01,4.707e-01,5.149e-01,5.540e-01,5.401e-01,6.224e-01,7.653e-01,7.827e-01,7.940e-01,1.076e+00,1.041e+00,6.923e-01,6.404e-01,5.779e-01,5.508e-01,6.053e-01,5.841e-01,6.443e-01,6.286e-01,4.259e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v8_mumu_data18_metlt125_ee2j = {4.547e-01,4.490e-01,4.627e-01,4.596e-01,4.939e-01,4.640e-01,4.773e-01,5.224e-01,5.450e-01,5.759e-01,6.161e-01,7.599e-01,7.870e-01,8.187e-01,1.017e+00,1.022e+00,6.668e-01,6.671e-01,6.104e-01,5.370e-01,5.549e-01,5.756e-01,5.584e-01,4.454e-01,2.731e-01,1.000e+00,};

// vector<float> sf_Vpt_b0_v9_data16_metlt125_eq0j{sf_Vpt_b0_v8_data16_metlt125_eq0j}, sf_Vpt_b0_v9_ee_data16_metlt125_eq0j{sf_Vpt_b0_v8_ee_data16_metlt125_eq0j}, sf_Vpt_b0_v9_mumu_data16_metlt125_eq0j{sf_Vpt_b0_v8_mumu_data16_metlt125_eq0j};
// vector<float> sf_Vpt_b1_v9_data16_metlt125_eq1j{sf_Vpt_b1_v8_data16_metlt125_eq1j}, sf_Vpt_b1_v9_ee_data16_metlt125_eq1j{sf_Vpt_b1_v8_ee_data16_metlt125_eq1j}, sf_Vpt_b1_v9_mumu_data16_metlt125_eq1j{sf_Vpt_b1_v8_mumu_data16_metlt125_eq1j};
// vector<float> sf_Vpt_b0_v9_data17_metlt125_eq0j{sf_Vpt_b0_v8_data17_metlt125_eq0j}, sf_Vpt_b0_v9_ee_data17_metlt125_eq0j{sf_Vpt_b0_v8_ee_data17_metlt125_eq0j}, sf_Vpt_b0_v9_mumu_data17_metlt125_eq0j{sf_Vpt_b0_v8_mumu_data17_metlt125_eq0j};
// vector<float> sf_Vpt_b1_v9_data17_metlt125_eq1j{sf_Vpt_b1_v8_data17_metlt125_eq1j}, sf_Vpt_b1_v9_ee_data17_metlt125_eq1j{sf_Vpt_b1_v8_ee_data17_metlt125_eq1j}, sf_Vpt_b1_v9_mumu_data17_metlt125_eq1j{sf_Vpt_b1_v8_mumu_data17_metlt125_eq1j};
// vector<float> sf_Vpt_b0_v9_data18_metlt125_eq0j{sf_Vpt_b0_v8_data18_metlt125_eq0j}, sf_Vpt_b0_v9_ee_data18_metlt125_eq0j{sf_Vpt_b0_v8_ee_data18_metlt125_eq0j}, sf_Vpt_b0_v9_mumu_data18_metlt125_eq0j{sf_Vpt_b0_v8_mumu_data18_metlt125_eq0j};
// vector<float> sf_Vpt_b1_v9_data18_metlt125_eq1j{sf_Vpt_b1_v8_data18_metlt125_eq1j}, sf_Vpt_b1_v9_ee_data18_metlt125_eq1j{sf_Vpt_b1_v8_ee_data18_metlt125_eq1j}, sf_Vpt_b1_v9_mumu_data18_metlt125_eq1j{sf_Vpt_b1_v8_mumu_data18_metlt125_eq1j};
// vector<float> sferr_Vpt_b0_v9_data16_metlt125_eq0j{sferr_Vpt_b0_v8_data16_metlt125_eq0j}, sferr_Vpt_b0_v9_ee_data16_metlt125_eq0j{sferr_Vpt_b0_v8_ee_data16_metlt125_eq0j}, sferr_Vpt_b0_v9_mumu_data16_metlt125_eq0j{sferr_Vpt_b0_v8_mumu_data16_metlt125_eq0j};
// vector<float> sferr_Vpt_b1_v9_data16_metlt125_eq1j{sferr_Vpt_b1_v8_data16_metlt125_eq1j}, sferr_Vpt_b1_v9_ee_data16_metlt125_eq1j{sferr_Vpt_b1_v8_ee_data16_metlt125_eq1j}, sferr_Vpt_b1_v9_mumu_data16_metlt125_eq1j{sferr_Vpt_b1_v8_mumu_data16_metlt125_eq1j};
// vector<float> sferr_Vpt_b0_v9_data17_metlt125_eq0j{sferr_Vpt_b0_v8_data17_metlt125_eq0j}, sferr_Vpt_b0_v9_ee_data17_metlt125_eq0j{sferr_Vpt_b0_v8_ee_data17_metlt125_eq0j}, sferr_Vpt_b0_v9_mumu_data17_metlt125_eq0j{sferr_Vpt_b0_v8_mumu_data17_metlt125_eq0j};
// vector<float> sferr_Vpt_b1_v9_data17_metlt125_eq1j{sferr_Vpt_b1_v8_data17_metlt125_eq1j}, sferr_Vpt_b1_v9_ee_data17_metlt125_eq1j{sferr_Vpt_b1_v8_ee_data17_metlt125_eq1j}, sferr_Vpt_b1_v9_mumu_data17_metlt125_eq1j{sferr_Vpt_b1_v8_mumu_data17_metlt125_eq1j};
// vector<float> sferr_Vpt_b0_v9_data18_metlt125_eq0j{sferr_Vpt_b0_v8_data18_metlt125_eq0j}, sferr_Vpt_b0_v9_ee_data18_metlt125_eq0j{sferr_Vpt_b0_v8_ee_data18_metlt125_eq0j}, sferr_Vpt_b0_v9_mumu_data18_metlt125_eq0j{sferr_Vpt_b0_v8_mumu_data18_metlt125_eq0j};
// vector<float> sferr_Vpt_b1_v9_data18_metlt125_eq1j{sferr_Vpt_b1_v8_data18_metlt125_eq1j}, sferr_Vpt_b1_v9_ee_data18_metlt125_eq1j{sferr_Vpt_b1_v8_ee_data18_metlt125_eq1j}, sferr_Vpt_b1_v9_mumu_data18_metlt125_eq1j{sferr_Vpt_b1_v8_mumu_data18_metlt125_eq1j};

// vector<float> sf_Vpt_b1_v9_data16_metlt125_ee2j = {1.282e-02,1.497e-02,1.679e-02,1.959e-02,2.256e-02,2.376e-02,2.839e-02,2.952e-02,3.311e-02,3.804e-02,4.033e-02,4.352e-02,4.580e-02,5.153e-02,5.663e-02,6.455e-02,7.033e-02,7.961e-02,8.324e-02,9.403e-02,1.047e-01,1.379e-01,};
// vector<float> sf_Vpt_b1_v9_data17_metlt125_ee2j = {1.496e-02,1.707e-02,1.944e-02,2.285e-02,2.589e-02,2.879e-02,3.155e-02,3.463e-02,3.952e-02,4.167e-02,4.681e-02,4.753e-02,5.071e-02,5.773e-02,6.422e-02,7.041e-02,7.786e-02,8.123e-02,9.289e-02,1.058e-01,1.143e-01,1.454e-01,};
// vector<float> sf_Vpt_b1_v9_data18_metlt125_ee2j = {1.483e-02,1.715e-02,1.869e-02,2.277e-02,2.661e-02,2.797e-02,3.118e-02,3.425e-02,3.772e-02,4.144e-02,4.534e-02,4.995e-02,5.255e-02,5.621e-02,6.360e-02,6.759e-02,7.695e-02,8.480e-02,9.149e-02,1.042e-01,1.169e-01,1.346e-01,};
// vector<float> sferr_Vpt_b1_v9_data16_metlt125_ee2j = {2.084e-04,2.644e-04,3.219e-04,4.066e-04,5.116e-04,5.072e-04,4.904e-04,5.531e-04,6.452e-04,7.146e-04,8.080e-04,9.260e-04,1.034e-03,8.982e-04,1.065e-03,1.254e-03,1.269e-03,1.697e-03,1.867e-03,2.003e-03,2.231e-03,4.035e-03,};
// vector<float> sferr_Vpt_b1_v9_data17_metlt125_ee2j = {2.427e-04,3.097e-04,3.894e-04,5.077e-04,6.348e-04,6.529e-04,5.599e-04,6.649e-04,8.002e-04,7.799e-04,9.366e-04,1.025e-03,1.177e-03,1.044e-03,1.138e-03,1.320e-03,1.368e-03,1.727e-03,1.785e-03,2.069e-03,2.146e-03,3.678e-03,};
// vector<float> sferr_Vpt_b1_v9_data18_metlt125_ee2j = {2.500e-04,3.261e-04,3.860e-04,5.302e-04,6.972e-04,7.807e-04,9.664e-04,6.616e-04,7.667e-04,7.578e-04,8.861e-04,1.060e-03,1.183e-03,9.997e-04,1.028e-03,1.119e-03,1.189e-03,1.559e-03,1.573e-03,1.836e-03,1.817e-03,2.888e-03,};
// vector<float> sf_Vpt_b1_v9_ee_data16_metlt125_ee2j = {4.923e-03,5.850e-03,6.341e-03,7.439e-03,8.644e-03,9.169e-03,1.139e-02,1.121e-02,1.316e-02,1.514e-02,1.615e-02,1.655e-02,1.846e-02,2.065e-02,2.267e-02,2.559e-02,2.937e-02,3.126e-02,3.293e-02,3.835e-02,4.241e-02,5.486e-02,};
// vector<float> sf_Vpt_b1_v9_ee_data17_metlt125_ee2j = {5.885e-03,6.640e-03,7.696e-03,9.018e-03,1.020e-02,1.136e-02,1.253e-02,1.375e-02,1.613e-02,1.669e-02,1.890e-02,1.956e-02,2.105e-02,2.325e-02,2.633e-02,2.899e-02,3.268e-02,3.381e-02,3.717e-02,4.447e-02,5.000e-02,6.286e-02,};
// vector<float> sf_Vpt_b1_v9_ee_data18_metlt125_ee2j = {5.899e-03,6.924e-03,7.487e-03,9.238e-03,1.093e-02,1.165e-02,1.251e-02,1.423e-02,1.527e-02,1.687e-02,1.930e-02,2.134e-02,2.178e-02,2.340e-02,2.601e-02,2.870e-02,3.257e-02,3.665e-02,3.837e-02,4.469e-02,5.024e-02,6.099e-02,};
// vector<float> sferr_Vpt_b1_v9_ee_data16_metlt125_ee2j = {1.058e-04,1.333e-04,1.565e-04,1.947e-04,2.428e-04,2.528e-04,2.778e-04,3.021e-04,3.615e-04,4.152e-04,4.690e-04,5.193e-04,5.968e-04,5.123e-04,6.198e-04,7.503e-04,7.776e-04,1.000e-03,1.102e-03,1.247e-03,1.379e-03,2.450e-03,};
// vector<float> sferr_Vpt_b1_v9_ee_data17_metlt125_ee2j = {1.147e-04,1.431e-04,1.796e-04,2.306e-04,2.846e-04,3.035e-04,2.892e-04,3.394e-04,4.152e-04,4.228e-04,5.044e-04,5.580e-04,6.393e-04,5.487e-04,6.417e-04,7.656e-04,7.960e-04,9.959e-04,1.069e-03,1.277e-03,1.378e-03,2.329e-03,};
// vector<float> sferr_Vpt_b1_v9_ee_data18_metlt125_ee2j = {1.141e-04,1.486e-04,1.736e-04,2.372e-04,3.117e-04,3.521e-04,4.187e-04,3.302e-04,3.755e-04,3.938e-04,4.701e-04,5.581e-04,6.062e-04,5.103e-04,5.607e-04,6.509e-04,6.824e-04,8.988e-04,9.363e-04,1.106e-03,1.155e-03,1.880e-03,};
// vector<float> sf_Vpt_b1_v9_mumu_data16_metlt125_ee2j = {7.900e-03,9.124e-03,1.044e-02,1.216e-02,1.392e-02,1.460e-02,1.700e-02,1.831e-02,1.995e-02,2.290e-02,2.418e-02,2.698e-02,2.735e-02,3.088e-02,3.396e-02,3.896e-02,4.096e-02,4.835e-02,5.031e-02,5.568e-02,6.232e-02,8.307e-02,};
// vector<float> sf_Vpt_b1_v9_mumu_data17_metlt125_ee2j = {9.077e-03,1.043e-02,1.174e-02,1.383e-02,1.568e-02,1.743e-02,1.901e-02,2.089e-02,2.339e-02,2.498e-02,2.791e-02,2.797e-02,2.966e-02,3.448e-02,3.789e-02,4.141e-02,4.518e-02,4.742e-02,5.571e-02,6.137e-02,6.425e-02,8.259e-02,};
// vector<float> sf_Vpt_b1_v9_mumu_data18_metlt125_ee2j = {8.931e-03,1.023e-02,1.120e-02,1.353e-02,1.568e-02,1.632e-02,1.868e-02,2.002e-02,2.245e-02,2.457e-02,2.604e-02,2.861e-02,3.077e-02,3.281e-02,3.758e-02,3.889e-02,4.438e-02,4.815e-02,5.312e-02,5.954e-02,6.665e-02,7.361e-02,};
// vector<float> sferr_Vpt_b1_v9_mumu_data16_metlt125_ee2j = {1.459e-04,1.818e-04,2.232e-04,2.788e-04,3.467e-04,3.502e-04,3.532e-04,4.059e-04,4.649e-04,5.258e-04,5.916e-04,6.893e-04,7.509e-04,6.502e-04,7.811e-04,9.427e-04,9.326e-04,1.272e-03,1.393e-03,1.514e-03,1.687e-03,3.055e-03,};
// vector<float> sferr_Vpt_b1_v9_mumu_data17_metlt125_ee2j = {1.604e-04,2.044e-04,2.526e-04,3.278e-04,4.078e-04,4.265e-04,3.848e-04,4.542e-04,5.387e-04,5.476e-04,6.506e-04,7.052e-04,8.024e-04,7.163e-04,8.032e-04,9.452e-04,9.666e-04,1.221e-03,1.334e-03,1.521e-03,1.573e-03,2.695e-03,};
// vector<float> sferr_Vpt_b1_v9_mumu_data18_metlt125_ee2j = {1.607e-04,2.064e-04,2.445e-04,3.306e-04,4.289e-04,4.753e-04,6.002e-04,4.279e-04,5.028e-04,5.117e-04,5.812e-04,6.891e-04,7.790e-04,6.540e-04,7.117e-04,7.828e-04,8.262e-04,1.066e-03,1.129e-03,1.305e-03,1.341e-03,2.077e-03,};

vector<float> sf_Vaeta_flatpt_v9_data16_metlt125_ee2j = {7.827e-01,8.133e-01,8.080e-01,8.014e-01,8.692e-01,8.689e-01,8.785e-01,9.533e-01,9.695e-01,9.496e-01,1.080e+00,1.250e+00,1.375e+00,1.471e+00,1.822e+00,1.893e+00,1.271e+00,1.075e+00,9.360e-01,8.339e-01,8.564e-01,7.880e-01,7.714e-01,6.801e-01,4.801e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v9_data17_metlt125_ee2j = {7.503e-01,7.596e-01,7.362e-01,7.699e-01,7.891e-01,8.022e-01,8.035e-01,9.013e-01,8.894e-01,9.259e-01,1.064e+00,1.265e+00,1.255e+00,1.301e+00,1.833e+00,1.708e+00,1.149e+00,1.041e+00,9.532e-01,8.691e-01,9.990e-01,9.617e-01,1.049e+00,1.098e+00,7.724e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v9_data18_metlt125_ee2j = {7.542e-01,7.571e-01,7.472e-01,7.919e-01,8.350e-01,7.673e-01,8.160e-01,8.718e-01,9.094e-01,9.486e-01,1.023e+00,1.220e+00,1.317e+00,1.356e+00,1.674e+00,1.751e+00,1.146e+00,1.122e+00,9.980e-01,9.377e-01,9.453e-01,9.854e-01,1.038e+00,8.251e-01,5.823e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v9_ee_data16_metlt125_ee2j = {8.356e-01,8.478e-01,8.527e-01,8.232e-01,8.723e-01,8.783e-01,9.059e-01,9.887e-01,9.697e-01,9.260e-01,1.114e+00,1.223e+00,1.355e+00,1.388e+00,1.767e+00,1.795e+00,1.217e+00,1.047e+00,9.291e-01,7.947e-01,8.047e-01,7.949e-01,7.975e-01,6.817e-01,4.977e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v9_ee_data17_metlt125_ee2j = {7.542e-01,7.806e-01,7.406e-01,8.232e-01,8.089e-01,8.217e-01,8.052e-01,9.309e-01,8.550e-01,9.619e-01,1.058e+00,1.263e+00,1.212e+00,1.286e+00,1.741e+00,1.710e+00,1.079e+00,1.010e+00,9.382e-01,8.496e-01,9.740e-01,9.440e-01,1.030e+00,1.084e+00,8.057e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v9_ee_data18_metlt125_ee2j = {7.478e-01,7.399e-01,7.245e-01,8.070e-01,8.323e-01,7.637e-01,8.171e-01,8.386e-01,8.725e-01,9.214e-01,9.773e-01,1.160e+00,1.231e+00,1.280e+00,1.574e+00,1.680e+00,1.134e+00,1.093e+00,1.001e+00,9.593e-01,9.744e-01,1.009e+00,1.165e+00,8.925e-01,6.785e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v9_mumu_data16_metlt125_ee2j = {7.486e-01,7.911e-01,7.791e-01,7.873e-01,8.673e-01,8.628e-01,8.608e-01,9.304e-01,9.693e-01,9.648e-01,1.058e+00,1.267e+00,1.389e+00,1.525e+00,1.859e+00,1.957e+00,1.306e+00,1.092e+00,9.404e-01,8.592e-01,8.898e-01,7.836e-01,7.546e-01,6.791e-01,4.689e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v9_mumu_data17_metlt125_ee2j = {7.477e-01,7.455e-01,7.332e-01,7.344e-01,7.759e-01,7.891e-01,8.024e-01,8.815e-01,9.123e-01,9.019e-01,1.067e+00,1.266e+00,1.284e+00,1.311e+00,1.895e+00,1.707e+00,1.195e+00,1.062e+00,9.632e-01,8.821e-01,1.016e+00,9.735e-01,1.061e+00,1.107e+00,7.503e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v9_mumu_data18_metlt125_ee2j = {7.586e-01,7.690e-01,7.629e-01,7.813e-01,8.367e-01,7.697e-01,8.151e-01,8.949e-01,9.351e-01,9.675e-01,1.054e+00,1.262e+00,1.377e+00,1.409e+00,1.743e+00,1.801e+00,1.155e+00,1.143e+00,9.956e-01,9.226e-01,9.251e-01,9.689e-01,9.503e-01,7.783e-01,5.152e-01,1.000e+00,};

vector<float> sf_Vpt_b0_v9_data16_metlt125_eq0j = {1.686e-02,1.847e-02,2.091e-02,2.232e-02,2.487e-02,2.696e-02,3.363e-02,3.604e-02,3.432e-02,3.207e-02,};
vector<float> sf_Vpt_b1_v9_data16_metlt125_eq1j = {1.878e-02,2.160e-02,2.348e-02,2.728e-02,2.936e-02,3.306e-02,3.640e-02,3.978e-02,4.305e-02,4.815e-02,4.893e-02,5.270e-02,5.474e-02,5.947e-02,6.338e-02,7.020e-02,7.436e-02,7.900e-02,8.668e-02,8.850e-02,9.578e-02,1.175e-01,};
vector<float> sf_Vpt_b1_v9_data16_metlt125_ee2j = {9.960e-03,1.163e-02,1.289e-02,1.502e-02,1.737e-02,1.804e-02,2.133e-02,2.182e-02,2.418e-02,2.798e-02,2.895e-02,3.206e-02,3.302e-02,3.618e-02,3.911e-02,4.490e-02,4.785e-02,5.300e-02,5.386e-02,6.185e-02,6.736e-02,8.002e-02,};
vector<float> sf_Vpt_b0_v9_data17_metlt125_eq0j = {1.778e-02,1.938e-02,2.177e-02,2.615e-02,2.588e-02,3.561e-02,4.583e-02,6.919e-02,1.001e-01,1.357e-01,};
vector<float> sf_Vpt_b1_v9_data17_metlt125_eq1j = {1.996e-02,2.336e-02,2.621e-02,2.981e-02,3.249e-02,3.583e-02,4.004e-02,4.226e-02,4.512e-02,4.913e-02,5.127e-02,5.467e-02,5.885e-02,6.235e-02,6.655e-02,7.298e-02,8.074e-02,8.774e-02,8.993e-02,9.349e-02,1.044e-01,1.108e-01,};
vector<float> sf_Vpt_b1_v9_data17_metlt125_ee2j = {1.194e-02,1.351e-02,1.531e-02,1.789e-02,2.010e-02,2.206e-02,2.398e-02,2.580e-02,2.910e-02,3.069e-02,3.454e-02,3.454e-02,3.649e-02,4.172e-02,4.524e-02,4.853e-02,5.356e-02,5.464e-02,6.225e-02,7.037e-02,7.384e-02,8.860e-02,};
vector<float> sf_Vpt_b0_v9_data18_metlt125_eq0j = {1.836e-02,2.063e-02,2.275e-02,2.640e-02,2.922e-02,3.370e-02,4.593e-02,7.491e-02,1.003e-01,1.425e-01,};
vector<float> sf_Vpt_b1_v9_data18_metlt125_eq1j = {2.076e-02,2.323e-02,2.678e-02,2.925e-02,3.278e-02,3.567e-02,3.851e-02,4.304e-02,4.658e-02,5.130e-02,5.250e-02,5.574e-02,5.822e-02,6.279e-02,6.829e-02,7.193e-02,7.935e-02,8.167e-02,8.956e-02,9.427e-02,9.834e-02,9.764e-02,};
vector<float> sf_Vpt_b1_v9_data18_metlt125_ee2j = {1.156e-02,1.328e-02,1.433e-02,1.741e-02,2.036e-02,2.107e-02,2.317e-02,2.504e-02,2.780e-02,3.012e-02,3.255e-02,3.611e-02,3.693e-02,3.923e-02,4.388e-02,4.583e-02,5.222e-02,5.639e-02,6.034e-02,6.560e-02,7.163e-02,7.672e-02,};
vector<float> sferr_Vpt_b0_v9_data16_metlt125_eq0j = {1.820e-04,2.598e-04,3.835e-04,5.295e-04,7.475e-04,6.378e-04,8.067e-04,1.697e-03,2.478e-03,4.991e-03,};
vector<float> sferr_Vpt_b1_v9_data16_metlt125_eq1j = {1.618e-04,2.138e-04,2.648e-04,3.516e-04,4.272e-04,4.684e-04,4.232e-04,5.141e-04,6.125e-04,6.778e-04,7.603e-04,9.004e-04,1.029e-03,8.942e-04,1.077e-03,1.266e-03,1.321e-03,1.744e-03,2.052e-03,2.158e-03,2.479e-03,4.600e-03,};
vector<float> sferr_Vpt_b1_v9_data16_metlt125_ee2j = {1.726e-04,2.177e-04,2.620e-04,3.289e-04,4.131e-04,4.104e-04,4.076e-04,4.539e-04,5.253e-04,5.928e-04,6.591e-04,7.659e-04,8.427e-04,7.177e-04,8.502e-04,1.021e-03,1.018e-03,1.341e-03,1.451e-03,1.602e-03,1.760e-03,2.997e-03,};
vector<float> sferr_Vpt_b0_v9_data17_metlt125_eq0j = {2.443e-04,3.533e-04,5.270e-04,8.634e-04,1.059e-03,1.173e-03,1.567e-03,4.902e-03,9.451e-03,2.089e-02,};
vector<float> sferr_Vpt_b1_v9_data17_metlt125_eq1j = {1.999e-04,2.736e-04,3.549e-04,4.679e-04,5.820e-04,6.157e-04,5.491e-04,6.496e-04,7.552e-04,7.794e-04,9.039e-04,1.061e-03,1.261e-03,1.081e-03,1.166e-03,1.402e-03,1.525e-03,2.030e-03,1.995e-03,2.266e-03,2.529e-03,4.127e-03,};
vector<float> sferr_Vpt_b1_v9_data17_metlt125_ee2j = {2.006e-04,2.534e-04,3.162e-04,4.090e-04,5.062e-04,5.193e-04,4.551e-04,5.305e-04,6.325e-04,6.290e-04,7.531e-04,8.151e-04,9.265e-04,8.202e-04,9.000e-04,1.042e-03,1.074e-03,1.334e-03,1.419e-03,1.641e-03,1.694e-03,2.799e-03,};
vector<float> sferr_Vpt_b0_v9_data18_metlt125_eq0j = {2.390e-04,3.606e-04,5.269e-04,8.057e-04,1.188e-03,1.388e-03,1.398e-03,4.456e-03,7.415e-03,1.671e-02,};
vector<float> sferr_Vpt_b1_v9_data18_metlt125_eq1j = {2.065e-04,2.670e-04,3.612e-04,4.529e-04,5.878e-04,7.263e-04,8.937e-04,6.386e-04,7.511e-04,7.587e-04,8.618e-04,1.008e-03,1.168e-03,1.031e-03,1.037e-03,1.159e-03,1.255e-03,1.600e-03,1.698e-03,1.948e-03,1.918e-03,2.926e-03,};
vector<float> sferr_Vpt_b1_v9_data18_metlt125_ee2j = {2.005e-04,2.593e-04,3.037e-04,4.146e-04,5.438e-04,6.001e-04,7.317e-04,5.109e-04,5.960e-04,5.939e-04,6.856e-04,8.211e-04,8.946e-04,7.500e-04,7.903e-04,8.677e-04,9.169e-04,1.182e-03,1.218e-03,1.382e-03,1.393e-03,2.124e-03,};

vector<float> sf_Vpt_b0_v9_ee_data16_metlt125_eq0j = {6.354e-03,7.071e-03,8.124e-03,8.692e-03,9.542e-03,1.051e-02,1.340e-02,1.483e-02,1.539e-02,1.603e-02,};
vector<float> sf_Vpt_b1_v9_ee_data16_metlt125_eq1j = {7.152e-03,8.297e-03,8.961e-03,1.047e-02,1.115e-02,1.274e-02,1.411e-02,1.554e-02,1.734e-02,1.865e-02,1.975e-02,2.085e-02,2.197e-02,2.385e-02,2.467e-02,2.833e-02,2.917e-02,3.203e-02,3.535e-02,3.569e-02,3.831e-02,4.512e-02,};
vector<float> sf_Vpt_b1_v9_ee_data16_metlt125_ee2j = {3.867e-03,4.488e-03,4.881e-03,5.785e-03,6.645e-03,6.987e-03,8.640e-03,8.526e-03,9.403e-03,1.115e-02,1.152e-02,1.216e-02,1.327e-02,1.427e-02,1.567e-02,1.768e-02,1.996e-02,2.123e-02,2.152e-02,2.526e-02,2.886e-02,3.138e-02,};
vector<float> sf_Vpt_b0_v9_ee_data17_metlt125_eq0j = {6.968e-03,7.570e-03,8.479e-03,1.018e-02,9.934e-03,1.432e-02,1.755e-02,2.558e-02,3.272e-02,4.524e-02,};
vector<float> sf_Vpt_b1_v9_ee_data17_metlt125_eq1j = {7.829e-03,9.171e-03,1.040e-02,1.177e-02,1.294e-02,1.411e-02,1.633e-02,1.708e-02,1.844e-02,1.964e-02,2.087e-02,2.233e-02,2.415e-02,2.565e-02,2.754e-02,3.081e-02,3.187e-02,3.747e-02,3.634e-02,3.974e-02,4.278e-02,4.898e-02,};
vector<float> sf_Vpt_b1_v9_ee_data17_metlt125_ee2j = {4.702e-03,5.258e-03,5.993e-03,7.094e-03,7.961e-03,8.726e-03,9.469e-03,1.005e-02,1.187e-02,1.209e-02,1.372e-02,1.450e-02,1.482e-02,1.678e-02,1.845e-02,2.022e-02,2.277e-02,2.187e-02,2.524e-02,2.974e-02,3.187e-02,3.695e-02,};
vector<float> sf_Vpt_b0_v9_ee_data18_metlt125_eq0j = {7.287e-03,8.227e-03,8.928e-03,1.073e-02,1.172e-02,1.329e-02,1.901e-02,2.860e-02,4.393e-02,6.405e-02,};
vector<float> sf_Vpt_b1_v9_ee_data18_metlt125_eq1j = {8.259e-03,9.368e-03,1.069e-02,1.179e-02,1.324e-02,1.447e-02,1.561e-02,1.734e-02,1.897e-02,2.136e-02,2.209e-02,2.321e-02,2.415e-02,2.658e-02,2.852e-02,3.034e-02,3.361e-02,3.508e-02,3.807e-02,3.963e-02,4.229e-02,4.555e-02,};
vector<float> sf_Vpt_b1_v9_ee_data18_metlt125_ee2j = {4.599e-03,5.358e-03,5.725e-03,7.018e-03,8.328e-03,8.852e-03,9.336e-03,1.038e-02,1.131e-02,1.231e-02,1.382e-02,1.533e-02,1.509e-02,1.608e-02,1.804e-02,1.932e-02,2.198e-02,2.441e-02,2.556e-02,2.793e-02,3.082e-02,3.511e-02,};
vector<float> sferr_Vpt_b0_v9_ee_data16_metlt125_eq0j = {8.775e-05,1.252e-04,1.841e-04,2.525e-04,3.481e-04,3.246e-04,4.520e-04,1.002e-03,1.594e-03,3.474e-03,};
vector<float> sferr_Vpt_b1_v9_ee_data16_metlt125_eq1j = {7.793e-05,1.019e-04,1.243e-04,1.631e-04,1.949e-04,2.245e-04,2.289e-04,2.771e-04,3.369e-04,3.810e-04,4.362e-04,5.077e-04,5.832e-04,5.030e-04,6.098e-04,7.610e-04,7.795e-04,1.044e-03,1.226e-03,1.337e-03,1.526e-03,2.757e-03,};
vector<float> sferr_Vpt_b1_v9_ee_data16_metlt125_ee2j = {9.075e-05,1.121e-04,1.316e-04,1.641e-04,2.020e-04,2.116e-04,2.373e-04,2.582e-04,2.981e-04,3.512e-04,3.893e-04,4.380e-04,4.964e-04,4.162e-04,5.062e-04,6.174e-04,6.333e-04,8.138e-04,8.785e-04,1.006e-03,1.131e-03,1.834e-03,};
vector<float> sferr_Vpt_b0_v9_ee_data17_metlt125_eq0j = {1.127e-04,1.614e-04,2.369e-04,3.811e-04,4.638e-04,5.572e-04,7.650e-04,2.373e-03,4.488e-03,1.097e-02,};
vector<float> sferr_Vpt_b1_v9_ee_data17_metlt125_eq1j = {9.160e-05,1.234e-04,1.599e-04,2.075e-04,2.582e-04,2.790e-04,2.811e-04,3.282e-04,3.856e-04,4.148e-04,4.845e-04,5.657e-04,6.697e-04,5.718e-04,6.586e-04,8.225e-04,8.525e-04,1.181e-03,1.203e-03,1.415e-03,1.573e-03,2.666e-03,};
vector<float> sferr_Vpt_b1_v9_ee_data17_metlt125_ee2j = {9.732e-05,1.202e-04,1.482e-04,1.910e-04,2.335e-04,2.488e-04,2.413e-04,2.764e-04,3.391e-04,3.476e-04,4.138e-04,4.631e-04,5.132e-04,4.458e-04,5.207e-04,6.239e-04,6.468e-04,7.752e-04,8.704e-04,1.031e-03,1.091e-03,1.764e-03,};
vector<float> sferr_Vpt_b0_v9_ee_data18_metlt125_eq0j = {1.068e-04,1.603e-04,2.291e-04,3.572e-04,5.170e-04,5.897e-04,6.872e-04,2.069e-03,4.084e-03,1.012e-02,};
vector<float> sferr_Vpt_b1_v9_ee_data18_metlt125_eq1j = {9.154e-05,1.187e-04,1.575e-04,1.979e-04,2.554e-04,3.154e-04,3.860e-04,3.042e-04,3.613e-04,3.912e-04,4.466e-04,5.153e-04,5.928e-04,5.264e-04,5.676e-04,6.680e-04,7.191e-04,9.221e-04,1.021e-03,1.168e-03,1.225e-03,1.951e-03,};
vector<float> sferr_Vpt_b1_v9_ee_data18_metlt125_ee2j = {9.389e-05,1.210e-04,1.396e-04,1.886e-04,2.471e-04,2.781e-04,3.250e-04,2.634e-04,3.028e-04,3.202e-04,3.758e-04,4.452e-04,4.719e-04,3.938e-04,4.481e-04,5.173e-04,5.412e-04,7.050e-04,7.472e-04,8.517e-04,8.963e-04,1.409e-03,};

vector<float> sf_Vpt_b0_v9_mumu_data16_metlt125_eq0j = {1.051e-02,1.140e-02,1.279e-02,1.363e-02,1.533e-02,1.645e-02,2.023e-02,2.121e-02,1.893e-02,1.603e-02,};
vector<float> sf_Vpt_b1_v9_mumu_data16_metlt125_eq1j = {1.163e-02,1.330e-02,1.452e-02,1.681e-02,1.821e-02,2.032e-02,2.228e-02,2.424e-02,2.571e-02,2.950e-02,2.918e-02,3.185e-02,3.277e-02,3.562e-02,3.871e-02,4.188e-02,4.519e-02,4.697e-02,5.132e-02,5.281e-02,5.747e-02,7.236e-02,};
vector<float> sf_Vpt_b1_v9_mumu_data16_metlt125_ee2j = {6.093e-03,7.137e-03,8.010e-03,9.237e-03,1.072e-02,1.105e-02,1.269e-02,1.329e-02,1.478e-02,1.683e-02,1.743e-02,1.990e-02,1.974e-02,2.192e-02,2.345e-02,2.722e-02,2.789e-02,3.177e-02,3.234e-02,3.660e-02,3.851e-02,4.863e-02,};
vector<float> sf_Vpt_b0_v9_mumu_data17_metlt125_eq0j = {1.081e-02,1.181e-02,1.329e-02,1.597e-02,1.594e-02,2.129e-02,2.828e-02,4.360e-02,6.737e-02,9.048e-02,};
vector<float> sf_Vpt_b1_v9_mumu_data17_metlt125_eq1j = {1.213e-02,1.419e-02,1.581e-02,1.805e-02,1.956e-02,2.172e-02,2.371e-02,2.518e-02,2.668e-02,2.949e-02,3.040e-02,3.235e-02,3.470e-02,3.670e-02,3.901e-02,4.217e-02,4.888e-02,5.027e-02,5.359e-02,5.375e-02,6.158e-02,6.185e-02,};
vector<float> sf_Vpt_b1_v9_mumu_data17_metlt125_ee2j = {7.240e-03,8.248e-03,9.316e-03,1.080e-02,1.214e-02,1.334e-02,1.451e-02,1.575e-02,1.723e-02,1.860e-02,2.082e-02,2.005e-02,2.167e-02,2.494e-02,2.679e-02,2.831e-02,3.079e-02,3.277e-02,3.700e-02,4.063e-02,4.197e-02,5.165e-02,};
vector<float> sf_Vpt_b0_v9_mumu_data18_metlt125_eq0j = {1.107e-02,1.240e-02,1.383e-02,1.567e-02,1.749e-02,2.041e-02,2.692e-02,4.631e-02,5.634e-02,7.843e-02,};
vector<float> sf_Vpt_b1_v9_mumu_data18_metlt125_eq1j = {1.250e-02,1.386e-02,1.609e-02,1.746e-02,1.955e-02,2.120e-02,2.291e-02,2.570e-02,2.761e-02,2.994e-02,3.041e-02,3.254e-02,3.407e-02,3.621e-02,3.977e-02,4.160e-02,4.574e-02,4.659e-02,5.149e-02,5.464e-02,5.605e-02,5.210e-02,};
vector<float> sf_Vpt_b1_v9_mumu_data18_metlt125_ee2j = {6.958e-03,7.926e-03,8.606e-03,1.039e-02,1.203e-02,1.222e-02,1.383e-02,1.466e-02,1.650e-02,1.781e-02,1.873e-02,2.078e-02,2.185e-02,2.316e-02,2.584e-02,2.652e-02,3.025e-02,3.198e-02,3.478e-02,3.767e-02,4.081e-02,4.161e-02,};
vector<float> sferr_Vpt_b0_v9_mumu_data16_metlt125_eq0j = {1.259e-04,1.775e-04,2.582e-04,3.547e-04,5.011e-04,4.414e-04,5.800e-04,1.230e-03,1.782e-03,3.474e-03,};
vector<float> sferr_Vpt_b1_v9_mumu_data16_metlt125_eq1j = {1.110e-04,1.448e-04,1.790e-04,2.354e-04,2.861e-04,3.174e-04,3.043e-04,3.669e-04,4.317e-04,4.986e-04,5.492e-04,6.535e-04,7.408e-04,6.413e-04,7.928e-04,9.427e-04,9.919e-04,1.290e-03,1.510e-03,1.640e-03,1.887e-03,3.536e-03,};
vector<float> sferr_Vpt_b1_v9_mumu_data16_metlt125_ee2j = {1.220e-04,1.528e-04,1.843e-04,2.276e-04,2.844e-04,2.876e-04,2.964e-04,3.342e-04,3.875e-04,4.412e-04,4.900e-04,5.775e-04,6.210e-04,5.311e-04,6.325e-04,7.764e-04,7.569e-04,1.010e-03,1.093e-03,1.218e-03,1.312e-03,2.302e-03,};
vector<float> sferr_Vpt_b0_v9_mumu_data17_metlt125_eq0j = {1.600e-04,2.310e-04,3.428e-04,5.571e-04,6.896e-04,7.611e-04,1.077e-03,3.448e-03,7.146e-03,1.630e-02,};
vector<float> sferr_Vpt_b1_v9_mumu_data17_metlt125_eq1j = {1.304e-04,1.770e-04,2.270e-04,2.987e-04,3.683e-04,3.980e-04,3.669e-04,4.344e-04,5.028e-04,5.421e-04,6.219e-04,7.254e-04,8.565e-04,7.309e-04,8.172e-04,9.912e-04,1.103e-03,1.413e-03,1.487e-03,1.665e-03,1.905e-03,3.015e-03,};
vector<float> sferr_Vpt_b1_v9_mumu_data17_metlt125_ee2j = {1.343e-04,1.693e-04,2.090e-04,2.667e-04,3.282e-04,3.439e-04,3.190e-04,3.723e-04,4.343e-04,4.526e-04,5.365e-04,5.668e-04,6.514e-04,5.747e-04,6.485e-04,7.551e-04,7.688e-04,9.779e-04,1.067e-03,1.217e-03,1.258e-03,2.100e-03,};
vector<float> sferr_Vpt_b0_v9_mumu_data18_metlt125_eq0j = {1.523e-04,2.281e-04,3.350e-04,4.990e-04,7.387e-04,8.688e-04,8.993e-04,2.996e-03,4.846e-03,1.143e-02,};
vector<float> sferr_Vpt_b1_v9_mumu_data18_metlt125_eq1j = {1.307e-04,1.669e-04,2.261e-04,2.809e-04,3.629e-04,4.460e-04,5.480e-04,4.143e-04,4.849e-04,4.994e-04,5.631e-04,6.600e-04,7.642e-04,6.634e-04,7.066e-04,8.103e-04,8.707e-04,1.101e-03,1.214e-03,1.403e-03,1.420e-03,2.093e-03,};
vector<float> sferr_Vpt_b1_v9_mumu_data18_metlt125_ee2j = {1.305e-04,1.662e-04,1.950e-04,2.625e-04,3.389e-04,3.673e-04,4.578e-04,3.375e-04,3.975e-04,4.086e-04,4.607e-04,5.478e-04,6.076e-04,5.066e-04,5.583e-04,6.212e-04,6.527e-04,8.272e-04,8.859e-04,1.005e-03,1.036e-03,1.538e-03,};

vector<float> sf_Vaeta_flatpt_v10_data16_metlt125_ee2j = {7.806e-01,8.104e-01,8.068e-01,7.990e-01,8.664e-01,8.667e-01,8.757e-01,9.503e-01,9.675e-01,9.469e-01,1.078e+00,1.245e+00,1.373e+00,1.467e+00,1.821e+00,1.886e+00,1.264e+00,1.071e+00,9.335e-01,8.316e-01,8.553e-01,7.865e-01,7.692e-01,6.796e-01,4.789e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v10_data17_metlt125_ee2j = {7.484e-01,7.584e-01,7.340e-01,7.670e-01,7.874e-01,8.006e-01,8.014e-01,8.994e-01,8.873e-01,9.242e-01,1.061e+00,1.262e+00,1.251e+00,1.298e+00,1.828e+00,1.703e+00,1.147e+00,1.038e+00,9.510e-01,8.672e-01,9.962e-01,9.589e-01,1.045e+00,1.096e+00,7.711e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v10_data18_metlt125_ee2j = {7.525e-01,7.557e-01,7.450e-01,7.903e-01,8.328e-01,7.644e-01,8.137e-01,8.693e-01,9.076e-01,9.448e-01,1.020e+00,1.217e+00,1.313e+00,1.354e+00,1.670e+00,1.747e+00,1.143e+00,1.119e+00,9.954e-01,9.349e-01,9.424e-01,9.833e-01,1.036e+00,8.235e-01,5.809e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v10_ee_data16_metlt125_ee2j = {8.325e-01,8.424e-01,8.512e-01,8.202e-01,8.699e-01,8.744e-01,9.044e-01,9.837e-01,9.664e-01,9.244e-01,1.114e+00,1.217e+00,1.353e+00,1.381e+00,1.764e+00,1.790e+00,1.211e+00,1.044e+00,9.263e-01,7.930e-01,8.028e-01,7.949e-01,7.953e-01,6.817e-01,4.955e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v10_ee_data17_metlt125_ee2j = {7.523e-01,7.792e-01,7.392e-01,8.193e-01,8.069e-01,8.203e-01,8.023e-01,9.288e-01,8.528e-01,9.608e-01,1.055e+00,1.261e+00,1.210e+00,1.284e+00,1.736e+00,1.705e+00,1.077e+00,1.009e+00,9.382e-01,8.490e-01,9.697e-01,9.408e-01,1.026e+00,1.081e+00,8.034e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v10_ee_data18_metlt125_ee2j = {7.463e-01,7.388e-01,7.216e-01,8.055e-01,8.299e-01,7.612e-01,8.153e-01,8.358e-01,8.692e-01,9.185e-01,9.738e-01,1.156e+00,1.229e+00,1.277e+00,1.570e+00,1.679e+00,1.132e+00,1.088e+00,9.977e-01,9.553e-01,9.729e-01,1.007e+00,1.163e+00,8.898e-01,6.768e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v10_mumu_data16_metlt125_ee2j = {7.471e-01,7.897e-01,7.781e-01,7.854e-01,8.641e-01,8.618e-01,8.572e-01,9.288e-01,9.683e-01,9.615e-01,1.055e+00,1.264e+00,1.386e+00,1.523e+00,1.858e+00,1.948e+00,1.299e+00,1.089e+00,9.381e-01,8.566e-01,8.892e-01,7.811e-01,7.524e-01,6.783e-01,4.682e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v10_mumu_data17_metlt125_ee2j = {7.458e-01,7.445e-01,7.304e-01,7.321e-01,7.743e-01,7.875e-01,8.008e-01,8.798e-01,9.102e-01,8.997e-01,1.065e+00,1.263e+00,1.279e+00,1.307e+00,1.889e+00,1.702e+00,1.193e+00,1.057e+00,9.595e-01,8.794e-01,1.014e+00,9.708e-01,1.058e+00,1.106e+00,7.496e-01,1.000e+00,};
vector<float> sf_Vaeta_flatpt_v10_mumu_data18_metlt125_ee2j = {7.567e-01,7.674e-01,7.614e-01,7.797e-01,8.348e-01,7.666e-01,8.125e-01,8.926e-01,9.343e-01,9.632e-01,1.052e+00,1.259e+00,1.372e+00,1.407e+00,1.740e+00,1.795e+00,1.152e+00,1.140e+00,9.937e-01,9.207e-01,9.212e-01,9.665e-01,9.483e-01,7.773e-01,5.141e-01,1.000e+00,};

vector<float> sf_Vpt_b0_v10_data16_metlt125_eq0j = {1.679e-02,1.837e-02,2.080e-02,2.222e-02,2.471e-02,2.680e-02,3.339e-02,3.592e-02,3.417e-02,3.134e-02,};
vector<float> sf_Vpt_b1_v10_data16_metlt125_eq1j = {1.874e-02,2.155e-02,2.342e-02,2.722e-02,2.931e-02,3.299e-02,3.628e-02,3.966e-02,4.298e-02,4.802e-02,4.882e-02,5.262e-02,5.463e-02,5.930e-02,6.326e-02,7.011e-02,7.410e-02,7.888e-02,8.645e-02,8.811e-02,9.554e-02,1.170e-01,};
vector<float> sf_Vpt_b1_v10_data16_metlt125_ee2j = {9.933e-03,1.160e-02,1.286e-02,1.496e-02,1.732e-02,1.800e-02,2.129e-02,2.181e-02,2.413e-02,2.790e-02,2.886e-02,3.195e-02,3.290e-02,3.612e-02,3.902e-02,4.476e-02,4.773e-02,5.288e-02,5.359e-02,6.174e-02,6.719e-02,7.970e-02,};
vector<float> sf_Vpt_b0_v10_data17_metlt125_eq0j = {1.776e-02,1.936e-02,2.176e-02,2.627e-02,2.578e-02,3.559e-02,4.570e-02,6.933e-02,1.004e-01,1.370e-01,};
vector<float> sf_Vpt_b1_v10_data17_metlt125_eq1j = {1.996e-02,2.334e-02,2.618e-02,2.982e-02,3.251e-02,3.583e-02,4.005e-02,4.220e-02,4.510e-02,4.911e-02,5.139e-02,5.461e-02,5.881e-02,6.239e-02,6.660e-02,7.300e-02,8.067e-02,8.777e-02,9.005e-02,9.357e-02,1.043e-01,1.106e-01,};
vector<float> sf_Vpt_b1_v10_data17_metlt125_ee2j = {1.195e-02,1.350e-02,1.528e-02,1.788e-02,2.012e-02,2.206e-02,2.397e-02,2.579e-02,2.904e-02,3.067e-02,3.456e-02,3.458e-02,3.650e-02,4.175e-02,4.518e-02,4.850e-02,5.350e-02,5.474e-02,6.227e-02,7.028e-02,7.386e-02,8.873e-02,};
vector<float> sf_Vpt_b0_v10_data18_metlt125_eq0j = {1.839e-02,2.061e-02,2.274e-02,2.642e-02,2.920e-02,3.369e-02,4.597e-02,7.519e-02,9.903e-02,1.412e-01,};
vector<float> sf_Vpt_b1_v10_data18_metlt125_eq1j = {2.075e-02,2.323e-02,2.676e-02,2.927e-02,3.279e-02,3.566e-02,3.850e-02,4.307e-02,4.656e-02,5.134e-02,5.247e-02,5.569e-02,5.818e-02,6.278e-02,6.826e-02,7.196e-02,7.929e-02,8.164e-02,8.967e-02,9.408e-02,9.832e-02,9.757e-02,};
vector<float> sf_Vpt_b1_v10_data18_metlt125_ee2j = {1.156e-02,1.329e-02,1.434e-02,1.744e-02,2.034e-02,2.105e-02,2.311e-02,2.502e-02,2.781e-02,3.017e-02,3.253e-02,3.608e-02,3.694e-02,3.924e-02,4.383e-02,4.576e-02,5.225e-02,5.643e-02,6.039e-02,6.558e-02,7.161e-02,7.666e-02,};
vector<float> sferr_Vpt_b0_v10_data16_metlt125_eq0j = {1.813e-04,2.586e-04,3.819e-04,5.274e-04,7.434e-04,6.349e-04,8.028e-04,1.694e-03,2.472e-03,4.931e-03,};
vector<float> sferr_Vpt_b1_v10_data16_metlt125_eq1j = {1.614e-04,2.134e-04,2.642e-04,3.509e-04,4.266e-04,4.676e-04,4.223e-04,5.130e-04,6.119e-04,6.766e-04,7.591e-04,8.994e-04,1.027e-03,8.925e-04,1.075e-03,1.265e-03,1.319e-03,1.742e-03,2.049e-03,2.153e-03,2.476e-03,4.589e-03,};
vector<float> sferr_Vpt_b1_v10_data16_metlt125_ee2j = {1.723e-04,2.173e-04,2.615e-04,3.278e-04,4.121e-04,4.096e-04,4.070e-04,4.537e-04,5.246e-04,5.918e-04,6.579e-04,7.643e-04,8.409e-04,7.169e-04,8.491e-04,1.019e-03,1.017e-03,1.340e-03,1.447e-03,1.601e-03,1.757e-03,2.990e-03,};
vector<float> sferr_Vpt_b0_v10_data17_metlt125_eq0j = {2.445e-04,3.537e-04,5.277e-04,8.707e-04,1.058e-03,1.175e-03,1.569e-03,4.930e-03,9.512e-03,2.110e-02,};
vector<float> sferr_Vpt_b1_v10_data17_metlt125_eq1j = {2.001e-04,2.735e-04,3.549e-04,4.687e-04,5.828e-04,6.164e-04,5.498e-04,6.492e-04,7.553e-04,7.802e-04,9.070e-04,1.060e-03,1.261e-03,1.083e-03,1.167e-03,1.403e-03,1.525e-03,2.031e-03,1.999e-03,2.269e-03,2.531e-03,4.128e-03,};
vector<float> sferr_Vpt_b1_v10_data17_metlt125_ee2j = {2.011e-04,2.536e-04,3.159e-04,4.093e-04,5.075e-04,5.198e-04,4.555e-04,5.312e-04,6.318e-04,6.293e-04,7.541e-04,8.169e-04,9.275e-04,8.217e-04,9.002e-04,1.042e-03,1.074e-03,1.338e-03,1.421e-03,1.641e-03,1.695e-03,2.804e-03,};
vector<float> sferr_Vpt_b0_v10_data18_metlt125_eq0j = {2.402e-04,3.609e-04,5.278e-04,8.093e-04,1.189e-03,1.392e-03,1.402e-03,4.493e-03,7.344e-03,1.661e-02,};
vector<float> sferr_Vpt_b1_v10_data18_metlt125_eq1j = {2.067e-04,2.673e-04,3.614e-04,4.538e-04,5.887e-04,7.269e-04,8.942e-04,6.398e-04,7.517e-04,7.603e-04,8.623e-04,1.008e-03,1.169e-03,1.032e-03,1.037e-03,1.161e-03,1.256e-03,1.601e-03,1.702e-03,1.947e-03,1.919e-03,2.928e-03,};
vector<float> sferr_Vpt_b1_v10_data18_metlt125_ee2j = {2.009e-04,2.598e-04,3.042e-04,4.163e-04,5.442e-04,5.999e-04,7.304e-04,5.111e-04,5.972e-04,5.956e-04,6.859e-04,8.209e-04,8.960e-04,7.509e-04,7.906e-04,8.680e-04,9.184e-04,1.184e-03,1.220e-03,1.383e-03,1.395e-03,2.125e-03,};

vector<float> sf_Vpt_b0_v10_ee_data16_metlt125_eq0j = {6.328e-03,7.035e-03,8.091e-03,8.664e-03,9.456e-03,1.046e-02,1.332e-02,1.483e-02,1.539e-02,1.531e-02,};
vector<float> sf_Vpt_b1_v10_ee_data16_metlt125_eq1j = {7.136e-03,8.279e-03,8.940e-03,1.044e-02,1.113e-02,1.271e-02,1.406e-02,1.551e-02,1.733e-02,1.861e-02,1.970e-02,2.083e-02,2.193e-02,2.377e-02,2.462e-02,2.833e-02,2.910e-02,3.194e-02,3.524e-02,3.550e-02,3.819e-02,4.480e-02,};
vector<float> sf_Vpt_b1_v10_ee_data16_metlt125_ee2j = {3.855e-03,4.476e-03,4.869e-03,5.773e-03,6.622e-03,6.964e-03,8.601e-03,8.520e-03,9.386e-03,1.113e-02,1.149e-02,1.212e-02,1.321e-02,1.423e-02,1.562e-02,1.762e-02,1.994e-02,2.114e-02,2.135e-02,2.514e-02,2.877e-02,3.128e-02,};
vector<float> sf_Vpt_b0_v10_ee_data17_metlt125_eq0j = {6.956e-03,7.562e-03,8.478e-03,1.022e-02,9.873e-03,1.431e-02,1.744e-02,2.574e-02,3.298e-02,4.567e-02,};
vector<float> sf_Vpt_b1_v10_ee_data17_metlt125_eq1j = {7.829e-03,9.160e-03,1.038e-02,1.178e-02,1.295e-02,1.411e-02,1.633e-02,1.706e-02,1.843e-02,1.964e-02,2.093e-02,2.227e-02,2.414e-02,2.566e-02,2.757e-02,3.082e-02,3.184e-02,3.749e-02,3.638e-02,3.974e-02,4.286e-02,4.895e-02,};
vector<float> sf_Vpt_b1_v10_ee_data17_metlt125_ee2j = {4.709e-03,5.260e-03,5.981e-03,7.079e-03,7.986e-03,8.732e-03,9.469e-03,1.005e-02,1.184e-02,1.209e-02,1.373e-02,1.451e-02,1.483e-02,1.679e-02,1.840e-02,2.021e-02,2.279e-02,2.190e-02,2.527e-02,2.965e-02,3.185e-02,3.704e-02,};
vector<float> sf_Vpt_b0_v10_ee_data18_metlt125_eq0j = {7.298e-03,8.222e-03,8.934e-03,1.076e-02,1.170e-02,1.328e-02,1.902e-02,2.881e-02,4.368e-02,6.275e-02,};
vector<float> sf_Vpt_b1_v10_ee_data18_metlt125_eq1j = {8.257e-03,9.370e-03,1.068e-02,1.180e-02,1.324e-02,1.446e-02,1.561e-02,1.737e-02,1.896e-02,2.138e-02,2.205e-02,2.319e-02,2.413e-02,2.658e-02,2.853e-02,3.035e-02,3.359e-02,3.505e-02,3.810e-02,3.956e-02,4.220e-02,4.546e-02,};
vector<float> sf_Vpt_b1_v10_ee_data18_metlt125_ee2j = {4.602e-03,5.365e-03,5.730e-03,7.026e-03,8.317e-03,8.846e-03,9.313e-03,1.038e-02,1.131e-02,1.234e-02,1.381e-02,1.529e-02,1.508e-02,1.606e-02,1.802e-02,1.926e-02,2.199e-02,2.440e-02,2.558e-02,2.793e-02,3.080e-02,3.502e-02,};
vector<float> sferr_Vpt_b0_v10_ee_data16_metlt125_eq0j = {8.751e-05,1.248e-04,1.835e-04,2.519e-04,3.457e-04,3.234e-04,4.503e-04,1.002e-03,1.594e-03,3.392e-03,};
vector<float> sferr_Vpt_b1_v10_ee_data16_metlt125_eq1j = {7.780e-05,1.017e-04,1.240e-04,1.629e-04,1.946e-04,2.242e-04,2.284e-04,2.768e-04,3.367e-04,3.804e-04,4.355e-04,5.074e-04,5.826e-04,5.020e-04,6.090e-04,7.610e-04,7.784e-04,1.042e-03,1.224e-03,1.333e-03,1.524e-03,2.747e-03,};
vector<float> sferr_Vpt_b1_v10_ee_data16_metlt125_ee2j = {9.057e-05,1.119e-04,1.314e-04,1.639e-04,2.015e-04,2.111e-04,2.367e-04,2.581e-04,2.978e-04,3.508e-04,3.886e-04,4.371e-04,4.951e-04,4.156e-04,5.054e-04,6.162e-04,6.330e-04,8.120e-04,8.748e-04,1.004e-03,1.129e-03,1.831e-03,};
vector<float> sferr_Vpt_b0_v10_ee_data17_metlt125_eq0j = {1.128e-04,1.616e-04,2.374e-04,3.839e-04,4.625e-04,5.580e-04,7.646e-04,2.391e-03,4.528e-03,1.108e-02,};
vector<float> sferr_Vpt_b1_v10_ee_data17_metlt125_eq1j = {9.168e-05,1.234e-04,1.598e-04,2.079e-04,2.587e-04,2.792e-04,2.814e-04,3.283e-04,3.857e-04,4.153e-04,4.864e-04,5.652e-04,6.699e-04,5.726e-04,6.596e-04,8.231e-04,8.529e-04,1.182e-03,1.205e-03,1.417e-03,1.576e-03,2.668e-03,};
vector<float> sferr_Vpt_b1_v10_ee_data17_metlt125_ee2j = {9.756e-05,1.204e-04,1.481e-04,1.910e-04,2.344e-04,2.493e-04,2.416e-04,2.768e-04,3.387e-04,3.480e-04,4.145e-04,4.641e-04,5.138e-04,4.464e-04,5.205e-04,6.245e-04,6.478e-04,7.770e-04,8.718e-04,1.031e-03,1.091e-03,1.768e-03,};
vector<float> sferr_Vpt_b0_v10_ee_data18_metlt125_eq0j = {1.073e-04,1.604e-04,2.297e-04,3.594e-04,5.169e-04,5.912e-04,6.892e-04,2.090e-03,4.069e-03,1.000e-02,};
vector<float> sferr_Vpt_b1_v10_ee_data18_metlt125_eq1j = {9.162e-05,1.188e-04,1.575e-04,1.983e-04,2.559e-04,3.155e-04,3.865e-04,3.049e-04,3.616e-04,3.921e-04,4.465e-04,5.154e-04,5.928e-04,5.271e-04,5.683e-04,6.690e-04,7.198e-04,9.224e-04,1.022e-03,1.168e-03,1.225e-03,1.950e-03,};
vector<float> sferr_Vpt_b1_v10_ee_data18_metlt125_ee2j = {9.406e-05,1.213e-04,1.399e-04,1.892e-04,2.472e-04,2.781e-04,3.245e-04,2.636e-04,3.034e-04,3.213e-04,3.760e-04,4.447e-04,4.724e-04,3.940e-04,4.481e-04,5.172e-04,5.421e-04,7.057e-04,7.485e-04,8.526e-04,8.971e-04,1.408e-03,};

vector<float> sf_Vpt_b0_v10_mumu_data16_metlt125_eq0j = {1.046e-02,1.133e-02,1.271e-02,1.355e-02,1.526e-02,1.634e-02,2.008e-02,2.109e-02,1.878e-02,1.603e-02,};
vector<float> sf_Vpt_b1_v10_mumu_data16_metlt125_eq1j = {1.160e-02,1.327e-02,1.448e-02,1.677e-02,1.818e-02,2.028e-02,2.222e-02,2.416e-02,2.565e-02,2.942e-02,2.912e-02,3.179e-02,3.269e-02,3.553e-02,3.864e-02,4.178e-02,4.500e-02,4.694e-02,5.121e-02,5.262e-02,5.735e-02,7.220e-02,};
vector<float> sf_Vpt_b1_v10_mumu_data16_metlt125_ee2j = {6.079e-03,7.119e-03,7.993e-03,9.185e-03,1.069e-02,1.103e-02,1.269e-02,1.329e-02,1.474e-02,1.677e-02,1.738e-02,1.983e-02,1.969e-02,2.189e-02,2.340e-02,2.714e-02,2.779e-02,3.174e-02,3.224e-02,3.660e-02,3.842e-02,4.843e-02,};
vector<float> sf_Vpt_b0_v10_mumu_data17_metlt125_eq0j = {1.080e-02,1.180e-02,1.328e-02,1.605e-02,1.591e-02,2.128e-02,2.826e-02,4.359e-02,6.741e-02,9.135e-02,};
vector<float> sf_Vpt_b1_v10_mumu_data17_metlt125_eq1j = {1.213e-02,1.418e-02,1.580e-02,1.805e-02,1.955e-02,2.173e-02,2.372e-02,2.513e-02,2.667e-02,2.947e-02,3.046e-02,3.234e-02,3.467e-02,3.673e-02,3.902e-02,4.218e-02,4.882e-02,5.028e-02,5.367e-02,5.383e-02,6.149e-02,6.170e-02,};
vector<float> sf_Vpt_b1_v10_mumu_data17_metlt125_ee2j = {7.246e-03,8.244e-03,9.303e-03,1.080e-02,1.213e-02,1.333e-02,1.450e-02,1.574e-02,1.720e-02,1.858e-02,2.082e-02,2.007e-02,2.167e-02,2.497e-02,2.678e-02,2.829e-02,3.071e-02,3.284e-02,3.700e-02,4.063e-02,4.200e-02,5.169e-02,};
vector<float> sf_Vpt_b0_v10_mumu_data18_metlt125_eq0j = {1.110e-02,1.238e-02,1.381e-02,1.566e-02,1.750e-02,2.041e-02,2.695e-02,4.638e-02,5.535e-02,7.843e-02,};
vector<float> sf_Vpt_b1_v10_mumu_data18_metlt125_eq1j = {1.250e-02,1.386e-02,1.608e-02,1.747e-02,1.955e-02,2.120e-02,2.288e-02,2.571e-02,2.760e-02,2.996e-02,3.042e-02,3.249e-02,3.406e-02,3.620e-02,3.973e-02,4.160e-02,4.570e-02,4.659e-02,5.157e-02,5.453e-02,5.612e-02,5.211e-02,};
vector<float> sf_Vpt_b1_v10_mumu_data18_metlt125_ee2j = {6.958e-03,7.929e-03,8.607e-03,1.042e-02,1.203e-02,1.221e-02,1.380e-02,1.464e-02,1.650e-02,1.783e-02,1.872e-02,2.079e-02,2.186e-02,2.317e-02,2.582e-02,2.650e-02,3.026e-02,3.203e-02,3.481e-02,3.764e-02,4.081e-02,4.164e-02,};
vector<float> sferr_Vpt_b0_v10_mumu_data16_metlt125_eq0j = {1.255e-04,1.767e-04,2.571e-04,3.531e-04,4.992e-04,4.395e-04,5.773e-04,1.226e-03,1.774e-03,3.474e-03,};
vector<float> sferr_Vpt_b1_v10_mumu_data16_metlt125_eq1j = {1.107e-04,1.445e-04,1.786e-04,2.349e-04,2.857e-04,3.169e-04,3.037e-04,3.661e-04,4.312e-04,4.978e-04,5.485e-04,6.527e-04,7.398e-04,6.403e-04,7.920e-04,9.415e-04,9.895e-04,1.290e-03,1.508e-03,1.637e-03,1.885e-03,3.532e-03,};
vector<float> sferr_Vpt_b1_v10_mumu_data16_metlt125_ee2j = {1.218e-04,1.526e-04,1.840e-04,2.266e-04,2.838e-04,2.872e-04,2.964e-04,3.341e-04,3.870e-04,4.403e-04,4.893e-04,5.763e-04,6.201e-04,5.307e-04,6.318e-04,7.751e-04,7.555e-04,1.009e-03,1.091e-03,1.218e-03,1.311e-03,2.297e-03,};
vector<float> sferr_Vpt_b0_v10_mumu_data17_metlt125_eq0j = {1.602e-04,2.313e-04,3.432e-04,5.619e-04,6.897e-04,7.627e-04,1.080e-03,3.460e-03,7.179e-03,1.646e-02,};
vector<float> sferr_Vpt_b1_v10_mumu_data17_metlt125_eq1j = {1.305e-04,1.769e-04,2.271e-04,2.990e-04,3.687e-04,3.985e-04,3.674e-04,4.340e-04,5.029e-04,5.425e-04,6.238e-04,7.257e-04,8.566e-04,7.321e-04,8.181e-04,9.920e-04,1.103e-03,1.414e-03,1.490e-03,1.668e-03,1.905e-03,3.014e-03,};
vector<float> sferr_Vpt_b1_v10_mumu_data17_metlt125_ee2j = {1.345e-04,1.694e-04,2.089e-04,2.671e-04,3.287e-04,3.442e-04,3.192e-04,3.726e-04,4.340e-04,4.527e-04,5.371e-04,5.681e-04,6.521e-04,5.758e-04,6.490e-04,7.555e-04,7.684e-04,9.804e-04,1.068e-03,1.218e-03,1.259e-03,2.104e-03,};
vector<float> sferr_Vpt_b0_v10_mumu_data18_metlt125_eq0j = {1.531e-04,2.281e-04,3.354e-04,5.006e-04,7.402e-04,8.714e-04,9.022e-04,3.014e-03,4.786e-03,1.143e-02,};
vector<float> sferr_Vpt_b1_v10_mumu_data18_metlt125_eq1j = {1.309e-04,1.671e-04,2.262e-04,2.815e-04,3.634e-04,4.465e-04,5.481e-04,4.150e-04,4.853e-04,5.004e-04,5.638e-04,6.598e-04,7.646e-04,6.641e-04,7.068e-04,8.113e-04,8.712e-04,1.101e-03,1.217e-03,1.402e-03,1.422e-03,2.095e-03,};
vector<float> sferr_Vpt_b1_v10_mumu_data18_metlt125_ee2j = {1.306e-04,1.665e-04,1.953e-04,2.636e-04,3.392e-04,3.671e-04,4.570e-04,3.376e-04,3.982e-04,4.095e-04,4.608e-04,5.482e-04,6.087e-04,5.073e-04,5.586e-04,6.219e-04,6.538e-04,8.292e-04,8.875e-04,1.005e-03,1.038e-03,1.541e-03,};

const vector<float> sf_pt_phvsel_llg = {0.774,1.09,1.222,1.437,1.201,1.002,0.853,0.717,0.788,1.225,0.834,0.724,0.873,0.853,0.842,0.468,0.491,0.476,0.499,};

std::pair<float,float> getBosonPtScale(float Vpt, float Veta, int year, string jetcat, vector<float> ptbin0, vector<float> ptbin1, bool is_data, bool extendEEphoton2j, bool doBosonEtaReweight) {
  float scale = 1.0;
  float sferr = 0.0;
  int icat0 = std::upper_bound(ptbin0.begin(), ptbin0.end(), Vpt) - ptbin0.begin() - 1;
  int icat1 = std::upper_bound(ptbin1.begin(), ptbin1.end(), Vpt) - ptbin1.begin() - 1;
  if (icat0 >= ptbin0.size()-1) icat0 = ptbin0.size()-2;
  if (icat1 >= ptbin1.size()-1) icat1 = ptbin1.size()-2;

  if (year == 2018) {
    if (is_data) {
      if (extendEEphoton2j && jetcat == "_ge2j") {
        scale = sf_Vpt_b1_v10_data18_metlt125_ee2j.at(icat1);
        sferr = sferr_Vpt_b1_v10_data18_metlt125_ee2j.at(icat1);
        if (doBosonEtaReweight) {
          int ibin = (fabs(Veta) >= 2.5)? 24 : int(fabs(Veta) / 0.1);
          scale *= sf_Vaeta_flatpt_v10_data18_metlt125_ee2j.at(ibin);
        }
      }
      else if (jetcat == "_eq0j") {
        scale = sf_Vpt_b0_v10_data18_metlt125_eq0j.at(icat0);
        sferr = sferr_Vpt_b0_v10_data18_metlt125_eq0j.at(icat0);
      }
      else if (jetcat == "_eq1j") {
        scale = sf_Vpt_b1_v10_data18_metlt125_eq1j.at(icat1);
        sferr = sferr_Vpt_b1_v10_data18_metlt125_eq1j.at(icat1);
      }
      else if (jetcat == "_eq2j") scale = sf_Vpt_b1_v4_data18_metlt125_eq2j.at(icat1);
    } else {
      if (extendEEphoton2j && jetcat == "_eq2j") {
        scale = sf_Vpt_b1_v4_closure18_metlt125_ee2j.at(icat1);
        if (doBosonEtaReweight) {
          int ibin = (fabs(Veta) >= 2.5)? 0 : int((Veta + 2.5) / 0.1);
          scale *= sf_Veta_flatpt_v4_closure18_metlt125_ee2j.at(ibin);
        }
      }
      else if (jetcat == "_eq0j") scale = sf_Vpt_b0_v4_closure18_metlt125_eq0j.at(icat0);
      else if (jetcat == "_eq1j") scale = sf_Vpt_b1_v4_closure18_metlt125_eq1j.at(icat1);
      else if (jetcat == "_eq2j") scale = sf_Vpt_b1_v4_closure18_metlt125_eq2j.at(icat1);
    }
  } else if (year == 2017) {
    if (is_data) {
      if (extendEEphoton2j && jetcat == "_ge2j") {
        scale = sf_Vpt_b1_v10_data17_metlt125_ee2j.at(icat1);
        sferr = sferr_Vpt_b1_v10_data17_metlt125_ee2j.at(icat1);
        if (doBosonEtaReweight) {
          int ibin = (fabs(Veta) >= 2.5)? 24 : int(fabs(Veta) / 0.1);
          scale *= sf_Vaeta_flatpt_v10_data17_metlt125_ee2j.at(ibin);
        }
      }
      else if (jetcat == "_eq0j") {
        scale = sf_Vpt_b0_v10_data17_metlt125_eq0j.at(icat0);
        sferr = sferr_Vpt_b0_v10_data17_metlt125_eq0j.at(icat0);
      }
      else if (jetcat == "_eq1j") {
        scale = sf_Vpt_b1_v10_data17_metlt125_eq1j.at(icat1);
        sferr = sferr_Vpt_b1_v10_data17_metlt125_eq1j.at(icat1);
      }
      else if (jetcat == "_eq2j") scale = sf_Vpt_b1_v4_data17_metlt125_eq2j.at(icat1);
    } else {
      if (extendEEphoton2j && jetcat == "_eq2j") {
        scale = sf_Vpt_b1_v4_closure17_metlt125_ee2j.at(icat1);
        if (doBosonEtaReweight) {
          int ibin = (fabs(Veta) >= 2.5)? 0 : int((Veta + 2.5) / 0.1);
          scale *= sf_Veta_flatpt_v4_closure17_metlt125_ee2j.at(ibin);
        }
      }
      else if (jetcat == "_eq0j") scale = sf_Vpt_b0_v4_closure17_metlt125_eq0j.at(icat0);
      else if (jetcat == "_eq1j") scale = sf_Vpt_b1_v4_closure17_metlt125_eq1j.at(icat1);
      else if (jetcat == "_eq2j") scale = sf_Vpt_b1_v4_closure17_metlt125_eq2j.at(icat1);
    }
  } else if (year == 2016) {
    if (is_data) {
      if (extendEEphoton2j && jetcat == "_ge2j") {
        scale = sf_Vpt_b1_v10_data16_metlt125_ee2j.at(icat1);
        sferr = sferr_Vpt_b1_v10_data16_metlt125_ee2j.at(icat1);
        if (doBosonEtaReweight) {
          int ibin = (fabs(Veta) >= 2.5)? 24 : int(fabs(Veta) / 0.1);
          scale *= sf_Vaeta_flatpt_v10_data16_metlt125_ee2j.at(ibin);
        }
      }
      else if (jetcat == "_eq0j") {
        scale = sf_Vpt_b0_v10_data16_metlt125_eq0j.at(icat0);
        sferr = sferr_Vpt_b0_v10_data16_metlt125_eq0j.at(icat0);
      }
      else if (jetcat == "_eq1j") {
        scale = sf_Vpt_b1_v10_data16_metlt125_eq1j.at(icat1);
        sferr = sferr_Vpt_b1_v10_data16_metlt125_eq1j.at(icat1);
      }
      else if (jetcat == "_eq2j") scale = sf_Vpt_b1_v4_data16_metlt125_eq2j.at(icat1);
    } else {
      if (extendEEphoton2j && jetcat == "_eq2j") {
        scale = sf_Vpt_b1_v4_closure16_metlt125_ee2j.at(icat1);
        if (doBosonEtaReweight) {
          int ibin = (fabs(Veta) >= 2.5)? 0 : int((Veta + 2.5) / 0.1);
          scale *= sf_Veta_flatpt_v4_closure16_metlt125_ee2j.at(ibin);
        }
      }
      else if (jetcat == "_eq0j") scale = sf_Vpt_b0_v4_closure16_metlt125_eq0j.at(icat0);
      else if (jetcat == "_eq1j") scale = sf_Vpt_b1_v4_closure16_metlt125_eq1j.at(icat1);
      else if (jetcat == "_eq2j") scale = sf_Vpt_b1_v4_closure16_metlt125_eq2j.at(icat1);
    }
  }
  return std::make_pair(scale,sferr);
}

std::tuple<float,float,float,float> getBosonPtScale2(float Vpt, float Veta, int year, int njet, vector<float> ptbin0, vector<float> ptbin1, bool doBosonEtaReweight) {
  float sfval_ee(1.0), sfval_mumu(1.0);
  float sferr_ee(0.0), sferr_mumu(0.0);
  int icat0 = std::upper_bound(ptbin0.begin(), ptbin0.end(), Vpt) - ptbin0.begin() - 1;
  int icat1 = std::upper_bound(ptbin1.begin(), ptbin1.end(), Vpt) - ptbin1.begin() - 1;
  if (icat0 >= ptbin0.size()-1) icat0 = ptbin0.size()-2;
  if (icat1 >= ptbin1.size()-1) icat1 = ptbin1.size()-2;

  // if (!is_data || !extendEEphoton2j) {
  //   throw std::invalid_argument("[getBosonPtScale2] >> Have to run uneer data and ee2j option.");
  // }

  if (year == 2018) {
    if (njet >= 2) {
      sfval_ee = sf_Vpt_b1_v10_ee_data18_metlt125_ee2j.at(icat1);
      sferr_ee = sferr_Vpt_b1_v10_ee_data18_metlt125_ee2j.at(icat1);
      sfval_mumu = sf_Vpt_b1_v10_mumu_data18_metlt125_ee2j.at(icat1);
      sferr_mumu = sferr_Vpt_b1_v10_mumu_data18_metlt125_ee2j.at(icat1);
      if (doBosonEtaReweight) {
        int ibin = (fabs(Veta) >= 2.5)? 24 : int(fabs(Veta) / 0.1);
        sfval_ee *= sf_Vaeta_flatpt_v10_ee_data18_metlt125_ee2j.at(ibin);
        sfval_mumu *= sf_Vaeta_flatpt_v10_mumu_data18_metlt125_ee2j.at(ibin);
      }
    }
    else if (njet == 0) {
      sfval_ee = sf_Vpt_b0_v10_ee_data18_metlt125_eq0j.at(icat0);
      sferr_ee = sferr_Vpt_b0_v10_ee_data18_metlt125_eq0j.at(icat0);
      sfval_mumu = sf_Vpt_b0_v10_mumu_data18_metlt125_eq0j.at(icat0);
      sferr_mumu = sferr_Vpt_b0_v10_mumu_data18_metlt125_eq0j.at(icat0);
    }
    else if (njet == 1) {
      sfval_ee = sf_Vpt_b1_v10_ee_data18_metlt125_eq1j.at(icat1);
      sferr_ee = sferr_Vpt_b1_v10_ee_data18_metlt125_eq1j.at(icat1);
      sfval_mumu = sf_Vpt_b1_v10_mumu_data18_metlt125_eq1j.at(icat1);
      sferr_mumu = sferr_Vpt_b1_v10_mumu_data18_metlt125_eq1j.at(icat1);
    }
  } else if (year == 2017) {
    if (njet >= 2) {
      sfval_ee = sf_Vpt_b1_v10_ee_data17_metlt125_ee2j.at(icat1);
      sferr_ee = sferr_Vpt_b1_v10_ee_data17_metlt125_ee2j.at(icat1);
      sfval_mumu = sf_Vpt_b1_v10_mumu_data17_metlt125_ee2j.at(icat1);
      sferr_mumu = sferr_Vpt_b1_v10_mumu_data17_metlt125_ee2j.at(icat1);
      if (doBosonEtaReweight) {
        int ibin = (fabs(Veta) >= 2.5)? 24 : int(fabs(Veta) / 0.1);
        sfval_ee *= sf_Vaeta_flatpt_v10_ee_data17_metlt125_ee2j.at(ibin);
        sfval_mumu *= sf_Vaeta_flatpt_v10_mumu_data17_metlt125_ee2j.at(ibin);
      }
    }
    else if (njet == 0) {
      sfval_ee = sf_Vpt_b0_v10_ee_data17_metlt125_eq0j.at(icat0);
      sferr_ee = sferr_Vpt_b0_v10_ee_data17_metlt125_eq0j.at(icat0);
      sfval_mumu = sf_Vpt_b0_v10_mumu_data17_metlt125_eq0j.at(icat0);
      sferr_mumu = sferr_Vpt_b0_v10_mumu_data17_metlt125_eq0j.at(icat0);
    }
    else if (njet == 1) {
      sfval_ee = sf_Vpt_b1_v10_ee_data17_metlt125_eq1j.at(icat1);
      sferr_ee = sferr_Vpt_b1_v10_ee_data17_metlt125_eq1j.at(icat1);
      sfval_mumu = sf_Vpt_b1_v10_mumu_data17_metlt125_eq1j.at(icat1);
      sferr_mumu = sferr_Vpt_b1_v10_mumu_data17_metlt125_eq1j.at(icat1);
    }
  } else if (year == 2016) {
    if (njet >= 2) {
      sfval_ee = sf_Vpt_b1_v10_ee_data16_metlt125_ee2j.at(icat1);
      sferr_ee = sferr_Vpt_b1_v10_ee_data16_metlt125_ee2j.at(icat1);
      sfval_mumu = sf_Vpt_b1_v10_mumu_data16_metlt125_ee2j.at(icat1);
      sferr_mumu = sferr_Vpt_b1_v10_mumu_data16_metlt125_ee2j.at(icat1);
      if (doBosonEtaReweight) {
        int ibin = (fabs(Veta) >= 2.5)? 24 : int(fabs(Veta) / 0.1);
        sfval_ee *= sf_Vaeta_flatpt_v10_ee_data16_metlt125_ee2j.at(ibin);
        sfval_mumu *= sf_Vaeta_flatpt_v10_mumu_data16_metlt125_ee2j.at(ibin);
      }
    }
    else if (njet == 0) {
      sfval_ee = sf_Vpt_b0_v10_ee_data16_metlt125_eq0j.at(icat0);
      sferr_ee = sferr_Vpt_b0_v10_ee_data16_metlt125_eq0j.at(icat0);
      sfval_mumu = sf_Vpt_b0_v10_mumu_data16_metlt125_eq0j.at(icat0);
      sferr_mumu = sferr_Vpt_b0_v10_mumu_data16_metlt125_eq0j.at(icat0);
    }
    else if (njet == 1) {
      sfval_ee = sf_Vpt_b1_v10_ee_data16_metlt125_eq1j.at(icat1);
      sferr_ee = sferr_Vpt_b1_v10_ee_data16_metlt125_eq1j.at(icat1);
      sfval_mumu = sf_Vpt_b1_v10_mumu_data16_metlt125_eq1j.at(icat1);
      sferr_mumu = sferr_Vpt_b1_v10_mumu_data16_metlt125_eq1j.at(icat1);
    }
  }
  return std::make_tuple(sfval_ee,sfval_mumu,sferr_ee,sferr_mumu);
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
  if (year == 2016) sfval *= 35.9/36.33;

  return std::make_pair(sfval, sferr);
}

// void getPoissonCountingConfidenceInterval_Frequentist(double sw_total, double swsq_total, double CL, double& vlow, double& vhigh){
//   double const quant = (1. - CL) / 2.;
//   double const count = (swsq_total<=0. ? (sw_total==0. ? 0. : sw_total) : std::pow(sw_total, 2)/swsq_total);
//   vlow = (count == 0. ? 0. : ROOT::Math::chisquared_quantile(quant, 2. * count) / 2.);
//   vhigh = ROOT::Math::chisquared_quantile_c(quant, 2 * (count + 1.)) / 2.;
//   if (count>0.){
//     vlow *= sw_total/count;
//     vhigh *= sw_total/count;
//   }
// }

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
