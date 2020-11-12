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

inline void conditionalizeHist(TH2* h2d, bool invert = false) {
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

void linkHist(string hnew, string hexist, std::map<std::string, TH1*> &allhistos)
{
  // Could be useful when mutiple ratio hists sharing a common denominator
  if (allhistos.count(hnew)) return;
  auto iter = allhistos.find(hexist);
  if (iter == allhistos.end()) throw std::logic_error("linkHist(): Histogram "+hexist+" need to be plotted first");
  allhistos.insert( std::pair<std::string, TH1*>(hnew, iter->second) );
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

  // vector<float> ptbin0;
  // for (int i = 55; i < 160; i += 5) ptbin0.emplace_back(i);
  // for (int i = 160; i < 200; i += 10) ptbin0.emplace_back(i);
  // ptbin0.emplace_back(250); ptbin0.emplace_back(450);

  // vector<float> ptbin1;
  // for (int i = 55; i < 160; i += 5) ptbin1.emplace_back(i);
  // for (int i = 160; i < 300; i += 10) ptbin1.emplace_back(i);
  // for (int i = 300; i <= 350; i += 25) ptbin1.emplace_back(i);
  // ptbin1.emplace_back(400); ptbin1.emplace_back(600);

  vector<float> ptbin0 = {55, 60, 65, 70, 75, 80, 90, 105, 130, 160, 200, 250, 450};

  vector<float> ptbin1;
  for (int i = 55;  i < 120; i += 5)  ptbin1.emplace_back(i);
  for (int i = 120; i < 150; i += 10) ptbin1.emplace_back(i);
  for (int i = 150; i < 180; i += 15) ptbin1.emplace_back(i);
  for (int i = 180; i < 220; i += 20) ptbin1.emplace_back(i);
  for (int i = 220; i < 270; i += 25) ptbin1.emplace_back(i);
  ptbin1.emplace_back(270); ptbin1.emplace_back(300); ptbin1.emplace_back(350);
  ptbin1.emplace_back(425); ptbin1.emplace_back(600);

  return std::make_tuple(ptbin0, ptbin1);
}


// Transfer factors from elec to photon
// binned in pt as vector<float> tnpbin1 = { 55,    82.5,    100,    135,    175,    190,    220,   270,   450};
const vector<float> sfval_elpt_tophCR_v2_16 = {0.0320, 0.0260, 0.0184, 0.0136, 0.0154, 0.0095, 0.013,  0.0142,};
const vector<float> sfval_elpt_tophCR_v2_17 = {0.0352, 0.0289, 0.0242, 0.0196, 0.0214, 0.0183, 0.0155, 0.0177,};
const vector<float> sfval_elpt_tophCR_v2_18 = {0.0334, 0.0276, 0.024,  0.0202, 0.0209, 0.021,  0.0168, 0.0153,};

const vector<float> sferr_elpt_tophCR_v2_16 = {0.0004, 0.0009, 0.0008, 0.0011, 0.0025, 0.0016, 0.0021, 0.0025};
const vector<float> sferr_elpt_tophCR_v2_17 = {0.0004, 0.0008, 0.0008, 0.0011, 0.0027, 0.0022, 0.0020, 0.0024};
const vector<float> sferr_elpt_tophCR_v2_18 = {0.0003, 0.0007, 0.0007, 0.0010, 0.0022, 0.0019, 0.0018, 0.0019};

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
 
 
std::pair<float,float> getElecToGammaRate(float Vpt, int year) {
  // Apply the single electron to photon fake rate measured from tag and probe
  const vector<float> tnpbin1 = {55, 82.5, 100, 135, 175, 190, 220, 270, 450};
  const vector<float> tnpbin2 = {55, 82.5, 100, 135, 180, 190, 450};
  int icat = std::min(std::upper_bound(tnpbin1.begin(), tnpbin1.end(), Vpt) - tnpbin1.begin() - 1, 7L);
  if (icat < 0 || Vpt < 54) {
    cout << "[getElecToGammaRate] >> WARNING: Input Vpt= " << Vpt << ", icat= " << icat << endl;
    return std::make_pair(1, 0);
  }
  float scale = 1, sferr = 0;
  if (year == 2016) {
    scale = sfval_elpt_tophCR_v2_16.at(icat);
    sferr = sferr_elpt_tophCR_v2_16.at(icat);
  }
  if (year == 2017) {
    scale = sfval_elpt_tophCR_v2_17.at(icat);
    sferr = sferr_elpt_tophCR_v2_17.at(icat);
  }
  if (year == 2018) {
    scale = sfval_elpt_tophCR_v2_18.at(icat);
    sferr = sferr_elpt_tophCR_v2_18.at(icat);
  }

  // Apply the single-photon/single-electron trigger weight
  const vector<float> trigptbin0 = {55, 82.5, 100, 135, 180, 220, 450};
  const vector<float> trigptbin1 = {55, 100, 220, 450};
  int trigcat0 = std::min(std::upper_bound(trigptbin0.begin(), trigptbin0.end(), Vpt) - trigptbin0.begin() - 1, 5L);
  int trigcat1 = std::min(std::upper_bound(trigptbin1.begin(), trigptbin1.end(), Vpt) - trigptbin1.begin() - 1, 2L);
  if (year == 2016) {
    scale *= trigeff_ph_data_barrel_v1_16.at(trigcat0) / trigeff_el_data_barrel_v1_16.at(trigcat1);
    sferr = sqrt(sferr*sferr + pow(tefferr_ph_data_barrel_v1_16.at(trigcat0), 2));
  }
  if (year == 2017) {
    scale *= trigeff_ph_data_barrel_v1_17.at(trigcat0) / trigeff_el_data_barrel_v1_17.at(trigcat1);
    sferr = sqrt(sferr*sferr + pow(tefferr_ph_data_barrel_v1_17.at(trigcat0), 2));
  }
  if (year == 2018) {
    scale *= trigeff_ph_data_barrel_v1_18.at(trigcat0) / trigeff_el_data_barrel_v1_18.at(trigcat1);
    sferr = sqrt(sferr*sferr + pow(tefferr_ph_data_barrel_v1_18.at(trigcat0), 2));
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

// vector<float> boson_aeta_metlt125_eq2j_18 = {1.567e-02,1.526e-02,1.662e-02,1.749e-02,1.858e-02,1.633e-02,1.713e-02,2.057e-02,1.916e-02,2.136e-02,2.209e-02,2.648e-02,2.876e-02,2.723e-02,3.645e-02,3.295e-02,2.063e-02,2.127e-02,1.813e-02,1.594e-02,1.591e-02,1.733e-02,1.970e-02,1.436e-02,1.004e-02,1.000e+00,};

const vector<float> sf_pt_phvsel_llg = {0.774,1.09,1.222,1.437,1.201,1.002,0.853,0.717,0.788,1.225,0.834,0.724,0.873,0.853,0.842,0.468,0.491,0.476,0.499,};

float getBosonPtScale(float Vpt, float Veta, int year, string jetcat, vector<float> ptbin0, vector<float> ptbin1, bool is_data, bool extendEEphoton2j, bool doBosonEtaReweight) {
  float scale = 1.0;
  int icat0 = std::upper_bound(ptbin0.begin(), ptbin0.end(), Vpt) - ptbin0.begin() - 1;
  int icat1 = std::upper_bound(ptbin1.begin(), ptbin1.end(), Vpt) - ptbin1.begin() - 1;
  if (icat0 >= ptbin0.size()-1) icat0 = ptbin0.size()-2;
  if (icat1 >= ptbin1.size()-1) icat1 = ptbin1.size()-2;

  if (year == 2018) {
    if (is_data) {
      if (extendEEphoton2j && jetcat == "_eq2j") {
        scale = sf_Vpt_b1_v6_data18_metlt125_ee2j.at(icat1);
        if (doBosonEtaReweight) {
          // int ibin = (fabs(Veta) >= 2.5)? 0 : int((Veta + 2.5) / 0.1);
          int ibin = (fabs(Veta) >= 2.5)? 24 : int(fabs(Veta) / 0.1);
          scale *= sf_Vaeta_flatpt_v6_data18_metlt125_ee2j.at(ibin);
        }
      }
      else if (jetcat == "_eq0j") scale = sf_Vpt_b0_v6_data18_metlt125_eq0j.at(icat0);
      else if (jetcat == "_eq1j") scale = sf_Vpt_b1_v6_data18_metlt125_eq1j.at(icat1);
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
      if (extendEEphoton2j && jetcat == "_eq2j") {
        scale = sf_Vpt_b1_v6_data17_metlt125_ee2j.at(icat1);
        if (doBosonEtaReweight) {
          // int ibin = (fabs(Veta) >= 2.5)? 0 : int((Veta + 2.5) / 0.1);
          int ibin = (fabs(Veta) >= 2.5)? 24 : int(fabs(Veta) / 0.1);
          scale *= sf_Vaeta_flatpt_v6_data17_metlt125_ee2j.at(ibin);
        }
      }
      else if (jetcat == "_eq0j") scale = sf_Vpt_b0_v6_data17_metlt125_eq0j.at(icat0);
      else if (jetcat == "_eq1j") scale = sf_Vpt_b1_v6_data17_metlt125_eq1j.at(icat1);
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
      if (extendEEphoton2j && jetcat == "_eq2j") {
        scale = sf_Vpt_b1_v6_data16_metlt125_ee2j.at(icat1);
        if (doBosonEtaReweight) {
          // int ibin = (fabs(Veta) >= 2.5)? 0 : int((Veta + 2.5) / 0.1);
          int ibin = (fabs(Veta) >= 2.5)? 24 : int(fabs(Veta) / 0.1);
          scale *= sf_Vaeta_flatpt_v6_data16_metlt125_ee2j.at(ibin);
        }
      }
      else if (jetcat == "_eq0j") scale = sf_Vpt_b0_v6_data16_metlt125_eq0j.at(icat0);
      else if (jetcat == "_eq1j") scale = sf_Vpt_b1_v6_data16_metlt125_eq1j.at(icat1);
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
  return scale;
}

float getNjetSFfromLLG(int njet, int year) {
  if (year == 2018) {
    if (njet == 0) return 1.0;
    if (njet == 1) return 1.0;
    if (njet == 2) return 1.0;
  }
  if (year == 2017) {
    if (njet == 0) return 1.0;
    if (njet == 1) return 1.0;
    if (njet == 2) return 1.0;
  }
  if (year == 2016) {
    if (njet == 0) return 1.07;
    if (njet == 1) return 1.00;
    if (njet == 2) return 0.82;
  }

  return 1.0;
}

}

#endif  // HZZLOOPER_UTILITIES_H
