#include "IsolationTools.h"

using namespace tas;

bool passMultiIso(int id, int idx, float cutMiniIso, float cutPtRatio, float cutPtRel) {
    float miniiso = abs(id)==11 ?  Electron_miniPFRelIso_all()[idx] : Muon_miniPFRelIso_all()[idx];
    if (miniiso > cutMiniIso) return false;
    // short circuit if we only are requiring miniiso (e.g., veto)
    if (cutPtRatio < 0 && cutPtRel < 0) return true;
    int closestJetIdx = abs(id) == 11 ? Electron_jetIdx()[idx] : Muon_jetIdx()[idx];
    float ptratio = 1.;
    float ptrel = 0.;
    if (closestJetIdx < 0) {
        ptratio = 1.;
        ptrel = 0.;
    } else {
        // NOTE missing JEC L2L3 (or L1) to correctly compute the closest jet p4
        LorentzVector closestJet = Jet_p4()[closestJetIdx];
        LorentzVector lepton = abs(id) == 11 ? Electron_p4()[idx] : Muon_p4()[idx];
        ptratio = getPtRatio(id,idx);
        ptrel = computePtRel(lepton,closestJet);
    }
    return (ptratio > cutPtRatio) || (ptrel > cutPtRel);
}

float computePtRel(const LorentzVector& lepp4, const LorentzVector& jetp4) {
    if (jetp4.pt()<=0.) return 0.;
    LorentzVector jp4 = jetp4;
    jp4 -= lepp4;
    float dot = lepp4.Vect().Dot( jp4.Vect() );
    float ptrel = lepp4.P2() - dot*dot/jp4.P2();
    ptrel = ptrel>0 ? pow(ptrel,0.5) : 0.0;
    return ptrel;
}

float getPtRel(int id, int idx) {
    // NOTE missing JEC L2L3 (or L1) to correctly compute the closest jet p4
    int closestJetIdx = abs(id) == 11 ? Electron_jetIdx()[idx] : Muon_jetIdx()[idx];
    if (closestJetIdx >= 0){
        LorentzVector closestJet = Jet_p4()[closestJetIdx];
        LorentzVector lepton = abs(id) == 11 ? Electron_p4()[idx] : Muon_p4()[idx];
        return computePtRel(lepton,closestJet);
    } else {
        return 0.;
    }
}

float getPtRatio(int id, int idx) {
    // NOTE missing JEC L2L3 (or L1) to correctly compute the closest jet p4
    int closestJetIdx = abs(id) == 11 ? Electron_jetIdx()[idx] : Muon_jetIdx()[idx];
    if (closestJetIdx >= 0){
        float leppt = abs(id) == 11 ? Electron_pt()[idx] : Muon_pt()[idx];
        return leppt/Jet_pt()[closestJetIdx];
    }
    return 1.;
}
