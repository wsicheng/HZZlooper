#ifndef ISOLATIONTOOLS_H
#define ISOLATIONTOOLS_H
#include "Nano.h"


float getPtRel(int id, int idx);
float getPtRatio(int id, int idx);
float computePtRel(const LorentzVector& lepp4, const LorentzVector& jetp4);
bool passMultiIso(int id, int idx, float cutMiniIso, float cutPtRatio, float cutPtRel);

#endif
