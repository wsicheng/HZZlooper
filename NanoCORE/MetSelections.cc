#include <iostream>

#include "Math/VectorUtil.h"
#include "MetSelections.h"
#include "Config.h"

using namespace tas;


bool passesMETfilters(bool isData){

    if (gconf.year == 2016) {
        if(isData) {
        }
    } else if (gconf.year == 2017) {
        if (!Flag_globalTightHalo2016Filter()) return false;
        if (!Flag_HBHENoiseFilter()) return false;
        if (!Flag_HBHENoiseIsoFilter()) return false;
        if (!Flag_EcalDeadCellTriggerPrimitiveFilter()) return false;
        if (!Flag_eeBadScFilter()) return false;
        if (!Flag_goodVertices()) return false;
        if (!Flag_chargedHadronTrackResolutionFilter()) return false;
        if (!Flag_muonBadTrackFilter()) return false;
        if(isData) {
        }
    } else if (gconf.year == 2018) {
        if(isData) {
        }
    } else {
        std::cout << "uh, I don't know what year this is" << std::endl;
    }

    // otherwise, good
    return true;
}

