# MT2PlotUtils.py
# various plotting utilities

def GetColor(sample):
    if 'ttbar'   in sample: return 870  # kAzure+10
    if 'singleT' in sample: return 625  # kRed-7
    if 'Vjets'   in sample: return 798  # kOrange-2
    if 'rare'    in sample: return 614  # kMagenta-2

    if 'wjets_lnu_0j'   in sample.lower(): return 625  # kRed-7
    if 'tt2l'    in sample.lower(): return 625  # kRed-7
    if 'tt_'     in sample.lower(): return 625  # kRed-7
    if 'zz2l'    in sample.lower(): return 861  # kAzure+1
    if 'qqzz'    in sample.lower(): return 861  # kAzure+1
    if 'ww2l'    in sample.lower(): return 798  # kOrange-2
    if 'ttw'     in sample.lower(): return 844  # kTeal+4
    if 'dy_m10'  in sample.lower(): return 844  # kTeal+4
    if 'dy_m50'  in sample.lower(): return 842  # kTeal+2
    if 'dy'      in sample.lower(): return 842  # kTeal+2
    if 'tt1l'    in sample.lower(): return 861  # kAzure+1
    if 'wjets'   in sample.lower(): return 798  # kOrange-2
    if "st_tw"   in sample.lower(): return 614  # kMagenta-2
    if "qcd"     in sample.lower(): return 921  # kGray+1
    if "zgjets"  in sample.lower(): return 425  # kCyan-7
    if "zgto2nu" in sample.lower(): return 425  # kCyan-7
    if "zgto2lg" in sample.lower(): return 425  # kCyan-7
    if "topg"    in sample.lower(): return 625  # kRed-7
    if "tgjets"  in sample.lower(): return 625  # kRed-7
    if "wg"      in sample.lower(): return 876  # kViolet-4
    if "gjets"   in sample.lower(): return 901  # kPink+1  
    if "zjets"   in sample.lower(): return 842  # kTeal+2
    if "qqwz"    in sample.lower(): return 876  # kViolet-4
    if "wlnu"    in sample.lower(): return 798  # kOrange-2


    if 'data'    in sample: return 14   # Gray
    if 'lostlep' in sample: return 866  # kAzure+6
    if 'FromW'   in sample: return 796  # kOrange-4
    if 'FromTop' in sample: return 625  # kRed-7
    if 'ZToNuNu' in sample: return 613  # kMagenta-3

    if "2lep"    in sample: return 866  # kAzure+6
    if "1lepW"   in sample: return 796  # kOrange-4
    if "1lepTop" in sample: return 625  # kRed-7
    if "Znunu"   in sample: return 872  # kViolet-8
    if 'others'  in sample: return 625  # kRed-7
    if 'Others'  in sample: return 625  # kRed-7

    if 'ZG'       == sample: return 796  # kOrange-4
    if 'ZZ'       == sample: return 861  # kAzure+1
    if 'WZ'       == sample: return 425  # kCyan-7
    if 'nonres'   == sample: return 872  # kViolet-8
    if 'diboson'  == sample: return 796  # kOrange-4
    if 'triboson' == sample: return 872  # kViolet-8

    if "bevt"  in sample: return 842  # kTeal+2
    if "cevt"  in sample: return 798  # kOrange-2
    if "lfevt" in sample: return 425  # kCyan-7


    return 870 # for everything else

def GetSampleName(sample):
    names = {
        "lostlepFromCRs": "Lost Lepton",
        "lostlep": "Lost Lepton",
        "tt2l" : "t#bar{t}#rightarrow2 #font[12]{l}" ,
        "tt1l" : "t#bar{t}#rightarrow1 #font[12]{l}" ,
        "ttbar" : "t#bar{t}+Jets" ,
        "Vjets": "V+Jets",
        "singleT": "single top",
        "rare": "t#bar{t}V+VV",
        "TTZToLLNuNu": "t#bar{t}Z(#font[12]{#nu#nu})+Jets",
        "TT2L2Nu" : "t#bar{t}#rightarrow2#font[12]{l}2#font[12]{#nu}" ,
        "WW2L2Nu" : "WW#rightarrow2#font[12]{l}2#font[12]{#nu}" ,
        "ZZ2L2Nu" : "ZZ#rightarrow2#font[12]{l}2#font[12]{#nu}" ,
        "DY_" : "DY ",
        "HTbinned" : " ",
        "ST_tW" : "tW#rightarrow2#font[12]{l}2#font[12]{#nu}",
        "_16" : " 16" ,
        "_17" : " 17" ,
        "_18" : " 18" ,
        "_run2" : "" ,
        "TTTo2L2Nu": "t#bar{t}#rightarrow2#font[12]{l}2#font[12]{#nu}",
        "ggHZZ_all" : "gg#rightarrowZZ total (#Gamma_{H}=#Gamma_{H}^{SM})",
        "ggHZZ_sigonly" : "gg#rightarrowZZ sig. (#Gamma_{H}=#Gamma_{H}^{SM})",
        # "allData" : "allData" ,
        # "allBkg " : "allBkg" ,
        "_mrs0" : " no METRes corr." ,
        "_mrs2" : " w/ METRes corr." ,
        "2lep" : "t#bar{t}/tW#rightarrow2 #font[12]{l}",
        "1lepW" : "1#font[12]{l} from W",
        "1lepTop" : "1#font[12]{l} from top",
        "Znunu" : "Z#rightarrow #nu#nu",
        "bevt"  : "Events with b-quarks",
        "cevt"  : "Events with c-quarks",
        "lfevt" : "Events with LF",
        "_datadriven" : "",
    }

    # use the above name if defined, otherwise use sample itself
    for k, v in names.items():
        if k in sample:
            sample = sample.replace(k, v)

    # return names.get(sample,sample)
    return sample

def GetVarName(var):
    names = {"ht": "H_{T}",
             "met": "E_{T}^{miss}",
             "mt2": "M_{T2}",
             "mt" : "M_{T}",
             "metbins": "#slash{E}_{T}",
             "mt2bins": "M_{T2}",
             "nJet30": "N(jet)",
             "nBJet20": "N(b jet)",
             "njets": "N(jet)",
             "nbjets": "N(b-jet)",
             "leppt": "p_{T}(lep)",
             "lepeta": "#eta(lep)",
             "lep1pt": "p_{T}(lep1)",
             "lep1eta": "#eta(lep1)",
             "lep2pt": "p_{T}(lep2)",
             "lep2eta": "#eta(lep2)",
             "nlepveto": "N(veto lep)",
             "zllmass": "m_{#font[12]{ll}}",
             "gammaPt": "P_{T}(#gamma)",
             "gammaEta": "#eta(#gamma)",
             "jet1pt": "Leading jet p_{T}",
             "E_{T}^{miss}" : "p_{T}^{miss}",
             }

    if "SigmaIetaIeta" in var:
        return "Photon #sigma_{i#eta i#eta}"

    # use the above name if defined, otherwise use var itself
    return names.get(var,var)

def GetUnit(vn):
    noUnit = ["nJet","nBJet","eta","nlep","drMin","SigmaIetaIeta"]
    for s in noUnit:
        if s.lower() in vn.lower():
            return None

    return "GeV"

def GetCRName(cr, outn=None):
    names = {
        "OffShell" : "",
        "_llg"   : " ee/#mu#mu+#gamma",
        "_eeg"   : " ee+#gamma",
        "_lle"   : " #font[12]{ll}+e",
        "_llmu"  : " #font[12]{ll}+#mu",
        "_mumug" : " #mu#mu+#gamma",
        "_gamma" : " Single-#gamma",
        "_eq0j" : "N_{j}=0",
        "_eq1j" : "N_{j}=1",
        "_eq2j" : "N_{j}=2",
        "_ge2j" : "N_{j}#geq2",
    }

    cr = cr.split('_')
    cr = '_'+', _'.join(reversed(cr))

    # use the above name if defined, otherwise use cr itself
    for k, v in names.items():
        if k in cr:
            cr = cr.replace(k, v)

    cr = cr.replace(', _', '')
    if 'metcrtd' in outn: cr += ", p_{T}^{miss} smeared"
    # elif 'nometres' not in outn: cr += ", p_{T}^{miss} smeared"

    return cr

def GetSubtitles(dirname, vn):

    subtitle = ["p_{T}^{#font[12]{ll}}#geq 55GeV, METSEL", "76 < m_{#font[12]{ll}} < 106 GeV", ""]
    if "_gamma" in vn:
        subtitle = ["p_{T}^{#gamma}#geq 55GeV, METSEL", "#Delta#phi(#gamma,p_{T}^{miss})>1.0", ""]
    elif "_llg" in vn:
        subtitle = ["p_{T}^{#gamma}#geq 55GeV, METSEL", "Analyssis #Delta#phi(#gamma,p_{T}^{miss}+p_{T}^{#font[12]{ll}}) sel.", "76<M_{#font[12]{ll}}<106 GeV" ]
    elif "_mumug" in vn:
        subtitle = ["p_{T}^{#gamma}#geq 55GeV, METSEL", "Analyssis #Delta#phi(#gamma,p_{T}^{miss}+p_{T}^{#mu#mu}) sel.", "76<M_{#mu#mu}<106 GeV" ]
    elif "_eeg" in vn:
        subtitle = ["p_{T}^{#gamma}#geq 55GeV, METSEL", "Analyssis #Delta#phi(#gamma,p_{T}^{miss}+p_{T}^{ee}) sel.", "76<M_{ee}<106 GeV" ]

    if "_eq0j" in vn:
        subtitle[1] = "Analysis #Delta#phi sel."
    elif "_eq1j" in vn:
        subtitle[1] = "Analysis #Delta#phi sel., 0b"
    elif "_eq2j" in vn:
        subtitle[1] = "Analysis #Delta#phi sel., 0b"
    elif "_ge2j" in vn:
        subtitle[1] = "Analysis #Delta#phi sel., 0b"


    if "_fullMET" in vn:
        subtitle[0] = subtitle[0].replace("METSEL", "p_{T}^{miss} > 0 GeV")
    elif "_met_lt125" in vn:
        subtitle[0] = subtitle[0].replace("METSEL", "0 < p_{T}^{miss} < 125 GeV")
    elif "_metlt125" in vn:
        subtitle[0] = subtitle[0].replace("METSEL", "p_{T}^{miss} < 125 GeV")
    elif "_metlt80" in vn:
        subtitle[0] = subtitle[0].replace("METSEL", "0 < p_{T}^{miss} < 80 GeV")
    elif "_met80to125" in vn:
        subtitle[0] = subtitle[0].replace("METSEL", "80 < p_{T}^{miss} < 125 GeV")
    elif "_met140to200" in vn:
        subtitle[0] = subtitle[0].replace("METSEL", "p_{T}^{miss}: [140, 200) GeV")
    elif "_metge125" in vn:
        subtitle[0] = subtitle[0].replace("METSEL", "p_{T}^{miss} #geq 125 GeV")
    elif "_metge200" in vn:
        subtitle[0] = subtitle[0].replace("METSEL", "p_{T}^{miss} #geq 200 GeV")
    elif "_metlt50" in vn:
        subtitle[0] = subtitle[0].replace("METSEL", "0 < p_{T}^{miss} < 50 GeV")
    elif "_metge50" in vn:
        subtitle[0] = subtitle[0].replace("METSEL", "p_{T}^{miss} > 50 GeV")
    elif "_final" in vn and '_ge2j' in vn:
        subtitle[0] = subtitle[0].replace("METSEL", "p_{T}^{miss} #geq140 GeV")
    elif "_final" in vn:
        subtitle[0] = subtitle[0].replace("METSEL", "p_{T}^{miss} #geq125 GeV")
    else:
        subtitle[0] = subtitle[0].replace("METSEL", "")
       
    if "_llg" in vn:
        subtitle[0] = subtitle[0].replace("p_{T}^{miss}", "p_{T}^{miss}+p_{T}^{#font[12]{ll}}")

    return subtitle

def RebinAll(h_bkg_vec, r):
    for h in h_bkg_vec:
        h.Rebin(r)
