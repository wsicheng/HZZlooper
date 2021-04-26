import os, sys
import ROOT
# run in batch mode to suppress plot windows
ROOT.gROOT.SetBatch(1)

from MT2PlotMaker import *

def makePlots3lCR():

    bkgnames = [ 'ZGJets', 'DY', 'other_llg' ]
    srNames = ['OffShell', ]
    exts = ['pdf', 'png']

    base_plots = [
        # ("njets",True,None,None),
        ("met",True,(0,560.),None,4),
        # ("rlmet",True,(0,560.),None,4),
        ("rlmtZZ",True,(150,650),None,4),
        ("mtWZ",True,(150,650),None,4),
        ("mWZ",True,(150,650),None,4),
        ("mtl3",True,None,None,2),
        # # ("mllmin",True,None,None,2),

        ("ptll",True,None,None,4),
        # ("mll",True,None,None,2),
        ("boson_pt",True,(55,855),None,4),
        ("boson_mass",True,None,None,3),

        # ("llid",True,None,None,1),
        ("lepZ1pt",True,None,None,2),
        ("lepZ2pt",True,None,None,2),
        ("lepWpt",True,None,None,2),

        ("dphi_boson_met",True,(1.0,3.2),None),
        ("dphi_lepW_met",True,None,None),
        ("dphi_ZW",True,None,None),
    ]

    jet_plots = [
        ("min_dphijmet",True,(0.2,3.2),None),
        ("min_dphiVjet",True,None,None),
        ("dphi_lljets_met",True,(2.0,3.2),None),
    ]

    # year = '2018'
    year = 'run2'
    input_dir = '../trilepLooper/output/v4_08_3lCR2_{}'.format(year)
    output_dir = 'plots{}_3lCR_Apr20'.format(year[2:])
    # bkg_set = bkgnames
    bkg_set = [ 'WZ', 'DY', 'ZG', 'ttZ', 'triboson', 'Others' ] # 2016
    dataname = 'data_3lskim'.format(year)
    
    systset = ['EW',]
    metsufs = ['_fullMET', '_final' ]
    ljsuffs = []
    for jsuf in ['_eq0j', '_eq1j', '_ge2j']:
    # for jsuf in [ '']:
        for lsuf in ['_lle', '_llmu']:
            ljsuffs.append(jsuf+lsuf)

    for ljsuf in ljsuffs:
        outdir = output_dir
        for srn in srNames:
            # multsuf = True if '_ll' in ljsuf else False
            multsuf = False
            sr = srn if multsuf else srn+ljsuf
            # slst = [ljsuf.replace('_ll', '_ee'), ljsuf.replace('_ll', '_mumu')] if multsuf else None
            slst = None
            suf, tag = ('', ljsuf) if multsuf else (ljsuf, '')
            for msuf in metsufs:
                plot_set = list(base_plots) if '0j' in suf else base_plots+jet_plots
                MT2PlotMaker(input_dir, bkg_set, dataname, sr, plot_set, output_dir, exts, suffix=msuf+suf, tag=tag, systset=systset, scaleMC=False)

if __name__ == '__main__':

    makePlots3lCR()
