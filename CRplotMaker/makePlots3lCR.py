import os, sys
import ROOT
# run in batch mode to suppress plot windows
ROOT.gROOT.SetBatch(1)

from MT2PlotMaker import *

def makePlots3lCR(year='run2'):

    bkgnames = [ 'ZGJets', 'DY', 'other_llg' ]
    srNames = ['OffShell', ]
    exts = ['pdf', 'png']

    base_plots = [
        ("njets",True,None,None),
        ("met",True,(0,560.),None,4),
        # ("rlmet",True,(0,560.),None,4),
        ("rlmtZZ",True,(150,650),None,4),
        ("mtWZ",True,(150,1000),None,4),
        ("mWZ",True,(150,1000),None,4),
        ("mtl3",True,None,None,2),
        # ("mllmin",True,None,None,2),
        ("mtWZbins",True,None,None),

        ("ptll",True,None,None,4),
        ("mll",True,(50, 150),None,2),
        ("boson_pt",True,(55,855),None,4),
        ("boson_mass",True,None,None,3),

        ("mllclose",True,None,None,4),
        ("mllx3",True,None,None,4),
        ("mllsfx3",True,None,None,4),

        # ("llid",True,None,None,1),
        ("lepZ1pt",True,None,None),
        ("lepZ2pt",True,None,None),
        ("lepWpt",True,None,None),
        ("lep1pt",True,None,None),
        ("lep2pt",True,None,None),
        ("lep3pt",True,None,None),

        ("dphi_boson_met",True,(1.0,3.2),None),
        ("dphi_lepW_met",True,None,None),
        ("dphi_ZW",True,None,None),
    ]

    jet_plots = [
        ("min_dphijmet",True,(0.2,3.2),None),
        ("min_dphiVjet",True,None,None),
        ("dphi_lljets_met",True,(2.0,3.2),None),
    ]

    input_dir = '../trilepLooper/output/v5_00_3lCR_{}'.format(year)
    output_dir = 'plots{}_3lCR_May10'.format(year[2:])
    bkg_set = [ 'WZ', 'DY', 'ZG', 'ttZ', 'triboson', 'Others' ]
    dataname = 'data_3lskim'.format(year)
    
    systset = ['EW',]
    metsufs = ['_final', ]
    ljsuffs = []
    for jsuf in ['_eq0j', '_eq1j', '_ge2j', '']:
        for lsuf in ['_lle', '_llmu']:
            ljsuffs.append(jsuf+lsuf)

    for ljsuf in ljsuffs:
        outdir = output_dir
        for srn in srNames:
            multsuf = False
            sr = srn if multsuf else srn+ljsuf
            slst = None
            suf, tag = ('', ljsuf) if multsuf else (ljsuf, '')
            for msuf in metsufs:
                plot_set = list(base_plots) if '0j' in suf else base_plots+jet_plots
                MT2PlotMaker(input_dir, bkg_set, dataname, sr, plot_set, output_dir, exts, suffix=msuf+suf, tag=tag, scaleMC=False)
                # MT2PlotMaker(input_dir, bkg_set, dataname, sr, plot_set, output_dir, exts, suffix=msuf+suf, tag=tag, systset=systset, scaleMC=False)

if __name__ == '__main__':

    # makePlots3lCR('2016')
    # makePlots3lCR('2017')
    # makePlots3lCR('2018')

    makePlots3lCR('run2')
