#! /usr/bin/env python

import os, sys
from multiprocessing.dummy import Pool as ThreadPool
import yaml
import argparse


dir17 = '/pnfs/iihe/cms/store/group/HZZ2l2nu/Production/2019-10-07_2017/DDF/'

ds17 = [
    'DoubleMuon',
    'SingleMuon',
    'DoubleEG',
    'MuonEG',
    'ZZTo2L2Nu_13TeV_powheg_pythia8',
    'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
    'WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8',
    'WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8',
    'WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8',
    'WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8',
    'WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8',
    'TTZToLLNuNu_M-10_TuneCP5_PSweights_13TeV-amcatnlo-pythia8',
    'ST_t-channel_top_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
    'TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8',
    'WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8',
    'WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8',
    'ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
    'GluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8',
    'GluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8',
    'TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8',
    'WZTo3LNu_13TeV-powheg-pythia8',
    'WWTo2L2Nu_NNPDF31_TuneCP5_PSweights_13TeV-powheg-pythia8',
    'TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8',
    'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8',
    'WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8',
    'WZZ_TuneCP5_13TeV-amcatnlo-pythia8',
    'ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8',
    'TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8',
    'WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8',
    'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8',
    'DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8',
    'WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8',
    'ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
    'SingleElectron',
    'ST_s-channel_4f_leptonDecays_TuneCP5_PSweights_13TeV-amcatnlo-pythia8',
    'WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8',
    'ZZTo4L_13TeV_powheg_pythia8',
]

def run_HZZlooper(args):
    
    print '  ./runHZZlooper '+args
    os.system('./runHZZlooper '+args)


if __name__ == '__main__':

    print 'Test!'

    parser = argparse.ArgumentParser('Run the simple HZZlooper')
    parser.add_argument('outdir', nargs='?', default='output/temp')

    args = parser.parse_args()

    outdir = args.outdir

    with open(dir17+'DoubleMuon'+'.yaml', 'r') as ddf:
        ds = yaml.safe_load(ddf)
        ff = ds['files'][0]
        print ff

        pool = ThreadPool(8)

        arglist = []
        arglist.append('{} {} {}'.format(ff, ds['stem'], outdir))

        pool.imap_unordered(run_HZZlooper, arglist)


    print 'Test end!'
        
