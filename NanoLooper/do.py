#! /usr/bin/env python
from __future__ import print_function, division

import os, sys, math, glob
from multiprocessing.dummy import Pool as ThreadPool
import yaml
import argparse
from metis.Sample import DBSSample, DirectorySample

# Control variables
verbose = 2
poolsize = 8
dryrun = False

# File location
proddir = '/pnfs/iihe/cms/store/group/HZZ2l2nu/Production/'
dir17 = proddir+'2019-11-01_2017/DDF/'
dir16p2 = proddir+'2019-09-05_2016/DDF/'
dir16p1 = proddir+'2019-08-16_2016/DDF/'

proddir = '/hadoop/cms/store/user/sicheng/NanoAOD'
dir16v7 = proddir+'/mc2016_v7'

# Dataset names
test17 = [
    'DoubleMuon',
    'DoubleEG',
    'MuonEG',
]

data17 = [
    'DoubleMuon',
    'SingleMuon',
    'DoubleEG',
    'SingleElectron',
    'MuonEG',
]
 
mc17 = [
    'ZZTo2L2Q',
    'ZZTo2L2Nu',
    'ZZTo4L',
    'GluGluToContinToZZTo2e2nu',
    'GluGluToContinToZZTo2mu2nu',
    'WWTo2L2Nu',
    'WZTo3LNu',
    'WZTo2L2Q',
    'TT',
    'DYJetsToLL_M-50',
    'TTWJetsToLNu',
    'TTZToLLNuNu_M-10',
    'ST_s-channel',
    'ST_t-channel_top',
    'ST_t-channel_antitop',
    'ST_tW_top',
    'ST_tW_antitop',
    'WWW',
    'WWZ',
    'WZZ',
    'ZZZ',
    'WJetsToLNu_LO',
    'WJetsToLNu_HT-100To200',
    'WJetsToLNu_HT-200To400',
    'WJetsToLNu_HT-400To600',
    'WJetsToLNu_HT-600To800',
    'WJetsToLNu_HT-800To1200',
    'WJetsToLNu_HT-1200To2500',
    'WJetsToLNu_HT-2500ToInf',
]

data16p1 = [
    'SingleElectron',
    'SingleMuon',
    'DoubleEG',
    'DoubleMuon',
    'MuonEG',
]

ds17 = data17 + mc17

dsinfos = {
    # diBoson
    # '/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM|3.05|1.00262|1',
    # 'DYJetsToLL_M-50' : '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20_ext2-v1/NANOAODSIM|130939668', # total 193215674
    # "GJets_HT-40To100" : "/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
    # "GJets_HT-100To200" : "/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_4cores5k_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
    # "GJets_HT-200To400" : "/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
    # "GJets_HT-400To600" : "/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
    # "GJets_HT-600ToInf" : "/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20_ext1-v1/NANOAODSIM",

    # "GJets_DR-0p4_HT-600ToInf" : "/GJets_DR-0p4_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
    # "GJets_DR-0p4_HT-400To600" : "/GJets_DR-0p4_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
    # "GJets_DR-0p4_HT-200To400" : "/GJets_DR-0p4_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
    # "GJets_DR-0p4_HT-100To200" : "/GJets_DR-0p4_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",

    # 'data_egamma_2018A' : '/EGamma/Run2018A-02Apr2020-v1/NANOAOD',
    # 'data_egamma_2018B' : '/EGamma/Run2018B-02Apr2020-v1/NANOAOD',
    # 'data_egamma_2018C' : '/EGamma/Run2018C-02Apr2020-v1/NANOAOD',
    # 'data_egamma_2018D' : '/EGamma/Run2018D-02Apr2020-v1/NANOAOD',

    # 'data_egamma_2018A' : '/EGamma/Run2018A-Nano25Oct2019-v1/NANOAOD',
    # 'data_egamma_2018B' : '/EGamma/Run2018B-Nano25Oct2019-v1/NANOAOD',
    # 'data_egamma_2018C' : '/EGamma/Run2018C-Nano25Oct2019-v1/NANOAOD',
    # 'data_egamma_2018D' : '/EGamma/Run2018D-Nano25Oct2019-v1/NANOAOD',

    # 'data_photon_2016B' : "/SinglePhoton/Run2016B_ver2-Nano25Oct2019_ver2-v1/NANOAOD",
    'data_photon_2016C' : "/SinglePhoton/Run2016C-Nano25Oct2019-v1/NANOAOD",
    'data_photon_2016D' : "/SinglePhoton/Run2016D-Nano25Oct2019-v1/NANOAOD",
    'data_photon_2016E' : "/SinglePhoton/Run2016E-Nano25Oct2019-v1/NANOAOD",
    'data_photon_2016F' : "/SinglePhoton/Run2016F-Nano25Oct2019-v1/NANOAOD",
    'data_photon_2016G' : "/SinglePhoton/Run2016G-Nano25Oct2019-v1/NANOAOD",
    'data_photon_2016H' : "/SinglePhoton/Run2016H-Nano25Oct2019-v1/NANOAOD",

    # 'data_photon_2016B' : "/SinglePhoton/Run2016B-02Apr2020_ver2-v1/NANOAOD",
    # 'data_photon_2016C' : "/SinglePhoton/Run2016C-02Apr2020-v1/NANOAOD",
    # 'data_photon_2016D' : "/SinglePhoton/Run2016D-02Apr2020-v1/NANOAOD",
    # 'data_photon_2016E' : "/SinglePhoton/Run2016E-02Apr2020-v1/NANOAOD",
    # 'data_photon_2016F' : "/SinglePhoton/Run2016F-02Apr2020-v1/NANOAOD",
    # 'data_photon_2016G' : "/SinglePhoton/Run2016G-02Apr2020-v1/NANOAOD",
    # 'data_photon_2016H' : "/SinglePhoton/Run2016H-02Apr2020-v1/NANOAOD",

}

merge_map = {
    'Wjets' : [ 'WJetsToLNu_HT-100To200', 'WJetsToLNu_HT-200To400',
                'WJetsToLNu_HT-400To600', 'WJetsToLNu_HT-600To800', 'WJetsToLNu_HT-800To1200',
                'WJetsToLNu_HT-1200To2500', 'WJetsToLNu_HT-2500ToInf',], 
    # 'WJetsToLNu_LO', 'WJetsToLNu',
    'ZZto2L2X' : ['ZZTo2L2Nu', 'ZZTo2L2Q', 'ZZTo4L',],
    'WZ' : ['WZTo3LNu', 'WZTo2L2Q', ],
    'contin2L2Nu' : ['GluGluToContinToZZTo2e2nu', 'GluGluToContinToZZTo2mu2nu',],
    'ggZZto2L2Nu' : ['GGToZZTo2E2Nu_BSI', 'GGToZZTo2Mu2Nu_BSI',], # what's the difference?
    # 'ttbar' : ['TT',],
    'DYjets' : ['DYJetsToLL_M-50',],
    'ttV' : ['TTWJetsToLNu', 'TTZToLLNuNu_M-10',],
    'singleTop' : ['ST_s-channel', 'ST_t-channel_top', 'ST_t-channel_antitop',
                   'ST_tW_top', 'ST_tW_antitop',],
    'triBoson' : ['WWW', 'WWZ', 'WZZ', 'ZZZ',],
    'ggH' : ['GGToHToZZTo2E2Nu', 'GGToHToZZTo2Mu2Nu',],
    'gjets' : ['GJets_HT-40To100', 'GJets_HT-100To200', 'GJets_HT-200To400',
               'GJets_HT-400To600', 'GJets_HT-600ToInf',],
    'qcd' : [ 'QCD_HT50to100', 'QCD_HT100to200', 'QCD_HT200to300', 'QCD_HT300to500',
              'QCD_HT500to700', 'QCD_HT700to1000', 'QCD_HT1000to1500', 'QCD_HT1500to2000', 
              'QCD_HT2000toInf',],
    'gammaX' : ['TTGJets', 'TGJets', 'ZGTo2LG', 'ZNuNuGJets', 'WGToLNuG', 'ZGTo2NuG', 'ZGTo2NuG_PtG-130',],

    'allData' : [ 'SinglePhoton', 'SingleElectron', 'SingleMuon',
                  'DoubleEG', 'DoubleMuon', 'MuonEG' ],
    'TopnW' : [ 'TT', 'Wjets', 'WWTo2L2Nu',],
    'ZZ' : ['ZZto2L2X', 'ggZZto2L2Nu'],

}


# Merge the files from all input dataset into a single one, for overlap removal consideration
def getFileListData(dsloc, dslist, dataname='data'):
    ds_flists = {}
    ds_flists[dataname] = []

    if dsloc[-1] != '/': ds_loc += '/'
    for dsname in dslist:
        with open(dsloc+dsname+'.yaml', 'r') as ddf:
            ds = yaml.safe_load(ddf)
            ds_flists[dataname] += (ds['files'], 0)

    return ds_flists
    
def getSampList(dsloc):
    if dsloc[-1] != '/': dsloc += '/'
    return [ fn.split('/')[-1][:-5] for fn in glob.glob(dsloc+'*.yaml')]

def getFileList(dsloc, dslist=None):
    if dsloc[-1] != '/': dsloc += '/'
    if dslist == None:
        dslist = getSampList(dsloc)

    ds_flists = {}
    for dsname in dslist:
        with open(dsloc+dsname+'.yaml', 'r') as ddf:
            ds = yaml.safe_load(ddf)
            if 'num_events' in ds:
                ds_flists[dsname] = (ds['files'], ds['num_events'])
            else:
                ds_flists[dsname] = (ds['files'], 0)

    return ds_flists

def getDBSFileList(dsinfos, dslist=None):
    if dslist == None:
        dslist = dsinfos.keys()
    ds_flists = {}
    for dsname in dslist:
        dsinfo = dsinfos[dsname].split('|')
        dataset = dsinfo[0]
        sample = DBSSample(dataset=dataset)
        nevts = dsinfo[1] if len(dsinfo) > 1 else sample.get_nevents()
        ds_flists[dsname] = ([ef.get_name() for ef in sample.get_files()], nevts)

    return ds_flists

def getDirFileList(findir, dslist=None):
    ds_flists = {}
    
    dins = glob.glob(findir+'/*') 
    dsel = dins if dslist==None else [ d for d in dins if any(k in d for k in dslist) ]
                
    for d in dsel:
        with open(d+'/info.yml') as finfo:
            info = yaml.safe_load(finfo)
        fins = glob.glob(d+'/*.root') 
        # dsname = info.get('dsname', d.replace(findir,'').replace('__','/')+'/NANOAODSIM')
        dsname = info.get('dsname', reduce(sum, [ k for k in dslist if k in d ]) if dslist!=None else d.replace(findir+'/',''))
        nevts = info.get('nevents', -1)
        ds_flists[dsname] = (fins, nevts)

    return ds_flists

def runHZZlooper(args):
    
    if verbose >= 3:
        print( '  ./runHZZlooper '+args)

    if not dryrun:
        return os.system('./runHZZlooper '+args)

    return 0

def mergeOutputHists(outdir, suf='', lst_samp=None):
    # os.system('pushd {}'.format(outdir))
    pwd = os.getcwd()
    os.chdir(outdir)

    for target, finlist in merge_map.items():
        if lst_samp != None and target not in lst_samp: continue
        fins = [ f+'*.root' for f in finlist ]
        fins = [ glob.glob(f+'*.root') for f in finlist]
        # print( fins, sum(fins, []))
        if len(fins) == 0: continue
        print( 'hadd -f {}.root {}'.format(target+suf, ' '.join(sum(fins, []))) )
        if not dryrun:
            os.system('hadd -f {}.root {} > /dev/null'.format(target+suf, ' '.join(sum(fins, []))))

    os.chdir(pwd)

if __name__ == '__main__':

    print( 'Test!' )

    parser = argparse.ArgumentParser('Run the simple HZZlooper')
    parser.add_argument('outdir', nargs='?', default='temp')
    parser.add_argument('-v', '--verbose', default=1)
    parser.add_argument('-j', '--ncore', default=12)
    parser.add_argument('-mo', '--merge_only', default=False)
    parser.add_argument('-dr', '--dryrun', action='store_true', default=False)
    parser.add_argument('-fg', '--nolog', action='store_true', default=False)
    parser.add_argument('--merge', action='store_true', default=False)
    parser.add_argument('--nice', default=10)

    args = parser.parse_args()

    outdir = 'output/'+args.outdir
    dryrun = args.dryrun
    verbose = args.verbose
    poolsize = args.ncore
    os.nice(args.nice)

    if args.merge_only != False:
        lst_samp = args.merge_only.split(',') if len(args.merge_only) > 0 else None
        if len(lst_samp) == 1 and lst_samp[0] == 'all': lst_samp = None
        mergeOutputHists(outdir, lst_samp=lst_samp)
        exit(0)

    ## Define amples to run over

    # samp17_flists = getFileList(dir17, data17+mc17)

    # mc16p1 = [ x for x in mc16p1 if x not in mc16p2 ]
    samp16_flists = getFileList(dir16p1)
    samp16_flists.update( getFileList(dir16p2) )

    # data_flists = data17_flists
    data_flists = {}
    samp_flists = {}

    # samp_flists = samp16_flists
    samp_flists = getDBSFileList(dsinfos)
    # samp_flists = getFileList(dir16p1, [ 'MuonEG',])
    # samp_flists = getDirFileList(dir16v7, [ 'ZGTo2',])

    rv = os.system('make -j 12')
    if rv != 0: exit()          # quit if make is not successful

    os.system('mkdir -p {}/logs'.format(outdir))

    pool = ThreadPool(poolsize)
    njobs_total = 0
    njobs_done = 0

    arglist = []

    combine_data = False
    # The samples that has to be merged to avoid overlap
    if combine_data:
        for dsname in ['data17',]:
            njobs_total += 1

            infilestr = ','.join(data_flists[dsname][0])
            argstr = '{0} {1} {2}'.format(infilestr, dsname, outdir)
            if not args.nolog: argstr += ' >& {}/logs/log_{}.txt'.format(outdir, dsname) 
            arglist.append(argstr)

            if verbose >= 1:
                print( '>>> Running data sample {} with 1 jobs!'.format(dsname) )

    # The samples that can be separated
    for dsname, (flist, nevt) in samp_flists.items():
        nfile_per_job = 10
        njobs = int(math.ceil(len(flist) / nfile_per_job))
        njobs_total += njobs

        for i in range(njobs):
            infilestr = ','.join(flist[ i*nfile_per_job : (i+1)*nfile_per_job ])
            argstr = '{0} {1}_{3} {2} {4}'.format(infilestr, dsname, outdir, i, nevt)
            if not args.nolog: argstr += ' >& {}/logs/log_{}_{}.txt'.format(outdir, dsname, i)
            arglist.append(argstr)

        if verbose >= 1:
            print( '>>> Running sample {} with {} jobs!'.format(dsname, njobs) )

    for rs in pool.imap_unordered(runHZZlooper, arglist):
        if rs == 0: njobs_done += 1
        else: print( '>>> The {}-th job return with {}!'.format(njobs_done, rs) )

        if verbose >= 1:
            print( '>>> Finish running {}/{} jobs!'.format(njobs_done, njobs_total) )

            
    if njobs_done >= njobs_total-1:
        print( 'All {} job done!'.format(njobs_done) )
        if args.merge:
            # mergeOutputHists(outdir,'',['allData',])
            mergeOutputHists(outdir)
    else:
        print( '{} jobs finished out of {}!'.format(njobs_done, njobs_total) )
        

    print( 'Test end!' )
        
