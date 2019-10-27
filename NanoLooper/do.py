#! /usr/bin/env python
from __future__ import print_function, division

import os, sys, math, glob
from multiprocessing.dummy import Pool as ThreadPool
import yaml
import argparse

# Control variables
verbose = 2
poolsize = 8
dryrun = False

# File location
proddir = '/pnfs/iihe/cms/store/group/HZZ2l2nu/Production/'
dir17 = proddir+'2019-10-07_2017/DDF/'
dir16p2 = proddir+'2019-09-05_2016/DDF/'
dir16p1 = proddir+'2019-08-16_2016/DDF/'

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

mc16p1 = [
    'ZZTo2L2Nu',
    'ZZTo2L2Q',
    'ZZTo4L',
    'GluGluToContinToZZTo2e2nu',
    'GluGluToContinToZZTo2mu2nu',
    'WZTo3LNu',
    'WZTo2L2Q',
    'TT',
    'DYJetsToLL_M-50',
    'WJetsToLNu',
    'TTWJetsToLNu',
    'TTZToLLNuNu_M-10',
    'ST_s-channel',
    'ST_t-channel_top',
    'ST_t-channel_antitop',
    'ST_tW_top',
    'ST_tW_antitop',
    'WWW',
    'WWZ',
    'WWTo2L2Nu',
    'WZZ',
    'ZZZ',
    'GGToHToZZTo2E2Nu',
    'GGToZZTo2E2Nu_BSI',
    'GGToHToZZTo2Mu2Nu',
    'GGToZZTo2Mu2Nu_BSI',
]

mc16p2 = [
    'WJetsToLNu_LO',
    'WJetsToLNu_HT-100To200',
    'WJetsToLNu_HT-200To400',
    'WJetsToLNu_HT-400To600',
    'WJetsToLNu_HT-600To800',
    'WJetsToLNu_HT-800To1200',
    'WJetsToLNu_HT-1200To2500',
    'WJetsToLNu_HT-2500ToInf',
    'SinglePhoton',
    'GJets_HT-40To100',
    'GJets_HT-100To200',
    'GJets_HT-200To400',
    'GJets_HT-400To600',
    'GJets_HT-600ToInf',
    'QCD_HT50to100',
    'QCD_HT100to200',
    'QCD_HT200to300',
    'QCD_HT300to500',
    'QCD_HT500to700',
    'QCD_HT700to1000',
    'QCD_HT1000to1500',
    'QCD_HT1500to2000',
    'QCD_HT2000toInf',
    'TTGJets',
    'TGJets',
    'ZGTo2LG',
    'ZNuNuGJets',
    'WGToLNuG',
    'ZGTo2NuG',
    'ZGTo2NuG_PtG-130',
]

ds17 = data17 + mc17

merge_map = {
    'Wjets' : [ 'WJetsToLNu_HT-100To200', 'WJetsToLNu_HT-200To400',
    'WJetsToLNu_HT-400To600', 'WJetsToLNu_HT-600To800', 'WJetsToLNu_HT-800To1200',
    'WJetsToLNu_HT-1200To2500', 'WJetsToLNu_HT-2500ToInf',], 
    'allData' : [ 'SinglePhoton', 'SingleElectron', 'SingleMuon',
                  'DoubleEG', 'DoubleMuon', 'MuonEG' ],
    'diboson' : ['ZZTo2L2Nu', 'ZZTo2L2Q', 'ZZTo4L', 'WZTo3LNu', 'WZTo2L2Q', 'WWTo2L2Nu',],
    'ggTo2L2Nu' : ['GluGluToContinToZZTo2e2nu', 'GluGluToContinToZZTo2mu2nu',],
    'ttbar' : ['TT',],
    'DYjets' : ['DYJetsToLL_M-50',],
    # 'WJetsToLNu_LO',
    # 'WJetsToLNu',
    # 'TTWJetsToLNu',
    # 'TTZToLLNuNu_M-10',
    'signleTop' : ['ST_s-channel', 'ST_t-channel_top', 'ST_t-channel_antitop',
                   'ST_tW_top', 'ST_tW_antitop',],
    'triBoson' : ['WWW', 'WWZ', 'WZZ', 'ZZZ',],
    # 'GGToHToZZTo2E2Nu',
    # 'GGToZZTo2E2Nu_BSI',
    # 'GGToHToZZTo2Mu2Nu',
    # 'GGToZZTo2Mu2Nu_BSI',
    'gjets' : ['GJets_HT-40To100', 'GJets_HT-100To200', 'GJets_HT-200To400',
               'GJets_HT-400To600', 'GJets_HT-600ToInf',],
    'qcd' : [ 'QCD_HT50to100', 'QCD_HT100to200', 'QCD_HT200to300', 'QCD_HT300to500',
              'QCD_HT500to700', 'QCD_HT700to1000', 'QCD_HT1000to1500', 'QCD_HT1500to2000', 
              'QCD_HT2000toInf',],
    # 'TTGJets',
    # 'TGJets',
    # 'ZGTo2LG',
    # 'ZNuNuGJets',
    # 'WGToLNuG',
    # 'ZGTo2NuG',
    # 'ZGTo2NuG_PtG-130',
}


# Merge the files from all input dataset into a single one, for overlap removal consideration
def getFileListData(dsloc, dslist, dataname='data'):
    ds_flists = {}
    ds_flists[dataname] = []

    if dsloc[-1] != '/': ds_loc += '/'
    for dsname in dslist:
        with open(dsloc+dsname+'.yaml', 'r') as ddf:
            ds = yaml.safe_load(ddf)
            ds_flists[dataname] += ds['files']

    return ds_flists
    
def getFileList(dsloc, dslist):
    ds_flists = {}

    if dsloc[-1] != '/': dsloc += '/'
    for dsname in dslist:
        with open(dsloc+dsname+'.yaml', 'r') as ddf:
            ds = yaml.safe_load(ddf)
            ds_flists[dsname] = ds['files']

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
    parser.add_argument('--nomerge', action='store_true', default=False)
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

    dsloc = dir17
    data17_flists = getFileListData(dsloc, data17, 'data17')
    # mc17_flists = getFileList(dsloc, mc17)
    # data_flists = getFileListData(dsloc, test17, 'data17')

    mc16p1 = [ x for x in mc16p1 if x not in mc16p2 ]
    samp16_flists = getFileList(dir16p1, mc16p1+data16p1) 
    samp16_flists.update(getFileList(dir16p2, mc16p2))

    # data_flists = data17_flists
    data_flists = {}
    samp_flists = {}

    samp_flists = samp16_flists
    # samp_flists = getFileList(dir16p1,data16p1)

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

            infilestr = ','.join(data_flists[dsname])
            argstr = '{0} {1} {2}'.format(infilestr, dsname, outdir)
            if not args.nolog: argstr += ' >& {}/logs/log_{}.txt'.format(outdir, dsname) 
            arglist.append(argstr)

            if verbose >= 1:
                print( '>>> Running data sample {} with 1 jobs!'.format(dsname) )

    # The samples that can be separated
    for dsname, flist in samp_flists.items():
        nfile_per_job = 20
        njobs = int(math.ceil(len(flist) / nfile_per_job))
        njobs_total += njobs

        for i in range(njobs):
            infilestr = ','.join(flist[ i*nfile_per_job : (i+1)*nfile_per_job ])
            argstr = '{0} {1}_{3} {2}'.format(infilestr, dsname, outdir, i)
            if not args.nolog: argstr += ' >& {}/logs/log_{}.txt'.format(outdir, dsname)
            arglist.append(argstr)

        if verbose >= 1:
            print( '>>> Running sample {} with {} jobs!'.format(dsname, njobs) )

    for rs in pool.imap_unordered(runHZZlooper, arglist):
        if rs == 0: njobs_done += 1
        else: print( '>>> Finish running {}/{} jobs!'.format(njobs_done, njobs_total) )

        if verbose >= 1:
            print( '>>> Finish running {}/{} jobs!'.format(njobs_done, njobs_total) )

            
    if njobs_done == njobs_total:
        print( 'All {} job done!'.format(njobs_done) )
        if not args.nomerge:
            # mergeOutputHists(outdir,'',['allData',])
            mergeOutputHists(outdir)
    else:
        print( '{} jobs finished out of {}!'.format(njobs_done, njobs_total) )
        

    print( 'Test end!' )
        
