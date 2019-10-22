#! /usr/bin/env python
from __future__ import print_function, division

import os, sys
from multiprocessing.dummy import Pool as ThreadPool
import yaml
import argparse

# Control variables
verbose = 3
poolsize = 8
dryrun = False

# File location
proddir = '/pnfs/iihe/cms/store/group/HZZ2l2nu/Production/'
dir17 = proddir+'2019-10-07_2017/DDF/'
dir16 = proddir+'2019-09-05_2016/DDF/'
# dir16 = proddir+'2019-08-16_2016/DDF/'

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


ds17 = data17 + mc17

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
    
def getFileListMC(dsloc, dslist):
    ds_flists = {}

    if dsloc[-1] != '/': dsloc += '/'
    for dsname in dslist:
        with open(dsloc+dsname+'.yaml', 'r') as ddf:
            ds = yaml.safe_load(ddf)
            ds_flists[dsname] = ds['files']

    return ds_flists

def run_HZZlooper(args):
    
    if verbose >= 3:
        print( '  ./runHZZlooper '+args)

    if not dryrun:
        return os.system('./runHZZlooper '+args)

    return 0

if __name__ == '__main__':

    print( 'Test!' )

    parser = argparse.ArgumentParser('Run the simple HZZlooper')
    parser.add_argument('outdir', nargs='?', default='output/temp')
    parser.add_argument('-dr', '--dryrun', action='store_true', default=False)
    parser.add_argument('-fg', '--nolog', action='store_true', default=False)

    args = parser.parse_args()

    outdir = args.outdir
    dryrun = args.dryrun

    dsloc = dir17

    # data_flists = getFileListData(dsloc, data17, 'data17')
    # mc_flists = getFileListMC(dsloc, mc17)
    data_flists = getFileListData(dsloc, test17, 'data17')
    mc_flists = {}

    os.system('make -j 12')
    os.system('mkdir -p {}/logs'.format(outdir))

    pool = ThreadPool(poolsize)
    njobs_total = 0
    njobs_done = 0

    arglist = []

    for dsname in ['data17',]:
        njobs_total += 1

        infilestr = ','.join(data_flists[dsname])
        argstr = '{0} {1} {2}'.format(infilestr, dsname, outdir)
        if not args.nolog: argstr += ' >& {}/logs/log_{}.txt'.format(outdir, dsname) 
        arglist.append(argstr)

        if verbose >= 1:
            print( '>>> Running data sample {} with 1 jobs!'.format(dsname) )

    # The MC samples can be separated
    for dsname, flist in mc_flists.items():
        nfile_per_job = 20
        njobs = int(len(flist) / nfile_per_job) + 1
        njobs_total += njobs

        for i in range(njobs):
            infilestr = ','.join(flist[ i*nfile_per_job : (i+1)*nfile_per_job ])
            argstr = '{0} {1}_{3} {2}'.format(infilestr, dsname, outdir, i)
            if args.nolog: argstr += ' >& {}/logs/log_{}.txt'.format(outdir, dsname)
            arglist.append(argstr)

        if verbose >= 1:
            print( '>>> Running sample {} with {} jobs!'.format(dsname, njobs) )

    for rs in pool.imap_unordered(run_HZZlooper, arglist):
        if rs == 0: njobs_done += 1
            
    if njobs_done == njobs_total:
        print( 'All {} job done!'.format(njobs_done) )
    else:
        print( '{} jobs finished out of {}!'.format(njobs_done, njobs_total) )
        

    print( 'Test end!' )
        
