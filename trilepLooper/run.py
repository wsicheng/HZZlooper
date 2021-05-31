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
osoutdir = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output";
# osoutdir = "root://cmsxrootd.fnal.gov//store/user/usarica/Offshell_2L2Nu/Worker/output";

skimver  = "v4_08";
skimdate = "/SkimTrees/210315/AK4Jets_WithPUJetId_NoTightLeptonJetId_ParticleCleaned";

skimver  = "v5_00";
skimver  = "v5_01";
skimdate = "/SkimTrees/210504/AK4Jets_WithPUJetId_NoTightLeptonJetId_ParticleCleaned";

datasuf  = "/PFMET_WithXY_WithPartMomCorr_P4Preserved/";
mcsuf    = "/PFMET_WithXY_NoJER_WithPartMomCorr_P4Preserved_ResolutionCorrected/";

def runOSlooper(args):
    if verbose >= 3:
        print( '  ./runOSlooper '+args)
    if not dryrun:
        return os.system('./runOSlooper '+args)
    return 0

def mergeOutputHists(outdir, merge_map, suf='', lst_samp=None, clean_unmerge=False):
    pwd = os.getcwd()
    os.chdir(outdir)

    for target, finlist in merge_map.items():
        if lst_samp != None and target not in lst_samp: continue
        fins = [ f+'*.root' for f in finlist ]
        fins = [ glob.glob(f+'*.root') for f in finlist]
        if len(fins) == 0: continue
        print( 'hadd -f {}.root {}'.format(target+suf, ' '.join(sum(fins, []))) )
        # if not dryrun:
        os.system('hadd -f {}.root {} > /dev/null'.format(target+suf, ' '.join(sum(fins, []))))
        if clean_unmerge:
            os.system('rm {} '.format( ' '.join(sum(fins, []))))
    os.chdir(pwd)


if __name__ == '__main__':

    print( 'Test!' )
    parser = argparse.ArgumentParser('Run the simple HZZlooper')
    parser.add_argument('outdir', nargs='?', default='temp')
    parser.add_argument('-v', '--verbose', default=1)
    parser.add_argument('-j', '--ncore', default=12)
    parser.add_argument('-d', '--indir', default=None)
    parser.add_argument('-sl', '--samp', default=None)
    parser.add_argument('-mo', '--merge_only', default=False)
    parser.add_argument('-dr', '--dryrun', action='store_true', default=False)
    parser.add_argument('-fg', '--nolog', action='store_true', default=False)
    parser.add_argument('--syst', default='nominal')
    parser.add_argument('--merge', action='store_true', default=False)
    parser.add_argument('--nice', default=10)
    parser.add_argument('--test', action='store_true', default=False)

    args = parser.parse_args()

    os.system('mkdir -p output')
    outdir = 'output/'+args.outdir
    dryrun = args.dryrun
    verbose = args.verbose
    poolsize = args.ncore
    os.nice(args.nice)

    ## Define amples to run over

    with open('sample_lists.yml', 'r') as flst:
        samplists = yaml.safe_load(flst)

    with open('sample_filelists.yml', 'r') as flst:
        filelists = yaml.safe_load(flst)

    if args.merge_only != False:
        lst_samp = args.merge_only.split(',') if len(args.merge_only) > 0 else None
        if len(lst_samp) == 1 and lst_samp[0] == 'all': lst_samp = None
        mergeOutputHists(outdir, samplists['merge_map'], lst_samp=lst_samp)
        exit(0)

    samp_flists_16 = {}
    samp_flists_17 = {}
    samp_flists_18 = {}
    samp_flists_16.update(samplists['data_2016'])
    samp_flists_17.update(samplists['data_2017'])
    samp_flists_18.update(samplists['data_2018'])
    samp_flists_16.update(samplists['bkg_3lskim_16'])
    samp_flists_17.update(samplists['bkg_3lskim_17'])
    samp_flists_18.update(samplists['bkg_3lskim_18'])

    yrsamp_flists = {
        "2016:3LEvents:3lCR" : samp_flists_16,
        "2017:3LEvents:3lCR" : samp_flists_17,
        "2018:3LEvents:3lCR" : samp_flists_18,
    }

    if args.samp != None: samp_flist = samplists[samp]

    rv = os.system('make -j 12')
    if rv != 0: exit()          # quit if make is not successful

    pool = ThreadPool(poolsize)
    njobs_total = 0
    njobs_done = 0

    arglist = []

    systypes = ["Nominal", ]
    if args.syst == 'all':
        systypes = ["Nominal", "JECUp", "JECDn", "METUp", "METDn", "JERUp", "JERDn"]
    elif args.syst == 'nomonly':
        systypes = ["NominalOnly", ]

    # The samples that can be separated
    for yrcfg, dslists in yrsamp_flists.items():
        year = yrcfg.split(':')[0]
        skimtype = yrcfg.split(':')[1]
        outtag = yrcfg.split(':')[2]
        skimdir = osoutdir+'/'+skimtype+skimdate
        tag = '{}_{}_{}'.format(skimver, outtag, year)
        for dsname, spec in dslists.items():
            indir = skimdir+mcsuf+year
            if dsname.startswith('Run') or dsname.startswith('EGamma_'):
                indir = skimdir+datasuf+dsname.replace('Run', '').replace('EGamma_', 'DefaultLeptons/')

            logdir = 'logs/logs_{}_{}'.format(outdir, year[2:])
            os.system('mkdir -p {}'.format(logdir))
            os.system('cp ScanChain.cc {}'.format(logdir)) # back up the state of the looper

            flist = ""
            if spec == None: spec = ""
            elif type(spec) == list:
                flist = ','.join(spec)
            for systype in systypes:
                argstr = '{0} {1} {2} {3} {4} {5}'.format(indir, dsname, tag, systype, spec, flist)
                logname = '{}/log_{}_{}.txt'.format(logdir, dsname, systype)
                if not args.nolog:
                    argstr += ' >& {}'.format(logname)
                arglist.append(argstr)
                njobs_total += 1
                if verbose >= 1:
                    print( '>>> Running over sample {} with log at {}'.format(dsname, logname) )

    ijob = 0
    njobs_done = 0
    for rs in pool.imap(runOSlooper, arglist):
        if rs == 0:
            njobs_done += 1
        else:
            print( '>>> The {}-th job return with {}! See {}'.format(ijob, rs, arglist[ijob].split()[-1]) )
        if verbose >= 1:
            print( '>>> Finish running {}/{} jobs!'.format(njobs_done, njobs_total) )
        ijob += 1

    if njobs_done >= njobs_total-1:
        print( 'All {} job done!'.format(njobs_done) )
        if args.merge:
            mergeOutputHists(outdir, samplists['merge_map'])
    else:
        print( '{} jobs finished out of {}!'.format(njobs_done, njobs_total) )

    print( 'Finished!' )
