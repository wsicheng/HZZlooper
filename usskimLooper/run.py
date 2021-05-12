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
skimdir_2l = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Skims/DileptonSkims/2018/200326/"
skimdir_ph = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Skims/GammaJetsSkims/2018/200326"

prodbase = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Production/200906_{year}"
workbase = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/Skims/200906_{year}/"
skimbase = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Skims/200906_{year}"
skimbd17 = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Skims/201221_2017/"

effbase = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/SinglePhotonTriggerEfficiencies/"

# skimver  = "v4_04";
# skimdate = "/SkimTrees/201128/AK4Jets_WithPUJetId_NoTightLeptonJetId_ParticleCleaned";
# mcsuf    = "/PFMET_WithXY_WithJER_WithPartMomCorr_P4Preserved_ResolutionCorrected/";

# skimver = "v4_05";
# skimver = "v4_06";
# skimver = "v4_07";
skimver  = "v4_08";
skimdate = "/SkimTrees/210315/AK4Jets_WithPUJetId_NoTightLeptonJetId_ParticleCleaned";
skimbase = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Skims/210313_{year}"

skimver  = "v5_00";
skimdate = "/SkimTrees/210504/AK4Jets_WithPUJetId_NoTightLeptonJetId_ParticleCleaned";

osoutdir = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output";
datasuf  = "/PFMET_WithXY_WithPartMomCorr_P4Preserved/";
mcsuf    = "/PFMET_WithXY_NoJER_WithPartMomCorr_P4Preserved_ResolutionCorrected/";

# skimver  = "v4_09";
# skimdate = "/SkimTrees/210329/AK4Jets_WithPUJetId_NoTightLeptonJetId_ParticleCleaned";
# mcsuf    = "/PFMET_WithXY_NoJER_WithPartMomCorr_P4Preserved_ResolutionUncorrected/";

def runOSlooper(args):
    if verbose >= 3:
        print( '  ./runOSlooper '+args)
    if not dryrun:
        return os.system('./runOSlooper '+args)
    return 0

def mergeOutputHists(outdir, merge_map, suf='', lst_samp=None, clean_unmerge=False):
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
    # parser.add_argument('--nomerge', action='store_true', default=False)
    parser.add_argument('--merge', action='store_true', default=False)
    parser.add_argument('--nice', default=10)
    parser.add_argument('--test', action='store_true', default=False)

    args = parser.parse_args()

    outdir = 'output/'+args.outdir
    dryrun = args.dryrun
    verbose = args.verbose
    poolsize = args.ncore
    os.nice(args.nice)

    ## Define amples to run over

    with open('sample_lists.yml', 'r') as flst:
        samplists = yaml.safe_load(flst)

    if args.merge_only != False:
        lst_samp = args.merge_only.split(',') if len(args.merge_only) > 0 else None
        if len(lst_samp) == 1 and lst_samp[0] == 'all': lst_samp = None
        mergeOutputHists(outdir, samplists['merge_map'], lst_samp=lst_samp)
        exit(0)

    # data_flists = data17_flists
    data_flists = {}
    samp_flists = {}
    samp_flists_16 = {}
    samp_flists_17 = {}
    samp_flists_18 = {}
    samp_flists_16.update(samplists['data_2016'])
    samp_flists_17.update(samplists['data_2017'])
    samp_flists_18.update(samplists['data_2018'])
    samp_flists_16.update(samplists['bkg_phskim_16'])
    samp_flists_17.update(samplists['bkg_phskim_17'])
    samp_flists_18.update(samplists['bkg_phskim_18'])
    # samp_flists_16.update(samplists['bkg_llgskim_16'])
    # samp_flists_17.update(samplists['bkg_llgskim_17'])
    # samp_flists_18.update(samplists['bkg_llgskim_18'])
    # samp_flists_16.update({'ZGTo2NuG' : '',})
    # samp_flists_17.update({'ZGTo2NuG' : '',})
    # samp_flists_18.update({'ZGTo2NuG' : '',})

    # samp_flists_17.update(samplists['bkg_wgtolnu_1718'])
    # samp_flists_18.update(samplists['bkg_wgtolnu_1718'])

    indir = skimbase

    yrsamp_flists = {
        "2016:SinglePhotonEvents:phCR_all2jsel_rwgtd" : samp_flists_16,
        "2017:SinglePhotonEvents:phCR_all2jsel_rwgtd" : samp_flists_17,
        "2018:SinglePhotonEvents:phCR_all2jsel_rwgtd" : samp_flists_18,
        # "2016:SingleLeptonEvents:slCR_all2jsel_rwgtd" : samplists['egamma_2016'],
        # "2017:SingleLeptonEvents:slCR_all2jsel_rwgtd" : samplists['egamma_2017'],
        # "2018:SingleLeptonEvents:slCR_all2jsel_rwgtd" : samplists['egamma_2018'],
        # "2016:LLGEvents:llgCR" : samp_flists_16,
        # "2017:LLGEvents:llgCR" : samp_flists_17,
        # "2018:LLGEvents:llgCR" : samp_flists_18,
        # "2016:DileptonEvents:2lSR_all2jsel" : samplists['data_2016'],
        # "2017:DileptonEvents:2lSR_all2jsel" : samplists['data_2017'],
        # "2018:DileptonEvents:2lSR_all2jsel" : samplists['data_2018'],
        # "2016:SinglePhotonEvents:phCR_all2jsel_rwgtd_raweta" : samplists['data_2016'],
        # "2017:SinglePhotonEvents:phCR_all2jsel_rwgtd_raweta" : samplists['data_2017'],
        # "2018:SinglePhotonEvents:phCR_all2jsel_rwgtd_raweta" : samplists['data_2018'],
        # "2016:LLGEvents:temp" : samplists['mc_tester'],
        # "2016:SinglePhotonEvents:temp" : samplists['data_tester'],
        # "2018:SinglePhotonEvents:phCR_all2jsel_rwgtd_wgtsbkup" : samp_flists_18,
        # "2018:SinglePhotonEvents:phCR_ee2j_test" : samplists['data_tester'],
        # "2018:SinglePhotonEvents:phCR_all2jsel_rwgtd_tmp" : samplists['mc_tester'],
    }

    if args.indir != None: indir = args.indir
    if args.samp != None: samp_flist = samplists[samp]

    if not indir.endswith('/'): indir += '/'

    rv = os.system('make -j 12')
    if rv != 0: exit()          # quit if make is not successful

    pool = ThreadPool(poolsize)
    njobs_total = 0
    njobs_done = 0

    arglist = []

    systypes = ["Nominal", "JECUp", "JECDn", "METUp", "METDn", "JERUp", "JERDn"]
    # systypes = ["Nominal", ]

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
            os.system('cp ScanChain.cc {}'.format(logdir))

            if '{year}' in indir:
                idir = indir.format(year=year)
                idir2 = workbase.format(year=year) # backup indir
            if True:
                if spec == None: spec = ""
                for systype in systypes:
                    argstr = '{0} {1} {2} {3} {4}'.format(indir, dsname, tag, systype, spec)
                    logname = '{}/log_{}_{}.txt'.format(logdir, dsname, systype)
                    if not args.nolog: argstr += ' >& {}'.format(logname)
                    arglist.append(argstr)
                    njobs_total += 1
                    if verbose >= 1:
                        print( '>>> Running over sample {} with log at {}'.format(dsname, logname) )
            else:
                flist = glob.glob('{}/{}*.root'.format(indir, dsname)) 

                nfile_per_job = 10
                # if sumwgt < 0: nfile_per_job = 1000
                njobs = int(math.ceil(len(flist) / nfile_per_job))
                njobs_total += njobs
                for i in range(njobs):
                    infilestr = ','.join(flist[ i*nfile_per_job : (i+1)*nfile_per_job ])
                    outname = dsname if njobs == 1 else '{}_{}'.format(dsname, i)
                    argstr = '{0} {1} {2} {3} {4}'.format(indir, infilestr, tag, systype, spec)
                    if not args.nolog: argstr += ' >& {}/log_{}.txt'.format(logdir, outname)
                    arglist.append(argstr)
                    # print( argstr )
                if verbose >= 1:
                    print( '>>> Running sample {} with {} jobs!'.format(dsname, njobs) )
                
    ijob = 0
    njobs_done = 0
    for rs in pool.imap_unordered(runOSlooper, arglist):
        if rs == 0:
            njobs_done += 1
        else:
            # print( '>>> The {}-th job return with {}! See {}'.format(ijob, rs, arglist[ijob].split()[-1]) )
            print( '>>> The {}-th job return with {}! '.format(ijob, rs) ) # unordered
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
        
