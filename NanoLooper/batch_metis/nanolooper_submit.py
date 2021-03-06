from __future__ import print_function

import os
import time
import itertools
import json
import traceback

from metis.Sample import DBSSample, DirectorySample
from metis.CMSSWTask import CMSSWTask, CondorTask
from metis.StatsParser import StatsParser
from metis.Utils import send_email,do_cmd


if __name__ == "__main__":


    dsinfos = {
        # diBoson
        # '/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM|3.05|1.00262|1',
        # 'DYJetsToLL_M-50' : '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20_ext2-v1/NANOAODSIM|130939668',

        'ZGTo2NuG' : '/ZGTo2NuG_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
        'ZGTo2NuG_PtG-130' : '/ZGTo2NuG_PtG-130_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
        'ZNuNuGJets_MonoPhoton_PtG-40to130' : '/ZNuNuGJets_MonoPhoton_PtG-40to130_TuneCP5_13TeV-madgraph/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
        'ZNuNuGJets_MonoPhoton_PtG-130' : '/ZNuNuGJets_MonoPhoton_PtG-130_TuneCP5_13TeV-madgraph/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
        'ZGTo2LG_PtG-130' : '/ZGTo2LG_PtG-130_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',

        'GJets_HT-40To100' : '/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM|9371355',
        'GJets_HT-100To200' : '/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_4cores5k_102X_upgrade2018_realistic_v20-v1/NANOAODSIM|9798176',
        'GJets_HT-200To400' : '/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM|19062809',
        'GJets_HT-400To600' : '/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM|4655985',
        'GJets_HT-600ToInf' : '/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20_ext1-v1/NANOAODSIM|4981121',

        'EGamma_2018A' : '/EGamma/Run2018A-Nano25Oct2019-v1/NANOAOD',
        'EGamma_2018B' : '/EGamma/Run2018B-Nano25Oct2019-v1/NANOAOD',
        'EGamma_2018C' : '/EGamma/Run2018C-Nano25Oct2019-v1/NANOAOD',
        'EGamma_2018D' : '/EGamma/Run2018D-Nano25Oct2019-v1/NANOAOD',

        # 'data_2018A_singlemu' : '/SingleMuon/Run2018A-Nano25Oct2019-v1/NANOAOD',
    }

    publish_to_dis = False

    ntuple_tasks = []
    for outname, dsinfo in dsinfos.items():
        infos = dsinfo.split("|");
        if len(infos) > 0: dsname = infos[0].strip()
        if len(infos) > 1: nevts = float(infos[1].strip())
        else: nevts = -1

        cmsswver = "CMSSW_10_1_0"
        scramarch = "slc6_amd64_gcc700"
        tag = "v1_6"

        tarfile = "tarfiles/input_"+tag+".tar.gz"

        sample = DBSSample( dataset=dsname,)
        task = CondorTask(
            sample = sample,
            # events_per_output = 450e3,
            files_per_output = 50 if 'data' in outname else 5,
            outdir_name = "HZZlooperOutput",
            output_name = outname+".root",
            tag = tag,
            scram_arch = scramarch,
            # global_tag = "", # if global tag blank, one from DBS is used
            executable = "condor_executable.sh",
            arguments = "{} ".format(sample.get_nevents() if nevts < 0 else nevts),
            cmssw_version = cmsswver,
            tarfile = tarfile,
            # condor_submit_params = {"sites": "UCSB"},
            # condor_submit_params = {"sites": "UAF,T2_US_UCSD"},
            condor_submit_params = {"use_xrootd": True},
            # publish_to_dis = publish_to_dis,
            no_load_from_backup = True,
            # recopy_inputs = True,
        )
        ntuple_tasks.append(task)

    for i in range(100):

        total_summary = {}
        # for dsname in dataset_names:
        for task in ntuple_tasks:
            
            try:
                task.process()
            except:
                traceback_string = traceback.format_exc()
                print( "Runtime error:\n{0}".format(traceback_string) )
                # send_email(subject="metis error", body=traceback_string)

            if task.complete():
                hdpdir = task.get_outputdir()
                outdir = '../output/{}_{}'.format(task.outdir_name, task.tag)
                cmdstr = '[[ -d {1} ]] || cp -r {0} {1}'.format( hdpdir, outdir )
                print( cmdstr )
                os.system(cmdstr)
            # if not task.complete():
            #     f = open("unfinished.txt","a+")
            #     print ""
            #     for output in task.get_uncompleted_outputs():
            #         f.write("for ouput file <"+output.get_name()+">:\n")
            #         for inp in task.get_inputs_for_output(output):
            #             f.write(inp.get_name()+'\n')
            #     f.write("\n\n")
            #     f.close()

            total_summary[task.get_sample().get_datasetname()] = task.get_task_summary()

        StatsParser(data=total_summary, webdir="~/public_html/dump/metis_HZZlooper/").do()

        time.sleep(60.*30)

