#!/bin/bash

RunningJobs=`jobs | grep runHZZlooper`
if [[ $RunningJobs == "" ]]; then
    make -j 12 || return $?
else
    echo "Running looper jobs exist, will not re-build the looper."
    read -p "Do append jobs? (use Ctrl-c to break if not) " dummy
    unset dummy
fi

LOGDIR=logs

# skimtype=merged

OUTDIR18=output/samp18_$vsuf
OUTDIR17=output/samp17_$vsuf
OUTDIR16=output/samp16_$vsuf
OUTDIRR2=output/combRun2_$vsuf

run18dat=0
run18bkg=0
run18sig=1

run17dat=0
run17bkg=0
run17sig=1

run16dat=0
run16bkg=0
run16sig=1

function runLooperJobs {
    # [[ -d ${OUTDIR} ]] && rm -r ${OUTDIR}
    mkdir -p ${OUTDIR}; mkdir -p ${LOGDIR}
    cp do.sh ScanChain.cc HZZSelections.cc ${LOGDIR}
    for SAMPLE in ${Samples[@]}; do
        ./runHZZlooper ${INDIR} ${SAMPLE} ${OUTDIR}
        echo ./runHZZlooper ${INDIR} ${SAMPLE} ${OUTDIR} '>&' ${LOGDIR}/log_${SAMPLE}.txt
        # eval "nohup nice -n 10 ./runHZZlooper ${INDIR} ${SAMPLE} ${OUTDIR} >& ${LOGDIR}/log_${SAMPLE}.txt &"
    done
}

########################
# 2018 Data

INDIR=/pnfs/iihe/cms/store/group/HZZ2l2nu/Production/2019-10-07_2017/stageout/

# Temporary test
INDIR=/pnfs/iihe/cms/store/group/HZZ2l2nu/Production/2019-10-07_2017/stageout/ZZTo4L_13TeV_powheg_pythia8/crab_ZZTo4L_p0_ZZ2l2vPruner-MC/191007_091528/0000/
OUTDIR=output/temp
LOGDIR=$OUTDIR/logs

declare -a Samples=( ZZTo4L_p0_ZZ2l2vPruner-MC_1 )

runLooperJobs
