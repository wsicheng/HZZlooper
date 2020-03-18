#!/bin/bash

# This is a baby making condor executable for CondorTask of ProjectMetis. Passed in arguments are:
# arguments = [outdir, outname_noext, inputs_commasep, index, cmssw_ver, scramarch, self.arguments]

OUTPUTDIR=$1
OUTPUTNAME=$2
INPUTFILENAMES=$3
IFILE=$4
CMSSW_VERSION=$5
SCRAM_ARCH=$6
NEVENTS=$7
EXTRARG=$8

OUTPUTNAME=$(echo $OUTPUTNAME | sed 's/\.root//')

echo -e "\n--- begin header output ---\n" #                     <----- section division
echo "OUTPUTDIR: $OUTPUTDIR"
echo "OUTPUTNAME: $OUTPUTNAME"
echo "INPUTFILENAMES: $INPUTFILENAMES"
echo "IFILE: $IFILE"
echo "CMSSW_VERSION: $CMSSW_VERSION"
echo "SCRAM_ARCH: $SCRAM_ARCH"

echo "hostname: $(hostname)"
echo "uname -a: $(uname -a)"
echo "system: $(cat /etc/system-release)"
echo "time: $(date +%s)"
echo "args: $@"

echo -e "\n--- end header output ---\n" #                       <----- section division
ls -ltrha
echo ----------------------------------------------

# Unpack the passed in tarfile
tar -xzf package.tar.gz
cd input
ls -ltrha
echo ----------------------------------------------

# Setup Enviroment
export SCRAM_ARCH=$SCRAM_ARCH
source /cvmfs/cms.cern.ch/cmsset_default.sh
pushd /cvmfs/cms.cern.ch/$SCRAM_ARCH/cms/cmssw/$CMSSW_VERSION/src/ > /dev/null
eval `scramv1 runtime -sh`
popd > /dev/null

# The output name is the sample name for stop baby

# if [ ! -z ${EXE_ARGS[1]} ]; then
#   NEVENTTOTAL=${EXE_ARGS[1]}
# fi

echo "Running BabyMaker:"
echo "  ./runHZZlooper $INPUTFILENAMES ${OUTPUTNAME}_${IFILE} ./ $NEVENTS $EXTRARG"
./runHZZlooper $INPUTFILENAMES ${OUTPUTNAME}_${IFILE} ./ $NEVENTS $EXTRARG

echo ----------------------------------------------
ls -ltrha
echo ----------------------------------------------

echo -e "\n--- end running ---\n" #                             <----- section division

# Copy back the output file

if [[ $(hostname) == "uaf"* ]]; then
    mkdir -p ${OUTPUTDIR}
    echo cp ${OUTPUTNAME}_${IFILE}.root ${OUTPUTDIR}/${OUTPUTNAME}_${IFILE}.root
    cp ${OUTPUTNAME}_${IFILE}.root ${OUTPUTDIR}/${OUTPUTNAME}_${IFILE}.root
else
    export LD_PRELOAD=/usr/lib64/gfal2-plugins//libgfal_plugin_xrootd.so # needed in cmssw versions later than 9_3_X
    gfal-copy -p -f -t 4200 --verbose file://`pwd`/${OUTPUTNAME}_${IFILE}.root gsiftp://gftp.t2.ucsd.edu${OUTPUTDIR}/${OUTPUTNAME}_${IFILE}.root --checksum ADLER32
fi

echo -e "\n--- cleaning up ---\n" #                             <----- section division
cd ..
rm -r package.tar.gz input/
