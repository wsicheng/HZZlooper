#! /bin/bash

METIS_PATH=~/working/ProjectMetis
BATCH_DIR=$PWD
INPUT_DIR=..
tag=${1}

# Checkout ProjectMetis, if necessary
if [ -z $METIS_PATH ] && [ ! -d ProjectMetis ]; then
    git clone git@github.com:aminnj/ProjectMetis.git
elif [ ! -d ProjectMetis ]; then
    ln -s $METIS_PATH ProjectMetis
fi

# Setup the environment for ProjectMetis
cd ProjectMetis
. setup.sh
cd ..

# Make the input tarfile

[[ `ls input/` ]] && echo "[setup] Cleaning input/" && rm -r input/*
mkdir -p input
echo "[setup] Making input tarfile... This assumes that the babymaker has been compiled and tested."
pushd $INPUT_DIR > /dev/null
cp -rL *.so config runHZZlooper testrun.sh $BATCH_DIR/input
popd > /dev/null
tar -czf input.tar.gz input
[[ ! -z $tag ]] && echo "[setup] Creating tag $tag" && cp input.tar.gz tarfiles/input_$tag.tar.gz

# if [ ! -f merge_scripts.tar.gz ]; then
#     echo "[setup] Making tarfile for merge scripts."
#     tar -czf merge_scripts.tar.gz mergeHadoopFiles.C ../skimBaby.C
# fi

# Flag for successful environment setup
CONDOR_ENV_READY=true
