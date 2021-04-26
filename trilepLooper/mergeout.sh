
skimtype=$1
dir=$2

if [ $skimtype == "3lskim" ]; then
    cd $dir;
    rename _Nominal.root .root *
    hadd -f data_3lskim.root Run201??.root  &> /dev/null
    hadd -f DY.root DYJetsToLL_M-50.root &> /dev/null
    hadd -f WZ.root WZTo3LNu_mllmin01.root WZTo2*.root &> /dev/null
    hadd -f WW.root WWTo2L2Nu.root &> /dev/null
    hadd -f ttW.root TTWJetsToLNu.root &> /dev/null
    hadd -f triboson.root W??.root Z??.root  &> /dev/null
    hadd -f ZZ.root ZZTo*.root &> /dev/null
    hadd -f ZH.root ZH_*.root &> /dev/null
    [[ -f ZGTo2LG.root ]] && hadd -f ZG.root ZGTo2LG.root &> /dev/null
    [[ -f ZLLGJets.root ]] && hadd -f ZG.root ZLLGJets.root &> /dev/null
    hadd -f top.root TTTo2L2Nu.root TTGJets.root &> /dev/null
    hadd -f ttZ.root TTZTo*.root &> /dev/null
    hadd -f Others.root WW.root TTGJets.root WGToLNuG.root ZZ.root ZH.root ttZ.root tZq.root  &> /dev/null
    hadd -f allBkgs_3l.root ZG.root DY.root Others.root &> /dev/null
    cd -;
fi

if [ $skimtype == "3lrun2" ]; then
    mkdir -p ${dir}_run2
    declare -a Samples=( DY ZGTo2LG WZ WW ttW ttZ top triboson ZZ ZG  Others )
    hadd -f ${dir}_run2/data_3lskim.root ${dir}_201?/data_3lskim.root &> /dev/null
    for sample in ${Samples[@]}; do
        hadd -f ${dir}_run2/$sample.root ${dir}_201?/$sample.root &> /dev/null
    done
fi
