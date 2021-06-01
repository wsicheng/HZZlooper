
skimtype=$1
dir=$2

if [ $skimtype == "3lskim" ]; then
    cd $dir;
    rename _Nominal.root .root *
    hadd -f data_3lskim.root Run201??.root  &> /dev/null
    hadd -f DY.root DYJetsToLL_M-50.root &> /dev/null
    hadd -f WZ.root WZTo3LNu_mllmin01.root WZTo2*.root &> /dev/null
    hadd -f WW.root WWTo2L2Nu.root &> /dev/null
    hadd -f Triboson.root W??.root Z??.root  &> /dev/null
    hadd -f ZZ.root ZZTo*.root &> /dev/null
    hadd -f ZH.root ZH_*.root &> /dev/null
    [[ -f ZGTo2LG.root ]] && hadd -f ZG.root ZGTo2LG.root &> /dev/null
    [[ -f ZLLGJets.root ]] && hadd -f ZG.root ZLLGJets.root &> /dev/null
    # hadd -f ttZ.root TTZTo*.root &> /dev/null
    # hadd -f top.root TTZTo*.root TTTo2L2Nu.root TTGJets.root &> /dev/null
    rename TTTo2L2Nu TTbar *
    hadd -f tVX.root TTZTo*.root TTGJets.root TTWJetsToLNu.root tZq.root &> /dev/null
    # hadd -f Others.root WW.root TTGJets.root WGToLNuG.root ZZ.root ZH.root ttZ.root tZq.root  &> /dev/null
    # hadd -f Others.root WW.root TTGJets.root WGToLNuG.root ZZ.root ttZ.root Triboson.root tZq.root  &> /dev/null
    hadd -f Others.root WW.root WGToLNuG.root Triboson.root  &> /dev/null
    # hadd -f allBkgs_3l.root ZG.root DY.root Others.root &> /dev/null
    cd -;
fi

if [ $skimtype == "3lrun2" ]; then
    mkdir -p ${dir}_run2
    declare -a Samples=( DY ZGTo2LG WZ WW TTbar tVX Triboson ZZ ZG ZH Others )
    hadd -f ${dir}_run2/data_3lskim.root ${dir}_201?/data_3lskim.root &> /dev/null
    for sample in ${Samples[@]}; do
        hadd -f ${dir}_run2/$sample.root ${dir}_201?/$sample.root &> /dev/null
    done
fi

if [ $skimtype == "3ltree" ]; then
    cd $dir;
    hadd -f trileptree_data_Nominal.root trileptree_Run201??_Nominal.root &> /dev/null
    # rm finaltree_Run201??_*.root
    rename _Tune _ALT *WZTo3LNu_Tune*.root
    rename ZLLGJets ZGJets_ll *ZLLGJets*.root
    rename ZGTo2LG ZGJets_ll *ZGTo2LG*.root
    # rm trileptree_ZGJets_ll_PtG-130_*.root
    # for s in *WGToLNuG_01J_*; do hadd -f finaltree_qqWG_lnu${s:23} ${s/_01J_/_} $s &> /dev/null && rm ${s/_01J_/_} $s; done
    for s in *WGToLNuG_01J_*; do mv $s trileptree_qqWG_lnu${s:23} && rm ${s/_01J_/_}; done
    rename WGToLNuG qqWG_lnu *.root
    cd -;
fi
