
skimtype=$1
dir=$2

if [ $skimtype == "llskim" ]; then
    cd $dir;
    rename _Nominal.root .root *
    [ -f WGToLNuG_01J.root ] && mv WGToLNuG_01J.root WGToLNuG.root
    hadd -f triboson.root WWW.root WWZ.root ZZZ.root WZG.root WZZ.root &> /dev/null;
    # rm WWW.root WWZ.root ZZZ.root WZG.root WZZ.root;
    # hadd -f ZZ.root ZZTo2L2Nu.root &> /dev/null
    # hadd -f diboson.root WWToLNuQQ.root WZTo1L3Nu.root WZTo1L1Nu2Q.root WZTo3LNu.root WZTo2L2Q.root ZLLGJets.root ZZTo2L2Q.root ZZTo2Q2Nu.root ZZTo4L.root WGToLNuG.root &> /dev/null
    # rm WWToLNuQQ.root WZTo1L3Nu.root WZTo1L1Nu2Q.root WZTo3LNu.root WZTo2L2Q.root ZLLGJets.root ZZTo2L2Q.root ZZTo2Q2Nu.root ZZTo4L.root WGToLNuG.root
    # hadd -f Others.root TTGJets.root TGJets.root QCD*.root tZq_ll.root ST.root TTWJetsToLNu.root WJetsToLNu.root triboson.root &> /dev/null
    # rm TTGJets.root TGJets.root QCD_HT.root tZq_ll.root ST.root TTWJetsToLNu.root WJetsToLNu.root
    hadd -f nonres.root TTTo2L2Nu.root TTZToLLNuNu.root WWTo2L2Nu.root &> /dev/null
    hadd -f DY.root DYJetsToLL_*.root &> /dev/null
    hadd -f ZZ.root ZZTo2L2Nu.root ZZTo2L2Q.root ZZTo2Q2Nu.root ZZTo4L.root &> /dev/null
    hadd -f WZ.root WZTo1L3Nu.root WZTo1L1Nu2Q.root WZTo3LNu.root WZTo2L2Q.root &> /dev/null
    hadd -f Others.root WWToLNuQQ.root TTGJets.root TGJets.root QCD*.root tZq_ll.root ST.root TTWJetsToLNu.root WJetsToLNu.root WGToLNuG.root ZLLGJets.root triboson.root &> /dev/null
    # hadd -f diboson.root WWToLNuQQ.root WZTo1L3Nu.root WZTo1L1Nu2Q.root WZTo3LNu.root WZTo2L2Q.root   &> /dev/null
    cd -;
fi

if [ $skimtype == "phskim" ]; then
    cd $dir;
    rename _Nominal.root .root *
    hadd -f data_phskim.root Run201??.root  &> /dev/null
    [ -f QCD_HT.root ] && mv QCD_HT.root QCD.root
    [[ -f WGToLNuG_01J.root ]] && [[ -f WGToLNuG.root ]] && hadd -f qqWG.root WGToLNuG*.root &> /dev/null
    # [ -f WGToLNuG_01J.root ] && mv WGToLNuG_01J.root WGToLNuG.root
    hadd -f Others.root TGJets.root TTGJets.root TTJets.root WZG.root ZJetsToNuNu.root  &> /dev/null
    # hadd -f subtractor.root ZGTo2NuG.root WlnuFromCR.root WGToLNuG.root Others.root  &> /dev/null
    hadd -f subtractor.root ZGTo2NuG.root WJetsToLNu.root qqWG.root Others.root  &> /dev/null
    hadd -f allbkg.root GJets.root subtractor.root &> /dev/null
    cd -;
fi

if [ $skimtype == "llgskim" ]; then
    cd $dir;
    rename _Nominal.root .root *
    hadd -f data_llgskim.root Run201??.root  &> /dev/null
    [ -f DYJetsToLL_*.root ] && hadd -f DY.root DYJetsToLL_*.root &> /dev/null
    hadd -f Others.root WZ*.root ZZ*.root  &> /dev/null
    [ -f ZLLGJets.root ] && hadd -f ZGTo2LG.root ZLLGJets.root ZGTo2LG_PtG-130.root &> /dev/null
    hadd -f allBkgs_llg.root ZGTo2LG.root DY.root Others.root &> /dev/null
    cd -;
fi

if [ $skimtype == "slskim" ]; then
    phdir=${dir/_slCR/_phCR}
    echo $phdir
    cd $dir;
    rename _Nominal.root .root *
    hadd -f data_slskim.root EGamma*.root  &> /dev/null
    cd -;
    cd $phdir
    ln -s ${dir/output/..}/data_slskim.root WlnuFromCR.root
    cd -;
fi

if [ $skimtype == "phtree" ]; then
    cd $dir;
    declare -a Systs=(Nominal TransferFactorUp TransferFactorDn )
    for syst in ${Systs[@]}; do
        hadd -f finaltree_data_$syst.root finaltree_Run201??_$syst.root &> /dev/null
    done
    # rm finaltree_Run201??_*.root
    rename ZGTo2NuG ZGJets_nunu *.root
    # rename WGToLNuG_01J qqWG_lnu_ALT *.root
    for s in *WGToLNuG_01J_*; do hadd -f finaltree_qqWG_lnu${s:22} ${s/_01J_/_} $s &> /dev/null && rm ${s/_01J_/_} $s; done
    rename WGToLNuG qqWG_lnu *.root
    for i in finaltree_GJets_*.root; do [[ -f ${i/GJets_/GJets_HT_} ]] && rm ${i/GJets_/GJets_HT_}; done
    rename GJets GJets_HT finaltree_GJets*.root
    rename lnu_01J lnu_ALT *.root
    rename ZJetsToNuNu ZJets_nunu_HT *.root
    rename _HT_HT _HT *.root
    # rename lnu_01J lnu *.root
    cd -;
fi

if [ $skimtype == "sltree" ]; then
    cd $dir;
    hadd -f finaltree_WJetsToLNu_datadriven_Nominal.root finaltree_EGamma_*_Nominal.root &> /dev/null
    declare -a Systs=(Nominal TransferFactorUp TransferFactorDn ElecToPhotonScaleUp ElecToPhotonScaleDn PhoTriggerEffUp PhoTriggerEffDn ElecTriggerEffUp ElecTriggerEffDn )
    for syst in ${Systs[@]}; do
        hadd -f finaltree_WJetsToLNu_datadriven_$syst.root finaltree_EGamma_*_$syst.root &> /dev/null
    done
    # hadd -f finaltree_WJetsToLNu_datadriven_ElecToPhotonScaleUp.root finaltree_EGamma_*_ElecToPhotonScaleUp.root &> /dev/null
    # hadd -f finaltree_WJetsToLNu_datadriven_ElecToPhotonScaleDn.root finaltree_EGamma_*_ElecToPhotonScaleDn.root &> /dev/null
    cd -;
fi

if [ $skimtype == "fthist" ]; then
    cd $dir;
    rename _M-50 "" *.root
    hadd -f Others.root TGJets.root TTGJets.root TTJets.root WZG.root  &> /dev/null
    hadd -f subtractor.root ZGJets_nunu.root WJetsToLNu.root qqWG_lnu.root Others.root  &> /dev/null
    hadd -f allbkg.root GJets_HT.root subtractor.root
    cd -;
fi

if [ $skimtype == "ftrun2" ]; then
    mkdir -p ${dir}_run2
    declare -a Samples=( data GJets_HT ZGJets_nunu qqWG_lnu WJetsToLNu WZG TGJets TTGJets TTJets ZZTo2L2Nu DYJetsToLL )
    for sample in ${Samples[@]}; do
        hadd -f ${dir}_run2/$sample.root ${dir}_201?/$sample*.root &> /dev/null
    done
    # rename _M-50 "" *.root
    cd ${dir}_run2
    hadd -f Others.root TGJets.root TTGJets.root TTJets.root WZG.root  &> /dev/null
    hadd -f subtractor.root ZGJets_nunu.root WJetsToLNu.root qqWG_lnu.root Others.root  &> /dev/null
    hadd -f allbkg.root GJets_HT.root subtractor.root
    cd -;
fi

if [ $skimtype == "llgrun2" ]; then
    mkdir -p ${dir}_run2
    declare -a Samples=( DY ZGTo2LG WZTo2L2Q WZTo3LNu ZZTo2L2Nu ZZTo2L2Q Others )
    hadd -f ${dir}_run2/data_llgskim.root ${dir}_201?/data_llgskim.root &> /dev/null
    for sample in ${Samples[@]}; do
        hadd -f ${dir}_run2/$sample.root ${dir}_201?/$sample.root &> /dev/null
    done
fi

if [ $skimtype == "phrun2" ]; then
    mkdir -p ${dir}_run2
    hadd -f ${dir}_run2/data_phskim.root ${dir}_201?/data_phskim.root &> /dev/null
    # declare -a Samples=( GJets ZGTo2NuG WGToLNuG ZJetsToNuNu TGJets TTGJets TTJets WZG QCD_HT WJetsToLNu WlnuFromCR Others)
    declare -a Samples=( GJets ZGTo2NuG qqWG ZJetsToNuNu TGJets TTGJets TTJets WZG QCD_HT WJetsToLNu WlnuFromCR Others)
    for sample in ${Samples[@]}; do
        hadd -f ${dir}_run2/$sample.root ${dir}_201?/$sample*.root &> /dev/null
    done
fi

if [ $skimtype == "llrun2" ]; then
    mkdir -p ${dir}_run2
    declare -a Samples=( data GJets_HT ZGJets_nunu qqWG_lnu WJetsToLNu WZG TGJets TTGJets TTJets ZZTo2L2Nu DYJetsToLL )
    for sample in ${Samples[@]}; do
        hadd -f ${dir}_run2/$sample.root ${dir}_201?/$sample*.root &> /dev/null
    done
    cd $dir;
    rename _Nominal.root .root *
    [ -f WGToLNuG_01J.root ] && mv WGToLNuG_01J.root WGToLNuG.root
    hadd -f triboson.root WWW.root WWZ.root ZZZ.root WZG.root WZZ.root &> /dev/null;
    # rm WWW.root WWZ.root ZZZ.root WZG.root WZZ.root;
    # hadd -f ZZ.root ZZTo2L2Nu.root &> /dev/null
    # hadd -f diboson.root WWToLNuQQ.root WZTo1L3Nu.root WZTo1L1Nu2Q.root WZTo3LNu.root WZTo2L2Q.root ZLLGJets.root ZZTo2L2Q.root ZZTo2Q2Nu.root ZZTo4L.root WGToLNuG.root &> /dev/null
    # rm WWToLNuQQ.root WZTo1L3Nu.root WZTo1L1Nu2Q.root WZTo3LNu.root WZTo2L2Q.root ZLLGJets.root ZZTo2L2Q.root ZZTo2Q2Nu.root ZZTo4L.root WGToLNuG.root
    # hadd -f Others.root TTGJets.root TGJets.root QCD*.root tZq_ll.root ST.root TTWJetsToLNu.root WJetsToLNu.root triboson.root &> /dev/null
    # rm TTGJets.root TGJets.root QCD_HT.root tZq_ll.root ST.root TTWJetsToLNu.root WJetsToLNu.root
    hadd -f nonres.root TTTo2L2Nu.root TTZToLLNuNu.root WWTo2L2Nu.root &> /dev/null
    hadd -f DY.root DYJetsToLL_*.root &> /dev/null
    hadd -f ZZ.root ZZTo2L2Nu.root ZZTo2L2Q.root ZZTo2Q2Nu.root ZZTo4L.root &> /dev/null
    hadd -f WZ.root WZTo1L3Nu.root WZTo1L1Nu2Q.root WZTo3LNu.root WZTo2L2Q.root &> /dev/null
    hadd -f Others.root WWToLNuQQ.root TTGJets.root TGJets.root QCD*.root tZq_ll.root ST.root TTWJetsToLNu.root WJetsToLNu.root WGToLNuG.root ZLLGJets.root triboson.root &> /dev/null
    # hadd -f diboson.root WWToLNuQQ.root WZTo1L3Nu.root WZTo1L1Nu2Q.root WZTo3LNu.root WZTo2L2Q.root   &> /dev/null
    cd -;
fi
