
{
    gSystem->Exec("mkdir -p plots");

    gROOT->ProcessLine(".L SkimTree.cc+");
    gROOT->ProcessLine(".L ScanChain.C+");

    TString skimver = "v2_07_metcrtd_2018";

    TString tagsuf = "_zgest3";
    TString tag = skimver+tagsuf;

    // gSystem->Exec(Form("mkdir -p output/%s", tag.Data()));

    // ScanChain(skimdir+skimver, "DY", tag);
    // ScanChain(skimdir+skimver, "GJets", tag);
    // ScanChain(skimdir+skimver, "QCD", tag);
    // ScanChain(skimdir+skimver, "ZJets", tag);

    // ScanChain(skimdir+skimver, "qqWG", tag);
    // ScanChain(skimdir+skimver, "WZG", tag);
    // ScanChain(skimdir+skimver, "TGJets", tag);
    // ScanChain(skimdir+skimver, "TTGJets", tag);

    // ScanChain(skimdir+skimver, "TT_2l2nu", tag);
    // ScanChain(skimdir+skimver, "ST", tag);
    // ScanChain(skimdir+skimver, "qqZZ", tag);
    // ScanChain(skimdir+skimver, "qqWZ", tag);

    // ScanChain(skimdir+skimver, "ZGJets_nunu_pTG_40", tag);
    // ScanChain(skimdir+skimver, "ZGJets_nunu_pTG_130", tag);
    // ScanChain(skimdir+skimver, "ZGJets_nunu_nlo_incl", tag);
    // ScanChain(skimdir+skimver, "ZGJets_nunu_nlo_pTG_130", tag);
    // ScanChain(skimdir+skimver, "ZGJets_ll_pTG_40", tag);
    // ScanChain(skimdir+skimver, "ZGJets_ll_pTG_130", tag);
    // ScanChain(skimdir+skimver, "ZGJets_ll_nlo_incl", tag);
    // ScanChain(skimdir+skimver, "ZGJets_ll_nlo_pTG_130", tag);

    // ScanChain(skimdir+skimver, "2018A_ph", tag);
    // ScanChain(skimdir+skimver, "2018B_ph", tag);
    // ScanChain(skimdir+skimver, "2018C_ph", tag);
    // ScanChain(skimdir+skimver, "2018D_ph", tag);

    // ScanChain(skimdir+skimver, "2018A_ll", tag);
    // ScanChain(skimdir+skimver, "2018B_ll", tag);
    // ScanChain(skimdir+skimver, "2018C_ll", tag);
    // ScanChain(skimdir+skimver, "2018D_ll", tag);

    skimver = "v2_07_metcrtd_2017";
    tag = skimver+tagsuf;

    // gSystem->Exec(Form("mkdir -p output/%s", tag.Data()));

    // ScanChain(skimdir+skimver, "DY", tag);
    // ScanChain(skimdir+skimver, "GJets", tag);
    // ScanChain(skimdir+skimver, "QCD", tag);
    // ScanChain(skimdir+skimver, "ZJets", tag);

    // ScanChain(skimdir+skimver, "qqWG", tag);
    // ScanChain(skimdir+skimver, "WZG", tag);
    // ScanChain(skimdir+skimver, "TGJets", tag);
    // ScanChain(skimdir+skimver, "TTGJets", tag);
    // ScanChain(skimdir+skimver, "TT_2l2nu", tag);

    // ScanChain(skimdir+skimver, "ST", tag);
    // ScanChain(skimdir+skimver, "qqZZ", tag);
    // ScanChain(skimdir+skimver, "qqWZ", tag);

    // ScanChain(skimdir+skimver, "ZGJets_nunu_pTG_40", tag);
    // ScanChain(skimdir+skimver, "ZGJets_nunu_pTG_130", tag);
    // ScanChain(skimdir+skimver, "ZGJets_nunu_nlo_incl", tag);
    // ScanChain(skimdir+skimver, "ZGJets_nunu_nlo_pTG_130", tag);
    // ScanChain(skimdir+skimver, "ZGJets_ll_pTG_40", tag);
    // ScanChain(skimdir+skimver, "ZGJets_ll_pTG_130", tag);
    // ScanChain(skimdir+skimver, "ZGJets_ll_nlo_incl", tag);
    // ScanChain(skimdir+skimver, "ZGJets_ll_nlo_pTG_130", tag);

    // ScanChain(skimdir+skimver, "2017B_ph", tag);
    // ScanChain(skimdir+skimver, "2017C_ph", tag);
    // ScanChain(skimdir+skimver, "2017D_ph", tag);
    // ScanChain(skimdir+skimver, "2017E_ph", tag);
    // ScanChain(skimdir+skimver, "2017F_ph", tag);
    // ScanChain(skimdir+skimver, "2017G_ph", tag);

    // ScanChain(skimdir+skimver, "2017B_ll", tag);
    // ScanChain(skimdir+skimver, "2017C_ll", tag);
    // ScanChain(skimdir+skimver, "2017D_ll", tag);
    // ScanChain(skimdir+skimver, "2017E_ll", tag);
    // ScanChain(skimdir+skimver, "2017F_ll", tag);
    // ScanChain(skimdir+skimver, "2017G_ll", tag);

    skimver = "v2_07_metcrtd_2016";
    tag = skimver+tagsuf;

    // gSystem->Exec(Form("mkdir -p output/%s", tag.Data()));

    // ScanChain(skimdir+skimver, "DY", tag);
    // ScanChain(skimdir+skimver, "GJets", tag);
    // ScanChain(skimdir+skimver, "QCD", tag);
    // ScanChain(skimdir+skimver, "ZJets", tag);

    // ScanChain(skimdir+skimver, "qqWG", tag);
    // ScanChain(skimdir+skimver, "WZG", tag);
    // ScanChain(skimdir+skimver, "TGJets", tag);
    // ScanChain(skimdir+skimver, "TTGJets", tag);

    // ScanChain(skimdir+skimver, "TT_2l2nu", tag);
    // ScanChain(skimdir+skimver, "ST", tag);
    // ScanChain(skimdir+skimver, "qqZZ", tag);
    // ScanChain(skimdir+skimver, "qqWZ", tag);

    // ScanChain(skimdir+skimver, "ZGJets_nunu_pTG_40", tag);
    // ScanChain(skimdir+skimver, "ZGJets_nunu_pTG_130", tag);
    // ScanChain(skimdir+skimver, "ZGJets_ll_pTG_40", tag);
    // ScanChain(skimdir+skimver, "ZGJets_ll_pTG_130", tag);

    // ScanChain(skimdir+skimver, "ZGJets_nunu_nlo_pTG_130", tag);
    // ScanChain(skimdir+skimver, "ZGJets_nunu_nlo_incl", tag);
    // ScanChain(skimdir+skimver, "ZGJets_ll_nlo_incl", tag);
    // ScanChain(skimdir+skimver, "ZGJets_ll_nlo_pTG_130", tag);

    // ScanChain(skimdir+skimver, "2016B_ph", tag);
    // ScanChain(skimdir+skimver, "2016C_ph", tag);
    // ScanChain(skimdir+skimver, "2016D_ph", tag);
    // ScanChain(skimdir+skimver, "2016E_ph", tag);
    // ScanChain(skimdir+skimver, "2016F_ph", tag);
    // ScanChain(skimdir+skimver, "2016G_ph", tag);
    // ScanChain(skimdir+skimver, "2016H_ph", tag);

    ScanChain(skimdir+skimver, "2016B_ll", tag);
    ScanChain(skimdir+skimver, "2016C_ll", tag);
    ScanChain(skimdir+skimver, "2016D_ll", tag);
    ScanChain(skimdir+skimver, "2016E_ll", tag);
    ScanChain(skimdir+skimver, "2016F_ll", tag);
    ScanChain(skimdir+skimver, "2016G_ll", tag);
    ScanChain(skimdir+skimver, "2016H_ll", tag);

    // skimver = "v2_06_metcrtd_2016";
    // tag = skimver+"_zgest2_rhen";
    // gSystem->Exec(Form("mkdir -p output/%s", tag.Data()));

    skimdir = "~/working/OffShellAnalyses/loop1/CMSSW_10_2_18/src/CMS3/AnalysisTree/test/output/DYestTest/SkimTrees/v2_05_skim2l_2016";
    // ScanChain(skimdir, "qqWG", tag);
    // ScanChain(skimdir, "WZG", tag);
    // ScanChain(skimdir, "TGJets", tag);
    // ScanChain(skimdir, "TTGJets", tag);
    // ScanChain(skimdir, "ZGJets_ll_nlo_incl", tag);
    // ScanChain(skimdir, "ZGJets_ll_nlo_pTG_130", tag);

    skimdir = "~/working/OffShellAnalyses/loop1/CMSSW_10_2_18/src/CMS3/AnalysisTree/test/output/DYestTest/SkimTrees/v2_06_skimph_rhen_2016";
    // ScanChain(skimdir, "ZGJets_nunu_nlo_incl", tag);
    // ScanChain(skimdir, "ZGJets_nunu_nlo_pTG_130", tag);

}

