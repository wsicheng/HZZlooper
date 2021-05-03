
{
    gSystem->Exec("mkdir -p plots");

    gROOT->ProcessLine(".L SkimTree.cc+");
    gROOT->ProcessLine(".L ScanChain.C+");

    TString skimdir = "/hadoop/cms/store/user/sicheng/Offshell_2L2Nu/Worker/output/DYestTest/SkimTrees/";
    TString skimver = "v2_08_metcrtd_2018";
    // TString skimver = "v2_10_fullsamp_2018";

    // TString tagsuf = "_phtest_allncol";
    TString tagsuf = "_zgest5";
    TString tag = skimver+tagsuf;

    // gSystem->Exec(Form("mkdir -p output/%s", tag.Data()));

    ScanChain(skimdir+skimver, "DY", tag);
    ScanChain(skimdir+skimver, "GJets", tag);
    ScanChain(skimdir+skimver, "QCD", tag);
    ScanChain(skimdir+skimver, "ZJets", tag);
    // ScanChain(skimdir+skimver, "WJets", tag);

    ScanChain(skimdir+skimver, "qqWG", tag);
    ScanChain(skimdir+skimver, "WZG", tag);
    ScanChain(skimdir+skimver, "TGJets", tag);
    ScanChain(skimdir+skimver, "TTGJets", tag);
    gSystem->Exec(Form("hadd -f output/%s/others.root output/%s/WZG.root output/%s/*TGJets.root", tag.Data(), tag.Data(), tag.Data()));

    ScanChain(skimdir+skimver, "TT_2l2nu", tag);
    ScanChain(skimdir+skimver, "ST", tag);
    ScanChain(skimdir+skimver, "qqZZ", tag);
    ScanChain(skimdir+skimver, "qqWZ", tag);

    ScanChain(skimdir+skimver, "ZGJets_nunu_pTG_40", tag);
    ScanChain(skimdir+skimver, "ZGJets_nunu_pTG_130", tag);
    ScanChain(skimdir+skimver, "ZGJets_ll_pTG_40", tag);
    ScanChain(skimdir+skimver, "ZGJets_ll_pTG_130", tag);

    ScanChain(skimdir+skimver, "ZGJets_nunu_nlo_incl", tag);
    ScanChain(skimdir+skimver, "ZGJets_nunu_nlo_pTG_130", tag);
    ScanChain(skimdir+skimver, "ZGJets_ll_nlo_incl", tag);
    ScanChain(skimdir+skimver, "ZGJets_ll_nlo_pTG_130", tag);

    // // ScanChain(skimdir+skimver, "2018A_ph", tag);
    // // ScanChain(skimdir+skimver, "2018B_ph", tag);
    // // ScanChain(skimdir+skimver, "2018C_ph", tag);
    // // ScanChain(skimdir+skimver, "2018D_ph", tag);
    // // gSystem->Exec(Form("hadd -f output/%s/data_2018_phskim.root output/%s/2018?_ph.root", tag.Data(), tag.Data()));

    ScanChain(skimdir+skimver, "2018A_ll", tag);
    ScanChain(skimdir+skimver, "2018B_ll", tag);
    ScanChain(skimdir+skimver, "2018C_ll", tag);
    ScanChain(skimdir+skimver, "2018D_ll", tag);
    gSystem->Exec(Form("hadd -f output/%s/data_2018_llskim.root output/%s/2018?_ll.root", tag.Data(), tag.Data()));

    skimver2 = "v2_11_fullsamp_2018";
    // ScanChain(skimdir+skimver2, "WJets_lnu_0j", tag);
    // ScanChain(skimdir+skimver2, "WJets_lnu_1j", tag);
    // ScanChain(skimdir+skimver2, "WJets_lnu_2j", tag);
    // gSystem->Exec(Form("hadd -f output/%s/WJets.root output/%s/WJets_lnu_?j.root", tag.Data(), tag.Data()));

    tag = skimver+"_phtest_tslt1";

    // ScanChain(skimdir+skimver2, "egamma_2018A", tag);
    // ScanChain(skimdir+skimver2, "egamma_2018B", tag);
    // ScanChain(skimdir+skimver2, "egamma_2018C", tag);
    // ScanChain(skimdir+skimver2, "egamma_2018D", tag);
    // gSystem->Exec(Form("hadd -f output/%s/data_2018_llskim.root output/%s/egamma_2018?.root", tag.Data(), tag.Data()));

    // skimver = "v2_10_fullsamp_2017";
    skimver = "v2_08_metcrtd_2017";
    tag = skimver+tagsuf;

    gSystem->Exec(Form("mkdir -p output/%s", tag.Data()));

    ScanChain(skimdir+skimver, "DY", tag);
    ScanChain(skimdir+skimver, "GJets", tag);
    ScanChain(skimdir+skimver, "QCD", tag);
    ScanChain(skimdir+skimver, "ZJets", tag);

    ScanChain(skimdir+skimver, "qqWG", tag);
    ScanChain(skimdir+skimver, "WZG", tag);
    ScanChain(skimdir+skimver, "TGJets", tag);
    ScanChain(skimdir+skimver, "TTGJets", tag);
    ScanChain(skimdir+skimver, "TT_2l2nu", tag);

    ScanChain(skimdir+skimver, "ST", tag);
    ScanChain(skimdir+skimver, "qqZZ", tag);
    ScanChain(skimdir+skimver, "qqWZ", tag);

    ScanChain(skimdir+skimver, "ZGJets_nunu_pTG_40", tag);
    ScanChain(skimdir+skimver, "ZGJets_nunu_pTG_130", tag);
    ScanChain(skimdir+skimver, "ZGJets_ll_pTG_40", tag);
    ScanChain(skimdir+skimver, "ZGJets_ll_pTG_130", tag);

    ScanChain(skimdir+skimver, "ZGJets_nunu_nlo_incl", tag);
    ScanChain(skimdir+skimver, "ZGJets_nunu_nlo_pTG_130", tag);
    ScanChain(skimdir+skimver, "ZGJets_ll_nlo_incl", tag);
    ScanChain(skimdir+skimver, "ZGJets_ll_nlo_pTG_130", tag);

    // ScanChain(skimdir+skimver, "2017B_ph", tag);
    // ScanChain(skimdir+skimver, "2017C_ph", tag);
    // ScanChain(skimdir+skimver, "2017D_ph", tag);
    // ScanChain(skimdir+skimver, "2017E_ph", tag);
    // ScanChain(skimdir+skimver, "2017F_ph", tag);
    // gSystem->Exec(Form("hadd -f output/%s/data_2017_phskim.root output/%s/2017?_ph.root", tag.Data(), tag.Data()));

    ScanChain(skimdir+skimver, "2017B_ll", tag);
    ScanChain(skimdir+skimver, "2017C_ll", tag);
    ScanChain(skimdir+skimver, "2017D_ll", tag);
    ScanChain(skimdir+skimver, "2017E_ll", tag);
    ScanChain(skimdir+skimver, "2017F_ll", tag);
    gSystem->Exec(Form("hadd -f output/%s/data_2017_llskim.root output/%s/2017?_ll.root", tag.Data(), tag.Data()));


    skimver2 = "v2_11_fullsamp_2017";

    // ScanChain(skimdir+skimver2, "WJets_lnu_0j", tag);
    // ScanChain(skimdir+skimver2, "WJets_lnu_1j", tag);
    // ScanChain(skimdir+skimver2, "WJets_lnu_2j", tag);
    // gSystem->Exec(Form("hadd -f output/%s/WJets.root output/%s/WJets_lnu_?j.root", tag.Data(), tag.Data()));

    // tag = skimver+"_phtest_tslt1";
    // ScanChain(skimdir+skimver2, "photon_2017B", tag);
    // ScanChain(skimdir+skimver2, "photon_2017C", tag);
    // ScanChain(skimdir+skimver2, "photon_2017D", tag);
    // ScanChain(skimdir+skimver2, "photon_2017E", tag);
    // ScanChain(skimdir+skimver2, "photon_2017F", tag);
    // gSystem->Exec(Form("hadd -f output/%s/data_2017_phskim.root output/%s/photon_2017?.root", tag.Data(), tag.Data()));

    skimver = "v2_08_metcrtd_2016";
    // skimver = "v2_10_fullsamp_2016";
    tag = skimver+tagsuf;

    gSystem->Exec(Form("mkdir -p output/%s", tag.Data()));

    // ScanChain(skimdir+skimver, "GJets", tag);
    // ScanChain(skimdir+skimver, "QCD", tag);
    // ScanChain(skimdir+skimver, "ZJets", tag);

    ScanChain(skimdir+skimver, "qqWG", tag);
    ScanChain(skimdir+skimver, "WZG", tag);
    ScanChain(skimdir+skimver, "TGJets", tag);
    ScanChain(skimdir+skimver, "TTGJets", tag);
    gSystem->Exec(Form("hadd -f output/%s/others.root output/%s/WZG.root output/%s/*TGJets.root", tag.Data(), tag.Data(), tag.Data()));

    ScanChain(skimdir+skimver, "ZGJets_nunu_pTG_40", tag);
    ScanChain(skimdir+skimver, "ZGJets_nunu_pTG_130", tag);
    ScanChain(skimdir+skimver, "ZGJets_ll_pTG_40", tag);
    ScanChain(skimdir+skimver, "ZGJets_ll_pTG_130", tag);

    ScanChain(skimdir+skimver, "ZGJets_nunu_nlo_pTG_130", tag);
    ScanChain(skimdir+skimver, "ZGJets_nunu_nlo_incl", tag);
    ScanChain(skimdir+skimver, "ZGJets_ll_nlo_incl", tag);
    ScanChain(skimdir+skimver, "ZGJets_ll_nlo_pTG_130", tag);
    gSystem->Exec(Form("hadd -f output/%s/ZGJets.root output/%s/ZGJets_*_nlo*.root", tag.Data(), tag.Data()));

    ScanChain(skimdir+skimver, "DY", tag);
    ScanChain(skimdir+skimver, "TT_2l2nu", tag);
    ScanChain(skimdir+skimver, "ST", tag);
    ScanChain(skimdir+skimver, "qqZZ", tag);
    ScanChain(skimdir+skimver, "qqWZ", tag);

    // ScanChain(skimdir+skimver, "2016B_ph", tag);
    // ScanChain(skimdir+skimver, "2016C_ph", tag);
    // ScanChain(skimdir+skimver, "2016D_ph", tag);
    // ScanChain(skimdir+skimver, "2016E_ph", tag);
    // ScanChain(skimdir+skimver, "2016F_ph", tag);
    // ScanChain(skimdir+skimver, "2016G_ph", tag);
    // ScanChain(skimdir+skimver, "2016H_ph", tag);
    // gSystem->Exec(Form("hadd -f output/%s/data_2016_phskim.root output/%s/2016?_ph.root", tag.Data(), tag.Data()));

    ScanChain(skimdir+skimver, "2016B_ll", tag);
    ScanChain(skimdir+skimver, "2016C_ll", tag);
    ScanChain(skimdir+skimver, "2016D_ll", tag);
    ScanChain(skimdir+skimver, "2016E_ll", tag);
    ScanChain(skimdir+skimver, "2016F_ll", tag);
    ScanChain(skimdir+skimver, "2016G_ll", tag);
    ScanChain(skimdir+skimver, "2016H_ll", tag);
    gSystem->Exec(Form("hadd -f output/%s/data_2016_llskim.root output/%s/2016?_ll.root", tag.Data(), tag.Data()));

    // tagsuf += "_allncol_conv";

    // skimver = "v2_11_fullsamp_2016";
    // tag = skimver+tagsuf;
    skimver2 = "v2_11_fullsamp_2016";

    tag = skimver+"_phtest_tslt1";

    // ScanChain(skimdir+skimver2, "photon_2016B", tag);
    // ScanChain(skimdir+skimver2, "photon_2016C", tag);
    // ScanChain(skimdir+skimver2, "photon_2016D", tag);
    // ScanChain(skimdir+skimver2, "photon_2016E", tag);
    // ScanChain(skimdir+skimver2, "photon_2016F", tag);
    // ScanChain(skimdir+skimver2, "photon_2016G", tag);
    // ScanChain(skimdir+skimver2, "photon_2016H", tag);
    // gSystem->Exec(Form("hadd -f output/%s/data_2016_phskim.root output/%s/photon_2016?.root", tag.Data(), tag.Data()));

    // ScanChain(skimdir+skimver2, "WJets", tag);

    // tag = skimver+tagsuf;



    // skimver = "v2_06_metcrtd_2016";
    // tag = skimver+"_zgest2_rhen";
    // gSystem->Exec(Form("mkdir -p output/%s", tag.Data()));

    skimdir = "/home/users/sicheng/working/OffShellAnalyses/loop2/CMSSW_10_2_18/src/CMS3/AnalysisTree/test/output/DYestTest/SkimTrees/v2_11_skimph_2016";
    // ScanChain(skimdir, "photon_2016B", tag);
    // ScanChain(skimdir, "photon_2016C", tag);
    // ScanChain(skimdir, "photon_2016D", tag);
    // ScanChain(skimdir, "photon_2016E", tag);
    // ScanChain(skimdir, "photon_2016F", tag);
    // ScanChain(skimdir, "photon_2016G", tag);
    // ScanChain(skimdir, "photon_2016H", tag);

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

