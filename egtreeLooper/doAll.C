
{
    // gSystem->Exec("mkdir -p plots");
    gSystem->Exec("mkdir -p output");

    gROOT->ProcessLine(".L EGTree.cc+");
    gROOT->ProcessLine(".L ScanChain.C+");

    // TString skimdir = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/SinglePhotonTriggerEfficiencies/SkimTrees/200624/WithPUJetId_NoTightLeptonJetId/";
    // TString skimdir = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/SinglePhotonEGTnP/SkimTrees/201007";
    // TString skimdir = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/SinglePhotonEGTnP/SkimTrees/201128";
    TString skimdir = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/SinglePhotonEGTnP/SkimTrees/210107/AK4Jets_WithPUJetId_NoTightLeptonJetId_ParticleCleaned/PFMET_WithXY_WithPartMomCorr_P4Preserved";

    // gSystem->Exec(Form("mkdir -p output/%s", tag.Data()));
    ScanChain(skimdir, "EGamma", "2018A");
    // ScanChain(skimdir, "EGamma", "2018B");
    // ScanChain(skimdir, "EGamma", "2018C");
    // ScanChain(skimdir, "EGamma", "2018D");

    // ScanChain(skimdir, "SingleElectron", "2017B");
    // ScanChain(skimdir, "SingleElectron", "2017C");
    // ScanChain(skimdir, "SingleElectron", "2017D");
    // ScanChain(skimdir, "SingleElectron", "2017E");
    // ScanChain(skimdir, "SingleElectron", "2017F");

    // ScanChain(skimdir, "SingleElectron", "2016B");
    // ScanChain(skimdir, "SingleElectron", "2016C");
    // ScanChain(skimdir, "SingleElectron", "2016D");
    // ScanChain(skimdir, "SingleElectron", "2016E");
    // ScanChain(skimdir, "SingleElectron", "2016F");
    // ScanChain(skimdir, "SingleElectron", "2016G");
    // ScanChain(skimdir, "SingleElectron", "2016H");

    // ScanChain(skimdir, "DoubleEG", "2017B");
    // ScanChain(skimdir, "DoubleEG", "2017C");
    // ScanChain(skimdir, "DoubleEG", "2017D");
    // ScanChain(skimdir, "DoubleEG", "2017E");
    // ScanChain(skimdir, "DoubleEG", "2017F");

    // ScanChain(skimdir, "DoubleEG", "2016B");
    // ScanChain(skimdir, "DoubleEG", "2016C");
    // ScanChain(skimdir, "DoubleEG", "2016D");
    // ScanChain(skimdir, "DoubleEG", "2016E");
    // ScanChain(skimdir, "DoubleEG", "2016F");
    // ScanChain(skimdir, "DoubleEG", "2016G");
    // ScanChain(skimdir, "DoubleEG", "2016H");

    // ScanChain(skimdir, "DYJetsToLL", "2016");
    // ScanChain(skimdir, "DYJetsToLL", "2017");
    // ScanChain(skimdir, "DYJetsToLL", "2018");

    // skimdir = "/home/users/sicheng/working/OffShellAnalyses/loop2/CMSSW_10_2_18/src/CMS3/AnalysisTree/test/output/SinglePhotonTriggerEfficiencies/SkimTrees/200622/WithPUJetId_NoTightLeptonJetId/";

    // ScanChain(skimdir, "DYJetsToLL", "2017");
    // ScanChain(skimdir, "DoubleEG", "2017B");
    // ScanChain(skimdir, "DoubleEG", "2017C");
    // ScanChain(skimdir, "DoubleEG", "2017D");
    // ScanChain(skimdir, "DoubleEG", "2017E");
    // ScanChain(skimdir, "DoubleEG", "2017F");


}

