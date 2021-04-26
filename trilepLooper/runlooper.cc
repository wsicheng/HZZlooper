#include <iostream>
#include <string>
#include "TChain.h"
#include "TString.h"
#include "TSystem.h"

using namespace std;

// int ScanChain(TChain* ch, string sample, string outdir, string tn = "Events", double sumwgt = 1.0, string portion = "");  // Header for ScanChain.cc
int ScanChain(TChain* ch, TString sample, TString tag, TString systype = "", TString specifiers = "", TString extrargs = "");

TString parseArg(const TString& input, TString arg, const TString dfval="") {
  if (!arg.EndsWith("=")) arg += "=";
  if (input.Contains(arg)) {
    int sidx = input.Index(arg) + arg.Sizeof() - 1;
    int eidx = input.Index(",", sidx);
    if (eidx < 0) eidx = input.Sizeof();
    return input(sidx, eidx-sidx);
  } else {
    return dfval;
  }
}

int main(int argc, char** argv)
{

  if (argc < 4) {
    cout << "   Usage: runLooper <input_dir> <sample> <tag> [systype] [specifiers] [extra]" << endl;
    return 1;
  }

  TString indir(argv[1]);
  TString sample(argv[2]);
  TString tag(argv[3]);
  TString systype = "Nominal";
  if (argc > 4) systype = TString(argv[4]);

  TString specifier = "";
  if (argc > 5) specifier = TString(argv[5]);

  TString extrargs;
  if (argc > 6) extrargs = TString(argv[6]);

  TString treename = parseArg(extrargs, "treename", "SkimTree");
  if (treename == "Events") treename = "cms3ntuple/Events";
  TChain *ch = new TChain(treename);
  if (sample.Contains(".root")) {
    vector<TString> vecInFiles;
    TString filestr(sample);
    while (filestr.Contains(',')) {
      TString fn = filestr(0, filestr.Index(','));
      vecInFiles.push_back( fn );
      filestr.Remove(0, fn.Length()+1);
    }
    vecInFiles.push_back( filestr );
    for (TString file : vecInFiles) {
      cout << "[runHZZlooper] >> Adding file " << file << " to be process." << endl;
      ch->Add(file);
    }
  } else {
    TString files_in = Form("%s/%s*.root", indir.Data(), sample.Data());
    cout << ">> Adding " << files_in << " into the chain." << endl;
    ch->Add(files_in);

    if (ch->GetEntries() == 0) {
      cout << "ERROR: no entries in chain for " << files_in << endl;
      return 2;
    }
  }

  ScanChain(ch, sample, tag, systype, specifier, extrargs);

  return 0;
}
