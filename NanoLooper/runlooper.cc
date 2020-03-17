#include <iostream>
#include <string>
#include "TChain.h"
#include "TString.h"

using namespace std;

int ScanChain(TChain *ch, string sample, string outdir, int nEventsSample = -1);  // Header for ScanChain.cc

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
    cout << "   Usage: runStopLooper <input_files> <sample> <output_dir> [nevt_samp]" << endl;
    return 1;
  }

  string infile_str(argv[1]);
  string sample(argv[2]);
  string output_dir(argv[3]);
  string nevt_samp;
  if (argc > 4) nevt_samp = string(argv[4]);
  if (argc > 5) string extrargs(argv[5]);

  bool usexrootd = true;  // for current usage

  TChain *ch = new TChain("Events");
  // TString infile = Form("%s/%s*.root", input_dir.c_str(), sample.c_str());

  vector<TString> vecInFiles;
  TString filestr(infile_str);
  while (filestr.Contains(',')) {
    TString fn = filestr(0, filestr.Index(','));
    vecInFiles.push_back( fn );
    filestr.Remove(0, fn.Length()+1);
  }
  vecInFiles.push_back( filestr );

  for (TString file : vecInFiles) {
    if (usexrootd && file.BeginsWith("/store"))
      file = "root://cmsxrootd.fnal.gov/"+file;
    cout << "[runHZZlooper] >> Adding file " << file << " to be process." << endl;
    ch->Add(file);
  }

  if (ch->GetEntries() == 0) {
    cout << "ERROR: no entries in chain. filenames were: " << infile_str << endl;
    return 2;
  }

  int nEvtSamp = (argc > 4)? stoi(nevt_samp) : -1;
  ScanChain(ch, sample, output_dir, nEvtSamp);

  return 0;
}
