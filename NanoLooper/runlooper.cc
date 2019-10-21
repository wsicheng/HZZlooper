#include <iostream>
#include <string>
#include "TChain.h"
#include "TString.h"

using namespace std;

int ScanChain(TChain *ch, string sample, string outdir);  // Header for ScanChain.cc

int main(int argc, char** argv)
{

  if (argc < 4) {
    cout << "   Usage: runStopLooper <input_dir> <sample> <output_dir>" << endl;
    return 1;
  }

  string input_dir(argv[1]);
  string sample(argv[2]);
  string output_dir(argv[3]);

  TChain *ch = new TChain("Events");
  TString infile = Form("%s/%s*.root", input_dir.c_str(), sample.c_str());

  cout << ">>> infile = " << infile << endl;
  ch->Add(infile);

  if (ch->GetEntries() == 0) {
    cout << "ERROR: no entries in chain. filename was: " << infile << endl;
    return 2;
  }

  ScanChain(ch, sample, output_dir);

  return 0;
}
