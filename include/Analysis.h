#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "ReadHistFile.h"
#include "Cuts.h"
#include "UsefulTypes.h"

#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TMath.h"

#include <string>
#include <map>
#include <algorithm>

class Analysis {
public:
  Bool_t Setup(std::string filename);

protected:
  ReadHistFile file;
  Cuts cuts;
  std::vector<Int_t> runlist;
  Bool_t status;
};

Bool_t Analysis::Setup(std::string filename)
{
  status = true;
  std::cout << "Analysis::Setup() -- Setup analysis: " << std::endl;
  if (file.ReadFile(filename) == 0)
    {
      status = false;
      std::cout << "Analysis::Setup() -- No hits read." << std::endl;
    }
  runlist = file.GetRunList();
  if (runlist.size() == 0)
    {
      status = false;
      std::cout << "Analysis::Setup() -- No runs read." << std::endl;
    }
  if (!(cuts.PrepareCuts("/home/mthiesse/PurityAnalysis/cuts.conf",runlist)))
    {
      status = false;
      std::cout << "Analysis::Setup() -- Cuts not set up." << std::endl;
    }
  return status;
}


#endif
