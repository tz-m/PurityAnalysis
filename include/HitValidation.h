#ifndef HITVALIDATION_H
#define HITVALIDATION_H

#include "ReadHistFile.h"
#include "PedCheck.h"

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TTreeReader.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TMath.h"

#include <string>
#include <map>
#include <algorithm>

class HitValidation {
public:
  struct AvgHits {
    std::vector<std::pair<Float_t,Float_t> > allavg; //for all wires
    std::map<Int_t,std::map<Int_t,std::vector<std::pair<Float_t, Float_t> > > > chanavg; //for each tpc, channel separately
  };

  HitValidation();
  void ImportData(std::string filename);
  void SetPedestalLimits(Float_t pedmeanmin, Float_t pedmeanmax, Float_t pedrmsmin, Float_t pedrmsmax);
  bool GapWire(Int_t channelnumber);
  void CalculateAvgHitsPerWire(AvgHits & ah, Bool_t MC);
  void GraphAvgHitsWires();
  void GraphChargeDistributions();
  void GraphAmplitudeDistribution();

private:
  std::vector<HitInfo> hitmap;
  std::map<Int_t,PedCheck> pedcheckmap;
  std::map<Int_t,Float_t> chanpedrmsmap;

};

HitValidation::HitValidation() {
}

/*
   void HitValidation::ImportData(std::string filename)
   {
   ReadHistFile datafile;
   unsigned int numhits = datafile.ReadFile(filename);
   for (unsigned int i = 0; i < numhits; i++)
    {
      HitInfo * hi = datafile.GetHitInfo(i);
      if (!GapWire(hi->channel))
        {
          if (pedcheck.find(*run) == pedcheck.end())
            {
              PedCheck pc(*run);
              pc.setMinMaxPed(-4000,4000);
              pc.setMinMaxRMS(0,400);
              pedcheck.emplace(std::pair<Int_t,PedCheck>(*run,pc));
            }
          if (chanpedrms.find(*channel) == chanpedrms.end())
            {
              chanpedrms.emplace(std::pair<Int_t,Float_t>(*channel,*rms));
            }
        }
    }
   }
 */
bool HitValidation::GapWire(Int_t channelnumber)
{
  if (channelnumber == 288 || channelnumber == 399 ||
      channelnumber == 400 || channelnumber == 511 ||
      channelnumber == 800 || channelnumber == 911 ||
      channelnumber == 912 || channelnumber == 1023 ||
      channelnumber == 1312 || channelnumber == 1423 ||
      channelnumber == 1424 || channelnumber == 1535 ||
      channelnumber == 1824 || channelnumber == 1935 ||
      channelnumber == 1936 || channelnumber == 2047)
    {
      return true;
    }
  return false;
}

#endif
