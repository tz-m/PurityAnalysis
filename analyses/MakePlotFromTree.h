#ifndef MAKEPLOTFROMTREE_H
#define MAKEPLOTFROMTREE_H

#include "Analysis.h"

class MakePlotFromTree : public Analysis
{
public:
  void run();

  void reconfigureCuts();

private:
  Float_t loamp;
  Float_t hiamp;
};

void MakePlotFromTree::run()
{
  reconfigureCuts();
  TCanvas * canv1 = new TCanvas("mpft","MakePlotFromTree",2000,1600);
  TH1F * amp = new TH1F("amp","Hit Amplitude",300,loamp,hiamp);
  for (auto const & hititr : *(file.GetHitMap()))
    {
      const types::HitInfo * hit = &(hititr.second);
      if (cuts.ChannelPass(hit) &&
          cuts.CounterPass(hit) &&
          cuts.HitPass(hit))
        {
          amp->Fill(hit->amplitude);
        }
    }
  canv1->cd();
  amp->Draw();
  amp->GetXaxis()->SetTitle("(ADC)");
}

void MakePlotFromTree::reconfigureCuts()
{
  loamp = cuts.GetAnalysisCutFloat("MakePlotFromTree","MinVal");
  hiamp = cuts.GetAnalysisCutFloat("MakePlotFromTree","MaxVal");
}

#endif
