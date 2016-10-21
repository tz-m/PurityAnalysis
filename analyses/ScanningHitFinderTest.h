#ifndef SCANHITTEST_H
#define SCANHITTEST_H

#include "Analysis.h"

class ScanningHitFinderTest : public Analysis
{
public:
  void run();

private:
};

void ScanningHitFinderTest::run()
{
  gStyle->SetOptFit(1);
  TF1 * gaus = new TF1("gaus","([0]/([2]*sqrt(2*3.1415926)))*exp(-0.5*(x-[1])*(x-[1])/([2]*[2]))+[3]+x*[4]",0,5300);

  gaus->SetNpx(10000);
  gaus->SetParLimits(2,1,14);
  for (auto const & hititr : *(file->GetHitMap()))
    {
      const types::HitInfo * hit = &(hititr.second);
      if (cuts.ChannelPass(hit)) // &&
      //cuts.CounterPass(hit) &&
      //cuts.HitPass(hit))
        {
          if (!hit->assumedhit) continue;
          TGraph * signal = new TGraph();
          Int_t n=0;
          for (Int_t i = 0; i < hit->signalsize; ++i)
            {
              if (i > hit->begintick-300 && i < hit->endtick+300)
                {
                  signal->SetPoint(n,i,hit->signal[i]);
                  ++n;
                }
            }
          std::cout << "peaktick=" << hit->peaktick << std::endl;
          gaus->SetParameter(1,hit->peaktick);
          gaus->SetParLimits(1,hit->peaktick-5,hit->peaktick+5);
          //signal->Fit(gaus);
          TLine * hitlow = new TLine(hit->begintick,-1000,hit->begintick,1000);
          TLine * hithigh = new TLine(hit->endtick,-1000,hit->endtick,1000);
          TCanvas * canv = new TCanvas("canv","canv",1800,900);
          canv->cd();
          signal->Draw("al");
          signal->Fit(gaus,"B");
          hitlow->Draw("same");
          hithigh->Draw("same");
          if (hit->assumedhit) signal->SetTitle("assumed hit");
          else signal->SetTitle("found hit");
          canv->WaitPrimitive();

          delete canv;
          delete signal;

        }
    }
}

#endif
