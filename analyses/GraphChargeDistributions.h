#ifndef GRAPHCHARGEDISTRIBUTIONS_H
#define GRAPHCHARGEDISTRIBUTIONS_H

#include "Analysis.h"

class GraphChargeDistributions : public Analysis
{
public:
  void run();
};

void GraphChargeDistributions::run()
{
  gStyle->SetOptStat(11);
  TCanvas * canv1 = new TCanvas("gcds1","GraphChargeDistributions1",2000,1600);
  TPad * pad1 = new TPad("pad1","",0,0,1,1);
  TPad * pad2 = new TPad("pad2","",0,0,1,1);
  pad2->SetFillStyle(4000);
  pad1->Draw();
  pad1->cd();
  TH1F * shortd = new TH1F("shortd","Hit Charge Distribution",200,0,20000);
  TH1F * longd = new TH1F("longd","Hit Charge Distribution",200,0,20000);
  for (auto const & hititr : *(file.GetHitMap()))
    {
      const types::HitInfo hit = hititr.second;
      if (cuts.ChannelPass(hit) && cuts.CounterPass(hit) && cuts.HitPass(hit))
        {
          if (hit.c1==30) shortd->Fill(hit.integral/hit.segmentlength);
          if (hit.c1==36) longd->Fill(hit.integral/hit.segmentlength);
        }
    }
  Double_t maxval = std::max(shortd->GetBinContent(shortd->GetMaximumBin()),longd->GetBinContent(longd->GetMaximumBin()));
  shortd->SetLineColor(kBlue);
  shortd->SetLineWidth(2);
  shortd->SetAxisRange(0,maxval*1.1,"Y");
  shortd->Draw();
  shortd->GetXaxis()->SetTitle("Summed ADC / Effective Track Length on Wire");
  pad1->Update();
  TPaveStats * ps1 = (TPaveStats*)shortd->GetListOfFunctions()->FindObject("stats");
  ps1->SetX1NDC(0.4); ps1->SetX2NDC(0.6);
  ps1->SetTextColor(kBlue);
  pad1->Modified();
  canv1->cd();
  Double_t ymin = 0;
  Double_t ymax = maxval*1.1;
  Double_t dy = (ymax-ymin)/0.8;
  Double_t xmin = 0;
  Double_t xmax = 20000;
  Double_t dx = (xmax-xmin)/0.8;
  pad2->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);
  pad2->Draw();
  pad2->cd();
  longd->SetLineColor(kRed);
  longd->SetLineWidth(2);
  longd->Draw("][sames");
  longd->GetXaxis()->SetTitle("Summed ADC / Effective Track Length on Wire");
  pad2->Update();
  TPaveStats * ps2 = (TPaveStats*)longd->GetListOfFunctions()->FindObject("stats");
  ps2->SetX1NDC(0.65); ps2->SetX2NDC(0.85);
  ps2->SetTextColor(kRed);
  TGaxis * axis = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,510,"+L");
  axis->SetLabelFont(42);
  axis->SetLabelSize(0.035);
  axis->Draw();
  TLegend * leg = new TLegend(0.55,0.6,0.85,0.75);
  leg->AddEntry(shortd,"Short Drift (<20cm)","l");
  leg->AddEntry(longd,"Long Drift (>200cm)","l");
  leg->Draw();

}



#endif
