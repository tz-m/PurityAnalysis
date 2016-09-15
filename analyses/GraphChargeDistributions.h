#ifndef GRAPHCHARGEDISTRIBUTIONS_H
#define GRAPHCHARGEDISTRIBUTIONS_H

#include "Analysis.h"

class GraphChargeDistributions : public Analysis
{
public:
  void run();
  void reconfigureCuts();
private:
  TCanvas * DrawCanvas(TH1F * shortd, TH1F * longd, std::string canvname);

  Float_t mincharge;
  Float_t maxcharge;
};

void GraphChargeDistributions::reconfigureCuts()
{
  mincharge = cuts.GetAnalysisCutFloat("GraphChargeDistributions","MinVal");
  std::cout << "Minimum charge window = " << mincharge << std::endl;
  maxcharge = cuts.GetAnalysisCutFloat("GraphChargeDistributions","MaxVal");
  std::cout << "Maximum charge window = " << maxcharge << std::endl;
}


void GraphChargeDistributions::run()
{
  reconfigureCuts();
  TH1F * shortd_real = new TH1F("shortd_real","Hit Charge Distribution -- Real Hits",200,mincharge,maxcharge);
  TH1F * longd_real = new TH1F("longd_real","Hit Charge Distribution -- Real Hits",200,mincharge,maxcharge);
  TH1F * shortd_fake = new TH1F("shortd_fake","Hit Charge Distribution -- Fake Hits",200,mincharge,maxcharge);
  TH1F * longd_fake = new TH1F("longd_fake","Hit Charge Distribution -- Fake Hits",200,mincharge,maxcharge);
  std::cout << "Filling Histograms" << std::endl;
  for (auto const & hititr : *(file.GetHitMap()))
    {
      const types::HitInfo * hit = &(hititr.second);
      if (cuts.ChannelPass(hit) && cuts.CounterPass(hit))
        {
          if (cuts.HitPass(hit))
            {
              if (hit->c1==30) shortd_real->Fill(hit->integral/hit->segmentlength);
              if (hit->c1==36) longd_real->Fill(hit->integral/hit->segmentlength);
            }
          if (!(cuts.HitPass(hit)))
            {
              if (hit->c1==30) shortd_fake->Fill(hit->integral/hit->segmentlength);
              if (hit->c1==36) longd_fake->Fill(hit->integral/hit->segmentlength);
            }
        }
    }
  std::cout << "Finished filling histograms" << std::endl;
  TCanvas * canv1 = DrawCanvas(shortd_real,longd_real,"real");
  canv1->Draw();
  TCanvas * canv2 = DrawCanvas(shortd_fake,longd_fake,"fake");
  canv2->Draw();
}

TCanvas * GraphChargeDistributions::DrawCanvas(TH1F * shortd, TH1F * longd, std::string canvname)
{
  gStyle->SetOptStat(11);
  TCanvas * canv = new TCanvas(canvname.c_str(),"GraphChargeDistributions",2000,1600);
  TPad * pad1 = new TPad("pad1","",0,0,1,1);
  TPad * pad2 = new TPad("pad2","",0,0,1,1);
  pad2->SetFillStyle(4000);
  pad1->Draw();
  pad1->cd();
  Double_t maxval = std::max(shortd->GetBinContent(shortd->GetMaximumBin()),longd->GetBinContent(longd->GetMaximumBin()));
  shortd->SetLineColor(kBlue);
  shortd->SetLineWidth(2);
  shortd->SetAxisRange(0,maxval*1.1,"Y");
  shortd->Draw();
  shortd->GetXaxis()->SetTitle("Summed ADC / Effective Track Length on Wire (ADC/cm)");
  pad1->Update();
  TPaveStats * ps1 = (TPaveStats*)shortd->GetListOfFunctions()->FindObject("stats");
  ps1->SetX1NDC(0.4); ps1->SetX2NDC(0.6);
  ps1->SetTextColor(kBlue);
  pad1->Modified();
  canv->cd();
  Double_t ymin = 0;
  Double_t ymax = maxval*1.1;
  Double_t dy = (ymax-ymin)/0.8;
  Double_t xmin = mincharge;
  Double_t xmax = maxcharge;
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
  return canv;
}



#endif
