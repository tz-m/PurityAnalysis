#include "HitValidation.h"

void HitValidation::CalculateAvgHitsPerWire(AvgHits & ah, Bool_t MC)
{
  ah.allavg.clear();
  ah.allavg.resize((MC) ? 6 : 7); // 7 for data, 6 for MC (because geometry differs)
  ah.chanavg.clear();

  std::map<UInt_t,std::map<Int_t, std::vector<Int_t> > > cpairs;
  std::map<UInt_t,Float_t> counterx;
  for (auto & hit : hits)
    {
      //if (!(hit.countercut)) continue;
      //if (!(hit.fitsuccess)) continue;
      //if (!(hit.fitrealhit)) continue;
      if (!(((hit.c1>=6 && hit.c1<=15) || (hit.c1>=28 && hit.c1<=37)) && (hit.c1%22==hit.c2 || hit.c2%22==hit.c1))) continue;
      if (std::find(cpairs[hit.c1][hit.run].begin(),cpairs[hit.c1][hit.run].end(),hit.event) == cpairs[hit.c1][hit.run].end()) cpairs[hit.c1][hit.run].push_back(hit.event);

      if (MC) {
          if (hit.c1==6 || hit.c1==28) counterx[hit.c1] = -33.5;
          if (hit.c1==7 || hit.c1==29) counterx[hit.c1] = -4;
          if (hit.c1==8 || hit.c1==30) counterx[hit.c1] = 25;
          if (hit.c1==9 || hit.c1==31) counterx[hit.c1] = 55;
          if (hit.c1==10 || hit.c1==32) counterx[hit.c1] = 85;
          if (hit.c1==11 || hit.c1==33) counterx[hit.c1] = 115;
          if (hit.c1==12 || hit.c1==34) counterx[hit.c1] = 145;
          if (hit.c1==13 || hit.c1==35) counterx[hit.c1] = 175;
          if (hit.c1==14 || hit.c1==36) counterx[hit.c1] = 205;
          if (hit.c1==15 || hit.c1==37) counterx[hit.c1] = 235;
        } else {
          counterx[hit.c1] = hit.c1x;
        }
    }

  std::map<UInt_t,Int_t> cpairEvents;
  for (auto const & cpair : cpairs)
    {
      std::cout << "Counter pair " << cpair.first;
      Int_t count = 0;
      for (auto const & run : cpair.second)
        {
          count += run.second.size();
        }
      cpairEvents[cpair.first] = count;
      std::cout << " contains " << count << " events." << std::endl;
    }

  std::map<UInt_t,Float_t> cpairHits;
  std::map<UInt_t,std::map<Int_t,std::map<Int_t,Float_t> > >cpairChanHits;
  for (auto & hit : hits)
    {
      if (!(hit.countercut)) continue;
      if (!(hit.fitsuccess)) continue;
      if (!(hit.fitrealhit)) continue;
      if (fabs(hit.peaktime-hit.peaktimeFilter)>5) continue;
      if (!(((hit.c1>=6 && hit.c1<=15) || (hit.c1>=28 && hit.c1<=37)) && (hit.c1%22==hit.c2 || hit.c2%22==hit.c1))) continue;
      bool goodped = pedcheck[hit.run].check(hit.channel);
      if (!MC && !goodped) continue;
      //if (hit.amplitude < 25) continue;
      cpairHits[hit.c1] += (cpairEvents[hit.c1] >= 1 ? 1.0/cpairEvents[hit.c1] : 0);
      cpairChanHits[hit.c1][hit.tpc][hit.channel] += (cpairEvents[hit.c1] >= 1 ? 1.0/cpairEvents[hit.c1] : 0);
    }

  for (auto const & cpair : cpairHits)
    {
      if (cpair.first >= 30 && cpair.first <= ((MC) ? 35 : 36))
        {
          ah.allavg.push_back(std::pair<Float_t,Float_t>(counterx[cpair.first],cpair.second));
          for (auto const & cpchtpc : cpairChanHits[cpair.first])
            {
              for (auto const & cpch : cpchtpc.second)
                {
                  ah.chanavg[cpchtpc.first][cpch.first].push_back(std::pair<Float_t,Float_t>(counterx[cpair.first],cpch.second));
                }
            }
        }

      std::cout << "Avg NHits in cpair " << cpair.first << " = " << cpair.second << std::endl;
    }

}

void HitValidation::SetPedestalLimits(Float_t pedmeanmin, Float_t pedmeanmax, Float_t pedrmsmin, Float_t pedrmsmax)
{
  for (auto & ped : pedcheck)
    {
      ped.second.setMinMaxRMS(pedrmsmin,pedrmsmax);
      ped.second.setMinMaxPed(pedmeanmin,pedmeanmax);
    }
}

void HitValidation::GraphAvgHitsWires()
{
  Bool_t MC = false;
  TCanvas * canv1 = new TCanvas("canv1","avghits",2000,1600);
  TCanvas * canv2 = new TCanvas("canv2","channelInfo",2000,1600);
  AvgHits ah;
  TMultiGraph * mg = new TMultiGraph();
  TGraphErrors * effrms = new TGraphErrors();
  Int_t i_effrms = 0;
  //TH1F * chanrms = new TH1F("chanrms","Channel RMS",200,0,0);
  TF1 * line = new TF1("line","pol1",0,250);
  CalculateAvgHitsPerWire(ah,MC);
  for (auto const & chtpc : ah.chanavg)
    {
      if (chtpc.first % 2 == 0) continue;
      if (chtpc.first == 3 || chtpc.first == 5) continue;
      for (auto const & ch : chtpc.second)
        {
          TGraph * gr = new TGraph(ch.second.size());
          gr->SetLineColor(chtpc.first+2);
          Int_t i = 0;
          std::vector<double> chsecond;
          for (auto const & ahc : ch.second)
            {
              chsecond.push_back(ahc.second);
              gr->SetPoint(i,ahc.first,ahc.second);
              i++;
            }
          effrms->SetPoint(i_effrms,chanpedrms[ch.first],TMath::Mean(chsecond.size(),chsecond.data()));
          effrms->SetPointError(i_effrms,0,TMath::RMS(chsecond.size(),chsecond.data()));
          i_effrms++;
          //if (gr->GetN()==((MC) ? 6 : 7)) gr->Fit("line","q");
          if (gr->GetN()==((MC) ? 6 : 7) && line->GetParameter(0) < 10.1 /* && line->GetParameter(1)<0.0001*/)
            {
              //if (chanpedrms[ch.first] > 6 && chanpedrms[ch.first] < 20)
              //if (line->GetParameter(0)<0.1)
              {
                //std::cout << ch.first << std::endl;
                mg->Add(gr);
                //chanrms->Fill(chanpedrms[ch.first]);
              }
              //mg->Add(gr);
            }
        }
    }
  canv1->cd();
  mg->SetTitle("Channel Efficiency");
  mg->Draw("al");
  mg->GetXaxis()->SetTitle("Drift Distance (cm)");
  mg->GetYaxis()->SetTitle("Mean # hits per event");
  mg->Draw("al");
  TGraph * tpc0 = new TGraph();
  tpc0->SetLineColor(2);
  tpc0->SetLineWidth(2);
  TGraph * tpc1 = new TGraph();
  tpc1->SetLineColor(3);
  tpc1->SetLineWidth(2);
  TGraph * tpc2 = new TGraph();
  tpc2->SetLineColor(4);
  tpc2->SetLineWidth(2);
  TGraph * tpc3 = new TGraph();
  tpc3->SetLineColor(5);
  tpc3->SetLineWidth(2);
  TGraph * tpc4 = new TGraph();
  tpc4->SetLineColor(6);
  tpc4->SetLineWidth(2);
  TGraph * tpc5 = new TGraph();
  tpc5->SetLineColor(7);
  tpc5->SetLineWidth(2);
  TGraph * tpc6 = new TGraph();
  tpc6->SetLineColor(8);
  tpc6->SetLineWidth(2);
  TGraph * tpc7 = new TGraph();
  tpc7->SetLineColor(9);
  tpc7->SetLineWidth(2);
  TLegend * leg = new TLegend(0.85,0.7,0.95,0.95);
  //leg->AddEntry(tpc0,"TPC0","l");
  leg->AddEntry(tpc1,"TPC1","l");
  //leg->AddEntry(tpc2,"TPC2","l");
  leg->AddEntry(tpc3,"TPC3","l");
  //leg->AddEntry(tpc4,"TPC4","l");
  leg->AddEntry(tpc5,"TPC5","l");
  //leg->AddEntry(tpc6,"TPC6","l");
  leg->AddEntry(tpc7,"TPC7","l");
  leg->Draw();
  canv2->cd();
  effrms->SetMarkerStyle(20);
  effrms->SetMarkerSize(1);
  effrms->Draw("AP");
  effrms->SetTitle("Channel Efficiency vs. Noise");
  effrms->GetXaxis()->SetTitle("RMS Channel Noise (ADC)");
  effrms->GetYaxis()->SetTitle("Mean # hits per event");
}

void HitValidation::GraphChargeDistributions()
{
  gStyle->SetOptStat(11);
  TCanvas * canv1 = new TCanvas("canv1","chargedistributions",2000,1600);
  TPad * pad1 = new TPad("pad1","",0,0,1,1);
  TPad * pad2 = new TPad("pad2","",0,0,1,1);
  pad2->SetFillStyle(4000);
  pad1->Draw();
  pad1->cd();
  TH1F * shortd = new TH1F("shortd","Hit Charge Distribution",200,0,20000);
  TH1F * longd = new TH1F("longd","Hit Charge Distribution",200,0,20000);
  for (auto const & hit : hits)
    {
      if (!hit.fitrealhit && hit.c1==30 /*&& hit.amplitude>25*/) shortd->Fill(hit.integral/hit.segmentlength);
      if (!hit.fitrealhit && hit.c1==36 /*&& hit.amplitude>25*/) longd->Fill(hit.integral/hit.segmentlength);
    }
  double maxval = std::max(shortd->GetBinContent(shortd->GetMaximumBin()),longd->GetBinContent(longd->GetMaximumBin()));
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

void HitValidation::GraphAmplitudeDistribution()
{
  TCanvas * canv1 = new TCanvas("canv1","amplitude",2000,1600);
  TH1F * amp = new TH1F("amp","Hit Amplitude",300,0,600);
  for (auto const & hit : hits)
    {
      if (hit.fitrealhit) amp->Fill(hit.amplitude);
    }
  canv1->cd();
  amp->Draw();
  amp->GetXaxis()->SetTitle("(#ADC)");
}

int run(std::string filename)
{
  std::cout << "Filename = " << filename << std::endl;
  HitValidation hv;
  hv.ImportData(filename);

  //hv.GraphAmplitudeDistribution();

  //hv.GraphChargeDistributions();

  hv.GraphAvgHitsWires();

  //std::vector<std::pair<Float_t,Float_t> > avghits;
  //hv.CalculateAvgHitsPerWire(avghits,450);
  //for (auto const & a : avghits)
  //  {
  //    std::cout << "drift distance = " << a.first << "  nhits = " << a.second << std::endl;
  //  }
  return 0;
}

int main(int argc, char** argv)
{
  std::string filename = argv[1];
  TApplication app("hitvalidation",&argc,argv);
  app.ExitOnException();
  run(filename);
  app.Run();
  return 0;
}
