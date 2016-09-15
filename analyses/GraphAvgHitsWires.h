#ifndef GRAPHAVGHITSWIRES_H
#define GRAPHAVGHITSWIRES_H

#include "Analysis.h"

class GraphAvgHitsWires : public Analysis
{
public:
  void run();
  void reconfigureCuts();
private:
  Bool_t MC;
};

void GraphAvgHitsWires::reconfigureCuts()
{
  MC = cuts.GetAnalysisCutBool("GraphAvgHitsWires","MC");
}

void GraphAvgHitsWires::run()
{
  reconfigureCuts();

  types::EfficiencyGraph allavg;
  types::TPCEfficiencyMap chanavg;

  allavg.resize((MC) ? 6 : 7);

  std::map<UInt_t, std::map<Int_t, std::vector<Int_t> > > neventspercpair;
  std::map<UInt_t, Float_t> counterx;

  std::cout << "GraphAvgHitsWires::run() -- Counting " << std::endl;
  for (auto const & hititr : *(file.GetHitMap()))
    {
      const types::HitInfo * hit = &(hititr.second);
      if (cuts.ChannelPass(hit) && cuts.CounterPass(hit))
        {
          if (std::find(neventspercpair[hit->c1][hit->run].begin(),
                        neventspercpair[hit->c1][hit->run].end(),
                        hit->event) ==
              neventspercpair[hit->c1][hit->run].end())
            {
              neventspercpair[hit->c1][hit->run].push_back(hit->event);
            }
          if (MC) {
              if (hit->c1==6 || hit->c1==28) counterx[hit->c1] = -33.5;
              if (hit->c1==7 || hit->c1==29) counterx[hit->c1] = -4;
              if (hit->c1==8 || hit->c1==30) counterx[hit->c1] = 25;
              if (hit->c1==9 || hit->c1==31) counterx[hit->c1] = 55;
              if (hit->c1==10 || hit->c1==32) counterx[hit->c1] = 85;
              if (hit->c1==11 || hit->c1==33) counterx[hit->c1] = 115;
              if (hit->c1==12 || hit->c1==34) counterx[hit->c1] = 145;
              if (hit->c1==13 || hit->c1==35) counterx[hit->c1] = 175;
              if (hit->c1==14 || hit->c1==36) counterx[hit->c1] = 205;
              if (hit->c1==15 || hit->c1==37) counterx[hit->c1] = 235;
            } else {
              counterx[hit->c1] = hit->c1x;
            }
        }
    }

  std::map<UInt_t,Int_t> cpairEvents;
  for (auto const & cpair : neventspercpair)
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
  std::map<UInt_t,std::map<Int_t,std::map<Int_t,Float_t> > > cpairChanHits;
  for (auto const & hititr : *(file.GetHitMap()))
    {
      const types::HitInfo * hit = &(hititr.second);
      if (cuts.ChannelPass(hit) && cuts.CounterPass(hit) && cuts.HitPass(hit))
        {
          cpairHits[hit->c1] += (cpairEvents[hit->c1] >= 1 ? 1.0/cpairEvents[hit->c1] : 0);
          cpairChanHits[hit->c1][hit->tpc][hit->channel] += (cpairEvents[hit->c1] >= 1 ? 1.0/cpairEvents[hit->c1] : 0);
        }
    }

  for (auto const & cpair : cpairHits)
    {
      if (cpair.first >= 30 && cpair.first <= ((MC) ? 35 : 36))
        {
          allavg.push_back(std::make_pair(counterx[cpair.first],cpair.second));
          for (auto const & cpchtpc : cpairChanHits[cpair.first])
            {
              for (auto const & cpch : cpchtpc.second)
                {
                  chanavg[cpchtpc.first][cpch.first].push_back(std::make_pair(counterx[cpair.first],cpch.second));
                }
            }
          std::cout << "Avg N hits in cpair " << cpair.first << " = " << cpair.second << std::endl;
        }
    }

  TCanvas * canv1 = new TCanvas("canv1","avghits",2000,1600);
  TCanvas * canv2 = new TCanvas("canv2","channelInfo",2000,1600);
  TMultiGraph * mg = new TMultiGraph();
  TGraphErrors * effrms = new TGraphErrors();
  Int_t i_effrms = 0;
  //TF1 * line = new TF1("line","pol1",0,250);

  std::map<Int_t,Float_t> chanpedmeanrms;
  std::map<Int_t,Float_t> chanpedrmsrms;

  std::cout << "Making some plots" << std::endl;
  for (auto const & chtpc : chanavg)
    {
      if (chtpc.first % 2 == 0) continue;
      //if (chtpc.first == 3 || chtpc.first == 5) continue;
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
          if (chanpedmeanrms.find(ch.first) == chanpedmeanrms.end())
            chanpedmeanrms[ch.first] = cuts.GetPedCheckPtr()->PedMeanRms(ch.first);
          if (chanpedrmsrms.find(ch.first) == chanpedrmsrms.end())
            chanpedrmsrms[ch.first] = cuts.GetPedCheckPtr()->PedRmsRms(ch.first);
          Float_t pedrms = chanpedmeanrms[ch.first];
          Float_t pedrmserr = chanpedrmsrms[ch.first];
          // cuts.GetPedCheckPtr()->PedRms(,ch.first);
          //std::cout << "Channel " << ch.first << "  run 15600  pedrms " << pedrms << std::endl;
          effrms->SetPoint(i_effrms,pedrms,TMath::Mean(chsecond.size(),chsecond.data()));
          effrms->SetPointError(i_effrms,pedrmserr,TMath::RMS(chsecond.size(),chsecond.data()));
          i_effrms++;
          //if (gr->GetN()==((MC) ? 6 : 7)) gr->Fit("line","q");
          if (gr->GetN()==((MC) ? 6 : 7) /*&& line->GetParameter(0) < 10.1*/)
            {
              mg->Add(gr);
            }
        }
    }
  std::cout << "Drawing plots" << std::endl;
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

#endif
