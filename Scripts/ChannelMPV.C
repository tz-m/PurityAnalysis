#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TMultiGraph.h"
#include <map>

void Put(std::map<int,int> & mymap, int chan)
{
  if (mymap.find(chan) == mymap.end())
    {
      mymap[chan] = 1;
    }
  else
    {
      mymap[chan] += 1;
    }
}

void ChannelMPV()
{
  TFile * datafile = TFile::Open("/home/mthiesse/PurityAnalysis/PurityAnalysis_data.root","READ");
  TTreeReader datareader("hits",datafile);
  TTreeReaderValue<int> datachannel(datareader,"channel");
  TTreeReaderValue<float> datadqdx(datareader,"dqdx");
  TTreeReaderValue<bool> dataassumedhit(datareader,"assumedhit");
  TTreeReaderValue<int> databinnum(datareader,"binnum");

  std::vector<int> skippeddatachannels;
  std::vector<int> skippedsimchannels;

  std::map<int,int> numtot_data;
  std::map<int,int> numass_data;

  std::map<int,TH1F*> datahists;
  while (datareader.Next())
    {
      if (*databinnum > 3) continue;
      int chan = *datachannel;
      if (datahists.find(chan) == datahists.end())
        {
          datahists[chan] = new TH1F(TString::Format("ch%i",chan),TString::Format("Channel %i",chan),100,-200,10000);
          datahists[chan]->GetXaxis()->SetTitle("dQ/dx [ADC/cm]");
        }
      datahists[chan]->Fill(*datadqdx);

      Put(numtot_data,chan);
      if (*dataassumedhit) Put(numass_data,chan);
    }

  std::map<int,std::pair<double,double> > channelMPVsData;
  std::map<int,std::pair<double,double> > channelMPVsErrData;
  for (auto hist : datahists)
    {
      if (hist.second->GetEntries()<50)
        {
          skippeddatachannels.push_back(hist.first);
          continue;
        }
      double minfactor = 0.5;
      double maxfactor = 0.6;
      if (hist.second->GetEntries()<1500)
        {
          minfactor = 0.3;
          maxfactor = 0.4;
          hist.second->Rebin();
        }
      if (hist.second->GetEntries()<300)
        {
          minfactor = 0.1;
          maxfactor = 0.1;
          hist.second->Rebin();
        }

      float maxval = hist.second->GetBinContent(hist.second->GetMaximumBin());
      float min = hist.second->GetBinCenter(hist.second->FindFirstBinAbove(minfactor*maxval));
      float max = hist.second->GetBinCenter(hist.second->FindLastBinAbove(maxfactor*maxval));
      TF1 * gaus = new TF1("gaus","gaus",min,max);
      gaus->SetParLimits(1,-1000,10000);
      gaus->SetParameter(1,hist.second->GetBinCenter(hist.second->GetMaximumBin()));
      hist.second->Fit("gaus","RQ");

      channelMPVsData[hist.first] = std::make_pair(gaus->GetParameter(1),gaus->GetParameter(2));
      channelMPVsErrData[hist.first] = std::make_pair(gaus->GetParError(1),gaus->GetParError(2));

      if (hist.first == 462)
        {
          gStyle->SetOptFit(1);
          TCanvas * canv = new TCanvas("c","",2000,1600);
          hist.second->Draw();
          hist.second->SetTitle("");
          canv->Update();
          canv->SaveAs("Channel462Data.png");
        }
    }

  delete datafile;

  ///////////////////////////////////////////

  bool oldgain = true;

  TFile * simfile = TFile::Open("/home/mthiesse/PurityAnalysis/PurityAnalysis_4500us_mcscale1.0.root","READ");
  TTreeReader simreader("hits",simfile);
  TTreeReaderValue<int> simchannel(simreader,"channel");
  TTreeReaderValue<float> simdqdx(simreader,"dqdx");
  TTreeReaderValue<bool> simassumedhit(simreader,"assumedhit");
  TTreeReaderValue<int> simbinnum(simreader,"binnum");

  std::map<int,int> numtot_sim;
  std::map<int,int> numass_sim;

  std::map<int,TH1F*> simhists;
  while (simreader.Next())
    {
      if (*simbinnum > 3) continue;
      int chan = *simchannel;
      if (simhists.find(chan) == simhists.end())
        {
          simhists[chan] = new TH1F(TString::Format("ch%i",chan),TString::Format("Channel %i",chan),200,-200,10000);
          simhists[chan]->GetXaxis()->SetTitle("dQ/dx [ADC/cm]");
        }
      simhists[chan]->Fill(*simdqdx);

      Put(numtot_sim,chan);
      if (*simassumedhit) Put(numass_sim,chan);
    }

  std::map<int,std::pair<double,double> > channelMPVsSim;
  std::map<int,std::pair<double,double> > channelMPVsErrSim;
  for (auto hist : simhists)
    {
      if (hist.second->GetEntries()<50)
        {
          skippedsimchannels.push_back(hist.first);
          continue;
        }
      double minfactor = 0.5;
      double maxfactor = 0.6;
      if (hist.second->GetEntries()<1500)
        {
          minfactor = 0.3;
          maxfactor = 0.4;
          hist.second->Rebin();
        }

      if (hist.second->GetEntries()<300)
        {
          minfactor = 0.1;
          maxfactor = 0.1;
          hist.second->Rebin();
        }


      float maxval = hist.second->GetBinContent(hist.second->GetMaximumBin());
      float min = hist.second->GetBinCenter(hist.second->FindFirstBinAbove(minfactor*maxval));
      float max = hist.second->GetBinCenter(hist.second->FindLastBinAbove(maxfactor*maxval));
      TF1 * gaus = new TF1("gaus","gaus",min,max);
      gaus->SetParLimits(1,-1000,10000);
      gaus->SetParameter(1,hist.second->GetBinCenter(hist.second->GetMaximumBin()));
      hist.second->Fit("gaus","RQ");

      channelMPVsSim[hist.first] = std::make_pair(gaus->GetParameter(1),gaus->GetParameter(2));
      channelMPVsErrSim[hist.first] = std::make_pair(gaus->GetParError(1),gaus->GetParError(2));


      if (hist.first == 462)
        {
          gStyle->SetOptFit(1);
          TCanvas * canv = new TCanvas("c","",2000,1600);
          hist.second->Draw();
          hist.second->SetTitle("");
          canv->Update();
          canv->SaveAs(TString::Format("Channel462Sim_%sgain.png",oldgain ? "old" : "new"));
        }

    }

  delete simfile;

  /////////////////////////////////

  TFile * outfile = TFile::Open("ChannelGains.root","RECREATE");
  TTree * outtree = new TTree("channelgains","Relative gains of channels with signal");
  double datampv;
  double datampverr;
  double datawidth;
  double datawidtherr;
  double simmpv;
  double simmpverr;
  double simwidth;
  double simwidtherr;
  int channel;
  double changain;
  double changainerr;
  outtree->Branch("datampv",&datampv,"datampv/D");
  outtree->Branch("datampverr",&datampverr,"datampverr/D");
  outtree->Branch("datawidth",&datawidth,"datawidth/D");
  outtree->Branch("datawidtherr",&datawidtherr,"datawidtherr/D");
  outtree->Branch("simmpv",&simmpv,"simmpv/D");
  outtree->Branch("simmpverr",&simmpverr,"simmpverr/D");
  outtree->Branch("simwidth",&simwidth,"simwidth/D");
  outtree->Branch("simwidtherr",&simwidtherr,"simwidtherr/D");
  outtree->Branch("channel",&channel,"channel/I");
  outtree->Branch("changain",&changain,"changain/D");
  outtree->Branch("changainerr",&changainerr,"changainerr/D");

  TH1D * hdata = new TH1D("hdata","",100,1000,4000);
  hdata->GetXaxis()->SetTitle("MPV dQ/dx");
  hdata->GetYaxis()->SetTitle("Channels");
  TH1D * hsim = new TH1D("hsim","",100,1000,4000);
  hsim->GetXaxis()->SetTitle("MPV dQ/dx");
  hsim->GetYaxis()->SetTitle("Channels");
  TH2D * hcomp = new TH2D("chanMPV",";Data MPV dQ/dx;Sim MPV dQ/dx",20,1000,4000,20,1000,4000);
  TH1D * hgain = new TH1D("hgain",";Data MPV dQ/dx / Sim MPV dQ/dx;",100,0.6,3);
  TGraphErrors * changaingraph = new TGraphErrors();
  int nchan = 0;
  TGraphErrors * gr_ass = new TGraphErrors();
  int nchanass = 0;
  std::vector<double> gains;
  std::vector<double> gainerrs;
  std::vector<double> invgainerrs;
  for (int chan = 0; chan < 2500; ++chan)
    {
      if (numtot_data.find(chan) != numtot_data.end() && numtot_sim.find(chan) != numtot_sim.end())
        {
          int ass_data = 0;
          if (numass_data.find(chan) != numass_data.end())
            {
              ass_data = numass_data[chan];
            }
          double frac_data = (double)ass_data / (double)numtot_data[chan];
          double fracerr_data = (ass_data==0) ? 0 : frac_data*sqrt((1/((double)ass_data))+(1/((double)numtot_data[chan])));

          int ass_sim = 0;
          if (numass_sim.find(chan) != numass_sim.end())
            {
              ass_sim = numass_sim[chan];
            }
          double frac_sim = (double)ass_sim / (double)numtot_sim[chan];
          double fracerr_sim = (ass_sim==0) ? 0 : frac_sim*sqrt((1/((double)ass_sim))+(1/((double)numtot_sim[chan])));

          gr_ass->SetPoint(nchanass,frac_sim,frac_data);
          gr_ass->SetPointError(nchanass,fracerr_sim,fracerr_data);
          ++nchanass;
        }
      if (channelMPVsData.find(chan) != channelMPVsData.end() && channelMPVsSim.find(chan) != channelMPVsSim.end())
        {
          double gain = channelMPVsData[chan].first / channelMPVsSim[chan].first;
          double dataerr = TMath::Power(channelMPVsErrData[chan].first / channelMPVsData[chan].first,2);
          double simerr = TMath::Power(channelMPVsErrSim[chan].first / channelMPVsSim[chan].first,2);
          double gainerr = sqrt(dataerr + simerr) * gain;
          if (gain<0 || gainerr<0) continue;
          if (channelMPVsErrData[chan].first > 200 || channelMPVsErrData[chan].second > 400 || channelMPVsErrSim[chan].first > 200 || channelMPVsErrSim[chan].second > 400) continue;

          hcomp->Fill(channelMPVsData[chan].first,channelMPVsSim[chan].first);
          changaingraph->SetPoint(nchan,chan,gain);
          changaingraph->SetPointError(nchan,0,gainerr);
          ++nchan;

          hgain->Fill(gain,1/gainerr);

          gains.push_back(gain);
          gainerrs.push_back(gainerr);
          invgainerrs.push_back(1/gainerr);

          datampv = channelMPVsData[chan].first;
          datampverr = channelMPVsErrData[chan].first;
          datawidth = channelMPVsData[chan].second;
          datawidtherr = channelMPVsErrData[chan].second;
          simmpv = channelMPVsSim[chan].first;
          simmpverr = channelMPVsErrSim[chan].first;
          simwidth = channelMPVsSim[chan].second;
          simwidtherr = channelMPVsErrSim[chan].second;
          channel = chan;
          changain = gain;
          changainerr = gainerr;
          outtree->Fill();
        }
      if (channelMPVsData.find(chan) != channelMPVsData.end()) hdata->Fill(channelMPVsData[chan].first);
      if (channelMPVsSim.find(chan) != channelMPVsSim.end()) hsim->Fill(channelMPVsSim[chan].first);
    }

  double mediangain = TMath::Median(gains.size(),gains.data(),invgainerrs.data());
  std::cout << "median gain (unweighted) = " << TMath::Median(gains.size(),gains.data()) << "     median scale factor (weighted 1/gainerr) = " << TMath::Median(gains.size(),gains.data(),invgainerrs.data()) << "     median gain (weighted gainerr) = " << TMath::Median(gains.size(),gains.data(),gainerrs.data()) << std::endl;
  double mediangainstdev = TMath::StdDev(gains.size(),gains.data(),invgainerrs.data());
  std::cout << "stdev gain (unweighted) = " << TMath::StdDev(gains.begin(),gains.end()) << "     stdev gain (weighted 1/gainerr) = " << TMath::StdDev(gains.begin(),gains.end(),invgainerrs.begin()) << "     stdev gain (weighted gainerr) = " << TMath::StdDev(gains.begin(),gains.end(),gainerrs.begin()) << std::endl;


  TCanvas * canvhdata = new TCanvas("canvhdata","",2000,1600);
  hdata->Draw();
  canvhdata->Update();
  canvhdata->Write();
  canvhdata->SaveAs(TString::Format("canvhdata_%sgain.png",oldgain ? "old" : "new"));

  TCanvas * canvhsim = new TCanvas("canvhsim","",2000,1600);
  hsim->Draw();
  canvhsim->Update();
  canvhsim->Write();
  canvhsim->SaveAs(TString::Format("canvhsim_%sgain.png",oldgain ? "old" : "new"));

  TCanvas * canvhcomp = new TCanvas("canvhcomp","",2000,1600);
  hcomp->Draw("colz2");
  gStyle->SetOptStat(0);
  canvhcomp->Update();
  canvhcomp->Write();
  canvhcomp->SaveAs(TString::Format("canvhcomp_%sgain.png",oldgain ? "old" : "new"));

  TCanvas * canvhgain = new TCanvas("canvhgain","",2000,1600);
  hgain->SetMarkerStyle(20);
  hgain->SetLineWidth(2);
  hgain->SetMarkerSize(2);
  hgain->Draw();
  hgain->GetYaxis()->SetTitleOffset(1.4);
  TPaveText * tpg = new TPaveText(0.4,0.8,0.9,0.9,"brNDC");
  tpg->AddText(TString::Format("Median scale factor (weighted by 1/error) = %.3f",mediangain));
  tpg->AddText(TString::Format("StdDev scale factor (weighted by 1/error) = %.3f",mediangainstdev));
  tpg->Draw();
  canvhgain->Update();
  canvhgain->Write();
  canvhgain->SaveAs(TString::Format("canvhgain_%sgain.png",oldgain ? "old" : "new"));

  TCanvas * canvchangain = new TCanvas("canvchangain","",2000,1600);
  changaingraph->SetMarkerStyle(20);
  changaingraph->SetMarkerSize(2);
  changaingraph->Draw("ape");
  changaingraph->SetTitle(";Channel ID;Data MPV dQ/dx / Sim MPV dQ/dx");
  TPaveText * tp = new TPaveText(0.25,0.8,0.75,0.9,"brNDC");
  tp->AddText(TString::Format("Median scale factor (weighted by 1/error) = %.3f",mediangain));
  tp->AddText(TString::Format("StdDev scale factor (weighted by 1/error) = %.3f",mediangainstdev));
  tp->Draw();
  canvchangain->Update();
  canvchangain->Write();
  canvchangain->SaveAs(TString::Format("canvchangain_%sgain.png",oldgain ? "old" : "new"));

  TCanvas * canvass = new TCanvas("canvass","",2000,1600);
  gr_ass->SetMarkerStyle(20);
  gr_ass->SetMarkerSize(2);
  gr_ass->SetMarkerColor(kBlue);
  gr_ass->Draw("APE");
  gr_ass->GetXaxis()->SetLimits(0,1);
  gr_ass->GetYaxis()->SetRangeUser(0,1);
  gr_ass->SetTitle(";Assumed hit fraction, Simulation;Assumed hit fraction, Data");
  TPaveText * tpass = new TPaveText(0.55,0.35,0.85,0.45,"brNDC");
  tpass->AddText(TString::Format("Correlation coefficient = %.3f",gr_ass->GetCorrelationFactor()));
  tpass->Draw();
  canvass->Update();
  canvass->Write();
  canvass->SaveAs(TString::Format("canvass_%sgain.png",oldgain ? "old" : "new"));

  for (auto c : skippedsimchannels)
    {
      std::cout << "Skipped SimChannel " << c << std::endl;
    }
  for (auto c : skippeddatachannels)
    {
      std::cout << "Skipped DataChannel " << c << std::endl;
    }
  std::cout << "Total # included wires = " << gains.size() << std::endl;
  outtree->Write();
  outfile->Close();
}
