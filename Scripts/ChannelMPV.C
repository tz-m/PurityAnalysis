void ChannelMPV()
{
  TFile * datafile = TFile::Open("/home/mthiesse/PurityAnalysis/PurityAnalysis_data_allruns_newchanmap_origgains.root","READ");
  TTreeReader datareader("hits",datafile);
  TTreeReaderValue<int> datachannel(datareader,"channel");
  TTreeReaderValue<float> datadqdx(datareader,"dqdx");

  std::vector<int> skippeddatachannels;
  std::vector<int> skippedsimchannels;

  std::map<int,TH1F*> datahists;
  while (datareader.Next())
    {
      if (datahists.find(*datachannel) == datahists.end())
        {
          datahists[*datachannel] = new TH1F(TString::Format("ch%i",*datachannel),TString::Format("Channel %i",*datachannel),100,-200,10000);
          datahists[*datachannel]->GetXaxis()->SetTitle("dQ/dx (ADC/cm)");
        }
      datahists[*datachannel]->Fill(*datadqdx);
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
      hist.second->Fit("gaus","R");

      channelMPVsData[hist.first] = std::make_pair(gaus->GetParameter(1),gaus->GetParameter(2));
      channelMPVsErrData[hist.first] = std::make_pair(gaus->GetParError(1),gaus->GetParError(2));

      //TCanvas * canv = new TCanvas("c","",2000,1600);
      //hist.second->Draw();
      //canv->Update();
      //cin.get();

    }

  delete datafile;

  ///////////////////////////////////////////

  TFile * simfile = TFile::Open("/home/mthiesse/PurityAnalysis/PurityAnalysis_newgain_3ms_mcscale1.0.root","READ");
  TTreeReader simreader("hits",simfile);
  TTreeReaderValue<int> simchannel(simreader,"channel");
  TTreeReaderValue<float> simdqdx(simreader,"dqdx");

  std::map<int,TH1F*> simhists;
  while (simreader.Next())
    {
      if (simhists.find(*simchannel) == simhists.end())
        {
          simhists[*simchannel] = new TH1F(TString::Format("ch%i",*simchannel),TString::Format("Channel %i",*simchannel),200,-200,10000);
          simhists[*simchannel]->GetXaxis()->SetTitle("dQ/dx (ADC/cm)");
        }
      simhists[*simchannel]->Fill(*simdqdx);
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
      hist.second->Fit("gaus","R");

      channelMPVsSim[hist.first] = std::make_pair(gaus->GetParameter(1),gaus->GetParameter(2));
      channelMPVsErrSim[hist.first] = std::make_pair(gaus->GetParError(1),gaus->GetParError(2));

      //TCanvas * canv = new TCanvas("c","",2000,1600);
      //hist.second->Draw();
      //canv->Update();
      //cin.get();
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
  TH2D * hcomp = new TH2D("chanMPV",";Data MPV dQ/dx;Sim MPV dQ/dx",40,1000,4000,40,1000,4000);
  TGraphErrors * changaingraph = new TGraphErrors();
  int nchan = 0;
  for (int chan = 0; chan < 2500; ++chan)
    {
      if (channelMPVsData.find(chan) != channelMPVsData.end() && channelMPVsSim.find(chan) != channelMPVsSim.end())
        {
          hcomp->Fill(channelMPVsData[chan].first,channelMPVsSim[chan].first);
          double gain = channelMPVsData[chan].first / channelMPVsSim[chan].first;
          changaingraph->SetPoint(nchan,chan,gain);
          double dataerr = TMath::Power(channelMPVsErrData[chan].first / channelMPVsData[chan].first,2);
          double simerr = TMath::Power(channelMPVsErrSim[chan].first / channelMPVsSim[chan].first,2);
          double gainerr = sqrt(dataerr + simerr) * gain;
          changaingraph->SetPointError(nchan,0,gainerr);
          ++nchan;

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

  TCanvas * canvhdata = new TCanvas("canvhdata","",2000,1600);
  hdata->Draw();
  canvhdata->Update();
  canvhdata->Write();
  canvhdata->SaveAs("canvhdata.png");

  TCanvas * canvhsim = new TCanvas("canvhsim","",2000,1600);
  hsim->Draw();
  canvhsim->Update();
  canvhsim->Write();
  canvhsim->SaveAs("canvhsim.png");

  TCanvas * canvhcomp = new TCanvas("canvhcomp","",2000,1600);
  hcomp->Draw("colz2");
  canvhcomp->Update();
  canvhcomp->Write();
  canvhcomp->SaveAs("canvhcomp.png");

  TCanvas * canvchangain = new TCanvas("canvchangain","",2000,1600);
  changaingraph->SetMarkerStyle(20);
  changaingraph->SetMarkerSize(2);
  changaingraph->Draw("ape");
  canvchangain->Update();
  canvchangain->Write();
  canvchangain->SaveAs("canvchangain.png");

  for (auto c : skippedsimchannels)
    {
      std::cout << "Skipped SimChannel " << c << std::endl;
    }
  for (auto c : skippeddatachannels)
    {
      std::cout << "Skipped DataChannel " << c << std::endl;
    }
  outtree->Write();
  outfile->Close();
}
