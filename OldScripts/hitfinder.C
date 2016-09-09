{
  TFile * file = TFile::Open("noZS_reco_hist.root","READ");
  TTreeReader reader("robusthits/RobustHitFinder",file);

  TTreeReaderValue<Int_t> run(reader,"run");
  TTreeReaderValue<Int_t> event(reader,"event");
  TTreeReaderValue<UInt_t> c1(reader,"c1");
  TTreeReaderValue<UInt_t> c2(reader,"c2");
  TTreeReaderValue<UInt_t> trignum(reader,"trignum");
  TTreeReaderValue<Float_t> c1x(reader,"c1x");
  TTreeReaderValue<Float_t> c1y(reader,"c1y");
  TTreeReaderValue<Float_t> c1z(reader,"c1z");
  TTreeReaderValue<Float_t> c2x(reader,"c2x");
  TTreeReaderValue<Float_t> c2y(reader,"c2y");
  TTreeReaderValue<Float_t> c2z(reader,"c2z");
  TTreeReaderValue<Int_t> channel(reader,"channel");
  TTreeReaderValue<Int_t> tpc(reader,"tpc");
  TTreeReaderValue<Int_t> signalsize(reader,"signalsize");
  TTreeReaderValue<Float_t> baseline(reader,"baseline");
  TTreeReaderValue<Float_t> rms(reader,"rms");
  TTreeReaderValue<Float_t> baselineFilter(reader,"baselineFilter");
  TTreeReaderValue<Float_t> rmsFilter(reader,"rmsFilter");
  TTreeReaderValue<Float_t> integral(reader,"integral");
  TTreeReaderValue<Float_t> integralFilter(reader,"integralFilter");
  TTreeReaderValue<Float_t> sigmaintegral(reader,"sigmaintegral");
  TTreeReaderValue<Float_t> sigmaintegralFilter(reader,"sigmaintegralFilter");
  TTreeReaderValue<Float_t> amplitude(reader,"amplitude");
  TTreeReaderValue<Float_t> amplitudeFilter(reader,"amplitudeFilter");
  TTreeReaderValue<Float_t> peaktick(reader,"peaktick");
  TTreeReaderValue<Float_t> peaktickFilter(reader,"peaktickFilter");
  TTreeReaderValue<Float_t> peaktime(reader,"peaktime");
  TTreeReaderValue<Float_t> peaktimeFilter(reader,"peaktimeFilter");
  TTreeReaderValue<Int_t> begintick(reader,"begintick");
  TTreeReaderValue<Int_t> endtick(reader,"endtick");
  TTreeReaderValue<Int_t> width(reader,"width");
  TTreeReaderValue<Float_t> hitx(reader,"hitx");
  TTreeReaderValue<Float_t> hity(reader,"hity");
  TTreeReaderValue<Float_t> hitz(reader,"hitz");
  TTreeReaderValue<Float_t> hiterrxlo(reader,"hiterrxlo");
  TTreeReaderValue<Float_t> hiterrxhi(reader,"hiterrxhi");
  TTreeReaderValue<Float_t> hiterrylo(reader,"hiterrylo");
  TTreeReaderValue<Float_t> hiterryhi(reader,"hiterryhi");
  TTreeReaderValue<Float_t> hiterrzlo(reader,"hiterrzlo");
  TTreeReaderValue<Float_t> hiterrzhi(reader,"hiterrzhi");
  TTreeReaderValue<Float_t> perpdist(reader,"perpdist");
  TTreeReaderValue<Float_t> hitt(reader,"hitt");
  TTreeReaderValue<Float_t> driftdist(reader,"driftdist");
  TTreeReaderValue<Bool_t> countercut(reader,"countercut");
  TTreeReaderValue<Float_t> fitconstant(reader,"fitconstant");
  TTreeReaderValue<Float_t> fitconstanterr(reader,"fitconstanterr");
  TTreeReaderValue<Float_t> fitlinear(reader,"fitlinear");
  TTreeReaderValue<Float_t> fitlinearerr(reader,"fitlinearerr");
  TTreeReaderValue<Float_t> fitquadratic(reader,"fitquadratic");
  TTreeReaderValue<Float_t> fitquadraticerr(reader,"fitquadraticerr");
  TTreeReaderValue<Float_t> fitchi2(reader,"fitchi2");
  TTreeReaderValue<Float_t> fitsumsqrresidual(reader,"fitsumsqrresidual");
  TTreeReaderValue<Float_t> fitndf(reader,"fitndf");
  TTreeReaderValue<Float_t> fitmle(reader,"fitmle");
  TTreeReaderValue<Bool_t> fitsuccess(reader,"fitsuccess");
  TTreeReaderValue<Bool_t> fitrealhit(reader,"fitrealhit");
  TTreeReaderValue<Float_t> segmentlength(reader,"segmentlength");

  Int_t nentries = reader.GetEntries(true);
  
  TH1F * snr = new TH1F("snr","Change in signal-to-noise in hit finding",100,0,0);
  
  TH1F * allchan = new TH1F("allchan","All channels",1000,0,4500);
  std::map<Int_t,TH1F*> channelhits;
  while (reader.Next())
    {
      if (channelhits.find(*channel) == channelhits.end())
	{
	  channelhits[*channel] = new TH1F(TString::Format("%d",*channel),TString::Format("Channel %d",*channel),1000,0,*signalsize);
	}
    }

  reader.SetEntry(0);
  
  while (reader.Next())
    {
      allchan->Fill(*hitPeakTime);
      channelhits[*channel]->Fill(*hitPeakTime);
      //if (*event != 2) continue;
      //if (*channel < 1000) continue;
      //if (*rms < 4) continue;
      TGraph * gr = new TGraph();
      TGraph * grcopy = new TGraph();
      TGraph * grfilt = new TGraph();
      Int_t npts = 0, nptscopy = 0;
      std::cout << "run=" << *run << "  event=" << *event << "  channel=" << *channel << "  size=" << *signalSize << "  integral=" << *hitIntegral << "  PeakTime=" << *hitPeakTime << "  width=" << *hitEndTime-*hitStartTime << "  amp=" << *hitPeakAmplitude << std::endl;
      for (size_t i_dig = 0; i_dig < *signalsize; ++i_dig)
	{
	  gr->SetPoint(npts,i_dig,signalVec[i_dig]);
	  grfilt->SetPoint(npts,i_dig,signalfilterVec[i_dig]);
	  if (i_dig < *hitStartTime || i_dig > *hitEndTime)
	    {
	      nptscopy++;
	      grcopy->SetPoint(nptscopy,i_dig,signalVec[i_dig]);
	    }
	  npts++;
	}
      std::cout << "baseline=" << *baseline << "  rms=" << *rms << std::endl;
      std::cout << "filtamp/filtrms=" << (*fhitPeakAmplitude - *fbaseline) / *frms << "   amp/rms=" << (*hitPeakAmplitude - *baseline) / *rms << std::endl;
      snr->Fill(((*fhitPeakAmplitude - *fbaseline)/(*frms))/((*hitPeakAmplitude - *baseline)/(*rms)));
      
      TCanvas * canv = new TCanvas("c","canv",2000,2000);
      canv->Divide(1,2);
      canv->cd(1);
      gr->Draw("al");
      canv->Update();
      TLine * hitbegin = new TLine(*hitStartTime,gPad->GetUymin(),*hitStartTime,gPad->GetUymax());
      hitbegin->SetLineColor(3);
      hitbegin->SetLineWidth(2);
      TLine * hitend = new TLine(*hitEndTime,gPad->GetUymin(),*hitEndTime,gPad->GetUymax());
      hitend->SetLineColor(3);
      hitend->SetLineWidth(2);
      TLine * base = new TLine(0,*baseline,4500,*baseline);
      base->SetLineColor(2);
      base->SetLineWidth(2);
      TLine * rmslow = new TLine(0,*baseline-*rms,4500,*baseline-*rms);
      rmslow->SetLineColor(9);
      rmslow->SetLineWidth(2);
      TLine * rmshi = new TLine(0,*baseline+*rms,4500,*baseline+*rms);
      rmshi->SetLineColor(9);
      rmshi->SetLineWidth(2);
      hitbegin->Draw();
      hitend->Draw();
      rmslow->Draw();
      rmshi->Draw();
      base->Draw();
      canv->cd(2);
      grfilt->Draw("al");
      canv->Update();
      hitbegin->Draw();
      hitend->Draw();
      TLine * fbase = new TLine(0,*fbaseline,4500,*fbaseline);
      fbase->SetLineColor(2);
      fbase->SetLineWidth(2);
      TLine * frmshi = new TLine(0,*fbaseline+*frms,4500,*fbaseline+*frms);
      frmshi->SetLineColor(9);
      frmshi->SetLineWidth(2);
      TLine * frmslow = new TLine(0,*fbaseline-*frms,4500,*fbaseline-*frms);
      frmslow->SetLineColor(9);
      frmslow->SetLineWidth(2);
      fbase->Draw();
      frmshi->Draw();
      frmslow->Draw();
      canv->WaitPrimitive();
      delete canv;
      
    }
  /*
  TCanvas * canv3 = new TCanvas("canv3","canv3",2000,1600);
  canv3->cd();
  allchan->Draw();
  */
  TCanvas * canv4 = new TCanvas("canv4","canv4",2000,1600);
  canv4->cd();
  snr->Draw();
  /*
  for (const auto &c : channelhits)
    {
      TCanvas * canv2 = new TCanvas("canv2","canv2",2000,1600);
      canv2->cd();
      c.second->Draw();
      canv2->WaitPrimitive();
      delete canv2;
    }
  */
}
