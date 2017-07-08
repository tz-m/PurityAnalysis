void Put(std::map<UInt_t, std::map<UInt_t, Int_t> > & thismap, UInt_t c1, UInt_t c2)
{
  if (thismap.find(c1) == thismap.end())
    {
      thismap[c1][c2] = 1;
    }
  else
    {
      if (thismap[c1].find(c2) == thismap[c1].end())
        {
          thismap[c1][c2] = 1;
        }
      else
        {
          thismap[c1][c2] += 1;
        }
    }
}

void DoTrackEfficiency(const char* filename, const char* suffix)
{
  std::map<UInt_t, std::map<UInt_t, Int_t> > all;
  std::map<UInt_t, std::map<UInt_t, Int_t> > good;

  TH2D * trackeff = new TH2D( "trackeff", "", 10, 27.5, 37.5, 10, 5.5, 15.5 );
  TH2D * trackefferr = new TH2D( "trackefferr", "", 10, 27.5, 37.5, 10, 5.5, 15.5 );

  TFile * file = TFile::Open(filename,"READ");

  TTreeReader allreader("robusthit/AllTree",file);
  TTreeReaderValue<UInt_t> allc1(allreader,"c1");
  TTreeReaderValue<UInt_t> allc2(allreader,"c2");
  TTreeReaderValue<UInt_t> alltrignum(allreader,"trignum");
  TTreeReaderValue<Int_t> allrun(allreader,"run");
  TTreeReaderValue<Int_t> allsubrun(allreader,"subrun");
  TTreeReaderValue<Int_t> allevent(allreader,"event");

  while (allreader.Next())
    {
      Put(all,*allc1,*allc2);
    }

  TTreeReader goodreader("robusthit/GoodTree",file);
  TTreeReaderValue<UInt_t> goodc1(goodreader,"c1");
  TTreeReaderValue<UInt_t> goodc2(goodreader,"c2");
  TTreeReaderValue<UInt_t> goodtrignum(goodreader,"trignum");
  TTreeReaderValue<Int_t> goodrun(goodreader,"run");
  TTreeReaderValue<Int_t> goodsubrun(goodreader,"subrun");
  TTreeReaderValue<Int_t> goodevent(goodreader,"event");

  while (goodreader.Next())
    {
      Put(good,*goodc1,*goodc2);
    }

  for (auto all1 : all)
    {
      for (auto all2 : all1.second)
        {
          UInt_t thisc1 = all1.first;
          UInt_t thisc2 = all2.first;
          Double_t thisall = all2.second;
          Double_t thisgood = good[thisc1][thisc2];
          Double_t thiseff = thisgood/thisall;
          Double_t thisefferr = thiseff*sqrt((1/thisgood)+(1/thisall));
          trackeff->Fill(thisc1,thisc2,thiseff);
          trackefferr->Fill(thisc1,thisc2,thisefferr);
        }
    }

  gStyle->SetPaintTextFormat( ".3f" );

  TCanvas * canvtrackeff = new TCanvas( "canvtrackeff", "", 2400, 1600 );
  trackeff->SetContour( 100 );
  trackeff->SetMinimum( 0.0 );
  trackeff->SetMaximum( 1.0 );
  trackeff->SetNdivisions( 10, "xy" );
  trackeff->Draw( "colz" );
  TH2D * trackefftext = ( TH2D* )trackeff->Clone( "trackefftext" );
  trackefftext->SetBarOffset( 0.2 );
  trackefftext->Draw( "same text" );
  trackefferr->SetBarOffset( -0.2 );
  trackefferr->SetMarkerColor( kRed );
  trackefferr->Draw( "same text" );
  trackeff->SetStats( false );
  trackeff->SetXTitle( "Counter ID -- West" );
  trackeff->SetYTitle( "Counter ID -- East" );
  trackeff->SetZTitle( "Track Finding Efficiency" );
  canvtrackeff->SetTopMargin( 0.04 );
  canvtrackeff->SetRightMargin( 0.14 );
  canvtrackeff->Update();
  canvtrackeff->SaveAs( TString::Format("/home/mthiesse/PurityAnalysis/Scripts/png/canvtrackeff_%s.png",suffix) );
  delete canvtrackeff;
  delete trackefftext;

}

void TrackEfficiency()
{
  DoTrackEfficiency("/media/mthiesse/Dell Portable Hard Drive/PurityData/robust_oldgain_3500us_mcscale1.0_hist.root","sim");
  DoTrackEfficiency("/home/mthiesse/PurityAnalysis/DataFiles/robustreco_prepostoutage_oldgain_hist.root","data");
}
