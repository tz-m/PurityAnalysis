Double_t CauchyDens(Double_t *x, Double_t *par)
{
  Double_t pi = TMath::Pi();
  Double_t mean = par[0];
  Double_t fwhm = par[1];

  Double_t arg = x[0]-mean;
  Double_t top = 0.5*fwhm;
  Double_t bot = pi*(arg*arg+top*top);

  Double_t func = top/bot;
  return func;
}

Double_t CauchyPeak(Double_t *x, Double_t *par)
{
  Double_t height = par[2];
  Double_t func = height*CauchyDens(x,par);
  return func;
}

const UInt_t nscale = 16;
const UInt_t npur   = 3;

void PurityEfficiency()
{
  std::ofstream systematics;
  systematics.open( "systematics.txt", ios::out );

  Double_t mcscale[nscale] = {0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};
  //Double_t mcscale[nscale] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0};
  Int_t elife[npur] = {2000,2500,3000};
  //Int_t elife[npur] = {1, 2, 3, 5, 8};

  std::map<Int_t, Double_t> counterx;
  counterx[28] = counterx[6] = -50;
  counterx[29] = counterx[7] = -21;
  counterx[30] = counterx[8] = 15;
  counterx[31] = counterx[9] = 44;
  counterx[32] = counterx[10] = 74.5;
  counterx[33] = counterx[11] = 105.5;
  counterx[34] = counterx[12] = 136;
  counterx[35] = counterx[13] = 166.5;
  counterx[36] = counterx[14] = 196.5;
  counterx[37] = counterx[15] = 225.5;

  std::map<Int_t, TGraphErrors*> purgraphs; // purity v noise
  std::map<Int_t, TGraphErrors*> effgraphs; // efficiency v noise
  std::map<Int_t, TGraphErrors*> cpurgraphs; // charge purity v noise
  std::map<Int_t, TGraphErrors*> ceffgraphs; // charge efficiency v noise
  //std::map<Int_t, TGraphErrors*> cratiographs; // charge ratio v noise

  std::map<Int_t, std::map<UInt_t, TH1D*> > chargeresolution;
  std::map<Int_t, std::map<UInt_t, TH1D*> > chargeresidual;

  // same as dd graphs, except use opposite counter information, rather than binned drift distance
  std::map<Int_t, std::map<UInt_t, TGraphErrors*> > purgraphsopp; // purity v drift distance opposite counters
  std::map<Int_t, std::map<UInt_t, TGraphErrors*> > effgraphsopp; // efficiency v drift distance
  std::map<Int_t, std::map<UInt_t, TGraphErrors*> > cpurgraphsopp; // charge purity v drift distance
  std::map<Int_t, std::map<UInt_t, TGraphErrors*> > ceffgraphsopp; // charge efficiency v drift distance
  //std::map<Int_t, std::map<UInt_t, TGraphErrors*> > cratiographsopp; // charge ratio v drift distance

  std::map<Int_t, std::map<UInt_t, TH2D*> > purhisttrig;
  std::map<Int_t, std::map<UInt_t, TH2D*> > effhisttrig;
  std::map<Int_t, std::map<UInt_t, TH2D*> > cpurhisttrig;
  std::map<Int_t, std::map<UInt_t, TH2D*> > ceffhisttrig;
  //std::map<Int_t, std::map<UInt_t, TH2D*> > ratiohisttrig;

  std::map<Int_t, std::map<UInt_t, TH2D*> > purerrhisttrig;
  std::map<Int_t, std::map<UInt_t, TH2D*> > efferrhisttrig;
  std::map<Int_t, std::map<UInt_t, TH2D*> > cpurerrhisttrig;
  std::map<Int_t, std::map<UInt_t, TH2D*> > cefferrhisttrig;
  //std::map<Int_t, std::map<UInt_t, TH2D*> > ratioerrhisttrig;

  TProfile * compareQ = new TProfile("compareQ",";Q_{MC};(Q_{MC}-Q_{Reco})/Q_{MC}",1000,200,50000,"s");
  TProfile * scaleQ = new TProfile("scaleQ",";Q_{MC};Q_{Reco}/Q_{MC}",1000,200,50000,"s");
  TProfile * effQ = new TProfile("effQ",";Q_{MC};Q_{Reco|MC}/Q_{MC}",1000,200,50000,"s");

  for ( UInt_t ipur = 0; ipur < npur; ++ipur )
    {
      purgraphs[elife[ipur]] = new TGraphErrors();
      effgraphs[elife[ipur]] = new TGraphErrors();
      cpurgraphs[elife[ipur]] = new TGraphErrors();
      ceffgraphs[elife[ipur]] = new TGraphErrors();
      //cratiographs[elife[ipur]] = new TGraphErrors();

      Int_t numactual = 0;
      for ( UInt_t iscale = 0; iscale < nscale; ++iscale )
        {

          purgraphsopp[elife[ipur]][iscale] = new TGraphErrors();
          effgraphsopp[elife[ipur]][iscale] = new TGraphErrors();
          cpurgraphsopp[elife[ipur]][iscale] = new TGraphErrors();
          ceffgraphsopp[elife[ipur]][iscale] = new TGraphErrors();
          //cratiographsopp[elife[ipur]][iscale] = new TGraphErrors();

          chargeresolution[elife[ipur]][iscale] = new TH1D( TString::Format( "chgresol_%u_%.1f", elife[ipur], mcscale[iscale] ), ";(Q_{reco} - Q_{MC}) / Q_{MC}; # Entries", 500, -1, 1 );
          //chargeresolution[elife[ipur]][iscale]->Sumw2();
          chargeresidual[elife[ipur]][iscale] = new TH1D( TString::Format( "chgresid_%u_%.1f", elife[ipur], mcscale[iscale] ), ";Q_{reco}-Q_{MC};# Entries", 500, -1000, 1000 );
          //chargeresidual[elife[ipur]][iscale]->Sumw2();

          purhisttrig[elife[ipur]][iscale] = new TH2D( TString::Format( "purhisttrig_%u_%.1f", elife[ipur], mcscale[iscale] ), "", 10, 27.5, 37.5, 10, 5.5, 15.5 );
          effhisttrig[elife[ipur]][iscale] = new TH2D( TString::Format( "effhisttrig_%u_%.1f", elife[ipur], mcscale[iscale] ), "", 10, 27.5, 37.5, 10, 5.5, 15.5 );
          cpurhisttrig[elife[ipur]][iscale] = new TH2D( TString::Format( "cpurhisttrig_%u_%.1f", elife[ipur], mcscale[iscale] ), "", 10, 27.5, 37.5, 10, 5.5, 15.5 );
          ceffhisttrig[elife[ipur]][iscale] = new TH2D( TString::Format( "ceffhisttrig_%u_%.1f", elife[ipur], mcscale[iscale] ), "", 10, 27.5, 37.5, 10, 5.5, 15.5 );
          ratiohisttrig[elife[ipur]][iscale] = new TH2D( TString::Format( "ratiohisttrig_%u_%.1f", elife[ipur], mcscale[iscale] ), "", 10, 27.5, 37.5, 10, 5.5, 15.5 );

          purerrhisttrig[elife[ipur]][iscale] = new TH2D( TString::Format( "purerrhisttrig_%u_%.1f", elife[ipur], mcscale[iscale] ), "", 10, 27.5, 37.5, 10, 5.5, 15.5 );
          efferrhisttrig[elife[ipur]][iscale] = new TH2D( TString::Format( "efferrhisttrig_%u_%.1f", elife[ipur], mcscale[iscale] ), "", 10, 27.5, 37.5, 10, 5.5, 15.5 );
          cpurerrhisttrig[elife[ipur]][iscale] = new TH2D( TString::Format( "cpurerrhisttrig_%u_%.1f", elife[ipur], mcscale[iscale] ), "", 10, 27.5, 37.5, 10, 5.5, 15.5 );
          cefferrhisttrig[elife[ipur]][iscale] = new TH2D( TString::Format( "cefferrhisttrig_%u_%.1f", elife[ipur], mcscale[iscale] ), "", 10, 27.5, 37.5, 10, 5.5, 15.5 );
          ratioerrhisttrig[elife[ipur]][iscale] = new TH2D( TString::Format( "ratioerrhisttrig_%u_%.1f", elife[ipur], mcscale[iscale] ), "", 10, 27.5, 37.5, 10, 5.5, 15.5 );

          TString filename = TString::Format( "/media/mthiesse/Dell Portable Hard Drive/PurityData/robust_newgain_%uus_mcscale%.1f_hist.root",elife[ipur],mcscale[iscale]);
          //TString filename = TString::Format( "robustmcana_%ums_%.1f.root", elife[ipur], mcscale[iscale] );
          std::cout << "Filename = " << filename << std::endl;

          TFile * file = TFile::Open( filename, "READ" );
          if ( !file || !file->IsOpen() ) continue;

          TTreeReader reader( "robustmcana/mcanaevents", file );
          TTreeReaderValue<UInt_t> c1( reader, "c1" );
          TTreeReaderValue<UInt_t> c2( reader, "c2" );
          TTreeReaderValue<Double_t> purity( reader, "purity" );
          TTreeReaderValue<Double_t> efficiency( reader, "efficiency" );
          TTreeReaderValue<Double_t> chargepurity( reader, "chargepurity" );
          TTreeReaderValue<Double_t> chargeefficiency( reader, "chargeefficiency" );
          TTreeReaderValue<Double_t> chargeratio( reader, "chargeratio" );
          TTreeReaderValue<Int_t> tp(reader,"tp");
          TTreeReaderValue<Int_t> fp(reader,"fp");
          TTreeReaderValue<Int_t> fn(reader,"fn");
          TTreeReaderValue<Double_t> tpc(reader,"tpc");
          TTreeReaderValue<Double_t> fpc(reader,"fpc");
          TTreeReaderValue<Double_t> fnc(reader,"fnc");

          std::vector<Double_t> purvec;
          std::vector<Double_t> effvec;
          std::vector<Double_t> cpurvec;
          std::vector<Double_t> ceffvec;
          std::vector<Double_t> cratiovec;

          std::map<UInt_t, std::map<UInt_t, std::vector<Double_t> > > pur_by_trigger;
          std::map<UInt_t, std::map<UInt_t, std::vector<Double_t> > > eff_by_trigger;
          std::map<UInt_t, std::map<UInt_t, std::vector<Double_t> > > cpur_by_trigger;
          std::map<UInt_t, std::map<UInt_t, std::vector<Double_t> > > ceff_by_trigger;
          std::map<UInt_t, std::map<UInt_t, std::vector<Double_t> > > ratio_by_trigger;

          while ( reader.Next() )
            {
              Double_t tp_val = *tp;
              Double_t fp_val = *fp;
              Double_t fn_val = *fn;
              Double_t tpc_val = *tpc;
              Double_t fpc_val = *fpc;
              Double_t fnc_val = *fnc;
              Double_t eventpurity = tp_val / (tp_val + fp_val);
              Double_t eventefficiency = tp_val / (tp_val + fn_val);
              Double_t eventchgpur = tpc_val / (tpc_val + fpc_val);
              Double_t eventchgeff = tpc_val / (tpc_val + fnc_val);

              if ( *purity>-0.5 && *purity<2.0 )
                {
                  purvec.push_back( *purity );
                  pur_by_trigger[*c1][*c2].push_back( *purity );
                }
              if ( *efficiency>0 && *efficiency<=1 )
                {
                  effvec.push_back( *efficiency );
                  eff_by_trigger[*c1][*c2].push_back( *efficiency );
                }
              if ( *chargepurity>-0.5 && *chargepurity<2 )
                {
                  cpurvec.push_back( *chargepurity );
                  cpur_by_trigger[*c1][*c2].push_back( *chargepurity );
                }
              if ( *chargeefficiency>0 && *chargeefficiency<=1 )
                {
                  ceffvec.push_back( *chargeefficiency );
                  ceff_by_trigger[*c1][*c2].push_back( *chargeefficiency );
                }
              if ( *chargeratio > 0 && *chargeratio < 10 )
                {
                  cratiovec.push_back( *chargeratio );
                  ratio_by_trigger[*c1][*c2].push_back( *chargeratio );
                }
            }

          Double_t meanpur = TMath::Mean( purvec.size(), purvec.data() );
          Double_t errpur = TMath::RMS( purvec.size(), purvec.data() )/sqrt( purvec.size() );
          Double_t meaneff = TMath::Mean( effvec.size(), effvec.data() );
          Double_t erreff = TMath::RMS( effvec.size(), effvec.data() )/sqrt( effvec.size() );
          Double_t meancpur = TMath::Mean( cpurvec.size(), cpurvec.data() );
          Double_t errcpur = TMath::RMS( cpurvec.size(), cpurvec.data() )/sqrt( cpurvec.size() );
          Double_t meanceff = TMath::Mean( ceffvec.size(), ceffvec.data() );
          Double_t errceff = TMath::RMS( ceffvec.size(), ceffvec.data() )/sqrt( ceffvec.size() );
          Double_t meancratio = TMath::Mean( cratiovec.size(), cratiovec.data() );
          Double_t errcratio = TMath::RMS( cratiovec.size(), cratiovec.data() )/sqrt( cratiovec.size() );

          std::cout << "Elife=" << elife[ipur] << "  mcscale=" << mcscale[iscale] << "  numactual=" << numactual << std::endl;
          std::cout << "meanpur=" << meanpur << "  errpur=" << errpur << "  meaneff=" << meaneff << "  erreff=" << erreff << std::endl;
          std::cout << "meancpur=" << meancpur << "  errcpur=" << errcpur << "  meanceff=" << meanceff << "  errceff=" << errceff << std::endl;
          std::cout << "meancratio=" << meancratio << "  errcratio=" << errcratio << std::endl;

          purgraphs[elife[ipur]]->SetPoint( numactual, mcscale[iscale], meanpur );
          purgraphs[elife[ipur]]->SetPointError( numactual, 0, errpur );
          effgraphs[elife[ipur]]->SetPoint( numactual, mcscale[iscale], meaneff );
          effgraphs[elife[ipur]]->SetPointError( numactual, 0, erreff );
          cpurgraphs[elife[ipur]]->SetPoint( numactual, mcscale[iscale], meancpur );
          cpurgraphs[elife[ipur]]->SetPointError( numactual, 0, errcpur );
          ceffgraphs[elife[ipur]]->SetPoint( numactual, mcscale[iscale], meanceff );
          ceffgraphs[elife[ipur]]->SetPointError( numactual, 0, errceff );
          cratiographs[elife[ipur]]->SetPoint( numactual, mcscale[iscale], meancratio );
          cratiographs[elife[ipur]]->SetPointError( numactual, 0, errcratio );

          ++numactual;


          Int_t numopp = 0;
          for ( auto c1it : pur_by_trigger )
            {
              for ( auto c2it : c1it.second )
                {
                  meanpur = TMath::Mean( c2it.second.size(), c2it.second.data() );
                  errpur = TMath::RMS( c2it.second.size(), c2it.second.data() )/sqrt( c2it.second.size() );
                  purhisttrig[elife[ipur]][iscale]->Fill( c1it.first, c2it.first, meanpur );
                  purerrhisttrig[elife[ipur]][iscale]->Fill( c1it.first, c2it.first, errpur );
                  std::cout << "c1=" << c1it.first << " c2=" << c2it.first << " meanpur=" << meanpur << " +/- " << errpur << std::endl;
                  systematics << TString::Format( "%ums %.1fscale %u %u pur %f %f\n", elife[ipur], mcscale[iscale], c1it.first, c2it.first, meanpur, errpur );
                  if ( c1it.first - c2it.first == 22 && c2it.first != 6 && c2it.first != 7 && c2it.first != 15 )
                    {
                      purgraphsopp[elife[ipur]][iscale]->SetPoint( numopp, counterx[c2it.first], meanpur );
                      purgraphsopp[elife[ipur]][iscale]->SetPointError( numopp, 14, errpur );
                      numopp++;
                    }
                }
            }
          numopp = 0;
          for ( auto c1it : eff_by_trigger )
            {
              for ( auto c2it : c1it.second )
                {
                  meaneff = TMath::Mean( c2it.second.size(), c2it.second.data() );
                  erreff = TMath::RMS( c2it.second.size(), c2it.second.data() )/sqrt( c2it.second.size() );
                  effhisttrig[elife[ipur]][iscale]->Fill( c1it.first, c2it.first, meaneff );
                  efferrhisttrig[elife[ipur]][iscale]->Fill( c1it.first, c2it.first, erreff );
                  std::cout << "c1=" << c1it.first << " c2=" << c2it.first << " meaneff=" << meaneff << " +/- " << erreff << std::endl;
                  systematics << TString::Format( "%ums %.1fscale %u %u eff %f %f\n", elife[ipur], mcscale[iscale], c1it.first, c2it.first, meaneff, erreff );
                  if ( c1it.first - c2it.first == 22 && c2it.first != 6 && c2it.first != 7 && c2it.first != 15 )
                    {
                      effgraphsopp[elife[ipur]][iscale]->SetPoint( numopp, counterx[c2it.first], meaneff );
                      effgraphsopp[elife[ipur]][iscale]->SetPointError( numopp, 14, erreff );
                      numopp++;
                    }
                }
            }
          numopp = 0;
          for ( auto c1it : cpur_by_trigger )
            {
              for ( auto c2it : c1it.second )
                {
                  meancpur = TMath::Mean( c2it.second.size(), c2it.second.data() );
                  errcpur = TMath::RMS( c2it.second.size(), c2it.second.data() )/sqrt( c2it.second.size() );
                  cpurhisttrig[elife[ipur]][iscale]->Fill( c1it.first, c2it.first, meancpur );
                  cpurerrhisttrig[elife[ipur]][iscale]->Fill( c1it.first, c2it.first, errcpur );
                  std::cout << "c1=" << c1it.first << " c2=" << c2it.first << " meancpur=" << meancpur << " +/- " << errcpur << std::endl;
                  systematics << TString::Format( "%ums %.1fscale %u %u cpur %f %f\n", elife[ipur], mcscale[iscale], c1it.first, c2it.first, meancpur, errcpur );
                  if ( c1it.first - c2it.first == 22 && c2it.first != 6 && c2it.first != 7 && c2it.first != 15 )
                    {
                      cpurgraphsopp[elife[ipur]][iscale]->SetPoint( numopp, counterx[c2it.first], meancpur );
                      cpurgraphsopp[elife[ipur]][iscale]->SetPointError( numopp, 14, errcpur );
                      numopp++;
                    }
                }
            }
          numopp = 0;
          for ( auto c1it : ceff_by_trigger )
            {
              for ( auto c2it : c1it.second )
                {
                  meanceff = TMath::Mean( c2it.second.size(), c2it.second.data() );
                  errceff = TMath::RMS( c2it.second.size(), c2it.second.data() )/sqrt( c2it.second.size() );
                  ceffhisttrig[elife[ipur]][iscale]->Fill( c1it.first, c2it.first, meanceff );
                  cefferrhisttrig[elife[ipur]][iscale]->Fill( c1it.first, c2it.first, errceff );
                  std::cout << "c1=" << c1it.first << " c2=" << c2it.first << " meanceff=" << meanceff << " +/- " << errceff << std::endl;
                  systematics << TString::Format( "%ums %.1fscale %u %u ceff %f %f\n", elife[ipur], mcscale[iscale], c1it.first, c2it.first, meanceff, errceff );
                  if ( c1it.first - c2it.first == 22 && c2it.first != 6 && c2it.first != 7 && c2it.first != 15 )
                    {
                      ceffgraphsopp[elife[ipur]][iscale]->SetPoint( numopp, counterx[c2it.first], meanceff );
                      ceffgraphsopp[elife[ipur]][iscale]->SetPointError( numopp, 14, errceff );
                      numopp++;
                    }
                }
            }
          numopp = 0;
          for ( auto c1it : ratio_by_trigger )
            {
              for ( auto c2it : c1it.second )
                {
                  meancratio = TMath::Mean( c2it.second.size(), c2it.second.data() );
                  errcratio = TMath::RMS( c2it.second.size(), c2it.second.data() )/sqrt( c2it.second.size() );
                  ratiohisttrig[elife[ipur]][iscale]->Fill( c1it.first, c2it.first, meancratio );
                  ratioerrhisttrig[elife[ipur]][iscale]->Fill( c1it.first, c2it.first, errcratio );
                  std::cout << "c1=" << c1it.first << " c2=" << c2it.first << " meancratio=" << meancratio << " +/- " << errcratio << std::endl;
                  systematics << TString::Format( "%ums %.1fscale %u %u cratio %f %f\n", elife[ipur], mcscale[iscale], c1it.first, c2it.first, meancratio, errcratio );
                  if ( c1it.first - c2it.first == 22 && c2it.first != 6 && c2it.first != 7 && c2it.first != 15 )
                    {
                      cratiographsopp[elife[ipur]][iscale]->SetPoint( numopp, counterx[c2it.first], meancratio );
                      cratiographsopp[elife[ipur]][iscale]->SetPointError( numopp, 14, errcratio );
                      numopp++;
                    }
                }
            }

          gStyle->SetPaintTextFormat( ".3f" );

          TCanvas * canvpurhisttrig = new TCanvas( TString::Format( "pur_%ums_%.1f", elife[ipur], mcscale[iscale] ), "pur", 2000, 1600 );
          purhisttrig[elife[ipur]][iscale]->SetContour( 100 );
          purhisttrig[elife[ipur]][iscale]->SetMinimum( 0.0 );
          purhisttrig[elife[ipur]][iscale]->SetMaximum( 1.0 );
          purhisttrig[elife[ipur]][iscale]->SetNdivisions( 10, "xy" );
          purhisttrig[elife[ipur]][iscale]->Draw( "colz" );
          TH2D * purtext = ( TH2D* )purhisttrig[elife[ipur]][iscale]->Clone( "purtext" );
          purtext->SetBarOffset( 0.2 );
          purtext->Draw( "same text" );
          purerrhisttrig[elife[ipur]][iscale]->SetBarOffset( -0.2 );
          purerrhisttrig[elife[ipur]][iscale]->SetMarkerColor( kRed );
          purerrhisttrig[elife[ipur]][iscale]->Draw( "same text" );
          purhisttrig[elife[ipur]][iscale]->SetStats( false );
          purhisttrig[elife[ipur]][iscale]->SetXTitle( "Counter ID -- West" );
          purhisttrig[elife[ipur]][iscale]->SetYTitle( "Counter ID -- East" );
          purhisttrig[elife[ipur]][iscale]->SetZTitle( "Purity TP/( TP+FP )" );
          canvpurhisttrig->SetTopMargin( 0.04 );
          canvpurhisttrig->SetRightMargin( 0.14 );
          canvpurhisttrig->Update();
          canvpurhisttrig->SaveAs( TString::Format( "/home/mthiesse/PurityAnalysis/Scripts/png/canvpurhisttrig_%ums_%.1f.png", elife[ipur], mcscale[iscale] ) );
          delete canvpurhisttrig;
          delete purtext;

          TCanvas * canveffhisttrig = new TCanvas( TString::Format( "eff_%ums_%.1f", elife[ipur], mcscale[iscale] ), "eff", 2000, 1600 );
          effhisttrig[elife[ipur]][iscale]->SetContour( 100 );
          effhisttrig[elife[ipur]][iscale]->SetMinimum( 0.0 );
          effhisttrig[elife[ipur]][iscale]->SetMaximum( 1.0 );
          effhisttrig[elife[ipur]][iscale]->SetNdivisions( 10, "xy" );
          effhisttrig[elife[ipur]][iscale]->Draw( "colz" );
          TH2D * efftext = ( TH2D* )effhisttrig[elife[ipur]][iscale]->Clone( "efftext" );
          efftext->SetBarOffset( 0.2 );
          efftext->Draw( "same text" );
          efferrhisttrig[elife[ipur]][iscale]->SetBarOffset( -0.2 );
          efferrhisttrig[elife[ipur]][iscale]->SetMarkerColor( kRed );
          efferrhisttrig[elife[ipur]][iscale]->Draw( "same text" );
          effhisttrig[elife[ipur]][iscale]->SetStats( false );
          effhisttrig[elife[ipur]][iscale]->SetXTitle( "Counter ID -- West" );
          effhisttrig[elife[ipur]][iscale]->SetYTitle( "Counter ID -- East" );
          effhisttrig[elife[ipur]][iscale]->SetZTitle( "Efficiency TP/( TP+FN )" );
          canveffhisttrig->SetTopMargin( 0.04 );
          canveffhisttrig->SetRightMargin( 0.14 );
          canveffhisttrig->Update();
          canveffhisttrig->SaveAs( TString::Format( "/home/mthiesse/PurityAnalysis/Scripts/png/canveffhisttrig_%ums_%.1f.png", elife[ipur], mcscale[iscale] ) );
          delete canveffhisttrig;
          delete efftext;

          TCanvas * canvcpurhisttrig = new TCanvas( TString::Format( "cpur_%ums_%.1f", elife[ipur], mcscale[iscale] ), "cpur", 2000, 1600 );
          cpurhisttrig[elife[ipur]][iscale]->SetContour( 100 );
          cpurhisttrig[elife[ipur]][iscale]->SetMinimum( 0.0 );
          cpurhisttrig[elife[ipur]][iscale]->SetMaximum( 1.0 );
          cpurhisttrig[elife[ipur]][iscale]->SetNdivisions( 10, "xy" );
          cpurhisttrig[elife[ipur]][iscale]->Draw( "colz" );
          TH2D * cpurtext = ( TH2D* )cpurhisttrig[elife[ipur]][iscale]->Clone( "cpurtext" );
          cpurtext->SetBarOffset( 0.2 );
          cpurtext->Draw( "same text" );
          cpurerrhisttrig[elife[ipur]][iscale]->SetBarOffset( -0.2 );
          cpurerrhisttrig[elife[ipur]][iscale]->SetMarkerColor( kRed );
          cpurerrhisttrig[elife[ipur]][iscale]->Draw( "same text" );
          cpurhisttrig[elife[ipur]][iscale]->SetStats( false );
          cpurhisttrig[elife[ipur]][iscale]->SetXTitle( "Counter ID -- West" );
          cpurhisttrig[elife[ipur]][iscale]->SetYTitle( "Counter ID -- East" );
          cpurhisttrig[elife[ipur]][iscale]->SetZTitle( "Charge Purity ChgTP/( ChgTP+ChgFP )" );
          canvcpurhisttrig->SetTopMargin( 0.04 );
          canvcpurhisttrig->SetRightMargin( 0.14 );
          canvcpurhisttrig->Update();
          canvcpurhisttrig->SaveAs( TString::Format( "/home/mthiesse/PurityAnalysis/Scripts/png/canvcpurhisttrig_%ums_%.1f.png", elife[ipur], mcscale[iscale] ) );
          delete canvcpurhisttrig;
          delete cpurtext;

          TCanvas * canvceffhisttrig = new TCanvas( TString::Format( "ceff_%ums_%.1f", elife[ipur], mcscale[iscale] ), "ceff", 2000, 1600 );
          ceffhisttrig[elife[ipur]][iscale]->SetContour( 100 );
          ceffhisttrig[elife[ipur]][iscale]->SetMinimum( 0.0 );
          ceffhisttrig[elife[ipur]][iscale]->SetMaximum( 1.0 );
          ceffhisttrig[elife[ipur]][iscale]->SetNdivisions( 10, "xy" );
          ceffhisttrig[elife[ipur]][iscale]->Draw( "colz" );
          TH2D * cefftext = ( TH2D* )ceffhisttrig[elife[ipur]][iscale]->Clone( "cefftext" );
          cefftext->SetBarOffset( 0.2 );
          cefftext->Draw( "same text" );
          cefferrhisttrig[elife[ipur]][iscale]->SetBarOffset( -0.2 );
          cefferrhisttrig[elife[ipur]][iscale]->SetMarkerColor( kRed );
          cefferrhisttrig[elife[ipur]][iscale]->Draw( "same text" );
          ceffhisttrig[elife[ipur]][iscale]->SetStats( false );
          ceffhisttrig[elife[ipur]][iscale]->SetXTitle( "Counter ID -- West" );
          ceffhisttrig[elife[ipur]][iscale]->SetYTitle( "Counter ID -- East" );
          ceffhisttrig[elife[ipur]][iscale]->SetZTitle( "Charge Efficiency ChgTP/( ChgTP+ChgFN )" );
          canvceffhisttrig->SetTopMargin( 0.04 );
          canvceffhisttrig->SetRightMargin( 0.14 );
          canvceffhisttrig->Update();
          canvceffhisttrig->SaveAs( TString::Format( "/home/mthiesse/PurityAnalysis/Scripts/png/canvceffhisttrig_%ums_%.1f.png", elife[ipur], mcscale[iscale] ) );
          delete canvceffhisttrig;
          delete cefftext;

          TCanvas * canvratiohisttrig = new TCanvas( TString::Format( "ratio_%ums_%.1f", elife[ipur], mcscale[iscale] ), "ratio", 2000, 1600 );
          ratiohisttrig[elife[ipur]][iscale]->SetContour( 100 );
          ratiohisttrig[elife[ipur]][iscale]->SetMinimum( 0.0 );
          ratiohisttrig[elife[ipur]][iscale]->SetMaximum( 3.0 );
          ratiohisttrig[elife[ipur]][iscale]->SetNdivisions( 10, "xy" );
          ratiohisttrig[elife[ipur]][iscale]->Draw( "colz" );
          TH2D * ratiotext = ( TH2D* )ratiohisttrig[elife[ipur]][iscale]->Clone( "ratiotext" );
          ratiotext->SetBarOffset( 0.2 );
          ratiotext->Draw( "same text" );
          ratioerrhisttrig[elife[ipur]][iscale]->SetBarOffset( -0.2 );
          ratioerrhisttrig[elife[ipur]][iscale]->SetMarkerColor( kRed );
          ratioerrhisttrig[elife[ipur]][iscale]->Draw( "same text" );
          ratiohisttrig[elife[ipur]][iscale]->SetStats( false );
          ratiohisttrig[elife[ipur]][iscale]->SetXTitle( "Counter ID -- West" );
          ratiohisttrig[elife[ipur]][iscale]->SetYTitle( "Counter ID -- East" );
          ratiohisttrig[elife[ipur]][iscale]->SetZTitle( "Reco Hit Charge / MC Hit Charge" );
          canvratiohisttrig->SetTopMargin( 0.04 );
          canvratiohisttrig->SetRightMargin( 0.14 );
          canvratiohisttrig->Update();
          canvratiohisttrig->SaveAs( TString::Format( "/home/mthiesse/PurityAnalysis/Scripts/png/canvratiohisttrig_%ums_%.1f.png", elife[ipur], mcscale[iscale] ) );
          delete canvratiohisttrig;
          delete ratiotext;

          //break;
          //continue;

          TTreeReader readerhits( "robustmcana/mcanahits", file );
          TTreeReaderValue<Bool_t> foundBoth( readerhits, "foundBoth" );
          TTreeReaderValue<Double_t> RecoQ( readerhits, "RecoQ" );
          TTreeReaderValue<Double_t> MCQ( readerhits, "MCQ" );

          while (readerhits.Next())
            {
              Bool_t fb = *foundBoth;
              Double_t recoq = *RecoQ;
              Double_t mcq = *MCQ;
              Double_t resid = mcq-recoq;
              Double_t resol = resid/mcq;
              Double_t ratio = recoq/mcq;
              if (fb) chargeresolution[elife[ipur]][iscale]->Fill( resol );
              if (fb) chargeresidual[elife[ipur]][iscale]->Fill( resid );
              if (mcq > 0 && mcq < 100000 && recoq > 0)
                {
                  if (fb && resol>-1 && resol<1) compareQ->Fill( mcq, resol );
                  if (fb) scaleQ->Fill( mcq, ratio );
                  effQ->Fill( mcq, (Double_t)(fb) );
                }
            }

          TF1 *func = new TF1("cauchy",CauchyDens,-1000,1000,2);

          Double_t par[2];
          par[0] = 200; par[1] = 1000;
          func->SetParameters(par);
          func->SetParLimits(0,-1000,1000);
          func->SetParLimits(1,0,3000);
          func->SetParNames("Mean","FWHM");
          func->SetNpx(10000);

          TCanvas * canvchgresol = new TCanvas(TString::Format( "canvchgresol_%ums_%.1f", elife[ipur], mcscale[iscale] ),"",2000,1600 );
          chargeresolution[elife[ipur]][iscale]->Draw();
          Double_t maxbincontent = chargeresolution[elife[ipur]][iscale]->GetBinContent(chargeresolution[elife[ipur]][iscale]->GetMaximumBin());
          Double_t min = chargeresolution[elife[ipur]][iscale]->GetBinCenter(chargeresolution[elife[ipur]][iscale]->FindFirstBinAbove(0.65*maxbincontent));
          Double_t max = chargeresolution[elife[ipur]][iscale]->GetBinCenter(chargeresolution[elife[ipur]][iscale]->FindLastBinAbove(0.65*maxbincontent));
          TF1 * fitresolution = new TF1("fitresolution",CauchyPeak,min,max,3);
          Double_t fitresolpar[3];
          fitresolpar[0] = chargeresolution[elife[ipur]][iscale]->GetBinCenter(chargeresolution[elife[ipur]][iscale]->GetMaximumBin());
          fitresolpar[1] = 0.5;
          fitresolpar[2] = chargeresolution[elife[ipur]][iscale]->GetEntries();
          fitresolution->SetParameters(fitresolpar);
          fitresolution->SetParLimits(0,-1,1);
          fitresolution->SetParLimits(1,0.001,1);
          fitresolution->SetParNames("Mean","FWHM","Height");
          //chargeresolution[elife[ipur]][iscale]->Fit(fitresolution,"R");
          //gStyle->SetOptFit(1);
          canvchgresol->Update();
          canvchgresol->SaveAs( TString::Format( "/home/mthiesse/PurityAnalysis/Scripts/png/canvchgresol_%ums_%.1f.png", elife[ipur], mcscale[iscale] ) );
          delete canvchgresol;

          TCanvas * canvchgresid = new TCanvas(TString::Format( "canvchgresid_%ums_%.1f", elife[ipur], mcscale[iscale] ),"",2000,1600 );
          chargeresidual[elife[ipur]][iscale]->Draw();
          Double_t maxbincontent2 = chargeresidual[elife[ipur]][iscale]->GetBinContent(chargeresidual[elife[ipur]][iscale]->GetMaximumBin());
          Double_t min2 = chargeresidual[elife[ipur]][iscale]->GetBinCenter(chargeresidual[elife[ipur]][iscale]->FindFirstBinAbove(0.65*maxbincontent2));
          Double_t max2 = chargeresidual[elife[ipur]][iscale]->GetBinCenter(chargeresidual[elife[ipur]][iscale]->FindLastBinAbove(0.65*maxbincontent2));
          TF1 * fitresidual = new TF1("fitresidual",CauchyPeak,min2,max2,3);
          Double_t fitresidpar[3];
          fitresidpar[0] = chargeresidual[elife[ipur]][iscale]->GetBinCenter(chargeresidual[elife[ipur]][iscale]->GetMaximumBin());
          fitresidpar[1] = par[1];
          fitresidpar[2] = chargeresidual[elife[ipur]][iscale]->GetEntries();
          fitresidual->SetParameters(fitresidpar);
          fitresidual->SetParLimits(0,-1000,1000);
          fitresidual->SetParLimits(1,0,3000);
          fitresidual->SetParNames("Mean","FWHM","Height");
          chargeresidual[elife[ipur]][iscale]->Fit(fitresidual,"R");
          gStyle->SetOptFit(1);
          canvchgresid->Update();
          canvchgresid->SaveAs( TString::Format( "/home/mthiesse/PurityAnalysis/Scripts/png/canvchgresid_%ums_%.1f.png", elife[ipur], mcscale[iscale] ) );
          delete canvchgresid;

          delete file;
          //break;
        }
      //break;
    }

  compareQ->Print();
  TCanvas * canvcompareQ = new TCanvas( "canvcompareQ", "", 2000, 1600 );
  canvcompareQ->cd();
  compareQ->SetMarkerStyle(20);
  compareQ->Draw();
  canvcompareQ->Update();
  canvcompareQ->SaveAs("/home/mthiesse/PurityAnalysis/Scripts/canvcompareQ.png");
  compareQ->SaveAs("/home/mthiesse/PurityAnalysis/Scripts/compareQ.root");
  //delete canvcompareQ;

  scaleQ->Print();
  TCanvas * canvscaleQ = new TCanvas( "canvscaleQ", "", 2000, 1600 );
  canvscaleQ->cd();
  scaleQ->SetMarkerStyle(20);
  scaleQ->Draw();
  canvscaleQ->Update();
  canvscaleQ->SaveAs("/home/mthiesse/PurityAnalysis/Scripts/canvscaleQ.png");
  scaleQ->SaveAs("/home/mthiesse/PurityAnalysis/Scripts/scaleQ.root");
  //delete canvscaleQ;

  effQ->Print();
  TCanvas * canveffQ = new TCanvas( "canveffQ", "", 2000, 1600 );
  canveffQ->cd();
  effQ->SetMarkerStyle(20);
  effQ->Draw();
  canveffQ->Update();
  canveffQ->SaveAs("/home/mthiesse/PurityAnalysis/Scripts/canveffQ.png");
  effQ->SaveAs("/home/mthiesse/PurityAnalysis/Scripts/effQ.root");
  //delete canveffQ;

  //////////////////////////////////////////////

  TCanvas * purvnoise = new TCanvas( "purvnoise", "purvnoise", 2000, 1600 );
  TMultiGraph * purmg = new TMultiGraph();
  TLegend * purleg = new TLegend( 0.7, 0.1, 0.9, 0.3 );
  for ( auto igr : purgraphs )
    {
      if ( igr.second->GetN()==0 ) continue;
      igr.second->SetLineColor( igr.first );
      igr.second->SetLineWidth( 3 );
      igr.second->SetMarkerStyle( 21 );
      igr.second->SetMarkerSize( 2 );
      igr.second->SetMarkerColor( igr.first );
      purmg->Add( igr.second, "LP" );
      purleg->AddEntry( igr.second, TString::Format( "%ims Lifetime", igr.first ), "lep" );
    }
  purmg->Draw( "A" );
  purmg->GetXaxis()->SetTitle( "Simulated Signal-to-Noise / 35-ton Signal-to-Noise" );
  purmg->GetYaxis()->SetTitle( "Purity, TP/( TP+FP )" );
  purleg->Draw();
  purvnoise->Update();
  purvnoise->SaveAs( "/home/mthiesse/PurityAnalysis/Scripts/png/purvnoise.png" );
  delete purvnoise;

  ///////////////////

  TCanvas * effvnoise = new TCanvas( "effvnoise", "effvnoise", 2000, 1600 );
  TMultiGraph * effmg = new TMultiGraph();
  TLegend * effleg = new TLegend( 0.7, 0.1, 0.9, 0.3 );
  for ( auto igr : effgraphs )
    {
      if ( igr.second->GetN()==0 ) continue;
      igr.second->SetLineColor( igr.first );
      igr.second->SetLineWidth( 3 );
      igr.second->SetMarkerStyle( 21 );
      igr.second->SetMarkerSize( 2 );
      igr.second->SetMarkerColor( igr.first );
      effmg->Add( igr.second, "LP" );
      effleg->AddEntry( igr.second, TString::Format( "%ims Lifetime", igr.first ), "lep" );
    }
  effmg->Draw( "A" );
  effmg->GetXaxis()->SetTitle( "Simulated Signal-to-Noise / 35-ton Signal-to-Noise" );
  effmg->GetYaxis()->SetTitle( "Efficiency, TP/( TP+FN )" );
  effleg->Draw();
  effvnoise->Update();
  effvnoise->SaveAs( "/home/mthiesse/PurityAnalysis/Scripts/png/effvnoise.png" );
  delete effvnoise;

  ////////////////////////////

  TCanvas * cpurvnoise = new TCanvas( "cpurvnoise", "cpurvnoise", 2000, 1600 );
  TMultiGraph * cpurmg = new TMultiGraph();
  TLegend * cpurleg = new TLegend( 0.7, 0.1, 0.9, 0.3 );
  for ( auto igr : cpurgraphs )
    {
      if ( igr.second->GetN()==0 ) continue;
      igr.second->SetLineColor( igr.first );
      igr.second->SetLineWidth( 3 );
      igr.second->SetMarkerStyle( 21 );
      igr.second->SetMarkerSize( 2 );
      igr.second->SetMarkerColor( igr.first );
      cpurmg->Add( igr.second, "LP" );
      cpurleg->AddEntry( igr.second, TString::Format( "%ims Lifetime", igr.first ), "lep" );
    }
  cpurmg->Draw( "A" );
  cpurmg->GetXaxis()->SetTitle( "Simulated Signal-to-Noise / 35-ton Signal-to-Noise" );
  cpurmg->GetYaxis()->SetTitle( "Charge Purity, ChgTP/( ChgTP+ChgFP )" );
  cpurleg->Draw();
  cpurvnoise->Update();
  cpurvnoise->SaveAs( "/home/mthiesse/PurityAnalysis/Scripts/png/cpurvnoise.png" );
  delete cpurvnoise;

  ///////////////////

  TCanvas * ceffvnoise = new TCanvas( "ceffvnoise", "ceffvnoise", 2000, 1600 );
  TMultiGraph * ceffmg = new TMultiGraph();
  TLegend * ceffleg = new TLegend( 0.7, 0.1, 0.9, 0.3 );
  for ( auto igr : ceffgraphs )
    {
      if ( igr.second->GetN()==0 ) continue;
      igr.second->SetLineColor( igr.first );
      igr.second->SetLineWidth( 3 );
      igr.second->SetMarkerStyle( 21 );
      igr.second->SetMarkerSize( 2 );
      igr.second->SetMarkerColor( igr.first );
      ceffmg->Add( igr.second, "LP" );
      ceffleg->AddEntry( igr.second, TString::Format( "%ims Lifetime", igr.first ), "lep" );
    }
  ceffmg->Draw( "A" );
  ceffmg->GetXaxis()->SetTitle( "Simulated Signal-to-Noise / 35-ton Signal-to-Noise" );
  ceffmg->GetYaxis()->SetTitle( "Charge Efficiency, ChgTP/( ChgTP+ChgFN )" );
  ceffleg->Draw();
  ceffvnoise->Update();
  ceffvnoise->SaveAs( "/home/mthiesse/PurityAnalysis/Scripts/png/ceffvnoise.png" );
  delete ceffvnoise;

  ///////////////////////////

  TCanvas * cratiovnoise = new TCanvas( "cratiovnoise", "cratiovnoise", 2000, 1600 );
  TMultiGraph * cratiomg = new TMultiGraph();
  TLegend * cratioleg = new TLegend( 0.7, 0.1, 0.9, 0.3 );
  for ( auto igr : cratiographs )
    {
      if ( igr.second->GetN()==0 ) continue;
      igr.second->SetLineColor( igr.first );
      igr.second->SetLineWidth( 3 );
      igr.second->SetMarkerStyle( 21 );
      igr.second->SetMarkerSize( 2 );
      igr.second->SetMarkerColor( igr.first );
      cratiomg->Add( igr.second, "LP" );
      cratioleg->AddEntry( igr.second, TString::Format( "%ims Lifetime", igr.first ), "lep" );
    }
  cratiomg->Draw( "A" );
  cratiomg->GetXaxis()->SetTitle( "Simulated Signal-to-Noise / 35-ton Signal-to-Noise" );
  cratiomg->GetYaxis()->SetTitle( "Charge Ratio, Q/Q" );
  cratioleg->Draw();
  cratiovnoise->Update();
  cratiovnoise->SaveAs( "/home/mthiesse/PurityAnalysis/Scripts/png/cratiovnoise.png" );
  delete cratiovnoise;

/////////////////////////////



  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////

  std::map<Int_t, TCanvas*> purcanvmap;

  Int_t palette[16];
  Double_t Red[] = {0., 0.0, 1.0, 1.0, 1.0};
  Double_t Green[] = {0., 0.0, 0.0, 1.0, 1.0};
  Double_t Blue[]   = {0., 1.0, 0.0, 0.0, 1.0};
  Double_t Length[] = {0., .25, .50, .75, 1.0};
  Int_t FI = TColor::CreateGradientColorTable( 5, Length, Red, Green, Blue, 16 );
  for ( int i=0; i<16; i++ ) palette[i] = FI+i;

  for ( auto igrp : purgraphsopp )
    {
      TString tit = TString::Format( "purvdrift_%ims", igrp.first );
      purcanvmap[igrp.first] = new TCanvas( tit.Data(), tit.Data(), 2000, 1600 );
      purcanvmap[igrp.first]->cd();
      TLegend * purleg = new TLegend( 0.1, 0.1, 0.25, 0.48 );
      purleg->SetHeader( TString::Format( "eLifetime=%ims", igrp.first ) );
      TIter next( purleg->GetListOfPrimitives() );
      TLegendEntry *first = ( TLegendEntry* )next();
      first->SetTextColor( 2 );
      first->SetTextSize( 0.025 );
      Int_t col = 0;
      Int_t ngraphs = 0;
      TMultiGraph * mg = new TMultiGraph();
      for ( auto igr : igrp.second )
        {
          if ( igr.second->GetN()==0 ) continue;
          igr.second->SetLineColor( palette[col] );
          igr.second->SetLineWidth( 3 );
          if ( fabs( mcscale[igr.first]-1.0 )<0.01 ) { igr.second->SetMarkerStyle( 29 ); igr.second->SetMarkerSize( 4 ); }
          else { igr.second->SetMarkerStyle( 21 ); igr.second->SetMarkerSize( 2 ); }
          igr.second->SetMarkerColor( palette[col] );
          mg->Add( igr.second, "LP" );
          purleg->AddEntry( igr.second, TString::Format( "MCScale=%.1f", mcscale[igr.first] ), "lep" );
          ++ngraphs;
          ++col;
        }
      if ( ngraphs!=0 )
        {
          mg->Draw( "A" );
          mg->GetXaxis()->SetTitle( "Mean Drift Distance ( cm )" );
          mg->GetYaxis()->SetTitle( "Purity, TP/( TP+FP )" );
          purleg->Draw();
          purcanvmap[igrp.first]->Update();
          purcanvmap[igrp.first]->SaveAs( TString::Format( "/home/mthiesse/PurityAnalysis/Scripts/png/%s.png", tit.Data() ) );
          delete purcanvmap[igrp.first];
        }
      else
        {
          delete purcanvmap[igrp.first];
        }
    }

  ///////////////////////////

  std::map<Int_t, TCanvas*> effcanvmap;

  for ( auto igrp : effgraphsopp )
    {
      TString tit = TString::Format( "effvdrift_%ims", igrp.first );
      effcanvmap[igrp.first] = new TCanvas( tit.Data(), tit.Data(), 2000, 1600 );
      effcanvmap[igrp.first]->cd();
      TLegend * effleg = new TLegend( 0.85, 0.6, 0.98, 0.98 );
      effleg->SetHeader( TString::Format( "eLifetime=%ims", igrp.first ) );
      TIter next( effleg->GetListOfPrimitives() );
      TLegendEntry *first = ( TLegendEntry* )next();
      first->SetTextColor( 2 );
      first->SetTextSize( 0.025 );
      Int_t col = 0;
      Int_t ngraphs = 0;
      TMultiGraph * mg = new TMultiGraph();
      for ( auto igr : igrp.second )
        {
          if ( igr.second->GetN()==0 ) continue;
          igr.second->SetLineColor( palette[col] );
          igr.second->SetLineWidth( 3 );
          if ( fabs( mcscale[igr.first]-1.0 )<0.01 ) { igr.second->SetMarkerStyle( 29 ); igr.second->SetMarkerSize( 4 ); }
          else { igr.second->SetMarkerStyle( 21 ); igr.second->SetMarkerSize( 2 ); }
          igr.second->SetMarkerColor( palette[col] );
          mg->Add( igr.second, "LP" );
          effleg->AddEntry( igr.second, TString::Format( "MCScale=%.1f", mcscale[igr.first] ), "lep" );
          ++ngraphs;
          ++col;
        }
      if ( ngraphs!=0 )
        {
          mg->Draw( "A" );
          mg->GetXaxis()->SetTitle( "Mean Drift Distance ( cm )" );
          mg->GetYaxis()->SetTitle( "Efficiency, TP/( TP+FN )" );
          effleg->Draw();
          effcanvmap[igrp.first]->Update();
          effcanvmap[igrp.first]->SaveAs( TString::Format( "/home/mthiesse/PurityAnalysis/Scripts/png/%s.png", tit.Data() ) );
          delete effcanvmap[igrp.first];
        }
      else
        {
          delete effcanvmap[igrp.first];
        }
    }

  //////////////////////////////////

  std::map<Int_t, TCanvas*> cpurcanvmap;

  for ( auto igrp : cpurgraphsopp )
    {
      TString tit = TString::Format( "cpurvdrift_%ims", igrp.first );
      cpurcanvmap[igrp.first] = new TCanvas( tit.Data(), tit.Data(), 2000, 1600 );
      cpurcanvmap[igrp.first]->cd();
      TLegend * purleg = new TLegend( 0.1, 0.1, 0.25, 0.48 );
      purleg->SetHeader( TString::Format( "eLifetime=%ims", igrp.first ) );
      TIter next( purleg->GetListOfPrimitives() );
      TLegendEntry *first = ( TLegendEntry* )next();
      first->SetTextColor( 2 );
      first->SetTextSize( 0.025 );
      Int_t col = 0;
      Int_t ngraphs = 0;
      TMultiGraph * mg = new TMultiGraph();
      for ( auto igr : igrp.second )
        {
          if ( igr.second->GetN()==0 ) continue;
          igr.second->SetLineColor( palette[col] );
          igr.second->SetLineWidth( 3 );
          if ( fabs( mcscale[igr.first]-1.0 )<0.01 ) { igr.second->SetMarkerStyle( 29 ); igr.second->SetMarkerSize( 4 ); }
          else { igr.second->SetMarkerStyle( 21 ); igr.second->SetMarkerSize( 2 ); }
          igr.second->SetMarkerColor( palette[col] );
          mg->Add( igr.second, "LP" );
          purleg->AddEntry( igr.second, TString::Format( "MCScale=%.1f", mcscale[igr.first] ), "lep" );
          ++ngraphs;
          ++col;
        }
      if ( ngraphs!=0 )
        {
          mg->Draw( "A" );
          mg->GetXaxis()->SetTitle( "Mean Drift Distance ( cm )" );
          mg->GetYaxis()->SetTitle( "Charge Purity, ChgTP/( ChgTP+ChgFP )" );
          purleg->Draw();
          cpurcanvmap[igrp.first]->Update();
          cpurcanvmap[igrp.first]->SaveAs( TString::Format( "/home/mthiesse/PurityAnalysis/Scripts/png/%s.png", tit.Data() ) );
          delete cpurcanvmap[igrp.first];
        }
      else
        {
          delete cpurcanvmap[igrp.first];
        }
    }

  ///////////////////////////

  std::map<Int_t, TCanvas*> ceffcanvmap;

  for ( auto igrp : ceffgraphsopp )
    {
      TString tit = TString::Format( "ceffvdrift_%ims", igrp.first );
      ceffcanvmap[igrp.first] = new TCanvas( tit.Data(), tit.Data(), 2000, 1600 );
      ceffcanvmap[igrp.first]->cd();
      TLegend * effleg = new TLegend( 0.85, 0.6, 0.98, 0.98 );
      effleg->SetHeader( TString::Format( "eLifetime=%ims", igrp.first ) );
      TIter next( effleg->GetListOfPrimitives() );
      TLegendEntry *first = ( TLegendEntry* )next();
      first->SetTextColor( 2 );
      first->SetTextSize( 0.025 );
      Int_t col = 0;
      Int_t ngraphs = 0;
      TMultiGraph * mg = new TMultiGraph();
      for ( auto igr : igrp.second )
        {
          if ( igr.second->GetN()==0 ) continue;
          igr.second->SetLineColor( palette[col] );
          igr.second->SetLineWidth( 3 );
          if ( fabs( mcscale[igr.first]-1.0 )<0.01 ) { igr.second->SetMarkerStyle( 29 ); igr.second->SetMarkerSize( 4 ); }
          else { igr.second->SetMarkerStyle( 21 ); igr.second->SetMarkerSize( 2 ); }
          igr.second->SetMarkerColor( palette[col] );
          mg->Add( igr.second, "LP" );
          effleg->AddEntry( igr.second, TString::Format( "MCScale=%.1f", mcscale[igr.first] ), "lep" );
          ++ngraphs;
          ++col;
        }
      if ( ngraphs!=0 )
        {
          mg->Draw( "A" );
          mg->GetXaxis()->SetTitle( "Mean Drift Distance ( cm )" );
          mg->GetYaxis()->SetTitle( "Charge Efficiency, ChgTP/( ChgTP+ChgFN )" );
          effleg->Draw();
          ceffcanvmap[igrp.first]->Update();
          ceffcanvmap[igrp.first]->SaveAs( TString::Format( "/home/mthiesse/PurityAnalysis/Scripts/png/%s.png", tit.Data() ) );
          delete ceffcanvmap[igrp.first];
        }
      else
        {
          delete ceffcanvmap[igrp.first];
        }
    }

  //////////////////////////////////////////

  std::map<Int_t, TCanvas*> ratiocanvmap;

  for ( auto igrp : cratiographsopp )
    {
      TString tit = TString::Format( "ratiovdrift_%ims", igrp.first );
      ratiocanvmap[igrp.first] = new TCanvas( tit.Data(), tit.Data(), 2000, 1600 );
      ratiocanvmap[igrp.first]->cd();
      TLegend * effleg = new TLegend( 0.85, 0.6, 0.98, 0.98 );
      effleg->SetHeader( TString::Format( "eLifetime=%ims", igrp.first ) );
      TIter next( effleg->GetListOfPrimitives() );
      TLegendEntry *first = ( TLegendEntry* )next();
      first->SetTextColor( 2 );
      first->SetTextSize( 0.025 );
      Int_t col = 0;
      Int_t ngraphs = 0;
      TMultiGraph * mg = new TMultiGraph();
      for ( auto igr : igrp.second )
        {
          if ( igr.second->GetN()==0 ) continue;
          igr.second->SetLineColor( palette[col] );
          igr.second->SetLineWidth( 3 );
          if ( fabs( mcscale[igr.first]-1.0 )<0.01 ) { igr.second->SetMarkerStyle( 29 ); igr.second->SetMarkerSize( 4 ); }
          else { igr.second->SetMarkerStyle( 21 ); igr.second->SetMarkerSize( 2 ); }
          igr.second->SetMarkerColor( palette[col] );
          mg->Add( igr.second, "LP" );
          effleg->AddEntry( igr.second, TString::Format( "MCScale=%.1f", mcscale[igr.first] ), "lep" );
          ++ngraphs;
          ++col;
        }
      if ( ngraphs!=0 )
        {
          mg->Draw( "A" );
          mg->GetXaxis()->SetTitle( "Mean Drift Distance ( cm )" );
          mg->GetYaxis()->SetTitle( "Charge Ratio, Qreco/Qsim" );
          effleg->Draw();
          ratiocanvmap[igrp.first]->Update();
          ratiocanvmap[igrp.first]->SaveAs( TString::Format( "/home/mthiesse/PurityAnalysis/Scripts/png/%s.png", tit.Data() ) );
          delete ratiocanvmap[igrp.first];
        }
      else
        {
          delete ratiocanvmap[igrp.first];
        }
    }
  systematics.close();
}
