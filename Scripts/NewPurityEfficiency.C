#include <map>
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TLegendEntry.h"
#include "TAxis.h"
#include "TF1.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

void FillTrigMap(std::map<UInt_t,std::map<UInt_t,Double_t> > & trigmap, UInt_t c1, UInt_t c2, Double_t val)
{
  if (trigmap.find(c1) == trigmap.end())
    {
      trigmap[c1][c2] = val;
    }
  else
    {
      if (trigmap[c1].find(c2) == trigmap[c1].end())
        {
          trigmap[c1][c2] = val;
        }
      else
        {
          trigmap[c1][c2] += val;
        }
    }
}

void PurEffErr(Double_t tp, Double_t fp, Double_t fn, Int_t N, Double_t & pur, Double_t & eff, Double_t & purerr, Double_t & efferr)
{
  pur = tp/(tp+fp);
  eff = tp/(tp+fn);
  if (N != -1)
    {
      purerr = sqrt(pur*(1-pur)/(N));
      efferr = sqrt(eff*(1-eff)/(N));
    }
  else
    {
      purerr = sqrt(pur*(1-pur)/(tp+fp));
      efferr = sqrt(eff*(1-eff)/(tp+fn));
    }
}

const UInt_t nscale = 16;
const UInt_t npur   = 5;

struct effchg {
  Double_t mcq;
  Double_t recoq;
  Bool_t fb;
  Double_t mct;
};

void NewPurityEfficiency()
{
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");

  Double_t mcscale[nscale] = {0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};
  Int_t elife[npur] = {2500,3000,3500,4000,4500};

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

  std::map<Int_t, TEfficiency*> purgraphs; // purity v noise
  std::map<Int_t, TEfficiency*> effgraphs; // efficiency v noise

  std::map<Int_t, std::map<UInt_t, TH1D*> > effhist;
  std::map<Int_t, std::map<UInt_t, TH1D*> > purhist;
  std::vector<Int_t> interesting_elife = {2500,3000,3500,4000,4500}; // select subset of elife[npur]
  std::vector<UInt_t> interesting_mcscale = {0,5,10,15}; // select subset of mcscale[nscale] with items in vector the indices
  //std::vector<Int_t> interesting_elife = {2500};
  //std::vector<UInt_t> interesting_mcscale = {5};

  std::map<Int_t, std::map<UInt_t, TEfficiency*> > purdrift; // purity v drift distance
  std::map<Int_t, std::map<UInt_t, TEfficiency*> > effdrift; // efficiency v drift distance

  std::map<Int_t, TEfficiency*> purdriftlife;
  std::map<Int_t, TEfficiency*> effdriftlife;

  std::map<Int_t, std::map<UInt_t, TEfficiency*> > purhisttrig;
  std::map<Int_t, std::map<UInt_t, TEfficiency*> > effhisttrig;
  std::map<Int_t, std::map<UInt_t, TH2D*> > purerrhisttrig;
  std::map<Int_t, std::map<UInt_t, TH2D*> > efferrhisttrig;

  //std::map<UInt_t, TProfile*> effcharge; // efficiency v MCQ
  //std::map<UInt_t, TProfile*> effcharge_stdev;
  std::map<UInt_t, TEfficiency*> effcharge_correct;
  //std::map<UInt_t, TEfficiency*> effcharge_short;
  //std::map<UInt_t, TEfficiency*> effcharge_long;
  std::map<UInt_t, TH1D*> numeffcharge;
  std::map<UInt_t, std::vector< effchg > > effchargedata;
  std::map<UInt_t, TProfile*> residual; // charge residual v MCQ
  std::map<UInt_t, std::map<Int_t, TH1D*> > resolution;

  std::map<UInt_t, TH1D*> reshist;

  TProfile * mm = new TProfile("mm",";Q_{MC} [ADC];Q_{Reco} [ADC]",200,0,10000,"e");
  TProfile * mm_stdev = new TProfile("mm_stdev","",200,0,10000,"s");

  TH1D * effhistbig = new TH1D("effhistbig","efficiency, Q > 6000 ADC/cm",200,0,1);

  for ( UInt_t ipur = 0; ipur < npur; ++ipur )
    {
      purgraphs[elife[ipur]] = new TEfficiency(TString::Format("purgraphs_%u",elife[ipur]),";something1;efficiency",16,0.45,2.05);
      effgraphs[elife[ipur]] = new TEfficiency(TString::Format("effgraphs_%u",elife[ipur]),";MCscale;Efficiency, #epsilon_{S}",16,0.45,2.05);

      if ((std::find(interesting_elife.begin(),interesting_elife.end(),elife[ipur]) != interesting_elife.end()))
        {
          purdriftlife[elife[ipur]] = new TEfficiency(TString::Format("purdriftlife_%u",elife[ipur]),";Drift Distance [cm];Purity, #phi",7,0,225);
          effdriftlife[elife[ipur]] = new TEfficiency(TString::Format("effdriftlife_%u",elife[ipur]),";Drift Distance [cm];Efficiency, #epsilon_{S}",7,0,225);
        }

      for ( UInt_t iscale = 0; iscale < nscale; ++iscale )
        {
          Bool_t do_this_elife = (std::find(interesting_elife.begin(),interesting_elife.end(),elife[ipur]) != interesting_elife.end());
          Bool_t do_this_mcscale = (std::find(interesting_mcscale.begin(),interesting_mcscale.end(),iscale) != interesting_mcscale.end());
          if (do_this_mcscale)
            {
              if (do_this_elife)
                {
                  effhist[elife[ipur]][iscale] = new TH1D(TString::Format("effhist_%u_%.1f",elife[ipur],mcscale[iscale]),"",200,0,1);
                  purhist[elife[ipur]][iscale] = new TH1D(TString::Format("purhist_%u_%.1f",elife[ipur],mcscale[iscale]),"",200,0,1);

                  purdrift[elife[ipur]][iscale] = new TEfficiency(TString::Format("purdrift_%u_%.1f",elife[ipur],mcscale[iscale]),"",7,0,225);
                  effdrift[elife[ipur]][iscale] = new TEfficiency(TString::Format("effdrift_%u_%.1f",elife[ipur],mcscale[iscale]),"",7,0,225);

                  purhisttrig[elife[ipur]][iscale] = new TEfficiency( TString::Format( "purhisttrig_%u_%.1f", elife[ipur], mcscale[iscale] ), "", 10, 27.5, 37.5, 10, 5.5, 15.5 );
                  effhisttrig[elife[ipur]][iscale] = new TEfficiency( TString::Format( "effhisttrig_%u_%.1f", elife[ipur], mcscale[iscale] ), "", 10, 27.5, 37.5, 10, 5.5, 15.5 );
                  purerrhisttrig[elife[ipur]][iscale] = new TH2D( TString::Format( "purerrhisttrig_%u_%.1f", elife[ipur], mcscale[iscale] ), "", 10, 27.5, 37.5, 10, 5.5, 15.5 );
                  efferrhisttrig[elife[ipur]][iscale] = new TH2D( TString::Format( "efferrhisttrig_%u_%.1f", elife[ipur], mcscale[iscale] ), "", 10, 27.5, 37.5, 10, 5.5, 15.5 );

                }
              /*
                 if (effcharge.find(iscale) == effcharge.end())
                 {
                 effcharge[iscale] = new TProfile(TString::Format("effcharge_%.1f",mcscale[iscale]),";dQ_{MC}/dx [ADC/cm];Efficiency, #epsilon_{S}", 200,0,8000,"e");
                 effcharge_stdev[iscale] = new TProfile(TString::Format("effcharge_stdev_%.1f",mcscale[iscale]),";dQ_{MC}/dx [ADC/cm];Efficiency, #epsilon_{S}",200,0,8000,"s");
                 }
               */
              if (residual.find(iscale) == residual.end())
                {
                  residual[iscale] = new TProfile(TString::Format("residual_%.1f",mcscale[iscale]),";Q_{MC} [ADC];(Q_{Reco}-Q_{MC}) / Q_{MC}",100,0,8000,"e");
                  for (Int_t i=1; i<=100; ++i)
                    {
                      resolution[iscale][i] = new TH1D(TString::Format("resolution_%.1f_bin%i",mcscale[iscale],i),"",500,-2,2);
                    }
                }
              if (reshist.find(iscale) == reshist.end())
                reshist[iscale] = new TH1D(TString::Format("reshist_%.1f",mcscale[iscale]),"",500,-2,2);
            }

          TString filename = TString::Format( "/media/mthiesse/Dell Portable Hard Drive/PurityData/robust_oldgain_%uus_mcscale%.1f_histredo.root",elife[ipur],mcscale[iscale]);
          //TString filename = "/home/mthiesse/Sandbox/anahist_2_1.0.root";
          //TString filename = "/home/mthiesse/Sandbox/anahist_largermlesaccut.root";
          //TString filename = TString::Format("/media/mthiesse/Dell Portable Hard Drive/PurityData/anahist_2_%uus_%.1f.root",elife[ipur],mcscale[iscale]);
          std::cout << "Filename = " << filename << std::endl;

          TFile * file = TFile::Open( filename, "READ" );
          if ( !file || !file->IsOpen() ) continue;

          TTreeReader reader( "robustmcana/mcanaevents", file );
          TTreeReaderValue<UInt_t> c1( reader, "c1" );
          TTreeReaderValue<UInt_t> c2( reader, "c2" );
          TTreeReaderValue<Int_t> tp(reader,"tp");
          TTreeReaderValue<Int_t> fp(reader,"fp");
          TTreeReaderValue<Int_t> fn(reader,"fn");
          TTreeReaderValue<Int_t> totalMCHits(reader,"totalMCHits");
          TTreeReaderValue<Int_t> matchedMCHits(reader,"matchedMCHits");
          TTreeReaderValue<Int_t> totalRecoHits(reader,"totalRecoHits");
          TTreeReaderValue<Int_t> matchedRecoHits(reader,"matchedRecoHits");

          Double_t total_tp = 0, total_fp = 0, total_fn = 0;

          std::map<UInt_t, std::map<UInt_t, Double_t> > totMC_by_trig;
          std::map<UInt_t, std::map<UInt_t, Double_t> > matchMC_by_trig;
          std::map<UInt_t, std::map<UInt_t, Double_t> > totReco_by_trig;
          std::map<UInt_t, std::map<UInt_t, Double_t> > matchReco_by_trig;
          std::map<UInt_t, std::map<UInt_t, Double_t> > nevents_by_trig;

          Int_t N = 0;
          Float_t countTotalMC = 0, countMatchedMC = 0;
          Float_t countTotalReco = 0, countMatchedReco = 0;
          Int_t nevents = 0;
          while ( reader.Next() )
            {
              countTotalMC += *totalMCHits;
              countMatchedMC += *matchedMCHits;
              countTotalReco += *totalRecoHits;
              countMatchedReco += *matchedRecoHits;
              ++nevents;

              if (do_this_elife && do_this_mcscale)
                {
                  FillTrigMap(totMC_by_trig,*c1,*c2,*totalMCHits);
                  FillTrigMap(matchMC_by_trig,*c1,*c2,*matchedMCHits);
                  FillTrigMap(totReco_by_trig,*c1,*c2,*totalRecoHits);
                  FillTrigMap(matchReco_by_trig,*c1,*c2,*matchedRecoHits);
                  FillTrigMap(nevents_by_trig,*c1,*c2,1);
                  if (*matchedRecoHits > 0) purhist[elife[ipur]][iscale]->Fill((double)(*matchedRecoHits) / (double)(*totalRecoHits));
                  if (*matchedMCHits > 0) effhist[elife[ipur]][iscale]->Fill((double)(*matchedMCHits) / (double)(*totalMCHits));
                }

/*
              Double_t tp_val = *tp;
              Double_t fp_val = *fp;
              Double_t fn_val = *fn;

              if (total_tp >= 0 && total_fp >= 0 && total_fn >= 0)
                {
                  total_tp += tp_val;
                  total_fp += fp_val;
                  total_fn += fn_val;
   ++N;
                }

              if (do_this_elife && do_this_mcscale)
                {
                  Double_t pur_val, eff_val, errpur_val, erreff_val;
                  PurEffErr(tp_val,fp_val,fn_val,-1,pur_val,eff_val,errpur_val,erreff_val);

                  if (eff_val > 0) effhist[elife[ipur]][iscale]->Fill(eff_val);
                  if (pur_val > 0) purhist[elife[ipur]][iscale]->Fill(pur_val);

                  FillTrigMap(totMC_by_trig,*c1,*c2,*totalMCHits);
                  FillTrigMap(fp_by_trig,*c1,*c2,fp_val);
                  FillTrigMap(fn_by_trig,*c1,*c2,fn_val);
                }
 */
            }

          Int_t bin = purgraphs[elife[ipur]]->FindFixBin(mcscale[iscale]);
          purgraphs[elife[ipur]]->SetTotalEvents(bin,nevents);
          purgraphs[elife[ipur]]->SetPassedEvents(bin,countMatchedReco*nevents/countTotalReco);
          effgraphs[elife[ipur]]->SetTotalEvents(bin,nevents);
          effgraphs[elife[ipur]]->SetPassedEvents(bin,countMatchedMC*nevents/countTotalMC);

          std::cout << "Nevents= " << nevents << "  bin=" << bin << "  countTotalReco=" << countTotalReco << "  countMatchedReco=" << countMatchedReco << "  countTotalMC=" << countTotalMC << "  countMatchedMC=" << countMatchedMC << std::endl;

          std::cout << "matchedMC(full)=" << countMatchedMC << "  matchedMC(trunc)=" << (countMatchedMC/countTotalMC)*nevents << std::endl;

          if (mcscale[iscale]-1.1<1e-5)
            {
              for (auto c1it : totMC_by_trig )
                {
                  for (auto c2it : c1it.second)
                    {
                      if (matchMC_by_trig[c1it.first].find(c2it.first) == matchMC_by_trig[c1it.first].end()) continue;
                      if (totReco_by_trig[c1it.first].find(c2it.first) == totReco_by_trig[c1it.first].end()) continue;
                      if (matchReco_by_trig[c1it.first].find(c2it.first) == matchReco_by_trig[c1it.first].end()) continue;
                      Double_t totMC_trig = totMC_by_trig[c1it.first][c2it.first];
                      Double_t matchMC_trig = matchMC_by_trig[c1it.first][c2it.first];
                      Double_t totReco_trig = totReco_by_trig[c1it.first][c2it.first];
                      Double_t matchReco_trig = matchReco_by_trig[c1it.first][c2it.first];
                      Double_t nevents_trig = nevents_by_trig[c1it.first][c2it.first];
                      if ( c1it.first - c2it.first == 22 && c2it.first != 6 && c2it.first != 7 && c2it.first != 15 )
                        {
                          Int_t bin = purdriftlife[elife[ipur]]->FindFixBin(counterx[c2it.first]);
                          purdriftlife[elife[ipur]]->SetTotalEvents(bin,nevents_trig);
                          purdriftlife[elife[ipur]]->SetPassedEvents(bin,matchReco_trig*nevents_trig/totReco_trig);
                          effdriftlife[elife[ipur]]->SetTotalEvents(bin,nevents_trig);
                          effdriftlife[elife[ipur]]->SetPassedEvents(bin,matchMC_trig*nevents_trig/totMC_trig);
                        }
                    }
                }
            }


          if (do_this_elife && do_this_mcscale)
            {
              Int_t numopp = 0;
              for ( auto c1it : totMC_by_trig )
                {
                  for ( auto c2it : c1it.second )
                    {
                      if (matchMC_by_trig[c1it.first].find(c2it.first) == matchMC_by_trig[c1it.first].end()) continue;
                      if (totReco_by_trig[c1it.first].find(c2it.first) == totReco_by_trig[c1it.first].end()) continue;
                      if (matchReco_by_trig[c1it.first].find(c2it.first) == matchReco_by_trig[c1it.first].end()) continue;
                      Double_t totMC_trig = totMC_by_trig[c1it.first][c2it.first];
                      Double_t matchMC_trig = matchMC_by_trig[c1it.first][c2it.first];
                      Double_t totReco_trig = totReco_by_trig[c1it.first][c2it.first];
                      Double_t matchReco_trig = matchReco_by_trig[c1it.first][c2it.first];
                      Double_t nevents_trig = nevents_by_trig[c1it.first][c2it.first];
                      Int_t binc = purhisttrig[elife[ipur]][iscale]->FindFixBin(c1it.first,c2it.first);
                      purhisttrig[elife[ipur]][iscale]->SetTotalEvents(binc,nevents_trig);
                      purhisttrig[elife[ipur]][iscale]->SetPassedEvents(binc,matchReco_trig*nevents_trig/totReco_trig);
                      effhisttrig[elife[ipur]][iscale]->SetTotalEvents(binc,nevents_trig);
                      effhisttrig[elife[ipur]][iscale]->SetPassedEvents(binc,matchMC_trig*nevents_trig/totMC_trig);
                      purerrhisttrig[elife[ipur]][iscale]->SetBinContent(binc,sqrt(pow(purhisttrig[elife[ipur]][iscale]->GetEfficiencyErrorLow(binc),2)+pow(purhisttrig[elife[ipur]][iscale]->GetEfficiencyErrorUp(binc),2)));
                      efferrhisttrig[elife[ipur]][iscale]->SetBinContent(binc,sqrt(pow(effhisttrig[elife[ipur]][iscale]->GetEfficiencyErrorLow(binc),2)+pow(effhisttrig[elife[ipur]][iscale]->GetEfficiencyErrorUp(binc),2)));
                      if ( c1it.first - c2it.first == 22 && c2it.first != 6 && c2it.first != 7 && c2it.first != 15 )
                        {
                          Int_t bin = purdriftlife[elife[ipur]]->FindFixBin(counterx[c2it.first]);
                          purdrift[elife[ipur]][iscale]->SetTotalEvents(bin,nevents_trig);
                          purdrift[elife[ipur]][iscale]->SetPassedEvents(bin,matchReco_trig*nevents_trig/totReco_trig);
                          effdrift[elife[ipur]][iscale]->SetTotalEvents(bin,nevents_trig);
                          effdrift[elife[ipur]][iscale]->SetPassedEvents(bin,matchMC_trig*nevents_trig/totMC_trig);
                        }
                    }
                }
/*
              gStyle->SetPaintTextFormat( ".3f" );

              TCanvas * canvpurhisttrig = new TCanvas( TString::Format( "pur_%uus_%.1f", elife[ipur], mcscale[iscale] ), "pur", 2000, 1600 );
              TH2D * thispurhisttrig = (TH2D*)(purhisttrig[elife[ipur]][iscale]->CreateHistogram());
              thispurhisttrig->SetContour( 100 );
              thispurhisttrig->SetMinimum( 0.0 );
              thispurhisttrig->SetMaximum( 1.0 );
              thispurhisttrig->SetNdivisions( 10, "xy" );
              thispurhisttrig->Draw( "colz" );
              TH2D * purtext = ( TH2D* )thispurhisttrig->Clone( "purtext" );
              purtext->SetBarOffset( 0.2 );
              purtext->Draw( "same text" );
              purerrhisttrig[elife[ipur]][iscale]->SetBarOffset( -0.2 );
              purerrhisttrig[elife[ipur]][iscale]->SetMarkerColor( kRed );
              purerrhisttrig[elife[ipur]][iscale]->Draw( "same text" );
              thispurhisttrig->SetStats( false );
              thispurhisttrig->SetXTitle( "Counter ID -- West" );
              thispurhisttrig->SetYTitle( "Counter ID -- East" );
              thispurhisttrig->SetZTitle( "Purity, #phi" );
              canvpurhisttrig->SetTopMargin( 0.04 );
              canvpurhisttrig->SetRightMargin( 0.14 );
              canvpurhisttrig->Update();
              canvpurhisttrig->SaveAs( TString::Format( "/home/mthiesse/PurityAnalysis/Scripts/png/canvpurhisttrig_%uus_%.1f.png", elife[ipur], mcscale[iscale] ) );
              delete canvpurhisttrig;
              delete purtext;

              TCanvas * canveffhisttrig = new TCanvas( TString::Format( "eff_%uus_%.1f", elife[ipur], mcscale[iscale] ), "eff", 2000, 1600 );
              TH2D * thiseffhisttrig = (TH2D*)(effhisttrig[elife[ipur]][iscale]->CreateHistogram());
              thiseffhisttrig->SetContour( 100 );
              thiseffhisttrig->SetMinimum( 0.0 );
              thiseffhisttrig->SetMaximum( 1.0 );
              thiseffhisttrig->SetNdivisions( 10, "xy" );
              thiseffhisttrig->Draw( "colz" );
              TH2D * efftext = ( TH2D* )thiseffhisttrig->Clone( "efftext" );
              efftext->SetBarOffset( 0.2 );
              efftext->Draw( "same text" );
              efferrhisttrig[elife[ipur]][iscale]->SetBarOffset( -0.2 );
              efferrhisttrig[elife[ipur]][iscale]->SetMarkerColor( kRed );
              efferrhisttrig[elife[ipur]][iscale]->Draw( "same text" );
              thiseffhisttrig->SetStats( false );
              thiseffhisttrig->SetXTitle( "Counter ID -- West" );
              thiseffhisttrig->SetYTitle( "Counter ID -- East" );
              thiseffhisttrig->SetZTitle( "Efficiency, #epsilon_{S}" );
              canveffhisttrig->SetTopMargin( 0.04 );
              canveffhisttrig->SetRightMargin( 0.14 );
              canveffhisttrig->Update();
              canveffhisttrig->SaveAs( TString::Format( "/home/mthiesse/PurityAnalysis/Scripts/png/canveffhisttrig_%uus_%.1f.png", elife[ipur], mcscale[iscale] ) );
              delete canveffhisttrig;
              delete efftext;
 */
            }

          if (!do_this_mcscale || !do_this_elife) continue;

          TTreeReader readerhits( "robustmcana/mcanahits", file );
          TTreeReaderValue<Bool_t> foundBoth( readerhits, "foundBoth" );
          TTreeReaderValue<Double_t> RecoQ( readerhits, "RecoQ" );
          TTreeReaderValue<Double_t> RecoT( readerhits, "RecoT" );
          TTreeReaderValue<Double_t> MCQ( readerhits, "MCQ" );
          TTreeReaderValue<Double_t> MCT( readerhits, "MCT");
          TTreeReaderValue<Int_t> channel( readerhits, "channel");
          TTreeReaderValue<UInt_t> c1new( readerhits, "c1");
          TTreeReaderValue<UInt_t> c2new( readerhits, "c2");

          std::vector<std::pair<Double_t,Double_t> > effcharge_map;

          while (readerhits.Next())
            {
              /*
                 if (*channel == 288 || *channel == 399 ||
                 *channel == 400 || *channel == 511 ||
                 *channel == 800 || *channel == 911 ||
                 *channel == 912 || *channel == 1023 ||
                 *channel == 1312 || *channel == 1423 ||
                 *channel == 1424 || *channel == 1535 ||
                 *channel == 1824 || *channel == 1935 ||
                 *channel == 1936 || *channel == 2047) continue;
               */

              if (*channel <= 400 ||
                  (*channel >= 511 && *channel <= 912) ||
                  (*channel >= 1023 && *channel <= 1424) ||
                  (*channel >= 1535 && *channel <= 1936) ||
                  *channel >= 2047) continue;

              //if (*RecoT < 200 || *RecoT > 2000) continue;


              if ( (*c1new == 28 && (*c2new == 6 || *c2new == 7 || *c2new == 8 || *c2new == 9)) ||
                   (*c1new == 29 && (*c2new == 6 || *c2new == 7 || *c2new == 8)) ||
                   (*c1new == 30 && (*c2new == 6 || *c2new == 7)) ||
                   (*c1new == 31 && (*c2new == 6)) ||
                   (*c1new == 36 && (*c2new == 15)) ||
                   (*c1new == 37 && (*c2new == 15 || *c2new == 14))) continue;


              Bool_t fb = *foundBoth;
              Double_t recoq = *RecoQ / mcscale[iscale];
              //if (recoq < 0) continue;
              Double_t mcq = *MCQ / mcscale[iscale];
              Double_t resid = recoq-mcq;
              Double_t resol = resid/mcq;
              //effcharge[iscale]->Fill( mcq, (Double_t)fb );
              //effcharge_stdev[iscale]->Fill(mcq, (Double_t)fb);
              effchg e;
              e.mcq = mcq;
              e.recoq = recoq;
              e.mct = *MCT;
              e.fb = fb;
              effchargedata[iscale].push_back(e);
              //if (mcq>0 && mcq<8000) effcharge_map.push_back(std::make_pair(mcq,(Double_t)fb));
              if (fb && mcq>0)
                {
                  if (resol>-2 && resol<2)
                    {
                      residual[iscale]->Fill(mcq,resol);
                      Int_t bin = residual[iscale]->FindBin(mcq);
                      if (bin >=1 && bin <= 100) resolution[iscale][bin]->Fill(resol);
                    }
                  mm->Fill(mcq,recoq);
                  mm_stdev->Fill(mcq,recoq);
                  if (mcq>1000 && mcq<1500) reshist[iscale]->Fill(resol);
                }
            }
/*
          std::sort(effcharge_map.begin(),effcharge_map.end());
          std::vector<Double_t> xbins;
          Int_t numinbin = effcharge_map.size()/400;
          for (UInt_t i=0; i<effcharge_map.size(); ++i)
            {
              if (i%numinbin==0) xbins.push_back(effcharge_map[i].first);
            }
          std::cout << "xbins.size=" << xbins.size() << std::endl;
          std::cout << "[";
          for (auto x : xbins)
            {
              std::cout << x << ",";
            }
          std::cout << "]" << std::flush << std::endl;
          return;
 */
          delete file;
        }
    }

//////////////////////////////////////////////
  TGaxis::SetMaxDigits(3);
//////////////////////////////////////////////

  TCanvas * purvnoise = new TCanvas( "purvnoise", "purvnoise", 2000, 1600 );
  purvnoise->SetMargin(0.12,0.02,0.11,0.02); //l,r,b,t
  TMultiGraph * purmg = new TMultiGraph();
  TLegend * purleg = new TLegend( 0.55, 0.15, 0.9, 0.5 );
  Int_t color = 53;
  for ( auto igr : purgraphs )
    {
      TGraphAsymmErrors * purgraph = igr.second->CreateGraph();
      //if ( purgraph->GetN()==0 ) continue;
      purgraph->SetLineColor( color );
      purgraph->SetLineWidth(4);
      purgraph->SetMarkerStyle( 21 );
      purgraph->SetMarkerSize( 2 );
      purgraph->SetMarkerColor( color );
      if (purgraph->GetN() != 0) purmg->Add( purgraph );
      purleg->AddEntry( purgraph, TString::Format( "%i#mus Lifetime", igr.first ), "lep" );
      color+=10;
    }
  purmg->Draw( "alpe" );
  purmg->SetTitle( ";MCscale;Purity, #phi" );
  //purmg->GetXaxis()->SetTitle("MCscale");
  //purmg->GetYaxis()->SetTitle("Purity, #phi");
  purmg->GetYaxis()->SetTitleOffset(1.2);
  purleg->Draw();
  purvnoise->Update();
  purvnoise->SaveAs( "/home/mthiesse/PurityAnalysis/Scripts/png/purvnoise.png" );
  delete purvnoise;

  ///////////////////

  TCanvas * effvnoise = new TCanvas( "effvnoise", "effvnoise", 2000, 1600 );
  effvnoise->SetMargin(0.12,0.02,0.11,0.02); //l,r,b,t
  TMultiGraph * effmg = new TMultiGraph();
  TLegend * effleg = new TLegend( 0.55, 0.15, 0.9, 0.5 );
  color = 53;
  for ( auto igr : effgraphs )
    {
      TGraphAsymmErrors * effgraph = igr.second->CreateGraph();
      //if ( effgraph->GetN()==0 ) continue;
      effgraph->SetLineColor( color);
      effgraph->SetLineWidth(4);
      effgraph->SetMarkerStyle( 21 );
      effgraph->SetMarkerSize( 2 );
      effgraph->SetMarkerColor( color );
      effmg->Add( effgraph, "LPE" );
      effleg->AddEntry( effgraph, TString::Format( "%i#mus Lifetime", igr.first ), "lep" );
      color+=10;
    }
  effmg->Draw( "A" );
  effmg->SetTitle( ";MCscale;Efficiency, #epsilon_{S}" );
  effvnoise->Update();
  effmg->GetYaxis()->SetTitleOffset(1.2);
  effleg->Draw();
  effvnoise->Update();
  effvnoise->SaveAs( "/home/mthiesse/PurityAnalysis/Scripts/png/effvnoise.png" );
  delete effvnoise;

  ////////////////////////////

  TCanvas * canvpurdriftlife = new TCanvas("canvpurdriftlife","",2000,1600);
  canvpurdriftlife->SetMargin(0.12,0.02,0.11,0.02); //l,r,b,t
  TMultiGraph * purdriftlifemg = new TMultiGraph();
  TLegend * purdriftlifeleg = new TLegend(0.12,0.12,0.47,0.47);
  color = 53;
  for (auto igr : purdriftlife)
    {
      TGraphAsymmErrors * pdl = igr.second->CreateGraph();
      //if (igr.second->GetN()==0) continue;
      pdl->SetLineColor(color);
      pdl->SetLineWidth(4);
      pdl->SetMarkerStyle(21);
      pdl->SetMarkerSize(2);
      pdl->SetMarkerColor(color);
      purdriftlifemg->Add(pdl,"LPE");
      purdriftlifeleg->AddEntry(pdl,TString::Format("%i#mus Lifetime",igr.first),"lep");
      color +=10;
    }
  purdriftlifemg->Draw("A");
  purdriftlifemg->SetTitle(";Drift Distance [cm];Purity, #phi");
  canvpurdriftlife->Update();
  purdriftlifemg->GetYaxis()->SetTitleOffset(1.2);
  purdriftlifeleg->Draw();
  canvpurdriftlife->Update();
  canvpurdriftlife->SaveAs("/home/mthiesse/PurityAnalysis/Scripts/png/canvpurdriftlife.png");
  delete canvpurdriftlife;

  //////////////////////////////////////////////////////////////////

  TCanvas * canveffdriftlife = new TCanvas("canveffdriftlife","",2000,1600);
  canveffdriftlife->SetMargin(0.12,0.02,0.11,0.02); //l,r,b,t
  TMultiGraph * effdriftlifemg = new TMultiGraph();
  TLegend * effdriftlifeleg = new TLegend(0.12,0.12,0.47,0.47);
  color = 53;
  for (auto igr : effdriftlife)
    {
      TGraphAsymmErrors * edl = igr.second->CreateGraph();
      //if (edl->GetN()==0) continue;
      edl->SetLineColor(color);
      edl->SetLineWidth(4);
      edl->SetMarkerStyle(21);
      edl->SetMarkerSize(2);
      edl->SetMarkerColor(color);
      effdriftlifemg->Add(edl,"LPE");
      effdriftlifeleg->AddEntry(edl,TString::Format("%i#mus Lifetime",igr.first),"lep");
      color +=10;
    }
  effdriftlifemg->Draw("A");
  effdriftlifemg->SetTitle(";Drift Distance [cm];Efficiency, #epsilon_{S}");
  canveffdriftlife->Update();
  effdriftlifemg->GetYaxis()->SetTitleOffset(1.2);
  effdriftlifeleg->Draw();
  canveffdriftlife->Update();
  canveffdriftlife->SaveAs("/home/mthiesse/PurityAnalysis/Scripts/png/canveffdriftlife.png");
  delete canveffdriftlife;

  //////////////////////////////////////////////////////////////////

  //std::map<Int_t, TCanvas*> purcanvmap;

  UInt_t numcols = interesting_mcscale.size();
  Int_t palette[numcols];
  Double_t Red[] =    {0.831, 0.0,   0.0,   0.133, 1.0};
  Double_t Green[] =  {0.349, 0.267, 1.0,   1.0,   0.933};
  Double_t Blue[]   = {0.329, 1.0,   0.800, 0.0,   0.0};
  Double_t Length[] = {0., .25, .50, .75, 1.0};
  Int_t FI = TColor::CreateGradientColorTable( 5, Length, Red, Green, Blue, numcols );
  for ( UInt_t i=0; i<numcols; i++ ) palette[i] = FI+i;
/*
   for ( auto igrp : purhist )
    {
      TString tit = TString::Format( "purhist_%ius", igrp.first );
      purcanvmap[igrp.first] = new TCanvas( tit.Data(), tit.Data(), 2000, 1600 );
      purcanvmap[igrp.first]->cd();
      TLegend * purleg = new TLegend( 0.1, 0.5, 0.35, 0.9 );
      purleg->SetHeader( TString::Format( "eLifetime=%i#mus", igrp.first ) );
      TIter next( purleg->GetListOfPrimitives() );
      TLegendEntry *first = ( TLegendEntry* )next();
      first->SetTextColor( 2 );
      //first->SetTextSize( 0.025 );
      Int_t col = 0;
      Int_t ngraphs = 0;
      TMultiGraph * mg = new TMultiGraph();
      for ( auto igr : igrp.second )
        {
          TGraph * gr = new TGraph(igr.second);
          if ( gr->GetN()==0 ) continue;
          gr->SetLineColor( palette[col] );
          gr->SetLineWidth(4);
          mg->Add( gr, "L" );
          purleg->AddEntry( gr, TString::Format( "MCScale=%.1f", mcscale[igr.first] ), "l" );
   ++ngraphs;
   ++col;
        }
      if ( ngraphs!=0 )
        {
          mg->Draw( "A" );
          mg->SetTitle( ";Purity, #phi;" );
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
 */
///////////////////////////////
/*
   std::map<Int_t, TCanvas*> effcanvmap;

   for ( auto igrp : effhist )
    {
      TString tit = TString::Format( "effhist_%ius", igrp.first );
      effcanvmap[igrp.first] = new TCanvas( tit.Data(), tit.Data(), 2000, 1600 );
      effcanvmap[igrp.first]->cd();
      TLegend * effleg = new TLegend( 0.4, 0.58, 0.65, 0.98 );
      effleg->SetHeader( TString::Format( "eLifetime=%i#mus", igrp.first ) );
      TIter next( effleg->GetListOfPrimitives() );
      TLegendEntry *first = ( TLegendEntry* )next();
      first->SetTextColor( 2 );
      //first->SetTextSize( 0.025 );
      Int_t col = 0;
      Int_t ngraphs = 0;
      TMultiGraph * mg = new TMultiGraph();
      for ( auto igr : igrp.second )
        {
          TGraph * gr = new TGraph(igr.second);
          if ( gr->GetN()==0 ) continue;
          gr->SetLineColor( palette[col] );
          gr->SetLineWidth(4);
          mg->Add( gr, "L" );
          effleg->AddEntry( gr, TString::Format( "MCScale=%.1f", mcscale[igr.first] ), "l" );
   ++ngraphs;
   ++col;
        }
      if ( ngraphs!=0 )
        {
          mg->Draw( "A" );
          mg->SetTitle( ";Efficiency, #epsilon_{S};" );
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
 */
///////////////////////////
/*
   std::map<Int_t, TCanvas*> purdriftcanvmap;

   for ( auto igrp : purdrift )
    {
      TString tit = TString::Format( "purvdrift_%ius", igrp.first );
      purdriftcanvmap[igrp.first] = new TCanvas( tit.Data(), tit.Data(), 2200, 1600 );
      purdriftcanvmap[igrp.first]->SetRightMargin(0.15);
      purdriftcanvmap[igrp.first]->SetLeftMargin(0.09);
      purdriftcanvmap[igrp.first]->cd();
      TLegend * purleg = new TLegend( 0.78, 0.5, 0.999, 0.98 );
      purleg->SetHeader( TString::Format( "eLifetime=%i#mus", igrp.first ) );
      TIter next( purleg->GetListOfPrimitives() );
      TLegendEntry *first = ( TLegendEntry* )next();
      first->SetTextColor( 2 );
      //first->SetTextSize( 0.025 );
      Int_t col = 0;
      Int_t ngraphs = 0;
      TMultiGraph * mg = new TMultiGraph();
      for ( auto igr : igrp.second )
        {
          TGraphAsymmErrors * pd = igr.second->CreateGraph();
          //if ( pd->GetN()==0 ) continue;
          pd->SetLineColor( palette[col] );
          pd->SetLineWidth(4);
          if ( fabs( mcscale[igr.first]-1.0 )<0.01 ) { pd->SetMarkerStyle( 29 ); pd->SetMarkerSize( 4 ); }
          else { pd->SetMarkerStyle( 21 ); pd->SetMarkerSize( 2 ); }
          pd->SetMarkerColor( palette[col] );
          mg->Add( pd, "LP" );
          purleg->AddEntry( pd, TString::Format( "MCScale=%.1f", mcscale[igr.first] ), "lep" );
   ++ngraphs;
   ++col;
        }
      if ( ngraphs!=0 )
        {
          mg->Draw( "A" );
          mg->SetTitle( ";Mean Drift Distance [cm];Purity, #phi" );
          purleg->Draw();
          purdriftcanvmap[igrp.first]->Update();
          purdriftcanvmap[igrp.first]->SaveAs( TString::Format( "/home/mthiesse/PurityAnalysis/Scripts/png/%s.png", tit.Data() ) );
          delete purdriftcanvmap[igrp.first];
        }
      else
        {
          delete purdriftcanvmap[igrp.first];
        }
    }
 */
//////////////////////////
/*
   std::map<Int_t, TCanvas*> effdriftcanvmap;

   for ( auto igrp : effdrift )
    {
      TString tit = TString::Format( "effvdrift_%ius", igrp.first );
      effdriftcanvmap[igrp.first] = new TCanvas( tit.Data(), tit.Data(), 2200, 1600 );
      effdriftcanvmap[igrp.first]->SetRightMargin(0.15);
      effdriftcanvmap[igrp.first]->SetLeftMargin(0.09);
      effdriftcanvmap[igrp.first]->cd();
      TLegend * effleg = new TLegend( 0.78, 0.5, 0.999, 0.98 );
      effleg->SetHeader( TString::Format( "eLifetime=%i#mus", igrp.first ) );
      TIter next( effleg->GetListOfPrimitives() );
      TLegendEntry *first = ( TLegendEntry* )next();
      first->SetTextColor( 2 );
      //first->SetTextSize( 0.025 );
      Int_t col = 0;
      Int_t ngraphs = 0;
      TMultiGraph * mg = new TMultiGraph();
      for ( auto igr : igrp.second )
        {
          TGraphAsymmErrors * ed = igr.second->CreateGraph();
          //if ( ed->GetN()==0 ) continue;
          ed->SetLineColor( palette[col] );
          ed->SetLineWidth(4);
          if ( fabs( mcscale[igr.first]-1.0 )<0.01 ) { ed->SetMarkerStyle( 29 ); ed->SetMarkerSize( 4 ); }
          else { ed->SetMarkerStyle( 21 ); ed->SetMarkerSize( 2 ); }
          ed->SetMarkerColor( palette[col] );
          mg->Add( ed, "LP" );
          effleg->AddEntry( ed, TString::Format( "MCScale=%.1f", mcscale[igr.first] ), "lep" );
   ++ngraphs;
   ++col;
        }
      if ( ngraphs!=0 )
        {
          mg->Draw( "A" );
          mg->SetTitle( ";Mean Drift Distance [cm];Efficiency, #epsilon_{S}" );
          effleg->Draw();
          effdriftcanvmap[igrp.first]->Update();
          effdriftcanvmap[igrp.first]->SaveAs( TString::Format( "/home/mthiesse/PurityAnalysis/Scripts/png/%s.png", tit.Data() ) );
          delete effdriftcanvmap[igrp.first];
        }
      else
        {
          delete effdriftcanvmap[igrp.first];
        }
    }
 */
///////////////////////////////////

  std::map<Int_t, std::vector<Double_t> > binedge;
  std::map<Int_t, std::vector<Double_t> > bininvcontent;

  TCanvas * resolutioncanv = new TCanvas("resolution","resolution",2000,1600);
  resolutioncanv->SetBottomMargin(0.11);
  resolutioncanv->SetMargin(0.11,0.1,0.11,0.02); //l,r,b,t
  TH1F * empty = resolutioncanv->DrawFrame(-100,0.,8100,150);
  empty->SetXTitle("Q_{MC} [ADC]");
  empty->SetYTitle("Charge Resolution, #sigma_{q} [%]");
  TLegend * resolutionleg = new TLegend( 0.58, 0.569, 0.9, 0.9 );
  TMultiGraph * resolutionmg = new TMultiGraph();
  Int_t col = 0;
  Int_t ngraphs = 0;
  Int_t nbinsless20entries = 0;
  std::vector<UInt_t> iscaleless20;
  std::vector<Int_t> ibinless20;
  for (auto ig : resolution)
    {
      TGraphErrors * gr = new TGraphErrors();
      //std::cout << "mcscale=" << ig.first << "   [" << std::endl;
      int pt = 0;
      for (auto igr : ig.second)
        {
          Double_t mcq = (*(residual.begin())).second->GetBinCenter(igr.first);
          Double_t histmax = igr.second->GetBinContent(igr.second->GetMaximumBin());
          Double_t binlow = igr.second->GetBinCenter(igr.second->FindFirstBinAbove(0.5*histmax));
          Double_t binhigh = igr.second->GetBinCenter(igr.second->FindLastBinAbove(0.5*histmax));
          TF1 * gaus = new TF1("gaus","gaus",binlow,binhigh);
          igr.second->Fit("gaus","RQ0");
          Double_t res = gaus->GetParameter(2);
          Double_t reserr = gaus->GetParError(2);
          if (igr.second->GetEntries()<20) {
              ++nbinsless20entries;
              iscaleless20.push_back(ig.first);
              ibinless20.push_back(igr.first);
            }
          //if (res > 0 && igr.second->GetEntries()>1)
          if (binhigh>binlow && igr.second->GetEntries()>1)
            {
              binedge[ig.first].push_back(igr.second->GetBinLowEdge(pt));
              bininvcontent[ig.first].push_back(1/res);
              gr->SetPoint(pt,mcq,res*100);
              gr->SetPointError(pt,0,reserr*100);
              //gr->SetPointError(pt,0,2*pow((binhigh-binlow)*100,4)/(igr.second->GetEntries()-1));
              ++pt;
            }
          //std::cout << "   mcq=" << mcq << "  N=" << igr.second->GetEntries() << "  FWHM=" << binhigh-binlow << "  gauswidth=" << res << " +/- " << reserr << std::endl;
        }
      //std::cout << "                                   ]" << std::endl;
      gr->SetLineColor(palette[col]);
      gr->SetMarkerStyle(20);
      gr->SetLineWidth(4);
      gr->SetMarkerSize(2);
      gr->SetMarkerColor(palette[col]);
      resolutionmg->Add(gr,"EP");
      resolutionleg->AddEntry(gr,TString::Format("MCScale=%.1f", mcscale[ig.first]),"ep");
      ++ngraphs;
      ++col;
    }
  if (ngraphs != 0)
    {
      resolutionmg->Draw("ep");
      resolutionleg->Draw();
      resolutioncanv->Update();
      resolutioncanv->SaveAs( "/home/mthiesse/PurityAnalysis/Scripts/png/resolution.png");
    }
  //delete resolutioncanv;

  ////////////////////////////////////////
  // Do some magic to get exponentially increasing bin widths (at least according to the residual plot...)
  for (UInt_t iscale = 0; iscale < interesting_mcscale.size(); ++iscale)
    {
      //std::vector<Double_t> binc = bininvcontent[interesting_mcscale[iscale]];
      //Int_t nbins = binc.size();
      Int_t nbins = 200;
      //Double_t sum = 0;
      //for (int i = 0; i < nbins; ++i)
      //  {
      //    sum += binc[i];
      //std::cout << "binc[" << i << "]=" << binc[i] << std::endl;
      //  }
      //std::cout << "   sum = " << sum << std::endl;
      std::vector<Double_t> xbins(nbins);
      Double_t lowedge = 0;
      Double_t hiedge = 8000;
      Double_t thislow = lowedge;
      Double_t dx = (hiedge-lowedge)/nbins;
      for (int i = 0; i < nbins; ++i)
        {
          //std::cout << "binc[" << i << "]=" << binc[i];
          //binc[i] = binc[i]*8000/sum;
          //std::cout << "     NEW binc[" << i << "]=" << binc[i];
          //std::cout << "     BinLowEdge = " << lowedge << std::endl;
          xbins[i] = thislow;
          //lowedge += binc[i];
          thislow += dx; //200 = number of bins
        }
      //effcharge[interesting_mcscale[iscale]] = new TProfile(TString::Format("effcharge_%.1f",mcscale[interesting_mcscale[iscale]]),";Q_{MC} [ADC];Efficiency, #epsilon_{S}",xbins.size()-1,xbins.data(),0,8000,"e");
      //effcharge_stdev[interesting_mcscale[iscale]] = new TProfile(TString::Format("effcharge_stdev_%.1f",mcscale[interesting_mcscale[iscale]]),";Q_{MC} [ADC];Efficiency, #epsilon_{S}",xbins.size()-1,xbins.data(),0,8000,"s");
      effcharge_correct[interesting_mcscale[iscale]] = new TEfficiency(TString::Format("effcharge_correct_%.1f",mcscale[interesting_mcscale[iscale]]),";Q_{MC} [ADC]; Efficiency, #epsilon_{S}",xbins.size()-1,xbins.data());
      //effcharge_short[interesting_mcscale[iscale]] = new TEfficiency(TString::Format("effcharge_short_%.1f",mcscale[interesting_mcscale[iscale]]),";Q_{MC} [ADC]; Efficiency, #epsilon_{S}",xbins.size()-1,xbins.data());
      //effcharge_long[interesting_mcscale[iscale]] = new TEfficiency(TString::Format("effcharge_long_%.1f",mcscale[interesting_mcscale[iscale]]),";Q_{MC} [ADC]; Efficiency, #epsilon_{S}",xbins.size()-1,xbins.data());
      numeffcharge[interesting_mcscale[iscale]] = new TH1D(TString::Format("numeffcharge_%.1f",mcscale[interesting_mcscale[iscale]]),";Q_{MC} [ADC]; # Entries",600,lowedge,hiedge);

      effcharge_correct[interesting_mcscale[iscale]]->SetStatisticOption(TEfficiency::kBUniform);
      //effcharge_short[interesting_mcscale[iscale]]->SetStatisticOption(TEfficiency::kBUniform);
      //effcharge_long[interesting_mcscale[iscale]]->SetStatisticOption(TEfficiency::kBUniform);
      effcharge_correct[interesting_mcscale[iscale]]->SetConfidenceLevel(0.68);
      //effcharge_short[interesting_mcscale[iscale]]->SetConfidenceLevel(0.951);
      //effcharge_long[interesting_mcscale[iscale]]->SetConfidenceLevel(0.951);

      for (size_t i = 0; i < effchargedata[interesting_mcscale[iscale]].size(); ++i)
        {
          effchg e = effchargedata[interesting_mcscale[iscale]][i];
          Bool_t foundboth = e.fb;
          if (foundboth && e.recoq<0) foundboth = false;
          //if (e.recoq > 1e5 || e.recoq < -1e3) continue;
          //effcharge[interesting_mcscale[iscale]]->Fill(e.mcq,e.fb);
          //effcharge_stdev[interesting_mcscale[iscale]]->Fill(e.mcq,e.fb);
          effcharge_correct[interesting_mcscale[iscale]]->Fill(e.fb,e.mcq);
          //if (e.mct < 500) effcharge_short[interesting_mcscale[iscale]]->Fill(e.fb,e.mcq);
          //if (e.mct > 3500) effcharge_long[interesting_mcscale[iscale]]->Fill(e.fb,e.mcq);
          if (e.mcq > 200) numeffcharge[interesting_mcscale[iscale]]->Fill(e.mcq);
        }
      //effcharge[interesting_mcscale[iscale]]->Rebin();
      //effcharge_stdev[interesting_mcscale[iscale]]->Rebin();

      numeffcharge[interesting_mcscale[iscale]]->Scale(1/numeffcharge[interesting_mcscale[iscale]]->GetBinContent(numeffcharge[interesting_mcscale[iscale]]->GetMaximumBin()));

    }

  ////////////////////////////////////////

  TCanvas * effchargecanv = new TCanvas("effcharge","effcharge",2000,1600);
  effchargecanv->SetBottomMargin(0.11);
  effchargecanv->SetMargin(0.11,0.1,0.11,0.02); //l,r,b,t
  TLegend * effchargeleg = new TLegend( 0.508, 0.2, 0.89, 0.506 );
  col = 0;
  ngraphs = 0;
  TMultiGraph * effchargemg = new TMultiGraph();
  for ( auto igr : numeffcharge )
    {
      TGraph * grnum = new TGraph(igr.second);
      grnum->SetFillColor(kGray);
      effchargemg->Add(grnum,"B");
      effchargeleg->AddEntry(grnum,"Histogram of Q_{MC}","F");
      break;
    }
  for ( auto igr : effcharge_correct )
    {
      TGraphAsymmErrors * grerr = igr.second->CreateGraph();
      grerr->SetLineColor( palette[col] );
      grerr->SetLineWidth( 4 );
      grerr->SetMarkerStyle(20);
      grerr->SetMarkerColor(palette[col]);
      grerr->SetMarkerSize(1.5);
      effchargemg->Add(grerr,"EP");
      effchargeleg->AddEntry( grerr, TString::Format( "MCScale=%.1f", mcscale[igr.first] ), "e3p" );
      ++ngraphs;
      ++col;
    }
  if ( ngraphs!=0 )
    {
      effchargecanv->cd();
      effchargemg->Draw("A");
      effchargemg->SetTitle( ";Q_{MC} [ADC];Efficiency, #epsilon_{S}" );
      effchargeleg->Draw();
      effchargecanv->Update();
      effchargecanv->SaveAs( "/home/mthiesse/PurityAnalysis/Scripts/png/effcharge.png");
    }
  //delete effchargecanv;


////////////////////////////////////////////////

  TCanvas * residualcanv = new TCanvas("residual","residual",2000,1600);
  residualcanv->SetBottomMargin(0.11);
  TLegend * residualleg = new TLegend(0.6,0.1,0.9,0.35);
  col = 0;
  ngraphs = 0;
  TMultiGraph * residualmg = new TMultiGraph();
  for ( auto igr : residual )
    {
      TGraphErrors * gr = new TGraphErrors(igr.second);
      gr->SetLineColor( palette[col] );
      gr->SetLineWidth(4);
      gr->SetMarkerStyle(20);
      gr->SetMarkerColor(palette[col]);
      gr->SetMarkerSize(1);
      residualmg->Add(gr,"EP");
      residualleg->AddEntry( gr, TString::Format( "MCScale=%.1f", mcscale[igr.first] ), "ep" );
      ++ngraphs;
      ++col;
    }
  if ( ngraphs!=0 )
    {
      residualmg->Draw( "A" );
      residualmg->SetTitle( ";Q_{MC} [ADC];Charge Residual, R_{q}" );
      residualleg->Draw();
      residualcanv->Update();
      residualcanv->SaveAs( "/home/mthiesse/PurityAnalysis/Scripts/png/residual.png");
    }
  //delete residualcanv;

//////////////////////////////////

  gStyle->SetOptStat(0);
  TCanvas * migrationcanv = new TCanvas("migrationcanv","migrationcanv",2000,1600);
  mm->SetLineColor(9);
  mm->SetLineWidth(4);
  mm->SetMarkerStyle(20);
  mm->SetMarkerColor(9);
  mm->SetMarkerSize(1);
  mm->Draw("ep");
  mm_stdev->SetFillColor(9);
  mm_stdev->SetFillColorAlpha(9,0.3);
  mm_stdev->Draw("E3 same");
  migrationcanv->Update();
  migrationcanv->SaveAs("/home/mthiesse/PurityAnalysis/Scripts/png/migrationmatrix.png");
  delete migrationcanv;

  TGraph * mm_graph = new TGraph(mm);
  std::cout << "Migration Matrix: correlation = " << mm_graph->GetCorrelationFactor() << std::endl;

  ///////////////////////////////

  TCanvas * reshistcanv = new TCanvas("reshist","reshist",2000,1600);
  reshistcanv->cd();
  TLegend * reshistleg = new TLegend( 0.65, 0.65, 0.9, 0.9 );
  col = 0;
  ngraphs = 0;
  TMultiGraph * reshistmg = new TMultiGraph();
  for ( auto igr : reshist )
    {
      TGraphErrors * gr = new TGraphErrors(igr.second);
      if ( gr->GetN()==0 ) continue;
      gr->SetLineColor( palette[col] );
      gr->SetLineWidth(4);
      reshistmg->Add( gr, "EP" );
      reshistleg->AddEntry( gr, TString::Format( "MCScale=%.1f", mcscale[igr.first] ), "ep" );
      ++ngraphs;
      ++col;
    }
  if ( ngraphs!=0 )
    {
      reshistmg->Draw( "A" );
      reshistmg->SetTitle( ";Charge residual, R_{q};");
      reshistleg->Draw();
      reshistcanv->Update();
      reshistcanv->SaveAs( "/home/mthiesse/PurityAnalysis/Scripts/png/reshist.png");
    }
  //delete reshistcanv;


  std::cout << "nbinsless20entries=" << nbinsless20entries << std::endl;
  for (int i = 0; i < nbinsless20entries; ++i)
    {
      std::cout << "iscale=" << iscaleless20[i] << "  ibin=" << ibinless20[i] << std::endl;
    }

}
