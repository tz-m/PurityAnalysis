#ifndef PURITYANALYSIS_H
#define PURITYANALYSIS_H

#include "Analysis.h"

#include "THStack.h"
#include "TFitResult.h"
#include "TSystem.h"
#include "TImage.h"
#include "TPaveLabel.h"
#include "TVirtualFitter.h"
#include "TColor.h"
#include "TFrame.h"
#include "TProfile.h"
#include "TUnfoldDensity.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "RooGlobalFunc.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooArgSet.h"
#include "RooAddition.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooMinuit.h"
#include "RooAbsPdf.h"
#include "RooPolynomial.h"
#include "RooProduct.h"
#include "Math/VavilovAccurate.h"
#include "RooFunctorBinding.h"
#include "Math/Functor.h"

struct ChargeHistogram
{
  TH1F * hist;
  TH1F * histuw;
  Float_t binlow;
  Float_t bincenter;
  Float_t binhigh;
  Float_t binwidth;
};

struct FitResult
{
  Int_t binnum;
  Float_t landmpv;
  Float_t landwidth;
  Float_t gaussmean;
  Float_t gausswidth;
  Float_t gaussmean_l;
  Float_t gaussmean_v;
  Float_t gausswidth_l;
  Float_t gausswidth_v;
  Float_t vavkappa;
  Float_t vavbeta2;
  Float_t vavmean;
  Float_t vavwidth;
  Float_t vavamp;
  Float_t landmpverr;
  Float_t landwidtherr;
  Float_t gaussmeanerr;
  Float_t gausswidtherr;
  Float_t gaussmean_lerr;
  Float_t gausswidth_lerr;
  Float_t gaussmean_verr;
  Float_t gausswidth_verr;
  Float_t vavkappaerr;
  Float_t vavbeta2err;
  Float_t vavmeanerr;
  Float_t vavwidtherr;
  Float_t vavamperr;
  Float_t bincenter;
  Float_t binwidtherr;
  Float_t lxg_chi2ndf;
  Float_t vxg_chi2ndf;
  Float_t pg_chi2ndf;
  Float_t lxglow;
  Float_t lxghigh;
  Bool_t fitsuccess;
};


struct Vavilov_Func {
  Vavilov_Func() {
  }

  ROOT::Math::VavilovAccurate pdf;

  double operator() ( const double *x ) {
    return x[5]*(  pdf.Pdf(  ( x[0]-x[3] )/x[4], x[1], x[2] ) );
  }
};

struct HitSave {
  Float_t hitt;
  Float_t dqdx;
  Float_t weight;
  Float_t unscaleddqdx;
};


class PurityAnalysis : public Analysis {
public:
  PurityAnalysis();
  void run();

private:

  int getBinNumber( Float_t t, const std::vector<ChargeHistogram> & tbh );

};

PurityAnalysis::PurityAnalysis()
{
}

int PurityAnalysis::getBinNumber( Float_t t, const std::vector<ChargeHistogram> & tbh )
{
  for ( UInt_t i_hist = 0; i_hist < tbh.size(); i_hist++ )
    {
      if ( t >= tbh[i_hist].binlow && t < tbh[i_hist].binhigh )
        {
          return ( int )i_hist;
        }
    }
  return -1;
}

void PurityAnalysis::run()
{
  bool DoPull = cuts.GetAnalysisCutBool( "PurityAnalysis", "DoPull" );
  bool DoResidual = cuts.GetAnalysisCutBool( "PurityAnalysis", "DoResidual" );
  bool DoLxG = cuts.GetAnalysisCutBool( "PurityAnalysis", "DoLxG" );
  bool DoVxG = cuts.GetAnalysisCutBool( "PurityAnalysis", "DoVxG" );
  bool DoPG = cuts.GetAnalysisCutBool( "PurityAnalysis", "DoPG" );
  float Nbinsfloat = cuts.GetAnalysisCutFloat( "PurityAnalysis", "Nbins" );
  std::string fileoutname = cuts.GetAnalysisCutString( "PurityAnalysis", "FileOut");

  UInt_t Nbins = ( UInt_t )Nbinsfloat;
  Float_t drifttimemax = 2012;
  Float_t ADCcutofflow, ADCcutoffhigh;
  ADCcutofflow = cuts.GetAnalysisCutFloat( "PurityAnalysis", "ADCCutoffLow" );
  ADCcutoffhigh = cuts.GetAnalysisCutFloat( "PurityAnalysis", "ADCCutoffHigh" );

  TFile * scalefile = TFile::Open("scaleQ.root","READ");
  TProfile * scaleQ = (TProfile*)scalefile->Get("scaleQ");
  scaleQ->SetDirectory(0);
  scalefile->Close();

  TFile * comparefile = TFile::Open("compareQ.root","READ");
  TProfile * compareQ = (TProfile*)comparefile->Get("compareQ");
  compareQ->SetDirectory(0);
  comparefile->Close();

  TFile * efffile = TFile::Open("effQ.root","READ");
  TProfile * effQ = (TProfile*)efffile->Get("effQ");
  effQ->SetDirectory(0);
  efffile->Close();

  Double_t mcscale = cuts.GetAnalysisCutFloat("PurityAnalysis","mcscale");

  TFile * fileout = TFile::Open( fileoutname.c_str(), "RECREATE" );
  Int_t runnum;
  Int_t event;
  Float_t hitt;
  Int_t binnum;
  Bool_t foundrealhit;
  Bool_t assumedhit;
  Float_t dqdx;
  Int_t channel;
  Int_t tpc;
  TTree * hitstree = new TTree( "hits", "Hit Information" );
  hitstree->Branch( "runnum", &runnum, "runnum/I" );
  hitstree->Branch( "event", &event, "event/I" );
  hitstree->Branch( "hitt", &hitt, "hitt/F" );
  hitstree->Branch( "binnum", &binnum, "binnum/I" );
  hitstree->Branch( "foundrealhit", &foundrealhit, "foundrealhit/O" );
  hitstree->Branch( "assumedhit", &assumedhit, "assumedhit/O" );
  hitstree->Branch( "dqdx", &dqdx, "dqdx/F" );
  hitstree->Branch( "channel", &channel, "channel/I" );
  hitstree->Branch( "tpc", &tpc, "tpc/I" );

  TH1F * found = new TH1F( "found", "Hits Found", 200, ADCcutofflow, ADCcutoffhigh );
  TH1F * assumed = new TH1F( "assumed", "Assumed Hits", 200, ADCcutofflow, ADCcutoffhigh );

  RooRealVar charge( "charge", "dQ/dx ( ADC/cm )", ADCcutofflow, ADCcutoffhigh );
  RooRealVar tdrift( "tdrift", "tdrift", 0, drifttimemax );
  RooRealVar weight( "weight", "weight", -10, 10 );
  RooDataSet * chgdata = new RooDataSet( "chgdata", "chgdata", RooArgSet( tdrift, charge ), RooFit::StoreError( RooArgSet( charge ) ) /* , RooFit::WeightVar( "weight" )*/ );
  chgdata->Print();

  std::vector<ChargeHistogram> timeBinHist( Nbins );
  for ( UInt_t i = 0; i < Nbins; i++ )
    {
      timeBinHist[i].binlow = i*( drifttimemax/Nbins );
      timeBinHist[i].binhigh = ( i+1 )*( drifttimemax/Nbins );
      timeBinHist[i].bincenter = ( timeBinHist[i].binlow + timeBinHist[i].binhigh )/2.0;
      timeBinHist[i].binwidth = timeBinHist[i].binhigh - timeBinHist[i].binlow;
      timeBinHist[i].hist = new TH1F( TString::Format( "h%i", i ), TString::Format( "RobustHitFinder, Bin %i, %.2f <= t < %.2f ( us )", i, timeBinHist[i].binlow, timeBinHist[i].binhigh ), 200, ADCcutofflow, ADCcutoffhigh );
      timeBinHist[i].hist->Sumw2();
      timeBinHist[i].histuw = new TH1F( TString::Format( "huw%i", i ), TString::Format( "RobustHitFinder, Bin %i, %.2f <= t < %.2f ( us )", i, timeBinHist[i].binlow, timeBinHist[i].binhigh ), 200, ADCcutofflow, ADCcutoffhigh );
      timeBinHist[i].histuw->Sumw2();
    }

  Float_t secondbinstart = timeBinHist[1].binlow;
  Float_t penultimatebinend = timeBinHist[Nbins-2].binhigh;

  std::map<Int_t, std::vector<HitSave> > hitmap;
  std::map<UInt_t, std::pair<Int_t, Int_t> > fabinmap;
  for ( UInt_t i = 0; i < Nbins; i++ )
    {
      fabinmap[i] = std::make_pair( 0, 0 );
    }

  for ( auto const & hititr : *( file->GetHitMap() ) )
    {
      const types::HitInfo * hit = &( hititr.second );
      if ( hit->run >= 15483 && hit->run <= 15502 ) continue;
      if ( hit->run == 15593 ) continue;
      if ( hit->run >= 15634 && hit->run <= 15643 ) continue;
      if ( hit->run >= 15664 && hit->run <= 15674 ) continue;
      if ( hit->run == 15823 || hit->run == 15914 || hit->run == 15940 || hit->run == 15980 || hit->run == 16025 || hit->run == 16083 || hit->run == 16147 || hit->run == 16586 ) continue;
      if ( hit->run >= 16593 && hit->run <= 16596 ) continue;
      if ( hit->channel == 566 || hit->channel == 885 || hit->channel == 1547 ) continue;
      if ( !( hit->fitsuccess ) ) continue;
      //if ( hit->tpc != 1 ) continue;//if you uncomment this, then must re-consider the fitndf cut below
      if ( hit->width > 400 ) continue;
      if ( hit->fitchi2 < -9999 || hit->fitchi2 > 10 ) continue;
      if ( hit->integral < -500e4 ) continue;
      if ( hit->rms > 40 ) continue;
      if ( hit->rms < 10 ) continue;
      if ( fabs( hit->baseline ) > 20 ) continue;
      if ( hit->fitconstant<-49 ) continue;
      if ( hit->fitconstant>249 ) continue;
      if ( hit->fitndf<50 ) continue;  // basically cut on 50 hits per track
      if ( fabs( hit->fitquadratic )>0.000199 ) continue;  // skip if fit parameter is maximal
      if ( cuts.ChannelPass( hit ) &&
           cuts.CounterPass( hit ) )
        {
          runnum = hit->run;
          event = hit->event;
          hitt = hit->hitt;
          binnum = getBinNumber( hitt, timeBinHist );
          foundrealhit = !hit->assumedhit && hit->fitrealhit;
          assumedhit = hit->assumedhit;
          Double_t hitintegral = hit->integral;
          dqdx =  hitintegral / hit->segmentlength;
          channel = hit->channel;
          tpc = hit->tpc;
          hitstree->Fill();

          if ( foundrealhit )
            {
              found->Fill( dqdx );
              fabinmap[binnum].first = 1 + fabinmap[binnum].first;
            }
          if ( assumedhit )
            {
              assumed->Fill( dqdx );
              fabinmap[binnum].second = 1 + fabinmap[binnum].second;
            }
          if ( hit->assumedhit || hit->fitrealhit )
            {
              HitSave hs;
              hs.hitt = hitt;
              hs.unscaleddqdx = dqdx;
              hs.dqdx = dqdx;
              hs.weight = 1; // / (compareQ->GetBinError(compareQ->FindBin(hitintegral)) * effQ->GetBinContent(effQ->FindBin(hitintegral)));

              // / effQ->GetBinContent(effQ->FindBin(hit->integral / scaleQ->GetBinContent(scaleQ->FindBin(hit->integral)) )); //compareQ->GetBinError(compareQ->FindBin(hit->integral / scaleQ->GetBinContent(scaleQ->FindBin(hit->integral)))); // / (compareQ->GetBinError(compareQ->FindBin(hit->integral /* / scaleQ->GetBinContent(scaleQ->FindBin(hit->integral)) */ )) * effQ->GetBinContent(effQ->FindBin(hit->integral /* / scaleQ->GetBinContent(scaleQ->FindBin(hit->integral))*/ )) );
              hitmap[channel].push_back( hs );

              //std::cout << "HitT: " << hitt << "  dQdx: " << dqdx << "  Weight: " << hs.weight << std::endl;
            }
          //if ( hit->fitrealhit ) hitmap[channel].push_back( std::make_pair( hitt, dqdx ) );
        }
    }
  hitstree->Write();

  TCanvas * canvFAratio = new TCanvas( "canvFAratio", "FAratio vs drift", 2000, 1600 );
  TPad *pad = new TPad( "pad", "", 0, 0, 1, 1 );
  pad->SetFillColor( 0 );
  pad->Draw();
  pad->cd();
  TH1F *hr = pad->DrawFrame( 0, 0, drifttimemax, 1 );
  hr->SetXTitle( "Drift Time ( us )" );
  hr->SetYTitle( "#Assumed / ( #Assumed + #Found )" );
  pad->GetFrame()->SetFillColor( 0 );
  pad->GetFrame()->SetBorderSize( 12 );
  TGraph * FAratio = new TGraph();
  TGraph * sumFA = new TGraph();
  Double_t maxsum = 0, minsum = 99999;
  for ( UInt_t b = 0; b < Nbins; ++b )
    {
      std::cout << "fabinmap[" << b << "].first=" << fabinmap[b].first << "  second=" << fabinmap[b].second << std::endl;
      if ( fabinmap[b].first != 0 )
        {
          FAratio->SetPoint( b, timeBinHist[b].bincenter, ( double )fabinmap[b].second / ( ( double )fabinmap[b].second+( double )fabinmap[b].first ) );
          Double_t sum = ( double )fabinmap[b].second+( double )fabinmap[b].first;
          sumFA->SetPoint( b, timeBinHist[b].bincenter, sum );
          if ( sum > maxsum ) maxsum = sum;
          if ( sum < minsum ) minsum = sum;
        }
    }
  FAratio->Draw( "lp" );
  //FAratio->GetXaxis()->SetTitle( "Drift Time ( us )" );
  //FAratio->GetYaxis()->SetTitle( "# Assumed / ( # Assumed + # Found )" );
  FAratio->SetMarkerStyle( 21 );
  FAratio->SetMarkerSize( 2 );
  //FAratio->SetMarkerColor( 2 );
  FAratio->SetLineWidth( 2 );
  canvFAratio->cd();
  TPad *overlay = new TPad( "overlay", "", 0, 0, 1, 1 );
  overlay->SetFillStyle( 4000 );
  overlay->SetFillColor( 0 );
  overlay->SetFrameFillStyle( 4000 );
  overlay->Draw();
  overlay->cd();
  Double_t xmin1 = 0;
  Double_t ymin1 = minsum*0.9;
  Double_t xmax1 = drifttimemax;
  Double_t ymax1 = maxsum*1.1;
  TH1F *hframe = overlay->DrawFrame( xmin1, ymin1, xmax1, ymax1 );
  hframe->GetXaxis()->SetLabelOffset( 99 );
  hframe->GetYaxis()->SetLabelOffset( 99 );
  sumFA->Draw( "lp" );
  sumFA->SetMarkerStyle( 21 );
  sumFA->SetMarkerSize( 2 );
  sumFA->SetMarkerColor( 4 );
  sumFA->SetLineWidth( 2 );
  TGaxis *axis1 = new TGaxis( xmax1, ymin1, xmax1, ymax1, ymin1, ymax1, 510, "+L" );
  axis1->SetTitle( "Total # of hits" );
  axis1->SetTitleOffset( 1.1 );
  axis1->SetLineColor( 4 );
  axis1->SetLabelColor( 4 );
  axis1->Draw();
  canvFAratio->Update();
  canvFAratio->Write();
  canvFAratio->SaveAs( "canvFAratio.png" );
  delete canvFAratio;

  //return;

  gStyle->SetOptStat( 11 );
  TCanvas * canvtest = new TCanvas( "fahist", "Found and Assumed Hits", 2000, 1600 );
  TPad * pad1 = new TPad( "pad1", "", 0, 0, 1, 1 );
  TPad * pad2 = new TPad( "pad2", "", 0, 0, 1, 1 );
  pad2->SetFillStyle( 4000 );
  pad1->Draw();
  pad1->cd();
  Double_t maxval = std::max( found->GetBinContent( found->GetMaximumBin() ), assumed->GetBinContent( assumed->GetMaximumBin() ) );
  found->SetLineColor( kRed );
  found->SetLineWidth( 2 );
  found->SetAxisRange( 0, maxval*1.1, "Y" );
  found->Draw();
  found->GetXaxis()->SetTitle( "Charge" );
  pad1->Update();
  TPaveStats * ps1 = ( TPaveStats* )found->GetListOfFunctions()->FindObject( "stats" );
  ps1->SetX1NDC( 0.4 ); ps1->SetX2NDC( 0.6 );
  ps1->SetTextColor( kRed );
  pad1->Modified();
  canvtest->cd();
  Double_t ymin = 0;
  Double_t ymax = maxval*1.1;
  Double_t dy = ( ymax-ymin )/0.8;
  Double_t xmin = ADCcutofflow;
  Double_t xmax = ADCcutoffhigh;
  Double_t dx = ( xmax-xmin )/0.8;
  pad2->Range( xmin-0.1*dx, ymin-0.1*dy, xmax+0.1*dx, ymax+0.1*dy );
  pad2->Draw();
  pad2->cd();
  assumed->SetLineColor( kBlue );
  assumed->SetLineWidth( 2 );
  assumed->Draw( "][sames" );
  assumed->GetXaxis()->SetTitle( "Charge" );
  pad2->Update();
  TPaveStats * ps2 = ( TPaveStats* )assumed->GetListOfFunctions()->FindObject( "stats" );
  ps2->SetX1NDC( 0.65 ); ps2->SetX2NDC( 0.85 );
  ps2->SetTextColor( kBlue );
  TGaxis * axis = new TGaxis( xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L" );
  axis->SetLabelFont( 42 );
  axis->SetLabelSize( 0.035 );
  axis->Draw();
  TLegend * leg = new TLegend( 0.55, 0.6, 0.85, 0.75 );
  leg->AddEntry( found, "Found Hits", "l" );
  leg->AddEntry( assumed, "Assumed Hits", "l" );
  leg->Draw();
  canvtest->Update();
  canvtest->Write();
  canvtest->SaveAs( "fahist.png" );
  delete canvtest;


  for ( auto & i_chan : hitmap )
    {
      for ( auto & i_hit : i_chan.second )
        {
          if ( i_hit.dqdx < ADCcutoffhigh && i_hit.dqdx > ADCcutofflow )
            {
              charge = ( Float_t )( i_hit.dqdx );
              tdrift = ( Float_t )( i_hit.hitt );
              weight = ( Float_t )( i_hit.weight );
              chgdata->add( RooArgSet( tdrift, charge ) /* , (Float_t)(i_hit.weight)*/ );
              for ( UInt_t i_hist = 0; i_hist < timeBinHist.size(); i_hist++ )
                {
                  if ( i_hit.hitt >= timeBinHist[i_hist].binlow && i_hit.hitt < timeBinHist[i_hist].binhigh )
                    {
                      timeBinHist[i_hist].hist->Fill( i_hit.dqdx );
                      timeBinHist[i_hist].histuw->Fill( i_hit.unscaleddqdx );
                      break;
                    }
                }
            }
        }
    }
  chgdata->Print();

  std::vector<FitResult> results( Nbins );

  Vavilov_Func vav;
  ROOT::Math::Functor func( vav, 6 );


  RooRealVar kappa( "VavilovKappa", "vavilov kappa", 0.1, 0, 10 );
  RooRealVar beta2( "VavilovBeta2", "vavilov beta2", 0.4, 0, 1 );
  RooRealVar mv( "VavilovMean", "vavilov mean", 2800, -1000, 20000 );
  RooRealVar sv( "VavilovSigma", "vavilov sigma", 1000, 0, 10000 );
  RooRealVar av( "VavilovAmplitude", "vavilov amp", 1e4, 1e0, 1e9 );
  RooAbsPdf * vavilov = new RooFunctorPdfBinding( "vavbind", "vavbind", func, RooArgList( charge, kappa, beta2, mv, sv, av ) );
  RooRealVar mgv( "mgv", "mean gauss vav", 0 );
  RooRealVar sgv( "sgv", "sigma gauss vav", 50, 0, 2000 );
  RooGaussian gaussvav( "gaussvav", "gaussvav", charge, mgv, sgv );
  RooFFTConvPdf vxg( "vxg", "vavilov ( x ) gauss", charge, *vavilov, gaussvav );
  RooPolynomial pol0_1( "pol0_1", "pol0", charge, RooArgList() );
  RooRealVar vxg_yield( "BkgCoef", "BkgCoef", 0.5, 0, 1 );
  RooAddPdf cvxg( "cvxg", "const + ( vavilov ( x ) gauss )", RooArgList( vxg, pol0_1 ), RooArgList( vxg_yield ) );

  Double_t mlscaleVal = 1000;
  Double_t slscaleVal = 100;
  Double_t sgscaleVal = 100;
  RooRealVar mlvar( "LandMPV/1000", "mean landau var", 2000/mlscaleVal, ADCcutoffhigh/mlscaleVal, ADCcutoffhigh/mlscaleVal );
  RooRealVar mlscale("LandMPVScale","scale mean landau", mlscaleVal);
  RooProduct ml("LandMPV","mlvar*mlscale",RooArgSet(mlvar,mlscale));
  RooRealVar slvar( "LandWidth/100", "sigma landau var", 200/slscaleVal, 60/slscaleVal, 2000/slscaleVal );
  RooRealVar slscale("LandWidthScale","scale sigma landau", slscaleVal);
  RooProduct sl("LandWidth","slvar*slscale",RooArgSet(slvar,slscale));
  RooLandau landau( "lx", "lx", charge, ml, sl );
  RooRealVar mg( "mg", "mean gauss", 0 );
  RooRealVar sgvar( "GaussWidth/100", "sigma gauss var", 200/sgscaleVal, 60/sgscaleVal, 2000/sgscaleVal );
  RooRealVar sgscale( "GaussWidthScale","scale sigma gauss",sgscaleVal);
  RooProduct sg("GaussWidth","sgvar*sgscale",RooArgSet(sgvar,sgscale));
  RooGaussian gauss( "gauss", "gauss", charge, mg, sg );
  RooFFTConvPdf lxg( "lxg", "landau ( x ) gauss", charge, landau, gauss );
  RooPolynomial pol0_2( "pol0_2", "pol0", charge, RooArgList() );
  RooRealVar lxg_yield( "BkgCoef", "BkgCoef", 0 );
  RooAddPdf clxg( "clxg", "const + ( landau ( x ) gauss )", RooArgList( pol0_2, lxg ), RooArgList( lxg_yield ) );

  RooRealVar pgmean( "PeakGaussMean", "peak gauss mean", 2800, -1000, 200000 );
  RooRealVar pgwidth( "PeakGaussWidth", "peak gauss width", 1000, 1, 200000 );
  RooGaussian pg( "PeakGaussian", "peak gaussian", charge, pgmean, pgwidth );
  RooPolynomial pol0_3( "pol0_3", "pol0", charge, RooArgList() );
  RooRealVar pg_yield( "BkgCoef", "BkgCoef", 0.5, 0, 1 );
  RooAddPdf cpg( "cpg", "const + peak gaussian", RooArgList( pol0_3, pg ), RooArgList( pg_yield ) );

  TGraph * hbin_lxg = new TGraph( Nbins );
  hbin_lxg->SetTitle( "Most Probable dQ/dx in L( x )g Fit" );
  Double_t minMPVlxg = 99999;
  Double_t maxMPVlxg = -99999;

  TGraphErrors * hbin_vxg = new TGraphErrors( Nbins );
  hbin_vxg->SetTitle( "Most Probable dQ/dx in V( x )g Fit" );
  Double_t minMPVvxg = 99999;
  Double_t maxMPVvxg = -99999;

  TGraphErrors * hbin_pg = new TGraphErrors( Nbins );
  hbin_pg->SetTitle( "Peak dQ/dx in PeakGauss Fit" );
  Double_t minMPVpg = 99999;
  Double_t maxMPVpg = -99999;

  TGraphErrors * gw_lxg = new TGraphErrors( Nbins );
  gw_lxg->SetTitle( "Width of Gaussian in L( x )g Convolution" );
  TGraphErrors * gw_vxg = new TGraphErrors( Nbins );
  gw_vxg->SetTitle( "Width of Gaussian in V( x )g Convolution" );
  TGraphErrors * lw = new TGraphErrors( Nbins );
  lw->SetTitle( "Width of Landau in L( x )g Convolution" );
  TGraphErrors * vw = new TGraphErrors( Nbins );
  vw->SetTitle( "Width of Vavilov in V( x )g Convolution" );
  TGraphErrors * kv = new TGraphErrors( Nbins );
  kv->SetTitle( "kappa of Vavilov in V( x )g Convolution" );
  TGraphErrors * bv = new TGraphErrors( Nbins );
  bv->SetTitle( "beta2 of Vavilov in V( x )g Convolution" );
  TGraphErrors * gw_pg = new TGraphErrors( Nbins );
  gw_pg->SetTitle( "Width of Gaussian Fit to Peak" );

  TGraph * lxgchi2 = new TGraph( Nbins );
  lxgchi2->SetTitle( "Chi2/NDF of L( x )g Best Fit" );

  for ( UInt_t i = 0; i < Nbins; i++ )
    {
      if ( timeBinHist[i].hist->GetEntries() > 0 )
        {
          char shortname[100];
          sprintf( shortname, "dQdx_bin%d", i );

          char cut[100];
          sprintf( cut, "tdrift > %f && tdrift < %f", timeBinHist[i].binlow, timeBinHist[i].binhigh );
          RooDataSet* roodata = ( RooDataSet* ) chgdata->reduce( RooFit::SelectVars( RooArgSet(charge) ), RooFit::Cut( cut ), RooFit::Name( shortname ), RooFit::Title( shortname ) );

          TH1D * clonehist = (TH1D*)roodata->createHistogram("roodata_clone",charge,RooFit::Binning(400,ADCcutofflow,ADCcutoffhigh));

          Float_t maxbin = clonehist->GetBinCenter( clonehist->GetMaximumBin() );
          //Float_t maxbin = 2000;

          Float_t maxvalue = clonehist->GetBinContent( clonehist->GetMaximumBin() );
          Float_t cutoffheightlow = maxvalue*0.66;
          Float_t cutoffheighthigh = maxvalue*0.66;
          Float_t lowval = clonehist->GetBinCenter( clonehist->FindFirstBinAbove( cutoffheightlow ) );
          Float_t highval = clonehist->GetBinCenter( clonehist->FindLastBinAbove( cutoffheighthigh ) );

          std::cout << "maxbin=" << maxbin << "  maxvalue=" << maxvalue << "  cutoffheight=( " << cutoffheightlow << ", " << cutoffheighthigh << " )  lowval=" << lowval << "  highval=" << highval << std::endl;

/*
          TCanvas * canvtest = new TCanvas("canvtest","",1600,1600);
          clonehist->Draw();
          clonehist->SetLineColor(kRed);
          timeBinHist[i].histuw->Draw("same");
          timeBinHist[i].histuw->SetLineColor(kBlue);
          canvtest->Update();
          //std::cin.get();
 */

          pgmean.setVal( maxbin );
          pgmean.setMin( lowval );
          pgmean.setMax( highval );
          pgwidth.setVal( 1000 );
          pgwidth.setMin( 1 );
          pgwidth.setMax( 10000 );

          Float_t lxgfitlow = maxvalue*0.3;
          Float_t lxgfithigh = maxvalue*0.5;
          Float_t lxglowval = clonehist->GetBinCenter( clonehist->FindFirstBinAbove( lxgfitlow ) );
          Float_t lxghighval = clonehist->GetBinCenter( clonehist->FindLastBinAbove( lxgfithigh ) );

          // Setup component pdfs
          // --------------------

          //setup observable
          //charge.setRange( shortname, maxbin-2000, maxbin+3000 );
          charge.setRange( shortname, ADCcutofflow, ADCcutoffhigh );

          //setup landau( t, ml, sl )
          mlvar.setVal( maxbin/mlscaleVal );
          mlvar.setMin( (maxbin-2000)/mlscaleVal ); //0 );
          mlvar.setMax( (maxbin+2000)/mlscaleVal ); //0 );
          slvar.setVal( 500/slscaleVal );
          //sl.setMin( 10 );
          //sl.setMax( 2000 );

          //setup gauss( t, mg, sg )
          //mg.setVal( 1 );
          //mg.setMin( -100 );
          //mg.setMax( 100 );
          sgvar.setVal( 600/slscaleVal );
          //sg.setMin( 10 );
          //sg.setMax( 2000 );

          // construct convolution pdf
          // -------------------------

          // set num bins to be used for FFT sampling
          //charge.setBins( 10000, "fft" );

          // fit convoluted pdf to binHist
          // -----------------------------



          if ( DoLxG )
            {
              RooFitResult * fr = lxg.fitTo( *roodata, RooFit::Save(), RooFit::PrintLevel( 1 ), RooFit::Range( lxglowval, lxghighval ), /*RooFit::SumW2Error( kTRUE ),*/ RooFit::Minimizer("Minuit2"), RooFit::Minos(kTRUE), RooFit::Offset(kTRUE) );
              RooPlot * lxgframe = charge.frame( RooFit::Title( TString::Format( "%s L( x )g", timeBinHist[i].hist->GetTitle() ) ) );
              RooDataHist * binhistuw = new RooDataHist(timeBinHist[i].histuw->GetName(),timeBinHist[i].histuw->GetTitle(),RooArgList(charge),timeBinHist[i].histuw);
              binhistuw->plotOn( lxgframe, RooFit::DrawOption("e2"), RooFit::FillColor(kBlue), RooFit::FillStyle(3001), RooFit::LineColor(kBlack), RooFit::MarkerColor(kBlack) );
              roodata->plotOn( lxgframe, RooFit::Binning( 200 ), RooFit::DataError(RooAbsData::SumW2), RooFit::LineColor(kBlack), RooFit::MarkerColor(kBlack) );
              lxg.plotOn( lxgframe, RooFit::Name( "lxg" ), RooFit::LineColor( kRed ), RooFit::Range( lxglowval, lxghighval ) );
              results[i].lxg_chi2ndf = lxgframe->chiSquare();

              std::cout << "status = " << fr->status() << "  covQual = " << fr->covQual() << "  numInvalidNLL = " << fr->numInvalidNLL() << "  EDM = " << fr->edm() << "  minNll = " << fr->minNll() << std::endl;

              TPaveLabel * tlxg = new TPaveLabel( 0.7, 0.83, 0.99, 0.9, Form( " #chi^{2}/NDF = %f", lxgframe->chiSquare() ), "NDC" );
              tlxg->SetFillColor( 0 );
              tlxg->SetBorderSize( 1 );
              tlxg->SetTextAlign( 12 );
              lxgframe->addObject( tlxg );
              lxgframe->getAttText()->SetTextSize( 0.30 );
              roodata->statOn( lxgframe, RooFit::Layout( 0.7, 0.99, 0.8 ) );
              lxgframe->getAttText()->SetTextSize( 0.019 );
              lxg.paramOn( lxgframe, RooFit::Layout( 0.7, 0.99, 0.6 ) );
              lxgframe->getAttText()->SetTextSize( 0.019 );

              char buff[100];
              sprintf( buff, "lxg_bin%d", i );
              TCanvas* canvlxg2 = new TCanvas( buff, buff, 1600, 1600 );
              gPad->SetLeftMargin( 0.15 );
              lxgframe->GetYaxis()->SetTitleOffset( 1.4 );
              lxgframe->Draw();
              canvlxg2->Write();
              canvlxg2->SaveAs( TString::Format( "%s.png", buff ) );
              /*
                 gSystem->ProcessEvents();
                 TImage *imglxg = TImage::Create();
                 imglxg->FromPad( canvlxg2 );
                 char newbuff[100];
                 sprintf( newbuff, "%s.png", buff );
                 imglxg->WriteImage( newbuff );
               */

              if ( DoResidual )
                {
                  char buffresid[100];
                  sprintf( buffresid, "bin%d_resid_lxg", i );
                  TCanvas* cresidlxg = new TCanvas( buffresid, buffresid, 1600, 1600 );
                  gPad->SetLeftMargin( 0.15 );
                  RooPlot* frameresidlxg = charge.frame( RooFit::Title( "lxg_residual" ) );
                  RooHist* hresidlxg = lxgframe->residHist();
                  frameresidlxg->addPlotable( hresidlxg, "P" );
                  frameresidlxg->GetYaxis()->SetTitleOffset( 1.4 );
                  frameresidlxg->Draw();
                  cresidlxg->Write();
                  gSystem->ProcessEvents();

                  delete frameresidlxg;
                  delete cresidlxg;
                }

              if ( DoPull )
                {
                  char buffpull[100];
                  sprintf( buffpull, "bin%d_pull_lxg", i );
                  TCanvas* cpulllxg = new TCanvas( buffpull, buffpull, 1600, 1600 );
                  gPad->SetLeftMargin( 0.15 );
                  RooPlot* framepulllxg = charge.frame( RooFit::Title( "lxg_pull" ) );
                  RooHist* hpull2lxg = lxgframe->pullHist();
                  TH1F* hpulllxg = new TH1F( "hpull_lxg", "pull_lxg", 100, 0, 0 );
                  for ( Int_t ibin=0; ibin<hpull2lxg->GetN(); ibin++ )
                    {
                      Double_t x, y;
                      hpull2lxg->GetPoint( ibin, x, y );
                      hpulllxg->Fill( y );
                    }
                  hpulllxg->Draw();
                  cpulllxg->Write();
                  gSystem->ProcessEvents();

                  delete cpulllxg;
                  delete hpulllxg;
                  delete framepulllxg;
                }

              results[i].landmpv = ml.getVal()/mcscale;
              results[i].landwidth = sl.getVal()/mcscale;
              results[i].gaussmean_l = mg.getVal()/mcscale;
              results[i].gausswidth_l = sg.getVal()/mcscale;
              results[i].landmpverr = (mlvar.getError()/mcscale)*mlscaleVal;
              results[i].landwidtherr = (slvar.getError()/mcscale)*slscaleVal;
              results[i].gaussmean_lerr = mg.getError()/mcscale;
              results[i].gausswidth_lerr = (sgvar.getError()/mcscale)*sgscaleVal;
              results[i].lxglow = lxglowval/mcscale;
              results[i].lxghigh = lxghighval/mcscale;
              //results[i].lxg_chi2ndf = lxgframe->chiSquare( "lxg", shortname );

              if ( minMPVlxg > results[i].landmpv ) minMPVlxg = results[i].landmpv;
              if ( maxMPVlxg < results[i].landmpv ) maxMPVlxg = results[i].landmpv;

              results[i].fitsuccess = false;
              if ( fr->covQual() == 3 && fr->edm() < 1 && results[i].lxg_chi2ndf < 3)
                {
                  results[i].fitsuccess = true;
                  hbin_lxg->SetPoint( i, timeBinHist[i].bincenter, results[i].landmpv );
                  //hbin_lxg->SetPointError( i, 0, results[i].landmpverr );
                }
              gw_lxg->SetPoint( i, timeBinHist[i].bincenter, results[i].gausswidth_l );
              gw_lxg->SetPointError( i, 0, results[i].gausswidth_lerr );
              lw->SetPoint( i, timeBinHist[i].bincenter, results[i].landwidth );
              lw->SetPointError( i, 0, results[i].landwidtherr );

              lxgchi2->SetPoint( i, timeBinHist[i].bincenter, results[i].lxg_chi2ndf );

              delete canvlxg2;
              //delete imglxg;
              delete lxgframe;
            }

          if ( DoVxG )
            {
              cvxg.fitTo( *roodata, RooFit::Save(), RooFit::PrintLevel( -1 ), RooFit::Range( lxglowval, lxghighval ) );
              RooPlot * vxgframe = charge.frame( RooFit::Title( TString::Format( "%s V( x )g", timeBinHist[i].hist->GetTitle() ) ) );
              roodata->plotOn( vxgframe, RooFit::Binning( 200 ) );
              cvxg.plotOn( vxgframe, RooFit::LineColor( kGreen ) /*, RooFit::ShiftToZero()*/, RooFit::Name( "vxg" ), RooFit::Range( lxglowval, lxghighval ) );
              //results[i].vxg_chi2ndf = vxgframe->chiSquare();

              TPaveLabel * tvxg = new TPaveLabel( 0.7, 0.83, 0.99, 0.9, Form( " #chi^{2}/NDF = %f", vxgframe->chiSquare() ), "NDC" );
              tvxg->SetFillColor( 0 );
              tvxg->SetBorderSize( 1 );
              tvxg->SetTextAlign( 12 );
              tvxg->SetTextSize( 0.019 );
              vxgframe->addObject( tvxg );
              vxgframe->getAttText()->SetTextSize( 0.3 );
              roodata->statOn( vxgframe, RooFit::Layout( 0.7, 0.99, 0.8 ) );
              vxgframe->getAttText()->SetTextSize( 0.019 );
              cvxg.paramOn( vxgframe, RooFit::Layout( 0.7, 0.99, 0.6 ) );
              vxgframe->getAttText()->SetTextSize( 0.019 );

              char buff[100];
              sprintf( buff, "vxg_bin%d", i );
              TCanvas* canvvxg2 = new TCanvas( buff, buff, 1600, 1600 );
              gPad->SetLeftMargin( 0.15 );
              vxgframe->GetYaxis()->SetTitleOffset( 1.4 );
              vxgframe->Draw();
              canvvxg2->Write();
              gSystem->ProcessEvents();
              TImage *imgvxg = TImage::Create();
              imgvxg->FromPad( canvvxg2 );
              char newbuff[100];
              sprintf( newbuff, "%s.png", buff );
              imgvxg->WriteImage( newbuff );

              if ( DoResidual )
                {
                  char buffresid[100];
                  sprintf( buffresid, "bin%d_resid_vxg", i );
                  TCanvas* cresidvxg = new TCanvas( buffresid, buffresid, 1600, 1600 );
                  gPad->SetLeftMargin( 0.15 );
                  RooPlot* frameresidvxg = charge.frame( RooFit::Title( "vxg_residual" ) );
                  RooHist* hresidvxg = vxgframe->residHist();
                  frameresidvxg->addPlotable( hresidvxg, "P" );
                  frameresidvxg->GetYaxis()->SetTitleOffset( 1.4 );
                  frameresidvxg->Draw();
                  cresidvxg->Write();
                  gSystem->ProcessEvents();

                  delete cresidvxg;
                  delete frameresidvxg;
                }

              if ( DoPull )
                {
                  char buffpull[100];
                  sprintf( buffpull, "bin%d_pull_vxg", i );
                  TCanvas* cpullvxg = new TCanvas( buffpull, buffpull, 1600, 1600 );
                  gPad->SetLeftMargin( 0.15 );
                  RooPlot* framepullvxg = charge.frame( RooFit::Title( "vxg_pull" ) );
                  RooHist* hpull2vxg = vxgframe->pullHist();
                  TH1F* hpullvxg = new TH1F( "hpull_vxg", "pull_vxg", 100, 0, 0 );
                  for ( Int_t ibin=0; ibin<hpull2vxg->GetN(); ibin++ )
                    {
                      Double_t x, y;
                      hpull2vxg->GetPoint( ibin, x, y );
                      hpullvxg->Fill( y );
                    }
                  hpullvxg->Draw();
                  cpullvxg->Write();
                  gSystem->ProcessEvents();

                  delete cpullvxg;
                  delete hpullvxg;
                  delete framepullvxg;
                }

              results[i].gaussmean_v = mgv.getVal();
              results[i].gausswidth_v = sgv.getVal();
              results[i].vavkappa = kappa.getVal();
              results[i].vavbeta2 = beta2.getVal();
              results[i].vavmean = mv.getVal();
              results[i].vavwidth = sv.getVal();
              results[i].vavamp = av.getVal();
              results[i].gaussmean_verr = mgv.getError();
              results[i].gausswidth_verr = sgv.getError();
              results[i].vavkappaerr = kappa.getError();
              results[i].vavbeta2err = beta2.getError();
              results[i].vavmeanerr = mv.getError();
              results[i].vavwidtherr = sv.getError();
              results[i].vavamperr = av.getError();
              results[i].vxg_chi2ndf = vxgframe->chiSquare( "cvxg", shortname );

              if ( minMPVvxg > results[i].vavmean ) minMPVvxg = results[i].vavmean;
              if ( maxMPVvxg < results[i].vavmean ) maxMPVvxg = results[i].vavmean;

              //if ( results[i].vxg_chi2ndf < 2 )
              {
                hbin_vxg->SetPoint( i, timeBinHist[i].bincenter, results[i].vavmean );
                hbin_vxg->SetPointError( i, 0, results[i].vavmeanerr );
              }
              gw_vxg->SetPoint( i, timeBinHist[i].bincenter, results[i].gausswidth_v );
              gw_vxg->SetPointError( i, 0, results[i].gausswidth_verr );
              vw->SetPoint( i, timeBinHist[i].bincenter, results[i].vavwidth );
              vw->SetPointError( i, 0, results[i].vavwidtherr );
              kv->SetPoint( i, timeBinHist[i].bincenter, results[i].vavkappa );
              kv->SetPointError( i, 0, results[i].vavkappaerr );
              bv->SetPoint( i, timeBinHist[i].bincenter, results[i].vavbeta2 );
              bv->SetPointError( i, 0, results[i].vavbeta2err );

              delete canvvxg2;
              delete vxgframe;
              delete imgvxg;
            }

          if ( DoPG )
            {
              cpg.fitTo( *roodata, RooFit::Save(), RooFit::PrintLevel( -1 ), RooFit::Range( lowval, highval ), RooFit::SumW2Error( true ) );
              RooPlot * pgframe = charge.frame( RooFit::Title( TString::Format( "%s PeakGaussian", timeBinHist[i].hist->GetTitle() ) ) );
              roodata->plotOn( pgframe, RooFit::Binning( 200 ) );
              cpg.plotOn( pgframe, RooFit::Name( "pgb" ), RooFit::LineColor( kBlue ), RooFit::Range( lowval, highval ) /*, RooFit::ShiftToZero() , RooFit::VLines()*/ );
              //results[i].pg_chi2ndf = pgframe->chiSquare();

              TPaveLabel * tpg = new TPaveLabel( 0.7, 0.83, 0.99, 0.9, Form( " #chi^{2}/NDF = %f", pgframe->chiSquare() ), "NDC" );
              tpg->SetFillColor( 0 );
              tpg->SetBorderSize( 1 );
              tpg->SetTextAlign( 12 );
              tpg->SetTextSize( 0.3 );
              pgframe->addObject( tpg );
              pgframe->getAttText()->SetTextSize( 0.3 );
              roodata->statOn( pgframe, RooFit::Layout( 0.7, 0.99, 0.8 ) );
              pgframe->getAttText()->SetTextSize( 0.019 );
              cpg.paramOn( pgframe, RooFit::Layout( 0.7, 0.99, 0.6 ) );
              pgframe->getAttText()->SetTextSize( 0.019 );

              char buff[100];
              sprintf( buff, "pg_bin%d", i );
              TCanvas* canvpg2 = new TCanvas( buff, buff, 1600, 1600 );
              gPad->SetLeftMargin( 0.15 );
              pgframe->GetYaxis()->SetTitleOffset( 1.4 );
              pgframe->Draw();
              canvpg2->Write();
              canvpg2->SaveAs( TString::Format( "%s.png", buff ) );
              /*
                 gSystem->ProcessEvents();
                 TImage *imgpg = TImage::Create();
                 imgpg->FromPad( canvpg2 );
                 char newbuff[100];
                 sprintf( newbuff, "%s.png", buff );
                 imgpg->WriteImage( newbuff );
               */

              if ( DoResidual )
                {
                  char buffresid[100];
                  sprintf( buffresid, "bin%d_resid_pg", i );
                  TCanvas* cresidpg = new TCanvas( buffresid, buffresid, 1600, 1600 );
                  gPad->SetLeftMargin( 0.15 );
                  RooPlot* frameresidpg = charge.frame( RooFit::Title( "pg_residual" ) );
                  RooHist* hresidpg = pgframe->residHist();
                  frameresidpg->addPlotable( hresidpg, "P" );
                  frameresidpg->GetYaxis()->SetTitleOffset( 1.4 );
                  frameresidpg->Draw();
                  cresidpg->Write();
                  gSystem->ProcessEvents();

                  delete cresidpg;
                  delete frameresidpg;

                }

              if ( DoPull )
                {
                  char buffpull[100];
                  sprintf( buffpull, "bin%d_pull_pg", i );
                  TCanvas* cpullpg = new TCanvas( buffpull, buffpull, 1600, 1600 );
                  gPad->SetLeftMargin( 0.15 );
                  RooPlot* framepullpg = charge.frame( RooFit::Title( "pg_pull" ) );
                  RooHist* hpull2pg = pgframe->pullHist();
                  TH1F* hpullpg = new TH1F( "hpull_pg", "pull_pg", 100, 0, 0 );
                  for ( Int_t ibin=0; ibin<hpull2pg->GetN(); ibin++ )
                    {
                      Double_t x, y;
                      hpull2pg->GetPoint( ibin, x, y );
                      hpullpg->Fill( y );
                    }
                  hpullpg->Draw();
                  cpullpg->Write();
                  gSystem->ProcessEvents();

                  delete cpullpg;
                  delete hpullpg;
                  delete framepullpg;
                }

              results[i].gaussmean = pgmean.getVal();
              results[i].gausswidth = pgwidth.getVal();
              results[i].gaussmeanerr = pgmean.getError();
              results[i].gausswidtherr = pgwidth.getError();
              results[i].pg_chi2ndf = pgframe->chiSquare( "cpg", shortname );

              if ( minMPVpg > results[i].gaussmean ) minMPVpg = results[i].gaussmean;
              if ( maxMPVpg < results[i].gaussmean ) maxMPVpg = results[i].gaussmean;

              //if ( results[i].pg_chi2ndf < 2 )
              {
                hbin_pg->SetPoint( i, timeBinHist[i].bincenter, results[i].gaussmean );
                hbin_pg->SetPointError( i, 0, results[i].gaussmeanerr );
              }
              gw_pg->SetPoint( i, timeBinHist[i].bincenter, results[i].gausswidth );
              gw_pg->SetPointError( i, 0, results[i].gausswidtherr );

              delete canvpg2;
              delete pgframe;
              //delete imgpg;
            }

          results[i].binnum = i;
          results[i].bincenter = timeBinHist[i].bincenter;
          results[i].binwidtherr = ( timeBinHist[i].binwidth )/sqrt( 12.0 );

        }
    }

  if ( DoLxG )
    {
      std::cout << "Landau( x )Gauss: " << std::endl;
      for ( UInt_t i = 0; i < Nbins; i++ )
        {
          std::cout << "     Bin " << i << ": ml=" << results[i].landmpv << " ( +/- " << results[i].landmpverr << " )  sl=" << results[i].landwidth << " ( +/- " << results[i].landwidtherr << " )  sg_l=" << results[i].gausswidth_l << " ( +/- " << results[i].gausswidth_lerr << " )  chi2ndf=" << results[i].lxg_chi2ndf << std::endl;
        }
    }
  if ( DoVxG )
    {
      std::cout << "Vavilov( x )Gauss: " << std::endl;
      for ( UInt_t i = 0; i < Nbins; i++ )
        {
          std::cout << "     Bin " << i << ": mv=" << results[i].vavmean << " ( +/- " << results[i].vavmeanerr << " )  sv=" << results[i].vavwidth << " ( +/- " << results[i].vavwidtherr << " )  sg_v=" << results[i].gausswidth_v << " ( +/- " << results[i].gausswidth_verr << " )  chi2ndf=" << results[i].vxg_chi2ndf << std::endl;
        }
    }
  if ( DoPG )
    {
      std::cout << "PeakFit Gauss: " << std::endl;
      for ( UInt_t i = 0; i < Nbins; i++ )
        {
          std::cout << "     Bin " << i << ": mean=" << results[i].gaussmean << " ( +/- " << results[i].gaussmeanerr << " )  width=" << results[i].gausswidth << " ( +/- " << results[i].gausswidtherr << " )  chi2ndf=" << results[i].pg_chi2ndf << std::endl;
        }
    }

  TF1 * expo = new TF1( "expo", "[0]*exp( -x/[1] )", secondbinstart, penultimatebinend );
  expo->SetParNames( "dQdx0", "eLifetime" );
  expo->SetParameters( 3000, 3000 );

  gStyle->SetOptStat( 0 );
  gStyle->SetOptFit( 1 );
  Int_t trans_red = TColor::GetColorTransparent( kRed, 0.35 );

  if ( DoLxG )
    {
      TCanvas * canvlxg = new TCanvas( "canvMPV", "LandauMPV", 1600, 900 );
      canvlxg->cd();
      hbin_lxg->Fit( "expo", "R" );
      hbin_lxg->Draw( "ap" );
      hbin_lxg->SetMarkerStyle( kFullDotLarge );
      hbin_lxg->SetMarkerSize( 2 );

      TH1F * hint_lxg = new TH1F( "hint_lxg", "68#% confidence band", 100, 0, drifttimemax );
      ( TVirtualFitter::GetFitter() )->GetConfidenceIntervals( hint_lxg, 0.68 );
      hint_lxg->SetStats( kFALSE );
      hint_lxg->SetFillColor( trans_red );
      //hint_lxg->Draw( "e3 same" );

      hbin_lxg->GetYaxis()->SetRangeUser( minMPVlxg-100, maxMPVlxg+100 );
      hbin_lxg->GetXaxis()->SetTitle( "Drift Time ( #mu s )" );
      hbin_lxg->GetYaxis()->SetTitle( "dQ/dx ( ADC/cm )" );
      canvlxg->Update();
      canvlxg->Write();
      canvlxg->SaveAs( "canvMPV.png" );

      TGraphErrors * lifetimedrift = new TGraphErrors();
      lifetimedrift->SetTitle( "Fitted Electron Lifetime vs. Fiducial Drift Cut" );
      lifetimedrift->SetMarkerStyle( 20 );
      lifetimedrift->SetMarkerColor( kBlue );
      lifetimedrift->GetXaxis()->SetTitle( "Drift Time ( us )" );
      lifetimedrift->GetYaxis()->SetTitle( "eLifetime ( us )" );
      TF1 * expotrunc = new TF1( "expotrunc", "[0]*exp( -x/[1] )", 0, 1 );
      expotrunc->SetParNames( "dQdx0", "eLifetime" );
      expotrunc->SetParameters( 3000, 3000 );
      UInt_t nbinwidth = 4;
      for ( UInt_t i_bin = 1; i_bin < Nbins-nbinwidth; ++i_bin )
        {
          expotrunc->SetRange( timeBinHist[i_bin].binlow, timeBinHist[i_bin+nbinwidth].binhigh );
          hbin_lxg->Fit( "expotrunc", "RQ" );
          lifetimedrift->SetPoint( i_bin-1, timeBinHist[i_bin+nbinwidth].bincenter, expotrunc->GetParameter( 1 ) );
          lifetimedrift->SetPointError( i_bin-1, 0, expotrunc->GetParError( 1 ) );
        }
      TCanvas * canvlifetimedrift_lxg = new TCanvas( "cltd_lxg", "cltd_lxg", 1600, 900 );
      canvlifetimedrift_lxg->cd();
      lifetimedrift->Draw( "alpe" );
      lifetimedrift->GetXaxis()->SetTitle( "Drift Distance Cut ( cm )" );
      lifetimedrift->GetYaxis()->SetTitle( "Elifetime ( us )" );
      canvlifetimedrift_lxg->Update();
      canvlifetimedrift_lxg->Write();
      canvlifetimedrift_lxg->SaveAs( "cltd_lxg.png" );

      TCanvas * canvgw_lxg = new TCanvas( "canvgw_lxg", "GaussWidth of Lxg", 1000, 800 );
      canvgw_lxg->cd();
      gw_lxg->Draw( "ap" );
      gw_lxg->SetMarkerStyle( kFullDotLarge );
      gw_lxg->SetMarkerSize( 2 );
      gw_lxg->GetXaxis()->SetTitle( "Drift Time ( #mu s )" );
      gw_lxg->GetYaxis()->SetTitle( "Gauss Width ( ticks )" );
      canvgw_lxg->Update();
      canvgw_lxg->Write();
      canvgw_lxg->SaveAs( "canvgw_lxg.png" );

      TCanvas * canvlw = new TCanvas( "canvlw", "LandauWidth", 1600, 900 );
      canvlw->cd();
      lw->Draw( "ap" );
      lw->SetMarkerStyle( kFullDotLarge );
      lw->SetMarkerSize( 2 );
      lw->GetXaxis()->SetTitle( "Drift Time ( #mu s )" );
      lw->GetYaxis()->SetTitle( "Landau Width ( ticks )" );
      canvlw->Update();
      canvlw->Write();
      canvlw->SaveAs( "canvlw.png" );

      TCanvas * canvlxgchi2 = new TCanvas( "canvlxgchi2", "LxG Chi2/NDF", 1000, 800 );
      canvlxgchi2->cd();
      lxgchi2->Draw( "ap" );
      lxgchi2->SetMarkerStyle( kFullDotLarge );
      lxgchi2->SetMarkerSize( 2 );
      lxgchi2->GetXaxis()->SetTitle( "Drift Time ( #mu s )" );
      lxgchi2->GetYaxis()->SetTitle( "Chi2/NDF of L( x )g Fit" );
      canvlxgchi2->Update();
      canvlxgchi2->Write();
      canvlxgchi2->SaveAs( "canvlxgchi2.png" );

    }

  if ( DoVxG )
    {
      TCanvas * canvvxg = new TCanvas( "canvMPVvxg", "VavilovMPV", 1600, 900 );
      canvvxg->cd();
      hbin_vxg->Fit( "expo", "R" );
      hbin_vxg->Draw( "ap" );
      hbin_vxg->SetMarkerStyle( kFullDotLarge );
      hbin_vxg->SetMarkerSize( 2 );

      TH1F * hint_vxg = new TH1F( "hint_vxg", "95#% confidence band", 100, 0, drifttimemax );
      ( TVirtualFitter::GetFitter() )->GetConfidenceIntervals( hint_vxg, 0.68 );
      hint_vxg->SetStats( kFALSE );
      hint_vxg->SetFillColor( trans_red );
      hint_vxg->Draw( "e3 same" );

      hbin_vxg->GetYaxis()->SetRangeUser( minMPVvxg-100, maxMPVvxg+100 );
      hbin_vxg->GetXaxis()->SetTitle( "Drift Time ( #mu s )" );
      hbin_vxg->GetYaxis()->SetTitle( "dQ/dx ( ADC/cm )" );
      canvvxg->Update();
      canvvxg->Write();
      canvvxg->SaveAs( "canvMPVvxg.png" );

      TCanvas * canvgw_vxg = new TCanvas( "canvgw_vxg", "GaussWidth of Vxg", 1600, 900 );
      canvgw_vxg->cd();
      gw_vxg->Draw( "ap" );
      gw_vxg->SetMarkerStyle( kFullDotLarge );
      gw_vxg->SetMarkerSize( 2 );
      gw_vxg->GetXaxis()->SetTitle( "Drift Time ( #mu s )" );
      gw_vxg->GetYaxis()->SetTitle( "Gauss Width ( ticks )" );
      canvgw_vxg->Update();
      canvgw_vxg->Write();
      canvgw_vxg->SaveAs( "canvgw_vxg.png" );

      TCanvas * canvvw = new TCanvas( "canvvw", "VavilovWidth", 1600, 900 );
      canvvw->cd();
      vw->Draw( "ap" );
      vw->SetMarkerStyle( kFullDotLarge );
      vw->SetMarkerSize( 2 );
      vw->GetXaxis()->SetTitle( "Drift Time ( #mu s )" );
      vw->GetYaxis()->SetTitle( "Vavilov Width ( ticks )" );
      canvvw->Update();
      canvvw->Write();
      canvvw->SaveAs( "canvvw.png" );

      TCanvas * canvkv = new TCanvas( "canvkv", "VavilovKappa", 1600, 900 );
      canvkv->cd();
      kv->Draw( "ap" );
      kv->SetMarkerStyle( kFullDotLarge );
      kv->SetMarkerSize( 2 );
      kv->GetXaxis()->SetTitle( "Drift Time ( #mu s )" );
      kv->GetYaxis()->SetTitle( "Vavilov Kappa" );
      canvkv->Update();
      canvkv->Write();
      canvkv->SaveAs( "canvkv.png" );

      TCanvas * canvbv = new TCanvas( "canvbv", "VavilovBeta2", 1600, 900 );
      canvbv->cd();
      bv->Draw( "ap" );
      bv->SetMarkerStyle( kFullDotLarge );
      bv->SetMarkerSize( 2 );
      bv->GetXaxis()->SetTitle( "Drift Time ( #mu s )" );
      bv->GetYaxis()->SetTitle( "Vavilov Beta2" );
      canvbv->Update();
      canvbv->Write();
      canvbv->SaveAs( "canvbv.png" );
    }

  if ( DoPG )
    {
      TCanvas * canvpg = new TCanvas( "canvpg", "PeakGaussMean", 1600, 900 );
      canvpg->cd();
      hbin_pg->Fit( "expo", "R" );
      hbin_pg->Draw( "ap" );
      hbin_pg->SetMarkerStyle( kFullDotLarge );
      hbin_pg->SetMarkerSize( 2 );

      TH1F * hint_pg = new TH1F( "hint_pg", "95#% confidence band", 100, 0, drifttimemax );
      ( TVirtualFitter::GetFitter() )->GetConfidenceIntervals( hint_pg, 0.68 );
      hint_pg->SetStats( kFALSE );
      hint_pg->SetFillColor( trans_red );
      hint_pg->Draw( "e3 same" );

      hbin_pg->GetYaxis()->SetRangeUser( minMPVpg-100, maxMPVpg+100 );
      hbin_pg->GetXaxis()->SetTitle( "Drift Time ( #mu s )" );
      hbin_pg->GetYaxis()->SetTitle( "dQ/dx ( ADC/cm )" );
      canvpg->Update();
      canvpg->Write();
      canvpg->SaveAs( "canvpg.png" );

      TCanvas * canvgw_pg = new TCanvas( "canvgw_pg", "PeakGaussWidth", 1600, 900 );
      canvgw_pg->cd();
      gw_pg->Draw( "ap" );
      gw_pg->SetMarkerStyle( kFullDotLarge );
      gw_pg->SetMarkerSize( 2 );
      gw_pg->GetXaxis()->SetTitle( "Drift Time ( #mu s )" );
      gw_pg->GetYaxis()->SetTitle( "Gauss Width ( ticks )" );
      canvgw_pg->Update();
      canvgw_pg->Write();
      canvgw_pg->SaveAs( "canvgw_pg.png" );

      TGraphErrors * lifetimedrift = new TGraphErrors();
      lifetimedrift->SetTitle( "Fitted Electron Lifetime vs. Fiducial Drift Cut" );
      lifetimedrift->SetMarkerStyle( 20 );
      lifetimedrift->SetMarkerColor( kBlue );
      lifetimedrift->GetXaxis()->SetTitle( "Drift Time ( us )" );
      lifetimedrift->GetYaxis()->SetTitle( "eLifetime ( us )" );
      TF1 * expotrunc = new TF1( "expotrunc", "[0]*exp( -x/[1] )", 0, 1 );
      expotrunc->SetParNames( "dQdx0", "eLifetime" );
      expotrunc->SetParameters( 3000, 3000 );
      UInt_t nbinwidth=4;
      for ( UInt_t i_bin = 1; i_bin < Nbins-nbinwidth; ++i_bin )
        {
          expotrunc->SetRange( timeBinHist[i_bin].binlow, timeBinHist[i_bin+nbinwidth].binhigh );
          hbin_pg->Fit( "expotrunc", "RQ" );
          lifetimedrift->SetPoint( i_bin-1, timeBinHist[i_bin+nbinwidth].bincenter, expotrunc->GetParameter( 1 ) );
          lifetimedrift->SetPointError( i_bin-1, 0, expotrunc->GetParError( 1 ) );
        }
      TCanvas * canvlifetimedrift_pg = new TCanvas( "cltd_pg", "cltd_pg", 1600, 900 );
      canvlifetimedrift_pg->cd();
      lifetimedrift->Draw( "alpe" );
      lifetimedrift->GetXaxis()->SetTitle( "Drift Distance Cut ( cm )" );
      lifetimedrift->GetYaxis()->SetTitle( "Elifetime ( us )" );
      canvlifetimedrift_pg->Update();
      canvlifetimedrift_pg->Write();
      canvlifetimedrift_pg->SaveAs( "cltd_pg.png" );
    }

  TTree * resulttree = new TTree( "results", "Results of fits" );
  Float_t landmpv;
  Float_t landwidth;
  Float_t gaussmean_l;
  Float_t gausswidth_l;
  Float_t landmpverr;
  Float_t landwidtherr;
  Float_t gaussmean_lerr;
  Float_t gausswidth_lerr;
  Float_t lxg_chi2ndf;
  Float_t lxglow;
  Float_t lxghigh;
  Bool_t fitsuccess;
  resulttree->Branch( "landmpv", &landmpv, "landmpv/F" );
  resulttree->Branch( "landwidth", &landwidth, "landwidth/F" );
  resulttree->Branch( "gaussmean_l", &gaussmean_l, "gaussmean_l/F" );
  resulttree->Branch( "gausswidth_l", &gausswidth_l, "gausswidth_l/F" );
  resulttree->Branch( "landmpverr", &landmpverr, "landmpverr/F" );
  resulttree->Branch( "landwidtherr", &landwidtherr, "landwidtherr/F" );
  resulttree->Branch( "gaussmean_lerr", &gaussmean_lerr, "gaussmean_lerr/F" );
  resulttree->Branch( "gausswidth_lerr", &gausswidth_lerr, "gausswidth_lerr/F" );
  resulttree->Branch( "lxg_chi2ndf", &lxg_chi2ndf, "lxg_chi2ndf/F" );
  resulttree->Branch( "lxglow", &lxglow, "lxglow/F" );
  resulttree->Branch( "lxghigh", &lxghigh, "lxghigh/F" );
  resulttree->Branch( "fitsuccess", &fitsuccess, "fitsuccess/O");

  Float_t gaussmean_v;
  Float_t gausswidth_v;
  Float_t vavkappa;
  Float_t vavbeta2;
  Float_t vavmean;
  Float_t vavwidth;
  Float_t vavamp;
  Float_t gaussmean_verr;
  Float_t gausswidth_verr;
  Float_t vavkappaerr;
  Float_t vavbeta2err;
  Float_t vavmeanerr;
  Float_t vavwidtherr;
  Float_t vavamperr;
  Float_t vxg_chi2ndf;
  resulttree->Branch( "gaussmean_v", &gaussmean_v, "gaussmean_v/F" );
  resulttree->Branch( "gausswidth_v", &gausswidth_v, "gausswidth_v/F" );
  resulttree->Branch( "vavkappa", &vavkappa, "vavkappa/F" );
  resulttree->Branch( "vavbeta2", &vavbeta2, "vavbeta2/F" );
  resulttree->Branch( "vavmean", &vavmean, "vavmean/F" );
  resulttree->Branch( "vavwidth", &vavwidth, "vavwidth/F" );
  resulttree->Branch( "vavamp", &vavamp, "vavamp/F" );
  resulttree->Branch( "gaussmean_verr", &gaussmean_verr, "gaussmean_verr/F" );
  resulttree->Branch( "gausswidth_verr", &gausswidth_verr, "gausswidth_verr/F" );
  resulttree->Branch( "vavkappaerr", &vavkappaerr, "vavkappaerr/F" );
  resulttree->Branch( "vavbeta2err", &vavbeta2err, "vavbeta2err/F" );
  resulttree->Branch( "vavmeanerr", &vavmeanerr, "vavmeanerr/F" );
  resulttree->Branch( "vavwidtherr", &vavwidtherr, "vavwidtherr/F" );
  resulttree->Branch( "vavamperr", &vavamperr, "vavamperr/F" );
  resulttree->Branch( "vxg_chi2ndf", &vxg_chi2ndf, "vxg_chi2ndf/F" );

  Float_t gaussmean;
  Float_t gausswidth;
  Float_t gaussmeanerr;
  Float_t gausswidtherr;
  Float_t bincenter;
  Float_t binwidtherr;
  Float_t pg_chi2ndf;
  resulttree->Branch( "binnum", &binnum, "binnum/I" );
  resulttree->Branch( "gaussmean", &gaussmean, "gaussmean/F" );
  resulttree->Branch( "gausswidth", &gausswidth, "gausswidth/F" );
  resulttree->Branch( "gaussmeanerr", &gaussmeanerr, "gaussmeanerr/F" );
  resulttree->Branch( "gausswidtherr", &gausswidtherr, "gausswidtherr/F" );
  resulttree->Branch( "bincenter", &bincenter, "bincenter/F" );
  resulttree->Branch( "binwidtherr", &binwidtherr, "binwidtherr/F" );
  resulttree->Branch( "pg_chi2ndf", &pg_chi2ndf, "pg_chi2ndf/F" );

  for ( auto const & r : results )
    {
      binnum = r.binnum;
      landmpv = r.landmpv;
      landwidth = r.landwidth;
      gaussmean = r.gaussmean;
      gausswidth = r.gausswidth;
      gaussmean_l = r.gaussmean_l;
      gausswidth_l = r.gausswidth_l;
      gaussmean_v = r.gaussmean_v;
      gausswidth_v = r.gausswidth_v;
      vavkappa = r.vavkappa;
      vavbeta2 = r.vavbeta2;
      vavmean = r.vavmean;
      vavwidth = r.vavwidth;
      vavamp = r.vavamp;
      landmpverr = r.landmpverr;
      landwidtherr = r.landwidtherr;
      gaussmeanerr = r.gaussmeanerr;
      gausswidtherr = r.gausswidtherr;
      gaussmean_lerr = r.gaussmean_lerr;
      gausswidth_lerr = r.gausswidth_lerr;
      gaussmean_verr = r.gaussmean_verr;
      gausswidth_verr = r.gausswidth_verr;
      vavkappaerr = r.vavkappaerr;
      vavbeta2err = r.vavbeta2err;
      vavmeanerr = r.vavmeanerr;
      vavwidtherr = r.vavwidtherr;
      vavamperr = r.vavamperr;
      bincenter = r.bincenter;
      binwidtherr = r.binwidtherr;
      lxg_chi2ndf = r.lxg_chi2ndf;
      vxg_chi2ndf = r.vxg_chi2ndf;
      pg_chi2ndf = r.pg_chi2ndf;
      lxglow = r.lxglow;
      lxghigh = r.lxghigh;
      fitsuccess = r.fitsuccess;
      resulttree->Fill();
    }

  std::cout << "Found Number = " << found->GetEntries() << "  Assumed Number = " << assumed->GetEntries() << std::endl;
  resulttree->Write();
  fileout->Close();
  return;
}

#endif
