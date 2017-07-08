#include "TMultiGraph.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TF1.h"
#include "TF2.h"
#include "TFitResult.h"
#include "TLine.h"
#include "TVirtualFitter.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TGraph2DErrors.h"
#include "TGraph2D.h"
#include "TMatrixD.h"
#include "TStyle.h"

const UInt_t nscale = 16;
const UInt_t npur   = 5;

void PlotMPVExpos()
{
  Double_t adctoelectrons = (1/2.808)*(1/14.0)*(1/7.615)*(1/1.602e-4);

  Double_t mcscale[nscale] = {0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};
  Int_t elife[npur] = {2500,3000,3500,4000,4500};
  Bool_t oldgain = true;
  Double_t interestingscale = 1.2;
  Int_t interestingelife = 4000;
  Double_t dqdx0 = 1.153;

  gStyle->SetLabelSize(0.05,"xyz");
  //gStyle->SetLegendTextSize(0.05);
  gStyle->SetTitleSize(0.05,"xyz");

  std::map<Int_t, std::map<UInt_t, TGraphErrors* > > mpvgraphs;
  std::map<Int_t, std::map<UInt_t, TGraphErrors* > > mpvgraphspaper;

  TMultiGraph * dqdx0mg = new TMultiGraph();

  TMultiGraph * mgmpvpaper = new TMultiGraph();
  TLegend * legpaper = new TLegend(0.249,0.75,0.998,0.99);

  TMultiGraph * mg = new TMultiGraph(); //elifeVnoise
  TMultiGraph * mgsyst = new TMultiGraph(); //elifeVnoise systematics

  Int_t palette[nscale];
  Double_t Red[] =    {0.831, 0.0,   0.0,   0.133, 1.0};
  Double_t Green[] =  {0.349, 0.267, 1.0,   1.0,   0.933};
  Double_t Blue[]   = {0.329, 1.0,   0.800, 0.0,   0.0};
  Double_t Length[] = {0., .25, .50, .75, 1.0};
  Int_t FI = TColor::CreateGradientColorTable( 5, Length, Red, Green, Blue, nscale );
  for (unsigned int i=0; i<nscale; i++ ) palette[i] = FI+i;

  Int_t color2500=2,color3000=3,color3500=4,color4000=6,color4500=12;

  Float_t fitstart = 100;
  Float_t fitend = 1000;

  TGraphErrors * grdata = new TGraphErrors();

  TGraphErrors * grelife = new TGraphErrors();
  Int_t numelife = 0;

  Double_t elife_data = 0;
  Double_t elife_data_err = 0;

  Double_t dqdx0data=0, dqdx0errdata=0;

  std::map<Int_t,std::vector<Double_t> > dqdx0values;

  TGraph2D * gr2d = new TGraph2D();
  Int_t numgr2d = 0;

  //TGraph2D * h2 = new TGraph2D();

  TMatrixD X(nscale*npur,6);
  TMatrixD Xerr(nscale*npur,6);
  TMatrixD Y(nscale*npur,1);

  TFile * datafile = TFile::Open( "/home/mthiesse/PurityAnalysis/PurityAnalysis_data.root", "READ" );
  if ( datafile && datafile->IsOpen() )
    {
      TTreeReader datareader( "results", datafile );
      TTreeReaderValue<Float_t> landmpvdata( datareader, "landmpv" );
      TTreeReaderValue<Float_t> landwidthdata( datareader, "landwidth" );
      TTreeReaderValue<Float_t> gaussmean_ldata( datareader, "gaussmean_l" );
      TTreeReaderValue<Float_t> gausswidth_ldata( datareader, "gausswidth_l" );
      TTreeReaderValue<Float_t> landmpverrdata( datareader, "landmpverr" );
      TTreeReaderValue<Float_t> landwidtherrdata( datareader, "landwidtherr" );
      TTreeReaderValue<Float_t> gaussmean_lerrdata( datareader, "gaussmean_lerr" );
      TTreeReaderValue<Float_t> gausswidth_lerrdata( datareader, "gausswidth_lerr" );
      TTreeReaderValue<Float_t> lxg_chi2ndfdata( datareader, "lxg_chi2ndf" );
      TTreeReaderValue<Float_t> lxglowdata( datareader, "lxglow" );
      TTreeReaderValue<Float_t> lxghighdata( datareader, "lxghigh" );
      TTreeReaderValue<Float_t> bincenterdata( datareader, "bincenter" );
      TTreeReaderValue<Bool_t> fitsuccessdata( datareader, "fitsuccess" );

      Int_t num = 0;
      while (datareader.Next())
        {
          if (*fitsuccessdata && *lxg_chi2ndfdata<5 && *bincenterdata<2012 && *bincenterdata>0)
            {
              grdata->SetPoint(num,*bincenterdata,*landmpvdata);
              //grdata->SetPointError(num,0,*landmpverrdata);
              grdata->SetPointError(num,91.45/sqrt(12),*landmpverrdata);
              ++num;
            }
        }

      if (num != 0)
        {
          TF1 * expo = new TF1( "expo", "[0]*exp( -x/[1] )", fitstart, fitend );
          expo->SetParNames( "dQdx0", "eLifetime" );
          expo->SetParameters( 3000, 3000 );
          expo->SetParLimits(0,0,10000);
          expo->SetParLimits(1,1000,100000);
          expo->SetLineWidth(4);

          TFitResultPtr r = grdata->Fit("expo","SRBQ");
          std::cout << "Data: " << r->Parameter(1) << " " << r->ParError(1) << std::endl;

          elife_data = r->Parameter(1);
          elife_data_err = r->ParError(1);

          grdata->SetLineWidth(4);
          grdata->SetMarkerStyle(20);
          grdata->SetMarkerSize(2);
          grdata->SetMarkerColor(1);
          mgmpvpaper->Add(grdata,"p");
          legpaper->AddEntry(grdata,"35-ton Data","pel");
          legpaper->AddEntry(expo,"Exponential Fit to Data","l");

          dqdx0data = r->Parameter(0);
          dqdx0errdata = r->ParError(0);

          TGraphErrors * data = new TGraphErrors();
          data->SetPoint(0,1, expo->GetParameter(1));
          data->SetPointError(0,0, expo->GetParError(1));
          data->SetMarkerStyle(29);
          data->SetMarkerSize(5);
          data->SetLineWidth(4);
          data->SetLineColor(1);
          data->SetMarkerColor(1);
          data->SetTitle("35-ton Data");
          mg->Add(data);

        }
      delete datafile;
    }

  Int_t irow = -1;
  for ( UInt_t ipur = 0; ipur < npur; ++ipur )
    {
      TGraphErrors * elifeVnoise = new TGraphErrors();
      Int_t numelifeVnoise = 0;
      elifeVnoise->SetTitle(TString::Format("Sim %u#mus",elife[ipur]));

      TGraphErrors * elifeVnoisesyst = new TGraphErrors();
      Int_t numelifeVnoisesyst = 0;
      elifeVnoisesyst->SetTitle(TString::Format("Sim %u#mus",elife[ipur]));

      TGraphErrors * dqdx0gr = new TGraphErrors();
      Int_t numdqdx0 = 0;
      dqdx0gr->SetTitle(TString::Format("Sim %u#mus",elife[ipur]));

      TCanvas * canvexpo = new TCanvas("canvexpo","",2000,1600);
      TMultiGraph * mgmpv = new TMultiGraph();
      //TH1F * hpad = new TH1F("hpad","",100,0,1050);
      //hpad->SetMinimum(1000);
      //hpad->SetMaximum(5000);
      //hpad->SetStats(false);
      mgmpv->SetTitle(";Drift time [#mus];Most Probable dQ/dx [ADC/cm]");
      //hpad->Draw();
      TLegend * leg = new TLegend(0.78,0.68,0.98,0.98);
      leg->SetHeader(TString::Format("Simulated eLifetime = %u#mus",elife[ipur]),"C");
      //Double_t maxval = -99999;
      //Double_t minval = 99999;

      for ( UInt_t iscale = 0; iscale < nscale; ++iscale )
        {
          TString filename = TString::Format( "PurityAnalysis_%uus_mcscale%.1f.root", elife[ipur], mcscale[iscale] );

          TFile * file = TFile::Open( filename, "READ" );
          if ( !file || !file->IsOpen() ) continue;

          TTreeReader reader( "results", file );
          TTreeReaderValue<Float_t> landmpv( reader, "landmpv" );
          TTreeReaderValue<Float_t> landwidth( reader, "landwidth" );
          TTreeReaderValue<Float_t> gaussmean_l( reader, "gaussmean_l" );
          TTreeReaderValue<Float_t> gausswidth_l( reader, "gausswidth_l" );
          TTreeReaderValue<Float_t> landmpverr( reader, "landmpverr" );
          TTreeReaderValue<Float_t> landwidtherr( reader, "landwidtherr" );
          TTreeReaderValue<Float_t> gaussmean_lerr( reader, "gaussmean_lerr" );
          TTreeReaderValue<Float_t> gausswidth_lerr( reader, "gausswidth_lerr" );
          TTreeReaderValue<Float_t> lxg_chi2ndf( reader, "lxg_chi2ndf" );
          TTreeReaderValue<Float_t> lxglow( reader, "lxglow" );
          TTreeReaderValue<Float_t> lxghigh( reader, "lxghigh" );
          TTreeReaderValue<Float_t> bincenter( reader, "bincenter" );
          TTreeReaderValue<Bool_t> fitsuccess( reader, "fitsuccess" );

          mpvgraphs[elife[ipur]][iscale] = new TGraphErrors();
          mpvgraphspaper[elife[ipur]][iscale] = new TGraphErrors();
          Int_t num = 0;
          while (reader.Next())
            {
              if (*fitsuccess && *lxg_chi2ndf<10 && *bincenter<2012 && *bincenter>0)
                {
                  //if (*landmpv * adctoelectrons / mcscale[iscale] > maxval) maxval = *landmpv * adctoelectrons / mcscale[iscale];
                  //if (*landmpv * adctoelectrons / mcscale[iscale] < minval) minval = *landmpv * adctoelectrons / mcscale[iscale];
                  mpvgraphs[elife[ipur]][iscale]->SetPoint(num,(*bincenter),(*landmpv) /* (oldgain ? originaloffset : residualoffset)*/);
                  mpvgraphspaper[elife[ipur]][iscale]->SetPoint(num,(*bincenter),(*landmpv) /* (oldgain ? originaloffset : residualoffset)*/);
                  mpvgraphspaper[elife[ipur]][iscale]->SetPointError(num,0 /*91.45/sqrt(12)*/,*landmpverr + *landmpv*fabs(dqdx0-mcscale[iscale]));
                  mpvgraphs[elife[ipur]][iscale]->SetPointError(num,0 /*91.45/sqrt(12)*/,*landmpverr);
                  ++num;
                }
            }

          if (num == 0) continue;

          TF1 * expo = new TF1( "expo", "[0]*exp( -x/[1] )", fitstart, fitend );
          expo->SetParNames( "dQdx0", "eLifetime" );
          expo->SetParameters( 3000, 3000 );
          expo->SetParLimits(0,1000,10000);
          expo->SetParLimits(1,1000,10000);

          TFitResultPtr r = mpvgraphs[elife[ipur]][iscale]->Fit("expo","SRNQ");
          std::cout << "Sim (elife=" << elife[ipur] << ", scale=" << mcscale[iscale] << "): " << r->Parameter(1) << " " << r->ParError(1) << ", " << 100*r->ParError(1)/(r->Parameter(1)-elife[ipur]) << "%" << std::endl;

          //hpad->SetMinimum(minval*0.98);
          //hpad->SetMaximum(maxval*1.1);

          mpvgraphs[elife[ipur]][iscale]->SetLineWidth(2);
          mpvgraphs[elife[ipur]][iscale]->SetLineColor(palette[iscale]);
          mpvgraphs[elife[ipur]][iscale]->SetMarkerStyle(20);
          mpvgraphs[elife[ipur]][iscale]->SetMarkerSize(3);
          mpvgraphs[elife[ipur]][iscale]->SetMarkerColor(palette[iscale]);
          mgmpv->Add(mpvgraphs[elife[ipur]][iscale],"lp");
          leg->AddEntry(mpvgraphs[elife[ipur]][iscale],TString::Format("MCscale=%.1f",mcscale[iscale]),"lp");

          if (elife[ipur] == interestingelife && mcscale[iscale] == interestingscale)
            {
/*
              TGraphErrors * hint_lxg = new TGraphErrors(mpvgraphs[elife[ipur]][iscale]->GetN());
              for (int i = 0; i < mpvgraphs[elife[ipur]][iscale]->GetN(); ++i)
                {
                  hint_lxg->SetPoint(i,mpvgraphs[elife[ipur]][iscale]->GetX()[i],0);
                }
              ( TVirtualFitter::GetFitter() )->GetConfidenceIntervals( hint_lxg, 0.95 );
              hint_lxg->SetLineColor( kRed );
              hint_lxg->SetLineWidth(4);
              hint_lxg->SetMarkerStyle(20);
              hint_lxg->SetMarkerSize(3);
              hint_lxg->SetMarkerColor( kRed);
              hint_lxg->SetFillColor(kRed);
              hint_lxg->SetFillColorAlpha(kRed,0.35);
              mgmpvpaper->Add(hint_lxg, "e3" );
 */

              TGraphErrors * comparegr = (TGraphErrors*)mpvgraphspaper[elife[ipur]][iscale]->Clone(TString::Format("%u#mus Lifetime, MCscale=%.1f",interestingelife,interestingscale));
              comparegr->SetLineColor(kBlue);
              comparegr->SetLineWidth(4);
              comparegr->SetFillColor(kGreen);
              comparegr->SetFillColorAlpha(kGreen,0.35);
              mgmpvpaper->Add(comparegr,"le3");
              legpaper->AddEntry(comparegr,TString::Format("Simulation (eLife=%u#mus, MCscale=%.1f)",interestingelife,interestingscale),"l");
              legpaper->AddEntry(comparegr,"Uncertainty due to Mis-modelling","f");
            }

          if (expo->GetParameter(1) < 10000 && expo->GetParameter(1) > 1000 && expo->GetParError(1) < 500)
            {
              elifeVnoise->SetPoint(numelifeVnoise,(mcscale[iscale])/dqdx0 /*(oldgain ? originaloffset : residualoffset)*/,expo->GetParameter(1));
              elifeVnoise->SetPointError(numelifeVnoise,/*(oldgain ? originaloffseterr : residualoffseterr)*/ 0,expo->GetParError(1));
              ++numelifeVnoise;

              //h2->SetBinContent(h2->FindBin(mcscale[iscale]/1.15,elife[ipur]),h2->GetBinContent(h2->FindBin(mcscale[iscale]/1.15,elife[ipur]))+expo->GetParameter(1));

              Double_t syst = expo->GetParError(1)*100/fabs(expo->GetParameter(1)-elife[ipur]);
              if (syst < 1000 && syst > -1000) { elifeVnoisesyst->SetPoint(numelifeVnoisesyst,(mcscale[iscale])/dqdx0,syst); ++numelifeVnoisesyst; }

              dqdx0gr->SetPoint(numdqdx0,expo->GetParameter(0),mcscale[iscale]/dqdx0);
              dqdx0gr->SetPointError(numdqdx0,expo->GetParError(0),0);
              ++numdqdx0;

              dqdx0values[iscale].push_back(expo->GetParameter(0));

              gr2d->SetPoint(numgr2d,mcscale[iscale]/dqdx0,expo->GetParameter(1),elife[ipur]);
              //gr2d->SetPointError(numgr2d,0.05,expo->GetParError(1),0);
              ++numgr2d;

              Float_t corr = dqdx0;
              Float_t correrr = 0.05;
              X[irow][0] = 1; X[irow][1] = mcscale[iscale]/corr; X[irow][2] = expo->GetParameter(1);
              X[irow][3] = (mcscale[iscale]/corr)*expo->GetParameter(1); X[irow][4] = (mcscale[iscale]/corr)*(mcscale[iscale]/corr);
              X[irow][5] = expo->GetParameter(1)*expo->GetParameter(1);
              Y[irow][0] = elife[ipur];
              Xerr[irow][0] = 0; Xerr[irow][1] = correrr; Xerr[irow][2] = expo->GetParError(1);
              Xerr[irow][3] = (mcscale[iscale]/corr)*expo->GetParameter(1)*sqrt(pow(correrr/(mcscale[iscale]/corr),2)+pow(expo->GetParError(1)/expo->GetParameter(1),2));
              Xerr[irow][4] = (mcscale[iscale]/corr)*2*correrr;
              Xerr[irow][5] = expo->GetParameter(1)*2*expo->GetParError(1);
              ++irow;
            }
          if (mcscale[iscale] == interestingscale)
            {
              grelife->SetPoint(numelife,elife[ipur],expo->GetParameter(1));
              grelife->SetPointError(numelife,0,expo->GetParError(1));
              ++numelife;
            }

          delete file;
        }
      Int_t thiscolor = 1;
      if (elife[ipur]==2500) thiscolor = color2500;
      else if (elife[ipur]==3000) thiscolor = color3000;
      else if (elife[ipur]==3500) thiscolor = color3500;
      else if (elife[ipur]==4000) thiscolor = color4000;
      else if (elife[ipur]==4500) thiscolor = color4500;

      mgmpv->Add(grdata,"p0");
      leg->AddEntry(grdata,"35-ton Data","p");
      mgmpv->Draw("a");
      mgmpv->GetYaxis()->SetTitleOffset(1.5);
      leg->Draw();
      canvexpo->Update();
      canvexpo->SaveAs(TString::Format("PlotMPVExpos_%u.png",elife[ipur]));
      elifeVnoise->SetMarkerStyle(22);
      elifeVnoise->SetMarkerSize(3);
      elifeVnoise->SetMarkerColor(thiscolor);
      elifeVnoise->SetLineColor(thiscolor);
      elifeVnoise->SetLineWidth(4);
      if (elifeVnoise->GetN() != 0) mg->Add(elifeVnoise);
      elifeVnoisesyst->SetMarkerStyle(20);
      elifeVnoisesyst->SetMarkerSize(3);
      elifeVnoisesyst->SetMarkerColor(thiscolor);
      elifeVnoisesyst->SetLineColor(thiscolor);
      elifeVnoisesyst->SetLineWidth(4);
      if (elifeVnoisesyst->GetN() != 0) mgsyst->Add(elifeVnoisesyst);

      dqdx0gr->SetMarkerStyle(20);
      dqdx0gr->SetMarkerSize(2);
      dqdx0gr->SetMarkerColor(thiscolor);
      //if (dqdx0gr->GetN() != 0) dqdx0mg->Add(dqdx0gr);
    }

  ////////////////////////////////////////////

  X.ResizeTo(irow,6);
  Y.ResizeTo(irow,1);
  Xerr.ResizeTo(irow,6);
  TMatrixD XTX(X,TMatrixD::kTransposeMult,X);
  TMatrixD XTXI = XTX.Invert();
  TMatrixD XTY(X,TMatrixD::kTransposeMult,Y);
  TMatrixD b = XTXI*XTY;
  TMatrixD Err = Y+Xerr*b-X*b;
  TMatrixD SSE(Err,TMatrixD::kTransposeMult,Err);
  Double_t MSE = SSE[0][0]/(irow-6);
  TMatrixD COV = XTXI*MSE;

  TMatrixD Xh(6,1);
  Xh[0][0] = 1; Xh[1][0] = 1; Xh[2][0] = elife_data; Xh[3][0] = elife_data;
  Xh[4][0] = 1; Xh[5][0] = elife_data*elife_data;
  TMatrixD XhErr(6,1);
  XhErr[0][0] = 0; XhErr[1][0] = 0; XhErr[2][0] = elife_data_err; XhErr[3][0] = elife_data_err;
  XhErr[4][0] = 0; XhErr[5][0] = elife_data*2*elife_data_err;
  TMatrixD predictedlife(Xh,TMatrixD::kTransposeMult,b);

  TMatrixD errlife(XhErr,TMatrixD::kTransposeMult,b);
  TMatrixD intermediateErr(XhErr,TMatrixD::kTransposeMult,XTXI*XhErr);
  Double_t seerr = sqrt(intermediateErr[0][0]*MSE);
  Double_t prederr = sqrt(MSE+seerr*seerr);
  TMatrixD intermediate(Xh,TMatrixD::kTransposeMult,XTXI*Xh);
  Double_t se = sqrt(intermediate[0][0]*MSE);
  Double_t pred = sqrt(MSE+se*se);

  std::cout << "predictedlife: " << predictedlife[0][0] << std::endl;
  std::cout << "combined errors: " << sqrt(pred*pred+prederr*prederr) << std::endl;

  std::cout << "Fitting 2d function" << std::endl;
  TCanvas * canvgr2d = new TCanvas("canvgr2d","",2000,1600);
  TF2 * poly2 = new TF2("poly2","[0]+[1]*x+[2]*y+[3]*x*y+[4]*x*x+[5]*y*y",0,2,2000,5000);
  poly2->SetParameter(0,10);
  poly2->SetParameter(1,0.1);
  poly2->SetParameter(2,1);
  std::cout << "gr2d points " << gr2d->GetN() << std::endl;
  TH2F *grint2 = new TH2F("grint2","",50,0,2,50,2000,5000);
  //for (int i = 0; i < gr2d->GetN(); ++i)
  //  {
  //    grint2->SetBinContent(i,gr2d->GetX()[i],gr2d->GetY()[i],gr2d->GetZ()[i]);
  //  }
  gr2d->Fit("poly2","BMER");
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint2,0.68);

  poly2->SetNpx(30);
  poly2->SetNpy(30);
  poly2->SetFillColor(kBlue);
  poly2->Draw("surf4");

  //grint2->SetNpx(30);
  //grint2->SetNpy(30);
  grint2->SetMarkerStyle(24);
  grint2->SetMarkerSize(0.7);
  grint2->SetMarkerColor(kRed);
  grint2->SetLineColor(kRed);
  grint2->SetFillColor(kRed);
  grint2->SetFillColorAlpha(kRed,0.35);

  gr2d->Draw("P0 same");
  //grint2->Draw("e3 same");
  canvgr2d->Update();
  std::cout << "Done fiting 2d function" << std::endl;
  std::cout << "Value of function at data: scale=1, observed=" << elife_data << "  = " << poly2->Eval(1,elife_data) << " +/- " << grint2->GetBinError(grint2->FindBin(1,elife_data)) << " in bin " << grint2->FindBin(1,elife_data) << " content " << grint2->GetBinContent(grint2->FindBin(1,elife_data)) << std::endl;

  //TCanvas * canvh2 = new TCanvas("canvh2","",2000,1600);
  //h2->Draw("colz");
  //canvh2->Update();

  TCanvas * canvelifeVnoise = new TCanvas("canvelifeVnoise","",2000,1600);
  canvelifeVnoise->SetLeftMargin(0.14);
  mg->Draw("ap");
  mg->GetXaxis()->SetTitle("MCscale");
  mg->GetYaxis()->SetTitle("Observed #tau [ #mus]");
  mg->GetYaxis()->SetTitleOffset(1.5);
  TLegend * thisleg = canvelifeVnoise->BuildLegend(0.421,0.804,0.931,0.99,"","lpe");
  thisleg->SetNColumns(2);
  Float_t linelowx = (0.5/1.149)-0.05;
  Float_t linehix = (2.0/1.149)+0.05;
  TLine * line3000 = new TLine(linelowx,3000,linehix,3000);
  line3000->SetLineStyle(kDashed);
  line3000->SetLineColor(color3000);
  line3000->SetLineWidth(5);
  //line3000->Draw();
  TLine * line3500 = new TLine(linelowx,3500,linehix,3500);
  line3500->SetLineStyle(kDashed);
  line3500->SetLineColor(color3500);
  line3500->SetLineWidth(5);
  //line3500->Draw();
  TLine * line2500 = new TLine(linelowx,2500,linehix,2500);
  line2500->SetLineStyle(kDashed);
  line2500->SetLineColor(color2500);
  line2500->SetLineWidth(5);
  //line2500->Draw();
  TLine * line4000 = new TLine(linelowx,4000,linehix,4000);
  line4000->SetLineStyle(kDashed);
  line4000->SetLineColor(color4000);
  line4000->SetLineWidth(5);
  //line4000->Draw();
  TLine * line4500 = new TLine(linelowx,4500,linehix,4500);
  line4500->SetLineStyle(kDashed);
  line4500->SetLineColor(color4500);
  line4500->SetLineWidth(5);
  //line4500->Draw();
  canvelifeVnoise->Update();
  canvelifeVnoise->SaveAs("canvelifeVnoise.png");

  TCanvas * canvdqdx0 = new TCanvas("canvdqdx0","",2000,1600);
  TGraphErrors * dqdx0average = new TGraphErrors();
  dqdx0average->SetMarkerStyle(20);
  dqdx0average->SetMarkerSize(2);
  dqdx0average->SetMarkerColor(4);
  dqdx0average->SetLineWidth(4);
  dqdx0average->SetLineColor(2);
  Int_t idqdx = 0;
  for (auto dqdx0scale : dqdx0values)
    {
      for (auto dqdx0 : dqdx0scale.second)
        {
          dqdx0average->SetPoint(idqdx,dqdx0,mcscale[dqdx0scale.first]);
          ++idqdx;
        }
      std::cout << "MCscale=" << mcscale[dqdx0scale.first] << "  dQ/dx_0=" << TMath::Mean(dqdx0scale.second.begin(),dqdx0scale.second.end()) << std::endl;
      //dqdx0average->SetPointError(idqdx,TMath::StdDev(dqdx0scale.second.begin(),dqdx0scale.second.end()),0);
      //++idqdx;
    }
  TF1 * pol2 = new TF1("pol2","pol2",1000,5000);
  pol2->SetParameters(-1.29e-1,5.65e-4,-2e-8);
  pol2->SetParLimits(0,-1,1);
  pol2->SetParLimits(1,-1,1);
  pol2->SetParLimits(2,-1e-3,-1e-12);
  pol2->SetLineWidth(4);
  dqdx0average->Fit("pol2");
  dqdx0mg->Add(dqdx0average);
  dqdx0mg->Draw("ap");
  dqdx0mg->SetTitle(";dQ/dx0 [ADC/cm];Simulated MCscale");
  //canvdqdx0->BuildLegend(0.1,0.7,0.3,0.9,"","lpe");
  Float_t val = pol2->Eval(dqdx0data);
  Float_t unc = fabs(val-pol2->Eval(dqdx0data-dqdx0errdata));
  TPaveText * valtext = new TPaveText(0.3,0.15,0.89,0.25,"brNDC");
  valtext->AddText(TString::Format("MCscale at data dQ/dx0 = %.3f #pm %.3f",val,unc));
  //valtext->Draw();
  canvdqdx0->Update();
  TLine * dqdx0dataline = new TLine(dqdx0data,canvdqdx0->GetUymin(),dqdx0data,canvdqdx0->GetUymax());
  dqdx0dataline->SetLineStyle(1);
  dqdx0dataline->SetLineWidth(4);
  dqdx0dataline->SetLineColor(4);
  dqdx0dataline->Draw();
  valtext->Draw();
  TLine * dqdx0simline = new TLine(canvdqdx0->GetUxmin(),val,canvdqdx0->GetUxmax(),val);
  dqdx0simline->SetLineStyle(9);
  dqdx0simline->SetLineWidth(4);
  dqdx0simline->SetLineColor(6);
  dqdx0simline->Draw();
  TLegend * dqdx0leg = new TLegend(0.15,0.8,0.673,0.98);
  dqdx0leg->AddEntry(dqdx0average,"Simulation","pl");
  //dqdx0leg->Add(pol2,"Quadratic Interpolation","l");
  dqdx0leg->AddEntry(dqdx0dataline,"Data dQ/dx0","l");
  dqdx0leg->AddEntry(dqdx0simline,"Interpolated Best dQ/dx0","l");
  dqdx0leg->Draw();
  canvdqdx0->Update();
  canvdqdx0->SaveAs("canvdqdx0.png");

  TCanvas * canvelifeVnoisesyst = new TCanvas("canvelifeVnoisesyst","",2000,1600);
  canvelifeVnoisesyst->SetLeftMargin(0.11);
  mgsyst->Draw("alp");
  mgsyst->GetXaxis()->SetTitle("MCscale");
  mgsyst->GetYaxis()->SetTitle("eLifetime Error [%]");
  mgsyst->GetYaxis()->SetTitleOffset(1.5);
  canvelifeVnoisesyst->BuildLegend(0.7,0.7,0.9,0.9,"","lpe");
  canvelifeVnoisesyst->Update();
  canvelifeVnoisesyst->SaveAs("canvelifeVnoisesyst.png");

  TCanvas * canvmpvpaper = new TCanvas("canvmpvpaper","",2000,1600);
  canvmpvpaper->cd();
  canvmpvpaper->SetLeftMargin(0.14);
  mgmpvpaper->Draw("a");
  mgmpvpaper->GetYaxis()->SetTitle("Most Probable dQ/dx [ADC/cm]");
  mgmpvpaper->GetYaxis()->SetTitleOffset(1.4);
  mgmpvpaper->GetXaxis()->SetTitle("Drift Time [ #mus]");
  legpaper->Draw();
  canvmpvpaper->Update();
  canvmpvpaper->SaveAs("canvmpvpaper.png");

  TCanvas * canvelife = new TCanvas("canvelife","",2000,1600);
  canvelife->cd();
  grelife->SetMarkerStyle(20);
  grelife->SetMarkerSize(2);
  grelife->Draw("ape");
  grelife->SetTitle(";True lifetime [ #mus]; Observed lifetime [ #mus]");
  grelife->GetYaxis()->SetTitleOffset(1.4);
  TF1 * datalife = new TF1("datalife","pol0",2400, 4600);
  datalife->SetParameter(0,elife_data);
  datalife->SetParError(0,elife_data_err);
  datalife->SetLineWidth(4);
  datalife->SetLineColor(kBlue);
  datalife->SetFillColor(kBlue);
  datalife->SetFillColorAlpha(kBlue,0.35);
  datalife->Draw("le3 same");
  TBox * datalifeerr = new TBox(2400,elife_data-elife_data_err,4600,elife_data+elife_data_err);
  datalifeerr->SetFillColor(kBlue);
  datalifeerr->SetFillColorAlpha(kBlue,0.35);
  datalifeerr->Draw();
  if (grelife->GetN() != 0) grelife->Fit("pol2","Q");
  TF1 * f = grelife->GetFunction("pol2");
  TH1F * grelifeband = new TH1F("grelifeband","",400,2400,4600);
  ( TVirtualFitter::GetFitter() )->GetConfidenceIntervals( grelifeband, 0.68 );
  grelifeband->SetFillColor(kRed);
  grelifeband->SetFillColorAlpha(kRed,0.35);
  grelifeband->SetLineColor(kRed);
  grelifeband->SetLineWidth(2);
  grelifeband->Draw("le3 same");

  Float_t minx=-1,miny=99999;
  Float_t minx_bb=-1,miny_bb=99999;
  Float_t minx_bt=-1,miny_bt=99999;
  Float_t minx_tb=-1,miny_tb=99999;
  Float_t minx_tt=-1,miny_tt=99999;
  for (Int_t i = 0; i < grelifeband->GetNbinsX(); ++i)
    {
      Float_t val = TMath::Abs(grelifeband->GetBinContent(i) - datalife->Eval(grelifeband->GetBinCenter(i)));
      if (val < miny) { miny = val; minx = grelifeband->GetBinCenter(i); }
      val = TMath::Abs(grelifeband->GetBinContent(i)-grelifeband->GetBinError(i) - (datalife->Eval(grelifeband->GetBinCenter(i))-elife_data_err));
      if (val < miny_bb) { miny_bb = val; minx_bb = grelifeband->GetBinCenter(i); }
      val = TMath::Abs(grelifeband->GetBinContent(i)-grelifeband->GetBinError(i) - (datalife->Eval(grelifeband->GetBinCenter(i))+elife_data_err));
      if (val < miny_bt) { miny_bt = val; minx_bt = grelifeband->GetBinCenter(i); }
      val = TMath::Abs(grelifeband->GetBinContent(i)+grelifeband->GetBinError(i) - (datalife->Eval(grelifeband->GetBinCenter(i))-elife_data_err));
      if (val < miny_tb) { miny_tb = val; minx_tb = grelifeband->GetBinCenter(i); }
      val = TMath::Abs(grelifeband->GetBinContent(i)+grelifeband->GetBinError(i) - (datalife->Eval(grelifeband->GetBinCenter(i))+elife_data_err));
      if (val < miny_tt) { miny_tt = val; minx_tt = grelifeband->GetBinCenter(i); }
    }
  std::cout << "minx=" << minx << "  miny=" << miny << std::endl;
  std::cout << "minx_bb=" << minx_bb << "  miny_bb=" << miny_bb << std::endl;
  std::cout << "minx_bt=" << minx_bt << "  miny_bt=" << miny_bt << std::endl;
  std::cout << "minx_tb=" << minx_tb << "  miny_tb=" << miny_tb << std::endl;
  std::cout << "minx_tt=" << minx_tt << "  miny_tt=" << miny_tt << std::endl;
  Float_t true_lifetime_data = minx;
  Float_t true_lifetime_data_errup = TMath::Max(minx_bb,TMath::Max(minx_bt,TMath::Max(minx_tb,minx_tt)))-minx;
  Float_t true_lifetime_data_errlo = minx-TMath::Min(minx_bb,TMath::Min(minx_bt,TMath::Min(minx_tb,minx_tt)));
  std::cout << "Lifetime: " << minx << " +" << TMath::Max(minx_bb,TMath::Max(minx_bt,TMath::Max(minx_tb,minx_tt)))-minx << " -" << minx-TMath::Min(minx_bb,TMath::Min(minx_bt,TMath::Min(minx_tb,minx_tt))) << std::endl;

  TGraphErrors * placeholder = new TGraphErrors();
  placeholder->SetLineWidth(4);
  placeholder->SetLineColor(kBlue);
  placeholder->SetFillColor(kBlue);
  placeholder->SetFillColorAlpha(kBlue,0.35);
  placeholder->Draw("e3l");
  TLegend * legelife = new TLegend(0.15,0.8,0.7,0.95);
  legelife->AddEntry(grelife,"Simulated data","ep");
  legelife->AddEntry(placeholder,"35-ton data (with error)","fl");
  legelife->AddEntry(grelifeband,"2D polynomial fit (inc. 68\% confidence interval)","fl");
  legelife->Draw();
  TPaveText * tp = new TPaveText(0.5,0.2,0.85,0.3,"brNDC");
  tp->AddText(TString::Format("True Data Lifetime: %.f ^{+%.f}_{-%.f} #mus",true_lifetime_data,true_lifetime_data_errup,true_lifetime_data_errlo));
  tp->Draw();
  canvelife->Update();
  canvelife->SaveAs("canvelife.png");

}
