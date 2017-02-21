#include "TMultiGraph.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TLine.h"

const UInt_t nscale = 22;
const UInt_t npur   = 1;

void PlotMPVExpos()
{
  Double_t adctoelectrons = (1/2.808)*(1/14.0)*(1/7.615)*(1/1.602e-4);

  //Double_t mcscale[nscale] = {1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0};
  Double_t mcscale[nscale] = {0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
  Int_t elife[npur] = {2};
  //Int_t elife[npur] = {1, 2, 3, 5, 8};

  std::map<Int_t, std::map<UInt_t, TGraph* > > mpvgraphs;

  TMultiGraph * mg = new TMultiGraph();

  Int_t palette[nscale];
  Double_t Red[] = {0., 0.0, 1.0, 1.0, 1.0};
  Double_t Green[] = {0., 0.0, 0.0, 1.0, 1.0};
  Double_t Blue[]   = {0., 1.0, 0.0, 0.0, 1.0};
  Double_t Length[] = {0., .25, .50, .75, 1.0};
  Int_t FI = TColor::CreateGradientColorTable( 5, Length, Red, Green, Blue, nscale );
  for ( int i=0; i<nscale; i++ ) palette[i] = FI+i;

  for ( UInt_t ipur = 0; ipur < npur; ++ipur )
    {
      TGraphErrors * elifeVnoise = new TGraphErrors();
      Int_t numelifeVnoise = 0;
      elifeVnoise->SetTitle(TString::Format("Simulated %ums",elife[ipur]));

      TCanvas * canvexpo = new TCanvas("canvexpo","",2000,1600);
      TH1F * hpad = new TH1F("hpad","",100,0,1050);
      hpad->SetMinimum(0);
      hpad->SetMaximum(5000);
      hpad->SetStats(false);
      hpad->GetXaxis()->SetTitle("Drift time (us)");
      hpad->GetYaxis()->SetTitle("Most Probable dQ/dx (electrons/cm)");
      hpad->Draw();
      TLegend * leg = new TLegend(0.7,0.6,0.9,0.9);
      leg->SetHeader(TString::Format("Simulated eLifetime = %ums",elife[ipur]),"C");
      Double_t maxval = -99999;
      Double_t minval = 99999;

      for ( UInt_t iscale = 0; iscale < nscale; ++iscale )
        {
          TString filename = TString::Format( "PurityAnalysis_newgain_%ums_mcscale%.1f.root", elife[ipur], mcscale[iscale] );

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

          mpvgraphs[elife[ipur]][iscale] = new TGraph();
          Int_t num = 0;
          while (reader.Next())
            {
              if (*fitsuccess && *lxg_chi2ndf<5)
                {
                  if (*landmpv * adctoelectrons > maxval) maxval = *landmpv * adctoelectrons;
                  if (*landmpv * adctoelectrons < minval) minval = *landmpv * adctoelectrons;
                  mpvgraphs[elife[ipur]][iscale]->SetPoint(num,*bincenter,*landmpv * adctoelectrons);
                  //mpvgraphs[elife[ipur]][iscale]->SetPoint(num,*bincenter,*landmpverr * adctoelectrons - 1000*elife[ipur]);
                  ++num;
                }
            }

          Float_t fitstart = 100;
          Float_t fitend = 1000;
          TF1 * expo = new TF1( "expo", "[0]*exp( -x/[1] )", fitstart, fitend );
          expo->SetParNames( "dQdx0", "eLifetime" );
          expo->SetParameters( 3000, 3000 );

          TFitResultPtr r = mpvgraphs[elife[ipur]][iscale]->Fit("expo","SR");
          std::cout << r->Parameter(1) << " " << r->ParError(1) << std::endl;

          hpad->SetMinimum(minval*0.95);
          hpad->SetMaximum(maxval*1.1);

          mpvgraphs[elife[ipur]][iscale]->SetLineWidth(2);
          mpvgraphs[elife[ipur]][iscale]->SetLineColor(palette[iscale]);
          mpvgraphs[elife[ipur]][iscale]->SetMarkerStyle(20);
          mpvgraphs[elife[ipur]][iscale]->SetMarkerSize(3);
          mpvgraphs[elife[ipur]][iscale]->SetMarkerColor(palette[iscale]);
          mpvgraphs[elife[ipur]][iscale]->Draw("lp same");
          leg->AddEntry(mpvgraphs[elife[ipur]][iscale],TString::Format("MCscale=%.1f",mcscale[iscale]),"lp");

          elifeVnoise->SetPoint(numelifeVnoise,mcscale[iscale],expo->GetParameter(1));
          elifeVnoise->SetPointError(numelifeVnoise,0,expo->GetParError(1));
          ++numelifeVnoise;

          delete file;
        }
      leg->Draw();
      canvexpo->Update();
      canvexpo->SaveAs("PlotMPVExpos.png");
      elifeVnoise->SetMarkerStyle(20);
      elifeVnoise->SetMarkerSize(3);
      elifeVnoise->SetMarkerColor(ipur+2);
      elifeVnoise->SetLineColor(ipur+2);
      elifeVnoise->SetLineWidth(3);
      mg->Add(elifeVnoise);
    }

  TGraphErrors * data = new TGraphErrors();
  data->SetPoint(0,1.0,5544.43);
  data->SetPointError(0,0,445.513);
  data->SetMarkerStyle(29);
  data->SetMarkerSize(5);
  data->SetLineWidth(3);
  data->SetLineColor(6);
  data->SetMarkerColor(6);
  data->SetTitle("35-ton Data");
  mg->Add(data);
  TCanvas * canvelifeVnoise = new TCanvas("canvelifeVnoise","",2000,1600);
  canvelifeVnoise->SetLeftMargin(0.11);
  mg->Draw("alp");
  mg->GetXaxis()->SetTitle("Sim SNR / 35-ton SNR");
  mg->GetYaxis()->SetTitle("Measured eLifetime (us)");
  mg->GetYaxis()->SetTitleOffset(1.5);
  canvelifeVnoise->BuildLegend(0.7,0.7,0.9,0.9,"","lpe");
  TLine * line3 = new TLine(0,3000,10,3000);
  line3->SetLineStyle(kDashed);
  line3->SetLineColor(kBlack);
  line3->SetLineWidth(3);
  line3->Draw();
  TLine * line5 = new TLine(0,5000,10,5000);
  line5->SetLineStyle(kDashed);
  line5->SetLineColor(kBlack);
  line5->SetLineWidth(3);
  line5->Draw();
  TLine * line2 = new TLine(0,2000,10,2000);
  line2->SetLineStyle(kDashed);
  line2->SetLineColor(kBlack);
  line2->SetLineWidth(3);
  line2->Draw();
  canvelifeVnoise->Update();
  canvelifeVnoise->SaveAs("canvelifeVnoise.png");
}
