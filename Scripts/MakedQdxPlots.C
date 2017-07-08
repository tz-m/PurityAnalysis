#include <algorithm>
#include "TString.h"
#include "TH1.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"
#include <sstream>
#include "TGaxis.h"

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
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooMinuit.h"
#include "RooAbsPdf.h"
#include "RooPolynomial.h"
#include "RooProduct.h"
#include "RooFormulaVar.h"

struct Histogram {
  TH1F * hist_data;
  TH1F * hist_sim;
  TH1F * hist_mc;
  double binlow;
  double binhigh;
  double bincenter;
  double binwidth;
  UInt_t binnum;

  double mpv_data;
  double mpv_data_err;
  double mpv_sim;
  double mpv_sim_err;
  double mpv_mc;
  double mpv_mc_err;
};

int getBinNumber( Float_t t, const std::vector<Histogram> & tbh )
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

void MakedQdxPlots()
{
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetLegendTextSize(0.05);

  double mcscale = 1.0;
  double elife = 4000;

  TString datafilename = "/home/mthiesse/PurityAnalysis/PurityAnalysis_data.root";
  TString simfilename = TString::Format("/home/mthiesse/PurityAnalysis/PurityAnalysis_%.fus_mcscale%.1f.root",elife,mcscale);

  TGaxis::SetMaxDigits(3);

  Double_t ADCcutofflow = 400;
  Double_t ADCcutoffhigh = 6000;
  Double_t ADCcutofflowMCT = 800;
  Double_t ADCcutoffhighMCT = 3000;
  std::vector<Histogram> timeBinHist(22);
  for (size_t i = 0; i < 22; ++i)
    {
      timeBinHist[i].binnum = i;
      timeBinHist[i].binlow = i*( 2012/22 );
      timeBinHist[i].binhigh = ( i+1 )*( 2012/22 );
      timeBinHist[i].bincenter = ( timeBinHist[i].binlow + timeBinHist[i].binhigh )/2.0;
      timeBinHist[i].binwidth = timeBinHist[i].binhigh - timeBinHist[i].binlow;
      timeBinHist[i].hist_data = new TH1F( TString::Format( "data_h%zu", i ),"", 200, ADCcutofflow, ADCcutoffhigh );
      timeBinHist[i].hist_sim = new TH1F( TString::Format( "sim_h%zu", i ),"", 200, ADCcutofflow, ADCcutoffhigh );
      timeBinHist[i].hist_mc = new TH1F( TString::Format( "mc_h%zu", i ),"", 600, ADCcutofflow, ADCcutoffhigh );
    }

  TFile * datafile = TFile::Open(datafilename,"READ");

  TTreeReader dataresults("results",datafile);
  TTreeReaderValue<Float_t> landmpv_data(dataresults,"landmpv");
  TTreeReaderValue<Float_t> landmpverr_data(dataresults,"landmpverr");
  TTreeReaderValue<Bool_t> fitsuccess_data(dataresults,"fitsuccess");
  TTreeReaderValue<Int_t> binnum2_data(dataresults,"binnum");

  while (dataresults.Next())
    {
      if (*fitsuccess_data && *binnum2_data != -1)
        {
          timeBinHist[*binnum2_data].mpv_data = *landmpv_data;
          timeBinHist[*binnum2_data].mpv_data_err = *landmpverr_data;
        }
    }

  TTreeReader datareader("hits",datafile);
  TTreeReaderValue<Float_t> hitt_data(datareader,"hitt");
  TTreeReaderValue<Float_t> dqdx_data(datareader,"dqdx");
  TTreeReaderValue<Int_t> channel_data(datareader,"channel");
  TTreeReaderValue<Bool_t> assumedhit_data(datareader,"assumedhit");
  TTreeReaderValue<Bool_t> fitrealhit_data(datareader,"foundrealhit");
  TTreeReaderValue<Int_t> tpc_data(datareader,"tpc");
  TTreeReaderValue<Int_t> binnum_data(datareader,"binnum");
  TTreeReaderValue<Bool_t> ismctruth_data(datareader,"ismctruth");

  while (datareader.Next())
    {
      if (*binnum_data == -1) continue;
      if (*assumedhit_data || *fitrealhit_data) timeBinHist[*binnum_data].hist_data->Fill(*dqdx_data);
    }

  delete datafile;

  TFile * simfile = TFile::Open(simfilename,"READ");

  TTreeReader simresults("results",simfile);
  TTreeReaderValue<Float_t> landmpv_sim(simresults,"landmpv");
  TTreeReaderValue<Float_t> landmpverr_sim(simresults,"landmpverr");
  TTreeReaderValue<Bool_t> fitsuccess_sim(simresults,"fitsuccess");
  TTreeReaderValue<Int_t> binnum2_sim(simresults,"binnum");

  while (simresults.Next())
    {
      if (*fitsuccess_sim && *binnum2_sim != -1)
        {
          timeBinHist[*binnum2_sim].mpv_sim = *landmpv_sim;
          timeBinHist[*binnum2_sim].mpv_sim_err = *landmpverr_sim;
        }
    }

  TTreeReader mcresults("resultsMCT",simfile);
  TTreeReaderValue<Float_t> landmpv_mc(mcresults,"landmpv");
  TTreeReaderValue<Float_t> landmpverr_mc(mcresults,"landmpverr");
  TTreeReaderValue<Bool_t> fitsuccess_mc(mcresults,"fitsuccess");
  TTreeReaderValue<Int_t> binnum2_mc(mcresults,"binnum");

  while (mcresults.Next())
    {
      if (*fitsuccess_mc && *binnum2_mc != -1)
        {
          timeBinHist[*binnum2_mc].mpv_mc = *landmpv_mc;
          timeBinHist[*binnum2_mc].mpv_mc_err = *landmpverr_mc;
        }
    }

  TTreeReader simreader("hits",simfile);
  TTreeReaderValue<Float_t> hitt_sim(simreader,"hitt");
  TTreeReaderValue<Float_t> dqdx_sim(simreader,"dqdx");
  TTreeReaderValue<Int_t> channel_sim(simreader,"channel");
  TTreeReaderValue<Bool_t> assumedhit_sim(simreader,"assumedhit");
  TTreeReaderValue<Bool_t> fitrealhit_sim(simreader,"foundrealhit");
  TTreeReaderValue<Int_t> tpc_sim(simreader,"tpc");
  TTreeReaderValue<Int_t> binnum_sim(simreader,"binnum");
  TTreeReaderValue<Bool_t> ismctruth_sim(simreader,"ismctruth");

  while (simreader.Next())
    {
      if (*binnum_sim == -1) continue;
      if (*ismctruth_sim) timeBinHist[*binnum_sim].hist_mc->Fill(*dqdx_sim+0.149*timeBinHist[*binnum_sim].mpv_mc);
      else timeBinHist[*binnum_sim].hist_sim->Fill(*dqdx_sim+0.149*timeBinHist[*binnum_sim].mpv_sim);
    }

  delete simfile;


  TMultiGraph * mg = new TMultiGraph();
  TGraphErrors * simmc = new TGraphErrors();
  simmc->SetName("simmc");
  simmc->SetMarkerStyle(20);
  simmc->SetMarkerSize(2);
  simmc->SetMarkerColor(kRed);
  TGraphErrors * datamc = new TGraphErrors();
  datamc->SetName("datamc");
  datamc->SetMarkerStyle(20);
  datamc->SetMarkerSize(2);
  datamc->SetMarkerColor(kBlue);
  Int_t point = 0;

  TGraph * mctruthmpv = new TGraph();

  TFile * fileout = TFile::Open("dqdxscale.root","RECREATE");

  std::stringstream binc;
  std::stringstream datampv;
  std::stringstream datampverr;
  std::stringstream simmpv;
  std::stringstream simmpverr;
  std::stringstream mcmpv;
  std::stringstream mcmpverr;

  for (auto hist : timeBinHist)
    {
      float maxval_data = hist.hist_data->GetBinContent(hist.hist_data->GetMaximumBin());
      int maxbin_data = hist.hist_data->GetBinCenter(hist.hist_data->GetMaximumBin());
      float min_data = hist.hist_data->GetBinCenter(hist.hist_data->FindFirstBinAbove(0.3*maxval_data));
      float max_data = hist.hist_data->GetBinCenter(hist.hist_data->FindLastBinAbove(0.5*maxval_data));
      TF1 * gaus_data = new TF1("gaus_data","gaus",min_data,max_data);
      //gaus_data->SetParLimits(1,500,7000);
      gaus_data->SetParameter(1,hist.hist_data->GetBinCenter(hist.hist_data->GetMaximumBin()));
      gaus_data->SetParameter(2,500);
      //hist.hist_data->Fit("gaus_data","RBQ");
      //hist.mpv_data = gaus_data->GetParameter(1);
      //hist.mpv_data_err = gaus_data->GetParError(1);

      float maxval_sim = hist.hist_sim->GetBinContent(hist.hist_sim->GetMaximumBin());
      float min_sim = hist.hist_sim->GetBinCenter(hist.hist_sim->FindFirstBinAbove(0.7*maxval_sim));
      float max_sim = hist.hist_sim->GetBinCenter(hist.hist_sim->FindLastBinAbove(0.7*maxval_sim));
      TF1 * gaus_sim = new TF1("gaus_sim","gaus",min_sim,max_sim);
      //gaus_sim->SetParLimits(1,500,7000);
      gaus_sim->SetParameter(1,hist.hist_sim->GetBinCenter(hist.hist_sim->GetMaximumBin()));
      gaus_sim->SetParameter(2,500);
      //hist.hist_sim->Fit("gaus_sim","RBQ");
      //hist.mpv_sim = gaus_sim->GetParameter(1);
      //hist.mpv_sim_err = gaus_sim->GetParError(1);

      float maxval_mc = hist.hist_mc->GetBinContent(hist.hist_mc->GetMaximumBin());
      float min_mc = hist.hist_mc->GetBinCenter(hist.hist_mc->FindFirstBinAbove(0.7*maxval_mc));
      float max_mc = hist.hist_mc->GetBinCenter(hist.hist_mc->FindLastBinAbove(0.7*maxval_mc));
      TF1 * gaus_mc = new TF1("gaus_mc","gaus",min_mc,max_mc);
      //gaus_mc->SetParLimits(1,500,7000);
      gaus_mc->SetParameter(1,hist.hist_mc->GetBinCenter(hist.hist_mc->GetMaximumBin()));
      gaus_mc->SetParameter(2,250);
      //hist.hist_mc->Fit("gaus_mc","RBQ");
      //hist.mpv_mc = gaus_mc->GetParameter(1);
      //hist.mpv_mc_err = gaus_mc->GetParError(1);

      double simmcratio = hist.mpv_sim/hist.mpv_mc;
      double simerr2 = TMath::Power(hist.mpv_sim_err/hist.mpv_sim,2);
      double mcerr2 = TMath::Power(hist.mpv_mc_err/hist.mpv_mc,2);
      simmc->SetPoint(point,hist.bincenter,simmcratio);
      simmc->SetPointError(point,0,simmcratio*sqrt(simerr2+mcerr2));

      double datamcratio = hist.mpv_data/hist.mpv_mc;
      double dataerr2 = TMath::Power(hist.mpv_data_err/hist.mpv_data,2);
      datamc->SetPoint(point,hist.bincenter,datamcratio);
      datamc->SetPointError(point,0,datamcratio*sqrt(dataerr2+mcerr2));

      if (hist.mpv_mc > 1) mctruthmpv->SetPoint(point,hist.bincenter,hist.mpv_mc);

      ++point;

      binc << hist.bincenter << ",";
      datampv << hist.mpv_data << ",";
      datampverr << hist.mpv_data_err << ",";
      simmpv << hist.mpv_sim << ",";
      simmpverr << hist.mpv_sim_err << ",";
      mcmpv << hist.mpv_mc << ",";
      mcmpverr << hist.mpv_mc_err << ",";

      std::cout << "BinCenter = " << hist.bincenter << std::endl;
      //std::cout << "MPV_Sim / MPV_MC = " << hist.mpv_sim / hist.mpv_mc << "  MPV_Data / MPV_MC = " << hist.mpv_data / hist.mpv_mc << std::endl;
      std::cout << "MPV_Sim = " << hist.mpv_sim << " MPV_Data = " << hist.mpv_data << " MPV_MC = " << hist.mpv_mc << std::endl;
      std::cout << "MPV_Sim_Error = " << hist.mpv_sim_err << " MPV_Data_Error = " << hist.mpv_data_err << " MPV_MC_Error = " << hist.mpv_mc_err << std::endl;

      //hist.hist_data->Scale(1/hist.hist_data->GetBinContent(hist.hist_data->GetMaximumBin()),"nosw2");
      //hist.hist_sim->Scale(1/hist.hist_sim->GetBinContent(hist.hist_sim->GetMaximumBin()),"nosw2");
      //hist.hist_mc->Scale(1/hist.hist_mc->GetBinContent(hist.hist_mc->GetMaximumBin()),"nosw2");
      hist.hist_data->Scale(1/hist.hist_data->Integral("width"),"nosw2");
      hist.hist_sim->Scale(1/hist.hist_sim->Integral("width"),"nosw2");
      hist.hist_mc->Scale(1/hist.hist_mc->Integral("width"),"nosw2");

      double maxamp = std::max(hist.hist_mc->GetBinContent(hist.hist_mc->GetMaximumBin()),std::max(hist.hist_data->GetBinContent(hist.hist_data->GetMaximumBin()),hist.hist_sim->GetBinContent(hist.hist_sim->GetMaximumBin())));

      TCanvas * canv = new TCanvas(TString::Format("canv_bin%i",getBinNumber(hist.bincenter,timeBinHist)),"",2000,1600);
      hist.hist_data->SetLineColor(kBlue);
      hist.hist_data->SetLineWidth(2);
      hist.hist_data->SetStats(false);
      hist.hist_data->SetFillColor(kBlue);
      hist.hist_data->SetFillColorAlpha(kBlue,0.35);
      hist.hist_data->Draw();
      hist.hist_data->GetYaxis()->SetTitle("Normalized Entries [cm/ADC]");
      hist.hist_data->GetXaxis()->SetTitle("dQ/dx [ADC/cm]");
      //hist.hist_data->GetYaxis()->SetLabelOffset(99);
      hist.hist_data->GetYaxis()->SetRangeUser(0,maxamp*1.05);
      hist.hist_sim->SetLineColor(kGreen);
      hist.hist_sim->SetLineWidth(2);
      hist.hist_sim->SetFillColor(kGreen);
      hist.hist_sim->SetFillColorAlpha(kGreen,0.35);
      hist.hist_sim->Draw("same");
      hist.hist_mc->SetLineColor(kCyan);
      hist.hist_mc->SetLineWidth(2);
      hist.hist_mc->SetFillColor(kCyan);
      hist.hist_mc->SetFillColorAlpha(kCyan,0.5);
      hist.hist_mc->Draw("same");
      TLegend * leg = new TLegend(0.45,0.65,0.95,0.95);
      leg->SetHeader(TString::Format("Drift Bin: %.f<t<%.f #mus",hist.binlow,hist.binhigh),"C");
      leg->AddEntry(hist.hist_data,"35-ton Data","fl");
      leg->AddEntry(hist.hist_sim,"35-ton Simulation","fl");
      leg->AddEntry(hist.hist_mc,"MCTruth","fl");
      leg->Draw();
      canv->Update();
      canv->SaveAs(TString::Format("canv_bin%i.png",getBinNumber(hist.bincenter,timeBinHist)));
      if (getBinNumber(hist.bincenter,timeBinHist) == 1) canv->SaveAs("short_dqdx.png");
      if (getBinNumber(hist.bincenter,timeBinHist) == 10) canv->SaveAs("middle_dqdx.png");
      if (getBinNumber(hist.bincenter,timeBinHist) == 19) canv->SaveAs("long_dqdx.png");
      delete canv;
    }

  std::cout << "Double_t binc[22] = " << binc.str() << ";" << std::endl;
  std::cout << "Double_t datampv[22] = " << datampv.str() << ";" << std::endl;
  std::cout << "Double_t datampverr[22] = " << datampverr.str() << ";" << std::endl;
  std::cout << "Double_t simmpv[22] = " << simmpv.str() << ";" << std::endl;
  std::cout << "Double_t simmpverr[22] = " << simmpverr.str() << ";" << std::endl;
  std::cout << "Double_t mcmpv[22] = " << mcmpv.str() << ";" << std::endl;
  std::cout << "Double_t mcmpverr[22] = " << mcmpverr.str() << ";" << std::endl;


  TCanvas * canvratios = new TCanvas("canvratios","",2000,1600);
  mg->Add(simmc,"pe");
  mg->Add(datamc,"pe");
  mg->Draw("a");
  mg->GetXaxis()->SetTitle("Drift Time (#mus)");
  //mg->GetYaxis()->SetTitle("MPV / MC Truth MPV");
  TLegend * leg = new TLegend(0.15,0.7,0.6,0.9);
  leg->AddEntry(simmc,"Simulation MPV / MCTruth MPV","pe");
  leg->AddEntry(datamc,"Data MPV / MCTruth MPV","pe");
  leg->Draw();
  canvratios->Update();
  canvratios->SaveAs("datasimmc_ratios.png");

  simmc->Write();
  datamc->Write();
  fileout->Close();

  TF1 * expo = new TF1( "expo", "[0]*exp( -x/[1] )", 0, 2012 );
  expo->SetParNames( "dQdx0", "eLifetime" );
  expo->SetParameters( 3000, 3000 );

  gStyle->SetOptStat( 0 );
  gStyle->SetOptFit( 1 );

  TCanvas * canvmcmpv = new TCanvas("canvmcmpv","",2000,1600);
  mctruthmpv->SetMarkerStyle(20);
  mctruthmpv->SetMarkerSize(2);
  mctruthmpv->Draw("ap");
  mctruthmpv->Fit("expo","Q");
  canvmcmpv->Update();
  canvmcmpv->SaveAs("mctruthmpv.png");
}
