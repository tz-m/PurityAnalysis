#ifndef PURITYANALYSIS_H
#define PURITYANALYSIS_H

#include "Analysis.h"

#include "THStack.h"
#include "TFitResult.h"
#include "TSystem.h"
#include "TImage.h"
#include "TPaveLabel.h"
#include "TVirtualFitter.h"
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
#include "Math/VavilovAccurate.h"
#include "RooFunctorBinding.h"
#include "Math/Functor.h"

struct ChargeHistogram
{
  TH1F * hist;
  Float_t binlow;
  Float_t bincenter;
  Float_t binhigh;
  Float_t binwidth;
};

struct FitResult
{
  Int_t binnum;
  Float_t meanlandau;
  Float_t widthlandau;
  Float_t meangauss;
  Float_t meangaussvav;
  Float_t widthgauss;
  Float_t widthgaussvav;
  Float_t kappavav;
  Float_t beta2vav;
  Float_t meanvav;
  Float_t sigmavav;
  Float_t ampvav;
  Float_t meanlandauerr;
  Float_t widthlandauerr;
  Float_t meangausserr;
  Float_t widthgausserr;
  Float_t meangaussvaverr;
  Float_t widthgaussvaverr;
  Float_t kappavaverr;
  Float_t beta2vaverr;
  Float_t meanvaverr;
  Float_t sigmavaverr;
  Float_t ampvaverr;
  Float_t bincenter;
  Float_t binwidtherr;
  Float_t lxg_chi2ndf;
  Float_t vxg_chi2ndf;
};


struct Vavilov_Func {
  Vavilov_Func() {
  }

  ROOT::Math::VavilovAccurate pdf;

  double operator() (const double *x) {
    return x[5]*( pdf.Pdf( (x[0]-x[3])/x[4],x[1],x[2]));
  }
};


class PurityAnalysis : public Analysis {
public:
  PurityAnalysis();
  void run();

  //void reconfigureCuts();

private:

  int getBinNumber(Float_t t,const std::vector<ChargeHistogram> & tbh);

};

PurityAnalysis::PurityAnalysis()
{
}

int PurityAnalysis::getBinNumber(Float_t t,const std::vector<ChargeHistogram> & tbh)
{
  for (UInt_t i_hist = 0; i_hist < tbh.size(); i_hist++)
    {
      if (t >= tbh[i_hist].binlow && t < tbh[i_hist].binhigh)
        {
          return (int)i_hist;
        }
    }
  return -1;
}

void PurityAnalysis::run()
{
  //reconfigureCuts();
  TFile * fileout = TFile::Open("PurityAnalysis.root","RECREATE");

  Int_t runnum;
  //Int_t subrun;
  Int_t event;
  Float_t hitt;
  Int_t binnum;
  Bool_t foundrealhit;
  Bool_t assumedhit;
  Float_t dqdx;
  Int_t channel;
  Int_t tpc;
  TTree * hitstree = new TTree("hits","Hit Information");
  hitstree->Branch("runnum",&runnum,"runnum/I");
  //hitstree->Branch("subrun",&subrun,"subrun/I");
  hitstree->Branch("event",&event,"event/I");
  hitstree->Branch("hitt",&hitt,"hitt/F");
  hitstree->Branch("binnum",&binnum,"binnum/I");
  hitstree->Branch("foundrealhit",&foundrealhit,"foundrealhit/O");
  hitstree->Branch("assumedhit",&assumedhit,"assumedhit/O");
  hitstree->Branch("dqdx",&dqdx,"dqdx/F");
  hitstree->Branch("channel",&channel,"channel/I");
  hitstree->Branch("tpc",&tpc,"tpc/I");

  UInt_t Nbins = 22;
  Float_t drifttimemax = 2012;
  Float_t ADCcutofflow = -3500;
  Float_t ADCcutoffhigh = 20000;

  TH1F * found = new TH1F("found","Hits Found",200,ADCcutofflow,ADCcutoffhigh);
  TH1F * assumed = new TH1F("assumed","Assumed Hits",200,ADCcutofflow,ADCcutoffhigh);

  RooRealVar charge("charge","Hit Integral (ADC*tick)",ADCcutofflow,ADCcutoffhigh);
  RooRealVar tdrift("tdrift","tdrift",0,drifttimemax);
  RooDataSet * chgdata = new RooDataSet("chgdata","chgdata",RooArgSet(tdrift,charge));

  std::vector<ChargeHistogram> timeBinHist(Nbins);
  for (UInt_t i = 0; i < Nbins; i++)
    {
      timeBinHist[i].binlow = i*(drifttimemax/Nbins);
      timeBinHist[i].binhigh = (i+1)*(drifttimemax/Nbins);
      timeBinHist[i].bincenter = (timeBinHist[i].binlow + timeBinHist[i].binhigh)/2.0;
      timeBinHist[i].binwidth = timeBinHist[i].binhigh - timeBinHist[i].binlow;
      timeBinHist[i].hist = new TH1F(TString::Format("h%i",i),TString::Format("RobustHitFinder, Bin %i, %.2f <= t < %.2f (us)",i,timeBinHist[i].binlow,timeBinHist[i].binhigh),100,0,drifttimemax);
    }

  Float_t secondbinstart = timeBinHist[1].binlow;
  Float_t penultimatebinend = timeBinHist[Nbins-2].binhigh;

  TH1F * hlong = new TH1F("hlong","Hit Integral (ADC)",1000,ADCcutofflow,ADCcutoffhigh);

  std::map<Int_t, std::vector<std::pair<Float_t,Float_t> > > hitmap;

  for (auto const & hititr : *(file->GetHitMap()))
    {
      const types::HitInfo * hit = &(hititr.second);
      if (hit->run >= 15483 && hit->run <= 15502) continue;
      if (hit->run == 15593) continue;
      if (hit->run >= 15634 && hit->run <= 15643) continue;
      if (hit->run >= 15664 && hit->run <= 15674) continue;
      if (hit->run == 15823 || hit->run == 15914 || hit->run == 15940 || hit->run == 15980 || hit->run == 16025 || hit->run == 16083 || hit->run == 16147 || hit->run == 16586) continue;
      if (hit->run >= 16593 && hit->run <= 16596) continue;
      if (hit->channel == 566 || hit->channel == 885 || hit->channel == 1547) continue;
      if (!(hit->fitsuccess)) continue;
      //if (hit->segmentlength > 2 || hit->segmentlength < 0.5) continue;
      if (hit->tpc != 1) continue;
      //if (hit->width > 100) continue;
      //if (hit->pedmean < 200 || hit->pedmean > 1500) continue;
      //if (hit->pedrms < 7 || hit->pedrms > 30) continue;
      if (cuts.ChannelPass(hit) &&
          cuts.CounterPass(hit)) //&&
      //cuts.HitPass(hit))
        {
          runnum = hit->run;
          //subrun = hit->subrun;
          event = hit->event;
          hitt = hit->hitt;
          binnum = getBinNumber(hitt,timeBinHist);
          foundrealhit = !hit->assumedhit && hit->fitrealhit;
          assumedhit = hit->assumedhit;
          dqdx = hit->integral / hit->segmentlength;
          channel = hit->channel;
          tpc = hit->tpc;
          hitstree->Fill();

          if (foundrealhit) found->Fill(dqdx);
          if (assumedhit) assumed->Fill(dqdx);
          if (assumedhit || hit->fitrealhit) hitmap[channel].push_back(std::make_pair(hitt,dqdx));
        }
    }
  hitstree->Write();

  TCanvas * ctest = new TCanvas("ctest","ctest",2000,1600);
  ctest->cd();
  THStack * hs = new THStack("hs","Found and Assumed Hits");
  found->SetFillColor(kRed);
  assumed->SetFillColor(kBlue);
  hs->Add(found);
  hs->Add(assumed);
  hs->Draw();
  ctest->Write();

  gStyle->SetOptStat(11);
  TCanvas * canvtest = new TCanvas("fahist","Found and Assumed Hits",2000,1600);
  TPad * pad1 = new TPad("pad1","",0,0,1,1);
  TPad * pad2 = new TPad("pad2","",0,0,1,1);
  pad2->SetFillStyle(4000);
  pad1->Draw();
  pad1->cd();
  Double_t maxval = std::max(found->GetBinContent(found->GetMaximumBin()),assumed->GetBinContent(assumed->GetMaximumBin()));
  found->SetLineColor(kRed);
  found->SetLineWidth(2);
  found->SetAxisRange(0,maxval*1.1,"Y");
  found->Draw();
  found->GetXaxis()->SetTitle("Charge");
  pad1->Update();
  TPaveStats * ps1 = (TPaveStats*)found->GetListOfFunctions()->FindObject("stats");
  ps1->SetX1NDC(0.4); ps1->SetX2NDC(0.6);
  ps1->SetTextColor(kRed);
  pad1->Modified();
  canvtest->cd();
  Double_t ymin = 0;
  Double_t ymax = maxval*1.1;
  Double_t dy = (ymax-ymin)/0.8;
  Double_t xmin = ADCcutofflow;
  Double_t xmax = ADCcutoffhigh;
  Double_t dx = (xmax-xmin)/0.8;
  pad2->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);
  pad2->Draw();
  pad2->cd();
  assumed->SetLineColor(kBlue);
  assumed->SetLineWidth(2);
  assumed->Draw("][sames");
  assumed->GetXaxis()->SetTitle("Charge");
  pad2->Update();
  TPaveStats * ps2 = (TPaveStats*)assumed->GetListOfFunctions()->FindObject("stats");
  ps2->SetX1NDC(0.65); ps2->SetX2NDC(0.85);
  ps2->SetTextColor(kBlue);
  TGaxis * axis = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,510,"+L");
  axis->SetLabelFont(42);
  axis->SetLabelSize(0.035);
  axis->Draw();
  TLegend * leg = new TLegend(0.55,0.6,0.85,0.75);
  leg->AddEntry(found,"Found Hits","l");
  leg->AddEntry(assumed,"Assumed Hits","l");
  leg->Draw();
  canvtest->Update();
  canvtest->Write();

  for (auto & i_chan : hitmap)
    {
      for (auto & i_hit : i_chan.second)
        {
          if (i_hit.second < ADCcutoffhigh && i_hit.second > ADCcutofflow)
            {
              hlong->Fill(i_hit.second);

              charge = (Float_t)(i_hit.second);
              tdrift = (Float_t)(i_hit.first);
              chgdata->add(RooArgSet(tdrift,charge));
              for (UInt_t i_hist = 0; i_hist < timeBinHist.size(); i_hist++)
                {
                  if (i_hit.first >= timeBinHist[i_hist].binlow && i_hit.first < timeBinHist[i_hist].binhigh)
                    {
                      timeBinHist[i_hist].hist->Fill(i_hit.second);
                      break;
                    }
                }
            }
        }
    }



  std::vector<FitResult> results(Nbins);

  Vavilov_Func vav;
  ROOT::Math::Functor func(vav,6);

  RooRealVar kappa("VavilovKappa","vavilov kappa",0.1,0,10);
  RooRealVar beta2("VavilovBeta2","vavilov beta2",0.4,0,1);
  RooRealVar mv("VavilovMean","vavilov mean",2800,-1000,20000);
  RooRealVar sv("VavilovSigma","vavilov sigma",1000,0,10000);
  RooRealVar av("VavilovAmplitude","vavilov amp",1e4,1e0,1e9);
  RooAbsPdf * vavilov = new RooFunctorPdfBinding("vavbind","vavbind",func,RooArgList(charge,kappa,beta2,mv,sv,av));
  RooRealVar mgv("mgv","mean gauss vav",0);
  RooRealVar sgv("sgv","sigma gauss vav",50,0,2000);
  RooGaussian gaussvav("gaussvav","gaussvav",charge,mgv,sgv);
  RooRealVar mgv2("mgv2","mean gauss vav 2",5000,2800,8000);
  RooRealVar sgv2("sgv2","sigma gauss vav 2",50,0,2000);
  RooGaussian deltavav("deltavav","deltavav",charge,mgv2,sgv2);

  RooRealVar ml("LandMPV","mean landau",2800,-1000,20000);
  RooRealVar sl("LandWidth","sigma landau",1000,0,10000);
  RooLandau landau("lx","lx",charge,ml,sl);
  RooRealVar mg("mg","mean gauss",0);
  RooRealVar sg("GaussWidth","sigma gauss",50,0,2000);
  RooGaussian gauss("gauss","gauss",charge,mg,sg);
  RooRealVar mg2("mg2","mean gauss delta",5000,2800,8000);
  RooRealVar sg2("sg2","sigma gauss delta",50,0,2000);
  RooGaussian delta("delta","delta",charge,mg2,sg2);
  RooFFTConvPdf lxg("lxg","landau (x) gauss",charge,landau,gauss);
  RooFFTConvPdf vxg("vxg","vavilov (x) gauss",charge,*vavilov,gaussvav);

  RooRealVar coef1("coef1","Coefficient delta",0.1,0,0.7);
  RooAddPdf lxgxg("lxgxg","(lxg)+g",lxg,delta,coef1);
  RooRealVar coef2("coef2","Coefficient delta 2",0.1,0,0.7);
  RooAddPdf vxgxg("vxgxg","(vxg)+g",vxg,deltavav,coef2);

  TGraphAsymmErrors * hbin = new TGraphAsymmErrors(Nbins);
  hbin->SetTitle("Most Probable dQ/dx");
  Double_t minMPV = 99999;
  Double_t maxMPV = -99999;

  TGraphAsymmErrors * gw = new TGraphAsymmErrors(Nbins);
  gw->SetTitle("Width of Gaussian in Convolution");
  TGraphAsymmErrors * lw = new TGraphAsymmErrors(Nbins);
  lw->SetTitle("Width of Landau in Convolution");

  for (UInt_t i = 0; i < Nbins; i++)
    {
      if (timeBinHist[i].hist->GetEntries() > 0)
        {
          char shortname[100];
          sprintf(shortname,"dQdx_bin%d",i);

          Float_t maxbin = timeBinHist[i].hist->GetBinCenter(timeBinHist[i].hist->GetMaximumBin());
          //Float_t maxbin = 2000;

          // Setup component pdfs
          // --------------------

          //setup observable
          //charge.setRange(shortname,maxbin-2000,maxbin+3000);
          charge.setRange(shortname,ADCcutofflow,ADCcutoffhigh);

          //setup landau(t,ml,sl)
          ml.setVal(maxbin);
          ml.setMin(maxbin-2000); //0);
          ml.setMax(maxbin+2000); //0);
          sl.setVal(850);
          sl.setMin(60);
          sl.setMax(2000);

          //setup gauss(t,mg,sg)
          //mg.setVal(1);
          //mg.setMin(-100);
          //mg.setMax(100);
          sg.setVal(550);
          sg.setMin(10);
          sg.setMax(2000);

          // construct convolution pdf
          // -------------------------

          // set num bins to be used for FFT sampling
          charge.setBins(10000,"fft");

          // fit convoluted pdf to binHist
          // -----------------------------

          char cut[100];
          sprintf(cut,"tdrift > %f && tdrift < %f",timeBinHist[i].binlow,timeBinHist[i].binhigh);
          RooDataSet* roodata = (RooDataSet*) chgdata->reduce(RooFit::SelectVars(charge),RooFit::Cut(cut),RooFit::Name(shortname),RooFit::Title(shortname));
          //RooDataHist* roodata = new RooDataHist("roodata","roodata histogram",RooArgSet(charge),*roodataset);


          //RooAbsReal* nll = lxg.createNLL(*roodata,RooFit::NumCPU(4));
          //RooMinuit m(*nll);
          //m.migrad();
          //m.hesse();
          //RooFitResult * rooresult = m.save();
          //RooFitResult * rooresult =
          vxgxg.fitTo(*roodata,RooFit::Save(),RooFit::PrintLevel(-1),RooFit::Range(shortname));
          lxgxg.fitTo(*roodata,RooFit::Save(),RooFit::PrintLevel(-1),RooFit::Range(shortname));
          RooPlot * frame = charge.frame(RooFit::Title(timeBinHist[i].hist->GetTitle()));
          roodata->plotOn(frame,RooFit::Binning(200));
          lxgxg.plotOn(frame,RooFit::Name("lxg"),RooFit::LineColor(kBlue),RooFit::ShiftToZero());
          results[i].lxg_chi2ndf = frame->chiSquare();
          lxg.plotOn(frame,RooFit::LineColor(kBlue),RooFit::ShiftToZero(),RooFit::LineStyle(kDashed));
          vxgxg.plotOn(frame,RooFit::Name("vxg"),RooFit::LineColor(kRed),RooFit::ShiftToZero());
          results[i].vxg_chi2ndf = frame->chiSquare();
          vxg.plotOn(frame,RooFit::LineColor(kRed),RooFit::ShiftToZero(),RooFit::LineStyle(kDashed));

          TPaveLabel * tl = new TPaveLabel(0.6,0.82,0.99,0.9,Form(" #chi^{2}/NDF = %f",frame->chiSquare()),"NDC");
          tl->SetFillColor(0);
          tl->SetBorderSize(1);
          tl->SetTextAlign(12);
          frame->addObject(tl);
          frame->getAttText()->SetTextSize(0.32);
          roodata->statOn(frame,RooFit::Layout(0.6,0.99,0.8));
          frame->getAttText()->SetTextSize(0.025);
          vxgxg.paramOn(frame,RooFit::Layout(0.6,0.99,0.6));
          frame->getAttText()->SetTextSize(0.025);
          lxgxg.paramOn(frame,RooFit::Layout(0.6,0.99,0.3));
          frame->getAttText()->SetTextSize(0.025);



          char buff[100];
          sprintf(buff,"bin%d",i);
          TCanvas* c = new TCanvas(buff,buff,1600,1600);
          gPad->SetLeftMargin(0.15);
          frame->GetYaxis()->SetTitleOffset(1.4);
          frame->Draw();
          c->Write();
          /*
             gSystem->ProcessEvents();
             TImage *img = TImage::Create();
             img->FromPad(c);
             char newbuff[100];
             sprintf(newbuff,"%s.png",buff);
             img->WriteImage(newbuff);
           */

          char buffnew[100];
          sprintf(buffnew,"bin%d_resid",i);
          TCanvas* cresid = new TCanvas(buffnew,buffnew,1600,1600);
          gPad->SetLeftMargin(0.15);
          RooPlot* frame2 = charge.frame(RooFit::Title("residual"));
          RooHist* hresid = frame->residHist();
          frame2->addPlotable(hresid,"P");
          frame2->GetYaxis()->SetTitleOffset(1.4);
          frame2->Draw();
          cresid->Write();
          /*
             gSystem->ProcessEvents();
             TImage *img2 = TImage::Create();
             img2->FromPad(cresid);
             char buffnew2[100];
             sprintf(buffnew2,"%s.png",buffnew);
             img2->WriteImage(buffnew2);
           */

          char buffpull[100];
          sprintf(buffpull,"bin%d_pull",i);
          TCanvas* cpull = new TCanvas(buffpull,buffpull,1600,1600);
          gPad->SetLeftMargin(0.15);
          RooPlot* frame3 = charge.frame(RooFit::Title("pull"));
          RooHist* hpull2 = frame->pullHist();
          TH1F* hpull = new TH1F("hpull","pull",100,0,0);
          for (Int_t ibin=0; ibin<hpull2->GetN(); ibin++)
            {
              Double_t x,y;
              hpull2->GetPoint(ibin,x,y);
              hpull->Fill(y);
            }
          hpull->Draw();
          cpull->Write();
          /*
             gSystem->ProcessEvents();
             TImage *img3 = TImage::Create();
             img3->FromPad(cpull);
             char buffpull2[100];
             sprintf(buffpull2,"%s.png",buffpull);
             img3->WriteImage(buffpull2);
           */



          results[i].binnum = i;
          results[i].meanlandau = ml.getVal();
          results[i].widthlandau = sl.getVal();
          results[i].meangauss = mg.getVal();
          results[i].widthgauss = sg.getVal();
          results[i].meangaussvav = mgv.getVal();
          results[i].widthgaussvav = sgv.getVal();
          results[i].kappavav = kappa.getVal();
          results[i].beta2vav = beta2.getVal();
          results[i].meanvav = mv.getVal();
          results[i].sigmavav = sv.getVal();
          results[i].ampvav = av.getVal();
          results[i].meanlandauerr = ml.getError();
          results[i].widthlandauerr = sl.getError();
          results[i].meangausserr = mg.getError();
          results[i].widthgausserr = sg.getError();
          results[i].meangaussvaverr = mgv.getError();
          results[i].widthgaussvaverr = sgv.getError();
          results[i].kappavaverr = kappa.getError();
          results[i].beta2vaverr = beta2.getError();
          results[i].meanvaverr = mv.getError();
          results[i].sigmavaverr = sv.getError();
          results[i].ampvaverr = av.getError();
          results[i].bincenter = timeBinHist[i].bincenter;
          results[i].binwidtherr = (timeBinHist[i].binwidth)/sqrt(12.0);


          //lxg_chi2ndf = frame->chiSquare("lxg",shortname);
          //vxg_chi2ndf = frame->chiSquare("vxg",shortname);
          //resulttree->Fill();

          if (minMPV > results[i].meanlandau) minMPV = results[i].meanlandau;
          if (maxMPV < results[i].meanlandau) maxMPV = results[i].meanlandau;

          hbin->SetPoint(i,timeBinHist[i].bincenter,ml.getVal());
          hbin->SetPointError(i,(timeBinHist[i].binwidth)/sqrt(12.0), (timeBinHist[i].binwidth)/sqrt(12.0),ml.getError(),ml.getError());
          gw->SetPoint(i,timeBinHist[i].bincenter,sg.getVal());
          gw->SetPointError(i,(timeBinHist[i].binwidth)/sqrt(12.0),(timeBinHist[i].binwidth)/sqrt(12.0),sg.getError(),sg.getError());
          lw->SetPoint(i,timeBinHist[i].bincenter,sl.getVal());
          lw->SetPointError(i,(timeBinHist[i].binwidth)/sqrt(12.0),(timeBinHist[i].binwidth)/sqrt(12.0),sl.getError(),sl.getError());

          //c->WaitPrimitive();
          delete c;
          delete cresid;
          delete cpull;
          delete frame;
          delete frame2;
          delete frame3;
          //delete img;
          //delete img2;
          //delete img3;

        }
    }

  for (UInt_t i = 0; i < Nbins; i++)
    {
      std::cout << "Bin " << i << ": ml=" << results[i].meanlandau << "  sl=" << results[i].widthlandau << "  sg=" << results[i].widthgauss << "  chi2ndf=" << results[i].lxg_chi2ndf << std::endl;
      std::cout << "      : mv=" << results[i].meanvav << "  sv=" << results[i].sigmavav << "  sgv=" << results[i].widthgaussvav << " chi2ndf=" << results[i].vxg_chi2ndf << std::endl;
    }

  TF1 * expo = new TF1("expo","[0]*exp(-x/[1])",secondbinstart,penultimatebinend);
  expo->SetParNames("dQdx0","eLifetime");
  expo->SetParameters(3000,3000);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  TCanvas * canv2 = new TCanvas("canv2","canv2",1600,800);
  canv2->cd();
  hbin->Fit("expo","R");
  hbin->Draw("ap");
  hbin->SetMarkerStyle(kFullDotLarge);
  hbin->SetMarkerSize(2);

  TH1F * hint = new TH1F("hint","95#% confidence band",100,0,drifttimemax);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint,0.68);
  hint->SetStats(kFALSE);
  hint->SetFillColorAlpha(kRed,0.35);
  hint->Draw("e3 same");

  hbin->GetYaxis()->SetRangeUser(minMPV-100,maxMPV+100);
  hbin->GetXaxis()->SetTitle("Drift Time (#mu s)");
  hbin->GetYaxis()->SetTitle("dQ/dx (ADC/cm)");

  canv2->Write();

  TCanvas * canv3 = new TCanvas("canv3","canv3",1600,800);
  canv3->cd();
  gw->Draw("ap");
  gw->SetMarkerStyle(kFullDotLarge);
  gw->SetMarkerSize(2);
  gw->GetXaxis()->SetTitle("Drift Time (#mu s)");

  canv3->Write();

  TCanvas * canv4 = new TCanvas("canv4","canv4",1600,800);
  canv4->cd();
  lw->Draw("ap");
  lw->SetMarkerStyle(kFullDotLarge);
  lw->SetMarkerSize(2);
  lw->GetXaxis()->SetTitle("Drift Time (#mu s)");

  canv4->Write();
/*
   TCanvas * canv3 = new TCanvas("canv3","canv3",1600,800);
   canv3->cd();
   hlong->Draw();
 */

  TTree * resulttree = new TTree("results","Results of fits");
  Float_t meanlandau;
  Float_t widthlandau;
  Float_t meangauss;
  Float_t widthgauss;
  Float_t meangaussvav;
  Float_t widthgaussvav;
  Float_t kappavav;
  Float_t beta2vav;
  Float_t meanvav;
  Float_t sigmavav;
  Float_t ampvav;
  Float_t meanlandauerr;
  Float_t widthlandauerr;
  Float_t meangausserr;
  Float_t widthgausserr;
  Float_t meangaussvaverr;
  Float_t widthgaussvaverr;
  Float_t kappavaverr;
  Float_t beta2vaverr;
  Float_t meanvaverr;
  Float_t sigmavaverr;
  Float_t ampvaverr;
  Float_t bincenter;
  Float_t binwidtherr;
  Float_t lxg_chi2ndf;
  Float_t vxg_chi2ndf;
  resulttree->Branch("binnum",&binnum,"binnum/I");
  resulttree->Branch("meanlandau",&meanlandau,"meanlandau/F");
  resulttree->Branch("widthlandau",&widthlandau,"widthlandau/F");
  resulttree->Branch("meangauss",&meangauss,"meangauss/F");
  resulttree->Branch("widthgauss",&widthgauss,"widthgauss/F");
  resulttree->Branch("meangaussvav",&meangaussvav,"meangaussvav/F");
  resulttree->Branch("widthgaussvav",&widthgaussvav,"widthgaussvav/F");
  resulttree->Branch("kappavav",&kappavav,"kappavav/F");
  resulttree->Branch("beta2vav",&beta2vav,"beta2vav/F");
  resulttree->Branch("meanvav",&meanvav,"meanvav/F");
  resulttree->Branch("sigmavav",&sigmavav,"sigmavav/F");
  resulttree->Branch("ampvav",&ampvav,"ampvav/F");
  resulttree->Branch("meanlandauerr",&meanlandauerr,"meanlandauerr/F");
  resulttree->Branch("widthlandauerr",&widthlandauerr,"widthlandauerr/F");
  resulttree->Branch("meangausserr",&meangausserr,"meangausserr/F");
  resulttree->Branch("widthgausserr",&widthgausserr,"widthgausserr/F");
  resulttree->Branch("meangaussvaverr",&meangaussvaverr,"meangaussvaverr/F");
  resulttree->Branch("widthgaussvaverr",&widthgaussvaverr,"widthgaussvaverr/F");
  resulttree->Branch("kappavaverr",&kappavaverr,"kappavaverr/F");
  resulttree->Branch("beta2vaverr",&beta2vaverr,"beta2vaverr/F");
  resulttree->Branch("meanvaverr",&meanvaverr,"meanvaverr/F");
  resulttree->Branch("sigmavaverr",&sigmavaverr,"sigmavaverr/F");
  resulttree->Branch("ampvaverr",&ampvaverr,"ampvaverr/F");
  resulttree->Branch("bincenter",&bincenter,"bincenter/F");
  resulttree->Branch("binwidtherr",&binwidtherr,"binwidtherr/F");
  resulttree->Branch("lxg_chi2ndf",&lxg_chi2ndf,"lxg_chi2ndf/F");
  resulttree->Branch("vxg_chi2ndf",&vxg_chi2ndf,"vxg_chi2ndf/F");

  for (auto const & r : results)
    {
      binnum = r.binnum;
      meanlandau = r.meanlandau;
      widthlandau = r.widthlandau;
      meangauss = r.meangauss;
      widthgauss = r.widthgauss;
      meangaussvav = r.meangaussvav;
      widthgaussvav = r.widthgaussvav;
      kappavav = r.kappavav;
      beta2vav = r.beta2vav;
      meanvav = r.meanvav;
      sigmavav = r.sigmavav;
      ampvav = r.ampvav;
      meanlandauerr = r.meanlandauerr;
      widthlandauerr = r.widthlandauerr;
      meangausserr = r.meangausserr;
      widthgausserr = r.widthgausserr;
      meangaussvaverr = r.meangaussvaverr;
      widthgaussvaverr = r.widthgaussvaverr;
      kappavaverr = r.kappavaverr;
      beta2vaverr = r.beta2vaverr;
      meanvaverr = r.meanvaverr;
      sigmavaverr = r.sigmavaverr;
      ampvaverr = r.ampvaverr;
      bincenter = r.bincenter;
      binwidtherr = r.binwidtherr;
      lxg_chi2ndf = r.lxg_chi2ndf;
      vxg_chi2ndf = r.vxg_chi2ndf;
      resulttree->Fill();
    }

  std::cout << "Found Number = " << found->GetEntries() << "  Assumed Number = " << assumed->GetEntries() << std::endl;


  resulttree->Write();
  fileout->Close();
  return;
}

#endif
