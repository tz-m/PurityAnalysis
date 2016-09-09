//#include "/home/mthiesse/Documents/PurityAnalysis/badchannel/PedCheck.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TTreeReader.h"

#include <string>
#include <map>
#include <algorithm>

void SNR(std::string fname)
{
  
  TFile * file = TFile::Open(fname.c_str(),"READ");
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

  TH1F * snr = new TH1F("snr","Ratio of signal-to-noise, pre- and post-filtering",200,0,0);
  
  while (reader.Next())
    {
      if (!(*countercut)) continue;
      if (!(*fitsuccess)) continue;
      if (!(*fitrealhit)) continue;
      if (!(((*c1>=6 && *c1<=15) || (*c1>=28 && *c1<=37)) && (*c1%22==*c2 || *c2%22==*c1))) continue;
      float snrfilt = *amplitudeFilter / *rmsFilter;
      float snrraw = *amplitude / *rms;
      snr->Fill(snrfilt/snrraw);
    }

  TCanvas * canv1 = new TCanvas("c1","c1",2000,1600);
  canv1->cd();
  snr->Draw();

}
