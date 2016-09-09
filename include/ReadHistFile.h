#ifndef READHISTFILE_H
#define READHISTFILE_H

#include "TTree.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTreeReader.h"

#include <string>
#include <map>
#include <algorithm>

class HitInfo {
public:
  Int_t run;
  Int_t event;
  UInt_t c1;
  UInt_t c2;
  UInt_t trignum;
  Float_t c1x;
  Float_t c1y;
  Float_t c1z;
  Float_t c2x;
  Float_t c2y;
  Float_t c2z;
  Int_t channel;
  Int_t tpc;
  Int_t signalsize;
  Float_t baseline;
  Float_t rms;
  Float_t baselineFilter;
  Float_t rmsFilter;
  Float_t pedmean;
  Float_t pedrms;
  Float_t integral;
  Float_t integralFilter;
  Float_t sigmaintegral;
  Float_t sigmaintegralFilter;
  Float_t amplitude;
  Float_t amplitudeFilter;
  Float_t peaktick;
  Float_t peaktickFilter;
  Float_t peaktime;
  Float_t peaktimeFilter;
  Int_t begintick;
  Int_t endtick;
  Int_t width;
  Float_t hitx;
  Float_t hity;
  Float_t hitz;
  Float_t hiterrxlo;
  Float_t hiterrxhi;
  Float_t hiterrylo;
  Float_t hiterryhi;
  Float_t hiterrzlo;
  Float_t hiterrzhi;
  Float_t perpdist;
  Float_t hitt;
  Float_t driftdist;
  Bool_t countercut;
  Float_t fitconstant;
  Float_t fitconstanterr;
  Float_t fitlinear;
  Float_t fitlinearerr;
  Float_t fitquadratic;
  Float_t fitquadraticerr;
  Float_t fitchi2;
  Float_t fitsumsqrresidual;
  Float_t fitndf;
  Float_t fitmle;
  Bool_t fitsuccess;
  Bool_t fitrealhit;
  Float_t segmentlength;
};

class ReadHistFile {

public:

  ReadHistFile(std::string filename);

  HitInfo * GetHitInfo(unsigned int i) {
    return hits[i];
  }

private:
  std::map<unsigned int, HitInfo> hits;
};

HitValidation::HitValidation(std::string filename)
{
  TFile * file = TFile::Open(filename.c_str(),"READ");
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
  //TTreeReaderValue<Float_t> pedmean(reader,"pedmean");
  //TTreeReaderValue<Float_t> pedrms(reader,"pedrms");
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

  unsigned int i = 0;
  while (reader.Next())
    {
      HitInfo hi;
      hi.run = *run;
      hi.event = *event;
      hi.c1 = *c1;
      hi.c2 = *c2;
      hi.trignum = *trignum;
      hi.c1x = *c1x;
      hi.c1y = *c1y;
      hi.c1z = *c1z;
      hi.c2x = *c2x;
      hi.c2y = *c2y;
      hi.c2z = *c2z;
      hi.channel = *channel;
      hi.tpc = *tpc;
      hi.signalsize = *signalsize;
      hi.baseline = *baseline;
      hi.rms = *rms;
      hi.baselineFilter = *baselineFilter;
      hi.rmsFilter = *rmsFilter;
      //hi.pedmean = *pedmean;
      //hi.pedrms = *pedrms;
      hi.integral = *integral;
      hi.integralFilter = *integralFilter;
      hi.sigmaintegral = *sigmaintegral;
      hi.sigmaintegralFilter = *sigmaintegralFilter;
      hi.amplitude = *amplitude;
      hi.amplitudeFilter = *amplitudeFilter;
      hi.peaktick = *peaktick;
      hi.peaktickFilter = *peaktickFilter;
      hi.peaktime = *peaktime;
      hi.peaktimeFilter = *peaktimeFilter;
      hi.begintick = *begintick;
      hi.endtick = *endtick;
      hi.width = *width;
      hi.hitx = *hitx;
      hi.hity = *hity;
      hi.hitz = *hitz;
      hi.hiterrxlo = *hiterrxlo;
      hi.hiterrxhi = *hiterrxhi;
      hi.hiterrylo = *hiterrylo;
      hi.hiterryhi = *hiterryhi;
      hi.hiterrzlo = *hiterrzlo;
      hi.hiterrzhi = *hiterrzhi;
      hi.perpdist = *perpdist;
      hi.hitt = *hitt;
      hi.driftdist = *driftdist;
      hi.countercut = *countercut;
      hi.fitconstant = *fitconstant;
      hi.fitconstanterr = *fitconstanterr;
      hi.fitlinear = *fitlinear;
      hi.fitlinearerr = *fitlinearerr;
      hi.fitquadratic = *fitquadratic;
      hi.fitquadraticerr = *fitquadraticerr;
      hi.fitchi2 = *fitchi2;
      hi.fitsumsqrresidual = *fitsumsqrresidual;
      hi.fitndf = *fitndf;
      hi.fitmle = *fitmle;
      hi.fitsuccess = *fitsuccess;
      hi.fitrealhit = *fitrealhit;
      hi.segmentlength = *segmentlength;
      hits.emplace(std::make_pair(i,hi));
      i++;
    }
}

#endif
