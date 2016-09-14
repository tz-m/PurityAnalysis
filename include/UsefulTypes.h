#ifndef USEFULTYPES_H
#define USEFULTYPES_H

namespace types {

  typedef std::vector<std::pair<Float_t, Float_t> > EfficiencyGraph; //vector of x position and efficiency
  typedef std::map<Int_t,EfficiencyGraph> ChannelEfficiencyMap;
  typedef std::map<Int_t,ChannelEfficiencyMap> TPCEfficiencyMap;

  struct HitInfo {
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

  typedef std::map<UInt_t,HitInfo> HitMap;

};

#endif
