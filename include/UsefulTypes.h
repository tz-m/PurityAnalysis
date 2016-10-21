#ifndef USEFULTYPES_H
#define USEFULTYPES_H

namespace types {

  typedef std::vector<std::pair<Float_t, Float_t> > EfficiencyGraph; //vector of x position and efficiency
  typedef std::map<Int_t,EfficiencyGraph> ChannelEfficiencyMap;
  typedef std::map<Int_t,ChannelEfficiencyMap> TPCEfficiencyMap;
  typedef std::map<Int_t,TPCEfficiencyMap> RunEfficiencyMap;

  struct HitInfo {
    Int_t run;
    Int_t event;
    Double_t t0;
    UInt_t c1;
    UInt_t c2;
    UInt_t trignum;
    Float_t c1x;
    Float_t c1y;
    Float_t c1z;
    Float_t c2x;
    Float_t c2y;
    Float_t c2z;
    Float_t distancecut;
    Int_t channel;
    Int_t wire;
    Int_t tpc;
    Int_t signalsize;
    std::vector<Float_t> signal;
    std::vector<Float_t> signalFilter;
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
    Bool_t assumedhit;
    Int_t numGoodHitsChan;
    Int_t nwiresTPC0;
    Int_t nwiresTPC1;
    Int_t nwiresTPC2;
    Int_t nwiresTPC3;
    Int_t nwiresTPC4;
    Int_t nwiresTPC5;
    Int_t nwiresTPC6;
    Int_t nwiresTPC7;
    Float_t prebaseline;
    Float_t postbaseline;
    Float_t prebaserms;
    Float_t postbaserms;
    Int_t trackid;
    Int_t numtrajpts;
    Double_t tracklength;
    Bool_t isOnTrack;
    Double_t dqdxatpt;
    Int_t prevStartTick;
    Int_t prevEndTick;
    Float_t prevPeakTime;
    Float_t prevSigmaPeakTime;
    Float_t prevRMS;
    Float_t prevPeakAmplitude;
    Float_t prevSigmaPeakAmplitude;
    Float_t prevSummedADC;
    Float_t prevIntegral;
    Float_t prevSigmaIntegral;
  };

  typedef std::map<UInt_t,HitInfo> HitMap;

  typedef std::map<std::string,Float_t> CutMap_F;
  typedef std::map<std::string,CutMap_F> AnalysisCuts_F;
  typedef std::map<std::string,Bool_t> CutMap_B;
  typedef std::map<std::string,CutMap_B> AnalysisCuts_B;

};

#endif
