#include "/home/mthiesse/Documents/BadChannel/PedCheck.h"
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
#include "TStyle.h"

#include <string>
#include <map>
#include <algorithm>

void HitValidation(std::string fname)
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

  TH1F * hitxpos = new TH1F("hitxpos","#splitline{Reconstructed hit locations, normalized by events in each muon counter pair}{RobustHitFinder threshold = 2.5*RMS}",200,-51,230);
  gStyle->SetTitleAlign(23);

  std::map<UInt_t,std::map<Int_t, std::vector<Int_t> > > cpairs;
  std::map<UInt_t,Float_t> counterx;
  std::map<Int_t,PedCheck> pedestalcheck;
  while (reader.Next())
    {
      if (!(*countercut)) continue;
      if (!(*fitsuccess)) continue;
      if (!(*fitrealhit)) continue;
      if (!(((*c1>=6 && *c1<=15) || (*c1>=28 && *c1<=37)) && (*c1%22==*c2 || *c2%22==*c1))) continue;
      if (std::find(cpairs[*c1][*run].begin(),cpairs[*c1][*run].end(),*event) == cpairs[*c1][*run].end()) cpairs[*c1][*run].push_back(*event);
      if (pedestalcheck.find(*run) == pedestalcheck.end())
	{
	  PedCheck pc(*run);
	  pc.setMinMaxPed(150,2000);
	  pc.setMinMaxRMS(6,40);
	  pedestalcheck.emplace(std::pair<Int_t,PedCheck>(*run,pc));
	}

      // if SIMULATION
      /*
      if (*c1==6 || *c1==28) counterx[*c1] = -33.5;
      if (*c1==7 || *c1==29) counterx[*c1] = -4;
      if (*c1==8 || *c1==30) counterx[*c1] = 25;
      if (*c1==9 || *c1==31) counterx[*c1] = 55;
      if (*c1==10 || *c1==32) counterx[*c1] = 85;
      if (*c1==11 || *c1==33) counterx[*c1] = 115;
      if (*c1==12 || *c1==34) counterx[*c1] = 145;
      if (*c1==13 || *c1==35) counterx[*c1] = 175;
      if (*c1==14 || *c1==36) counterx[*c1] = 205;
      if (*c1==15 || *c1==37) counterx[*c1] = 235;
      */
      
      // if DATA
      counterx[*c1] = *c1x;
    }

  std::map<UInt_t,Int_t> cpairEvents;
  for (auto const & cpair : cpairs)
    {
      std::cout << "Counter pair " << cpair.first;
      Int_t count = 0;
      for (auto const & run : cpair.second)
        {
          count += run.second.size();
        }
      cpairEvents[cpair.first] = count;
      std::cout << " contains " << count << " events." << std::endl;
    }

  reader.SetEntry(0);

  std::map<UInt_t,Float_t> cpairHits;
  while (reader.Next())
    {
      if (!(*countercut)) continue;
      if (!(*fitsuccess)) continue;
      if (!(*fitrealhit)) continue;
      if (!(((*c1>=6 && *c1<=15) || (*c1>=28 && *c1<=37)) && (*c1%22==*c2 || *c2%22==*c1))) continue;
      bool goodped = pedestalcheck[*run].check(*channel);
      if (!goodped) continue;
      Float_t bincontent = hitxpos->GetBinContent(hitxpos->FindBin(*hitx));
      bincontent += (cpairEvents[*c1] >= 1 ? 1.0/cpairEvents[*c1] : 0);
      hitxpos->SetBinContent(hitxpos->FindBin(*hitx),bincontent);
      cpairHits[*c1] += (cpairEvents[*c1] >= 1 ? 1.0/cpairEvents[*c1] : 0);
    }

  TGraph * grsum = new TGraph(cpairHits.size()-3);
  int i = 0;
  for (auto const & cpair : cpairHits)
    {
      if (cpair.first == 28 || cpair.first == 29 || cpair.first == 37) continue;
      grsum->SetPoint(i,counterx[cpair.first],cpair.second);
      i++;
      std::cout << "NHits in cpair " << cpair.first << " = " << cpair.second << std::endl;
    }
  TCanvas * canv1 = new TCanvas("canv1","canv1",2000,1600);
  TPad *pad1 = new TPad("hitxposhist","",0,0,1,1);
  TPad *pad2 = new TPad("graphintegral","",0,0,1,1);
  pad2->SetFillStyle(4000);
  pad1->Draw();
  pad1->cd();
  TPaveText * label = new TPaveText(0.45,0.85,0.55,0.9,"brNDC");
  label->AddText("Data");
  hitxpos->Draw();
  label->Draw();
  hitxpos->GetXaxis()->SetTitle("Drift distance (cm)");
  pad1->Update();
  gStyle->SetOptStat(0);
  //TPaveStats *ps1 = (TPaveStats*)hitxpos->GetListOfFunctions()->FindObject("stats");
  //ps1->SetX1NDC(0.65); ps1->SetX2NDC(0.85);
  pad1->Modified();
  canv1->cd();
  Double_t ymin = 0;
  Double_t ymax = 450;
  Double_t dy = (ymax-ymin)/0.8;
  Double_t xmin = -51;
  Double_t xmax = 230;
  Double_t dx = (xmax-xmin)/0.8;
  pad2->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);
  pad2->Draw();
  pad2->cd();
  grsum->SetLineColor(kRed);
  grsum->SetLineWidth(2);
  grsum->SetMarkerStyle(20);
  grsum->SetMarkerSize(2);
  grsum->SetMarkerColor(kRed);
  grsum->Draw("][samelp");
  pad2->Update();
  grsum->GetListOfFunctions()->Print();
  TGaxis *axis = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,50510,"+L");
  axis->SetLabelColor(kRed);
  axis->SetTitle("Total hits per counter pair / Events per counter pair");
  axis->SetTitleSize(0.03);
  axis->SetTitleOffset(1.6);
  axis->Draw();
  gROOT->GetListOfCanvases()->Draw();
  canv1->WaitPrimitive();


}
