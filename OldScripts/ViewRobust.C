#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TLine.h"

#include <string>
#include <map>
#include <algorithm>

void ViewRobust()
{
  TFile * file = TFile::Open("/home/mthiesse/Documents/PurityAnalysis/DataFiles/full_100-4400_split_minprep_reco_1.5_robust_hist.root","READ");
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
  Int_t nevents = 0;
  
  std::map<Int_t,std::vector<Int_t> > runevent;
  while (reader.Next())
    {
      if (*c1==5 || *c2==5 || *c1==43 || *c2==43) continue;
      if (std::count(runevent[*run].begin(),runevent[*run].end(),*event)==0)
	{
	  nevents++;
	  runevent[*run].push_back(*event);
	}
    }
  reader.SetEntry(0);

  TCanvas * canv1 = new TCanvas("canv1","canv1",1200,800);
  for (auto i_run = runevent.begin(); i_run != runevent.end(); i_run++)
    {
      for (auto i_evt = (i_run->second).begin(); i_evt != (i_run->second).end(); i_evt++)
	{
	  std::cout << "Run=" << (i_run->first) << "  Event=" << *i_evt << std::endl;
	  TGraphAsymmErrors * gr_all = new TGraphAsymmErrors();//all
	  TGraphAsymmErrors * gr_cut = new TGraphAsymmErrors();//counter cut
	  TGraphAsymmErrors * gr_sel = new TGraphAsymmErrors();//selected
	  Int_t n_pts_all = 0;
	  Int_t n_pts_cut = 0;
	  Int_t n_pts_sel = 0;
	  std::vector<std::pair<double,double> > data;
	  int thisc1=-1,thisc2=-1;
	  double thisc1z=-1,thisc2z=-1;
	  double thisc1x=-1,thisc2x=-1;
	  double thisconstant=-1,thisconstanterr=-1;
	  double thislinear=-1,thislinearerr=-1;
	  double thisquadratic=-1,thisquadraticerr=-1;
	  double thischi2=-1,thisndf=-1,thissumsqrresidual=-1;
	  double thissuccess=-1;
	  double thiserror=1;

	  double minx=99999,maxx=-99999;
	  double miny=99999,maxy=-99999;
	  
	  bool thisNS = false;
	  
	  reader.SetEntry(0);
	  while (reader.Next())
	    {
	            if (*c1==5 || *c2==5 || *c1==43 || *c2==43) continue;
	      if (*event == *i_evt && *run == (i_run->first))
		{
		  if ((*c1 >= 16 && *c1 <= 27) || (*c2 >= 16 && *c2 <= 27)) thisNS =true;
		  else thisNS = false;

		  //if (*ransac_success == 0) continue;
		  
		      if (*hitx < -400) continue;

		      if (thisNS)
			{
			  if (minx>*hitx) minx=*hitx;
			  if (maxx<*hitx) maxx=*hitx;
			  if (miny>*hitz) miny=*hitz;
			  if (maxy<*hitz) maxy=*hitz;
			}
		      else
			{
			  if (minx>*hitz) minx=*hitz;
			  if (maxx<*hitz) maxx=*hitz;
			  if (miny>*hitx) miny=*hitx;
			  if (maxy<*hitx) maxy=*hitx;
			}
		      
		      if (*fitrealhit)
			{
			  if (thisNS)
			    {
			      gr_sel->SetPoint(n_pts_sel,*hitx,*hitz);
			      gr_sel->SetPointError(n_pts_sel,*hiterrxlo,*hiterrxhi,*hiterrzlo,*hiterrzhi);
			    }
			  else
			    {
			      gr_sel->SetPoint(n_pts_sel,*hitz,*hitx);
			      gr_sel->SetPointError(n_pts_sel,*hiterrzlo,*hiterrzhi,*hiterrxlo,*hiterrxhi);
			    }
			  n_pts_sel++;
			}
		      else if (!(*countercut))
			{
		          if (thisNS)
			    {
			      gr_cut->SetPoint(n_pts_cut,*hitx,*hitz);
			      gr_cut->SetPointError(n_pts_cut,*hiterrxlo,*hiterrxhi,*hiterrzlo,*hiterrzhi);
			    }
			  else
			    {
			      gr_cut->SetPoint(n_pts_cut,*hitz,*hitx);
			      gr_cut->SetPointError(n_pts_cut,*hiterrzlo,*hiterrzhi,*hiterrxlo,*hiterrxhi);
			    }
			  n_pts_cut++;
		        }
		      else
			{
		          if (thisNS)
			    {
			      gr_all->SetPoint(n_pts_all,*hitx,*hitz);
			      gr_all->SetPointError(n_pts_all,*hiterrxlo,*hiterrxhi,*hiterrzlo,*hiterrzhi);
			    }
			  else
			    {
			      gr_all->SetPoint(n_pts_all,*hitz,*hitx);
			      gr_all->SetPointError(n_pts_all,*hiterrzlo,*hiterrzhi,*hiterrxlo,*hiterrxhi);
			    }
			  n_pts_all++;
		    }
		  thisc1 = *c1; thisc2 = *c2;
		  if (thisNS)
		    {
		      thisc1z = *c1z; thisc2z = *c2z;
		      thisc1x = *c1x; thisc2x = *c2x;
		    }
		  else
		    {
		      thisc1z = *c1x; thisc2z = *c2x;
		      thisc1x = *c1z; thisc2x = *c2z;
		    }
		  thisconstant = *fitconstant;
		  thisconstanterr = *fitconstanterr;
		  thislinear = *fitlinear;
		  thislinearerr = *fitlinearerr;
		  thisquadratic = *fitquadratic;
		  thisquadraticerr = *fitquadraticerr;
		  thischi2 = *fitchi2;
		  thissumsqrresidual = *fitsumsqrresidual;
		  thisndf = *fitndf;
		  thissuccess = *fitsuccess;
		  thiserror = *fitmle;
		}
	    }
		      
	  std::cout << "RANSAC result = " << thissuccess << std::endl;
	  std::cout << "     Chi2/NDF = " << thischi2/thisndf << std::endl;
	  std::cout << "     SSR/NDF  = " << thissumsqrresidual/thisndf << std::endl;
	  std::cout << "     Error    = " << thiserror << std::endl;
	  std::cout << "     Success Code = " << thissuccess << std::endl;
	  std::cout << " const=" << thisconstant << " (err=" << thisconstanterr << ")" << std::endl;
	  std::cout << " linear=" << thislinear << " (err=" << thislinearerr << ")" << std::endl;
	  std::cout << " quad=" << thisquadratic << " (err=" << thisquadraticerr << ")" << std::endl;

	  if (gr_all->GetN() + gr_cut->GetN() + gr_sel->GetN() == 0) continue;
	  
	  //if (thissumsqrresidual/thisndf > 2 || thischi2/thisndf > 2) continue;
	  //if (thissuccess != 1 || thissuccess != -1) continue;
	  
	  TF1 * model = new TF1("model","pol2",-50,250);
	  model->SetNpx(10000);
	  model->SetParameters(thisconstant,thislinear,thisquadratic);
	  
	  gr_all->SetMarkerStyle(8);
	  gr_all->SetMarkerSize(1);
	  gr_all->SetMarkerColor(4);
	  std::cout << "gr_all->N() = " << gr_all->GetN() << std::endl;
	  
	  gr_cut->SetMarkerStyle(8);
	  gr_cut->SetMarkerSize(1);
	  gr_cut->SetMarkerColor(3);
	  std::cout << "gr_cut->N() = " << gr_cut->GetN() << std::endl;

	  gr_sel->SetMarkerStyle(8);
	  gr_sel->SetMarkerSize(2);
	  gr_sel->SetMarkerColor(2);
	  std::cout << "gr_sel->N() = " << gr_sel->GetN() << std::endl;

	  TLine * l = new TLine(thisc1x,thisc1z,thisc2x,thisc2z);

	  canv1->cd();
	  if (gr_all->GetN() > 0)
	    {
	      gr_all->Draw("EAP");
	      gr_all->GetXaxis()->SetLimits(minx-5,maxx+5);
	      gr_all->GetHistogram()->SetMinimum(miny-5);
	      gr_all->GetHistogram()->SetMaximum(maxy+5);
	    }
	  if (gr_cut->GetN() > 0) gr_cut->Draw("PE SAME");
	  l->Draw("SAME");
	  if (gr_sel->GetN() > 0) gr_sel->Draw("PE SAME");
	  model->Draw("L SAME");
	  canv1->Update();
	  canv1->WaitPrimitive();
	}
    }
}
