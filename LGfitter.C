#include "string.h"
#include "TMinuit.h"
#include "TROOT.h"


Double_t LandFun(Double_t *x, Double_t *par) {

  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density (sigma) 
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Landau Amplitude
  //
  //In the Landau distribution (represented by the CERNLIB approximation), 
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter
  
  // Numeric constants
  Double_t mpshift  = -0.22278298;       // Landau maximum location


  Double_t mpc;
  mpc=par[1]-mpshift*par[0];

  Double_t temp;
  temp = par[2]*(TMath::Landau(x[0],par[1],par[0]));

  return(temp);
}


Double_t langaufun(Double_t *x, Double_t *par) {

  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation), 
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location

  // Control constants - need to be tuned for each specific application
  Double_t np = 10000.0;      // number of convolution steps
  Double_t sc =   8.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;


  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp-xlow) / np;

  if (par[0]==0) sum=0;
  else {
  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}



TF1 *LGfitter(TH1 *h,  Double_t lbound, Double_t rbound, Double_t *fitparams, Double_t *fiterrors, Double_t *covmat, Int_t ief)
{

  // Fit the histogram h to a Landau function and pass back the parameters
  //    and their errors.  
  //  Use TMath::Landau, but correct for MPV offset
  //  Initializes the fit parameters to reasonable values
  //  Pass back the fit parameters and their errors
  //  Return the fit function

  //  Note: to get this to work correctly, you need to tune the 
  //    two constants "control parameters" in langaufun
  //    to match your application (sc to match the gaus width and 
  //    np to accomodate the histogram binning)!!

  //gStyle->SetOptFit(12);

  //  Fit histogram to Landau/Gaussian conv
  Char_t FunName[100]; sprintf(FunName,"Fitfcn_%s",h->GetName());
  TF1 *ffit = new TF1(FunName,langaufun,lbound,rbound,4);
  cout << "LW error" << fiterrors[0] << endl;
  ffit->SetParameters(fitparams[0],fitparams[1],fitparams[2],fitparams[3]);
  ffit->SetParError(0,fiterrors[0]);
  if (fiterrors[0]==0) ffit->FixParameter(0,fitparams[0]);
  ffit->SetParError(1,fiterrors[1]);
  ffit->SetParError(2,fiterrors[2]);
  ffit->SetParError(3,fiterrors[3]);
  ffit->SetParNames("Width","MP","Amp","Sigma");  
  // If the bins are large w.r.t. to the rising slope, you may 
  //     need to use the I option when fitting.  i.e. "IMLEVR"
  TFitResultPtr r;
  // removed M,E options as a test
  if (ief)   r = h->Fit(FunName,"LVRE");
  else  r = h->Fit(FunName,"LQRE");

  // Check fit status 
  TString test =  gMinuit->fCstatu.Data();
  Double_t a = ffit->GetParameter(0);
  fiterrors[0] = -1000.0;
  fiterrors[1] = -1000.0;
  fiterrors[2] = -1000.0;
  fiterrors[3] = -1000.0;
  fitparams[0] = -1000.0;
  fitparams[1] = -1000.0;
  fitparams[2] = -1000.0;
  fitparams[3] = -1000.0;
  covmat[0] = -1000.0;
  covmat[1] = -1000.0;
  covmat[2] = -1000.0;  
  covmat[3] = -1000.0;  
  //  if (ii<20) return(ffit);
  cout << "here  " << test << endl;
  if (test.BeginsWith("SUCC")) {   //successful fit 

  // Get Fit Parameters, their errors and cov matrix
  fitparams[0] = h->GetFunction(FunName)->GetParameter(0);
  fitparams[1] = h->GetFunction(FunName)->GetParameter(1);
  fitparams[2] = h->GetFunction(FunName)->GetParameter(2);
  fitparams[3] = h->GetFunction(FunName)->GetParameter(3);
  fiterrors[0] = h->GetFunction(FunName)->GetParError(0);
  fiterrors[1] = h->GetFunction(FunName)->GetParError(1);
  fiterrors[2] = h->GetFunction(FunName)->GetParError(2);
  fiterrors[3] = h->GetFunction(FunName)->GetParError(3);
  TVirtualFitter *fitter = TVirtualFitter::GetFitter();
  TMatrixD matrix(4,4, fitter->GetCovarianceMatrix());
  covmat[0] = fitter->GetCovarianceMatrixElement(0, 1);
  covmat[1] = fitter->GetCovarianceMatrixElement(0, 2);
  covmat[2] = fitter->GetCovarianceMatrixElement(1, 2);
  covmat[3] = fitter->GetCovarianceMatrixElement(0, 3);
  cout << "cov int " << covmat[3] << endl;
  //missing covariance terms here !
  }

  return(ffit);

}

