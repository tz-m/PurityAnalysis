{
  Double_t  x[22] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21};
  Double_t  ex[22] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t  y[22] = {2518,2449,2341,2316,2286,2240,2154,2120,2080,2068,2027,2009,1950,1932,1871,1884,1817,1800,1835,1732,1739,1751};
  Double_t  ey[22] = {10.7,8.2,6.1,8.0,8.8,9.1,6.0,7.5,10.6,11.9,12.9,14.6,10.4,13.0,8.2,12.9,14.5,8.6,13.9,12.1,12.8,127.4};

  for (int i = 0; i < 22; ++i)
    {
      x[i] = x[i]*91.45+45.725;
    }

  TGraphErrors * gr = new TGraphErrors(22,x,y,ex,ey);
  
  TF1 * expo = new TF1( "expo", "[0]*exp( -x/[1] )", 100, 900 );

  expo->SetParNames( "dQdx0", "eLifetime" );
  expo->SetParameters( 3000, 3000 );

  gStyle->SetOptStat( 0 );
  gStyle->SetOptFit( 1 );
  
  TCanvas * canv = new TCanvas("canv","",1600,900);
  gr->Draw("APE");
  gr->SetMarkerStyle(kFullDotLarge);
  gr->SetMarkerSize(1);
  gr->Fit("expo","R");
  canv->SaveAs("michelleMPV.png");
  
		       
}
