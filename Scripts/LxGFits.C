void LxGFits()
{
  TFile * file = TFile::Open("/home/mthiesse/PurityAnalysis/PurityAnalysis.root","READ");
  if (!file || !file->IsOpen()) return;

  TTreeReader reader("results",file);
  TTreeReaderValue<Float_t> landmpv(reader,"landmpv");
  TTreeReaderValue<Float_t> landwidth(reader,"landwidth");
  TTreeReaderValue<Float_t> gaussmean_l(reader,"gaussmean_l");
  TTreeReaderValue<Float_t> gausswidth_l(reader,"gausswidth_l");
  TTreeReaderValue<Float_t> landmpverr(reader,"landmpverr");
  TTreeReaderValue<Float_t> landwidtherr(reader,"landwidtherr");
  TTreeReaderValue<Float_t> gaussmean_lerr(reader,"gaussmean_lerr");
  TTreeReaderValue<Float_t> gausswidth_lerr(reader,"gausswidth_lerr");
  TTreeReaderValue<Float_t> bincenter(reader,"bincenter");
  TTreeReaderValue<Float_t> lxg_chi2ndf(reader,"lxg_chi2ndf");
  TTreeReaderValue<Float_t> lxglow(reader,"lxglow");
  TTreeReaderValue<Float_t> lxghigh(reader,"lxghigh");

  RooRealVar charge("charge","dQdx (ADC/cm)",200,5000);
  charge.setBins(10000,"fft");
  RooRealVar ml("LandMPV","mean landau",2800,-1000,200000);
  RooRealVar sl("LandWidth","sigma landau",1000,0,100000);
  RooLandau landau("lx","lx",charge,ml,sl);
  RooRealVar mg("mg","mean gauss",0);
  RooRealVar sg("GaussWidth","sigma gauss",50,0,20000);
  RooGaussian gauss("gauss","gauss",charge,mg,sg);
  RooFFTConvPdf lxg("lxg","landau (x) gauss",charge,landau,gauss);

  Int_t palette[22];
  Double_t Red[] = {0.,0.0,1.0,1.0,1.0};
  Double_t Green[] = {0.,0.0,0.0,1.0,1.0};
  Double_t Blue[]   = {0., 1.0, 0.0, 0.0, 1.0};
  Double_t Length[] = {0., .25, .50, .75, 1.0};
  Int_t FI = TColor::CreateGradientColorTable(5, Length, Red, Green, Blue, 22);
  for (int i=0; i<22; i++) palette[i] = FI+i;
  
  RooPlot * lxgframe = charge.frame(RooFit::Title("Landau(x)Gauss Best Fit"));
  TLegend * tlxg = new TLegend(0.8,0.3,0.95,0.95);

  Int_t col = 0;
  while (reader.Next())
    {
      ml.setVal(*landmpv);
      sl.setVal(*landwidth);
      mg.setVal(*gaussmean_l);
      sg.setVal(*gausswidth_l);
      lxg.plotOn(lxgframe,RooFit::Name(TString::Format("%.1f",*bincenter)),RooFit::LineColor(palette[col]),RooFit::Range(*lxglow,*lxghigh));
      tlxg->AddEntry(lxgframe->findObject(TString::Format("%.1f",*bincenter)),TString::Format("Bin @ %.1f",*bincenter),"L");
      col++;
    }

  TCanvas * canv = new TCanvas("lxg","lxg",2000,1200);
  lxgframe->GetYaxis()->SetTitleOffset(1.4);
  lxgframe->Draw();
  tlxg->Draw();
  canv->Update();
  canv->SaveAs("LxGFits.png");
  
}
