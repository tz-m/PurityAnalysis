void MigrationMatrix()
{
  TH2D * mm = new TH2D("mm",";chg(reco);chg(MC)",400,200,8000,200,0,10000);

  TFile * file = TFile::Open("/home/mthiesse/PurityAnalysis/Scripts/ana_mixer_hist.root","READ");
  TTreeReader reader("robustmcana/mcanahits",file);
  TTreeReaderValue<Bool_t> foundBoth(reader,"foundBoth");
  TTreeReaderValue<Double_t> RecoQ(reader,"RecoQ");
  TTreeReaderValue<Double_t> MCQ(reader,"MCQ");
  while (reader.Next())
    {
      if (*foundBoth && *RecoQ > 200 && *RecoQ < 8000 && *MCQ > 0 && *MCQ < 10000)
	mm->Fill(*RecoQ,*MCQ);
      else
	mm->Fill(-1,*MCQ);
    }
  gStyle->SetOptStat(0);
  TCanvas * canv = new TCanvas("canv","",2000,1600);
  mm->Draw("colz");
  canv->SaveAs("migrationmatrix.png");
  mm->SaveAs("MigrationMatrix.root");
}
