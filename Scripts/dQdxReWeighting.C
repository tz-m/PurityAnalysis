void dQdxReWeighting()
{
  std::vector<Float_t> qreco;
  std::vector<Float_t> qmc;

  TProfile * recoprof = new TProfile("recoprof",";Drift Distance;dQ/dx",22,0,2012,-1000,20000,"s");
  recoprof->SetLineColor(4);
  recoprof->SetLineWidth(2);
  recoprof->SetStats(false);
  TProfile * mcprof = new TProfile("mcprof",";Drift Distance;dQ/dx",22,0,2012,-1000,20000,"s");
  mcprof->SetStats(false);
  mcprof->SetLineColor(3);
  mcprof->SetLineWidth(2);

  TFile * recofile = TFile::Open("/home/mthiesse/PurityAnalysis/PurityAnalysis_data.root","READ");

  //TTreeReader recohitsreader("hits",recofile);
  //TTreeReaderValue<Float_t> recohitt(recohitsreader,"hitt");
  //TTreeReaderValue<Float_t> recodqdx(recohitsreader,"dqdx");

  TTreeReader recoresultreader("results",recofile);
  TTreeReaderValue<Float_t> recocenter(recoresultreader,"bincenter");
  TTreeReaderValue<Float_t> recompv(recoresultreader,"landmpv");

  while (recoresultreader.Next())
    {
      recoprof->Fill(*recocenter,*recompv);
    }
  recoprof->SetMinimum(1700);

  TFile * mcfile = TFile::Open("/home/mthiesse/PurityAnalysis/PurityAnalysis_sim_3ms_1.0.root","READ");

  //TTreeReader mcreader("hits",mcfile);
  //TTreeReaderValue<Float_t> mchitt(mcreader,"hitt");
  //TTreeReaderValue<Float_t> mcdqdx(mcreader,"dqdx");

  TTreeReader mcresultreader("results",mcfile);
  TTreeReaderValue<Float_t> mccenter(mcresultreader,"bincenter");
  TTreeReaderValue<Float_t> mcmpv(mcresultreader,"landmpv");

  while (mcresultreader.Next())
    {
      mcprof->Fill(*mccenter,*mcmpv);
    }

  TProfile * copyrecoprof = (TProfile*)recoprof->Clone();
  copyrecoprof->SetStats(false);
  copyrecoprof->Divide(mcprof);
  copyrecoprof->SetMinimum(-0.25);
  copyrecoprof->SetMaximum(2.5);
  copyrecoprof->GetYaxis()->SetTitle("(MPV dQ_{Data}/dx) / (MPV dQ_{MC}/dx)");
  copyrecoprof->SetMarkerColor(1);
  copyrecoprof->SetLineWidth(1);
  copyrecoprof->SetMarkerStyle(20);

  TCanvas * c = new TCanvas("c","",2000,1600);
  c->Divide(1,2);
  c->cd(1);
  recoprof->Draw();
  mcprof->Draw("same");
  c->cd(2);
  copyrecoprof->Draw();
  TLine * l = new TLine(0,1,2012,1);
  l->SetLineStyle(kDashed);
  l->Draw();
  c->Update();



}
