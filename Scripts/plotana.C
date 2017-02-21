{
  Double_t noise[4] = {0.5,1.0,1.5,2.0};
  Double_t noiseerr[4] = {0,0,0,0};
  Double_t noise3[3] = {0.5,1.0,1.5};
  Double_t noiseerr3[3] = {0,0,0};
  
  Double_t p1[4] = {0.5965,0.8209,0.9055,0.9386};
  Double_t pe1[4] = {0.4198,0.3219,0.2351,0.1842};
  Double_t p2[4] = {0.7809,0.9414,0.9646,0.9708};
  Double_t pe2[4] = {0.3405,0.1781,0.1275,0.1069};
  Double_t p3[4] = {0.8474,0.9567,0.971,0.9728};
  Double_t pe3[4] = {0.2912,0.1478,0.1061,0.1016};
  Double_t p5[3] = {0.8886,0.9628,0.9723};
  Double_t pe5[3] = {0.252,0.132,0.1029};
  Double_t p8[4] = {0.909,0.9662,0.9727,0.9742};
  Double_t pe8[4] = {0.2274,0.1231,0.1022,0.09806};

  Double_t e1[4] = {0.05339,0.1754,0.3017,0.3922};
  Double_t ee1[4] = {0.07873,0.1831,0.2476,0.2723};
  Double_t e2[4] = {0.09195,0.2986,0.4761,0.5746};
  Double_t ee2[4] = {0.08871,0.1873,0.2348,0.2502};
  Double_t e3[4] = {0.1168,0.3643,0.5455,0.6276};
  Double_t ee3[4] = {0.09185,0.1879,0.2348,0.2535};
  Double_t e5[3] = {0.1436,0.4259,0.5951};
  Double_t ee5[3] = {0.09535,0.1915,0.24};
  Double_t e8[4] = {0.1627,0.46,0.617,0.6721};
  Double_t ee8[4] = {0.09884,0.1983,0.2473,0.2657};
  
  Double_t mcc1[4] = {0.04336,0.2368,0.3764,0.4681};
  Double_t mcce1[4] = {0.2179,0.2502,0.2571,0.2506};
  Double_t mcc2[4] = {0.1344,0.3856,0.552,0.643};
  Double_t mcce2[4] = {0.1946,0.1877,0.1775,0.1638};
  Double_t mcc3[4] = {0.1799,0.4504,0.6164,0.6926};
  Double_t mcce3[4] = {0.18,0.1648,0.1511,0.1493};
  Double_t mcc5[3] = {0.2212,0.5059,0.6604};
  Double_t mcce5[3] = {0.1673,0.1472,0.1432};
  Double_t mcc8[4] = {0.2479,0.5373,0.6817,0.736};
  Double_t mcce8[4] = {0.1596,0.1409,0.1437,0.1514};
  
  TGraphErrors * purity1 = new TGraphErrors(4,noise,p1,noiseerr,pe1);
  purity1->SetMarkerStyle(21);
  purity1->SetMarkerColor(4);
  purity1->SetLineColor(4);
  purity1->SetLineWidth(3);
  purity1->SetMarkerSize(2);
  TGraphErrors * purity2 = new TGraphErrors(4,noise,p2,noiseerr,pe2);
  purity2->SetMarkerStyle(21);
  purity2->SetMarkerColor(5);
  purity2->SetLineColor(5);
  purity2->SetLineWidth(3);
  purity2->SetMarkerSize(2);
  TGraphErrors * purity3 = new TGraphErrors(4,noise,p3,noiseerr,pe3);
  purity3->SetMarkerStyle(21);
  purity3->SetMarkerColor(6);
  purity3->SetLineColor(6);
  purity3->SetLineWidth(3);
  purity3->SetMarkerSize(2);
  TGraphErrors * purity5 = new TGraphErrors(3,noise3,p5,noiseerr3,pe5);
  purity5->SetMarkerStyle(21);
  purity5->SetMarkerColor(7);
  purity5->SetLineColor(7);
  purity5->SetLineWidth(3);
  purity5->SetMarkerSize(2);
  TGraphErrors * purity8 = new TGraphErrors(4,noise,p8,noiseerr,pe8);
  purity8->SetMarkerStyle(21);
  purity8->SetMarkerColor(3);
  purity8->SetLineColor(3);
  purity8->SetLineWidth(3);
  purity8->SetMarkerSize(2);

  TCanvas * canvpurity = new TCanvas("canvpurity","Purity",2000,1600);
  purity1->Draw("ALP");
  purity1->SetTitle("Reconstruction Purity");
  purity1->GetXaxis()->SetTitle("Signal-To-Noise Scale Factor");
  purity1->GetYaxis()->SetTitle("TP/(TP+FP)");
  purity2->Draw("sameLP");
  purity3->Draw("sameLP");
  purity5->Draw("sameLP");
  purity8->Draw("sameLP");

  TLegend * legpurity = new TLegend(0.7,0.1,0.9,0.3);
  legpurity->AddEntry(purity1,"1ms lifetime","lep");
  legpurity->AddEntry(purity2,"2ms lifetime","lep");
  legpurity->AddEntry(purity3,"3ms lifetime","lep");
  legpurity->AddEntry(purity5,"5ms lifetime","lep");
  legpurity->AddEntry(purity8,"8ms lifetime","lep");
  legpurity->Draw();

  /////////////////////////////////////////////////////////////////////////

  TGraphErrors * efficiency1 = new TGraphErrors(4,noise,e1,noiseerr,ee1);
  efficiency1->SetMarkerStyle(21);
  efficiency1->SetMarkerColor(4);
  efficiency1->SetLineColor(4);
  efficiency1->SetLineWidth(3);
  efficiency1->SetMarkerSize(2);
  TGraphErrors * efficiency2 = new TGraphErrors(4,noise,e2,noiseerr,ee2);
  efficiency2->SetMarkerStyle(21);
  efficiency2->SetMarkerColor(5);
  efficiency2->SetLineColor(5);
  efficiency2->SetLineWidth(3);
  efficiency2->SetMarkerSize(2);
  TGraphErrors * efficiency3 = new TGraphErrors(4,noise,e3,noiseerr,ee3);
  efficiency3->SetMarkerStyle(21);
  efficiency3->SetMarkerColor(6);
  efficiency3->SetLineColor(6);
  efficiency3->SetLineWidth(3);
  efficiency3->SetMarkerSize(2);
  TGraphErrors * efficiency5 = new TGraphErrors(3,noise3,e5,noiseerr3,ee5);
  efficiency5->SetMarkerStyle(21);
  efficiency5->SetMarkerColor(7);
  efficiency5->SetLineColor(7);
  efficiency5->SetLineWidth(3);
  efficiency5->SetMarkerSize(2);
  TGraphErrors * efficiency8 = new TGraphErrors(4,noise,e8,noiseerr,ee8);
  efficiency8->SetMarkerStyle(21);
  efficiency8->SetMarkerColor(3);
  efficiency8->SetLineColor(3);
  efficiency8->SetLineWidth(3);
  efficiency8->SetMarkerSize(2);

  TCanvas * canvefficiency = new TCanvas("canvefficiency","Efficiency",2000,1600);
  efficiency1->Draw("ALP");
  efficiency1->SetTitle("Reconstruction Efficiency");
  efficiency1->GetXaxis()->SetTitle("Signal-To-Noise Scale Factor");
  efficiency1->GetYaxis()->SetTitle("TP/(TP+FN)");
  efficiency2->Draw("sameLP");
  efficiency3->Draw("sameLP");
  efficiency5->Draw("sameLP");
  efficiency8->Draw("sameLP");

  TLegend * legefficiency = new TLegend(0.7,0.1,0.9,0.3);
  legefficiency->AddEntry(efficiency1,"1ms lifetime","lep");
  legefficiency->AddEntry(efficiency2,"2ms lifetime","lep");
  legefficiency->AddEntry(efficiency3,"3ms lifetime","lep");
  legefficiency->AddEntry(efficiency5,"5ms lifetime","lep");
  legefficiency->AddEntry(efficiency8,"8ms lifetime","lep");
  legefficiency->Draw();

  //////////////////////////////////////////////////////////////////////////////////

  TGraphErrors * mccoef1 = new TGraphErrors(4,noise,mcc1,noiseerr,mcce1);
  mccoef1->SetMarkerStyle(21);
  mccoef1->SetMarkerColor(4);
  mccoef1->SetLineColor(4);
  mccoef1->SetLineWidth(3);
  mccoef1->SetMarkerSize(2);
  TGraphErrors * mccoef2 = new TGraphErrors(4,noise,mcc2,noiseerr,mcce2);
  mccoef2->SetMarkerStyle(21);
  mccoef2->SetMarkerColor(5);
  mccoef2->SetLineColor(5);
  mccoef2->SetLineWidth(3);
  mccoef2->SetMarkerSize(2);
  TGraphErrors * mccoef3 = new TGraphErrors(4,noise,mcc3,noiseerr,mcce3);
  mccoef3->SetMarkerStyle(21);
  mccoef3->SetMarkerColor(6);
  mccoef3->SetLineColor(6);
  mccoef3->SetLineWidth(3);
  mccoef3->SetMarkerSize(2);
  TGraphErrors * mccoef5 = new TGraphErrors(3,noise3,mcc5,noiseerr3,mcce5);
  mccoef5->SetMarkerStyle(21);
  mccoef5->SetMarkerColor(7);
  mccoef5->SetLineColor(7);
  mccoef5->SetLineWidth(3);
  mccoef5->SetMarkerSize(2);
  TGraphErrors * mccoef8 = new TGraphErrors(4,noise,mcc8,noiseerr,mcce8);
  mccoef8->SetMarkerStyle(21);
  mccoef8->SetMarkerColor(3);
  mccoef8->SetLineColor(3);
  mccoef8->SetLineWidth(3);
  mccoef8->SetMarkerSize(2);

  TCanvas * canvmccoef = new TCanvas("canvmccoef","MCC",2000,1600);
  mccoef1->Draw("ALP");
  mccoef1->SetTitle("Matthews Correlation Coefficient");
  mccoef1->GetXaxis()->SetTitle("Signal-To-Noise Scale Factor");
  mccoef1->GetYaxis()->SetTitle("(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))");
  mccoef2->Draw("sameLP");
  mccoef3->Draw("sameLP");
  mccoef5->Draw("sameLP");
  mccoef8->Draw("sameLP");

  TLegend * legmccoef = new TLegend(0.7,0.1,0.9,0.3);
  legmccoef->AddEntry(mccoef1,"1ms lifetime","lep");
  legmccoef->AddEntry(mccoef2,"2ms lifetime","lep");
  legmccoef->AddEntry(mccoef3,"3ms lifetime","lep");
  legmccoef->AddEntry(mccoef5,"5ms lifetime","lep");
  legmccoef->AddEntry(mccoef8,"8ms lifetime","lep");
  legmccoef->Draw();
  
}
 
