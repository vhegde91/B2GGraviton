void PDF_ScaleuncPlotter(TString iFname, TString catName){
  TFile *f=new TFile(iFname);
  if(iFname.Contains("2016")) intLumi=35.8;
  else if(iFname.Contains("2017")) intLumi=41.5;
  else if(iFname.Contains("2018")) intLumi=59.5;
 
  TLegend *l = new TLegend(0.55,0.78,0.85,0.88); l->SetLineColor(0);
  TString name1,name2,name1_2d,name2_2d,name1_3d,name2_3d,newName,legName;
  if(catName=="VBF HP"){ name1="MTvBinvbf_0";          name2="MTvBinvbf_4"; legName="#alpha_{VBF HP}";}
  if(catName=="VBF LP"){ name1="MTvBinvbf_1";          name2="MTvBinvbf_5"; legName="#alpha_{VBF LP}";}
  if(catName=="ggF HP"){ name1="MTvBinggF_2";          name2="MTvBinggF_6"; legName="#alpha_{ggF HP}";}
  if(catName=="ggF LP"){ name1="MTvBinggF_3";          name2="MTvBinggF_7"; legName="#alpha_{ggF LP}";}

  TH1D *h = (TH1D*)f->Get(name1);
  TH1D *hd= (TH1D*)f->Get(name2);
  h->Divide(hd);
  l->AddEntry(h,"Statistical (2018 MC)","f");
  float err=0;
  for(int i=0;i<=h->GetNbinsX()+1;i++){
    err = h->GetBinError(i)/h->GetBinContent(i);
    h->SetBinContent(i,1);
    h->SetBinError(i,err);
  }
  //---------- for scale syst
  // name1_2d.push_back("MTvBinvbf_Sclidx0"); name2_2d.push_back("MTvBinvbf_Sclidx4"); rebin.push_back(1); legName.push_back("#alpha_{VBF HP}");
  // name1_2d.push_back("MTvBinvbf_Sclidx1"); name2_2d.push_back("MTvBinvbf_Sclidx5"); rebin.push_back(1); legName.push_back("#alpha_{VBF LP}");
  // name1_2d.push_back("MTvBinggF_Sclidx2"); name2_2d.push_back("MTvBinggF_Sclidx6"); rebin.push_back(1); legName.push_back("#alpha_{ggF HP}");
  // name1_2d.push_back("MTvBinggF_Sclidx3"); name2_2d.push_back("MTvBinggF_Sclidx7"); rebin.push_back(1); legName.push_back("#alpha_{ggF LP}");

  //---------- for PDF syst
  // name1_2d.push_back("MTvBinvbf_PDFidx0"); name2_2d.push_back("MTvBinvbf_PDFidx4"); rebin.push_back(1); legName.push_back("#alpha_{VBF HP}");
  // name1_2d.push_back("MTvBinvbf_PDFidx1"); name2_2d.push_back("MTvBinvbf_PDFidx5"); rebin.push_back(1); legName.push_back("#alpha_{VBF LP}");
  // name1_2d.push_back("MTvBinggF_PDFidx2"); name2_2d.push_back("MTvBinggF_PDFidx6"); rebin.push_back(1); legName.push_back("#alpha_{ggF HP}");
  // name1_2d.push_back("MTvBinggF_PDFidx3"); name2_2d.push_back("MTvBinggF_PDFidx7"); rebin.push_back(1); legName.push_back("#alpha_{ggF LP}");

  double x1,y1, x2,y2, x3,y3, x4,y4;
  x1 = h->GetBinLowEdge(1);
  x4 = h->GetBinLowEdge(h->GetNbinsX()+1);

  if(catName=="VBF HP"){
    y1 = 0.01;
    x2 = 700.1;
    y2 = 0.01;
    x3 = 700.1;
    y3 = 0.02;
    y4 = 0.02;
  }
  if(catName=="VBF LP"){
    y1 = 0.02;
    x2 = 1200.1;
    y2 = 0.005;
    x3 = 1200.1;
    y3 = 0.005;
    y4 = 0.02;
  }

  if(catName=="ggF HP"){
    y1 = 0.02;
    x2 = 700;
    y2 = 0.01;
    x3 = 700.1;
    y3 = 0.02;
    x4 = 3200;
    y4 = 0.02;
  }
  if(catName=="ggF LP"){
    y1 = 0.02;
    x2 = 1300.;
    y2 = 0.005;
    x3 = 1300.1;
    y3 = 0.005;
    x4 = 3200;
    y4 = 0.02;
  }


  vector<double> x,y,ex,eyl,eyh;
  double a=0.05;
  for(int i=1;i<=h->GetNbinsX()+1;i++){
    x.push_back(h->GetBinLowEdge(i));
    y.push_back(1.0);
    if(h->GetBinLowEdge(i) < x2){
      eyl.push_back(0);
      a = y1 + ((y2-y1)/(x2-x1))*(h->GetBinLowEdge(i)-x1);
      eyh.push_back(a);
    }
    else{
      eyh.push_back(0);
      a = y3 + ((y4-y3)/(x4-x3))*(h->GetBinLowEdge(i)-x3);
      eyl.push_back(a);
    }
    ex.push_back(0);
  }

  TGraphAsymmErrors *g=new TGraphAsymmErrors(x.size(),&(x[0]),&(y[0]),&(ex[0]),&(ex[0]),&(eyl[0]),&(eyh[0]));
  g->SetFillColor(kYellow);
  g->SetLineColor(0);
  g->GetXaxis()->SetRangeUser(h->GetBinLowEdge(1),h->GetBinLowEdge(h->GetNbinsX()));
  g->GetYaxis()->SetRangeUser(0.9,1.1);
  g->Draw("a3");
  g->SetTitle(";mT (GeV);Normalized uncertainty");
  l->AddEntry(g,"#mu_{R}, #mu_{F} scale (2016/2017/2018)","f");
  
  //h->SetFillStyle(3001);
  h->SetLineColor(0);
  h->SetFillColor(kGray);
  h->Draw("E2 same");
  l->Draw();

  TPaveText *tx = new TPaveText(0.35,0.8,0.45,0.87,"NDC");
  tx->SetFillColor(0);
  tx->SetShadowColor(0);
  tx->AddText(catName);
  tx->Draw();
  //gPad->SetGridx(0);
  //gPad->SetGridy(0);

}
