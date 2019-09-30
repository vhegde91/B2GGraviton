void systUncHistMaker(){
  TFile *fout = new TFile("SystUnc.root","recreate");
  TH1D *h_pdf[4], *h_scale[4];
  vector<TString> name = {"VBF_HP","VBF_LP","ggF_HP","ggF_LP"};
  
  vector<double> ggfEdges = {500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2250,2400,2550,2700,2900,3200};
  vector<double> vbfEdges = {500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1850,2100,2350};

  h_pdf[0] = new TH1D("PDFunc_"+name[0],"PDF unc for 137/fb for "+name[0],vbfEdges.size()-1,&(vbfEdges)[0]);
  h_pdf[1] = new TH1D("PDFunc_"+name[1],"PDF unc for 137/fb for "+name[1],vbfEdges.size()-1,&(vbfEdges)[0]);

  h_pdf[2] = new TH1D("PDFunc_"+name[2],"PDF unc for 137/fb for "+name[2],ggfEdges.size()-1,&(ggfEdges)[0]);
  h_pdf[3] = new TH1D("PDFunc_"+name[3],"PDF unc for 137/fb for "+name[3],ggfEdges.size()-1,&(ggfEdges)[0]);

  h_scale[0] = new TH1D("Scaleunc_"+name[0],"Scale unc for 137/fb for "+name[0],vbfEdges.size()-1,&(vbfEdges)[0]);
  h_scale[1] = new TH1D("Scaleunc_"+name[1],"Scale unc for 137/fb for "+name[1],vbfEdges.size()-1,&(vbfEdges)[0]);

  h_scale[2] = new TH1D("Scaleunc_"+name[2],"Scale unc for 137/fb for "+name[2],ggfEdges.size()-1,&(ggfEdges)[0]);
  h_scale[3] = new TH1D("Scaleunc_"+name[3],"Scale unc for 137/fb for "+name[3],ggfEdges.size()-1,&(ggfEdges)[0]);

  float pdfUncValues[4] = {1.027, 1.027, 1.013, 1.013};
  for(int i=0;i<name.size();i++){
    for(int j=1;j<=h_pdf[i]->GetNbinsX();j++){
      h_pdf[i]->SetBinContent(j,pdfUncValues[i]);}}
  for(int i=0;i<name.size();i++)
    h_pdf[i]->Write();
  //-------------------------
  double x1,y1, x2,y2, x3,y3, x4,y4;
  for(int i=0;i<name.size();i++){
    x1 = h_scale[i]->GetBinLowEdge(1);
    x4 = h_scale[i]->GetBinLowEdge(h_scale[i]->GetNbinsX()+1);
    if(name[i]=="VBF_HP"){
      y1 = 0.01;
      x2 = 700.1;
      y2 = 0.01;
      x3 = 700.1;
      y3 = 0.02;
      y4 = 0.02;
    }
    if(name[i]=="VBF_LP"){
      y1 = 0.02;
      x2 = 1200.1;
      y2 = 0.005;
      x3 = 1200.1;
      y3 = 0.005;
      y4 = 0.02;
    }
    if(name[i]=="ggF_HP"){
      y1 = 0.02;
      x2 = 700;
      y2 = 0.01;
      x3 = 700.1;
      y3 = 0.02;
      x4 = 3200;
      y4 = 0.02;
    }
    if(name[i]=="ggF_LP"){
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
    for(int j=1;j<=h_scale[i]->GetNbinsX()+1;j++){
      //      x.push_back(h_scale[i]->GetBinLowEdge(j));
      y.push_back(1.0);
      if(h_scale[i]->GetBinLowEdge(j) < x2){
	eyl.push_back(0);
	a = y1 + ((y2-y1)/(x2-x1))*(h_scale[i]->GetBinLowEdge(j)-x1);
	eyh.push_back(a);
	h_scale[i]->SetBinContent(j,1+a);
      }
      else{
	eyh.push_back(0);
	a = y3 + ((y4-y3)/(x4-x3))*(h_scale[i]->GetBinLowEdge(j)-x3);
	eyl.push_back(a);
	h_scale[i]->SetBinContent(j,1-a);
      }
      //      ex.push_back(0);
    }
  }
  //-------------------------
  //  fout->cd();
  for(int i=0;i<name.size();i++)
    h_scale[i]->Write();
  
}
