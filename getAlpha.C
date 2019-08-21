#include<iostream>
#include<iomanip>
#include"TH1.h"
#include"TROOT.h"
#include"TH2.h"
#include"TFile.h"
#include"TDirectory.h"
#include"TTree.h"
#include"TBrowser.h"
#include"TF1.h"
#include<string>
#include<vector>
#include"TGraphErrors.h"
#include"TGraph.h"
#include"TLegend.h"
#include"TLatex.h"
#include"TCanvas.h"
#include"THStack.h"
#include"TStyle.h"

//int col[10]={kOrange,kBlue,kTeal+9,kGray+1,kCyan,kOrange-9,kYellow+2,kRed,kMagenta+2,kMagenta};  //Specify Colors
int col[11]={kBlack,kBlue,kRed,kTeal+9,kBlack,kCyan,kOrange,kGray+1,kMagenta+2,kYellow+2,kMagenta};  //Specify Colors
TString name;
TString nameX;
double sr_Integral=0,cr_Integral=0,sr_IntrErr=0,cr_IntrErr=0,tf_err=0;
void inclOverflow(TH2D*);
TH2D rebin2DHist(TH2D*,int,double*,int,double*);
TString getCatName(TString);
TLatex textOnTop,intLumiE;
double intLumi=35.9;
bool savepdf=0;

void getAlpha(TString iFname){
  TH1::SetDefaultSumw2(1);
  gStyle->SetOptStat(0);
  
  TFile *f = new TFile(iFname);
  name="alpha_"+iFname;
  TFile *fout=new TFile(name,"RECREATE");
  if(iFname.Contains("2016")) intLumi=35.8;
  else if(iFname.Contains("2017")) intLumi=41.5;
  else if(iFname.Contains("2018")) intLumi=59.5;

  
  vector<TString> name1,name2,name1_2d,name2_2d,name1_3d,name2_3d,newName;
  vector<int> rebin;
  vector<TString> legName;
  name1.push_back("MTvBinvbf_0");          name2.push_back("MTvBinvbf_4");       rebin.push_back(1); legName.push_back("#alpha_{VBF HP}");
  name1.push_back("MTvBinvbf_1");          name2.push_back("MTvBinvbf_5");       rebin.push_back(1); legName.push_back("#alpha_{VBF LP}");
  name1.push_back("MTvBinggF_2");          name2.push_back("MTvBinggF_6");       rebin.push_back(1); legName.push_back("#alpha_{ggF HP}");
  name1.push_back("MTvBinggF_3");          name2.push_back("MTvBinggF_7");       rebin.push_back(1); legName.push_back("#alpha_{ggF LP}");

  //---------- for scale syst
  // name1_2d.push_back("MTvBinvbf_Sclidx0"); name2_2d.push_back("MTvBinvbf_Sclidx4"); rebin.push_back(1); legName.push_back("#alpha_{VBF HP}");
  // name1_2d.push_back("MTvBinvbf_Sclidx1"); name2_2d.push_back("MTvBinvbf_Sclidx5"); rebin.push_back(1); legName.push_back("#alpha_{VBF LP}");
  // name1_2d.push_back("MTvBinggF_Sclidx2"); name2_2d.push_back("MTvBinggF_Sclidx6"); rebin.push_back(1); legName.push_back("#alpha_{ggF HP}");
  // name1_2d.push_back("MTvBinggF_Sclidx3"); name2_2d.push_back("MTvBinggF_Sclidx7"); rebin.push_back(1); legName.push_back("#alpha_{ggF LP}");

  //---------- for PDF syst
  name1_2d.push_back("MTvBinvbf_PDFidx0"); name2_2d.push_back("MTvBinvbf_PDFidx4"); rebin.push_back(1); legName.push_back("#alpha_{VBF HP}");
  name1_2d.push_back("MTvBinvbf_PDFidx1"); name2_2d.push_back("MTvBinvbf_PDFidx5"); rebin.push_back(1); legName.push_back("#alpha_{VBF LP}");
  name1_2d.push_back("MTvBinggF_PDFidx2"); name2_2d.push_back("MTvBinggF_PDFidx6"); rebin.push_back(1); legName.push_back("#alpha_{ggF HP}");
  name1_2d.push_back("MTvBinggF_PDFidx3"); name2_2d.push_back("MTvBinggF_PDFidx7"); rebin.push_back(1); legName.push_back("#alpha_{ggF LP}");
  int nVars = 9;
  if(name1_2d.size() && name1_2d[0].Contains("PDF")) nVars = 100;
  
  TLegend *leg[legName.size()];
  TCanvas *c_dphi[name1.size()+2*name1_2d.size()];
  for(int i=0;i<name1.size();i++){
    name = name1[i]+name2[i];
    c_dphi[i]=new TCanvas(name,name,1500,800);
    TH1D *h_histG,*h_histE,*h_histGcopy;
    c_dphi[i]->cd();
    h_histG=(TH1D*)f->FindObjectAny(name1[i]);
    h_histE=(TH1D*)f->FindObjectAny(name2[i]);
    
    if(h_histG && h_histE){
      h_histG->Rebin(rebin[i]);
      h_histE->Rebin(rebin[i]);
      //	h_histG->SetBinContent(h_histG->GetNbinsX(),h_histG->GetBinContent(h_histG->GetNbinsX())+h_histG->GetBinContent(h_histG->GetNbinsX()+1));
      h_histGcopy=(TH1D*)h_histG->Clone("alpha");
      //	h_histE->SetBinContent(h_histE->GetNbinsX(),h_histE->GetBinContent(h_histE->GetNbinsX())+h_histE->GetBinContent(h_histE->GetNbinsX()+1));
      h_histGcopy->Divide(h_histE);
      h_histGcopy->SetLineWidth(2);
      h_histGcopy->SetTitle(";m_{T} (GeV);#alpha");
      h_histGcopy->Draw("e1");
      h_histGcopy->SetLineColor(kBlue);
	
      leg[i] = new TLegend(0.72,0.78,0.86,0.86);
      leg[i]->AddEntry(h_histGcopy,legName[i],"epl");
      leg[i]->SetBorderSize(0);
      leg[i]->Draw();
	
    }
  }
  
  TH2D *h2_histG,*h2_histE,*h2_histGcopy;
  int canvasIdx = name1.size();
  TH1D *h_varRatio[name1_2d.size()][nVars],*h_temp;
  TPaveText *tx[name1_2d.size()];
  for(int i=0;i<name1_2d.size();i++){
    name = name1_2d[i]+name2_2d[i];
    c_dphi[canvasIdx]=new TCanvas(name,name,1500,800);//c_dphi[i]->Divide(4,2);
    c_dphi[canvasIdx]->cd();
    canvasIdx++; 
    h2_histG=(TH2D*)f->FindObjectAny(name1_2d[i]);
    inclOverflow(h2_histG);
    h2_histE=(TH2D*)f->FindObjectAny(name2_2d[i]);
    inclOverflow(h2_histE);

    tx[i] = new TPaveText(0.15,0.8,0.25,0.87,"NDC");
    tx[i]->SetFillColor(0);
    tx[i]->SetShadowColor(0);
    tx[i]->AddText(getCatName(h2_histG->GetTitle()));

    if(h2_histG && h2_histE){
      name = "Variations_"+name2_2d[i];
      h2_histGcopy=(TH2D*)h2_histG->Clone(name);
      h2_histGcopy->Divide(h2_histE);
      h2_histGcopy->Draw("colz");
      gPad->Update();
      if(h2_histGcopy) h2_histGcopy->Write();
      nameX = "temp";
      h_temp = (TH1D*)h2_histGcopy->ProjectionX(nameX,1,1);
      for(int p=1;p<=h_temp->GetNbinsX();p++){
	h_temp->SetBinError(p,0);
      }
      for(int k=0;k<nVars;k++){
	nameX = to_string(k) + "wrt0_" + name1_2d[i];
	h_varRatio[i][k] = (TH1D*)h2_histGcopy->ProjectionX(nameX,k+1,k+1);
	//	h_varRatio[i][k]->Add(h_temp,-1);
	h_varRatio[i][k]->Divide(h_temp);
	if(nameX.Contains("PDF")) h_varRatio[i][k]->SetTitle(";m_{T} (GeV);PDF variation wrt nominal");
	else h_varRatio[i][k]->SetTitle(";m_{T} (GeV);#mu_{R}, #mu_{F} variation wrt nominal");
	//	h_varRatio[i][k]->GetYaxis()->SetRangeUser(0.7499,1.2499);
	h_varRatio[i][k]->GetYaxis()->SetRangeUser(0.9,1.1);
	//	h_varRatio[i][k]->Scale(100);
	h_varRatio[i][k]->Write();
      }
    }
  }
  for(int i=0;i<name1_2d.size();i++){
    //    nameX = to_string(k) + "wrt1_" + name1_2d[i];
    nameX = "VariationRatio"+name1_2d[i];
    c_dphi[name1.size()+name1_2d.size()+i]=new TCanvas(nameX,nameX,1500,800);
    c_dphi[name1.size()+name1_2d.size()+i]->cd();
    for(int k=0;k<nVars;k++){
      if(nVars<10) h_varRatio[i][k]->SetLineColor(col[k]);
      else h_varRatio[i][k]->SetLineColor(k+1);
      if(k==0){
	h_varRatio[i][k]->SetFillStyle(3001);
	h_varRatio[i][k]->SetFillColor(kGray);
	h_varRatio[i][k]->Draw("E2SAME");
      }
      //h_varRatio[i][k]->Draw("histe");
      else h_varRatio[i][k]->Draw("hist same");
      tx[i]->Draw();
      //      if(i==name1_2d.size()-1){
      char name2[100];
      textOnTop.SetTextSize(0.05);
      intLumiE.SetTextSize(0.05);
      textOnTop.DrawLatexNDC(0.12,0.91,"CMS #it{#bf{Simulation Preliminary}}");
      sprintf(name2,"#bf{%0.1f fb^{-1} (13 TeV)}",intLumi);
      intLumiE.DrawLatexNDC(0.7,0.91,name2);
      //    }
    }
  }
  
  gStyle->SetTextSize(2);
  fout->cd();
  //if(h2_histGcopy) h2_histGcopy->Write();
  
  gStyle->SetTextSize(2);
  fout->cd();
  // nameX = to_string(k) + "wrt1_" + name1_2d[i];
  // h_varRatio[k-1] = (TH1D*)h2_histGcopy->ProjectionX(nameX,k,k);
  // h_varRatio[k-1]->Divide((TH1D*)h2_histGcopy->ProjectionX(nameX,1,1));
  // c_dphi[canvasIdx]=new TCanvas(nameX,nameX,1500,800);
  // c_dphi[canvasIdx]->cd();
  // h_varRatio[k-1]->Draw();
  if(savepdf){
    for(int i=0;i<name1.size()+2*name1_2d.size();i++){
      name = c_dphi[i]->GetName();
      if(iFname.Contains("2016")) name = name+"_2016.pdf";
      if(iFname.Contains("2017")) name = name+"_2017.pdf";
      if(iFname.Contains("2018")) name = name+"_2018.pdf";
      c_dphi[i]->SaveAs(name);
    }
  }

}



void inclOverflow(TH2D* h){
  int nx = h->GetNbinsX(), ny = h->GetNbinsY();
  TH1D *h1 = (TH1D*)h->ProjectionX("h1",ny,ny+1);
  TH1D *h2 = (TH1D*)h->ProjectionY("h2",nx,nx+1);

  for(int i=0;i<=nx+1;i++){
    h->SetBinContent(i,ny,h1->GetBinContent(i));
    h->SetBinError(i,ny,h1->GetBinError(i));
  }
  for(int i=0;i<=ny+1;i++){
    h->SetBinContent(nx,i,h2->GetBinContent(i));
    h->SetBinError(nx,i,h2->GetBinError(i));  
  }
}

TString getCatName(TString s){
  if(s.Contains("VBF_HP")) return "VBF HP";
  if(s.Contains("VBF_LP")) return "VBF LP";
  if(s.Contains("ggF_HP")) return "ggF HP";
  if(s.Contains("ggF_LP")) return "ggF LP";
  return s;
}
