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
int col[11]={kBlack,kBlue,kCyan,kMagenta,kRed,kOrange,kGray+1,kMagenta+2,kYellow+2,kMagenta,kCyan};  //Specify Colors
TString name;
TString nameX;

void inclOverflow(TH2D*);
TH2D rebin2DHist(TH2D*,int,double*,int,double*);
TString getCatName(TString);
TLatex textOnTop,intLumiE;
double intLumi=35.9;
bool savepdf=1;

void compAlpha(){
  TH1::SetDefaultSumw2(1);
  gStyle->SetOptStat(0);
  const int nFiles = 5;
  TString iFname[nFiles];
  TFile *f[nFiles];
  vector<TString> legName;
  
  iFname[0] = "VJets_MC2016.root"; legName.push_back("Nominal");
  iFname[1] = "JECup48_VJets_MC2016.root"; legName.push_back("JEC up");
  iFname[2] = "JECdown48_VJets_MC2016.root"; legName.push_back("JEC down");
  iFname[3] = "JERup48_VJets_MC2016.root"; legName.push_back("JER up");
  iFname[4] = "JERdown48_VJets_MC2016.root"; legName.push_back("JER down");

  for(int i=0;i<nFiles;i++){
    f[i] = new TFile(iFname[i]);
  }
      
  name="alpha_comp_"+iFname[0]+iFname[1];
  //  TFile *fout=new TFile(name,"RECREATE");
  if(iFname[0].Contains("2016")) intLumi=35.8;
  else if(iFname[0].Contains("2017")) intLumi=41.5;
  else if(iFname[0].Contains("2018")) intLumi=59.5;

  
  vector<TString> name1,name2;
  vector<TString> catName;
  name1.push_back("MTvBinvbf_0");          name2.push_back("MTvBinvbf_4");        catName.push_back("VBF HP");
  name1.push_back("MTvBinvbf_1");          name2.push_back("MTvBinvbf_5");        catName.push_back("VBF LP");
  name1.push_back("MTvBinggF_2");          name2.push_back("MTvBinggF_6");        catName.push_back("ggF HP");
  name1.push_back("MTvBinggF_3");          name2.push_back("MTvBinggF_7");        catName.push_back("ggF LP");
 
  TLegend *leg = new TLegend(0.72,0.5,0.86,0.86);
  TCanvas *c_dphi[name1.size()];
  TPad *p_top[name1.size()];
  TPad *p_bot[name1.size()];
  TLatex tx[name1.size()];
  
  TH1D *h_sr[nFiles][name1.size()],*h_sb[nFiles][name1.size()],*h_alpha[nFiles][name1.size()],*h_ratio[nFiles-1][name1.size()];
  for(int i=0;i<name1.size();i++){
    name = name1[i]+name2[i];
    c_dphi[i]=new TCanvas(name,name,1500,800);
    p_top[i]=new TPad(name+"_top",name+"_top",0,0.5,1,1);
    p_bot[i]=new TPad(name+"_bot",name+"_bot",0,0.0,1,0.5);
    p_top[i]->Draw();p_top[i]->SetGridx();p_top[i]->SetGridy();
    p_top[i]->SetBottomMargin(0);
    p_bot[i]->SetTopMargin(0);
    p_bot[i]->SetBottomMargin(0.3);
    p_bot[i]->Draw();p_bot[i]->SetGridx();p_bot[i]->SetGridy();
    c_dphi[i]->cd();    p_top[i]->cd();
    for(int p=0;p<nFiles;p++){
      if(p==0){
	h_sr[0][i] = (TH1D*)f[0]->Get(name1[i]);
	h_sb[0][i] = (TH1D*)f[0]->Get(name2[i]);
	name = "alpha0"+to_string(i);
	h_alpha[0][i] = (TH1D*)h_sr[0][i]->Clone(name);
	h_alpha[0][i]->Divide(h_sb[0][i]);
	h_alpha[0][i]->SetLineColor(col[0]);
	h_alpha[0][i]->SetMinimum(0.01);
	h_alpha[0][i]->SetTitle(";mT (GeV);#alpha");
	h_alpha[0][i]->Draw("e1");
	h_alpha[0][i]->GetYaxis()->SetTitleOffset(0.45);
	h_alpha[0][i]->GetYaxis()->SetTitleSize(0.1);
	h_alpha[0][i]->GetYaxis()->SetLabelSize(0.1);
	h_alpha[0][i]->GetYaxis()->SetNdivisions(505);
      }
      else{
	h_sr[p][i] = (TH1D*)f[p]->Get(name1[i]);
	h_sb[p][i] = (TH1D*)f[p]->Get(name2[i]);
	name = "alpha"+to_string(p)+to_string(i);
	h_alpha[p][i] = (TH1D*)h_sr[p][i]->Clone(name);
	h_alpha[p][i]->Divide(h_sb[p][i]);
	h_alpha[p][i]->SetLineColor(col[p]);
	h_alpha[p][i]->SetLineWidth(2);
	h_alpha[p][i]->SetTitle(0);
	c_dphi[i]->cd();    p_top[i]->cd();
	h_alpha[p][i]->Draw("e1 same");
	c_dphi[i]->cd();    p_bot[i]->cd();
	h_ratio[p-1][i]=(TH1D*)h_alpha[p][i]->Clone("ratio"+name);
	h_ratio[p-1][i]->Divide(h_alpha[0][i]);
	h_ratio[p-1][i]->SetTitle(0); name=name1[i];
	h_ratio[p-1][i]->GetXaxis()->SetTitleOffset(0.96);
	h_ratio[p-1][i]->GetXaxis()->SetTitleSize(0.1);
	h_ratio[p-1][i]->GetXaxis()->SetLabelSize(0.1);
	
	h_ratio[p-1][i]->GetYaxis()->SetTitleOffset(0.4);
	h_ratio[p-1][i]->GetYaxis()->SetTitleSize(0.1);
	h_ratio[p-1][i]->GetYaxis()->SetLabelSize(0.1);
	h_ratio[p-1][i]->GetYaxis()->SetNdivisions(505);
	h_ratio[p-1][i]->SetMaximum(1.22);
	h_ratio[p-1][i]->SetMinimum(0.78);
	h_ratio[p-1][i]->SetTitle(";mT (GeV);Nominal/Variation");
	if(p<2) h_ratio[p-1][i]->Draw("hist");
	else h_ratio[p-1][i]->Draw("hist same");
      }
      if(i==0){
	leg->AddEntry(h_alpha[p][i],legName[p],"epl");
	leg->SetBorderSize(0);
      }
      
    }
    c_dphi[i]->cd();    p_top[i]->cd();
    leg->Draw();
    tx[i].SetTextSize(0.1);
    tx[i].DrawLatexNDC(0.15,0.1,catName[i]);
    char name3[100];
    textOnTop.SetTextSize(0.08);
    intLumiE.SetTextSize(0.08);
    textOnTop.DrawLatexNDC(0.12,0.91,"CMS #it{#bf{Simulation Preliminary}}");
    sprintf(name3,"#bf{%0.1f fb^{-1} (13 TeV)}",intLumi);
    intLumiE.DrawLatexNDC(0.7,0.91,name3);
  }
  for(int i=0;i<name1.size() && savepdf;i++){
    name = catName[i];
    name.ReplaceAll(" ","_");
    if(iFname[0].Contains("2016")) name = name+"_JEscaleResUnc_2016.pdf";
    if(iFname[0].Contains("2017")) name = name+"_JEscaleResUnc_2017.pdf";
    if(iFname[0].Contains("2018")) name = name+"_JEscaleResUnc_2018.pdf";
    c_dphi[i]->SaveAs(name);
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
