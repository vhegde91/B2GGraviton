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

const int nfiles=4;    //Specify no. of files
TFile *f[nfiles];
//int col[10]={kBlack,kPink-2,kOrange,kYellow,kGreen,kTeal+9,kPink+1,kCyan,kBlue,kRed};  //Specify Colors
int col[10]={kBlack,kCyan,kBlue,kTeal+9,kGreen,kTeal+9,kPink+1,kCyan,kBlue,kRed};  //Specify Colors

//int col[10]={kBlack,kBlue,kOrange,kYellow,kGreen,kCyan,kPink+1,kTeal+9,kBlue,kRed};  //Specify Colors
//int col[10]={kBlack,kPink-2,kGreen,kYellow,kBlue,kCyan,kPink+1,kCyan,kBlue,kRed};  //Specify Colors

TString name;
bool saveCanvas=0;
void setLastBinAsOverFlow(TH1D*);
void setMyRange(TH1D*,double,double);
void setMyRange(THStack*,double,double);
TString getLegName(TString);
TString getXaxisName(TString);
double data_Integral=0,mc_Integral=0;
TLatex textOnTop,intLumiE;
double intLumi=35.9;
void stackedDataVsMC(TString iFname){
  TH1::SetDefaultSumw2(1);
  gStyle->SetOptStat(0);
  //  gStyle->SetOptStat("nemri"); 
  f[0] = new TFile(iFname);
  if(iFname.Contains("2016")){
    intLumi=35.8;
    f[1] = new TFile("TTJets_MC2016.root");
    f[2] = new TFile("WJetsToLNu_HT_MC2016.root");
    f[3] = new TFile("ZJetsToNuNu_HT_MC2016.root");
  }
  else if(iFname.Contains("2017")){
    intLumi=41.5;
    f[1] = new TFile("TTJets_MC2017.root");
    f[2] = new TFile("WJetsToLNu_HT_MC2017.root");
    f[3] = new TFile("ZJetsToNuNu_HT_MC2017.root");
  }
  else if(iFname.Contains("2018")){
    intLumi=59.5;
    f[1] = new TFile("TTJets_MC2018.root");
    f[2] = new TFile("WJetsToLNu_HT_MC2018.root");
    f[3] = new TFile("ZJetsToNuNu_HT_MC2018.root");
  }
  
  vector<TString> name1;
  vector<int> rebin;
  vector<double> xLow,xHigh;
  name1.push_back("MET_6");  rebin.push_back(5);   xLow.push_back(0); xHigh.push_back(1500);
  name1.push_back("AK8J1Pt_6");   rebin.push_back(5);  xLow.push_back(-300000); xHigh.push_back(300000);
  name1.push_back("AK8J1Eta_6");  rebin.push_back(5);   xLow.push_back(-300000); xHigh.push_back(3000);
  //  name1.push_back("AK8J1Phi_6");  rebin.push_back(2);   xLow.push_back(-300000); xHigh.push_back(3000);
  name1.push_back("AK8J1Mass_6");  rebin.push_back(2);   xLow.push_back(0); xHigh.push_back(3000);
  name1.push_back("AK8J1Tau21_6");  rebin.push_back(1);   xLow.push_back(0); xHigh.push_back(3000);
  name1.push_back("MTvBinggF_6");  rebin.push_back(1);   xLow.push_back(0); xHigh.push_back(300000);
  name1.push_back("MTvBinvbf_6");  rebin.push_back(1);   xLow.push_back(0); xHigh.push_back(300000);
  name1.push_back("dPhi1_6");  rebin.push_back(1);   xLow.push_back(-300000); xHigh.push_back(3000);
  name1.push_back("dPhi4_6");  rebin.push_back(1);   xLow.push_back(-300000); xHigh.push_back(3000);

  TLegend *legend[name1.size()];//=new TLegend(0.6, 0.90,  0.98, 0.45);
  TCanvas *c_cA[name1.size()];
  TPad *p_top[name1.size()];
  TPad *p_bot[name1.size()];
  THStack *hs_hist[name1.size()];
  TH1D *h_sum[name1.size()];
  // gStyle->SetTitleYSize(0.08);
  // gStyle->SetTitleYOffset(0.68);
  for(int i=0;i<name1.size();i++){
    name=name1[i];//+name2[i];
    c_cA[i]=new TCanvas(name,name,1500,800);//c_cA[i]->Divide(4,2);
    p_top[i]=new TPad(name+"_top",name+"_top",0,0.4,1,1);
    p_bot[i]=new TPad(name+"_bot",name+"_bot",0,0.0,1,0.4);
    p_top[i]->Draw();p_top[i]->SetGridx();p_top[i]->SetGridy();p_top[i]->SetLogy();
    p_top[i]->SetBottomMargin(0);
    p_bot[i]->SetTopMargin(0);
    p_bot[i]->SetBottomMargin(0.3);
    p_bot[i]->Draw();p_bot[i]->SetGridx();p_bot[i]->SetGridy();    
    name=name1[i]+"_Stack";
    hs_hist[i] = new THStack(name,name);
    //legend[i]=new TLegend(0.55, 0.88,  0.85, 0.63);
    legend[i]=new TLegend(0.746996,0.5417207,0.8638184,0.8767834);
    legend[i]->SetBorderSize(0);
    legend[i]->SetFillColor(0);
  }
  //cout<<getLegName(f[0]->GetName());
  TH1D *h_histG,*h_histE,*h_histGcopy;
  for(int i=0;i<name1.size();i++){
    c_cA[i]->cd();
    name=name1[i];
    h_histG=(TH1D*)f[0]->FindObjectAny(name);
    h_histG->Rebin(rebin[i]);
    setLastBinAsOverFlow(h_histG);
    setMyRange(h_histG,xLow[i],xHigh[i]);
    if(xLow[i] > -5000 && xHigh[i] < 5000) h_histG->GetXaxis()->SetRangeUser(xLow[i],xHigh[i]);
    if(name1[i].Contains("dPhi1")) data_Integral=h_histG->Integral();
    c_cA[i]->cd();p_top[i]->cd();
    p_top[i]->SetTickx();p_top[i]->SetTicky();
    h_histG->SetLineColor(kBlack);
    h_histG->SetLineWidth(2);
    h_histG->SetMarkerStyle(20);
    h_histG->SetMarkerColor(h_histG->GetLineColor());
    // h_histG->SetTitle(";;Events");
    // h_histG->GetYaxis()->SetTitleOffset(0.60);    
    // h_histG->GetYaxis()->SetTitleSize(0.08);
    // h_histG->GetYaxis()->SetLabelSize(0.08);
    //h_histG->Draw();
    legend[i]->AddEntry(h_histG,"Data","lep");
    //gPad->Update();
    // gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
      
    for(int p=1;p<nfiles;p++){  
      name=name1[i];
      h_histE=(TH1D*)f[p]->FindObjectAny(name);
      h_histE->Rebin(rebin[i]);
      //      h_histE->Scale(29.0/457.519);
      setLastBinAsOverFlow(h_histE);
      setMyRange(h_histE,xLow[i],xHigh[i]);
      if(xLow[i] > -5000 && xHigh[i] < 5000) h_histE->GetXaxis()->SetRangeUser(xLow[i],xHigh[i]);
      if(name1[i].Contains("dPhi1")){
        mc_Integral+=h_histE->Integral();
        cout<<f[p]->GetName()<<" # events "<<h_histE->Integral()<<endl;
      }
      h_histE->SetLineColor(kBlack);
      h_histE->SetLineWidth(1);
      //      h_histE->SetMarkerStyle(21);
      h_histE->SetMarkerColor(col[p]);
      h_histE->SetFillColor(col[p]);
     
      hs_hist[i]->Add(h_histE);
      if(p==1){
	name=name1[i]+"_Sum";
	h_sum[i] = (TH1D*)h_histE->Clone(name);
      }
      else h_sum[i]->Add(h_histE);
      legend[i]->AddEntry(h_histE,getLegName(f[p]->GetName()),"f");
    }
    c_cA[i]->cd();    p_top[i]->cd();
    hs_hist[i]->SetMinimum(0.8);
    hs_hist[i]->SetMaximum(5*hs_hist[i]->GetMaximum());
    hs_hist[i]->Draw("HIST");
    
    h_sum[i]->SetFillColor(1);
    h_sum[i]->SetFillStyle(3013);
    h_sum[i]->Draw("e2 same");

    name = h_histG->GetName();
    if(xLow[i] > -5000 && xHigh[i] < 5000) hs_hist[i]->GetXaxis()->SetRangeUser(xLow[i],xHigh[i]);
    hs_hist[i]->GetYaxis()->SetNdivisions(5);
    //    if(name.Contains("MuPt")) setMyRange(hs_hist[i],xLow[i],xHigh[i]);
    h_histG->Draw("e1 same");
    hs_hist[i]->SetTitle(";;Events");
    hs_hist[i]->GetYaxis()->SetTitleOffset(0.50);    
    hs_hist[i]->GetYaxis()->SetTitleSize(0.09);
    hs_hist[i]->GetYaxis()->SetLabelSize(0.09);
    
    c_cA[i]->Modified();
    c_cA[i]->Update();
    //    legend[i]->SetNColumns(3);
    legend[i]->Draw();

    TH1D *h_numr=(TH1D*)h_histG->Clone();
    h_numr->Divide(h_sum[i]);
    h_numr->SetLineColor(kBlack);
    h_numr->SetMarkerColor(kBlack);
    h_numr->SetTitle(0); name=name1[i];
    h_numr->GetXaxis()->SetTitle(getXaxisName(name));
    h_numr->GetXaxis()->SetTitleOffset(0.96);
    h_numr->GetXaxis()->SetTitleSize(0.13);
    h_numr->GetXaxis()->SetLabelSize(0.13);

    h_numr->GetYaxis()->SetTitle("#frac{Data}{MC}");
    h_numr->GetYaxis()->SetTitleOffset(0.35);
    h_numr->GetYaxis()->SetTitleSize(0.13);
    h_numr->GetYaxis()->SetLabelSize(0.14);
    h_numr->GetYaxis()->SetNdivisions(505);
    h_numr->SetMaximum(1.99);
    h_numr->SetMinimum(0.01);
    c_cA[i]->cd();    p_bot[i]->cd();
    p_bot[i]->SetTickx();p_bot[i]->SetTicky();
    //    c_cB->cd(i+1);    p_bot[i]->cd();
    h_numr->Draw("e1");

    c_cA[i]->cd();    p_top[i]->cd(); gPad->RedrawAxis();
    char name2[100];
    textOnTop.SetTextSize(0.08);
    intLumiE.SetTextSize(0.08);
    textOnTop.DrawLatexNDC(0.12,0.91,"CMS #it{#bf{Preliminary}}");
    intLumiE.SetTextSize(0.08);
    sprintf(name2,"#bf{%0.1f fb^{-1} (13 TeV)}",intLumi);
    intLumiE.DrawLatexNDC(0.7,0.91,name2);
    if(saveCanvas){name=name1[i]+".pdf";c_cA[i]->SaveAs(name);}
    
  }
  
  cout<<"Integral in Data "<<data_Integral<<endl
      <<"Integral in MC "<<mc_Integral<<endl;
}

TString getLegName(TString fname){
  if(fname.Contains("TTJets")) return "t#bar{t}";
  else if(fname.Contains("WJets")) return "W+Jets";
  else if(fname.Contains("ZJets")) return "Z+Jets";

  else return fname;
}
TString getXaxisName(TString axname){
  if(axname.Contains("nHadJets")) return "Jets";
  else if(axname.Contains("ST")) return "ST(GeV)";
  else if(axname.Contains("BTags")) return "b tags";
  else if(axname.Contains("ElePt")) return "e pT(GeV)";
  else if(axname.Contains("MuPt")) return "#mu pT(GeV)";
  else if(axname.Contains("MuEta")) return "#mu #eta";
  else if(axname.Contains("PhotonPt")) return "#gamma pT(GeV)";
  else if(axname.Contains("PhotonEta")) return "#gamma #eta";
  else if(axname.Contains("PhotonPhi")) return "#gamma #Phi";
  else if(axname.Contains("mT")) return "mT(#gamma,p_{T}^{miss})(GeV)";
  else if(axname.Contains("dR_MuPho")) return "#DeltaR(#mu,#gamma)";
  else if(axname.Contains("AllSBin")) return "Bin Number";
  else if(axname.Contains("nVtx")) return "No. of primary vertices";
  else if(axname.Contains("dPhi_METjet1") || axname.Contains("dphi1_METjet1")) return "#Delta#Phi_{1}";
  else if(axname.Contains("dPhi_METjet2") || axname.Contains("dphi2_METjet2")) return "#Delta#Phi_{2}";
  else if(axname.Contains("dPhi_METBestPhoton") ) return "#Delta#Phi(MET,#gamma)";
  else if(axname.Contains("dPhi_Muon_Photon") ) return "#Delta#Phi(#mu,#gamma)";
  else if(axname.Contains("QMut") || axname.Contains("Qmut")) return "QMult";
  else if(axname.Contains("MT_Mu")) return "mT(#mu,p_{T}^{miss})(GeV)";
  else if(axname.Contains("MET")) return "p_{T}^{miss}(GeV)";
  else return axname;

}

void setMyRange(TH1D *h1,double xLow,double xHigh){
  double err=0;
  if(xHigh > 13000) return;
  if(xLow < 13000) return;
  int nMax=h1->FindBin(xHigh);
  h1->SetBinContent(nMax,h1->IntegralAndError(nMax,h1->GetNbinsX(),err));
  h1->SetBinError(nMax,err);
  for(int i=nMax+1;i<=h1->GetNbinsX()+1;i++){
    h1->SetBinContent(i,0);
    h1->SetBinError(i,0);
  }
  //  h1->GetXaxis()->SetRangeUser(xLow,xHigh);
}
void setMyRange(THStack *h1,double xLow,double xHigh){
  // double err=0;
  // if(xHigh > 13000) return;
  // int nMax=h1->FindBin(xHigh);
  // h1->SetBinContent(nMax,h1->IntegralAndError(nMax,h1->GetNbinsX(),err));
  // h1->SetBinError(nMax,err);
  // for(int i=nMax+1;i<=h1->GetNbinsX()+1;i++){
  //   h1->SetBinContent(i,0);
  //   h1->SetBinError(i,0);
  // }
  //h1->GetXaxis()->SetRangeUser(xLow,xHigh);
}

void setLastBinAsOverFlow(TH1D* h_hist){
  double lastBinCt =h_hist->GetBinContent(h_hist->GetNbinsX()),overflCt =h_hist->GetBinContent(h_hist->GetNbinsX()+1);
  double lastBinErr=h_hist->GetBinError(h_hist->GetNbinsX()),  overflErr=h_hist->GetBinError(h_hist->GetNbinsX()+1);
  
  if(lastBinCt!=0 && overflCt!=0)
    lastBinErr = (lastBinCt+overflCt)* (sqrt( ((lastBinErr/lastBinCt)*(lastBinErr/lastBinCt)) + ((overflErr/overflCt)*(overflErr/overflCt)) ) );
  
  else if(lastBinCt==0 && overflCt!=0)
    lastBinErr = overflErr;
  else if(lastBinCt!=0 && overflCt==0)
    lastBinErr = lastBinErr;
  else lastBinErr=0;

  lastBinCt = lastBinCt+overflCt;
  h_hist->SetBinContent(h_hist->GetNbinsX(),lastBinCt);
  h_hist->SetBinError(h_hist->GetNbinsX(),lastBinErr);
    
}



/*


	double lastBinCt =h_histG->GetBinContent(h_histG->GetNbinsX()),overflCt =h_histG->GetBinContent(h_histG->GetNbinsX()+1);
	double lastBinErr=h_histG->GetBinError(h_histG->GetNbinsX()),  overflErr=h_histG->GetBinError(h_histG->GetNbinsX()+1);
	lastBinEr = (lastBinCt+overflCt)* (sqrt( ((lastBinErr/lastBinCt)*(lastBinErr/lastBinCt)) + ((overflErr/overflCt)*(overflErr/overflCt)) ) );
	lastBinCt = lastBinCt+overflCt;
	h_histG->SetBinContent(h_histG->GetNbinsX(),h_histG->GetBinContent(h_histG->GetNbinsX())+h_histG->GetBinContent(h_histG->GetNbinsX()+1));
	h_histE->SetBinContent(h_histE->GetNbinsX(),h_histE->GetBinContent(h_histE->GetNbinsX())+h_histE->GetBinContent(h_histE->GetNbinsX()+1));
*/
