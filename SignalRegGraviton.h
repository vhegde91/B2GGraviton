#ifndef SignalRegGraviton_H
#define SignalRegGraviton_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "NtupleVariables.h"
#include "TH1F.h"
#include "TH2.h"
#include <TProfile.h>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TDirectory.h"

class SignalRegGraviton : public NtupleVariables{

 public:
  SignalRegGraviton(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data");
  ~SignalRegGraviton();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *,const char *);
  void     BookHistogram(const char *);
  void print(Long64_t);

  //Variables defined
  bool isMC=true;
  double wt=0,lumiInfb=35.815165;
  vector<double> ggfEdges = {500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2250,2400,2550,2700,2900,3200};
  vector<double> vbfEdges = {500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1850,2100,2350};

  const static int nCategories = 8;  
  vector<TString> categName = {"0: SR_VBF_HP",
			       "1: SR_VBF_LP",
			       "2: SR_ggF_HP",
			       "3: SR_ggF_LP",
			       "4: SB_VBF_HP",
			       "5: SB_VBF_LP",
			       "6: SB_ggF_HP",
			       "7: SB_ggF_LP",
			       "8: ND"};
  TH1D *h_EvtCategory;
  TH1D *h_vetos;
  TH1D *h_filters;
  TH1D *h_MET[nCategories];
  TH1D *h_MHT[nCategories];
  TH1D *h_HT[nCategories];

  TH1D *h_MT[nCategories];
  TH1D *h_MTvBinggf[nCategories];
  TH1D *h_MTvBinvbf[nCategories];

  TH1D *h_AK8J1Pt[nCategories];
  TH1D *h_AK8J1Eta[nCategories];
  TH1D *h_AK8J1Phi[nCategories];
  TH1D *h_AK8J1Mass[nCategories];
  TH1D *h_AK8J1Tau21[nCategories];
  TH1D *h_dhi1[nCategories];
  TH1D *h_dhi2[nCategories];
  TH1D *h_dhi3[nCategories];
  TH1D *h_dhi4[nCategories];
  TH1D *h_Jet1Pt[nCategories];
  TH1D *h_Jet1Eta[nCategories];
  TH1D *h_Jet1Phi[nCategories];

  TH2D *h2_MTggfPDF[nCategories];
  TH2D *h2_MTvbfPDF[nCategories];
  TH2D *h2_MTggfScl[nCategories];
  TH2D *h2_MTvbfScl[nCategories];
  TH1F *h_PDFNorm; TString nameOfHist1FromFile = "PDFNorm";

  TH1F *h_cutflow;
  TFile *oFile;
  
};
#endif

#ifdef SignalRegGraviton_cxx

void SignalRegGraviton::BookHistogram(const char *outFileName) {

  //  char hname[200], htit[200];
  //  double xlow = 0.0,  xhigh = 2000.0;
  //  int nbins = 2000;
  TString name,title;
 
  oFile = new TFile(outFileName, "recreate");
  TH1::SetDefaultSumw2(1);

  h_cutflow = new TH1F("CutFlow","cut flow",25,0,25);
  h_EvtCategory = new TH1D("EvtCategory","Event category",10,0,10);
  h_vetos = new TH1D("vetos","0 means e veto, mu veto and photon veto applied",5,0,5);
  h_filters = new TH1D("Filters","Filters: Bin1 : all nEvnts, other bins: filter pass/fail",10,0,10);

  for(int i=0;i<nCategories;i++){
    name  = "MET_"+to_string(i); title = "MET "+categName[i];
    h_MET[i] = new TH1D(name,title,200,0,2000);

    name  = "MHT_"+to_string(i); title = "MHT "+categName[i];
    h_MHT[i] = new TH1D(name,title,200,0,2000);

    name  = "HT_"+to_string(i); title = "HT "+categName[i];
    h_HT[i] = new TH1D(name,title,100,0,5000);

    name  = "MT_"+to_string(i); title = "m_{T}(MET,AK8J1) "+categName[i];
    h_MT[i] = new TH1D(name,title,300,0,3000);

    name  = "MTvBinggF_"+to_string(i); title = "ggF m_{T}(MET,AK8J1) "+categName[i];
    h_MTvBinggf[i] = new TH1D(name,title,ggfEdges.size()-1,&(ggfEdges[0]));

    name  = "MTvBinvbf_"+to_string(i); title = "vbf m_{T}(MET,AK8J1) "+categName[i];
    h_MTvBinvbf[i] = new TH1D(name,title,vbfEdges.size()-1,&(vbfEdges[0]));

    name  = "AK8J1Pt_"+to_string(i); title = "Pt of leading AK8Jet "+categName[i];
    h_AK8J1Pt[i] = new TH1D(name,title,200,0,2000);

    name  = "AK8J1Eta_"+to_string(i); title = "#eta_{AK8} "+categName[i];
    h_AK8J1Eta[i] = new TH1D(name,title,120,-6,6);

    name  = "AK8J1Phi_"+to_string(i); title = "#Phi_{AK8} "+categName[i];
    h_AK8J1Phi[i] = new TH1D(name,title,80,-4,4);

    name = "AK8J1Mass_"+to_string(i); title = "M_{AK8} "+categName[i];
    h_AK8J1Mass[i] = new TH1D(name,title,60,0,300);

    name = "AK8J1Tau21_"+to_string(i); title = "AK8J #tau_{21} "+categName[i];
    h_AK8J1Tau21[i] = new TH1D(name,title,20,0,1);

    name = "dPhi1_"+to_string(i); title = "#Delta#Phi(AK4J1,MET) "+categName[i];
    h_dhi1[i] = new TH1D(name,title,40,0,4);

    name = "dPhi2_"+to_string(i); title = "#Delta#Phi(AK4J2,MET) "+categName[i];
    h_dhi2[i] = new TH1D(name,title,40,0,4);

    name = "dPhi3_"+to_string(i); title = "#Delta#Phi(AK4J3,MET) "+categName[i];
    h_dhi3[i] = new TH1D(name,title,40,0,4);

    name = "dPhi4_"+to_string(i); title = "#Delta#Phi(AK4J4,MET) "+categName[i];
    h_dhi4[i] = new TH1D(name,title,40,0,4);

    name = "Jet1Pt_"+to_string(i); title = "Pt of leading jet "+categName[i];
    h_Jet1Pt[i] = new TH1D(name,title,200,0,2000);

    name = "Jet1Eta_"+to_string(i); title = "#eta of leading jet "+categName[i];
    h_Jet1Eta[i] = new TH1D(name,title,120,-6,6);

    name = "Jet1Phi_"+to_string(i); title = "#Phi of leading jet "+categName[i];
    h_Jet1Phi[i] = new TH1D(name,title,80,-4,4);

    name = "MTvBinvbf_PDFidx"+to_string(i); title = "x:MT VBF, y:PDF index "+categName[i];
    h2_MTvbfPDF[i] = new TH2D(name,title,vbfEdges.size()-1,&(vbfEdges[0]),110,0.,110.);

    name = "MTvBinggF_PDFidx"+to_string(i); title = "x:MT ggF, y:PDF index "+categName[i];
    h2_MTggfPDF[i] = new TH2D(name,title,ggfEdges.size()-1,&(ggfEdges[0]),110,0.,110.);

    name = "MTvBinvbf_Sclidx"+to_string(i); title = "x:MT VBF, y:Scale index "+categName[i];
    h2_MTvbfScl[i] = new TH2D(name,title,vbfEdges.size()-1,&(vbfEdges[0]),10,0.,10.);

    name = "MTvBinggF_Sclidx"+to_string(i); title = "x:MT ggF, y:Scale index "+categName[i];
    h2_MTggfScl[i] = new TH2D(name,title,ggfEdges.size()-1,&(ggfEdges[0]),10,0.,10.);
  }
}

SignalRegGraviton::SignalRegGraviton(const TString &inputFileList, const char *outFileName, const char* dataset) {
  string nameData=dataset;
  TChain *tree = new TChain("tree");
  if(nameData=="signal") tree = new TChain("TreeMaker2/PreSelection");
  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
  }

  if(nameData!="signalH") nameData="BG";
  if(nameData=="signalH") nameData="signal";
  cout<<"Treating the input files as "<<nameData<<" for setting tree branches"<<endl;
  NtupleVariables::Init(tree,nameData);

  BookHistogram(outFileName);
  
}

Bool_t SignalRegGraviton::FillChain(TChain *chain, const TString &inputFileList) {
  int itr=0;
  TFile *filePointer;
  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1) {
    infile >> buffer;
    if(!infile.good()) break;
    //    std::cout << "Adding tree from " << buffer.c_str() << std::endl;

    if(isMC){
      filePointer = TFile::Open(buffer.c_str());
      if(itr==0) h_PDFNorm = (TH1F*)filePointer->Get("PDFNorm");
      else h_PDFNorm->Add((TH1F*)filePointer->Get("PDFNorm"));
      itr++;
    }

    chain->Add(buffer.c_str());
  }
  std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
  return kTRUE;
}

Long64_t SignalRegGraviton::LoadTree(Long64_t entry) {
  // Set the environment to read one entry                                                                                          
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
  }
  return centry;
}

SignalRegGraviton::~SignalRegGraviton() { 

  if (!fChain) return;
  delete fChain->GetCurrentFile();
  oFile->cd();
  oFile->Write();
  if(isMC) h_PDFNorm->Write();
  oFile->Close();

}

#endif

