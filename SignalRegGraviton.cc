#define SignalRegGraviton_cxx
#include "SignalRegGraviton.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <cstring>
#include <string>
#include <fstream>
#include "btag/BTagCorrector.h"

using namespace std;

int main(int argc, char* argv[])
{

  if (argc < 2) {
    cerr << "Please give 3 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" << endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];

  SignalRegGraviton ana(inputFileList, outFileName, data);
  cout << "dataset " << data << " " << endl;

  ana.EventLoop(data,inputFileList);

  return 0;
}

void SignalRegGraviton::EventLoop(const char *data,const char *inputFileList) {
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;

  TString s_data=data, inFileList = inputFileList;
  Long64_t nbytes = 0, nb = 0;
  int decade = 0;
  Long64_t nEvtSurv = 0;
  int ak8J1Idx = -1;
  TLorentzVector metVec;  
  vector<TLorentzVector> myjets;
  int isSB=0;
  vector<int> categNum;
  bool applySDmassCorr = 1, applyMadHTless600 = 0, applyPUwt = 1, applybTagSFs = 1;
  TString puWtFileName;
  TFile *f_puwt;

  for(int i=0;i<categName.size();i++) h_EvtCategory->Fill(categName[i],0);
  h_cutflow->Fill("0",0);
  h_cutflow->Fill("Weighted",0);    
  h_cutflow->Fill("MET>200",0);
  h_cutflow->Fill("AK8L1JPt>200",0);
  h_cutflow->Fill("dPhiCuts",0);
  h_cutflow->Fill("photonVeto",0);
  h_cutflow->Fill("LVeto",0);
  h_cutflow->Fill("bVeto",0);
  h_cutflow->Fill("Filters",0);
  h_cutflow->Fill("HTRatioDphi",0);
  h_cutflow->Fill("PassedTrigger",0);
  h_cutflow->Fill("SBorSR_anyPurity",0);
  h_cutflow->Fill("SR_anyPurity",0);
  h_cutflow->Fill("SB_anyPurity",0);
  h_cutflow->Fill("SR+HP",0);
  h_cutflow->Fill("SR+LP",0);
  h_cutflow->Fill("SB+HP",0);
  h_cutflow->Fill("SB+LP",0);
  h_cutflow->Fill("isVBF",0);
  h_cutflow->Fill("HEMaffected",0);
  h_cutflow->Fill("NEvtsNoWtLeft",0);

  h_filters->Fill("TotEvnts",0);
  h_filters->Fill("globalSuperTightHalo2016Filter",0);
  h_filters->Fill("HBHENoiseFilter",0);
  h_filters->Fill("HBHEIsoNoiseFilter",0);
  h_filters->Fill("eeBadScFilter",0);
  h_filters->Fill("EcalDeadCellTriggerPrimitiveFilter",0);
  h_filters->Fill("BadChargedCandidateFilter",0);
  h_filters->Fill("BadPFMuonFilter",0);
  h_filters->Fill("NVtx>0",0);
  h_filters->Fill("JetID",0);
  h_filters->Fill("(MET/CaloMET<5.)",0);

  double p0=-100, p1=-100, p2=-100, deepCSVvalue = 0.0, tempVar = 0.;
  int dataRun = 0;
  if(s_data.Contains("MC_2016")){ dataRun = -2016;      lumiInfb = 35.815165; deepCSVvalue = 0.6321; p0 = 1.21277e+02; p1 = 8.77679e+01; p2 = 9.94172e-01;}
  else if(s_data.Contains("MC_2017")){ dataRun = -2017; lumiInfb = 41.486136; deepCSVvalue = 0.4941; p0 = 1.61724e+02; p1 = 6.91644e+01; p2 = 9.89446e-01;}
  else if(s_data.Contains("MC_2018")){ dataRun = -2018; lumiInfb = 59.546381; deepCSVvalue = 0.4184; p0 = 1.70454e+02; p1 = 6.64100e+01; p2 = 9.89298e-01;}

  else if(s_data.Contains("2016")){ dataRun = 2016; isMC = false; deepCSVvalue = 0.6321;}
  else if(s_data.Contains("2017")){ dataRun = 2017; isMC = false; deepCSVvalue = 0.4941;}
  else if(s_data.Contains("2018")){ dataRun = 2018; isMC = false; deepCSVvalue = 0.4184;}

  //----------PU reweighting
  puWtFileName = "PileupHistograms_"+to_string(abs(dataRun))+"_69mb_pm5.root";
  f_puwt = TFile::Open(puWtFileName,"READ");
  TH1* puhist = (TH1*)f_puwt->Get("pu_weights_central");
  //----------btag SFs
  int fListIndxOld=-1;
  vector<TString> inFileName;
  TString sampleName;
  string str1;
  ifstream runListFile(inputFileList);
  TFile *currFile;
  while (std::getline(runListFile, str1)) {
    inFileName.push_back(str1);
  }runListFile.close();
  cout<<"applying b-tag SFs? "<<applybTagSFs<<endl;
  BTagCorrector btagcorr;
  //----------mad HT cuts for tt_dilep and single_lep
  if(inFileList.Contains("TTJets_SingleLeptFromT") || inFileList.Contains("TTJets_DiLept")) applyMadHTless600 = true;

  if(dataRun>0) cout<<"Processing it as "<<dataRun<<" data"<<endl;
  else if(dataRun<0) cout<<"Processing it as "<<abs(dataRun)<<" MC"<<endl;
  else cout<<"No specific data/MC year"<<endl;
  cout<<"Apply PU reweighting on MC? "<<applyPUwt<<endl;
  cout<<"Apply SD mass corrections on data & MC? "<<applySDmassCorr<<endl;

  int jec2Use = 0;//-1 for JEC down, 0 for CV, 1 for JEC up
  int jer2Use = 0;//-1 for JER down, 0 for CV, 1 for JER up
  int jet2Vary = 48;//4: only AK4 jets, 8: only AK8 jets, 48 or 84: both AK4 and AK8.
  if(jet2Vary!=0 && (jec2Use!=0 || jer2Use!=0)){
    cout<<"Varying jets like:"<<jet2Vary<<" 4: only AK4 jets, 8: only AK8 jets, 48 or 84: both AK4 and AK8"<<endl;
    if(jec2Use!=0) cout<<"!!!!!!!!!! Applying JECs. -1 for JEC down, 0 for CV, 1 for JEC up. I am using "<<jec2Use<<" !!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    if(jer2Use!=0) cout<<"!!!!!!!!!! Applying JERs. -1 for JER down, 0 for CV, 1 for JER up. I am using "<<jer2Use<<" !!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
  }
  //----- SD mass corrections
  if(applySDmassCorr){
    if(abs(dataRun) == 2016 || abs(dataRun) == 2017 || abs(dataRun) == 2018){
      sdCorrFile = TFile::Open("puppiCorr.root");//https://github.com/cms-jet/PuppiSoftdropMassCorrections/tree/80X
      if(sdCorrFile->IsOpen()){
	puppisd_corrGEN      = (TF1*)sdCorrFile->Get("puppiJECcorr_gen");
	puppisd_corrRECO_cen = (TF1*)sdCorrFile->Get("puppiJECcorr_reco_0eta1v3");
	puppisd_corrRECO_for = (TF1*)sdCorrFile->Get("puppiJECcorr_reco_1v3eta2v5");
      }
      else{
	cout<<"!!!!! can't get SD mass corr file. Will not apply any corrn."<<endl;
	applySDmassCorr = false;
      }
    }
  }


  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    // ==============print number of events done == == == == == == == =
    double progress = 10.0 * jentry / (1.0 * nentries);
    int k = int (progress);
    if (k > decade)
      cout << 10 * k << " %" <<endl;
    decade = k;
    // cout<<"j:"<<jentry<<" fcurrent:"<<fCurrent<<endl;
    // ===============read this entry == == == == == == == == == == == 
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //    print(jentry);  
    if(dataRun==-2016 || dataRun==-2017) wt=Weight*1000.0*lumiInfb*NonPrefiringProb;
    else if(dataRun <=0) wt=Weight*1000.0*lumiInfb;
    else wt = 1.0;
    if(applyPUwt && isMC){
      tempVar = (puhist->GetBinContent(puhist->GetXaxis()->FindBin(min(TrueNumInteractions,puhist->GetBinLowEdge(puhist->GetNbinsX()+1)))));
      if(tempVar < 100) wt = wt * tempVar;
    }

    //-------------- btag SFs
    if(fListIndxOld!=fCurrent && isMC){ 
      fListIndxOld = fCurrent;
      sampleName = inFileName[fCurrent];
      if(applybTagSFs){
	// btagcorr.SetBtagSFunc(1);
	// btagcorr.SetMistagSFunc(1);
	currFile = TFile::Open(sampleName);
	btagcorr.SetEffs(currFile);
	if(dataRun==-2016) btagcorr.SetCalib("btag/DeepCSV_Moriond17_B_H_mod.csv");
	if(dataRun==-2017) btagcorr.SetCalib("btag/DeepCSV_94XSF_V3_B_F_mod.csv");
	if(dataRun==-2018) btagcorr.SetCalib("btag/DeepCSV_102XSF_V1_mod.csv");
      }
    }
    double corrbtag = 1.0;
    if(isMC && applybTagSFs){
      corrbtag = btagcorr.GetSimpleCorrection(Jets,Jets_hadronFlavor,Jets_HTMask,Jets_bJetTagDeepCSVBvsAll,deepCSVvalue);
      wt = wt * corrbtag;
    }

    if(isMC && applyMadHTless600 && madHT > 600) continue;

    h_cutflow->Fill("0",1);
    h_cutflow->Fill("Weighted",wt);
    h_madHT->Fill(madHT,wt);
    //--------------
    if(jec2Use!=0 || jer2Use!=0) changeJets(jec2Use, jer2Use, jet2Vary);
    // cout<<"HT values from tree: Nominal, JECup, JECDn, JERup, JERDn "<<HT<<" "<<HTJECup<<" "<<HTJECdown<<" "<<HTJERup<<" "<<HTJERdown<<endl;
    // float myHT=0;
    // for(int i=0;i<Jets->size();i++){
    //   if((*Jets)[i].Pt() > 30 && abs((*Jets)[i].Eta()) < 2.4) myHT=myHT+(*Jets)[i].Pt();
    // }
    // cout<<"myHT:"<<myHT<<endl;
    //----MET
    if(MET < 200) continue;
    h_cutflow->Fill("MET>200",wt);
    //----AK8Jets
    if(JetsAK8->size()==0) continue;
    ak8J1Idx = -1;
    for(int i=0;i<JetsAK8->size();i++){
      if(ak8J1Idx < 0) ak8J1Idx = 0;
      if((*JetsAK8)[ak8J1Idx].Pt() < (*JetsAK8)[i].Pt()) ak8J1Idx = i;
    }
    if(ak8J1Idx!=0) cout<<"!!!!! AK8Jets are not sorted!!!!"<<endl;
    if((*JetsAK8)[ak8J1Idx].Pt() <= 200 || (abs((*JetsAK8)[ak8J1Idx].Eta()) >=2.5)) continue;
    h_cutflow->Fill("AK8L1JPt>200",wt);

    //    if(abs((*JetsAK8)[ak8J1Idx].Eta()) >= 2.4) continue;
    float dphi1=4, dphi2=4, dphi3=4, dphi4=4;
    if(Jets->size() > 0 && (*Jets)[0].Pt() > 30 && abs((*Jets)[0].Eta()) < 6.0)
      dphi1 = (abs(DeltaPhi(METPhi,(*Jets)[0].Phi())));

    if(Jets->size() > 1 && (*Jets)[1].Pt() > 30 && abs((*Jets)[1].Eta()) < 6.0)
      dphi2 = (abs(DeltaPhi(METPhi,(*Jets)[1].Phi())));

    if(Jets->size() > 2 && (*Jets)[2].Pt() > 30 && abs((*Jets)[2].Eta()) < 6.0)
      dphi3 = (abs(DeltaPhi(METPhi,(*Jets)[2].Phi())));

    if(Jets->size() > 3 && (*Jets)[3].Pt() > 30 && abs((*Jets)[3].Eta()) < 6.0)
      dphi4 = (abs(DeltaPhi(METPhi,(*Jets)[3].Phi())));
    // print(jentry);
    // cout<<"METPhi:"<<METPhi<<" dphi1:"<<dphi1<<" dphi2:"<<dphi2<<" dphi3:"<<dphi3<<" dphi4:"<<dphi4<<endl;
    
    if(!(dphi1 > 0.5 && dphi2 > 0.5 && dphi3 > 0.5 && dphi4 > 0.5)) continue;
    h_cutflow->Fill("dPhiCuts",wt);
    //----Photon veto
    int nPhotons=0;
    for(int i=0;i<Photons->size();i++){
      if((*Photons)[i].Pt() > 100 && (*Photons_fullID)[i] && (!(*Photons_hasPixelSeed)[i]) ){ nPhotons++; break;}
    }
    if(nPhotons>0) continue;
    //    if(Photons->size()!=0) continue;
    h_cutflow->Fill("photonVeto",wt);
    //----Lepton veto
    if(NMuons!=0 || NElectrons!=0) continue;
    if(isoMuonTracks != 0 || isoElectronTracks != 0 || isoPionTracks != 0) continue;
    h_cutflow->Fill("LVeto",wt);
    //----b-tag veto
    if(BTagsDeepCSV!=0) continue;
    h_cutflow->Fill("bVeto",wt);

    h_filters->Fill("TotEvnts",1);
    h_filters->Fill("globalSuperTightHalo2016Filter",globalSuperTightHalo2016Filter);
    h_filters->Fill("HBHENoiseFilter",HBHENoiseFilter);
    h_filters->Fill("HBHEIsoNoiseFilter",HBHEIsoNoiseFilter);
    h_filters->Fill("eeBadScFilter",eeBadScFilter);
    h_filters->Fill("EcalDeadCellTriggerPrimitiveFilter",EcalDeadCellTriggerPrimitiveFilter);
    h_filters->Fill("BadChargedCandidateFilter",BadChargedCandidateFilter);
    h_filters->Fill("BadPFMuonFilter",BadPFMuonFilter);
    h_filters->Fill("NVtx>0",NVtx > 0);
    h_filters->Fill("JetID",JetID);
    h_filters->Fill("(MET/CaloMET<5.)",(MET/CaloMET < 5.));

    if(!(globalSuperTightHalo2016Filter==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadChargedCandidateFilter && BadPFMuonFilter && NVtx > 0 && JetID && (MET/CaloMET < 5.))) continue;
    h_cutflow->Fill("Filters",wt);

    //if(!HTRatioDPhiFilter) continue;
    float htratio = HT5/HT;
    if( !((htratio < 1.2) || ((dphi1) >= ((2.64*htratio)-2.64))) ) continue;
    h_cutflow->Fill("HTRatioDphi",wt);
    //--------------------------triggers
    if(!isMC){
      bool trgPass = false;
      TString trgName;
      for(int i=0;i<TriggerPass->size();i++){
	trgName = (*TriggerNames)[i];
	if(!(trgName.Contains("MET"))) continue;
	if((*TriggerPass)[i]==1 && (trgName.Contains("HLT_PFMET100_PFMHT100_IDTight_v") || trgName.Contains("HLT_PFMET110_PFMHT110_IDTight_v") ||
				    trgName.Contains("HLT_PFMET120_PFMHT120_IDTight_v") || trgName.Contains("HLT_PFMET130_PFMHT130_IDTight_v") ||
				    trgName.Contains("HLT_PFMET140_PFMHT140_IDTight_v") || trgName.Contains("HLT_PFMET90_PFMHT90_IDTight_v")))
	  trgPass = true;
      }
      if(trgPass) h_cutflow->Fill("PassedTrigger",wt);
      else continue;
    }
    else wt = wt*( (TMath::Erf((MET - p0) / p1) + 1) / 2. * p2);//apply trigger eff weights to MC

    bool HEMaffected = false;
    if(dataRun==2018 && RunNum >=319077){//for data 2018
      for(int i=0;i<Jets->size();i++){
	if((*Jets)[i].Pt() < 30) continue;
	if( (*Jets)[i].Eta() >= -3.20 && (*Jets)[i].Eta() <= -1.2 && 
	    (*Jets)[i].Phi() >= -1.77 && (*Jets)[i].Phi() <= -0.67 &&
	    (abs(DeltaPhi(METPhi,(*Jets)[i].Phi())) < 0.5) ){HEMaffected = true; break;}
      }
    }
    if(dataRun==-2018){//for MC 2018
      if( (EvtNum % 1000 > 1000*21.0/59.6) && !passHEMjetVeto(30.)) HEMaffected = true;
    }
    //--------------------------end of triggers
    if(applySDmassCorr) applySDmassCorrAllAK8();//SD mass corrections

    if(((*JetsAK8_softDropMass)[ak8J1Idx] > 65) && ((*JetsAK8_softDropMass)[ak8J1Idx] < 105)){ isSB=0; h_cutflow->Fill("SR_anyPurity",wt);}
    else if( (((*JetsAK8_softDropMass)[ak8J1Idx] > 30) && ((*JetsAK8_softDropMass)[ak8J1Idx] < 65)) || (((*JetsAK8_softDropMass)[ak8J1Idx] > 135) && ((*JetsAK8_softDropMass)[ak8J1Idx] < 300) )){ isSB=1; h_cutflow->Fill("SB_anyPurity",wt);}
    else isSB=-1;
    if(isSB==-1) continue;
    h_cutflow->Fill("SBorSR_anyPurity",wt);
    if(!isSB && !isMC) continue;
    //----HP or LP
    int purityIdx = 0;
    if(((*JetsAK8_NsubjettinessTau2)[ak8J1Idx]/(*JetsAK8_NsubjettinessTau1)[ak8J1Idx]) < 0.35) purityIdx = 2;
    else if(((*JetsAK8_NsubjettinessTau2)[ak8J1Idx]/(*JetsAK8_NsubjettinessTau1)[ak8J1Idx]) < 0.75) purityIdx = 1;
    else{ purityIdx = 0; continue;}
    if(!isSB && purityIdx==2) h_cutflow->Fill("SR+HP",wt);
    if(!isSB && purityIdx==1) h_cutflow->Fill("SR+LP",wt);
    if(isSB && purityIdx==2) h_cutflow->Fill("SB+HP",wt);
    if(isSB && purityIdx==1) h_cutflow->Fill("SB+LP",wt);
    // sortTLorVec(&myjets);
    //----V(J)
    //----MT
    double mt = sqrt(2*(*JetsAK8)[ak8J1Idx].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*JetsAK8)[ak8J1Idx].Phi()))));
    //----VBF
    bool isVBF=true;
    int vbfJ1Idx=-1, vbfJ2Idx=-1;
    TLorentzVector vbfJpair;
    bool ak8JLooksLikeZ=(((*JetsAK8_softDropMass)[ak8J1Idx] > 30 && (*JetsAK8_softDropMass)[ak8J1Idx] < 300) &&
			 ((*JetsAK8_NsubjettinessTau2)[ak8J1Idx]/(*JetsAK8_NsubjettinessTau1)[ak8J1Idx]) < 0.75);
    for(int i=0;i<Jets->size();i++){//choosing VBF pair of jets
      if((*Jets)[i].Pt() < 30) continue;
      if( ((*Jets)[i].DeltaR((*JetsAK8)[ak8J1Idx]) < 0.8) && ak8JLooksLikeZ) continue;
      for(int j=i+1;j<Jets->size();j++){
	if((*Jets)[j].Pt() < 30) continue;
	if( ((*Jets)[j].DeltaR((*JetsAK8)[ak8J1Idx]) < 0.8) && ak8JLooksLikeZ) continue;
	if(vbfJpair.M() < (((*Jets)[i]+(*Jets)[j]).M()) ){
	  vbfJpair = (*Jets)[i]+(*Jets)[j];
	  vbfJ1Idx = i; vbfJ2Idx = j;
	}
      }
    }
    if(vbfJpair.M() < 500) isVBF=false;
    if(vbfJ1Idx >=0 && vbfJ2Idx >=0){
      if(((*Jets)[vbfJ1Idx].Eta())*((*Jets)[vbfJ2Idx].Eta()) > 0) isVBF=false;
      if(abs((*Jets)[vbfJ1Idx].Eta()-(*Jets)[vbfJ2Idx].Eta()) < 4) isVBF=false;
    }
    else isVBF=false;
    if(isVBF) h_cutflow->Fill("isVBF",wt);
    //---------------------------------
    categNum.resize(0);
    if(!isSB){
      if(isVBF && purityIdx==2) categNum.push_back(0);
      if(isVBF && purityIdx==1) categNum.push_back(1);
      if(!isVBF && purityIdx==2) categNum.push_back(2);
      if(!isVBF && purityIdx==1) categNum.push_back(3);
    }
    else{
      if(isVBF && purityIdx==2) categNum.push_back(4);
      if(isVBF && purityIdx==1) categNum.push_back(5);
      if(!isVBF && purityIdx==2) categNum.push_back(6);
      if(!isVBF && purityIdx==1) categNum.push_back(7);
    }

    if(categNum.size()==0 || categNum.size() > 1){ cout<<"Could not IDfy category"<<endl; break;} //at max categs can be VBF or nonVBF. So size <=1.
    if(categNum[0] >=4 && categNum[0] <=7) categNum.push_back(8);//if categNum[0] == 4,5,6,7 then it is Baseline + SB.

    h_vetos->Fill(NMuons+NElectrons+BTags+nPhotons,wt);    

    if(mt < 400) continue;
    if(HEMaffected){
      h_cutflow->Fill("HEMaffected",wt);
      continue;
    }
    
    for(int iCategory=0;iCategory<categNum.size();iCategory++){
      h_EvtCategory->Fill(categName[categNum[iCategory]],wt);
      h_MET[categNum[iCategory]]->Fill(MET,wt);
      h_MHT[categNum[iCategory]]->Fill(MHT,wt);
      h_HT[categNum[iCategory]]->Fill(HT,wt);
      h_METPhi[categNum[iCategory]]->Fill(METPhi,wt);

      h_MT[categNum[iCategory]]->Fill(mt,wt);
      if(mt > ggfEdges[ggfEdges.size()-1]) 
      	h_MTvBinggf[categNum[iCategory]]->Fill(0.5*(ggfEdges[ggfEdges.size()-1]+ggfEdges[ggfEdges.size()-2]),wt);
      else h_MTvBinggf[categNum[iCategory]]->Fill(mt,wt);

      if(mt > vbfEdges[vbfEdges.size()-1]) 
      	h_MTvBinvbf[categNum[iCategory]]->Fill(0.5*(vbfEdges[vbfEdges.size()-1]+vbfEdges[vbfEdges.size()-2]),wt);
      else h_MTvBinvbf[categNum[iCategory]]->Fill(mt,wt);

      h_AK8J1Pt[categNum[iCategory]]->Fill((*JetsAK8)[ak8J1Idx].Pt(),wt);
      h_AK8J1Eta[categNum[iCategory]]->Fill((*JetsAK8)[ak8J1Idx].Eta(),wt);
      h_AK8J1Phi[categNum[iCategory]]->Fill((*JetsAK8)[ak8J1Idx].Phi(),wt);
      h_AK8J1Mass[categNum[iCategory]]->Fill((*JetsAK8_softDropMass)[ak8J1Idx],wt);
      h_AK8J1Tau21[categNum[iCategory]]->Fill(((*JetsAK8_NsubjettinessTau2)[ak8J1Idx])/((*JetsAK8_NsubjettinessTau1)[ak8J1Idx]),wt);
      h_dhi1[categNum[iCategory]]->Fill(dphi1,wt);
      h_dhi2[categNum[iCategory]]->Fill(dphi2,wt);
      h_dhi3[categNum[iCategory]]->Fill(dphi3,wt);
      h_dhi4[categNum[iCategory]]->Fill(dphi4,wt);
      if(myjets.size() >0) {
	h_Jet1Pt[categNum[iCategory]]->Fill(myjets[0].Pt(),wt);
	h_Jet1Eta[categNum[iCategory]]->Fill(myjets[0].Eta(),wt);
	h_Jet1Phi[categNum[iCategory]]->Fill(myjets[0].Phi(),wt);
      }
      h_NVtx[categNum[iCategory]]->Fill(NVtx,wt);
      if(dataRun <=0){
	for(int i=0;i<PDFweights->size();i++){
	  h2_MTggfPDF[categNum[iCategory]]->Fill(mt,i,(*PDFweights)[i]*wt);
	  h2_MTvbfPDF[categNum[iCategory]]->Fill(mt,i,(*PDFweights)[i]*wt);
	}
	for(int i=0;i<ScaleWeights->size();i++){
	  h2_MTggfScl[categNum[iCategory]]->Fill(mt,i,(*ScaleWeights)[i]*wt);
	  h2_MTvbfScl[categNum[iCategory]]->Fill(mt,i,(*ScaleWeights)[i]*wt);
	}
      }
    }//categories
    nEvtSurv++;
    h_cutflow->Fill("NEvtsNoWtLeft",1);
  } // loop over entries
  cout<<"zzzzzzzzzzzzzz No. of entries survived: "<<nEvtSurv<<endl;
}

void SignalRegGraviton::applySDmassCorrAllAK8(){
  float genCorr  = 1.;
  float recoCorr = 1.;
  float totalWeight = 1.;

  TLorentzVector subjetSum;
  uncorrAK8SubjSum.resize(0);
  SDmassCorrFac.resize(0);

  for(int i=0;i< JetsAK8->size();i++){
    for(int j=0;j < (*JetsAK8_subjets)[i].size();j++){
      if(j==0) subjetSum = (*JetsAK8_subjets)[i][j];
      else subjetSum = subjetSum + (*JetsAK8_subjets)[i][j];
      //      cout<<"subj pt,eta,phi,M:"<<(*JetsAK8_subjets)[i][j].Pt()<<", "<<(*JetsAK8_subjets)[i][j].Eta()<<", "<<(*JetsAK8_subjets)[i][j].Phi()<<" ,"<<(*JetsAK8_subjets)[i][j].M()<<endl;
    }
    uncorrAK8SubjSum.push_back(subjetSum);
    SDmassCorrFac.push_back(1);
    
    genCorr =  puppisd_corrGEN->Eval( (*JetsAK8)[i].Pt() );
    if( abs((*JetsAK8)[i].Eta())  <= 1.3 ) 
      recoCorr = puppisd_corrRECO_cen->Eval( (*JetsAK8)[i].Pt() );
    else 
      recoCorr = puppisd_corrRECO_for->Eval( (*JetsAK8)[i].Pt() );
    
    SDmassCorrFac.push_back(genCorr * recoCorr);
    (*JetsAK8_softDropMass)[i] = (subjetSum.M())*genCorr*recoCorr;
    //    cout<<"genCorr:"<<genCorr<<" recoCorr:"<<recoCorr<<" subjetSum.M():"<<subjetSum.M()<<endl;
  }
}


void SignalRegGraviton::print(Long64_t jentry){
  //cout<<endl;
  TLorentzVector v1,photo;
  cout<<"------------------------------------------------------------"<<endl;
  cout<<"MomMass:"<<SusyMotherMass<<" Kid Mass:"<<SusyLSPMass<<endl;
  for(int i=0;i<GenParticles->size();i++){  
    //    cout<<EvtNum<<" "<<jentry<<" "<<GenParticles->size()<<" "<<i<<" PdgId:"<<(*GenParticles_PdgId)[i]<<" parentId:"<<(*GenParticles_ParentId)[i]<<" parentIndx:"<<(*GenParticles_ParentIdx)[i]<<" Status:"<<(*GenParticles_Status)[i]<<"\tPx :"<<(*GenParticles)[i].Px()<<" Py :"<<(*GenParticles)[i].Py()<<" Pz :"<<(*GenParticles)[i].Pz()<<" E: "<<(*GenParticles)[i].Energy()<<" M:"<<(*GenParticles)[i].M()<<endl;
    cout<<EvtNum<<" "<<jentry<<" "<<GenParticles->size()<<" "<<i<<" PdgId:"<<(*GenParticles_PdgId)[i]<<" parentId:"<<(*GenParticles_ParentId)[i]<<" parentIndx:"<<(*GenParticles_ParentIdx)[i]<<" Status:"<<(*GenParticles_Status)[i]<</*"\tPx:"<<(*GenParticles)[i].Px()<<" Py:"<<(*GenParticles)[i].Py()<<" Pz:"<<(*GenParticles)[i].Pz()<<*/"\tPt:"<<(*GenParticles)[i].Pt()<<" Eta:"<<(*GenParticles)[i].Eta()<<" Phi:"<<(*GenParticles)[i].Phi()<<" E:"<<(*GenParticles)[i].Energy()<<endl;

  }
  for(int i=0;i<Photons->size();i++){
    double dR=0;//DeltaR( bestPhoton.Eta(),bestPhoton.Phi(),(*Photons)[i].Eta(),(*Photons)[i].Phi() );
    //cout<<jentry<<" i:"<<i<<" phoSize:"<<Photons->size()<<" Pt:"<<bestPhoton.Pt()<<" eta:"<<bestPhoton.Eta()<<" phi:"<<bestPhoton.Phi()<<" otherP:"<<(*Photons)[i].Pt()<<" eta:"<<(*Photons)[i].Eta()<<" phi:"<<(*Photons)[i].Phi()<<" dR:"<<dR<<endl;
  }
  for(int i=0;i<Jets->size();i++){
    cout<<"JetPt:"<<(*Jets)[i].Pt()<<" JetEta:"<<(*Jets)[i].Eta()<<" JetPhi:"<<(*Jets)[i].Phi()<<endl;
  }
  cout<<"MHTPhi:"<<MHTPhi<<" DPhi1:"<<DeltaPhi1<<" DeltaPhi2:"<<DeltaPhi2<<" DeltaPhi3:"<<DeltaPhi3<<" DeltaPhi4:"<<DeltaPhi4<<endl;
}

bool SignalRegGraviton::passHEMjetVeto(double ptThresh) {
  for (int p = 0; p < Jets->size(); p++){
    if (-3.2 <= Jets->at(p).Eta() && Jets->at(p).Eta() <= -1.2 &&
	-1.77 <= Jets->at(p).Phi() && Jets->at(p).Phi() <= -0.67 &&
	Jets->at(p).Pt() > ptThresh && abs(DeltaPhi(Jets->at(p).Phi(),METPhi)) < 0.5)
      return false;
  }
  return true;
}

void SignalRegGraviton::changeJets(int jec2Use, int jer2Use, int jet2Vary){
  if((jec2Use*jer2Use)!=0) return;
  if(jet2Vary!=4 && jet2Vary!=8 && jet2Vary!=48 && jet2Vary!=84) return;
  
  if(jet2Vary==4 || jet2Vary==48 || jet2Vary==84){
    TLorentzVector iJet;
    vector<TLorentzVector> jets;
    if(jec2Use==-1){ BTagsDeepCSV = BTagsDeepCSVJECdown; MET = (*METDown)[1]; METPhi = (*METPhiDown)[1];}
    else if(jec2Use== 1){ BTagsDeepCSV = BTagsDeepCSVJECup; MET = (*METUp)[1]; METPhi = (*METPhiUp)[1];}
    else if(jer2Use==-1){ BTagsDeepCSV = BTagsDeepCSVJERdown; MET = (*METDown)[0]; METPhi = (*METPhiDown)[0];}
    else if(jer2Use== 1){ BTagsDeepCSV = BTagsDeepCSVJERup; MET = (*METUp)[0]; METPhi = (*METPhiUp)[0];}
  
    //---- looking at https://github.com/kpedro88/Analysis/blob/da0bb3e24e04768d3c205b45cbd74b18e638133c/KCode/KSkimmerVariators.h#L409-L436
    const auto& JetsUnc_origIndex = (jec2Use== 1)?*JetsJECup_origIndex:(jec2Use==-1)?*JetsJECdown_origIndex:(jer2Use==1)?*JetsJERup_origIndex:(jer2Use==-1)?*JetsJERdown_origIndex:*Jets_origIndex; //last one is a dummy value
    const auto& JetsUnc_jerFactor = (jec2Use== 1)?*JetsJECup_jerFactor:(jec2Use==-1)?*JetsJECdown_jerFactor:*Jets_jerFactor; //last one is a dummy value
    const auto& Jets_unc = (jec2Use== 1)?*Jets_jecUnc:(jec2Use==-1)?*Jets_jecUnc:(jer2Use==1)?*Jets_jerFactorUp:(jer2Use==-1)?*Jets_jerFactorDown:*Jets_jecFactor; //last one is a dummy value
    
    // if(EvtNum==207632 || EvtNum == 2055901){
    vector<int> newIndex(Jets_origIndex->size(),-1);
    for(unsigned k = 0; k < Jets_origIndex->size(); ++k){
      //reverse the index vector
      newIndex[(*Jets_origIndex)[k]] = k;
    }
    for(unsigned j = 0; j < JetsUnc_origIndex.size(); ++j){
      //Jets[Unc]_origIndex is sorted in the final order after uncertainty variation is applied
      //go up to common ancestor, then down to central smeared collection
      int i = newIndex[JetsUnc_origIndex[j]];
      //undo central smearing, apply JEC unc, redo smearing w/ new smearing factor
      if(jec2Use== 1)       jets.push_back((*Jets)[i]*(1./(*Jets_jerFactor)[i])*(1+Jets_unc[i])*JetsUnc_jerFactor[j]);
      else if(jec2Use== -1) jets.push_back((*Jets)[i]*(1./(*Jets_jerFactor)[i])*(1-Jets_unc[i])*JetsUnc_jerFactor[j]);
      else if(jer2Use==  1) jets.push_back((*Jets)[i]*(1./(*Jets_jerFactor)[i])*Jets_unc[i]);
      else if(jer2Use== -1) jets.push_back((*Jets)[i]*(1./(*Jets_jerFactor)[i])*Jets_unc[i]);
    }
    sortTLorVec(&jets);
    if(Jets->size()!=jets.size()){
      cout<<"oooooo.... Jets size:"<<Jets->size()<<" jets size:"<<jets.size()<<endl;
      return;
    }
    for(int i=0;i<jets.size();i++){
      //      if(EvtNum==207632 || EvtNum == 2055901){
      (*Jets)[i] = jets[i];
    }
  }//AK4 jets
  //---------------For AK8 jets
  if(jet2Vary==8 || jet2Vary==48 || jet2Vary==84){
    TLorentzVector iJet;
    vector<TLorentzVector> jets;
    const auto& JetsAK8Unc_origIndex = (jec2Use== 1)?*JetsAK8JECup_origIndex:(jec2Use==-1)?*JetsAK8JECdown_origIndex:(jer2Use==1)?*JetsAK8JERup_origIndex:(jer2Use==-1)?*JetsAK8JERdown_origIndex:*JetsAK8_origIndex; //last one is a dummy value
    const auto& JetsAK8Unc_jerFactor = (jec2Use== 1)?*JetsAK8JECup_jerFactor:(jec2Use==-1)?*JetsAK8JECdown_jerFactor:*JetsAK8_jerFactor; //last one is a dummy value
    const auto& JetsAK8_unc = (jec2Use== 1)?*JetsAK8_jecUnc:(jec2Use==-1)?*JetsAK8_jecUnc:(jer2Use==1)?*JetsAK8_jerFactorUp:(jer2Use==-1)?*JetsAK8_jerFactorDown:*JetsAK8_jecFactor; //last one is a dummy value
    
    // if(EvtNum==207632 || EvtNum == 2055901){
    vector<int> newIndexAK8(JetsAK8_origIndex->size(),-1);
    for(unsigned k = 0; k < JetsAK8_origIndex->size(); ++k){
      //reverse the index vector
      newIndexAK8[(*JetsAK8_origIndex)[k]] = k;
    }
    for(unsigned j = 0; j < JetsAK8Unc_origIndex.size(); ++j){
      //JetsAK8[Unc]_origIndex is sorted in the final order after uncertainty variation is applied
      //go up to common ancestor, then down to central smeared collection
      int i = newIndexAK8[JetsAK8Unc_origIndex[j]];
      //undo central smearing, apply JEC unc, redo smearing w/ new smearing factor
      if(jec2Use== 1)       jets.push_back((*JetsAK8)[i]*(1./(*JetsAK8_jerFactor)[i])*(1+JetsAK8_unc[i])*JetsAK8Unc_jerFactor[j]);
      else if(jec2Use== -1) jets.push_back((*JetsAK8)[i]*(1./(*JetsAK8_jerFactor)[i])*(1-JetsAK8_unc[i])*JetsAK8Unc_jerFactor[j]);
      else if(jer2Use==  1) jets.push_back((*JetsAK8)[i]*(1./(*JetsAK8_jerFactor)[i])*JetsAK8_unc[i]);
      else if(jer2Use== -1) jets.push_back((*JetsAK8)[i]*(1./(*JetsAK8_jerFactor)[i])*JetsAK8_unc[i]);
    }
    sortTLorVec(&jets);
    if(JetsAK8->size()!=jets.size()){
      cout<<"oooooo.... JetsAK8 size:"<<JetsAK8->size()<<" jets AK8 size:"<<jets.size()<<endl;
      return;
    }
    //get corrections for AK8 mass using Pt corrections applied.
    vector<double> ak8SDmass;
    for(int j=0;j<jets.size();j++){
      double mindr = 1000.;
      int mindrIdx = -1;
      for(int J=0;J<JetsAK8->size();J++){
	if(jets[j].DeltaR((*JetsAK8)[J]) < mindr){
	  mindr = jets[j].DeltaR((*JetsAK8)[J]);
	  mindrIdx = J;
	}
      }
      ak8SDmass.push_back(((*JetsAK8_softDropMass)[mindrIdx]) * (jets[j].Pt() / (*JetsAK8)[mindrIdx].Pt()));
    }
    for(int i=0;i<jets.size();i++){
      (*JetsAK8)[i] = jets[i];
      (*JetsAK8_softDropMass)[i] = ak8SDmass[i];
    }
  }
  //---------------------------
}


    //////////////////////////
    // cout<<"Before:"<<endl;
    // for(int J=0;J<JetsAK8->size();J++)
    //   cout<<"AK8 Pt: "<<(*JetsAK8)[J].Pt()<<" Eta: "<<(*JetsAK8)[J].Eta()<<" Phi: "<<" Mass: "<<(*JetsAK8)[J].M()<<" SD mass: "<<(*JetsAK8_softDropMass)[J]<<endl;
    /////////////////////////////

    // cout<<"After:"<<endl;
    // for(int j=0;j<JetsAK8->size();j++)
    //   cout<<"AK8 Pt: "<<(*JetsAK8)[j].Pt()<<" Eta: "<<(*JetsAK8)[j].Eta()<<" Phi: "<<" Mass: "<<(*JetsAK8)[j].M()<<" SD mass: "<<ak8SDmass[j]<<endl;
    // cout<<"**************"<<endl;
