#include "TTree.h"
#include "TNtuple.h"
#include "TEventList.h"
#include "TMath.h"
#include <iostream>
#include "TFile.h"
#include <cmath>

float EcalR=129.0;

void Rematch();

  unsigned int run, event; 
//  Int_t  HLTPhoton175, HLTPhoton160, HLTMET, HLTPFMET170, passGeneral;
  vector<int> *SubHits, *SatSubHits, *ClustEB, *ClustEE;
  vector<float> *TIso, *Ndof;
vector<float> *RZPar0, *RZPar1, *RZPar2, *RZErr2,  *XYPar0, *XYPar1, *XYPar2;
  vector<double> *DistEB, *seedFracEB, *EEB, *EtaEB, *PhiEB, *HIsoEB;
  vector<double> *DistEE, *seedFracEE, *EEE, *EtaEE, *PhiEE, *HIsoEE;
  UInt_t  HLTPhoton175, HLTPhoton160, HLTMET, HLTPFMET170, passGeneral;
void Trim(TTree *OldTree, string FileName, bool Refit=false, TEventList *Passed=NULL){
  TFile OutFile(FileName.c_str(), "recreate");
  TNtuple NewTree("monopoles", "monopoles", "run:event:Dist:SubHits:SatSubHits:seedFrac:E:TIso:HIso:XYPar0:XYPar1:XYPar2:RZPar0:RZPar1:RZPar2:RZErr2:Eta:Phi:HLTPhoton175:HLTPhoton160:HLTMET:HLTPFMET170:passGeneral");

  OldTree->SetBranchAddress("run", &run);
  OldTree->SetBranchAddress("event", &event);
  /*OldTree->SetBranchAddress("passtrig_HLTPhoton175",&HLTPhoton175);
  OldTree->SetBranchAddress("passtrig_HLTPhoton160",&HLTPhoton160);
  OldTree->SetBranchAddress("passtrig_HLTMET",&HLTMET);
  OldTree->SetBranchAddress("passtrig_HLTPFMET170",&HLTPFMET170);
  OldTree->SetBranchAddress("passGeneral",&passGeneral);
  */
  OldTree->SetBranchAddress("Track_clustMatchEB", &ClustEB);
  OldTree->SetBranchAddress("Track_clustDistEB", &DistEB);
  OldTree->SetBranchAddress("Track_clustMatchEE", &ClustEE);
  OldTree->SetBranchAddress("Track_clustDistEE", &DistEE);
  OldTree->SetBranchAddress("Track_Ndof", &Ndof);
  OldTree->SetBranchAddress("Track_SubHits", &SubHits);
  OldTree->SetBranchAddress("Track_SatSubHits", &SatSubHits);
  OldTree->SetBranchAddress("Track_Iso", &TIso);
  OldTree->SetBranchAddress("Track_XYPar0", &XYPar0);
  OldTree->SetBranchAddress("Track_XYPar1", &XYPar1);
  OldTree->SetBranchAddress("Track_XYPar2", &XYPar2);
  OldTree->SetBranchAddress("Track_RZPar0", &RZPar0);
  OldTree->SetBranchAddress("Track_RZPar1", &RZPar1);
  OldTree->SetBranchAddress("Track_RZPar2", &RZPar2);
  OldTree->SetBranchAddress("Track_RZErr2", &RZErr2);

  OldTree->SetBranchAddress("egComb_E", &EEB);
  OldTree->SetBranchAddress("egComb_frac51", &seedFracEB);
  OldTree->SetBranchAddress("egComb_eta", &EtaEB);
  OldTree->SetBranchAddress("egComb_phi", &PhiEB);
  OldTree->SetBranchAddress("egComb_hcalIso", &HIsoEB);

  OldTree->SetBranchAddress("eeComb_E", &EEE);
  OldTree->SetBranchAddress("eeComb_frac51", &seedFracEE);
  OldTree->SetBranchAddress("eeComb_eta", &EtaEE);
  OldTree->SetBranchAddress("eeComb_phi", &PhiEE);
  OldTree->SetBranchAddress("eeComb_hcalIso", &HIsoEE);
  OldTree->SetBranchAddress("passtrig_HLTPhoton175",&HLTPhoton175);
  OldTree->SetBranchAddress("passtrig_HLTPhoton160",&HLTPhoton160);
  OldTree->SetBranchAddress("passtrig_HLTMET",&HLTMET);
  OldTree->SetBranchAddress("passtrig_HLTPFMET170",&HLTPFMET170);
  OldTree->SetBranchAddress("passGeneral",&passGeneral);
  float x[23];

  int NEntries=OldTree->GetEntries();
  for(int entry=0; entry<NEntries; entry++){
    if(Passed != NULL && !Passed->Contains(entry)) continue;

//    cout << entry << endl;
    OldTree->GetEntry(entry);

    for(unsigned int track=0; track < TIso->size(); track++){
//      cout << track << '\t' << Clust->at(track) << '\t' << E->size() << endl;
      if(Ndof->at(track)<=3) continue;
      if(Refit) Rematch(); //Eta, Phi, XYPar0, XYPar1, XYPar2, RZPar0, RZPar1, RZPar2, Clust, Dist);
      if(ClustEB->at(track)==-1 && ClustEE->at(track)==-1) continue;

x[0]=run;
x[1]=event;
x[3]=SubHits->at(track); 
x[4]=SatSubHits->at(track); 
x[7]=TIso->at(track);
x[9]=XYPar0->at(track);
x[10]=XYPar1->at(track);
x[11]=XYPar2->at(track); 
x[12]=RZPar0->at(track);
x[13]=RZPar1->at(track);
x[14]=RZPar2->at(track);
x[15]=RZErr2->at(track);
x[18]=HLTPhoton175;
x[19]=HLTPhoton160;
x[20]=HLTMET;
x[21]=HLTPFMET170;
x[22]=passGeneral;



      if(DistEB->at(track) < DistEE->at(track)){
	x[2]=DistEB->at(track);x[5]=seedFracEB->at(ClustEB->at(track)), x[6]=EEB->at(ClustEB->at(track)); x[8]=HIsoEB->at(ClustEB->at(track)); x[16]=EtaEB->at(ClustEB->at(track)); x[17]=PhiEB->at(ClustEB->at(track));
      }else{
	x[2]=DistEE->at(track); x[5]=seedFracEE->at(ClustEE->at(track)), x[6]=EEE->at(ClustEE->at(track)); x[8]=HIsoEE->at(ClustEE->at(track)); x[16]=EtaEE->at(ClustEE->at(track)); x[17]=PhiEE->at(ClustEE->at(track));




      }
      NewTree.Fill(x);
    }
    if(entry%1000==0) cout << entry << " / " << NEntries << endl;
  }

  NewTree.Write();
  //OutFile.Write();
  //OutFile.Close();
}

void Rematch(){ //vector<double> *Eta, vector<double> *Phi, vector<float> *XYPar0, vector<float> *XYPar1, vector<float> *XYPar2, vector<float> *RZPar0, vector<float> *RZPar1, vector<float> *RZPar2, vector<int> *Clust, vector<double> *Dist){
  for(int i=0; i<XYPar0->size(); i++){

    //    if (E > 50){

    // Calculate expected Eta, Phi for ECAL cluster
    float ThisZ = RZPar0->at(i) + EcalR*RZPar1->at(i) + EcalR*EcalR*RZPar2->at(i);
    float ThisEta = TMath::ASinH(ThisZ / EcalR);
    float ThisPhi = XYPar1->at(i) - TMath::ASin((EcalR*EcalR-XYPar0->at(i)*XYPar2->at(i))/(2*EcalR*(XYPar2->at(i)-XYPar0->at(i))));

    float MinDR = 999;
    int BestEBCluster=-1, BestEECluster=-1;

    
    for(int j=0; j < EtaEB->size(); j++){


      float DEta = ThisEta-(*EtaEB)[j];
      float DPhi = ThisPhi-(*PhiEB)[j];
      while(DPhi < -3.14159) DPhi += 2*3.14159;
      while(DPhi >  3.14159) DPhi -= 2*3.14159;

      float ThisDR = sqrt(pow(DEta,2) + pow(DPhi,2));
      if(ThisDR < MinDR){
	MinDR = ThisDR;
	BestEBCluster = j;
      }
    }
    for(int j=0; j < EtaEE->size(); j++){
      float DEta = ThisEta-(*EtaEE)[j];
      float DPhi = ThisPhi-(*PhiEE)[j];
      while(DPhi < -3.14159) DPhi += 2*3.14159;
      while(DPhi >  3.14159) DPhi -= 2*3.14159;

      float ThisDR = sqrt(pow(DEta,2) + pow(DPhi,2));
      if(ThisDR < MinDR){
	MinDR = ThisDR;
	BestEECluster = j;
	BestEBCluster = -1;
      }
    }
    (*ClustEB)[i] = BestEBCluster;
    (*ClustEE)[i] = BestEECluster;
    (*DistEB)[i] = BestEBCluster >= 0 ? MinDR : 999;
    (*DistEE)[i] = BestEECluster >= 0 ? MinDR : 999;
  }
}

//}
