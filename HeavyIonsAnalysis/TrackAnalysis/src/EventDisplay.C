#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TAttLine.h"
#include "THelix.h"
#include "TView.h"
#include <string>
#include "TAttPad.h"
#include "TMath.h"
#include "TVector3.h"
#include "TView3D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"

#include <iostream>
#include <vector>

// Dataformat
//#include "include/smartJetName.h"
//#include "include/particleData.h"
//#include "include/jetData.h"
//#include "include/boostedEvtData.h"
//#include "include/eventData.h"

//nasty global variables for nextEvent() function
std::string currentFile;
int currentEvtIndx;
float currentPtCut;

float dR(float eta1, float phi1, float eta2, float phi2){
  return TMath::Power(TMath::Power(TMath::ACos(TMath::Cos(phi1-phi2)),2)+TMath::Power(eta1-eta2,2),0.5);
}

void formatHelix(THelix * h, float pz, float pt, int color = 0, bool doWTA = false){
  if(color==0) h->SetLineColor(kBlue);
  if(color==1) h->SetLineColor(kWhite);
  if(color==2) h->SetLineColor(kRed);
  if(color==3) h->SetLineColor(6);
  if(color==4) h->SetLineColor(7);
  if(color==5) h->SetLineColor(kGreen);
  h->SetLineWidth(1);
  
  float rangeBound = 1.0;
  if(!doWTA && pt<1.0 && TMath::Abs(pz)<0.5) rangeBound = 2*TMath::Abs(pz);
  //if(!doWTA && pt<1.5 && TMath::Abs(pz)<0.5) rangeBound = TMath::Abs(pz);
  //if(!doWTA && pt>1.5) rangeBound = rangeBound*pz/pt*10;
  //if(!doWTA && pt>5) rangeBound = 0.2;
  h->SetRange(0,rangeBound);
  if(pz<0) h->SetRange(-rangeBound,0);
  
  //if(pz>0)  h->SetRange(0,1, kHelixZ);
  //else      h->SetRange(-1,0,kHelixZ);
  //h->SetRange(-2,2,kLabX);
  //h->SetRange(-2,2,kLabY);

}

void EventDisplay(std::string inputFile = "/storage1/users/aab9/13TeV_ppGen_Nov24/Pythia/pthat800/output_MC_5.root" ,
                  int eventIndx = 595, 
		  float ptCut = 0, 
		  bool doWTA = false,
		  bool verbose = 0){
  currentFile=inputFile;
  currentEvtIndx=eventIndx;
  currentPtCut=ptCut;

  TFile * f 	= TFile::Open(inputFile.c_str(),"read");
  TTree * t 	= (TTree*)f->Get("analyzer/trackTree");


  // Setup branches
  std::vector< std::vector< float > > * pPt =  0;  
  std::vector< std::vector< float > > * pEta = 0;  
  std::vector< std::vector< float > > * pPhi = 0;  
  std::vector< std::vector< int > >   * pChg = 0;  

  t->SetBranchAddress("genDau_pt", &pPt);
  t->SetBranchAddress("genDau_eta",&pEta);
  t->SetBranchAddress("genDau_phi",&pPhi);
  t->SetBranchAddress("genDau_chg",&pChg);

  // Take the event of interest
  t->GetEntry(eventIndx);

  /* Draw Thrust and WTA axis
  TVector3 thrust = TVector3(0,0,0);
  thrust.SetMagThetaPhi(10,event.TTheta,event.TPhi);
  TVector3 wta1 = TVector3(0,0,0);
  if(WTAjet.nref>0) wta1.SetMagThetaPhi(10,2*TMath::ATan(TMath::Exp(-WTAjet.jteta[0])),WTAjet.jtphi[0]);
  TVector3 wta2 = TVector3(0,0,0);
  if(WTAjet.nref>1) wta2.SetMagThetaPhi(10,2*TMath::ATan(TMath::Exp(-WTAjet.jteta[1])),WTAjet.jtphi[1]);
  */

  std::cout << (pPt->at(0)).at(0) << " " << (pEta->at(0)).at(0) << " " << (pPhi->at(0)).at(0) << " " << (pChg->at(0)).at(0)  << std::endl;

  TCanvas *c = new TCanvas("c","Event Display",1000, 1000);
  c->Divide(2,2);
  c->SetFillColor(kGray+3);
  TCanvas* c1 = (TCanvas*)c->GetPad(1);//new TCanvas("AzimuthalView","AzimuthalView",600,600);
  c1->SetFillColor(kBlack);  
  TView *view = TView::CreateView(1);
  view->SetRange(-1,-1,-1,1,1,1);
  view->TopView(c1);

  TCanvas* c2 = (TCanvas*)c->GetPad(2); //new TCanvas("TopView","TopView",600,600);
  c2->SetFillColor(kBlack);  
  TView *view2 = TView::CreateView(1);
  view2->SetRange(-1,-1,-1,1,1,1);
  view2->SideView(c2);
  
  TCanvas* c3 = (TCanvas*)c->GetPad(3); //new TCanvas("AzimuthalThrustView","AzimuthalThrustView",600,600);
  c3->SetFillColor(kBlack);  
  TView *view3 = TView::CreateView(1);
  view3->SetRange(-1,-1,-1,1,1,1);
//  view3->RotateView(thrust.Phi()*180/TMath::Pi(),thrust.Theta()*180/TMath::Pi(),c3);
  view3->FrontView(c3);
  view3->ZoomView(c3,10);
  
  TCanvas* c4 = (TCanvas*)c->GetPad(4); //new TCanvas("AzimuthalThrustView","AzimuthalThrustView",600,600);
  c4->SetFillColor(kGray+3);  
  TView *view4 = TView::CreateView(1);
  view4->ZoomView(c4,1);
  c4->Update();
  
  TLegend *leg = new TLegend(0.1,0.1,0.9,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextColor(kWhite);
  


  TCanvas* cA = new TCanvas("cA","Azimuthal View",1000,1000);//new TCanvas("AzimuthalView","AzimuthalView",600,600);
  cA->SetFillColor(kBlack);  
  TView *viewA = TView::CreateView(1);
  viewA->SetRange(-1,-1,-1,1,1,1);
  viewA->TopView(c1);

  THelix * helix[1000];

  /*
  if(!doWTA){
    helix[999] = new THelix(0,0,0,thrust.Px(),thrust.Py(),thrust.Pz(),0.00000001);
    helix[999]->SetRange(-2,2);
    helix[999]->SetLineColor(8);
    helix[999]->SetLineWidth(3);
    c1->cd();
    helix[999]->Draw();
    c2->cd();
    helix[999]->Draw();
    c3->cd();
    helix[999]->Draw();
    
    helix[998] = new THelix(0,0,0,wta1.Px(),wta1.Py(),wta1.Pz(),0.00000001);
    if(wta1.Pz()<0) helix[998]->SetRange(-1,0);
    if(wta1.Pz()>=0) helix[998]->SetRange(0,1);
    helix[998]->SetLineColor(38);
    helix[998]->SetLineWidth(3);
    c1->cd();
    helix[998]->Draw();
    c2->cd();
    helix[998]->Draw();
    c3->cd();
    helix[998]->Draw();
    
    helix[997] = new THelix(0,0,0,wta2.Px(),wta2.Py(),wta2.Pz(),0.000001);
    if(wta2.Pz()<0) helix[997]->SetRange(-1,0);
    if(wta2.Pz()>=0) helix[997]->SetRange(0,1);
    helix[997]->SetLineColor(38);
    helix[997]->SetLineWidth(3);
    c1->cd();
    helix[997]->Draw();
    c2->cd();
    helix[997]->Draw();
    c3->cd();
    helix[997]->Draw();
    leg->AddEntry(helix[999],"ALEPH Archived Data Event Display","t");
    leg->AddEntry(helix[999],"Thrust Axis","l");
    leg->AddEntry(helix[998],"WTA Jet Axis 1","l");
    leg->AddEntry(helix[997],"WTA Jet Axis 2","l");
    
    helix[992] = new THelix(0,0,0,1,1,1,0);
    helix[993] = new THelix(0,0,0,1,1,1,0);
    helix[994] = new THelix(0,0,0,1,1,1,0);
    helix[995] = new THelix(0,0,0,1,1,1,0);
    helix[996] = new THelix(0,0,0,1,1,1,0);
    formatHelix(helix[993],1,1,1,0);
    leg->AddEntry(helix[993],"Tracks in Leading Jet (#Delta R<0.8)","l");
    formatHelix(helix[994],1,1,2,0);
    leg->AddEntry(helix[994],"Tracks in Subleading Jet (#Delta R<0.8)","l");
    formatHelix(helix[995],1,1,3,0);
    leg->AddEntry(helix[995],"Tracks in Third Jet (#Delta R<0.8)","l");
    formatHelix(helix[996],1,1,4,0);
    leg->AddEntry(helix[996],"Tracks in Fourth Jet (#Delta R<0.8)","l");
    formatHelix(helix[992],1,1,0,0);
    leg->AddEntry(helix[992],"Other Tracks","l");
       
    
  }
  */


  // Draw charged particles
  int nHelix = 0;
  for(int j = 0; j<(pPt->size()); j++){
 // std::cout << pPt->size() << std::endl;
  for(int i = 0; i<(pPt->at(j).size()); i++){
 //   std::cout << pPt->at(j).size() << std::endl;
    if((pChg->at(j)).at(i) == 0) continue;
    if((pPt->at(j)).at(i) < ptCut) continue;

    //jet tracks
    int trackColor = 2;
    if(j>0) trackColor = 4;
    //if(jet.jtpt[0]*TMath::CosH(jet.jteta[0])>0 && dR(jet.jteta[0],jet.jtphi[0],particle.eta[i],particle.phi[i])<0.8) trackColor = 1;
    //else if(jet.jtpt[1]*TMath::CosH(jet.jteta[1])>0 && dR(jet.jteta[1],jet.jtphi[1],particle.eta[i],particle.phi[i])<0.8) trackColor = 2;
    //else if(jet.jtpt[2]*TMath::CosH(jet.jteta[2])>0 && dR(jet.jteta[2],jet.jtphi[2],particle.eta[i],particle.phi[i])<0.8) trackColor = 3;
    //else if(jet.jtpt[3]*TMath::CosH(jet.jteta[3])>0 && dR(jet.jteta[3],jet.jtphi[3],particle.eta[i],particle.phi[i])<0.8) trackColor = 4;
    //if(doWTA && boosted.pt[i]<0.01) trackColor = 5;
    //if (verbose) std::cout <<  jet.jteta[0] << " " << jet.jtphi[0] <<" " <<  particle.eta[i] << " " << particle.phi[i] << std::endl;
    //if (verbose) std::cout <<  dR(jet.jteta[0],jet.jtphi[0],particle.eta[i],particle.phi[i]) << std::endl;

    Float_t px,py,pz,pt;
    float phi = (pPhi->at(j)).at(i);
    pt = (pPt->at(j)).at(i);
    px = pt*TMath::Cos(phi);
    py = pt*TMath::Sin(phi);
    pz = pt*TMath::SinH( (pEta->at(j)).at(i) );

    

   // if(i>15) continue;
    std::cout << j << " "<<  i<<" " << pt << " " << (pEta->at(j)).at(i) << " " << (pPhi->at(j)).at(i) << " " << px<<" "  << py<<" " << pz <<  std::endl; 

    /*
    if(i%4==0) helix[nHelix] = new THelix(0,0,0,10,10,5,1 * 1.5);
    if(i%4==1) helix[nHelix] = new THelix(0,0,0,10,10,-5,1 * 1.5);
    if(i%4==2) helix[nHelix] = new THelix(0,0,0,10,10,5,-1 * 1.5);
    if(i%4==3) helix[nHelix] = new THelix(0,0,0,10,10,-5,-1 * 1.5);

    if(i%4==0 || i%4==2) formatHelix(helix[nHelix],5,14.14,trackColor,doWTA);
    if(i%4==1 || i%4==3) formatHelix(helix[nHelix],-5,14.14,trackColor,doWTA);
    */ 

    helix[nHelix] = new THelix(0,0,0,px,py,pz,(pChg->at(j)).at(i) * 0.2);
    formatHelix(helix[nHelix],pz,pt,trackColor,doWTA);
    c1->cd();
    helix[nHelix]->Draw("same");
    c2->cd();
    helix[nHelix]->Draw("same");
    c3->cd();
    helix[nHelix]->Draw("same");
    cA->cd();
    helix[nHelix]->Draw("same");
    nHelix++;
  }
  }

  c1->cd();
  TLatex *l1 = new TLatex(-0.97,0.90,"Azimuthal View");
  l1->SetTextColor(kWhite);
  l1->Draw();
  c2->cd();
  TLatex *l2 = new TLatex(-0.97,0.90,"Side View");
  l2->SetTextColor(kWhite);
  l2->Draw();
  c3->cd();
//  TLatex *l3 = new TLatex(-0.97,0.90,"Azimuthal Thrust View");
  TLatex *l3 = new TLatex(-0.97,0.90,"Front View");
  l3->SetTextColor(kWhite);
  l3->Draw();

  c4->cd();
  leg->Draw();

  cA->cd();
//  TLatex *l3 = new TLatex(-0.97,0.90,"Azimuthal Thrust View");
  TLatex *lA1 = new TLatex(0.3,0.35,"p_{T} = 939 GeV");
  lA1->SetTextColor(kRed+1);
  lA1->Draw("same");
  TLatex *lA2 = new TLatex(0.3,0.2,"N_{ch} = 27");
  lA2->SetTextColor(kRed+1);
  lA2->Draw("same");
  
  TLatex *lA3 = new TLatex(-0.9,-0.5,"p_{T} = 880 GeV");
  lA3->SetTextColor(7);
  lA3->Draw("same");
  TLatex *lA4 = new TLatex(-0.9,-0.65,"N_{ch} = 112");
  lA4->SetTextColor(7);
  lA4->Draw("same");
  

  TLatex *lA5 = new TLatex(0.4,0.9,"PYTHIA 8");
  lA5->SetTextColor(kWhite);
  lA5->Draw("same");
  TLatex *lA6 = new TLatex(0.4,0.78,"13 TeV pp");
  lA6->SetTextColor(kWhite);
  lA6->Draw("same");
  TLatex *lA7 = new TLatex(0.4,0.66,"anti-k_{t} R=0.8");
  lA7->SetTextColor(kWhite);
  lA7->Draw("same");
  cA->SaveAs("../plots/EventDisplay.pdf");
  cA->SaveAs("../plots/EventDisplay.png");

}

int nextEvent(){
  EventDisplay(currentFile,currentEvtIndx+1,currentPtCut);
  return currentEvtIndx;
}
