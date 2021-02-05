#include <string>
#include "include/Timer.h"
#include <TStyle.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TProfile.h"
#include "TGraph.h"
#include <vector>
#include "math.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TCutG.h"
#include "TCanvas.h"
#include <numeric>
#include <iomanip>
#include <vector>
#include <iostream>
#include <fstream>
#include "include/coordinateTools.h"

void analyze( std::vector< std::string> files){

  const int nqHats = 8;
  const float qHatBoundaries[nqHats+1] = {470, 600, 800, 1000, 1400, 1800, 2400, 3200, 13000};
  //units in pb
  float xs[nqHats] = {552.1, 156.5, 26.28, 7.47, 0.6484, 0.08743, 0.005236, 0.0001357};
  TH1D * qHatHist = new TH1D("qhatHit",";;qHat",nqHats,qHatBoundaries);

  float genQScale = 0;
  std::vector< float > * genJetPt = 0;
  std::vector< float > * genJetEta = 0;
  std::vector< float > * genJetPhi = 0;
  std::vector< int > * genJetChargedMultiplicity = 0;
  std::vector< std::vector< int > > * genDau_chg = 0;
  std::vector< std::vector< int > > * genDau_pid = 0;
  std::vector< std::vector< float > > * genDau_pt = 0;
  std::vector< std::vector< float > > * genDau_eta = 0;
  std::vector< std::vector< float > > * genDau_phi = 0;


  std::cout << "calculating weighting factors for qhat bins" << std::endl;
  for(unsigned int f = 0; f<files.size(); f++){
    //for testing
    //if(f%10!=0) continue;
    //if(f%10==0) std::cout << f << "/" << files.size() << std::endl;
    std::cout << f << "/" << files.size() << std::endl;

    TFile * inputFile = TFile::Open(files.at(f).c_str(),"read");
    TTree * t = (TTree*)inputFile->Get("analyzer/trackTree");

    t->SetBranchAddress("genQScale",&genQScale);
    for( int i = 0; i<t->GetEntries(); i++){
      t->GetEntry(i);

      qHatHist->Fill(genQScale);    
    }
  
    inputFile->Close();
  } 
  qHatHist->Print("All");
  
  const int nMultBins = 5;
  int multLow[nMultBins] =  {0,   0 , 30 , 50 , 70};
  int multHigh[nMultBins] = {999, 30 ,50 , 70,  999};

  //analysis loop
  TFile * output = TFile::Open("PythiaOutput.root","recreate");
  TH1D * pthat = new TH1D("pthat",";#hat{q};#sigma (pb)",303,470,3500);
  TH1D * leadingJetPt = new TH1D("leadingJetPt",";Leading p_{T}^{gen};#sigma (pb)",50,500,1500);
  TH1D * jetPt = new TH1D("JetPt",";p_{T}^{gen};#sigma (pb)",50,100,1500);
  TH1D * genJetChargedMultiplicity_h = new TH1D("genJetChargedMultiplicity",";gen Charged Multiplicity;#sigma (pb)",70,0,140);
  TH2D * multVsPt = new TH2D("multVsPt",";p_{T}^{gen};genChargedMultiplicity",50,500,1500,50,0,100);
  TH1D * mult = new TH1D("mult",";n_{ch};Arbitrary Units",70,0,140);
  TH1D * ptStar_h = new TH1D("ptStar",";j_{T};#sigma (pb)",50,0,10);
  TH1D * ptStar_0_h = new TH1D("ptStar_0",";j_{T};#sigma (pb)",50,0,10);
  TH1D * ptStar_30_h = new TH1D("ptStar_30",";j_{T};#sigma (pb)",50,0,10);
  TH1D * ptStar_50_h = new TH1D("ptStar_50",";j_{T};#sigma (pb)",50,0,10);
  TH1D * ptStar_70_h = new TH1D("ptStar_70",";j_{T};#sigma (pb)",50,0,10);
  TH1D * ptStar_ks_h = new TH1D("ptStar_ks",";j_{T};#sigma (pb)",50,0,10);
  TH1D * ptStar_phi_h = new TH1D("ptStar_phi",";j_{T};#sigma (pb)",50,0,10);
  TH1D * ptStar_lambda_h = new TH1D("ptStar_lambda",";j_{T};#sigma (pb)",50,0,10);
  TH1D * ptStar_pi_h[nMultBins];
  TH1D * ptStar_k_h[nMultBins];
  TH1D * ptStar_p_h[nMultBins];
  for(int m = 0; m<nMultBins; m++){
    ptStar_pi_h[m] = new TH1D(Form("ptStar_pi_%d_%d",multLow[m],multHigh[m] ),";j_{T};#sigma (pb)",100,0,20);
    ptStar_k_h[m] = new  TH1D(Form("ptStar_p_%d_%d",multLow[m],multHigh[m] ),";j_{T};#sigma (pb)",100,0,20);
    ptStar_p_h[m] = new  TH1D(Form("ptStar_k_%d_%d",multLow[m],multHigh[m] ),";j_{T};#sigma (pb)",100,0,20);
  }
  TH1D * mtStar_pi_h[nMultBins];
  TH1D * mtStar_k_h[nMultBins];
  TH1D * mtStar_p_h[nMultBins];
  for(int m = 0; m<nMultBins; m++){
    mtStar_pi_h[m] = new TH1D(Form("mtStar_pi_%d_%d",multLow[m],multHigh[m] ),";m_{T};#sigma (pb)",100,0,20);
    mtStar_k_h[m] = new  TH1D(Form("mtStar_p_%d_%d",multLow[m],multHigh[m] ),";m_{T};#sigma (pb)",100,0,20);
    mtStar_p_h[m] = new  TH1D(Form("mtStar_k_%d_%d",multLow[m],multHigh[m] ),";m_{T};#sigma (pb)",100,0,20);
  }

  float etaStarJetWeightSum = 0;
  TH1D * etaStar_h = new TH1D("etaStar",";#eta^{*};#sigma (pb)",50,0,10);
  TH1D * etaStar_0_h = new TH1D("etaStar_0",";#eta^{*};#sigma (pb)",50,0,10);
  TH1D * etaStar_30_h = new TH1D("etaStar_30",";#eta^{*};#sigma (pb)",50,0,10);
  TH1D * etaStar_50_h = new TH1D("etaStar_50",";#eta^{*};#sigma (pb)",50,0,10);
  TH1D * etaStar_70_h = new TH1D("etaStar_70",";#eta^{*};#sigma (pb)",50,0,10);
  TH1D * etaStar_p_h = new TH1D("etaStar_p",";#eta^{*};#sigma (pb)",50,0,10);
  TH1D * etaStar_k_h = new TH1D("etaStar_k",";#eta^{*};#sigma (pb)",50,0,10);
  TH1D * etaStar_pi_h = new TH1D("etaStar_pi",";#eta^{*};#sigma (pb)",50,0,10);
  TH1D * phiStar_h = new TH1D("phiStar",";#phi^{*};#sigma (pb)",50,-TMath::Pi(),TMath::Pi());
 
  const int jtMultBins = 4;
  float multJtBins[jtMultBins+1] =  {0,   30 , 50 , 70 , 90};
  TH1D * avgJt_pi = new TH1D("avgJt_pi",";mult;<j_{T}>", jtMultBins, multJtBins);
  TH1D * avgJt_p = new TH1D("avgJt_p",";mult;<j_{T}>", jtMultBins, multJtBins);
  TH1D * avgJt_k = new TH1D("avgJt_k",";mult;<j_{T}>", jtMultBins, multJtBins);
  TH1D * avgMt_pi = new TH1D("avgMt_pi",";mult;<m_{T}>", jtMultBins, multJtBins);
  TH1D * avgMt_p = new TH1D("avgMt_p",";mult;<m_{T}>", jtMultBins, multJtBins);
  TH1D * avgMt_k = new TH1D("avgMt_k",";mult;<m_{T}>", jtMultBins, multJtBins);

  for(unsigned int f = 0; f<files.size(); f++){
    //for testing
    //if(f%10!=0) continue;

    //if(f%10==0) std::cout << f << "/" << files.size() << std::endl;
    std::cout << f << "/" << files.size() << std::endl;

    TFile * inputFile = TFile::Open(files.at(f).c_str(),"read");
    TTree * t = (TTree*)inputFile->Get("analyzer/trackTree");

    t->SetBranchAddress("genQScale",&genQScale);
    t->SetBranchAddress("genJetPt",&genJetPt);
    t->SetBranchAddress("genJetEta",&genJetEta);
    t->SetBranchAddress("genJetPhi",&genJetPhi);
    t->SetBranchAddress("genJetChargedMultiplicity",&genJetChargedMultiplicity); 
    t->SetBranchAddress("genDau_chg",&genDau_chg); 
    t->SetBranchAddress("genDau_pid",&genDau_pid); 
    t->SetBranchAddress("genDau_pt",&genDau_pt); 
    t->SetBranchAddress("genDau_eta",&genDau_eta); 
    t->SetBranchAddress("genDau_phi",&genDau_phi); 
    //event loop
    for(int i = 0; i<t->GetEntries(); i++){
      t->GetEntry(i);
      //weight by xsection/total number of gen events in the pthat bin
      pthat->Fill(genQScale, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)) );    

      //jet spectra
      if(genJetPt->size()==0) continue;
      if(genJetChargedMultiplicity->size()==0) continue;
      leadingJetPt->Fill(genJetPt->at(0), xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)) );
      for(int j = 0; j<genJetPt->size(); j++){
        jetPt->Fill(genJetPt->at(j), xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)) );

        if(genJetPt->at(j) > 500){
          etaStarJetWeightSum += xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale));
          genJetChargedMultiplicity_h->Fill(genJetChargedMultiplicity->at(j), xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)) );
          mult->Fill(genJetChargedMultiplicity->at(j), xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)) );
          multVsPt->Fill(genJetPt->at(j), genJetChargedMultiplicity->at(j), xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)) );

          for(int k = 0; k<(genDau_chg->at(j)).size(); k++){
            //exception for lambdas, phi, kshort
            if( (genDau_chg->at(j)).size() == 0 ) continue;
            if( (genDau_chg->at(j)).at(k) == 0 ) continue;
            float ptStar = ptWRTJet(genJetPt->at(j), genJetEta->at(j), genJetPhi->at(j), (genDau_pt->at(j)).at(k), (genDau_eta->at(j)).at(k), (genDau_phi->at(j)).at(k));
            float etaStar = etaWRTJet(genJetPt->at(j), genJetEta->at(j), genJetPhi->at(j), (genDau_pt->at(j)).at(k), (genDau_eta->at(j)).at(k), (genDau_phi->at(j)).at(k));
            float phiStar = phiWRTJet(genJetPt->at(j), genJetEta->at(j), genJetPhi->at(j), (genDau_pt->at(j)).at(k), (genDau_eta->at(j)).at(k), (genDau_phi->at(j)).at(k));
            float mtStar = 0;
            if( TMath::Abs((genDau_pid->at(j)).at(k)) == 211 ) mtStar = TMath::Sqrt(ptStar*ptStar+0.1395*0.1395) - 0.1395;  
            if( TMath::Abs((genDau_pid->at(j)).at(k)) == 2212) mtStar = TMath::Sqrt(ptStar*ptStar+0.938*0.938) - 0.938;
            if( TMath::Abs((genDau_pid->at(j)).at(k)) == 321 ) mtStar = TMath::Sqrt(ptStar*ptStar+0.494*0.494) - 0.494;
 
            ptStar_h->Fill(ptStar, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)));
            if( genJetChargedMultiplicity->at(j) < 30) ptStar_0_h->Fill(ptStar, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)));
            if( genJetChargedMultiplicity->at(j) >= 30 && genJetChargedMultiplicity->at(j) <50) ptStar_30_h->Fill(ptStar, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)));
            if( genJetChargedMultiplicity->at(j) >= 50 && genJetChargedMultiplicity->at(j) < 70) ptStar_50_h->Fill(ptStar, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)));
            if( genJetChargedMultiplicity->at(j) >= 70) ptStar_70_h->Fill(ptStar, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)));
            if( TMath::Abs((genDau_pid->at(j)).at(k)) == 211 ){
              for(int m = 0; m<nMultBins; m++){
                if(genJetChargedMultiplicity->at(j)>= multLow[m] && genJetChargedMultiplicity->at(j)<multHigh[m]) ptStar_pi_h[m]->Fill(ptStar, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)));
                if(genJetChargedMultiplicity->at(j)>= multLow[m] && genJetChargedMultiplicity->at(j)<multHigh[m]) mtStar_pi_h[m]->Fill(mtStar, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)));
              }
            }
            if( TMath::Abs((genDau_pid->at(j)).at(k)) == 2212 ){
              for(int m = 0; m<nMultBins; m++){
                if(genJetChargedMultiplicity->at(j)>= multLow[m] && genJetChargedMultiplicity->at(j)<multHigh[m]) ptStar_p_h[m]->Fill(ptStar, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)));
                if(genJetChargedMultiplicity->at(j)>= multLow[m] && genJetChargedMultiplicity->at(j)<multHigh[m]) mtStar_p_h[m]->Fill(mtStar, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)));
              }
            }
            if( TMath::Abs((genDau_pid->at(j)).at(k)) == 321 ){
              for(int m = 0; m<nMultBins; m++){
                if(genJetChargedMultiplicity->at(j)>= multLow[m] && genJetChargedMultiplicity->at(j)<multHigh[m]) ptStar_k_h[m]->Fill(ptStar, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)));
                if(genJetChargedMultiplicity->at(j)>= multLow[m] && genJetChargedMultiplicity->at(j)<multHigh[m]) mtStar_k_h[m]->Fill(mtStar, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)));
              }
            }
            etaStar_h->Fill(etaStar, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)));
            if( genJetChargedMultiplicity->at(j) < 30) etaStar_0_h->Fill(etaStar, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)));
            if( genJetChargedMultiplicity->at(j) >= 30 && genJetChargedMultiplicity->at(j) < 50) etaStar_30_h->Fill(etaStar, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)));
            if( genJetChargedMultiplicity->at(j) >= 50 && genJetChargedMultiplicity->at(j) < 70) etaStar_50_h->Fill(etaStar, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)));
            if( genJetChargedMultiplicity->at(j) >= 70) etaStar_70_h->Fill(etaStar, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)));
            if( TMath::Abs((genDau_pid->at(j)).at(k)) == 211 ) etaStar_pi_h->Fill(etaStar, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)));
            if( TMath::Abs((genDau_pid->at(j)).at(k)) == 2212 ) etaStar_p_h->Fill(etaStar, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)));
            if( TMath::Abs((genDau_pid->at(j)).at(k)) == 321 ) etaStar_k_h->Fill(etaStar, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)));
            phiStar_h->Fill(phiStar, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)));

          }
        }
      }
    }


    inputFile->Close();
  } 

  pthat->Print("All");
  

  gStyle->SetOptStat(0);
  TCanvas * c1 = new TCanvas("c1","c1",600.600);
  c1->SetLeftMargin(0.2);
  pthat->Draw("p");
  c1->SetLogy();
  c1->SaveAs("plots/pthat.png");
  
  jetPt->Draw("p");
  leadingJetPt->SetMarkerColor(kRed);
  leadingJetPt->SetLineColor(kRed);
  leadingJetPt->Draw("same p");
  c1->SaveAs("plots/jetSpectra.png");

  genJetChargedMultiplicity_h->Draw("p");
  c1->SaveAs("plots/genJetChargedMult.png");

  ptStar_h->Draw("p");
  c1->SaveAs("plots/ptStar.png");
  etaStar_h->Draw("p");
  c1->SaveAs("plots/etaStar.png");
  phiStar_h->Draw("p");
  c1->SetLogy(0);
  c1->SaveAs("plots/phiStar.png");
  
  output->cd();
  TH1D * etaStar_h_dNdEta = (TH1D*) etaStar_h->Clone("etaStar_h_dNdEta");
  etaStar_h_dNdEta->SetDirectory(output);
  etaStar_h_dNdEta->Scale(1.0/etaStarJetWeightSum);

  c1->SetLogy();
  float ptStarEntries = ptStar_h->Integral();
  float etaStarEntries = etaStar_h->Integral();
  float ptStarEntries_0 = ptStar_0_h->Integral();
  float ptStarEntries_30 = ptStar_30_h->Integral();
  float ptStarEntries_50 = ptStar_50_h->Integral();
  float ptStarEntries_70 = ptStar_70_h->Integral();
  float ptStarEntries_p[nMultBins];
  float ptStarEntries_k[nMultBins];
  float ptStarEntries_pi[nMultBins];
  for(int m = 0; m<nMultBins; m++){
    ptStarEntries_p[m] = ptStar_p_h[m]->Integral();
    ptStarEntries_k[m] = ptStar_k_h[m]->Integral();
    ptStarEntries_pi[m] = ptStar_pi_h[m]->Integral();
  }
  float etaStarEntries_0 = etaStar_0_h->Integral();
  float etaStarEntries_30 = etaStar_30_h->Integral();
  float etaStarEntries_50 = etaStar_50_h->Integral();
  float etaStarEntries_70 = etaStar_70_h->Integral();
  float etaStarEntries_p = etaStar_p_h->Integral();
  float etaStarEntries_k = etaStar_k_h->Integral();
  float etaStarEntries_pi = etaStar_pi_h->Integral();

  ptStar_h->Scale(1.0/ptStarEntries); 
  etaStar_h->Scale(1.0/etaStarEntries); 
  ptStar_0_h->Scale(1.0/ptStarEntries_0); 
  ptStar_30_h->Scale(1.0/ptStarEntries_30); 
  ptStar_50_h->Scale(1.0/ptStarEntries_50); 
  ptStar_70_h->Scale(1.0/ptStarEntries_70); 

  for(int m = 0; m<nMultBins; m++){
    ptStar_p_h[m]->Scale(1.0/ptStarEntries_p[m]); 
    ptStar_k_h[m]->Scale(1.0/ptStarEntries_k[m]); 
    ptStar_pi_h[m]->Scale(1.0/ptStarEntries_pi[m]); 
  }
  etaStar_0_h->Scale(1.0/etaStarEntries_0); 
  etaStar_30_h->Scale(1.0/etaStarEntries_30); 
  etaStar_50_h->Scale(1.0/etaStarEntries_50); 
  etaStar_70_h->Scale(1.0/etaStarEntries_70); 
  etaStar_p_h->Scale(1.0/etaStarEntries_p); 
  etaStar_k_h->Scale(1.0/etaStarEntries_k); 
  etaStar_pi_h->Scale(1.0/etaStarEntries_pi); 

  ptStar_h->Draw("p");
  ptStar_30_h->SetLineColor(kRed);
  ptStar_50_h->SetLineColor(kBlue);
  ptStar_70_h->SetLineColor(kGreen);
  ptStar_30_h->Draw("same p");
  ptStar_50_h->Draw("same p");
  ptStar_70_h->Draw("same p");
  c1->SaveAs("plots/jtVsMult.png");

  ptStar_h->SetLineColor(kBlack);
  ptStar_0_h->SetLineColor(kBlack);
  ptStar_30_h->SetLineColor(kBlack);
  ptStar_50_h->SetLineColor(kBlack);
  ptStar_70_h->SetLineColor(kBlack);
  for(int m = 0; m<nMultBins; m++){
    if(m==0) ptStar_h->Draw("p");
    if(m==1) ptStar_0_h->Draw("p");
    if(m==2) ptStar_30_h->Draw("p");
    if(m==3) ptStar_50_h->Draw("p");
    if(m==4) ptStar_70_h->Draw("p");
    ptStar_p_h[m]->SetLineColor(kRed);
    ptStar_k_h[m]->SetLineColor(kBlue);
    ptStar_pi_h[m]->SetLineColor(kGreen);
    ptStar_p_h[m]->Draw("same p");
    ptStar_k_h[m]->Draw("same p");
    ptStar_pi_h[m]->Draw("same p");
    c1->SaveAs(Form("plots/jtVspid_%d_%d.png",multLow[m], multHigh[m]));
  }
  
  avgJt_pi->SetBinContent(1,ptStar_pi_h[1]->GetMean()); 
  avgJt_pi->SetBinContent(2,ptStar_pi_h[2]->GetMean()); 
  avgJt_pi->SetBinContent(3,ptStar_pi_h[3]->GetMean()); 
  avgJt_pi->SetBinContent(4,ptStar_pi_h[4]->GetMean()); 
  avgJt_pi->SetBinError(1,ptStar_pi_h[1]->GetMeanError()); 
  avgJt_pi->SetBinError(2,ptStar_pi_h[2]->GetMeanError()); 
  avgJt_pi->SetBinError(3,ptStar_pi_h[3]->GetMeanError()); 
  avgJt_pi->SetBinError(4,ptStar_pi_h[4]->GetMeanError()); 
  avgJt_p->SetBinContent(1,ptStar_p_h[1]->GetMean()); 
  avgJt_p->SetBinContent(2,ptStar_p_h[2]->GetMean()); 
  avgJt_p->SetBinContent(3,ptStar_p_h[3]->GetMean()); 
  avgJt_p->SetBinContent(4,ptStar_p_h[4]->GetMean()); 
  avgJt_p->SetBinError(1,ptStar_p_h[1]->GetMeanError()); 
  avgJt_p->SetBinError(2,ptStar_p_h[2]->GetMeanError()); 
  avgJt_p->SetBinError(3,ptStar_p_h[3]->GetMeanError()); 
  avgJt_p->SetBinError(4,ptStar_p_h[4]->GetMeanError()); 
  avgJt_k->SetBinContent(1,ptStar_k_h[1]->GetMean()); 
  avgJt_k->SetBinContent(2,ptStar_k_h[2]->GetMean()); 
  avgJt_k->SetBinContent(3,ptStar_k_h[3]->GetMean()); 
  avgJt_k->SetBinContent(4,ptStar_k_h[4]->GetMean()); 
  avgJt_k->SetBinError(1,ptStar_k_h[1]->GetMeanError()); 
  avgJt_k->SetBinError(2,ptStar_k_h[2]->GetMeanError()); 
  avgJt_k->SetBinError(3,ptStar_k_h[3]->GetMeanError()); 
  avgJt_k->SetBinError(4,ptStar_k_h[4]->GetMeanError()); 
  avgJt_k->SetMarkerColor(kBlue); 
  avgJt_pi->SetMarkerColor(kGreen); 
  avgJt_p->SetMarkerColor(kRed); 
  avgJt_k->SetLineColor(kBlue); 
  avgJt_pi->SetLineColor(kGreen); 
  avgJt_p->SetLineColor(kRed); 

 
  avgMt_pi->SetBinContent(1,mtStar_pi_h[1]->GetMean()); 
  avgMt_pi->SetBinContent(2,mtStar_pi_h[2]->GetMean()); 
  avgMt_pi->SetBinContent(3,mtStar_pi_h[3]->GetMean()); 
  avgMt_pi->SetBinContent(4,mtStar_pi_h[4]->GetMean()); 
  avgMt_pi->SetBinError(1,mtStar_pi_h[1]->GetMeanError()); 
  avgMt_pi->SetBinError(2,mtStar_pi_h[2]->GetMeanError()); 
  avgMt_pi->SetBinError(3,mtStar_pi_h[3]->GetMeanError()); 
  avgMt_pi->SetBinError(4,mtStar_pi_h[4]->GetMeanError()); 
  avgMt_p->SetBinContent(1,mtStar_p_h[1]->GetMean()); 
  avgMt_p->SetBinContent(2,mtStar_p_h[2]->GetMean()); 
  avgMt_p->SetBinContent(3,mtStar_p_h[3]->GetMean()); 
  avgMt_p->SetBinContent(4,mtStar_p_h[4]->GetMean()); 
  avgMt_p->SetBinError(1,mtStar_p_h[1]->GetMeanError()); 
  avgMt_p->SetBinError(2,mtStar_p_h[2]->GetMeanError()); 
  avgMt_p->SetBinError(3,mtStar_p_h[3]->GetMeanError()); 
  avgMt_p->SetBinError(4,mtStar_p_h[4]->GetMeanError()); 
  avgMt_k->SetBinContent(1,mtStar_k_h[1]->GetMean()); 
  avgMt_k->SetBinContent(2,mtStar_k_h[2]->GetMean()); 
  avgMt_k->SetBinContent(3,mtStar_k_h[3]->GetMean()); 
  avgMt_k->SetBinContent(4,mtStar_k_h[4]->GetMean()); 
  avgMt_k->SetBinError(1,mtStar_k_h[1]->GetMeanError()); 
  avgMt_k->SetBinError(2,mtStar_k_h[2]->GetMeanError()); 
  avgMt_k->SetBinError(3,mtStar_k_h[3]->GetMeanError()); 
  avgMt_k->SetBinError(4,mtStar_k_h[4]->GetMeanError()); 
  avgMt_k->SetMarkerColor(kBlue); 
  avgMt_pi->SetMarkerColor(kGreen); 
  avgMt_p->SetMarkerColor(kRed); 
  avgMt_k->SetLineColor(kBlue); 
  avgMt_pi->SetLineColor(kGreen); 
  avgMt_p->SetLineColor(kRed); 

  etaStar_h->Draw("p");
  etaStar_30_h->SetLineColor(kRed);
  etaStar_50_h->SetLineColor(kBlue);
  etaStar_70_h->SetLineColor(kGreen);
  etaStar_30_h->Draw("same p");
  etaStar_50_h->Draw("same p");
  etaStar_70_h->Draw("same p");
  c1->SaveAs("plots/etaVsMult.png");
  etaStar_h->Draw("p");
  etaStar_p_h->SetLineColor(kRed);
  etaStar_k_h->SetLineColor(kBlue);
  etaStar_pi_h->SetLineColor(kGreen);
  etaStar_p_h->Draw("same p");
  etaStar_k_h->Draw("same p");
  etaStar_pi_h->Draw("same p");
  c1->SaveAs("plots/etaVspid.png");

  multVsPt->Draw("colz");
  c1->SetLogz();
  c1->SaveAs("plots/multVsJetPt.png"); 


  c1->SetLogx(0);
  c1->SetLogy(0);
  c1->SetLogz(0);
  avgJt_k->GetYaxis()->SetRangeUser(0,2.5);
  avgJt_k->Draw("p");
  avgJt_p->Draw("p same");
  avgJt_pi->Draw("p same");
  c1->SaveAs("plots/avgJt_vsMult.png");
  
  avgMt_k->GetYaxis()->SetRangeUser(0.5,1.5);
  avgMt_k->Draw("p");
  avgMt_p->Draw("p same");
  avgMt_pi->Draw("p same");
  c1->SaveAs("plots/avgMt_vsMult.png");

  //inclusive dNdeta plot
  etaStar_h_dNdEta->Scale(1.0/etaStar_h_dNdEta->GetBinWidth(etaStar_h_dNdEta->FindBin(1)));
  etaStar_h_dNdEta->GetYaxis()->SetTitle("#frac{1}{N_{jet}} #frac{dN_{ch}}{d#eta*}");
  etaStar_h_dNdEta->GetXaxis()->SetTitle("#eta*");
  etaStar_h_dNdEta->SetMarkerColor(kBlack);
  etaStar_h_dNdEta->SetLineColor(kBlack);
  etaStar_h_dNdEta->GetYaxis()->CenterTitle();
  etaStar_h_dNdEta->GetYaxis()->SetRangeUser(0.001,20);
  etaStar_h_dNdEta->GetXaxis()->CenterTitle();
  etaStar_h_dNdEta->SetMarkerStyle(8);
  etaStar_h_dNdEta->Draw("p");
  c1->SetLogy();
  c1->SaveAs("plots/dNdEta_jet.png"); 
  c1->SetLogy(0);
  etaStar_h_dNdEta->GetYaxis()->SetRangeUser(0,7);
  c1->SaveAs("plots/dNdEta_jet_linear.png"); 
  etaStar_h_dNdEta->Write();

  c1->SetLogy();
  mult->SetMarkerColor(kBlack);
  mult->SetLineColor(kBlack);
  mult->GetYaxis()->CenterTitle();
  mult->GetYaxis()->SetRangeUser(0.000000001,70);
  mult->GetYaxis()->SetTitle("#sigma");
  mult->GetXaxis()->CenterTitle();
  mult->SetMarkerStyle(8);
  mult->Draw("p");
  c1->SaveAs("plots/multCrossSection.png");  
  
  mult->Scale(1.0/etaStarJetWeightSum);
  mult->SetMarkerColor(kBlack);
  mult->SetLineColor(kBlack);
  mult->GetYaxis()->CenterTitle();
  mult->GetYaxis()->SetRangeUser(0.000000001,0.3);
  mult->GetYaxis()->SetTitle("Normalized to Unity");
  mult->GetXaxis()->CenterTitle();
  mult->SetMarkerStyle(8);
  mult->Draw("p");

  c1->SaveAs("plots/multNormalizedToUnity.png");  

  output->Write();

  output->Close();
}

//Code enters execution here
int main(int argc, const char* argv[])
{
  if(argc != 2)
  {
    std::cout << "Usage: Z_mumu_Channel <fileList>" << std::endl;
    return 1;
  }  


  //read input parameters
  std::string fList = argv[1];
  std::string buffer;
  std::vector<std::string> listOfFiles;
  std::ifstream inFile(fList.data());


  //read the file list and spit it into a vector of strings based on how the parallelization is to be done
  //each vector is a separate subset of the fileList based on the job number
  if(!inFile.is_open())
  {
    std::cout << "Error opening jet file. Exiting." <<std::endl;
    return 1;
  }
  else
  {
    int line = 0;
    while(true)
    {
      inFile >> buffer;
      if(inFile.eof()) break;
      listOfFiles.push_back(buffer);
      line++;
    }
  }

  analyze(listOfFiles);

  return 0; 
}
