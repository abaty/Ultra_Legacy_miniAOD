#include <string>
#include "include/Timer.h"
#include <TStyle.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
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

void analyze( std::vector< std::string> files, std::vector< std::string> files2, bool redoWeightHist, bool isHerwig){
  TH1::SetDefaultSumw2();


  //pthat combination stufg
  const int nqHats = 8;
  const float qHatBoundaries[nqHats+1] = {470, 600, 800, 1000, 1400, 1800, 2400, 3200, 13000};
      //units in pb
  float xs[nqHats] = {552.1, 156.5, 26.28, 7.47, 0.6484, 0.08743, 0.005236, 0.0001357};
  const int nMultBins = 5;
  int multBinLow[nMultBins] = {0,0,35,70,100};
  int multBinHigh[nMultBins] = {999,35,70,100,999};
  
  const int nMultBinsFine = 23;
  float multBinsFine[nMultBinsFine+1] = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,140};
  
  const int nMultBinsFine2 = 12;
  float multBinsFine2[nMultBinsFine2+1] = {0,10,20,30,40,50,60,70,80,90,100,120,140};

  TH1D * qHatHist; 
  TH1D * qHatHist2; 

  TFile * weights; 
  if(!isHerwig){ 
    if(redoWeightHist){
      weights = TFile::Open("rootFiles/weights.root","recreate");  
      qHatHist= new TH1D("qhatHist",";;qHat",nqHats,qHatBoundaries);
      qHatHist2= new TH1D("qhatHist2",";;qHat2",nqHats,qHatBoundaries);
    } else{
      weights = TFile::Open("rootFiles/weights.root","read");  
      qHatHist= (TH1D*) weights->Get("qhatHist");
      qHatHist2= (TH1D*) weights->Get("qhatHist2");
      qHatHist->Print("All");
      qHatHist2->Print("All");
    }
  }
  
  //branch definitions    
  float genWeight = 0;
  float genQScale = 0;
  std::vector< float > * genJetPt = 0;
  std::vector< float > * genJetEta = 0;
  std::vector< float > * genJetPhi = 0;
  std::vector< int > * genJetChargedMultiplicity = 0;
  std::vector< std::vector< int > > * genDau_chg = 0;
  std::vector< std::vector< int > > * genDau_pid = 0;
  std::vector< std::vector< int > > * genDau_mom = 0;
  std::vector< std::vector< float > > * genDau_pt = 0;
  std::vector< std::vector< float > > * genDau_eta = 0;
  std::vector< std::vector< float > > * genDau_phi = 0;

  if(redoWeightHist && !isHerwig){
    std::cout << "calculating weighting factors for qhat bins" << std::endl;
    for(unsigned int f = 0; f<files.size(); f++){
      std::cout << f << "/" << files.size() << std::endl;

      TFile * inputFile = TFile::Open(files.at(f).c_str(),"read");
      TTree * t = (TTree*)inputFile->Get("analyzer/trackTree");

      t->SetBranchAddress("genQScale",&genQScale);
      for( int i = 0; i<t->GetEntries(); i++){
        genQScale = 0;
        t->GetEntry(i);

        qHatHist->Fill(genQScale);    
      }
    
      inputFile->Close();
    } 
    qHatHist->Print("All");
    
    //repeat for filtered samples (copy pasted :( )
    std::cout << "calculating weighting factors for qhat bins" << std::endl;
    for(unsigned int f = 0; f<files2.size(); f++){
      std::cout << f << "/" << files2.size() << std::endl;

      TFile * inputFile = TFile::Open(files2.at(f).c_str(),"read");
      TTree * t = (TTree*)inputFile->Get("analyzer/trackTree");

      t->SetBranchAddress("genQScale",&genQScale);
      for( int i = 0; i<t->GetEntries(); i++){
        genQScale = 0;
        t->GetEntry(i);

        qHatHist2->Fill(genQScale);    
      }
    
      inputFile->Close();
    } 
    qHatHist2->Print("All");
    weights->Write();
  }
 
  TFile * output;
  if(!isHerwig)  output = TFile::Open("rootFiles/PythiaOutput.root","recreate");
  else           output = TFile::Open("rootFiles/HerwigOutput.root","recreate");
  //TH1D * pthat = new TH1D("pthat",";#hat{q};#sigma (pb)",303,470,3500);
  TH1D * multiplicity = new TH1D("multiplicity",";n_{ch};Normalized to Unity",70,0,140);
  TH1D * multiplicity500 = new TH1D("multiplicity500",";n_{ch};Normalized to Unity",70,0,140);
  TH1D * multiplicity800 = new TH1D("multiplicity800",";n_{ch};Normalized to Unity",70,0,140);

  TProfile * avgJtVsMult500[4];
  TProfile * avgMtVsMult500[4];
  TProfile * avgYieldVsMult500[8];
  for(int i = 0; i<8; i++){
    if(i<4){
      avgJtVsMult500[i] = new  TProfile(Form("avgJtVsMult500_%d",i),";n_{ch};<j_{t}>",nMultBinsFine,multBinsFine);
      avgMtVsMult500[i] = new  TProfile(Form("avgMtVsMult500_%d",i),";n_{ch};<m_{t}>",nMultBinsFine,multBinsFine);
    }
    avgYieldVsMult500[i] = new  TProfile(Form("avgYieldVsMult500_%d",i),";n_{ch};<Yield per jet>",nMultBinsFine,multBinsFine);
  }
  TH1D * avgYieldVsMult500Dummy = new TH1D("avgYieldVsMult500Dummy",";n_{ch};<Yield per jet>",nMultBinsFine,multBinsFine);

  //0 = inclusive, 1 = proton-initiated, 2 = shower photon, 3 = pi0 decays , 4 = other hadron decays
  TH1D * gammaJt500[nMultBins][5];
  float gammaJt500_w[nMultBins][5] = {0};
  
  TH1D * dNdEta500[nMultBins];
  TH1D * dNdEta800[nMultBins];
  TH2D * jtVsEta500[nMultBins];
  float dNdEta500_w[nMultBins] = {0};
  float dNdEta800_w[nMultBins] = {0};
  for(int i = 0; i<nMultBins; i++){
    dNdEta500[i] = new TH1D(Form("dNdEta500_%d_%d",multBinLow[i],multBinHigh[i]),";#eta*;#frac{dN}{d#eta*}",80,0,8);
    dNdEta800[i] = new TH1D(Form("dNdEta800_%d_%d",multBinLow[i],multBinHigh[i]),";#eta*;#frac{dN}{d#eta*}",80,0,8);
    jtVsEta500[i] = new TH2D(Form("jtVsEta500_%d_%d",multBinLow[i],multBinHigh[i]),";#eta*;j_{t}",40,0,8,50,0,10);
    for(int j = 0; j<5; j++){
      gammaJt500[i][j] = new TH1D(Form("gammaJt500_%d_%d_%d",multBinLow[i],multBinHigh[i],j),";j_{t};#frac{1}{N_{jet}}#frac{dN}{dj_{T}}",20,0,10);
    }
  }


  //file loop
  for(int x = 0; x<2; x++){
    std::vector< std::string > fileIter;
    if(x==0) fileIter = files;
    if(x==1) fileIter = files2;
    for(unsigned int f = 0; f<fileIter.size(); f++){
      //for testing
      //if(f%10!=0) continue;

      //if(f%10==0) std::cout << f << "/" << files.size() << std::endl;
      std::cout << f << "/" << fileIter.size() << std::endl;

      TFile * inputFile = TFile::Open(fileIter.at(f).c_str(),"read");
      TTree * t = (TTree*)inputFile->Get("analyzer/trackTree");
   
      //branch definitions
      if(isHerwig) t->SetBranchAddress("genWeight",&genWeight);
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
      t->SetBranchAddress("genDau_mom",&genDau_mom);   
 
      //event loop
      for(int i = 0; i<t->GetEntries(); i++){
        //zero out vectors and get entry
        genQScale = 0;
        genWeight = 0;
        *genJetPt = std::vector< float >();
        *genJetEta = std::vector< float >();
        *genJetPhi = std::vector< float >();
        *genJetChargedMultiplicity = std::vector< int >();
        *genDau_chg = std::vector< std::vector< int > >();
        *genDau_pid = std::vector< std::vector< int > >();
        *genDau_mom = std::vector< std::vector< int > >();
        *genDau_pt =  std::vector< std::vector< float > >();
        *genDau_eta = std::vector< std::vector< float > >();
        *genDau_phi = std::vector< std::vector< float > >();
        t->GetEntry(i);

        //hard cut at pthat of 1000 to prevent low-stats anomalies from higher pthats
        if(!isHerwig && genQScale>999.9) continue;

        //weight by xsection/total number of gen events in the pthat bin
        float w = 1;
        if(!isHerwig){
          if(x==0){
            w = xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale));
          }
          if(x==1) {
            w = xs[qHatHist2->FindBin(genQScale) - 1 ] / qHatHist2->GetBinContent(qHatHist2->FindBin(genQScale)); 
          }
        } else{
          w = genWeight; 
        }

        //jet spectra
        if(genJetPt->size()==0) continue;
        if(genJetChargedMultiplicity->size()==0) continue;
        
        for(int j = 0; j<genJetPt->size(); j++){
          if( TMath::Abs( genJetEta->at(j) ) > 2) continue;
          if( (genJetChargedMultiplicity->at(j) <  60 && x==0 ) || 
              (genJetChargedMultiplicity->at(j) >= 60 && x==1) ){

            //multiplicity distributions
            multiplicity->Fill( genJetChargedMultiplicity->at(j), w);
            if(genJetPt->at(j) > 500){
              multiplicity500->Fill(genJetChargedMultiplicity->at(j),w);
             
              for(int l = 0; l<nMultBins; l++){
                if(genJetChargedMultiplicity->at(j) >= multBinLow[l] && genJetChargedMultiplicity->at(j) < multBinHigh[l]){
                  dNdEta500_w[l] += w;
                  for(int m = 0; m<5; m++) gammaJt500_w[l][m] += w;
                }
              }
            }
            if(genJetPt->at(j) > 800){
              multiplicity800->Fill(genJetChargedMultiplicity->at(j),w);
              
              for(int l = 0; l<nMultBins; l++){
                if(genJetChargedMultiplicity->at(j) >= multBinLow[l] && genJetChargedMultiplicity->at(j) < multBinHigh[l]) dNdEta800_w[l] += w;
              }
            } 
  
            //daughter loop
            float yieldCounter[8] = {0};
            for(int k = 0; k<(genDau_chg->at(j)).size(); k++){

              float ptStar = ptWRTJet(genJetPt->at(j), genJetEta->at(j), genJetPhi->at(j), (genDau_pt->at(j)).at(k), (genDau_eta->at(j)).at(k), (genDau_phi->at(j)).at(k));
              
              //photon stuff first
              if(( genDau_pid->at(j)).at(k) == 22 && (genDau_pt->at(j)).at(k) > 2.0 ){
                if(genJetPt->at(j) > 500){
                  
                  int momPidIndx = -9999;
                  if(      (genDau_mom->at(j)).at(k) == 2212 ) momPidIndx = 1;
                  else if( TMath::Abs((genDau_mom->at(j)).at(k)) < 10 || (genDau_mom->at(j)).at(k) == 21 ) momPidIndx = 2;
                  else if( (genDau_mom->at(j)).at(k) == 111 ) momPidIndx = 3;
                  else momPidIndx = 4;
                  
                  for(int l = 0; l<nMultBins; l++){
                    if(genJetChargedMultiplicity->at(j) >= multBinLow[l] && genJetChargedMultiplicity->at(j) < multBinHigh[l]){
                      gammaJt500[l][0]->Fill(ptStar, w);
                      gammaJt500[l][momPidIndx]->Fill(ptStar, w);
                    }
                  }
                }  
              }
  
              //neutral yield counters
              if( TMath::Abs((genDau_pid->at(j)).at(k)) == 310 ) yieldCounter[3]++;//k0s
              if( TMath::Abs((genDau_pid->at(j)).at(k)) == 3122 ) yieldCounter[4]++;//lambda
              //if( TMath::Abs((genDau_pid->at(j)).at(k)) == 111 ) yieldCounter[7]++;//pi0
              
              //now charged stuff
              if( (genDau_chg->at(j)).at(k) == 0 ) continue;

              float etaStar = etaWRTJet(genJetPt->at(j), genJetEta->at(j), genJetPhi->at(j), (genDau_pt->at(j)).at(k), (genDau_eta->at(j)).at(k), (genDau_phi->at(j)).at(k));
              float phiStar = phiWRTJet(genJetPt->at(j), genJetEta->at(j), genJetPhi->at(j), (genDau_pt->at(j)).at(k), (genDau_eta->at(j)).at(k), (genDau_phi->at(j)).at(k));
              float mtStar = 0;
              int pidIndx = -1;
              if( TMath::Abs((genDau_pid->at(j)).at(k)) == 2212){
                mtStar = TMath::Sqrt(ptStar*ptStar+0.938*0.938) - 0.938;
                pidIndx = 1;
              }
              if( TMath::Abs((genDau_pid->at(j)).at(k)) == 321 ){
                mtStar = TMath::Sqrt(ptStar*ptStar+0.494*0.494) - 0.494;
                pidIndx = 2;
              }
              if( TMath::Abs((genDau_pid->at(j)).at(k)) == 211 ){
                mtStar = TMath::Sqrt(ptStar*ptStar+0.1395*0.1395) - 0.1395;  
                pidIndx = 3;
              }
           
              //eta distributions
              if(genJetPt->at(j) > 500){
                for(int l = 0; l<nMultBins; l++){
                  if(genJetChargedMultiplicity->at(j) >= multBinLow[l] && genJetChargedMultiplicity->at(j) < multBinHigh[l]){
                    dNdEta500[l]->Fill(  etaStar  ,w);
                    jtVsEta500[l]->Fill(etaStar,ptStar,w);
                  }
                }
                //jt and mt profiles
                avgJtVsMult500[0]->Fill( genJetChargedMultiplicity->at(j)  , ptStar, w);
                avgMtVsMult500[0]->Fill( genJetChargedMultiplicity->at(j)  , mtStar, w);
                if(pidIndx>-1){
                  avgJtVsMult500[pidIndx]->Fill( genJetChargedMultiplicity->at(j)  , ptStar, w);
                  avgMtVsMult500[pidIndx]->Fill( genJetChargedMultiplicity->at(j)  , mtStar, w);
                }

                //charged yield counters
                if( TMath::Abs((genDau_pid->at(j)).at(k)) == 211 ) yieldCounter[0]++;//ch. pion
                if( TMath::Abs((genDau_pid->at(j)).at(k)) == 2212 ) yieldCounter[1]++;//proton
                if( TMath::Abs((genDau_pid->at(j)).at(k)) == 321 ) yieldCounter[2]++;//ch. kaon
                if( TMath::Abs((genDau_pid->at(j)).at(k)) == 3312 ) yieldCounter[5]++;//cascade
                if( TMath::Abs((genDau_pid->at(j)).at(k)) == 3334 ) yieldCounter[6]++;//omega
 
              }
              if(genJetPt->at(j) > 800){
                for(int l = 0; l<nMultBins; l++){
                  if(genJetChargedMultiplicity->at(j) >= multBinLow[l] && genJetChargedMultiplicity->at(j) < multBinHigh[l]) dNdEta800[l]->Fill(  etaStar  ,w);
                }
              }
            }//close daughter loop

            //fill yield counters
            for(int k = 0; k<8; k++){
              avgYieldVsMult500[k]->Fill(genJetChargedMultiplicity->at(j) ,yieldCounter[k],w);           
            }                      
          }//close if statement
        }//close jet  loop
      }//close event loop
      inputFile->Close();
    }//close file loop
  }//close x loop
 
  float scale = multiplicity->Integral();
  multiplicity->Scale(1.0/scale);
  scale = multiplicity500->Integral();
  multiplicity500->Scale(1.0/scale);
  scale = multiplicity800->Integral();
  multiplicity800->Scale(1.0/scale);


  for(int l = 0; l<nMultBins; l++){
    float binWidth = dNdEta500[l]->GetXaxis()->GetBinWidth(1);
    dNdEta500[l]->Scale(1.0/(binWidth*dNdEta500_w[l]) ); 
    
    binWidth = dNdEta800[l]->GetXaxis()->GetBinWidth(1);
    dNdEta800[l]->Scale(1.0/(binWidth*dNdEta800_w[l]) ); 

    for(int i = 0; i<5; i++){
      binWidth = gammaJt500[l][i]->GetXaxis()->GetBinWidth(1);
      gammaJt500[l][i]->Scale(1.0/(binWidth*gammaJt500_w[l][i]) ); 
    }
  }

  output->Write();
  output->Close();
  

 
  /* 
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
  */

  
  if(!isHerwig)  weights->Close();
}

//Code enters execution here
int main(int argc, const char* argv[])
{
  if(argc != 5)
  {
    std::cout << "Usage: Z_mumu_Channel <fileList1> <fileListFiltered> <redo Weight hist> <isHerwig>" << std::endl;
    return 1;
  }  


  //read input parameters
  std::string fList = argv[1];
  std::string fList2 = argv[2];
  std::string buffer;
  std::vector<std::string> listOfFiles;
  std::ifstream inFile(fList.data());
  
  std::string buffer2;
  std::vector<std::string> listOfFiles2;
  std::ifstream inFile2(fList2.data());

  bool redoWeightHist = (bool)( std::stoi(argv[3]) );
  bool isHerwig = (bool)( std::stoi(argv[4]) );

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

  if(!inFile2.is_open())
  {
    std::cout << "Error opening jet file. Exiting." <<std::endl;
    return 1;
  }
  else
  {
    int line = 0;
    while(true)
    {
      inFile2 >> buffer2;
      if(inFile2.eof()) break;
      listOfFiles2.push_back(buffer2);
      line++;
    }
  }

  analyze(listOfFiles, listOfFiles2, redoWeightHist, isHerwig);

  return 0; 
}
