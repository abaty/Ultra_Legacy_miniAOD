#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMarker.h"


void GenPlotter(){
  TH1::SetDefaultSumw2();
  gStyle->SetLegendBorderSize(0);

  //multiplicity
  TFile * f = TFile::Open("../rootFiles/PythiaOutput.root","read");
  TH1D * multiplicity500 = (TH1D*) f->Get("multiplicity500");
  TH1D * multiplicity800 = (TH1D*) f->Get("multiplicity800");

  TCanvas * c1 = new TCanvas("c1","c1",1000,800);
  c1->SetLogy();
  c1->SetLeftMargin(0.2); 
  c1->SetTickx(1);
  c1->SetTicky(1);

  multiplicity500->SetMarkerColor(kBlack);
  multiplicity500->SetLineColor(kBlack);
  multiplicity500->SetMarkerStyle(8);
  multiplicity500->GetYaxis()->CenterTitle();
  multiplicity500->GetXaxis()->CenterTitle();
  multiplicity500->GetXaxis()->SetTitle("N_{ch}");
  multiplicity500->SetStats(0);

  multiplicity500->GetYaxis()->SetRangeUser(0.00000001,0.1);
  multiplicity500->Draw("p");

  multiplicity800->SetMarkerColor(kRed+1);
  multiplicity800->SetMarkerStyle(21);
  multiplicity800->SetLineColor(kRed+1);
  multiplicity800->Draw("same p");

  TLegend * l = new TLegend(0.25,0.2,0.6,0.5);
  l->SetFillStyle(0);
  l->AddEntry((TObject*)0,"PYTHIA 8 pp 13 TeV","");
  l->AddEntry((TObject*)0,"Anti k_{T} R=0.8","");
  l->AddEntry(multiplicity500,"jet p_{T} > 500 GeV","p");
  l->AddEntry(multiplicity800,"jet p_{T} > 800 GeV","p");
  l->Draw("same");  


  c1->SaveAs("../plots/multiplicity500.png");
  c1->SaveAs("../plots/multiplicity500.pdf");
  c1->SaveAs("../plots/multiplicity500.C");
 
  //dNdEta
  TH1D * dNdEta500_0_999 = (TH1D*) f->Get("dNdEta500_0_999");
  TH1D * dNdEta800_0_999 = (TH1D*) f->Get("dNdEta800_0_999");
  TH1D * dNdEta500_0_35 = (TH1D*) f->Get("dNdEta500_0_35");
  TH1D * dNdEta800_0_35 = (TH1D*) f->Get("dNdEta800_0_35");
  TH1D * dNdEta500_100_999 = (TH1D*) f->Get("dNdEta500_100_999");
  TH1D * dNdEta800_100_999 = (TH1D*) f->Get("dNdEta800_100_999");
 
  dNdEta500_0_999->SetMarkerColor(kBlack);
  dNdEta500_0_999->SetLineColor(kBlack);
  dNdEta500_0_999->SetMarkerStyle(8);
  dNdEta500_0_999->GetYaxis()->CenterTitle();
  dNdEta500_0_999->GetXaxis()->CenterTitle();
  dNdEta500_0_999->GetXaxis()->SetTitle("#eta*");
  dNdEta500_0_999->GetYaxis()->SetTitle("#frac{1}{N_{jet}} #frac{N_{ch}}{#eta*}");
  dNdEta500_0_999->SetStats(0);

  dNdEta500_0_999->GetYaxis()->SetRangeUser(0.01,100);
  dNdEta500_0_999->Draw("p");
  
  dNdEta800_0_999->SetMarkerColor(kBlack);
  dNdEta800_0_999->SetLineColor(kBlack);
  dNdEta800_0_999->SetMarkerStyle(24);
  dNdEta800_0_999->Draw("p same");
  
  dNdEta500_0_35->SetMarkerColor(kRed+1);
  dNdEta500_0_35->SetLineColor(kRed+1);
  dNdEta500_0_35->SetMarkerStyle(21);
  dNdEta500_0_35->Draw("p same");
  
  dNdEta800_0_35->SetMarkerColor(kRed+1);
  dNdEta800_0_35->SetLineColor(kRed+1);
  dNdEta800_0_35->SetMarkerStyle(25);
  dNdEta800_0_35->Draw("p same");
  
  dNdEta500_100_999->SetMarkerColor(kBlue);
  dNdEta500_100_999->SetLineColor(kBlue);
  dNdEta500_100_999->SetMarkerStyle(34);
  dNdEta500_100_999->Draw("p same");
  
  dNdEta800_100_999->SetMarkerColor(kBlue);
  dNdEta800_100_999->SetLineColor(kBlue);
  dNdEta800_100_999->SetMarkerStyle(28);
  dNdEta800_100_999->Draw("p same");
  
  TLegend * l2 = new TLegend(0.25,0.2,0.6,0.5);
  l2->SetFillStyle(0);
  l2->AddEntry((TObject*)0,"PYTHIA 8 pp 13 TeV","");
  l2->AddEntry((TObject*)0,"Anti k_{T} R=0.8","");
  l2->AddEntry(dNdEta500_0_999,"p_{T} > 500 GeV","p");
  l2->AddEntry(dNdEta800_0_999,"p_{T} > 800 GeV","p");
  l2->AddEntry(dNdEta500_0_35,"p_{T} > 500 GeV, N_{ch}<35","p");
  l2->AddEntry(dNdEta800_0_35,"p_{T} > 800 GeV, N_{ch}<35","p");
  l2->AddEntry(dNdEta500_100_999,"p_{T} > 500 GeV, N_{ch}>100","p");
  l2->AddEntry(dNdEta800_100_999,"p_{T} > 800 GeV, N_{ch}>100","p");
  l2->Draw("same");  
  
  c1->SaveAs("../plots/dNdEta.png");
  c1->SaveAs("../plots/dNdEta.pdf");
  c1->SaveAs("../plots/dNdEta.C");


  //mean m_{T}
  c1->SetLogy(0);
  TH1D * mt500_0 = (TH1D*) f->Get("avgMtVsMult500_0");
  TH1D * mt500_1 = (TH1D*) f->Get("avgMtVsMult500_1");
  TH1D * mt500_2 = (TH1D*) f->Get("avgMtVsMult500_2");
  TH1D * mt500_3 = (TH1D*) f->Get("avgMtVsMult500_3");


  mt500_0->SetMarkerColor(kBlack);
  mt500_0->SetLineColor(kBlack);
  mt500_0->SetMarkerStyle(8);
  mt500_0->GetYaxis()->CenterTitle();
  mt500_0->GetXaxis()->CenterTitle();
  mt500_0->GetXaxis()->SetTitle("N_{ch}");
  mt500_0->GetYaxis()->SetTitle("<m_{T}>");
  mt500_0->SetStats(0);

  mt500_0->GetYaxis()->SetRangeUser(0,4);
  mt500_0->Draw("p");


  mt500_1->SetMarkerColor(kRed+1);
  mt500_1->SetLineColor(kRed+1);
  mt500_1->SetMarkerStyle(24);
  mt500_1->Draw("p same");

  mt500_2->SetMarkerColor(kBlue);
  mt500_2->SetLineColor(kBlue);
  mt500_2->SetMarkerStyle(21);
  mt500_2->Draw("p same");
  
  mt500_3->SetMarkerColor(kGreen+2);
  mt500_3->SetLineColor(kGreen+2);
  mt500_3->SetMarkerStyle(25);
  mt500_3->Draw("p same");

  TLegend * l3 = new TLegend(0.22,0.58,0.57,0.88);
  l3->SetFillStyle(0);
  l3->AddEntry((TObject*)0,"PYTHIA 8 pp 13 TeV","");
  l3->AddEntry((TObject*)0,"Anti k_{T} R=0.8","");
  l3->AddEntry((TObject*)0,"jet p_{T} > 500 GeV","");
  l3->AddEntry(mt500_0,"Charged Hadrons","p");
  l3->AddEntry(mt500_1,"p","p");
  l3->AddEntry(mt500_2,"K^{#pm}","p");
  l3->AddEntry(mt500_3,"#pi^{#pm}","p");
  l3->Draw("same");

  c1->SaveAs("../plots/avgMtVsMult.png");
  c1->SaveAs("../plots/avgMtVsMult.pdf");
  c1->SaveAs("../plots/avgMtVsMult.C");

  
  //total yield ratios
  c1->SetLogy();
  TH1D * yield[7];
  yield[0] = (TH1D*) f->Get("avgYieldVsMult500_0");
  yield[1] = (TH1D*) f->Get("avgYieldVsMult500_1");
  yield[2] = (TH1D*) f->Get("avgYieldVsMult500_2");
  yield[3] = (TH1D*) f->Get("avgYieldVsMult500_3");
  yield[4] = (TH1D*) f->Get("avgYieldVsMult500_4");
  yield[5] = (TH1D*) f->Get("avgYieldVsMult500_5");
  yield[6] = (TH1D*) f->Get("avgYieldVsMult500_6");
  TH1D * dummy = (TH1D*) f->Get("avgYieldVsMult500Dummy");  

  TH1D * yieldRatio[7];
  for(int i = 0; i<7; i++){
    yieldRatio[i] = (TH1D*)dummy->Clone(Form("yieldRatio_%d",i));
    for(int j = 0; j<yield[0]->GetNbinsX()+2; j++){
      yieldRatio[i]->SetBinContent(j, yield[i]->GetBinContent(j)/yield[0]->GetBinContent(j));
      float relError1 = yield[0]->GetBinError(j)/yield[0]->GetBinContent(j);
      float relError2 = yield[i]->GetBinError(j)/yield[i]->GetBinContent(j);
      float relErrorNet = TMath::Sqrt( relError1*relError1 + relError2*relError2 );
      yieldRatio[i]->SetBinError(j, relErrorNet * yieldRatio[i]->GetBinContent(j) );
    }
    yieldRatio[i]->Print("All");
  }

  yieldRatio[1]->SetMarkerColor(kBlack);
  yieldRatio[1]->SetLineColor(kBlack);
  yieldRatio[1]->SetMarkerStyle(8);
  yieldRatio[1]->GetYaxis()->CenterTitle();
  yieldRatio[1]->GetXaxis()->CenterTitle();
  yieldRatio[1]->GetXaxis()->SetTitle("N_{ch}");
  yieldRatio[1]->GetYaxis()->SetTitle("#frac{<Yield>}{<Yield_{#pi^{#pm}}>}");
  yieldRatio[1]->SetStats(0);

  yieldRatio[1]->GetYaxis()->SetRangeUser(0.00001,200);
  yieldRatio[1]->Draw("p");
  
  yieldRatio[2]->SetMarkerColor(kRed+1);
  yieldRatio[2]->SetLineColor(kRed+1);
  yieldRatio[2]->SetMarkerStyle(24);
  yieldRatio[2]->Draw("p same");

  yieldRatio[3]->SetMarkerColor(kBlue);
  yieldRatio[3]->SetLineColor(kBlue);
  yieldRatio[3]->SetMarkerStyle(21);
  yieldRatio[3]->Draw("p same");
  
  yieldRatio[4]->SetMarkerColor(kGreen+2);
  yieldRatio[4]->SetLineColor(kGreen+2);
  yieldRatio[4]->SetMarkerStyle(25);
  yieldRatio[4]->Draw("p same");

  yieldRatio[5]->SetMarkerColor(kViolet);
  yieldRatio[5]->SetLineColor(kViolet);
  yieldRatio[5]->SetMarkerStyle(34);
  yieldRatio[5]->Draw("p same");
  
  yieldRatio[6]->SetMarkerColor(kCyan+1);
  yieldRatio[6]->SetLineColor(kCyan+1);
  yieldRatio[6]->SetMarkerStyle(28);
  yieldRatio[6]->Draw("p same");

  TLegend * l4 = new TLegend(0.22,0.63,0.57,0.88);
  l4->SetFillStyle(0);
  l4->AddEntry((TObject*)0,"PYTHIA 8 pp 13 TeV","");
  l4->AddEntry((TObject*)0,"Anti k_{T} R=0.8","");
  l4->AddEntry((TObject*)0,"jet p_{T} > 500 GeV","");
  l4->AddEntry(yieldRatio[1],"p","p");
  l4->AddEntry(yieldRatio[2],"K^{#pm}","p");
  l4->AddEntry(yieldRatio[3],"K^{0}_{S}","p");
  l4->AddEntry(yieldRatio[4],"#Lambda + #bar{#Lambda}","p");
  l4->AddEntry(yieldRatio[5],"#Xi^{-} + #bar{#Xi^{-}}","p");
  l4->AddEntry(yieldRatio[6],"#Omega^{-} + #bar{#Omega^{-}}","p");
  l4->Draw("same");

  c1->SaveAs("../plots/yieldRatiosVsMult.png");
  c1->SaveAs("../plots/yieldRatiosVsMult.pdf");
  c1->SaveAs("../plots/yieldRatiosVsMult.C");
 
  c1->SetLogy(0);
  c1->SetLogz();

  //photons
  c1->SetLogy();
  TH1D * gamma_0_35[5];
  TH1D * gamma_100_999[5];
  for(int i = 0; i<5; i++){
   gamma_0_35[i] = (TH1D*) f->Get(Form("gammaJt500_0_35_%d",i));
   gamma_100_999[i] = (TH1D*) f->Get(Form("gammaJt500_100_999_%d",i));
  }
  gamma_0_35[0]->Print("All");
  gamma_0_35[0]->SetMarkerColor(kBlack);
  gamma_0_35[0]->SetLineColor(kBlack);
  gamma_0_35[0]->SetMarkerStyle(8);
  gamma_0_35[0]->GetYaxis()->CenterTitle();
  gamma_0_35[0]->GetXaxis()->CenterTitle();
  gamma_0_35[0]->GetXaxis()->SetTitle("j_{T}");
  gamma_0_35[0]->GetYaxis()->SetTitle("#frac{1}{N_{jet}} #frac{dN_{#gamma}}{dj_{T}}");
  gamma_0_35[0]->GetYaxis()->SetTitleOffset(2.0);
  gamma_0_35[0]->SetStats(0);

  gamma_0_35[0]->GetYaxis()->SetRangeUser(0.00001,2000);
  gamma_0_35[0]->Draw("p");
  
  gamma_0_35[1]->SetMarkerColor(kRed+1);
  gamma_0_35[1]->SetLineColor(kRed+1);
  gamma_0_35[1]->SetMarkerStyle(24);
  gamma_0_35[1]->Draw("p same");

  gamma_0_35[2]->SetMarkerColor(kBlue);
  gamma_0_35[2]->SetLineColor(kBlue);
  gamma_0_35[2]->SetMarkerStyle(21);
  gamma_0_35[2]->Draw("p same");
  
  gamma_0_35[3]->SetMarkerColor(kGreen+2);
  gamma_0_35[3]->SetLineColor(kGreen+2);
  gamma_0_35[3]->SetMarkerStyle(25);
  gamma_0_35[3]->Draw("p same");

  gamma_0_35[4]->SetMarkerColor(kViolet);
  gamma_0_35[4]->SetLineColor(kViolet);
  gamma_0_35[4]->SetMarkerStyle(34);
  gamma_0_35[4]->Draw("p same");

  TLegend * l5 = new TLegend(0.52,0.63,0.87,0.88);
  l5->SetFillStyle(0);
  l5->AddEntry((TObject*)0,"PYTHIA 8 pp 13 TeV","");
  l5->AddEntry((TObject*)0,"Anti k_{T} R=0.8","");
  l5->AddEntry((TObject*)0,"jet p_{T} > 500 GeV, N_{ch}<35","");
  l5->AddEntry(gamma_0_35[0],"Inclusive #gamma","p");
  l5->AddEntry(gamma_0_35[1],"Beam","p");
  l5->AddEntry(gamma_0_35[2],"Parton Shower","p");
  l5->AddEntry(gamma_0_35[3],"#pi^{0} decay","p");
  l5->AddEntry(gamma_0_35[4],"Other hadronic decay","p");
  l5->Draw("same");

  c1->SaveAs("../plots/gamma500_0_35.png");
  c1->SaveAs("../plots/gamma500_0_35.pdf");
  c1->SaveAs("../plots/gamma500_0_35.C");
  
  gamma_100_999[0]->Print("All");
  gamma_100_999[0]->SetMarkerColor(kBlack);
  gamma_100_999[0]->SetLineColor(kBlack);
  gamma_100_999[0]->SetMarkerStyle(8);
  gamma_100_999[0]->GetYaxis()->CenterTitle();
  gamma_100_999[0]->GetXaxis()->CenterTitle();
  gamma_100_999[0]->GetXaxis()->SetTitle("j_{T}");
  gamma_100_999[0]->GetYaxis()->SetTitle("#frac{1}{N_{jet}} #frac{dN_{#gamma}}{dj_{T}}");
  gamma_100_999[0]->GetYaxis()->SetTitleOffset(2.0);
  gamma_100_999[0]->SetStats(0);

  gamma_100_999[0]->GetYaxis()->SetRangeUser(0.00001,2000);
  gamma_100_999[0]->Draw("p");
  
  gamma_100_999[1]->SetMarkerColor(kRed+1);
  gamma_100_999[1]->SetLineColor(kRed+1);
  gamma_100_999[1]->SetMarkerStyle(24);
  gamma_100_999[1]->Draw("p same");

  gamma_100_999[2]->SetMarkerColor(kBlue);
  gamma_100_999[2]->SetLineColor(kBlue);
  gamma_100_999[2]->SetMarkerStyle(21);
  gamma_100_999[2]->Draw("p same");
  
  gamma_100_999[3]->SetMarkerColor(kGreen+2);
  gamma_100_999[3]->SetLineColor(kGreen+2);
  gamma_100_999[3]->SetMarkerStyle(25);
  gamma_100_999[3]->Draw("p same");

  gamma_100_999[4]->SetMarkerColor(kViolet);
  gamma_100_999[4]->SetLineColor(kViolet);
  gamma_100_999[4]->SetMarkerStyle(34);
  gamma_100_999[4]->Draw("p same");

  TLegend * l6 = new TLegend(0.52,0.63,0.87,0.88);
  l6->SetFillStyle(0);
  l6->AddEntry((TObject*)0,"PYTHIA 8 pp 13 TeV","");
  l6->AddEntry((TObject*)0,"Anti k_{T} R=0.8","");
  l6->AddEntry((TObject*)0,"jet p_{T} > 500 GeV, N_{ch}#geq100, #gamma p_{T} > 2 GeV","");
  l6->AddEntry(gamma_100_999[0],"Inclusive #gamma","p");
  l6->AddEntry(gamma_100_999[1],"Beam","p");
  l6->AddEntry(gamma_100_999[2],"Parton Shower","p");
  l6->AddEntry(gamma_100_999[3],"#pi^{0} decay","p");
  l6->AddEntry(gamma_100_999[4],"Other hadronic decay","p");
  l6->Draw("same");

  c1->SaveAs("../plots/gamma500_100_999.png");
  c1->SaveAs("../plots/gamma500_100_999.pdf");
  c1->SaveAs("../plots/gamma500_100_999.C");

  
  //jtVsEta
 
  TH2D * jtVEta_0_35 = (TH2D*) f->Get("jtVsEta500_0_35");
  TH2D * jtVEta_100_999 = (TH2D*) f->Get("jtVsEta500_100_999");
  
  jtVEta_0_35->SetStats(0); 
  jtVEta_0_35->GetXaxis()->CenterTitle(); 
  jtVEta_0_35->GetYaxis()->CenterTitle(); 
  jtVEta_0_35->GetXaxis()->SetTitle("#eta*"); 
  jtVEta_0_35->GetYaxis()->SetTitle("j_{T}"); 
  jtVEta_0_35->Draw("colz");
  c1->SaveAs("../plots/jtVsEta_0_35.png");
  c1->SaveAs("../plots/jtVsEta_0_35.pdf");
  c1->SaveAs("../plots/jtVsEta_0_35.C");


  jtVEta_100_999->SetStats(0); 
  jtVEta_100_999->GetXaxis()->CenterTitle(); 
  jtVEta_100_999->GetYaxis()->CenterTitle(); 
  jtVEta_100_999->GetXaxis()->SetTitle("#eta*"); 
  jtVEta_100_999->GetYaxis()->SetTitle("j_{T}"); 
  jtVEta_100_999->Draw("colz");
  c1->SaveAs("../plots/jtVsEta_100_999.png");
  c1->SaveAs("../plots/jtVsEta_100_999.pdf");
  c1->SaveAs("../plots/jtVsEta_100_999.C");

}
