#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <TChain.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TDirectory.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>

/////////////////////////////////////////////////////////////////////////////////
/* Calculate the average upper limit of an experiment with no signal 
   and expected Poissonian background */

#include <TFeldmanCousins.h>
#include <TMath.h>

using namespace std;

Double_t MeanUL(Double_t Nbackground)
{
  Double_t CL=0.9; // set confidence level
  TFeldmanCousins FC(CL);
  FC.SetMuMax(100.0); // maximum value of signal to calculate the tables! 
                      // increase it for greater values of Nbackground!
  FC.SetMuStep(0.1);
  Double_t PoisProb=0., mu90=0., ul=0.;//, ll=0.; 
  
  //for (Int_t Nobserved=0; Nobserved<200; Nobserved++) // decrease this value for Nbackground<20
  //for (Int_t Nobserved=0; Nobserved<100; Nobserved++) // decrease this value for Nbackground<20
  /* * * * * * * * * * NON FUNZIONA CORRETTAMENTE PER N>20 * * * * * * * * * */
  for (Int_t Nobserved=0; Nobserved<200; Nobserved++) // increase this value for Nbackground>20
    {
      ul = FC.CalculateUpperLimit(Nobserved, Nbackground);
      //ll = FC.GetLowerLimit();
      PoisProb = TMath::Poisson(Nobserved,Nbackground);
      
      mu90 = mu90 + (ul * PoisProb);
      //      cout << "\t" << Nobserved << "\t" << mu90 << endl; // DEBUG!
    }
  
  return mu90;
}

/////////////////////////////////////////////////////////////////////////////////



//Int_t MRF(const TH1F& signal, const TH1F& background) {
//Int_t MRF_calc(TH1* signal, TH1* background) {
Int_t MRF_calc() {

  /*  
      TH1F *astrorepet = new TH1F("astrorepet","cosmic signal",80,6,10);
      astrorepet->GetXaxis()->SetTitle("log_{10}#rho");
      astrorepet->GetYaxis()->SetTitle("events above log_{10}#rho");
      astrorepet->SetLineColor(2);
      astrorepet->SetFillStyle(3004);
      astrorepet->SetFillColor(kRed-7);
      //astrorepet->GetYaxis()->SetRangeUser(1e-3,1e3);
      //astrorepet->GetXaxis()->SetRangeUser(2,8);
      TH1F *atmorepet = new TH1F("atmorepet","atmospheric neutrinos",80,6,10);
      atmorepet->SetLineStyle(7);
      atmorepet->SetLineWidth(2);
      //TH1F *promrepet = new TH1F("promrepet","prompt neutrinos",300,1,4);
      //promrepet->SetLineStyle(7);
      //promrepet->SetLineWidth(2);
      //promrepet->SetLineColor(4);
  

      Double_t bin_sig=0.,  bin_bkg=0.;
      Int_t Nbins = signal->GetNbinsX();

      for (Int_t i=Nbins; i>=0; i--) {
      bin_sig = signal->Integral(i,Nbins,"");
      cout << bin_sig << endl;
      astrorepet->SetBinContent(i,bin_sig);
      }
  
      Nbins = background->GetNbinsX();
      for (Int_t i=Nbins; i>=0; i--) {
      bin_bkg = background->Integral(i,Nbins,"");
      cout << bin_bkg << endl;
      atmorepet->SetBinContent(i,bin_bkg);
      }
    
  */
  
  TChain *snueCC = new TChain("t1");
  TChain *snueNC = new TChain("t1");
  TChain *snumuCC = new TChain("t1");
  TChain *snumuNC = new TChain("t1");
  
  TChain *bn = new TChain("t1");
  TChain *bm = new TChain("t1");

  snueCC->Add("input/*nueCC*.root");
  snueNC->Add("input/*nueNC*.root");
  snumuCC->Add("input/*numuCC*.root");
  snumuNC->Add("input/*numuNC*.root");
  
  bn->Add("input/*nu*.root");
  bm->Add("input/mupage.root");

  TH1D *hsnueCC = new TH1D("hsnueCC", "h", 30, 3, 6);
  TH1D *hsnueNC = new TH1D("hsnueNC", "h", 30, 3, 6);
  TH1D *hsnumuCC = new TH1D("hsnumuCC", "h", 30, 3, 6);
  TH1D *hsnumuNC = new TH1D("hsnumuNC", "h", 30, 3, 6);
  TH1D *hbn = new TH1D("hbn", "h", 30, 3, 6);
  TH1D *hbm = new TH1D("hbm", "h", 30, 3, 6);

  new TCanvas;

  //  *(w2*(4.8*(1e-7))*(pow(nu_e,-2.3))*(1e4)*(0.5))

  //  *w2*(pow(nu_e, -2.))*0.5e-8*1e4")
  snueCC->Draw("log10(t_en)>>hsnueCC", "1.3*(bdt_c>0.39)*(w2*(4.8*(1e-7))*(pow(nu_e,-2.3))*(1e4)*(0.5))");
  snueNC->Draw("log10(t_en)>>hsnueNC", "2.5*(bdt_c>0.39)*(w2*(4.8*(1e-7))*(pow(nu_e,-2.3))*(1e4)*(0.5))");
  snumuCC->Draw("log10(t_en)>>hsnumuCC", "1.1*(bdt_c>0.39)*(w2*(4.8*(1e-7))*(pow(nu_e,-2.3))*(1e4)*(0.5))");
  snumuNC->Draw("log10(t_en)>>hsnumuNC", "(bdt_c>0.39)*(w2*(4.8*(1e-7))*(pow(nu_e,-2.3))*(1e4)*(0.5))");
  
  bn->Draw("log10(t_en)>>hbn", "(bdt_c>0.39)*w_atmo", "same");

  //  bm->Draw("log10(t_en)>>hbm", "(bdt_c>0.39)*w_atmo", "same");


  double p0 =  1.25144e+04;
  double p1 = -3.66461e-01;
  double p2 = -1.43237e-01;
  double p3 = 2.72018e+00;
  double p4 = 5.62589e-01;

  TF2 *f2 = new TF2("f2","[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])",0, 1, 3, 7);
  f2->SetParameters(p0,p1,p2,p3,p4);
  TH2D *hm2d = new TH2D("hm2d", "hm2d", 100, 0, 1, 30, 3, 6);
  hm2d->FillRandom("f2", 10000000);
  hm2d->Scale(p0/10000000);
  
  //4.1 10-6

  Float_t ca=0, cc=0;
  TH1D *atmorepet = new TH1D("atmorepet", "ch", 30, 3, 6);
  TH1D *astrorepet = new TH1D("astrorepet", "ch", 30, 3, 6);

  for(Int_t i = 30; i > 0; i--)
    {
      ca += hbn->GetBinContent(i);

      for (int j = 39; j < 100; j++) {
	ca += hm2d->GetBinContent(j, i);
      }

      
      cc += hsnueCC->GetBinContent(i);
      cc += hsnueNC->GetBinContent(i);
      cc += hsnumuCC->GetBinContent(i);
      cc += hsnumuNC->GetBinContent(i);

      //      out << cc << ' ' << ca << endl;

      atmorepet->SetBinContent(i, ca);
      astrorepet->SetBinContent(i, cc);
      
    }

  //  cout << atmorepet->GetBinContent(1) << endl;
  //  cout << astrorepet->GetBinContent(1) << endl;

  TCanvas *ccum = new TCanvas("ccum","ccum");
  ccum->SetLogy();
  astrorepet->Draw();
  atmorepet->Draw("same");


  TLegend *leg = new TLegend(0.7,0.75,0.9,0.9,NULL,"brNDC");  
  TLegendEntry *entry=leg->AddEntry("astrorepet","Signal #nu","L"); 
  entry=leg->AddEntry("atmorepet","Atms #nu","L"); 
  leg->Draw();

  
  
  ////////////////////////////////////////////////////////////////////////////////////
  //                         MODEL REJECTON FACTOR                                  //
  ////////////////////////////////////////////////////////////////////////////////////
  
  TH1F *hMRF = new TH1F("hMRF","Model Rejection Factor", 30, 3, 6);

  Int_t Nbins = atmorepet->GetNbinsX();
  Double_t bin_sig=0.,  bin_bg=0., mu_90=0., MRF=0., ii=0.;

  //for (Int_t i=1; i<=Nbins; i++) {
  //for (Int_t i=7; i<=15; i++) {
  for (Int_t i=1; i<=30; i++) {

    bin_sig = astrorepet->GetBinContent(i);
    bin_bg = atmorepet->GetBinContent(i);
    if(bin_bg > 0)
      {
	mu_90 = MeanUL(bin_bg);
	MRF = mu_90/bin_sig;
	
	hMRF->SetBinContent(i,MRF);

	ii=i;
	cout << ii << ' ' <<  bin_bg << ' ' <<  bin_sig << ' ' <<  mu_90 << ' ' << MRF << endl;
      }
  }

  TCanvas *cMRF = new TCanvas("cMRF","cMRF");
  cMRF->SetLogy();
  hMRF->Draw();
  


  return 1;
}
