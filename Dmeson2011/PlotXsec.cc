#include "../Bmeson/Bplusdsigmadpt_1ptbins.h"
#include "iostream"
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
using namespace std;
void PlotXsec(){
    TString input = "D0_PbPb_acc_eff_ptbin_14_ybin_6_prompt_FONLLweight_cent-0to100_dataptshape_y1_Ncollweight1.root";
    TFile *inf = new TFile(input.Data());
  	TH1D* d0accxeff_pt = (TH1D*)inf->Get("d0accxeff_pt");
	//for(int i = 1; i <= d0accxeff_pt->GetNbinsX(); i++){  
	//	cout<<d0accxeff_pt->GetBinLowEdge(i)<<endl;
	//	cout<<d0accxeff_pt->GetBinLowEdge(i)+d0accxeff_pt->GetBinWidth(i)<<endl;;
	//	cout<<"---"<<endl;
	//}
	//return;
	string fonlldata="../fonllData/D_pp_pt_rap1_2p75_pt0to200.dat";
	double diffXsec = 0;
    double BRchain=3.93e-2+1.399e-4;
    double Fraction=0.557;
	
	const int BINS = 14;
	double apt[BINS];
	double asigma[BINS];
	double aptl[BINS];
	double aerrorl[BINS];
	double aerrorh[BINS];
	double aerrorl2[BINS];
	double aerrorh2[BINS];
	double aerrorl3[BINS];
	double aerrorh3[BINS];
	double genB[BINS];
	double effB[BINS];
	double recoB[BINS];
    double bins[BINS+1] = {1.5, 2.5, 3.5, 4.5, 5.5, 7, 9, 11, 13, 16, 20, 28, 40, 60, 100};
	double effD[BINS] = {0.0004110786, 0.004868644, 0.01467014, 0.03677896, 0.07008137, 0.1185297, 0.1670131, 0.2038537, 0.2419492, 0.2935908, 0.3535325, 0.4175788, 0.4175788, 0.4175788};
	double effDerr[BINS] = {8.195642e-05, 0.0004130658 , 0.0009189658, 0.001791014, 0.002483458, 0.003386794, 0.004592759, 0.005642633, 0.005676826, 0.006227245, 0.00599039, 0.007349736, 0.007349736, 0.007349736};
	double RAA[BINS] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
	for(int i = 0; i < BINS; i++){
		Bplusdsigmadpt_1ptbins(false, diffXsec, bins[i], bins[i+1], fonlldata);
		cout<<"diffXsec: "<<diffXsec<<endl;
		apt[i] = (bins[i]+bins[i+1])/2;
		asigma[i] = diffXsec;
		aptl[i] = (bins[i+1]-bins[i])/2;
		aerrorl[i] = 0;
		aerrorh[i] = 0;

		//genB[i] = 0.5 * 2 * 166e-6*10 * 208*208 * diffXsec * (bins[i+1]-bins[i]);
		genB[i] = RAA[i] * 2 * 166e-6* 208*208 * diffXsec * (bins[i+1]-bins[i]);
		aerrorl2[i] = 0;
		aerrorh2[i] = 0;
		cout<<"genB: "<<genB[i]<<endl;
		
        recoB[i] = genB[i]*effD[i]*BRchain*Fraction;	
		aerrorl3[i] = recoB[i]*effDerr[i];
		aerrorh3[i] = recoB[i]*effDerr[i];
		cout<<"aerrorl3: "<<aerrorl3[i]<<endl;
		cout<<"recoB: "<<recoB[i]<<endl;
	}
    TCanvas* cr = new TCanvas("cr","cr",600,500);
    gStyle->SetOptTitle(1);
    cr->SetLogy();

    TH2F* hempty=new TH2F("hempty","",1,0,100,10.,10,10000000);
	hempty->SetTitle("average B meson dsigma/dpt");
    hempty->Draw();
    TGraphAsymmErrors* dXsec = new TGraphAsymmErrors(BINS, apt, asigma, aptl, aptl, aerrorl, aerrorh);
	dXsec->SetLineWidth(3);
	dXsec->SetLineColor(4);
	dXsec->Draw("p same");

    TH2F* hempty2=new TH2F("hempty2","",1,0,100,10.,1e5,1e11);
	hempty2->SetTitle("# of produced D0 meson");
	hempty2->GetXaxis()->SetTitle("D0 pt");
    hempty2->Draw();
    TGraphAsymmErrors* GenB = new TGraphAsymmErrors(BINS, apt, genB, aptl, aptl, aerrorl2, aerrorh2);
	GenB->SetLineWidth(3);
	GenB->SetLineColor(4);
	GenB->Draw("p same");
	cr->SaveAs("Plots/GenD.png");

    TH2F* hempty3=new TH2F("hempty3","",1,0,100,10.,1e3,1e8);
	hempty3->SetTitle("# of reconstructed D0");
	hempty3->GetXaxis()->SetTitle("D0 pt");
    hempty3->Draw();
    TGraphAsymmErrors* RecoB = new TGraphAsymmErrors(BINS, apt, recoB, aptl, aptl, aerrorl3, aerrorh3);
	RecoB->SetLineWidth(3);
	RecoB->SetLineColor(4);
	RecoB->Draw("p same");
	cr->SaveAs("Plots/RecoD.png");
}
