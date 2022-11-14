#include "particletype.h"
#include "resonancetype.h"
#include "particle.h"

#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TMath.h"
#include "TCanvas.h"


void analysis(){
    TFile* file = new TFile("particles.root");
    file -> ls();
    TH1F* hParticleType = (TH1F*) file -> Get("hParticleType");
    TH1F* hAzimuthAngle = (TH1F*) file -> Get("hAzimuthAngle");
    TH1F* hPolarAngle = (TH1F*) file -> Get("hPolarAngle");
    TH1F* hMomentum = (TH1F*) file -> Get("hMomentum");
    TH1F* hTransverseMomentum = (TH1F*) file -> Get("hTransverseMomentum");
    TH1F* hEnergy = (TH1F*) file -> Get("hEnergy");
    TH1F* hInvMass = (TH1F*) file -> Get("hInvMass");
    TH1F* hInvMassDisc = (TH1F*) file -> Get("hInvMassDisc");
    TH1F* hInvMassConc = (TH1F*) file -> Get("hInvMassConc");
    TH1F* hInvMassPkd = (TH1F*) file -> Get("hInvMassPkd");
    TH1F* hInvMassPkc = (TH1F*) file -> Get("hInvMassPkc");
    TH1F* hInvMassResonance = (TH1F*) file -> Get("hInvMassResonance");
    
    //number of particles
    TString s = "i";
    for(int i=1; i<8; ++i){
        std::cout << "bin " << s+i << " content of types of particle: " << hParticleType -> GetBinContent(i) << " +/- " << hParticleType -> GetBinError(i) << '\n';
    }

    //angles
    TF1* fAzimuthAngle = new TF1("f1", "[0]", 0, 2*TMath::Pi());
    fAzimuthAngle -> SetParameter(0,130); //aggiungere valore per il parametro
    hAzimuthAngle -> Fit(fAzimuthAngle);
    std::cout << "Azimuth angle: " << '\n';
    std::cout << "f(x) = " << fAzimuthAngle->GetParameter(0) << '\n';
    std::cout << "chisquare/NDF: " << fAzimuthAngle->GetChisquare() << " / " << fAzimuthAngle->GetNDF() << '\n';
    std::cout << "probability: " << fAzimuthAngle->GetProb() << '\n';
    std::cout << "___________________________" << '\n';

    TF1* fPolarAngle = new TF1("f2", "[0]", 0, TMath::Pi());
    fAzimuthAngle -> SetParameter(0,130); //aggiungere valore per il parametro
    hPolarAngle -> Fit(fPolarAngle);
    std::cout << "Polar angle: " << '\n';
    std::cout << "f(x) = " << fPolarAngle->GetParameter(0) << '\n';
    std::cout << "chisquare/NDF: " << fPolarAngle->GetChisquare() << " / " << fPolarAngle->GetNDF() << '\n';
    std::cout << "probability: " << fPolarAngle->GetProb() << '\n';
    std::cout << "___________________________" << '\n';


    //momentum
    TF1* fMom = new TF1("f3", "exp(-x/[0])", 0,5); //????
    fMom -> SetParameter(0,1); //avevamo fissato 1 come media
    hMomentum -> Fit(fMom);
    std::cout << "Momentum: " << '\n';
    std::cout << "mean: " << fMom->GetParameter(0) << '\n';
    std::cout << "chisquare/NDF: " << fMom->GetChisquare() << " / " << fMom->GetNDF() << '\n';
    std::cout << "probability: " << fMom->GetProb() << '\n';
    std::cout << "___________________________" << '\n';


    //InvMass
    TH1F* hres1 = new TH1F("hres1", "division 1,2", 1E5, 0, 5);
    hres1 -> Add(hInvMassDisc, hInvMassConc, 1,-1);
    
    TH1F* hres2 = new TH1F("hres2", "division 3,4", 1E5, 0, 5);
    hres2 -> Add(hInvMassPkd, hInvMassPkc, 1,-1);

    //fitting
    TF1* fgaus1 = new TF1("f4", "[0]*exp(-(x-[1])^2/2*[2]^2)", 0,5);
    fgaus1 -> SetParameter(0,10); //impostare
    fgaus1 -> SetParameter(1,0.89166);
    fgaus1 -> SetParameter(2,0.050);

    hres1 -> Fit(fgaus1);
    std::cout << "hres1: " << '\n';
    std::cout << "maximum: " << fgaus1->GetParameter(0) << '\n';
    std::cout << "mean: " << fgaus1->GetParameter(1) << '\n';
    std::cout << "sigma: " << fgaus1->GetParameter(2) << '\n';
    std::cout << "chisquare/NDF: " << fgaus1->GetChisquare() << " / " << fgaus1->GetNDF() << '\n';
    std::cout << "probability: " << fgaus1->GetProb() << '\n';
    std::cout << "___________________________" << '\n';

    TF1* fgaus2 = new TF1("f4", "[0]*exp(-(x-[1])^2/2*[2]^2)", 0,5);
    fgaus2 -> SetParameter(0,4000); //impostare
    fgaus2 -> SetParameter(1,0.89166);
    fgaus2 -> SetParameter(2,0.050);
    hres2 -> Fit(fgaus2);
    std::cout << "hres2: " << '\n';
    std::cout << "maximum: " << fgaus2->GetParameter(0) << '\n';
    std::cout << "mean: " << fgaus2->GetParameter(1) << '\n';
    std::cout << "sigma: " << fgaus2->GetParameter(2) << '\n';
    std::cout << "chisquare/NDF: " << fgaus2->GetChisquare() << " / " << fgaus2->GetNDF() << '\n';
    std::cout << "probability: " << fgaus2->GetProb() << '\n';
    std::cout << "___________________________" << '\n';



    //drawing histos of InvMass
    TCanvas* cRes = new TCanvas("cRes", "InvMass", 600,600);
    cRes -> Divide(2,1);

    cRes -> cd(1);
    hres1 -> Draw("HIST");
    fgaus1 -> Draw("SAME");

    cRes -> cd(2);
    hres2 -> Draw("HIST");
    fgaus2 -> Draw("SAME");


    //adding cosmetics and printing
    hres1 -> GetYaxis() -> SetTitle("entries");
    hres1 -> GetXaxis() -> SetTitle("InvMass");
    //hres1 -> GetYaxis() -> SetRangeUser();
    hres1 -> SetFillColor(kBlue);

    fgaus1 -> SetLineColor(kRed);

    hres2 -> GetYaxis() -> SetTitle("entries");
    hres2 -> GetXaxis() -> SetTitle("InvMass");
    //hres2 -> GetYaxis() -> SetRangeUser();
    hres2 -> SetFillColor(kGreen);

    fgaus2 -> SetLineColor(kRed);



    file -> Close();
}