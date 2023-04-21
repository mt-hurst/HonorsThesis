//Script written 11/30/22 by Tristan Hurst
//
//Plotting individual detector responses of 
//Showermax and pion detector to demonstrate
//distinguishability of pions from electrons

#include <TF1.h>
void shmx_pion_separate(const TString& files)
{
TChain* P = new TChain("T"); //tchain for pion
TChain* E = new TChain("T"); //tchain for moll


if (files == "moll" || files == "pion" || files == "both"){
    if (files == "moll" || files == "both"){
        E->Add("/w/halla-scshelf2102/moller12gev/mthurst/remoll/rootfiles/small_moll_noblock.root");//symbolic link to moller files
    }
    if (files == "pion" || files == "both"){
        P->Add("/w/halla-scshelf2102/moller12gev/mthurst/remoll/rootfiles/small_pion_noblock.root");//symbolic link to pion files
    }
}
if (files != "moll" && files != "pion" && files != "both"){
    P->Add(files);
}

Double_t rate_e = 0;
std::vector<remollEventParticle_t>* parts_e = 0;
std::vector<remollGenericDetectorHit_t>* hits_e = 0;

Double_t rate_p = 0;
std::vector<remollEventParticle_t>* parts_p = 0;
std::vector<remollGenericDetectorHit_t>* hits_p = 0;

// Define some branches of the Tree (which is "T")
E->SetBranchAddress("rate", &rate_e);
E->SetBranchAddress("hit", &hits_e);
E->SetBranchAddress("part", &parts_e);

P->SetBranchAddress("rate", &rate_p);
P->SetBranchAddress("hit", &hits_p);
P->SetBranchAddress("part", &parts_p);
//Start historgram structure 
TH1F* hr29_e = new TH1F("hr29_e","Moll Electron hits of det 29",200,1080,1180);
hr29_e->GetXaxis()->SetTitle("Radial Dist(mm)");


TH1F* hshmx_moll_rate = new TH1F("hshmx_moll_rate","Moll Shmx Rate", 100,0,200);
hshmx_moll_rate->GetXaxis()->SetTitle("SHMX Layer Hits");

TH1F* hshmx_pion_rate = new TH1F("hshmx_pion_rate","Pion Shmx Rate", 100,0,200);
hshmx_pion_rate->GetXaxis()->SetTitle("SHMX Layer Hits");

TH1F* hpd_moll_rate = new TH1F("hpd_moll_rate","Moll PD Rate", 100,0,2000);
hpd_moll_rate->GetXaxis()->SetTitle("PD Opt Hits");

TH1F* hpd_pion_rate = new TH1F("hpd_pion_rate","Pion PD Rate", 100,0,2000);
hpd_pion_rate->GetXaxis()->SetTitle("PD Opt Hits");


TH2F* hshmx_e_pion = new TH2F("hshmx_e_pion","SHMX Layer hits vs Opt Photon Hits",50,0,200,50,0,2000);
hshmx_e_pion->GetXaxis()->SetTitle("SHMX Layer Hits");hshmx_e_pion->GetYaxis()->SetTitle("Pion PMT Opt Hits");

TH2F* hshmx_e_pion_trig = new TH2F("hshmx_e_pion_trig","SHMX Layer hits vs Opt Photon Hits (w/ Trig)",50,0,200,50,0,2000);
hshmx_e_pion_trig->GetXaxis()->SetTitle("SHMX Layer Hits");hshmx_e_pion_trig->GetYaxis()->SetTitle("Pion PMT Opt Hits");

TH2F* hshmx_e_pion_rate = new TH2F("hshmx_e_pion_rate","SHMX Layer hits vs Opt Photon Hits (Rate)",50,0,200,50,0,2000);
hshmx_e_pion_rate->GetXaxis()->SetTitle("SHMX Layer Hits [Hz/mA]");hshmx_e_pion_rate->GetYaxis()->SetTitle("Pion PMT Opt Hits [Hz/mA]");

TH2F* hshmx_e_pion_trig_rate = new TH2F("hshmx_e_pion_trig_rate","SHMX Layer hits vs Opt Photon Hits (w/ Trig)",50,0,200,50,0,2000);
hshmx_e_pion_trig_rate->GetXaxis()->SetTitle("SHMX Layer Hits[Hz/mA]");hshmx_e_pion_trig_rate->GetYaxis()->SetTitle("Pion PMT Opt Hits [Hz/mA]");


TH2F* hshmx_p_pion = new TH2F("hshmx_p_pion","SHMX Layer hits vs Opt Photon Hits",50,0,200,50,0,2000);
hshmx_p_pion->GetXaxis()->SetTitle("SHMX Layer Hits");hshmx_p_pion->GetYaxis()->SetTitle("Pion PMT Opt Hits");

TH2F* hshmx_p_pion_trig = new TH2F("hshmx_p_pion_trig","SHMX Layer hits vs Opt Photon Hits (w/ Trig)",50,0,200,50,0,2000);
hshmx_p_pion_trig->GetXaxis()->SetTitle("SHMX Layer Hits");hshmx_p_pion_trig->GetYaxis()->SetTitle("Pion PMT Opt Hits");

TH2F* hshmx_p_pion_rate = new TH2F("hshmx_p_pion_rate","SHMX Layer hits vs Opt Photon Hits (Rate)",50,0,200,50,0,2000);
hshmx_p_pion_rate->GetXaxis()->SetTitle("SHMX Layer Hits [Hz/mA]");hshmx_p_pion_rate->GetYaxis()->SetTitle("Pion PMT Opt Hits [Hz/mA]");

TH2F* hshmx_p_pion_trig_rate = new TH2F("hshmx_p_pion_trig_rate","SHMX Layer hits vs Opt Photon Hits (w/ Trig)",50,0,200,50,0,2000);
hshmx_p_pion_trig_rate->GetXaxis()->SetTitle("SHMX Layer Hits[Hz/mA]");hshmx_p_pion_trig_rate->GetYaxis()->SetTitle("Pion PMT Opt Hits [Hz/mA]");

//End Histogram structure

//Define constants and flags
bool opthit = false; // flag for if optical photon hit the PMT of Pion, 8000
bool trig29 = false; //flag for if particles exit the pion det and hit the slice of 29

const double pi = 3.14159265358979323846;
double mincut = -pi;
double maxcut = pi;
//begin moll integers

int q_shmx_e = 0; //# of charged particle hits for each event in SHMX 
int pmt_opt_e = 0; // # of opt photons on pion PMT
int q_shmx_trig_e = 0; //same as q_shmx but with trigger scintillator slice hit
int pmt_opt_trig_e = 0; //same as pmt_opt_trig but with trigger scintillator slice hit
//end moll event integers

//pion integers
int q_shmx_p = 0; //# of charged particle hits for each event in SHMX
int pmt_opt_p = 0; // # of opt photons on pion PMT
int q_shmx_trig_p = 0; //same as q_shmx but with trigger scintillator slice hit
int pmt_opt_trig_p = 0; //same as pmt_opt_trig but with trigger scintillator slice hit
//end pion event integers

//Begin trigger cuts 
double minrad = 1110;//default 1110
double maxrad = 1150;//default 1150
//End Trigger cuts

double ratio_pion_moll_rate = 0;
double pion_rate = 1; //filled with dummy value
double moll_rate = 0.5; //filled with dummy value to prevent divide by 0 error

double ratio_pion_moll_rate_mollregion = 0;
double pion_rate_mollregion =1;
double moll_rate_mollregion =1;

double ratio_pion_moll_rate_all = 0;
double pion_rate_all =1;
double moll_rate_all =1;


//End defining constants and flags

//Begin Loop over different radial extents of Trigger Scintillator
int n = 1;
//initialize 4 data sets, first is the width,
//the other three are for holding the values of the rate ratios in
//different detector response regions
double_t width[n], ratio_vals[n],ratio_vals_moll[n],ratio_vals_all[n];


for (int i = 0 ; i<n ; i++){
/*minrad = 1125-i*5 ;
maxrad = 1135+i*5 ;*/
width[i] = maxrad - minrad;


//begin Moll event processing
if (files == "moll" || files =="both"){
    //Loop over all Moller Events
    //E->GetEntries()
 for (size_t iev = 0; iev < E->GetEntries(); iev++){
    E->GetEntry(iev);
    trig29 = false;
    opthit = false;
    q_shmx_e = 0;
    q_shmx_trig_e = 0;
    pmt_opt_e = 0;
    pmt_opt_trig_e = 0;
    //first loop over hits, defining flags
    for (size_t ihit = 0; ihit< hits_e->size(); ihit++){
        remollGenericDetectorHit_t hit = hits_e->at(ihit);
     if (hit.det == 29 && (hit.r > minrad && hit.r <maxrad) && (hit.pid == 11 ||hit.pid ==-11 ||hit.pid ==211 ||hit.pid ==-211 ||hit.pid == 13 ||hit.pid == -13 ||hit.pid == 2212)){
        trig29 = true;
     }
     if (hit.det == 8000 && hit.pid == 0){
        opthit = true;
     }
    }//end first loop over hits

    //second loop over hits
    for (size_t ihit = 0; ihit<hits_e->size();ihit++){
        remollGenericDetectorHit_t hit = hits_e->at(ihit);
/*
        if (opthit && (hit.det == 29) && (hit.pid ==11 || hit.pid == -11)){
            hr29_e->Fill(hit.r,rate_e);
        }
*/

        if (opthit && hit.pid == 0){
            pmt_opt_e = pmt_opt_e + 1;
        }
        if (trig29 && opthit && hit.pid == 0){
            pmt_opt_trig_e = pmt_opt_trig_e + 1;
        }
        if ( hit.det > 7001 && hit.det < 7279 && hit.pid != 22 && hit.pid != 0){
            q_shmx_e = q_shmx_e + 1;
        }
        if (trig29 && hit.det > 7001 && hit.det < 7279 && hit.pid != 22 && hit.pid != 0 ){
            q_shmx_trig_e = q_shmx_trig_e + 1;
        }

    }//end second loop over hits 
    if (q_shmx_e != 0 && pmt_opt_e != 0){
    hshmx_e_pion-> Fill(q_shmx_e,pmt_opt_e);
    hshmx_e_pion_rate-> Fill(q_shmx_e,pmt_opt_e,rate_e);
    hshmx_moll_rate->Fill(q_shmx_e,rate_e);
    hpd_moll_rate->Fill(pmt_opt_e,rate_e);
    }
    if (q_shmx_trig_e != 0 && pmt_opt_trig_e != 0){
    hshmx_e_pion_trig->Fill(q_shmx_trig_e,pmt_opt_trig_e);
    hshmx_e_pion_trig_rate->Fill(q_shmx_trig_e,pmt_opt_trig_e,rate_e);
    }
        if ((iev%(E->GetEntries()/10)) ==0){cout<<"Moll Events Finished : "<<round((double)iev*100/(double)E->GetEntries())<<"\%"<<endl;}

    }//end loop over moll events
}//end moll event processing
    moll_rate = hshmx_e_pion_trig_rate->Integral(0,15,0,40);
    moll_rate_mollregion= hshmx_e_pion_trig_rate->Integral(15,100,0,20);
    moll_rate_all = hshmx_e_pion_trig_rate->Integral(0,100,0,100);


//begin Pion event processing

if (files == "pion" || files =="both"){
    //Loop over all Pion Events
 for (size_t iev = 0; iev < P->GetEntries(); iev++){
    P->GetEntry(iev);
    trig29 = false;
    opthit = false;

    q_shmx_p = 0;
    q_shmx_trig_p = 0;
    pmt_opt_p = 0;
    pmt_opt_trig_p = 0;

    //first loop over hits, defining flags
    for (size_t ihit = 0; ihit< hits_p->size(); ihit++){
        remollGenericDetectorHit_t hit = hits_p->at(ihit);
     if (hit.det == 29 && (hit.r > minrad && hit.r <maxrad) && (hit.pid == 11 ||hit.pid ==-11 ||hit.pid ==211 ||hit.pid ==-211 ||hit.pid == 13 ||hit.pid == -13 ||hit.pid == 2212)){
        trig29 = true;
     }
     if (hit.det == 8000 && hit.pid == 0){
        opthit = true;
     }
    }//end first loop over hits

    //second loop over hits
    for (size_t ihit = 0; ihit<hits_p->size();ihit++){
        remollGenericDetectorHit_t hit = hits_p->at(ihit);
/*
        if (opthit && (hit.det == 29) && (hit.pid ==11 || hit.pid == -11)){
            hr29_e->Fill(hit.r,rate_p);
        }
*/

        if (opthit && hit.pid == 0){
            pmt_opt_p = pmt_opt_p + 1;
        }
        if (trig29 && opthit && hit.pid == 0){
            pmt_opt_trig_p = pmt_opt_trig_p + 1;
        }
        if ( hit.det > 7001 && hit.det < 7279 && hit.pid != 22 && hit.pid != 0){
            q_shmx_p = q_shmx_p + 1;
        }
        if (trig29 && hit.det > 7001 && hit.det < 7279 && hit.pid != 22 && hit.pid != 0 ){
            q_shmx_trig_p = q_shmx_trig_p + 1;
        }

    }//end second loop over hits 
    if (q_shmx_p != 0 && pmt_opt_p != 0){
    hshmx_p_pion-> Fill(q_shmx_p,pmt_opt_p);
    hshmx_p_pion_rate-> Fill(q_shmx_p,pmt_opt_p,rate_p);
    hshmx_pion_rate->Fill(q_shmx_p,rate_p);
    hpd_pion_rate->Fill(pmt_opt_p,rate_p);
    }
    if (q_shmx_trig_p != 0 && pmt_opt_trig_p != 0){
    hshmx_p_pion_trig->Fill(q_shmx_trig_p,pmt_opt_trig_p);
    hshmx_p_pion_trig_rate->Fill(q_shmx_trig_p,pmt_opt_trig_p,rate_p);
    }
    if (iev%(P->GetEntries()/10) ==0){cout<<"Moll Events Finished : "<<round((double)iev*100/(double)P->GetEntries())<<"\%"<<endl;
    }//end loop over Pion events
    pion_rate = hshmx_p_pion_trig_rate->Integral(0,15,0,40);
    pion_rate_mollregion= hshmx_p_pion_trig_rate->Integral(15,100,0,20);
    pion_rate_all = hshmx_p_pion_trig_rate->Integral(0,100,0,100);

    }

}//end moll event processing
ratio_pion_moll_rate = (double)pion_rate /(double)moll_rate;
ratio_pion_moll_rate_mollregion = (double)pion_rate_mollregion/(double)moll_rate_mollregion;
ratio_pion_moll_rate_all = (double)pion_rate_all/(double)moll_rate_all;

cout<<"Annulus Width : "<< maxrad - minrad << "mm"<<endl;
cout<< "Pion rate (Pion Region): "<<pion_rate<<"   Moll rate (Pion Region): "<< moll_rate << endl;
cout<<"Ratio of Pion Rate/Moll Rate : "<< ratio_pion_moll_rate<<endl;

cout<< "Pion rate (Moll Region): "<<pion_rate_mollregion<<"   Moll rate (Moll Region): "<< moll_rate_mollregion << endl;
cout<<"Ratio of Pion Rate/Moll Rate : "<< ratio_pion_moll_rate_mollregion<<endl;

cout<< "Pion rate (All Region): "<<pion_rate_all<<"   Moll rate (all Region): "<< moll_rate_all << endl;
cout<<"Ratio of Pion Rate/Moll Rate : "<< ratio_pion_moll_rate_all<<endl;

ratio_vals[i]= ratio_pion_moll_rate;
ratio_vals_moll[i]= ratio_pion_moll_rate_mollregion;
ratio_vals_all[i] = ratio_pion_moll_rate_all;
hshmx_e_pion_trig_rate->Reset("ICESM");
hshmx_p_pion_trig_rate->Reset("ICESM");
}
TGraph *gr1 = new TGraph (n, width, ratio_vals);
gr1->GetXaxis()->SetTitle("Width of Annulus (mm)");
gr1->SetTitle("Pion Sensitive Region: Rate vs Width");
gr1->GetYaxis()->SetTitle("Ratio of Mu/e");

TGraph *gr2 = new TGraph (n,width,ratio_vals_moll);
gr2->GetXaxis()->SetTitle("Width of Annulus (mm)");
gr2->SetTitle("Moll Sensitive Region: Rate vs Width");
gr2->GetYaxis()->SetTitle("Ratio of Mu/e");


TGraph *gr3 = new TGraph (n,width,ratio_vals_all);
gr3->GetXaxis()->SetTitle("Width of Annulus (mm)");
gr3->SetTitle("Total Region: Rate vs Width");
gr3->GetYaxis()->SetTitle("Ratio of Mu/e");

TCanvas *c1 = new TCanvas("c1","Annulus Width vs Ratio of Mu/e");
c1->Divide(1,3);
c1->cd(1);
gr1->Draw("AC*");
c1->cd(2);
gr2->Draw("AC*");
c1->cd(3);
gr3->Draw("AC*");
c1->cd();



//c1->SaveAs("/home/mthurst/SHMX_Layers/Scint_Images/PiMu_rate_vs_width.pdf");
TCanvas *p1 = new TCanvas("p1", "Overlaid Graphs SHMX");
p1->Divide(1,1);
p1->cd(1);
hshmx_pion_rate->Draw();
hshmx_moll_rate->SetLineColor(kRed);
hshmx_moll_rate->Draw("same");

TCanvas *p2 = new TCanvas("p2", "Overlaid Graphs PD");
p2->Divide(1,1);
p2->cd(1);
hpd_pion_rate->Draw();
hpd_moll_rate->SetLineColor(kRed);
hpd_moll_rate->Draw("same");

/*
TCanvas *c2 = new TCanvas("c2", "Radial Hits on Det 29");
c2->Divide(1,1);
c2->cd(1);
hr29_e->Draw();
*/
}
