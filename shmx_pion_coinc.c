// Script Written by Tristan Hurst  June 16th  2022
//
// Goal: Plotting # of charged particles in SHMX layers vs # of photons seen in pionPMT
//
#include <TF1.h>
void shmx_pion_coinc(const TString& files)
{
 TChain* T = new TChain("T");
 T->Add(files);
/*
 T->Add("shmx_pion.root");
 T->Add("shmx_moll.root");
*/
 Double_t rate = 0;
 std::vector<remollEventParticle_t>* parts = 0;
 std::vector<remollGenericDetectorHit_t>* hits = 0;

//

 int q_shmx = 0; //# of charged particle hits for each event in SHMX
 int pmt_opt = 0; // # of opt photons on pion PMT

 int q_shmx_trig = 0; //same as q_shmx but with trigger scintillator slice hit
 int pmt_opt_trig = 0; //same as pmt_opt_trig but with trigger scintillator slice hit

 // Define some branches of the Tree (which is "T")
 T->SetBranchAddress("rate", &rate);
 T->SetBranchAddress("hit", &hits);
 T->SetBranchAddress("part", &parts);

 //Define some histograms
 // these first are all 1D histograms of floating point quantities (i.e. TH1F type)


 TH1F* h331_e_pz = new TH1F("h331_e_pz","z Momenta of 331 hits",100,0,10000);
 h331_e_pz -> GetXaxis()->SetTitle("Z momenta [MeV]");//using for debugging

 TH1F* h8001_e = new TH1F("h8001_e", "Radius of Det 8001 Electron hits",100,1080,1180);
 h8001_e ->GetXaxis()->SetTitle("Radius of Hit [mm]");

 TH1F* h8001_p = new TH1F("h8001_p", "Radius of Det 8001 Pion hits",100,1080,1180);
 h8001_p ->GetXaxis()->SetTitle("radius of Hit [mm]");

 TH1F* hr29_e = new TH1F("hr29_e","Hit radius at detector plane 29 (pion det)",100,1080,1180);
 hr29_e->GetXaxis()->SetTitle("r [mm], 20 mm / bin"); hr29_e->GetYaxis()->SetTitle("Hits");

 TH1F* hr29_pion = new TH1F("hr29_pion","Hit radius at detector plane 29 (pion det) raw hits",100,1080,1180);
 hr29_pion->GetXaxis()->SetTitle("r [mm], 20 mm / bin"); hr29_pion->GetYaxis()->SetTitle("hits");


 TH1F* h29_ep = new TH1F("h29_ep","All Electrons in Slice",100,0,100);
 h29_ep->GetXaxis()->SetTitle("Momentum [MeV]");

 TH1F* h29_pp = new TH1F("h29_pp","All Pions in Slice",100,0,10000);
 h29_pp->GetXaxis()->SetTitle("Momentum [MeV]");

 TH1F* h29_ep_pmt = new TH1F("h29_ep_pmt","Opthit Electrons in Slice",100,0,100);
 h29_ep_pmt->GetXaxis()->SetTitle("Momentum [MeV]");

 TH1F* h29_pp_pmt = new TH1F("h29_pp_pmt","Opthit Pions in Slice",100,0,10000);
 h29_pp_pmt->GetXaxis()->SetTitle("Momentum [MeV]");



 TH1F* h8001_ep = new TH1F("h8001_ep","All Electrons in 8001",100,0,1000);
 h8001_ep->GetXaxis()->SetTitle("Momentum [MeV]");

 TH1F* h8001_pp = new TH1F("h8001_pp","All Pions in 8001",100,0,10000);
 h8001_pp->GetXaxis()->SetTitle("Momentum [MeV]");

 TH1F* h8001_ep_pmt = new TH1F("h8001_ep_pmt","Opthit Electrons in 8001",100,0,1000);
 h8001_ep_pmt->GetXaxis()->SetTitle("Momentum [MeV]");

 TH1F* h8001_pp_pmt = new TH1F("h8001_pp_pmt","Opthit Pions in 8001",100,0,10000);
 h8001_pp_pmt->GetXaxis()->SetTitle("Momentum [MeV]");


 // these next 2D scatterplots  of floating point quantities (i.e. TH2F type)

 TH2F* h28_xy = new TH2F("h28_xy","Detector 28 charged hits y vs x ",40,-2000,2000,40,-2000,2000);
 TH2F* h29_xy_e = new TH2F("h29_xy_e","Detector 29 Electron hits y vs x ",60,-2000,2000,60,-2000,2000);
 TH2F* h29_xy_p = new TH2F("h29_xy_p","Detector 29 Pion hits y vs x ",60,-2000,2000,60,-2000,2000);

 TH2F* hshmx_pion = new TH2F("hshmx_pion","SHMX Layer hits vs Opt Photon Hits",50,0,200,50,0,2000);
 hshmx_pion->GetXaxis()->SetTitle("SHMX Layer Hits");hshmx_pion->GetYaxis()->SetTitle("Pion PMT Opt Hits");

 TH2F* hshmx_pion_trig = new TH2F("hshmx_pion_trig","SHMX Layer hits vs Opt Photon Hits (w/ Trig)",50,0,200,50,0,2000);
 hshmx_pion_trig->GetXaxis()->SetTitle("SHMX Layer Hits");hshmx_pion_trig->GetYaxis()->SetTitle("Pion PMT Opt Hits");

 TH2F* hshmx_pion_rate = new TH2F("hshmx_pion_rate","SHMX Layer hits vs Opt Photon Hits (Rate)",50,0,100,50,0,1000);
 hshmx_pion_rate->GetXaxis()->SetTitle("SHMX Layer Hits [Hz/mA]");hshmx_pion_rate->GetYaxis()->SetTitle("Pion PMT Opt Hits [Hz/mA]");

 TH2F* hshmx_pion_trig_rate = new TH2F("hshmx_pion_trig_rate","SHMX Layer hits vs Opt Photon Hits (w/ Trig)",50,0,100,50,0,1000);
 hshmx_pion_trig_rate->GetXaxis()->SetTitle("SHMX Layer Hits[Hz/mA]");hshmx_pion_trig_rate->GetYaxis()->SetTitle("Pion PMT Opt Hits [Hz/mA]");


 bool opthit = false; // flag for if optical photon hit the PMT of Pion, 8000
 bool trig29 = false; //flag for if particles exit the pion det and hit the slice of 29
const double pi = 3.14159265358979323846;
const double mincut = -pi;
const double maxcut = pi;
 const double minrad = 1110;//min radius of Dstream trigger cut
 const double maxrad = 1150;//max radius of Dstream trigger cut
 // Loop over all events
 //T-> GetEntries(); use this for a full run
 //iev < 1000-10000 for a quick run
 for (size_t iev = 0; iev < T->GetEntries(); iev++) {
   T->GetEntry(iev);
   trig29 = false;
   opthit = false;
   q_shmx = 0;
   pmt_opt = 0;
   q_shmx_trig = 0;
   pmt_opt_trig = 0;


   // Process hits: loop over all the hits in this event
   for (size_t ihit = 0; ihit < hits->size(); ihit++) {
     remollGenericDetectorHit_t hit = hits->at(ihit);
     // implementing a debugging phi restriction
     if (hit.ph < mincut || hit.ph > maxcut ){continue;}
     //end debugging phi restriction
     // fill histograms weighted by rate variable, and raw (as thrown)
     if (hit.det == 29 && (hit.r > minrad && hit.r <maxrad && (hit.pid == 11 ||hit.pid ==-11 ||hit.pid ==211 ||hit.pid ==-211 ||hit.pid == 13 ||hit.pid == -13 ||hit.pid == 2212))){
       trig29 = true;
     }

     if (hit.det ==8000 && hit.pid == 0){
       opthit = true;//this flags hits in events that generated light on the PMT of the Pion Detector
       // pmt_opt = pmt_opt + 1;
     }

     if (opthit && hit.det == 29 && (hit.pid == 11 || hit.pid == -11)){
       hr29_e-> Fill(hit.r);
     }
     if (opthit && hit.det == 29 && (hit.pid == 211 || hit.pid == -211)){
       hr29_pion-> Fill(hit.r);
     }

   }  // end loop over hits

  // debugging comment
  //cout << iev <<"||event ||" <<pmt_opt << "||pmt hits || "<<q_shmx << "||  charged layer hits"<<endl;

   // Process hits again
   for (size_t ihit = 0; ihit < hits->size(); ihit++) {
    remollGenericDetectorHit_t hit = hits->at(ihit);
    // implementing a debugging phi restriction
    if (hit.ph < mincut || hit.ph > maxcut ){continue;}
    //end debugging phi restriction
    if (opthit && (hit.r > minrad && hit.r <maxrad) ) {
      if (hit.det == 29 &&(hit.pid == 11||hit.pid==-11))  {h29_xy_e->Fill(hit.x,hit.y);}
      if (hit.det == 29 &&((hit.pid == 211 || hit.pid == -211)||(hit.pid == 13 ||hit.pid ==-13))) {h29_xy_p->Fill(hit.x,hit.y);}
    }

    if (opthit && hit.pid == 0){
      pmt_opt = pmt_opt + 1;
    }
    if (trig29 && opthit && hit.pid == 0){
      pmt_opt_trig = pmt_opt_trig + 1;
    }
    if ( hit.det > 7001 && hit.det < 7279 && hit.pid != 22 && hit.pid != 0){
      q_shmx = q_shmx + 1;
    }
    if (trig29 && hit.det > 7001 && hit.det < 7279 && hit.pid != 22 && hit.pid != 0 ){
      q_shmx_trig = q_shmx_trig + 1;
    }
    // if (trig29 && hit.det ==331 && (hit.pid ==11 || hit.pid == -11) && hit.r > 1035 && hit.r<1190){
    //    h331_e_pz->Fill(hit.pz);
    // }
  }  // end 2nd loop over hits
  if (q_shmx != 0 && pmt_opt != 0){
  hshmx_pion-> Fill(q_shmx,pmt_opt);
  hshmx_pion_rate-> Fill(q_shmx,pmt_opt,rate);
  }
  if (q_shmx_trig != 0 && pmt_opt_trig != 0){
  hshmx_pion_trig->Fill(q_shmx_trig,pmt_opt_trig);
  hshmx_pion_trig_rate->Fill(q_shmx_trig,pmt_opt_trig,rate);
  }

 }    //end loop over events
//----------------------------------------------------------------------
  // Draw and save the histograms into png files in images subdirectory

 TStyle *st1 = new TStyle("st1","my style");
 st1->SetOptStat(1111111);
 st1->SetStatY(1);
 st1->SetStatX(0.9);
 gROOT->SetStyle("st1");

/*TCanvas *d1 = new TCanvas("d1", "SHMX hits vs PMT opt hits");
// gStyle->SetOptStat("ouneMRi");
d1->Divide(1,2);
d1->UseCurrentStyle();
d1->SetHighLightColor(0);
d1->cd(1);
hshmx_pion->Draw("colz");
d1->cd(2);
hshmx_pion_trig->Draw("colz");*/

TCanvas *d2 = new TCanvas("d2", "SHMX hits vs PMT opt hits (Rate)");
d2->Divide(1,2);
d2->UseCurrentStyle();
d2->SetHighLightColor(0);
d2->cd(1);
hshmx_pion_rate->Draw("colz");
Double_t rate_coinc_region = hshmx_pion_rate->Integral(0,15,0,25);
cout << "Rate of Coinc ="<< rate_coinc_region << "Hz/mA" <<endl;
d2->cd(2);
hshmx_pion_trig_rate->Draw("colz");
Double_t rate_coinc_region_trig = hshmx_pion_trig_rate->Integral(0,15,0,25);
cout << "Rate of Coinc ="<< rate_coinc_region_trig << "Hz/mA" <<endl;

cout << "Ratio of Trigger to Open = " << rate_coinc_region_trig / rate_coinc_region << "Hz/mA" <<endl;
/*
TCanvas *d3 = new TCanvas("d3", "Polar Plot: 8001 and 29");
d3->Divide(2,1);
d2->UseCurrentStyle();
d2->SetHighLightColor(0);
d3->cd(1);
h29_xy_e->Draw("colz");
d3->cd(2);
h29_xy_p->Draw("colz");
d3->cd();*/
}
