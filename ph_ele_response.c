//Script written on 9/19/2022 by Tristan Hurst
//
//Goal: determine the # of photo electrons per particle
//hitting 8000
//Det 8000: Pion det PMT
//Det 8001: Pion det body(lucite)
//
//

#include<TF1.h>

void ph_ele_response(const TString& files, int energycut){
    TChain* T = new TChain("T");
    if (files == "pion"){
    T->Add("~/SHMX_Layers/shmx_pion.root");}
    if (files == "moll"){
    T->Add("~/SHMX_Layers/shmx_moll.root");}

Double_t rate = 0;
std::vector<remollEventParticle_t>* parts = 0;
std::vector<remollGenericDetectorHit_t>* hits = 0;

T->SetBranchAddress("rate", &rate);
T->SetBranchAddress("hit", &hits);
T->SetBranchAddress("part", &parts);


const double pi = 3.14159265358979323846;
TH1F* h_8001_phi_pion = new TH1F("h_8001_phi_pion","#pi's, #mu 's, and, e 's Full Phi Extent Raw",200, -pi, pi);
h_8001_phi_pion->GetXaxis()->SetTitle("Phi (radians)"); 

TH1F* h_8001_phi_pion_pes = new TH1F("h_8001_phi_pion_pes","#pi's, #mu 's, and, e 's Full Phi Extent rated by PEs",200, -pi, pi);
h_8001_phi_pion_pes->GetXaxis()->SetTitle("Phi (radians)"); 

TH1F* h_8001_phi_pion_stack = new TH1F("h_8001_phi_pion_stack","#pi's, #mu 's, and, e 's Stacked Phi Extent Raw",15, 0, pi/14);
h_8001_phi_pion_stack->GetXaxis()->SetTitle("Phi (mod pi/14) (radians)"); 

TH1F* h_8001_phi_pion_pes_stack = new TH1F("h_8001_phi_pion_pes_stack","#pi's, #mu 's, and, e 's Stacked Phi Extent rated by PEs",15, 0, pi/14);
h_8001_phi_pion_pes_stack->GetXaxis()->SetTitle("Phi (mod pi/14) (radians)"); 



TH1F* h_8001_r_pion = new TH1F("h_8001_r_pion","Pion radial extent in Lucite",100,1100,1200);
h_8001_r_pion->GetXaxis()->SetTitle("Radius (mm)");

TH1F* h_8001_r_pion_pes = new TH1F("h_8001_r_pion_pes","Pion radial extent in Lucite weight by PEs",100,1100,1200);
h_8001_r_pion->GetXaxis()->SetTitle("Radius (mm)");



TAxis* a = h_8001_phi_pion->GetXaxis();
a->SetNdivisions(-502);
a->ChangeLabel(1,-1,-1,-1,-1,-1,"-#pi");
a->ChangeLabel(-1,-1,-1,-1,-1,-1,"#pi");
TAxis* b = h_8001_phi_pion_pes->GetXaxis();
b->SetNdivisions(-502);
b->ChangeLabel(1,-1,-1,-1,-1,-1,"-#pi");
b->ChangeLabel(-1,-1,-1,-1,-1,-1,"#pi");
TAxis* c = h_8001_phi_pion_stack->GetXaxis();
c->SetNdivisions(-502);
c->ChangeLabel(1,-1,-1,-1,-1,-1,"0");
c->ChangeLabel(-1,-1,-1,-1,-1,-1,"#pi/14");
TAxis* d = h_8001_phi_pion_pes_stack->GetXaxis();
d->SetNdivisions(-502);
d->ChangeLabel(1,-1,-1,-1,-1,-1,"0");
d->ChangeLabel(-1,-1,-1,-1,-1,-1,"#pi/14");


TH2F* h_lucite_pes = new TH2F("h_lucite_pes","PEs per charged part hit",10,0,10,50,0,1000);
h_lucite_pes->GetXaxis()->SetTitle("Charged Hits on Lucite"); h_lucite_pes->GetYaxis()->SetTitle("Det 8001 Photon Hits");

bool opthit = false; // flag for if optical photon hit the PMT of Pion, 8000
bool trig29 = false; //flag for if particles exit the pion det and hit the slice of 29

double mincut = -pi;
double maxcut = pi;


int total_events = T->GetEntries();

int counter = 1;

int lucitehits = 0; //# of chraged particle hits on 8001
int pmt_pes = 0; //# of photoelectrons on pmt cathode

const double minrad = 1110;//min radius of Dstream trigger cut
const double maxrad = 1150;//max radius of Dstream trigger cut

//T->GetEntries()
//20000
for (size_t iev = 0; iev < T->GetEntries();iev++){
    T->GetEntry(iev);
    lucitehits = 0;
    pmt_pes = 0;
    counter = 0;
    trig29 = false;
    opthit = false;

    for (size_t ihit = 0 ; ihit < hits->size(); ihit++){
        remollGenericDetectorHit_t hit = hits->at(ihit);
        
        if (hit.det == 8000 && hit.pid == 0){
            opthit = true;
            pmt_pes++;
        }
        if (hit.det == 29 && (hit.r < maxrad && hit.r>minrad) && (hit.pid == 11 ||hit.pid ==-11 ||hit.pid ==211 ||hit.pid ==-211 ||hit.pid == 13 ||hit.pid == -13 ||hit.pid == 2212)){
            trig29 = true;
        }

    }
    
    for (size_t ihit = 0 ; ihit < hits->size(); ihit++){
        remollGenericDetectorHit_t hit = hits->at(ihit);

//||hit.pid ==11 ||hit.pid==-11
        if (opthit && hit.p/MeV > energycut && hit.det == 8001 && (hit.pid ==211 || hit.pid ==-211||hit.pid == 13 ||hit.pid == -13||hit.pid ==11 ||hit.pid==-11) ){
        //-0.11211
            if (counter == 0){
            h_8001_phi_pion_stack->Fill(abs(fmod((hit.ph-0.11211),(pi/14))),1);
            h_8001_phi_pion_pes_stack->Fill(abs(fmod((hit.ph-0.11211),(pi/14))),pmt_pes);
            
            h_8001_phi_pion->Fill(hit.ph,1);
            h_8001_phi_pion_pes->Fill(hit.ph,pmt_pes);
            }

            h_8001_r_pion->Fill(hit.r);
            h_8001_r_pion_pes->Fill(hit.r,pmt_pes);
            //h_8001_phi_pion->Fill(hit.ph);
            counter++;

        }
        /* 
        if (opthit && hit.det == 8000 && hit.pid == 0){
            pmt_pes++;
        }*/
        if (opthit && hit.det ==8001 && (/*hit.pid == 11 ||hit.pid ==-11 ||*/hit.pid ==211 ||hit.pid ==-211 ||hit.pid == 13 ||hit.pid == -13 ||hit.pid == 2212)){
            lucitehits++;
        }
    }
    if(lucitehits != 0 && pmt_pes != 0){
    h_lucite_pes->Fill(lucitehits,pmt_pes);}
    
    if (iev%(T->GetEntries()/10) ==0){cout<<"Events Finished : "<<round((double)iev*100/(double)total_events)<<"\%"<<endl;}
}

//cout << abs(fmod(-pi,(double)3))<<endl;
/*
TCanvas* c1 = new TCanvas("c1", "Pion Lucite PMT PEs");
gStyle->SetOptStat("ouneMRi");
gStyle->SetStatX(0.9);
c1->Divide(1,1);
c1->cd(1);
h_lucite_pes->Draw("colz");
c1->cd();
*/
char energycut_str[20];
sprintf(energycut_str, "%d", energycut);
TString s (energycut_str);
TString gen (files);
TCanvas* c2 = new TCanvas("c2", "Pion hits on PD Phi",1080,720);
gStyle->SetOptStat("ouneMRi");
gStyle->SetStatX(1);
c2->SetTitle("Pion Generator hits >"+ s +"MeV on Lucite that produced light in Pion PMT");
c2->Divide(2,2);
c2->cd(1);
h_8001_phi_pion->Draw();
c2->cd(2);
h_8001_phi_pion_pes->Draw();
c2->cd(3);
h_8001_phi_pion_stack->Draw();
c2->cd(4);
h_8001_phi_pion_pes_stack->Draw();
c2->cd();
//c2->SaveAs("/home/mthurst/PionDetResponse/PionDet_Images/"+gen+"Generator hits >"+ s +"MeV on Lucite that produced light in Pion PMT.pdf");

TCanvas* p1 = new TCanvas("p1", "Pion Detector PE response",1080,720);
gStyle->SetOptStat("ouneMRi");
gStyle->SetStatX(1);
p1->SetTitle(gen+"Generator hits >"+s+"MeV which produced light in Pion PMT");
p1->Divide(2,1);
p1->cd(1);
h_8001_phi_pion_stack->Draw();
p1->cd(2);
h_8001_phi_pion_pes_stack->Draw();
p1->cd();


/*
TCanvas* c3 = new TCanvas("c3","XY Plot 8001");
//gStyle->SetOptStat("");
c3->Divide(1,1);
c3->cd(1);
h_x_y_8001->Draw("colz");*/
}