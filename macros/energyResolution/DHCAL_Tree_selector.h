//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jan 19 11:52:16 2011 by ROOT version 5.26/00
// from TTree DHCAL/DHCAL
// found on file: hitsDHCAL_TB.root
//////////////////////////////////////////////////////////





#ifndef DHCAL_Tree_selector_h
#define DHCAL_Tree_selector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TFrame.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TNtuple.h>
#include <iostream>
#include <fstream>
#include <Riostream.h>
#include <string>
#include <TStyle.h>
#include <TLegend.h>
#include <TText.h>
#include <TTree.h>
#include <TBranch.h>


using namespace std;
const Int_t kMaxDHCALEvent = 2000000;
TFile *f=new TFile("output.root","RECREATE");
TTree *sdhcal = new TTree("sdhcal", "all info");
TCanvas *c1 = new TCanvas("c1","",300,200,900,900);
TCanvas *c2 = new TCanvas("c2","",200,200,900,900);
Float_t pi1,pi2,pi3,RMAX,XGM,YGM;
Int_t ncounthr1,ncounthr2,ncounthr3,ncount,Nedge,zcount,NKFIN,KMAX,NKMAX,Nhole,Ndeb,Nfin,Ncern,Ndoub,Nlarg,nhough1,nhough2,nhough3,Nend, ncl,nis1,nis2,nis3, Ndir, nNdir;
Float_t thr1m,thr2m,thr3m,XGMu,YGMu;
Int_t numu,ncountpm,zcountm;
Int_t dt,Neutcosm,ntrk;
Int_t flag,kp;
Float_t dx,dy,xt,yt,zt,xp,yp;
Int_t ncountdh1, ncountdh2, ncountdh3;
Int_t ncountdl1, ncountdl2, ncountdl3;
unsigned long long SpillEvtTime;
//std::vector<double> densities;
Int_t densities[kMaxDHCALEvent];
Int_t densities_1[kMaxDHCALEvent];

TBranch *branche_1=sdhcal->Branch("ncounthr1", &ncounthr1);
TBranch *branche_2=sdhcal->Branch("ncounthr2", &ncounthr2);
TBranch *branche_3=sdhcal->Branch("ncounthr3", &ncounthr3);
TBranch *branche_4=sdhcal->Branch("pi1", &pi1);
TBranch *branche_5=sdhcal->Branch("pi2", &pi2);
TBranch *branche_6=sdhcal->Branch("pi3", &pi3);
TBranch *branche_7=sdhcal->Branch("ncount", &ncount);
TBranch *branche_8=sdhcal->Branch("Nedge", &Nedge);
TBranch *branche_9=sdhcal->Branch("zcount", &zcount);
TBranch *branche_10=sdhcal->Branch("NKFIN", &NKFIN);
TBranch *branche_11=sdhcal->Branch("KMAX", &KMAX);
TBranch *branche_12=sdhcal->Branch("NKMAX", &NKMAX);
TBranch *branche_13=sdhcal->Branch("RMAX", &RMAX);
TBranch *branche_14=sdhcal->Branch("Nhole", &Nhole);
TBranch *branche_15=sdhcal->Branch("Ndeb", &Ndeb);
TBranch *branche_16=sdhcal->Branch("Nfin", &Nfin);
TBranch *branche_17=sdhcal->Branch("XGM", &XGM);
TBranch *branche_18=sdhcal->Branch("YGM", &YGM);
TBranch *branche_19=sdhcal->Branch("Ncern", &Ncern);
TBranch *branche_20=sdhcal->Branch("Ndoub", &Ndoub);
TBranch *branche_21=sdhcal->Branch("Nlarg", &Nlarg);
TBranch *branche_22=sdhcal->Branch("numu", &numu);
TBranch *branche_23=sdhcal->Branch("thr1m", &thr1m);
TBranch *branche_24=sdhcal->Branch("thr2m", &thr2m);
TBranch *branche_25=sdhcal->Branch("thr3m", &thr3m);
TBranch *branche_27=sdhcal->Branch("XGMu", &XGMu);
TBranch *branche_28=sdhcal->Branch("YGMu", &YGMu);
TBranch *branche_29=sdhcal->Branch("dt", &dt);
TBranch *branche_30=sdhcal->Branch("dx", &dx);
TBranch *branche_31=sdhcal->Branch("dy", &dy);
TBranch *branche_32=sdhcal->Branch("Neutcosm", &Neutcosm);
TBranch *branche_33=sdhcal->Branch("nhough1", &nhough1);
TBranch *branche_34=sdhcal->Branch("nhough2", &nhough2);
TBranch *branche_35=sdhcal->Branch("nhough3", &nhough3);
TBranch *branche_36=sdhcal->Branch("Nend", &Nend);
TBranch *branche_37=sdhcal->Branch("ncl", &ncl);
TBranch *branche_38=sdhcal->Branch("nis1", &nis1);
TBranch *branche_39=sdhcal->Branch("nis2", &nis2);
TBranch *branche_40=sdhcal->Branch("nis3", &nis3);

TBranch *branche_41=sdhcal->Branch("Ndir", &Ndir);
TBranch *branche_42=sdhcal->Branch("nNdir", &nNdir);
TBranch *branche_43=sdhcal->Branch("TimeInSpill", &SpillEvtTime );

TBranch *branche_44=sdhcal->Branch("ncountdl1", &ncountdl1);
TBranch *branche_45=sdhcal->Branch("ncountdl2", &ncountdl2);
TBranch *branche_46=sdhcal->Branch("ncountdl3", &ncountdl3);
TBranch *branche_47=sdhcal->Branch("ncountdh1", &ncountdh1);
TBranch *branche_48=sdhcal->Branch("ncountdh2", &ncountdh2);
TBranch *branche_49=sdhcal->Branch("ncountdh3", &ncountdh3);

float ref1=79.82,ref2=9.205,ref3=0.8575;

TH2F H2("H2","",96,0,1,96,0,100);
TH1F HT1("HT1","",2001,0,2000);
TH1F HT2("HT2","",5001,0,500);
TH1F HT3("HT3","",5110,-0.1,5);
TH1F HT4("HT4","",2001,0,2000);
TH1F HT5("HT5","",51,0,50);
TH1F HT6("HT6","",101,0,100);

TH1F HTime[50];
TH1F HNhit[50];
TH1F HNoise[50];


//TNtuple ntup1("ntup1","","ncounthr1:ncounthr2:ncounthr3:pi1:pi2:pi3:ncount:Ned//ge:zcount:NKFIN:KMAX:NKMAX:RMAX:Nhole:Ndeb");
//:Nfin:XGM:YGM:cern:doub:Nlarg");
//TNtuple ntup2("ntup2","","nmu:thr1m:thr2m:thr3m:ncountm:zcountm");

//TNtuple ntup("ntup","","flag:ntrk[50]:xt:yt:zt:xp:yp:kp");

 
ofstream out;
int Ntrack=0;   

void CaliceStyle()
{
  /*CALICE style for figure: use in a ROOT macro like this:*/
  //gROOT->ProcessLine(".L ~/RootStuff/CaliceStyle.C");
  //CaliceStyle();

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTitleFont(42);
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetStatFont(42);

  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetCanvasColor(kWhite);  
  gStyle->SetOptStat(0); /*don't show statistics box*/
  gStyle->SetTitleSize(0.05, "xyz"); 
  gStyle->SetLegendBorderSize(1);

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);

  gROOT->ForceStyle();
}



class DHCAL_Tree_selector : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           DHCALEvent_;
   //UInt_t          DHCALEvent_fUniqueID[kMaxDHCALEvent];   //[DHCALEvent_]
   //UInt_t          DHCALEvent_fBits[kMaxDHCALEvent];   //[DHCALEvent_]
   Double_t        DHCALEvent_fX[kMaxDHCALEvent];   //[DHCALEvent_]
   Double_t        DHCALEvent_fY[kMaxDHCALEvent];   //[DHCALEvent_]
   Double_t        DHCALEvent_fZ[kMaxDHCALEvent];   //[DHCALEvent_]
   UInt_t          DHCALEvent_fDif_id[kMaxDHCALEvent];   //[DHCALEvent_]
   //Int_t           DHCALEvent_fAsic_id[kMaxDHCALEvent];   //[DHCALEvent_]
   //Int_t           DHCALEvent_fChan_id[kMaxDHCALEvent];   //[DHCALEvent_]
   Double_t        DHCALEvent_fTime[kMaxDHCALEvent];   //[DHCALEvent_]
   Int_t           DHCALEvent_fThr0[kMaxDHCALEvent];   //[DHCALEvent_]
   Int_t           DHCALEvent_fThr1[kMaxDHCALEvent];   //[DHCALEvent_]
   Int_t           DHCALEvent_fThr2[kMaxDHCALEvent];   //[DHCALEvent_]
   Int_t           Nframes;
   Int_t           NEvent;
   unsigned long long triggerTime;
   Int_t           eventTime;
// List of branches
   TBranch        *b_DHCALEvent_;   //!
   //TBranch        *b_DHCALEvent_fUniqueID;   //!
	 // TBranch        *b_DHCALEvent_fBits;   //!
   TBranch        *b_DHCALEvent_fX;   //!
   TBranch        *b_DHCALEvent_fY;   //!
   TBranch        *b_DHCALEvent_fZ;   //!
   TBranch        *b_DHCALEvent_fDif_id;   //!
   //TBranch        *b_DHCALEvent_fAsic_id;   //!
   //TBranch        *b_DHCALEvent_fChan_id;   //!
   TBranch        *b_DHCALEvent_fTime;   //!
   TBranch        *b_DHCALEvent_fThr0;   //!
   TBranch        *b_DHCALEvent_fThr1;   //!
   TBranch        *b_DHCALEvent_fThr2;   //!
   TBranch        *b_Nframes;   //!
   TBranch        *b_NEvent;   //!
   TBranch        *b_TriggerTime;   //!
   TBranch        *b_EventTime;   //!

   DHCAL_Tree_selector(TTree * /*tree*/ =0) { }
   virtual ~DHCAL_Tree_selector() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();


   //mes variables
   Int_t Ntoto;
   Double_t LessT;
   Double_t MaxT;
   unsigned long long _prevTime;
   unsigned long long timeRef;
   map<string, int> myCounters;
 
 
   ClassDef(DHCAL_Tree_selector,0);
};

#endif

#ifdef DHCAL_Tree_selector_cxx
void DHCAL_Tree_selector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("DHCALEvent", &DHCALEvent_, &b_DHCALEvent_);
	 // fChain->SetBranchAddress("DHCALEvent.fUniqueID", DHCALEvent_fUniqueID, &b_DHCALEvent_fUniqueID);
	 // fChain->SetBranchAddress("DHCALEvent.fBits", DHCALEvent_fBits, &b_DHCALEvent_fBits);
   fChain->SetBranchAddress("DHCALEvent.fX", DHCALEvent_fX, &b_DHCALEvent_fX);
   fChain->SetBranchAddress("DHCALEvent.fY", DHCALEvent_fY, &b_DHCALEvent_fY);
   fChain->SetBranchAddress("DHCALEvent.fZ", DHCALEvent_fZ, &b_DHCALEvent_fZ);
   fChain->SetBranchAddress("DHCALEvent.fDif_id", DHCALEvent_fDif_id, &b_DHCALEvent_fDif_id);
   //fChain->SetBranchAddress("DHCALEvent.fAsic_id", DHCALEvent_fAsic_id, &b_DHCALEvent_fAsic_id);
   //fChain->SetBranchAddress("DHCALEvent.fChan_id", DHCALEvent_fChan_id, &b_DHCALEvent_fChan_id);
   fChain->SetBranchAddress("DHCALEvent.fTime", DHCALEvent_fTime, &b_DHCALEvent_fTime);
   fChain->SetBranchAddress("DHCALEvent.fThr0", DHCALEvent_fThr0, &b_DHCALEvent_fThr0);
   fChain->SetBranchAddress("DHCALEvent.fThr1", DHCALEvent_fThr1, &b_DHCALEvent_fThr1);
   fChain->SetBranchAddress("DHCALEvent.fThr2", DHCALEvent_fThr2, &b_DHCALEvent_fThr2);
   fChain->SetBranchAddress("Nframes", &Nframes, &b_Nframes);
   fChain->SetBranchAddress("NEvent", &NEvent, &b_NEvent);
   fChain->SetBranchAddress("TriggerTime", &triggerTime, &b_TriggerTime);
   fChain->SetBranchAddress("eventTime", &eventTime, &b_EventTime);

}

Bool_t DHCAL_Tree_selector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


///

// Try to cut some noisy regions

Bool_t  CutNoise(Double_t x_x, Double_t y_y, Double_t z_z)
{
  
  if (x_x<100 && x_x>90 && y_y<40 && y_y>20 && z_z>120) 
    {
      return kTRUE;
    } else if (x_x<65 && x_x>55 && y_y<55 && y_y>45 && z_z<115 && z_z>105)
    {
      return kTRUE;
    } else if (x_x<65 && x_x>55 && y_y<65 && y_y>55 && z_z>125)
    {
      return kTRUE;
    } else if (x_x<55 && x_x>45 && y_y<100 && y_y>90 && z_z>120)
    {
      return kTRUE;
    } else if (y_y>97)
    {
      return kTRUE;
    } else if (x_x>97)
    {
      return kTRUE;
    } else if (x_x<23 && x_x>13 && y_y<10 && y_y>0 && z_z<125 && z_z>115)
    {
      return kTRUE;
    } else if (x_x<23 && x_x>13 && y_y<45 && y_y>35 && z_z<130 && z_z>120)
    {
      return kTRUE;
    } else if (x_x<45 && x_x>35 && y_y<25 && y_y>15 && z_z<130 && z_z>120)
    {
      return kTRUE;
    } else if (x_x<65 && x_x>55 && y_y<65 && y_y>55 && z_z>120)
    {
      return kTRUE;
    } else if (x_x<75 && x_x>65 && y_y<55 && y_y>45 && z_z<130 && z_z>120)
    {
      return kTRUE;
    } else if (x_x<75 && x_x>65 && y_y<45 && y_y>35 && z_z>125)
    {
      return kTRUE;
    } else if (x_x<85 && x_x>75 && y_y<65 && y_y>55 && z_z>125)
    {
      return kTRUE;
    } else if (x_x<15 && x_x>5 && y_y<5 && y_y>0 && z_z<130 && z_z>120)
    {
      return kTRUE;
    } else if (x_x<15 && x_x>5 && y_y<55 && y_y>45 && z_z<120 && z_z>110)
    {
      return kTRUE;
    } else if (x_x<45 && x_x>35 && y_y<95 && y_y>85 && z_z<125 && z_z>115)
    {
      return kTRUE;
    } else if (x_x<40 && x_x>30 && y_y<50 && y_y>40 && z_z>125)
    {
      return kTRUE;
    } else if (x_x<40 && x_x>30 && y_y<15 && y_y>5 && z_z>120)
    {
      return kTRUE;
    } else if (x_x<85 && x_x>75 && y_y<55 && y_y>45 && z_z<120 && z_z>110)
    {
      return kTRUE;
    } else if (x_x<95 && x_x>85 && y_y<95 && y_y>85 && z_z<10 && z_z>0)
    {
      return kTRUE;
    } else
    {
      return kFALSE;
    }
 }


///
#endif // #ifdef DHCAL_Tree_selector_cxx
