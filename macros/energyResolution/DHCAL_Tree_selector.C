#define DHCAL_Tree_selector_cxx
// The class definition in DHCAL_Tree_selector.h has been generated automatically
// by the ROOT utilityg TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For moredraw information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("DHCAL_Tree_selector.C")
// Root > T->Process("DHCAL_Tree_selector.C","some options")
// Root > T->Process("DHCAL_Tree_selector.C+")
//
#include "DHCAL_Tree_selector.h"
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TH3.h>
#include <TStyle.h>
#include <TProfile2D.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TFile.h>
#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>
#include <math.h>
#include <TPrincipal.h>
#include <TVector3.h>
#include <TMatrixD.h>

// DHCAL->Process("DHCAL_Tree_selector.C++")
std::string ou_txt_file = "/gridgroup/ilc/petr/data-Nov2012/corr-715751.txt";
std::string ou_txt_file1 = "/gridgroup/ilc/petr/data-Nov2012/calibfile_thr.txt";
std::string ou_txt_file2 = "/gridgroup/ilc/petr/data-Nov2012/corr-avet.txt";
std::string ou_txt_file3 = "/gridgroup/ilc/petr/data-Nov2012/corr-lin-2014.txt-726403";

double vect_co[50][12][12];
double vect1[13];
double vect_cot[47];
double vect_col[3];

void DHCAL_Tree_selector::Begin(TTree * /*tree*/)
{



  cout << "message avant de commencer" << endl;
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

       
  TString option = GetOption();
  std::ifstream file(ou_txt_file.c_str());
  if(file){
    while(1){
      if(!file.good()) break;
      int iplate,asic_i,asic_j;
      double ef,mu;
      double co;
      
      file >> iplate
           >> asic_i
           >> asic_j
           >> ef
           >> mu
           >> co;
      //std::cout << "co == " << co << "  //  plate == " << iplate << std::endl;
      //   vect_co[iplate][asic_i][asic_j] = co;
      vect_co[iplate][asic_i][asic_j] = co*mu;
    }
  }else
    std::cerr << "ton fichier est pourri mon pote !!!" << std::endl;
	

	
  std::ifstream file1(ou_txt_file1.c_str());
  if(file1){
    while(1){
      if(!file1.good()) break;
      double par1,par2,par3,par4,par5,par22,par23,par24,par25,par32,par33,par34,par35;
      
      file1 >> par1
	    >> par2
	    >> par3
	    >> par4
	    >> par5
	    >> par22
	    >> par23
	    >> par24
	    >> par25
            >> par32
	    >> par33
	    >> par34
	    >> par35;

      cout << "Calibration block thr1:" << par1 << " " << par2 << " " << par3  << " " << par4 << " " << par5 << endl;
      cout << "Calibration block thr2:" << par1 << " " << par22 << " " << par23  << " " << par24 << " " << par25 << endl;
      cout << "Calibration block thr3:" << par1 << " " << par32 << " " << par33  << " " << par34 << " " << par35 << endl;
      vect1[1]=par1;vect1[2]=par2;vect1[3]=par3;vect1[4]=par4;vect1[5]=par5;vect1[6]=par22;vect1[7]=par23;vect1[8]=par24;vect1[9]=par25;vect1[10]=par32;vect1[11]=par33;vect1[12]=par34;vect1[13]=par35;

    }
  }else
    std::cerr << " No calib. open !!! " << std::endl;

  // Average per layer correction from AP:
  std::ifstream file2(ou_txt_file2.c_str());
  if(file2){
    while(1){
      if(!file2.good()) break;
      int iplate2;
      double co2;
      
      file2 >> iplate2
	    >> co2;
      //     std::cout << "??????????????????? corr == " << co2 << "  //  plate == " << iplate2 << std::endl;
      vect_cot[iplate2] = co2;

    }
  }else
    std::cerr << "No file found !!!" << std::endl;  
 
  // Linear fit correction:
  std::ifstream file3(ou_txt_file3.c_str());
  if(file3){
    while(1){
      if(!file3.good()) break;
      double par1,par2,par3,par4,par5,par6;
      
      file3 >> par1
	    >> par2
	    >> par3;

      cout << "Linear fit calibration:" << par1 << " " << par2 << " " << par3  <<endl;
      vect_col[1]=par1;vect_col[2]=par2;vect_col[3]=par3;

    }
  }else
    std::cerr << " No calib. open !!! " << std::endl;
       


  myCounters["thr1"]=0;
  myCounters["thr2"]=0;
  myCounters["thr3"]=0;
  myCounters["anaevent"]=0;
  myCounters["mevent"]=0;
  myCounters["time"]=0;
  myCounters["distance"]=0;
  
  for(Int_t kpll=0; kpll< 50;kpll++)
    {
      /*	  
	TH1F 	*HTime("%d",kpll) = new TH1F(Form("HTime%d",kpll),"",1001,0,100); 
	TH1F	*HNhit("%d",kpll)= new TH1F(Form("HNhit%d",kpll),"",1001,0,1000); 
	TH1F 	*HNoise("%d",kpll)= new TH1F(Form("HNoise%d",kpll),"",101,0,10); 
      */	
		
      HTime[kpll].SetBins(1001,0,400);
      HNhit[kpll].SetBins(1001,0,4000);	 
      HNoise[kpll].SetBins(101,0,40);	 
		 
    }
  // yacine add
 
 
  //  for(int i=0;i < 50; i++)
  //  {
  //    h_density_asic[i] =  TProfile2D(Form("h_density_asic_%i",i),"",12,0,96,12,0,96);
  //  }




 
}

void DHCAL_Tree_selector::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).
  TString option = GetOption();
}

Bool_t DHCAL_Tree_selector::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either DHCAL_Tree_selector::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.
  //TFile f;

  
	// cout << "message avant de commencer1" << endl;  
  Float_t XGMBef=1000;
  Float_t YGMBef=1000;

  float XIN=20, YIN=20;  
  Int_t NZ = 55; 
  Int_t NT = 3500;
  Int_t NEPH=0;
	Int_t NHITPH=0;
  Int_t IBINBef=-100000;
  Int_t IBIN0=-10;
  Double_t XTX[NT],YTY[NT],ZTZ[NT],THR1[NT],THR2[NT],THR3[NT];
  Int_t TIME[NT],DIF[NT];
  Int_t TAG1[NT],TAG2[NT],TAGX[NT],TAGY[NT],TAGZ[NT], TTAG[NT];
  Int_t NATAGX[NT],NBTAGX[NT],NATAGY[NT],NBTAGY[NT],NATAGZ[NT],NBTAGZ[NT],TAGH[NT],TAGN[NT];
  Int_t Nho = 30;
   
  Int_t nclust[NT];
  Double_t xcl[NT];
  Double_t ycl[NT];
  Double_t zcl[NT];
  Int_t nhc[NT];
  Int_t sett[NT];
  Int_t tagcl[NT];
  ncl=0;

  Int_t NH=250;
  Int_t vec[NT][NH];


  Int_t irx,iry,irz,irmax=140,ith,ithmax=60,cellmax=14000; 
  Int_t IRX[Nho],IRY[irmax],ITH[ithmax];
  Int_t houghx[ithmax][irmax],houghy[ithmax][irmax], houghz[ithmax][irmax];
  Int_t houghxn[ithmax][irmax][Nho];
  Int_t houghyn[ithmax][irmax][Nho];
  Int_t houghzn[ithmax][irmax][Nho];
  // Double_t densities[kMaxDHCALEvent];

  float Pi=3.14,Pi2=1.57;
	//cout << "message avant de commencer2" << endl;

  GetEntry(entry);
	//cout << "message avant de commencer3" << endl;
  double E_corr = 0;
  double E_corr1 = 0;
  double E_corr2 = 0;
  double E_corr3 = 0;  
  bool draw=false;
  //	bool draw=true;
  bool detail=true;
	// bool detail=false;
  bool  waiting = false;
  bool effcorr = false;
  bool  Ef=true;
  bool elec=false;
  bool  cern=false;
  bool doub=false;
  bool deb=false;
  bool fin=false;
  bool neut=false;
  bool pion=false;
  bool muon=false;
  bool hough=true;
  Int_t npion=0;
  Int_t noise=0;
  bool calibr=false;
  bool calibl=true;
  bool clustering=true;

  Double_t XG[NZ],YG[NZ],ZG[NZ],MAXX[NZ],MAXY[NZ],MINX[NZ],MINY[NZ],EXX[NZ],EYY[NZ],EZ[NZ],ncountp[NZ],ncountpp[NZ],ncountp1[NZ],ncountp2[NZ],ncountp3[NZ],XT[NZ],YT[NZ],ZT[NZ],E_corrt[NZ];
  Double_t E_corre=0;

  bool PP[NZ];  
  bool PSELECT[NZ];  

  //  for(int i=0 ; i<kMaxDHCALEvent ; i++) densities[i] = -1;

  Double_t ZP[NZ];
  CaliceStyle();
  //TCanvas *c1 = new TCanvas("c1","",200,200,800,800);




  if(NEvent >0)
    {

      
			NEPH=0;
			NHITPH=0;
      //if(NEvent>100) continue;

			/*
        if( (triggerTime-_prevTime) > 100*TMath::Power(10.0,6) ){
      if( (triggerTime-_prevTime) > 100*pow(10.0,6.0) ){
	timeRef=triggerTime;
	  
	//	  cout << "SPIIIIIIIIIIIIIIIL: " << timeRef << endl;
      }
      SpillEvtTime=triggerTime-timeRef;
      _prevTime=triggerTime;

      //	c1->Clear();
      //c1->Divide(1,1);

      //        cout << "Je processe l'event number " << NEvent;
      //  cout << " qui a " << Nframes << " frames" << endl;
      //  cout << " qui a " << DHCALEvent_ << " framesMax" << endl;
	  

      //  	if (SpillEvtTime<20000000) { // AP: use only first 2 slots in time
			*/
      if( (triggerTime-_prevTime)*200 > 5*pow(10.0,9.0) ){
	timeRef=triggerTime;
	  
	//	  cout << "SPIIIIIIIIIIIIIIIL: " << timeRef << endl;
      }
      SpillEvtTime=triggerTime-timeRef;
      _prevTime=triggerTime;

      for(Int_t kpn=0; kpn< NZ;kpn++) 
	{
	  myCounters[Form("timemax%d",kpn)] = 0;
	  myCounters[Form("timemin%d",kpn)] = 10000000;
	  myCounters[Form("hits%d",kpn)] = 0;
	}
	  
      Int_t fframe =0;

      bool found=false;
      // 	const
 
 
      Int_t IBINP=-5;		 
 	
      /* 	
	 gROOT->SetStyle("Plain"); 
	 gROOT->ForceStyle(); 
	 gStyle->SetLabelFont(22,"X");
	 gStyle->SetLabelFont(22,"Y");
	 gStyle->SetTitleFont(22,"X");
	 gStyle->SetTitleFont(22,"Y");
	 gStyle->SetTitleFont(22," ");
	 gStyle->SetTextFont(22);
      */
      const Int_t NFixed= 10; 
      const Int_t NFT = 90 ;//150 minimum hit number;			
      const Int_t NFixedd= 3000;
      const Double_t dis = 3.; 
      const Int_t dlim=3;
      const Int_t ZFixed=7;//3;
      const Int_t ZFixedd= 40;
      //	cout << " Nframes1" << Nframes << endl;    
      //bool ok[Nframes]; 
      Int_t NTlim=3;//4
      Int_t NHMU=4;//4
      Int_t NHF=4;//4
      Int_t NHFF=2;//4
      Long64_t  maxtime=00;
      Long64_t  mintime=10e10;

      Int_t houghlim=5;
      Int_t Nhitl=10;
      Int_t Nisot=10;
      Int_t Nisol1=0;
      Int_t Nisol2=30;
      for (Long64_t kframe=0; kframe<Nframes;kframe++)
	{
	  // AP: exclude bad events for runs 716250(e) and 716305(pi): 
		//        if (NEvent==11105 || NEvent==11978) continue;
	  // AP: exclude bad events for runs 715553:
	  //	  if (NEvent==4774) continue;

	  if(int(DHCALEvent_fTime[kframe])>=maxtime)  maxtime=(DHCALEvent_fTime[kframe]);	
		 
		
	  if(int(DHCALEvent_fTime[kframe])<mintime)  mintime=(DHCALEvent_fTime[kframe]);	 

		//	cout << mintime << "   " << maxtime<<endl;
	  for(Int_t kpl=0; kpl< NZ;kpl++)
	    {
	      if( fabs(DHCALEvent_fZ[kframe]-2.8*kpl)<.1 &&DHCALEvent_fTime[kframe] >= myCounters[Form("timemax%d",kpl)]) myCounters[Form("timemax%d",kpl)] = DHCALEvent_fTime[kframe];
	      if( fabs(DHCALEvent_fZ[kframe]-2.8*kpl)<.1 &&DHCALEvent_fTime[kframe] <= myCounters[Form("timemin%d",kpl)]) myCounters[Form("timemin%d",kpl)] = DHCALEvent_fTime[kframe];
	      if( fabs(DHCALEvent_fZ[kframe]-2.8*kpl)<.1) myCounters[Form("hits%d",kpl)]++;
			 
			 
			 
	    }
		 
		 
		 
	}

      // if(maxtime>100000) 
      // {
      //cout << "maxtime"<< " " << maxtime*2./10000.<<endl;
      // maxtime=100000;
      // }

  
      {   
	TH1F  H("H","",maxtime+1,0,maxtime);
	for (Long64_t iframe=0; iframe< Nframes; iframe++)
    
	  {
	    H.Fill(DHCALEvent_fTime[iframe]);
		
	  }
	//				 H->SetMarkerColor(4);
	//			 H->SetMarkerStyle(20);
	//	c1->cd();
	//	 H->Draw();
	//getchar();
	//c1->Update();
	//c1->GetFrame()->SetFillColor(21);
	//c1->GetFrame()->SetBorderSize(12);
	//c1->Modified();

	Int_t NB0=0, NB1p=0, NB1m=0;

	IBINBef=IBIN0;
	for (Int_t  ibin=0; ibin< (maxtime+1); ibin++)
	  {
	    NB0= H.GetBinContent(ibin);
	    if(ibin<maxtime) Int_t NB1p= H.GetBinContent(ibin+1);
	    if(ibin>0)       Int_t NB1m= H.GetBinContent(ibin-1);
	    if(NB0>NFixed&&NB0>NB1p&&NB0>=NB1m&&ibin>IBINP+5)
	      {


		//5 time slots to protect against 
		{
		  IBIN0 = ibin;

		  IBINP=IBIN0+NTlim;
		}


 

	
	  
		TH2F H10("H10","",280,-135,145,101,0,100);
		TH2F H11("H11","",280,-135,145,101,0,100);
		TH2F H12("H12","",280,-135,145,101,0,100);

						 
		TH2F H20("H20","",280,-135,145,101,0,100);
		TH2F H21("H21","",280,-135,145,101,0,100);
		TH2F H22("H22","",280,-135,145,101,0,100);


		TH2F H12H("H12H","",280,-135,145,101,0,100);
		TH2F H22H("H22H","",280,-135,145,101,0,100);





						 
		TH3F H3D0("H3D0","",280,-135,145,101,0,100,101,0,100);  
		TH3F H3D1("H3D1","",280,-135,145,101,0,100,101,0,100);          
		TH3F H3D2("H3D2","",280,-135,145,101,0,100,101,0,100);          
		TH3F H3DH("H3DH","",280,-135,145,101,0,100,101,0,100);          

		TH3F H3DI("H3DI","",280,-135,145,101,0,100,101,0,100);          
	
	  
		//		 cout <<  "time slot" <<ibin<< "               NEvent" << NEvent<< endl;
		Nedge =0;
		ncount =0;
		Int_t ncountEd =0;

		ncounthr1 =0;
		ncounthr2 =0;
		ncounthr3 =0;
		Nlarg=0;
		float ncounthrc1 =0;
		float ncounthrc2 =0;
		float ncounthrc3 =0;

		float ncounthrm1 =0;
		float ncounthrm2 =0;
		float ncounthrm3 =0;




		zcount=0;
		Double_t zt =0;
		if(IBIN0>0)
		  {



		    ncount =0;

		    ncounthr1 =0;
		    ncounthr2 =0;
		    ncounthr3 =0;

		    ncounthrc1 =0;
		    ncounthrc2 =0;
		    ncounthrc3 =0;

		    ncounthrm1 =0;
		    ncounthrm2 =0;
		    ncounthrm3 =0;

		    nhough1=0;
		    nhough2=0;
		    nhough3=0;

		    nis1=0;
		    nis2=0;
		    nis3=0;
 
		    zcount=0;
					 						 
		    for(Int_t j=0; j<NZ; j++) 
		      {  


						if( j < 5)
							{
        ZP[j]=-130+j*9;
							}

			//PS G.Chambres seules				
			//	ZP[j]=j*9.;
						if(j>=5)
							{
						ZP[j]=(j-5)*2.8;
			// ZP[j]=j*9.-8;

              }
			PP[j]=false;
			/* 
			   if(j==0) ZP[j]=0; 
			   if(j==1) ZP[j]=3.3; 
			   if(j==2) ZP[j]=6.6;
			   if(j==3) ZP[j]=10;
			   if(j==4) ZP[j]=13;
			*/
		 
			XG[j]=0; 
			YG[j]=0; 
			MINX[j]=100;
			MINY[j]=100;
			MAXX[j]=0;
			MAXY[j]=0;
			ncountp[j]=0;
			ncountp1[j]=0;
			ncountp2[j]=0;
			ncountp3[j]=0;
			ncountpp[j]=0;

			XT[j]=0;
			YT[j]=0;
			EXX[j]=1;
			EYY[j]=1;

			//		h_density_asic[j].Reset(); // AP: Used to reset a hist cont.: we do not want Gb sizes of added stuff!

		      }  

		    //*******************
		    //	    for(Int_t j=0;j<NZ;j++) {
		    //  h_density_asic[j].Reset(); // AP: Used to reset a hist cont.: we do not want Gb sizes of added stuff!
		    // }
		    //	********************//

		    Int_t jk=0;
	
		    for (Int_t mframe=0; mframe< Nframes; mframe++)
		      {			     
					
				
					 
				
					
			if(fabs(DHCALEvent_fTime[mframe]-IBIN0)<NTlim)
			  //if scintillator
			  //		IBIN0= 10;
			  //			if(fabs(DHCALEvent_fTime[mframe]-IBIN0)<3)
           
			  {
                 
			    bool topo;
			    if(ncount>NT) continue;
			    if (DHCALEvent_fDif_id[mframe]>=205&&DHCALEvent_fDif_id[mframe]<=216) continue; // new changes 27/02/14
			    //    if(DHCALEvent_fDif_id[mframe]==125) cern=true;
			    for (Int_t mmframe= IBIN0-15; mmframe <IBIN0+15; mmframe++)
			      {
				if(DHCALEvent_fDif_id[mmframe]==1) cern=true;
			      }
			    Double_t xx= DHCALEvent_fX[mframe];
			    Double_t yy= DHCALEvent_fY[mframe];
			    Double_t zz= DHCALEvent_fZ[mframe];
			    Int_t time=DHCALEvent_fTime[mframe];
			    Int_t difid=DHCALEvent_fDif_id[mframe];
			    topo=false;
			    //	if(fabs(xx-50)<10&&fabs(yy-50)<10)ncountEd++;
			    if(fabs(xx-50)<41&&fabs(yy-50)<41)topo=true;
			    Int_t Thr0=DHCALEvent_fThr0[mframe];
			    Int_t Thr1=DHCALEvent_fThr1[mframe];
			    Int_t Thr2=DHCALEvent_fThr2[mframe];
			    
			    //  ncount++;
			    if(!topo)Nedge++;
			    //					{
            
			    if( ! CutNoise(xx, yy, zz))
			      {   
				// yacine contribution
				//std::cout << "yacine est la .... " << std::endl;

				// Protection against z=0 noise (already done in CutNoise):
				/*
				bool append=false; // new changes 27/02/14
				for(int ix=0; ix<NT; ix++){
				  if(xx==XTX[ix]&&yy==YTY[ix]&&zz==ZTZ[ix]){
				    append=true;
				    break;
				  }
				}
				if(append==true) continue;
				*/
				ncount++;
			    
			       
				if(Thr0>0)ncounthr1++;
				if(Thr1>0)ncounthr2++;
				if(Thr2>0)ncounthr3++;


				XTX[jk] = xx;
				YTY[jk] = yy;
				ZTZ[jk] = zz;
				TIME[jk] = time;
				DIF[jk] = difid;

				THR1[jk]=0;
				THR2[jk]=0;
				THR3[jk]=0;
				TAGH[jk]=0;
				TAGN[jk]=0;

				if(Thr0)THR1[jk]=1;
				if(Thr1)THR2[jk]=1;
				if(Thr2)THR3[jk]=1;
				jk++;
				for(Int_t k=0;k<NZ;k++) 
				  { 
				    ZG[k]=ZP[k];
				    if(fabs(zz-ZP[k])< 1)
				      {
					//	ZG[k]=zz;
					PP[k]=true;

										
					XG[k]=xx+XG[k];
					if(MINX[k]>xx)MINX[k]=xx;
					if(MAXX[k]<xx)MAXX[k]=xx;
					YG[k]=yy+YG[k];
					if(MINY[k]>yy)MINY[k]=yy;
					if(MAXY[k]<yy)MAXY[k]=yy;
					ncountp[k]++;
					if(Thr0>0)ncountp1[k]++;
					if(Thr1>0)ncountp2[k]++;
					if(Thr2>0)ncountp3[k]++;

                      
				      }

										
				  }

										
										
				//		cout <<" x  y  z "<< ""<<xx<<"  "<<yy<<" "<<zz<<endl;



				if(draw)
				  {

				    /*									            H3D0->Fill(zz,xx,yy);
				      if(Thr1>0) H3D1->Fill(zz,xx,yy);
				      if(Thr2>0) H3D2->Fill(zz,xx,yy);

				    */



				    if(Thr1<1&&Thr2<1)H10.Fill(zz,xx);	
				    if(Thr2>0) H12.Fill(zz,xx);
				    if(Thr1>0&&Thr2<1) H11.Fill(zz,xx);
	  
								



									
				    if(Thr1<1&&Thr2<1) H20.Fill(zz,yy);								
				    if(Thr2>0) H22.Fill(zz,yy);
				    if(Thr1>0&&Thr2<1) H21.Fill(zz,yy);
	  
	  
	  
				    if(Thr1<1&&Thr2<1) H3D0.Fill(zz,xx,yy);								
				    if(Thr2>0) H3D2.Fill(zz,xx,yy);
				    if(Thr1>0&&Thr2<1) H3D1.Fill(zz,xx,yy);


				  }

			      }//alexey										
			    //					}//topo
								
								
								
								
								
								
								
								
								
									 
                    
								
			  }
						

		      }

		  }

		   


		//for Remi
		if(hough)
		  {

		    for(Int_t kt=0;kt<ithmax;kt++)
		      {

			for(Int_t kr=0;kr<irmax;kr++)
			  {


			    houghx[kt][kr]=0;
			    houghy[kt][kr]=0;
			    houghz[kt][kr]=0;
			    for(Int_t kn=0;kn<Nho;kn++)

			      {



				houghxn[kt][kr][kn]=0;
				houghyn[kt][kr][kn]=0;
				houghzn[kt][kr][kn]=0;
			      }



			  }
               

		      }
					
		  }

		// end for remi

		////clustering
			
		for(Int_t klj=0; klj<ncount;klj++)

		  {
		    tagcl[klj]=0;
		    sett[klj]=-1;
		    nhc[ncl]=0;	
		  }

		Double_t nrcl=0;
				ncl=-1;
		Double_t xclmin=1000000;
		Double_t yclmin=1000000;
		Double_t zclmin=1000000;
		Double_t xclmax=0;
		Double_t yclmax=0;
		Double_t zclmax=0;

		for(Int_t kkj=0; kkj<ncount;kkj++)
		  {

				//		    if (DHCALEvent_fDif_id[kkj]>=205&&DHCALEvent_fDif_id[kkj]<=216) continue;

  
		    if(tagcl[kkj]==0)
		      {
			//intialisation du cluster
			ncl++;

			xcl[ncl]=XTX[kkj];
			ycl[ncl]=YTY[kkj];
			zcl[ncl]=ZTZ[kkj];
			sett[kkj]=ncl;
			tagcl[kkj]=1;
			nhc[ncl]=1;
			//définis vec
			vec[ncl][1]=kkj;
  



			Int_t NHP=0;
      Int_t NHA=1;

				while(NHP!=NHA)
					{

					NHP=nhc[ncl];
			for(Int_t mkj=0; mkj<ncount;mkj++)

			  {

					//			    if (DHCALEvent_fDif_id[mkj]>=205&&DHCALEvent_fDif_id[mkj]<=216) continue;


			    for(Int_t ncc=1; ncc<=nhc[ncl];ncc++)
			      {

							//							cout<<"here I am baby" << ncc << "    " << nhc[ncl] <<endl;
							if(tagcl[mkj]==0&&fabs(XTX[mkj]-XTX[vec[ncl][ncc]])<2&&fabs(YTY[mkj]-YTY[vec[ncl][ncc]])<2&&fabs(ZTZ[mkj]-ZTZ[vec[ncl][ncc]])<.5&&sett[vec[ncl][ncc]]!=sett[mkj]&&nhc[ncl]<NH)
				  {

      
				    tagcl[mkj]=1;
 
				    nhc[ncl]++;
            NHA=nhc[ncl];

						//						cout<< "test me " <<  NHP << " against " << NHA<<endl;
				    sett[mkj]=sett[kkj];
				    xcl[ncl]=xcl[ncl]+XTX[mkj];
				    ycl[ncl]=ycl[ncl]+YTY[mkj];
				    zcl[ncl]=zcl[ncl]+ZTZ[mkj];
				    vec[ncl][nhc[ncl]]=mkj;
				    

				    //	 	cout<< "xx=" <<xcl[ncl]<<endl;
					 				    continue;
				  }
		   
						}
			  }//end mkj

					}

			xcl[ncl]=xcl[ncl]/nhc[ncl];
			ycl[ncl]=ycl[ncl]/nhc[ncl];
			zcl[ncl]=zcl[ncl]/nhc[ncl];

			// AP: protection against big cluster positions and NANs

			if (xcl[ncl]<-200 || xcl[ncl]>200) continue;
			if (ycl[ncl]<-200 || ycl[ncl]>200) continue;
			//			if (zcl[ncl]<-200 || zcl[ncl]>200) continue;

			if (xcl[ncl]!=xcl[ncl]) continue;
			if (ycl[ncl]!=ycl[ncl]) continue;
			if (zcl[ncl]!=zcl[ncl]) continue;

			/**
			 // AP: some more cluster vars. Not really working yet!

			 if (xcl[ncl]>xclmax) xclmax=xcl[ncl];
			 if (ycl[ncl]>yclmax) yclmax=ycl[ncl];
			 if (zcl[ncl]>zclmax) zclmax=zcl[ncl];
			 if (xcl[ncl]<xclmin) xclmin=xcl[ncl];
			 if (ycl[ncl]<yclmin) yclmin=ycl[ncl];
			 if (zcl[ncl]<zclmin) zclmin=zcl[ncl];
  

			 rcl = sqrt(xcl[ncl]*xcl[ncl] + ycl[ncl]*ycl[ncl]); // AP: cluster radius
			 nrcl += rcl;
			 rclav = nrcl/ncl; // AP: average clust. radius
			 zcla = zcl[ncl]; // AP: z of cluster
			 clw = sqrt((xclmax-xclmin)*(xclmax-xclmin)+(yclmax-yclmin)*(yclmax-yclmin)); // AP: clusters width
			 clwav = clw/ncl; // AP: relative clust. width
			 cll = fabs(zclmax - zclmin); // AP: clusters length
			 clv = 0.33*3.14*clw*clw*0.25*cll; // AP: clusters volume assuming a cone shower shape

			 //			cout << "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH=" << rcl << endl; 
			 **/

			//	cout<< "ncl=  "<< ncl<< " nhit= " <<  nhc[ncl]<< "  x= " << xcl[ncl] << "   y=" << ycl[ncl] << "  z= " << zcl[ncl]  <<  endl;
		      } //if kkkj
		    
		  }// endkkj

		// cout << "number of luster    "<< ncl<<endl;





		/////////
		//for remi
		if(hough&&ncount>40)
		  {

 

		    //creation des histo bidimentionnels (rx,theta) et (ry,theta) 


		    for(Int_t jlj=1; jlj<ncl;jlj++)

		      {

			Int_t numbt=0;
			Int_t numbl=0;
			Int_t nth1=0;
			Int_t nth2=0;
			Int_t nth3=0;

			TAG1[jlj]=0;
			TAG2[jlj]=0;
			TAGX[jlj]=0;
			TAGY[jlj]=0;
			TAGZ[jlj]=0;
			TTAG[jlj]=0;
			NATAGX[jlj]=0;
			NATAGY[jlj]=0;
			NBTAGX[jlj]=0;
			NBTAGY[jlj]=0;
			NATAGZ[jlj]=0;
			NBTAGZ[jlj]=0;
          
			for(Int_t jjj=1; jjj<ncl+1;jjj++)

              
			  {

              
			    if(fabs(xcl[jjj]-xcl[jlj])>1&&fabs(xcl[jjj]-xcl[jlj])<9&&fabs(ycl[jjj]-ycl[jlj])>1&&fabs(ycl[jjj]-ycl[jlj])<9&&fabs(zcl[jjj]-zcl[jlj])<1)
			      {
				numbt++;
                  
			      }

			    if(fabs(zcl[jjj]-zcl[jlj])<9&&fabs(zcl[jjj]-zcl[jlj])>=0&&fabs(xcl[jjj]-xcl[jlj])<9&&fabs(ycl[jjj]-ycl[jlj])<9&&jjj!=jlj)numbl++;
			    //  if(fabs(zcl[jjj]-zcl[jlj])==0&&fabs(xcl[jjj]-xcl[jlj])<10&&fabs(ycl[jjj]-ycl[jlj])<10)numbl++;
               
          
			  }
			if(numbt<Nisot&&numbl>Nisol1&&numbl<Nisol2&&nhc[jlj]<Nhitl)TAG1[jlj]=1;
					 
			Double_t xxx=xcl[jlj];
			Double_t yyy=ycl[jlj];
			Double_t zzz=zcl[jlj];
							



         



			if(TAG1[jlj]==1) 
						
			  //remplir 
			  {
								
			    for (Int_t ith=0; ith<ithmax;ith++)

			      {
				irx=abs((int)(1*(xxx*cos(-Pi2+ith*(Pi/ithmax))+zzz*sin(-Pi2+ith*(Pi/ithmax)))));
				//AP: protection
				if (irx<-200 || irx>139) continue;						
				
				houghx[ith][irx]= 	houghx[ith][irx] +1;
				houghxn[ith][irx][houghx[ith][irx]]=jlj;
										

				//											houghz[ith][irz]= houghz[ith][irz]+1;  ;
				//	houghzn[ith][irz][houghz[ith][irz]]=jlj;  

				//	cout<< "ir=  "<< ith  << "   " << irx << "   " << iry << "  "<<	houghx[ith][irx]<< "     "<<	houghy[ith][iry]<<endl;

			      }  



							
			  }




		      }



		    Ndir=0;
		    nNdir=0;
		    Int_t TAGTX=0;
		    Int_t TAGTY=0;
		    Int_t TAGTZ=0;

		    for (Int_t ai=0;ai<ithmax;ai++)
		      {

			for (Int_t bi=0;bi<irmax;bi++)

			  {
          
			    bool comp=true;
			    if(houghx[ai][bi]>houghlim)
			      {
				
				Int_t maxl=houghx[ai][bi];
					
				for (Int_t lx=-0;lx<=0;lx++)
				  {
				    for (Int_t ly=-0;ly<=0;ly++)

				      {
					if((ai+lx)>=0&&(ai+lx)<=ithmax&&(bi+ly)>=0&&(bi+ly)<=irmax)
 
					  {
					    if(maxl<houghx[ai+lx][bi+ly]&&(lx!=0&&ly!=0)) comp=false;
					  }

						 
				      }



				  }

			      }

					
			    if(houghx[ai][bi]>houghlim&&comp)
			      {

				TAGTX=1;
				//nNdir=houghx[ai][bi];

				for(Int_t ktt=0;ktt<ithmax;ktt++)
				  {


				    //intialisation


				    for(Int_t krr=0;krr<irmax;krr++)
				      {


					houghy[ktt][krr]=0;
					houghz[ktt][krr]=0;

					for(Int_t knn=0;knn<Nho;knn++)

					  {
										
					    houghyn[ktt][krr][knn]=0;
					    houghzn[ktt][krr][knn]=0;
										
					  }



				      }

				  }
				//intialisation end

               

							






				for (Int_t ci=1;ci<houghx[ai][bi]+1;ci++)
				  {
				    //	cout << " cellx" <<  houghx[ai][bi]<< "    "<<  houghxn[ai][bi][ci]<<endl;
				    TAGX[houghxn[ai][bi][ci]]=1;

                 

				    Int_t jhj= houghxn[ai][bi][ci];



				    Double_t  xxv=xcl[jhj];

				    Double_t  yyv=ycl[jhj];

				    Double_t  zzv=zcl[jhj];

				    Int_t nunu=0;

				    TAG2[jhj]=0;


				    for (Int_t ei=1;ei<houghx[ai][bi]+1;ei++)
				      {

					

					Int_t jhj1=  houghxn[ai][bi][ei];
 
 
					Double_t  xxv1=xcl[jhj1];

					Double_t  yyv1=ycl[jhj1];

					Double_t  zzv1=zcl[jhj1];





					if( fabs(zzv-zzv1)<12&&fabs(zzv-zzv1)>1&&fabs(yyv-yyv1)<9&&fabs(xxv-xxv1)<9/*&&jhj!=jhj1*/)nunu++;
				

				      }

				    if((nunu>4)||(nunu>2&&((zzv>44*2.8)||(zzv<5*2.8))  )) TAG2[jhj]=1;

				    if(TAG2[jhj]==1)
				      {



					for (Int_t ithh=0; ithh<ithmax;ithh++)

					  {

					    Double_t yyyy=ycl[houghxn[ai][bi][ci]];
					    Double_t zzzz=zcl[houghxn[ai][bi][ci]];
					    Double_t xxxx=xcl[houghxn[ai][bi][ci]];
										
					    iry=abs((int)(1*(yyyy*cos(-Pi2+ithh*Pi/ithmax)+zzzz*sin(-Pi2+ithh*Pi/ithmax))));
					    irz=abs((int)(1*(xxxx*cos(-Pi2+ithh*Pi/ithmax)+yyyy*sin(-Pi2+ithh*Pi/ithmax))));
					    //AP: protection
					    if (iry<-200 || iry>139) continue;						
					    //   if (irz<-200 || irz>200) continue;						
																		
					    //	    cout << "IRZZZZZZZZZZZZZZ = " << irz << endl;

					    houghy[ithh][iry]= (houghy[ithh][iry])+1;  ;
					    houghyn[ithh][iry][houghy[ithh][iry]]=houghxn[ai][bi][ci]; 



					    houghz[ithh][irz]= (houghz[ithh][irz])+1;  ;
					    houghzn[ithh][irz][houghz[ithh][irz]]=houghxn[ai][bi][ci]; 
					  }

				      }

				  }//ci


 
				TAGTY=0;
				TAGTZ=0;
 
				for (Int_t aii=0;aii<ithmax;aii++)
				  {

				    for (Int_t bii=0;bii<irmax;bii++)

				      {






					bool compp=true;
					if(houghy[aii][bii]>houghlim)
					  {
				
					    Int_t maxll=houghy[aii][bii];
					
					    for (Int_t llx=-0;llx<=0;llx++)
					      {
						for (Int_t lly=-0;lly<=0;lly++)

						  {
						    if((aii+llx)>=0&&(aii+llx)<=ithmax&&(bii+lly)>=0&&(bii+lly)<=irmax)
 
						      {
							if(maxll<houghy[aii+llx][bii+lly]&&(llx!=0&&lly!=0)) compp=false;
						      }

						 
						  }



					      }

					  }


					if(houghy[aii][bii]>houghlim&&compp)
					  {



					    Ndir++;

					    TAGTY=1;
					    //	    cout << "Ndir evol "<<Ndir<< " aii "<< aii << " bii  " << bii<< " nomb "<< houghy[aii][bii]<<endl;
					    // cout << "Ndir evol "<<Ndir<< " ai "<< ai << " bi  " << bi<< " nombx "<< houghx[ai][bi]<<endl;
					    //   nNdir=houghy[aii][bii];




					    for (Int_t dii=1;dii<(houghy[aii][bii])+1;dii++)
					      {
						TAGY[houghyn[aii][bii][dii]]=1;
						TTAG[houghyn[aii][bii][dii]]=Ndir;	
       

					      }

					  }


				      }

				  }



				for (Int_t fii=0;fii<ithmax;fii++)
				  {

				    for (Int_t gii=0;gii<irmax;gii++)

				      {

					if(houghz[fii][gii]>houghlim)
					  {
					    TAGTZ=1;




					    nNdir=houghz[fii][gii];



 
					    for (Int_t hii=1;hii<(houghz[fii][gii])+1;hii++)
					      {
						TAGZ[houghzn[fii][gii][hii]]=1;


					      }

					  }


				      }

				  }




				//									NATAGX[houghxn[ai][bi][ci]]=ai;
				//									NBTAGX[houghxn[ai][bi][ci]]=bi;

			      }


			    //			if(TAGTX==1&&TAGTY==1&&TAGTZ==1)
			    //	{

			    //							Ndir++;

			    //	cout <<"-------Ndir =  "<<Ndir<<"----nNdir="<<nNdir<<endl;    
			    //	}
				
					

				











			  }

		      }






		  }//hough


		////////////





				

		if(hough&&ncount>40)
		  {

		    for(Int_t jmj=1;jmj<=ncl;jmj++)
		      {


			Double_t  xxh=xcl[jmj];

			Double_t  yyh=ycl[jmj];

			Double_t  zzh=zcl[jmj];

					
					
					 


			Int_t nfinal;
			 

	
			if((TAGX[jmj]==1&&TAGY[jmj]==1&&TAGZ[jmj]==1))
			  {


			    nfinal=0;
			    for( Int_t lll=1;lll<=ncl;lll++)

			      {
				float ddd=fabs(xxh-xcl[lll])+fabs(yyh-ycl[lll]);

				if(ddd<9&&fabs(zcl[jmj]-zcl[lll])<9&&TAGX[lll]==1&&TAGY[lll]==1&&TAGZ[lll]==1&&TTAG[lll]==TTAG[jmj]/*&&lll!=jmj*/)
				  { nfinal++;
	
				    //		cout << " 1 " <<jmj << "  2  "<< lll<<endl;
				  }
			      }



								   


			    if
			      (TAGX[jmj]==1&&TAGY[jmj]==1&&TAGZ[jmj]==1&&(nfinal>1))

			      {

		              	Double_t xxf= xcl[jmj];
				Double_t yyf= ycl[jmj];
				Double_t zzf= zcl[jmj];

				//			cout << jmj << "  Ndir="<< TTAG[jmj]   <<" nfinal " << nfinal <<endl;
				nfinal=0;
			        H12H.Fill(zzf,xxf);
				H22H.Fill(zzf,yyf);
				H3DH.Fill(zzf,xxf,yyf);
				for(Int_t nth=1;nth<=nhc[jmj];nth++)

				  {

				    TAGH[vec[jmj][nth]]=1;
				    if(THR1[vec[jmj][nth]]==1) nhough1++;
				    if(THR2[vec[jmj][nth]]==1) nhough2++;
				    if(THR3[vec[jmj][nth]]==1) nhough3++;

							

				  }


			      }



			    //	for(Int_t hhh=1; hhh<=nhc[jmj];hhh++)
			    {
			      //	Double_t xxf= XTX[vec[jmj][hhh]];
			      //	Double_t yyf= YTY[vec[jmj][hhh]];
			      //	Double_t zzf= ZTZ[vec[jmj][hhh]];
			      //			cout << nhc[jmj]<< "   "  << xxf<< "   "<< yyf<<"   "<<zzf<<endl;

			    }
	
			  }

		      }




		  }





		// end for rémi
		  
		// Density part:
		// 8+9+9 cells max
		 ncountdh1 = 0;
		 ncountdh2 = 0;
		 ncountdh3 = 0;
		 ncountdl1 = 0;
		 ncountdl2 = 0;
		 ncountdl3 = 0;
		double xNeighbor = 1.5;
		double yNeighbor = 1.5;
		double zNeighbor = 3.1;
		
		bool tag = false;
		for( Int_t h1=0 ; h1<ncount; h1++ ) {
		  densities[h1]=0;
		  //		  if (NEvent==59) {
		  /*
		  if (NEvent==77) {
		    
		    if (h1==0) cout << "time 1 hit 1 = " << TIME[0] << endl;
		    if (h1==(ncount-1)) cout << "time last hit 1 = " << TIME[ncount-1] << endl;
		    
		  }
		  */
	
		  for( int h2=0 ; h2<ncount ; h2++ ) {
		    if(h1==h2) continue;
		    if (DHCALEvent_fDif_id[h1]>=205&&DHCALEvent_fDif_id[h1]<=216) continue;
		    if (DHCALEvent_fDif_id[h2]>=205&&DHCALEvent_fDif_id[h2]<=216) continue;

		    if( fabs(XTX[h1] - XTX[h2] ) < xNeighbor
			&& fabs(YTY[h1] - YTY[h2] ) < yNeighbor
			&& fabs(ZTZ[h1] - ZTZ[h2] ) < zNeighbor ) {
		      
		      /*
		      if (XTX[h1]==XTX[h2]&&YTY[h1]==YTY[h2]&&ZTZ[h1]==ZTZ[h2]){
			  // &&TIME[h1]==TIME[h2]) {
			tag = true;
			break;
			cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
			cout << "event = " << NEvent << endl;
			cout << "nhit = " << ncount << endl;
			cout << "hit1 = " << h1 << " hit2 = " << h2 << endl;
			cout << "x1 = " << XTX[h1] << " x2 = " << XTX[h2] <<  endl;
			cout << "y1 = " << YTY[h1] << " y2 = " << YTY[h2] <<  endl;
			cout << "z1 = " << ZTZ[h1] << " z2 = " << ZTZ[h2] <<  endl;	
			cout << "t1 = " << TIME[h1] << " t2 = " << TIME[h2] <<  endl;	
			cout << "dif1 = " << DIF[h1] << " dif2 = " << DIF[h2] <<  endl;	
			cout << "bcid = " << triggerTime <<  endl;
			cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

		      }
		      if (tag==false)
		      */
		      densities[h1]++;
		    
		    }
		    
		  }
		  //  if (tag==true) continue;
		  //cout << "densities=" << densities[h1] << " density=" << density << endl;
		  //	  h_density_asic[int(DHCALEvent_fZ[h1]/2.8)].Fill(int(DHCALEvent_fX[h1]/(1.04125)),int(DHCALEvent_fY[h1]/(1.04125)),densities[h1]);
		  
		  //		  cout << "density : " << densities[h1] << endl;
		  // std::cout << "density:" << densities[h1]/63. << "; number of hit is: " << ncount << std::endl;
		  if ((densities[h1]/27.)>0.3) {
		    if( DHCALEvent_fThr0[h1] ) ncountdh1++;
		    if( DHCALEvent_fThr1[h1] ) ncountdh2++;
		    if( DHCALEvent_fThr2[h1] ) ncountdh3++;
		  } 
		  else {
		    if( DHCALEvent_fThr0[h1] ) ncountdl1++;
		    if( DHCALEvent_fThr1[h1] ) ncountdl2++;
		    if( DHCALEvent_fThr2[h1] ) ncountdl3++;   
		  }
		}

		// isolation
	       
		
		for(Int_t mlm=0;mlm<ncount+1;mlm++)
		  {

		    Int_t Nis=0;
		    Int_t NAis=0;
		    if(TAGH[mlm]==0)

		      {

			for(Int_t lml=0;lml<ncount;lml++)

			  {

			    if(fabs(XTX[mlm]-XTX[lml])<10&&fabs(YTY[mlm]-YTY[lml])<10&&fabs(ZTZ[mlm]-ZTZ[lml])<10&&TAGH[lml]==0) Nis++;
			    if(fabs(XTX[mlm]-XTX[lml])<5&&fabs(YTY[mlm]-YTY[lml])<5&&fabs(ZTZ[mlm]-ZTZ[lml])<5&&TAGH[lml]==1&&lml!=mlm) NAis++;




			  }
            
     


			

			if(Nis<9&&NAis<1) 
			  {
			    TAGN[mlm]=1;
			    if(THR1[mlm]==1)nis1++;  
			    if(THR2[mlm]==1)nis2++;  
			    if(THR3[mlm]==1)nis3++;  

			    H3DI.Fill(ZTZ[mlm],XTX[mlm],YTY[mlm]);
			  }

		      }

			
			
		    //	    if(TAGN[mlm]==0) // AP 09/05/14, commented to use wider range
		      {

			for (Int_t ko=0;ko<NZ;ko++)
			  {
			    if(fabs(ZTZ[mlm]-ZP[ko])< 1)
			      {
											 										
				//	if(MINX[ko]>XTX[mlm])MINX[ko]=XTX[mlm];
				//	if(MAXX[ko]<XTX[mlm])MAXX[ko]=XTX[mlm];
											
				//		if(MINY[ko]>YTY[mlm])MINY[ko]=YTY[mlm];
				//	if(MAXY[ko]<YTY[mlm])MAXY[ko]=YTY[mlm];
											           
			      }

			  }
										
		      }		



		  }








		   
		zcount =0;	   
		noise = 0;
		neut=false;
		npion=0;
		Int_t  nelec=0;
		elec=false;
		Int_t First=0;
		Ndeb=0;
		Nend=0;
		Nfin=50;
		Nhole=0;
   




		XGM=0;
		YGM=0;
		Int_t KGM=0;

		XGMu=0;
		YGMu=0;
		Int_t KGMu=0;


		Int_t KGPi=0;
		Int_t  XGPi=0;
		Int_t YGPi=0;



		Float_t MAXXT=0;
		Float_t MAXYT=0;
		RMAX=0;
		KMAX=0;
		NKMAX=0;
		NKFIN=0;

		Neutcosm=0;
		if(ncountp[0]+ncountp[1]+ncountp[2]+ncountp[3]+ncountp[4]<5)Neutcosm=1; 

		
		// Dense part (densities>0.2)
		// Not dense: densities=0-0.2
	       	       		
		// Asic density part 
		
	       
		
		for(Int_t k=0;k<NZ;k++) 
		  {
		    if(ncountp[k]>0)
		      { 
			XG[k]=XG[k]/ncountp[k];
			YG[k]=YG[k]/ncountp[k];
			 


			
				 
			
			if(k<50&&k>=0&&ncountp[k]>0&&Neutcosm==0)
			  {
			    KGM++;
			    XGM=XG[k]+XGM;
			    YGM=YG[k]+YGM;
			  }

			if(KGM>0&&ncount/KGM>NHMU)
			  {
			    if(k<10&&k>1&&ncountp[k]>0)
			      {
				KGPi++;
				XGPi=XG[k]+XGPi;
				YGPi=YG[k]+YGPi;
			      }
			  }
        
			if(KGM>0&&ncount/KGM<NHMU)
			  {
			    if(k<10&&k>1&&ncountp[k]>0)
			      {
				KGMu++;
				XGMu=XG[k]+XGMu;
				YGMu=YG[k]+YGMu;
			      }
			  }


		      }





		    if(ncountp[k]>100) noise++;

		 
				
		    if(ncountp[k]>NHF&&ncountp[k+1]>NHF&&ncountp[k+2]>NHF&&ncountp[k+3]>NHF&&Ndeb==0) Ndeb=k;

  
		    if(k>=40) Nend=Nend+ncountp[k];

		    if(ncountp[k]>NHFF&&ncountp[k+1]<NHFF&&ncountp[k+2]<NHFF&&ncountp[k+3]<NHFF&&k<47&&Nfin==50) Nfin=k;
       
		    if(ncountp[k]>1) NKFIN=k;

      	

		    if(ncountp[k]<1&&k<46&&(ncountp[k+1]>4|| ncountp[k+2]>4||ncountp[k+3]>4))Nhole++;

		    if(PP[k])zcount++;

		    if((MAXX[k]-MINX[k])>5||(MAXY[k]-MINY[k])>5)Nlarg++;
		    if((MAXX[k]-MINX[k])>MAXXT) MAXXT=MAXX[k]-MINX[k];
		    if((MAXY[k]-MINY[k])>MAXYT) MAXYT=MAXY[k]-MINY[k];
		    if(ncountp[k]*ncountp[k]*(MAXY[k]-MINY[k])*(MAXY[k]-MINY[k])+(MAXY[k]-MINY[k])*(MAXY[k]-MINY[k])>RMAX*RMAX&&k>0)

		      {
			//					cout <<"semme"<<RMAX<<endl;
			RMAX =ncountp[k]*sqrt((MAXX[k]-MINX[k])*(MAXX[k]-MINX[k])+(MAXY[k]-MINY[k])*(MAXY[k]-MINY[k]));
			KMAX=k;
			NKMAX=ncountp[k]; 
		      }  


		  }

		RMAX=RMAX/NKMAX;
		if(ncountp[47]>NHF) Nfin=47;
		if(ncountp[48]>NHF) Nfin=48;
		if(ncountp[49]>NHF) Nfin=49;
		if(KGM>0)     XGM=XGM/KGM;
		if(KGM>0)     YGM=YGM/KGM;
 

		if(KGPi>0)     XGPi=XGPi/KGPi;
		if(KGPi>0)     YGPi=YGPi/KGPi;

		if(KGMu>0)     XGMu=XGMu/KGMu;
		if(KGMu>0)     YGMu=YGMu/KGMu;
 

		Int_t Nmu=0;
		
		for(Int_t mima=0; mima<10;mima++)  
		  {
		    
		    //		cout << "  here it"<< ((MAXX[mima]-MINX[mima])*(MAXX[mima]-MINX[mima]) +   (MAXY[mima]-MINY[mima])*(MAXY[mima]-MINY[mima]) ) <<endl;
		    
		    
		    if( (((MAXX[mima]-MINX[mima])*(MAXX[mima]-MINX[mima]) +   (MAXY[mima]-MINY[mima])*(MAXY[mima]-MINY[mima]) )>9)&&PP[mima]  )Nmu++; 
		  }
		
		if(Nmu>4) doub=true;
		
	// if((MAXX[5]-MINX[5])>5&&(MAXX[6]-MINX[6])>5&&(MAXX[7]-MINX[7])>5&&(MAXX[8]-MINX[8])>5&&(MAXX[9]-MINX[9])>5) doub=true;


       


       

		/*
		   
		if(nelec>9)elec=true;
		if(npion>25) pion=true;
		if(Ndeb<3) deb=true;
		if(Nfin>2) fin=true;
		*/

		   		  
		//selection 

		// muon 

		  
		   
		   
	
		if(Ef) 
		  //effiency study
		  {	   

		    for(Int_t kE=0;kE<NZ;kE++) 

		      {
			Int_t kj=0; 
			for(Int_t k=0;k<NZ;k++) 
			  { 
       
			
			    // if(ncountp[k]>0) XG[k]=XG[k]/ncountp[k];
			    // if(ncountp[k]>0) YG[k]=YG[k]/ncountp[k];


			    if(k!=kE&&fabs(k-kE)<20)				
			      {
				if(((MAXX[k]-MINX[k])*(MAXX[k]-MINX[k])+(MAXY[k]-MINY[k])*(MAXY[k]-MINY[k]))<dis*dis&&ncountp[k]<3&&PP[k])
				  //if(PP[k])
				  {
				    //							cout << "testssssssssss"<< k << " " <<((MAXX[k]-MINX[k])*(MAXX[k]-MINX[k])+(MAXY[k]-MINY[k])*(MAXY[k]-MINY[k])) <<endl;
				    PSELECT[kj]=true;
				    XT[kj]=XG[k] ;
				    YT[kj]=YG[k] ;
				    ZT[kj]=ZG[k] ;
				    EXX[kj]=fabs(MAXX[k]-MINX[k]);
				    if(EXX[kj]<1) EXX[kj]=1.;
				    EXX[kj]= EXX[kj]/sqrt(12);
				    EXX[kj]= EXX[kj]/sqrt(12);
				    EYY[kj]=fabs(MAXY[k]-MINY[k]);
				    if(EYY[kj]<1) EYY[kj]=1.;
				    EYY[kj]= EYY[kj]/sqrt(12);
				    EZ[kj]=1.;

				    kj++;
			
				  }
			      }
			  }

			// cout << "kj"<< kj <<endl;
			//attentionnnnnnnnnnnnnnnnnnnnnnn	 
			if(kj>=6) 


			  {
	
	


			    //	c2->Clear();					 	
  
			    //	c2->Divide(1,2);

			    //c2->cd(1);


			    TGraphErrors* grxz = new TGraphErrors(kj,ZT,XT,EZ,EXX);

			    grxz->SetLineColor(2);
			    grxz->SetLineWidth(4);
			    grxz->SetMarkerColor(4);
			    grxz->SetMarkerStyle(21);
			    grxz->GetXaxis()->SetTitle("Z title");
			    grxz->GetYaxis()->SetTitle("X title");
	 
	
			    grxz->Fit("pol1","Q");
			    TF1 *myfitxz = (TF1*) grxz->GetFunction("pol1");
			    Double_t pxz0 = myfitxz->GetParameter(0);           
			    Double_t  pxz1 = myfitxz->GetParameter(1);
			    Double_t  kxz = myfitxz->GetChisquare();



			    //				if(kxz<10) grxz->Draw("AP");

			    //c2->cd(2);

			    grxz->Delete();
	
		
			    TGraphErrors* gryz = new TGraphErrors(kj,ZT,YT,EZ,EYY);



			    gryz->SetLineColor(2);
			    gryz->SetLineWidth(4);
			    gryz->SetMarkerColor(4);
			    gryz->SetMarkerStyle(21);
			    gryz->GetXaxis()->SetTitle("Z title");
			    gryz->GetYaxis()->SetTitle("Y title");



			    gryz->Fit("pol1","Q");
			    TF1 *myfityz = (TF1*) gryz->GetFunction("pol1");
			    Double_t pyz0 = myfityz->GetParameter(0);           
			    Double_t  pyz1 = myfityz->GetParameter(1);
			    Double_t  kyz = myfityz->GetChisquare();

			    //		if(kyz<10) gryz->Draw("AP");
		


			    gryz->Delete();

			    //		 c2->Update();
			    //	c2->GetFrame()->SetFillColor(21);
			    //	c2->GetFrame()->SetBorderSize(12);
			    //	c2->Modified();
		

 
			    ////////  for(Int_t k=0;k<NZ;k++) {
	    
			    bool candid = false;
			    Double_t Xexp = pxz0+pxz1*ZG[kE];
			    Double_t Yexp = pyz0+pyz1*ZG[kE];

			    //							cout << "here0   "<< XGM<< "      "<<YGM << endl;



			    if ((fabs(Xexp-XGPi)<XIN&&fabs(Yexp-YGPi)<YIN)&&kyz<10&&kxz<10) 
								
			      //	cout << "here0   "<< ZG[kE]<< "   " <<  pyz0 <<" " <<pyz1 << " " << kyz << endl;
			      //        if (kyz<10&&kxz<10) 

			      {
				myCounters[Form("track%d",kE)]++;
				//**	cout  << "here1   "<< kE << "   " <<  myCounters[Form("track%d",kE)] << endl;
				candid=true;


	

				//					 	cout << "here   "<< ZG[kE]<< "   " <<  pxz0 <<" " <<pxz1 << " " << kxz << endl;

			      }
			    bool  dd=false;
			    ncountpp[kE]=0;


			
			    for(Int_t jkj=0; jkj<ncount;jkj++)

			      {

					

				//					cout << "here0   "<< XTX[jkj]<< " "<<  YTY[jkj] << " "<<  endl;
				if(PP[kE]&&fabs(ZG[kE]-ZTZ[jkj])<.5&&fabs(Xexp-XTX[jkj])<dlim&&fabs(Yexp-YTY[jkj])<dlim)
				  {
				    dd=true;
				    ncountpp[kE]++;
				    /*						
				      H2.Fill(XTX[jkj],YTY[jkj]);
				      xt=XTX[jkj];
				      yt=YTY[jkj];
				      zt=ZG[kE];
								
				      ntup.Fill(ntrk,xt,yt,zt);
				    */						 
							
								
								
				  }



			      }
			    //					cout << "here   "<< PP[kE]<< " "<<  candid << " "<< dd  << endl;
								    if(dd&&candid)
			      {

         
				myCounters[Form("trackf%d",kE)]++;
				myCounters[Form("multip%d",kE)]=myCounters[Form("multip%d",kE)]+ncountpp[kE];

						}
										//										cout << "here2 "  <<kE << " " <<myCounters[Form("trackf%d",kE)] << " "    << myCounters[Form("track%d",kE)] << " " << Double_t(myCounters[Form("trackf%d",kE)])/ myCounters[Form("track%d",kE)] << "  " << Double_t(myCounters[Form("multip%d",kE)])/myCounters[Form("track%d",kE)]  << endl;
			



			      //end of condition

			  }


		      }



		  }








		//end of efficiency study

 

          
        	///tracer l'histo ici

		if(draw)
		  {
		    /*
		      H3D0->SetMarkerColor(4);//bleu
		      H3D1->SetMarkerColor(3); //vert
		      H3D2->SetMarkerColor(2); //rouge

		      H3D0->SetMarkerStyle(20);
		      H3D1->SetMarkerStyle(20);
		      H3D2->SetMarkerStyle(20);

		      H3D2->SetMarkerSize(.5);
		      H3D1->SetMarkerSize(.5);
		      H3D0->SetMarkerSize(.5);

		    */
	  
	  
		    H3D0.SetMarkerColor(3);
		    H3D1.SetMarkerColor(3);//4
		    H3D2.SetMarkerColor(3);//4

		    H3DH.SetMarkerColor(1);

		    H3DI.SetMarkerColor(2);

	  
		    H3D0.SetMarkerStyle(20);
		    H3D1.SetMarkerStyle(20);
		    H3D2.SetMarkerStyle(20);

		    H3DH.SetMarkerStyle(20);

		    H3DI.SetMarkerStyle(20);
	  
		    H3D0.SetMarkerSize(0.40);
		    H3D1.SetMarkerSize(0.40);
		    H3D2.SetMarkerSize(0.40);

		    H3DH.SetMarkerSize(0.7);
	  
		    H3DI.SetMarkerSize(0.8);


					 
		    H10.SetMarkerColor(3);
		    H11.SetMarkerColor(4);
		    H12.SetMarkerColor(2);


		    H12H.SetMarkerColor(1);

		    H10.SetMarkerStyle(20);
		    H11.SetMarkerStyle(20);
		    H12.SetMarkerStyle(20);

		    H12H.SetMarkerStyle(20);


		    H10.SetMarkerSize(0.70);
		    H11.SetMarkerSize(0.70);
		    H12.SetMarkerSize(0.70);

		    H12H.SetMarkerSize(0.75);

		    H20.SetMarkerColor(3);
		    H21.SetMarkerColor(4);
		    H22.SetMarkerColor(2);


		    H22H.SetMarkerColor(1);

				



		    H20.SetMarkerStyle(20);
		    H21.SetMarkerStyle(20);
		    H22.SetMarkerStyle(20);

		    H22H.SetMarkerStyle(20);

		    H20.SetMarkerSize(0.70);
		    H21.SetMarkerSize(0.70);
		    H22.SetMarkerSize(0.70);

		    H22H.SetMarkerSize(0.75);


				 

		    //H1->SetMarkerStyle(20);
		  }
		//selection
	


		Float_t  mu1,mu2,mu3;//pi1,pi2,pi3;

	   
		// if(ncount>=NFixed&&zcount>=ZFixed&&npion>2&&!doub&&!neut&&noise==0)
		// if(ncount>=NFixed&&zcount>=ZFixed&&!noise)
		// if(ncount>=NFixed&&zcount>=ZFixed&&First>10&&First<100&&!doub)
		// if(ncount>=NFixed&&zcount>=ZFixed&&doub)
		//---------
		//electrons
		//  if(ncount>=NFT&&zcount>=ZFixed&&!doub&&ncount/zcount>2.5&&elec&&!doub&&!pion)
		//----------
		//muons


	

		if(zcount>ZFixed&&ncount/zcount<2.2&&ncount>NFixed&&ncount<200&&Neutcosm==0)

		  {


		    myCounters["mevent"]= myCounters["mevent"]+1;

 
		    for (Int_t kE=0; kE<NZ;kE++)
		      {
	
			Double_t effm =0;
			effm =Double_t(myCounters[Form("trackf%d",kE)])/Double_t(myCounters[Form("track%d",kE)]); 


			Double_t multipm=1;
			multipm =Double_t(myCounters[Form("multip%d",kE)])/Double_t(myCounters[Form("trackf%d",kE)]);

			if(effm>0)
			  {

			    ncounthrm1= (ncountp1[kE]/effm/multipm)+ncounthrm1;
			    ncounthrm2= (ncountp2[kE]/effm/multipm)+ncounthrm2;
			    ncounthrm3= (ncountp3[kE]/effm/multipm)+ncounthrm3;

			  }


			if(effm<=0)
			  {
			    ncounthrm1= (ncountp1[kE])+ncounthrm1;
			    ncounthrm2= (ncountp2[kE])+ncounthrm2;
			    ncounthrm3= (ncountp3[kE])+ncounthrm3;

			  }

		      }      
		    myCounters["mu1"]= myCounters["mu1"]+ncounthrm1;
		    myCounters["mu2"]= myCounters["mu2"]+ncounthrm2;
		    myCounters["mu3"]= myCounters["mu3"]+ncounthrm3;

	 
		    if(myCounters["mevent"]>0)
						 
		      {
			mu1=float(myCounters["mu1"])/ myCounters["mevent"];
			mu2=float(myCounters["mu2"])/ myCounters["mevent"];
			mu3=float(myCounters["mu3"])/ myCounters["mevent"];

			//			cout << "mevent000= " << myCounters["mevent"] <<  "   mu1= "<<mu1 << "  mu2= "<< mu2<< "  mu3= "<<mu3<<endl; 
			//		cout << "x= " << XGMu<< " y=" << YGMu<<endl; 

		      }
		    //nmu:thrm:thr1m:thr2m:ncountm:zcountm
	
		    Int_t nmu=myCounters["mevent"];
		    //				ntup2.Fill(nmu,ncounthrm1,ncounthrm2,ncounthrm3,ncount,zcount);
 

		  }

		//-----	
		//pions


		Ndoub=0;

		Ncern=0;
		if(cern)Ncern=1;
		if(doub)Ndoub=1;
		//	cout<<"count "<<ncount<< "zcount" << zcount<< endl;
		//		if(zcount>ZFixed&&ncount>NFixed&&ncount>150&&ncount/zcount>2&&noise<1)
		if(zcount>ZFixed&&ncount>NFixed&&ncount>40&&ncount/zcount>2&&noise<1)
		//// imads	if(zcount>ZFixed&&ncount>NFixed&&ncount>200&&noise<1&&Ndoub<1)
		// if ( (Ndeb>3||(zcount>30&&Ndeb>1))&&Nedge/ncount<.05&&(ncount/zcount>2&&RMAX/zcount>.05)&&Ndeb<15&&fabs(XGM-50)<35&&fabs(YGM-50)<35&&Nhole<1&&Ndoub<1&&noise<1)



		  {	

		    NEPH++;
				 NHITPH=NHITPH+ncount;
		    //				cout<<"count="<<ncount<< " nc1="<< ncounthr1 << " nc2=" <<ncounthr2 << " nc3="<< ncounthr3 << ", zcount= " << zcount<< ", time slot "<< ibin << endl;
		    // cout<<"count "<<ncount<< "zcount" << zcount<< " Nlarg=" << Nlarg<<endl;
						



		    for(Int_t kEE1=0; kEE1<NZ;kEE1++)
		      {
			
			Double_t eff =0;
			if(myCounters[Form("track%d",kEE1)]>0)  eff =Double_t(myCounters[Form("trackf%d",kEE1)])/Double_t(myCounters[Form("track%d",kEE1)]); 

			Double_t multip=1;
			if(myCounters[Form("track%d",kEE1)]>0) multip =Double_t(myCounters[Form("multip%d",kEE1)])/Double_t(myCounters[Form("trackf%d",kEE1)]);
 








 if(NEvent%10==0&&detail)
			 {
				 Float_t fac=1.;
								 
				 
				 cout << "number of tracks " <<kEE1 << " " <<Double_t(myCounters[Form("trackf%d",kEE1)]) << " " << Double_t(myCounters[Form("track%d",kEE1)]) << endl;
				 
				 cout << "Eff  "  <<kEE1 << " " <<fac*Double_t(myCounters[Form("trackf%d",kEE1)])/Double_t(myCounters[Form("track%d",kEE1)]) << endl;
				 
				 cout << "multip  " << kEE1 << "    "  <<Double_t(myCounters[Form("multip%d",kEE1)])/Double_t(myCounters[Form("trackf%d",kEE1)]) << endl;
				 
				 
				 
				 
			
			 cout <<"layer" << kEE1 << " ";
			 cout << "Nhits="<< myCounters[Form("hits%d",kEE1)] << " ";
			 cout << "Time="<< (myCounters[Form("timemax%d",kEE1)]-myCounters[Form("timemin%d",kEE1)])*2/10000. << " ms " << endl;
			 if(myCounters[Form("timemax%d",kEE1)]-myCounters[Form("timemin%d",kEE1)]>0) cout << "noise="<< 542.5 * myCounters[Form("hits%d",kEE1)]/(myCounters[Form("timemax%d",kEE1)]-myCounters[Form("timemin%d",kEE1)])<<endl;


			 }
			




			if(eff>0)
			  {
					 
			  
			    ncounthrc1= (ncountp1[kEE1]/eff)+ncounthrc1;
			    ncounthrc2= (ncountp2[kEE1]/eff)+ncounthrc2;
			    ncounthrc3= (ncountp3[kEE1]/eff)+ncounthrc3;
			 

			  }

			//				 	 cout << "eff  " << eff << " multip3" << multip3<< endl; 

			if( myCounters[Form("trackf%d",kEE1)]>30)
			  {
			    myCounters[Form("trackf%d",kEE1)]=0;
			    myCounters[Form("track%d",kEE1)]=0;
			    myCounters[Form("multip%d",kEE1)]=0;
          
			  }

				 

	 
			 
				 
	 
			if(eff<=0)   
			  {
					 
			    ncounthrc1= ncountp1[kEE1]+ncounthrc1;
			    ncounthrc2= ncountp2[kEE1]+ncounthrc2;
			    ncounthrc3= ncountp3[kEE1]+ncounthrc3;
			 
			  }

 
			/*
			  if(myCounters[Form("trackf%d",KEE)]>100)
			  {
			  myCounters[Form("track%d",KEE)]=0;
			  myCounters[Form("multip%d",KEE)]=0;
			  myCounters[Form("trackf%d",KEE)]=0;
			  }
 

			*/ 

		      }

 if((NEvent%2)==0)
	 {
			 cout << "beam intensity="<< 5*1e6*Double_t(NEPH)/(maxtime-mintime)<< endl;

			 cout <<  " Number of particles =  " <<  NEPH <<  "     mintime=  " << mintime << " maxtime= " << maxtime << endl;

			 if(NEPH>0) cout <<  " Average number of hits=  " <<  NHITPH/NEPH << endl;

	 }
 
		    //cout << "test" << ncounthrc1<< "   "<<ncounthrc2<< "  "<<ncounthrc3<<endl;
		    //if(mu1>0) ncounthrc1=ncounthrc1/mu1;
		    //if(mu2>0) ncounthrc2=ncounthrc2/mu2;
		    //if(mu3>0) ncounthrc3=ncounthrc3/mu3;
 
 

		    myCounters["anaevent"]= myCounters["anaevent"]+1;
		    myCounters["thr1"]= myCounters["thr1"]+ncounthrc1;
		    myCounters["thr2"]= myCounters["thr2"]+ncounthrc2;
		    myCounters["thr3"]= myCounters["thr3"]+ncounthrc3;
		    //  cout<<" anaevent = " << myCounters["thr3"]<<endl;
							
		    if(myCounters["anaevent"]>0)
		      { 


			if(myCounters["mevent"]>0)
						 
			  {
			    mu1=float(myCounters["mu1"])/ myCounters["mevent"];
			    mu2=float(myCounters["mu2"])/ myCounters["mevent"];
			    mu3=float(myCounters["mu3"])/ myCounters["mevent"];

			    //		cout << "mevent= " << myCounters["mevent"] <<  "   mu1= "<<mu1 << "  mu2= "<< mu2<< "  mu3= "<<mu3<<endl; 
			  }

	


			// pi1=(float(myCounters["thr1"])/ myCounters["anaevent"]);
			//	pi2=(float(myCounters["thr2"])/ myCounters["anaevent"]);
			//	pi3=(float(myCounters["thr3"])/ myCounters["anaevent"]);

			pi1=ncounthrc1;
			pi2=ncounthrc2;
			pi3=ncounthrc3;

			/* 
			   if(myCounters["mevent"]>100)
			   {

			   myCounters["thr1"]=0;
			   myCounters["thr2"]=0;
			   myCounters["thr3"]=0;
			   myCounters["mu1"]=0;
			   myCounters["mu2"]=0;
			   myCounters["mu3"]=0;
			   myCounters["mevent"]=0;
			   myCounters["anaevent"]=0;
			   }
			*/
			// yacine add : nhit for each asics

		       
			/**
			   for(int ilayer = 0; ilayer < 47; ilayer++){
			   for(int asic_i=0; asic_i < 12; asic_i++){
			   for(int asic_j=0; asic_j < 12; asic_j++){

			   yacine_density_asic = h_density_asic[ilayer].GetBinContent(asic_i+1,asic_j+1); // AP: density var.
			   //			      HistDensity.Fill(yacine_density_asic);
			   }
			   }
			   }
			**/

			//******************
			// Calibration from str line fit for run 715751:
			//	ncounthr1 = (ncounthr1-(-1.17971*SpillEvtTime*200/(pow(10.0,9)))); // time in sec.
			//		ncounthr2 = (ncounthr2-(-2.02072*SpillEvtTime*200/(pow(10.0,9))));
			//	ncounthr3 = (ncounthr3-(-1.90554*SpillEvtTime*200/(pow(10.0,9))));

		       

			if (calibl) {
			  ncounthr1 = (ncounthr1-(vect_col[1]*SpillEvtTime*200/(pow(10.0,9)))); // time in sec.
			  ncounthr2 = (ncounthr2-(vect_col[2]*SpillEvtTime*200/(pow(10.0,9)))); // time in sec.
			  ncounthr3 = (ncounthr3-(vect_col[3]*SpillEvtTime*200/(pow(10.0,9)))); // time in sec.
			  //		  ncountdh1 = (ncountdh1-(vect_col[1]*SpillEvtTime*200/(pow(10.0,9)))); // time in sec.
			  //  ncountdh2 = (ncountdh2-(vect_col[2]*SpillEvtTime*200/(pow(10.0,9)))); // time in sec.
			  //  ncountdh3 = (ncountdh3-(vect_col[3]*SpillEvtTime*200/(pow(10.0,9)))); // time in sec.
			  //  ncountdl1 = (ncountdl1-(vect_col[4]*SpillEvtTime*200/(pow(10.0,9)))); // time in sec.
			  //  ncountdl2 = (ncountdl2-(vect_col[5]*SpillEvtTime*200/(pow(10.0,9)))); // time in sec.
			  //  ncountdl3 = (ncountdl3-(vect_col[6]*SpillEvtTime*200/(pow(10.0,9)))); // time in sec.
			}

			// AP: try to use a spill time calibration:
			if (calibr) {
			  //	if (SpillEvtTime<10000000) {
			  if (SpillEvtTime>=10000000&&SpillEvtTime<20000000) {
			    ncounthr1 = ncounthr1 * vect1[2];
			    ncounthr2 = ncounthr2 * vect1[6];
			    ncounthr3 = ncounthr3 * vect1[10];
			    //	    nhough1 = nhough1 * vect1[2];
			    // nhough2 = nhough2 * vect1[6];
			    //  nhough3 = nhough3 * vec1[10];
			  }
			  
			  if (SpillEvtTime>=20000000&&SpillEvtTime<30000000) {
			    ncounthr1 = ncounthr1 * vect1[3];
			    ncounthr2 = ncounthr2 * vect1[7];
			    ncounthr3 = ncounthr3 * vect1[11];
			    //  nhough1 = nhough1 * vect1[3];
			    //  nhough2 = nhough2 * vect1[7];
			    //  nhough3 = nhough3 * vect1[11];
			  }
			  
			  if (SpillEvtTime>=30000000&&SpillEvtTime<40000000) {
			    ncounthr1 = ncounthr1 * vect1[4];
			    ncounthr2 = ncounthr2 * vect1[8];
			    ncounthr3 = ncounthr3 * vect1[12];
			    //   nhough1 = nhough1 * vect1[4];
			    //  nhough2 = nhough2 * vect1[8];
			    //  nhough3 = nhough3 * vect1[12];
			  }
			  
			  if (SpillEvtTime>=40000000&&SpillEvtTime<50000000) {
			    ncounthr1 = ncounthr1 * vect1[5];
			    ncounthr2 = ncounthr2 * vect1[9];
			    ncounthr3 = ncounthr3 * vect1[13];
			    //  nhough1 = nhough1 * vect1[5];
			    //  nhough2 = nhough2 * vect1[9];
			    //  nhough3 = nhough3 * vect1[13];
			  }
			  
			}

		       
			Int_t nhitt = ncounthr1+ncounthr2+ncounthr3;

		 



			double aa1=0.0407388; double bb1 =0.0657203;double cc1 =0.107314;double aa2 =3.38241e-05; double bb2=-.72138e-05;double cc2=1.28111e-05;double aa3=-1.58769e-08; double bb3=-2.0741e-10;double cc3 =3.23606e-08;
			float Erec=(aa1+nhitt*aa2+nhitt*nhitt*aa3)*(ncounthr1)+(bb1+nhitt*bb2+nhitt*nhitt*bb3)*(ncounthr2)+(cc1+nhitt*cc2+nhitt*nhitt*cc3)*(ncounthr3);

			float ErecH=(aa1+nhitt*aa2+nhitt*nhitt*aa3)*(ncounthr1-nhough1-nis1)+(bb1+nhitt*bb2+nhitt*nhitt*bb3)*(ncounthr2-nhough2-nis2)+(cc1+nhitt*cc2+nhitt*nhitt*cc3)*(ncounthr3-nhough3-nis3)+.06*(nhough1+nhough2+nhough3)+.065*(nis1+nis2+nis3);



			float Erecc=(aa1+nhitt*aa2+nhitt*nhitt*aa3)*pi1+(bb1+nhitt*bb2+nhitt*nhitt*bb3)*pi2+(cc1+nhitt*cc2+nhitt*nhitt*cc3)*pi3;


			//cout<< " ratio1 = " <<  ncounthr1/pi1 << " ratio3= " << ncounthr3/pi3<<endl;
			//cout << "ENERGYYYYYY=  "<< Erec<< "  EnerH = "  << ErecH<< " ratio="<<ErecH/Erec<< endl;					
			//cout << "ENERGYCORR=  "<< Erecc<<endl;					
			float Ereccc=Erec;
			if(ncounthr1/pi1<ncounthr3/pi3)

			  {

			    Ereccc = Erec*(ncounthr3/pi3)/(ncounthr1/pi1);
			  }

			//cout << "ENERGYCORR2fois=  "<< Ereccc<<endl;					



		
			//cout << "anaevent=  " << myCounters["anaevent"] <<  "    pi1=  "<<pi1 << "    pi2= "<< pi2<< "   pi3= "<<pi3<<endl; 


			//cout << "anaevent=  " << myCounters["anaevent"] <<  "    ncounthr1=  "<<ncounthr1 << "    ncounthr2= "<< ncounthr2<< "   ncounthr3= "<<ncounthr3<<endl; 


			//cout << "anaevent=  " << myCounters["anaevent"] <<  "    nhough1=  "<<nhough1 << "    nhough2= "<< nhough2<< "   nhough3= "<<nhough3<<endl; 

			//cout << "anaevent=  " << myCounters["anaevent"] <<  "    nis1=  "<<nis1 << "    nis2= "<< nis2<< "   nis3= "<<nis3<<endl; 

			//				cout <<  "XGM =" <<XGM << "  XGMu= " << XGMu << " YGM= "  <<  YGM << " YGMu= "<< YGMu <<endl;     

			//	cout << "pi1n= "<<ncounthrc1 << "  pi2n= "<< ncounthrc2<< "  pi3n= "<<ncounthrc3<<endl; 



	 


		

		
		
   

			//cout << "anaevent=  " << myCounters["anaevent"] <<  "    Ndeb=  "<<Ndeb << "    Ndoub= "<< Ndoub<< "   Ncern= "<<Ncern<< "  Nfin=" << Nfin<< "Nend= " << Nend<< " ncl=  "<< ncl<< " RMAX= " << RMAX<< " KMAX = " << KMAX*2.8<<endl; 

			//		ntup1.Fill(ncounthr1,ncounthr2,ncounthr3,pi1,pi2,pi3,ncount,Nedge,zcount,NKFIN,KMAX,NKMAX,RMAX,Nhole,Ndeb);
			//,Nfin,XGM,YGM,Ncern,dNdoub,Nlarg);
    
			//		ntup1.Fill(nfin,ndeb,thr,ncounthr1,ncounthr2,ncounthr3,thrp,ncounthrc1,ncounthrc2,ncounthrc3,pi1,pi2,pi3,ncount,zcount);
 
			//	cout<<"hereeeeeeeeee"<<IBIN0-IBINBef<<endl;
 
			dt=IBIN0-IBINBef;
			dx=fabs(XGM-XGMBef);
			dy=fabs(YGM-YGMBef);
			//cout<<"hereeeeeeeeee"<<XGM<< " "<<XGMBef<<endl;
			//cout<<"hereeeeeeeeee"<<dt<< " "<<dx<<" "<<dy<<endl;
			sdhcal->Fill();

			XGMBef=XGM;
			YGMBef=YGM;

			float ratio = ncount/zcount;
			HT1.Fill(ncounthrc1);
			HT2.Fill(ncounthrc2);
			HT3.Fill(ncounthrc3);
			HT4.Fill(ncount);
			HT5.Fill(zcount);
			HT6.Fill(ratio);


		      }
		    cern=false;
		    		    
		    //  cout << "ncount = " << ncount << endl;
		    //  cout << "Ndeb = " << Ndeb << endl;
		    //  cout << "zcount = " << zcount << endl;
		    //  cout << "Nedge = " << Nedge << endl;
		    //  cout << "XGM = " << XGM << endl;
		    // cout << "YGM = " << YGM << endl;
		    // cout << "Neutcosm = " << Neutcosm << endl;
		    // cout << "Nlarg = " << Nlarg << endl;
		    // cout << "****************************" << endl;
		    

				//+ if(draw&&ncount>70&&(Ndeb>5||(zcount>30&&Ndeb>-1))&&float(Nedge)/ncount<.02&&(float(ncount)/zcount>2.2)&&(Ndeb<50)&&abs(XGM-50)<20&&abs(YGM-50)<20&&Neutcosm<1&&(float(Nlarg)/zcount>.3&&zcount>44)&&dt>10000&&doub==1) // draw option
		    //  if(draw&&(Ndeb>5||(zcount>30&&Ndeb>-1))&&float(Nedge)/ncount<.02&&(float(ncount)/zcount>2.2)&&(float(ncount)/zcount<2.6)&&(Ndeb<50)&&fabs(XGM-50)<20&&fabs(YGM-50)<20&&Neutcosm<1&&(zcount<42||float(Nlarg)/zcount>.3)&&float(Nlarg)/zcount>.2) // draw option
				//		  if (draw&&(float(ncount)/zcount<2.2))
		      if (draw&&zcount>7)
		      {
    
			c1->cd();
   
				     	
			H3D0.Draw();
			H3D1.Draw("same");
			H3D2.Draw("same");
			H3DH.Draw("same");

			//		 H3DI.Draw("same");

					 
			c1->Update(); 	
			c1->Modified();	        

			//			c1->SaveAs("F80Sep.pdf");
					



    


			//cout << "time slot" << IBIN0 <<endl;
			//getchar();
	  
	  
			c2->Clear();					 			  	
			c2->Divide(1,2);
					 
	  
			//	  c2->cd(1);
			//	  H3D0.Draw();
			// H3D1.Draw("same");
			//H3D2.Draw("same");
	  
	  
	  
			c2->cd(1);
			H10.Draw();
			H11.Draw("same");
			H12.Draw("same");
			H12H.Draw("same");

			H10.GetXaxis()->SetTitle("Z (cm)");
			H10.GetYaxis()->SetTitle("X (cm)");


  
			TText *test1=new TText(10,85," CALICE Preliminary");
			test1->Draw();


 


			//  TLegend *legend = new TLegend(0.553482,0.653319,0.950072,0.948683,NULL,"brNDC");
  
			gStyle->SetOptStat(0);
			//  legend->SetTextAlign(22);
			// legend->SetTextSize(0.055);
			//   legend->AddEntry(H10,"Hit1","LPF");
			// legend->AddEntry(H11,"Hit2","LPF");
			//legend->AddEntry(H12,"Hit3","LPF");


  
			c2->cd(2);

			H20.Draw();
			H21.Draw("same");
			H22.Draw("same");
			H22H.Draw("same");

			//						 TText *test2=new TText(10,85," CALICE Preliminary");
			// test2->Draw();

						
			//H1->Draw(); 
			H20.GetXaxis()->SetTitle("Z (cm)");
			H20.GetYaxis()->SetTitle("Y (cm)");


			c2->Update(); 	
			c2->Modified();
			sleep(10);
	  
		  
			//
			//	 c2->SaveAs("muon80GeV.pdf");
					
			if(waiting) getchar();
   
			/*	 

			H3D0->Delete();
			H3D1->Delete();
			H3D2->Delete();
			*/

 
			/*
			  H10->Delete();
			  H11->Delete();
			  H12->Delete();

	
			  H20->Delete();
			  H21->Delete();
			  H22->Delete();
			*/
		      }

		    //   ncount=0;



					    
		  } //fin selection pion, 
					 
 


	      }
		
			 
	  }
	H.Delete();
	HT1.Delete();
	HT2.Delete();
	HT3.Delete();
      }
	
      //	} // 2 time slots














	   
    }

  // cout << "DHCALEvent_ : " << DHCALEvent_ << endl;
  // cout << "ncount : " << ncount << endl;

  //  delete [] densities;

  return kTRUE;
}
	
	
void DHCAL_Tree_selector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

}

void DHCAL_Tree_selector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  float_t fac =1.;
	
  for(Int_t k=0;k<50;k++)
    {
		
		
		
	 
						
		
		
		
		
      //  cout << k << "  "<< "track"<<k<<" " <<  myCounters[Form("track%d",k)]<< "trackf"<<k<< " "<<  myCounters[Form("trackf%d",k)] <<endl;

      Double_t eff =fac*Double_t(myCounters[Form("trackf%d",k)])/Double_t(myCounters[Form("track%d",k)]); 
      //  cout <<" EFF" << k<<"= "<< eff <<endl;
      //  cout << "Delta"<<k<<"= " << sqrt(eff*(1.-eff)/Double_t(myCounters[Form("track%d",k)]))<<endl;
      // cout << "multip" << k << "      " <<  Double_t(myCounters[Form("multip%d",k)])/Double_t(myCounters[Form("trackf%d",k)]) << endl;
		
		

      //      cout << "eff"<<k<<" " <<  myCounters[Form("eff%d",k)]<< "err"<<k<< " "<<  myCounters[Form("err%d",k)] <<endl;
 
    
		
    }
  //       	for (map<string,int>::iterator it=myCounters.begin(); it != myCounters.end(); it++)

  //	cout << it->first << " " << it->second << endl;

  //  	out.close();
	
  //		H2->Delete();
  //sdhcal->Write();
  f->Write();
  // f->Delete();
}



