#include <TRandom.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include <string>
#include <map.h>
#include <iostream>
#include <fstream>
#include <Riostream.h>
#include <iostream>
#include <fstream>
#include <TF1.h>
#include <TLegend.h>

TNtuple* event  = new TNtuple("event","event ","N1:N2:N3");


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



void minimisation_3seuilsTI3compress_HT() 
{   
	CaliceStyle();
   ofstream calibration; calibration.open("calibration.txt");       
   TString outputFileName= "minim.root";
   TFile *target  = new TFile(outputFileName,"RECREATE" );	
   TH1F** GeVH=new TH1F * [120];
   TH1F** N_Hit=new TH1F * [120];
   float Min_a1[120]; float Min_b1[120]; float Min_c1[120];float Min_a2[120]; float Min_b2[120]; float Min_c2[120];float Min_a3[120];float Min_b3[120];float Min_c3[120];float Min_d1[120];float Min_e1[120];

   for(int i=0;i<120;i++){Min_a1[i]=0;Min_b1[i]=0;Min_c1[i]=0;Min_a2[i]=0;Min_b2[i]=0;Min_c2[i]=0;Min_a3[i]=0;Min_b3[i]=0;Min_c3[i]=0; Min_d1[i]=0; Min_e1[i]=0;}
   float cal1[20]; float cal2[20]; float cal3[20];
   char file_name[1000];
   Int_t Ene[20];
   /*  	 
   Ene[0]=5;    
   Ene[1]=10;    
   Ene[2]=20;
   Ene[3]=60;    
   Ene[4]=80;    
   */ 	
   
   Ene[0]=20;    
   Ene[1]=30;    
   Ene[2]=40;
   Ene[3]=50;    
   Ene[4]=60;    
   Ene[5]=70;    
   Ene[6]=80;    
   	  
   /*  
   Ene[0]=10;    
   Ene[1]=30;
   Ene[2]=50;    
   Ene[3]=70;    
   */
   Int_t first_plate=0;  
   Int_t ppoint =6;

      ofstream coool("temporaire.txt");



	 for(kk=first_plate;kk<=ppoint;kk++)
   {
		 

  
		 for( Int_t jj=0;jj<1;jj++)

			 		 {
						 sprintf(file_name,"/gridgroup/ilc/petr/data-2014/output-TBcal-%dGeV.root",Ene[kk]); cout<< kk<<"----------------"<<file_name<<" -----  Energy= "<<Ene[kk]<<endl;
      TFile input(file_name);
      //--------------------------------------------------- 
      TTree *sdhcal = (TTree *)  input.Get("sdhcal");
      //---------------------------------------------------      
         
			Int_t ncount,ncounthr1,ncounthr2,ncounthr3,Ndeb,Nedge,Ndoub,Nfin, Nhole, zcount, ncount, NKMAX, dt, Neutcosm, NKFIN,nis1, nis2, nis3,nhough1,nhough2,nhough3, Nend,ncl , Nlarg ;
			Float_t RMAX,XGM,YGM, pi1, pi2, pi3, pi11, pi21, pi31;

			
      sdhcal->SetBranchAddress("ncounthr1", &ncounthr1);   sdhcal->SetBranchAddress("ncounthr2", &ncounthr2);   sdhcal->SetBranchAddress("ncounthr3", &ncounthr3);  sdhcal->SetBranchAddress("Ndeb", &Ndeb);  sdhcal->SetBranchAddress("Nedge", &Nedge);  sdhcal->SetBranchAddress("Ndoub", &Ndoub);  sdhcal->SetBranchAddress("RMAX", &RMAX);   sdhcal->SetBranchAddress("XGM", &XGM); sdhcal->SetBranchAddress("YGM", &YGM); sdhcal->SetBranchAddress("Nfin", &Nfin);   sdhcal->SetBranchAddress("zcount", &zcount);    sdhcal->SetBranchAddress("ncount", &ncount);    sdhcal->SetBranchAddress("NKMAX", &NKMAX);     sdhcal->SetBranchAddress("dt", &dt);    sdhcal->SetBranchAddress("Neutcosm", &Neutcosm);      sdhcal->SetBranchAddress("NKFIN", &NKFIN);   sdhcal->SetBranchAddress("Nend", &Nend);   sdhcal->SetBranchAddress("ncl", &ncl);  sdhcal->SetBranchAddress("nis1", &nis1);  sdhcal->SetBranchAddress("nis2", &nis2); sdhcal->SetBranchAddress("nis3", &nis3); sdhcal->SetBranchAddress("nhough1", &nhough1); sdhcal->SetBranchAddress("nhough2", &nhough2);sdhcal->SetBranchAddress("nhough3", &nhough3);  sdhcal->SetBranchAddress("pi1", &pi1); sdhcal->SetBranchAddress("pi2", &pi2); sdhcal->SetBranchAddress("pi3", &pi3); sdhcal->SetBranchAddress("Nlarg", &Nlarg); 



      cout <<sdhcal->GetEntries()<<" "<<Ene[kk]<<endl;
    

     
      for(int i=1;i<=sdhcal->GetEntries();i++)
      {
         sdhcal->GetEntry(i);
				 //				 cout<<zcountp<< "   "<<Ene[kk]<<endl;
	 


	 Float_t nhitt= (ncounthr1)+ncounthr2+ncounthr3;
				 //   Float_t nhitt= ncounthr1-0*nis1+(ncounthr2-0*nis2)+(ncounthr3-0*nis3);

	 //	 if(Ndeb>3&&Nedge/ncount<.05&&ncount/zcount>2&&Ndoub<1&&RMAX/zcount>.2&&Ndeb<45&&abs(XGM-50)<15&&abs(YGM-50)<15&&Nhole<1)
	 if( 

			//(Ndeb>4||(zcount>30&&Ndeb>1))&&Nedge/ncount<.02&&ncount/zcount>2&&RMAX/zcount>.2&&Ndoub<1&&(Ndeb<20)&&abs(XGM-50)<30&&abs(YGM-50)<30&&Nhole<1

			//(Ndeb>5||(zcount>30&&Ndeb>-1&&Ndeb<6))&&float(Nedge)/ncount<.02&&(float(ncount)/zcount>2&&RMAX/zcount>.2)&&Ndoub<1&&(Ndeb<10)&&abs(XGM-50)<30&&abs(YGM-50)<30&&Nhole<1&&Neutcosm<1&&float(Nlarg)/zcount>.2

			//	 (Ndeb>4||(zcount>30&&Ndeb>-1))&&float(Nedge)/ncount<.02&&(Ndeb<50)&&abs(XGM-50)<20&&abs(YGM-50)<20&&Neutcosm<1&&(float(Nlarg)/zcount>.25)&&dt>10000&&Ndoub<1
	 //&&float(Nend)/ncount<.05


			// (Ndeb>4||(zcount>30&&Ndeb>-1))&&float(Nedge)/ncount<.02&&(Ndeb<50)&&abs(XGM-50)<20&&abs(YGM-50)<20&&Neutcosm<1&&(float(Nlarg)/zcount>.25)&&dt>10000&&Ndoub<1
		 


	    //	    (nhough1>5)&&float(Nedge)/ncount<.02&&(float(ncount)/zcount>2.2)&&(Ndeb<50)&&abs(XGM-50)<50&&abs(YGM-50)<50&&Neutcosm<1&&float(Nlarg)/zcount>.2&&Ndoub<1&&Ndeb>20
	    //    (Ndeb>4||(zcount>30&&Ndeb>-1))&&float(Nedge)/ncount<.02&&(float(ncount)/zcount>2.2)&&(Ndeb<50)&&(zcount<42||float(Nlarg)/zcount>.2)&&float(Nlarg)/zcount>.2&&Ndoub<1&&float(Nlarg)/zcount<.95

	    //	    (Ndeb>4||(zcount>30&&Ndeb>-1))&&float(Nedge)/ncount<.02&&(float(ncount)/zcount>2.2)&&(Ndeb<50)&&abs(XGM-50)<20&&abs(YGM-50)<20&&Neutcosm<1&&(zcount<42||float(Nlarg)/zcount>.2)&&float(Nlarg)/zcount>.2&&Ndoub<1&&float(Nlarg)/zcount<.95

	    (Ndeb>9||(zcount>35&&Ndeb>-1))&&float(Nedge)/ncount<.02&&(float(ncount)/zcount>2.2)&&(Ndeb<55)&&abs(XGM-50)<20&&abs(YGM-50)<20&&Neutcosm<1&&(zcount<47||float(Nlarg)/zcount>.25)&&float(Nlarg)/zcount>.2&&float(Nlarg)/zcount<.95&&Ndoub<1

	    )
	   coool<< ncounthr1-nhough1<<"  "<<ncounthr2-nhough2<<"  "<<ncounthr3-nhough3<<"  "<<nhough1+nhough2+nhough3<<"   "<<Ene[kk]<<endl;


       
			 }

      
			 }
					 }
	 coool.close();

	 cout<<"End of first step"<<endl;
	 int koko = minimizer(1,1,1,1,1,1,1,1,1,1);


   	  
	 Ene[0]=20;    
	 Ene[1]=30;    
	 Ene[2]=40;
	 Ene[3]=50;    
	 Ene[4]=60;    
	 Ene[5]=70;    
	 Ene[6]=80;    
   	  

	 int kkk=0;

	 first_plate=0;  
	 ppoint =6;
	 
for(kkk=first_plate;kkk<=ppoint;kkk++)
   {
          char hName[200];
		 for( Int_t jjj=0;jjj<1;jjj++)

			 		 {


				 sprintf(file_name,"/gridgroup/ilc/petr/data-2014/output-TBcal-%dGeV.root",Ene[kkk]); cout<<"----------------"<<file_name<<" -----  Energy= "<<Ene[kkk]<<endl;
      TFile input(file_name);
      //--------------------------------------------------- 
      TTree *sdhcal = (TTree *)  input.Get("sdhcal");
      //---------------------------------------------------      
         
    sdhcal->SetBranchAddress("ncounthr1", &ncounthr1);   sdhcal->SetBranchAddress("ncounthr2", &ncounthr2);   sdhcal->SetBranchAddress("ncounthr3", &ncounthr3);  sdhcal->SetBranchAddress("Ndeb", &Ndeb);  sdhcal->SetBranchAddress("Nedge", &Nedge);  sdhcal->SetBranchAddress("Ndoub", &Ndoub);  sdhcal->SetBranchAddress("RMAX", &RMAX);   sdhcal->SetBranchAddress("XGM", &XGM); sdhcal->SetBranchAddress("YGM", &YGM); sdhcal->SetBranchAddress("Nfin", &Nfin);  sdhcal->SetBranchAddress("zcount", &zcount);    sdhcal->SetBranchAddress("ncount", &ncount);    sdhcal->SetBranchAddress("NKMAX", &NKMAX);     sdhcal->SetBranchAddress("dt", &dt);    sdhcal->SetBranchAddress("Neutcosm", &Neutcosm);     sdhcal->SetBranchAddress("NKFIN", &NKFIN);   sdhcal->SetBranchAddress("Nend", &Nend);   sdhcal->SetBranchAddress("ncl", &ncl);  sdhcal->SetBranchAddress("nis1", &nis1);  sdhcal->SetBranchAddress("nis2", &nis2); sdhcal->SetBranchAddress("nis3", &nis3); sdhcal->SetBranchAddress("nhough1", &nhough1); sdhcal->SetBranchAddress("nhough2", &nhough2);sdhcal->SetBranchAddress("nhough3", &nhough3); sdhcal->SetBranchAddress("pi1", &pi1); sdhcal->SetBranchAddress("pi2", &pi2); sdhcal->SetBranchAddress("pi3", &pi3);  sdhcal->SetBranchAddress("Nlarg", &Nlarg);  

 
   
    sprintf(hName,"%s%d%s","toto",Ene[kkk],"GeV");
      
	//	if(jjj==0)TH1F *temp = new TH1F("temp", "temp", 560, 0, 140.); temp->SetDirectory(0); GeVH[kkk]=temp;
    if(jjj==0)TH1F *temp = new TH1F(hName, "temp", 1120, 0, 140.); temp->SetDirectory(0); GeVH[kkk]=temp;
    if(jjj==0)  TH1F *temp2 = new TH1F("temp2", "temp2", 2000, 0, 2000.); temp2->SetDirectory(0); N_Hit[kkk]=temp2;
            

		
      
      ifstream ff; ff.open("temp2.txt");
      float aa1,aa2,aa3,bb1,bb2,bb3,cc1,cc2,cc3,dd1,ee1;  
      while(ff>>aa1>>aa2>>aa3>>bb1>>bb2>>bb3>>cc1>>cc2>>cc3>>dd1){Min_a1[kkk]=aa1; Min_b1[kkk]=bb1; Min_c1[kkk]=cc1;Min_a2[kkk]=aa2; Min_b2[kkk]=bb2; Min_c2[kkk]=cc2; Min_a3[kkk]=aa3; Min_b3[kkk]=bb3; Min_c3[kkk]=cc3;Min_d1[kkk]=dd1;}

      calibration<<aa1<<"  "<<aa2<<"  "<<aa3<<"  "<<bb1<<"  "<<bb2<<"  "<<bb3<<"   "<<cc1<<"  "<<cc2<<"  "<<cc3<<"  "<<dd1<<"  "<<Ene[kk]<<endl;   
			ff.close();
		


			for(int ii=1;ii<=sdhcal->GetEntries();ii++)
			  {       
			    
			    
			    sdhcal->GetEntry(ii);

			      

			
			    Float_t nhitt=  ncounthr1+ncounthr2+ncounthr3;
			    //			 Float_t nhitt=  (ncounthr1-0*nis1)+(ncounthr2-0*nis2)+(ncounthr3-0*nis3);
			    //Float_t nhitt=  pi1+pi2+pi3;
			    Float_t eff=1.;
			    if(nis1>0) eff=ncounthr1/nis1; 

		
					
if(
	 //	 (Ndeb>4||(zcount>30&&Ndeb>1))&&Nedge/ncount<.02&&(ncount/zcount>2&&RMAX/zcount>.2)&&Ndoub<1&&(Ndeb<11)&&abs(XGM-50)<20&&abs(YGM-50)<20&&Nhole<1&&Neutcosm<1

	 //	 	 (Ndeb>5||(zcount>30&&Ndeb>-1))&&float(Nedge)/ncount<.02&&(float(ncount)/zcount>1.2&&RMAX/zcount>.2)&&Ndoub<1&&(Ndeb<21)&&abs(XGM-50)<20&&abs(YGM-50)<20&&Nhole<1&&Neutcosm<1&&float(Nlarg)/zcount>.2&&dt>10000&&float(Nend)/ncount<.1


	 //(Ndeb>5||(zcount>30&&Ndeb>-1&&Ndeb<15))&&(float(ncount)/zcount>1.2)&&(Ndeb<15)&&abs(XGM-50)<20&&abs(YGM-50)<20&&Neutcosm<1&&float(Nlarg)/zcount>.2&&Ndoub<1
	 // (Ndeb>5||(zcount>30&&Ndeb>-1))&&float(Nedge)/ncount<.02&&(float(ncount)/zcount>1.2)&&(Ndeb<10)&&abs(XGM-50)<20&&abs(YGM-50)<20&&Neutcosm<1&&float(Nlarg)/zcount>.2&&dt>10000&&float(Nend)/ncount<.15
  
 //	 (Ndeb>4||(zcount>30&&Ndeb>-1))&&float(Nedge)/ncount<.02&&(float(ncount)/zcount>1.2)&&Ndoub<1&&(Ndeb<6)&&abs(XGM-50)<20&&abs(YGM-50)<20&&Nhole<1&&Neutcosm<1&&Ndeb<10



 


	 //(Ndeb>4||(zcount>30&&Ndeb>-1))&&float(Nedge)/ncount<.02&&(float(ncount)/zcount>2)&&(Ndeb<20)&&abs(XGM-50)<50&&abs(YGM-50)<50&&Neutcosm<1&&float(Nlarg)/zcount>.2&&Ndoub<1

   //	 (nhough1>5)&&float(Nedge)/ncount<.02&&(float(ncount)/zcount>2.2)&&(Ndeb<50)&&abs(XGM-50)<50&&abs(YGM-50)<50&&Neutcosm<1&&float(Nlarg)/zcount>.2&&Ndeb>20
   //  (Ndeb>4||(zcount>30&&Ndeb>-1))&&float(Nedge)/ncount<.02&&(float(ncount)/zcount>2.2)&&(Ndeb<50)&&abs(XGM-50)<20&&abs(YGM-50)<20&&Neutcosm<1&&(zcount<42||float(Nlarg)/zcount>.2)&&float(Nlarg)/zcount>.2&&float(Nlarg)/zcount<.95

   (Ndeb>9||(zcount>35&&Ndeb>-1))&&float(Nedge)/ncount<.02&&(float(ncount)/zcount>2.2)&&(Ndeb<55)&&abs(XGM-50)<20&&abs(YGM-50)<20&&Neutcosm<1&&(zcount<47||float(Nlarg)/zcount>.25)&&float(Nlarg)/zcount>.2&&float(Nlarg)/zcount<.95&&Ndoub<1

 
)
		

	 

	{ 

    float  corr=1;
		if(ncounthr3/pi3> ncounthr1/pi1) corr = (ncounthr3/pi3)/( ncounthr1/pi1);
		float Erec=(aa1+nhitt*aa2+nhitt*nhitt*aa3)*(ncounthr1-nhough1)+(bb1+nhitt*bb2+nhitt*nhitt*bb3)*(ncounthr2-nhough2)+(cc1+nhitt*cc2+nhitt*nhitt*cc3)*(ncounthr3-nhough3)+dd1*(nhough1+nhough2+nhough3);
//Erec=Erec*corr;
         GeVH[kkk]->Fill(Erec);

			 }
				 if(



						Ndeb>2&&Nedge/ncount<.05&&ncount/zcount>2&&Ndoub<1&&RMAX/zcount>.2&&Ndeb<50&&abs(XGM-50)<15&&abs(YGM-50)<15 /*&&Nhole<3*/

 )
					 { float tot=nhitt;
         N_Hit[kkk]->Fill(tot);
					 }
      }

		 
		
					 }
   }


 cout<<aa1<<"  "<<aa2<<"  "<<aa3<<"  "<<bb1<<"  "<<bb2<<"  "<<bb3<<"   "<<cc1<<"  "<<cc2<<"  "<<cc3<<" " <<dd1<<"  "<<Ene[kkk]<<endl;   

 cout<<"double bestX1="<<aa1<<"; double bestY1 ="<<bb1<<";double bestZ1 ="<<cc1<<";double bestX2 ="<<aa2<<"; double bestY2="<<bb2<<";double bestZ2="<<cc2<<";double bestX3="<<aa3<<"; double bestY3="<<bb3<<";double bestZ3 ="<<cc3<<";"<<";double bestW1 ="<<dd1<<";"<< endl;

   cout<<"==========================  Lecture des .root terminee"<<endl;
   
   const int ttt=ppoint+1-first_plate;
   float Mean[ttt];float RMS[ttt];float energie[ttt];float z[ttt];float z_error[ttt];float z_h_error[ttt];
   float rapport[ttt];float z[ttt]; float energie[ttt]; float z_h[ttt];
   float Mean_nhit[ttt];float RMS_nhit[ttt]; float aa1[ttt]; float bb1[ttt]; float cc1[ttt]; float Er[ttt];
   
	 //   int ee=Ene[first_plate];
	 gStyle->SetOptStat(111111);
 TCanvas *cv3 = new TCanvas("cv3", "cv3",1000,600); cv3.Divide(1,1);   
 cv3.cd(1);


TF1 *ff3 = new TF1("ff3","[0]+[1]*x+[2]*x*x",0,1400);
 ff3->GetHistogram()->GetXaxis()->SetTitle("Nhit");
 ff3->GetHistogram()->GetYaxis()->SetTitle("Calibration coefficients");
	 ff3->SetParameters(cc1,cc2,cc3);
 ff3->SetLineColor(2);ff3->Draw(""); 
 ff3->GetHistogram()->GetXaxis()->SetTitle("Nhit");
 ff3->GetHistogram()->GetYaxis()->SetTitle("Calibration coefficients");
 ff3->GetYaxis()->SetRangeUser(0.,.3);
   TF1 *ff1 = new TF1("ff1","[0]+[1]*x+[2]*x*x",0,1400);
	 ff1->SetParameters(aa1,aa2,aa3);
 ff1->SetLineColor(3);ff1->Draw("same"); 

TF1 *ff2 = new TF1("ff2","[0]+[1]*x+[2]*x*x",0,1400);
	 ff2->SetParameters(bb1,bb2,bb3);
 ff2->GetXaxis()->SetTitle("Nhit");
ff2->SetLineColor(4);ff2->Draw("same"); 

TF1 *ff4 = new TF1("ff4","0.045*p0",0,1400);
 ff4->GetXaxis()->SetTitle("Nhit");
ff4->SetLineColor(5);ff4->Draw("same"); 

   test=new TText(50,.25," ");
   test->Draw();
   gStyle->SetOptTitle(0);
   TLegend *legend= new TLegend(0.6,0.7,0.95,0.95,NULL,"brNDC");
   legend->SetTextAlign(22);
   legend->SetTextSize(0.055);
   legend->AddEntry(ff1,"#alpha","L");
   legend->AddEntry(ff2,"#beta","L");
   legend->AddEntry(ff3,"#gamma","L");
   legend->AddEntry(ff4,"C","L");
   

	 //   legend->SetHeader("Calibration coefficients");
   legend->SetFillColor(kWhite);
   legend->Draw();


 cv3->Print("evolution.jpg");


   TCanvas *cv1 = new TCanvas("cv1", "cv1",1000,1000); 
	 cv1.Divide(3,4); 
	 cv1.cd(1);

gStyle->SetOptFit(1);
gStyle->SetOptStat(111);

  // AP: for Sept systematics:
 double sigma_CB[20];
 sigma_CB[0]=1.24;sigma_CB[1]=1.82;sigma_CB[2]=2.43;sigma_CB[3]=2.65;sigma_CB[4]=3.28;sigma_CB[5]=3.74;
 sigma_CB[6]=4.4;sigma_CB[7]=5.09;sigma_CB[8]=5.37;sigma_CB[9]=6.15;sigma_CB[10]=5.52;sigma_CB[11]=8.06;
 double dE_CB[20];
 dE_CB[0]=0.54;dE_CB[1]=0.18;dE_CB[2]=0.64;dE_CB[3]=0.2;dE_CB[4]=0.79;dE_CB[5]=-0.5;dE_CB[6]=1.34;dE_CB[7]=1.23;
 dE_CB[8]=0;dE_CB[9]=1.;dE_CB[10]=-1.73;dE_CB[11]=-2.3;
 double E_CB[20];
 E_CB[0]=5.54;E_CB[1]=7.68;E_CB[2]=10.64;E_CB[3]=15.2;E_CB[4]=20.79;E_CB[5]=24.50;E_CB[6]=31.34;E_CB[7]=41.23;
 E_CB[8]=50.;E_CB[9]=61.;E_CB[10]=68.27;E_CB[11]=77.70;

   for(int ii=0;ii<ppoint+1-first_plate;ii++)
	 // for(int ii=0;ii<11;ii++)
		 {  
			 int jj=ii+first_plate;    
			 cout<<"------------------> Fit de l'energie "<<Ene[jj]<<endl;
			 // added by yacine 
			 TF1 *fitfun = new TF1("fitfun","gaus");
			 // GeVH[jj]->Fit("fitfun","Q");
			 float meany =GeVH[jj]->GetMean();
       float rmsy =GeVH[jj]->GetRMS();
			 			 	  GeVH[jj]->Fit("fitfun","Q","",meany-rmsy,meany+rmsy);
								// GeVH[jj]->Fit("fitfun","Q","",meany-10,meany+20);

								float centre=fitfun->GetParameter(1);
								float ecart=fitfun->GetParameter(2);

							  GeVH[jj]->Fit("fitfun","Q","",centre-1.5*ecart,centre+1.5*ecart);

								//	if(centre>=35) GeVH[jj]->Fit("fitfun","Q","",centre-2.*ecart,centre+2*ecart);

// GeVH[jj]->Fit("fitfun","Q");
			
			//
      //GeVH[jj]->Fit("gaus","Q");

			 Mean[ii]=fitfun->GetParameter(1);RMS[ii]=fitfun->GetParameter(2);energie[ii]=Ene[jj]; 
		


	 if(energie[ii]==8) energie[ii]=(7.5/8)*energie[ii]; 		 	 
	 rapport[ii]=(Mean[ii]-energie[ii])/energie[ii];   
     z[ii]=RMS[ii]/energie[ii];
		 z_h[ii]=RMS[ii]/Mean[ii];
			 
			 //	TH1F *tmptmp = new TH1F("tmptmpt","temptemp",141,0,140);tmptmp = GeVH[ii];
			 

			 cout <<"mean= " <<Mean[ii]<< "  rms="<<RMS[ii]<< "  rapport="<<rapport[ii]<<endl; 
			 //	 Er[ii]= fitfun->GetParError(1);
			 // AP: systematics:
			 Er[ii]= fitfun->GetParError(1);
			 double meanErr=fitfun->GetParError(1);
			 double se=fitfun->GetParError(2);
			 double sigmaErr=fitfun->GetParError(2) + (sigma_CB[ii]-RMS[ii]);
			 double err_reso = sqrt((sigmaErr/Mean[ii])**2 + (meanErr*RMS[ii]/Mean[ii]**2)**2);
			 double dE_G = Mean[ii] - energie[ii];
			 double dE = dE_CB[ii] - dE_G;
			 double meanErrTot=sqrt(meanErr**2 + (Mean[ii]-E_CB[ii])**2); // Fit and sys errors included
			 double err_mean = sqrt((meanErr/energie[ii])**2 + (0.01*Mean[ii]/energie[ii])**2 + dE**2);
			 z_h_error[ii]=z[ii]*((fitfun->GetParError(2)/fitfun->GetParameter(2))+(fitfun->GetParError(1)/fitfun->GetParameter(1))); 
			 z_error[ii]=(fitfun->GetParError(2)/energie[ii]) ;
			 


			 GeVH[jj]->SetLineWidth(2);GeVH[jj]->(3001);GeVH[jj]->SetFillColor(ii+2);
			 
			 N_Hit[jj]->Fit("gaus","Q");Mean_nhit[ii]=gaus->GetParameter(1);RMS_nhit[ii]=gaus->GetParameter(2);
			 //		 z_h[ii]=RMS_nhit[ii]/Mean_nhit[ii];

			 
			 cout <<"meanH= " <<Mean_nhit[ii]<< "  rmsH="<<RMS_nhit[ii]<< "  rapport="<<z_h[ii]<<endl;
			 cout << "\033[1;31mErrors for the energy= \033[0m"; cout << energie[ii]; 
			 cout << ": Err on resol= " << err_reso*100 << " , Err on Mean= " << err_mean*100 << " , Err on Ereco= " << meanErrTot << endl;
			 //      cout << "se= " << se << endl;
			 //	 cout << "dE_G= " << dE_G << "  dE_CB[ii]= " << dE_CB[ii] << endl; 
			 //	 cout << "sigmaErr= " << sigmaErr << " MeanErr[ii]= " << meanErr << " Ene[jj]= " << Ene[jj] << endl;
			 
			 delete fitfun;
		 }
   
	 // ee=first_plate;
	 for(int ii=0;ii<ppoint+1-first_plate;ii++)
	 
		 {
			 cv1.cd(ii+1);    
			 int jj=ii+first_plate;		
	 
			 GeVH[jj]->Draw();
			 target->cd();
                         GeVH[jj]->Write(); 

                 }
	 
	   target->Write();
           target->Close();  
	   cv1->Print("Ene-spec.jpg");
	 
   cv1->Delete();
	 
	  TCanvas *cv2 = new TCanvas("cv2", "cv2",1000,600); cv2.Divide(2,1);   
 
   cv2.cd(1);
   TGraphErrors *gr1 = new TGraphErrors(ttt,energie,Mean,0,Er); 
   gr1->Draw("A*"); gr1->SetTitle("Linearity (Mean_fit : energy)");
   TF1 *fa1 = new TF1("fa1","x",0,100); fa1->SetLineColor(2);fa1->Draw("same"); 
   test=new TText(5,70," Very Preliminary");
   test->Draw();
   
 gStyle->SetOptStat(1);
   cv2.cd(2);
   TGraphErrors *gr2 = new TGraphErrors(ttt,energie,z_h,0,z_h_error); 
	 gr2->GetYaxis()->SetRangeUser(0.,.4); 
   gr2->Draw("A*"); gr2->SetTitle("Resolution (sigma_fit/Energy) : energy");
   gr2->GetXaxis()->SetTitle("Energy(GeV)"); gr2->GetYaxis()->SetTitle("Resolution(%)"); gStyle->SetOptFit(kTRUE);
	   TF1 *fun_res1 = new TF1("fun_res1","sqrt(  (([0]/100.)**2/x)+(([1]/100.)**2))",1.,100.);
	 //  TF1 *fun_res1 = new TF1("fun_res1","sqrt(  (([0]/100.)**2/x))",1.,100.);
		  fun_res1->SetLineColor(1); fun_res1->SetLineStyle(1); fun_res1->SetLineWidth(2.2); gr2->Fit("fun_res1","r");
		 double p0_res1=fun_res1->GetParameter(0); double p1_res1=fun_res1->GetParameter(1); double p2_res1=fun_res1->GetParameter(2);  
       

test1=new TText(5,.35,"  Very Preliminary");
test1->Draw();


       
	 //	  TGraphErrors *gr3 = new TGraphErrors(ppoint,energie,z_h,0,0); gr3->SetMarkerColor(2); gr3->SetMarkerStyle(20);gr3->Draw("p*"); 
   
 
   cv2->Print("Ene-Res.jpg");
	 //	  cv1->Delete();
   //===================================================================================================
   //===================================================================================================
   //===================================================================================================
}

//=============================================================


double minimizer(double n1, double n2, double n3, double n4, double n5, double n6, double n7,double n8, double n9, double n10) 
{
   TFitter* minimizer = new TFitter(9);
   double p1 = -1; minimizer->ExecuteCommand("SET PRINTOUT",&p1,9);   
   minimizer->SetFCN(minuitFunction);
   minimizer->SetParameter(0,"X1",9,2,0,0);
   minimizer->SetParameter(1,"Y1",9,2,0,0);
   minimizer->SetParameter(2,"Z1",9,2,0,0);
   minimizer->SetParameter(3,"X2",9,2,0,0);
   minimizer->SetParameter(4,"Y2",9,2,0,0);
   minimizer->SetParameter(5,"Z2",9,2,0,0);
   minimizer->SetParameter(6,"X3",9,2,0,0);
   minimizer->SetParameter(7,"Y3",9,2,0,0);
   minimizer->SetParameter(8,"Z3",9,2,0,0);
   minimizer->SetParameter(9,"W1",9,2,0,0);
  
	 // minimizer->ExecuteCommand("SIMPLEX",0,0);
	 // minimizer->ExecuteCommand("MIGRAD",0,0);



	


 



// c'est l'option de base
   //double bestX1=0.0392888; double bestY1 =0.0748402;double bestZ1 =0.113314;double bestX2 =3.39241e-05; double bestY2=-3.6138e-06;double bestZ2=1.28111e-05;double bestX3=-1.68769e-08; double bestY3=-1.50741e-10;double bestZ3 =3.23606e-08;;double bestW1 =0.033;









// better resolution
//   double bestX1=0.0404481; double bestY1 =0.076236;double bestZ1 =0.114274;double bestX2 =3.63877e-05; double bestY2=-2.28733e-06;double bestZ2=1.40615e-05;double bestX3=-2.04628e-08; double bestY3=-3.49567e-09;double bestZ3 =2.80403e-08;;double bestW1 =0.045;


  // new pars1:

   //double bestX1=0.0405337; double bestY1 =0.0771617;double bestZ1 =0.114501;double bestX2 =3.70663e-05; double bestY2=-1.8391e-06;double bestZ2=1.42794e-05;double bestX3=-2.17112e-08; double bestY3=-4.80447e-09;double bestZ3 =2.67082e-08;;double bestW1 =0.05;

   // new pars2:

   //  double bestX1=0.0408761; double bestY1 =0.0780871;double bestZ1 =0.115239;double bestX2 =3.77442e-05; double bestY2=-1.39171e-06;double bestZ2=1.44963e-05;double bestX3=-2.30762e-08; double bestY3=-5.83356e-09;double bestZ3 =2.56976e-08;;double bestW1 =0.045;

   // Sept from new pars2, best so far:
   //double bestX1=0.0395046; double bestY1 =0.0781585;double bestZ1 =0.115296;double bestX2 =3.84302e-05; double bestY2=-5.96698e-06;double bestZ2=1.59277e-05;double bestX3=-2.30225e-08; double bestY3=-7.04039e-09;double bestZ3 =2.39387e-08;;double bestW1 =0.04;

   // New cuts 2 (by hand):
   //double bestX1=0.0379046; double bestY1 =0.0752585;double bestZ1 =0.129396;double bestX2 =3.44302e-05; double bestY2=4.03302e-06;double bestZ2=2.39277e-05;double bestX3=-2.30225e-08; double bestY3=-7.04039e-09;double bestZ3 =2.59387e-08;;double bestW1 =0.042;
   //   double bestX1=0.0395046; double bestY1 =0.0781585;double bestZ1 =0.1;double bestX2 =3.84302e-05; double bestY2=-5.96698e-06;double bestZ2=1.59277e-05;double bestX3=-2.30225e-08; double bestY3=-7.04039e-09;double bestZ3 =2.39387e-08;;double bestW1 =0.04;

   double bestX1=0.0301259; double bestY1 =0.072341;double bestZ1 =0.122424;double bestX2 =4.36049e-05; double bestY2=-2.60055e-06;double bestZ2=1.74974e-05;double bestX3=-2.23545e-08; double bestY3=-7.97219e-09;double bestZ3 =1.18452e-08;;double bestW1 =0.04;



   
   // 2 slots:

   //   double bestX1=0.0360361; double bestY1 =0.0820871;double bestZ1 =0.127339;double bestX2 =3.37442e-05; double bestY2=8.60829e-06;double bestZ2=2.44963e-05;double bestX3=-2.30762e-08; double bestY3=-5.83356e-09;double bestZ3 =2.56976e-08;;double bestW1 =0.045;

   // HT 1.76:

   //double bestX1=0.0408765; double bestY1 =0.0774474;double bestZ1 =0.114558;double bestX2 =3.77515e-05; double bestY2=-1.38254e-06;double bestZ2=1.45073e-05;double bestX3=-2.17934e-08; double bestY3=-6.21398e-09;double bestZ3 =2.48111e-08;;double bestW1 =0.045;

   // HT new:

   // double bestX1=0.0404311; double bestY1 =0.0775187;double bestZ1 =0.114615;double bestX2 =3.84332e-05; double bestY2=-9.30463e-07;double bestZ2=1.47298e-05;double bestX3=-2.25027e-08; double bestY3=-8.08667e-09;double bestZ3 =2.2282e-08;;double bestW1 =0.0449809;



   double minimum = myFunc(bestX1,bestY1,bestZ1,bestX2,bestY2,bestZ2,bestX3,bestY3,bestZ3,bestW1);




	 
 for (int i=0; i<1; ++i) 
   {
      //cout<<"--------------------  "<<i<<endl;
		 /*
				 double tryX1 = gRandom->Uniform(0.02,.08);
      double tryY1 = gRandom->Uniform(0.02,0.2);
      double tryZ1 = gRandom->Uniform(0.,0.5);

    	double tryX2 = gRandom->Uniform(-.001,.001);
      double tryY2 = gRandom->Uniform(-.001,.001);
      double tryZ2 = gRandom->Uniform(-.001,.001);

      double tryX3 = gRandom->Uniform(-.00001,.00001);
      double tryY3 = gRandom->Uniform(-.00001,.00001);
      double tryZ3 = gRandom->Uniform(-.00001,.00001);


		 */

		  double tryX1 = bestX1;
      double tryY1 = bestY1;
      double tryZ1 = bestZ1;

    	double tryX2 = bestX2;
      double tryY2 = bestY2;
      double tryZ2 = bestZ2;

      double tryX3 = bestX3;
      double tryY3 = bestY3;
      double tryZ3 = bestZ3;


      double tryW1 = bestW1;
 
      minimizer->SetParameter(0,"X1",tryX1,.0001,0,0);
      minimizer->SetParameter(1,"Y1",tryY1,.0001,0,0);
      minimizer->SetParameter(2,"Z1",tryZ1,.0001,0,0);



      minimizer->SetParameter(3,"X2",tryX2,.00001,0,0);
      minimizer->SetParameter(4,"Y2",tryY2,.00001,0,0);
      minimizer->SetParameter(5,"Z2",tryZ2,.00001,0,0);
   
      minimizer->SetParameter(6,"X3",tryX3,.000001,0,0);  
      minimizer->SetParameter(7,"Y3",tryY3,.000001,0,0);
      minimizer->SetParameter(8,"Z3",tryZ3,.000001,0,0);

      minimizer->SetParameter(9,"W1",tryW1,.01,0,0);
    



      minimizer->ExecuteCommand("SIMPLEX",0,0);
			 //			 	 minimizer->ExecuteCommand("MIGRAD",0,0);
      double newX1 = minimizer->GetParameter(0); 
      double newY1 = minimizer->GetParameter(1); 
      double newZ1 = minimizer->GetParameter(2);

      double newX2 = minimizer->GetParameter(3); 
      double newY2 = minimizer->GetParameter(4); 
      double newZ2 = minimizer->GetParameter(5);

      double newX3 = minimizer->GetParameter(6);
      double newY3 = minimizer->GetParameter(7);
      double newZ3 = minimizer->GetParameter(8);
      
      double newW1 = minimizer->GetParameter(9);
       
      cout <<"i=  "<<i<<endl;

	
			// if(newX1>0&&newY1>0&&newZ1>0&&newX1>0&&newY1>0&&newZ1>0)
      double newMin = myFunc(newX1, newY1, newZ1, newX2, newY2, newZ2, newX3,newY3,newZ3, newW1);
      if (newMin < minimum) 
	{
	  minimum = newMin;
	  bestX1 = newX1; bestY1 = newY1; bestZ1 = newZ1;
	  bestX2 = newX2; bestY2 = newY2; bestZ2 = newZ2;
          bestX3=newX3;   bestY3=newY3; bestZ3=newZ3; 
          bestW1=newW1;
	  cout<<bestX1<<"  "<<bestY1<<"  "<<bestZ1<<"  "<<newMin<<"  "<<i<<endl;
	  bestX2 = newX2; bestY2 = newY2; bestZ2 = newZ2;
	  cout<<bestX2<<"  "<<bestY2<<"  "<<bestZ2<<"  "<<i<<endl;
	  cout<<bestX3<<"  "<<bestY3<<"  "<<bestZ3<<"  "<<i<<endl;
	  cout<<bestW1<<"  "<<i<<endl; 
	}
   }  
   ofstream temp2("temp2.txt");
   temp2<<bestX1<<"  "<<bestX2<<"  "<<bestX3<<"  "<<bestY1<<"  "<<bestY2<<"  "<<bestY3<<"  "<<bestZ1<<" "<<bestZ2<<" " <<bestZ3<<" "<<bestW1<<endl;
   temp2.close();
   return newMin;
}

double myFunc(double x1, double y1, double z1,double x2, double y2, double z2, double x3, double y3, double z3, double w1) 
{


	//tou=10000;
	// if(x1<0||x2<0||y1<0||y2<0||z1<0||z2<0) return tou;
 
	
	
   int conteur =0; 
   float n1,n2,n3,zc,nc,nh,nd,ee;  
   int N_1[2000000]; int N_2[2000000]; int N_3[2000000]; int E[2000000];  int N_C[2000000]; int N_H[2000000]; int N_D[2000000];
   ifstream f; f.open("temporaire.txt");  
   while(f>>n1>>n2>>n3>>nh>>ee){N_1[conteur]=n1; N_2[conteur]=n2; N_3[conteur]=n3;N_C[conteur]=n1+n2+n3+nh; N_H[conteur]=nh; E[conteur]=ee;conteur++;}
		
   float test=0.; 
   float energy=0.;
	
   for(int i=0;i<conteur;i++)
   {

     int mm=E[i]; int N1=N_1[i]; int N2=N_2[i]; int N3=N_3[i]; int NC=N_C[i]+N_H[i]; int NH=N_H[i];
     
     //     if(mm==20)mm=16;
     // if(mm==30)mm=26;
     //     if(mm==80)mm=90;
     
     energy= (x1+x2*NC+x3*NC*NC)*N1+(y1+y2*NC+y3*NC*NC)*N2+(z1+z2*NC+z3*NC*NC)*N3 + w1* NH ;      
     if (mm>0) test=test+((mm-energy)*(mm-energy)/mm);   
   }
	 
	 //	 cout << test/(conteur+1) << "    " << RES123<<endl;
	 return  (test/(conteur+1));
	
      
	 }

void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg) 
{
  result = myFunc(par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9]);
}
//=======================================================================================
//=======================================================================================



/*
cp  /group/ilc/mirabito/SPS2012/showers_715695_0.root 8_GeVL0.root
$ cp /group/ilc/mirabito/SPS2012/showers_715693_0.root 10_GeVL0.root
$ cp /group/ilc/mirabito/SPS2012/showers_715693_3.root 10_GeVL3.root
$ cp /group/ilc/mirabito/SPS2012/showers_715675_0.root 20_GeVL0.root
$ cp /group/ilc/mirabito/SPS2012/showers_715675_3.root 20_GeVL3.root
$ cp /group/ilc/mirabito/SPS2012/showers_715671_0.root 30_GeVL0.root
$ cp /group/ilc/mirabito/SPS2012/showers_715671_3.root 30_GeVL3.root
$ cp /group/ilc/mirabito/SPS2012/showers_715651_0.root 40_GeVL0.root
$ cp /group/ilc/mirabito/SPS2012/showers_715651_3.root 40_GeVL3.root
$ cp /group/ilc/mirabito/SPS2012/showers_715596_0.root 50_GeVL0.root 
$ cp /group/ilc/mirabito/SPS2012/showers_715596_3.root 50_GeVL3.root 
$ cp /group/ilc/mirabito/SPS2012/showers_715531_0.root 60_GeVL0.root 
$ cp /group/ilc/mirabito/SPS2012/showers_715531_3.root 60_GeVL3.root 
$ cp /group/ilc/mirabito/SPS2012/showers_715493_0.root 70_GeVL0.root
$ cp /group/ilc/mirabito/SPS2012/showers_715493_3.root 70_GeVL3.root
$ cp /group/ilc/mirabito/SPS2012/showers_715491_3.root 80_GeVL3.root
$ cp /group/ilc/mirabito/SPS2012/showers_715491_0.root 80_GeVL0.root
$ cp /group/ilc/mirabito/SPS2012/showers_715492_0.root 90_GeVL0.root
*/

/*
common parameters for May and August 
double bestX1=0.0420046; double bestY1 =0.05344;double bestZ1 =0.258335;double bestX2 =2.57049e-05; double bestY2=7.77407e-05;double bestZ2=-7.89366e-05;double bestX3=-6.22608e-09; double bestY3=-9.53908e-08;double bestZ3 =-9.74462e-08;
*/
/*
very good for ugust up to 100
double bestX1=0.0469224; double bestY1 =0.0417256;double bestZ1 =0.252199;double bestX2 =-5.85709e-06; double bestY2=0.000127489;double bestZ2=-4.41272e-05;double bestX3=1.10296e-08; double bestY3=-6.78472e-08;double bestZ3 =-8.52935e-08;

last
double bestX1=0.0312548; double bestY1 =0.0707213;double bestZ1 =0.20901;double bestX2 =3.9031e-05; double bestY2=2.40972e-05;double bestZ2=-4.98226e-06;double bestX3=-1.7052e-08; double bestY3=-1.27052e-08;double bestZ3 =-2.42474e-08;

*/
