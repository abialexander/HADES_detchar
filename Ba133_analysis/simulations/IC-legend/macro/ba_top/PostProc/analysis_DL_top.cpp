//{
//gInterpreter->AddIncludePath("/lfs/l1/legend/users/bianca/IC_geometry/IC-gerda/general_geometry/macro/PostProcessing/");

#include "T4SimGeometry.h"
//}


void analysis_DL_top(double dlt,  double fccd, string output){
	string name="IC160A";
        double radius=39.9;
        double height=65.4;
        double grooveOuterRadius=15.5;
        double grooveInnerRadius=11.3;
        double grooveDepth=2.0;
        double coneRadius=37.75;
        double coneHeight=20.1;
        double boreRadius=4.65;
        double boreDepth=33.7;
        double a_smear=0.35;
        double b_smear=1.99e-3;

        TF1 *f_smear = new TF1("f_smear","sqrt([0]+[1]*x)",0,120);
        f_smear->SetParameter(0,a_smear);
        f_smear->SetParameter(1,b_smear);

	
	T4SimHPGe *detector= new T4SimHPGe(name,
   	fccd, dlt, radius,  height,
        grooveOuterRadius , grooveInnerRadius,
        grooveDepth,
        coneRadius ,
        coneHeight ,
        boreRadius ,
 	boreDepth );

	cout << "start" << endl;

	TFile *f = TFile::Open("top_4z_1x/energy_spectrum_IC160A_top_HS_run0003_4z_1x.root");
	//TFile *f = TFile::Open("energy_spectrum_IC160A_top_HS_run0003.root");
   	TTree *t = (TTree*)f->Get("new_tree");

        vector<double> *energy_hit = new vector<double> (200);
        vector<double> *x = new vector<double> (200);
        vector<double> *y = new vector<double> (200);
        vector<double> *z = new vector<double> (200);
        
        vector<double> *new_energy_hit = new vector<double> (200);
	vector<double> *r = new vector<double> (200);
	vector<double> *new_z = new vector<double> (200);

        t->SetBranchAddress("Edep", &energy_hit);
        t->SetBranchAddress("x", &x);
        t->SetBranchAddress("y", &y);
        t->SetBranchAddress("z", &z);

        int entries= t->GetEntries();

	TFile* output_file = new TFile(("top_4z_1x/energy_spectrum_am_run0003_4z_1x_"+output+".root").c_str(),"recreate");
	//TFile* output_file = new TFile(("energy_spectrum_am_run0003_"+output+".root").c_str(),"recreate");
	TH1D* dlt_histo_energy= new TH1D ("dlt_histo_energy", "dlt_histo_energy", 24000, 0, 120);
	dlt_histo_energy->SetTitle(" ; energy [keV]; counts");

	TTree *tree_dlt = new TTree ("tree_dlt", "tree_dlt");
	tree_dlt->Branch("new_energy_hit", &new_energy_hit);
	tree_dlt->Branch("r",&r);
	tree_dlt->Branch("new_z",&new_z);
	
	//FILE *ff = fopen("energy_spectrum_run0003_02dlt.txt","w");

	double efficiency;
	int n=0;
	int p=0;
	int l=0;
	int m=0;
        for (int i=0; i<entries; i++){
		t->GetEntry(i);
		int hits=energy_hit->size();
		double total_energy=0;
	      	double smear_total_energy=0;
		for (int j=0; j<hits; j++){
			r->at(j) = sqrt(pow(x->at(j),2)+pow(y->at(j),2));
			new_z->at(j)=(z->at(j)-7.0);
			//if ((z->at(j))<7)  cout << z->at(j) << endl;
			
			efficiency =detector->GetChargeCollectionEfficiency(r->at(j), new_z->at(j));
			new_energy_hit->at(j)=energy_hit->at(j)*efficiency;
			total_energy+=new_energy_hit->at(j);
			
			double fHeight=detector->GetHeight();
			double fRadius=detector->GetRadius();
			//if ((new_z->at(j))<0)  cout << new_z->at(j) << endl;
			if ((r->at(j))>fRadius)  cout << r->at(j) << endl;
			
			
			if((r->at(j))>fRadius || (new_z->at(j))>fHeight || (new_z->at(j))<0)  p++;		
			if (efficiency==0)  n++;
			if (efficiency==1)  l++;
		
			m++;
     	 		//tree_dlt->Fill();
		}
		if (total_energy==0) continue;
			 smear_total_energy=gRandom->Gaus(total_energy*1000,(f_smear->Eval(total_energy*1000))/2.355);
			 dlt_histo_energy->Fill(smear_total_energy);
			 //fprintf(ff,"%f\n",smear_total_energy);
		
	 }
	
	dlt_histo_energy->Write();
	//tree_dlt->Write();


	cout << "n " << n << endl;
	cout << "l " << l << endl;
	cout << "p " << p << endl;
	cout << "n+l " << n+l << endl;
	cout << "m " << m << endl;
		




}



