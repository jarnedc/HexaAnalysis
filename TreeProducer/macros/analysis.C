#include <iostream>

#include "classes/HQClass.h"

using namespace std;


double covariance(vector<vector<double> >* &cov, int idx, int i, int j){
    int dim = 3;                   //vtx covariance is 3x3
    if(cov->size() == 16) dim = 4; //track covarinance is 4x4
    double covij = cov->at(idx).at(i*dim+j);
    return covij;
}

double vector_size(vector<double> a){
 return pow(a[0]*a[0]+a[1]*a[1]+a[2]*a[2],0.5);
}

double vector_product_size(vector<double> a,vector<double> b){
  return pow(a[0]*b[0]+a[1]*b[1]+a[2]*b[2],0.5);

}

double vector_sum_size(vector<double> a,vector<double> b){
  vector<double> temp;
  temp.push_back(a[0]+b[0]);
  temp.push_back(a[1]+b[1]);
  temp.push_back(a[2]+b[2]);
  return vector_size(temp);
}

double vector_sum_size3(vector<double> a,vector<double> b, vector<double> c){
  vector<double> temp;
  temp.push_back(a[0]+b[0]+c[0]);
  temp.push_back(a[1]+b[1]+c[1]);
  temp.push_back(a[2]+b[2]+c[2]);
  return vector_size(temp);
}

double Energy(vector<double> impuls, double mass){
  return pow(vector_size(impuls)*vector_size(impuls)+mass*mass,0.5);
}

void plotpdgids(){
    TChain *tree = new TChain("tree/HexaQAnalysis");
    tree->Add("../test/LambdaGun_GENSIM_IIDD_3122_05to10GeV_tree.root");
    string analysisPlots  = "pdgids.root";
    TFile *fAnalysisPlots =new TFile(analysisPlots.c_str(),"RECREATE");
    fAnalysisPlots->cd();
    TH1F *h_pdgids = new TH1F("pdgids", "pdgids lambda decay (Gev)", 10000,-5000, 5000);
    
    HQClass hqhand;
    hqhand.Init(tree);

    Long64_t nEntries = tree->GetEntries();   
    
    for (int i=0; i<nEntries; i++){
   	hqhand.GetEntry(i);
        for(int genp = 0; genp < hqhand.gen_pdgid->size(); genp++){
             h_pdgids->Fill(hqhand.gen_pdgid->at(genp));
        }

    }
    h_pdgids->Write();
}

void LambdaInvMass(){
    TChain *tree = new TChain("tree/HexaQAnalysis");
    tree->Add("../test/LambdaGun_GENSIM_IIDD_3122_05to10GeV_tree.root");
   // tree->Add("/user/jdeclerc/SSearch/HexaAnalysis/CMSSW_8_0_21/src/HexaAnalysis/TreeProducer/test/MC_tree_Lambda.root");
 

    string analysisPlots  = "analysisPlots_Lambda_sample.root";
    TFile *fAnalysisPlots =new TFile(analysisPlots.c_str(),"RECREATE");
    fAnalysisPlots->cd();
    TH1F *h_lambda_mass = new TH1F("lambda_mass", "mass Lambda (Gev)", 100, 1.07, 1.16);

    HQClass hqhand;
    hqhand.Init(tree);

    Long64_t nEntries = tree->GetEntries();
    cout<<nEntries<<endl;

    
   
    for (int i=0; i<nEntries; i++){
    cout << "----------------------------------" << endl;
        hqhand.GetEntry(i);
        vector<Double_t> p_proton;
	vector<Double_t> p_pion;
	vector<Double_t> p_lambda;
	cout << "nEvent: " << hqhand.nEvent << endl;
	int charge = 0;
	for(int gent=0; gent<hqhand.nTrack; gent++){
	    charge = hqhand.track_charge->at(gent);
	    cout << "charge: " << charge << endl;
	    cout << "pt: " << hqhand.track_pt->at(gent) << endl;
//	    cout << "fromPV" << hqhand.track_fromPV->at(genp) << endl;
	    if(charge == 1){ //proton
		p_proton.push_back(hqhand.track_px->at(gent));
		p_proton.push_back(hqhand.track_py->at(gent));
		p_proton.push_back(hqhand.track_pz->at(gent));	
		cout << "proton found" << endl;
		cout << "++++++++++++" << endl; 
	    }
	    else if(charge == -1){ //pi minus
		p_pion.push_back(hqhand.track_px->at(gent));
		p_pion.push_back(hqhand.track_py->at(gent));
		p_pion.push_back(hqhand.track_pz->at(gent));	
		cout << "pi found" << endl;
		cout << "++++++++++++" << endl;
	    }

	    else continue;
	}
	if(p_proton.size() != 3 || p_pion.size() != 3) continue;
	Double_t energy_proton = Energy(p_proton,m_proton);
	Double_t energy_pi_charged = Energy(p_pion, m_pi_charged);
	Double_t p_proton_pi_charged = vector_product_size(p_proton,p_pion);
	Double_t m_inv = pow(pow(energy_proton+energy_pi_charged,2)-pow(vector_sum_size(p_proton,p_pion),2),0.5);
        cout << m_inv << endl;
	h_lambda_mass->Fill(m_inv);


    }
    h_lambda_mass->Write();
    fAnalysisPlots->Write();   
}

void LambdaXiInvMass(){
    TChain *tree = new TChain("tree/HexaQAnalysis");
    //tree->Add("../test/MC_tree_Lambda2.root");
    tree->Add("../test/AntiXiGun_GENSIM_IIDD_min3312_05to10GeV_tree.root");
 
    string analysisPlots  = "analysisPlots_Xi_sample.root";
    TFile *fAnalysisPlots =new TFile(analysisPlots.c_str(),"RECREATE");
    fAnalysisPlots->cd();
    TH1F *h_lambda_mass = new TH1F("lambda_mass", "Lambda mass (Gev)", 100, 0., 2.);
    TH1F *h_Xi_mass_lambda_m_const = new TH1F("h_Xi_mass_lambda_m_const", "Xi mass (no lambda mass constraint) (Gev)", 100, 0., 2.);
    TH1F *h_Xi_mass_no_lambda_m_const = new TH1F("h_Xi_mass_no_lambda_m_const", "mass (lambda mass constraint) (Gev)", 100, 0., 2.);
    TH1F *h_nrPosParticles = new TH1F("nrPosParticles", "nrPosParticles", 10, 0, 10);
    HQClass hqhand;
    hqhand.Init(tree);

    Long64_t nEntries = tree->GetEntries();
    cout<<nEntries<<endl;

    
   
    for (int i=0; i<nEntries; i++){
    cout << "----------------------------------" << endl;
        hqhand.GetEntry(i);
        vector<Double_t> p_proton;
	vector<vector<Double_t>> p_positive; //save the momenta of the positive particles. These can be pions from the Xi or Lambda decay or the Xi itself
	vector<Double_t> p_pion_lambda;
	vector<Double_t> p_pion_xi;
	vector<Double_t> p_lambda;
	vector<Double_t> p_Xi;
	cout << "nEvent: " << hqhand.nEvent << endl;
	int charge = 0;
	int nrPosParticles = 0;
	for(int genp=0; genp<hqhand.nTrack; genp++){
	    charge = hqhand.track_charge->at(genp);
	    cout << "charge: " << charge << endl;
	    cout << "pt: " << hqhand.track_pt->at(genp) << endl;
	    if(charge == -1){ //anti-proton
		p_proton.push_back(hqhand.track_px->at(genp));
		p_proton.push_back(hqhand.track_py->at(genp));
		p_proton.push_back(hqhand.track_pz->at(genp));	
		cout << "anti-proton found" << endl;
		cout << "++++++++++++" << endl; 
	    }
	    else if(charge == 1){ //positive particle
		p_positive.push_back(vector<double>());
	  	p_positive[nrPosParticles].push_back(hqhand.track_px->at(genp));
		p_positive[nrPosParticles].push_back(hqhand.track_py->at(genp));
		p_positive[nrPosParticles].push_back(hqhand.track_pz->at(genp));
		cout << "positive particle found" << endl;
		nrPosParticles++;
		cout << "++++++++++++" << endl;
	
	    }
	    else continue;
	}
	h_nrPosParticles->Fill(nrPosParticles);	
	if(p_proton.size() != 3 || nrPosParticles != 2 ) continue;
	if(vector_size(p_positive[0]) > vector_size(p_positive[1])){
		 cout << "1" << endl;
		 p_pion_xi.push_back(p_positive[0][0]);
		 p_pion_xi.push_back(p_positive[0][1]);
		 p_pion_xi.push_back(p_positive[0][2]);
		 p_pion_lambda.push_back(p_positive[1][0]);
                 p_pion_lambda.push_back(p_positive[1][1]);
                 p_pion_lambda.push_back(p_positive[1][2]);
	}
	else{
		 cout << "2" << endl;
		 p_pion_xi.push_back(p_positive[1][0]);
                 p_pion_xi.push_back(p_positive[1][1]);
                 p_pion_xi.push_back(p_positive[1][2]);
                 p_pion_lambda.push_back(p_positive[0][0]);
                 p_pion_lambda.push_back(p_positive[0][1]);
                 p_pion_lambda.push_back(p_positive[0][2]);	
	}
	Double_t energy_proton = Energy(p_proton,m_proton);
	Double_t energy_pion_xi = Energy(p_pion_xi, m_pi_charged);
	Double_t energy_pion_lambda = Energy(p_pion_lambda, m_pi_charged);
	Double_t p_proton_pi_charged = vector_product_size(p_proton,p_pion_lambda);
	Double_t m_inv_lambda = pow(pow(energy_proton+energy_pion_lambda,2)-pow(vector_sum_size(p_proton,p_pion_lambda),2),0.5);
        cout << "m_inv_lambda: " <<  m_inv_lambda << endl;
	h_lambda_mass->Fill(m_inv_lambda);
	p_lambda.push_back(p_proton[0]+p_pion_lambda[0]);
	p_lambda.push_back(p_proton[1]+p_pion_lambda[1]);
	p_lambda.push_back(p_proton[2]+p_pion_lambda[2]);
	Double_t energy_lambda = Energy(p_lambda,m_lambda);	
	Double_t m_inv_Xi_Lambda_mass_const =  pow(pow(energy_lambda+energy_pion_xi,2)-pow(vector_sum_size(p_lambda,p_pion_xi),2),0.5);
	Double_t m_inv_Xi_no_Lambda_mass_const = pow(pow(energy_proton+energy_pion_xi+energy_pion_lambda,2)-pow(vector_sum_size3(p_proton,p_pion_xi,p_pion_lambda),2),0.5);
	cout << "m_inv_Xi: " << m_inv_Xi_Lambda_mass_const << endl;
	cout << "m_inv_Xi2: " << m_inv_Xi_no_Lambda_mass_const << endl;
	h_Xi_mass_lambda_m_const->Fill(m_inv_Xi_Lambda_mass_const);
	h_Xi_mass_no_lambda_m_const->Fill(m_inv_Xi_no_Lambda_mass_const);
    }
    h_lambda_mass->Write();
    h_Xi_mass_lambda_m_const->Write();
    h_Xi_mass_no_lambda_m_const->Write();
    h_nrPosParticles->Write();
    fAnalysisPlots->Write();   
}

int analysis(){
   plotpdgids();
   //calculate the LambdaInvMass
   LambdaInvMass();
   LambdaXiInvMass();
   return 0;

}
