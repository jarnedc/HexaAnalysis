//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 18 14:25:27 2017 by ROOT version 6.06/01
// from TTree HexaQAnalysis/tree
// found on file: MC_tree.root
//////////////////////////////////////////////////////////

#ifndef HQClass_h
#define HQClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

//Header file for the classes stored in the TTree if any.
#include <vector>
#include <TMath.h>
#include <TH2F.h>


//constants
const Double_t m_proton = 0.9382720813;
const Double_t m_pi_charged = 0.13957018;
const Double_t m_lambda = 1.115683;

class HQClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           nEvent;
   Int_t           nRun;
   Int_t           nLumi;
   Int_t           vtx_N;
   Int_t           vtx_N_stored;
   vector<double>  *vtx_normalizedChi2;
   vector<int>     *vtx_ndof;
   vector<int>     *vtx_nTracks;
   vector<double>  *vtx_d0;
   vector<double>  *vtx_x;
   vector<double>  *vtx_y;
   vector<double>  *vtx_z;
   vector<vector<double> > *vtx_covariance;
   vector<bool>    *vtx_isFake;
   vector<bool>    *vtx_isValid;
   Int_t           nTrack_stored;
   Int_t           nTrack;
   vector<double>  *track_pt;
   vector<double>  *track_px;
   vector<double>  *track_py;
   vector<double>  *track_pz;
   vector<double>  *track_eta;
   vector<double>  *track_phi;
   vector<double>  *track_normalizedChi2;
   vector<int>     *track_ndof;
   vector<double>  *track_ptError;
   vector<double>  *track_dzError;
   vector<double>  *track_dz;
   vector<int>     *track_fromPV;
   vector<int>     *track_purity;
   vector<int>     *track_nhits;
   vector<int>     *track_nPixHits;
   vector<double>  *track_d0;
   vector<int>     *track_charge;
   vector<double>  *track_dxy;
   vector<vector<double> > *track_covariance;
   vector<double>  *gen_px;
   vector<double>  *gen_py;
   vector<double>  *gen_pz;
   vector<double>  *gen_p;
   vector<double>  *gen_pt;
   vector<double>  *gen_eta;
   vector<double>  *gen_phi;
   vector<double>  *gen_mass;
   vector<double>  *gen_energy;
   vector<int>     *gen_charge;
   vector<int>     *gen_pdgid;
   vector<int>     *gen_status;
   
   Int_t           HLT_PFJet450;
   Double_t        pswgt_singlejet_450;

   // List of branches
   TBranch        *b_nEvent;   //!
   TBranch        *b_nRun;   //!
   TBranch        *b_nLumi;   //!
   TBranch        *b_vtx_N;   //!
   TBranch        *b_vtx_N_stored;   //!
   TBranch        *b_vtx_normalizedChi2;   //!
   TBranch        *b_vtx_ndof;   //!
   TBranch        *b_vtx_nTracks;   //!
   TBranch        *b_vtx_d0;   //!
   TBranch        *b_vtx_x;   //!
   TBranch        *b_vtx_y;   //!
   TBranch        *b_vtx_z;   //!
   TBranch        *b_vtx_covariance;   //!
   TBranch        *b_vtx_isFake;   //!
   TBranch        *b_vtx_isValid;   //!
   TBranch        *b_nTrack_stored;   //!
   TBranch        *b_nTrack;   //!
   TBranch        *b_track_pt;   //!
   TBranch        *b_track_px;   //!
   TBranch        *b_track_py;   //!
   TBranch        *b_track_pz;   //!
   TBranch        *b_track_eta;   //!
   TBranch        *b_track_phi;   //!
   TBranch        *b_track_normalizedChi2;   //!
   TBranch        *b_track_ndof;   //!
   TBranch        *b_track_ptError;   //!
   TBranch        *b_track_dzError;   //!
   TBranch        *b_track_dz;   //!
   TBranch        *b_track_fromPV;   //!
   TBranch        *b_track_purity;   //!
   TBranch        *b_track_nhits;   //!
   TBranch        *b_track_nPixHits;   //!
   TBranch        *b_track_d0;   //!
   TBranch        *b_track_charge;   //!
   TBranch        *b_track_dxy;   //!
   TBranch        *b_track_covariance;   //!
   TBranch        *b_gen_px;   //!
   TBranch        *b_gen_py;   //!
   TBranch        *b_gen_pz;   //!
   TBranch        *b_gen_p;   //!
   TBranch        *b_gen_pt;   //!
   TBranch        *b_gen_eta;   //!
   TBranch        *b_gen_phi;   //!
   TBranch        *b_gen_mass;   //!
   TBranch        *b_gen_energy;   //!
   TBranch        *b_gen_charge;   //!
   TBranch        *b_gen_pdgid;   //!
   TBranch        *b_gen_status;   //!   
   TBranch        *b_HLT_PFJet450;   //!
   TBranch        *b_pswgt_singlejet_450;   //!

   HQClass(TTree *tree=0);
   virtual ~HQClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};


HQClass::HQClass(TTree *tree) : fChain(0) 
{
   Init(tree);
}

HQClass::~HQClass()
{
   if (!fChain) return;
   //delete fChain->GetCurrentFile();
}

Int_t HQClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HQClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void HQClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   vtx_normalizedChi2 = 0;
   vtx_ndof = 0;
   vtx_nTracks = 0;
   vtx_d0 = 0;
   vtx_x = 0;
   vtx_y = 0;
   vtx_z = 0;
   vtx_covariance = 0;
   vtx_isFake = 0;
   vtx_isValid = 0;
   track_pt = 0;
   track_px = 0;
   track_py = 0;
   track_pz = 0;
   track_eta = 0;
   track_phi = 0;
   track_normalizedChi2 = 0;
   track_ndof = 0;
   track_ptError = 0;
   track_dzError = 0;
   track_dz = 0;
   track_fromPV = 0;
   track_purity = 0;
   track_nhits = 0;
   track_nPixHits = 0;
   track_d0 = 0;
   track_charge = 0;
   track_dxy = 0;
   track_covariance = 0;
   gen_px = 0;
   gen_py = 0;
   gen_pz = 0;
   gen_p = 0;
   gen_pt = 0;
   gen_eta = 0;
   gen_phi = 0;
   gen_mass = 0;
   gen_energy = 0;
   gen_charge = 0;
   gen_pdgid = 0;
   gen_status = 0;   
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nEvent", &nEvent, &b_nEvent);
   fChain->SetBranchAddress("nRun", &nRun, &b_nRun);
   fChain->SetBranchAddress("nLumi", &nLumi, &b_nLumi);
   fChain->SetBranchAddress("vtx_N", &vtx_N, &b_vtx_N);
   fChain->SetBranchAddress("vtx_N_stored", &vtx_N_stored, &b_vtx_N_stored);
   fChain->SetBranchAddress("vtx_normalizedChi2", &vtx_normalizedChi2, &b_vtx_normalizedChi2);
   fChain->SetBranchAddress("vtx_ndof", &vtx_ndof, &b_vtx_ndof);
   fChain->SetBranchAddress("vtx_nTracks", &vtx_nTracks, &b_vtx_nTracks);
   fChain->SetBranchAddress("vtx_d0", &vtx_d0, &b_vtx_d0);
   fChain->SetBranchAddress("vtx_x", &vtx_x, &b_vtx_x);
   fChain->SetBranchAddress("vtx_y", &vtx_y, &b_vtx_y);
   fChain->SetBranchAddress("vtx_z", &vtx_z, &b_vtx_z);
   fChain->SetBranchAddress("vtx_covariance", &vtx_covariance, &b_vtx_covariance);
   fChain->SetBranchAddress("vtx_isFake", &vtx_isFake, &b_vtx_isFake);
   fChain->SetBranchAddress("vtx_isValid", &vtx_isValid, &b_vtx_isValid);
   fChain->SetBranchAddress("nTrack_stored", &nTrack_stored, &b_nTrack_stored);
   fChain->SetBranchAddress("nTrack", &nTrack, &b_nTrack);
   fChain->SetBranchAddress("track_pt", &track_pt, &b_track_pt);
   fChain->SetBranchAddress("track_px", &track_px, &b_track_px);
   fChain->SetBranchAddress("track_py", &track_py, &b_track_py);
   fChain->SetBranchAddress("track_pz", &track_pz, &b_track_pz);
   fChain->SetBranchAddress("track_eta", &track_eta, &b_track_eta);
   fChain->SetBranchAddress("track_phi", &track_phi, &b_track_phi);
   fChain->SetBranchAddress("track_normalizedChi2", &track_normalizedChi2, &b_track_normalizedChi2);
   fChain->SetBranchAddress("track_ndof", &track_ndof, &b_track_ndof);
   fChain->SetBranchAddress("track_ptError", &track_ptError, &b_track_ptError);
   fChain->SetBranchAddress("track_dzError", &track_dzError, &b_track_dzError);
   fChain->SetBranchAddress("track_dz", &track_dz, &b_track_dz);
   fChain->SetBranchAddress("track_fromPV", &track_fromPV, &b_track_fromPV);
   fChain->SetBranchAddress("track_purity", &track_purity, &b_track_purity);
   fChain->SetBranchAddress("track_nhits", &track_nhits, &b_track_nhits);
   fChain->SetBranchAddress("track_nPixHits", &track_nPixHits, &b_track_nPixHits);
   fChain->SetBranchAddress("track_d0", &track_d0, &b_track_d0);
   fChain->SetBranchAddress("track_charge", &track_charge, &b_track_charge);
   fChain->SetBranchAddress("track_dxy", &track_dxy, &b_track_dxy);
   fChain->SetBranchAddress("track_covariance", &track_covariance, &b_track_covariance);
   fChain->SetBranchAddress("gen_px", &gen_px, &b_gen_px);
   fChain->SetBranchAddress("gen_py", &gen_py, &b_gen_py);
   fChain->SetBranchAddress("gen_pz", &gen_pz, &b_gen_pz);
   fChain->SetBranchAddress("gen_p", &gen_p, &b_gen_p);
   fChain->SetBranchAddress("gen_pt", &gen_pt, &b_gen_pt);
   fChain->SetBranchAddress("gen_eta", &gen_eta, &b_gen_eta);
   fChain->SetBranchAddress("gen_phi", &gen_phi, &b_gen_phi);
   fChain->SetBranchAddress("gen_mass", &gen_mass, &b_gen_mass);
   fChain->SetBranchAddress("gen_energy", &gen_energy, &b_gen_energy);
   fChain->SetBranchAddress("gen_charge", &gen_charge, &b_gen_charge);
   fChain->SetBranchAddress("gen_pdgid", &gen_pdgid, &b_gen_pdgid);
   fChain->SetBranchAddress("gen_status", &gen_status, &b_gen_status);   
   fChain->SetBranchAddress("HLT_PFJet450", &HLT_PFJet450, &b_HLT_PFJet450);
   fChain->SetBranchAddress("pswgt_singlejet_450", &pswgt_singlejet_450, &b_pswgt_singlejet_450);
   Notify();
}

Bool_t HQClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HQClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HQClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef HQClass_cxx
