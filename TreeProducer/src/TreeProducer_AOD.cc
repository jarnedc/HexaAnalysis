#include "HexaAnalysis/TreeProducer/interface/TreeProducer_AOD.h"

//
// constructors and destructor
//
TreeProducer_AOD::TreeProducer_AOD(const edm::ParameterSet& pset):
  _trigResultsTag(pset.getParameter<edm::InputTag>("triggerResults")),
  _vertexCollectionTag(pset.getParameter<edm::InputTag>("vertexCollection")),
  _trackCollectionTag(pset.getParameter<edm::InputTag>("trackCollection")),
  _trigResultsToken(consumes<edm::TriggerResults>(_trigResultsTag)),
  _vertexCollectionToken(consumes<vector<reco::Vertex> >(_vertexCollectionTag)),
  _trackCollectionToken(consumes<vector<reco::Track> >(_trackCollectionTag)),
  _isData(pset.getUntrackedParameter<bool>("isData")),
  m_partons(consumes<vector<reco::GenParticle> >(pset.getParameter<edm::InputTag>("genCollection"))),
  hltPrescale_(pset, consumesCollector(), *this)
{
 triggerNames_ = pset.getParameter<std::vector<std::string> > ("triggerName");
 
}


TreeProducer_AOD::~TreeProducer_AOD()
{
}

//
// member functions
//

// ------------ method called for each event  ------------
void
TreeProducer_AOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // Initialize branches
  Init();

  // HANDLES //
  // Get collections
  edm::Handle<edm::TriggerResults> H_trig;//, H_trig1, H_trig2;
  iEvent.getByToken(_trigResultsToken, H_trig);

  edm::Handle<vector<reco::Vertex> > H_vert;
  iEvent.getByToken(_vertexCollectionToken, H_vert);

  edm::Handle<vector<reco::Track> > H_track;
  iEvent.getByToken(_trackCollectionToken , H_track);

  edm::Handle<vector<reco::GenParticle> > H_partons;
  if(!_isData) iEvent.getByToken(m_partons, H_partons);


  // Check validity
  if(!H_trig.isValid()) {
    if(verbose>0) cout << "Missing collection : " << _trigResultsTag << " ... skip entry !" << endl;
    return;
  }

  if(!H_vert.isValid()) {
    if(verbose>0) cout << "Missing collection : " << _vertexCollectionTag << " ... skip entry !" << endl;
    return;
  }

  if(!H_track.isValid()) {
    if(verbose>0) cout << "Missing collection : " << _trackCollectionTag << " ... skip entry !" << endl;
    return;
  }

  if(!_isData){
   if(!H_partons.isValid()) {
     if(verbose>0) cout << "Missing collection : genParticles ... skip entry !" << endl;
     return;
   }
  }
  // GLOBAL EVENT INFORMATIONS //
  _nRun   = iEvent.id().run();
  _nLumi  = iEvent.luminosityBlock();
  _nEvent = iEvent.id().event();

  // VERTICES //

  // select the primary vertex as the one with higest sum of (pt)^2 of tracks
  PrimaryVertexSorter PVSorter;
  std::vector<reco::Vertex> sortedVertices = PVSorter.sortedList( *(H_vert.product()) );

  for( std::vector<reco::Vertex>::const_iterator PV = sortedVertices.begin(); PV != sortedVertices.end(); ++PV){
    _vtx_normalizedChi2.push_back(PV->normalizedChi2());
    _vtx_ndof.push_back(PV->ndof());
    _vtx_nTracks.push_back(PV->nTracks());
    _vtx_tracksSize.push_back(PV->tracksSize());
    _vtx_d0.push_back(PV->position().Rho());
    _vtx_x.push_back(PV->x());
    _vtx_y.push_back(PV->y());
    _vtx_z.push_back(PV->z());
    _vtx_isValid.push_back(PV->isValid());
    _vtx_isFake.push_back(PV->isFake());

    std::vector<double> cov;
    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
        cov.push_back(PV->covariance(i,j));
      }
    }
    _vtx_covariance.push_back(cov);

  } // for loop on primary vertices
  _vtx_N = H_vert->size();


  // TRACKS //
    vector<reco::TrackRef> trackRef;
    for (vector<reco::Track>::const_iterator theTrack = H_track->begin(); theTrack != H_track->end(); ++theTrack){
            reco::TrackRef ref(H_track, theTrack - H_track->begin());
                    trackRef.push_back(ref);
    }

    reco::RecoPtSorter<reco::TrackRef> trackSorter;
    std::sort( trackRef.begin(), trackRef.end(), trackSorter);

    ///FIXME: High pt cuts only for testing crab - Avoid large size files
   // if(trackRef.size()==0) return;
   // if(trackRef[0]->pt()<5.) return;

    for (size_t i = 0; i < trackRef.size(); i++) {
           // if(trackRef[1]->pt()<0.5) return;

            _track_purity.push_back(trackRef[i]->highPurity);
            _track_Nhits.push_back(trackRef[i]->numberOfValidHits());
            _track_NpixHits.push_back(trackRef[i]->hitPattern().numberOfValidPixelHits());
            if (trackRef[i]->vz() == 0) _track_fromPV.push_back(1);
            else _track_fromPV.push_back(0);
            _track_pt.push_back(trackRef[i]->pt());
            _track_px.push_back(trackRef[i]->px());
            _track_py.push_back(trackRef[i]->py());
            _track_pz.push_back(trackRef[i]->pz());
            _track_eta.push_back(trackRef[i]->eta());
            _track_phi.push_back(trackRef[i]->phi());
            _track_normalizedChi2.push_back(trackRef[i]->normalizedChi2());
            _track_ndof.push_back(trackRef[i]->ndof());
            _track_ptError.push_back(trackRef[i]->ptError());
            _track_dzError.push_back(trackRef[i]->dzError());
            _track_dz.push_back(trackRef[i]->dz());
            _track_dxy.push_back(trackRef[i]->dxy());
            _track_d0.push_back(trackRef[i]->d0());
            _track_charge.push_back(trackRef[i]->charge());

            std::vector<double> cov;
            for(int k=0; k<4; k++){
              for(int j=0; j<4; j++){
                cov.push_back(trackRef[i]->covariance(k,j));
              }
            }
            _track_covariance.push_back(cov);

    }
    _nTrack = trackRef.size();


    if(!_isData){

        for (vector<reco::GenParticle>::const_iterator thepartons = H_partons->begin();
        thepartons != H_partons->end(); ++thepartons){

        //   if(thepartons->fromHardProcessFinalState()){
            _genp_px.push_back(thepartons->px());
            _genp_py.push_back(thepartons->py());
            _genp_pz.push_back(thepartons->pz());
            _genp_p.push_back(thepartons->p());
            _genp_pt.push_back(thepartons->pt());
            _genp_eta.push_back(thepartons->eta());
            _genp_phi.push_back(thepartons->phi());
            _genp_mass.push_back(thepartons->mass());
            _genp_energy.push_back(thepartons->energy());
            _genp_charge.push_back(thepartons->charge());
            _genp_pdgid.push_back(thepartons->pdgId());
            _genp_status.push_back(thepartons->status());
            
          // }//if HardProcess
        }///for partons
    }

//TRIGGER//
if (triggerPathsMap[triggerPathsVector[0]] != -1 && H_trig->accept(triggerPathsMap[triggerPathsVector[0]])) _singlejet_450 = 1;
// sanity check
assert(H_trig->size() == hltConfig_.size());
//------ loop over all trigger names ---------
for(unsigned itrig=0;itrig<triggerNames_.size() ;itrig++) {
  int preL1(-1);
  int preHLT(-1);

  if (triggerIndex_[itrig] < hltConfig_.size()) {
    ///In detail get prescale info from hltConfig_
    std::pair<std::vector<std::pair<std::string,int> >,int> detailedPrescaleInfo = hltPrescale_.prescaleValuesInDetail(iEvent, iSetup, triggerNames_[itrig]);
    preHLT = detailedPrescaleInfo.second ;

    // save l1 prescale values in standalone vector
    std::vector <int> l1prescalevals;
    for( size_t varind = 0; varind < detailedPrescaleInfo.first.size(); varind++ ){
      l1prescalevals.push_back(detailedPrescaleInfo.first.at(varind).second);
    }

    //find and save minimum l1 prescale of any ORed L1 that seeds the HLT
    std::vector<int>::iterator result = std::min_element(std::begin(l1prescalevals), std::end(l1prescalevals));
    size_t minind = std::distance(std::begin(l1prescalevals), result);
    // sometimes there are no L1s associated with a HLT. In that case, this branch stores -1 for the l1prescale
    preL1 = minind < l1prescalevals.size() ? l1prescalevals.at(minind) : -1 ;//commented for 76X

    //prescales
    if(triggerNames_[itrig].find("HLT_PFJet450_v") != std::string::npos)   _pswgt_singlejet_450 = preL1 * preHLT;

   }
  }

  _tree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void
TreeProducer_AOD::beginJob()
{
        // Initialize when class is created
        edm::Service<TFileService> fs ;
        _tree = fs->make <TTree>("HexaQAnalysis","tree");

        // Declare tree's branches
        // Event
        _tree->Branch("nEvent",&_nEvent,"nEvent/I");
        _tree->Branch("nRun",&_nRun,"nRun/I");
        _tree->Branch("nLumi",&_nLumi,"nLumi/I");
        //
        // Vertices
        _tree->Branch("vtx_N",&_vtx_N,"vtx_N/I");
        _tree->Branch("vtx_N_stored",&_vtx_N_stored,"vtx_N_stored/I");
        _tree->Branch("vtx_normalizedChi2",&_vtx_normalizedChi2);
        _tree->Branch("vtx_ndof",&_vtx_ndof);
        _tree->Branch("vtx_nTracks",&_vtx_nTracks);
        _tree->Branch("vtx_d0",&_vtx_d0);
        _tree->Branch("vtx_x",&_vtx_x);
        _tree->Branch("vtx_y",&_vtx_y);
        _tree->Branch("vtx_z",&_vtx_z);
        _tree->Branch("vtx_covariance",&_vtx_covariance);
        _tree->Branch("vtx_isFake",&_vtx_isFake);
        _tree->Branch("vtx_isValid",&_vtx_isValid);
        //
	// Tracks
        _tree->Branch("nTrack_stored",&_nTrack_stored,"nTrack_stored/I");
        _tree->Branch("nTrack",&_nTrack,"nTrack/I");
        _tree->Branch("track_pt",&_track_pt);
        _tree->Branch("track_px",&_track_px);
        _tree->Branch("track_py",&_track_py);
        _tree->Branch("track_pz",&_track_pz);
        _tree->Branch("track_eta",&_track_eta);
	_tree->Branch("track_phi",&_track_phi);
	_tree->Branch("track_normalizedChi2",&_track_normalizedChi2);
	_tree->Branch("track_ndof",&_track_ndof);
	_tree->Branch("track_ptError",&_track_ptError);
	_tree->Branch("track_dzError",&_track_dzError);
	_tree->Branch("track_dz",&_track_dz);
	_tree->Branch("track_fromPV",&_track_fromPV);
	_tree->Branch("track_purity",&_track_purity);
	_tree->Branch("track_nhits",&_track_Nhits);
	_tree->Branch("track_nPixHits",&_track_NpixHits);
        _tree->Branch("track_d0",&_track_d0);
        _tree->Branch("track_charge",&_track_charge);
        _tree->Branch("track_dxy",&_track_dxy);
	_tree->Branch("track_covariance",&_track_covariance);

    ///GenParticles
	_tree->Branch("gen_px",&_genp_px);
	_tree->Branch("gen_py",&_genp_py);
	_tree->Branch("gen_pz",&_genp_pz);
        _tree->Branch("gen_p",&_genp_p);
        _tree->Branch("gen_pt",&_genp_pt);
	_tree->Branch("gen_eta",&_genp_eta);
	_tree->Branch("gen_phi",&_genp_phi);
	_tree->Branch("gen_mass",&_genp_mass);
        _tree->Branch("gen_energy",&_genp_energy);
        _tree->Branch("gen_charge",&_genp_charge);
        _tree->Branch("gen_pdgid",&_genp_pdgid);
        _tree->Branch("gen_status",&_genp_status);
        
        //Trigger
        _tree->Branch("HLT_PFJet450", &_singlejet_450);
        //prescales
        _tree->Branch("pswgt_singlejet_450", &_pswgt_singlejet_450, "pswgt_singlejet_450/D");

}

// ------------ method called once each job just after ending the event loop  ------------
void
TreeProducer_AOD::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void
TreeProducer_AOD::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  triggerPathsVector.push_back("HLT_PFJet450_v");

  bool changedConfig = false;
  hltConfig_.init(iRun, iSetup, _trigResultsTag.process(), changedConfig);

  for (size_t i = 0; i < triggerPathsVector.size(); i++) {
    triggerPathsMap[triggerPathsVector[i]] = -1;
  }

  for(size_t i = 0; i < triggerPathsVector.size(); i++){
    TPRegexp pattern(triggerPathsVector[i]);
    for(size_t j = 0; j < hltConfig_.triggerNames().size(); j++){
      std::string pathName = hltConfig_.triggerNames()[j];
      if(TString(pathName).Contains(pattern)){
        triggerPathsMap[triggerPathsVector[i]] = j;
      }
    }
  }


  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,"HLT",changed) && hltPrescale_.init(iRun, iSetup, "HLT", changed) ) {
    if (changed) {
      // check if trigger names in (new) config
      cout<<"New trigger menu found !!!"<<endl;
      triggerIndex_.clear();
      const unsigned int n(hltConfig_.size());
      for(unsigned itrig=0;itrig<triggerNames_.size();itrig++) {
        triggerIndex_.push_back(hltConfig_.triggerIndex(triggerNames_[itrig]));
        //cout<<triggerNames_[itrig]<<" "<<triggerIndex_[itrig]<<" ";
        if (triggerIndex_[itrig] >= n)
          cout<<"does not exist in the current menu"<<endl;
        else
          cout<<"exists"<<endl;
      }
      cout << "Available TriggerNames are: " << endl;
      //if (true)
      //  hltConfig_.dump("Triggers");
    }
  }
  else {
    cout << "TreeProducer_AOD::analyze: config extraction failure with process name HLT" << endl;
}

}

// ------------ method called when ending the processing of a run  ------------
void
TreeProducer_AOD::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
TreeProducer_AOD::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
TreeProducer_AOD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TreeProducer_AOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void
TreeProducer_AOD::Init()
{
  _nEvent = _nRun = _nLumi = 0;

  //Trigger
  _singlejet_450 = 0;
  //prescales
  _pswgt_singlejet_450 = 1;

  // Vertices
  _vtx_N = 0;
  _vtx_N_stored = 0;
  _vtx_normalizedChi2.clear();
  _vtx_ndof.clear();
  _vtx_nTracks.clear();
  _vtx_d0.clear();
  _vtx_x.clear();
  _vtx_y.clear();
  _vtx_z.clear();
  _vtx_covariance.clear();
  _vtx_isFake.clear();
  _vtx_isValid.clear();


  //Tracks
  _nTrack = 0;
  _nTrack_stored = 0;
  _track_eta.clear();
  _track_fromPV.clear();
  _track_ndof.clear();
  _track_Nhits.clear();
  _track_normalizedChi2.clear();
  _track_NpixHits.clear();
  _track_phi.clear();
  _track_pt.clear();
  _track_px.clear();
  _track_py.clear();
  _track_pz.clear();  
  _track_ptError.clear();
  _track_dzError.clear();
  _track_dz.clear();
  _track_d0.clear();
  _track_charge.clear();
  _track_covariance.clear();
  _track_dxy.clear();
  _track_purity.clear();

  //GenParticles
  _genp_px.clear();
  _genp_py.clear();
  _genp_pz.clear();
  _genp_p.clear();
  _genp_pt.clear();  
  _genp_eta.clear();
  _genp_phi.clear();
  _genp_mass.clear();
  _genp_energy.clear();
  _genp_charge.clear();
  _genp_pdgid.clear();
  _genp_status.clear();
  

}

//define this as a plug-in
DEFINE_FWK_MODULE(TreeProducer_AOD);
