// -*- C++ -*-
//
// Package:    stage/ntuple
// Class:      ntuple
// 
/**\class ntuple ntuple.cc stage/ntuple/plugins/ntuple.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Collard Caroline
//         Created:  Mon, 14 Jan 2019 15:48:08 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"
#include "RecoTracker/DeDx/interface/DeDxTools.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/MuonSegment.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"


#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLink.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
//#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"

#include "HSCP_codeFromAnalysis.h"

#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPIsolation.h"

#include "TH1.h"
#include <TTree.h>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.
//

const int nMaxTrack = 10000;
const int nMaxDeDxH = 100000;
const int nMaxStrip = 100000;
const int nMaxStripprim = 100000;
const int nMaxSimHit = 100000;
const int nMaxGen = 1000;
const int nMaxMuon = 10000;
const int nMaxHSCP = 10000;


class ntuple : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ntuple(const edm::ParameterSet&);
      ~ntuple();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
       std::string m_format;
       edm::EDGetTokenT<reco::VertexCollection> m_primaryVertexTag;
       edm::EDGetTokenT<reco::TrackCollection>  m_tracksTag;
       edm::EDGetTokenT<std::vector<pat::IsolatedTrack>>  m_isotracksTag;
       edm::EDGetTokenT< susybsm::HSCParticleCollection  > m_hscpcand;
       edm::EDGetTokenT< susybsm::HSCPIsolationValueMap  > m_hscpiso;
//       edm::EDGetTokenT<std::vector<reco::PFCandidate> >  m_tracksMiniAOD;
       edm::EDGetTokenT<reco::DeDxHitInfoAss>   m_dedxTag;
       edm::EDGetTokenT<edm::Association<std::vector<reco::DeDxHitInfo> > >  m_dedxMiniTag;
       edm::EDGetTokenT<edm::ValueMap<int>> m_dedxPrescaleTag;
       edm::EDGetTokenT< edm::DetSetVector<StripDigiSimLink> > m_stripSimLink;
       edm::EDGetTokenT< std::vector<reco::GenParticle> > genParticlesToken_;
       edm::EDGetTokenT<reco::MuonCollection>  m_muonTag;
       edm::EDGetTokenT<reco::MuonTimeExtraMap>  m_muontimecb;
       edm::EDGetTokenT<reco::MuonTimeExtraMap>  m_muontimedt;
       edm::EDGetTokenT<reco::MuonTimeExtraMap>  m_muontimecsc;


//       edm::EDGetTokenT< edm::ValueMap<reco::DeDxData> > dEdxTrackToken_;

       bool m_runOnGS;
       int printOut_;
//       edm::InputTag trackTags_; //used to select what tracks to read from configuration file
//
       std::vector< edm::EDGetTokenT<CrossingFrame<PSimHit> > > cfTokens_;
       std::vector< edm::EDGetTokenT<std::vector<PSimHit> > > simHitTokens_;
       typedef std::pair<unsigned int, unsigned int> simHitCollectionID;
       typedef std::pair<simHitCollectionID, unsigned int> simhitAddr;
       typedef std::map<simHitCollectionID, std::vector<PSimHit> > simhit_collectionMap;
       simhit_collectionMap SimHitCollMap_;


       TH3F* dEdxTemplatesUncorr = NULL;
       TH3F* dEdxTemplatesCorr = NULL;


//       TH1D * histo; 
       TTree *smalltree;
       int      tree_runNumber ;
       int      tree_event ;
       int      tree_npv;

       int      tree_genpart;
       int      tree_gen_pdg[nMaxGen];
       float    tree_gen_pt[nMaxGen];
       float    tree_gen_eta[nMaxGen];
       float    tree_gen_phi[nMaxGen];
       float    tree_gen_mass[nMaxGen];
       bool     tree_gen_isHardProcess[nMaxGen];
       int      tree_gen_status[nMaxGen];
       int      tree_gen_moth_pdg[nMaxGen];
       int      tree_gen_ndaughter[nMaxGen];
       int      tree_gen_daughter_pdg[nMaxGen];

       int      tree_ntracks ;
       float    tree_track_pt[nMaxTrack];
       float    tree_track_pterr[nMaxTrack];
       float    tree_track_p[nMaxTrack];
       float    tree_track_eta[nMaxTrack];
       float    tree_track_phi[nMaxTrack];
       float    tree_track_charge[nMaxTrack];
       float    tree_track_chi2[nMaxTrack];
       int      tree_track_nvalidhits[nMaxTrack];
       int      tree_track_npixhits[nMaxTrack];
       int      tree_track_missing[nMaxTrack];
       float    tree_track_validfrac[nMaxTrack];
       float    tree_track_validlast[nMaxTrack];
       int      tree_track_qual[nMaxTrack];
//       float    tree_track_dedx_harmonic2[nMaxTrack];
       int      tree_track_index_hit[nMaxTrack];
       int      tree_track_nhits[nMaxTrack];
       int      tree_track_prescale[nMaxTrack];
       float    tree_track_ih_ampl[nMaxTrack];
       float    tree_track_ih_ampl_corr[nMaxTrack];
       float    tree_track_ias_ampl[nMaxTrack];
       float    tree_track_ias_ampl_corr[nMaxTrack];

       int      tree_dedxhits ;
       uint32_t tree_dedx_detid[nMaxDeDxH];
       int      tree_dedx_subdetid[nMaxDeDxH];
       int      tree_dedx_modulgeom[nMaxDeDxH];
       float    tree_dedx_charge[nMaxDeDxH];
       float    tree_dedx_pathlength[nMaxDeDxH];
       float    tree_dedx_posx[nMaxDeDxH];
       float    tree_dedx_posy[nMaxDeDxH];
       float    tree_dedx_posz[nMaxDeDxH];
       bool     tree_dedx_strip[nMaxDeDxH];
       bool     tree_dedx_pixel[nMaxDeDxH];
       bool     tree_dedx_insideTkMod[nMaxDeDxH];
    
       int      tree_sclus_index_strip[nMaxDeDxH];
       int      tree_sclus_nstrip[nMaxDeDxH];
       int      tree_sclus_firstsclus[nMaxDeDxH];
       float    tree_sclus_barycenter[nMaxDeDxH];
       float    tree_sclus_charge[nMaxDeDxH];
       float    tree_sclus_errorclus[nMaxDeDxH];
       bool     tree_sclus_ismerged[nMaxDeDxH];
       bool     tree_sclus_sat254[nMaxDeDxH];
       bool     tree_sclus_sat255[nMaxDeDxH];
       bool     tree_sclus_shape[nMaxDeDxH];
       int      tree_sclus_index_strip_corr[nMaxDeDxH];
       int      tree_sclus_nstrip_corr[nMaxDeDxH];
       float    tree_sclus_charge_corr[nMaxDeDxH];
       bool     tree_sclus_clusclean[nMaxDeDxH];
       bool     tree_sclus_clusclean2[nMaxDeDxH];
  
       int      tree_nstrips;
       int      tree_strip_ampl[nMaxStrip];
       int      tree_nstrips_corr;
       int      tree_strip_ampl_corr[nMaxStripprim];

       int      tree_sclus_index_simhit[nMaxDeDxH];
       int      tree_sclus_nsimhit[nMaxDeDxH];
       float    tree_sclus_eloss[nMaxDeDxH];
       int      tree_nsimhits;
       int      tree_simhit_pid[nMaxSimHit];
       int      tree_simhit_process[nMaxSimHit];
       float    tree_simhit_p[nMaxSimHit];
       float    tree_simhit_eloss[nMaxSimHit];
       float    tree_simhit_tof[nMaxSimHit];
       float    tree_simhit_segment[nMaxSimHit];
       float    tree_simhit_xentry[nMaxSimHit];
       float    tree_simhit_yentry[nMaxSimHit];
       float    tree_simhit_zentry[nMaxSimHit];
       float    tree_simhit_xexit[nMaxSimHit];
       float    tree_simhit_yexit[nMaxSimHit];
       float    tree_simhit_zexit[nMaxSimHit];
       
       int      tree_nmuons ;
       float    tree_muon_pt[nMaxMuon];
       float    tree_muon_ptSA[nMaxMuon];
       float    tree_muon_ptIT[nMaxMuon];
       float    tree_muon_p[nMaxMuon];
       float    tree_muon_eta[nMaxMuon];
       float    tree_muon_phi[nMaxMuon];
       float    tree_muon_comb_inversebeta[nMaxMuon];
       float    tree_muon_comb_inversebetaerr[nMaxMuon];
       int      tree_muon_comb_tofndof[nMaxMuon];
       float    tree_muon_comb_vertextime[nMaxMuon];
       float    tree_muon_dt_inversebeta[nMaxMuon];
       float    tree_muon_dt_inversebetaerr[nMaxMuon];
       int      tree_muon_dt_tofndof[nMaxMuon];
       float    tree_muon_dt_vertextime[nMaxMuon];
       float    tree_muon_csc_inversebeta[nMaxMuon];
       float    tree_muon_csc_inversebetaerr[nMaxMuon];
       int      tree_muon_csc_tofndof[nMaxMuon];
       float    tree_muon_csc_vertextime[nMaxMuon];

       int      tree_hscp ;
       int      tree_hscp_gen_id[nMaxHSCP];
       float    tree_hscp_gen_dr[nMaxHSCP];
       int      tree_hscp_track_idx[nMaxHSCP];
       int      tree_hscp_muon_idx[nMaxHSCP];
       float    tree_hscp_iso_tk[nMaxHSCP];
       float    tree_hscp_iso_ecal[nMaxHSCP];
       float    tree_hscp_iso_hcal[nMaxHSCP];


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ntuple::ntuple(const edm::ParameterSet& iConfig)
/*
 :
  trackTags_(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))

{
*/
{
   m_format   = iConfig.getParameter<std::string>("format_file");
   m_primaryVertexTag   = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexColl"));
   m_tracksTag = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"));
   m_isotracksTag = consumes<std::vector<pat::IsolatedTrack>>(iConfig.getParameter<edm::InputTag>("isotracks"));
//   m_tracksMiniAOD = consumes<std::vector<reco::PFCandidate> >(iConfig.getParameter<edm::InputTag>("PFCand"));
   m_hscpcand = consumes<susybsm::HSCParticleCollection >(iConfig.getParameter<edm::InputTag>("collectionHSCP"));
   m_hscpiso = consumes< susybsm::HSCPIsolationValueMap  >(iConfig.getParameter<edm::InputTag>("isoHSCP"));
   m_dedxTag = consumes<reco::DeDxHitInfoAss>(iConfig.getParameter<edm::InputTag>("dedx"));
   m_dedxMiniTag = consumes<edm::Association<std::vector<reco::DeDxHitInfo> >>(iConfig.getParameter<edm::InputTag>("MiniDedx"));
   m_dedxPrescaleTag = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("dEdxHitInfoPrescale"));
   m_runOnGS = iConfig.getParameter<bool>("runOnGS");
   printOut_ = iConfig.getUntrackedParameter<int>("printOut");
   m_stripSimLink = consumes< edm::DetSetVector<StripDigiSimLink> >(iConfig.getParameter<edm::InputTag>("stripSimLinks"));
   std::vector<std::string> trackerContainers(iConfig.getParameter<std::vector<std::string> >("ROUList"));
   cfTokens_.reserve(trackerContainers.size());
   simHitTokens_.reserve(trackerContainers.size());
   for(auto const& trackerContainer : trackerContainers) {
      cfTokens_.push_back(consumes<CrossingFrame<PSimHit> >(edm::InputTag("mix", trackerContainer)));
      simHitTokens_.push_back(consumes<std::vector<PSimHit> >(edm::InputTag("g4SimHits", trackerContainer)));
   }
   genParticlesToken_ = consumes< std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("GenPart"));
//   dEdxTrackToken_ = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));

   //loadDeDxTemplates
   dEdxTemplatesUncorr = loadDeDxTemplate ("/opt/sbg/data/safe1/cms/ccollard/HSCP/CMSSW_9_4_3/src/stage/ntuple/test/minbias_template_uncorr_iter1.root", true);
   dEdxTemplatesCorr = loadDeDxTemplate ("/opt/sbg/data/safe1/cms/ccollard/HSCP/CMSSW_9_4_3/src/stage/ntuple/test/minbias_template_corr_iter1.root", true);
   m_muonTag = consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
   m_muontimecb = consumes< reco::MuonTimeExtraMap >(iConfig.getParameter<edm::InputTag>("muonTOF"));
   m_muontimedt = consumes< reco::MuonTimeExtraMap >(iConfig.getParameter<edm::InputTag>("muonTDT"));
   m_muontimecsc = consumes< reco::MuonTimeExtraMap >(iConfig.getParameter<edm::InputTag>("muonTCSC"));


   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;
   //histo = fs->make<TH1D>("charge" , "Charges" , 200 , -2 , 2 );

   smalltree = fs->make<TTree>("ttree", "ttree");
   smalltree -> Branch ( "runNumber", &tree_runNumber ) ;
   smalltree -> Branch ( "event",     &tree_event ) ;
   smalltree -> Branch ( "npv",     &tree_npv ) ;


   smalltree -> Branch ( "ngenpart",  &tree_genpart) ;
   smalltree -> Branch ( "gen_pdg",   tree_gen_pdg,    "gen_pdg[ngenpart]/I");
   smalltree -> Branch ( "gen_pt",    tree_gen_pt,     "gen_pt[ngenpart]/F");
   smalltree -> Branch ( "gen_eta",   tree_gen_eta,    "gen_eta[ngenpart]/F");
   smalltree -> Branch ( "gen_phi",   tree_gen_phi,    "gen_phi[ngenpart]/F");
   smalltree -> Branch ( "gen_mass",   tree_gen_mass,    "gen_mass[ngenpart]/F");
   smalltree -> Branch ( "gen_isHardProcess",   tree_gen_isHardProcess,    "gen_isHardProcess[ngenpart]/O");
   smalltree -> Branch ( "gen_status",   tree_gen_status,    "gen_status[ngenpart]/I");
   smalltree -> Branch ( "gen_moth_pdg",   tree_gen_moth_pdg,    "gen_moth_pdg[ngenpart]/I");
   smalltree -> Branch ( "gen_ndaughter",   tree_gen_ndaughter,    "gen_ndaughter[ngenpart]/I");
   smalltree -> Branch ( "gen_daughter_pdg",   tree_gen_daughter_pdg,    "gen_daughter_pdg[ngenpart]/I");
  

   smalltree -> Branch ( "ntracks",              &tree_ntracks ) ;
   smalltree -> Branch ( "track_pt",             tree_track_pt,             "track_pt[ntracks]/F" );
   smalltree -> Branch ( "track_pterr",          tree_track_pterr,          "track_pterr[ntracks]/F" );
   smalltree -> Branch ( "track_p",              tree_track_p,              "track_p[ntracks]/F"  );
   smalltree -> Branch ( "track_eta",            tree_track_eta,            "track_eta[ntracks]/F" );
   smalltree -> Branch ( "track_phi",            tree_track_phi,            "track_phi[ntracks]/F" );
   smalltree -> Branch ( "track_charge",         tree_track_charge,         "track_charge[ntracks]/F" );
   smalltree -> Branch ( "track_chi2",           tree_track_chi2,           "track_chi2[ntracks]/F" );
   smalltree -> Branch ( "track_nvalidhits",     tree_track_nvalidhits,     "track_nvalidhits[ntracks]/I" );
   smalltree -> Branch ( "track_npixhits",       tree_track_npixhits,       "track_npixhits[ntracks]/I" );
   smalltree -> Branch ( "track_missing",        tree_track_missing,        "track_missing[ntracks]/I" );
   smalltree -> Branch ( "track_validfraction",  tree_track_validfrac,      "track_validfraction[ntracks]/F" );
   smalltree -> Branch ( "track_validlast",      tree_track_validlast,      "track_validlast[ntracks]/F" );
   smalltree -> Branch ( "track_qual",           tree_track_qual,           "track_qual[ntracks]/I" );
//   smalltree -> Branch ( "track_dedx_harmonic2", tree_track_dedx_harmonic2, "track_dedx_harmonic2[ntracks]/F" );
   smalltree -> Branch ( "track_index_hit",      tree_track_index_hit,      "track_index_hit[ntracks]/I" );
   smalltree -> Branch ( "track_nhits",          tree_track_nhits,          "track_nhits[ntracks]/I"  );
   smalltree -> Branch ( "track_prescale",       tree_track_prescale,       "track_prescale[ntracks]/I" );
   smalltree -> Branch ( "track_ih_ampl",        tree_track_ih_ampl,        "track_ih_ampl[ntracks]/F" );
   smalltree -> Branch ( "track_ih_ampl_corr",   tree_track_ih_ampl_corr,   "track_ih_ampl_corr[ntracks]/F" );
   smalltree -> Branch ( "track_ias_ampl",       tree_track_ias_ampl,       "track_ias_ampl[ntracks]/F" );
   smalltree -> Branch ( "track_ias_ampl_corr",  tree_track_ias_ampl_corr,  "track_ias_ampl_corr[ntracks]/F" );

   smalltree -> Branch ( "ndedxhits",        &tree_dedxhits ) ;
   smalltree -> Branch ( "dedx_detid",       tree_dedx_detid,      "dedx_detid[ndedxhits]/i" );
   smalltree -> Branch ( "dedx_subdetid",    tree_dedx_subdetid,   "dedx_subdetid[ndedxhits]/I" );
   smalltree -> Branch ( "dedx_modulgeom",   tree_dedx_modulgeom,  "dedx_modulgeom[ndedxhits]/I" );
   smalltree -> Branch ( "dedx_charge",      tree_dedx_charge,     "dedx_charge[ndedxhits]/F" );
   smalltree -> Branch ( "dedx_pathlength",  tree_dedx_pathlength, "dedx_pathlength[ndedxhits]/F" );
   smalltree -> Branch ( "dedx_posx",        tree_dedx_posx,       "dedx_posx[ndedxhits]/F" );
   smalltree -> Branch ( "dedx_posy",        tree_dedx_posy,       "dedx_posy[ndedxhits]/F" );
   smalltree -> Branch ( "dedx_posz",        tree_dedx_posz,       "dedx_posz[ndedxhits]/F"  );
   smalltree -> Branch ( "dedx_isstrip",     tree_dedx_strip,      "dedx_isstrip[ndedxhits]/O" );
   smalltree -> Branch ( "dedx_ispixel",     tree_dedx_pixel,      "dedx_ispixel[ndedxhits]/O" );
   smalltree -> Branch ( "dedx_insideTkMod", tree_dedx_insideTkMod,"dedx_insideTkMod[ndedxhits]/O" );

   smalltree -> Branch ( "sclus_firstsclus",  tree_sclus_firstsclus, "sclus_firstsclus[ndedxhits]/I" );
   smalltree -> Branch ( "sclus_barycenter",  tree_sclus_barycenter, "sclus_barycenter[ndedxhits]/F" );
   smalltree -> Branch ( "sclus_charge",      tree_sclus_charge,     "sclus_charge[ndedxhits]/F" );
   smalltree -> Branch ( "sclus_errorclus",   tree_sclus_errorclus,  "sclus_errorclus[ndedxhits]/F"  );
   smalltree -> Branch ( "sclus_ismerged",    tree_sclus_ismerged,   "sclus_ismerged[ndedxhits]/O" );
   smalltree -> Branch ( "sclus_index_strip", tree_sclus_index_strip,"sclus_index_strip[ndedxhits]/I"  );
   smalltree -> Branch ( "sclus_nstrip",      tree_sclus_nstrip,     "sclus_nstrip[ndedxhits]/I"  );
   smalltree -> Branch ( "sclus_sat254",      tree_sclus_sat254,     "sclus_sat254[ndedxhits]/O" );
   smalltree -> Branch ( "sclus_sat255",      tree_sclus_sat255,     "sclus_sat255[ndedxhits]/O" );
   smalltree -> Branch ( "sclus_shape",       tree_sclus_shape,      "sclus_shape[ndedxhits]/O" );
   smalltree -> Branch ( "sclus_index_strip_corr", tree_sclus_index_strip_corr,"sclus_index_strip_corr[ndedxhits]/I"  );
   smalltree -> Branch ( "sclus_nstrip_corr", tree_sclus_nstrip_corr,"sclus_nstrip_corr[ndedxhits]/I"  );
   smalltree -> Branch ( "sclus_charge_corr", tree_sclus_charge_corr,"sclus_charge_corr[ndedxhits]/F" );
   smalltree -> Branch ( "sclus_clusclean",   tree_sclus_clusclean,  "sclus_clusclean[ndedxhits]/O" );
   smalltree -> Branch ( "sclus_clusclean2",   tree_sclus_clusclean2,  "sclus_clusclean2[ndedxhits]/O" );

   smalltree -> Branch ( "sclus_index_simhit", tree_sclus_index_simhit,"sclus_index_simhit[ndedxhits]/I"  );
   smalltree -> Branch ( "sclus_nsimhit",      tree_sclus_nsimhit,     "sclus_nsimhit[ndedxhits]/I"  );
   smalltree -> Branch ( "sclus_eloss",        tree_sclus_eloss,       "sclus_eloss[ndedxhits]/F"  );

   smalltree -> Branch ( "nstrips",    &tree_nstrips ) ;
   smalltree -> Branch ( "strip_ampl", tree_strip_ampl, "strip_ampl[nstrips]/I");
   smalltree -> Branch ( "nstrips_corr",    &tree_nstrips_corr ) ;
   smalltree -> Branch ( "strip_ampl_corr", tree_strip_ampl_corr, "strip_ampl_corr[nstrips_corr]/I");

   smalltree -> Branch ( "nsimhits",       &tree_nsimhits ) ;
   smalltree -> Branch ( "simhit_pid",     tree_simhit_pid,     "simhit_pid[nsimhits]/I");
   smalltree -> Branch ( "simhit_process", tree_simhit_process, "simhit_process[nsimhits]/I");
   smalltree -> Branch ( "simhit_p",       tree_simhit_p,       "simhit_p[nsimhits]/F");
   smalltree -> Branch ( "simhit_eloss",   tree_simhit_eloss,   "simhit_eloss[nsimhits]/F");
   smalltree -> Branch ( "simhit_tof",     tree_simhit_tof,     "simhit_tof[nsimhits]/F");
   smalltree -> Branch ( "simhit_segment", tree_simhit_segment, "simhit_segment[nsimhits]/F");
   smalltree -> Branch ( "simhit_xentry",  tree_simhit_xentry,  "simhit_xentry[nsimhits]/F");
   smalltree -> Branch ( "simhit_yentry",  tree_simhit_yentry,  "simhit_yentry[nsimhits]/F");
   smalltree -> Branch ( "simhit_zentry",  tree_simhit_zentry,  "simhit_zentry[nsimhits]/F");
   smalltree -> Branch ( "simhit_xexit",   tree_simhit_xexit,   "simhit_xexit[nsimhits]/F");
   smalltree -> Branch ( "simhit_yexit",   tree_simhit_yexit,   "simhit_yexit[nsimhits]/F");
   smalltree -> Branch ( "simhit_zexit",   tree_simhit_zexit,   "simhit_zexit[nsimhits]/F");

   smalltree -> Branch ( "nmuons",              &tree_nmuons ) ;
   smalltree -> Branch ( "muon_pt",             tree_muon_pt,             "muon_pt[nmuons]/F" );
   smalltree -> Branch ( "muon_ptSA",           tree_muon_ptSA,           "muon_ptSA[nmuons]/F" );
   smalltree -> Branch ( "muon_ptIT",           tree_muon_ptIT,           "muon_ptIT[nmuons]/F" );
   smalltree -> Branch ( "muon_p",              tree_muon_p,              "muon_p[nmuons]/F"  );
   smalltree -> Branch ( "muon_eta",            tree_muon_eta,            "muon_eta[nmuons]/F" );
   smalltree -> Branch ( "muon_phi",            tree_muon_phi,            "muon_phi[nmuons]/F" );
   smalltree -> Branch ( "muon_comb_inversebeta",    tree_muon_comb_inversebeta,    "muon_comb_inversebeta[nmuons]/F" );
   smalltree -> Branch ( "muon_comb_inversebetaerr", tree_muon_comb_inversebetaerr, "muon_comb_inversebetaerr[nmuons]/F" );
   smalltree -> Branch ( "muon_comb_tofndof",        tree_muon_comb_tofndof,        "muon_comb_tofndof[nmuons]/I" );
   smalltree -> Branch ( "muon_comb_vertextime",     tree_muon_comb_vertextime,     "muon_comb_vertextime[nmuons]/I" );
   smalltree -> Branch ( "muon_dt_inversebeta",      tree_muon_dt_inversebeta,      "muon_dt_inversebeta[nmuons]/F" );
   smalltree -> Branch ( "muon_dt_inversebetaerr",   tree_muon_dt_inversebetaerr,   "muon_dt_inversebetaerr[nmuons]/F" );
   smalltree -> Branch ( "muon_dt_tofndof",          tree_muon_dt_tofndof,          "muon_dt_tofndof[nmuons]/I" );
   smalltree -> Branch ( "muon_dt_vertextime",       tree_muon_dt_vertextime,       "muon_dt_vertextime[nmuons]/I" );
   smalltree -> Branch ( "muon_csc_inversebeta",     tree_muon_csc_inversebeta,     "muon_csc_inversebeta[nmuons]/F" );
   smalltree -> Branch ( "muon_csc_inversebetaerr",  tree_muon_csc_inversebetaerr,  "muon_csc_inversebetaerr[nmuons]/F" );
   smalltree -> Branch ( "muon_csc_tofndof",         tree_muon_csc_tofndof,         "muon_csc_tofndof[nmuons]/I" );
   smalltree -> Branch ( "muon_csc_vertextime",      tree_muon_csc_vertextime,      "muon_csc_vertextime[nmuons]/I" );

   smalltree -> Branch ( "nhscp",              &tree_hscp ) ;
   smalltree -> Branch ( "hscp_gen_id",        tree_hscp_gen_id,             "hscp_gen_id[nhscp]/I" );
   smalltree -> Branch ( "hscp_gen_dr",        tree_hscp_gen_dr,             "hscp_gen_dr[nhscp]/F" );
   smalltree -> Branch ( "hscp_track_idx",     tree_hscp_track_idx,          "hscp_track_idx[nhscp]/I" );
   smalltree -> Branch ( "hscp_muon_idx",      tree_hscp_muon_idx,           "hscp_muon_idx[nhscp]/I" );
   smalltree -> Branch ( "hscp_iso_tk",        tree_hscp_iso_tk,             "hscp_iso_tk[nhscp]/F" );
   smalltree -> Branch ( "hscp_iso_ecal",      tree_hscp_iso_ecal,           "hscp_iso_ecal[nhscp]/F" );
   smalltree -> Branch ( "hscp_iso_hcal",      tree_hscp_iso_hcal,           "hscp_iso_hcal[nhscp]/F" );
}


ntuple::~ntuple()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ntuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace reco;

    using reco::TrackCollection;


    EventID myEvId = iEvent.id();
    tree_runNumber = myEvId.run();
    tree_event = myEvId.event();

    
    edm::Handle<reco::VertexCollection> primaryVertex ;
    iEvent.getByToken(m_primaryVertexTag,primaryVertex);
    tree_npv = primaryVertex->size();


    edm::Handle<reco::TrackCollection> trackCollectionHandle;
//    edm::Handle<std::vector<reco::PFCandidate> > pCands;
    edm::Handle<std::vector<pat::IsolatedTrack>>  IsotrackCollectionHandle;
    if (m_format != "miniAOD" ) {iEvent.getByToken(m_tracksTag,trackCollectionHandle); }
//    else {iEvent.getByToken(m_tracksMiniAOD,pCands);}
    else {iEvent.getByToken(m_isotracksTag,IsotrackCollectionHandle);}


    // inspired by /opt/sbg/data/safe1/cms/ccollard/HSCP/CMSSW_9_4_3/src/SUSYBSMAnalysis-HSCP/test/UsefulScripts/DeDxStudy/DeDxStudy
    edm::Handle<reco::DeDxHitInfoAss> dedxCollH;
    edm::Handle<edm::Association<std::vector<reco::DeDxHitInfo> > > dedxMiniH;
    if (m_format != "miniAOD" ) { 
            iEvent.getByToken(m_dedxTag,dedxCollH);
            if (!dedxCollH.isValid()) std::cout << " access problem to dedxHitInfo collection " << std::endl;
    } 
    else {
            iEvent.getByToken(m_dedxMiniTag,dedxMiniH);
            if (!dedxMiniH.isValid()) std::cout << " access problem to dedxHitInfo collection " << std::endl;
    }

    // inspired by https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/plugins/PATIsolatedTrackProducer.cc
    edm::Handle<edm::ValueMap<int>> dedxHitInfoPrescale;
    if (m_format != "miniAOD" ) {
     iEvent.getByToken(m_dedxPrescaleTag, dedxHitInfoPrescale);
    }


    // inspired by SUSYBSMAnalysis/HSCP/test/AnalysisCode/Analysis_Step1_EventLoop.C
    edm::Handle<susybsm::HSCParticleCollection> hscpCollH;
    edm::Handle<susybsm::HSCPIsolationValueMap> isoHSCPColl;
    if (m_format != "miniAOD" ) {
      iEvent.getByToken(m_hscpcand,hscpCollH);
      if(!hscpCollH.isValid()){printf("HSCParticle Collection NotFound\n");}
      iEvent.getByToken(m_hscpiso,isoHSCPColl);
      if(!isoHSCPColl.isValid()){printf("HSCPIsolation Collection NotFound\n");}
    }




    edm::Handle<reco::MuonCollection> muonCollectionHandle;
    if (m_format != "miniAOD" ) {
     iEvent.getByToken(m_muonTag,muonCollectionHandle);
    }


    edm::Handle<reco::MuonTimeExtraMap> TOFCollH;
//    iEvent.getByLabel("muons", "combined",TOFCollH);
//    if (m_format != "miniAOD" ) iEvent.getByLabel("muons","combined",TOFCollH);
    if (m_format != "miniAOD" ) {
     iEvent.getByToken(m_muontimecb,TOFCollH);
     if (!TOFCollH.isValid()) { std::cout << " access problem to TOF collection " << std::endl; }
    }
    edm::Handle<reco::MuonTimeExtraMap> TOFDTCollH;
    if (m_format != "miniAOD" ) {
      iEvent.getByToken(m_muontimedt,TOFDTCollH);
      if (!TOFDTCollH.isValid()) std::cout << " access problem to DT TOF collection " << std::endl;
    }
    edm::Handle<reco::MuonTimeExtraMap> TOFCSCCollH;
    if (m_format != "miniAOD" ) {
     iEvent.getByToken(m_muontimecsc,TOFCSCCollH);
     if (!TOFCSCCollH.isValid()) std::cout << " access problem to CSC TOF collection " << std::endl;
    }


    // inspired by https://github.com/cms-sw/cmssw/blob/master/Validation/TrackerRecHits/test/StripClusterMCanalysis.cc
    edm::Handle< edm::DetSetVector<StripDigiSimLink> >  stripdigisimlink;
    iEvent.getByToken(m_stripSimLink, stripdigisimlink);

    SimHitCollMap_.clear();
    for(auto const& cfToken : cfTokens_) {
      edm::Handle<CrossingFrame<PSimHit> > cf_simhit;
      int Nhits = 0;
      if (iEvent.getByToken(cfToken, cf_simhit)) {
      std::unique_ptr<MixCollection<PSimHit> > thisContainerHits(new MixCollection<PSimHit>(cf_simhit.product()));
      for (auto const& isim : *thisContainerHits) {
        DetId theDet(isim.detUnitId());
        edm::EDConsumerBase::Labels labels;
        iEvent.labelsForToken(cfToken, labels);
        std::string trackerContainer(labels.productInstance);
        if (printOut_ && Nhits==0) std::cout << "  trackerContainer " << trackerContainer << std::endl;
        unsigned int tofBin = StripDigiSimLink::LowTof;
        if (trackerContainer.find(std::string("HighTof")) != std::string::npos) tofBin = StripDigiSimLink::HighTof;
        simHitCollectionID theSimHitCollID = std::make_pair(theDet.subdetId(), tofBin);
        SimHitCollMap_[theSimHitCollID].push_back(isim);
        ++Nhits;
       }
       if (printOut_ > 0) std::cout << "simHits from crossing frames; map size = " << SimHitCollMap_.size()
                                   << ", Hit count = " << Nhits << ", " << sizeof(SimHitCollMap_)
                                   << " bytes" << std::endl;
       }
    }


    for(auto const& simHitToken : simHitTokens_) {
     edm::Handle<std::vector<PSimHit> > simHits;
     int Nhits = 0;
     if(iEvent.getByToken(simHitToken, simHits)) {
      for (auto const& isim : *simHits) {
        DetId theDet(isim.detUnitId());
        edm::EDConsumerBase::Labels labels;
        iEvent.labelsForToken(simHitToken, labels);
        std::string trackerContainer(labels.productInstance);
        if (printOut_>0 && Nhits==0) std::cout << "  trackerContainer " << trackerContainer << std::endl;
        unsigned int tofBin = StripDigiSimLink::LowTof;
        if (trackerContainer.find(std::string("HighTof")) != std::string::npos) tofBin = StripDigiSimLink::HighTof;
        simHitCollectionID theSimHitCollID = std::make_pair(theDet.subdetId(), tofBin);
        SimHitCollMap_[theSimHitCollID].push_back(isim);
        ++Nhits;
      }
      if (printOut_ > 0) std::cout << "simHits from hard-scatter collection; map size = " << SimHitCollMap_.size()
                                   << ", Hit count = " << Nhits << ", " << sizeof(SimHitCollMap_)
                                   << " bytes" << std::endl;
     }
    }


    // inspired by SUSYBSMAnalysis-HSCP/plugins/HSCPValidator.cc
    
/*
    Handle<ValueMap<DeDxData> >          dEdxTrackHandle;
    iEvent.getByToken(dEdxTrackToken_, dEdxTrackHandle);
    const ValueMap<DeDxData> dEdxTrack = *dEdxTrackHandle.product();
    if (!dEdxTrackHandle.isValid()) std::cout << " access problem to DeDxData ValueMap " << std::endl;
*/

    // inspired by Stephanie (mail of Feb 6, 2019) 
    // and https://github.com/jozzez1/cmssw/blob/Run2HSCP16_v4/SUSYBSMAnalysis/HSCP/test/UsefulScripts/DeDxStudy/DeDxStudy.C#L698
    // https://github.com/CMS-HSCP/SUSYBSMAnalysis-HSCP/blob/master/plugins/BigNtuplizer.cc
    
    tree_genpart=0;
//    if (m_runOnGS) {
    edm::Handle< std::vector<reco::GenParticle> > GenColl;
    iEvent.getByToken(genParticlesToken_, GenColl);


    int n_genp = GenColl->size();
    if (n_genp>0) {
     for(int i=0;i< n_genp ;++i){
      const reco::GenParticle* genCand = &(*GenColl)[i];
/*
      if(genCand->pt()<5)continue;
      if(genCand->status()!=1)continue;
      int AbsPdg=abs(genCand->pdgId());
      if(AbsPdg<1000000 && AbsPdg!=17)continue;
*/

      if (genCand->pt()<5) continue;
/*
      int absPdg = abs(genCand->pdgId());
      if (absPdg>1000000 || absPdg==17) {
          std::cout << " particle " << absPdg << "     status "  << genCand->status()          << std::endl;
      } 
*/
      if (tree_genpart<nMaxGen) {
       tree_gen_pdg[tree_genpart]=genCand->pdgId();
       tree_gen_pt[tree_genpart]=genCand->pt();
       tree_gen_eta[tree_genpart]=genCand->eta();
       tree_gen_phi[tree_genpart]=genCand->phi();
       tree_gen_mass[tree_genpart]=genCand->mass();
       tree_gen_isHardProcess[tree_genpart]=genCand->isHardProcess();
       tree_gen_status[tree_genpart]=genCand->status();
       if (genCand->numberOfMothers()>0) tree_gen_moth_pdg[tree_genpart]=genCand->mother()->pdgId();
       else tree_gen_moth_pdg[tree_genpart]=-9999;
       tree_gen_ndaughter[tree_genpart]=genCand->numberOfDaughters();
       if (genCand->numberOfDaughters()>0) tree_gen_daughter_pdg[tree_genpart]=genCand->daughter(0)->pdgId();
       else tree_gen_daughter_pdg[tree_genpart]=-9999;
       tree_genpart++;
      }
     }
    } // if n_genp
//    } // if runOnGS

     //
    tree_ntracks=0;
    tree_dedxhits=0;
    tree_nstrips=0;
    tree_nstrips_corr=0;
    tree_nsimhits=0;
    tree_nmuons=0;
     
    std::vector<reco::TrackRef> KeepTrackRefVec;
    unsigned int ntrack_in_coll=0;
    if (m_format != "miniAOD" ) {ntrack_in_coll=trackCollectionHandle->size();}
//    else {ntrack_in_coll=pCands->size();} 
    else {ntrack_in_coll=IsotrackCollectionHandle->size();} 

//    for(unsigned int c=0;c<trackCollectionHandle->size();c++){
    for(unsigned int c=0;c<ntrack_in_coll;c++){
//      reco::TrackRef track;

      float pt_tr=0;
      float pterr_tr=0;
      float p_tr=0;
      float eta_tr=0;
      float phi_tr=0;
      float charge_tr=0;
      float chi2_tr=0;
      int nh_tr=0;
      int nhpix_tr=0;
      int nmisstilL_tr=0;
      float frac_tr=-1;
      float frac2_tr=-1;
      int qual_tr=0;
      const reco::DeDxHitInfo* dedxHits = nullptr;
      reco::DeDxHitInfoRef dedxHitsRef;

      if (m_format != "miniAOD" ) {  // AOD on RECO General Tracks

          reco::TrackRef  track = reco::TrackRef( trackCollectionHandle.product(), c ); 
          if (printOut_ > 0) std::cout << " track with pT =  " << track->pt() << std::endl;
          pt_tr=track->pt();
          pterr_tr=track->ptError();
          p_tr=track->p();
          eta_tr=track->eta();
          phi_tr=track->phi();
          charge_tr=track->charge();
          chi2_tr=track->chi2()/track->ndof();
          nh_tr=track->numberOfValidHits();
          nhpix_tr=track->hitPattern().numberOfValidPixelHits();
          nmisstilL_tr=track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_INNER_HITS) + track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);
          frac_tr=track->validFraction();
          frac2_tr=track->found()<=0?-1:track->found() / float(track->found() + nmisstilL_tr);
          qual_tr=track->qualityMask();
          dedxHitsRef = dedxCollH->get(track.key());
          if(!dedxHitsRef.isNull())dedxHits = &(*dedxHitsRef);
          if(!dedxHitsRef.isNull() && tree_ntracks < nMaxTrack) KeepTrackRefVec.push_back(track);
      }
      else {   // miniAOD on PAT Isolated Tracks
//          auto &track = (*IsotrackCollectionHandle)[c];
          edm::Ref<std::vector<pat::IsolatedTrack> > track = edm::Ref<std::vector<pat::IsolatedTrack> >( IsotrackCollectionHandle, c );
          if (printOut_ > 0) std::cout << " track with pT =  " << track->pt() << std::endl;
          pt_tr=track->pt();
          p_tr=track->p();
          eta_tr=track->eta();
          phi_tr=track->phi();
//          chi2_tr=track.chi2()/track.ndof();
//          nh_tr=track.numberOfValidHits();
          dedxHitsRef = (*dedxMiniH)[track];
          if(!dedxHitsRef.isNull())dedxHits = &(*dedxHitsRef);
      }
//     for(reco::TrackCollection::const_iterator track = trackCollectionHandle->begin(); track != trackCollectionHandle->end(); ++track)
      //basic track quality cuts
      ////if(track.isNull())continue;
      //if(track->chi2()/track->ndof()>5 )continue;
      //if(track->found()<8) continue;
      //if(track->pt() < 5) continue;

      if (tree_ntracks < nMaxTrack) {


/*
       tree_track_pt[tree_ntracks]= track->pt();
       tree_track_p[tree_ntracks]= track->p();
       tree_track_eta[tree_ntracks]= track->eta();
       tree_track_phi[tree_ntracks]= track->phi();
       tree_track_chi2[tree_ntracks]= track->chi2()/track->ndof();
       tree_track_nvalidhits[tree_ntracks]= track->numberOfValidHits();
*/
       tree_track_pt[tree_ntracks]= pt_tr;
       tree_track_pterr[tree_ntracks]= pterr_tr;
       tree_track_p[tree_ntracks]=  p_tr;
       tree_track_eta[tree_ntracks]= eta_tr;
       tree_track_phi[tree_ntracks]= phi_tr;
       tree_track_charge[tree_ntracks]= charge_tr; 
       tree_track_chi2[tree_ntracks]= chi2_tr; 
       tree_track_nvalidhits[tree_ntracks]= nh_tr;
       tree_track_npixhits[tree_ntracks]= nhpix_tr;
       tree_track_missing[tree_ntracks]= nmisstilL_tr;
       tree_track_validfrac[tree_ntracks]= frac_tr;
       tree_track_validlast[tree_ntracks]= frac2_tr;
       tree_track_qual[tree_ntracks]= qual_tr;


       //load dEdx informations
       if(!dedxHits)continue;

       if (m_format != "miniAOD" ) {  // AOD on RECO General Tracks
        tree_track_prescale[tree_ntracks] = (*dedxHitInfoPrescale)[dedxHitsRef];
       }
       else {
        tree_track_prescale[tree_ntracks] = -1;
       }

       // load the dedx estimator
       /*
       if (dEdxTrackHandle.isValid()) {
         if (printOut_ > 0) std::cout << " valid dEdxTrackHandle " << std::endl;
         if (printOut_ > 0) std::cout << " DeDX esimtator  "<< dEdxTrack[track].dEdx()
                 << "  " << dEdxTrack[track].dEdxError() 
                 << "  " << dEdxTrack[track].numberOfSaturatedMeasurements()
                 << "  " << dEdxTrack[track].numberOfMeasurements() << std::endl;
         // existe aussi un dedxPixelHarmonic2 et dedxTruncated40 qui sont stockés
         tree_track_dedx_harmonic2[tree_ntracks]= dEdxTrack[track].dEdx();
       }
       else tree_track_dedx_harmonic2[tree_ntracks]= -10;
       */

       //hit level dEdx information (only done for MIPs)

       if (printOut_ > 0) std::cout << " with " << dedxHits->size() << " dedxHits info " << std::endl;

       tree_track_nhits[tree_ntracks]= dedxHits->size();
       tree_track_index_hit[tree_ntracks]=tree_dedxhits;

       // values for 2016 MC --> which values for 2017???
       //double dEdxSF [2] = { 1.09711, 1.09256 };  // 0 : Strip SF, 1 : Pixel to Strip SF
       // new SF Profile
//       double dEdxSF [2] = { 1., 1.2025 };  // 0 : Strip SF, 1 : Pixel to Strip SF
//       double dEdxSF_corr [2] = { 1., 1.5583 };  // 0 : Strip SF, 1 : Pixel to Strip SF
       double dEdxSF [2] = { 1., 1. };  // 0 : Strip SF, 1 : Pixel to Strip SF
       double dEdxSF_corr [2] = { 1., 1. };  // 0 : Strip SF, 1 : Pixel to Strip SF

       // Hits, SF, Templates, usePixel, useClusterCleaning, reverseProb, useTrunc,  TrackerGains, useStrips, mustBeInside, MaxStripNOM, 
       // correctFEDSat, CrossTalkInv, dropLowerValue, ErrorDeDx
       tree_track_ih_ampl[tree_ntracks]=0;
       // correction inverseXtalk = 0 --> take the raw amplitudes of the cluster;
       reco::DeDxData dedxMObjTmp = computedEdx(dedxHits, dEdxSF, NULL, true, true, false , false, NULL, true, true, 99, false, 0, 0.15,  NULL);
       reco::DeDxData* dedxMObj = dedxMObjTmp.numberOfMeasurements()>0?&dedxMObjTmp:NULL;
       if(dedxMObj) tree_track_ih_ampl[tree_ntracks]=dedxMObj->dEdx();

       tree_track_ih_ampl_corr[tree_ntracks]=0;
       // correction inverseXtalk = 1 --> modify the amplitudes based on xtalk for non-sat cluster + correct for saturation
       reco::DeDxData dedxMObjTmp2 = computedEdx(dedxHits, dEdxSF_corr, NULL, true, true, false , false, NULL, true, true, 99, false, 1, 0.15,  NULL);
       reco::DeDxData* dedxMObj2 = dedxMObjTmp2.numberOfMeasurements()>0?&dedxMObjTmp2:NULL;
       if(dedxMObj2) tree_track_ih_ampl_corr[tree_ntracks]=dedxMObj2->dEdx();

       tree_track_ias_ampl[tree_ntracks]=0;
       tree_track_ias_ampl_corr[tree_ntracks]=0;
       if (dEdxTemplatesUncorr) {
       // correction inverseXtalk = 0 --> take the raw amplitudes of the cluster;
       reco::DeDxData dedxSObjTmp = computedEdx(dedxHits, dEdxSF, dEdxTemplatesUncorr, true, true, false , false, NULL, true, true, 99, false, 0, 0.0,  NULL);
       reco::DeDxData* dedxSObj = dedxSObjTmp.numberOfMeasurements()>0?&dedxSObjTmp:NULL;
       if(dedxSObj) tree_track_ias_ampl[tree_ntracks]=dedxSObj->dEdx();
       }
       if (dEdxTemplatesCorr) {
       // correction inverseXtalk = 1 --> modify the amplitudes based on xtalk for non-sat cluster + correct for saturation
       reco::DeDxData dedxSObjTmp2 = computedEdx(dedxHits, dEdxSF_corr, dEdxTemplatesCorr, true, true, false , false, NULL, true, true, 99, false, 1, 0.0,  NULL);
       reco::DeDxData* dedxSObj2 = dedxSObjTmp2.numberOfMeasurements()>0?&dedxSObjTmp2:NULL;
       if(dedxSObj2) tree_track_ias_ampl_corr[tree_ntracks]=dedxSObj2->dEdx();
       }

    

       for(unsigned int h=0;h< dedxHits->size();h++){
          if (tree_dedxhits<nMaxDeDxH) {
             DetId detid(dedxHits->detId(h));
             if (printOut_ > 0) std::cout << " DetId " << (uint32_t) detid << std::endl;
             tree_dedx_detid[tree_dedxhits]=(uint32_t) detid;
             tree_dedx_subdetid[tree_dedxhits]=detid.subdetId();
             if (detid.subdetId()<3) tree_dedx_modulgeom[tree_dedxhits]=15;
             else {
               SiStripDetId SSdetId(detid);
               tree_dedx_modulgeom[tree_dedxhits]=SSdetId.moduleGeometry();
             }
             tree_dedx_insideTkMod[tree_dedxhits]=isHitInsideTkModule(dedxHits->pos(h), detid, detid.subdetId()>=3?dedxHits->stripCluster(h):NULL);
             tree_dedx_charge[tree_dedxhits]=dedxHits->charge(h);
             tree_dedx_pathlength[tree_dedxhits]=dedxHits->pathlength(h);
             const LocalPoint position_hit=dedxHits->pos(h);
             tree_dedx_posx[tree_dedxhits]=position_hit.x();
             tree_dedx_posy[tree_dedxhits]=position_hit.y();
             tree_dedx_posz[tree_dedxhits]=position_hit.z();

             if (printOut_ > 0) std::cout << " charge of " << h << " = " << dedxHits->charge(h) << std::endl;     
             if (printOut_ > 0) std::cout << " path length = " << dedxHits->pathlength(h) << std::endl;
             if (printOut_ > 0) std::cout << " ChargeOverPathlength " << dedxHits->charge(h)/dedxHits->pathlength(h) << std::endl;
             if (printOut_ > 0) std::cout << " DetId " << detid.subdetId() << std::endl;

             if (detid.subdetId()>=3) {
               tree_dedx_strip[tree_dedxhits]=true;
               tree_dedx_pixel[tree_dedxhits]=false;
 
               tree_sclus_firstsclus[tree_dedxhits]=(int) dedxHits->stripCluster(h)->firstStrip();
               tree_sclus_barycenter[tree_dedxhits]=dedxHits->stripCluster(h)->barycenter();
               tree_sclus_charge[tree_dedxhits]=dedxHits->stripCluster(h)->charge();
               tree_sclus_errorclus[tree_dedxhits]=dedxHits->stripCluster(h)->getSplitClusterError();
               tree_sclus_ismerged[tree_dedxhits]=dedxHits->stripCluster(h)->isMerged();
               tree_sclus_index_strip[tree_dedxhits]=tree_nstrips;
               tree_sclus_index_strip_corr[tree_dedxhits]=tree_nstrips_corr;
               tree_sclus_charge_corr[tree_dedxhits]=0;
               tree_sclus_sat254[tree_dedxhits]=false;
               tree_sclus_sat255[tree_dedxhits]=false;
               tree_sclus_shape[tree_dedxhits]=DeDxTools::shapeSelection(*(dedxHits->stripCluster(h)));

               std::vector <uint8_t> amplis = dedxHits->stripCluster(h)->amplitudes();
               std::vector <int> amps = convert(amplis);
               if (printOut_ > 0) std::cout << " amps.size() "<< amps.size() << std::endl;


               tree_sclus_nstrip[tree_dedxhits]=amps.size();
               for (unsigned int iclu=0;iclu<amps.size();iclu++) {

                 if (printOut_ > 0) std::cout << "   amplitude (cluster" << iclu <<") = " <<  amps[iclu] <<std::endl;
                 if ( amps[iclu]>=254) tree_sclus_sat254[tree_dedxhits]=true;
                 if ( amps[iclu]==255) tree_sclus_sat255[tree_dedxhits]=true;
                 if (tree_nstrips< nMaxStrip) {
                  tree_strip_ampl[tree_nstrips]= amps[iclu];
                  tree_nstrips++;
                 } // end if MaxStrips
                 else {
                   std::cout << "Limit reached for tree_nstrips "<< tree_nstrips << std::endl;
                 } // end if MaxStrips
               } // end loop strips
               
               std::vector <int> ampsprim = CrossTalkInv(amps, 0.10, 0.04, true);
               tree_sclus_nstrip_corr[tree_dedxhits]=ampsprim.size();
               for (unsigned int iclu=0;iclu<ampsprim.size();iclu++) {
                 if (printOut_ > 0) std::cout << " corrected amplitude (cluster" << iclu <<") = " << (int) ampsprim[iclu] <<std::endl;
                 tree_sclus_charge_corr[tree_dedxhits]+=ampsprim[iclu];
                 if (tree_nstrips_corr< nMaxStripprim) {
                  tree_strip_ampl_corr[tree_nstrips_corr]= ampsprim[iclu];
                  tree_nstrips_corr++;
                 } // end if nMaxStripprim
                 else {
                   std::cout << "Limit reached for tree_nstrips_corr "<< tree_nstrips_corr << std::endl;
                 } // end if nMaxStripprim
               } // end loop strips

               tree_sclus_clusclean[tree_dedxhits]= clusterCleaning(amps, 0);
               tree_sclus_clusclean2[tree_dedxhits]= clusterCleaning(ampsprim, 1);



               // search for the associated simHit
               if (m_runOnGS) {
//                edm::DetSetVector<StripDigiSimLink>::const_iterator isearch = stripdigisimlink->find(detid);
                edm::DetSetVector<StripDigiSimLink>::const_iterator isearch = stripdigisimlink->find(dedxHits->detId(h));
                
                // Look for a digisimlink matching this cluster
                if(isearch != stripdigisimlink->end()) {
                  edm::DetSet<StripDigiSimLink> link_detset = (*isearch);
                  float clusEloss = 0;
                  std::vector<simhitAddr> CFaddr;
                  std::vector<unsigned int> hitProcess;
                  std::vector<int> hitPID;
                  std::vector<float> trackCharge;
                  std::vector<float> hitPmag;
                  std::vector<float> hitPathLength;
                  tree_sclus_index_simhit[tree_dedxhits]=tree_nsimhits;
                  tree_sclus_nsimhit[tree_dedxhits]=0;
                  tree_sclus_eloss[tree_dedxhits]=0;
                  for(edm::DetSet<StripDigiSimLink>::const_iterator linkiter = link_detset.data.begin(), linkEnd = link_detset.data.end();
                      linkiter != linkEnd; ++linkiter) {
                      int theChannel = linkiter->channel();
                      if( theChannel >= tree_sclus_firstsclus[tree_dedxhits]  && theChannel < tree_sclus_firstsclus[tree_dedxhits]+tree_sclus_nstrip[tree_dedxhits] ) { 

                        bool booly=false;
                        if ((int)detid==369120293 && tree_sclus_firstsclus[tree_dedxhits] ==73) booly=true;
                        if ((int)detid==369120294 && tree_sclus_firstsclus[tree_dedxhits] ==689) booly=true;
                        if ((int)detid==369120313 && tree_sclus_firstsclus[tree_dedxhits] ==505) booly=true;
                        if ((int)detid==369120314 && tree_sclus_firstsclus[tree_dedxhits] ==272) booly=true;
                        if ((int)detid==369120333 && tree_sclus_firstsclus[tree_dedxhits] ==95) booly=true;
                        if (booly) {
                            std::cout << "Stage  detID " << (int)detid << " pos "  << linkiter->CFposition() << std::endl;
                         }


                         // This digisimlink points to a strip in the current cluster
                         int stripIdx = theChannel - tree_sclus_firstsclus[tree_dedxhits] ;
                         if (printOut_ > 0) std::cout << "association cluster-simhit" << std::endl;
                         if (printOut_ > 0) std::cout <<  "channel = " << linkiter->channel() << " TrackID = " << linkiter->SimTrackId() 
                                   <<  " EventID = " <<  linkiter->eventId().rawId() << " TofBin = " <<  linkiter->TofBin()
                                   <<  " CFPos = " << linkiter->CFposition() << " fraction = " << linkiter->fraction()
                                   <<  " stripIdx = " << stripIdx << std::endl;
                                   // << " amp = " << amps[stripIdx] << std::endl;
                         unsigned int currentCFPos = linkiter->CFposition();
                         unsigned int tofBin = linkiter->TofBin();
                         simHitCollectionID theSimHitCollID = std::make_pair(detid.subdetId(), tofBin);
                         simhitAddr currentAddr = std::make_pair(theSimHitCollID, currentCFPos);
                         bool newHit = true;
                         if (std::find(CFaddr.begin(), CFaddr.end(), currentAddr) != CFaddr.end()) newHit = false;
                         if (newHit) {
                           simhit_collectionMap::const_iterator it = SimHitCollMap_.find(theSimHitCollID);
                           if (it!= SimHitCollMap_.end()) {
                              if (currentCFPos < (it->second).size()) {
                                const PSimHit& theSimHit = (it->second)[currentCFPos];
                                CFaddr.push_back(currentAddr);
                                hitProcess.push_back(theSimHit.processType());
                                hitPID.push_back(theSimHit.particleType());
                                hitPmag.push_back(theSimHit.pabs());
                                Local3DPoint entry = theSimHit.entryPoint();
                                Local3DPoint exit = theSimHit.exitPoint();
                                Local3DVector segment = exit - entry;
                                hitPathLength.push_back(segment.mag());
                                clusEloss += theSimHit.energyLoss();  // Add up the contributions of all simHits to this cluster
                                if (printOut_ > 0) std::cout << "SimHit " <<  int(CFaddr.size()) << ", process = " <<  theSimHit.processType() 
                                          << ", PID = " << theSimHit.particleType() << ", p = " << theSimHit.pabs() << ", Eloss = " << theSimHit.energyLoss()
                                          << ", segment = " << segment.mag() << std::endl;
//                                          << ", segment = " << segment.mag() << ", str segment = " << modPathLength << std::endl;
                                if (booly) std::cout << "Stage      "  << stripIdx << " " << entry.x() << "  " << entry.y() << "  "  << theSimHit.energyLoss() << std::endl;
                                tree_sclus_nsimhit[tree_dedxhits]+=1;
                                if (tree_nsimhits< nMaxSimHit) {
                                    tree_simhit_pid[tree_nsimhits]=theSimHit.particleType();
                                    tree_simhit_process[tree_nsimhits]=theSimHit.processType();
                                    tree_simhit_p[tree_nsimhits]=theSimHit.pabs();
                                    tree_simhit_eloss[tree_nsimhits]=theSimHit.energyLoss();
                                    tree_sclus_eloss[tree_dedxhits]+=theSimHit.energyLoss();
                                    tree_simhit_tof[tree_dedxhits]=theSimHit.timeOfFlight();
                                    tree_simhit_segment[tree_nsimhits]=segment.mag();
                                    tree_simhit_xentry[tree_nsimhits]=entry.x();
                                    tree_simhit_yentry[tree_nsimhits]=entry.y();
                                    tree_simhit_zentry[tree_nsimhits]=entry.z();
                                    tree_simhit_xexit[tree_nsimhits]=exit.x();
                                    tree_simhit_yexit[tree_nsimhits]=exit.y();
                                    tree_simhit_zexit[tree_nsimhits]=exit.z();
                                    tree_nsimhits++;
                                } // end if MaxSimHits
                                else {
                                    std::cout << "Limit reached for tree_nsimhits "<< tree_nsimhits << std::endl;
                                } // end if MaxSimHits
                              } else {
                                   if (printOut_ > 0) std::cout << "currentCFPos " << currentCFPos << " is out of range for " << (it->second).size() << std::endl;
                              }
                          }
                        }  // if (newHit)

                      } // if theChannel
                   } // end for linkiter
                   if (tree_sclus_nsimhit[tree_dedxhits]==0) tree_sclus_index_simhit[tree_dedxhits]=-1;
                } // if isearch
                else {
                 tree_sclus_index_simhit[tree_dedxhits]=-1;
                 tree_sclus_nsimhit[tree_dedxhits]=0;
                 tree_sclus_eloss[tree_dedxhits]=0;
                }
               } // if runOnGS
               else {
                 tree_sclus_index_simhit[tree_dedxhits]=-1;
                 tree_sclus_nsimhit[tree_dedxhits]=0;
                 tree_sclus_eloss[tree_dedxhits]=0;
               }
               // end search for associated simHit
                   //
             }
             else {
               tree_dedx_strip[tree_dedxhits]=false;
               tree_dedx_pixel[tree_dedxhits]=true;
               tree_sclus_firstsclus[tree_dedxhits]=-1;
               tree_sclus_barycenter[tree_dedxhits]=-1;
               tree_sclus_charge[tree_dedxhits]=-1;
               tree_sclus_errorclus[tree_dedxhits]=-1;
               tree_sclus_ismerged[tree_dedxhits]=0;
               tree_sclus_index_strip[tree_dedxhits]=-1;
               tree_sclus_nstrip[tree_dedxhits]=0;
               tree_sclus_index_simhit[tree_dedxhits]=-1;
               tree_sclus_nsimhit[tree_dedxhits]=0;
               tree_sclus_eloss[tree_dedxhits]=0;
               tree_sclus_sat254[tree_dedxhits]=0;
               tree_sclus_sat255[tree_dedxhits]=0;
               tree_sclus_shape[tree_dedxhits]=0;
               tree_sclus_clusclean[tree_dedxhits]=0;
               tree_sclus_charge_corr[tree_dedxhits]=-1;
               tree_sclus_nstrip_corr[tree_dedxhits]=0;
               tree_sclus_clusclean2[tree_dedxhits]=0;

             } // end if detID
             tree_dedxhits++;
             
          } //end if MaxDeDx
          else {
                   std::cout << "Limit reached for tree_dedxhits "<< tree_dedxhits << std::endl;
          } // end if MaxDeDx
        } // end loop dEdx
        tree_ntracks++;
      } // end if MaxTracks 
      else {
                   std::cout << "Limit reached for tree_ntracks "<< tree_ntracks << std::endl;
      } // end if MaxTracks




    } // end loop TrackCollection


    std::vector<reco::MuonRef> KeepMuonRefVec;
    if (m_format != "miniAOD" )  {
    for(unsigned int c=0;c<muonCollectionHandle->size();c++){
      reco::MuonRef muon  = reco::MuonRef( muonCollectionHandle, c );

      if (tree_nmuons < nMaxMuon) {
       KeepMuonRefVec.push_back(muon);
       if (printOut_ > 0) std::cout << " muon with pT =  " << muon->pt() << std::endl;
       tree_muon_pt[tree_nmuons]= muon->pt();
       if (muon->isStandAloneMuon()) {
             tree_muon_ptSA[tree_nmuons] = muon->standAloneMuon()->pt();
       }
       else {
            tree_muon_ptSA[tree_nmuons]= -1;
       }
       if (muon->innerTrack().isNonnull()) tree_muon_ptIT[tree_nmuons]= muon->innerTrack()->pt();
       else tree_muon_ptIT[tree_nmuons]= -1;
       tree_muon_p[tree_nmuons]= muon->p();
       tree_muon_eta[tree_nmuons]= muon->eta();
       tree_muon_phi[tree_nmuons]= muon->phi();

       const reco::MuonTimeExtra* tof = NULL;
       const reco::MuonTimeExtra* dttof = NULL;
       const reco::MuonTimeExtra* csctof = NULL;

       if (TOFCollH.isValid()) {
       tof  = &TOFCollH->get(muon.key()); 
       tree_muon_comb_inversebeta[tree_nmuons] = tof->inverseBeta();
       tree_muon_comb_inversebetaerr[tree_nmuons] = tof->inverseBetaErr();
       tree_muon_comb_tofndof[tree_nmuons] = tof->nDof();
       tree_muon_comb_vertextime[tree_nmuons] = tof->timeAtIpInOut();
       }
       else {
       tree_muon_comb_inversebeta[tree_nmuons] = -10;
       tree_muon_comb_inversebetaerr[tree_nmuons] = -10;
       tree_muon_comb_tofndof[tree_nmuons] = 0;
       tree_muon_comb_vertextime[tree_nmuons] = -10;
       }
       if (TOFDTCollH.isValid()) {
       dttof = &TOFDTCollH->get(muon.key());  
       tree_muon_dt_inversebeta[tree_nmuons] = dttof->inverseBeta();
       tree_muon_dt_inversebetaerr[tree_nmuons] = dttof->inverseBetaErr();
       tree_muon_dt_tofndof[tree_nmuons] = dttof->nDof();
       tree_muon_dt_vertextime[tree_nmuons] = dttof->timeAtIpInOut();
       }
       else {
       tree_muon_dt_inversebeta[tree_nmuons] = -10;
       tree_muon_dt_inversebetaerr[tree_nmuons] = -10;
       tree_muon_dt_tofndof[tree_nmuons] = 0;
       tree_muon_dt_vertextime[tree_nmuons] = -10;
       }
       if (TOFCSCCollH.isValid()) {
       csctof = &TOFCSCCollH->get(muon.key());
       tree_muon_csc_inversebeta[tree_nmuons] = csctof->inverseBeta();
       tree_muon_csc_inversebetaerr[tree_nmuons] = csctof->inverseBetaErr();
       tree_muon_csc_tofndof[tree_nmuons] = csctof->nDof();
       tree_muon_csc_vertextime[tree_nmuons] = csctof->timeAtIpInOut();
       }
       else {
       tree_muon_csc_inversebeta[tree_nmuons] = -10;
       tree_muon_csc_inversebetaerr[tree_nmuons] = -10;
       tree_muon_csc_tofndof[tree_nmuons] = 0;
       tree_muon_csc_vertextime[tree_nmuons] = -10;
       }
       tree_nmuons++;

      } // end if MaxMuon
    } // end loop MuonCollectoon
    } // end if AOD


     //loop on HSCP candidates
     tree_hscp=0;
     if (m_format != "miniAOD" && hscpCollH.isValid()) {

     const susybsm::HSCParticleCollection& hscpColl = *hscpCollH;
     for(unsigned int c=0;c<hscpColl.size();c++){
       //define alias for important variable
       susybsm::HSCParticle hscp  = hscpColl[c];
       reco::MuonRef  muon  = hscp.muonRef();
 
       //For TOF only analysis use updated stand alone muon track.
       //Otherwise use inner tracker track      
       reco::TrackRef track = hscp.trackRef();

       if (printOut_ > 0) std::cout << "HSCP candidate " << c << std::endl;
       if (track.isNull() && printOut_ > 0 ) std::cout << "probleme track isNull"<< std::endl;
       if (muon.isNull() && printOut_ > 0 ) std::cout << "probleme muon isNull"<< std::endl;

       if (tree_hscp<nMaxHSCP) {



       if (n_genp>0 && !track.isNull()) {
        double RMin = 9999;
        int idxG=-1;
        for(int i=0;i< n_genp ;++i){
         const reco::GenParticle* genCand = &(*GenColl)[i];
         if(genCand->pt()<5) continue;
         if(genCand->status()!=1)continue;
         int AbsPdg=abs(genCand->pdgId());
         if(AbsPdg<1000000 && AbsPdg!=17)continue;
         double dR = deltaR(track->eta(), track->phi(), genCand->eta(), genCand->phi());
         if(dR<RMin){ 
           RMin=dR;
           idxG=genCand->pdgId();
         }
        }
        tree_hscp_gen_id[tree_hscp]=idxG;
        tree_hscp_gen_dr[tree_hscp]=RMin;
       }
       else {
        tree_hscp_gen_id[tree_hscp]=-1;
        tree_hscp_gen_dr[tree_hscp]=-1;
       }

       float Mindelta=0.001;
       int tr_index=-1;
       if (KeepTrackRefVec.size()>0 && !track.isNull()) {
         if (printOut_ > 0 ) std::cout << "KeepTrackRefVec.size()" << KeepTrackRefVec.size() << std::endl;
         for(unsigned int tr=0;tr<KeepTrackRefVec.size();tr++){
            if( fabs( (1.0/KeepTrackRefVec[tr]->pt())-(1.0/track->pt())) < Mindelta) {
                  Mindelta=fabs( (1.0/KeepTrackRefVec[tr]->pt())-(1.0/track->pt()));
                  tr_index=tr;
            } 
         }
         if (printOut_ > 0  && tr_index>-1)   std::cout << "association HSCP - track " << KeepTrackRefVec[tr_index]->pt() << "  "  << track->pt() 
                          << "   " << KeepTrackRefVec[tr_index]->eta() << "  "  << track->eta() 
                          << "   " << KeepTrackRefVec[tr_index]->phi() << "  "  << track->phi() << std::endl;
       }
       tree_hscp_track_idx[tree_hscp]=tr_index;

       float Mindelta2=0.001;
       int mu_index=-1;
       if (KeepMuonRefVec.size()>0 && !muon.isNull()) {
         if (printOut_ > 0 ) std::cout << "KeepMuonRefVec.size()" << KeepMuonRefVec.size() << std::endl;
         for(unsigned int tr=0;tr<KeepMuonRefVec.size();tr++){
            if( fabs( (1.0/KeepMuonRefVec[tr]->pt())-(1.0/muon->pt())) < Mindelta2) {
                  Mindelta2=fabs( (1.0/KeepMuonRefVec[tr]->pt())-(1.0/muon->pt()));
                  mu_index=tr;
            } 
         }
         if (printOut_ > 0  && mu_index>-1  )   std::cout << "association HSCP - muon " << KeepMuonRefVec[mu_index]->pt() << "  "  << muon->pt() 
                          << "   " << KeepMuonRefVec[mu_index]->eta() << "  "  << muon->eta() 
                          << "   " << KeepMuonRefVec[mu_index]->phi() << "  "  << muon->phi() << std::endl;
       }
       tree_hscp_muon_idx[tree_hscp]=mu_index;


       if (!track.isNull() && isoHSCPColl.isValid()) {
           const ValueMap<susybsm::HSCPIsolation>& IsolationMap = *isoHSCPColl.product();
           susybsm::HSCPIsolation hscpIso = IsolationMap.get((size_t)track.key());   
           if (printOut_ > 0) std::cout << "HSCP iso " << hscpIso.Get_TK_SumEt() << "  "  << hscpIso.Get_ECAL_Energy() 
                                        << "     "  << hscpIso.Get_HCAL_Energy() << "  for HSCP p " << track->p()  << std::endl;
           tree_hscp_iso_tk[tree_hscp]=hscpIso.Get_TK_SumEt();
           tree_hscp_iso_ecal[tree_hscp]=hscpIso.Get_ECAL_Energy();
           tree_hscp_iso_hcal[tree_hscp]=hscpIso.Get_HCAL_Energy();
       }
       else {
           tree_hscp_iso_tk[tree_hscp]=-10;
           tree_hscp_iso_ecal[tree_hscp]=-10;
           tree_hscp_iso_hcal[tree_hscp]=-10;
       }

//       if(TypeMode!=3) track = hscp.trackRef();
//       else {
//         if(muon.isNull()) continue;
//         track = muon->standAloneMuon();
//       }
       //skip events without track
//       if(track.isNull())continue;
       
       //require a track segment in the muon system
       //if(TypeMode>1 && TypeMode!=5 && (muon.isNull() || !muon->isStandAloneMuon()))continue;
      
       //Apply a scale factor to muon only analysis to account for differences seen in data/MC preselection efficiency
       //For eta regions where Data > MC no correction to be conservative
       //if(!isData && TypeMode==3 && scaleFactor(track->eta())<RNG->Uniform(0, 1)) continue;
     
        tree_hscp++;
       }

     } // end loop HSCP 

    } // end if HSCP collection valid and !miniAOD






    smalltree -> Fill();

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
ntuple::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ntuple::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ntuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ntuple);
