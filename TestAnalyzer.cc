// -*- C++ -*-
//
// Package:    ana/TestAnalyzer
// Class:      TestAnalyzer
// 
/**\class TestAnalyzer TestAnalyzer.cc ana/TestAnalyzer/plugins/TestAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yuxiang Guo
//         Created:  Tue, 12 May 2015 13:17:39 GMT
//
//


// system include files
#include <memory>
// user include files
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <iostream>
#include <cmath>
#include<TCanvas.h>
using namespace edm;
using namespace reco;
using namespace std;
//
// class declaration
//

class TestAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TestAnalyzer(const edm::ParameterSet&);
      ~TestAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
	TH1F *histWpt;
	TH1F *histZpt;
	TH1F *histZeta;
	TH1F *histWeta;
	TH1F *histWphi;
	TH1F *histempt;
	TH1F *histemeta;
	TH1F *histemphi;
	TH1F *histmupt;
	TH1F *histmuallpt;
	TH1F *histmualleta;
	TH1F *histmueta;
	TH1F *histmuphi;
	TH1F *histmugenpt;
	TH1F *histmugeneta;
	TH1F *histmugenphi;
	TH1F *histtaupt;
	TH1F *histtaueta;
	TH1F *histtauphi;
	TH1F *histelimiate;
	TH1F *histZmass;
	TH1F *histDmass;
	TH1F *histinZmass;
	TH1F *histDpt;
	TH1F *histZchild;
	TH1F *histtauchild;  
	TH1F *histtaumass;   
	TH1F *histtaumchild; 
	TH1F *histgmass; 
	TH1F *hist4mmass; 
	TH1F *hist4mpt;
	TH1F *histreco4mmass; 
	TCanvas *ptcomparison;
	double tau_px;
	double tau_py;
	double tau_pz;
	double tau_energy;
	double tau_pt;
	double tau_eta;
	double tau_phi;
	double tau_deltaR;
	double in_z;
	double px_m;         //the variable in order to calculate the invariant mass
	double py_m;
	double pz_m;
	double e_p[50];
	double e_m;	     //e for energy
	double pt_p[50];	//p is on generator level m is on muon level
	double pt_m[50];  
	double gen_muon_eta[50]; // the selected 3 muons
	double gen_muon_pt[50];
	double gen_muon_phi[50];
	double reco_pt[50];
	double reco_eta[50];
	double reco_phi[50];
	TH1F *histmutauptratio;
	TH1F *histmutaudeltaR;
	TH1F *histmutaupt;
	TH1F *histmuongenallpt;
	TH1F *histmuongenalleta;
	TH1F *histmumatchpt;
	TH1F *histmulastpt;
	TH1F *histmulasteta;
	TH1F *histmumatcheta;
	TH1F *histmuongencutpt;
	TH1F *histmuongencuteta;
	TH1F *histmaxptgen;
	TH1F *histminptgen;
	TH1F *histmaxptreco;
	TH1F *histminptreco;
	TH1F *histrecomuongroupeta;
	TH1F *histgenmuongroupeta;
  	TH1F *histrecomuongrouppt;
	TH1F *histrecomuongroupphi;
  	TH1F *histgenmuongrouppt;
  	TH1F *histgenmuongroupphi;
  	TH1F *histcount;
	TH1F *histgenaccept;
	TH1F *histrecoaccept;
	TH1F *histgenFpt;
	TH1F *histgenSpt;
	TH1F *histgenTpt;
	TH1F *histgenFeta;
	TH1F *histgenSeta;
	TH1F *histgenTeta;
	TH1F *histgenmaxdeltaR;
	TH1F *histgenmindeltaR;
	TH1F *histrecomaxISO;
	TH1F *histgenmaxrecodeltaR;
	TH1F *histgenminrecodeltaR;
	TH1F *histrecoFeta;
	TH1F *histrecoSeta;
	TH1F *histrecoTeta;
	TH1F *histrecoFpt;
	TH1F *histrecoSpt;
	TH1F *histrecoTpt;
	TH1F *histgenrecodeltaRmax;
	TH1F *histgenrecodeltaRmin;
	TH1F *histrecodeltaRaccept;
	TH1F *histrecomatchingaccept;
	TH1F *histrecopTaccept;
	TH1F *histrecoetaaccept;
	TH1F *histrecopfaccept;
	TH1F *histgenothermupttau;
	TH1F *histgenothermupthar;
	TH1F *histgenothermudeltaRmax;
	TH1F *histgenothermudeltaRmin;
	TH1F *histgenothermunumber;
	TH1F *histgenothermumother;	
	TH1F *histMETpt;	
	TH1F *histMETinvariantmass;	
	TH1F *histMETphi;	
			//virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
			//
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
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
TestAnalyzer::TestAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
	Service<TFileService> fs;
	histmaxptreco=fs->make<TH1F>("leading_pt_reco","leading_muon_pt_reco",100,0,100);
	histmaxptgen=fs->make<TH1F>("leading_pt_gen","leading_muon_pt_gen",100,0,100);
	histminptgen=fs->make<TH1F>("least_pt_gen","leading_muon_pt_gen",100,0,100);
	histminptreco=fs->make<TH1F>("least_pt_reco","leading_muon_pt_reco",100,0,100);
	histZpt=fs->make<TH1F>("Z_pt_gen","Z pt",100,0,100);
	histZeta=fs->make<TH1F>("Z_eta_gen","Z eta",100,-15,15);
	histWpt=fs->make<TH1F>("Wpt_gen","W+ pt",10,0,10);
	histWeta=fs->make<TH1F>("Weta_gen","W+ eta",100,-15,15);
	histWphi=fs->make<TH1F>("Wphi_gen","W+ phi",100,-10,10);
	histempt=fs->make<TH1F>("empt_gen","e+ pt",100,10,100);
	histemeta=fs->make<TH1F>("emeta_gen","e+ eta",100,-10,10);
	histemphi=fs->make<TH1F>("emphi_gen","e+ phi",100,-10,10);
	histmupt=fs->make<TH1F>("mupt_reco_from_tau","mu pt_from_tau_matching",100,0,100);
	histmutauptratio=fs->make<TH1F>("mu_tau_pt_ratio","mu_tau_pt_ratio",100,0,1);
	histmutaudeltaR=fs->make<TH1F>("mu_tau_deltaR","mu_tau_pt_deltaR",100,0,6);
	histmutaupt=fs->make<TH1F>("mu_tau_combinedpT","mu_tau_combinedpT",100,0,100);
	histmupt=fs->make<TH1F>("mupt_reco_from_tau","mu pt_from_tau_matching",100,0,100);
	histmupt=fs->make<TH1F>("mupt_reco_from_tau","mu pt_from_tau_matching",100,0,100);
     	histmueta=fs->make<TH1F>("mueta_reco","mu eta_from_tau_matching",40,-10,10);
	histmuphi=fs->make<TH1F>("muphi_reco","mu phi_from_tau_matching",100,-10,10);
	histmuallpt=fs->make<TH1F>("mupt_reco_pt","mu_reco_all_pt",100,0,100);
	histmualleta=fs->make<TH1F>("mupt_reco_eta","mu_reco_all_eta",100,-10,10);
	histmumatchpt=fs->make<TH1F>("mupt_reco_matching_pt","mu_reco_matching_pt",100,0,100);
	histmumatcheta=fs->make<TH1F>("mupt_reco_matching_eta","mu_reco_matching_eta",100,-3,3);
	histmulastpt=fs->make<TH1F>("mupt_reco_matching&ISO_pt","mu_reco_matching&ISO&PF_pt",100,0,100);
	histmulasteta=fs->make<TH1F>("mupt_reco_matching&ISO_eta","mu_reco_matching&ISO&PF_eta",100,-3,3);
	histmugenpt=fs->make<TH1F>("mupt_gen_from_tau","mu_gen_from_tau_pt",100,0,100);
	histmugeneta=fs->make<TH1F>("mueta_gen_from_tau","mu_gen_from_tau_eta",100,-10,10);
	histmugenphi=fs->make<TH1F>("muphi_gen_from_tau","mu_gen_from_tau_ phi",100,-10,10);
	histtaupt=fs->make<TH1F>("tau_gen_pt","tau_gen_pt",200,0,200);
	histtaueta=fs->make<TH1F>("tau_gen_eta","tau_gen_eta",100,-15,15);
	histtauphi=fs->make<TH1F>("tau_gen_phi","tau_gen_phi",100,-10,10);
	histelimiate=fs->make<TH1F>("elmiate","elimiate",10,0,10);
	histZmass=fs->make<TH1F>("Z_mass_gen","Z_mass_gen",100,20,160);
	histDmass=fs->make<TH1F>("invariance_mass_reco","invariance_mass_reco",100,1.65,1.85);
	histinZmass=fs->make<TH1F>("invariance_mass_Z","invariance_mass_Z",100,0,100);
	histDpt=fs->make<TH1F>("Decay_pt","Decay_pt",100,0,100);
	histZchild=fs->make<TH1F>("child_of_Z_gen","Z_child",60,-30,30);	
	histtauchild=fs->make<TH1F>("child_of_tau_gen","child_of_tau_13_gen",60,-30,30);
	histtaumass=fs->make<TH1F>("tau_mass_gen","tau mass_gen",100,1.7,1.8);
	histtaumchild=fs->make<TH1F>("tau+_child_gen","tau+_child_gen",60,-30,30);
	histgmass=fs->make<TH1F>("invariance_mass_gen","invariance_mass_gen",100,1.5,2);
	hist4mmass=fs->make<TH1F>("invariance_mass_4m","invariance_mass_4m",100,0,100);
	hist4mpt=fs->make<TH1F>("4th_mu_pt","4th_mu_pt",100,0,100);
	histreco4mmass=fs->make<TH1F>("invariance_mass_reco_4m","invariance_mass_reco_4m",100,0,100);
	ptcomparison=fs->make<TCanvas>("Pt_Comparision","Pt_Comparision",800,1000);
	histmuongencutpt=fs->make<TH1F>("muon_gen_pt,eta,cut_pt","muon_gen_pt,eta_cut_pt",100,0,100);
	histmuongencuteta=fs->make<TH1F>("muon_gen_pt,eta,cut,eta","muon_gen_pt,eta_cut,eta",100,-5,5);
	histmuongenallpt=fs->make<TH1F>("hist_muon_gen_all_pt","hist_muon_gen_all_pt",100,0,100);
	histmuongenalleta=fs->make<TH1F>("hist_muon_gen_all_eta","hist_muon_gen_all_eta",100,-10,10);
	histrecomuongrouppt=fs->make<TH1F>("3_muon_reco_pt","3_muon_reco_pt",100,100,100);
	histrecomuongroupeta=fs->make<TH1F>("3_muon_reco_eta","3_muon_reco_eta",100,-5,5);
	histrecomuongroupphi=fs->make<TH1F>("3_muon_reco_phi","3_muon_reco_phi",100,-10,10);
	histgenmuongrouppt=fs->make<TH1F>("3_muon_gen_pt","3_muon_gen_pt",100,100,100);
	histgenmuongroupeta=fs->make<TH1F>("3_muon_gen_eta","3_muon_gen_eta",100,-5,5);
	histgenmuongroupphi=fs->make<TH1F>("3_muon_gen_phi","3_muon_gen_phi",100,-10,10);
	histcount=fs->make<TH1F>("no.muon.each.event","no.muon.each.event",100,-10,10);
	histgenaccept=fs->make<TH1F>("No_of_accpeted_muon_gen","number of accept muon_gen",4,0,4);
	histrecoaccept=fs->make<TH1F>("No_of_accpeted_muon_reco","number of accept muon_reco",10,0,10);
	histgenFpt=fs->make<TH1F>("Gen_Muon_least_pt","Gen_Muon_least_Pt",100,0,100);
	histgenSpt=fs->make<TH1F>("Gen_Muon_mid_pt","Gen_Muon_mid_Pt",100,0,100);
	histgenTpt=fs->make<TH1F>("Gen_Muon_leading_pt","Gen_Muon_leading_Pt",100,0,100);
	histgenFeta=fs->make<TH1F>("Gen_Muon_least_eta","Gen_Muon_least_eta",100,-10,10);
	histgenSeta=fs->make<TH1F>("Gen_Muon_mid_eta","Gen_Muon_mid_eta",100,-10,10);
	histgenTeta=fs->make<TH1F>("Gen_Muon_leading_eta","Gen_Muon_leading_eta",100,-10,10);
	histgenmaxdeltaR=fs->make<TH1F>("Gen_Muon_max_deltaR","Gen_Muon_max_deltaR",100,0,0.5);
	histgenmindeltaR=fs->make<TH1F>("Gen_Muon_min_deltaR","Gen_Muon_min_deltaR",100,0,0.5);
	histrecomaxISO=fs->make<TH1F>("reco_Muon_max_ISO","reco_Muon_max_ISO",100,0,1);
	//histmuongenallpt=fs->make<TH1F>("hist_muon_gen_all_pt","hist_muon_gen_all_pt",0,100);
	histgenmaxrecodeltaR=fs->make<TH1F>("reco_Muon_max_deltaR","reco_Muon_max_deltaR",100,0,0.5);
	histgenminrecodeltaR=fs->make<TH1F>("reco_Muon_min_deltaR","reco_Muon_min_deltaR",100,0,0.5);
	histrecoFpt=fs->make<TH1F>("reco_Muon_least_pt","reco_muon_least_Pt",100,0,100);
	histrecoSpt=fs->make<TH1F>("reco_Muon_mid_pt","reco_muon_mid_Pt",100,0,100);
	histrecoTpt=fs->make<TH1F>("reco_Muon_leading_pt","reco_Muon_leading_Pt",100,0,100);
	histrecoFeta=fs->make<TH1F>("reco_Muon_least_eta","reco_Muon_least_eta",100,-10,10);
	histrecoTeta=fs->make<TH1F>("reco_Muon_mid_eta","reco_Muon_mid_eta",100,-10,10);
	histMETpt=fs->make<TH1F>("MET_pt","MET_pt",100,0,100);
	histMETphi=fs->make<TH1F>("MET_phi_4th_mu","MET_phi_4th_mu",100,-3.2,3.2);
	histMETinvariantmass=fs->make<TH1F>("MET_included_invariant_mass","MET_included_invariant_mass",100,0,200);
	histrecoTeta=fs->make<TH1F>("reco_Muon_mid_eta","reco_Muon_mid_eta",100,-10,10);
	histgenrecodeltaRmax=fs->make<TH1F>("gen_reco_deltaR_max","gen_reco_deltaR_max",1000,0,0.1);
	histgenrecodeltaRmin=fs->make<TH1F>("gen_reco_deltaR_min","gen_reco_deltaR_min",1000,0,0.1);
	histrecodeltaRaccept=fs->make<TH1F>("reco_deltaR_accept","reco_deltaR_accept",6,0,6);
	histrecomatchingaccept=fs->make<TH1F>("reco_matching_accept","reco_matching_accept",6,0,6);
	histrecopTaccept=fs->make<TH1F>("reco_pT_accept","reco_pT_accept",6,0,6);
	histrecoetaaccept=fs->make<TH1F>("reco_eta_accept","reco_eta_accept",6,0,6);
	histrecopfaccept=fs->make<TH1F>("reco_pf_accept","reco_pf_accept",6,0,6);
	histgenothermupttau=fs->make<TH1F>("gen_other_mu_pt_from_tau","gen_other_mu_pt_from_tau",100,0,100);
	histgenothermupthar=fs->make<TH1F>("gen_other_mu_pt_from_other","gen_other_mu_pt_from_tau",100,0,100);
	histgenothermudeltaRmax=fs->make<TH1F>("gen_other_mu_deltaR_TAU","gen_other_mu_deltaR_TAU",100,0,5);
	histgenothermudeltaRmin=fs->make<TH1F>("gen_other_mu_deltaR_HARD","gen_other_mu_deltaR_HAR",100,0,5);
	histgenothermumother=fs->make<TH1F>("gen_mu_mother","gen_mu_mother",40,-20,20);
	histgenothermunumber=fs->make<TH1F>("gen_other_mu_number","gen_other_mu_number",10,0,10);
}


TestAnalyzer::~TestAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)


}



//
// member functions
//

// ------------ method called for each event  ------------
void
TestAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;
   // using namespace pat;
   int count=0; //this is the count for 3 muons from tau
   bool left[50];
   Handle<GenParticleCollection> genParticles;
 //  Handle<MuonCollection> muons;
 //  Handle<reco::PFMETCollection> pfMET;
	 Handle<pat::METCollection> pfMET;
	//iEvent.getByLabel("caloMet", met);
	 // Handle<View<CaloMET>> caloMEThandle;
   
	 iEvent.getByLabel("slimmedMETs",pfMET);
 	 Handle<pat::MuonCollection> muons;
   iEvent.getByLabel("slimmedMuons", muons); 
	 iEvent.getByLabel("prunedGenParticles", genParticles);
 //  iEvent.getByLabel("genParticles", genParticles);
//   iEvent.getByLabel("muons", muons);
	 float met_px=(*(pfMET->begin())).px();
	 float met_py=(*(pfMET->begin())).py();
	 //float met_mass=(*(pfMET->begin())).mass();
	 float met_energy=(*(pfMET->begin())).energy();
 	 float met_pt=(*(pfMET->begin())).pt();
 	 float met_phi=(*(pfMET->begin())).phi();
   int noOfmuons=muons->size();
   int NoGenMuon=0;
   float  deltaR=0;
   double max_pt_gen=0;
   double max_pt_reco=0;
   double min_pt_gen=2000;
   double min_pt_reco=2000;
   double px_g=0;
   double py_g=0;
   double pz_g=0;
   double energy_g=0;
   double mass_g=0;
   double pt_gen[50];
   double eta_gen[50];
   double phi_gen[50];
   double charge_gen[50];
   int accept_reco_1=noOfmuons;
   int accept_reco_2=noOfmuons;
   double mass=0;
   int positive=1;
   double reco_px[50];
   double reco_py[50];
   double reco_pz[50];
   double reco_energy[50];
   double gen_other_pt[50];
   double gen_other_px[50];
   double gen_other_py[50];
   double gen_other_pz[50];
   double gen_other_energy[50];
   double gen_other_eta[50];
   double gen_other_phi[50];
   long gen_other_mother[50];
   int count_other_muon=0;
	for (int i=0;i<30;++i)  //initial the selection arrary;
	{
	left[i]=true;
	reco_pt[i]=0.0;
	reco_eta[i]=0.0;
	reco_phi[i]=0.0;
	}	


        for (int i=0;i<10;i++)  //initiallize
		{
		gen_muon_eta[i]=0.0;
		gen_muon_phi[i]=0.0;
		}
		e_m=0.0;
		px_m=0.0;
		py_m=0.0;
		pz_m=0.0;
		
	for(size_t i = 0; i < genParticles->size(); ++ i) //starting from Z analysis_gen 
	{
		const GenParticle & p = (*genParticles)[i];  //geting the partilce as p
	   			if(p.pdgId()==23) //Z is 23
					{
							int final=1;
							int k = p.numberOfDaughters();
							for(int j = 0; j < k; ++ j) 
								{
       					const Candidate * d = p.daughter(j);
       					int dauId = d->pdgId();
								if (dauId == 23)
									{
										final=0; //not final
										break;
									}		
								}
							if (final==1)
							{
								histZmass->Fill(p.mass());
								histZpt->Fill(p.pt());
								histZeta->Fill(p.eta());
								//int k = p.numberOfDaughters();
               	//cout<<"number of Child of Z"<<k<<endl;
							}		
					}	//end of if if pdg.ID=z
	}//end of Gennerator Loop 1.
//cout<<"fuck"<<endl;
for(size_t i = 0; i < genParticles->size(); ++ i) //starting from tau analysis_gen 
	{
		  int final=1;
			int count_muon=0;
			const GenParticle & p = (*genParticles)[i];  //geting the partilce as p
	   			if(p.pdgId()==15 || p.pdgId()== -15 ) //tau is 15
					{
							int k = p.numberOfDaughters();
	   					int   n=p.numberOfMothers();
							for (int j=0; j<n; ++j)
							{
								const Candidate * mon=p.mother(j);
								int monId=mon->pdgId();
								if (monId != 23 && monId!=p.pdgId())
								{
									final=0;
									break;
								}
							}
							for (int j=0;j<k;++ j)
							{
								const Candidate * d = p.daughter(j);
       					int dauId = d->pdgId();
								if (dauId==p.pdgId())
								{
									final=0;
								}
							}
							if (final == 1)
							{
								histZchild->Fill(p.pdgId());
							}
							for(int j = 0; j < k; ++ j) 
								{
       					const Candidate * d = p.daughter(j);
       					int dauId = d->pdgId();
								if (dauId != 13 && dauId != -13)
								{
									final=0;
									break;
								}
								if (dauId == 13 || dauId == -13)
								{
									count_muon++;
								}
								}
							
							if (p.pdgId()==-15 && count_muon==3)
							{
								positive=0;
							}	
						histWpt->Fill(k);	
            if (final==1 && count_muon==3)
							{
						//		int k = p.numberOfDaughters();
								histtaupt->Fill(p.pt());
								histtaueta->Fill(p.eta());			
								histtaumass->Fill(p.mass());
								} //end of if (final==1 && count_muon==3)
				}	//end of if if pdg.ID
	}//end of Gennerator Loop 2.
//cout<<"fuck"<<endl;
positive=3;
tau_px=0;
tau_py=0;
tau_pz=0;
tau_energy=0;
if (positive == 1) 
{
NoGenMuon=0;
for(size_t i = 0; i < genParticles->size(); ++ i) //starting from muon analysis_gen for muon
	{
		  int final=1;
			const GenParticle & p = (*genParticles)[i];  //geting the partilce as p
	//   	const Candidate * mon=p.mother();
	   	int   n=p.numberOfMothers();
			int   k=p.numberOfDaughters();
	
				if(p.pdgId()==13 || p.pdgId()==-13 ) //muon is 13
					{
				
						histmuongenallpt->Fill(p.pt());
						histmuongenalleta->Fill(p.eta());
						for (int j=0; j<n; ++j)
							{
								const Candidate * mon=p.mother(j);
								int monId=mon->pdgId();
                const Candidate * grandmon=mon->mother(0);
								int grandmonId=grandmon->pdgId();
								if (monId != 15 || k>0)
								{
									final=0;
								}
								if (grandmonId!=23 && grandmonId!=15)
								{
									final=0;
								}
							}
						int child=1;
		
						for (int j=0;j<k;++j)
						{
							const Candidate * child=p.daughter(j);
						  if (child->pdgId() == p.pdgId())	
							{
								child=0;
								break;
							}
						}
				
						if (final ==0 && child ==1)
						{
						 	gen_other_pt[count_other_muon]=p.pt();
						 	gen_other_eta[count_other_muon]=p.eta();
						 	gen_other_phi[count_other_muon]=p.phi();
						 	gen_other_px[count_other_muon]=p.px();
						 	gen_other_py[count_other_muon]=p.py();
						 	gen_other_pz[count_other_muon]=p.pz();
						 	gen_other_energy[count_other_muon]=p.energy();
							gen_other_mother[count_other_muon]=p.mother(0)->pdgId();
							count_other_muon++;		
							histgenothermumother->Fill(p.mother(0)->pdgId());
							if (p.mother(0)->pdgId()==-15)
							{
								tau_px=p.mother(0)->px();		
								tau_py=p.mother(0)->py();		
								tau_pz=p.mother(0)->pz();		
								tau_energy=p.mother(0)->energy();		
								tau_pt=p.mother(0)->pt();
								tau_eta=p.mother(0)->eta();
								tau_phi=p.mother(0)->phi();
							}
						}
						if (final ==1 && n==1)
						{
								histmugenpt->Fill(p.pt());
								histmugeneta->Fill(p.eta());
								histmugenphi->Fill(p.phi());
								pt_gen[count]=p.pt();
								eta_gen[count]=p.eta();
								phi_gen[count]=p.phi();
								charge_gen[count]=p.pdgId();
								count++;
								if (p.pt()>3.0 && p.eta()>-2.4 && p.eta()<2.4)
								{	
								histmuongencutpt->Fill(p.pt());
								histmuongencuteta->Fill(p.eta());
								gen_muon_eta[NoGenMuon]=p.eta();
                gen_muon_pt[NoGenMuon]=p.pt();
								gen_muon_phi[NoGenMuon]=p.phi();
		//						gen_muon_charge[NoGenMuon]=p.pdgId();
								px_g=px_g+p.px();
								py_g=py_g+p.py();
								pz_g=pz_g+p.pz();
								energy_g=energy_g+p.energy();
								NoGenMuon++;
					  		}
						}
					}
    
	} //end of gen_loop_3
}
//cout<<"fuck 8"<<endl;
if (positive == 0)  // analysis for the tau+ -15
{
		NoGenMuon=0;
    for(size_t i = 0; i < genParticles->size(); ++ i) //starting from muon analysis_gen for tau+
    {
        int final=1;
        const GenParticle & p = (*genParticles)[i];  //geting the partilce as p
        //   	const Candidate * mon=p.mother();
        int   n=p.numberOfMothers();
				int   k=p.numberOfDaughters();
        if(p.pdgId()==13 || p.pdgId()==-13 ) //muon is 13
        {
            histmuongenallpt->Fill(p.pt());
            histmuongenalleta->Fill(p.eta());
            for (int j=0; j<n; ++j)
            {
                const Candidate * mon=p.mother(j);
                int monId=mon->pdgId();
                const Candidate * grandmon=mon->mother(0);
                if (grandmon->pdgId() != 23 && grandmon->pdgId() != -15)
								{
								 final=0;
								}
								if (monId != -15 || k>0)
                {
                    final=0;
                }
            
						}
						int child=1;
						for (int j=0;j<k;++j)
						{
							const Candidate * child=p.daughter(j);
						  if (child->pdgId() == p.pdgId())	
							{
								child=0;
								break;
							}
						}
						if (final == 0 && child ==1)
						{
						 	gen_other_pt[count_other_muon]=p.pt();
						 	gen_other_eta[count_other_muon]=p.eta();
						 	gen_other_phi[count_other_muon]=p.phi();
							gen_other_mother[count_other_muon]=p.mother(0)->pdgId();
							gen_other_px[count_other_muon]=p.px();
						 	gen_other_py[count_other_muon]=p.py();
						 	gen_other_pz[count_other_muon]=p.pz();
						 	gen_other_energy[count_other_muon]=p.energy();
							count_other_muon++;		
							histgenothermumother->Fill(p.mother(0)->pdgId());
							if (p.mother(0)->pdgId()==15)
							{
								tau_px=p.mother(0)->px();		
								tau_py=p.mother(0)->py();		
								tau_pz=p.mother(0)->pz();		
								tau_energy=p.mother(0)->energy();		
								tau_pt=p.mother(0)->pt();
								tau_eta=p.mother(0)->eta();
								tau_phi=p.mother(0)->phi();

							}
						}
            if (final ==1 && n==1)
            {
                pt_gen[count]=p.pt();
                eta_gen[count]=p.eta();
                phi_gen[count]=p.phi();
                charge_gen[count]=p.pdgId();
                count++;
                histmugenpt->Fill(p.pt());
                histmugeneta->Fill(p.eta());
                histmugenphi->Fill(p.phi());
								
                if (p.pt()>3.0 && p.eta()>-2.4 && p.eta()<2.4)
                {
                    histmuongencutpt->Fill(p.pt());
                    histmuongencuteta->Fill(p.eta());
                    gen_muon_pt[NoGenMuon]=p.pt();
                    gen_muon_eta[NoGenMuon]=p.eta();
                    gen_muon_phi[NoGenMuon]=p.phi();
      //              gen_muon_charge[NoGenMuon]=p.pdgId();
                    px_g=px_g+p.px();
                    py_g=py_g+p.py();
                    pz_g=pz_g+p.pz();
                    energy_g=energy_g+p.energy();
                    NoGenMuon++;
                }
            }
        }
    } //end of gen_loop_4
	//cout<<"fuck99"<<endl;
}
//cout<<"fuck 1999"<<endl;
int reco_muon_count=0;
//cout<<"fuck 1999"<<endl;
histgenaccept->Fill(NoGenMuon);
//cout<<"fuck 1999"<<endl;
histgenothermunumber->Fill(count_other_muon);
//cout<<"fuck 1999"<<endl;
//histcount->Fill(count);
//cout<<count<<endl;
//other muon_gen_analysis with gen from 3 tau
if (NoGenMuon==3)
{
in_z=sqrt(((energy_g+tau_energy)*(energy_g+tau_energy))-(((px_g+tau_px)*(px_g+tau_px))+((py_g+tau_py)*(py_g+tau_py))+((pz_g+tau_pz)*(pz_g+tau_pz))));
histinZmass->Fill(in_z);
}

for (int i=0;i<count_other_muon;++i)
{
//	cout<<"fuck 1999"<<endl;
	
	//double max_deltaR_other=0;
	double min_deltaR_other=1000;
	for (int j=0;j<3;j++)
	{
			double delta_eta=gen_other_eta[i]-gen_muon_eta[j];
			double delta_phi=gen_other_phi[i]-gen_muon_phi[j];
			if (delta_phi < -M_PI) { delta_phi=delta_phi+2*M_PI;}
			if (delta_phi >  M_PI) { delta_phi=delta_phi-2*M_PI;}
			double deltaR=sqrt((delta_eta*delta_eta)+(delta_phi*delta_phi));
		//	if (max_deltaR_other <deltaR)
		//	{
		//		max_deltaR_other=deltaR;
		//	}
			if (min_deltaR_other >deltaR)
			{
				min_deltaR_other=deltaR;
			}
	}
	if (gen_other_mother[i]==15 || gen_other_mother[i]== -15)
	{
	histgenothermudeltaRmax->Fill(min_deltaR_other);
	histgenothermupttau->Fill(gen_other_pt[i]);
	//histmutaupt->Fill(tau_pt);
	histmutauptratio->Fill(gen_other_pt[i]/tau_pt);
	double delta_eta=gen_other_eta[i]-tau_eta;
	double delta_phi=gen_other_phi[i]-tau_phi;
	if (delta_phi < -M_PI) { delta_phi=delta_phi+2*M_PI;}
	if (delta_phi >  M_PI) { delta_phi=delta_phi-2*M_PI;}
	double deltaR=sqrt((delta_eta*delta_eta)+(delta_phi*delta_phi));
	histmutaudeltaR->Fill(deltaR);
	double combinedpt=(gen_other_px[i]+tau_px)*(gen_other_px[i]+tau_px)+(gen_other_py[i]+tau_py)*(gen_other_py[i]+tau_py);
	combinedpt=sqrt(combinedpt);
	histmutaupt->Fill(combinedpt);
	if (gen_other_pt[i]>3.5 && abs(gen_other_eta[i])<2.4 && NoGenMuon==3)
		{
	 	double px_4=px_g+gen_other_px[i];
	 	double py_4=py_g+gen_other_py[i];
	 	double pz_4=pz_g+gen_other_pz[i];
		double p_4=px_4*px_4+py_4*py_4+pz_4*pz_4;
		double en_4=energy_g+gen_other_energy[i];
		hist4mmass->Fill(sqrt((en_4*en_4)-p_4));		
		}

	
	}
	else{
	histgenothermudeltaRmin->Fill(min_deltaR_other);
	histgenothermupthar->Fill(gen_other_pt[i]);
	}
}
//cout<<count<<"  "<<positive<<"  "<<NoGenMuon<<endl;
if (count>3)
{
//cout<<"fuck"<<count<<endl;
}
if (NoGenMuon==3) // nothing happen if the gen_didn't pass//
{
float max_deltaR=0;
float min_deltaR=10;
for (int i=0;i<3;i++)
	{
		for (int j=i+1;j<3;j++)
		{
			double delta_eta=gen_muon_eta[i]-gen_muon_eta[j];
			double delta_phi=gen_muon_phi[i]-gen_muon_phi[j];
			if (delta_phi < -M_PI) { delta_phi=delta_phi+2*M_PI;}
			if (delta_phi >  M_PI) { delta_phi=delta_phi-2*M_PI;}
			deltaR=sqrt((delta_eta*delta_eta)+(delta_phi*delta_phi));
			if (max_deltaR <deltaR)
			{
				max_deltaR=deltaR;
			}
		
			if (min_deltaR >deltaR)
			{
				min_deltaR=deltaR;
			}
		}
	}
histgenmaxdeltaR->Fill(max_deltaR);
histgenmindeltaR->Fill(min_deltaR);
histminptgen->Fill(min_pt_gen);
histmaxptgen->Fill(max_pt_gen);
double p_g=(px_g*px_g)+(py_g*py_g)+(pz_g*pz_g);
mass_g=sqrt((energy_g*energy_g)-p_g);
histgmass->Fill(mass_g);
for (int i=0;i<NoGenMuon;i++)
{
	histgenmuongrouppt->Fill(gen_muon_pt[i]);
	histgenmuongroupeta->Fill(gen_muon_eta[i]);
	histgenmuongroupphi->Fill(gen_muon_phi[i]);
}
}  //end of if NoofGenmuon==3

double maxISO=0;
int count_reco=0;
double deltaRmax=0;
double deltaRmin=0.2;
//int recodeltaRcut=0;
int recopTcut=0;
int recoetacut=0;
int recopfcut=0;
int recomatchingcut=0;
for(size_t i = 0; i < muons->size(); ++ i) //analysis on Muon_Reco
	{
		const Muon & p = (*muons)[i];  //geting the partilce as p
		//float pt=p.pt();
		//float eta=p.eta();
		bool iden=false;
		
//		if (pt >10000.0)

	//iden=true;	

	if (left[i])
		{
			histmuallpt->Fill(p.pt());	
			histmualleta->Fill(p.eta());
//		  pt_reco[count_reco]=p.pt();
	//	  eta_reco[count_reco]=p.eta();
			count_reco++;
		}
	if (left[i]==true && p.pt()<min_pt_reco)
		{
			min_pt_reco=p.pt();
		}
	if (left[i]==true && p.pt()>max_pt_reco)
		{
			  max_pt_reco=p.pt();
		}


		if (p.pt()<3 && left[i]==true)
		{
			left[i]=false;
			noOfmuons--;
			recopTcut++;
		}
		if (abs(p.eta())>2.4 && left[i]==true)
		{
			left[i]=false;
			noOfmuons--;
			recoetacut++;
		}
		if (left[i] && accept_reco_1>2 )
		{
			histmumatchpt->Fill(p.pt());	
			histmumatcheta->Fill(p.eta());
		}	

		if (p.isPFMuon()==false && left[i]==true)
		{
			left[i]=false;
			recopfcut++;
			noOfmuons--;
		}
		

		for (int j=0;j<-1;++j)   //partical identification
		{
			double delta_eta=p.eta()-eta_gen[j];
			double delta_phi=p.phi()-phi_gen[j];
			if (delta_phi < -M_PI) { delta_phi=delta_phi+2*M_PI;}
			if (delta_phi >  M_PI) { delta_phi=delta_phi-2*M_PI;}
			double pt_match=abs((p.pt()-pt_gen[j]))/p.pt();
			deltaR=sqrt((delta_eta*delta_eta)+(delta_phi*delta_phi));
			if (deltaR<0.05 && p.pdgId()==charge_gen[j] && pt_match<0.1)
			{
				if (deltaR>deltaRmax)
				{
					deltaRmax=deltaR;
				}
			
				if (deltaR<deltaRmin)
				{
					deltaRmin=deltaR;
				}
			}	
			if (deltaR<0.05 && p.pdgId()==charge_gen[j] && pt_match<0.1)
			//if (delta_eta <0.01 && delta_phi<0.01)
			{
				iden=true;				
			}
		}
	iden=true;
	if (iden == false && left[i]==true)
		{
			left[i]=false;
			recomatchingcut++;
			noOfmuons--;
		}

		

		if (left[i] == false)
		{
			accept_reco_1--;
		}
		double iso1 = (p.pfIsolationR04().sumChargedHadronPt+ TMath::Max(p.pfIsolationR04().sumNeutralHadronEt+p.pfIsolationR04().sumPhotonEt-0.5*p.pfIsolationR04().sumPUPt, 0.))/p.pt();//isolation parmeter
		if (iso1>0.4 && left[i])
		{
			left[i]=false;
			recopfcut++;
			noOfmuons--;
		}
	
	if (iso1>maxISO && left[i]) 
		{
			maxISO=iso1;
		}	
	
		if (left[i] == false)
		{
			accept_reco_2--;
		}
	 if (left[i] && accept_reco_2>2)
	{
			histmulastpt->Fill(p.pt());	
			histmulasteta->Fill(p.eta());
	}
		
		if (left[i]==true)
		{

					reco_px[reco_muon_count]=p.px();
					reco_py[reco_muon_count]=p.py();
					reco_pz[reco_muon_count]=p.pz();
					reco_energy[reco_muon_count]=p.energy();
					reco_pt[reco_muon_count]=p.pt();
					reco_eta[reco_muon_count]=p.eta();
					reco_phi[reco_muon_count]=p.phi();
					reco_muon_count++;
					histmupt->Fill(p.pt());					
					histmueta->Fill(p.eta());					
					histmuphi->Fill(p.phi());					
					histgenrecodeltaRmax->Fill(deltaRmax);
					histgenrecodeltaRmin->Fill(deltaRmin);

			}
	} //end of muon loop

//plot the acceptance
histrecopTaccept->Fill(muons->size()-recopTcut);
histrecoetaaccept->Fill(muons->size()-recopTcut-recoetacut);
histrecopfaccept->Fill(muons->size()-recopTcut-recoetacut-recopfcut);
histrecomatchingaccept->Fill(muons->size()-recomatchingcut-recopTcut-recoetacut-recopfcut);
histrecoaccept->Fill(noOfmuons);
int reco_deltaR_cut_after=0;
double f_mu_pt[10];
double f_mu_px[10];
double f_mu_py[10];
double f_mu_pz[10];
double f_mu_energy[10];
double f_mu_phi[10];
int f_mu_number=0;
for (int i=0;i<noOfmuons;i++)
		{
			//double sum_deltaR=0;	
			int btb=0;
			int deltaR_count=0;
			for (int j=0;j<noOfmuons;j++)
			{
				double delta_eta=reco_eta[i]-reco_eta[j];
				double delta_phi=reco_phi[i]-reco_phi[j];
				if (delta_phi < -M_PI) { delta_phi=delta_phi+2*M_PI;}
				if (delta_phi >  M_PI) { delta_phi=delta_phi-2*M_PI;}
			 	double	deltaR=sqrt((delta_eta*delta_eta)+(delta_phi*delta_phi));
				//sum_deltaR=sum_deltaR+deltaR;
				if (deltaR>0.25)
				{
					deltaR_count++;
				//	sum_deltaR=sum_deltaR+deltaR;
				}
				if (deltaR>1.57 && deltaR<4.71)
				{
				  btb++;
				}
			}
			if (deltaR_count>1)
			{
				if (btb>=3)
				{
			 	f_mu_pt[f_mu_number]=reco_pt[i];
			 	f_mu_px[f_mu_number]=reco_px[i];
			 	f_mu_py[f_mu_number]=reco_py[i];
				f_mu_pz[f_mu_number]=reco_pz[i];
				f_mu_phi[f_mu_number]=reco_phi[i];
				f_mu_energy[f_mu_number]=reco_energy[i];
				f_mu_number++;
				}
				reco_deltaR_cut_after++;
				reco_pt[i]=0.0;
			}
		}
double f_max_pt=0;
for (int i=0;i<f_mu_number;i++)
{
	if (f_mu_pt[i]>f_max_pt)
	{
		f_mu_pt[0]=f_mu_pt[i];
		f_mu_px[0]=f_mu_px[i];
		f_mu_py[0]=f_mu_py[i];
		f_mu_pz[0]=f_mu_pz[i];
		f_mu_phi[0]=f_mu_phi[i];
		f_mu_energy[0]=f_mu_energy[i];
		f_max_pt=f_mu_pt[i];
	}
}

//noOfmuons=noOfmuons+reco_deltaR_cut_after;
if (noOfmuons >= 3) //selection of the largest 3 pT muon in Reco
 {	
		for (int i=0;i<noOfmuons;++i)
		{
			for (int j=i+1;j<noOfmuons;++j)
			{
				if (reco_pt[j]>reco_pt[i])
				{
					double temp_px=reco_px[i];
					double temp_py=reco_py[i];
					double temp_pz=reco_pz[i];
					double temp_e=reco_energy[i];
					double temp_pt=reco_pt[i];
					double temp_eta=reco_eta[i];
					double temp_phi=reco_phi[i];
					reco_px[i]=reco_px[j];
					reco_py[i]=reco_py[j];
					reco_pz[i]=reco_pz[j];
					reco_energy[i]=reco_energy[j];
					reco_pt[i]=reco_pt[j];
					reco_eta[i]=reco_eta[j];
					reco_phi[i]=reco_phi[j];
					reco_px[j]=temp_px;
					reco_py[j]=temp_py;
					reco_pz[j]=temp_pz;
					reco_energy[j]=temp_e;
					reco_pt[j]=temp_pt;
					reco_eta[j]=temp_eta;
					reco_phi[j]=temp_phi;
				}
			}	
		}
} 


		
noOfmuons=noOfmuons-reco_deltaR_cut_after;
histrecodeltaRaccept->Fill(muons->size()-recomatchingcut-reco_deltaR_cut_after-recopTcut-recoetacut-recopfcut);
//histrecodeltaRaccept->Fill(noOfmuons);
histrecoaccept->Fill(noOfmuons);

	//plot the kinematic of reco muons
	for (int i=0;i<3;++i)
	{
		px_m=px_m+reco_px[i];
		py_m=py_m+reco_py[i];
		pz_m=pz_m+reco_pz[i];
		e_m=e_m+reco_energy[i];
		double p=(px_m*px_m)+(py_m*py_m)+(pz_m*pz_m);
		mass=sqrt((e_m*e_m)-(p));
	}
if (noOfmuons>3) { noOfmuons=3;	}
if (noOfmuons==3)
	{
		for (int i=0; i<reco_muon_count; i++)
		{
		histrecomaxISO->Fill(maxISO);
		histrecomuongrouppt->Fill(reco_pt[i]);					
		histrecomuongroupeta->Fill(reco_eta[i]);					
		histrecomuongroupphi->Fill(reco_phi[i]);					
		}
		
		for (int i=0; i<3; i++)
		{
				for (int j=i+1;j<3;j++)
				{
					if (reco_pt[i]>reco_pt[j])
					{
						double temp=reco_pt[j];
						reco_pt[j]=reco_pt[i];
						reco_pt[i]=temp;
					}
				}
		}
		histrecoFpt->Fill(reco_pt[0]);
		histrecoSpt->Fill(reco_pt[1]);
		histrecoTpt->Fill(reco_pt[2]);
			
		float max_deltaR=0;
		float min_deltaR=10;
		for (int i=0;i<3;i++)
		{
			for (int j=i+1;j<3;j++)
			{
				double delta_eta=reco_eta[i]-reco_eta[j];
				double delta_phi=reco_phi[i]-reco_phi[j];
				if (delta_phi < -M_PI) { delta_phi=delta_phi+2*M_PI;}
				if (delta_phi >  M_PI) { delta_phi=delta_phi-2*M_PI;}
			 	double	deltaR=sqrt((delta_eta*delta_eta)+(delta_phi*delta_phi));
				if (max_deltaR <deltaR)
				{
					max_deltaR=deltaR;
				}
				if (min_deltaR >deltaR)
				{
				min_deltaR=deltaR;
				}
			}
		}
		histgenmaxrecodeltaR->Fill(max_deltaR);
		histgenminrecodeltaR->Fill(min_deltaR);
		histDmass->Fill(mass);
		if (f_mu_number >0)
		{
		hist4mpt->Fill(f_mu_pt[0]);
		}
		if (f_mu_number>0 && f_mu_pt[0]>3.5)
		{
		px_m=px_m+f_mu_px[0];
		py_m=py_m+f_mu_py[0];
		pz_m=pz_m+f_mu_pz[0];
		e_m=e_m+f_mu_energy[0];
		double p=(px_m*px_m)+(py_m*py_m)+(pz_m*pz_m);
		mass=sqrt((e_m*e_m)-(p));
		histreco4mmass->Fill(mass);
		histMETpt->Fill(met_pt);
		met_px=met_px+px_m;
		met_py=met_py+py_m;
		met_energy=met_energy+e_m;
	  p=(met_px*met_px)+(met_py*met_py)+(pz_m*pz_m);
		mass=sqrt((met_energy*met_energy)-(p));
		histMETinvariantmass->Fill(mass);
		double delta_phi=f_mu_phi[0]-met_phi;
		if (delta_phi < -M_PI) { delta_phi=delta_phi+2*M_PI;}
		if (delta_phi >  M_PI) { delta_phi=delta_phi-2*M_PI;}
		histMETphi->Fill(delta_phi);
		}
} // if no of muon=3


histmaxptreco->Fill(max_pt_reco);
histminptreco->Fill(min_pt_reco);

//kinematic analysis on Gen level
          
        
       
       
}


// ------------ method called once each job just before starting event loop  ------------
void 
TestAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TestAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
TestAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
TestAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
TestAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
TestAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TestAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TestAnalyzer);
