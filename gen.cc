// -*- C++ -*-
//
// Package:    gen/gen
// Class:      gen
// 
/**\class gen gen.cc gen/gen/plugins/gen.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yuxiang Guo
//         Created:  Wed, 30 Sep 2015 23:43:08 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
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
#include <TCanvas.h>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
using namespace edm; 
using namespace reco;
using namespace std;
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class declaration
//

class gen : public edm::EDAnalyzer {
   public:
      explicit gen(const edm::ParameterSet&);
      ~gen();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      TH1F *histWpt;
      TH1F *histZchild;
      TH1F *histtaudecaymode;
			
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
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
gen::gen(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   Service<TFileService> fs;
   histWpt=fs->make<TH1F>("Wpt_gen","W+ pt",10,0,10);
   histZchild=fs->make<TH1F>("Z_tau_child","Z_tau_gen",40,-20,20);
   histtaudecaymode=fs->make<TH1F>("tau_decay_mode","tau_decay_mode",5,1,5);
}


gen::~gen()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
gen::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;
   Handle<GenParticleCollection> genParticles;
   iEvent.getByLabel("genParticles", genParticles);
   for(size_t i = 0; i < genParticles->size(); ++ i) //starting from tau analysis_gen 
	{
		  int final=1;
			int decay_mode_tau=0; //1 for 3mu,2 for 1mu 3 for 1 electron 4 for hardronic for others
			int count_muon=0;
			int count_electron=0;
			const GenParticle & p = (*genParticles)[i];  //geting the partilce as p
	   			if(p.pdgId()==15 || p.pdgId()== -15 ) //tau is 15
					{
							int k = p.numberOfDaughters();
	   					int   n=p.numberOfMothers();
							for (int j=0; j<n; ++j)
							{
								const Candidate * mon=p.mother(j);
								int monId=mon->pdgId();
								if (monId != 23 && monId!=p.pdgId()) //check if tau from Z
								{
									final=0;
									break;
								}
							}
							for (int j=0;j<k;++ j) //check whether they duplicate another tau
							{
								const Candidate * d = p.daughter(j);
       					int dauId = d->pdgId();
								if (dauId==p.pdgId())
								{
									final=0;
									break;
								}
								if (dauId==13 || dauId==-13)
								{
									count_muon++;
								}
								
								if (dauId==11 || dauId==-11)
								{
									count_electron++;
								}
							}
							if (final == 1)
							{
								histZchild->Fill(p.pdgId());
								if (count_muon == 3)
								{
									decay_mode_tau=1;
								} 		
								if (count_muon == 1)
								{
									decay_mode_tau=2;
								} 		
								if (count_electron == 1)
								{
									decay_mode_tau=3;
								} 		
								if (count_muon == 0 && count_electron==0)
								{
									decay_mode_tau=4;
								} 		
								histtaudecaymode->Fill(decay_mode_tau);
							}
							
							if (p.pdgId()==-15 && count_muon==3)
							{
								positive=0;
							}	
            
						if (final==1 && count_muon==3)
							{
						//		int k = p.numberOfDaughters();
						//		histtaupt->Fill(p.pt());
						//		histtaueta->Fill(p.eta());			
						//		histtaumass->Fill(p.mass());
							} //end of if (final==1 && count_muon==3)
				}	//end of if if pdg.ID
	}//end of Tau Analysis.
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
} //end if positive 
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
gen::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
gen::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
gen::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
gen::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
gen::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
gen::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
gen::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(gen);
