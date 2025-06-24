#include <iostream>
#include <vector>
#include <cmath>
#include "fastjet/ClusterSequence.hh"
#include "TH1.h"
#include "TMath.h"
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TFile.h"
#include "Pythia8/Pythia.h"

using namespace fastjet;
using namespace Pythia8;
using namespace std;

const Int_t nDeltaRBinsEEC = 13;
Float_t deltaRBinsEEC[nDeltaRBinsEEC + 1];

//!!!Don't forget to add to Makefile

//Function to compute Delta R
float deltaR(const PseudoJet& p1, const PseudoJet& p2) {
  float dphi = std::abs(p1.phi() - p2.phi());
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  float deta = p1.eta() - p2.eta();
  return std::sqrt(deta * deta + dphi * dphi);
}

 
int main(int argc, char* argv[]) {
    //Check if a .cmnd file is provided + allows for the .cmnd file to be accepted as an argument
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <config.cmnd>" << endl;
        return 1;
    } 
  
    TFile *out = new TFile("pythia_pp_eec_jun23.root", "RECREATE");
    out->cd();
    
    cout << "created pythia_pp_eec_jun23.root!" << endl;
    
    Pythia pythia;
    
    if (!pythia.readFile(argv[1])) {
        cerr << "Error: Could not read file.cmnd!" << endl;
        return 1;
    }
    
    pythia.init();
    
    float jet_radius = 0.4;
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, jet_radius);

        
    //declare power vectors
    vector<float> v1 = {0.25, 0.5, 0.75, 0.5, 1, 1, 0.5, 2, 1.25, 1.5, 1, 2, 2};
    vector<float> v2 = {0.25, 0.5, 0.75, 1, 0.5, 1, 2, 0.5, 1.25, 1.5, 2, 1, 2};
    
    //Declare-create histogram variables
    string histname;
    string histtitle;
    TH1F* hists[v1.size()];
    
    for(int i = 0; i < v1.size(); i++){
            histname = "E" + to_string(v1[i]) + "_E" + to_string(v2[i]) + ", " ;
            histtitle = "E" + to_string(v1[i]) + " E" + to_string(v2[i]) + " histogram" ;
            hists[i] = new TH1F(histname.c_str(), histtitle.c_str(), 20, -8, 1);
            cout << histname << "created" << endl;
        
            }

    //Event loop
    for (int iEvent = 0; iEvent < 1000; ++iEvent) { // 1000 events for now
        if (!pythia.next()) continue;
        
        vector<fastjet::PseudoJet> event;
        
        // pT cut > 120 & < 140. GeV
        for (int i = 0; i < pythia.event.size(); ++i) {
            Particle& p = pythia.event[i];
            if (p.isFinal() && p.pT() > 0.3) {
	            event.push_back( PseudoJet(p.px(), p.py(), p.pz(), p.e()));
	            }
            }
        
        fastjet::ClusterSequence cs(event, jet_def);
        vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());

        // Require one jet
        if (jets.size() < 1) 
            continue;
        
        //Jet pT cut
        if(jets.at(0).pt() < 120 || jets.at(0).pt() < 140)
            continue;
        
        vector<fastjet::PseudoJet> jet_constituents = jets.at(0).constituents();
            
        for(int k = 0; k < v1.size(); k++){
            for (size_t i = 0; i < jet_constituents.size(); ++i) {
                for (size_t j = i + 1; j < jet_constituents.size(); ++j) {
                    
                    if(jet_constituents.at(i).pt() > 1 && jet_constituents.at(j).pt() > 1){
                        float eec = pow(jet_constituents.at(i).pt(), v1[k]) * pow(jet_constituents.at(j).pt(), v2[k]);
                        //float eec = jet_constituents.at(i).pt() * jet_constituents.at(j).pt();
                    
                        float delr = jet_constituents.at(i).delta_R(jet_constituents.at(j));
	                    //float z = (1 - TMath::Cos(delr))/2; 
	                
                        // Epow Epow filling
                        hists[k]->Fill(TMath::Log(delr), eec);
	                    }
                    }
                }
            }
    } //Event loop close
    
    out->Write();
  
    //Pythia cleanup
    pythia.stat();
    out->Close();
    return 0;    
        
}//Main close