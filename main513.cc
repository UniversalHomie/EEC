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

//!!!Don't forget to add to Makefile

//Function to compute Delta R
float deltaR(const PseudoJet& p1, const PseudoJet& p2) {
  float dphi = std::abs(p1.phi() - p2.phi());
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  float deta = p1.eta() - p2.eta();
  return std::sqrt(deta * deta + dphi * dphi);
}

 
int main(int argc, char* argv[]) {

    //Shower types loop open
    for (int showercount = 1; showercount <= 3; showercount++){

        //Check if a .cmnd file is provided + allows for the .cmnd file to be accepted as an argument
        if (argc < 2) {
            cerr << "Usage: " << argv[0] << " <config.cmnd>" << endl;
            return 1;
        } 
    
        //Shower type stuff
        string showty_st = "PartonShowers:model = ";
        string showertype = showty_st + to_string(showercount);
        vector<string> showerlist = {"standard","vincia","dire"};
        cout << showertype << endl;
        cout << "shower type is " << showerlist[showercount - 1] << endl;
    
        string date = "jul21_";
        string filename = "eec_" + date + showerlist[showercount - 1] + ".root";
        TString TFilename(filename);
      
        TFile *out = new TFile(TFilename, "RECREATE");
        out->cd();
        
        cout << "created " + filename + ".root!" << endl;
        
        Pythia pythia;
        
        
        pythia.readString(showertype);
        
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
        
        //declare pt cut vector
        vector<int> pcut = {1,2,3,4,5};
        
       
        //Declare-create histogram variables
        string histname;
        string histtitle;
        string pcutname;
        string pcuttitle;
        TH1F* hists[v1.size()];
        TH1F* charged_hist[0];
        TH1F* pcut_hists[pcut.size()];
        TH1F* pcutregent[v1.size()][pcut.size()];
    
        for(int i = 0; i < v1.size(); i++){
            histname = "Charged_E" + to_string(v1[i]) + "_E" + to_string(v2[i]) + ", " ;
            histtitle = "Charged E" + to_string(v1[i]) + " E" + to_string(v2[i]) + " histogram" ;
            hists[i] = new TH1F(histname.c_str(), histtitle.c_str(), 60, -8, 1);
            cout << histname << "created" << endl;
            
            for(int j = 0; j < pcut.size(); j++){
                pcutname = "Pt_" + to_string(pcut[j]) + "_Charged_E" + to_string(v1[i]) + "_E" + to_string(v2[i]) + ", " ;
                pcuttitle = "Pt_" + to_string(pcut[j]) + "Charged E" + to_string(v1[i]) + " E" + to_string(v2[i]) + " histogram" ;
                pcutregent[i][j] = new TH1F(pcutname.c_str(), pcuttitle.c_str(), 60, -8, 1);
                cout << pcutname << "created" << endl;
            
                }
            }
    
        charged_hist[0] = new TH1F("charged-delr", "charged delta r", 20, -2, 2);
    
        //Event loop
        for (int iEvent = 0; iEvent < 6000; ++iEvent) { // 6000 events 
            if (!pythia.next()) continue;
            
            vector<fastjet::PseudoJet> event;
            
            vector<fastjet::PseudoJet> neutral_particles;
            vector<fastjet::PseudoJet> charged_particles;
            
            // pT cut > 120 & < 140. GeV
            for (int i = 0; i < pythia.event.size(); ++i) {
                Particle& p = pythia.event[i];
                if (p.isFinal() && p.pT() > 0.3) {
    	            event.push_back( PseudoJet(p.px(), p.py(), p.pz(), p.e()));
    	            }
    	        if(p.isCharged()){
    	            charged_particles.push_back( PseudoJet(p.px(), p.py(), p.pz(), p.e()));
    	        }
    	        if(p.isNeutral()){
                    neutral_particles.push_back( PseudoJet(p.px(), p.py(), p.pz(), p.e()));
    	        }
                }
            
            fastjet::ClusterSequence cs(event, jet_def);
            vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
            fastjet::ClusterSequence cs1(neutral_particles, jet_def);
            vector<fastjet::PseudoJet> neutral_jets = fastjet::sorted_by_pt(cs1.inclusive_jets());
            fastjet::ClusterSequence cs2(charged_particles, jet_def);
            vector<fastjet::PseudoJet> charged_jets = fastjet::sorted_by_pt(cs2.inclusive_jets());
            
            // Require atleast one jet
            if (jets.size() < 1) 
                continue;
            
            //Jet pT cut
            if(jets.at(0).pt() < 120 || jets.at(0).pt() < 140)
                continue;
            
            vector<fastjet::PseudoJet> jet_constituents = jets.at(0).constituents();
            vector<fastjet::PseudoJet> neutral_constituents = neutral_jets.at(0).constituents();
            vector<fastjet::PseudoJet> charged_constituents = charged_jets.at(0).constituents();
            
            //Declaration for particles within a certain radius "_rad"
            
            vector<fastjet::PseudoJet> charged_rad;
            for (size_t c_i = 0; c_i < charged_particles.size(); c_i++){
                if (jets.at(0).delta_R(charged_particles.at(c_i)) < 0.8){
                    charged_rad.push_back(charged_particles.at(c_i));
                }
            }
            
            vector<fastjet::PseudoJet> neutral_rad;
            for (size_t c_i = 0; c_i < neutral_particles.size(); c_i++){
                if (jets.at(0).delta_R(neutral_particles.at(c_i)) < 0.8){
                    neutral_rad.push_back(neutral_particles.at(c_i));
                }
            }
    
                
            for(int k = 0; k < v1.size(); k++){
                for (size_t i = 0; i < charged_rad.size(); ++i) {
                    for (size_t j = i + 1; j < charged_rad.size(); ++j) {
                        
                        if (charged_rad.at(i).pt() > 1 && charged_rad.at(j).pt() > 1){
                            float eec = pow(charged_rad.at(i).pt(), v1[k]) * pow(charged_rad.at(j).pt(), v2[k]);
                            //float eec = charged_rad.at(i).pt() * charged_rad.at(j).pt();
                            
                        
                            float delr = charged_rad.at(i).delta_R(charged_rad.at(j));
    	                    //float z = (1 - TMath::Cos(delr))/2; 
    	                
                            // Epow Epow filling
                            hists[k]->Fill(TMath::Log(delr), eec);
                            
                            const double& crpti = charged_rad.at(i).pt();
                            const double& crptj = charged_rad.at(j).pt();
                            if (crpti < 2 && crptj < 2){
                                pcutregent[k][0]->Fill(TMath::Log(delr), eec);
                                }
                            if (crpti > 2 && crpti < 3 && crptj > 2 && crptj < 3){
                                pcutregent[k][1]->Fill(TMath::Log(delr), eec);
                                }
                            if (crpti > 3 && crpti < 4 && crptj > 3 && crptj < 4){
                                pcutregent[k][2]->Fill(TMath::Log(delr), eec);
                                }
                            if (crpti > 4 && crpti < 5 && crptj > 4 && crptj < 5){
                                pcutregent[k][3]->Fill(TMath::Log(delr), eec);
                                }
                            if (charged_rad.at(i).pt() > 5 && charged_rad.at(j).pt() > 5){
                                pcutregent[k][4]->Fill(TMath::Log(delr), eec);
                                }
                            
    	                    }
                        }
                    } 
                }
    
            for (int cj_i = 0; cj_i < charged_jets.size(); ++cj_i) {
                for (int nj_i = 0; nj_i < neutral_jets.size(); ++nj_i) {
                    float chargedelr = deltaR(charged_jets.at(cj_i), neutral_jets.at(nj_i));
                    charged_hist[0]->Fill(TMath::Log(chargedelr));
                    //charged_hist[0]->Fill(chargedelr);
                }
            }
        } //Event loop close
    
        
        out->Write();
      
        //Pythia cleanup
        pythia.stat();
        out->Close();
            
               } //Shower types loop close 
        return 0;
}//Main close