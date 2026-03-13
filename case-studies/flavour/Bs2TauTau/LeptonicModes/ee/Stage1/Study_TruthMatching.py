import ROOT as r
from termcolor import colored
import json
import numpy as np

r.gInterpreter.Declare('''
    string Translate_PDG(ROOT::VecOps::RVec<int> Daughters_ind, ROOT::VecOps::RVec<int> MC_PDG){
        
        string result ("");

        if (Daughters_ind.at(0) == -1){
            result = "No CA";
            return result;
        }
        
        int i (0);
        do{ 
            bool Mother (true);
            int j (0);
            string subdecay ("");
            do {
                int PDG (std::abs(MC_PDG.at(Daughters_ind.at(i+j))));
                if (PDG < 100) {
                    if (PDG == 11) subdecay += "E";
                    else if (PDG == 13) subdecay += "Mu";
                    else if (PDG == 15) subdecay += "Tau";
                    else if (PDG == 12) subdecay += "Nue";
                    else if (PDG == 14) subdecay += "Numu";
                    else if (PDG == 16) subdecay += "Nutau";
                    else if (PDG == 22) subdecay += "Gamma";
                    else subdecay += std::to_string(PDG);
                }
                else {
                    if (PDG == 130) subdecay += "K"; //KL
                    else if (PDG < 300) subdecay += "Pi";
                    else if (PDG == 333) subdecay += "Phi";
                    else if (PDG < 400) subdecay += "K";
                    else if (PDG == 443) subdecay += "J/Psi";
                    else if (PDG < 500) subdecay += "D";
                    else if (PDG < 600) subdecay += "B";
                    else if (PDG == 10433) subdecay += "D";
                    else if (PDG == 10411) subdecay += "D";
                    else if (PDG == 20433) subdecay += "D";
                    else subdecay += std::to_string(PDG);
                }
                if (Mother){
                    subdecay += "->";
                    Mother = false;
                }
                ++j;
            } while (Daughters_ind.at(i+j) != -2);
            result += subdecay + "|";
            i = i+j+1;
        } while (i < Daughters_ind.size());

        return result;
    }

    ROOT::VecOps::RVec<int> Decay_Chain(ROOT::VecOps::RVec<int> Daughters_ind){
        
        ROOT::VecOps::RVec<int> result;

        if (Daughters_ind.at(0) == -1){
            result.push_back(-1);
            return result;
        }

        int i (0);
        do{ 
            ROOT::VecOps::RVec<int> temp;
            int j (0);
            do {
                temp.push_back(Daughters_ind.at(i+j));
                ++j;
            } while (Daughters_ind.at(i+j) != -2);
            if (temp.size() > 2){
                for (size_t k=0; k<temp.size(); ++k){
                    result.push_back(temp.at(k));
                }
                result.push_back(-2);
            }
            i = i+j+1;
        } while (i < Daughters_ind.size());

        return result;

    }
''')


#Get the histograms in a numpy form -----------------------------------
def Load_RDF(Mode):
        
    Path = "/afs/cern.ch/work/t/tomonnar/public/Bs2TauTau/ee/eede/"

    Links = {"bb":"p8_ee_Zbb_ecm91",
             "cc":"p8_ee_Zcc_ecm91",
             "ss":"p8_ee_Zss_ecm91",
             "ud":"p8_ee_Zcc_ecm91",
             "sig":"p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau",
            }

    #Load the rdf properly
    rdf = r.RDataFrame("events",Path+Links[Mode]+".root")

    return rdf



#=========================================================================

rdf = Load_RDF("sig")
rdf.Filter("Selected>0").Display(["nevent","TM_RECO_electronwithDuplicates_ind","TM_RECO_positronwithDuplicates_ind"]).Print()


#print(f"Number of events : {rdf.Count().GetValue()}")
#print(f"Number of events with more than 2 TM electrons : {rdf.Filter('n_TM_electrons > 2').Count().GetValue()}")
#print(f"Number of events with more than 1 TM e+ : {rdf.Filter('n_TM_positron > 1').Count().GetValue()}")
#print(f"Number of events with more than 1 TM e- : {rdf.Filter('n_TM_electron > 1').Count().GetValue()}")
#print(f"Number of events with more than 3 TM electrons : {rdf.Filter('n_TM_electrons > 3').Count().GetValue()}")


#rdf = rdf.Define("MC_dielectron_CADaughters","Translate_PDG(MC_dielectron_CADaughters_ind,MC_PDG)")
#rdf = rdf.Define("MC_dielectron_CADaughters_ind_reduced","Decay_Chain(MC_dielectron_CADaughters_ind)")
#rdf = rdf.Define("MC_dielectron_CADaughters_reduced","Translate_PDG(MC_dielectron_CADaughters_ind_reduced,MC_PDG)")

#rdf.Filter("n_TM_electrons > 2").Display(["MCRecoAsso","MCRecoAssoo","TM_MC_positron_ind","TM_RECO_positron_ind","Electroo","Photoo"],4,100).Print()
#rdf.Display(["MCRecoAsso","MCRecoAssoo","TM_MC_positron_ind","TM_RECO_positron_ind","Electroo","Photoo"],4,100).Print()

#rdf.Display(["MC_dimuon_CADaughters_reduced"],30).Print()