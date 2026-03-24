import ROOT as r
from termcolor import colored
import json
import numpy as np
import matplotlib.pyplot as plt


r.gInterpreter.Declare('''
    #include <string>

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
                    if (PDG == 11) subdecay += "e";
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

    string Plain_PDG(ROOT::VecOps::RVec<int> Daughters_ind, ROOT::VecOps::RVec<int> MC_PDG){
        
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
                subdecay += std::to_string(MC_PDG.at(Daughters_ind.at(i+j)));
                if (Mother){
                    subdecay += "->";
                    Mother = false;
                }
                else subdecay += ",";
                ++j;
            } while (Daughters_ind.at(i+j) != -2);
            result += subdecay + "|";
            i = i+j+1;
        } while (i < Daughters_ind.size());

        return result;
    }


    ROOT::VecOps::RVec<int> Find_MuSubdecays(ROOT::VecOps::RVec<int> Decays, ROOT::VecOps::RVec<int> MC_PDG, int pm){
        
        ROOT::VecOps::RVec<int> results;
        
        //No Common Ancestor case (to be refined with the daughters ind)
        if (Decays.at(0) == -1){
            results.push_back(-1);
        }
        
        else {
            
            int i (0);
            do{ 
                int j (0);
                do {
                    if (MC_PDG.at(Decays.at(i+j)) == pm*(-13) || MC_PDG.at(Decays.at(i+j)) == pm*(-11)){ //Start a new loop to save the PDG values of the muon/electron mother and sister
                        int k (0);
                        do {
                            results.push_back(MC_PDG.at(Decays.at(i+k)));
                            ++k;
                        } while (Decays.at(i+k) != -2);
                        j=k-1; //skip the decay since already scanned
                    }
                    ++j;
                } while (Decays.at(i+j) != -2);
                i = i+j+1;
            } while (i < Decays.size());
        }
        return results;
    }

    int Find_Categories(ROOT::VecOps::RVec<int> MuFamily){

        //No Common ancestor case
        if (MuFamily.size() == 1 && MuFamily.at(0) == -1){
            return 0;
        }

        else {
            if (500 < std::abs(MuFamily.at(0)) && std::abs(MuFamily.at(0)) < 600) return 1; //B mother
            else if ((400 < std::abs(MuFamily.at(0)) && std::abs(MuFamily.at(0)) < 500) || std::abs(MuFamily.at(0)) == 10433 || std::abs(MuFamily.at(0)) == 10411 || std::abs(MuFamily.at(0)) == 20433) return 2; //D or excited D state mother
            else if (300 < std::abs(MuFamily.at(0)) && std::abs(MuFamily.at(0)) < 400) return 3; //K mother
            else if (std::abs(MuFamily.at(0)) == 15) return 4; //Tau mother
            else return 5; //Other to be refined if most cases
        }

    }

''')

#Get the histograms in a numpy form -----------------------------------
def Load_RDF(Mode):
        
    Path = "./"

    Links = {"bb":"p8_ee_Zbb_ecm91",
             "cc":"p8_ee_Zcc_ecm91",
             "ss":"p8_ee_Zss_ecm91",
             "ud":"p8_ee_Zcc_ecm91",
             "sig":"p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau",
            }

    #Load the rdf properly
    rdf = r.RDataFrame("events",Path+Links[Mode]+".root")
    rdf = rdf.Define("MC_dilepton_CADaughters_PDG","Plain_PDG(MC_dilepton_CADaughters_ind,MC_PDG)")
    rdf = rdf.Define("MC_dilepton_CADaughters","Translate_PDG(MC_dilepton_CADaughters_ind,MC_PDG)")
    rdf = rdf.Define("MC_dilepton_CADaughters_ind_reduced","Decay_Chain(MC_dilepton_CADaughters_ind)")
    rdf = rdf.Define("MC_dilepton_CADaughters_reduced","Translate_PDG(MC_dilepton_CADaughters_ind_reduced,MC_PDG)")
    rdf = rdf.Define("MC_dilepton_CADaughters_reduced_PDG","Plain_PDG(MC_dilepton_CADaughters_ind_reduced,MC_PDG)")
    
    rdf = rdf.Define("MC_dileptonplus_Family","Find_MuSubdecays(MC_dilepton_CADaughters_ind_reduced,MC_PDG,+1)")
    rdf = rdf.Define("MC_dileptonminus_Family","Find_MuSubdecays(MC_dilepton_CADaughters_ind_reduced,MC_PDG,-1)")
    rdf = rdf.Define("MC_dilepton_BkgCat","10*Find_Categories(MC_dileptonplus_Family) + Find_Categories(MC_dileptonminus_Family)")

    return rdf



#=========================================================================

rdf = Load_RDF("bb")
print("===== B->D->K =====\n")
rdf.Filter("MC_dilepton_BkgCat == 12 || MC_dilepton_BkgCat == 21").Display("MC_dilepton_CADaughters_reduced_PDG",30).Print()
print("\n===== Unkown ======\n")
rdf.Filter("MC_dilepton_BkgCat > 50").Display("MC_dilepton_CADaughters_reduced_PDG",30).Print()

heights, edges = np.histogram(rdf.AsNumpy(["MC_dilepton_BkgCat"])["MC_dilepton_BkgCat"],bins=61,range=(-0.5,60.5))

fig, ax = plt.subplots()

ax.stairs(heights,edges)

fig.savefig("BKGCat.pdf")

#c = r.TCanvas("c","c")
#h = rdf.Histo1D("MC_dimuon_BkgCat")
#h.Draw()
#c.SaveAs("BkgCat.pdf")

#for Mode in ["bb","cc","ss","ud","sig"]:
#    rdf = Load_RDF(Mode)
#
#    #rdf.Display(["MC_dimuon_CADaughters_reduced"],30).Print()
#
#    h = rdf.Histo1D("MC_dimuon_BkgCat")
#    h.Draw()
    
