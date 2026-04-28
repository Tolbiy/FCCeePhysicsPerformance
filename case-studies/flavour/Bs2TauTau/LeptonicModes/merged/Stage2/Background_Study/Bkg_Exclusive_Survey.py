import ROOT as r
from termcolor import colored
import json
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amssymb}'

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
            else if (5000 < std::abs(MuFamily.at(0)) && std::abs(MuFamily.at(0)) < 6000) return 5; //b-baryon mother
            else if (4000 < std::abs(MuFamily.at(0)) && std::abs(MuFamily.at(0)) < 5000) return 6; //c-baryon mother
            else return 7;
        }

    }

    int Find_SuperCat(int Cat){
        if (Cat == 12 or Cat == 21) return 0; //B->D->K 
        else if (Cat == 56 or Cat == 65) return 1; //Lb->Lc->L 
        else if (Cat == 14 or Cat == 24 or Cat == 41 or Cat == 42) return 2; //B->TaunuD or D->TaunuK
        else if (Cat == 54 or Cat == 64 or Cat == 45 or Cat == 46) return 3; //Lb->TaunuLc or Lc->TaunuL
        else if (Cat == 22) return 4; //B->DD
        else if (Cat == 44) return 5; //B->TauD
        else if (Cat == 0) return 6; //Independant B/D
        else return 7; //Code effect or UNKNOWN
    }

    string NOCA_case_PDG(ROOT::VecOps::RVec<int> Mother_ind, ROOT::VecOps::RVec<int> MCPDG){

        string result ("");

        if (Mother_ind.at(0) == -1) result = "Has CA or no dilepton";
        else if (Mother_ind.at(0) == -3) result = "No or multiparents";
        else {
            for (size_t i=0; i<Mother_ind.size(); ++i){
                if (Mother_ind.at(i) == -2) break;
                else result += std::to_string(MCPDG.at(Mother_ind.at(i)))+"/";
            }
        }
        return result;
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
    rdf = rdf.Define("MC_dilepton_nFamily","MC_dileptonplus_Family.size()+MC_dileptonminus_Family.size()")
    #print(f"Total amount = {rdf.Count().GetValue()}")
    #rdf = rdf.Filter("MC_dilepton_nFamily == 8")
    #print(f"No additional lepton amount = {rdf.Count().GetValue()}")
    rdf = rdf.Define("MC_dilepton_BkgCat","10*Find_Categories(MC_dileptonplus_Family) + Find_Categories(MC_dileptonminus_Family)")
    rdf = rdf.Define("MC_dilepton_BkgSupCat","Find_SuperCat(MC_dilepton_BkgCat)")

    rdf = rdf.Define("MC_dilep1_NoCAMother_PDG","NOCA_case_PDG(MC_dilep1_NoCADaughters_ind,MC_PDG)")
    rdf = rdf.Define("MC_dilep2_NoCAMother_PDG","NOCA_case_PDG(MC_dilep2_NoCADaughters_ind,MC_PDG)")
    
    rdf = rdf.Define("lp_ta","lepton_thrustangles.at(dilepton_plus_ind)")
    rdf = rdf.Define("lm_ta","lepton_thrustangles.at(dilepton_minus_ind)")

    return rdf



#=========================================================================

rdf = Load_RDF("bb")

print("\n===== B->D->K =====\n")
rdf.Filter("MC_dilepton_BkgCat == 12 || MC_dilepton_BkgCat == 21").Display("MC_dilepton_CADaughters_reduced_PDG",30).Print()

print("\n===== Tau semileptonic (mesonic) =====\n")
rdf.Filter("MC_dilepton_BkgCat == 14 || MC_dilepton_BkgCat == 24 || MC_dilepton_BkgCat == 41 || MC_dilepton_BkgCat == 42").Display("MC_dilepton_CADaughters_reduced_PDG",30).Print()

print("\n===== Double D =====\n")
rdf.Filter("MC_dilepton_BkgCat == 22").Display("MC_dilepton_CADaughters_reduced_PDG",30).Print()

print("\n===== Double B =====\n")
rdf.Filter("MC_dilepton_BkgCat == 11").Display("MC_dilepton_CADaughters_reduced_PDG",30).Print()

print("\n===== Lb->Lc->L =====\n")
rdf.Filter("MC_dilepton_BkgCat == 56 || MC_dilepton_BkgCat == 65").Display("MC_dilepton_CADaughters_reduced_PDG",30).Print()

print("\n===== Tau semileptonic (baryonic) =====\n")
rdf.Filter("MC_dilepton_BkgCat == 54 || MC_dilepton_BkgCat == 64 || MC_dilepton_BkgCat == 45 || MC_dilepton_BkgCat == 46").Display("MC_dilepton_CADaughters_reduced_PDG",30).Print()

print("\n===== Double Lc =====\n")
rdf.Filter("MC_dilepton_BkgCat == 66").Display("MC_dilepton_CADaughters_reduced_PDG",30).Print()

print("\n===== Double Lb =====\n")
rdf.Filter("MC_dilepton_BkgCat == 55").Display("MC_dilepton_CADaughters_reduced_PDG",30).Print()

print("\n===== Double Tau =====\n")
rdf.Filter("MC_dilepton_BkgCat == 44").Display("MC_dilepton_CADaughters_reduced_PDG",30).Print()

print("\n===== No Common Ancestor ======\n")
rdf.Filter("MC_dilepton_BkgCat == 0").Display(["MC_dilep1_NoCAMother_PDG","MC_dilep2_NoCAMother_PDG"],30).Print()
rdf.Filter("MC_dilepton_BkgCat == 0").Display(["MC_dilep1_NoCADaughters_ind","MC_dilep2_NoCADaughters_ind"],30).Print()
rdf.Filter("MC_dilepton_BkgCat == 0").Display(["lp_ta","lm_ta"],30).Print()

print("\n===== Unkown ======\n")
rdf.Filter("MC_dilepton_BkgCat > 70 || MC_dilepton_BkgCat == 17 || MC_dilepton_BkgCat == 27 || MC_dilepton_BkgCat == 37 || MC_dilepton_BkgCat == 47 || MC_dilepton_BkgCat == 57 || MC_dilepton_BkgCat == 67").Display("MC_dilepton_CADaughters_reduced_PDG",30).Print()

Tot = rdf.Count().GetValue()
heights, edges = np.histogram(rdf.AsNumpy(["MC_dilepton_BkgSupCat"])["MC_dilepton_BkgSupCat"],bins=8,range=(-0.5,7.5))
heights2, edges2 = np.histogram(rdf.AsNumpy(["MC_dilepton_BkgCat"])["MC_dilepton_BkgCat"],bins=77,range=(-0.5,76.5))
heights = 100*heights/Tot

fig, ax = plt.subplots()
ax.stairs(heights2,edges2)
ax.set_xlim([-1,77])
#ax.set_ylim([0,10])
fig.savefig("BKGCat.pdf")

fig2, ax2 = plt.subplots()
#Decays = [
#    r"$\boldsymbol{B}\rightarrow \boldsymbol{D}\rightarrow K$",
#    r"$\boldsymbol{\Lambda_b}\rightarrow\boldsymbol{\Lambda_c}\rightarrow{\Lambda}$",
#    r"$B/D\rightarrow\boldsymbol{\tau}\nu_{\tau} \boldsymbol{D/K}$",
#    r"$B\rightarrow \boldsymbol{DD}$",
#    r"$B\rightarrow D(\rightarrow\boldsymbol{\tau}\nu_{\tau}K)\boldsymbol{\tau}$",
#    r"\textrm{Independent }$\boldsymbol{B}/\boldsymbol{D}$",
#    r"TBD",
#]

heights_BM = np.array([heights[0]+heights[1],heights[2]+heights[3],heights[4],heights[5],heights[6],heights[7]])
heights_M = np.array([heights[0],heights[2],heights[4],heights[5],heights[6],heights[7]])

Decaysn = [
    r"\textrm{Cascade without }$\tau$",
    r"\textrm{Cascade with }$\tau$",
    r"\textrm{Double }$D$",
    r"\textrm{Double }$\tau$",
    r"\textrm{Independent}",
    r"TBD",
]

ax2.grid(color='grey', linestyle='--', linewidth=0.5, alpha=0.5)
ax2.bar(Decaysn, heights_BM, color="lightsteelblue",label=r"\textrm{Baryonic}")
ax2.bar(Decaysn, heights_M, color="steelblue",label=r"\textrm{Mesonic}")
ax2.set_ylabel(r"\textrm{Ratio (out of }$"+f"{Tot}"+r"$\textrm{ events) [}$\%$\textrm{]}",size="x-large")
ax2.set_xticks(np.arange(0,6,1),labels=Decaysn,ha="right",size="large")
plt.xticks(rotation=35)
ax2.set_xlim([-0.5,4.5])

ax2.legend(fontsize="x-large")

fig2.tight_layout()
fig2.savefig("SupCat.pdf")

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
    
