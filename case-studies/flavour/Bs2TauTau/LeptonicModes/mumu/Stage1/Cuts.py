import ROOT as r
from termcolor import colored
import json
import numpy as np

#Function to find the smallest (and thus biggest cosine) opening angle -> chosen has the dimuon coming from Bs2TauTau
r.gInterpreter.Declare('''
    float GetMin(ROOT::VecOps::RVec<float> list){
        if (list.size() == 0) return -2.0;
        else if (list.size() == 1) return list[0];
        else{
            float result (list[0]);
            for (size_t i=1; i<list.size();++i){
                if (result < list[i]) result = list[i];
            }
            return result;
        }
    }
''')


#Function that finds the dimuon system based on muon opening and thrust angles, hemisphere emission and charges
#Explicitly we need :
# - An opening angle whose cosine > 0 (emitted on the same side)
# - The thrust angles of each muon needs to be positive (emitted in the signal hemisphere)
# - The two muons need to be of opposite charge
r.gInterpreter.Declare('''
    //All in one function
    int good_dimuon(ROOT::VecOps::RVec<int> dimuon_ind,
                    ROOT::VecOps::RVec<float> muon_OAs,
                    ROOT::VecOps::RVec<int> muon_charges,
                    ROOT::VecOps::RVec<float> muon_thrustangles){

        if (muon_OAs.size() == 1 && std::abs(muon_OAs[0]+2.0) < 1e-4) return 0;
        else{
            if (GetMin(muon_OAs) < 0.0) return 0;
            else {
                if (muon_thrustangles[dimuon_ind.at(0)] < 0.0 && muon_thrustangles[dimuon_ind.at(1)] < 0.0) return 0;
                else {
                    if (muon_charges[dimuon_ind.at(0)]*muon_charges[dimuon_ind.at(1)] > 0) return 0;
                    else return 1;
                }
            }
        }
    }

    //Split the cuts to study them individually

    int Check_dimuon_Presence(ROOT::VecOps::RVec<float> muon_OAs){
        if (muon_OAs.size() == 1 && std::abs(muon_OAs[0]+2.0) < 1e-4) return 0;
        else return 1;
    }

    int Check_dimuon_SameSide(ROOT::VecOps::RVec<float> muon_OAs){
        if (GetMin(muon_OAs) < 0.0) return 0;
        else return 1;
    }

    int Check_dimuon_SigHemi(ROOT::VecOps::RVec<int> dimuon_ind, ROOT::VecOps::RVec<float> muon_thrustangles){
        if (muon_thrustangles[dimuon_ind.at(0)] < 0.0 && muon_thrustangles[dimuon_ind.at(1)] < 0.0) return 0;
        else return 1;
    }

    int Check_dimuon_Charges(ROOT::VecOps::RVec<int> dimuon_ind, ROOT::VecOps::RVec<int> muon_charges){
        if (muon_charges[dimuon_ind.at(0)]*muon_charges[dimuon_ind.at(1)] > 0) return 0;
        else return 1;
    }
''')

#Create the cut argument from the cut list --------------------------------

def Build_FilterArg(CutList):
    
    FullList = ""
    
    if len(CutList) == 1:
        FullList = CutList[0]
    
    else:
        for Cut in CutList:
            if Cut == CutList[0]:
                FullList = Cut
            else:
                FullList += f" && {Cut}"

    return FullList


#Get the histograms in a numpy form -----------------------------------
def Load_RDF(Mode):
        
    Path = "/eos/experiment/fcc/ee/analyses_storage/flavor/Bs2TauTau/flatNtuples/winter2023/analysis_stage1_mumu_noFilter/"

    Links = {"bb":"p8_ee_Zbb_ecm91",
             "cc":"p8_ee_Zcc_ecm91",
             "ss":"p8_ee_Zss_ecm91",
             "ud":"p8_ee_Zcc_ecm91",
             "sig":"p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau",
            }

    #Load the rdf properly (avoid segfault from Tree going out of scope)
    filenames = r.std.vector('string')()
    if Mode == "sig":
        for i in np.arange(0,20,1):
            filenames.push_back(Path+Links[Mode]+f"/chunk_{i}.root")
    else:
        for i in np.arange(0,10,1):
            filenames.push_back(Path+Links[Mode]+f"/chunk_{i}.root")
    rdf = r.RDataFrame("events",filenames)
    rdf = rdf.Define("has_good_dimuon","good_dimuon(dimuon_ind, muon_OpeningAngle, muon_charge, muon_thrustangles)")
    rdf = rdf.Define("has_dimuon","Check_dimuon_Presence(muon_OpeningAngle)")
    rdf = rdf.Define("has_dimuon_SameSide","Check_dimuon_SameSide(muon_OpeningAngle)")
    rdf = rdf.Define("has_dimuon_SigHemi","Check_dimuon_SigHemi(dimuon_ind, muon_thrustangles)")
    rdf = rdf.Define("has_dimuon_OppositeCharges","Check_dimuon_Charges(dimuon_ind, muon_charge)")

    return rdf

#=======================================================================================

#Storage Values
Ntot = {}
NCut = {}

NTot = 0
NCuts = 0

#The Cut list to be studied
CutList = ["n_muons > 1",
           "EVT_ThrustEmin_E < 38",
           "recoEmiss_e > 10",
           "EVT_ThrustEmin_Eneutral < 10",

           #dimuon system requirement
           "has_dimuon > 0", #/!\ Always required if want to study any other dimuon cuts
           "has_dimuon_SameSide > 0",
           #"has_dimuon_SigHemi > 0",
           "has_dimuon_OppositeCharges > 0"
          ] 


#prevents adding same lines in case calling multiple time the scipt 
Previous = {}
AlreadyDone = False
with open("Cuts.json","r") as ofile:
    Previous = json.load(ofile)
    
    for Cuts in Previous.keys():
        if Build_FilterArg(CutList) == Cuts:
            AlreadyDone = True
            print("Already studied")
            exit()

print("\n== Cuts studied ==\n")
for cut in CutList:
    print(f"{cut}")
print("\n")

for mode in ["bb","cc","ss","ud","sig"]:

    print(f"_____________ {mode} _______________")

    #Load the number of cuts before any cuts
    Ntot[mode] = Previous["Before"][mode]
    NTot += Ntot[mode]

    rdf = Load_RDF(mode)
    rdf = rdf.Filter(Build_FilterArg(CutList))
    NCut[mode] = rdf.Count().GetValue()
    NCuts += NCut[mode]
    
    #Efficiencies of how many got rejected so red should be close to 100% and blue close to 0%
    print(f"Cuts eff = {round(100*NCut[mode]/Ntot[mode],5)} %")

print("\n-----------------------------------------")
print(f"=====> Total eff = {round(100*NCuts/NTot,5)} %")
print("-----------------------------------------\n")


#Store these numbers
with open("Cuts.json","w") as ofile:
    Previous[Build_FilterArg(CutList)] = NCut
    json.dump(Previous, ofile)

     
            
