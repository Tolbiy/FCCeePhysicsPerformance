import ROOT as r
from termcolor import colored
import json
import numpy as np

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
        
    Path = "/eos/experiment/fcc/ee/analyses_storage/flavor/Bs2TauTau/flatNtuples/winter2023/analysis_stage1_Leptons_noFilter/"

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
    #rdf = rdf.Define("has_good_dimuon","good_dimuon(dimuon_ind, muon_OpeningAngle, muon_charge, muon_thrustangles)")
    #rdf = rdf.Define("has_dimuon","Check_dimuon_Presence(muon_OpeningAngle)")
    #rdf = rdf.Define("has_dimuon_SameSide","Check_dimuon_SameSide(muon_OpeningAngle)")
    #rdf = rdf.Define("has_dimuon_SigHemi","Check_dimuon_SigHemi(dimuon_ind, muon_thrustangles)")
    #rdf = rdf.Define("has_dimuon_OppositeCharges","Check_dimuon_Charges(dimuon_ind, muon_charge)")

    return rdf

#=======================================================================================

#Storage Values
Ntot = {}
NCut = {}

NTot = 0
NCuts = 0

#The Cut list to be studied
CutList = ["n_lepton > 1",
           "EVT_ThrustEmin_E < 38",
           "recoEmiss_e > 10",
           "EVT_ThrustEmin_Eneutral < 10",

           #dimuon system requirement
           "has_dilepton > 0", #/!\ Always required if want to study any other dimuon cuts
           "has_dilepton_SameSide > 0",
           "has_dilepton_SigHemi > 0",
           "has_dilepton_OppositeCharges > 0",
           "has_dilepton_Vertex == 0",
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

     
            
