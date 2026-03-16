import ROOT as r
import numpy as np
import sys

#sys.path.insert(1,"../Utils")
#import Functions_ForPreTreat

Links = {"bb": "p8_ee_Zbb_ecm91",
         "cc": "p8_ee_Zcc_ecm91",
         "ss": "p8_ee_Zss_ecm91",
         "ud": "p8_ee_Zud_ecm91",
         "sig":"p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau"
        }

def RDF_Treatment(rdf):

    for var in ["px","py","pz","phi","eta","energy","mass","charge","PDG","thrustangles"]:
        rdf = rdf.Define(f"plus_{var}",f"lepton_{var}.at(dilepton_plus_ind)").Define(f"minus_{var}",f"lepton_{var}.at(dilepton_minus_ind)")
    rdf = rdf.Define("Opening_Angle","plus_px*minus_px+plus_py*minus_py+plus_pz*minus_pz/sqrt(plus_px*plus_px+plus_py*plus_py+plus_pz*plus_pz)/sqrt(minus_px*minus_px+minus_py*minus_py+minus_pz*minus_pz)")

    return rdf

    

def Load_RDF(mode,amount):

    Path = "/eos/experiment/fcc/ee/analyses_storage/flavor/Bs2TauTau/flatNtuples/winter2023/analysis_stage1_Leptons_withCuts/"

    if mode == "sig": 
        
        NF = 10
        filenames = r.std.vector('string')()
        for chunk in np.arange(0,NF,1):
            filenames.push_back(Path+Links[mode]+f"/chunk_{chunk}.root")
        rdf = r.RDataFrame("events",filenames)
        rdf = RDF_Treatment(rdf)
        return rdf
    
    else: 
        
        NF = 10
        filenames = r.std.vector('string')()
        for chunk in np.arange(0,NF,1):
            filenames.push_back(Path+Links[mode]+f"/chunk_{chunk}.root")
        rdf = r.RDataFrame("events",filenames)
        rdf = RDF_Treatment(rdf)
        return rdf


def Do_RDF_PreTreatment(amount):

    VarList = ["plus_px","plus_py","plus_pz","plus_phi","plus_eta","plus_energy","plus_mass","plus_thrustangles",
               "minus_px","minus_py","minus_pz","minus_phi","minus_eta","minus_energy","minus_mass","minus_thrustangles",
               "Opening_Angle","dilepton_case"
               ]

    print(f"Variables selected for training the BDT:\n{VarList}")

    columns = r.std.vector('string')()
    for var in VarList:
        columns.push_back(var)

    RDFs = {}
    for mode in ["sig","bb"]:
        RDFs[mode] = Load_RDF(mode,amount)
        RDFs[mode].Snapshot("events",f"PreTreatedData/{mode}/Naive.root",columns)

#=========================================================================

Do_RDF_PreTreatment(1)

