import ROOT as r
import numpy as np
import matplotlib.pyplot as plt
import json
import sys

from matplotlib import rc
rc('font',**{'family':'serif','serif':['Roman']})
rc('text', usetex=True)

Links = {"bb": "p8_ee_Zbb_ecm91",
         "cc": "p8_ee_Zcc_ecm91",
         "ss": "p8_ee_Zss_ecm91",
         "ud": "p8_ee_Zud_ecm91",
         "sig":"p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau"
        }

#----------------------
#----------------------

def Add_BDT(rdf,BDTModel,varlist):

    for vari in ["px","py","pz","phi","eta","energy","mass","charge","PDG","thrustangles"]:
        rdf = rdf.Define(f"plus_{vari}",f"lepton_{vari}.at(dilepton_plus_ind)").Define(f"minus_{vari}",f"lepton_{vari}.at(dilepton_minus_ind)")
    rdf = rdf.Define("Opening_Angle","(plus_px*minus_px+plus_py*minus_py+plus_pz*minus_pz)/sqrt(plus_px*plus_px+plus_py*plus_py+plus_pz*plus_pz)/sqrt(minus_px*minus_px+minus_py*minus_py+minus_pz*minus_pz)")
    rdf = rdf.Define("EAsymm","(EVT_ThrustEmax_E-EVT_ThrustEmin_E)/(EVT_ThrustEmin_E+EVT_ThrustEmax_E)")
    rdf = rdf.Redefine("dilepton_case","float(dilepton_case)")

    rdf = rdf.Define("MVAVec2", r.computeModel1,varlist).Define("MVA2", "MVAVec2.at(0)")

    return rdf

#----------------------

def Load_Snap_RDF(mode,BDTModel,varlist):

    Path = "/eos/experiment/fcc/ee/analyses_storage/flavor/Bs2TauTau/flatNtuples/winter2023/analysis_stage1_Leptons_withCuts/"
        
    NF = 100
    if mode == "sig": 
        NF = 20
    
    filenames = r.std.vector('string')()
    for chunk in np.arange(0,NF,1):
        filenames.push_back(Path+Links[mode]+f"/chunk_{chunk}.root")
    rdf = r.RDataFrame("events",filenames)
    rdf = Add_BDT(rdf,BDTModel,varlist)

    columns = r.std.vector('string')()
    for var in varlist:
        columns.push_back(var)
    columns.push_back("MVA2")
    rdf.Snapshot("events",f"{mode}/{BDTModel}.root",columns)

#----------------------

def Do_AddBDT(BDTName,VarList):

    for mode in ["sig","bb","cc"]: #Restrained to these while waiting for the rest of the samples to be produced
        print(f"Loading and Saving {mode} mode")
        Load_Snap_RDF(mode,BDTName,VarList)


#========================================================================

#BDT columns
ModelsDict = {"Naive": [#Muon related
                        "plus_px","plus_py","plus_pz","plus_phi","plus_eta","plus_energy","plus_mass","plus_thrustangles",
                        "minus_px","minus_py","minus_pz","minus_phi","minus_eta","minus_energy","minus_mass","minus_thrustangles","Opening_Angle","dilepton_case",

                        #Event level
                        "EVT_ThrustEmax_E","EVT_ThrustEmin_E","EVT_ThrustEmax_Echarged","EVT_ThrustEmin_Echarged","EVT_ThrustEmax_Eneutral","EVT_ThrustEmin_Eneutral",
                        "EVT_ThrustEmax_N","EVT_ThrustEmin_N","EVT_ThrustEmax_Ncharged","EVT_ThrustEmin_Ncharged","EVT_ThrustEmax_Nneutral","EVT_ThrustEmin_Nneutral",
                        "recoEmiss_thrustangle","recoEmiss_e","EVT_Thrust_Mag","EVT_Thrust_X","EVT_Thrust_Y","EVT_Thrust_Z","EAsymm",

                        #Vertex related
                        "EVT_ThrustEmin_NDV","EVT_ThrustEmax_NDV","EVT_dPV2DVmin","EVT_dPV2DVmax","EVT_dPV2DVave",
                        "EVT_NtracksPV","EVT_NVertex",
                       ],
              "Naive_withMoreData": [#Muon related
                        "plus_px","plus_py","plus_pz","plus_phi","plus_eta","plus_energy","plus_mass","plus_thrustangles",
                        "minus_px","minus_py","minus_pz","minus_phi","minus_eta","minus_energy","minus_mass","minus_thrustangles","Opening_Angle","dilepton_case",

                        #Event level
                        "EVT_ThrustEmax_E","EVT_ThrustEmin_E","EVT_ThrustEmax_Echarged","EVT_ThrustEmin_Echarged","EVT_ThrustEmax_Eneutral","EVT_ThrustEmin_Eneutral",
                        "EVT_ThrustEmax_N","EVT_ThrustEmin_N","EVT_ThrustEmax_Ncharged","EVT_ThrustEmin_Ncharged","EVT_ThrustEmax_Nneutral","EVT_ThrustEmin_Nneutral",
                        "recoEmiss_thrustangle","recoEmiss_e","EVT_Thrust_Mag","EVT_Thrust_X","EVT_Thrust_Y","EVT_Thrust_Z","EAsymm",

                        #Vertex related
                        "EVT_ThrustEmin_NDV","EVT_ThrustEmax_NDV","EVT_dPV2DVmin","EVT_dPV2DVmax","EVT_dPV2DVave",
                        "EVT_NtracksPV","EVT_NVertex",
                       ]
}


BDT_Name = "Naive_withMoreData"
#Get the model to be added to the dataframe
r.gInterpreter.ProcessLine(f'''
    TMVA::Experimental::RBDT<> bdt("{BDT_Name}", "/afs/cern.ch/work/t/tomonnar/public/FCCeePhysicsPerformance/case-studies/flavour/Bs2TauTau/LeptonicModes/merged/Stage2/BDT/TrainTest/Train_Results/Models/{BDT_Name}.root");
    computeModel1 = TMVA::Experimental::Compute<44, float>(bdt);
''')
Do_AddBDT(BDT_Name,ModelsDict[BDT_Name])