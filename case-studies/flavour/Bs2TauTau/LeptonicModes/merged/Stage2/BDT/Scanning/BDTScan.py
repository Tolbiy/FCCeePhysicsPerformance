import ROOT as r
import numpy as np
import matplotlib.pyplot as plt
import json
import sys

from matplotlib import rc
rc('font',**{'family':'serif','serif':['Roman']})
rc('text', usetex=True)

#Get the model to be added to the dataframe
r.gInterpreter.ProcessLine('''
TMVA::Experimental::RBDT<> bdt("Naive", "/afs/cern.ch/work/t/tomonnar/public/FCCeePhysicsPerformance/case-studies/flavour/Bs2TauTau/LeptonicModes/merged/Stage2/BDT/TrainTest/Train_Results/Models/Naive.root");
computeModel1 = TMVA::Experimental::Compute<44, float>(bdt);
''')


Links = {"bb": "p8_ee_Zbb_ecm91",
         "cc": "p8_ee_Zcc_ecm91",
         "ss": "p8_ee_Zss_ecm91",
         "ud": "p8_ee_Zud_ecm91",
         "sig":"p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau"
        }

def Add_BDT(rdf):

    for var in ["px","py","pz","phi","eta","energy","mass","charge","PDG","thrustangles"]:
        rdf = rdf.Define(f"plus_{var}",f"lepton_{var}.at(dilepton_plus_ind)").Define(f"minus_{var}",f"lepton_{var}.at(dilepton_minus_ind)")
    rdf = rdf.Define("Opening_Angle","plus_px*minus_px+plus_py*minus_py+plus_pz*minus_pz/sqrt(plus_px*plus_px+plus_py*plus_py+plus_pz*plus_pz)/sqrt(minus_px*minus_px+minus_py*minus_py+minus_pz*minus_pz)")
    rdf = rdf.Define("EAsymm","(EVT_ThrustEmax_E-EVT_ThrustEmin_E)/(EVT_ThrustEmin_E+EVT_ThrustEmax_E)")
    rdf = rdf.Redefine("dilepton_case","float(dilepton_case)")

    BDTvars = [#Muon related
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

    rdf = rdf.Define("MVAVec2", r.computeModel1,BDTvars).Define("MVA2", "MVAVec2.at(0)")

    return rdf


def Load_RDF(mode):

    Path = "/eos/experiment/fcc/ee/analyses_storage/flavor/Bs2TauTau/flatNtuples/winter2023/analysis_stage1_Leptons_withCuts/"
        
    NF = 100
    if mode == "sig": 
        NF = 20 
    
    filenames = r.std.vector('string')()
    for chunk in np.arange(0,NF,1):
        filenames.push_back(Path+Links[mode]+f"/chunk_{chunk}.root")
    rdf = r.RDataFrame("events",filenames)
    rdf = Add_BDT(rdf)
    return rdf

def BDT_Scan(rdf):

    After = []
    cuts = np.linspace(0.66,0.99,20)

    for cut in cuts:

        rdf = rdf.Filter(f"MVA2 > {cut}")
        After.append(rdf.Count().GetValue())
        if int(100*cut)%25 == 0:
            print(f"  Treated {int(100*cut)}% of the scan")

    return After


def Do_Scan():

    Afters = {}
    Befores = {}
    Eff = {}
    for mode in ["bb","cc","sig"]:
        print(f"->Loading {mode}")
        rdf = Load_RDF(mode)
        print("  Getting Before and After BDT cuts numbers")
        Befores[mode] = rdf.Count().GetValue()
        Afters[mode] = BDT_Scan(rdf)
        Eff[mode] = [Afters[mode][i]/Befores[mode] for i in np.arange(0,len(Afters[mode]),1)]

    print(Afters)
    print(Befores)
    print(Eff)

    #Store these numbers
    with open("Naive_66_99.json","w") as ofile:
        json.dump({"Before":Befores,"After_Scan":Afters,"Efficiencies":Eff}, ofile)


def Do_Plots():

    with open("Naive_66_99.json","r") as ofile:
        BDT_Eff = json.load(ofile)

    NZ = 6e12 #Number of Z (To be refined depending on FCC scenarios)
    Bs_had = 0.1 #Ratio of b->Bs hadronisation (given by Xunwu, to be refined)
    Bs2TauTau_BR = 7.73e-7 #Highest SM estimate for Bs->TauTau (from the LHCb paper that searched for it)
    Tau23Pi = 0.0931 #Tau->3pi estimate (PDG)
    Tau2l_PDG = 0.1737**2 + 0.1785**2 + 2*0.1737*0.1785 #Tau->lnunu from PDG values
    Tau2l_Sim = 1250154/10000000 #Tau->lnunu from sim, keep in mind the num here includes the hasPV

    #BR of Z->qq
    BR_qq = {}
    BR_qq["bb"] = 0.1512
    BR_qq["cc"] = 0.1203
    BR_qq["ss"] = 0.1560
    BR_qq["ud"] = 0.6991 - BR_qq["bb"] - BR_qq["cc"] - BR_qq["ss"]

    #------------------------------------

    #bkg eff preBDT (hasPV, Stage1 cuts)
    PreBDTEff = {}
    PreBDTEff["bb"] = 2105721/438738637
    PreBDTEff["cc"] = 6802/499786495

    #sig eff preBDT (hasPV, Stage1 cuts from MC decay selections)
    PreBDTEff["sig"] = 522843/1250154 #Warning the sig denom is with hasPV, hence the number should be close but not exact

    #------------------------------------

    #Get the expected number of events until Stage1
    NexpS1 = {}
    NexpS1["sig"] = NZ*BR_qq["bb"]*Bs_had*Bs2TauTau_BR*Tau2l_Sim*PreBDTEff["sig"]
    NexpS1["bb"] = NZ*BR_qq["bb"]*PreBDTEff["bb"]
    NexpS1["cc"] = NZ*BR_qq["cc"]*PreBDTEff["cc"]

    #-----------------------------------

    #Find the Nexp after Stage2 by looking at the BDT scan
    NexpS2 = {}
    NexpS2["sig"] = [NexpS1["sig"]*BDT_Eff["Efficiencies"]["sig"][i] for i in np.arange(0,len(BDT_Eff["Efficiencies"]["sig"]),1)]
    NexpS2["bb"] = [NexpS1["bb"]*BDT_Eff["Efficiencies"]["bb"][i] for i in np.arange(0,len(BDT_Eff["Efficiencies"]["bb"]),1)]
    NexpS2["cc"] = [NexpS1["cc"]*BDT_Eff["Efficiencies"]["cc"][i] for i in np.arange(0,len(BDT_Eff["Efficiencies"]["cc"]),1)]
    
    #------------------------------------

    #Add the error due to MC simulation limited numbers (background only)
    TotAfter = [BDT_Eff["After_Scan"]["bb"][i] + BDT_Eff["After_Scan"]["cc"][i] for i in np.arange(0,len(BDT_Eff["Efficiencies"]["bb"]),1)]
    TotAfterPlus = [TotAfter[i] + np.sqrt(TotAfter[i]) for i in np.arange(0,len(BDT_Eff["Efficiencies"]["bb"]),1)]
    TotAfterMinus = [TotAfter[i] - np.sqrt(TotAfter[i]) for i in np.arange(0,len(BDT_Eff["Efficiencies"]["bb"]),1)]
    TotBefore = BDT_Eff["Before"]["bb"] + BDT_Eff["Before"]["cc"]

    EffS2Plus = [TotAfterPlus[i]/TotBefore for i in np.arange(0,len(BDT_Eff["Efficiencies"]["bb"]),1)]
    EffS2Minus = [TotAfterMinus[i]/TotBefore for i in np.arange(0,len(BDT_Eff["Efficiencies"]["bb"]),1)]
    
    NexpS2Err = {}
    NexpS2Err["Plus"] = [(NexpS1["bb"] + NexpS1["cc"])*EffS2Plus[i] for i in np.arange(0,len(BDT_Eff["Efficiencies"]["bb"]),1)]
    NexpS2Err["Minus"] = [(NexpS1["bb"] + NexpS1["cc"])*EffS2Minus[i] for i in np.arange(0,len(BDT_Eff["Efficiencies"]["bb"]),1)]

    #-----------------------------------

    #Scan S/sqrt(B)
    S = np.array(NexpS2["sig"])
    
    B = np.array(NexpS2["bb"]) + np.array(NexpS2["cc"])
    BPlus = np.array(NexpS2Err["Plus"])
    BMinus = np.array(NexpS2Err["Minus"])

    print(B)
    print(BPlus)
    print(BMinus)
    
    FoM = S/np.sqrt(B)
    FoMMinus = S/np.sqrt(BPlus)
    FoMPlus = S/np.sqrt(BMinus)

    x = np.linspace(0.66,0.99,20)

    fig, ax = plt.subplots()
    ax.grid(color='grey', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.plot(x,S/NexpS1["sig"],ls="--",color="tab:blue",label=r"$\textrm{Normalised Signal Events }S$")
    ax.plot(x,B/(NexpS1["bb"]+NexpS1["cc"]),ls="--",color="tab:orange",label=r"$\textrm{Normalised Background Events }B$")
    ax.plot(x,FoM,ls="-",color="tab:purple",label=r"$\textrm{Figure of Merit: }\frac{S}{\sqrt{B}}$")
    
    print(FoM)
    print(f"Maximum FoM at Stage 2 BDT cut = {x[np.argmax(FoM[:-1])]}")

    ax.fill_between(x,FoMMinus,FoM,alpha=0.2,color="tab:purple")
    ax.fill_between(x,FoMPlus,FoM,alpha=0.2,color="tab:purple")

    ax.set_xlabel("Stage 2 BDT cuts")
    ax.set_ylabel("Normalised events")
    ax.set_xlim([x[0],x[-1]])
    ax.set_ylim([1e-4,1])
    ax.set_yscale("log")
    ax.legend()
    fig.savefig("Naive_66_99.pdf")

Do_Scan()
Do_Plots()
        
        


    


