import ROOT as r
import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amssymb}'

#Dict for shortened dict key(global var -> put it in a config file)
Links = {"bb": "p8_ee_Zbb_ecm91",
         "cc": "p8_ee_Zcc_ecm91",
         "ss": "p8_ee_Zss_ecm91",
         "ud": "p8_ee_Zud_ecm91",
         "sig":"p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau",
         "exc":"p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU"
        }

#Get the model to be added to the dataframe
r.gInterpreter.ProcessLine('''
TMVA::Experimental::RBDT<> bdt("Naive_withMoreData", "/afs/cern.ch/work/t/tomonnar/public/FCCeePhysicsPerformance/case-studies/flavour/Bs2TauTau/LeptonicModes/merged/Stage2/BDT/TrainTest/Train_Results/Models/Naive_withMoreData.root");
computeModel1 = TMVA::Experimental::Compute<44, float>(bdt);
''')

#Get the histograms in a numpy form -----------------------------------
def Load_AsNumpy(Mode,var,Nbins,withEdge):
        
    Path = "/eos/experiment/fcc/ee/analyses_storage/flavor/Bs2TauTau/flatNtuples/winter2023/analysis_stage1_Leptons_noFilter/"

    #Load the rdf properly (avoid segfault from Tree going out of scope)
    filenames = r.std.vector('string')()
    if Mode == "sig":
        for i in np.arange(0,20,1):
            filenames.push_back(Path+Links[Mode]+f"/chunk_{i}.root")
    else:
        for i in np.arange(0,10,1):
            filenames.push_back(Path+Links[Mode]+f"/chunk_{i}.root")
    rdf = r.RDataFrame("events",filenames)
    print(f"{Mode}: Number of events = {rdf.Count().GetValue()}")
    rdf = rdf.Filter("n_lepton > 1 && EVT_ThrustEmin_E < 38 && recoEmiss_e > 10 && EVT_ThrustEmin_Eneutral < 10 && has_dilepton > 0 && has_dilepton_SameSide > 0 && has_dilepton_SigHemi > 0 && has_dilepton_OppositeCharges > 0 && has_dilepton_Vertex == 0")
    print(f"{Mode}: Number of selected events = {rdf.Count().GetValue()}")
    

    for vari in ["px","py","pz","phi","eta","energy","mass","charge","PDG","thrustangles"]:
        rdf = rdf.Define(f"plus_{vari}",f"lepton_{vari}.at(dilepton_plus_ind)").Define(f"minus_{vari}",f"lepton_{vari}.at(dilepton_minus_ind)")
    rdf = rdf.Define("Opening_Angle","(plus_px*minus_px+plus_py*minus_py+plus_pz*minus_pz)/sqrt(plus_px*plus_px+plus_py*plus_py+plus_pz*plus_pz)/sqrt(minus_px*minus_px+minus_py*minus_py+minus_pz*minus_pz)")
    rdf = rdf.Define("EAsymm","(EVT_ThrustEmax_E-EVT_ThrustEmin_E)/(EVT_ThrustEmin_E+EVT_ThrustEmax_E)")
    rdf = rdf.Redefine("dilepton_case","float(dilepton_case)")

    #Add BDT
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

    rdf = rdf.Define("MVAVec2", r.computeModel1,BDTvars).Define("MVA2_withMoreData", "MVAVec2.at(0)")


    ##### To Add Any RDF Treatment Needed #####

    ###########################################

    hist, edges = np.histogram(rdf.AsNumpy([var])[var],bins=Nbins[0],range=(Nbins[1],Nbins[2]))

    if withEdge:
        return hist, edges
    else:
        return hist

#Draw the resolution sig + bkg ------------------------------
def Plot(hlist,edges,var):

    #To normalise and compare
    Norm = {}
    for mode in hlist.keys():
        Norm[mode] = np.sum(hlist[mode])

    fig, ax = plt.subplots()
    ax.grid(color='grey', linestyle='--', linewidth=0.5, alpha=0.5)

    ax.stairs(hlist["ud"]/Norm["bkg"] + hlist["ss"]/Norm["bkg"] + hlist["cc"]/Norm["bkg"] + hlist["bb"]/Norm["bkg"],  edges,lw=0,fill=True,color='slategrey',label=r"$Z^0\rightarrow u\overline{d}$")
    ax.stairs(hlist["ss"]/Norm["bkg"] + hlist["cc"]/Norm["bkg"] + hlist["bb"]/Norm["bkg"],                            edges,lw=0,fill=True,color='mediumseagreen',label=r"$Z^0\rightarrow s\overline{s}$")
    ax.stairs(hlist["cc"]/Norm["bkg"] + hlist["bb"]/Norm["bkg"],                                                      edges,lw=0,fill=True,color='goldenrod',label=r"$Z^0\rightarrow c\overline{c}$")
    ax.stairs(hlist["bb"]/Norm["bkg"],                                                                                edges,lw=0,fill=True,color='steelblue',label=r"$Z^0\rightarrow b\overline{b}$")
    
    ax.stairs(hlist["bkg"]/Norm["bkg"],edges,ec='black',ls='-',lw=2,label=r"$\textrm{Total background}$")
    #ax.stairs(hlist["exc"]/Norm["exc"],edges,ec='maroon',ls='--',lw=2,label=r"$\textrm{Signal (Exclusive)}$")
    ax.stairs(hlist["sig"]/Norm["sig"],edges,ec='red',ls='-',lw=2,label=r"$\textrm{Signal}$")

    #Set axis range
    ax.set_xlim([edges[0],edges[-1]])

    #Labels
    ax.set_xlabel(r"$\textrm{"+f"{var}"+r"}$",size="large")
    ax.set_ylabel(r"$\textrm{Normalised events}$",size="large")

    #Legend (Count the number of bkg modes drawn to write it down properly)
    ax.legend()
        
    #Save
    fig.savefig(f'{var}.pdf')

#=======================================================================================================

def Do_Plots(var,bins):

    for v in var:
        hmode = {}
        hmode["sig"], edges = Load_AsNumpy("sig",v,bins,True)
        for mode in ["bb","cc","ss","ud"]:
            hmode[mode] = Load_AsNumpy(mode,v,bins,False)
        hmode["bkg"] = hmode["bb"]+hmode["cc"]+hmode["ss"]+hmode["ud"]
        Plot(hmode,edges,v)


#=======================================================================================================

parser = argparse.ArgumentParser()
parser.add_argument("NBins", type=int)
parser.add_argument("Low",   type=float)
parser.add_argument("High",  type=float)
parser.add_argument("List",   nargs='+')
args=parser.parse_args()
Do_Plots(args.List,(args.NBins,args.Low,args.High))