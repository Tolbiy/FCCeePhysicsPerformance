import ROOT as r
import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amssymb}'

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

r.gInterpreter.Declare('''
    float Compare_MinMuon_TMMuon_OA(float min, float TM){
        float result;
        if (std::abs(min+2.0)<1e-4 || std::abs(TM+2.0)<1e-4) return -3.0;
        else return min-TM;
    }
''')

r.gInterpreter.Declare('''
    int Check_dimuon(ROOT::VecOps::RVec<float> OAs){
        if (OAs.size() == 1 && std::abs(OAs[0]+2.0) < 1e-4) return 0;
        else return 1;
    }
''')

#Function that finds the dimuon system based on muon opening and thrust angles, hemisphere emission and charges
#Explicitly we need :
# - An opening angle whose cosine > 0 (emitted on the same side)
# - The thrust angles of each muon needs to be positive (emitted in the signal hemisphere)
# - The two muons need to be of opposite charge
r.gInterpreter.Declare('''
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
''')


#Dict for shortened dict key(global var -> put it in a config file)
Links = {"bb": "p8_ee_Zbb_ecm91",
         "cc": "p8_ee_Zcc_ecm91",
         "ss": "p8_ee_Zss_ecm91",
         "ud": "p8_ee_Zud_ecm91",
         "sig":"p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau",
         "exc":"p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU"
        }


#Get the histograms in a numpy form -----------------------------------
def Load_AsNumpy(Mode,var,Nbins,withEdge):
        
    Path = "/eos/experiment/fcc/ee/analyses_storage/flavor/Bs2TauTau/flatNtuples/winter2023/analysis_stage1_mumu_noFilter/"

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
    #rdf = rdf.Filter("n_muons > 1")


    ##### To Add Any RDF Treatment Needed #####

    ###########################################

    print(f"{Mode}: Number of events = {rdf.Count().GetValue()}")
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
    fig.savefig(f'Plots/{var}.pdf')

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