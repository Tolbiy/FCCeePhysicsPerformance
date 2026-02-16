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
    fig.savefig(f'Plots/{var}.pdf')

#=======================================================================================================

def Do_Plots(var,bins):

    for v in var:
        hmode = {}
        hmode["sig"], edges = Load_AsNumpy("sig",v,bins,True)
        for mode in ["bb","cc","ss","ud"]:
            hmode[mode] = Load_AsNumpy("bb",v,bins,False)
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