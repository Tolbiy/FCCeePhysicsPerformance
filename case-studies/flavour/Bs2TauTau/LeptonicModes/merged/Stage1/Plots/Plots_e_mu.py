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
        
    Path_m = "/eos/experiment/fcc/ee/analyses_storage/flavor/Bs2TauTau/flatNtuples/winter2023/analysis_stage1_mumu_noFilter/"
    Path_e = "/eos/experiment/fcc/ee/analyses_storage/flavor/Bs2TauTau/flatNtuples/winter2023/analysis_stage1_ee_noFilter/"

    #Load the rdf properly (avoid segfault from Tree going out of scope)
    filenames = r.std.vector('string')()
    if "sig" in Mode:
        for i in np.arange(0,20,1):
            if "_m" in Mode:
                filenames.push_back(Path_m+Links["sig"]+f"/chunk_{i}.root")
            else:
                filenames.push_back(Path_e+Links["sig"]+f"/chunk_{i}.root")
    else:
        for i in np.arange(0,10,1):
            filenames.push_back(Path_m+Links[Mode]+f"/chunk_{i}.root")
    rdf = r.RDataFrame("events",filenames)
    print(f"{Mode}: Number of events = {rdf.Count().GetValue()}")
    #rdf = rdf.Filter("n_lepton > 1 && EVT_ThrustEmin_E < 38 && recoEmiss_e > 10 && EVT_ThrustEmin_Eneutral < 10 && has_dilepton > 0 && has_dilepton_SameSide > 0 && has_dilepton_SigHemi > 0 && has_dilepton_OppositeCharges > 0 && has_dilepton_Vertex == 0")
    #print(f"{Mode}: Number of selected events = {rdf.Count().GetValue()}")
    ### /!\ Removed cuts once the filtered samples are produced

    #Compute the dilepton properties
    #for prop in ["energy","px","py","pz"]:    
    #    rdf = rdf.Define(f"dileptonplus_{prop}",f"lepton_{prop}.at(dilepton_plus_ind)")
    #    rdf = rdf.Define(f"dileptonminus_{prop}",f"lepton_{prop}.at(dilepton_minus_ind)")
    #rdf = rdf.Define("dileptonplus_p4","TLorentzVector(dileptonplus_px,dileptonplus_py,dileptonplus_pz,dileptonplus_energy)")
    #rdf = rdf.Define("dileptonminus_p4","TLorentzVector(dileptonminus_px,dileptonminus_py,dileptonminus_pz,dileptonminus_energy)")
    #rdf = rdf.Define("dilepton_Vismass","(dileptonplus_p4+dileptonminus_p4).M()")
    #rdf = rdf.Define("dileptonplus_pt","sqrt(dileptonplus_px*dileptonplus_px + dileptonplus_py*dileptonplus_py)")
    #rdf = rdf.Define("dileptonminus_pt","sqrt(dileptonminus_px*dileptonminus_px + dileptonminus_py*dileptonminus_py)")
    #rdf = rdf.Define("dileptonplus_p","sqrt(dileptonplus_px*dileptonplus_px + dileptonplus_py*dileptonplus_py + dileptonplus_pz*dileptonplus_pz)")
    #rdf = rdf.Define("dileptonminus_p","sqrt(dileptonminus_px*dileptonminus_px + dileptonminus_py*dileptonminus_py + dileptonplus_pz*dileptonplus_pz)")

    #Compute the collinear mass
    #rdf = rdf.Define("recoEmiss_p4",  "TLorentzVector(recoEmiss_px, recoEmiss_py, recoEmiss_pz, recoEmiss_e)")
    #rdf = rdf.Define("recoEmiss_p","recoEmiss_p4.P()")
    #rdf = rdf.Define("ratio_Emiss","(recoEmiss_px * dileptonplus_px + recoEmiss_py * dileptonplus_py + recoEmiss_pz * dileptonplus_pz) / (recoEmiss_p4.P() * dileptonplus_p4.P())")
    #rdf = rdf.Define("dilepton_Fullmass","(dileptonplus_p4+dileptonminus_p4+ratio_Emiss*recoEmiss_p4).M()")
    #rdf = rdf.Define("LpxLm_x",    "dileptonplus_py*dileptonminus_pz-dileptonplus_pz*dileptonminus_py")
    #rdf = rdf.Define("LpxLm_y",    "dileptonplus_pz*dileptonminus_px-dileptonplus_px*dileptonminus_pz")
    #rdf = rdf.Define("LpxLm_z",    "dileptonplus_px*dileptonminus_py-dileptonplus_py*dileptonminus_px")
    #rdf = rdf.Define("LpxEm_x",    "dileptonplus_py*recoEmiss_pz-dileptonplus_pz*recoEmiss_py")
    #rdf = rdf.Define("LpxEm_y",    "dileptonplus_pz*recoEmiss_px-dileptonplus_px*recoEmiss_pz")
    #rdf = rdf.Define("LpxEm_z",    "dileptonplus_px*recoEmiss_py-dileptonplus_py*recoEmiss_px")
    #rdf = rdf.Define("EmxLm_x",    "recoEmiss_py*dileptonminus_pz-recoEmiss_pz*dileptonminus_py")
    #rdf = rdf.Define("EmxLm_y",    "recoEmiss_pz*dileptonminus_px-recoEmiss_px*dileptonminus_pz")
    #rdf = rdf.Define("EmxLm_z",    "recoEmiss_px*dileptonminus_py-recoEmiss_py*dileptonminus_px")
    #rdf = rdf.Define("LpxLm_M2",   "LpxLm_x * LpxLm_x + LpxLm_y * LpxLm_y + LpxLm_z * LpxLm_z")
    #rdf = rdf.Define("denom",      "(LpxLm_x + EmxLm_x)*(LpxLm_x + LpxEm_x) + (LpxLm_y + EmxLm_y)*(LpxLm_y + LpxEm_y) + (LpxLm_z + EmxLm_z)*(LpxLm_z + LpxEm_z)")
    #rdf = rdf.Define("dilepton_Collmass",    "dilepton_Vismass/sqrt(LpxLm_M2/denom)")

    #rdf = rdf.Define("EAsymm","(EVT_ThrustEmax_E-EVT_ThrustEmin_E)/(EVT_ThrustEmin_E+EVT_ThrustEmax_E)")
    #Find garbage particles around candidates (since bkg is mostly in the form B->D->K, should be dealt with by reuiring no vertex)

    #Track related variables
    #rdf = rdf.Define("dileptonplus_DOCA","Compute_DOCA(dileptonplus_px,dileptonminus_py,dileptonplus_pz,dileptonminus_px,dileptonminus_py,dileptonminus_pz)")

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
    ax.stairs(hlist["sig_e"]/Norm["sig_e"],edges,ec='firebrick',ls='-',lw=2,label=r"$\textrm{Signal (electrons)}$")
    ax.stairs(hlist["sig_m"]/Norm["sig_m"],edges,ec='darkorange',ls='-',lw=2,label=r"$\textrm{Signal (muons)}$")

    #Set axis range
    ax.set_xlim([edges[0],edges[-1]])

    #Labels
    ax.set_xlabel(r"$\textrm{Neutral signal hemisphere energy [GeV]}$",size="large")
    ax.set_ylabel(r"$\textrm{Normalised events}$",size="large")

    #Legend (Count the number of bkg modes drawn to write it down properly)
    ax.legend()
        
    #Save
    fig.savefig(f'{var}_ForPres.pdf')

#=======================================================================================================

def Do_Plots(var,bins):

    for v in var:
        hmode = {}
        hmode["sig_m"], edges = Load_AsNumpy("sig_m",v,bins,True)
        hmode["sig_e"] = Load_AsNumpy("sig_e",v,bins,False)
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