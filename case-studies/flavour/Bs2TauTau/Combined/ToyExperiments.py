import numpy as np
import ROOT as r
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import argparse
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tqdm import tqdm

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amssymb}'
plt.rcParams["figure.figsize"] = (9,9)

#- Physics number --------------------------------------------------------------

#Some numbers to evaluate S and B
NZ = 6e12 #Number of Z (To be refined depending on FCC scenarios)
Bs_had = 0.1 #Ratio of b->Bs hadronisation (given by Xunwu, to be refined)
Bs2TauTau_BR = 7.73e-7 #SM estimate for Bs->TauTau (from the LHCb paper that searched for it)

#Tau->3pinu
Tau23pi = 0.0931**2

#Tau->lnunu BR
Tau2l_PDG = 0.1737**2 + 0.1785**2 + 2*0.1737*0.1785 #from PDG values
Tau2l_Sim = 1250154/10000000 #from sim, keep in mind the num here includes the hasPV

#- BR of Z->qq ----------------------------------------------------------------
BR_qq = {}
BR_qq["bb"] = 0.1512
BR_qq["cc"] = 0.1203
BR_qq["ss"] = 0.1560
BR_qq["ud"] = 0.6991 - BR_qq["bb"] - BR_qq["cc"] - BR_qq["ss"]

#- Efficiencies of tau2l and tau23pi up to BDT training --------------------------------------------

#bkg eff preBDT (hasPV, Stage1 cuts)
PreBDTEff = {}


PreBDTEff["bb_l"] = 2105721/438738637
PreBDTEff["cc_l"] = 6802/499786495
PreBDTEff["ss_l"] = 1416/489770989
PreBDTEff["ud_l"] = 1338/492658654

PreBDTEff["bb_3pi"] = 97114/434383092
PreBDTEff["cc_3pi"] = 2310/494686495
PreBDTEff["ss_3pi"] = 709/499842440
PreBDTEff["ud_3pi"] = 414/499786495

#sig eff preBDT (hasPV, Stage1 cuts from MC decay selections)
PreBDTEff["sig_l"] = 522843/1250154 #Warning the sig denom is with hasPV, hence the number should be close but not exact
PreBDTEff["sig_3pi"] = 4544951/10000000

#- Lumi scaling -------------------------------------------------------------------

LumiScale = {}

#split because of corruption
LumiScale["bb_l"]  = NZ/438738637
LumiScale["cc_l"]  = NZ/499786495
LumiScale["ss_l"]  = NZ/489770989
LumiScale["ud_l"]  = NZ/492658654

LumiScale["bb_3pi"]  = NZ/434383092
LumiScale["cc_3pi"]  = NZ/494686495
LumiScale["ss_3pi"]  = NZ/499842440
LumiScale["ud_3pi"]  = NZ/497658654

LumiScale["sig_l"] = NZ/1250154
LumiScale["sig_3pi"] = NZ/10000000

#==================================================================================================================================

#Load the TTree for each mode, get the rdf, get the numpy histogram
def Load_Files_Hists(var,DataPath,Nbins):

    rdf = r.RDataFrame("events",DataPath)
    hm = np.histogram(rdf.AsNumpy([var])[var],bins=Nbins[0],range=(Nbins[1],Nbins[2]))

    return hm

#------------------------------------------------------------------------------------------------

#Create the Sig and Bkg shapes to be fitted to the data
def Make_Shapes(hlist,Vars):

    SigShape = {}
    BkgShape = {}

    for var in Vars:
        
        #Scale the different such that there relative shape agrees
        htemp = {}

        for mode in ["bb","cc","ss","ud"]:
            htemp[mode+"_l"] = hlist[mode+"_l_"+var][0]*LumiScale[mode+"_l"]*BR_qq[mode]
            htemp[mode+"_3pi"] = hlist[mode+"_3pi_"+var][0]*LumiScale[mode+"_3pi"]*BR_qq[mode]
        htemp["sig_l"] = hlist["sig_l_"+var][0]*LumiScale["sig_l"]*BR_qq["bb"]*Bs_had*Bs2TauTau_BR*Tau2l_Sim
        htemp["sig_3pi"] = hlist["sig_3pi_"+var][0]*LumiScale["sig_3pi"]*BR_qq["bb"]*Bs_had*Bs2TauTau_BR*Tau23pi
        
        #Fuse the bin edges (I just add them one behind the other by removing the first element of the second list of edges and shifting them to the end of the first list of edges)
        #Not sure how this will affect the fit but shouldn't since the only thing we care about is the bin height
        Edge2 = hlist["sig_3pi_"+var][1][1:]
        Shift = hlist["sig_l_"+var][1][-1]-hlist["sig_l_"+var][1][0]
        for i in range(len(Edge2)):
            Edge2[i] += Shift
        #print(hlist["sig_l_"+var][1])
        #print(Edge2)
        FusedEdges = np.concatenate((hlist["sig_l_"+var][1],Edge2))

        #Sum the backgrounds together to get a single shape and fuse the different decay histo into a single one for combined fit
        FusedBkg = np.concatenate((htemp["bb_l"]+htemp["cc_l"]+htemp["ss_l"]+htemp["ud_l"],htemp["bb_3pi"]+htemp["cc_3pi"]+htemp["ss_3pi"]+htemp["ud_3pi"]))
        BkgShape[var] = (FusedBkg,FusedEdges)
        
        #Same in signal
        FusedSig = np.concatenate((htemp["sig_l"],htemp["sig_3pi"]))
        SigShape[var] = (FusedSig,FusedEdges)

    return SigShape, BkgShape

#------------------------------------------------------------------------------------------------

#Scale the initial distr and Poisson vary it to simulate data to be fitted
def Make_Data(hlist,Vars,bkg_dummyEff,sig_dummyEff,seed):

    #Store the various sim data histo
    SimData = {}    

    for var in Vars:

        htemp = {}
        for mode in ["bb","cc","ss","ud"]:
            htemp[mode+"_l"] = np.rint(hlist[mode+"_l_"+var][0]*LumiScale[mode+"_l"]*BR_qq[mode]*bkg_dummyEff) #Round to the int, still a float type hope it works fine with the Poisson toy
            htemp[mode+"_3pi"] = np.rint(hlist[mode+"_3pi_"+var][0]*LumiScale[mode+"_3pi"]*BR_qq[mode]*bkg_dummyEff)
        htemp["sig_l"] = np.rint(hlist["sig_l_"+var][0]*LumiScale["sig_l"]*BR_qq["bb"]*Bs_had*Bs2TauTau_BR*Tau2l_Sim*sig_dummyEff)
        htemp["sig_3pi"] = np.rint(hlist["sig_3pi_"+var][0]*LumiScale["sig_3pi"]*BR_qq["bb"]*Bs_had*Bs2TauTau_BR*Tau23pi*sig_dummyEff)

        if seed == 10: #Print only at the beginning for info
            print(f"Total Tau2L Data Background = {np.sum(htemp['bb_l']+htemp['cc_l']+htemp['ss_l']+htemp['ud_l'])}")
            print(f"Total Tau23Pi Data Background = {np.sum(htemp['bb_3pi']+htemp['cc_3pi']+htemp['ss_3pi']+htemp['ud_3pi'])}")
            print(f"Total Tau2L Data Signal = {np.sum(htemp['sig_l'])}")
            print(f"Total Tau23Pi Data Signal = {np.sum(htemp['sig_3pi'])}")

        #Fuse the bin edges (I just add them one behind the other by removing the first element of the second list of edges and shifting them to the end of the first list of edges)
        #Not sure how this will affect the fit but shouldn't since the only thing we care about is the bin height
        Edge2 = hlist["sig_3pi_"+var][1][1:]
        Shift = hlist["sig_l_"+var][1][-1]-hlist["sig_l_"+var][1][0]
        for i in range(len(Edge2)):
            Edge2[i] += Shift
        FusedEdges = np.concatenate((hlist["sig_l_"+var][1],Edge2))

        #Sum the backgrounds together to get a single shape and fuse the different decay histo into a single one for combined fit
        FusedBkg = np.concatenate((htemp["bb_l"]+htemp["cc_l"]+htemp["ss_l"]+htemp["ud_l"],htemp["bb_3pi"]+htemp["cc_3pi"]+htemp["ss_3pi"]+htemp["ud_3pi"]))

        #Same in signal
        FusedSig = np.concatenate((htemp["sig_l"],htemp["sig_3pi"]))

        #Make the toy by Poisson varying it (give again the bin edges to keep the same data structure)
        np.random.seed(seed) #for reproducibility
        SimData[var] = (np.random.poisson(FusedSig+FusedBkg),FusedEdges)

    return SimData

#------------------------------------------------------------------------------------------------

#The fit performed, Some troubles defining it such that we do have Bkg and Sig yields (both should be applied separatly with a constaint such that Ysig + Ybkg = Total events)
def Yields(bin_center,SigShape,BkgShape,yields):
    
    TotShape = []
    for i in np.arange(0,len(bin_center),1):
        TotShape.append(yields[0]*SigShape[i]+yields[1]*BkgShape[i])
    return TotShape


#Wrapper to start fitting
def Fitter(Data,SigShape,BkgShape,Vars):

    popt = {}
    pcov = {}
    Chi2 = {}

    param = [1.0,1.0] #Nsig,Nbkg

    for var in Vars:

        #Get the bin centers
        bin_center = []
        for i in np.arange(0,len(Data[var][1])-1,1):
            bin_center.append((Data[var][1][i+1]+Data[var][1][i])/2.0)

        #Actual fitting, the lambda expression is to tell curve_fit that only the param have to be fitted while the shape should not
        popt[var], pcov[var] = curve_fit(lambda bin_center, *param: Yields(bin_center,SigShape[var][0],BkgShape[var][0],param), bin_center, Data[var][0], p0=param, sigma=np.sqrt(Data[var][0]), absolute_sigma=True)

        #Compute Chi2
        Chi2[var] = 0.0
        for i in np.arange(0,len(Data[var][1])-1,1):
            Chi2[var] += (Data[var][0][i]-(SigShape[var][0][i]*popt[var][0]+BkgShape[var][0][i]*popt[var][1]))**2/np.sqrt(Data[var][0][i])**2/(len(Data[var][1])-1-2)

    return popt, pcov, Chi2

#------------------------------------------------------------------------------------------------

def Draw_Toys(ValDist,Vars,BDTNames):

    for var in Vars:
        
        fig, axs = plt.subplots(2,2)
        for ax in axs.flat:
            ax.grid(color='grey', linestyle='--', linewidth=0.5, alpha=0.5)

        axs[0,0].hist(ValDist[var]["Best_Sig"],bins=20,color="firebrick")
        axs[0,0].set_xlabel(r"$<N_{\rm sig}>$",fontsize="x-large")
        axs[0,0].set_ylabel(r"\textrm{Occurences}",fontsize="xx-large")
        axs[0,1].hist(ValDist[var]["Sigma_Sig"],bins=20,color="firebrick")
        axs[0,1].set_xlabel(r"$\sigma_{N_{\rm sig}}$",fontsize="x-large")
        #axs[0,1].set_ylabel(r"\textrm{Occurences}",fontsize="large")
        axs[1,0].hist(ValDist[var]["Best_Bkg"],bins=20,color="steelblue")
        axs[1,0].set_xlabel(r"$<N_{\rm bkg}>$",fontsize="x-large")
        axs[1,0].set_ylabel(r"\textrm{Occurences}",fontsize="xx-large")
        axs[1,1].hist(ValDist[var]["Sigma_Bkg"],bins=20,color="steelblue")
        axs[1,1].set_xlabel(r"$\sigma_{N_{\rm bkg}}$",fontsize="x-large")
        #axs[1,1].set_ylabel(r"\textrm{Occurences}",fontsize="large")

        axs[0,0].text(0.5,1.05,r"\textrm{Best Fit Value}",size="xx-large",transform=axs[0,0].transAxes,ha="center",va="center")
        axs[0,1].text(0.5,1.05,r"\textrm{Fit Uncertainty}",size="xx-large",transform=axs[0,1].transAxes,ha="center",va="center")
        #axs[0,0].text(-0.2,0.5,r"\textrm{Signal Distributions}",size="xx-large",transform=axs[0,0].transAxes,ha="center",va="center",rotation="vertical")
        #axs[1,0].text(-0.2,0.5,r"\textrm{Background Distributions}",size="xx-large",transform=axs[1,0].transAxes,ha="center",va="center",rotation="vertical")

        axs[0,0].text(0.02,0.98,r"$\mu="+f"{round(np.array(ValDist[var]['Best_Sig']).mean(),2)}"+r"$"+"\n"+r"$\sigma="+f"{round(np.array(ValDist[var]['Best_Sig']).std(),2)}"+r"$",size="large",transform=axs[0,0].transAxes,ha="left",va="top")
        axs[0,1].text(0.02,0.98,r"$\mu="+f"{round(np.array(ValDist[var]['Sigma_Sig']).mean(),2)}"+r"$"+"\n"+r"$\sigma="+f"{round(np.array(ValDist[var]['Sigma_Sig']).std(),2)}"+r"$",size="large",transform=axs[0,1].transAxes,ha="left",va="top")
        axs[1,0].text(0.02,0.98,r"$\mu="+f"{round(np.array(ValDist[var]['Best_Bkg']).mean(),2)}"+r"$"+"\n"+r"$\sigma="+f"{round(np.array(ValDist[var]['Best_Bkg']).std(),2)}"+r"$",size="large",transform=axs[1,0].transAxes,ha="left",va="top")
        axs[1,1].text(0.02,0.98,r"$\mu="+f"{round(np.array(ValDist[var]['Sigma_Bkg']).mean(),2)}"+r"$"+"\n"+r"$\sigma="+f"{round(np.array(ValDist[var]['Sigma_Bkg']).std(),2)}"+r"$",size="large",transform=axs[1,1].transAxes,ha="left",va="top")


        axs[0,0].text(1.1,1.2,r"\textrm{"+f"{var}"+r" Toys Results (}$N_{\rm toys}="+f"{len(ValDist[var]['Best_Sig'])}"+r"$\textrm{)}",size="xx-large",transform=axs[0,0].transAxes,ha="center",va="center")

        fig.savefig(f"{BDTNames[0]}-{BDTNames[1]}_{var}.pdf")    
            
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Collect info, scale and make toys, fits, draw according to a dumy eff to amplify the signal
def Do_Toys(Vars,BDTNames,Bins,BkgEff,SigEff,NToy):

    tau23pi_datapath = "/afs/cern.ch/work/t/tomonnar/public/Bs2TauTau/Final/Datasets/"
    tau2l_datapath = "/afs/cern.ch/work/t/tomonnar/public/FCCeePhysicsPerformance/case-studies/flavour/Bs2TauTau/LeptonicModes/merged/Final/Template_Fitting/Data/"

    h = {}

    print("Loading Tau23pi...")
    for var in Vars:
        h["sig_3pi_"+var] = Load_Files_Hists(var,tau23pi_datapath+BDTNames[0]+"/sig.root",Bins)
        h["bb_3pi_"+var] = Load_Files_Hists(var,tau23pi_datapath+BDTNames[0]+"/bb.root",Bins)
        h["cc_3pi_"+var] = Load_Files_Hists(var,tau23pi_datapath+BDTNames[0]+"/cc.root",Bins)
        h["ss_3pi_"+var] = Load_Files_Hists(var,tau23pi_datapath+BDTNames[0]+"/ss.root",Bins)
        h["ud_3pi_"+var] = Load_Files_Hists(var,tau23pi_datapath+BDTNames[0]+"/ud.root",Bins)
    
    print("Loading Tau2L...")
    for var in Vars:
        h["sig_l_"+var] = Load_Files_Hists(var,tau2l_datapath+"sig/"+BDTNames[1]+".root",Bins)
        h["bb_l_"+var] = Load_Files_Hists(var,tau2l_datapath+"bb/"+BDTNames[1]+".root",Bins)
        h["cc_l_"+var] = Load_Files_Hists(var,tau2l_datapath+"cc/"+BDTNames[1]+".root",Bins)
        h["ss_l_"+var] = Load_Files_Hists(var,tau2l_datapath+"ss/"+BDTNames[1]+".root",Bins)
        h["ud_l_"+var] = Load_Files_Hists(var,tau2l_datapath+"ud/"+BDTNames[1]+".root",Bins)

    print("Getting the Sig and Bkg shapes...")
    Sig, Bkg = Make_Shapes(h,Vars)
    
    Values = {}
    for var in Vars:
        Values[var] = {"Best_Sig":[],"Sigma_Sig":[],"Best_Bkg":[],"Sigma_Bkg":[]}
    print("Start producing and fitting toys...\n") 
    for i in tqdm(np.arange(0,NToy,1)):
        Data = Make_Data(h,Vars,BkgEff,SigEff,10*i)
        popt, pcov, Chi2 = Fitter(Data,Sig,Bkg,Vars)
        
        for var in Vars:
            Values[var]["Best_Sig"].append(int(np.sum(popt[var][0]*Sig[var][0])))
            Values[var]["Best_Bkg"].append(int(np.sum(popt[var][1]*Bkg[var][0])))

            Values[var]["Sigma_Sig"].append(int(np.sum(np.sqrt(pcov[var][0][0])*Sig[var][0])))
            Values[var]["Sigma_Bkg"].append(int(np.sum(np.sqrt(pcov[var][1][1])*Bkg[var][0])))


    print("\nDrawing...")
    Draw_Toys(Values,Vars,BDTNames)

#====================================================================================================================================================


#MassRes = 0.2 # in GeV Arbitrary value, need to evaluate properly with pion mass maybe
#MLow = 0
#MHigh = 6

parser = argparse.ArgumentParser()
parser.add_argument("NToys",  type=int)
parser.add_argument("BkgEff", type=float)
parser.add_argument("SigEff", type=float) #To amplify the sig for drawing purposes only (between 2-5)
parser.add_argument("MassRes",type=float)
parser.add_argument("MLow",   type=float)
parser.add_argument("MHigh",  type=float)
parser.add_argument("BDTName_3pi",type=str)
parser.add_argument("BDTName_l",type=str)
parser.add_argument("List",   nargs='+')
args=parser.parse_args()

Nbin = int((args.MHigh-args.MLow)/args.MassRes)

Do_Toys(args.List,(args.BDTName_3pi,args.BDTName_l),(Nbin,args.MLow,args.MHigh),args.BkgEff,args.SigEff,args.NToys)

