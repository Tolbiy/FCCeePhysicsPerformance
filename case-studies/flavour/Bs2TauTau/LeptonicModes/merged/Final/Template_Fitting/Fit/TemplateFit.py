import numpy as np
import ROOT as r
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import argparse
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amssymb}'

#- Physics number --------------------------------------------------------------

#Some numbers to evaluate S and B
NZ = 6e12 #Number of Z (To be refined depending on FCC scenarios)
Bs_had = 0.1 #Ratio of b->Bs hadronisation (given by Xunwu, to be refined)
Bs2TauTau_BR = 7.73e-7 #SM estimate for Bs->TauTau (from the LHCb paper that searched for it)

#Tau->lnunu BR
Tau2l_PDG = 0.1737**2 + 0.1785**2 + 2*0.1737*0.1785 #from PDG values
Tau2l_Sim = 1250154/10000000 #from sim, keep in mind the num here includes the hasPV

#- BR of Z->qq ----------------------------------------------------------------
BR_qq = {}
BR_qq["bb"] = 0.1512
BR_qq["cc"] = 0.1203
BR_qq["ss"] = 0.1560
BR_qq["ud"] = 0.6991 - BR_qq["bb"] - BR_qq["cc"] - BR_qq["ss"]

#- Efficiencies up to BDT training --------------------------------------------

#bkg eff preBDT (hasPV, Stage1 cuts)
PreBDTEff = {}
PreBDTEff["bb"] = 2105721/438738637
PreBDTEff["cc"] = 6802/499786495
#sig eff preBDT (hasPV, Stage1 cuts from MC decay selections)
PreBDTEff["sig"] = 522843/1250154 #Warning the sig denom is with hasPV, hence the number should be close but not exact

#- Lumi scaling ----------------------------------------------------------------

LumiScale = {}
LumiScale["sig"] = NZ/1250154 #if only a fraction of the signal, the 10M should be changed to the fraction of total number (for instance if 10 files instead of 20, change it to 5M)
LumiScale["bb"]  = NZ/438738637
LumiScale["cc"]  = NZ/499786495

#==================================================================================================================================

#Load the TTree for each mode, get the rdf, get the numpy histogram
def Load_Files_Hists(Vars,Nbins,BDTName):

    #dict to be returned with the numpy histo
    hnp = {}

    #Search for the files
    pathtofile = "../Data/"

    #Split per mode
    for mode in ["sig","bb","cc"]:
        rdf = r.RDataFrame("events",pathtofile+mode+"/"+BDTName+".root")
        #rdf = rdf.Filter("MVA2 > 0.6")
        for var in Vars:
            hnp[mode+"_"+var] = np.histogram(rdf.AsNumpy([var])[var],bins=Nbins[0],range=(Nbins[1],Nbins[2]))

    return hnp

#------------------------------------------------------------------------------------------------

#Create the Sig and Bkg shapes to be fitted to the data
def Make_Shapes(hlist,Vars):

    SigShape = {}
    BkgShape = {}

    for var in Vars:
        
        #Scale the different such that there relative shape agrees
        htemp = {}
        for mode in ["bb","cc"]:
            htemp[mode] = hlist[mode+"_"+var][0]*LumiScale[mode]*BR_qq[mode]
        htemp["sig"] = hlist["sig_"+var][0]*LumiScale["sig"]*BR_qq["bb"]*Bs_had*Bs2TauTau_BR*Tau2l_Sim #Scaling the sig is technically not necessary
        
        #Temporary solution to fit a realistic dataset (with the expected total number of events) while not all samples have been produced
        #Let's amplify the bb bkg (main source of the shape in any case) such that the background accounts for the total amount of Z
        #The ratio done here are not perfect: 1% stat evaluation + not exactly the same amount of ss and ud as bb
        htemp["ss"] = BR_qq["ss"]/BR_qq["bb"]*3.754e-6/4.906895e-3*htemp["bb"]
        htemp["ud"] = BR_qq["ss"]/BR_qq["bb"]*1.3868e-5/4.906895e-3*htemp["bb"]
        
        #Sum the backgrounds together to get a single shape
        BkgShape[var] = (htemp["bb"]+htemp["cc"]+htemp["ss"]+htemp["ud"],hlist["sig_"+var][1])
        
        #The signal
        SigShape[var] = (htemp["sig"],hlist["sig_"+var][1])

    return SigShape, BkgShape

#------------------------------------------------------------------------------------------------

#Scale the initial distr and Poisson vary it to simulate data to be fitted
def Make_Data(hlist,Vars,bkg_dummyEff,sig_dummyEff):

    #Store the various sim data histo
    SimData = {}    

    for var in Vars:

        htemp = {}
        for mode in ["bb","cc"]:
            htemp[mode] = np.rint(hlist[mode+"_"+var][0]*LumiScale[mode]*BR_qq[mode]*bkg_dummyEff) #Round to the int, still a float type hope it works fine with the Poisson toy
        htemp["sig"] = np.rint(hlist["sig_"+var][0]*LumiScale["sig"]*BR_qq["bb"]*Bs_had*Bs2TauTau_BR*Tau2l_Sim*sig_dummyEff)

        #Temporary solution to fit a realistic dataset (with the expected total number of events) while not all samples have been produced
        #Let's amplify the bb bkg (main source of the shape in any case) such that the background accounts for the total amount of Z
        #The ratio done here are not perfect: 1% stat evaluation + not exactly the same amount of ss and ud as bb
        htemp["ss"] = np.rint(BR_qq["ss"]/BR_qq["bb"]*3.754e-6/4.906895e-3*htemp["bb"])
        htemp["ud"] = np.rint(BR_qq["ss"]/BR_qq["bb"]*1.3868e-5/4.906895e-3*htemp["bb"])

        print(f"Total Data Background = {np.sum(htemp['bb']+htemp['cc']+htemp['ss']+htemp['ud'])}")
        print(f"Total Data Signal = {np.sum(htemp['sig'])}")

        #Sum the histo together to get the simulated data distribution
        htemps = htemp["sig"]+htemp["bb"]+htemp["cc"]+htemp["ss"]+htemp["ud"]

        #Make the toy by Poisson varying it (give again the bin edges to keep the same data structure)
        np.random.seed(12) #for reproducibility
        SimData[var] = (np.random.poisson(htemps),hlist["sig_"+var][1])

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

        print(f"--- {var} FIT RESULTS ---")
        print(f"\nTotal Fitted Signal Events = {int(np.sum(popt[var][0]*SigShape[var][0]))} pm {int(np.sum(np.sqrt(pcov[var][0][0])*SigShape[var][0]))}")
        print(f"Total Fitted Background Events = {int(np.sum(popt[var][1]*BkgShape[var][0]))} pm {int(np.sum(np.sqrt(pcov[var][1][1])*BkgShape[var][0]))}\n")
        
        print(f"Signal Scale = {popt[var][0]}")
        print(f"Background Scale = {popt[var][1]}")
        print(f"Chi2/ndf = {Chi2[var]}")
        print("---------------------------")

    return popt, pcov, Chi2

#------------------------------------------------------------------------------------------------

def Draw_Hists(Data,SigShape,BkgShape,Fitted_Val,Vars,BDTName):

    for var in Vars:

        fig, ax = plt.subplots()
        
        #Get the bin centers
        bin_center = []
        for i in np.arange(0,len(Data[var][1])-1,1):
            bin_center.append((Data[var][1][i+1]+Data[var][1][i])/2.0)

        #Actually draw
        #ax.set_yscale('log', base=10) not that useful since to see something for the signal the shape must not span too much order of magnitude
        ax.grid(color='grey', linestyle='--', linewidth=0.5, alpha=0.5)
        ax.errorbar(bin_center,Data[var][0],yerr=np.sqrt(Data[var][0]),marker='.',linestyle='',color='black')
        ax.stairs(Fitted_Val[var][1]*BkgShape[var][0] + Fitted_Val[var][0]*SigShape[var][0],BkgShape[var][1],fill=True,color='steelblue',ec='black',ls='-',lw=1) #Kind of shitty way to stack but works for now
        ax.stairs(Fitted_Val[var][0]*SigShape[var][0],SigShape[var][1],fill=True,color='firebrick')

        #Compute and display pulls
        Pulls = (Fitted_Val[var][1]*BkgShape[var][0]+Fitted_Val[var][0]*SigShape[var][0] - Data[var][0])/np.sqrt(Data[var][0])
        divider = make_axes_locatable(ax)
        ax_pulls = divider.append_axes("bottom", 0.6, pad=0.05, sharex=ax)
        ax.xaxis.set_tick_params(labelbottom=False,bottom=False)
        ax_pulls.plot([0.0,1.0],[0.0,0.0],color="black",lw=1,ls="--")
        ax_pulls.errorbar(bin_center,Pulls,yerr=1.0,marker='.',linestyle='',color="black")
        ax_pulls.set_ylim([-5,5])
        ax_pulls.set_xlim([Data[var][1][0],Data[var][1][-1]])
        ax_pulls.set_xlabel(r"$\textrm{BDT Score}$",size="large")
        ax_pulls.set_ylabel(r"$\textrm{Pulls}$",size="large")

        #Set axis range
        ax.set_xlim([Data[var][1][0],Data[var][1][-1]])
        ax.set_yscale("log")
        ax.set_ylim([1.01e0,1e8])

        #Labels
        ax.set_xlabel(r"$\textrm{BDT Score}$",size="large")
        ax.set_ylabel(r"$\textrm{Events}$",size="large")
        #Collinear mass
        #ax.set_xlabel(r"$\textrm{Collinear Mass [GeV]}$",size="large")
        #ax.set_ylabel(r"$\textrm{Events / ("+f"{Data[var][1][1]-Data[var][1][0]}"+r" GeV)}$",size="large")
        ax.text(0.489,1.015,r"\textrm{FCC-ee Simulation (IDEA Delphes)}",size="large",transform=ax.transAxes)

        #Legend
        patch1 = mpatches.Patch(color='firebrick', label=r"$B_s^0\rightarrow\tau^+\tau^-(\tau\rightarrow\ell\nu_{\tau}\nu_{\ell})$")
        patch2 = mpatches.Patch(color='steelblue', label=r"$Z^0\rightarrow q\overline{q}\textrm{, Bkg}$")
        patch3 = Line2D([0], [0], marker='.', color='black', label=r"$\textrm{data}$")
        patch4 = Line2D([0,0], [0,1], label=r'$\textrm{Total fit}$', color='k')
        ax.legend(handles=[patch4,patch1,patch2,patch3], frameon=True, framealpha=1, fancybox=True, edgecolor='lightgrey', loc="center left", bbox_to_anchor=(0.02,0.22))

        plt.savefig(f"{BDTName}_{var}.pdf")    
            
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Collect info, scale and make toys, fits, draw according to a dumy eff to amplify the signal
def Do_MassFit(Vars,BDTName,Bins,BkgEff,SigEff):

    print("Loading...")
    h = Load_Files_Hists(Vars,Bins,BDTName)
    print("Getting the Sig and Bkg shapes...")
    Sig, Bkg = Make_Shapes(h,Vars)
    print("Getting the Simulated data...\n")
    Data = Make_Data(h,Vars,BkgEff,SigEff)

    print("\nFitting Sig+Bkg shape to Data...\n")
    popt, pcov, Chi2 = Fitter(Data,Sig,Bkg,Vars)
    print("\nDrawing...")
    Draw_Hists(Data,Sig,Bkg,popt,Vars,BDTName)

#====================================================================================================================================================


#MassRes = 0.2 # in GeV Arbitrary value, need to evaluate properly with pion mass maybe
#MLow = 0
#MHigh = 6

parser = argparse.ArgumentParser()
parser.add_argument("BkgEff", type=float)
parser.add_argument("SigEff", type=float) #To amplify the sig for drawing purposes only (between 2-5)
parser.add_argument("MassRes",type=float)
parser.add_argument("MLow",   type=float)
parser.add_argument("MHigh",  type=float)
parser.add_argument("BDTName",type=str)
parser.add_argument("List",   nargs='+')
args=parser.parse_args()

Nbin = int((args.MHigh-args.MLow)/args.MassRes)

Do_MassFit(args.List,args.BDTName,(Nbin,args.MLow,args.MHigh),args.BkgEff,args.SigEff)

