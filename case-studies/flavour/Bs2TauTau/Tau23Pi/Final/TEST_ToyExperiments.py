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

#-----

#Some numbers to evaluate S and B
NZ = 6e12 #Number of Z (To be refined depending on FCC scenarios)
Bs_had = 0.1 #Ratio of b->Bs hadronisation (given by Xunwu, to be refined)
Bs2TauTau_BR = 7.73e-7 #Highest SM estimate for Bs->TauTau (from the LHCb paper that searched for it)
Tau23Pi = 0.0931 #Tau->3pi estimate (PDG) 

#BR of Z->qq
BR_qq = {}
BR_qq["bb"] = 0.1512
BR_qq["cc"] = 0.1203
BR_qq["ss"] = 0.1560
BR_qq["ud"] = 0.6991 - BR_qq["bb"] - BR_qq["cc"] - BR_qq["ss"]

#Luminostiy Scaling
LumiScale = {}
LumiScale["sig"] = NZ/10000000 #if only a fraction of the signal, the 10M should be changed to the fraction of total number (for instance if 10 files instead of 20, change it to 5M)
LumiScale["bb"]  = NZ/434383092
LumiScale["cc"]  = NZ/494686495
LumiScale["ss"]  = NZ/499842440
LumiScale["ud"]  = NZ/497658654

#-----

#Load the TTree for each mode, get the rdf, get the numpy histogram
def Load_Files_Hists(Vars,Nbins):

    #dict to be returned with the numpy histo
    hnp = {}

    #Search for the files
    pathtofile = "/afs/cern.ch/work/t/tomonnar/public/Bs2TauTau/Final/Datasets/Baseline_NoBug_BDT/"

    #Split per mode
    for mode in ["sig","bb","cc","ss","ud"]:
        rdf = r.RDataFrame("events",pathtofile+mode+".root")
        for var in Vars:
            hnp[mode+"_"+var] = np.histogram(rdf.AsNumpy([var])[var],bins=Nbins[0],range=(Nbins[1],Nbins[2]))

    return hnp


#Create the Sig and Bkg shapes to be fitted to the data
def Make_Shapes(hlist,Vars):

    SigShape = {}
    BkgShape = {}

    for var in Vars:
        
        #Scale the different such that there relative shape agrees
        htemp = {}
        for mode in ["bb","cc","ss","ud"]:
            htemp[mode] = hlist[mode+"_"+var][0]*LumiScale[mode]*BR_qq[mode]
        htemp["sig"] = hlist["sig_"+var][0]*LumiScale["sig"]*BR_qq["bb"]*Bs_had*Bs2TauTau_BR*Tau23Pi**2 #Scaling the sig is technically not necessary
        
        #Sum the backgrounds together to get a single shape
        BkgShape[var] = (htemp["bb"]+htemp["cc"]+htemp["ss"]+htemp["ud"],hlist["sig_"+var][1])
        
        #The signal
        SigShape[var] = (htemp["sig"],hlist["sig_"+var][1])

    return SigShape, BkgShape


#Scale the initial distr and Poisson vary it to simulate data to be fitted
def Make_Data(hlist,Vars,Dummy_Eff,ToDraw):

    #Store the various sim data histo
    SimData = {}    

    for var in Vars:

        #Lumi and eff scale the histo
        #The preeff (HasPV == 1) is embeded in the number of events in the histo -> no need to include it
        #EVT_hemEmin_NTau23PiCandidates == 2 is also embedded because of the way we compute the collinear mass
        #Also assumes for now that the range of the histo includes 100% of the events, might correct for that later (add a mass range eff)

        print(hlist["sig_"+var][0])

        htemp = {}
        for mode in ["bb","cc","ss","ud"]:
            htemp[mode] = np.rint(hlist[mode+"_"+var][0]*LumiScale[mode]*BR_qq[mode]*Dummy_Eff) #Round to the int, still a float type hope it works fine with the Poisson toy
        htemp["sig"] = np.rint(hlist["sig_"+var][0]*LumiScale["sig"]*BR_qq["bb"]*Bs_had*Bs2TauTau_BR*Tau23Pi**2*ToDraw)

        print(htemp["sig"])

        print(f"Total Data Background = {np.sum(htemp['bb']+htemp['cc']+htemp['ss']+htemp['ud'])}")
        print(f"Total Data Signal = {np.sum(htemp['sig'])}")

        #Sum the histo together to get the simulated data distribution
        htemps = htemp["sig"]+htemp["bb"]+htemp["cc"]+htemp["ss"]+htemp["ud"]

        #Make the toys by Poisson varying it (give again the bin edges to keep the same data structure)
        Toys = []
        for i in np.arange(0,1000,1):
            np.random.seed(1+50*i)
            Toys.append((np.random.poisson(htemps),hlist["sig_"+var][1]))
        SimData[var] = Toys

    return SimData


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
        for j in np.arange(0,len(Data[var][0][1])-1,1):
            bin_center.append((Data[var][0][1][j+1]+Data[var][0][1][j])/2.0)

        Toys_Mean = []
        Toys_Sigma = []
        Toys_Chi2 = []

        for i in np.arange(0,len(Data[var]),1):

            #Actual fitting, the lambda expression is to tell curve_fit that only the param have to be fitted while the shape should not
            popt_Toy, pcov_Toy = curve_fit(lambda bin_center, *param: Yields(bin_center,SigShape[var][0],BkgShape[var][0],param), bin_center, Data[var][i][0], p0=param, sigma=np.sqrt(Data[var][i][0]), absolute_sigma=True)

            #Compute Chi2
            Chi2_Toy = 0.0
            for k in np.arange(0,len(Data[var][i][1])-1,1):
                Chi2_Toy += (Data[var][i][0][k]-(SigShape[var][0][k]*popt_Toy[0]+BkgShape[var][0][k]*popt_Toy[1]))**2/np.sqrt(Data[var][i][0][k])**2/(len(Data[var][i][1])-1-2)

            Toys_Mean.append((np.sum(popt_Toy[0]*SigShape[var][0]),np.sum(popt_Toy[1]*BkgShape[var][0])))
            Toys_Sigma.append((np.sum(np.sqrt(pcov_Toy[0][0])*SigShape[var][0]),np.sum(np.sqrt(pcov_Toy[1][1])*BkgShape[var][0])))
            Toys_Chi2.append(Chi2_Toy)

        popt[var] = Toys_Mean
        pcov[var] = Toys_Sigma
        Chi2[var] = Toys_Chi2             

    return popt, pcov, Chi2


def Draw_Toys_Mean(Fitted_Mean,Fitted_Chi2,Vars):

    for var in Vars:

        NSigs = np.array([Fitted_Mean[var][i][0] for i in np.arange(0,1000,1)])
        NBkgs = np.array([Fitted_Mean[var][i][1] for i in np.arange(0,1000,1)])
        
        fig_sig, ax_sig = plt.subplots()
        ax_sig.grid(color='grey', linestyle='--', linewidth=0.5, alpha=0.5)
        
        ax_sig.hist(NSigs, 50)

        #Labels
        ax_sig.set_xlabel(r"$<N_{sig}>$",size="large")
        ax_sig.set_ylabel(r"$\textrm{Number of Toys}$",size="large")
        ax_sig.text(0.489,1.015,r"\textrm{FCC-ee Simulation (IDEA Delphes)}",size="large",transform=ax_sig.transAxes)

        fig_sig.savefig(f"Toys_{var}_NSig_2.pdf")

        #-------

        fig_bkg, ax_bkg = plt.subplots()
        ax_bkg.grid(color='grey', linestyle='--', linewidth=0.5, alpha=0.5)
        
        ax_bkg.hist(NBkgs, 50)

        #Labels
        ax_bkg.set_xlabel(r"$<N_{bkg}>$",size="large")
        ax_bkg.set_ylabel(r"$\textrm{Number of Toys}$",size="large")
        ax_bkg.text(0.489,1.015,r"\textrm{FCC-ee Simulation (IDEA Delphes)}",size="large",transform=ax_bkg.transAxes)

        fig_bkg.savefig(f"Toys_{var}_NBkg_2.pdf")      


def Draw_Toys_Sigma(Fitted_Sigma,Vars):

    for var in Vars:

        Sigma_Sigs = np.array([Fitted_Sigma[var][i][0] for i in np.arange(0,1000,1)])
        Sigma_Bkgs = np.array([Fitted_Sigma[var][i][1] for i in np.arange(0,1000,1)])
        
        fig_sig, ax_sig = plt.subplots()
        ax_sig.grid(color='grey', linestyle='--', linewidth=0.5, alpha=0.5)
        
        ax_sig.hist(Sigma_Sigs)

        #Labels
        ax_sig.set_xlabel(r"$\sigma(N_{sig})$",size="large")
        ax_sig.set_ylabel(r"$\textrm{Number of Toys}$",size="large")
        ax_sig.text(0.489,1.015,r"\textrm{FCC-ee Simulation (IDEA Delphes)}",size="large",transform=ax_sig.transAxes)

        fig_sig.savefig(f"Toys_{var}_SigmaSig_2.pdf")

        #-------

        fig_bkg, ax_bkg = plt.subplots()
        ax_bkg.grid(color='grey', linestyle='--', linewidth=0.5, alpha=0.5)
        
        ax_bkg.hist(Sigma_Bkgs)

        #Labels
        ax_bkg.set_xlabel(r"$\sigma(N_{bkg})$",size="large")
        ax_bkg.set_ylabel(r"$\textrm{Number of Toys}$",size="large")
        ax_bkg.text(0.489,1.015,r"\textrm{FCC-ee Simulation (IDEA Delphes)}",size="large",transform=ax_bkg.transAxes)

        fig_bkg.savefig(f"Toys_{var}_SigmaBkg_2.pdf")     
            
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Collect info, scale and make toys, fits, draw according to a dumy eff to amplify the signal
def Do_MassFit(Vars,Bins,Dummy_Eff,Draw_Eff):

    print("Loading...")
    h = Load_Files_Hists(Vars,Bins)
    print("Getting the Sig and Bkg shapes...")
    Sig, Bkg = Make_Shapes(h,Vars)
    print("Getting the Simulated data...")
    Data = Make_Data(h,Vars,Dummy_Eff,Draw_Eff)

    print("Fitting Sig+Bkg shape to Data...")
    popt, pcov, Chi2 = Fitter(Data,Sig,Bkg,Vars)
    print("Drawing...")
    Draw_Toys_Mean(popt,Chi2,Vars)
    Draw_Toys_Sigma(pcov,Vars)

#====================================================================================================================================================


#MassRes = 0.2 # in GeV Arbitrary value, need to evaluate properly with pion mass maybe
#MLow = 0
#MHigh = 6

parser = argparse.ArgumentParser()
parser.add_argument("Eff",    type=float)
parser.add_argument("DrawEff",type=float) #To amplify the sig for drawing purposes only (between 2-5)
parser.add_argument("MassRes",type=float)
parser.add_argument("MLow",   type=float)
parser.add_argument("MHigh",  type=float)
parser.add_argument("List",   nargs='+')
args=parser.parse_args()

Nbin = int((args.MHigh-args.MLow)/args.MassRes)

Do_MassFit(args.List,(Nbin,args.MLow,args.MHigh),args.Eff,args.DrawEff)