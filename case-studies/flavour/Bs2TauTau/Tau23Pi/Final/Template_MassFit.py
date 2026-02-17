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

#Properly evaluated Stage1 (Preeff (hasPV == 1) * Stage1 cuts) not necessary here since we work with the cutted out dataset but given explicitly
#Warning does not include the removal of the phantom events where no diTau are found 
EffStage1 = {}
EffStage1["sig"] = 4544951/10000000
EffStage1["bb"]  = 97114/434383092
EffStage1["cc"]  = 2310/494686495
EffStage1["ss"]  = 709/499842440
EffStage1["ud"]  = 414/497658654 

#Luminostiy Scaling
LumiScale = {}
LumiScale["sig"] = NZ/10000000 #if only a fraction of the signal, the 10M should be changed to the fraction of total number (for instance if 10 files instead of 20, change it to 5M)
LumiScale["bb"]  = NZ/434383092
LumiScale["cc"]  = NZ/494686495
LumiScale["ss"]  = NZ/499842440
LumiScale["ud"]  = NZ/497658654

#Our efficiencies from Stage 1 and Preff, not to be used here because of the implicit cuts in Stage 1 and collinear mass computation
#eff_Plain = {}# =  has PV == 1                  *           #Stage 1
#eff_Plain["bb"]  = 4356115/(4*500000+6*400000)*             0.0002408
#eff_Plain["cc"]  = 5088774/(600000+8*500000+494560)*        5.1e-06
#eff_Plain["ss"]  = 5060852/(492297+600000+7*500000+471230)* 1.4e-06
#eff_Plain["ud"]  = 4975473/(600000+7*500000+470159+408276)* 6e-07
#eff_Plain["sig"] = 889007/(9*100000)*                       0.4673304

#The efficiencies for the collinear mass, without EVT_hemEmin_NTau23PiCandidates == 2
#eff = {}                     
#eff["bb_mDiTau_collinear3D"] = 0.1629636
#eff["cc_mDiTau_collinear3D"] = 0.0740768
#eff["ss_mDiTau_collinear3D"] = 0.0283642
#eff["ud_mDiTau_collinear3D"] = 0.0112108
#eff["sig_mDiTau_collinear3D"]= 0.579622

#The total amount of events per file to scale the them to the lumi
#NSim = {}
#NSim["bb"] = 4*500000+6*400000
#NSim["cc"] = 600000+8*500000+494560
#NSim["ss"] = 492297+600000+7*500000+471230
#NSim["ud"] = 600000+7*500000+470159+408276
#NSim["sig"]= 9*100000

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
        #rdf = rdf.Filter("MVA2 > 0.6")
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

        #Make the toy by Poisson varying it (give again the bin edges to keep the same data structure)
        np.random.seed(12) #for reproducibility
        SimData[var] = (np.random.poisson(htemps),hlist["sig_"+var][1])

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
        for i in np.arange(0,len(Data[var][1])-1,1):
            bin_center.append((Data[var][1][i+1]+Data[var][1][i])/2.0)

        #Actual fitting, the lambda expression is to tell curve_fit that only the param have to be fitted while the shape should not
        popt[var], pcov[var] = curve_fit(lambda bin_center, *param: Yields(bin_center,SigShape[var][0],BkgShape[var][0],param), bin_center, Data[var][0], p0=param, sigma=np.sqrt(Data[var][0]), absolute_sigma=True)

        #Compute Chi2
        Chi2[var] = 0.0
        for i in np.arange(0,len(Data[var][1])-1,1):
            Chi2[var] += (Data[var][0][i]-(SigShape[var][0][i]*popt[var][0]+BkgShape[var][0][i]*popt[var][1]))**2/np.sqrt(Data[var][0][i])**2/(len(Data[var][1])-1-2)


        print(f"\nTotal Signal Events = {int(np.sum(popt[var][0]*SigShape[var][0]))} pm {int(np.sum(np.sqrt(pcov[var][0][0])*SigShape[var][0]))}")
        print(f"Total Background Events = {int(np.sum(popt[var][1]*BkgShape[var][0]))} pm {int(np.sum(np.sqrt(pcov[var][1][1])*BkgShape[var][0]))}\n")

    return popt, pcov, Chi2


def Draw_Hists(Data,SigShape,BkgShape,Fitted_Val,Vars):

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
        ax.set_ylim([1.01e0,1e9])

        #Labels
        ax.set_xlabel(r"$\textrm{BDT Score}$",size="large")
        ax.set_ylabel(r"$\textrm{Events}$",size="large")
        #Collinear mass
        #ax.set_xlabel(r"$\textrm{Collinear Mass [GeV]}$",size="large")
        #ax.set_ylabel(r"$\textrm{Events / ("+f"{Data[var][1][1]-Data[var][1][0]}"+r" GeV)}$",size="large")
        ax.text(0.489,1.015,r"\textrm{FCC-ee Simulation (IDEA Delphes)}",size="large",transform=ax.transAxes)

        #Legend
        patch1 = mpatches.Patch(color='firebrick', label=r"$B_s^0\rightarrow\tau^+\tau^-(\tau\rightarrow 3\pi)$")
        patch2 = mpatches.Patch(color='steelblue', label=r"$Z^0\rightarrow q\overline{q}\textrm{, Bkg}$")
        patch3 = Line2D([0], [0], marker='.', color='black', label=r"$\textrm{data}$")
        patch4 = Line2D([0,0], [0,1], label=r'$\textrm{Total fit}$', color='k')
        ax.legend(handles=[patch4,patch1,patch2,patch3], frameon=True, framealpha=1, fancybox=True, edgecolor='lightgrey', loc="center left", bbox_to_anchor=(0.6,0.85))

        plt.savefig(f"MassFit_NoBug_{var}.pdf")    
            
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
    print(popt,pcov,Chi2)
    print("Drawing...")
    Draw_Hists(Data,Sig,Bkg,popt,Vars)

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

