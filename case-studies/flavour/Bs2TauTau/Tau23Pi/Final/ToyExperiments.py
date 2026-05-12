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
Tau23Pi = 0.0931 #Tau->3pi estimate (PDG)

#- BR of Z->qq ----------------------------------------------------------------
BR_qq = {}
BR_qq["bb"] = 0.1512
BR_qq["cc"] = 0.1203
BR_qq["ss"] = 0.1560
BR_qq["ud"] = 0.6991 - BR_qq["bb"] - BR_qq["cc"] - BR_qq["ss"]

#- Efficiencies up to BDT training --------------------------------------------

#bkg eff preBDT (hasPV, Stage1 cuts)
PreBDTEff = {}
PreBDTEff["bb"] = 97114/434383092
PreBDTEff["cc"] = 2310/494686495
PreBDTEff["ss"] = 709/499842440
PreBDTEff["ud"] = 414/499786495
PreBDTEff["sig"] = 4544951/10000000

#- Lumi scaling ----------------------------------------------------------------

LumiScale = {}
LumiScale["bb"]  = NZ/434383092
LumiScale["cc"]  = NZ/494686495
LumiScale["ss"]  = NZ/499842440
LumiScale["ud"]  = NZ/497658654
LumiScale["sig"] = NZ/10000000

#==================================================================================================================================

#Load the TTree for each mode, get the rdf, get the numpy histogram
def Load_Files_Hists(Vars,Nbins,BDTName):

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

#------------------------------------------------------------------------------------------------

#Create the Sig and Bkg shapes to be fitted to the data
def Make_Shapes(hlist,Vars):

    SigShape = {}
    BkgShape = {}

    for var in Vars:
        
        #Scale the different such that there relative shape agrees
        htemp = {}
        for mode in ["bb","cc","ss","ud"]:
            htemp[mode] = hlist[mode+"_"+var][0]*LumiScale[mode]*BR_qq[mode]
        htemp["sig"] = hlist["sig_"+var][0]*LumiScale["sig"]*BR_qq["bb"]*2*Bs_had*Bs2TauTau_BR*Tau23Pi**2 #Scaling the sig is technically not necessary
        
        #Temporary solution to fit a realistic dataset (with the expected total number of events) while not all samples have been produced
        #Let's amplify the bb bkg (main source of the shape in any case) such that the background accounts for the total amount of Z
        #The ratio done here are not perfect: 1% stat evaluation + not exactly the same amount of ss and ud as bb
        #htemp["ss"] = BR_qq["ss"]/BR_qq["bb"]*3.754e-6/4.906895e-3*htemp["bb"]
        #htemp["ud"] = BR_qq["ss"]/BR_qq["bb"]*1.3868e-5/4.906895e-3*htemp["bb"]
        
        #Sum the backgrounds together to get a single shape
        BkgShape[var] = (htemp["bb"]+htemp["cc"]+htemp["ss"]+htemp["ud"],hlist["sig_"+var][1])
        
        #The signal
        SigShape[var] = (htemp["sig"],hlist["sig_"+var][1])

    return SigShape, BkgShape

#------------------------------------------------------------------------------------------------

#Scale the initial distr and Poisson vary it to simulate data to be fitted
def Make_Data(hlist,Vars,bkg_dummyEff,sig_dummyEff,seed):

    #Store the various sim data histo
    SimData = {}    

    for var in Vars:

        htemp = {}
        for mode in ["bb","cc","ss","ud"]:
            htemp[mode] = np.rint(hlist[mode+"_"+var][0]*LumiScale[mode]*BR_qq[mode]*bkg_dummyEff) #Round to the int, still a float type hope it works fine with the Poisson toy
        htemp["sig"] = np.rint(hlist["sig_"+var][0]*LumiScale["sig"]*BR_qq["bb"]*2*Bs_had*Bs2TauTau_BR*Tau23Pi**2*sig_dummyEff)

        #Temporary solution to fit a realistic dataset (with the expected total number of events) while not all samples have been produced
        #Let's amplify the bb bkg (main source of the shape in any case) such that the background accounts for the total amount of Z
        #The ratio done here are not perfect: 1% stat evaluation + not exactly the same amount of ss and ud as bb
        #htemp["ss"] = np.rint(BR_qq["ss"]/BR_qq["bb"]*3.754e-6/4.906895e-3*htemp["bb"])
        #htemp["ud"] = np.rint(BR_qq["ss"]/BR_qq["bb"]*1.3868e-5/4.906895e-3*htemp["bb"])

        if seed == 10: #Print only at the beginning for info
            print(f"Total Data Background = {np.sum(htemp['bb']+htemp['cc']+htemp['ss']+htemp['ud'])}")
            print(f"Total Data Signal = {np.sum(htemp['sig'])}")

        #Sum the histo together to get the simulated data distribution
        htemps = htemp["sig"]+htemp["bb"]+htemp["cc"]+htemp["ss"]+htemp["ud"]

        #Make the toy by Poisson varying it (give again the bin edges to keep the same data structure)
        np.random.seed(seed) #for reproducibility
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

    return popt, pcov, Chi2

#------------------------------------------------------------------------------------------------

def Draw_Toys(ValDist,Vars,BDTName):

    for var in Vars:

        #The expected value (true number of signal events used for toy gen) (to be sepecialised per var)
        Nexp = {}
        Nexp["sig"] = 503
        Nexp["bkg"] = 17568598

        #Compute the pulls
        Pulls_sig = (np.array(ValDist[var]["Best_Sig"])-Nexp["sig"])/np.array(ValDist[var]['Best_Sig']).std()
        Pulls_bkg = (np.array(ValDist[var]["Best_Bkg"])-Nexp["bkg"])/np.array(ValDist[var]['Best_Bkg']).std()
        
        fig, axs = plt.subplots(2,2)
        for ax in axs.flat:
            ax.grid(color='grey', linestyle='--', linewidth=0.5, alpha=0.5)

        hsig, edgsig = np.histogram(Pulls_sig, bins=19, range=[-5,5])
        hbkg, edgbkg = np.histogram(Pulls_bkg, bins=19, range=[-5,5])

        axs[0,0].stairs(hsig/len(ValDist[var]["Best_Sig"]),edges=edgsig,color="firebrick",fill=True)
        axs[0,1].hist(ValDist[var]["Sigma_Sig"],bins=20,color="firebrick")
        axs[1,0].stairs(hbkg/len(ValDist[var]["Best_Bkg"]),edges=edgbkg,color="steelblue",fill=True)
        axs[1,1].hist(ValDist[var]["Sigma_Bkg"],bins=20,color="steelblue")

        axs[0,0].set_xlabel(r"$\textrm{Signal Pulls}$",size="x-large")
        axs[0,1].set_xlabel(r"$\sigma_{N_{\rm sig}}$",size="x-large")
        axs[1,0].set_xlabel(r"$\textrm{Background Pulls}$",size="x-large")
        axs[1,1].set_xlabel(r"$\sigma_{N_{\rm bkg}}$",size="x-large")

        axs[0,0].text(0.5,1.05,r"\textrm{Best Fitted Value}",size="xx-large",transform=axs[0,0].transAxes,ha="center",va="center")
        axs[0,1].text(0.5,1.05,r"\textrm{Fit Uncertainty}",size="xx-large",transform=axs[0,1].transAxes,ha="center",va="center")
        #axs[0,0].text(-0.2,0.5,r"\textrm{Occurences}",size="xx-large",transform=axs[0,0].transAxes,ha="center",va="center",rotation="vertical")
        #axs[1,0].text(-0.2,0.5,r"\textrm{Occurences}",size="xx-large",transform=axs[1,0].transAxes,ha="center",va="center",rotation="vertical")
        axs[0,0].set_ylabel(r"$\textrm{Normalised Count}$",size="xx-large")
        axs[1,0].set_ylabel(r"$\textrm{Normalised Count}$",size="xx-large")

        axs[0,0].text(0.02,0.98,r"$N_{\rm sig}="+f"{int(np.array(ValDist[var]['Best_Sig']).mean())}"+r"\pm"+f"{int(np.array(ValDist[var]['Best_Sig']).std())}"+r"$",size="large",transform=axs[0,0].transAxes,ha="left",va="top")
        axs[0,1].text(0.02,0.98,r"$\mu="+f"{round(np.array(ValDist[var]['Sigma_Sig']).mean(),2)}"+r"$"+"\n"+r"$\sigma="+f"{round(np.array(ValDist[var]['Sigma_Sig']).std(),2)}"+r"$",size="large",transform=axs[0,1].transAxes,ha="left",va="top")
        axs[1,0].text(0.02,0.98,r"$N_{\rm bkg}="+f"{int(np.array(ValDist[var]['Best_Bkg']).mean())}"+r"\pm"+f"{int(np.array(ValDist[var]['Best_Bkg']).std())}"+r"$",size="large",transform=axs[1,0].transAxes,ha="left",va="top")
        axs[1,1].text(0.02,0.98,r"$\mu="+f"{round(np.array(ValDist[var]['Sigma_Bkg']).mean(),2)}"+r"$"+"\n"+r"$\sigma="+f"{round(np.array(ValDist[var]['Sigma_Bkg']).std(),2)}"+r"$",size="large",transform=axs[1,1].transAxes,ha="left",va="top")


        axs[0,0].text(1.1,1.2,r"\textrm{"+f"{var}"+r" Toys Results (}$N_{toys}="+f"{len(ValDist[var]['Best_Sig'])}"+r"$\textrm{)}",size="xx-large",transform=axs[0,0].transAxes,ha="center",va="center")

        fig.savefig(f"LumiFix_{BDTName}_{var}.pdf")    
            
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Collect info, scale and make toys, fits, draw according to a dumy eff to amplify the signal
def Do_Toys(Vars,BDTName,Bins,BkgEff,SigEff,NToy):

    print("Loading...")
    h = Load_Files_Hists(Vars,Bins,BDTName)
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
    Draw_Toys(Values,Vars,BDTName)

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
parser.add_argument("BDTName",type=str)
parser.add_argument("List",   nargs='+')
args=parser.parse_args()

Nbin = int((args.MHigh-args.MLow)/args.MassRes)

Do_Toys(args.List,args.BDTName,(Nbin,args.MLow,args.MHigh),args.BkgEff,args.SigEff,args.NToys)

