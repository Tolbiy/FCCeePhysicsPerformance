import sys,os, argparse
import json
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.utils.class_weight import compute_sample_weight
from sklearn.metrics import roc_curve, auc
#from root_pandas import read_root
import uproot
import ROOT
import joblib
import glob
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

from matplotlib import rc
rc('font',**{'family':'serif','serif':['Roman']})
rc('text', usetex=True)

#Local code
#from userConfig import loc, mode, train_vars, train_vars_vtx, mode_names
#import plotting
#import utils as ut

def TrainTest_Samples(vars):

    path = "/afs/cern.ch/work/t/tomonnar/public/Bs2TauTau/Stage2/BDT/PreTreatedFiles/Baseline_NoBug/"
    modes = ["sig","bb","cc","ss","ud"]
    dfs = {}
    train = {}
    test = {}   

    print("Loading PreTreated DF...")
    for mode in modes:
        

        if mode == "sig":

            dfs[mode] = uproot.open(f"{path}{mode}.root:events").arrays(library="pd").sample(n=100000,random_state=12) #Select 100k events to have as much sig as background
            dfs[mode] = dfs[mode][vars]
            dfs[mode]["label"] = 1

        else:

            dfs[mode] = uproot.open(f"{path}{mode}.root:events").arrays(library="pd")
            dfs[mode] = dfs[mode][vars]
            dfs[mode]["label"] = 0

        

        train[mode], test[mode] = train_test_split(dfs[mode], test_size=0.2, random_state=911)

    return train, test

#_____________________________________________________________________________________________________________________________________________

def Train(train):

    print("Start Training")
    #Regroupe all modes to train the BDT
    train_tot = pd.concat([train[mode] for mode in ["sig","bb","cc","ss","ud"]])

    vars_list = list(train_tot.columns.values)[:-1]
    print(f"Variables used: {vars_list}")

    #Split into class label (y) and training vars (x)
    y = train_tot["label"]
    x = train_tot[vars_list]

    y = y.to_numpy()
    x = x.to_numpy()

    #Sample weights to balance the classes
    weights = compute_sample_weight(class_weight='balanced', y=y)

    #BDT
    config_dict = {
            "n_estimators": 400,
            "learning_rate": 0.3,
            "max_depth": 3,
            }

    bdt = xgb.XGBClassifier(n_estimators=config_dict["n_estimators"],
                            max_depth=config_dict["max_depth"],
                            learning_rate=config_dict["learning_rate"],
                            )

    #Fit the model
    print("Training model")
    bdt.fit(x, y, sample_weight=weights)

    feature_importances = pd.DataFrame(bdt.feature_importances_,
                                     index = vars_list,
                                     columns=['importance']).sort_values('importance',ascending=False)

    print("Feature importances")
    print(feature_importances)
    feature_importances.to_json("Feature/feature_importances_NoMass_NoBug.json")

    print("Writting BDT model")
    #Write it for additional testing
    joblib.dump(bdt, f"/afs/cern.ch/work/t/tomonnar/public/Bs2TauTau/Stage2/BDT/Models/xgb_bdt_Baseline_NoMass_NoBug.joblib")
    
    #Write the model to a ROOT file on EOS, for application elsewhere in FCCAnalyses
    ROOT.TMVA.Experimental.SaveXGBoost(bdt, "Bs2TauTau_Stage2_BDT", f"/afs/cern.ch/work/t/tomonnar/public/Bs2TauTau/Stage2/BDT/Models/xgb_bdt_Baseline_NoMass_NoBug.root", num_inputs=len(vars_list)) 
    #To add it to the dataset column will need these columns -> add them in the Stage2 script

#_____________________________________________________________________________________________________________________________________________  
    
def Test(train,test):    
    
    print("Start Testing")

    vars_list = list(train["sig"].columns.values)[:-1]
    print(f"Variables used: {vars_list}")
    
    #Get the correlation matrix
    fig2, ax2 = plt.subplots(figsize=(8,5), dpi=80)
    fig2.colorbar(ax2.matshow(pd.concat([train[mode][vars_list] for mode in ["sig","bb","cc","ss","ud"]]).corr(),vmin=-1.0,vmax=1.0),label="Correlation")
    ax2.matshow(pd.concat([train[mode][vars_list] for mode in ["sig","bb","cc","ss","ud"]]).corr(),vmin=-1.0,vmax=1.0)
    ax2.set_xticks(ticks=np.arange(0,len(vars_list),1),labels=vars_list,rotation=90,size="small")
    ax2.set_yticks(ticks=np.arange(0,len(vars_list),1),labels=vars_list,size="small")
    fig2.savefig("Feature/Correlation_Baseline_NoBug.pdf")

    bdt = joblib.load(f"/afs/cern.ch/work/t/tomonnar/public/Bs2TauTau/Stage2/BDT/Models/xgb_bdt_Baseline_NoBug.joblib")
    
    #Split into class label (y) and training vars (x)
    #Regroupe all modes to train the BDT
    train_tot = pd.concat([train[mode] for mode in ["sig","bb","cc","ss","ud"]])
    y = train_tot["label"]
    x = train_tot[vars_list]
    y = y.to_numpy()
    x = x.to_numpy()

    #Create ROC curves
    decisions = bdt.predict_proba(x)[:,1]

    # Compute ROC curves and area under the curve
    fpr, tpr, thresholds = roc_curve(y, decisions)
    roc_auc = auc(fpr, tpr)

    fig, ax = plt.subplots(figsize=(8,8))
    plt.plot(fpr, tpr, lw=1.5, color="k", label='ROC (area = %0.3f)'%(roc_auc))
    plt.plot([0., 1.], [0., 1.], linestyle="--", color="k", label='50/50')
    plt.xlim(0.,1.)
    plt.ylim(0.,1.)
    plt.ylabel('True positive rate',fontsize=30)
    plt.xlabel('False positive rate',fontsize=30)
    ax.tick_params(axis='both', which='major', labelsize=25)
    plt.legend(loc="lower left",fontsize=20)
    plt.grid()
    plt.tight_layout()
    fig.savefig(f"ROC/Stage2_BDT_Baseline_NoBug.pdf")



    #Train-Test comparison

    for mode in ["sig","bb","cc","ss","ud"]:
        train[mode]["BDT"] = bdt.predict_proba(train[mode][vars_list]).tolist()
        train[mode]["BDT"] = train[mode]["BDT"].apply(lambda x: x[1])

        test[mode]["BDT"] = bdt.predict_proba(test[mode][vars_list]).tolist()
        test[mode]["BDT"] = test[mode]["BDT"].apply(lambda x: x[1])

    train["bkg"] = pd.concat([train[mode] for mode in ["bb","cc","ss","ud"]])
    test["bkg"] = pd.concat([test[mode] for mode in ["bb","cc","ss","ud"]])

    BDTHist = {}
    BDTHist["train_sig"] = np.histogram(train["sig"]["BDT"].to_numpy(),bins=50,range=[0,1])
    BDTHist["train_bkg"] = np.histogram(train["bkg"]["BDT"].to_numpy(),bins=50,range=[0,1])
    BDTHist["test_sig"] = np.histogram(test["sig"]["BDT"].to_numpy(),bins=50,range=[0,1])
    BDTHist["test_bkg"] = np.histogram(test["bkg"]["BDT"].to_numpy(),bins=50,range=[0,1])

    NormHeight = {}
    NormHeight["train_sig"] = np.array(BDTHist["train_sig"][0])/len(train["sig"]["BDT"].to_numpy())
    NormHeight["train_bkg"] = np.array(BDTHist["train_bkg"][0])/len(train["bkg"]["BDT"].to_numpy())
    NormHeight["test_sig"] = np.array(BDTHist["test_sig"][0])/len(test["sig"]["BDT"].to_numpy())
    NormHeight["test_bkg"] = np.array(BDTHist["test_bkg"][0])/len(test["bkg"]["BDT"].to_numpy())

    fig, ax = plt.subplots()
    ax.grid(color='grey', linestyle='--', linewidth=0.5, alpha=0.5)

    ax.stairs(NormHeight["train_sig"],BDTHist["train_sig"][1],ls="-",color="firebrick",label="Sig. training")
    ax.stairs(NormHeight["test_sig"],BDTHist["train_sig"][1],ls="--",color="firebrick",label="Sig. testing")
    ax.stairs(NormHeight["train_bkg"],BDTHist["train_sig"][1],ls="-",color="steelblue",label="Bkg. training")
    ax.stairs(NormHeight["test_bkg"],BDTHist["train_sig"][1],ls="--",color="steelblue",label="Bkg. testing")

    ax.set_xlabel(r"$\textrm{BDT Score}$")
    ax.set_ylabel(r"$\textrm{Normalised Counts}$")
    ax.set_yscale("log")
    ax.legend()
    fig.savefig("Overtrain/Pres_Baseline_NoBug.pdf")

    '''
    cuts = np.linspace(0.0,1.0,100)
    Eff = {}
    
    #Collect the efficiencies
    for mode in ["sig","bb","cc","ss","ud"]:

            eff_train = []
            eff_test = []
            for cut in cuts: 
                eff_train.append(len(train[mode].query(f"BDT > {cut}").index)/len(train[mode].index))
                eff_test.append(len(test[mode].query(f"BDT > {cut}").index)/len(test[mode].index))
            Eff[f"{mode}_train"] = eff_train
            Eff[f"{mode}_test"] = eff_test

    #Start drawing

    colors = {"sig":"firebrick",
              "bb": "steelblue",
              "cc": "goldenrod",
              "ss": "mediumseagreen",
              "ud": "slategrey",
             }

    linstyle = {"train":"-","test":"--"}

    labels = {"sig":r"$B_s^0\rightarrow \tau^+\tau^- (\tau\rightarrow3\pi^{\pm}\nu_{\tau})$",
              "bb": r"$Z^0\rightarrow b\overline{b}$",
              "cc": r"$Z^0\rightarrow c\overline{c}$",
              "ss": r"$Z^0\rightarrow s\overline{s}$",
              "ud": r"$Z^0\rightarrow u\overline{d}$",
             }

    fig, ax = plt.subplots()
    ax.grid(color='grey', linestyle='--', linewidth=0.5, alpha=0.5)
    
    for key in Eff.keys():
        mode, cat = key.split("_")
        ax.plot(cuts,Eff[key],ls=linstyle[cat],color=colors[mode],label=labels[mode]+r"\textrm{ "+f"{cat}"+r"}$")
    
    ax.set_xlabel("Stage 2 BDT cuts")
    ax.set_ylabel("Efficiencies")
    ax.set_xlim([0.0,1.0])
    ax.set_ylim([1e-3,1.05])
    ax.set_yscale("log")

    #Legend
    #ax.legend(frameon=True, framealpha=1, fancybox=True, edgecolor='lightgrey', loc="center left", bbox_to_anchor=(0.1,0.2),ncol=2)

    #Legend
    testlab = Line2D([0,0], [0,1], label=r'$\textrm{Test}$', ls="--", color='k')
    trainlab = Line2D([0,0], [0,1], label=r'$\textrm{Train}$', ls="-", color='k')
    siglab = Line2D([0,0],[0,1], label=labels["sig"], ls="-",color=colors["sig"])
    bblab = Line2D([0,0],[0,1], label=labels["bb"], ls="-",color=colors["bb"])
    cclab = Line2D([0,0],[0,1], label=labels["cc"], ls="-",color=colors["cc"])
    sslab = Line2D([0,0],[0,1], label=labels["ss"], ls="-",color=colors["ss"])
    udlab = Line2D([0,0],[0,1], label=labels["ud"], ls="-",color=colors["ud"])
    handlist = [trainlab,testlab,siglab,bblab,cclab,sslab,udlab]
    ax.legend(handles=handlist, frameon=True, framealpha=1, fancybox=True, edgecolor='lightgrey', loc="lower left")

    fig.savefig("Overtrain/TrainTest_Baseline_NoBug.pdf")
    '''
    




def main():
    parser = argparse.ArgumentParser(description='Train xgb model for Bs -> Tau Tau (Tau -> 3pi) Stage 2 (focuses on identified di-Tau systems)')
    parser.add_argument("--TrainTest", "-t", choices=["Train","Test"],required=False,help="Train: Train and Test BDT, Test: Test only",default="Train")
    args = parser.parse_args()

    #Select a subset of the variables
    VarSet = [
              "diTau_Angles","diTauPlus_IP","diTauPlus_IPV","diTauMinus_IP","diTauMinus_IPV","diTauMinus_rho1mass","diTauMinus_rho2mass",
              "diTauPlus_rho1mass","diTauPlus_rho2mass","diTauPlus_mass","diTauMinus_mass","diTauPlus_Lifetime","diTauMinus_Lifetime","Bs_Lifetime","Bs_IPV","mDiTau_Vis","DeltaM"
             ]

    train, test = TrainTest_Samples(VarSet)

    if args.TrainTest == "Train":
        Train(train)
        Test(train,test)
    elif args.TrainTest == "Test":
        Test(train,test)

if __name__ == '__main__':
    main()