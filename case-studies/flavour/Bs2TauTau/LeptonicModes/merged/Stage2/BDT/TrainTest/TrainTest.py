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

from matplotlib import rc
rc('font',**{'family':'serif','serif':['Roman']})
rc('text', usetex=True)

#Local code
#from userConfig import loc, mode, train_vars, train_vars_vtx, mode_names
#import plotting
#import utils as ut

def TrainTest_Samples(vars):

    path = "../PreTreatedData"
    modes = ["sig","bb","cc"]
    dfs = {}
    train = {}
    test = {}   

    print("Loading PreTreated DF...")
    for mode in modes:

        dfs[mode] = uproot.open(f"{path}/{mode}/Naive_withMoreData.root:events").arrays(library="pd")#.sample(n=10000,random_state=12)
        dfs[mode] = dfs[mode][vars]
        if mode == "sig":
            dfs[mode]["label"] = 1
        else:
            dfs[mode]["label"] = 0


        train[mode], test[mode] = train_test_split(dfs[mode], test_size=0.2, random_state=911)

    return train, test

#_____________________________________________________________________________________________________________________________________________

def Train(train):

    print("Start Training")
    #Regroupe all modes to train the BDT
    train_tot = pd.concat([train[mode] for mode in ["sig","bb","cc"]])

    vars_list = list(train_tot.columns.values)[:-1]
    print(f"Variables used:\n {vars_list}")

    #Split into class label (y) and training vars (x)
    y = train_tot["label"]
    x = train_tot[vars_list]

    y = y.to_numpy()
    x = x.to_numpy()

    #Sample weights to balance the classes
    weights = compute_sample_weight(class_weight='balanced', y=y)

    #BDT
    config_dict = {
            "n_estimators": 1000,
            "learning_rate": 0.3,
            "max_depth": 5,
            "min_child_weight": 50,
            }

    bdt = xgb.XGBClassifier(n_estimators=config_dict["n_estimators"],
                            max_depth=config_dict["max_depth"],
                            learning_rate=config_dict["learning_rate"],
                            min_child_weight=config_dict["min_child_weight"]
                            )

    #Fit the model
    print("Training model")
    bdt.fit(x, y, sample_weight=weights)

    feature_importances = pd.DataFrame(bdt.feature_importances_,
                                     index = vars_list,
                                     columns=['importance']).sort_values('importance',ascending=False)

    print("Feature importances")
    print(feature_importances)
    feature_importances.to_json("Train_Results/Feature/Naive_withMoreData.json")

    #Create ROC curves
    decisions = bdt.predict_proba(x)[:,1]

    # Compute ROC curves and area under the curve
    fpr, tpr, thresholds = roc_curve(y, decisions)
    roc_auc = auc(fpr, tpr)

    fig, ax = plt.subplots(figsize=(8,8))
    plt.plot(tpr, 1-fpr, lw=1.5, color="k", label='ROC (area = %0.3f)'%(roc_auc))
    plt.plot([0.45, 1.], [0.45, 1.], linestyle="--", color="k", label='50/50')
    plt.xlim(0.45,1.)
    plt.ylim(0.45,1.)
    plt.ylabel('Background rejection',fontsize=30)
    plt.xlabel('Signal efficiency',fontsize=30)
    ax.tick_params(axis='both', which='major', labelsize=25)
    plt.legend(loc="upper left",fontsize=20)
    plt.grid()
    plt.tight_layout()
    fig.savefig(f"Train_Results/ROC/Naive_withMoreData.pdf")

    print("Writting BDT model")
    #Write it for additional testing
    joblib.dump(bdt, f"Train_Results/Models/Naive_withMoreData.joblib")
    
    #Write the model to a ROOT file on EOS, for application elsewhere in FCCAnalyses
    ROOT.TMVA.Experimental.SaveXGBoost(bdt, "Naive_withMoreData", f"Train_Results/Models/Naive_withMoreData.root", num_inputs=len(vars_list)) 
    #To add it to the dataset column will need these columns -> add them in the Stage2 script

#_____________________________________________________________________________________________________________________________________________  
    
def Test(train,test):    
    
    print("Start Testing")

    vars_list = list(train["sig"].columns.values)[:-1]
    print(f"Variables used: {vars_list}")
    
    #Get the correlation matrix
    fig2, ax2 = plt.subplots(figsize=(8,5), dpi=80)
    fig2.colorbar(ax2.matshow(pd.concat([train[mode][vars_list] for mode in ["sig","bb"]]).corr(),vmin=-1.0,vmax=1.0),label="Correlation")
    ax2.matshow(pd.concat([train[mode][vars_list] for mode in ["sig","bb"]]).corr(),vmin=-1.0,vmax=1.0)
    ax2.set_xticks(ticks=np.arange(0,len(vars_list),1),labels=vars_list,rotation=90,size="small")
    ax2.set_yticks(ticks=np.arange(0,len(vars_list),1),labels=vars_list,size="small")
    fig2.savefig("Train_Results/Correlation/Naive_withMoreData.pdf")

    bdt = joblib.load(f"Train_Results/Models/Naive_withMoreData.joblib")
    
    #Train-Test comparison

    for mode in ["sig","bb","cc"]:
        train[mode]["BDT"] = bdt.predict_proba(train[mode][vars_list]).tolist()
        train[mode]["BDT"] = train[mode]["BDT"].apply(lambda x: x[1])

        test[mode]["BDT"] = bdt.predict_proba(test[mode][vars_list]).tolist()
        test[mode]["BDT"] = test[mode]["BDT"].apply(lambda x: x[1])

    
    cuts = np.linspace(0.0,1.0,100)
    Eff = {}
    
    #Collect the efficiencies
    for mode in ["sig","bb","cc"]:

            eff_train = []
            eff_test = []
            for cut in cuts: 
                eff_train.append(len(train[mode].query(f"BDT > {cut}").index)/len(train[mode].index))
                eff_test.append(len(test[mode].query(f"BDT > {cut}").index)/len(test[mode].index))
            Eff[f"{mode}_train"] = eff_train
            Eff[f"{mode}_test"] = eff_test

    #Start drawing

    colors = {"sig":"royalblue",
              "bb": "indigo",
              "cc": "darkorchid",
              "ss": "darkviolet",
              "ud": "mediumorchid",
             }

    linstyle = {"train":"-","test":"--"}

    labels = {"sig":r"$B_s^0\rightarrow \tau^+\tau^- (\tau\rightarrow\ell\nu_{\tau}\nu_{\ell})",
              "bb": r"$Z^0\rightarrow b\overline{b}",
              "cc": r"$Z^0\rightarrow c\overline{c}",
              "ss": r"$Z^0\rightarrow s\overline{s}",
              "ud": r"$Z^0\rightarrow u\overline{d}",
             }

    fig, ax = plt.subplots()
    ax.grid(color='grey', linestyle='--', linewidth=0.5, alpha=0.5)
    
    for key in Eff.keys():
        mode, cat = key.split("_")
        ax.plot(cuts,Eff[key],ls=linstyle[cat],color=colors[mode],label=labels[mode]+r"\textrm{ "+f"{cat}"+r"}$")
    
    ax.set_xlabel("Stage 2 BDT cuts")
    ax.set_ylabel("Efficiencies")
    ax.set_xlim([0.0,1.0])
    ax.set_ylim([1e-4,1.05])
    ax.set_yscale("log")

    #Legend
    ax.legend(frameon=True, framealpha=1, fancybox=True, edgecolor='lightgrey', loc="center left", bbox_to_anchor=(0.1,0.2),ncol=2)

    fig.savefig("Train_Results/Overtrain/Naive_withMoreData.pdf")
    




def main():
    parser = argparse.ArgumentParser(description='Train xgb model for Bs -> Tau Tau (Tau -> lnunu (l=e/mu)) Stage 2 (focuses on identified dilepton systems)')
    parser.add_argument("--TrainTest", "-t", choices=["Train","Test"],required=False,help="Train: Train and Test BDT, Test: Test only",default="Train")
    args = parser.parse_args()

    #Select a subset of the variables
    VarSet = [#Muon related
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

    train, test = TrainTest_Samples(VarSet)

    if args.TrainTest == "Train":
        Train(train)
        Test(train,test)
    elif args.TrainTest == "Test":
        Test(train,test)

if __name__ == '__main__':
    main()