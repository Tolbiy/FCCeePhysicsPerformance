#from termcolor import colored
import json
import numpy as np

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

#preeff compute (efficiency coming from applying the FCCAnalyses layer which requires some Tau23PiCandidates)
#WARNING: There is no such thing as an implicit cuts due to the FCCAnalyses requiring Tau23PiCandidates: if none are find, the list is simply empty and the event is kept as such
#The events 'unselected' by the FCCAnalyses layer are those where the simulation failed for a reason or another (means that we should probably also disregard them)


#These preeff are to be taken into account, it corresponds to event with no PV (implicitly cutted out in the Stage1 script)
PreEff = {}
PreEff["bb"] = 4356115/(4*500000+6*400000)
PreEff["cc"] = 5088774/(600000+8*500000+494560)
PreEff["ss"] = 5060852/(492297+600000+7*500000+471230)
PreEff["ud"] = 4975473/(600000+7*500000+470159+408276)
PreEff["Bs2TauTau"] = 987681/(10*100000)
PreEff["TAUHADNU"] = 889007/(9*100000)

#A dict to link two ways of writting the bkg modes
Link = {"p8_ee_Zbb_ecm91":"bb",
        "p8_ee_Zcc_ecm91":"cc",
        "p8_ee_Zss_ecm91":"ss",
        "p8_ee_Zud_ecm91":"ud",
        "p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau":"Bs2TauTau",
        "p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU":"TAUHADNU"}

#Modes to be scanned
modes = ["p8_ee_Zbb_ecm91",
         "p8_ee_Zcc_ecm91",
         "p8_ee_Zss_ecm91",
         "p8_ee_Zud_ecm91",
         "p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU",
         #"p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau"
         ]

#print("\nPre-efficiencies (finding Tau23Pi candidates):")
#for i in modes:
#    print(f"{i}: {PreEff[Link[i]]}")

Previous = {}

with open("Cuts_Studied.json","r") as ofile:
    Previous = json.load(ofile)


print("\n---------Recaping-------------\n")
for CutsListed in Previous.keys():

    if CutsListed == "Before":
        continue
    
    print("-------------------------------")
    print(CutsListed.replace(" && ","\n"))
    print("-------------------------------\n")
    
    Nsig_Before = 0
    Nsig_After = 0
    Nbkg_Before = 0
    Nbkg_After = 0

    Nsig_exp = 0
    Nbkg_exp = 0

    for mode in modes:
        
        #ADD THE PREEFF (SMALL SELECTION BY REQUIRING TAU23PI CANDIDATES)
        eff = 1-(Previous['Before'][mode]-Previous[CutsListed][mode])/Previous['Before'][mode]
        #if eff > 0.01:
        #    Toprint = f"{100*round(eff,7)} %"
        #else:
        #  effnt = f"{round(eff,7)}"
            

        if mode == "p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU" or mode == "p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau":
            
            N_exp = int(NZ*BR_qq['bb']*2*Bs_had*Bs2TauTau_BR*Tau23Pi**2*eff*PreEff[Link[mode]])
            print(f"{mode} efficiency: "+f"{eff}"+f" ---> Nevents = {N_exp}")
            Nsig_Before += Previous["Before"][mode]
            Nsig_After  += Previous[CutsListed][mode]
            Nsig_exp += N_exp

        #elif mode not in Previous[CutsListed].keys():
        #    NBkg_Before += 0
        #    NBkg_After += 0

        else:
            
            N_exp = int(NZ*BR_qq[Link[mode]]*eff*PreEff[Link[mode]])
            print(f"{mode} efficiency: "+f"{eff}"+f" ---> Nevents = {N_exp:e}")
            Nbkg_Before += Previous["Before"][mode]
            Nbkg_After  += Previous[CutsListed][mode]
            Nbkg_exp += N_exp
    
    print("\nSignal efficiency: ",f"{round(100-100*(Nsig_Before-Nsig_After)/Nsig_Before,5)} %")
    print("Background efficiency: ",f"{round(100-100*(Nbkg_Before-Nbkg_After)/Nbkg_Before,5)} %")
    print(f"\nSignificance S/sqrt(B) = {round(Nsig_exp/2/np.sqrt(Nbkg_exp),3)}\n") #Divide by 2 Nsig because the two signal samples are actually the same
