import numpy as np

#Some numbers to evaluate S and B
NZ = 6e12 #Number of Z (To be refined depending on FCC scenarios)
Bs_had = 0.1 #Ratio of b->Bs hadronisation (given by Xunwu, to be refined)
Bs2TauTau_BR = 7.73e-7 #SM estimate for Bs->TauTau (from the LHCb paper that searched for it)

Tau2e = 0.1785**2

#- BR of Z->qq ----------------------------------------------------------------
BR_qq = {}
BR_qq["bb"] = 0.1512
BR_qq["cc"] = 0.1203
BR_qq["ss"] = 0.1560
BR_qq["ud"] = 0.6991 - BR_qq["bb"] - BR_qq["cc"] - BR_qq["ss"]

Total = {}
Total["sig"] = 320059
Total["bb"]  = 4356115
Total["cc"]  = 5088774
Total["ss"]  = 5060852
Total["ud"]  = 4975473

Pass = {}
Pass["sig"] = 129481
Pass["bb"]  = 5276
Pass["cc"]  = 16
Pass["ss"]  = 4
Pass["ud"]  = 3

PreBDTEff = {}
for mode in Total.keys():
    PreBDTEff[mode] = Pass[mode]/Total[mode]

PreBDTEffUnc = {}
for mode in Total.keys():
    PreBDTEffUnc[mode] = np.sqrt(PreBDTEff[mode]*(1-PreBDTEff[mode])/Total[mode])
    print(f"{mode}: {PreBDTEff[mode]} pm {PreBDTEffUnc[mode]}")

print("\n")

LumiScale = {}
for mode in Total.keys():
    if mode == "sig":
        LumiScale[mode] = NZ*BR_qq["bb"]*Bs_had*Bs2TauTau_BR*Tau2e*PreBDTEff[mode]
    else:
        LumiScale[mode] = NZ*BR_qq[mode]*PreBDTEff[mode]

LumiScaleUnc = {}
for mode in Total.keys():
    if mode == "sig":    
        LumiScaleUnc[mode] = NZ*BR_qq["bb"]*Bs_had*Bs2TauTau_BR*Tau2e*PreBDTEffUnc[mode]
    else:
        LumiScaleUnc[mode] = NZ*BR_qq[mode]*PreBDTEffUnc[mode]
    print(f"{mode}: {LumiScale[mode]} pm {LumiScaleUnc[mode]}")