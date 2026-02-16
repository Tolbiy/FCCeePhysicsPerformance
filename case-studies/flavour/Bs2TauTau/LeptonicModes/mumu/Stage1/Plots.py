import ROOT as r
import numpy as np



def Load_RDF(mode):

    Path = "/eos/experiment/fcc/ee/analyses_storage/flavor/Bs2TauTau/flatNtuples/winter2023/analysis_stage1_mumu_noFilter/"

    if mode == "sig": 
        
        NF = 1
        filenames = r.std.vector('string')()
        for chunk in np.arange(0,20,1):
            filenames.push_back(Path+f"p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau/chunk_{chunk}.root")
        rdf = r.RDataFrame("events",filenames)
        #rdf = RDF_Treatment(rdf)
        return rdf
    
    else: 
        
        filenames = r.std.vector('string')()
        for chunk in np.arange(0,10,1):
            filenames.push_back(Path+f"p8_ee_Zbb_ecm91/chunk_{chunk}.root")
        rdf = r.RDataFrame("events",filenames)
        #rdf = RDF_Treatment(rdf)
        return rdf



def Do_Plots():

    sig = Load_RDF("sig")
    bb = Load_RDF("bb")

    print(f"Nsig = {sig.Count().GetValue()}")
    print(f"Nbb = {bb.Count().GetValue()}")


#=======================================================================================================

Do_Plots()