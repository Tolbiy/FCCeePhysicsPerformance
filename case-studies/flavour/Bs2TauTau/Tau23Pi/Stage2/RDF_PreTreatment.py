import ROOT as r
import numpy as np
import sys

sys.path.insert(1,"../Utils")
import Functions_ForPreTreat

Links = {"bb": "p8_ee_Zbb_ecm91",
         "cc": "p8_ee_Zcc_ecm91",
         "ss": "p8_ee_Zss_ecm91",
         "ud": "p8_ee_Zud_ecm91",
         "sig":"p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU"
        }

def RDF_Treatment(rdf):

    #Compute the Tau momenta angle and ID them and find the di-Tau, diTau for check only stage1 cuts requires only two tau candidates (split the pAngles column by Taus charge, or more usefully by TauCand to avoid mixing the two selections ways)
    rdf = rdf.Define("Tau23PiCandidates_pAngles","Compute_Momenta_Angles(EVT_NTau23Pi,Tau23PiCandidates_px,Tau23PiCandidates_py,Tau23PiCandidates_pz)")\
             .Define("Tau23PiCandidates_pAngles_ID","ID_Angles(EVT_NTau23Pi)")\
             .Define("diTau","Find_diTau(EVT_NTau23Pi,Tau23PiCandidates_pAngles,Tau23PiCandidates_pAngles_ID,Tau23PiCandidates_q)")\
             .Define("n_diTau","diTau.size()")\
             .Filter("n_diTau > 1") #Get rid of the the phantom events (where no diTau where found)

    #Get the vertex position and momentum of the diTaus (disentangled by charge)
    rdf = rdf.Define("diTauPlus_x", "if (Tau23PiCandidates_q[diTau[0]] > 0) return Tau23PiCandidates_x[diTau[0]];  else return Tau23PiCandidates_x[diTau[1]];")\
             .Define("diTauPlus_y", "if (Tau23PiCandidates_q[diTau[0]] > 0) return Tau23PiCandidates_y[diTau[0]];  else return Tau23PiCandidates_y[diTau[1]];")\
             .Define("diTauPlus_z", "if (Tau23PiCandidates_q[diTau[0]] > 0) return Tau23PiCandidates_z[diTau[0]];  else return Tau23PiCandidates_z[diTau[1]];")\
             .Define("diTauPlus_px","if (Tau23PiCandidates_q[diTau[0]] > 0) return Tau23PiCandidates_px[diTau[0]]; else return Tau23PiCandidates_px[diTau[1]];")\
             .Define("diTauPlus_py","if (Tau23PiCandidates_q[diTau[0]] > 0) return Tau23PiCandidates_py[diTau[0]]; else return Tau23PiCandidates_py[diTau[1]];")\
             .Define("diTauPlus_pz","if (Tau23PiCandidates_q[diTau[0]] > 0) return Tau23PiCandidates_pz[diTau[0]]; else return Tau23PiCandidates_pz[diTau[1]];")\
             .Define("diTauPlus_p", "if (Tau23PiCandidates_q[diTau[0]] > 0) return Tau23PiCandidates_p[diTau[0]];  else return Tau23PiCandidates_p[diTau[1]];")\
             .Define("diTauPlus_rho1mass", "if (Tau23PiCandidates_q[diTau[0]] > 0) return Tau23PiCandidates_rho1mass[diTau[0]];  else return Tau23PiCandidates_rho1mass[diTau[1]];")\
             .Define("diTauPlus_rho2mass", "if (Tau23PiCandidates_q[diTau[0]] > 0) return Tau23PiCandidates_rho2mass[diTau[0]];  else return Tau23PiCandidates_rho2mass[diTau[1]];")\
             .Define("diTauPlus_mass","if (Tau23PiCandidates_q[diTau[0]] > 0) return Tau23PiCandidates_mass[diTau[0]];  else return Tau23PiCandidates_mass[diTau[1]];")\
             .Define("diTauPlus_chi2","if (Tau23PiCandidates_q[diTau[0]] > 0) return Tau23PiCandidates_chi2[diTau[0]];  else return Tau23PiCandidates_chi2[diTau[1]];")\
             \
             .Define("diTauMinus_x", "if (Tau23PiCandidates_q[diTau[0]] < 0) return Tau23PiCandidates_x[diTau[0]];  else return Tau23PiCandidates_x[diTau[1]];")\
             .Define("diTauMinus_y", "if (Tau23PiCandidates_q[diTau[0]] < 0) return Tau23PiCandidates_y[diTau[0]];  else return Tau23PiCandidates_y[diTau[1]];")\
             .Define("diTauMinus_z", "if (Tau23PiCandidates_q[diTau[0]] < 0) return Tau23PiCandidates_z[diTau[0]];  else return Tau23PiCandidates_z[diTau[1]];")\
             .Define("diTauMinus_px","if (Tau23PiCandidates_q[diTau[0]] < 0) return Tau23PiCandidates_px[diTau[0]]; else return Tau23PiCandidates_px[diTau[1]];")\
             .Define("diTauMinus_py","if (Tau23PiCandidates_q[diTau[0]] < 0) return Tau23PiCandidates_py[diTau[0]]; else return Tau23PiCandidates_py[diTau[1]];")\
             .Define("diTauMinus_pz","if (Tau23PiCandidates_q[diTau[0]] < 0) return Tau23PiCandidates_pz[diTau[0]]; else return Tau23PiCandidates_pz[diTau[1]];")\
             .Define("diTauMinus_p", "if (Tau23PiCandidates_q[diTau[0]] < 0) return Tau23PiCandidates_p[diTau[0]];  else return Tau23PiCandidates_p[diTau[1]];")\
             .Define("diTauMinus_rho1mass", "if (Tau23PiCandidates_q[diTau[0]] < 0) return Tau23PiCandidates_rho1mass[diTau[0]];  else return Tau23PiCandidates_rho1mass[diTau[1]];")\
             .Define("diTauMinus_rho2mass", "if (Tau23PiCandidates_q[diTau[0]] < 0) return Tau23PiCandidates_rho2mass[diTau[0]];  else return Tau23PiCandidates_rho2mass[diTau[1]];")\
             .Define("diTauMinus_mass","if (Tau23PiCandidates_q[diTau[0]] < 0) return Tau23PiCandidates_mass[diTau[0]];  else return Tau23PiCandidates_mass[diTau[1]];")\
             .Define("diTauMinus_chi2","if (Tau23PiCandidates_q[diTau[0]] < 0) return Tau23PiCandidates_chi2[diTau[0]];  else return Tau23PiCandidates_chi2[diTau[1]];")\
             \
             .Define("diTau_Angles","Compute_Momenta_Angles(diTauPlus_x,diTauPlus_y,diTauPlus_z,diTauMinus_x,diTauMinus_y,diTauMinus_z)")
             
    
    #Find the RECO (Visible) vertex and get its "quality", given by the length of the line on which the Visible vertex is found (not sure yet how to properly use it to assess the actual quality, would need a reference distance)
    rdf = rdf.Define("Bs2TauTau",  "Compute_BsVisibleVertex(diTauPlus_x,diTauPlus_y,diTauPlus_z,diTauMinus_x,diTauMinus_y,diTauMinus_z,diTauPlus_px,diTauPlus_py,diTauPlus_pz,diTauMinus_px,diTauMinus_py,diTauMinus_pz)")\
             .Define("Bs2TauTau_x","Bs2TauTau[0]")\
             .Define("Bs2TauTau_y","Bs2TauTau[1]")\
             .Define("Bs2TauTau_z","Bs2TauTau[2]")

    #Get the PV coordinates 
    rdf = rdf.Define("PV_x","for (size_t i=0;i<Vertex_isPV.size();++i) if (Vertex_isPV[i]) return Vertex_x[i];")\
             .Define("PV_y","for (size_t i=0;i<Vertex_isPV.size();++i) if (Vertex_isPV[i]) return Vertex_y[i];")\
             .Define("PV_z","for (size_t i=0;i<Vertex_isPV.size();++i) if (Vertex_isPV[i]) return Vertex_z[i];")

    #Get the IPs of both taus
    rdf = rdf.Define("diTauPlus_IPP","Compute_IPP(diTauPlus_x,diTauPlus_y,diTauPlus_z,Bs2TauTau_x,Bs2TauTau_y,Bs2TauTau_z,PV_x,PV_y,PV_z)")\
             .Define("diTauPlus_IPx","diTauPlus_IPP[0]-diTauPlus_x")\
             .Define("diTauPlus_IPy","diTauPlus_IPP[1]-diTauPlus_y")\
             .Define("diTauPlus_IPz","diTauPlus_IPP[2]-diTauPlus_z")\
             .Define("diTauPlus_IP", "Get_Norm(diTauPlus_IPx,diTauPlus_IPy,diTauPlus_IPz)")\
             .Define("diTauMinus_IPP","Compute_IPP(diTauMinus_x,diTauMinus_y,diTauMinus_z,Bs2TauTau_x,Bs2TauTau_y,Bs2TauTau_z,PV_x,PV_y,PV_z)")\
             .Define("diTauMinus_IPx","diTauMinus_IPP[0]-diTauMinus_x")\
             .Define("diTauMinus_IPy","diTauMinus_IPP[1]-diTauMinus_y")\
             .Define("diTauMinus_IPz","diTauMinus_IPP[2]-diTauMinus_z")\
             .Define("diTauMinus_IP", "Get_Norm(diTauMinus_IPx,diTauMinus_IPy,diTauMinus_IPz)")

    #Get the Flight Distance and Distance
    rdf = rdf.Define("diTauPlus_FlightDistance", "sqrt(pow(diTauPlus_x-Bs2TauTau_x,2) + pow(diTauPlus_y-Bs2TauTau_y,2) + pow(diTauPlus_z-Bs2TauTau_z,2))")\
             .Define("diTauPlus_Lifetime",       "1.7769/diTauPlus_p*diTauPlus_FlightDistance/299792458e3")\
             .Define("diTauMinus_FlightDistance","sqrt(pow(diTauMinus_x-Bs2TauTau_x,2) + pow(diTauMinus_y-Bs2TauTau_y,2) + pow(diTauMinus_z-Bs2TauTau_z,2))")\
             .Define("diTauMinus_Lifetime",      "1.7769/diTauMinus_p*diTauMinus_FlightDistance/299792458e3")

    #Get the Tau Vertex IP (Distance between vertex and z-oriented line passing by PV, EXACT) and the "usual" Tau IP (Distance between Taus flight dir and PV, BIASED)
    rdf = rdf.Define("diTauPlus_IPVP","Compute_IPP(diTauPlus_x,diTauPlus_y,diTauPlus_z,PV_x,PV_y,PV_z+1.0,PV_x,PV_y,PV_z)")\
             .Define("diTauPlus_IPVx","diTauPlus_IPVP[0]-diTauPlus_x")\
             .Define("diTauPlus_IPVy","diTauPlus_IPVP[1]-diTauPlus_y")\
             .Define("diTauPlus_IPVz","diTauPlus_IPVP[2]-diTauPlus_z")\
             .Define("diTauPlus_IPV","Get_Norm(diTauPlus_IPVx,diTauPlus_IPVy,diTauPlus_IPVz)")\
             .Define("diTauMinus_IPVP","Compute_IPP(diTauMinus_x,diTauMinus_y,diTauMinus_z,PV_x,PV_y,PV_z+1.0,PV_x,PV_y,PV_z)")\
             .Define("diTauMinus_IPVx","diTauMinus_IPVP[0]-diTauMinus_x")\
             .Define("diTauMinus_IPVy","diTauMinus_IPVP[1]-diTauMinus_y")\
             .Define("diTauMinus_IPVz","diTauMinus_IPVP[2]-diTauMinus_z")\
             .Define("diTauMinus_IPV","Get_Norm(diTauMinus_IPVx,diTauMinus_IPVy,diTauMinus_IPVz)")

    rdf = rdf.Define("diTauPlus_IPFP","Compute_IPP(PV_x,PV_y,PV_z, diTauPlus_x,diTauPlus_y,diTauPlus_z, diTauPlus_x+diTauPlus_px,diTauPlus_y+diTauPlus_py,diTauPlus_z+diTauPlus_pz)")\
             .Define("diTauPlus_IPFx","diTauPlus_IPFP[0]-PV_x")\
             .Define("diTauPlus_IPFy","diTauPlus_IPFP[1]-PV_y")\
             .Define("diTauPlus_IPFz","diTauPlus_IPFP[2]-PV_z")\
             .Define("diTauPlus_IPF","Get_Norm(diTauPlus_IPFx,diTauPlus_IPFy,diTauPlus_IPFz)")\
             .Define("diTauMinus_IPFP","Compute_IPP(PV_x,PV_y,PV_z, diTauMinus_x,diTauMinus_y,diTauMinus_z, diTauMinus_x+diTauMinus_px,diTauMinus_y+diTauMinus_py,diTauMinus_z+diTauMinus_pz)")\
             .Define("diTauMinus_IPFx","diTauMinus_IPFP[0]-PV_x")\
             .Define("diTauMinus_IPFy","diTauMinus_IPFP[1]-PV_y")\
             .Define("diTauMinus_IPFz","diTauMinus_IPFP[2]-PV_z")\
             .Define("diTauMinus_IPF","Get_Norm(diTauMinus_IPFx,diTauMinus_IPFy,diTauMinus_IPFz)")

    #Get the Bs IP (vertex)
    rdf = rdf.Define("Bs_IPVP","Compute_IPP(Bs2TauTau_x,Bs2TauTau_y,Bs2TauTau_z,PV_x,PV_y,PV_z+1.0,PV_x,PV_y,PV_z)")\
             .Define("Bs_IPVx","Bs_IPVP[0]-Bs2TauTau_x")\
             .Define("Bs_IPVy","Bs_IPVP[1]-Bs2TauTau_y")\
             .Define("Bs_IPVz","Bs_IPVP[2]-Bs2TauTau_z")\
             .Define("Bs_IPV","Get_Norm(Bs_IPVx,Bs_IPVy,Bs_IPVz)")

    #Bs properties from RECO technique, Only visible momenta
    rdf = rdf.Define("Bs_FlightDistance","Get_Norm(Bs2TauTau_x-PV_x,Bs2TauTau_y-PV_y,Bs2TauTau_z-PV_z)")\
             .Define("Bs_px","diTauMinus_px+diTauPlus_px")\
             .Define("Bs_py","diTauMinus_py+diTauPlus_py")\
             .Define("Bs_pz","diTauMinus_pz+diTauPlus_pz")\
             .Define("Bs_p","Get_Norm(Bs_px,Bs_py,Bs_pz)")\
             .Define("Bs_Lifetime","5.3669/Bs_p*Bs_FlightDistance/299792458e3")

    #Get the various mass variable, but not the collinear mass on
    rdf = rdf.Define("DeltaM","mDiTau_Vis-diTauPlus_mass-diTauMinus_mass")


    return rdf

    

def Load_RDF(mode,amount):

    Path = "/eos/experiment/fcc/ee/analyses_storage/flavor/Bs2TauTau/flatNtuples/winter2023/analysis_stage1_withSimpleCut/"

    if mode == "sig": 
        
        #Get only one with 250k events, reduce it afterward in pandas
        NF = 1
        filenames = r.std.vector('string')()
        filenames.push_back(Path+Links[mode]+f"/chunk_0.root")
        rdf = r.RDataFrame("events",filenames)
        rdf = RDF_Treatment(rdf)
        return rdf
    
    else: 
        
        NF = 100
        filenames = r.std.vector('string')()
        for chunk in np.arange(0,NF,1):
            if (mode == "bb" and chunk == 51) or (mode == "cc" and chunk == 6) or (mode == "ud" and chunk == 97): continue #Remove some empty files
            filenames.push_back(Path+Links[mode]+f"/chunk_{chunk}.root")
        rdf = r.RDataFrame("events",filenames)
        rdf = RDF_Treatment(rdf)
        return rdf


def Do_RDF_PreTreatment(amount):

    VarList = ["diTau_Angles",
               "diTauPlus_IP","diTauPlus_IPV",
               "diTauMinus_IP","diTauMinus_IPV",
               
               "diTauPlus_FlightDistance","diTauMinus_FlightDistance",
               "diTauPlus_Lifetime","diTauMinus_Lifetime",
               
               "Bs_Lifetime","Bs_FlightDistance",
               "Bs_IPV",
               "mDiTau_Vis","DeltaM",
               ]

    print(f"Variables selected for training the BDT:\n{VarList}")

    columns = r.std.vector('string')()
    for var in VarList:
        columns.push_back(var)

    RDFs = {}
    for mode in Links.keys():
        RDFs[mode] = Load_RDF(mode,amount)
        RDFs[mode].Snapshot("events",f"/afs/cern.ch/work/t/tomonnar/public/Bs2TauTau/Stage2/BDT/PreTreatedFiles/Baseline_NoBug/{mode}.root",columns)

#=========================================================================

Do_RDF_PreTreatment(1)

