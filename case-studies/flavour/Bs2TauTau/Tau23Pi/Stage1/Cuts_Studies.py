import ROOT as r
#from termcolor import colored
import json

r.gInterpreter.Declare('''
    ROOT::VecOps::RVec<float> Cut_OnLists_float(ROOT::VecOps::RVec<float> List, float Threshold1, float Threshold2){
        ROOT::VecOps::RVec<float> newList;
        for (size_t i=0; i<List.size(); ++i){
            if (List.at(i) > Threshold1 && List.at(i) < Threshold2) {
                newList.push_back(List.at(i));
            }
        }
        return newList;
    }
''')

r.gInterpreter.Declare('''
    bool Cut_Tau_Mass(ROOT::VecOps::RVec<float> List, float Mmin, float Mmax, int Nwanted){
        int counter = 0;
        for (size_t i=0; i<List.size(); ++i){
            if (List.at(i) > Mmin && List.at(i) < Mmax){
                counter += 1;
            }
        }
        if (counter >= Nwanted){
            return true;
        }
        else{
            return false;
        }
    }
''')

r.gInterpreter.Declare('''
    float Compute_CosTheta(float axis_x, float axis_y, float axis_z, float x, float y, float z){
        float num = axis_x*x + axis_y*y + axis_z*z;
        float den = sqrt(axis_x*axis_x+axis_y*axis_y+axis_z*axis_z)*sqrt(x*x+y*y+z*z);
        return num/den;
    }
''')


r.gInterpreter.Declare('''
    bool CloseEnough(ROOT::VecOps::RVec<float> in, float cut_low, float cut_high){
        int Counter = 0;
        for (size_t i=0; i < in.size(); ++i){
            if (in[i] > cut_low and in[i] < cut_high) Counter += 1;
        }
        if (Counter == 0) return false;
        else return true;
    }
''')

r.gInterpreter.Declare('''
    bool OppositeSign_CloseEnough(ROOT::VecOps::RVec<float> angles, float threshold, ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> ID_Angles, ROOT::VecOps::RVec<int> Charges){
        bool result(false);
        for (size_t i = 0; i < angles.size(); ++i){
            if (angles[i] > threshold){
                if (Charges[ID_Angles[i][0]] + Charges[ID_Angles[i][1]] == 0){
                    result = true;
                }
            }
        }
        return result;
    }
''')

r.gInterpreter.Declare('''
    TLorentzVector build_p4(float px, float py, float pz, float mass) {
        TLorentzVector p4;
        p4.SetXYZM(px, py, pz, mass);
        return p4;
    }
''') 


def Build_FilterArg(CutList):
    
    FullList = ""
    
    if len(CutList) == 1:
        FullList = CutList[0]
    
    else:
        for Cut in CutList:
            if Cut == CutList[0]:
                FullList = Cut
            else:
                FullList += f" && {Cut}"

    return FullList


#-----------------------------------------------------------------------------------------------------------------------------------

################################# Warning: we may be reaching the statistical limit for conducting these studies ##################################################################

pathtofile = "/eos/experiment/fcc/ee/analyses_storage/flavor/Bs2TauTau/flatNtuples/winter2023/analysis_stage1_noFilter_fullP4/"

#Modes to be scanned
modes = ["p8_ee_Zbb_ecm91",
         "p8_ee_Zcc_ecm91",
         "p8_ee_Zss_ecm91",
         "p8_ee_Zud_ecm91",
         "p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU",
         "p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau"]

#Storage Values
Ntot = {}
NCut = {}

NTot = 0
NCuts = 0

#The Cut list to be studied
CutList = ["EVT_ThrustEmin_NTau23PiCand > 1", 
           "recoEmiss_e > 3.5",
           #"EVT_E < 87.0",
           #"EVT_Eneutral < 30.0",
           "EVT_ThrustEmin_Eneutral < 10",
           "EVT_ThrustEmin_Nneutral <= 13",
           #"EVT_ThrustEmin_p < 38.",
           "EVT_ThrustEmin_E < 38.",
           #"EVT_MVA1 > 0.9"
           #"EVT_ThrustEmin_NDV < 4",
           #"EVT_ThrustEmin_N <= 28",
           #"EVT_ThrustEmin_Echarged > 21."
           #"Cut_Tau_Mass(Tau23PiCandidates_mass,0.6,1.6,1)",
           #"Tau23PiCandidates_q.size() == 2", #Lost a lot of signal (25%->42%) while not much has changed for background
           #"Tau23PiCandidates_q.at(0) + Tau23PiCandidates_q.at(1) == 0",
           #"Bs_Mass > 2.3",
           #"Bs_Mass < 5.",
           #"Delta_M > 0.5",
           #"(EVT_ThrustEmin_NTau23PiCand == 2 || EVT_ThrustEmax_NTau23PiCand == 2)",
           #"EVT_ThrustEmax_Eneutral < 10",
           #"EVT_ThrustEmax_Nneutral <= 13",
           #"EVT_ThrustEmax_E < 38."
           #"EVT_NTau23Pi > 1",
           #"EVT_NVertex > 1",
           #"EVT_dPV2DVmax > 1",
           #"CloseEnough(Tau23PiCandidates_vertex_angles,0.95)",
           #"OppositeSign_CloseEnough(Tau23PiCandidates_p_angles,0.9,Tau23PiCandidates_vertex_angID,Tau23PiCandidates_q)",
           #"CloseEnough(Tau23PiCandidates_p_angles,0.9,0.99)"
           "diTau_p4DotProd > 2.5"
          ] 
           
           #Tau23PiCandidates_rho1mass.at(0) < 0.9
           #"EVT_ThrustEmin_Eneutral < 10",
           #"EVT_ThrustEmin_Nneutral <= 12", #not much gain since Eneutral is strongly correlated with Nneutral
           #"Tau23PiCandidates_rho1mass.size() > 0",
           #"Tau23PiCandidates_rho2mass.size() > 0",
           #"Tau23PiCandidates_q.size() == 2", #Lost a lot of signal (25%->42%) while not much has changed for background
           #"Tau23PiCandidates_q.at(0) + Tau23PiCandidates_q.at(1) == 0"


#prevents adding same lines in case calling multiple time the scipt 
Previous = {}
AlreadyDone = False
with open("Cuts_Studied.json","r") as ofile:
    Previous = json.load(ofile)
    
    for Cuts in Previous.keys():
        if Build_FilterArg(CutList) == Cuts:
            AlreadyDone = True
            print("Already studied")
            exit()

print("\n== Cuts studied ==\n")
for cut in CutList:
    print(f"{cut}")
print("\n")

for mode in modes:

    print(f"_____________ {mode} _______________")

    #Load the number of cuts before any cuts
    Ntot[mode] = Previous["Before"][mode]
    NTot += Ntot[mode]

    #load all the chunks
    c = r.TChain("Chunks")
    for chunk in [0,1,2,3,4,5,6,7,8,9]:
        c.Add(pathtofile+mode+f"/chunk_{chunk}.root/events")
    rdf = r.RDataFrame(c)

    #rdf = rdf.Define("EVT_ThrustEmin_p","sqrt(pow(EVT_ThrustEmin_px,2)+pow(EVT_ThrustEmin_py,2)+pow(EVT_ThrustEmin_pz,2))")\
    #         .Define("EVT_ThrustEmax_p","sqrt(pow(EVT_ThrustEmax_px,2)+pow(EVT_ThrustEmax_py,2)+pow(EVT_ThrustEmax_pz,2))")

    #Find close Tau23PiCandidates vertex (The CloseEnough function has a tunable threshold)
    rdf = rdf.Define("Tau23PiCandidates_vertex_x","ROOT::VecOps::RVec<float> result; for (size_t i=0; i<Tau23PiCandidates_vertex.size(); ++i) result.push_back(Vertex_x[Tau23PiCandidates_vertex[i]]); return result;")\
             .Define("Tau23PiCandidates_vertex_y","ROOT::VecOps::RVec<float> result; for (size_t i=0; i<Tau23PiCandidates_vertex.size(); ++i) result.push_back(Vertex_y[Tau23PiCandidates_vertex[i]]); return result;")\
             .Define("Tau23PiCandidates_vertex_z","ROOT::VecOps::RVec<float> result; for (size_t i=0; i<Tau23PiCandidates_vertex.size(); ++i) result.push_back(Vertex_z[Tau23PiCandidates_vertex[i]]); return result;")\
             .Define("Tau23PiCandidates_vertex_angles","ROOT::VecOps::RVec<float> result; for (size_t i=0;i<EVT_NTau23Pi;++i) {for (size_t j=i+1;j<EVT_NTau23Pi;++j) result.push_back(Compute_CosTheta(Tau23PiCandidates_vertex_x[i],Tau23PiCandidates_vertex_y[i],Tau23PiCandidates_vertex_z[i],Tau23PiCandidates_vertex_x[j],Tau23PiCandidates_vertex_y[j],Tau23PiCandidates_vertex_z[j]));} return result;")\
             \
             .Define("Tau23PiCandidates_p_angles","ROOT::VecOps::RVec<float> result; for (size_t i=0;i<EVT_NTau23Pi;++i) {for (size_t j=i+1;j<EVT_NTau23Pi;++j) result.push_back(Compute_CosTheta(Tau23PiCandidates_px[i],Tau23PiCandidates_py[i],Tau23PiCandidates_pz[i],Tau23PiCandidates_px[j],Tau23PiCandidates_py[j],Tau23PiCandidates_pz[j]));} return result;")\
             \
             .Define("EVT_E","EVT_ThrustEmin_E + EVT_ThrustEmax_E")\
             .Define("EVT_Echarged","EVT_ThrustEmin_Echarged + EVT_ThrustEmax_Echarged")\
             .Define("EVT_Eneutral","EVT_ThrustEmin_Eneutral + EVT_ThrustEmax_Eneutral")\
             .Define("EVT_N","EVT_ThrustEmin_N + EVT_ThrustEmax_N")\
             .Define("EVT_Ncharged","EVT_ThrustEmin_Ncharged + EVT_ThrustEmax_Ncharged")\
             .Define("EVT_Nneutral","EVT_ThrustEmin_Nneutral + EVT_ThrustEmax_Nneutral")\
             \
             .Define("TauCandPlus_x","float result; if (TauCand1_q > 0) return TauCand1_x; else return TauCand2_x;")\
             .Define("TauCandPlus_y","float result; if (TauCand1_q > 0) return TauCand1_y; else return TauCand2_y;")\
             .Define("TauCandPlus_z","float result; if (TauCand1_q > 0) return TauCand1_z; else return TauCand2_z;")\
             .Define("TauCandPlus_px","float result; if (TauCand1_q > 0) return TauCand1_px; else return TauCand2_px;")\
             .Define("TauCandPlus_py","float result; if (TauCand1_q > 0) return TauCand1_py; else return TauCand2_py;")\
             .Define("TauCandPlus_pz","float result; if (TauCand1_q > 0) return TauCand1_pz; else return TauCand2_pz;")\
             .Define("TauCandPlus_p","float result; if (TauCand1_q > 0) return TauCand1_p; else return TauCand2_p;")\
             .Define("TauCandPlus_m","float result; if (TauCand1_q > 0) return TauCand1_m; else return TauCand2_m;")\
             .Define("TauCandMinus_x","float result; if (TauCand1_q < 0) return TauCand1_x; else return TauCand2_x;")\
             .Define("TauCandMinus_y","float result; if (TauCand1_q < 0) return TauCand1_y; else return TauCand2_y;")\
             .Define("TauCandMinus_z","float result; if (TauCand1_q < 0) return TauCand1_z; else return TauCand2_z;")\
             .Define("TauCandMinus_px","float result; if (TauCand1_q < 0) return TauCand1_px; else return TauCand2_px;")\
             .Define("TauCandMinus_py","float result; if (TauCand1_q < 0) return TauCand1_py; else return TauCand2_py;")\
             .Define("TauCandMinus_pz","float result; if (TauCand1_q < 0) return TauCand1_pz; else return TauCand2_pz;")\
             .Define("TauCandMinus_p","float result; if (TauCand1_q < 0) return TauCand1_p; else return TauCand2_p;")\
             .Define("TauCandMinus_m","float result; if (TauCand1_q < 0) return TauCand1_m; else return TauCand2_m;")\
             \
             .Define("TauCandPlus_p4","build_p4(TauCandPlus_px,TauCandPlus_py,TauCandPlus_pz,TauCandPlus_m)")\
             .Define("TauCandMinus_p4","build_p4(TauCandMinus_px,TauCandMinus_py,TauCandMinus_pz,TauCandMinus_m)")\
             .Define("diTau_p4DotProd","TauCandPlus_p4.Dot(TauCandMinus_p4)")
             #.Define("Tau23PiCandidates_vertex_angID", "ID_Angles(Tau23PiCandidates_vertex_angles,EVT_NTau23Pi)")\


    #flatten the lists, cut and count the survivors
    #rdf = rdf.Redefine("Tau23PiCandidates_rho1mass","Cut_OnLists_float(Tau23PiCandidates_rho1mass,0.55,1.0)")
    #rdf = rdf.Redefine("Tau23PiCandidates_rho2mass","Cut_OnLists_float(Tau23PiCandidates_rho2mass,0.55,1.0)")
    #rdf = rdf.Define("EVT_ThrustEmin_p","sqrt(pow(EVT_ThrustEmin_px,2.0) + pow(EVT_ThrustEmin_py,2.0) + pow(EVT_ThrustEmin_pz,2.0))")
    
    #rdf = rdf.Filter("Tau23PiCandidates_q.size() == 2 && Tau23PiCandidates_q.at(0) + Tau23PiCandidates_q.at(1) == 0")
    #rdf = rdf.Define("Tau23PiCandidates_firstTau_pion1p",  "Tau23PiCandidates_pion1p.at(0)")\
    #         .Define("Tau23PiCandidates_firstTau_pion2p",  "Tau23PiCandidates_pion2p.at(0)")\
    #         .Define("Tau23PiCandidates_firstTau_pion3p",  "Tau23PiCandidates_pion3p.at(0)")\
    #         \
    #         .Define("Tau23PiCandidates_secondTau_pion1p",  "Tau23PiCandidates_pion1p.at(1)")\
    #         .Define("Tau23PiCandidates_secondTau_pion2p",  "Tau23PiCandidates_pion2p.at(1)")\
    #         .Define("Tau23PiCandidates_secondTau_pion3p",  "Tau23PiCandidates_pion3p.at(1)")\
    #         \
    #         .Define("Tau23PiCandidates_firstTau_pion1E",  "sqrt(pow(Tau23PiCandidates_firstTau_pion1p,2) + pow(0.13957,2))")\
    #         .Define("Tau23PiCandidates_firstTau_pion2E",  "sqrt(pow(Tau23PiCandidates_firstTau_pion2p,2) + pow(0.13957,2))")\
    #         .Define("Tau23PiCandidates_firstTau_pion3E",  "sqrt(pow(Tau23PiCandidates_firstTau_pion3p,2) + pow(0.13957,2))")\
    #         \
    #         .Define("Tau23PiCandidates_secondTau_pion1E",  "sqrt(pow(Tau23PiCandidates_secondTau_pion1p,2) + pow(0.13957,2))")\
    #         .Define("Tau23PiCandidates_secondTau_pion2E",  "sqrt(pow(Tau23PiCandidates_secondTau_pion2p,2) + pow(0.13957,2))")\
    #         .Define("Tau23PiCandidates_secondTau_pion3E",  "sqrt(pow(Tau23PiCandidates_secondTau_pion3p,2) + pow(0.13957,2))")\
    #         \
    #         .Define("Tau23PiCandidates_firstTau_pion1px",  "Tau23PiCandidates_pion1px.at(0)")\
    #         .Define("Tau23PiCandidates_firstTau_pion1py",  "Tau23PiCandidates_pion1py.at(0)")\
    #         .Define("Tau23PiCandidates_firstTau_pion1pz",  "Tau23PiCandidates_pion1pz.at(0)")\
    #         .Define("Tau23PiCandidates_firstTau_pion2px",  "Tau23PiCandidates_pion2px.at(0)")\
    #         .Define("Tau23PiCandidates_firstTau_pion2py",  "Tau23PiCandidates_pion2py.at(0)")\
    #         .Define("Tau23PiCandidates_firstTau_pion2pz",  "Tau23PiCandidates_pion2pz.at(0)")\
    #         .Define("Tau23PiCandidates_firstTau_pion3px",  "Tau23PiCandidates_pion3px.at(0)")\
    #         .Define("Tau23PiCandidates_firstTau_pion3py",  "Tau23PiCandidates_pion3py.at(0)")\
    #         .Define("Tau23PiCandidates_firstTau_pion3pz",  "Tau23PiCandidates_pion3pz.at(0)")\
    #         \
    #         .Define("Tau23PiCandidates_secondTau_pion1px",  "Tau23PiCandidates_pion1px.at(1)")\
    #         .Define("Tau23PiCandidates_secondTau_pion1py",  "Tau23PiCandidates_pion1py.at(1)")\
    #         .Define("Tau23PiCandidates_secondTau_pion1pz",  "Tau23PiCandidates_pion1pz.at(1)")\
    #         .Define("Tau23PiCandidates_secondTau_pion2px",  "Tau23PiCandidates_pion2px.at(1)")\
    #         .Define("Tau23PiCandidates_secondTau_pion2py",  "Tau23PiCandidates_pion2py.at(1)")\
    #         .Define("Tau23PiCandidates_secondTau_pion2pz",  "Tau23PiCandidates_pion2pz.at(1)")\
    #         .Define("Tau23PiCandidates_secondTau_pion3px",  "Tau23PiCandidates_pion3px.at(1)")\
    #         .Define("Tau23PiCandidates_secondTau_pion3py",  "Tau23PiCandidates_pion3py.at(1)")\
    #         .Define("Tau23PiCandidates_secondTau_pion3pz",  "Tau23PiCandidates_pion3pz.at(1)")\
    #         \
    #         .Define("Bs_Mass","sqrt(pow(Tau23PiCandidates_firstTau_pion1E+Tau23PiCandidates_firstTau_pion2E+Tau23PiCandidates_firstTau_pion3E+Tau23PiCandidates_secondTau_pion1E+Tau23PiCandidates_secondTau_pion2E+Tau23PiCandidates_secondTau_pion3E,2) - "+\
    #                                "pow(Tau23PiCandidates_firstTau_pion1px+Tau23PiCandidates_firstTau_pion2px+Tau23PiCandidates_firstTau_pion3px+Tau23PiCandidates_secondTau_pion1px+Tau23PiCandidates_secondTau_pion2px+Tau23PiCandidates_secondTau_pion3px,2) - "+\
    #                                "pow(Tau23PiCandidates_firstTau_pion1py+Tau23PiCandidates_firstTau_pion2py+Tau23PiCandidates_firstTau_pion3py+Tau23PiCandidates_secondTau_pion1py+Tau23PiCandidates_secondTau_pion2py+Tau23PiCandidates_secondTau_pion3py,2) - "+\
    #                                "pow(Tau23PiCandidates_firstTau_pion1pz+Tau23PiCandidates_firstTau_pion2pz+Tau23PiCandidates_firstTau_pion3pz+Tau23PiCandidates_secondTau_pion1pz+Tau23PiCandidates_secondTau_pion2pz+Tau23PiCandidates_secondTau_pion3pz,2))")\
    #         \
    #         .Define("Delta_M","Bs_Mass - Tau23PiCandidates_mass.at(0) - Tau23PiCandidates_mass.at(1)")\
    
    rdf = rdf.Filter(Build_FilterArg(CutList))
    NCut[mode] = rdf.Count().GetValue()
    NCuts += NCut[mode]

    if mode == "p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU" or mode == "p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau":
        color = "blue"
    else:
        color = "red"
    
    #Efficiencies of how many got rejected so red should be close to 100% and blue close to 0%
    print(f"Cuts eff = {round(100*(Ntot[mode]-NCut[mode])/Ntot[mode],5)} %")

print("\n-----------------------------------------")
print(f"=====> Total eff = {round(100*(NTot-NCuts)/NTot,5)} %")
print("-----------------------------------------\n")


#Store these numbers
with open("Cuts_Studied.json","w") as ofile:
    Previous[Build_FilterArg(CutList)] = NCut
    json.dump(Previous, ofile)

     
            
