import ROOT as r
import numpy as np
import sys

#sys.path.insert(1,"../Utils")
#import Functions_ForPreTreat

r.gInterpreter.Declare('''

    ROOT::VecOps::RVec<int> Find_Vertex2Daughters(ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>>recoPDG, ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>>recoq, int PDG1, int PDG2){
    
        ROOT::VecOps::RVec<int> results;

        for (size_t i=0; i<recoPDG.size(); ++i){
            if (recoPDG.at(i).size() == 2){ //exactly two tracks
                int AbsPDG0 (std::abs(recoPDG.at(i).at(0)));
                int AbsPDG1 (std::abs(recoPDG.at(i).at(1)));

                int whoPDG1 (-1);
                int whoPDG2 (-1);
                
                //Check the id of both part separately
                if (AbsPDG0 == PDG1){
                    whoPDG1 = 0;
                    if (AbsPDG1 == PDG2) whoPDG2 = 1;
                }
                else if (AbsPDG1 == PDG1){
                    whoPDG1 = 1;
                    if (AbsPDG0 == PDG2) whoPDG2 = 0;
                }

                if(whoPDG1 + whoPDG2 == 1){
                    //The two expected particles are present
                    if (recoq.at(i).at(0)*recoq.at(i).at(1) < 0){ 
                        //They are of opposite charge (works with pion and protons, tbc with others)
                        results.push_back(i);
                    }
                }
            }
        }

        return results;
    }

    //Select the candidate of interest
    int Select_Vertex2Daughters_ind(ROOT::VecOps::RVec<int> Vertex_ind, ROOT::VecOps::RVec<float> vertex_x, ROOT::VecOps::RVec<float> vertex_y){
        if (Vertex_ind.size() == 0) return -1;
        else if (Vertex_ind.size() == 1) return Vertex_ind.at(0);
        
        //Select smallest radius position if multiple vertices
        else {
            float rmin ( sqrt( vertex_x.at(0)*vertex_x.at(0) + vertex_y.at(0)*vertex_y.at(0) ) );
            int min_ind (0);
            for (size_t i=1;i<Vertex_ind.size();++i){
                float ri ( sqrt( vertex_x.at(i)*vertex_x.at(i) + vertex_y.at(i)*vertex_y.at(i) ) );
                if ( ri < rmin ) rmin = ri;
                min_ind = i;
            }
            return min_ind;
        }
    
    }

    
    //-----------------------------------------------------------------------------------------------------------------
    
    
    
    float Find_Vertex2Daughters_MassDistrib(ROOT::VecOps::RVec<float> recoMass, ROOT::VecOps::RVec<int> vertexcand_ind){

        //Treat only the events with one Id strange vertices (no ambiguity)
        if (vertexcand_ind.size() == 1) return recoMass.at(vertexcand_ind.at(0));
        else return -999.9;
    }

    float Compute_Vertex2Daughters_r(ROOT::VecOps::RVec<float> x, ROOT::VecOps::RVec<float> y, ROOT::VecOps::RVec<int> vertexcand_ind){
    
        //Treat only the events with one Id strange vertices (no ambiguity)
        if (vertexcand_ind.size() == 1) return sqrt(x.at(vertexcand_ind.at(0))*x.at(vertexcand_ind.at(0)) + y.at(vertexcand_ind.at(0))*y.at(vertexcand_ind.at(0)));
        else return -999.9;
    }

    float Compute_Vertex2Daughters_d2PV(ROOT::VecOps::RVec<float> x, 
                                        ROOT::VecOps::RVec<float> y, 
                                        ROOT::VecOps::RVec<float> z, 
                                        ROOT::VecOps::RVec<int> vertexcand_ind,
                                        ROOT::VecOps::RVec<float> x_all,
                                        ROOT::VecOps::RVec<float> y_all,
                                        ROOT::VecOps::RVec<float> z_all,
                                        ROOT::VecOps::RVec<int> isPV){
    
        //Treat only the events with one Id strange vertices (no ambiguity)
        if (vertexcand_ind.size() == 1){
            float xp (0.0);
            float yp (0.0);
            float zp (0.0);
            bool FoundPV (false);
            for (size_t i=0; i<isPV.size(); ++i){
                if (isPV.at(i) == 1){
                    xp = x_all.at(i);
                    yp = y_all.at(i);
                    zp = z_all.at(i);
                    FoundPV = true;
                }
            }
            if (not FoundPV) return -999.9;
            float xv (x.at(vertexcand_ind.at(0)));
            float yv (y.at(vertexcand_ind.at(0)));
            float zv (z.at(vertexcand_ind.at(0)));
            return sqrt( (xv-xp)*(xv-xp) + (yv-yp)*(yv-yp) + (zv-zp)*(zv-zp) );
        }
        else return -999.9;
    }

    float Find_Vertex2Daughters_PV_DOCA(float x, ROOT::VecOps::RVec<float> px, float y, ROOT::VecOps::RVec<float> py, float z, ROOT::VecOps::RVec<float> pz, float pvx, float pvy, float pvz){

        //Get the total vertex momentum
        float pxtot (0.0);
        float pytot (0.0);
        float pztot (0.0);
        for (size_t i=0;i<px.size();++i){
            pxtot += px.at(i);
            pytot += py.at(i);
            pztot += pz.at(i);
        }

        //Compute the DOCA of (Vertex-Momentum,PV)
        float VPVx (x-pvx);
        float VPVy (y-pvy);
        float VPVz (z-pvz);

        float VPVxMx (VPVy*pztot-VPVz*pytot);
        float VPVxMy (VPVz*pxtot-VPVx*pztot);
        float VPVxMz (VPVx*pytot-VPVy*pxtot);

        return sqrt(VPVxMx*VPVxMx + VPVxMy*VPVxMy + VPVxMz*VPVxMz)/sqrt(pxtot*pxtot + pytot*pytot + pztot*pztot);

    }

''')

Links = {"bb": "p8_ee_Zbb_ecm91",
         "cc": "p8_ee_Zcc_ecm91",
         "ss": "p8_ee_Zss_ecm91",
         "ud": "p8_ee_Zud_ecm91",
         "sig":"p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau"
        }

def RDF_Treatment(rdf):

    for var in ["px","py","pz","phi","eta","energy","mass","charge","PDG","thrustangles"]:
        rdf = rdf.Define(f"plus_{var}",f"lepton_{var}.at(dilepton_plus_ind)").Define(f"minus_{var}",f"lepton_{var}.at(dilepton_minus_ind)")
    rdf = rdf.Define("Opening_Angle","(plus_px*minus_px+plus_py*minus_py+plus_pz*minus_pz)/sqrt(plus_px*plus_px+plus_py*plus_py+plus_pz*plus_pz)/sqrt(minus_px*minus_px+minus_py*minus_py+minus_pz*minus_pz)")
    rdf = rdf.Define("EAsymm","(EVT_ThrustEmax_E-EVT_ThrustEmin_E)/(EVT_ThrustEmin_E+EVT_ThrustEmax_E)")

    #------ The strange hadron vertices info --------------------

    rdf = rdf.Define("Vertex_sighemi_2pi_ind","Find_Vertex2Daughters(Vertex_sighemi_RECO_PDG,Vertex_sighemi_RECO_charge,211,211)")
    rdf = rdf.Define("Vertex_sighemi_ppi_ind","Find_Vertex2Daughters(Vertex_sighemi_RECO_PDG,Vertex_sighemi_RECO_charge,211,2212)")
    rdf = rdf.Define("nVertex_sighemi_2pi","Vertex_sighemi_2pi_ind.size()")
    rdf = rdf.Define("nVertex_sighemi_ppi","Vertex_sighemi_ppi_ind.size()")
    
    rdf = rdf.Define("TheVertex_sighemi_2pi","Select_Vertex2Daughters_ind(Vertex_sighemi_2pi_ind,Vertex_sighemi_x,Vertex_sighemi_y)")
    rdf = rdf.Define("TheVertex_sighemi_ppi","Select_Vertex2Daughters_ind(Vertex_sighemi_ppi_ind,Vertex_sighemi_x,Vertex_sighemi_y)")
    
    rdf = rdf.Define("Vertex_sighemi_2pi_mass","Find_Vertex2Daughters_MassDistrib(Vertex_sighemi_mass,Vertex_sighemi_2pi_ind)")
    rdf = rdf.Define("Vertex_sighemi_ppi_mass","Find_Vertex2Daughters_MassDistrib(Vertex_sighemi_mass,Vertex_sighemi_ppi_ind)")
    rdf = rdf.Define("TheVertex_sighemi_2pi_mass","if (TheVertex_sighemi_2pi>-1) return Vertex_sighemi_mass.at(TheVertex_sighemi_2pi); else return float(-999.9);")
    rdf = rdf.Define("TheVertex_sighemi_ppi_mass","if (TheVertex_sighemi_ppi>-1) return Vertex_sighemi_mass.at(TheVertex_sighemi_ppi); else return float(-999.9);")

    rdf = rdf.Define("Vertex_sighemi_2pi_r","Compute_Vertex2Daughters_r(Vertex_sighemi_x,Vertex_sighemi_y,Vertex_sighemi_2pi_ind)")
    rdf = rdf.Define("Vertex_sighemi_ppi_r","Compute_Vertex2Daughters_r(Vertex_sighemi_x,Vertex_sighemi_y,Vertex_sighemi_ppi_ind)")
    
    rdf = rdf.Define("TheVertex_sighemi_2pi_x","if (TheVertex_sighemi_2pi>-1) return Vertex_sighemi_x.at(TheVertex_sighemi_2pi); else return float(-999.9);")
    rdf = rdf.Define("TheVertex_sighemi_ppi_x","if (TheVertex_sighemi_ppi>-1) return Vertex_sighemi_x.at(TheVertex_sighemi_ppi); else return float(-999.9);")
    rdf = rdf.Define("TheVertex_sighemi_2pi_y","if (TheVertex_sighemi_2pi>-1) return Vertex_sighemi_y.at(TheVertex_sighemi_2pi); else return float(-999.9);")
    rdf = rdf.Define("TheVertex_sighemi_ppi_y","if (TheVertex_sighemi_ppi>-1) return Vertex_sighemi_y.at(TheVertex_sighemi_ppi); else return float(-999.9);")
    rdf = rdf.Define("TheVertex_sighemi_2pi_z","if (TheVertex_sighemi_2pi>-1) return Vertex_sighemi_z.at(TheVertex_sighemi_2pi); else return float(-999.9);")
    rdf = rdf.Define("TheVertex_sighemi_ppi_z","if (TheVertex_sighemi_ppi>-1) return Vertex_sighemi_z.at(TheVertex_sighemi_ppi); else return float(-999.9);")

    rdf = rdf.Define("TheVertex_sighemi_2pi_r","sqrt(TheVertex_sighemi_2pi_x*TheVertex_sighemi_2pi_x + TheVertex_sighemi_2pi_y*TheVertex_sighemi_2pi_y)")
    rdf = rdf.Define("TheVertex_sighemi_ppi_r","sqrt(TheVertex_sighemi_ppi_x*TheVertex_sighemi_ppi_x + TheVertex_sighemi_ppi_y*TheVertex_sighemi_ppi_y)")

    rdf = rdf.Define("PV_x_","Vertex_x [Vertex_isPV > 0]")
    rdf = rdf.Define("PV_x","PV_x_.at(0)")
    rdf = rdf.Define("PV_y_","Vertex_y [Vertex_isPV > 0]")
    rdf = rdf.Define("PV_y","PV_y_.at(0)")
    rdf = rdf.Define("PV_z_","Vertex_z [Vertex_isPV > 0]")
    rdf = rdf.Define("PV_z","PV_z_.at(0)")
    
    rdf = rdf.Define("TheVertex_sighemi_2pi_d2PV_x","TheVertex_sighemi_2pi_x-PV_x")
    rdf = rdf.Define("TheVertex_sighemi_ppi_d2PV_x","TheVertex_sighemi_ppi_x-PV_x")
    rdf = rdf.Define("TheVertex_sighemi_2pi_d2PV_y","TheVertex_sighemi_2pi_y-PV_y")
    rdf = rdf.Define("TheVertex_sighemi_ppi_d2PV_y","TheVertex_sighemi_ppi_y-PV_y")
    rdf = rdf.Define("TheVertex_sighemi_2pi_d2PV_z","TheVertex_sighemi_2pi_z-PV_z")
    rdf = rdf.Define("TheVertex_sighemi_ppi_d2PV_z","TheVertex_sighemi_ppi_z-PV_z")
    
    rdf = rdf.Define("TheVertex_sighemi_2pi_d2PV","sqrt(TheVertex_sighemi_2pi_d2PV_x*TheVertex_sighemi_2pi_d2PV_x + TheVertex_sighemi_2pi_d2PV_y*TheVertex_sighemi_2pi_d2PV_y + TheVertex_sighemi_2pi_d2PV_z*TheVertex_sighemi_2pi_d2PV_z)")
    rdf = rdf.Define("TheVertex_sighemi_ppi_d2PV","sqrt(TheVertex_sighemi_ppi_d2PV_x*TheVertex_sighemi_ppi_d2PV_x + TheVertex_sighemi_ppi_d2PV_y*TheVertex_sighemi_ppi_d2PV_y + TheVertex_sighemi_ppi_d2PV_z*TheVertex_sighemi_ppi_d2PV_z)")

    rdf = rdf.Define("TheVertex_sighemi_2pi_px","if (TheVertex_sighemi_2pi>-1) return Vertex_sighemi_RECO_px.at(TheVertex_sighemi_2pi); else return ROOT::VecOps::RVec<float>({-999.9,-999.9});")
    rdf = rdf.Define("TheVertex_sighemi_ppi_px","if (TheVertex_sighemi_ppi>-1) return Vertex_sighemi_RECO_px.at(TheVertex_sighemi_ppi); else return ROOT::VecOps::RVec<float>({-999.9,-999.9});")
    rdf = rdf.Define("TheVertex_sighemi_2pi_py","if (TheVertex_sighemi_2pi>-1) return Vertex_sighemi_RECO_py.at(TheVertex_sighemi_2pi); else return ROOT::VecOps::RVec<float>({-999.9,-999.9});")
    rdf = rdf.Define("TheVertex_sighemi_ppi_py","if (TheVertex_sighemi_ppi>-1) return Vertex_sighemi_RECO_py.at(TheVertex_sighemi_ppi); else return ROOT::VecOps::RVec<float>({-999.9,-999.9});")
    rdf = rdf.Define("TheVertex_sighemi_2pi_pz","if (TheVertex_sighemi_2pi>-1) return Vertex_sighemi_RECO_pz.at(TheVertex_sighemi_2pi); else return ROOT::VecOps::RVec<float>({-999.9,-999.9});")
    rdf = rdf.Define("TheVertex_sighemi_ppi_pz","if (TheVertex_sighemi_ppi>-1) return Vertex_sighemi_RECO_pz.at(TheVertex_sighemi_ppi); else return ROOT::VecOps::RVec<float>({-999.9,-999.9});")

    rdf = rdf.Define("TheVertex_sighemi_2pi_fromPV","Find_Vertex2Daughters_PV_DOCA(TheVertex_sighemi_2pi_x,TheVertex_sighemi_2pi_px,TheVertex_sighemi_2pi_y,TheVertex_sighemi_2pi_py,TheVertex_sighemi_2pi_z,TheVertex_sighemi_2pi_pz,PV_x,PV_y,PV_z)")
    rdf = rdf.Define("TheVertex_sighemi_ppi_fromPV","Find_Vertex2Daughters_PV_DOCA(TheVertex_sighemi_ppi_x,TheVertex_sighemi_ppi_px,TheVertex_sighemi_ppi_y,TheVertex_sighemi_ppi_py,TheVertex_sighemi_ppi_z,TheVertex_sighemi_ppi_pz,PV_x,PV_y,PV_z)")

    return rdf

    

def Load_RDF(mode,amount):

    Path = "/eos/experiment/fcc/ee/analyses_storage/flavor/Bs2TauTau/flatNtuples/winter2023/analysis_stage1_Leptons_withCuts_April26/"

    if mode == "sig": 
        
        NF = 15
        filenames = r.std.vector('string')()
        for chunk in np.arange(0,NF,1):
            filenames.push_back(Path+Links[mode]+f"/chunk_{chunk}.root")
        rdf = r.RDataFrame("events",filenames)
        rdf = RDF_Treatment(rdf)
        return rdf
    
    elif mode == "bb":

        NF = 19
        filenames = r.std.vector('string')()
        for chunk in np.arange(0,NF,1):
            filenames.push_back(Path+Links[mode]+f"/chunk_{chunk}.root")
        rdf = r.RDataFrame("events",filenames)
        rdf = RDF_Treatment(rdf)
        return rdf

    else: 
        
        NF = 100
        filenames = r.std.vector('string')()
        for chunk in np.arange(0,NF,1):
            filenames.push_back(Path+Links[mode]+f"/chunk_{chunk}.root")
        rdf = r.RDataFrame("events",filenames)
        rdf = RDF_Treatment(rdf)
        return rdf


def Do_RDF_PreTreatment(amount):

    VarList = [#Muon related
               "plus_px","plus_py","plus_pz","plus_phi","plus_eta","plus_energy","plus_mass","plus_thrustangles",
               "minus_px","minus_py","minus_pz","minus_phi","minus_eta","minus_energy","minus_mass","minus_thrustangles","Opening_Angle","dilepton_case",

               #Event level
               "EVT_ThrustEmax_E","EVT_ThrustEmin_E","EVT_ThrustEmax_Echarged","EVT_ThrustEmin_Echarged","EVT_ThrustEmax_Eneutral","EVT_ThrustEmin_Eneutral",
               "EVT_ThrustEmax_N","EVT_ThrustEmin_N","EVT_ThrustEmax_Ncharged","EVT_ThrustEmin_Ncharged","EVT_ThrustEmax_Nneutral","EVT_ThrustEmin_Nneutral",
               "recoEmiss_thrustangle","recoEmiss_e","EVT_Thrust_Mag","EVT_Thrust_X","EVT_Thrust_Y","EVT_Thrust_Z","EAsymm",

               #Vertex related
               "EVT_ThrustEmin_NDV","EVT_ThrustEmax_NDV","EVT_dPV2DVmin","EVT_dPV2DVmax","EVT_dPV2DVave",
               "EVT_NtracksPV","EVT_NVertex",

               #Strange vertices related
               "nVertex_sighemi_2pi","nVertex_sighemi_ppi",
               "TheVertex_sighemi_2pi_mass","TheVertex_sighemi_ppi_mass",
               "TheVertex_sighemi_2pi_r","TheVertex_sighemi_ppi_r","TheVertex_sighemi_2pi_d2PV","TheVertex_sighemi_ppi_d2PV",
               "TheVertex_sighemi_2pi_fromPV","TheVertex_sighemi_ppi_fromPV",
               ]

               #Same as EVT_NVertex
               #"Vertex_n",

               #Not simple columns
               #"Vertex_d2PV","Vertex_mass","Vertex_ntrk","Vertex_chi2",

    print(f"Variables selected for training the BDT:\n{VarList}")

    columns = r.std.vector('string')()
    for var in VarList:
        columns.push_back(var)

    RDFs = {}
    for mode in ["sig","bb","cc","ss","ud"]: #READD AFTER PRODUCITON
        RDFs[mode] = Load_RDF(mode,amount)
        RDFs[mode].Snapshot("events",f"{mode}/StrangeHadronBkg.root",columns)

#=========================================================================

Do_RDF_PreTreatment(1)

