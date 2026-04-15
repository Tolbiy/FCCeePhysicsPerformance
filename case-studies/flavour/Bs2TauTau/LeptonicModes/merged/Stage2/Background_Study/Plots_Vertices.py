import ROOT as r
import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amssymb}'

#Dict for shortened dict key(global var -> put it in a config file)
Links = {"bb": "p8_ee_Zbb_ecm91",
         "cc": "p8_ee_Zcc_ecm91",
         "ss": "p8_ee_Zss_ecm91",
         "ud": "p8_ee_Zud_ecm91",
         "sig":"p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau",
         "exc":"p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU"
        }

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

#Get the model to be added to the dataframe
#r.gInterpreter.ProcessLine('''
#TMVA::Experimental::RBDT<> bdt("Naive_withMoreData", "/afs/cern.ch/work/t/tomonnar/public/FCCeePhysicsPerformance/case-studies/flavour/Bs2TauTau/LeptonicModes/merged/Stage2/BDT/TrainTest/Train_Results/Models/Naive_withMoreData.root");
#computeModel1 = TMVA::Experimental::Compute<44, float>(bdt);
#''')

#Get the histograms in a numpy form -----------------------------------
def Load_AsNumpy(Mode,var,Nbins,withEdge):
        
    Path = "/eos/experiment/fcc/ee/analyses_storage/flavor/Bs2TauTau/flatNtuples/winter2023/analysis_stage1_Leptons_withCuts_bkgStudies/"

    #Load the rdf properly (avoid segfault from Tree going out of scope)
    filenames = r.std.vector('string')()
    if Mode == "sig":
        for i in np.arange(0,20,1):
            filenames.push_back(Path+Links[Mode]+f"/chunk_{i}.root")
    else:
        for i in np.arange(0,10,1):
            filenames.push_back(Path+Links[Mode]+f"/chunk_{i}.root")
    rdf = r.RDataFrame("events",filenames)
    print(f"{Mode}: Number of events = {rdf.Count().GetValue()}")
 
    #rdf = rdf.Filter("n_lepton > 1 && EVT_ThrustEmin_E < 38 && recoEmiss_e > 10 && EVT_ThrustEmin_Eneutral < 10 && has_dilepton > 0 && has_dilepton_SameSide > 0 && has_dilepton_SigHemi > 0 && has_dilepton_OppositeCharges > 0 && has_dilepton_Vertex == 0")
    #print(f"{Mode}: Number of selected events = {rdf.Count().GetValue()}")
    

    #for vari in ["px","py","pz","phi","eta","energy","mass","charge","PDG","thrustangles"]:
    #    rdf = rdf.Define(f"plus_{vari}",f"lepton_{vari}.at(dilepton_plus_ind)").Define(f"minus_{vari}",f"lepton_{vari}.at(dilepton_minus_ind)")
    #rdf = rdf.Define("Opening_Angle","(plus_px*minus_px+plus_py*minus_py+plus_pz*minus_pz)/sqrt(plus_px*plus_px+plus_py*plus_py+plus_pz*plus_pz)/sqrt(minus_px*minus_px+minus_py*minus_py+minus_pz*minus_pz)")
    #rdf = rdf.Define("EAsymm","(EVT_ThrustEmax_E-EVT_ThrustEmin_E)/(EVT_ThrustEmin_E+EVT_ThrustEmax_E)")
    #rdf = rdf.Redefine("dilepton_case","float(dilepton_case)")

    #Add BDT
    #BDTvars = [#Muon related
    #           "plus_px","plus_py","plus_pz","plus_phi","plus_eta","plus_energy","plus_mass","plus_thrustangles",
    #           "minus_px","minus_py","minus_pz","minus_phi","minus_eta","minus_energy","minus_mass","minus_thrustangles","Opening_Angle","dilepton_case",
#
    #           #Event level
    #           "EVT_ThrustEmax_E","EVT_ThrustEmin_E","EVT_ThrustEmax_Echarged","EVT_ThrustEmin_Echarged","EVT_ThrustEmax_Eneutral","EVT_ThrustEmin_Eneutral",
    #           "EVT_ThrustEmax_N","EVT_ThrustEmin_N","EVT_ThrustEmax_Ncharged","EVT_ThrustEmin_Ncharged","EVT_ThrustEmax_Nneutral","EVT_ThrustEmin_Nneutral",
    #           "recoEmiss_thrustangle","recoEmiss_e","EVT_Thrust_Mag","EVT_Thrust_X","EVT_Thrust_Y","EVT_Thrust_Z","EAsymm",
#
    #           #Vertex related
    #           "EVT_ThrustEmin_NDV","EVT_ThrustEmax_NDV","EVT_dPV2DVmin","EVT_dPV2DVmax","EVT_dPV2DVave",
    #           "EVT_NtracksPV","EVT_NVertex",
    #           ]

    #rdf = rdf.Define("MVAVec2", r.computeModel1,BDTvars).Define("MVA2_withMoreData", "MVAVec2.at(0)")

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
    
    #rdf.Display(["PV_x","TheVertex_sighemi_2pi_x","TheVertex_sighemi_ppi_x"],30).Print()

    #rdf = rdf.Define("Vertex_sighemi_2pi_d2PV","Compute_Vertex2Daughters_d2PV(Vertex_sighemi_x,Vertex_sighemi_y,Vertex_sighemi_z,Vertex_sighemi_2pi_ind,Vertex_x,Vertex_y,Vertex_z,Vertex_isPV)")
    #rdf = rdf.Define("Vertex_sighemi_ppi_d2PV","Compute_Vertex2Daughters_d2PV(Vertex_sighemi_x,Vertex_sighemi_y,Vertex_sighemi_z,Vertex_sighemi_ppi_ind,Vertex_x,Vertex_y,Vertex_z,Vertex_isPV)")

    #Select the smallest signal hemisphere vertex thrust angle

    ##### To Add Any RDF Treatment Needed #####

    ###########################################

    hist, edges = np.histogram(rdf.AsNumpy([var])[var],bins=Nbins[0],range=(Nbins[1],Nbins[2]))

    if withEdge:
        return hist, edges
    else:
        return hist

#Draw the resolution sig + bkg ------------------------------
def Plot(hlist,edges,var):

    #To normalise and compare
    Norm = {}
    for mode in hlist.keys():
        Norm[mode] = np.sum(hlist[mode])

    fig, ax = plt.subplots()
    ax.grid(color='grey', linestyle='--', linewidth=0.5, alpha=0.5)

    #ax.stairs(hlist["ud"]/Norm["bkg"] + hlist["ss"]/Norm["bkg"] + hlist["cc"]/Norm["bkg"] + hlist["bb"]/Norm["bkg"],  edges,lw=0,fill=True,color='slategrey',label=r"$Z^0\rightarrow u\overline{d}$")
    #ax.stairs(hlist["ss"]/Norm["bkg"] + hlist["cc"]/Norm["bkg"] + hlist["bb"]/Norm["bkg"],                            edges,lw=0,fill=True,color='mediumseagreen',label=r"$Z^0\rightarrow s\overline{s}$")
    #ax.stairs(hlist["cc"]/Norm["bkg"] + hlist["bb"]/Norm["bkg"],                                                      edges,lw=0,fill=True,color='goldenrod',label=r"$Z^0\rightarrow c\overline{c}$")
    ax.stairs(hlist["bb"]/Norm["bkg"],                                                                                edges,lw=0,fill=True,color='steelblue',label=r"$Z^0\rightarrow b\overline{b}$")
    
    ax.stairs(hlist["bkg"]/Norm["bkg"],edges,ec='black',ls='-',lw=2,label=r"$\textrm{Total background}$")
    #ax.stairs(hlist["exc"]/Norm["exc"],edges,ec='maroon',ls='--',lw=2,label=r"$\textrm{Signal (Exclusive)}$")
    ax.stairs(hlist["sig"]/Norm["sig"],edges,ec='red',ls='-',lw=2,label=r"$\textrm{Signal}$")

    #Set axis range
    ax.set_xlim([edges[0],edges[-1]])

    #Labels
    ax.set_xlabel(r"$\textrm{"+f"{var}"+r"}$",size="large")
    ax.set_ylabel(r"$\textrm{Normalised events}$",size="large")

    #Legend (Count the number of bkg modes drawn to write it down properly)
    ax.legend()
        
    #Save
    fig.savefig(f'{var}.pdf')

#=======================================================================================================

def Do_Plots(var,bins):

    for v in var:
        hmode = {}
        hmode["sig"], edges = Load_AsNumpy("sig",v,bins,True)
        for mode in ["bb",]: #"cc","ss","ud"
            hmode[mode] = Load_AsNumpy(mode,v,bins,False)
        hmode["bkg"] = hmode["bb"]#+hmode["cc"]+hmode["ss"]+hmode["ud"]
        Plot(hmode,edges,v)


#=======================================================================================================

parser = argparse.ArgumentParser()
parser.add_argument("NBins", type=int)
parser.add_argument("Low",   type=float)
parser.add_argument("High",  type=float)
parser.add_argument("List",   nargs='+')
args=parser.parse_args()
Do_Plots(args.List,(args.NBins,args.Low,args.High))