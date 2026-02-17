import ROOT as r


#Some analyzers needed in general (to be included in Functions.h when FCCAnalyses Stage 1)

#General Purpose--------------------------------------------------------------------------------------------

#Compute the distance between two points
r.gInterpreter.Declare('''
    ROOT::VecOps::RVec<float> Compute_CosTheta(float axis_x, float axis_y, float axis_z, ROOT::VecOps::RVec<float> x, ROOT::VecOps::RVec<float> y, ROOT::VecOps::RVec<float> z){
        ROOT::VecOps::RVec<float> results; 
        for(size_t i=0; i<x.size(); ++i) 
            {
                float num = axis_x*x[i] + axis_y*y[i] + axis_z*z[i];
                float den = sqrt(axis_x*axis_x+axis_y*axis_y+axis_z*axis_z)*sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
                results.push_back(num/den);
            } 
            return results;
    }
''')

#Compute the cosine of the angle between two vectors
r.gInterpreter.Declare('''
    float Compute_CosTheta(float axis_x, float axis_y, float axis_z, float x, float y, float z){
        float num = axis_x*x + axis_y*y + axis_z*z;
        float den = sqrt(axis_x*axis_x+axis_y*axis_y+axis_z*axis_z)*sqrt(x*x+y*y+z*z);
        return num/den;
    }
''')

#Get the maximum of a list
r.gInterpreter.Declare('''
    int Find_Max(ROOT::VecOps::RVec<float> List){
        int Max(0);
        for (size_t i=1;i<List.size();++i){
            if (List[i] > List[Max]) Max = i;
        }
        return Max;
    }
''')

#Compute cross-product
r.gInterpreter.Declare('''
    float Cross_Product_comp(int comp, float vx, float vy, float vz, float wx, float wy, float wz, int i){
    
        if (i == 1)      return vy*wz-vz*wy;
        else if (i == 2) return vz*wx-vx*wz;
        else if (i == 3) return vx*wy-vy*wx;
        else             return float(0.0);
    
    }
''')

#Compute the dot-product
r.gInterpreter.Declare('''
    float Compute_DotProduct (float Vx, float Vy, float Vz, float Wx, float Wy, float Wz){
        return Vx*Wx + Vy*Wy + Vz*Wz;
    }
''')

#Compute norm of a vector
r.gInterpreter.Declare('''

    float Get_Norm(float x, float y, float z){

        return sqrt(pow(x,2) + pow(y,2) + pow(z,2));

    }

''')

#RECO-----------------------------------------------------------------------------------------------------

#To compute the angle between both candidates for each pair of candidates
r.gInterpreter.Declare('''
    ROOT::VecOps::RVec<float> Compute_Momenta_Angles(int NTauCand, ROOT::VecOps::RVec<float> Tau_x, ROOT::VecOps::RVec<float> Tau_y, ROOT::VecOps::RVec<float> Tau_z){

        ROOT::VecOps::RVec<float> result; 
        for (size_t i=0;i<NTauCand;++i){
            for (size_t j=i+1;j<NTauCand;++j) result.push_back(Compute_CosTheta(Tau_x[i],Tau_y[i],Tau_z[i],Tau_x[j],Tau_y[j],Tau_z[j]));
        } return result;
    }
''')

#An alternative to get a single number in the angles column
r.gInterpreter.Declare('''
    float Compute_Momenta_Angles(float Tau1_x, float Tau1_y, float Tau1_z, float Tau2_x, float Tau2_y, float Tau2_z){

        return Compute_CosTheta(Tau1_x,Tau1_y,Tau1_z,Tau2_x,Tau2_y,Tau2_z);
    }
''')

#To trace back which Tau Candidates are link to which angles
r.gInterpreter.Declare('''
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> ID_Angles(int NTauCand){
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> result;
        for (size_t i = 0; i<NTauCand; ++i){
            for (size_t j = i+1; j<NTauCand; ++j){
                ROOT::VecOps::RVec<int> Pass;
                Pass.push_back(i);
                Pass.push_back(j);
                result.push_back(Pass);
            }
        }
        return result;
    }
''')


#Get the diTau system indices based on opening angle and charge
r.gInterpreter.Declare('''
    ROOT::VecOps::RVec<int> Find_diTau (int NTau23Pi, ROOT::VecOps::RVec<float> Taus_MomentaAngles, ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> IDs, ROOT::VecOps::RVec<int> Taus_q){

        ROOT::VecOps::RVec<int> diTau;

        if (NTau23Pi<2){  //Discarded, cannot do anything without at least 2 candidates
            return diTau;
        }

        else if (NTau23Pi==2){  //Check if in the same hemisphere and opposite charge, otherwise discarded
            if (Taus_MomentaAngles[0] < 0 || Taus_q[0]*Taus_q[1] > 0){
                return diTau;
            }
            else {
                diTau.push_back(IDs[0][0]);
                diTau.push_back(IDs[0][1]);
                return diTau;
            }
        }
        
        else{  //Check all the angles, Selects the biggest cosine then check if opposite sign if same sign move on to the next angle until cos < 0
            
            bool swapped; //Order the angles list and rearange the angle IDs list in the same way (easier to select the angle in case the q are not opposite)
            do{
                swapped = false;
                for (int i=1; i<NTau23Pi; ++i){
                    if (Taus_MomentaAngles[i-1] > Taus_MomentaAngles[i]){
                        
                        float tempa (Taus_MomentaAngles[i]);
                        ROOT::VecOps::RVec<int> tempi (IDs[i]);

                        Taus_MomentaAngles[i]   = Taus_MomentaAngles[i-1];
                        Taus_MomentaAngles[i-1] = tempa;

                        IDs[i]   = IDs[i-1];
                        IDs[i-1] = tempi;

                        swapped = true;
                    }
                }
            } while (swapped);
            
            for (int i=NTau23Pi-1;i>=0;--i){
                if (Taus_MomentaAngles[i] < 0 || Taus_q[IDs[i][0]]*Taus_q[IDs[i][1]] > 0) continue; //Could check directly that the maximum cosine is negative and get out of the function, the maximum cosine will always be positive since starting from 3 candidates two are always found on the same side
                else{
                    diTau.push_back(IDs[i][0]);
                    diTau.push_back(IDs[i][1]);
                    break;
                }
            }
            return diTau;
        }
    }
''')

#To Compute Bs visible vertex, that is using the reco tau dv and infering the flight direction from visible momenta, the visible vertex is the midpoint of the line joining the 2 skewed lines
r.gInterpreter.Declare('''
    ROOT::VecOps::RVec<float> Compute_BsVisibleVertex(float Ver1x, float Ver1y, float Ver1z, float Ver2x, float Ver2y, float Ver2z, float p1x, float p1y, float p1z, float p2x, float p2y, float p2z){
        ROOT::VecOps::RVec<float> Bs_Vertex;

        float n_x = p1y*p2z - p1z*p2y;
        float n_y = p1z*p2x - p1x*p2z;
        float n_z = p1x*p2y - p1y*p2x; 

        float t1 = ((p2y*n_z-p2z*n_y)*(Ver2x-Ver1x) + (p2z*n_x-p2x*n_z)*(Ver2y-Ver1y) + (p2x*n_y-p2y*n_x)*(Ver2z-Ver1z))/(n_x*n_x + n_y*n_y + n_z*n_z);
        float t2 = ((p1y*n_z-p1z*n_y)*(Ver2x-Ver1x) + (p1z*n_x-p1x*n_z)*(Ver2y-Ver1y) + (p1x*n_y-p1y*n_x)*(Ver2z-Ver1z))/(n_x*n_x + n_y*n_y + n_z*n_z);

        float Pt1x = t1*p1x + Ver1x;
        float Pt1y = t1*p1y + Ver1y;
        float Pt1z = t1*p1z + Ver1z;

        float Pt2x = t2*p2x + Ver2x;
        float Pt2y = t2*p2y + Ver2y;
        float Pt2z = t2*p2z + Ver2z;

        Bs_Vertex.push_back(Pt2x + (Pt1x-Pt2x)/2);
        Bs_Vertex.push_back(Pt2y + (Pt1y-Pt2y)/2);
        Bs_Vertex.push_back(Pt2z + (Pt1z-Pt2z)/2);
        return Bs_Vertex;
    }
''')

#To evaluate the quality of the RECO Bs vertex, based on the length of the line of minimal distance between the two skewed lines
r.gInterpreter.Declare('''
    float Quality_BsVisibleVertex(float Ver1x, float Ver1y, float Ver1z, float Ver2x, float Ver2y, float Ver2z, float p1x, float p1y, float p1z, float p2x, float p2y, float p2z){
        float result;

        float n_x = p1y*p2z - p1z*p2y;
        float n_y = p1z*p2x - p1x*p2z;
        float n_z = p1x*p2y - p1y*p2x; 

        float t1 = ((p2y*n_z-p2z*n_y)*(Ver2x-Ver1x) + (p2z*n_x-p2x*n_z)*(Ver2y-Ver1y) + (p2x*n_y-p2y*n_x)*(Ver2z-Ver1z))/(n_x*n_x + n_y*n_y + n_z*n_z);
        float t2 = ((p1y*n_z-p1z*n_y)*(Ver2x-Ver1x) + (p1z*n_x-p1x*n_z)*(Ver2y-Ver1y) + (p1x*n_y-p1y*n_x)*(Ver2z-Ver1z))/(n_x*n_x + n_y*n_y + n_z*n_z);

        float Pt1x = t1*p1x + Ver1x;
        float Pt1y = t1*p1y + Ver1y;
        float Pt1z = t1*p1z + Ver1z;

        float Pt2x = t2*p2x + Ver2x;
        float Pt2y = t2*p2y + Ver2y;
        float Pt2z = t2*p2z + Ver2z;

        result = sqrt(pow(Pt1x-Pt2x,2) + pow(Pt1y-Pt2y,2) + pow(Pt1z-Pt2z,2));
        return result;
    }
''')

#IP computation  ----------------------------------------------------------------------------------------------------------

#Get the point on the flight direction that is the closest to the tau->3pi candidates
r.gInterpreter.Declare('''
    ROOT::VecOps::RVec<float> Compute_IPP(float Vx, float Vy, float Vz, float SVx, float SVy, float SVz, float PVx, float PVy, float PVz){
        ROOT::VecOps::RVec<float> result;
        ROOT::VecOps::RVec<float> FlightDir = {SVx-PVx,SVy-PVy,SVz-PVz};
        ROOT::VecOps::RVec<float> ToProject = {Vx-PVx,Vy-PVy,Vz-PVz};

        float t = Compute_DotProduct(ToProject[0],ToProject[1],ToProject[2],FlightDir[0],FlightDir[1],FlightDir[2])/pow(Get_Norm(FlightDir[0],FlightDir[1],FlightDir[2]),2);
        result.push_back(PVx+t*FlightDir[0]);
        result.push_back(PVy+t*FlightDir[1]);
        result.push_back(PVz+t*FlightDir[2]);
        return result;
    }
''')

#==========================================================================================================================================

#The functions that creates the columns used by the BDT in the current dataframe
#Proceed in two stage, first properly ID's the di-Tau system, should work in most cases but adds a check on opposite charge tau's candidates, which might reject some events since the check in stage 1 is that we have exactly two tau's in the signal hemi