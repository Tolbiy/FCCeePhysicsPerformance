#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include <string>
#include <vector>

void Compute_PreBDTeff() {
    
    std::vector<string> Modes = {"sig","bb","cc","ss","ud"};

    for (size_t m=0; m<Modes.size(); ++m){
        int NF (100);
        if (Modes.at(m) == "sig") NF = 20;
        
        int Processed (0);
        int Selected (0);

        for (size_t i=0; i<NF; ++i)
        {
            string filename ("/eos/experiment/fcc/ee/analyses_storage/flavor/Bs2TauTau/flatNtuples/winter2023/analysis_stage1_Leptons_withCuts/");
            if (Modes.at(m) == "sig") filename += "p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau/chunk_"+std::to_string(i)+".root";
            else filename += "p8_ee_Z"+Modes.at(m)+"_ecm91/chunk_"+std::to_string(i)+".root";
            
            TFile *file = TFile::Open(filename.c_str(),"READ");
            if (!file || file->IsZombie()) {
                std::cerr << "Error opening file #" << i << " in mode " << Modes.at(m) << endl;
                //exit(-1);
            }
            else{
                TParameter<int> *Processedi = nullptr;
                file->GetObject("eventsProcessed",Processedi);
                //std::cout << Processedi->GetVal() << "/";
                Processed += Processedi->GetVal();

                TParameter<int> *Selectedi = nullptr;
                file->GetObject("eventsSelected",Selectedi);
                //std::cout << Selectedi->GetVal() << std::endl;
                Selected += Selectedi->GetVal();
            }
        }

        std::cout << "_________ Modes = " << Modes.at(m) << " __________" << std::endl;
        std::cout << "Processed: " << Processed << std::endl;
        std::cout << "Selected:  " << Selected << std::endl;
        std::cout << "___________________________" << std::endl;

    }
}

