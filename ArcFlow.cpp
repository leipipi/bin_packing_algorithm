#include "Job.h"
#include "Machine.h"
#include "MathFunctions.h"
#include "gurobi_c++.h"
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#define FIRSTSHELLSCRIPT "\
#/bin/bash \n\
vpsolver/bin/vbp2afg instan.vbp graph.afg\n\
vpsolver/bin/afg2mps graph.afg model.mps\n\
"
#define SECONDSHELLSCRIPT "\
#/bin/bash \n\
vpsolver/bin/vbpsol graph.afg vars.sol \n\
"
//gurobi_cl ResultFile=vars.sol model.mps\n\

extern "C"
int planJobs(int &counter, std::vector <Machine> **out_machines, std::vector <Job> *jq, std::vector <Machine> * _mq, float _epsilon) {
    
    int dims = jq->at(0).getDimensions();
    float* tmp_jv;
    float* tmp_mv;
    float* total = new float[dims];
    for(unsigned i = 0; i < dims; i++)
    {
        total[i] = 0;
    }
    for(int i = 0; i < jq->size(); i++) {
        tmp_jv = jq->at(i).getKPIVec();
        for(unsigned j = 0; j < dims; j++)
        {
            total[j] += tmp_jv[j];
        }
    }
    MachineGenerator mach_gen(dims);
    std::vector <Machine>* mq = mach_gen.generateMachines(0,1.0);

    //Change input format into file instan.vbp
    std::ofstream destFile("instan.vbp");
    destFile<<dims<<"\n";
    for(unsigned j = 0; j < dims; j++){
        if(j != dims-1)
            destFile<<"100001"<<" ";
        else
            destFile<<"100001"<<"\n";
    }
    destFile<<jq->size()<<"\n";
    
    for(int i = 0; i < jq->size(); i++) {
        tmp_jv = jq->at(i).getKPIVec();
        for(unsigned j = 0; j < dims; j++){
            destFile<<int(100000*tmp_jv[j])<<" ";
        }
        destFile<<"1"<<"\n";
    }
    destFile.close();
    //Shell script
    system(FIRSTSHELLSCRIPT);
    //Solve mps model, clean Gurobi's print
    try {
        GRBEnv *env = 0;
        env = new GRBEnv();
        //std::string para = "OutputFlag";
        //std::string newVal = "0";
        //env->set(para, newVal);
        env->set(GRB_IntParam_OutputFlag, 0);
        GRBModel model = GRBModel(*env, "model.mps");
        model.set(GRB_IntParam_OutputFlag,0);
        if (model.get(GRB_IntAttr_IsMIP) == 0) {
            throw GRBException("Model is not a MIP");
        }
        model.optimize();
        int error;
        //Write solution into file vars.sol
        model.write("vars.sol");
        //Solve the solution
        system(SECONDSHELLSCRIPT);
        //Read result file and assign jobs
        std::ifstream infile( "result.txt" );
        std::string line;
        int count = 0;
        int i = 0;
        while (std::getline(infile, line))
        {
            if(count!=0)
            {
                mq->push_back(mach_gen.generateMachine(mq->size() + 1, 1.0));
                std::istringstream iss(line);
                int jobid;
                while(iss >> jobid)
                {
                    mq->at(count-1).assignJob(&(jq->at(jobid-1)));
                }
            }
            count++;
        }
        infile.close();
        if(mq->size() > _mq->size()) {
            counter++;
        }
        *out_machines = mq;
        return 1;
    } catch(GRBException e) {
        printf("Error code = %d\n", e.getErrorCode());
        printf("GurobiError: %s\n", e.getMessage().c_str());
        return 1;
    }
}


