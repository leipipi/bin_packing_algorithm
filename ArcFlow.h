#include "Job.h"
#include "Machine.h"
#include "MathFunctions.h"
#include "Sorting.h"
#include "gurobi_c++.h"
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#define FIRSTMVPSHELLSCRIPT "\
#/bin/bash \n\
vpsolver/bin/vbp2afg instan.mvp graph.afg\n\
vpsolver/bin/afg2mps graph.afg model.mps\n\
"
#define CLEARFILESCRIPT "\
#/bin/bash \n\
rm -fr instan.mvp\n\
rm -fr graph.afg\n\
rm -fr model.mps\n\
rm -fr vars.sol\n\
rm -fr result.txt\n\
"
#define FIRSTVBPSHELLSCRIPT "\
#/bin/bash \n\
vpsolver/bin/vbp2afg instan.vbp graph.afg\n\
vpsolver/bin/afg2mps graph.afg model.mps\n\
"
#define SECONDSHELLSCRIPT "\
#/bin/bash \n\
vpsolver/bin/vbpsol graph.afg vars.sol \n\
"

class ArcFlow {
public:

    static bool CheckInput(std::vector <Machine> * _mq, std::vector <Job> *jq) {
        float* tmp_jv;
        float* tmp_mv;

        int dims = jq->at(0).getDimensions();
        for(unsigned i = 0; i < jq->size(); i++) {
            int flag = 1;
            tmp_jv = jq->at(i).getKPIVec();
            for(unsigned j = 0; j < _mq->size(); j++) {
                tmp_mv = _mq->at(j).getRemainingCapacity();
                //      printf("mq_ %d remaining capactiry is %f  %f  %f\n",j,tmp_mv[0],tmp_mv[1],tmp_mv[2]);
                if(MathFunctions::isVectorDifferencePositive(tmp_mv, tmp_jv, dims))
                    flag = 0;
            }
            if(flag)
                return false;
        }
        return true;
    }

    static void CreateMVPInputFile(std::vector <Job> *jq, std::vector <Machine> * mq) {

        int dims = jq->at(0).getDimensions();
        std::ofstream destFile("instan.mvp");
        destFile << dims << "\n";
        destFile << "\n";
        destFile << mq->size() << "\n";
        int i = 0;
        float* tmp_mv;
        //machine input
        for(unsigned i = 0; i < mq->size(); i++) {
            tmp_mv = mq->at(i).getRemainingCapacity();
            for(unsigned j = 0; j < dims; j++) {
                //    printf("machine capacity dimens 1 is %d   ",int(tmp_mv[j]*10000));
                destFile << int(tmp_mv[j]*10000) << " ";
            }
             bool flag = true;
            // check whether this is a new machine
            for(unsigned k = 0; k < dims; k++){
                if(tmp_mv[k] < 0.9998){
                    flag = false;
                    break;
                }
            }
            if(flag)
                destFile << "1000" << " " << "1" << "\n";
            else
                destFile << "1" << " " << "1" << "\n";
        //    destFile << "1" << " " << "1" << "\n";
        }

        destFile << "\n";
        destFile << jq->size() << "\n";
        destFile << "\n";
        float* tmp_jv;
        for(int i = 0; i < jq->size(); i++) {
            destFile << "1" << " " << "1" << "\n";
            tmp_jv = jq->at(i).getKPIVec();
            for(unsigned j = 0; j < dims; j++) {
                if(j == dims - 1)
                    destFile << int(10000 * tmp_jv[j]) << "\n";
                else
                    destFile << int(10000 * tmp_jv[j]) << " ";
            }
            destFile << "\n";
        }
        destFile.close();
    }

    static void CreateVBPInputFile(std::vector <Job> *jq) {
        int dims = jq->at(0).getDimensions();
        std::ofstream destFile("instan.vbp");
        destFile << dims << "\n";
        int i = 0;
        while(i < dims - 1) {
            destFile << "10000" << " ";
            i++;
        }
        destFile << "10000" << "\n";
        destFile << jq->size() << "\n";
        float* tmp_jv;
        for(int i = 0; i < jq->size(); i++) {
            tmp_jv = jq->at(i).getKPIVec();
            for(unsigned j = 0; j < dims; j++) {
                if(j == dims - 1)
                    destFile << int(10000 * tmp_jv[j]) << " " << "1" << "\n";
                else
                    destFile << int(10000 * tmp_jv[j]) << " ";
            }
        }
        destFile.close();
    }

    static int SolveMPSmodel() {
        GRBEnv *env = 0;
        env = new GRBEnv();
        env->set("OutputFlag", "0");
        GRBModel model = GRBModel(*env, "model.mps");
        model.set(GRB_IntParam_OutputFlag, 0);
        if(model.get(GRB_IntAttr_IsMIP) == 0) {
            throw GRBException("Model is not a MIP");
        }
        model.optimize();
        try {
            model.write("vars.sol");
        } catch(GRBException e) {
//            printf("GRBException in SolveMPSModel\n");
            return 0;
        }
        printf("MPS model solved\n");
        return 1;
    }

    static int AssignVBPJobs(std::vector <Machine> * mq, std::vector <Job> *jq) {
        std::ifstream infile("result.txt");
        std::string line;
        int _mqsize = mq->size();
        int count = 0;
        while(std::getline(infile, line)) {
            if(count != 0) {
                std::istringstream iss(line);
                int jobid;
                if(count > _mqsize) {
                    return 0;
                }
                while(iss >> jobid) {
                    mq->at(count - 1).assignJob(&(jq->at(jobid - 1)));
                }
            }
            count++;
        }
        infile.close();
        return 1;
    }

    static void AssignMVPJobs(std::vector <Machine> * mq, std::vector <Job> *jq) {
        std::ifstream infile("result.txt");
        std::string line;
        while(std::getline(infile, line)) {
            std::istringstream iss(line);
            int jobid;
            int machineid;
            iss >> machineid;
            while(iss >> jobid) {
                mq->at(machineid).assignJob(&(jq->at(jobid - 1)));
            }
        }
        infile.close();
    }

    //All machines should be empty

    static int ArcFlowToSolveVBP(std::vector <Machine> * _mq, std::vector <Job> *jq) {
        //     CreateMVPInputFile(jq, _mq);
        //   printf("CreateVBPInputFile\n");
        CreateVBPInputFile(jq);
        //   printf("FIRSTVBPSHELLSCRIPT\n");
        system(FIRSTVBPSHELLSCRIPT);
        //   printf("SolveMPSmodel\n");
        SolveMPSmodel();
        printf("SECONDSHELLSCRIPT\n");
        system(SECONDSHELLSCRIPT);
        printf("finish second shell\n");
        if(AssignVBPJobs(_mq, jq)) {
            printf("Finish VBP\n");
            return 1;
        } else {
            printf("fail to assign\n");
            return 0;
        }
    }

    static int ArcFlowToSolveMVP(std::vector <Machine> * _mq, std::vector <Job> *jq) {

        system(CLEARFILESCRIPT);
        bool flag = false;
        flag = CheckInput(_mq, jq);

        if(flag) {
            //    printf("create MVP file\n ");
            CreateMVPInputFile(jq, _mq);
            //    printf("Finish create MVP file\n ");
            //     printf("start first shell\n ");
            system(FIRSTMVPSHELLSCRIPT);
            //    printf("Finish first shell\n ");
            if(SolveMPSmodel()) {
                //     printf("Start last shell\n");
                system(SECONDSHELLSCRIPT);
                //     printf("Finish last shell\n");
                //     printf("Start Assign MVP jobs\n");
                AssignMVPJobs(_mq, jq);
                //     printf("Finish Assign MVP jobs\n");
                return 1;
            } else {
                printf("MPS model cannot be solved\n");
                return 0;
            }
        } else {
            printf("input cannot be solved\n");
            return 2;
        }
    }
};

