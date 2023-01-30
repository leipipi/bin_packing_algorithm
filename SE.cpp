#include "Job.h"
#include "Machine.h"
#include "MathFunctions.h"
#include <vector>
#include <cmath>
#include <algorithm>
/**
 * Simulated Evolution by SS17.
 */
struct CompareFuncBo{
    bool operator()(Machine &a, Machine &b) {
        return a.getSumOfLoad() < b.getSumOfLoad();
    }
};
struct CompareFunc {
    
    bool operator()(Job const &a, Job const &b) {
        return a.getSumOfComponents() > b.getSumOfComponents();
    }
};
float getGi(float* mv, float* jv, int dimensions) {
    float *nmv = new float[dimensions];
    for(unsigned i = 0; i < dimensions; i++)
    {
        nmv[i] = mv[i];
    }
    float smv = 0.;
    float sjv = 0.;
    float gi;
    MathFunctions::VectorAddition(nmv, jv, dimensions);
    for(unsigned i = 0; i < dimensions; ++i) {
        sjv += jv[i]/nmv[i];
    }
    delete[] nmv;
    gi = sjv/dimensions;
    return gi;
};
extern "C"
int planJobs(int &counter, std::vector <Machine> **out_machines, std::vector <Job> *_jq, std::vector <Machine> * _mq, float _epsilon) {
    std::vector <Job*>* aJQ = new std::vector <Job*> ();
    int SizeofMachine = 100000;
    std::sort(_jq->begin(), _jq->end(), CompareFunc());
    for(int i = 0; i < _jq->size(); i++) {
        aJQ->push_back(&(_jq->at(i)));
    }
    int maxSelection = int(0.4*aJQ->size());
    int dimensions = aJQ->at(0)->getDimensions();
    MachineGenerator mach_gen(dimensions);
    std::vector <Machine>* mq = mach_gen.generateMachines(1, 1.0);
    // initial assignment FFD
    for(unsigned i = 0; i < aJQ->size(); ++i) {
        for(unsigned j = 0; j < mq->size(); ++j) {
            if(mq->at(j).assignJob((aJQ->at(i)))) {
                break;
            }
            else {
                if(j == mq->size() - 1) {
                    mq->push_back(mach_gen.generateMachine(mq->size() + 1, 1.0));
                    mq->at(j).assignJob((aJQ->at(i)));
                }
            }
        }
    }
    int initialSize = mq->size();
    int T = 100;
    // start SE
    while(T>0){
        //Evaluation
        float **Gi;
        float* tmp_jv;
        float* tmp_mv;
        int mq_size = mq->size();
        Gi = new float*[mq_size];
        //Define Gi
        for (int i = 0; i < mq->size(); ++i)
        {
            tmp_mv = mq->at(i).getRemainingCapacity();
            //define goodness function for all items
            Gi[i] = new float[mq->at(i).assignedJobs->size()];
            for(int j = 0;j < mq->at(i).assignedJobs->size(); ++j){
                tmp_jv = mq->at(i).assignedJobs->at(j).getKPIVec();
                Gi[i][j] = getGi(tmp_mv, tmp_jv, dimensions);
            }
        }
         std::vector <Job>* JQ = new std::vector <Job> ();
        int numOfSelected = 0;
        //Selection
        int itr = 0;
        while(numOfSelected < maxSelection){
            for (int i = mq->size()-1; i >= 0; --i)
            {
                float * tmpp;
                int mqIsize = mq->at(i).assignedJobs->size()-1;
                for(int j = mqIsize; j >= 0; --j){
                    float probability = rand()%100;
                    probability = probability/100;
                    if(probability <= (1-Gi[i][j])){
                        numOfSelected++;
                        JQ->push_back((mq->at(i).assignedJobs->at(j)));
                        mq->at(i).removeJob(j);
                    }
                }
                if(mq->at(i).assignedJobs->size() == 0)
                {
                    mq->erase(mq->begin()+i);
                }
            }
            itr++;
        }
        //Allocation
        std::sort(JQ->begin(), JQ->end(), CompareFunc());
        if(mq->size() == 0)
            mq->push_back(mach_gen.generateMachine(mq->size() + 1, 1.0));
        std::sort(mq->begin(), mq->end(), CompareFuncBo());
        for(unsigned i = 0; i < JQ->size(); ++i) {
            for(unsigned j = 0; j < mq->size(); ++j) {
                if(mq->at(j).assignJob(&(JQ->at(i)))) {
                    break;
                }
                else {
                    if(j == mq->size() - 1) {
                        mq->push_back(mach_gen.generateMachine(mq->size() + 1, 1.0));
                        mq->at(j).assignJob(&(JQ->at(i)));
                    }
                }
            }
        }
        --T;
        for(int i = mq_size-1; i >=0; i--)
        {
             delete [] Gi[i];
        }
        delete [] Gi;
        delete JQ;
        if(mq->size() < SizeofMachine){
            SizeofMachine = mq->size();
        }
    }
    for(unsigned i = 0; i < aJQ->size(); i++) {
	    aJQ->at(i)->deleteData();
    }
    delete aJQ;
    //failed
    if(SizeofMachine > _mq->size()){
    	counter++;
    }
    *out_machines = mq;
    return 1;
}


