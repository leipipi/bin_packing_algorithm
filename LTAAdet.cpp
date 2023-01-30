#include "Job.h"
#include "JobConfiguration.h"
#include "Machine.h"
#include "MathFunctions.h"
#include "LTAA.h"
#include <cmath>
#include <random>

#include <iostream>
#include <chrono>
#include <ctime> 




extern "C"
int planJobs(int &_counter, std::vector<Machine> **out_machines, std::vector<Job> *jq, std::vector<Machine>* _mq, float _epsilon) {
    
//    auto start1 = std::chrono::system_clock::now();
    
    // Dimensions of vectors
    int dimensions = jq->at(0).getDimensions();
    // Round all jobs; rounded values are stored in job
    for(int i=0; i<jq->size(); i++) {
        jq->at(i).roundLTAA(_epsilon);
    }
    // Categorize vectors into zero vectors, small vectors and large vectors
    std::vector<Job>* null_vecs = new std::vector <Job>();
    std::vector<std::vector<Job>*>* S = new std::vector<std::vector<Job>*>();
    std::vector<std::vector<Job>*>* L = new std::vector<std::vector<Job>*>();
    LTAA::categorizeJobs(jq, null_vecs, S, L);
//    LTAA::printVectors(S, L, null_vecs);
    // Generate job configurations
    std::vector<JobConfiguration>* configs = new std::vector<JobConfiguration>();
    LTAA::generateConfigurations(L, configs, _epsilon, dimensions);
    // MILP
    if(!LTAA::solveMILP(configs, S, L, _mq, _epsilon, false)) {
        _counter++;
    }
    // Assign zero vectors
    if(!null_vecs->empty()) {
        for(unsigned i=0; i<null_vecs->size(); i++) {
            for(unsigned j=0; j<_mq->size(); j++) {
                if(_mq->at(j).assignJobLTAA(&(null_vecs->at(i)), false)) {
                    break;
                }
                if(_mq->size() - j == 1) {
                    return false;
                }
            }
        }
    }
    // Clean up
    delete configs;
    for(int i=0; i<S->size(); i++) {
        delete S->at(i);
    }
    delete S;
    for(int i=0; i<L->size(); i++) {
        delete L->at(i);
    }
    delete L;
    delete null_vecs;
    // Create output and return
    std::vector<Machine>* mq_clone = new std::vector<Machine>(*_mq);
    *out_machines = mq_clone;
#ifdef DEBUG
    printf("out_machines size: %lu\n", (*out_machines)->size());
#endif
    return 1;
}

