#include "Pani.h"


extern "C"
int planJobs(int& counter, std::vector<Machine>** out_machines, std::vector<Job>* _jq, std::vector<Machine>* _mq, float _epsilon) {
    bool successful = Pani::runHybrid(counter, out_machines, _jq, _mq, _epsilon, Pani::DP);
    if(!successful) {
        counter++;
    }
    return 1;
}

