#include "Kou.h"


extern "C"
int planJobs(int &counter, std::vector<Machine> **out_machines, std::vector<Job> *jq, std::vector<Machine> *mq, float _eps) {
    bool successful = Kou::runKouFFD(counter, out_machines, jq, mq, _eps, Sorting::NO_SORTING);
    if(!successful) {
        counter++;
        return 0;
    }
    return 1;
}


