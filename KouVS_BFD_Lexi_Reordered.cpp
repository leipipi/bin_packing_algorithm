#include "Kou.h"


extern "C"
int planJobs(int &counter, std::vector<Machine> **out_machines, std::vector<Job> *jq, std::vector<Machine> *mq, float _eps) {
    bool successful = Kou::runKouBFD(counter, out_machines, jq, mq, _eps, Sorting::LEXI_REORDERED);
    if(!successful) {
        counter++;
        return 0;
    }
    return 1;
}

