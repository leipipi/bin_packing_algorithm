#include "Kou.h"

/**
 * Vector Scheduling. Algorithm by Kou and Markowsky.
 * First fit; jobs sorted lexicographically in descending order:
 * b >= a if b-a = 0 or if the first component != 0 of b-a is positive.
 * @param counter
 * @param out_machines
 * @param jq
 * @param mq
 * @return 
 */
extern "C"
int planJobs(int &counter, std::vector<Machine> **out_machines, std::vector<Job> *jq, std::vector<Machine> *mq, float _eps) {
    bool successful = Kou::runKouFFD(counter, out_machines, jq, mq, _eps, Sorting::LEXI);
    if(!successful) {
        counter++;
        return 0;
    }
    return 1;
}

