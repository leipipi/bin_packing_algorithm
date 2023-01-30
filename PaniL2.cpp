#include <iostream>
#include <sstream>
#include "Pani.h"
#include "Kou.h"


extern "C"
int planJobs(int &counter, std::vector <Machine> **out_machines, std::vector<Job> *_jq, std::vector<Machine> * _mq, float _epsilon) {

    bool successful = Pani::runPanigrahy(counter, out_machines, _jq, _mq, _epsilon, Pani::L2);
    if(!successful) {
        counter++;
    }
    return 1;
}    
    
