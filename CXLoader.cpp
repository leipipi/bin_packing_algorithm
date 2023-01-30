/* 
The interval limits can be set in config file by setting the values: interval_min and interval_max
*/

#include <stdio.h>
#include <unistd.h>
#include <vector>
#include <algorithm>
#include <random>

#include "Job.h"
#include "Environment.h"

extern "C"
std::vector<Job>* loadJobs(const char* _path, JobMetaData* _jobdata, Environment* _env, int _runs) {
    std::vector<Job>* jobVec = new std::vector<Job>();

#ifdef RANDOM
    std::random_device device;
    std::mt19937 generator(device());
#else
    char* seed = getenv("RAND_SEED");
    std::mt19937 generator;
    if(!seed)
        generator = std::mt19937(1234);
    else
        generator = std::mt19937(atoi(seed));
#endif

    std::uniform_int_distribution<int> distBurst(_jobdata->lowerBurst, _jobdata->upperBurst);
    std::uniform_int_distribution<int> distArrive(0, _jobdata->duration);
    std::uniform_real_distribution<float> noise(_env->intervalMin, _env->intervalMax);

    // use vectors to generate jobs

    std::vector <float> dummy;

    for(unsigned i = 0; i < _runs * _jobdata->numberOfJobs; i++) {
        dummy.clear();
    	for(int i = 0; i < _jobdata->dimensions; i++) {
            dummy.push_back(noise(generator));
        }
        jobVec->push_back(Job(i + 1, &dummy, distBurst(generator), distArrive(generator)));
    }

    return jobVec;
}

