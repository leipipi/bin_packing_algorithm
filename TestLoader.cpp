
#include <stdio.h>
#include <unistd.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <random>
#include <limits>

#include "Job.h"
#include "Environment.h"

class TestVectorCreator {
public:

    TestVectorCreator() {
    }

    std::vector<std::vector<float>*>* createVectors() {
        std::vector<std::vector<float>*>* queue = new std::vector<std::vector<float>*>();
        int n = 55;
        /*float arr[n][3] = {0.15, 0.6, 0.5, 
            0.1, 0.3, 0.4,
            0.2, 0.6, 0.65,
            0.5, 0.2, 0.4,
            0.2, 0.2, 0.2,
            0.3, 0.2, 0.6,
            0.001, 0.001, 0.001,
            0.0001, 0.0001, 0.0001,
            0.7, 0.6, 0.2,
            0.3, 0.7, 0.4,
            0.3, 0.1, 0.4,
            0.2, 0.4, 0.2,
            0.7, 0.1, 0.2,
            0.2, 0.4, 0.2,
            0.000098642, 0.000222, 0.00003,
            0.3, 0.3, 0.0,
            0.000065821, 0.00007733, 0.0001898,
            0.1, 0.0, 0.2};*/
        /*float arr[n][3] = {0.332119, 0.453951, 0.452049,
            0.111607, 0.062492, 0.015700,
            0.224890, 0.481214, 0.400552,
            0.379173, 0.412910, 0.296092,
            0.470631, 0.219494, 0.441325,
            0.290518, 0.045938, 0.290334,
            0.217204, 0.372747, 0.467408,
            0.392509, 0.203318, 0.412532,
            0.828791, 0.372807, 0.065226,
            0.243506, 0.320526, 0.442990,
            0.359168, 0.472996, 0.111308,
            0.296973, 0.146231, 0.436123,
            0.817396, 0.396736, 0.365139,
            0.435872, 0.289045, 0.223725,
            0.227363, 0.297880, 0.406907,
            0.256110, 0.474134, 0.230483,
            0.396314, 0.081267, 0.448638,
            0.197765, 0.251269, 0.281567,
            0.025069, 0.433326, 0.490609,
            0.458429, 0.273608, 0.160627,
            0.360489, 0.197126, 0.224832,
            0.208729, 0.465938, 0.229472,
            0.443136, 0.046245, 0.373267,
            0.311908, 0.466956, 0.117453,
            0.287870, 0.107752, 0.410989,
            0.062129, 0.457587, 0.335977,
            0.495643, 0.120738, 0.050751,
            0.097365, 0.389048, 0.378559,
            0.314809, 0.150661, 0.244059,
            0.204715, 0.456759, 0.057876,
            0.014581, 0.220436, 0.451318,
            0.155608, 0.489759, 0.038612,
            0.287609, 0.176897, 0.115770,
            0.241661, 0.029004, 0.265823,
            0.022072, 0.422766, 0.349821,
            0.234729, 0.278154, 0.120408,
            0.012540, 0.382196, 0.341582,
            0.271430, 0.212852, 0.104247,
            0.088180, 0.485753, 0.076546,
            0.045264, 0.040954, 0.427246,
            0.039307, 0.497529, 0.100301,
            0.208004, 0.178701, 0.116471,
            0.188189, 0.190253, 0.128679,
            0.039307, 0.497529, 0.100301,
            0.015879, 0.495200, 0.108503};*/
        float arr[n][3] = {0.99, 0.99, 0.0,
            0.99, 0.99, 0.0,
            0.99, 0.99, 0.0,
            0.99, 0.99, 0.0,
            0.99, 0.99, 0.0,
            0.0, 0.0, 0.01949,
            0.0, 0.0, 0.01948, 
            0.0, 0.0, 0.01947, 
            0.0, 0.0, 0.01946, 
            0.0, 0.0, 0.01945, 
            0.0, 0.0, 0.01944, 
            0.0, 0.0, 0.01943, 
            0.0, 0.0, 0.01942, 
            0.0, 0.0, 0.01941, 
            0.0, 0.0, 0.01940, 
            0.0, 0.0, 0.01939,
            0.0, 0.0, 0.01938, 
            0.0, 0.0, 0.01937, 
            0.0, 0.0, 0.01936, 
            0.0, 0.0, 0.01935, 
            0.0, 0.0, 0.01934, 
            0.0, 0.0, 0.01933, 
            0.0, 0.0, 0.01932, 
            0.0, 0.0, 0.01931, 
            0.0, 0.0, 0.01930, 
            0.0, 0.0, 0.01929,
            0.0, 0.0, 0.01928, 
            0.0, 0.0, 0.01927, 
            0.0, 0.0, 0.01926, 
            0.0, 0.0, 0.01925, 
            0.0, 0.0, 0.01924, 
            0.0, 0.0, 0.01923, 
            0.0, 0.0, 0.01922, 
            0.0, 0.0, 0.01921, 
            0.0, 0.0, 0.01920, 
            0.0, 0.0, 0.01919,
            0.0, 0.0, 0.01918, 
            0.0, 0.0, 0.01917, 
            0.0, 0.0, 0.01916, 
            0.0, 0.0, 0.01915, 
            0.0, 0.0, 0.01914, 
            0.0, 0.0, 0.01913, 
            0.0, 0.0, 0.01912, 
            0.0, 0.0, 0.01911, 
            0.0, 0.0, 0.01910, 
            0.0, 0.0, 0.01909,
            0.0, 0.0, 0.01908, 
            0.0, 0.0, 0.01907, 
            0.0, 0.0, 0.01906, 
            0.0, 0.0, 0.01905, 
            0.0, 0.0, 0.01904, 
            0.0, 0.0, 0.01903, 
            0.0, 0.0, 0.01902, 
            0.0, 0.0, 0.01901, 
            0.0, 0.0, 0.01900};


        for(int i = 0; i < n; i++) {
            std::vector<float>* vec = new std::vector<float>();
            for(int j = 0; j < 3; j++) {
                vec->push_back(arr[i][j]);
            }
            queue->push_back(vec);
        }

        return queue;
    }
};

extern "C"
std::vector <Job>* loadJobs(const char* _path, JobMetaData* _jobdata, Environment* _env, int _runs) {
    std::vector <Job>* jobVec = new std::vector <Job>();

    std::random_device device;
    std::mt19937 generator(device());

    std::uniform_int_distribution<int> distBurst(_jobdata->lowerBurst, _jobdata->upperBurst);
    std::uniform_int_distribution<int> distArrive(0, _jobdata->duration);

    TestVectorCreator creator = TestVectorCreator();

    std::vector<std::vector<float>*>* vectors = creator.createVectors();

    // use vectors to generate jobs
    for(unsigned i = 0; i < vectors->size(); i++) {
        jobVec->push_back(Job(i + 1, vectors->at(i), distBurst(generator), distArrive(generator)));
    }
    for(unsigned i = 0; i < vectors->size(); i++) {
        delete vectors->at(i);
    }
    delete vectors;

    return jobVec;
}

