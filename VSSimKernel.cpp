#include <cstdio>
#include <cstdlib>
#include <dlfcn.h>
#include <vector>
#include <random>
#include <ostream>
#include <chrono>
#include "Job.h"
#include "Environment.h"
#include "Machine.h"
#include "ConfigReader.h"
#include "CSVWriter.h"
#include "Timer.h"

int main(int argc, char** argv) {
    if(argc != 2) {
        fprintf(stderr, "usage: %s <path to config file> \nexit\n", argv[0]);
        return 1;
    }

    ConfigReader cr = ConfigReader(argv[1]);
    Environment* env = cr.getEnvironment();
    JobMetaData* jmd = cr.getJobMetaData();

    printf("Opening loader: %s ...\n", (env->loader).c_str());
    void* handle = dlopen((env->loader).c_str(), RTLD_LAZY);
    if(!handle) {
        fprintf(stderr, "Cannot open library: %s\n", dlerror());
        return 1;
    }

    printf("Loading symbol loadJobs ...\n");
    typedef std::vector <Job>* (*JobLoader_t)(const char*, JobMetaData*, Environment*, int);

    // Reset errors    
    dlerror();
    JobLoader_t jl = (JobLoader_t) dlsym(handle, "loadJobs");

    const char *dlsym_error = dlerror();
    if(dlsym_error) {
        fprintf(stderr, "Cannot load symbol 'loadJobs': %s\n", dlsym_error);
        dlclose(handle);
        return 1;
    }

    printf("Calling loadJobs ...\n");
    std::vector <Job>* jv = jl(env->file.c_str(), jmd, env, env->runs); // close the library    

    //printf("Closing library...\n");
    printf("    Jobs loaded? ... ");
    if(!jv) {
        printf("No jobs loaded! ...\n");
    } else {
        printf("%lu jobs loaded! ...\n", jv->size());
    }
    dlclose(handle);

    printf("Opening algorithm: %s ...\n", (env->algorithm).c_str());
    handle = dlopen((env->algorithm).c_str(), RTLD_LAZY);
    if(!handle) {
        fprintf(stderr, "Cannot open algorithm %s: %s\n", (env->algorithm).c_str(), dlerror());
        return 1;
    }
    //load symbols
    //int planJobs(int &counter, std::vector <Machine> *out_machines, std::vector <Job> *jq, std::vector <Machine> * mq)
    printf("Loading symbol planJobs...\n");
    typedef int (*JobPlaner_t)(int&, std::vector<Machine>**, std::vector<Job>*, std::vector<Machine>*, float);
    dlerror();
    JobPlaner_t jp = (JobPlaner_t) dlsym(handle, "planJobs");

    dlsym_error = dlerror();
    if(dlsym_error) {
        fprintf(stderr, "Cannot load symbol 'planJobs': %s\n", dlsym_error);
        dlclose(handle);
        return 1;
    }

    printf("Starting random number generator (randomJobSelection = %d) ...\n", env->randomJobSelection);
#ifdef RANDOM
    std::random_device device;
    std::mt19937 generator(device());
#else
    char* seedStr = getenv("RAND_SEED");
    std::mt19937 generator;
    int seed = 1234;
    if(seedStr) {
        seed = atoi(seedStr);
    }
    generator = std::mt19937(seed);
#endif
    std::uniform_int_distribution<int> uni;

    if(env->randomJobSelection) {
        uni = std::uniform_int_distribution<int>(0, jv->size() - 1);
    }

    int fails;
    int active = true;
    float capacity = env->machineSize;

    printf("Running experiments (jobs = %d, ..., %d; step = %d; runs = %d) ...\n", env->startJobs, jmd->numberOfJobs, env->stepWidth, env->runs);
    float avgUsedMachine;

    std::vector<double> runtimeVec = std::vector<double>();
    std::vector<double> successVec = std::vector<double>();

    std::size_t pos = (env->algorithm).find(".");
    std::string algoName = env->algorithm;
    if(pos != std::string::npos) {
        algoName = (env->algorithm).substr(0, pos - 1);
    }
    pos = (env->algorithm).find("LTAA");
    if(pos != std::string::npos) {
        capacity += env->epsilon;
    }
    pos = (env->algorithm).find("Chekuri");
    if(pos != std::string::npos) {
        capacity += env->epsilon;
    }
#ifdef DEBUG
    printf("Capacity is %f\n", capacity);
#endif
    std::stringstream ss;
    ss << algoName;
    ss << "-Jobs:" << env->startJobs << "_" << jmd->numberOfJobs << "_" << env->stepWidth << "_" << env->runs;
    ss << "-Mach:" << env->machines << "-Dims:" << jmd->dimensions;

    std::vector<std::vector < std::string>> machineLoad1 = std::vector<std::vector < std::string >> ();
    std::vector<std::vector < std::string>> machineLoad2 = std::vector<std::vector < std::string >> ();

    //using namespace std::chrono;

    Timer timer = Timer();


#ifdef SAVEJOBS
    char jobFileName[128];
    memset(jobFileName, 0, 128);
    char* jobAlgName = (char*) algoName.c_str();
    while(jobAlgName[0] == '.' || jobAlgName[0] == '/')
        jobAlgName++;
    sprintf(jobFileName, "JobsInProcess-%s-%d.txt", jobAlgName, getpid());
    std::ofstream jobofs(jobFileName);
    if(!jobofs.is_open()) {
        fprintf(stderr, "ERROR: Cannot open file %s\n", jobFileName);
    } else {
        printf("Write to file: %s\n", jobFileName);
    }
#endif

    for(int i = env->startJobs; i <= jmd->numberOfJobs && active; i += env->stepWidth) {
        printf("Experiment with %d jobs and %d machines ...\n", i, env->machines);
        MachineGenerator machineGen(jmd->dimensions);
        std::vector<Machine>* machineQueue = machineGen.generateMachines(env->machines, capacity);
        fails = 0;
        avgUsedMachine = 0.0;
        double averageTime = 0.0;
        int j;
        for(j = 0; j < env->runs && active; j++) {
            for(unsigned i = 0; i < machineQueue->size(); i++) {
                machineQueue->at(i).reset(capacity);
            }
            std::vector<Job>* jq = new std::vector<Job>();
            for(int k = 0; k < i; k++) {
                int pos = env->randomJobSelection ? uni(generator) : k;
                Job tmp = jv->at(pos);
                jq->push_back(tmp);
            }
#ifdef SAVEJOBS
            Job* aTmpJobPtr = NULL;
            float* aTmpKPI = NULL;
            for(unsigned ii = 0; ii < jq->size(); ii++) {
                aTmpJobPtr = &(jq->at(ii));
                aTmpKPI = aTmpJobPtr->getKPIVec();
                for(int jj = 0; jj < aTmpJobPtr->getDimensions(); jj++) {
                    jobofs << aTmpKPI[jj] << ";";
                }
                jobofs << "\n";
            }
            jobofs << "----\n";
            jobofs.flush();
#endif

            std::vector<Machine>* outMachines = 0;
            timer.start();
            jp(fails, &outMachines, jq, machineQueue, env->epsilon);
            timer.stop();
            averageTime += timer.getElapsedTimeInMicroSec();
            avgUsedMachine += outMachines->size();
            machineLoad1.push_back(std::vector <std::string>());
            machineLoad2.push_back(std::vector <std::string>());
            std::vector <std::string>* mload1 = &machineLoad1.back();
            std::vector <std::string>* mload2 = &machineLoad2.back();
            for(int k = 0; k < outMachines->size(); k++) {
                std::stringstream sss1;
                std::stringstream sss2;
                sss1 << outMachines->at(k).getLoad()[0];
                sss2 << outMachines->at(k).assignedJobs->size();
                for(int j = 1; j < outMachines->at(k).getDimensions(); j++) {
                    sss1 << ";" << outMachines->at(k).getLoad()[j];
                }
                mload1->push_back(sss1.str());
                mload2->push_back(sss2.str());
                //std::cout << sss.str() << std::endl;
            }
#ifdef SAVEJOBS
            jobofs << j << " -- " << fails << "\n";
            jobofs.flush();
#endif

            delete jq;
            delete outMachines;
        }
        runtimeVec.push_back(averageTime / (j));
        successVec.push_back((((double) env->runs - fails)) / ((double) j)* 100.0);

        if(fails == env->runs) { // none of the runs was successful
            active = false;
        }

#ifdef DEBUG
        printf("        Machines:\n");
        for(std::vector<Machine>::iterator m = machineQueue->begin(); m != machineQueue->end(); ++m) {
            printf("            ");
            (*m).printLoad(false);
        }
#endif
        avgUsedMachine /= env->runs;
        printf("    Result for %d jobs and %d machines (avg used machines: %3.2f): %3.2f%% success ...\n",
                i, env->machines, avgUsedMachine, (100. - ((float) fails) * 100. / ((float) env->runs)));

        delete machineQueue;

#ifndef RANDOM
        generator = std::mt19937(seed + env->stepWidth);
#endif
    }

#ifdef SAVEJOBS
    jobofs.close();
#endif
    srand48(getpid());
    double shift = (2 * drand48() - 1.0) * 0.2;

    for(unsigned i = 0; i < successVec.size(); i++) {
        successVec[i] += shift;
    }

#ifdef SAVEJOBS
    std::stringstream sss, tss, mss1, mss2;
    int theSeed = -1;
#ifndef RANDOM
    theSeed = seed;
#endif
    sss << "SuccessesTab-" << env->startJobs << "-" << jmd->numberOfJobs << "-" << env->stepWidth << "-S" << theSeed << ".csv";
    tss << "TimingTab-" << env->startJobs << "-" << jmd->numberOfJobs << "-" << env->stepWidth << "-S" << theSeed << ".csv";
    mss1 << "Machineload1-" << env->startJobs << "-" << jmd->numberOfJobs << "-" << env->stepWidth << "-S" << theSeed << ".csv";
    mss2 << "Machineload2-" << env->startJobs << "-" << jmd->numberOfJobs << "-" << env->stepWidth << "-S" << theSeed << ".csv";
//    CSVWriter cw = CSVWriter(ss.str(), successVec, runtimeVec, env->startJobs, jmd->numberOfJobs, env->stepWidth, sss.str(), tss.str(), true);
#endif
    return 0;
}
