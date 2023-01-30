
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

/** Fields of GWF as defined at website (http://gwa.ewi.tudelft.nl/fileadmin/pds/trace-archives/grid-workloads-archive/docs/TheGridWorkloadFormat_v001.pdf). */
struct GWFData {
    long JobID;
    long SubmitTime;
    long WaitTime;
    long RunTime;
    long NProcs; // Number of Allocated Processors
    float AverageCPUTimeUsed;
    float UsedMemory;
    long ReqNProcs; // Requested Number of Processors
    long ReqTime; // Requested Time
    float ReqMemory; // Requested Memory
    long Status;
    char UserID[64];
    char GroupID[64];
    char ExecutableID[256];
    char QueueID[64];
    char PartitionID[64];
    char OrigSiteID[64];
    char LastRunSiteID[64];
    char JobStructure[16];
    char JobStructureParams[64];
    float UsedNetwork;
    float UsedLocalDiskSpace;
    char UsedResources[128];
    char ReqPlatform[64];
    float ReqNetwork;
    float ReqLocalDiskSpace;
    char ReqResources[128];
    long VOID; // Virtual Organization ID
    long ProjectID;
};

class GWFParser {
public:

    GWFParser() {
    }

    std::vector<std::vector<float>*>* parse(const char* _path, Environment* _env, int _num) {
        std::vector<std::vector<float>*>* queue = new std::vector<std::vector<float>*>();

        int count = 0;
        int linecount = 0;
        std::ifstream file(_path);
        if(!file.is_open()) {
            fprintf(stderr, "ERROR: Cannot open job file: %s - exit\n", _path);
            exit(1);
        }
        std::string line;

        float minCPUs = std::numeric_limits<float>::max();
        float maxCPUs = -std::numeric_limits<float>::max();
        float minRAM = std::numeric_limits<float>::max();
        float maxRAM = -std::numeric_limits<float>::max();
        float minLDS = std::numeric_limits<float>::max();
        float maxLDS = -std::numeric_limits<float>::max();
        float minNet = std::numeric_limits<float>::max();
        float maxNet = -std::numeric_limits<float>::max();
        float minTime = std::numeric_limits<float>::max();
        float maxTime = -std::numeric_limits<float>::max();
        float minRCPU = std::numeric_limits<float>::max();
        float maxRCPU = -std::numeric_limits<float>::max();
        float minRRAM = std::numeric_limits<float>::max();
        float maxRRAM = -std::numeric_limits<float>::max();
        float minRLDS = std::numeric_limits<float>::max();
        float maxRLDS = -std::numeric_limits<float>::max();
        float minRNet = std::numeric_limits<float>::max();
        float maxRNet = -std::numeric_limits<float>::max();
        float minRTim = std::numeric_limits<float>::max();
        float maxRTim = -std::numeric_limits<float>::max();
        float val;

        while(std::getline(file, line)) {
            linecount++;
            GWFData in;

            if(isdigit(line[0])) {
                //                     1  2   3   4   5   6  7  8   9   10 11  12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28  29
                sscanf(line.c_str(), "%ld %ld %ld %ld %ld %f %f %ld %ld %f %ld %s %s %s %s %s %s %s %s %s %f %f %s %s %f %f %s %ld %ld", &in.JobID, &in.SubmitTime, &in.WaitTime, &in.RunTime, &in.NProcs, &in.AverageCPUTimeUsed, &in.UsedMemory, &in.ReqNProcs, &in.ReqTime, &in.ReqMemory, &in.Status, in.UserID, in.GroupID, in.ExecutableID, in.QueueID, in.PartitionID, in.OrigSiteID, in.LastRunSiteID, in.JobStructure, in.JobStructureParams, &in.UsedNetwork, &in.UsedLocalDiskSpace, in.UsedResources, in.ReqPlatform, &in.ReqNetwork, &in.ReqLocalDiskSpace, in.ReqResources, &in.VOID, &in.ProjectID);

                std::vector <float>* vec = new std::vector <float>();

                if(in.NProcs >= 0) {
					if(_env->cpus > 0 && in.NProcs > _env->cpus) continue;
                    vec->push_back(((float) in.NProcs) / ((float) _env->cpus));
                    if(in.NProcs < minCPUs) minCPUs = in.NProcs;
                    if(in.NProcs > maxCPUs) maxCPUs = in.NProcs;
                } else {
                    vec->push_back(0.1);
                }
                if(in.UsedMemory >= 0) {
					if(_env->ram > 0 && in.UsedMemory > _env->ram) continue;
                    vec->push_back(((float) in.UsedMemory) / ((float) _env->ram));
                    if(in.UsedMemory < minRAM) minRAM = in.UsedMemory;
                    if(in.UsedMemory > maxRAM) maxRAM = in.UsedMemory;
                } else {
                    vec->push_back(0.1);
                }

                if(in.UsedLocalDiskSpace >= 0) {
					if(_env->maxLocalDisk > 0 && in.UsedLocalDiskSpace > _env->maxLocalDisk) continue;
                    vec->push_back(((float) in.UsedLocalDiskSpace));
                    if(in.UsedLocalDiskSpace < minLDS) minLDS = in.UsedLocalDiskSpace;
                    if(in.UsedLocalDiskSpace > maxLDS) maxLDS = in.UsedLocalDiskSpace;
                } else {
                    vec->push_back(0.1);
                }

                if(in.UsedNetwork >= 0) {
					if(_env->maxNetwork > 0 && in.UsedNetwork > _env->maxNetwork) continue;
                    vec->push_back(((float) in.UsedNetwork));
                    if(in.UsedNetwork < minNet) minNet = in.UsedNetwork;
                    if(in.UsedNetwork > maxNet) maxNet = in.UsedNetwork;
                } else {
                    vec->push_back(0.1);
                }

                if(in.AverageCPUTimeUsed > 0) {
                    vec->push_back(((float) in.AverageCPUTimeUsed));
                    if(in.AverageCPUTimeUsed < minTime) minTime = in.AverageCPUTimeUsed;
                    if(in.AverageCPUTimeUsed > maxTime) maxTime = in.AverageCPUTimeUsed;
                } else {
                    vec->push_back(0.1);
                }


                if(in.ReqNProcs >= 0) {
                    vec->push_back(((float) in.ReqNProcs) / ((float) _env->cpus));
                    if(in.ReqNProcs < minRCPU) minRCPU = in.ReqNProcs;
                    if(in.ReqNProcs > maxRCPU) maxRCPU = in.ReqNProcs;
                } else {
                    vec->push_back(0.1);
                }

                if(in.ReqMemory >= 0) {
                    vec->push_back(((float) in.ReqMemory) / ((float) _env->ram));
                    if(in.ReqMemory < minRRAM) minRRAM = in.ReqMemory;
                    if(in.ReqMemory > maxRRAM) maxRRAM = in.ReqMemory;
                } else {
                    vec->push_back(0.1);
                }

                if(in.ReqLocalDiskSpace >= 0) {
                    vec->push_back(((float) in.ReqLocalDiskSpace));
                    if(in.ReqLocalDiskSpace < minRLDS) minRLDS = in.ReqLocalDiskSpace;
                    if(in.ReqLocalDiskSpace > maxRLDS) maxRLDS = in.ReqLocalDiskSpace;
                } else {
                    vec->push_back(0.1);
                }

                if(in.ReqNetwork >= 0) {
                    vec->push_back(((float) in.ReqNetwork));
                    if(in.ReqNetwork < minRNet) minRNet = in.ReqNetwork;
                    if(in.ReqNetwork > maxRNet) maxRNet = in.ReqNetwork;
                } else {
                    vec->push_back(0.1);
                }

                if(in.ReqTime >= 0) {
                    vec->push_back(((float) in.ReqTime));
                    if(in.ReqTime < minRTim) minRTim = in.ReqTime;
                    if(in.ReqTime > maxRTim) maxRTim = in.ReqTime;
                } else {
                    vec->push_back(0.1);
                }

                queue->push_back(vec);
            }
        }
        file.close();


        float* d;
        for(int i = 0; i < queue->size(); i++) {
            d = queue->at(i)->data();
			if(_env->cpus < 0)
            	if(d[0] > 1.0) d[0] /= maxCPUs;
			else
            	if(d[0] > 1.0) d[0] /= _env->cpus;
			
			if(_env->ram < 0)
            	if(d[1] > 1.0) d[1] /= maxRAM;
			else
            	if(d[1] > 1.0) d[1] /= _env->ram;

			if(_env->maxLocalDisk < 0)
            	if(d[2] > 1.0) d[2] /= maxLDS;
            else
				if(d[2] > 1.0) d[2] /= _env->maxLocalDisk;

			if(_env->maxNetwork < 0)
            	if(d[3] > 1.0) d[3] /= maxNet;
			else
            	if(d[3] > 1.0) d[3] /= _env->maxNetwork;

            if(d[4] > 1.0) d[4] /= maxTime;
            if(d[5] > 1.0) d[5] /= maxRCPU;
            if(d[6] > 1.0) d[6] /= maxRRAM;
            if(d[7] > 1.0) d[7] /= maxRLDS;
            if(d[8] > 1.0) d[8] /= maxRNet;
            if(d[9] > 1.0) d[9] /= maxRTim;
            /*if(d[0] > 1.0) d[0] = 1.0;
            if(d[1] > 1.0) d[1] = 1.0;
            if(d[2] > 1.0) d[2] = 1.0;
            if(d[3] > 1.0) d[3] = 1.0;
            if(d[4] > 1.0) d[4] = 1.0;
            if(d[5] > 1.0) d[5] = 1.0;
            if(d[6] > 1.0) d[6] = 1.0;
            if(d[7] > 1.0) d[7] = 1.0;
            if(d[8] > 1.0) d[8] = 1.0;
            if(d[9] > 1.0) d[9] = 1.0;*/
        }

        printf("min CPU: %f max CPU: %f min RAM: %f max RAM: %f \n", minCPUs, maxCPUs, minRAM, maxRAM);
        return queue;
    }


};

extern "C"
std::vector <Job>* loadJobs(const char* _path, JobMetaData* _jobdata, Environment* _env, int _runs) {
    std::vector <Job>* jobVec = new std::vector <Job>();

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

    GWFParser parser = GWFParser();

    // generate factor* m_num vectors from swf file
    std::vector < std::vector < float >* >* vectors = parser.parse(_path, _env, _runs * _jobdata->numberOfJobs);
    if(_jobdata->dimensions > vectors->at(0)->size()) {
        std::vector < float >* tmp;
        for(int i = 0; i < vectors->size(); i++) {
            tmp = vectors->at(i);
            int pos = 0;
            while(tmp->size() < _jobdata->dimensions) {
                tmp->push_back(tmp->at(pos));
                pos++;
            }
        }
    }

    if(_jobdata->dimensions < vectors->at(0)->size()) {
        for(int i = 0; i < vectors->size(); i++) {
            vectors->at(i)->resize(_jobdata->dimensions);
        }
    }

    unsigned dims = vectors->at(0)->size();
    unsigned elem = vectors->size();
    printf("vectors size1: %lu\n", vectors->size());
    printf("vectors size2: %lu\n", vectors->at(0)->size());

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

