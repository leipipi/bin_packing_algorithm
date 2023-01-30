#ifndef PANI_H_
#define PANI_H_

#include "Job.h"
#include "Machine.h"
#include "MathFunctions.h"
#include "Sorting.h"
#include <vector>
#include <cmath>
#include <list>
#include <limits.h>
#include <stdlib.h> 

class Pani {
public:

    static const int LINF = 0;
    static const int L1 = 1;
    static const int L2 = 2;
    static const int DP = 3;

    static char* getName(int id) {
        switch(id) {
            case LINF:
                return strdup("LINF");
            case L1:
                return strdup("L1");
            case L2:
                return strdup("L2");
            case DP:
                return strdup("DP");
            default:
                return strdup("Error");
        }
    }
    
    struct pair {
        int machine;
        int job;
        float value;
    };

    static float* averageDemandExp(std::vector<Job>* _jobs, float _factor) {
        // Stop if there are no jobs
        if(_jobs->empty()) {
            return 0;
        }
        // Create float array vec
        int dimensions = _jobs->at(0).getDimensions();
        float* vec = new float[dimensions];
        // Initialize all coordinates in vec with 0
        for(unsigned d = 0; d < dimensions; d++) {
            vec[d] = 0;
        }
        // Add all job vectors to vec
        for(unsigned i = 0; i < _jobs->size(); i++) {
            MathFunctions::VectorAddition(vec, _jobs->at(i).getKPIVec(), dimensions);
        }
        // Divide vec by number of jobs to get average
        MathFunctions::VecDiv(vec, static_cast<float> (_jobs->size()), dimensions);
        // Replace i-th coordinate in vec by exp(eps * vec[i]) and return result
        for(int d = 0; d < dimensions; d++) {
            vec[d] = std::exp(_factor * vec[d]);
        }
        return vec;
    }

 static float applyNormPercentage(float* avgDemandExp, float* difference, float* jv, float* remCap, float* linfInOut, int dimensions, const int type) {
        if(type == DP) {
            return -MathFunctions::DotProd(avgDemandExp, jv, remCap, dimensions);
        }
        // Compute difference of remaining capacity and job vector
        for(int k = 0; k < dimensions; k++) {
            difference[k] = remCap[k] - jv[k];
        }
        float percentage = 0.0;
        switch(type) {
            case L1:
            {
                percentage = MathFunctions::DotProd(avgDemandExp, difference, dimensions)/MathFunctions::DotProd(avgDemandExp, remCap, dimensions);
           //    printf("percentage is %f", percentage);
                return percentage;
            }
            case L2:
            {
                percentage = MathFunctions::DotProd(avgDemandExp, difference, difference, dimensions)/MathFunctions::DotProd(avgDemandExp, remCap, remCap, dimensions);
          //      printf("percentage is %f", percentage);
                return percentage;
            }
            case LINF:
            {
                percentage = MathFunctions::Linf(MathFunctions::VecCoordProd(linfInOut, avgDemandExp, difference, dimensions), dimensions)/MathFunctions::Linf(MathFunctions::VecCoordProd(linfInOut, avgDemandExp, remCap, dimensions), dimensions);
           //     printf("percentage is %f", percentage);
                return percentage;
            }
            default:
                printf("ERROR: type needs to be specified!\n");
                return FLT_MIN;
        }
    }

    static float applyNorm(float* avgDemandExp, float* difference, float* jv, float* remCap, float* linfInOut, int dimensions, const int type) {
        if(type == DP) {
            return -MathFunctions::DotProd(avgDemandExp, jv, remCap, dimensions);
        }
        // Compute difference of remaining capacity and job vector
        for(int k = 0; k < dimensions; k++) {
            difference[k] = remCap[k] - jv[k];
        }

        switch(type) {
            case L1:
                return MathFunctions::DotProd(avgDemandExp, difference, dimensions);
            case L2:
                return MathFunctions::DotProd(avgDemandExp, difference, difference, dimensions);
            case LINF:
                return MathFunctions::Linf(MathFunctions::VecCoordProd(linfInOut, avgDemandExp, difference, dimensions), dimensions);
            default:
                printf("ERROR: type needs to be specified!\n");
                return FLT_MIN;
        }
    }
    
    struct CompareFuncWeight {
        bool operator()(Job const &a, Job const &b) {
            return a.getWeight() < b.getWeight();
        }
    };

    static void updateBestJobs(std::vector<Machine>* mq, int indexOfLastJobAssigned) {
        for(int i = 0; i < mq->size(); i++) {
            int bestJobIndex = mq->at(i).getBestJobIndex();
            if(bestJobIndex >= indexOfLastJobAssigned) {
                if(bestJobIndex > indexOfLastJobAssigned) {
                    mq->at(i).setBestJobIndex(bestJobIndex - 1);
                } else {
                    mq->at(i).setBestJobIndex(-1);
                }
            }
        }
    }


    static bool runPanigrahy(int& counter, std::vector<Machine>** out_machines, std::vector<Job>* _jq, std::vector<Machine>* _mq, float _epsilon, const int type) {
        // Make a copy of _jq (not sure if this actually achieves anything)
        std::vector<Job*>* aJQ = new std::vector<Job*>();
        for(int i = 0; i < _jq->size(); i++) {
            aJQ->push_back(&(_jq->at(i)));
        }

        // Create machine vector with first machine
        int dimensions = aJQ->at(0)->getDimensions();
        MachineGenerator mach_gen(dimensions);
        std::vector<Machine>* mq = mach_gen.generateMachines(1, 1.0);
        float dp;
        float* tmp_jv;
        float* tmp_mv;
        // Create Panigrahy weight vector 
        float* avgDemandExp = averageDemandExp(_jq, _epsilon);
        float* difference = new float[dimensions];
        float* linfInOut = new float[dimensions];

        while(!aJQ->empty()) {
            int job_index = -1;
            float min = FLT_MAX;
            tmp_mv = mq->back().getRemainingCapacity();

            for(unsigned i = 0; i < aJQ->size(); i++) {
                tmp_jv = aJQ->at(i)->getKPIVec();
                if(!MathFunctions::isVectorDifferencePositive(tmp_mv, tmp_jv, dimensions)) {
                    continue;
                }
                // Assess job by norm (LINF, L1, L2) or dot product (DP)
                dp = applyNorm(avgDemandExp, difference, tmp_jv, tmp_mv, linfInOut, dimensions, type);
                // If job is better than current minimum, update job_index 
                if(dp < min) {
                    job_index = i;
                    min = dp;
                }
            }

            if(job_index != -1) {
                mq->back().assignJob(aJQ->at(job_index));
                aJQ->erase(aJQ->begin() + job_index);
            } else {
                mq->push_back(mach_gen.generateMachine(mq->size() + 1, 1.0));
            }
        }
        delete aJQ;
        delete [] avgDemandExp;
        delete [] difference;
        delete [] linfInOut;
        *out_machines = mq;
        if(mq->size() > _mq->size()) {
            return false;
        }
        return true;
    }

    /**
     * DONE: If best job for bin is not assigned in round, it will be the best job for that bin in the next round again!
     * DONE: Machines are marked as full; it would be better to remove them from the vector.
     * TODO: If jobs were sorted by "size" (which depends on norm), the first job that cannot be assigned to any active bin can immediately be assigned to new (empty) bin
     * TODO: A job could remember in how many bins it fits; when a bin is added / altered, this information needs to be updated for every active job
     * 
     * TODO: Alternative: maintain an m x n bit matrix storing for every job-machine pair if the job can be placed on the machine; bitsets would allow quick checks if a job can be assigned to any machine.
     * In addition store a matrix with the weight of assigning job i to machine j. However, there are complications, and it might not be faster.
     */
    static bool runHybrid(int& counter, std::vector<Machine>** out_machines, std::vector<Job>* _jq, std::vector<Machine>* _mq, float _epsilon, int type) {
        // Make a copy of _jq (not sure if this actually achieves anything)
        std::vector<Job*>* aJQ = new std::vector<Job*>();
        for(int i = 0; i < _jq->size(); i++) {
            aJQ->push_back(&(_jq->at(i)));
        }
        // Create machine vector with first machine
        int dimensions = aJQ->at(0)->getDimensions();
        MachineGenerator mach_gen(dimensions);
        std::vector<Machine>* mq = mach_gen.generateMachines(1, 1.0);
        std::vector<Machine>* mqFull = mach_gen.generateMachines(0, 1.0);
        float dp;
        int indexOfLastJobAssigned = -1;
        float* tmp_jv;
        float* tmp_mv = mq->back().getRemainingCapacity();
        // Create Panigrahy weight vector 
        float* avgDemandExp = averageDemandExp(_jq, _epsilon);
        float* difference = new float[dimensions];
        float* linfInOut = new float[dimensions]; // only for LINF
        float* remCapOfEmptyMachine = new float[dimensions];
        for(int i = 0; i < dimensions; i++) {
            remCapOfEmptyMachine[i] = tmp_mv[i];
        }
        // Assign all jobs in a loop
        while(!aJQ->empty()) {
            int job_index = -1;
            float min = FLT_MAX;

            for(int j = 0; j < aJQ->size(); j++) {
                tmp_jv = aJQ->at(j)->getKPIVec();
                bool couldBeAssigned = false;
                for(int i = 0; i < mq->size(); i++) {
                    tmp_mv = mq->at(i).getRemainingCapacity();
                    if(MathFunctions::isVectorDifferencePositive(tmp_mv, tmp_jv, dimensions)) {
                        couldBeAssigned = true;
                        break;
                    }
                }
                // If jobs do not fit into any machine, the best of them will be placed in a new machine
                if(!couldBeAssigned) {
                    tmp_jv = aJQ->at(j)->getKPIVec();
                    dp = applyNorm(avgDemandExp, difference, tmp_jv, remCapOfEmptyMachine, linfInOut, dimensions, type);
                    if(dp < min) {
                        job_index = j;
                        min = dp;
                    }
                }
            }
            // If at least one job could not be fit into any machine, place the largest of them in a new machine
            if(job_index != -1) {
                // Create new machine
                mq->push_back(mach_gen.generateMachine(mq->size() + 1, 1.0));
                // Assign best job to new bin and remove it from job queue
                mq->back().setBestJobIndex(job_index);
                mq->back().assignJob(aJQ->at(job_index));
                aJQ->erase(aJQ->begin() + job_index);
                indexOfLastJobAssigned = job_index;
            } else {
                // Loop over all machines
                for(int i = (mq->size())-1; i >= 0; i--) {
                    if(mq->at(i).getBestJobIndex() != -1) {
                        continue;
                    }
                    //printf("XXXXXX machine: %d, jq size: %ld", i, aJQ->size());
                    int job_index = -1;
                    float min = FLT_MAX;
                    // Loop over all jobs
                    for(int j = 0; j < aJQ->size(); j++) {
                        tmp_jv = aJQ->at(j)->getKPIVec();
                        tmp_mv = mq->at(i).getRemainingCapacity();
                        // If job does not fit, go to next one
                        if(!MathFunctions::isVectorDifferencePositive(tmp_mv, tmp_jv, dimensions)) {
                            continue;
                        }
                        dp = applyNorm(avgDemandExp, difference, tmp_jv, tmp_mv, linfInOut, dimensions, type);

                        if(dp < min) {
                            job_index = j;
                            min = dp;
                        }
                    }
                    if(job_index == -1) {
                        //printf(", machine removed\n");
                        mqFull->push_back(mq->at(i));
                        mq->erase(mq->begin() + i); 
                    } else {
                        //printf(", best job: %d\n", job_index);
                        mq->at(i).setBestJobIndex(job_index);
//                        if(type == DP) {
//                            // Best machine is determined by L2 here; one could try L1 or LINF
//                            tmp_jv = aJQ->at(job_index)->getKPIVec();
//                            tmp_mv = mq->at(i).getRemainingCapacity();
//                            for(int k = 0; k < dimensions; k++) {
//                                difference[k] = tmp_mv[k] - tmp_jv[k];
//                            }
//                            mq->at(i).setValueOfBestJob(MathFunctions::DotProd(avgDemandExp, difference, difference, dimensions));
//                        } else {
                            mq->at(i).setValueOfBestJob(min);
//                        }
                    }
                }
                
                if(type == DP) {
                    // Idea: 
                    // (1) If two machines want the same job, we prefer the machine with the smaller dot product because it is a better fit
                    // (2) If two machines want different jobs, we prefer the machine with the larger dot product because we want to fill lower-loaded bins first
                    // This way we hopefully create more fuller bins among which the later jobs can "choose"
                    // Implementation:
                    // Apply rule 1 and remove machines, until each machine in the list wants a different job
                    // Then apply rule 2 and pick the best pair
                    // List of job-machine pairs; for each job, the best machine is the one that minimises dot product (or maximises negative dot product)
                    std::list<pair> pairs;
                    // Go through list of machines, retrieve best job for it and add pair to list, unless there is already a better pair with same machine
                    for(int i = 0; i < mq->size(); i++) {
                        int jobIndex = mq->at(i).getBestJobIndex();
                        float value = mq->at(i).getValueOfBestJob();
                        bool found = false;
                        for(pair p : pairs) { // TODO: could be done more efficiently by keeping list sorted by job index
                            if(p.job == jobIndex) {
                                if(p.value < value) {
                                    p.machine = i;
                                    p.job = jobIndex;
                                    p.value = value;
                                }
                                found = true;
                                break;
                            }
                        }
                        if(!found) {
                            pair p;
                            p.machine = i;
                            p.job = jobIndex;
                            p.value = value;
                            pairs.push_back(p);
                        }
                    }
                    // Find best pair; now the best pair is the one that maximises the dot product
                    pair bestPair;
                    bestPair.value = 0.;
                    for(pair p : pairs) {
                        if(p.value < bestPair.value) {
                            bestPair = p;
                        }
                    }
                    // Assign job to machine
                    int bestMachine = bestPair.machine;
                    int chosenJobIndex = bestPair.job;
                    mq->at(bestMachine).assignJob(aJQ->at(chosenJobIndex));
                    aJQ->erase(aJQ->begin() + chosenJobIndex);
                    indexOfLastJobAssigned = chosenJobIndex;

                } else {
                    int bestMachine = -1;
                    float bestValue = FLT_MAX;
                    for(int i = 0; i < mq->size(); i++) {
                        if(mq->at(i).getValueOfBestJob() < bestValue) {
                            bestValue = mq->at(i).getValueOfBestJob();
                            bestMachine = i;
                        }
                    }
                    int chosenJobIndex = mq->at(bestMachine).getBestJobIndex();
                    mq->at(bestMachine).assignJob(aJQ->at(chosenJobIndex));
                    aJQ->erase(aJQ->begin() + chosenJobIndex);
                    indexOfLastJobAssigned = chosenJobIndex;
                }
            }

            for(int i = 0; i < mq->size(); i++) {
                int bestJobIndex = mq->at(i).getBestJobIndex();
                if(bestJobIndex >= indexOfLastJobAssigned) {
                    if(bestJobIndex > indexOfLastJobAssigned) { 
                        mq->at(i).setBestJobIndex(bestJobIndex - 1);
                    } else {
                        mq->at(i).setBestJobIndex(-1);
                    }
                }
            }
        }
        // Clean up and return success / failure
        bool success = mq->size() + mqFull->size() <= _mq->size();
        for(int i = 0; i < mq->size(); i++) {
            mqFull->push_back(mq->at(i));
        }
        *out_machines = mqFull;
        delete mq;
        delete aJQ;
        delete [] avgDemandExp;
        delete [] difference;
        delete [] linfInOut;
        delete [] remCapOfEmptyMachine;
        return success;
    }
    static bool runHybrid2(int& counter, std::vector<Machine>** out_machines, std::vector<Job>* _jq, std::vector<Machine>* _mq, float _epsilon, int type) {
        // Make a copy of _jq (not sure if this actually achieves anything)
        std::vector<Job*>* aJQ = new std::vector<Job*>();
        for(int i = 0; i < _jq->size(); i++) {
            aJQ->push_back(&(_jq->at(i)));
        }
        // Create machine vector with first machine
        int dimensions = aJQ->at(0)->getDimensions();
        MachineGenerator mach_gen(dimensions);
        std::vector<Machine>* mq = mach_gen.generateMachines(1, 1.0);
        std::vector<Machine>* mqFull = mach_gen.generateMachines(0, 1.0);
        float dp;
        int indexOfLastJobAssigned = -1;
        float* tmp_jv;
        
        
        
        
        
        
     
        //   MachineGenerator mach_gen(dimensions);
        
        float* total = new float[dimensions];
        
        //analyze standard deviation:
        
        for(unsigned j = 0; j < dimensions; j++)
        {
            total[j] = 0;
        }
        for(int i = 0; i < _jq->size(); i++) {
            
            //  printf("Job id is %d\n",_jq->at(i).getID());
            tmp_jv = _jq->at(i).getKPIVec();
            for(unsigned j = 0; j < dimensions; j++)
            {
                total[j] += tmp_jv[j];
            }
            
        }
    //    printf("total 0 is %f\n",total[0]);
        for(unsigned i = 0; i < dimensions; i++){
            if(total[i]>=_mq->size()){
           //     counter++;
       //         printf("xzz\n");
                *out_machines = mq;
                return false;
            }
        }
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        float* tmp_mv = mq->back().getRemainingCapacity();
        // Create Panigrahy weight vector 
        float* avgDemandExp = averageDemandExp(_jq, _epsilon);
        float* difference = new float[dimensions];
        float* linfInOut = new float[dimensions]; // only for LINF
        float* remCapOfEmptyMachine = new float[dimensions];
        for(int i = 0; i < dimensions; i++) {
            remCapOfEmptyMachine[i] = tmp_mv[i];
        }

        // Assign all jobs in a loop
        while(!aJQ->empty()) {
            int job_index = -1;
            float min = FLT_MAX;
            
            //check whether a new bin should be open
            for(int j = 0; j < aJQ->size(); j++) {
                tmp_jv = aJQ->at(j)->getKPIVec();
                bool couldBeAssigned = false;
                for(int i = 0; i < mq->size(); i++) {
                    tmp_mv = mq->at(i).getRemainingCapacity();
                    if(MathFunctions::isVectorDifferencePositive(tmp_mv, tmp_jv, dimensions)) {
                        couldBeAssigned = true;
                        break;
                    }
                }
                
                // If jobs do not fit into any machine, the best of them will be placed in a new machine
                if(!couldBeAssigned) {
                    tmp_jv = aJQ->at(j)->getKPIVec();
                    dp = applyNorm(avgDemandExp, difference, tmp_jv, remCapOfEmptyMachine, linfInOut, dimensions, type);
                    if(dp < min) {
                        job_index = j;
                        min = dp;
                    }
                }
            }
            // If at least one job could not be fit into any machine, place the largest of them in a new machine
            if(job_index != -1) {
                // Create new machine
                mq->push_back(mach_gen.generateMachine(mq->size() + 1, 1.0));
                // Assign best job to new bin and remove it from job queue
                mq->back().setBestJobIndex(job_index);
                mq->back().assignJob(aJQ->at(job_index));
                aJQ->erase(aJQ->begin() + job_index);
                indexOfLastJobAssigned = job_index;
            } else {
                // Loop over all machines to send their proposal and
                for(int i = (mq->size())-1; i >= 0; i--) {
                    if(mq->at(i).getBestJobIndex() != -1) {
                        continue;
                    }
                    int job_index = -1;
                    float min = FLT_MAX;
                    // Loop over all jobs
                    for(int j = 0; j < aJQ->size(); j++) {
                        tmp_jv = aJQ->at(j)->getKPIVec();
                        tmp_mv = mq->at(i).getRemainingCapacity();
                        // If job does not fit, go to next one
                        if(!MathFunctions::isVectorDifferencePositive(tmp_mv, tmp_jv, dimensions)) {
                            continue;
                        }
                     //   printf("\n ajq size is %ld ",aJQ->size());
                        //righthere
                        dp = applyNormPercentage(avgDemandExp, difference, tmp_jv, tmp_mv, linfInOut, dimensions, type);
                        // use the percentage of the dp compared with the total remaining capacity as a parameter.
                        if(dp < min) {
                            job_index = j;
                            min = dp;
                        }
                    }
                    if(job_index == -1) {
                        //full machine
                        mqFull->push_back(mq->at(i));
                        mq->erase(mq->begin() + i); 
                    } else {
                        mq->at(i).setBestJobIndex(job_index);
                        if(type == DP) {
                            // Best machine is determined by L2 here; one could try L1 or LINF
                            tmp_jv = aJQ->at(job_index)->getKPIVec();
                            tmp_mv = mq->at(i).getRemainingCapacity();
                            for(int k = 0; k < dimensions; k++) {
                                difference[k] = tmp_mv[k] - tmp_jv[k];
                            }
                            mq->at(i).setValueOfBestJob(MathFunctions::DotProd(avgDemandExp, difference, difference, dimensions));
                        } else {
                            mq->at(i).setValueOfBestJob(min);
                        }
                    }
                }
                int bestMachine = -1;
                float bestValue = FLT_MAX;
                for(int i = 0; i < mq->size(); i++) {
                    if(mq->at(i).getValueOfBestJob() < bestValue) {
                        bestValue = mq->at(i).getValueOfBestJob();
                        bestMachine = i;
                    }
                }
                //get the best proposal and assign
                int chosenJobIndex = mq->at(bestMachine).getBestJobIndex();
                mq->at(bestMachine).assignJob(aJQ->at(chosenJobIndex));
                aJQ->erase(aJQ->begin() + chosenJobIndex);
                indexOfLastJobAssigned = chosenJobIndex;
            }
//update bestjob in the machines
            for(int i = 0; i < mq->size(); i++) {
                int bestJobIndex = mq->at(i).getBestJobIndex();
                if(bestJobIndex >= indexOfLastJobAssigned) {
                    if(bestJobIndex > indexOfLastJobAssigned) { 
                        mq->at(i).setBestJobIndex(bestJobIndex - 1);
                    } else {
                        mq->at(i).setBestJobIndex(-1);
                    }
                }
            }
        }
        // Clean up and return success / failure
        bool success = mq->size() + mqFull->size() <= _mq->size();
        for(int i = 0; i < mq->size(); i++) {
            mqFull->push_back(mq->at(i));
        }
        *out_machines = mqFull;
        delete mq;
        delete aJQ;
        delete [] avgDemandExp;
        delete [] difference;
        delete [] linfInOut;
        delete [] remCapOfEmptyMachine;
        return success;
    }
    /**
     * This algorithm assumes that all coordinates of all bins have the same capacity. 
     * This algorithm assumes that no job on its own exceeds the machine capacity.
     * 
     * 
     * DONE: If best job for bin is not assigned in round, it will be the best job for that bin in the next round again!
     * DONE: Machines are marked as full; it would be better to remove them from the vector.
     * TODO: If jobs were sorted by "size" (which depends on norm), the first job that cannot be assigned to any active bin can immediately be assigned to new (empty) bin
     * TODO: A job could remember in how many bins it fits; when a bin is added / altered, this information needs to be updated for every active job
     * 
     * TODO: Alternative: maintain an m x n bit matrix storing for every job-machine pair if the job can be placed on the machine; bitsets would allow quick checks if a job can be assigned to any machine.
     * In addition store a matrix with the weight of assigning job i to machine j. However, there are complications, and it might not be faster.
     */
    static bool runAlternativeHybrid(int& counter, std::vector<Machine>** out_machines, std::vector<Job>* _jq, std::vector<Machine>* _mq, float _epsilon, int type) {
        // Create machine vector with first machine
        int dimensions = _jq->at(0).getDimensions();
        MachineGenerator mach_gen(dimensions);
        std::vector<Machine>* mq = mach_gen.generateMachines(1, 1.0);
        std::vector<Machine>* mqFull = mach_gen.generateMachines(0, 1.0);
        float dp;
        float* tmp_jv;
        float* tmp_mv = mq->back().getRemainingCapacity();
        // Create Panigrahy weight vector 
        float* avgDemandExp = averageDemandExp(_jq, _epsilon);
        float* difference = new float[dimensions];
        float* linfInOut = new float[dimensions]; // only for LINF
        float* remCapOfEmptyMachine = new float[dimensions];
        for(int i = 0; i < dimensions; i++) {
            remCapOfEmptyMachine[i] = tmp_mv[i];
        }

        // Make a copy of _jq and set weight and initialise counter for each job
        std::vector<Job>* aJQ = new std::vector<Job>();
        for(int i = 0; i < _jq->size(); i++) {
            aJQ->push_back(_jq->at(i));
            // The weight of a job is the value one gets applying the norm / dot product to the job and an empty machine
            dp = applyNorm(avgDemandExp, difference, aJQ->back().getKPIVec(), remCapOfEmptyMachine, linfInOut, dimensions, type);
            aJQ->back().setWeight(dp);
            // Counter counts number of machines the job can be fitted in; since there is one empty machine at the start, it is initiated with 1
            aJQ->back().setCounter(1);
        }
        // Sort jobs by weight 
        std::sort(aJQ->begin(), aJQ->end(), CompareFuncWeight());


        // Assign all jobs in a loop
        while(!aJQ->empty()) {
            // Loop over all machines
            for(int i = (mq->size()) - 1; i >= 0; i--) {
                if(mq->at(i).getBestJobIndex() != -1) {
                    continue;
                }
                int job_index = -1;
                float min = FLT_MAX;
                // Loop over all jobs
                for(int j = 0; j < aJQ->size(); j++) {
                    tmp_jv = aJQ->at(j).getKPIVec();
                    tmp_mv = mq->at(i).getRemainingCapacity();
                    // If job does not fit, go to next one
                    if(!MathFunctions::isVectorDifferencePositive(tmp_mv, tmp_jv, dimensions)) {
                        continue;
                    }
                    dp = applyNorm(avgDemandExp, difference, tmp_jv, tmp_mv, linfInOut, dimensions, type);

                    if(dp < min) {
                        job_index = j;
                        min = dp;
                    }
                }
                if(job_index == -1) {
                    mqFull->push_back(mq->at(i));
                    mq->erase(mq->begin() + i);
                } else {
                    mq->at(i).setBestJobIndex(job_index);
                    if(type == DP) {
                        // Best machine is determined by L2 here; one could try L1 or LINF
                        tmp_jv = aJQ->at(job_index).getKPIVec();
                        tmp_mv = mq->at(i).getRemainingCapacity();
                        for(int k = 0; k < dimensions; k++) {
                            difference[k] = tmp_mv[k] - tmp_jv[k];
                        }
                        mq->at(i).setValueOfBestJob(MathFunctions::DotProd(avgDemandExp, difference, difference, dimensions));
                    } else {
                        mq->at(i).setValueOfBestJob(min);
                    }
                }
            }
            int bestMachine = -1;
            float bestValue = FLT_MAX;
            for(int i = 0; i < mq->size(); i++) {
                if(mq->at(i).getValueOfBestJob() < bestValue) {
                    bestValue = mq->at(i).getValueOfBestJob();
                    bestMachine = i;
                }
            }
            int chosenJobIndex = mq->at(bestMachine).getBestJobIndex();
            // Save chosen job vector for updating counters of other jobs
            MathFunctions::copyVector(mq->at(bestMachine).getRemainingCapacity(), linfInOut, dimensions);
            // Assign chosen job to chosen machine and remove it from job list
            mq->at(bestMachine).assignJob(&(aJQ->at(chosenJobIndex)));
            aJQ->erase(aJQ->begin() + chosenJobIndex);
            updateBestJobs(mq, chosenJobIndex);

            tmp_mv = mq->at(bestMachine).getRemainingCapacity();
            int indexOfFirstZeroJob = -1;
            for(int j = 0; j < aJQ->size(); j++) {
                // If job j cannot be be fitted into machine, but could prior to assignment, then decrement its counter
                if(!MathFunctions::isVectorDifferencePositive(tmp_mv, aJQ->at(j).getKPIVec(), dimensions)) {
                    if(MathFunctions::isVectorDifferencePositive(linfInOut, aJQ->at(j).getKPIVec(), dimensions)) {
                        aJQ->at(j).decrementCounter();
                        if(indexOfFirstZeroJob == -1 && aJQ->at(j).isCounterZero()) {
                            indexOfFirstZeroJob = j;
                        }
                    }
                }
            }
            while(indexOfFirstZeroJob != -1) {
                // Create new machine
                mq->push_back(mach_gen.generateMachine(mq->size() + 1, 1.0));
                // Assign best job to new bin and remove it from job queue
                mq->back().setBestJobIndex(indexOfFirstZeroJob);
                mq->back().assignJob(&(aJQ->at(indexOfFirstZeroJob)));
                aJQ->erase(aJQ->begin() + indexOfFirstZeroJob);
                updateBestJobs(mq, indexOfFirstZeroJob);
                indexOfFirstZeroJob = -1;
                tmp_mv = mq->back().getRemainingCapacity();
                for(int j = 0; j < aJQ->size(); j++) {
                    if(MathFunctions::isVectorDifferencePositive(tmp_mv, aJQ->at(j).getKPIVec(), dimensions)) {
                        aJQ->at(j).incrementCounter();
                    } else if(indexOfFirstZeroJob == -1 && aJQ->at(j).isCounterZero()) {
                        indexOfFirstZeroJob = j;
                    }
                }
            }
        }
        // Clean up and return success / failure
        bool success = mq->size() + mqFull->size() <= _mq->size();
        for(int i = 0; i < mq->size(); i++) {
            mqFull->push_back(mq->at(i));
        }
        *out_machines = mqFull;
        delete mq;
        delete aJQ;
        delete [] avgDemandExp;
        delete [] difference;
        delete [] linfInOut;
        delete [] remCapOfEmptyMachine;
        return success;
    }



};

#endif
