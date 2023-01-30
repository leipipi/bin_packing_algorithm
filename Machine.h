#ifndef MACHINE_H_
#define MACHINE_H_

#include <vector>
#include "Job.h"

class Machine {
public:
    Machine(const Machine& _m);
    Machine(int _id, int _dimensions, float _capacity);
    Machine& operator=(const Machine& _m);
    virtual ~Machine();

    float* getLoad();
    float getSumOfLoad();
    int getDimensions();
    void reset(float _cap);
    bool assignJob(Job* _job);
    bool assignJobLTAA(Job* _job, bool small);
    bool doesRoundedVectorFit(const float* roundedVector);
    bool doesVectorFit(const float* vector);
    void assignJobOverload(Job* _job);
    // TODO: We should get rid of setVector(...) which is used only once
    void setVector(const std::vector<float> vec);
    float* getRemainingCapacity();
    std::vector<Job>* getSmallJobs();
    void removeSmallJobs();
    //assigned ID is the position of job in the assignedJobs vector
    void removeJob(int assignedID);
    Job popBackJobInit();
    Job popJob(int index);
    bool isFlaggedAsFull();
    void flagAsFull();
    void setBestJobIndex(int index);
    int getBestJobIndex();
    void setValueOfBestJob(float value);
    float getValueOfBestJob();
    void printLoad(bool printRemCap);
    void printMachine();
    float* getJobKPI(int index);
    int getSizeOfAssignedJobs();
    // TODO: assignedJobs should be private, too
    std::vector<Job>* assignedJobs;
    
private:
    int id;
    int dimensions;    
    float capacity;
    float* load;
    float* roundedLoad;
    float* remainingCapacity;
    float* roundedRemainingCapacity;
    std::vector<Job>* smallJobs;
    bool flaggedAsFull;
    int bestJobIndex;
    float valueOfBestJob;
};



class MachineGenerator {
public:
    MachineGenerator(int _dims);

    std::vector<Machine>* generateMachines(int _numberOfMachines, float _capacity);
    Machine generateMachine(int _id, float _capacity);

private:
    int dimensions;
};

#endif
