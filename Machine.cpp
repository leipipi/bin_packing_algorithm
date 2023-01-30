#include "Machine.h"
#include "MathFunctions.h"
#include <cstring>


Machine::Machine(const Machine& _m) {
    id = _m.id;
    dimensions = _m.dimensions;
    capacity = _m.capacity;
    assignedJobs = new std::vector<Job>();
    smallJobs = new std::vector<Job>();
    load = new float[dimensions];
    roundedLoad = new float[dimensions];
    remainingCapacity = new float[dimensions];
    roundedRemainingCapacity = new float[dimensions];
    memcpy(load, _m.load, dimensions * sizeof(float));
    memcpy(roundedLoad, _m.roundedLoad, dimensions * sizeof(float));
    memcpy(remainingCapacity, _m.remainingCapacity, dimensions * sizeof(float));
    memcpy(roundedRemainingCapacity, _m.roundedRemainingCapacity, dimensions * sizeof(float));
    for(unsigned i=0; i<_m.assignedJobs->size(); i++) {
        assignedJobs->push_back(_m.assignedJobs->at(i));
    }
    for(unsigned i=0; i<_m.smallJobs->size(); i++) {
        smallJobs->push_back(_m.smallJobs->at(i));
    }
    flaggedAsFull = _m.flaggedAsFull;
    bestJobIndex = _m.bestJobIndex;
    valueOfBestJob = _m.valueOfBestJob;
}

Machine::Machine(int _id, int _dimensions, float _capacity) {
    id = _id;
    dimensions = _dimensions;
    capacity = _capacity;
    load = new float[dimensions];
    roundedLoad = new float[dimensions];
    remainingCapacity = new float[dimensions];
    roundedRemainingCapacity = new float[dimensions];
    for(int i=0; i<dimensions; i++) {
        load[i] = 0.0;
        roundedLoad[i] = 0.0;
        remainingCapacity[i] = capacity;
        roundedRemainingCapacity[i] = capacity;
        remainingCapacity[i] = capacity;
    }
    assignedJobs = new std::vector<Job>();
    smallJobs = new std::vector<Job>();
    flaggedAsFull = false;
    bestJobIndex = -1;
    valueOfBestJob = -1.;
}

Machine& Machine::operator=(const Machine& _m) {
    if(this == &_m)
        return *this;
    if(dimensions)
    {
        delete [] load;
        delete [] roundedLoad;
        delete [] remainingCapacity;
        delete [] roundedRemainingCapacity;
        delete assignedJobs;
        delete smallJobs;
    }
    id = _m.id;
    dimensions = _m.dimensions;
    capacity = _m.capacity;
    load = 0;
    roundedLoad = 0;
    remainingCapacity = 0;
    roundedRemainingCapacity = 0;
    assignedJobs = new std::vector <Job>();
    smallJobs = new std::vector <Job>();
    if(_m.dimensions) {
        load = new float[dimensions];
        roundedLoad = new float[dimensions];
        remainingCapacity = new float[dimensions];
        roundedRemainingCapacity = new float[dimensions];
        memcpy(load, _m.load, dimensions * sizeof(float));
        memcpy(roundedLoad, _m.roundedLoad, dimensions * sizeof(float));
        memcpy(remainingCapacity, _m.remainingCapacity, dimensions * sizeof(float));
        memcpy(roundedRemainingCapacity, _m.roundedRemainingCapacity, dimensions * sizeof(float));
        for(unsigned i=0; i<_m.assignedJobs->size(); i++) {
            assignedJobs->push_back(_m.assignedJobs->at(i));
        }
        for(unsigned i=0; i<_m.smallJobs->size(); i++) {
            smallJobs->push_back(_m.smallJobs->at(i));
        }
    }
    flaggedAsFull = _m.flaggedAsFull;
    bestJobIndex = _m.bestJobIndex;
    valueOfBestJob = _m.valueOfBestJob;
    return *this;
}

Machine::~Machine() {
    delete [] load;
    delete [] roundedLoad;
    delete [] remainingCapacity;
    delete [] roundedRemainingCapacity;
    delete assignedJobs;
    delete smallJobs;
}

void Machine::reset(float _cap) {
    capacity = _cap;
    for(int i=0; i<dimensions; i++) {
        load[i] = 0.0;
        roundedLoad[i] = 0.0;
        remainingCapacity[i] = _cap;
        roundedRemainingCapacity[i] = _cap;
    }
    if(assignedJobs) {
        assignedJobs->clear();
    }
    if(smallJobs) {
        smallJobs->clear();
    }
    flaggedAsFull = false;
    bestJobIndex = -1;
    valueOfBestJob = -1.;
}

bool Machine::assignJob(Job* _job) {
    float* tmp = _job->getKPIVec();
    if(doesVectorFit(tmp)) {
        assignedJobs->push_back(*_job);
        for(int i=0; i<dimensions; i++) {
            load[i] += tmp[i];
            remainingCapacity[i] -= tmp[i];
        }
        return true;
    } else {
        return false;
    }
}
 float* Machine::getJobKPI(int index){
     float* tmp = assignedJobs->at(index).getKPIVec();
     return tmp;
}
bool Machine::assignJobLTAA(Job* _job, bool small) {
    float* roundedVector = _job->getRoundedVector();
    // TODO: The if statement is probably not needed for large jobs, but maybe for small ones
    // TODO: Make sure that the right vectors are compared; (1, ..., 1) or (1+eps, ..., 1+eps)?
    if(small || doesRoundedVectorFit(roundedVector)) {
        assignedJobs->push_back(*_job);
        float* vector = _job->getKPIVec();
        for(int i=0; i<dimensions; i++) {
            load[i] += vector[i];
            roundedLoad[i] += roundedVector[i];
            remainingCapacity[i] -= vector[i];
            roundedRemainingCapacity[i] -= roundedVector[i];
        }
        if(small) {
            smallJobs->push_back(*_job);
        }
        return true;
    } else {
        return false;
    }
}

bool Machine::doesRoundedVectorFit(const float* roundedVector) {
    for(int i=0; i<dimensions; i++) {
        if(roundedLoad[i] + roundedVector[i] > capacity + 0.00001) { // + 0.00001 because of rounding effects 
            return false;
        }
    }
    return true;
};

bool Machine::doesVectorFit(const float* vector) {
    for(int i=0; i<dimensions; i++) {
        if(load[i] + vector[i] > capacity + 0.00001) { // + 0.00001 because of rounding effects 
            return false;
        }
    }
    return true;
};

// TODO: Only used in ChekuriVS; can we get rid of it?
void Machine::assignJobOverload(Job* _job) {
    float* tmp = _job->getKPIVec();
    for(unsigned i=0; i<dimensions; i++) {
        load[i] += tmp[i];
        //roundedLoad[i] += roundedVector[i];
        remainingCapacity[i] -= tmp[i];
        //roundedRemainingCapacity[i] -= roundedVector[i];
    }
    assignedJobs->push_back(*_job);
}


// TODO: Only used in LTAA; can we get rid of it?
void Machine::removeSmallJobs() {
    for(unsigned i=0; i<smallJobs->size(); i++) {
        this->popBackJobInit();
    }
    smallJobs->clear();
}
Job Machine::popJob(int index){
    Job job;
    if(!assignedJobs->empty()) {
        job = assignedJobs->at(index);
        float* vector = assignedJobs->at(index).getKPIVec();
        float* roundedVector = assignedJobs->at(index).getRoundedVector();
        for(unsigned i=0; i<dimensions; i++) {
            load[i] -= vector[i];
            roundedLoad[i] -= roundedVector[i];
            remainingCapacity[i] += vector[i];
            roundedRemainingCapacity[i] += roundedVector[i];
        }
        assignedJobs->erase(assignedJobs->begin()+index);
    }
    return job;
}
// TODO: Only used in LTAA; can we get rid of it?
Job Machine::popBackJobInit() {
    Job job;
    if(!assignedJobs->empty()) {
        job = assignedJobs->back();
        float* vector = assignedJobs->back().getKPIVec();
        float* roundedVector = assignedJobs->back().getRoundedVector();
        for(unsigned i=0; i<dimensions; i++) {
            load[i] -= vector[i];
            roundedLoad[i] -= roundedVector[i];
            remainingCapacity[i] += vector[i];
            roundedRemainingCapacity[i] += roundedVector[i];
        }
        assignedJobs->pop_back();
    }
    return job;
}

std::vector<Job>* Machine::getSmallJobs() {
    return smallJobs;
}

int Machine::getDimensions() {
    return dimensions;
}

float* Machine::getLoad() {
    return load;
}

float Machine::getSumOfLoad(){
    float sum;
    for(unsigned i=0; i<dimensions; i++) {
        sum += load[i];
    }
    return sum;
}

float* Machine::getRemainingCapacity() {
    return remainingCapacity;
}

// TODO: Be careful with this function because it does, for example, not adapt assignedJobs
void Machine::setVector(const std::vector<float> vec) {
    for(unsigned i=0; i<dimensions; i++) {
        load[i] = capacity - vec[i];
        remainingCapacity[i] = vec[i];
    }
}

void Machine::removeJob(int assignedID){
//    printf("size: %d id: %d\n", assignedJobs->size(), assignedID);
    float* tmp_jv = assignedJobs->at(assignedID).getKPIVec();
    for(int i=0; i<dimensions; i++) {
        load[i] -= tmp_jv[i];
        remainingCapacity[i] += tmp_jv[i];
    }
    assignedJobs->erase(assignedJobs->begin() + assignedID);
}

bool Machine::isFlaggedAsFull() {
    return flaggedAsFull;
}

void Machine::flagAsFull() {
    flaggedAsFull = true;
}

void Machine::setBestJobIndex(int index) {
    bestJobIndex = index;
}
int Machine::getSizeOfAssignedJobs() {
    return assignedJobs->size();
}
int Machine::getBestJobIndex() {
    return bestJobIndex;
}

void Machine::setValueOfBestJob(float value) {
    valueOfBestJob = value;
}

float Machine::getValueOfBestJob() {
    return valueOfBestJob;
}

void Machine::printLoad(bool printRemCap) {
    if(printRemCap) {
        printf("Machine id: %d; dimensions: %d; remaining capacity: [%5.4f", id, dimensions, remainingCapacity[0]);
    } else {
        printf("Machine id: %d; dimensions: %d; load: [%5.4f", id, dimensions, load[0]);
    }
    for(int i=1; i<dimensions; i++) {
        printf(", %5.4f", printRemCap ? remainingCapacity[i] : load[i]);
    }
    printf("]\n");
}

void Machine::printMachine() {
    printf("    Machine %d:\n", id);
    printf("      Load: %f", load[0]);
    for(int d=1; d<dimensions; d++) {
        printf(", %f", load[d]);
    }
    printf("      Remaining capacity: %f", remainingCapacity[0]);
    for(int d=1; d<dimensions; d++) {
        printf(", %f", remainingCapacity[d]);
    }
    printf("\n      All jobs:\n");
    for(int i=0; i<assignedJobs->size(); i++) {
        printf("        ");
        assignedJobs->at(i).print();
    }
    printf("      Small jobs:\n");
    for(int i=0; i<smallJobs->size(); i++) {
        printf("        ");
        smallJobs->at(i).print();
    }
}









MachineGenerator::MachineGenerator(int _dimensions) {
    dimensions = _dimensions; 
}

/**
 * This function creates a vector of machines. 
 * @param _numberOfMachines number of machines to be created.
 * @param _capacity capacity in each dimension, usually 1.
 * @return vector of machines with IDs 1, 2, ..., numberOfMachines.
 */
std::vector<Machine>* MachineGenerator::generateMachines(int _numberOfMachines, float _capacity) {
    std::vector<Machine>* mv = new std::vector<Machine>();
    for(int i=0; i<_numberOfMachines; i++) {
        mv->push_back(generateMachine(i+1, _capacity));
    }
    return mv;
}

/**
 * This function creates and returns a single machine with given ID and capacity.
 * @param _id ID of machine.
 * @param _capacity capacity for each dimension.
 * @return machine.
 */
Machine MachineGenerator::generateMachine(int _id, float _capacity) {
    return (Machine(_id, dimensions, _capacity));
}

