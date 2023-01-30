//Search keywords
// could related to statistics
#include "Job.h"
#include "Machine.h"
#include "MathFunctions.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
/**
 * Local Search Randomization.
 */

bool cmp(std::pair<int, float>& a,
         std::pair<int, float>& b)
{
    return a.second < b.second;
}
// Function to sort the map according
// to value in a (key-value) pairs
void sort(std::map<int, float>& M)
{
    std::vector<std::pair<int, float> > A;
    for (auto& it : M) {
        A.push_back(it);
    }
    sort(A.begin(), A.end(), cmp);
}
struct CompareFuncJQ{
    
    bool operator()(Job const &a, Job const &b) {
        return a.getSumOfComponents() > b.getSumOfComponents();
    }
};
int findBiggestJob(std::vector < Job >* _jobs, int type) {
    float max = 0.0;
    int index = -1;
    for(unsigned i = 0; i < _jobs->size(); i++){
        float sum = 0.0;
        if(type == 1)
            sum = _jobs->at(i).getSumOfComponents();
        else if(type == 2)
            sum = _jobs->at(i).getSumOfSquares();
        else if(type == 3)
            sum = _jobs->at(i).getMaxComponent();
        if(sum > max){
            max = sum;
            index = i;
        }
    }
    return index;
}
bool seachJID(std::vector < int >_mj, int _nj) {
    for(unsigned i = 0; i < _mj.size(); ++i) {
        if(_nj == _mj[i])
            return true;
    }
    return false;
}
extern "C"
int planJobs(int &counter, std::vector <Machine> **out_machines, std::vector <Job> *_jq, std::vector <Machine> * _mq, float _epsilon) {
   
    std::vector <Job>* aJQ = new std::vector <Job> ();
    float* tmp_jv;
    int dims = _jq->at(0).getDimensions();
    MachineGenerator mach_gen(dims);
    std::vector <Machine>* mq = mach_gen.generateMachines(1,1.0);
    float* total = new float[dims];
    for(unsigned j = 0; j < dims; j++)
    {
        total[j] = 0;
    }
    for(int i = 0; i < _jq->size(); i++) {
        aJQ->push_back((_jq->at(i)));
        tmp_jv = _jq->at(i).getKPIVec();
        for(unsigned j = 0; j < dims; j++)
        {
            total[j] += tmp_jv[j];
        }
    }
    for(unsigned i = 0; i < dims; i++){
        if(total[i] >= _mq->size()){
            counter++;
        *out_machines = mq;
        return 1;
        }
    }
    float maxOfTotal = 0;
    for(unsigned i = 0; i < dims; i++){
        if(total[i] > maxOfTotal)
            maxOfTotal = total[i];
    }
    float Avg = 0;
    for(unsigned j = 0; j < dims; j++)
    {
        Avg +=total[j];
    }
    Avg = Avg/(dims*(_jq->size()));
    float stdvation = 0;
    for(int i = 0; i < _jq->size(); i++) {
        tmp_jv = _jq->at(i).getKPIVec();
        float jAvg = 0;
        for(unsigned j = 0; j < dims; j++)
        {
            jAvg += tmp_jv[j];
        }
        jAvg = jAvg / dims;
        stdvation += pow((jAvg-Avg),2);
    }
    stdvation = std::sqrt(stdvation/(_jq->size()));
   int iniRound = int(((_mq->size()/Avg)*stdvation*10)*(1+Avg+stdvation));
    if(iniRound < 100)
        iniRound = 100;
    if(iniRound > 1000)
        iniRound = 1000;
    delete []total;
    int round ;
    float p=0.1;
    int goodd=0;
    int badd=0;
    float* tmp_mv;
    float* tmp_kpi;
    double dp;
    float BC;
    int BCID;
    int count = 0;
    float NC;
    float JC;
    int machineSize = 0;
    long jqSize = aJQ->size();
    bool quitFlag = false;
    int bin = 0;
    // start algorithm
    //terminate condition is reaching maximum machine size
    while(mq->size() <= _mq->size()) {
        bin++;
        float max = 0.0;
        float maxID = -1;
        //arrange random jobs into machine
        int randomSelect;
        int count = 0;
        mq->back().assignJob(&(aJQ->at(findBiggestJob(aJQ,3))));
        max = aJQ->at(findBiggestJob(aJQ,3)).getMaxComponent();
        aJQ->erase(aJQ->begin()+findBiggestJob(aJQ,3));
        float* initialResult;
        float* bestResult;
        float* finalResult;
        int timesOfT = 0;
        //randomly assign
        while((aJQ->size())!=0) {
            int job_index = -1;
            randomSelect = rand() % aJQ->size();
            tmp_jv = aJQ->at(randomSelect).getKPIVec();
            //error data
            for(unsigned k = 0; k < dims; k++){
                if(tmp_jv[k] > 1.0)
                    printf("ERROR: value to big!!! (%f)\n", tmp_jv[k]);
            }
            tmp_mv = mq->back().getRemainingCapacity();
            if(MathFunctions::isVectorDifferencePositive(tmp_mv, tmp_jv, dims)) {
                mq->back().assignJob(&(aJQ->at(randomSelect)));
                aJQ->erase(aJQ->begin()+randomSelect);
            }
            else{
                count++;
            }
            if(count > int(aJQ->size())) { // could related to statistics
                break;
                //              aJID.push_back(job_index);
            }
        }
        initialResult = mq->back().getRemainingCapacity();
        float initial = MathFunctions::DotProd(initialResult,initialResult,dims);
        float best = MathFunctions::DotProd(initialResult,initialResult,dims);
        if(aJQ->size()==0)
            break;
        int success = 10;
        int totalJob = aJQ->size() + mq->back().getSizeOfAssignedJobs();
        std::vector<int> bestAssignment;
        while(success>0)
        {
            success--;
            timesOfT++;
            round = iniRound;
            while(round>0){
                bool swap = false;
                int i = rand() % (aJQ->size());
                round--;
                tmp_mv = mq->back().getRemainingCapacity();
                tmp_jv = aJQ->at(i).getKPIVec();
                float min_good = 0;
                float start = 0;
                float min_bad = 100;
                int aJID_index_good=-1;
                int aJID_index_bad=-1;
                float costfunction;
                for(unsigned j = 0; j < dims; j++){
                    min_good += pow((tmp_mv[j]),2);
                }
                start = min_good;
                std::map<int,float>toBeSwap;
                if(mq->back().getSizeOfAssignedJobs()>1){
                    for (unsigned k = 1; k < mq->back().getSizeOfAssignedJobs();++k){
                        costfunction = 0;
                        tmp_kpi = mq->back().getJobKPI(k);
                        MathFunctions::VectorAddition(tmp_mv, tmp_kpi, dims);
                            // if the select item is able to be assigned into the bin
                        if(MathFunctions::isVectorDifferencePositive(tmp_mv, tmp_jv, dims)){
                            for(unsigned j = 0; j < dims; j++){
                                float difference =tmp_mv[j] - tmp_jv[j];
                                costfunction += pow(difference,2);
                            }
                            if(costfunction <= min_good){
                                if(costfunction != min_good){
                                    min_good = costfunction;
                                    aJID_index_good = k;
                                }else
                                {
                                    MathFunctions::VectorSubtraction(tmp_mv, tmp_kpi, dims);
                                    continue;
                                }
                            }
                            if(aJID_index_good == -1){
                                toBeSwap.insert(std::pair<int,float>(k,costfunction));
                            }
                        }
                        MathFunctions::VectorSubtraction(tmp_mv, tmp_kpi, dims);
                    }
                }
                if(aJID_index_good != -1){
                    Job removedJob;
                    removedJob = mq->back().popJob(aJID_index_good);
                    tmp_jv = removedJob.getKPIVec();
                    aJQ->push_back(removedJob);
                    tmp_mv = mq->back().getRemainingCapacity();
                    if(mq->back().assignJob(&(aJQ->at(i)))){
                        aJQ->erase(aJQ->begin() + i);
                        success = true;
                        swap = true;
                    }
                }
                    // replace BC by NC with p
                else if(toBeSwap.size()!=0){
                    sort(toBeSwap);
                    int candidateJob = 0;
                    std::map<int,float>::iterator it =toBeSwap.begin();
                    for(unsigned w = 0; w < candidateJob; w++)
                        it++;
                    //probability is the ratio of the original remainning capacity divided by
                    //the est to the swap
                    
			p = (min_good/it->second);
			p = 0.05;
                    float randomNumber =float((rand() % 100)) / 100.0;
                    if(randomNumber < p){
                        Job removedJob;
                        removedJob = mq->back().popJob(it->first);
                        aJQ->push_back(removedJob);
                        if(mq->back().assignJob(&(aJQ->at(i)))){
                            aJQ->erase(aJQ->begin() + i);
                            success = true;
                            swap = true;
                        }
                    }
                }
                count = 0;
                if(swap)
                {
                    //random assignment
                    while((aJQ->size())!=0) {
                        randomSelect = rand() % aJQ->size();
                        tmp_jv = aJQ->at(randomSelect).getKPIVec();
                            //error data
                        for(unsigned k = 0; k < dims; k++){
                            if(tmp_jv[k] > 1.0)
                                printf("ERROR: value to big!!! (%f)\n", tmp_jv[k]);
                        }
                        tmp_mv = mq->back().getRemainingCapacity();
                        if(MathFunctions::isVectorDifferencePositive(tmp_mv, tmp_jv, dims)) {
                            if(mq->back().assignJob(&(aJQ->at(randomSelect))))
                                aJQ->erase(aJQ->begin() + randomSelect);
                        }
                        else{
                            count++;
                        }
                        if(count > int(aJQ->size())) {
                            break;
                        }
                    }
                }
                if(aJQ->size() == 0 || aJQ->size() == 1)
                {
                    quitFlag = true;
                    break;
                }
            }
            if(quitFlag == true)
                break;
            tmp_mv = mq->back().getRemainingCapacity();
            float current = MathFunctions::DotProd(tmp_mv,tmp_mv,dims);
            if(current < best){
                bestAssignment.clear();
                for(unsigned a = 0; a < mq->back().getSizeOfAssignedJobs(); a++){
                    bestAssignment.push_back(mq->back().assignedJobs->at(a).getID());
                }
                best = current;
            }
        }
        totalJob =  mq->back().getSizeOfAssignedJobs();
        tmp_mv = mq->back().getRemainingCapacity();
        float current = MathFunctions::DotProd(tmp_mv,tmp_mv,dims);
        if(current == best){
            goodd++;
        }
        else{
            badd++;
        }
        // open new machine
        if(aJQ->size()>0 && (mq->size() < _mq->size()))
        {
            mq->push_back(mach_gen.generateMachine(mq->size() + 1, 1.0));
            machineSize++;
        }
        else
        {
            break;
        }
    }
    if(aJQ->size()!= 0)
    {
        std::sort(aJQ->begin(), aJQ->end(), CompareFuncJQ());
        for(unsigned i = 0; i < aJQ->size(); i++){
            for(unsigned j = 0; j < mq->size(); j++)
                if(mq->at(j).assignJob(&aJQ->at(i)))
                    break;
        }
    }
    if(aJQ->size()>0) {
        counter++;
    }
    delete aJQ;
    *out_machines = mq;
    return 1;
}





