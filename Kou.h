#ifndef KOU_H_
#define KOU_H_

#include "Job.h"
#include "Machine.h"
#include "MathFunctions.h"
#include "Sorting.h"
#include <vector>
#include <algorithm>

class Kou {
public:
   
    static bool runKouFFD(int &counter, std::vector<Machine> **out_machines, std::vector<Job> *jq, std::vector<Machine> *mq, float _eps, int type) {
        std::vector<Machine>* mq_clone = new std::vector<Machine>(*mq);
        int d = jq->at(0).getDimensions();
        Sorting::sort(jq, d, type);
        
        for(unsigned i = 0; i < jq->size(); i++) {
            for(unsigned j = 0; j < mq_clone->size(); j++) {
                if(mq_clone->at(j).assignJob(&(jq->at(i)))) {
                    break;
                } else {
                    if(j == mq_clone->size() - 1) {
                        *out_machines = mq_clone;
                        return false;
                    }
                }
            }
        }
        
        *out_machines = mq_clone;
        return true;
    }

    static bool runKouBFD(int &counter, std::vector<Machine> **out_machines, std::vector<Job> *jq, std::vector<Machine> *mq, float _eps, int type) {
        std::vector<Machine>* mq_clone = new std::vector<Machine>(*mq);
        int d = jq->at(0).getDimensions();
        Sorting::sort(jq, d, type);

        for(unsigned i = 0; i < jq->size(); i++) {
            // find the best fit for current job
            float min_sum = std::numeric_limits<float>::max();
            int index = -1;
            float* v = (jq->at(i)).getKPIVec();
            // Go through array of machines to find best fit
            for(unsigned j = 0; j < mq_clone->size(); j++) {
                float* remainder = (mq_clone->at(j)).getRemainingCapacity();
                float tmp_sum = 0.;
                for(unsigned k = 0; k < d; k++) {
                    float diff = remainder[k] - v[k];
                    if(diff < -0.000001) { // job does not fit
                        tmp_sum = std::numeric_limits<float>::max();
                        break;
                    } else { // component fits
                        tmp_sum += diff;
                    }
                }
                if(tmp_sum < min_sum) { // machine j is new best fit
                    index = j;
                    min_sum = tmp_sum;
                }
            }
            // If no machine found, return failure
            if(index == -1) { 
                *out_machines = mq_clone;
                return false;
            }
            // Assign job to best machine
            mq_clone->at(index).assignJob(&(jq->at(i)));
        }
        
        *out_machines = mq_clone;
        return true;
    }
};


#endif
