#ifndef SORTING_H_
#define SORTING_H_

#include "Job.h"
#include "Machine.h"
#include "MathFunctions.h"
#include <vector>
#include <algorithm>

class Sorting {
public:

    static const int NO_SORTING = 0;
    static const int LINF = 1;
    static const int SUM = 2;
    static const int L2 = 3;
    static const int LEXI = 4;
    static const int LEXI_REORDERED = 5;

    struct CompareFuncLinf {
        bool operator()(Job const &a, Job const &b) {
            return a.getMaxComponent() > b.getMaxComponent();
        }
    };

    struct CompareFuncSum {
        bool operator()(Job const &a, Job const &b) {
            return a.getSumOfComponents() > b.getSumOfComponents();
        }
    };
    
    struct CompareFuncLexi {
        bool operator()(Job const &a, Job const &b) {
            return a.compareLexicographicallyTo(b) >= 0;
        }
    };

    struct CompareFuncLexiReordered {
        int* order;
        bool operator()(Job const &a, Job const &b) {
            return a.compareLexicographicallyTo(b, order) >= 0;
        }
    };

    struct CompareFuncL2 {
        bool operator()(Job const &a, Job const &b) {
            return a.getSumOfSquares() > b.getSumOfSquares();
        }
    };
    
    static int* determineOrder(int* order, std::vector<Job>* jq) {
        // Sum up all job vectors
        int d = jq->at(0).getDimensions();
        float sum[d] = {0};
        for(unsigned i = 0; i < jq->size(); i++) {
            float* v = jq->at(i).getKPIVec();
            for(unsigned j = 0; j < d; j++) {
                sum[j] += v[j];
            }
        }
        // Save order (order[0] = index of largest sum, ..., index[d-1] = index of smallest sum)
        for(unsigned i = 0; i < d; i++) {
            int index = 0;
            for(unsigned j = 1; j < d; j++) {
                if(sum[j] > sum[j - 1]) {
                    index = j;
                }
            }
            order[i] = index;
            sum[index] = -1;
        }
        return order;
    }
    
    static void sort(std::vector<Job> *jq, int dimensions, int type) {
        switch(type) {
            case NO_SORTING:
                break;
            case LINF:
                std::sort(jq->begin(), jq->end(), CompareFuncLinf());
                break;
            case SUM:
                std::sort(jq->begin(), jq->end(), CompareFuncSum());
                break;
            case LEXI: 
                std::sort(jq->begin(), jq->end(), CompareFuncLexi());
                break;
            case LEXI_REORDERED: {
                int order[dimensions] = {0};
                determineOrder(order, jq);
                CompareFuncLexiReordered cf;
                cf.order = order;
                std::sort(jq->begin(), jq->end(), cf);
                break;
            }
            case L2: 
                std::sort(jq->begin(), jq->end(), CompareFuncL2());
            default:
                break;
        }
    }
};


#endif
