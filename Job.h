#ifndef JOB_H_
#define JOB_H_

#include <vector>
#include <cstring>
#include <cstdio>
#include <cfloat>
#include <cmath>

class JobMetaData {
public:
    float lowerBound;
    float upperBound;
    int lowerBurst;
    int upperBurst;
    int duration;
    int dimensions;
    int numberOfJobs;
    int arrivalSpan;

    JobMetaData() {
        lowerBound = 0;
        upperBound = 0;
        lowerBurst = 0;
        upperBurst = 0;
        duration = 0;
        dimensions = 0;
        numberOfJobs = 0;
        arrivalSpan = 0;
    }

    JobMetaData(int _dimensions, int _numberOfJobs = 1, float _lowerBound = 0, float _upperBound = 0, int _lowerBurst = 0, int _upperBurst = 0, int _duration = 0, int _arrivalSpan = 0) {
        lowerBound = _lowerBound;
        upperBound = _upperBound;
        lowerBurst = _lowerBurst;
        upperBurst = _upperBurst;
        duration = _duration;
        dimensions = _dimensions;
        numberOfJobs = _numberOfJobs;
        arrivalSpan = _arrivalSpan;
    }
};

/**
 * TODO: If guaranteed that job never changes, one could precompute sum of components / squares and (rounded) max component.
 */
class Job {
public:

    Job() {
        id = -1;
        kpiVec = 0;
        roundedVector = 0;
        bursttime = -1;
        arrivaltime = -1;
        dimensions = 0;
        large = false;
        zero = false;
        roundedVectorInitialized = false;
        weight = -1;
        counter = 0;
    };

    Job(const Job& _j) {
        id = _j.id;
        dimensions = _j.dimensions;
        kpiVec = 0;
        roundedVector = 0;
        roundedVectorInitialized = _j.roundedVectorInitialized;
        if(dimensions) {
            kpiVec = new float[dimensions];
            memcpy(kpiVec, _j.kpiVec, dimensions * sizeof (float));
            if(_j.roundedVectorInitialized) {
                roundedVector = new float[dimensions];
                memcpy(roundedVector, _j.roundedVector, dimensions * sizeof (float));
            }
        }
        bursttime = _j.bursttime;
        arrivaltime = _j.arrivaltime;
        large = false;
        zero = false;
        weight = _j.weight;
        counter = _j.counter;
    };

    Job(int _id, std::vector<float>* _vec, int _bursttime, int _arrivaltime = 0) {
        id = _id;
        dimensions = 0;
        if (_vec && _vec->size()) {
            dimensions = _vec->size();
            kpiVec = new float[_vec->size()];
            roundedVector = new float[_vec->size()];
            memcpy(kpiVec, &(_vec->at(0)), dimensions * sizeof (float));
            memcpy(roundedVector, &(_vec->at(0)), dimensions * sizeof (float));
        } else {
            kpiVec = 0;
            roundedVector = 0;
        }
        bursttime = _bursttime;
        arrivaltime = _arrivaltime;
        large = false;
        zero = false;
        roundedVectorInitialized = true;
        weight = -1;
        counter = 0;
    };

    void deleteData() {
        if(kpiVec) {
            delete[] kpiVec;
	    kpiVec = 0;
        }
        if(roundedVector) {
            delete[] roundedVector;
	    roundedVector = 0;
        }
    };

    virtual ~Job() {
        if(kpiVec) {
            delete[] kpiVec;
        }
        if(roundedVector) {
            delete[] roundedVector;
        }
    };

    Job& operator=(const Job& _j) {
        if(this == &_j)
            return *this;
        if(kpiVec)
        {
            delete[] kpiVec;
        }
        if(roundedVector) {
            delete[] roundedVector;
        }
        id = _j.id;
        dimensions = _j.dimensions;
        kpiVec = 0;
        roundedVector = 0;
        if (dimensions) {
            kpiVec = new float[dimensions];
            memcpy(kpiVec, _j.kpiVec, dimensions * sizeof (float));
            if(_j.roundedVectorInitialized) {
                roundedVector = new float[dimensions];
                memcpy(roundedVector, _j.roundedVector, dimensions * sizeof (float));
            }
        }
        bursttime = _j.bursttime;
        arrivaltime = _j.arrivaltime;
        large = _j.large;
        zero = _j.zero;
        roundedVectorInitialized = _j.roundedVectorInitialized;
        weight = _j.weight;
        counter = _j.counter;
        return *this;
    }

    float* getKPIVec() const {
        return kpiVec;
    }

    float* getRoundedVector() const {
        return roundedVector;
    }

    bool isLarge() const {
        return large;
    }

    bool isZero() const {
        return zero;
    }

    int getID() const {
        return id;
    }

    void print() {
        printf("Job id: %d; dims: %d; vector: [%f (%f)", id, dimensions, kpiVec[0], roundedVector[0]);
        for(int i=1; i<dimensions; i++) {
            printf(", %f (%f)", kpiVec[i], roundedVector[i]);
        }
        printf("]\n");
    }

    int getDimensions() const {
        return dimensions;
    }

    float getSumOfComponents() const {
        float res = 0.;
        for(unsigned i=0; i<dimensions; i++) {
            res += kpiVec[i];
        }
        return res;
    }
    
    float getSumOfSquares() const {
        float res = 0.;
        for(unsigned i=0; i<dimensions; i++) {
            res += kpiVec[i] * kpiVec[i];
        }
        return res;
    }

    float getMaxComponent() const {
        float max = -FLT_MAX;
        for(unsigned i=0; i<dimensions; i++) {
            if (max < kpiVec[i]) {
                max = kpiVec[i];
            }
        }
        return max;
    }

    float getMaxRoundedComponent() const {
        float max = -FLT_MAX;
        for(unsigned i=0; i<dimensions; i++) {
            if(max < roundedVector[i]) {
                max = roundedVector[i];
            }
        }
        return max;
    }
    
    void setWeight(float w) {
        weight = w;
    }

    float getWeight() const {
        return weight;
    }
    
    void setCounter(int c) {
        counter = c;
    }
    
    void incrementCounter() {
        counter++;
    }
    
    void decrementCounter() {
        counter--;
    }
    
    bool isCounterZero() const {
        return counter == 0;
    }

//    bool isCounterOne() const {
//        return counter == 1;
//    }

    int compareLexicographicallyTo(const Job& _j) const {
        float* v = _j.getKPIVec();
        for(unsigned i=0; i<dimensions; i++) {
            float d = kpiVec[i] - v[i];
            if(d > 0.000001) {
                return 1;
            }
            if(d < 0.000001) {
                return -1;
            }
        }
        return 0;
    }

    int compareLexicographicallyTo(const Job& _j, int* order) const {
        float* v = _j.getKPIVec();
        for(unsigned i=0; i<dimensions; i++) {
            int ii = order[i];
            float d = kpiVec[ii] - v[ii];
            if(d > 0.000001) {
                return 1;
            }
            if(d < 0.000001) {
                return -1;
            }
        }
        return 0;
    }

    void roundChekuri(float _delta, float _epsilon) {
        // get maximum vector component
        float v_max = getMaxComponent();
        // determine if job is "large"
        large = v_max > _delta;
        // make roundedVector a copy of kpiVec, but round every small vector component (< delta * v_max) down to 0
        for(unsigned i=0; i<dimensions; i++) {
            if(kpiVec[i] <= _delta * v_max) {
                roundedVector[i] = 0;
            } else {
                roundedVector[i] = kpiVec[i];
            }
        }
        // if job is small, we are done
        if(!large) {
            return;
        }
        // round each coordinate of large vectors down to left end point of interval they are in
        for(unsigned i=0; i<dimensions; i++) {
            if(kpiVec[i] != 0) {
                roundedVector[i] = intervalRounding((_delta * _delta), kpiVec[i], _epsilon);
            }
        }
        roundedVectorInitialized = true;
    }
    
    void roundLTAA(float _epsilon) {
        double eps_divided_by_d = _epsilon / static_cast<float>(dimensions);
        double eps3_divided_by_d = _epsilon * _epsilon * eps_divided_by_d;
        // delta = epsilon^4 / d^2   (not called delta in paper)
        double delta = eps3_divided_by_d * eps_divided_by_d;
        // Round every vector coordinate down to nearest value in {0, delta, delta*(1+eps), delta*(1+eps)^2, ..., 1}
        for(unsigned i=0; i<dimensions; i++) {
            // TODO: MathFunctions::IntervalRounding(..., false, ...); was called before; but using true should not make a difference
            roundedVector[i] = kpiVec[i] < delta ? 0 : intervalRounding(delta, kpiVec[i], _epsilon);
        }
        // Set coordinates to 0 that are small compared to largest coordinate in the vector
        float max = getMaxRoundedComponent();
        float limit = max * eps_divided_by_d;
        for(unsigned i=0; i<dimensions; i++) {
            if(roundedVector[i] < limit) {
                roundedVector[i] = 0;
            }
        }
        // Compute if vector is large or small
        large = max >= eps3_divided_by_d;
        zero = max == 0;
        roundedVectorInitialized = true;
    }

    /**
     * This function returns the left end point of the interval containing num.
     * That is, it returns left if left < num <= (1+_eps)*num, otherwise it 
     * calls itself with (1+eps)*left and returns the result of this call.
     * TODO: Interval rounding could be implemented more efficiently.
     * @param left left end point of interval.
     * @param num number for which the fitting interval is determined.
     * @param _eps used to determine right end point of interval which is (1+_eps)*left.
     * @return left end point of interval that contains num.
     */
    double intervalRounding(float left, float num, float epsilon) {
        float right = (1.0 + epsilon) * left;
        if(left < num && num <= right) {
            return left;
        } else {
            return intervalRounding(right, num, epsilon);
        }
    }


private:
    int id;
    float* kpiVec;
    float* roundedVector;
    int bursttime;
    int arrivaltime;
    int dimensions;
    float weight;
    int counter;
    bool large;
    bool zero;
    bool roundedVectorInitialized;
};


#endif

