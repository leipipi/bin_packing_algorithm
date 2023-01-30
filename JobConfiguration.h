#ifndef JOBCONFIGURATION_H_
#define JOBCONFIGURATION_H_

#include <vector>
#include <cstring>
#include <cstdio>
#include <cfloat>
#include <cmath>

/**
 * This class is only used by the LTAA algorithms. A job configuration is a multiset B of
 * rounded large jobs (called types). The respective jobs added up must not exceed (1, 1, ..., 1).
 * The vector types represent this multiset. types_max is a vector indicating the maximum number
 * of jobs for each type in the problem instance; types_max is therefore the same in each 
 * JobConfiguration. The profile f represents the remaining capacity which can be filled with small jobs.
 * The load is the sum of the jobs of the configuration.
 */
class JobConfiguration {
private:
    // Dimensions of job vectors
    int dimensions;
    // Epsilon of approximation scheme; bins can be overloaded up to (1+eps, ..., 1+eps)
    float epsilon;
    // Vector indicating the number of jobs of each type (multiset B in the paper)
    std::vector<int> types;
    // Vector indicating the maximum number of jobs of each type in the problem instance
    std::vector<int> types_max;
    // Load of all the rounded jobs
    float* load;
    // Profile f, i.e. the remaining capacity which can be used for small jobs
    float* f; 

public:

    /**
     * Constructor that creates configuration with an empty multiset B of large jobs.
     * @param _epsilon epsilon.
     * @param _dimensions dimensions of job vectors.
     * @param _types_max vector storing for each type t the maximum number of jobs of type t.
     */
    JobConfiguration(float _epsilon, int _dimensions, const std::vector<int> &_types_max) {
        epsilon = _epsilon;
        dimensions = _dimensions;
        // For each type, set the maximum number of jobs of that type in the instance
        types_max = std::vector<int>(_types_max);
        // Create load vector
        load = new float[dimensions];
        std::fill(load, load + dimensions, 0.0);
        // Create profile f
        f = new float[dimensions];
        setF();
        // For each type, set the number of jobs of that type in this configuration to 0
        types = std::vector<int>(types_max.size(), 0);
    }

    JobConfiguration(float _epsilon, int _dimensions, const std::vector<int> &_types, const std::vector<int> &_types_max, float* _load) {
        epsilon = _epsilon;
        dimensions = _dimensions;
        types = std::vector<int>(_types);
        types_max = std::vector<int>(_types_max);
        f = new float[dimensions];
        load = new float[dimensions];
        memcpy(load, _load, dimensions * sizeof(float));
    }

    JobConfiguration(const JobConfiguration& _jc) {
        epsilon = _jc.epsilon;
        dimensions = _jc.dimensions;
        types = std::vector<int>(_jc.types);
        types_max = std::vector<int>(_jc.types_max);
        f = new float[dimensions];
        memcpy(f, _jc.f, dimensions * sizeof(float));
        load = new float[dimensions];
        memcpy(load, _jc.load, dimensions * sizeof(float));
    }

    virtual ~JobConfiguration() {
        delete [] f;
    }

    JobConfiguration& operator=(const JobConfiguration& _jc) {
        epsilon = _jc.epsilon;
        dimensions = _jc.dimensions;
        types = std::vector<int>(_jc.types);
        types_max = std::vector<int>(_jc.types_max);
        f = new float[dimensions];
        memcpy(f, _jc.f, dimensions * sizeof (float));
        load = new float[dimensions];
        memcpy(load, _jc.load, dimensions * sizeof(float));
        return *this;
    }

    void print() {
        printf("Configuration: d = %d, eps = %f\n", dimensions, epsilon);
        for(int i=0; i<types.size(); i++) {
            printf("  Type: %d, types_max: %d\n", types[i], types_max[i]);
        }
        printf("  Profile: ");
        for(int i=0; i<dimensions; i++) {
            printf(" %f,", f[i]);
        }
        printf("\n");
    }

    void setF() {
        std::fill(f, f + dimensions, 0.0);

        for(unsigned i=0; i<dimensions; i++) {
            unsigned int exp = 0;
            float last = 0.0;
            while(load[i] + f[i] <= 1.0 + epsilon) {
                last = f[i];
                if(f[i] < 1.0 && std::pow(1 + epsilon, exp) * epsilon > 1.0) {
                    f[i] = 1.0;
                } else {
                    f[i] = std::pow(1.0 + epsilon, exp) * epsilon;
                    exp++;
                }
            }
            f[i] = last;
        }
    }

    float* getF() {
        return f;
    }

    float* getLoad() {
        return load;
    }

    std::vector<int> getTypes() {
        return types;
    }

    std::vector<int> getTypesMax() {
        return types_max;
    }

    int getDimensions() {
        return dimensions;
    }
};


#endif

