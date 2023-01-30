#ifndef LTAA_H_
#define LTAA_H_

#include "Job.h"
#include "JobConfiguration.h"
#include <string>
#include <algorithm>
#include <random>

#include <ilcplex/ilocplex.h>
#include <vector>

class LTAA {

private:
    
    static void generateConfigurationsRecursively(std::vector<std::vector<Job>*>* jobs, std::vector<JobConfiguration>* configs, JobConfiguration* config, int typeIndex, int jobIndex, float epsilon, int dimensions) {
        // Load of new configuration
        float* newLoad = new float[dimensions];
        // Load of old configuration that is to be extended
        float* load = config->getLoad();
        // Vector storing maximum number of jobs for each type; just passed on
        std::vector<int> types_max = config->getTypesMax();
        // Vector storing number of jobs for each type in old configuration
        std::vector<int> types = config->getTypes();
        // Loop over types and try to extend old configuration by one job of type t
        for(int t=typeIndex; t<jobs->size(); t++) {
            // Load of job of type t (we can always pick the first job because all jobs of the same type have the same load)
            float* jobLoad = jobs->at(t)->at(0).getRoundedVector();
            // Check if adding the new job exceeds the limits of a machine of size (1, 1, ..., 1)
            bool newJobFits = true;
            for(int d=0; d<dimensions; d++) {
                newLoad[d] = load[d] + jobLoad[d];
                if(newLoad[d] > 1.00000001) {
                    newJobFits = false;
                    break;
                }
            }
            // If the job fits, create a new configuration by extending the old configuration with this job
            if(newJobFits) {
                // Add the type of the new job by incrementing the type counter in the types array
                std::vector<int> newTypes(types);
                newTypes[t]++;
                // Create and add configuration
                JobConfiguration newConfig(epsilon, dimensions, newTypes, types_max, newLoad);
                configs->push_back(newConfig);
                // Continue with next job of the same type or, if there is none left, with the first job of the next type
                if(jobIndex + 1 < jobs->at(t)->size()) {
                    generateConfigurationsRecursively(jobs, configs, &newConfig, t, jobIndex+1, epsilon, dimensions);
                } else {
                    generateConfigurationsRecursively(jobs, configs, &newConfig, t+1, 0, epsilon, dimensions);
                }
            }
        }
        delete [] newLoad;
    }

public:

    /**
     * This function gets the set of large jobs and creates all job configurations. A job
     * configuration is a subset of (rounded) large jobs that would into a machine of size (1, ..., 1).
     * This function creates the first configuration containing no jobs and then calls 
     * generateConfigurationsRecursively(...) to create all other job configurations recursively.
     * @param _jobs vector of vector of large jobs; each subvector contains jobs of the same type.
     * @param _configs set of configurations; when new configurations are created, they are added to it.
     * @param _epsilon epsilon of the approximation algorithm.
     * @param dimensions dimensions of the job vectors.
     * @return false if there are no large jobs, otherwise true.
     */
    static bool generateConfigurations(std::vector<std::vector<Job>*>* _jobs, std::vector<JobConfiguration>* _configs, float _epsilon, int dimensions) {
        // Check if there are any large jobs
        if(_jobs->empty()) { 
            // Create configuration with no jobs and return
            std::vector<int> tmp_types(1, 0);
            JobConfiguration config(_epsilon, dimensions, tmp_types);
            _configs->push_back(config);
            return false;
        }
        // For each type i of job vector store the number of occurrences in types_max[i]
        std::vector<int> types_max(_jobs->size(), 0);
        for(unsigned i=0; i<_jobs->size(); i++) {
            types_max[i] = _jobs->at(i)->size();
        }
        // Create first configuration with no jobs
        JobConfiguration config(_epsilon, dimensions, types_max);
        // Recursively create all valid configurations and add them to configs
        generateConfigurationsRecursively(_jobs, _configs, &config, 0, 0, _epsilon, dimensions);
        // Add first configuration (with no jobs) at last position 
        _configs->push_back(config);
        // Compute the profile for every configuration
        for(int i=0; i<_configs->size(); i++) {
            _configs->at(i).setF();
        }
        return true;
    }

    /**
     * This function checks if two small job vectors are of the same type. The job vectors a and b are
     * of the same type if a / ||a||_max = b / ||b||_max.
     * @param _a first small job.
     * @param _b second small job.
     * @param _dimensions dimensions of vectors.
     * @return true if and only if _a and _b are of the same small type.
     */
    static bool compareSmallTypes(float* _a, float* _b, int _dimensions) {
        float aLinf = MathFunctions::Linf(_a, _dimensions);
        float bLinf = MathFunctions::Linf(_b, _dimensions);
        for(unsigned i=0; i<_dimensions; i++) {
            if(!MathFunctions::areFloatsEqual(_a[i] / aLinf, _b[i] / bLinf)) {
                return false;
            }
        }
        return true;
    }

    static void printVectors(std::vector<std::vector<Job>*>* S, std::vector<std::vector<Job>*>* L, std::vector<Job>* null_vecs) {
        printVectors(S, L);
        int dimensions = null_vecs->empty() ? 0 : null_vecs->at(0).getDimensions();
        for(int i=0; i<null_vecs->size(); i++) {
            printf("N job %.3d     vec: ( ", i);
            for(int k=0; k<dimensions; k++) {
                printf("%f ", null_vecs->at(i).getKPIVec()[k]);
                printf("[%f], ", null_vecs->at(i).getRoundedVector()[k]);
            }
            printf(")\n");
        }
    }

    static void printVectors(std::vector<std::vector<Job>*>* S, std::vector<std::vector<Job>*>* L) {
        int dimensions = S->empty() ? L->empty() ? 0 : L->at(0)->at(0).getDimensions() : S->at(0)->at(0).getDimensions();
        for(int i=0; i<S->size(); i++) {
            for(int j=0; j<S->at(i)->size(); j++) {
                printf("S job %.3d %.3d id: %d vec: (", i, j, S->at(i)->at(j).getID());
                for(int k=0; k<dimensions; k++) {
                    printf("%f ", S->at(i)->at(j).getKPIVec()[k]);
                    printf("[%f], ", S->at(i)->at(j).getRoundedVector()[k]);
                }
                printf(")\n");
            }
        }
        for(int i=0; i<L->size(); i++) {
            for(int j=0; j<L->at(i)->size(); j++) {
                printf("L job %.3d %.3d id: %d vec: (", i, j, L->at(i)->at(j).getID());
                for(int k=0; k<dimensions; k++) {
                    printf("%f ", L->at(i)->at(j).getKPIVec()[k]);
                    printf("[%f], ", L->at(i)->at(j).getRoundedVector()[k]);
                }
                printf(")\n");
            }
        }
    }
    
    
    /**
     * This function goes through the job queue jq and sorts its jobs into three categories. Jobs (0, ..., 0) are put in null_vecs.
     * Small (large) vectors are placed in S (L) which is a vector of vectors. Each ob the subvectors stores jobs of the same type. 
     * @param jq all jobs.
     * @param null_vecs all jobs with only zero coordinates.
     * @param S all small jobs (except for zero jobs).
     * @param L all large jobs.
     */
    static void categorizeJobs(std::vector<Job> *jq, std::vector<Job>* null_vecs, std::vector<std::vector<Job>*>* S, std::vector<std::vector<Job>*>* L) {
        int dimensions = jq->at(0).getDimensions();
        for(unsigned i=0; i<jq->size(); i++) {
            if(jq->at(i).isZero()) {
                null_vecs->push_back(jq->at(i));
                continue;
            }
            float* roundedVec = jq->at(i).getRoundedVector();
            bool found = false;
            if(jq->at(i).isLarge()) {
                for(unsigned j=0; j<L->size(); j++) {
                    if(MathFunctions::areVectorsEqual(L->at(j)->at(0).getRoundedVector(), roundedVec, dimensions)) {
                        L->at(j)->push_back(jq->at(i));
                        found = true;
                        break;
                    }
                }
                if(!found) {
                    std::vector<Job>* temp = new std::vector <Job>();
                    temp->push_back(jq->at(i));
                    L->push_back(temp);
                }
            } else {
                for(unsigned j=0; j<S->size(); j++) {
                    if(compareSmallTypes(S->at(j)->at(0).getRoundedVector(), roundedVec, dimensions)) {
                        S->at(j)->push_back(jq->at(i)); 
                        found = true;
                        break;
                    }
                }
                if(!found) {
                    std::vector<Job>* temp = new std::vector <Job>();
                    temp->push_back(jq->at(i));
                    S->push_back(temp);
                }
            }
        }
    }

    
    static void prepareProfilesForMILP(std::vector<JobConfiguration>* _configs, std::vector<float*>* f_vec, std::vector<int>* f_conf_indx) {
        f_vec->push_back(_configs->at(0).getF());
        f_conf_indx->push_back(0);
        // Fill f_vec and f_conf_indx
        for(unsigned i=1; i<_configs->size(); i++) {
            bool found = false;
            // Check if profile has been seen before
            for(unsigned j=0; j<f_vec->size(); j++) {
                if(MathFunctions::areVectorsEqual(f_vec->at(j), _configs->at(i).getF(), _configs->at(i).getDimensions())) {
                    // Store index of previous occurrence
                    f_conf_indx->push_back(j);
                    found = true;
                    break;
                }
            }
            if(!found) {
                // Add index of new profile in f_vec to vector of indices
                f_conf_indx->push_back(f_vec->size());
                // Store new profile 
                f_vec->push_back(_configs->at(i).getF());
            }
        }
    }

    
    static void assignLargeJobsDeterministic(std::vector<JobConfiguration>* _configs, std::vector<std::vector<Job>*>* _L, std::vector<Machine>* _machines, std::vector<int>* f_conf_indx, std::vector<int>* chosen_configs, std::vector<int>* mach_f) {
        if(_L->empty()) {
            return;
        }
        // Assign chosen configurations to machines
        // Note: the number of configurations can be less than the number of machines
        for(unsigned i=0; i<chosen_configs->size(); i++) {
            std::vector<int> tmp_types = _configs->at(chosen_configs->at(i)).getTypes();
            for(unsigned j=0; j<tmp_types.size(); j++) {
                int job_count = tmp_types[j];
                while(job_count > 0) {
                    if(!_L->at(j)->empty()) {
                        _machines->at(i).assignJobLTAA(&(_L->at(j)->back()), false);
                        _L->at(j)->pop_back();
                        job_count--;
                    } else {
                        break;
                    }
                }
            }
            // Save identifiers of chosen configurations to retrieve profiles when assigning small jobs
            mach_f->at(i) = f_conf_indx->at(chosen_configs->at(i));
        }
    }

    
    static bool assignSmallJobsDeterministic(IloCplex* cplex, IloArray<IloNumVarArray>* y, std::vector<std::vector<Job>*>* _S, std::vector<Machine>* _machines, std::vector<int>* mach_f, float _epsilon, int dimensions) {
        // If there are no small jobs, return
        if(_S->empty()) {
            return true;
        }
        std::vector<std::vector<int>> types_mach;
        int mach_count = 0;
        bool check;
        
        for(unsigned i=0; i<_S->size(); i++) {
            check = false;
            for(unsigned j=0; j<mach_f->size(); j++) {
                if(cplex->getValue((*y)[i][mach_f->at(j)]) > 0) { 
                    if(!check) {
                        types_mach.push_back(std::vector<int>(1, j));
                        check = true;
                        continue;
                    }
                    types_mach[i].push_back(j);
                }
            }
            if(!check) {
                types_mach.push_back(std::vector<int>(1, mach_count % static_cast<int>(_machines->size())));
                mach_count++;
            }
        }
        for(unsigned i=0; i<_S->size(); i++) {
            std::vector<Job>* tmp_jobs = _S->at(i);
            for(unsigned j=0; j<tmp_jobs->size(); j++) {
                double min_phi = DBL_MAX;
                int mach_index = 0;
                std::vector<int> type_mach = types_mach[i];
                for(unsigned k=0; k<type_mach.size(); k++) {
                    double tmp_phi = phiBansal(&(_machines->at(type_mach[k])), &(tmp_jobs->at(j)), dimensions, _epsilon);
                    if(tmp_phi < min_phi) {
                        min_phi = tmp_phi;
                        mach_index = type_mach[k];
                    }
                }
                _machines->at(mach_index).assignJobLTAA(&(tmp_jobs->at(j)), true);
            }
        }
        // distribute jobs of overloaded machines
        for(unsigned i=0; i<_machines->size(); i++) {
            float small_job_sum[dimensions];
            std::fill(small_job_sum, small_job_sum + dimensions, 0.0);
            std::vector<Job>* small_jobs = _machines->at(i).getSmallJobs();

            for(unsigned j=0; j<small_jobs->size(); j++) {
                MathFunctions::VectorAddition(small_job_sum, small_jobs->at(j).getRoundedVector(), dimensions);
            }
            if(MathFunctions::isVectorNegative(_machines->at(i).getRemainingCapacity(), dimensions)) {
                _machines->at(i).removeSmallJobs();
                bool mach_found;
                while(!small_jobs->empty()) {
                    mach_found = false;
                    for(unsigned j=0; j<_machines->size(); j++) {
                        if(_machines->at(j).assignJobLTAA(&(small_jobs->back()), true)) {
                            mach_found = true;
                            small_jobs->pop_back();
                            if(small_jobs->empty()) {
                                break;
                            }
                        }
                    }
                    // could not pack overloaded jobs
                    if(!mach_found) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    static double phiBansal(Machine* _machine, Job* _job, int _dimensions, float _epsilon) {
        // calculate phi as described in the paper
        // TODO: I could not find it in the paper.
        double phi = 0.0;
        double lambda = (1.0 / 2) * std::log2(std::pow(_dimensions, 3) / std::pow(2, 3));
        std::vector<Job>* smallJobs = _machine->getSmallJobs();
        float* jobVector = _job->getRoundedVector();
        double L;
        double f;

        for(int i=0; i<_dimensions; i++) {
            L = 0.0;
            f = 0.0;
            L += jobVector[i];

            if(!smallJobs->empty()) {
                for(unsigned j=0; j<smallJobs->size(); j++) {
                    L += smallJobs->at(j).getRoundedVector()[i];
                }
            }

            unsigned int ex = 0;
            double last = 0.0;
            double dif = _machine->getLoad()[i] - L;

            while(dif + f <= 1.0 + _epsilon) {
                last = f;
                if(f < 1 && std::pow(1.0 + _epsilon, ex) * _epsilon > 1) {
                    f = 1;
                } else {
                    f = std::pow(1.0 + _epsilon, ex) * _epsilon;
                    ex++;
                }
            }
            f = last;

            phi += std::pow(2, lambda * (L - f));
        }
        return phi;
    }

    static void assignLargeJobsRandomized(std::vector<JobConfiguration>* _configs, std::vector<std::vector<Job>*>* _L, std::vector<Machine>* _machines, std::vector<int>* chosen_configs, std::vector<float*>* mach_f, float _epsilon, int dimensions) {
        // TODO: What is done here?
        // mach_f includes the profiles of all machines given by the assigned config
        // first fill up mach_f with the empty config profiles
        float max_f[dimensions];
        std::fill(max_f, max_f + dimensions, 0.0);
        float null_vec[dimensions];
        std::fill(null_vec, null_vec + dimensions, 0.0);
        for(unsigned i=0; i<dimensions; i++) {
            unsigned int exp = 0;
            float last = 0.0;
            while(null_vec[i] + max_f[i] <= 1.0 + _epsilon) {
                last = max_f[i];
                if(max_f[i] < 1 && std::pow(1.0 + _epsilon, exp) * _epsilon > 1) {
                    max_f[i] = 1;
                } else {
                    max_f[i] = std::pow(1.0 + _epsilon, exp) * _epsilon;
                    exp++;
                }
            }
            max_f[i] = last;
        }
        // TODO: What is done here?
        for(unsigned i=0; i<_machines->size(); i++) {
            mach_f->push_back(new float[dimensions]);
            for(int j=0; j<dimensions; j++) {
                (*mach_f)[i][j] = max_f[j];
            }
        }
        // If there are no large jobs, return
        if(_L->empty()) {
            return;
        }
        // Assign large jobs
        for(unsigned i=0; i<chosen_configs->size(); i++) {
            std::vector<int> tmp_types = _configs->at(chosen_configs->at(i)).getTypes();
            for(unsigned j=0; j<tmp_types.size(); j++) {
                int job_count = tmp_types[j];
                while(job_count > 0) {
                    if(!_L->at(j)->empty()) {
                        _machines->at(i).assignJobLTAA(&(_L->at(j)->back()), false);
                        _L->at(j)->pop_back();
                        job_count--;
                    } else {
                        break;
                    }
                }
            }
            // note the profile of each machine
            float* t = _configs->at(chosen_configs->at(i)).getF();
            for(int j=0; j<dimensions; j++) {
                (*mach_f)[i][j] = t[j];
            }
        }
    }

    static void printArray(float* arr) {
        for(int i=0; i<3; i++) {
            printf(" .. %f ", arr[i]);
        }
        printf("\n ");
    }
    
    
    static bool assignSmallJobsRandomized(IloCplex* cplex, IloArray<IloNumVarArray>* y, std::vector<std::vector<Job>*>* _S, std::vector<Machine>* _machines, std::vector<float*>* f_vec, std::vector<float*>* mach_f, int dimensions) {
        // If no small jobs, return
        if(_S->empty()) {
            return true;
        }
        // Setup random number generator
#ifdef RANDOM
        std::random_device device;
        std::mt19937 generator(device());
#else
        char* seed = getenv("RAND_SEED");
        std::mt19937 generator;
        if(!seed) {
            generator = std::mt19937(1234);
        } else {
            generator = std::mt19937(atoi(seed));
        }
#endif
        // Randomized algorithm for assigning the small jobs
        for(unsigned t=0; t<_S->size(); t++) {
            for(unsigned j=0; j<_S->at(t)->size(); j++) {
                // beta(f, t) = y_{f, t} / sum_{g in F} y_{g, t}
                // where t is a small type and F the set of profiles
                // y(f, t) sums up the largest coordinates of the jobs from J(f, t)
                // where J(f, t) is the set of small jobs assigned to configurations having profile f
                float betaDenominator = 0.0;
                // Keep track of profiles g with y(f, t) > 0
                std::vector<int> profile_indices;
                // Loop over all profile vectors
                for(unsigned g=0; g<f_vec->size(); g++) {
                    // Add y_{g, t} to sum
                    float tmp_y = cplex->getValue((*y)[t][g]);
                    if(tmp_y > 0.0) {
                        betaDenominator += tmp_y;
                        profile_indices.push_back(g);
                    }
                }
                // Define uniform distribution over interval [0, sum_{g in F} y_{g, t})
                std::uniform_real_distribution<float> beta(0.0, betaDenominator);
                // Compute random value in interval 
                float rand_beta = beta(generator);
                // Index of randomly chosen profile
                int index_f = 0;
                // Choose profile at random depending on the weights y_{f, t}
                for(unsigned k=0; k<profile_indices.size(); k++) {
                    // Get first positive y_{f, t}
                    float tmp_y = cplex->getValue((*y)[t][profile_indices[k]]);
                    // If random number < y_{f, t}, choose profile
                    if(rand_beta < tmp_y) {
                        index_f = profile_indices[k];
                        break;
                    }
                    // Reduce random number by y_{f, t}
                    rand_beta -= tmp_y;
                }

                // Collect the indices of the machines with the same profile
                std::vector<int> current_f_indices;
                for(unsigned l=0; l<mach_f->size(); l++) {
                    if(MathFunctions::areVectorsEqual(mach_f->at(l), (*f_vec)[index_f], dimensions)) {
                        current_f_indices.push_back(l);
                    }
                }
                // If several machines have profile, choose machine uniformly at random and assign job
                if(current_f_indices.size() > 1) {
                    std::uniform_int_distribution<int> uni(0, current_f_indices.size() - 1);
                    _machines->at(current_f_indices[uni(generator)]).assignJobLTAA(&(_S->at(t)->at(j)), true);
                } else {
                    // If only one machine has profile, assign job to that machine; otherwise return false
                    if(!current_f_indices.empty()) {
                        _machines->at(current_f_indices[0]).assignJobLTAA(&(_S->at(t)->at(j)), true);
                    } else {
                        return false;
                    }
                }
            }
        }

        // Redistribute jobs of every overloaded machine
        for(unsigned i=0; i<_machines->size(); i++) {
            // Check if machine i is overloaded
            if(MathFunctions::isVectorNegative(_machines->at(i).getRemainingCapacity(), dimensions)) {
                // Take all small jobs from machine i
                std::vector<Job>* smallJobs = _machines->at(i).getSmallJobs();
                _machines->at(i).removeSmallJobs();

                bool machineFound;
                while(!smallJobs->empty()) {
                    machineFound = false;
                    for(unsigned j=0; j<_machines->size(); j++) {
                        if(_machines->at(j).assignJobLTAA(&smallJobs->back(), true)) {
                            machineFound = true;
                            smallJobs->pop_back();
                            if(smallJobs->empty()) {
                                break;
                            }
                        }
                    }
                    // Could not pack overloaded jobs
                    if(!machineFound) {
                        return false;
                    }
                }
            }
        }
        return true;
    }
    
    
    static bool solveMILP(std::vector<JobConfiguration>* _configs, std::vector<std::vector<Job>*>* _S, std::vector<std::vector<Job>*>* _L, std::vector<Machine>* _machines, float _epsilon, bool randomized) {
        // Dimension of job vectors
        int dimensions = _configs->at(0).getDimensions();
        // Create CPLEX environment
        IloEnv env;
        // Create IloModel which is a set of constraints / objectives
        IloModel model(env);
        // Define t_big[i] as the total number of large jobs of type i
        std::vector<int> t_big;
        if(!_configs->empty()) {
            t_big = _configs->at(0).getTypesMax();
        }
        // Vector of (different) profiles 
        std::vector<float*> f_vec;
        // Vector that stores the index of profile in f_vec for each configuration
        // That is, f_vec.at(f_conf_indx.at(i)) would be the profile of configuration i
        std::vector<int> f_conf_indx;
        // Fill the vectors f_vec and f_conf_indx
        prepareProfilesForMILP(_configs, &f_vec, &f_conf_indx);
        // Define variables x_C (number of machines with jobs according to configuration C)
        IloNumVarArray x(env, _configs->size(), 0, IloInfinity, ILOINT);
        // Define variables y_{f, t} (largest coordinate of type t * size of J(f, t) for each small job type and profile)
        IloArray<IloNumVarArray> y(env, _S->size());
        for(unsigned i=0; i<_S->size(); i++) {
            y[i] = IloNumVarArray(env, f_vec.size(), 0, IloInfinity, ILOFLOAT);
        }
        // C1: sum_{c \in C} x_c * n(c, t) >= n(t)  for all t \in T_big
        // x_c     -- #machines that are assigned jobs according to c
        // n(c, t) -- number of big jobs of type t in configuration c
        // n(t)    -- total number of big jobs of type t
        for(unsigned t=0; t<t_big.size(); t++) {
            IloExpr v(env);
            for(unsigned c=0; c<_configs->size(); c++) {
                // Get vector of types in configuration
                std::vector<int> tmp_types = _configs->at(c).getTypes();
                // Add x_C * n(C, t) to sum
                v += x[c] * tmp_types[t];
            }
            // Add C1 to model 
            model.add(v >= t_big[t]);
            v.end();
        }
        // C2: sum_{f \in F} y_{f, t} >= a(t)  for all t \in T_small
        // f        -- profile
        // J(f, t)  -- Set of small jobs of type t assigned to configurations with profile f
        // y_{f, t} -- largest coordinate of type t * size of J(f, t)
        // a(t)     -- largest coordinate of type t * number of jobs of type t 
        for(unsigned i=0; i<_S->size(); i++) {
            IloExpr v(env);
            // TODO: Is this simply MathFunctions::Linf(_S->at(i)->at(0).getRoundedVector(), dimensions) * _S->at(i)->size()?
            float f2 = 0;
            for(unsigned j=0; j<_S->at(i)->size(); j++) {
                v += MathFunctions::Linf(_S->at(i)->at(j).getRoundedVector(), dimensions);
            }
            // Add C2 to model 
            model.add(IloSum(y[i]) >= v);
            v.end();
        }
        // TODO: Check this
        // C3
        for(IloInt i=0; ((unsigned) i)<f_vec.size(); i++) {
            for(IloInt d=0; d<dimensions; d++) {
                IloExpr v(env);
                IloExpr w(env);
                for(unsigned tsm=0; tsm<_S->size(); tsm++) {
                    // get vector t, S[tsm][0] is never empty
                    // pick the first vector, they are all identical
                    float* t = _S->at(tsm)->at(0).getRoundedVector();
                    float aLinf = MathFunctions::Linf(t, dimensions);
                    if(aLinf != 0) {
                        //MathFunctions::VecDiv(t, aLinf, dimensions); // Should be wrong. Why should we divide t by aLinf twice??
                        w += y[tsm][i] * t[d] / aLinf;
                    } else {
                        w += 0; // TODO: What does that do?
                    }
                }
                for(unsigned c=0; c<_configs->size(); c++) {
                    // Add x_C to sum if configuration C has current f as its profile
                    if(MathFunctions::areVectorsEqual(_configs->at(c).getF(), f_vec[i], dimensions)) {
                        v += x[c];
                    }
                }
                model.add(w <= f_vec[i][d] * v);
                v.end();
                w.end();
            }
        }
        // Objective: minimize the sum of machines used
        IloObjective objective = IloMinimize(env, IloSum(x));
        model.add(objective);
        // Solve problem (solution stored in cplex)
        IloCplex cplex(model);
        cplex.setOut(env.getNullStream());
        if(!cplex.solve()) { 
            env.end();
            return false;
        }
        // Retrieve indices of chosen configurations (from cplex)
        std::vector<int> chosen_configs;
        for(unsigned i=0; i<_configs->size(); i++) {
            int tmp = cplex.getValue(x[i]);
            for(int j=0; j<tmp; j++) {
                chosen_configs.push_back(i);
            }
        }
        // Check if there are enough machines
        if(chosen_configs.size() > _machines->size()) {
            return false;
        }
        bool successful;
        if(randomized) {
            // TODO: What is mach_f?
            std::vector<float*> mach_f;
            // Assign large jobs (according to chosen configurations)
            assignLargeJobsRandomized(_configs, _L, _machines, &chosen_configs, &mach_f, _epsilon, dimensions);
            // Assign small jobs (according to cplex solution) which can fail
            successful = assignSmallJobsRandomized(&cplex, &y, _S, _machines, &f_vec, &mach_f, dimensions);
            // Clean up
            for(int i=0; i<mach_f.size(); i++) {
                delete [] mach_f[i];
            }
        } else {
            // Vector with entry for each machine storing the identifier of its profile
            // It initializes every entry with the index of the last configuration (which is the one with no jobs)
            std::vector<int> machine_f(static_cast<int>(_machines->size()), f_conf_indx.back());
            // Assign large jobs (according to chosen configurations)
            assignLargeJobsDeterministic(_configs, _L, _machines, &f_conf_indx, &chosen_configs, &machine_f);
            // Assign small jobs (according to cplex solution) which can fail
            successful = assignSmallJobsDeterministic(&cplex, &y, _S, _machines, &machine_f, _epsilon, dimensions);
        }
        env.end();
        return successful;
    }
};

#endif
