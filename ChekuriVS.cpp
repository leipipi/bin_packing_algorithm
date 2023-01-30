#include "Machine.h"
#include "MathFunctions.h"
#include <vector>
#include <algorithm>

#include <ilcplex/ilocplex.h>
#include <thread>

// typedef for CPLEX matrix variables
typedef IloArray<IloNumArray> FloatMatrix;
typedef IloArray<IloNumArray> Float3DMatrix;
typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<NumVarMatrix> NumVar3DMatrix;

/**
 * Vector b is added to vector a.
 * @param a first vector.
 * @param b second vector.
 * @param _dim dimension of vectors.
 */
void vecAdd(float* a, float* b, int _dim) {
    for(unsigned i = 0; i < _dim; i++) {
        a[i] += b[i];
    }
}

/**
 * Vector b is subtracted from vector a, unless a would have negative components.
 * @param a first vector.
 * @param b second vector.
 * @param _dim dimension of vectors.
 * @return true if and if subtraction was done.
 */
bool vecDifPosCheck(float* a, const float* b, int _dim) {
    for(unsigned i = 0; i < _dim; i++) {
        if(a[i] - b[i] < 0) {
            return false;
        }
    }

    // TODO: Mit einem Job funktioniert es nicht, sonst schon! Und es werden die gerundeten Jobs in den Maschinen addiert.
    // TODO: Wenn man die naechsten drei Zeilen entfernt, scheint es immer zu funktionieren, aber die normalen Jobs werden in den Maschinen addiert. 
    // TODO: Das Problem ist irgendwo in PackL.

    for(unsigned i = 0; i < _dim; i++) {
        a[i] -= b[i];
    }
    return true;
}

bool packS(std::vector<Job> jobs, std::vector<Machine> &machines, int _dim) {
    // declare CPLEX environment and model
    IloEnv env;
    IloModel model(env);

    // set up values of job vector coordinates
    FloatMatrix p(env, jobs.size());
    for(IloInt i = 0; ((unsigned) i) < jobs.size(); i++) {
        float* tmp_vec = jobs[i].getKPIVec();
        p[i] = IloNumArray(env, _dim);
        for(IloInt j = 0; j < _dim; j++) {
            p[i][j] = tmp_vec[j];
        }
    }

    // set up values of machine vector coordinates
    FloatMatrix b(env, machines.size());
    for(IloInt i = 0; ((unsigned) i) < machines.size(); i++) {
        float* tmp_vec = machines[i].getRemainingCapacity();
        b[i] = IloNumArray(env, _dim);
        for(IloInt j = 0; j < _dim; j++) {
            b[i][j] = tmp_vec[j];
        }
    }

    // set up matrix for x_ij variables
    NumVarMatrix x(env, jobs.size());
    for(IloInt i = 0; ((unsigned) i) < jobs.size(); i++) {
        x[i] = IloNumVarArray(env, machines.size(), 0, IloInfinity, ILOFLOAT);
    }

    // C1
    for(IloInt i = 0; ((unsigned) i) < jobs.size(); i++) {
        model.add(IloSum(x[i]) == 1);
    }

    // C2
    for(IloInt j = 0; ((unsigned) j) < machines.size(); j++) {
        for(IloInt k = 0; k < _dim; k++) {
            IloExpr v(env);
            for(IloInt i = 0; ((unsigned) i) < jobs.size(); i++) {
                v += x[i][j] * p[i][k];
            }
            model.add(v <= b[j][k]);
            v.end();
        }
    }

    // optimize the problem and obtain solution.
    IloCplex cplex(model);
    if(!cplex.solve()) {
        env.end();
        return false;
    }

    // assign jobs to machines depending on results of x
    for(IloInt i = (jobs.size() - 1); i >= 0; i--) {
        for(IloInt j = 0; ((unsigned) j) < machines.size(); j++) {
            double x_val = cplex.getValue(x[i][j]);
            if(x_val == 1.0) {
                if(machines[j].assignJob(&(jobs[i]))) {
                    jobs.erase(jobs.begin() + i);
                } else {
                    return false;
                }
            }
        }
    }

    // if there are still jobs left, devide jobs into 
    // m partitions and assign to machines
    if(!jobs.empty()) {
        int set_size = std::ceil(static_cast<double> (jobs.size()) / static_cast<double> (machines.size()));
        if(set_size > _dim) {
            return false;
        }
        int mach_count = 0;
        while(!jobs.empty()) {
            for(int i = 0; i < set_size; i++) {
                machines[mach_count].assignJobOverload(&(jobs.back()));
                jobs.pop_back();
                if(jobs.empty()) {
                    break;
                }
            }
            mach_count++;
        }
    }

    env.end();
    return true;
}

void feasPackRecursion(std::vector<int> &types, int rand_left, int rand_right, std::vector<std::vector<Job>> &jobs, std::vector<std::vector<int>> &feas_pack, float* machine_vec, int _dim) {
    if(rand_right - rand_left == -1) {
        // check if selected jobs would fit into machine of height "mach_vec"
        float sum_job[_dim] = {};

        std::vector<int> types_copy = types;
        for(unsigned i = 0; i < types_copy.size(); i++) {
            while(types_copy[i] > 0) {
                vecAdd(sum_job, jobs[i][0].getRoundedVector(), _dim);
                --types_copy[i];
            }
        }
        if(vecDifPosCheck(machine_vec, sum_job, _dim)) {
            feas_pack.push_back(types);
        }
        return;
    }
    while(types[rand_left] >= 1) {
        feasPackRecursion(types, rand_left + 1, rand_right, jobs, feas_pack, machine_vec, _dim);
        --types[rand_left];
        if(types[rand_left] == 0) {
            return;
        }
    }
    if(types[rand_left] == 0) {
        feasPackRecursion(types, rand_left + 1, rand_right, jobs, feas_pack, machine_vec, _dim);
    }
}

/*
"jobs" contains jobs listed by their job type. For exmple jobs[0] will return a 
JQUEUE of jobs that are from the first Type. 

generateFeasPack returns a 2d vector that contains vectors which point out
what amount of vectors of each type will fit into a machine.

For example: Let vec be a std::vector<std::vector<int>> and vec[0] looks like this (1,2,0,5).
This means that a set consisting of one job of type 1, two of type 2
and five of type 4 will fit into a machine.
 */
std::vector<std::vector<int>> generateFeasPack(std::vector<std::vector<Job>> &jobs, float* machine_vec, int _dim) {
    // "feas_pack" contains vectors that describe all possible 
    // subsets of "jobs" that would fit in a machine of height "machine_vec".
    std::vector<std::vector<int>> feas_pack;
    for(int i = jobs.size(); i > 0; --i) {
        // use macro string to iterate over all combinations of job types
        std::string s(i, '1');
        s.resize(jobs.size(), '0');
        do {
            //cout << "Str: " << s << endl;
            // types contains the number of vectors of each type for the current iteration
            std::vector<int> types(jobs.size(), 0);
            int rand_left = 0;
            int rand_right = 0;
            for(unsigned j = 0; j < jobs.size(); j++) {
                if(s[j] == '1') {
                    if(jobs[j].empty()) {
                        goto next_generateFeasPack;
                    }
                    if(rand_left == 0 && j != 0) {
                        rand_left = j;
                    }
                    if(j > ((unsigned) rand_right)) {
                        rand_right = j;
                    }
                    types[j] = jobs[j].size();
                }
            }
            // iterate over all possible vector amounts of the chosen types
            feasPackRecursion(types, rand_left, rand_right, jobs, feas_pack, machine_vec, _dim);
next_generateFeasPack:
            ;
        } while(std::prev_permutation(s.begin(), s.end()));
    }
    return feas_pack;
}

/*
Try to pack all jobs from "jobs" into the machines from "machines".
"step" is the index of the machine that is currently been viewed.
 */
bool packL(std::vector<std::vector<Job>> &jobs, std::vector<Machine> &machines, int step, int _dim) {
    // If all jobs are packed, return true
    for(unsigned i = 0; i < jobs.size(); i++) {
        if(!jobs[i].empty()) {
            break;
        }
        if(i == jobs.size() - 1) {
            return true;
        }
    }
    // If all machines are full, return false 
    if(machines.size() == step) {
        return false;
    }
    // Backup jobs and machine
    std::vector<std::vector<Job>> jobs_backup = jobs;
    Machine machine_backup = machines[step];
    
    // iterate over all combinations of job sets that will fit into a machine
    std::vector<std::vector<int>> feas_pack = generateFeasPack(jobs, machines[step].getRemainingCapacity(), _dim);
    for(unsigned j = 0; j < feas_pack.size(); j++) {
        for(unsigned k = 0; k < feas_pack[j].size(); k++) {
            while(feas_pack[j][k] > 0) {
                machines[step].assignJob(&(jobs[k].back()));
                jobs[k].pop_back();
                --feas_pack[j][k];
            }
        }

        // try next machine
        if(packL(jobs, machines, step + 1, _dim)) {
            return true;
        } else {
            machines[step] = machine_backup;
            jobs = jobs_backup;
        }
    }
    return false;
}

/**
 * Recursive function that creates all capacity configurations.
 * @param vec partial capacity configuration; this function creates all possible continuations.
 * @param numbers numbers used for creating capacity configurations.
 * @param vectors vector of all capacity configurations.
 * @param mq_size number of machines.
 * @param _eps bins have size 1+eps instead of 1.
 * @param index current vector index.
 */
void generateCapacityVectors(std::vector<float> vec, std::vector<int> numbers, std::vector<std::vector<float>> &vectors, int mq_size, float _eps, int index) {
    // generate m many of each capacity configuration vector
    // this eases the generation of the machine configurations later
    if(((unsigned) index) == vec.size()) {
        for(int i = 0; i < mq_size; i++) {
            vectors.push_back(vec);
        }
        return;
    }

    for(unsigned i = 0; i < numbers.size(); i++) {
        vec[index] = numbers[i] * _eps;
        generateCapacityVectors(vec, numbers, vectors, mq_size, _eps, index + 1);
    }
}

/**
 * Every distinct integer tuple (a_1, a_2, ..., a_d) where 0 <= ai <= ceil(1/eps),
 * is a capacity configuration. This function creates all of them and each mq_size
 * times. For two machines and eps = 0.5, the output would be {{1, 1, 1}, {1, 1, 1},
 * {1, 1, 0.5}, {1, 1, 0.5}, ..., {0.5, 1, 1}, {0.5, 1, 1}, ..., {0, 0, 0.5}, {0, 0, 0.5},
 * {0, 0, 0}, {0, 0, 0}}.
 * 
 * @param mq_size number of machines.
 * @param _dim dimension of job vectors.
 * @param _eps bins have size 1+eps instead of 1.
 * @return capacity configurations. 
 */
std::vector<std::vector<float>> generateCapacityConfigs(int mq_size, int _dim, float _eps) {
    std::vector<std::vector<float>> vecs;
    std::vector<int> numbers;

    // numbers = {ceil(1.0/eps), ..., 2, 1, 0}; e.g., for eps = 0.5, numbers = {2, 1, 0}
    for(int i = std::ceil(1.0 / _eps); i >= 0; i--) {
        numbers.push_back(i);
    }
    std::vector<float> vec(_dim, 0.0);
    generateCapacityVectors(vec, numbers, vecs, mq_size, _eps, 0);

    return vecs;
}

/**
 * Function for debugging. It prints a vector of capacity configurations.
 * @param cc vector of vectors of floats.
 */
void print(std::vector<std::vector<float>> cc) {
    for(std::vector<std::vector<float>>::iterator it = cc.begin(); it != cc.end(); it++) {
        for(std::vector<float>::iterator a = it->begin(); a != it->end(); a++) {
            printf("  %4.2f", *a);
        }
        printf("\n");
    }
}

bool checkL_(std::vector<std::vector<Job>> L_) {
    int i = 0;
    for(std::vector<std::vector<Job>>::iterator it = L_.begin(); it != L_.end(); it++) {
        printf("CHECK %d ", i++);
        for(std::vector<Job>::iterator a = it->begin(); a != it->end(); a++) {
            a->print();
        }
        printf("\n");
    }
    return true;
}


/**
 * PTAS-Algorithm by Chekuri and Khanna for Vector Scheduling.
 * In Nikolay's thesis, this algorithm is called PTASChK.
 * Like in LTAA, the bins have a size 1+epsilon instead of 1.
 * @param counter
 * @param out_machines
 * @param jq
 * @param mq
 * @return 
 */
extern "C"
int planJobs(int &counter, std::vector<Machine> **out_machines, std::vector<Job> *jq, std::vector<Machine>* mq, float _eps) {

    // Dimension of jobs.
    int d = (jq->at(0)).getDimensions();
    // delta is used for preprocessing the jobs / vectors p_j = (p_j^1, ..., p_j^d). If p_j^i <= delta * ||p_j||, then p_j^i := 0.
    float delta = _eps / d;
    // Set of small jobs.
    std::vector<Job> S;
    // Set of sets of large jobs with vectors rounded. (L_ is L' in the paper.)
    std::vector<std::vector<Job>> L_;
    // TODO: out_machines are returned as empty vector
    *out_machines = new std::vector<Machine>();

    // Preprocessing: round coordinates of jobs and partition jobs in large and small ones
    for(std::vector<Job>::iterator it = jq->begin(); it != jq->end(); it++) {
        (*it).roundChekuri(delta, _eps);
        if((*it).isLarge()) {
            // store large jobs in a list depending on their types
            if(L_.empty()) {
                L_.push_back(std::vector<Job>(1, (*it)));
            } else {
                for(unsigned j = 0; j < L_.size(); j++) {
                    if(L_[j][0].getRoundedVector() == (*it).getRoundedVector()) {
                        L_[j].push_back(*it);
                        break;
                    }
                    if(L_.size() - j == 1) {
                        L_.push_back(std::vector<Job>(1, (*it)));
                        break;
                    }
                }
            }
        } else {
            S.push_back((*it));
        }
    }

    // Create all capacity configurations and each of them mq_size many times
    std::vector<std::vector<float>> capacity_configs = generateCapacityConfigs(mq->size(), d, _eps);
    // Create string of form "110000000000000000000000000000000000000000000000000000" 
    // of length capacity_configs.size() starting with as many 1's as there are machines
    std::string s(mq->size(), '1');
    s.resize(capacity_configs.size(), '0');
    // Iterate over all machine configurations and try to pack S and L_ in one of them
    do {
        int machine_index = 0;
        // The string has as many characters as there are capacity configurations. If the character is "1",
        // the respective capacity configurations is used. The number of 1's in the string equals the number
        // of machines. The ith capacity configuration used is assigned to the ith machine.
        for(unsigned i=0; i<capacity_configs.size(); i++) {
            if(s[i] == '1') {
                mq->at(machine_index).setVector(capacity_configs[i]);
                machine_index++;
            }
        }
        
        if(!L_.empty()) {
//            checkL_(L_);
            if(!packL(L_, *mq, 0, d)) {
                // If packing the large jobs fails, continue with the next combination of capacity configurations.
                continue;
            }
            // If packing the large jobs is successful, continue with small jobs.
        } else {
            if(packS(S, *mq, d)) {
                // If there no large jobs and packing the small jobs is successful, return true.
                // TODO: *out_machines = new std::vector<Machine>(*mq);
                return true;
            } else {
                // If there no large jobs and packing the small jobs fails, increase fail counter and return false.
                break;
            }
        }
        if(!S.empty()) {
            if(packS(S, *mq, d)) {
                // If packing the large and small jobs is successful, return true.
                // TODO: *out_machines = new std::vector<Machine>(*mq);
                return true;
            }
            // If packing the small jobs fails, continue with next combination of capacity configurations.
        } else {
            // If packing the large jobs is successful and there are no small jobs, return true.
            // TODO: *out_machines = new std::vector<Machine>(*mq);
            return true;
        }
    } while(std::prev_permutation(s.begin(), s.end()));
    // If packing fails for every configuration, increase fail counter and return false.
    counter++;
    return false;
}




