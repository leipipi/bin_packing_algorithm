#include "ConfigReader.h"
#include <cstdio>
#include <fstream>
#include <sstream>

ConfigReader::ConfigReader(const char* _path) {
    env.loader = "";
    env.algorithm = "";
    env.online = false;
    env.stepWidth = 1;
    env.testType = 'm';
    jmd.dimensions = 3;
    env.epsilon = 0.01;
    jmd.numberOfJobs = 1000;
    env.machines = 10;
    env.machineSize = 1.0;
    env.runs = 10;
    env.startJobs = 1;
    env.file = "";
    env.cpus = 1;
    env.ram = 1;
	env.maxLocalDisk = -1;
	env.maxNetwork = -1;
    env.randomJobSelection = false;
    jmd.lowerBound = 0.1;
    jmd.upperBound = 0.3;
    jmd.lowerBurst = 500;
    jmd.upperBurst = 3000;
    jmd.arrivalSpan = 1500;
	env.intervalMin = 0.01;
	env.intervalMax = 0.25;

    readConfigFile(_path);
}

void ConfigReader::readConfigFile(const char* _path) {
    std::ifstream file(_path);
    if(!file.is_open()) {
        fprintf(stderr, "ERROR: Cannot open config file %s - exit\n", _path);
        exit(1);
    }
    std::string line;

    while(std::getline(file, line)) {
        std::istringstream iss(line);
        if(line.empty() || line[0] == '#')
            continue;
        std::string word1;
        std::string word2;

        iss >> word1;
        iss >> word2;

        if(word1 == "loader") {
            env.loader = word2;
        } else if(word1 == "algorithm") {
            env.algorithm = word2;
        } else if(word1 == "online") {
            if(word2 == "true")
                env.online = true;
            else
                env.online = false;
        } else if(word1 == "random_job_selection") {
            if(word2 == "true")
                env.randomJobSelection = true;
            else
                env.randomJobSelection = false;
        } else if(word1 == "step_width") {
            env.stepWidth = atoi(word2.c_str());
        } else if(word1 == "max_ram") {
            env.ram = atoi(word2.c_str());
        } else if(word1 == "max_cpus") {
            env.cpus = atoi(word2.c_str());
        } else if(word1 == "max_local_disk") {
            env.maxLocalDisk = atoi(word2.c_str());
        } else if(word1 == "max_network") {
            env.maxNetwork = atoi(word2.c_str());
        } else if(word1 == "test_type") {
            env.testType = word2[0];
        } else if(word1 == "interval_min") {
            env.intervalMin = atof(word2.c_str());
        } else if(word1 == "interval_max") {
            env.intervalMax = atof(word2.c_str());
        } else if(word1 == "dimensions") {
            jmd.dimensions = atoi(word2.c_str());
        } else if(word1 == "epsilon") {
            env.epsilon = atof(word2.c_str());
        } else if(word1 == "jobs") {
            jmd.numberOfJobs = atoi(word2.c_str());
        } else if(word1 == "machines") {
            env.machines = atoi(word2.c_str());
        } else if(word1 == "machine_size") {
            env.machineSize = atof(word2.c_str());
        } else if(word1 == "runs") {
            env.runs = atoi(word2.c_str());
        } else if(word1 == "start_jobs") {
            env.startJobs = atoi(word2.c_str());
        } else if(word1 == "file") {
            env.file = word2;
        } else if(word1 == "lower_bound") {
            jmd.lowerBound = atof(word2.c_str());
        } else if(word1 == "upper_bound") {
            jmd.upperBound = atof(word2.c_str());
        } else if(word1 == "lower_burst") {
            jmd.lowerBurst = atoi(word2.c_str());
        } else if(word1 == "upper_burst") {
            jmd.upperBurst = atoi(word2.c_str());
        } else if(word1 == "arrival_span") {
            jmd.arrivalSpan = atoi(word2.c_str());
        } else {
            fprintf(stderr, "Error: Unknown option: %s - exit\n", word1.c_str());
            exit(1);
        }
    }
}

Environment* ConfigReader::getEnvironment() {
    return &env;
}

JobMetaData* ConfigReader::getJobMetaData() {
    return &jmd;
}

