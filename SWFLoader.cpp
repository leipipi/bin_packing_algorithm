
#include <stdio.h>
#include <unistd.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <random>

#include "Job.h"
#include "Environment.h"

class SWFParser {
public:

    SWFParser() {
    }

    std::vector < std::vector < float >* >* parse(const char* _path, int _num) {
        std::vector < float > machine = std::vector < float > ();

        std::vector < int > fields = std::vector <int> ();
        std::vector < std::vector < float >* >* queue = new std::vector < std::vector < float >* >();

        int count = 0;
        int linecount = 0;
        std::ifstream file(_path);
        if(!file.is_open()) {
            fprintf(stderr, "ERROR: Cannor open job file: %s - exit\n", _path);
            exit(1);
        }
        std::string line;

        while(std::getline(file, line)) {
            linecount++;

            std::istringstream iss(line);
            std::string word;

            if(line[0] == '#') {
                int count = 1;
                while(iss >> word) {
                    if(word != "#") {
                        std::replace(word.begin(), word.end(), '|', ' ');
                        std::istringstream iss2(word);
                        if(count == 1) {
                            int temp;
                            while(iss2 >> temp)
                                fields.push_back(temp);
                        }
                        if(count == 2) {
                            float temp;
                            while(iss2 >> temp)
                                machine.push_back(temp);
                        }
                        ++count;
                    }
                }
                continue;
            }

            // skip the header
            if(line[0] == ';')
                continue;

            // check if demanded amount of vectors is exceeded
            if(count < _num) {
                // copy line data to vector
                std::vector <float> line_vec = std::vector <float> ();
                while(iss >> word) {
                    line_vec.push_back(std::stod(word.c_str()));
                }

                // get relevant entries and normalize
                std::vector <float>* vector = new std::vector <float> ();
                for(unsigned i = 0; i < fields.size(); ++i) {
                    double tmp_data = line_vec[fields[i] - 1] / machine[i];
                    if(tmp_data <= 0.0 || tmp_data > 1.0)
                        goto next;
                    vector->push_back(tmp_data);
                }
                queue->push_back(vector);
                count++;
            } else {
                break;
            }
next:
            ;
        }
        file.close();
        return queue;
    }


};

extern "C"
std::vector <Job>* loadJobs(const char* _path, JobMetaData* _jobdata, Environment* _env, int _runs) {
    std::vector <Job>* jobVec = new std::vector <Job>();

#ifdef RANDOM
    std::random_device device;
    std::mt19937 generator(device());
#else
    char* seed = getenv("RAND_SEED");
    std::mt19937 generator;
    if(!seed)
        generator = std::mt19937(1234);
    else
        generator = std::mt19937(atoi(seed));
#endif

    std::uniform_int_distribution<int> distBurst(_jobdata->lowerBurst, _jobdata->upperBurst);
    std::uniform_int_distribution<int> distArrive(0, _jobdata->duration);

    SWFParser parser = SWFParser();

    // generate factor* m_num vectors from swf file
    std::vector < std::vector < float >* >* vectors = parser.parse(_path, _runs * _jobdata->numberOfJobs);
    printf("vectors size1: %lu\n", vectors->size());
    printf("vectors size2: %lu\n", vectors->at(0)->size());

    // use vectors to generate jobs
    for(unsigned i = 0; i < vectors->size(); i++) {
        jobVec->push_back(Job(i + 1, vectors->at(i), distBurst(generator), distArrive(generator)));
    }
    for(unsigned i = 0; i < vectors->size(); i++) {
        delete vectors->at(i);
    }
    delete vectors;

    return jobVec;
}

