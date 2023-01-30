#ifndef CONFIGREADER_H_
#define CONFIGREADER_H_

#include "Job.h"
#include "Environment.h"
#include <string>

class ConfigReader {
public:
    ConfigReader(const char* _path);

    Environment* getEnvironment();
    JobMetaData* getJobMetaData();
private:
    void readConfigFile(const char* _path);
    Environment env;
    JobMetaData jmd;

};
#endif

