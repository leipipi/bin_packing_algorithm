#include <stdio.h>
#include <unistd.h>
#include <vector>

#include "Job.h"

extern "C"
std::vector <Job>* loadJobs(const char* _path, JobMetaData* _jobdata, int _runs) {
	std::vector <Job>* jobVec = 0;
	printf ("start job loading ...\n");
	sleep(1);
	printf ("jobs loaded.\n");

	return jobVec;
}
