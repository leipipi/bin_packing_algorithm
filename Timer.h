#ifndef TIMER_H_DEF
#define TIMER_H_DEF

#include <sys/time.h>


class Timer {
public:
    Timer();
    ~Timer();

    void   start();
    void   stop();
    double getElapsedTime();
    double getElapsedTimeInSec();
    double getElapsedTimeInMilliSec();
    double getElapsedTimeInMicroSec();

private:
    double startTimeInMicroSec;
    double endTimeInMicroSec;
    int    stopped;
    timeval startCount;
    timeval endCount;
};

#endif
