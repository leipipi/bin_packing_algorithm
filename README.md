####Introduction:

This is the simulator developed to run experiments with vector scheduling and vector bin packing problems. 
Existing codes for algorithms are implemented for the paper "Analysis of Heuristics for Vector Scheduling and Vector Bin Packing". 
Users can implement their own algorithms and test them in the simulator.
Several useful .h files could be helpful for users to write their own algorithms:
Machine.h: Functions to create bins or machines.
Mathfunctions.h: Functions for doing mathematical operations with vectors.
Sorting.h: Various sorting methods for vectors are defined.

####Set Up:

To build the simulator, one is encouraged to run either ./buildwithoutArcFlow.sh or ./buildwithArcFlow.sh.
To run the first bash script, Gurobi and Cplex are required to build the exact and approximation algorithms.
There are no specific requirements for running the second script.

####Useage:

After the simulator has been built into the system, one could test it by running the following:
./VSSimKernel config
to run the PaniDP algorithm with a pre-defined setting in the configuration file.
One can also run the TestRunner.sh defined in the TestRunner folder to run multiple algorithms together.
The TestRunner.sh will generate the config file and run the VSSimKernel automatically based on the defined settings. 
The default config file and TestRunner.sh are commented on with explanations for different parameters that should be helpful for understanding.

####Data Set:

With the default CXLoader, job vectors can be generated uniformly in the given range.
Users only need to provide the configuration file's lower and upper bound numbers to define the range
