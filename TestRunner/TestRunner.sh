#/bin/bash

##### Data loader
LOADER=../CXLoader.so    	# ../C0Loader.so ../C6Loader.so ../CXLoader.so ../SWFLoader.so ../GWFLoader.so
## required for CXLoader
INTERVAL_MIN=0.01
INTERVAL_MAX=0.25
## required for SWFLoader and GWFLoader
JOB_FILE=../Datasets/anon_jobs.gwf

##### Algorithm
ALGO=(../CNS.so ../KouVS_BFD_Lexi.so ../KouVS_BFD_LInf.so ../KouVS_BFD_Sum.so ../KouVS_BF.so ../KouVS_FFD_Lexi.so ../KouVS_FFD_LInf.so ../KouVS_FFD_Sum.so ../KouVS_FF.so ../PaniDP.so ../PaniL1.so ../PaniL2.so ../PaniLinf.so ../LS.so ../HybridDP.so ../SE.so)
ONLINE=false				# NOT USED YET
STEP_WIDTH=2
TESTTYPE=maxMachineLoad		# NOT USED YET
DIMENSIONS=3
EPSILON=0.1
RUNS=5
START_JOB=1
RANDOM_JOB_SELECTION=true

##### Pseudo-random number generator seed
SEED=4321

##### Jobs and environment
JOBS=200
MACHINES=20
MACHINE_SIZE=1.0
MAX_CPUS=-1
MAX_RAM=-1
MAX_LDS=-1
MAX_NET=-1

### Parameters for online algorithms
LOWERBOUND=0.1
UPPERBOUND=0.3
LOWERBURST=500
UPPERBURST=3000
ARRIVAL_SPAN=1500




for i in ${ALGO[@]}
do
	echo loader $LOADER > config.txt
	echo algorithm $i  >> config.txt
	echo online $ONLINE >> config.txt 
	echo step_width $STEP_WIDTH >> config.txt 
	echo test_type $TEST_TYPE >> config.txt 
	echo dimensions $DIMENSIONS >> config.txt 
	echo epsilon $EPSILON >> config.txt 
	echo jobs $JOBS >> config.txt 
	echo machines $MACHINES >> config.txt 
	echo machine_size $MACHINE_SIZE >> config.txt 
	echo runs $RUNS >> config.txt 
	echo start_jobs $START_JOB >> config.txt 
	echo file $JOB_FILE >> config.txt 
	echo lower_bound $LOWER_BOUND >> config.txt 
	echo upper_bound $UPPER_BOUND >> config.txt 
	echo lower_burst $LOWER_BURST >> config.txt 
	echo upper_burst $UPPER_BURST >> config.txt 
	echo arrival_span $ARRIVAL_SPAN >> config.txt 
	echo random_job_selection $RANDOM_JOB_SELECTION >> config.txt 
	echo max_ram $MAX_RAM >> config.txt 
	echo max_cpus $MAX_CPUS >> config.txt 
	echo max_local_disk $MAX_LDS >> config.txt 
	echo max_network $MAXNET >> config.txt 
	echo interval_min $INTERVAL_MIN >> config.txt 
	echo interval_max $INTERVAL_MAX >> config.txt 
	echo Currently running: $ALGO
	RAND_SEED=$SEED ../VSSimKernel config.txt
done
