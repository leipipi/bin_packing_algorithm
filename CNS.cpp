 #include <unistd.h>
#include <chrono>
#include "Job.h"
#include "Machine.h"
#include "MathFunctions.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <iostream>
/**
 * CNS by BV16.
 */
int limit = 20;
int WTF=0;
int totaljob;
int WTFID=-1;
//#define DEBUG 0
float strangM;
float strangTC;
int iter = 0;
int **m;
int Swapiter = 0;
int Totaliter =0;
std::vector<Job> * TCP = new std::vector <Job>();
std::vector<Job> * TC_tmp = new std::vector <Job>();
std::vector<Job> * mq_tmp = new std::vector <Job>();
std::vector<Job> * TCP_tmp = new std::vector <Job>();
std::vector<Job> * S = new std::vector <Job>();
std::vector<Job> * S_tmp = new std::vector <Job>();
std::vector<Job> * T = new std::vector <Job>();
std::vector<Job> * TT = new std::vector <Job>();
std::vector<Job> * P = new std::vector <Job>();
std::vector<Job> * Q = new std::vector <Job>();
std::vector<Job> S_permant;
	template <class T>
int getArrayLen(T& array)
{
	return (sizeof(array) / sizeof(array[0]));
}
void cleanJobVector(std::vector<Job>*S){
	for(int i = S->size()-1; i >= 0; i--){
		S->erase(S->begin()+i);
	}
	S->shrink_to_fit();
}
float Max(float* a, int dimens)
{
	float max = 0.;
	for(unsigned i = 0 ; i < dimens; i++)
	{
		if(a[i]>max)
			max = a[i];
	}
	return max;
}
float Min(float* a, int dimens)
{
	float min = 100.;
	for(unsigned i = 0 ; i < dimens; i++)
	{
		if(a[i]<min)
			min = a[i];
	}
	return min;
}
float WeightofPair(Machine * MQ, int * a, int dimens)
{
	float* tmp_jb;
	float* tmp_sum = new float[dimens];
	float weight = 0;
	tmp_jb = MQ->assignedJobs->at(a[0]).getKPIVec();
	for(unsigned i = 0; i<dimens; i++)
	{
		tmp_sum[i] = tmp_jb[i];
	}
	if(a[1]>0)
	{
		tmp_jb = MQ->assignedJobs->at(a[1]).getKPIVec();
		MathFunctions::VectorAddition(tmp_sum ,tmp_jb,dimens);
	}
	weight = Max(tmp_sum, dimens);
	delete []tmp_sum;
	return weight;
}
int SizeofItems(int * a)
{
	int size;
	if(a[1]>0)
		size = 2;
	else
		size = 1;
	if(a[0] ==-1 && a[1] == -1)
		size = 0;
	return size;
}
int SizeofPair(int size)
{
	size = size*(size-1)/2 + size;
	return size;
}
int ** GeneratePair(int size)
{
	int inisize = size;
	size = size*(size-1)/2 + size;
	delete [] m;
	int **m=new int*[size];
	int first = 0;
	int second = -1;
	for(unsigned i = 0; i < size; ++i)
	{
		m[i] = new int[2];
		if(first == inisize)
		{
			second = 1;
			first = 0;
			m[i][0] = first;
			m[i][1] = second;
			++second;
			continue;
		}
		if(second == -1){
			m[i][0] = first;
			m[i][1] = second;
			++first;
		}
		else
		{
			if(second < inisize)
			{
				m[i][0] = first;
				m[i][1] = second;
				++second;
			}
			else
			{
				++first;
				second = first + 1;
				m[i][0] = first;
				m[i][1] = second;
				++second;
			}
		}
	}
	return m;
}
//set everything in tabu as -1
int ** resetTabu(std::vector<Machine>*mq)
{
	iter = 0;
	int  size;
	size = mq->size();
	int **m=new int*[size];
	for(unsigned i = 0; i < size; i++ )
	{
		m[i] = new int[mq->at(i).assignedJobs->size()];
		for(unsigned j = 0; j < mq->at(i).assignedJobs->size(); j++)
		{
			m[i][j] = -1;
		}
	}
	return m;
}
bool isnotTabu(int * S , int * Tabu)
{
	if(Tabu[S[0]] >=iter || (S[1]>0 && Tabu[S[1]]>iter))
	{
		return false;
	}
	return true;
}
//simplized tabulist
int * updateTabu(int * S , int * T, int * Tabu, Machine m)
{
	iter = iter +1;
	int sizeofTabu = getArrayLen(Tabu);
	std::vector<int> tmp_Tabulist;
	std::vector<int> tmp_Tabu;
	for(unsigned i = 0; i < sizeofTabu; i++)
	{
		if(Tabu[i]!=0)
		{
			tmp_Tabulist.push_back(i);
			tmp_Tabu.push_back(Tabu[i]);
		}
		Tabu[i] = 0;
	}
	int sizeoft = SizeofItems(T);
	int sizeofs = SizeofItems(S);
	int sizeofm = m.assignedJobs->size();
	int dimensions = m.getDimensions();
	// update same item
	for(unsigned i = 0; i < sizeofm; i++ )
	{
		if(MathFunctions::areVectorsEqual(m.assignedJobs->at(i).getKPIVec(),m.assignedJobs->back().getKPIVec(),dimensions)||(sizeoft == 2 && MathFunctions::areVectorsEqual(m.assignedJobs->at(i).getKPIVec(),m.assignedJobs->at(sizeofm-2).getKPIVec(),dimensions)))
			Tabu[i] = iter + 1;
	}
	// update exsit tabu item
	for(unsigned i = 0; i < tmp_Tabulist.size(); i++)
	{
		if(sizeofs == 2 && tmp_Tabulist[i] > S[1])
			Tabu[tmp_Tabulist[i]-2] += tmp_Tabu[i];
		else if(tmp_Tabulist[i]>S[0])
			Tabu[tmp_Tabulist[i]-1] += tmp_Tabu[i];
		else if(tmp_Tabulist[i]<S[0])
			Tabu[tmp_Tabulist[i]] += tmp_Tabu[i];
	}
	return Tabu;
}
float Distance(float*tmp,int dimens)
{

	float total=0.;
	for (unsigned i = 0; i < dimens; i++){
		total +=pow(tmp[i],2);
	}
	total = sqrt(total);
	return total;
}
//packset
int iteration = 0;
float distance(std::vector<Job>*S,int dimens)
{
	float * tmp;
	float tmp_S[dimens] = {0.};
	float * tmp_jv;
	float total=0;
	if(S->size()>0){
		tmp = S->at(0).getKPIVec();
		for(int i = 0; i<dimens ; i++)
		{
			tmp_S[i] += tmp[i];
		}
		for(unsigned i = 1; i < S->size(); i++)
		{
			tmp_jv = S->at(i).getKPIVec();
			MathFunctions::VectorAddition(tmp_S,tmp_jv,dimens);
		}
		for (unsigned i = 0; i < dimens; i++){
			total +=pow(tmp_S[i],2);
		}
		total = sqrt(total);
	}
	return total;
}
//packset
bool compare(float * tmp, int i, int dimens)
{
	for(unsigned j = 0; j < dimens; j++)
	{
		if(tmp[j]>i)
			return false;
	}
	return true;
}
float* loadV(std::vector<Job>* S,int dimens)
{
	float * tmp;
	float * tmp_r = new float[dimens];
	for(unsigned i = 0; i < dimens; i++)
	{
		tmp_r[i] =0.;
	}
	if(S->size()!=0)
	{
		tmp = S->at(0).getKPIVec();
		for(unsigned i = 0; i < dimens; i++)
		{
			tmp_r[i] += tmp[i];
		}
	}
	for(unsigned i = 1 ; i < S->size(); i++)
	{
		MathFunctions::VectorAddition(tmp_r,S->at(i).getKPIVec(),dimens);
	}
	return tmp_r;
}
int emptysize;
int eemptysize;
int emptyiter;
int compareiter;
int nonemptiter;
int neemptysize;
bool packset(std::vector<Job> *Input, std::vector<Job> *Output,int dimens)
{
	std::vector<Job> *Input_ = new std::vector<Job>();
	std::vector<Job> *Q = new std::vector<Job>();
	Input_->assign(Input->begin(),Input->end());
	iter++;
	int size;
	if(Input_->empty())
	{
		emptyiter++;
		if(distance(Output,dimens)>distance(P,dimens))
		{
			++iteration;
			P->assign(Output->begin(),Output->end());

		}
	}
	else
	{
		nonemptiter++;
		float * tmp_job;
		float * tmp_T = new float[dimens];
		for(unsigned i = 0; i < dimens; i++)
		{
			tmp_T[i] = 0.;
		}
		Job job0=Input_->at(0);
		tmp_job =loadV(Output,dimens);
		for(unsigned i = 0; i < dimens; i++)
		{
			tmp_T[i] += tmp_job[i];
		}
		delete[] tmp_job;
		tmp_job = job0.getKPIVec();
		MathFunctions::VectorAddition(tmp_T,tmp_job,dimens);
		Input_->erase(Input_->begin());
		if(compare(tmp_T,1,dimens))
		{
			compareiter++;
			Q->assign(Output->begin(),Output->end());
			Q->push_back(job0);
			if(Input_->size()==0){
				emptysize++;
			}
			packset(Input_,Q,dimens);
		}
		delete[] tmp_T;
		if(Input_->size()==0){
            eemptysize++;
        }
        else{neemptysize++;}
		packset(Input_,Output,dimens);
	}
	for(int i = Input_->size()-1; i >=0; i--){
		Input_->erase(Input_->begin()+i);
	}
	for(int i = Q->size()-1; i >=0; i--){
		Q->erase(Q->begin()+i);
	}
	delete Input_;
	delete Q;
	return true;
}

extern "C"
int planJobs(int &counter, std::vector <Machine> **out_machines, std::vector <Job> *_jq, std::vector <Machine> * _mq, float _epsilon) {
	std::vector <Job*>* aJQ = new std::vector <Job*> ();

	float* tmp_jv;
	int dimensions = _jq->at(0).getDimensions();
	float* LBvector = new float[dimensions];
	for(unsigned j = 0; j < dimensions; j++){
		LBvector[j] = 0.;
	}
	for(unsigned i = 0; i < _jq->size(); i++)
	{
		tmp_jv = _jq->at(i).getKPIVec();
		for(unsigned j = 0; j < dimensions; j++){
			LBvector[j] += tmp_jv[j];
		}
		aJQ->push_back(&(_jq->at(i)));
	}
	int LB = int(MathFunctions::max(LBvector,dimensions)) +1;
	printf("lower bound is %d\n",LB);
	//define variables
	MachineGenerator mach_gen(dimensions);
	std::vector <Machine>* mq = mach_gen.generateMachines(1, 1.0);
	std::vector <Machine>* TC = mach_gen.generateMachines(1, 3.0);
	std::vector <Machine>* TCtmp = mach_gen.generateMachines(1, 3.0);
	float* tmp_mv;
	// initialized solution
	for(unsigned i = 0; i < aJQ->size(); ++i)
	{
		tmp_jv = aJQ->at(i)->getKPIVec();
		for(unsigned j = 0; j < mq->size();++j )
		{
			tmp_mv = mq->at(j).getRemainingCapacity();
			if(!MathFunctions::isVectorDifferencePositive(tmp_mv, tmp_jv, dimensions))
			{
				if(j == mq->size()-1)
				{
					mq->push_back(mach_gen.generateMachine(mq->size() + 1, 1.0));
					mq->at(j+1).assignJob(aJQ->at(i));
					break;
				}
				//   continue;
			}
			else
			{
				mq->at(j).assignJob(aJQ->at(i));
				break;
			}
		}
	}
	limit = mq->size() - LB;
	int JobNum =0;
	JobNum = 0;
	int mqsize = mq->size();
	if(mqsize<=3)
	{
		for(unsigned i = 0; i < mq->size(); i++)
		{
			for(unsigned j = 0; j < mq->at(i).assignedJobs->size(); j++)
			{
				S_permant.push_back(mq->at(i).assignedJobs->at(j));
				S->push_back(mq->at(i).assignedJobs->at(j));
			}
		}
		S_tmp->assign(S->begin(),S->end());
		for(unsigned i = 0 ; i < mqsize; i++)
		{
			mq->pop_back();
		}
		for(int k = 0; k < S_permant.size();k++)
		{
			Job tmpjob = S_permant[k];
			float * tmpjobv;
			tmpjobv = S_permant[k].getKPIVec();
		}
		while(S->size()>0)
		{
			srand(getpid());
			int a = rand();
			std::mt19937 g(a);
            std::shuffle(std::begin(*S_tmp), std::end(*S_tmp), g);
			if(S_tmp->size()>20){
				while(S_tmp->size()>20)
					S_tmp->erase(S_tmp->begin()+S_tmp->size()-1);
			}
			iter = 0;
			iteration = 0;
			packset(S_tmp,T,dimensions);
			if(P->size()!=0)
			{
				mq->push_back(mach_gen.generateMachine(mq->size() + 1, 1.0));
				for(unsigned i = 0; i< P->size();i++)
				{
					mq->back().assignJob(&P->at(i));
				}
				for(int k = 0; k < P->size();k++)
				{
					float * tmpjobv;
					tmpjobv = P->at(k).getKPIVec();
				}
			}
			tmp_mv = mq->back().getRemainingCapacity();
			T->clear();
			for(unsigned j = 0; j < P->size(); j++)
			{
				for(int i = S->size()-1; i >= 0; i--)
				{
					if(S->at(i).getID()==P->at(j).getID())
					{
						S->erase(S->begin()+i);
						break;
					}
				}
			}
			if(S->size()>0)
			{
				S_tmp->assign(S->begin(),S->end());
			}
			P->clear();
		}
	}
	if(mqsize>3)
	{
		float totalTime = 0;
#ifdef DEBUG	
	    printf("mqsize start with %ld\n", mq->size());
#endif	
//		limit += 50;
	while((mq->size() > LB) && (totalTime < 3600) && (limit > 0)){
            printf("LB is %d, current size of bins is %ld\n",LB,mq->size());
	    auto begin =  std::chrono::steady_clock::now();
			float Time = 0;
			limit--;
#ifdef DEBUG
			printf("limit is %d\n",limit);
#endif
			int inimqsize = mq->size();
			int deleteList[3];
			if(mq->size()>3){
				JobNum =0;
				int BigFlag = 0;
				int mainiter =10;
				while(mainiter)
				{
			
					totaljob = 0;
					for(unsigned i = 0; i < mq->size(); i++){
                                		totaljob+= mq->at(i).assignedJobs->size();
                        		}
#ifdef DEBUG
                			printf("limit %d mainiter %d total job is %d mqsize is %ld\n ",limit,mainiter,totaljob,mq->size());
#endif
					srand(getpid());
                    int a = rand();
                    std::mt19937 g(a);
                    std::shuffle(std::begin(*S_tmp), std::end(*S_tmp), g);
					mqsize = mq->size();
					int MQSIZETABU = mq->size();
					for(int i = (mqsize-1); i >= (mqsize-3); --i)
					{
						for(int j = 0; j < mq->at(i).assignedJobs->size(); j++)
						{
							tmp_jv =mq->at(i).assignedJobs->at(j).getKPIVec();
							TC->back().assignJob(&mq->at(i).assignedJobs->at(j));
						}
						mq->erase(mq->begin()+i);
					}
					BigFlag = 0;
					int ** Tabu;
					Tabu = resetTabu(mq);
					iter = 0;
					tmp_mv = TC->back().getRemainingCapacity();
					float optimal = distance(TC->back().assignedJobs,dimensions);
					optimal = aJQ->size()*optimal - TC->back().assignedJobs->size();
					int jobnumber = 0;
					jobnumber += TC->back().assignedJobs->size();
					for(unsigned i = 0; i < mq->size(); i++){
						jobnumber += mq->at(i).assignedJobs->size();
					}
					float ooo = optimal;
					int tabuiter = _mq->size() * _jq->size();
                    // store pair for trash can
					int S[2];
					int T[2];
					float lastobject = 0.;
					while(tabuiter)
					{
#ifdef DEBUG
printf("tabuiter is %d\n",tabuiter);
#endif
						S[0]= -1;
						S[1]= -1;
						T[0]= -1;
						T[1]= -1;
						int b = -1;
						float min = 100000;
						//store pair for bin
						int **TCpair = GeneratePair(TC->back().assignedJobs->size());
						int SizeofTCPair = SizeofPair(TC->back().assignedJobs->size());
//						printf("Get pair to swap\n");
						for(unsigned i = 0; i < mq->size(); ++i)
						{
							int **Binpair = GeneratePair(mq->at(i).assignedJobs->size());
							int SizeofBinPair = SizeofPair(mq->at(i).assignedJobs->size());
							for(unsigned j = 0; j < SizeofBinPair; ++j)
							{
#ifdef DEBUG	
												printf("pair %d in machine %d\n",j,i);
#endif							
								if(isnotTabu(Binpair[j],Tabu[i]))
								{
									for(unsigned k = 0; k < SizeofTCPair; ++k)
									{
										float tmp_b[dimensions];
										float tmp_tc[dimensions];
										float tmp_tcw[dimensions];
										tmp_mv = mq->at(i).getRemainingCapacity();
										float tmp_mw[dimensions];
										for(unsigned i = 0; i < dimensions; i++ )
										{
											tmp_mw[i] = tmp_mv[i];
										}
										tmp_mv = mq->at(i).assignedJobs->at(Binpair[j][0]).getKPIVec();
										for(unsigned i = 0; i < dimensions; i++ )
										{
											tmp_b[i] = tmp_mv[i];
										}
										if(Binpair[j][1]>0)
										{
											MathFunctions::VectorAddition(tmp_b,mq->at(i).assignedJobs->at(Binpair[j][1]).getKPIVec(),dimensions);
										}
										//machine capacity increase with weight of pair
										MathFunctions::VectorAddition(tmp_mw,tmp_b,dimensions);
										tmp_mv = TC->back().assignedJobs->at(TCpair[k][0]).getKPIVec();
										for(unsigned i = 0; i < dimensions; i++ ){
											tmp_tc[i] = tmp_mv[i];
										}
										if(TCpair[k][1]>0)
										{ MathFunctions::VectorAddition(tmp_tc,TC->back().assignedJobs->at(TCpair[k][1]).getKPIVec(),dimensions);}
										tmp_mv = TC->back().getRemainingCapacity();
										for(unsigned i = 0; i < dimensions; i++ ){
											tmp_tcw[i] = tmp_mv[i];
										}
										//TC capacity increase with weight of TCpair
										MathFunctions::VectorAddition(tmp_tcw,tmp_tc,dimensions);
										//check legal swap
										if(MathFunctions::isVectorDifferencePositive(tmp_mw, tmp_tc, dimensions)&&MathFunctions::isVectorDifferencePositive(tmp_tcw, tmp_b, dimensions))
										{
											int SizeSub = SizeofItems(Binpair[j]) - SizeofItems(TCpair[k]);
											//if machine's pair number bigger than TC's
											if(SizeSub >= 0)
											{
												float WeightSub=5;
												Totaliter++;
												Machine* Mp = &mq->at(i);
												Machine* TCp = &TC->back();
												//L2 norm of Bin pair and TC pair
												float WeightofM = WeightofPair(Mp,Binpair[j],dimensions);
												float WeightofTC = WeightofPair(TCp,TCpair[k],dimensions);
												//check whether weight of Pair equal to values got from above.
												WeightSub =WeightofM - WeightofTC;
												float current = aJQ->size()*WeightSub - SizeSub;
												if(current <= min)
												{
													//changed
													S[1]=-1;
													T[1]=-1;
													min = current;
													b = i;
													S[0] = Binpair[j][0];
													if(SizeofItems(Binpair[j])==2)
														S[1] = Binpair[j][1];
													T[0] = TCpair[k][0];
													if(SizeofItems(TCpair[k])==2)
														T[1] = TCpair[k][1];
												}
											}
										}
									}
								}
							}
							for(int i = 0; i < SizeofBinPair; i++)	
								delete[] Binpair[i];
							delete[] Binpair;
						}
						if(b == -1){ 
#ifdef DEBUG
						printf("No avalibale swap\n");
#endif						
						break;
						}
						else
						{
#ifdef DEBUG
							printf("Try to swap\n");
#endif						
							tmp_mv = TC->back().getLoad();
							float tmp_weight[dimensions];
							for(unsigned i = 0; i < dimensions; i++ ){
								tmp_weight[i] = tmp_mv[i];
							}
							MathFunctions::VectorSubtraction(tmp_weight,TC->back().assignedJobs->at(T[0]).getKPIVec(),dimensions);
							//second job
							if(T[1]>0)
							{
								MathFunctions::VectorSubtraction(tmp_weight,TC->back().assignedJobs->at(T[1]).getKPIVec(),dimensions);
							}
							//deal with pair in machine
							MathFunctions::VectorAddition(tmp_weight,mq->at(b).assignedJobs->at(S[0]).getKPIVec(),dimensions);
							if(S[1]>0)
							{
								MathFunctions::VectorAddition(tmp_weight,mq->at(b).assignedJobs->at(S[1]).getKPIVec(),dimensions);
							}
							// get distance after swap
							int SizeSub = getArrayLen(S)- getArrayLen(T);
							float object = 0;
							object = aJQ->size()*Max(tmp_weight,dimensions) - SizeSub - TC->back().assignedJobs->size();
							if(object!=lastobject)
							{
								lastobject = object;
							}
							else{
#ifdef DEBUG
								printf("failed to swap, break\n");
#endif
								break;
							}
							//swap
							if(object < optimal)
							{
								Job tmp_fj;
								Job tmp_sj;
								Job tmp_ftc;
								Job tmp_stc;
								Job * tmp_fj_tmp;
								Job * tmp_sj_tmp;
								Job * tmp_ftc_tmp;
                                Job * tmp_stc_tmp;
								optimal = object;
								tmp_fj=mq->at(b).assignedJobs->at(S[0]);
								tmp_fj_tmp=&tmp_fj;
								int s1flag = 0;
								if(S[1]>0)
								{
									s1flag =1;
									tmp_sj =mq->at(b).assignedJobs->at(S[1]);
									tmp_sj_tmp=&tmp_sj;
									mq->at(b).removeJob(S[1]);
								}
								mq->at(b).removeJob(S[0]);

								tmp_ftc=TC->back().assignedJobs->at(T[0]);
								tmp_ftc_tmp=&tmp_ftc;
								int t1flag=0;
								if(T[1]>0)
								{
									t1flag=1;
									tmp_stc=TC->back().assignedJobs->at(T[1]);
									TC->back().removeJob(T[1]);
									tmp_stc_tmp=&tmp_stc;
								}
								TC->back().removeJob(T[0]);
								mq->at(b).assignJob(tmp_ftc_tmp);
								if(T[1]>0)
								{
									mq->at(b).assignJob(tmp_stc_tmp);
								}
								TC->back().assignJob(tmp_fj_tmp);
								if(S[1]>0)
								{
									TC->back().assignJob(tmp_sj_tmp);
								}
								Swapiter++;
								updateTabu(S,T,Tabu[b],mq->at(b));
								tmp_mv = TC->back().getRemainingCapacity();
								totaljob = 0;	
								for(unsigned i = 0; i < mq->size(); i++){
                                					totaljob+= mq->at(i).assignedJobs->size();
                        					}
								totaljob += TC->back().assignedJobs->size();
#ifdef DEBUG
                			printf("222222 limit %d mainiter %d total job is %d mqsize is %ld\n ",limit,mainiter,totaljob,mq->size());
#endif
							}
                            for(unsigned i = 0; i < SizeofTCPair; i++)
                                delete[] TCpair[i];
                            delete[] TCpair;
						}
#ifdef DEBUG	
					printf("Finish Tabu\n");
#endif
						if(Min(TC->back().getRemainingCapacity(),dimensions) > 1.0)
						{
#ifdef DEBUG
							printf("Try to pack TC into two\n");
#endif
							int a = mq->size();
							tmp_mv=TC->back().getRemainingCapacity();
							cleanJobVector(TCP_tmp);
							cleanJobVector(TT);
							for(unsigned j = 0; j < TC->back().assignedJobs->size(); j++)
							{
								TCP_tmp->push_back(TC->back().assignedJobs->at(j));
							}
							for(int j = 0; j < TCP_tmp->size(); j++)
							{
								tmp_jv =TCP_tmp->at(j).getKPIVec();
							}
							TCP->assign(TCP_tmp->begin(),TCP_tmp->end());
							while(TCP->size()>0)
							{	
				//				printf("Start packset in TrytopackTC\n");
								iteration = 0;
								std::shuffle(std::begin(*TCP_tmp), std::end(*TCP_tmp), g);
								if(TCP_tmp->size()>20){
									while(TCP_tmp->size() > 20){
										TCP_tmp->erase(TCP_tmp->begin()+TCP_tmp->size()-1);
									}
								}
								packset(TCP_tmp,TT,dimensions);
								cleanJobVector(TCP_tmp);
								if(P->size()!=0)
								{
									mq->push_back(mach_gen.generateMachine(mq->size() + 1, 1.0));
									for(unsigned i = 0; i< P->size();i++)
									{
										mq->back().assignJob(&P->at(i));
									}
								}
								tmp_mv=mq->back().getRemainingCapacity();
								cleanJobVector(TT);
								for(unsigned j = 0; j < P->size(); j++)
								{
									for(int i = TCP->size()-1; i >= 0; i--)
									{
										if(TCP->at(i).getID()==P->at(j).getID())
										{
											TCP->erase(TCP->begin()+i);
											break;
										}
									}
								}
								if(TCP->size()>0)
								{
									TCP_tmp->assign(TCP->begin(),TCP->end());
								}
								
								cleanJobVector(P);
							}
							if((mq->size()-a)>2){
								totaljob = 0;	
								for(unsigned i = 0; i < mq->size(); i++){
                                    totaljob+= mq->at(i).assignedJobs->size();
                                }
#ifdef DEBUG
                			printf("333333 limit %d mainiter %d total job is %d mqsize is %ld\n ",limit,mainiter,totaljob,mq->size());
#endif 
								while(mq->size()>a){
									mq->erase(mq->begin()+(mq->size()-1));
								}
								--tabuiter;
							}else{
								TC->erase(TC->begin());
								TC->push_back(mach_gen.generateMachine(TC->size() + 1, 3.0));
								BigFlag = 1;
								totaljob = 0;	
								for(unsigned i = 0; i < mq->size(); i++){
                                					totaljob+= mq->at(i).assignedJobs->size();
                        					}
#ifdef DEBUG
                			printf("444444 limit %d mainiter %d total job is %d mqsize is %ld\n ",limit,mainiter,totaljob,mq->size());
#endif 
								break;
							}
						}
						else
							--tabuiter;
					}
					if(BigFlag == 1){break;}
						//start decent
#ifdef DEBUG
					printf("start decent\n");
#endif
					
					totaljob = 0;	
					for(unsigned i = 0; i < mq->size(); i++){
                           				totaljob+= mq->at(i).assignedJobs->size();
                        		}
					totaljob += TC->back().assignedJobs->size();
#ifdef DEBUG
                			printf("555555 limit %d mainiter %d total job is %d mqsize is %ld\n ",limit,mainiter,totaljob,mq->size());
#endif 
					JobNum =0;
					for(unsigned i = 0; i < sizeof(Tabu)/sizeof(Tabu[0]); i++)
						delete[] Tabu[i];
					delete[] Tabu;
					int mqsize = mq->size();
					int PSize= 0;
					int totalJobSize = 0;
                    cleanJobVector(TCP_tmp);
                    cleanJobVector(TT);
					std::shuffle(std::begin(*mq), std::end(*mq), g);
					while(mqsize)
					{
						int machineJobSize = mq->at(0).assignedJobs->size();
						int remainMqSize = 0;
						remainMqSize = _jq->size() - machineJobSize - TC->back().assignedJobs->size();
						int TCJobSize = TC->back().assignedJobs->size();
						totalJobSize = machineJobSize + TCJobSize;
						for(int i = mq_tmp->size()-1; i >=0;i--){
							mq_tmp->erase(mq_tmp->begin()+i);
						}
						for(int i = TC_tmp->size()-1; i >=0;i--){
							TC_tmp->erase(TC_tmp->begin()+i);
						}
						for(unsigned j = 0; j < mq->at(0).assignedJobs->size(); j++)
						{
							mq_tmp->push_back(mq->at(0).assignedJobs->at(j));
							TCP_tmp->push_back(mq->at(0).assignedJobs->at(j));
						}
						for(unsigned j = 0; j < TC->back().assignedJobs->size(); j++)
						{
							TC_tmp->push_back(TC->back().assignedJobs->at(j));
							TCP_tmp->push_back(TC->back().assignedJobs->at(j));
						}
						TCP->assign(TCP_tmp->begin(),TCP_tmp->end());
						iter = 0;
						iteration = 0;
						eemptysize = 0;
						emptysize = 0;
						emptyiter = 0;
						nonemptiter = 0;
						neemptysize = 0;
						compareiter = 0;
						std::shuffle(std::begin(*TCP_tmp), std::end(*TCP_tmp), g);
                        if(TCP_tmp->size()>20){
                            while(TCP_tmp->size() > 20){
                                    TCP_tmp->erase(TCP_tmp->begin()+TCP_tmp->size()-1);
                            }
                        }
						PSize = P->size();
						if(P->size()!=0){ printf("BBBBBBUUUUUGGGG\n");
							return 0;}
						packset(TCP_tmp,TT,dimensions);
#ifdef DEBUG
					//	printf("limit is %d tabuiter is %d iter iis %d neemptysize is %d eemptysize is %d emptysize is %d emptyiter is %d nonemptiter is %d compare iter is %d \n",limit,tabuiter,iter,neemptysize,eemptysize,emptysize,emptyiter,nonemptiter,compareiter);
#endif
						cleanJobVector(TCP_tmp);
						mq->erase(mq->begin());
						TC->erase(TC->begin());
						PSize = P->size();
						//Assign jobs into machine
						if(P->size()!=0)
						{
							mq->push_back(mach_gen.generateMachine(mq->size() + 1, 1.0));
							for(unsigned i = 0; i< P->size();i++)
							{
								mq->back().assignJob(&P->at(i));												  }
						}
						if(PSize != mq->back().assignedJobs->size()){
							printf("Here it is \n\n\n\n");
							return 0;
						}
						for(unsigned g = 0; g < P->size(); g++)
						{
							for(int k = TCP->size()-1; k >= 0; k--)
							{
								if(TCP->at(k).getID()==P->at(g).getID())
								{
									TCP->erase(TCP->begin()+k);
									break;
								}
							}
						}
						if((TCP->size() + PSize) != totalJobSize){
							printf("It IS Here \n\n\n");
							return 0;
							}
						if(TCP->size()!=0)
						{
							TC->push_back(mach_gen.generateMachine(TC->size() + 1, 3.0));
							for(unsigned h = 0; h< TCP->size();h++)
							{
								TC->back().assignJob(&TCP->at(h));
							}
						}
						else
						{
							TC->push_back(mach_gen.generateMachine(TC->size() + 1, 3.0));
						}
						if(TCP->size()!=TC->back().assignedJobs->size()){
							//printf("ATENTION ATENTION ATENTION JOBS CAN NOT BE PACKED IN TO TC, START REVERSE\n");
							mq->erase(mq->begin()+(mq->size()-1));
							mq->push_back(mach_gen.generateMachine(mq->size() + 1, 1.0));
							TC->erase(TC->begin()+(TC->size()-1));
							TC->push_back(mach_gen.generateMachine(TC->size() + 1, 3.0));
							for(unsigned i = 0; i < mq_tmp->size(); i++){
								mq->back().assignJob(&mq_tmp->at(i));
							}	
							for(unsigned i = 0; i < TC_tmp->size(); i++){
								TC->back().assignJob(&TC_tmp->at(i));
							}
						}
                        cleanJobVector(TT);
                        cleanJobVector(P);
						mqsize--;
						int remainMqSize_copy = 0;
						remainMqSize_copy = _jq->size() - mq->back().assignedJobs->size() - TC->back().assignedJobs->size();
						if(remainMqSize_copy != remainMqSize){
							printf("remainMq size is not equal original size is %d, new size is %d TCP size is %ld, TC size is %ld mqjob size is %ld P size is %d\n",remainMqSize,remainMqSize_copy,TCP->size(),TC->back().assignedJobs->size(), mq->back().assignedJobs->size(), PSize);
							return 0;
						}
						totaljob = 0;	
						for(unsigned i = 0; i < mq->size(); i++){
                            totaljob+= mq->at(i).assignedJobs->size();
                        }
						totaljob += TC->back().assignedJobs->size();
						if(totaljob != _jq->size()){
							printf("totajob is %d P size is %d TC size is %ld totalTrash size is %d, last machine size is %ld \n",totaljob,PSize,TC->back().assignedJobs->size(),totalJobSize,mq->back().assignedJobs->size());
							return 0;
						}
					}
#ifdef DEBUG
                			printf("666666 limit %d mainiter %d total job is %d mqsize is %ld\n ",limit,mainiter,totaljob,mq->size());
#endif
					cleanJobVector(TCP_tmp);
					for(unsigned j = 0; j < TC->back().assignedJobs->size(); j++)
					{
						TCP_tmp->push_back(TC->back().assignedJobs->at(j));
					}
					for(int j = 0; j < TCP_tmp->size(); j++)
					{
						tmp_jv =TCP_tmp->at(j).getKPIVec();
					}
					TCP->assign(TCP_tmp->begin(),TCP_tmp->end());
					while(TCP->size()>0){
						std::shuffle(std::begin(*TCP_tmp), std::end(*TCP_tmp), g);
                	    if(TCP_tmp->size()>20){
                            while(TCP_tmp->size() > 20){
                                TCP_tmp->erase(TCP_tmp->begin()+TCP_tmp->size()-1);
                            }
                        }
						packset(TCP_tmp,TT,dimensions);
						cleanJobVector(TCP_tmp);
						if(P->size()!=0)
						{
							mq->push_back(mach_gen.generateMachine(mq->size() + 1, 1.0));
							for(unsigned i = 0; i< P->size();i++)
							{
								mq->back().assignJob(&P->at(i));
							}
						}
						tmp_mv=mq->back().getRemainingCapacity();							//Q->clear();
						cleanJobVector(TT);
						for(unsigned j = 0; j < P->size(); j++)
						{
							for(int i = TCP->size()-1; i >= 0; i--)
							{
								if(TCP->at(i).getID()==P->at(j).getID())
								{
									TCP->erase(TCP->begin()+i);
									break;
								}
							}
						}
						if(TCP->size()>0)
						{
							TCP_tmp->assign(TCP->begin(),TCP->end());
						}
						cleanJobVector(P);
					}

					//	printf("After packed TC Machine size is %d\n",mq->size());
					totaljob = 0;	
					for(unsigned i = 0; i < mq->size(); i++){
                           				totaljob+= mq->at(i).assignedJobs->size();
                        		}
#ifdef DEBUG
                			printf("777777 limit %d mainiter %d total job is %d mqsize is %ld\n ",limit,mainiter,totaljob,mq->size());
#endif 
					TC->erase(TC->begin());
					TC->push_back(mach_gen.generateMachine(TC->size() + 1, 3.0));
					--mainiter;
				}
			}
            auto end =  std::chrono::steady_clock::now();
            Time = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
            totalTime += std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
            printf("Time is %f Total elapsed time is %f\n",Time,totalTime);
            if(totalTime > 3600){
                printf("33333333Total elapsed time is %f\n",totalTime);
                break;
			}
		}
	}
	totaljob = 0;
	for(unsigned i = 0; i < mq->size(); i++){
		totaljob+= mq->at(i).assignedJobs->size();
	}
	mqsize = mq->size();
	delete aJQ;
	if(mq->size() > _mq->size()) {
		counter++;
	}
	*out_machines = mq;
	return 1;
}



