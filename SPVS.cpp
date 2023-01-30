#include "Job.h"
#include "Machine.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <map>


bool isOneFree(int n,bool*P)
{
    for(int i=0; i<n; i++)
    {
        if(*(P+i)==true)
            return true;
    }
    return false;
}

int* stableMarriage(int n,int**MP,int**WP)
{
    bool isMachineFree[n];
    bool isJobFree[n];
    bool **isMachineProposed = new bool*[n];
    for(int j=0;j<n;j++)
        isMachineProposed[j]=new bool[n];
    int *match = new int[n];
    for(int i=0; i<n; i++)
    {
        isMachineFree[i]=true;
        isJobFree[i]=true;
        for(int j=0; j<n; j++)
            isMachineProposed[i][j]=false;
        match[i]=-1;
    }
    while(isOneFree(n,isMachineFree))
    {
        int indexM;
        for(int i=0;i<n;i++)
        {
            if(isMachineFree[i]==true)
            {
                indexM=i;
                break;
            }
        }
        int indexJ;
        for(int i=0; i<n; i++)
        {
            int w = MP[indexM][i];
            if(isMachineProposed[indexM][w] == false)
            {
                indexJ = w;
                break;
            }
        }
        isMachineProposed[indexM][indexJ]=true;
        if(isJobFree[indexJ])
        {
            isMachineFree[indexM]=false;
            isJobFree[indexJ]=false;
            match[indexM]=indexJ;
        }
        else
        {
            int indexRival;
            for(int i=0; i<n; i++)
            {
                if(match[i]==indexJ)
                {
                    indexRival = i;
                    break;
                }
            }
            int pM,pRival;
            for(int i=0; i<n; i++)
            {
                if(WP[indexJ][i]==indexM)
                    pM=i;
                if(WP[indexJ][i]==indexRival)
                    pRival=i;
            }
            if(pM<pRival)
            {
                isMachineFree[indexM]=false;
                isMachineFree[indexRival]=true;
                isJobFree[indexJ]=false;
                match[indexM]=indexJ;
            }
        }
    }
    for(int j=0;j<n;j++)
        delete [] isMachineProposed[j];
    delete [] isMachineProposed;
    return match;
}

class SRMJob;
class SRMJobPrio {
public:
    SRMJob* job;
    float prio;
    
    SRMJobPrio() {
        job = 0;
        prio = -1;
    }
    SRMJobPrio(SRMJob* _j, float _p) {
        job = _j;
        prio = _p;
    }
    SRMJobPrio(const SRMJobPrio& _j) {
        job = _j.job;
        prio = _j.prio;
    }
    SRMJobPrio& operator=(const SRMJobPrio& _j) {
        job = _j.job;
        prio = _j.prio;
        return *this;
    }
    bool operator<(const SRMJobPrio& _j) const {
        return prio < _j.prio;
    }
};

class SRMJob {
public:
    static int KPIS;
    static int ID;
    int id;
    std::vector <SRMJobPrio> prefs;
    float* kpis;
    int favorite;
    Job* job;
    
    SRMJob() {
        id = ID++;
        kpis = new float[KPIS];
        favorite = -1;
        job = 0;
        for(int i = 0; i < KPIS; i++)
            kpis[i] = drand48();
        prefs = std::vector <SRMJobPrio> ();
    }
    
    SRMJob(Job* _j) {
        id = ID++;
        kpis = new float[KPIS];
        favorite = -1;
        job = _j;
        for(int i = 0; i < KPIS; i++)
            kpis[i] = job->getKPIVec()[i];
        prefs = std::vector <SRMJobPrio> ();
    }
    
    SRMJob(const SRMJob& _j) {
        id = _j.id;
        favorite = _j.favorite;
        job = _j.job;
        kpis = new float[KPIS];
        memcpy(kpis, _j.kpis, KPIS*sizeof(float));
        prefs = std::vector <SRMJobPrio> (_j.prefs);
    }
    
    ~SRMJob() {
        delete [] kpis;
    }
    
    SRMJob& operator=(const SRMJob& _j) {
        id = _j.id;
        favorite = _j.favorite;
        memcpy(kpis, _j.kpis, KPIS*sizeof(float));
        prefs = std::vector <SRMJobPrio> (_j.prefs);
        return *this;
    }
    
    bool addJob(SRMJob* _j) {
        int flag = 1;
        for(unsigned i = 0; i < KPIS; i++)
        {
            if(kpis[i] + _j->kpis[i] > 1)
            {
                flag = 0;
                printf("jobs is so big\n");
                return false;
                break;
            }
        }
        if(flag)
        {
            float locScalar = 0;
            //get the distance of origin vector
            for(int i = 0; i < KPIS; i++)
                locScalar += kpis[i] * kpis[i];
            locScalar = std::sqrt(locScalar);
            
            float extScalar = 0;
            //get the distance of new vector
            for(int i = 0; i < KPIS; i++)
                extScalar += _j->kpis[i] * _j->kpis[i];
            extScalar = std::sqrt(extScalar);
            
            float tmpPrio = 0;
            for(int i = 0; i < KPIS; i++)
                tmpPrio += (kpis[i]*_j->kpis[i]);
            tmpPrio += tmpPrio/(locScalar*extScalar);
            //      printf("TmpPrio is %f\n",tmpPrio);
            //check
            //     std::cout<<i<<"change is"<<change<<std::endl;
            prefs.push_back(SRMJobPrio(_j, fabs(tmpPrio)));
            return true;
        }
        return false; // LN: my compiler complained, that is why I added this line
    }
    
    void sortPrios() {
        std::sort(prefs.begin(), prefs.end());
    }
    
    int makeProposal(){
        SRMJob* fav = prefs[0].job;
        if(fav->favorite == -1) {
            fav->favorite = id;
            return -1;
        } else if(fav->favorite == id) {
            return -1;
        }else {
            //here
            return fav->checkForBetterRoommate(id);
        }
    }
    
    void removeCandidate(int _id) {
        for(std::vector<SRMJobPrio>::iterator it = prefs.begin(); it != prefs.end(); it++) {
            if((*it).job->id == _id) {
                prefs.erase(it);
                return;
            }
        }
    }
    
    void removeFavorite() {
        prefs.erase(prefs.begin());
    }
    
    int checkForBetterRoommate(int _id) {
        int oldFavorite = favorite;
        int currID = 0;
        favorite = -1;
        int delFirst = -1;
        for(int i = 0; i < prefs.size(); i++) {
            currID = prefs[i].job->id;
            if(favorite == -1) {
                if(currID == _id) {
                    favorite = _id;
                    removeCandidate(oldFavorite);
                    return oldFavorite;
                } else if (currID == oldFavorite) {
                    favorite = oldFavorite;
                    removeCandidate(_id);
                    return _id;
                }
            }
        }
        
        return delFirst;
    }
    
    void reducePreflists(int _minCandidate) {
        std::vector <SRMJobPrio>::iterator it;
        for(it = prefs.begin(); it != prefs.end(); it++) {
            if((*it).job->id == _minCandidate)
                break;
        }
        if(it == prefs.end()) {
            printf("OHOH - not found\n");
            return;
        }
        it++;
        if(it != prefs.end()) {
            std::vector <SRMJobPrio>::iterator it2 = it;
            for(; it2 != prefs.end(); it2++) {
                (*it2).job->removeCandidate(id);
            }
        }
        prefs.erase(it, prefs.end());
    }
    
    void print() {
        if(prefs.empty())
            return;
        
        printf("id: %03d{{%03d, %.6f}", id, prefs[0].job->id, prefs[0].prio);
        for(int i = 1; i < prefs.size(); i++)
            printf(", {%03d, %.6f}", prefs[i].job->id, prefs[i].prio);
        printf("}\n");
    }
};

int SRMJob::KPIS = 3;
int SRMJob::ID = 0;

bool findJobInVector(SRMJob* _j, std::vector<SRMJob*>* _vec) {
    if(!_vec || !_j)
        return false;
    for(int i = 0; i < _vec->size(); i++) {
        if(_vec->at(i) == _j)
            return true;
    }
    return false;
}


bool buildRoommates(std::vector <SRMJob*> &_jobs) {
    for(int i = 0; i < _jobs.size(); i++) {
        for(int j = 0; j < _jobs.size(); j++) {
            if(i == j)
                continue;
            if(!_jobs[i]->addJob(_jobs[j]))
                return false;
        }
        _jobs[i]->sortPrios();
    }
    
    // step one: make proposals
    for (int i = 0; i < _jobs.size(); ) {
        //    printf("Start make proposal\n");
        int res = _jobs[i]->makeProposal();
        if(res >= 0) {
            //       printf("Remove favorite\n");
            _jobs[res]->removeFavorite();
            i = res;
        } else {
            i++;
        }
    }
    
    //step two: remove candidates worse than current proposal
    for (int i = 0; i < _jobs.size(); i++) {
        SRMJob* jobsFavorite = _jobs[_jobs[i]->favorite];
        jobsFavorite->reducePreflists(i);
    }
    //step three:
    bool thereAreLoops = true;
    while(thereAreLoops) {
        bool phaseThreeRequired = false;
        int firstUnfinished = 0;
        for (int i = 0; i < _jobs.size(); i++) {
            if(_jobs[i]->prefs.size() > 1) {
                phaseThreeRequired = true;
                firstUnfinished = i;
                break;
            }
        }
        if(phaseThreeRequired) {
            int counter = 0;
            std::vector <SRMJob*> rowR = std::vector <SRMJob*>();
            std::vector <SRMJob*> rowQ = std::vector <SRMJob*>();
            SRMJob* candidate = _jobs[firstUnfinished];
            SRMJob* val = 0;
            while(counter < _jobs.size()) {
                if(findJobInVector(candidate, &rowR))
                    break;
                val = candidate->prefs[candidate->prefs.size()-1].job;
                rowR.push_back(candidate);
                rowQ.push_back(val);
                candidate = val->prefs[val->prefs.size()-1].job;
                counter++;
            }
        } else {
            thereAreLoops = false;
        }
    }
    return true;
}
void buildpreference(float* man, std::vector<float*>&woman, int prefs[],int dimens)
{
  //  int dimens = 3;
    std::map<float,int>priorimap;
    std::vector<float> priorivector;
    for(unsigned i = 0; i < woman.size(); i++)
    {
        float locScalar = 0;
        //get the distance of origin vector
        for(int j = 0; j < dimens; j++)
            locScalar += man[j] * man[j];
        locScalar = std::sqrt(locScalar);
        float extScalar = 0;
        //get the distance of new vector
        for(int j = 0; j < dimens; j++)
            extScalar += woman[i][j] * woman[i][j];
        extScalar = std::sqrt(extScalar);
        float tmpPrio = 0;
        for(int j = 0; j < dimens; j++)
            tmpPrio += man[j]*woman[i][j];
        tmpPrio = tmpPrio/(locScalar*extScalar);
        float change = float(i*3)/float(10000);
        //     std::cout<<i<<"change is"<<change<<std::endl;
        tmpPrio += change;
        priorimap.insert(std::pair<float,int>(tmpPrio,i));
        priorivector.push_back(tmpPrio);
    }
    sort(priorivector.begin(), priorivector.end(),std::greater<float>());
    for(unsigned i = 0; i < woman.size(); i++)
    {
        prefs[i] = priorimap[priorivector[i]];
    }
}
struct CompareFunc {
    
    bool operator()(Job const &a, Job const &b) {
        return a.getSumOfComponents() > b.getSumOfComponents();
    }
};

extern "C"
int planJobs(int &counter, std::vector <Machine> **out_machines, std::vector <Job> *jq, std::vector <Machine> * mq, float _epsilon) {
    
    int machinesize =mq->size();
    int jqsize = jq->size();
    int dimens =jq->at(0).getDimensions();
    std::sort(jq->begin(), jq->end(), CompareFunc());
    SRMJob::ID = 0;
    //Build Roommate stage
    if(jq->size() >= 2*mq->size()) {
     //   printf("Start SRM\n");
        SRMJob::KPIS = jq->at(0).getDimensions();
        std::vector <SRMJob*> srmjobs = std::vector <SRMJob*>();
        for(unsigned i = 0; i < 2*mq->size(); i++) {
            srmjobs.push_back(new SRMJob(&jq->at(i)));
        }
        if(buildRoommates(srmjobs))
        {
            std::map<int,int>roommate;
            for(unsigned i = 0; i < srmjobs.size(); i++)
            {
                roommate.insert(std::pair<int,int>(i,srmjobs[i]->prefs[0].job->id));
            }
            int j = 0;
            for(unsigned i = 0; i < 2*mq->size(); i++)
            {
                if(roommate[i] != 0)
                {
                    roommate.erase(roommate[i]);
                    mq->at(j).assignJob((srmjobs[i]->job));
                    mq->at(j).assignJob((srmjobs[roommate[i]]->job));
                    j++;
                }
            }
            for(int i = 2*mq->size() - 1; i >= 0; i--)
            {
		jq->at(i).deleteData();
                jq->erase(jq->begin()+i);
            }

        }
        else
        {
            for(unsigned i = 0; i < 2*mq->size(); ++i) {
                for(unsigned j = 0; j < mq->size(); ++j) {
                    if(mq->at(j).assignJob(&(jq->at(i)))) {
                        break;
                    } else {
                        if(j == mq->size() - 1) {
                            std::vector<Machine>* mq_clone = new std::vector<Machine>(*mq);
			    std::vector<Machine>* tmp = *out_machines;
                            *out_machines = mq_clone;
			    // REPAIRED
			    if(tmp != NULL)
				    delete tmp;
                            ++counter;
	    		    for(unsigned i = 0; i < srmjobs.size(); i++)
		    	        delete srmjobs[i];
                            return 0;
                        }
                    }
                }
            }
            for(int i = 2*mq->size() - 1; i >= 0; i--)
            {
                jq->erase(jq->begin()+i);
            }
        }
	for(unsigned i = 0; i < srmjobs.size(); i++)
		delete srmjobs[i];
   //      printf("Finish SRM\n");
        // stable marriage
        while(jq->size() >= mq->size())
        {
            // man is machine, woman is job
            std::vector<float*> machineset;
            std::vector<float*> jobset;
            for(unsigned i = 0; i < mq->size(); i++)
            {
                float* tmp_mv = mq->at(i).getLoad();
                float* tmp_jv = jq->at(i).getKPIVec();
                machineset.push_back(tmp_mv);
                jobset.push_back(tmp_jv);
            }
            int mList[mq->size()][mq->size()];
            int wList[mq->size()][mq->size()];
            for(unsigned i = 0; i < mq->size(); i++ )
            {
                buildpreference(machineset[i], jobset, mList[i], dimens);
                buildpreference(jobset[i], machineset, wList[i], dimens);
            }
            int **man = new int*[machinesize];
            int **woman = new int*[machinesize];
            for(unsigned i = 0; i < machinesize; i++)
            {
                man[i] = new int[machinesize];
                woman[i] = new int[machinesize];
            }
            for(unsigned i = 0; i < mq->size(); i++)
            {
                for(unsigned j = 0; j < mq->size(); j++)
                {
                    man[i][j] = mList[i][j];
                    woman[i][j] = wList[i][j];
                }
            }
            //printf("Start Marriage\n");
            int *assignment = stableMarriage(machinesize,man,woman);
            //printf("Good\n");
            std::vector<int> updateJob;
            for(unsigned i = 0; i < machinesize; i++)
            {
                if( mq->at(assignment[i]).assignJob(&(jq->at(i))))
                {
            //        printf("canbe assigned \n");
                    updateJob.push_back(i);
                }
                else{
            //        printf("job is big \n");
                }
            }
            if(updateJob.size()!= 0)
            {
                for(int i = updateJob.size()-1; i >= 0; i--)
                {
                    jq->erase(jq->begin()+updateJob[i]);
                }
            }
            else
            {
                for(unsigned i = 0; i < mq->size(); ++i) {
                    for(unsigned j = 0; j < mq->size(); ++j) {
                        if(mq->at(j).assignJob(&(jq->at(i)))) {
                            break;
                        } else {
                            if(j == mq->size() - 1) {
                                std::vector<Machine>* mq_clone = new std::vector<Machine>(*mq);
                                std::vector<Machine>* tmp = *out_machines;
                                *out_machines = mq_clone;
				// REPAIR
				if(tmp != NULL)
					delete tmp;
                                ++counter;
				
				for(unsigned i = 0; i < machinesize; i++)
            			{
                			delete [] man[i];
                			delete [] woman[i];
            			}
				delete [] man;
				delete [] woman;
	    			delete [] assignment;
                                return 0;
                            }
                        }
                    }
                }
                for(int i = mq->size() - 1; i >= 0; i--)
                {
                    jq->erase(jq->begin()+i);
                }
		// REPAIR
		for(unsigned i = 0; i < machinesize; i++)
            	{
                	delete [] man[i];
                	delete [] woman[i];
            	}
		delete [] man;
		delete [] woman;
            }
	    delete [] assignment;
        }
      //    printf("Finish Marriage\n");
    }
    if(jq->size()!=0)
    {
        for(unsigned i = 0; i < jq->size(); ++i) {
            for(unsigned j = 0; j < mq->size(); ++j) {
                if(mq->at(j).assignJob(&(jq->at(i)))) {
                    break;
                } else {
                    if(j == mq->size() - 1) {
                        std::vector<Machine>* mq_clone = new std::vector<Machine>(*mq);
                        std::vector<Machine>* tmp = *out_machines;
                       	*out_machines = mq_clone;
			// REPAIR
			if(tmp != NULL)
				delete tmp;
                        ++counter;
                        return 0;
                    }
                }
            }
        }
    }
    std::vector<Machine>* mq_clone = new std::vector<Machine>(*mq);
    std::vector<Machine>* tmp = *out_machines;
    *out_machines = mq_clone;
    // REPAIR
    if(tmp != NULL)
         delete tmp;
    return 1;
}


