//=================================
// include guard
#ifndef __BlockAverage_H_INCLUDED__
#define __BlockAverage_H_INCLUDED__

//=================================
// forward declared dependencies
class RunningStats;
//=================================
// included dependencies


//=================================
// the actual class
class BlockAverage{
	
	public:
	

	void Initialize(unsigned nblock, unsigned nperblock){
		nBlocks=nblock;
		nperBlock=nperblock;
		block= new RunningStats[nBlocks];
		cblock=0;
		npushed=0;
		return;
	}
	void Push(long double val){
		//push the value to running stat block
		if (cblock<nBlocks){
			block[cblock].Push(val);
		
			//Increment the counter
			++npushed;
			//check the count
			if (npushed==nperBlock){
				//increment to the next block
				++cblock;
				npushed=0;
			}
		}
		return;
		
	}

	
	

	//Get the average of the block averages
	long double Average(void){
		long double sum=0.0;
		//only loop over filled blocks
		for (unsigned int i = 0; i < cblock; ++i)
		{
			sum+=block[i].Mean();
		}
		long double N = (long double)(cblock);
		long double avg = sum/N;
		return avg;
	}
	
	long double StandardError(void){
		RunningStats Temp;
		//only loop over filled blocks
		for (unsigned int i = 0; i < cblock; ++i)
		{
			Temp.Push(block[i].Mean());
		}
		long double sigma = Temp.StandDev();
		long double serr = 2.0*sigma/sqrt(cblock);
		return serr;

	}
	void Reset(void){
		for(unsigned long long i=0;i<nBlocks;++i){
			
			block[i].Reset();
			cblock=0;
			npushed=0;
		}
		return;
	}
	
	

	

	private:
	unsigned npushed;
	unsigned cblock;
	unsigned nBlocks;
	unsigned nperBlock;
	RunningStats* block;
	
};


#endif // __MYCLASS_H_INCLUDED__ 
