//=================================
// include guard
#ifndef __HistogramAvg_H_INCLUDED__
#define __HistogramAvg_H_INCLUDED__

//=================================
// forward declared dependencies
class RunningStats;
//=================================
// included dependencies


//=================================
// the actual class
class HistogramAvg{
	
	public:
	

	void Initialize(unsigned long long Nbins, long double Taile, long double Tope){
		nBins=Nbins;
		
		sBin= new RunningStats[nBins];
		Tail = Taile;
		Top = Tope;
		incr = (Top-Tail)/((long double)(nBins));
		return;
	}
	void Push(long double val, long double at){
		
		if(at>=Tail && at < Top){
			unsigned long long Ind = (unsigned long long)((at - Tail)/incr);
		//cout << "value " << val << " is being pushed into bin number " << Ind << endl;
			
			sBin[Ind].Push(val);
			
		}
		return;
		
	}

	
	void Dump(string filename){
		ofstream outfile( filename.c_str());
		
		for(unsigned long long i=0;i<nBins;i++){
			long double pop = Tail+((long double)(i+1)*incr); 
			outfile << pop << " " << sBin[i].Mean() << " " << sBin[i].StandDev() << endl;
		}
		outfile.close();
		return;
	}

	
	long double GetIncrement(void){
		return(incr);
	}
	
	void Reset(void){
		for(unsigned long long i=0;i<nBins;++i){
			
			sBin[i].Reset();
		}
		return;
	}
	
	

	

	private:
	unsigned long long nBins;
	RunningStats* sBin;
	long double Tail, Top;
	long double incr;
};


#endif // __MYCLASS_H_INCLUDED__ 
