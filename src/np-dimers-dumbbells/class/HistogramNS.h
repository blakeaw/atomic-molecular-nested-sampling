//=================================
// include guard
#ifndef __HistogramNS_H_INCLUDED__
#define __HistogramNS_H_INCLUDED__

//=================================
// forward declared dependencies

//=================================
// included dependencies


//=================================
// the actual class
class HistogramNS{
	
	public:
	

	HistogramNS(): nBins(0), nObs(0), Tail(0.0), Top(0.0), incr(0.0), medIndex(0), NSfrac(0.5){}

	void Initialize(unsigned long long Nbins, long double Taile, long double Tope){
		nBins=Nbins;
		Bin= new unsigned long long[nBins];
	//	dBin= new long double[nBins];
		
		Tail = Taile;
		Top = Tope;
		incr = (Top-Tail)/((long double)(nBins));
		for(unsigned long long i=0;i<nBins;++i){
			Bin[i] = 0ULL;
		//	dBin[i] = 0.0;
			
		}
		return;
	}
	void Push(long double val){
		
		if(val>=Tail && val < Top){
			unsigned long long Ind = (unsigned long long)((val - Tail)/incr);
		//cout << "value " << val << " is being pushed into bin number " << Ind << endl;
			Bin[Ind]++;
			
		}
		return;
		
	}

	void PushObs(long double Obs, unsigned oIndex, long double En){
		
		if(En>=Tail && En < Top){
			unsigned long long Ind = (unsigned long long)((En - Tail)/incr);
			oBin[Ind][oIndex] += Obs;
			
		}
		return;

	}
	void Normalize(void){
				
		long double Total =0.0;
		
		for(unsigned long long i=0;i<nBins;++i){
			Total += ((long double)(Bin[i]))*incr;
			
		}
	//	cout << "Total is " << Total;
		for(unsigned long long i=0;i<nBins;++i){
			dBin[i] = ((long double)(Bin[i]))/((long double)(Total));
		//	cout << "dBin[i] " << dBin[i] << endl;
			
		}
		return;
	}

	void Dump(string filename){
		ofstream outfile( filename.c_str());
		
		for(unsigned long long i=0;i<nBins;i++){
			long double pop = Tail+((double)(i+1)*incr); 
			outfile << pop << " " << Bin[i] << endl;
		}
		outfile.close();
		return;
	}

	long double GetNormal(int i){
		return(dBin[i]);
	}
	long double GetIncrement(void){
		return(incr);
	}
	//Need to run Normalize before calling Median
	long double Median(void){
		long double median;
		unsigned long long m_n=1;
		long double area=0.0;
		for(unsigned long long i=0; i<nBins;++i){
			area += incr*dBin[i];
			//cout << "incr: " << incr << " dbin[i] " << dBin[i] << " i: " << i << endl;
			if(area >=NSfrac){
				
				m_n+=i;
				medIndex=i;
				break;
			}
		}
	
		median=Tail+(((long double)(m_n))*incr);
		
		return(median);
	}
	
	void Reset(void){
		for(unsigned long long i=0;i<nBins;++i){
			Bin[i] = 0ULL;
		//	dBin[i] = 0.0;
			for (unsigned int j = 0; j < nObs; ++j)
			{
				oBin[i][j]=0.0;
			}
			
		}
		return;
	}

	void SmartReset(void){
	//	cout << "Initial reset medIndex: "<<medIndex<<endl;
		
		if (medIndex==0)
		{
			this->Reset();
		}
		else{
			unsigned long long rindex;
			if(nBins-medIndex<=20){
				rindex = nBins;
			}
			else{
				rindex = medIndex+10;
			}
		//	cout<<"Going to use rindex of: "<<rindex<<endl;
			for(unsigned long long i=0;i<rindex;++i){
				Bin[i] = 0ULL;
				//dBin[i] = 0.0;
				for (unsigned int j = 0; j < nObs; ++j)
				{
					oBin[i][j]=0.0;
				}
				
			}
		}
		return;
	}
	
	void Renormalize(void){
				
		long double Total =0.0;
		
		for(unsigned long long i=0;i<nBins;++i){
			Total += ((long double)(Bin[i]))*incr;
			dBin[i] = 0.0;
		}
		for(unsigned long long i=0;i<nBins;++i){
			dBin[i] = ((long double)(Bin[i]))/(Total);
			
		}
		return;
	}

	void DumpNormal(string filename){
		ofstream outfile( filename.c_str());
		for(unsigned long long i=0;i<nBins;i++){
			
				long double pop = Tail+((long double)(i+1)*incr); 
				outfile << pop << " " << dBin[i] << endl;
			
		}
		outfile.close();
		return;
	}

	void InitializeObservables(unsigned Nobs){
		nObs = Nobs;
		oBin = new long double*[nBins];
		for (unsigned long long  i = 0; i < nBins; ++i)
		{
			oBin[i] = new long double[Nobs];
			for (unsigned j = 0; j < Nobs; ++j)
			{
				oBin[i][j] = 0.0;
			}
		}
	
	}

	long double GetObsAvg(unsigned oIndex){
		
		long double Total = 0.0;
		unsigned long long Number=0;
		for (unsigned long long i = medIndex+1; i < nBins; ++i)
		{
			Total += oBin[i][oIndex];
			Number += Bin[i];
		}
		long double Oavg = Total / (long double)(Number);
		if(Number == 0){
			Oavg = 0.0;
		}
		return(Oavg);
	}

	void ResetObs(void){
		for(unsigned long long i=0;i<nBins;++i){
			
			for (unsigned int j = 0; j < nObs; ++j)
			{
				oBin[i][j]=0.0;
			}
			
		}
		return;
	}

	void SetNSfrac(long double f){
		NSfrac = f;
		return;
	}

	long double NestedEnergy(void){
		long double Total =0.0;
		unsigned long zerocount=0;
		unsigned long long topBin;
		if(medIndex==0){
			topBin = nBins;
		}
		else{
			topBin = medIndex+1;
		}
		
		for(unsigned long long i=0;i<topBin;++i){
			long double current = ((long double)(Bin[i]))*incr;
			Total += current;
			if (Bin[i]==0 && i==zerocount)
			{
				++zerocount;
			}
			
		}
//		cout<<"Zerocount is: "<<zerocount<<endl;
//		exit(0);	
		long double nse;
		unsigned long long m_n=1;
		long double area=0.0;
		
		if (zerocount>=10)
		{
			zerocount-=10;
		}
		for(unsigned long long i=zerocount; i<topBin;++i){
			long double curr=((long double)(Bin[i]))/(Total);
			area += incr*curr;
			
			if(area >=NSfrac){
				
				m_n+=i;
				medIndex=i;
				break;
			}
		}
	
		nse=Tail+(((long double)(m_n))*incr);
	//	Top = nse;
	//	incr = (Top-Tail)/((long double)(nBins));
		return(nse);
	}
	private:
	unsigned long long nBins;
	unsigned long long* Bin;
	long double* dBin;
	long double** oBin;
	unsigned nObs;
	long double Tail, Top;
	long double incr;
	unsigned long long medIndex;
	long double NSfrac;
};


#endif // __MYCLASS_H_INCLUDED__ 
