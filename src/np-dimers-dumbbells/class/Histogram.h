//=================================
// include guard
#ifndef __Histogram_H_INCLUDED__
#define __Histogram_H_INCLUDED__

//=================================
// forward declared dependencies

//=================================
// included dependencies


//=================================
// the actual class
class Histogram{
	
	public:
	

	void Initialize(unsigned long long Nbins, long double Taile, long double Tope){
		nBins=Nbins;
		Bin= new unsigned long long[nBins];
		dBin= new long double[nBins];
		Tail = Taile;
		Top = Tope;
		incr = (Top-Tail)/((long double)(nBins));
		for(unsigned long long i=0;i<nBins;++i){
			Bin[i] = 0ULL;
			dBin[i] = 0.0;
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

	void Dump(ofstream& fileobj){
		
		
		for(unsigned long long i=0;i<nBins;i++){
			long double pop = Tail+((double)(i+1)*incr); 
			fileobj << pop << " " << Bin[i] << endl;
		}
		fileobj<<" "<<endl;
		//outfile.close();
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
			if(area >=0.60){
				
				m_n+=i;
				break;
			}
		}
	
		median=Tail+(((long double)(m_n))*incr);
		
		return(median);
	}
	
	void Reset(void){
		for(unsigned long long i=0;i<nBins;++i){
			Bin[i] = 0ULL;
			dBin[i] = 0.0;
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

	private:
	unsigned long long nBins;
	unsigned long long* Bin;
	long double* dBin;
	long double Tail, Top;
	long double incr;
};


#endif // __MYCLASS_H_INCLUDED__ 
