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
	
	// Initialize empty histogram -- bin increments stored in Bin, dBin is for normalized values
	void Initialize(unsigned long long Nbins, long double Taile, long double Tope){
		// Number of bins
		nBins=Nbins;
		// to store bin increments
		Bin= new unsigned long long int[nBins];
		// used in normalization to store normalized values
		dBin= new long double[nBins];
		// Lower bound
		Tail = Taile;
		// upper bound
		Top = Tope;
		// Bin width
		incr = (Top-Tail)/((long double)(nBins));
		for(unsigned long long i=0;i<nBins;++i){
			Bin[i] = 0U;
			dBin[i] = 0.0;
		}
		return;
	}

	void Initialize(long double bw, long double Taile, long double Tope){
		nBins=(unsigned long long)((Tope-Taile)/bw);
		nBins+=1;
		Bin= new unsigned long long[nBins];
		dBin= new long double[nBins];
		Tail = Taile;
		Top = Tope;
		incr = bw;
		for(unsigned long long i=0;i<nBins;++i){
			Bin[i] = 0ULL;
			dBin[i] = 0.0;
			
		}
		return;
	}
	// Add values to the histogram
	void Push(long double val){
		
		if(val>=Tail && val < Top){
			unsigned long long Ind = (unsigned long long)((val - Tail)/incr);
//			
			Bin[Ind]++;
			
		}
		return;
		
	}
	// Normalize the histogram - integral normalization
	// Normalized values are stored in dBin
	void Normalize(void){
				
		long double Total =0.0;
		
		for(unsigned long long i=0;i<nBins;++i){
			Total += ((long double)(Bin[i]))*incr;
			
		}
		
		for(unsigned long long i=0;i<nBins;++i){
			dBin[i] = ((long double)(Bin[i]))/((long double)(Total));
			
		}
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

	// Creates a new file and dumps the current contents of the histogram Bins to it 
	void Dump(string filename){
		ofstream outfile( filename.c_str());
		
		for(unsigned long long i=0;i<nBins;i++){
			long double pop = Tail+((double)(i+1)*incr); 
			outfile << pop << " " << Bin[i] << endl;
		}
		outfile.close();
		return;
	}
	void DumpNorm(ofstream& fileobj){
		
		
		for(unsigned long long i=0;i<nBins;i++){
			long double pop = Tail+((double)(i+1)*incr); 
			fileobj << pop << " " << dBin[i] << endl;
		}
		fileobj<<" "<<endl;
		//outfile.close();
		return;
	}
	// Returns the Normalized bin value at bin i
	long double GetNormal(int i){
		return(dBin[i]);
	}
	// Retuens the bin width
	long double GetIncrement(void){
		return(incr);
	}
	//Need to run Normalize before calling Median
	// Median integrates the normalized values up to a fraction of 1/2
	// and returns the associated bin value
	long double Median(void){
		long double median;
		unsigned long long m_n=1;
		long double area=0.0;
		for(unsigned long long i=0; i<nBins;++i){
			area += incr*dBin[i];
			
			if(area >=0.50){
				
				m_n+=i;
				break;
			}
		}
	
		median=Tail+(((long double)(m_n))*incr);
		
		return(median);
	}
	// Reset all bin to zero
	void Reset(void){
		for(unsigned long long i=0;i<nBins;++i){
			Bin[i] = 0U;
			dBin[i] = 0.0;
		}
		return;
	}
	// Renormalization - can be used after the first instance of 
	// of Normalize (e.g. in a loop) -- currently the same as 
	// Normalize
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
	// Creates a new file and dumps the current contents of the histogram dBins to it 
	void DumpNormal(string filename){
		ofstream outfile( filename.c_str());
		for(unsigned long long i=0;i<nBins;i++){
			
				long double pop = Tail+((long double)(i+1)*incr); 
				outfile << pop << " " << dBin[i] << endl;
			
		}
		outfile.close();
		return;
	}
	
	unsigned long long GetNumberOfBins(void){

		return nBins;
	}

	private:
	unsigned long long int nBins;
	unsigned long long int* Bin;
	long double* dBin;
	long double Tail, Top;
	long double incr;
};


#endif // __MYCLASS_H_INCLUDED__ 
