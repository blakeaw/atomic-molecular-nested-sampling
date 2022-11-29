//=================================
// include guard
#ifndef __Histogram_H_INCLUDED__
#define __Histogram_H_INCLUDED__

//=================================
// forward declared dependencies

//=================================
// included dependencies
#define MAXBIN 500000

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
		//pushcount
		pushcount=0;
		for(unsigned long long i=0;i<nBins;++i){
			Bin[i] = 0U;
			dBin[i] = 0.0;
		}
		return;
	}

	void Initialize(long double bw, long double Taile, long double Tope){
		nBins=(unsigned long long)((Tope-Taile)/bw);
		nBins+=1;
        Tail = Taile;
		Top = Tope;
		incr = bw;
        if (nBins==0){
            std::cout<<"Initializing range for Histogram leads to zero bins!"<<endl;
        }
        else if(nBins>MAXBIN){

            std::cout<<"Warning: Initializing range for Histogram leads to nBins larger than global max: MAXBIN !"<<endl;
            std::cout<<" Resetting nBins to MAXBIN: "<< MAXBIN <<" and adjusting bin width accordinly."<<endl;
            nBins= MAXBIN ;
            incr = (Top-Tail)/((long double)(nBins));
            std::cout<<" New bin width is: "<<incr<<endl;
        }
		Bin= new unsigned long long[nBins];
		dBin= new long double[nBins];
		
		//pushcount
		pushcount=0;
		for(unsigned long long i=0;i<nBins;++i){
			Bin[i] = 0ULL;
			dBin[i] = 0.0;
			
		}
		return;
	}
	// Add values to the histogram
	void Push(long double val){
		
		if(val>=Tail && val < Top){
			++pushcount;
			unsigned long long Ind = (unsigned long long)((val - Tail)/incr);
//			
			Bin[Ind]++;
			
		}
		return;
		
	}
    	// Add values to the histogram -- if val is greater than the top, it 
        // is added to the top bin
	void PushOverTop(long double val){
		
		if(val>=Tail && val < Top){
			++pushcount;
			unsigned long long Ind = (unsigned long long)((val - Tail)/incr);
//			
			Bin[Ind]++;
			
		}
        if (val>Top){
            ++pushcount;
            Bin[nBins-1]++;
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
		// Returns the Un-Normalized bin value at bin i
	long double GetValueInBin(int i){
		return(Bin[i]);
	}
	// Returns the Normalized bin value at bin i
	long double GetNormal(int i){
		return(dBin[i]);
	}
	// Retuens the bin width
	long double GetIncrement(void){
		return(incr);
	}

	unsigned long long GetNumberOfBins(void){

		return nBins;
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
				//cout<<"area: "<<area<<endl;
				m_n+=i;
				break;
			}
		}
	
		median=Tail+(((long double)(m_n))*incr);
		
		return(median);
	}

	long double IntegrateToFraction(long double f){
		long double Valf;
		unsigned long long m_n=1;
		long double area=0.0;
		for(unsigned long long i=0; i<nBins;++i){
			area += incr*dBin[i];
			
			if(area >=f){
				//cout<<"area: "<<area<<endl;
				m_n+=i;
				break;
			}
		}
	
		Valf=Tail+(((long double)(m_n))*incr);
		
		return(Valf);
	}
	// Reset all bin to zero
	void Reset(void){
		for(unsigned long long i=0;i<nBins;++i){
			Bin[i] = 0U;
			dBin[i] = 0.0;
			pushcount=0;
		}
		return;
	}
	// Renormalization - can be used to normalize using 
	// dbin counts if weighted samples were added. 
	void Renormalize(void){
				
		long double Total =0.0;
		
		for(unsigned long long i=0;i<nBins;++i){
			Total += ((long double)(dBin[i]))*incr;
			//dBin[i] = 0.0;
		}
		if(pushcount==0){
			return;
		}
		for(unsigned long long i=0;i<nBins;++i){
			dBin[i]/= Total;
			
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

	void DumpNorm(ofstream& fileobj){
		
		
		for(unsigned long long i=0;i<nBins;i++){
			long double pop = Tail+((double)(i+1)*incr); 
			fileobj << pop << " " << dBin[i] << endl;
		}
		fileobj<<" "<<endl;
		//outfile.close();
		return;
	}

	~Histogram(){
		if(Bin){
			delete [] Bin;
		}
		if(dBin){
			delete [] dBin;
		}

	}

	unsigned long long GetPushCount(void){
		return pushcount;

	}


	unsigned long long GetPushBinFor(long double val){
		
		if(val>=Tail && val < Top){
			unsigned long long Ind = (unsigned long long)((val - Tail)/incr);
//			
			return Ind;
			
		}
		else{
			return 0;
		}
	}

	long double IntegralFractionAt(long double Value){
		Normalize();
		unsigned long long Ind = (unsigned long long)((Value - Tail)/incr);
		long double frac=0.0;
		if(Ind+1>nBins){
			frac=1.0;
		}
		else{
			for (unsigned int i = 0; i < Ind+1; ++i)
			{
				frac+=dBin[i]*incr;
			}
		}
		return(frac);
	}

	void MultiplyBinNormalByFactor(long double factor, unsigned binnum){
		dBin[binnum]*=factor;

	}

	void AddWithWeight(Histogram& other, long double weight){
		if (other.GetNumberOfBins()!=nBins){
			cout<<"Unable to add histograms with differing bin counts"<<endl;
			return;
		}
		for (unsigned int i = 0; i < nBins; ++i)
		{
			//if (other.GetNormal(i)>0.0){
				dBin[i]+=other.GetNormal(i)*weight;
			//}
			
		}

	}

	void ZeroBelow(long double val){
		unsigned long long Ind = (unsigned long long)((val - Tail)/incr);
		for (unsigned int i = 0; i < Ind+1; ++i)
		{
			dBin[i]=0.0;
		}
		
	}

	void ZeroAbove(long double val){
		unsigned long long Ind = (unsigned long long)((val - Tail)/incr);
		for (unsigned int i = Ind; i < nBins; ++i)
		{
			dBin[i]=0.0;
		}
		
	}

	long double NormTotal(void){
		long double Total=0.0;
		for(unsigned long long i=0;i<nBins;++i){
			Total += ((long double)(Bin[i]))*incr;
		}
		return Total;	
	}
	void NormalizeBy(long double Value){
		
		for(unsigned long long i=0;i<nBins;++i){
			dBin[i] = ( ((long double)(Bin[i])) )/(Value);
			
		}
	}

	void PushToNormWithWeight(long double val, long double weight){

		if(val>=Tail && val < Top){
			++pushcount;
			unsigned long long Ind = (unsigned long long)((val - Tail)/incr);
//			
			dBin[Ind]+=weight;
			
		}
		return;
	}

	void AddWithWeightToUpTo(Histogram& other, long double weight, long double Val){
		
		unsigned long long Ind = (unsigned long long)((Val - Tail)/incr);
		if (Ind>nBins){
			Ind=nBins;
		}
		for (unsigned int i = 0; i < Ind; ++i)
		{
			//if (other.GetNormal(i)>0.0){
			long double bval = ((long double)(i)*incr)+Tail;
			long double w = (long double)(Bin[i])*weight;
				other.PushToNormWithWeight(bval, w);
			//}

		}

	}
	long double TotalWeight(void){
			long double Total=0.0;
		for (unsigned int i = 0; i < nBins; ++i)
		{
			Total+=dBin[i];
		}
		return Total;
	}
	private:
	unsigned long long int nBins;
	unsigned long long int* Bin;
	long double* dBin;
	long double Tail, Top;
	long double incr;
	unsigned long long pushcount;
};


#endif // __MYCLASS_H_INCLUDED__ 
