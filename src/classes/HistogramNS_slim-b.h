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
	

	HistogramNS(): nBins(0), nObs(0), Tail(0.0), Top(0.0), incr(0.0), medIndex(0), NSfrac(0.5), prevTot(0.0), zcprev(0), obin(0){}

	void Initialize(unsigned long long Nbins, double Taile, double Tope){
		nBins=Nbins;
		Bin= new unsigned long long[nBins];
		prevTop = nBins;
		Tail = Taile;
		Top = Tope;
		incr = (Top-Tail)/((double)(nBins));
		for(unsigned long long i=0;i<nBins;++i){
			Bin[i] = 0ULL;
			
			
		}
		return;
	}

	void Initialize(double bw, double Taile, double Tope){
		nBins=(unsigned long long)((Tope-Taile)/bw);
		nBins+=1;
		Bin= new unsigned long long[nBins];
		prevTop = nBins;
		Tail = Taile;
		Top = Tope;
		incr = bw;
		for(unsigned long long i=0;i<nBins;++i){
			Bin[i] = 0ULL;
			
			
		}
		return;
	}
	void Push(double val){
		
		if(val>Tail && val < Top){
			unsigned long long Ind = (unsigned long long)((val - Tail)/incr);
		//cout << "value " << val << " is being pushed into bin number " << Ind << endl;
			Bin[Ind]++;
			
		}
		return;
		
	}

	void PushObs(double Obs, unsigned oIndex, double En){
		
		if(En>Tail && En < Top){
			unsigned long long Ind = (unsigned long long)((En - Tail)/incr);
			oBin[Ind][oIndex] += Obs;
			
		}
		return;

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
				rindex = medIndex+20;
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
	

	void InitializeObservables(unsigned Nobs){
		nObs = Nobs;
		oBin = new double*[nBins];
		for (unsigned long long  i = 0; i < nBins; ++i)
		{
			oBin[i] = new double[Nobs];
			for (unsigned j = 0; j < Nobs; ++j)
			{
				oBin[i][j] = 0.0;
			}
		}
		obin=1;		
		return;
	}

	

	void SetNSfrac(double f){
		NSfrac = f;
		return;
	}

	double NestedEnergy(void){
		double Total =0.0;
		unsigned long zerocount=0;
		
		if(medIndex==0){
			prevTop = nBins;
		}
		else{
			prevTop = medIndex+1;
		}
		
		for(unsigned long long i=0;i<prevTop;++i){
			double current = ((double)(Bin[i]))*incr;
			Total += current;
			if (Bin[i]==0 && i==zerocount)
			{
				++zerocount;
			}
			
		}
//		cout<<"Zerocount is: "<<zerocount<<endl;
//		exit(0);	
		double nse;
		unsigned long long m_n=1;
		double area=0.0;
		
		if (zerocount>=10)
		{
			zerocount-=10;
		}
		prevTot=Total;
		zcprev=zerocount;
		for(unsigned long long i=zerocount; i<prevTop;++i){
			double curr=((double)(Bin[i]))/(Total);
			area += incr*curr;
			
			if(area >=NSfrac){
				
				m_n+=i;
				medIndex=i;
				break;
			}
		}
	
		nse=Tail+(((double)(m_n))*incr);
		//Top = nse;
		//incr = (Top-Tail)/((double)(nBins));
		cout<<"nested value at area "<<area<<endl;
		return(nse);
	}

	void DumpNested(string& name){
		ofstream outfile(name.c_str());
		//cout<<"zcprev "<<zcprev<<" prevTot "<<prevTot<<endl;
		outfile<< setiosflags(ios::fixed) << setprecision(15);
		for(unsigned long long i=zcprev; i<prevTop;i+=1){
			double curr= ((double)(Bin[i]))/(prevTot);
			curr*=incr;
			double value = Tail+incr*(i+1);
			outfile<<value<<" "<<curr<<endl;
			//cout<<"value "<<value<<" curr "<<curr<<" total "<<prevTot<<endl;
		}

		outfile.close();
		return;
	}

	void Dump(ofstream& fileobj){
		
		
		for(unsigned long long i=0;i<nBins;i++){
			double pop = Tail+((double)(i+1)*incr); 
			fileobj << pop << " " << Bin[i] << endl;
		}
		fileobj<<" "<<endl;
		//outfile.close();
		return;
	}

	unsigned long long GetNumberOfBins(void){

		return nBins;
	}

	double GetObsAvg(unsigned oIndex){
		
		double Total = 0.0;
		unsigned long long Number=0;
		//bool failflag=0;
		for (unsigned long long i = medIndex+1; i < prevTop; ++i)
		{
			Total += oBin[i][oIndex];
			Number += Bin[i];
		}
		double Oavg = Total / (double)(Number);
	//	if(Oavg<0.0001 && Oavg>-0.0001){
	//		cout<<"Flag1!!!---Number is: "<<Number<<" and Total is: "<<Total<<endl;
	//		cout<<"----medIndex is: "<<medIndex<<" and prevTop: "<<prevTop<<endl;
	//		cout<<"----oIndex is: "<<oIndex<<" and Oavg is: "<<Oavg<<endl;
	//		failflag=1;
	//	}
		if(Number == 0){
			Oavg = 0.0;
			cout<<"Flag2!!!---Number is: "<<Number<<" and Total is: "<<Total<<endl;
			cout<<"----medIndex is: "<<medIndex<<" and prevTop: "<<prevTop<<endl;
			//failflag=1;
			
			unsigned long long topbin = prevTop;
			if(prevTop==nBins){
				topbin-=1;
			}
			Total = oBin[medIndex][oIndex]+oBin[topbin][oIndex];
			Number = Bin[medIndex]+Bin[topbin];
			Oavg = Total / (double)(Number);
		}
//		if(failflag==1){
//			unsigned long long topbin = prevTop;
//			if(prevTop==nBins){
//				topbin-=1;
//			}
//			Total = oBin[medIndex][oIndex]+oBin[topbin][oIndex];
//			Number = Bin[medIndex]+Bin[topbin];
//			Oavg = Total / (double)(Number);
		//	cout<<"--Oavg values were adjusted to compensate---"<<endl;
		//	cout<<"---Number is now: "<<Number<<" and Total is now: "<<Total<<endl;
		//	cout<<"----medIndex is now: "<<medIndex<<" and prevTop is now: "<<prevTop<<endl;
		//	cout<<"----oIndex is now: "<<oIndex<<" and Oavg is now: "<<Oavg<<endl;
//		}
		
		
		return(Oavg);
	}
	/* Calculates the average of all the observables simultaneously - reads in a reference variable and 
       will save the final average in the ref variable
	   As is will need to update if new/more observables are added or if some are removed
	*/
	void GetObsAvgAll(double& ot0, double& ot1, double& ot2, double& ot3, double& ot4, double& ot5, double& ot6, double& ot7, double& ot8, double& ot9 ){
		
		double t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;
		unsigned long long Number=0;
		t0=0.0, t1=0.0, t2=0.0, t3=0.0, t4=0.0, t5=0.0, t6=0.0;
		t7=0.0, t8=0.0, t9=0.0;
		//bool failflag=0;
		for (unsigned long long i = medIndex+1; i < prevTop; ++i)
		{
			t0 += oBin[i][0];
			t1 += oBin[i][1];
			t2 += oBin[i][2];
			t3 += oBin[i][3];
			t4 += oBin[i][4];
			t5 += oBin[i][5];
			t6 += oBin[i][6];
			t7 += oBin[i][7];
			t8 += oBin[i][8];
			t9 += oBin[i][9];
			Number += Bin[i];
		}
		double Num = (double)(Number);
		ot0 = t0 / Num;
		ot1 = t1 / Num;
		ot2 = t2 / Num;
		ot3 = t3 / Num;
		ot4 = t4 / Num;
		ot5 = t5 / Num;
		ot6 = t6 / Num;
		ot7 = t7 / Num;
		ot8 = t8 / Num;
		ot9 = t9 / Num;

		if(Number == 0){
			//Oavg = 0.0;
			//cout<<"Flag2!!!---Number is: "<<Number<<" and Total is: "<<Total<<endl;
			cout<<"Flag2!!!---"<<endl;
			cout<<"----medIndex is: "<<medIndex<<" and prevTop: "<<prevTop<<endl;
			//failflag=1;
			
			unsigned long long topbin = prevTop;
			if(prevTop==nBins){
				topbin-=1;
			}
			t0 = oBin[medIndex][0]+oBin[topbin][0];
			t1 = oBin[medIndex][1]+oBin[topbin][1];
			t2 = oBin[medIndex][2]+oBin[topbin][2];
			t3 = oBin[medIndex][3]+oBin[topbin][3];
			t4 = oBin[medIndex][4]+oBin[topbin][4];
			t5 = oBin[medIndex][5]+oBin[topbin][5];
			t6 = oBin[medIndex][6]+oBin[topbin][6];
			t7 = oBin[medIndex][7]+oBin[topbin][7];
			t8 = oBin[medIndex][8]+oBin[topbin][8];
			t9 = oBin[medIndex][9]+oBin[topbin][9];
			Number = Bin[medIndex]+Bin[topbin];
			Num = (double)(Number);
			ot0 = t0 / (double)(Num);
			ot1 = t1 / (double)(Num);
			ot2 = t2 / (double)(Num);
			ot3 = t3 / (double)(Num);
			ot4 = t4 / (double)(Num);
			ot5 = t5 / (double)(Num);
			ot6 = t6 / (double)(Num);
			ot7 = t7 / (double)(Num);
			ot8 = t8 / (double)(Num);
			ot9 = t9 / (double)(Num);
		}

		
		
		return;
	}
	void GetObsAvgAll(double& ot0, double& ot1, double& ot2, double& ot3, double& ot4, double& ot5, double& ot6){
		
		double t0, t1, t2, t3, t4, t5, t6;
		unsigned long long Number=0;
		t0=0.0, t1=0.0, t2=0.0, t3=0.0, t4=0.0, t5=0.0, t6=0.0;
		
		//bool failflag=0;
		for (unsigned long long i = medIndex+1; i < prevTop; ++i)
		{
			t0 += oBin[i][0];
			t1 += oBin[i][1];
			t2 += oBin[i][2];
			t3 += oBin[i][3];
			t4 += oBin[i][4];
			t5 += oBin[i][5];
			t6 += oBin[i][6];
			
			Number += Bin[i];
		}
		double Num = (double)(Number);
		ot0 = t0 / Num;
		ot1 = t1 / Num;
		ot2 = t2 / Num;
		ot3 = t3 / Num;
		ot4 = t4 / Num;
		ot5 = t5 / Num;
		ot6 = t6 / Num;
	

		if(Number == 0){
			//Oavg = 0.0;
			//cout<<"Flag2!!!---Number is: "<<Number<<" and Total is: "<<Total<<endl;
			cout<<"Flag2!!!---"<<endl;
			cout<<"----medIndex is: "<<medIndex<<" and prevTop: "<<prevTop<<endl;
			//failflag=1;
			
			unsigned long long topbin = prevTop;
			if(prevTop==nBins){
				topbin-=1;
			}
			t0 = oBin[medIndex][0]+oBin[topbin][0];
			t1 = oBin[medIndex][1]+oBin[topbin][1];
			t2 = oBin[medIndex][2]+oBin[topbin][2];
			t3 = oBin[medIndex][3]+oBin[topbin][3];
			t4 = oBin[medIndex][4]+oBin[topbin][4];
			t5 = oBin[medIndex][5]+oBin[topbin][5];
			t6 = oBin[medIndex][6]+oBin[topbin][6];
			
			Number = Bin[medIndex]+Bin[topbin];
			Num = (double)(Number);
			ot0 = t0 / (double)(Num);
			ot1 = t1 / (double)(Num);
			ot2 = t2 / (double)(Num);
			ot3 = t3 / (double)(Num);
			ot4 = t4 / (double)(Num);
			ot5 = t5 / (double)(Num);
			ot6 = t6 / (double)(Num);
			
		}

		
		
		return;
	}
	

	private:
	unsigned long long nBins;
	unsigned long long* Bin;
	double** oBin;
	unsigned nObs;
	double Tail, Top;
	double incr;
	unsigned long long medIndex;
	unsigned long long prevTop;
	double NSfrac;
	double prevTot;
	unsigned long long zcprev;
	bool obin;

};


#endif // __MYCLASS_H_INCLUDED__ 
