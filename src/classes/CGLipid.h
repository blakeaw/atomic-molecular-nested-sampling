//=================================
// include guard
#ifndef __CGLipid_H_INCLUDED__
#define __CGLipid_H_INCLUDED__

//=================================
// forward declared dependencies

//class Bead;
//class HarmonicBond;
//class FeneBond;

//=================================
// included dependencies
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/Bead.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/HarmonicBond.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/FeneBond.h"

//=================================
// the actual class
class CGLipid {
	public:
		const unsigned short nbeads;
		const unsigned short nfbonds;
		const unsigned short nhbonds;
		long double* COM;
		long double sigma;
		long double epsilon;	
		long double sigmahh, sigmaht, sigmatt;
		unsigned index;

		Bead * bead;
		FeneBond * fbond;
		HarmonicBond * hbond;
//		ConfigCGL * config;		
		//Constructor
		CGLipid(): nbeads(3), nfbonds(2), nhbonds(1), sigma(7.3), epsilon(0.596) {
			
			//b values
			//head head
			sigmahh = 0.95*sigma;
			//head tail
			sigmaht = 0.95*sigma;
			//tail tail
			sigmatt = sigma;

			bead = new Bead[nbeads];
			fbond = new FeneBond[nfbonds];
			hbond = new HarmonicBond[nhbonds];
			COM = new long double[3];
			//Initialize values of COM to zero
			for	(short i =0; i<3 ; ++i){
				COM[i] = 0.0;
			}
			//Fene Bond between head-tail1 and tail1-tail2
			for(short i=0; i<nfbonds; ++i){
				fbond[i].Bindex1=i;
				fbond[i].Bindex2=i+1;
				
			}

			for (short i = 0; i < nfbonds; ++i)
			{
				fbond[i].SetMolec(this);
				long double tsig;
				if(i == 0){
					tsig = sigmaht;
				}
				else{
					tsig = sigmatt;
				}
				fbond[i].SetParams(epsilon, (30.0*epsilon)/(sigma*sigma), 1.5*sigma, tsig);	
			}
			//Harmonic bond between head (0) and 2nd tail (2) bead
			hbond[0].Bindex1=0;
			hbond[0].Bindex2=2;
			hbond[0].SetMolec(this);
			hbond[0].SetParams((10.0*epsilon)/(sigma*sigma), 4.0*sigma);

			for (unsigned int i = 0; i < nbeads; ++i)
			{
				bead[i].SetMolec(this);
				if(i==0){
					bead[i].type = 0;
				}
				else{
					bead[i].type = 1;
				}
			}
		}
		//Destructor
		~CGLipid(){
			delete [] bead;
			delete [] fbond;
			delete [] hbond;
			delete [] COM;
		}
		
	
		//Member Function Prototypes
		void CalcCom(void);
		long double ComDistFrom(CGLipid& one);
		void EquateLipidCoord(const CGLipid& other);
};

#endif // __Dimer_H_INCLUDED__ 
