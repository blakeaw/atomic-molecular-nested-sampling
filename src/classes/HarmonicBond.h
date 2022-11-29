//=================================
// include guard
#ifndef __HarmonicBond_H_INCLUDED__
#define __HarmonicBond_H_INCLUDED__

//=================================
// forward declared dependencies

//class CGLipid;
//=================================
// included dependencies
#include <stdlib.h>
#include <math.h>
#include <cmath>

//=================================
// the actual class
class HarmonicBond {
	public:
		
		long double length;
		long double constant;
		int Bindex1;
		int Bindex2;
		CGLipid * molec;
		
		//Assign Default values for values
		HarmonicBond(): length(1.0), constant(0.0){}
	
		void SetMolec(CGLipid * m);

		void SetParams(long double k, long double l);

		long double CalcPotential(void);
		long double CalcPotential(long double r);
};

#endif // __Bond_H_INCLUDED__ 
