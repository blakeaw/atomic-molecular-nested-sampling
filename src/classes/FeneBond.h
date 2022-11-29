//=================================
// include guard
#ifndef __FeneBond_H_INCLUDED__
#define __FeneBond_H_INCLUDED__

//=================================
// forward declared dependencies
class CGLipid;

//=================================
// included dependencies


//=================================
// the actual class
class FeneBond {
	public:
		
		long double divlength;
		long double constant;
		long double epsilon;
		long double b;
		int Bindex1;
		int Bindex2;
		CGLipid * molec;
		
		//Assign Default values for values
		FeneBond(): divlength(1.0), constant(0.0), epsilon(1.0), b(1.0){};

		void SetMolec(CGLipid * m);

		void SetParams(long double e, long double k, long double dl, long double bb);

		long double CalcPotential(void);
		
		long double CalcPotential(long double r);
};

#endif // __FeneBond_H_INCLUDED__ 
