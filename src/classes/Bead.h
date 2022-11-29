//=================================
// include guard
#ifndef __Bead_H_INCLUDED__
#define __Bead_H_INCLUDED__

//=================================
// forward declared dependencies
class CGLipid;

//=================================
// included dependencies
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/CGLipid.h"

//=================================
// the actual class
class Bead {
	public:
//		string type;
		long double x, y, z;
		long double mass;
		// 0 -> head bead ; 1 -> tail bead
		bool type;
		CGLipid * molec;
		
		//Assign Default values for values
		Bead(): x(0.0), y(0.0), z(0.0), mass(1.0){};
	
		void SetMolec(CGLipid * m);
};

#endif // __Bead_H_INCLUDED__ 
