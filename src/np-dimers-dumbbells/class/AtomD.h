//=================================
// include guard
#ifndef __AtomD_H_INCLUDED__
#define __AtomD_H_INCLUDED__

//=================================
// forward declared dependencies


//=================================
// included dependencies
#include <string>

//=================================
// the actual class
class Atom {
	public:
		string type;
		long double x, y, z;
		long double Ae;
		long double Ac;
		long double rad;
		long double mass;
		long double rsmall;
		long double sigc;
		//charge tag: 1 => positive, 0 => negative
		bool ctag;
		//Assign Default values for values
		Atom(): x(0.0), y(0.0), z(0.0), Ae(1.0), Ac(1.0), rad(1.0), mass(1.0), rsmall(1.0), sigc(50.0), ctag(0){};
};

#endif // __Atom_H_INCLUDED__ 
