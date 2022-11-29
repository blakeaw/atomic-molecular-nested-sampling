//=================================
// include guard
#ifndef __Atom_H_INCLUDED__
#define __Atom_H_INCLUDED__

//=================================
// forward declared dependencies


//=================================
// included dependencies


//=================================
// the actual class
class Atom {
	public:
		
		long double x, y, z;
		long double Eu, dE;
		Atom(): Eu(0.0), dE(0.0) {}
};

#endif // __Atom_H_INCLUDED__ 
