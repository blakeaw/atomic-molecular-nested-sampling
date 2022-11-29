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
		//int Aindex;
		long double x, y, z;
		long double mass;
		//charge tag - 0 negative, 1 positive
		bool ctag;
//		Atom(Atom another) {};
		Atom(): x(0.0), y(0.0), z(0.0), mass(1.0), ctag(1){}
		//Copy Constructor
		Atom(const Atom& other){
			x=other.x;
			y=other.y;
			z=other.z;
			mass=other.mass;
			ctag=other.ctag;
		}
		void Equate(const Atom& other){
			x=other.x;
			y=other.y;
			z=other.z;
			mass=other.mass;
			ctag=other.ctag;
		}
};

#endif // __Atom_H_INCLUDED__ 
