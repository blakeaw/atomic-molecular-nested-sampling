//=================================
// include guard
#ifndef __Angle_H_INCLUDED__
#define __Angle_H_INCLUDED__

//=================================
// forward declared dependencies


//=================================
// included dependencies


//=================================
// the actual class
class Angle {
	public:
		
		long double angle;
		long double constant;
		unsigned Bindex1;
		unsigned Bindex2;
		unsigned Bindex3;
		
		
		//Assign Default values for values
		Angle(): angle(90.0), constant(0.0){};
};

#endif // __Angle_H_INCLUDED__ 
