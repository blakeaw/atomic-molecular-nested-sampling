//=================================
// include guard
#ifndef __HardOpenCylinder_H_INCLUDED__
#define __HardOpenCylinder_H_INCLUDED__

//=================================
// forward declared dependencies
class Atom;
//=================================
// included dependencies


//=================================
// the actual class
class HardOpenCylinder {
	public:
		long double* center;
		long double length;
		long double radius;
		//Constructor - no arguments
		HardOpenCylinder(void):length(1.0), radius(1.0){
			center = new long double[3];
			for (unsigned int i = 0; i < 3; ++i)
			{
				center[i]=0.0;
			}
			
		}
		//Destructor
		~HardOpenCylinder(){
			
			if (center){
				
				delete[] center;
				
			}
			
		}
		//Constructor- takes length (l) and radius (r)
		HardOpenCylinder(long double l, long double r){
			center = new long double[3];
			for (unsigned int i = 0; i < 3; ++i)
			{
				center[i]=0.0;
			}
			length=l;
			radius=r;
		}
		//returns 1 if inside and 0 if not
		bool CheckIfInsideYZAndWrapX(Atom& a){
				long double x,y,z;
				x = a.x;
				y = a.y;
				z= a.z;
				//first check y and z to see if inside
				long double Rs = (y-center[1])*(y-center[1])+(z-center[2])*(z-center[2]);
				if (Rs>radius*radius){
					return 0;
				}
				//y and z inside check x for wrapping
				long double xmin, xmax;
				xmin = center[0]-length/2.0;
				xmax = center[0]+length/2.0;
				if (x<xmin){
					long double dx = xmin - x;
					a.x = xmax-dx;
				}
				else if (x>xmax){
					long double dx = x - xmax;
					a.x = xmin+dx;
				}
				
					return 1;
				
				
		}
		

			//returns 1 if inside and 0 if not
		bool CheckIfInsideYZ(const Atom& a){
				long double y,z;
				y = a.y;
				z= a.z;
				//first check y and z to see if inside
				long double Rs = (y-center[1])*(y-center[1])+(z-center[2])*(z-center[2]);
				if (Rs>radius*radius){
					return 0;
				}
				else {
					return 1;
				}
		}

		long double ComputeVolume(void){
			long double pi=3.14159265359;
			return (pi*radius*radius*length);
		}
};

#endif // __HardOpenCylinder_H_INCLUDED__ 
