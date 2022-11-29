//=================================
// include guard
#ifndef __HardBox_H_INCLUDED__
#define __HardBox_H_INCLUDED__

//=================================
// forward declared dependencies

class Atom;
//=================================
// included dependencies


//=================================
// the actual class

//3d Hard box. Used to encapsulate atoms and serve as boundary. 
class HardBox {
	public:
		
		long double* length;
		long double* center;
		// constructor - empty
		HardBox() {
			length = new long double[3];
			center = new long double[3];
			for (unsigned int i = 0; i < 3; ++i)
			{
				length[i]=1.0;
				center[i]=0.0;
			}
		}
	
	//Destructor
		~HardBox(){
			
			if (center){
				
				delete[] center;
				
			}
			if (length){
				delete[] length;
			}
			
		}
		

	void SetLengths(long double x, long double y, long double z){
			length[0]=x;
			length[1]=y;
			length[2]=z;
		
	}

	void SetCenter(long double x, long double y, long double z){
			center[0]=x;
			center[1]=y;
			center[2]=z;
		
	}
	bool CheckBoundary(const Atom& a){
				long double ax, ay, az;
				ax = a.x, ay = a.y, az = a.z;
				long double dx, dy, dz;
				dx = ax-center[0], dy = ay-center[1], dz = az-center[2];
				//boundaries
				long double xh = center[0]+length[0]/2.0;
				long double xl = center[0]-length[0]/2.0;
				if (dx>xh || dx<xl){
					return 0;
				}
				long double yh = center[1]+length[1]/2.0;
				long double yl = center[1]-length[1]/2.0;
				if (dy>yh || dy<yl){
					return 0;
				}
				long double zh = center[2]+length[2]/2.0;
				long double zl = center[2]-length[2]/2.0;
				if (dz>zh || dz<zl){
					return 0;
				}
				
				// all boundaries passed so return 1 for success - Atom is inside the box
				return 1;
		}
		
		
};

#endif // __Sphere_H_INCLUDED__ 
