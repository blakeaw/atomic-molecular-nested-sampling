//=================================
// include guard
#ifndef __SphereWall_H_INCLUDED__
#define __SphereWall_H_INCLUDED__

//=================================
// forward declared dependencies

class Config;
//=================================
// included dependencies


//=================================
// the actual class

//Spherical LJ wall. Used to encapsulate atoms and serve as boundary. 
class SphereWall {
	public:
		long double xc, yc, zc;
		long double radius;
		long double sigma, epsilon;
		SphereWall(): xc(0.0), yc(0.0), zc(0.0), radius(7.6), sigma(1.0), epsilon(0.000010){}


		
		
		long double CalcPotential(const Atom& a){
				bool check = CheckBoundary(a);
				if(check == 1){
				long double ax, ay, az;
				ax = a.x, ay = a.y, az = a.z;
				long double asig = a.sig;
				long double aeps = a.eps;
				long double sig = (sigma+asig)/2.0;
				long double eps = sqrt(epsilon*aeps);
				long double dx, dy, dz;
				dx = ax-xc, dy = ay-yc, dz = az-zc;
				long double p = sqrt(dx*dx+dy*dy+dz*dz);
				long double d = radius - p;
				long double potent;
				long double sor = sig/d;
				potent = 4.0*eps*(pow(sor, 12)-pow(sor, 6));
				return(potent);
				}
				else{
					long double potent = 100000.0;
					return(potent);
				}
		}

		long double TotalPotential(const Config& c){

			unsigned long long N = c.natoms;
			long double Etot = 0.0;
			for (unsigned int i = 0; i < N; ++i)
			{
				Etot += CalcPotential(c.atom[i]);
			}
				return(Etot);
	

		}
	

			long double CalcPotentialOC(const Atom& a, const long double& ox, const long double& oy, const long double& oz){

				long double ax, ay, az;
				ax = ox, ay = oy, az = oz;
				long double asig = a.sig;
				long double aeps = a.eps;
				long double sig = (sigma+asig)/2.0;
				long double eps = sqrt(epsilon*aeps);
				long double dx, dy, dz;
				dx = ax-xc, dy = ay-yc, dz = az-zc;
				long double p = sqrt(dx*dx+dy*dy+dz*dz);
				long double d = radius - p;
				long double potent;
				long double sor = sig/d;
				potent = 4.0*eps*(pow(sor, 12)-pow(sor, 6));
				return(potent);			


		}

		long double DeltaE(const Atom& a, const long double& ox, const long double& oy, const long double& oz){
				long double En, Eo;
				Eo = CalcPotentialOC(a, ox, oy, oz);
				En = CalcPotential(a);
				long double dE = En - Eo;
				return(dE);
		}

	
		
		private:

			bool CheckBoundary(const Atom& a){
				long double ax, ay, az;
				ax = a.x, ay = a.y, az = a.z;
				long double dx, dy, dz;
				dx = ax-xc, dy = ay-yc, dz = az-zc;
				long double p = sqrt(dx*dx+dy*dy+dz*dz);
				if(p>=radius){
					return 0;
				}
				else{
					return 1;
				}
		}
		
		
};

#endif // __Sphere_H_INCLUDED__ 
