//=================================
// include guard
#ifndef __ContinuumSphere_H_INCLUDED__
#define __ContinuumSphere_H_INCLUDED__

#define PI 3.14159265359
//=================================
// forward declared dependencies
//class Atom

//=================================
// included dependencies


//=================================
// the actual class

//2-dimensional continuum surface
class ContinuumSphere {

public:
	long double x,y,z;
	long double eps, sig;
	long double density;

	ContinuumSphere(): x(0.0),y(0.0),z(0.0), eps(1.0), R(1.0) {}

	ContinuumSphere(long double e, long double r, long double d){
			x=0.0;
			y=0.0;
			z=0.0;
			eps = e;
			R = r;
			density = d;
	}



	long double CalcPotential(const Atom& a){
				
				//long double az;
				//az = a.z;
				long double asig = a.sig;
				long double aeps = a.eps;
				//long double csig = (sig+asig)/2.0;
				//long double rc = 12.0*csig;
				//long double rcs = rc*rc;
				//long double ceps = sqrt(eps*aeps);
				long double dx, dy, dz;
				dx = a.x-x;
				dy = a.y-y;
				dz = a.z-z;
				long double rs = dx*dx+dy*dy+dz*dz;
				long double r = sqrt(rs);
				long double potent=0.0;
				
					
					//long double sor = (csig*csig)/p;
					//potent = 4.0*ceps*(pow(sor, 6)-pow(sor, 3));
					//potent = 8.0*PI*ceps*density*( (2.0/5.0)*pow(csig, 13.0)*pow(p, -10.0) - pow(csig, 6.0)*pow(p, -4.0) );
				 potent= (16.0*pow(R,3)*density*eps*PI*pow(asig,6)*(15.0*pow(pow(R,2) - pow(r,2),6) - (5.0*pow(R,6) + 45.0*pow(R,4)*pow(r,2) + 63.0*pow(R,2)*pow(r,4) + 15.0*pow(r,6))*pow(asig,6)))/(45.0*pow(pow(R,2) - pow(r,2),9));
				
				//cout<<"p: "<<p<<" potent: "<<potent<<" az "<<az<<" dz "<<dz<<" WallC: "<<WallCoord<<endl;
				return(potent);
				
				
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

				//long double az;
				//az = a.z;
				long double asig = a.sig;
				long double aeps = a.eps;
				//long double csig = (sig+asig)/2.0;
				//long double rc = 12.0*csig;
				//long double rcs = rc*rc;
				//long double ceps = sqrt(eps*aeps);
				long double dx, dy, dz;
				dx = ox-x;
				dy = oy-y;
				dz = oz-z;
				long double rs = dx*dx+dy*dy+dz*dz;
				long double r = sqrt(rs);
				long double potent=0.0;
				
					
					//long double sor = (csig*csig)/p;
					//potent = 4.0*ceps*(pow(sor, 6)-pow(sor, 3));
					//potent = 8.0*PI*ceps*density*( (2.0/5.0)*pow(csig, 13.0)*pow(p, -10.0) - pow(csig, 6.0)*pow(p, -4.0) );
				 potent= (16.0*pow(R,3)*density*eps*PI*pow(asig,6)*(15.0*pow(pow(R,2) - pow(r,2),6) - (5.0*pow(R,6) + 45.0*pow(R,4)*pow(r,2) + 63.0*pow(R,2)*pow(r,4) + 15.0*pow(r,6))*pow(asig,6)))/(45.0*pow(pow(R,2) - pow(r,2),9));
				
				//cout<<"p: "<<p<<" potent: "<<potent<<" az "<<az<<" dz "<<dz<<" WallC: "<<WallCoord<<endl;
				return(potent);


		}

		long double DeltaE(const Atom& a, const long double& ox, const long double& oy, const long double& oz){
				long double En, Eo;
				Eo = CalcPotentialOC(a, ox, oy, oz);
				En = CalcPotential(a);
				long double dE = En - Eo;
				return(dE);
		}

							
};

#endif // __2dContinuumSurface_H_INCLUDED__ 
