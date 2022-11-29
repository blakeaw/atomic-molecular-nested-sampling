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
	long double R;
	ContinuumSphere(): x(0.0),y(0.0),z(0.0), eps(1.0), R(1.0), density(1.0) {}

	ContinuumSphere(long double e, long double r, long double d){
			x=0.0;
			y=0.0;
			z=0.0;
			eps = e;
			R = r;
			density = d;
	}

	//Copy Construct
	ContinuumSphere(const ContinuumSphere& other){
			x=other.x;
			y=other.y;
			z=other.z;
			eps = other.eps;
			R = other.R;
			density = other.density;

	}
	//Destructor
	~ContinuumSphere(){};
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
				// u += (16.0*pow(a,3)*density*epsilon*M_PI*pow(sigma,6)*(15.0*pow(pow(a,2) - pow(R,2),6) - (5.0*pow(a,6) + 45.0*pow(a,4)*pow(R,2) + 63.0*pow(a,2)*pow(R,4) + 15.0*pow(R,6))*pow(sigma,6)))/(45.0*pow(pow(a,2) - pow(R,2),9));
				
				 potent= (16.0*pow(R,3)*density*eps*PI*pow(asig,6)*(15.0*pow(pow(R,2) - pow(r,2),6) - (5.0*pow(R,6) + 45.0*pow(R,4)*pow(r,2) + 63.0*pow(R,2)*pow(r,4) + 15.0*pow(r,6))*pow(asig,6))) / (45.0*pow(pow(R,2) - pow(r,2),9));		
//				long double pref = 9.0*PI*eps*density*pow(R, 3);
//				long double pa = pow(asig, 9)*(3.0*pow(R, 4) + 42.0*R*R*r*r+35.0*pow(r, 4))/(35.0*r*pow(r*r*-R*R, 6));
//				long double pb = pow(asig, 6)/pow(r*r*-R*R, 3);
//				potent = pref*(pa - pb);
				//cout<<"p: "<<p<<" potent: "<<potent<<" az "<<az<<" dz "<<dz<<" WallC: "<<WallCoord<<endl;
				return(potent);
				
				
		}


		long double CalcPotentialDerivative(const Atom& a){
				
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
				long double dudR=0.0;
				
				 dudR += (16.0*density*eps*PI*pow(R,2)*pow(asig,6)*(-5.0*pow(pow(R,2) - pow(r,2),6)*(pow(R,2) + pow(r,2)) +
       (5.0*pow(R,4) + 10.0*pow(R,2)*pow(r,2) + pow(r,4))*(pow(R,4) + 10.0*pow(R,2)*pow(r,2) + 5.0*pow(r,4))*pow(asig,6)))/(5.0*pow(pow(R,2) - pow(r,2),10));
				
				return(dudR);
				
				
		}
	long double TotalPotential(const Config& c){
			//return 0.0;
			unsigned long long N = c.natoms;
			long double Etot = 0.0;
			for (unsigned int i = 0; i < N; ++i)
			{
				Etot += CalcPotential(c.atom[i]);
			}
				return(Etot);
	

		}

	long double TotalPotentialDerivative(const Config& c){
			//return 0.0;
			unsigned long long N = c.natoms;
			long double Etot = 0.0;
			for (unsigned int i = 0; i < N; ++i)
			{
				Etot += CalcPotentialDerivative(c.atom[i]);
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
//				long double pref = 9.0*PI*eps*density*pow(R, 3);
//				long double pa = pow(asig, 9)*(3.0*pow(R, 4) + 42.0*R*R*r*r+35.0*pow(r, 4))/(35.0*r*pow(r*r*-R*R, 6));
//				long double pb = pow(asig, 6)/pow(r*r*-R*R, 3);
//				potent = pref*(pa - pb);
				//cout<<"p: "<<p<<" potent: "<<potent<<" az "<<az<<" dz "<<dz<<" WallC: "<<WallCoord<<endl;
				return(potent);


		}

		long double DeltaE(const Atom& a, const long double& ox, const long double& oy, const long double& oz){
				//return 0.0;
				long double En, Eo;
				Eo = CalcPotentialOC(a, ox, oy, oz);
				En = CalcPotential(a);
				long double dE = En - Eo;
				return(dE);
		}

		long double Dist(const Atom& a){
				long double dx, dy, dz;
				dx = a.x-x;
				dy = a.y-y;
				dz = a.z-z;
				long double rs = dx*dx+dy*dy+dz*dz;
				long double r = sqrt(rs);
				return(r);
		}
		
		long double ComputeVolume(void){
			long double V = (4.0/3.0)*PI*R*R*R;
			return V;
		}
		long double ComputeSurfaceArea(void){
			long double S = (4.0)*PI*R*R;
			return S;
		}
		bool CheckDist(const Atom& a){
			long double dx, dy, dz;
				dx = a.x-x;
				dy = a.y-y;
				dz = a.z-z;
				long double rs = dx*dx+dy*dy+dz*dz;
				long double r = sqrt(rs);
				if (r>R){
					return 1;
				}
				else{
					return 0;
				}
		}	

			bool CheckDistAll(const Config& c){
				long double dx, dy, dz;
				for (unsigned int i = 0; i < c.natoms; ++i)
				{
					dx = c.atom[i].x-x;
					dy = c.atom[i].y-y;
					dz = c.atom[i].z-z;
					long double rs = dx*dx+dy*dy+dz*dz;
					long double r = sqrt(rs);
					if (r<R){
						return 0;
					}
				
				}
		
				return 1;
				
		}	
         void WrapCoord(long double& boxx, long double& boxy, long double& boxz){
            long double xc = this->x;
			long double yc = this->y;
			long double zc = this->z;

			long double hboxx = boxx/2.0;
			while(xc > hboxx || xc < -hboxx) {	
				if(xc > hboxx){
					xc = xc - boxx;
				}
				else if(xc < -hboxx){
					xc = xc + boxx;
				}
			}	


	
			long double hboxy = boxy/2.0;
			while(yc > hboxy || yc < -hboxy){
				if(yc > hboxy){
					yc = yc - boxy;
				}
				else if(yc < -hboxy){
					yc = yc + boxy;
				}
			}


	        //cout<<"zc : "<<zc<<" boxz: "<<boxz<<endl;
			long double hboxz = boxz/2.0;
            //cout<<"hboxz: "<<hboxz<<endl;
			while( (zc > hboxz) || (zc < -hboxz)){
              //  cout<<"entering wrap loop..."<<endl;				
                if(zc > hboxz){
					zc-= boxz;
				}
				else if(zc < (-hboxz)){
					zc += boxz;
				}
			}
            //cout<<"resetting the coordinates to wrapped values..."<<endl;
            //cout<<"xc "<<xc<<" yc "<<yc<<" zc "<<zc<<endl;
			this->x = xc;
			this->y = yc;
			this->z = zc;
			//cout<<"about to return from wrapping..."<<endl;	
		    return;
         }		
        void CopyCoord(ContinuumSphere& other){
            x=other.x;
            y=other.y;
            z=other.z;
            return;
        }		
};

#endif // __2dContinuumSurface_H_INCLUDED__ 
