//=================================
// class definition
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/FeneBond.h"

//=================================
//Dependencies
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/CGLipid.h"


//=================================
// the class functions



void FeneBond::SetMolec(CGLipid * m){
			molec = m;
			return;
		}

void FeneBond::SetParams(long double e, long double k, long double dl, long double bb){
			epsilon = e;
			constant = k;
			divlength = dl;
			b = bb;
			return;
		}

long double FeneBond::CalcPotential(void){
			long double x1 = molec->bead[Bindex1].x;
			long double y1 = molec->bead[Bindex1].y;
			long double z1 = molec->bead[Bindex1].z;
			long double x2 = molec->bead[Bindex2].x;
			long double y2 = molec->bead[Bindex2].y;
			long double z2 = molec->bead[Bindex2].z;
			
			long double dx, dy, dz;
			dx = x2 - x1;
			dy = y2 - y1;
			dz = z2 - z1;
			long double r = sqrt(dx*dx+dy*dy+dz*dz);
			long double rc = pow(2.0, 1.0/6.0)*b;
			long double E1 = 0.5*constant*divlength*divlength*log(1.0 - (r/divlength)*(r/divlength));
			long double E2=0.0;
			if (r<=rc){
				long double r2 = r*r;
				long double b2 = b*b;
				E2 = 4.0*epsilon*(pow(b2/r2, 6)-pow(b2/r2, 3) + 0.25);
			} 
			//cout << "Inside Fene Calc potential; r: " << r << " E1: " << E1 << " E2: " << E2 << endl;
		//	cout << "Bindex1: " << Bindex1 << " Bindex2: " << Bindex2 << endl;
			long double E = E1+E2;
			return(E);

		}
		
long double FeneBond::CalcPotential(long double r){
			
			long double rc = pow(2.0, 1.0/6.0)*b;
			long double E1 = 0.5*constant*divlength*divlength*log(1.0 - (r/divlength)*(r/divlength));
			long double E2;
			if (r<=rc){
				const long double b2 = b*b;
				const long double r2 = r*r;
				E2 = 4.0*epsilon*(pow(b2/r2, 6)-pow(b2/r2, 3) + 0.25);
			} 
			else{
				E2 = 0.0;
			}
			long double E = E1+E2;
			return(E);

		}
