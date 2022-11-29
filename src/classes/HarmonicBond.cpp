//=================================
// class definition
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/HarmonicBond.h"

//=================================
//Dependencies
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/CGLipid.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/Bead.h"

//=================================
// the class functions
void HarmonicBond::SetMolec(CGLipid * m){
			molec = m;
			return;
		}

void HarmonicBond::SetParams(long double k, long double l){
			constant = k;
			length = l;
			return;
		}

long double HarmonicBond::CalcPotential(void){
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
			long double E = 0.5*constant*(r - length)*(r - length);
			return(E);

		}
long double HarmonicBond::CalcPotential(long double r){
			
			long double E = 0.5*constant*(r - length)*(r - length);
			return(E);

		}
