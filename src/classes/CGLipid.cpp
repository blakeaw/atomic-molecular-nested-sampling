	
//=================================
// class definition
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/CGLipid.h"

//=================================
//Dependencies
//#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/Bead.h"
//#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/HarmonicBond.h"
//#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/FeneBond.h"

//=================================
// the class functions

void CGLipid::CalcCom(void){
			
			long double M;
			long double Rx, Ry, Rz;
			long double sumrx=0.0, sumry=0.0, sumrz=0.0;
			for(unsigned int i = 0;i<nbeads;++i){
				sumrx += bead[i].x*bead[i].mass;
				sumry += bead[i].y*bead[i].mass;
				sumrz += bead[i].z*bead[i].mass;
				M += bead[i].mass;
			}
			Rx = sumrx/M;
			Ry = sumry/M;
			Rz = sumrz/M;

			COM[0]=Rx;
			COM[1]=Ry;
			COM[2]=Rz;
			return;

		}

long double CGLipid::ComDistFrom(CGLipid& one){
			long double dx, dy, dz;
			long double rs, r;
			this->CalcCom();
			one.CalcCom();
			dx = COM[0]-one.COM[0];
			dy = COM[1]-one.COM[1];
			dz = COM[2]-one.COM[2];
		
			rs = (dx*dx)+(dy*dy)+(dz*dz);
			
			r = sqrt(rs);
		
			return(r);


		}
void CGLipid::EquateLipidCoord(const CGLipid& other){
			unsigned short nb = other.nbeads;
			for (unsigned int i = 0; i < nb; ++i)
			{
				long double otx, oty, otz;
				otx = other.bead[i].x;
				oty = other.bead[i].y;
				otz = other.bead[i].z;
				bead[i].x = otx;
				bead[i].y = oty;
				bead[i].z = otz;
				
			}
			return;
		}
