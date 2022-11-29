//=================================
// include guard
#ifndef __ConfigCGL_H_INCLUDED__
#define __ConfigCGL_H_INCLUDED__

//=================================
// forward declared dependencies
//class CGLipid;
//=================================
// included dependencies
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/CGLipid.h"

//=================================
// the actual class
class ConfigCGL {
	public:
		
		unsigned int nlipids;
		long double boxx, boxy, boxz;
		long double* COM;
	

		CGLipid * lipid;
		//Constructor
		ConfigCGL(): boxx(0.0), boxy(0.0), boxz(0.0) {
 			
			COM = new long double[3];
			for	(short i =0; i<3 ; ++i){
				COM[i] = 0.0;
			}
			
		}

		//Copy Constructor
		ConfigCGL(const ConfigCGL& other){
			long double x, y, z;
			long double bx, by, bz;
			
			unsigned nlip=other.nlipids;
			bx=other.boxx;
			by=other.boxy;
			bz=other.boxz;
			//Initialize same number of lipids in the copy
			nlipids=nlip;
			InitializeLipids(nlip);
			//Copy box sizes
			boxx=bx;
			boxy=by;
			boxz=bz;
			//copy COM values
			COM = new long double[3];
			for	(short i =0; i<3 ; ++i){
				COM[i] = other.COM[i];
			}
			//Copy coodinates
				unsigned N = nlipids*3;
				for(unsigned int i = 0; i<N;++i){
					unsigned lindex = i/3;
					unsigned bindex = i - lindex*3;
					x = other.lipid[lindex].bead[bindex].x;
					y = other.lipid[lindex].bead[bindex].y;
					z = other.lipid[lindex].bead[bindex].z;
					
					
					
					lipid[lindex].bead[bindex].x = x;
					lipid[lindex].bead[bindex].y = y;
					lipid[lindex].bead[bindex].z = z;
					
					
				}
				return;
			
		}
		//Destructor
		~ConfigCGL(){
			delete [] lipid;
			delete [] COM;
		}

		void InitializeLipids(unsigned number);
		void SetBox(long double x, long double y, long double z);
						
		void EquateCoord(const ConfigCGL& other);
		
		void CalcCom(void);

		long double ComDist(CGLipid& one);
		
		bool BoxCheck(void);
		bool BoxCheck(const unsigned& lindex);	

	
};

#endif // __ConfigCGL_H_INCLUDED__ 
