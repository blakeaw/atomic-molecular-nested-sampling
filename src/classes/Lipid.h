//=================================
// include guard
#ifndef __Lipid_H_INCLUDED__
#define __Lipid_H_INCLUDED__

//=================================
// forward declared dependencies
class Bead;

//=================================
// included dependencies


//=================================
// the actual class
class CGLipid {
	public:
		const unsigned short nbeads;
		const unsigned short nbonds;
		long double* COM;	

		Atom * atom;
		Dimer(): nbeads(3), nbonds(2) {
			bead = new Bead[natoms];
			COM = new long double[3];
			for	(short i =0; i<3 ; ++i){
				COM[i] = 0.0;
			}
		}

		~Dimer(){
			delete [] atom;
			delete [] COM;
		}

		void CalcCom(void){
			
			long double M = atom[0].mass + atom[1].mass;
			long double Rx, Ry, Rz;
			long double sumrx=0.0, sumry=0.0, sumrz=0.0;
			for(unsigned int i = 0;i<natoms;++i){
				sumrx += atom[i].x*atom[i].mass;
				sumry += atom[i].y*atom[i].mass;
				sumrz += atom[i].z*atom[i].mass;
				
			}
			Rx = sumrx/M;
			Ry = sumry/M;
			Rz = sumrz/M;

			COM[0]=Rx;
			COM[1]=Ry;
			COM[2]=Rz;
			return;

		}

		long double ComDistFrom(Dimer& one){
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
};

#endif // __Dimer_H_INCLUDED__ 
