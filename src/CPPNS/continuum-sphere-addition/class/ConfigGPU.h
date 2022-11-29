//=================================
// include guard
#ifndef __ConfigGPU_H_INCLUDED__
#define __ConfigGPU_H_INCLUDED__

//=================================
// forward declared dependencies
//class AtomGPU;
//class Config;
//=================================
// included dependencies


//=================================
// the actual class
class ConfigGPU {
	public:
		
		//values
		int natoms;
		double* COM;
		double *box;
		double *ax, *ay, *az;
	

	//Constructor
	__device__ ConfigGPU(int na) {
			natoms= na;
			// Allocate host
			box = new double[3];
			COM = new double[3];
			ax = new double[na];
			ay = new double[na];
			az = new double[na];
				
		}

	//Copy Constructor

//	__device__ ConfigGPU(const Config& other){
//			natoms = other.natoms;
//			atom = new Atom[natoms];
//			COM = new long double[3];

//			//box size
//			boxx = other.boxx;
//			boxy = other.boxy;
//			boxz = other.boxz;
//			
//			//copy com
//			for (unsigned int i = 0; i < 3; ++i)
//			{
//				COM[i] = other.COM[i];
//			}
//			//copy coordinates and properties
//			for (unsigned int i = 0; i < natoms; ++i)
//			{
//				//call atom equate
//				atom[i].Equate(other.atom[i]);
//			}
//			
//			return;

//	}
		//destructor
	__device__	~ConfigGPU(){
			delete [] box;
			delete [] COM;
			delete [] ax;
			delete [] ay;
			delete [] az;
			
		}

	

		

};

#endif // __Config_H_INCLUDED__ 
