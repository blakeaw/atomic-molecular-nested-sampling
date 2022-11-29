//=================================
// include guard
#ifndef __Config_H_INCLUDED__
#define __Config_H_INCLUDED__

//=================================
// forward declared dependencies
class Atom;
//=================================
// included dependencies


//=================================
// the actual class
class Config {
	public:
		//int Cindex;
		unsigned int natoms;
		//long double boxx, boxy, boxz;
		Atom * atom;
		
		Config(): natoms(a_natoms_a) {
 			atom = new Atom[a_natoms_a];
		}

		~Config(){
			delete [] atom;
		}
				
		void Equate(const Config& other){
			
			for(unsigned int i = 0; i<natoms;++i){
				atom[i].x = other.atom[i].x;
				atom[i].y = other.atom[i].y;
				atom[i].z = other.atom[i].z;
				atom[i].Eu = other.atom[i].Eu;
				atom[i].dE = other.atom[i].dE;
			}
		}
			
		void UpdateAtomEnergy(void){
				for (unsigned int i = 0; i < natoms; ++i)
				{
					atom[i].Eu+=atom[i].dE;
				}
			
		}
		
		
};

#endif // __Config_H_INCLUDED__ 
