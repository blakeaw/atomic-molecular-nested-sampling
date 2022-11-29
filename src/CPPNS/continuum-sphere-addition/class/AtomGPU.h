//=================================
// include guard
#ifndef __AtomGPU_H_INCLUDED__
#define __AtomGPU_H_INCLUDED__

//=================================
// forward declared dependencies
class Atom;

//=================================
// included dependencies


//=================================
// the actual class
class AtomGPU {
	public:
		//int Aindex;
		double x, y, z;
		double mass;
//		Atom(Atom another) {};
	__device__	AtomGPU(): x(0.0), y(0.0), z(0.0), mass(1.0){};
		//Copy Constructor
	__device__	AtomGPU(const AtomGPU& other){
			x=other.x;
			y=other.y;
			z=other.z;
			mass=other.mass;
		}
	__device__	void Equate(const AtomGPU& other){
			x=other.x;
			y=other.y;
			z=other.z;
			mass=other.mass;
		}

	__host__ void EquateCPU(const Atom& other){
				double omass = (double)(other.mass);
				double ox, oy, oz;
				ox = (double)(other.x);
				oy = (double)(other.y);
				oz = (double)(other.z);
				x=ox;
				y=oy;
				z=oz;
				mass=omass;
//				cudaMemcpy(mass, omass, sizeof(double), cudaMemcpyHostToDevice);
//				cudaMemcpy(x, ox, sizeof(double), cudaMemcpyHostToDevice);
//				cudaMemcpy(y, oy, sizeof(double), cudaMemcpyHostToDevice);
//				cudaMemcpy(z, oz, sizeof(double), cudaMemcpyHostToDevice);
				
			}

	__host__ void UpdateCoord(double nx, double ny, double nz){
//			cudaMemcpy(x, nx, sizeof(double), cudaMemcpyHostToDevice);
//			cudaMemcpy(y, ny, sizeof(double), cudaMemcpyHostToDevice);
//			cudaMemcpy(z, nz, sizeof(double), cudaMemcpyHostToDevice);
			x=nx, y=ny, z=nz;
	}
};

#endif // __Atom_H_INCLUDED__ 
