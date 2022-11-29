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
		unsigned int maxAtoms;
		unsigned int natoms;
		long double* COM;
		long double boxx, boxy, boxz;
		long double boxold;
		Atom * atom;
		
		
		//Config(): natoms(a_natoms_a), boxx(1.0), boxy(1.0), boxz(1.0) {
 		//	atom = new Atom[a_natoms_a];
		//	COM = new long double[3];
		//}

		Config(unsigned na):maxAtoms(10000), boxx(0.0), boxy(0.0), boxz(0.0){
			natoms = na;
			//boxx = 0.0, boxy = 0.0, boxz = 0.0;
			atom = new Atom[maxAtoms];
			COM = new long double[3];
	
		}

		Config(void):maxAtoms(1000), boxx(0.0), boxy(0.0), boxz(0.0){
			//boxx = 0.0, boxy = 0.0, boxz = 0.0;
		//	atom=NULL;
			natoms=0;
			atom = new Atom[maxAtoms];
			COM = new long double[3];
	
		}

		~Config(){
			
			if (atom){
				
				delete[] atom;
				
			}
			if (COM){
			
				delete[] COM;
			}
			
		}

		//Copy Constructor

		Config(const Config& other){
			natoms = other.natoms;
			atom = new Atom[natoms];
			COM = new long double[3];

			//box size
			boxx = other.boxx;
			boxy = other.boxy;
			boxz = other.boxz;
			
			//copy com
			for (unsigned int i = 0; i < 3; ++i)
			{
				COM[i] = other.COM[i];
			}
			//copy coordinates and properties
			for (unsigned int i = 0; i < natoms; ++i)
			{
				//call atom equate
				atom[i].Equate(other.atom[i]);
			}
			
			return;
		}
				
		void Equate(const Config& other){
			
			for(unsigned int i = 0; i<natoms;++i){
				atom[i].x = other.atom[i].x;
				atom[i].y = other.atom[i].y;
				atom[i].z = other.atom[i].z;
			}
		}
		
		void FullEquate(const Config& other){
			
			//box size
			boxx = other.boxx;
			boxy = other.boxy;
			boxz = other.boxz;
			natoms=other.natoms;
			//copy com
			for (unsigned int i = 0; i < 3; ++i)
			{
				COM[i] = other.COM[i];
			}
			//copy coordinates and properties
			for (unsigned int i = 0; i < natoms; ++i)
			{
				//call atom equate
				atom[i].Equate(other.atom[i]);
			}
			
			return;
		}
		void CopyCoordsUnequalAtoms(const Config& other){
				// test which has larger number of atoms
				unsigned n;
				if (natoms>other.natoms){
					n=other.natoms;
				}
				else{
					n=natoms;
				}
				for (unsigned int i = 0; i < n; ++i)
				{
					atom[i].x = other.atom[i].x;
					atom[i].y = other.atom[i].y;
					atom[i].z = other.atom[i].z;
				}
			return;
		}

		void SwapAtomCoords(Atom& one, Atom& two){
			
			long double tx,ty,tz;
			tx=one.x, ty=one.y, tz=one.z;
			one.x=two.x, one.y=two.y, one.z=two.z;
			two.x=tx, two.y=ty, two.z=tz;
			return;			
		}

		void DeleteAtoms(void){
		
			if (atom){
				delete [] atom;
				atom = NULL;
			}
		
			return;
		}
		void CalcCom(void){
			//cout<<"natoms "<<natoms<<endl;
			long double M = 0.0;
			long double Rx, Ry, Rz;
			long double sumrx=0.0, sumry=0.0, sumrz=0.0;
			
			
			for(unsigned int i = 0;i<natoms;++i){
				
				long double mass = atom[i].mass;
				sumrx += atom[i].x*mass;
				sumry += atom[i].y*mass;
				sumrz += atom[i].z*mass;
				M += mass;
				
			}
			Rx = sumrx/M;
			Ry = sumry/M;
			Rz = sumrz/M;

			COM[0]=Rx;
			COM[1]=Ry;
			COM[2]=Rz;
			return;

		}
	
		long double RadGyr(void){
			CalcCom();
			long double sumgyrx = 0.0, sumgyry = 0.0, sumgyrz = 0.0;
			long double M = 0.0;
			long double Rsq, R;
			long double sumgyr=0.0;
//			cout << " M " << M << endl;
			for(unsigned long y = 0;y<natoms; ++y){
				long double rx, ry, rz;
				rx = atom[y].x;
				ry = atom[y].y;
				rz = atom[y].z;
				long double mass = atom[y].mass;
				sumgyrx = (rx - COM[0])*(rx - COM[0]);
				sumgyry = (ry - COM[1])*(ry - COM[1]);
				sumgyrz = (rz - COM[2])*(rz - COM[2]);
				sumgyr += (sumgyrx + sumgyry + sumgyrz)*mass;
				M+= mass;		
			}

			Rsq = sumgyr;
		
			R = sqrt(Rsq/M);
	
			return(R);	
		}
		
		void ScaleCoordinates(long double& Vscale){

			
			long double Cscale = pow(Vscale, 1.0/3.0);
			for (unsigned i = 0; i < natoms; ++i)
			{
				atom[i].x *= Cscale;
				atom[i].y *= Cscale;
				atom[i].z *= Cscale;
			}


			return;

		}

		long double ComDist(unsigned Aindex){
			long double dx, dy, dz;
			dx = atom[Aindex].x - COM[0];
			dy = atom[Aindex].y - COM[1];
			dz = atom[Aindex].z - COM[2];
			long double r, rsq;
			rsq = (dx*dx) + (dy*dy) + (dz*dz);
			r = sqrt(rsq);

			return(r);
}
		

	bool CheckComDist(unsigned Aindex, long double dist){
			this->CalcCom();
			long double adist = ComDist(Aindex);
			if (adist>dist){
				return 0;
			}
			else{
				return 1;
			}
}
	bool CheckComDistAll(long double dist){
			this->CalcCom();
			bool flag = 1;
			for (unsigned int i = 0; i < natoms; ++i)
			{
				long double r = ComDist(i);
				if (r>dist){
					flag=0;
					break;
				}
			}
			

			return flag;
}

	bool CheckComDistAllInsert(long double dist, long double x, long double y, long double z){
			//compute current com
			this->CalcCom();
			//compute perturbed com by adding the ghost particle
			long double N = (long double)(natoms);
			long double ncx, ncy, ncz;
			ncx= COM[0]*(N);
			ncy= COM[1]*(N);
			ncz= COM[2]*(N);
			ncx+=x;
			ncy+=y;
			ncz+=z;
			ncx/=(N+1.0);
			ncy/=(N+1.0);
			ncz/=(N+1.0);
			
			bool flag = 1;
			//check old particles against perturbed com
			for (unsigned int i = 0; i < natoms; ++i)
			{
				long double dx, dy, dz;
				dx = atom[i].x - ncx;
				dy = atom[i].y - ncy;
				dz = atom[i].z - ncz;
				long double r = sqrt(dx*dx+dy*dy+dz*dz) ;
				if (r>dist){
					flag=0;
					break;
				}
			}
			//check ghost
			long double dx, dy, dz;
				dx = x - ncx;
				dy = y - ncy;
				dz = z - ncz;
				long double r = sqrt(dx*dx+dy*dy+dz*dz) ;
				if (r>dist){
					flag=0;
				}

			return flag;
}

	bool CheckComDistAllDelete(long double dist, unsigned Ai){
			//compute current com
			this->CalcCom();
			//compute perturbed com by subtracting the test delete particle
			long double N = (long double)(natoms);
			long double ncx, ncy, ncz;
			ncx= COM[0]*(N);
			ncy= COM[1]*(N);
			ncz= COM[2]*(N);
			ncx-=atom[Ai].x;
			ncy-=atom[Ai].y;
			ncz-=atom[Ai].z;
			ncx/=(N-1.0);
			ncy/=(N-1.0);
			ncz/=(N-1.0);
			
			bool flag = 1;
			//check old particles (except delete test) against perturbed com
			for (unsigned int i = 0; i < natoms; ++i)
			{
				if (i!=Ai){
					long double dx, dy, dz;
					dx = atom[i].x - ncx;
					dy = atom[i].y - ncy;
					dz = atom[i].z - ncz;
					long double r = sqrt(dx*dx+dy*dy+dz*dz) ;
					if (r>dist){
						flag=0;
						break;
					}
				}
				
			}
			
			return flag;
}
	void SetBox(long double& bx, long double& by, long double& bz){

		boxold=boxx;
		boxx=bx;
		boxy=by;
		boxz=bz;	
		return;
	}
	
	void InitializeAtoms(unsigned nat){
		natoms = nat;
	//	atom = new Atom[nat];
		
	}

	void WrapCoordinates(void){

		for (unsigned int i = 0; i < natoms; ++i)
		{
			long double xc = atom[i].x;
			long double yc = atom[i].y;
			long double zc = atom[i].z;

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


	
			long double hboxz = boxz/2.0;
			while(zc > hboxz || zc < -hboxz){
				if(zc > hboxz){
					zc = zc - boxz;
				}
				else if(zc < -hboxz){
					zc = zc + boxz;
				}
			}

			atom[i].x = xc;
			atom[i].y = yc;
			atom[i].z = zc;
		}
		
		return;

}	

	void WrapAtom(const unsigned& aindex){

			long double xc = atom[aindex].x;
			long double yc = atom[aindex].y;
			long double zc = atom[aindex].z;

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


	
			long double hboxz = boxz/2.0;
			while(zc > hboxz || zc < -hboxz){
				if(zc > hboxz){
					zc = zc - boxz;
				}
				else if(zc < -hboxz){
					zc = zc + boxz;
				}
			}

			atom[aindex].x = xc;
			atom[aindex].y = yc;
			atom[aindex].z = zc;
			
			return;
	}	
		
	void ReCenterZero(void){
		for (unsigned int i = 0; i < natoms; ++i)
		{
			atom[i].x += -COM[0];
			atom[i].y += -COM[1];
			atom[i].z +=-COM[2];
		}
	}
};

#endif // __Config_H_INCLUDED__ 
