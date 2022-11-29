//=================================
// include guard
#ifndef __ConfigD_H_INCLUDED__
#define __ConfigD_H_INCLUDED__

//=================================
// forward declared dependencies
class Dimer;
//=================================
// included dependencies


//=================================
// the actual class
class Config {
	public:
		
		unsigned int ndimers;
		long double boxx, boxy, boxz;
		long double* COM;
		long double WallZ;
		long double kappa;

		Dimer * dimer;
		
	//	Config(): ndimers(n_dimers_n) {
 	//		dimer = new Dimer[ndimers];
	//		COM = new long double[3];
	//		for	(short i =0; i<3 ; ++i){
	//			COM[i] = 0.0;
	//		}
			
	//	}
		Config(){}


		Config(unsigned nd): ndimers(nd){
			dimer = new Dimer[ndimers];
			COM = new long double[3];
			for	(short i =0; i<3 ; ++i){
				COM[i] = 0.0;
			}
			
		}
		
		~Config(){
			delete [] dimer;
			delete [] COM;
		}
		
		void InitializeDimers(unsigned nd){
			ndimers = nd;
			dimer = new Dimer[nd];
			COM = new long double[3];
			for	(short i =0; i<3 ; ++i){
				COM[i] = 0.0;
			}
			
		}
		void SetBox(long double x, long double y, long double z){
				boxx=x;
				boxy=y;
				boxz=z;

		}

		void SetKappa(long double k){
			kappa = k;
			return;
		}
			
		//Copy Constructor

		Config(const Config& other){

			long double x, y, z, rad, Ae, mass;
			long double sc, rsm;
			bool ctag;
			string type;
			ndimers=other.ndimers;
			dimer = new Dimer[ndimers];
			boxx=other.boxx;
			boxy=other.boxz;
			boxz=other.boxz;
			COM = new long double[3];
			for	(short i =0; i<3 ; ++i){
				COM[i] = other.COM[i];
			}
			WallZ=other.WallZ;
			kappa=other.kappa;
			for(unsigned int k=0; k<ndimers;++k){
				unsigned N = dimer[0].natoms;
				long double BL = other.dimer[k].Lb;
				dimer[k].SetBond(BL);
				for(unsigned int i = 0; i<N;++i){
					x = other.dimer[k].atom[i].x;
					y =  other.dimer[k].atom[i].y;
					z = other.dimer[k].atom[i].z;
					rad = other.dimer[k].atom[i].rad;
					ctag = other.dimer[k].atom[i].ctag;
					mass = other.dimer[k].atom[i].mass;
					type = other.dimer[k].atom[i].type;
					Ae = other.dimer[k].atom[i].Ae;
					sc = other.dimer[k].atom[i].sigc;
					rsm = other.dimer[k].atom[i].rsmall;					

					dimer[k].atom[i].x = x;
					dimer[k].atom[i].y = y;
					dimer[k].atom[i].z = z;
					dimer[k].atom[i].rad = rad;
					dimer[k].atom[i].ctag = ctag ;
					dimer[k].atom[i].type = type;
					dimer[k].atom[i].mass = mass;
					dimer[k].atom[i].Ae = Ae;
					dimer[k].atom[i].sigc = sc;
					dimer[k].atom[i].rsmall = rsm;
				}
			}
		}


		

	
		void Equate(const Config& other){
			long double x, y, z, rad, Ae, mass;
			long double sc, rsm;
			bool ctag;
			string type;
			unsigned int nd=ndimers;
		
			for(unsigned int k=0; k<nd;++k){
				unsigned N = dimer[0].natoms;
				dimer[k].Lb = other.dimer[k].Lb;
				for(unsigned int i = 0; i<N;++i){
					x = other.dimer[k].atom[i].x;
					y =  other.dimer[k].atom[i].y;
					z = other.dimer[k].atom[i].z;
					rad = other.dimer[k].atom[i].rad;
					ctag = other.dimer[k].atom[i].ctag;
					mass = other.dimer[k].atom[i].mass;
					type = other.dimer[k].atom[i].type;
					Ae = other.dimer[k].atom[i].Ae;
					sc = other.dimer[k].atom[i].sigc;
					rsm = other.dimer[k].atom[i].rsmall;


					dimer[k].atom[i].x = x;
					dimer[k].atom[i].y = y;
					dimer[k].atom[i].z = z;
					dimer[k].atom[i].rad = rad;
					dimer[k].atom[i].ctag = ctag ;
					dimer[k].atom[i].type = type;
					dimer[k].atom[i].mass = mass;
					dimer[k].atom[i].Ae = Ae;
					dimer[k].atom[i].sigc = sc;
					dimer[k].atom[i].rsmall = rsm;
				}
			}
		}
			
		void CalcCom(void){
			
			long double M = 0.0;
			long double Rx, Ry, Rz;
			long double sumrx=0.0, sumry=0.0, sumrz=0.0;
			unsigned Di;
			unsigned short Ai;
			unsigned long long N=ndimers*2;
			for(unsigned int i = 0;i<N;++i){
				Di = i/2;
				Ai = i - (Di*2);
				long double mass = dimer[Di].atom[Ai].mass;
				sumrx += dimer[Di].atom[Ai].x*mass;
				sumry += dimer[Di].atom[Ai].y*mass;
				sumrz += dimer[Di].atom[Ai].z*mass;
				M += mass;
//				if (sumrx!=sumrx){
//					cout<<" x : "<<dimer[Di].atom[Ai].x<<" mass: "<<mass<<endl;
//				}
				
			}
			Rx = sumrx/M;
			Ry = sumry/M;
			Rz = sumrz/M;

			COM[0]=Rx;
			COM[1]=Ry;
			COM[2]=Rz;
			//cout<<"COMx "<<COM[0]<<" COMy: "<<COM[0]<<" M "<<M<<endl;
			return;

		}

		long double ComDist(Dimer& one){
			long double dx, dy, dz;
			long double rs, r;
			one.CalcCom();
			dx = COM[0]-one.COM[0];
			dy = COM[1]-one.COM[1];
			dz = COM[2]-one.COM[2];
		//	cout << "dx is " << dx << " dy is " << dy << " dz is " << dz << endl;
			rs = (dx*dx)+(dy*dy)+(dz*dz);
			//cout << "distsquared is " << rs << endl;
			r = sqrt(rs);
		//	cout << "dist is " << r << endl;
			return(r);


		}
		
		bool BoxCheck(void){
			bool btag=0;
			unsigned Di;
			unsigned short Ai;
			unsigned long long N=ndimers*2;
			for(unsigned int i = 0;i<N;++i){
				Di = i/2;
				Ai = i - (Di*2);
				long double x = dimer[Di].atom[Ai].x;
				if(x>boxx/2.0 || x<-boxx/2.0){
					btag = 1;
					break;
				}
				long double y = dimer[Di].atom[Ai].y;
				if(y>boxy/2.0 || y<-boxy/2.0){
					btag = 1;
					break;
				}
				long double z = dimer[Di].atom[Ai].z;
				if(z>boxz/2.0 || z<-boxz/2.0){
					btag = 1;
					break;
				}
			}
			
			return(btag);
		}

		bool BoxCheck(const unsigned& dindex){
			bool btag=0;
			
			unsigned long long N=2;
			for(unsigned int i = 0;i<N;++i){
				
				long double x = dimer[dindex].atom[i].x;
				if(x>boxx/2.0 || x<-boxx/2.0){
					btag = 1;
					break;
				}
				long double y = dimer[dindex].atom[i].y;
				if(y>boxy/2.0 || y<-boxy/2.0){
					btag = 1;
					break;
				}
				long double z = dimer[dindex].atom[i].z;
				if(z>boxz/2.0 || z<-boxz/2.0){
					btag = 1;
					break;
				}
			}
			
			return(btag);

		}
};

#endif // __ConfigD_H_INCLUDED__ 
