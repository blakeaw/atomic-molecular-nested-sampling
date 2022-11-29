//=================================
// class definition
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/ConfigCGL.h"

//=================================
//Dependencies
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/CGLipid.h"


//=================================
// the class functions




void ConfigCGL::InitializeLipids(unsigned number){
			nlipids = number;
			lipid = new CGLipid[nlipids];
			for (unsigned int i = 0; i < nlipids; ++i)
			{
				//lipid[i].SetConfig(this);
				lipid[i].index = i;
			}
			return;
		}
void ConfigCGL::SetBox(long double x, long double y, long double z){
				boxx=x;
				boxy=y;
				boxz=z;

		}

						
void ConfigCGL::EquateCoord(const ConfigCGL& other){
			long double x, y, z;
			
		
				unsigned N = nlipids*3;
				for(unsigned int i = 0; i<N;++i){
					unsigned lindex = i/3;
					unsigned bindex = i - lindex*3;
					x = other.lipid[lindex].bead[bindex].x;
					y =  other.lipid[lindex].bead[bindex].y;
					z = other.lipid[lindex].bead[bindex].z;
					
					
					
					lipid[lindex].bead[bindex].x = x;
					lipid[lindex].bead[bindex].y = y;
					lipid[lindex].bead[bindex].z = z;
					
					
				}
				return;
			
		}
		
void ConfigCGL::CalcCom(void){
			
			long double M = 0.0;
			long double Rx, Ry, Rz;
			long double sumrx=0.0, sumry=0.0, sumrz=0.0;
			unsigned lindex; 
			unsigned bindex;
			unsigned long long N=nlipids*3;
			for(unsigned int i = 0;i<N;++i){
				lindex = i/3;
				bindex = i - lindex*3;
				long double mass = lipid[lindex].bead[bindex].mass;
				sumrx += lipid[lindex].bead[bindex].x*mass;
				sumry += lipid[lindex].bead[bindex].y*mass;
				sumrz += lipid[lindex].bead[bindex].z*mass;
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

long double ConfigCGL::ComDist(CGLipid& one){
			long double dx, dy, dz;
			long double rs, r;
			one.CalcCom();
			dx = COM[0]-one.COM[0];
			dy = COM[1]-one.COM[1];
			dz = COM[2]-one.COM[2];
		
			rs = (dx*dx)+(dy*dy)+(dz*dz);
			
			r = sqrt(rs);
		
			return(r);


		}
		
bool ConfigCGL::BoxCheck(void){
			bool btag=0;
			unsigned lindex;
			unsigned short bindex;
			unsigned long long N=nlipids*3;
			for(unsigned int i = 0;i<N;++i){
				lindex = i/3;
				bindex = i - (lindex*3);
				long double x = lipid[lindex].bead[bindex].x;
				if(x>boxx/2.0 || x<-boxx/2.0){
					btag = 1;
					
				}
				long double y = lipid[lindex].bead[bindex].y;
				if(y>boxy/2.0 || y<-boxy/2.0){
					btag = 1;
				}
				long double z = lipid[lindex].bead[bindex].z;
				if(z>boxz/2.0 || z<-boxz/2.0){
					btag = 1;
				}
			}
			
			return(btag);
		}

bool ConfigCGL::BoxCheck(const unsigned& lindex){
		bool btag=0;
			
		unsigned long long N=lipid[lindex].nbeads;
		for(unsigned int i = 0;i<N;++i){
				
			long double x = lipid[lindex].bead[i].x;
			if(x>boxx/2.0 || x<-boxx/2.0){
				btag = 1;
				break;
			}
			long double y = lipid[lindex].bead[i].y;
			if(y>boxy/2.0 || y<-boxy/2.0){
				btag = 1;
				break;
			}
			long double z = lipid[lindex].bead[i].z;
			if(z>boxz/2.0 || z<-boxz/2.0){
				btag = 1;
				break;
			}
		}
			
		return(btag);

}
