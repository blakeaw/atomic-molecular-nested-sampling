//=================================
// include guard
#ifndef __MCiterator_H_INCLUDED__
#define __MCiterator_H_INCLUDED__

//=================================
// forward declared dependencies
class Config;
class Dimer;
class dsfmt_t;
long double TotalPairEnergy(const Config& c);
long double DeltaE(const Config& c, const unsigned& Dindex, const unsigned long long& Ar, const long double& ox1, const long double& oy1, const long double& oz1, const long double& ox2, const long double& oy2, const long double& oz2);
void RandomRotation2d(const long double& u, const long double& v, const long double& stepa, Dimer& D, const unsigned long& Ar);
double dsfmt_genrand_close_open(dsfmt_t*);
//=================================
// included dependencies
#ifdef _OPENMP
 #include <omp.h>
#endif
//#include <sstream>
//#include "/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/src/class/AtomD.h"
//#include "/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/src/class/Dimer.h"
//#include "/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/src/class/ConfigD.h"
//#include "/net/uu/nm/cm/bxw109120/src/MersenneTwist/dSFMT-src-2.2.1/dSFMT.c" 
//=================================
// the actual function
void MCiterator(dsfmt_t& dsfmt, Config* config, int medcindex, long double** En, int MCit, int Samples, long double& Emedian){

	
 	#pragma omp parallel for num_threads(2) shared(config, En, medcindex)
	for (unsigned int j = medcindex+1; j <Samples ; ++j)
	{

			#pragma omp flush
//			#pragma omp critical
//			{
//				cout<<"j: "<<j<<endl;
//			}
			//sample index-- upper interval		
			int sind = (int)(En[j][1]);
			//sample index -- lower interval
			int iind = j-medcindex-1;
		
			int bind = (int)(En[iind][1]);
//				#pragma omp critical
//			{
//			cout<<"j "<<j<<" iind "<<iind<<" sind "<<sind<<" bind "<<bind<<endl;
//			}
			// replace upper with lower	
		//	config[sind].Equate(config[bind]);
			
			
			
			//acceptance ratio counters
			int attempttrans=0;
			int accepttrans=0;
			//acceptance ratio counters
			int attemptrot=0;
			int acceptrot=0;
//			long double steplocal = ;
//			long double stepalocal = stepa;
			//copy conf bind into thread local mem
			Config Tlocal(config[bind]);
			long double Ec = TotalPairEnergy(Tlocal);

			long double steplocal = Tlocal.boxx*0.9;
			long double stepalocal = 1.0/12.0;
			//now launch short markov chain to make different
			for (unsigned int m = 0; m < MCit; ++m)
			{
				//Store Initial Energy value
				long double Et1b = Ec;
				//Trial Move -----

				// randomly select an atom index	
				unsigned int Arandom;
				#pragma omp critical
				{
				Arandom = (unsigned int)(dsfmt_genrand_close_open(&dsfmt)*Tlocal.natoms);
				}
				//Get the proper dimer index and atom index within the dimer
				unsigned long Dr = Arandom/2;
				unsigned long long Ar = Arandom - (Dr*2);
				//Store the original coordinates of randomly selected atom
				long double origx, origy, origz;
				origx = Tlocal.dimer[Dr].atom[Ar].x;
				origy = Tlocal.dimer[Dr].atom[Ar].y;
				origz = Tlocal.dimer[Dr].atom[Ar].z;
				//Now store coordinates of other atom in the dimer
				unsigned long long Ar2;
				if (Ar==0)
				{
					Ar2=1;
				}
				else{
					Ar2=0;
				}
				//Store original coordinates
				long double origx2, origy2, origz2;
				origx2 = Tlocal.dimer[Dr].atom[Ar2].x;
				origy2 = Tlocal.dimer[Dr].atom[Ar2].y;
				origz2 = Tlocal.dimer[Dr].atom[Ar2].z;
		
				long double ranprob = dsfmt_genrand_close_open(&dsfmt);
				bool rottag = 0;
				//Only rotate 30% of trial moves
				if(ranprob<=0.30){
					//Pick random x, y, and z values
					long double u, v;
					#pragma omp critical
					{
					u = dsfmt_genrand_close_open(&dsfmt);
					v = dsfmt_genrand_close_open(&dsfmt);
					}
					//rotates about Ar
					RandomRotation2d(u, v, stepalocal, Tlocal.dimer[Dr], Ar);
					long double xn = Tlocal.dimer[Dr].atom[Ar2].x;
//					if (xn!=xn){
//						cout<<"xn is "<<xn<<endl;
//						exit(0);
//					}
					++attemptrot;
					rottag=1;
				}
				else{
					long double movex, movey;
					/* set random movements of the coordinates
				   of the selected atom */
					#pragma omp critical
					{
					movex = (dsfmt_genrand_close_open(&dsfmt)-0.5)*steplocal;
					movey = (dsfmt_genrand_close_open(&dsfmt)-0.5)*steplocal;
					}
					// Implement those moves
					Tlocal.dimer[Dr].atom[Ar].x += movex;
					Tlocal.dimer[Dr].atom[Ar].y += movey;
				
					Tlocal.dimer[Dr].atom[Ar2].x += movex;
					Tlocal.dimer[Dr].atom[Ar2].y += movey;
					++attempttrans;
			
				}	
				//check boundary
				bool dccheck = Tlocal.BoxCheck(Dr);
				//boundary passes
				if (!(dccheck)){
					//Check change in energy for pairwise interactions
					long double dEp=DeltaE(Tlocal, Dr, Ar, origx, origy, origz, origx2, origy2, origz2);
					//Add difference to current total potential energy
					//cout<<"boundary passed. dE is : "<<dEp<<endl;
					Ec+=dEp;
				}
			
				if(dccheck || !(Ec < Emedian)){
					/*Failed criteria--reset atom coordinates and remove
				  the difference in energy */ 
					Tlocal.dimer[Dr].atom[Ar].x = origx;
					Tlocal.dimer[Dr].atom[Ar].y = origy;
				
					Tlocal.dimer[Dr].atom[Ar2].x = origx2;
					Tlocal.dimer[Dr].atom[Ar2].y = origy2;
			
					Ec = Et1b;
				
				}
				else{
					//cout<<"trial move accepted new energy is : "<<Ec<<endl;
					if (rottag){
						++acceptrot;
					}
					else{
						++accepttrans;
					}
				}

				if (m>0 && m%50==0){
						//check acceptance ratios and adjust
					long double art = (long double)(accepttrans)/(long double)(attempttrans);
					long double arr = (long double)(acceptrot)/(long double)(attemptrot);
					//cout<<"art: "<<art<<" arr: "<<arr<<endl;
					//cout<<"stepa: "<<stepa<<endl;
					if (art<0.10){
					//	#pragma omp atomic
						steplocal*=0.5;
						
					}
					else if (art<0.33){
						long double msc = art/0.33;
					//	#pragma omp atomic
						steplocal*=msc;
						
					}
					else if(art>0.51){
						long double msc = art/0.51;
						//#pragma omp atomic
						steplocal*=msc;
					
					}
					accepttrans=0;
					attempttrans=0;
					if (arr<0.10){
						stepalocal*=0.5;
						
					}
					else if (arr<0.33){
						long double msc = arr/0.33;
						stepalocal*=msc;
						
					}
					else if(arr>0.51 && stepalocal<1.0){
						long double msc = arr/0.51;
						stepalocal*=msc;
						if (stepalocal>1.0){
							stepalocal=1.0;
						}
					}
					acceptrot=0;
					attemptrot=0;
					//cout<<"stepa: "<<stepa<<endl;
				}
			} //End Markov Chain
			//Update energy in array
		//	cout<<"j : "<<j<<" Ec: "<<Ec<<endl;
			En[j][0]=Ec;
			//copy Tlocal back to config sind
			config[sind].Equate(Tlocal);
			//#pragma omp flush(config)
			//long double art = (long double)(taccept)/(long double)(MCit);
			//cout<<"j "<<j<<" total acceptance ratio: "<<art<<endl;
			
	} //end upper interval loop
	
};

#endif // __MCiterator_H_INCLUDED__ 
