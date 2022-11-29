/* Nested Sampling Routine--Collected nested energies.
Configured to run for lj particles with sigma=1.0 and epsilon=1.0.
The Monte Carlo trial move consist of moving a single atom 
 The boundary condition is set to a cubic box with min image periodicity; 
---Contains conditions to modulate step sizes under each new median energy.

Modifications


9-16-14 Added output of a potential energy histogram from each NS step

Author:Blake Wilson  blake.wilson@utdallas.edu
-------------------------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string>
#include <fstream>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <time.h>
#include <sstream>
// OpenMP Library
#ifdef _OPENMP
 #include <omp.h>
#endif
//RNG Fast Mersenne Twist
#include "/net/uu/nm/cm/bxw109120/src/MersenneTwist/dSFMT-src-2.2.1/dSFMT.c" 



using namespace std;

 // Classes
#include "/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/src/class/AtomD.h"
#include "/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/src/class/Dimer.h"
#include "/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/src/class/ConfigD.h"
#include "/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/src/class/HistogramNS_slim.h"
#include "/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/src/class/RunningStats.h"
#include "/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/src/class/Histogram.h"

//#include "./src/class/CenterofMass.h"


// prototype functions

long double TotalPairEnergy(const Config& c);
//long double DeltaE(const Config& c, const unsigned& Dindex, const unsigned long long& Ar, long double ox1, long double oy1, long double oz1, long double ox2, long double oy2, long double oz2);
long double DeltaE(const Config& c, const unsigned& Dindex, const unsigned long long& Ar, const long double& ox1, const long double& oy1, const long double& oz1, const long double& ox2, const long double& oy2, const long double& oz2);
long double PairSChargeShell(const Atom& one, const Atom& two, const long double& r, const long double& kappa);
long double PairLJColloid(const Atom& one, const Atom& two, const long double& r);
long double RDist(const Atom& one, const Atom& two);
long double RDisto(const Atom& one, const long double& xo, const long double& yo, const long double& zo);
long double RadGyr(Config& c);
long double LjWallZ(const Atom& one, const long double& WallZ );
long double LjWallZp(const long double& oz, const long double& rad, const long double& WallZ);
long double TotalWallEnergyZ(const Config& c);
long double dEWallZ(const Dimer& d, long double oz1, long double r1, long double oz2, long double r2, const long double& WallZ);
long double AverageZorientation(const Config& c);
void Zdensity(const Config& c, Histogram& gofz);
long double NPvecProd(const Config& c);
void RandomRotation2d(const long double& u, const long double& v, const long double& stepa, Dimer& D, const unsigned long& Ar);
long double AverageZwallDist(Config& c);
string itos( int Number );
void CountContacts(long double& Npp, long double& Npn, long double& Nnn, Config& c, long double& packparm);
inline long double fastPow(long double a, long double b);
long double NPvecProdNeighbor(const Config& c);
bool CompareDoubles2 (long double A, long double B, long double tolerance);
void dquickSort(long double arr[][2], int left, int right);

int main(int argc, char *argv[])
{

	//Get the local time for the logfile
	time_t rawtime;
 	struct tm * timeinfo;
 	time (&rawtime);
 	timeinfo = localtime (&rawtime);

	// Parameters that need to be set

	// the number of atoms to have in the system
	const unsigned short ndimers = 100;
	const unsigned natoms = ndimers*2; 
	
	// the number of nested sampling iterations to run 
	const unsigned int Iterations = 20000;
	//Number of samples to use--should be odd
	const unsigned int Samples = 201;
	//Total number of MC style trial moves under each NS iteration
	const unsigned int MCit = natoms*10;
	//const unsigned short replicatenum = n_rep_n ;
	const unsigned short replicatenum = 0 ;

	
	//Dimer Properties
	const long double rad1 = 1.0; // Positive charged
	const bool ctag1 = 1;
	long double Ae1 = 1.0;
	const long double Ac1 = 1.0;
	const long double mass1 = 1.0;
	string type1 = "P";
	const long double rad2 = 1.0; // Negative charged
	const bool ctag2 = 0;
	long double Ae2 = 1.0;
	const long double Ac2 = 1.0;
	const long double mass2 = 1.0;
	string type2 = "N";	
	const long double BondLength = rad1+rad2;
	//long double sigconv = 50.0;
	// kappa value - inverse screening length	
	long double rmin = pow(2.0, 1.0/50.0)*(rad1+rad2);
	//long double minscale = k_val_k;
	long double minscale = 0.5;
	long double kappa;
	// determin kappa based on electro decay term scaling at the colloid potential
	// minimum i.e. e^(-kappa*rmin)
	if (minscale>0.0){
		kappa = - log(minscale)/rmin;
	}
	else{
		//turn off electrostat - infinite screening
		kappa=1.0;
		Ae1=0.0;
		Ae2=0.0;
	}
	//const long double kappa = 1.0 ;
	const long double density = 0.040;
	// Hard Cubic box size
//	const long double box = ( (long double)(ndimers) )*(rad1+rad2)/2.0 ;
	const long double box = sqrt( (1.0/density)*( (long double)(ndimers) ) );
		//Translation step
	const long double bsfac = 0.900;
	long double step = bsfac*box;
	//Rotational step;
	long double stepa = 1.0/12.0;
	//Pi
	const long double pi = 3.14159265359;
	
	//Enthalpy bounds
	const long double Etop = 200.0;
	const long double Elow = (-100.0)*(long double)(ndimers);


		//run descriptor
	string nds = itos(ndimers);
	string rundescriptor = "cube_d"+nds+"_test_xxntestnxx_3";
	
	
//------------------------------------------------------------

	string rundescript = rundescriptor;

	
	string t_id = itos(replicatenum);
	rundescript+="_rep"+t_id;

	cout<<"Will run replica " <<replicatenum<<" of system:"<<endl;
	cout <<"---- " <<ndimers<<" dimers and a total of " << natoms << " atoms" << endl;
	cout<<"density set to: "<<density<<" and box is : "<<box<<endl;
		//Name of logfile
	string logfilename = "oMMC_" + rundescript + ".log";
	//Name of Emedian output file
	string eoutname = "Epoints_int_" + rundescript + ".dat";
	//Name of Configuration output file
	string coutname = "ConfigsDimer_int_" + rundescript + ".xyz";
	//Name of Observable output file
	//string outfilename = "zOutput_int_" + rundescript + ".dat";
	//Name of NS values file for observables
	string nsoboutname = "zNSobs_int_" + rundescript + ".dat";
	
	ofstream logfile(logfilename.c_str());
	ofstream eoutfile(eoutname.c_str());
	
	//ofstream outfile(outfilename.c_str());
	ofstream coordsfile(coutname.c_str());
	ofstream nsobout(nsoboutname.c_str());
	eoutfile << setiosflags(ios::fixed) << setprecision(15);
	nsobout << setiosflags(ios::fixed) << setprecision(15);
cout << "Nested sampling iterations is set to " << Iterations << " with MC iterations set to " << MCit << endl;
cout<<"Etop: "<<Etop<<" Elow: "<<Elow<<endl;	

logfile << " " << endl;
	logfile << " " << endl;
	logfile << "Local time and date: " << asctime(timeinfo) << endl;
	logfile << "Parameters: " << endl;
	logfile << "Nested Sampling Energy Integration Fraction: " << 0.5 << endl;
	logfile << "R1: " << rad1 << " R2: " << rad2 << " ctag1: " << ctag1 << " ctag2 " << ctag2 << endl;
	logfile << "Ae1: " << Ae1 << " Ae2: " << Ae2 << " mass1: " << mass1 << " mass2: " << mass2 << endl;
	logfile<<"Ac1: "<<Ac1<<" Ac2: "<<Ac2<<endl;
	logfile << "kappa: " << kappa << endl;
	//logfile<<"MCitpp: "<<MCitpp<<"ConfigsOutFreq: "<<ConfigsOutFreq<<endl;

	// Initialize the MT RNG object
	dsfmt_t dsfmt;
	//seed -using time
	unsigned long long seed = time(0);
	
	logfile<<"RNG seed: "<<seed<<endl;
	logfile<<endl;
	// Initialize the RNG object with the seed
	
	dsfmt_init_gen_rand(&dsfmt, seed);


//Initialize configuration objects
	

//if (Samples%2==0){
//	++Samples;
//	cout<<"Samples was even and has been reset to "<<Samples<<endl;
//}

//Initialize the configuration objects
Config * config  = new Config[Samples];
// Initialize an array to store the Energy values
long double En[Samples][2];
for (unsigned int i = 0; i < Samples; ++i)
{
	config[i].InitializeDimers(ndimers);
	config[i].SetBox(box, box, box);
	config[i].kappa=kappa;
}

cout<<"Assigning Initial coordinates for "<<Samples<<" configurations.."<<endl;
logfile<<"Assigning Initial coordinates.."<<endl;
//Randomly place atoms in box - but filter to configuration with Enthalpy below Htop
for (unsigned int k = 0; k < Samples; ++k)
{
	
	logfile << "Dimer Bond Length: " << BondLength << endl;
	logfile<<endl;
	//Assign dimer parameters
	logfile<<"Assigning particle parameters..."<<endl;
	for(unsigned int jj = 0;jj<ndimers;jj++){
			
			config[k].dimer[jj].SetBond(BondLength);
			config[k].dimer[jj].atom[0].type = type1;
			config[k].dimer[jj].atom[0].rad = rad1;
			config[k].dimer[jj].atom[0].mass = mass1;	
			config[k].dimer[jj].atom[0].ctag = ctag1;	
			config[k].dimer[jj].atom[0].Ae = Ae1;	
			config[k].dimer[jj].atom[0].Ac = Ac1;			
			config[k].dimer[jj].atom[1].type = type2;
			config[k].dimer[jj].atom[1].rad = rad2;
			config[k].dimer[jj].atom[1].mass = mass2;
			config[k].dimer[jj].atom[1].ctag = ctag2;	
			config[k].dimer[jj].atom[1].Ae = Ae2;	
			config[k].dimer[jj].atom[1].Ac = Ac2;
	}
	//Assign an Initial Random Configuration 
	logfile<<"Assigning an initial random configuration in the box..."<<endl;
	bool atag = 0;
	long double Et=0.0;
	while(atag == 0){
		atag=1;
		//Assign coordinates	
		for (unsigned int jj = 0; jj < ndimers; ++jj)
		{
			long double ranx, rany;
			ranx = (dsfmt_genrand_close_open(&dsfmt)-0.5);
			rany = (dsfmt_genrand_close_open(&dsfmt)-0.5);
			//ranz = (dsfmt_genrand_close_open(&dsfmt)-0.5);
			long double x, y;
			x = ranx*(box);
			y = rany*(box);
			//z = ranz*(boxzz);
			config[k].dimer[jj].atom[0].x=x;
			config[k].dimer[jj].atom[0].y=y;
			config[k].dimer[jj].atom[0].z= 0.0;
			long double u;
			u = dsfmt_genrand_close_open(&dsfmt);
			//v = dsfmt_genrand_close_open(&dsfmt);		
			long double theta = u*2.0*pi;
			//long double phi = acos((2.0*v) - 1.0);
			long double BondLength = config[k].dimer[jj].Lb;	
			//long double dz = BondLength*cos(phi);
			//long double r = sqrt(BondLength*BondLength-dz*dz);
			long double r = BondLength;
			long double dx = r*cos(theta);
			long double dy = r*sin(theta);
			config[k].dimer[jj].atom[1].x = x+dx;
			config[k].dimer[jj].atom[1].y = y+dy;
			config[k].dimer[jj].atom[1].z = 0.0;
			bool flag = 0;
			for (unsigned int i = 0; i < jj; ++i)
			{
				long double rnn, rnp, rpp, rpn;
				rpp = RDist(config[k].dimer[jj].atom[0], config[k].dimer[i].atom[0]);
				rpn = RDist(config[k].dimer[jj].atom[0], config[k].dimer[i].atom[1]); 
				rnp = RDist(config[k].dimer[jj].atom[1], config[k].dimer[i].atom[0]);
				rnn = RDist(config[k].dimer[jj].atom[1], config[k].dimer[i].atom[1]);
				long double tdist = r;
				if (rpp<tdist || rpn<tdist || rnp< tdist || rnn<tdist){
					flag=1;
					break;
				}
			}
			if (flag){
				--jj;
			}
		} // end assigment for
		
		config[k].CalcCom();
		//Test boundary condition
		if(config[k].BoxCheck() == 1){
			atag = 0;
			
		}
		else{
			TotalPairEnergy(config[k]);
			
			if(Et>Etop-0.0010){
				atag=0;
			}	

		}
		//long double Et = TotalPairEnergy(config);//+TotalWallEnergyZ(config);
//		if (Et<1000.00){
//			cout<<"Et "<<Et<<endl;
//		}
		
			
	} // end assignment whil
	cout<<"Et "<<Et<<" sample : "<<k<<endl;
	En[k][0]=Et;
	En[k][1]=k;
}

cout<<"Initial coordinates assigned!"<<endl;
logfile<<"Initial coordinates assigned!"<<endl;

//Track the radius of gyration and gyration density
	RunningStats Rgstat;
	RunningStats Dstat;
	// Average Energy and Average Energy^2
	RunningStats Estat;
	RunningStats Esqstat;
	//Orientation OPs
	RunningStats Dorstat;

	//hcp packing order param;
	RunningStats Pparmstat;
	//N-N contacts
	RunningStats PPstat;
	//N-P contacts
	RunningStats NPstat;
	//P-P contacts
	RunningStats NNstat; 

//Sort Energies and get median

dquickSort(En, 0, Samples-1);	
//define the median index
int medcindex= (Samples/2);

//Define the median energy
long double Emedian = En[medcindex][0];

// to track convergence
long double Emprev=Emedian;
// out put initial values to files
eoutfile<<Emedian<<endl;
int medcind = (int)(En[medcindex][1]);
//observables to output

for (unsigned int eo = 0; eo < Samples; ++eo)
{
		//Calculate Observables				
	long double Rg = RadGyr(config[eo]);
	long double dynVol = (4.0/3.0)*pi*(Rg*Rg*Rg);
	long double dynDens = ((long double)(config[eo].ndimers)*2.0)/dynVol;
	long double Dor= NPvecProdNeighbor(config[eo]);
	long double Packparm;
	long double Npp, Npn, Nnn;
	Npp=0.0, Npn=0.0, Nnn=0.0;
	CountContacts(Npp, Npn, Nnn, config[eo], Packparm);




	//Push observables for Trajectory Averages
	long double Eetloc=En[eo][0];
	Estat.Push(Eetloc);
	Esqstat.Push(Eetloc*Eetloc);
	Rgstat.Push(Rg);
	Dstat.Push(dynDens);
	Dorstat.Push(Dor);
	Pparmstat.Push(Packparm);
	PPstat.Push(Npp);
	NPstat.Push(Npn);
	NNstat.Push(Nnn);
		
}
//Trajectory averages
long double Estatmean=Estat.Mean();
long double Esqstatmean=Esqstat.Mean();
long double Rgstatmean=Rgstat.Mean();
long double Dstatmean=Dstat.Mean();
long double Dorstatmean=Dorstat.Mean();
long double Pparmstatmean=Pparmstat.Mean();
long double PPstatmean=PPstat.Mean();
long double NPstatmean=NPstat.Mean();
long double NNstatmean=NNstat.Mean();
//output to file
nsobout << 0 << " " << Emedian << " " << Estatmean << " " << Esqstatmean << " ";
nsobout << Rgstatmean << " " << Dstatmean << " " << " " << Dorstatmean << " " << Pparmstatmean<< " ";
nsobout<<PPstatmean<<" "<<NPstatmean<<" "<<NNstatmean<<endl;
// Observables-- logfile output - same as screen output
		logfile << "* Trajectory Average Values:" << endl;
		logfile << "Energy: " << Estatmean << " and Energy Squared: " << Esqstatmean << endl;
		logfile << "Radius of gyration: " << Rgstatmean << endl;
		logfile << "Gyration Density: " << Dstatmean << endl;
		logfile << "Dimer-Dimer(NP-NP) Orientation Parameter " << Dorstatmean << endl;
		logfile << "HCP packing parameter: " << Pparmstatmean << endl;
		logfile << "Contacts; PP: "<<PPstatmean<<" NP: "<<NPstatmean<<" NN: "<<NNstatmean<<endl;
		// Observables-- logfile output - same as screen output
		cout << "* Trajectory Average Values:" << endl;
		cout << "Energy: " << Estatmean << " and Energy Squared: " << Esqstatmean << endl;
		cout << "Radius of gyration: " << Rgstatmean << endl;
		cout << "Gyration Density: " << Dstatmean << endl;
		cout << "Dimer-Dimer(NP-NP) Orientation Parameter " << Dorstatmean << endl;
		cout << "HCP packing parameter: " << Pparmstatmean << endl;
		cout << "Contacts; PP: "<<PPstatmean<<" NP: "<<NPstatmean<<" NN: "<<NNstatmean<<endl;
//reset running stat objects
Rgstat.Reset();
		Dstat.Reset();
		Dorstat.Reset();
		Estat.Reset();
		Esqstat.Reset();
		Pparmstat.Reset();
		PPstat.Reset();
		NPstat.Reset();
		NNstat.Reset();
//start the NS loop
//output median config
coordsfile<<natoms<<endl;
coordsfile<<"charged shell lj particles with Energy "<<Emedian<<endl;
for(unsigned long long hg = 0; hg<natoms; hg++){
	unsigned long Di = hg/2;
	unsigned long long Ai = hg - (2*Di);
	
	coordsfile<<setiosflags(ios::fixed)<<setprecision(15);
	string ttype = config[medcind].dimer[Di].atom[Ai].type;
	long double tx, ty, tz;
	tx = config[medcind].dimer[Di].atom[Ai].x;
	ty = config[medcind].dimer[Di].atom[Ai].y;
	tz = config[medcind].dimer[Di].atom[Ai].z;
	coordsfile<<ttype.c_str()<<" "<<tx<<" "<<ty<<" "<<tz<< endl;
}

for(unsigned int ii=1;ii<Iterations;++ii){
	cout << "Iteration " << ii << "......." << endl;
	logfile << "Iteration " << ii << "......." << endl;
	cout<<"Emedian is : "<<Emedian<<endl;
	logfile<<"Emedian is : "<<Emedian<<endl;
	medcind = (int)(En[medcindex][1]);


//first loop over samples with E>Em and replace with ones E<Em
	#pragma omp for
	for (unsigned int j = medcindex+1; j <Samples ; ++j)
	{
			//sample index-- upper interval		
			int sind = (int)(En[j][1]);
			//sample index -- lower interval
			int bind = (int)(En[j-medcindex-1][1]);
			// replace upper with lower	
			config[sind].Equate(config[bind]);
			
			long double Ec = TotalPairEnergy(config[sind]);
			//acceptance ratio counters
			int attempttrans=0;
			int accepttrans=0;
			//acceptance ratio counters
			int attemptrot=0;
			int acceptrot=0;
		
			//now launch short markov chain to make different
			for (unsigned int m = 0; m < MCit; ++m)
			{
				//Store Initial Energy value
				long double Et1b = Ec;
				//Trial Move -----

				// randomly select an atom index	
				unsigned int Arandom = (unsigned int)(dsfmt_genrand_close_open(&dsfmt)*natoms);
				//Get the proper dimer index and atom index within the dimer
				unsigned long Dr = Arandom/2;
				unsigned long long Ar = Arandom - (Dr*2);
				//Store the original coordinates of randomly selected atom
				long double origx, origy, origz;
				origx = config[sind].dimer[Dr].atom[Ar].x;
				origy = config[sind].dimer[Dr].atom[Ar].y;
				origz = config[sind].dimer[Dr].atom[Ar].z;
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
				origx2 = config[sind].dimer[Dr].atom[Ar2].x;
				origy2 = config[sind].dimer[Dr].atom[Ar2].y;
				origz2 = config[sind].dimer[Dr].atom[Ar2].z;
		
				long double ranprob = dsfmt_genrand_close_open(&dsfmt);
				bool rottag = 0;
				//Only rotate 30% of trial moves
				if(ranprob<=0.30){
					//Pick random x, y, and z values
					long double u, v;
					u = dsfmt_genrand_close_open(&dsfmt);
					v = dsfmt_genrand_close_open(&dsfmt);
			
					//rotates about Ar
					RandomRotation2d(u, v, stepa, config[sind].dimer[Dr], Ar);
					long double xn = config[sind].dimer[Dr].atom[Ar2].x;
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
					movex = (dsfmt_genrand_close_open(&dsfmt)-0.5)*step;
					movey = (dsfmt_genrand_close_open(&dsfmt)-0.5)*step;
			
					// Implement those moves
					config[sind].dimer[Dr].atom[Ar].x += movex;
					config[sind].dimer[Dr].atom[Ar].y += movey;
				
					config[sind].dimer[Dr].atom[Ar2].x += movex;
					config[sind].dimer[Dr].atom[Ar2].y += movey;
					++attempttrans;
			
				}	
				//check boundary
				bool dccheck = config[sind].BoxCheck(Dr);
				//boundary passes
				if (!(dccheck)){
					//Check change in energy for pairwise interactions
					long double dEp=DeltaE(config[sind], Dr, Ar, origx, origy, origz, origx2, origy2, origz2);
					//Add difference to current total potential energy
					//cout<<"boundary passed. dE is : "<<dEp<<endl;
					Ec+=dEp;
				}
			
				if(dccheck || !(Ec < Emedian)){
					/*Failed criteria--reset atom coordinates and remove
				  the difference in energy */ 
					config[sind].dimer[Dr].atom[Ar].x = origx;
					config[sind].dimer[Dr].atom[Ar].y = origy;
				
					config[sind].dimer[Dr].atom[Ar2].x = origx2;
					config[sind].dimer[Dr].atom[Ar2].y = origy2;
			
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
						step*=0.5;
						
					}
					else if (art<0.33){
						long double msc = art/0.33;
						step*=msc;
						
					}
					else if(art>0.51){
						long double msc = art/0.51;
						step*=msc;
					
					}
					accepttrans=0;
					attempttrans=0;
					if (arr<0.10){
						stepa*=0.5;
						
					}
					else if (arr<0.33){
						long double msc = arr/0.33;
						stepa*=msc;
						
					}
					else if(arr>0.51 && stepa<1.0){
						long double msc = arr/0.51;
						stepa*=msc;
						if (stepa>1.0){
							stepa=1.0;
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
			//long double art = (long double)(taccept)/(long double)(MCit);
			//cout<<"j "<<j<<" total acceptance ratio: "<<art<<endl;

	} //end upper interval loop

	//Now put the median energy conf on a Markov Chain enforcing E<Emed
	long double Ec = Emedian;
	//acceptance ratio counters
	//acceptance ratio counters
			int attempttrans=0;
			int accepttrans=0;
			//acceptance ratio counters
			int attemptrot=0;
			int acceptrot=0;
	for (unsigned int m = 0; m < MCit; ++m)
	{
		//Store Initial Energy value
				long double Et1b = Ec;
				//Trial Move -----

				// randomly select an atom index	
				unsigned int Arandom = (unsigned int)(dsfmt_genrand_close_open(&dsfmt)*natoms);
				//Get the proper dimer index and atom index within the dimer
				unsigned long Dr = Arandom/2;
				unsigned long long Ar = Arandom - (Dr*2);
				//Store the original coordinates of randomly selected atom
				long double origx, origy, origz;
				origx = config[medcind].dimer[Dr].atom[Ar].x;
				origy = config[medcind].dimer[Dr].atom[Ar].y;
				origz = config[medcind].dimer[Dr].atom[Ar].z;
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
				origx2 = config[medcind].dimer[Dr].atom[Ar2].x;
				origy2 = config[medcind].dimer[Dr].atom[Ar2].y;
				origz2 = config[medcind].dimer[Dr].atom[Ar2].z;
		
				long double ranprob = dsfmt_genrand_close_open(&dsfmt);
				bool rottag = 0;
				//Only rotate 30% of trial moves
				if(ranprob<=0.30){
					//Pick random x, y, and z values
					long double u, v;
					u = dsfmt_genrand_close_open(&dsfmt);
					v = dsfmt_genrand_close_open(&dsfmt);
			
					//rotates about Ar
					RandomRotation2d(u, v, stepa, config[medcind].dimer[Dr], Ar);
					++attemptrot;
					rottag=1;
				}
				else{
					long double movex, movey;
					/* set random movements of the coordinates
				   of the selected atom */
					movex = (dsfmt_genrand_close_open(&dsfmt)-0.5)*step;
					movey = (dsfmt_genrand_close_open(&dsfmt)-0.5)*step;
			
					// Implement those moves
					config[medcind].dimer[Dr].atom[Ar].x += movex;
					config[medcind].dimer[Dr].atom[Ar].y += movey;
				
					config[medcind].dimer[Dr].atom[Ar2].x += movex;
					config[medcind].dimer[Dr].atom[Ar2].y += movey;
					++attempttrans;
			
				}	
				//check boundary
				bool dccheck = config[medcind].BoxCheck(Dr);
				//boundary passes
				if (!(dccheck)){
					//Check change in energy for pairwise interactions
					long double dEp=DeltaE(config[medcind], Dr, Ar, origx, origy, origz, origx2, origy2, origz2);
					//Add difference to current total potential energy
					Ec+=dEp;
				}
			
						
				
		
				
				if(dccheck || !(Ec < Emedian)){
					/*Failed criteria--reset atom coordinates and remove
				  the difference in energy */ 
					config[medcind].dimer[Dr].atom[Ar].x = origx;
					config[medcind].dimer[Dr].atom[Ar].y = origy;
				
					config[medcind].dimer[Dr].atom[Ar2].x = origx2;
					config[medcind].dimer[Dr].atom[Ar2].y = origy2;
			
						Ec = Et1b;
				
				}
				else{
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
					if (art<0.10){
						step*=0.5;
						
					}
					else if (art<0.33){
						long double msc = art/0.33;
						step*=msc;
						
					}
					else if(art>0.51){
						long double msc = art/0.51;
						step*=msc;
					
					}
					accepttrans=0;
					attempttrans=0;
					if (arr<0.10){
						stepa*=0.5;
						
					}
					else if (arr<0.33){
						long double msc = arr/0.33;
						stepa*=msc;
						
					}
					else if(arr>0.51 && stepa<1.0){
						long double msc = arr/0.51;
						stepa*=msc;
						if(stepa>1.0){
							stepa=1.0;
						}
					}
					acceptrot=0;
					attemptrot=0;
				}
			

		} //End Markov Chain
		//test to make sure it found something with E<Emed
		if(!(Ec < Emedian)){
				//if not then set to lower interval config and launch new Markov Chain
				// randomly select a lower interval conf
				int randind = (int)(dsfmt_genrand_close_open(&dsfmt)*medcindex);
				int rcind = (int)(En[randind][1]);
				config[medcind].Equate(config[rcind]);
				Ec = En[randind][0];
				attempttrans=0;
				accepttrans=0;
				attemptrot=0;
				acceptrot=0;
				for (unsigned int m = 0; m < MCit; ++m)
				{
						//Store Initial Energy value
					long double Et1b = Ec;
					//Trial Move -----

					// randomly select an atom index	
					unsigned int Arandom = (unsigned int)(dsfmt_genrand_close_open(&dsfmt)*natoms);
					//Get the proper dimer index and atom index within the dimer
					unsigned long Dr = Arandom/2;
					unsigned long long Ar = Arandom - (Dr*2);
					//Store the original coordinates of randomly selected atom
					long double origx, origy, origz;
					origx = config[medcind].dimer[Dr].atom[Ar].x;
					origy = config[medcind].dimer[Dr].atom[Ar].y;
					origz = config[medcind].dimer[Dr].atom[Ar].z;
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
					origx2 = config[medcind].dimer[Dr].atom[Ar2].x;
					origy2 = config[medcind].dimer[Dr].atom[Ar2].y;
					origz2 = config[medcind].dimer[Dr].atom[Ar2].z;
		
					long double ranprob = dsfmt_genrand_close_open(&dsfmt);
					bool rottag = 0;
					//Only rotate 30% of trial moves
					if(ranprob<=0.30){
						//Pick random x, y, and z values
						long double u, v;
						u = dsfmt_genrand_close_open(&dsfmt);
						v = dsfmt_genrand_close_open(&dsfmt);
			
						//rotates about Ar
						RandomRotation2d(u, v, stepa, config[medcind].dimer[Dr], Ar);
						++attemptrot;
						rottag=1;
					}
					else{
						long double movex, movey;
						/* set random movements of the coordinates
					   of the selected atom */
						movex = (dsfmt_genrand_close_open(&dsfmt)-0.5)*step;
						movey = (dsfmt_genrand_close_open(&dsfmt)-0.5)*step;
			
						// Implement those moves
						config[medcind].dimer[Dr].atom[Ar].x += movex;
						config[medcind].dimer[Dr].atom[Ar].y += movey;
				
						config[medcind].dimer[Dr].atom[Ar2].x += movex;
						config[medcind].dimer[Dr].atom[Ar2].y += movey;
						++attempttrans;
			
					}	
					//check boundary
					bool dccheck = config[medcind].BoxCheck(Dr);
					//boundary passes
					if (!(dccheck)){
						//Check change in energy for pairwise interactions
						long double dEp=DeltaE(config[medcind], Dr, Ar, origx, origy, origz, origx2, origy2, origz2);
						//Add difference to current total potential energy
						Ec+=dEp;
					}
			
						
				
		
					
					if(dccheck || !(Ec < Emedian)){
						/*Failed criteria--reset atom coordinates and remove
					  the difference in energy */ 
						config[medcind].dimer[Dr].atom[Ar].x = origx;
						config[medcind].dimer[Dr].atom[Ar].y = origy;
				
						config[medcind].dimer[Dr].atom[Ar2].x = origx2;
						config[medcind].dimer[Dr].atom[Ar2].y = origy2;
			
							Ec = Et1b;
				
					}
					else{
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
						if (art<0.10){
							step*=0.5;
						
						}
						else if (art<0.33){
							long double msc = art/0.33;
							step*=msc;
						
						}
						else if(art>0.51){
							long double msc = art/0.51;
							step*=msc;
					
						}
						accepttrans=0;
						attempttrans=0;
						if (arr<0.10){
							stepa*=0.5;
						
						}
						else if (arr<0.33){
							long double msc = arr/0.33;
							stepa*=msc;
						
						}
						else if(arr>0.51){
							long double msc = arr/0.51;
							stepa*=msc;
					
						}
						acceptrot=0;
						attemptrot=0;
					}
				} //End Markov Chain
			
			}//End energy check
	//	cout<<" medcindex: "<<medcindex<<" Ec "<<Ec<<endl;	
		//Update energy in array
		En[medcindex][0]=Ec;

//Collected new sample set so.. 
//Sort Energies and get new median
//cout<<"Enmid before sort: "<<En[medcindex][0]<<endl;
dquickSort(En, 0, Samples-1);	
//cout<<"Enmid after sort: "<<En[medcindex][0]<<endl;
Emedian = En[medcindex][0];

//output
eoutfile<<Emedian<<endl;
medcind = (int)(En[medcindex][1]);
//observables to output

for (unsigned int eo = 0; eo < Samples; ++eo)
{
		//Calculate Observables				
	long double Rg = RadGyr(config[eo]);
	long double dynVol = (4.0/3.0)*pi*(Rg*Rg*Rg);
	long double dynDens = ((long double)(config[eo].ndimers)*2.0)/dynVol;
	long double Dor= NPvecProdNeighbor(config[eo]);
	long double Packparm;
	long double Npp, Npn, Nnn;
	Npp=0.0, Npn=0.0, Nnn=0.0;
	CountContacts(Npp, Npn, Nnn, config[eo], Packparm);




	//Push observables for Trajectory Averages
	long double Eetloc=En[eo][0];
	Estat.Push(Eetloc);
	Esqstat.Push(Eetloc*Eetloc);
	Rgstat.Push(Rg);
	Dstat.Push(dynDens);
	Dorstat.Push(Dor);
	Pparmstat.Push(Packparm);
	PPstat.Push(Npp);
	NPstat.Push(Npn);
	NNstat.Push(Nnn);
		
}
//Trajectory averages
long double Estatmean=Estat.Mean();
long double Esqstatmean=Esqstat.Mean();
long double Rgstatmean=Rgstat.Mean();
long double Dstatmean=Dstat.Mean();
long double Dorstatmean=Dorstat.Mean();
long double Pparmstatmean=Pparmstat.Mean();
long double PPstatmean=PPstat.Mean();
long double NPstatmean=NPstat.Mean();
long double NNstatmean=NNstat.Mean();
//output to file
nsobout << 0 << " " << Emedian << " " << Estatmean << " " << Esqstatmean << " ";
nsobout << Rgstatmean << " " << Dstatmean << " " << " " << Dorstatmean << " " << Pparmstatmean<< " ";
nsobout<<PPstatmean<<" "<<NPstatmean<<" "<<NNstatmean<<endl;
// Observables-- logfile output - same as screen output
		logfile << "* Trajectory Average Values:" << endl;
		logfile << "Energy: " << Estatmean << " and Energy Squared: " << Esqstatmean << endl;
		logfile << "Radius of gyration: " << Rgstatmean << endl;
		logfile << "Gyration Density: " << Dstatmean << endl;
		logfile << "Dimer-Dimer(NP-NP) Orientation Parameter " << Dorstatmean << endl;
		logfile << "HCP packing parameter: " << Pparmstatmean << endl;
		logfile << "Contacts; PP: "<<PPstatmean<<" NP: "<<NPstatmean<<" NN: "<<NNstatmean<<endl;
		// Observables-- logfile output - same as screen output
		cout << "* Trajectory Average Values:" << endl;
		cout << "Energy: " << Estatmean << " and Energy Squared: " << Esqstatmean << endl;
		cout << "Radius of gyration: " << Rgstatmean << endl;
		cout << "Gyration Density: " << Dstatmean << endl;
		cout << "Dimer-Dimer(NP-NP) Orientation Parameter " << Dorstatmean << endl;
		cout << "HCP packing parameter: " << Pparmstatmean << endl;
		cout << "Contacts; PP: "<<PPstatmean<<" NP: "<<NPstatmean<<" NN: "<<NNstatmean<<endl;
//reset running stat objects
Rgstat.Reset();
		Dstat.Reset();
		Dorstat.Reset();
		Estat.Reset();
		Esqstat.Reset();
		Pparmstat.Reset();
		PPstat.Reset();
		NPstat.Reset();
		NNstat.Reset();

//output median config
coordsfile<<natoms<<endl;
coordsfile<<"charged shell lj particles with Energy "<<Emedian<<endl;
for(unsigned long long hg = 0; hg<natoms; hg++){
	unsigned long Di = hg/2;
	unsigned long long Ai = hg - (2*Di);
	
	coordsfile<<setiosflags(ios::fixed)<<setprecision(15);
	string ttype = config[medcind].dimer[Di].atom[Ai].type;
	long double tx, ty, tz;
	tx = config[medcind].dimer[Di].atom[Ai].x;
	ty = config[medcind].dimer[Di].atom[Ai].y;
	tz = config[medcind].dimer[Di].atom[Ai].z;
	coordsfile<<ttype.c_str()<<" "<<tx<<" "<<ty<<" "<<tz<< endl;
}

//test for convergence
if(ii>1 && CompareDoubles2(Emedian, Emprev, 1.0e-6)){
	cout<<"Emedian: "<<Emedian<<endl;
	cout<<"System has converged. Ending simulation."<<endl;
	
	break;
} 
Emprev=Emedian;
} // End NS loop

logfile.close();
eoutfile.close();
nsobout.close();
//coutfile.close();
//stepoutfile.close();
logfile<<"Simulation Complete!"<<endl;
logfile<<"Output files are: "<<endl;
logfile<<"Nested Energies: "<<eoutname<<endl;
//logfile<<"Configuration Samples: "<<coutname<<endl;
logfile<<"Average Observable Values: "<<nsoboutname<<endl;

cout<<"Simulation Complete!"<<endl;
cout<<"Output files are: "<<endl;
cout<<"Logfile: "<<logfilename<<endl;
cout<<"Nested Energies: "<<eoutname<<endl;
//cout<<"Configuration Samples: "<<coutname<<endl;
cout<<"Average Observable Values: "<<nsoboutname<<endl;

cout<<endl;
cout << " " << endl;
cout << " Thank You and Have a Nice Day!" << endl;
cout << " " << endl;
cout << "... End of Line ..." << endl;
	return EXIT_SUCCESS;



}

// Define Functions

//Total Energy- sum of all pair interactions
long double TotalPairEnergy(const Config& c){
	long double Etot = 0.0;
	long double r;
	long double kappa = c.kappa;
	unsigned long long N = c.ndimers*2;
	//cout << "N " << N << endl;
		for(unsigned long long  pp = 0;pp<N-1;pp++){
			unsigned long Di1 = pp/2;
			unsigned long long Ai1 = pp - (Di1*2);
		//	cout << "pp " << pp;
			for(unsigned long long uu = pp+1;uu<N;uu++){
				unsigned long Di2 = uu/2;
				unsigned long long Ai2 = uu - (Di2*2);
			//	cout << " uu " << uu << endl;
				if(Di2==Di1){

					Etot += 0.0;
				}
				else{
					r = RDist(c.dimer[Di1].atom[Ai1], c.dimer[Di2].atom[Ai2]);
					Etot += PairSChargeShell(c.dimer[Di1].atom[Ai1], c.dimer[Di2].atom[Ai2], r, kappa)+PairLJColloid(c.dimer[Di1].atom[Ai1], c.dimer[Di2].atom[Ai2], r);
				
				}	
			}
		}
//	Etot *= 60.0;
	//cout << "Etot: " << Etot << endl;
	return(Etot);
}

// calculates the change in energy associated with one moved dimer
long double DeltaE(const Config& c, const unsigned& Dindex, const unsigned long long& Ar, const long double& ox1, const long double& oy1, const long double& oz1, const long double& ox2, const long double& oy2, const long double& oz2){
	long double Etot1=0.0;	
	long double Etoto1=0.0;
	long double Etot2=0.0;
	long double Etoto2=0.0;
	long double dE;
	unsigned long long Ar2;
	long double r1, r2, r1o, r2o;
	long double kappa = c.kappa;
	if (Ar==0)
	{
		Ar2=1;
	}	
	else{
		Ar2=0;
	}
	unsigned long long N=c.ndimers*2;
	for (unsigned int i = 0; i < N; ++i)
	{
		unsigned long Di = i/2;
		unsigned long long Ai = i - (Di*2);
	//	cout << "Di: " << Di << " Ai: " << Ai << endl;
		if(Di!=Dindex){
			r1 = RDist(c.dimer[Dindex].atom[Ar], c.dimer[Di].atom[Ai]);
			r2 = RDist(c.dimer[Dindex].atom[Ar2], c.dimer[Di].atom[Ai]);
			r1o = RDisto(c.dimer[Di].atom[Ai], ox1, oy1, oz1);
			r2o = RDisto(c.dimer[Di].atom[Ai], ox2, oy2, oz2);
			Etot1 += PairSChargeShell(c.dimer[Di].atom[Ai], c.dimer[Dindex].atom[Ar], r1, kappa)+PairLJColloid(c.dimer[Di].atom[Ai], c.dimer[Dindex].atom[Ar], r1);
			Etoto1 += PairSChargeShell(c.dimer[Di].atom[Ai], c.dimer[Dindex].atom[Ar], r1o, kappa)+PairLJColloid(c.dimer[Di].atom[Ai], c.dimer[Dindex].atom[Ar], r1o);
			Etot2 += PairSChargeShell(c.dimer[Di].atom[Ai], c.dimer[Dindex].atom[Ar2], r2, kappa)+PairLJColloid(c.dimer[Di].atom[Ai], c.dimer[Dindex].atom[Ar2], r2);
			Etoto2 +=PairSChargeShell(c.dimer[Di].atom[Ai], c.dimer[Dindex].atom[Ar2], r2o, kappa)+PairLJColloid(c.dimer[Di].atom[Ai], c.dimer[Dindex].atom[Ar2], r2o);
					
		}	
		else{
			Etot1 += 0.0;
			Etoto1 += 0.0;
			Etot2 += 0.0;
			Etoto2 += 0.0;
		}
	}
	dE = (Etot1+Etot2) - (Etoto1+Etoto2);
	//cout << "In DeltaE funct, dE is : " << dE << endl;
	//if(dE != dE){
		
	//}
//	cout<<"dE "<<dE<<endl;
	return( dE ); 
	

}

long double PairSChargeShell(const Atom& one, const Atom& two, const long double& r, const long double& kappa){
	
	//long double pi = 3.14159265359;
	long double potent;
	long double R1, R2, Ae1, Ae2, Ae;
	bool ctag1, ctag2;
//	long double k4 = pow(kappa, 4);
	ctag1 = (bool)(one.ctag), ctag2 = (bool)(two.ctag);
	R1 = one.rad, R2 = two.rad;
	Ae1 = one.Ae, Ae2 = two.Ae;
//	kappa =1.00;
	potent=0.0;
	if(ctag1 == ctag2){
		Ae = sqrt(Ae1*Ae2);
		//Ae = Ae1;
	}
	else{
		Ae = -sqrt(Ae1*Ae2);
	}
	long double rc = (4.0*(R1+R2))/kappa;

	//if(r<R1+R2){
	//	potent = 100000.0;
	//}
	 if(r<rc){
	 potent = (Ae/(r*kappa))*exp(-kappa*r);
	}
	else{
		potent=0.000000;
	}
//	cout<<"rc: "<<rc<<" kappa: "<<kappa<<" r: "<<r<<" Ae: "<<Ae<<endl;
//cout<<"charged pair: "<<potent<<endl;
	return(potent);
}

long double PairLJColloid(const Atom& one, const Atom& two, const long double& r){
	
//long double pi = 3.14159265359;
	long double potent;
	//long double R1s, R2s;
	long double Ac1, Ac2;
	//long double Acc, F;
	long double R1, R2;
	//long double rsm = one.rsmall;
	//long double rcon=one.sigc;
	string type1 = one.type;
	string type2 = two.type;
	R1 = one.rad, R2 = two.rad;
	Ac1 = one.Ac, Ac2 = one.Ac;
	long double eps = sqrt(Ac1*Ac2);
	long double rc = 2.0*(R1+R2);
	long double D = R1+R2;
	
//	if (r<0.5*D){
//		potent = 10000.0/r;
//	}
//	else
 if(r<rc){
		long double dor = D/r;
		long double pw50 = pow(dor, 50);
		potent = 4.0*eps*( pw50*pw50 - pw50 );
	}
	else{
		potent=0.0;
	}
//	cout << "PC potent: " << potent << endl;
	//if(potent != potent){
		//exit(0);
	//}
	//cout<<"colloid pair: "<<potent<<endl;
	return(potent);
}

//Calculates the distance between two atoms--wraps coordinates and gets minimum image
long double RDist(const Atom& one, const Atom& two){

	long double dx, dy, dz, r, rs;
	long double x1, y1, z1, x2, y2, z2;
	x1 = one.x, y1=one.y, z1=one.z;
	x2=two.x, y2=two.y, z2=two.z;
	
	dx = x1-x2;
	
	dy = y1-y2;
	
	dz = z1-z2;
	
	rs= (pow(dx, 2))+(pow(dy, 2))+(pow(dz, 2));
	r = sqrt(rs);
	return(r);



}
//Calculates the distance between an atom and the old coordinates of another
long double RDisto(const Atom& one, const long double& xo, const long double& yo, const long double& zo){

	long double dx, dy, dz, r, rs;
	long double x1, y1, z1, x2, y2, z2;
	x1 = one.x, y1=one.y, z1=one.z;
	x2 = xo, y2 = yo, z2 = zo;
	
	dx = x1-x2;
	
	dy = y1-y2;
	
	dz = z1-z2;
	
	rs= (pow(dx, 2))+(pow(dy, 2))+(pow(dz, 2));
	r = sqrt(rs);
	return(r);



}


long double RadGyr(Config& c){
	c.CalcCom();
	long double sumgyrx = 0.0, sumgyry = 0.0, sumgyrz = 0.0;
	long double M = 0.0;
	long double Rsq, R;
//	double ri, rcm, sumgyr=0.0;
	long double sumgyr=0.0;
//	cout << " M " << M << endl;
	long double comx, comy, comz;
	comx=c.COM[0];
	comy=c.COM[1];
	comz=c.COM[2];
	unsigned long long N = c.ndimers*2;
	for(unsigned int y = 0;y<N; ++y){
		unsigned long Di = y/2;
		unsigned Ai = y - (Di*2);
		double rx, ry, rz;
		rx = c.dimer[Di].atom[Ai].x;
		ry = c.dimer[Di].atom[Ai].y;
		rz = c.dimer[Di].atom[Ai].z;
		
//		ri = sqrt(rx*rx+ry*ry+rz*rz);
//		rcm = sqrt(com.x*com.x + com.y*com.y + com.z*com.z);
//		sumgyr += (ri-rcm)*(ri-rcm);
		long double mc = c.dimer[Di].atom[Ai].mass;
		M += mc;
		sumgyrx = (rx - comx)*(rx - comx);
		sumgyry = (ry - comy)*(ry - comy);
		sumgyrz = (rz - comz)*(rz - comz);
		sumgyr += (sumgyrx + sumgyry + sumgyrz)*mc;		
	}
//	cout << " sumgyrx " << sumgyrx ;
//	cout << " sumgyry " << sumgyry ;
//	cout << " sumgyrz " << sumgyrz << endl;
//	sumgyrx /= M;
//	sumgyry /= M;
//	sumgyrz /= M;
//	cout << " sumgyrx/M " << sumgyrx ;
//	cout << " sumgyry/M " << sumgyry ;
//	cout << " sumgyrz/M " << sumgyrz << endl;

//	Rsq = (sumgyrx*sumgyrx)+(sumgyry*sumgyry)+(sumgyrz*sumgyrz);
//	Rsq = sumgyrx + sumgyry + sumgyrz;
	Rsq = sumgyr;
//	cout << " Rsq " << Rsq << " M "<<endl;
//	Rsq = sumgyr;
	R = sqrt(Rsq/M);
//	cout << " RoG " << R << endl;
	if (R!=R){
		cout << " Rsq " << Rsq << " M "<<M<<endl;
		cout << " RoG " << R << endl;
		cout<<"comx: "<<comx<<" comy: "<<comy<<" comz: "<<comz<<endl;
		exit(0);
	}
	return(R);
}

long double LjWallZ(const Atom& one, const long double& WallZ ){

	long double Zwall, dz;
	//long double pi = 3.14159265359;
	long double rs;
	long double zc = one.z;
	long double rr = one.rad;
	long double potent, eps, sigma;
	// Energy Unit Reducer = {(Faraday's Constant*1 V)} * (1 Kcal/4184J) = Kcal/mol
	long double reducer = (96352.0/4184.0)/10.0;
	eps = 20.0*rr, sigma = rr+2.0;
	long double rc = 4.0*sigma;
	long double rcs = pow(rc, 2);
	Zwall = WallZ;
	dz = Zwall - zc;
	rs = pow(dz, 2.0);
	long double sigmasq = pow(sigma, 2);
	if ( rs < rcs){

		potent = 4.0*eps*((pow(sigmasq/rs, 6))-(pow(sigmasq/rs, 3)));
	}
	else{
		potent = 0.0;
	}
	
	potent /= reducer;
	return(potent);

} 

long double LjWallZp(const long double& oz, const long double& rad, const long double& WallZ){

	long double Zwall, dz;
	//long double pi = 3.14159265359;
//	long double zc = oz;
	long double rr = rad;
	long double rs;
	long double potent, eps, sigma;
	// Energy Unit Reducer = {(Faraday's Constant*1 V)} * (1 Kcal/4184J) = Kcal/mol
	long double reducer = (96352.0/4184.0)/10.0;
	eps = 20.0*rr, sigma = rr+2.0;
	long double rc = 4.0*sigma;
	long double rcs = pow(rc, 2);
	Zwall = WallZ;
	dz = Zwall - oz;
	rs = pow(dz, 2);
	long double sigmasq = pow(sigma, 2);
	if ( rs < rcs){

		potent = 4.0*eps*((pow(sigmasq/rs, 6))-(pow(sigmasq/rs, 3)));
	}
	else{
		potent = 0.0;
	}
	
	potent /= reducer;
	return(potent);

} 

long double TotalWallEnergyZ(const Config& c){
	unsigned long N = c.ndimers*2;
	long double wEtot=0.0;
	long double WallZ = c.WallZ;
	for (unsigned long i = 0; i < N; ++i)
	{
		unsigned long Di = i/2;
		unsigned long long Ai = i-(Di*2);
		wEtot+=LjWallZ(c.dimer[Di].atom[Ai], WallZ);
		
	}
	return(wEtot);
}

long double dEWallZ(const Dimer& d, long double oz1, long double r1, long double oz2, long double r2, const long double& WallZ){

	long double Eo=0.0, En=0.0;
	//unsigned short N = d.natoms;
	
	
	
	Eo= LjWallZp(oz1, r1, WallZ)+LjWallZp(oz2, r2, WallZ);
	En = LjWallZ(d.atom[0], WallZ) + LjWallZ(d.atom[1], WallZ);
	
	long double dE = En-Eo;
	return(dE);


}

long double AverageZorientation(const Config& c){
	unsigned long N=c.ndimers;
	long double Avgdz=0.0;
	long double BL = c.dimer[0].Lb;
	for (unsigned long i = 0; i < N; ++i)
	{
		//unsigned long Di = i/2;
		//Neg
		long double z1 = c.dimer[i].atom[1].z;
		//Pos
		long double z2 = c.dimer[i].atom[0].z;
		Avgdz += (z2-z1);		
	}
	Avgdz /= (BL*(long double)(N));
	//if(Avgdz!=Avgdz){
	//	cout<<"AverageZorientation Avgdz is: "<<Avgdz<<endl;
	//	cout<<"BL: "<<BL<<" N: "<<N<<endl;
	//}
	return(Avgdz);

}

void Zdensity(const Config& c, Histogram& gofz){
		unsigned long N = c.ndimers*2;
		for (unsigned long i = 0; i < N; ++i)
		{
			unsigned long Di = i/2;
			unsigned long long Ai = i-(Di*2);
			long double cz = c.dimer[Di].atom[Ai].z;
			gofz.Push(cz);
		}
		

}

long double NPvecProd(const Config& c){

	unsigned long N = c.ndimers;
	long double BL = c.dimer[0].Lb;
	long double BLs = BL*BL;
	long double Avg = 0.0;
	unsigned long count = 0;
	for (unsigned long i = 0; i < N-1; ++i)
	{
		long double dx1, dy1, dz1;
		dx1 = c.dimer[i].atom[0].x - c.dimer[i].atom[1].x;
		dy1 = c.dimer[i].atom[0].y - c.dimer[i].atom[1].y;
		dz1 = c.dimer[i].atom[0].z - c.dimer[i].atom[1].z;
		
		for (unsigned long j = i+1; j < N; ++j)
		{
			long double dx2, dy2, dz2;
			dx2 = c.dimer[j].atom[0].x - c.dimer[j].atom[1].x;
			dy2 = c.dimer[j].atom[0].y - c.dimer[j].atom[1].y;
			dz2 = c.dimer[j].atom[0].z - c.dimer[j].atom[1].z;
			
			Avg += (dx1*dx2+dy1*dy2+dz1*dz2)/(BLs);
			++count;
		}
	}
	Avg /= (long double)(count);
	//if(Avg!=Avg){
	//	cout<<"NPvecProd Avg is: "<<Avg<<endl;
	//	cout<<"BL: "<<BL<<" N: "<<N<<" count "<<count<<endl;
	//}
	return(Avg);		

}

long double NPvecProdNeighbor(const Config& c){

	unsigned long N = c.ndimers;
	long double BL = c.dimer[0].Lb;
	long double BLs = BL*BL;
	long double Avg = 0.0;
	unsigned long count = 0;
	long double avgdiam = c.dimer[0].atom[0].rad + c.dimer[0].atom[1].rad;
	long double cutoff = 3.0*avgdiam;

	for (unsigned long i = 0; i < N-1; ++i)
	{
		
		
		for (unsigned long j = i+1; j < N; ++j)
		{
			
			long double comdist = c.dimer[i].ComDistFrom(c.dimer[j]);
			if (comdist<cutoff){
				long double dx1, dy1, dz1;
				dx1 = c.dimer[i].atom[0].x - c.dimer[i].atom[1].x;
				dy1 = c.dimer[i].atom[0].y - c.dimer[i].atom[1].y;
				dz1 = c.dimer[i].atom[0].z - c.dimer[i].atom[1].z;

				long double dx2, dy2, dz2;
				dx2 = c.dimer[j].atom[0].x - c.dimer[j].atom[1].x;
				dy2 = c.dimer[j].atom[0].y - c.dimer[j].atom[1].y;
				dz2 = c.dimer[j].atom[0].z - c.dimer[j].atom[1].z;
			
				Avg += (dx1*dx2+dy1*dy2+dz1*dz2)/(BLs);
				++count;

			}

			
		}
	}
	if (count==0){
		++count;
	}
	Avg /= (long double)(count);
	//if(Avg!=Avg){
	//	cout<<"NPvecProd Avg is: "<<Avg<<endl;
	//	cout<<"BL: "<<BL<<" N: "<<N<<" count "<<count<<endl;
	//}
	return(Avg);		

}
void RandomRotation2d(const long double& u, const long double& v, const long double& stepa, Dimer& D, const unsigned long& Ar){
		long double theta;
        long double pi=3.14159265358979323;
		long double rho = D.Lb;
		long double dt, dp;
		unsigned long long Ar2;
		if (Ar==0)
		{
			Ar2=1;
		}
		else{
			Ar2=0;
		}
			
		dt=(2.0*pi)*stepa;
		//cout << "rho " << rho << "stepa " << stepa << " dt " << dt << " dp " << dp << endl;
		long double x1, y1, z1, x2, y2, z2;
		x1 = D.atom[Ar].x;
		y1 = D.atom[Ar].y;
		
		x2 = D.atom[Ar2].x;
		y2 = D.atom[Ar2].y;
	
	//	cout << "x1 " << x1 << " y1 " << y1 << " z1 " << z1 << endl;
	//	cout << "x2 " << x2 << " y2 " << y2 << " z2 " << z2 << endl;
		long double xo, yo, zo;
		xo=x2-x1; 
		yo=y2-y1;
		
		//long double phio = acos(zo/rho);
		long double ro = rho;
		long double thetao;//=acos(xo/ro);
		if (yo<0.0 && xo<0.0)
		{
			thetao=-asin(yo/ro)+pi;
		}
		else if(yo<0.0 && xo > 0.0){
			thetao = asin(yo/ro)+2.0*pi;
		}
		else{
			 thetao=acos(xo/ro);
		}
	//	cout << "xo: " << xo << " yo: " << yo << " zo: " << zo << endl;
	//	cout << "u " << u << " v " << v << endl;
	//	cout <<" Thetao: " << thetao << endl;
        theta=thetao+(u-0.5)*dt;
		if (theta>2.0*pi)
	{
		theta -= 2.0*pi;
	}
		else if(theta<0.0){
			theta+=2.0*pi;
		}

		  		
		//cout<<"theta "<<theta<<endl;
		long double xr, yr, rr;
	//	zr = rho*cos(phi);
		rr = rho;
		xr = rr*cos(theta);
		yr= rr*sin(theta);
		D.atom[Ar2].x = x1 + xr;
		D.atom[Ar2].y = y1 + yr;
//		D.atom[Ar2].z = z1 + zr;
	//	cout << "xr: " << xr << " yr: " << yr << " rr " << rr << endl; 
		//cout << "Ar2.x " << D.atom[Ar2].x << " Ar2.y " << D.atom[Ar2].y << " Ar2.z " << D.atom[Ar2].z << endl;
	//		cout << "----------------" << endl;
	if (x1+xr!=x1+xr){
			cout << "xr: " << xr << " yr: " << yr << " rr " << rr << endl; 
		cout << "Ar2.x " << D.atom[Ar2].x << " Ar2.y " << D.atom[Ar2].y << " Ar2.z " << D.atom[Ar2].z << endl;
		cout<<"theta "<<theta<<endl;
		cout <<" Thetao: " << thetao << endl;
		cout << "u " << u << " v " << v << endl;
		cout << "xo: " << xo << " yo: " << yo << " zo: " << zo << endl;
		cout << "x1 " << x1 << " y1 " << y1 << " z1 " << z1 << endl;
		cout << "x2 " << x2 << " y2 " << y2 << " z2 " << z2 << endl;
		cout << "rho " << rho << "stepa " << stepa << " dt " << dt << " dp " << dp << endl;
	}
		return;

}

long double AverageZwallDist(Config& c){

	unsigned Nd = c.ndimers;
	long double N = (long double)(Nd);
	long double WallZ = c.WallZ;
	long double sum = 0.0;
	for (unsigned int i = 0; i <Nd ; ++i)
	{
		c.dimer[i].CalcCom();
		long double diff = c.dimer[i].COM[2] - WallZ;
		sum += diff;
	}
	return(sum/N);
} 

string itos( int Number ){
     ostringstream ss;
     ss << Number;
     return ss.str();
 }

void CountContacts(long double& Npp, long double& Npn, long double& Nnn, Config& c, long double& packparm){

	long double r;
	long double R1, R2, rc;
	unsigned long long N = c.ndimers*2;
	long double packnum[N];
	for (unsigned int i = 0; i < N; ++i)
	{
		packnum[i]=0.0;
	}
	//cout << "N " << N << endl;
		for(unsigned long long  pp = 0;pp<N-1;pp++){
			unsigned long Di1 = pp/2;
			unsigned long long Ai1 = pp - (Di1*2);
			//cout << "pp " << pp;
			for(unsigned long long uu = pp+1;uu<N;uu++){
				unsigned long Di2 = uu/2;
				unsigned long long Ai2 = uu - (Di2*2);
			//	cout << " uu " << uu << endl;
				//cout<<"Di1: "<<Di1<<" Di2: "<<Di2<<" Ai1: "<<Ai1<<" Ai2: "<<Ai2<<endl;
				if(Di2!=Di1){

					R1 = c.dimer[Di1].atom[Ai1].rad;
					R2 = c.dimer[Di2].atom[Ai2].rad;
					rc = (R1+R2)*pow(202.0/51.0, 1.0/50.0);
					
					r = RDist(c.dimer[Di1].atom[Ai1], c.dimer[Di2].atom[Ai2]);
					//cout<<"r1: "<<R1<<" r2: "<<R2<<" rc: "<<rc<<" r: "<<r<<endl;
					//Check proper distance
					if (r < rc){
						packnum[pp]+=1.0;
						packnum[uu]+=1.0;
						//Check pp
						if (Ai1==0 && Ai2==0){
							Npp+=1.0;
						//	cout<<"increment Npp!"<<endl;
						}
						//check np
						else if (Ai1==1 && Ai2==0){
							Npn+=1.0;
							//cout<<"increment Npn!"<<endl;
						}
						//check pn
						else if (Ai1==0 && Ai2==1){
							Npn+=1.0;
							//cout<<"increment Npn!"<<endl;
						}
						// nn
						else{
							Nnn+=1.0;
							//cout<<"increment Nnn!"<<endl;
						}
						
					}
					
				
				}	
			}
		}
		long double avgpacknum=0.0;
		for (unsigned int i = 0; i < N; ++i)
		{
			avgpacknum+=packnum[i]/5.0;
		}
		avgpacknum/=(long double)(N);
		packparm=avgpacknum;
	return;

} 

bool CompareDoubles2 (long double A, long double B, long double tolerance) 
{
	long double EPSILON = tolerance;
   long double diff = A - B;
   return (diff < EPSILON) && (-diff < EPSILON);
}

void dquickSort(long double arr[][2], int left, int right) {

      int i = left, j = right;

      long double tmp1, tmp2;

      long double pivot = arr[(left + right) / 2][0];

 

      /* partition */

      while (i <= j) {

            while (arr[i][0] < pivot)

                  i++;

            while (arr[j][0] > pivot)

                  j--;

            if (i <= j) {

                  tmp1 = arr[i][0];
		  tmp2 = arr[i][1];

                  arr[i][0] = arr[j][0];
		  arr[i][1] = arr[j][1];

                  arr[j][0] = tmp1;
		  arr[j][1] = tmp2;

                  i++;

                  j--;

            }

      };

 

      /* recursion */

      if (left < j)

            dquickSort(arr, left, j);

      if (i < right)

            dquickSort(arr, i, right);

}
