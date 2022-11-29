/* Nested Sampling Routine--Collected nested median energys.
Configured for NP dimer system with Radii r1, r2, mass m1, m2, and bondlength Lb, where r1<=r2.
Uses combination of Colloid potential and screened electrostatics.
Has an LJ12-6 wall in the z-plane.

params and units:
 R1;p1 (distance units), R2;p2 (distance units), sigma (distance units), Ac (energy units), Ae (energy*distance units), kappa (1/distance units) 

 Using reduced units:
 R1 is reduced by sigma, but the other distance values are reduced by R1

 R1' = R1/sigma , R1'' = R1' / R1' = 1.0 , R2' = R2/sigma, R2'' = R2' / R1' = R2/R1 = 2.0
 distance' = distance/sigma, distance'' = distance'/R1' = distance/R1, Ac' = Ac/Ac,
 Ae' = Ae/(Ac*R1'*sigma), kappa' = kappa*sigma, kappa'' = kappa'*R1' = kappa*R1

 Note: R1os = R1' ; the value of sigma in the colloid potential should be replaced with 1/R1os, 
 i.e. sigma^6 => (1/R1os6)^6, to get correct unit reduction


The Monte Carlo trial move consist of moving a single dimer. Performs a translation move each time
and a rotation move only 30% of the time.
The boundary condition is set to hard cubic box for particle coordinates (i.e. does not directly account for radius of particle)
Contains a simple routine to adjust step size under each new median energy.

Note: Uses OpenMP multi-threading to run replicates in parallel. 
      Uses nested layer of parallelization to split the sampling 
	  loop (over MCit) across nprocs threads.
Must add the -fopenmp flag when compiling with g++
Calls some omp functions that require the -fopenmp flag for compilation
	omp_set_num_threads()
	omp_set_nested() -- use omp_set_nested(1) to allow nested parallelization

Safe Compile: 
 c++ -O2 -fopenmp prog.cpp
Faster Execution Speed Compiles: 
 c++ -O3 -fopenmp prog.cpp
 c++ -O2 -ffast-math -fcx-fortran-rules -fopenmp prog.cpp
 c++ -O3 -ffast-math -fcx-fortran-rules -fopenmp prog.cpp  

	
Modifications
 2/4/14
  -Fixed all -Wall warnings
   Used gprof profiler data to see that pair potential functions were eating up large amounts of call time %
  -Rewrote some calculations in the PairColloid function that have multiple calls in the computation to evaluate beforehand
  - Added cutoffs to both pair potentials (Note colloid potential dies fast so is safe for cutoff (rc=3.0*(R1+R2). For 
    charge shell checked potential for 
    r1=1 r2=1 and Ae1=5000 and Ae2=5000 and kappa=0.1 at same charge, which out of values used is one of the highest combos and
    should have one the larger cutoffs. 
    The potential dies well to zero (to second decimal or less) around r=60, so cutoff rc=(4.0*(R1+R2))/kappa 
	(which here would be 80.0) should 
    be more than sufficient. However, may need to adjust if going to even larger values of Ae.
  - PairChargeShell updated function input of r and kappa to const ref
  - PairColloid updated fucntion input of r to const ref	

  Summary: these modifications lead to an overall approximate 2.5x execution speed increase
  versus pre-modified code both compiled with c++ -O2 -fopenmp prog.cpp doing 1 million trial
  moves with 8 dimers; single rep with one cpu

 2/6/14
 -Added new member function to calculate NS average values of observables all at one (HistogramNS_slim::GetObsAvgAll
  Updated main code to call this instead of calling GetObsAvg for each separately.

  Summary: Makes execution over larger numbers of NS iterations run a little bit faster

 2/24/14
	-Added CountContacts function and updated code to do counting of contacts
	-Updated HistogramNS_slim with new GetObsAvgAll function to accomodate 3 new contact count variables
   Progfile updated to: pNest_dimers_parallel_replicates_nestedomp_3_pre.cpp
 	-Updated to a new Dump function in Histogram.h to take ofstream object with already opened file
	- Added Histogramming of Rog and Dor observables at each NS step
   Progfile updated to: pNest_dimers_parallel_replicates_nestedomp_4_pre.cpp

 3/19/14
	-Added explicit definitions of Ac (colloid energy prefactor) - using geometrical mixing; added Ac var to Atom 		 class from AtomD.h; Updated PairColloid function to use geometrical mixing to determin Acc. 

 3/19/14
	-Fully separate types of trial moves
	-Add counting and tracking for each type of trial move and modulation for each
	 Progfile updated to: pNest_dimers_parallel_replicates_nestedomp_5_pre.cpp
4-23-14
	-Split the OpenMP critical section into separate bits under each respective conditional
	-Updated computation of MCitpp and Configsoutfreq

6-25-14
	-Updated pairwise interaction functions - LJcolloid to LJ 100-50; Screened Electro - Yukawa Screened Function
		see Soft Matter, 2014, 10, 4479
	Progfile updated to: pNest_dimers_parallel_replicates_nestedomp_6_pre.cpp
	-Added fastPow function, but not using currently
6-26-14
	-Added in local walker minima tracking and usage
7-30-14
	-Previously converted to simulate 2d in xy
	- stripped out the outer level of parallelism. No longer parallel over number of
	  replicates and removes nested parallel. 

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
//#include <format.h>

// OpenMP Library
#ifdef _OPENMP
 #include <omp.h>
#endif
//RNG Fast Mersenne Twist- http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/
#include "/net/uu/nm/cm/bxw109120/src/MersenneTwist/dSFMT-src-2.2.1/dSFMT.c" 

// Specify the number of dimers
//#define n_dimers_n 40

using namespace std;

 // Classes -- cl
#include "/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/src/class/AtomD.h"
#include "/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/src/class/Dimer.h"
#include "/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/src/class/ConfigD.h"
#include "/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/src/class/HistogramNS_slim.h"
#include "/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/src/class/RunningStats.h"
#include "/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/src/class/Histogram.h"



// prototype functions -- pf

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
void RandomRotation(const long double& u, const long double& v, const long double& stepa, Dimer& D, const unsigned long& Ar);
long double AverageZwallDist(Config& c);
string itos( int Number );
void CountContacts(long double& Npp, long double& Npn, long double& Nnn, Config& c);
inline long double fastPow(long double a, long double b);
long double NPvecProdNeighbor(const Config& c);

//long double inline __declspec (naked) __fastcall sqrt14(long double n);

//Main Function -- bm
int main(int argc, char *argv[])
{
	//Get the local time for the logfile
	time_t rawtime;
 	struct tm * timeinfo;
 	time (&rawtime);
 	timeinfo = localtime (&rawtime);
	// Parameters that need to be set

	
	// the number of atoms to have in the system
	const unsigned short ndimers = 10;
	const unsigned natoms = ndimers*2; 
	
	// the number of iterations to run on collecting new samples
	const unsigned long int Iterations = 2000;
	//Number of Markov style iterations to adjust configuration
	const unsigned long long MCit = 4000000;
	//Number of Trial Moves To Initially not collect data--allow for equilibration
	const unsigned long Eqit = 1000000;
	/* Number of replicate sims to run in parallel 
	   -- for single run (i.e. no parallel runs) set to 1 */
	const unsigned short replicatenum = 1 ;
	//Number of threads to split the sampling loop across
	const unsigned short nprocs= 2 ;
	//Nested Sampling Integration Fraction
	const long double NSfrac = 0.50;
	//Dimer Properties
	const long double rad1 = 1.0; // Positive charged
	const bool ctag1 = 1;
	const long double Ae1 = 1.0;
	const long double Ac1 = 1.0;
	const long double mass1 = 1.0;
	string type1 = "P";
	const long double rad2 = 1.0; // Negative charged
	const bool ctag2 = 0;
	const long double Ae2 = 1.0;
	const long double Ac2 = 1.0;
	const long double mass2 = 1.0;
	string type2 = "N";	
	const long double BondLength = rad1+rad2;
	//long double sigconv = 50.0;
	// kappa value - inverse screening length	
	const long double kappa = 0.005 ;
	// Hard Cubic box size
	const long double box = ( (long double)(ndimers) )*(rad1+rad2)*2.0 ;
	// Set for initial coordinates
	long double boxxx = box, boxyy = box;
	//z-coordinate to place interacting surface
	const long double WallZ = -box/2.0;
	const long double pi = 3.14159265359;		
	//Movement Step sizelong double Eavgns = Ehist.GetObsAvg(0);
	const long double stepguess = 0.90*box;	
	long double step = stepguess;
	const long double stepaguess = 1.0/12.0;
	long double stepa = stepaguess;
	//Energy Histogram Properties
	unsigned long long nbins;
	if(ndimers<=10){
		nbins=6000000;
	}
	else {
		nbins = 6000000+(((ndimers-10)*500000));
	}
	const long double Elow = (-10.0)*(long double)(ndimers);
	const long double Ehigh = 200.0;
	//Observable histogram properties
	// Rog
	const long double Roglow = 1.0;
	const long double Roghigh = box/1.5;
	const unsigned Rognbins = 1000;
	// Dor
	const long double Dorlow = -1.0;
	const long double Dorhigh =1.0;
	const unsigned Dorbins = 200;
	//
		
	//collect step
	const unsigned collectstepguess = ndimers*5;
	unsigned collectstep = collectstepguess;
	//run descriptor
	string nds = itos(ndimers);
	string rundescriptor = "cube_d"+nds+"_test_24";
	
	//----End of Parameter setting-----

	long double rsmall;
	if(rad1<rad2){
		rsmall=rad1;
	}
	else if(rad2<rad1){
		rsmall=rad2;
	}
	else {
		rsmall=rad1;
	}
	unsigned long long MCitpp = (MCit+Eqit*nprocs)/nprocs;
	unsigned int ConfigsOutFreq = ((MCitpp-Eqit))/2;

cout<<"Will run replica " <<replicatenum<<" of system:"<<endl;
cout <<"---- " <<ndimers<<" dimers and a total of " << natoms << " atoms" << endl;




	//Output Files
	string rundescript = rundescriptor;
	unsigned short th_id;
	
	string t_id = itos(replicatenum);
	rundescript+="_rep"+t_id;
	//Name of logfile
	string logfilename = "oMMC_" + rundescript + ".log";
	//Name of Emedian output file
	string eoutname = "Epoints_int_" + rundescript + ".dat";
	//Name of Configuration output file
	string coutname = "ConfigsDimer_int_" + rundescript + ".xyz";
	//Name of Observable output file
	string outfilename = "zOutput_int_" + rundescript + ".dat";
	//Name of NS values file for observables
	string nsoboutname = "zNSobs_int_" + rundescript + ".dat";
	//Name of Observable Histrogram Distribution
	string rogdistname = "zNS_dist_rg_int_"+rundescript+".dat";
	string dodistname = "zNS_dist_do_int_"+rundescript+".dat";
	
		cout << "Running Replicate: "<<th_id<<endl;
		cout << "Name of logfile is : " << logfilename.c_str() << endl;
	
	//#pragma omp threadprivate(coordsfile)
   //Initialize output file streams
	
	ofstream logfile(logfilename.c_str());
	ofstream eoutfile(eoutname.c_str());
	
	ofstream outfile(outfilename.c_str());
	ofstream nsobout(nsoboutname.c_str());	
	ofstream rgout(rogdistname.c_str());
	ofstream doout(dodistname.c_str());
	ofstream coordsfile;
	//#pragma omp threadprivate(coordsfile)
	coordsfile.open(coutname.c_str());
	
	logfile << " " << endl;
	logfile << "Local time and date: " << asctime(timeinfo) << endl;
	logfile << "Parameters: " << endl;
	logfile << "Nested Sampling Energy Integration Fraction: " << NSfrac << endl;
	logfile << "R1: " << rad1 << " R2: " << rad2 << " ctag1: " << ctag1 << " ctag2 " << ctag2 << endl;
	logfile << "Ae1: " << Ae1 << " Ae2: " << Ae2 << " mass1: " << mass1 << " mass2: " << mass2 << endl;
	logfile<<"Ac1: "<<Ac1<<" Ac2: "<<Ac2<<endl;
	logfile << "kappa: " << kappa << endl;
	logfile << "box: " << box << " WallZ: " << WallZ << " stepguess: " << stepguess << " stepaguess: " << stepaguess << endl;
	logfile << "Energy Histogram Bins: " << nbins << " Elow: " << Elow << " Ehigh: " << Ehigh << endl;
	logfile<<"Radius of Gyration Histogram Bins: "<<Rognbins<<" Roglow: "<<Roglow<<" Roghigh: "<<Roghigh<<endl;
	logfile<<"Dimer Orientation Param Histogram Bins: "<<Dorbins<<" Dorlow: "<<Dorlow<<" Dorhigh: "<<Dorhigh<<endl;
	logfile << "collectstepguess: " << collectstepguess << endl;
	logfile << "Number of Equilibration Steps: " << Eqit << endl;
	logfile<<"ConfigsOutFreq: "<<ConfigsOutFreq<<" MCitpp: "<<MCitpp<<endl;
	logfile<<"Replicate number is set to : "<<replicatenum<<" and number of procs for sampling parallelization is: "<<nprocs<<endl;
	logfile << "System running with " << ndimers << " dimers and has a total of " << natoms << " atoms" << endl;

	cout << "Running...."<<endl;
	

	// Initialize the MT RNG object
	dsfmt_t dsfmt;	
	//Define a seed
	//time(0) is the clock value in seconds
	unsigned long long seed = time(0)+(th_id*2)+(th_id*14);
	logfile<<"FMT Seed value: "<<seed<<endl;
	// Initialize the RNG object with the seed
	
	dsfmt_init_gen_rand(&dsfmt, seed);

	//Create a Configuration Object
	Config confn(ndimers);
	confn.SetBox(box, box, box);
	confn.WallZ = WallZ;
	confn.SetKappa(kappa);
	
	logfile << "Dimer Bond Length: " << BondLength << endl;
	logfile<<endl;
	//Assign dimer parameters
	logfile<<"Assigning particle parameters..."<<endl;
	for(unsigned int jj = 0;jj<ndimers;jj++){
			
			confn.dimer[jj].SetBond(BondLength);
			confn.dimer[jj].atom[0].type = type1;
			confn.dimer[jj].atom[0].rad = rad1;
			confn.dimer[jj].atom[0].mass = mass1;	
			confn.dimer[jj].atom[0].ctag = ctag1;	
			confn.dimer[jj].atom[0].Ae = Ae1;	
			confn.dimer[jj].atom[0].Ac = Ac1;			
			confn.dimer[jj].atom[1].type = type2;
			confn.dimer[jj].atom[1].rad = rad2;
			confn.dimer[jj].atom[1].mass = mass2;
			confn.dimer[jj].atom[1].ctag = ctag2;	
			confn.dimer[jj].atom[1].Ae = Ae2;	
			confn.dimer[jj].atom[1].Ac = Ac2;
	}
	//Assign an Initial Random Configuration 
	logfile<<"Assigning an initial random configuration in the box..."<<endl;
	bool atag = 0;
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
			x = ranx*(boxxx);
			y = rany*(boxyy);
			//z = ranz*(boxzz);
			confn.dimer[jj].atom[0].x=x;
			confn.dimer[jj].atom[0].y=y;
			confn.dimer[jj].atom[0].z= WallZ+rad1;
			long double u;
			u = dsfmt_genrand_close_open(&dsfmt);
			//v = dsfmt_genrand_close_open(&dsfmt);		
			long double theta = u*2.0*pi;
			//long double phi = acos((2.0*v) - 1.0);
			long double BondLength = confn.dimer[jj].Lb;	
			//long double dz = BondLength*cos(phi);
			//long double r = sqrt(BondLength*BondLength-dz*dz);
			long double r = BondLength;
			long double dx = r*cos(theta);
			long double dy = r*sin(theta);
			confn.dimer[jj].atom[1].x = x+dx;
			confn.dimer[jj].atom[1].y = y+dy;
			confn.dimer[jj].atom[1].z = WallZ+rad2;
			
		} // end assigment for
		
		confn.CalcCom();
		//Test boundary condition
		if(confn.BoxCheck() == 1){
			atag = 0;
		}
		long double Et = TotalPairEnergy(confn);//+TotalWallEnergyZ(confn);
		if(Et>Ehigh-0.0010 || Et<Elow+0.0010){
				atag=0;
		}		
	} // end assignment while
	//Assign the energy of the initial configuration
	long double Et1 = TotalPairEnergy(confn);//+TotalWallEnergyZ(confn);
	Config Bottom(confn);
	// local min confs
	Config * lgmin  = new Config[nprocs];

	// Initialize an array to store the lgmin Energy values
	long double Enmin[nprocs]={0.0};

	logfile << "Energy of Initial Configuration is: " << Et1 << endl;
	//Initialize the Energy Histogram	
	HistogramNS Ehist;
	Ehist.Initialize(nbins, Elow, Ehigh);
	Ehist.InitializeObservables(10);
	Ehist.SetNSfrac(NSfrac);
	//Initialize Observable histograms
	Histogram Rghist;
	Rghist.Initialize(Rognbins, Roglow, Roghigh);
	Histogram Dohist;
	Dohist.Initialize(Dorbins, Dorlow, Dorhigh);
	//Track the radius of gyration and gyration density
	RunningStats Rgstat;
	RunningStats Dstat;
	// Average Energy and Average Energy^2
	RunningStats Estat;
	RunningStats Esqstat;
	//Orientation OPs
	RunningStats Dorstat;
	RunningStats Zorstat;
	//Zwall particel distance
	RunningStats Zdstat;
	//N-N contacts
	RunningStats PPstat;
	//N-P contacts
	RunningStats NPstat;
	//P-P contacts
	RunningStats NNstat; 
	/*Begin iterating down phase space in 0.5 intervals collecting
	MC trajectories and histogramming energies under each new median energy*/
	logfile << "Nested sampling iterations is set to " << Iterations << " with MC iterations set to " << MCit << endl;
	/*Define the first median energy--typically set
  	to the upper limit of the Energy histogram */
	long double Emedian = Ehigh;

	long double Emprev = Emedian;
	//Used to track last accepted energy for histogram
	for (unsigned int np = 0; np < nprocs; ++np)
	{
		lgmin[np].InitializeDimers(ndimers);
		lgmin[np].Equate(Bottom);
		Enmin[np]=Et1;
	}	

	Et1 = TotalPairEnergy(confn);//+TotalWallEnergyZ(confn);
	logfile << setiosflags(ios::fixed) << setprecision(15);
//Variable to track the Bottom Config
	long double Elowprev = Et1;


// Begin The Nested Sampling Iterations -- biter
	for (unsigned int i = 0; i < Iterations; ++i){
	
		logfile << "Iteration " << i << "......." << endl;
		//#pragma omp flush
		
	 if(i>0){
		Emprev = Emedian;
		Emedian = Ehist.NestedEnergy();
		logfile << setiosflags(ios::fixed) << setprecision(15) <<"*** Emedian: " << Emedian << endl;
		eoutfile << setiosflags(ios::fixed) << setprecision(15) << i << " " << Emedian << endl;

		long double Eavgns;
		long double Esqavgns;
		long double Rgns;
		long double Drgns;
		long double Zorns;
		long double Dorns;
		long double Zdns;
		long double PPns, NPns, NNns;
		Ehist.GetObsAvgAll(Eavgns, Esqavgns, Rgns, Drgns, Zorns, Dorns, Zdns, PPns, NPns, NNns);
		//Trajectory averages
		long double Estatmean=Estat.Mean();
		long double Esqstatmean=Esqstat.Mean();
		long double Rgstatmean=Rgstat.Mean();
		long double Dstatmean=Dstat.Mean();
		long double Zorstatmean=Zorstat.Mean();
		long double Dorstatmean=Dorstat.Mean();
		long double Zdstatmean=Zdstat.Mean();
		long double PPstatmean=PPstat.Mean();
		long double NPstatmean=NPstat.Mean();
		long double NNstatmean=NNstat.Mean();
	
	// Observables-- logfile output - same as screen output
		logfile << "* Trajectory Average Values:" << endl;
		logfile << "Energy: " << Estatmean << " and Energy Squared: " << Esqstatmean << endl;
		logfile << "Radius of gyration: " << Rgstatmean << endl;
		logfile << "Gyration Density: " << Dstatmean << endl;
		logfile << "Z-orientation Parameter: " << Zorstatmean << endl;
		logfile << "Dimer-Dimer(NP-NP) Orientation Parameter " << Dorstatmean << endl;
		logfile << "Zwall distance: " << Zdstatmean << endl;
		logfile << "Contacts; PP: "<<PPstatmean<<" NP: "<<NPstatmean<<" NN: "<<NNstatmean<<endl;
		logfile << "* NS Average Values (between NS energy at n and n-1):" << endl;
		logfile << "Energy: " << Eavgns << " and Energy Squared: " << Esqavgns << endl;
		logfile << "Radius of Gyration: " << Rgns << endl;
		logfile << "Gyration Density: " << Drgns << endl;
		logfile << "Z-orientation Parameter: " << Zorns << endl;
		logfile << "Dimer-Dimer(NP-NP) Orientation Parameter: " << Dorns << endl;
		logfile << "Zwall distance: " << Zdns << endl;
		logfile << "Contacts; PP: "<<PPns<<" NP: "<<NPns<<" NN: "<<NNns<<endl;
		// Observables-- output to files
		// Trajectory Averages
		outfile << setiosflags(ios::fixed) << setprecision(15) << i << " " << Emedian << " " << Estatmean << " " << Esqstatmean << " ";
		//outfile << Rgstatmean << " " << Dstatmean << " " << Zorstat.Mean() << " " << Dorstat.Mean() << endl;
		outfile << Rgstatmean << " " << Dstatmean << " " << Zorstatmean << " " << Dorstatmean << " " << Zdstatmean<< " ";
		outfile<<PPstatmean<<" "<<NPstatmean<<" "<<NNstatmean<<endl;
		// NS Averages
		nsobout << setiosflags(ios::fixed) << setprecision(15) << i << " " << Emedian << " " << Eavgns << " " << Esqavgns << " ";
		nsobout << Rgns << " " << Drgns << " " << Zorns << " " << Dorns << " " << Zdns << " ";
		nsobout<<PPns<<" "<<NPns<<" "<<NNns<<endl;
		//Output Observable histograms
		Rghist.Dump(rgout);
		Dohist.Dump(doout);


		// Reset Histogram and RunningStats
		Ehist.SmartReset();
		Rgstat.Reset();
		Dstat.Reset();
		Zorstat.Reset();
		Dorstat.Reset();
		Estat.Reset();
		Esqstat.Reset();
		Zdstat.Reset();
		PPstat.Reset();
		NPstat.Reset();
		NNstat.Reset();
		Rghist.Reset();
		Dohist.Reset();
		if (abs(Emprev-Emedian)<0.01){
			logfile<<"Convergence criteria met. Exiting simulation."<<endl;
			break;
		}
	
	//Set the new starting config
		for (unsigned int np = 0; np < nprocs; ++np)
		{	
			if (Enmin[np]>Emedian){
				lgmin[np].Equate(Bottom);
				Enmin[np]=Elowprev;
			
			}
		}

		logfile << "* Energy of new initial configuration is " << Et1 << endl;
	} //end output if
	cout<<"Elowprev: "<<Elowprev<<" and Enest is "<<Emedian<<endl;
	long double percent = 0.0;
	unsigned long long rotcount=0;
	unsigned long long transcount=0;
	unsigned long long srotcount=0;
	unsigned long long stranscount=0;

	// Parallell
	
	omp_set_num_threads(nprocs);
	//#pragma omp parallel shared(coordsfile, Ehist, Estat, Esqstat, Rgstat, Dstat, Zorstat, Dorstat, Zdstat, PPstat, NPstat, NNstat)
	#pragma omp parallel
	{
		//#pragma omp flush
		//local thread number
		 unsigned tid = omp_get_thread_num();
		Config Tlocal(confn);
		long double Eetloc=0.0;
		long double Etloc=0.0;	
		 //set starting conf
		 if (i==0){
		 	Tlocal.Equate(Bottom);
			Eetloc=Elowprev;
			Etloc=Elowprev;
		 }
		 else{
			Tlocal.Equate(lgmin[tid]);
			Eetloc = Enmin[tid];
			Etloc = Enmin[tid];
			
		 }
		Etloc = TotalPairEnergy(Tlocal);
		cout<<"--------starting energy is "<<Etloc<<endl;
		long double percentlocal=0.0;
		unsigned long long  rotcountlocal = 0;
		unsigned long long  transcountlocal = 0;
		unsigned long long  srotcountlocal = 0;
		unsigned long long  stranscountlocal = 0;
		//Begin the MC trajectory under current Emedian -- bmc
		for (unsigned long long int j = 0; j < MCitpp; ++j){
		
			//Save the Energy before trial move	
			long double Et1s = Etloc;
		
			// randomly select an atom index	
		
			unsigned int Arandom = (unsigned int)(dsfmt_genrand_close_open(&dsfmt)*natoms);
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
			//tag which trial move: rotation or translation
			bool mtag;
			//Only rotate 30% of trial moves
			if(ranprob<=0.30){
				//Pick random x, y, and z values
				long double u, v;
				u = dsfmt_genrand_close_open(&dsfmt);
				v = dsfmt_genrand_close_open(&dsfmt);
			
				//rotates about Ar
				RandomRotation(u, v, stepa, Tlocal.dimer[Dr], Ar);
				mtag=0;
				++rotcountlocal;
			}
			else{
				long double movex, movey;
				/* set random movements of the coordinates
			   of the selected atom */
				movex = (dsfmt_genrand_close_open(&dsfmt)-0.5)*step;
				movey = (dsfmt_genrand_close_open(&dsfmt)-0.5)*step;
			//	movez = (dsfmt_genrand_close_open(&dsfmt)-0.5)*step*0.5;
				// Implement those moves
				Tlocal.dimer[Dr].atom[Ar].x += movex;
				Tlocal.dimer[Dr].atom[Ar].y += movey;
				//Tlocal.dimer[Dr].atom[Ar].z += movez;
				Tlocal.dimer[Dr].atom[Ar2].x += movex;
				Tlocal.dimer[Dr].atom[Ar2].y += movey;
				//Tlocal.dimer[Dr].atom[Ar2].z += movez;
				mtag=1;
				++transcountlocal;
			
			}	

		
			long double dE = 0.00000000000000;

			bool tagg = 1;
			bool dccheck = Tlocal.BoxCheck(Dr);
		
			bool wallcheck = 0;
//			if(nz1<=WallZ+r1+0.100 || nz2 <= WallZ+r2+0.100){
//				
//				wallcheck = 1;
//			}
		
			bool hardsuccess = 0;
	 		if(wallcheck == 0 && dccheck==0){
				hardsuccess = 1;
			
				long double dEpp = DeltaE(Tlocal, Dr, Ar, origx, origy, origz, origx2, origy2, origz2);
	
			//	long double dEww = dEWallZ(Tlocal.dimer[Dr], origz, r1, origz2, r2, WallZ);
		
				dE = dEpp;// + dEww;
			//	cout<<"dE "<<dE<<" Etloc "<<Etloc<<" Etloc+dE "<<Etloc+dE<<endl;
				Etloc+=dE;
				//long double Eneww = TotalPairEnergy(Tlocal);
			//	cout<<"New energy is: "<<Eneww<<endl;			
			}
			
			bool softsuccess = 0;
			
			if(Etloc>Elow && Etloc<Emedian){
				softsuccess = 1;
			}
		
			
			if (hardsuccess == 1 && softsuccess == 1){
			
				tagg=0;
				percentlocal += 1.0;
				if (mtag==0){
					++srotcountlocal;
				}
				else{
					++stranscountlocal;
				}
		
				Eetloc = Etloc;
			
				if (Etloc<Enmin[tid]){
				
			
					Enmin[tid]=Etloc;
					lgmin[tid].Equate(Tlocal);
				
				}
				if(Etloc<Elowprev){
					#pragma omp critical
					{
						Bottom.Equate(Tlocal);
						Elowprev = Etloc;
					
					}//end critical
				}
						
			}// end trial move acceptance if
			/*Failed criteria--reset atom coordinates and remove
			  the difference in energy */ 
			if(tagg == 1){
			
				Tlocal.dimer[Dr].atom[Ar].x = origx;
				Tlocal.dimer[Dr].atom[Ar].y = origy;
			//	Tlocal.dimer[Dr].atom[Ar].z = origz;
				Tlocal.dimer[Dr].atom[Ar2].x = origx2;
				Tlocal.dimer[Dr].atom[Ar2].y = origy2;
			//	Tlocal.dimer[Dr].atom[Ar2].z = origz2;
			
				if(hardsuccess == 1){
					Etloc=Et1s;
				}
			
			} // end failed criteria if
		
			//dE = 0.0;
		
			if(j>Eqit && j%(collectstep)==0){

					//Calculate Observables				
					long double Rg = RadGyr(Tlocal);
					long double dynVol = (4.0/3.0)*pi*(Rg*Rg*Rg);
					long double dynDens = ((long double)(Tlocal.ndimers)*2.0)/dynVol;
					long double Zor = AverageZorientation(Tlocal);
					long double Dor= NPvecProdNeighbor(Tlocal);
					long double Zdist = AverageZwallDist(Tlocal);
					long double Npp, Npn, Nnn;
					Npp=0.0, Npn=0.0, Nnn=0.0;
					CountContacts(Npp, Npn, Nnn, Tlocal);
					//Flush relevant constructs
					#pragma omp flush(Ehist, Estat, Esqstat, Rgstat)
					#pragma omp flush(Dstat, Zorstat, Dorstat, Zdstat)
					#pragma omp flush(PPstat, NPstat, NNstat, Rghist, Dohist)

					#pragma omp critical
					{
					//	cout<<"############# Pushing energy: "<<Eetloc<< " Etloc is "<<Etloc<<endl;
						Ehist.Push(Eetloc);
						//Push observables for Trajectory Averages
						Estat.Push(Eetloc);
						Esqstat.Push(Eetloc*Eetloc);
						Rgstat.Push(Rg);
						Dstat.Push(dynDens);
						Zorstat.Push(Zor);
						Dorstat.Push(Dor);
						Zdstat.Push(Zdist);
						PPstat.Push(Npp);
						NPstat.Push(Npn);
						NNstat.Push(Nnn);
						//Push observables for NS Averages
						Ehist.PushObs(Eetloc, 0, Eetloc);
						Ehist.PushObs(Eetloc*Eetloc, 1, Eetloc);
						Ehist.PushObs(Rg, 2, Eetloc);
						Ehist.PushObs(dynDens, 3, Eetloc);
						Ehist.PushObs(Zor, 4, Eetloc);
						Ehist.PushObs(Dor, 5, Eetloc);
						Ehist.PushObs(Zdist, 6, Eetloc);
						Ehist.PushObs(Npp, 7, Eetloc);
						Ehist.PushObs(Npn, 8, Eetloc);
						Ehist.PushObs(Nnn, 9, Eetloc);
						// Observable Histograms
						Rghist.Push(Rg);
						Dohist.Push(Dor);
	
					
				}//End Critical	
			
			} //end value pushing if
		
			/* Output coordinates from
	  		 MC trajectory in an xyz file format */
			if(j>Eqit && j%ConfigsOutFreq ==0){
				//#pragma omp barrier
				//#pragma omp flush(coordsfile)
				#pragma omp critical(outer)
				{
				//	unsigned thrid = omp_get_thread_num();
				//	cout<<"outputting to coordfile from thread: "<<thrid<<endl;
					unsigned long nna = natoms;
					coordsfile<<nna<<endl;
					//coordsfile.flush();
					coordsfile<<"charged shell lj particles with Energy "<<Etloc<<endl;
					//coordsfile.flush();
					for(unsigned long long hg = 0; hg<natoms; hg++){
						unsigned long Di = hg/2;
						unsigned long long Ai = hg - (2*Di);
						
						coordsfile<<setiosflags(ios::fixed)<<setprecision(15);
						string ttype = Tlocal.dimer[Di].atom[Ai].type;
						long double tx, ty, tz;
						tx = Tlocal.dimer[Di].atom[Ai].x;
						ty = Tlocal.dimer[Di].atom[Ai].y;
						tz = Tlocal.dimer[Di].atom[Ai].z;
						coordsfile<<ttype.c_str()<<" "<<tx<<" "<<ty<<" "<<tz<< endl; 
						//coordsfile.flush();
					}//end atom loop
				
					
			 	 }//End Critical
			}// end coordinate out if

		} //End of MCit Loop -- emc
		#pragma omp barrier

		#pragma omp critical
		{
		 rotcount+=rotcountlocal;
		 transcount+=transcountlocal;
		 srotcount+=srotcountlocal;
		 stranscount+=stranscountlocal;
		}

	}//End Nested parallel
	
//	Et1 = Elowprev;
//			#pragma omp critical
//			{
//				unsigned long nna = ndimers;
//				coordsfile << nna << endl;
//				coordsfile << "charged shell lj particles with Energy " << Et1 << endl;
//				for(unsigned long long hg = 0; hg<natoms; hg++){
//					unsigned long Di = hg/2;
//					unsigned long long Ai = hg - (2*Di);
//	
//					coordsfile<<setiosflags(ios::fixed)<<setprecision(15)<<Bottom.dimer[Di].atom[Ai].type.c_str()<<" "<<Bottom.dimer[Di].atom[Ai].x<<" "<<Bottom.dimer[Di].atom[Ai].y<<" "<<Bottom.dimer[Di].atom[Ai].z<<endl; 
//				}
//			}//end critical
	
	
//Calc percent
percent /= (long double)(MCit)/100.0;
long double rper = (long double)(srotcount)/(long double)(rotcount);
long double tper = (long double)(stranscount)/(long double)(transcount);
//Adjust step size if needed
if (rper<0.30){
	stepa*=0.98;
}
if (tper<0.30){
	step*=0.95;
}

//cout<<"Elowprev: "<<Elowprev<<endl;
logfile<<"rper: "<<rper<<" stepa: "<<stepa<<" tper: "<<tper<<" step: "<<step<<endl;
if (rper<0.15 && tper<0.15){
	logfile << "!!!*** Configuration is frozen. Exiting simulation. ***!!!" << endl;
		break;
}

logfile << endl;	

} //End of Iteration Loop -- eiter

eoutfile.close();
coordsfile.close();
outfile.close();
nsobout.close();

cout << "Rep "<<th_id<<":"<<endl;
cout << "Output files are: " << endl;
cout << eoutname.c_str() << endl;
cout << coutname.c_str() << endl;
cout << outfilename.c_str() << endl;
//cout << gofzofname.c_str() << endl;
cout << nsoboutname.c_str() << endl;

logfile << "Output files are: " << endl;
logfile << eoutname.c_str() << endl;
logfile << coutname.c_str() << endl;
logfile << outfilename.c_str() << endl;
//logfile << gofzofname.c_str() << endl;
logfile << nsoboutname.c_str() << endl;

cout << " " << endl;
cout << " Thank You and Have a Nice Day!" << endl;
cout << " " << endl;
cout << "... End of Line ..." << endl;


return EXIT_SUCCESS;

} // End of Main -- eom

// Define Functions -- df

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
	
	if (r<D){
		potent = 10000.0;
	}
	else if(r<rc){
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
		M += c.dimer[Di].atom[Ai].mass;
		sumgyrx = (rx - comx)*(rx - comx);
		sumgyry = (ry - comy)*(ry - comy);
		sumgyrz = (rz - comz)*(rz - comz);
		sumgyr += sumgyrx + sumgyry + sumgyrz;		
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
//	cout << " Rsq " << Rsq << endl;
//	Rsq = sumgyr;
	R = sqrt(Rsq/M);
//	cout << " RoG " << R << endl;
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
void RandomRotation(const long double& u, const long double& v, const long double& stepa, Dimer& D, const unsigned long& Ar){
		long double theta, phi;
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
		phi=0.0;
		
		dt=(2.0*pi)*stepa;
		dp = pi*stepa;
		//cout << "rho " << rho << "stepa " << stepa << " dt " << dt << " dp " << dp << endl;
		long double x1, y1, z1, x2, y2, z2;
		x1 = D.atom[Ar].x;
		y1 = D.atom[Ar].y;
		z1 = D.atom[Ar].z;
		x2 = D.atom[Ar2].x;
		y2 = D.atom[Ar2].y;
		z2 = D.atom[Ar2].z;
		//cout << "x1 " << x1 << " y1 " << y1 << " z1 " << z1 << endl;
		//cout << "x2 " << x2 << " y2 " << y2 << " z2 " << z2 << endl;
		long double xo, yo, zo;
		xo=x2-x1; 
		yo=y2-y1;
		zo=z2-z1;
		//long double phio = acos(zo/rho);
		long double ro = sqrt(rho*rho-zo*zo);
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
		//cout << "xo: " << xo << " yo: " << yo << " zo: " << zo << endl;
		//cout << "u " << u << " v " << v << endl;
		//cout << "Phio: " << phio << " Thetao: " << thetao << endl;
        theta=thetao+(u-0.5)*dt;
		if (theta>2.0*pi)
	{
		theta -= 2.0*pi;
	}
		else if(theta<0.0){
			theta+=2.0*pi;
		}
//		long double phip = phio+dp;
//		long double phim = phio-dp;
//			if(phip<=pi && phim>=0.0){
//			long double cpm=cos(phim);
//			long double cpp=cos(phip);
//			phi= acos( cpm+v*(cpp-cpm) );

//		} 
//		else if (phip>pi)
//		{
//			long double cpm=cos(phim);
//			phi= acos( cpm + v*(-1.0-cpm) );
//		}
//		else if (phim<0.0)
//		{
//			long double cpp=cos(phip);
//			phi = acos( 1.0 + v*(cpp-1.0));
//		}
        // cout << "Phi: " << phi << " Theta: " << theta << endl;  
		  		
	
		long double xr, yr, rr;
	//	zr = rho*cos(phi);
		rr = rho;
		xr = rr*cos(theta);
		yr= rr*sin(theta);
		D.atom[Ar2].x = x1 + xr;
		D.atom[Ar2].y = y1 + yr;
//		D.atom[Ar2].z = z1 + zr;
		//cout << "xr: " << xr << " yr: " << yr << " zr: " << zr << " rr " << rr << endl; 
		//cout << "Ar2.x " << D.atom[Ar2].x << " Ar2.y " << D.atom[Ar2].y << " Ar2.z " << D.atom[Ar2].z << endl;
	//		cout << "----------------" << endl;
	
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

void CountContacts(long double& Npp, long double& Npn, long double& Nnn, Config& c){

	long double r;
	long double R1, R2, rc;
	unsigned long long N = c.ndimers*2;
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
					rc = (R1+R2)+1.3;
					
					r = RDist(c.dimer[Di1].atom[Ai1], c.dimer[Di2].atom[Ai2]);
					//cout<<"r1: "<<R1<<" r2: "<<R2<<" rc: "<<rc<<" r: "<<r<<endl;
					//Check proper distance
					if (r < rc){
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

	return;

} 

inline long double fastPow(long double a, long double b) {
	union {
		long double d;
		int x[2];
	} u = { a };
	u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
	u.x[0] = 0;
	return u.d;

}

//string itos( int Number ){
//	string ss = static_cast<ostringstream*>( &(ostringstream() << Number) )->str();
     
//     return ss;
// }

//string itos(int Number){
//	string s;
//	s = FormatInt(((Number & 1) ? Number : -Number)).c_str();
//	return(s);
//}

//long double inline __declspec (naked) __fastcall sqrt14(long double n)
//{
//	_asm fld qword ptr [esp+4]
//	_asm fsqrt
//	_asm ret 8
//} 
//long double inline __declspec attribute_(naked) __fastcall sqrt14(long double n)
//{
//	_asm fld qword ptr [esp+4]
//	_asm fsqrt
//	_asm ret 8
//} 







