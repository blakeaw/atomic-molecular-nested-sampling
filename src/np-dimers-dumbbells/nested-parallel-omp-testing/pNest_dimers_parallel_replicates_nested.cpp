/* Nested Sampling Routine--Collected nested median energys.
Configured for NP dimer system with Radii r1, r2, mass m1, m2, and bondlength Lb, where r1<r2.
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
and a rotation move only 50% of the time.
The boundary condition is set to hard cubic box for particle coordinates (i.e. does not directly account for radius of particle)
Contains a simple routine to adjust step size under each new median energy.

Note: Uses OpenMP multi-threading to run replicates in parallel. 
Must add the -fopenmp flag when compiling with g++
Calls some omp functions that require the -fopenmp flag for compilation



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
//RNG Fast Mersenne Twist- http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/
#include "/net/uu/nm/cm/bxw109120/src/MersenneTwist/dSFMT-src-2.2.1/dSFMT.c" 

// Specify the number of dimers
#define n_dimers_n 10

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
long double PairSChargeShell(const Atom& one, const Atom& two, long double r, long double kappa);
long double PairLJColloid(const Atom& one, const Atom& two, long double r);
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
	const unsigned short ndimers = n_dimers_n;
	const unsigned natoms = ndimers*2; 
	
	// the number of iterations to run on collecting new samples
	const unsigned long int Iterations = 300;
	//Number of Markov style iterations to adjust configuration
	const unsigned long long MCit = 4000000;
	//Number of Trial Moves To Initially not collect data--allow for equilibration
	const unsigned long Eqit = 100;
	/* Number of replicate sims to run in parallel 
	   -- for single run (i.e. no parallel runs) set to 1 */
	const unsigned short numreplicates = 2;
	const unsigned short nprocs=4;
	//Nested Sampling Integration Fraction
	const long double NSfrac = 0.50;
	//Dimer Properties
	const long double rad1 = 1.0; // Positive charged
	const bool ctag1 = 1;
	const long double Ae1 = 1000.0;
	const long double mass1 = 1.0;
	string type1 = "P";
	const long double rad2 = 2.0; // Negative charged
	const bool ctag2 = 0;
	const long double Ae2 = 5000.0;
	const long double mass2 = 2.0;
	string type2 = "N";	
	const long double BondLength = rad1+rad2;
	// kappa value - inverse screening length	
	const long double kappa = 0.5000;
	// Used in place of pbc to adjust atom x,y coordinates
	const long double box = ((long double)(ndimers)*(rad1+rad2)*2.0)+rad1+rad2;
	// Set for initial coordinates
	long double boxxx = box, boxyy = box, boxzz = box/8.0;
	const long double WallZ = -box/2.0;
	const long double pi = 3.14159265359;		
	//Movement Step size
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
	
	const long double Elow = (-100.0)*(long double)(ndimers);
	const long double Ehigh = 200.0;
	const unsigned collectstepguess = ndimers*5;
	unsigned collectstep = collectstepguess;
	
	string rundescriptor = "cube_d2_test_mtest";
	
	//----End of Parameter setting-----


	//cout << "Name of logfile is : " << logfilename.c_str() << endl;
	//note: ConfigsOutFreq should be large enough to collect only a few
	// configs from each nested sampling iteration
	unsigned int ConfigsOutFreq = ((MCit/nprocs)-Eqit)/2; 


cout<<"Will run " <<numreplicates<<" replicate(s) of system:"<<endl;
cout <<"---- " <<ndimers<<" dimers and a total of " << natoms << " atoms" << endl;
omp_set_num_threads(numreplicates);
omp_set_nested(1);

#pragma omp parallel
{
	//Output Files
	string rundescript = rundescriptor;
	unsigned short th_id;
	th_id = omp_get_thread_num();
	
	string t_id = itos(th_id);
	rundescript+="_rep"+t_id;
	//Name of logfile
	string logfilename = "oMMC_" + rundescript + ".log";
	//Name of Emedian output file
	string eoutname = "Epoints_int_" + rundescript + ".dat";
	//Name of Configuration output file
	string coutname = "ConfigsDimer_int_" + rundescript + ".xyz";
	//Name of Observable output file
	string outfilename = "zOutput_int_" + rundescript + ".dat";
	//Name of Particle distribution file
	string gofzofname = "zGofz_int_" + rundescript + ".dat";
	//Name of NS values file for observables
	string nsoboutname = "zNSobs_int_" + rundescript + ".dat";
	#pragma omp critical(output)
	{
	cout << "Running Replicate with thread id: "<<th_id<<endl;
	cout << "Name of logfile is : " << logfilename.c_str() << endl;
	}
   //Initialize output file streams
	ofstream logfile(logfilename.c_str());
	ofstream eoutfile(eoutname.c_str());
	
	ofstream coordsfile(coutname.c_str());
	
	ofstream outfile(outfilename.c_str());
	ofstream nsobout(nsoboutname.c_str());	

	logfile << " " << endl;
	logfile << "Local time and date: " << asctime(timeinfo) << endl;
	logfile << "Parameters: " << endl;
	logfile << "Nested Sampling Energy Integration Fraction: " << NSfrac << endl;
	logfile << "R1: " << rad1 << " R2: " << rad2 << " ctag1: " << ctag1 << " ctag2 " << ctag1 << endl;
	logfile << "Ae1: " << Ae1 << " Ae2: " << Ae2 << " mass1: " << mass1 << " mass2: " << mass2 << endl;
	logfile << "kappa: " << kappa << endl;
	logfile << "box: " << box << " WallZ: " << WallZ << " stepguess: " << stepguess << " stepaguess: " << stepaguess << endl;
	logfile << "Energy Histogram Bins: " << nbins << " Elow: " << Elow << " Ehigh: " << Ehigh << endl;
	logfile << "collectstepguess: " << collectstepguess << endl;
	logfile << "Number of Equilibration Steps: " << Eqit << endl;
	logfile<<"ConfigsOutFreq: "<<ConfigsOutFreq<<endl;

	//
	logfile << "System running with " << ndimers << " dimers and has a total of " << natoms << " atoms" << endl;

	#pragma omp single
	{
		cout << "Running...."<<endl;
	}

	// Initialize the MT RNG object
	dsfmt_t dsfmt;	
	//Define a seed
	//time(0) is the clock value in seconds
	unsigned long long seed = time(0)+th_id;
	logfile<<"FMT Seed value: "<<seed<<endl;
	// Initialize the RNG object with the seed
	
	dsfmt_init_gen_rand(&dsfmt, seed);

    


	//Create a Configuration Object
	Config confn;
	confn.SetBox(box, box, box);
	confn.WallZ = WallZ;
	confn.SetKappa(kappa);
	Config Bottom;
	Bottom.SetBox(box, box, box);
	Bottom.WallZ = WallZ;
	Bottom.SetKappa(kappa);

	logfile << "Dimer Bond Length: " << confn.dimer[0].Lb << endl;
	//Assign dimer parameters
	for(unsigned int jj = 0;jj<ndimers;jj++){
			
			
			confn.dimer[jj].SetBond(BondLength);
			confn.dimer[jj].atom[0].type = type1;
			confn.dimer[jj].atom[0].rad = rad1;
			confn.dimer[jj].atom[0].mass = mass1;	
			confn.dimer[jj].atom[0].ctag = ctag1;	
			confn.dimer[jj].atom[0].Ae = Ae1;				
			confn.dimer[jj].atom[0].rsmall=rsmall;			
			confn.dimer[jj].atom[0].sigc=sigconv;	

			confn.dimer[jj].atom[1].type = type2;
			confn.dimer[jj].atom[1].rad = rad2;
			confn.dimer[jj].atom[1].mass = mass2;
			confn.dimer[jj].atom[1].ctag = ctag2;	
			confn.dimer[jj].atom[1].Ae = Ae2;	
			confn.dimer[jj].atom[1].rsmall=rsmall;
			confn.dimer[jj].atom[1].sigc=sigconv;
			//cout << "Assigned dimer number " << jj << endl;
	}
	//Assign an Initial Random Configuration 
	bool atag = 0;
	while(atag == 0){
		atag=1;
		//Assign coordinates	
		for (unsigned int jj = 0; jj < ndimers; ++jj)
		{
			long double ranx, rany, ranz;
			ranx = (dsfmt_genrand_close_open(&dsfmt)-0.5);
			rany = (dsfmt_genrand_close_open(&dsfmt)-0.5);
			ranz = (dsfmt_genrand_close_open(&dsfmt)-0.5);
	//		cout << " ranz: " << ranz << endl;
			long double x, y, z;
			x = ranx*(boxxx);
			y = rany*(boxyy);
			z = ranz*(boxzz);
	//		cout << " z: " << z << endl;
			confn.dimer[jj].atom[0].x=x;
			confn.dimer[jj].atom[0].y=y;
			confn.dimer[jj].atom[0].z=z;
			long double u, v;
			u = dsfmt_genrand_close_open(&dsfmt);
			v = dsfmt_genrand_close_open(&dsfmt);		
			long double theta = u*2.0*pi;
			long double phi = acos((2.0*v) - 1.0);
			long double BondLength = confn.dimer[jj].Lb;	
			long double dz = BondLength*cos(phi);
			long double r = sqrt(BondLength*BondLength-dz*dz);
			long double dx = r*cos(theta);
			long double dy = r*sin(theta);
		//	cout << "dx " << dx << " dy " << dy << " dz " << dz << endl;
			confn.dimer[jj].atom[1].x = x+dx;
			confn.dimer[jj].atom[1].y = y+dy;
			confn.dimer[jj].atom[1].z = z+dz;
			
		}
		
		confn.CalcCom();
		//Test boundary condition
		if(confn.BoxCheck() == 1){
			atag = 0;
		}
		long double Et = TotalPairEnergy(confn)+TotalWallEnergyZ(confn);
		//cout << "Et: " << Et << endl;
			if(Et>Ehigh-0.0010 || Et<Elow+0.0010){
				atag=0;
			}		
	}
	//Assign the energy of the initial configuration
	long double Et1 = TotalPairEnergy(confn)+TotalWallEnergyZ(confn);
//	long double Et1 = TotalPairEnergy(confn);
	Bottom.Equate(confn);
	//cout << "Energy of Initial Configuration is: " << Et1 << endl;
	logfile << "Energy of Initial Configuration is: " << Et1 << endl;
	//Initialize the Energy Histogram	
	HistogramNS Ehist;
	Ehist.Initialize(nbins, Elow, Ehigh);
	Ehist.InitializeObservables(7);
	Ehist.SetNSfrac(NSfrac);
	//Track the random angles statistics
	//RunningStats ThetaStat;
	//RunningStats PhiStat;
	//Histogram Thist;
	//Thist.Initialize(10000, 0.0, 360.0);
	//Histogram Phist;
	//Phist.Initialize(10000, 0.0, 180.0);

	//Track the radius of gyration and gyration density
	RunningStats Rgstat;
	RunningStats Dstat;
	// Average Energy and Average Energy^2
	RunningStats Estat;
	RunningStats Esqstat;
	//Z-wall particle distribution
	Histogram Gofz;
	Gofz.Initialize(100000, WallZ, box/2.0 );
	//Orientation OPs
	RunningStats Dorstat;
	RunningStats Zorstat;
	//Zwall particel distance
	RunningStats Zdstat;
	
	RunningStats dEstat;

/*Begin iterating down phase space in 0.5 intervals collecting
MC trajectories and histogramming energies under each new median energy*/
//cout << "Nested sampling iterations is set to " << Iterations << " with MC iterations set to " << MCit << endl;
logfile << "Nested sampling iterations is set to " << Iterations << " with MC iterations set to " << MCit << endl;
/*Define the first median energy--typically set
  to the upper limit of the Energy histogram */
long double Emedian = Ehigh;
//Used to track last accepted energy for histogram
long double Ee = Et1;
long double dEe = 0.0;
long double dXe = 0.0;
long double Xe = 0.0;

Et1 = TotalPairEnergy(confn)+TotalWallEnergyZ(confn);
//cout << setiosflags(ios::fixed) << setprecision(15);
logfile << setiosflags(ios::fixed) << setprecision(15);

long double Tempapprox;
// Begin The Nested Sampling Iterations -- biter
for (unsigned int i = 0; i < Iterations; ++i){
	

	//cout << "Iteration " << i << "......." << endl;
	logfile << "Iteration " << i << "......." << endl;

	dEe = 0.0, dXe = 0.0, Xe = 0.0;
	long double Emednext;

 if(i>0){
	Emedian = Ehist.NestedEnergy();
//	Emednext = Ehist.PredictNextNestedEnergy();
 	//cout << setiosflags(ios::fixed) << setprecision(15) <<"*** Emedian: " << Emedian << endl;
	logfile << setiosflags(ios::fixed) << setprecision(15) <<"*** Emedian: " << Emedian << endl;
	eoutfile << setiosflags(ios::fixed) << setprecision(15) << i << " " << Emedian << endl;
//	//cout << "The next Nested Energy is predicted to be: " << Emednext << endl;
//cout << "Average energy change: " << dEstat.Mean() << endl;	
//	cout << "Approximate Temperature is: " << Tempapprox << endl;	

	//NS averages
	long double Eavgns = Ehist.GetObsAvg(0);
	long double Esqavgns = Ehist.GetObsAvg(1);
	long double Rgns = Ehist.GetObsAvg(2);
	long double Drgns = Ehist.GetObsAvg(3);
	long double Zorns = Ehist.GetObsAvg(4);
	long double Dorns = Ehist.GetObsAvg(5);
	long double Zdns = Ehist.GetObsAvg(6);

	// Observables-- screen output
		//Trajectory averages
	//cout << "* Trajectory Average Values:" << endl;
	//cout << "Energy: " << Estat.Mean() << " and Energy Squared: " << Esqstat.Mean() << endl;
	//cout << "Radius of gyration: " << Rgstat.Mean() << endl;
	//cout << "Gyration Density: " << Dstat.Mean() << endl;
	//cout << "Z-orientation Parameter: " << Zorstat.Mean() << endl;
	//cout << "Dimer-Dimer(NP-NP) Orientation Parameter " << Dorstat.Mean() << endl;
	//cout << "Zwall distance: " << Zdstat.Mean() << endl;
	
		//NS averages
	//cout << "* NS Average Values (between NS energy at n and n-1):" << endl;
	//cout << "Energy: " << Eavgns << " and Energy Squared: " << Esqavgns << endl;
	//cout << "Radius of Gyration: " << Rgns << endl;
	//cout << "Gyration Density: " << Drgns << endl;
	//cout << "Z-orientation Parameter: " << Zorns << endl;
	//cout << "Dimer-Dimer(NP-NP) Orientation Parameter: " << Dorns << endl;
	//cout << "Zwall distance: " << Zdns << endl;

	// Observables-- logfile output - same as screen output
	logfile << "* Trajectory Average Values:" << endl;
	logfile << "Energy: " << Estat.Mean() << " and Energy Squared: " << Esqstat.Mean() << endl;
	logfile << "Radius of gyration: " << Rgstat.Mean() << endl;
	logfile << "Gyration Density: " << Dstat.Mean() << endl;
	logfile << "Z-orientation Parameter: " << Zorstat.Mean() << endl;
	logfile << "Dimer-Dimer(NP-NP) Orientation Parameter " << Dorstat.Mean() << endl;
	logfile << "Zwall distance: " << Zdstat.Mean() << endl;
	logfile << "* NS Average Values (between NS energy at n and n-1):" << endl;
	logfile << "Energy: " << Eavgns << " and Energy Squared: " << Esqavgns << endl;
	logfile << "Radius of Gyration: " << Rgns << endl;
	logfile << "Gyration Density: " << Drgns << endl;
	logfile << "Z-orientation Parameter: " << Zorns << endl;
	logfile << "Dimer-Dimer(NP-NP) Orientation Parameter: " << Dorns << endl;
	logfile << "Zwall distance: " << Zdns << endl;
	// Observables-- output to files
	// Trajectory Averages
	outfile << setiosflags(ios::fixed) << setprecision(15) << i << " " << Emedian << " " << Estat.Mean() << " " << Esqstat.Mean() << " ";
	//outfile << Rgstat.Mean() << " " << Dstat.Mean() << " " << Zorstat.Mean() << " " << Dorstat.Mean() << endl;
	outfile << Rgstat.Mean() << " " << Dstat.Mean() << " " << Zorstat.Mean() << " " << Dorstat.Mean() << " " << Zdstat.Mean() << endl;
	// NS Averages
	nsobout << setiosflags(ios::fixed) << setprecision(15) << i << " " << Emedian << " " << Eavgns << " " << Esqavgns << " ";
	nsobout << Rgns << " " << Drgns << " " << Zorns << " " << Dorns << " " << Zdns << endl;
	// Reset Histogram and RunningStats
	Ehist.SmartReset();
//	Ehist.ResetObs();
	Rgstat.Reset();
	Dstat.Reset();
	Zorstat.Reset();
	Dorstat.Reset();
	Estat.Reset();
	Esqstat.Reset();
	Zdstat.Reset();
	dEstat.Reset();
	//Set the new starting config
	confn.Equate(Bottom);
	//Output the energy to screen and logfile
	//cout << "* Energy of new initial configuration is " << Et1 << endl;
	logfile << "* Energy of new initial configuration is " << Et1 << endl;
}


	//if (i==Iterations)
	//{
		//break;
	//}
	//Variable to track the Bottom Config
	long double Elowprev = 100000.0;
	//Used to track trial move success
	long double percent = 0.0;
	long double Et1prev = Et1;
	unsigned long long MCitpp = MCit/nprocs;
//	unsigned long long MCitpp = MCit;
	omp_set_num_threads(nprocs);
	#pragma omp parallel
	{
	Config Tlocal(confn);
	//Tlocal.SetBox(box, box, box);
	//Tlocal.Equate(confn);
	long double Eetloc;
	long double Etloc=Et1;	
	long double percentlocal=0.0;
	
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
		long double origx, origy, origz, origr;
		origx = Tlocal.dimer[Dr].atom[Ar].x;
		origy = Tlocal.dimer[Dr].atom[Ar].y;
		origz = Tlocal.dimer[Dr].atom[Ar].z;
		Tlocal.dimer[Dr].CalcCom();
		long double ocomx, ocomy, ocomz;
		ocomx = Tlocal.dimer[Dr].COM[0];
		ocomy = Tlocal.dimer[Dr].COM[1];
		ocomz = Tlocal.dimer[Dr].COM[2];
		
		long double movex, movey, movez;
		/* set random movements of the coordinates
		   of the selected atom */
		movex = (dsfmt_genrand_close_open(&dsfmt)-0.5)*step;
		movey = (dsfmt_genrand_close_open(&dsfmt)-0.5)*step;
		movez = (dsfmt_genrand_close_open(&dsfmt)-0.5)*step*0.5;
		// Implement those moves
		Tlocal.dimer[Dr].atom[Ar].x += movex;
		Tlocal.dimer[Dr].atom[Ar].y += movey;
		Tlocal.dimer[Dr].atom[Ar].z += movez;
		
		//Now apply a move to the other atom of the dimer
		unsigned long long Ar2;
		if (Ar==0)
		{
			Ar2=1;
		}
		else{Tlocal.dimer[Dr].CalcCom();
			Ar2=0;
		}
		//Store original coordinates
		long double origx2, origy2, origz2, origr2;
		origx2 = Tlocal.dimer[Dr].atom[Ar2].x;
		origy2 = Tlocal.dimer[Dr].atom[Ar2].y;
		origz2 = Tlocal.dimer[Dr].atom[Ar2].z;
		
				

		long double ranprob = dsfmt_genrand_close_open(&dsfmt);
		//Only rotate 50% of trial moves
		
		if(ranprob<=0.25){
			//Pick random x, y, and z values
			long double u, v;
			u = dsfmt_genrand_close_open(&dsfmt);
			v = dsfmt_genrand_close_open(&dsfmt);
			Tlocal.dimer[Dr].atom[Ar2].x += movex;
			Tlocal.dimer[Dr].atom[Ar2].y += movey;
			Tlocal.dimer[Dr].atom[Ar2].z += movez;
			RandomRotation(u, v, stepa, Tlocal.dimer[Dr], Ar);
		}
		else{
			Tlocal.dimer[Dr].atom[Ar2].x += movex;
			Tlocal.dimer[Dr].atom[Ar2].y += movey;
			Tlocal.dimer[Dr].atom[Ar2].z += movez;

		}
		long double nz1, nz2;
		nz1=Tlocal.dimer[Dr].atom[Ar].z;
		nz2=Tlocal.dimer[Dr].atom[Ar2].z;
		long double r1, r2;
		r1 = Tlocal.dimer[Dr].atom[Ar].rad;
		r2 = Tlocal.dimer[Dr].atom[Ar2].rad;
		// Calculate the difference in total potential energy
		//cout << "Emedian: " << Emedian << endl;
		Tlocal.dimer[Dr].CalcCom();
		long double ncomx, ncomy, ncomz;
		ncomx = Tlocal.dimer[Dr].COM[0];
		ncomy = Tlocal.dimer[Dr].COM[1];
		ncomz = Tlocal.dimer[Dr].COM[2];
		long double dXx, dXy, dXz;
		dXx = ncomx - ocomx;
		dXy = ncomy - ocomy;
		dXz = ncomz - ocomz;
		long double Xn = sqrt(ncomx*ncomx+ncomy*ncomy+ncomz*ncomz);
		
		long double dE = 0.00000000000000;
	//	long double dX = sqrt(dXx*dXx+dXy*dXy+dXz*dXz);
		long double dX = Xn - sqrt(ocomx*ocomx+ocomy*ocomy+ocomz*ocomz);
	//cout << "Et1 is " << Et1 << " and dE is " << dE << endl;
		//long double Ea = TotalPairEnergy(confn)+TotalWallEnergyZ(confn);
		//Selection Criteria-- energy less than Emedian is the Nested Sampling criterion
		
	//	cout << "Dist is " << dist << endl;
		bool tagg = 1;
		bool dccheck = Tlocal.BoxCheck(Dr);
		
		bool wallcheck = 0;
		if(nz1<=WallZ+r1+0.100 || nz2 <= WallZ+r2+0.100){
			//dE += 1000000.0;
			wallcheck = 1;
		//	Ea += 1000000.0;
		//	cout << "Wall Fail!!!" << endl;
		}
		
		bool hardsuccess = 0;
 		if(wallcheck == 0 && dccheck==0){
			hardsuccess = 1;
			
			long double dEpp = DeltaE(Tlocal, Dr, Ar, origx, origy, origz, origx2, origy2, origz2);
	
			long double dEww = dEWallZ(Tlocal.dimer[Dr], origz, r1, origz2, r2, WallZ);
		
			dE = dEpp + dEww;
			//Check for nan--if nan set dE high to throw out trial move
			if(dE!=dE){
			//	cout<<"---Error-- dE is: "<<dE<<endl;
				dE=1000000.0;
			}
			//Add difference to current total potential energy
			Etloc+=dE;
			
		}
	
	//	unsigned short th_idb;
	//	th_idb = omp_get_thread_num();
		//cout<<"Thread id: "<<th_idb << " dE is now " << dE << " and Etloc plus dE is " << Etloc << endl;
	//	cout << endl;
		bool softsuccess = 0;
		
		if(Etloc>Elow && Etloc<Emedian){
			softsuccess = 1;
		}
		#pragma omp critical
		{
			
		if (hardsuccess == 1 && softsuccess == 1){
		//	cout << "Trial move was accepted! " << endl;
			tagg=0;
			percentlocal += 1.0;
		//	cout << "Et1 " << Et1 << endl;;
			
			//Ee=Et1;
		//	Ee = TotalPairEnergy(confn)+TotalWallEnergyZ(confn);
			Eetloc = Etloc;
			dEe = dE;	
			dXe = dX;	
			Xe = Xn;	
		//	cout << "Et1: " << Et1 << " Ee: " << Ee << endl;
		
			if(Etloc<Elowprev){
				Bottom.Equate(Tlocal);
				Elowprev = Etloc;
				//cout << " Elowprev is now : " << Elowprev << endl;
				if(Elowprev > Emedian){
					//cout << "!!!!Error: Elowprev is greater than Emedian!--- exiting..." << endl;
					logfile << "!!!!Error: Elowprev is greater than Emedian!--- exiting..." << endl;
					exit(0);
				}
				
			}
			if(abs(Eetloc)-abs(Etloc)>0.00000000001){
				//cout << "Ee and Et1 are no longer equal!!!!---exiting..." << endl;
				logfile << "Ee and Et1 are no longer equal!!!!---exiting..." << endl;
				exit(0);
			}
			
		}
		//cout << endl;
		/*Failed criteria--reset atom coordinates and remove
		  the difference in energy */ 
		if(tagg == 1){
			
			Tlocal.dimer[Dr].atom[Ar].x = origx;
			Tlocal.dimer[Dr].atom[Ar].y = origy;
			Tlocal.dimer[Dr].atom[Ar].z = origz;
			Tlocal.dimer[Dr].atom[Ar2].x = origx2;
			Tlocal.dimer[Dr].atom[Ar2].y = origy2;
			Tlocal.dimer[Dr].atom[Ar2].z = origz2;
			if(hardsuccess == 1){
				Etloc=Et1s;
			//	Et1 = TotalPairEnergy(confn)+TotalWallEnergyZ(confn);
			}
			//cout << "in fail dE: " << dE << " Et1: " << Et1 << endl;
		}
	//	cout << " Et1: " << Et1 << endl;
		dE = 0.0;
		
		if(j>Eqit && j%(collectstep)==0){
			//#pragma omp critical
			 // {			
				//if(Eetloc>Emedian){
					//cout << "!!!Error: Pushed energy: " << Ee << " at MCit " << j << " is greater than Emed " << Emedian << endl;
					//cout << "----Skipping the push!!" << endl;
				//	logfile << "!!!Error: Pushed energy: " << Eetloc << " at MCit " << j << " is greater than Emed " << Emedian << endl;
				//	logfile << "----Skipping the push!!" << endl;
					//break;
			//	}
				//cout << "Thread id: "<<th_idb<<" Pushing energy: " << Eetloc << endl;
				//Add Energy to Histogram
				
				Ehist.Push(Eetloc);
				
			//	long double xdedx = Xe*(abs(dEe)/dXe);
				long double xdedx = abs(Xe*(dEe/dXe));

				//Calculate Observables				
				long double Rg = RadGyr(Tlocal);
				long double dynVol = (4.0/3.0)*pi*(Rg*Rg*Rg);
				long double dynDens = ((long double)(Tlocal.ndimers)*2.0)/dynVol;
				long double Zor = AverageZorientation(Tlocal);
				long double Dor= NPvecProd(Tlocal);
				long double Zdist = AverageZwallDist(Tlocal);

				//Push observables for Trajectory Averages
				Estat.Push(Eetloc);
				Esqstat.Push(Eetloc*Eetloc);
				Rgstat.Push(Rg);
				Dstat.Push(dynDens);
				Zdensity(Tlocal, Gofz);
				Zorstat.Push(Zor);
				Dorstat.Push(Dor);
				Zdstat.Push(Zdist);
				//Push observables for NS Averages
				Ehist.PushObs(Eetloc, 0, Eetloc);
				Ehist.PushObs(Eetloc*Eetloc, 1, Eetloc);
				Ehist.PushObs(Rg, 2, Eetloc);
				Ehist.PushObs(dynDens, 3, Eetloc);
				Ehist.PushObs(Zor, 4, Eetloc);
				Ehist.PushObs(Dor, 5, Eetloc);
				Ehist.PushObs(Zdist, 6, Eetloc);

				dEstat.Push(xdedx);
			//}//End Critical	
			
		}
		
		/* Output coordinates from
  		 MC trajectory in an xyz file format */
		if(j>Eqit && j%ConfigsOutFreq ==0){
			
			coordsfile << natoms << endl;
			coordsfile << "charged shell lj particles with Energy " << Etloc << endl;
			for(unsigned long long hg = 0; hg<natoms; hg++){
				unsigned long Di = hg/2;
				unsigned long long Ai = hg - (2*Di);
		
				coordsfile << setiosflags(ios::fixed) << setprecision(15) << Tlocal.dimer[Di].atom[Ai].type.c_str() << " " << Tlocal.dimer[Di].atom[Ai].x << " " << Tlocal.dimer[Di].atom[Ai].y << " " << Tlocal.dimer[Di].atom[Ai].z << endl; 
			}
		 //End Critical
		}

	  }//End Critical

	} //End of MCit Loop -- emc
	#pragma omp barrier

	#pragma omp atomic
	 percent+=percentlocal;

}//End Nested parallel
	
	Et1 = Elowprev;
	
	
	coordsfile << natoms << endl;
			coordsfile << "charged shell lj particles with Energy " << Et1 << endl;
			for(unsigned long long hg = 0; hg<natoms; hg++){
				unsigned long Di = hg/2;
				unsigned long long Ai = hg - (2*Di);
	
				coordsfile << setiosflags(ios::fixed) << setprecision(15) << Bottom.dimer[Di].atom[Ai].type.c_str() << " " << Bottom.dimer[Di].atom[Ai].x << " " << Bottom.dimer[Di].atom[Ai].y << " " << Bottom.dimer[Di].atom[Ai].z << endl; 
			}

	
	
//Calc percent
percent /= (long double)(MCit)/100.0;

//Adjust step size if needed
if (percent < 45.0)
{
	step *= 0.80;
	stepa *= 0.90;
	if(percent<20.0&&percent>10.0){
		step *= 0.50;
		stepa *= 0.80;
	}
	if(percent<10.0){
		step *= 0.25;
		stepa *= 0.75;
	}
}
if (step < 0.01)
{
	step = stepguess/5.0;
}
if (stepa < stepaguess/100.0){
	stepa = stepaguess/5.0;
}

if (percent>65.0)
{
	step *= 1.5;
	stepa *= 1.5;
}

//cout << "Percent: " << percent << " step: " << step << " stepa: " << stepa << endl;
logfile << "Percent: " << percent << " step: " << step << " stepa: " << stepa << endl;
if((unsigned)(Elowprev) == 100000){
		//cout << "!!!*** Configuration is frozen. Exiting simulation. ***!!!" << endl;
		logfile << "!!!*** Configuration is frozen. Exiting simulation. ***!!!" << endl;
		break;
	}
if(percent<0.10){
	collectstep=(unsigned)(1.0/percent);
}

if(percent<0.000001){
	//cout << "!!!*** Configuration is trapped in minima or otherwise prevented from sufficient sampling. ***!!!" << endl;
	//cout << "!-----Ending Nested Sampling Iterations------!" << endl;
	logfile << "!!!*** Configuration is trapped in minima or otherwise prevented from sufficient sampling. ***!!!" << endl;
	logfile << "!-----Ending Nested Sampling Iterations------!" << endl;
	break;
}
//cout << endl;
logfile << endl;	
//Approximate the Temperature
//Tempapprox = abs(Emedian-Elowprev)/2.0;
} //End of Iteration Loop -- eiter

Gofz.Renormalize();
Gofz.DumpNormal(gofzofname.c_str());
eoutfile.close();
coordsfile.close();
outfile.close();
nsobout.close();
#pragma omp critical(filenames)
{
cout << "Rep "<<th_id<<":"<<endl;
cout << "Output files are: " << endl;
cout << eoutname.c_str() << endl;
cout << coutname.c_str() << endl;
cout << outfilename.c_str() << endl;
cout << gofzofname.c_str() << endl;
cout << nsoboutname.c_str() << endl;
}
logfile << "Output files are: " << endl;
logfile << eoutname.c_str() << endl;
logfile << coutname.c_str() << endl;
logfile << outfilename.c_str() << endl;
logfile << gofzofname.c_str() << endl;
logfile << nsoboutname.c_str() << endl;

} // End of parallel
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
	return( dE ); 
	

}

long double PairSChargeShell(const Atom& one, const Atom& two, long double r, long double kappa){
	
	long double pi = 3.14159265359;
	long double potent;
	long double R1, R2, Ae1, Ae2, Ae;
	bool ctag1, ctag2;
//	long double k4 = pow(kappa, 4);
	ctag1 = (bool)(one.ctag), ctag2 = (bool)(two.ctag);
	R1 = one.rad, R2 = two.rad;
	Ae1 = one.Ae, Ae2 = two.Ae;
//	kappa =1.00;
	if(ctag1 == ctag2){
		//Ae = sqrt(Ae1*Ae2);
		Ae = Ae1;
	}
	else{
		Ae = -sqrt(Ae1*Ae2);
	}

	if(r<R1+R2){
		potent = 100000.0;
	}
	else{
	potent = (Ae/(r))*(exp(-kappa*(R1+R2+r)))*(-1.0+exp(-2.0*kappa*R1))*(-1.0+exp(-2.0*kappa*R2));
	}

	//	cout<<"R1: "<<R1<<" R2: "<<R2<<endl;
		//cout<<"ctag1: "<<one.ctag<<" ctag2: "<<two.ctag<<endl;
	//	cout<<"ctag1: "<<ctag1<<" ctag2: "<<ctag2<<endl;
	//	cout<<"Ae1: "<<Ae1<<" Ae2: "<<Ae2<<endl;
	//	cout << "Ae "<<Ae<<" kappa "<<kappa<<endl;
	//	cout << "Screened Pair Potential is: " << potent << " at distance: " << r << endl;
	//	cout<<endl;
//	cout << "PSCS potent: " << potent << endl;
//	if(potent != potent){
	//	cout << "R1: " << R1 << " R2: " << R2 << " r: " << r << " Ae: " << Ae << endl;
		//exit(0);
	//}
	
	
	return(potent);
}

long double PairLJColloid(const Atom& one, const Atom& two, long double r){
	
	long double pi = 3.14159265359;
	long double potent, Pa, Pb;
	long double R1s, R2s;
	long double Acc, F;
	long double R1, R2;
	R1 = one.rad, R2 = two.rad;
	//Reduced Units
	Acc = 1.0;
	F = 50.0;
	R1s = pow(R1, 2);
	R2s = pow(R2, 2);
	
	
	if(r<R1+R2+0.01){
		
		potent = 100000.0;
	}
	else{
		Pa = -(Acc/6.0) * ( ( (2.0*R1*R2)/(pow(r, 2)-pow(R1+R2, 2)) ) + ( (2.0*R1*R2)/(pow(r, 2)-pow(R1-R2, 2)) ) + log((pow(r, 2)-pow(R1+R2, 2))/(pow(r, 2)-pow(R1-R2, 2))) );

		Pb = ((Acc)/(37800.0*r*pow(F, 6))) * ( ( (pow(r, 2)-7.0*r*(R1+R2)+6.0*(R1s+7.0*R1*R2+R2s))/(pow(r-R1-R2, 7))) + ( (pow(r, 2)+7.0*r*(R1+R2)+6.0*(R1s+7.0*R1*R2+R2s))/(pow(r+R1+R2, 7))) - ( (pow(r, 2)+7.0*r*(R1-R2)+6.0*(R1s-7.0*R1*R2+R2s))/(pow(r+R1-R2, 7))) - ( (pow(r, 2)-7.0*r*(R1-R2)+6.0*(R1s-7.0*R1*R2+R2s))/(pow(r-R1+R2, 7))) ) ;
		potent = Pa+Pb;
	}
//	cout << "PC potent: " << potent << endl;
	//if(potent != potent){
		//exit(0);
	//}
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
//Calculates the distance between an atom and the old coordinates of another--wraps coordinates and gets minimum image
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

long double LjWallZp(const long double& oz, const long double& rad, const long double& WallZ){

	long double Zwall, dz;
	//long double pi = 3.14159265359;
	long double zc = oz;
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

	long double Eo=0.0, En=0.0, zn;
	unsigned short N = d.natoms;
	
	
	
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
		long double phio = acos(zo/rho);
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
		long double phip = phio+dp;
		long double phim = phio-dp;
			if(phip<=pi && phim>=0.0){
			long double cpm=cos(phim);
			long double cpp=cos(phip);
			phi= acos( cpm+v*(cpp-cpm) );

		} 
		else if (phip>pi)
		{
			long double cpm=cos(phim);
			phi= acos( cpm + v*(-1.0-cpm) );
		}
		else if (phim<0.0)
		{
			long double cpp=cos(phip);
			phi = acos( 1.0 + v*(cpp-1.0));
		}
        // cout << "Phi: " << phi << " Theta: " << theta << endl;  
		  		
	
		long double xr, yr, zr, rr;
		zr = rho*cos(phi);
		rr = sqrt(rho*rho-zr*zr);
		xr = rr*cos(theta);
		yr= rr*sin(theta);
		D.atom[Ar2].x = x1 + xr;
		D.atom[Ar2].y = y1 + yr;
		D.atom[Ar2].z = z1 + zr;
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











