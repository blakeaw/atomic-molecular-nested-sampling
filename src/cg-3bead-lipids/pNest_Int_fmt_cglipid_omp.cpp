/* Nested Sampling Routine--Collected nested energies.
Configured for 3 bead lipid model.

params and units:
sigma(distance) epsilon(energy)

 Using reduced units:
 d'=d/sigma  E'=E/epsilon


The Monte Carlo trial move consist of moving a lipid. There are three different trial moves.
	1) Center of mass translation of lipid molecule
	2) Rotation of the head bead or tail bead 2 about the central tail bead
	3) Translation of a single bead within a lipid
Each trial move has equal probability (~0.33) of being executed. 

Note: Has some multi-threading using OpenMP. To activate multi-threading compiling with g++ add the -fopenmp tag
To run without OpenMP multi-threading simply leave out the -fopenmp compiler tag when compiling.

Note: compiling with
	  c++ -O2 -ffast-math -fopenmp prog -o name
 and running with 10,000 lipids gives ~0.8% speed increase. Smaller numbers of lipids actually run slower
with current OMP pragmas.

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
// OpenMP Library
#ifdef _OPENMP
 #include <omp.h>
#endif
//RNG Fast Mersenne Twist- http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/
#include "/net/uu/nm/cm/bxw109120/src/MersenneTwist/dSFMT-src-2.2.1/dSFMT.c" 


using namespace std;

 // Classes -- cl
//#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/ConfigCGL.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/ConfigCGL.cpp"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/FeneBond.cpp"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/HarmonicBond.cpp"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/Bead.cpp"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/CGLipid.cpp"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/Histogram.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/HistogramNS.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Lipid/src/class/RunningStats.h"



// prototype functions -- pf

long double TotalPairEnergy(const ConfigCGL& c);
long double PairCGLipid(const Bead& one, const Bead& two, long double r);
long double TotalBondEnergy(const ConfigCGL& c);
long double DeltaEtype1(const ConfigCGL& c, const CGLipid& moved, const CGLipid& orig);
long double DeltaEtype2a3(const ConfigCGL& c, const CGLipid& moved, const CGLipid& orig, unsigned short bindex);
long double RDist(const Bead& one, const Bead& two);


void RandomRotation(const long double& u, const long double& v, const long double& stepa, CGLipid& L, const unsigned long& Br2);


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
	const unsigned short nlipids = 100;
	const unsigned nbeads = nlipids*3; 
	
	// the number of iterations to run on collecting new samples
	unsigned long int Iterations = 2;
	//Number of Markov style iterations to adjust configuration
	unsigned long long MCit = 200000;
	//Number of Trial Moves To Initially not collect data--allow for equilibration
	unsigned long Eqit = 1;
	//Nested Sampling Integration Fraction
	long double NSfrac = 0.50;
	
	// Used in place of pbc to adjust atom x,y coordinates
//	long double box = ((long double)(nlipids)*(rad1+rad2)*2.0)+rad1+rad2;
	long double box = 10000000.0;
	// Set for initial coordinates
	long double boxxx = box, boxyy = box, boxzz = box;
	
	long double pi = 3.14159265359;		
	//Movement Step size
	long double stepguess = 0.90*box;	
	long double step = stepguess;
	long double stepaguess = 1.0/12.0;
	long double stepa = stepaguess;
	//Energy Histogram Properties
	unsigned long long nbins = 6000000;
	long double Elow = (-2000.0)*(long double)(nlipids);
	long double Ehigh = 200.0;
	unsigned collectstepguess = nlipids*5;
	unsigned collectstep = collectstepguess;
	//Output Files
	string rundescript = "cube_l2_test_omp";
	//Name of logfile
	string logfilename = "oMMC_" + rundescript + ".log";
	//Name of Emedian output file
	string eoutname = "Epoints_int_" + rundescript + ".dat";
	//Name of Configuration output file
	string coutname = "Configs_int_" + rundescript + ".xyz";
	//Name of Observable output file
	string outfilename = "zOutput_int_" + rundescript + ".dat";
	
	
	//Name of NS values file for observables
	string nsoboutname = "zNSobs_int_" + rundescript + ".dat";
	cout << "Name of logfile is : " << logfilename.c_str() << endl;
	//note: ConfigsOutFreq should be large enough to collect only a few
	// configs from each nested sampling iteration
	unsigned int ConfigsOutFreq = (MCit-Eqit)/5; 

//----End of Parameter setting-----

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
	logfile << "Energy Histogram Bins: " << nbins << " Elow: " << Elow << " Ehigh: " << Ehigh << endl;
	logfile << "collectstepguess: " << collectstepguess << endl;
	logfile << "Number of Equilibration Steps: " << Eqit << endl;

	#pragma omp parallel sections
	{
		#pragma omp section
		{
		 	cout<<"There are "<<omp_get_num_threads()<<" threads"<<endl;
		}
	}
	cout << "System running with " << nlipids << " lipids and has a total of " << nbeads << " beads" << endl;
	logfile << "System running with " << nlipids << " lipids and has a total of " << nbeads << " beads" << endl;


	// Initialize the MT RNG object
	dsfmt_t dsfmt;	
	//Define a seed
	//time(0) is the clock value in seconds
	unsigned long long seed = time(0);
	
	// Initialize the RNG object with the seed
	
	dsfmt_init_gen_rand(&dsfmt, seed);

    


	//Create a Configuration Object
	ConfigCGL confn;
	confn.SetBox(box, box, box);
	confn.InitializeLipids(nlipids);
	ConfigCGL Bottom;
	Bottom.SetBox(box, box, box);
	Bottom.InitializeLipids(nlipids);
	
	
	//Assign an Initial Random Configuration 
	bool atag = 0;
	while(atag == 0){
		atag=1;
		//Assign coordinates	
		for (unsigned int jj = 0; jj < nlipids; ++jj)
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
//			cout << "assignind bead 0 to x: " << x << " y: " << y << " z: " << z << endl;
	//		cout << " z: " << z << endl;
			confn.lipid[jj].bead[0].x=x;
			confn.lipid[jj].bead[0].y=y;
			confn.lipid[jj].bead[0].z=z;
			long double u, v;
			u = dsfmt_genrand_close_open(&dsfmt);
			v = dsfmt_genrand_close_open(&dsfmt);		
			long double theta = u*2.0*pi;
			long double phi = acos((2.0*v) - 1.0);
			long double BondLength = confn.lipid[0].sigma;	
			long double dz = BondLength*cos(phi);
			long double r = sqrt(BondLength*BondLength-dz*dz);
			long double dx = r*cos(theta);
			long double dy = r*sin(theta);
		
			confn.lipid[jj].bead[1].x = x+dx;
			confn.lipid[jj].bead[1].y = y+dy;
			confn.lipid[jj].bead[1].z = z+dz;
			x +=dx;
			y+=dy;
			z+=dz;
//			cout << "assignind bead 1 to x: " << x << " y: " << y << " z: " << z << endl;
			u = dsfmt_genrand_close_open(&dsfmt);
			v = dsfmt_genrand_close_open(&dsfmt);		
			theta = u*2.0*pi;
			phi = acos((2.0*v) - 1.0);
			dz = BondLength*cos(phi);
			r = sqrt(BondLength*BondLength-dz*dz);
			dx = r*cos(theta);
			dy = r*sin(theta);
		
			confn.lipid[jj].bead[2].x = x+dx;
			confn.lipid[jj].bead[2].y = y+dy;
			confn.lipid[jj].bead[2].z = z+dz;
			x +=dx;
			y+=dy;
			z+=dz;
//			cout << "assignind bead 2 to x: " << x << " y: " << y << " z: " << z << endl;
		}
		
		confn.CalcCom();
		//Test boundary condition
		if(confn.BoxCheck() == 1){
			atag = 0;
		}
		long double Et = TotalPairEnergy(confn)+TotalBondEnergy(confn);
		//	exit(0);
		//cout << "Et: " << Et << endl;
			if(Et>Ehigh-0.0010 || Et<Elow+0.0010){
				atag=0;
			}		
	}
	//Assign the energy of the initial configuration
	long double Et1 = TotalPairEnergy(confn)+TotalBondEnergy(confn);
//	long double Et1 = TotalPairEnergy(confn);
	Bottom.EquateCoord(confn);
	cout << "Energy of Initial Configuration is: " << Et1 << endl;
	logfile << "Energy of Initial Configuration is: " << Et1 << endl;

	//Initialize the Energy Histogram	
	HistogramNS Ehist;
	Ehist.Initialize(nbins, Elow, Ehigh);
	Ehist.InitializeObservables(2);
	Ehist.SetNSfrac(NSfrac);


	
	// Average Energy and Average Energy^2
	RunningStats Estat;
	RunningStats Esqstat;
	

/*Begin iterating down phase space in 0.5 intervals collecting
MC trajectories and histogramming energies under each new median energy*/
cout << "Nested sampling iterations is set to " << Iterations << " with MC iterations set to " << MCit << endl;
logfile << "Nested sampling iterations is set to " << Iterations << " with MC iterations set to " << MCit << endl;
/*Define the first median energy--typically set
  to the upper limit of the Energy histogram */
long double Emedian = Ehigh;
//Used to track last accepted energy for histogram
long double Ee = Et1;



Et1 = TotalPairEnergy(confn)+TotalBondEnergy(confn);
cout << setiosflags(ios::fixed) << setprecision(15);
logfile << setiosflags(ios::fixed) << setprecision(15);

CGLipid OClip;

// Begin The Nested Sampling Iterations -- biter
for (unsigned int i = 0; i < Iterations; ++i){
	

	cout << "Iteration " << i << "......." << endl;
	logfile << "Iteration " << i << "......." << endl;

	

	if(i==1){

 		Ehist.Normalize();
 		Emedian = Ehist.Median();

		
	
	}
	if(i>1){
		Ehist.Renormalize();
		Emedian = Ehist.Median();
		
	}
 if(i>0){
 	cout << setiosflags(ios::fixed) << setprecision(15) <<"*** Emedian: " << Emedian << endl;
	logfile << setiosflags(ios::fixed) << setprecision(15) <<"*** Emedian: " << Emedian << endl;
	eoutfile << setiosflags(ios::fixed) << setprecision(15) << i << " " << Emedian << endl;

	//NS averages
	long double Eavgns, Esqavgns;
	#pragma omp parallel sections num_threads(2)
	{
		#pragma omp section
		{
			Eavgns = Ehist.GetObsAvg(0);
		}
		#pragma omp section
		{
			Esqavgns = Ehist.GetObsAvg(1);
		}
	}
	
	// Observables-- screen output
		//Trajectory averages
	cout << "* Trajectory Average Values:" << endl;
	cout << "Energy: " << Estat.Mean() << " and Energy Squared: " << Esqstat.Mean() << endl;
	
	
		//NS averages
	cout << "* NS Average Values (between NS energy at n and n-1):" << endl;
	cout << "Energy: " << Eavgns << " and Energy Squared: " << Esqavgns << endl;
	

	// Observables-- logfile output - same as screen output
	logfile << "* Trajectory Average Values:" << endl;
	logfile << "Energy: " << Estat.Mean() << " and Energy Squared: " << Esqstat.Mean() << endl;
	
	logfile << "* NS Average Values (between NS energy at n and n-1):" << endl;
	logfile << "Energy: " << Eavgns << " and Energy Squared: " << Esqavgns << endl;
	
	// Observables-- output to files
	// Trajectory Averages
	outfile << setiosflags(ios::fixed) << setprecision(15) << i << " " << Emedian << " " << Estat.Mean() << " " << Esqstat.Mean() <<  endl;
	// NS Averages
	nsobout << setiosflags(ios::fixed) << setprecision(15) << i << " " << Emedian << " " << Eavgns << " " << Esqavgns << " " << endl;

	// Reset Histogram and RunningStats
	#pragma omp parallel sections num_threads(2)
	{
		#pragma omp section 
		{
			Ehist.Reset();
		}
		#pragma omp section
		{
			Ehist.ResetObs();
		}
	//Set the new starting config
		#pragma omp section
		{
			confn.EquateCoord(Bottom);
		}
	}
	Estat.Reset();
	Esqstat.Reset();
	

	//Output the energy to screen and logfile
	cout << "* Energy of new initial configuration is " << Et1 << endl;
	logfile << "* Energy of new initial configuration is " << Et1 << endl;
}


	
	//Variable to track the Bottom Config
	long double Elowprev = 100000.0;
	//Used to track trial move success
	long double percent = 0.0;
	//long double Et1prev = Et1;
	//Begin the MC trajectory under current Emedian -- bmc
	for (unsigned long long int j = 0; j < MCit; ++j){
		
		//Save the Energy before trial move	
		long double Et1s = Et1;
		
		long double dE = 0.00000000000000;
	//Trial Moves

		// 1) randomly select an lipid index	
		unsigned int Lrandom = (unsigned int)(dsfmt_genrand_close_open(&dsfmt)*nlipids);
		//Assign a storage Lipid to save current coordinates
	
		OClip.EquateLipidCoord(confn.lipid[Lrandom]);

		// 2) pick the type of trial move

		// there 3 types of trial moves
		unsigned short ntrials = 3;
		// randomly chose one
		unsigned short ttype = (unsigned short)(ntrials*dsfmt_genrand_close_open(&dsfmt));
		
		

		// First type: ttype = 0; center of mass translation of lipid
		if (ttype == 0)
		{
			long double movex, movey, movez;
		/* set random movements of the coordinates
		   of the selected atom */
			movex = (dsfmt_genrand_close_open(&dsfmt)-0.5)*step;
			movey = (dsfmt_genrand_close_open(&dsfmt)-0.5)*step;
			movez = (dsfmt_genrand_close_open(&dsfmt)-0.5)*step*0.5;
			// Implement those moves
			unsigned short nb = confn.lipid[Lrandom].nbeads;
			for (unsigned int b = 0; b < nb; ++b)
			{
				confn.lipid[Lrandom].bead[b].x += movex;
				confn.lipid[Lrandom].bead[b].y += movey;
				confn.lipid[Lrandom].bead[b].z += movez;
			}
			dE = DeltaEtype1(confn, confn.lipid[Lrandom], OClip);
		}
		/* Second type: ttype = 1; random rotation
			of one of the connected beads
			about  the central bead */
		else if(ttype == 1){
			unsigned short ranrindex;
			long double ranp = dsfmt_genrand_close_open(&dsfmt);
			if (ranp<0.5)
			{
				ranrindex = 0;
			}
			else if (ranp>0.5)
			{
				ranrindex = 2;
			}
			long double u, v;
			u = dsfmt_genrand_close_open(&dsfmt);
			v = dsfmt_genrand_close_open(&dsfmt);
			
			RandomRotation(u, v, stepa, confn.lipid[Lrandom], ranrindex);
			dE = DeltaEtype2a3(confn, confn.lipid[Lrandom], OClip, ranrindex);
			
		}
		/* Third type: ttype = 2; random translation
			of one of the beads in the lipid */
		else if (ttype == 2)
		{
			long double movex, movey, movez;
		/* set random movements of the coordinates
		   of the selected atom */
			movex = (dsfmt_genrand_close_open(&dsfmt)-0.5)*step;
			movey = (dsfmt_genrand_close_open(&dsfmt)-0.5)*step;
			movez = (dsfmt_genrand_close_open(&dsfmt)-0.5)*step*0.5;
			// Implement those moves
			unsigned short rb = (unsigned short)(confn.lipid[Lrandom].nbeads*dsfmt_genrand_close_open(&dsfmt));
			
				confn.lipid[Lrandom].bead[rb].x += movex;
				confn.lipid[Lrandom].bead[rb].y += movey;
				confn.lipid[Lrandom].bead[rb].z += movez;
			
			dE =DeltaEtype2a3(confn, confn.lipid[Lrandom], OClip, rb);
			
		}
		
		
		
		
		
		
	
		//Selection Criteria-- energy less than Emedian is the Nested Sampling criterion
		
	
		bool tagg = 1;
		bool dccheck = confn.BoxCheck();
		
		
		
		bool hardsuccess = 0;
 		if(dccheck==0){
			hardsuccess = 1;
		
			//Add difference to current total potential energy
			Et1+=dE;
			
		}
	
		bool softsuccess = 0;
		
		if(Et1>Elow && Et1<Emedian){
			softsuccess = 1;
		}

		if (hardsuccess == 1 && softsuccess == 1){
		//	cout << "Trial move was accepted! " << endl;
			tagg=0;
			percent += 1.0;
		//	cout << "Et1 " << Et1 << endl;;
			
			//Ee=Et1;
		//	Ee = TotalPairEnergy(confn)+TotalWallEnergyZ(confn);
			Ee = Et1;
			
			if(Et1<Elowprev){
				Bottom.EquateCoord(confn);
				Elowprev = Et1;
				//cout << " Elowprev is now : " << Elowprev << endl;
				if(Elowprev > Emedian){
					cout << "!!!!Error: Elowprev is greater than Emedian!--- exiting..." << endl;
					logfile << "!!!!Error: Elowprev is greater than Emedian!--- exiting..." << endl;
					exit(0);
				}
				
			}
			if(abs(Ee)-abs(Et1)>0.00000000001){
				cout << "Ee and Et1 are no longer equal!!!!---exiting..." << endl;
				logfile << "Ee and Et1 are no longer equal!!!!---exiting..." << endl;
				exit(0);
			}
			
		}
		//cout << endl;
		/*Failed criteria--reset atom coordinates and remove
		  the difference in energy */ 
		if(tagg == 1){
			
			confn.lipid[Lrandom].EquateLipidCoord(OClip);
			if(hardsuccess == 1){
				Et1=Et1s;
			//	Et1 = TotalPairEnergy(confn)+TotalWallEnergyZ(confn);
			}
			//cout << "in fail dE: " << dE << " Et1: " << Et1 << endl;
		}
	//	cout << " Et1: " << Et1 << endl;
		dE = 0.0;
		if(j>Eqit && j%(collectstep)==0){
						
				if(Ee>Emedian){
					cout << "!!!Error: Pushed energy: " << Ee << " at MCit " << j << " is greater than Emed " << Emedian << endl;
					cout << "----Skipping the push!!" << endl;
					logfile << "!!!Error: Pushed energy: " << Ee << " at MCit " << j << " is greater than Emed " << Emedian << endl;
					logfile << "----Skipping the push!!" << endl;
					break;
				}
			//	cout << "Pushing energy: " << Ee << endl;
				//Add Energy to Histogram
				
				Ehist.Push(Ee);
				
			

				//Calculate Observables				
				

				//Push observables for Trajectory Averages
				Estat.Push(Ee);
				Esqstat.Push(Ee*Ee);
				
				//Push observables for NS Averages
				Ehist.PushObs(Ee, 0, Ee);
				Ehist.PushObs(Ee*Ee, 1, Ee);
			
		}
		
		/* Output coordinates from
  		 MC trajectory in an xyz file format */
		if(j>Eqit && j%ConfigsOutFreq ==0){
			coordsfile << nbeads << endl;
			coordsfile << "charged shell lj particles with Energy " << Ee << endl;
			for(unsigned long long hg = 0; hg<nbeads; hg++){
				unsigned long Di = hg/3;
				unsigned long long Ai = hg - (3*Di);
		
				coordsfile << setiosflags(ios::fixed) << setprecision(15) << confn.lipid[Di].bead[Ai].type << " " << confn.lipid[Di].bead[Ai].x << " " << confn.lipid[Di].bead[Ai].y << " " << confn.lipid[Di].bead[Ai].z << endl; 
			}
		}
		
	} //End of MCit Loop -- emc


	
	Et1 = Elowprev;
	
	
	coordsfile << nbeads << endl;
			coordsfile << "charged shell lj particles with Energy " << Et1 << endl;
			for(unsigned long long hg = 0; hg<nbeads; hg++){
				unsigned long Di = hg/3;
				unsigned long long Ai = hg - (3*Di);
	
				coordsfile << setiosflags(ios::fixed) << setprecision(15) << Bottom.lipid[Di].bead[Ai].type << " " << Bottom.lipid[Di].bead[Ai].x << " " << Bottom.lipid[Di].bead[Ai].y << " " << Bottom.lipid[Di].bead[Ai].z << endl; 
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

cout << "Percent: " << percent << " step: " << step << " stepa: " << stepa << endl;
logfile << "Percent: " << percent << " step: " << step << " stepa: " << stepa << endl;
if((unsigned)(Elowprev) == 100000){
		cout << "!!!*** Configuration is frozen. Exiting simulation. ***!!!" << endl;
		logfile << "!!!*** Configuration is frozen. Exiting simulation. ***!!!" << endl;
		break;
	}
if(percent<0.10){
	collectstep=(unsigned)(1.0/percent);
}

if(percent<0.000001){
	cout << "!!!*** Configuration is trapped in minima or otherwise prevented from sufficient sampling. ***!!!" << endl;
	cout << "!-----Ending Nested Sampling Iterations------!" << endl;
	logfile << "!!!*** Configuration is trapped in minima or otherwise prevented from sufficient sampling. ***!!!" << endl;
	logfile << "!-----Ending Nested Sampling Iterations------!" << endl;
	break;
}
cout << endl;
logfile << endl;	
//Approximate the Temperature
//Tempapprox = abs(Emedian-Elowprev)/2.0;
} //End of Iteration Loop -- eiter


eoutfile.close();
coordsfile.close();
outfile.close();
nsobout.close();

cout << "Output files are: " << endl;
cout << eoutname.c_str() << endl;
cout << coutname.c_str() << endl;
cout << outfilename.c_str() << endl;
cout << nsoboutname.c_str() << endl;

logfile << "Output files are: " << endl;
logfile << eoutname.c_str() << endl;
logfile << coutname.c_str() << endl;
logfile << outfilename.c_str() << endl;
logfile << nsoboutname.c_str() << endl;

cout << " " << endl;
cout << " Thank You and Have a Nice Day!" << endl;
cout << " " << endl;
cout << "... End of Line ..." << endl;


return EXIT_SUCCESS;

} // End of Main -- eom

// Define Functions -- df

//Total Energy- sum of all pair interactions
long double TotalPairEnergy(const ConfigCGL& c){
	long double Etot = 0.0;
	long double r;
	
	unsigned long long Nl = c.nlipids;
//	#pragma omp parallel shared(Etot)
	{
		long double Et = 0.0;
//		#pragma omp for
		for (unsigned int i = 0; i < Nl-1; ++i)
		{
			for (unsigned int j = i+1; j < Nl; ++j)
			{
				unsigned short Nb = c.lipid[i].nbeads;
				for (unsigned int b = 0; b < Nb; ++b)
				{
					unsigned short Nb2 = c.lipid[j].nbeads;
					for (unsigned int p = 0; p < Nb2; ++p)
					{
						r = RDist(c.lipid[i].bead[b], c.lipid[j].bead[p]);
						
						
							Et += PairCGLipid(c.lipid[i].bead[b], c.lipid[j].bead[p], r);
						
					}
				}
			}		
		}
//		#pragma omp atomic
		Etot += Et;
	}	
//	cout << "Total pair energy: " << Etot << endl;
	return(Etot);
}

long double PairCGLipid(const Bead& one, const Bead& two, long double r){
	long double pi=3.14159265358979323;
	long double E = 0.0;
	long double t1 = one.type;
	long double t2 = two.type;
	long double b, e;
	long double wc = 1.6*one.molec->sigma;
	long double rc;
	bool itag;

	if (t1==0 && t2==0)
	{
		b = one.molec->sigmahh;
		itag = 0;
	}
	if (t1==0 && t2==1)
	{
		b = one.molec->sigmaht;
		itag = 0;
	}
	if (t1==1 && t2==0)
	{
		b = one.molec->sigmaht;
		itag = 0;
	}
	if (t1==1 && t2==1)
	{
		b = one.molec->sigmatt;
		itag = 1;
	}
	
	rc = pow(2.0, 1.0/6.0)*b;
	e = one.molec->epsilon;
	if (r<=rc)
	{
		E += 4.0*e*(pow(b/r, 12) - pow(b/r, 6) + 0.25);
	}

	if (itag==1)
	{
		if (r<rc)
		{
			E += -e;
		}
		else if (r>= r && r<= rc+wc)
		{
			E += -e*cos((pi*(r-rc))/(2.0*wc))*cos((pi*(r-rc))/(2.0*wc));
		}
		else if (r> rc+wc){
			E += 0.0;
		}
	}
//	cout << "Inside pair CGLipid E: " << E << " r: " << r << " eps: " << e << " b " << b << " itag " << itag << " rc " << rc << " wc " << wc << endl;
	return(E);
	
}

long double TotalBondEnergy(const ConfigCGL& c){
	long double Et = 0.0;
	unsigned Nl = c.nlipids;
//	#pragma omp parallel shared(Et)
	{
		long double E = 0.0;
//		#pragma omp for
		for (unsigned int i = 0; i < Nl; ++i)
		{
			unsigned short nfb = c.lipid[i].nfbonds;
			for (unsigned int j = 0; j < nfb; ++j)
			{
				
				E += c.lipid[i].fbond[j].CalcPotential();
			}
//			cout << "Bond E after fene: " << E << " lipid " << i << endl;
			unsigned short nhb = c.lipid[i].nhbonds;
			for (unsigned int k = 0; k < nhb; ++k)
			{
				
				E += c.lipid[i].hbond[k].CalcPotential();
			}
//			cout << "Bond E after harmonic: " << E << " lipid " << i << endl;
		}
//		#pragma omp atomic
		Et += E;
	}
//	cout << "Total bond energy: " << E << endl;
	return(Et);

}

long double DeltaEtype1(const ConfigCGL& c, const CGLipid& moved, const CGLipid& orig){
	long double Etotp = 0.0;
	long double Etotop = 0.0;
	long double r;
	long double ro;
	unsigned mindex = moved.index;	
	
	unsigned long long Nl = c.nlipids;
	#pragma omp parallel shared(Etotp, Etotop)
	{
		long double Etot = 0.0;
		long double Etoto = 0.0;
	 	#pragma omp for
		for (unsigned int j = 0; j < Nl; ++j)
		{	
			if (j!=mindex)
			{
				unsigned short Nb = c.lipid[mindex].nbeads;
				for (unsigned int b = 0; b < Nb; ++b)
				{
					unsigned short Nb2 = c.lipid[j].nbeads;
					for (unsigned int p = 0; p < Nb2; ++p)
					{
						r = RDist(c.lipid[mindex].bead[b], c.lipid[j].bead[p]);	
						ro = RDist(orig.bead[b], c.lipid[j].bead[p]);
						
						{
							Etot += PairCGLipid(c.lipid[mindex].bead[b], c.lipid[j].bead[p], r);
							Etoto += PairCGLipid(c.lipid[mindex].bead[b], c.lipid[j].bead[p], ro);
						}
					}
				}
			}		
		}	
		#pragma omp atomic
		Etotp += Etot;
		#pragma omp atomic
		Etotop += Etoto;
	}
		
		long double dE = Etotp - Etotop;
		return(dE);
	
}

long double DeltaEtype2a3(const ConfigCGL& c, const CGLipid& moved, const CGLipid& orig, unsigned short bindex){
	long double Etotp = 0.0;
	long double Etotop = 0.0;
	
	long double r;
	long double ro;
	unsigned mindex = moved.index;	
	unsigned Nl = c.nlipids;
	#pragma omp parallel shared(Etotp, Etotop) 
	{
		long double Etot = 0.0;
		long double Etoto = 0.0;
		//Pair Interactions
		#pragma omp for
		for (unsigned int j = 0; j < Nl; ++j)
			{	
				if (j!=mindex)
				{
					
						unsigned short Nb2 = c.lipid[j].nbeads;
						for (unsigned int p = 0; p < Nb2; ++p)
						{
							r = RDist(c.lipid[mindex].bead[bindex], c.lipid[j].bead[p]);	
							ro = RDist(orig.bead[bindex], c.lipid[j].bead[p]);
							Etot += PairCGLipid(c.lipid[mindex].bead[bindex], c.lipid[j].bead[p], r);
							Etoto += PairCGLipid(c.lipid[mindex].bead[bindex], c.lipid[j].bead[p], ro);
						}
				}		
			}	
	
		//Bond Interactions
		unsigned nfb = moved.nfbonds;
		unsigned nhb = moved.nhbonds;
		//Fene bonds
		#pragma omp for
		for (unsigned int i = 0; i < nfb; ++i)
		{
			unsigned bi1 = moved.fbond[i].Bindex1;
			unsigned bi2 = moved.fbond[i].Bindex2;
			unsigned otb;
			if (bindex == bi1)
			{
				otb = bi2;
				r = RDist(c.lipid[mindex].bead[bindex], c.lipid[mindex].bead[otb]);	
				ro = RDist(orig.bead[bindex], c.lipid[mindex].bead[otb]);
				Etot += moved.fbond[i].CalcPotential(r);
				Etoto += moved.fbond[i].CalcPotential(ro);
			}
			else if(bindex == bi2){
				otb = bi1;
				r = RDist(c.lipid[mindex].bead[bindex], c.lipid[mindex].bead[otb]);	
				ro = RDist(orig.bead[bindex], c.lipid[mindex].bead[otb]);
				Etot += moved.fbond[i].CalcPotential(r);
				Etoto += moved.fbond[i].CalcPotential(ro);
			}
		
		}
		// Harmonic bonds
		#pragma omp for
		for (unsigned int i = 0; i < nhb; ++i)
		{
			unsigned bi1 = moved.hbond[i].Bindex1;
			unsigned bi2 = moved.hbond[i].Bindex2;
			unsigned otb;
			if (bindex == bi1)
			{
				otb = bi2;
				r = RDist(c.lipid[mindex].bead[bindex], c.lipid[mindex].bead[otb]);	
				ro = RDist(orig.bead[bindex], c.lipid[mindex].bead[otb]);
				Etot += moved.hbond[i].CalcPotential(r);
				Etoto += moved.hbond[i].CalcPotential(ro);
			}
			else if(bindex == bi2){
				otb = bi1;
				r = RDist(c.lipid[mindex].bead[bindex], c.lipid[mindex].bead[otb]);	
				ro = RDist(orig.bead[bindex], c.lipid[mindex].bead[otb]);
				Etot += moved.hbond[i].CalcPotential(r);
				Etoto += moved.hbond[i].CalcPotential(ro);
			}
		}
		#pragma omp atomic
		Etotp += Etot;
		#pragma omp atomic
		Etotop += Etoto;
	}//End omp parallel
	long double dE = Etotp - Etotop;
	return(dE);
}




//Calculates the distance between two atoms--wraps coordinates and gets minimum image
long double RDist(const Bead& one, const Bead& two){

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


void RandomRotation(const long double& u, const long double& v, const long double& stepa, CGLipid& L, const unsigned long& Br2){
		long double theta, phi;
        long double pi=3.14159265358979323;
		long double rho;
		long double dt, dp;
		
		unsigned long long Br1= 1;
		rho = RDist(L.bead[Br1], L.bead[Br2]);
		
		
		dt=(2.0*pi)*stepa;
		dp = pi*stepa;
		//cout << "rho " << rho << "stepa " << stepa << " dt " << dt << " dp " << dp << endl;
		long double x1, y1, z1, x2, y2, z2;
		x1 = L.bead[Br1].x;
		y1 = L.bead[Br1].y;
		z1 = L.bead[Br1].z;
		x2 = L.bead[Br2].x;
		y2 = L.bead[Br2].y;
		z2 = L.bead[Br2].z;
		//cout << "x1 " << x1 << " y1 " << y1 << " z1 " << z1 << endl;
		//cout << "x2 " << x2 << " y2 " << y2 << " z2 " << z2 << endl;
		long double xo, yo, zo;
		xo=x2-x1; 
		yo=y2-y1;
		zo=z2-z1;
		long double phio = acos(zo/rho);
		long double ro = sqrt(rho*rho-zo*zo);
		long double thetao=acos(xo/ro);
		if (yo<0.0 && xo<0.0)
		{
			thetao=-asin(yo/ro)+pi;
		}
		if(yo<0.0 && xo > 0.0){
			thetao = asin(yo/ro)+2.0*pi;
		}
		
		//cout << "xo: " << xo << " yo: " << yo << " zo: " << zo << endl;
		//cout << "u " << u << " v " << v << endl;
		//cout << "Phio: " << phio << " Thetao: " << thetao << endl;
        theta=thetao+(u-0.5)*dt;
		if (theta>2.0*pi)
	{
		theta -= 2.0*pi;
	}
		if(theta<0.0){
			theta+=2.0*pi;
		}
		long double phip = phio+dp;
		long double phim = phio-dp;
		if(phip<=pi && phim>=0.0){
			phi= acos( cos(phim)+v*(cos(phip)-cos(phim)) );

		} 
		else if (phip>pi)
		{
			phi= acos( cos(phim) + v*(-1.0-cos(phim)) );
		}
		else if (phim<0.0)
		{
			phi = acos( 1.0 + v*(cos(phip)-1.0));
		}
        // cout << "Phi: " << phi << " Theta: " << theta << endl;  
		  		
	
		long double xr, yr, zr, rr;
		zr = rho*cos(phi);
		rr = sqrt(rho*rho-zr*zr);
		xr = rr*cos(theta);
		yr= rr*sin(theta);
		L.bead[Br2].x = x1 + xr;
		L.bead[Br2].y = y1 + yr;
		L.bead[Br2].z = z1 + zr;
		//cout << "xr: " << xr << " yr: " << yr << " zr: " << zr << " rr " << rr << endl; 
		//cout << "Ar2.x " << D.bead[Ar2].x << " Ar2.y " << D.bead[Ar2].y << " Ar2.z " << D.bead[Ar2].z << endl;
	//		cout << "----------------" << endl;
	
		return;

}













