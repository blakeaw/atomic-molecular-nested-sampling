/* Nested Sampling Routine--Collected nested enthapies.
Configured to run for lj particles with sigma=1.0 and epsilon=1.0.
The Monte Carlo trial move consist of moving a single atom and volume deformations
 The boundary condition is set to a cubic box ; 
---Contains conditions to modulate step sizes under each new median enthalpy.

Modifications
3-28-14
	-Added in OpenMP Multithreading to break sampling loop over multiple procs/threads
	Updated from v2 to v3.

3-28-14
	-Switched from spherical COM box condition to cubic box - so does uniform xyz box size 
 		volume deformations (i.e. maintains a cubic box during volume deformation)
	-Coordinates are wrapped into box at trial move stage
	-Added minimum image periodicity in pair energy funtions

4-2-14
	-Updated volume deformation trial moves to compute the volume deformation probability (vprob)
		before it does anything else and then use continue statement to move to next trial move
		if the deformation is rejected by that criterion
	Updated from v4 to v4.1 -- runs a little faster. speed increase is magnified with more volume 
		deformation trial moves and smaller acceptance of deformation trial moves

4-7-14
	-In sampling loop, removed the large critical section and added two small critical sections
	inside the corresponding if statements
	-Added calculation and output of radius of gyration
	Updated from v4.1 to v4.2 -- Even with added overhead of computing the 
		rog still runs a little faster than previous version due to the critical section modification

4-25-14
	-Update for spherical rms radius boundary
	Updated from v4.2 to rmsr_v1.0 -- new boundary condition 

4-25-14
	-Converted to MC code 
	Updated to pMMC_Enthalpy_rmsr_v1.0.cpp

5-1-14
	-Added logfile output

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

#define PI 3.14159265359

// Specify the number of atoms
#define a_natoms_a a_atoms_a


using namespace std;

 // Classes
#include "./src/class/Atom.h"
#include "./src/class/Config.h"
#include "./src/class/RunningStats.h"
#include "./src/class/HistogramNS_slim.h"
//#include "./src/class/Histogram.h"

//#include "./src/class/CenterofMass.h"


// prototype functions

long double PairEnergy(const Atom& one, const Atom& two);
long double PairEnergyO(const Atom& one, const long double& ox, const long double& oy, const long double& oz);
long double TotalEnergy(const Config& c);
long double DeltaE(const Config& c, const unsigned& Aindex, const long double& ox, const long double& oy, const long double& oz);
long double ConfProb(long double E, long double Beta);
string itos( int Number );

int main(int argc, char *argv[])
{

	//Get the local time for the logfile
	time_t rawtime;
 	struct tm * timeinfo;
 	time (&rawtime);
 	timeinfo = localtime (&rawtime);

	// Parameters that need to be set

	// the number of atoms to have in the system
	const unsigned int natoms = a_natoms_a;
	//Number of Temperature Iterations
	const unsigned Tloop = 30;
	//Temp Range
	const long double Tstart = 1.0;
	const long double Tend = 0.05;
	//Total number of MC style trial moves under each NS iteration
	const unsigned int MCit = 32000000;
	//Number of procs to split Sampling loop (MCit)
	const unsigned short nprocs = n_procs_n ;
	// Number of initial MC steps to allow equilibration
	const unsigned Eqit = 1000000;
	// Step interval to collect data after equilibration
	const unsigned Cit = 100;
	
	
	long double rad = 10.0;
	
	//Translation/Volume Def step -- for rms radius condition they are the same
	const long double bsfac = 0.25;
 	long double stepguess= bsfac*rad;
    long double step = stepguess; 
	
	//Boltzmann constant
	const long double kb = 1.0;
	const long double pi = 3.14159265359;
	//Pressure
	const long double Press = p_press_p ;
	const long double rmax = 3.0;
	const long double vfact = (4.0/3.0)*pi*pow(5.0/3.0, 3.0/2.0);
	const long double Vmax = vfact*pow(rmax, 3.0);
	//const long double Htop = Press*Vmax;
	//const long double Hlow = -500.0;
 //	long double Vexc = 0.5*(4.0/3.0)*pi*(long double)(natoms);
	long double Vexc = 1.0;

	//run descriptor
	string nds = itos(natoms);
	string rundescript = "MC_Ent_n"+nds+"_rmsr_pp_press_p_reprrreprr";
	//Name of file to read in initial coordinates
	string filename = "LJa_atoms_a_GlobalMin.txt";	
//------------------------------------------------------------

	unsigned long long MCitpp = (MCit+Eqit*nprocs)/nprocs;
	unsigned int ConfigsOutFreq = ((MCitpp)-Eqit)/10;

	//Name of logfile
	string logfilename = "oMMC_" + rundescript + ".log";
	
	//Name of Configuration output file
	string coutname = "Coord_" + rundescript + ".xyz";
	//Name of Observable output file
	string outfilename = "zObservables_" + rundescript + ".dat";

	
	ofstream logfile(logfilename.c_str());
	//ofstream houtfile(houtname.c_str());
	ofstream obsoutfile(outfilename.c_str());
	//string coutfilename = "Coord_NS_Ent_a0.xyz";
	ofstream coutfile(coutname.c_str());

	//ofstream houtfile("Hpoints_MC_Ent_rmsr_a0.dat");
	//ofstream obsoutfile("zObservables_MC_Ent_rmsr_a4.dat");
	//string coutfilename = "Coord_MC_Ent_rmsr_a4.xyz";
	//ofstream coutfile(coutfilename.c_str());
//	ofstream stepoutfile("zStepOut_int_test_a1.dat");	
	// configs from each nested sampling iteration
	 

cout<<" MC iterations set to " << MCit << endl;
cout<<" Vmax: "<<Vmax<< " Vexc: "<<Vexc<<endl;	
	logfile << " " << endl;
	logfile << "Local time and date: " << asctime(timeinfo) << endl;
	logfile<<"Constant Pressure Monte Carlo Simulation of LJ particles with a spherical root mean squared radius boundary."<<endl;
	logfile << "Parameters: " << endl;
	logfile<<"Temperature Iterations: "<<Tloop<<" Number of Trial Moves (MCit): "<<MCit<<endl;
	logfile<<"Tstart: "<<Tstart<<" Tend: "<<Tend<<endl;
	logfile<<"natoms: "<<natoms<<" nprocs: "<<nprocs<<" Eqit: "<<Eqit<<" Cit: "<<Cit<<endl;
	logfile<<"Press: "<<Press<<endl;
	logfile<<"bsfac: "<<bsfac<<" initial guess: "<<stepguess<<endl;
	logfile<<"rmax: "<<rmax<<" Vmax: "<<Vmax<< " Vexc: "<<Vexc<<endl;	
	logfile<<"MCitpp: "<<MCitpp<<"ConfigsOutFreq: "<<ConfigsOutFreq<<endl;

	// Initialize the MT RNG object
	dsfmt_t dsfmt;
	
	// Need to initialize member mti as NN+1
	//Random1.mti = NN+1;
	
	//Define a seed
	//unsigned long long seed = 12983456;
	
	unsigned long long seed = time(0);
	//unsigned long long seed = 9737238;
	logfile<<"RNG seed: "<<seed<<endl;
	logfile<<endl;
	// Initialize the RNG object with the seed
	//Random1.initialize(seed2);
	dsfmt_init_gen_rand(&dsfmt, seed);

	//HistogramNS Hhist;
	//Hhist.Initialize(1000000, Hlow, Htop);
	RunningStats Estat;
	RunningStats Esqstat;
	RunningStats Vstat;
	RunningStats Hstat;
	RunningStats Rgstat;
	RunningStats Hsqstat;
	RunningStats Cpstat;
//	RunningStats UHstat;
//	RunningStats VHstat;
//Initialize configuration objects
	
Config trial1;
Config Bottom;

ifstream infile;
infile.open(filename.c_str());
// Returns an error if the file does not exist
if (!infile) {
	cerr << "Unable to open file " << filename << endl;
	cout << "Unable to open file " << filename.c_str() << endl;
}
unsigned int trace = 0;
cout<<"Reading in inital coordinates from: "<<filename<<endl;
logfile<<"Reading in inital coordinates from: "<<filename<<endl;
while (! infile.eof()){
	infile >> Bottom.atom[trace].x >> Bottom.atom[trace].y >> Bottom.atom[trace].z;
	trace++;
	if(trace == natoms){
		break;
	}
}
infile.close();
rad = Bottom.RadGyr();
long double radbottom = rad;
//Bottom.WrapCoordinates();
trial1.Equate(Bottom);
//trial1.SetBox(box, box, box);
step=bsfac*rad;	
long double Egmin = TotalEnergy(trial1);
long double Vinitial = vfact*pow(rad, 3.0);
long double Hstart = Egmin + Press*Vinitial;
long double Et1 = Egmin;
long double Ht1 = Hstart;

//long double Hnest = Htop;
//long double Htlowprev = Ht1;
//long double Vlowprev = Vinitial;
//long double rlowprev = rad;
//long double Elowprev = Et1;
unsigned long long cfailcount = 0;
//long double Hnprev=Hnest;
long double dT = (Tend - Tstart)/ (long double)(Tloop);

for(unsigned int ii=0;ii<=Tloop;++ii){

long double Tcurr = Tstart + (long double)(ii)*dT;
long double Beta = 1.0/Tcurr;

 if(ii>0){
 
//Hnest=Hhist.NestedEnergy();
//Hhist.Normalize();
//Hnest=Hhist.Median();
long double Tprev = Tcurr-dT;
long double Vavg = Vstat.Mean();
long double Eavg = Estat.Mean();
long double Havg = Hstat.Mean();
long double Rgavg = Rgstat.Mean();
long double Hsqavg = Hsqstat.Mean();
long double Cp = Cpstat.Mean();
long double cpp = (Havg*Havg - Hsqavg)/(Tprev*Tprev);
cout<<" Average Entalpy is "<<Havg<<" and Average Volume is: "<<Vavg<<" and Average Energy: "<<Eavg<<" Cp: "<<Cp<<endl;
cout<<"Final Instantaneous Cp: "<<cpp<<endl;
logfile<<" Average Entalpy is "<<Havg<<" and Average Volume is: "<<Vavg<<" and Average Energy: "<<Eavg<<" Cp: "<<Cp<<endl;
logfile<<"Final Instantaneous Cp: "<<cpp<<endl;
obsoutfile<<Tprev<<" "<<Havg<<" "<<Eavg<<" "<<Vavg<<" "<<Rgavg<<" "<<Cp<<endl;
//obsoutfile<<Eavg<<" "<<Vavg<<" "<<Rgavg<<endl;
//Hhist.SmartReset();
//Hhist.Reset();
Estat.Reset();
Esqstat.Reset();
Vstat.Reset();
Hstat.Reset();
Rgstat.Reset();
Hsqstat.Reset();
Cpstat.Reset();
//UHstat.Reset();
//VHstat.Reset();
trial1.Equate(Bottom);
Et1=Egmin;
Ht1=Hstart;
rad=radbottom;
//trial1.SetBox(box, box, box);
step=bsfac*rad*pow(0.90, cfailcount);
}
cout << "## Iteration " << ii << " Current Temp: "<<Tcurr<<" Beta: "<<Beta<<" ......." << endl;
Vinitial = vfact*pow(rad, 3.0);
cout<<"Initial Entalpy: "<<Ht1<<" Initial Energy: "<<Et1<<" Volume: "<<Vinitial<< " rms radius: "<<rad<<endl;
logfile << "## Iteration " << ii << " Current Temp: "<<Tcurr<<" Beta: "<<Beta<<" ......." << endl;
logfile<<"Initial Entalpy: "<<Ht1<<" Initial Energy: "<<Et1<<" Volume: "<<Vinitial<< " rms radius: "<<rad<<endl;
//	long double Ee = 0.0;
//	long double Cv=0.0;
//	long double rog=0.0;
//	long double He = 0.0;
//	long double Ve = 0.0;
//	unsigned long long acceptcount = 0;
	//unsigned long long vtries = 0;
	//unsigned long long vsuccess = 0;
	unsigned long long ctries = 0;
	unsigned long long csuccess = 0;
	//long double Hlowprev = Ht1;
	
	//cout<<"Moving to parallel section with nprocs: "<<nprocs<<endl;
	omp_set_num_threads(nprocs);
	#pragma omp parallel 
	{

	 //Define a thread local configuration
	 Config Tlocal;
	 Tlocal.Equate(trial1);
	// Tlocal.SetBox(box, box, box);
	// Config Tlocold;

	//Define thread local variables
	long double Et1local = Et1;
	long double Helocal = 0.0;
	long double Velocal = 0.0;
	long double Eelocal = 0.0;
	long double Ht1local = Ht1;
	long double rlocal = rad;
	 
	//unsigned long long vtriesloc = 0;
	//unsigned long long vsuccessloc = 0;
	unsigned long long ctriesloc = 0;
	unsigned long long csuccessloc = 0;
	
  //Perform markov moves on the configuration
    for(unsigned int n = 0; n<MCitpp+Eqit; ++n){
		
		
		//Store Initial Et1 value
		long double Et1b = Et1local;
		long double Ht1b = Ht1local;
		//Trial Move -----
		// randomly select an atom index	
			unsigned int Arandom = (unsigned int)(dsfmt_genrand_close_open(&dsfmt)*natoms);
		
			long double origx, origy, origz;
			origx = Tlocal.atom[Arandom].x;
			origy = Tlocal.atom[Arandom].y;
			origz = Tlocal.atom[Arandom].z;
		
		// Calculate the difference in total potential energy
		//long double dH = 0.0;
		long double dE = 0.0;
		//long double dV = 0.0;
		long double Vscale = 1.0;
		bool vdeftag = 0;
		long double Vcurr = vfact*pow(rlocal, 3.0);
		long double rscale = 1.0;
		//Coordinate move
		long double movex, movey, movez;
			/* set random movements of the coordinates
			   of the selected atom */
			movex = (dsfmt_genrand_close_open(&dsfmt)-0.5)*step;
			movey = (dsfmt_genrand_close_open(&dsfmt)-0.5)*step;
			movez = (dsfmt_genrand_close_open(&dsfmt)-0.5)*step;
			// Implement those moves
			Tlocal.atom[Arandom].x += movex;
			Tlocal.atom[Arandom].y += movey;
			Tlocal.atom[Arandom].z += movez;
		//	Tlocal.WrapAtom(Arandom);
			
			
			//Tlocal.CalcCom();
			++ctriesloc;
				long double rold = rlocal;
				//Determin new rms radius
				long double rnew = Tlocal.RadGyr();
				//Get new volume based rms r
				long double Vnew = vfact*pow(rnew, 3.0);
				Vscale = Vnew/Vcurr;
				long double dV = Vnew - Vcurr;
				long double ranprob = dsfmt_genrand_close_open(&dsfmt);
				long double vprob = 1.0;
				if (Vscale<1.0){
					long double N = (long double)(natoms);
					vprob = pow(Vscale, N);
				}
				
				if (vprob>ranprob){
					//continue;
					
				
					vdeftag=1;
				//Tlocold.Equate(Tlocal);
					rscale = pow(Vscale, 1.0/3.0);
				
					rlocal = rold*rscale;
					step=bsfac*rlocal*pow(0.90, cfailcount);
					Vcurr = Vnew;
					//Check change in energy for pairwise interactions
					dE=DeltaE(Tlocal, Arandom, origx, origy, origz);
				}
		//Selection Criteria--Metropolis Criteria
		bool tagg = 1;
		
		
		//Add difference to current total potential energy
		Et1local+=dE;
		//Calculate the change in Enthalpy from trial moves
		//dH += dE + Press*dV;
		//cout <<" dE: " << dE << endl;
		//Ht1 += dH;
		Ht1local = Et1local + Press*Vcurr;
		long double dH = Ht1local - Ht1b;
		//cout<<"NS "<<ii<<" Ht1: "<<Ht1local<<" Et1: "<<Et1local<<" Vcurr: "<<Vcurr<<endl;
		//if (condition){
	//		
	//	}
		long double prob = ConfProb(dH, Beta);
		long double rp = dsfmt_genrand_close_open(&dsfmt);
		

		if(prob > rp && vdeftag==1){
			//cout<<"NS "<<ii<<" Trial move was accepted! " << endl;
			tagg=0;
			Eelocal=Et1local;
			Helocal = Ht1local;
			Velocal = Vcurr;
			//cout<<"NS "<<ii<<" He: "<<Helocal<<" Ee: "<<Eelocal<<" Ve: "<<Velocal<<endl;
			
			//++acceptcount;
						
			
			++csuccessloc;			
			
		}


				/*Failed criteria--r(dsfmt_genrand_close_open(&dsfmt)-0.5)*vstep;eset atom coordinates and remove
		  the difference in energy */ 
		if(tagg == 1){
			
			rlocal/=rscale;
			Tlocal.atom[Arandom].x = origx;
			Tlocal.atom[Arandom].y = origy;
			Tlocal.atom[Arandom].z = origz;
			step=bsfac*rlocal*pow(0.80, cfailcount);
			Et1local = Et1b;
			Ht1local = Ht1b;
			
		}

			if(n>Eqit && n%(Cit)==0){
			//	cout << "!!!!!Going to push values!!!!.................--------------------------!!!!!" << endl;
			//	cout << "Ee: " << Eelocal << " He: " << Helocal << " Ve: " << Velocal << endl;	
			//	cout << "He*He: " << He*He << endl;	
				#pragma omp critical
				{
					Estat.Push(Eelocal);
					Esqstat.Push(Eelocal*Eelocal);
					//Hhist.Push(Helocal);
					//long double rog = Tlocal.RadGyr();
					Rgstat.Push(rlocal);
					Hstat.Push(Helocal);
					Vstat.Push(Velocal);
					Hsqstat.Push(Helocal*Helocal);
				//	UHstat.Push(Eelocal*Helocal);
				//	VHstat.Push(Velocal*Helocal);
					long double Ha = Hstat.Mean();
					long double Hsa = Hsqstat.Mean();
					//long double Ea = Estat.Mean();
					//long double Va = Vstat.Mean();
					//long double UHa = UHstat.Mean();
					//long double VHa = VHstat.Mean();
					long double delH = Hsa-(Ha*Ha);
					delH/=(Tcurr*Tcurr);
					long double cp = delH;
					Cpstat.Push(cp);
				}
			
		}


		//Output coordinates of sample n
	   if(n>Eqit && n % (ConfigsOutFreq) == 0) {
		#pragma omp critical
		{
		coutfile << Tlocal.natoms << endl;
		coutfile <<"lj 17 H: "<<Ht1local<<" E: "<<Et1local<<" V: "<<Vcurr<<endl;
		for(unsigned long hg = 0; hg<Tlocal.natoms; ++hg){
			coutfile << setiosflags(ios::fixed) << setprecision(15) << "lj " << Tlocal.atom[hg].x << " " << Tlocal.atom[hg].y << " " << Tlocal.atom[hg].z << endl; 
		}
		}
	   }	

	

	//	cout << endl;
	
		
	


		
    } //End Monte Carlo loop

	#pragma omp barrier
	
	#pragma omp critical
	{ 
	//	vsuccess+=vsuccessloc;
		csuccess+=csuccessloc;
	//	vtries+=vtriesloc;
		ctries+=ctriesloc;

	}

}//End parallel

//	long double vafrac = (long double)(vsuccess)/(long double)(vtries);
	long double cafrac = (long double)(csuccess)/(long double)(ctries);
//	if (vafrac<0.30){
//		vstep*=0.90;
//	}
	if (cafrac<0.30){
		++cfailcount;
	}

cout<<" cfrac: "<<cafrac<<" cfailcount: "<<cfailcount<<endl;
logfile<<" cfrac: "<<cafrac<<" cfailcount: "<<cfailcount<<endl;
} // End T loop

//houtfile.close();
obsoutfile.close();
coutfile.close();
//stepoutfile.close();
logfile<<"Simulation Complete!"<<endl;
logfile<<"Output files are: "<<endl;
logfile<<"Configuration Samples: "<<coutname<<endl;
logfile<<"Average Observable Values: "<<outfilename<<endl;
cout<<"Simulation Complete!"<<endl;
cout<<"Output files are: "<<endl;
cout<<"Logfile: "<<logfilename<<endl;
cout<<"Configuration Samples: "<<coutname<<endl;
cout<<"Average Observable Values: "<<outfilename<<endl;
cout<<endl;


cout << " " << endl;
cout << " Thank You and Have a Nice Day!" << endl;
cout << " " << endl;
cout << "... End of Line ..." << endl;
	return EXIT_SUCCESS;



}

// Define Functions

// Total potential energy calculator, uses only the lj pot
long double TotalEnergy(const Config& c){
	
	long double Etot = 0.0;
	unsigned int natoms = c.natoms;

		for(unsigned int pp = 0;pp<natoms-1;++pp){

			for(unsigned int uu = pp+1;uu<natoms;++uu){
				
				
				Etot += PairEnergy(c.atom[pp], c.atom[uu]);
					
			}
		}
//	Etot *= 60.0;
	return(Etot);
}

// potential energy function
long double PairEnergy(const Atom& one, const Atom& two){
	
	
	long double dx, dy, dz, rs;
	long double potent, eps, sigma;
	eps = 1.0, sigma = 1.0;

	long double rc = 3.0*sigma;
	long double rcs = rc*rc;
	
	dx = one.x - two.x;
	
	dy = one.y - two.y;
	
	dz = one.z - two.z;
	rs= (pow(dx, 2))+(pow(dy, 2))+(pow(dz, 2));
	if (rs<rcs){
		potent = 4.0*eps*((pow(sigma/rs, 6))-(pow(sigma/rs, 3)));
	}
	else{
		potent=0.0;
	}
	
	return(potent);
}

//Uses the o position values (xo, yo, zo) from atom one
long double PairEnergyO(const Atom& one, const long double& ox, const long double& oy, const long double& oz){
	
	
	long double dx, dy, dz, rs;
	long double potent, eps, sigma;
	
	eps = 1.0, sigma = 1.0;
	long double rc = 3.0*sigma;
	long double rcs = rc*rc;
	dx = one.x - ox;
	
	dy = one.y - oy;
	
	dz = one.z - oz;

	
	rs= (pow(dx, 2))+(pow(dy, 2))+(pow(dz, 2));
	if (rs<rcs){
		potent = 4.0*eps*((pow(sigma/rs, 6))-(pow(sigma/rs, 3)));
	}
	else{
		potent=0.0;
	}
	
	return(potent);
}

// calculates the change in energy associated with one moved atom
long double DeltaE(const Config& c, const unsigned& Aindex, const long double& ox, const long double& oy, const long double& oz){
	long double Etot=0.0;	
	long double Etoto=0.0;
	long double dE;
	
		
	
	for (unsigned int i = 0; i < c.natoms; ++i)
	{
		if(i!=Aindex){
			Etot += PairEnergy(c.atom[Aindex], c.atom[i]);
			Etoto += PairEnergyO(c.atom[i], ox, oy, oz);
					
		}	
	}
	dE = Etot - Etoto;
	//cout << "In DeltaE funct, dE is : " << dE << endl;
	return( dE ); 
	

}

long double ConfProb(long double E, long double Beta){
//	cout << "Beta " << Beta << endl;
	long double P = exp(-Beta*E);
	return(P);
}

string itos( int Number ){
     ostringstream ss;
     ss << Number;
     return ss.str();
 }
