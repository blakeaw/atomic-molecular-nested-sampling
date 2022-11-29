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

5-1-14
	-Added logfile output

5-16-14
	-Added thread local walker minima tracking and usage
	Updated from v1.0 to v1.1 - should have enhanced sampling over previous version

6-4-14
	-Converted to spherical com boundary
	Updated from rmsr_v1.1 to comr_v1.0

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



using namespace std;

 // Classes
#include "./src/class/Atom.h"
#include "./src/class/Config.h"
#include "./src/class/RunningStats.h"
#include "./src/class/HistogramNS_slim.h"
#include "./src/class/Histogram.h"

//#include "./src/class/CenterofMass.h"


// prototype functions

long double PairEnergy(const Atom& one, const Atom& two);
long double PairEnergyO(const Atom& one, const long double& ox, const long double& oy, const long double& oz);
long double TotalEnergy(const Config& c);
long double DeltaE(const Config& c, const unsigned& Aindex, const long double& ox, const long double& oy, const long double& oz);
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
	const unsigned int natoms = 17;
	
	// the number of nested sampling iterations to run 
	const unsigned int Iterations = 100;
	//Total number of MC style trial moves under each NS iteration
	const unsigned int MCit = 2000000;
	//Number of procs to split Sampling loop (MCit)
	const unsigned short nprocs = 1;
	// Number of initial MC steps to allow equilibration
	const unsigned Eqit = 500000;
	//Step Interval to perform volume deformation
	const unsigned Vint = 4*natoms;
	// Step interval to collect data after equilibration
	const unsigned Cit = 5*Vint;
	//Boltzmann constant
	const long double kb = 1.0;
	const long double pi = 3.14159265359;
	//Pressure
	const long double Press = 1.0;
	
	 long double Htop = 500.0;
	 long double Hlow = -100.0;
	 unsigned long long nhbins = 1000000;
	const long double vfact = (4.0/3.0)*pi;
	long double Vexc = 1.0;
	const long double Vmax = (0.5*Htop) / Press;
	/*Spherical COM boundary condition radius */
	long double dconstraint;
	dconstraint = pow(Vmax/vfact, 1.0/3.0);
	cout<<"Vmax: "<<Vmax<<" dconstraint: "<<dconstraint<<" vact: "<<vfact<<" Htop: "<<Htop<<endl;
	//Translation step
 	long double stepguess= 2.0*dconstraint;
    long double step = stepguess; 
	
	//Volume deformation step
	long double vstepguess = 0.1;
	long double vstep = 0.0;	


	//run descriptor
	string nds = itos(natoms);
	string rundescript = "NS_Ent_n"+nds+"_comr_a4";
	//Name of file to read in initial coordinates
	string filename = "LJ17_GlobalMin.txt";	
//------------------------------------------------------------

	unsigned long long MCitpp = (MCit+Eqit*nprocs)/nprocs;
	unsigned int ConfigsOutFreq = ((MCitpp)-Eqit)/5;

	//Name of logfile
	string logfilename = "oMMC_" + rundescript + ".log";
	//Name of Emedian output file
	string houtname = "Hpoints_" + rundescript + ".dat";
	//Name of Configuration output file
	string coutname = "Coord_" + rundescript + ".xyz";
	//Name of Observable output file
	string outfilename = "zObservables_" + rundescript + ".dat";
	string hhistname = "zHhist_"+rundescript;
	
	ofstream logfile(logfilename.c_str());
	ofstream houtfile(houtname.c_str());
	ofstream obsoutfile(outfilename.c_str());
	//string coutfilename = "Coord_NS_Ent_a0.xyz";
	ofstream coutfile(coutname.c_str());

	logfile << " " << endl;
	logfile << "Local time and date: " << asctime(timeinfo) << endl;
	logfile<<"Constant Pressure Nested Sampling Simulation of LJ particles with a spherical root mean squared radius boundary."<<endl;
	logfile << "Parameters: " << endl;
	logfile<<"NS Iterations: "<<Iterations<<" Number of Trial Moves (MCit): "<<MCit<<endl;
	logfile<<"natoms: "<<natoms<<" nprocs: "<<nprocs<<" Eqit: "<<Eqit<<" Cit: "<<Cit<<endl;
	logfile<<"Press: "<<Press<<endl;
	logfile<<"vstepguess: "<<vstepguess<<" initial guess: "<<stepguess<<endl;
	logfile<<"dconstraint: "<<dconstraint<<" Htop: "<<Htop<<" Hlow: "<<Hlow<<" nhbins: "<<nhbins<<" Vmax: "<<Vmax<< " Vexc: "<<Vexc<<endl;	
	logfile<<"MCitpp: "<<MCitpp<<"ConfigsOutFreq: "<<ConfigsOutFreq<<endl;
//	ofstream stepoutfile("zStepOut_int_test_a1.dat");	
	// configs from each nested sampling iteration
	 
	

cout << "Nested sampling iterations is set to " << Iterations << " with MC iterations set to " << MCit << endl;
cout<<"Htop: "<<Htop<<" Hlow: "<<Hlow<<" Vmax: "<<Vmax<< " Vexc: "<<Vexc<<endl;	
	

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
	long double bwd = 1.0e-4;
	HistogramNS Hhist;
	Hhist.Initialize(bwd, Hlow, Htop);
	RunningStats Estat;
	RunningStats Esqstat;
	RunningStats Vstat;
	RunningStats Hstat;
	RunningStats Rgstat;
//Initialize configuration objects
	Histogram Hohist;	

	Hohist.Initialize(bwd*1000.0, Hlow, Htop);
Config trial1(natoms);
Config Bottom(natoms);
// local min confs
Config * lgmin  = new Config[nprocs];

// Initialize an array to store the lgmin Energy values
long double Hnmin[nprocs]={Htop};
long double Enmin[nprocs]={0.0};
// store lgmin box size
long double lgminr[nprocs]={0.0};
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
//rad = Bottom.RadGyr();
//Bottom.WrapCoordinates();
Bottom.CalcCom();
Bottom.ReCenterZero();
trial1.Equate(Bottom);
//trial1.SetBox(box, box, box);
//step=bsfac*rad;	
long double Egmin = TotalEnergy(trial1);
long double Vinitial = vfact*pow(dconstraint, 3.0);
long double Hstart = Egmin + Press*Vinitial;
long double Et1 = Egmin;
long double Ht1 = Hstart;
Hstart=Egmin;
Ht1=Et1;
long double Hnest = Htop;
long double Htlowprev = Ht1;
long double Vlowprev = Vinitial;
long double dclowprev = dconstraint;
long double Elowprev = Et1;
unsigned long long cfailcount = 0;
long double Hnprev=Hnest;

for (unsigned int np = 0; np < nprocs; ++np)
{
	lgmin[np].InitializeAtoms(natoms);
	lgmin[np].Equate(Bottom);
//	lgmin[np].SetBox(box, box, box);
	//Hnmin[np]=Htlowprev;
	Hnmin[np]=Htop;
	Enmin[np]=Elowprev;
	lgminr[np]=dclowprev;
	//cout<<"min "<<np<<" has "<<lgmin[np].natoms<<" atoms"<<endl;
}	


for(unsigned int ii=0;ii<Iterations;++ii){
cout << "Iteration " << ii << "......." << endl;
logfile << "Iteration " << ii << "......." << endl;
 if(ii>0){
 
Hnest=Hhist.NestedEnergy();
string hname = hhistname+"_iter"+itos(ii)+".dat";
Hohist.Normalize();
Hohist.DumpNormal(hname);
Hohist.Reset();
//Hhist.Normalize();
//Hnest=Hhist.Median();
long double Vavg = Vstat.Mean();
long double Eavg = Estat.Mean();
long double Havg = Hstat.Mean();
long double Rgavg = Rgstat.Mean();
cout<<"Nested Entalpy is: "<<Hnest<<" and Average Entalpy is "<<Havg<<" and Average Volume is: "<<Vavg<<" and Average Energy: "<<Eavg<<endl;
logfile<<"Nested Entalpy is: "<<Hnest<<" and Average Entalpy is "<<Havg<<" and Average Volume is: "<<Vavg<<" and Average Energy: "<<Eavg<<endl;
//obsoutfile<<Havg<<" "<<Eavg<<" "<<Vavg<<" "<<Rgavg<<endl;
obsoutfile<<Eavg<<" "<<Vavg<<" "<<Rgavg<<endl;
Hhist.SmartReset();
//Hhist.Reset();
Estat.Reset();
Esqstat.Reset();
Vstat.Reset();
Hstat.Reset();
Rgstat.Reset();
//trial1.Equate(Bottom);
//Et1=Elowprev;
//Ht1=Htlowprev;
//rad=rlowprev;
//trial1.SetBox(box, box, box);
//step=bsfac*rad*pow(0.85, cfailcount);
if (abs(Hnest-Hnprev)<0.0001){
	break;
}
Hnprev=Hnest;


//Make sure starting conf is lower enthalpy than new nested value - if not set to global found minima
for (unsigned int np = 0; np < nprocs; ++np)
{
	if (Hnmin[np]>Hnest){
		lgmin[np].Equate(Bottom);
		Hnmin[np]=Htlowprev;
		Enmin[np]=Elowprev;
		lgminr[np]=dclowprev;
		//lgmin[np].SetBox(blowprev, blowprev, blowprev);
	}
}


}
houtfile<<Hnest<<endl;


//	long double Ee = 0.0;
//	long double Cv=0.0;
//	long double rog=0.0;
//	long double He = 0.0;
//	long double Ve = 0.0;
//	unsigned long long acceptcount = 0;
	unsigned long long vtries = 0;
	unsigned long long vsuccess = 0;
	unsigned long long ctries = 0;
	unsigned long long csuccess = 0;
	//long double Hlowprev = Ht1;
	
	//cout<<"Moving to parallel section with nprocs: "<<nprocs<<endl;
	omp_set_num_threads(nprocs);
	#pragma omp parallel 
	{

	  //local thread number
	 unsigned tid = omp_get_thread_num();
	 //Define a thread local configuration
	 Config Tlocal(natoms);
	long double Et1local;
	long double Ht1local;
	long double dclocal;
	Config Vlocal(natoms);
	 //set starting conf
	 if (ii==0){
	 	Tlocal.Equate(Bottom);
		//Tlocal.SetBox(blowprev, blowprev, blowprev);
		Et1local=Elowprev;
		Ht1local=Htlowprev;
		dclocal = dclowprev;
		
	 }
	 else{
		Tlocal.Equate(lgmin[tid]);
		dclocal = lgminr[tid];
		//Tlocal.SetBox(blocal, blocal, blocal);
		Et1local = Enmin[tid];
		Ht1local = Hnmin[tid];
		//cout<<"initial values thread "<<tid<<"  box: "<<lgminbox[tid]<<" Et1: "<<Enmin[tid]<<" Ht1: "<<Hnmin[tid]<<endl;
	 }
	Vinitial = vfact*pow(dclocal, 3.0);
	step=0.98*dclocal*pow(0.85, cfailcount);
//long double Vinitial = pow(blocal, 3);
	#pragma omp critical
	{
		cout<<"thread walker: "<<tid<<" Initial Entalpy: "<<Ht1local<<" Initial Energy: "<<Et1local<<" Volume: "<<Vinitial<<" dclocal: "<<dclocal<<endl;
		logfile<<"thread walker: "<<tid<<" Initial Entalpy: "<<Ht1local<<" Initial Energy: "<<Et1local<<" Volume: "<<Vinitial<<endl;
	}
	//Define thread local variables
	//long double Et1local = Et1;
	long double Helocal = 0.0;
	long double Velocal = 0.0;
	long double Eelocal = 0.0;
	//long double Ht1local = Ht1;
	//long double rlocal = rad;
	 
	unsigned long long vtriesloc = 0;
	unsigned long long vsuccessloc = 0;
	unsigned long long ctriesloc = 0;
	unsigned long long csuccessloc = 0;
	
  //Perform markov moves on the configuration
    for(unsigned int n = 0; n<MCitpp; ++n){
		
		
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
		long double Vcurr = vfact*pow(dclocal, 3.0);
		long double rscale = 1.0;
		bool dflag = 1;
		//Volume deformation move
		if(n%Vint==0){
				++vtriesloc;
				long double rold = dclocal;
				long double lnVnew = log(Vcurr) + (dsfmt_genrand_close_open(&dsfmt)-0.5)*vstep;
				long double Vnew = exp(lnVnew);
				while (Vnew<Vexc){
					lnVnew = log(Vcurr) + (dsfmt_genrand_close_open(&dsfmt)-0.5)*vstep;
					Vnew = exp(lnVnew);
				}
			//	dV = Vnew - Vcurr;
				
				
				
				Vscale = Vnew/Vcurr;
				long double ranprob = dsfmt_genrand_close_open(&dsfmt);
				long double vprob = 1.0;
				if (Vscale<1.0){
					long double N = (long double)(natoms);
					vprob = pow(Vscale, N);
				}
				if (vprob>ranprob){
					long double r3 = Vnew/vfact;
					dclocal = pow(r3, 1.0/3.0);
					rscale = dclocal/rold;
					Vlocal.Equate(Tlocal);
					Tlocal.ScaleCoordinates(Vscale);
					long double Enew = TotalEnergy(Tlocal);
			//	dE += Enew - Et1;
					Et1local = Enew;
					vdeftag = 1;
					step = 0.98*dclocal*pow(0.95, cfailcount);
					Vcurr = Vnew;
					
				}
		}
		else{
		//Coordinate change		

			
			
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
			//Check change in energy for pairwise interactions
			
			Tlocal.CalcCom();
		
			dflag = Tlocal.CheckComDistAll(dclocal);
			
		
			if(dflag){
				dE=DeltaE(Tlocal, Arandom, origx, origy, origz);
				Et1local+=dE;
			}
			++ctriesloc;
		}
		
		//Selection Criteria--Metropolis Criteria
		bool tagg = 1;
		
		
		//Add difference to current total potential energy
	//	Et1local+=dE;
		//Calculate the change in Enthalpy from trial moves
		//dH += dE + Press*dV;
		//cout <<" dE: " << dE << endl;
		//Ht1 += dH;
		//Ht1local = Et1local + Press*Vcurr;
		Ht1local=Et1local;
		//cout<<"NS "<<ii<<" Ht1: "<<Ht1local<<" Et1: "<<Et1local<<" Vcurr: "<<Vcurr<<endl;
		//if (condition){
	//		
	//	}
		if(Ht1local < Hnest && dflag){
			//cout<<"NS "<<ii<<" Trial move was accepted! " << endl;
			tagg=0;
			Eelocal=Et1local;
			Helocal = Ht1local;
			Velocal = Vcurr;
			//cout<<"NS "<<ii<<" He: "<<Helocal<<" Ee: "<<Eelocal<<" Ve: "<<Velocal<<endl;
		//	if (Eelocal<-500.0){
		//		exit(0);
		//	}
			//++acceptcount;

			if (Helocal<Hnmin[tid]){
				
				Hnmin[tid]=Helocal;
				Enmin[tid]=Eelocal;
				lgminr[tid]=dclocal;
				lgmin[tid].Equate(Tlocal);
				//lgmin[tid].SetBox(blocal, blocal, blocal);
				//cout<<"local min for thread: "<<tid<<" setting!"<<endl;
				//cout<<"Hnmin: "<<Hnmin[tid]<<" Helocal "<<Helocal<<endl;
				//Vnmin[tid]=Velocal;
			}
			if (Helocal<Htlowprev){

				#pragma omp critical
				{
					Htlowprev=Helocal;
					Vlowprev=Velocal;
					dclowprev=dclocal;
					Elowprev=Eelocal;
					Bottom.Equate(Tlocal);
					//Bottom.SetBox(blocal, blocal, blocal);
				}
			}
			
			
			if (vdeftag==1){
				++vsuccessloc;
			}
			else{
				++csuccessloc;			
			}	
			
		}


				/*Failed criteria--r(dsfmt_genrand_close_open(&dsfmt)-0.5)*vstep;eset atom coordinates and remove
		  the difference in energy */ 
		if(tagg == 1){
			
			if(vdeftag==1){
				long double Unscale = 1.0/Vscale;
				//Tlocal.ScaleCoordinates(Unscale);
				Tlocal.Equate(Vlocal);
				dclocal/=rscale;
				step = 0.98*dclocal*pow(0.95, cfailcount);
				Vcurr = vfact*pow(dclocal, 3);
			}
			else{
				Tlocal.atom[Arandom].x = origx;
				Tlocal.atom[Arandom].y = origy;
				Tlocal.atom[Arandom].z = origz;
			}
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
					Hhist.Push(Helocal);
					long double rog = Tlocal.RadGyr();
					Rgstat.Push(rog);
					Hstat.Push(Helocal);
					Vstat.Push(Velocal);
					Hohist.Push(Eelocal);
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
		vsuccess+=vsuccessloc;
		csuccess+=csuccessloc;
		vtries+=vtriesloc;
		ctries+=ctriesloc;

	}

}//End parallel

	long double vafrac = (long double)(vsuccess)/(long double)(vtries);
	long double cafrac = (long double)(csuccess)/(long double)(ctries);
	if (vafrac<0.30){
		vstep*=0.90;
	}
	if (cafrac<0.30){
		++cfailcount;
	}
	cout<<"vfrac: "<<vafrac<<" cfrac: "<<cafrac<<" cfailcount: "<<cfailcount<<endl;
	cout<<"vstep: "<<vstep<<endl;
//cout<<" cfrac: "<<cafrac<<" cfailcount: "<<cfailcount<<endl;
cout<<"Outputting lowest enthalpy configuration.."<<endl;
logfile<<" cfrac: "<<cafrac<<" cfailcount: "<<cfailcount<<endl;
logfile<<"Outputting lowest enthalpy configuration.."<<endl;
		coutfile << trial1.natoms << endl;
		coutfile <<"lj 17 H: "<<Htlowprev<<" E: "<<Elowprev<<" V: "<<Vlowprev<<endl;
		for(unsigned long hg = 0; hg<trial1.natoms; ++hg){
			coutfile << setiosflags(ios::fixed) << setprecision(15) << "lj " << Bottom.atom[hg].x << " " << Bottom.atom[hg].y << " " << Bottom.atom[hg].z << endl; 
		}

} // End NS loop

houtfile.close();
obsoutfile.close();
coutfile.close();
//stepoutfile.close();
logfile<<"Simulation Complete!"<<endl;
logfile<<"Output files are: "<<endl;
logfile<<"Nested Enthalpies: "<<houtname<<endl;
logfile<<"Configuration Samples: "<<coutname<<endl;
logfile<<"Average Observable Values: "<<outfilename<<endl;
cout<<"Simulation Complete!"<<endl;
cout<<"Output files are: "<<endl;
cout<<"Logfile: "<<logfilename<<endl;
cout<<"Nested Enthalpies: "<<houtname<<endl;
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

string itos( int Number ){
     ostringstream ss;
     ss << Number;
     return ss.str();
 }
