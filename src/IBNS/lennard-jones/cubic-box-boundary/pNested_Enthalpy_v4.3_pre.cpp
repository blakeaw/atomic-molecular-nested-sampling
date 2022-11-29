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

5-14-14
	-Added logfile output
	-Added convergence criteria exit

5-16-14
	-Added thread local walker minima tracking and usage
	Updated from v4.2 to v4.3 - should have enhanced sampling over previous version

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
// Specify the number of atoms
//#define a_natoms_a 17


using namespace std;

 // Classes
#include "/net/uu/nm/cm/bxw109120/NestedSampling/ConstPress/EnthalpyNested/ljcluster/src/class/Atom.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/ConstPress/EnthalpyNested/ljcluster/src/class/Config.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/ConstPress/EnthalpyNested/ljcluster/src/class/RunningStats.h"
//#include "/net/uu/nm/cm/bxw109120/NestedSampling/ConstPress/EnthalpyNested/ljcluster/src/class/HistogramNS_slim.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/ConstPress/EnthalpyNested/ljcluster/src/class/Histogram.h"
//#include "./src/class/Histogram.h"

//#include "./src/class/CenterofMass.h"


// prototype functions

long double PairEnergy(const Atom& one, const Atom& two, const long double& box);
long double PairEnergyO(const Atom& one, const long double& ox, const long double& oy, const long double& oz, const long double& box);
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
	const unsigned int natoms = a_atoms_a ;
	
	// the number of nested sampling iterations to run 
	const unsigned int Iterations = 5000 ;
	//Total number of MC style trial moves under each NS iteration
	const unsigned int MCit = 1000000 ;
	//Number of procs to split Sampling loop (MCit)
	const unsigned short nprocs = n_procs_n ;
	// Number of initial MC steps to allow equilibration
	const unsigned Eqit = 40000*natoms ;
	//Step Interval to perform volume deformation
	const unsigned Vint = 1;
	// Step interval to collect data after equilibration
	const unsigned Cit = 10*Vint;
	/*Spherical COM boundary condition radius */
//	long double dconstraint = 6.60;
	long double box = 14.0;
	//Dump nested sampling distribution - 0=no, 1=yes
	bool dnest = 0;
	//Translation step
	const long double bsfac = 0.98;
 	long double stepguess= bsfac*box;
    long double step = stepguess; 
	
	//Volume deformation step
	long double vstepguess = 5.0;
	long double vstep = vstepguess;	
	//Boltzmann constant
	const long double kb = 1.0;
	const long double pi = 3.14159265359;
	//Pressure
	const long double Press = p_press_p ;
	
	const long double Htop = 800.0;
	const long double Hlow = -1.0;
	const long double binwidth = 1.0e-4;
	const long double Vmax = (Htop/Press);
	const long double bmax = pow(Vmax, 1.0/3.0);
	box = 0.8*bmax;
 //	long double Vexc = 0.5*(4.0/3.0)*pi*(long double)(natoms);
	long double Vexc = 1.0e-3;
   
	//run descriptor
	string nds = itos(natoms);
	string rundescript = "NS_Ent_n"+nds+"_cbb_pp_press_p_rep_rrreprr";
	//Name of file to read in initial coordinates
	string filename = "LJa_atoms_a_GlobalMin.txt";	
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

	
	ofstream logfile(logfilename.c_str());
	ofstream houtfile(houtname.c_str());
	ofstream obsoutfile(outfilename.c_str());
	//string coutfilename = "Coord_NS_Ent_a0.xyz";
	ofstream coutfile(coutname.c_str());
	unsigned long long nbins =(unsigned long long)( (Htop-Hlow)/binwidth );
	logfile << " " << endl;
	logfile << "Local time and date: " << asctime(timeinfo) << endl;
	logfile<<"Constant Pressure Nested Sampling Simulation of LJ particles with a periodic cubic boundary using minimum image interaction."<<endl;
	logfile << "Parameters: " << endl;
	logfile<<"NS Iterations: "<<Iterations<<" Number of Trial Moves (MCit): "<<MCit<<endl;
	logfile<<"natoms: "<<natoms<<" nprocs: "<<nprocs<<" Eqit: "<<Eqit<<" Cit: "<<Cit<<endl;
	logfile<<"Press: "<<Press<<endl;
	logfile<<"bsfac: "<<bsfac<<" initial guess: "<<stepguess<<endl;
	logfile<<"bmax: "<<bmax<<" Htop: "<<Htop<<" Hlow: "<<Hlow<<" binwidth: "<<binwidth<<" nbins: "<<nbins<<" Vmax: "<<Vmax<< " Vexc: "<<Vexc<<endl;	
	logfile<<"MCitpp: "<<MCitpp<<"ConfigsOutFreq: "<<ConfigsOutFreq<<endl;

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
	
	// Initialize the RNG object with the seed
	//Random1.initialize(seed2);
	dsfmt_init_gen_rand(&dsfmt, seed);

	Histogram Hhist;
	Hhist.Initialize(binwidth, Hlow, Htop);
	//Hhist.SetNSfrac(0.99);
	RunningStats Estat;
	RunningStats Esqstat;
	RunningStats Vstat;
	RunningStats Hstat;
	RunningStats Rgstat;
//Initialize configuration objects
	
Config trial1(natoms);
// Global min conf
Config Bottom(natoms);
// local min confs
Config * lgmin  = new Config[nprocs];

// Initialize an array to store the lgmin Energy values
long double Hnmin[nprocs]={Htop};
long double Enmin[nprocs]={0.0};
// store lgmin box size
long double lgminbox[nprocs]={0.0};
//string filename = "LJ50_GlobalMin.txt";	
ifstream infile;
infile.open(filename.c_str());
// Returns an error if the file does not exist
if (!infile) {
	cerr << "Unable to open file " << filename << endl;
//	cout << "Unable to open file " << filename.c_str() << endl;
	cerr<<"Exiting simulation..."<<endl;
	exit(0);
}
unsigned int trace = 0;
while (! infile.eof()){
	infile >> Bottom.atom[trace].x >> Bottom.atom[trace].y >> Bottom.atom[trace].z;
	trace++;
	if(trace == natoms){
		break;
	}
}
infile.close();
Bottom.SetBox(box, box, box);
//Bottom.WrapCoordinates();
trial1.Equate(Bottom);
trial1.SetBox(box, box, box);

long double Egmin = TotalEnergy(trial1);
long double Hstart = Egmin + Press*pow(box, 3);
long double Et1 = Egmin;
long double Ht1 = Hstart;

long double Hnest = Htop;
long double Htlowprev = Ht1;
long double Vlowprev = pow(box, 3);
long double blowprev = box;
long double Elowprev = Et1;
unsigned long long cfailcount = 0;
long double Hnprev=Hnest;

for (unsigned int np = 0; np < nprocs; ++np)
{
	lgmin[np].InitializeAtoms(natoms);
	lgmin[np].Equate(Bottom);
	lgmin[np].SetBox(box, box, box);
	//Hnmin[np]=Htlowprev;
	Hnmin[np]=Htop;
	Enmin[np]=Elowprev;
	lgminbox[np]=blowprev;
	//cout<<"min "<<np<<" has "<<lgmin[np].natoms<<" atoms"<<endl;
}	

for(unsigned int ii=0;ii<Iterations;++ii){
cout << "Iteration " << ii << "......." << endl;
logfile << "Iteration " << ii << "......." << endl;
 if(ii>0){
 
//Hnest=Hhist.NestedEnergy();
Hhist.Normalize();
Hnest=Hhist.Median();
long double Vavg = Vstat.Mean();
long double Eavg = Estat.Mean();
long double Havg = Hstat.Mean();
long double Rgavg = Rgstat.Mean();
cout<<"Nested Entalpy is: "<<Hnest<<" and Average Entalpy is "<<Havg<<" and Average Volume is: "<<Vavg<<" and Average Energy: "<<Eavg<<endl;
logfile<<"Nested Entalpy is: "<<Hnest<<" and Average Entalpy is "<<Havg<<" and Average Volume is: "<<Vavg<<" and Average Energy: "<<Eavg<<endl;
//obsoutfile<<Havg<<" "<<Eavg<<" "<<Vavg<<" "<<Rgavg<<endl;
obsoutfile<<Eavg<<" "<<Vavg<<" "<<Rgavg<<endl;
if (dnest==1){
	string fname = "dNestDist_"+rundescript+"_iter_"+itos(ii-1)+".dat";
	cout<<"Dumping Nested Distribution to file: "<<fname.c_str()<<" for iteration "<<ii-1<<endl;
///	Hhist.DumpNested(fname);
}
//Hhist.SmartReset();
Hhist.Reset();
Estat.Reset();
Esqstat.Reset();
Vstat.Reset();
Hstat.Reset();
Rgstat.Reset();
//trial1.Equate(Bottom);
//Et1=Elowprev;
//Ht1=Htlowprev;
//box=blowprev;
//trial1.SetBox(box, box, box);
//step=bsfac*box*pow(0.90, cfailcount);
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
		lgminbox[np]=blowprev;
		lgmin[np].SetBox(blowprev, blowprev, blowprev);
	}
}


}
houtfile<<Hnest<<endl;
//long double Vinitial = pow(box, 3);
//cout<<"Initial Entalpy: "<<Ht1<<" Initial Energy: "<<Et1<<" Volume: "<<Vinitial<<endl;
//logfile<<"Initial Entalpy: "<<Ht1<<" Initial Energy: "<<Et1<<" Volume: "<<Vinitial<<endl;
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
	long double Hlowprev = Ht1;
	

	omp_set_num_threads(nprocs);
	#pragma omp parallel 
	{
	
	 //local thread number
	 unsigned tid = omp_get_thread_num();
	 //Define a thread local configuration
	 Config Tlocal(natoms);
	long double Et1local;
	long double Ht1local;
	long double blocal;

	 //set starting conf
	 if (ii==0){
	 	Tlocal.Equate(Bottom);
		Tlocal.SetBox(blowprev, blowprev, blowprev);
		Et1local=Elowprev;
		Ht1local=Htlowprev;
		blocal = blowprev;
		
	 }
	 else{
		Tlocal.Equate(lgmin[tid]);
		blocal = lgminbox[tid];
		Tlocal.SetBox(blocal, blocal, blocal);
		Et1local = Enmin[tid];
		Ht1local = Hnmin[tid];
		//cout<<"initial values thread "<<tid<<"  box: "<<lgminbox[tid]<<" Et1: "<<Enmin[tid]<<" Ht1: "<<Hnmin[tid]<<endl;
	 }
	//Hlowprev = Ht1;
	// Tlocal.Equate(trial1);
	// Tlocal.SetBox(box, box, box);
	// Config Tlocold;
	step=bsfac*box*pow(0.90, cfailcount);	
	long double Vinitial = pow(blocal, 3);
	#pragma omp critical
	{
		cout<<"thread walker: "<<tid<<" Initial Entalpy: "<<Ht1local<<" Initial Energy: "<<Et1local<<" Volume: "<<Vinitial<<endl;
		logfile<<"thread walker: "<<tid<<" Initial Entalpy: "<<Ht1local<<" Initial Energy: "<<Et1local<<" Volume: "<<Vinitial<<endl;
	}
	//Define thread local variables
//	long double Et1local = Et1;
	long double Helocal = 0.0;
	long double Velocal = 0.0;
	long double Eelocal = 0.0;
	//long double Ht1local = Ht1;
	//long double blocal = box;
	 
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
		long double Vcurr = pow(blocal, 3);
		long double bscale = 1.0;
		//Volume deformation move
		if(n%Vint==0){
				++vtriesloc;
				long double bold = blocal;
				//long double lnVnew = log(Vcurr) + (dsfmt_genrand_close_open(&dsfmt)-0.5)*vstep;
				//long double Vnew = exp(lnVnew);
				long double Vnew = Vcurr + (dsfmt_genrand_close_open(&dsfmt)-0.5)*vstep;
				while (Vnew<Vexc){
				//	lnVnew = log(Vcurr) + (dsfmt_genrand_close_open(&dsfmt)-0.5)*vstep;
					//Vnew = exp(lnVnew);
					Vnew = Vcurr + (dsfmt_genrand_close_open(&dsfmt)-0.5)*vstep;
				}
			//	dV = Vnew - Vcurr;
				Vscale = Vnew/Vcurr;
				long double ranprob = dsfmt_genrand_close_open(&dsfmt);
//			//	long double vprob = 1.0/Vscale;
				long double vprob = 1.0;
				if (Vscale<1.0){
					long double N = (long double)(natoms);
					//N-=1.0;
					//N/=2.0;
					vprob = pow(Vscale, N);
					//vprob = exp(N*log(Vscale));
				}
				if (vprob>ranprob){
					//continue;
					
				

				//Tlocold.Equate(Tlocal);
					bscale = pow(Vscale, 1.0/3.0);
				
					blocal = bold*bscale;
				//rscale = dclocal/rold;
					Tlocal.SetBox(blocal, blocal, blocal);
				
					Tlocal.ScaleCoordinates(Vscale);
					long double Enew = TotalEnergy(Tlocal);
				//	cout<<"NS "<<ii<<" Will make v-def! Ecurr: "<<Et1local<<" Vcurr: "<<Vcurr<<" Enew: "<<Enew<<" Vnew: "<<Vnew<<endl;
				//	dE += Enew - Et1;
					Et1local = Enew;
					vdeftag = 1;
					step = bsfac*blocal*pow(0.90, cfailcount);
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
			Tlocal.WrapAtom(Arandom);
			//Check change in energy for pairwise interactions
			dE+=DeltaE(Tlocal, Arandom, origx, origy, origz);
			//Tlocal.CalcCom();
			
		
			++ctriesloc;
		}
		

		//Selection Criteria--Metropolis Criteria
		bool tagg = 1;
		
		
		//Add difference to current total potential energy
		Et1local+=dE;
		//Calculate the change in Enthalpy from trial moves
		//dH += dE + Press*dV;
		//cout << "dH: " << dH << " dE: " << dE << " dV: " << dV << endl;
		//Ht1 += dH;
		Ht1local = Et1local + Press*Vcurr;
		//cout<<"NS "<<ii<<" Ht1: "<<Ht1local<<" Et1: "<<Et1local<<" Vcurr: "<<Vcurr<<endl;
		//if (condition){
	//		
	//	}
		if(Ht1local < Hnest){
			//cout<<"NS "<<ii<<" Trial move was accepted! " << endl;
			tagg=0;
			Eelocal=Et1local;
			Helocal = Ht1local;
			Velocal = Vcurr;
			//cout<<"NS "<<ii<<" He: "<<Helocal<<" Ee: "<<Eelocal<<" Ve: "<<Velocal<<endl;
			//++acceptcount;
			if (Helocal<Hnmin[tid]){
				
				Hnmin[tid]=Helocal;
				Enmin[tid]=Eelocal;
				lgminbox[tid]=blocal;
				lgmin[tid].Equate(Tlocal);
				lgmin[tid].SetBox(blocal, blocal, blocal);
				//cout<<"local min for thread: "<<tid<<" setting!"<<endl;
				//cout<<"Hnmin: "<<Hnmin[tid]<<" Helocal "<<Helocal<<endl;
				//Vnmin[tid]=Velocal;
			}
			if (Helocal<Htlowprev){

				#pragma omp critical
				{
					Htlowprev=Helocal;
					Vlowprev=Velocal;
					blowprev=blocal;
					Elowprev=Eelocal;
					Bottom.Equate(Tlocal);
					Bottom.SetBox(blocal, blocal, blocal);
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
				Tlocal.ScaleCoordinates(Unscale);
				//Tlocal.Equate(Tlocold);
				blocal/=bscale;
				step = bsfac*blocal*pow(0.90, cfailcount);
				Vcurr = pow(blocal, 3);
				Tlocal.SetBox(blocal, blocal, blocal);
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
				}
			
		}


		//Output coordinates of sample n
	   if(n>Eqit && n % (ConfigsOutFreq) == 0) {
		coutfile << Tlocal.natoms << endl;
		coutfile <<"lj 17 H: "<<Ht1local<<" E: "<<Et1local<<" V: "<<Vcurr<<endl;
		for(unsigned long hg = 0; hg<Tlocal.natoms; ++hg){
			coutfile << setiosflags(ios::fixed) << setprecision(15) << "lj " << Tlocal.atom[hg].x << " " << Tlocal.atom[hg].y << " " << Tlocal.atom[hg].z << endl; 
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
	logfile<<"vfrac: "<<vafrac<<" cfrac: "<<cafrac<<" cfailcount: "<<cfailcount<<endl;
	logfile<<"vstep: "<<vstep<<endl;	
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
	return 0.0;
	long double Etot = 0.0;
	unsigned int natoms = c.natoms;
	long double box = c.boxx;
		for(unsigned int pp = 0;pp<natoms-1;++pp){

			for(unsigned int uu = pp+1;uu<natoms;++uu){
				
				
				Etot += PairEnergy(c.atom[pp], c.atom[uu], box);
					
			}
		}
//	Etot *= 60.0;
	return(Etot);
}

// potential energy function
long double PairEnergy(const Atom& one, const Atom& two, const long double& box){
	
	
	long double dx, dy, dz, rs;
	long double potent, eps, sigma;
	eps = 1.0, sigma = 1.0;
	long double boxx=box;
	long double boxy=box;
	long double boxz=box;
//	eps = sqrt(one.eps*two.eps), sigma = (one.sig+two.sig)/2.0;
//	double rc = 3.0*sigma;
	
	dx = one.x - two.x;
	
	dy = one.y - two.y;
	
	dz = one.z - two.z;
	//Minimum image -- coordinates must be pre-wrapped
	if(abs(dx)> boxx/2.0){
		dx = boxx - abs(one.x) - abs(two.x);
	}
	if(abs(dy)>boxy/2.0){
		dy = boxy - abs(one.y) - abs(two.y);
	}
	if(abs(dz)>boxz/2.0){
		dz = boxz - abs(one.z) - abs(two.z);
	}
	
	rs= (pow(dx, 2))+(pow(dy, 2))+(pow(dz, 2));
	
	potent = 4.0*eps*((pow(sigma/rs, 6))-(pow(sigma/rs, 3)));
	
	return(potent);
}

//Uses the o position values (xo, yo, zo) from atom one
long double PairEnergyO(const Atom& one, const long double& ox, const long double& oy, const long double& oz, const long double& box){
	
	
	long double dx, dy, dz, rs;
	long double potent, eps, sigma;
	long double boxx=box;
	long double boxy=box;
	long double boxz=box;
	eps = 1.0, sigma = 1.0;
//	double rc = 3.0*sigma;
	
	dx = one.x - ox;
	
	dy = one.y - oy;
	
	dz = one.z - oz;

	//Minimum image -- coordinates must be pre-wrapped 
	if(abs(dx)> boxx/2.0){
		dx = boxx - abs(one.x) - abs(ox);
	}
	
	if(abs(dy)>boxy/2.0){
		dy = boxy - abs(one.y) - abs(oy);
	}
	
	if(abs(dz)>boxz/2.0){
		dz = boxz - abs(one.z) - abs(oz);
	}
	
	rs= (pow(dx, 2))+(pow(dy, 2))+(pow(dz, 2));
	
	potent = 4.0*eps*((pow(sigma/rs, 6))-(pow(sigma/rs, 3)));
	
	return(potent);
}

// calculates the change in energy associated with one moved atom
long double DeltaE(const Config& c, const unsigned& Aindex, const long double& ox, const long double& oy, const long double& oz){
	return 0.0;	
	long double Etot=0.0;	
	long double Etoto=0.0;
	long double dE;
	long double box = c.boxx;
		
	for (unsigned int i = 0; i < c.natoms; ++i)
	{
		if(i!=Aindex){
			Etot += PairEnergy(c.atom[Aindex], c.atom[i], box);
			Etoto += PairEnergyO(c.atom[i], ox, oy, oz, box);
					
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
