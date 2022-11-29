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

4-10-14

	-Update to sample lj gas interacting with continuum solid surface;
		Surface in the center of box in xy plane. Volume deformations
		updated to just adjustment of box size in z coordinate
	Updated from v4.2 to v5.0 -- for sampling a different system
	-Added logfile output

4-23-14
	-Added another volume deformation in which the particle coordinates are not scaled
	Updated from v5.0 to v5.1 -- additional volume trial move

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

#include <boost/math/special_functions/spherical_harmonic.hpp>


using namespace std;

 // Classes
#include "/net/uu/nm/cm/bxw109120/NestedSampling/NSTAM/code/src/class/Atom.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/NSTAM/code/src/class/Config.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/NSTAM/code/src/class/RunningStats.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/NSTAM/code/src/class/HistogramNS_slim.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/NSTAM/code/src/class/Histogram.h"

//#include "./src/class/CenterofMass.h"


// prototype functions

long double PairEnergy(const Atom& one, const Atom& two, const long double& boxx, const long double& boxy, const long double& boxz);
long double PairEnergyO(const Atom& one, const long double& ox, const long double& oy, const long double& oz, const long double& boxx, const long double& boxy, const long double& boxz);
long double TotalEnergy(const Config& c);
long double DeltaE(const Config& c, const unsigned& Aindex, const long double& ox, const long double& oy, const long double& oz);
string itos( int Number );
long double RDist(const Atom& one, const Atom& two);
long double RDistMinImage(const Atom& one, const Atom& two, long double bx, long double by, long double bz);
long double LocalDensity(const Config& c);
long double TestArea(const Config& c, long double scale, long double Ec);
long double ClusterQ4(Config& c);

int main(int argc, char *argv[])
{

	//Get the local time for the logfile
	time_t rawtime;
 	struct tm * timeinfo;
 	time (&rawtime);
 	timeinfo = localtime (&rawtime);

	// Parameters that need to be set

	// the number of atoms to have in the system
	const unsigned int natoms = 38;
	
	// the number of nested sampling iterations to run 
	const unsigned int Iterations = 4;
	//Total number of MC style trial moves under each NS iteration
	const unsigned int Nsamples = 20000;
	//Number of procs to split Sampling loop (MCit)
	const unsigned short nprocs = 2;
	//number of nested energies back additional walkers should be from each other 
	int nskip = 2;
	// Number of initial MC steps to allow equilibration
	const unsigned Eqit = 0;
	// Step interval to collect data after equilibration
	const unsigned Cit = natoms*10;
	long double NSfrac = 0.5;
	// Swap probability
	const long double Pex = 0.10;
	//Box sizes
	long double boxx = 3.15;//*1.1;
	long double boxy = 3.15;//*1.1;
	long double boxz = 3.15;///(1.1*1.1);
	long double rcon = 2.25;
	// Test area scaling value
	long double scale = 1.00001;
	long double abox = (boxx+boxy+boxz)/3.0;
	//Translation step
	const long double bsfac = 0.2500;
 	
    long double stepx = boxx*bsfac;
	long double stepy = boxy*bsfac;
	long double stepz = boxz*bsfac;
	
	//Pi
	//const long double pi = 3.14159265359;
	
	//Enthalpy bounds
	const long double Etop = 20000.0;
	const long double Elow = -200.0;
	long double Binwidth = 1.0e-2;

		//run descriptor
	string nds = itos(natoms);
	string rundescript = "NS_comr_n"+nds+"_a5";
	//coordinates to read in
//	string infilename = "coord_n256_lowE.xyz";
//	ifstream infile;
//	infile.open(infilename.c_str());
//	// Returns an error if the file does not exist
//	if (!infile) {
//		cerr << "Unable to open file " << infilename.c_str() << endl;
//		exit(0);
//	}	
	
//------------------------------------------------------------
	
	//Name of logfile
	string logfilename = "oMMC_" + rundescript + ".log";
	//Name of Emedian output file
	string houtname = "Epoints_" + rundescript + ".dat";
	//Name of Configuration output file
	string coutname = "Coord_" + rundescript + ".xyz";
	//Name of Observable output file
	string outfilename = "zObservables_" + rundescript + ".dat";

	string qfoutname = "zNestDist_Q4_"+rundescript+".dat";
	ofstream qfout(qfoutname.c_str());
	ofstream logfile(logfilename.c_str());
	ofstream houtfile(houtname.c_str());
	ofstream obsoutfile(outfilename.c_str());
	//string coutfilename = "Coord_NS_Ent_a0.xyz";
	ofstream coutfile(coutname.c_str());
//	ofstream stepoutfile("zStepOut_int_test_a1.dat");	
	// configs from each nested sampling iteration
	 
//	unsigned long long MCitpp = MCit+Eqit;
//	unsigned int ConfigsOutFreq = ((MCitpp-Eqit))/2;

cout << "Nested sampling iterations is set to " << Iterations << " with NSamples set to " << Nsamples << endl;
cout<<"Etop: "<<Etop<<" Elow: "<<Elow<<endl;	

logfile << " " << endl;
	logfile << "Local time and date: " << asctime(timeinfo) << endl;
	logfile<<"Constant Volume Nested Sampling of lennard jones system with TAM"<<endl;
	logfile << "Parameters: " << endl;
	logfile<<"NS Iterations: "<<Iterations<<" Number of Samples to Collect (Nsamples): "<<Nsamples<<endl;
	logfile<<"NS fraction : "<<NSfrac<<endl;
	logfile<<"TAM scaling factor: "<<scale<<endl;
	logfile<<"natoms: "<<natoms<<" nprocs: "<<nprocs<<" Eqit: "<<Eqit<<" Cit: "<<Cit<<endl;
	logfile<<" boxx: "<<boxx<<" boxy: "<<boxy<<" boxz: "<<boxz<<endl;
	logfile<<"bsfac: "<<bsfac<<" initial guess stepx: "<<stepx<<" stepy: "<<stepy<<" "<<stepz<<endl;
	logfile<<" Etop: "<<Etop<<" Elow: "<<Elow<<endl;	
//	logfile<<"MCitpp: "<<MCitpp<<"ConfigsOutFreq: "<<ConfigsOutFreq<<endl;
	logfile<<"nskip: "<<nskip<<" Pex: "<<Pex<<endl;
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
	//Initialize histogram and running stat objects
	HistogramNS Ehist;
	Ehist.Initialize(Binwidth, Elow, Etop);
	Ehist.SetNSfrac(NSfrac);
	RunningStats Estat;
	RunningStats Esqstat;
	RunningStats LDstat;
	RunningStats dUstat;
	RunningStats dUsstat;
	RunningStats dUcstat;
//	RunningStats Rgstat;
	Histogram Qfhist;
	Qfhist.Initialize(Binwidth/10.0, -0.1, 1.0);
//Initialize configuration objects
	
Config trial1(natoms);
Config Bottom(natoms);
// local min confs
Config * lgmin  = new Config[nprocs];
// local min confs
Config * Tlocal  = new Config[nprocs];


Bottom.SetBox(boxx, boxy, boxz);

trial1.SetBox(boxx, boxy, boxz);


cout<<"Assigning Initial coordinates.."<<endl;
logfile<<"Assigning Initial coordinates.."<<endl;
//Randomly place atoms in box - but filter to configuration with Enthalpy below Htop

long double Econf;
bool rccheck=0;
do
{
	unsigned track=0;
	for (unsigned int i = 0; i < natoms; ++i)
	{
		long double rx = (dsfmt_genrand_close_open(&dsfmt)-0.5)*boxx;
		long double ry = (dsfmt_genrand_close_open(&dsfmt)-0.5)*boxy;
		long double rz = (dsfmt_genrand_close_open(&dsfmt)-0.5)*boxz;
		trial1.atom[i].x = rx;
		trial1.atom[i].y = ry;
		trial1.atom[i].z = rz;
		bool flag=0;
		//cout<<"assigned atom: "<<i<<endl;
		if (track<2*natoms){
			//cout<<"checking distance..."<<endl;
			for (unsigned int j = 0; j < i; ++j)
			{
				long double rdist = RDist(trial1.atom[i], trial1.atom[j]);
				//cout<<"distance from atom "<<i<<" to atom "<<j<<" is "<<rdist<<endl;
				if (rdist<1.0){
				//	cout<<"rdist less than closest allowed"<<endl;
					flag=1;
					break;
				}
			}
		}
		if (flag){
			--i;
			++track;
			//cout<<"reseting i back from "<<i+1<<" to "<<i<<" and track is "<<track<<endl;
		}
	}
	Econf = TotalEnergy(trial1);
	rccheck = trial1.CheckComDistAll(rcon);
	//cout<<Econf<<" "<<!rccheck<<endl;
	if (!rccheck){
		Econf=Etop+100.0;
	}
} while (Econf>Etop);

//Read coordinates from xyz file
//string junk;
//long double x,y,z;
//infile>>junk;
//for (unsigned int i = 0; i < 4; ++i)
//{
//	infile>>junk;
//}
//for (unsigned int i = 0; i < natoms; ++i)
//{
//	infile>>junk>>x>>y>>z;
//	trial1.atom[i].x=x;
//	trial1.atom[i].y=y;
//	trial1.atom[i].z=z;
//}

logfile<<"Initial coordinates assigned!"<<endl;
Bottom.Equate(trial1);

cout<<"Initial coordinates assigned!"<<endl;	
long double Egmin = TotalEnergy(trial1);
long double Enest = Etop;
long double Enestloc[nprocs];
long double Et1 = Egmin;
long double Enmin[nprocs];
for (unsigned int np = 0; np < nprocs; ++np)
{
	lgmin[np].InitializeAtoms(natoms);
	lgmin[np].Equate(Bottom);
	lgmin[np].SetBox(boxx, boxy, boxz);
	Tlocal[np].InitializeAtoms(natoms);
	Tlocal[np].SetBox(boxx, boxy, boxz);
	Enestloc[np]=Etop;
	Enmin[np]=Egmin;
}
// Output starting coordinates
coutfile << trial1.natoms << endl;
coutfile <<"lj "<<natoms<<" E: "<<Egmin<<endl;
for(unsigned long hg = 0; hg<trial1.natoms; ++hg){
		coutfile << setiosflags(ios::fixed) << setprecision(15) << "lj " << Bottom.atom[hg].x << " " << Bottom.atom[hg].y << " " << Bottom.atom[hg].z << endl; 
		}

long double Elowprev = Et1;
unsigned long long cfailcount = 0;
long double Enprev=Enest;

for(unsigned int ii=0;ii<Iterations;++ii){
cout << "Iteration " << ii << "......." << endl;
logfile << "Iteration " << ii << "......." << endl;
 if(ii>0){
 
Enest=Ehist.NestedEnergy();

string oname = "zNestDist_Q4_"+rundescript+"_iter_"+itos(ii-1)+".dat";
Qfhist.Normalize();
Qfhist.DumpProbability(oname);
Qfhist.DumpNorm(qfout);
Qfhist.Reset();
//Hhist.Normalize();
//Hnest=Hhist.Median();
long double Eavg = Estat.Mean();
long double Esqavg = Esqstat.Mean();
long double LDavg = LDstat.Mean();
long double dUavg = dUstat.Mean();
long double dUsavg = dUsstat.Mean();
long double dUcavg = dUcstat.Mean();
//long double Rgavg = Rgstat.Mean();
cout<<"Nested Energy is: "<<Enest<<" and Average Energy is "<<Eavg<<" and Average Energy Squared: "<<Esqavg<<endl;
logfile<<"Nested Energy is: "<<Enest<<" and Average Energy is "<<Eavg<<" and Average Energy Squared: "<<Esqavg<<endl;
cout<<"Averag Local Density: "<<LDavg<<endl;
logfile<<"Average Local Density: "<<LDavg<<endl;
cout<<"Average dU: "<<dUavg<<endl;
logfile<<"Average dU: "<<dUavg<<endl;
obsoutfile<<Eavg<<" "<<Esqavg<<" "<<LDavg<<" "<<dUavg<<" "<<dUsavg<<" "<<dUcavg<<endl;
Ehist.SmartReset();
//Hhist.Reset();
Estat.Reset();
Esqstat.Reset();
LDstat.Reset();
dUstat.Reset();
dUsstat.Reset();
dUcstat.Reset();
//Rgstat.Reset();
trial1.Equate(Bottom);
Et1=Elowprev;

//trial1.SetBox(boxx, boxy, boxz);
abox = (boxx+boxy+boxz)/3.0;

stepx = bsfac*boxx*pow(0.90, cfailcount);
stepy = bsfac*boxy*pow(0.90, cfailcount);
stepz = bsfac*boxz*pow(0.90, cfailcount);
if (abs(Enest-Enprev)<0.0001){
	break;
}
long double Nsk = 1.0+4.0*exp(-abs(Enest-Enprev));
cout<<"Nsk: "<<Nsk<<endl;
Enprev=Enest;
if (ii%nskip==0){
	

// set so that index 0 is the current target median
long double Enpb[nprocs];
		for (unsigned int i = 1; i < nprocs; ++i)
		{
			Enpb[i]=Enestloc[i-1];
		}
		Enpb[0]=Enest;
		for (unsigned int i = 0; i < nprocs; ++i)
		{
			Enestloc[i]=Enpb[i];
		}
}
else{
	Enestloc[0]=Enest;
}

//		for (unsigned int i = 0; i < nprocs; ++i)
//		{
//			cout<<"walker "<<i<<" Enest[i]: "<<Enestloc[i]<<endl;
//		}
	//Make sure starting conf is lower enthalpy than new nested value - if not set to global found minima
for (unsigned int np = 0; np < nprocs; ++np)
{
	cout<<"Walker "<<np<<" Enest "<<Enestloc[np]<<" Enmin : "<<Enmin[np]<<endl;
	if (Enmin[np]>Enestloc[np]){
		lgmin[np].Equate(Bottom);
		Enmin[np]=Elowprev;
	
	}
}
}
houtfile<<Enest<<endl;

cout<<"Initial Energy: "<<Et1<<endl;
logfile<<"Initial Energy: "<<Et1<<endl;

	//unsigned long long acceptcount = 0;
	
	unsigned long long ctries = 0;
	unsigned long long csuccess = 0;
	//long double Elowprev = Et1;
	unsigned long long nsamp=0;
	long double Et1local[nprocs];
	long double Et1b[nprocs];
	omp_set_num_threads(nprocs);
	#pragma omp parallel 
	{
	unsigned tid = omp_get_thread_num();
	 //Define a thread local configuration
	// Config Tlocal(natoms);
	 Tlocal[tid].Equate(lgmin[tid]);
	 //Tlocal.SetBox(boxx, boxy, boxz);
	// Config Tlocold(natoms);
	bool lgflag=1;
	//Define thread local variables
	 Et1local[tid] = Enmin[tid];
	
	long double Eelocal = 0.0;
		 
	unsigned long long ctriesloc = 0;
	unsigned long long csuccessloc = 0;

	unsigned long long n=0;
  //Perform markov moves on the configuration
    while(nsamp<Nsamples){
//		if(n%10==0){
//			#pragma omp critical
//			{
//			cout<<"Thread "<<tid<<" Trial move "<<n<<endl;
//			}
//		}
		#pragma omp flush
		//Store Initial Et1 value
		Et1b[tid] = Et1local[tid];
		long double ran[5];
		#pragma omp critical(rangen)
		{
			//cout<<"tid "<<tid<<" generating random numbers"<<endl;
			for (unsigned int r = 0; r < 5; ++r)
			{
				ran[r]=dsfmt_genrand_close_open(&dsfmt);
			}
			

		}
		//move tag - 0=swap, 1=translation
		bool mtag = 0;
		//vars needed for translation
		long double nx, ny, nz;
		unsigned int Arandom = (unsigned int)(ran[0]*natoms);
		//attemp swap with tid+1
		//bool eflag=1;
		//no exchange do translation
		//if(eflag){
		//Trial Move -----
		// randomly select an atom index	
//		#pragma omp critical
//		{
//		cout<<"tid "<<tid<<" attempting translation"<<endl;
//		}
			mtag=1;
			
		
			
			nx = Tlocal[tid].atom[Arandom].x;
			ny = Tlocal[tid].atom[Arandom].y;
			nz = Tlocal[tid].atom[Arandom].z;
	
		//Coordinate change		

			
			long double movex, movey, movez;
			/* set random movements of the coordinates
			   of the selected atom */
			movex = (ran[1]-0.5)*stepx;
			movey =(ran[2]-0.5)*stepy;
			movez = (ran[3]-0.5)*stepz;
			// Implement those moves
			Tlocal[tid].atom[Arandom].x+=movex;
			Tlocal[tid].atom[Arandom].y+=movey;
			Tlocal[tid].atom[Arandom].z+=movez;
			//Tlocal[tid].CalcCom();
			bool dcheck= Tlocal[tid].CheckComDistAll(rcon);
			long double dEp=0.0;
			if (dcheck){
//				#pragma omp critical
//				{
//				cout<<"tid "<<tid<<" boundary success. Compute dE"<<endl;
//				}
				//Check change in energy for pairwise interactions
				dEp=DeltaE(Tlocal[tid], Arandom, nx, ny, nz);
			}
			else{
				dEp=1.0e10;

			}
		//	cout<<"dcheck "<<dcheck<<" dEp "<<dEp<<endl;
		
		
			++ctriesloc;
		
			//Selection Criteria--Metropolis Criteria
			//bool ttag = 1;
		
		
			//Add difference to current total potential energy
			Et1local[tid]+=dEp;
//				#pragma omp critical
//			{
//			cout<<"tid "<<tid<<"  translation complete"<<endl;
//			}
		//}
		
		if(Et1local[tid] < Enestloc[tid]){
			
			if (mtag){
//				Tlocal[tid].atom[Arandom].x=nx;
//				Tlocal[tid].atom[Arandom].y=ny;
//				Tlocal[tid].atom[Arandom].z=nz;
				++csuccessloc;
			}
			//tagg=0;
		//	Eelocal=Et1local;
//			Tlocal.CalcCom();
//			Tlocal.ReCenterZero();
			#pragma omp flush(Et1local, Elowprev)
			if (Et1local[tid]<Elowprev){
				
				#pragma omp critical(bottom)
				{
					
					Elowprev=Et1local[tid];
					Bottom.Equate(Tlocal[tid]);
					//Bottom.SetBox(boxx, boxy, blocal);
					#pragma omp flush(Bottom, Elowprev)
				}
				
			}
			#pragma omp flush(lgflag, Enmin)
			if (Et1local[tid]<Enmin[tid]){
				lgmin[tid].Equate(Tlocal[tid]);
				Enmin[tid]=Et1local[tid];
				//lgflag=0;
			}
							
			
		}
		else{
			
			Et1local[tid] = Et1b[tid];
			if (mtag){
				Tlocal[tid].atom[Arandom].x=nx;
				Tlocal[tid].atom[Arandom].y=ny;
				Tlocal[tid].atom[Arandom].z=nz;
			}
			
		}

			if(n>Eqit && Et1local[tid]<Enest && n%(Cit)==0){
//					#pragma omp critical
//				{
//				cout<<"tid "<<tid<<" pushing values"<<endl;
//				}
			//	cout << "!!!!!Going to push values!!!!.................--------------------------!!!!!" << endl;
			//	cout << "Ee: " << Eelocal << " He: " << Helocal << " Ve: " << Velocal << endl;	
			//	cout << "He*He: " << He*He << endl;	
			//	long double LD = LocalDensity(Tlocal[tid]);
			//	long double dU = TestArea(Tlocal[tid], scale, Eelocal);
				long double LD=ClusterQ4(Tlocal[tid]);
				long double dU = 0.0;
				
				#pragma omp critical(pushval)
				{
					#pragma omp flush(Ehist, Estat, Esqstat, nsamp)
					Ehist.Push(Et1local[tid]);
					Estat.Push(Et1local[tid]);
					Esqstat.Push(0.0);
					//if (LD<0.3){
						LDstat.Push(LD);
					//}
					Qfhist.Push(LD);
					dUstat.Push(dU);
					dUsstat.Push(dU*dU);
					dUcstat.Push(dU*dU*dU);
					
					++nsamp;
					#pragma omp flush(Ehist, Estat, Esqstat, nsamp)
				}
				//#pragma omp flush
		}
//		if (n%10000==0){
//			Tlocal[tid].CalcCom();
//			Tlocal[tid].ReCenterZero();
//		}

//		//Output coordinates of sample n
//	   if(n>Eqit && (n-Eqit)%(ConfigsOutFreq)==0) {
//			#pragma omp critical 
//			{
//				//cout<<"Output coordinates to file at trial move: "<<n<<endl;
//				coutfile << Tlocal.natoms << endl;
//				coutfile <<"lj "<<Tlocal.natoms<<" E: "<<Et1local<<endl;
//				for(unsigned long hg = 0; hg<Tlocal.natoms; ++hg){
//					coutfile << setiosflags(ios::fixed) << setprecision(15) << "lj " << Tlocal.atom[hg].x << " " << Tlocal.atom[hg].y << " " << Tlocal.atom[hg].z << endl; 
//				}
//			}
//		//cout<<" -----------------------------"<<endl;

//	   }	

	

	//	cout << endl;
	
		
	

	++n;

		//check for exchange
		//synchronise threads
		#pragma omp barrier
		//critical only one thread at a time
		#pragma omp critical(exchange)
		{
		//cout<<"tid "<<tid<<" checking for exchange"<<endl;
		if (tid<nprocs-1 && ran[4]<Pex){
			//cout<<"tid "<<tid<<" attempting exchange"<<endl;
			//eflag=0;
			//#pragma omp critical(exchange)
			//{
				#pragma omp flush(Tlocal, Et1local, Et1b)
				
				//long double Etidp1=TotalEnergy(Tlocal[tid+1]);
				  long double Etidp1 = Et1local[tid+1];
				if (Etidp1<Enestloc[tid]){
					//store temp of 
					Config Swap(Tlocal[tid+1]);
					//replace 
					Tlocal[tid+1].Equate(Tlocal[tid]);
					// replace with temp
					Tlocal[tid].Equate(Swap);
					
					Et1local[tid+1]=Et1local[tid];
					//Et1b[tid+1]=Et1local[tid+1];
					Et1local[tid]=Etidp1;
					//cout<<"tid "<<tid<<" exchange successful"<<endl;
				}
				#pragma omp flush(Tlocal, Et1local, Et1b)
			//}
			
		}
		}
		//make sure synchronised again before continuing
		#pragma omp barrier
		
    } //End Monte Carlo loop

	#pragma omp barrier
	
	#pragma omp critical(trialadd)
	{ 
		
		csuccess+=csuccessloc;
		ctries+=ctriesloc;

	}


}//End parallel

	
	long double cafrac = (long double)(csuccess)/(long double)(ctries);
	
	if (cafrac<0.30){
		++cfailcount;
	}
	cout<<" cfrac: "<<cafrac<<" cfailcount: "<<cfailcount<<endl;
	logfile<<" cfrac: "<<cafrac<<" cfailcount: "<<cfailcount<<endl;
	logfile<<" stepx: "<<stepx<<" stepy: "<<stepy<<" stepz: "<<stepz<<endl;
		
//		coutfile << trial1.natoms << endl;
//		coutfile <<"lj 17 H: "<<Htlowprev<<" E: "<<Elowprev<<" V: "<<Vlowprev<<endl;
//		for(unsigned long hg = 0; hg<trial1.natoms; ++hg){
//			coutfile << setiosflags(ios::fixed) << setprecision(15) << "lj " << Bottom.atom[hg].x << " " << Bottom.atom[hg].y << " " << Bottom.atom[hg].z << endl; 
//		}
//	if (cafrac<0.10){
//		cout<<"config frozen. exiting!"<<endl;
//		exit(0);
//	}

} // End NS loop

logfile.close();
houtfile.close();
obsoutfile.close();
coutfile.close();
//stepoutfile.close();
logfile<<"Simulation Complete!"<<endl;
logfile<<"Output files are: "<<endl;
logfile<<"Nested Energies: "<<houtname<<endl;
logfile<<"Configuration Samples: "<<coutname<<endl;
logfile<<"Average Observable Values: "<<outfilename<<endl;
cout<<"Simulation Complete!"<<endl;
cout<<"Output files are: "<<endl;
cout<<"Logfile: "<<logfilename<<endl;
cout<<"Nested Energies: "<<houtname<<endl;
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
	long double bx, by, bz;
	bx = c.boxx, by = c.boxy, bz = c.boxz;
		for(unsigned int pp = 0;pp<natoms-1;++pp){

			for(unsigned int uu = pp+1;uu<natoms;++uu){
				
				
				Etot += PairEnergy(c.atom[pp], c.atom[uu], bx, by, bz);
					
			}
		}
//	Etot *= 60.0;
	return(Etot);
}

// potential energy function
long double PairEnergy(const Atom& one, const Atom& two, const long double& boxx, const long double& boxy, const long double& boxz){
	
	
	long double dx, dy, dz, rs;
	long double potent, eps, sigma;
	eps = 1.0, sigma = 1.0;
	long double bx=boxx;
	long double by=boxy;
	long double bz=boxz;
//	eps = sqrt(one.eps*two.eps), sigma = (one.sig+two.sig)/2.0;
	long double rc = 3.0*sigma;
	long double rcs = rc*rc;
	dx = one.x - two.x;
	
	dy = one.y - two.y;
	
	dz = one.z - two.z;
	rs= (pow(dx, 2))+(pow(dy, 2))+(pow(dz, 2));
	potent = 4.0*eps*((pow(sigma/rs, 6))-(pow(sigma/rs, 3)));

	
	return(potent);
}

//Uses the o position values (xo, yo, zo) from atom one
long double PairEnergyO(const Atom& one, const long double& ox, const long double& oy, const long double& oz, const long double& boxx, const long double& boxy, const long double& boxz){
	
	
	long double dx, dy, dz, rs;
	long double potent, eps, sigma;
	long double bx=boxx;
	long double by=boxy;
	long double bz=boxz;
	eps = 1.0, sigma = 1.0;
	long double rc = 3.0*sigma;
	long double rcs = rc*rc;
	dx = one.x - ox;
	
	dy = one.y - oy;
	
	dz = one.z - oz;



	rs= (pow(dx, 2))+(pow(dy, 2))+(pow(dz, 2));
	
		potent = 4.0*eps*((pow(sigma/rs, 6))-(pow(sigma/rs, 3)));
	

	return(potent);
}

// calculates the change in energy associated with one moved atom
long double DeltaE(const Config& c, const unsigned& Aindex, const long double& ox, const long double& oy, const long double& oz){
	long double Etot=0.0;	
	long double Etoto=0.0;
	long double dE;
	long double bx, by, bz;
	bx = c.boxx, by = c.boxy, bz = c.boxz;
		
	
	for (unsigned int i = 0; i < c.natoms; ++i)
	{
		if(i!=Aindex){
			Etot += PairEnergy(c.atom[Aindex], c.atom[i], bx, by, bz);
			Etoto += PairEnergyO(c.atom[i], ox, oy, oz, bx, by, bz);
					
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

//Calculates the distance between two atoms--
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
	//cout<<"dx: "<<dx<<" dy: "<<dy<<" dz: "<<dz<<" r "<<r<<endl;
	return(r);



}
long double RDistMinImage(const Atom& one, const Atom& two, long double bx, long double by, long double bz){

	long double dx, dy, dz, r, rs;
	long double x1, y1, z1, x2, y2, z2;
	x1 = one.x, y1=one.y, z1=one.z;
	x2=two.x, y2=two.y, z2=two.z;
	
	dx = x1-x2;
	
	dy = y1-y2;
	
	dz = z1-z2;
	//Minimum image -- coordinates must be pre-wrapped
	if(abs(dx)> bx/2.0){
		dx = bx - abs(one.x) - abs(two.x);
	}
	if(abs(dy)>by/2.0){
		dy = by - abs(one.y) - abs(two.y);
	}
	if(abs(dz)>bz/2.0){
		dz = bz - abs(one.z) - abs(two.z);
	}
	rs= (pow(dx, 2))+(pow(dy, 2))+(pow(dz, 2));
	r = sqrt(rs);
	//cout<<"dx: "<<dx<<" dy: "<<dy<<" dz: "<<dz<<" r "<<r<<endl;
	return(r);



}

// Total potential energy calculator, uses only the lj pot
long double LocalDensity(const Config& c){
	
	long double r;
	unsigned int natoms = c.natoms;
	long double bx, by, bz;
	bx = c.boxx, by = c.boxy, bz = c.boxz;
	long double cutoff = 2.0;
	long double volume = (4.0/3.0)*3.14*pow(cutoff, 3);
	long double avgden = 0.0;
	long double dc = 0.0;
		for(unsigned int pp = 0;pp<natoms-1;++pp){
			long double count = 0.0;
			for(unsigned int uu = pp+1;uu<natoms;++uu){
				
				
				r= RDistMinImage(c.atom[pp], c.atom[uu], bx, by, bz);
				if (r<cutoff){
					count+=1.0;
				}	
			}
		 long double dense = count/volume;
		 avgden +=dense;
		dc+=1.0;
		}
	long double finavg = avgden/dc;
//	Etot *= 60.0;
	return(finavg);
}

long double TestArea(const Config& c, long double scale, long double Ec){
	Config Test(c);
//	long double scale = 1.0001;
	Test.ScaleCoordinates(scale, scale, 1.0/(scale*scale));
	Test.boxx*=scale;
	Test.boxy*=scale;
	Test.boxz/=(scale*scale);
	//long double Ec = TotalEnergy(c);
	long double Et = TotalEnergy(Test);
	long double dE = Et-Ec;
	return dE;

}

long double ClusterQ4(Config& c){
	c.CalcCom();
	long double rmin=1e10;
	int mina;
	//find the atom closest to COM
	for (unsigned int i = 0; i < c.natoms; ++i)
	{
		long double dist = c.ComDist(i);
			if (dist<rmin){
				rmin=dist;
				mina=i;
			}
	}
	//cout<<"mina "<<mina<<endl;
	//now find the closest 12 neighbors
	int neighbors[12];
	long double minprev = 1e10;
	for (unsigned int i = 0; i < 12; ++i)
	{
		long double min = 1e10;
		
		for (unsigned int j = 0; j < c.natoms; ++j)
		{
			if (i==0 && j!=mina){
				long double rdist = RDist(c.atom[mina], c.atom[j]);
				if (rdist<min){
					min=rdist;
					neighbors[i]=j;
					//cout<<"neighbors["<<i<<"] "<<neighbors[i]<<" j "<<j<<endl;
				}
			}
			else if(i>0 && j!=mina){
				long double rdist = RDist(c.atom[mina], c.atom[j]);
				//cout<<"rdist: "<<rdist<<" min "<<min<<" minprev "<<minprev<<endl;
				if (rdist<min && rdist>minprev){
					min=rdist;
					neighbors[i]=j;
					//cout<<"neighbors["<<i<<"] "<<neighbors[i]<<" j "<<j<<endl;
				}
			}
		}
		minprev=min;
	}
	//define the core cluster
	Config Core(13);
	Core.atom[0].Equate(c.atom[mina]);
	for (unsigned int i = 1; i < 13; ++i)
	{
		int ind = neighbors[i-1];
		//cout<<"ind "<<ind<<" i-1 "<<i-1<<endl;
		Core.atom[i].Equate(c.atom[ind]);
	}
	
	Core.CalcCom();
	//Identify bonds and calculate their angles
	int maxbonds = (13*13-13)/2;
	int nbonds=0;
	long double angles[maxbonds][2];
	for (unsigned int i = 0; i < Core.natoms-1; ++i)
	{
		for (unsigned int j = i+1; j < Core.natoms; ++j)
		{
			long double dist = RDist(Core.atom[i], Core.atom[j]);
			if(dist<1.39){
				++nbonds;
				Config two(2);
				two.atom[0].Equate(Core.atom[i]);
				two.atom[1].Equate(Core.atom[j]);
				two.CalcCom();
				long double vector[3];
				vector[0]=two.COM[0]-Core.COM[0];
				vector[1]=two.COM[1]-Core.COM[1];
				vector[2]=two.COM[2]-Core.COM[2];
				long double Rsq = vector[0]*vector[0]+vector[1]*vector[1];
				long double rhosq = Rsq + vector[2]*vector[2];
				long double R = sqrt(Rsq);
				long double rho = sqrt(rhosq);
				long double theta = acos(vector[0]/R);
				long double phi = asin(vector[2]/rho);
				//cout<<"theta "<<theta<<" phi "<<phi<<endl;
//				if (theta!=theta){
//					cout<<"vector[0] "<<vector[0]<<" R "<<R<<endl;
//					theta=3.14/2.0;
//				}
//				if (phi!=phi){
//					//cout<<"vector[1] "<<vector[1]<<" rho "<<rho<<endl;
//					phi=3.14/2.0;
//				}
				angles[nbonds-1][0]=theta;
				angles[nbonds-1][1]=phi;
			}
		}
	}
//	cout<<"nbonds "<<nbonds<<endl;
	if (nbonds>0){
		
	
	//now get Q4 sum
	long double sumo=0.0;
	for (int i = -4; i <5 ; ++i)
	{
		long double sumin = 0.0;
		for (unsigned int j = 0; j < nbonds; ++j)
		{
			long double val ;
		   val = boost::math::spherical_harmonic_r(4, i, angles[j][0], angles[j][1]);
			//cout<<"val "<<endl;
			if (val!=val){
				cout<<"val "<<val<<endl;
				exit(0);
			}
			sumin+=val;
 
		}
		//cout<<"sumin "<<sumin<<endl;
		sumin/= (long double)(nbonds);
		sumo+=sumin*sumin;
	}
	//cout<<"sumo "<<sumo<<endl;
	long double Q4 = sqrt((4.0*3.14/9.0)*sumo);
	//cout<<"Q4 "<<Q4<<endl;
//	if (Q4!=Q4){
//		Q4=1.0;
//	}
	//cout<<"nbonds "<<nbonds<<endl;
	//exit(0);
	//cout<<"returning Q4 "<<Q4<<endl;
	return Q4;
	}
	else{
		return 1.0;
	}
}
