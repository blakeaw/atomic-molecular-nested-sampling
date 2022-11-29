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
#include "/net/uu/nm/cm/bxw109120/src/SimulationCpp/class/Atom.h"
#include "/net/uu/nm/cm/bxw109120/src/SimulationCpp/class/ConfigDynamic.h"
#include "/net/uu/nm/cm/bxw109120/src/SimulationCpp/class/RunningStats.h"
#include "/net/uu/nm/cm/bxw109120/src/SimulationCpp/class/HistogramNS_slim.h"
#include "/net/uu/nm/cm/bxw109120/src/SimulationCpp/class/HardOpenCylinder.h"
//#include "./src/class/Histogram.h"

//#include "./src/class/CenterofMass.h"


// prototype functions

long double PairEnergy(const Atom& one, const Atom& two, const long double& boxx);
long double PairEnergyO(const Atom& one, const long double& ox, const long double& oy, const long double& oz, const long double& boxx);
long double TotalEnergy(const Config& c);
long double DeltaE(const Config& c, const unsigned& Aindex, const long double& ox, const long double& oy, const long double& oz);
string itos( int Number );
long double DeltaEInsert(const Config& c, const long double& ix, const long double& iy, const long double& iz);
long double DeltaEDelete(const Config& c, const unsigned& Dindex);
long double ConfProb(long double E, long double Beta);

int main(int argc, char *argv[])
{

	//Get the local time for the logfile
	time_t rawtime;
 	struct tm * timeinfo;
 	time (&rawtime);
 	timeinfo = localtime (&rawtime);

	// Parameters that need to be set

	// the number of atoms to have in the system initially
	const unsigned int natoms = 17;
	

	//Total number of MC style trial moves under each NS iteration
	const unsigned int MCit = 2000000;
	//Number of procs to split Sampling loop (MCit)
	const unsigned short nprocs = 1;
	// Number of initial MC steps to allow equilibration
	const unsigned Eqit = 4000000;
	//Step Interval to attempt insertion/deletion move
	const short Vint = 10;
	// Step interval to collect data after equilibration
	const unsigned Cit = 500;
	//Boltzmann constant
	const long double kb = 1.0;
	const long double pi = 3.14159265359;
	//Chemical Potential
	const long double ChemPo = 0.010;
	//Temperature
	const long double Temp = 0.4;
	//Hard Open Cylinder Params
	const long double Length = 10.0;
	const long double Radius = 5.0;
	//Nested Histogram Params
	const long double Htop = 100.0;
	const long double Hlow = -2000.0;
	const long double binwidth = 1.0e-3;
	//spherical volume prefactor
	const long double vfact = (4.0/3.0)*pi;
	/*Spherical COM boundary condition radius */
	const long double rdist = 2.60;
	//Planck reduced - material specific
		/* argon - 0.185
		neon - 0.563
		krypton - 0.104
		xenon - 0.065
		---Ref: Quasi Harmonic Lattice Dynamics and Molecular Dynamics calculations
		for the Lennard-Jones solids
		*/
	const long double hred = 0.185;
	//Translation step 
 	long double stepguess= 0.4e-1*rdist;
    long double step = stepguess; 

	//run descriptor
	string nds = itos(natoms);
	string rundescript = "MC_Ent_n"+nds+"_comr_b0_T0.4";
	//Name of file to read in initial coordinates
	string filename = "LJ17_GlobalMin.txt";	
//------------------------------------------------------------
	cout<<"rdist: "<<rdist<<" Htop: "<<Htop<<" Hlow: "<<Hlow<<" binwidth: "<<binwidth<<endl;
	unsigned long long MCitpp = (MCit+Eqit*nprocs)/nprocs;
	unsigned int ConfigsOutFreq = ((MCitpp)-Eqit)/20;

	//Name of logfile
	string logfilename = "oMMC_" + rundescript + ".log";
	//Name of Emedian output file
	//string houtname = "Hpoints_" + rundescript + ".dat";
	//Name of Configuration output file
	string coutname = "Coord_" + rundescript + ".xyz";
	//Name of Observable output file
	string outfilename = "zObservables_" + rundescript + ".dat";

	
	ofstream logfile(logfilename.c_str());
	//ofstream houtfile(houtname.c_str());
	ofstream obsoutfile(outfilename.c_str());
	//string coutfilename = "Coord_NS_Ent_a0.xyz";
	ofstream coutfile(coutname.c_str());

	logfile << " " << endl;
	logfile << "Local time and date: " << asctime(timeinfo) << endl;
	logfile<<"Constant Chemical Potential MMMC Simulation of LJ particles with a hard spherical boundary centered at the configuration COM."<<endl;
	logfile << "Parameters: " << endl;
	logfile<<" Number of Trial Moves (MCit): "<<MCit<<endl;
	logfile<<"natoms: "<<natoms<<" nprocs: "<<nprocs<<" Eqit: "<<Eqit<<" Cit: "<<Cit<<endl;
	logfile<<"ChemPo: "<<ChemPo<<endl;
	logfile<<" initial guess: "<<stepguess<<endl;
	logfile<<" Htop: "<<Htop<<" Hlow: "<<Hlow<<" binwidth: "<<binwidth<<endl;	
	logfile<<"MCitpp: "<<MCitpp<<"ConfigsOutFreq: "<<ConfigsOutFreq<<endl;
	logfile<<"Cylinder Parms:"<<endl;
	logfile<<"Length: "<<Length<<" Radius: "<<Radius<<endl;
//	ofstream stepoutfile("zStepOut_int_test_a1.dat");	
	// configs from each nested sampling iteration
	 
	long double Beta = 1.0/Temp;
	long double box = Length;
cout <<" MC iterations set to " << MCit << endl;
cout<<"Htop: "<<Htop<<" Hlow: "<<Hlow<<endl;	
cout<<"Temperature: "<<Temp<<endl;

	// Initialize the MT RNG object
	dsfmt_t dsfmt;
	
	// Need to initialize member mti as NN+1
	//Random1.mti = NN+1;
	//thermal wavelength
	long double tw = sqrt( pow(hred, 2)/(2.0*Temp*pi));
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
	//Hhist.Initialize(binwidth, Hlow, Htop);
	RunningStats Estat;
	RunningStats Esqstat;
	RunningStats Nstat;
	RunningStats Hstat;
	RunningStats Rgstat;
//Initialize configuration objects
	
Config trial1(natoms);
Config Bottom(natoms);
HardOpenCylinder boundary(Length, Radius);
// local min confs
Config * lgmin  = new Config[nprocs];

// Initialize an array to store the lgmin Energy values
long double Hnmin[nprocs]={Htop};
long double Enmin[nprocs]={0.0};
// store lgmin n 
int lgminn[nprocs]={0};
ifstream infile;
infile.open(filename.c_str());
// Returns an error if the file does not exist
if (!infile) {
	cerr << "Unable to open file " << filename << endl;
	cout << "Unable to open file " << filename.c_str() << endl;
	exit(0);
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
long double srg = Bottom.RadGyr();
cout<<"Initial Radius of Gyration: "<<srg<<endl;
Bottom.SetBox(box, box, box);
//if (!Bottom.CheckComDistAll(rdist)){
//	cout<<"Radius is too small for initial coordinates."<<endl;
//	cout<<"exiting..."<<endl;
//	exit(1);
//}
trial1.Equate(Bottom);
trial1.SetBox(box, box, box);
//step=bsfac*rad;	
long double Egmin = TotalEnergy(trial1);
cout<<"Egmin "<<Egmin<<endl;
long double Ninitial = (long double)(natoms);
long double Hstart = Egmin - ChemPo*Ninitial;
long double Et1 = Egmin;
long double Ht1 = Hstart;
cout<<"Einitial "<<Egmin<<" Ninitial "<<Ninitial<<" Hstart "<<Hstart<<endl;
//exit(0);

	
long double V = boundary.ComputeVolume();








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
	//	unsigned tid = 0;
	 //Define a thread local configuration
	 Config Tlocal;
	Tlocal.SetBox(box, box, box);
	long double Et1local;
	long double Ht1local;
	int nlocal;

		nlocal = Bottom.natoms;
		Tlocal.InitializeAtoms(nlocal);
	 	Tlocal.Equate(Bottom);
		//Tlocal.SetBox(blowprev, blowprev, blowprev);
		Et1local=Egmin;
		Ht1local=Hstart;
		

	#pragma omp critical
	{
		cout<<"thread walker: "<<tid<<" Initial Enthalpy: "<<Ht1local<<" Initial Energy: "<<Et1local<<" nlocal: "<<nlocal<<endl;
		logfile<<"thread walker: "<<tid<<" Initial Enthalpy: "<<Ht1local<<" Initial Energy: "<<Et1local<<" Initial N: "<<nlocal<<endl;
	}
//	Et1local=TotalEnergy(Tlocal);
//	Ht1local=Et1local - ( (long double)(Tlocal.natoms))*ChemPo;
//	cout<<"Starting values. E: "<<Et1local<<" H "<<Ht1local<<endl;
	//Define thread local variables
	//long double Et1local = Et1;
	long double Helocal = 0.0;
	long double Nelocal = 0.0;
	long double Eelocal = 0.0;
	//long double Ht1local = Ht1;
	//long double rlocal = rad;
	 
	unsigned long long vtriesloc = 0;
	unsigned long long vsuccessloc = 0;
	unsigned long long ctriesloc = 0;
	unsigned long long csuccessloc = 0;
	
  //Perform markov moves on the configuration
    for(unsigned int n = 0; n<MCitpp; ++n){
		//cout<<endl;
		//cout<<"trial move "<<n<<endl;
		//Store Initial Et1 value
		long double Et1b = Et1local;
		long double Ht1b = Ht1local;
		//Trial Move -----
		// randomly select an atom index	
			unsigned int Arandom = (unsigned int)(dsfmt_genrand_close_open(&dsfmt)*nlocal);
		
			long double origx, origy, origz;
			origx = Tlocal.atom[Arandom].x;
			origy = Tlocal.atom[Arandom].y;
			origz = Tlocal.atom[Arandom].z;
			long double nx, ny, nz;
		// Calculate the difference in total potential energy
		//long double dH = 0.0;
		long double dE = 0.0;
	
		short vdeftag = 0;
		bool tflag = 0;
		long double idp = 1.0;
		//Volume deformation move
		if(n%Vint==0){
			//cout<<" insertion/deletion attempt"<<endl;
			++vtriesloc;
			long double rp = dsfmt_genrand_close_open(&dsfmt);
			if (rp<0.5){
				
				//particle insertion	
			//	rp = dsfmt_genrand_close_open(&dsfmt);
				long double vp = V / ( (long double)(nlocal+1)*pow(tw, 3) );
				//if (rp<vp){
					idp = vp;
					//cout<<"particle insertion"<<endl;
					vdeftag = 1;
					//pick coords for ghost atom
				
					long double u, v;
					u = dsfmt_genrand_close_open(&dsfmt);
					v = dsfmt_genrand_close_open(&dsfmt);		
					long double theta = u*2.0*pi;
//					long double phi = acos((2.0*v) - 1.0);
//					long double Rho = rdist*dsfmt_genrand_close_open(&dsfmt);
//					long double dz = Rho*cos(phi);
					long double r = Radius*dsfmt_genrand_close_open(&dsfmt);
					long double dz = r*cos(theta);
					long double dy = r*sin(theta);
					long double dx = Length*(dsfmt_genrand_close_open(&dsfmt)-0.5);
				//	Tlocal.CalcCom();
					nx = boundary.center[0] + dx;
					ny = boundary.center[1] + dy;
					nz = boundary.center[2] + dz;
//					//check boundary
//					bool bdist = boundary.;
//			
//					if(!bdist){
//						//boundary failed
//						dE += 1000000.0;
//			
//					}
//					else{
					//boundary passed
					//Check change in energy for pairwise interactions
						dE+=DeltaEInsert(Tlocal, nx, ny, nz);
					//}
					//cout<<"Insertion dE: "<<dE<<endl;
				//}
			}
			else{
				//particle deletion
				//rp = dsfmt_genrand_close_open(&dsfmt);
				long double vp = ((long double)(nlocal)* pow(tw, 3))/V;
				//if (rp<vp){
					idp=vp;
					vdeftag=2;
//					//check boundary
//					bool bdist = Tlocal.CheckComDistAllDelete(rdist, Arandom);
//			
//					if(!bdist){
//						//boundary failed
//						dE += 1000000.0;
//			
//					}
//					else{
						//boundary passed
						//Check change in energy for pairwise interactions
						dE+=DeltaEDelete(Tlocal, Arandom);
					//}
					//cout<<"Deletion dE: "<<dE<<endl;
					//cout<<"particle deletion"<<endl;
				//}
			}
		}
		else{
		//Coordinate change		

			tflag=1;
			//cout<<"translation move"<<endl;
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
			
			//check boundary
			bool bdist = boundary.CheckIfInsideYZAndWrapX(Tlocal.atom[Arandom]);
			
			if(!bdist){
				//boundary failed
				dE += 1000000.0;
			
			}
			else{
				//boundary passed
				//Check change in energy for pairwise interactions
				dE+=DeltaE(Tlocal, Arandom, origx, origy, origz);
			}
			++ctriesloc;
		}
		
		//Selection Criteria--Metropolis Criteria
		bool tagg = 1;
		
		
		//Add difference to current total potential energy
		Et1local+=dE;
		//Calculate the change in Enthalpy from trial moves
		//dH += dE + Press*dV;
		//cout <<" dE: " << dE << endl;
		//Ht1 += dH;
		long double Hp1 = ConfProb(Ht1local, Beta);
		long double Ncurr = (long double)(nlocal);
		Ht1local = Et1local - ChemPo*Ncurr;
		//cout<<"Ht1local: "<<Ht1local<<endl;
		long double Hrp = dsfmt_genrand_close_open(&dsfmt);
		long double Hp2 = ConfProb(Ht1local, Beta);
		long double Hrat = (Hp2/Hp1)*idp;
		if(Hrp < Hrat ){
			//cout<<"trial move accepted! vdeftag is: "<<vdeftag<<endl;
			tagg=0;
			Eelocal=Et1local;
			Helocal = Ht1local;
			
		
			if (vdeftag!=0){
				if (vdeftag==1){
					//update the config for particle insertion
					Config temp(Tlocal);
					
					Tlocal.DeleteAtoms();
					++nlocal;
					Tlocal.InitializeAtoms(nlocal);
					Tlocal.CopyCoordsUnequalAtoms(temp);
					Tlocal.atom[nlocal-1].SetCoord(nx, ny, nz);
									
				}
				else{
					// update the config for particle deletion
					Config temp(Tlocal);
					
					Tlocal.DeleteAtoms();
					--nlocal;
					Tlocal.InitializeAtoms(nlocal);
					temp.SwapAtomCoords(temp.atom[Arandom], temp.atom[nlocal]);
					Tlocal.CopyCoordsUnequalAtoms(temp);
				}
			}
			Nelocal=(long double)(nlocal);
			if (nlocal!=Tlocal.natoms){
				cout<<"numbers of atoms do not match"<<endl;
				exit(0);
			}
					
			
			if (vdeftag>0){
				++vsuccessloc;
			}
			else if(tflag) {
				++csuccessloc;			
			}	
			
		}


				/*Failed criteria--r(dsfmt_genrand_close_open(&dsfmt)-0.5)*vstep;eset atom coordinates and remove
		  the difference in energy */ 
		if(tagg == 1){
			//cout<<"trial move failed. vdeftag: "<<vdeftag<<" nlocal" <<nlocal<<" T.n "<<Tlocal.natoms<<" Aran: "<<Arandom<<" tflag: "<<tflag<<endl;
			if(vdeftag==0 && tflag){
			//	cout<<"T.aAr.x "<<Tlocal.atom[Arandom].x<<" origx "<<origx<<endl;
				Tlocal.atom[Arandom].x = origx;
				Tlocal.atom[Arandom].y = origy;
				Tlocal.atom[Arandom].z = origz;
			//	cout<<"T.aAr.x is now: "<<Tlocal.atom[Arandom].x<<endl;
			}
			Et1local = Et1b;
			Ht1local = Ht1b;
			
		}

			if(n>Eqit && n%(Cit)==0){
			//	cout << "!!!!!Going to push values!!!!.................--------------------------!!!!!" << endl;
			//	cout << "Ee: " << Eelocal << " He: " << Helocal << " Ve: " << Velocal << endl;	
			//	cout << "He*He: " << He*He << endl;	
			//	cout<<"pusing value Helocal: "<<Helocal<<endl;
				#pragma omp flush
				#pragma omp critical
				{
					Estat.Push(Eelocal);
					Esqstat.Push(Eelocal*Eelocal);
					//Hhist.Push(Helocal);
					//long double rog = Tlocal.RadGyr();
					long double rho = Nelocal/V;
					Rgstat.Push(rho);
					Hstat.Push(Helocal);
					Nstat.Push(Nelocal);
				}
			
		}
		
//		if (n>0 && n%10000){
//			Tlocal.CalcCom();
//			Tlocal.ReCenterZero();
//		}

		//Output coordinates of sample n
	   if(n>Eqit && n % (ConfigsOutFreq) == 0) {
		//cout<<"outputting coordinates"<<endl;
		#pragma omp critical
		{
		coutfile << Tlocal.natoms << endl;
		coutfile <<"lj H: "<<Ht1local<<" E: "<<Et1local<<endl;
		for(unsigned long hg = 0; hg<Tlocal.natoms; ++hg){
			coutfile << setiosflags(ios::fixed) << setprecision(15) << "lj " << Tlocal.atom[hg].x << " " << Tlocal.atom[hg].y << " " << Tlocal.atom[hg].z << endl; 
		}
		}
	   }	

	

	//	cout << endl;
	
		
	


		
    } //End Monte Carlo loop
	//cout<<"end MC loop"<<endl;
	#pragma omp barrier
	
	#pragma omp critical
	{ //	cout<<"update global success value"<<endl;
		vsuccess+=vsuccessloc;
		csuccess+=csuccessloc;
		vtries+=vtriesloc;
		ctries+=ctriesloc;
	}

	

}//End parallel
//cout<<"out of parallel region"<<endl;
	long double vafrac = (long double)(vsuccess)/(long double)(vtries);
	long double cafrac = (long double)(csuccess)/(long double)(ctries);
	

	cout<<"vfrac: "<<vafrac<<" cfrac: "<<cafrac<<" step: "<<step<<endl;
//	cout<<"vstep: "<<vstep<<endl;
//cout<<" cfrac: "<<cafrac<<" cfailcount: "<<cfailcount<<endl;

// cout<<"bottom output"<<endl;

long double Vavg = Nstat.Mean();
long double Eavg = Estat.Mean();
long double Havg = Hstat.Mean();
long double Rgavg = Rgstat.Mean();
cout<<" and Average Enthalpy is "<<Havg<<" and Average Particles is: "<<Vavg<<" and Average Energy: "<<Eavg<<endl;
logfile<<" and Average Enthalpy is "<<Havg<<" and Average Particles is: "<<Vavg<<" and Average Energy: "<<Eavg<<endl;
//obsoutfile<<Havg<<" "<<Eavg<<" "<<Vavg<<" "<<Rgavg<<endl;
obsoutfile<<Temp<<" "<<Havg<<" "<<Eavg<<" "<<Vavg<<" "<<Rgavg<<endl;




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
	long double bx = c.boxx;
		for(unsigned int pp = 0;pp<natoms-1;++pp){

			for(unsigned int uu = pp+1;uu<natoms;++uu){
				
				
				Etot += PairEnergy(c.atom[pp], c.atom[uu], bx);
					
			}
		}
//	Etot *= 60.0;
	return(Etot);
}

// potential energy function
long double PairEnergy(const Atom& one, const Atom& two, const long double& boxx){
	
	
	long double dx, dy, dz, rs;
	long double potent, eps, sigma;
	eps = 1.0, sigma = 1.0;
	long double bx=boxx;
//	long double by=boxy;
//	long double bz=boxz;
//	eps = sqrt(one.eps*two.eps), sigma = (one.sig+two.sig)/2.0;
	long double rc = 3.0*sigma;
	long double rcs = rc*rc;
	dx = one.x - two.x;
	
	dy = one.y - two.y;
//	
	dz = one.z - two.z;
	//Minimum image -- coordinates must be pre-wrapped
	if(abs(dx)> bx/2.0){
		dx = bx - abs(one.x) - abs(two.x);
	}
//	if(abs(dy)>by/2.0){
//		dy = by - abs(one.y) - abs(two.y);
//	}
//	if(abs(dz)>bz/2.0){
//		dz = bz - abs(one.z) - abs(two.z);
//	}
	
	rs= (pow(dx, 2))+(pow(dy, 2))+(pow(dz, 2));
	if (rs<rcs){
		potent = 4.0*eps*((pow(sigma/rs, 6))-(pow(sigma/rs, 3))) - 4.0*eps*((pow(sigma/rcs, 6))-(pow(sigma/rcs, 3)));
	}
	else{
		potent=0.0;
	}	
	
	return(potent);
}

//Uses the o position values (xo, yo, zo) from atom one
long double PairEnergyO(const Atom& one, const long double& ox, const long double& oy, const long double& oz, const long double& boxx){
	
	
	long double dx, dy, dz, rs;
	long double potent, eps, sigma;
	long double bx=boxx;

	eps = 1.0, sigma = 1.0;
	long double rc = 3.0*sigma;
	long double rcs = rc*rc;
	dx = one.x - ox;
	
	dy = one.y - oy;
	
	dz = one.z - oz;

	//Minimum image -- coordinates must be pre-wrapped 
	if(abs(dx)> bx/2.0){
		dx = bx - abs(one.x) - abs(ox);
	}
	
//	if(abs(dy)>by/2.0){
//		dy = by - abs(one.y) - abs(oy);
//	}
//	
//	if(abs(dz)>bz/2.0){
//		dz = bz - abs(one.z) - abs(oz);
//	}

	rs= (pow(dx, 2))+(pow(dy, 2))+(pow(dz, 2));
	if (rs<rcs){
		potent = 4.0*eps*((pow(sigma/rs, 6))-(pow(sigma/rs, 3))) - 4.0*eps*((pow(sigma/rcs, 6))-(pow(sigma/rcs, 3)));
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
	long double bx = c.boxx;
		
	
	for (unsigned int i = 0; i < c.natoms; ++i)
	{
		if(i!=Aindex){
			Etot += PairEnergy(c.atom[Aindex], c.atom[i], bx);
			Etoto += PairEnergyO(c.atom[i], ox, oy, oz, bx);
					
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

// calculates the change in energy associated inserting one atom
long double DeltaEInsert(const Config& c, const long double& ix, const long double& iy, const long double& iz){
		
	long double Etoto=0.0;
	long double dE;
	long double bx = c.boxx;
		
	
	for (unsigned int i = 0; i < c.natoms; ++i)
	{
		
			
			Etoto += PairEnergyO(c.atom[i], ix, iy, iz, bx);
					
			
	}
	dE = Etoto;
	//cout << "In DeltaE funct, dE is : " << dE << endl;
	return( dE ); 
	

}

// calculates the change in energy associated with one deleted atom
long double DeltaEDelete(const Config& c, const unsigned& Dindex){
	long double Etot=0.0;	

	long double dE;
	long double bx = c.boxx;
		
	
	for (unsigned int i = 0; i < c.natoms; ++i)
	{
		if(i!=Dindex){
			Etot += PairEnergy(c.atom[Dindex], c.atom[i], bx);
			
					
		}	
	}
	dE = -Etot;
	//cout << "In DeltaE funct, dE is : " << dE << endl;
	return( dE ); 
	

}

long double ConfProb(long double E, long double Beta){
//	cout << "Beta " << Beta << endl;
	long double P = exp(-Beta*E);
	return(P);
}
