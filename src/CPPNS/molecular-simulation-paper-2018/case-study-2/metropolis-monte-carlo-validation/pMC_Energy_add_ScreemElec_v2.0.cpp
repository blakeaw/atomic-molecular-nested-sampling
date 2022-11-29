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
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Additive/code/src/class/Atom.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Additive/code/src/class/Config.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/NSTAM/code/src/class/RunningStats.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/NSTAM/code/src/class/HistogramNS_slim.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/NSTAM/code/src/class/Histogram.h"
#include "/net/uu/nm/cm/bxw109120/src/SimulationCpp/class/HardBox.h"
#include "/net/uu/nm/cm/bxw109120/src/SimulationCpp/class/BlockAverage.h"
//#include "./src/class/CenterofMass.h"


// prototype functions

long double PairEnergy(const Atom& one, const Atom& two, long double bx, long double by, long double bz);
long double PairEnergyO(const Atom& one, const Atom& two, const long double& ox, const long double& oy, const long double& oz, long double bx, long double by, long double bz);
long double TotalEnergy(const Config& c);
long double TotalEnergyTwo(const Config& c1, const Config& c2);
long double DeltaE(const Config& c, const unsigned& Aindex, const long double& ox, const long double& oy, const long double& oz);
string itos( int Number );
long double RDist(const Atom& one, const Atom& two);
long double RDistMinImage(const Atom& one, const Atom& two, long double bx, long double by, long double bz);
long double LocalDensity(const Config& c);
long double LocalDensityTwo(const Config& ca, const Config& cb);
long double TestArea(const Config& c, long double scale, long double Ec);
long double ConfProb(long double E, long double Beta);
long double PairEnergyES(const Atom& one, const Atom& two, long double bx, long double by, long double bz);
long double TotalEnergyES(const Config& c);
long double TotalEnergyTwoES(const Config& c1, const Config& c2);
long double PairEnergyESO(const Atom& one, const Atom& two, const long double& ox, const long double& oy, const long double& oz, long double bx, long double by, long double bz);
long double DeltaEES(const Config& c, const unsigned& Aindex, const long double& ox, const long double& oy, const long double& oz);
long double DeltaEEStwo(const Config& c, const unsigned& Aindex, const long double& ox, const long double& oy, const long double& oz, const Config& s);
long double DeltaEtwo(const Config& c, const unsigned& Aindex, const long double& ox, const long double& oy, const long double& oz, const Config& s);

int main(int argc, char *argv[])
{

	//Get the local time for the logfile
	time_t rawtime;
 	struct tm * timeinfo;
 	time (&rawtime);
 	timeinfo = localtime (&rawtime);

	// Parameters that need to be set

	// the number of atoms to have in the system
	
	const unsigned int na = 100;
	const unsigned int nb = 100;
	const unsigned int natoms = na+nb;
	//starting c value
	const long double Cstart = 0.0;
	const long double Cend = 2.0;
	//temperature
	const long double Temp = 0.7;
	// the number of nested sampling iterations to run 
	const unsigned int Cloop = 30;
	//Total number of MC style trial moves
	const unsigned int MCit = 20000000;
	//Number of procs to split Sampling loop (MCit)
	const unsigned short nprocs = 1;
	// Number of initial MC steps to allow equilibration
	const unsigned Eqit = 80000000;
	// Step interval to collect data after equilibration
	const unsigned Cit = natoms*10;
	//long double NSfrac = 0.5;
	//Box sizes
	long double boxx = 6.0;//*1.1;
	long double boxy = 6.0;//*1.1;
	long double boxz = 7.0;///(1.1*1.1);
	// Test area scaling value
	long double scale = 1.00001;
	long double sca = 1.01;
//	boxx*=sca;
//	boxy*=sca;
//	boxz/=(sca*sca);
	long double abox = (boxx+boxy+boxz)/3.0;
	//Translation step
	const long double bsfac = 0.900;
 	
    long double stepx = boxx*bsfac;
	long double stepy = boxy*bsfac;
	long double stepz = boxz*bsfac;
	

	//Pi
	//const long double pi = 3.14159265359;
	
	//Enthalpy bounds
	const long double Etop = 8000.0;
	const long double Elow = -2000.0;
	long double Binwidth = 1.0e-3;

		//run descriptor
	string nds = itos(natoms);
	string rundescript = "MC_n"+nds+"_testb";
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
	long double dC = (Cend - Cstart)/ (long double)(Cloop-1);
	//Name of logfile
	string logfilename = "oMMC_" + rundescript + ".log";
//	//Name of Emedian output file
//	string houtname = "Epoints_" + rundescript + ".dat";
	//Name of Configuration output file
	string coutname = "Coord_" + rundescript + ".xyz";
	//Name of Observable output file
	string outfilename = "zObservables_" + rundescript + ".dat";
	string ehistoutname = "zEabhist_"+rundescript;
	
	ofstream logfile(logfilename.c_str());
	//ofstream houtfile(houtname.c_str());
	ofstream obsoutfile(outfilename.c_str());
	//string coutfilename = "Coord_NS_Ent_a0.xyz";
	ofstream coutfile(coutname.c_str());
	ofstream eoutfile(ehistoutname.c_str());
//	ofstream stepoutfile("zStepOut_int_test_a1.dat");	
	// configs from each nested sampling iteration
	 
	unsigned long long MCitpp = (MCit+Eqit*nprocs)/nprocs;
	unsigned int ConfigsOutFreq = ((MCitpp-Eqit))/30;

//cout << "Nested sampling iterations is set to " << Iterations << " with MC iterations set to " << MCit << endl;
//cout<<"Etop: "<<Etop<<" Elow: "<<Elow<<endl;	

logfile << " " << endl;
	logfile << "Local time and date: " << asctime(timeinfo) << endl;
	logfile<<"NVT Monte Carlo of lennard jones system with variable ab epsilon"<<endl;
	logfile << "Parameters: " << endl;
	logfile<<"Parameter Iterations: "<<Cloop<<" Number of Trial Moves (MCit): "<<MCit<<endl;
//	logfile<<"NS fraction : "<<NSfrac<<endl;
//	logfile<<"TAM scaling factor: "<<scale<<endl;
	logfile<<" na: "<<na<<" nb: "<<nb<<endl;
	logfile<<"natoms: "<<natoms<<" nprocs: "<<nprocs<<" Eqit: "<<Eqit<<" Cit: "<<Cit<<endl;
	logfile<<" boxx: "<<boxx<<" boxy: "<<boxy<<" boxz: "<<boxz<<endl;
	logfile<<"bsfac: "<<bsfac<<" initial guess stepx: "<<stepx<<" stepy: "<<stepy<<" "<<stepz<<endl;
	//logfile<<" Etop: "<<Etop<<" Elow: "<<Elow<<endl;	
	logfile<<"MCitpp: "<<MCitpp<<"ConfigsOutFreq: "<<ConfigsOutFreq<<endl;
	logfile<<"Temperature: "<<Temp<<endl;

	cout<<" na: "<<na<<" nb: "<<nb<<endl;
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
//	HistogramNS Ehist;
//	Ehist.Initialize(Binwidth, Elow, Etop);
//	Ehist.SetNSfrac(NSfrac);
	BlockAverage Estat;
	Estat.Initialize(20, 5000);
	RunningStats Esqstat;
	RunningStats LDstat;
	RunningStats dUstat;
	RunningStats dUsstat;
	RunningStats dUcstat;
	BlockAverage AAdstat;
	BlockAverage ABdstat;
	AAdstat.Initialize(20, 5000);
	ABdstat.Initialize(20, 5000);
//	RunningStats Rgstat;
//	Histogram Eohist;
//	long double ldlow = Elow;
//	long double ldhigh = Etop;
//	//Eohist.Initialize(Binwidth*1000.0, ldlow, ldhigh);
	
	//logfile<<"number of NS bins: "<<Ehist.GetNumberOfBins()<<endl;
//	logfile<<"Potential Energy Histogram Params: "<<endl;
//	logfile<<"number of bins: "<<Eohist.GetNumberOfBins()<<" bindwidth: "<<Eohist.GetIncrement()<<" Elow: "<<Elow<<" Etop: "<<Etop<<endl;
//Initialize configuration objects

long double Beta = 1.0/Temp;
	
Config trial1a(na);
Config trial1b(nb);
Config Bottoma(na);
Config Bottomb(nb);
Bottoma.SetBox(boxx, boxy, boxz);
Bottomb.SetBox(boxx, boxy, boxz);
trial1a.SetBox(boxx, boxy, boxz);
trial1b.SetBox(boxx, boxy, boxz);
//Config Lconfa(na);
//Lconfa.SetBox(boxx, boxy, boxz);
//Config Lconfb(nb);
//Lconfb.SetBox(boxx, boxy, boxz);
for (unsigned int i = 0; i < na; ++i)
{
	trial1a.atom[i].ctag=0;
}
for (unsigned int i = 0; i < nb; ++i)
{
	trial1b.atom[i].ctag=1;
}
cout<<"Assigning Initial coordinates.."<<endl;
logfile<<"Assigning Initial coordinates.."<<endl;
//Randomly place atoms in box - but filter to configuration with Enthalpy below Htop

long double Econf;

do
{

	//set a
	unsigned track=0;
	for (unsigned int i = 0; i < na; ++i)
	{
		long double rx = (dsfmt_genrand_close_open(&dsfmt)-0.5)*boxx;
		long double ry = (dsfmt_genrand_close_open(&dsfmt)-0.5)*boxy;
		long double rz = (dsfmt_genrand_close_open(&dsfmt))*(-boxz/2.0);
		trial1a.atom[i].x = rx;
		trial1a.atom[i].y = ry;
		trial1a.atom[i].z = rz;
		bool flag=0;
		//cout<<"assigned atom: "<<i<<endl;
		if (track<2*na){
			//cout<<"checking distance..."<<endl;
			for (unsigned int j = 0; j < i; ++j)
			{
				long double rdist = RDistMinImage(trial1a.atom[i], trial1a.atom[j], boxx, boxy, boxz);
				//cout<<"distance from atom "<<i<<" to atom "<<j<<" is "<<rdist<<endl;
				if (rdist<0.8){
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
	//set b
	track=0;
	for (unsigned int i = 0; i < nb; ++i)
	{
		long double rx = (dsfmt_genrand_close_open(&dsfmt)-0.5)*boxx;
		long double ry = (dsfmt_genrand_close_open(&dsfmt)-0.5)*boxy;
		long double rz = (dsfmt_genrand_close_open(&dsfmt))*boxz/2.0;
		trial1b.atom[i].x = rx;
		trial1b.atom[i].y = ry;
		trial1b.atom[i].z = rz;
		bool flag=0;
		//cout<<"assigned atom: "<<i<<endl;
		if (track<2*nb){
			//cout<<"checking distance..."<<endl;
			for (unsigned int j = 0; j < i; ++j)
			{
				long double rdist = RDistMinImage(trial1b.atom[i], trial1b.atom[j], boxx, boxy, boxz);
				//cout<<"distance from atom "<<i<<" to atom "<<j<<" is "<<rdist<<endl;
				if (rdist<0.8){
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

	Econf = (TotalEnergyES(trial1a)+TotalEnergyES(trial1b)+TotalEnergyTwoES(trial1a, trial1b))*Cstart;
	//Econf = TotalEnergy(trial1);

//	cout<<Econf<<endl;
	
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
//HardBox boundary;
//boundary.SetLengths(boxx, boxy, boxz);

cout<<"Initial coordinates assigned!"<<endl;
logfile<<"Initial coordinates assigned!"<<endl;
Bottoma.FullEquate(trial1a);
Bottomb.FullEquate(trial1b);	
long double Egmin = Cstart*(TotalEnergyES(trial1a)+TotalEnergyES(trial1b)+TotalEnergyTwoES(trial1a, trial1b));

long double Et1 = Egmin;

cout<<"Initial local densities; aa: "<<LocalDensity(trial1a)<<" bb: "<<LocalDensity(trial1b)<<" ab "<<LocalDensityTwo(trial1a, trial1b)<<" ba "<<LocalDensityTwo(trial1b, trial1a)<<endl;
// Output starting coordinates
coutfile << natoms << endl;
coutfile <<"lj "<<natoms<<" E: "<<Egmin<<endl;
for(unsigned long hg = 0; hg<trial1a.natoms; ++hg){
		coutfile << setiosflags(ios::fixed) << setprecision(15) << "a " << Bottoma.atom[hg].x << " " << Bottoma.atom[hg].y << " " << Bottoma.atom[hg].z << endl; 
		}
for(unsigned long hg = 0; hg<trial1b.natoms; ++hg){
		coutfile << setiosflags(ios::fixed) << setprecision(15) << "b " << Bottomb.atom[hg].x << " " << Bottomb.atom[hg].y << " " << Bottomb.atom[hg].z << endl; 
		}
long double Enest = Etop;
long double Elowprev = Et1;
long double Elconf = Et1;
unsigned long long cfailcount = 0;
long double Enprev=Enest;

for(unsigned int ii=0;ii<Cloop;++ii){
cout << "Iteration " << ii << "......." << endl;
logfile << "Iteration " << ii << "......." << endl;
long double Ccurr = Cstart + dC*(long double)(ii);
//if (ii>0){
	Et1 = (TotalEnergyES(trial1a)+TotalEnergyES(trial1b)+TotalEnergyTwoES(trial1a, trial1b))*Ccurr;
//}
//houtfile<<Enest<<endl;

cout<<"Initial Energy: "<<Et1<<" Ccurr "<<Ccurr<<endl;
logfile<<"Initial Energy: "<<Et1<<" Ccurr "<<Ccurr<<endl;

	//unsigned long long acceptcount = 0;
	
	unsigned long long ctries = 0;
	unsigned long long csuccess = 0;
	//long double Elowprev = Et1;
	bool coflag = 1;

//	omp_set_num_threads(nprocs);
	#pragma omp parallel 
	{

	 //Define a thread local configuration
	 Config Tlocala(na);
	 Tlocala.FullEquate(trial1a);
	 Tlocala.SetBox(boxx, boxy, boxz);
	 Config Tlocalb(nb);
	 Tlocalb.FullEquate(trial1b);
	 Tlocalb.SetBox(boxx, boxy, boxz);
	// Config Tlocold(natoms);

//Define thread local variables
	long double Et1local = Et1;
	
	long double Eelocal = 0.0;
	long double Eljcrossloc = TotalEnergyTwo(Tlocala, Tlocalb);	 
	cout<<"Eljcrossloc start "<<Eljcrossloc<<endl;
	unsigned long long ctriesloc = 0;
	unsigned long long csuccessloc = 0;
	long double Elja = TotalEnergy(Tlocala);
		long double Eljb = TotalEnergy(Tlocalb);
  //Perform markov moves on the configuration
    for(unsigned int n = 0; n<MCitpp; ++n){
		
//	 Elja = TotalEnergy(Tlocala);
//		 Eljb = TotalEnergy(Tlocalb);
//		Eljcrossloc = TotalEnergyTwo(Tlocala, Tlocalb);	 
		//Store Initial Et1 value
		long double Et1b = Et1local;
	//	long double Eljcb = Eljcrossloc;
		long double Eljab=Elja;
		long double Eljbb=Eljb;
		//Trial Move -----
		// randomly select an atom index	
			unsigned int Arandom = (unsigned int)(dsfmt_genrand_close_open(&dsfmt)*natoms);
			//phase flag: 0 for a and 1 for b
			bool pflag;
			if (Arandom/na == 1){
				pflag=1;
				Arandom -= na;
			}
			else{
				pflag=0;
			}
			long double origx, origy, origz;
			if (pflag==0){
				origx = Tlocala.atom[Arandom].x;
				origy = Tlocala.atom[Arandom].y;
				origz = Tlocala.atom[Arandom].z;
			}
			else if(pflag==1){
					origx = Tlocalb.atom[Arandom].x;
				origy = Tlocalb.atom[Arandom].y;
				origz = Tlocalb.atom[Arandom].z;
			}
		
		// Calculate the difference in total potential energy
		//long double dH = 0.0;
		long double dE = 0.0;
		//long double dV = 0.0;
		


		//Coordinate change		

			
			long double movex, movey, movez;
			/* set random movements of the coordinates
			   of the selected atom */
			movex = (dsfmt_genrand_close_open(&dsfmt)-0.5)*stepx;
			movey = (dsfmt_genrand_close_open(&dsfmt)-0.5)*stepy;
			movez = (dsfmt_genrand_close_open(&dsfmt)-0.5)*stepz;
			bool bflag = 1;
			// Implement those moves
			if (pflag==0){
				Tlocala.atom[Arandom].x += movex;
				Tlocala.atom[Arandom].y += movey;
				Tlocala.atom[Arandom].z += movez;
				//bflag = boundary.CheckBoundary(Tlocala.atom[Arandom]);
				Tlocala.WrapAtom(Arandom);
			}
			else if(pflag==1){
				Tlocalb.atom[Arandom].x += movex;
				Tlocalb.atom[Arandom].y += movey;
				Tlocalb.atom[Arandom].z += movez;
				//bflag = boundary.CheckBoundary(Tlocalb.atom[Arandom]);
				Tlocalb.WrapAtom(Arandom);
			}
			++ctriesloc;
//			Tlocal.CalcCom();
//			Tlocal.ReCenterZero();
//			Tlocal.WrapCoordinates();
			//Check change in energy for pairwise interactions
			bool Tflag = 0;
			if (bflag){
				
				long double dEpljc=0.0;
				long double dEpes=0.0;
				if (pflag==0){
					dE = DeltaE(Tlocala, Arandom, origx, origy, origz);
					dEpljc= DeltaEtwo(Tlocala, Arandom, origx, origy, origz, Tlocalb);
					dEpes=(DeltaEES(Tlocala, Arandom, origx, origy, origz))*Ccurr;
					dEpes+=(DeltaEEStwo(Tlocala, Arandom, origx, origy, origz, Tlocalb))*Ccurr;
				}
				else if (pflag==1){
					dE = DeltaE(Tlocalb, Arandom, origx, origy, origz);
					dEpljc= DeltaEtwo(Tlocalb, Arandom, origx, origy, origz, Tlocala);
					dEpes=(DeltaEES(Tlocalb, Arandom, origx, origy, origz))*Ccurr;
					dEpes+=(DeltaEEStwo(Tlocalb, Arandom, origx, origy, origz, Tlocala))*Ccurr;
				}
			
				//long double dEp=DeltaE(Tlocal, Arandom, origx, origy, origz);
			//	long double Enew = TotalEnergy(Tlocal);
				//dE = Enew - Et1b;
				//dE = dEp;
				//long double dEpljc;
//				cout<<"pflag "<<pflag<<endl;
//				Eljcrossloc+=dEpljc;
//				cout<<"a Eljcross "<<Eljcrossloc<<" dEpljc "<<dEpljc<<endl;
//				Eljcrossloc = TotalEnergyTwo(Tlocala, Tlocalb);
//				dEpljc = Eljcrossloc - Eljcb;
//				cout<<"b Eljcross "<<Eljcrossloc<<" dEpljc "<<dEpljc<<endl;
//				
//				if (pflag){
//					Eljb+=dE;
//				}
//				else{
//					Elja+=dE;
//				}
//				cout<<"Eljab "<<Eljab<<" Eljbb "<<Eljbb<<endl;
//				cout<<"dE "<<dE<<" Elja "<<Elja<<" Eljb "<<Eljb<<endl;
//				Elja=TotalEnergy(Tlocala);
//				Eljb=TotalEnergy(Tlocalb);
//				long double dEa = Elja-Eljab;
//				long double dEb = Eljb-Eljbb;
//				cout<<"dEa "<<dEa<<" dEb "<<dEb<<" Elja "<<Elja<<" Eljb "<<Eljb<<endl;
//				if (n==10){
//							exit(0);
//				}
		
				//Et1local= (TotalEnergyES(Tlocala)+TotalEnergyES(Tlocalb)+TotalEnergyTwoES(Tlocala, Tlocalb))*Ccurr;			
				//dEpes*=Ccurr;
				Et1local+=dEpes;
				//long double dEpelc = Et1local-Et1b;
				dE+=dEpljc+dEpes;
				//Thermal criterion
				long double Prand = dsfmt_genrand_close_open(&dsfmt);
				long double Prat = ConfProb(dE, Beta);
				if (Prat>Prand){
					Tflag=1;
				}


				
			}
			
		
		//Selection Criteria--Metropolis Criteria
		bool tagg = 1;
		
		
		//Add difference to current total potential energy
			
			//cout<<"pflag "<<pflag<<" dE "<<dE<<" Ar "<<Arandom<<" Prand "<<Prand<<" Prat "<<Prat<<" Et1b "<<Et1b<<" Et1local "<<Et1local<<" Enest "<<Enest<<endl;
			
		if(bflag && Tflag){
			
			
			tagg=0;
			Eelocal=Et1local;
//			Tlocal.CalcCom();
//			Tlocal.ReCenterZero();
			if (Eelocal<Elowprev){

				#pragma omp critical
				{
					
					Elowprev=Eelocal;
					Bottoma.FullEquate(Tlocala);
					Bottomb.FullEquate(Tlocalb);
					//Bottom.SetBox(boxx, boxy, blocal);
				}
			}
//			else if (coflag){
//				if (Eelocal>Elowprev && Eelocal<Enest-1.0){
//					Elconf = Eelocal;
//					Lconfa.Equate(Tlocala);
//					Lconfb.Equate(Tlocalb);
//					coflag=0;
//				}
//			}
			
			
				++csuccessloc;			
			
		}


				/*Failed criteria--r(dsfmt_genrand_close_open(&dsfmt)-0.5)*vstep;eset atom coordinates and remove
		  the difference in energy */ 
		if(tagg == 1){
			
				if (pflag==0){
					Tlocala.atom[Arandom].x = origx;
				Tlocala.atom[Arandom].y = origy;
				Tlocala.atom[Arandom].z = origz;
				}
				else if(pflag==1){
					Tlocalb.atom[Arandom].x = origx;
				Tlocalb.atom[Arandom].y = origy;
				Tlocalb.atom[Arandom].z = origz;

				}
				
			
			Et1local = Et1b;
			//Eljcrossloc=Eljcros
			
		}

			if(n>Eqit && (n%(Cit))==0){
			//	cout << "!!!!!Going to push values!!!!.................--------------------------!!!!!" << endl;
			//	cout << "Ee: " << Eelocal << " He: " << Helocal << " Ve: " << Velocal << endl;	
			//	cout << "He*He: " << He*He << endl;	
				long double LD = LocalDensity(Tlocala);
				long double dU = LocalDensityTwo(Tlocala, Tlocalb);
//				long double dU = TestArea(Tlocal, scale, Eelocal);
				//long double LD =0.0;
				//long double dU = 0.0;
				#pragma omp critical
				{
				//	Ehist.Push(Eelocal);
					Estat.Push(Eelocal);
					Esqstat.Push(Eelocal*Eelocal);
					LDstat.Push(LD);
					dUstat.Push(dU);
					dUsstat.Push(dU*dU);
					dUcstat.Push(dU*dU*dU);
					AAdstat.Push(LD);
					ABdstat.Push(dU);
					//Eohist.Push(Eelocal);
				}
			
		}

		if (n>0 && n%2000==0){
							//check acceptance ratios and adjust
						long double art = (long double)(csuccessloc)/(long double)(ctriesloc);
						//cout<<"n "<<n<<" art "<<art<<endl;
						if (art<0.10){
							stepx*=0.5;
							stepy*=0.5;
							stepz*=0.5;
						}
						else if (art<0.33){
							long double msc = art/0.33;
							stepx*=msc;
							stepy*=msc;
							stepz*=msc;
						
						}
						else if(art>0.51){
							long double msc = art/0.51;
							stepx*=msc;	
							stepy*=msc;
							stepz*=msc;
					
						}
						csuccessloc=0;
						ctriesloc=0;
						
					}
		//Output coordinates of sample n
	   if(n>Eqit && (n-Eqit)%(ConfigsOutFreq)==0) {
			#pragma omp critical 
			{
				//cout<<"Output coordinates to file at trial move: "<<n<<endl;
				coutfile << natoms << endl;
				coutfile <<"lj "<<natoms<<" E: "<<Et1local<< " eab: "<<Ccurr<<endl;
				for(unsigned long hg = 0; hg<Tlocala.natoms; ++hg){
					coutfile << setiosflags(ios::fixed) << setprecision(15) << "a " << Tlocala.atom[hg].x << " " << Tlocala.atom[hg].y << " " << Tlocala.atom[hg].z << endl; 
				}
				for(unsigned long hg = 0; hg<Tlocalb.natoms; ++hg){
					coutfile << setiosflags(ios::fixed) << setprecision(15) << "b " << Tlocalb.atom[hg].x << " " << Tlocalb.atom[hg].y << " " << Tlocalb.atom[hg].z << endl; 
				}
			}
		//cout<<" -----------------------------"<<endl;

	   }	

	

	//	cout << endl;
	
		
	


		
    } //End Monte Carlo loop

	#pragma omp barrier
	Eljcrossloc=TotalEnergyTwo(Tlocala, Tlocalb);
	cout<<"Eljcross end "<<Eljcrossloc<<endl;
//	#pragma omp critical
//	{ 
//		
//		csuccess+=csuccessloc;
//		ctries+=ctriesloc;

//	}

}//End parallel

	
//	long double cafrac = (long double)(csuccess)/(long double)(ctries);
//	
//	if (cafrac<0.30){
//		++cfailcount;
//	}
//	cout<<" cfrac: "<<cafrac<<" cfailcount: "<<cfailcount<<endl;
//	logfile<<" cfrac: "<<cafrac<<" cfailcount: "<<cfailcount<<endl;
//	logfile<<" stepx: "<<stepx<<" stepy: "<<stepy<<" stepz: "<<stepz<<endl;
cout<<" stepx: "<<stepx<<" stepy: "<<stepy<<" stepz: "<<stepz<<endl;
//		coutfile << trial1.natoms << endl;
//		coutfile <<"lj 17 H: "<<Htlowprev<<" E: "<<Elowprev<<" V: "<<Vlowprev<<endl;
//		for(unsigned long hg = 0; hg<trial1.natoms; ++hg){
//			coutfile << setiosflags(ios::fixed) << setprecision(15) << "lj " << Bottom.atom[hg].x << " " << Bottom.atom[hg].y << " " << Bottom.atom[hg].z << endl; 
//		}
//	if (cafrac<0.10){
//		cout<<"config frozen. exiting!"<<endl;
//		exit(0);
//	}
 
//Hhist.Normalize();
//Hnest=Hhist.Median();
long double Eavg = Estat.Average();
long double Esqavg = Esqstat.Mean();
long double LDavg = LDstat.Mean();
long double dUavg = dUstat.Mean();
long double dUsavg = dUsstat.Mean();
long double dUcavg = dUcstat.Mean();
//long double Rgavg = Rgstat.Mean();
cout<<" Average Energy is "<<Eavg<<" and Average Energy Squared: "<<Esqavg<<endl;
logfile<<" Average Energy is "<<Eavg<<" and Average Energy Squared: "<<Esqavg<<endl;
cout<<"Average Local Density: "<<LDavg<<endl;
logfile<<"Average Local Density: "<<LDavg<<endl;
cout<<"Average dU: "<<dUavg<<endl;
logfile<<"Average dU: "<<dUavg<<endl;
obsoutfile<<Ccurr<<" "<<Eavg<<" "<<Estat.Uncertainty()<<" "<<LDavg<<" "<<dUavg<<" "<<dUsavg<<" "<<dUcavg<<" "<<AAdstat.Average()<<" "<<AAdstat.Uncertainty()<<" "<<ABdstat.Average()<<" "<<ABdstat.Uncertainty()<<endl;
//string eoname = ehistoutname+"_it_"+itos(ii-1)+".dat";
//Eohist.Normalize();
//Eohist.DumpNormal(eoname);
//Eohist.Reset();


Estat.Reset();
Esqstat.Reset();
LDstat.Reset();
dUstat.Reset();
dUsstat.Reset();
dUcstat.Reset();
AAdstat.Reset();
ABdstat.Reset();
//Rgstat.Reset();
//if (Elconf<Enest){
//	trial1a.Equate(Lconfa);
//	trial1b.Equate(Lconfb);
//	Et1 = Elconf;
//	cout<<"setting new start conf to Lconf"<<endl;
//}
//else{
	trial1a.FullEquate(Bottoma);
	trial1b.FullEquate(Bottomb);
	//Et1=Elowprev;
	cout<<"setting new start conf to Bottom"<<endl;

//trial1.SetBox(boxx, boxy, boxz);
//abox = (boxx+boxy+boxz)/3.0;

//stepx = bsfac*boxx*pow(0.90, cfailcount);
//stepy = bsfac*boxy*pow(0.90, cfailcount);
//stepz = bsfac*boxz*pow(0.90, cfailcount);





} // End NS loop

logfile.close();
//houtfile.close();
obsoutfile.close();
coutfile.close();
//stepoutfile.close();
logfile<<"Simulation Complete!"<<endl;
logfile<<"Output files are: "<<endl;
//logfile<<"Nested Energies: "<<houtname<<endl;
logfile<<"Configuration Samples: "<<coutname<<endl;
logfile<<"Average Observable Values: "<<outfilename<<endl;
logfile<<"Energy Histograms : "<<ehistoutname<<endl;
cout<<"Simulation Complete!"<<endl;
cout<<"Output files are: "<<endl;
cout<<"Logfile: "<<logfilename<<endl;
//cout<<"Nested Energies: "<<houtname<<endl;
cout<<"Configuration Samples: "<<coutname<<endl;
cout<<"Average Observable Values: "<<outfilename<<endl;
cout<<"Energy Histograms : "<<ehistoutname<<endl;
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
long double TotalEnergyTwo(const Config& c1, const Config& c2){
	
	long double Etot = 0.0;
	unsigned int na = c1.natoms;
	unsigned int nb = c2.natoms;
//	cout<<"nb "<<nb<<" na "<<na<<endl;
	long double bx, by, bz;
	bx = c1.boxx, by = c1.boxy, bz = c1.boxz;
		for(unsigned int pp = 0;pp<na;++pp){

			for(unsigned int uu = 0;uu<nb;++uu){
				
				long double PE = PairEnergy(c1.atom[pp], c2.atom[uu], bx, by, bz);
				//cout<<"PE "<<PE<<endl;
				Etot += PE;
					
			}
		}
//	Etot *= 60.0;
	//Etot/=4.0;
	return(Etot);
}
// potential energy function
long double PairEnergy(const Atom& one, const Atom& two, long double bx, long double by, long double bz){
	
	
	long double dx, dy, dz, rs;
	long double potent, eps, sigma;
	if (one.ctag==two.ctag){
		eps = 1.0, sigma = 1.0;
	}
	else{
		eps = 0.1, sigma = 1.0;
	}
	long double rc=3.0;
	long double rcs=rc*rc;
	dx = one.x - two.x;
	
	dy = one.y - two.y;
	
	dz = one.z - two.z;
	//cout<<"bx "<<bx<<" by "<<by<<" bz "<<bz<<endl;
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
	if (rs<rcs){
		potent = 4.0*eps*((pow(sigma/rs, 6))-(pow(sigma/rs, 3)))  -  4.0*eps*((pow(sigma/rcs, 6))-(pow(sigma/rcs, 3)));
	}
	else{
		potent=0.0;
	}

		

	return(potent);
}

//Uses the o position values (xo, yo, zo) from atom one
long double PairEnergyO(const Atom& one, const Atom& two, const long double& ox, const long double& oy, const long double& oz, long double bx, long double by, long double bz){
	
	
	long double dx, dy, dz, rs;
	long double potent, eps, sigma;
	if (one.ctag==two.ctag){
		eps = 1.0, sigma = 1.0;
	}
	else{
		eps = 0.1, sigma = 1.0;
	}
	long double rc=3.0;
	long double rcs=rc*rc;
	//cout<<"bx "<<bx<<" by "<<by<<" bz "<<bz<<endl;
//	long double rc = 3.0*sigma;
//	long double rcs = rc*rc;
	dx = one.x - ox;
	
	dy = one.y - oy;
	
	dz = one.z - oz;
//Minimum image -- coordinates must be pre-wrapped
	if(abs(dx)> bx/2.0){
		dx = bx - abs(one.x) - abs(ox);
	}
	if(abs(dy)>by/2.0){
		dy = by - abs(one.y) - abs(oy);
	}
	if(abs(dz)>bz/2.0){
		dz = bz - abs(one.z) - abs(oz);
	}

	rs= (pow(dx, 2))+(pow(dy, 2))+(pow(dz, 2));

	if (rs<rcs){
		potent = 4.0*eps*((pow(sigma/rs, 6))-(pow(sigma/rs, 3)))  -  4.0*eps*((pow(sigma/rcs, 6))-(pow(sigma/rcs, 3)));
	}
	else{
		potent=0.0;
	}
	//cout<<"rs "<<rs<<" rcs "<<rcs<<" potent "<<potent<<endl;
	return potent;
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
			Etoto += PairEnergyO(c.atom[i], c.atom[Aindex], ox, oy, oz, bx, by, bz);
					
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
		for(unsigned int pp = 0;pp<natoms;++pp){
			long double count = 0.0;
			for(unsigned int uu = 0;uu<natoms;++uu){
				
				if (uu!=pp){

					r= RDistMinImage(c.atom[pp], c.atom[uu], bx, by, bz);
				if (r<cutoff){
					count+=1.0;
				}

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
long double LocalDensityTwo(const Config& ca, const Config& cb){
	
	long double r;
	unsigned int na = ca.natoms;
	unsigned int nb = cb.natoms;
	long double bx, by, bz;
	bx = ca.boxx, by = ca.boxy, bz = ca.boxz;
	long double cutoff = 2.0;
	long double volume = (4.0/3.0)*3.14*pow(cutoff, 3);
	long double avgden = 0.0;
	long double dc = 0.0;
		for(unsigned int pp = 0;pp<na;++pp){
			long double count = 0.0;
			for(unsigned int uu = 0;uu<nb;++uu){
				
				
				r= RDistMinImage(ca.atom[pp], cb.atom[uu], bx, by, bz);
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

long double ConfProb(long double E, long double Beta){
//	cout << "Beta " << Beta << endl;
	long double P = exp(-Beta*E);
	return(P);
}

// potential energy function
long double PairEnergyES(const Atom& one, const Atom& two, long double bx , long double by, long double bz){
	
	
	long double dx, dy, dz, rs;
	long double potent, A;
	long double rc = 3.0;
	long double kappa = 1.0;
	if (one.ctag==two.ctag){
		A=1.0;
	}
	else{
		A=-1.0;
	}
	
	dx = one.x - two.x;
	
	dy = one.y - two.y;
	
	dz = one.z - two.z;
	
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
	long double r = sqrt(rs);
	if (r<rc){
		potent = (A/(r*kappa))*exp(-kappa*r) - (A/(rc*kappa))*exp(-kappa*rc);
	}
	else{
		potent=0.0;
	}	
	
	return(potent);
}

// potential energy function
long double PairEnergyESO(const Atom& one, const Atom& two, const long double& ox, const long double& oy, const long double& oz, long double bx, long double by, long double bz){
	
	
	long double dx, dy, dz, rs;
	long double potent, A;
	long double kappa = 1.0;
	long double rc=3.0;
	if (one.ctag==two.ctag){
		A=1.0;
	}
	else{
		A=-1.0;
	}
	//cout<<"bx "<<bx<<" by "<<by<<" bz "<<bz<<endl;
	dx = one.x - ox;
	
	dy = one.y - oy;
	
	dz = one.z - oz;
	
	//Minimum image -- coordinates must be pre-wrapped
	if(abs(dx)> bx/2.0){
		dx = bx - abs(one.x) - abs(ox);
	}
	if(abs(dy)>by/2.0){
		dy = by - abs(one.y) - abs(oy);
	}
	if(abs(dz)>bz/2.0){
		dz = bz - abs(one.z) - abs(oz);
	}
	
	rs= (pow(dx, 2))+(pow(dy, 2))+(pow(dz, 2));
	long double r = sqrt(rs);
	if (r<rc){
		potent = (A/(r*kappa))*exp(-kappa*r) - (A/(rc*kappa))*exp(-kappa*rc);
	}
	else{
		potent=0.0;
	}	
	
	return(potent);
}

// Total potential energy calculator, uses only the lj pot
long double TotalEnergyES(const Config& c){
	
	long double Etot = 0.0;
	unsigned int natoms = c.natoms;

		for(unsigned int pp = 0;pp<natoms-1;++pp){

			for(unsigned int uu = pp+1;uu<natoms;++uu){
				
				
				Etot += PairEnergyES(c.atom[pp], c.atom[uu], c.boxx, c.boxy, c.boxz);
					
			}
		}
//	Etot *= 60.0;
	return(Etot);
}
// calculates the change in energy associated with one moved atom
long double DeltaEES(const Config& c, const unsigned& Aindex, const long double& ox, const long double& oy, const long double& oz){
	long double Etot=0.0;	
	long double Etoto=0.0;
	long double dE;
	long double bx, by, bz;
	bx = c.boxx, by = c.boxy, bz = c.boxz;
		
	
	for (unsigned int i = 0; i < c.natoms; ++i)
	{
		if(i!=Aindex){
			Etot += PairEnergyES(c.atom[Aindex], c.atom[i], bx, by, bz);
			Etoto += PairEnergyESO(c.atom[i], c.atom[Aindex], ox, oy, oz, bx, by, bz);
					
		}	
	}
	dE = Etot - Etoto;
	//cout << "In DeltaE funct, dE is : " << dE << endl;
	return( dE ); 
	

}
// calculates the change in energy associated with one moved atom
long double DeltaEEStwo(const Config& c, const unsigned& Aindex, const long double& ox, const long double& oy, const long double& oz, const Config& s){
	long double Etot=0.0;	
	long double Etoto=0.0;
	long double dE;
	long double bx, by, bz;
	bx = c.boxx, by = c.boxy, bz = c.boxz;
		
	
	for (unsigned int i = 0; i < s.natoms; ++i)
	{
		//if(i!=Aindex){
			Etot += PairEnergyES(c.atom[Aindex], s.atom[i], bx, by, bz);
			Etoto += PairEnergyESO(s.atom[i], c.atom[Aindex], ox, oy, oz, bx, by, bz);
					
		//}	
	}
	dE = Etot - Etoto;
	//cout << "In DeltaE funct, dE is : " << dE << endl;
	return( dE ); 
	

}
long double TotalEnergyTwoES(const Config& c1, const Config& c2){
	
	long double Etot = 0.0;
	unsigned int na = c1.natoms;
	unsigned int nb = c2.natoms;

		for(unsigned int pp = 0;pp<na;++pp){

			for(unsigned int uu = 0;uu<nb;++uu){
				
				
				Etot += PairEnergyES(c1.atom[pp], c2.atom[uu], c1.boxx, c1.boxy, c1.boxz);
					
			}
		}
//	Etot *= 60.0;
	//Etot/=4.0;
	return(Etot);
}

long double DeltaEtwo(const Config& c, const unsigned& Aindex, const long double& ox, const long double& oy, const long double& oz, const Config& s){
	long double Etot=0.0;	
	long double Etoto=0.0;
	long double dE;
	long double bx, by, bz;
	bx = s.boxx, by = s.boxy, bz = s.boxz;
		
	//cout<<"s.natoms "<<s.natoms<<" bx "<<bx<<" by "<<by<<" bz "<<bz<<endl;
	for (unsigned int i = 0; i < s.natoms; ++i)
	{
		//if(i!=Aindex){
			Etot += PairEnergy(s.atom[i], c.atom[Aindex], bx, by, bz);
			Etoto += PairEnergyO(s.atom[i], c.atom[Aindex], ox, oy, oz, bx, by, bz);
					
		//}	
	}
	dE = Etot - Etoto;
	//cout << "In DeltaE funct, dE is : " << dE << endl;
	return( dE ); 
	

}
