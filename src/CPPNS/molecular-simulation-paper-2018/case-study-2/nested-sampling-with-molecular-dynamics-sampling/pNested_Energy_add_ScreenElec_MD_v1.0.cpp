/* Nested Sampling Routine--Collected nested energies.
Configured to run for lj particles with sigma=1.0 and epsilon=1.0.
The Monte Carlo trial move consist of moving a single atom 
 The boundary condition is set to a cubic box with min image periodicity; 
---Contains conditions to modulate step sizes under each new median energy.

Modifications




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
// Values for RNG
#define NN 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define LM 0x7FFFFFFFULL /* Least significant 31 bits */


using namespace std;

 // Classes
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Additive/code/src/class/Atom.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/Additive/code/src/class/Config.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/NSTAM/code/src/class/RunningStats.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/NSTAM/code/src/class/HistogramNS_slim.h"
//#include "/net/uu/nm/cm/bxw109120/NestedSampling/NSTAM/code/src/class/Histogram.h"
#include "/net/uu/nm/cm/bxw109120/src/SimulationCpp/class/Histogram.h"
#include "/net/uu/nm/cm/bxw109120/src/SimulationCpp/class/HardBox.h"
//#include "./src/class/CenterofMass.h"
//MD
class AtomMD {
	public:
		double x, y, z;  // position
		double v_x, v_y, v_z; //velocity (angstrom/fs)
		double f_x, f_y, f_z; //force
		double sigma , epsilon ; // sigma (Angstrom) & epsilon (kcal/mol)
		double mass; // mass (kg/mol)
        bool ctype;
	AtomMD (): sigma(1.0), epsilon(1.0), mass (1.0){ 
                    x=0.0,y=0.0,z=0.0;
                    v_x=0.0,v_y=0.0,v_z=0.0;
                    f_x=0.0,f_y=0.0,f_z=0.0;
               }
           };
           
//class Frame which contains an array of atoms!
class Frame {
	public:
		int natoms; //number of atoms
          double boxx, boxy, boxz; //box dimensions!
          double rCutOff, rSkin; //values of cutoff and skin!
		double t; //time
		AtomMD *atom;
		
//        Frame(int na){
//            natoms=na;
//            atom = new AtomMD[na]; //builds an array of atoms 
//            rCutOff (2.5);
//        }//string eoname = ehistoutname+"_it_"+itos(ii-1)+".dat";
//Histogram Eohist;
//long double Eohistlow=Elist[0];
//long double Eohisthigh=Elist[Nsamples-1];
//if (Eohisthigh>100.0){
//    Eohisthigh=100.0;
//}
//long double Eohistbinw = (Eohisthigh-Eohistlow)/250.0;
//cout<<"Elow: "<<Eohistlow<<" Ehigh: "<<Eohisthigh<<" bw "<<Eohistbinw<<endl;
//Eohist.Initialize(Eohistbinw, Eohistlow, Eohisthigh);
//for (unsigned int i = 0; i < Nsamples; ++i)
//{
//    //Eohist.PushWithWeight(Elist[i], exp(-Beta*Elist[i]);
//    Eohist.PushOverTop(Elist[i]);
//}
//Eohist.Normalize();
//Eohist.DumpNormal(eoname);
//Eohist.Reset();

        Frame(): rCutOff(2.5){};
        
//		Frame(): natoms(n_atoms_n), boxx(5.0), boxy(5.0), boxz(7.0),rCutOff (2.5) {
// 			atom = new AtomMD[natoms]; //builds an array of atoms 
//		}
        
		~Frame(){
			delete [] atom;
		}

        void  InitializeAtoms(int na){
            natoms=na;
            atom = new AtomMD[na]; //builds an array of atoms 
            rCutOff=3.0;
        }

        void SetBox(long double bx, long double by, long double bz){

            boxx=bx,boxy=by,boxz=bz;
            return;
        }
          void equate (Frame &qqq) {
                 for ( int i = 0; i < natoms; i++)
                 { 
                      atom[i].x = qqq.atom[i].x;
                      atom[i].y = qqq.atom[i].y;
                      atom[i].z = qqq.atom[i].z; 
                      atom[i].v_x = qqq.atom[i].v_x;
                      atom[i].v_y = qqq.atom[i].v_y;	
	                 atom[i].v_z = qqq.atom[i].v_z;	
                      atom[i].f_x = qqq.atom[i].f_x;	
                      atom[i].f_y = qqq.atom[i].f_y;
                      atom[i].f_z = qqq.atom[i].f_z;	
                 }
               }
          };


class MTRandomNum {
      public:
	unsigned long long int mt[NN];
        int mti;
	MTRandomNum() : mti(NN+1){}
	// initializes with a seed
        void initialize(unsigned long long int seed){
        	mt[0] = seed;
   		 for (mti=1; mti<NN; mti++){ 
      			  mt[mti] =  (unsigned long long int)(6364136223846793005ULL * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti);
		}
	
		return;
        }
	// called to generate random number
	// on [0,1] 
	long double generate(void){
		 int i;
   		 unsigned long long x;
   		 static unsigned long long mag01[2]={0ULL, MATRIX_A};

    		if (mti >= NN) { /* generate NN words at one time */	

        	/* if initialize has not been called, */
        	/* a default initial seed is used     */
        		if (mti == NN+1){
        		    initialize(5489ULL);
			} 

        		for (i=0;i<NN-MM;i++) {
        		    x = (mt[i]&UM)|(mt[i+1]&LM);
        		    mt[i] = mt[i+MM] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        		}
        		for (;i<NN-1;i++) {
        		    x = (mt[i]&UM)|(mt[i+1]&LM);
        		    mt[i] = mt[i+(MM-NN)] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        		}
        		x = (mt[NN-1]&UM)|(mt[0]&LM);
        		mt[NN-1] = mt[MM-1] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
	
        		mti = 0;
    		}
  
   		 x = mt[mti++];

    		x ^= (x >> 29) & 0x5555555555555555ULL;
    		x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
    		x ^= (x << 37) & 0xFFF7EEE000000000ULL;
    		x ^= (x >> 43);
		// [0,1] real
//		return (x >> 11) * (1.0/9007199244640991.0);
		// [0,1) real
//		return (x >> 11) * (1.0/9007199254740992.0);
		// (0,1) real 
		return ((x >> 12) + 0.5) * (1.0/4503599627370496.0);
	
	}
};
// Prototype Functions
	double PotentialEnergy (Frame &qqq);
	void ForceCalculation (Frame &qqq);
	void verlet (Frame &qqq, Frame &nnn, double timestep);  
     double kinetic (Frame &qqq);
     void pbc (double &dx, double &dy, double &dz, double boxx, double boxy, double boxz);
    void correct (Frame &qqq);
     double gauss (MTRandomNum &obj, double x);
     void printlog (Frame &qqq, int i);
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
long double TotalEnergyES(const Frame& c);
long double PairEnergyES(const AtomMD& one, const AtomMD& two, long double bx, long double by, long double bz);
void AddBiasingForceES(const Frame& c, long double k, long double Uo, long double Uc);
void AddBiasingForceESexpo(const Frame& c, long double Uo);
void PairForceES(const AtomMD& one, const AtomMD& two, long double bx , long double by, long double bz, long double Fv[3]);
void PairForceLJ(const AtomMD& one, const AtomMD& two, long double bx , long double by, long double bz, long double Fv[3]);

int main(int argc, char *argv[])
{

	//Get the local time for the logfile
	time_t rawtime;
 	struct tm * timeinfo;
 	time (&rawtime);
 	timeinfo = localtime (&rawtime);

	// Parameters that need to be set

	// the number of atoms to have in the system
	
	const unsigned int na = 50;
	const unsigned int nb = 50;
	const unsigned int natoms = na+nb;
	// the number of nested sampling iterations to run 
	const unsigned int Iterations = 5;
	const unsigned int Nsamples = 100;
	//Total number of MC style trial moves under each NS iteration
	const unsigned int MCit =2000000;
	//Number of procs to split Sampling loop (MCit)
	const unsigned short nprocs = 1;
	// Number of initial MC steps to allow equilibration
	int Eqit = 1000;
	// Step interval to collect data after equilibration
	const unsigned Cit = natoms*2;
	long double NSfrac = 0.75;
	//Box sizes
	long double boxx = 8.0;//*1.1;
	long double boxy = 8.0;//*1.1;
	long double boxz = 8.2;///(1.1*1.1);
	long double Wd = 16.0;
    long double kharm=1.0;
	// Test area scaling value
	long double scale = 1.00001;
	long double sca = 1.01;
    //how many steps between MD printlogs
    int sout = 5000;
//	boxx*=sca;
//	boxy*=sca;
//	boxz/=(sca*sca);
	long double abox = (boxx+boxy+boxz)/3.0;
	//Translation step
	const long double bsfac = 0.0400;
 	
    long double stepx = boxx*bsfac;
	long double stepy = boxy*bsfac;
	long double stepz = boxz*bsfac;
	
	long double Temp = 0.70;
	//Pi
	//const long double pi = 3.14159265359;
	
	//Enthalpy bounds
	const long double Etop = 100000.0;
	const long double Elow = -1000.0;
	long double Binwidth = 1.0e-3;

		//run descriptor
	string nds = itos(natoms);
	string rundescript = "NS_D_n"+nds+"_binaryLJSE_kappa1_t1_rep6";
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
	string ehistoutname = "zEabhist_"+rundescript;
	
	ofstream logfile(logfilename.c_str());
	ofstream houtfile(houtname.c_str());
	ofstream obsoutfile(outfilename.c_str());
	//string coutfilename = "Coord_NS_Ent_a0.xyz";
	ofstream coutfile(coutname.c_str());
	ofstream eoutfile(ehistoutname.c_str());
//	ofstream stepoutfile("zStepOut_int_test_a1.dat");	
	// configs from each nested sampling iteration
	 
	unsigned long long MCitpp = (MCit+Eqit*nprocs)/nprocs;
	unsigned int ConfigsOutFreq =5000;

cout << "Nested sampling iterations is set to " << Iterations << " with MC iterations set to " << MCit << endl;
cout<<"Etop: "<<Etop<<" Elow: "<<Elow<<endl;	

logfile << " " << endl;
	logfile << "Local time and date: " << asctime(timeinfo) << endl;
	logfile<<"Constant Volume Nested Sampling of binary lennard jones system with charge. Charge is screened electrostatics with kappa=1. The cutoff for all interactions is rc=3.0"<<endl;
	logfile<<" Diffusive sampling with weight: "<<Wd<<endl;
	logfile << "Parameters: " << endl;
	logfile<<"NS Iterations: "<<Iterations<<" Number of Trial Moves (MCit): "<<MCit<<endl;
	logfile<<"NS fraction : "<<NSfrac<<endl;
	logfile<<"TAM scaling factor: "<<scale<<endl;
	logfile<<" na: "<<na<<" nb: "<<nb<<endl;
	logfile<<"natoms: "<<natoms<<" nprocs: "<<nprocs<<" Eqit: "<<Eqit<<" Cit: "<<Cit<<endl;
	logfile<<" boxx: "<<boxx<<" boxy: "<<boxy<<" boxz: "<<boxz<<endl;
	logfile<<"bsfac: "<<bsfac<<" initial guess stepx: "<<stepx<<" stepy: "<<stepy<<" "<<stepz<<endl;
	logfile<<" Etop: "<<Etop<<" Elow: "<<Elow<<endl;	
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
	HistogramNS Ehist;
    //Histogram Ehist;
	Ehist.Initialize(Binwidth, Elow, 300.0);
	Ehist.SetNSfrac(NSfrac);
	RunningStats Estat;
	RunningStats Esqstat;
	RunningStats LDstat;
	RunningStats dUstat;
	RunningStats dUsstat;
	RunningStats dUcstat;
//	RunningStats Rgstat;
	Histogram Eohist;
	long double ldlow = -200.0;
	long double ldhigh = 200.0;
	Eohist.Initialize(Binwidth*100.0, ldlow, ldhigh);
	Histogram Eohistfull;
	Eohistfull.Initialize(Binwidth*100.0, ldlow, ldhigh);

	logfile<<"number of NS bins: "<<Ehist.GetNumberOfBins()<<endl;
	logfile<<"Potential Energy Histogram Params: "<<endl;
	logfile<<"number of bins: "<<Eohist.GetNumberOfBins()<<" bindwidth: "<<Eohist.GetIncrement()<<" Elow: "<<Elow<<" Etop: "<<ldhigh<<endl;
    
//    Histogram Eshist;
//	long double sdlow = -200.0;
//	long double sdhigh = 200.0;
//	Eshist.Initialize(Binwidth*100.0, ldlow, ldhigh);
//	Histogram Eohistfull;
//	Eohistfull.Initialize(Binwidth*100.0, ldlow, ldhigh);

long double Elist[Iterations]={0.0};
Elist[0]=Etop;
unsigned Nelist=1;

//string eoname = ehistoutname+"_it_"+itos(ii-1)+".dat";
//Histogram Eohist;
//long double Eohistlow=Elist[0];
//long double Eohisthigh=Elist[Nsamples-1];
//if (Eohisthigh>100.0){
//    Eohisthigh=100.0;
//}
//long double Eohistbinw = (Eohisthigh-Eohistlow)/250.0;
//cout<<"Elow: "<<Eohistlow<<" Ehigh: "<<Eohisthigh<<" bw "<<Eohistbinw<<endl;
//Eohist.Initialize(Eohistbinw, Eohistlow, Eohisthigh);
//for (unsigned int i = 0; i < Nsamples; ++i)
//{
//    //Eohist.PushWithWeight(Elist[i], exp(-Beta*Elist[i]);
//    Eohist.PushOverTop(Elist[i]);
//}
//Eohist.Normalize();
//Eohist.DumpNormal(eoname);
//Eohist.Reset();
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
Config Lconfa(na);
Lconfa.SetBox(boxx, boxy, boxz);
Config Lconfb(nb);
Lconfb.SetBox(boxx, boxy, boxz);

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

long double Econf=Etop+100.0;
long double Econf2=Etop+100.0;
long double Econf3=Etop+100.0;
long double Econf4=Etop+100.0;
while (Econf3+Econf4+Econf2>Etop)
{

	//set a
	unsigned track=0;
	for (unsigned int i = 0; i < na; ++i)
	{
		long double rx = (dsfmt_genrand_close_open(&dsfmt)-0.5)*boxx*0.99;
		long double ry = (dsfmt_genrand_close_open(&dsfmt)-0.5)*boxy*0.99;
		long double rz = ((dsfmt_genrand_close_open(&dsfmt))*(-boxz/2.0)*0.90 )-0.1;
		trial1a.atom[i].x = rx;
		trial1a.atom[i].y = ry;
		trial1a.atom[i].z = rz;
		bool flag=0;
		//cout<<"assigned atom: "<<i<<endl;
		if (track<1000*na){
			//cout<<"checking distance..."<<endl;
			for (unsigned int j = 0; j < i; ++j)
			{
				long double rdist = RDistMinImage(trial1a.atom[i], trial1a.atom[j], boxx, boxy, boxz);
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
	//set b
	track=0;
	for (unsigned int i = 0; i < nb; ++i)
	{
		long double rx = (dsfmt_genrand_close_open(&dsfmt)-0.5)*boxx*0.99;
		long double ry = (dsfmt_genrand_close_open(&dsfmt)-0.5)*boxy*0.99;
		long double rz = ((dsfmt_genrand_close_open(&dsfmt))*(boxz/2.0)*0.9)+0.1;
		trial1b.atom[i].x = rx;
		trial1b.atom[i].y = ry;
		trial1b.atom[i].z = rz;
		bool flag=0;
		//cout<<"assigned atom: "<<i<<endl;
		if (track<1000*nb){//string eoname = ehistoutname+"_it_"+itos(ii-1)+".dat";
//Histogram Eohist;
//long double Eohistlow=Elist[0];
//long double Eohisthigh=Elist[Nsamples-1];
//if (Eohisthigh>100.0){
//    Eohisthigh=100.0;
//}
//long double Eohistbinw = (Eohisthigh-Eohistlow)/250.0;
//cout<<"Elow: "<<Eohistlow<<" Ehigh: "<<Eohisthigh<<" bw "<<Eohistbinw<<endl;
//Eohist.Initialize(Eohistbinw, Eohistlow, Eohisthigh);
//for (unsigned int i = 0; i < Nsamples; ++i)
//{
//    //Eohist.PushWithWeight(Elist[i], exp(-Beta*Elist[i]);
//    Eohist.PushOverTop(Elist[i]);
//}
//Eohist.Normalize();
//Eohist.DumpNormal(eoname);
//Eohist.Reset();
			//cout<<"checking distance..."<<endl;
			for (unsigned int j = 0; j < i; ++j)
			{
				long double rdist = RDistMinImage(trial1b.atom[i], trial1b.atom[j], boxx, boxy, boxz);
				//cout<<"distance from atom "<<i<<" to atom "<<j<<" is "<<rdist<<endl;
				if (rdist<0.99){
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
	trial1a.WrapCoordinates();
	trial1b.WrapCoordinates();
	Econf = TotalEnergyES(trial1a)+TotalEnergyES(trial1b)+TotalEnergyTwoES(trial1a, trial1b);
	Econf3 = TotalEnergy(trial1a);
	Econf4 = TotalEnergy(trial1b);
	Econf2=TotalEnergyTwo(trial1a, trial1b);
	cout<<Econf<<" "<<Econf2<<" "<<Econf3<<" "<<Econf4<<endl;
//	if (Econf2>Etop || Econf3>Etop || Econf4>Etop){
//		Econf=Etop+100.0;
//	}
	
} 

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
cout<<"Initial coordinates assigned!"<<endl;
logfile<<"Initial coordinates assigned!"<<endl;
Bottoma.FullEquate(trial1a);
Bottomb.FullEquate(trial1b);	
long double Egmin = TotalEnergyES(Bottoma)+TotalEnergyES(Bottomb)+TotalEnergyTwoES(Bottoma, Bottomb);
//Egmin*=0.1;
cout<<"Egmin "<<Egmin<<endl;
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
long double Enest = Etop+5000.0;
//long double Enest = -5.0;
long double Elowprev = Et1;
long double Elconf = Et1;
unsigned long long cfailcount = 0;
long double Enprev=Enest;
cout<< setiosflags(ios::fixed) << setprecision(10);

//benchmarking 
RunningStats Timer;

  //Initialize MD Copy Configs into an MD frame
      Frame *frame;
      frame = new Frame[2]; //builds 2 frames! one as current frame and the other one as next frame!
       frame[0].InitializeAtoms(natoms);
        frame[1].InitializeAtoms(natoms);
        frame[0].SetBox(boxx, boxy, boxz);
        frame[1].SetBox(boxx, boxy, boxz);
      double timestep = 0.001000; //timestep (femtosecond)
//using the random number generator!
	MTRandomNum random;
	unsigned long long seed1 = time(0); //returns a new random number each time we run the code!
	random.initialize(seed1);
for(int ii=0;ii<Iterations;++ii){
cout << "Iteration " << ii << "......." << endl;
logfile << "Iteration " << ii << "......." << endl;
//if (ii==0){
//	Eqit=200000;
//}
//else{
//	Eqit=0;
//}

 if(ii>0){

//Eqit=1000000; 
//Eqit/=5;
Enest=Ehist.NestedEnergy();
//Ehist.Renormalize();
//Enest=Ehist.IntegrateToFraction(NSfrac);
Elist[Nelist]=Enest;
++Nelist;

if (ii>1){
    
    cout<<"kharm was: "<<kharm<<endl;
//    //recompute the value of kharm 
    long double bprob = 0.05000;
    long double deltaEm = Enest-Elist[Nelist-2];
    kharm = -(2.0*log(bprob))/( Beta*deltaEm*deltaEm);
    //kharm=0.5;
    cout<<"kharm adjusted to :"<<kharm<<endl;
}
 

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
cout<<"Average Local Density: "<<LDavg<<endl;
logfile<<"Average Local Density: "<<LDavg<<endl;
cout<<"Average dU: "<<dUavg<<endl;
logfile<<"Average dU: "<<dUavg<<endl;
obsoutfile<<Eavg<<" "<<Esqavg<<" "<<LDavg<<" "<<dUavg<<" "<<dUsavg<<" "<<dUcavg<<endl;
if ( (ii-1) % 5 ==0){
    

string eoname = ehistoutname+"_it_"+itos(ii-1)+".dat";
Eohist.Normalize();
Eohist.DumpNormal(eoname);


string eonamef = ehistoutname+"_full_it_"+itos(ii-1)+".dat";
Eohistfull.Normalize();
Eohistfull.DumpNormal(eonamef);

}
Eohist.Reset();
Eohistfull.Reset();
Ehist.SmartReset();
//Ehist.Reset();
//Hhist.Reset();
Estat.Reset();
Esqstat.Reset();
LDstat.Reset();
dUstat.Reset();
dUsstat.Reset();
dUcstat.Reset();
//Rgstat.Reset();
if (Elconf<Enest){
	trial1a.FullEquate(Lconfa);
	trial1b.FullEquate(Lconfb);
	Et1 = Elconf;
	cout<<"setting new start conf to Lconf"<<endl;
}
else{
	trial1a.FullEquate(Bottoma);
	trial1b.FullEquate(Bottomb);
	Et1=Elowprev;
	cout<<"setting new start conf to Bottom"<<endl;
}

//trial1.SetBox(boxx, boxy, boxz);
abox = (boxx+boxy+boxz)/3.0;

stepx = bsfac*boxx*pow(0.90, cfailcount);
stepy = bsfac*boxy*pow(0.90, cfailcount);
stepz = bsfac*boxz*pow(0.90, cfailcount);
if (abs(Enest-Enprev)<0.01){
	break;
}
Enprev=Enest;
houtfile<<Enest<<endl;

}
cout<<"Eqit "<<Eqit<<endl;

cout<<"Initial Energy: "<<Et1<<endl;
logfile<<"Initial Energy: "<<Et1<<endl;
	//start time for benchmark	 
	clock_t startTime = clock();
	//unsigned long long acceptcount = 0;
	
	unsigned long long ctries = 0;
	unsigned long long csuccess = 0;
	//long double Elowprev = Et1;
	bool coflag = 1;
	long double Emax=-1000.0;
	//omp_set_num_threads(nprocs);
	#pragma omp parallel 
	{

	 //Define a thread local configuration
	 Config Tlocala(na);
	 Tlocala.FullEquate(trial1a);
	 Tlocala.SetBox(boxx, boxy, boxz);
	Tlocala.WrapCoordinates();
	 Config Tlocalb(nb);
	 Tlocalb.FullEquate(trial1b);
	 Tlocalb.SetBox(boxx, boxy, boxz);
	Tlocalb.WrapCoordinates();
	// Config Tlocold(natoms);
	
  
      frame[0].t = 0.0; // sets initial time to zero!
        for (unsigned int i = 0; i < na; ++i)
        {
            frame[0].atom[i].ctype=0;
            frame[1].atom[i].ctype=0;
            frame[0].atom[i].x = Tlocala.atom[i].x;
            frame[0].atom[i].y = Tlocala.atom[i].y;
            frame[0].atom[i].z = Tlocala.atom[i].z;
            frame[0].atom[i].v_x=0.0;
            frame[0].atom[i].v_y=0.0;
            frame[0].atom[i].v_z=0.0;
            //cout<<"Tlocal a x "<<Tlocala.atom[i].x<<" frame "<<frame[0].atom[i].x<< " i "<<i<<endl;
        }
        for (unsigned int i = na; i < natoms; ++i)
        {
             frame[0].atom[i].ctype=1;
            frame[1].atom[i].ctype=1;
            frame[0].atom[i].x = Tlocalb.atom[i-na].x;
            frame[0].atom[i].y = Tlocalb.atom[i-na].y;
            frame[0].atom[i].z = Tlocalb.atom[i-na].z;
            frame[0].atom[i].v_x=0.0;
            frame[0].atom[i].v_y=0.0;
            frame[0].atom[i].v_z=0.0;
           // cout<<"Tlocal b x "<<Tlocalb.atom[i-na].x<<" frame "<<frame[0].atom[i].x<< " i "<<i<<" i-na "<<i-na<<endl;
        }
          // exit(0);
        ForceCalculation (frame[0]); //calculates forces for the first frame!
        long double potentialEnergy = PotentialEnergy (frame[0]);
        cout<<"starting frame potential energy: "<<potentialEnergy<<endl;
        long double Eljtot = TotalEnergy(Tlocala)+TotalEnergy(Tlocalb)+TotalEnergyTwo(Tlocala, Tlocalb);
        cout<<"starting frame potential energy via MC: "<<Eljtot<<endl;
      //
	//Define thread local variables
	//long double Et1local = Et1;
	long double Et1local = TotalEnergyES(Tlocala)+TotalEnergyES(Tlocalb)+TotalEnergyTwoES(Tlocala, Tlocalb);
     //Et1local*=0.1;
	cout<<"Starting ES energy is "<<Et1local<<" from frame0 "<<TotalEnergyES(frame[0])<<endl;
	long double Eelocal = 0.0;
	long double Eljcrossloc = TotalEnergyTwo(Tlocala, Tlocalb);	
	cout<<"Starting LJcross energy "<<Eljcrossloc<<endl; 
	unsigned long long ctriesloc = 0;
	unsigned long long csuccessloc = 0;
	unsigned int nsamp = 0;
	unsigned long int n = 0;
    long double bias=0.0;
        //update the forces with the biasing potential
        if (ii>0){
           
            if (Et1local>Enest){
//               bias = -kharm*(Et1local-Enest)/(3.0*(long double)(natoms));
//                 for (unsigned int i = 0; i < natoms; ++i)
//             {
//                 frame[0].atom[i].f_x+=bias;
//                 frame[0].atom[i].f_y+=bias;
//                 frame[0].atom[i].f_z+=bias;
//             } 
                AddBiasingForceES(frame[0], kharm, Enest, Et1local);
              //  AddBiasingForceESexpo(frame[0], Enest);
            }
             
                     
        }
    long double Et1local1 = Et1local;
  //Perform markov moves on the configuration
    while(nsamp<Nsamples){
		
//		long double term=0.0;
//        if(n<Eqit){
//           term=(long double)(Eqit)-(long double)(n);

//        }
       // kharm=25.0*exp(-term/ (long double)(Eqit/4));
		//Store Initial Et1 value
		long double Et1b = Et1local;
		long double Eljcb = Eljcrossloc;

        // turns next frame[1] to new current frame[0] after the loop is done!
//          if (n>0)
//          {
//          	frame[0].equate(frame[1]);
//          }
         
        // updates the position, velocity and force for the next frame - verlet thermostat!
        frame[1].t = frame[0].t + timestep; // time evolution
     
    
      double C1, C2, gamma, param, beta;
      gamma = 0.1;
      beta = (1.0/Temp);
      param = -gamma * (timestep/2.0);
      C1 = exp (param);
      C2 = sqrt ((1.0 - pow (C1, 2.0)) * (frame[0].atom[0].mass/beta));
   // 1st Step
     for ( int i = 0; i < frame[0].natoms; i++)
         {
    	    frame[0].atom[i].v_x = C1 * frame[0].atom[i].v_x + (C2/frame[0].atom[i].mass) * gauss(random, 1.0);
         frame[0].atom[i].v_y = C1 * frame[0].atom[i].v_y + (C2/frame[0].atom[i].mass) * gauss(random, 1.0);
         frame[0].atom[i].v_z = C1 * frame[0].atom[i].v_z + (C2/frame[0].atom[i].mass) * gauss(random, 1.0);
      
         frame[1].atom[i].x = frame[0].atom[i].x + (frame[0].atom[i].v_x) * timestep + (pow (timestep, 2.0)/(2.0*frame[0].atom[i].mass)) * frame[0].atom[i].f_x;
         frame[1].atom[i].y = frame[0].atom[i].y + (frame[0].atom[i].v_y) * timestep + (pow (timestep, 2.0)/(2.0*frame[0].atom[i].mass)) * frame[0].atom[i].f_y;
         frame[1].atom[i].z = frame[0].atom[i].z + (frame[0].atom[i].v_z) * timestep + (pow (timestep, 2.0)/(2.0*frame[0].atom[i].mass)) * frame[0].atom[i].f_z;
         }   
  // Updates force
        //correct wraps the coordinates
       correct (frame[1]);
          ForceCalculation (frame[1]);
//         //copy frame1 coords into configs and recompute the screened electro energy
//        for (unsigned int i = 0; i < na; ++i)
//        {

//            Tlocala.atom[i].x = frame[1].atom[i].x;
//            Tlocala.atom[i].y = frame[1].atom[i].y;
//            Tlocala.atom[i].z = frame[1].atom[i].z;
//           
//        }
//        for (unsigned int i = na; i < natoms; ++i)
//        {
//            Tlocalb.atom[i-na].x = frame[1].atom[i].x;
//            Tlocalb.atom[i-na].y = frame[1].atom[i].y;
//            Tlocalb.atom[i-na].z = frame[1].atom[i].z;
//         
//        }
       // Et1local = TotalEnergyES(Tlocala)+TotalEnergyES(Tlocalb)+TotalEnergyTwoES(Tlocala, Tlocalb);
        Et1local=Et1local1;
        Et1local1 = TotalEnergyES(frame[1]);
        if (ii>0){
            //long double bias=0.0;
           // bias=0.0;
            if (Et1local1>Enest){
//               bias = -kharm*(Et1local-Enest)/(3.0*(long double)(natoms));
//                 for (unsigned int i = 0; i < natoms; ++i)
//             {
//                 frame[0].atom[i].f_x+=bias;
//                 frame[0].atom[i].f_y+=bias;
//                 frame[0].atom[i].f_z+=bias;
//             } 
                //AddBiasingForceES(frame[0], kharm, Enest);
                
                AddBiasingForceES(frame[1], kharm, Enest, Et1local1);
            }
             
                     
        }
     
       
     for (int i=0; i<frame[0].natoms; i++)
         {
          frame[1].atom[i].v_x = frame[0].atom[i].v_x + (timestep/2.0) * (frame[0].atom[i].f_x + frame[1].atom[i].f_x);
          frame[1].atom[i].v_y = frame[0].atom[i].v_y + (timestep/2.0) * (frame[0].atom[i].f_y + frame[1].atom[i].f_y);
          frame[1].atom[i].v_z = frame[0].atom[i].v_z + (timestep/2.0) * (frame[0].atom[i].f_z + frame[1].atom[i].f_z);
      
          frame[1].atom[i].v_x = C1 * frame[1].atom[i].v_x + (C2/frame[1].atom[i].mass) * gauss(random, 1.0);
          frame[1].atom[i].v_y = C1 * frame[1].atom[i].v_y + (C2/frame[1].atom[i].mass) * gauss(random, 1.0);
          frame[1].atom[i].v_z = C1 * frame[1].atom[i].v_z + (C2/frame[1].atom[i].mass) * gauss(random, 1.0);
          //reset the next frame zero
          frame[0].atom[i].x = frame[1].atom[i].x;
          frame[0].atom[i].y = frame[1].atom[i].y;
          frame[0].atom[i].z = frame[1].atom[i].z; 
          frame[0].atom[i].v_x = frame[1].atom[i].v_x;
          frame[0].atom[i].v_y = frame[1].atom[i].v_y;	
          frame[0].atom[i].v_z = frame[1].atom[i].v_z;	
          frame[0].atom[i].f_x = frame[1].atom[i].f_x;	
          frame[0].atom[i].f_y = frame[1].atom[i].f_y;
          frame[0].atom[i].f_z = frame[1].atom[i].f_z;


        }      
//          // verlet (frame[0], frame[1], timestep);
//         frame[1].t = frame[0].t + timestep; // time evolution
//   // 1st Step
//     for ( int i = 0; i < frame[0].natoms; i++)
//       {
//    	   frame[1].atom[i].x = frame[0].atom[i].x + frame[0].atom[i].v_x * timestep + (pow (timestep, 2.0)/(2.0 * frame[0].atom[i].mass)) * frame[0].atom[i].f_x;
//           frame[1].atom[i].y = frame[0].atom[i].y + frame[0].atom[i].v_y * timestep + (pow (timestep, 2.0)/(2.0 * frame[0].atom[i].mass)) * frame[0].atom[i].f_y;
//           frame[1].atom[i].z = frame[0].atom[i].z + frame[0].atom[i].v_z * timestep + (pow (timestep, 2.0)/(2.0 * frame[0].atom[i].mass)) * frame[0].atom[i].f_z;
//       //   cout<<"i "<<i<<" 1st step frame1 x: "<<frame[1].atom[i].x<<" frame0 x "<<frame[0].atom[i].x<<" v_x "<<frame[0].atom[i].v_x<<" f_x "<<frame[0].atom[i].f_x<<endl;
//   // Half-step
//           frame[1].atom[i].v_x = frame[0].atom[i].v_x + (timestep / (2.0 * frame[0].atom[i].mass)) * frame[0].atom[i].f_x;
//           frame[1].atom[i].v_y = frame[0].atom[i].v_y + (timestep / (2.0 * frame[0].atom[i].mass)) * frame[0].atom[i].f_y; 
//           frame[1].atom[i].v_z = frame[0].atom[i].v_z + (timestep / (2.0 * frame[0].atom[i].mass)) * frame[0].atom[i].f_z;
//     
//         }
//        // exit(0);
//  // Updates force
//       correct(frame[1]);
//       ForceCalculation (frame[1]);
//        //copy frame1 coords into configs and recompute the screened electro energy
//        for (unsigned int i = 0; i < na; ++i)
//        {

//            Tlocala.atom[i].x = frame[1].atom[i].x;
//            Tlocala.atom[i].y = frame[1].atom[i].y;
//            Tlocala.atom[i].z = frame[1].atom[i].z;
//           
//        }
//        for (unsigned int i = na; i < natoms; ++i)
//        {
//            Tlocalb.atom[i-na].x = frame[1].atom[i].x;
//            Tlocalb.atom[i-na].y = frame[1].atom[i].y;
//            Tlocalb.atom[i-na].z = frame[1].atom[i].z;
//         
//        }
//        Et1local = TotalEnergyES(Tlocala)+TotalEnergyES(Tlocalb)+TotalEnergyTwoES(Tlocala, Tlocalb);
//        long double bias;
//            if (Et1local>Enest){
//                cout<<"adding bias to forces!"<<endl;
//               bias = kharm*(Et1local-Enest)/(3.0*(long double)(natoms));
//                 for (unsigned int i = 0; i < natoms; ++i)
//             {
//                 frame[1].atom[i].f_x+=bias;
//                 frame[1].atom[i].f_x+=bias;
//                 frame[1].atom[i].f_x+=bias;
//             } 
//            }
// // Second half-step
//     for ( int i = 0; i < frame[0].natoms; i++)
//         {
//  	      frame[1].atom[i].v_x = frame[1].atom[i].v_x + (timestep/(2.0 * frame[0].atom[i].mass)) * frame[1].atom[i].f_x;
//           frame[1].atom[i].v_y = frame[1].atom[i].v_y + (timestep/(2.0 * frame[0].atom[i].mass)) * frame[1].atom[i].f_y;
//           frame[1].atom[i].v_z = frame[1].atom[i].v_z + (timestep/(2.0 * frame[0].atom[i].mass)) * frame[1].atom[i].f_z;
//         } 

		
		
		 //integration complete - copy frame0 coords into configs and recompute the screened electro energy
      //  cout<<"copying coordinates into configs after integration.."<<endl;
//        for (unsigned int i = 0; i < na; ++i)
//        {

//            Tlocala.atom[i].x = frame[0].atom[i].x;
//            Tlocala.atom[i].y = frame[0].atom[i].y;
//            Tlocala.atom[i].z = frame[0].atom[i].z;
//          // cout<<"Tlocal a x "<<Tlocala.atom[i].x<<" frame "<<frame[0].atom[i].x<< " i "<<i<<endl;
//        }
//        for (unsigned int i = na; i < natoms; ++i)
//        {
//            Tlocalb.atom[i-na].x = frame[0].atom[i].x;
//            Tlocalb.atom[i-na].y = frame[0].atom[i].y;
//            Tlocalb.atom[i-na].z = frame[0].atom[i].z;
//         //     cout<<"Tlocal b x "<<Tlocalb.atom[i-na].x<<" frame "<<frame[0].atom[i].x<< " i "<<i<<" i-na "<<i-na<<endl;
//         
//        }
        //wrap coordinates
      //  cout<<"Wrapping coordinates.. "<<endl;
      //  Tlocala.WrapCoordinates();
      //  Tlocalb.WrapCoordinates();
        //compute new energy
      //  cout<<"computing new ES energy..."<<endl;
      //  Et1local = TotalEnergyES(Tlocala)+TotalEnergyES(Tlocalb)+TotalEnergyTwoES(Tlocala, Tlocalb);

      //     Et1local=TotalEnergyES(frame[0]);

	//	cout<<"Et1local "<<Et1local<<" Et1b "<<Et1b<<endl;
       // potentialEnergy = PotentialEnergy (frame[0]);
      //  cout<<"frame potential energy: "<<potentialEnergy<<endl;
      //  Eljtot = TotalEnergy(Tlocala)+TotalEnergy(Tlocalb)+TotalEnergyTwo(Tlocala, Tlocalb);
     //   cout<<"frame potential energy via MC: "<<Eljtot<<endl;
       // exit(0)	;
		if (n%sout==0){
            if (Et1local>Enest){
                bias = kharm*(Et1local-Enest);
            }
		    else{
                bias=0.0;
            }
            printlog(frame[1], n);
             cout<<"Current Screened ES energy: "<<Et1local<<" Current nested value: "<<Enest<<" nsamp "<<nsamp<<endl;
             cout<<"Harmonic bias prefactor: "<<bias<< " Current kharm: "<<kharm<<endl;
		}
//			if (Et1local<Elowprev){

//				#pragma omp critical
//				{
//					
//					Elowprev=Et1local;
//					Bottoma.FullEquate(Tlocala);
//					Bottomb.FullEquate(Tlocalb);
//					//Bottom.SetBox(boxx, boxy, blocal);
//				}
//			}
//			else if (coflag){
//				if (Et1local>Elowprev && Et1local<Enest-5.0){
//					Elconf = Et1local;
//					Lconfa.FullEquate(Tlocala);
//					Lconfb.FullEquate(Tlocalb);
//					coflag=0;
//				}
//			}
			if (n>Eqit && Et1local>Emax){
				Emax=Et1local;
			}
			
			if(n==Eqit){
               cout<<"Equilibration period passed....."<<endl; 
            }	

			if(n>Eqit && ((n-Eqit)%Cit)==0){
			//	cout << "!!!!!Going to push values!!!!.................--------------------------!!!!!" << endl;
			//	cout << "Ee: " << Eelocal << " He: " << Helocal << " Ve: " << Velocal << endl;	
			//	cout << "He*He: " << He*He << endl;	

            Eohistfull.Push(Et1local);
			if (Et1local < Enest){
                //copy from frame to conf object for LD calcs
				for (unsigned int i = 0; i < na; ++i)
                    {

                        Tlocala.atom[i].x = frame[0].atom[i].x;
                        Tlocala.atom[i].y = frame[0].atom[i].y;
                        Tlocala.atom[i].z = frame[0].atom[i].z;
                      // cout<<"Tlocal a x "<<Tlocala.atom[i].x<<" frame "<<frame[0].atom[i].x<< " i "<<i<<endl;
                    }
                    for (unsigned int i = na; i < natoms; ++i)
                    {
                        Tlocalb.atom[i-na].x = frame[0].atom[i].x;
                        Tlocalb.atom[i-na].y = frame[0].atom[i].y;
                        Tlocalb.atom[i-na].z = frame[0].atom[i].z;
                     //     cout<<"Tlocal b x "<<Tlocalb.atom[i-na].x<<" frame "<<frame[0].atom[i].x<< " i "<<i<<" i-na "<<i-na<<endl;
                     
                    }
                    if (Et1local<Elowprev){

				    //#pragma omp critical
				    //{
					
					        Elowprev=Et1local;
					        Bottoma.FullEquate(Tlocala);
					        Bottomb.FullEquate(Tlocalb);
					        //Bottom.SetBox(boxx, boxy, blocal);
				       // }
			        }
			        else if (coflag){
				        if (Et1local>Elowprev && Et1local<Enest-5.0){
					        Elconf = Et1local;
					        Lconfa.FullEquate(Tlocala);
					        Lconfb.FullEquate(Tlocalb);
					        coflag=0;
				        }
			        }
				//cout<<"pushing value : "<<Et1local<<" Enest "<<Enest<<endl;
				long double LD = LocalDensity(Tlocala);
				long double dU = LocalDensityTwo(Tlocala, Tlocalb);
//				long double dU = TestArea(Tlocal, scale, Eelocal);
				//long double LD =0.0;
				//long double dU = 0.0;
				//#pragma omp critical
				//{
					
					Ehist.Push(Et1local);
                	Estat.Push(Et1local);
					Esqstat.Push(Et1local*Et1local);
					LDstat.Push(LD);
					dUstat.Push(dU);
					dUsstat.Push(dU*dU);
					dUcstat.Push(dU*dU*dU);
					Eohist.Push(Et1local);
                   
					++nsamp;
				//	cout<<"energy "<<Et1local<<" accepted for Enest "<<Enest<<" nsamp now "<<nsamp<<endl;
				//}
			}
		}


//		//Output coordinates of sample n
//	   if( (n)%(ConfigsOutFreq)==0) {
//			#pragma omp critical 
//			{
//				//cout<<"Output coordinates to file at trial move: "<<n<<endl;
//				coutfile << natoms << endl;
//				coutfile <<"lj "<<natoms<<" E: "<<Et1local<<endl;
//				for(unsigned long hg = 0; hg<frame[0].natoms; ++hg){
//					coutfile << setiosflags(ios::fixed) << setprecision(10) << frame[0].atom[hg].ctype<<" "<< frame[0].atom[hg].x << " " << frame[0].atom[hg].y << " " << frame[0].atom[hg].z << endl; 
//				}
////				for(unsigned long hg = 0; hg<Tlocalb.natoms; ++hg){
////					coutfile << setiosflags(ios::fixed) << setprecision(15) << "b " << Tlocalb.atom[hg].x << " " << Tlocalb.atom[hg].y << " " << Tlocalb.atom[hg].z << endl; 
////				}
//			}
		//cout<<" -----------------------------"<<endl;

	  // }	

	

	//	cout << endl;
	
		
	


	++n;		
    } //End Monte Carlo loop

	

}//End parallel

	
	
	if (ii==0){
		houtfile<<Emax<<endl;
		cout<<"Maximum Energy found: "<<Emax<<endl;
	}	


	// to compute its execution duration in runtime
	// sampling execution time in seconds
	long double tsec= ( clock() - startTime ) / (long double)(CLOCKS_PER_SEC);
	
	if (ii>0)
	{
		Timer.Push(tsec);
	}
	long double tsecavg = Timer.Mean();
	cout<<"Benchmarking: "<<tsec<<" seconds for sampling loop"<<endl;
	cout<<"Benchmarking: Average second per sampling loop is "<<tsecavg<<endl; 
	long double iterperhour;
	if (ii>1)
	{
		iterperhour= 3600.0/tsecavg;
	}
	else{
		tsecavg=tsec;
		iterperhour = 3600.0/tsec;
	}
 
	cout<<"Benchmarking: "<<iterperhour<<" sampling loops per hour "<<endl;
	logfile<<"Benchmarking: "<<tsec<<" seconds for sampling loop"<<endl;
	logfile<<"Benchmarking: Average second per sampling loop is "<<tsecavg<<endl; 
	logfile<<"Benchmarking: "<<iterperhour<<" sampling loops per hour "<<endl;

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
logfile<<"Energy Histograms : "<<ehistoutname<<endl;
cout<<"Simulation Complete!"<<endl;
cout<<"Output files are: "<<endl;
cout<<"Logfile: "<<logfilename<<endl;
cout<<"Nested Energies: "<<houtname<<endl;
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
	long double rc=2.5;
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
	long double rc = 2.5;
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
	long double rc=2.5;
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
		
	
	for (unsigned int i = 0; i < c.natoms; ++i)
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
	bx = c.boxx, by = c.boxy, bz = c.boxz;
		
	
	for (unsigned int i = 0; i < c.natoms; ++i)
	{
		//if(i!=Aindex){
			Etot += PairEnergy(c.atom[Aindex], s.atom[i], bx, by, bz);
			Etoto += PairEnergyO(s.atom[i], c.atom[Aindex], ox, oy, oz, bx, by, bz);
					
		//}	
	}
	dE = Etot - Etoto;
	//cout << "In DeltaE funct, dE is : " << dE << endl;
	return( dE ); 
	

}

// Potential energy function: returns potential energy!

double PotentialEnergy (Frame &qqq){
   int n;
   n = qqq.natoms;
   double r2, dx, dy, dz, sigma, epsilon, sigma2,rCutOff2, u;

   u = 0.0;
   for (int i = 0; i < n - 1 ; i++)
  	  {
 
             for (int j = i + 1; j < n; j++)
               {  
	            dx = qqq.atom[i].x - qqq.atom[j].x;
                 dy = qqq.atom[i].y - qqq.atom[j].y;
                 dz = qqq.atom[i].z - qqq.atom[j].z;
                 pbc(dx, dy, dz, qqq.boxx, qqq.boxy, qqq.boxz);
    	            r2 = pow (dx , 2.0) + pow (dy, 2.0) + pow (dz , 2.0);
                 rCutOff2 = pow (qqq.rCutOff, 2.0);
              // test cutoff
                 if (r2 < rCutOff2)
                    {
  		    // sigma = (qqq.atom[i].sigma + qqq.atom[j].sigma) / 2.0;
                sigma=1.0;
                if (qqq.atom[i].ctype!=qqq.atom[j].ctype){
                    epsilon=0.1;
                }
                else{
                    epsilon=1.0;
                }
    		     sigma2 = sigma * sigma;
   		    // epsilon = sqrt (qqq.atom[i].epsilon * qqq.atom[j].epsilon);
   		     u += 4.0 * epsilon * (pow ((sigma2 / r2), 6.0) - pow ((sigma2 / r2), 3.0) ) - 4.0 * epsilon * (pow ((sigma2 / rCutOff2), 6.0) - pow ((sigma2 / rCutOff2), 3.0) );
                    }
              }	
  	}		
return u;
 
}

//Force Calculation function

void ForceCalculation (Frame &qqq) {
	int n;
	n = qqq.natoms;
	double r2, dx, dy, dz, sigma, epsilon, sigma2,rCutOff2, u, Fx, Fy, Fz;


 //Dump old force values and reset to zero

   for(int i = 0; i < n; ++i)
         {
	qqq.atom[i].f_x = 0.0;
	qqq.atom[i].f_y = 0.0;
	qqq.atom[i].f_z = 0.0;
         }

//Calculating the new force values

       for (int i = 0; i < n - 1 ; i++)
            {
 
             for (int j = i + 1; j < n; j++)
                 {  
	            dx = qqq.atom[i].x- qqq.atom[j].x;
                 dy = qqq.atom[i].y- qqq.atom[j].y;
                 dz = qqq.atom[i].z- qqq.atom[j].z;
                pbc(dx, dy, dz, qqq.boxx, qqq.boxy, qqq.boxz);
    	           r2 = pow (dx , 2.0) + pow (dy, 2.0) + pow (dz , 2.0);
                 rCutOff2 = pow (qqq.rCutOff, 2.0);
          // test cutoff 
                 if (r2 < rCutOff2)
                    {
    		      //sigma = (qqq.atom[i].sigma + qqq.atom[j].sigma) / 2.0;
                sigma=1.0;
                epsilon=1.0;
                if (qqq.atom[i].ctype!=qqq.atom[j].ctype){
                    epsilon=0.1;
                }
               
   		      sigma2 = sigma * sigma;
//   		      epsilon = sqrt (qqq.atom[i].epsilon * qqq.atom[j].epsilon);

  		      u = 4.0 * epsilon * (12.0 * pow ((sigma2 / r2), 6.0) - 6.0 * pow ((sigma2 / r2), 3.0));
  		      Fx = u * dx * (1.0/r2) ;//* 4.18400e-7;
  		      Fy = u * dy * (1.0/r2) ;//* 4.18400e-7;
  		      Fz = u * dz * (1.0/r2) ;//* 4.18400e-7;
//                if (abs(Fx)>2.0){
//                    long double factor = abs(Fx)/2.0;
//                    Fx/=factor;
//                }
//                   if (abs(Fy)>2.0){
//                    long double factor = abs(Fy)/2.0;
//                    Fy/=factor;
//                }
//                   if (abs(Fz)>2.0){
//                    long double factor = abs(Fz)/2.0;
//                    Fz/=factor;
//                }
  		        qqq.atom[i].f_x +=  Fx;
 	             qqq.atom[j].f_x += -Fx;
 	             qqq.atom[i].f_y +=  Fy;
 	             qqq.atom[j].f_y += -Fy;
 	             qqq.atom[i].f_z +=  Fz;
 	             qqq.atom[j].f_z += -Fz; 
                //cout<<"dx "<<dx<<" Fx "<<Fx<<" atom "<< i<<" f_x "<<qqq.atom[i].f_x<<" j "<<j<<" f_x "<<qqq.atom[j].f_x<<" eps "<<epsilon<< " r "<<sqrt(r2)<<endl;
                   } 
                 }	
              }		
}

// Verlet function: integrates equations of motion!

void verlet (Frame &qqq, Frame &nnn, double timestep) {
       
      nnn.t = qqq.t + timestep; // time evolution
   // 1st Step
     for ( int i = 0; i < qqq.natoms; i++)
       {
    	      nnn.atom[i].x = qqq.atom[i].x + qqq.atom[i].v_x * timestep + (pow (timestep, 2.0)/(2.0 * qqq.atom[i].mass)) * qqq.atom[i].f_x;
           nnn.atom[i].y = qqq.atom[i].y + qqq.atom[i].v_y * timestep + (pow (timestep, 2.0)/(2.0 * qqq.atom[i].mass)) * qqq.atom[i].f_y;
           nnn.atom[i].z = qqq.atom[i].z + qqq.atom[i].v_z * timestep + (pow (timestep, 2.0)/(2.0 * qqq.atom[i].mass)) * qqq.atom[i].f_z;
   // Half-step
           nnn.atom[i].v_x = qqq.atom[i].v_x + (timestep / (2.0 * qqq.atom[i].mass)) * qqq.atom[i].f_x;
           nnn.atom[i].v_y = qqq.atom[i].v_y + (timestep / (2.0 * qqq.atom[i].mass)) * qqq.atom[i].f_y; 
           nnn.atom[i].v_z = qqq.atom[i].v_z + (timestep / (2.0 * qqq.atom[i].mass)) * qqq.atom[i].f_z;
     
         }
  // Updates force
     
       ForceCalculation (nnn);

 // Second half-step
     for ( int i = 0; i < qqq.natoms; i++)
         {
  	      nnn.atom[i].v_x = nnn.atom[i].v_x + (timestep/(2.0 * qqq.atom[i].mass)) * nnn.atom[i].f_x;
           nnn.atom[i].v_y = nnn.atom[i].v_y + (timestep/(2.0 * qqq.atom[i].mass)) * nnn.atom[i].f_y;
           nnn.atom[i].v_z = nnn.atom[i].v_z + (timestep/(2.0 * qqq.atom[i].mass)) * nnn.atom[i].f_z;
         } 
}
  

void pbc (double &dx, double &dy, double &dz, double boxx, double boxy, double boxz){
	int pbc_int (double);
	int xInt = pbc_int(dx/boxx);
	int yInt = pbc_int(dy/boxy);
	int zInt = pbc_int(dz/boxz);

	dx = dx - boxx * (double)(xInt);
	dy = dy - boxy * (double)(yInt);
	dz = dz - boxz * (double)(zInt);
}

// pbc_int function : rounds a real number to the nearest integer!

int pbc_int (double x) {
  int i;
  double y;
  i = (int)(x);
  y = x-(double)(i);
     if (y > 0.5)  {i++;}
     if (y < -0.5) {i--;}      
     return i;
}

// GAUSSIAN FUNCTION -------------------------------------------------------------------------------------------------------------------------------------------------------------//
// Gaussian function: returns a value taken from a Gaussian distribution with zero mean and standard deviation x as input!
double gauss  (MTRandomNum &obj, double x) {

   static int iset = 0;
   static float gset;
   float fac, r, v1, v2;
 
     if (iset ==0) {
                    do {
                         v1=2.0*obj.generate()-1.0;
                         v2=2.0*obj.generate()-1.0;
                         r = v1*v1+v2*v2;
                      } while (r>=1.0 || r == 0.0);

          fac = sqrt (-2.0 * log(r)/r);
          gset = v1*fac;
          iset=1;
          return v2*fac*x;
                 } else {
                     iset=0;
                   return gset*x;

                        }
}
// WRAPPING FUNCTION -------------------------------------------------------------------------------------------------------------------------------------------------------------//
// wrapping up the coordinates
void correct (Frame &qqq){
	
	for (int i = 0; i < qqq.natoms; i++) 
	{
	while (qqq.atom[i].x > qqq.boxx/2.0) {qqq.atom[i].x = qqq.atom[i].x - qqq.boxx;}
	while (qqq.atom[i].y > qqq.boxy/2.0) {qqq.atom[i].y = qqq.atom[i].y - qqq.boxy;}
	while (qqq.atom[i].z > qqq.boxz/2.0) {qqq.atom[i].z = qqq.atom[i].z - qqq.boxz;}
	while (qqq.atom[i].x < -qqq.boxx/2.0) {qqq.atom[i].x = qqq.atom[i].x + qqq.boxx;}
	while (qqq.atom[i].y < -qqq.boxy/2.0) {qqq.atom[i].y = qqq.atom[i].y + qqq.boxy;}
	while (qqq.atom[i].z < -qqq.boxz/2.0) {qqq.atom[i].z = qqq.atom[i].z + qqq.boxz;}	
     }
}

// KINETIC FUNCTION------------------------------------------------------------------------------------------------------------------------------------------------------------//    

double kinetic (Frame &qqq) {
  double k, v2;
  
  k = 0.0;
  for (int i = 0; i<qqq.natoms; i++)
     { 
       v2 = pow (qqq.atom[i].v_x, 2.0) + pow (qqq.atom[i].v_y, 2.0) + pow (qqq.atom[i].v_z, 2.0);
       k += 0.5 * qqq.atom[i].mass * v2;
     }
  return k;
}

// PRINTING DATA FUNCTION --------------------------------------------------------------------------------------------------------------------------------------------------------//

void printlog (Frame &qqq, int i) {

       // ofstream out ("Data.log", ios::app);
        
        double temperature;
	   double potentialEnergy; 	
        double kineticEnergy;  
        double totalEnergy;
        
        
        potentialEnergy = PotentialEnergy(qqq);
        kineticEnergy = kinetic(qqq);
        temperature = (2.0/(3.0 * qqq.natoms)) * kineticEnergy; 
        totalEnergy = kineticEnergy + potentialEnergy;
         
      // prints out the potential energy, kinetic energy and total energy for each timestep to the output file!     
       //out << i << "     " << potentialEnergy << "       " <<  kineticEnergy << "       "  << totalEnergy <<"     "<< temperature << endl; 
       cout <<"step "<< i << " PE " << potentialEnergy << " KE " <<  kineticEnergy << " TE  "  << totalEnergy <<" Temp  "<< temperature << endl; 

}

long double TotalEnergyES(const Frame& c){
	
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

// potential energy function
long double PairEnergyES(const AtomMD& one, const AtomMD& two, long double bx , long double by, long double bz){
	
	
	long double dx, dy, dz, rs;
	long double potent, A;
	long double rc = 2.5;
	long double kappa = 1.0;
	if (one.ctype==two.ctype){
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

void AddBiasingForceES(const Frame& c, long double k, long double Uo, long double Uc){
	

	unsigned int natoms = c.natoms;
     long double pre = k*(Uc-Uo);
      if (abs(pre)>50.0){
                   // cout<<" prefac for bias "<<pre<<" Eij: "<<Eij<<" Uo "<<Uo<<endl;
                    long double factor = abs(pre)/50.0;
                    pre/=factor;
      }
     // pre=1.0;
		for(unsigned int pp = 0;pp<natoms-1;++pp){

			for(unsigned int uu = pp+1;uu<natoms;++uu){
				
				
				//long double Eij = PairEnergyES(c.atom[pp], c.atom[uu], c.boxx, c.boxy, c.boxz);
				long double Fi[3];
                PairForceES(c.atom[pp], c.atom[uu], c.boxx, c.boxy, c.boxz, Fi);
                long double Flj[3];
               // PairForceLJ(c.atom[pp], c.atom[uu], c.boxx, c.boxy, c.boxz, Flj);
              
//                while (abs(pre)>2.0){
//                    
//                    pre*=0.5;
//                }
                long double Fbix,Fbiy,Fbiz,Fbjx,Fbjy,Fbjz;
	            Fbix = pre*Fi[0];
                Fbiy = pre*Fi[1];
                Fbiz = pre*Fi[2];
                Fbjx = pre*(-Fi[0]);
                Fbjy = pre*(-Fi[1]);
                Fbjz = pre*(-Fi[2]);
              //  cout<<" pre "<<pre<<" Fbix "<<Fbix<<" F_ix" <<c.atom[pp].f_x<<" F(LJ)_ix "<<Flj[0]<<endl;
//                c.atom[pp].f_x= 0.0;
//                c.atom[pp].f_y= 0.0;
//                c.atom[pp].f_z= 0.0;
//                c.atom[uu].f_x= 0.0;
//                c.atom[uu].f_y= 0.0;
//                c.atom[uu].f_z= 0.0;

//                c.atom[pp].f_x*= 0.5;
//                c.atom[pp].f_y*= 0.5;
//                c.atom[pp].f_z*= 0.5;
//                c.atom[uu].f_x*= 0.5;
//                c.atom[uu].f_y*= 0.5;
//                c.atom[uu].f_z*= 0.5;

                c.atom[pp].f_x+= Fbix;
                c.atom[pp].f_y+= Fbiy;
                c.atom[pp].f_z+= Fbiz;
                c.atom[uu].f_x+= Fbjx;
                c.atom[uu].f_y+= Fbjy;
                c.atom[uu].f_z+= Fbjz;
			}
		}
//	Etot *= 60.0;
	return;
}
void AddBiasingForceESexpo(const Frame& c, long double Uo){
	

	unsigned int natoms = c.natoms;

		for(unsigned int pp = 0;pp<natoms-1;++pp){

			for(unsigned int uu = pp+1;uu<natoms;++uu){
				
				
				long double Eij = PairEnergyES(c.atom[pp], c.atom[uu], c.boxx, c.boxy, c.boxz);
				long double Fi[3];
                PairForceES(c.atom[pp], c.atom[uu], c.boxx, c.boxy, c.boxz, Fi);
                long double pre = -exp(-Eij/Uo)/Uo;
                long double Fbix,Fbiy,Fbiz,Fbjx,Fbjy,Fbjz;
	            Fbix = pre*(Fi[0]);
                Fbiy = pre*(Fi[1]);
                Fbiz = pre*(Fi[2]);
                Fbjx = pre*(-Fi[0]);
                Fbjy = pre*(-Fi[1]);
                Fbjz = pre*(-Fi[2]);
                c.atom[pp].f_x+= Fbix;
                c.atom[pp].f_y+= Fbiy;
                c.atom[pp].f_z+= Fbiz;
                c.atom[uu].f_x+= Fbjx;
                c.atom[uu].f_y+= Fbjy;
                c.atom[uu].f_z+= Fbjz;
			}
		}
//	Etot *= 60.0;
	return;
}

// potential energy function
void PairForceES(const AtomMD& one, const AtomMD& two, long double bx , long double by, long double bz, long double Fv[3]){
	
	
	long double dx, dy, dz, rs;
	long double potent, A;
    long double Fx, Fy, Fz;
	long double rc = 2.5;
	long double kappa = 1.0;
	if (one.ctype==two.ctype){
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
       long double Fpre = A*exp(-r)*( pow(r,-3.0) + pow(r, -2.0) ); 
		Fx = Fpre*dx;
        Fy = Fpre*dy;
        Fz=Fpre*dz;
	}
	else{
		Fx=0.0;
        Fy=0.0;
        Fz=0.0;
	}	
	
    Fv[0]=Fx;
    Fv[1]=Fy;
    Fv[2]=Fz;
	return;
}

void PairForceLJ(const AtomMD& one, const AtomMD& two, long double bx , long double by, long double bz, long double Fv[3]){
	
	
	long double dx, dy, dz, rs;
	long double potent, A;
    long double Fx, Fy, Fz;
	long double rc = 2.5;
	
	
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
       long double Fpre = 4.0*( -12.0/pow(r, 13) + 6.0/pow(r, 7) ); 
		Fx = Fpre*dx;
        Fy = Fpre*dy;
        Fz=Fpre*dz;
	}
	else{
		Fx=0.0;
        Fy=0.0;
        Fz=0.0;
	}	
	
    Fv[0]=Fx;
    Fv[1]=Fy;
    Fv[2]=Fz;
	return;
}
