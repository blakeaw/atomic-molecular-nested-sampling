/* Nested Sampling Routine--Collected nested energies.
Configured to run for lj particles with sigma=1.0 and epsilon=1.0.
The Monte Carlo trial move consist of moving a single atom 
 The boundary condition is set to a cubic box with min image periodicity; 
---Contains conditions to modulate step sizes under each new median energy.

Modifications

12-17-14 Creation - adapated from pNested_Energy_TAM_v1.0.cpp

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
#include "/net/uu/nm/cm/bxw109120/NestedSampling/ContinuumSphere/code/src/class/Atom.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/ContinuumSphere/code/src/class/Config.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/ContinuumSphere/code/src/class/RunningStats.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/ContinuumSphere/code/src/class/HistogramNS_slim.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/ContinuumSphere/code/src/class/Histogram.h"
#include "/net/uu/nm/cm/bxw109120/NestedSampling/ContinuumSphere/code/src/class/ContinuumSphere.h"

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
long double ConfProb(long double E, long double Beta);
long double TestAreaIso(const Config& c, const ContinuumSphere&, long double scale, long double Ec);
int main(int argc, char *argv[])
{

	//Get the local time for the logfile
	time_t rawtime;
 	struct tm * timeinfo;
 	time (&rawtime);
 	timeinfo = localtime (&rawtime);

	// Parameters that need to be set

	// the number of atoms to have in the system
	const unsigned int natoms = 100;
	
	//Temp bound
	long double Tstart = 1.4;
	long double Tend = 2.0;
	long double Tloop = 4;
	//Total number of MC style trial moves under each NS iteration
	const unsigned int MCit = 8000000;
	//Number of procs to split Sampling loop (MCit)
	const unsigned short nprocs = 1;
	// Number of initial MC steps to allow equilibration
	const unsigned Eqit = 80000*natoms;
	
	// Step interval to collect data after equilibration
	const unsigned Cit = natoms*200;
	//Box sizes
	long double boxx = 12.0;//*1.1;
	long double boxy = 12.0;//*1.1;
	long double boxz = 12.0;///(1.1*1.1);
	// Test area scaling value
	long double scale = 1.0001;
	long double sca = 1.01;

	

	const long double bsfac = 0.900;
 	
    long double stepx = boxx*bsfac;
	long double stepy = boxy*bsfac;
	long double stepz = boxz*bsfac;
	//Sphere Quantities
	long double sdensity = 0.3;
	long double seps = 1.0;
	long double srad = 2.0;

	//Pi
	const long double pi = 3.14159265359;
	
	//Enthalpy bounds
	const long double Etop = 200.0;
//	const long double Elow = -200.0;
//	long double Binwidth = 1.0e-4;

		//run descriptor
	string nds = itos(natoms);
	string rundescript = "MC_Iso_n"+nds+"_a2";
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
	//string ehistoutname = "zLDhist_"+rundescript+".dat";
	
	ofstream logfile(logfilename.c_str());
	ofstream houtfile(houtname.c_str());
	ofstream obsoutfile(outfilename.c_str());
	//string coutfilename = "Coord_NS_Ent_a0.xyz";
	ofstream coutfile(coutname.c_str());
	//ofstream eoutfile(ehistoutname.c_str());
//	ofstream stepoutfile("zStepOut_int_test_a1.dat");	
	// configs from each nested sampling iteration
	 
	unsigned long long MCitpp = (MCit+Eqit*nprocs)/nprocs;
	unsigned int ConfigsOutFreq = ((MCitpp-Eqit))/2;

cout << "Temperature iterations is set to " <<Tloop  << " from Temps "<<Tstart<<" - "<<Tend<<endl;
cout<<" MC iterations set to " << MCit << endl;
//cout<<"Etop: "<<Etop<<" Elow: "<<Elow<<endl;	

logfile << " " << endl;
	logfile << "Local time and date: " << asctime(timeinfo) << endl;
	logfile<<"Constant Volume and Constant Temperature Nested Sampling of lennard jones system with continuum sphere fixed at origin"<<endl;
	logfile << "Parameters: " << endl;
	logfile<<" Number of Trial Moves (MCit): "<<MCit<<endl;
	//logfile<<"NS fraction : "<<NSfrac<<endl;
//	logfile<<"TAM scaling factor: "<<scale<<endl;
	logfile<<"natoms: "<<natoms<<" nprocs: "<<nprocs<<" Eqit: "<<Eqit<<" Cit: "<<Cit<<endl;
	logfile<<" boxx: "<<boxx<<" boxy: "<<boxy<<" boxz: "<<boxz<<endl;
	logfile<<"bsfac: "<<bsfac<<" initial guess stepx: "<<stepx<<" stepy: "<<stepy<<" "<<stepz<<endl;
	//logfile<<" Etop: "<<Etop<<" Elow: "<<Elow<<endl;	
	logfile<<"MCitpp: "<<MCitpp<<"ConfigsOutFreq: "<<ConfigsOutFreq<<endl;
	//logfile<<"Temperature: "<<Temp<<endl;
	// Initialize the MT RNG object
	dsfmt_t dsfmt;
	
	// Need to initialize member mti as NN+1
	//Random1.mti = NN+1;
	//long double Beta = 1.0/Temp;
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
	RunningStats Estat;
	RunningStats Esqstat;
	RunningStats LDstat;
	RunningStats dUstat;
	RunningStats dUsstat;
	RunningStats dUcstat;
//	RunningStats Rgstat;
	Histogram Eohist;
	long double ldlow = 0.0;
	long double ldhigh = 5.0;
	//Eohist.Initialize(Binwidth, ldlow, ldhigh);
	
//logfile<<"number of NS bins: "<<Ehist.GetNumberOfBins()<<endl;
	//logfile<<"Potential Energy Histogram Params: "<<endl;
	//logfile<<"number of bins: "<<Eohist.GetNumberOfBins()<<" bindwidth: "<<Eohist.GetIncrement()<<" Elow: "<<Elow<<" Etop: "<<Etop<<endl;

	ContinuumSphere CS(seps, srad, sdensity);
	logfile<<"Cont Sphere parms: "<<endl;
	logfile<<"Epsilon: "<<seps<<" density: "<<sdensity<<" radius "<<srad<<endl;

//Initialize configuration objects
	
Config trial1(natoms);
Config Bottom(natoms);
Bottom.SetBox(boxx, boxy, boxz);

trial1.SetBox(boxx, boxy, boxz);
//Config Lconf(natoms);
//Lconf.SetBox(boxx, boxy, boxz);


cout<<"Assigning Initial coordinates.."<<endl;
logfile<<"Assigning Initial coordinates.."<<endl;
//Randomly place atoms in box - but filter to configuration with Enthalpy below Htop

long double Econf;

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
//		if (track<2*natoms){
//			//cout<<"checking distance..."<<endl;
//			for (unsigned int j = 0; j < i; ++j)
//			{
//				long double rdist = RDist(trial1.atom[i], trial1.atom[j]);
//				//cout<<"distance from atom "<<i<<" to atom "<<j<<" is "<<rdist<<endl;
//				if (rdist<0.8){
//				//	cout<<"rdist less than closest allowed"<<endl;
//					flag=1;
//					break;
//				}
//			}
//		}

		if (track<2*natoms){
			long double dcs = CS.Dist(trial1.atom[i]);
			if (dcs<srad){
				flag=1;
			}
		}
		if (flag){
			--i;
			++track;
			//cout<<"reseting i back from "<<i+1<<" to "<<i<<" and track is "<<track<<endl;
		}
	}
	Econf = TotalEnergy(trial1)+CS.TotalPotential(trial1);

	//cout<<Econf<<endl;
	
} while (Econf>100.0);

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
Bottom.Equate(trial1);
	
long double Egmin = TotalEnergy(trial1)+CS.TotalPotential(trial1);

long double Et1 = Egmin;
long double dU = 0.0;
// Output starting coordinates
//coutfile << trial1.natoms << endl;
//coutfile <<"lj "<<natoms<<" E: "<<Egmin<<endl;
//for(unsigned long hg = 0; hg<trial1.natoms; ++hg){
//		coutfile << setiosflags(ios::fixed) << setprecision(15) << "lj " << Bottom.atom[hg].x << " " << Bottom.atom[hg].y << " " << Bottom.atom[hg].z << endl; 
//		}
long double Enest = Etop;
long double Elowprev = Et1;
//long double Elconf = Et1;
long double dT;
if (Tloop>1){
	dT = (Tend - Tstart)/ (long double)(Tloop-1);
}
else{
	dT=0.0;
}
unsigned long long cfailcount = 0;
long double Enprev=Enest;
long double Gvol = pow(boxx, 3)-CS.ComputeVolume();
cout<<"Gas volume: "<<Gvol<<endl;
long double Tprev;
for(unsigned int ii=0;ii<Tloop;++ii){
cout << "Iteration " << ii << "......." << endl;
logfile << "Iteration " << ii << "......." << endl;
long double Temp = Tstart + dT*(long double)(ii);
long double Beta = 1.0/Temp;
cout<<"Temperature : "<<Temp<<endl;
if (ii>0){
	long double fact = Temp/Tprev;
	stepx*=fact;
	stepy*=fact;
	stepz*=fact;
}


cout<<"Initial Energy: "<<Et1<<endl;
logfile<<"Initial Energy: "<<Et1<<endl;

	//unsigned long long acceptcount = 0;
	
	unsigned long long ctries = 0;
	unsigned long long csuccess = 0;

	//long double Elowprev = Et1;
	bool coflag = 1;
	
//	omp_set_num_threads(nprocs);
//	#pragma omp parallel 
//	{

	 //Define a thread local configuration
	 Config Tlocal(natoms);
	 Tlocal.Equate(trial1);
	 Tlocal.SetBox(boxx, boxy, boxz);
	//Config Velocal(natoms);
	// Config Tlocold(natoms);

	//Define thread local variables
	long double Et1local = Et1;
	long double Eelocal = 0.0;

	unsigned long long ctriesloc = 0;
	unsigned long long csuccessloc = 0;

  //Perform markov moves on the configuration
    for(unsigned int n = 0; n<MCitpp; ++n){
		
		
		//Store Initial Et1 value
		long double Et1b = Et1local;
	
		//Trial Move -----
		// randomly select an atom index	
			unsigned int Arandom = (unsigned int)(dsfmt_genrand_close_open(&dsfmt)*natoms);
		
			long double origx, origy, origz;
			origx = Tlocal.atom[Arandom].x;
			origy = Tlocal.atom[Arandom].y;
			origz = Tlocal.atom[Arandom].z;
		
		// Calculate the difference in total potential energy

		long double dE = 0.0;
	
		bool csflag = 0;
		bool tflag = 1;
		//Coordinate change		

			
			long double movex, movey, movez;
			/* set random movements of the coordinates
			   of the selected atom */
			movex = (dsfmt_genrand_close_open(&dsfmt)-0.5)*stepx;
			movey = (dsfmt_genrand_close_open(&dsfmt)-0.5)*stepy;
			movez = (dsfmt_genrand_close_open(&dsfmt)-0.5)*stepz;
			// Implement those moves
			Tlocal.atom[Arandom].x += movex;
			Tlocal.atom[Arandom].y += movey;
			Tlocal.atom[Arandom].z += movez;
			Tlocal.WrapAtom(Arandom);
//			Tlocal.CalcCom();
//			Tlocal.ReCenterZero();
//			Tlocal.WrapCoordinates();
			//Check change in energy for pairwise interactions
			
			long double dEp=0.0;
			long double dEs=0.0;
			csflag = CS.CheckDist(Tlocal.atom[Arandom]);
			//cout<<"dcs "<<dcs<<endl;
			if (csflag){
				dEp=DeltaE(Tlocal, Arandom, origx, origy, origz);
				dEs=CS.DeltaE(Tlocal.atom[Arandom], origx, origy, origz);
				//csflag=1;
				//cout
			}
			//DeltaE(Tlocal, Arandom, origx, origy, origz)+CS.DeltaE(Tlocal.atom[Arandom], origx, origy, origz);
		//	long double Enew = TotalEnergy(Tlocal);
			//dE = Enew - Et1b;
			dE = dEs+dEp;
			//long double Pt1 = ConfProb(Et1b, Beta);
			//Add difference to current total potential energy
			Et1local+=dE;
			//long double Ptrial2 = ConfProb(Et1local, Beta);
			long double Prat = ConfProb(dE, Beta);
			long double Prand = dsfmt_genrand_close_open(&dsfmt);
			if (Prand>Prat){
				tflag=0;
			}
			++ctriesloc;
		
		//Selection Criteria--Metropolis Criteria
		bool tagg = 1;
		
		
		//Add difference to current total potential energy
	//	Et1local+=dE;
		
		
		if(csflag && tflag){
			
			
			tagg=0;
			Eelocal=Et1local;
		
			
					++csuccessloc;
			
			//dUelocal=dUlocal;
//			Tlocal.CalcCom();
//			Tlocal.ReCenterZero();
			if (Et1local<Elowprev){

				#pragma omp critical
				{
					
					Elowprev=Eelocal;
					Bottom.Equate(Tlocal);
					
				}
			}
//			else if (coflag){
//				if (Eelocal>Elowprev && Eelocal<Enest-10.0){
//					Elconf = Eelocal;
//					Lconf.Equate(Tlocal);
//					Lconf.SetBox(boxx, boxx, boxx);
//					coflag=0;
//				}
//			}
			
			
						
			
		}


				/*Failed criteria--r(dsfmt_genrand_close_open(&dsfmt)-0.5)*vstep;eset atom coordinates and remove
		  the difference in energy */ 
		if(tagg == 1){
			
			Et1local = Et1b;
			
			
				Tlocal.atom[Arandom].x = origx;
				Tlocal.atom[Arandom].y = origy;
				Tlocal.atom[Arandom].z = origz;
			
			
			
		}

			if(n>Eqit && n%(Cit)==0){
			//	cout << "!!!!!Going to push values!!!!.................--------------------------!!!!!" << endl;
			//	cout << "Ee: " << Eelocal << " He: " << Helocal << " Ve: " << Velocal << endl;	
			//	cout << "He*He: " << He*He << endl;	
			//	long double LD = LocalDensity(Tlocal);
			//	long double dU = TestAreaIso(Tlocal, CS, scale, Eelocal);
                long double dU = CS.TotalPotentialDerivative(Tlocal);
				long double edU = exp(-Beta*dU);
				#pragma omp critical
				{
				//	Ehist.Push(exp(-Beta*dU));
					Estat.Push(Eelocal);
					Esqstat.Push(Eelocal*Eelocal);
					LDstat.Push(dU);
					dUstat.Push(edU);
					dUsstat.Push(dU);
//					dUcstat.Push(dU*dU*dU);
					//Eohist.Push(LD);
				}
			
		}



			if (n>0 && n%(natoms*1000)==0){
							//check acceptance ratios and adjust
						long double art = (long double)(csuccessloc)/(long double)(ctriesloc);
						cout<<"n "<<n<<" art "<<art<<endl;
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
				coutfile << Tlocal.natoms+1 << endl;
				coutfile <<"lj+CS "<<Tlocal.natoms+1<<" E: "<<Et1local<<endl;
				coutfile<<"CS 0.0 0.0 0.0"<<endl;
				for(unsigned long hg = 0; hg<Tlocal.natoms; ++hg){
					coutfile << setiosflags(ios::fixed) << setprecision(15) << "lj " << Tlocal.atom[hg].x << " " << Tlocal.atom[hg].y << " " << Tlocal.atom[hg].z << endl; 
				}
			}
		//cout<<" -----------------------------"<<endl;

	   }	

	

	//	cout << endl;
	
		
	


		
    } //End Monte Carlo loop

	#pragma omp barrier
	
//	#pragma omp critical
//	{ 
//		
//		csuccess+=csuccessloc;
//		ctries+=ctriesloc;
//			asuccess+=asuccessloc;
//		atries+=atriesloc;

//	}

//}//End parallel

	
//	long double cafrac = (long double)(csuccess)/(long double)(ctries);
//	
//	if (cafrac<0.30){
//		++cfailcount;
//	}

//	cout<<" cfrac: "<<cafrac<<" cfailcount: "<<cfailcount<<endl;
//	
//	logfile<<" cfrac: "<<cafrac<<" cfailcount: "<<cfailcount<<endl;
//	logfile<<" stepx: "<<stepx<<" stepy: "<<stepy<<" stepz: "<<stepz<<endl;
		
//		coutfile << trial1.natoms << endl;
//		coutfile <<"lj 17 H: "<<Htlowprev<<" E: "<<Elowprev<<" V: "<<Vlowprev<<endl;
//		for(unsigned long hg = 0; hg<trial1.natoms; ++hg){
//			coutfile << setiosflags(ios::fixed) << setprecision(15) << "lj " << Bottom.atom[hg].x << " " << Bottom.atom[hg].y << " " << Bottom.atom[hg].z << endl; 
//		}
//	if (cafrac<0.10){
//		cout<<"config frozen. exiting!"<<endl;
//		exit(0);
//	}




long double Eavg = Estat.Mean();
long double Esqavg = Esqstat.Mean();
long double LDavg = LDstat.Mean();
long double dUavg = dUstat.Mean();
long double dUsavg = dUsstat.Mean();
//long double dUcavg = dUcstat.Mean();
//long double Rgavg = Rgstat.Mean();
cout<<" Average Energy is "<<Eavg<<" Average Energy Squared: "<<Esqavg<<endl;
logfile<<"Average Energy is "<<Eavg<<" and Average Energy Squared: "<<Esqavg<<endl;
cout<<"Average dU: "<<LDavg<<endl;
logfile<<"Average dU: "<<LDavg<<endl;
cout<<"Average e^dU: "<<dUavg<<endl;
logfile<<"Average e^dU: "<<dUavg<<endl;
cout<<"Average dUs: "<<dUsavg<<endl;
logfile<<"Average dUs: "<<dUsavg<<endl;
//obsoutfile<<Eavg<<" "<<Esqavg<<" "<<LDavg<<" "<<dUavg<<" "<<dUsavg<<" "<<dUcavg<<endl;
obsoutfile<<Eavg<<" "<<Esqavg<<" "<<LDavg<<" "<<dUavg<<" "<<dUsavg<<" "<<endl;

Estat.Reset();
Esqstat.Reset();
LDstat.Reset();
dUstat.Reset();
dUsstat.Reset();

//gamma
long double Vi = CS.ComputeVolume();
long double sai = CS.ComputeSurfaceArea();
	long double dA = scale*sai-sai;
	long double Rf = sqrt( (dA)/(4.0*pi)+pow(CS.R, 2) );
	long double Vf = (4.0/3.0)*pi*Rf*Rf*Rf;
	long double dV = Vf-Vi;
	long double bxi = boxx;
	long double Vti = pow(bxi, 3);
	long double Vtf = Vti+dV;
	long double Vfrac=Vtf/Vti;
	long double Vterm = (long double)(natoms+1)*Temp*log(Vfrac);
	long double gamma = -Temp*log(dUavg)+Vterm;
	gamma/=dA;
	cout<<"Temp: "<<Temp<<" gamma "<<gamma<<" dA "<<dA<<" Vfrac "<<Vfrac<<" Vterm "<<Vterm<<" dUavg "<<dUavg<<" dV "<<dV<<" gamma/2 "<<gamma/2.0<<endl;
	houtfile<<Temp<<" "<<gamma<<endl;
//dUcstat.Reset();
//Rgstat.Reset();
//if (Elconf<Enest){
//	trial1.Equate(Lconf);
//	Et1 = Elconf;
//	dU = dUconf;
//	cout<<"setting new start conf to Lconf"<<endl;
//}
//else{
	trial1.Equate(Bottom);

	Et1=Elowprev;

	cout<<"setting new start conf to Bottom"<<endl;

	Tprev=Temp;



} // End Temp loop

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
//logfile<<"Energy Histograms : "<<ehistoutname<<endl;
cout<<"Simulation Complete!"<<endl;
cout<<"Output files are: "<<endl;
cout<<"Logfile: "<<logfilename<<endl;
cout<<"Nested Energies: "<<houtname<<endl;
cout<<"Configuration Samples: "<<coutname<<endl;
cout<<"Average Observable Values: "<<outfilename<<endl;
//cout<<"Energy Histograms : "<<ehistoutname<<endl;
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
	long double rc = 2.5*sigma;
	long double rcs = rc*rc;
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
	if (rs<rcs){
		potent = 4.0*eps*((pow(sigma/rs, 6))-(pow(sigma/rs, 3))) - 4.0*eps*((pow(sigma/rcs, 6))-(pow(sigma/rcs, 3)));
	}
	else{
		potent=0.0;
	}	
	
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
	long double rc = 2.5*sigma;
	long double rcs = rc*rc;
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

long double ConfProb(long double E, long double Beta){
//	cout << "Beta " << Beta << endl;
	long double P = exp(-Beta*E);
	return(P);
}
long double TestAreaIso(const Config& c, const ContinuumSphere& s, long double scale, long double Ec){
	Config Test(c);
	ContinuumSphere Tests(s);
	long double pi = 3.14159265359;
	long double Vi = Tests.ComputeVolume();
	long double sai = Tests.ComputeSurfaceArea();
	long double dA = scale*sai-sai;
	long double Rf = sqrt( (dA)/(4.0*pi)+pow(s.R, 2) );
	Tests.R = Rf;
	long double Vf = Tests.ComputeVolume();
	long double dV = Vf-Vi;
	long double bxi = Test.boxx;
	long double Vti = pow(bxi, 3);
	long double Vtf = Vti+dV;
	long double Vfrac=Vtf/Vti;
//	cout<<"dA "<<dA<<" dV "<<dV<<endl;
	Test.ScaleCoordinates(Vfrac);
	long double bsc = pow(Vfrac, 1.0/3.0);
	//bxi*=bsc;
//	cout<<"bsc "<<bsc<<" Vf "<<Vfrac<<endl;
	//Test.SetBox(bxi, bxi, bxi);
	//csflag = CS.CheckDistAll(Tlocal);
//	long double scale = 1.0001;
	//Test.ScaleCoordinates(scale, scale, scale));
	Test.boxx*=bsc;
	Test.boxy*=bsc;
	Test.boxz*=bsc;
	//long double Ec = TotalEnergy(c);
	long double Et = TotalEnergy(Test)+Tests.TotalPotential(Test);
	long double dE = Et-Ec;
	return dE;

}
