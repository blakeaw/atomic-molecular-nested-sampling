#include <stdlib.h>
#include <string>
#include <fstream>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cmath>
//#include "/net/uu/nm/cm/bxw109120/lib/gsl/gsl-1.16/gsl_math.h"
//#include <quadmath.h>
//#include <ctime>

using namespace std;

// Define Classes

// Prototype Functions
double CastDouble(double input);
unsigned long long factorial(unsigned long long n);
unsigned long long Factorial(unsigned long long n);
// Main Function 
int main(int argc, char *argv[])
{
///
cout<<"Launching c++ Program: pNestProc_Energy_TAM_v1.0"<<endl;

///
//Input File Name
string filename = "Epoints_NS_n100_f0.dat";
string filename2 = "zObservables_NS_n100_f0.dat";
//Output File Name
string outfilename = "zNestOut_n100_f0.dat";
// Define the the value of K used in simulation
//int Ksample = 700;
// Define the number of atoms in simulation
long double N = 60.0;
long double V = pow(12.0, 3.0);
//long double bx = 5.0/1.587;
//long double by = 5.0/1.587;
long double bx = 6.7;
long double by = 6.7;
//long double bx = 5.0;
//long double by = 5.0;
long double scale = 1.00001;
//long double sca = 1.01;
//bx*=sca;
//by*=sca;
long double dA = bx*by*(scale*scale - 1.0);
cout<<"Parms set:"<<endl;
cout<<"bx "<<bx<<" by "<<by<<" scale "<<scale<<" dA "<<dA<<endl;
//double Ksampled = CastDouble(Ksample);
// Calculates the alpha based on Ksample
//double alpha = Ksampled/(Ksampled + 1.0);

long double kb = 1.0; 
long double Beta;
long double Temp;

//Planck reduced - material specific
/* argon - 0.185
neon - 0.563
krypton - 0.104
xenon - 0.065
---Ref: Quasi Harmonic Lattice Dynamics and Molecular Dynamics calculations
for the Lennard-Jones solids
*/
long double h = 0.185;
long double m = 1.0;

// Initialize an array to store the Energy values read in from file



//string junk;
ifstream infile;
infile.open(filename.c_str());
// Returns an error if the file does not exist
if (!infile) {
	cerr << "Unable to open file " << filename.c_str() << endl;
	exit(0);
}

long double Hn[200000] = {0};
long double tex, junk;
int i = 0;
cout<<"Reading file: "<<filename.c_str()<<endl;
// Reads the file line by line until end 
while (! infile.eof()){
	//infile>>junk>>tex;
	infile>>tex;
	Hn[i] = tex;
	++i;
	if(i==200000-1){
		break;
	}
}
infile.close();

ifstream infile2;
infile2.open(filename2.c_str());
// Returns an error if the file does not exist
if (!infile2) {
	cerr << "Unable to open file " << filename2.c_str() << endl;
}
const unsigned long ii = i;
long double Obs[ii][2];

//long double tex1, tex2, tex3, tex4, tex5, tex6, tex7;
long double tex1, tex2, tex3, tex4, tex5, tex6;
//long double tex1, tex2, tex3;
//int i = 0;
int kuku=0;
//long double Vlow=1000000000000.0;
cout<<"Reading file: "<<filename2.c_str()<<endl;
// Reads the file line by line until end 
while (! infile2.eof()){
//	infile2 >> junk >> junk >> tex1 >> tex2 >> tex3 >> tex4 >> tex5 >> tex6 >> tex7;
//	infile2 >> tex1 >> tex2 >> tex3 >> tex4 >> tex5 >> tex6;
	infile2 >> tex1 >> tex2;
	Obs[kuku][0] = tex1; //Currently: Mean Energy
	Obs[kuku][1] = tex2; //Currently: Mean Energy Squared
//	Obs[kuku][2] = tex3; //Currently: Local Density
//	Obs[kuku][3] = tex4; //Currently: dU
//	Obs[kuku][4] = tex5; //Currently: dU^2
//	Obs[kuku][5] = tex6; //Currently: dU^3
//	Obs[kuku][6] = tex7; //Currently: Average Dimer-Surface Separation
//	cout<<"t1: "<<tex1<<" t2 "<<tex2<<" t3: "<<tex3<<" t4 "<<tex4<<" t5: "<<tex5<<" t6 "<<tex6<<endl;
//	cout<<" obs[2] "<<Obs[kuku][2]<<endl;
	++kuku;
	if(kuku == i){
		break;
	}
}
infile2.close();
//exit(0);

//exit(0);
long double kT ;
long double Ck;
//Set the number of Temperature points.
// should be a multiple of 10
int Tloop = 4000;
long double Ts = 0.01;
long double Te = 2.0;
long double dT = (Te - Ts)/(long double)(Tloop);
ofstream outfile(outfilename.c_str());
outfile << "# kT T Z U C Ck" << endl;

for(int T=0;T<=Tloop;++T){
	
	Temp = Ts + dT * (long double)(T);
	kT = kb*Temp;
	Beta = 1.0/kT;
	// Estimate the partition function
	int n ;
	//long double np1;
	long double Zest = 0.0 ;
//	long double Zunscale = 0.0;
	long double nd;
	long double He = 0.0;
	long double H;
	long double C;
	long double Ce = 0.0;
	long double Wn;	
	long double Hcurr;
	long double S = 0.0;
	long double lwmax = 100.0;
	long double lw;
	long double Obs1 = 0.0;
	long double Obs2 = 0.0;
//	long double Obs3 = 0.0;
//	long double Obs4 = 0.0;
//	long double Obs5 = 0.0;
//	long double Obs6 = 0.0;
//	long double Obs7 = 0.0;
	long double Obs1curr;
	long double Obs2curr;
//	long double Obs3curr;
//	long double Obs4curr;
//	long double Obs5curr;
//	long double Obs6curr;
//	long double Obs7curr;
	
	for(n=1;n<i-2;++n){
		
		if(n==0){
			Hcurr =Hn[n];
		}
		else{
			Hcurr = 0.50*(Hn[n]+Hn[n-1]);
			Obs1curr = 0.50*(Obs[n][0]+Obs[n-1][0]);
			Obs2curr = 0.50*(Obs[n][1]+Obs[n-1][1]);
//		    Obs3curr = 0.50*(Obs[n][2]+Obs[n-1][2]);
//			Obs4curr = 0.50*(Obs[n][3]+Obs[n-1][3]);
//		    Obs5curr = 0.50*(Obs[n][4]+Obs[n-1][4]);
//			Obs6curr = 0.50*(Obs[n][5]+Obs[n-1][5]);
		  //	cout<<"obs[2] "<<Obs[n][5]<<" obs3cur "<<Obs3curr<<endl;
		}
	//	nd = CastDouble(n);
		//np1 = nd + 1.0;
//		cout<<"n: "<<n<<" obs1curr: "<<Obs1curr<<" obs2curr: "<<Obs2curr<<endl;
//		lw = (-((long double)(n))*log(2.0)) - (Beta*Hcurr);
//		lw = log( pow(0.9, n-1)-pow(0.9, n) )- (Beta*Hcurr);
			
		while(lw>S+lwmax){
			Zest /= exp(lwmax);
			He /= exp(lwmax);
			Ce /= exp(lwmax);
			Obs1 /= exp(lwmax);
			Obs2 /= exp(lwmax);
//			Obs3 /= exp(lwmax);
//			Obs4 /= exp(lwmax);
//			Obs5 /= exp(lwmax);
//			Obs6 /= exp(lwmax);
		//	Obs7 /= exp(lwmax);	
			
			S+=lwmax;	
			
		}
//		long double Wnunscale=exp(lw);
//		Zunscale+=Wnunscale;
		lw += -S;
		
		Wn = exp(lw);	
		//cout << "Wn: " << Wn << endl;
		Zest +=  Wn;
		He += Wn*Hcurr;
		Ce += Wn*Hcurr*Hcurr;
		
		Obs1 += Wn*Obs1curr;
		Obs2 += Wn*Obs2curr;
//		Obs3 += Wn*Obs3curr;
//		Obs4 += Wn*Obs4curr;
//		Obs5 += Wn*Obs5curr;
//		Obs6 += Wn*Obs6curr;
//	//	Obs7 += Wn*Obs7curr;
		//Wn = pow(2, -(n+1));
		//Zest +=  Wn* exp(-Beta*Ecurr);
		//Ue += Wn * (Ecurr * exp(-Beta*Ecurr));
		//Ce += Wn *Ecurr*Ecurr*exp(-Beta*Ecurr);	
		//cout<<"Obs1: "<<Obs1<<" Obs2: "<<Obs2<<endl;		


	}
	Obs1 = Obs1/Zest;
	Obs2 = Obs2/Zest;
//	Obs3 = Obs3/Zest;
//	Obs4 = Obs4/Zest;
//	Obs5 = Obs5/Zest;
//	Obs6 = Obs6/Zest;

	H = ((1.0/Zest)*He);
	
	C = -( (1.0/(Zest*Zest))*(1.0/(kb*Temp*Temp))*(He*He) ) + ( (1.0/(Zest*kb*Temp*Temp))*Ce ) ;

	long double cvn = C/N;
	unsigned long long na = (unsigned long long)(N);
	unsigned long long nafac = Factorial(na);
	//long double lcnzp = log(1.0/( nafac*pow( h*h*Beta/(2.0*3.14*m), 3.0*N/2.0) ));
	//long double lcnzp = log(cnzp);
	long double G = -kb*Temp*( log(Zest) + S );
//	cout<<" G: "<<G<<" lcnzp: "<<lcnzp<<" nafac: "<<nafac<<" na: "<<na<< " nafac2 "<<Factorial(na)<<endl;
	long double mu = G/N;
	//cout<<std::scientific << setiosflags(ios::fixed)<< setprecision(15)<<"T: "<<Temp<<" du : "<<Obs4<<" du^2: "<<Obs5<<" du^3:"<<Obs6<<endl;
	//	cout<<std::scientific << setiosflags(ios::fixed)<< setprecision(15)<<Temp<<" "<<Obs4<<" "<<Obs5<<" "<<Obs6<<endl;
//	long double gamma = Obs4-0.5*Beta*( Obs5 - Obs4*Obs4);
//	gamma += (1.0/6.0)*Beta*Beta*( Obs6 - 3.0*Obs5*Obs4+2.0*Obs6*Obs6);
//	gamma/=dA;
//	gamma/=2.0;
	//Output to file
//	cout<<"obs1 "<<Obs1<<" obs3 "<<Obs3<<" obs4 "<<Obs4<<endl;
	outfile<<std::scientific << setiosflags(ios::fixed)<< setprecision(15)<<Temp<<" "<<Zest<<" "<<H<<" "<<C<<" "<<cvn<<" "<<G<<" "<<mu<<" "<<Obs1<<" "<<Obs2<<" "<<endl;
}

outfile.close();
cout<<"Generated output file: "<<outfilename.c_str()<<endl;
cout<<"Program Complete."<<endl;
return EXIT_SUCCESS;

}

// Define Functions

//Implicitly casts the input value to variable 
// type double
double CastDouble(double input){
	return(input);
}
unsigned long long factorial(unsigned long long n)
{
  cout<<"n "<<n<<endl;
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

unsigned long long Factorial(unsigned long long n)
{
	if (n==0 || n==1){
		return 1;
	}
	else{
		unsigned long long factor = 1;
		for (unsigned long long i = 2; i <= n; ++i)
		{
			factor*=i;
		}
		return factor;
	}
}

