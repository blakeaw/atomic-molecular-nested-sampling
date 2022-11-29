#include <stdlib.h>
#include <string>
#include <fstream>
#include <math.h>
#include <iostream>
#include <iomanip>
//#include <quadmath.h>
//#include <ctime>

using namespace std;

// Define Classes

// Prototype Functions
double CastDouble(double input);

// Main Function 
int main(int argc, char *argv[])
{
//Input File Name
//string filename = "Hpoints_NS_Ent_na_atoms_a_ttagt_pp_press_p_rep_rrreprr.dat";
//string filename2 = "zObservables_NS_Ent_na_atoms_a_ttagt_pp_press_p_rep_rrreprr.dat";
////Output File Name
//string outfilename = "zNest_Cv_na_atoms_a_ttagt_pp_press_p_rep_rrreprr.dat";
string filename = "Hpoints_NS_Ent_n1_hoc_f3.dat";
string filename2 = "zObservables_NS_Ent_n1_hoc_f3.dat";
//Output File Name
string outfilename = "zNestOut_n1_hoc_f2.dat";
// Define the the value of K used in simulation
//int Ksample = 700;
// Define the number of atoms in simulation
//long double N = (long double)(a_atoms_a) ;

//double Ksampled = CastDouble(Ksample);
// Calculates the alpha based on Ksample
//double alpha = Ksampled/(Ksampled + 1.0);

long double kb = 1.0; 
long double Beta;
long double Temp;
//long double Press = p_press_p ;
//Press = 1.0;
//Planck reduced - material specific
/* argon - 0.185
neon - 0.563
krypton - 0.104
xenon - 0.065
---Ref: Quasi Harmonic Lattice Dynamics and Molecular Dynamics calculations
for the Lennard-Jones solids
*/
long double hr = 0.185;
// define pi
const long double pi = 3.14159265359;
// Initialize an array to store the Energy values read in from file

//string junk;
ifstream infile;
infile.open(filename.c_str());
// Returns an error if the file does not exist
if (!infile) {
	cerr << "Unable to open file " << filename.c_str() << endl;
}

const int i = 244;
long double Hn[i] = {0};
long double tex;


// Reads the file line by line until end 
//while (!infile.eof()){
//	infile >> tex;
//	Hn[i] = tex;
//	++i;
//	if(i==200000-1){
//		break;
//	}
//	
//}
for (unsigned int j = 0; j < i; ++j)
{
	infile >> tex;
//	cout<<"tex "<<tex<<endl;
	Hn[j] = tex;
}

infile.close();
cout<<"read in "<<i<<" nested values "<<endl;
ifstream infile2;
infile2.open(filename2.c_str());
// Returns an error if the file does not exist
if (!infile2) {
	cerr << "Unable to open file " << filename2.c_str() << endl;
}
const unsigned long ii = i;
long double Obs[ii][4];

//long double tex1, tex2, tex3, tex4, tex5, tex6, tex7;
long double tex1, tex2, tex3, tex4;
//int i = 0;
int kuku=0;
long double Vlow=1000000000000.0;
// Reads the file line by line until end 
for (unsigned int j = 0; j < i; ++j)
{
//	infile2 >> junk >> junk >> tex1 >> tex2 >> tex3 >> tex4 >> tex5 >> tex6 >> tex7;
	infile2 >> tex1 >> tex2 >> tex3 >> tex4;
	Obs[j][0] = tex1; //Currently: Mean Energy
	Obs[j][1] = tex2; //Currently: Mean Particle Count
	Obs[j][2] = tex3; //Currently: Density
	Obs[j][3] = tex4; //Currently: <N^2>
//	Obs[kuku][4] = tex5; //Currently: Surface Orientation Param
//	Obs[kuku][5] = tex6; //Currently: NP-NP Orientation Param
//	Obs[kuku][6] = tex7; //Currently: Average Dimer-Surface Separation
	//cout<<"t1: "<<tex1<<" t2 "<<tex2<<" tex4 "<<tex4<<endl;
	
	
}
infile2.close();
//exit(0);


long double kT ;
long double Ck;
//Set the number of Temperature points.
// should be a multiple of 10
int Tloop = 4000;
long double Ts = 0.05;
long double Te = 2.0;
long double dT = (Te - Ts)/(long double)(Tloop);
ofstream outfile(outfilename.c_str());
outfile << "# kT T Z U C Ck" << endl;

for(int T=0;T<=Tloop;++T){
	
	Temp = Ts + dT * (long double)(T);
	kT = kb*Temp;
	Beta = 1.0/kT;
//	long double tw = sqrt( (Beta*pow(hr, 2))/(2.0*pi));
//	long double tw = 1.0;
	long double tw=sqrt(Beta);
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
	long double lwmax = 200.0;
	long double lw;
	long double Obs1 = 0.0;
	long double Obs2 = 0.0;
	long double Obs3 = 0.0;
	long double Obs4 = 0.0;
//	long double Obs5 = 0.0;
//	long double Obs6 = 0.0;
//	long double Obs7 = 0.0;
	long double Obs1curr;
	long double Obs2curr;
	long double Obs3curr;
	long double Obs4curr;
//	long double Obs5curr;
//	long double Obs6curr;
//	long double Obs7curr;
	
	for(n=1;n<i-1;++n){
		
		if(n==0){
			Hcurr =Hn[n];
		}
		else{
			Hcurr = 0.50*(Hn[n]+Hn[n-1]);
			Obs1curr = 0.50*(Obs[n][0]+Obs[n-1][0]);
			Obs2curr = 0.50*(Obs[n][1]+Obs[n-1][1]);
			Obs3curr = 0.50*(Obs[n][2]+Obs[n-1][2]);
			Obs4curr = 0.50*(Obs[n][3]+Obs[n-1][3]);
			//Hcurr = Obs1curr+Press*Obs2curr;
		}
	//	cout<<"Hcurr "<<Hcurr<<endl;
	//	nd = CastDouble(n);
		//np1 = nd + 1.0;
	//	cout<<"n: "<<n<<" obs1curr: "<<Obs1curr<<" obs2curr: "<<Obs2curr<<endl;
		lw = (-((long double)(n))*log(2.0)) - (Beta*Hcurr) - 3.0*Obs2curr*log(tw);//+N*log(Obs2curr);
	//	lw = (-((long double)(n))*log(2.0)) - (Beta*Hcurr);
	//	lw = (-((long double)(n))*log(2.0)) - (Beta*Hcurr) - 3.0*Obs2curr*log(Beta)/2.0;
		//cout<<"lw "<<lw<<" 3.0*Obs2curr*log(Beta)/2.0 "<<3.0*Obs2curr*log(Beta)/2.0<<endl;	
		while(lw>S+lwmax){
			Zest /= exp(lwmax);
			He /= exp(lwmax);
			Ce /= exp(lwmax);
			Obs1 /= exp(lwmax);
			Obs2 /= exp(lwmax);
			Obs3 /= exp(lwmax);
			Obs4 /= exp(lwmax);
		//	Obs5 /= exp(lwmax);
		//	Obs6 /= exp(lwmax);
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
		Obs3 += Wn*Obs3curr;
		Obs4 += Wn*Obs4curr;
	//	Obs5 += Wn*Obs5curr;
	//	Obs6 += Wn*Obs6curr;
	//	Obs7 += Wn*Obs7curr;
		//Wn = pow(2, -(n+1));
		//Zest +=  Wn* exp(-Beta*Ecurr);
		//Ue += Wn * (Ecurr * exp(-Beta*Ecurr));
		//Ce += Wn *Ecurr*Ecurr*exp(-Beta*Ecurr);	
		//cout<<"Obs1: "<<Obs1<<" Obs2: "<<Obs2<<endl;		


	}
		Obs1 = Obs1/Zest;
	Obs2 = Obs2/Zest;
	Obs3 = Obs3/Zest;
	Obs4 = Obs4/Zest;
	H = ((1.0/Zest)*He);

	C = -( (1.0/(Zest*Zest))*(1.0/(kb*Temp*Temp))*(He*He) ) + ( (1.0/(Zest*kb*Temp*Temp))*Ce ) ;

	long double cvn = C/Obs2;


	long double G = -kb*Temp*( log(Zest) + S);

	long double mu = G/Obs2;
	long double Entropy = (H-G)/(Temp);
	
	long double rho = Obs3;

	long double en = Obs1/Obs2;
	long double hn = H/Obs2;
	long double Sn = Entropy/Obs2;
	long double Ndev = Obs4 - Obs2*Obs2;
	//cout<<"<N^2> : "<<Obs4<<" <N>: "<< Obs2 <<" Ndev: "<<Ndev<<endl;
outfile<<std::scientific << setiosflags(ios::fixed)<< setprecision(15)<<Temp<<" "<<Zest<<" "<<H<<" "<<C<<" "<<Obs1<<" "<<Obs2<< " "<<Obs3<<" "<<G<<" "<<mu<<" "<<hn<<" "<<cvn<<" "<<en<<" "<<rho<<" "<<Entropy<<" "<<Sn<< " "<< Ndev<<endl;

}

outfile.close();

return EXIT_SUCCESS;

}

// Define Functions

//Implicitly casts the input value to variable 
// type double
double CastDouble(double input){
	return(input);
}


