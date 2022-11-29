#include <stdlib.h>
#include <string>
#include <fstream>
#include <math.h>
#include <iostream>
#include <iomanip>
//#include <ctime>

using namespace std;

// Define Classes

// Prototype Functions
double CastDouble(double input);

// Main Function 
int main(int argc, char *argv[])
{
//Input File Name
string filename = "Epoints_int_cube_d10_test_itagi_reprrtagrr.dat";
string filename2 = "zNSobs_int_cube_d10_test_itagi_reprrtagrr.dat";
//Output File Name
string outfilename = "zCv_int_cube_d10_test_itagi_reprrtagrr.dat";
string outfilename2 = "zNestOut_int_cube_d10_test_itagi_reprrtagrr.dat";
long double kb = 1.0;
long double Tstart = 0.0001;
long double Tfinish = 10.0;
unsigned long Tloop = 2000;
long double Beta;
long double Temp;


// Initialize an array to store the Energy values read in from file

//string junk;
ifstream infile;
infile.open(filename.c_str());
// Returns an error if the file does not exist
if (!infile) {
	cerr << "Unable to open file: " << filename.c_str() << endl;
//	cout << "Unable to open file: " << filename.c_str() << endl;
	cout << "Exiting..." << endl;
	exit(0);
}
cout << "Extracting Data from input files...." << endl;
long double En[200000] = {0};
long double tex;
unsigned long  i = 0;
long double junk;
// Reads the file line by line until end 
while (! infile.eof()){
	infile >> junk >> tex;
	En[i] = tex;
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
long double Obs[ii][10];

long double tex1, tex2, tex3, tex4, tex5, tex6, tex7, tex8, tex9, tex10;
//int i = 0;
int kuku=0;
// Reads the file line by line until end 
while (! infile2.eof()){
	infile2 >> junk >> junk >> tex1 >> tex2 >> tex3 >> tex4 >> tex5 >> tex6 >> tex7 >> tex8 >> tex9 >> tex10;
	//cout << tex1 << " "<<tex2<<" "<<tex3<<" "<<tex4<<" "<<tex5<<endl;
	Obs[kuku][0] = tex1; //Currently: Mean Energy
	Obs[kuku][1] = tex2; //Currently: Mean (Energy)^2
	Obs[kuku][2] = tex3; //Currently: Radius of Gyration
	Obs[kuku][3] = tex4; //Currently: Gyration Density
	Obs[kuku][4] = tex5; //Currently: Surface Orientation Param
	Obs[kuku][5] = tex6; //Currently: NP-NP Orientation Param
	Obs[kuku][6] = tex7; //Currently: Average Dimer-Surface Separation
	Obs[kuku][7] = tex8; //Currently: PP contacts
	Obs[kuku][8] = tex9; //Currently: PN contacts
	Obs[kuku][9] = tex10; //Currently: NN contacts
//	cout<<"pp: "<<tex8<<" np: "<<tex9<<" nn: "<<tex10<<endl;
//	cout<<"Obskuku8: "<<Obs[kuku][8]<<endl;
	++kuku;
	if(kuku == i){
		break;
	}
}
infile2.close();
long double kT ;
//long double Ck;
//exit(0);

unsigned int T;
ofstream outfile(outfilename.c_str());
outfile << "# kT T Z U C Ck Cvb" << endl;
ofstream outfile2(outfilename2.c_str());
//outfile2 << "# E   Esq   Cv " << endl;

long double dT = (Tfinish-Tstart)/(long double)(Tloop);
cout << "Beginning the Temperature Loop...." << endl;
for(T=1;T<=Tloop;++T){
	
//	Temp = CastDouble(T)/CastDouble(Tloop);
	Temp = Tstart + dT*(long double)(T); // K
	kT = kb*Temp; // Unitless
	Beta = 1.0/kT;
	// Estimate the partition function
	unsigned long n ;
	//long double np1;
	long double Zest = 0.0 ;
	long double nd;
	long double Ue = 0.0;
	long double U;
	long double C;
	long double Ce = 0.0;
	long double Wn;	
	long double Ecurr;
	long double S = 0.0;
	long double lw;
	long double Obs1 = 0.0;
	long double Obs2 = 0.0;
	long double Obs3 = 0.0;
	long double Obs4 = 0.0;
	long double Obs5 = 0.0;
	long double Obs6 = 0.0;
	long double Obs7 = 0.0;
	long double Obs8 = 0.0;
	long double Obs9 = 0.0;
	long double Obs10 = 0.0;
	long double Obs1curr;
	long double Obs2curr;
	long double Obs3curr;
	long double Obs4curr;
	long double Obs5curr;
	long double Obs6curr;
	long double Obs7curr;
	long double Obs8curr;
	long double Obs9curr;
	long double Obs10curr;

	for(n=1;n<i;++n){

		
		//	Ecurr = 0.5*(En[n]+En[n-1]);
		//	Ecurr = 0.5*(Obs[n][0]+Obs[n-1][0]);
			Ecurr = Obs[n][0];
			Obs1curr = Obs[n][0];
			Obs2curr = Obs[n][1];
			Obs3curr = Obs[n][2];
			Obs4curr = Obs[n][3];
			Obs5curr = Obs[n][4];
			Obs6curr = Obs[n][5];
			Obs7curr = Obs[n][6];
			Obs8curr = Obs[n][7];
			Obs9curr = Obs[n][8];
			Obs10curr = Obs[n][9];
			//Ecurr = 0.5*( En[n]+En[n+1] ) ;
			//Ecurr = Obs1curr;
			
		
		
	//	Wn = pow(2, -(n+1)) * exp(-Beta*Ecurr);
		//cout << "Ecurr: " << Ecurr << " Beta: " << Beta << " kb " << kb << " Temp " << Temp <<  endl;
	lw = (-((long double)(n))*log(2.0)) - (Beta*Ecurr);
//		lw = (((long double)(n))*log(9.0/10.0)) - (Beta*Ecurr);
//		lw = (((long double)(n))*log(0.75)) - (Beta*Ecurr);
//		lw = (((long double)(n))*log(0.60)) - (Beta*Ecurr);
//		lw = (((long double)(n))*log(0.10)) - (Beta*Ecurr);
//		lw = (((long double)(n))*log(0.25)) - (Beta*Ecurr);
		long double lwmax = 50.0;
		while(lw>S+lwmax){
			Zest /= exp(lwmax);
			Ue /= exp(lwmax);
			Ce /= exp(lwmax);
			Obs1 /= exp(lwmax);
			Obs2 /= exp(lwmax);
			Obs3 /= exp(lwmax);
			Obs4 /= exp(lwmax);
			Obs5 /= exp(lwmax);
			Obs6 /= exp(lwmax);
			Obs7 /= exp(lwmax);	
			Obs8 /= exp(lwmax);
			Obs9 /= exp(lwmax);	
			Obs10 /= exp(lwmax);		
			
			S+=lwmax;	
			
		}
		lw += -S;
		
		Wn = exp(lw);	
		//cout << "Wn: " << Wn << endl;
		Zest +=  Wn;
		Ue += Wn*Ecurr;
		Ce += Wn*Ecurr*Ecurr;
		
		Obs1 += Wn*Obs1curr;
		Obs2 += Wn*Obs2curr;
		Obs3 += Wn*Obs3curr;
		Obs4 += Wn*Obs4curr;
		Obs5 += Wn*Obs5curr;
		Obs6 += Wn*Obs6curr;
		Obs7 += Wn*Obs7curr;	
		Obs8 += Wn*Obs8curr;	
		Obs9 += Wn*Obs9curr;	
		Obs10 += Wn*Obs10curr;	
	//	cout<<"obs9curr: "<<Obs9curr<<" Obs9 "<<Obs9<<" Obsn8: "<<Obs[n][8]<<endl;
	//exit(0);
	}

	//cout << "**********************************************************************************" << endl;
	// Internal Energy
	long double UoZ = Ue/Zest;
	long double CoZ = Ce/Zest;

	U = UoZ;
	// Heat Capacity
	
	//C =  (1.5*N*kb) - ( (1.0/(Zest*Zest))*(1.0/(kb*Temp*Temp))*(Ue*Ue) ) + ( (1.0/(Zest*kb*Temp*Temp))*Ce ) ;
	C = -( (1.0/(kb*Temp*Temp))*(UoZ*UoZ) ) + ( (1.0/(kb*Temp*Temp))*CoZ ) ;
	
//	cout<<"Obs9a: "<<Obs9<<endl;
	
	//Observables
	Obs1 = Obs1/Zest;
	Obs2 = Obs2/Zest;
	Obs3 = Obs3/Zest;
	Obs4 = Obs4/Zest;
	Obs5 /=Zest;
	Obs6 /=Zest;
	Obs7 /=Zest;
	Obs8 /=Zest;
	Obs9 /=Zest;
	Obs10 /=Zest;
//	cout<<"Obs9b: "<<Obs9<<endl;
	long double Cvb = Obs2 - (Obs1*Obs1);
	Cvb *= (Beta/Temp);
//outfile2 << std::scientific << setiosflags(ios::fixed) << setprecision(15) << Temp << " " << Obs1 << " " << Obs2 << " " << Cvb <<  endl; 	
	

	//Output to file
	outfile << std::scientific << setiosflags(ios::fixed) << setprecision(15);
	 outfile << Temp << " " << Zest << " " << U << " " << C <<  endl;

	outfile2 << std::scientific << setiosflags(ios::fixed) << setprecision(15);
	outfile2 << Temp << " " << Obs1 << " " << Obs2 << " " << Cvb << " " << Obs3;
	outfile2 << " " << Obs4 << " " << Obs5 << " " << Obs6 << " " << Obs7 << " ";
	outfile2<<Obs8<<" "<<Obs9<<" "<<Obs10<<endl;
}
cout << "Processing complete!" << endl;
cout << "Check output files:" << endl;
cout << "    " << outfilename.c_str() << endl;
cout << "    " << outfilename2.c_str() << endl;
cout << " " << endl;

outfile.close();
outfile2.close();
return EXIT_SUCCESS;

}

// Define Functions

//Implicitly casts the input value to variable 
// type double
double CastDouble(double input){
	return(input);
}


