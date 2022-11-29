#include <stdlib.h>
#include <string>
#include <fstream>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cmath>
//#include <ctime>

using namespace std;

// Define Classes

// Prototype Functions
double CastDouble(double input);
unsigned long long factorial(unsigned long long n);
// Main Function 
int main(int argc, char *argv[])
{

/*
Note: Files should have the same number of
temperature points across the same range.
i.e. all temperature values should match.
*/

//Input File Name -- binary system
string filename = "aNest_ab2_G.dat";
//number of columns
unsigned ncol1 = 4;
//column number that free energy is in - starts at 1
unsigned lcol1 = 2;
//column number that temperature is in - starts at 1
unsigned tcol1 = 1;
//column number that error is in - starts at 1
unsigned acol1 = 4;
string filename3 = "aNest_ab2_A.dat";
//number of columns
unsigned ncol3 = 4;
//column number that area is in - starts at 1
unsigned lcol3 = 2;
//column number that temperature is in - starts at 1
unsigned tcol3 = 1;
//column number that error is in - starts at 1
unsigned acol3 = 4;
string filename2 = "aNest_a0_G.dat";
unsigned ncol2 = 4;
unsigned lcol2 = 2;
//column number that temperature is in - starts at 1
unsigned tcol2 = 1;
//column number that error is in - starts at 1
unsigned acol2 = 4;
//Output File Name
string outfilename = "zNest_Gamma_ab2_a0.dat";

// Define the number of atoms in simulation
long double N = 100.0 ;
const long double pi = 3.14159265359;

//Define dimension of slab
long double sx1, sy1, sx2, sy2;
sx1 = 5.0, sy1 = 5.0;
sx2 = 6.7*1.01, sy2 = 6.7*1.01; 
long double P=0.1;
long double a1 = sx1*sy1*2.0;
long double a2 = sx2*sy2*2.0;
long double rad = 2.0;
long double Vb = 90000.0;
long double area = a1;
long double Va = Vb;
long double Vab = 90000.0;
long double Na = 101.0;
long double Nb = 101.0;
long double Nab = 201.0;
//Nb=3.2;
long double kb = 1.0;
long double h =  0.185;
// Initialize an array to store the Energy values read in from file

//string junk;
ifstream infile;
infile.open(filename.c_str());
// Returns an error if the file does not exist
if (!infile) {
	cerr << "Unable to open file " << filename.c_str() << endl;
	exit(0);
}

long double G1[200000][3];
long double tex;
unsigned i = 0;
cout<<"Will read in file: "<<filename<<endl;
// Reads the file line by line until end 
while (! infile.eof()){
	for (unsigned int ii = 0; ii < ncol1; ++ii)
	{
		infile >> tex;
		//cout<<"tex: "<<tex<<endl;
		if (ii == tcol1-1 ){
			G1[i][0] = tex;
			//cout<<"tex: "<<tex<<" ";
		}
		if (ii == lcol1-1){
			G1[i][1] = tex;
			
		//	cout<<"tex: "<<tex<<endl;
		}
		if (ii == acol1-1){
			G1[i][2] = tex;
			++i;
			//cout<<"tex: "<<tex<<endl;
		}
		if(i==200000-1){
		break;
	}
	}
		
}
infile.close();
//exit(0);

ifstream infile2;
infile2.open(filename2.c_str());
// Returns an error if the file does not exist
if (!infile2) {
	cerr << "Unable to open file " << filename2.c_str() << endl;
	exit(0);
}
const unsigned iii = i;

long double G2[iii][3];

cout<<"Will read in file: "<<filename2<<endl;
// Reads the file line by line until end 
for (unsigned int j = 0; j < iii; ++j)
{
	for (unsigned int ii = 0; ii < ncol1; ++ii)
	{
		infile2 >> tex;
		//cout<<"tex: "<<tex<<endl;
		if (ii == tcol2-1 ){
			G2[j][0] = tex;
			//cout<<"tex: "<<tex<<" ";
		}
		if (ii == lcol2-1){
			G2[j][1] = tex;
			
		//	cout<<"tex: "<<tex<<endl;
		}
		if (ii == acol2-1){
			G2[j][2] = tex;
	
			//cout<<"tex: "<<tex<<endl;
		}
	
	}
		
}
infile2.close();
//exit(0);


//string junk;
ifstream infile3;
infile3.open(filename3.c_str());
// Returns an error if the file does not exist
if (!infile3) {
	cerr << "Unable to open file " << filename3.c_str() << endl;
	exit(0);
}

long double G3[iii][3];


cout<<"Will read in file: "<<filename3<<endl;
// Reads the file line by line until end 
for (unsigned int j = 0; j < iii; ++j)
{
	
	for (unsigned int ii = 0; ii < ncol1; ++ii)
	{
		infile3 >> tex;
		//cout<<"tex: "<<tex<<endl;
		if (ii == tcol3-1 ){
			G3[j][0] = tex;
			//cout<<"tex: "<<tex<<" ";
		}
		if (ii == lcol3-1){
			G3[j][1] = tex;
			
		//	cout<<"tex: "<<tex<<endl;
		}
		if (ii == acol3-1){
			G3[j][2] = tex;
			
			//cout<<"tex: "<<tex<<endl;
		}
		
	}
		
}
infile3.close();



ofstream outfile(outfilename.c_str());
outfile << "# T G1 A1 G2 A2 dG dA Gamma" << endl;
cout<<"Opening file: "<<outfilename<<" for writing.."<<endl;
//Nab-=1.0;
cout<<"Beggining the Loop.."<<endl;
for (unsigned int bb = 0;bb < i-1; ++bb)
{
//	cout<<"T: "<<G1[bb][0]<<" g1: "<<G1[bb][1]<<" g2: "<<G2[bb][1]<<endl;
	long double T = G1[bb][0];
	long double g1 = G1[bb][1];
	long double g2 = G2[bb][1];
	long double dG =  g1-2.0*g2;
	long double dA = G3[bb][1];
	long double eg1 = G1[bb][2];
	long double eg2 = G2[bb][2];
	long double ear = G3[bb][2];
	//cout<<"dA "<<dA<<endl;
	long double gamma = dG/dA;
	long double gab = g1;
	long double ga = g2;
	long double gb = ga;
	long double dg = gab - (ga+gb);
	long double dgv = -T*( Nab*log(Vab) - Na*log(Va) - Nb*log(Vb) );
	long double lh = 3.0*Nab*log(h);	
	long double lnb = 1.5*log(1.0/(T*2.0*pi));
	long double dg2 = dg - lh - lnb;
	long double dgv2 = -T*( log( pow(Vab, Nab)/(pow(Va, Na)*pow(Vb, Nb)) ) );
	long double dgVo = T*log(T/P) - 2.0*T*log(T/P);
//	long double dgv = -T*( Nab*log(Vab) - Na*log(Va) - Nb*log(area) );
	//long double dgv = -T*(Na*log(Vab)+Nb*log(Vab)+log(log(Vab))-Na*log(Va)-Nb*log(Vb)-log(log(Vb)));
	//long double dgb = dg;
	long double dg3 = dg-dgv;
	dg += dgv;
	
	gamma = (dg+dgVo)/(2.0*dA);
	
	long double error = abs(gamma)*sqrt( (eg1*eg1+2.0*eg2*eg2)/pow(g1-2.0*g2, 2) + pow(ear/dA, 2)  );
	//gamma/=T;
	//gamma*=58.8;
	//gamma/=10.0;
	//gamma/=2.0;
	//gamma*=-1.0;
	//long double 
	//
	//long double gamma2 = dG - 
	outfile<<T<<" "<<gab<<" "<<dg<<" "<<ga<<" "<<area<<" "<<dG<<" "<<dA<<" "<<gamma<<" "<<error<<endl;
	//cout<<T<<" "<<g1<<" "<<a1<<" "<<g2<<" "<<a2<<" "<<dG/dA<<" "<<dA<<" "<<gamma<<" "<<gamma/2.0<<endl;
	//cout<<" dG "<<dG<<" dg "<<dg<<" va "<<Na*log(Va)<<" Vb "<<Nb*log(Vb)<<" Vab "<<Nab*log(Vab)<<" dgv "<<dgv<<" dgv2 "<<dgv2<<endl;
	//cout<<"lh "<<T<<" lnb "<<lnb<< " dg2 "<<dg2<<" dg3 "<<dg3<<" dg3/dA "<<dg3/dA<<" dg3/(2dA) "<<dg3/(2.0*dA)<<" dgVo "<<dgVo<<" (dg+dgVo)/2dA "<<(dg+dgVo)/(2.0*dA)<<endl;
}

outfile.close();
cout<<"Complete!"<<endl;
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
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
