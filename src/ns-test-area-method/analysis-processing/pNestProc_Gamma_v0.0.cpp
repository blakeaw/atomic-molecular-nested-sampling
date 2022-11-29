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

/*
Note: Files should have the same number of
temperature points across the same range.
i.e. all temperature values should match.
*/

//Input File Name
string filename = "a0_clean.dat";
//number of columns
unsigned ncol1 = 10;
//column number that free energy is in - starts at 1
unsigned lcol1 = 6;
//column number that temperature is in - starts at 1
unsigned tcol1 = 1;
string filename2 = "b0_clean.dat";
unsigned ncol2 = 10;
unsigned lcol2 = 6;
//column number that temperature is in - starts at 1
unsigned tcol2 = 1;
//Output File Name
string outfilename = "zNest_Gamma256_a0_b0.dat";

// Define the number of atoms in simulation
long double N = 256.0 ;


//Define dimension of slab
long double sx1, sy1, sx2, sy2;
sx1 = 6.7, sy1 = 6.7;
sx2 = 6.7*1.01, sy2 = 6.7*1.01; 

long double a1 = sx1*sy1*2.0;
long double a2 = sx2*sy2*2.0;

// Initialize an array to store the Energy values read in from file

//string junk;
ifstream infile;
infile.open(filename.c_str());
// Returns an error if the file does not exist
if (!infile) {
	cerr << "Unable to open file " << filename.c_str() << endl;
	exit(0);
}

long double G1[200000][2];
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
			++i;
		//	cout<<"tex: "<<tex<<endl;
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
const unsigned ii = i;
long double G2[ii][2];

//long double tex1, tex2, tex3, tex4, tex5, tex6, tex7;


cout<<"Will read in file: "<<filename2<<endl;
// Reads the file line by line until end 
for (unsigned int j = 0; j < i; ++j)
{
	
	for (unsigned int jj = 0; jj < ncol2; ++jj)
	{
		infile2 >> tex;
		if (jj+1 == tcol2){
			G2[j][0] = tex;
		}
		if (jj+1 == lcol2){
			G2[j][1] = tex;
			
		}
		
	}
			
}

infile2.close();
//exit(0);


ofstream outfile(outfilename.c_str());
outfile << "# T G1 A1 G2 A2 dG dA Gamma" << endl;
cout<<"Opening file: "<<outfilename<<" for writing.."<<endl;

cout<<"Beggining the Loop.."<<endl;
for (unsigned int bb = 0;bb < i-1; ++bb)
{
//	cout<<"T: "<<G1[bb][0]<<" g1: "<<G1[bb][1]<<" g2: "<<G2[bb][1]<<endl;
	long double T = G1[bb][0];
	long double g1 = G1[bb][1];
	long double g2 = G2[bb][1];
	long double dG = g2 - g1;
	long double dA = a2 - a1;
	long double gamma = dG/dA;

	outfile<<T<<" "<<g1<<" "<<a1<<" "<<g2<<" "<<a2<<" "<<dG<<" "<<dA<<" "<<gamma<<endl;
	cout<<T<<" "<<g1<<" "<<a1<<" "<<g2<<" "<<a2<<" "<<dG<<" "<<dA<<" "<<gamma<<endl;
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


