#include <stdlib.h>
#include <string>
#include <fstream>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>

//#include <ctime>

using namespace std;

// Define Classes
//class Array3d{

//	public:
	
//	Array3d():{




//}
// Prototype Functions
double CastDouble(double input);
string dtos(long double Number);

// Main Function 
int main(int argc, char *argv[])
{
//Input File Name
string filename = "zNSobs_int_cube_d8_test_itagi_reprrtagrr.dat";
string filename2 = "zNS_dist_otago_int_cube_d8_test_itagi_reprrtagrr.dat";
//Output File Name
//string outfilename = "zCv_int_cube_d8_test_itagi_reprrtagrr.dat";
string outfilenamepre = "zNS_freeorder_otago_int_cube_d8_test_itagi_reprrtagrr_";
long double kb = 1.0;
long double Tstart = 0.5;
long double Tfinish = 5.0;
unsigned long Tloop = 5;
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
	infile >> junk>>junk>>tex>>junk>>junk>>junk>>junk>>junk>>junk>>junk>>junk>>junk;
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

const unsigned NumElements = n_el_n ;
unsigned DEPTH = 2;
long double ***Obs3darray;


  // Allocate memory
  Obs3darray = new long double**[i];
  for (unsigned j = 0; j < i; ++j) {
    Obs3darray[j] = new long double*[NumElements];

    for (unsigned k = 0; k < NumElements; ++k)
      Obs3darray[j][k] = new long double[DEPTH];
  }
//Read Histogram Data in from file
for (unsigned int b = 0; b < i; ++b)
{
	long double tex1, tex2;
	for (unsigned int l = 0; l < NumElements; ++l)
	{
		infile2 >> tex1 >> tex2;
		//obs value
		Obs3darray[b][l][0]=tex1;
		// histogram count
		Obs3darray[b][l][1]=tex2;
	//	cout<<"l: "<<l<<" tex1: "<<tex1<<" tex2: "<<tex2<<endl;
	}
	
}

infile2.close();


//Normalize
for (unsigned int b = 0; b < i; ++b)
{
	//Get Total
	long double Total=0.0;
	for (unsigned int l = 0; l < NumElements; ++l)
	{
		
		
		// histogram count
		long double val = Obs3darray[b][l][1];
		Total+=val;
	}
	//Divide by Total
	for (unsigned int l = 0; l < NumElements; ++l)
	{
		
		
		// histogram count
		long double val = Obs3darray[b][l][1];
		//new normalized value
		long double nval = val/Total;
		//replace in array
		Obs3darray[b][l][1]=nval;
	}
}



long double kT ;
//long double Ck;


unsigned int T;
//ofstream outfile(outfilename.c_str());
//outfile << "# kT T Z U C Ck Cvb" << endl;
//ofstream outfile2(outfilename2.c_str());
//outfile2 << "# E   Esq   Cv " << endl;

long double dT = (Tfinish-Tstart)/(long double)(Tloop);
cout << "Beginning the Temperature Loop...." << endl;
for(T=0;T<=Tloop;++T){
	
//	Temp = CastDouble(T)/CastDouble(Tloop);
	Temp = Tstart + dT*(long double)(T); // K
	cout<<"Temp: "<<Temp<<endl;
	kT = kb*Temp; // Unitless
	Beta = 1.0/kT;
	// Estimate the partition function
	unsigned long n ;
	//long double np1;
	
	
	long double Ce = 0.0;
	long double Wn;	
	long double Ecurr;

	long double lw;
	long double Zprof[NumElements] = {0.0};

	for(n=1;n<i;++n){

		

			Ecurr = En[n];
			
		
		
	
	lw = (-((long double)(n))*log(2.0)) - (Beta*Ecurr);

		
		Wn = exp(lw);	
		//cout << "Wn: " << Wn << endl;
		//Zest +=  Wn;
		
		for (unsigned int p = 0; p < NumElements; ++p)
			{
						
				Zprof[p]+=Wn*abs(Obs3darray[n][p][0]);
			}
			
		
	
	//exit(0);
	} //end loop nested iterations
	
	//cout << "**********************************************************************************" << endl;
	long double Aprof[NumElements] = {0.0};
		for (unsigned int p = 0; p < NumElements; ++p)
			{
						
				Aprof[p]=-Temp*log(Zprof[p]);
			}

		string cTemp = dtos(Temp);
		string ofname = outfilenamepre+"T"+cTemp+".dat";
		
		ofstream outfile(ofname.c_str());
		for (unsigned int p = 0; p < NumElements; ++p)
			{
				long double cObs = Obs3darray[1][p][0];	
				long double cA = Aprof[p];		
				outfile<<cObs<<" "<<cA<<endl;
			}
		outfile.close();
		
} // end loop Temp





cout << "Processing complete!" << endl;
//cout << "Check output files:" << endl;
//cout << "    " << outfilename.c_str() << endl;
//cout << "    " << outfilename2.c_str() << endl;
//cout << " " << endl;

//outfile.close();
//outfile2.close();
return EXIT_SUCCESS;

}

// Define Functions

//Implicitly casts the input value to variable 
// type double
double CastDouble(double input){
	return(input);
}

string dtos(long double Number ){
     ostringstream ss;
     ss << Number;
     return ss.str();
 }

