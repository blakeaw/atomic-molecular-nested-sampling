/*
Program to read in temperature vs Cv values from a file
and average the Cv values
**NOTE: Each set should be catted to a single file and 
the number of lines for an individual set (before concatenation)
must be specified correctly(variable -> numpts). If there
are any comments in the Cv outputs they need to be removed before/during
catting.
**Also note that this program assumes that the same
number of data points was collected for each set
and that the temperature values are the same for each 
set
*/

#include <stdlib.h>
#include <string>
#include <fstream>
#include <math.h>
#include <iostream>
#include <iomanip>
//#include <ctime>

using namespace std;

// Define Classes

// Class to collect running statistics
class RunningStats {

	public:
		RunningStats(): n(0), Mn_old(0.0), Mn_new(0.0), Sn_old(0.0), Sn_new(0.0) {}
		void Push(long double val){
			++n;
			if(n == 1){
				Mn_old = val;
				Sn_old = 0.0;
			}
			else{
				Mn_new = Mn_old + (val - Mn_old)/(double(n));
				Sn_new = Sn_old + (val - Mn_old)*(val-Mn_new);
				
				Mn_old = Mn_new;
				Sn_old = Sn_new;
				
			}
		}

	 long double Mean() const {
            return (n > 0) ? Mn_new : 0.0;
        }

        long double Variance() const {
            return ( (n > 1) ? Sn_new/(n - 1) : 0.0 );
        }

        long double StandDev() const
        {
            return sqrt( Variance() );
        }

	void Reset(void){
		n=0;
		return;
	}

    private:
        unsigned int n;
        long double Mn_old, Mn_new, Sn_old, Sn_new;

};
// Prototype Functions


// Main Function
int main(int argc, char *argv[])
{

//Parameters to set
//Number of sets to avg 
const unsigned numsets =3;
//Number of data points in each set
const unsigned numpts = 2000;
//Number of columns in the file/set
const unsigned numcol = 6;
//Column number that temp is stored. Starts at 1
const unsigned tcol = 1;
//Column number that Cv is stored. Starts at 1
const unsigned ccol = 4;
//Input file
string infilename = "i_file_i";
//Output file
string outfilename = "aNest_Cv_avg_rscalek_val_k.dat";

//to store temperature values
long double Temp[numpts];
//to store heat capacities
long double Cv[numpts][numsets];

//Read in the file
cout << endl;
ifstream infile;
infile.open(infilename.c_str());
if (!infile) {
	cerr << "!! Unable to open file " << infilename.c_str() << endl;
	exit(0);
}
cout << "# Opening: " << infilename.c_str() << " as the input...." << endl;
//Get the values from the file and store them
//while (! infile.eof()){
cout << "# Extracting Data from input file..." << endl;
for(unsigned i =0;i<numsets;++i){
	for(unsigned k=0;k<numpts;++k){
		long double input=0.0;
		for(unsigned j=0;j<numcol;++j){
			infile >> input;
			//cout << "input: " << input << endl;
			if(i == 0 && j==tcol-1){
				//cout << "Temp " << input << " read in!" << endl;
				Temp[k]=input;
			}
			if(j == ccol-1){
				//cout << "Heat Capacity " << input << " read in!" << endl;
				Cv[k][i]=input;
			}
		}
	}

}
//}
infile.close();
//Do the averaging of the heat capacities
RunningStats Cvstat;
ofstream outfile(outfilename.c_str());
outfile  << setiosflags(ios::fixed) << setprecision(15);
cout << "# Averaging the Cv values..." << endl;
for(unsigned i=0;i<numpts;++i){
	for(unsigned j=0;j<numsets;++j){
		Cvstat.Push(Cv[i][j]);
	}
	outfile << Temp[i] << " " << Cvstat.Mean() << " " << Cvstat.StandDev() << endl;
	Cvstat.Reset();

}
outfile.close();
cout << "Program Complete. Check the output file: " << outfilename.c_str() << endl;
cout << " " << endl;


return EXIT_SUCCESS;

}

// Define Functions
