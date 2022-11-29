#include <stdlib.h>
#include <string>
#include <fstream>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <cmath>

using namespace std;

// Define Classes

// Prototype Functions


// Main Function
int main(int argc, char *argv[])
{

string infilename = "./MCresults/LJ17_p0.001.dat";
string infilename2 = "./lj17/comr/aNest_n17_comr_p0.001_cpn.dat";

int ncol = 17;
int tcol = 1;
int dcol = 16;
int ecol = 17;
int npoints = 17;

int ncol2 = 4;
int tcol2 = 1;
int	dcol2 = 2;
int ecol2 = 3;
int npoints2 = 2001;



ifstream infile;
infile.open(infilename.c_str());
// Returns an error if the file does not exist
if (!infile) {
	cerr << "Unable to open file " << infilename.c_str() << endl;
	exit(0);
}
ifstream infile2;
infile2.open(infilename2.c_str());
// Returns an error if the file does not exist
if (!infile2) {
	cerr << "Unable to open file " << infilename2.c_str() << endl;
	exit(0);
}

long double dat[npoints][3];
long double dat2[npoints2][3];

long double input;
cout<<" Reading in file one: "<<infilename.c_str()<<".."<<endl;
for (unsigned int i = 0; i < npoints; ++i)
	{
		for (unsigned int j = 0; j < ncol; ++j)
		{
			infile>>input;
			if (j==tcol-1){
				dat[i][0]=input;
			}
			else if (j==dcol-1){
				dat[i][1]=input;
			}
			else if (j==ecol-1){
				dat[i][2]=input;
			}
		}
	}
infile.close();
cout<<" Reading in file two: "<<infilename2.c_str()<<".."<<endl;

	for (unsigned int i = 0; i < npoints2; ++i)
	{
		for (unsigned int j = 0; j < ncol2; ++j)
		{
			infile2>>input;
			//cout<<"input: "<<input<<endl;
			if (j==tcol2-1){
				//cout<<"T: "<<input<<endl;
				dat2[i][0]=input;
			}
			else if (j==dcol2-1){
				//cout<<"val: "<<input<<endl;
				dat2[i][1]=input;
			}
			else if(j==ecol2-1){
				dat2[i][2]=input;
			}
		}
	}

infile2.close();



cout<<" Going to compute ratios.."<<endl;
long double pdmaxprev = 0.0;
long double pdmaxprev2 = 0.0;
long double pdsum = 0.0;
long double perrsum=0.0;
long double pratsum=0.0;
long double tme = 0.0;
long double tme2=0.0;
for (unsigned int i = 0; i < npoints; ++i)
{
	
	long double Tcurr = dat[i][0];
	long double Valcurr = dat[i][1];
	long double Errcurr = dat[i][2];
	long double Tcurr2;
	long double Valcurr2;
	long double ecurr2;
	//find closest value in data set 2
	long double diffprev = 1000.0;
	
	for (unsigned int j = 0; j < npoints2; ++j)
	{
		long double tcurr = dat2[j][0];
		//cout<<"tcurr: "<<tcurr<<endl;
		long double diff = abs(tcurr - Tcurr);
		if (diff<diffprev){
			Tcurr2=tcurr;
			Valcurr2=dat2[j][1];
			ecurr2=2.0*dat2[j][2];
			diffprev=diff;	
		//	cout << setiosflags(ios::fixed) << setprecision(15)<<"Tcurr: "<<Tcurr<<" tcurr: "<<tcurr<<" diff: "<<diff<<endl;		
		}
		
	}
	//cout<<"Tcurr: "<<Tcurr<<" Tcurr2: "<<Tcurr2<<endl;
	//cout<<"Val1: "<<Valcurr<<" Val2: "<<Valcurr2<<endl;
	//compare
	long double diff = abs(Valcurr2-Valcurr);
	long double pdiff = (diff/Valcurr)*100.0;
	long double perr = (Errcurr/ecurr2);
	long double prat = pdiff/perr;
	//cout<<" Tcurr: "<<Tcurr<<" Tcurr2: "<<Tcurr2<<" diff: "<<diff<<" pdiff: "<<pdiff<<" perr: "<<perr<<" prat: "<<prat<<endl;
	if (Errcurr>pdmaxprev){
		pdmaxprev=Errcurr;
		tme=Tcurr;
	}
	if (ecurr2>pdmaxprev2){
		pdmaxprev2=ecurr2;
		tme2=Tcurr2;
	}
	pdsum+=pdiff;
	pratsum+=prat;
	perrsum+=perr;
}
cout<<"maximum PT error : "<<pdmaxprev<<" at T="<<tme<<endl;
cout<<"maximum NS error : "<<pdmaxprev2<<" at T="<<tme2<<endl;
long double maxrat = pdmaxprev/pdmaxprev2;
long double perravg = perrsum/(long double)(npoints);

cout<<"maximum error ratio: "<<maxrat<<endl;
cout<<"average error ratio (PT/NS): "<<perravg<<endl;

return EXIT_SUCCESS;

}

// Define Functions
