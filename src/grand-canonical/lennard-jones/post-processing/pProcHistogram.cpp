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
string filename = "Hpoints_NS_Ent_n17_comr_d7.dat";
string filename2 = "zNhist_NS_Ent_n17_comr_d7.dat";
//Output File Name
//string outfilename = "zCv_int_cube_d8_test_itagi_reprrtagrr.dat";
string outfilenamepre = "zOut_Nhist_NS_d7_";
string outtdname = "zOut_td_Nhist_NS_d7.dat";
ofstream outtd(outtdname.c_str());
long double kb = 1.0;
long double Tstart = 0.4;
long double Tfinish = 2.0;
unsigned long Tloop = 17;
long double binwidth = 1.0;
const unsigned NumElements = 1001 ;
long double Beta;
long double Temp;
long double chempo = 0.05;

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
	//infile >> junk>>tex;;
	infile>>tex;
 //   cout<<"tex "<<tex<<endl;
	En[i] = tex;
	++i;
	if(i==200000-1){
		break;
	}
}
cout<<"i "<<i<<endl;
infile.close();
--i;
--i;
//exit(0);
ifstream infile2;
infile2.open(filename2.c_str());
// Returns an error if the file does not exist
if (!infile2) {
	cerr << "Unable to open file " << filename2.c_str() << endl;
}


unsigned DEPTH = 2;
long double ***Obs3darray;


  // Allocate memory
  Obs3darray = new long double**[i];
  for (unsigned j = 0; j < i; ++j) {
    Obs3darray[j] = new long double*[NumElements];

    for (unsigned k = 0; k < NumElements; ++k)
      Obs3darray[j][k] = new long double[DEPTH];
  }
////Initialize values of zero in array
//for (unsigned int b = 0; b < i; ++b)
//{
//	
//	for (unsigned int l = 0; l < NumElements; ++l)
//	{
//		
//		//obs value
//		Obs3darray[b][l][0]=0.0;
//		// histogram count
//		Obs3darray[b][l][1]=0.0;
//		//cout<<"l: "<<l<<" tex1: "<<tex1<<" tex2: "<<tex2<<endl;
//	}
//	
//}
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
		Obs3darray[b][l][1]=(long double)(tex2);
	//	cout<<"l: "<<l<<" tex1: "<<tex1<<" tex2: "<<tex2<<" obsbl1 "<<Obs3darray[b][l][1]<<endl;
	}
	
}
//exit(0);
infile2.close();

cout<<"Normalizing.."<<endl;
//Normalize
for (unsigned int b = 0; b < i; ++b)
{
	//Get Total
	long double Total=0.0;
	for (unsigned int l = 0; l < NumElements; ++l)
	{
		
		
		// histogram count
		long double val = Obs3darray[b][l][1];
		Total+=val*binwidth;
		//cout<<"Total: "<<Total<<" val "<<val<<" binwidth "<<binwidth<<" obs "<<Obs3darray[b][l][1]<<endl;
	}
	//cout<<"Total: "<<Total<<endl;
	//Divide by Total
	for (unsigned int l = 0; l < NumElements; ++l)
	{
		
		
		// histogram count
		long double val = Obs3darray[b][l][1];
		//new normalized value
		long double nval = val/Total;
		//replace in array
		Obs3darray[b][l][1]=nval;
		//cout<<"nval "<<nval<<endl;
		if(nval!=nval){
			cout<<"b "<<b<<" l "<<l<<" nval: "<<nval<<" Total: "<<Total<<" val "<<val<<endl;
			exit(0);
		}
	}
}

//exit(0);

long double kT ;
//long double Ck;


unsigned int T;
//ofstream outfile(outfilename.c_str());
//outfile << "# kT T Z U C Ck Cvb" << endl;
//ofstream outfile2(outfilename2.c_str());
//outfile2 << "# E   Esq   Cv " << endl;

long double dT = (Tfinish-Tstart)/(long double)(Tloop-1);
cout << "Beginning the Temperature Loop...." << endl;
for(T=0;T<Tloop;++T){
	
//	Temp = CastDouble(T)/CastDouble(Tloop);
	Temp = Tstart + dT*(long double)(T); // K
	cout<<"Temp: "<<Temp<<endl;
	kT = kb*Temp; // Unitless
	Beta = 1.0/kT;
	//cout<<"Temp "<<Temp<<" kb "<<kb<<" kT "<<kT<<" Beta "<<Beta<<endl;
	// Estimate the partition function
	unsigned long n ;
	//long double np1;
	
	long double Zest=0.0;
	long double Ce = 0.0;
	long double Wn;	
	long double Ecurr;
	long double S = 0.0;
	long double lwmax = 100.0;
	long double lw;
	long double Zprof[NumElements] = {0.0};

	for(n=0;n<i;++n){

		

			//Ecurr = En[n];
		if(n==0){
			Ecurr =En[n];
		}
		else{
			Ecurr = 0.50*(En[n]+En[n-1]);
		
		}
	
		lw = (-((long double)(n))*log(2.0)) - (Beta*Ecurr);
	//	cout<<"lw "<<lw<<" Ecurr" <<Ecurr<<" Beta "<<Beta<<endl;
		while(lw>S+lwmax){
			Zest /= exp(lwmax);
		
			for (unsigned int p = 0; p < NumElements; ++p)
			{
				Zprof[p]/=exp(lwmax);
			}
			
			S+=lwmax;	
			
		}

		lw += -S;
		
		Wn = exp(lw);	
		//cout << "Wn: " << Wn << endl;
		Zest +=  Wn;
		//cout<<"zest "<<Zest<<" Wn "<<Wn<<" lw "<<lw<<endl;
		if (Zest!=Zest){
			exit(1);
		}
		for (unsigned int p = 0; p < NumElements; ++p)
			{
				//cout<<"Zprofp "<<Zprof[p]<<" p "<<p<<endl;	
				if (n==0){
					Zprof[p]+=Wn*Obs3darray[n][p][1];
				}				
				else{
					Zprof[p]+=Wn*0.5*(Obs3darray[n][p][1]+Obs3darray[n-1][p][1]);
				}
				//cout<<"reweighted Zprofp "<<Zprof[p]<<" p "<<p<<" Wn "<<Wn<<" val "<<Obs3darray[n][p][1]<<endl;		
				
			}
			
		
	
	//exit(0);
	} //end loop nested iterations
	
	//cout << "**********************************************************************************" << endl;
//	long double Aprof[NumElements] = {0.0};
//		for (unsigned int p = 0; p < NumElements; ++p)
//			{
//						
//				Aprof[p]=-Temp*log(Zprof[p]);
//			}
		for (unsigned int p = 0; p < NumElements; ++p)
		{
						
				Zprof[p]/=Zest;
			//	cout<<"zest "<<Zest<<" Zprofp "<<Zprof[p]<<endl;
		}
		string cTemp = dtos(Temp);
		string ofname = outfilenamepre+"T"+cTemp+".dat";
		
		ofstream outfile(ofname.c_str());
//		long double Total = 0.0;
//		for (unsigned int ne = 0; ne < NumElements; ++ne)
//		{
//					Total+=Zprof[ne]*binwidth;
//		}
//	//	Total*=binwidth;
//		for (unsigned int ne = 0; ne < NumElements; ++ne)
//		{
//			Zprof[ne]/=Total;
//		}
		for (unsigned int p = 0; p < NumElements; ++p)
			{
				long double cObs = Obs3darray[1][p][0];	
			//	long double cA = Aprof[p];		
				outfile<<cObs<<" "<<Zprof[p]<<endl;
				outtd<<chempo<<" "<<Temp<<" "<<" "<<cObs<<" "<<Zprof[p]<<endl;
			}
		//outtd<<endl;
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

