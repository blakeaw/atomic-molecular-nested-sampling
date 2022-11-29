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
void quickSort(long double* arr, unsigned long left, unsigned long right, unsigned long indexsort[]);
void POVRayOut(string* mat, const unsigned long& ccurr, const unsigned long Eindsort[], const unsigned& na, const unsigned& nb, const unsigned long long& cls, string inname, const long double& bz);

// Main Function
int main(int argc, char *argv[])
{

const unsigned maxConf = 5000;
//Format tag -- should be xyz
string format = ".xyz";
/*Name of file -- without file type i.e. 
  coord.xyz -> coord    */ 
string ifname = "ConfigsDimer_int_cube_d8_test_24_rep2";

string ofname = ifname+"_sorted"+format;
//number of comment strings before energy on comment line;
unsigned short numb = 6;
//number of atoms/beads in each configuration
unsigned numparticles = 2* 8;
//box size in z direction--used for POVRayOut function
long double boxz=(long double)(numparticles)*3.0+1.5;

numb+=1;
unsigned numa = numparticles*4;
ifname+=format;
ifstream infilea;
infilea.open(ifname.c_str());

if (!infilea) {
	cerr << "Unable to open input file: " << ifname.c_str() << endl;
//	cout << "Unable to open file: " << filename.c_str() << endl;
	cout << "Exiting..." << endl;
	exit(0);
}

string junk;
unsigned long long numelements =0;
const unsigned maxElements = (numb+numa)*maxConf;
cout << "Determinging number of elements and configurations..."<<endl;
while(! infilea.eof()){
	infilea >> junk;
//	cout << "junk: " <<junk<<endl;
	++numelements;
	if(numelements>=maxElements){
		cout <<maxConf<<" configuration limit reached. Stopping at current point. Increase maximum number of configurations to read more"<<endl;
		break;
	}
}
infilea.close();
unsigned long numconfig = numelements/(numb+numa+1);
numconfig -=1;
numelements -= (numb+numa+1);
cout << "Will process " << numconfig<<" configurations.."<<endl;


string* all;
all = new string[numelements];
long double Energies[numconfig];


ifstream infile;

infile.open(ifname.c_str());
// Returns an error if the file does not exist
if (!infile) {
	cerr << "Unable to open input file: " << ifname.c_str() << endl;
//	cout << "Unable to open file: " << filename.c_str() << endl;
	cout << "Exiting..." << endl;
	exit(0);
}
cout << std::scientific << setiosflags(ios::fixed) << setprecision(15);

ofstream outfile(ofname.c_str());
string inputs;
int tracker = 0;
cout <<"Reading in input..."<<endl;
for (unsigned long long i = 0; i < numelements; ++i)
{
	infile >> inputs;	
	all[i] = inputs;
	if(tracker==numb){
		unsigned index = i/(numa+numb+1);
		Energies[index]=atof(inputs.c_str());
		//cout << "inputs "<<inputs<<" index "<<index<<" Energy "<<Energies[index]<<endl;
	}
	
	++tracker;
	if (tracker%(numa+numb+1)==0)
	{
		tracker=0;
	}
}

infile.close();

unsigned long Eindexsort[numconfig];
for (unsigned int i = 0; i < numconfig; ++i)
{
	Eindexsort[i]=i;
}
//long double * Epoint
cout<<"Sorting...."<<endl;
quickSort(Energies, 0, (numconfig-1), Eindexsort);
cout<<"Finished sorting..."<<endl;
//exit(0);
unsigned long long cels = numa+numb+1;
//cout << "cels "<<cels<<endl;
cout<<"Outputting sorted values to file; Energies high to low..."<<endl;
for (unsigned long long i = numconfig-1; i >= 0; --i)
{
	unsigned long long cind = Eindexsort[i];
	unsigned long long sel = (cels)*cind;
//	cout<<"cind "<<cind<<" sel "<<sel<<endl;
	//Number of atoms/beads
	outfile<<all[sel].c_str()<<endl;
	//Comment line
	for (unsigned  j = 0; j < numb; ++j)
	{
		unsigned long long jind=sel+1+j;
		outfile<<all[jind].c_str()<<" ";
	}
	outfile<<endl;
	//Type and Coordinates
	for (unsigned j = 1; j <= numa; ++j)
	{
		
			unsigned long long nind=sel+numb+j;
		
			outfile<<all[nind].c_str()<<" ";
	//		cout<<"nind "<<nind<<" all "<<all[nind]<<endl;
			if(j%4==0){
				outfile<<endl;
			}
	}
	if (i==0){
		
		POVRayOut(all, i, Eindexsort, numa, numb, cels, ifname, boxz);
		outfile.close();
		cout<<" "<<endl;
		cout<<"Config Sorting Complete!"<<endl;
		exit(1);
	}
}

outfile.close();
//for (unsigned int i = 0; i < numconfig; ++i)
//{
	//cout <<"i "<<i<<" Eindex[i] "<<Eindexsort[i]<<endl;
//}
return EXIT_SUCCESS;

}

// Define Functions

void quickSort(long double* arr, unsigned long left, unsigned long right, unsigned long indexsort[]) {

      unsigned long i = left, j = right;

     long double tmp;
	  unsigned long temp;
      long double pivot = arr[(left + right) / 2];

 

      /* partition */

      while (i <= j) {

            while (arr[i] < pivot){

                  i++;
			}
            while (arr[j] > pivot){

                  j--;
			}
            if (i <= j) {

                  tmp = arr[i];

                  arr[i] = arr[j];

                  arr[j] = tmp;
				  temp=indexsort[i];
				  indexsort[i]=indexsort[j];
				  indexsort[j]=temp;
				

                  i++;

                  j--;

            }

      };

 

      /* recursion */

      if (left < j){

            quickSort(arr, left, j, indexsort);
	 }

      if (i < right){

            quickSort(arr, i, right, indexsort);
	}
	//for (unsigned int k = 0; k < right; ++k)
	//{
	//	cout<<"Energy: "<<arr[k]<<" oindex "<<indexsort[k]<<endl;
	//}

	return;
}

void POVRayOut(string* mat, const unsigned long& ccurr, const unsigned long Eindsort[], const unsigned& na, const unsigned& nb, const unsigned long long& cls, string inname, const long double& bz){

	long double floor=-bz/2.0;
	long double top=bz;	
	
	string povname = inname+"_MinConf.pov";
	ofstream poutfile(povname.c_str());
	cout<<"POVRayOut; Creating file "<<povname<<endl;
	
	poutfile<<"// Render Try: povray -H1000 -W1075 -I"<<povname<<" -O"<<povname<<".png +P +A"<<endl;
	poutfile<<"// "<<endl;

	poutfile<<"#  macro blake_sphere (P1, R1, C1)"<<endl;
	poutfile<<"   sphere {"<<endl;
	poutfile<<"     P1, R1"<<endl;
	poutfile<<"     texture {"<<endl;
	poutfile<<"       pigment { color C1 }"<<endl;
	poutfile<<"       finish { phong 0 phong_size 250 specular 0.3 ambient 0.20 diffuse 0.5  }"<<endl;
	poutfile<<"       normal{ bumps 0.1 omega 1.0 scale 0.002 }"<<endl;
	poutfile<<"  	}"<<endl;
	poutfile<<"  }"<<endl;
	poutfile<<"#end"<<endl;
	poutfile<<endl; 
	poutfile<<"#  macro bond (P1, P2)"<<endl;
	poutfile<<"   cylinder {"<<endl;
	poutfile<<"     P1, P2,"<<endl;
	poutfile<<"     0.300"<<endl;
	poutfile<<"     open"<<endl;
	poutfile<<"      texture { pigment{ color rgb<0.1, 0.1, 0.1> } }"<<endl;
	poutfile<<"   }"<<endl;
	poutfile<<"#end"<<endl;
	poutfile<<" "<<endl;
	
	poutfile<<"camera {"<<endl;
  poutfile<<"    orthographic"<<endl;
  poutfile<<"    location <0.0000, 0.0000, "<<top<<">"<<endl;
  poutfile<<"    look_at <0.0000, 0.0000, "<<floor<<">"<<endl;
 // poutfile<<"    up <0.0000, 3.0000, 0.0000>"<<endl;
 // poutfile<<"    right <2.4065, 0.0000, 0.0000>"<<endl;
poutfile<<"}"<<endl;
poutfile<<"light_source { "<<endl;
poutfile<<"  <-0.1000, 0.1000, "<<top<<">"<<endl;
poutfile<<"  color rgb<1.000, 1.000, 1.000> "<<endl;
poutfile<<"  parallel "<<endl;
poutfile<<"  point_at <0.0, 0.0, "<<floor<<">"<<endl;
poutfile<<"}"<<endl;
poutfile<<"light_source { "<<endl;
poutfile<<"  <20.0000, 20.0000, 10.0000> "<<endl;
poutfile<<"  color rgb<1.000, 1.000, 1.000> "<<endl;
poutfile<<"  parallel "<<endl;
poutfile<<"  point_at <0.0, 0.0, 0.0> "<<endl;
poutfile<<"}"<<endl;
poutfile<<"background {"<<endl;
poutfile<<"  color rgb<1.000, 1.000, 1.000>"<<endl;
poutfile<<"}"<<endl;
poutfile<<endl;

	unsigned long long cind = Eindsort[ccurr];
	unsigned long long sel = (cls)*cind;
	long double R;
	string C;
	bool flag = 0;
	string point1, point2;
	//Type and Coordinates
	for (unsigned j = 1; j <= na; j+=4)
	{
			
			unsigned long long nind=sel+nb+j;
						

			
				if (mat[nind]=="P"){
					R=1.0;
					C="rgbt<0.000,0.900,0.040,0.350>";
					flag=0;
					point1="<"+mat[nind+1]+","+mat[nind+2]+","+mat[nind+3]+">";
				}
				else{
					R=1.359;
					C="rgbt<0.010,0.040,0.930,0.350>";
					flag=1;
					point2="<"+mat[nind+1]+","+mat[nind+2]+","+mat[nind+3]+">";
				}
			
			
				string line;
				line = "blake_sphere(<"+mat[nind+1]+","+mat[nind+2]+","+mat[nind+3]+">,";
				poutfile<<line<<R<<","<<C<<")"<<endl;
				if(flag==1){
					string bline = "bond("+point1+","+point2+")";
					poutfile<<bline<<endl;
				}
			
	}

poutfile<<"// Output from POVRayOut function of pConfigSort_rep_minprov_pre.cpp"<<endl;
poutfile<<"// Edit values in POVRayOut function to change base parameters(color, camera, etc.)"<<endl;
poutfile.close();
cout<<"POVRayOut Complete!"<<endl;
	
}
