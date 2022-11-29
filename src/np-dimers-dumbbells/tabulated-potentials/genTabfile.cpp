#include <stdlib.h>
#include <string>
#include <fstream>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <ctime>

using namespace std;
// Script used to generate energies (Colloid, Charge, and Combo)

// Define Classes

// Prototype Functions

long double PairColloid(long double r, long double rs, long double R1, long double R2);
long double PairChargeShell(long double r, long double rs, long double R1, long double R2, long double kappa, long double Ae);
long double ForceColloid(long double r, long double rs, long double R1, long double R2);
long double ForceSE(long double r, long double rs, long double R1, long double R2, long double kappa, long double Ae);


// Main Function
int main(int argc, char *argv[])
{
unsigned Rloop = 2000;
long double Rstart, Rend, r, rs;
long double pColloid, pSE, fColloid, fSE;
long double pCombo, fCombo;
long double R1, R2, Ae, kappa;
long double Ae1, Ae2;
//Radii of particles--stays the same
R1 = 1.0;
R2 = 1.0;

// inverse screening distance -- set to specific value
kappa= 1.0;

//Prefactor for Screened Electrostatics--
Ae1 = 1000.0;
Ae2 = 1000.0;
//Determine Total Prefactor -- Ae
//set chargetag - 0:same charge ; 1:opposite charges
bool chargetag =0;

if(chargetag == 0){
	//Same charge (e.g. positive-positive) -- geometrical mixing 
	Ae = sqrt(Ae1*Ae2);
}
if(chargetag == 1){
	//Opposite Charges (e.g. positive-negative) -- negative geometrical mixing
	Ae = -sqrt(Ae1*Ae2);
}

Rstart = R1+R2+0.00000000001;
Rend = 10.0;
cout<<"Ae: "<<Ae<<endl;
ofstream outfile2("Colloid_test_R1R1.tab");
ofstream outfile3("Charge_r1r1Ae5000_5000kappa0.10.tab");
ofstream outfile4("Combo_bbr1r1Ae4000_4000kappa0.50.tab");

outfile2 << "# Colloid r, U, F " << endl;
outfile3 << "# ChargeShell r, U, F " << endl;
outfile4 << "# Combo of Colloid and Charge " << endl;
outfile4 << "START" << endl;
outfile4 << "N " << Rloop << " R " << Rstart << " " << Rend << endl;
outfile4 << " " << endl;

for(unsigned R=0;R<Rloop;R++){
	r = ( Rstart + ((Rend-Rstart)/((long double)Rloop))*R );
	rs = pow(r, 2);
	unsigned Rp1 = R+1;
	pColloid = PairColloid(r, rs, R1, R2);
	pSE = PairChargeShell(r, rs, R1, R2, kappa, Ae);
	pCombo = pColloid+pSE;
	fColloid = -ForceColloid(r, rs, R1, R2);
	fSE = -ForceSE(r, rs, R1, R2, kappa, Ae);
	fCombo = fColloid+fSE;
//	fCombo *= -1.0;

	outfile2 << std::scientific << setiosflags(ios::fixed) << setprecision(10) << r << " " << pColloid << " " << fColloid << endl;
	
	outfile3 << std::scientific << setiosflags(ios::fixed) << setprecision(10) << r << " " << pSE << " " << fSE << endl;

	
	outfile4 << std::scientific << setiosflags(ios::fixed) << setprecision(10) << Rp1 << " " << r << " " << pCombo << " " << fCombo << endl;
}

outfile2.close();
outfile3.close();
outfile4.close();
return EXIT_SUCCESS;

}

// Define Functions


long double PairColloid(long double r, long double rs, long double R1, long double R2){
	
	long double pi = 3.14159265359;
	long double potent, Pa, Pb;
	long double R1s, R2s;
	long double Acc,sigma, F;
	//Acc = 39.5;
	//sigma = 1.0;
	//Reduced Units
	Acc = 1.21;
	F = 20.0;
	R1s = pow(R1, 2);
	R2s = pow(R2, 2);
	
	
	
		Pa = -(Acc/6.0) * ( ( (2.0*R1*R2)/(pow(r, 2)-pow(R1+R2, 2)) ) + ( (2.0*R1*R2)/(pow(r, 2)-pow(R1-R2, 2)) ) + log((pow(r, 2)-pow(R1+R2, 2))/(pow(r, 2)-pow(R1-R2, 2))) );

		Pb = ((Acc)/(37800.0*r*pow(F, 6))) * ( ( (pow(r, 2)-7.0*r*(R1+R2)+6.0*(R1s+7.0*R1*R2+R2s))/(pow(r-R1-R2, 7))) + ( (pow(r, 2)+7.0*r*(R1+R2)+6.0*(R1s+7.0*R1*R2+R2s))/(pow(r+R1+R2, 7))) - ( (pow(r, 2)+7.0*r*(R1-R2)+6.0*(R1s-7.0*R1*R2+R2s))/(pow(r+R1-R2, 7))) - ( (pow(r, 2)-7.0*r*(R1-R2)+6.0*(R1s-7.0*R1*R2+R2s))/(pow(r-R1+R2, 7))) ) ;
		potent = Pa+Pb;
	
	
	return(potent);
}

long double PairChargeShell(long double r, long double rs, long double R1, long double R2, long double kappa, long double Ae){
	
	long double potent;

	
	potent = (Ae/(r))*(exp(-kappa*(R1+R2+r)))*(1.0+exp(-2.0*kappa*R1))*(1.0+exp(-2.0*kappa*R2));
	
	
	return(potent);
}

long double ForceColloid(long double r, long double rs, long double R1, long double R2){
	long double dea, deb, debb, dec, ded, dee, def, deg, deh;
	long double dei, dej, dek, del, dem, den, dp, der, de2, f2;

	long double Acc=1.0, F=50.0;			

	dea=Acc/(37800.0*pow(r, 2)*pow(F, 6));
	deb=(403200.0*R1*R1*R1*R2*R2*R2*r*r*r*pow(F, 6))/(pow(R1*R1-2*R1*R2+R2*R2-r*r, 2));
	debb=pow((R1*R1+2*R1*R2+R2*R2-r*r), 2);
	dec=(6.0*(R1*R1-7.0*R1*R2+R2*R2)+7.0*(R1-R2)*r+r*r)/pow(R1-R2+r, 7);
	ded=(6.0*(R1*R1+7.0*R1*R2+R2*R2)-7.0*(R1+R2)*r+r*r)/pow(R1+R2-r, 7);
	dee=(6.0*(R1*R1+7.0*R1*R2+R2*R2)+7.0*(R1+R2)*r+r*r)/pow(R1+R2+r, 7);
	def=(6.0*R1*R1+6.0*R2*R2+7.0*R2*r+r*r-7.0*R1*(6.0*R2+r))/pow(-R1+R2+r, 7);
	deg=r;
	deh=(7.0*R1+7.0*R2-2*r)/pow(R1+R2-r, 7);
	dei=(-7.0*R1+7.0*R2-2*r)/pow(R1-R2+r, 7);
	dej=(-7.0*R1+7.0*R2+2*r)/pow(R1-R2-r, 7);
	dek=(7.0*R1+7.0*R2+2.0*r)/pow(R1+R2+r, 7);
	del=(7.0*(6.0*(R1*R1-7.0*R1*R2+R2*R2)+7.0*(R1-R2)*r+r*r))/pow(R1-R2+r, 8);
	dem=(7.0*(6.0*(R1*R1+7.0*R1*R2+R2*R2)-7.0*(R1+R2)*r+r*r))/pow(R1+R2-r, 8);
	den=(7.0*(6.0*(R1*R1+7.0*R1*R2+R2*R2)+7.0*(R1+R2)*r+r*r))/pow(R1+R2+r, 8);
	dp=(7.0*(6.0*R1*R1+6.0*R2*R2+7.0*R2*r+r*r-7.0*R1*(6.0*R2+r)))/pow(-R1+R2+r, 8);
	
	der=(deh+dei+dej+dek+del-dem-den+dp);
        de2=dea*(deb/debb+dec+ded-dee+def+deg*der);
        f2=-de2;
	return(f2);

}

long double ForceSE(long double r, long double rs, long double R1, long double R2, long double kappa, long double Ae){
		long double j1, j2, j3, jj;
		long double de, f;

		j1= exp(-kappa*(R1+R2+r));
		j2 = (-1.0+exp(-2.0*kappa*R1));
		j3 = (-1.0+exp(-2.0*kappa*R2));
		jj = j1*j2*j3;
		de = -Ae*((jj/rs) + (jj/r)*kappa);
		f = -de;

		return(f);
}

