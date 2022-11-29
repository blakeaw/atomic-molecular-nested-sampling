//=================================
// include guard
#ifndef __BoxNonOrthogonal_H_INCLUDED__
#define __BoxNonOrthogonal_H_INCLUDED__

//=================================
// forward declared dependencies

class Atom;
//=================================
// included dependencies
#include <math.h>
#include <cmath>



//=================================
// the actual class

//3d box that can be nonorthogonal. Used to encapsulate atoms and serve as boundary. 
class BoxNonOrthogonal {
	public:
		//box matrix
		long double length[3][3];
		//inverse box matrix
		long double linverse[3][3];
		//origin 
		long double origin[3];
		// constructor - no args
		BoxNonOrthogonal() {
			
			for (unsigned int i = 0; i < 3; ++i)
			{
				origin[i]=0.0;
			}
			for (unsigned int i = 0; i < 3; ++i)
			{
				for (unsigned int j = 0; j < 3; ++j)
				{
					length[i][j]=0.0;
				}
			}
			//start orthogonal
			length[0][0]=1.0;
			length[1][1]=1.0;
			length[2][2]=1.0;
			ComputeBoxInverse();
		}
	
	//Destructor
//		~BoxNonOrthogonal(){
//			
////			if (center){
////				
////				delete[] center;
////				
////			}
////			if (length){
////				delete[] length;
////			}
//			
//		}
	//copy constructor
	BoxNonOrthogonal(BoxNonOrthogonal& other){
			for (unsigned int i = 0; i < 3; ++i)
			{
				origin[i]=other.origin[i];
			}
			for (unsigned int i = 0; i < 3; ++i)
			{
				for (unsigned int j = 0; j < 3; ++j)
				{
					length[i][j]=other.length[i][j];
				}
			}
			ComputeBoxInverse();
	}	

	//copy

	void Copy(BoxNonOrthogonal& other){

			for (unsigned int i = 0; i < 3; ++i)
			{
				origin[i]=other.origin[i];
			}
			for (unsigned int i = 0; i < 3; ++i)
			{
				for (unsigned int j = 0; j < 3; ++j)
				{
					length[i][j]=other.length[i][j];
				}
			}
			ComputeBoxInverse();
	}
	
	void ComputeBoxInverse(){

		long double a, b,c,d,e,f,g,h,i;
		a=length[1][1]*length[2][2]-length[1][2]*length[2][1];
		b=length[0][2]*length[2][1]-length[2][2]*length[0][1];
		c=length[0][1]*length[1][2]-length[1][1]*length[0][2];
		d=length[1][2]*length[2][0]-length[1][0]*length[2][2];
		e=length[0][0]*length[2][2]-length[0][2]*length[2][0];
		f=length[0][2]*length[1][0]-length[0][0]*length[1][2];
		g=length[1][0]*length[2][1]-length[1][1]*length[2][0];
		h=length[0][1]*length[2][0]-length[0][0]*length[2][1];
		i=length[0][0]*length[1][1]-length[0][1]*length[1][0];
		long double detinv = 1.0/ComputeVolume();
		linverse[0][0]=a*detinv;
		linverse[0][1]=b*detinv;
		linverse[0][2]=c*detinv;
		linverse[1][0]=d*detinv;
		linverse[1][1]=e*detinv;
		linverse[1][2]=f*detinv;
		linverse[2][0]=g*detinv;
		linverse[2][1]=h*detinv;
		linverse[2][2]=i*detinv;
		return;
		
	}
	//multiply a vector by the inverse box matrix, i.e. linverse*vec, and replace vec with outcome
	//Used to convert coordinates into box scaled coordinates
	void MultiplyCoordVecByBoxInverse(long double mat[]){
		long double x,y,z;
		//ComputeBoxInverse();
		x=linverse[0][0]*mat[0]+linverse[0][1]*mat[1]+linverse[0][2]*mat[2];
		y=linverse[1][0]*mat[0]+linverse[1][1]*mat[1]+linverse[1][2]*mat[2];
		z=linverse[2][0]*mat[0]+linverse[2][1]*mat[1]+linverse[2][2]*mat[2];
		mat[0]=x;
		mat[1]=y;
		mat[2]=z;
		return;
	}
	//multiply a vector by the box matrix, i.e. length*vec, and replace vec with outcome
	//used to convert coordinates from box scaled to unscaled
	void MultiplyCoordVecByBox(long double mat[]){
		long double x,y,z;
		x=length[0][0]*mat[0]+length[0][1]*mat[1]+length[0][2]*mat[2];
		y=length[1][0]*mat[0]+length[1][1]*mat[1]+length[1][2]*mat[2];
		z=length[2][0]*mat[0]+length[2][1]*mat[1]+length[2][2]*mat[2];
		mat[0]=x;
		mat[1]=y;
		mat[2]=z;
		return;
	}

	bool CheckForNegativeElement(){

		for (unsigned int i = 0; i < 3; ++i)
			{
				for (unsigned int j = 0; j < 3; ++j)
				{
					if(length[i][j]<0.0){
						
						return 1;

					}
				}
			}
		return 0;
	}
	//Enforce lammps style triclinic box
	bool CheckForNegativeElementLammps(){
		int ifer = 0;
		for (unsigned int i = 0; i < 3; ++i)
			{
				for (unsigned int j = 0; j < 3; ++j)
				{
					if (i==0 && j==1){
						++ifer;
					}
					if(i==0 && j==2){
						++ifer;
					}
					else if(i==1 && j==2){
						++ifer;
					}
					else if(length[i][j]<0.0){
						return 1;						

					}
				}
			}
		return 0;
	}
	void CoutElements(){
		for (unsigned int i = 0; i < 3; ++i)
		{
			for (unsigned int j = 0; j < 3; ++j)
			{
				cout<<length[i][j]<<" ";
			}
			cout<<endl;
		}
	}
	void CoutInverseElements(){
		ComputeBoxInverse();
		for (unsigned int i = 0; i < 3; ++i)
		{
			for (unsigned int j = 0; j < 3; ++j)
			{
				cout<<linverse[i][j]<<" ";
			}
			cout<<endl;
		}
	}
	void CoutDimensions(){
		long double a,b,c;
		a=GetA();
		b=GetB();
		c=GetC();
		cout<<"box dimensions; a: "<<a<<" b: "<<b<<" c: "<<c<<endl;
	}
	long double ComputeVolume(){
			//volume is determinant of box matrix: det(length)
			long double a = length[0][0]*( length[1][1]*length[2][2]- length[1][2]*length[2][1] );
			long double b = length[0][1]*( length[1][0]*length[2][2]- length[1][2]*length[2][0] );
			long double c = length[0][2]*( length[1][0]*length[2][1]- length[1][1]*length[2][0] );
			long double V = a-b+c;
			return V;
	}
	void AddToBox(long double mat[3][3]){
		for (unsigned int i = 0; i < 3; ++i)
		{
			for (unsigned int j = 0; j < 3; ++j)
			{
				length[i][j]+=mat[i][j];
			}
		}
		ComputeBoxInverse();
		return;
	}	

	void SetOrthoLengths(long double x, long double y, long double z){
			length[0][0]=x;
			length[1][1]=y;
			length[2][2]=z;
			ComputeBoxInverse();
			return;
	}

	void SetOrigin(long double x, long double y, long double z){
			origin[0]=x;
			origin[1]=y;
			origin[2]=z;
			return;
	}

	//compute the components or box matrix [a b c]
	long double GetA(void){
		long double a = sqrt(pow(length[0][0], 2)+pow(length[1][0],2)+pow(length[2][0],2));
		return a;

	}
	long double GetB(void){
		long double b = sqrt(pow(length[0][1], 2)+pow(length[1][1],2)+pow(length[2][1],2));
		return b;

	}
		long double GetC(void){
		long double c = sqrt(pow(length[0][2], 2)+pow(length[1][2],2)+pow(length[2][2],2));
		return c;

	}
	//coordinate wrapping
	void WrapCoordinates(long double& x , long double& y, long double& z){
		
		//store coordinates in a vector
		long double cvec[3];
		cvec[0]=x;
		cvec[1]=y;
		cvec[2]=z;
		//convert to scaled coordinates
		ComputeBoxInverse();
		MultiplyCoordVecByBoxInverse(cvec);
		//Wrap
		//x
		
		while (cvec[0]>1.0){
			cvec[0]-=1.0;
		}
		while(cvec[0]<0.0){
			cvec[0]+=1.0;
		}
		//y
		
		while (cvec[1]>1.0){
			cvec[1]-=1.0;
		}
		while(cvec[1]<0.0){
			cvec[1]+=1.0;
		}
		//z
		
		while(cvec[2]>1.0){
			cvec[2]-=1.0;
		}
		while(cvec[2]<0.0){
			cvec[2]+=1.0;
		}
		//Convert back to unscaled
		MultiplyCoordVecByBox(cvec);
		//set new value
		x=cvec[0];
		y=cvec[1];
		z=cvec[2];
		return;
	}
		
	void WrapAtom(Atom& a){
		long double ax, ay, az;
		ax=a.x;
		ay=a.y;
		az=a.z;
		WrapCoordinates(ax, ay, az);
		a.x=ax;
		a.y=ay;
		a.z=az;
		return;

	}
		
	void ScaleCoordinates(long double factor, Config& c){
		long double cvec[3];
		ComputeBoxInverse();
		for (unsigned int i = 0; i < c.natoms; ++i)
		{
			cvec[0]=c.atom[i].x;
			cvec[1]=c.atom[i].y;
			cvec[2]=c.atom[i].z;
			MultiplyCoordVecByBoxInverse(cvec);
			cvec[0]*=factor;
			cvec[1]*=factor;
			cvec[2]*=factor;
			MultiplyCoordVecByBox(cvec);
			c.atom[i].x=cvec[0];
			c.atom[i].y=cvec[1];
			c.atom[i].z=cvec[2];
		}

	}

	void ScaleCoordinates(long double fx,long double fy,long double fz, Config& c){
		long double cvec[3];
		ComputeBoxInverse();
		for (unsigned int i = 0; i < c.natoms; ++i)
		{
			cvec[0]=c.atom[i].x;
			cvec[1]=c.atom[i].y;
			cvec[2]=c.atom[i].z;
			MultiplyCoordVecByBoxInverse(cvec);
			cvec[0]*=fx;
			cvec[1]*=fy;
			cvec[2]*=fz;
			MultiplyCoordVecByBox(cvec);
			c.atom[i].x=cvec[0];
			c.atom[i].y=cvec[1];
			c.atom[i].z=cvec[2];
		}

	}
	void ScaleCoordinatesByMatrix(Config& c, BoxNonOrthogonal& other){
		ComputeBoxInverse();
		other.ComputeBoxInverse();
		long double scalingmatrix[3][3];
		for (unsigned int i = 0; i < 3; ++i)
		{
			for (unsigned int j = 0; j < 3; ++j)
			{
				scalingmatrix[i][j]=0.0;
			}
		}
		for (int row = 0; row < 3; row++) {
		    for (int col = 0; col < 3; col++) {
		        // Multiply the row of A by the column of B to get the row, column of product.
		        for (int inner = 0; inner < 3; inner++) {
		            scalingmatrix[row][col] += length[row][inner] * other.linverse[inner][col];
		        }
		        //std::cout << scalingmatrix[row][col] << "  ";
		    }
		  //  std::cout << "\n";
		}
		for (unsigned int i = 0; i < c.natoms; ++i)
		{
			long double mat[3];
			mat[0]=c.atom[i].x;
			mat[1]=c.atom[i].y;
			mat[2]=c.atom[i].z;
//			cout<<"mat before scaling matrix: x "<<mat[0]<<" y "<<mat[1]<<" z "<<mat[2]<<endl;
			//other.MultiplyCoordVecByBoxInverse(mat);
			long double x,y,z;
			x=scalingmatrix[0][0]*mat[0]+scalingmatrix[0][1]*mat[1]+scalingmatrix[0][2]*mat[2];
			y=scalingmatrix[1][0]*mat[0]+scalingmatrix[1][1]*mat[1]+scalingmatrix[1][2]*mat[2];
			z=scalingmatrix[2][0]*mat[0]+scalingmatrix[2][1]*mat[1]+scalingmatrix[2][2]*mat[2];
			mat[0]=x;
			mat[1]=y;
			mat[2]=z;
//			cout<<"mat after scaling matrix: x "<<mat[0]<<" y "<<mat[1]<<" z "<<mat[2]<<endl;
//			exit(0);
			//MultiplyCoordVecByBox(mat);
			c.atom[i].x=mat[0];
			c.atom[i].y=mat[1];
			c.atom[i].z=mat[2];
		}

			
    
	}


		
};


#endif // __BoxNonOrthogonal_H_INCLUDED__ 
