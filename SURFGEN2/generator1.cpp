#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

int main () 
{
        double a=1.0; double b=1.0; double c=1.0; //keep a=1 always, change b and c
	int nx=50; int ny=nx*b/a+2; int nz=nx*c/a+2;
	printf("ni,nj,nk=%d, %d, %d \n", nx,ny,nz);
	int sizex=2*nx+1; int sizey=2*ny+1; int sizez=2*nz+1;
	
	double rho0=160.0;
	double rho1=1.0;
	double LL=1.0;
	double mean=rho0;
	double rho;
	double rho_th=0.4;
	double t1;
	
	//double data[size*size*size]; int index;
	//double test=1.5*n+1.5; 

  ofstream myfile,myfile1;
  myfile.open ("oblate_c_by_a_1_50.dat"); myfile.precision(15);myfile1.open ("plot.dat"); myfile1.precision(15);
	myfile<<"#"<<sizex-1<<"\t"<<sizey-1<<"\t"<<sizez-1<<"\t"<<LL<<"\t"<<mean<<endl;
	for(int i=0;i<sizex-1;i++)
	{
		for(int j=0;j<sizey-1;j++)
		{
			for(int k=0;k<sizez-1;k++)
			{
				t1=(i-nx)*(i-nx)/(a*a)+(j-ny)*(j-ny)/(b*b)+(k-nz)*(k-nz)/(c*c);
				rho=rho0/sqrt(t1); //rho decreases as we move away from centre
				//rho=sqrt(t1);	
				if (i==nx && j==ny && k==nz){ rho=rho0/0.0001;} //rho is very high at the centre
				if(i==0 && j==14 && k==14) {cout<<"rho= "<<rho<<endl;}
				//index=k+size*j+size*size*i;
				//data[index]=rho;
				myfile<<rho<<endl;
			       // if(rho>=rho_th){myfile1<<i<<"\t"<<j<<"\t"<<k<<endl;}					
			
			}		
		}
	}

  myfile.close();myfile1.close();
	//cout<<data[141195]<<endl;
	cout<<endl<<"sizex= "<<sizex<<"\t sizey= "<<sizey<<"\t sizez= "<<sizez<<"\t tot= "<<sizex*sizey*sizez<<endl;
	cout<<" rho_min= "<<rho0*a/nx<<endl; // minimum value of rho_th can be given to fit the ellipse in the box
  return 0;
}
