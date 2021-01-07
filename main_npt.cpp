#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <time.h>
//#include "conj_grad.h"
#include "soft_sphere_npt.h"

using namespace std;

int main(int argc, char *argv[])
{
    int rand();
    char buffer[64];
    fstream input;
    input.open(argv[1]);
    BOX=stod(argv[2]);
    p_applied=stod(argv[3]);
    double a,b,c,d,e;
	srand( (unsigned)time(NULL) );
	double ke=0,lmx=0,lmy=0;
    while(input>>a>>b>>c>>d)
    {
        X[2*(int(a)-1)]=b;
        X[2*(int(a)-1)+1]=c;

        V[2*(int(a)-1)]  =double(rand())/RAND_MAX-0.5;
        V[2*(int(a)-1)+1]=double(rand())/RAND_MAX-0.5;
		lmx=lmx+V[2*(int(a)-1)];
		lmy=lmy+V[2*(int(a)-1)+1];
		ke =ke +0.5*(V[2*(int(a)-1)]*V[2*(int(a)-1)]+V[2*(int(a)-1)+1]*V[2*(int(a)-1)+1]);
      	//cout<<V[2*(int(a)-1)]<<"\t"<<V[2*(int(a)-1)+1]<<"\t"<<RAD[int(a)-1]<<"\n";
		//cout<<V[2*int(a)-1]*V[2*int(a)-1]+V[2*(int(a)-1)+1]*V[2*(int(a)-1)+1]<<"\n";
        //ATOMz[int(a)-1]=e;
        RAD[int(a)-1]=d;
      //if(b==1)
      //    RAD[int(a)-1]=0.5;
      //if(b==2)
      //    RAD[int(a)-1]=0.7;
      	//cout<<X[2*(int(a)-1)]<<"\t"<<X[2*(int(a)-1)+1]<<"\t"<<RAD[int(a)-1]<<"\n";
		//cout<<(rand())/RAND_MAX<<"\n";
    }
	lmx=lmx/N;
	lmy=lmy/N;
	ke =ke*1./N;
	double fs=sqrt((T)/ke);
	//cout<<fs<<"\n";
	ke=0.;
	pressure_virial1=0.;;
	for(int i=0; i<N; i++)
	{
		V[2*i]  =(V[2*i]  -lmx)*fs;
		V[2*i+1]=(V[2*i+1]-lmy)*fs;
		ke =ke +0.5*(V[2*i]*V[2*i]+V[2*i+1]*V[2*i+1]);
		pressure_virial1 = pressure_virial1 +(V[2*i]*V[2*i]+V[2*i+1]*V[2*i+1]);
      	//cout<<V[2*i]<<"\t"<<V[2*i+1]<<"\n";
	}
	//cout<<ke/N<<"\t"<<N<<"\t"<<ke<<"\n";;
	//return 0;
	
	int ndim=2*N;
//	for(int i=0;i <N; i++)
//	{
//		X[i*2] =ATOMx[i];
//		X[i*2+1]=ATOMy[i];
//		//X[i*3+2]=ATOMz[i];
//	}	
	//X[ndim-1]=BOX;
	make_list();
	cout<<std::setprecision(16);
    cout<<"this "<<energy_force(X,G,ndim)/2000<<"\n";
	pressure_virial=(pressure_virial1-pressure_virial2)/(2.*BOX*BOX);
//  double XX=0.;
//  for(int i=0; i<N; i++)
//  {
//  	//cout<<i<<"\t"<<G[2*i+1]<<"\n";
//  	XX=XX+G[i]*G[i];
//  	//G[i]=GRAD[i];
//  }
	double t_accum=0.;
	double delta_t=0.0001;
	for(int t=0; t<10000000; t++)
	{
		integrate_NPT_SLLOD(X,G,delta_t);
		if(t%5000==0)
		{
			cout<<t*delta_t<<"\n";
			write_config(t,X,"NPT_shear_low_shear");
		}
	}
//////calculate_gradient(X,G,ndim);
////cout<<XX<<"\n";
	//minimize(ndim,X,100000,energy_force,write_config);
}
