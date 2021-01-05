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
	double ke=0,lmx=0,lmy=0,lmz=0;
    while(input>>a>>b>>c>>d>>e)
    {
		//cout<<a<<"\t"<<b<<"\t"<<c<<"\t"<<d<<"\t"<<e<<"\n";
        X[3*(int(a)-1)]=b;
        X[3*(int(a)-1)+1]=c;
        X[3*(int(a)-1)+2]=d;

        V[3*(int(a)-1)]  =double(rand())/RAND_MAX-0.5;
        V[3*(int(a)-1)+1]=double(rand())/RAND_MAX-0.5;
        V[3*(int(a)-1)+2]=double(rand())/RAND_MAX-0.5;
		lmx=lmx+V[3*(int(a)-1)];
		lmy=lmy+V[3*(int(a)-1)+1];
		lmz=lmz+V[3*(int(a)-1)+2];
		ke =ke +0.5*(V[3*(int(a)-1)]*V[3*(int(a)-1)]+V[3*(int(a)-1)+1]*V[3*(int(a)-1)+1]+V[3*(int(a)-1)+2]*V[3*(int(a)-1)+2]);
      	//cout<<V[2*(int(a)-1)]<<"\t"<<V[2*(int(a)-1)+1]<<"\t"<<RAD[int(a)-1]<<"\n";
		//cout<<V[2*int(a)-1]*V[2*int(a)-1]+V[2*(int(a)-1)+1]*V[2*(int(a)-1)+1]<<"\n";
        //ATOMz[int(a)-1]=e;
        RAD[int(a)-1]=e;
      //if(b==1)
      //    RAD[int(a)-1]=0.5;
      //if(b==2)
      //    RAD[int(a)-1]=0.7;
      	//cout<<X[3*(int(a)-1)]<<"\t"<<X[3*(int(a)-1)+1]<<"\t"<<X[3*(int(a)-1)+2]<<"\t"<<RAD[int(a)-1]<<"\n";
		//cout<<(rand())/RAND_MAX<<"\n";
    }
	lmx=lmx/N;
	lmy=lmy/N;
	lmz=lmz/N;
	ke =ke*1./N;
	double fs=sqrt((T)/(2.*ke/3.));
	//cout<<fs<<"\n";
	ke=0.;
	pressure_virial1=0.;;
	for(int i=0; i<N; i++)
	{
		V[3*i]  =(V[3*i]  -lmx)*fs;
		V[3*i+1]=(V[3*i+1]-lmy)*fs;
		V[3*i+2]=(V[3*i+2]-lmz)*fs;
		ke =ke +0.5*(V[3*i]*V[3*i]+V[3*i+1]*V[3*i+1]+V[3*i+2]*V[3*i+2]);
		pressure_virial1 = pressure_virial1 +(V[3*i]*V[3*i]+V[3*i+1]*V[3*i+1]+V[3*i+2]*V[3*i+2]);
      	//cout<<V[2*i]<<"\t"<<V[2*i+1]<<"\n";
	}
	//cout<<ke/N<<"\t"<<N<<"\t"<<ke<<"\n";;
	//return 0;
	
	int ndim=3*N;
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
	pressure_virial=(pressure_virial1-pressure_virial2)/(3*BOX*BOX*BOX);
//  double XX=0.;
//  for(int i=0; i<N; i++)
//  {
//  	//cout<<i<<"\t"<<G[2*i+1]<<"\n";
//  	XX=XX+G[i]*G[i];
//  	//G[i]=GRAD[i];
//  }
	double t_accum=0.;
	double delta_t=0.0001;
	double g_max=0.035;
	//double g_dot=0.0001;
	double tau=2000;
	double g_dot=(2.*3.14159265359)*g_max/tau;
	//double tau=(2.*3.14159265359)*g_max/g_dot;
	double omega=(2.*3.14159265359)/tau;
	cout<<g_dot<<"\n";
	//double g_max=(g_dot*tau)/2.;
	//cout<<g_max<<"\n";
	for(int t=0; t*delta_t<1000*tau; t++)
	{
		
		double shear_rate = g_dot*cos(omega*t*delta_t);
		double STRAIN = g_max*sin(omega*t*delta_t);
		integrate_NPT_SLLOD(X,G,delta_t,shear_rate,STRAIN);
		if(t%5000==0)
		{
			cout<<t*delta_t<<"\n";
			write_config(t,X,"NPT_shear_cyclic_shear");
		}
	  	if((t%5000)==0)
  		{	
      		cout<<t*delta_t<<"\n";
      		write_config(t,X,"NVT_cyclic_shear_");
  		}
  		if(t%int(0.5*tau/delta_t)==0)
  		{
      		cout<<t*delta_t<<"\n";
      		write_config(int(t/int(0.5*tau/delta_t)),X,"NVT_CYC_STROB_");
  		}

	}
//////calculate_gradient(X,G,ndim);
////cout<<XX<<"\n";
	//minimize(ndim,X,100000,energy_force,write_config);
}
