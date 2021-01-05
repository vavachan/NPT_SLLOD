#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <iomanip>

using namespace std;

// double BOX=2*7.2005488787721159;
// double BOX=2*26.2107109945647245;
double BOX;
int N=2000;
double epsilon=1.0;
double magnitude=0.;
double p_applied=0.002;
double KB=1.0;
double Q=0.005;
double iQ=1./Q;
double T=0.0001;
double STRAIN=0.0;
double zeta=0.;
double alpha=0.;
double tilt=STRAIN*BOX;
double PEnergy=0.;
double KEnergy=0.;
int dim=3;

//double *ATOMx = new  double[N];
//double *ATOMy = new  double[N];
// double *ATOMz = new  double[N];
double *RAD   = new  double[N];
//static  double *DIST= new  double[N];
static int *VLIST= new int [N*(N+1)];

double stress_tensor[3][3];
double pressure_virial2=0.;
double pressure_virial1=0.;
double pressure_virial=0.;

//static  double *tATOMx = new  double[N];
//static  double *tATOMy = new  double[N];
//static  double *tATOMz = new  double[N];
// gradient is calculated and stored in this 

//double *CONJ_D= new  double[2*N+1];
static double *prevG= new  double[dim*N];
double *G= new  double[dim*N];
double *X= new  double[dim*N];
double *V= new  double[dim*N];
static double *tX= new  double[dim*N];

double deltadrv=0.5;
double deltadrvsq=deltadrv*deltadrv;


void make_list()
{
	for(int i=0; i<dim*N; i++)
	{
		tX[i]=-999;
		//tATOMy[i]=-999;
//		tATOMz[i]=ATOMz[i];
	}
	//cout<<tATOMx[0]<<"\t"<<ATOMx[0]<<" htis \n";
}

double energy_force( double X[], double GRAD [], int ndim)
{
	
	//BOX=X[ndim-1];
	clock_t begin=clock();
	double delta_rx,delta_ry,delta_rz;
	double *dR = new double[dim];
	double delta_r;
	double pot_e=0;
    double skin_depth=1.4+deltadrv;
    double skin_depth_sq=skin_depth*skin_depth;
    int nnei=0;
    int nlistbeg=0;
	double iBOX=1./BOX;
	
	double maxdisp=0.;
	double maxdisp2=0.;
	//cout<<tATOMx[0]<<"\t"<<X[0]<<"\n";
	for(int i=0; i<N; i++)
	{
		for(int d=0; d<dim;d++)
		{
			GRAD[i*dim+d] = 0.;
			dR[d]=tX[i*dim+d]-X[i*dim+d];
		}
		delta_r=0.;
		for(int d=0; d<dim;d++)
		{
			if(tilt and d==0)
				dR[d]=dR[d]-tilt*lroundl(dR[d+1]*iBOX);
			dR[d]=dR[d]-BOX*lroundl(dR[d]*iBOX);
			delta_r=delta_r+dR[d]*dR[d];
		}
		
		delta_r=sqrt(delta_r);

		if(delta_r > maxdisp)
		{
			maxdisp2=maxdisp;
			maxdisp=delta_r;
			
		} 
		else if(delta_r > maxdisp2)
		{
			maxdisp2=delta_r;
		}
		//cout<<ATOMx[i]<<"\n";
	}	
	//cout<<"( energyforce ) \t"<<maxdisp<<"\t"<<maxdisp2<<"\t"<<maxdisp+maxdisp2<<"\t"<<deltadrv<<"\n" ;
	if(maxdisp+maxdisp2 > deltadrv )
	{
		//cout<<"inside ( energyforce ) \t"<<maxdisp<<"\t"<<maxdisp2<<"\t"<<maxdisp+maxdisp2<<"\t"<<deltadrv<<"\n" ;
		//cout<<tATOMx[0]<<"\t"<<X[0]<<"\n";
		//make_list();
		for(int i=0; i<N; i++)
		{
			nnei=0;
			
			for(int j=0 ; j<N; j++)
			{
				if ( i != j )
				{
					for(int d=0; d<dim;d++)
					{
						dR[d]=X[i*dim+d]-X[j*dim+d];
					}
					delta_r=0.;
					for(int d=0; d<dim;d++)
					{
						if(tilt and d==0)
							dR[d]=dR[d]-tilt*lroundl(dR[d+1]*iBOX);
						dR[d]=dR[d]-BOX*lroundl(dR[d]*iBOX);
						delta_r=delta_r+dR[d]*dR[d];
					}
					
					//delta_r=sqrt(delta_r);
					if(delta_r< skin_depth_sq)
					{
						nnei=nnei+1;
						VLIST[nlistbeg+nnei]=j;
					}
				}
			}
			VLIST[nlistbeg]=nnei;
			nlistbeg=nlistbeg+nnei+1;
		}
		for(int i=0; i<dim*N; i++)
		{
			tX[i]=X[i];
		}
	}
	nlistbeg=0;

	for(int d=0; d<dim;d++)
	{
		for(int d1=d; d1<dim;d1++)
		{
			stress_tensor[d][d1]=0.;
		}
	}
	pressure_virial2=0.;
	for(int i=0; i<N; i++)
	{
		int jj;
		//cout<<i<<"\n";
		for(jj=nlistbeg+1; jj<nlistbeg+VLIST[nlistbeg]+1; jj++)
		{
			int	j=VLIST[jj];
			if(j>i)
			{
				//delta_rx=(X[2*i]-X[2*j]);
				for(int d=0; d<dim;d++)
				{
					dR[d]=X[i*dim+d]-X[j*dim+d];
				}
				delta_r=0.;
				for(int d=0; d<dim;d++)
				{
					if(tilt and d==0)
						dR[d]=dR[d]-tilt*lroundl(dR[d+1]*iBOX);
					dR[d]=dR[d]-BOX*lroundl(dR[d]*iBOX);
					delta_r=delta_r+dR[d]*dR[d];
				}
				
				//delta_r=sqrt(delta_r);
				double sigma=RAD[i]+RAD[j];
				double sigmasq=sigma*sigma;
				//cout<<i<<"\t"<<j<<"\t"<<delta_r<<"\n";
				if(delta_r<sigmasq)
				{
					double isigma=1./sigma;
					delta_r=sqrt(delta_r);
					double vh=(1-delta_r*isigma);
					pot_e=pot_e+epsilon*vh*vh;
					//dwreturn epsilon*2.*(1-r/sigma)*(-1./sigma);
					double U=epsilon*2.*vh*(-1.*isigma)/delta_r;

					for(int d=0; d<dim;d++)
					{
						GRAD[i*dim+d]   = GRAD[i*dim+d]   + U*dR[d];
						GRAD[j*dim+d]   = GRAD[j*dim+d]   - U*dR[d];
						
						for(int d1=d; d1<dim;d1++)
						{
							stress_tensor[d][d1]=stress_tensor[d][d1]+(U*dR[d])*dR[d1];
						}
						pressure_virial2 = pressure_virial2+U*dR[d]*dR[d];
					}
				}
			}
		}
		nlistbeg=jj;//nlistbeg+jj+1;
	}
	for(int d=0; d<dim;d++)
	{
		for(int d1=d; d1<dim;d1++)
		{
			stress_tensor[d][d1]=stress_tensor[d][d1]/N;
		}
	}

	clock_t end=clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	//cout<<time_spent<<" in energy \n";
	return pot_e;
}
void NPT_mVV( double X[], double GRAD[], double delta_t, double v_e, double shear_rate, double STRAIN)
{
	double delta_t4=delta_t/4.;
	double delta_t2=delta_t/2.;
	double delta_t8=delta_t/8.;
	
	for( int i=0; i < N; i++)
	{
		for( int d=0; d< dim ; d++)
		{
			V[dim*i+d]  =V[dim*i+d]  + delta_t2*(-1.*GRAD[dim*i+d]);
		//V[2*i+1]=V[2*i+1]+delta_t2*(-1.*GRAD[2*i+1]);
		}
	}
	double scaley1=exp(v_e*delta_t2);;
	double scaley2=exp(v_e*delta_t4);;
	double sin_hypy=sinh(delta_t4*v_e);

	double scalex1=exp(v_e*delta_t);;
	double scalex2=exp(v_e*delta_t2);;
	double sin_hypx=sinh(delta_t2*v_e);

	if(dim==2)
	{
		for( int i=0; i < N; i++)
		{
			//X[2*i+1] = X[2*i+1]*scaley1 + 4.*(V[2*i+1])*scaley2*sin_hypy;
			X[2*i+1] = X[2*i+1]*scaley1 + 2.*(V[2*i+1]/v_e)*scaley2*sin_hypy;
			X[2*i]   = X[2*i]*scalex1   + 2.*((V[2*i]+shear_rate*X[2*i+1])/v_e)*scalex2*sin_hypx;
			//X[2*i+1] = X[2*i+1]*scaley1 + 2.*V[2*i+1]*scaley2*sin_hypy;
			X[2*i+1] = X[2*i+1]*scaley1 + 2.*(V[2*i+1]/v_e)*scaley2*sin_hypy;
			//X[2*i+1] = X[2*i+1]*scaley1 + 2.*V[2*i+1]*scaley2*sin_hypy;
		}
	}
	if(dim==3)
	{
		for( int i=0; i < N; i++)
		{
			//X[2*i+1] = X[2*i+1]*scaley1 + 4.*(V[2*i+1])*scaley2*sin_hypy;
			X[3*i+1] = X[3*i+1]*scaley1 + 2.*(V[3*i+1]/v_e)*scaley2*sin_hypy;
			X[3*i]   = X[3*i]*scalex1   + 2.*((V[3*i]+shear_rate*X[3*i+1])/v_e)*scalex2*sin_hypx;
			X[3*i+2] = X[3*i+2]*scalex1 + 2.*((V[3*i+2]/v_e)*scalex2*sin_hypx);
			//3[2*i+1] = 3[2*i+1]*scaley1 + 2.*V[2*i+1]*scaley2*sin_hypy;
			X[3*i+1] = X[3*i+1]*scaley1 + 2.*(V[3*i+1]/v_e)*scaley2*sin_hypy;
			//X[2*i+1] = X[2*i+1]*scaley1 + 2.*V[2*i+1]*scaley2*sin_hypy;
		}
	}

	BOX=BOX*exp(delta_t*v_e);

	tilt=STRAIN*BOX;
	PEnergy=energy_force(X,GRAD,2*N);
	KEnergy=0.;
	
	double Vol=1.;
	for( int d=0; d< dim ; d++)
	{
		Vol=Vol*BOX;
	}
	pressure_virial1=0.;
	for(int i=0; i <N; i++)
	{
		for( int d=0; d< dim ; d++)
		{
			V[dim*i+d]  =V[dim*i+d]  + delta_t2*(-1.*GRAD[dim*i+d]);
		//V[2*i+1]=V[2*i+1]+delta_t2*(-1.*GRAD[2*i+1]);
			pressure_virial1=pressure_virial1+(V[dim*i+d]*V[dim*i+d]);
			KEnergy=KEnergy+0.5*(V[dim*i+d]*V[dim*i+d]);
		}

		//V[2*i]  =V[2*i]  +delta_t2*(-1.*GRAD[2*i]  );
		//V[2*i+1]=V[2*i+1]+delta_t2*(-1.*GRAD[2*i+1]);
	}
	pressure_virial=(pressure_virial1-pressure_virial2)/(dim*Vol);
		
}


void NH_NPT_SLLOD(double X[], double GRAD[], double delta_t, double &eta, double &eta_e, double &v_eta, double &v_e , double &v_eta_e, double W, double Q,double Qe, double &KE, double shear_rate)
{

	double L=dim*N;
	double G_eta_e=W*v_eta*v_eta-T;
	double G=2.*KE-L*T;
	double Vol=1.;// BOX*BOX*BOX;
	
	for( int d=0; d< dim ; d++)
	{
		Vol=Vol*BOX;
	}
	double Ge=dim*Vol*(pressure_virial-p_applied)+(dim/L)*(2.*KE);

	double delta_t4=delta_t/4.;
	double delta_t2=delta_t/2.;
	double delta_t8=delta_t/8.;

	v_eta = v_eta + delta_t4*G/Q;

	v_eta_e = v_eta_e + delta_t4*G_eta_e/Qe;

	v_e = v_e*exp(-1.*delta_t8*v_eta_e);
	v_e = v_e + delta_t4*Ge/W;
	v_e = v_e*exp(-1.*delta_t8*v_eta_e);
	
	eta = eta + delta_t2*v_eta;
	eta_e = eta_e +delta_t2*v_eta_e;

	double scalex= exp(-1.*delta_t8*(v_eta +(1+dim/L)*v_e));
	double scaley= exp(-1.*delta_t2*(v_eta +(1+dim/L)*v_e));
	double scalez= exp(-1.*delta_t2*(v_eta +(1+dim/L)*v_e));
	pressure_virial1=0.;
	KE=0.;
	if(dim==2)
	{
		for(int i=0; i <N; i++)
		{
			V[2*i] = V[2*i]*scalex;
			V[2*i] = V[2*i] - shear_rate*delta_t4*V[2*i+1];
			V[2*i] = V[2*i]*scalex;

			V[2*i+1] = V[2*i+1]*scaley;

			V[2*i] = V[2*i]*scalex;
			V[2*i] = V[2*i] - shear_rate*delta_t4*V[2*i+1];
			V[2*i] = V[2*i]*scalex;

			pressure_virial1=pressure_virial1+(V[2*i]*V[2*i]+V[2*i+1]*V[2*i+1]);
			KE=KE+0.5*(V[2*i]*V[2*i]+V[2*i+1]*V[2*i+1]);
		}
	}
	if(dim==3)
	{
		for(int i=0; i <N; i++)
		{
			V[3*i] = V[3*i]*scalex;
			V[3*i] = V[3*i] - shear_rate*delta_t4*V[3*i+1];
			V[3*i] = V[3*i]*scalex;

			V[3*i+1] = V[3*i+1]*scaley;
			V[3*i+2] = V[3*i+2]*scalez;

			V[3*i] = V[3*i]*scalex;
			V[3*i] = V[3*i] - shear_rate*delta_t4*V[3*i+1];
			V[3*i] = V[3*i]*scalex;

			pressure_virial1=pressure_virial1+(V[3*i]*V[3*i]+V[3*i+1]*V[3*i+1]+V[3*i+2]*V[3*i+2]);
			KE=KE+0.5*(V[3*i]*V[3*i]+V[3*i+1]*V[3*i+1]+V[3*i+2]*V[3*i+2]);
		}
	}
	pressure_virial=(pressure_virial1-pressure_virial2)/(dim*Vol);

	Ge=dim*Vol*(pressure_virial-p_applied)+(dim/L)*(2.*KE);
	
	v_e = v_e*exp(-1.*delta_t8*v_eta_e);
	v_e = v_e + delta_t4*Ge/W;
	v_e = v_e*exp(-1.*delta_t8*v_eta_e);
	
	G=2.*KE-L*T;
	G_eta_e=W*v_eta*v_eta-T;
	
	
	v_eta_e = v_eta_e + delta_t4*G_eta_e/Qe;

	v_eta = v_eta + delta_t4*G/Q;
	
}

void integrate_NPT_SLLOD(double X[], double GRAD[], double delta_t, double shear_rate, double STRAIN)
{
	static double eta=0.;
	static double eta_e=0.;
	static double v_eta=0.;
	static double v_eta_e=0.;
	static double v_e=0.;
	//double shear_rate=0.0001;
	double W=10.;
	double Q=10.;
	double Qe=10.;
	//STRAIN=STRAIN+shear_rate*delta_t;
	NH_NPT_SLLOD(X,GRAD,delta_t, eta, eta_e, v_eta, v_e, v_eta_e, W,Q,Qe, KEnergy, shear_rate);
	NPT_mVV(X,GRAD,delta_t,v_e,shear_rate,STRAIN);
	NH_NPT_SLLOD(X,GRAD,delta_t, eta, eta_e, v_eta, v_e, v_eta_e, W,Q,Qe, KEnergy, shear_rate);
	//cout<<STRAIN<<"\t"<<PEnergy<<"\t"<<KEnergy<<"\t"<<KEnergy/(N)<<"\t"<<T<<"\t"<<stress_tensor[0][0]<<"\t"<<stress_tensor[0][1]<<"\t"<<stress_tensor[1][1]<<"\t"<<pressure_virial<<"\t"<<p_applied<<"\t"<<3.14159265359*1000*(0.5*0.5+0.7*0.7)/BOX*BOX<<"\n";
}
void write_config(int k,  double X[], char name[])
{
////for(int i=0;i <N; i++)
////{
////	ATOMx[i]=X[i*3];
////	ATOMy[i]=X[i*3+1];
////	ATOMz[i]=X[i*3+2];
////}	
	fstream input(std::string(name)+"_Config"+"_"+std::to_string(k)+".xyz", fstream::out);
	input<<N<<"\n";
	input<<"\n";
	input<<std::setprecision(16);
	//double pressure_virial=0.;
	for(int i=0; i<N; i++)
	{
		if(tilt)
			X[3*i]=(X[3*i]-(tilt*lroundf(X[3*i+1]/BOX)));
		X[3*i]  =(X[3*i]  -(BOX*lroundf(X[3*i]/BOX)));
		X[3*i+1]=(X[3*i+1]-(BOX*lroundf(X[3*i+1]/BOX)));
		X[3*i+2]=(X[3*i+2]-(BOX*lroundf(X[3*i+2]/BOX)));

		input<<X[3*i]<<"\t"<<X[3*i+1]<<"\t"<<X[3*i+2]<<"\t"<<RAD[i]<<"\t"<<V[3*i]<<"\t"<<V[3*i+1]<<"\t"<<V[3*i+2]<<"\n";
		//pressure_virial=pressure_virial-(V[2*i]*V[2*i]+V[2*i+1]*V[2*i+1]);
	}
	fstream stat("stat_"+std::string(name)+".dat", fstream::out | fstream::app);
	double Vol=1.;
	for( int d=0; d< dim ; d++)
	{
		Vol=Vol*BOX;
	}
	//cout<<k<<"\t"<<"\t"<<PEnergy<<"\t"<<KEnergy<<"\t"<<KEnergy/(N)<<"\t"<<T<<"\n";I/
	//cout<<k<<"\t"<<tilt/BOX<<"\t"<<PEnergy<<"\t"<<KEnergy<<"\t"<<(2./3.)*KEnergy/(N)<<"\t"<<T<<"\t"<<stress_tensor[0][0]<<"\t"<<stress_tensor[0][1]<<"\t"<<stress_tensor[1][1]<<"\t"<<pressure_virial<<"\t"<<p_applied<<"\t"<<3.14159265359*(N/2.)*(0.5*0.5+0.7*0.7)/(BOX*BOX)<<"\n";;
	//cout<<k<<"\t"<<tilt/BOX<<"\t"<<PEnergy<<"\t"<<KEnergy<<"\t"<<2./3.*KEnergy/(N)<<"\t"<<T<<"\t"<<stress_tensor[0][0]<<"\t"<<stress_tensor[0][1]<<"\t"<<stress_tensor[1][1]<<"\t"<<(pressure_virial)/(2./3.*KEnergy/(N))*(N/(BOX*BOX*BOX))<<"\t"<<pressure_virial2<<"\t"<<pressure_virial<<"\t"<<p_applied<<"\t"<<(4./3.)*3.14159265359*(N/2.)*(0.5*0.5*0.5+0.7*0.7*0.7)/(BOX*BOX*BOX)<<"\n";;
	stat<<k<<"\t"<<tilt/BOX<<"\t"<<PEnergy<<"\t"<<KEnergy<<"\t"<<2./3.*KEnergy/(N)<<"\t"<<T<<"\t"<<stress_tensor[0][0]<<"\t"<<stress_tensor[0][1]<<"\t"<<stress_tensor[1][1]<<"\t"<<(pressure_virial)/((2./3.*KEnergy/(N))*(N/(BOX*BOX*BOX)))<<"\t"<<-1.*pressure_virial2<<"\t"<<pressure_virial<<"\t"<<p_applied<<"\t"<<(4./3.)*3.14159265359*(N/2.)*(0.5*0.5*0.5+0.7*0.7*0.7)/(BOX*BOX*BOX)<<"\n";;
	
}
