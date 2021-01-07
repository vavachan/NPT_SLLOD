#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <iomanip>

using namespace std;

// double BOX=2*7.2005488787721159;
// double BOX=2*26.2107109945647245;
double BOX;
int N=1000;
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

//double *ATOMx = new  double[N];
//double *ATOMy = new  double[N];
// double *ATOMz = new  double[N];
double *RAD   = new  double[N];
//static  double *DIST= new  double[N];
static int *VLIST= new int [N*(N+1)];

double stress_tensor[2][2];
double pressure_virial2=0.;
double pressure_virial1=0.;
double pressure_virial=0.;

//static  double *tATOMx = new  double[N];
//static  double *tATOMy = new  double[N];
//static  double *tATOMz = new  double[N];
// gradient is calculated and stored in this 

//double *CONJ_D= new  double[2*N+1];
static double *prevG= new  double[2*N+1];
double *G= new  double[2*N+1];
double *X= new  double[2*N+1];
double *V= new  double[2*N+1];
static double *tX= new  double[2*N+1];

double deltadrv=0.5;
double deltadrvsq=deltadrv*deltadrv;


void make_list()
{
	for(int i=0; i<2*N; i++)
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
		//ATOMx[i]=X[i*2];
		//ATOMy[i]=X[i*2+1];
		GRAD[i*2] = 0.;
		GRAD[i*2+1] = 0.;
//		ATOMz[i]=X[i*3+2];

		delta_rx=(tX[2*i]-X[i*2]);
		delta_ry=(tX[2*i+1]-X[i*2+1]);
		//cout<<i<<"\t"<<delta_rx<<"\n";
//		delta_rz=tATOMz[i]-ATOMz[i];
		    
		if(tilt)
	    	delta_rx=(delta_rx-(tilt*lroundl(delta_ry*iBOX)));
	    delta_rx=(delta_rx-(BOX*lroundl(delta_rx*iBOX)));
  	    delta_ry=(delta_ry-(BOX*lroundl(delta_ry*iBOX)));
  //	    delta_rz=(delta_rz-(BOX*lroundl(delta_rz/BOX)));
		
		delta_r=delta_rx*delta_rx+delta_ry*delta_ry;//+delta_rz*delta_rz;
		
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
					delta_rx=(X[2*i]-X[2*j]);
					delta_ry=(X[2*i+1]-X[2*j+1]);
					//delta_rz=ATOMz[i]-ATOMz[j];

					if(tilt)
						delta_rx=(delta_rx-(tilt*lroundl(delta_ry*iBOX)));
					delta_rx=(delta_rx-(BOX*lroundl(delta_rx*iBOX)));
					delta_ry=(delta_ry-(BOX*lroundl(delta_ry*iBOX)));
					//delta_rz=(delta_rz-(BOX*lroundl(delta_rz/BOX)));

					delta_r=delta_rx*delta_rx+delta_ry*delta_ry;//+delta_rz*delta_rz;
					
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
		for(int i=0; i<2*N; i++)
		{
			tX[i]=X[i];
		}
	}
	nlistbeg=0;
	stress_tensor[0][0]=0.;
	stress_tensor[0][1]=0.;
	stress_tensor[1][1]=0.;
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
				delta_rx=(X[2*i]-X[2*j]);
				delta_ry=(X[2*i+1]-X[2*j+1]);
				//delta_rx=ATOMx[i]-ATOMx[j];
				//delta_ry=ATOMy[i]-ATOMy[j];
//				delta_rz=ATOMz[i]-ATOMz[j];

				if(tilt)
					delta_rx=(delta_rx-(tilt*lroundl(delta_ry*iBOX)));
				delta_rx=(delta_rx-(BOX*lroundl(delta_rx*iBOX)));
				delta_ry=(delta_ry-(BOX*lroundl(delta_ry*iBOX)));
//					delta_rz=(delta_rz-(BOX*lroundl(delta_rz/BOX)));

				delta_r=delta_rx*delta_rx+delta_ry*delta_ry;//+delta_rz*delta_rz;
				
				//delta_r=pow(delta_r,0.5);
				double sigma=RAD[i]+RAD[j];
				double sigmasq=sigma*sigma;
				//cout<<i<<"\t"<<j<<"\t"<<delta_r<<"\n";
				if(delta_r<sigmasq)
				{
					//pot_e=pot_e+U_ij(delta_r,RAD[i]+RAD[j]);
					double isigma=1./sigma;
					delta_r=sqrt(delta_r);
					double vh=(1-delta_r*isigma);
					pot_e=pot_e+epsilon*vh*vh;
					//dwreturn epsilon*2.*(1-r/sigma)*(-1./sigma);
					double U=epsilon*2.*vh*(-1.*isigma)/delta_r;

					GRAD[i*2]   = GRAD[i*2]   + U*delta_rx;
					GRAD[i*2+1] = GRAD[i*2+1] + U*delta_ry;

					GRAD[j*2]   = GRAD[j*2]   - U*delta_rx;
					GRAD[j*2+1] = GRAD[j*2+1] - U*delta_ry;

					stress_tensor[0][0]=stress_tensor[0][0]+(U*delta_rx)*delta_rx;
					stress_tensor[0][1]=stress_tensor[0][1]+(U*delta_rx)*delta_ry;
					stress_tensor[1][1]=stress_tensor[1][1]+(U*delta_ry)*delta_ry;
					
					pressure_virial2 = pressure_virial2+(U*delta_rx)*delta_rx+(U*delta_ry)*delta_ry;
				}
			}
		}
		nlistbeg=jj;//nlistbeg+jj+1;
	}
	stress_tensor[0][0]=stress_tensor[0][0]/N;
	stress_tensor[0][1]=stress_tensor[0][1]/N;
	stress_tensor[1][1]=stress_tensor[1][1]/N;
	clock_t end=clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	//cout<<time_spent<<" in energy \n";
////for(int i=0; i<N; i++)
////{
////	GRAD[2*i]=GRAD[2*i]*BOX;
////	GRAD[2*i+1]=GRAD[2*i+1]*BOX;
////	GRAD[ndim-1]=GRAD[ndim-1]+(GRAD[2*i]/BOX*ATOMx[i]+GRAD[2*i+1]/BOX*ATOMy[i]);
////}
////GRAD[ndim-1]=GRAD[ndim-1]/BOX+2*BOX*p_applied;
	return pot_e;
}
void velocity_verlet(double X[], double GRAD[], double delta_t)
{
	double delta_t_sq=delta_t*delta_t;
	double linear_momentum=0.;
	double linear_momentum_x=0.;
	double linear_momentum_y=0.;
	//PEnergy=energy_force(X,GRAD,delta_t);
	cout<<std::setprecision(16);
	//cout<<KEnergy+PEnergy<<"\t"<<KEnergy<<"\t"<<PEnergy<<"\n";
	for(int i=0; i <N; i++)
	{
		//velocity update
		//cout<<V[2*i]<<"\n";
		V[2*i]=V[2*i]-0.5*delta_t*(GRAD[2*i]);
		V[2*i+1]=V[2*i+1]-0.5*delta_t*(GRAD[2*i+1]);

		//poition update
		X[2*i]=X[2*i]+delta_t*V[2*i];
		X[2*i+1]=X[2*i+1]+delta_t*V[2*i+1];

		
	}
	PEnergy=energy_force(X,GRAD,2*N);
	//cout<<PEnergy<<"\n";
	KEnergy=0.;
	for(int i=0; i <N; i++)
	{
		V[2*i]=V[2*i]-0.5*delta_t*(GRAD[2*i]);
		V[2*i+1]=V[2*i+1]-0.5*delta_t*(GRAD[2*i+1]);
		
		KEnergy=KEnergy+0.5*(V[2*i]*V[2*i]+V[2*i+1]*V[2*i+1]);
		//cout<<KEnergy<<"\n";
		linear_momentum_x=linear_momentum_x+V[2*i];
		linear_momentum_y=linear_momentum_y+V[2*i+1];
	}
	linear_momentum=linear_momentum_x*linear_momentum_x+linear_momentum_y*linear_momentum_y;
	
	cout<<KEnergy+PEnergy<<"\t"<<KEnergy<<"\t"<<PEnergy<<"\t"<<linear_momentum<<"\t"<<KEnergy/(N)<<"\n";
		//prevG[2*i]=GRAD[2*i];
		//prevG[2*i+1]=GRAD[2*i+1];
		

} 

void NPT_mVV( double X[], double GRAD[], double delta_t, double v_e, double shear_rate, double STRAIN)
{
	double delta_t4=delta_t/4.;
	double delta_t2=delta_t/2.;
	double delta_t8=delta_t/8.;
	
	for( int i=0; i < N; i++)
	{
		V[2*i]  =V[2*i]  +delta_t2*(-1.*GRAD[2*i]  );
		V[2*i+1]=V[2*i+1]+delta_t2*(-1.*GRAD[2*i+1]);
	}
	double scaley1=exp(v_e*delta_t2);;
	double scaley2=exp(v_e*delta_t4);;
	double sin_hypy=sinh(delta_t4*v_e);

	double scalex1=exp(v_e*delta_t);;
	double scalex2=exp(v_e*delta_t2);;
	double sin_hypx=sinh(delta_t2*v_e);

	for( int i=0; i < N; i++)
	{
		//X[2*i+1] = X[2*i+1]*scaley1 + 4.*(V[2*i+1])*scaley2*sin_hypy;
		X[2*i+1] = X[2*i+1]*scaley1 + 2.*(V[2*i+1]/v_e)*scaley2*sin_hypy;
		X[2*i]   = X[2*i]*scalex1   + 2.*((V[2*i]+shear_rate*X[2*i+1])/v_e)*scalex2*sin_hypx;
		//X[2*i+1] = X[2*i+1]*scaley1 + 2.*V[2*i+1]*scaley2*sin_hypy;
		X[2*i+1] = X[2*i+1]*scaley1 + 2.*(V[2*i+1]/v_e)*scaley2*sin_hypy;
		//X[2*i+1] = X[2*i+1]*scaley1 + 2.*V[2*i+1]*scaley2*sin_hypy;
	}

	BOX=BOX*exp(delta_t*v_e);

	tilt=STRAIN*BOX;
	PEnergy=energy_force(X,GRAD,2*N);
	KEnergy=0.;
	
	pressure_virial1=0.;
	for(int i=0; i <N; i++)
	{
		V[2*i]  =V[2*i]  +delta_t2*(-1.*GRAD[2*i]  );
		V[2*i+1]=V[2*i+1]+delta_t2*(-1.*GRAD[2*i+1]);
		pressure_virial1=pressure_virial1+(V[2*i]*V[2*i]+V[2*i+1]*V[2*i+1]);
		KEnergy=KEnergy+0.5*(V[2*i]*V[2*i]+V[2*i+1]*V[2*i+1]);
	}
	pressure_virial=(pressure_virial1-pressure_virial2)/(2*BOX*BOX);
		
}


void NH_NPT_SLLOD(double X[], double GRAD[], double delta_t, double &eta, double &eta_e, double &v_eta, double &v_e , double &v_eta_e, double W, double &KE, double shear_rate)
{

	double L=2*N;
	double G_eta_e=W*v_eta*v_eta-T;
	double G=2.*KE-L*T;
	double Vol=BOX*BOX;
	double Ge=2.*Vol*(pressure_virial-p_applied)+(2./L)*(2.*KE);

	double delta_t4=delta_t/4.;
	double delta_t2=delta_t/2.;
	double delta_t8=delta_t/8.;

	v_eta = v_eta + delta_t4*G;

	v_eta_e = v_eta_e + delta_t4*G_eta_e;

	v_e = v_e*exp(-1.*delta_t8*v_eta_e);
	v_e = v_e + delta_t4*Ge;
	v_e = v_e*exp(-1.*delta_t8*v_eta_e);
	
	eta = eta + delta_t2*v_eta;
	eta_e = eta_e +delta_t2*v_eta_e;

	double scalex= exp(-1.*delta_t8*(v_eta +(1+2./L)*v_e));
	double scaley= exp(-1.*delta_t2*(v_eta +(1+2./L)*v_e));
	pressure_virial1=0.;
	KE=0.;
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
	pressure_virial=(pressure_virial1-pressure_virial2)/(2*BOX*BOX);

	Ge=2.*Vol*(pressure_virial-p_applied)+(2./L)*(2.*KE);
	
	v_e = v_e*exp(-1.*delta_t8*v_eta_e);
	v_e = v_e + delta_t4*Ge;
	v_e = v_e*exp(-1.*delta_t8*v_eta_e);
	
	G=2.*KE-L*T;
	G_eta_e=W*v_eta*v_eta-T;
	
	
	v_eta_e = v_eta_e + delta_t4*G_eta_e;

	v_eta = v_eta + delta_t4*G;
	
}

void integrate_NPT_SLLOD(double X[], double GRAD[], double delta_t)
{
	static double eta=0.;
	static double eta_e=0.;
	static double v_eta=0.;
	static double v_eta_e=0.;
	static double v_e=0.;
	double shear_rate=0.0001;
	double W=10.;
	STRAIN=STRAIN+shear_rate*delta_t;
	NH_NPT_SLLOD(X,GRAD,delta_t, eta, eta_e, v_eta, v_e, v_eta_e, W, KEnergy, shear_rate);
	NPT_mVV(X,GRAD,delta_t,v_e,shear_rate,STRAIN);
	NH_NPT_SLLOD(X,GRAD,delta_t, eta, eta_e, v_eta, v_e, v_eta_e, W, KEnergy, shear_rate);
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
			X[2*i]=(X[2*i]-(tilt*lroundf(X[2*i+1]/BOX)));
		X[2*i]  =(X[2*i]  -(BOX*lroundf(X[2*i]/BOX)));
		X[2*i+1]=(X[2*i+1]-(BOX*lroundf(X[2*i+1]/BOX)));

		input<<X[2*i]<<"\t"<<X[2*i+1]<<"\t"<<RAD[i]<<"\t"<<V[2*i]<<"\t"<<V[2*i+1]<<"\n";
		//pressure_virial=pressure_virial-(V[2*i]*V[2*i]+V[2*i+1]*V[2*i+1]);
	}
	fstream stat("stat_"+std::string(name)+".dat", fstream::out | fstream::app);
	//cout<<k<<"\t"<<"\t"<<PEnergy<<"\t"<<KEnergy<<"\t"<<KEnergy/(N)<<"\t"<<T<<"\n";
	stat<<k<<"\t"<<"\t"<<PEnergy<<"\t"<<KEnergy<<"\t"<<KEnergy/(N)<<"\t"<<T<<"\t"<<stress_tensor[0][0]<<"\t"<<stress_tensor[0][1]<<"\t"<<stress_tensor[1][1]<<"\t"<<pressure_virial<<"\t"<<p_applied<<"\t"<<3.14159265359*(N/2.)*(0.5*0.5+0.7*0.7)/(BOX*BOX)<<"\n";;
}
