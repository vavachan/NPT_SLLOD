extern int N;
extern  double BOX;
extern  double epsilon;
extern  double magnitude;
extern  double p_applied;
extern double T;

//extern  double *tATOMz;

extern  double *RAD  ;
extern  double *CONJ_D;
extern  double *G;
extern  double *X;
extern  double *V;
extern  double delta ;

extern double pressure_virial;
extern double pressure_virial1;
extern double pressure_virial2;

//extern double PEnergy;
//extern double KEnergy;



double U_ij( double , double );

double der_U_ij( double , double );

double energy( double *,int );

double energy_force( double *,  double *, int);
void calculate_gradient( double *,  double *, int);

void write_config(int ,  double *, char *);

void make_list();
void integrate_NPT_SLLOD(double *, double *, double );
void integrate_NVT_SLLOD(double *, double *, double );
void integrate_NH(double *,  double *,double);
void velocity_verlet(double *,  double *,double);
void velocity_verlet_sllod(double *,double *,double);
void velocity_verlet_Gaussian(double *,double *,double);
void velocity_verlet_Nose_Hoover(double *,double *,double);
