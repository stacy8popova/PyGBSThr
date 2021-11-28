#include <iostream>
#include <fstream>
#include <complex>
#include <iomanip>
#include <time.h>
#include <stdlib.h>


#define complex complex<long double>
using namespace std;

#define do_once  {static int flag_once=0; if (flag_once==0) {flag_once=1;
#define end_do_once ;};}

int sgn(long double x) {if (x>0) return 1; if (x<0) return -1; return 0;}
long double sqr(long double x) {return x*x;}
//long double sqr(long double x) {return x*x;}

complex sqr(const complex & x) {return x*x;}
long double to0(long double r) {if (abs(r)>1e-10) return r; else return 0;}


complex I(0,1); long double Pi=4.l*atan(1.l);
inline long double norm2(complex x) {long double a=real(x), b=imag(x); return a*a+b*b;}


int initial_random=time(NULL);
int INT_RANDOM=initial_random;//  //from 1 to 2^31-1

int int_rnd()
{
//for (int j=1; j<3; j++)
{
	 int k=INT_RANDOM/127773;
	 INT_RANDOM=16807*(INT_RANDOM-k*127773)-2836*k;
	 if (INT_RANDOM<0) INT_RANDOM+=2147483647;
;}
	 return INT_RANDOM;
;}

int INT_RANDOM_AUX_VARIABLE=int_rnd()+int_rnd()+int_rnd()+int_rnd();  //randomising...

long double rnd ()
{
	return int_rnd()/2147483647.0
;}

int rnd (int k)
{
	int d=2147483647%k, d1=2147483647-d, r=int_rnd();
	if (r>=d1) return rnd(k);
	return r/(d1/k)
;}

complex rnd_gauss2()
{
  long double R=sqrt(-log(rnd())), phi=2*Pi*rnd(); return R*exp(I*phi);
  
//Gauss::   long double R=sqrt(-2.*log(rnd())), phi=2*Pi*rnd(); return R*exp(I*phi);
  
}

void swing(int & x, int & y) {int a=x; x=y; y=a;}
void swing(long double & x, long double & y) {long double a=x; x=y; y=a;}
void swing(complex & x, complex & y) {complex a=x; x=y; y=a;}



int min(int & x, int &y) {if (x<y) return x; else return y;}
int max(int & x, int &y) {if (x>y) return x; else return y;}

long double max_abs(const long double & x, const long double &y) {if (abs(x)>abs(y)) return x; else return y;}

int file_size(const char * file_name) {ifstream ff(file_name); int i=-1; long double x; do {ff>>x; i++;} while(!ff.eof()); return i; }

class function_of_time
{
  long double ** x, * x0; int nt;
public:
  function_of_time(const char * filename, int nv)
    {
      nt=0; {ifstream fs(filename); long double r; for (; !fs.eof(); nt++) fs>>r;} nt=(nt-1)/(4*nv+1)-2; //cout<<"nt="<<nt<<"\n"<<flush;
      x=new long double * [nt+1]; ifstream fs(filename); x0=new long double[4*nv+1]; for (int j=0; j<4*nv+1; j++) fs>>x0[j];
      for (int i=0; i<=nt; i++) {x[i]=new long double [4*nv+1]; for (int j=0; j<4*nv+1; j++) fs>>x[i][j];}
    }
  ~function_of_time() {for (int i=0; i<=nt; i++) delete [] x[i]; delete [] x;}
  complex value(long double t, int v)
  {
    if (t<0) return complex(x0[4*v+1],x0[4*v+2]);
    long double dt=(x[nt][0]-x[0][0])/nt;
    int i=int((t-x[0][0])/dt); if (i>=nt) i=nt-1; if (i<0) i=0;
    long double t1=t-x[0][0]-i*dt;
    complex a0(x[i][1+4*v], x[i][2+4*v]), a1(x[i][3+4*v], x[i][4+4*v]),b0(x[i+1][1+4*v], x[i+1][2+4*v]), b1(x[i+1][3+4*v], x[i+1][4+4*v]);
    complex a2=(3.l*(b0-a0)-(b1+2.l*a1)*dt)/(dt*dt), a3=(2.l*(a0-b0)+(b1+a1)*dt)/(dt*dt*dt);
    return a0+(a1+(a2+a3*t1)*t1)*t1;
  }
};


int ** new_int2(int n1, int n2)
{
  int ** r=new int *[n1];
  for (int i=0; i<n1; i++)  r[i]=new int [n2];   
  return r;
}

void delete_int2(int ** &r, int n1, int n2)
{
  for (int i=0; i<n1; i++)  delete [] r[i];   
  delete [] r;
}


long double ** new_double2(int n1, int n2)
{
  long double ** r=new long double *[n1];
  for (int i=0; i<n1; i++)  r[i]=new long double [n2];   
  return r;
}

void delete_double2(long double ** &r, int n1, int n2)
{
  for (int i=0; i<n1; i++)  delete [] r[i];   
  delete [] r;
}

long double ** new_ldouble2(int n1, int n2)
{
  long double ** r=new long double *[n1];
  for (int i=0; i<n1; i++)  r[i]=new long double [n2];   
  return r;
}

void delete_ldouble2(long double ** &r, int n1, int n2)
{
  for (int i=0; i<n1; i++)  delete [] r[i];   
  delete [] r;
}


complex ** new_complex2(int n1, int n2)
{
  complex ** r=new complex *[n1];
  for (int i=0; i<n1; i++)  r[i]=new complex [n2];   
  return r;
}

void delete_complex2(complex ** &r, int n1, int n2)
{
  for (int i=0; i<n1; i++)  delete [] r[i];   
  delete [] r;
}


int ***  new_int3(int n1, int n2, int n3)
{
  int *** r;
  r=new int **[n1];
  for (int i=0; i<n1; i++) 
  {
    r[i]=new int *[n2];
    for (int j=0; j<n2; j++) r[i][j]=new int [n3];
  }
  return r;
}

void delete_int3(int *** &r, int n1, int n2, int n3)
{
  for (int i=0; i<n1; i++)  
  {
    for (int j=0; j<n2; j++) delete [] r[i][j];
    delete [] r[i];  
  }
  delete [] r;
}



complex *** new_complex3(int n1, int n2, int n3)
{
  complex *** r=new complex **[n1];
  for (int i=0; i<n1; i++) 
  {
    r[i]=new complex *[n2];
    for (int j=0; j<n2; j++) r[i][j]=new complex [n3];
  }
  return r;
}

void delete_complex3(complex *** &r, int n1, int n2, int n3)
{
  for (int i=0; i<n1; i++)  
  {
    for (int j=0; j<n2; j++) delete [] r[i][j];
    delete [] r[i];  
  }
  delete [] r;
}


int **** new_int4(int n1, int n2, int n3, int n4)
{
  int **** r=new int ***[n1];
  for (int i=0; i<n1; i++) 
  {
    r[i]=new int **[n2];
    for (int j=0; j<n2; j++) 
    {
      r[i][j]=new int * [n3];
      for (int l=0; l<n3; l++) r[i][j][l]=new int [n4];
    }
  }
  return r;
}


complex **** new_complex4(int n1, int n2, int n3, int n4)
{
  complex **** r=new complex ***[n1];
  for (int i=0; i<n1; i++) 
  {
    r[i]=new complex **[n2];
    for (int j=0; j<n2; j++) 
    {
      r[i][j]=new complex * [n3];
      for (int l=0; l<n3; l++) r[i][j][l]=new complex [n4];
    }
  }
  return r;
}


#define for1(i, N) for (int i=0; i<N; i++) 
#define for2(i, j, N) for (int i=0; i<N; i++) for (int j=0; j<N; j++)
#define for3(i, j, k, N) for (int i=0; i<N; i++) for (int j=0; j<N; j++) for (int k=0; k<N; k++)

