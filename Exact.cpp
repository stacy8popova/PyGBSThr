#include "headers.h" 
#include "matrix.h"

const char* in_name="in_out/Submatrix.dat";


int N=sqrt(file_size(in_name)/2);



bool const calc_sectors=true;  // choose "false" if you need the final value of probability, without calculating sectors, this works faster

int Nnu2=10*N; long double dnu2=2*Pi/Nnu2;

complex * stat = new complex [N+1];
ifstream m_str(in_name); ofstream ou("in_out/Exact.dat");


complex ** M=new_complex2(N,N), ** M0=new_complex2(N,N), ** A=new_complex2(N,N);;

complex ** Z=new_complex2(N+1, Nnu2), ** ZN=new_complex2(N+1, Nnu2);;



int n_ph(int i)
{
    int s=0; for1 (n, N) {s+=i%2; i/=2;} return s;
}

long double factorial(int x)
{
    long double s=1; for1 (i, x) s*=(i+1.); return s;
}




int main()

{   
    ou<<setprecision(20);
    for1(n,N+1) {stat[n]=0; for1 (j, Nnu2) {Z[n][j]=0;ZN[n][j]=0;}}
  
  
  for2(n1, n2, N) {long double x, y; m_str>>x; m_str>>y; M[n1][n2]=x+I*y;} //cout<<M[0][0]<<"\n";
     
    
    
   // for2(n1,n2,N)    {        M[n1][n2]=0;        for1(n_, N/2) M[n1][n2]+=M0[n1][n_]*M0[n2][n_]*.44;}
  


    

    for1(i, 1<<N)
    {
        //cout<<i<<"  "<<flush;
        int n1_=0, i1_=i; for1 (n1, N)
        {
            if (i1_%2==1)
            {
                int n2_=0, i2_=i; for1 (n2, N)
                {
                    if (i2_%2==1)
                    {
                        M0[n1_][n2_]=M[n1][n2];
                        n2_++;
                    }
                    i2_/=2;
                }
                n1_++;
            }
            i1_/=2; 
        }
                
        
        int n_ph_i=n_ph(i); static int n_ph_old=0; if (n_ph_i>n_ph_old) // {cout<<n_ph_i<<" "<<flush; n_ph_old=n_ph_i;}
        for2(n1,n2, n_ph_i) A[n1][n2]=0.; 
        for1(n,n_ph_i) A[n][n]=0;
        for1(n1, n_ph_i) for (int n2=n1; n2<n_ph_i; n2++) {for1(n_, n_ph_i) A[n1][n2]+=conj(M0[n1][n_])*M0[n_][n2]; A[n2][n1]=conj(A[n1][n2]);}
        
        if (n_ph_i==0) for1(j, Nnu2) Z[i][j]=1;
        else
        {
            static long double * lambda = new long double [N];
            if (calc_sectors) EigenJacobi(A,M0,lambda, n_ph_i);
            //if (n_ph_i==N) {for1(n, N) cout<<4.*lambda[n]<<"  "; cout<<"\n";}
            static complex * e4_nu2=new complex [Nnu2]; do_once  for1(j, Nnu2) e4_nu2[j]=4.l*exp(I*(j*dnu2)); end_do_once;
            if (calc_sectors) for1(j, Nnu2)
            {
                //long double nu2=dnu2*j;
                complex r=1.; for1(n, n_ph_i) r/=sqrt(1.l-e4_nu2[j]*lambda[n]); Z[n_ph_i][j]+=r; //if (n_ph_i==N) ou_debug<<nu2<<"  "<<real(r)<<"  "<<imag(r)<<"\n";
                
            }
        }
        
        
        
        for2(n1,n2, n_ph_i) A[n1][n2]*=-4.l; for1(n, n_ph_i) A[n][n]+=1.l;    
        
        if (n_ph_i==0) stat[n_ph_i]+=1.; else stat[n_ph_i]+=1.l/sqrt(det(A,n_ph_i));  
    }
    
    //cout<<"\n";

    //cout<<stat[N]<<" = "<<Z[N][0]<<"\n";
    
    
    
    for1(k, N)
    {
        for(int n=k+1; n<=N; n++) {stat[n]-=stat[k]*factorial(N-k)/(factorial(n-k)*factorial(N-n)); for1(j,Nnu2) Z[n][j]-=Z[k][j]*factorial(N-k)/(factorial(n-k)*factorial(N-n));}
    }
    
    
    //cout<<stat[N]<<" = "<<Z[N][0]<<"\n";
    
    
    for1(n, N+1) for1(j, Nnu2)
    {
        for1(k, Nnu2) ZN[n][j]+=Z[n][k]*exp(-I*(j*k*dnu2))/(1.l*Nnu2);
    }
    
    long double s=0; for1(n, N+1) s+=real(stat[n]);
    for1(n, N+1) {ou<<n<<"  "<<real(stat[n])/s; if (calc_sectors) for1(j,  Nnu2) {ou<<"  "<<real(ZN[n][j])/s;} ou<<"\n"<<flush; }
     
    
    
    
//    for1(n1, n_ph(i)) {for1(n2, n_ph(i)) cout<<M0[n1][n2]<<"  "; cout<<"\n"; }
    
    return 0;
}
