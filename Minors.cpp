#include "headers.h"
#include "matrix.h"
const char* in_name="in_out/Submatrix.dat";
bool enable_checks=false;

int N=sqrt(file_size(in_name)/2), Nnu=N*10; ifstream a_str(in_name); ofstream ou1("in_out/Minors0-1.dat"), ou2 ("in_out/Minors2.dat"), ou3 ("in_out/Minors3.dat"), ou4 ("in_out/Minors4.dat");

complex ** M=new_complex2(N,N), ** MM=new_complex2(N,N), ** U=new_complex2(N,N), ** MU=new_complex2(N,N), ** MMU=new_complex2(N,N);
long double * lambda=new long double [N];
complex ** A=new_complex2(N,N), ** MA=new_complex2(N,N), ** MMA=new_complex2(N,N), ** MAM=new_complex2(N,N), ** MMAM=new_complex2(N,N), ** MMAMM=new_complex2(N,N), ** AM=new_complex2(N,N), ** AMM=new_complex2(N,N), ** MAMM=new_complex2(N,N);

complex sqrt_det(complex ** & b, int N)
{
    
     
 	complex s=1;  
    
	for (int i=0; i<N-1; i++)
	{
		if (norm2(b[i][i])<=0) return 0;
		for (int j=i+1; j<N; j++)
		{
		complex f=b[j][i]/b[i][i];
		for (int k=i; k<N; k++) b[j][k]-=f*b[i][k];
		;}
	;}

 	
   if (abs(b[N-1][N-1])<=0) return 0;
   {for (int i=0; i<N; i++) s*=sqrt(b[i][i]);}
   return s;
;}


void makeMMU()
{
    
    for2(i,j, N){MM[i][j]=0; for1(l, N) MM[i][j]+=conj(M[l][i])*M[l][j];}
    EigenJacobi(MM,U,lambda, N); for1(i, N) for1(j, i) swing(U[i][j], U[j][i]);
    
    for2(i,j,N) {MMU[i][j]=0; for1(l,N) MMU[i][j]+=MM[i][l]*U[l][j];}
    for2(i,j,N) {MU[i][j]=0; for1(l,N) MU[i][j]+=M[i][l]*U[l][j];}    
}


complex make_products(complex r)
{
    complex d=1;
    for2(i,j, N) A[i][j]=0; for1 (i, N) {A[i][i]=1.l/(1.l-r*lambda[i]); d/=sqrt(1.l-r*lambda[i]);}
    
    for2(i,j, N) {MA[i][j]=0; for1(l,N) MA[i][j]+=MU[i][l]*A[l][l]*conj(U[j][l]);}
    for2(i,j, N) {MAM[i][j]=0; for1(l,N) MAM[i][j]+=MU[i][l]*A[l][l]*conj(MU[j][l]);}
    
    for2(i,j, N) {MMA[i][j]=0; for1(l,N) MMA[i][j]+=MMU[i][l]*A[l][l]*conj(U[j][l]);}
    for2(i,j, N) {MMAM[i][j]=0; for1(l,N) MMAM[i][j]+=MMU[i][l]*A[l][l]*conj(MU[j][l]);}
    for2(i,j, N) {MMAMM[i][j]=0; for1(l,N) MMAMM[i][j]+=MMU[i][l]*A[l][l]*conj(MMU[j][l]);}
    for2(i,j, N) {MA[i][j]=0; for1(l,N) MA[i][j]+=MU[i][l]*A[l][l]*conj(U[j][l]);}
    for2(i,j, N) {AM[i][j]=0; for1(l,N) AM[i][j]+=U[i][l]*A[l][l]*conj(MU[j][l]);}
    for2(i,j, N) {AMM[i][j]=0; for1(l,N) AMM[i][j]+=U[i][l]*A[l][l]*conj(MMU[j][l]);}    
    for2(i,j, N) {MAMM[i][j]=0; for1(l,N) MAMM[i][j]+=MU[i][l]*A[l][l]*conj(MMU[j][l]);}
    

    for2(i,j, N) {A[i][j]=0; for1(l,N) A[i][j]+=U[i][l]*(1.l/(1.l-r*lambda[l]))*conj(U[j][l]);}

    
    return d;
}

complex rat_Z(int h, complex r)
{
    static complex **b=new_complex2(3,3);
     
    b[0][0]=A[h][h]; b[1][0]=MA[h][h]; b[2][0]=MMA[h][h];
    b[1][1]=MAM[h][h]; b[2][1]=MMAM[h][h]; b[2][2]=MMAMM[h][h];
    b[0][1]=AM[h][h];b[0][2]=AMM[h][h];b[1][2]=MAMM[h][h];
    

    static complex ** a=new_complex2(3,3); for2(i,j,3) a[i][j]=0;
    a[0][0]=MM[h][h]-conj(M[h][h])*M[h][h]; 
    a[0][1]=conj(M[h][h]); a[1][0]=M[h][h]; 
    a[0][2]=-1; a[1][1]=-1; a[2][0]=-1; 
 

    
    static complex **ab1=new_complex2(3,3);
    for2(k,l,3)
    {
        ab1[k][l]=0;
        for1(i, 3) ab1[k][l]-=r*a[k][i]*b[i][l];
        if (k==l) ab1[k][l]+=1;
    }
    
    return 1.l/sqrt_det(ab1,3); 
}


complex rat_Z(int h1, int h2, complex r)
{
    static complex **b=new_complex2(6,6); for2(i,j,6) b[i][j]=0;
    
 
    b[0][0]=A[h1][h1]; b[0][1]=A[h1][h2]; b[1][0]=A[h2][h1]; b[1][1]=A[h2][h2];
    b[2][0]=MA[h1][h1];b[2][1]=MA[h1][h2];b[3][0]=MA[h2][h1];b[3][1]=MA[h2][h2];
    b[4][0]=MMA[h1][h1];b[4][1]=MMA[h1][h2];b[5][0]=MMA[h2][h1];b[5][1]=MMA[h2][h2];
    b[2][2]=MAM[h1][h1];b[2][3]=MAM[h1][h2];b[3][2]=MAM[h2][h1];b[3][3]=MAM[h2][h2];
    b[4][2]=MMAM[h1][h1];b[4][3]=MMAM[h1][h2];b[5][2]=MMAM[h2][h1];b[5][3]=MMAM[h2][h2];
    b[4][4]=MMAMM[h1][h1];b[4][5]=MMAMM[h1][h2];b[5][4]=MMAMM[h2][h1];b[5][5]=MMAMM[h2][h2];
    b[0][2]=AM[h1][h1];b[0][3]=AM[h1][h2];b[1][2]=AM[h2][h1];b[1][3]=AM[h2][h2];
    b[0][4]=AMM[h1][h1];b[0][5]=AMM[h1][h2];b[1][4]=AMM[h2][h1];b[1][5]=AMM[h2][h2];
    b[2][4]=MAMM[h1][h1];b[2][5]=MAMM[h1][h2];b[3][4]=MAMM[h2][h1];b[3][5]=MAMM[h2][h2];
    
    static complex **a=new_complex2(6,6); for2(i,j,6) a[i][j]=0;
    
    a[0][0]=MM[h1][h1]-conj(M[h1][h1])*M[h1][h1]-conj(M[h1][h2])*M[h2][h1];
    a[0][1]=MM[h1][h2]-conj(M[h1][h1])*M[h1][h2]-conj(M[h1][h2])*M[h2][h2];
    a[1][0]=MM[h2][h1]-conj(M[h2][h1])*M[h1][h1]-conj(M[h2][h2])*M[h2][h1];
    a[1][1]=MM[h2][h2]-conj(M[h2][h1])*M[h1][h2]-conj(M[h2][h2])*M[h2][h2];    
    a[0][2]=conj(M[h1][h1]); a[1][2]=conj(M[h1][h2]); a[0][3]=conj(M[h2][h1]);a[1][3]=conj(M[h2][h2]); a[3][1]=M[h2][h2]; a[3][0]=M[h2][h1]; a[2][1]=M[h1][h2]; a[2][0]=M[h1][h1]; 
    a[0][4]=-1; a[1][5]=-1; a[2][2]=-1; a[3][3]=-1; a[5][1]=-1;a[4][0]=-1;
    
    
    static complex **ab1=new_complex2(6,6);
    for2(k,l,6)
    {
        ab1[k][l]=0;
        for1(i, 6) ab1[k][l]-=r*a[k][i]*b[i][l];
        if (k==l) ab1[k][l]+=1;
    }
    
    return 1.l/sqrt_det(ab1,6); 
}

complex rat_Z(int h1, int h2, int h3, complex r)
{
    static complex **b=new_complex2(9,9); for2(i,j,9) b[i][j]=0;
    static int hh[3]; hh[0]=h1, hh[1]=h2; hh[2]=h3;
    
    for2(i,j,3)
    {
        b[i][j]=A[hh[i]][hh[j]];
        b[3+i][j]=MA[hh[i]][hh[j]];
        b[6+i][j]=MMA[hh[i]][hh[j]];
        b[3+i][3+j]=MAM[hh[i]][hh[j]];
        b[6+i][3+j]=MMAM[hh[i]][hh[j]];
        b[6+i][6+j]=MMAMM[hh[i]][hh[j]];
        b[i][3+j]=AM[hh[i]][hh[j]];
        b[i][6+j]=AMM[hh[i]][hh[j]];
        b[3+i][6+j]=MAMM[hh[i]][hh[j]];
    }
    
    
    static complex **a=new_complex2(9,9); for2(i,j,9) a[i][j]=0;
    for2(i,j,3)
    {
        a[i][j]=MM[hh[i]][hh[j]]-conj(M[hh[i]][h1])*M[h1][hh[j]]-conj(M[hh[i]][h2])*M[h2][hh[j]]-conj(M[hh[i]][h3])*M[h3][hh[j]];
        a[i][3+j]=conj(M[hh[j]][hh[i]]); a[3+i][j]=M[hh[i]][hh[j]];
    }
    a[0][6]=-1; a[1][7]=-1; a[2][8]=-1; a[3][3]=-1; a[4][4]=-1; a[5][5]=-1; a[8][2]=-1; a[7][1]=-1; a[6][0]=-1;
    
    static complex **ab1=new_complex2(9,9);
    for2(k,l,9)
    {
        ab1[k][l]=0;
        for1(i, 9) ab1[k][l]-=r*a[k][i]*b[i][l];
        if (k==l) ab1[k][l]+=1;
    }
    
    return 1.l/sqrt_det(ab1,9); 
}


complex rat_Z(int h1, int h2, int h3, int h4, complex r)
{
    static complex **b=new_complex2(12,12); for2(i,j,12) b[i][j]=0;
    static int hh[4]; hh[0]=h1, hh[1]=h2; hh[2]=h3; hh[3]=h4;
    
    for2(i,j,4)
    {
        b[i][j]=A[hh[i]][hh[j]];
        b[4+i][j]=MA[hh[i]][hh[j]];
        b[8+i][j]=MMA[hh[i]][hh[j]];
        b[4+i][4+j]=MAM[hh[i]][hh[j]];
        b[8+i][4+j]=MMAM[hh[i]][hh[j]];
        b[8+i][8+j]=MMAMM[hh[i]][hh[j]];
        b[i][4+j]=AM[hh[i]][hh[j]];
        b[i][8+j]=AMM[hh[i]][hh[j]];
        b[4+i][8+j]=MAMM[hh[i]][hh[j]];
    }
    
    
    static complex **a=new_complex2(12,12); for2(i,j,12) a[i][j]=0;
    for2(i,j,4)
    {
        a[i][j]=MM[hh[i]][hh[j]]-conj(M[hh[i]][h1])*M[h1][hh[j]]-conj(M[hh[i]][h2])*M[h2][hh[j]]-conj(M[hh[i]][h3])*M[h3][hh[j]]-conj(M[hh[i]][h4])*M[h4][hh[j]];
        a[i][4+j]=conj(M[hh[j]][hh[i]]); a[4+i][j]=M[hh[i]][hh[j]];
    }
    a[0][8]=-1; a[1][9]=-1; a[2][10]=-1; a[3][11]=-1;  a[4][4]=-1; a[5][5]=-1; a[6][6]=-1; a[7][7]=-1; a[11][3]=-1; a[10][2]=-1; a[9][1]=-1; a[8][0]=-1;
    
    static complex **ab1=new_complex2(12,12);
    for2(k,l,12)
    {
        ab1[k][l]=0;
        for1(i, 12) ab1[k][l]-=r*a[k][i]*b[i][l];
        if (k==l) ab1[k][l]+=1;
    }
    
    return 1.l/sqrt_det(ab1,12); 
}



complex Z_direct(int h1, int h2, int h3, int h4, complex r)
{
    static complex ** MM0=new_complex2(N,N);
    for2(i,j, N) 
    {
        MM0[i][j]=0; if (i==j) MM0[i][j]=1;
        for1(l,N) if ((i-h1)*(j-h1)*(l-h1)!=0 && (i-h2)*(j-h2)*(l-h2)!=0 && (i-h3)*(j-h3)*(l-h3)!=0 && (i-h4)*(j-h4)*(l-h4)!=0) 
            MM0[i][j]-=r*conj(M[l][i])*M[l][j];
    }
    return 1.l/sqrt_det(MM0,N);    
}


int main()
{
    //cout<<N<<"\n";
//    complex ** M0=new_complex2 (N,N); 
    for2(n1, n2, N) {long double x, y; a_str>>x; a_str>>y; M[n1][n2]=2.l*(x+I*y);}
     
    
/*    for2(n1,n2,N)
    {
        M[n1][n2]=0;
        for1(n_, N/2) M[n1][n2]+=M0[n1][n_]*M0[n2][n_]*.88;
    }
*/
    ou1<<setprecision(24); ou2<<setprecision(24); ou3<<setprecision(24); ou4<<setprecision(24); 
  
    makeMMU();
    complex * z0=new complex [N]; for1(h, N) z0[h]=0;
    
   // int k0=31; complex Z0=0, Z0k=0;
    
    for1(nu, Nnu)
    {
        complex einu=exp(I*(nu*2*Pi/Nnu));
        
        complex z=make_products(einu); ou1<<nu*2*Pi/Nnu<<"  "<<real(z)<<"  "<<imag(z)<<"  "; ou2<<nu*2*Pi/Nnu<<"  "; ou3<<nu*2*Pi/Nnu<<"  "; ou4<<nu*2*Pi/Nnu<<"  ";
       // Z0+=z*exp(-2.l*(k0*nu*Pi/Nnu)*I);// cout<<z<<" "<<abs(z-Z_direct(-1,-1,-1,-1,einu))<<"   ";
        for1(h1,N) 
        {
            complex z1=z*rat_Z(h1, einu); ou1<<real(z1)<<"  "<<imag(z1)<<"  "; z0[h1]+=z1;
            //if (h1==N-1) 
        //    {Z0k+=z1*exp(-2.l*(k0*nu*Pi/Nnu)*I);}//cout<<z1<<"  "<<abs(z1-Z_direct(N-1,-1,-1,-1,einu))<< "\n";}
            for1(h2,h1) 
            {
                complex z2=z*rat_Z(h1, h2, einu); ou2<<real(z2)<<"  "<<imag(z2)<<"  ";
                for1(h3, h2) 
                {
                    complex z3=z*rat_Z(h1, h2, h3, einu); ou3<<real(z3)<<"  "<<imag(z3)<<"  ";
                    for1 (h4,h3) {complex z4=z*rat_Z(h1, h2, h3, h4, einu); ou4<<real(z4)<<"  "<<imag(z4)<<"  ";}
                }
            }
        }
        ou1<<"\n"; ou2<<"\n"; ou3<<"\n"; ou4<<"\n";
    }
    //for1(h, N) cout<<z0[h]<<"\n";
    
    if (!enable_checks) return 0;
    
    complex eiv=exp((2*Pi*rnd())*I); 
    int h1=rnd(N), h2, h3, h4; do h2=rnd(N); while (h2==h1); do h3=rnd(N); while ((h3==h1) || (h3==h2)); do h4=rnd(N); while ((h4==h1) || (h4==h2) || (h4==h3));  
    cout<<h1<<" "<<h2<<"  "<<h3<<"\n";
 
    makeMMU(); 
     cout<<"Z="<<make_products(1)<<" == "<<Z_direct(-1,-1,-1,-1,1)<<"\n"; 
     make_products(eiv);
     cout<<"Z1/Z="<<rat_Z(h1,eiv)<<" == "<<Z_direct(h1,-1,-1,-1,eiv)/Z_direct(-1,-1,-1,-1,eiv)<<"\n";
     //cout<<"Z2/Z="<<rat_Z(h1,h2,eiv)<<" == "<<Z_direct(h1,h2,-1,-1,eiv)/Z_direct(-1,-1,-1,-1,eiv)<<"\n";
     //cout<<"Z3/Z="<<rat_Z(h1,h2,h3,eiv)<<" == "<<Z_direct(h1,h2,h3,-1,eiv)/Z_direct(-1,-1,-1,-1,eiv)<<"\n";
     //cout<<"Z4/Z="<<rat_Z(h1,h2,h3,h4,eiv)<<" == "<<Z_direct(h1,h2,h3,h4,eiv)/Z_direct(-1,-1,-1,-1,eiv)<<"\n";
     
    // cout<<Z0k<<"  "<<Z0<<"  "<<1.l*N-Z0k/Z0<<"\n";
     

     
    return 0;
    
 
}
