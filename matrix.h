void EigenJacobi (long double ** h, long double ** a, long double * e, int n, long double accuracy=1e-10, int max_sweep=100)
{
  int flag; long double ** v=new long double * [n]; {for (int i=0; i<n; i++) v[i]=new long double [n];}
  
  {for (int i=0; i<n; i++) for (int j=0; j<n; j++) {a[i][j]=h[i][j]; if (i==j) v[i][j]=1; else v[i][j]=0;};}

	do
	{
		max_sweep--;
		flag=0;
		for (int pp=n-1; pp>=1; pp--)
		for (int q=0; q<n-pp; q++)
		{
			int p=pp+q;
			if (fabs(a[p][q])>accuracy)
			{
				flag=1;
				long double theta=(a[q][q]-a[p][p])/(2*a[p][q]);
				long double t=1/(abs(theta)+sqrt(theta*theta+1));if (theta<0) t=-t;
				long double c=1/sqrt(t*t+1), s=t*c, tau=s/(1+c);
				for (int r=0; r<n; r++)
				{
					if ((r!=p) && (r!=q))
					{
						long double a1=a[r][p]-s*(a[r][q]+tau*a[r][p]);
						long double a2=a[r][q]+s*(a[r][p]-tau*a[r][q]);
						a[r][p]=a1; a[p][r]=a1;
						a[r][q]=a2; a[q][r]=a2;


					;}
						long double v1=c*v[p][r]-s*v[q][r];
						long double v2=s*v[p][r]+c*v[q][r];
						v[p][r]=v1; v[q][r]=v2;
				;}
				long double ap=a[p][p]-t*a[p][q], aq=a[q][q]+t*a[p][q];
				a[p][p]=ap; a[q][q]=aq;
				a[p][q]=0;  a[q][p]=0;
			;}
		;}
	;}
	while( (flag>0) && (max_sweep>0) );
	
	
	for (int i=0; i<n; i++) e[i]=a[i][i];

/*
	for (int i=0; i<n-1; i++)
	for (int j=i+1; j<n; j++)
		if (e[j]<e[i])
		{
			long double a=e[j]; e[j]=e[i]; e[i]=a;
			long double * b=v[j]; v[j]=v[i];v[i]=b;
		;}
*/
	{for (int i=0; i<n; i++) for (int j=0; j<n; j++) a[i][j]=v[i][j];}
	


	for (int i=0; i<n; i++) //normalization
	{
		long double s=0;
		for (int j=0; j<n; j++) s+=a[i][j]*a[i][j];
		s=sqrt(s);
		if (fabs(s)>accuracy*accuracy) {for (int j=0; j<n; j++) a[i][j]/=s;}
		else {for (int j=0; j<n; j++) a[i][j]=1/sqrt(n);}
	;}
	
	{for (int i=0; i<n; i++) delete [] v[i]; delete [] v;}
;}



void EigenJacobi (complex ** h, complex ** a, long double * e, int n, long double accuracy=1e-10, int max_sweep=100)
{
    static complex ** v=new_complex2(n,n);
    static int nn=n; if (n>nn) {delete_complex2(v,nn,nn);  v=new_complex2(n,n); nn=n;}
    
    {for (int i=0; i<n; i++) for (int j=0; j<n; j++) {a[i][j]=h[i][j]; if (i==j) v[i][j]=1; else v[i][j]=0;};}
 
    for1(sweep, max_sweep)
    {
        int f=0;
        for1(p,n) for1(q,p)
//		for (int p=n-1; p>0; p--)
//		for (int q=0; q<n-p; q++)
        if (abs(a[p][q])>accuracy)
        {
            complex z=sqrt(I*conj(a[p][q])/abs(a[p][q])), x1,x2;
            for1(j,n) {x1=a[p][j]*z/sqrt(2.l); x2=a[q][j]*conj(z)/sqrt(2.l);a[p][j]=x1-x2; a[q][j]=x1+x2;}
//            for1(j,n) {x1=v[p][j]*z/sqrt(2); x2=v[q][j]*conj(z)/sqrt(2);v[p][j]=x1-x2; v[q][j]=x1+x2;}
            for1(j,n) {x1=a[j][p]*conj(z)/sqrt(2.l); x2=a[j][q]*z/sqrt(2.l);a[j][p]=x1-x2; a[j][q]=x1+x2;}
            for1(j,n) {x1=v[j][p]*conj(z)/sqrt(2.l); x2=v[j][q]*z/sqrt(2.l);v[j][p]=x1-x2; v[j][q]=x1+x2;}
             
  
            z=sqrt(conj(a[p][q])/abs(a[p][q]));
            for1(j,n) {x1=a[p][j]*z/sqrt(2.l); x2=a[q][j]*conj(z)/sqrt(2.l);a[p][j]=x1-x2; a[q][j]=x1+x2;}
//            for1(j,n) {x1=v[p][j]*z/sqrt(2); x2=v[q][j]*conj(z)/sqrt(2);v[p][j]=x1-x2; v[q][j]=x1+x2;}
            for1(j,n) {x1=a[j][p]*conj(z)/sqrt(2.l); x2=a[j][q]*z/sqrt(2.l);a[j][p]=x1-x2; a[j][q]=x1+x2;}
            for1(j,n) {x1=v[j][p]*conj(z)/sqrt(2.l); x2=v[j][q]*z/sqrt(2.l);v[j][p]=x1-x2; v[j][q]=x1+x2;}
            
            
            f=1;
            
        }
        if (f==0) break;
      
    }

    
	for1 (i,n) e[i]=real(a[i][i]);    
	for2(i,j,n) a[j][i]=v[i][j];

	


	for (int i=0; i<n; i++) //normalization
	{
		complex s=0;
		for (int j=0; j<n; j++) s+=norm2(a[i][j]);
		s=sqrt(s); //if (abs(a[i][0])!=0) s*=(a[i][0])/abs(a[i][0]);
		if (abs(s)>accuracy*accuracy) {for (int j=0; j<n; j++) a[i][j]/=s;}
		else {for (int j=0; j<n; j++) a[i][j]=1.l/sqrt(n);}
	;}

}



void rotate(long double ** &a, long double ** & v, int size)
{
	long double ** b=new long double * [size]; {for (int i=0; i<size; i++) b[i]=new long double [size];}
	
	{
		for (int i=0;i<size;i++)
		for (int j=0;j<size;j++)
		{
			 b[i][j]=0; 
			 for (int k=0; k<size; k++) b[i][j]+=a[i][k]*v[j][k]; //check index order!!!
		;}
	;}

	{
		for (int i=0;i<size;i++)
		for (int j=0;j<size;j++)
		{
			 a[i][j]=0; 
			 for (int k=0; k<size; k++) a[i][j]+=v[i][k]*b[k][j]; //check index order!!!
		;}
	;}

	
	{for (int i=0; i<size; i++) delete [] b[i];} delete [] b;
;}

void mult_left(complex ** a, complex **  r, int size) //a=r*a
{
  static int old_size=0;
  static complex ** y=NULL;
  if (size>old_size)
  {
    for (int i=0; i<old_size; i++) delete [] y[i]; delete [] y; y=new complex * [size]; for (int i=0; i<size; i++) y[i]=new complex [size]; old_size=size;
  }

  for (int i=0; i<size; i++) for (int j=0; j<size; j++) y[i][j]=0;
  #pragma omp parallel for
  for (int i=0; i<size; i++) for (int k=0; k<size; k++) 
      if (abs(r[i][k])>1e-15) for (int j=0; j<size; j++) y[i][j]+=r[i][k]*a[k][j];
  for (int i=0; i<size; i++) for (int j=0; j<size; j++) a[i][j]=y[i][j]; 
  
}

void mult_right(complex ** a, complex **  r, int size) //a=a*r
{
  static int old_size=0;
  static complex ** y=NULL;
  if (size>old_size)
  {
    for (int i=0; i<old_size; i++) delete [] y[i]; delete [] y; y=new complex * [size]; for (int i=0; i<size; i++) y[i]=new complex [size]; old_size=size;
  }

  for (int i=0; i<size; i++) for (int j=0; j<size; j++) y[i][j]=0;
  #pragma omp parallel for
  for (int j=0; j<size; j++) for (int k=0; k<size; k++) 
      if (abs(r[k][j])>1e-15) for (int i=0; i<size; i++) y[i][j]+=a[i][k]*r[k][j];
  for (int i=0; i<size; i++) for (int j=0; j<size; j++) a[i][j]=y[i][j]; 
  
}


#define c_number long double

void Inverse(c_number ** &a, int size) //Gauss with partial pivoting
{
	c_number ** r=new c_number * [size];
	{
		for (int i=0;i<size;i++)
		{
			r[i]=new c_number [size]; 
			for (int j=0; j<size; j++) r[i][j]=0.; r[i][i]=1.;
		;}
	;}
	
	for (int i=0; i<size; i++)
	{
		//pivoting
      {
		long double fmax=-1.; int jmax;
		for (int j=i; j<size; j++)
			if ( fmax<norm2(a[j][i]) ) {jmax=j; fmax=norm2(a[j][i]);}
		c_number * aux; 
		aux=a[i]; a[i]=a[jmax]; a[jmax]=aux;
		aux=r[i]; r[i]=r[jmax]; r[jmax]=aux;
		if (norm2(a[i][i])<=0.) {cout<<"LInverse!"; return;}
      ;}
		
		//main body
	
#pragma omp parallel for
		for (int j=0; j<i; j++)
		{
			c_number f=a[j][i]/a[i][i];a[j][i]=0;
			{
			for (int l=i+1; l<size; l++) a[j][l]-=a[i][l]*f; 
			for (int l=0; l<size; l++) r[j][l]-=r[i][l]*f;
			}
		;}
		{
			c_number f=1.l/a[i][i];
			{for (int l=i; l<size; l++) a[i][l]*=f;}
			{for (int l=0; l<size; l++) r[i][l]*=f;}
		;}
#pragma omp parallel for
		for (int j=i+1; j<size; j++)
		{
			c_number f=a[j][i]/a[i][i]; a[j][i]=0;
			{
			for (int l=i+1; l<size; l++) a[j][l]-=a[i][l]*f;
			for (int l=0; l<size; l++) r[j][l]-=r[i][l]*f;
			}
		;}
		
        ;}
	
	{for (int i=0; i<size; i++) for (int j=0; j<size; j++) a[i][j]=r[i][j];}
	{for (int i=0; i<size; i++) delete [] r[i]; delete [] r;}
;}

void Inverse(complex ** &a, int size) //Gauss with partial pivoting
{
	complex ** r=new complex * [size];
	{
		for (int i=0;i<size;i++)
		{
			r[i]=new complex [size]; 
			for (int j=0; j<size; j++) r[i][j]=0.; r[i][i]=1.;
		;}
	;}
	
	for (int i=0; i<size; i++)
	{
		//pivoting
      {
		long double fmax=-1.; int jmax;
		for (int j=i; j<size; j++)
			if ( fmax<norm2(a[j][i]) ) {jmax=j; fmax=norm2(a[j][i]);}
		complex * aux; 
		aux=a[i]; a[i]=a[jmax]; a[jmax]=aux;
		aux=r[i]; r[i]=r[jmax]; r[jmax]=aux;
		if (norm2(a[i][i])<=0.) {cout<<"LInverse!"; return;}
      ;}
		
		//main body
	
#pragma omp parallel for
		for (int j=0; j<i; j++)
		{
			complex f=a[j][i]/a[i][i];a[j][i]=0;
			{
			for (int l=i+1; l<size; l++) a[j][l]-=a[i][l]*f; 
			for (int l=0; l<size; l++) r[j][l]-=r[i][l]*f;
			}
		;}
		{
			complex f=1.l/a[i][i];
			{for (int l=i; l<size; l++) a[i][l]*=f;}
			{for (int l=0; l<size; l++) r[i][l]*=f;}
		;}
#pragma omp parallel for
		for (int j=i+1; j<size; j++)
		{
			complex f=a[j][i]/a[i][i]; a[j][i]=0;
			{
			for (int l=i+1; l<size; l++) a[j][l]-=a[i][l]*f;
			for (int l=0; l<size; l++) r[j][l]-=r[i][l]*f;
			}
		;}
		
        ;}
	
	{for (int i=0; i<size; i++) for (int j=0; j<size; j++) a[i][j]=r[i][j];}
	{for (int i=0; i<size; i++) delete [] r[i]; delete [] r;}
;}


void Fourier(complex * x, complex * y, int size, int dir) //dir==-1 => r to k; dir==1 => k to r
{
  long double r=dir*2*Pi/size;
  for (int j=0; j<size; j++) {y[j]=0; for (int l=0; l<size; l++) y[j]+=x[l]*exp(I*(r*j*l));}
  if (dir==1) for (int i=0; i<size; i++) y[i]/=1.l*size;
;}

void Fourier(complex * x, int size, int dir) //dir==-1 => r to k; dir==1 => k to r
{
  static complex * y=NULL; static int old_size=0;
  if (size>old_size) {delete[]y; y=new complex [size];old_size=size;}
  Fourier(x,y, size, dir);
  for (int i=0; i<size; i++) x[i]=y[i];
}

void kk_elements (complex ** x, complex *y, int size)  
{
  long double r=2*Pi/size;
  for (int k=0; k<size; k++)
  {
    y[k]=0;
    for (int i=0; i<size; i++)
    for (int j=0; j<size; j++) 
      y[k]+=x[i][j]*exp(I*(r*(i-j)));
  }
  for (int k=0; k<size; k++) y[k]/=1.l*size;
}


void Fourier (complex ** x, int size, int dir)
{
    static int old_size=0;
  static complex ** y=NULL;
  if (size>old_size)
  {
    for (int i=0; i<old_size; i++) delete [] y[i]; delete [] y; y=new complex * [size]; for (int i=0; i<size; i++) y[i]=new complex [size]; old_size=size;
  }
//long double s; for (int i=0; i<size; i++) for (int j=0; j<size; j++) s+=abs(sqr(x[i][j]));
  long double r=dir*2*Pi/size; 
  for (int i=0; i<size; i++) for (int j=0; j<size; j++) y[i][j]=exp(I*(r*i*j));
  mult_left(x,y,size);
  for (int i=0; i<size; i++) for (int j=0; j<size; j++)  y[j][i]=exp(-I*(r*i*j));
  mult_right(x,y,size);
  for (int i=0; i<size; i++) for (int j=0; j<size; j++) x[i][j]/=(1.l*size);
  
//  long double s2; for (int i=0; i<size; i++) for (int j=0; j<size; j++) s2+=abs(sqr(x[i][j]));  cout<<"\n"<<s<<"="<<s2<<"\n";
}


complex det(complex ** & b, int N)
{
    
     
 	complex s=1;  
    
	for (int i=0; i<N-1; i++)
	{
     //pivoting
      {
		long double fmax=0; int jmax=i;
		for (int j=i; j<N; j++)
		if ( fmax<abs(b[j][i]) ) {jmax=j; fmax=abs(b[j][i]);}
		if (jmax!=i)
		{
			complex * aux;
			aux=b[i]; b[i]=b[jmax]; b[jmax]=aux;
			s=-s;
		;}
      ;}
		//main body
		if (norm2(b[i][i])<=0) return 0;
		for (int j=i+1; j<N; j++)
		{
		complex f=b[j][i]/b[i][i];
		for (int k=i; k<N; k++) b[j][k]-=f*b[i][k];
		;}
	;}

 	
   if (abs(b[N-1][N-1])<=0) return 0;
	{for (int i=0; i<N; i++) s*=b[i][i];}
   return s;
;}

long double det(long double ** & b, int N)
{
    
     
 	long double s=1;  
    
	for (int i=0; i<N-1; i++)
	{
     //pivoting
      {
		long double fmax=0; int jmax=i;
		for (int j=i; j<N; j++)
		if ( fmax<abs(b[j][i]) ) {jmax=j; fmax=abs(b[j][i]);}
		if (jmax!=i)
		{
			long double * aux;
			aux=b[i]; b[i]=b[jmax]; b[jmax]=aux;
			s=-s;
		;}
      ;}
		//main body
		if (abs(b[i][i])<=0) return 0;
		for (int j=i+1; j<N; j++)
		{
		long double f=b[j][i]/b[i][i];
		for (int k=i; k<N; k++) b[j][k]-=f*b[i][k];
		;}
	;}

 	
   if (abs(b[N-1][N-1])<=0) return 0;
	{for (int i=0; i<N; i++) s*=b[i][i];}
   return s;
;}

