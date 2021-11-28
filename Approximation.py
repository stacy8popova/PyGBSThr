# -*- coding: utf-8 -*-
"""
@author: A. Popova 
"""


import numpy as np
old_settings = np.seterr(all='ignore') 


def fact(x):
    res = 1
    for i in range(int(x)):
        res *= (i+1.)
    return res

def n_combinations(nmodes, n_moments):
    
    return round(fact(nmodes)/(fact(nmodes-n_moments)*fact(n_moments))) 
    
def in_Z(data_minors_n, pn, Nu):
    
    Z_v_n = np.zeros((pn, Nu),dtype=np.complex128)
    
    for j in range(Nu):
        for n in range(0,2*pn,2):
            Z_v_n[n//2,j] =  data_minors_n[j,int(1+n)] + 1j*data_minors_n[j,int(2+n)] 
            
    return Z_v_n 

def fft_Z(Z_v_n, pn, Nu):
    
    Z_v_nf = np.zeros((pn, Nu), dtype = np.complex128)
    
    Z_v_nf = np.fft.fft(Z_v_n)/Nu
    
    return Z_v_nf


def gauss_fun(x, *args):
    
    for c in args:
        c = args
    
    if len(c) == 3:    
    
        res = c[0]*np.exp(-(x - c[1])**2/(2*c[2])) 
        
    if len(c) == 4:
        
        res = c[0]*np.exp(-(x - c[1])**2/(2*c[2])) * np.exp(+ c[3]*(x - c[1])**3/(6*c[2]**3)) 
        
    if len(c) == 5:
        
        res = c[0]*np.exp(-(x - c[1])**2/(2*c[2])) * np.exp(+ c[3]*(x - c[1])**3/(6*c[2]**3)) * np.exp(+ c[4]*(x - c[1])**4/(8*c[2]**4))
    
        
    return res


data_ = np.genfromtxt(r'in_out/Submatrix.dat')
m = len(data_) # number of clicked detectors 

Nu = int(10*m)
dnu = 2*np.pi/Nu

data_minors = np.genfromtxt(r'in_out/Minors0-1.dat')
    
Z_v_0 = np.zeros((Nu),dtype=np.complex128)

Z_v_1 = np.zeros((m, Nu),dtype=np.complex128)

for j in range(Nu):
     Z_v_0[j] =  data_minors[j,1] + 1j*data_minors[j,2]
        
for j in range(Nu):
    for n in range(0,2*m,2):
        Z_v_1[n//2,j] =  data_minors[j,int(3+n)] + 1j*data_minors[j,int(4+n)] 


Z_v_0f =  np.zeros((Nu), dtype = np.complex128)
Z_v_1f = np.zeros((m, Nu), dtype = np.complex128)

Z_v_0f = np.fft.fft(Z_v_0)/Nu
Z_v_1f = np.fft.fft(Z_v_1)/Nu


mean_ = np.zeros(Nu)
disp_ = np.zeros(Nu)
m3_ = np.zeros(Nu)
m4_ =np.zeros(Nu)

Mu1_2 = np.zeros(Nu)
Mu2_2 = np.zeros(Nu)
A_2 = np.zeros(Nu)

Mu1_3 = np.zeros(Nu)
Mu2_3 = np.zeros(Nu)
Mu3_3 = np.zeros(Nu)
A_3 = np.zeros(Nu)

Mu1_4 = np.zeros(Nu)
Mu2_4 = np.zeros(Nu)
Mu3_4 = np.zeros(Nu)
Mu4_4 = np.zeros(Nu)
A_4 = np.zeros(Nu)



data_mom = np.genfromtxt(r'in_out/Moments.dat')

for j in range(Nu):
    mean_[j] = data_mom[j, 0]
    disp_[j] = data_mom[j, 1]
    m3_[j] = data_mom[j, 2]
    m4_[j] = data_mom[j, 3]

# Approximation 

mu0 = np.zeros(Nu)
mu1 = np.zeros(Nu)
mu2 = np.zeros(Nu)
mu3 = np.zeros(Nu)
mu4 = np.zeros(Nu)

for nu in range(Nu): 
    mu0[nu] = (Z_v_0f[nu]/Z_v_0[0]).real
    mu1[nu] = mean_[nu]
    mu2[nu] = disp_[nu]
    mu3[nu] = m3_[nu]
    mu4[nu] = m4_[nu]

n_cut = int(m+1)

# 2 order 


for nu in range(Nu):

    A_2[nu] =  mu0[nu] 
    Mu1_2[nu] = mu1[nu] 
    Mu2_2[nu] = mu2[nu] 


    for z in range(300):
        s0 = 0
        s1 = 0 
        s2 = 0 
        s3 = 0
        s4 = 0

        for j in range(n_cut):

            s0 +=  gauss_fun(j, A_2[nu], Mu1_2[nu], Mu2_2[nu])
            s1 +=  gauss_fun(j, A_2[nu], Mu1_2[nu], Mu2_2[nu])* j
            s2 +=  gauss_fun(j, A_2[nu], Mu1_2[nu], Mu2_2[nu])* j**2
            
        if s0==s0:

            mu0_ = s0 
            mu1_ = s1/s0 
            mu2_ = s2/s0 - mu1_**2
            
            # To improve accuracy you can vary 'step_i' 
            
            step_0 =  0.1
            step_1 =  0.1
            step_2 =  0.1

            A_2[nu] += step_0*(mu0[nu]  - mu0_ )
            Mu1_2[nu] += step_1*(mu1[nu] - mu1_) 
            Mu2_2[nu] += step_2*(mu2[nu] - mu2_)  
            
        else:
            
            A_2[nu] += 0
            Mu1_2[nu] += 0
            Mu2_2[nu] += 0
            
            
# 3 order 


for nu in range(Nu):

    A_3[nu] = mu0[nu]
    Mu1_3[nu] = mu1[nu] 
    Mu2_3[nu] = mu2[nu] 
    Mu3_3[nu] = 0
    diverg_tracking = []


    for z in range(500):
        s0 = 0
        s1 = 0 
        s2 = 0 
        s3 = 0

        for j in range(n_cut):

            s0 += gauss_fun(j, A_3[nu], Mu1_3[nu], Mu2_3[nu], Mu3_3[nu])
            s1 += gauss_fun(j, A_3[nu], Mu1_3[nu], Mu2_3[nu], Mu3_3[nu])* j
            s2 += gauss_fun(j, A_3[nu], Mu1_3[nu], Mu2_3[nu], Mu3_3[nu])* j**2
            s3 += gauss_fun(j, A_3[nu], Mu1_3[nu], Mu2_3[nu], Mu3_3[nu])* j**3
            
        
         # attention! s0 may diverge 
        
        diverg_tracking.append(s0)
        
        s0_mean = np.mean(diverg_tracking)
        
        if s0==s0 and 10**(-2) < s0/s0_mean < 10**(1):

            mu0_ = s0 
            mu1_ = s1/s0 
            mu2_ = s2/s0 - mu1_**2
            mu3_ = s3/s0 - 3*mu2_*mu1_ - mu1_**3
            
            # To improve accuracy you can vary 'step_i'  
            
            step_0 =  0.05
            step_1 =  0.1
            step_2 =  0.1
            step_3 =  0.05
            
            A_3[nu] += step_0*(mu0[nu] - mu0_)
            Mu1_3[nu] += step_1*(mu1[nu] - mu1_ )
            Mu2_3[nu] += step_2*(mu2[nu] - mu2_)  
            Mu3_3[nu] += step_3*(mu3[nu] - mu3_) 
            
        else:
            A_3[nu] += 0
            Mu1_3[nu] += 0
            Mu2_3[nu] += 0
            Mu3_3[nu] += 0 
        

for nu in range(Nu):

    A_4[nu] = mu0[nu]
    Mu1_4[nu] = mu1[nu] 
    Mu2_4[nu] = mu2[nu] 
    Mu3_4[nu] = 0 
    Mu4_4[nu] = 0 
    
    diverg_tracking = []

    for z in range(700):

        s0 = 0
        s1 = 0 
        s2 = 0 
        s3 = 0
        s4 = 0

        for j in range(n_cut):

            s0 +=  gauss_fun(j, A_4[nu], Mu1_4[nu], Mu2_4[nu], Mu3_4[nu], Mu4_4[nu])
            s1 +=  gauss_fun(j, A_4[nu], Mu1_4[nu], Mu2_4[nu], Mu3_4[nu], Mu4_4[nu])* j
            s2 +=  gauss_fun(j, A_4[nu], Mu1_4[nu], Mu2_4[nu], Mu3_4[nu], Mu4_4[nu])* j**2
            s3 +=  gauss_fun(j, A_4[nu], Mu1_4[nu], Mu2_4[nu], Mu3_4[nu], Mu4_4[nu])* j**3
            s4 +=  gauss_fun(j, A_4[nu], Mu1_4[nu], Mu2_4[nu], Mu3_4[nu], Mu4_4[nu])* j**4

        # attention! s0 may diverge 
        
        diverg_tracking.append(s0)
        
        s0_mean = np.mean(diverg_tracking)
        
        if s0==s0 and 10**(-2) < s0/s0_mean < 10**(1):
            mu0_ = s0 
            mu1_ = s1/s0 
            mu2_ = s2/s0 - mu1_**2
            mu3_ = s3/s0 - 3*mu2_*mu1_ - mu1_**3
            mu4_ = s4/s0 - 4*mu3_*mu1_ - 3*mu2_**2 - 6*mu2_*mu1_**2 - mu1_**4
            
            # To improve accuracy you can vary 'step_i' 
            step_0 =  0.1
            step_1 =  0.1
            step_2 =  0.1
            step_3 =  0.05
            step_4 =  0.008
            
            A_4[nu] +=  step_0*(mu0[nu] - mu0_) 
            Mu1_4[nu] +=  step_1*(mu1[nu] - mu1_)
            Mu2_4[nu] +=  step_2*(mu2[nu] - mu2_) 
            Mu3_4[nu] +=  step_3*(mu3[nu] - mu3_) 
            Mu4_4[nu] +=  step_4*(mu4[nu] - mu4_)
            
        else:
            A_4[nu] += 0
            Mu1_4[nu] += 0
            Mu2_4[nu] += 0 
            Mu3_4[nu] += 0 
            Mu4_4[nu] += 0
        

# Export 

with open(r"in_out/Approximation.dat", 'w') as ouf:
    ouf.write( str('nu') + '\t' + str('A2')+ '\t' +  str( 'M12')+ '\t'+ str('M22') +'\t'  + str('A3') +'\t' + str('M13') + '\t' +  str( 'M23') + '\t' +  str('M33') +'\t' + str( 'A4') + '\t' + str('M14') + '\t' + str( 'M24') +'\t'+ str('M34') + '\t' +str('M44') +'\n' ) 
    for k in range(Nu):
        ouf.write( str(k) + '\t' + str(A_2[k])+ '\t' +  str(Mu1_2[k])+ '\t'+ str(Mu2_2[k]) +'\t'  + str(A_3[k]) +'\t' + str(Mu1_3[k]) + '\t' +  str( Mu2_3[k]) + '\t' +  str(Mu3_3[k]) +'\t' + str( A_4[k]) + '\t' + str(Mu1_4[k]) + '\t' + str( Mu2_4[k]) +'\t'+ str( Mu3_4[k]) + '\t'+ str(Mu4_4[k]) +'\n' ) 
        
