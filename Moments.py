# -*- coding: utf-8 -*-
"""
@author: A. Popova
"""

import numpy as np

#old_settings = np.seterr(all='ignore')


def fact(x):
    res = 1
    for i in range(int(x)):
        res *= (i+1.)
    return res


data_ = np.genfromtxt(r'in_out/Submatrix.dat')
m = len(data_)


Nu = int(10*m)
dnu = 2*np.pi/Nu


M =  np.zeros((m, m), dtype=np.complex128)

real_part = []
imaginary_part = []

for i in range(m):
    for k in range(0,2*m,2):
        real_part.append(data_[i,k])
        
for i in range(m):
    for k in range(1,2*m+1,2):
        imaginary_part.append(data_[i,k])

for i in range(m**2):
    M[i//m,i%m] = real_part[i] + 1j*  imaginary_part[i]
    


# Import Minors

data_minors = np.genfromtxt(r'in_out/Minors0-1.dat')
data_minors2 = np.genfromtxt(r'in_out/Minors2.dat')
data_minors3 = np.genfromtxt(r'in_out/Minors3.dat')
data_minors4 = np.genfromtxt(r'in_out/Minors4.dat')




p2 = round(fact(m)/(fact(m-2)*2)) 
p3 = round(fact(m)/(fact(m - 3)*fact(3)))
p4 = round(fact(m)/(fact(m - 4)*fact(4)))


Z_v_0 = np.zeros((Nu),dtype=np.complex128)

Z_v_1 = np.zeros((m, Nu),dtype=np.complex128)

Z_v_2 = np.zeros((p2, Nu),dtype=np.complex128)

Z_v_3 = np.zeros((p3, Nu),dtype=np.complex128)

Z_v_4 = np.zeros((p4, Nu),dtype=np.complex128)

for j in range(Nu):
     Z_v_0[j] =  data_minors[j,1:2] + 1j*data_minors[j,2:3]
        
for j in range(Nu):
    for n in range(0,2*m,2):
        Z_v_1[n//2,j] =  data_minors[j,int(3+n)] + 1j*data_minors[j,int(4+n)] 

for j in range(Nu):
    for n in range(0,2*p2,2):
        Z_v_2[n//2,j] =  data_minors2[j,int(1+n)] + 1j*data_minors2[j,int(2+n)] 
        
for j in range(Nu):
    for n in range(0,2*p3,2):
        Z_v_3[n//2,j] =  data_minors3[j,int(1+(n))] + 1j*data_minors3[j,int(2+(n))] 
        
for j in range(Nu):
    for n in range(0,2*p4,2):
        Z_v_4[n//2,j] =  data_minors4[j,int(1+(n))] + 1j*data_minors4[j,int(2+(n))] 
        
        
Z_v_0f = np.fft.fft(Z_v_0)/Nu
Z_v_1f = np.fft.fft(Z_v_1)/Nu            
Z_v_2f = np.fft.fft(Z_v_2)/Nu
Z_v_3f = np.fft.fft(Z_v_3)/Nu
Z_v_4f = np.fft.fft(Z_v_4)/Nu 
            

# a range of computed sectors 

Nuk = Nu

mean_ = np.zeros(Nuk)
disp_ = np.zeros(Nuk)
m3_ = np.zeros(Nuk)
m4_ = np.zeros(Nuk)
m5_ = np.zeros(Nuk)

def moment_formula(n, *args):
    
    m = 0 
    
    for x in args:
            moments = x
    
    if n == 2:
        m = moments[0] + 2*moments[1] - moments[0]**2
        
    if n == 3:
        m = moments[0] + 6*moments[1] + 6*moments[2] - 3*mean_[nu]*(moments[0] + 2*moments[1]) +  2*moments[0]**3
       
    if n == 4:
        m_2 = moments[0] + 2*moments[1]
        
        m_3 = moments[0] + 6*moments[1] + 6*moments[2]
        
        m_4 = moments[0] + 14*moments[1] + 36*moments[2] + 24*moments[3]
        
        m =  m_4 - 4*m_3*moments[0]- 3*m_2**2 + 12*m_2*moments[0]**2 - 6*moments[0]**4
        
    return m


n_ij_v =  np.zeros(Nuk)
n_ijk_v = np.zeros(Nuk)
n_ijkl_v = np.zeros(Nuk)
n_ijklp_v = np.zeros(Nuk)

ind_2 = []
ind_3 = []
ind_4 = []


for i in range(m):
        for j in range(i+1, m):
            ind_2.append([i,j]) 

for i in range(m):
        for j in range(i+1, m):
            for k in range(j+1, m):
                ind_3.append([i,j,k]) 
                
for i in range(m):
    for j in range(i+1, m):
        for k in range(j+1, m):
            for l in range(k+1, m):
                ind_4.append([i,j,k,l]) 


for z in range(Nuk): 
    for j in range(m):
        mean_[z] += 1 - (Z_v_1f[j,z]/Z_v_0f[z]).real
        

for nu in range(Nuk):
    i_ = 0
    for i in range(m):
        for j in range(i+1, m):
            n_ij_v[nu] += 1 - (( Z_v_1f[j,nu] + Z_v_1f[i,nu] - Z_v_2f[i_,nu])/Z_v_0f[nu]).real
            i_ += 1
    disp_[nu] =  moment_formula(2, [mean_[nu], n_ij_v[nu]])
            

for nu in range(Nuk):
    i_= 0
    for i in range(m):
        for j in range(i+1, m):
            for k in range(j+1, m):
                
                z1 = ind_2.index([i,j])
                z2 = ind_2.index([i,k])
                z3 = ind_2.index([j,k])
                
                n_ijk_v[nu] += 1 - ((Z_v_1f[i,nu] + Z_v_1f[j,nu] + Z_v_1f[k,nu] - Z_v_2f[z1,nu] - Z_v_2f[z2,nu] - Z_v_2f[z3,nu] + Z_v_3f[i_,nu])/Z_v_0f[nu]).real
                i_ += 1 
                
    m3_[nu] = moment_formula(3, [mean_[nu], n_ij_v[nu], n_ijk_v[nu]])
    
for nu in range(Nuk): 
    i_= 0
    for i in range(m):
        for j in range(i+1, m):
            for k in range(j+1, m):
                for l in range(k+1, m):

                    z1 = ind_2.index([i,j])
                    z2 = ind_2.index([i,k])
                    z3 = ind_2.index([i,l])

                    z4 = ind_2.index([j,k])
                    z5 = ind_2.index([k,l])
                    z6 = ind_2.index([j,l])

                    h1 = ind_3.index([i,j,k])
                    h2 = ind_3.index([j,k,l])
                    h3 = ind_3.index([i,k,l])
                    h4 = ind_3.index([i,j,l])

                    n_ijkl_v[nu] += 1 - ((Z_v_1f[i,nu] + Z_v_1f[j,nu] + Z_v_1f[k,nu] + Z_v_1f[l,nu] - Z_v_2f[z1,nu] - Z_v_2f[z2,nu] - Z_v_2f[z3,nu] - Z_v_2f[z4,nu] - Z_v_2f[z5,nu] - Z_v_2f[z6,nu] + Z_v_3f[h1,nu] + Z_v_3f[h2,nu] +  Z_v_3f[h3,nu] + Z_v_3f[h4,nu] -  Z_v_4f[i_,nu])/Z_v_0f[nu]).real 
                    
                    i_ += 1 

    m4_[nu] = moment_formula(4, [mean_[nu], n_ij_v[nu], n_ijk_v[nu], n_ijkl_v[nu]])


# Export Moments

with open( r"in_out/Moments.dat", 'w') as ouf:
    for nu in range(Nuk):
        ouf.write( str(mean_[nu].real) + '\t' + str(disp_[nu].real) +'\t'  + str(m3_[nu].real) +'\t'  + str(m4_[nu].real) +'\t' ) 
        if nu < (Nu+1):
            ouf.write('\n')
            



