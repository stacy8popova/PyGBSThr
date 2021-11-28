# -*- coding: utf-8 -*-
"""
@author: Popova A.S. 
"""

import numpy as np
old_settings = np.seterr(all='ignore') 

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

def Z(M):
    
    M = M.conjugate().T@M
    
    II = np.eye(len(M))
    
    z = (np.linalg.det(II - 4*M))**(-0.5)
    
    return z.real 

# Here we find the normalization factor for the probability 

data_ = np.genfromtxt(r'in_out/Submatrix.dat')
m = len(data_)

M_sub =  np.zeros((m, m), dtype=np.complex128)

real_part = []
imaginary_part = []

for i in range(m):
    for k in range(0,2*m,2):
        real_part.append(data_[i,k])
        
for i in range(m):
    for k in range(1,2*m+1,2):
        imaginary_part.append(data_[i,k])

for i in range(m**2):
    M_sub[i//m,i%m] = real_part[i] + 1j*  imaginary_part[i]
    
data_GBS = np.genfromtxt(r'in_out/GBS_matrix.dat')
m_ = len(data_GBS)

M =  np.zeros((m_, m_), dtype=np.complex128)

real_part = []
imaginary_part = []

for i in range(m):
    for k in range(0,2*m,2):
        real_part.append(data_GBS[i,k])
        
for i in range(m):
    for k in range(1,2*m+1,2):
        imaginary_part.append(data_GBS[i,k])

for i in range(m**2):
    M[i//m,i%m] = real_part[i] + 1j*  imaginary_part[i]
    

norm = Z(M_sub)/Z(M)

# The import the exact result 
 
data_exact = np.genfromtxt(r'in_out/Exact.dat')

exact_p = data_exact[m,1]*norm

# The import the approximate coefficients of the distribution of each sector  
Nu = int(10*m)

A_2 = np.zeros(Nu)
Mu1_2 = np.zeros(Nu)
Mu2_2 = np.zeros(Nu)

A_3 = np.zeros(Nu)
Mu1_3 = np.zeros(Nu)
Mu2_3 = np.zeros(Nu)
Mu3_3 = np.zeros(Nu)

A_4 = np.zeros(Nu)
Mu1_4 = np.zeros(Nu)
Mu2_4 = np.zeros(Nu)
Mu3_4 = np.zeros(Nu)
Mu4_4 = np.zeros(Nu)

data_cum = np.genfromtxt(r'in_out/Approximation.dat')


for j in range(1,len(data_cum)-1):
    A_2[j] = data_cum[j,1]
    Mu1_2[j] = data_cum[j,2]
    Mu2_2[j] = data_cum[j,3]
    A_3[j] = data_cum[j,4]
    Mu1_3[j] = data_cum[j,5]
    Mu2_3[j] = data_cum[j,6]
    Mu3_3[j] = data_cum[j,7]
    A_4[j] = data_cum[j,8]
    Mu1_4[j] = data_cum[j,9]
    Mu2_4[j] = data_cum[j,10]
    Mu3_4[j] = data_cum[j,11]
    Mu4_4[j] = data_cum[j,12]

cut_off_ = np.zeros(Nu)
for i in range(Nu):
    if gauss_fun(m, A_4[i],Mu1_4[i], Mu2_4[i],Mu3_4[i],Mu4_4[i] ) == gauss_fun(m, A_4[i],Mu1_4[i], Mu2_4[i],Mu3_4[i],Mu4_4[i] ): # excludes 'nan'  
        cut_off_[i] = gauss_fun(m, A_4[i],Mu1_4[i], Mu2_4[i],Mu3_4[i],Mu4_4[i] )
        
k_max = list(cut_off_).index(np.max(cut_off_))


# You can vary 'accur' to obtain more precise results 
accur = 1000

# Choosing of k_0 and k_cut (cut off of the number of sectors) is heuristic

# Let's find the left cut off over sectors 
k_0 = 1

for j in range(int(Nu/10)):
    
    if  gauss_fun(m, A_4[j],Mu1_4[j], Mu2_4[j],Mu3_4[j], Mu4_4[j] ) < 10**(-12): # > 10**(-12) and A_2[i,j]!= 0 :
        k_0 = j              

# Let's find the right cut off over sectors  

p_4 = 0 

i = k_max

while cut_off_[k_max]/cut_off_[i] < accur and i < Nu - 1:
    
    p_4 += gauss_fun(m, A_4[k], Mu1_4[k], Mu2_4[k], Mu3_4[k], Mu4_4[k]) 

    i += 1

k_cut = i  

if 0 <= k_cut <= k_max + 3:
    for j in range(Nu):
        if  gauss_fun(m, A_4[j],Mu1_4[j], Mu2_4[j],Mu3_4[j],Mu4_4[j]) > gauss_fun(m, A_4[k_0],Mu1_4[k_0], Mu2_4[k_0],Mu3_4[k_0],Mu4_4[k_0]): # > 10**(-12) and A_2[i,j]!= 0 :
            k_cut = j
    
# The probability computation for different orders of approximation

p_2 = 0
p_3 = 0
p_4 = 0


for k in range(k_0, k_cut):
    p_2 += norm*gauss_fun(m, A_2[k], Mu1_2[k], Mu2_2[k]) 
    p_3 += norm*gauss_fun(m, A_3[k], Mu1_3[k], Mu2_3[k], Mu3_3[k])
    p_4 += norm*gauss_fun(m, A_4[k], Mu1_4[k], Mu2_4[k], Mu3_4[k], Mu4_4[k])



p_4_new = p_4

if p_4 > 10*p_2 or (p_4 == p_4)==False:
    while p_4_new > 10*p_2 and (p_4 == p_4)==False:
        k_cut -= 1 

        p_4_new = 0 
        for k in range(k_0, k_cut):
            p_4_new += norm*gauss_fun(m, A_4[k], Mu1_4[k], Mu2_4[k], Mu3_4[k], Mu4_4[k])
    p_4 = p_4_new
    
p_3 = 0 
for k in range(k_0, k_cut):
    p_3 += norm*gauss_fun(m, A_3[k], Mu1_3[k], Mu2_3[k], Mu3_3[k])   
# Warning! The 3rd order might give the 'nan' result  due to the function's kind.  

with open( "in_out/Result.dat", 'w') as ouf:
    ouf.write(  str('N\t') + str('K_0\t') + str('K_cut\t') +  str('Z\t') + str('Pexact          ') +  str( 'P2              ') + str('P3              ') + str( 'P4\n') ) 
    ouf.write( str(m) +  '\t' + str(k_0) + '\t' + str(k_cut) + '\t' + str( "%.3f" % Z(M_sub)) + '\t' + str("{:.3e}".format(exact_p))+ '\t' +  str("{:.3e}".format(p_2))+ '\t'+ str("{:.3e}".format(p_3)) +'\t'  + str("{:.3e}".format(p_4)) + '\n' ) 

