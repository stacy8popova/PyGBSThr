# -*- coding: utf-8 -*-
"""
@author: Popova A.S. 
"""


import random
import math
import numpy as np
import sys


#functions 

#input state
def A_ii(r, phi):
    
    y1 = math.tanh(r) 
    
    f = (-np.exp(1j*phi)/2)*y1
    
    return(f)

#interferometer matrix 
    

    
    
def full_ps_1(phi_list_1, ind,  i, m):
    
    V = np.eye(m, dtype=np.complex128)


    V[ind[0], ind[0]] = np.exp(1j*phi_list_1[i])
    V[ind[1], ind[1]] = 1 
    
    return(V)

def full_ps_2(phi_list_2, ind,  i, m):
    
    V = np.eye(m, dtype=np.complex128)


    V[ind[0], ind[0]] = 1
    V[ind[1], ind[1]] = np.exp(1j*phi_list_2[i]) 
    
    return(V)
    

def full_bm(h_list, i, ind, m):
    
    V = np.eye(m, dtype=np.complex128)

    V[ind[0], ind[0]] = np.cos(h_list[i])
    V[ind[0], ind[1]] = -np.sin(h_list[i])
    V[ind[1], ind[0]] = np.sin(h_list[i])
    V[ind[1], ind[1]] = np.cos(h_list[i])
                   
    return(V)

# Gaussian multi-mode matrix

def M_matrix(V, A, m):
    
    M = np.zeros((m,m), dtype=np.complex128)
    
    for k in range(m):
        for i in range(m):
            for j in range(m):
                M[i,j] += V[k,i]*V[k,j].conjugate()*A[k,k]
                
    return(M)


def fact(x):
    res = 1
    for i in range(int(x)):
        res *=(i+1.)
    return res


def C_nm(x,y):
    
    res = fact(int(x))/(fact(int(x-y))*fact(int(y)))
    return int(res) 

# Submatrix 

def red_mat(M_big, list_det): # [0, ..., n,k,l, ..., m-1] - numbers of clicked detectors     
       
    n = len(list_det)
    
    small_mat =  np.zeros((n, n), dtype = np.complex128)
    
    # [0,+,0,+] == [1,3]
    
    for i in range(n):
        for j in range(n):
            ind_i = list_det[i]
            ind_j = list_det[j]
            small_mat[i,j] = M_big[ind_i,ind_j]
    
    return(small_mat)


path = r'in_out'

m = input("\nNumber of modes: ")
m = int(m)

r = input("\nSqueezing parameter: ")

n = input("\nNumber of clicked detectors: ")
n = int(n)

if n > m:
    sys.exit('\nError: number of clicked detectors must be <= number of modes')
    
    
r_list = np.zeros(m) 

for i in range(m):
    if i < int(m/2):
        r_list[i] = r
    else:
        r_list[i] = 0

phi_list = np.zeros(m) 

for i in range(int(m/2)):
    phi_list[i] =  np.round(np.pi*random.random(),5) 
    
N_BS = 5*m**2 

N_PS = 2*N_BS



# Computing of a random interferometer matrix

phi_list_1 = [ np.round(2*np.pi*random.random(),5)  for i in range(N_BS)]
phi_list_2 = [ np.round(2*np.pi*random.random(),5)  for i in range(N_BS)]

psi_list = [ np.round(2*np.pi*random.random(),5)  for i in range(N_BS)]

ind_list = []

while len(ind_list) < (N_BS):
    k, j = random.randint(0,m-1), random.randint(0,m-1)
    if k != j:
        k_j = sorted([k, j])
        ind_list.append(k_j)

# Interferometer matrix

V =  np.eye(m, dtype=np.complex128)

for l in range(N_BS):
    
     V = V@full_ps_1(phi_list_1, ind_list[l],  l, m)@full_bm(psi_list, l, ind_list[l], m)@full_ps_2(phi_list_2, ind_list[l],  l, m)#
    
V = np.matrix(V, dtype=np.complex128)

# Generating of the Gaussian multi-mode matrix


A = np.zeros((m,m), dtype=np.complex128)

for i in range(m):
    A[i,i] = A_ii(r_list[i], phi_list[i])


M = M_matrix(V, A, m )

# Export data

with open( path + r'/Initial_state.dat', 'w') as ouf:
    
    ouf.write('N' + '\t'+str( 'r') +'\t'+str( 'phi') +'\t'+str( 'A_real')+'\t'+str('A_imag') +'\n')
        

    for k in range(m):
        
        ouf.write(str(k)+'\t'+str(r_list[k]) +'\t'+str(phi_list[k]) +'\t'+str((A_ii(r_list[k], phi_list[k]).real))+'\t'+str((A_ii(r_list[k], phi_list[k]).imag)))
        
        if k < (m+1):
    
            ouf.write('\n')
        
with open( path +'/Parmeters_of_interferometer.dat', 'w') as ouf:

    ouf.write('# N_in = ' + str(m)+ '\t' + str( 'N_bs = ')  + str(N_BS) + '\t' + str( 'N_ps = ') + str(N_PS) +'\n')
    ouf.write('[n1, n2]'  + '\t'+ str( 'phi_1')  + '\t'+ str( 'phi_2')  + '\t' +  str( 'alpha') + '\n')


    for z in range(N_BS):
        ouf.write(str(ind_list[z][0]) +'\t' +str(ind_list[z][1]) +'\t'+ str( phi_list_1[z] )  + '\t'+  str(  phi_list_2[z] ) + '\t' + str(  psi_list[z]) + '\n')

with open( path + "/GBS_matrix.dat", 'w') as ouf:
    
    for k in range(m):
        for j in range(m):
            ouf.write( str(M[k,j].real) + '\t' + str(M[k,j].imag) +'\t' ) 
            
        if k < (m+1):
            ouf.write('\n')
 
# n - number detectors cliked 

# N - total number of samples 

N = 1

list_det = []

for i in range(N):
    list_det.append(sorted(random.sample(range(m), n)))
    
for sp in range(N):

    sample = list(list_det[sp][:])

    R = red_mat(M, sample)

    
    with open( path + r"/Submatrix.dat", 'w') as ouf:  #Submatrix_"+str(sp)+".dat" if N>1
        
        for k in range(n):
            for j in range(n):
                ouf.write( str(R[k,j].real) + '\t' + str(R[k,j].imag) +'\t' ) 

            if k < (n+1):
                ouf.write('\n')
                
            
with open( path + r"/Sample.dat", 'w') as ouf: 
    
    ouf.write('The string of clicked detectors numbers' + '\n')

    for k in range(N):
        ouf.write(str(list_det[k][:]) + '\t') 
        if k < (N+1):
            ouf.write('\n')
