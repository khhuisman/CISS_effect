#!/usr/bin/env python
# coding: utf-8



# 1. Representation is made for any l
# 2. Onsite SOC $\vec{L}\cdot \vec{S}$ is calculated in cubic harmonic basis for $l=2$


# Output Module:
# - Onsite atomic spin-orbit coupling Hamiltonian in the cubic harmonic basis for l =1 or l=2.


# # Packages

# In[38]:


import numpy as np
import math
from sympy import *

#pritty print the np.array: B
def matprint(B):
    pprint(Matrix(B.tolist()))

#pritty print the list: B
def listprint(B):
    pprint(Matrix(B))


# In[39]:


import tinyarray

# define Pauli-matrices for convenience
sigma_0 = tinyarray.array([[1, 0], [0, 1]])
sigma_x = tinyarray.array([[0, 1], [1, 0]])
sigma_y = tinyarray.array([[0, -1j], [1j, 0]])
sigma_z = tinyarray.array([[1, 0], [0, -1]])


# # Representation $J_z, J_x, J_y$ 

# in $Y_{lm}$ basis ($\hbar =1 $): 
# \begin{align}
# J_z|l,m\rangle = m |l,m\rangle \\
# J^{\pm} = \sqrt{j(j+1) - m(m\pm 1) } |l,m \pm 1\rangle
# \end{align}

# In[40]:


    

def angular_momentum(jp=2):
    j = jp
    deg = int(2*j+1)

    J_p = np.zeros((deg,deg))


    list_mz = -1*np.arange(-j,j+1,1)
    J_z = np.array(np.diag(list_mz),dtype = complex)

    for m in range(0,deg-1):
        mp = m-j
        J_p[m,m+1] = np.sqrt(j*(j+1)- mp*(mp+1))

    J_m = np.transpose(J_p)
    J_x = 0.5*(J_p + J_m)
    J_y = -0.5*1j*(J_p - J_m)
    return [J_x,J_y,J_z]


# # Calculate $L\cdot S$  for $l =2$

# ####  Basis Transformation: $Y_{lm} \rightarrow  d_i$ (cubic harmonic basis)

# Basis transfomation matrix $U$: obtained from https://en.wikipedia.org/wiki/Cubic_harmonic

# In[51]:

def LdotS(l=2):
    
    L_vec = angular_momentum(l)
    # Caculate L*S
    LS_x = np.kron(L_vec[0],sigma_x)
    LS_y = np.kron(L_vec[1],sigma_y)
    LS_z = np.kron(L_vec[2],sigma_z)
    
    
    # The inner product of L.S in the Y_lm basis
    LS_sum_ylm = np.multiply(np.add(np.add(LS_x,LS_y),LS_z),0.5)
    
    if l == 1:
        #Basis transformation matrix: Y_lm -> [px,py,pz]
        Ul = [   [-0.5*np.sqrt(2), 0, 0.5*np.sqrt(2) ],
                 [(-1j)*0.5*np.sqrt(2), 0, (-1j)*0.5*np.sqrt(2)],
                 [0, 1, 0]
             ]
        #Spin DOF arre added:
        Uj = np.kron(Ul,sigma_0)
        #Hermitian conjugate
        Ujv = np.transpose(np.conjugate(Uj))


        # L*S term for p orbitals
        LS_cubic = np.round(np.dot(Uj,np.dot(LS_sum_ylm,Ujv)),15)

        return LS_cubic
    
    
    if l == 2:
        #Basis transformation matrix: Y_lm -> [d_xy,dx2my2,dxz,dyz,dz2]
        Ul = np.array([ 
                      [1j*0.5*np.sqrt(2), 0, 0, 0, -1j*0.5*np.sqrt(2)], 
                      [1*0.5*np.sqrt(2), 0, 0, 0, 1*0.5*np.sqrt(2)], 
                      [0, -1*0.5*np.sqrt(2), 0, 1*0.5*np.sqrt(2), 0], 
                      [0, -1j*0.5*np.sqrt(2), 0, -1j*0.5*np.sqrt(2), 0], 
                      [0, 0, 1, 0, 0] 
                  ])
        #Spin DOF arre added:
        Uj = np.kron(Ul,sigma_0)
        #Hermitian conjugate
        Ujv = np.transpose(np.conjugate(Uj))


        # L*S term for d orbitals
        LS_cubic = np.round(np.dot(Uj,np.dot(LS_sum_ylm,Ujv)),15)

        return LS_cubic
    
    if l > 2:
        print("Code not built for l = {}".format(l))
    if l == 0:
        print("No spin orbit coupling for l = {}".format(l))
        
 

   
    
    
    

 

