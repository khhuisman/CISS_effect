#!/usr/bin/env python
# coding: utf-8


##### Purpose of notebook
# With this notebook we can generate hopping matrices between the orbitals two gold atoms on different sites.
# 

# In[3]:


import numpy as np
from sympy import *


# In[4]:


sqrt3 = np.sqrt(3)
sqrt2 = np.sqrt(2)


# Hopping between s,p orbitals

# In[3]:


def Hsp(l,m,n,Vsss,Vsps,Vppp,Vpps):
    
    
    H_sp = [[Vsss, l*Vsps, m*Vsps, n*Vsps],
     [-(l*Vsps), (1 - l**2)*Vppp + l**2*Vpps, -(l*m*Vppp) + l*m*Vpps, -(l*n*Vppp) + l*n*Vpps],
     [-(m*Vsps), -(l*m*Vppp) + l*m*Vpps, (1 - m**2)*Vppp + m**2*Vpps, -(m*n*Vppp) + m*n*Vpps],
     [-(n*Vsps), -(l*n*Vppp) + l*n*Vpps, -(m*n*Vppp) + m*n*Vpps, (1 - n**2)*Vppp + n**2*Vpps]]
    
    
    
    return H_sp
        

# s,d orbitals

# In[4]:


def Hsd(l,m,n,Vsds):
    Hsd = [
            sqrt3*l*m*Vsds, 
            (sqrt3*(l**2 - m**2)*Vsds)/2, 
           sqrt3*l*n*Vsds, sqrt3*m*n*Vsds, 
           ((-l**2 - m**2)/2 + n**2)*Vsds
          ]
    return Hsd




# p,d orbitals

# In[8]:


def Hpd(l,m,n,Vpds,Vpdp):
    
    H_pd = [
            [
                (1 - 2*l**2)*m*Vpdp + sqrt3*l**2*m*Vpds, 
             l*(1 - l**2 + m**2)*Vpdp + (sqrt3*l*(l**2 - m**2)*Vpds)/2,
             (1 - 2*l**2)*n*Vpdp + sqrt3*l**2*n*Vpds, -2*l*m*n*Vpdp + sqrt3*l*m*n*Vpds, 
             -(sqrt3*l*n**2*Vpdp) + l*((-l**2 - m**2)/2 + n**2)*Vpds
            ],
            [
             l*(1 - 2*m**2)*Vpdp + sqrt3*l*m**2*Vpds, 
             -(m*(1 + l**2 - m**2)*Vpdp) + (sqrt3*m*(l**2 - m**2)*Vpds)/2, 
             -2*l*m*n*Vpdp + sqrt3*l*m*n*Vpds, 
             (1 - 2*m**2)*n*Vpdp + sqrt3*m**2*n*Vpds, 
             -(sqrt3*m*n**2*Vpdp) + m*((-l**2 - m**2)/2 + n**2)*Vpds
            ],
            [
                -2*l*m*n*Vpdp + sqrt3*l*m*n*Vpds,
              -((l**2 - m**2)*n*Vpdp) + (sqrt3*(l**2 - m**2)*n*Vpds)/2, 
                 l*(1 - 2*n**2)*Vpdp + sqrt3*l*n**2*Vpds, 
                 m*(1 - 2*n**2)*Vpdp + sqrt3*m*n**2*Vpds, 
                 sqrt3*(l**2 + m**2)*n*Vpdp + n*((-l**2 - m**2)/2 + n**2)*Vpds
            ]
        ]
    
    return H_pd


# In[ ]:



# d,d orbitals

# In[13]:


def Hdd(l,m,n,  Vdds = -0.04971, Vddp = 0.02624, Vddl = -0.00457):
    
    E_dd = [
            [(l**2*m**2 + n**2)*Vddl + (l**2 + m**2 - 4*l**2*m**2)*Vddp + 3*l**2*m**2*Vdds, 
             (l*m*(l**2 - m**2)*Vddl)/2 + 2*l*m*(-l**2 + m**2)*Vddp + (3*l*m*(l**2 - m**2)*Vdds)/2, 
             (-1 + l**2)*m*n*Vddl + (1 - 4*l**2)*m*n*Vddp + 3*l**2*m*n*Vdds, l*(-1 + m**2)*n*Vddl + l*(1 - 4*m**2)*n*Vddp + 3*l*m**2*n*Vdds, 
             sqrt3*((l*m*(1 + n**2)*Vddl)/2 - 2*l*m*n**2*Vddp + l*m*((-l**2 - m**2)/2 + n**2)*Vdds)
            ],
             [
              (l*m*(l**2 - m**2)*Vddl)/2 + 2*l*m*(-l**2 + m**2)*Vddp + (3*l*m*(l**2 - m**2)*Vdds)/2, 
              ((l**2 - m**2)**2/4 + n**2)*Vddl + (l**2 + m**2 - (l**2 - m**2)**2)*Vddp + (3*(l**2 - m**2)**2*Vdds)/4, 
                 -(l*(1 + (-l**2 + m**2)/2)*n*Vddl) + l*(1 - 2*(l**2 - m**2))*n*Vddp + (3*l*(l**2 - m**2)*n*Vdds)/2, 
                 m*(1 + (l**2 - m**2)/2)*n*Vddl - m*(1 + 2*(l**2 - m**2))*n*Vddp + (3*m*(l**2 - m**2)*n*Vdds)/2, 
                 sqrt3*(((l**2 - m**2)*(1 + n**2)*Vddl)/4 + (-l**2 + m**2)*n**2*Vddp + ((l**2 - m**2)*((-l**2 - m**2)/2 + n**2)*Vdds)/2)
             ],
             [
            (-1 + l**2)*m*n*Vddl + (1 - 4*l**2)*m*n*Vddp + 3*l**2*m*n*Vdds,
              -(l*(1 + (-l**2 + m**2)/2)*n*Vddl) + l*(1 - 2*(l**2 - m**2))*n*Vddp + (3*l*(l**2 - m**2)*n*Vdds)/2, 
              (m**2 + l**2*n**2)*Vddl + (l**2 + n**2 - 4*l**2*n**2)*Vddp + 3*l**2*n**2*Vdds, l*m*(-1 + n**2)*Vddl + l*m*(1 - 4*n**2)*Vddp + 3*l*m*n**2*Vdds, 
              sqrt3*(-(l*(l**2 + m**2)*n*Vddl)/2 + l*n*(l**2 + m**2 - n**2)*Vddp + l*n*((-l**2 - m**2)/2 + n**2)*Vdds)
             ],
             [l*(-1 + m**2)*n*Vddl + l*(1 - 4*m**2)*n*Vddp + 3*l*m**2*n*Vdds, 
              m*(1 + (l**2 - m**2)/2)*n*Vddl - m*(1 + 2*(l**2 - m**2))*n*Vddp + (3*m*(l**2 - m**2)*n*Vdds)/2, 
              l*m*(-1 + n**2)*Vddl + l*m*(1 - 4*n**2)*Vddp + 3*l*m*n**2*Vdds, 
              (l**2 + m**2*n**2)*Vddl + (m**2 + n**2 - 4*m**2*n**2)*Vddp + 3*m**2*n**2*Vdds, 
              sqrt3*(-(m*(l**2 + m**2)*n*Vddl)/2 + m*n*(l**2 + m**2 - n**2)*Vddp + m*n*((-l**2 - m**2)/2 + n**2)*Vdds)
             ],
             [sqrt3*((l*m*(1 + n**2)*Vddl)/2 - 2*l*m*n**2*Vddp + l*m*((-l**2 - m**2)/2 + n**2)*Vdds), 
              sqrt3*(((l**2 - m**2)*(1 + n**2)*Vddl)/4 + (-l**2 + m**2)*n**2*Vddp + ((l**2 - m**2)*((-l**2 - m**2)/2 + n**2)*Vdds)/2), 
              sqrt3*(-(l*(l**2 + m**2)*n*Vddl)/2 + l*n*(l**2 + m**2 - n**2)*Vddp + l*n*((-l**2 - m**2)/2 + n**2)*Vdds), 
              sqrt3*(-(m*(l**2 + m**2)*n*Vddl)/2 + m*n*(l**2 + m**2 - n**2)*Vddp + m*n*((-l**2 - m**2)/2 + n**2)*Vdds), 
              (3*(l**2 + m**2)**2*Vddl)/4 + 3*(l**2 + m**2)*n**2*Vddp + ((-l**2 - m**2)/2 + n**2)**2*Vdds
             ]
            ]
    return E_dd


# In[ ]:




# Full Hopping Hamiltonian

# In[21]:


def H_hop_spd(l,m,n,  Vsss = -0.0668,Vsps = 0.09721,
                  Vpps = 0.17866, Vppp = -0.01645, 
                  Vdds = -0.04971, Vddp = 0.02624, Vddl = -0.00457, 
                  Vsds = -0.04722, 
                  Vpdp = 0.01896, Vpds = -0.06399):
    
    
    Hhop = [[Vsss, l*Vsps, m*Vsps, n*Vsps, sqrt3*l*m*Vsds, (sqrt3*(l**2 - m**2)*Vsds)/2, sqrt3*l*n*Vsds, sqrt3*m*n*Vsds, ((-l**2 - m**2)/2 + n**2)*Vsds],
             [-(l*Vsps), (1 - l**2)*Vppp + l**2*Vpps, -(l*m*Vppp) + l*m*Vpps, -(l*n*Vppp) + l*n*Vpps, (1 - 2*l**2)*m*Vpdp + sqrt3*l**2*m*Vpds, l*(1 - l**2 + m**2)*Vpdp + (sqrt3*l*(l**2 - m**2)*Vpds)/2, (1 - 2*l**2)*n*Vpdp + sqrt3*l**2*n*Vpds, -2*l*m*n*Vpdp + sqrt3*l*m*n*Vpds, -(sqrt3*l*n**2*Vpdp) + l*((-l**2 - m**2)/2 + n**2)*Vpds],
             [-(m*Vsps), -(l*m*Vppp) + l*m*Vpps, (1 - m**2)*Vppp + m**2*Vpps, -(m*n*Vppp) + m*n*Vpps, l*(1 - 2*m**2)*Vpdp + sqrt3*l*m**2*Vpds, -(m*(1 + l**2 - m**2)*Vpdp) + (sqrt3*m*(l**2 - m**2)*Vpds)/2, -2*l*m*n*Vpdp + sqrt3*l*m*n*Vpds, (1 - 2*m**2)*n*Vpdp + sqrt3*m**2*n*Vpds, -(sqrt3*m*n**2*Vpdp) + m*((-l**2 - m**2)/2 + n**2)*Vpds],
             [-(n*Vsps), -(l*n*Vppp) + l*n*Vpps, -(m*n*Vppp) + m*n*Vpps, (1 - n**2)*Vppp + n**2*Vpps, -2*l*m*n*Vpdp + sqrt3*l*m*n*Vpds, (-l**2 + m**2)*n*Vpdp + (sqrt3*(l**2 - m**2)*n*Vpds)/2, l*(1 - 2*n**2)*Vpdp + sqrt3*l*n**2*Vpds, m*(1 - 2*n**2)*Vpdp + sqrt3*m*n**2*Vpds, sqrt3*(l**2 + m**2)*n*Vpdp + n*((-l**2 - m**2)/2 + n**2)*Vpds],
             [sqrt3*l*m*Vsds, (-1 + 2*l**2)*m*Vpdp - sqrt3*l**2*m*Vpds, -(l*(1 - 2*m**2)*Vpdp) - sqrt3*l*m**2*Vpds, 2*l*m*n*Vpdp - sqrt3*l*m*n*Vpds, (l**2*m**2 + n**2)*Vddl + (l**2 + m**2 - 4*l**2*m**2)*Vddp + 3*l**2*m**2*Vdds, (l*m*(l**2 - m**2)*Vddl)/2 + 2*l*m*(-l**2 + m**2)*Vddp + (3*l*m*(l**2 - m**2)*Vdds)/2, (-1 + l**2)*m*n*Vddl + (1 - 4*l**2)*m*n*Vddp + 3*l**2*m*n*Vdds, l*(-1 + m**2)*n*Vddl + l*(1 - 4*m**2)*n*Vddp + 3*l*m**2*n*Vdds, sqrt3*((l*m*(1 + n**2)*Vddl)/2 - 2*l*m*n**2*Vddp + l*m*((-l**2 - m**2)/2 + n**2)*Vdds)],
             [(sqrt3*(l**2 - m**2)*Vsds)/2, -(l*(1 - l**2 + m**2)*Vpdp) - (sqrt3*l*(l**2 - m**2)*Vpds)/2, m*(1 + l**2 - m**2)*Vpdp - (sqrt3*m*(l**2 - m**2)*Vpds)/2, (l**2 - m**2)*n*Vpdp - (sqrt3*(l**2 - m**2)*n*Vpds)/2, (l*m*(l**2 - m**2)*Vddl)/2 + 2*l*m*(-l**2 + m**2)*Vddp + (3*l*m*(l**2 - m**2)*Vdds)/2, ((l**2 - m**2)**2/4 + n**2)*Vddl + (l**2 + m**2 - (l**2 - m**2)**2)*Vddp + (3*(l**2 - m**2)**2*Vdds)/4, -(l*(1 + (-l**2 + m**2)/2)*n*Vddl) + l*(1 - 2*(l**2 - m**2))*n*Vddp + (3*l*(l**2 - m**2)*n*Vdds)/2, m*(1 + (l**2 - m**2)/2)*n*Vddl - m*(1 + 2*(l**2 - m**2))*n*Vddp + (3*m*(l**2 - m**2)*n*Vdds)/2, sqrt3*(((l**2 - m**2)*(1 + n**2)*Vddl)/4 + (-l**2 + m**2)*n**2*Vddp + ((l**2 - m**2)*((-l**2 - m**2)/2 + n**2)*Vdds)/2)],
             [sqrt3*l*n*Vsds, (-1 + 2*l**2)*n*Vpdp - sqrt3*l**2*n*Vpds, 2*l*m*n*Vpdp - sqrt3*l*m*n*Vpds, -(l*(1 - 2*n**2)*Vpdp) - sqrt3*l*n**2*Vpds, (-1 + l**2)*m*n*Vddl + (1 - 4*l**2)*m*n*Vddp + 3*l**2*m*n*Vdds, -(l*(1 + (-l**2 + m**2)/2)*n*Vddl) + l*(1 - 2*(l**2 - m**2))*n*Vddp + (3*l*(l**2 - m**2)*n*Vdds)/2, (m**2 + l**2*n**2)*Vddl + (l**2 + n**2 - 4*l**2*n**2)*Vddp + 3*l**2*n**2*Vdds, l*m*(-1 + n**2)*Vddl + l*m*(1 - 4*n**2)*Vddp + 3*l*m*n**2*Vdds, sqrt3*(-(l*(l**2 + m**2)*n*Vddl)/2 + l*n*(l**2 + m**2 - n**2)*Vddp + l*n*((-l**2 - m**2)/2 + n**2)*Vdds)],
             [sqrt3*m*n*Vsds, 2*l*m*n*Vpdp - sqrt3*l*m*n*Vpds, (-1 + 2*m**2)*n*Vpdp - sqrt3*m**2*n*Vpds, -(m*(1 - 2*n**2)*Vpdp) - sqrt3*m*n**2*Vpds, l*(-1 + m**2)*n*Vddl + l*(1 - 4*m**2)*n*Vddp + 3*l*m**2*n*Vdds, m*(1 + (l**2 - m**2)/2)*n*Vddl - m*(1 + 2*(l**2 - m**2))*n*Vddp + (3*m*(l**2 - m**2)*n*Vdds)/2, l*m*(-1 + n**2)*Vddl + l*m*(1 - 4*n**2)*Vddp + 3*l*m*n**2*Vdds, (l**2 + m**2*n**2)*Vddl + (m**2 + n**2 - 4*m**2*n**2)*Vddp + 3*m**2*n**2*Vdds, sqrt3*(-(m*(l**2 + m**2)*n*Vddl)/2 + m*n*(l**2 + m**2 - n**2)*Vddp + m*n*((-l**2 - m**2)/2 + n**2)*Vdds)],
             [((-l**2 - m**2)/2 + n**2)*Vsds, sqrt3*l*n**2*Vpdp - l*((-l**2 - m**2)/2 + n**2)*Vpds, sqrt3*m*n**2*Vpdp - m*((-l**2 - m**2)/2 + n**2)*Vpds, -(sqrt3*(l**2 + m**2)*n*Vpdp) - n*((-l**2 - m**2)/2 + n**2)*Vpds, sqrt3*((l*m*(1 + n**2)*Vddl)/2 - 2*l*m*n**2*Vddp + l*m*((-l**2 - m**2)/2 + n**2)*Vdds), sqrt3*(((l**2 - m**2)*(1 + n**2)*Vddl)/4 + (-l**2 + m**2)*n**2*Vddp + ((l**2 - m**2)*((-l**2 - m**2)/2 + n**2)*Vdds)/2), sqrt3*(-(l*(l**2 + m**2)*n*Vddl)/2 + l*n*(l**2 + m**2 - n**2)*Vddp + l*n*((-l**2 - m**2)/2 + n**2)*Vdds), sqrt3*(-(m*(l**2 + m**2)*n*Vddl)/2 + m*n*(l**2 + m**2 - n**2)*Vddp + m*n*((-l**2 - m**2)/2 + n**2)*Vdds), (3*(l**2 + m**2)**2*Vddl)/4 + 3*(l**2 + m**2)*n**2*Vddp + ((-l**2 - m**2)/2 + n**2)**2*Vdds]]
    return Hhop


# In[29]:


print("H_hop_spd(l,m,n,  Vsss,Vsps, Vpps, Vppp, Vdds, Vddp, Vddl, Vsds, Vpdp, Vpds)" )


# In[164]:


def H_Au_NN(l,m,n):
    
    """
    Input:
    - direct cosines l,m,n in the x,y,z direction for NN hopping on the fcc lattice. 
    - NN hopping on the fcc lattice means that the absoluate value of l,m,n can only take the value: 1/sqrt(2)
    - Slater Koster Tight Binding elements of gold taken from Pappaconstantopoulos.
    Ouput:
    - NN hopping matrix of gold.
    """

   

    HAuNN = H_hop_spd(l,m,n,  Vsss = -0.0668,Vsps = 0.09721,
                  Vpps = 0.17866, Vppp = -0.01645, 
                  Vdds = -0.04971, Vddp = 0.02624, Vddl = -0.00457, 
                  Vsds = -0.04722, 
                  Vpdp = 0.01896, Vpds = -0.06399 )
    
    
    return HAuNN
    
    


# In[165]:


def H_Au_NNN(l,m,n):
    """
    Input:
    - direct cosines l,m,n in the x,y,z direction for NNN hopping on the fcc lattice. 
    - This means that the absoluate value of l,m,n can only take the values: 1.
    - Slater Koster Tight Binding elements  of gold taken from Pappaconstantopoulos.
    Ouput:
    - NNN hopping matrix of gold.
    """



    HAu_NNN = H_hop_spd(l,m,n,  Vsss = 0.00277,  Vsps = 0.00261,
                          Vpps = 0.03707, Vppp = -0.01025, 
                          Vdds = -0.00305, Vddp = 0.0024, Vddl = -0.00057, 
                          Vsds = -0.00784, 
                          Vpds = -0.00762, Vpdp = 0.0047 )
    
    
    return HAu_NNN


# In[151]:


def generate_hoppings():
    
    
    
    H110 = np.array(H_Au_NN(1/sqrt2,1/sqrt2,0))
    H1m10 = np.array(H_Au_NN(1/sqrt2,-1/sqrt2,0))
    
    H101 = np.array(H_Au_NN(1/sqrt2,0,1/sqrt2))
    H10m1 = np.array(H_Au_NN(1/sqrt2,0,-1/sqrt2))
    
    H011 = np.array(H_Au_NN(0,1/sqrt2,1/sqrt2))
    H01m1 = np.array(H_Au_NN(0,1/sqrt2,-1/sqrt2))
    
    
    H001 = np.array(H_Au_NNN(0,0,1))
    H010 = np.array(H_Au_NNN(0,1,0))
    H100 = np.array(H_Au_NNN(1,0,0) )
    
    return [    H110, H1m10, 
                H101, H10m1,
                H011, H01m1,
                H001,H010,H100]




    

