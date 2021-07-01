#!/usr/bin/env python
# coding: utf-8

# # Helicene - Cuniberti et al.

#######################################################################

import sympy

import kwant
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


# For matrix support
import tinyarray

# define Pauli-matrices for convenience
sigma_0 = tinyarray.array([[1, 0], [0, 1]])
sigma_x = tinyarray.array([[0, 1], [1, 0]])
sigma_y = tinyarray.array([[0, -1j], [1j, 0]])
sigma_z = tinyarray.array([[1, 0], [0, -1]])
sigma_n = tinyarray.array([[0, 0], [0, 0]])


# # Helicene 3D

# 1. n: number of benzene rings
# 2. xi_p: Atomic Spin Orbit coupling of Carbon
# 3. tvec = number that multiplies vec (x,y) = (1,1) to translate the molecule in xy plane.
# 4. E_z = z component of External electric field.
# 5. r: radius of the inner helix.
# 6. b = pitch of the helices
# 7. chirality_left: Boolean, True/False corresponds to left-,right-handed version of hel.

# ### Hoppings

# In[20]:
# In[4]:
import sys
sys.path.insert(0, '/Users/khhuisman/Documents/Jupyter_notebooks/py_files')



import Helicene_Hamiltonian

def attachement_position(I,M,E1,E2):
    """
    Input:
    - labels for the sites on the inner (I), middle (M) and exterior (E1,E2) helix.
    Ouput:
    - x,y coordinates on the cubic lattice corresponding to the attachement position of the probe.
    """
    if E1 !=0 and  E2 ==0 and M ==0:
        attposxE = 1
        attposyE = 4
    
    if E2 !=0 and E1 ==0 and M ==0:
        attposxE = 2
        attposyE = 5
        
    if E1 ==0 and  E2 ==0:
        attposxE = 0
        attposyE = 0
        
    if E1 ==0 and  E2 ==0:
        attposxE = 0
        attposyE = 0
          
    if E2 !=0 and E1 !=0:
        print('attposE1,attposE2 must be different,attposE set to 0 by default')
        attposxE = 0
        attposyE = 0
        
    if E1 !=0 and  E2 == 0 and M !=0:
        
        print('Invalid Configuration,set M or E1 to zero')
       
        
    if E1 ==0 and  E2 !=0 and M !=0:
        print('Invalid Configuration,set M or E2 to zero')
        
    return 4 + 3*I + attposxE ,4 + 3*I + attposyE + 2*M





def make_helicene_buttiker(nrings,xi_p = 0.006,tvec=0,E_z=0,
                  r=1.4,b=3.6,zstart=0,chirality_left =True,
                 es=-18 , ep_sigma=-18, ep_pi=-10.5,u0=0,t0=1,
                           attposI=1,attposM=0,
                           attposE1=0,attposE2=0):
    
    
    """
    Input:
    - systems parameters of helicene
    - system parameters of the Büttiker probe
    
    Output:
    - kwant system of helicene with a Büttiker probe attached to one site of the helicene helix. 
    """
    
    L= 3*nrings + 1
    
    hel = kwant.Builder()
    lat = kwant.lattice.cubic(a=1,norbs =8)
    
    sigma_08 = np.identity(8)
    sigma_z8 = np.kron(np.identity(4),sigma_z) 
    
    
    # Onsite Hamiltonian
    H_onsite = Helicene_Hamiltonian.H_onsite_Carbon(xi_p, E_z, es, ep_sigma, ep_pi)
    
    
    # Set hoppings
    h_list = Helicene_Hamiltonian.H_hop_list(r,b,chirality_left)
    H_IIs = h_list[0]
    H_MEs = h_list[1]
    H_EMs = h_list[2]
    H_EEs = h_list[3]
    H_IMs = h_list[4]
    
    
#     #Create sites of:
        # Inner helix
    hel[(lat(i,i,zstart) for i in range(1+tvec,L+1+tvec,3))]   = H_onsite
        #Middle helix
    hel[(lat(i,2+i,zstart) for i in range(1+tvec,L+1+tvec,3))] = H_onsite
        # Exterior helix
    hel[(lat(i+j,3+i+j,zstart)   for i in range(2+tvec,L+1+tvec,3) 
                                 for j in range(0,2))] = H_onsite



    #Hoppings Inner Helix
    for i in range(1+tvec,L+1+tvec,3):
        if i>= 4+tvec:
            hel[lat(i-3,i-3,zstart),lat(i,i,zstart)]= H_IIs
    
    #Hopping Exterior Helix
    for i in range(1+tvec,L+tvec,3):
        hel[lat(i,2+i,zstart),lat(i+1,3+i+1,zstart)] = H_MEs
    
    for i in range(3+tvec,L+tvec,3):
        hel[lat(i,3+i,zstart),lat(i+1,2+i+1,zstart)] = H_EMs
    
    for i in range(2+tvec,L+tvec,3):
        hel[lat(i,3+i,zstart),lat(i+1,3+i+1,zstart)] = H_EEs
        
        
    #Hopping Middle-Interior Helix
    for i in range(1+tvec,L+1+tvec,3):
            hel[lat(i,i,zstart),lat(i,2+i,zstart)]=  H_IMs
            
            
            
    sigma_z8 = np.kron(np.identity(4),sigma_z)  
    
    ##Create sites of Buttiker probe
    H_buttiker =  u0*np.identity(8)     
    
    sym_buttiker = kwant.TranslationalSymmetry((0, 0,1))

    # attachement position coordinates of buttiker probe
    x_att,y_att = attachement_position(attposI,attposM,
                           attposE1,attposE2)
        

    lead_buttiker3 = kwant.Builder(sym_buttiker,conservation_law = sigma_z8)
    lead_buttiker3[lat(tvec + x_att,
                       tvec + y_att,
                       zstart + 2) ] = H_buttiker
    lead_buttiker3[kwant.builder.HoppingKind((0,0,1), lat, lat)] = -1*np.kron(np.diag([t0,t0,t0,t0]),np.identity(2))
    hel.attach_lead(lead_buttiker3)
    
    return hel




