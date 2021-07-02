#!/usr/bin/env python
# coding: utf-8

# # Import Packages

# In[3]:


import kwant
import math
import cmath
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sympy import *

import time

import numpy as np


# For matrix support
import tinyarray

# define Pauli-matrices for convenience
sigma_0 = tinyarray.array([[1, 0], [0, 1]])
sigma_x = tinyarray.array([[0, 1], [1, 0]])
sigma_y = tinyarray.array([[0, -1j], [1j, 0]])
sigma_z = tinyarray.array([[1, 0], [0, -1]])
sigma_n = tinyarray.array([[0, 0], [0, 0]])


# # Add path

# In[4]:


import sys
sys.path.insert(0, '/Users/khhuisman/Documents/Jupyter_notebooks/py_files')


# # Helicene 3D

# In[6]:


import HeliceneSystem


# 1. n: number of benzene rings
# 2. xi_p: Atomic Spin Orbit coupling of Carbon
# 3. tvec = number that multiplies vec (x,y) = (1,1) to translate the molecule in xy plane.
# 4. E_z = z component of External electric field.
# 5. r: radius of the inner helix.
# 6. b = pitch of the helices
# 7. chirality_left: Boolean, True/False corresponds to left-,right-handed version of hel.
# 8. es , ep_sigma, ep_pi: onsite energies of 2s,2$p_{x,y}$ and 2$p_z$ orbitals respectively of Carbon

# ### Handy Functions

# In[9]:


def nplist(vlist):
    varray = np.array(vlist)
    return varray

#convert units of matrix a from Rydberg to Electronvolt
def r2v(a):
    c_2eV = 13.605698066
    mat9x9= np.array(np.multiply(c_2eV,a),dtype =complex)
#     matrix =np.kron(mat9x9,sigma_0)
    return mat9x9

#convert units of matrix a from Rydberg to Electronvolt & add spin DOF

def r2v_soc(a):
    c_2eV = 13.605698066
    mat9x9= np.array(np.multiply(c_2eV,a),dtype =complex)
    matrix18x18 =np.kron(mat9x9,sigma_0)
    return matrix18x18

#convert units of matrix a from Hartree to Electronvolt & add spin DOF


def hart2ev_soc(a):
    c_2eV = 27.211399
    mat9x9= np.array(np.multiply(c_2eV,a),dtype =complex)
    matrix18x18 =np.kron(mat9x9,sigma_0)
    return matrix18x18

#convert units of matrix a from Hartree to Electronvolt


def hart2ev(a):
    c_2eV = 27.211399
    mat9x9= np.array(np.multiply(c_2eV,a),dtype =complex)
    return mat9x9


# Direct sum of matrix v1 & v2
def dis(v1,v2): 
    ab = np.zeros( np.add(v1.shape,v2.shape),dtype = complex )
    ab[:v1.shape[0],:v1.shape[1]]= v1
    ab[v1.shape[0]:,v1.shape[1]:]= v2
    return ab


# Fermi energy of gold in Electron volt/ Rydberg
def fermi_gold(eV= False):
    if eV == False:
        # energy in rydberg
        energy_fermi = 0.5380
        return energy_fermi
    
    if eV == True:
        # energy in rydberg
        energy_fermi = 0.5380*13.605698066
        return energy_fermi


# # Gold Leads

#  1. function that makes the lead
#  2. cut = 0,1,2. It ensures that we have ABC, CAB or BCA stacking of gold

# In[10]:


import GoldNormal


# 1. nlayers = number of gold sites with SOC.
# 2. Txz,Tyz = number of sites in xz,yz direction.
# 3. d_L     = distance between the leads. (simply moves around leads s.t.                               it does not intersect with the molecule).
# 3. xi_d   = onsite SOC parameter of 5d_orbitals of gold.
# 4. check_conjugate = boolean, takes conjugate of onsite Hamiltonian if True.
# 5. delta_e = energy that shifts all onsite energies of 6s,6p,5d orbitals.

# # Coupling Au-S-C

# In[11]:


import Coupling_Au_S_C


# # Add system

# In[12]:



# Color the different lattices
# Sulfur : Bright yellow
# Gold   : Goldlike yellow
# Helicene: inner,middle , exterior helix are coloured blue,red,green respectively.
def family_color(site):
    
    
    lat_fcc = kwant.lattice.general([( 0.5, 0.5, 0), 
                                 (0.5, 0, 0.5),
                                 ( 0,0.5,0.5)
                              ], norbs=18)
    lat_sulfur = kwant.lattice.cubic(a=1,norbs =6)


    site_family = site.family
    
    if site_family == lat_fcc:
        return 'gold'
    if site_family == lat_sulfur:
        return 'yellow'
    else:
        n1,n2,n3 = site[1]

        #Inner helix
        if n2 == n1:
            return 'blue'

        #Middle helix
        if n2 == 2+n1:
            return 'red'

        #Exterior Helix
        if n2 == 3+n1:
            return 'green'


# In[13]:


def make_system(nlayers=3, xi_d=0.6,xi_p = 0.006,
                Txz=3,Tyz=3, d_L= 50,
                 txz_left_1=0,tyz_left_1=0,
                txz_left_2=1,tyz_left_2=1,
                txz_left_3=0,tyz_left_3=1,
                plot_left = False,
                txz_right_1=0,tyz_right_1=0, 
                txz_right_2=1,tyz_right_2=1,
                txz_right_3=0,tyz_right_3=1,
                plot_right = False,
                nrings =6 ,tvec=3,E_z=0,r=1.4,zstart=0,b=3.6,
                chirality_left=True,
                u=-8.1,delta_e=-12.6,Es=-6.21,quantization_axis=sigma_z):
    
    
    # Create fcc lattice
    lat_fcc = kwant.lattice.general([( 0.5, 0.5, 0), 
                                 (0.5, 0, 0.5),
                                 ( 0,0.5,0.5)
                              ], norbs=18)
    
    
    L_m = int(3*nrings+1) # molecule's length
    lat = kwant.lattice.cubic(a=1,norbs =8) # lattice for helicene 
    lat_sulfur = kwant.lattice.cubic(a=1,norbs =6) # lattice for sulfur
    
    
    ### Create Gold Leads ###

    if chirality_left == True:
    #create leads
        sys_gold = GoldNormal.make_fcc_lead(nlayers, Txz,Tyz,d_L, 
                                            xi_d,delta_e,quantization_axis)
    if chirality_left == False:
        sys_gold = GoldNormal.make_fcc_lead_chiral_hop(nlayers, Txz,Tyz,d_L,
                                                       xi_d, delta_e,quantization_axis)
                                    

    ### create helicene molecule ###
    sys_hel = HeliceneSystem.make_helicene(nrings,xi_p,tvec ,E_z,r,b,zstart,chirality_left,
                                           u,u,u+7.5)
    
    
   
    
    #add gold system and helicene
    sys_gold.update(sys_hel)
    H_S_onsite = np.kron(np.diag([Es,Es,Es]) ,np.identity(2))
    
    
    ### Create Sulfur atoms ###
    # Sulfur atoms
    sys_gold[lat_sulfur(1+tvec-1,1+tvec-1,zstart+1)] = H_S_onsite  
    sys_gold[lat_sulfur(L_m + tvec+2,L_m+tvec+2,zstart+1)] = H_S_onsite  
    
    
    
    ### Au - S hopping ### 
    # Au{1,2,3} -S interface 
    
    HAU1  = Coupling_Au_S_C.HAUS(1,chirality_left)
    HAU2  = Coupling_Au_S_C.HAUS(2,chirality_left)
    HAU3  = Coupling_Au_S_C.HAUS(3,chirality_left)


    #LEFT
    sys_gold[lat_fcc(nlayers-1,txz_left_1,tyz_left_1),
             lat_sulfur(1+tvec-1,1+tvec-1,zstart+1)] = \
                HAU1
    sys_gold[lat_fcc(nlayers-1,txz_left_2,tyz_left_2),
             lat_sulfur(1+tvec-1,1+tvec-1,zstart+1)] = \
                HAU2
    sys_gold[lat_fcc(nlayers-1,txz_left_3,tyz_left_3),
             lat_sulfur(1+tvec-1,1+tvec-1,zstart+1)] = \
                HAU3
    
    #RIGHT
    sys_gold[lat_fcc( d_L,txz_right_1,tyz_right_1),
             lat_sulfur(L_m + tvec+ 2, L_m+tvec+2,zstart+1)] = \
            HAU1
    sys_gold[lat_fcc( d_L,txz_right_2,tyz_right_2),
             lat_sulfur(L_m + tvec+ 2, L_m+tvec+2,zstart+1)] = \
            HAU2
    sys_gold[lat_fcc( d_L,txz_right_3,tyz_right_3),
             lat_sulfur(L_m + tvec+2,L_m+tvec+2,zstart+1)] = \
            HAU3

    
     
    ### C-S Hopping ###
    HCS = Coupling_Au_S_C.HCS_mirrored(chirality_left)
    # S - C interface 
    sys_gold[lat(tvec+1,tvec+1,zstart),
             lat_sulfur(1+tvec-1,1+tvec-1,zstart+1)]= \
            HCS
    sys_gold[lat(L_m + tvec,L_m+tvec,zstart),
             lat_sulfur(L_m + tvec+2,L_m+tvec+2,zstart+1)]= \
            HCS
    
    
    if plot_left == True:
        GoldNormal.plot_attachement_sulfur(Txz,Tyz, 
                                              txz_left_1,tyz_left_1, 
                                              txz_left_2,tyz_left_2,
                                               txz_left_3,tyz_left_3)
    if plot_right == True:
        GoldNormal.plot_attachement_sulfur(Txz,Tyz, txz_right_1,tyz_right_1, 
                                                txz_right_2,tyz_right_2,
                                                txz_right_3,tyz_right_3)


                     
    
    return sys_gold





