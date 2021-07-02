#!/usr/bin/env python
# coding: utf-8

# Model taken from: "Chirality-Induced Spin Selectivity in a Coarse-Grained Tight-Binding Model for Helicene"
# author,year : Matthias Geyer et al. , 2019

#Input of module:
# - helix related paramters: radius r, pitch b
# - energies of 2p,2s orbitals of carbon
# - soc parameter of carbon: xip
# - slater-koster tight binding elemetns of carbon
# - chirality related boolean parameter: chirality_left. True/False relates to M/P Helicene.


# Output of Module:
# - Hopping matrices between the carbon sites of helicene
# - Onsite Hamiltonian of Helicene with spin-orbit coupling is also defined.



# In[1]:


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


# # Points on the Helix

# We have 3 helices. An Inner, Middle and Exterior helix.
# The inner helix has a constant radius r. 

# In doing x,y coordinates $\Phi_i \rightarrow - \Phi_i$ the chirality if the molecule can be flipped

# In[2]:


def R_atom_pos(Phi_i,r,b,chirality_left =True):
     # the position of the atom on the helix. 
    # With chirality_left the chirality of the molecule can be flipped

    if chirality_left == True:
        Ri_l = [r*np.cos(Phi_i ), r*np.sin(Phi_i ), b*Phi_i/(2*np.pi)]
        return Ri_l
    
    if chirality_left == False:
        Ri_r = [r*np.cos(Phi_i ), r*np.sin(-Phi_i ), b*Phi_i/(2*np.pi)]
        return Ri_r


# In[3]:


# r is the radius of the inner helix, r_M and r_E follow from it.
# b is the pitch of the helix.


    
def RI(i,r,b,chirality_left = True ):
    phi_0 = np.pi/3
    Phi_i = (i-1)*phi_0
    
    Ri = R_atom_pos(Phi_i,r,b,chirality_left)
    return Ri

def RE(i,r,b,chirality_left = True):
    phi_0 = np.pi/3
    phi_2 = 2*np.arcsin(1/(2*np.sqrt(7) ))
    phi_1 = (1/2)*(phi_0 - phi_2)
    Phi_i = np.floor(i/3)*phi_2 + phi_1*(i - 1 - np.floor(i/3))
    
    Ri =  R_atom_pos(Phi_i,r,b,chirality_left)

    return Ri

def RM(i,r,b,chirality_left = True):
    phi_0 = np.pi/3
    phi_2 = 2*np.arcsin(1/(2*np.sqrt(7) ))
    phi_1 = (1/2)*(phi_0 - phi_2)
    
    Phi_i = np.floor(i/3)*phi_2 + phi_1*(i - 1 - np.floor(i/3))

    Ri =  R_atom_pos(Phi_i,r,b,chirality_left)
    return Ri

def REM(i,r,b,chirality_left = True):
    r_m = 2*r
    r_e = np.sqrt(7)*r
        
    if i % 3 == 1:
        rv = RM(i,r_m,b,chirality_left)
        return rv
    else:
        rv = RE(i,r_e,b,chirality_left)
        return rv
        
        


# # Onsite Hamiltonian

# In[4]:
#import path to modules here.
# import sys
# sys.path.insert(0,'path_to_modules')

import onsite_soc

# energies in eV



def H_onsite_Carbon(xi_p, E_z, es=-18 , ep_sigma=-18, ep_pi=-10.5):
    
    E_x = 0
    E_y = 0
    sigma_0 = np.identity(2)
    
    # Onsite energies
    Hos = np.array(
            np.kron( np.diag([es,ep_sigma,ep_sigma,ep_pi]),sigma_0),dtype = complex)
    
    # Atomic SOC
    H_SOC = np.array(np.zeros((8,8)),dtype =complex)
    HSOC_pp = xi_p*onsite_soc.LdotS(l=1)
    H_SOC[2:8,2:8] = HSOC_pp
    
    # Combining atomic SOC and onsite energies
    Honsite = Hos + H_SOC
    
    # Electric field
    HsE = np.kron([E_x,E_y,E_z],sigma_0)
    Honsite[0:2,2:8] = HsE
    Honsite[2:8,0:2] = np.transpose(HsE)
    
    return Honsite



# print(H_od(0.6,1)-np.conjugate(np.transpose(H_od(0.6,1))))


# # Helix Vectors: Ri, Rij, ni

# ### vector Rij

# In[5]:


# h_l is the helix label
def Rij(i,j,h_li='I',h_lj='I',r=1.4,b=3.6,chirality_left=True):
    
    if h_li == 'I' and h_lj == 'I':
        dRij = np.subtract(RI(i,r,b,chirality_left), RI(j,r,b,chirality_left))
        return dRij
    if h_li == 'E' and h_lj == 'E':
        dRij = np.subtract(REM(i,r,b,chirality_left), REM(j,r,b,chirality_left))
        return dRij
    if h_li == 'I' and h_lj == 'E':
         
        if j != 3*i -2:
            print("check site number: Rij")
        
        dRij = np.subtract(RI(i,r,b,chirality_left), REM(j,r,b,chirality_left))
        return dRij
    if h_li == 'E' and h_lj == 'I':
        
        if i != 3*j -2:
            print("check site number second: Rij")
            
        dRij = np.subtract(REM(i,r,b,chirality_left), RI(j,r,b,chirality_left))
        return dRij


# ### Length of Rij

# In[6]:


#length of Rij
def length_Rij(i,j,h_li ='I',h_lj ='I',r=1.4,b=3.6,chirality_left =True):
    
    dRij = Rij(i,j,h_li,h_lj,r,b,chirality_left)
    lenght_Rij = np.linalg.norm(dRij)
    return lenght_Rij


# ### norm vec

# In[7]:


def norm_vec(mu,i,h_li = 'I',chirality_left = True):
    phi_0 = np.pi/3
    phi_2 = 2*np.arcsin(1/(2*np.sqrt(7) ))
    phi_1 = (1/2)*(phi_0 - phi_2)


    if chirality_left == True:

        if h_li == 'I':
            Phi_i = (i-1)*phi_0

        if h_li == 'E':
             Phi_i = np.floor(i/3)*phi_2 + phi_1*(i - 1 - np.floor(i/3))


        if mu ==1:
            nx = [np.cos(Phi_i ), np.sin(Phi_i ),0]
            return nx
        if mu ==2:
            ny = [- np.sin(Phi_i),np.cos(Phi_i),0]
            return ny
        if mu ==3:
            nz = [0,0,1]
            return nz
        else:
            print('Check first input nvec')
    
    if chirality_left == False:
        if h_li == 'I':
            Phi_i = (i-1)*phi_0
    
        if h_li == 'E':
             Phi_i = np.floor(i/3)*phi_2 + phi_1*(i - 1 - np.floor(i/3))


        if mu ==1:
            nx = [np.cos(-Phi_i ), np.sin(-Phi_i ),0]
            return nx
        if mu ==2:
            ny = [- np.sin(-Phi_i),np.cos(-Phi_i),0]
            return ny
        if mu ==3:
            nz = [0,0,1]
            return nz
        else:
            print('Check first input nvec')




# # n_parallel, n_orth

# Vectors that indirectly depend on the helix label. 

# ### n_parallel

# In[8]:


def nvec_parallel(mu,i,j,h_li='I',h_lj='I',r=1.4,b=3.6,chirality_left = True):
    
    Rji = Rij(j,i,h_lj,h_li,r,b,chirality_left)
    length = length_Rij(j,i,h_lj,h_li ,r,b,chirality_left)

    if mu == 0:
        return np.multiply(Rji,1/length)
    
    if mu == 1 or mu ==2 or mu == 3:
        
        nv = norm_vec(mu,i,h_li,chirality_left)
        coef = np.inner(Rji,nv)/(length**2)
        n_vij = np.multiply(coef,Rji)

        return n_vij 
    else:
        print("check mu,nu: nvec_parallel") 
        


# ### nvec_orthogonal

# In[9]:


def nvec_orthogonal(mu,i,j,h_li='I',h_lj='I',r=1.4,b=3.6,chirality_left=True):
    
    if mu == 0:
        return [0,0,0]
    
    if mu == 1 or mu ==2 or mu == 3:
        nv = norm_vec(mu,i,h_li,chirality_left)
        n_v_para = nvec_parallel(mu,i,j,h_li,h_lj,r,b,chirality_left)
        
        n_v_orth = np.subtract(nv,n_v_para)
        return n_v_orth 
    else:
        print("check mu,nu: nvec_orthogonal")


# # Hopping Hamiltonian

# ### Slater-Koster parameters

# In[10]:




def V_sigma(mu,nu,Vss=-7.92, Vsp_sigma=8.08,Vpp_sigma=7.09):
    if mu == 0 and nu ==0:
        return Vss
    if mu == 0 and ( nu == 1 or nu == 2 or nu == 3):
        return Vsp_sigma
    if nu == 0 and ( mu == 1 or mu == 2 or mu == 3):
        return Vsp_sigma
    if (nu == 1 or nu == 2 or nu == 3) and ( mu == 1 or mu == 2 or mu == 3):
        return Vpp_sigma
    else:
        print("check mu,nu: V_sigma")

def V_pi(mu,nu,Vpp_pi=-3.44):
    if mu == 0 and nu ==0:
        return 0
    if mu == 0 and ( nu == 1 or nu == 2 or nu == 3):
        return 0
    if nu == 0 and ( mu == 1 or mu == 2 or mu == 3):
        return 0
    if (nu == 1 or nu == 2 or nu == 3) and ( mu == 1 or mu == 2 or mu == 3):
        return Vpp_pi
    else:
        print("check mu,nu: V_pi")


# ### Hopping Hamiltonians

# In[11]:


def H_parallel(mu,nu,i,j,h_li='I',h_lj='I',r=1.4,b=3.6,chirality_left=True):
    
    nv_mu_p = nvec_parallel(mu,i,j,h_li,h_lj,r,b,chirality_left)
    nv_nu_p = nvec_parallel(nu,j,i,h_lj,h_li,r,b,chirality_left)

    H_paralell = np.multiply(np.inner( nv_mu_p,nv_nu_p),V_sigma(mu,nu))

    return H_paralell


# In[12]:


def H_orthogonal(mu,nu,i,j,h_li='I',h_lj='I',r=1.4,b=3.6,chirality_left = True):
    
    nv_mu_p = nvec_orthogonal(mu,i,j,h_li,h_lj,r,b,chirality_left)
    nv_nu_p = nvec_orthogonal(nu,j,i,h_lj,h_li,r,b,chirality_left)

    H_orth = np.multiply(np.inner( nv_mu_p,nv_nu_p),V_pi(mu,nu))

    return H_orth


# In[13]:


def H_hop(mu,nu,i,j,h_li='I',h_lj='I',r=1.4,b=3.6,chirality_left =True):
    H_par = H_parallel(mu,nu,i,j,h_li,h_lj,r,b,chirality_left)
    H_orth = H_orthogonal(mu,nu,i,j,h_li,h_lj,r,b,chirality_left)
    H_hop = H_par + H_orth
    return H_hop


# ### Hopping Hamiltonians

#  Hopping between sites on the inner helix

# In[14]:


def H_ij_II(i,j,r=1.4,b=3.6,chirality_left= True):
    Hij = np.array([[np.round(H_hop(mu,nu,i,j,'I','I',r,b,chirality_left),6)
                for mu in range(0,4) ]
                for nu in range(0,4) ] )
    return Hij


# Hopping between sites on the Exterior helix

# In[15]:


def H_ij_EE(i,j,r=1.4,b=3.6,chirality_left=True):
    Hij = np.array([[np.round(H_hop(mu,nu,i,j,'E','E',r,b,chirality_left),6)
                for mu in range(0,4) ]
                for nu in range(0,4) ] )
    return Hij


# #### Hopping between middle and Inner Helix

# In[16]:


def H_ij_IM(i,j,r=1.4,b=3.6,chirality_left=True):
    Hij = np.array([[np.round(H_hop(mu,nu,i,j,'I','E',r,b,chirality_left),6)
                for mu in range(0,4) ]
                for nu in range(0,4) ] )
    return Hij

def H_ij_MI(i,j,r=1.4,b=3.6,chirality_left=True):
    Hij = np.array([[np.round(H_hop(mu,nu,i,j,'E','I',r,b,chirality_left),6)
                for mu in range(0,4) ]
                for nu in range(0,4) ] )
    return Hij




# # Kwant Helicene: constant hoppings system

# ### Colour function

# In[18]:


#Colour sites according to helix label: Interior, Middle, Exterior
def family_color(site):
        n1,n2,n3 = site[1]

        #Inner helix
        if n2 == n1:
            return 'black'

        #Middle helix
        if n2 == 2+n1:
            return 'blue'

        #Exterior Helix
        if n2 == 3+n1:
            return 'green'
        else:
            return 'yellow'


# ### Hopping Definitions

# In[19]:


import numpy as np
import sympy

def H_hop_list(r,b,chirality_left= True):
    sigma_0 = np.identity(2)
    H_IIs = np.kron(H_ij_II(1,2,r,b,chirality_left),sigma_0)
    H_MEs = np.kron(H_ij_EE(1,2,r,b,chirality_left),sigma_0)
    H_EMs = np.kron(H_ij_EE(3,4,r,b,chirality_left),sigma_0)
    H_EEs = np.kron(H_ij_EE(2,3,r,b,chirality_left),sigma_0)
    H_IMs = np.kron(H_ij_IM(2,4,r,b,chirality_left),sigma_0)
    
    return [H_IIs,H_MEs,H_EMs,H_EEs,H_IMs]








