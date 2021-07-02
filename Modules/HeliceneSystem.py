#!/usr/bin/env python
# coding: utf-8

# # Helicene - Cuniberti et al.

# In[1]:

import numpy as np
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


#import path to modules here.
# import sys
# sys.path.insert(0,'path_to_modules')

# Helicene Hamiltonian is imported
import Helicene_Hamiltonian

def make_helicene(n,xi_p = 0.006,tvec=0,E_z=0,
                  r=1.4,b=3.6,zstart=0,chirality_left =True,
                 es=-18 , ep_sigma=-18, ep_pi=-10.5):
    
    
    """
    Input:
    - n = number of benzene rings
    - xi_p = spin orbit coupling paramter of carbon
    - E_z = electric field in the z direction
    - r = radius of helix
    - b = pitch of helix
    - zstart = z coordinate in (x,y,z) space where first site of helicene will be placed
    - chirality_left = boolean that controls handedness of the helix
    - es,ep_sigma,ep_pi = onsite energies of 2s, 2p sigma and 2p pi orbitals
    
    Output
    - kwant system of the helicene molecule. 
    - Note: The sytem looks flat when plotted but the geometry is hidden in the hopping matrices.
    """
    
    L= 3*n+1
    
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
            
            
            
   
            
    return hel


# In[21]:

# kwant.plot(make_helicene(n=6,xi_p = 0.006,E_z=0,r=1.4,b=3.6,chirality_left = True));

# print("make_helicene(n,xi_p,tvec,E_z,r,b,zstart,chirality_left, es , ep_sigma, ep_pi)")
                  
                

def make_helicene_empty(n,tvec=0,
                  zstart=0):
    
    """ 
    Output
    - kwant system with empty sites everywhere and no hopping except for th edge sites to which we want to attach WBL leads
    """
    
    L= 3*n+1
    
    hel = kwant.Builder()
    lat = kwant.lattice.cubic(a=1,norbs =8)
    
    sigma_08 = np.identity(8)
    sigma_z8 = np.kron(np.identity(4),sigma_z) 
    
    
    # Onsite Hamiltonian
    H_onsite = sigma_08*0
    
    
    
    
    
#     #Create sites of:
        # Inner helix
    hel[(lat(i,i,zstart) for i in range(1+tvec,L+1+tvec,3))]   = H_onsite
        #Middle helix
    hel[(lat(i,2+i,zstart) for i in range(1+tvec,L+1+tvec,3))] = H_onsite
        # Exterior helix
    hel[(lat(i+j,3+i+j,zstart)   for i in range(2+tvec,L+1+tvec,3) 
                                 for j in range(0,2))] = H_onsite


            
    return hel

# Potential linear voltage drop/ step voltage drop

def potential_linear_drop(L,V,i,tvec):
    
    
        if L > 1 and V !=0 and i >=tvec:
            pot = -V/2 + (i-tvec-1)/(L-1)*V
            return pot
        if L > 1 and V ==0:
            pot = 0
            return pot
        
        
# Potential step voltage drop

       
def potential_step(L,V,i,tvec):


    if i == L+tvec:
        pot = V/2
        return pot
    if i == 1+tvec:
        pot = -V/2
        return pot
    if i+tvec!=1 and i!= L+tvec:
        pot = 0
        return pot
    
    
    
    
def make_helicene_linear_voltage(n,V=0,xi_p = 0.006,tvec=0,E_z=0,
                  r=1.4,b=3.6,zstart=0,chirality_left =True,
                 es=-18 , ep_sigma=-18, ep_pi=-10.5):
    
    """
    Input:
    - n = number of benzene rings
    - xi_p = spin orbit coupling paramter of carbon
    - E_z = electric field in the z direction
    - r = radius of helix
    - b = pitch of helix
    - zstart = z coordinate in (x,y,z) space where first site of helicene will be placed
    - chirality_left = boolean that controls handedness of the helix
    - es,ep_sigma,ep_pi = onsite energies of 2s, 2p sigma and 2p pi orbitals
    - V = bias voltage 
    
    
    Output:
    - kwant system of helicene molecule with linear bias voltage drop along the helix axis.
    """
    
    L= 3*n+1
    
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
    
    
    for i in range(1+tvec,L+1+tvec,3):
        hel[lat(i,i,zstart) ]   = H_onsite + potential_linear_drop(L,V,i,tvec)*sigma_08
        #Middle helix
    for i in range(1+tvec,L+1+tvec,3):
        hel[lat(i,2+i,zstart) ] = H_onsite + potential_linear_drop(L,V,i,tvec)*sigma_08
        # Exterior helix
    for i in range(2+tvec,L+1+tvec,3):
        for j in range(0,2):
            hel[lat(i+j,3+i+j,zstart)   ] = H_onsite \
                                            + potential_linear_drop(L,V,i,tvec)*sigma_08





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
            
            
            
   
            
    return hel   


def make_helicene_step_voltage(n,V=0,xi_p = 0.006,tvec=0,E_z=0,
                  r=1.4,b=3.6,zstart=0,chirality_left =True,
                 es=-18 , ep_sigma=-18, ep_pi=-10.5):
    
    
    """
    Input:
    - n = number of benzene rings
    - xi_p = spin orbit coupling paramter of carbon
    - E_z = electric field in the z direction
    - r = radius of helix
    - b = pitch of helix
    - zstart = z coordinate in (x,y,z) space where first site of helicene will be placed
    - chirality_left = boolean that controls handedness of the helix
    - es,ep_sigma,ep_pi = onsite energies of 2s, 2p sigma and 2p pi orbitals
    - V = bias voltage 
    
    
    Output:
    - kwant system of helicene molecule with step bias voltage drop at the edge sites.
    """
    
    L= 3*n+1
    
    
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
    
    
    for i in range(1+tvec,L+1+tvec,3):
        hel[lat(i,i,zstart) ]   = H_onsite + potential_step(L,V,i,tvec)*sigma_08
        #Middle helix
    for i in range(1+tvec,L+1+tvec,3):
        hel[lat(i,2+i,zstart) ] = H_onsite 
        # Exterior helix
    for i in range(2+tvec,L+1+tvec,3):
        for j in range(0,2):
            hel[lat(i+j,3+i+j,zstart)   ] = H_onsite 





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
            
            
            
   
            
    return hel   
    
    
     
def make_helicene_2d(nrings,xi_p = 0.006,E_z=0,r=1.4,b=3.6
                  ,chirality_left= True,
                    es=-18 , ep_sigma=-18, ep_pi=-10.5):
    
    """Input:
    - n = number of benzene rings
    - xi_p = spin orbit coupling paramter of carbon
    - E_z = electric field in the z direction
    - r = radius of helix
    - b = pitch of helix
    - zstart = z coordinate in (x,y,z) space where first site of helicene will be placed
    - chirality_left = boolean that controls handedness of the helix
    - es,ep_sigma,ep_pi = onsite energies of 2s, 2p sigma and 2p pi orbitals
    
    Output
    - kwant system of the helicene molecule on a 2d lattices (instead of 3d). 
    """
    
    
    #n is number of benzene rings in helicene molecule
    
    L = 3*nrings+1
    hel = kwant.Builder()
    lat = kwant.lattice.square(a=1,norbs =8)
    H_os = Helicene_Hamiltonian.H_onsite_Carbon(xi_p, E_z, es, ep_sigma, ep_pi)
    
    h_list = Helicene_Hamiltonian.H_hop_list(r,b,chirality_left)
    H_IIs = h_list[0]
    H_MEs = h_list[1]
    H_EMs = h_list[2]
    H_EEs = h_list[3]
    H_IMs = h_list[4]

    #     #Create system
    # Inner helix
    hel[(lat(i,0) for i in range(1,L+1,3))] = H_os
    #Middle helix
    hel[(lat(i,2) for i in range(1,L+1,3))] = H_os
#    # Exterior helix
    hel[(lat(i+j,3) for i in range(2,L+1,3) for j in range(0,2))] = H_os



    #Hoppings Inner Helix
    for i in range(1,L+1,3):
        if i>= 4:
            hel[lat(i-3,0),lat(i,0)]= \
                     H_IIs
    
    #Hopping Exterior Helix
    for i in range(1,L,3):
        hel[lat(i,2),lat(i+1,3)] = H_MEs
    
    for i in range(3,L,3):
        hel[lat(i,3),lat(i+1,2)] = H_EMs
    
    for i in range(2,L,3):
        hel[lat(i,3),lat(i+1,3)] = H_EEs
        
        
    #Hopping Middle-Interior Helix
    for i in range(1,L+1,3):
            hel[lat(i,0),lat(i,2)]= \
                     H_IMs
            
            
    
  
    return hel    
    
    
def make_helicene_buttiker_largeprobe(n,xi_p = 0.006,tvec=0,E_z=0,
                  r=1.4,b=3.6,zstart=0,chirality_left =True,
                 es=-18 , ep_sigma=-18, ep_pi=-10.5,u0=0,t0=1):
    
    """Input:
    - n = number of benzene rings
    - xi_p = spin orbit coupling paramter of carbon
    - E_z = electric field in the z direction
    - r = radius of helix
    - b = pitch of helix
    - zstart = z coordinate in (x,y,z) space where first site of helicene will be placed
    - chirality_left = boolean that controls handedness of the helix
    - es,ep_sigma,ep_pi = onsite energies of 2s, 2p sigma and 2p pi orbitals
    
    -u0 = onsite energy of the buttiker probe
    - t0 = the hopping paramter of the buttiker probe.
    
    Output
    - kwant system of the helicene molecule on a 3d lattices where every site is connected to a Buttiker probe. 
    """
    
    L= 3*n+1
    
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
    
    
    sys_buttiker = kwant.Builder(sym_buttiker,conservation_law = sigma_z8)
  
    # Inner helix
    sys_buttiker[(lat(i,i,zstart + 1) for i in range(7+tvec,L+1+tvec-6,3))]   = H_buttiker
    #Middle helix
    sys_buttiker[(lat(i,2+i,zstart + 1) for i in range(7+tvec,L+1+tvec-6,3))] = H_buttiker
    # Exterior helix
    sys_buttiker[(lat(i+j,3+i+j,zstart + 1)   for i in range(8+tvec,L+1+tvec-6,3) 
                  for j in range(0,2))]         = H_buttiker
    
    sys_buttiker[kwant.builder.HoppingKind((0,0,1), lat, lat)]  =  -t0*np.identity(8)
    hel.attach_lead(sys_buttiker)
   
            
    return hel

def H_onsite(e0,B):
    H2x2 = np.diag([e0 +B/2,e0-B/2])
    return H2x2


sigma_z8 =np.kron(np.identity(4),sigma_z)


def make_system_helicene(e0=0,B=0,t=1,
                nrings=6,xi_p = 0.3,E_z=0,
                      r=1.4,b=3.6,chirality_left =True,
                es=-18,ep_sigma=-18,ep_pi=-10.5):
    
    """
    Input:
    
    Helicene parameters:
    - n = number of benzene rings
    - xi_p = spin orbit coupling paramter of carbon
    - E_z = electric field in the z direction
    - r = radius of helix
    - b = pitch of helix
    - zstart = z coordinate in (x,y,z) space where first site of helicene will be placed
    - chirality_left = boolean that controls handedness of the helix
    - es,ep_sigma,ep_pi = onsite energies of 2s, 2p sigma and 2p pi orbitals
    
    
    Lead parameters:
    - e0 = onsite energy of the leads
    - B = magnetization of the left lead
    - t = hopping paramter of the lead
    
    Output
    - kwant system of the helicene molecule with 2 leads attached. One lead is magnetized. 
    """

    hels = make_helicene_2d(nrings,xi_p,E_z,r,b
                  ,chirality_left,es , ep_sigma, ep_pi)
                
    #ATTACH LEADS
    
    lat =  kwant.lattice.square(a=1,norbs = 8)
    
    H_left = np.kron(np.identity(4),H_onsite(e0,B))
    #LEFT lead (=lead 0)
    sym_left = kwant.TranslationalSymmetry((-1, 0))
    lead_left = kwant.Builder(sym_left,conservation_law = sigma_z8)
    
    #onsight energy
    lead_left[lat(-1, 0) ] = H_left
    lead_left[lat(-1, 2) ] = H_left
    
    # hopping terms
    lead_left[kwant.builder.HoppingKind((1, 0), lat, lat)] =  -t*np.identity(8)
    hels.attach_lead(lead_left)
    
    



    #RIGHT lead (=lead 1)
    sym_right = kwant.TranslationalSymmetry((1, 0))
    lead_right = kwant.Builder(sym_right,conservation_law = sigma_z8)
    
    #onsite hamiltonian
    H_right = np.kron(np.identity(4),H_onsite(e0,0))
    
    lead_right[lat(nrings+1, 0)] = H_right
    lead_right[lat(nrings+1, 2)] = H_right
    #hopping
    lead_right[kwant.builder.HoppingKind((1, 0), lat, lat)] = -t*np.identity(8)
    hels.attach_lead(lead_right)

    
    return hels




