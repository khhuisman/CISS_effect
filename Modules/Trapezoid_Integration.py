#!/usr/bin/env python
# coding: utf-8



# In[1]:


import kwant
import math
import cmath
from matplotlib import pyplot as plt
import numpy as np
import tinyarray


from scipy import optimize
import scipy
import time



# ## Function $T_{if}(E)$

# - We plot the function $T_{if}(E)$ in an energy window $[e_{0},e_{1},..,e_{N}]$. Where $e_i = e_0 + i*de$ and create a list of values T_if_list
# - We would like this list between two continuos bounds $\mu_l,\mu_h$
# - Example: For some values of lies between two energies:$ \mu_h$: $e_m <= \mu_h < e_m + de$. So we need to know what $e_m$ is.

# # Currents

# We want to  calculate the integral: 
# - $I =\int^{\mu_h}_{\mu_l} T_{if}(E) dE = I_h + I_m + I_t  $

# Where
# 1. $ I_h = \int^{\mu_h}_{e_m} T_{if}(E) dE $
# 2. $I_m = \int^{e_m}_{\tilde{e}_m} T_{if}(E) dE $
# 3. $I_t =  \int^{\tilde{e}_m}_{\mu_l} T_{if}(E) dE $

# #### Usefull functions

# given an energy h, and a list energies.
# 1. Returns the index $n$ & energy $e[n]$ for which: $e[n] <= h < e[n] + de$
# 2. If h > max(energies) $ = e_{max}$ ,h < min(energies) $ = e_{min}$  returns $ len(energies)-1 \ \&  \ e_{max}$, $ 0 \ \& \  e_{min}$ respectively
# 

# In[308]:



def findh_upper(h,energies):
    
    for n in range(len(energies)-1):
        en = energies[n]
        en1 = energies[n+1]
        
        if h > max(energies):
            index = len(energies)-1
            eupper = max(energies)
            return index ,eupper
        
        if h < min(energies):
            index = 0
            eupper = min(energies)
            return index ,eupper
        

        if h >= en and h < en1:
            index = n
            eupper = en
            return index,eupper
    


# Given an energy l, and a list energies.
# 1. Returns the index n & energy $e[n]$ for which: $e[n]-de <= l < e[n] $
# 2. If h > max(energies) $ = e_{max}$ ,h < min(energies) $ = e_{min}$  returns $ len(energies)-1 \ \&  \ e_{max}$, $ 0 \ \& \  e_{min}$ respectively
# 
# 
# 
# 

# In[309]:


def findh_lower(l,energies):
    
    for n in range(len(energies)):
        en = energies[n]
        en1 = energies[n-1]
        
        if l > max(energies):
            index = len(energies)-1
            elower = max(energies)
            return index ,elower
        
        if l < min(energies):
            index = 0
            elower = min(energies)
            return index ,elower

        if l <= en and l > en1:
            index = n
            elower = en
            return index,elower
    
# print(findh_lower(1.5,energiesB),
# energiesB )


# ## $I_m$ : findbounds & integrate

# Input: 
# 
# 1. Upper, lower bound $\mu_h$ ,$\mu_l$ respectively
# 2. A list of energies and transmissions $T_{list}$
# 
# Ouput:  
# 1. Index ih & energy (em) for which ${e}_m  < h < {e}_m +de $
# 2. Index il & energy (emt) for which $\tilde{e}_m -de < l < \tilde{e}_m $
# 3. $I_{trapz}$: Integral $T_{list}$ from   $\tilde{e}_m$ to $e_m$ 
# 
# 
# 

# In[261]:


# T_list is an array of transmissions, which depends on the energies in the list energies.
# muh,mul are the upper,lower bound of the integral.



def findbounds_integrate(muh,mul,energies,T_list):
    mus = [muh,mul]
    h = max(mus)
    l = min(mus)

    ih,em = findh_upper(h,energies)   # find index, energy of energies for which em < h < em + de
    il,emt = findh_lower(l,energies)  # find index, energy of energies for which emt-de < l < emt
    

    energies_emt_em = []

    # find all elements of energies which satisfy : emt <= energy <= em:
    for energy in energies:
        if energy <= em and energy >= emt:
            energies_emt_em.append(energy)

    T_trapz = [ T_list[i] for i in range(il,ih+1)] # Transmissions on interval [emt , em]
    
    #  integrate  T_list  on interval [emt , em]  taking the sign of the integral in account
    I_trapz  = np.trapz(T_trapz,energies_emt_em)*np.sign(muh-mul)
    
    
    if em == max(energies):
        h = em
    if emt == min(energies):
        l = emt

    return  I_trapz,ih,h,em, il,l,emt


# # Calculate rests at boundary

# # $I_h$

# Input:
# 
# - index: corresponding to: Tiflist[index] = Tif_energy(em)
# - Tiflist: list of tranmissions
# - Tif_energy(energy,*args): function for calculating transmission from i to f at energy for *args of system. 
# - Lower,upper bound e_m,h respectively
# 
# 

# Output:
# - Integral of the transmission $T_{if}$ over $[e_m,h]$ with trapezoid method.

# In[262]:


def Icurrent_upper(index,h,em,Tif_list,Tif_energy,*args):
    fb = Tif_energy(h,*args)
    fa = Tif_list[index]
    dx = (h-em)/2
    
    Iupper = (fb + fa)*dx
    
    return Iupper
    


# # $I_l$

# Input:
# 
# - index: corresponding to: Tiflist[index] = Tif_energy(emt)
# - Tiflist: list of tranmissions
# - Function: Tif_energy(energy,*args) for calculating transmission from lead i to lead f at energy for *args of system. 
# - Lower,upper bound $l, \tilde{e}_m$ respectively

# Output:
# - Integral of the transmission $T_{if}$ over $[l, \tilde{e}_m]$ with trapezoid method.

# In[263]:


def Icurrent_lower(index,l,emt,Tif_list,Tif_energy,*args):
    
    fb = Tif_list[index]
    
    fa = Tif_energy(l,*args)
    dx = (emt-l)/2
    
    Ilower = (fb + fa)*dx
    
    return Ilower


# ## Total integral: $I =\int^{\mu_h}_{\mu_l} T_{if}(E) dE = I_h + I_m + I_t  $

# In[ ]:


def calculate_integral(muh,mul,energies,Tif_list,Tif_energy,*args):
    
    """
    Input:
    - muh = upper bound integral
    - mul = lower bound integral
    - energies = list of energies over which we numerically integrate
    - Tif_list = list of transmission corresponding to the energies in the list "energies"
    - Tif_energy= function to calculate the transmission for a certain energy that lies between two elements of the list energies
    - *args = arguments of the function Tif_energy
    
    Output:
    -( numerical) integral of transmission: "Tif" from "muh" to "mul" over energy.
    """
    
    Im_if_trapz     ,ih,h,em     , il,l,emt = findbounds_integrate(muh,mul,energies,Tif_list)
    
    
    # Calculate left overs at the upper/lower bound
    Ih_if = Icurrent_upper(ih,h,em,Tif_list,Tif_energy, 
                                               *args)   
    Il_if = Icurrent_lower(il,l,emt,Tif_list,Tif_energy, 
                                               *args)  


    
    # Total integral including letfovers. Gets a negative sign if muh < mul
    I_if_total = (Ih_if + Il_if )*np.sign(muh-mul) + Im_if_trapz
    
    return I_if_total


# In[ ]:





# In[ ]:




