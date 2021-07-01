#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np


# In[2]:

#converts list to np.array
def nplist(vlist):
    varray = np.array(vlist)
    return varray

#convert units of matrix a from Rydberg to Electronvolt
def r2v(a):
    c_2eV = 13.605698066
    mat9x9= np.array(np.multiply(c_2eV,a),dtype =complex)
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
    
    
def fermi_gold_shifted(eV= False):
    
    if eV == True:
        # energy in rydberg
        energy_fermi = 0.5380*13.605698066 - 12.6
        return energy_fermi


# In[ ]:




