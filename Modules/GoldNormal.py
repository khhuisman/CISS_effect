#!/usr/bin/env python
# coding: utf-8

# # Import Packages

# In[1]:


import kwant
import math
import cmath
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import numpy as np


# For matrix support
import tinyarray


#import path to modules here.
# import sys
# sys.path.insert(0,'path_to_modules')

# define Pauli-matrices for convenience
sigma_0 = tinyarray.array([[1, 0], [0, 1]])
sigma_x = tinyarray.array([[0, 1], [1, 0]])
sigma_y = tinyarray.array([[0, -1j], [1j, 0]])
sigma_z = tinyarray.array([[1, 0], [0, -1]])
sigma_n = tinyarray.array([[0, 0], [0, 0]])


# ### Handy Functions

# In[2]:


def nplist(vlist):
    varray = np.array(vlist)
    return varray

def r2v(a):
    c_2eV = 13.605698066
    mat9x9= np.array(np.multiply(c_2eV,a),dtype =complex)
#     matrix =np.kron(mat9x9,sigma_0)
    return mat9x9

def r2v_soc(a):
    c_2eV = 13.605698066
    mat9x9= np.array(np.multiply(c_2eV,a),dtype =complex)
    matrix18x18 = np.kron(mat9x9,sigma_0)
    return matrix18x18

def direct_sum(v1,v2): 
    ab = np.zeros( np.add(v1.shape,v2.shape),dtype = complex )
    ab[:v1.shape[0],:v1.shape[1]]= v1
    ab[v1.shape[0]:,v1.shape[1]:]= v2
    return ab


# ### Onsite Hamiltonian

# In[3]:


ess   =  0.56220
ep = 1.27897
exzxz = 0.26097
eyzyz = 0.26097
exyxy = 0.26097 
exy2xy2 = 0.25309
ez2z2 = 0.25309


energies = [ess, ep,ep,ep, exyxy,exy2xy2,exzxz,eyzyz,ez2z2]

H_at_ev_soc = r2v_soc(np.diag(energies))


# ### Onsite SOC

# Onsite SOC is generated, only d orbitals feel SOC.

import onsite_soc

# In[5]:


def H_SO(xi_d=0.3):
    # SOC: L*S term is written in orbital basis: 
    # xy,x^2-y^2,xz,yz,z2
    
    H_soc_dorbitals = onsite_soc.LdotS(l=2)
    
    HSO_d  =  np.multiply(xi_d,H_soc_dorbitals) #Total L*S term for 5d orbitals
    
    # total onsite,atomic SOC in eV.
    H_SO_sp = np.zeros((8,8), dtype = complex ) 
    H_SO_ev = direct_sum(H_SO_sp,HSO_d) 
   
    return H_SO_ev


# # Hopping Elements: Gold

# Hhoptot is a list of all hopping elements on the gold lattice. 

# In[6]:


# Hhoptot


# In[7]:


Hhoptot = [[[-0.0668, 0.06873785019914429, 0.06873785019914429, 0, -0.04089371956670119, 0, 0, 0, 0.02361],
 [-0.06873785019914429, 0.08110500000000001, 0.097555, 0, -0.03918571216017389, 0.01340674457129694, 0, 0, 0.022623881464063588],
 [-0.06873785019914429, 0.097555, 0.08110500000000001, 0, -0.03918571216017389, -0.01340674457129694, 0, 0, 0.022623881464063588],
 [0, 0, 0, -0.01645, 0, 0, 0.01340674457129694, 0.01340674457129694, 0],
 [-0.04089371956670119, 0.03918571216017389, 0.03918571216017389, 0, -0.038424999999999994, 0, 0, 0, 0.019546193363414777],
 [0, -0.01340674457129694, 0.01340674457129694, 0, 0, 0.02624, 0, 0, 0],
 [0, 0, 0, -0.01340674457129694, 0, 0, 0.010835, 0.015405, 0],
 [0, 0, 0, -0.01340674457129694, 0, 0, 0.015405, 0.010835, 0],
 [0.02361, -0.022623881464063588, -0.022623881464063588, 0, 0.019546193363414777, 0, 0, 0, -0.015855]],
 [[-0.0668, -0.06873785019914429, 0.06873785019914429, 0, 0.04089371956670119, 0, 0, 0, 0.02361],
 [0.06873785019914429, 0.08110500000000001, -0.097555, 0, -0.03918571216017389, -0.01340674457129694, 0, 0, -0.022623881464063588],
 [-0.06873785019914429, -0.097555, 0.08110500000000001, 0, 0.03918571216017389, -0.01340674457129694, 0, 0, 0.022623881464063588],
 [0, 0, 0, -0.01645, 0, 0, -0.01340674457129694, 0.01340674457129694, 0],
 [0.04089371956670119, 0.03918571216017389, -0.03918571216017389, 0, -0.038424999999999994, 0, 0, 0, -0.019546193363414777],
 [0, 0.01340674457129694, 0.01340674457129694, 0, 0, 0.02624, 0, 0, 0],
 [0, 0, 0, 0.01340674457129694, 0, 0, 0.010835, -0.015405, 0],
 [0, 0, 0, -0.01340674457129694, 0, 0, -0.015405, 0.010835, 0],
 [0.02361, 0.022623881464063588, -0.022623881464063588, 0, -0.019546193363414777, 0, 0, 0, -0.015855]],
 [[-0.0668, 0.06873785019914429, -0.06873785019914429, 0, 0.04089371956670119, 0, 0, 0, 0.02361],
 [-0.06873785019914429, 0.08110500000000001, -0.097555, 0, 0.03918571216017389, 0.01340674457129694, 0, 0, 0.022623881464063588],
 [0.06873785019914429, -0.097555, 0.08110500000000001, 0, -0.03918571216017389, 0.01340674457129694, 0, 0, -0.022623881464063588],
 [0, 0, 0, -0.01645, 0, 0, 0.01340674457129694, -0.01340674457129694, 0],
 [0.04089371956670119, -0.03918571216017389, 0.03918571216017389, 0, -0.038424999999999994, 0, 0, 0, -0.019546193363414777],
 [0, -0.01340674457129694, -0.01340674457129694, 0, 0, 0.02624, 0, 0, 0],
 [0, 0, 0, -0.01340674457129694, 0, 0, 0.010835, -0.015405, 0],
 [0, 0, 0, 0.01340674457129694, 0, 0, -0.015405, 0.010835, 0],
 [0.02361, -0.022623881464063588, 0.022623881464063588, 0, -0.019546193363414777, 0, 0, 0, -0.015855]],
 [[-0.0668, -0.06873785019914429, -0.06873785019914429, 0, -0.04089371956670119, 0, 0, 0, 0.02361],
 [0.06873785019914429, 0.08110500000000001, 0.097555, 0, 0.03918571216017389, -0.01340674457129694, 0, 0, -0.022623881464063588],
 [0.06873785019914429, 0.097555, 0.08110500000000001, 0, 0.03918571216017389, 0.01340674457129694, 0, 0, -0.022623881464063588],
 [0, 0, 0, -0.01645, 0, 0, -0.01340674457129694, -0.01340674457129694, 0],
 [-0.04089371956670119, -0.03918571216017389, -0.03918571216017389, 0, -0.038424999999999994, 0, 0, 0, 0.019546193363414777],
 [0, 0.01340674457129694, -0.01340674457129694, 0, 0, 0.02624, 0, 0, 0],
 [0, 0, 0, 0.01340674457129694, 0, 0, 0.010835, 0.015405, 0],
 [0, 0, 0, 0.01340674457129694, 0, 0, 0.015405, 0.010835, 0],
 [0.02361, 0.022623881464063588, 0.022623881464063588, 0, 0.019546193363414777, 0, 0, 0, -0.015855]],
 [[-0.0668, 0.06873785019914429, 0., 0.06873785019914429, 0., -0.020446859783350594, -0.04089371956670119, 0., -0.011805],
 [-0.06873785019914429, 0.08110500000000001, 0., 0.097555, 0., -0.012889483794438477, -0.03918571216017389, 0., -0.022922522112824056],
 [0., 0., -0.01645, 0., 0.01340674457129694, 0., 0., 0.01340674457129694, 0.],
 [-0.06873785019914429, 0.097555, 0., 0.08110500000000001, 0., -0.026296228365735416, -0.03918571216017389, 0., 0.00029864064876047],
 [0., 0., -0.01340674457129694, 0., 0.010835, 0., 0., 0.015405, 0.],
 [-0.020446859783350594, 0.012889483794438477, 0., 0.026296228365735416, 0., -0.005331249999999999, -0.016927499999999998, 0., -0.018227669686152973],
 [-0.04089371956670119, 0.03918571216017389, 0., 0.03918571216017389, 0., -0.016927499999999998, -0.038424999999999994, 0., -0.009773096681707389],
 [0., 0., -0.01340674457129694, 0., 0.015405, 0., 0., 0.010835, 0.],
 [-0.011805, 0.022922522112824056, 0., -0.00029864064876047, 0., -0.018227669686152973, -0.009773096681707389, 0., 0.01571625]],
 [[-0.0668, -0.06873785019914429, 0., 0.06873785019914429, 0., -0.020446859783350594, 0.04089371956670119, 0., -0.011805],
 [0.06873785019914429, 0.08110500000000001, 0., -0.097555, 0., 0.012889483794438477, -0.03918571216017389, 0., 0.022922522112824056],
 [0., 0., -0.01645, 0., -0.01340674457129694, 0., 0., 0.01340674457129694, 0.],
 [-0.06873785019914429, -0.097555, 0., 0.08110500000000001, 0., -0.026296228365735416, 0.03918571216017389, 0., 0.00029864064876047],
 [0., 0., 0.01340674457129694, 0., 0.010835, 0., 0., -0.015405, 0.],
 [-0.020446859783350594, -0.012889483794438477, 0., 0.026296228365735416, 0., -0.005331249999999999, 0.016927499999999998, 0., -0.018227669686152973],
 [0.04089371956670119, 0.03918571216017389, 0., -0.03918571216017389, 0., 0.016927499999999998, -0.038424999999999994, 0., 0.009773096681707389],
 [0., 0., -0.01340674457129694, 0., -0.015405, 0., 0., 0.010835, 0.],
 [-0.011805, -0.022922522112824056, 0., -0.00029864064876047, 0., -0.018227669686152973, 0.009773096681707389, 0., 0.01571625]],
 [[-0.0668, 0.06873785019914429, 0., -0.06873785019914429, 0., -0.020446859783350594, 0.04089371956670119, 0., -0.011805],
 [-0.06873785019914429, 0.08110500000000001, 0., -0.097555, 0., -0.012889483794438477, 0.03918571216017389, 0., -0.022922522112824056],
 [0., 0., -0.01645, 0., 0.01340674457129694, 0., 0., -0.01340674457129694, 0.],
 [0.06873785019914429, -0.097555, 0., 0.08110500000000001, 0., 0.026296228365735416, -0.03918571216017389, 0., -0.00029864064876047],
 [0., 0., -0.01340674457129694, 0., 0.010835, 0., 0., -0.015405, 0.],
 [-0.020446859783350594, 0.012889483794438477, 0., -0.026296228365735416, 0., -0.005331249999999999, 0.016927499999999998, 0., -0.018227669686152973],
 [0.04089371956670119, -0.03918571216017389, 0., 0.03918571216017389, 0., 0.016927499999999998, -0.038424999999999994, 0., 0.009773096681707389],
 [0., 0., 0.01340674457129694, 0., -0.015405, 0., 0., 0.010835, 0.],
 [-0.011805, 0.022922522112824056, 0., 0.00029864064876047, 0., -0.018227669686152973, 0.009773096681707389, 0., 0.01571625]],
 [[-0.0668, -0.06873785019914429, 0., -0.06873785019914429, 0., -0.020446859783350594, -0.04089371956670119, 0., -0.011805],
 [0.06873785019914429, 0.08110500000000001, 0., 0.097555, 0., 0.012889483794438477, 0.03918571216017389, 0., 0.022922522112824056],
 [0., 0., -0.01645, 0., -0.01340674457129694, 0., 0., -0.01340674457129694, 0.],
 [0.06873785019914429, 0.097555, 0., 0.08110500000000001, 0., 0.026296228365735416, 0.03918571216017389, 0., -0.00029864064876047],
 [0., 0., 0.01340674457129694, 0., 0.010835, 0., 0., 0.015405, 0.],
 [-0.020446859783350594, -0.012889483794438477, 0., -0.026296228365735416, 0., -0.005331249999999999, -0.016927499999999998, 0., -0.018227669686152973],
 [-0.04089371956670119, -0.03918571216017389, 0., -0.03918571216017389, 0., -0.016927499999999998, -0.038424999999999994, 0., -0.009773096681707389],
 [0., 0., 0.01340674457129694, 0., 0.015405, 0., 0., 0.010835, 0.],
 [-0.011805, -0.022922522112824056, 0., 0.00029864064876047, 0., -0.018227669686152973, -0.009773096681707389, 0., 0.01571625]],
 [[-0.0668, 0., 0.06873785019914429, 0.06873785019914429, 0., 0.020446859783350594, 0., -0.04089371956670119, -0.011805],
 [0., -0.01645, 0., 0., 0.01340674457129694, 0., 0.01340674457129694, 0., 0.],
 [-0.06873785019914429, 0., 0.08110500000000001, 0.097555, 0., 0.012889483794438477, 0., -0.03918571216017389, -0.022922522112824056],
 [-0.06873785019914429, 0., 0.097555, 0.08110500000000001, 0., 0.026296228365735416, 0., -0.03918571216017389, 0.00029864064876047],
 [0., -0.01340674457129694, 0., 0., 0.010835, 0., 0.015405, 0., 0.],
 [0.020446859783350594, 0., -0.012889483794438477, -0.026296228365735416, 0., -0.005331249999999999, 0., 0.016927499999999998, 0.018227669686152973],
 [0., -0.01340674457129694, 0., 0., 0.015405, 0., 0.010835, 0., 0.],
 [-0.04089371956670119, 0., 0.03918571216017389, 0.03918571216017389, 0., 0.016927499999999998, 0., -0.038424999999999994, -0.009773096681707389],
 [-0.011805, 0., 0.022922522112824056, -0.00029864064876047, 0., 0.018227669686152973, 0., -0.009773096681707389, 0.01571625]],
 [[-0.0668, 0., -0.06873785019914429, 0.06873785019914429, 0., 0.020446859783350594, 0., 0.04089371956670119, -0.011805],
 [0., -0.01645, 0., 0., -0.01340674457129694, 0., 0.01340674457129694, 0., 0.],
 [0.06873785019914429, 0., 0.08110500000000001, -0.097555, 0., -0.012889483794438477, 0., -0.03918571216017389, 0.022922522112824056],
 [-0.06873785019914429, 0., -0.097555, 0.08110500000000001, 0., 0.026296228365735416, 0., 0.03918571216017389, 0.00029864064876047],
 [0., 0.01340674457129694, 0., 0., 0.010835, 0., -0.015405, 0., 0.],
 [0.020446859783350594, 0., 0.012889483794438477, -0.026296228365735416, 0., -0.005331249999999999, 0., -0.016927499999999998, 0.018227669686152973],
 [0., -0.01340674457129694, 0., 0., -0.015405, 0., 0.010835, 0., 0.],
 [0.04089371956670119, 0., 0.03918571216017389, -0.03918571216017389, 0., -0.016927499999999998, 0., -0.038424999999999994, 0.009773096681707389],
 [-0.011805, 0., -0.022922522112824056, -0.00029864064876047, 0., 0.018227669686152973, 0., 0.009773096681707389, 0.01571625]],
 [[-0.0668, 0., 0.06873785019914429, -0.06873785019914429, 0., 0.020446859783350594, 0., 0.04089371956670119, -0.011805],
 [0., -0.01645, 0., 0., 0.01340674457129694, 0., -0.01340674457129694, 0., 0.],
 [-0.06873785019914429, 0., 0.08110500000000001, -0.097555, 0., 0.012889483794438477, 0., 0.03918571216017389, -0.022922522112824056],
 [0.06873785019914429, 0., -0.097555, 0.08110500000000001, 0., -0.026296228365735416, 0., -0.03918571216017389, -0.00029864064876047],
 [0., -0.01340674457129694, 0., 0., 0.010835, 0., -0.015405, 0., 0.],
 [0.020446859783350594, 0., -0.012889483794438477, 0.026296228365735416, 0., -0.005331249999999999, 0., -0.016927499999999998, 0.018227669686152973],
 [0., 0.01340674457129694, 0., 0., -0.015405, 0., 0.010835, 0., 0.],
 [0.04089371956670119, 0., -0.03918571216017389, 0.03918571216017389, 0., -0.016927499999999998, 0., -0.038424999999999994, 0.009773096681707389],
 [-0.011805, 0., 0.022922522112824056, 0.00029864064876047, 0., 0.018227669686152973, 0., 0.009773096681707389, 0.01571625]],
 [[-0.0668, 0., -0.06873785019914429, -0.06873785019914429, 0., 0.020446859783350594, 0., -0.04089371956670119, -0.011805],
 [0., -0.01645, 0., 0., -0.01340674457129694, 0., -0.01340674457129694, 0., 0.],
 [0.06873785019914429, 0., 0.08110500000000001, 0.097555, 0., -0.012889483794438477, 0., 0.03918571216017389, 0.022922522112824056],
 [0.06873785019914429, 0., 0.097555, 0.08110500000000001, 0., -0.026296228365735416, 0., 0.03918571216017389, -0.00029864064876047],
 [0., 0.01340674457129694, 0., 0., 0.010835, 0., 0.015405, 0., 0.],
 [0.020446859783350594, 0., 0.012889483794438477, 0.026296228365735416, 0., -0.005331249999999999, 0., 0.016927499999999998, 0.018227669686152973],
 [0., 0.01340674457129694, 0., 0., 0.015405, 0., 0.010835, 0., 0.],
 [-0.04089371956670119, 0., -0.03918571216017389, -0.03918571216017389, 0., 0.016927499999999998, 0., -0.038424999999999994, -0.009773096681707389],
 [-0.011805, 0., -0.022922522112824056, 0.00029864064876047, 0., 0.018227669686152973, 0., -0.009773096681707389, 0.01571625]],
 [[0.00277, 0., 0., 0.00261, 0., 0., 0., 0., -0.00784],
 [0., -0.01025, 0., 0., 0., 0., 0.0047, 0., 0.],
 [0., 0., -0.01025, 0., 0., 0., 0., 0.0047, 0.],
 [-0.00261, 0., 0., 0.03707, 0., 0., 0., 0., -0.00762],
 [0., 0., 0., 0., -0.00057, 0., 0., 0., 0.],
 [0., 0., 0., 0., 0., -0.00057, 0., 0., 0.],
 [0., -0.0047, 0., 0., 0., 0., 0.0024, 0., 0.],
 [0., 0., -0.0047, 0., 0., 0., 0., 0.0024, 0.],
 [-0.00784, 0., 0., 0.00762, 0., 0., 0., 0., -0.00305]],
 [[0.00277, 0., 0., -0.00261, 0., 0., 0., 0., -0.00784],
 [0., -0.01025, 0., 0., 0., 0., -0.0047, 0., 0.],
 [0., 0., -0.01025, 0., 0., 0., 0., -0.0047, 0.],
 [0.00261, 0., 0., 0.03707, 0., 0., 0., 0., 0.00762],
 [0., 0., 0., 0., -0.00057, 0., 0., 0., 0.],
 [0., 0., 0., 0., 0., -0.00057, 0., 0., 0.],
 [0., 0.0047, 0., 0., 0., 0., 0.0024, 0., 0.],
 [0., 0., 0.0047, 0., 0., 0., 0., 0.0024, 0.],
 [-0.00784, 0., 0., -0.00762, 0., 0., 0., 0., -0.00305]],
 [[0.00277, 0, 0.00261, 0, 0, 0.006789639165669999, 0, 0, 0.00392],
 [0, -0.01025, 0, 0, 0.0047, 0, 0, 0, 0],
 [-0.00261, 0, 0.03707, 0, 0, 0.006599113576837422, 0, 0, 0.00381],
 [0, 0, 0, -0.01025, 0, 0, 0, 0.0047, 0],
 [0, -0.0047, 0, 0, 0.0024, 0, 0, 0, 0],
 [0.006789639165669999, 0, -0.006599113576837422, 0, 0, -0.00243, 0, 0, -0.001073871500692704],
 [0, 0, 0, 0, 0, 0, -0.00057, 0, 0],
 [0, 0, 0, -0.0047, 0, 0, 0, 0.0024, 0],
 [0.00392, 0, -0.00381, 0, 0, -0.001073871500692704, 0, 0, -0.00119]],
 [[0.00277, 0, -0.00261, 0, 0, 0.006789639165669999, 0, 0, 0.00392],
 [0, -0.01025, 0, 0, -0.0047, 0, 0, 0, 0],
 [0.00261, 0, 0.03707, 0, 0, -0.006599113576837422, 0, 0, -0.00381],
 [0, 0, 0, -0.01025, 0, 0, 0, -0.0047, 0],
 [0, 0.0047, 0, 0, 0.0024, 0, 0, 0, 0],
 [0.006789639165669999, 0, 0.006599113576837422, 0, 0, -0.00243, 0, 0, -0.001073871500692704],
 [0, 0, 0, 0, 0, 0, -0.00057, 0, 0],
 [0, 0, 0, 0.0047, 0, 0, 0, 0.0024, 0],
 [0.00392, 0, 0.00381, 0, 0, -0.001073871500692704, 0, 0, -0.00119]],
 [[0.00277, 0.00261, 0, 0, 0, -0.006789639165669999, 0, 0, 0.00392],
 [-0.00261, 0.03707, 0, 0, 0, -0.006599113576837422, 0, 0, 0.00381],
 [0, 0, -0.01025, 0, 0.0047, 0, 0, 0, 0],
 [0, 0, 0, -0.01025, 0, 0, 0.0047, 0, 0],
 [0, 0, -0.0047, 0, 0.0024, 0, 0, 0, 0],
 [-0.006789639165669999, 0.006599113576837422, 0, 0, 0, -0.00243, 0, 0, 0.001073871500692704],
 [0, 0, 0, -0.0047, 0, 0, 0.0024, 0, 0],
 [0, 0, 0, 0, 0, 0, 0, -0.00057, 0],
 [0.00392, -0.00381, 0, 0, 0, 0.001073871500692704, 0, 0, -0.00119]],
 [[0.00277, -0.00261, 0, 0, 0, -0.006789639165669999, 0, 0, 0.00392],
 [0.00261, 0.03707, 0, 0, 0, 0.006599113576837422, 0, 0, -0.00381],
 [0, 0, -0.01025, 0, -0.0047, 0, 0, 0, 0],
 [0, 0, 0, -0.01025, 0, 0, -0.0047, 0, 0],
 [0, 0, 0.0047, 0, 0.0024, 0, 0, 0, 0],
 [-0.006789639165669999, -0.006599113576837422, 0, 0, 0, -0.00243, 0, 0, 0.001073871500692704],
 [0, 0, 0, 0.0047, 0, 0, 0.0024, 0, 0],
 [0, 0, 0, 0, 0, 0, 0, -0.00057, 0],
 [0.00392, 0.00381, 0, 0, 0, 0.001073871500692704, 0, 0, -0.00119]]]


# # 110

# In[8]:


H_110_eV_soc = r2v_soc(Hhoptot[0])
H_m110_eV_soc = r2v_soc(Hhoptot[1])
H_1m10_eV_soc = r2v_soc(Hhoptot[2])
H_m1m10_eV_soc = r2v_soc(Hhoptot[3]) 


# # 101

# In[9]:


H_101_eV_soc = r2v_soc(Hhoptot[4])
H_m101_eV_soc = r2v_soc(Hhoptot[5])
H_10m1_eV_soc= r2v_soc(Hhoptot[6])
H_m10m1_eV_soc = r2v_soc(Hhoptot[7])


# # 011

# In[10]:


H_011_eV_soc = r2v_soc(Hhoptot[8])
H_0m11_eV_soc = r2v_soc(Hhoptot[9])
H_01m1_eV_soc = r2v_soc(Hhoptot[10])
H_0m1m1_eV_soc = r2v_soc(Hhoptot[11])


# # 001

# In[11]:


H_001_eV_soc = r2v_soc(Hhoptot[12])
H_00m1_eV_soc = r2v_soc(Hhoptot[13])


# # 010

# In[12]:


H_010_eV_soc = r2v_soc(Hhoptot[14])
H_0m10_eV_soc = r2v_soc(Hhoptot[15])


# # 100

# In[13]:


H_100_eV_soc = r2v_soc(Hhoptot[16])
H_m100_eV_soc = r2v_soc(Hhoptot[17])


# # Create Gold Leads

#  1. function that makes the left,right lead. These are blocks of shape nlayers*Txz*Tyz     with onsite atomic soc from d-orbtials $\xi_d$. $\frac{\text{d_L}}{2}$ is the distance in the xy direction between the blocks.
#  2. cut = 0,1,2. It ensures that we have ABC, CAB or BCA stacking of gold
# delta_e = shift in onsite energies of 6p,6s,5d orbitals to align system with its fermi-level
# In[14]:


def make_fcc_lead(nlayers, Txz,Tyz,d_L= 12, xi_d=0.6,delta_e=0,quantization_axis=sigma_z):
    
    
    """ 
     Input:
    - nlayers = number of gold atoms that feel SOC in xy direction
    - Txz, Tyz = dimensions of the plain gold atoms that feel SOC
    - d_L = distance between the two blocks of gold atoms.
    - xi_d = SOC parameter of 5d orbitals of gold
    - delta_e = shift in onsite energies of gold to align it to the proper fermi level.
    - quantization_axis = axis along which spins are quantized. 
    Ouput:
    - finalzed kwant system which contains two blocks of gold atoms where one block has a Buttiker probe attached. 
    - Each block has a lead attached. 
    """
    
    
    lat_fcc = kwant.lattice.general([( 0.5, 0.5, 0), 
                                 (0.5, 0, 0.5),
                                 ( 0,0.5,0.5)
                              ], norbs=18)
    
    
     
    sigma_18 = np.identity(18)
    sigmaz_18 = np.kron(np.identity(9),quantization_axis) 
    
    #n_layers = number of layers with strong SOC
    # Txy = Start value for the gold layer in xy plane.
     
    # Make a block of shape nlayers*Txz*Tyz with onsite soc xi_d.
    def sys_gold(nlayers,Txz,Tyz,xi_d,Txy):

        sys = kwant.Builder()  
        
        
        H_onsite = np.add(H_at_ev_soc,H_SO(xi_d)) + delta_e*np.identity(18)
#         

       


        #FCC lattice
        

        sys[(lat_fcc(d ,j,k) for d in range(Txy, Txy+nlayers) 
                             for j in range(0,Txz) 
                             for k in range(0,Tyz))] = H_onsite


        #Lattice hopping
        #     hopping in x direction
        sys[kwant.builder.HoppingKind((1,1,-1), lat_fcc, lat_fcc)] =             H_100_eV_soc
        #     hopping in y direction
        sys[kwant.builder.HoppingKind((1,-1,1), lat_fcc, lat_fcc)] =             H_010_eV_soc
        #     hopping in z direction
        sys[kwant.builder.HoppingKind((-1,1,1), lat_fcc, lat_fcc)] =             H_001_eV_soc



    # Diagonal hopping
        #xy direction
        sys[kwant.builder.HoppingKind((1,0,0), lat_fcc, lat_fcc)] =             H_110_eV_soc
        # x(-)y direction
        sys[kwant.builder.HoppingKind((0,1,-1), lat_fcc, lat_fcc)] =             H_1m10_eV_soc

        #hopping yz
        sys[kwant.builder.HoppingKind((0,0,1), lat_fcc, lat_fcc)] =             H_011_eV_soc
        #hopping y(-)z:V
        sys[kwant.builder.HoppingKind((1,-1,0), lat_fcc, lat_fcc)] =             H_01m1_eV_soc

        #hopping xz
        sys[kwant.builder.HoppingKind((0,1,0), lat_fcc, lat_fcc)] =              H_101_eV_soc
        #hopping x(-)z:V
        sys[kwant.builder.HoppingKind((1,0,-1), lat_fcc, lat_fcc)] =             H_10m1_eV_soc


        return sys 
    
    sys = sys_gold(nlayers,Txz,Tyz,xi_d,0)
    
    sys_reversed = sys_gold(nlayers,Txz,Tyz,xi_d,d_L)

    sys.update(sys_reversed)
    
    
    
     # Lead 
    lead1_sym =  kwant.TranslationalSymmetry((-0.5,-0.5,0))
    lead1 = kwant.Builder(lead1_sym, conservation_law = sigmaz_18)
    lead1[(lat_fcc(0,j,k)  
                         for j in range(0,Txz) 
                         for k in range(0,Tyz))] = H_at_ev_soc+ delta_e*np.identity(18)
          
            
            
    #Lead1 hopping
    #     hopping in x direction
    lead1[kwant.builder.HoppingKind((1,1,-1), lat_fcc, lat_fcc)] =         H_100_eV_soc
    #     hopping in y direction
    lead1[kwant.builder.HoppingKind((1,-1,1), lat_fcc, lat_fcc)] =         H_010_eV_soc
    #     hopping in z direction
    lead1[kwant.builder.HoppingKind((-1,1,1), lat_fcc, lat_fcc)] =         H_001_eV_soc

    
    
    # Diagonal hopping
    #xy direction
    lead1[kwant.builder.HoppingKind((1,0,0), lat_fcc, lat_fcc)] =         H_110_eV_soc
    # x(-)y direction
    lead1[kwant.builder.HoppingKind((0,1,-1), lat_fcc, lat_fcc)] =         H_1m10_eV_soc
    
    #hopping yz
    lead1[kwant.builder.HoppingKind((0,0,1), lat_fcc, lat_fcc)] =         H_011_eV_soc
    #hopping y(-)z:V
    lead1[kwant.builder.HoppingKind((1,-1,0), lat_fcc, lat_fcc)] =         H_01m1_eV_soc
    
    #hopping xz
    lead1[kwant.builder.HoppingKind((0,1,0), lat_fcc, lat_fcc)] =          H_101_eV_soc
    #hopping x(-)z:V
    lead1[kwant.builder.HoppingKind((1,0,-1), lat_fcc, lat_fcc)] =         H_10m1_eV_soc
                         
    
    sys.attach_lead(lead1)
    sys.attach_lead(lead1.reversed())
    
    return sys
    
# This system make_fcc_lead where y-> -y, this corresponds to changing all hoppings as: H_{l,m,n} -> H_{l,-m,n} 
def make_fcc_lead_chiral_hop(nlayers, Txz,Tyz,d_L= 12, xi_d=0.6,delta_e=0,quantization_axis=sigma_z):
    
    
    """ 
     Input:
    - nlayers = number of gold atoms that feel SOC in xy direction
    - Txz, Tyz = dimensions of the plain gold atoms that feel SOC
    - d_L = distance between the two blocks of gold atoms.
    - xi_d = SOC parameter of 5d orbitals of gold
    - delta_e = shift in onsite energies of gold to align it to the proper fermi level.
    - quantization_axis = axis along which spins are quantized. 
    Ouput:
    - finalzed kwant system which contains two blocks of gold atoms. 
    - Each block has a lead attached. 
    - The hopping matrices of gold have been changes:H_{l,m,n} -> H_{l,-m,n}
    """
    
    
    lat_fcc = kwant.lattice.general([( 0.5, 0.5, 0), 
                                 (0.5, 0, 0.5),
                                 ( 0,0.5,0.5)
                              ], norbs=18)
    
    
     
    sigma_18 = np.identity(18)
    sigmaz_18 = np.kron(np.identity(9),quantization_axis) 
    
    #n_layers = number of layers with strong SOC
    # Txy = Start value for the gold layer in xy plane.
     
    # Make a block of shape nlayers*Txz*Tyz with onsite soc xi_d.
    def sys_gold(nlayers,Txz,Tyz,xi_d,Txy):

        sys = kwant.Builder()  
        
        
        H_onsite = np.add(H_at_ev_soc,H_SO(xi_d)) + delta_e*np.identity(18)
        

       


        #FCC lattice
        

        sys[(lat_fcc(d ,j,k) for d in range(Txy, Txy+nlayers) 
                             for j in range(0,Txz) 
                             for k in range(0,Tyz))] = H_onsite


        #Lattice hopping
        #     hopping in x direction
        sys[kwant.builder.HoppingKind((1,1,-1), lat_fcc, lat_fcc)] =             H_100_eV_soc
        #     hopping in y direction
        sys[kwant.builder.HoppingKind((1,-1,1), lat_fcc, lat_fcc)] =             H_0m10_eV_soc # instead of H_010_ev_soc
        #     hopping in z direction
        sys[kwant.builder.HoppingKind((-1,1,1), lat_fcc, lat_fcc)] =             H_001_eV_soc



    # Diagonal hopping
        #xy direction
        sys[kwant.builder.HoppingKind((1,0,0), lat_fcc, lat_fcc)] =             H_1m10_eV_soc # instead of H_110_ev_soc
        # x(-)y direction
        sys[kwant.builder.HoppingKind((0,1,-1), lat_fcc, lat_fcc)] =             H_110_eV_soc # instead of H_1m10_ev_soc

        #hopping yz
        sys[kwant.builder.HoppingKind((0,0,1), lat_fcc, lat_fcc)] =             H_0m11_eV_soc # instead of H_011_ev_soc
        #hopping y(-)z:V
        sys[kwant.builder.HoppingKind((1,-1,0), lat_fcc, lat_fcc)] =             H_0m1m1_eV_soc  # instead of H_01m1_ev_soc

        #hopping xz
        sys[kwant.builder.HoppingKind((0,1,0), lat_fcc, lat_fcc)] =              H_101_eV_soc
        #hopping x(-)z
        sys[kwant.builder.HoppingKind((1,0,-1), lat_fcc, lat_fcc)] =             H_10m1_eV_soc


        return sys 
    
    sys = sys_gold(nlayers,Txz,Tyz,xi_d,0)
    
    sys_reversed = sys_gold(nlayers,Txz,Tyz,xi_d,d_L)

    sys.update(sys_reversed)
    
    
    
     # Lead 
    lead1_sym =  kwant.TranslationalSymmetry((-0.5,-0.5,0))
    lead1 = kwant.Builder(lead1_sym, conservation_law = sigmaz_18)
    lead1[(lat_fcc(0,j,k)  
                         for j in range(0,Txz) 
                         for k in range(0,Tyz))] = H_at_ev_soc+ delta_e*np.identity(18)
          
            
            
    #Lead1 hopping
    #     hopping in x direction
    lead1[kwant.builder.HoppingKind((1,1,-1), lat_fcc, lat_fcc)] =         H_100_eV_soc
    #     hopping in y direction
    lead1[kwant.builder.HoppingKind((1,-1,1), lat_fcc, lat_fcc)] =         H_0m10_eV_soc # instead of H_010_eV_soc
    #     hopping in z direction
    lead1[kwant.builder.HoppingKind((-1,1,1), lat_fcc, lat_fcc)] =         H_001_eV_soc

    
    
    # Diagonal hopping
    #xy direction
    lead1[kwant.builder.HoppingKind((1,0,0), lat_fcc, lat_fcc)] =         H_1m10_eV_soc 
    # x(-)y direction
    lead1[kwant.builder.HoppingKind((0,1,-1), lat_fcc, lat_fcc)] =         H_110_eV_soc
    
    #hopping yz
    lead1[kwant.builder.HoppingKind((0,0,1), lat_fcc, lat_fcc)] =         H_0m11_eV_soc
    #hopping y(-)z:V
    lead1[kwant.builder.HoppingKind((1,-1,0), lat_fcc, lat_fcc)] =         H_0m1m1_eV_soc
    
    #hopping xz
    lead1[kwant.builder.HoppingKind((0,1,0), lat_fcc, lat_fcc)] =          H_101_eV_soc
    #hopping x(-)z:V
    lead1[kwant.builder.HoppingKind((1,0,-1), lat_fcc, lat_fcc)] =         H_10m1_eV_soc
                         
    
    sys.attach_lead(lead1)
    sys.attach_lead(lead1.reversed())
    
    return sys




# In[55]:


def plot_sys_leads(nlayers, Txz=2,Tyz=2,d_L= 12, xi_d=0.6):
    
    """
    Input:
     Input:
    - nlayers = number of gold atoms that feel SOC in xy direction
    - Txz, Tyz = dimensions of the plain gold atoms that feel SOC
    - d_L = distance between the two blocks of gold atoms.
    - xi_d = SOC parameter of 5d orbitals of gold
   
    
    Output:
    - plot of the scattering region
    """
    
    sys_leads = make_fcc_lead( nlayers, Txz,Tyz,d_L, xi_d)
    kwant.plot(sys_leads,num_lead_cells=1 ,site_size=0.4, site_lw=0.01, hop_lw=0.05
              )
    sys_leads = sys_leads.finalized()


# In[54]:


# plot_sys_leads(nlayers=3)


# # In[ ]:


# print('make_fcc_lead(nlayers, Txz,Tyz,d_L, xi_d,check_conjugate,delta_e)')


# In[ ]:


from matplotlib import pyplot as plt


# Function that plots the surface atoms of gold leads in gold/yellow (at Au-molecule interface).
# Red/Blue point indicates the site to which left/right lead is attached.
# A Purple point means that the attachement position left/right lead is the same.
# A green point means has the same meaning of a purple point plus the fact that molecule is attached to the most central symmetric atom of the gold surface.
def plot_attachement(Txz=3,Tyz=3,
             txz_left=1,tyz_left=1,
             txz_right=2,tyz_right=2):


    list_xz = [j for j in range(Tyz) for k in range(Txz)]
    list_yz = [k for j in range(Tyz) for k in range(Txz)]

    nxz_sym = 0.5*(Txz-1)
    nyz_sym = 0.5*(Tyz-1)

   
    plt.subplots(figsize=(10, 5))

    sub1= plt.subplot(1, 2, 1)
    sub1.set_title(' Lead Attachement position')
    plt.scatter(list_yz, list_xz ,  color = 'gold')
    
    if tyz_left == tyz_right and txz_right == txz_left:
        
        if tyz_left == nyz_sym and txz_right == txz_left:
            plt.scatter(tyz_left, txz_left ,color = 'lightgreen')
        else:
            plt.scatter(tyz_left, txz_left ,color = 'purple')
    else:
        plt.scatter(tyz_left,  txz_left ,color = 'blue')
        plt.scatter(tyz_right, txz_right ,color = 'red')



    plt.tight_layout()
    plt.show()
    
# print('plot_attachement(Txz,Tyz, txz_left,tyz_left, txz_right,tyz_right)')


# Plot of the gold surface, a point corresponds to a gold atom.  The colour of the point indicates whether an atom is attached to sulfur:
# yellow: Au atom not attached to sulfur
# blue,red,green indicate that atom is attached and that hopping for gold atom 1,2,3 is used respectively.



def plot_attachement_sulfur(Txz=3,Tyz=3,
             txz_left_1=1,tyz_left_1=1,
             txz_left_2 =2,tyz_left_2=2,
             txz_left_3=0,tyz_left_3=0):

    # Make grid of gold surface:
    vecyz = np.array([0, 0.5,0.5])
    vecxz = np.array([0.5, 0, 0.5])
    
    x_list = []
    y_list = []
    z_list = []



    for k in range(Tyz):
        for i in range(Txz):
            rvec = np.add(np.multiply(k,vecyz), np.multiply(i,vecxz))
        
            
            z_list.append(rvec[2])
            y_list.append(rvec[1])
            x_list.append(rvec[0])

       
   
    plt.subplots(figsize=(10, 5))

    sub1= plt.subplot(1, 2, 1)
    sub1.set_title(' Lead Attachement position')
#     plt.scatter(x_list, y_list ,  color = 'gold')
    plt.scatter(x_list, z_list ,  color = 'gold')
#     plt.scatter(y_list, z_list ,  color = 'gold')



    
    rvec_1 = tyz_left_1*vecyz +  txz_left_1*vecxz
    rvec_2 = tyz_left_2*vecyz +  txz_left_2*vecxz
    rvec_3 = tyz_left_3*vecyz +  txz_left_3*vecxz
    
    plt.scatter(rvec_1[0],rvec_1[2] ,color = 'blue')
    plt.scatter(rvec_2[0], rvec_2[2] ,color = 'red')
    plt.scatter(rvec_3[0], rvec_3[2] ,color = 'green')


    plt.tight_layout()
    plt.show()


def sys_gold_bulk(nlayers,Txz,Tyz,xi_d,Txy):

        sys = kwant.Builder()  
        lat_fcc = kwant.lattice.general([
                                ( 0.5, 0.5, 0), 
                                 (0.5, 0, 0.5),
                                 ( 0,0.5,0.5)
                              ], norbs=18)
    
    
        sigma_z = tinyarray.array([[1, 0], [0, -1]])
        sigma_18 = np.identity(18)
        sigmaz_18 = np.kron(np.identity(9),sigma_z)
        H_onsite = np.add(H_at_ev_soc,H_SO(xi_d))
       


        #FCC lattice
        

        sys[(lat_fcc(d ,j,k) for d in range(Txy, Txy+nlayers) 
                             for j in range(0,Txz) 
                             for k in range(0,Tyz))] = H_onsite


        #Lattice hopping
        #     hopping in x direction
        sys[kwant.builder.HoppingKind((1,1,-1), lat_fcc, lat_fcc)] =             H_100_eV_soc
        #     hopping in y direction
        sys[kwant.builder.HoppingKind((1,-1,1), lat_fcc, lat_fcc)] =             H_010_eV_soc
        #     hopping in z direction
        sys[kwant.builder.HoppingKind((-1,1,1), lat_fcc, lat_fcc)] =             H_001_eV_soc



    # Diagonal hopping
        #xy direction
        sys[kwant.builder.HoppingKind((1,0,0), lat_fcc, lat_fcc)] =             H_110_eV_soc
        # x(-)y direction
        sys[kwant.builder.HoppingKind((0,1,-1), lat_fcc, lat_fcc)] =             H_1m10_eV_soc

        #hopping yz
        sys[kwant.builder.HoppingKind((0,0,1), lat_fcc, lat_fcc)] =             H_011_eV_soc
        #hopping y(-)z:V
        sys[kwant.builder.HoppingKind((1,-1,0), lat_fcc, lat_fcc)] =             H_01m1_eV_soc

        #hopping xz
        sys[kwant.builder.HoppingKind((0,1,0), lat_fcc, lat_fcc)] =              H_101_eV_soc
        #hopping x(-)z:V
        sys[kwant.builder.HoppingKind((1,0,-1), lat_fcc, lat_fcc)] =             H_10m1_eV_soc

        # Lead 
        lead1_sym =  kwant.TranslationalSymmetry((-0.5,-0.5,0))
        lead1 = kwant.Builder(lead1_sym, conservation_law = sigmaz_18)
        lead1[(lat_fcc(0,j,k)  
                             for j in range(0,Txz) 
                             for k in range(0,Tyz))] = H_at_ev_soc



        #Lead1 hopping
        #     hopping in x direction
        lead1[kwant.builder.HoppingKind((1,1,-1), lat_fcc, lat_fcc)] =         H_100_eV_soc
        #     hopping in y direction
        lead1[kwant.builder.HoppingKind((1,-1,1), lat_fcc, lat_fcc)] =         H_010_eV_soc
        #     hopping in z direction
        lead1[kwant.builder.HoppingKind((-1,1,1), lat_fcc, lat_fcc)] =         H_001_eV_soc



        # Diagonal hopping
        #xy direction
        lead1[kwant.builder.HoppingKind((1,0,0), lat_fcc, lat_fcc)] =         H_110_eV_soc
        # x(-)y direction
        lead1[kwant.builder.HoppingKind((0,1,-1), lat_fcc, lat_fcc)] =         H_1m10_eV_soc

        #hopping yz
        lead1[kwant.builder.HoppingKind((0,0,1), lat_fcc, lat_fcc)] =         H_011_eV_soc
        #hopping y(-)z:V
        lead1[kwant.builder.HoppingKind((1,-1,0), lat_fcc, lat_fcc)] =         H_01m1_eV_soc

        #hopping xz
        lead1[kwant.builder.HoppingKind((0,1,0), lat_fcc, lat_fcc)] =          H_101_eV_soc
        #hopping x(-)z:V
        lead1[kwant.builder.HoppingKind((1,0,-1), lat_fcc, lat_fcc)] =         H_10m1_eV_soc


        sys.attach_lead(lead1)
        sys.attach_lead(lead1.reversed())
        return sys 
    

def make_fcc_lead_LR(nlayers, Txz,Tyz, gamma, d_L= 12):
    
    
    lat_fcc = kwant.lattice.general([( 0.5, 0.5, 0), 
                                 (0.5, 0, 0.5),
                                 ( 0,0.5,0.5)
                              ], norbs=18)
    
    
     
    def gold_gamma_L(nlayers,Txz,Tyz,Txy,gamma):

        gold_system = kwant.Builder()  
        lat_fcc = kwant.lattice.general([( 0.5, 0.5, 0), 
                                 (0.5, 0, 0.5),
                                 ( 0,0.5,0.5)
                              ], norbs=18)
    
    
        #FCC lattice
        
        for i in range(Txy,Txy+nlayers):
            if i == Txy:
                gold_system[( lat_fcc(i,j,k)  
                                        for j in range(Txz)
                                        for k in range(Tyz))]  = gamma*np.identity(18)
            if i != Txy:
                gold_system[( lat_fcc(i,j,k)  
                                        for j in range(Txz)
                                        for k in range(Tyz))]  = 0*np.identity(18)

        return gold_system
    
    
    def gold_gamma_R(nlayers,Txz,Tyz,Txy,gamma):

        gold_system = kwant.Builder()  
        lat_fcc = kwant.lattice.general([( 0.5, 0.5, 0), 
                                 (0.5, 0, 0.5),
                                 ( 0,0.5,0.5)
                              ], norbs=18)
    
    
        #FCC lattice
        
        for i in range(Txy,Txy+nlayers):
            if i == Txy+nlayers-1:
                gold_system[( lat_fcc(i,j,k)  
                                        for j in range(Txz)
                                        for k in range(Tyz))]  = gamma*np.identity(18)
            if i != Txy+nlayers-1:
                gold_system[( lat_fcc(i,j,k)  
                                        for j in range(Txz)
                                        for k in range(Tyz))]  = 0*np.identity(18)

        return gold_system
    
    def gold_empty(nlayers,Txz,Tyz,Txy):

        sys_empty = kwant.Builder()  
        lat_fcc = kwant.lattice.general([(0.5, 0.5,   0), 
                                         (0.5, 0  , 0.5),
                                         (0  ,0.5 , 0.5)
                              ], norbs=18)
    
    
        #FCC lattice
        
        sys_empty[( lat_fcc(i,j,k) 
                                for i in range(Txy,Txy+nlayers)
                                for j in range(Txz)
                                for k in range(Tyz))]  = 0*np.identity(18)
            
                

        return sys_empty 
    
    #Make left gamma
    sys_left = gold_gamma_L(nlayers,Txz,Tyz,0,gamma)
    sys_left_empty = gold_empty(nlayers,Txz,Tyz,d_L)
    sys_left.update(sys_left_empty)
    
    #Make right gamma
    sys_right_empty      = gold_empty(nlayers,Txz,Tyz,0)
    sys_right = gold_gamma_R(nlayers,Txz,Tyz,d_L,gamma)
    sys_right.update(sys_right_empty)

    
    
    return sys_left,sys_right








    
