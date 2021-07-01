#!/usr/bin/env python
# coding: utf-8



# notebook for hopping elements between gold & sulfur, sulfur & Carbon


# In[6]:


import numpy as np


# In[7]:


import tinyarray

# define Pauli-matrices for convenience
sigma_0 = tinyarray.array([[1, 0], [0, 1]])
sigma_x = tinyarray.array([[0, 1], [1, 0]])
sigma_y = tinyarray.array([[0, -1j], [1j, 0]])
sigma_z = tinyarray.array([[1, 0], [0, -1]])
sigma_n = tinyarray.array([[0, 0], [0, 0]])


# In[8]:


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


# In[9]:




# # Coupling Au-S-C

# Onsite energy values of the sulfur atom. Values are taken from KC his notebook.
# H_S_onsite = r2v_soc( np.diag([0.30320,0.30320,0.30320]) )
# H_S_onsite = np.kron(np.diag([-7.16,-7.16,-7.16]) ,np.identity(2)) 

# ### Coupling Au-S

# #### Hopping Hamiltonian 

# The matrices H_Au1S_low,..  below represent the hopping elements between $3p_{x,y,z}$ on Sulfur and $6s$ orbitals on gold and are given in eV. Taken from masterproject K.Choi.

# In[10]:


H_Au1S_sp = np.kron([1.308,  -2.200, 3.471],np.identity(2))
H_Au2S_sp = np.kron([1.326, 2.219, 3.503],np.identity(2))
H_Au3S_sp = np.kron( [ -2.323, 0.004, 3.191],np.identity(2))


# The matrices below  the hopping hamiltonians between $3p_{x,y,z}$ orbitals of Sulfur and $5d:d_{xy},d_{x^2 -y^2}, d_{xz},d_{yz},d_{z^2} $ orbitals of gold atoms 1,2,3.  The are are given in Hartree.

# In[11]:


HAu1_dp= hart2ev_soc([
       [ 2.04692638068954E-02,1.38238785169911E-02 ,-4.17730509119422E-02],
       [-3.37979700449427E-02,-1.58424751625304E-02,-2.97437051891800E-02],
     [-1.86711092555072E-02,-3.62593777267764E-02, 2.57542298196935E-02],
      [-3.81824183136023E-02,3.13166149136391E-02,-4.85029693810615E-02],
      [2.28663606943767E-02,-4.25648203641807E-02,-2.50357036730428E-02],
]
)
HAu2_dp= hart2ev_soc(
    [
         [-2.03988248157826E-02,1.34874211085483E-02, 4.26898159494969E-02],
         [-3.42857899216452E-02 ,1.68335324578422E-02,-2.92632810103988E-02],
         [-1.88397792262509E-02,3.74501477536945E-02,2.68433554264235E-02],
         [3.91419200421737E-02, 3.07175573688418E-02,4.97204089522636E-02],
         [ 2.40496515823447E-02, 4.41384465197868E-02 , -2.47745832963384E-02]
    ]

)


HAu3_dp= hart2ev_soc([
    [-1.61524463103285E-04,3.98008119669623E-02 ,3.78942792263684E-04],            
    [-7.02170497629566E-03,-2.45189093201272E-04,4.70096524068124E-02],
    [4.78413468268273E-02,1.61140786665182E-04,-4.28564173790477E-02],
    [2.93121277446684E-04, -3.03209230445146E-02, -1.19166810526199E-04],
    [-3.38309477724411E-02,-1.11511474517445E-04,-2.61516364433633E-02],
]

    )


# ### Mirrored hopping y $\rightarrow$ -y

# 1. $5d:[d_{xy},d_{x^2 -y^2}, d_{xz},d_{yz},d_{z^2}] \rightarrow [-d_{xy},d_{x^2 -y^2}, d_{xz},-d_{yz},d_{z^2}]$  
# 2. $3p:[p_{x},p_{y},p_{z}] \rightarrow [p_{x},-p_{y},p_{z}]$ 
# 3. $s\rightarrow s$

# In[12]:


Myp = np.kron(np.diag([1,-1,1]),sigma_0)
Myd = np.kron(np.diag([-1,1,1,-1,1]),sigma_0)


# $sp$ orbitals

# In[13]:


H_Au1S_sp_my = np.dot(H_Au1S_sp,Myp)
H_Au2S_sp_my = np.dot(H_Au2S_sp,Myp)
H_Au3S_sp_my = np.dot(H_Au3S_sp,Myp)


# $dp$ orbitals

HAu1_dp_my = np.dot(Myd,np.dot(HAu1_dp,Myp))
HAu2_dp_my = np.dot(Myd,np.dot(HAu2_dp,Myp))
HAu3_dp_my = np.dot(Myd,np.dot(HAu3_dp,Myp))

# ### Combined hopping

# In[14]:


# Under the spatial operation: y -> -y, the s,p,d orbitals changes as indicated above
# chirality left changes the hoppings of sulfur accordingly.


# In[15]:


def HAUS(site_number_gold,chirality_left):
    """ 
    Input:
    -site_number_gold = label for the gold atom to which sulfur. Take possible values 1,2,3.
    - chirality_left = boolean. If False, returns hopping that is not mirrored in y. If True, returns hopping that is mirrored  in y.
    
    Ouput:
    - Hopping matrix between Sulfur atoms and gold atom 
    """
    
    H_hopr = np.zeros((18,6)) # define empty matrix
    
    
    if chirality_left == False:
        if site_number_gold== 1:
            H_hopr[0:2,0:6] = H_Au1S_sp
            H_hopr[2:8,0:6] = np.zeros((6,6))
            H_hopr[8:18,0:6] = np.real(HAu1_dp)
            return np.round(H_hopr,3)
        if site_number_gold== 2:
            H_hopr[0:2,0:6] = H_Au2S_sp
            H_hopr[2:8,0:6] = np.zeros((6,6))
            H_hopr[8:18,0:6]= np.real(HAu2_dp)
            return np.round(H_hopr,3)
        if site_number_gold== 3:
            H_hopr[0:2,0:6] = H_Au3S_sp
            H_hopr[2:8,0:6] = np.zeros((6,6))
            H_hopr[8:18,0:6] = np.real(HAu3_dp)
            return np.round(H_hopr,3)
        
    if chirality_left == True:
        if site_number_gold== 1:
            H_hopr[0:2,0:6] = H_Au1S_sp_my
            H_hopr[2:8,0:6] = np.zeros((6,6))
            H_hopr[8:18,0:6] = np.real(HAu1_dp_my)
            return np.round(H_hopr,3)
        if site_number_gold== 2:
            H_hopr[0:2,0:6] = H_Au2S_sp_my
            H_hopr[2:8,0:6] = np.zeros((6,6))
            H_hopr[8:18,0:6]= np.real(HAu2_dp_my)
            return np.round(H_hopr,3)
        if site_number_gold== 3:
            H_hopr[0:2,0:6] = H_Au3S_sp_my
            H_hopr[2:8,0:6] = np.zeros((6,6))
            H_hopr[8:18,0:6] = np.real(HAu3_dp_my)
            return np.round(H_hopr,3)



# # C-S Coupling

# Coupling between $3p_{x,y,z}$ orbitals on sulfur and $2s,2p_{x,y,z}$ orbitals on Carbon


HCS = np.round(hart2ev_soc(
[
[-1.21268947676568E-01, -1.18826938498954E-02,2.88950698095945E-01],
[-7.71747979573655E-02, 7.70837483746968E-03,-1.04764690721968E-01],
[7.63852653738964E-03,-1.40449988697671E-01,-9.78146785092129E-03 ],
[-1.11391052083447E-01,-1.05939774089218E-02,1.21069961179607E-01]
]),3)


# In[17]:

# transformation matrix
Mysp = np.kron(np.diag([1,1,-1,1]),sigma_0)

# coupling matrix under mirror operation y -> -y:
HCS_my = np.dot(np.dot(Mysp, HCS),Myp)


# In[18]:


def HCS_mirrored(chirality_left):
    
    """ 
    Input:
    - chirality_left = boolean. If False, returns hopping that is not mirrored in y. If True, returns hopping that is mirrored  in y.
    
    Ouput:
    - Hopping matrix between sulfur atoms and carbon atom (of helicene)
    """
     
    if chirality_left == False: 
        return HCS
    
    if chirality_left == True:
        return HCS_my
        


# In[ ]:





# In[ ]:




