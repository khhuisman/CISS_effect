#!/usr/bin/env python
# coding: utf-8

# # General system and one Buttiker lead added.

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
sigma_0 = tinyarray.array([[1, 0], [0, 1]])
sigma_x = tinyarray.array([[0, 1], [1, 0]])
sigma_y = tinyarray.array([[0, -1j], [1j, 0]])
sigma_z = tinyarray.array([[1, 0], [0, -1]])




# In[3]:


def list_tag_transmission_functions(nleads,systf,emin,emax,de):
    """ Input: 
            - finalized kwant system (systf) with nleads attached
            - energy window: [emin,emax,de] on which the transmission are calculated.  
        Ouput:  
            - 2 Lists of transmissions from lead i to j (Tij) and j to i (Tji)  where i != j:
                * One list for the + magnetization and one for the opposite - magnetization
                * IMPORTANT: Buttiker reciproical relations are used here: Tij(m) = Tji(-m)!
            - tag_list_ij,tag_list_ji are explicit labels for the calculated transmission 
                (to keep track of what transmission is calculted )"""
    

    list_transmissions_P = []
    list_transmissions_M = []

    for lead_i in range(nleads):
        energies, Tij_list_i,tag_list_ij_i  = CalculateCurrentN.Tij(systf,lead_i =lead_i,n_leads=nleads,
                                                  emin=emin,emax=emax, de=de,
                                                  plot_var=False,check_hermiticity=True)

        energies, Tji_list_i,tag_list_ji_i   = CalculateCurrentN.Tji(systf,lead_i =lead_i,n_leads=nleads,
                                                  emin=emin,emax=emax, de=de,
                                                  plot_var=False,check_hermiticity=True)



        list_transmissions_P.append([ Tij_list_i,tag_list_ij_i, Tji_list_i,tag_list_ji_i ])
        
        # Using the Buttiker reciprocity relations
        list_transmissions_M.append([ Tji_list_i,tag_list_ij_i, Tij_list_i,tag_list_ji_i ])
        
    return energies,list_transmissions_P,list_transmissions_M 



def list_tag_transmission_functions_nobuttiker(nleads,systf,emin,emax,de):
    """ Input: 
            - finalized kwant system (systf) with nleads attached
            - energy window: [emin,emax,de] on which the transmission are calculated.  
        Ouput:  
            - 2 Lists of transmissions from lead i to j (Tij) and j to i (Tji)  where i != j:
                * One list for the + magnetization and one for the opposite - magnetization
                * IMPORTANT: Buttiker reciproical relations are NOT used here: Tij(m) = Tji(-m)!
                             if the system breaks TRS, then  Tij(m) != Tji(-m)
            - tag_list_ij,tag_list_ji are explicit labels for the calculated transmission 
                (to keep track of what transmission is calculted )"""
    

    list_transmissions = []

    for lead_i in range(nleads):
        energies, Tij_list_i,tag_list_ij_i  = CalculateCurrentN.Tij(systf,lead_i =lead_i,n_leads=nleads,
                                                  emin=emin,emax=emax, de=de,
                                                  plot_var=False,check_hermiticity=True)

        energies, Tji_list_i,tag_list_ji_i   = CalculateCurrentN.Tji(systf,lead_i =lead_i,n_leads=nleads,
                                                  emin=emin,emax=emax, de=de,
                                                  plot_var=False,check_hermiticity=True)



        list_transmissions.append([ Tij_list_i,tag_list_ij_i, Tji_list_i,tag_list_ji_i ])
        
       
        
    return energies,list_transmissions 


# In[4]:


def list_tag_transmission_functions_prime(nleads,systf,emin,emax,de):
    """ Input: 
            - finalized kwant system (systf) with nleads attached
            - energy window: [emin,emax,de] on which the transmission are calculated.  
        Ouput:  
            - 2 Lists of transmissions from lead i to j (Tij) and j to i (Tji)  where i != j:
                * One list for the + magnetization and one for the opposite - magnetization
                * IMPORTANT: Buttiker reciproical relations are used here: Tij(m) = Tji(-m)!
            - tag_list_ij,tag_list_ji are explicit labels for the calculated transmission 
                (to keep track of what transmission is calculted )"""
    

    list_transmissions_P = []
    list_transmissions_M = []

    for lead_i in range(nleads):
        energies, Tij_list_i,tag_list_ij_i  = CalculateCurrentN.Tij(systf,lead_i =lead_i,n_leads=nleads,
                                                  emin=emin,emax=emax, de=de,
                                                  plot_var=False,check_hermiticity=True)

        energies, Tji_list_i,tag_list_ji_i   = CalculateCurrentN.Tji(systf,lead_i =lead_i,n_leads=nleads,
                                                  emin=emin,emax=emax, de=de,
                                                  plot_var=False,check_hermiticity=True)



        list_transmissions_P.append([  Tji_list_i,tag_list_ji_i,Tij_list_i,tag_list_ij_i ])
        
        # Using the Buttiker reciprocity relations
        list_transmissions_M.append([  Tij_list_i,tag_list_ji_i , Tji_list_i,tag_list_ij_i])
        
    return energies ,list_transmissions_P,list_transmissions_M


# In[5]:


def overlap_two_lists(list1,list2):
    "Returns the overlapping elements of list1 and list2"
    joinedlist = sorted(list(set(list1+list2)) )
    overlap_list = []

    for element in joinedlist:

        if element in list1 and element in list2:
            overlap_list.append(element)
            
            
    return overlap_list


def index_list_adjusted(IP_list_2,IM_list_2,V_list,nacc=7,plot_current=True):
    
    """ Input: - Buttiker probe currents for plus/minus magnetization 
              - list of voltages V_list,
              - accuracy for convergence 
       Ouput: - List of voltages for which both Buttiker currents converged, the default accuracy is 10**-7.
              - List of indices corresponding to returned voltage list."""
    
    if plot_current == True:
        plt.plot(V_list,IM_list_2,label = 'IB(-m)')
        plt.plot(V_list,IP_list_2,label = 'IB(+m)')
        plt.legend()
        plt.xlabel('Voltage ')
        plt.ylabel('Charge Current ')




        plt.show()
    
    index_list_p = []
    index_list_m = []
    
    
    # Values in IM_list_2,IP_list_2 that are to large (current IB did not converge ) are ignored
    # corresponding voltage are added to a list
    for i in range(len(IM_list_2)):
        
        IBm = np.round(IM_list_2[i],nacc)
        
        if IBm == 0:
            index_list_m.append(i)

        if IBm != 0:
            print(i,V_list[i],IM_list_2[i])
            
    for i in range(len(IP_list_2)):
        IBp = np.round(IP_list_2[i],nacc) 
        if IBp == 0:
            index_list_p.append(i)

        if IBp != 0:
            print(i,V_list[i],IP_list_2[i])
            
            
    # Find the overlapping set
    overlap_list = overlap_two_lists(index_list_p,index_list_m)
    
    V_list_adjusted = [ V_list[i] for i in overlap_list]
    
    return overlap_list, V_list_adjusted


# In[6]:



def calc_MR_list(IP_list_0,IM_list_0,V_list):
    
    """Input: 
    Currents into left lead for plus/minus magnetization of the lead, for the voltages in V_list.      
    Ouput: 
    MR for the voltages in V_list"""
    #At V = 0 we must have I(m) +I(-m) = 0
    # numerically we often find:  I(m) +I(-m) = 10**-12, which leads to division of a small number in the MR.
    # therefore I round I(m) +I(-m) with n:
#     if V_list[0] ==0:
        
#         T = IM_list_0[0] + IP_list_0[0]

#         if T !=0:
#             n = int(round(-np.log10(abs(T)))) -1
#         if T ==0:
#             n=16
#     elif V_list[0] !=0:
#         n =16

   
    
    n=12
    #Calculate MR"
    MR_list = []

    for k in range(len(V_list)):
        V = V_list[k]
        ILm = IM_list_0[k]
        ILp = IP_list_0[k]

        T = ILp + ILm
        
        
        if np.round(T,n) == 0:
            MR_list.append(0)
        if np.round(T,n) != 0:

            P = 100*(ILp - ILm)/T


            MR_list.append(P)
            
    return MR_list


# In[7]:


def func_list_weighted_transmissions(nleads,ef,V, muB_list,
                                list_tags_transmissions,
                                energies,beta_list=[10**11,10**11,10**11]):
    
    """
    Input:
        
        - beta_list = list of temperatured for lead 0,1,2,... nleads
        
        
    Ouput:
        - list of transmission weighted with temperature.
        
    """
    
    list_weighted_transmissions = []
    
    for lead_i in range(nleads):
        
        beta = beta_list[lead_i]
        
        
        weighted_transmissions_listi = CalculateCurrentN.weighted_transmissions(lead_i,nleads,
                           beta=beta,ef=ef,V=V,muB_list=muB_list,
                          list_tags_transmissions=list_tags_transmissions,energies=energies)
        
        list_weighted_transmissions.append(weighted_transmissions_listi)
        
    return list_weighted_transmissions


# #### Plot functions 

# In[8]:


import PlotFunctions


# In[9]:



# ### Test system


# In[14]:


import CalculateCurrentN





def muB_selfconsistent(ef,V_list,nleads,list_tags_transmissions,energies,system,print_bool,tol = 10**-17):

    ''' Input:
    - ef = fermi level/fermi energy
    - V_list = list of voltages, first value (V_list[0]) should be close to ZERO. 
            Important: Always start with a small first voltage to improve convergence speed 
    - nleads = number of leads (terminals) attached to the scattering region.
    - list_tags_transmissions = list of transmissions and tags
    - energies = energy window to integrate over
    - system = the finalized kwant system 
    - print_bool = weither or not to print the ouput of each self-consistency step.
    - tol = convergence criterium.
    Ouput:
    - muBs_sol_list = List of chemical potentials for buttiker probes 1,2,3... s.t. the current into each probe IB1,IB2,... is zero for voltage V at T =0.'''
    

    muBs_sol_list = []

    
    #intial guess for muB's 
    # make i
    

    for n in range(len(V_list)):
        V = V_list[n]
        if n == 0:
            muBs0_P_list = [ef for i in range(nleads-2)]
        if n != 0:
            muBs0_P_list = muBs_sol_list[n-1]


        print('--- V = {} --- '.format(V))
        sol = optimize.root(CalculateCurrentN.currents_B_probe, muBs0_P_list, method='hybr',args = (ef,V,nleads,
                                                                                list_tags_transmissions,energies,
                                                                                  system,print_bool) ,
                                                                        tol = tol)

        x_list_sol = sol.x
        muBs_sol_list.append(x_list_sol)

        print('--- mu2,mu3,... ={}'.format(x_list_sol))
        
        
    return muBs_sol_list



def currents_list(muBs_sol_list,V_list,ef,nleads,list_tags_transmissions,energies,
                                                                                  system,print_voltage=False):
    
    
    ''' Input:
    muBs_sol_list = List of chemical potentials for buttiker probes 1,2,3... s.t.
                    the current into each probe IB1,IB2,... is zero for voltage V
    V_list        = list of voltages
    ef            = fermi level
    nleads        = number of leads (terminals) attached to the scattering region.
    list_tags_transmissions = list of transmissions and tags
    energies      = energy window to integrate over
    system        = the finalized kwant system 
    print_voltage    = weither or not to print eacht voltage
    Ouput:
    - Current into left,right lead and buttiker probe(s)
    - plot to check that IL +IR = 0 and IBs = 0
    '''
    
    
    
    I_list_0 = []
    I_list_1 = []
    I_list_2 = []


    for n in range(len(V_list)):

        # Generate muLs
        muB_list = muBs_sol_list[n]
        V        = V_list[n]

        #Current into lead 0
        I0 = CalculateCurrentN.currents_leadi_probe(0,muB_list,ef,V,nleads,
                             energies,list_tags_transmissions,system,
                             print_bool=False)

        I_list_0.append(I0)

        #Current into lead 1
        I1 = CalculateCurrentN.currents_leadi_probe(1,muB_list,ef,V,nleads,
                             energies,list_tags_transmissions,system,
                             print_bool=False)

        I_list_1.append(I1)



        #Current into lead 2,3,..,nleads-2
        IB_list = CalculateCurrentN.currents_B_probe( muB_list, ef,V,nleads,list_tags_transmissions,energies,
                                                                                  system,False)

        I_list_2.append(IB_list)
        if print_voltage == True:
            print('V= {}'.format(V))
            
            
            
    plt.plot(V_list,I_list_1,label = 'IR')
    plt.plot(V_list,I_list_0,label = 'IL')
    plt.legend()



    plt.show()

    plt.plot(V_list,I_list_2,label = 'IB')
    plt.plot(V_list,np.add(I_list_1,I_list_0),label = 'IR + IL')
    plt.legend()
    plt.show()
        
        
    return I_list_0, I_list_1, I_list_2



