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


# In[2]:


#import path to modules here.
# import sys
# sys.path.insert(0,'path_to_modules')


# #### Functions

# In[3]:




############################################################################################################################
#                             """" The following code calculates the current T = 0 K"""""
############################################################################################################################


def tags_list(n_leads = 2, lead_i = 1):
    
    """ Input:
       - n_leads = number of leads attached to the scattering region.
       - lead_i = label for lead i. Can take the possible values 0, 1, 2, ..., n_leads -1
       
       Returns:
        - list of numbers: [0,1,2,..], which does not contain 'lead_i'
    """
    list_tags = []
    if lead_i >= n_leads:
        print('Invalid ')
        return []
    for n in range(n_leads):

        if n!= lead_i:
            list_tags.append(n)
    return list_tags


# # $T_i= \sum_{j\neq i}T_{ij}$

# In[8]:


# n_leads: Number of leads connected to the scattering area.
# 

def Tji(systemf,lead_i =0,n_leads=4,emin=-3,emax=3, de=0.1,plot_var=False,check_hermiticity=True):
    """" 
    Input: 
        - systemf = finalized kwant system
        - leadi = label for a lead
        - n_leads = number of leads attached to the scattering region.
        - energy window os given by [emin,emax] and 'de' is the gridsize
        - plot_var = variable that allows for direct plotting of all transmissions
        - check_hermiticity = boolean, it checks weither the Hamiltonian matrix of the finalized system is hermetian 
    
    Output:
     - list of transmissions from lead j = {1,2,3.., nleads}, into lead i. For j !=i (!!important!!):
     that is:
        # [T1i_sublist, T2i_sublist,....,Tji_sublist] with i!=j
     - Each Tji_sublist is a list of transmissions from lead j to lead i for the energies in the energywindow [emin,emax].
     - List of corresponding transmission labels: [[1,i],[2,i],....[nleads-1,i]] (for j !=i )"""
    
    # Call finalized system
    system =  systemf
    
    #Energy window
    
    
    
    if emin <= emax:
        energies=np.arange(emin,emax+de,de)
    if emin > emax:
        n = abs(int((emax-emin)/de)) + 1
        energies = np.linspace(emin,emax,n)
    
    #labels of the leads for which we want to calculate the transmission
    leads_j_tags = tags_list(n_leads, lead_i)
    
    # We need transmissions from {nleads} into i, but not for nlead = i
    j_leads = n_leads - 1
    Tji_list = [ [] for _ in range(j_leads) ]
    tag_list_ji = []
    
    
    for n in range(len(leads_j_tags)):
        k = leads_j_tags[n]
        tag_list_ji.append([k,lead_i])
        
        
    for energy in energies:
        smatrix_m = kwant.smatrix(system, energy,check_hermiticity=check_hermiticity)
        
        
        
        for n in range(len(leads_j_tags)):
            
            #For every lead k, 
            k = leads_j_tags[n]
            
            #the transmission from lead k to lead i is calculated
            Tki = smatrix_m.transmission(lead_i,k)
            
            # and added to the list:
            Tji_list[n].append(Tki)
            
    
            
           
            
            
                
    
       
    return energies, Tji_list,tag_list_ji

# n_leads: Number of leads connected to the scattering area.
# returns list of transmissions from lead i into leads j = {1,2,3.., nleads}. !!! With j !=i:
    # [Ti1_sublist, Ti2_sublist,....,Tij_sublist] met i!=j

def Tij(systemf,lead_i =0,n_leads=4,emin=-3,emax=3, de=0.1,plot_var=False,check_hermiticity=True):
    
    
    """" 
    Input: 
    - systemf = finalized kwant system
    - leadi = label for a lead
    - n_leads = number of leads attached to the scattering region.
    - energy window os given by [emin,emax] and 'de' is the gridsize
    - plot_var = variable that allows for direct plotting of all transmissions
    - check_hermiticity = boolean. Checks weither the Hamiltonian matrix of the finalized system is hermetian 
    
    Ouput:
     - list of transmissions from lead j = {1,2,3.., nleads}, into lead i. For j !=i (!!important!!):
       [Ti1_sublist, Ti2_sublist,....,Tij_sublist] with i!=j. 
    - Each Tij_sublist is a list of transmissions from lead i to lead j for the energies in the energywindow [emin,emax].
    
    - List of corresponding transmission labels: [[i,1],[i,2],....[i,nleads-1]] (for j !=i ) """

    # Call finalized system
    system =  systemf
    
    #Energy window
    
    
    
    
#     n = abs(int((emax-emin)/de)) + 1
#     energies = np.linspace(emin,emax,n)
    if emin <= emax:
        energies=np.arange(emin,emax+de,de)
    if emin > emax:
        n = abs(int((emax-emin)/de)) + 1
        energies = np.linspace(emin,emax,n)
    
    #labels of the leads for which we want to calculate the transmission
    leads_j_tags = tags_list(n_leads, lead_i)
    
    # We need transmissions from {nleads} into i, but not for nlead = i
    j_leads = n_leads - 1
    Tij_list = [ [] for _ in range(j_leads) ]
    tag_list_ij = [ ]
    
    
    for n in range(len(leads_j_tags)):
        k = leads_j_tags[n]
        tag_list_ij.append([lead_i,k])
     
    for energy in energies:
        smatrix_m = kwant.smatrix(system, energy,check_hermiticity=check_hermiticity)
        
        
        
        for n in range(len(leads_j_tags)):
            
            #For every lead k, 
            k = leads_j_tags[n]
            
            #the transmission from lead k to lead i is calculated
            Tik = smatrix_m.transmission(k,lead_i)
            
            # and added to the list:
            Tij_list[n].append(Tik)
            

       
    return energies, Tij_list ,tag_list_ij


# # Transmission as a function of energy

# # $T_{ij}(E)$

# In[11]:


# transmission from left/right into buttiker probe for a certain energy 
def Tij_energy(energy,leadi,leadj,
                systemf,check_hermiticity=True):
    
    """ 
    Input:
    - energy = the energy for which we want to evaluate the transmission
    - leadi = label for initial lead
    - leadj = label for final lead 
    - systemf = finalzed kwant system
    - check_hermiticity = boolean.  Checks weither the Hamiltonian matrix of the finalized system is hermetian  


    Ouput:
    Transmission from lead i to lead j at for energy (returns a number)
    """

    #System created
    system = systemf
    
    #Transmission Tif from lead initial i to lead final (j) 
    smatrix_m = kwant.smatrix(system, energy,check_hermiticity=check_hermiticity)
    Tij = smatrix_m.transmission(leadj,leadi)
    

    return Tij




# # $T_{ji}(E)$

# In[12]:


def Tji_energy(energy,leadj,leadi,
                systemf,check_hermiticity=True):
    
    """ 
    Input:
    - energy = the energy for which we want to evaluate the transmission
    - leadj = label for intial lead 
    - leadi = label for final lead
    - systemf = finalzed kwant system
    - check_hermiticity = boolean.  Checks weither the Hamiltonian matrix of the finalized system is hermetian  


    Ouput:
    Transmission from lead j to lead i at for energy (returns a number)
    """
     
    
    system = systemf
    
    #Transmission Tfi from lead j to  lead i
    smatrix_m = kwant.smatrix(system, energy,check_hermiticity=check_hermiticity)
    Tji = smatrix_m.transmission(leadi,leadj)

    
    
    return Tji


# # Numerical integration: Trapezoid method

# In[13]:


# Numerical integration of transmission over the energy is done with:
# IMPORTANT: We need transmission T(E) = 0 for E <= energies[0], E >= energies[imax]
# otherwise integrals with upper bounds b > energies[imax], 
# and lower bounds a < energies[0] give wrong results.
# The code in this notbook, is designed s.t. this is always True. But please check it when making another code.


# In[14]:
#plot transmission for lead i:
def plot_transmission_i(lead_i,list_tags_transmissions,energies):
    
    """
    Input:
    - lead_i = label of lead i
    - list_tags_transmissions = list of possible transmissions into lead i for all possible leads of the system 
    - energies = list of energies for which the transmissions list_tags_transmissions are calculated
    
    Output:
    - plot of all transmission into lead i
    """
    
    Tij_list_i,tag_list_ij_i, Tji_list_i,tag_list_ji_i = list_tags_transmissions[lead_i]


    for i in range(len(tag_list_ij_i)):
        plt.plot(energies,Tij_list_i[i],label = 'T{}'.format(tag_list_ij_i[i]))

    for i in range(len(tag_list_ji_i)):
        plt.plot(energies,Tji_list_i[i],label = 'T{}'.format(tag_list_ji_i[i]))


    plt.legend()
    plt.show()
    
    
# Module for calculating integerals:

import Trapezoid_Integration


def plot_reciprocity_i(lead_i,list_tags_transmissionsP,list_tags_transmissionsM,energies):
    
    TijP_list_i,tagP_list_ij_i, TjiP_list_i,tagP_list_ji_i = list_tags_transmissionsP[lead_i]
    TijM_list_i,tagM_list_ij_i, TjiM_list_i,tagM_list_ji_i = list_tags_transmissionsM[lead_i]

    for i in range(len(tagP_list_ij_i)):
        plt.plot(energies,np.subtract(TijP_list_i[i],TjiM_list_i[i]),label = 'TP{}-TM{}'.format(tagP_list_ij_i[i],tagM_list_ji_i[i]))

    for i in range(len(tagP_list_ji_i)):
        plt.plot(energies,np.subtract(TjiP_list_i[i],TijM_list_i[i]),label = 'TP{}-TM{}'.format(tagP_list_ji_i[i],tagM_list_ij_i[i]))


    plt.legend()
    plt.show()

import Trapezoid_Integration


# ## Current: $ I_i$

# In[15]:


# mu_j_list = list of chemical potential corresponding to the leads in lead T_ij_list_n
# mu_i = chemical potential of lead i.
# energies = list of energies to be integrated over
# systemf = finalized Kwant system.


def Ii_current(mu_i,mu_j_list,
               tag_list_ij,tag_list_ji,
               energies,
               Tij_list_n,Tji_list_n,
             systemf,nleads):
    
    """ 
    Input:
    - mu_i = chemical potential of lead i
    - mu_j_list = list of chemical potentials all other leads that are not lead i.
    - energies = list of energies over which the integral is performed
    - Tij_list_n = list of transmission from lead i to lead j.
    - Tji_list_n = list of transmission from lead j to lead i.
    
    Ouput:
    - current into lead i calculated with the Landauer Buttiker formula which does not assume TRS.
        #sum_{j \neq i} fi Tij -   sum_{j \neq i} fj Tji
    """

    if len(Tij_list_n)+1 != nleads:
        print('Check number of leads',print(len(Tij_list_n)+1),nleads)
    
    # Current into left lead is given by the integral over:
    #   sum_{j \neq i} fi Tij -   sum_{j \neq i} fj Tji
    
    
    # --- POSITIVE SIGN ---
    #   sum_{j \neq i} fi Tij
    Iij_list = []
    for n in range(len(tag_list_ij)):
        Tij_list = Tij_list_n[n]
        tags = tag_list_ij[n]
        leadi = tags[0]
        leadj = tags[1]
        
        
        Iij_new = Trapezoid_Integration.calculate_integral(mu_i,-10**8,energies,Tij_list,Tij_energy,
                                                           leadi,leadj, systemf)
        
        Iij_list.append(Iij_new)
        
        
    # --- NEGATIVE SIGN ---
    ## sum_{j \neq i} fj Tji
    Iji_list = []
    
    for n in range(len(tag_list_ji)):
        Tji_list = Tji_list_n[n]
        tags = tag_list_ji[n]
        
        
        muj = mu_j_list[n]
        
        leadj = tags[0]
        leadi = tags[1]
        

        
        Iji_new = Trapezoid_Integration.calculate_integral(muj,-10**8,energies,Tji_list,Tji_energy,
                                                           leadj,leadi, systemf)
        
        Iji_list.append(Iji_new)
    
    Icur_total = sum(Iij_list) - sum(Iji_list)
   
    return Icur_total


# # Chemical potentials

# #### List of chemical potentials for lead 0,1,2,... respectively.

# In[ ]:



def mu_list_generate(muB_list,ef,V,nleads,rho_L=1/2,rho_R=1/2):
    
    """ 
    Input
    - ef = fermi energy
    - V = bias voltage
    - nleads = number of leads attached to the scattering region.
    - muB_list = list of chemical potentials which for lead 2,3,4.. (the Buttiker probes)
    Ouput:
    - list of chemical potentials for all leads in ascending order:
      [mu0,mu1,mu2,...]
    
    """
    
    mu_list = []
    
    # The left,right lead have labels 0,1 respectively.
    mu0 = ef + V*rho_L
    mu1 = ef - V*rho_R
    
    mu_list.append(mu0)
    mu_list.append(mu1)
    
    for n in range(len(muB_list)):
        mu_list.append(muB_list[n])
        
    
    return mu_list


# # Current into Buttiker porbe

# In[18]:


def currents_B_probe(muB_list,ef,V,nleads,list_tags_transmissions,energies,systf,print_bool=False,rho_L=1/2,rho_R=1/2):
    
    
    """ 
    Input:
    - muB_list = list of chemical potentials which for lead 2,3,4.. (the Buttiker probes)
    - ef = fermi energy
    - V = bias voltage
    - energies = list of energies over which the integral is performed
    - list_tags_transmissions = list of possible transmissions into lead i for all possible leads of the system 
    - energies = list of energies for which the transmissions list_tags_transmissions are calculated
    - systf = finalzed kwant system
    - print_bool= boolean. If evaluated to True, then prints Output.
    Ouput:
    - list of currents into leads: 2,3,..:
        [I2,I3,...]
    """


    
    
    
    list_currents_B_probe = []
    
    
    # We only have Buttiker probes for n>=2
    for lead_i in range(2,nleads):
        
        #Pick out transmission for lead_i
        Tij_list,tag_list_ij, Tji_list,tag_list_ji = list_tags_transmissions[lead_i]
        
        # Generate mus for voltage V
        mu_list_object = mu_list_generate(muB_list,ef,V,nleads,rho_L,rho_R)

        # Generate tags j for lead i s.t. the chemical potential muj corresponds to the correct transmission
        tagsi_j_list = tags_list(nleads, lead_i)

        #mus for lead i
        mui_j_list = [mu_list_object[n] for n in tagsi_j_list]
        mui = mu_list_object[lead_i]

        Ii_cur = Ii_current(mui,mui_j_list,
                   tag_list_ij,tag_list_ji,
                   energies,
                   Tij_list,Tji_list,
                   systf,nleads)
        
        list_currents_B_probe.append(Ii_cur)
    
    if print_bool == True:
        print(list_currents_B_probe)
        
        
    return list_currents_B_probe




def currents_leadi_probe(lead_i,muB_list,ef,V,nleads,energies,list_tags_transmissions,systf,print_bool=False,
                         rho_L=1/2,rho_R=1/2):
    
    
    """ 
    Input:
    - lead_i = label for lead i.
    - muB_list = list of chemical potentials which for lead 2,3,4.. (the Buttiker probes)
    - ef = fermi energy
    - V = bias voltage
    - energies = list of energies over which the integral is performed
    - list_tags_transmissions = list of possible transmissions into lead i for all possible leads of the system 
    - energies = list of energies for which the transmissions list_tags_transmissions are calculated
    - systf = finalzed kwant system
    - print_bool= boolean. If evaluated to True, then prints Output.
    
    Ouput:
    - current into lead i:
    """
    
    Tij_list,tag_list_ij, Tji_list,tag_list_ji = list_tags_transmissions[lead_i]
    
    # generate chemical potentials for lead 0,1,2,...
    mu_list_object = mu_list_generate(muB_list,ef,V,nleads,rho_L,rho_R)


    # Generate tags j for lead i s.t. the chemical potential muj corresponds to the correct transmission
    tagsi_j_list = tags_list(nleads, lead_i)
    
    # Relevant mus for voltage V
    mui_j_list = [mu_list_object[n] for n in tagsi_j_list]
    mui = mu_list_object[lead_i]
    
    # Current into lead i  
    Ii_cur = Ii_current(mui,mui_j_list,
               tag_list_ij,tag_list_ji,
               energies,
               Tij_list,Tji_list,
               systf,nleads)
    


    return Ii_cur



## Current for finite temperature.


# Fermi - Dirac Function


def fermi_dirac(beta,mui,energy):
    if beta < 10**9:
        fd = 1/(np.exp(beta*(energy-mui) ) + 1 )
        return fd
    if beta >= 10**9:
        fd = np.heaviside(-(energy-mui),1) 
        return fd

# Weighted transmission

### $\overline{T}_{ij}$

# transmission from left/right into buttiker probe for a certain energy 
def Tbar_ij_energy(energy,beta,mu,leadi,leadj,
                systemf,check_hermiticity=True):

    #System created
    system = systemf
    
    #Transmission Tif from lead initial i to lead final (j) 
    smatrix_m = kwant.smatrix(system, energy,check_hermiticity=check_hermiticity)
    Tbarij = smatrix_m.transmission(leadj,leadi)*fermi_dirac(beta,mu,energy)
    

    return Tbarij

### $\overline{T}_{ji}$


def Tbar_ji_energy(energy,beta,mu,leadj,leadi,
                systemf,check_hermiticity=True):
     
    
    system = systemf
    
    #Transmission Tfi from lead j to  lead i
    smatrix_m = kwant.smatrix(system, energy,check_hermiticity=check_hermiticity)
    Tji = smatrix_m.transmission(leadi,leadj)*fermi_dirac(beta,mu,energy)

    
    
    return Tji


# Weigthing function



def weighted_transmissions(lead_i,nleads,
                           beta,ef,V,muB_list,
                          list_tags_transmissions,energies):
    
    """Input: a) list transmissions from lead i to lead j (Tij) and list of transmissions from j to i Tji 
                necessary to calculate the current into lead i
              b) Parameters V,ef,muB_list deterimine the chemical potentials of the leads 
              c) nleads: number of leads
              d) lead_i: A label for the lead. Says which transmissions Tij, Tji are converted
              e) beta: 1/(kB*T) with T the temperature of all the leads"""
             
    """ Ouput: - Transmissions that are weigthed with the fermi-dirac function at finite temperature."""
    
    Tij_list_i,tag_list_ij_i, Tji_list_i,tag_list_ji_i = list_tags_transmissions[lead_i]
    leads_j_tags = tags_list(nleads, lead_i)                                                                               


    j_leads = nleads - 1
    Tbar_ij_list = [ [] for _ in range(j_leads) ]
    Tbar_ji_list = [ [] for _ in range(j_leads) ]

    # tags for the leads that j!= i
    leads_j_tags = tags_list(nleads, lead_i)

    # generate chemical potentials for lead 0,1,2,...
    mu_list_object = mu_list_generate(muB_list,ef,V,nleads)

#     # Chemical potential for lead_i
#     mu_i = mu_list_object[lead_i]


    for n in range(len(energies)):


        energy = energies[n]

        # Calculate fermi-dirac functions for lead 0,1,2,...
        fd_list = [ fermi_dirac(beta,mu_i,energy) for mu_i in mu_list_object]

        for j in range(j_leads):

            #The transmission Tij are weigthed with FD of lead i
            Tbarij_i = Tij_list_i[j][n]*fd_list[lead_i]
            Tbar_ij_list[j].append(Tbarij_i)


        for k in range(j_leads):

            # tag for lead_j
            jtag = leads_j_tags[k]

            #The transmission Tji are weigthed with FD of lead j
            Tbarji_i = Tji_list_i[k][n]*fd_list[jtag]
            Tbar_ji_list[k].append(Tbarji_i)


    list_transmissions = [Tbar_ij_list,
                          tag_list_ij_i, 
                          Tbar_ji_list,tag_list_ji_i]
    
    
    return list_transmissions



############################################################################################################################
#                         """" The following code calculates the charge current at finite temperature."""""
############################################################################################################################
# Weighted current

# mu_j_list = list of chemical potential corresponding to the leads in lead T_ij_list_n
# mu_i = chemical potential of lead i.
# energies = list of energies to be integrated over
# systemf = finalized Kwant system.


def Ibari_current(beta,mu_i,mu_j_list,
               tag_list_ij,tag_list_ji,
               energies,
               Tij_list_n,Tji_list_n,
             systemf,nleads):

    if len(Tij_list_n) + 1 != nleads:
        print('Check number of leads',print(len(Tij_list_n)+1),nleads)
    
    # Current into left lead is given by the integral over:
    #   sum_{j \neq i} fi Tij -   sum_{j \neq i} fj Tji
    
    
    # --- POSITIVE ---
    #   sum_{j \neq i} fi Tij
    Iij_list = []
    for n in range(len(tag_list_ij)):
        Tij_list = Tij_list_n[n]
        tags = tag_list_ij[n]
        leadi = tags[0]
        leadj = tags[1]
        
        
        Iij_new = Trapezoid_Integration.calculate_integral(10**8,-10**8,energies,Tij_list,Tbar_ij_energy,
                                                           beta,mu_i,
                                                           leadi,leadj, systemf)
        
        Iij_list.append(Iij_new)
        
        
    
    ## sum_{j \neq i} fj Tji
    Iji_list = []
    
    for n in range(len(tag_list_ji)):
        Tji_list = Tji_list_n[n]
        tags = tag_list_ji[n]
        
        
        muj = mu_j_list[n]
        
        leadj = tags[0]
        leadi = tags[1]
        

        
        Iji_new = Trapezoid_Integration.calculate_integral(10**8,-10**8,energies,Tji_list,Tbar_ji_energy,
                                                           beta,muj,
                                                           leadj,leadi, systemf)
        
        Iji_list.append(Iji_new)
    
    Ibar_cur_total = sum(Iij_list) - sum(Iji_list)
   
    return Ibar_cur_total



# $I_i(T,V) $ Current into lead $i$: finite temperature



def current_bar_leadi(lead_i,beta,muB_list,ef,V,nleads,energies,list_tags_transmissions,systf,print_bool=False):
    
    # Weigh transmissions with Fermi-Dirac functions
    Tij_list,tag_list_ij, Tji_list,tag_list_ji = weighted_transmissions(lead_i,nleads,
                           beta,ef,V,muB_list,
                          list_tags_transmissions,energies)
    
    # generate chemical potentials for lead 0,1,2,...
    mu_list_object = mu_list_generate(muB_list,ef,V,nleads)


    # Generate tags j for lead i s.t. the chemical potential muj corresponds to the correct transmission
    tagsi_j_list = tags_list(nleads, lead_i)
    
    # Relevant mus for voltage V
    mui_j_list = [ mu_list_object[n] for n in tagsi_j_list]
    mui = mu_list_object[lead_i]
    
    # Current into lead i  
    Ii_cur = Ibari_current(beta,mui,mui_j_list,
               tag_list_ij,tag_list_ji,
               energies,
               Tij_list,Tji_list,
               systf,nleads)
    


    return Ii_cur



def currents_bar_B_probe(muB_list,beta,ef,V,nleads,list_tags_transmissions,energies,systf,print_bool=False):
    
    list_currents_bar_B_probe = []
    
    
   
    
    
    # We only have Buttiker probes for n>=2
    for lead_i in range(2,nleads):
        
        #Pick out WEIGHTED transmissions for lead_i
        Tij_list,tag_list_ij, Tji_list,tag_list_ji =  weighted_transmissions(lead_i,nleads,
                                                       beta,ef,V,muB_list,
                                                      list_tags_transmissions,energies)
        
        
        
        # Generate mus for voltage V
        mu_list_object = mu_list_generate(muB_list,ef,V,nleads)

        # Generate tags j for lead i s.t. the chemical potential muj corresponds to the correct transmission
        tagsi_j_list = tags_list(nleads, lead_i)

        #mus for lead i
        mui_j_list = [mu_list_object[n] for n in tagsi_j_list]
        mui = mu_list_object[lead_i]

        Ii_cur =  Ibari_current(beta,mui,mui_j_list,
                   tag_list_ij,tag_list_ji,
                   energies,
                   Tij_list,Tji_list,
                   systf,nleads)
        
        list_currents_bar_B_probe.append(Ii_cur)
    
    if print_bool == True:
        print(list_currents_bar_B_probe)
        
    return list_currents_bar_B_probe


def current_allmus(mu_list,lead_i,nleads,energies,list_tags_transmissions,systf,print_bool=False):
    
    Tij_list,tag_list_ij, Tji_list,tag_list_ji = list_tags_transmissions[lead_i]
    
    # generate chemical potentials for lead 0,1,2,...
    # mu_list = [mu0,mu1,mu2,....]


    # Generate tags j for lead i s.t. the chemical potential muj corresponds to the correct transmission
    tagsi_j_list = tags_list(nleads, lead_i)
    
    # Relevant mus for voltage V
    mui_j_list = [mu_list[n] for n in tagsi_j_list]
    mui = mu_list[lead_i]
    
    # Current into lead i  
    Ii_cur = Ii_current(mui,mui_j_list,
               tag_list_ij,tag_list_ji,
               energies,
               Tij_list,Tji_list,
               systf,nleads)
    


    return Ii_cur
