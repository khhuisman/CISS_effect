#!/usr/bin/env python
# coding: utf-8

# # Plot functions
# - To plot spin-polarization
# - To plot transmission 

# In[2]:


# syst: finalized kwant system with leads.
# emin,emax = minimum and maximum energy.
# de = the stepsize for each point plotted.
# e_steps = number of points plotted.


# ### Total conductance

# Plot the total conductance/transmission: $T_{LR} e^2/hbar = G_{LR}$

# In[1]:


import kwant
import math
import cmath
from matplotlib import pyplot as plt
import numpy as np

def plot_cond_total(syst, emin, emax, de,check_hermiticity=True):
    
    dataG = []
    energies = np.arange(emin,emax,de)

    for energy in energies:
            
            smatrix = kwant.smatrix(syst, energy,check_hermiticity =check_hermiticity)


            ## Conductance with kwant function transmission
            G = smatrix.transmission(0,1)
            

            # Add them to list
            dataG.append(G)
        

   
    plt.plot(energies, dataG, color ='blue' )
    plt.xlabel("E [eV] ")
    plt.ylabel("Transmission  ")
    plt.tight_layout()
    plt.show()
    
    return energies, dataG

def plot_cond_total_LR(syst, emin, emax, de,plot_var=False):
    
    dataG = []
    energies = np.arange(emin,emax+de,de)

    for energy in energies:
            
            smatrix = kwant.smatrix(syst, energy)

             # Transmission from lead L (0) to lead R (1)

            TLR = smatrix.transmission(1,0)
            

            # Add them to list
            dataG.append(TLR)
        

   
    if plot_var == True:
        plt.plot(energies, dataG, color ='blue' )
        plt.xlabel("E [eV] ")
        plt.ylabel("Transmission  ")
        plt.tight_layout()
        plt.show()
    
    return energies, dataG


def plot_cond_total_RL(syst, emin, emax, de,plot_var=False):
    
    dataG = []
    energies = np.arange(emin,emax+de,de)

    for energy in energies:
            
            smatrix = kwant.smatrix(syst, energy)

             # Transmission from lead R (1) to lead L (0)

            TRL = smatrix.transmission(0,1)
            

            # Add them to list
            dataG.append(TRL)
        

   
    if plot_var == True:
        plt.plot(energies, dataG, color ='blue' )
        plt.xlabel("E [eV] ")
        plt.ylabel("Transmission  ")
        plt.tight_layout()
        plt.show()
    
    return energies, dataG


def plot_cond_total_BL(syst, emin, emax, de,plot_var=False):
    
    dataG = []
    energies = np.arange(emin,emax+de,de)

    for energy in energies:
            
            smatrix = kwant.smatrix(syst, energy)

            ## # Transmission from lead B (2) to lead L (0)
            TBL = smatrix.transmission(0,2)
            

            # Add them to list
            dataG.append(TBL)
        

   
    if plot_var == True:
        plt.plot(energies, dataG, color ='blue' )
        plt.xlabel("E [eV] ")
        plt.ylabel("Transmission  ")
        plt.tight_layout()
        plt.show()
    
    return energies, dataG


def plot_cond_total_LB(syst, emin, emax, de,plot_var=False):
    
    dataG = []
    energies = np.arange(emin,emax+de,de)

    for energy in energies:
            
            smatrix = kwant.smatrix(syst, energy)

            ## # Transmission from lead L (0) to lead B (2)
            TLB = smatrix.transmission(2,0)
            

            # Add them to list
            dataG.append(TLB)
        

   
    if plot_var == True:
        plt.plot(energies, dataG, color ='blue' )
        plt.xlabel("E [eV] ")
        plt.ylabel("Transmission  ")
        plt.tight_layout()
        plt.show()
    
    return energies, dataG

def plot_cond_total_RB(syst, emin, emax, de,plot_var=False):
    
    dataG = []
    energies = np.arange(emin,emax+de,de)

    for energy in energies:
            
            smatrix = kwant.smatrix(syst, energy)

            ## # Transmission from lead R (1) to lead B (2)
            TRB = smatrix.transmission(2,1)
            

            # Add them to list
            dataG.append(TRB)
        

    if plot_var == True:
        plt.plot(energies, dataG, color ='blue' )
        plt.xlabel("E [eV] ")
        plt.ylabel("Transmission  ")
        plt.tight_layout()
        plt.show()
    
    return energies, dataG

def plot_cond_total_noplot(syst,energies,check_hermiticity=True):
    
    
    
    dataG = []
    

    for energy in energies:
            
            smatrix = kwant.smatrix(syst, energy,check_hermiticity=check_hermiticity)

            ## Conductance with kwant function transmission
            G = smatrix.transmission(0,1)
            

            # Add them to list
            dataG.append(G)

    
    return energies, dataG

def plot_cond_total_noplot_reversed(syst,energies,check_hermiticity=True):
    
    
    
    dataG = []
    

    for energy in energies:
            
            smatrix = kwant.smatrix(syst, energy,check_hermiticity=check_hermiticity)

            ## Conductance with kwant function transmission
            G = smatrix.transmission(1,0)
            

            # Add them to list
            dataG.append(G)

    
    return energies, dataG

def plot_reflection_total_noplot_00(syst,energies,check_hermiticity=True):
    
    
    
    dataG = []
    

    for energy in energies:
            
            smatrix = kwant.smatrix(syst, energy,check_hermiticity=check_hermiticity)

            ## Conductance with kwant function transmission
            R00 = smatrix.transmission(0,0)
            

            # Add them to list
            dataG.append(R00)

    
    return energies, dataG

def calc_cond_total(syst, emin, emax, de):
    
    dataG = []
    energies = np.arange(emin,emax,de)

    for energy in energies:
            
            smatrix = kwant.smatrix(syst, energy)

            ## Conductance with kwant function transmission
            G = smatrix.transmission(0,1)
            

            # Add them to list
            dataG.append(G)
        
    
    return energies, dataG



# ### Spin Conductance

# Plot transmission components: $T_{LR}^{\uparrow \uparrow}, T_{LR}^{\downarrow \downarrow}$ and $T_{LR}^{\uparrow \downarrow},  T_{LR}^{\downarrow \uparrow}$ as a function of energy

# In[1]:


import kwant
import math
import cmath
from matplotlib import pyplot
import numpy

def plot_conductances_LR(syst, emin, emax, de,plot_transmission=True):
    
    dataG00 = []
    dataG11 = []
    dataG01 = []
    dataG10 = []
    
    energies = np.arange(emin,emax,de)
    
    
    # Function below calculates:  
    # transmission from left lead (0) with spin b,
    # too the right (1) lead with spin a
    
    
    def TLR_ab(smatrix,b,a):
        Sm = smatrix
        Gab = Sm.transmission((1,a),(0,b))
        return Gab

    
    
    for energy in energies:
        
            
            smatrix = kwant.smatrix(syst, energy)
            

            ## Conductance with kwant function transmission
            G00 = TLR_ab(smatrix,0,0)
            G01 = TLR_ab(smatrix,0,1)
            G10 = TLR_ab(smatrix,1,0)
            G11 = TLR_ab(smatrix,1,1)

            # Add them to list
            dataG00.append(G00)
            dataG11.append(G11)
            dataG01.append(G01)
            dataG10.append(G10)
    
    if plot_transmission == True:
        plt.subplots(figsize=(10, 5))
        #Plot G00, G11
        sub1= pyplot.subplot(1, 2, 1)
        sub1.set_title(' G00(blue,solid), G11(red, dashed) ')
        plt.plot(energies, dataG00, color = 'blue' )
        plt.plot(energies, dataG11, color = 'red', linestyle = '--' )
        plt.xlabel("energy [t]")
        plt.ylabel("Conductance e^2/hbar")
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))


        #Plot G01, G10
        sub2 = pyplot.subplot(1, 2, 2)
        sub2.set_title('G01 (blue), G10 (dashed)')
        plt.plot(energies, dataG01, color ='blue' )
        plt.plot(energies, dataG10, color = 'red', linestyle = '--' , markersize = 20)
        plt.xlabel("energy [t]")
        plt.ylabel("Conductance e^2/hbar   ")
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        plt.tight_layout()
        plt.show()
    
    return energies, dataG00,dataG01,dataG10,dataG11

def plot_conductances_RL(syst, emin, emax, de,plot_transmission=True):
    
    dataG00 = []
    dataG11 = []
    dataG01 = []
    dataG10 = []
    
    energies = np.arange(emin,emax,de)
    
    
    # Function below calculates:  
    # transmission from right lead (1) with spin b,
    # too the left (0) lead with spin a
    
    
    def TRL_ab(smatrix,b,a):
        Sm = smatrix
        Gab = Sm.transmission((0,a),(1,b))
        return Gab

    
    
    for energy in energies:
        
            
    
            
            smatrix = kwant.smatrix(syst, energy)
            

            ## Conductance with kwant function transmission
            G00 = TRL_ab(smatrix,0,0)
            G01 = TRL_ab(smatrix,0,1)
            G10 = TRL_ab(smatrix,1,0)
            G11 = TRL_ab(smatrix,1,1)

            # Add them to list
            dataG00.append(G00)
            dataG11.append(G11)
            dataG01.append(G01)
            dataG10.append(G10)
    
    if plot_transmission == True:
        plt.subplots(figsize=(10, 5))
        #Plot G00, G11
        sub1= pyplot.subplot(1, 2, 1)
        sub1.set_title(' G00(blue,solid), G11(red, dashed) ')
        plt.plot(energies, dataG00, color = 'blue' )
        plt.plot(energies, dataG11, color = 'red', linestyle = '--' )
        plt.xlabel("energy [t]")
        plt.ylabel("Conductance e^2/hbar")
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))


        #Plot G01, G10
        sub2 = pyplot.subplot(1, 2, 2)
        sub2.set_title('G01 (blue), G10 (dashed)')
        plt.plot(energies, dataG01, color ='blue' )
        plt.plot(energies, dataG10, color = 'red', linestyle = '--' , markersize = 20)
        plt.xlabel("energy [t]")
        plt.ylabel("Conductance e^2/hbar   ")
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        plt.tight_layout()
        plt.show()
    
    return energies, dataG00,dataG01,dataG10,dataG11



def plot_conductances_LB(syst, emin, emax, de,plot_transmission=True):
    
    dataG00 = []
    dataG11 = []
    dataG01 = []
    dataG10 = []
    
    energies = np.arange(emin,emax,de)
    
    
    # Function below calculates:  
    # transmission from left lead (0) with spin b,
    # too the Buttiker probe (2) lead with spin a
    
    
    def TLB_ab(smatrix,b,a):
        Sm = smatrix
        Gab = Sm.transmission((2,a),(0,b))
        return Gab

    
    
    for energy in energies:
        
            
    
            
            smatrix = kwant.smatrix(syst, energy)
            

            ## Conductance with kwant function transmission
            G00 = TLB_ab(smatrix,0,0)
            G01 = TLB_ab(smatrix,0,1)
            G10 = TLB_ab(smatrix,1,0)
            G11 = TLB_ab(smatrix,1,1)

            # Add them to list
            dataG00.append(G00)
            dataG11.append(G11)
            dataG01.append(G01)
            dataG10.append(G10)
    
    if plot_transmission == True:
        plt.subplots(figsize=(10, 5))
        #Plot G00, G11
        sub1= pyplot.subplot(1, 2, 1)
        sub1.set_title(' G00(blue,solid), G11(red, dashed) ')
        plt.plot(energies, dataG00, color = 'blue' )
        plt.plot(energies, dataG11, color = 'red', linestyle = '--' )
        plt.xlabel("energy [t]")
        plt.ylabel("Conductance e^2/hbar")
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))


        #Plot G01, G10
        sub2 = pyplot.subplot(1, 2, 2)
        sub2.set_title('G01 (blue), G10 (dashed)')
        plt.plot(energies, dataG01, color ='blue' )
        plt.plot(energies, dataG10, color = 'red', linestyle = '--' , markersize = 20)
        plt.xlabel("energy [t]")
        plt.ylabel("Conductance e^2/hbar   ")
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        plt.tight_layout()
        plt.show()
    
    return energies, dataG00,dataG01,dataG10,dataG11


def plot_conductances_BL(syst, emin, emax, de,plot_transmission=True):
    
    dataG00 = []
    dataG11 = []
    dataG01 = []
    dataG10 = []
    
    energies = np.arange(emin,emax,de)
    
    
    # Function below calculates:  
    # transmission from buttiker lead (2) with spin b,
    # too the left lead probe (0) lead with spin a
    
    
    def TBL_ab(smatrix,b,a):
        Sm = smatrix
        Gab = Sm.transmission((0,a),(2,b))
        return Gab

    
    
    for energy in energies:
        
            
    
            
            smatrix = kwant.smatrix(syst, energy)
            

            ## Conductance with kwant function transmission
            G00 = TBL_ab(smatrix,0,0)
            G01 = TBL_ab(smatrix,0,1)
            G10 = TBL_ab(smatrix,1,0)
            G11 = TBL_ab(smatrix,1,1)

            # Add them to list
            dataG00.append(G00)
            dataG11.append(G11)
            dataG01.append(G01)
            dataG10.append(G10)
    
    if plot_transmission == True:
        plt.subplots(figsize=(10, 5))
        #Plot G00, G11
        sub1= pyplot.subplot(1, 2, 1)
        sub1.set_title(' G00(blue,solid), G11(red, dashed) ')
        plt.plot(energies, dataG00, color = 'blue' )
        plt.plot(energies, dataG11, color = 'red', linestyle = '--' )
        plt.xlabel("energy [t]")
        plt.ylabel("Conductance e^2/hbar")
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))


        #Plot G01, G10
        sub2 = pyplot.subplot(1, 2, 2)
        sub2.set_title('G01 (blue), G10 (dashed)')
        plt.plot(energies, dataG01, color ='blue' )
        plt.plot(energies, dataG10, color = 'red', linestyle = '--' , markersize = 20)
        plt.xlabel("energy [t]")
        plt.ylabel("Conductance e^2/hbar   ")
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        plt.tight_layout()
        plt.show()
    
    return energies, dataG00,dataG01,dataG10,dataG11


def plot_conductances_RR(syst, emin, emax, de,plot_transmission=True):
    
    dataG00 = []
    dataG11 = []
    dataG01 = []
    dataG10 = []
    
    energies = np.arange(emin,emax,de)
    
    
    # Function below calculates:  
    # transmission from right lead (0) with spin b,
    # too the right (0) lead with spin a
    
    
    def TRR_ab(smatrix,b,a):
        Sm = smatrix
        Gab = Sm.transmission((0,a),(0,b))
        return Gab

    
    
    for energy in energies:
        
            
    
            
            smatrix = kwant.smatrix(syst, energy)
            

            ## Conductance with kwant function transmission
            Guu = TRR_ab(smatrix,0,0)
            Gud = TRR_ab(smatrix,0,1)
            Gdu = TRR_ab(smatrix,1,0)
            Gdd = TRR_ab(smatrix,1,1)

            # Add them to list
            dataG00.append(Guu)
            dataG11.append(Gdd)
            dataG01.append(Gud)
            dataG10.append(Gdu)
    
    if plot_transmission == True:
        plt.subplots(figsize=(10, 5))
        #Plot G00, G11
        sub1= pyplot.subplot(1, 2, 1)
        sub1.set_title(' Ruu(blue,solid), Rdd(red, dashed) ')
        plt.plot(energies, dataG00, color = 'blue' )
        plt.plot(energies, dataG11, color = 'red', linestyle = '--' )
        plt.xlabel("energy [t]")
        plt.ylabel("Conductance e^2/hbar")
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))


        #Plot G01, G10
        sub2 = pyplot.subplot(1, 2, 2)
        sub2.set_title('Rud (blue), Rdu (dashed)')
        plt.plot(energies, dataG01, color ='blue' )
        plt.plot(energies, dataG10, color = 'red', linestyle = '--' , markersize = 20)
        plt.xlabel("energy [t]")
        plt.ylabel("Conductance e^2/hbar   ")
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        plt.tight_layout()
        plt.show()
    
    return energies, dataG00,dataG01,dataG10,dataG11

def plot_conductances_BB(syst, emin, emax, de,plot_transmission=True):
    
    dataG00 = []
    dataG11 = []
    dataG01 = []
    dataG10 = []
    
    energies = np.arange(emin,emax,de)
    
    
    # Function below calculates:  
    # transmission from buttiker lead (2) with spin b,
    # too the buttiker (2) lead with spin a
    
    
    def TBB_ab(smatrix,b,a):
        Sm = smatrix
        Gab = Sm.transmission((0,a),(0,b))
        return Gab

    
    
    for energy in energies:
        
            
    
            
            smatrix = kwant.smatrix(syst, energy)
            

            ## Conductance with kwant function transmission
            Guu = TBB_ab(smatrix,0,0)
            Gud = TBB_ab(smatrix,0,1)
            Gdu = TBB_ab(smatrix,1,0)
            Gdd = TBB_ab(smatrix,1,1)

            # Add them to list
            dataG00.append(Guu)
            dataG11.append(Gdd)
            dataG01.append(Gud)
            dataG10.append(Gdu)
    
    if plot_transmission == True:
        plt.subplots(figsize=(10, 5))
        #Plot G00, G11
        sub1= pyplot.subplot(1, 2, 1)
        sub1.set_title(' Ruu(blue,solid), Rdd(red, dashed) ')
        plt.plot(energies, dataG00, color = 'blue' )
        plt.plot(energies, dataG11, color = 'red', linestyle = '--' )
        plt.xlabel("energy [t]")
        plt.ylabel("Conductance e^2/hbar")
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))


        #Plot G01, G10
        sub2 = pyplot.subplot(1, 2, 2)
        sub2.set_title('Rud (blue), Rdu (dashed)')
        plt.plot(energies, dataG01, color ='blue' )
        plt.plot(energies, dataG10, color = 'red', linestyle = '--' , markersize = 20)
        plt.xlabel("energy [t]")
        plt.ylabel("Conductance e^2/hbar   ")
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        plt.tight_layout()
        plt.show()
    
    return energies, dataG00,dataG01,dataG10,dataG11





def plot_conductances_BR(syst, emin, emax, de,plot_transmission=True):
    
    dataG00 = []
    dataG11 = []
    dataG01 = []
    dataG10 = []
    
    energies = np.arange(emin,emax,de)
    
    
    # Function below calculates:  
    # transmission from buttiker lead (2) with spin b,
    # too the right lead (1) with spin a
    
    
    def TBR_ab(smatrix,b,a):
        Sm = smatrix
        Gab = Sm.transmission((1,a),(2,b))
        return Gab

    
    
    for energy in energies:
        
            
            smatrix = kwant.smatrix(syst, energy)
            

            ## Conductance with kwant function transmission
            G00 = TBR_ab(smatrix,0,0)
            G01 = TBR_ab(smatrix,0,1)
            G10 = TBR_ab(smatrix,1,0)
            G11 = TBR_ab(smatrix,1,1)

            # Add them to list
            dataG00.append(G00)
            dataG11.append(G11)
            dataG01.append(G01)
            dataG10.append(G10)
    
    if plot_transmission == True:
        plt.subplots(figsize=(10, 5))
        #Plot G00, G11
        sub1= pyplot.subplot(1, 2, 1)
        sub1.set_title(' G00(blue,solid), G11(red, dashed) ')
        plt.plot(energies, dataG00, color = 'blue' )
        plt.plot(energies, dataG11, color = 'red', linestyle = '--' )
        plt.xlabel("energy [t]")
        plt.ylabel("Conductance e^2/hbar")
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))


        #Plot G01, G10
        sub2 = pyplot.subplot(1, 2, 2)
        sub2.set_title('G01 (blue), G10 (dashed)')
        plt.plot(energies, dataG01, color ='blue' )
        plt.plot(energies, dataG10, color = 'red', linestyle = '--' , markersize = 20)
        plt.xlabel("energy [t]")
        plt.ylabel("Conductance e^2/hbar   ")
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        plt.tight_layout()
        plt.show()
    
    return energies, dataG00,dataG01,dataG10,dataG11

def plot_conductances_RB(syst, emin, emax, de,plot_transmission=True):
    
    dataG00 = []
    dataG11 = []
    dataG01 = []
    dataG10 = []
    
    energies = np.arange(emin,emax,de)
    
    
    # Function below calculates:  
    # transmission from right lead (1) with spin b,
    # too the buttiker lead (2) with spin a
    
    
    def TRB_ab(smatrix,b,a):
        Sm = smatrix
        Gab = Sm.transmission((2,a),(1,b))
        return Gab

    
    
    for energy in energies:
        
            
            smatrix = kwant.smatrix(syst, energy)
            

            ## Conductance with kwant function transmission
            G00 = TRB_ab(smatrix,0,0)
            G01 = TRB_ab(smatrix,0,1)
            G10 = TRB_ab(smatrix,1,0)
            G11 = TRB_ab(smatrix,1,1)

            # Add them to list
            dataG00.append(G00)
            dataG11.append(G11)
            dataG01.append(G01)
            dataG10.append(G10)
    
    if plot_transmission == True:
        plt.subplots(figsize=(10, 5))
        #Plot G00, G11
        sub1= pyplot.subplot(1, 2, 1)
        sub1.set_title(' G00(blue,solid), G11(red, dashed) ')
        plt.plot(energies, dataG00, color = 'blue' )
        plt.plot(energies, dataG11, color = 'red', linestyle = '--' )
        plt.xlabel("energy [t]")
        plt.ylabel("Conductance e^2/hbar")
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))


        #Plot G01, G10
        sub2 = pyplot.subplot(1, 2, 2)
        sub2.set_title('G01 (blue), G10 (dashed)')
        plt.plot(energies, dataG01, color ='blue' )
        plt.plot(energies, dataG10, color = 'red', linestyle = '--' , markersize = 20)
        plt.xlabel("energy [t]")
        plt.ylabel("Conductance e^2/hbar   ")
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        plt.tight_layout()
        plt.show()
    
    return energies, dataG00,dataG01,dataG10,dataG11



def plot_conductances_LL(syst, emin, emax, de,plot_transmission=True):
    
    dataG00 = []
    dataG11 = []
    dataG01 = []
    dataG10 = []
    
    energies = np.arange(emin,emax,de)
    
    
    # Function below calculates:  
    # transmission from left lead (1) with spin b,
    # too the left (1) lead with spin a
    
    
    def TLL_ab(smatrix,b,a):
        Sm = smatrix
        Gab = Sm.transmission((1,a),(1,b))
        return Gab

    
    
    for energy in energies:
        
            
            smatrix = kwant.smatrix(syst, energy)
            

            ## Conductance with kwant function transmission
            G00 = TLL_ab(smatrix,0,0)
            G01 = TLL_ab(smatrix,0,1)
            G10 = TLL_ab(smatrix,1,0)
            G11 = TLL_ab(smatrix,1,1)

            # Add them to list
            dataG00.append(G00)
            dataG11.append(G11)
            dataG01.append(G01)
            dataG10.append(G10)
    
    if plot_transmission == True:
        plt.subplots(figsize=(10, 5))
        #Plot G00, G11
        sub1= pyplot.subplot(1, 2, 1)
        sub1.set_title(' G00(blue,solid), G11(red, dashed) ')
        plt.plot(energies, dataG00, color = 'blue' )
        plt.plot(energies, dataG11, color = 'red', linestyle = '--' )
        plt.xlabel("energy [t]")
        plt.ylabel("Conductance e^2/hbar")
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))


        #Plot G01, G10
        sub2 = pyplot.subplot(1, 2, 2)
        sub2.set_title('G01 (blue), G10 (dashed)')
        plt.plot(energies, dataG01, color ='blue' )
        plt.plot(energies, dataG10, color = 'red', linestyle = '--' , markersize = 20)
        plt.xlabel("energy [t]")
        plt.ylabel("Conductance e^2/hbar   ")
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        plt.tight_layout()
        plt.show()
    
    return energies, dataG00,dataG01,dataG10,dataG11
# ### Spin Polarization

# Plot the spin polarization as a function of energy.

#  $ P =\frac{ T_{LR}^{\uparrow \uparrow} - T_{LR}^{\uparrow \downarrow} + T_{LR}^{\downarrow \uparrow} -T_{LR}^{\downarrow \downarrow}}{T_{LR}}$ with 
#  $T_{LR} = T_{LR}^{\uparrow \uparrow} +  T_{LR}^{\uparrow \downarrow} + T_{LR}^{\downarrow \uparrow} + T_{LR}^{\downarrow \downarrow}$

# In[4]:


def plot_pol(syst, emin, emax, de,check_hermiticity=True):
    
    dataP = []
    Tp_list = []
    energies = np.arange(emin,emax,de)

    # Function below calculates:  
    # transmission from left lead (1) with spin b,
    # too the right (0) lead with spin a
    def GLR_ba(smatrix,a,b):
        Gab = smatrix.transmission((0,a),(1,b))
        return Gab

    for energy in energies:
        

        smatrix = kwant.smatrix(syst, energy,check_hermiticity=check_hermiticity)

        # Total transmission:
        G_total = smatrix.transmission(0,1)

        ## Conductance with kwant function transmission
        GLRuu = GLR_ba(smatrix,0,0)
        GLRud = GLR_ba(smatrix,0,1)
        GLRdu = GLR_ba(smatrix,1,0)
        GLRdd = GLR_ba(smatrix,1,1)
        Gtot = GLRuu + GLRud + GLRdu + GLRdd 
        
        Tp = GLRuu + GLRdu - GLRud - GLRdd
        
        if G_total !=0:
            
            P = (GLRuu + GLRdu - GLRud - GLRdd)/(G_total)
        if G_total ==0:
            P=0
            
        dataP.append(P)
        Tp_list.append(Tp)
       
        

    
    plt.plot(energies, 100*np.array(dataP), color ='blue' )
    plt.xlabel("E [eV] ")
    plt.ylabel("% Polarization")
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    plt.tight_layout()
    plt.show()
    
    return energies, dataP,Tp_list


def plot_pol_LB(syst, emin, emax, de,check_hermiticity=True):
    
    dataP = []
    Tp_list = []
    energies = np.arange(emin,emax,de)

    # Function below calculates:  
    # transmission from left lead (1) with spin b,
    # too the right (0) lead with spin a
    def GLB_ba(smatrix,a,b):
        Gab = smatrix.transmission((0,a),(2,b))
        return Gab

    for energy in energies:
        

        smatrix = kwant.smatrix(syst, energy,check_hermiticity=check_hermiticity)

        # Total transmission:
        G_total = smatrix.transmission(0,2)

        ## Conductance with kwant function transmission
        GLRuu = GLB_ba(smatrix,0,0)
        GLRud = GLB_ba(smatrix,0,1)
        GLRdu = GLB_ba(smatrix,1,0)
        GLRdd = GLB_ba(smatrix,1,1)
        Gtot = GLRuu + GLRud + GLRdu + GLRdd 
        
        Tp = GLRuu + GLRdu - GLRud - GLRdd
        
        if G_total !=0:
            
            P = (GLRuu + GLRdu - GLRud - GLRdd)/(G_total)
        if G_total ==0:
            P=0
            
        dataP.append(P)
        Tp_list.append(Tp)
       
        

    
    plt.plot(energies, 100*np.array(dataP), color ='blue' )
    plt.xlabel("E [eV] ")
    plt.ylabel("% Polarization")
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    plt.tight_layout()
    plt.show()
    
    return energies, dataP,Tp_list


def calc_unnormalized_polarizations(syst, emin, emax, de,check_hermiticity=True):
    
    P_list_LR = []
    P_list_RL = []
    
    P_list_RR = []
    
    energies = np.arange(emin,emax,de)

    # Function below calculates:  
    # transmission from left lead (labeled 1) with spin b,
    # too the right lead (labeled 0) with spin a
    def GLR_ba(smatrix,a,b):
        Gab = smatrix.transmission((0,a),(1,b))
        return Gab
    
    # Function below calculates:  
    # transmission from right lead (labeled 0) with spin b,
    # too the left lead (labeled 1) with spin a
    def GRL_ba(smatrix,a,b):
        Gab = smatrix.transmission((1,a),(0,b))
        return Gab
    
    
    def GRR_ba(smatrix,a,b):
        Gab = smatrix.transmission((0,a),(0,b))
        return Gab

    for energy in energies:
        

        smatrix = kwant.smatrix(syst, energy,check_hermiticity=check_hermiticity)

        

        ## Conductance with kwant function transmission
        GLRuu = GLR_ba(smatrix,0,0)
        GLRud = GLR_ba(smatrix,0,1)
        GLRdu = GLR_ba(smatrix,1,0)
        GLRdd = GLR_ba(smatrix,1,1)
        
        
        P_LR = GLRuu + GLRdu - GLRud - GLRdd
        P_list_LR.append(P_LR)
        
        ## Conductance with kwant function transmission
        GRRuu = GRR_ba(smatrix,0,0)
        GRRud = GRR_ba(smatrix,0,1)
        GRRdu = GRR_ba(smatrix,1,0)
        GRRdd = GRR_ba(smatrix,1,1)
        
        
        P_RR = GRRuu + GRRdu  - GRRud - GRRdd 
        P_list_RR.append(P_RR)
        
        
#         ## Conductance with kwant function transmission
#         GRLuu = GRL_ba(smatrix,0,0)
#         GRLud = GRL_ba(smatrix,0,1)
#         GRLdu = GRL_ba(smatrix,1,0)
#         GRLdd = GRL_ba(smatrix,1,1)
        
        
#         P_RL = GRLuu + GRLdu - GRLud - GRLdd
#         P_list_RL.append(P_RL)
        
            
        
        
       
        

    
#     plt.plot(energies, 100*np.array(dataP), color ='blue' )
#     plt.xlabel("E [eV] ")
#     plt.ylabel("% Polarization")
#     plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
#     plt.tight_layout()
#     plt.show()
    
    return energies, P_list_LR,P_list_RR
    
    
def plot_updown_cond(syst, emin, emax, de):
    
    dataGup = []
    dataGdown = []
    
    
    energies = np.arange(emin,emax,de)
    
    def G_ab(smatrix,a,b):
        Sm = smatrix
        Gab = Sm.transmission((0,a),(1,b))
        return Gab

    
    
    for energy in energies:
        
      
            
        smatrix = kwant.smatrix(syst, energy)


        ## Conductance with kwant function transmission
        G00 = G_ab(smatrix,0,0)
        G01 = G_ab(smatrix,0,1)
        G10 = G_ab(smatrix,1,0)
        G11 = G_ab(smatrix,1,1)

        Gup = G00+G10
        Gdown = G11+G01

        # Add them to list
        dataGup.append(Gup)
        dataGdown.append(Gdown)
            

    plt.subplots(figsize=(10, 5))
    
    #Plot Gup, Gdown
    plt.plot(energies, dataGup, color ='blue',
              label = 'G$^\u2191$'
            )

    plt.plot(energies, dataGdown, color = 'red', label  = 'G$^\u2193$' )
    plt.xlabel("(E-E_F) [eV]")
    plt.ylabel("Transmission")
    
    
    plt.legend(loc='upper right')


    
    plt.tight_layout()
    plt.show()
    
    return energies,dataGup,dataGdown
    
    

def list_pol(syst, energies):
    
    dataP = []
    
    
    def G_lead_a2b(a,b):
        Gab = smatrix.transmission((0,a),(1,b))
        return Gab

    for energy in energies:

            smatrix = kwant.smatrix(syst, energy)
            
            # Total transmission:
            G_total = smatrix.transmission(0,1)

            ## Conductance with kwant function transmission
            G00 = G_lead_a2b(0,0)
            G01 = G_lead_a2b(0,1)
            G10 = G_lead_a2b(1,0)
            G11 = G_lead_a2b(1,1)
            Gtot = G00 + G01 + G10 + G11 
            
            if G_total !=0:
                P = (G00 - G01 + G10 - G11)/(G_total)
            if G_total ==0:
                P=0
           
            
        
            # Add them to list
            dataP.append(P)
            
    return dataP
# # Pure Polarization

# "pure polarization" $= T_{LR}^{\uparrow \uparrow} + T_{LR}^{\uparrow \downarrow} - T_{LR}^{\downarrow \uparrow} -T_{LR}^{\downarrow \downarrow}$

# In[3]:


def plot_pure_pol(syst, emin,emax, de):
    
    dataP = []
    energies = np.arange(emin,emax,de)

    
    
    def G_lead_a2b(a,b):
        Gab = smatrix.transmission((0,a),(1,b))
        return Gab

    for energy in energies:
            
            smatrix = kwant.smatrix(syst, energy)
            
            # Total transmission:
            G_total = smatrix.transmission(0,1)

            ## Conductance with kwant function transmission
            G00 = G_lead_a2b(0,0)
            G01 = G_lead_a2b(0,1)
            G10 = G_lead_a2b(1,0)
            G11 = G_lead_a2b(1,1)
            Gtot = G00 + G01 + G10 + G11 
            
            
            P = (G00 - G01 + G10- G11).real
            
           
            
        
            # Add them to list
            dataP.append(P)
        

    
    plt.plot(energies, dataP, color ='blue' )
    plt.xlabel("energy ")
    plt.ylabel("Pure Polarization (no transmission)")
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    plt.tight_layout()
    plt.show()

def list_cond_total(syst, emin, emax, de):
    
    dataG = []
    energies = np.arange(emin,emax,de)

    for energy in energies:
            
            smatrix = kwant.smatrix(syst, energy)

            ## Conductance with kwant function transmission
            G = smatrix.transmission(0,1)
            

            # Add them to list
            dataG.append(G)
            
    return energies,dataG



def plot_3_lists(list1,list2,list3,ymin=0,ymax=40,xmin=0,xmax=20):
    
    
    plt.clf()
    plt.plot(list1, list2 ,color = 'orange')
    plt.plot( list1, list3 ,color = 'blue')



    plt.xlabel("energy [eV]")
    plt.xlim(xmin,xmax)

    plt.ylim(ymin,ymax)
    plt.show()
        