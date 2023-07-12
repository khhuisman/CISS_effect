# CISS_effect

The Kwant package comes with a file called "leads.py". 
In that file the function "modes" is defined. 
We changed the variable "stabilization" from the default  "None: to (False,False).  Alternatively we can change it to (True,False), (False, True) or (True,True), we found that it did not matter for our calculation. 
By changing "stabilization" in this way the internal algorithms Kwant uses are changed such that it accurately calculates the transmission. 
When one does not change \textit{stabilization} this will result in a transmission that will diverge to unrealistic values for some energies. 
It is therefore crucial to change "stabilization".


In the folder 2T_noninteracting
We model a two-terminal junction consisting of gold,sulfur and helicene and calculate:
1. The spin-polarization for unmagnetized leads (UnmagnetizedSystem_Pz_SpinPolarization_git.ipynb)
2. The transmission and current for positive and negative magnetization B (MagnetizedSystem_Transmission_IVcurves.ipynb)


In the folder 2T_interacting, (GoldLeadMag_1probe_github.ipynb)
We add Büttiker voltage probes to the scattering region and calculate the current into left lead for positive and negative magnetization. 
From that we calculate i) the MR (the size of the effect) and ii) the odd and even part of the difference between the currents for opposite magnetizations of the lead.

In the folder Modules you can find the modules that are used in this calculation.
