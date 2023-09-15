# CISS_effect

IMPORTANT: The Kwant package comes with a file called "leads.py". 
In that file the function "modes" is defined. 
We changed the variable "stabilization" from the default  "None" to (False,False).  Alternatively we can change it to (True,False), (False, True) or (True,True), we found that it did not matter for our calculation. 
By changing "stabilization" in this way the internal algorithms Kwant uses are changed such that it accurately calculates the transmission. 
When one does not change ${ \it \text{stabilization}}$ this will result in a transmission that will diverge to unrealistic values for some energies. 
It is therefore crucial to change "stabilization".


In the folder "two terminal noninteracting"
We model a two-terminal junction consisting of gold,sulfur and helicene and calculate:
1. The spin-polarization for unmagnetized leads (UnmagnetizedSystem_Pz_SpinPolarization_git.ipynb)
2. The transmission and current for positive and negative magnetization B (MagnetizedSystem_Transmission_IVcurves.ipynb)


In the folder "two terminal noninteracting interacting", (GoldLeadMag_1probe_github.ipynb)
We add Büttiker voltage probes to the scattering region and calculate the current into left lead for positive and negative magnetization. 
From that we calculate i) the MR (the size of the effect) and ii) the odd and even part of the difference between the currents for opposite magnetizations of the lead.

In the folder "Modules" you can find the modules that are used in this calculation.

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
