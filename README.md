# CISS_effect

The Kwant package comes with a file called "lead.py". 
In that file the function "modes" is defined. 
We changed the variable "stabilization" from the default  "None: to (False,False).  Alternatively we can change it to (True,False), (False, True) or (True,True), we found that it did not matter for our calculation. 
By changing "stabilization" in this way the internal algorithms Kwant uses are changed such that it accurately calculates the transmission. 
When one does not change \textit{stabilization} this will result in a transmission that will diverge to unrealistic values for some energies. 
It is therefore crucial to change "stabilization".
