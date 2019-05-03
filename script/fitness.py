import pyvolve
import os
import ete3
import numpy
import AuxiliarFunctions as af
import numpy as np

	
fichier = open("codon_fitness.txt", "w") 


codon_fitness = np.random.normal(size = 61) 
for i in range(61):
	fichier.write("%s\n" % (str(codon_fitness[i])))    

fichier.close()
