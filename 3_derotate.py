#
#	Script for FRB polarization analysis
#
#								AB, April 2023

#	--------------------------	Import modules	---------------------------

import os, sys
import numpy as np
from globalpars import *
#from htrfns.dynamic_fns import *
#from htrfns.scintillation_fns import *
from htrfns.burstfns import *

def print_instructions():

	#	Print instructions to terminal
	
	print("\n            You probably need some assistance here!\n")
	print("\n Arguments are       --- <FRB name string> <number of channels> <f_avg> <t_avg> <RM>\n")	
	print("\n            Now let's try again!\n")
	
	return(0)

#	--------------------------	Read inputs	-------------------------------
if(len(sys.argv)<6):
	print_instructions()
	sys.exit()

frbname		=	sys.argv[1]					#	FRB name string (YYMMDDx)
nchan		=	int(sys.argv[2])			#	Number of fequency channels
ffac		=	int(sys.argv[3])			#	Averaging factor along frequency
tavgfac		=	int(sys.argv[4])			#	Averaging factor along time
frbrm0		=	float(sys.argv[5])			#	RM used for de-rotation

#	-------------------------	Do steps	-------------------------------

frbdm		=	-1.0
frblocms	=	-1
fmhz0		=	-1.0
	
frblist 	=	np.genfromtxt(frbcat)
for ifrb in frblist:
	if(str(int(ifrb[0]))==frbname):
		fmhz0		=	ifrb[1]	
		frbdm		=	str(ifrb[3])		

print("FRB {} F0 = {:.2f} MHz DM = {} RM = {:.2f}".format(frbname, fmhz0, frbdm, frbrm0))
	
#	[Step 0]	Undo Faraday rotation
print("\n[Step 0] Undoing Faraday rotation for RM = {} rad/m2".format(frbrm0))
unfarot(frbname, frbdm, nchan, ffac, tavgfac, fmhz0, frbrm0)
















































