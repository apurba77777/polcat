#
#	Script for FRB polarization analysis
#
#								AB, April 2023

#	--------------------------	Import modules	---------------------------

import os, sys
import numpy as np
from globalpars import *
from htrfns.dynamic_fns import *
from htrfns.scintillation_fns import *
from htrfns.burstfns import *

def print_instructions():

	#	Print instructions to terminal
	
	print("\n            You probably need some assistance here!\n")
	print("\n Arguments are       --- <FRB name string> <number of channels> <f_avg> <t_avg> <startms> <stopms> <startchan> <endchan>\n")	

	print("\n            Now let's try again!\n")
	
	return(0)

#	--------------------------	Read inputs	-------------------------------
if(len(sys.argv)<9):
	print_instructions()
	sys.exit()

frbname		=	sys.argv[1]					#	FRB name string (YYMMDDx)
nchan		=	int(sys.argv[2])			#	Number of fequency channels
ffac		=	int(sys.argv[3])			#	Averaging factor along frequency
tavgfac		=	int(sys.argv[4])			#	Averaging factor along time
startms		=	float(sys.argv[5])			#	Starting time from peak (ms)
stopms		=	float(sys.argv[6])			#	Stopping time from peak (ms)
startchan	=	int(sys.argv[7])
endchan		=	int(sys.argv[8])

#	-------------------------	Do steps	-------------------------------

frbdm		=	-1.0
frblocms	=	-1
fmhz0		=	-1.0
frbrm0		=	0.0
frbphi0		=	0.0
tpeakms		=	0.0

if(startchan < 0):
	startchan	=	0

if(endchan <= 0):
	endchan	=	nchan-1
	
frblist 	=	np.genfromtxt(frbcat)
for ifrb in frblist:
	if(str(int(ifrb[0]))==frbname):
		fmhz0		=	ifrb[1]
		frbdm		=	str(ifrb[3])		
		tpeakms		=	float(ifrb[14])			

print("FRB {} F0 = {:.2f} MHz DM = {} RM = {:.2f} Phi_0 = {:.2f} deg".format(frbname, fmhz0, frbdm, frbrm0, frbphi0))

print("\n[Step 1] Estimating RM")
(res_rmsynth, res_rmnest, res_rmtool, tpeak) = estimate_rm(frbname, frbdm, nchan, ffac, tavgfac, fmhz0, startms, stopms, 1.0e3, 1.0, startchan, endchan, tpeakms)
frbrm0		=	res_rmtool[0]
frbphi0		=	res_rmtool[2]	
#plot_polang(frbname, frbdm, nchan, ffac, tavgfac, fmhz0, tpeak, startms, stopms, frbphi0, frbrm0)

	
	















































