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
	print("\n Arguments are       --- <FRB name string> <number of channels> <f_avg> <t_avg> <t_left/ms> <t_right/ms> <mode>\n")	
	print(" Supported Modes are --- peakwidth (Get peak location and estimate width)")
	print("                         profile   (Plot pulse profile)")
	print("                         spec      (Plot spectra)")
	
	print("\n            Now let's try again!\n")
	
	return(0)

#	--------------------------	Read inputs	-------------------------------
if(len(sys.argv)!=8):
	print_instructions()
	sys.exit()

frbname		=	sys.argv[1]										#	FRB name string (YYMMDDx)
nchan		=	int(sys.argv[2])								#	Number of fequency channels
ffac		=	int(sys.argv[3])								#	Averaging factor along frequency
tavgfac		=	int(sys.argv[4])								#	Averaging factor along time
tbasems		=	[float(sys.argv[5]),float(sys.argv[6])]			#	Time baseline in ms
exmode		=	sys.argv[7]										#	What to do

#	-------------------------	Do steps	-------------------------------

frbdm		=	-1.0
frblocms	=	-1
fmhz0		=	-1.0
tpeakms		=	0.0
	
frblist 	=	np.genfromtxt(frbcat)
for ifrb in frblist:
	if(str(int(ifrb[0]))==frbname):
		fmhz0		=	ifrb[1]
		frbdm		=	ifrb[3]
		frblocms	=	int(ifrb[4])	
		tpeakms		=	float(ifrb[14])		

print("FRB {} F0 = {:.2f} MHz DM = {} Peak = {} ms".format(frbname, fmhz0, frbdm, tpeakms))

#	[Step 0]	Get peak location & Calculate FRB width
if (exmode=='peakwidth'):
	print("\n[Step 0] Peak location at time resolution = {:.6f} ms".format(tavgfac*nchan*Raw_time_res_ms))
	(tpeak,eqwms,w95ms,lw95ms,rw95ms,wms,lwms,rwms) = cal_width(frbname, frbdm, nchan, ffac, tavgfac, max([abs(tbasems[0]), abs(tbasems[1])]))
	print("\n we = {:.3f} ms \n w = {:.3f} ms start = {:.3f} ms end = {:.3f} ms\n w95 = {:.3f} ms start95 = {:.3f} ms end95 = {:.3f} ms\n".format(eqwms,wms,lwms,rwms,w95ms,lw95ms,rw95ms))

#	[Step 1]	Plot sub-band profiles
if (exmode=='profile'):
	print("\n[Step 0] Plotting profiles at time resolution = {:.6f} ms".format(tavgfac*nchan*Raw_time_res_ms))
	tpeak		=	sub_ts(frbname, frbdm, nchan, ffac, tavgfac,tbasems, fmhz0, tpeakms)

#	[Step 2]	Plot sub-component spectra
if (exmode=='spec'):
	print("\n[Step 1] Plotting spectra at time resolution = {:.6f} ms".format(tavgfac*nchan*Raw_time_res_ms))
	tpeak		=	sub_spec(frbname, frbdm, nchan, ffac, tavgfac,fmhz0,tbasems, tpeakms)



















































