#
#	Script for detailed polarization inspection across frequency
#
#								AB, June 2023

#	--------------------------	Import modules	---------------------------

import os, sys
import numpy as np
from globalpars import *
from htrfns.polspecfns import *

def print_instructions():

	#	Print instructions to terminal
	
	print("\n            You probably need some assistance here!\n")
	print("\n Arguments are     --- <FRB name string> <number of channels> <f_avg> <t_avg> <RM> <t_start/ms> <t_end/ms> <mode> <refang>\n")	
	print(" Supported Modes are --- polwave    (Polarization spectra vs wavelength)")
	print("                     --- anglam     (Polarization angles vs wavelength)")
	print("                     --- stkspec    (Stokes parameters vs frequency)")
	
	print("\n            Now let's try again!\n")
	
	return(0)

#	--------------------------	Read inputs	-------------------------------
if(len(sys.argv)<10):
	print_instructions()
	sys.exit()

frbname		=	sys.argv[1]					#	FRB name string (YYMMDDx)
nchan		=	int(sys.argv[2])			#	Number of fequency channels
ffac		=	int(sys.argv[3])			#	Averaging factor along frequency
tavgfac		=	int(sys.argv[4])			#	Averaging factor along time
frbrm0		=	float(sys.argv[5])			#	RM used for de-rotation
tbasems		=	[float(sys.argv[6]),float(sys.argv[7])]			#	Time baseline in ms
exmode		=	sys.argv[8]					#	What to do
refdeg		=	float(sys.argv[9])			#	RM used for de-rotation

#	-------------------------	Do steps	-------------------------------

frbdm		=	-1.0
frblocms	=	-1
fmhz0		=	-1.0
	
frblist 	=	np.genfromtxt(frbcat)
for ifrb in frblist:
	if(str(int(ifrb[0]))==frbname):
		fmhz0		=	ifrb[1]	
		frbdm		=	str(ifrb[3])			

print("FRB {} F0 = {:.2f} MHz DM = {} RM = {:.2f} time resolution = {:.6f} ms".format(frbname, fmhz0, frbdm, frbrm0,tavgfac*nchan*Raw_time_res_ms))

#	[Step 0]	Plot Polarization spectra vs wavelength
if (exmode=='polwave'):
	print("\n[Step 0] Plotting spectra after undoing Faraday rotation for RM = {} rad/m2".format(frbrm0))
	polambda(frbname, frbdm, nchan, ffac, tavgfac, fmhz0, frbrm0, tbasems)

#	[Step 0]	Plot angles vs wavelength
if (exmode=='anglam'):
	print("\n[Step 0] Plotting angles after undoing Faraday rotation for RM = {} rad/m2".format(frbrm0))
	gamdelam(frbname, frbdm, nchan, ffac, tavgfac, fmhz0, frbrm0, tbasems, refdeg)

#	[Step 0]	Plot Stokes parameters vs frequency
if (exmode=='stkspec'):
	print("\n[Step 0] Plotting spectra after undoing Faraday rotation for RM = {} rad/m2".format(frbrm0))
	stokespec(frbname, frbdm, nchan, ffac, tavgfac, fmhz0, frbrm0, tbasems)






































