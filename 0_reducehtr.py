#
#	Script for FRB polarization analysis
#
#								AB, April 2023

#	--------------------------	Import modules	---------------------------

import os, sys
import numpy as np
from globalpars import *
from htrfns.anfns import *

def print_instructions():

	#	Print instructions to terminal
	
	print("\n            You probably need some assistance here!\n")
	print("\n Arguments are       --- <FRB name string> <number of channels> <f_avg> <t_avg> <mode> \n")
	print(" Supported Modes are --- go        (do everything)")
	print("                         locate    (locate the FRB in time)")
	print("                         extract   (Extract zoomed data)")
	print("                         volds     (generate voltage dynamic spectra)")
	print("                         stksds    (generate I,Q,U,V dynamic spectra)")
	print("                         favg      (Average I,Q,U,V dynamic spectra along frequency)")
	print("                         avgds     (Average I,Q,U,V dynamic spectra along time)")
	print("                         genoise   (Subtract mean & generate noise spectra)")
	print("                         badchan   (Generate list of bad/noisy channels)")
	print("                         its       (Generate I time series)")
	print("                         getstks   (extract to stksds -- see above)")
	print("                         getit     (avgds to its -- see above)")
	
	print("\n            Now let's try again!\n")
	
	return(0)

#	--------------------------	Read inputs	-------------------------------
if(len(sys.argv)!=6):
	print_instructions()
	sys.exit()

frbname		=	sys.argv[1]					#	FRB name string (YYMMDDx)
nchan		=	int(sys.argv[2])			#	Number of fequency channels
ffac		=	int(sys.argv[3])			#	Averaging factor along frequency
tavgfac		=	int(sys.argv[4])			#	Averaging factor along time
exmode		=	sys.argv[5]					#	What to do	-----	go => do all steps

#	-------------------------	Do steps	-------------------------------

frbdm		=	-1.0
frblocms	=	-1
fmhz0		=	-1.0
	
frblist 	=	np.genfromtxt(frbcat)
for ifrb in frblist:
	if(str(int(ifrb[0]))==frbname):
		fmhz0		=	ifrb[1]
		frbdm		=	ifrb[3]
		frblocms	=	int(ifrb[4])		

print("FRB {} F0 = {:.2f} MHz DM = {} location = {} ms".format(frbname, fmhz0, frbdm, frblocms))

#	[Step 0]	Locate FRB in time series
if ((exmode=='locate') or (exmode=='go')):
	print("\n[Step 0] Locating FRB in time...\n")
	(frblocms,frbsnrms) =	locate_frb(frbname, frbdm)
	print("Located at {} ms at S/N = {:.2f}".format(frblocms,frbsnrms))

#	[Step 1]	Extract zoomed data
if ((exmode=='extract') or (exmode=='go') or (exmode=='getstks')):
	print("\n[Step 1] Extracting zoomed data...\n")
	zlen	=	extract_zoom(frbname, frbdm, frblocms, ZoomLenms/2.0)
	print("Zoomed votage data (len = {}) written to disk".format(zlen))

#	[Step 2]	Construct X and Y voltage dynamic spectra
if ((exmode=='volds') or (exmode=='go') or (exmode=='getstks')):
	print("\n[Step 2] Constructing voltage dynamic spectra for X and Y with {} channels...\n".format(nchan))
	dslen	=	voltage_ds(frbname,frbdm, nchan)
	print("Dynamic spectra (len = {}) written to disk".format(dslen))

#	[Step 3]	Construct dynamic spectra for I,Q,U,V
if ((exmode=='stksds') or (exmode=='go') or (exmode=='getstks')):
	print("\n[Step 3] Constructing dynamic spectra for I, Q, U, V with {} channels...\n".format(nchan))
	natresms	=	gen_stokes_ds(frbname, frbdm, nchan, fmhz0)
	print("Dynamic spectra for Stokes parameters written to disk with time resolution = {:.6f} ms".format(natresms))

#	[Step 4]	Average dynamic spectra for I,Q,U,V in frequency
if ((exmode=='favg') or (exmode=='go')):
	if(nchan%ffac == 0):
		print("\n[Step 4] Average dynamic spectra for I,Q,U,V in frequency by {}...\n".format(ffac))
		if(ffac>1):
			avgchans	=	favg_stokes(frbname, frbdm, nchan, ffac)
			print("Frequency averaged Stokes parameters with {} chanels and frequency resolution = {:.6f} MHz".format(avgchans,NatBWmhz/(nchan*ffac)))
		else:
			print("No frequency averaging required...")
	else:
		print("Wrong frequency averaging factor! Should be a factor of {}".format(nchan))

#	[Step 5]	Average dynamic spectra for I,Q,U,V in time
if ((exmode=='avgds') or (exmode=='go') or (exmode=='getit')):
	print("\n[Step 5] Average dynamic spectra for I,Q,U,V in time with {} channels...\n".format(nchan))
	if(tavgfac>1):
		avglen	=	avg_stokes_ds(frbname, frbdm, nchan, ffac, tavgfac)
		print("Time averaged Stokes parameters (len = {}) with time resolution = {:.6f} ms".format(avglen,tavgfac*nchan*Raw_time_res_ms))
	else:
		print("No time averaging required...")

#	[Step 6]	Generate noise spectra
if ((exmode=='genoise') or (exmode=='go') or (exmode=='getit')):
	print("\n[Step 6] Generating noise spectra for all Stokes parameters")
	gen_noise_spec(frbname, frbdm, nchan, ffac, tavgfac)

#	[Step 7]	Generate Bad channel list
if ((exmode=='badchan') or (exmode=='go') or (exmode=='getit')):
	print("\n[Step 7] Generating Bad channel list at time resolution = {:.6f} ms".format(tavgfac*nchan*Raw_time_res_ms))
	gen_bad_chan_list(frbname, frbdm, nchan, ffac, tavgfac)

#	[Step 8]	Generate I time series
if ((exmode=='its') or (exmode=='go') or (exmode=='getit')):
	print("\n[Step 8] Generating Stokes I time series with time resolution = {:.6f} ms".format(tavgfac*nchan*Raw_time_res_ms))
	gen_I_ts(frbname, frbdm, nchan, ffac, tavgfac)




















































