#
#	Investigate polarization properties of FRBs
#
#								AB, October 2023

#	--------------------------	Import modules	---------------------------

import os, sys
import numpy as np
from globalpars import *
from polfns.polprofile import *
from polfns.polreduce import *

def print_instructions():

	#	Print instructions to terminal
	
	print("\n            You probably need some assistance here!\n")
	print("\n Arguments are       --- <FRB name string> <number of channels> <f_avg> <t_avg> <RM> <t_left/ms> <t_right/ms> <sub_band> <refangdeg> <mode> <b_avg>\n")	
	print(" Supported Modes are --- pkl        (Process and pickle the outputs)")
	print("                         iquvt	   (Plot IQUV profiles and dynspec)")
	print("                         stkst	   (Plot IQUV profiles and dynspec separately)")
	print("                         stkspec	   (Plot fractional IQUV spectra)")
	print("                         polspec	   (Plot polarization fraction spectra)")
	print("                         fitcomp	   (Fit n gaussians)")
	print("                         ilvt	   (Plot ILV profiles)")
	print("                         ilvpat	   (Plot ILV & PA profiles)")
	print("                         ilvparm	   (Plot ILV, PA & RM profiles)")
	print("                         irm		   (Plot RM profile)")
	print("                         dpa		   (Plot PA profile and its rate of change)")
	print("                         fracparm   (Plot polarization fraction, PA & RM profiles)")
	print("                         fracpa	   (Plot polarization fraction & PA profiles)")
	print("                         haspec	   (Plot Poincare sphere frequency track in Hammer projection)")
	
	print("\n            Now let's try again!\n")
	
	return(0)

#	--------------------------	Read inputs	-------------------------------
if(len(sys.argv)<12):
	print_instructions()
	sys.exit()

frbname		=	sys.argv[1]										#	FRB name string (YYMMDDx)
nchan		=	int(sys.argv[2])								#	Number of fequency channels
ffac		=	int(sys.argv[3])								#	Averaging factor along frequency
tavgfac		=	int(sys.argv[4])								#	Averaging factor along time
frbrm0		=	float(sys.argv[5])								#	RM used for derotation
twinms		=	[float(sys.argv[6]),float(sys.argv[7])]			#	Time window in ms
subid		=	int(sys.argv[8])								#	Sub-band id
refdeg		=	float(sys.argv[9])								#	Reference angle in deg (for PSP derotation)
exmode		=	sys.argv[10]									#	What to do
bavgfac		=	int(sys.argv[11])								#	Averaging factor along time

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
		tpeakms		=	float(ifrb[14])		

print("FRB {} F0 = {:.2f} MHz DM = {} location = {} ms".format(frbname, fmhz0, frbdm, frblocms))

#	Process and pickle the outputs
if (exmode=='pkl'):
	print("\nPickling at time resolution = {:.6f} ms".format(tavgfac*nchan*Raw_time_res_ms))
	picklit(frbname, frbdm, nchan, ffac, tavgfac, fmhz0, frbrm0, twinms, tpeakms, subid)
else:
	print("\nPlotting at time resolution = {:.6f} ms".format(tavgfac*nchan*Raw_time_res_ms))
	plothemall(frbname, frbdm, nchan, ffac, tavgfac, fmhz0, frbrm0, twinms, tpeakms, subid, exmode, bavgfac)



































