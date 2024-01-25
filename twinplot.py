#
#	Script for plotting twin FRBs
#
#								AB, june 2023

#	--------------------------	Import modules	---------------------------

import os, sys
import numpy as np
from globalpars import *
from htrfns.twinfns import *

def print_instructions():

	#	Print instructions to terminal
	
	print("\n            You probably need some assistance here!\n")
	print("\n Arguments are -- <FRB name string> <number of channels> <f_avg> <t_avg> <t_left/ms> <t_right/ms> <FRB name string> <number of channels> <f_avg> <t_avg> <t_left/ms> <t_right/ms>\n")	
	
	print("\n            Now let's try again!\n")
	
	return(0)

#	--------------------------	Read inputs	-------------------------------
if(len(sys.argv)!=13):
	print_instructions()
	sys.exit()

frbname		=	sys.argv[1]										#	FRB name string (YYMMDDx)
nchan		=	int(sys.argv[2])								#	Number of fequency channels
ffac		=	int(sys.argv[3])								#	Averaging factor along frequency
tavgfac		=	int(sys.argv[4])								#	Averaging factor along time
tbasems		=	[float(sys.argv[5]),float(sys.argv[6])]			#	Time baseline in ms

frbname0	=	sys.argv[7]										#	FRB name string (YYMMDDx)
nchan0		=	int(sys.argv[8])								#	Number of fequency channels
ffac0		=	int(sys.argv[9])								#	Averaging factor along frequency
tavgfac0	=	int(sys.argv[10])								#	Averaging factor along time
tbasems0	=	[float(sys.argv[11]),float(sys.argv[12])]		#	Time baseline in ms

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

frbdm0		=	-1.0
frblocms0	=	-1
fmhz00		=	-1.0
tpeakms0	=	0.0
	
for ifrb in frblist:
	if(str(int(ifrb[0]))==frbname0):
		fmhz00		=	ifrb[1]
		frbdm0		=	ifrb[3]
		frblocms0	=	int(ifrb[4])	
		tpeakms0	=	float(ifrb[14])	

print("FRB {} F0 = {:.2f} MHz DM = {} Peak = {} ms".format(frbname, fmhz0, frbdm, tpeakms))
print("FRB {} F0 = {:.2f} MHz DM = {} Peak = {} ms".format(frbname0, fmhz00, frbdm0, tpeakms0))

#tpeak		=	twin_side_plot(frbname, frbdm, nchan, ffac, tavgfac,tbasems, fmhz0, tpeakms, frbname0, frbdm0, nchan0, ffac0, tavgfac0, tbasems0, fmhz00, tpeakms0)

tpeak		=	twin_same_plot(frbname, frbdm, nchan, ffac, tavgfac, 1.27, tpeakms, frbname0, frbdm0, nchan0, ffac0, tavgfac0, 0.81, tpeakms0)
#tpeak		=	twin_same_plot(frbname, frbdm, nchan, ffac, tavgfac, 1.00, tpeakms, frbname0, frbdm0, nchan0, ffac0, tavgfac0, 1.00, tpeakms0)
















































