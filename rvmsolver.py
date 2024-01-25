#
#	Script for Solving RVM parameters for two component FRBs
#
#								AB, October 2023

#	--------------------------	Import modules	---------------------------

import os, sys
import numpy as np

def print_instructions():

	#	Print instructions to terminal
	
	print("\n            Hi Stranger! You probably need some assistance here!\n")
	print("\n Arguments are       --- <Your name> <ta> <tb> <tres> <ra> <dra> <rb> <drb> <width>\n")	
	
	print("\n            Now let's try again!\n")
	
	return(0)

#	--------------------------	Read inputs	-------------------------------
if(len(sys.argv)<10):
	print_instructions()
	sys.exit()

yname	=	sys.argv[1]				#	Name
ta		=	float(sys.argv[2])		#	Fastest swing point in A in ms
tb		=	float(sys.argv[3])		#	Fastest swing point in B in ms
tres	=	float(sys.argv[4])		#	Time resolution in ms
ra		=	float(sys.argv[5])		#	Fastest swing rate in A in deg/ms
dra		=	float(sys.argv[6])		#	error in Fastest swing rate in A in deg/ms
rb		=	float(sys.argv[7])		#	Fastest swing rate in B in deg/ms
drb		=	float(sys.argv[8])		#	error in Fastest swing rate in B in deg/ms
pw		=	float(sys.argv[9])		#	Width in ms

print("\n                  Hi "+yname+"!\n")

dt		=	np.abs(tb-ta)
aa		=	(180.0/dt)/ra
ab		=	(180.0/dt)/rb

delt	=	2*tres
dela	=	aa*np.sqrt((delt/dt)**2 + (dra/ra)**2)
delb	=	ab*np.sqrt((delt/dt)**2 + (drb/rb)**2)

delalp	=	0.0

if(np.abs(aa+ab)>1.0):
	alp = 	90.0
	print(aa,ab)
	print("Caution! Assuming inclination = 90 deg")
else:
	alp		=	np.rad2deg(np.arccos(-(aa+ab)/2.0))
	delalp	=	np.rad2deg(np.sqrt(dela**2 + delb**2)/(2*np.sin(np.deg2rad(alp))))

theta	=	np.rad2deg(np.arctan(np.sqrt((4.0 - (aa+ab)**2)/(aa-ab)**2)))
thetb	=	180.0 - theta
deltheta=	np.rad2deg(np.sqrt((dela**2 + delb**2)*(4*(aa*aa - ab*ab)**2 + 4 - (aa+ab)**2)) / (4*(1.0 - aa*ab)))

ba		=	alp - theta
bb		=	alp - thetb
delb	=	np.sqrt(delalp**2 + deltheta**2)

alrad	=	np.deg2rad(alp)
thetrad	=	np.deg2rad(thetb)

rho		=	np.rad2deg(np.arccos( np.cos(alrad)*np.cos(thetrad)+ np.sin(alrad)*np.sin(thetrad)*np.cos(pw/(2*dt))))

delrho	=	np.rad2deg(np.sqrt((np.cos(alrad)*np.sin(thetrad)*np.deg2rad(deltheta))**2 + (np.sin(alrad)*np.cos(thetrad)*np.deg2rad(delalp))**2 + \
				(np.cos(alrad)*np.sin(thetrad)*np.cos(pw/(2*dt))*np.deg2rad(delalp))**2 + (np.sin(alrad)*np.cos(thetrad)*np.cos(pw/(2*dt))*np.deg2rad(deltheta))**2 + \
					(np.sin(alrad)*np.sin(thetrad)*np.sin(pw/(2*dt))*tres*np.sqrt(dt*dt + pw*pw)/(2*dt*dt))**2)/np.sin(np.deg2rad(rho)))

print("RVM parameters are here ...\n")
print("Inclination            = %.2f deg +/- %.2f"%(alp,delalp))
print("Magnetic obliquity     = %.2f deg / %.2f deg +/- %.2f"%(theta,thetb,deltheta))
print("Impact angles          = %.2f deg / %.2f deg +/- %.2f"%(ba,bb,delb))
print("beam width             = %.2f deg +/- %.2f"%(rho,delrho))





























































