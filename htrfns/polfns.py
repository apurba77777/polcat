#
#	Functions for FRB polarization analysis
#
#								AB, April 2023

#	Import modules

import os
import numpy as np
from RMtools_1D.do_RMsynth_1D import run_rmsynth
from RMtools_1D.do_RMclean_1D import run_rmclean
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
import matplotlib.pyplot as plt

#	---------------------------------------------------------------------------------

def pangdiffdeg(ang, ang0):
	ang			=	np.deg2rad(ang)
	ang0		=	np.deg2rad(ang0)
	dpang		=	np.rad2deg(np.arcsin(np.sin(ang-ang0)))
	return(dpang)

#	---------------------------------------------------------------------------------

def fr_rmtool(fghz, iquv, diquv):
	'''
	Function to determine RM using RM synthesis with RMtool
	
	Inputs	-	Frequencies in GHz
				I Q U spectrum
				I Q U noise spectrum					
	'''
	
	rmtdata	=	np.array([fghz*1.0e9, iquv[0], iquv[1], iquv[2], diquv[0], diquv[1], diquv[2]])
	
	rmd,rmad=	run_rmsynth(rmtdata, polyOrd=3, phiMax_radm2=1.0e3, dPhi_radm2=1.0, nSamples=100.0, weightType='variance', fitRMSF=False, noStokesI=False, phiNoise_radm2=1000000.0, \
						nBits=32, showPlots=True, debug=False, verbose=False, log=print, units='Jy/beam', prefixOut='prefixOut', saveFigures=None,fit_function='log')
	
	rmc		=	run_rmclean(rmd, rmad, 0.1, maxIter=1000, gain=0.1, nBits=32, showPlots=False, verbose=False, log=print)
	
	print(rmc[0])
	
	res		=	[rmc[0]['phiPeakPIfit_rm2'], rmc[0]['dPhiPeakPIfit_rm2'], rmc[0]['polAngle0Fit_deg'], rmc[0]['dPolAngle0Fit_deg']]
	
	return(res)	

#	---------------------------------------------------------------------------------




























































