#
#	Functions for fitting pulse profiles
#
#								AB, 2 november 2023
#
#	Function list
#
#	expg(x, t, aa, x0, w):
#		Gaussian convolved with exponential
#
#	ncompexp(x, n, t, pars):
#		n Gaussians convolved with the same exponential	
#
#	r2ad(yarr, farr, p):
#		Returns the adjusted R^2
#
#	gsn(x, a0, x0, w):
#		Gaussian
#
#	--------------------------	Import modules	---------------------------

import os, sys
import numpy as np
from scipy import special

#	-------------------------------------------------------------------------

def expg(x, t, aa, x0, w):

#	Gaussian convolved with exponential
	
	exg		=	aa * np.exp((x0-x) / t) * special.erfc((x0 + (w*w/t) - x) / (np.sqrt(2.0)*w))
	
	return(exg)

#	-------------------------------------------------------------------------

def ncompexp(x, n, t, *pars):

#	n Gaussians convolved with the same exponential	
		
	nexg	=	0.0
	for i in range(0,n):
		nexg	=	nexg + expg(x, t, pars[0][i*3], pars[0][i*3+1], pars[0][i*3+2])
	
	return(nexg)

#	-------------------------------------------------------------------------

def r2ad(yarr, farr, p):

#	Returns the adjusted R^2
	
	ymn		=	np.nanmean(yarr)
	res		=	np.sum((yarr - farr)**2)
	tvar	=	np.sum((yarr - ymn)**2)
	nnn		=	len(yarr)
	
	r2		=	1.0 - res*(nnn-1)/(tvar*(nnn-p-1))
	
	return(r2)

#	-------------------------------------------------------------------------

def gsn(x, a0, x0, w):

#	Gaussian 
	
	nexg	= 	a0*np.exp(-((x - x0)/w)**2)
	
	return(nexg)

#	-------------------------------------------------------------------------



































