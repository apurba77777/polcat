#
#	Functions describing models for polarization
#
#								AB, May 2023

#	Function list
#
#		chi_nonrel_rvm_ew01(tms, t0ms, tperms, chi0, theta0, inclin0):
#				PA for non-relativisitc RVM using EW 2001 model
#
#		ddt_chi_nonrel_rvm_ew01(tms, t0ms, tperms, theta0, inclin0):
#				Time derivative of PA for non-relativisitc RVM using EW 2001 model
#
#		chi_rvm_tong(tms, t0ms, tperms, chi0, theta0, inclin0, lam):
#				PA for RVM using Tong 2021 model
#
#		chi_rrvm_cr12(tms, t0ms, tperms, chi0, theta0, inclin0, rr):
#				PA for relativisitc RVM using CR 2012 model
#
#		chi_rrvm_ha01(tms, t0ms, tperms, chi0, theta0, inclin0, rr):
#				PA for relativisitc RVM using HA 2001 model
#
#		chi_srrvm_p20(tms, t0ms, tperms, chi0, theta0, inclin0, beta):
#				PA for relativisitc RVM using Poutaten 2020 model
#
#		chi_rvm_fr_l22(tms, t0ms, tperms, chi0, theta0, inclin0, fac):
#				PA for RVM using L 2022 model
#
#	--------------------------	Import modules	---------------------------

import os, sys
import numpy as np
import matplotlib.pyplot as plt
from globalpars import *
from scipy.optimize import curve_fit

#	-------------------------------------------------------------------------------

def chi_nonrel_rvm_ew01(tms, t0ms, tperms, chi0, theta0, inclin0):

	#	PA for non-relativisitc RVM using EW 2001 model
	
	theta		=	np.deg2rad(theta0)
	inclin		=	np.deg2rad(inclin0)
	phi			=	(tms - t0ms)*2*np.pi/tperms 
	schi		=	np.sin(theta)*np.sin(phi)
	cchi		=	np.sin(inclin)*np.cos(theta) - np.cos(inclin)*np.sin(theta)*np.cos(phi) 
	chideg		=	chi0 + np.rad2deg(np.arctan(schi/cchi))
	return (chideg)
	
#	-------------------------------------------------------------------------------

def ddt_chi_nonrel_rvm_ew01(tms, t0ms, tperms, theta0, inclin0):

	#	Time derivative of PA for non-relativisitc RVM using EW 2001 model
	
	theta		=	np.deg2rad(theta0)
	inclin		=	np.deg2rad(inclin0)
	phi			=	(tms - t0ms)*2*np.pi/tperms 
	num			=	np.sin(theta)*(np.sin(theta)*np.cos(inclin) - np.sin(inclin)*np.cos(theta)*np.cos(phi))
	den			=	(-np.sin(inclin)*np.cos(theta) + np.cos(inclin)*np.sin(theta)*np.cos(phi))**2 + (np.sin(theta)*np.sin(phi))**2
	sens		=	1.0 
	dcdt		=	np.rad2deg(2*np.pi*(sens*num/den)/tperms)
	return (dcdt)
	
#	-------------------------------------------------------------------------------

def chi_nonrel_rvm_ew01_wob(tms, t0ms, tperms, chi0, theta0, inclin0, dinc, nfac):

	#	PA for non-relativisitc RVM using EW 2001 model
	
	phi			=	(tms - t0ms)*2*np.pi/tperms 
	theta		=	np.deg2rad(theta0)
	inclin		=	np.deg2rad(inclin0) + dinc*np.sin(nfac*phi)
	schi		=	np.sin(theta)*np.sin(phi)
	cchi		=	 -np.sin(inclin)*np.cos(theta) + np.cos(inclin)*np.sin(theta)*np.cos(phi) 
	chideg		=	chi0 + np.rad2deg(np.arctan(schi/cchi))
	return (chideg)
	
#	-------------------------------------------------------------------------------

def chi_rvm_tong(tms, t0ms, tperms, chi0, theta0, inclin0, lam):

	#	PA for RVM using Tong 2021 model
	
	theta		=	np.deg2rad(theta0)
	inclin		=	np.deg2rad(inclin0)
	phi			=	(tms - t0ms)*2*np.pi/tperms 
	schi		=	np.sin(theta)*np.sin(phi)
	cchi		=	 -np.sin(inclin)*np.cos(theta) + np.cos(inclin)*np.sin(theta)*np.cos(phi) 
	chideg		=	chi0 + np.rad2deg(np.arctan(schi/cchi))
	twst		=	np.rad2deg(lam*(1.0 - (np.cos(theta)*np.cos(inclin) + np.sin(theta)*np.sin(inclin)*np.cos(phi))**2))
	return (chideg+twst)
	
#	-------------------------------------------------------------------------------

def chi_rrvm_cr12(tms, t0ms, tperms, chi0, theta0, inclin0, rr):

	#	PA for relativisitc RVM using CR 2012 model
	
	theta		=	np.deg2rad(theta0)
	inclin		=	np.deg2rad(inclin0)
	phi			=	(tms - t0ms)*2*np.pi/tperms 
	schi		=	np.sin(theta)*np.sin(phi) - 3*rr*np.sin(inclin)
	cchi		=	 -np.sin(inclin)*np.cos(theta) + np.cos(inclin)*np.sin(theta)*np.cos(phi) 
	chideg		=	chi0 + np.rad2deg(np.arctan(schi/cchi))
	return (chideg)
	
#	-------------------------------------------------------------------------------

def chi_rrvm_ha01(tms, t0ms, tperms, chi0, theta0, inclin0, rr):

	#	PA for relativisitc RVM using HA 2001 model
	
	theta		=	np.deg2rad(theta0)
	inclin		=	np.deg2rad(inclin0)
	phi			=	(tms - t0ms)*2*np.pi/tperms 
	schi		=	np.sin(theta)*np.sin(phi)
	cchi		=	 -np.sin(inclin)*np.cos(theta) + np.cos(inclin)*np.sin(theta)*np.cos(phi) 
	chideg		=	chi0 + np.rad2deg(np.arctan(schi/cchi))
	dchi		=	np.rad2deg((23.0/144.0)*rr*np.sin(theta)*np.sin(inclin)*(np.sin(phi)**2))
	return (chideg-dchi)
	
#	-------------------------------------------------------------------------------

def chi_srrvm_p20(tms, t0ms, tperms, chi0, theta0, inclin0, beta):

	#	PA for relativisitc RVM using Poutaten 2020 model
	
	theta		=	np.deg2rad(theta0)
	inclin		=	np.deg2rad(inclin0)
	phi			=	(tms - t0ms)*2*np.pi/tperms 
	schi		=	np.sin(theta)*np.sin(phi) + beta*(np.sin(inclin)*np.sin(theta) + np.cos(inclin)*np.cos(theta)*np.cos(phi))
	cchi		=	 -np.sin(inclin)*np.cos(theta) + np.cos(inclin)*np.sin(theta)*np.cos(phi) - beta*np.cos(theta)*np.sin(phi)
	chideg		=	chi0 + np.rad2deg(np.arctan(schi/cchi))
	return (chideg)
	
#	-------------------------------------------------------------------------------

def chi_rvm_fr_l22(tms, t0ms, tperms, chi0, theta0, inclin0, fac):

	#	PA for RVM using L 2022 model
	
	theta		=	np.deg2rad(theta0)
	inclin		=	np.deg2rad(inclin0)
	phi			=	(tms - t0ms)*2*np.pi/tperms 
	schi		=	np.sin(theta)*np.sin(phi)
	cchi		=	 -np.sin(inclin)*np.cos(theta) + np.cos(inclin)*np.sin(theta)*np.cos(phi) 
	tanchi		=	schi/cchi
	chideg		=	chi0 + np.rad2deg(0.5*np.arctan(2*fac*tanchi,(1.0 - tanchi**2)))
	return (chideg)
	
#	-------------------------------------------------------------------------------

def fit_chi_nonrel_rvm_ew01_wob(tmsarr, chiarr, chierr):

	#	Fit PA vs time for non-relativisitc RVM
	
	goodchi			=	np.where(np.isfinite(chiarr))
	tmsarr0			=	tmsarr[goodchi]
	chiarr0			=	chiarr[goodchi]
	chierr0			=	chierr[goodchi]
	
	fitopt, fitcov	=	curve_fit(chi_nonrel_rvm_ew01_wob, tmsarr0, chiarr0, sigma=chierr0)
	fiterr			=	np.sqrt(np.diag(fitcov))
	return (fitopt, fiterr)
	
#	--------------------------	Analysis functions	-------------------------------

def chilamfn(lm, grm, chi0deg, pindx):
	
	#pindx		=	2.0
	lm0			=	0.25
	chideg		=	chi0deg + np.rad2deg(grm*(lm**pindx - lm0**pindx))
	#chideg		=	np.rad2deg(np.arctan(np.deg2rad(chideg)))
	
	return(chideg)

#	-------------------------------------------------------------------------------

def chirmfit(lm2, rm, chi0deg):
	chideg		=	chi0deg + np.rad2deg(rm*lm2)
	return(chideg)

#	-------------------------------------------------------------------------------

def fit_chi_nonrel_rvm(tmsarr, chiarr, chierr):

	#	Fit PA vs time for non-relativisitc RVM
	
	goodchi			=	np.where(np.isfinite(chiarr))
	tmsarr0			=	tmsarr[goodchi]
	chiarr0			=	chiarr[goodchi]
	chierr0			=	chierr[goodchi]
	
	fitopt, fitcov	=	curve_fit(chi_nonrel_rvm, tmsarr0, chiarr0, sigma=chierr0, absolute_sigma=True, p0=(0.1, 10000.0, -1.0, 1.0, 1.0))
	fiterr			=	np.sqrt(np.diag(fitcov))
	return (fitopt, fiterr)
	
#	-------------------------------------------------------------------------------

def chi_str_rvm(tms, t0ms, omeg, chi0, theta, inclin, beta):

	#	PA for STR RVM
	
	phi			=	(tms - t0ms)*omeg*1.0e-3 
	schi		=	( np.sin(theta)*np.sin(phi) ) + beta*( np.sin(inclin)*np.sin(theta) + np.cos(inclin)*np.cos(theta)*np.cos(phi) ) 
	cchi		=	( -np.sin(inclin)*np.cos(theta) ) + ( np.cos(inclin)*np.sin(theta)*np.cos(phi) ) - beta*( np.cos(theta)*np.sin(phi) ) 
	return (chi0 + np.arctan(schi/cchi))
	
#	-------------------------------------------------------------------------------

def fit_chi_str_rvm(tmsarr, chiarr, chierr):

	#	Fit PA vs time for non-relativisitc RVM
	
	goodchi			=	np.where(np.isfinite(chiarr))
	tmsarr0			=	tmsarr[goodchi]
	chiarr0			=	chiarr[goodchi]
	chierr0			=	chierr[goodchi]
	
	fitopt, fitcov	=	curve_fit(chi_str_rvm, tmsarr0, chiarr0, sigma=chierr0, absolute_sigma=True, p0=(0.08, 28000.0, -0.6, -0.09, 1.3, 0.6))
	fiterr			=	np.sqrt(np.diag(fitcov))
	return (fitopt, fiterr)
	
#	-------------------------------------------------------------------------------

def chi_gtr_rvm(tms, t0ms, omeg, chi0, theta, inclin, beta, ssssi):

	#	PA for GTR RVM
	
	phi			=	(tms - t0ms)*omeg*1.0e-3 
	constb		=	sin(i)*sin(theta) + cos(i)*cos(theta)*cos(phi)
	constc		=	const*()
	consta		=	0
	
	schi		=	( np.sin(theta)*np.sin(phi) ) + beta*consta
	cchi		=	( -np.sin(inclin)*np.cos(theta) ) + ( np.cos(inclin)*np.sin(theta)*np.cos(phi) ) - beta*np.sin(phi)*constc  
	return (chi0 + np.arctan(schi/cchi))
	
#	-------------------------------------------------------------------------------


































































































