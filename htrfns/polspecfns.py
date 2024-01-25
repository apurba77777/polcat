#
#	Functions for FRB polarization spectra
#
#								AB, June 2023

#	Function list
#
#		polambda(frbname, dm, nchan, ffac, avgfac, fmhz0, rm0, tbasems):
#				Plot FRB polarization against wavelength	
#
#
#
#	--------------------------	Import modules	---------------------------

import os, sys
import numpy as np
from globalpars import *
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import FuncFormatter
from matplotlib import colors
from scipy.optimize import curve_fit

mpl.rcParams['pdf.fonttype']	= 42
mpl.rcParams['ps.fonttype'] 	= 42
mpl.rcParams['savefig.dpi'] 	= 600
mpl.rcParams['font.family'] 	= 'sans-serif'
mpl.rcParams['font.size']		= 7

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

def polambda(frbname, dm, nchan, ffac, avgfac, fmhz0, rm0, tbasems):

	#	Plot FRB polarization against wavelength
	
	datadir		=	htrdir+frbname+'/data/reduced/'
	odatadir	=	htrdir+frbname+'/data/outputs/'
	plotdir		=	htrdir+frbname+'/plots/'
		
	stksdsfile	=	"{}{}_stks_ds_unfarot_{}_{}_avg_{}_{}_rm_{:.2f}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac,rm0)
	noisefile	=	"{}{}_noise_spec_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)	
	chanqfile	=	"{}{}_chanq_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)	
	avgchan		=	int(round(nchan/ffac))
	
	if(os.path.exists(chanqfile)):
		chanq	=	np.load(chanqfile)
	else:
		chanq	=	np.ones(avgchan, dtype=int)
		print("No bad channel list found...")		
		
	tresms		=	avgfac*nchan*Raw_time_res_ms
	
	if(os.path.exists(stksdsfile) and os.path.exists(noisefile)):
		stksds		=	np.load(stksdsfile)
		noispec		=	np.load(noisefile)
		badcs		=	np.where(chanq<1)
		stksds[:,badcs]	=	np.nan		
		its			=	np.nanmean(stksds[0],axis=0)
		tpeak		=	np.argmax(its)
		avgchans	=	stksds.shape[1]
		chanwmhz	=	NatBWmhz/avgchans
		fmhzarr		=	np.arange(fmhz0-(NatBWmhz/2)+(chanwmhz/2), fmhz0+NatBWmhz/2, chanwmhz)
		dl2marr		=	1.0e-16*( (ccC/fmhzarr)**2 - (ccC/fmhz0)**2)
		#dl2marr	=	1.0e-16*((ccC/fmhzarr)**2) 
		lmarr		=	(1.0e-8*ccC/fmhzarr)
		lm2arr		=	(1.0e-8*ccC/fmhzarr)**2
		
		goodc		=	np.where(chanq>0)
		
		tstart		=	tbasems[0]
		tend		=	tbasems[1]			
		istart		=	int(tpeak + (tstart/tresms))
		iend		=	int(tpeak + (tend/tresms) - 1)
		
		ispec		=	np.nanmean(stksds[0,:,istart:iend+1], axis=1)
		vspec		=	np.nanmean(stksds[3,:,istart:iend+1], axis=1)
		qspec		=	np.nanmean(stksds[1,:,istart:iend+1], axis=1)
		uspec		=	np.nanmean(stksds[2,:,istart:iend+1], axis=1)			
		noispec0	=	noispec/np.sqrt(float(iend+1-istart))
		lspec		=	np.sqrt(uspec**2 + qspec**2)
		dlspec		=	np.sqrt((uspec*noispec0[2])**2 + (qspec*noispec0[1])**2)/lspec
		pspec		=	np.sqrt(lspec**2 + vspec**2)
		dpspec		=	np.sqrt((vspec*dlspec)**2 + (lspec*noispec0[3])**2)/pspec
		
		qfracspec	=	qspec/ispec
		ufracspec	=	uspec/ispec
		vfracspec	=	vspec/ispec
		dqfrac		=	np.sqrt((qspec*noispec0[0])**2 + (ispec*noispec0[1])**2)/(ispec**2)
		dufrac		=	np.sqrt((uspec*noispec0[0])**2 + (ispec*noispec0[2])**2)/(ispec**2)
		dvfrac		=	np.sqrt((vspec*noispec0[0])**2 + (ispec*noispec0[3])**2)/(ispec**2)
		
		lfracspec	=	lspec/ispec
		dlfrac		=	np.sqrt((lspec*noispec0[0])**2 + (ispec*dlspec)**2)/(ispec**2)
		pfracspec	=	pspec/ispec
		dpfrac		=	np.sqrt((pspec*noispec0[0])**2 + (ispec*dpspec)**2)/(ispec**2)
		
		phispec		=	np.rad2deg(0.5*np.arctan2(uspec,qspec))		
		dphispec	=	np.rad2deg(0.5*np.sqrt((uspec*noispec0[1])**2 + (qspec*noispec0[2])**2) / (uspec**2 + qspec**2))
		
		phispec		=	phispec - 180.0*np.rint((phispec-np.nanmedian(phispec))/180.0)
		
		popt,pcov	=	curve_fit(chirmfit, dl2marr[goodc], phispec[goodc], sigma=dphispec[goodc], absolute_sigma=True)
		perr		=	np.sqrt(np.diag(pcov))
		fitrm		=	[popt[0],perr[0]]
		fitchi		=	[popt[1],perr[1]]
		
		print("Fit results --- RM = %.2f +/- %.2f rad/m2  PA = %.2f +/- %.2f deg (at lambda = %.2f m)"%(fitrm[0],fitrm[1],fitchi[0],fitchi[1],1.0e-8*ccC/fmhz0))
		#print("Fit results --- RM = %.2f +/- %.2f rad/m2  PA = %.2f +/- %.2f deg (at lambda = %.2f m)"%(fitrm[0],fitrm[1],fitchi[0],fitchi[1],0.0))
		
		psispec		=	np.rad2deg(0.5*np.arctan2(vspec,lspec))		
		dpsispec	=	np.rad2deg(0.5*np.sqrt((vspec*dlspec)**2 + (lspec*noispec0[2])**2) / (vspec**2 + lspec**2))		
		
		popt2,pcov2	=	curve_fit(chirmfit, dl2marr[goodc], psispec[goodc], sigma=dpsispec[goodc], absolute_sigma=True)
		perr2		=	np.sqrt(np.diag(pcov2))
		fitrm2		=	[popt2[0],perr2[0]]
		fitchi2		=	[popt2[1],perr2[1]]
		
		print("Fit results --- SRM = %.2f +/- %.2f rad/m2  EA = %.2f +/- %.2f deg (at lambda = %.2f m)"%(fitrm2[0],fitrm2[1],fitchi2[0],fitchi2[1],1.0e-8*ccC/fmhz0))
		
		angrange	=	max(np.nanmax(psispec[goodc])-np.nanmin(psispec[goodc]),np.nanmax(phispec[goodc])-np.nanmin(phispec[goodc])) + 1.0*max(np.nanmedian(dpsispec[goodc]),np.nanmedian(dphispec[goodc]))
		
		'''								
		fig			=	plt.figure(figsize=(6,8))				
		ax1			=	fig.add_subplot(3,1,1)
		ax1.set_title("Time range  %.2f $-$ %.2f ms"%(tstart,tend))
		ax1.errorbar(lmarr[goodc], qfracspec[goodc], dqfrac[goodc], fmt='ko', capsize=2, lw=0.1, markersize=6, label='Q/I')
		ax1.errorbar(lmarr[goodc], ufracspec[goodc], dufrac[goodc], fmt='rs', capsize=2, lw=0.1, fillstyle='none', markersize=4, label='U/I')
		ax1.errorbar(lmarr[goodc], vfracspec[goodc], dvfrac[goodc], fmt='md', capsize=2, lw=0.1, markersize=5, label='V/I')
		ax1.set_ylim([-1.2,1.2])
		plt.legend(loc='upper right', ncol=3)
		ax1.set_ylabel(r"Stokes paramaters")'''
		'''		
		ax1			=	fig.add_subplot(2,1,1)
		ax1.set_title("Time range  %.2f $-$ %.2f ms"%(tstart,tend))
		ax1.errorbar(lmarr[goodc], pfracspec[goodc], dpfrac[goodc], fmt='ko', capsize=2, lw=0.1, markersize=6, label='P/I')
		ax1.errorbar(lmarr[goodc], lfracspec[goodc], dlfrac[goodc], fmt='rs', capsize=2, lw=0.1, fillstyle='none', markersize=4, label='L/I')
		ax1.errorbar(lmarr[goodc], vfracspec[goodc], dvfrac[goodc], fmt='md', capsize=2, lw=0.1, markersize=5, label='V/I')
		ax1.set_ylim([-1.2,1.2])
		plt.legend(loc='lower right', ncol=3)
		ax1.set_ylabel(r"Polarization fraction")
		'''
		fig			=	plt.figure(figsize=(6,5))
		ax2,ax1	=	fig.subplots(2,1,height_ratios=[1,1],sharex=True)
		ax2.set_title("Time range  %.2f $-$ %.2f ms"%(tstart,tend))
		
		ax2.errorbar(lm2arr[goodc], phispec[goodc], dphispec[goodc], fmt='b*', capsize=2, lw=0.5, fillstyle='none', markersize=5)
		ax2.plot(lm2arr, chirmfit(dl2marr, *popt), 'k--', lw=2)
		ax2.set_ylim([0.5*(np.nanmin(phispec[goodc]) + np.nanmax(phispec[goodc])) - angrange, 0.5*(np.nanmin(phispec[goodc]) + np.nanmax(phispec[goodc])) + angrange])	
		ax2.set_ylabel(r"$\chi$ (deg)")
		
		ax1.errorbar(lm2arr[goodc], psispec[goodc], dpsispec[goodc], fmt='b*', capsize=2, lw=0.5, fillstyle='none', markersize=5)
		ax1.plot(lm2arr, chirmfit(dl2marr, *popt2), 'k--', lw=2)
		ax1.set_ylim([0.5*(np.nanmin(psispec[goodc]) + np.nanmax(psispec[goodc])) - angrange, 0.5*(np.nanmin(psispec[goodc]) + np.nanmax(psispec[goodc])) + angrange])			
		ax1.set_ylabel(r"$\psi$ (deg)")
		ax1.set_xlabel(r"$\lambda^2$ (m$^2$)")	
								
		plt.tight_layout(h_pad=0, w_pad=0)			
		fig.align_ylabels()	
		#plt.savefig("{}{}_t_{}_to_{}_ms_dm_{}_chan_{}_avg_{}_{}.png".format(plotdir,frbname,tbasems[0],tbasems[1],dm,nchan,ffac,avgfac))
		#plt.savefig("{}{}_t_{}_to_{}_ms_dm_{}_chan_{}_avg_{}_{}.pdf".format(plotdir,frbname,tbasems[0],tbasems[1],dm,nchan,ffac,avgfac))
		plt.show()
		
		del(stksds)
		del(noispec)
	else:
		print("polambda ---- PANIC - File not found!")

	return(0)

#	-------------------------------------------------------------------------------

def stokespec(frbname, dm, nchan, ffac, avgfac, fmhz0, rm0, tbasems):

	#	Plot Stokes parameters against frequency
	
	datadir		=	htrdir+frbname+'/data/reduced/'
	odatadir	=	htrdir+frbname+'/data/outputs/'
	plotdir		=	htrdir+frbname+'/plots/'
		
	stksdsfile	=	"{}{}_stks_ds_unfarot_{}_{}_avg_{}_{}_rm_{:.2f}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac,rm0)
	noisefile	=	"{}{}_noise_spec_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)	
	chanqfile	=	"{}{}_chanq_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)	
	avgchan		=	int(round(nchan/ffac))
	
	if(os.path.exists(chanqfile)):
		chanq	=	np.load(chanqfile)
	else:
		chanq	=	np.ones(avgchan, dtype=int)
		print("No bad channel list found...")		
		
	tresms		=	avgfac*nchan*Raw_time_res_ms
	
	if(os.path.exists(stksdsfile) and os.path.exists(noisefile)):
		stksds		=	np.load(stksdsfile)
		noispec		=	np.load(noisefile)
		badcs		=	np.where(chanq<1)
		stksds[:,badcs]	=	np.nan		
		its			=	np.nanmean(stksds[0],axis=0)
		tpeak		=	np.argmax(its)
		avgchans	=	stksds.shape[1]
		chanwmhz	=	NatBWmhz/avgchans
		fmhzarr		=	np.arange(fmhz0-(NatBWmhz/2)+(chanwmhz/2), fmhz0+NatBWmhz/2, chanwmhz)
		dl2marr		=	1.0e-16*( (ccC/fmhzarr)**2 - (ccC/fmhz0)**2)
		#dl2marr	=	1.0e-16*((ccC/fmhzarr)**2) 
		lmarr		=	(1.0e-8*ccC/fmhzarr)
		lm2arr		=	(1.0e-8*ccC/fmhzarr)**2
		
		goodc		=	np.where(chanq>0)
		
		tstart		=	tbasems[0]
		tend		=	tbasems[1]			
		istart		=	int(tpeak + (tstart/tresms))
		iend		=	int(tpeak + (tend/tresms) - 1)
		
		ispec		=	np.nanmean(stksds[0,:,istart:iend+1], axis=1)
		vspec		=	np.nanmean(stksds[3,:,istart:iend+1], axis=1)
		qspec		=	np.nanmean(stksds[1,:,istart:iend+1], axis=1)
		uspec		=	np.nanmean(stksds[2,:,istart:iend+1], axis=1)			
		noispec0	=	noispec/np.sqrt(float(iend+1-istart))
		lspec		=	np.sqrt(uspec**2 + qspec**2)
		dlspec		=	np.sqrt((uspec*noispec0[2])**2 + (qspec*noispec0[1])**2)/lspec
		pspec		=	np.sqrt(lspec**2 + vspec**2)
		dpspec		=	np.sqrt((vspec*dlspec)**2 + (lspec*noispec0[3])**2)/pspec
		
		qfracspec	=	qspec/ispec
		ufracspec	=	uspec/ispec
		vfracspec	=	vspec/ispec
		dqfrac		=	np.sqrt((qspec*noispec0[0])**2 + (ispec*noispec0[1])**2)/(ispec**2)
		dufrac		=	np.sqrt((uspec*noispec0[0])**2 + (ispec*noispec0[2])**2)/(ispec**2)
		dvfrac		=	np.sqrt((vspec*noispec0[0])**2 + (ispec*noispec0[3])**2)/(ispec**2)
		
		lfracspec	=	lspec/ispec
		dlfrac		=	np.sqrt((lspec*noispec0[0])**2 + (ispec*dlspec)**2)/(ispec**2)
		pfracspec	=	pspec/ispec
		dpfrac		=	np.sqrt((pspec*noispec0[0])**2 + (ispec*dpspec)**2)/(ispec**2)
		
		phispec		=	np.rad2deg(0.5*np.arctan2(uspec,qspec))		
		dphispec	=	np.rad2deg(0.5*np.sqrt((uspec*noispec0[1])**2 + (qspec*noispec0[2])**2) / (uspec**2 + qspec**2))
		
		phispec		=	phispec - 180.0*np.rint((phispec-np.nanmedian(phispec))/180.0)
		
		popt,pcov	=	curve_fit(chirmfit, dl2marr[goodc], phispec[goodc], sigma=dphispec[goodc], absolute_sigma=True)
		perr		=	np.sqrt(np.diag(pcov))
		fitrm		=	[popt[0],perr[0]]
		fitchi		=	[popt[1],perr[1]]
		
		print("Fit results --- RM = %.2f +/- %.2f rad/m2  PA = %.2f +/- %.2f deg (at lambda = %.2f m)"%(fitrm[0],fitrm[1],fitchi[0],fitchi[1],1.0e-8*ccC/fmhz0))
		#print("Fit results --- RM = %.2f +/- %.2f rad/m2  PA = %.2f +/- %.2f deg (at lambda = %.2f m)"%(fitrm[0],fitrm[1],fitchi[0],fitchi[1],0.0))
		
		psispec		=	np.rad2deg(0.5*np.arctan2(vspec,lspec))		
		dpsispec	=	np.rad2deg(0.5*np.sqrt((vspec*dlspec)**2 + (lspec*noispec0[2])**2) / (vspec**2 + lspec**2))		
		
		popt2,pcov2	=	curve_fit(chirmfit, dl2marr[goodc], psispec[goodc], sigma=dpsispec[goodc], absolute_sigma=True)
		perr2		=	np.sqrt(np.diag(pcov2))
		fitrm2		=	[popt2[0],perr2[0]]
		fitchi2		=	[popt2[1],perr2[1]]
		
		print("Fit results --- SRM = %.2f +/- %.2f rad/m2  EA = %.2f +/- %.2f deg (at lambda = %.2f m)"%(fitrm2[0],fitrm2[1],fitchi2[0],fitchi2[1],1.0e-8*ccC/fmhz0))
									
		fig			=	plt.figure(figsize=(5.4,5.6))
		ax1,ax2		=	fig.subplots(2,1,height_ratios=[1.5,1],sharex=True)				
		ax1.set_title("FRB 20%sA - sub-burst A"%(frbname))		
		ax1.fill_between(fmhzarr/1.0e3,qfracspec-dqfrac,qfracspec+dqfrac,color='coral',alpha=0.5)
		ax1.fill_between(fmhzarr/1.0e3,ufracspec-dufrac,ufracspec+dufrac,color='lightgrey',alpha=0.6)
		ax1.fill_between(fmhzarr/1.0e3,vfracspec-dvfrac,vfracspec+dvfrac,color='aqua',alpha=0.6)
		ax1.plot(fmhzarr/1.0e3, qfracspec, 'rd', markersize=5, lw=0.5, label='$q$')
		ax1.plot(fmhzarr/1.0e3, vfracspec, 'b*', markersize=6, lw=0.5, label='$v$')
		ax1.plot(fmhzarr/1.0e3, ufracspec, 'ko', markersize=5,lw=0.5, label='$u$')
		
		ax1.set_ylim([-1.2,1.2])
		#plt.legend(loc='upper right', ncol=3)
		ax1.set_ylabel(r"Fractional Q, U, V")
		
		ax2.errorbar(fmhzarr[goodc]/1.0e3, phispec[goodc], dphispec[goodc], fmt='g*', capsize=2, lw=0.5, fillstyle='none', markersize=5)
		ax2.plot(fmhzarr/1.0e3, chirmfit(dl2marr, *popt), 'k--', lw=2, label='$\delta$RM = ($%.1f \pm %.1f$) rad m$^{-2}$'%(fitrm[0],fitrm[1]))
		ax2.set_ylim([-70,-70+80])
		ax2.set_yticks(np.arange(-60,-70+80, 20))
		ax2.set_ylabel(r"PA (deg)")
		ax2.set_xlabel(r"Frequency (GHz)")
		ax2.legend(loc='upper right')
												
		plt.tight_layout(h_pad=0, w_pad=0)			
		fig.align_ylabels()	
		plt.savefig("{}{}_t_{}_to_{}_stkspec_{}_chan_{}_avg_{}_{}.png".format(plotdir,frbname,tbasems[0],tbasems[1],dm,nchan,ffac,avgfac))
		plt.savefig("{}{}_t_{}_to_{}_stkspec_{}_chan_{}_avg_{}_{}.pdf".format(plotdir,frbname,tbasems[0],tbasems[1],dm,nchan,ffac,avgfac))
		plt.show()
		
		del(stksds)
		del(noispec)
	else:
		print("stokespec ---- PANIC - File not found!")

	return(0)

#	-------------------------------------------------------------------------------

def gamdelam(frbname, dm, nchan, ffac, avgfac, fmhz0, rm0, tbasems, refangdeg):

	#	Plot gamma and delta against wavelength
	
	datadir		=	htrdir+frbname+'/data/reduced/'
	odatadir	=	htrdir+frbname+'/data/outputs/'
	plotdir		=	htrdir+frbname+'/plots/'
		
	stksdsfile	=	"{}{}_stks_ds_unfarot_{}_{}_avg_{}_{}_rm_{:.2f}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac,rm0)
	noisefile	=	"{}{}_noise_spec_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)	
	chanqfile	=	"{}{}_chanq_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)	
	avgchan		=	int(round(nchan/ffac))
	
	if(os.path.exists(chanqfile)):
		chanq	=	np.load(chanqfile)
	else:
		chanq	=	np.ones(avgchan, dtype=int)
		print("No bad channel list found...")		
		
	tresms		=	avgfac*nchan*Raw_time_res_ms
	
	if(os.path.exists(stksdsfile) and os.path.exists(noisefile)):
		stksds		=	np.load(stksdsfile)
		noispec		=	np.load(noisefile)
		badcs		=	np.where(chanq<1)
		stksds[:,badcs]	=	np.nan		
		its			=	np.nanmean(stksds[0],axis=0)
		tpeak		=	np.argmax(its)
		avgchans	=	stksds.shape[1]
		chanwmhz	=	NatBWmhz/avgchans
		fmhzarr		=	np.arange(fmhz0-(NatBWmhz/2)+(chanwmhz/2), fmhz0+NatBWmhz/2, chanwmhz)
		dl2marr		=	1.0e-16*( (ccC/fmhzarr)**2 - (ccC/fmhz0)**2)
		#dl2marr	=	1.0e-16*((ccC/fmhzarr)**2) 
		lmarr		=	(1.0e-8*ccC/fmhzarr)
		lm2arr		=	(1.0e-8*ccC/fmhzarr)**2
		
		goodc		=	np.where(chanq>0)
		
		tstart		=	tbasems[0]
		tend		=	tbasems[1]			
		istart		=	int(tpeak + (tstart/tresms))
		iend		=	int(tpeak + (tend/tresms) - 1)
		
		ispec		=	np.nanmean(stksds[0,:,istart:iend+1], axis=1)
		vspec		=	np.nanmean(stksds[3,:,istart:iend+1], axis=1)
		qspec		=	np.nanmean(stksds[1,:,istart:iend+1], axis=1)
		uspec		=	np.nanmean(stksds[2,:,istart:iend+1], axis=1)			
		noispec0	=	noispec/np.sqrt(float(iend+1-istart))
		lspec		=	np.sqrt(uspec**2 + qspec**2)
		dlspec		=	np.sqrt((uspec*noispec0[2])**2 + (qspec*noispec0[1])**2)/lspec
		pspec		=	np.sqrt(lspec**2 + vspec**2)
		dpspec		=	np.sqrt((vspec*dlspec)**2 + (lspec*noispec0[3])**2)/pspec
		
		cosang			=	np.cos(2*refangdeg*np.pi/180.0)
		sinang			=	np.sin(2*refangdeg*np.pi/180.0)
		
		qsrc			=	qspec*cosang + uspec*sinang
		usrc			=	-qspec*sinang + uspec*cosang
		
		uvspec			=	np.sqrt(usrc**2 + vspec**2)		
		duvspec			=	np.sqrt((vspec*noispec0[3])**2 + (usrc*noispec0[2])**2)/uvspec
		
		lspec			=	np.sqrt(usrc**2 + qsrc**2)
		dlspec			=	np.sqrt((usrc*noispec0[2])**2 + (qsrc*noispec0[1])**2)/lspec
		
		phispec		=	np.rad2deg(0.5*np.arctan2(usrc,qsrc))		
		dphispec	=	np.rad2deg(0.5*np.sqrt((usrc*noispec0[1])**2 + (qsrc*noispec0[2])**2) / (usrc**2 + qsrc**2))		
		phispec		=	phispec - 180.0*np.rint((phispec-np.nanmedian(phispec))/180.0)
		
		psispec		=	np.rad2deg(0.5*np.arctan2(vspec,lspec))		
		dpsispec	=	np.rad2deg(0.5*np.sqrt((vspec*dlspec)**2 + (lspec*noispec0[3])**2) / (vspec**2 + lspec**2))		
		
		delspec			=	np.rad2deg(np.arctan2(vspec,usrc))		
		ddelspec		=	np.rad2deg(np.sqrt((usrc*noispec0[3])**2 + (vspec*noispec0[2])**2) / (usrc**2 + vspec**2))	
		delspec			=	delspec + 180.0
		
		gamspec			=	np.rad2deg(0.5*np.arctan2(uvspec,qsrc))		
		dgamspec		=	np.rad2deg(0.5*np.sqrt((qsrc*duvspec)**2 + (uvspec*noispec0[1])**2) / (qsrc**2 + uvspec**2))
		
		fig				=	plt.figure(figsize=(5,10))
		ax2,ax1,ax3,ax4	=	fig.subplots(4,1,height_ratios=[1,1,1,1],sharex=True)
		ax2.set_title("Time range  %.2f $-$ %.2f ms"%(tstart,tend))
		
		ax2.errorbar(lmarr[goodc], phispec[goodc], dphispec[goodc], fmt='b*', capsize=2, lw=0.5, fillstyle='none', markersize=5)
		ax2.set_ylim([np.nanmedian(phispec[goodc])-45, np.nanmedian(phispec[goodc])+45])
		ax2.set_ylabel(r"$\chi$' (deg)")
		
		ax1.errorbar(lmarr[goodc], psispec[goodc], dpsispec[goodc], fmt='b*', capsize=2, lw=0.5, fillstyle='none', markersize=5)
		ax1.set_ylim([np.nanmedian(psispec[goodc])-45, np.nanmedian(psispec[goodc])+45])
		ax1.set_ylabel(r"$\psi$' (deg)")
		
		ax3.errorbar(lmarr[goodc], gamspec[goodc], dgamspec[goodc], fmt='b*', capsize=2, lw=0.5, fillstyle='none', markersize=5)
		ax3.set_ylim([np.nanmedian(gamspec[goodc])-45, np.nanmedian(gamspec[goodc])+45])
		ax3.set_ylabel(r"$\gamma$ (deg)")
		
		ax4.errorbar(lmarr[goodc], delspec[goodc], ddelspec[goodc], fmt='b*', capsize=2, lw=0.5, fillstyle='none', markersize=5)
		ax4.set_ylim([np.nanmedian(delspec[goodc])-45, np.nanmedian(delspec[goodc])+45])
		ax4.set_ylabel(r"$\delta$ (deg)")
				
		ax4.set_xlabel(r"$\lambda$ (m)")	
								
		plt.tight_layout(h_pad=0, w_pad=0)			
		fig.align_ylabels()	
		plt.savefig("{}{}_t_{}_to_{}_ms_dm_{}_chan_{}_avg_{}_{}.png".format(plotdir,frbname,tbasems[0],tbasems[1],dm,nchan,ffac,avgfac))
		plt.savefig("{}{}_t_{}_to_{}_ms_dm_{}_chan_{}_avg_{}_{}.pdf".format(plotdir,frbname,tbasems[0],tbasems[1],dm,nchan,ffac,avgfac))
		plt.show()
		
		del(stksds)
		del(noispec)
	else:
		print("polambda ---- PANIC - File not found!")

	return(0)

#	-------------------------------------------------------------------------------






















































































