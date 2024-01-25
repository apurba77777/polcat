#
#	Plotting functions
#
#								AB, October 2023
#
#	Function list
#
#			plot_stkspec(plotdir,frbdat,ylim,fsize):
#					Plot fractional IQUV spectra
#
#			plot_polspec(plotdir,frbdat,ylim,fsize,palim):
#					Plot polarization fraction spectra
#
#			plot_spechammer(plotdir,frbdat,fsize):
#					Plot Poincare sphere frequency track in Hammer projection
#
#	--------------------------	Import modules	---------------------------

import os, sys
import numpy as np
from globalpars import *
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import FuncFormatter
import matplotlib.colors as mpc
import matplotlib.ticker as ticker
from scipy.optimize import curve_fit
from polfns.polmodelfns import *

mpl.rcParams['pdf.fonttype']	= 42
mpl.rcParams['ps.fonttype'] 	= 42
mpl.rcParams['savefig.dpi'] 	= 600
mpl.rcParams['font.family'] 	= 'sans-serif'
mpl.rcParams['font.size']		= 7

#frbhtr		=	namedtuple('frbhtr',['frbname','dm','nchan','ffac','tavgfac','rm','tresms','twinms','subband','refdeg','tmsarr','fmhzarr',\
		#							 'irms','qrms','urms','vrms',\	
		#							 'it','qt','ut','vt','pt','lt','epts','elts','qfract','ufract','vfract','eqfract','eufract','evfract','lfract','pfract','elfract','epfract',\
		#							 'chit','psit','echit','epsit','chispec','echispec','psispec','epsispec',\
		#							 'ispec','qspec','uspec','vspec','eispec','eqspec','euspec','evspec','lspec','pspec','elspec','epspec',\
		#							 'qfracspec','ufracspec','vfracspec','dqfracspec','dufracspec','dvfracspec','lfracspec','pfracspec','dlfracspec','dpfracspec', \
		#							 'ids','qds','uds','vds','irmspec','qrmspec','urmspec','vrmspec'])

#	----------------------------------------------------------------------------------------------------------

def plot_stkspec(plotdir,frbdat,ylim,fsize):

	#	Plot fractional IQUV spectra
	
	fig 	= plt.figure(figsize=(fsize[0], fsize[1]))
	ax 		= fig.add_axes([0.12, 0.12, 0.86,0.86])
	ax.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
	
	ax.axhline(c='c',ls='--')	
	ax.fill_between(frbdat.fmhzarr/1.0e3,frbdat.qfracspec-frbdat.dqfracspec,frbdat.qfracspec+frbdat.dqfracspec,color='coral',alpha=0.5)
	ax.fill_between(frbdat.fmhzarr/1.0e3,frbdat.ufracspec-frbdat.dufracspec,frbdat.ufracspec+frbdat.dufracspec,color='plum',alpha=0.6)
	ax.fill_between(frbdat.fmhzarr/1.0e3,frbdat.vfracspec-frbdat.dvfracspec,frbdat.vfracspec+frbdat.dvfracspec,color='aqua',alpha=0.6)
	ax.plot(frbdat.fmhzarr/1.0e3, frbdat.qfracspec, 'r-.',lw=2, label='Q/I')
	ax.plot(frbdat.fmhzarr/1.0e3, frbdat.ufracspec, 'm--',lw=2, label='U/I')
	ax.plot(frbdat.fmhzarr/1.0e3, frbdat.vfracspec, 'b-',lw=2, label='V/I')
	
	ax.set_ylim(ylim)
	ax.legend(loc='upper right', ncol=3)	
	ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
	ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
	ax.set_ylabel(r'Fraction')	
	ax.set_xlabel(r'$\nu$ (GHz)')
	
	plt.savefig("{}{}_stkspec_sb_{}_of_{}_rm_{:.2f}.pdf".format(plotdir,frbdat.frbname,frbdat.subband[0],frbdat.subband[1],frbdat.rm))
	plt.show()

	return(0)

#	----------------------------------------------------------------------------------------------------------

def plot_polspec(plotdir,frbdat,ylim,fsize,palim):

	#	Plot polarization fraction spectra
	
	fig 	= plt.figure(figsize=(fsize[0]/intocm, fsize[1]/intocm))
	ax 		= fig.add_axes([0.15, 0.48, 0.83,0.50])
	ax.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
	#ax.text(0.04,0.06,'e',fontsize=8,fontweight='bold',transform=ax.transAxes)
	
	ax.axhline(c='c',ls='--', lw=0.25)	
	ax.fill_between(frbdat.fmhzarr/1.0e3,frbdat.lfracspec-frbdat.dlfracspec,frbdat.lfracspec+frbdat.dlfracspec,color='coral',alpha=0.5)
	ax.fill_between(frbdat.fmhzarr/1.0e3,frbdat.pfracspec-frbdat.dpfracspec,frbdat.pfracspec+frbdat.dpfracspec,color='plum',alpha=0.6)
	ax.fill_between(frbdat.fmhzarr/1.0e3,frbdat.vfracspec-frbdat.dvfracspec,frbdat.vfracspec+frbdat.dvfracspec,color='aqua',alpha=0.6)
	ax.plot(frbdat.fmhzarr/1.0e3, frbdat.lfracspec, 'r:',lw=1, label='L/I')
	ax.plot(frbdat.fmhzarr/1.0e3, frbdat.pfracspec, 'm-',lw=1, label='P/I')
	ax.plot(frbdat.fmhzarr/1.0e3, frbdat.vfracspec, 'b--',lw=1, label='V/I')
	
	ax.set_ylim(ylim)
	ax.legend(loc='upper right', ncol=3)	
	ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
	ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
	ax.set_ylabel(r'Fraction')	
	ax.set_xticklabels([])
	ax.yaxis.set_label_coords(-0.11, 0.5)
	
	phi			=	frbdat.chispec
	dphi		=	frbdat.echispec
	badl		=	np.where((frbdat.lspec/frbdat.elspec) < 3.0)
	#nongoodphi	=	np.union1d(np.where(dphi > 50.0),badl)
	#phi[nongoodphi]		=	np.nan
	#dphi[nongoodphi]	=	np.nan
	
	goodc		=	np.isfinite(phi)	
	fmhz0		=	np.nanmedian(frbdat.fmhzarr)
	dl2marr		=	1.0e-16*( (ccC/frbdat.fmhzarr)**2 - (ccC/fmhz0)**2)
	
	popt,pcov	=	curve_fit(chirmfit, dl2marr[goodc], phi[goodc], sigma=dphi[goodc], absolute_sigma=True)
	perr		=	np.sqrt(np.diag(pcov))
	fitrm		=	[popt[0],perr[0]]
	fitchi		=	[popt[1],perr[1]]
		
	
	ax1 		= fig.add_axes([0.15, 0.12, 0.83,0.36])
	ax1.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
	#ax1.text(0.04,0.1,'f',fontsize=8,fontweight='bold',transform=ax1.transAxes)
	#ax1.fill_between(frbdat.fmhzarr/1.0e3,phi-dphi, phi+dphi,color='coral',alpha=0.6)
	ax1.errorbar(frbdat.fmhzarr/1.0e3, phi, dphi, fmt='ro', markersize=3,lw=0.5,capsize=2)
	ax1.plot(frbdat.fmhzarr/1.0e3, chirmfit(dl2marr, *popt), 'k--', lw=1, label='$\delta$RM = ($%.1f \pm %.1f$) rad m$^{-2}$'%(fitrm[0],fitrm[1]))
	ax1.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
	ax1.set_ylim([palim[0],palim[1]])
	ax1.yaxis.set_major_locator(ticker.MultipleLocator(palim[2]))
	ax1.set_xlabel(r'$\nu$ (GHz)')
	ax1.set_ylabel(r'PA (deg)')
	ax1.yaxis.set_label_coords(-0.11, 0.5)
	ax1.legend(loc='upper right')
	
	fig.align_ylabels()
	
	plt.savefig("{}{}_polspec_sb_{}_of_{}_rm_{:.2f}.pdf".format(plotdir,frbdat.frbname,frbdat.subband[0],frbdat.subband[1],frbdat.rm), transparent=True, format='pdf')
	plt.savefig("{}{}_polspec_sb_{}_of_{}_rm_{:.2f}.eps".format(plotdir,frbdat.frbname,frbdat.subband[0],frbdat.subband[1],frbdat.rm), transparent=True, format='eps')
	plt.show()

	return(0)

#	----------------------------------------------------------------------------------------------------------

def plot_spechammer(plotdir,frbdat,fsize):

	#	Plot Poincare sphere frequency track in Hammer projection
	
	fig 	= plt.figure(figsize=(fsize[0], fsize[1]))
	ax 		= fig.add_axes([0.12, 0.12, 0.86,0.86],projection="hammer")
	plt.grid(True)
	
	phi			=	frbdat.chispec
	dphi		=	frbdat.echispec
	badl		=	np.where((frbdat.lspec/frbdat.elspec) < 3.0)
	nongoodphi	=	np.union1d(np.where(dphi > 5.0),badl)
	phi[nongoodphi]		=	np.nan
	dphi[nongoodphi]	=	np.nan
	
	ax.errorbar(np.deg2rad(2*phi),np.deg2rad(2*phi),yerr=np.deg2rad(2*dphi),xerr=np.deg2rad(2*dphi),fmt='bo')
	
	'''ax.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
	
	ax.axhline(c='c',ls='--')	
	ax.fill_between(frbdat.fmhzarr/1.0e3,frbdat.qfracspec-frbdat.dqfracspec,frbdat.qfracspec+frbdat.dqfracspec,color='coral',alpha=0.5)
	ax.fill_between(frbdat.fmhzarr/1.0e3,frbdat.ufracspec-frbdat.dufracspec,frbdat.ufracspec+frbdat.dufracspec,color='plum',alpha=0.6)
	ax.fill_between(frbdat.fmhzarr/1.0e3,frbdat.vfracspec-frbdat.dvfracspec,frbdat.vfracspec+frbdat.dvfracspec,color='aqua',alpha=0.6)
	ax.plot(frbdat.fmhzarr/1.0e3, frbdat.qfracspec, 'r-.',lw=2, label='Q/I')
	ax.plot(frbdat.fmhzarr/1.0e3, frbdat.ufracspec, 'm--',lw=2, label='U/I')
	ax.plot(frbdat.fmhzarr/1.0e3, frbdat.vfracspec, 'b-',lw=2, label='V/I')
	
	ax.set_ylim(ylim)
	ax.legend(loc='upper right', ncol=3)	
	ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
	ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
	ax.set_ylabel(r'Fraction')	
	ax.set_xlabel(r'$\nu$ (GHz)')'''
	
	#plt.savefig("{}{}_stkspec_sb_{}_of_{}_rm_{:.2f}.pdf".format(plotdir,frbdat.frbname,frbdat.subband[0],frbdat.subband[1],frbdat.rm))
	plt.show()

	return(0)

#	----------------------------------------------------------------------------------------------------------














