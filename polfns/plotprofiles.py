#
#	Plotting functions
#
#								AB, October 2023
#
#	Function list
#
#			plot_iquvt(plotdir,frbdat,fmhz0,xlim,fsize,ttms):
#					Plot IQUV profiles and dynspecs
#
#			plot_iquvtsep(plotdir,frbdat,fmhz0,xlim,ylim,fsize,ttms):
#					Plot IQUV profiles and dynspecs separately
#
#			plot_ilvt(plotdir,frbdat,bfrbdat,xlim,fsize,ttms):
#					Plot ILV profiles
#
#			plot_ilvpat(plotdir,frbdat,bfrbdat,xlim,fsize, ttms, ttpa):
#					Plot ILV and PA profiles
#
#			plot_ilvparm(plotdir,rmfile,frbdat,bfrbdat,xlim,fsize,ttms,ttrm,ttpa):
#					Plot ILV, PA & RM profiles
#
#			plot_fracparm(plotdir,rmfile,frbdat,bfrbdat,xlim,fsize,ttms,ttrm):
#					Plot polarization fraction, PA & RM profiles
#
#			plot_fracpa(plotdir,rmfile,frbdat,bfrbdat,xlim,fsize,ttms,ttrm):
#					Plot polarization fraction & PA profiles
#
#			congaussfitter(plotdir,bfrbdat,xlim,fsize,ttms,ng,pol):
#					Fit profile with n Gaussians convolved with the same exponential
#
#			plot_irm(plotdir,rmfile,frbdat,bfrbdat,xlim,fsize,ttms,ttrm, baserm):
#					Plot I & RM profiles
#
#			plot_dpa(plotdir,rmfile,frbdat,bfrbdat,xlim,fsize,ttms,ntp):
#					Plot PA profile and dPA/dt
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
from scipy.optimize import minimize_scalar
from scipy.optimize import root_scalar
from polfns.polmodelfns import *
from polfns.promodels import *

mpl.rcParams['pdf.fonttype']	= 42
mpl.rcParams['ps.fonttype'] 	= 42
mpl.rcParams['savefig.dpi'] 	= 600
mpl.rcParams['font.family'] 	= 'sans-serif'
mpl.rcParams['font.size']		= 8

#frbhtr		=	namedtuple('frbhtr',['frbname','dm','nchan','ffac','tavgfac','rm','tresms','twinms','subband','refdeg','tmsarr','fmhzarr',\
		#							 'irms','qrms','urms','vrms',\	
		#							 'it','qt','ut','vt','pt','lt','epts','elts','qfract','ufract','vfract','eqfract','eufract','evfract','lfract','pfract','elfract','epfract',\
		#							 'chit','psit','echit','epsit','chispec','echispec','psispec','epsispec',\
		#							 'ispec','qspec','uspec','vspec','eispec','eqspec','euspec','evspec','lspec','pspec','elspec','epspec',\
		#							 'qfracspec','ufracspec','vfracspec','dqfracspec','dufracspec','dvfracspec','lfracspec','pfracspec','dlfracspec','dpfracspec', \
		#							 'ids','qds','uds','vds','irmspec','qrmspec','urmspec','vrmspec'])

#	----------------------------------------------------------------------------------------------------------

def plot_iquvt(plotdir,frbdat,fmhz0,xlim,fsize,ttms):

	#	Plot IQUV profiles and dynspecs
	
	avgchan		=	int(round(frbdat.nchan/frbdat.ffac))
	
	chanwmhz	=	NatBWmhz/avgchan
	fmhzarr		=	np.arange(fmhz0-(NatBWmhz/2)+(chanwmhz/2), fmhz0+NatBWmhz/2, chanwmhz)
	
	fig		=	plt.figure(figsize=(fsize[0]/intocm,fsize[1]/intocm))
	ax 		= 	fig.add_axes([0.10, 0.70, 0.88,0.28])
	ax.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
	
	ax.axhline(c='c',ls='--', lw=0.25)
	ax.plot(frbdat.tmsarr, frbdat.it/np.nanmax(frbdat.it), 'k-', lw=0.5, label='I')
	ax.plot(frbdat.tmsarr, frbdat.qt/np.nanmax(frbdat.it), 'r-', lw=0.5, label='Q')
	ax.plot(frbdat.tmsarr, frbdat.ut/np.nanmax(frbdat.it), 'm-', lw=0.5, label='U')
	ax.plot(frbdat.tmsarr, frbdat.vt/np.nanmax(frbdat.it), 'b-', lw=0.5, label='V')
	ax.text(0.02,0.9,'a',fontsize=8,fontweight='bold',transform=ax.transAxes)
	ax.set_ylim(ymax=1.1)
	ax.set_xlim(xlim)
	ax.legend(loc='upper right', ncol=2)
	ax.set_ylabel(r'Normalized flux density')
	ax.set_xticklabels([])
	ax.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
	ax.yaxis.set_label_coords(-0.07, 0.5)
		
	ax0 		= fig.add_axes([0.10, 0.54, 0.88,0.16])
	ax0.imshow(frbdat.ids, aspect='auto',cmap='seismic',interpolation="none",vmin=-np.nanmax(np.abs(frbdat.ids)),vmax=np.nanmax(np.abs(frbdat.ids)), \
		extent=([xlim[0],xlim[1],(fmhzarr[-1]-chanwmhz)/1.0e3, (fmhzarr[0]+chanwmhz)/1.0e3]))
	ax0.text(0.02,0.8,'b',fontsize=8,fontweight='bold',transform=ax0.transAxes)
	ax0.set_xticklabels([])
	ax0.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax0.set_ylabel(r'$\nu$ (GHz)')
	ax0.yaxis.set_label_coords(-0.07, 0.5)
	
	ax1 		= fig.add_axes([0.10, 0.38, 0.88,0.16])
	ax1.imshow(frbdat.qds, aspect='auto',cmap='seismic',interpolation="none",vmin=-np.nanmax(np.abs(frbdat.qds)),vmax=np.nanmax(np.abs(frbdat.qds)), \
		extent=([xlim[0],xlim[1],(fmhzarr[-1]-chanwmhz)/1.0e3, (fmhzarr[0]+chanwmhz)/1.0e3]))
	ax1.text(0.02,0.8,'c',fontsize=8,fontweight='bold',transform=ax1.transAxes)
	ax1.set_xticklabels([])
	ax1.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax1.set_ylabel(r'$\nu$ (GHz)')
	ax1.yaxis.set_label_coords(-0.07, 0.5)
	
	ax2 		= fig.add_axes([0.10, 0.22, 0.88,0.16])
	ax2.imshow(frbdat.uds, aspect='auto',cmap='seismic',interpolation="none",vmin=-np.nanmax(np.abs(frbdat.uds)),vmax=np.nanmax(np.abs(frbdat.uds)), \
		extent=([xlim[0],xlim[1],(fmhzarr[-1]-chanwmhz)/1.0e3, (fmhzarr[0]+chanwmhz)/1.0e3]))
	ax2.text(0.02,0.8,'d',fontsize=8,fontweight='bold',transform=ax2.transAxes)
	ax2.set_xticklabels([])
	ax2.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax2.set_ylabel(r'$\nu$ (GHz)')
	ax2.yaxis.set_label_coords(-0.07, 0.5)
	
	ax3 		= fig.add_axes([0.10, 0.06, 0.88,0.16])
	ax3.imshow(frbdat.vds, aspect='auto',cmap='seismic',interpolation="none",vmin=-np.nanmax(np.abs(frbdat.vds)),vmax=np.nanmax(np.abs(frbdat.vds)), \
		extent=([xlim[0],xlim[1],(fmhzarr[-1]-chanwmhz)/1.0e3, (fmhzarr[0]+chanwmhz)/1.0e3]))
	ax3.text(0.02,0.8,'e',fontsize=8,fontweight='bold',transform=ax3.transAxes)
	ax3.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax3.set_xlabel(r'Time (ms)')
	ax3.set_ylabel(r'$\nu$ (GHz)')
	ax3.yaxis.set_label_coords(-0.07, 0.5)
		
	plt.savefig("{}{}_iquv_sb_{}_of_{}_rm_{:.2f}.eps".format(plotdir,frbdat.frbname,frbdat.subband[0],frbdat.subband[1],frbdat.rm), transparent=True, format='eps')
	plt.savefig("{}{}_iquv_sb_{}_of_{}_rm_{:.2f}.pdf".format(plotdir,frbdat.frbname,frbdat.subband[0],frbdat.subband[1],frbdat.rm), transparent=True, format='pdf')
	plt.show()

	return(0)

#	----------------------------------------------------------------------------------------------------------

def plot_iquvtsep(plotdir,frbdat,fmhz0,xlim,ylim,fsize,ttms):

	#	Plot IQUV profiles and dynspecs separately
	
	avgchan		=	int(round(frbdat.nchan/frbdat.ffac))
	
	chanwmhz	=	NatBWmhz/avgchan
	fmhzarr		=	np.arange(fmhz0-(NatBWmhz/2)+(chanwmhz/2), fmhz0+NatBWmhz/2, chanwmhz)
	
	fig 	= plt.figure(figsize=(fsize[0], fsize[1]))
	ax 		= fig.add_axes([0.07, 0.50, 0.23,0.48])
	ax.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
	
	ax.axhline(c='c',ls='--')
	ax.plot(frbdat.tmsarr, frbdat.it/np.nanmax(frbdat.it), 'k-', label='I')
	ax.set_ylim(ylim)
	ax.set_xlim(xlim)
	ax.legend(loc='upper right')
	ax.set_ylabel(r'Normalized flux density')
	ax.set_xticklabels([])
	ax.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
	ax.yaxis.set_label_coords(-0.2, 0.5)
	
	ax5 		= fig.add_axes([0.30, 0.50, 0.23,0.48])
	ax5.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
	
	ax5.axhline(c='c',ls='--')
	ax5.plot(frbdat.tmsarr, frbdat.qt/np.nanmax(frbdat.it), 'r-', label='Q')
	ax5.set_ylim(ylim)
	ax5.set_xlim(xlim)
	ax5.legend(loc='upper right')
	ax5.set_xticklabels([])
	ax5.set_yticklabels([])
	ax5.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax5.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
	
	ax6 		= fig.add_axes([0.53, 0.50, 0.23,0.48])
	ax6.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
	
	ax6.axhline(c='c',ls='--')
	ax6.plot(frbdat.tmsarr, frbdat.ut/np.nanmax(frbdat.it), 'm-', label='U')
	ax6.set_ylim(ylim)
	ax6.set_xlim(xlim)
	ax6.legend(loc='upper right')
	ax6.set_xticklabels([])
	ax6.set_yticklabels([])
	ax6.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax6.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
	
	ax8 		= fig.add_axes([0.76, 0.50, 0.23,0.48])
	ax8.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
	
	ax8.axhline(c='c',ls='--')
	ax8.plot(frbdat.tmsarr, frbdat.vt/np.nanmax(frbdat.it), 'b-', label='V')
	ax8.set_ylim(ylim)
	ax8.set_xlim(xlim)
	ax8.legend(loc='upper right')
	ax8.set_xticklabels([])
	ax8.set_yticklabels([])
	ax8.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax8.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
	
	ax0 		= fig.add_axes([0.07, 0.10, 0.23,0.40])
	ax0.imshow(frbdat.ids, aspect='auto',cmap='seismic',interpolation="none",vmin=-np.nanmax(np.abs(frbdat.ids)),vmax=np.nanmax(np.abs(frbdat.ids)), \
		extent=([xlim[0],xlim[1],(fmhzarr[-1]-chanwmhz)/1.0e3, (fmhzarr[0]+chanwmhz)/1.0e3]))
	ax0.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax0.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
	ax0.set_ylabel(r'$\nu$ (GHz)')
	ax0.set_xlabel(r'Time (ms)')
	ax0.yaxis.set_label_coords(-0.2, 0.5)
	
	ax1 		= fig.add_axes([0.30, 0.10, 0.23,0.40])
	ax1.imshow(frbdat.qds, aspect='auto',cmap='seismic',interpolation="none",vmin=-np.nanmax(np.abs(frbdat.qds)),vmax=np.nanmax(np.abs(frbdat.qds)), \
		extent=([xlim[0],xlim[1],(fmhzarr[-1]-chanwmhz)/1.0e3, (fmhzarr[0]+chanwmhz)/1.0e3]))
	ax1.set_yticklabels([])
	ax1.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
	ax1.set_xlabel(r'Time (ms)')
	
	ax2 		= fig.add_axes([0.53, 0.10, 0.23,0.40])
	ax2.imshow(frbdat.uds, aspect='auto',cmap='seismic',interpolation="none",vmin=-np.nanmax(np.abs(frbdat.uds)),vmax=np.nanmax(np.abs(frbdat.uds)), \
		extent=([xlim[0],xlim[1],(fmhzarr[-1]-chanwmhz)/1.0e3, (fmhzarr[0]+chanwmhz)/1.0e3]))
	ax2.set_yticklabels([])
	ax2.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
	ax2.set_xlabel(r'Time (ms)')
	
	ax3 		= fig.add_axes([0.76, 0.10, 0.23,0.40])
	ax3.imshow(frbdat.vds, aspect='auto',cmap='seismic',interpolation="none",vmin=-np.nanmax(np.abs(frbdat.vds)),vmax=np.nanmax(np.abs(frbdat.vds)), \
		extent=([xlim[0],xlim[1],(fmhzarr[-1]-chanwmhz)/1.0e3, (fmhzarr[0]+chanwmhz)/1.0e3]))
	ax3.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax3.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
	ax3.set_yticklabels([])
	ax3.set_xlabel(r'Time (ms)')
		
	plt.savefig("{}{}_stks_sb_{}_of_{}_rm_{:.2f}.pdf".format(plotdir,frbdat.frbname,frbdat.subband[0],frbdat.subband[1],frbdat.rm))
	plt.show()

	return(0)

#	----------------------------------------------------------------------------------------------------------

def plot_ilvt(plotdir,frbdat,bfrbdat,xlim,fsize,ttms):

	#	Plot ILV profiles
	
	phits		=	frbdat.chit
	dphits		=	frbdat.echit
	badi		=	np.where((frbdat.it/frbdat.irms) < 3.0)
	badl		=	np.union1d(np.where((frbdat.lt/frbdat.elts) < 3.0),badi)
	nongoodphi	=	np.union1d(np.where(dphits > 5.0),badl)
	phits[nongoodphi]	=	np.nan
	dphits[nongoodphi]	=	np.nan
		
	bbadi		=	np.where((bfrbdat.it/bfrbdat.irms) < 3.0)
	bbadl		=	np.union1d(np.where((bfrbdat.lt/bfrbdat.elts) < 3.0),bbadi)
	blts		=	bfrbdat.lt
	blts[bbadl]	=	np.nan
	
	fig 	= plt.figure(figsize=(fsize[0], fsize[1]))
	ax 		= fig.add_axes([0.12, 0.10, 0.85,0.85])
	ax.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
	
	ax.axhline(c='c',ls='--')
	ax.plot(bfrbdat.tmsarr, bfrbdat.it/np.nanmax(bfrbdat.it), 'k-', label='I')
	ax.plot(bfrbdat.tmsarr, blts/np.nanmax(bfrbdat.it), 'r-', label='L',lw=1.0)
	ax.plot(bfrbdat.tmsarr, bfrbdat.vt/np.nanmax(bfrbdat.it), 'b-', label='V',lw=1.0)
	ax.set_ylim(ymax=1.1)
	ax.set_xlim(xlim)
	ax.legend(loc='upper right')
	ax.set_ylabel(r'Normalized flux density')
	ax.set_xlabel(r'Time (ms)')
	ax.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
	ax.yaxis.set_label_coords(-0.09, 0.5)
		
	plt.savefig("{}{}_ilv_sb_{}_of_{}_rm_{:.2f}.pdf".format(plotdir,frbdat.frbname,frbdat.subband[0],frbdat.subband[1],frbdat.rm))
	plt.show()

	return(0)

#	----------------------------------------------------------------------------------------------------------

def plot_ilvpat(plotdir,frbdat,bfrbdat,xlim,fsize,ttms,ttpa):

	#	Plot ILV and PA profiles
	
	phits		=	frbdat.chit
	dphits		=	frbdat.echit
	badi		=	np.where((frbdat.it/frbdat.irms) < 3.0)
	badl		=	np.union1d(np.where((frbdat.lt/frbdat.elts) < 3.0),badi)
	nongoodphi	=	np.union1d(np.where(dphits > 5.0),badl)
	phits[nongoodphi]	=	np.nan
	dphits[nongoodphi]	=	np.nan
		
	bbadi		=	np.where((bfrbdat.it/bfrbdat.irms) < 3.0)
	bbadl		=	np.union1d(np.where((bfrbdat.lt/bfrbdat.elts) < 3.0),bbadi)
	blts		=	bfrbdat.lt
	blts[bbadl]	=	np.nan
	
	fig 	= plt.figure(figsize=(fsize[0], fsize[1]))
	ax 		= fig.add_axes([0.12, 0.45, 0.85,0.53])
	ax.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
	
	ax.axhline(c='c',ls='--')
	ax.plot(bfrbdat.tmsarr, bfrbdat.it/np.nanmax(bfrbdat.it), 'k-', label='I')
	ax.plot(bfrbdat.tmsarr, blts/np.nanmax(bfrbdat.it), 'r-', label='L',lw=1.0)
	ax.plot(bfrbdat.tmsarr, bfrbdat.vt/np.nanmax(bfrbdat.it), 'b-', label='V',lw=1.0)
	#ax.set_ylim(ymax=0.45, ymin=-0.15)
	ax.set_xlim(xlim)
	ax.legend(loc='upper right')
	ax.set_ylabel(r'Normalized flux density')
	ax.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
	ax.set_xticklabels([])
	ax.yaxis.set_label_coords(-0.09, 0.5)
		
	ax1 	= fig.add_axes([0.12, 0.10, 0.85,0.35])
	ax1.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
	ax1.errorbar(frbdat.tmsarr, phits, dphits, fmt='bo', markersize=4, lw=0.5, capsize=2)
	
	'''fitind		=	np.where(np.isfinite(phits))	
	popt,pcov	=	curve_fit(chi_nonrel_rvm_ew01, frbdat.tmsarr[fitind], phits[fitind], sigma=dphits[fitind], absolute_sigma=True, p0=(0.0, 2.0*0.8, 70, 100.0, 120.0))
	perr		=	np.sqrt(np.diag(pcov))
	print(popt)
	print(perr)	
	ax1.plot(frbdat.tmsarr, chi_nonrel_rvm_ew01(frbdat.tmsarr, *popt), c='k', ls='--', lw=2)'''
	
	ax1.set_xlim(xlim)
	#ax1.set_ylim(ymax=-45)
	ax1.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax1.yaxis.set_major_locator(ticker.MultipleLocator(ttpa))
	ax1.set_xlabel(r'Time (ms)')
	ax1.set_ylabel(r'PA (deg)')
	ax1.yaxis.set_label_coords(-0.09, 0.5)
	
	plt.savefig("{}{}_ilvpa_sb_{}_of_{}_rm_{:.2f}.pdf".format(plotdir,frbdat.frbname,frbdat.subband[0],frbdat.subband[1],frbdat.rm))
	plt.show()

	return(0)

#	----------------------------------------------------------------------------------------------------------

def plot_ilvparm(plotdir,rmfile,frbdat,bfrbdat,xlim,fsize,ttms,ttrm,ttpa):

	#	Plot ILV, PA & RM profiles
	
	phits		=	frbdat.chit
	dphits		=	frbdat.echit
	badi		=	np.where((frbdat.it/frbdat.irms) < 3.0)
	badl		=	np.union1d(np.where((frbdat.lt/frbdat.elts) < 3.0),badi)
	nongoodphi	=	np.union1d(np.where(dphits > 5.0),badl)
	phits[nongoodphi]	=	np.nan
	dphits[nongoodphi]	=	np.nan
	
	bbadi		=	np.where((bfrbdat.it/bfrbdat.irms) < 3.0)
	bbadl		=	np.union1d(np.where((bfrbdat.lt/bfrbdat.elts) < 3.0),bbadi)
	blts		=	bfrbdat.lt
	blts[bbadl]	=	np.nan
	
	fig 	= plt.figure(figsize=(fsize[0]/intocm, fsize[1]/intocm))
	ax 		= fig.add_axes([0.12, 0.58, 0.87,0.40])
	ax.text(0.02,0.85,'A',fontsize=10,fontweight='bold',transform=ax.transAxes)
	ax.tick_params(axis="both",direction="out",bottom=True,right=False,top=False,left=True)
	
	ax.axhline(c='c',ls='--',lw=0.5)
	ax.plot(bfrbdat.tmsarr, bfrbdat.it/np.nanmax(bfrbdat.it), 'k-', label='I')
	ax.plot(bfrbdat.tmsarr, blts/np.nanmax(bfrbdat.it), 'r-', label='L',lw=1.0)
	ax.plot(bfrbdat.tmsarr, bfrbdat.vt/np.nanmax(bfrbdat.it), 'b-', label='V',lw=1.0)
	ax.set_ylim(ymax=1.1)
	ax.set_xlim(xlim)
	ax.legend(loc='upper right')
	ax.set_ylabel(r'Normalized flux density')
	ax.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
	ax.set_xticklabels([])
	ax.yaxis.set_label_coords(-0.09, 0.5)
		
	ax1 	= fig.add_axes([0.12, 0.30, 0.87,0.28])
	ax1.text(0.02,0.8,'B',fontsize=10,fontweight='bold',transform=ax1.transAxes)
	ax1.tick_params(axis="both",direction="out",bottom=True,right=False,top=False,left=True)
	ax1.errorbar(frbdat.tmsarr, phits, dphits, fmt='ro', markersize=2, lw=0.5, capsize=2)
	#ax1.plot(frbdat.tmsarr, chi_nonrel_rvm_ew01(frbdat.tmsarr, 0.026, 2.0*1.238, -38.5, 122.3, 107.8), c='k', ls='--', lw=2, zorder=2.5)
	#ax1.set_ylim([-64,-16])
	#ax1.set_ylim([16,64])
	ax1.set_xlim(xlim)
	ax1.set_xticklabels([])
	ax1.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax1.yaxis.set_major_locator(ticker.MultipleLocator(ttpa))
	ax1.set_ylabel(r'PA (deg)')
	ax1.yaxis.set_label_coords(-0.09, 0.5)
	
	rmt		=	np.loadtxt(rmfile)
	tstarts	=	rmt[1:,0]
	tends	=	rmt[1:,1]
	tcent	=	(tstarts+tends)/2.0
	etl		=	tcent-tstarts
	etr		=	tends-tcent
	rmv		=	rmt[1:,6]	#	rmt[1:,2]
	ermv	=	rmt[1:,7]	#	rmt[1:,3]
	
	ax2 	= fig.add_axes([0.12, 0.10, 0.87,0.20])
	ax2.text(0.02,0.76,'C',fontsize=10,fontweight='bold',transform=ax2.transAxes)
	ax2.tick_params(axis="both",direction="out",bottom=True,right=False,top=False,left=True)
	#ax2.axvspan(-0.0165,0.045, fc='lightgrey',alpha=0.5)
	#ax2.axvspan(1.235,1.355, fc='lightgrey',alpha=0.5)
	ax2.errorbar(tcent, rmv, yerr=ermv, xerr=(etl, etr), fmt='bo', lw= 0.5, markersize=2, capsize=2)
	ax2.set_xlim(xlim)
	ax2.set_ylim([7,31])
	#ax2.set_ylim([-24,24])
	ax2.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax2.yaxis.set_major_locator(ticker.MultipleLocator(ttrm))
	ax2.set_xlabel(r'Time (ms)')
	ax2.set_ylabel(r'RM (rad m$^{-2}$)')
	ax2.yaxis.set_label_coords(-0.09, 0.5)
	
	plt.savefig("{}{}_ilvparm_sb_{}_of_{}_rm_{:.2f}.pdf".format(plotdir,frbdat.frbname,frbdat.subband[0],frbdat.subband[1],frbdat.rm), transparent=True, format='pdf')
	plt.savefig("{}{}_ilvparm_sb_{}_of_{}_rm_{:.2f}.eps".format(plotdir,frbdat.frbname,frbdat.subband[0],frbdat.subband[1],frbdat.rm), transparent=True, format='eps')
	plt.show()

	return(0)

#	----------------------------------------------------------------------------------------------------------

def plot_fracparm(plotdir,rmfile,frbdat,bfrbdat,xlim,fsize,ttms,ttrm):

	#	Plot polarization fraction, PA & RM profiles
	
	fig 	= plt.figure(figsize=(fsize[0], fsize[1]))
	ax 		= fig.add_axes([0.12, 0.60, 0.85,0.38])
	ax.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
	
	pf		=	frbdat.pfract
	epf		=	frbdat.epfract
	lf		=	frbdat.lfract
	elf		=	frbdat.elfract
	vf		=	frbdat.vfract
	evf		=	frbdat.evfract
	badi	=	np.where((frbdat.it/frbdat.irms) < 5.0)
	pf[badi]	=	np.nan
	lf[badi]	=	np.nan
	vf[badi]	=	np.nan
	epf[badi]	=	np.nan
	elf[badi]	=	np.nan
	evf[badi]	=	np.nan
	
	ax.axhline(c='c',ls='--')
	ax.fill_between(frbdat.tmsarr,pf-epf,pf+epf,color='lightgray',alpha=0.5)
	ax.fill_between(frbdat.tmsarr,lf-elf,lf+elf,color='coral',alpha=0.6)
	ax.fill_between(frbdat.tmsarr,vf-evf,vf+evf,color='aqua',alpha=0.6)
	ax.plot(frbdat.tmsarr, pf, 'md',lw=2, fillstyle='none', label='P/I')
	ax.plot(frbdat.tmsarr, lf, 'ro',lw=2, label='L/I')
	ax.plot(frbdat.tmsarr, vf, 'b*',lw=2, label='V/I')	
	ax.plot(bfrbdat.tmsarr, bfrbdat.it/np.nanmax(bfrbdat.it),'k-', lw=1.0)
	
	ax.set_ylim(ymax=1.3)
	ax.set_xlim(xlim)
	ax.legend(loc='upper right')
	ax.set_ylabel(r'Fraction')
	ax.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
	ax.set_xticklabels([])
	ax.yaxis.set_label_coords(-0.09, 0.5)
	
	phits		=	frbdat.chit
	dphits		=	frbdat.echit
	badi		=	np.where((frbdat.it/frbdat.irms) < 3.0)
	badl		=	np.union1d(np.where((frbdat.lt/frbdat.elts) < 3.0),badi)
	nongoodphi	=	np.union1d(np.where(dphits > 5.0),badl)
	phits[nongoodphi]	=	np.nan
	dphits[nongoodphi]	=	np.nan
	
	ax1 	= fig.add_axes([0.12, 0.33, 0.85,0.27])
	ax1.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
	ax1.fill_between(frbdat.tmsarr,phits-dphits, phits+dphits,color='coral',alpha=0.6)
	ax1.plot(frbdat.tmsarr, phits,'r*', markersize=6)
	ax1.set_xlim(xlim)
	ax1.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax1.set_ylabel(r'PA (deg)')
	ax1.yaxis.set_label_coords(-0.09, 0.5)
	
	rmt		=	np.loadtxt(rmfile)
	tstarts	=	rmt[1:,0]
	tends	=	rmt[1:,1]
	tcent	=	(tstarts+tends)/2.0
	etl		=	tcent-tstarts
	etr		=	tends-tcent
	rmv		=	rmt[1:,6]	#	rmt[1:,2]
	ermv	=	rmt[1:,7]	#	rmt[1:,3]
	
	ax2 	= fig.add_axes([0.12, 0.10, 0.85,0.23])
	ax2.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
	ax2.errorbar(tcent, rmv, yerr=ermv, xerr=(etl, etr), fmt='md', lw= 0.5, markersize=5, capsize=2)
	ax2.set_xlim(xlim)
	ax2.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax2.yaxis.set_major_locator(ticker.MultipleLocator(ttrm))
	ax2.set_xlabel(r'Time (ms)')
	ax2.set_ylabel(r'RM (rad m$^{-2}$)')
	ax2.yaxis.set_label_coords(-0.09, 0.5)
	
	plt.savefig("{}{}_fracparm_sb_{}_of_{}_rm_{:.2f}.pdf".format(plotdir,frbdat.frbname,frbdat.subband[0],frbdat.subband[1],frbdat.rm))
	plt.show()

	return(0)

#	----------------------------------------------------------------------------------------------------------

def plot_fracpa(plotdir,rmfile,frbdat,bfrbdat,xlim,fsize,ttms,ttrm):

	#	Plot polarization fraction & PA profiles
	
	fig 	= plt.figure(figsize=(fsize[0], fsize[1]))
	ax 		= fig.add_axes([0.12, 0.45, 0.85,0.53])
	ax.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
	
	pf		=	frbdat.pfract
	epf		=	frbdat.epfract
	lf		=	frbdat.lfract
	elf		=	frbdat.elfract
	vf		=	frbdat.vfract
	evf		=	frbdat.evfract
	badi	=	np.where((frbdat.it/frbdat.irms) < 5.0)
	pf[badi]	=	np.nan
	lf[badi]	=	np.nan
	vf[badi]	=	np.nan
	epf[badi]	=	np.nan
	elf[badi]	=	np.nan
	evf[badi]	=	np.nan
	
	ax.axhline(c='c',ls='--')
	ax.fill_between(frbdat.tmsarr,pf-epf,pf+epf,color='lightgray',alpha=0.5)
	ax.fill_between(frbdat.tmsarr,lf-elf,lf+elf,color='coral',alpha=0.6)
	ax.fill_between(frbdat.tmsarr,vf-evf,vf+evf,color='aqua',alpha=0.6)
	ax.plot(frbdat.tmsarr, pf, 'md',lw=2, fillstyle='none', label='P/I')
	ax.plot(frbdat.tmsarr, lf, 'ro',lw=2, label='L/I')
	ax.plot(frbdat.tmsarr, vf, 'b*',lw=2, label='V/I')	
	ax.plot(bfrbdat.tmsarr, bfrbdat.it/np.nanmax(bfrbdat.it),'k-', lw=1.0)
	
	ax.set_ylim(ymax=1.3)
	ax.set_xlim(xlim)
	ax.legend(loc='upper right')
	ax.set_ylabel(r'Fraction')
	ax.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
	ax.set_xticklabels([])
	ax.yaxis.set_label_coords(-0.09, 0.5)
	
	phits		=	frbdat.chit
	dphits		=	frbdat.echit
	badi		=	np.where((frbdat.it/frbdat.irms) < 3.0)
	badl		=	np.union1d(np.where((frbdat.lt/frbdat.elts) < 3.0),badi)
	nongoodphi	=	np.union1d(np.where(dphits > 5.0),badl)
	phits[nongoodphi]	=	np.nan
	dphits[nongoodphi]	=	np.nan
	
	ax1 	= fig.add_axes([0.12, 0.10, 0.85,0.35])
	ax1.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
	ax1.fill_between(frbdat.tmsarr,phits-dphits, phits+dphits,color='coral',alpha=0.6)
	ax1.plot(frbdat.tmsarr, phits,'r*', markersize=6)
	ax1.set_xlim(xlim)
	ax1.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax1.set_xlabel(r'Time (ms)')
	ax1.set_ylabel(r'PA (deg)')
	ax1.yaxis.set_label_coords(-0.09, 0.5)
	
	plt.savefig("{}{}_fracpa_sb_{}_of_{}_rm_{:.2f}.pdf".format(plotdir,frbdat.frbname,frbdat.subband[0],frbdat.subband[1],frbdat.rm))
	plt.show()

	return(0)

#	----------------------------------------------------------------------------------------------------------

def congaussfitter(plotdir,bfrbdat,xlim,fsize,ttms,ng,pol):

	#	Fit profile with n Gaussians convolved with the same exponential
	
	xdat0	=	bfrbdat.tmsarr
	relind	=	np.where((xdat0 - xlim[0])*(xdat0 - xlim[1]) < 0.0)[0]
		
	if(pol=='L'):
		fdat0	=	bfrbdat.lt/np.nanmax(bfrbdat.it)
	else:
		fdat0	=	bfrbdat.it/np.nanmax(bfrbdat.it)
	
	xdat		=	xdat0[relind]
	fdat		=	fdat0[relind]
	inoise		=	bfrbdat.irms/np.nanmax(bfrbdat.it)
	
	nfitfn		=	lambda x, t, *pars : ncompexp(x, ng, t, pars)
	cnfitfn		=	lambda x, t, pars, f0 : ncompexp(x, ng, t, pars) - f0
	infitfn		=	lambda x, t, pars : -ncompexp(x, ng, t, pars)
	
	iniguess	=	np.ones(1+3*ng, dtype=float)
	
	popt,pcov	=	curve_fit(nfitfn, xdat, fdat, p0 = iniguess, bounds = (-2.0, 2.0), maxfev=50000)
	perr		=	np.sqrt(np.diag(pcov))
	
	adr2		=	r2ad(fdat, nfitfn(xdat, *popt), 1+3*ng)
	
	print("Fit parameters \n")
	print(popt)
	print(perr)
	print(" 1 - Adjusted R2 = %.5f"%(1.0 - adr2))
	
	#	Calculation of FWHM
	
	tpeak		=	minimize_scalar(infitfn, bracket=(xlim[0],xlim[1]), args=(popt[0],popt[1:]), tol = bfrbdat.tresms/10.0)
	lhalfpt		=	root_scalar(cnfitfn, args=(popt[0],popt[1:],nfitfn(tpeak.x, *popt)/2.0), bracket=(xlim[0],tpeak.x), xtol = bfrbdat.tresms/10.0)
	rhalfpt		=	root_scalar(cnfitfn, args=(popt[0],popt[1:],nfitfn(tpeak.x, *popt)/2.0), bracket=(tpeak.x,xlim[1]), xtol = bfrbdat.tresms/10.0)
	fitfwhm		=	(rhalfpt.root - lhalfpt.root)	
	
	print("Peak at %.6f ms"%(tpeak.x))
		
	fwhmarr		=	[]
	for kk in range(0,1):
		jklen		=	int(round(0.8*len(xdat)))
		indcs		=	np.random.choice(len(xdat),jklen,replace=False)
		popt0,pcov0	=	curve_fit(nfitfn, xdat[indcs], fdat[indcs], p0 = iniguess, bounds = (-2.0, 2.0), maxfev=100000)
		tpeak0		=	minimize_scalar(infitfn, bracket=(xlim[0],xlim[1]), args=(popt0[0],popt0[1:]), tol = bfrbdat.tresms/10.0)
		lhalfpt0	=	root_scalar(cnfitfn, args=(popt0[0],popt0[1:],nfitfn(tpeak0.x, *popt0)/2.0), bracket=(xlim[0],tpeak0.x), xtol = bfrbdat.tresms/10.0)
		rhalfpt0	=	root_scalar(cnfitfn, args=(popt0[0],popt0[1:],nfitfn(tpeak0.x, *popt0)/2.0), bracket=(tpeak0.x,xlim[1]), xtol = bfrbdat.tresms/10.0)
		fitfwhm0	=	(rhalfpt0.root - lhalfpt0.root)
		
		print("%d / %d FWHM = %.6f"%(kk,1234,fitfwhm0))
		fwhmarr.append(fitfwhm0)
	
	sdfwhm			=	1.48*np.nanmedian(np.abs(np.array(fwhmarr) - np.nanmedian(np.array(fwhmarr))))
	
	print("FWHM = %.6f +/- %.6f ms"%(fitfwhm,sdfwhm))
	#plt.hist(fwhmarr)
	#plt.show()
	
	#	Done now
	
	fig 	= plt.figure(figsize=(fsize[0]/intocm, fsize[1]/intocm))
	ax 		= fig.add_axes([0.12, 0.10, 0.85,0.63])
	ax.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
	ax.text(0.65, 0.9, r"$\tau = %d \pm %d \: \mu$s"%(int(round(popt[0]*1.0e3)),int(round(perr[0]*1.0e3))),transform=ax.transAxes)
	ax.text(0.65, 0.8, r"$1 - R^2 = %.5f $"%(1.0 - adr2),transform=ax.transAxes)
	ax.text(0.04,0.9,'b',fontsize=8,fontweight='bold',transform=ax.transAxes)
	ax.axhline(c='c',ls='--',lw=0.25)
	ax.plot(xdat, fdat, 'r-',lw=0.5)
	ax.plot(xdat, nfitfn(xdat, *popt), 'k-',lw=0.5)
	
	for i in range(0,ng):
		ax.plot(xdat, expg(xdat, popt[0], popt[1+3*i], popt[2+3*i], popt[3+3*i]), 'b:',lw=0.5)
	
	#ax.set_ylim(ymax=0.19)
	ax.set_xlim(xlim)
	#ax.legend(loc='upper right')
	ax.set_ylabel(r'Normalized flux density')
	ax.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
	ax.set_xlabel(r'Time (ms)')
	ax.yaxis.set_label_coords(-0.09, 0.5)
	
	ax1 	= fig.add_axes([0.12, 0.73, 0.85,0.25])
	ax1.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
	ax1.axhline(c='c',ls='--',lw=0.25)
	ax1.plot(xdat, (fdat - nfitfn(xdat, *popt))/inoise, 'r-', lw=0.5)
	ax1.text(0.04,0.8,'a',fontsize=8,fontweight='bold',transform=ax1.transAxes)
	ax1.set_xlim(xlim)
	ax1.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax1.set_xticklabels([])
	ax1.set_ylabel(r'Deviation / $\sigma$')
	ax1.yaxis.set_label_coords(-0.09, 0.5)
		
	plt.savefig("{}{}_pcomp_sb_{}_of_{}_rm_{:.2f}.eps".format(plotdir,bfrbdat.frbname,bfrbdat.subband[0],bfrbdat.subband[1],bfrbdat.rm), transparent=True, format='eps')
	plt.savefig("{}{}_pcomp_sb_{}_of_{}_rm_{:.2f}.pdf".format(plotdir,bfrbdat.frbname,bfrbdat.subband[0],bfrbdat.subband[1],bfrbdat.rm), transparent=True, format='pdf')
	plt.show()

	return(0)

#	----------------------------------------------------------------------------------------------------------

def plot_irm(plotdir,rmfile,frbdat,bfrbdat,xlim,fsize,ttms,ttrm, baserm):

	#	Plot I & RM profiles
	
	rmt		=	np.loadtxt(rmfile)
	tstarts	=	rmt[1:,0]
	tends	=	rmt[1:,1]
	tcent	=	(tstarts+tends)/2.0
	etl		=	tcent-tstarts
	etr		=	tends-tcent
	rmv		=	rmt[1:,6]	#	rmt[1:,2]
	ermv	=	rmt[1:,7]	#	rmt[1:,3]
	'''
	popt,pcov	=	curve_fit(gsn, tcent, rmv-baserm, sigma=ermv, absolute_sigma=True, p0=(0.0, tcent[np.argmax(np.abs(rmv-baserm))], 1.0))
	perr	=	np.sqrt(np.diag(pcov))
	print(popt)
	print(perr)
	'''
	fig 	= 	plt.figure(figsize=(fsize[0]/intocm, fsize[1]/intocm))
	ax 		= 	fig.add_axes([0.12, 0.12, 0.86,0.86])
	ax.text(0.04,0.9,'a',fontsize=8,fontweight='bold',transform=ax.transAxes)
	ax2		=	ax.twinx()	
	ax2.axhline(c='c',ls='--',lw=0.25)
	ax2.plot(bfrbdat.tmsarr, bfrbdat.it/np.nanmax(bfrbdat.it), 'b-',lw=0.5)
	ax2.set_ylim(ymax=1.1)
	ax2.set_yticks([])
		
	ax.errorbar(tcent, rmv, yerr=ermv, xerr=(etl, etr), fmt='ro', lw= 0.5, markersize=5, capsize=2)	
	
	#ax.plot(bfrbdat.tmsarr, baserm+gsn(bfrbdat.tmsarr, -9.464, 0.0106, 0.025), 'k--', lw=1.0)
	
	ax.plot(bfrbdat.tmsarr, baserm+gsn(bfrbdat.tmsarr, 17.97, 0.040, 0.047)+gsn(bfrbdat.tmsarr, -15.15, 1.290, 0.063), 'k--')
	
	ax.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax.yaxis.set_major_locator(ticker.MultipleLocator(ttrm))
	ax.set_xlim(xlim)
	#ax.set_ylim(ymax=6.0)
	ax.set_xlabel(r'Time (ms)')
	ax.set_ylabel(r'RM (rad m$^{-2}$)')
	ax.yaxis.set_label_coords(-0.08, 0.5)
	
	plt.savefig("{}{}_irm_sb_{}_of_{}_rm_{:.2f}.pdf".format(plotdir,frbdat.frbname,frbdat.subband[0],frbdat.subband[1],frbdat.rm), transparent=True, format='pdf')
	plt.savefig("{}{}_irm_sb_{}_of_{}_rm_{:.2f}.eps".format(plotdir,frbdat.frbname,frbdat.subband[0],frbdat.subband[1],frbdat.rm), transparent=True, format='eps')
	plt.show()

	return(0)

#	----------------------------------------------------------------------------------------------------------

def plot_dpa(plotdir,rmfile,frbdat,bfrbdat,xlim,fsize,ttms,ntp):

	#	Plot PA profile and dPA/dt
	
	print("Calculating slope from %d points"%(2*ntp+1))
	
	phits		=	frbdat.chit
	dphits		=	frbdat.echit
	badi		=	np.where((frbdat.it/frbdat.irms) < 3.0)
	badl		=	np.union1d(np.where((frbdat.lt/frbdat.elts) < 3.0),badi)
	nongoodphi	=	np.union1d(np.where(dphits > 5.0),badl)
	phits[nongoodphi]	=	np.nan
	dphits[nongoodphi]	=	np.nan
	
	dpadt			=	np.zeros(phits.shape, dtype=float)
	edpadt			=	np.zeros(phits.shape, dtype=float)	
	dpadt[:ntp]		=	np.nan
	edpadt[:ntp]	=	np.nan
	dpadt[-ntp:]	=	np.nan
	edpadt[-ntp:]	=	np.nan
		
	for ti in range(ntp,len(phits)-ntp):
		phi3	=	phits[ti-ntp:ti+ntp+1]
		dphi3	=	dphits[ti-ntp:ti+ntp+1]
		tarr3	=	frbdat.tmsarr[ti-ntp:ti+ntp+1]
		
		if(np.count_nonzero(np.isfinite(phi3))==(2*ntp+1)):
			popt,pcov	=	np.polyfit(tarr3, phi3, deg=1, w=1.0/dphi3, cov=True)
			perr		=	np.sqrt(np.diag(pcov))
			dpadt[ti]	=	popt[0]
			edpadt[ti]	=	perr[0]
		else:
			dpadt[ti]	=	np.nan
			edpadt[ti]	=	np.nan
	
	grat		=	np.isfinite(dpadt)	
	popt,pcov	=	curve_fit(gsn, frbdat.tmsarr[grat], dpadt[grat]/1.0e3, sigma=edpadt[grat]/1.0e3, absolute_sigma=True, p0=(-1.0, 0.0, 0.1))
	perr		=	np.sqrt(np.diag(pcov))
	print(popt)
	print(perr)
		
	fig 	= plt.figure(figsize=(fsize[0]/intocm, fsize[1]/intocm))
	ax 		= 	fig.add_axes([0.15, 0.48, 0.83,0.50])
	ax.text(0.04,0.9,'e',fontsize=8,fontweight='bold',transform=ax.transAxes)
	ax2		=	ax.twinx()	
	ax2.axhline(c='c',ls='--',lw=0.25)
	ax2.plot(bfrbdat.tmsarr, bfrbdat.it/np.nanmax(bfrbdat.it), 'c-', lw=0.5)
	ax2.set_ylim([-0.1, 1.1])
	ax2.set_yticks([])
	
	ax.errorbar(frbdat.tmsarr, phits,dphits, fmt='b*', markersize=5, lw = 0.5, capsize=2)
	
	ax.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax.set_xticklabels([])
	ax.set_xlim(xlim)
	ax.set_ylabel(r'PA (deg)')
	ax.yaxis.set_label_coords(-0.11, 0.5)	
	
	ax1 		= 	fig.add_axes([0.15, 0.10, 0.83,0.38])
	ax1.text(0.04,0.84,'f',fontsize=8,fontweight='bold',transform=ax1.transAxes)
	ax1.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
	
	ax1.errorbar(frbdat.tmsarr, dpadt/1.0e3, edpadt/1.0e3, fmt='ro', markersize=3, lw = 0.5, capsize=2)
	ax1.plot(frbdat.tmsarr, gsn(frbdat.tmsarr, *popt), 'k--',lw=1.0)
	
	ax1.set_xlim(xlim)
	#ax1.set_ylim(ymax=0.5)
	ax1.xaxis.set_major_locator(ticker.MultipleLocator(ttms))
	ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
	ax1.set_xlabel(r'Time (ms)')
	ax1.set_ylabel(r'Rate (deg / $\mu$s)')
	ax1.yaxis.set_label_coords(-0.11, 0.5)
	
	plt.savefig("{}{}_dpa_sb_{}_of_{}_rm_{:.2f}.pdf".format(plotdir,frbdat.frbname,frbdat.subband[0],frbdat.subband[1],frbdat.rm), transparent=True, format='pdf')
	plt.savefig("{}{}_dpa_sb_{}_of_{}_rm_{:.2f}.eps".format(plotdir,frbdat.frbname,frbdat.subband[0],frbdat.subband[1],frbdat.rm), transparent=True, format='eps')
	plt.show()

	return(0)

#	----------------------------------------------------------------------------------------------------------










