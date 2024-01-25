#
#	Functions for FRB polarization analysis
#
#								AB, April 2023

#	Function list
#
#		twin_side_plot(frbname, dm, nchan, ffac, avgfac, tbasems, fmhz0, tpeakms, frbname0, dm0, nchan0, ffac0, avgfac0, tbasems0, fmhz00, tpeakms0):
#			Plot two FRBs side by side
#
#		def twin_same_plot(frbname, dm, nchan, ffac, avgfac, perms, tpeakms, frbname0, dm0, nchan0, ffac0, avgfac0, perms0, tpeakms0):
#			Plot two FRBs the same graph
#
#	--------------------------	Import modules	---------------------------

import os, sys
import numpy as np
from globalpars import *
import matplotlib.pyplot as plt
import matplotlib.colors as mpc
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker as ticker

mpl.rcParams['pdf.fonttype']	= 42
mpl.rcParams['ps.fonttype'] 	= 42
mpl.rcParams['savefig.dpi'] 	= 600
mpl.rcParams['font.family'] 	= 'sans-serif'
mpl.rcParams['font.size']		= 8

#	--------------------------	Analysis functions	-------------------------------

def twin_side_plot(frbname, dm, nchan, ffac, avgfac, tbasems, fmhz0, tpeakms, frbname0, dm0, nchan0, ffac0, avgfac0, tbasems0, fmhz00, tpeakms0):

	#	Plot two FRBs side by side
	
	datadir		=	htrdir+frbname+'/data/reduced/'
	plotdir		=	htrdir+frbname+'/plots/'	
	datadir0	=	htrdir+frbname0+'/data/reduced/'
	plotdir0	=	htrdir+frbname0+'/plots/'
	
	chanqfile	=	"{}{}_chanq_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)			
	avgchan		=	int(round(nchan/ffac))	
	chanqfile0	=	"{}{}_chanq_{}_{}_avg_{}_{}.npy".format(datadir0,frbname0,dm0,nchan0,ffac0,avgfac0)			
	avgchan0	=	int(round(nchan0/ffac0))	
	
	if(os.path.exists(chanqfile)):
		chanq	=	np.load(chanqfile)
	else:
		chanq	=	np.ones(avgchan, dtype=int)
		print("twin_side_plot ---- No bad channel list found for "+frbname)	
	if(os.path.exists(chanqfile0)):
		chanq0	=	np.load(chanqfile0)
	else:
		chanq0	=	np.ones(avgchan0, dtype=int)
		print("twin_side_plot ---- No bad channel list found for "+frbname0)
		
	itsfile		=	"{}{}_I_ts_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)	
	stksdsfile	=	"{}{}_stks_ds_zeromean_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)	
	noisefile	=	"{}{}_noise_spec_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)
	itsfile0	=	"{}{}_I_ts_{}_{}_avg_{}_{}.npy".format(datadir0,frbname0,dm0,nchan0,ffac0,avgfac0)	
	stksdsfile0	=	"{}{}_stks_ds_zeromean_{}_{}_avg_{}_{}.npy".format(datadir0,frbname0,dm0,nchan0,ffac0,avgfac0)	
	noisefile0	=	"{}{}_noise_spec_{}_{}_avg_{}_{}.npy".format(datadir0,frbname0,dm0,nchan0,ffac0,avgfac0)
	
	tresms		=	avgfac*nchan*Raw_time_res_ms
	tpeak		=	tpeakms/tresms
	tresms0		=	avgfac0*nchan0*Raw_time_res_ms
	tpeak0		=	tpeakms0/tresms0
	
	if((os.path.exists(itsfile)) and (os.path.exists(stksdsfile)) and (os.path.exists(itsfile0)) and (os.path.exists(stksdsfile0))):
		its		=	np.load(itsfile)
		stksds	=	np.load(stksdsfile)
		rmspec	=	np.load(noisefile)
		badcs	=	np.where(chanq==0)
		stksds[:,badcs]		=	np.nan
		
		its0	=	np.load(itsfile0)
		stksds0	=	np.load(stksdsfile0)
		rmspec0	=	np.load(noisefile0)
		badcs0	=	np.where(chanq0==0)
		stksds0[:,badcs0]	=	np.nan
				
		chanwmhz	=	NatBWmhz/avgchan
		fmhzarr		=	np.arange(fmhz0-(NatBWmhz/2)+(chanwmhz/2), fmhz0+NatBWmhz/2, chanwmhz)
		chanwmhz0	=	NatBWmhz/avgchan0
		fmhzarr00	=	np.arange(fmhz00-(NatBWmhz/2)+(chanwmhz0/2), fmhz00+NatBWmhz/2, chanwmhz0)
		
		reltms	=	np.arange(-tpeak, len(its) -tpeak - 0.5, 1.0, dtype=float)*tresms
		reltms0	=	np.arange(-tpeak0, len(its0) -tpeak0 - 0.5, 1.0, dtype=float)*tresms0
				
		fig		=	plt.figure(figsize=(5,8))			
		
		ax0		=	fig.add_subplot(2,1,1)		
		plt.axhline(c='c', ls='--')
		plt.plot(reltms, its,'k-',label='FRB '+frbname)
		#plt.step(reltms, its,where='mid',c='k',ls='-',label='FRB '+frbname)
		plt.ylim(ymax=1.2*np.nanmax(its[int(tpeak + tbasems[0]/tresms): int(tpeak + tbasems[1]/tresms)]), ymin=-0.2*np.nanmax(its[int(tpeak + tbasems[0]/tresms): int(tpeak + tbasems[1]/tresms)]))
		plt.xlim((tbasems[0], tbasems[1]))
		#plt.xlabel('Time - Peak (ms)')
		plt.ylabel('Flux density (arbitrary unit)')
		ax0.legend(loc='upper right')
		
		ax0		=	fig.add_subplot(2,1,2)		
		plt.axhline(c='c', ls='--')
		plt.plot(reltms0, its0,'k-',label='FRB '+frbname0)
		#plt.step(reltms0, its0,where='mid',c='k',ls='-',label='FRB '+frbname0)
		plt.ylim(ymax=1.2*np.nanmax(its0[int(tpeak0 + tbasems0[0]/tresms0): int(tpeak0 + tbasems0[1]/tresms0)]),ymin=-0.2*np.nanmax(its0[int(tpeak0 + tbasems0[0]/tresms0): int(tpeak0 + tbasems0[1]/tresms0)]))
		plt.xlim((tbasems0[0], tbasems0[1]))
		plt.xlabel('Time (ms)')
		plt.ylabel('Flux density (arbitrary unit)')
		ax0.legend(loc='upper right')
		
		plt.tight_layout()			
		plt.savefig("{}{}_twin_ts_{}_{}_avg_{}_{}.png".format(plotdir,frbname,dm,nchan,ffac,avgfac))
		plt.savefig("{}{}_twin_ts_{}_{}_avg_{}_{}.pdf".format(plotdir,frbname,dm,nchan,ffac,avgfac))
		plt.show()
		
		fig		=	plt.figure(figsize=(5,8))			
		
		ax0		=	fig.add_subplot(2,1,1)		
		plt.imshow(stksds[0], cmap='plasma', aspect='auto', interpolation='none', norm=mpc.SymLogNorm(vmin=0.0, linthresh= 3.0*np.nanmin(rmspec[0])))
		plt.xlim((tpeak + tbasems[0]/tresms, tpeak + tbasems[1]/tresms))
		#plt.xlabel('Time - Peak (ms)')
		plt.ylabel('Frequency channel')
		
		ax0		=	fig.add_subplot(2,1,2)		
		plt.imshow(stksds0[0], cmap='plasma', aspect='auto', interpolation='none', norm=mpc.SymLogNorm(vmin=0.0, linthresh= 3.0*np.nanmin(rmspec0[0])))
		plt.xlim((tpeak0 + tbasems0[0]/tresms0, tpeak0 + tbasems0[1]/tresms0))
		plt.xlabel('Time (ms)')

		plt.ylabel('Frequency channel')
		
		plt.tight_layout()			
		plt.savefig("{}{}_twin_ds_{}_{}_avg_{}_{}.png".format(plotdir,frbname,dm,nchan,ffac,avgfac))
		plt.savefig("{}{}_twin_ds_{}_{}_avg_{}_{}.pdf".format(plotdir,frbname,dm,nchan,ffac,avgfac))
		plt.show()
		
		del(stksds)
		del(its)
		del(stksds0)
		del(its0)
	else:
		print("twin_side_plot ---- PANIC - File not found!")

	return(tpeak)

#	-------------------------------------------------------------------------------

def twin_same_plot(frbname, dm, nchan, ffac, avgfac, perms, tpeakms, frbname0, dm0, nchan0, ffac0, avgfac0, perms0, tpeakms0):

	#	Plot two FRBs the same graph
	
	datadir		=	htrdir+frbname+'/data/reduced/'
	plotdir		=	htrdir+frbname+'/plots/'	
	datadir0	=	htrdir+frbname0+'/data/reduced/'
	plotdir0	=	htrdir+frbname0+'/plots/'
	
	chanqfile	=	"{}{}_chanq_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)			
	avgchan		=	int(round(nchan/ffac))	
	chanqfile0	=	"{}{}_chanq_{}_{}_avg_{}_{}.npy".format(datadir0,frbname0,dm0,nchan0,ffac0,avgfac0)			
	avgchan0	=	int(round(nchan0/ffac0))	
	
	if(os.path.exists(chanqfile)):
		chanq	=	np.load(chanqfile)
	else:
		chanq	=	np.ones(avgchan, dtype=int)
		print("twin_side_plot ---- No bad channel list found for "+frbname)	
	if(os.path.exists(chanqfile0)):
		chanq0	=	np.load(chanqfile0)
	else:
		chanq0	=	np.ones(avgchan0, dtype=int)
		print("twin_side_plot ---- No bad channel list found for "+frbname0)
		
	itsfile		=	"{}{}_I_ts_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)	
	stksdsfile	=	"{}{}_stks_ds_zeromean_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)	
	noisefile	=	"{}{}_noise_spec_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)
	itsfile0	=	"{}{}_I_ts_{}_{}_avg_{}_{}.npy".format(datadir0,frbname0,dm0,nchan0,ffac0,avgfac0)	
	stksdsfile0	=	"{}{}_stks_ds_zeromean_{}_{}_avg_{}_{}.npy".format(datadir0,frbname0,dm0,nchan0,ffac0,avgfac0)	
	noisefile0	=	"{}{}_noise_spec_{}_{}_avg_{}_{}.npy".format(datadir0,frbname0,dm0,nchan0,ffac0,avgfac0)
	
	tresms		=	avgfac*nchan*Raw_time_res_ms
	tpeak		=	tpeakms/tresms
	tresms0		=	avgfac0*nchan0*Raw_time_res_ms
	tpeak0		=	tpeakms0/tresms0
	
	if((os.path.exists(itsfile)) and (os.path.exists(stksdsfile)) and (os.path.exists(itsfile0)) and (os.path.exists(stksdsfile0))):
		its		=	np.load(itsfile)
		stksds	=	np.load(stksdsfile)
		rmspec	=	np.load(noisefile)
		badcs	=	np.where(chanq==0)
		stksds[:,badcs]		=	np.nan
		
		its0	=	np.load(itsfile0)
		stksds0	=	np.load(stksdsfile0)
		rmspec0	=	np.load(noisefile0)
		badcs0	=	np.where(chanq0==0)
		stksds0[:,badcs0]	=	np.nan
				
		reltms	=	np.arange(-tpeak, len(its) -tpeak - 0.5, 1.0, dtype=float)*tresms
		reltms0	=	np.arange(-tpeak0, len(its0) -tpeak0 - 0.5, 1.0, dtype=float)*tresms0
				
		fig		=	plt.figure(figsize=(12.1/intocm,9.0/intocm))				
		ax2		=	fig.add_subplot(1,1,1)		
		
		ax2.axhline(c='grey', ls='--',lw=0.5)
		ax2.plot((reltms/perms)+7.5e-5, its/np.nanmax(its),'b-', lw=1, label='FRB 20'+frbname+'A')
		ax2.plot((reltms0/perms0)+0.0024691, its0/np.nanmax(its0),'r-', lw=1, label='FRB 20'+frbname0+'A')
		#plt.step(reltms, its,where='mid',c='k',ls='-',label='FRB '+frbname)
		ax2.set_ylim([-0.05,1.05])
		ax2.set_xlim([-0.35, 1.65])
		ax2.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
		ax2.set_ylabel('Normalized flux density')
		ax2.set_xlabel(r'$\rm{t}$ / $\rm{\Delta T}$')
		#ax2.set_xlabel(r'$t$ (ms)')
		ax2.legend(loc='upper right')
				
		plt.tight_layout()			
		plt.savefig("{}{}_twins_{}_{}_avg_{}_{}.eps".format(plotdir,frbname,dm,nchan,ffac,avgfac), transparent=True, format='eps')
		plt.savefig("{}{}_twins_{}_{}_avg_{}_{}.pdf".format(plotdir,frbname,dm,nchan,ffac,avgfac), transparent=True, format='pdf')
		plt.show()
				
		del(stksds)
		del(its)
		del(stksds0)
		del(its0)
	else:
		print("twin_same_plot ---- PANIC - File not found!")

	return(tpeak)

#	-------------------------------------------------------------------------------































