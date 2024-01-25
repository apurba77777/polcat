#
#	Functions for FRB polarization analysis
#
#								AB, April 2023

#	Function list
#
#			locate_frb(frbname, dm):
#						Locates FRB in the Raw Stokes I time series at 1ms resolution
#
#			extract_zoom(frbname, dm, locms, halflenms):
#						Extract X & Y complex time series zooming around the FRB
#
#			voltage_ds(frbname, dm, nchan):
#						Construct voltage dynamic spectra for X and Y
#
#			gen_stokes_ds(frbname, dm, nchan, fmhz0):
#						Construct dynamic spectra for I, Q, U, V
#
#			favg_stokes(frbname, dm, nchan, ffac):
#						Average dynamic spectra for I, Q, U, V along frequency axis
#	
#			avg_stokes_ds(frbname, dm, nchan, ffac, avgfac):
#						Average dynamic spectra for I, Q, U, V along time axis
#
#			gen_noise_spec(frbname, dm, nchan, ffac, avgfac):
#						Subtract mean & generate noise spectra
#
#			gen_bad_chan_list(frbname, dm, nchan, ffac, avgfac):
#						Generate a list of bad/noisy channels
#
#			gen_I_ts(frbname, dm, nchan, ffac, avgfac):
#						Generate I time series

#	--------------------------	Import modules	---------------------------

import os, sys
import numpy as np
from globalpars import *
from htrfns.dynamic_fns import *
from htrfns.scintillation_fns import *
import scipy.stats as st

#	--------------------------	Analysis functions	-------------------------------

def locate_frb(frbname, dm):

	#	Locates FRB in the Raw Stokes I time series at 1ms resolution
	
	datadir		=	htrdir+frbname+'/data/celebi/'
	
	locms		=	-1
	snrms		=	-1.0
	idsfile		=	"{}{}_{}_1ms_I.npy".format(datadir,frbname,dm)			#	From CELEBI output
	
	if(os.path.exists(idsfile)):
		idspec		=	np.load(idsfile)
		if(idspec.shape[0]!=CelebiNchan):
			print("locate_frb ---- WARNING - Wrong data shape {}".format(idspec.shape))
		itseries	=	np.nanmean(idspec,axis=0)
		itrms		=	1.48*np.nanmedian(np.abs(itseries))
		locms		=	np.argmax(itseries)
		snrms		=	itseries[locms]/itrms	
		#print(itrms,itseries[locms])		
		del(idspec)		
	else:
		print("locate_frb ---- PANIC - Location undetermined! File {} not found!".format(idsfile))
	
	return(locms,snrms)

#	-------------------------------------------------------------------------------

def extract_zoom(frbname, dm, locms, halflenms):

	#	Extract X & Y complex time series zooming around the FRB
	
	datadir		=	htrdir+frbname+'/data/celebi/'
	
	ptimenat	=	int(round(float(locms)/Raw_time_res_ms))
	lenat		=	int(round(halflenms/Raw_time_res_ms))
	
	xvtfile		=	"{}{}_X_t_{}.npy".format(datadir,frbname,dm)			#	From CELEBI output
	yvtfile		=	"{}{}_Y_t_{}.npy".format(datadir,frbname,dm)			#	From CELEBI output
		
	if(os.path.exists(xvtfile) and os.path.exists(yvtfile)):
		xvt		=	np.load(xvtfile,mmap_mode='r')
		yvt		=	np.load(yvtfile,mmap_mode='r')
		#print(xvt.shape,yvt.shape)	
		xzoom	=	xvt[ptimenat-lenat:ptimenat+lenat]	
		yzoom	=	yvt[ptimenat-lenat:ptimenat+lenat]	
		#print(xzoom.shape,yzoom.shape)	
		np.save("{}{}_x_t_zoom_{}.npy".format(datadir,frbname,dm),xzoom)
		np.save("{}{}_y_t_zoom_{}.npy".format(datadir,frbname,dm),yzoom)	
		del(xvt)
		del(yvt)	
	else:
		print("extract_zoom ---- PANIC - File(s) not found!")
	
	return(2*lenat)

#	-------------------------------------------------------------------------------

def voltage_ds(frbname, dm, nchan):

	#	Construct voltage dynamic spectra for X and Y
	
	datadir		=	htrdir+frbname+'/data/'
	
	xvfile		=	"{}celebi/{}_x_t_zoom_{}.npy".format(datadir,frbname,dm)			
	yvfile		=	"{}celebi/{}_y_t_zoom_{}.npy".format(datadir,frbname,dm)			
	
	dslen		=	0
	
	if(os.path.exists(xvfile) and os.path.exists(yvfile)):
		xvdat		=	np.load(xvfile,mmap_mode='r')
		yvdat		=	np.load(yvfile,mmap_mode='r')
		xvdspec		=	generate_dynspec(xvdat, nchan)
		yvdspec		=	generate_dynspec(yvdat, nchan)
		np.save("{}reduced/{}_x_vds_{}_{}.npy".format(datadir,frbname,dm,nchan),xvdspec.T)
		np.save("{}reduced/{}_y_vds_{}_{}.npy".format(datadir,frbname,dm,nchan),yvdspec.T)		
		dslen		=	xvdspec.shape[0]
		del(xvdat)
		del(yvdat)
		del(xvdspec)
		del(yvdspec)
	else:
		print("voltage_ds ---- PANIC - File(s) not found!")

	return(dslen)

#	-------------------------------------------------------------------------------

def gen_stokes_ds(frbname, dm, nchan, fmhz0):

	#	Construct dynamic spectra for I, Q, U, V
	
	datadir		=	htrdir+frbname+'/data/reduced/'
	polcalfile	=	htrdir+frbname+'/data/celebi/'+frbname+'_polcal.txt'
		
	xvdsfile	=	"{}{}_x_vds_{}_{}.npy".format(datadir,frbname,dm,nchan)			
	yvdsfile	=	"{}{}_y_vds_{}_{}.npy".format(datadir,frbname,dm,nchan)			
	
	if(os.path.exists(xvdsfile) and os.path.exists(yvdsfile) and os.path.exists(polcalfile)):
		xvds	=	np.load(xvdsfile)
		yvds	=	np.load(yvdsfile)
		pcals	=	np.loadtxt(polcalfile)
		stksds	=	np.zeros((4,xvds.shape[0],xvds.shape[1]), dtype=float)
		
		for cc in range(0,nchan):
			stkst	=	calculate_stokes_unnormalised(xvds[cc], yvds[cc])
			for i in range(0,4):
				stksds[i,cc]	=	stkst[i]	
		
		#np.save("{}{}_stks_ds_{}_{}_unpolcal.npy".format(datadir,frbname,dm,nchan),stksds)
		
		#	Apply pol cal
		
		tus		=	pcals[1]
		offrad	=	pcals[2]
		chanwmhz=	NatBWmhz/nchan
		fmhzarr	=	np.arange(fmhz0-(NatBWmhz/2)+(chanwmhz/2), fmhz0+NatBWmhz/2, chanwmhz)
		mixrads	=	(2*np.pi*fmhzarr*tus) + offrad 	
		
		plt.plot(fmhzarr,np.rad2deg(mixrads))
		plt.show()
		
		print("\n\n ************************ \n")
		print("Polcaling...\n")
		print("Central frequency = %.2f MHz MA = %.2f deg"%(fmhz0,np.rad2deg(2*np.pi*fmhz0*tus + offrad)))
		print("Mixing angles <= %.1f deg"%np.nanmax(np.abs(np.rad2deg(mixrads))))
		print("\n ************************ \n\n")
		
		sinmix	=	np.sin(mixrads)
		cosmix	=	np.cos(mixrads)
		
		for cc in range(0,nchan):
			ucor	=	stksds[2,cc]*cosmix[cc] - stksds[3,cc]*sinmix[cc]
			vcor	=	stksds[2,cc]*sinmix[cc] + stksds[3,cc]*cosmix[cc]
			stksds[2,cc]	=	ucor
			stksds[3,cc]	=	vcor	
		
		np.save("{}{}_stks_ds_{}_{}_avg_1_1.npy".format(datadir,frbname,dm,nchan),stksds)	
		#print(stksds.shape)
		del(stksds)
		del(xvds)
		del(yvds)
	else:
		print("gen_stokes_ds ---- PANIC - File(s) not found!")

	return(nchan*Raw_time_res_ms)

#	-------------------------------------------------------------------------------

def favg_stokes(frbname, dm, nchan, ffac):

	#	Average dynamic spectra for I, Q, U, V along frequency axis
	
	datadir		=	htrdir+frbname+'/data/reduced/'
		
	stksdsfile	=	"{}{}_stks_ds_{}_{}_avg_1_1.npy".format(datadir,frbname,dm,nchan)			
	avgchans	=	int(round(nchan/ffac))
	
	if(os.path.exists(stksdsfile)):
		stksds	=	np.load(stksdsfile)
		favgds	=	np.nanmean(np.reshape(stksds,(4,avgchans,ffac,-1)),axis=2)
		np.save("{}{}_stks_ds_{}_{}_avg_{}_1.npy".format(datadir,frbname,dm,nchan,ffac),favgds)		
		#print(stksds.shape,favgds.shape)
		del(stksds)
		del(favgds)
	else:
		print("favg_stokes ---- PANIC - File not found!")

	return(avgchans)

#	-------------------------------------------------------------------------------

def avg_stokes_ds(frbname, dm, nchan, ffac, avgfac):

	#	Average dynamic spectra for I, Q, U, V along time axis
	
	datadir		=	htrdir+frbname+'/data/reduced/'
		
	stksdsfile	=	"{}{}_stks_ds_{}_{}_avg_{}_1.npy".format(datadir,frbname,dm,nchan,ffac)			
	avglen		=	0
	avgchan		=	int(round(nchan/ffac))
	
	if(os.path.exists(stksdsfile)):
		stksds	=	np.load(stksdsfile,mmap_mode='r')
		natlen	=	stksds.shape[2]
		avglen	=	int(float(natlen)/avgfac)
		avgds	=	np.nanmean(np.reshape(stksds[:,:,int((natlen-avglen*avgfac)/2):int((natlen+avglen*avgfac)/2)],(4,avgchan,-1,avgfac)),axis=3)
		np.save("{}{}_stks_ds_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac),avgds)		
		#print(stksds.shape,avgds.shape)
		del(stksds)
		del(avgds)
	else:
		print("avg_stokes_ds ---- PANIC - File not found!")

	return(avglen)

#	-------------------------------------------------------------------------------
		
def gen_noise_spec(frbname, dm, nchan, ffac, avgfac):

	#	Generate noise spectra
	
	datadir		=	htrdir+frbname+'/data/reduced/'
	plotdir		=	htrdir+frbname+'/plots/'
		
	stksdsfile	=	"{}{}_stks_ds_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)			
	noiselen	=	0
	
	if(os.path.exists(stksdsfile)):
		if(ZoomLenms > MaxWms):
			stksds	=	np.load(stksdsfile)
			natlen	=	stksds.shape[2]
			exlen	=	int(MaxWms/(Raw_time_res_ms*nchan*avgfac))
			stksdum	=	np.copy(stksds)
			stksdum[:,:,int((natlen-exlen)/2):int((natlen+exlen)/2)]	= np.nan
			meanspec=	np.nanmean(stksdum, axis=2)
			for ss in range(0,stksds.shape[0]):
				for ci in range(0,stksds.shape[1]):
					stksds[ss,ci]	=	stksds[ss,ci] - meanspec[ss,ci]
					stksdum[ss,ci]	=	stksdum[ss,ci] - meanspec[ss,ci]
			rmspec	=	np.nanstd(stksdum, axis=2)
			np.save("{}{}_stks_ds_zeromean_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac),stksds)
			np.save("{}{}_meanbkg_spec_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac),meanspec)	
			np.save("{}{}_noise_spec_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac),rmspec)					
			fig		=	plt.figure(figsize=(5,8))
			ax		=	fig.add_subplot(2,1,1)
			if(NatBWmhz < 0.0):
				plt.imshow(stksds[0],origin='upper',cmap='Blues',interpolation='none',aspect='auto',vmin=0.5*np.nanmin(rmspec[0]), vmax=4*np.nanmin(rmspec[0]))
			else:
				plt.imshow(stksds[0],origin='lower',cmap='Blues',interpolation='none',aspect='auto',vmin=0.5*np.nanmin(rmspec[0]), vmax=4*np.nanmin(rmspec[0]))
			plt.ylabel('Channel (Frequency $ \longrightarrow $)')
			plt.xlabel('Time')
			ax		=	fig.add_subplot(2,1,2)
			plt.plot(rmspec[0], 'b-', label='I')
			plt.plot(rmspec[1], 'r-', label='Q')
			plt.plot(rmspec[2], 'c-', label='U')
			plt.plot(rmspec[3], 'k-', label='V')
			plt.xlabel('Channel')
			plt.ylabel('RMS noise (arbitrary unit)')
			plt.legend(loc='upper right')
			plt.tight_layout()
			#plt.show()
			plt.savefig("{}{}_noise_spec_{}_{}_avg_{}_{}.png".format(plotdir,frbname,dm,nchan,ffac,avgfac))
			plt.close()
			del(stksds)	
			del(stksdum)	
		else:
			print("gen_noise_spec ---- PANIC - ZoomLenms is smaller than MaxWms!")
	else:
		print("gen_noise_spec ---- PANIC - File not found!")

	return(noiselen)

#	-------------------------------------------------------------------------------

def gen_bad_chan_list(frbname, dm, nchan, ffac, avgfac):

	#	Generate a list of bad/noisy channels
	
	datadir		=	htrdir+frbname+'/data/reduced/'
		
	noisefile	=	"{}{}_noise_spec_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)			
	
	avgchan		=	int(round(nchan/ffac))	
	chanq		=	np.ones(avgchan, dtype=int)
	
	if(os.path.exists(noisefile)):
		
		#	Flag channels based on noise
		
		rmspec	=	np.load(noisefile)		
		lims	=	st.norm.interval(1.0-(1.0/float(avgchan)), loc=np.nanmedian(rmspec[0]), scale=1.48*np.nanmedian(np.abs(rmspec[0] - np.nanmedian(rmspec[0]))))
		
		for ci in range(0, avgchan):
			if(rmspec[0,ci]>lims[1]):
				chanq[ci]	=	0
		
		#	Flag edge channels
		
		edgech	=	int(np.fix(avgchan*EdgeFrac)+1)
		chanq[:edgech]				=	0
		chanq[(avgchan-edgech):]	=	0
		print("Flagged channels %d - %d and %d - %d\n"%(0,edgech-1, avgchan-edgech,avgchan-1))	
		
		np.save("{}{}_chanq_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac),chanq)
	
		del(rmspec)
	else:
		print("gen_bad_chan_list ---- PANIC - File not found!")

	return(0)

#	-------------------------------------------------------------------------------

def gen_I_ts(frbname, dm, nchan, ffac, avgfac):

	#	Generate I time series
	
	datadir		=	htrdir+frbname+'/data/reduced/'
	plotdir		=	htrdir+frbname+'/plots/'
		
	stksdsfile	=	"{}{}_stks_ds_zeromean_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)
	noisefile	=	"{}{}_noise_spec_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)	
	chanqfile	=	"{}{}_chanq_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)			
	
	avgchan		=	int(round(nchan/ffac))	
	
	if(os.path.exists(chanqfile)):
		chanq	=	np.load(chanqfile)
	else:
		chanq		=	np.ones(avgchan, dtype=int)
		print("gen_I_ts ---- No bad channel list found...")
	
	if(os.path.exists(stksdsfile) and os.path.exists(noisefile)):
		stksds	=	np.load(stksdsfile)
		rmspec	=	np.load(noisefile)
		badcs	=	np.where(chanq==0)
		stksds[:,badcs]	=	np.nan
		rmspec[:,badcs]	=	np.nan
		inoise	=	np.sqrt(np.nansum(rmspec[0]**2))/nchan
		#print(inoise)
		its		=	np.nanmean(stksds[0],axis=0)
		np.save("{}{}_I_ts_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac),its)
		fig		=	plt.figure(figsize=(10,10))
		ax0		=	fig.add_subplot(2,1,1)
		ax0.set_title("Time resolution = %.6f ms"%(avgfac*nchan*Raw_time_res_ms))
		plt.plot(its,'b-')
		plt.axhline(c='k',ls='-')
		plt.axhline(y=inoise, c='k',ls='--')
		plt.axhline(y=-inoise, c='k',ls='--')
		plt.axhline(y=5*inoise, c='k',ls=':')
		plt.axhline(y=-5*inoise, c='k',ls=':')		
		plt.xlim(0,len(its)-1)
		ax0		=	fig.add_subplot(2,1,2)
		if(NatBWmhz < 0.0):
			plt.imshow(stksds[0],origin='upper',interpolation='none',cmap='Blues',aspect='auto',vmin=0.5*np.nanmin(rmspec[0]), vmax=4*np.nanmin(rmspec[0]))
		else:
			plt.imshow(stksds[0],origin='lower',interpolation='none',cmap='Blues',aspect='auto',vmin=0.5*np.nanmin(rmspec[0]), vmax=4*np.nanmin(rmspec[0]))
		plt.ylabel('Channel (Frequency $ \longrightarrow $)')
		plt.show()	
		#plt.savefig("{}{}_stks_I_full_{}_{}_avg_{}_{}.png".format(plotdir,frbname,dm,nchan,ffac,avgfac))
		plt.close()
		del(stksds)
		del(its)
	else:
		print("gen_I_ts ---- PANIC - File(s) not found!")

	return(0)

#	-------------------------------------------------------------------------------











































