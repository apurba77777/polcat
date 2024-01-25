#
#	Functions for FRB polarization analysis
#
#								AB, April 2023

#	Function list
#
#			cal_width(frbname, dm, nchan, ffac, avgfac, tbasems):
#						Estimate FRB width
#
#			estimate_rm(frbname, dm, nchan, ffac, avgfac, fmhz0, lwms, rwms, phirange, dphi, startchan, endchan, tpeakms):
#						Estimate rotation measure
#
#			plot_polang(frbname, dm, nchan, ffac, avgfac, fmhz0, tpeakms, lwms, rwms, phi0, rm0):
#						Plot polarization angle vs frequency
#
#			unfarot(frbname, dm, nchan, ffac, avgfac, fmhz0, rm0):
#						Generate RM corrected dynamic spectrum 
#
#			sub_ts(frbname, dm, nchan, ffac, avgfac, tbasems, fmhz0, tpeakms):
#						Plot intensity profile for different subbands
#
#			sub_spec(frbname, dm, nchan, ffac, avgfac, fmhz0,tbasems, tpeakms):
#						Plot spectra for different components
#
#
#	--------------------------	Import modules	---------------------------

import os, sys
import matplotlib as mpl
import numpy as np
import matplotlib.ticker as ticker
from globalpars import *
from htrfns.dynamic_fns import *
from htrfns.scintillation_fns import *
from htrfns.polfns import *

mpl.rcParams['pdf.fonttype']	= 42
mpl.rcParams['ps.fonttype'] 	= 42
mpl.rcParams['savefig.dpi'] 	= 600
mpl.rcParams['font.family'] 	= 'sans-serif'
mpl.rcParams['font.size']		= 7

#	--------------------------	Analysis functions	-------------------------------

def cal_width(frbname, dm, nchan, ffac, avgfac, tbasems):

	#	Calculate FRB width
	
	datadir		=	htrdir+frbname+'/data/reduced/'
	plotdir		=	htrdir+frbname+'/plots/'
	
	chanqfile	=	"{}{}_chanq_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)			
	avgchan		=	int(round(nchan/ffac))	
	
	if(os.path.exists(chanqfile)):
		chanq	=	np.load(chanqfile)
	else:
		chanq		=	np.ones(avgchan, dtype=int)
		print("cal_width ---- No bad channel list found...")
		
	itsfile		=	"{}{}_I_ts_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)	
	stksdsfile	=	"{}{}_stks_ds_zeromean_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)	
	noisefile	=	"{}{}_noise_spec_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)
	tresms		=	avgfac*nchan*Raw_time_res_ms
	eqwms		=	0
	w95ms		=	0
	lw95ms		=	0
	rw95ms		=	0
	wms			=	0
	lwms		=	0
	rwms		=	0
	tpeak		=	0
	
	if((os.path.exists(itsfile)) and (os.path.exists(stksdsfile))):
		its		=	np.load(itsfile)
		stksds	=	np.load(stksdsfile)
		badcs	=	np.where(chanq==0)
		stksds[:,badcs]	=	np.nan
		rmspec	=	np.load(noisefile)
		itsdum	=	np.copy(its)
		natlen	=	len(its)
		exlen	=	int(MaxWms/(Raw_time_res_ms*nchan*avgfac))
		itsdum[int((natlen-exlen)/2):int((natlen+exlen)/2)]	= np.nan
		goodinds=	np.arange(natlen)[np.isfinite(itsdum)]
		dumean	=	np.nanmean(itsdum)
		its		=	its - dumean 
		itsdum	=	itsdum - dumean 
		irms	=	np.nanstd(itsdum)
		peak0	=	np.nanmax(its)
		tpeak	=	np.argmax(its)
		fl0		=	peak0
		dlfl	=	0.0
		drfl	=	0.0
		istart	=	tpeak
		iend	=	tpeak
		
		print("Peak found at %.3f ms"%(tpeak*tresms))
		print("Peak S/N = %.2f"%(its[tpeak]/irms))
		
		wleft	=	int(round(tpeak-MaxWms/tresms))
		wright	=	int(round(tpeak+MaxWms/tresms))
		
		while((np.nansum(its[wleft:istart])>5*dlfl) or (its[istart-1] > 2*irms)):
			istart 	=	istart - 1
			dlfl	=	irms*np.sqrt(float(tpeak-istart))
		while((np.nansum(its[iend+1:wright])>5*drfl) or (its[iend+1] > 2*irms)):
			iend 	=	iend + 1
			drfl	=	irms*np.sqrt(float(iend+1-tpeak))
		
		dfl0	=	irms*np.sqrt(float(iend+1-istart))
		while(np.nansum(its[wleft:istart])>2*dfl0):
			istart 	=	istart - 1
			dfl0	=	irms*np.sqrt(float(iend+1-istart))
		while(np.nansum(its[iend+1:wright])>2*dfl0):
			iend 	=	iend + 1
			dfl0	=	irms*np.sqrt(float(iend+1-istart))
		
		fl0		=	np.nansum(its[istart:iend+1])
		snrfl	=	fl0/(irms*np.sqrt(iend+1-istart))
		eqwms	=	tresms*fl0/peak0
		lwms	=	-tresms*(tpeak-istart)
		rwms	=	tresms*(iend+1-tpeak)
		wms		=	-lwms+rwms
		
		lfl		=	np.nansum(its[istart:tpeak]) + its[tpeak]/2
		rfl		=	np.nansum(its[tpeak+1:iend+1]) + its[tpeak]/2
		l95		=	tpeak
		r95		=	tpeak
		
		while((np.nansum(its[l95:tpeak])+its[tpeak]/2.0)<0.95*lfl):
			l95	=	l95 - 1
		while((np.nansum(its[tpeak+1:r95+1])+its[tpeak]/2.0)<0.95*rfl):
			r95	=	r95 + 1
		
		lw95ms	=	-tresms*(tpeak-l95)
		rw95ms	=	tresms*(r95+1-tpeak)
		w95ms	=	-lw95ms+rw95ms
		
		reltms	=	np.arange(-tpeak, len(its) -tpeak - 0.5, 1.0, dtype=float)*tresms		
		fig		=	plt.figure(figsize=(15,10))
		ax0		=	fig.add_subplot(2,2,1)
		ax0.set_title(r"Time resolution = %.6f ms, W$_{e}$ = %.2f ms, W$_{95}$ = %.2f ms, W = %.2f ms, S/N = %.1f"%(tresms, eqwms, w95ms, wms, snrfl))
		plt.plot(reltms, its,'b-')
		#plt.axvline(lw95ms, c='r', ls='--')
		#plt.axvline(rw95ms, c='r', ls='--')
		#plt.axvline(lwms, c='k', ls='--')
		#plt.axvline(rwms, c='k', ls='--')
		plt.xlim((-30*eqwms, 30*eqwms))
		plt.xlabel('Time - Peak (ms)')
		plt.ylabel('Flux density (arbitrary unit)')
		
		ax1		=	fig.add_subplot(2,2,2)
		ax1.set_title("Zoomed")
		plt.plot(reltms, its,'b-')
		#plt.axvline(lw95ms, c='r', ls='--')
		#plt.axvline(rw95ms, c='r', ls='--')
		#plt.axvline(lwms, c='k', ls='--')
		#plt.axvline(rwms, c='k', ls='--')
		plt.xlim((-tbasems, tbasems))
		#plt.xlim((2*lwms, 2*rwms))
		plt.xlabel('Time - Peak (ms)')
		#plt.ylabel('Flux density (arbitrary unit)')
		
		ax2		=	fig.add_subplot(2,2,3)
		if(NatBWmhz < 0.0):
			plt.imshow(stksds[0],origin='upper',aspect='auto',cmap='plasma',interpolation='none',vmin=0, vmax=5*np.nanmin(rmspec[0]))
		else:
			plt.imshow(stksds[0],origin='lower',aspect='auto',cmap='plasma',interpolation='none',vmin=0, vmax=5*np.nanmin(rmspec[0]))
		#plt.axvline(tpeak+lw95ms/tresms, c='r', ls='--')
		#plt.axvline(tpeak+rw95ms/tresms, c='r', ls='--')
		#plt.axvline(tpeak+lwms/tresms, c='k', ls='--')
		#plt.axvline(tpeak+rwms/tresms, c='k', ls='--')
		plt.xlim((tpeak-30*eqwms/tresms, tpeak+30*eqwms/tresms))
		plt.xlabel('Time (ms)')
		plt.ylabel('Channel (Frequency $ \longrightarrow $)')
		
		ax2		=	fig.add_subplot(2,2,4)
		if(NatBWmhz < 0.0):
			plt.imshow(stksds[0],origin='upper',aspect='auto',cmap='plasma',interpolation='none',vmin=0, vmax=5*np.nanmin(rmspec[0]))
		else:
			plt.imshow(stksds[0],origin='lower',aspect='auto',cmap='plasma',interpolation='none',vmin=0, vmax=5*np.nanmin(rmspec[0]))
		#plt.axvline(tpeak+lw95ms/tresms, c='r', ls='--')
		#plt.axvline(tpeak+rw95ms/tresms, c='r', ls='--')
		#plt.axvline(tpeak+lwms/tresms, c='k', ls='--')
		#plt.axvline(tpeak+rwms/tresms, c='k', ls='--')
		plt.xlim((tpeak-tbasems/tresms, tpeak+tbasems/tresms))
		plt.xlabel('Time (ms)')
		#plt.ylabel('Channel')
		#plt.savefig("{}{}_stks_I_w_{}_{}_avg_{}_{}.png".format(plotdir,frbname,dm,nchan,ffac,avgfac))
		plt.show()	
		del(its)
	else:
		print("cal_width ---- PANIC - File not found!")

	return(tpeak,eqwms,w95ms,lw95ms,rw95ms,wms,lwms,rwms)

#	-------------------------------------------------------------------------------

def sub_ts(frbname, dm, nchan, ffac, avgfac, tbasems, fmhz0, tpeakms):

	#	Plot intensity profile for different subbands
	
	datadir		=	htrdir+frbname+'/data/reduced/'
	plotdir		=	htrdir+frbname+'/plots/'
	
	chanqfile	=	"{}{}_chanq_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)			
	avgchan		=	int(round(nchan/ffac))	
	
	compfile	=	htrdir+frbname+'/subands.txt'
	colarr		=	['r', 'c', 'm', 'b']
	lsarr		=	['-', '-.', ':', '--']
	
	edgefile	=	htrdir+frbname+'/compedge.txt'
	
	if(os.path.exists(chanqfile)):
		chanq	=	np.load(chanqfile)
	else:
		chanq		=	np.ones(avgchan, dtype=int)
		print("sub_ts ---- No bad channel list found...")
		
	itsfile		=	"{}{}_I_ts_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)	
	stksdsfile	=	"{}{}_stks_ds_zeromean_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)	
	noisefile	=	"{}{}_noise_spec_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)
	tresms		=	avgfac*nchan*Raw_time_res_ms
	tpeak		=	tpeakms/tresms
	
	if((os.path.exists(itsfile)) and (os.path.exists(stksdsfile)) and (os.path.exists(compfile))):
		its		=	np.load(itsfile)
		stksds	=	np.load(stksdsfile)
		rmspec	=	np.load(noisefile)
		badcs	=	np.where(chanq==0)
		stksds[:,badcs]	=	np.nan
		
		subfracs=	np.loadtxt(compfile)
		#comedges=	np.loadtxt(edgefile)
		nsub	=	subfracs.shape[0]-1
		print("Found %d sub-bands...\n"%(nsub))
		
		isubts	=	[]
		subchs	=	[]
		startch	=	0
		for si in range(0,nsub):
			endch	=	min(avgchan, startch+int(subfracs[si+1]*avgchan))
			if(si == (nsub-1)):
				endch	=	avgchan-1
			#print(startch, endch)
			subchs.append([startch,endch])
			isubts.append(np.nanmean(stksds[0,startch:endch],axis=0))			
			startch	=	endch		
		
		chanwmhz=	NatBWmhz/avgchan
		fmhzarr	=	np.arange(fmhz0-(NatBWmhz/2)+(chanwmhz/2), fmhz0+NatBWmhz/2, chanwmhz)
		
		reltms	=	np.arange(-tpeak, len(its) -tpeak - 0.5, 1.0, dtype=float)*tresms		
		fig		=	plt.figure(figsize=(18.3/intocm,17.0/intocm))			
		ax2		=	fig.add_subplot(2,1,2)
		ax2.text(0.04,0.9,'b',fontsize=8,fontweight='bold',transform=ax2.transAxes)
		#for cs in comedges[2:-1]:
			#plt.axvline(cs,c='grey',ls='--', lw=0.8)
		plt.axhline(c='k', ls='--',lw=0.25)
		plt.plot(reltms, its / np.nanmax(its),'k-',lw=0.5,label='Full band')
		plt.xlim((tbasems[0], tbasems[1]))
		ax2.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
		plt.ylim([-0.1,1.1])
		plt.xlabel('Time (ms)')
		plt.ylabel('Normalized flux density')
		#plt.yticks([])
		#plt.ylabel('Flux density')
		ax2.legend(loc='upper right')
		
		ax0		=	fig.add_subplot(2,1,1)
		ax0.text(0.04,0.9,'a',fontsize=8,fontweight='bold',transform=ax0.transAxes)
		plt.axhline(c='k', ls='--', lw=0.25)
		for si in range(0,nsub):			
			plt.plot(reltms, isubts[si] / np.nanmax(its),c=colarr[si%nsub], ls='-',lw=0.5,label="%d MHz"%(np.nanmedian(fmhzarr[subchs[si][0]:subchs[si][1]])))
		
		plt.xlim((tbasems[0], tbasems[1]))	
		plt.ylabel('Normalized flux density')	
		ax0.legend(loc='upper right')
		#plt.xlabel('Time - Peak (ms)')	
		#plt.ylabel('Flux density')
		plt.xticks([])
		#plt.yticks([])
		
		plt.tight_layout(h_pad=0, w_pad=5)			
		plt.savefig("{}{}_sub_{}_bands_{}_{}_avg_{}_{}.eps".format(plotdir,frbname,nsub,dm,nchan,ffac,avgfac), transparent=True, format='eps')
		plt.savefig("{}{}_sub_{}_bands_{}_{}_avg_{}_{}.pdf".format(plotdir,frbname,nsub,dm,nchan,ffac,avgfac), transparent=True, format='pdf')
		plt.show()
		del(stksds)
		del(its)
	else:
		print("sub_ts ---- PANIC - File not found!")

	return(tpeak)

#	-------------------------------------------------------------------------------

def sub_spec(frbname, dm, nchan, ffac, avgfac, fmhz0,tbasems, tpeakms):

	#	Plot spectra for different components
	
	datadir		=	htrdir+frbname+'/data/reduced/'
	plotdir		=	htrdir+frbname+'/plots/'
	
	chanqfile	=	"{}{}_chanq_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)			
	avgchan		=	int(round(nchan/ffac))	
	
	compfile	=	htrdir+frbname+'/subands.txt'	
	edgefile	=	htrdir+frbname+'/complist.txt'
	edge2file	=	htrdir+frbname+'/compedge.txt'
	colarr		=	['k','c', 'r', 'm', 'g']
	markarr		=	['*','x','o','d']
	
	if(os.path.exists(chanqfile)):
		chanq	=	np.load(chanqfile)
	else:
		chanq		=	np.ones(avgchan, dtype=int)
		print("sub_spec ---- No bad channel list found...")
		
	itsfile		=	"{}{}_I_ts_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)	
	stksdsfile	=	"{}{}_stks_ds_zeromean_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)	
	noisefile	=	"{}{}_noise_spec_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)
	tresms		=	avgfac*nchan*Raw_time_res_ms
	tpeak		=	tpeakms/tresms
	
	if((os.path.exists(itsfile)) and (os.path.exists(stksdsfile)) and (os.path.exists(edgefile))):
		its		=	np.load(itsfile)
		stksds	=	np.load(stksdsfile)
		rmspec	=	np.load(noisefile)
		badcs	=	np.where(chanq==0)
		goodcs	=	np.where(chanq>0)
		stksds[:,badcs]	=	np.nan
		chanwmhz=	NatBWmhz/avgchan
		fmhzarr	=	np.arange(fmhz0-(NatBWmhz/2)+(chanwmhz/2), fmhz0+NatBWmhz/2, chanwmhz)
		
		subfracs=	np.loadtxt(compfile)
		comedges=	np.loadtxt(edgefile)
		coedges	=	np.loadtxt(edge2file)
		nband	=	subfracs.shape[0]-1
		nsub	=	comedges.shape[0]-2
		print("Found %d sub-bands and %d components...\n"%(nband,nsub))
		
		isubspec=	np.zeros((nsub,nband),dtype=float)
		isubnois=	np.zeros((nsub,nband),dtype=float)
		centfmhz=	np.zeros(nband,dtype=float)
		centms	=	np.zeros(nsub,dtype=float)
		for si in range(0,nsub):
			istart		=	int(tpeak + (comedges[si+1][0]/tresms))
			iend		=	int(tpeak + (comedges[si+1][1]/tresms) - 1)
			centms[si]	=	((comedges[si+1][0]+comedges[si+1][1])/2.0)
			startch	=	0
			for ci in range(0,nband):
				endch	=	min(avgchan, startch+int(subfracs[ci+1]*avgchan))
				if(ci == (nsub-1)):
					endch	=	avgchan
				isubspec[si,ci]	=	np.nanmean(stksds[0,startch:endch,istart:iend+1])
				isubnois[si,ci]	=	np.sqrt(np.nanmean(rmspec[0,startch:endch]**2))/np.sqrt(float(iend+1-istart))	
				centfmhz[ci]	=	np.nanmean(fmhzarr[startch:endch])					
				startch			=	endch	
								
		reltms	=	np.arange(-tpeak, len(its) -tpeak - 0.5, 1.0, dtype=float)*tresms		
		fig		=	plt.figure(figsize=(15,6))			
		ax2		=	fig.add_subplot(1,2,1)
		for cs in coedges[1:]:
			plt.axvline(cs,c='c',ls='--', lw=0.8)
		for si in range(0,nsub):
			plt.text(centms[si],1.01*np.nanmax(its),str(si+1),fontsize=10,color='c')
		plt.axhline(c='c', ls='--',lw=0.5)
		plt.plot(reltms, its,'b-',label='Full band')
		plt.xlim((tbasems[0], tbasems[1]))
		plt.xlabel('Time - Peak (ms)')
		plt.ylabel('Flux density (arbitrary unit)')
		#ax2.legend(loc='upper right')
		
		ax2		=	fig.add_subplot(1,2,2)
		for si in range(0,nsub):
			plt.errorbar(centfmhz, isubspec[si], isubnois[si], fmt=colarr[si%5]+markarr[si%4]+'-',markersize=10,capsize=2,label='C'+str(si+1))
		plt.yscale('log')
		plt.xlabel("Frequency (MHz)")
		ax2.legend(loc='upper left')
		
		plt.tight_layout(h_pad=0, w_pad=2)			
		plt.savefig("{}{}_sub_{}_spec_{}_{}_avg_{}_{}.png".format(plotdir,frbname,nsub,dm,nchan,ffac,avgfac))
		plt.show()
		del(stksds)
		del(its)
	else:
		print("sub_spec ---- PANIC - File not found!")

	return(tpeak)

#	-------------------------------------------------------------------------------

def estimate_rm(frbname, dm, nchan, ffac, avgfac, fmhz0, lwms, rwms, phirange, dphi, startchan, endchan, tpeakms):

	#	Estimate rotation measure
	
	datadir		=	htrdir+frbname+'/data/reduced/'
	plotdir		=	htrdir+frbname+'/plots/'
	
	chanqfile	=	"{}{}_chanq_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)			
	avgchan		=	int(round(nchan/ffac))
	
	if(os.path.exists(chanqfile)):
		chanq	=	np.load(chanqfile)
	else:
		chanq		=	np.ones(avgchan, dtype=int)
		print("cal_width ---- No bad channel list found...")
		
	stksdsfile	=	"{}{}_stks_ds_zeromean_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)
	noisefile	=	"{}{}_noise_spec_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)			
		
	tresms		=	avgfac*nchan*Raw_time_res_ms
	res_rmsynth	=	[0.0,0.0,0.0,0.0]
	res_rmnest	=	[0.0,0.0,0.0,0.0]
	res_rmtool	=	[0.0,0.0,0.0,0.0]
	
	if(os.path.exists(stksdsfile) and os.path.exists(noisefile)):
		stksds		=	np.load(stksdsfile)
		noispec		=	np.load(noisefile)
		goodc		=	np.where(chanq>0)[0]
		badcs		=	np.where(chanq==0)
		stksds[:,badcs]	=	np.nan
		its			=	np.nanmean(stksds[0],axis=0)
		tpeak		=	tpeakms/tresms
		istart		=	int(tpeak + (lwms/tresms))
		iend		=	int(tpeak + (rwms/tresms) - 1)
		avgchans	=	stksds.shape[1]
		chanwmhz	=	NatBWmhz/avgchans
		fmhzarr		=	np.arange(fmhz0-(NatBWmhz/2)+(chanwmhz/2), fmhz0+NatBWmhz/2, chanwmhz)
		lm2arr		=	(ccC*1.0e-8 / fmhzarr)**2
		lm20		=	(ccC*1.0e-8 / fmhz0)**2
		
		goodc0		=	goodc[startchan:endchan+1]
		
		ispec		=	np.nanmean(stksds[0,:,istart:iend+1], axis=1)
		vspec		=	np.nanmean(stksds[3,:,istart:iend+1], axis=1)
		qspec0		=	np.nanmean(stksds[1,:,istart:iend+1], axis=1)
		uspec0		=	np.nanmean(stksds[2,:,istart:iend+1], axis=1)
		noispec		=	noispec/np.sqrt(float(iend+1-istart))
								
		iqu			=	(ispec[goodc0],qspec0[goodc0],uspec0[goodc0])
		eiqu		=	(noispec[0][goodc0],noispec[1][goodc0],noispec[2][goodc0])
		
		iquv		=	(ispec[goodc0],qspec0[goodc0],uspec0[goodc0], vspec[goodc0])
		eiquv		=	(noispec[0][goodc0],noispec[1][goodc0],noispec[2][goodc0],noispec[3][goodc0])
		
		res_rmtool	=	fr_rmtool(fmhzarr[goodc0]/1.0e3, iquv, eiquv)
		
		del(stksds)
		del(noispec)
	else:
		print("estimate_rm ---- PANIC - File not found!")
	
	print("\nResults from RMtool (RM synthesis) \n")
	print("RM = %.2f +/- %.2f rad/m2   PolAng0 = %.2f +/- %.2f deg\n"%(res_rmtool[0], res_rmtool[1], res_rmtool[2], res_rmtool[3]))
	
	return(res_rmsynth, res_rmnest, res_rmtool, tpeak)

#	-------------------------------------------------------------------------------

def plot_polang(frbname, dm, nchan, ffac, avgfac, fmhz0, tpeakms, lwms, rwms, phi0, rm0):

	#	Plot polarization angle vs frequency
	
	datadir		=	htrdir+frbname+'/data/reduced/'
	plotdir		=	htrdir+frbname+'/plots/'
	
	chanqfile	=	"{}{}_chanq_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)			
	avgchan		=	int(round(nchan/ffac))
	
	if(os.path.exists(chanqfile)):
		chanq	=	np.load(chanqfile)
	else:
		chanq	=	np.ones(avgchan, dtype=int)
		print("No bad channel list found...")
		
	stksdsfile	=	"{}{}_stks_ds_zeromean_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)
	noisefile	=	"{}{}_noise_spec_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)			
		
	tresms		=	avgfac*nchan*Raw_time_res_ms
	tpeak		=	tpeakms/tresms
	istart		=	int(tpeak + (lwms/tresms))
	iend		=	int(tpeak + (rwms/tresms) - 1)
	
	if(os.path.exists(stksdsfile) and os.path.exists(noisefile)):
		stksds		=	np.load(stksdsfile)
		noispec		=	np.load(noisefile)
		badcs		=	np.where(chanq==0)
		stksds[:,badcs]	=	np.nan
		avgchans	=	stksds.shape[1]
		chanwmhz	=	NatBWmhz/avgchans
		fmhzarr		=	np.arange(fmhz0-(NatBWmhz/2)+(chanwmhz/2), fmhz0+NatBWmhz/2, chanwmhz)
		lm2arr		=	1.0e-4*((ccC*1.0e-6 / fmhzarr)**2)
		lm20		=	1.0e-4*((ccC*1.0e-6 / fmhz0)**2)
				
		ispec		=	np.nanmean(stksds[0,:,istart:iend+1], axis=1)
		vspec		=	np.nanmean(stksds[3,:,istart:iend+1], axis=1)
		qspec0		=	np.nanmean(stksds[1,:,istart:iend+1], axis=1)
		uspec0		=	np.nanmean(stksds[2,:,istart:iend+1], axis=1)
		noispec		=	noispec/np.sqrt(float(iend+1-istart))
		
		goodc			=	np.where(chanq>0)
		
		lspec		=	np.sqrt(qspec0**2 + uspec0**2)
		phispec		=	np.rad2deg(0.5*np.unwrap(np.arctan2(uspec0[goodc],qspec0[goodc]), period=2*np.pi))
		phimodel	=	np.rad2deg(np.deg2rad(phi0) + rm0*lm2arr) 
		dphispec	=	np.rad2deg(0.5*np.sqrt((uspec0*noispec[1])**2 + (qspec0*noispec[2])**2) / (uspec0**2 + qspec0**2))
		
		fig		=	plt.figure(figsize=(15,5))
		ax0		=	fig.add_subplot(1,2,1)
		plt.plot(fmhzarr[goodc], uspec0[goodc], 'g+', label='U')
		plt.plot(fmhzarr[goodc], qspec0[goodc], 'rs', label='Q')
		plt.plot(fmhzarr[goodc], vspec[goodc], 'co', label='V')
		plt.plot(fmhzarr[goodc], ispec[goodc], 'b*', label='I')
		plt.plot(fmhzarr[goodc], lspec[goodc], 'kx', label='L')
		plt.xlabel('Frequency (MHz)')
		plt.legend(loc='upper right')
					
		ax3		=	fig.add_subplot(1,2,2)
		plt.errorbar(fmhzarr[goodc], pangdiffdeg(phispec,phimodel[goodc]),dphispec[goodc],fmt='*', c='b', capsize=4, label='Observed - Model')
		#plt.plot(lm2arr[goodc], pangdiffdeg(phispec,0.0)[goodc],'k+', label='Observed')
		#plt.plot(lm2arr, phimodel,'cx', label='Model')
		plt.axhline(y=0, c='k', ls='--', lw=0.5)
		plt.ylim([-45.0,45.0])
		plt.legend(loc='upper right')
		plt.xlabel('Frequency (MHz)')
		plt.ylabel('$\Delta \phi$ (deg)')
				
		#plt.savefig("{}{}_polang_{}_{}_avg_{}_{}_rm_{}.png".format(plotdir,frbname,dm,nchan,ffac,avgfac,rm0))
		plt.show()
		
		del(stksds)
		del(noispec)
	else:
		print("plot_polang ---- PANIC - File not found!")

	return(0)

#	-------------------------------------------------------------------------------

def unfarot(frbname, dm, nchan, ffac, avgfac, fmhz0, rm0):

	#	Generate RM corrected dynamic spectrum 
	datadir		=	htrdir+frbname+'/data/reduced/'	
	stksdsfile	=	"{}{}_stks_ds_zeromean_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)			
		
	if(os.path.exists(stksdsfile)):
		stksds		=	np.load(stksdsfile)
		stks0		=	np.copy(stksds)
		avgchans	=	stksds.shape[1]
		chanwmhz	=	NatBWmhz/avgchans
		fmhzarr		=	np.arange(fmhz0-(NatBWmhz/2)+(chanwmhz/2), fmhz0+NatBWmhz/2, chanwmhz)
		lm2arr		=	(ccC*1.0e-8 / fmhzarr)**2
		lm20		=	(ccC*1.0e-8 / fmhz0)**2
		
		for ci in range (0,avgchans):
			rotang		=	-2*rm0*(lm2arr[ci]-lm20)
			stks0[1,ci]	=	stksds[1,ci]*np.cos(rotang) - stksds[2,ci]*np.sin(rotang) 
			stks0[2,ci]	=	stksds[2,ci]*np.cos(rotang) + stksds[1,ci]*np.sin(rotang) 
				
		np.save("{}{}_stks_ds_unfarot_{}_{}_avg_{}_{}_rm_{:.2f}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac,rm0), stks0)	
				
		del(stksds)
		del(stks0)		
	else:
		print("unfarot ---- PANIC - File not found!")

	return(0)

#	-------------------------------------------------------------------------------































