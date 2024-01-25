#
#	Functions for FRB polarization
#
#								AB, October 2023
#
#	Function list
#
#			picklit(frbname, frbdm, nchan, ffac, tavgfac, fmhz0, frbrm0, twinms, tpeakms, subid)
#					Process and pickle the outputs
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
from scipy.optimize import curve_fit

mpl.rcParams['pdf.fonttype']	= 42
mpl.rcParams['ps.fonttype'] 	= 42
mpl.rcParams['savefig.dpi'] 	= 600
mpl.rcParams['font.family'] 	= 'sans-serif'
mpl.rcParams['font.size']		= 7

#	------------------------------------------------------------------------------------------------------------

def picklit(frbname, dm, nchan, ffac, avgfac, fmhz0, rm0, tbasems, tpeakms, subid):

#	Process and pickle the outputs
		
	datadir		=	htrdir+frbname+'/data/reduced/'
	odatadir	=	htrdir+frbname+'/data/outputs/'
	compfile	=	htrdir+frbname+'/subands.txt'
	
	itsfile		=	"{}{}_I_ts_{}_{}_avg_{}_{}.npy".format(datadir,frbname,dm,nchan,ffac,avgfac)		
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
	tpeak		=	tpeakms/tresms
	
	if(os.path.exists(stksdsfile) and os.path.exists(noisefile) and os.path.exists(compfile)):
		
		its				=	np.load(itsfile)						
		stksds			=	np.load(stksdsfile)
		noispec			=	np.load(noisefile)
		badcs			=	np.where(chanq==0)
		stksds[:,badcs]	=	np.nan
		noispec[:,badcs]=	np.nan
		
		chanwmhz		=	NatBWmhz/avgchan
		fmhzarr			=	np.arange(fmhz0-(NatBWmhz/2)+(chanwmhz/2), fmhz0+NatBWmhz/2, chanwmhz)
		
		subfracs		=	np.loadtxt(compfile)
		nsub			=	subfracs.shape[0]-1
				
		startch			=	max(0, int(np.nansum(subfracs[:subid])*avgchan))
		endch			=	min(avgchan, startch+int(subfracs[subid]*avgchan))
		
		startmhz		=	np.nanmin(fmhzarr[startch:endch])
		endmhz			=	np.nanmax(fmhzarr[startch:endch])
		centmhz			=	np.nanmedian(fmhzarr[startch:endch])
		
		goodc			=	np.where(chanq>0)			
		ngoodc			=	np.count_nonzero(chanq[startch:endch])
			
		print("Sub-band %d (of %d)...\n"%(subid,nsub))
		print("Channels    %d -- %d (of %d) ----- %d usable -------- %.1f MHz (%.1f -- %.1f)\n"%(startch, endch, nchan, ngoodc, centmhz, startmhz, endmhz))
		
		reltms			=	np.arange(-tpeak, stksds.shape[2] -tpeak - 0.5, 1.0, dtype=float)*tresms
		
		itsub			=	np.nanmean(stksds[0,startch:endch],axis=0)					
		noistks			=	np.sqrt(np.nansum(noispec[:,startch:endch]**2,axis=1))/ngoodc
		#print(noistks)
		qtsub			=	np.nanmean(stksds[1,startch:endch],axis=0)	
		utsub			=	np.nanmean(stksds[2,startch:endch],axis=0)	
		vtsub			=	np.nanmean(stksds[3,startch:endch],axis=0)			
							
		lts				=	np.sqrt(utsub**2 + qtsub**2)
		lts				=	noistks[0]*np.sqrt((lts/noistks[0])**2 - 1.0)						
		elts			=	np.sqrt((qtsub*noistks[1])**2 + (utsub*noistks[2])**2)/lts
		pts				=	np.sqrt(lts**2 + vtsub**2)
		epts			=	np.sqrt((qtsub*noistks[1])**2 + (utsub*noistks[2])**2 + (vtsub*noistks[3])**2)/pts
				
		phits			=	np.rad2deg(0.5*np.arctan2(utsub,qtsub))		
		dphits			=	np.rad2deg(0.5*np.sqrt((utsub*noistks[1])**2 + (qtsub*noistks[2])**2) / (utsub**2 + qtsub**2))						
		psits			=	np.rad2deg(0.5*np.arctan2(vtsub,lts))		
		dpsits			=	np.rad2deg(0.5*np.sqrt((vtsub*elts)**2 + (lts*noistks[3])**2) / (vtsub**2 + lts**2))	
				
		vfrac			=	vtsub/itsub
		lfrac			=	lts/itsub
		pfrac			=	pts/itsub		
		qfrac			=	qtsub/itsub
		ufrac			=	utsub/itsub
		
		evfrac			=	np.abs(vfrac)*np.sqrt((noistks[3]/vtsub)**2 + (noistks[0]/itsub)**2)
		eqfrac			=	np.abs(qfrac)*np.sqrt((noistks[1]/qtsub)**2 + (noistks[0]/itsub)**2)
		eufrac			=	np.abs(ufrac)*np.sqrt((noistks[2]/utsub)**2 + (noistks[0]/itsub)**2)
		elfrac			=	np.abs(lfrac)*np.sqrt((elts/lts)**2 + (noistks[0]/itsub)**2)
		epfrac			=	np.abs(pfrac)*np.sqrt((epts/pts)**2 + (noistks[0]/itsub)**2)
				
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
		
		psispec		=	np.rad2deg(0.5*np.arctan2(vspec,lspec))		
		dpsispec	=	np.rad2deg(0.5*np.sqrt((vspec*dlspec)**2 + (lspec*noispec0[2])**2) / (vspec**2 + lspec**2))		
				
		#frbhtr		=	namedtuple('frbhtr',['frbname','dm','nchan','ffac','tavgfac','rm','tresms','twinms','subband','refdeg','tmsarr','fmhzarr',\
		#							 'irms','qrms','urms','vrms',\	
		#							 'it','qt','ut','vt','pt','lt','epts','elts','qfract','ufract','vfract','eqfract','eufract','evfract','lfract','pfract','elfract','epfract',\
		#							 'chit','psit','echit','epsit','chispec','echispec','psispec','epsispec',\
		#							 'ispec','qspec','uspec','vspec','eispec','eqspec','euspec','evspec','lspec','pspec','elspec','epspec',\
		#							 'qfracspec','ufracspec','vfracspec','dqfracspec','dufracspec','dvfracspec','lfracspec','pfracspec','dlfracspec','dpfracspec', \
		#							 'ids','qds','uds','vds','irmspec','qrmspec','urmspec','vrmspec'])
		
		anfrb		=	frbhtr(frbname,dm,nchan,ffac,avgfac,rm0,tresms,tbasems,[subid,nsub],0.0,reltms,fmhzarr,\
									noistks[0],noistks[1],noistks[2],noistks[3],\
									itsub,qtsub,utsub,vtsub,pts,lts,epts,elts,qfrac,ufrac,vfrac,eqfrac,eufrac,evfrac,lfrac,pfrac,elfrac,epfrac,\
									phits,psits,dphits,dpsits,phispec,dphispec,psispec,dpsispec,\
									ispec,qspec,uspec,vspec,noispec0[0],noispec0[1],noispec0[2],noispec0[3],lspec,pspec,dlspec,dpspec,\
									qfracspec,ufracspec,vfracspec,dqfrac,dufrac,dvfrac,lfracspec,pfracspec,dlfrac,dpfrac, \
									stksds[0,:,istart:iend+1],stksds[1,:,istart:iend+1],stksds[2,:,istart:iend+1],stksds[3,:,istart:iend+1],noispec[0],noispec[1],noispec[2],noispec[3])
		
		frbfile		=	open("{}{}_polpkl_{}_to_{}_{}_{}_avg_{}_{}_sb_{}_of_{}_rm_{:.2f}.pkl".format(odatadir,frbname,tbasems[0],tbasems[1],dm,nchan,ffac,avgfac,subid,nsub,rm0),'wb')
		
		pkl.dump(anfrb, frbfile)		
		frbfile.close()
		del(stksds)
		del(noispec)
	else:
		print("Pickling unsuccessful ---- PANIC - File not found!")

	return(tpeak)

#	-------------------------------------------------------------------------------








