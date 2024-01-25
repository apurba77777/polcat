#
#	Functions for FRB polarization
#
#								AB, October 2023
#
#	Function list
#
#		plothemall(frbname, dm, nchan, ffac, avgfac, fmhz0, rm0, tbasems, tpeakms, subid, mode, bavg)
#				Plot polarization properties
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
import matplotlib.colors as mpc
from scipy.optimize import curve_fit
from polfns.polmodelfns import *
from polfns.polreduce import *
from polfns.plotprofiles import *
from polfns.plotspec import *

#	----------------------------------------------------------------------------------------------------------

def plothemall(frbname, dm, nchan, ffac, avgfac, fmhz0, rm0, tbasems, tpeakms, subid, mode, bavg):

	#	Plot polarization properties
	
	odatadir	=	htrdir+frbname+'/data/outputs/'
	plotdir		=	htrdir+frbname+'/plots/'
	rmfile		=	htrdir+frbname+'/rm_time.txt'
	compfile	=	htrdir+frbname+'/subands.txt'
	
	avgchan		=	int(round(nchan/ffac))
	
	subfracs	=	np.loadtxt(compfile)
	nsub		=	subfracs.shape[0]-1
	
	frbhtr		=	"{}{}_polpkl_{}_to_{}_{}_{}_avg_{}_{}_sb_{}_of_{}_rm_{:.2f}.pkl".format(odatadir,frbname,tbasems[0],tbasems[1],dm,nchan,ffac,avgfac,subid,nsub,rm0)
	bfrbhtr		=	"{}{}_polpkl_{}_to_{}_{}_{}_avg_{}_{}_sb_{}_of_{}_rm_{:.2f}.pkl".format(odatadir,frbname,tbasems[0],tbasems[1],dm,nchan,ffac,bavg,subid,nsub,rm0)
		
	if(os.path.exists(frbhtr) and os.path.exists(bfrbhtr)):
		
		frbfile		=	open(frbhtr,'rb')
		bfrbfile	=	open(bfrbhtr,'rb')
		
		pdata	=	pkl.load(frbfile)
		p0data	=	pkl.load(bfrbfile)		
		frbfile.close()
		bfrbfile.close()
		
		#	All figure sizes are in cm
		
		if(mode=='iquvt'):			
			plot_iquvt(plotdir,pdata,fmhz0,tbasems,[18.3,17.0],0.2)
			
		elif(mode=='stkst'):			
			plot_iquvtsep(plotdir,pdata,fmhz0,tbasems,[-0.3,0.3],[12.0,6.0],1.0)
			
		elif(mode=='stkspec'):			
			plot_stkspec(plotdir,pdata,[-1.1,1.8],[6.0,4.5])
			
		elif(mode=='polspec'):			
			plot_polspec(plotdir,pdata,[-1.1,1.7],[8.9,8.5],[-90,90,20])
		
		elif(mode=='ilvt'):			
			plot_ilvt(plotdir,pdata,p0data,tbasems,[9.0,6.0],0.2)
		
		elif(mode=='fitcomp'):		
			congaussfitter(plotdir,pdata,tbasems,[8.9,11.0],0.1,4,"I")
			
		elif(mode=='ilvpat'):			
			plot_ilvpat(plotdir,pdata,p0data,tbasems,[6.0,6.0],0.01,10.0)
			
		elif(mode=='ilvparm'):			
			plot_ilvparm(plotdir,rmfile,pdata,p0data,tbasems,[12.1,10.0],0.2,10.0,10.0)
		
		elif(mode=='irm'):			
			plot_irm(plotdir,rmfile,pdata,p0data,tbasems,[18.3,8.5],0.2,5.0,-3.4)
			#plot_irm(plotdir,rmfile,pdata,p0data,tbasems,[8.9,8.0],0.2,5.0,0.0)
		
		elif(mode=='dpa'):			
			plot_dpa(plotdir,rmfile,pdata,p0data,tbasems,[8.9,8.5],0.1,2)
			
		elif(mode=='fracparm'):			
			plot_fracparm(plotdir,rmfile,pdata,p0data,tbasems,[6.0,6.0],0.2,5.0)
			
		elif(mode=='fracpa'):			
			plot_fracpa(plotdir,rmfile,pdata,p0data,tbasems,[6.0,6.0],0.2,5.0)
			
		elif(mode=='haspec'):			
			plot_spechammer(plotdir,pdata,[6.0,4.0])
		
		
	else:
		print("plothemall ---- PANIC - Pickle not found! Have you prepared it?")

	return(0)

#	-------------------------------------------------------------------------------
























