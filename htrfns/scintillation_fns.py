#
#	Functions related to scintillation and scattering
#
#											Based on scripts by MS
#
#	List of functions
#
#				twoScreenDist(tScatt, nuDc, dSO, nu):
#		    				Calculates the maximum product of screen-observer and
#    						screen-source distances in a two scattering screen scenario
#				    
#				autocorr(delta):
#    						to within some factor of normalisation and a 1 cell difference in delay this
#    						does exactly what scipy.signal.correlate(x,x) does    
#
#				autoCovariance(delta):
#							Calculates the auto-covariance function as detailed in Macquart et al 2019
#
#				lorentzian(nulag, nuDc, C):
#							a lorentzian of decorrelation bandwidth nuDc (HWHP)
#
#				fitLorentzian(auto, N, fName, subFreq, subBurst, burstFreq, burst, skippedACFCells, xlenSet = np.arange(20,200,1)):
#							fits a lorentzian to the autocorrrelation of an FRB and produces some plots of the results
#
#				subBandMakeNewFormat(FRB, Fx, Fy, f0, nChan, nSubChannels, xlenSet=np.arange(30,200,1), dcIndex=4.4, skippedACFCells=0):
#						    Constructs an FRB from voltages with nChan total frequency channels, broken into nSubChannels sub channels
#
#				decorrelationStudyNewFormat(FRB, Fx, Fy, f0, nChan, nSubChannels, xlenSet=np.arange(30,200,1), dcIndex=4.4, skippedACFCells=0):
#    						A higher level function that forms the required normalised FRB spectra from voltages then looks for scintillation,
#    						makes some diagnostic plots and saves some of the important data for extended data analysis
#
#				decorrelationStudyPremadeFormat(FRB, stokesI_burst, nChan, nSubChannels, burstStartNu, burstEndNu, f0, xlenSet=np.arange(30,200,1), dcIndex=4.4, skippedACFCells=0):
#							An alternate form of the above decorrStudyNewFormat function which takes in an already normalised 
#    						FRB dynamic spectrum rather than forming one from the unnormalised voltages
#

    
import htrfns.dynamic_fns as frf
import scipy.signal
import scipy.interpolate
import scipy.integrate
import scipy.optimize
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy import constants as const
import glob
import os

#	---------------------------------------------------------------------------------------------------

def twoScreenDist(tScatt, nuDc, dSO, nu):
    """
    Calculates the maximum product of screen-observer and
    screen-source distances in a two scattering screen scenario

       input:
       		tScatt 	= 	scattering time 
       		nuDc	= 	decorrelation bandwidth
       		dSO 	= 	distance between source and observer
       		nu 		= 	central frequency each observable has been measured at 

       output:
       		LxLg 	= 	product of distance between 1st screen and host Lx and 
              			distance between 2nd screen and observer Lg
    """
    
    LxLg = dSO**2/((2*np.pi*nu)**2*tScatt*(1/(2*np.pi*nuDc)))
    
    return(LxLg)
    
#	---------------------------------------------------------------------------------------------------

def autocorr(delta):
    """
    to within some factor of normalisation and a 1 cell difference in delay this
    does exactly what scipy.signal.correlate(x,x) does

    input:
    	delta - thing to autocorrelate

    output:
    	r = autocorrelation 
    """
    
    r = np.zeros(len(delta)*2)
    deltaPad = np.pad(delta, len(delta))
    for j in range(len(r)):
        r[j]=1/(len(delta))*np.sum(delta*deltaPad[j:j+len(delta)])
    
    return(r)

#	---------------------------------------------------------------------------------------------------

def autoCovariance(delta):
    """
    Calculates the auto-covariance function as detailed in Macquart et al 2019

       intput:
       delta = series you want to know the autocovariance of
    
       output:
       r = autocovariance as a function of series cell lag 
    """
    
    r = np.zeros(len(delta)*2)
    deltaPad = np.pad(delta, len(delta))
    for j in range(len(r)):
        r[j]=(np.mean((deltaPad[j:j+len(delta)]-np.mean(delta))*(delta-np.mean(delta))))/(np.mean(delta)**2)
    
    return(r)

#	---------------------------------------------------------------------------------------------------

def lorentzian(nulag, nuDc, C):
    #	a lorentzian of decorrelation bandwidth nuDc (HWHP)
    
    return(C/(1+nulag**2/(nuDc**2)))

#	---------------------------------------------------------------------------------------------------

def fitLorentzian(auto, N, fName, subFreq, subBurst, burstFreq, burst, skippedACFCells, xlenSet = np.arange(20,200,1)):
    """
    fits a lorentzian to the autocorrrelation of an FRB and produces some plots of the results

       input : 
       	auto 			= autocorrelation 
       	N 				= number of channels
       	fName 			= filename
       	subFreq 		= frequencies of hte subband corresponding to the autocorrelation 
       	subBurts 		= burst within the subband
       	burstFreq 		= frequencies of the whole band 
       	burst 			= burst in the whole band
       	skippedACFCells = number of cells skipped in ACF if any 
       	xlenset 		= set of lengths of ACF to be used in fitting process

    """
    
    xdata =336/N*np.arange(-int(auto.shape[0]/2),int(auto.shape[0]/2),1)
    optSetSuperSet = np.zeros([3, len(xlenSet)])
    midpoint = int(xdata.shape[0]/2)
    for i in range(len(xlenSet)):
         optSetTemp = scipy.optimize.curve_fit(lorentzian, xdata[midpoint+skippedACFCells+1:midpoint+1+skippedACFCells+xlenSet[i]], auto[midpoint+skippedACFCells+1:midpoint+1+skippedACFCells+xlenSet[i]], p0 =[1*336/N,np.amax(auto[midpoint+1+skippedACFCells:])], sigma=np.arange(1,len(xdata[midpoint:midpoint+1+xlenSet[i]])))
         optSetSuperSet[:2,i] = np.abs(optSetTemp[0])
         optSetSuperSet[2,i] = np.sum((lorentzian(xdata[midpoint+1+skippedACFCells:midpoint+1+skippedACFCells+xlenSet[i]], optSetTemp[0][0], optSetTemp[0][1])-auto[midpoint+1+skippedACFCells:midpoint+1+skippedACFCells+xlenSet[i]])**2/lorentzian(xdata[midpoint+1+skippedACFCells:midpoint+1+skippedACFCells+xlenSet[i]], optSetTemp[0][0], optSetTemp[0][1]))/xlenSet[i]

    np.save(fName+'optSetSuperSet', optSetSuperSet)
    fig  = plt.figure(tight_layout=True)
    spec = fig.add_gridspec(4,5)
    ax0 = fig.add_subplot(spec[0:2,0:3])
    ax1 = fig.add_subplot(spec[2,0:3])
    ax2 = fig.add_subplot(spec[3,0:3])
    ax3 = fig.add_subplot(spec[0:2,3:])
    ax4 = fig.add_subplot(spec[2:,3:])
    
    ax0.plot(np.arange(xlenSet[len(xlenSet)-1]*1.5), auto[midpoint+1:midpoint+1+len(np.arange(xlenSet[len(xlenSet)-1]*1.5))], color='black')
    ax0.scatter(xlenSet, auto[midpoint+skippedACFCells+1+xlenSet[0]:midpoint+skippedACFCells+2+xlenSet[len(xlenSet)-1]], color='red', label='End points')
    ax0.scatter(np.arange(1,skippedACFCells+1), auto[midpoint+1:midpoint+1+skippedACFCells], color='green', label='skipped')
    ax0.plot(np.arange(xlenSet[len(xlenSet)-1]*1.5), lorentzian(xdata[midpoint+1+skippedACFCells:midpoint+skippedACFCells+len(np.arange(xlenSet[len(xlenSet)-1]*1.5))+1], np.mean(optSetSuperSet[0,:]), np.mean(optSetSuperSet[1,:])), color='C0', label='<best fit>')
    BFIT = optSetSuperSet[2,:]==np.amin(optSetSuperSet[2,:])
    ax0.plot(np.arange(xlenSet[len(xlenSet)-1]*1.5), lorentzian(xdata[midpoint+1+skippedACFCells:midpoint+skippedACFCells+len(np.arange(xlenSet[len(xlenSet)-1]*1.5))+1], optSetSuperSet[0,BFIT], optSetSuperSet[1,BFIT]), color='C1', label='best fit')
      
    ax0.legend()
    ax0.set_xlabel('$\Delta$ cells')
    ax0.set_ylabel('ACF Power')
    ax1.plot(xlenSet, optSetSuperSet[2,:], label='$\chi^2$')
    ax1.set_ylabel('$\chi^2$')
    ax1.set_xlabel('Data length used in fit')
    ax2.plot(xlenSet, optSetSuperSet[0,:], label='$\\nu_{dc}$(MHz)')
    ax2.plot(xlenSet, np.ones(len(xlenSet))*np.mean(optSetSuperSet[0,:]), linestyle='dotted')
    ax2.set_xlabel('Data length used in fit')
    ax2.set_ylabel('$\\nu_{dc}(MHz)$')
    ax3.imshow(subBurst, extent=[0,len(subBurst[0,:]), np.amin(subFreq), np.amax(subFreq)], aspect='auto', interpolation='none')
    ax3.set_xlabel('Time Samples')
    ax3.set_ylabel('frequency (MHz)')
    subBurstIdent = np.repeat(np.expand_dims(((burstFreq>=np.amin(subFreq))*(burstFreq<=np.amax(subFreq))), 1), len(burst[0,:]), axis=1)
    ax4.imshow(burst, extent=[0,len(burst[0,:]), np.amin(burstFreq), np.amax(burstFreq)], aspect='auto', interpolation='none', cmap='gray') 
    ax4.imshow(burst*subBurstIdent, extent=[0,len(burst[0,:]), np.amin(burstFreq), np.amax(burstFreq)], aspect='auto', interpolation='none', cmap='viridis', alpha=0.6) 
    ax4.plot(np.arange(len(burst[0,:])+1), np.ones(len(burst[0,:])+1)*np.amin(subFreq), color='red')
    ax4.plot(np.arange(len(burst[0,:])+1), np.ones(len(burst[0,:])+1)*np.amax(subFreq), color='red')
    ax4.set_xlabel('Time Samples')
    ax4.set_ylabel('frequency (MHz)')
     
    fig.savefig(fName+'lorentzianFitDiagnostic.pdf', bbox_inches='tight')
 
#	--------------------------------------------------------------------------------------------------- 
    
def subBandMakeNewFormat(FRB, Fx, Fy, f0, nChan, nSubChannels, xlenSet=np.arange(30,200,1), dcIndex=4.4, skippedACFCells=0):
    """
    Constructs an FRB from voltages with nChan total frequency channels, broken into nSubChannels sub channels

       input:
       		FRB 			= string of the FRBs name, must correspond to the name of the FRB in burstCutout.cfg
       		Fx 				= shortened voltage time series of x polarisation data
       		Fy 				= shortened voltage time series of y polarisation data
       		nChan 			= number of channels to be formed
       		nSubChannels 	= number of subChannels the FRB should be divivded into for analysis after nChan channels are formed
       		xlenSet 		= set of lengths of ACF to be used when fitting the lorentzian to find scintillation 
       		dcIndex 		= expected frequency evolution power infex of scintillation decorrelation bandwidth size to be 
                 					used in a fiducial plot in the fitLorentzian function 
       		skippedACFCells = number of initial cells in ACF to skip, used to stop fit including spurious structures

       	output :
       		stokesI_burst 			= intensity dynamic spectra of burst, no subbands
       		frequencyRange_burst 	= frequency range of cropped burst dynamic spectra
       		stokesDyn[0] 			= intensity dynamic spectra of whole shortened voltage time series, contains more than burst 
       		frequencyRange 			= frequency range corresponding to stokesDyn[0]
    """
    
    FRB = str(FRB)
    direcs = glob.glob(str(FRB))
    subdirecs = glob.glob(str(FRB)+'/'+str(nSubChannels)+'subChannels')
    if len(direcs)==0:
        os.mkdir(str(FRB)+'/')
        os.mkdir(str(FRB)+'/UNNORMALISED') 
    if len(subdirecs)==0:
        os.mkdir(str(FRB)+'/'+str(nSubChannels)+'subChannels/')
    cutouts = np.loadtxt('burstCutout.cfg', delimiter=',', dtype='str')
    burstStartT = int(float(cutouts[cutouts[:,0]==FRB,1])*336/nChan)-1
    burstEndT = int(float(cutouts[cutouts[:,0]==FRB,2])*336/nChan)+1
    burstStartNu = int(float(cutouts[cutouts[:,0]==FRB,3])*nChan/336)-1
    burstEndNu = int(float(cutouts[cutouts[:,0]==FRB,4])*nChan/336)+1
    if burstStartNu<0:
        burstStartNu=0
    noiseT1 = int(burstStartT-0.2*burstStartT)
    noiseT2 = int(burstEndT+0.2*burstEndT)
    xdyn = frf.generate_dynspec(Fx,n=nChan)
    ydyn = frf.generate_dynspec(Fy,n=nChan)
    stokesDyn = calculate_stokes_unnormalised(xdyn, ydyn)
    frequencyRange = np.flip(np.linspace(f0.value-336/2+1.5, f0.value+336/2+0.5, nChan))
    stokesI_burstUN = stokesDyn[0][burstStartNu:burstEndNu, burstStartT:burstEndT]
    stokesI_noise1UN = stokesDyn[0][burstStartNu:burstEndNu, noiseT1-burstEndT+burstStartT:noiseT1]
    stokesI_noise2UN = stokesDyn[0][burstStartNu:burstEndNu, noiseT2:noiseT2+burstEndT-burstStartT]
    stokesI_burst = ((stokesI_burstUN.T - np.mean(stokesI_noise1UN,1))/np.std(stokesI_noise1UN,1)).T
    stokesI_noise1 = ((stokesI_noise1UN.T - np.mean(stokesI_noise1UN,1))/np.std(stokesI_noise1UN,1)).T
    stokesI_noise2 = ((stokesI_noise2UN.T - np.mean(stokesI_noise2UN,1))/np.std(stokesI_noise2UN,1)).T
    frequencyRange_burst = frequencyRange[burstStartNu:burstEndNu]
    np.save(str(FRB)+'/'+str(nSubChannels)+'subChannels/'+'stokesI_burst', stokesI_burst)
    np.save(str(FRB)+'/'+str(nSubChannels)+'subChannels/'+'stokesI_noise1', stokesI_noise1)
    np.save(str(FRB)+'/'+str(nSubChannels)+'subChannels/'+'stokesI_noise2', stokesI_noise2)
    np.save(str(FRB)+'/UNNORMALISED/stokesI_burst', stokesI_burstUN)
    np.save(str(FRB)+'/UNNORMALISED/stokesI_noise1', stokesI_noise1UN)
    np.save(str(FRB)+'/UNNORMALISED/stokesI_noise2', stokesI_noise2UN)
    if (len(stokesI_burst)/nSubChannels==int(len(stokesI_burst)/nSubChannels)):
        subBursts = np.zeros([int(len(stokesI_burst[:,0])/nSubChannels), len(stokesI_burst[0,:]), nSubChannels])
        subNoises = np.zeros([int(len(stokesI_noise1[:,0])/nSubChannels), len(stokesI_noise1[0,:]), nSubChannels, 2])
        subFreqRange = np.zeros([int(len(frequencyRange_burst)/nSubChannels), nSubChannels])
        for i in range(nSubChannels):
            subFreqRange[:,i]=frequencyRange_burst[int(i*len(frequencyRange_burst)/nSubChannels):int((i+1)*len(frequencyRange_burst)/nSubChannels)]
            subBursts[:, :, i]= stokesI_burst[int(i*len(stokesI_burst[:,0])/nSubChannels):int((i+1)*len(stokesI_burst[:,0])/nSubChannels), :]
            subNoises[:, :, i, 0]= stokesI_noise1[int(i*len(stokesI_noise1[:,0])/nSubChannels):int((i+1)*len(stokesI_noise1[:,0])/nSubChannels), :]
            subNoises[:, :, i, 1]= stokesI_noise2[int(i*len(stokesI_noise2[:,0])/nSubChannels):int((i+1)*len(stokesI_noise2[:,0])/nSubChannels), :]
            np.save(str(FRB)+'/'+str(nSubChannels)+'subChannels/'+'subFreqRange', subFreqRange)
            np.save(str(FRB)+'/'+str(nSubChannels)+'subChannels/'+'subBursts', subBursts)
            np.save(str(FRB)+'/'+str(nSubChannels)+'subChannels/'+'subNoises', subNoises)
    else:
        print('number of subchannels not consistent with number of channels, try again :', nSubChannels, len(stokesI_burst))
        subBursts = 0
        subNoises = 0
    
    return(stokesI_burst, frequencyRange_burst, stokesDyn[0], frequencyRange)    

#	---------------------------------------------------------------------------------------------------

def decorrelationStudyNewFormat(FRB, Fx, Fy, f0, nChan, nSubChannels, xlenSet=np.arange(30,200,1), dcIndex=4.4, skippedACFCells=0):
    """
    A higher level function that forms the required normalised FRB spectra from voltages then looks for scintillation,
    makes some diagnostic plots and saves some of the important data for extended data analysis

       input: 
       		FRB 			= string of the FRBs name, must correspond to the name of the FRB in burstCutout.cfg
       		Fx 				= shortened voltage time series of x polarisation data
       		Fy 				= shortened voltage time series of y polarisation data
       		nChan 			= number of channels to be formed
       		nSubChannels 	= number of subChannels the FRB should be divivded into for analysis after nChan channels are formed
       		xlenSet 		= set of lengths of ACF to be used when fitting the lorentzian to find scintillation 
       		dcIndex 		= expected frequency evolution power infex of scintillation decorrelation bandwidth size to be 
                 				used in a fiducial plot in the fitLorentzian function 
       		skippedACFCells = number of initial cells in ACF to skip, used to stop fit including spurious structures
 	"""
    
    cutouts = np.loadtxt('burstCutout.cfg', delimiter=',', dtype='str')
    burstStartT = int(float(cutouts[cutouts[:,0]==FRB,1])*336/nChan)-1
    burstEndT = int(float(cutouts[cutouts[:,0]==FRB,2])*336/nChan)+1
    burstStartNu = int(float(cutouts[cutouts[:,0]==FRB,3])*nChan/336)-1
    burstEndNu = int(float(cutouts[cutouts[:,0]==FRB,4])*nChan/336)+1
    if burstStartNu<0:
        burstStartNu=0
    noiseT1 = int(burstStartT-0.2*burstStartT)
    noiseT2 = int(burstEndT+0.2*burstEndT)
    stokesI_burst, frequencyRange_burst, stokesI, frequencyRange = subBandMakeNewFormat(FRB, Fx, Fy, f0, nChan, nSubChannels)
    subFreqRange = np.load(str(FRB)+'/'+str(nSubChannels)+'subChannels/'+'subFreqRange.npy')
    subBursts = np.load(str(FRB)+'/'+str(nSubChannels)+'subChannels/'+'subBursts.npy')
    subNoises = np.load(str(FRB)+'/'+str(nSubChannels)+'subChannels/'+'subNoises.npy')
    autos = np.zeros([int(len(frequencyRange_burst)*2/nSubChannels), nSubChannels])
    autosNoise = np.zeros([int(len(frequencyRange_burst)*2/nSubChannels), nSubChannels, 2])
    for i in range(nSubChannels):
        print(subBursts.shape)
        autos[:,i]= autocorr((np.sum(subBursts[:,:,i],1))-(np.sum(subBursts[:,:,i],1)).mean())
        np.save(str(FRB)+'/'+str(nSubChannels)+'subChannels/'+'autos', autos)
        autosNoise[:,i,0]= autocorr((np.sum(subNoises[:,:,i,0],1))-(np.sum(subNoises[:,:,i,0],1)).mean())
        autosNoise[:,i,1]= autocorr((np.sum(subNoises[:,:,i,1],1))-(np.sum(subNoises[:,:,i,1],1)).mean())
        fitLorentzian(autos[:,i], nChan, str(FRB)+'/'+str(nSubChannels)+'subChannels/'+'dCIndex'+str(dcIndex)+'_'+str(i+1)+'of'+str(nSubChannels)+'_'+str(FRB), subFreqRange[:,i], subBursts[:,:,i], frequencyRange_burst, stokesI_burst, skippedACFCells, xlenSet=xlenSet)


        fig,ax = plt.subplots(2)
        ax[0].plot((np.linspace(0, int(len(autos[:,i])/2), int(len(autos[:,i])/2))[1:])*336/nChan,autos[int(len(autos[:,i])/2)+1:,i], label='burst')
        ax[0].plot((np.linspace(0,int(len(autos[:,i])/2) , int(len(autosNoise[:,i,0])/2))[1:])*336/nChan,autosNoise[int(len(autosNoise[:,i,0])/2)+1:,i,0], label='Pre-burst noise')
        ax[0].plot((np.linspace(0,int(len(autos[:,i])/2) , int(len(autosNoise[:,i,1])/2))[1:])*336/nChan,autosNoise[int(len(autosNoise[:,i,1])/2)+1:,i,1], label='Post-burst noise')
        ax[0].set_xlabel('lag (MHz)')
        ax[0].set_ylabel('ACF (arbitrary units)')
        ax[0].legend()
        subburstStartNu = np.amin(subFreqRange[:,i])
        subburstEndNu = np.amax(subFreqRange[:,i])
        ax[1].imshow(stokesI, extent=[0,len(stokesI[0,:]), np.amin(frequencyRange), np.amax(frequencyRange)], interpolation='none', aspect='auto')
        ax[1].plot(np.linspace(burstStartT, burstEndT, 10), np.ones(10)*subburstStartNu, color='C0')
        ax[1].plot(np.linspace(burstStartT, burstEndT, 10), np.ones(10)*subburstEndNu, color='C0')
        ax[1].plot(np.ones(10)*burstStartT, np.linspace(subburstStartNu,subburstEndNu,10), color='C0')
        ax[1].plot(np.ones(10)*burstEndT, np.linspace(subburstStartNu,subburstEndNu,10), color='C0')
        ax[1].plot(np.linspace(noiseT1-(burstEndT-burstStartT), noiseT1, 10), np.ones(10)*subburstStartNu, color='C1')
        ax[1].plot(np.linspace(noiseT1-(burstEndT-burstStartT), noiseT1, 10), np.ones(10)*subburstEndNu, color='C1')
        ax[1].plot(np.ones(10)*(noiseT1-(burstEndT-burstStartT)), np.linspace(subburstStartNu,subburstEndNu,10), color='C1')
        ax[1].plot(np.ones(10)*(noiseT1), np.linspace(subburstStartNu,subburstEndNu,10), color='C1')
        ax[1].plot(np.linspace(noiseT2, noiseT2+(burstEndT-burstStartT), 10), np.ones(10)*subburstStartNu, color='C2')
        ax[1].plot(np.linspace(noiseT2, noiseT2+(burstEndT-burstStartT), 10), np.ones(10)*subburstEndNu, color='C2')
        ax[1].plot(np.ones(10)*(noiseT2), np.linspace(subburstStartNu,subburstEndNu,10), color='C2')
        ax[1].plot(np.ones(10)*(noiseT2+(burstEndT-burstStartT)), np.linspace(subburstStartNu,subburstEndNu,10), color='C2')

        fig.savefig(str(FRB)+'/'+str(nSubChannels)+'subChannels/'+'burst-noiseACFComp'+str(FRB)+'dCIndex'+str(dcIndex)+'_'+str(i)+'of'+str(nSubChannels)+'.pdf', bbox_inches='tight')
    optSetSuperSets = np.zeros([nSubChannels, 3, len(xlenSet)])
    for i in range(nSubChannels):
        optSetSuperSets[i,:,:] = np.load(str(FRB)+'/'+str(nSubChannels)+'subChannels/'+'dCIndex'+str(dcIndex)+'_'+str(i+1)+'of'+str(nSubChannels)+'_'+str(FRB)+'optSetSuperSet.npy')
    fig, ax = plt.subplots(2, tight_layout=True)
    for i in range(nSubChannels):
        ax[1].plot(xlenSet*336/nChan, optSetSuperSets[i,2,:], label='SC = '+str(i+1))
    for i in range(nSubChannels-1):
        ax[0].plot(xlenSet*336/nChan, optSetSuperSets[0,0,:]/optSetSuperSets[i+1,0,:], label='SC = '+str(i+1+1), color='C'+str(i))
        ax[0].plot(xlenSet*336/nChan, np.ones(len(optSetSuperSets[0,0,:]))*(np.mean(subFreqRange[:,0])/np.mean(subFreqRange[:,i+1]))**dcIndex, color='C'+str(i), linestyle='dashed')

    ax[0].set_xlabel('data length (MHz)')
    ax[1].set_ylabel('$\chi^2$')
    ax[0].set_ylabel('$\\nu_{dc,SC0}/\\nu_{dc,SCX}$')
    ax[0].legend()
    ax[1].legend()
    fig.savefig(str(FRB)+'/'+str(nSubChannels)+'subChannels/'+'nuDCExpectationDiagnostic'+str(FRB)+'dCIndex'+str(dcIndex)+'_'+str(nSubChannels)+'.pdf', bbox_inches='tight')

 

    fig, ax = plt.subplots(1)
    nudcFitSet = np.zeros([nSubChannels, 2])
    HPPFitSet = np.zeros([nSubChannels])
    for i in range(nSubChannels):
        nudcFitSet[i,0] = np.log10(np.mean(subFreqRange[:,i]))
        nudcFitSet[i,1] = np.log10(np.mean(optSetSuperSets[i,0,:]))
        ax.scatter(np.log10(np.mean(subFreqRange[:,i])), np.log10(np.mean(optSetSuperSets[i,0,:])), color='C'+str(i))
        interpTEMP = scipy.interpolate.interp1d(autos[int(len(autos[:,i])/2)+1:,i], (np.linspace(0, int(len(autos[:,i])/2), int(len(autos[:,i])/2))[1:])*336/nChan)
        HPPFitSet[i] = np.log10(interpTEMP(np.amax(autos[int(len(autos[:,i])/2)+1:,i])/2.0))
        ax.scatter(nudcFitSet[i,0], HPPFitSet[i], color='C'+str(i), marker='s')
    regSet = scipy.stats.linregress(nudcFitSet[:,0], nudcFitSet[:,1])
    showfitX = np.linspace(np.amin(nudcFitSet[:,0]), np.amax(nudcFitSet[:,0]),100)
    showfitY = regSet[0]*showfitX+regSet[1]
    ax.plot(showfitX, showfitY, label='Index = '+str(regSet[0]))
    regSetHPP = scipy.stats.linregress(nudcFitSet[:,0], HPPFitSet)
    showfitYHPP = regSetHPP[0]*showfitX+regSetHPP[1]
    ax.plot(showfitX, showfitYHPP, label='Index = '+str(regSetHPP[0]), linestyle='dashed')
    ax.legend()
    ax.set_xlabel('$\log_{10}\\nu$(MHz)')
    ax.set_ylabel('$\log_{10}\\nu_{dc}$(MHz)')
    fig.savefig(str(FRB)+'/'+str(nSubChannels)+'subChannels/'+'nuDCvsCentralFreqsAVG'+str(FRB)+'dCIndex'+str(dcIndex)+'_'+str(nSubChannels)+'.pdf', bbox_inches='tight')
    
    fig, ax = plt.subplots(1)
    nudcFitSet = np.zeros([nSubChannels, 2])
    for i in range(nSubChannels):
        nudcFitSet[i,0] = np.log10(np.mean(subFreqRange[:,i]))
        nudcFitSet[i,1] = np.log10(optSetSuperSets[i, 0, optSetSuperSets[i,2,:]==np.amin(optSetSuperSets[i,2,:])])
        ax.scatter(nudcFitSet[i,0], nudcFitSet[i,1])
        ax.scatter(nudcFitSet[i,0], HPPFitSet[i], color='C'+str(i), marker='s')
    regSet = scipy.stats.linregress(nudcFitSet[:,0], nudcFitSet[:,1])
    showfitX = np.linspace(np.amin(nudcFitSet[:,0]), np.amax(nudcFitSet[:,0]),100)
    showfitY = regSet[0]*showfitX+regSet[1]
    ax.plot(showfitX, showfitY, label='Index = '+str(regSet[0]))
    ax.plot(showfitX, showfitYHPP, label='Index = '+str(regSetHPP[0]), linestyle='dashed')
    ax.legend()
    ax.set_xlabel('$\log_{10}\\nu$')
    ax.set_ylabel('$\log_{10}\\nu_{dc}$')
    fig.savefig(str(FRB)+'/'+str(nSubChannels)+'subChannels/'+'nuDCvsCentralFreqsMIN'+str(FRB)+'dCIndex'+str(dcIndex)+'_'+str(nSubChannels)+'.pdf', bbox_inches='tight')

#	---------------------------------------------------------------------------------------------------

def decorrelationStudyPremadeFormat(FRB, stokesI_burst, nChan, nSubChannels, burstStartNu, burstEndNu, f0, xlenSet=np.arange(30,200,1), dcIndex=4.4, skippedACFCells=0):
    """
    An alternate form of the above decorrStudyNewFormat function which takes in an already normalised 
    FRB dynamic spectrum rather than forming one from the unnormalised voltages

       input:
       		FRB 			= string of the FRBs name, must correspond to the name of the FRB in burstCutout.cfg
       		stokesI_burst 	= normalised FRB dynamic spectrum
       		nChan 			= number of channels in stokesI_burst
       		nSubChannels 	= number of subChannels the FRB should be divivded into for analysis after nChan channels are formed
       		burstStartNu 	= starting frequency channel of stokesI_burst, cell number that the stokesI_burst was cropped from
       		burstEndNu 		= ending frequency channel of stokesI_burst, cell number that the stokesI_burst was cropped from
       		f0 				= central frequency of stokesI_burst observation, before cropping to isolate burst
       		xlenSet 		= set of lengths of ACF to be used when fitting the lorentzian to find scintillation 
       		dcIndex 		= expected frequency evolution power infex of scintillation decorrelation bandwidth size to be 
                 				used in a fiducial plot in the fitLorentzian function 
       		skippedACFCells = number of initial cells in ACF to skip, used to stop fit including spurious structures

       output:
    """
    
    FRB = str(FRB)
    direcs = glob.glob(str(FRB))
    subdirecs = glob.glob(str(FRB)+'/'+str(nSubChannels)+'subChannels')
    if len(direcs)==0:
        os.mkdir(str(FRB)+'/')
    if len(subdirecs)==0:
        os.mkdir(str(FRB)+'/'+str(nSubChannels)+'subChannels/')
    frequencyRange = np.flip(np.linspace(f0.value-336/2+1.5, f0.value+336/2+0.5, nChan))
    
    burstStartT = 0
    burstEndT = stokesI_burst.shape[1]
    frequencyRange_burst = frequencyRange[burstStartNu:burstEndNu]
    if (len(stokesI_burst)/nSubChannels==int(len(stokesI_burst)/nSubChannels)):
        subBursts = np.zeros([int(len(stokesI_burst[:,0])/nSubChannels), len(stokesI_burst[0,:]), nSubChannels])
        subFreqRange = np.zeros([int(len(frequencyRange_burst)/nSubChannels), nSubChannels])
        autos = np.zeros([int(len(frequencyRange_burst)*2/nSubChannels), nSubChannels])
        for i in range(nSubChannels):
            subFreqRange[:,i]=frequencyRange_burst[int(i*len(frequencyRange_burst)/nSubChannels):int((i+1)*len(frequencyRange_burst)/nSubChannels)]
            subBursts[:, :, i]= stokesI_burst[int(i*len(stokesI_burst[:,0])/nSubChannels):int((i+1)*len(stokesI_burst[:,0])/nSubChannels), :]
            autos[:,i]= autocorr((np.sum(subBursts[:,:,i],1))-(np.sum(subBursts[:,:,i],1)).mean())
            fitLorentzian(autos[:,i], nChan, str(FRB)+'/'+str(nSubChannels)+'subChannels/'+'dCIndex'+str(dcIndex)+'_'+str(i+1)+'of'+str(nSubChannels)+'_'+str(FRB), subFreqRange[:,i], subBursts[:,:,i], frequencyRange_burst, stokesI_burst, skippedACFCells, xlenSet=xlenSet)


            fig,ax = plt.subplots(2)
            ax[0].plot((np.linspace(0, int(len(autos[:,i])/2), int(len(autos[:,i])/2))[1:])*336/nChan,autos[int(len(autos[:,i])/2)+1:,i], label='burst')
            ax[0].set_xlabel('lag (MHz)')
            ax[0].set_ylabel('ACF (arbitrary units)')
            ax[0].legend()
            subburstStartNu = np.amin(subFreqRange[:,i])
            subburstEndNu = np.amax(subFreqRange[:,i])
            ax[1].imshow(stokesI_burst, extent=[0,len(stokesI_burst[0,:]), np.amin(frequencyRange), np.amax(frequencyRange)], interpolation='none', aspect='auto')
            ax[1].plot(np.linspace(burstStartT, burstEndT, 10), np.ones(10)*subburstStartNu, color='C0')
            ax[1].plot(np.linspace(burstStartT, burstEndT, 10), np.ones(10)*subburstEndNu, color='C0')
            ax[1].plot(np.ones(10)*burstStartT, np.linspace(subburstStartNu,subburstEndNu,10), color='C0')
            ax[1].plot(np.ones(10)*burstEndT, np.linspace(subburstStartNu,subburstEndNu,10), color='C0')

            fig.savefig(str(FRB)+'/'+str(nSubChannels)+'subChannels/'+'burst-noiseACFComp'+str(FRB)+'dCIndex'+str(dcIndex)+'_'+str(i)+'of'+str(nSubChannels)+'.pdf', bbox_inches='tight')
        optSetSuperSets = np.zeros([nSubChannels, 3, len(xlenSet)])
        for i in range(nSubChannels):
            optSetSuperSets[i,:,:] = np.load(str(FRB)+'/'+str(nSubChannels)+'subChannels/'+'dCIndex'+str(dcIndex)+'_'+str(i+1)+'of'+str(nSubChannels)+'_'+str(FRB)+'optSetSuperSet.npy')
        fig, ax = plt.subplots(2, tight_layout=True)
        for i in range(nSubChannels):
            ax[1].plot(xlenSet*336/nChan, optSetSuperSets[i,2,:], label='SC = '+str(i+1))
        for i in range(nSubChannels-1):
            ax[0].plot(xlenSet*336/nChan, optSetSuperSets[0,0,:]/optSetSuperSets[i+1,0,:], label='SC = '+str(i+1+1), color='C'+str(i))
            ax[0].plot(xlenSet*336/nChan, np.ones(len(optSetSuperSets[0,0,:]))*(np.mean(subFreqRange[:,0])/np.mean(subFreqRange[:,i+1]))**dcIndex, color='C'+str(i), linestyle='dashed')

        ax[0].set_xlabel('data length (MHz)')
        ax[1].set_ylabel('$\chi^2$')
        ax[0].set_ylabel('$\\nu_{dc,SC0}/\\nu_{dc,SCX}$')
        ax[0].legend()
        ax[1].legend()
        fig.savefig(str(FRB)+'/'+str(nSubChannels)+'subChannels/'+'nuDCExpectationDiagnostic'+str(FRB)+'dCIndex'+str(dcIndex)+'_'+str(nSubChannels)+'.pdf', bbox_inches='tight')

 

        fig, ax = plt.subplots(1)
        nudcFitSet = np.zeros([nSubChannels, 2])
        HPPFitSet = np.zeros([nSubChannels])
        for i in range(nSubChannels):
            nudcFitSet[i,0] = np.log10(np.mean(subFreqRange[:,i]))
            nudcFitSet[i,1] = np.log10(np.mean(optSetSuperSets[i,0,:]))
            ax.scatter(np.log10(np.mean(subFreqRange[:,i])), np.log10(np.mean(optSetSuperSets[i,0,:])), color='C'+str(i))
            interpTEMP = scipy.interpolate.interp1d(autos[int(len(autos[:,i])/2)+1:,i], (np.linspace(0, int(len(autos[:,i])/2), int(len(autos[:,i])/2))[1:])*336/nChan)
            HPPFitSet[i] = np.log10(interpTEMP(np.amax(autos[int(len(autos[:,i])/2)+1:,i])/2.0))
            ax.scatter(nudcFitSet[i,0], HPPFitSet[i], color='C'+str(i), marker='s')
        regSet = scipy.stats.linregress(nudcFitSet[:,0], nudcFitSet[:,1])
        showfitX = np.linspace(np.amin(nudcFitSet[:,0]), np.amax(nudcFitSet[:,0]),100)
        showfitY = regSet[0]*showfitX+regSet[1]
        ax.plot(showfitX, showfitY, label='Index = '+str(regSet[0]))
        regSetHPP = scipy.stats.linregress(nudcFitSet[:,0], HPPFitSet)
        showfitYHPP = regSetHPP[0]*showfitX+regSetHPP[1]
        ax.plot(showfitX, showfitYHPP, label='Index = '+str(regSetHPP[0]), linestyle='dashed')
        ax.legend()
        ax.set_xlabel('$\log_{10}\\nu$(MHz)')
        ax.set_ylabel('$\log_{10}\\nu_{dc}$(MHz)')
        fig.savefig(str(FRB)+'/'+str(nSubChannels)+'subChannels/'+'nuDCvsCentralFreqsAVG'+str(FRB)+'dCIndex'+str(dcIndex)+'_'+str(nSubChannels)+'.pdf', bbox_inches='tight')
    
        fig, ax = plt.subplots(1)
        nudcFitSet = np.zeros([nSubChannels, 2])
        for i in range(nSubChannels):
            nudcFitSet[i,0] = np.log10(np.mean(subFreqRange[:,i]))
            nudcFitSet[i,1] = np.log10(optSetSuperSets[i, 0, optSetSuperSets[i,2,:]==np.amin(optSetSuperSets[i,2,:])])
            ax.scatter(nudcFitSet[i,0], nudcFitSet[i,1])
            ax.scatter(nudcFitSet[i,0], HPPFitSet[i], color='C'+str(i), marker='s')
        regSet = scipy.stats.linregress(nudcFitSet[:,0], nudcFitSet[:,1])
        showfitX = np.linspace(np.amin(nudcFitSet[:,0]), np.amax(nudcFitSet[:,0]),100)
        showfitY = regSet[0]*showfitX+regSet[1]
        ax.plot(showfitX, showfitY, label='Index = '+str(regSet[0]))
        ax.plot(showfitX, showfitYHPP, label='Index = '+str(regSetHPP[0]), linestyle='dashed')
        ax.legend()
        ax.set_xlabel('$\log_{10}\\nu$')
        ax.set_ylabel('$\log_{10}\\nu_{dc}$')
        fig.savefig(str(FRB)+'/'+str(nSubChannels)+'subChannels/'+'nuDCvsCentralFreqsMIN'+str(FRB)+'dCIndex'+str(dcIndex)+'_'+str(nSubChannels)+'.pdf', bbox_inches='tight')

    else:
        print('number of subchannels not consistent with number of channels, try again :', nSubChannels, len(stokesI_burst))
        
#	---------------------------------------------------------------------------------------------------








