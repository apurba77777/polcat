import numpy as np

#	--------------------------------------	Useful functions related to dynamic spectra		------------
#
#																Based on scripts by MS
#
#	List of functions
#
#			generate_dynspec(t_ser, n):
#						Creates a dynamic spectrum at the highest time resolution from the given voltage time series.
#
#			calculate_stokes(x, y):
# 						lambda functions for each of the Stokes parameters
#
#			normalise(ds):
#						Normalises a dynamic spectrum along frequency channels to cancel out RFI
#
#			calculate_stokes_unnormalised(x, y):
# 						lambda functions for each of the Stokes parameters (un-normalized)
#
#			movingAverage(x, w):
#						A simple moving average function, series x is smoothed using a boxcar of width w
#
#	----------------------------------------------------------------------------------------------------

def generate_dynspec(t_ser, n):
	#	Creates a dynamic spectrum at the highest time resolution from the given voltage time series.
    #
    #	Inputs		input time series of voltages
    #				number of channels
    #
    #	Returns		dynamic spectrum of voltages
    
    dynspec = np.zeros((int(t_ser.shape[0] / n), n), dtype=np.complex64)
    for i in range(int(t_ser.shape[0] / n)):
        dynspec[i, :] = np.fft.fft(t_ser[i * n : (i + 1) * n])

    return(dynspec)
    
#	------------------------------------------------------------------------------------------

def calculate_stokes_unnormalised(x, y):
    # 	lambda functions for each of the Stokes parameters
    #
    #	Inputs		X & Y voltage time series
    #
    #	Returns		Stokes time series (I, Q, U, V)
        
    stokes = {
        "i": lambda x, y: np.abs(x) ** 2 + np.abs(y) ** 2,
        "q": lambda x, y: np.abs(x) ** 2 - np.abs(y) ** 2,
        "u": lambda x, y: 2 * np.real(np.conj(x) * y),
        "v": lambda x, y: 2 * np.imag(np.conj(x) * y),
    }
    stks = []

    for idx, stk in enumerate(["i", "q", "u", "v"]):
        par = stokes[stk](x, y)
        par_norm = par
        del par
        par = par_norm.transpose()
        del par_norm
        stks += [par]
    
    return(stks)

#	----------------------------------------------------------------------------------------








