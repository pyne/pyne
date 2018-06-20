"""Purpose:
This module is for spectrometry analysis
will have functions for general spectrum processing
"""
from warnings import warn
from pyne.utils import QAWarning


import copy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import splrep, sproot

warn(__name__ + " is not yet QA compliant.", QAWarning)


class PhSpectrum(object):
    """Pulse height spectrum class"""

    def __init__(self, spec_name='', start_chan_num=0, num_channels=0,
                 channels=None, ebin=None, counts=None):
        """Initialise Ph Spectrum variables"""
        self.channels = [] if channels is None else channels
        self.ebin = [] if ebin is None else ebin
        self.counts = [] if counts is None else counts
        self.spec_name = spec_name
        self.start_chan_num = start_chan_num
        self.num_channels = num_channels

def rect_smooth(spectrum, m):
    """Rectangular smoothing function.

    Parameters
    ----------
    spectrum: str
        a spectrum object
    m : int
        the smoothing width, must be an odd integer more than 3

    Returns
    -------
    smooth_spect: a spectrum object

    """

    if(m < 3):
        raise ValueError('Error:Smoothing width less than 3')
    if(m % 2 == 0):
        raise ValueError('Error:Smoothing width not odd')

    smooth_spec = copy.deepcopy(spectrum)
    smooth_spec.counts = [] #reset counts

    ext = int((m - 1.0) / 2.0)

    i = 0
    # 3 stages of loops a small one at start and end to deal
    # with end cases and a main one for bulk of spectrum
    while i < ext:
        smooth_spec.counts.append(spectrum.counts[i])
        i = i + 1
    i = ext
    while i < (len(spectrum.counts) - ext):
        j = i - ext
        sum_m = 0
        while j <= (i + ext):
            sum_m = sum_m + spectrum.counts[j]
            j = j + 1
        smooth_spec.counts.append(sum_m/m)
        i = i + 1
    while i < (len(spectrum.counts)):
        smooth_spec.counts.append(spectrum.counts[i])
        i = i + 1

    smooth_spec.spec_name = spectrum.spec_name + ' smoothed'
    return smooth_spec

def five_point_smooth(spec):
    """5 point smoothing function.

    Recommended for use in low statistics in
    G.W. Phillips , Nucl. Instrum. Methods 153 (1978), 449

    Parameters
    ----------
    spec: a spectrum object

    Returns
    -------
    smooth_spect: a spectrum object

    """
    smooth_spec = copy.deepcopy(spec)
    smooth_spec.counts = []
    smooth_spec.counts.append(spec.counts[0])
    smooth_spec.counts.append(spec.counts[1])
    spec_len = len(spec.counts)
    i = 2
    while i < spec_len - 2:
        val = (1.0 / 9.0) * (spec.counts[i - 2] +
                             spec.counts[i + 2] + (2 * spec.counts[i + 1]) +
                             (2 * spec.counts[i - 1]) + (3 * spec.counts[i]))
        smooth_spec.counts.append(val)
        i = i + 1
    smooth_spec.counts.append(spec.counts[i])
    smooth_spec.counts.append(spec.counts[i + 1])
    smooth_spec.spec_name = spec.spec_name + ' smoothed'
    return smooth_spec

def calc_bg(spec, c1, c2, m):
    """Returns background under a peak"""

    if c1 > c2:
       raise ValueError('c1 must be less than c2')
    if c1 < 0:
       raise ValueError('c1 must be positive number above 0')
    if c2 > max(spec.channels):
       raise ValueError('c2 must be less than max number of channels')

    if m == 1:
        low_sum = sum(spec.counts[c1 - 2:c1])
        high_sum = sum(spec.counts[c2:c2 + 2])
        bg = (low_sum + high_sum) * ((c2 - c1 + 1) / 6)
    else:
        raise ValueError('m is not set to a valud method id')

    return bg

def gross_count(spec, c1, c2):
    """Returns total number of counts in a spectrum between two channels"""

    if c1 > c2:
       raise ValueError('c1 must be less than c2')
    if c1 < 0:
       raise ValueError('c1 must be positive number above 0')
    if c2 > max(spec.channels):
       raise ValueError('c2 must be less than max number of channels')

    gc=sum(spec.counts[c1:c2])
    return gc

def net_counts(spec, c1, c2, m):
    """Calculates net counts between two channels"""
    bg = calc_bg(spec, c1, c2, m)
    gc = gross_count(spec, c1, c2)
    nc = gc - bg
    return nc

def net_area(spec,c1,c2):
    """Calculates the net area under a peak by subtracting background from total
    counts under the peak. 
    
    Parameters
    ----------
    spec : a spectrum object

    c1 : int
        First channel of the peak
    c2 : int
        Second channel of the peak
    Returns
    -------
    net_area: float
        net area under the peak that is not background
    Notes
    ----
    Counts list has to be iterated by the channel number, i.e. counts[300] is
    counts at the 300th channel.
    """
    p = sum(spec.counts[c1:c2])
    n = c2 - c1
    count1 = spec.counts[c1]
    count2 = spec.counts[c2]
    net_area = p - ((n/2.0)*(count1+count2))
    return net_area

def end_point_average_area(spec,c1,c2,var=5):
    """Calculates the net area under a peak by end point averaging the 
    background counts and subtracting the background from total
    counts under the peak. 
    
    Parameters
    ----------
    spec : a spectrum object
     
    c1 : int
        First channel of the peak
    c2 : int
        Second channel of the peak
    var : int
        amount of channels to collect intial and final count averages from.
    
    Returns
    -------
    end_point_average : float
        net area under the peak that is not background counts
    Notes
    ----
    Counts list has to be iterated by the channel number, i.e. counts[300] is
    counts at the 300th channel.
    var is preset to 5
    """
    p = sum(spec.counts[c1:c2])
    n = c2 - c1
    n1 = var
    n2 = var
    count1 = sum(spec.counts[c1-var:c1])
    count2 = sum(spec.counts[c2:c2+var])
    end_point_average = p - ((n/2.0)*((count1/n1)+(count2/n2)))
    return end_point_average

def gaussian_fit(spec,c1,c2):
    """Plots the gaussian fit to the raw data. 
    
    Parameters
    ----------
    spec : a spectrum object
     
    c1 : int
        First channel of the peak
    c2 : int
        Second channel of the peak
    
    Returns
    -------
    gaussian fit : plot or figure
        Gaussian fit of channels versus counts data
    
    """   
    x = spec.channels[c1:c2]
    y = spec.counts[c1:c2]
    n = len(x)
    mean = sum(x*y)/n                   
    sigma = np.sqrt(sum(y*(x-mean)**2)/n)
    amp = max(y)
    x_center = x.mean()
    x_peak = x[list(y).index(max(y))]
    def gauss(x,amp,x_center,sigma):
        return amp*np.exp(-(x-x_center)**2/(2*sigma**2))
    
    popt,pcov = curve_fit(gauss,x,y,p0=[amp,x_center, n/2])
    #po = # for [amp, cen, wid]
    plt.plot(x,y,'b*:',label='data')
    plt.plot(x,gauss(x,*popt),'ro:',label='fit')
    plt.legend()
    plt.title('Gaussian fit' )
    plt.xlabel('Channels')
    plt.ylabel('Counts')
    plt.show()
    
def fwhm(spec,c1,c2,k=3):
    """Gives FWHM of a peak. 
    
    Parameters
    ----------
    spec : a spectrum object
     
    c1 : int
        First channel of the peak
    c2 : int
        Second channel of the peak
    k : int
	The order of the spline fit 
    Returns
    -------
    fwhm : int or float
        FWHM of a peak.
    
    """   
    x = np.array(spec.channels[c1:c2])
    y = np.array(spec.counts[c1:c2])
    n = len(x)                          
    amp = max(y)
    x_center = x.mean()
    def gauss(x,amp,x_center,sigma):
        return amp*np.exp(-(x-x_center)**2/(2.0*sigma**2))  
    popt,pcov = curve_fit(gauss,x,y,p0=[amp,x_center, n/2.0])
    half_max = max(gauss(x,*popt))/2.0
    s = splrep(x, y - half_max, k=k)
    roots = sproot(s)
    fwhm = roots[-1] - roots[0]
    return fwhm

def resolution(fwhm,e_0):
    """Gives resolution of a peak. 
    
    Parameters
    ----------
    fwhm : float
        Full Width Half Maximum of peak
     e_0 : int
        Peak centroid energyt 
    Returns
    -------
    res : float
        Resolution of a peak.
    """   
    res = fwhm/float(e_0)
    return(res)