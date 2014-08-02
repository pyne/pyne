"""Purpose:

  This module is for spectrometry analysis
  will have functions for general spectrum processing

"""
from warnings import warn
from pyne.utils import VnVWarning

<<<<<<< HEAD
import copy

=======
warn(__name__ + " is not yet V&V compliant.", VnVWarning)
>>>>>>> upstream/staging

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
    a spectrum object,
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
    """ something breaks in this part, the original spectrum is reset as well as the new"""
    smooth_spec = copy.deepcopy(spectrum)
    smooth_spec.counts = [] #reset counts

    ext = (m - 1.0) / 2.0
    
    i = 0
    # 3 stages of loops a small one at start and end to deal 
    # with end cases and a main one for bulk of spectrum
    while i < ext:
        smooth_spec.counts.append(0)
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
    while len(smooth_spec.counts) < (len(spectrum.counts)):
        smooth_spec.counts.append(0)
	
    smooth_spec.spec_name=spectrum.spec_name +' smoothed'
    return smooth_spec

