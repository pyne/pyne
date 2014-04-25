"""Purpose:

  This module is for spectrometry analysis
  will have functions for general spectrum processing

"""
from warnings import warn
from pyne.utils import VnVWarning

warn(__name__ + " is not yet V&V compliant.", VnVWarning)

class PhSpectrum(object):
    """Pulse height spectrum class"""

    def __init__(self, spec_name="", start_chan_num=0, num_channels=0,
                 channels=None, ebin=None, counts=None):
        """Initialise Ph Spectrum variables"""
        self.channels = [] if channels is None else channels
        self.ebin = [] if ebin is None else ebin
        self.counts = [] if counts is None else counts
        self.spec_name = spec_name
        self.start_chan_num = start_chan_num
        self.num_channels = num_channels
