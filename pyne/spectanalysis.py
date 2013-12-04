"""Purpose:

  This module is for spectrometry analysis
  will have functions for general spectrum processing
   
  .. moduleauthor:: S Lilley

"""

import numpy as np


class PhSpectra():
    """ pulse height spectra class"""

    def __init__(self):
        """Initialise phSpectra variables"""
        self.channels = []
        self.ebin = []
        self.counts = []
        self.spec_name = ""
        self.start_chan_num = 0
        self.num_channels = 0


