"""Purpose:

  This module is for gamma spectrometry analysis
  Currently defines a GammaSpectrum, reads a .spe file
  Will in the future have functions for activity calculations

"""

import numpy as np
import spectanalysis


class GammaSpectrum(spectanalysis.PhSpectrum):
    """ GammaSpec class includes gamma specific variables"""

    def __init__(self, real_time = 0, live_time = 0, det_id = "", 
    det_descp = "", start_date = "", start_time = "", 
    calib_e_fit=None, calib_fwhm_fit=None, eres = 0, file_name = ""):
        """ define gamma spectrum specific variables """    
        super(GammaSpectrum, self).__init__()
        self.real_time = real_time
        self.live_time = live_time
        self.dead_time = real_time-live_time
        self.det_id = det_id
        self.det_descp = det_descp
        self.start_date = start_date
        self.start_time = start_time
        self.calib_e_fit = [] if calib_e_fit is None else calib_e_fit
        self.calib_fwhm_fit = [] if calib_fwhm_fit is None else calib_fwhm_fit
        self.eres = eres
        self.file_name = file_name

    def calc_ebins(self):
        """ Calculate the energy value for each channel."""
        channels = self.channels = np.asarray(self.channels, float)
        self.ebin = self.calib_e_fit[0] + (self.calib_e_fit[1]*channels) + \
        (self.calib_e_fit[2]*channels**2)
        
    def __str__(self):
        """ print debug information"""
        print_string= "Debug print of all header variables" \
        + "\n" + "The real time is:" + str(self.real_time) \
        + "\n" + "The live time is:" + str(self.live_time) \
        + "\n" + "The dead time is:" + str(self.dead_time) \
        + "\n" + "Detector ID:" + self.det_id \
        + "\n" + "Detector description:" + self.det_descp \
        + "\n" + "Start date:" + self.start_date \
        + "\n" + "Start time:" + self.start_time \
        + "\n" + "Start channel number:" + str(self.start_chan_num) \
        + "\n" + "Number of channels:" + str(self.num_channels) \
        + "\n" + "Energy calibration fit:" + str(self.calib_e_fit) \
        + "\n" + "FWHM calibration fit:" + str(self.calib_fwhm_fit) \
        + "\n" + "Spectrum:" + str(self.counts) \
        + "\n" + "File name:" + self.file_name
        return print_string


def read_spe_file(spec_file_path):
    """Reads a .spe file

    Parameters
    ----------
    spec_file_path : str 
        Path to spe file

    Returns
    -------
    spectrum : a GammaSpec object
        Contains all information from .spe file
        
    """

    spectrum = GammaSpectrum()
    spectrum.file_name = spec_file_path
   
    with open(spec_file_path, "r") as spec_file:
        full_file_text = spec_file.read()
    file_split = full_file_text.split("\n")
    spec_file.close()
    inspec = False
    
    # check version of .spe file matches currently supported version
    if (file_split[0] == '$SPEC_ID:'):
        raise RuntimeError('This type of spe file is not supported')

    for item in file_split:
        line = item.split(":")
        # processes the spectrum into 2 lists 1 for channel numbers
        # the other for counts
        if (inspec):
            if (len(line) > 1):
                spectrum.channels.append(int(line[0]))
                temp = line[1].strip()
                spectrum.counts.append(float(temp))

        if (line[0] == "Spectrum name"):
            spectrum.spec_name = line[1]
        elif (line[0] == "Detector ID"):
            spectrum.det_id = line[1].strip()
        elif (line[0] == "Detector description"):
            spectrum.det_descp = line[1].strip()
        elif (line[0] == "Real Time"):
            spectrum.real_time = float(line[1])
        elif (line[0] == "Live Time"):
            spectrum.live_time = float(line[1])
        elif (line[0] == "Acquisition start date"):
            spectrum.start_date = line[1].strip()
        elif (line[0] == "Acquisition start time"):
            spectrum.start_time = line[1].strip() + ":" + line[2] + ":"+line[3]
        elif (line[0] == "Starting channel number"):
            spectrum.start_chan_num = int(line[1].strip())
        elif (line[0] == "Number of channels"):
            spectrum.num_channels = int(line[1].strip())
        elif (line[0] == "Energy Fit"):
            temp = line[1].strip()
            temp = temp.split(" ")
            spectrum.calib_e_fit.append(float(temp[0]))
            spectrum.calib_e_fit.append(float(temp[2]))
            spectrum.calib_e_fit.append(float(temp[4]))
        elif (line[0] == "FWHM Fit"):
            temp = line[1].strip()
            temp = temp.split(" ")
            spectrum.calib_fwhm_fit.append(float(temp[0]))
            spectrum.calib_fwhm_fit.append(float(temp[2]))
            spectrum.calib_fwhm_fit.append(float(temp[4]))
        elif (line[0] == "SPECTRUM"):
            inspec = True

    spectrum.counts = np.array(spectrum.counts)
    spectrum.channels = np.array(spectrum.channels)
    # calculate additional parameters based on .spe file
    spectrum.dead_time=spectrum.real_time-spectrum.live_time
    spectrum.calc_ebins()
    return spectrum
