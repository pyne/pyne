"""This module is for gamma spectrometry analysis. Currently defines a 
GammaSpectrum class, reads a .spe file Will in the future have functions 
for activity calculations.
"""

from warnings import warn
from pyne.utils import VnVWarning

import numpy as np

from pyne import spectanalysis

warn(__name__ + " is not yet V&V compliant.", VnVWarning)

class GammaSpectrum(spectanalysis.PhSpectrum):
    """GammaSpec class includes gamma specific variables"""

    def __init__(self, real_time=0, live_time=0, det_id='',
                 det_descp='', start_date='', start_time='',
                 calib_e_fit=None, calib_fwhm_fit=None, eres=0, file_name=''):
        """Define gamma spectrum specific variables """
        super(GammaSpectrum, self).__init__()
        self.real_time = real_time
        self.live_time = live_time
        self.dead_time = real_time - live_time
        self.det_id = det_id
        self.det_descp = det_descp
        self.start_date = start_date
        self.start_time = start_time
        self.calib_e_fit = [] if calib_e_fit is None else calib_e_fit
        self.calib_fwhm_fit = [] if calib_fwhm_fit is None else calib_fwhm_fit
        self.eres = eres
        self.file_name = file_name

    def calc_ebins(self):
        """Calculate the energy value for each channel."""
        channels = self.channels = np.asarray(self.channels, float)
        self.ebin = self.calib_e_fit[0] + (self.calib_e_fit[1] * channels) + \
                    (self.calib_e_fit[2] * channels ** 2)

    def __str__(self):
        """Print debug information"""
        print_string = ('Debug print of all header variables\n'
                        'The real time is: {x.real_time}\n'
                        'The live time is: {x.live_time}\n'
                        'The dead time is: {x.dead_time}\n'
                        'Detector ID: {x.det_id}\n'
                        'Detector description: {x.det_descp}\n'
                        'Start date: {x.start_date}\n'
                        'Start time: {x.start_time}\n'
                        'Start channel number: {x.start_chan_num}\n'
                        'Number of channels: {x.num_channels}\n'
                        'Energy calibration fit: {x.calib_e_fit}\n'
                        'FWHM calibration fit: {x.calib_fwhm_fit}\n'
                        'Spectrum: {x.counts}\n'
                        'File name: {x.file_name}').format(x=self)
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

    with open(spec_file_path, 'r') as spec_file:
        full_file_text = spec_file.read()
    file_split = full_file_text.splitlines()
    spec_file.close()
    inspec = False

    # check version of .spe file matches currently supported version
    if (file_split[0] == '$SPEC_ID:'):
        raise RuntimeError('This type of spe file is not supported')

    for item in file_split:
        line = item.split(':')
        # processes the spectrum into 2 lists 1 for channel numbers
        # the other for counts
        if (inspec):
            if (len(line) > 1):
                spectrum.channels.append(int(line[0]))
                temp = line[1].strip()
                spectrum.counts.append(float(temp))

        if (line[0] == 'Spectrum name'):
            spectrum.spec_name = line[1]
        elif (line[0] == 'Detector ID'):
            spectrum.det_id = line[1].strip()
        elif (line[0] == 'Detector description'):
            spectrum.det_descp = line[1].strip()
        elif (line[0] == 'Real Time'):
            spectrum.real_time = float(line[1])
        elif (line[0] == 'Live Time'):
            spectrum.live_time = float(line[1])
        elif (line[0] == 'Acquisition start date'):
            spectrum.start_date = line[1].strip()
        elif (line[0] == 'Acquisition start time'):
            spectrum.start_time = line[1].strip() + ':' + line[2] + ':' + line[3]
        elif (line[0] == 'Starting channel number'):
            spectrum.start_chan_num = int(line[1].strip())
        elif (line[0] == 'Number of channels'):
            spectrum.num_channels = int(line[1].strip())
        elif (line[0] == 'Energy Fit'):
            temp = line[1].strip()
            temp = temp.split(' ')
            spectrum.calib_e_fit.append(float(temp[0]))
            spectrum.calib_e_fit.append(float(temp[2]))
            spectrum.calib_e_fit.append(float(temp[4]))
        elif (line[0] == 'FWHM Fit'):
            temp = line[1].strip()
            temp = temp.split(' ')
            spectrum.calib_fwhm_fit.append(float(temp[0]))
            spectrum.calib_fwhm_fit.append(float(temp[2]))
            spectrum.calib_fwhm_fit.append(float(temp[4]))
        elif (line[0] == 'SPECTRUM'):
            inspec = True

    spectrum.counts = np.array(spectrum.counts)
    spectrum.channels = np.array(spectrum.channels)
    # calculate additional parameters based on .spe file
    spectrum.dead_time = spectrum.real_time - spectrum.live_time
    spectrum.calc_ebins()
    return spectrum


def calc_e_eff(energy, eff_coeff, eff_fit=1):
    """Detector efficiency calculation

    Parameters
    ----------
    energy : float  
        Energy to calcuate det eff 
    eff_coeff : arr
        An array with the coefficients for the energy fit
        the length is not fixed, the length of the array determines the 
        number of terms in the expansion
    eff_fit : int
        Determines what type of fit to use

    Returns
    -------
    eff : float 
        Value of efficiency for the input energy using the selected fitting eqn

    """
    # eff_fit used to choose between calibration fit eqns
    # energy to be in MeV

    if eff_fit == 1:
        # eff_fit 1 uses series ao + a1(lnE)^1+ a2(lnE)^2+ ....
        log_eff = eff_coeff[0]
        i = 1
        while i < len(eff_coeff):
            log_eff = log_eff + (eff_coeff[i] *
                                 (np.power(np.log(energy), i)))
            i = i + 1
        eff = np.exp(log_eff)
    elif eff_fit == 2:
        # eff_fit 2 uses series a0 + a1(1/E)^1 + a2(1/E)^2+...
        log_eff = eff_coeff[0]
        i = 1
        while i < len(eff_coeff):
            log_eff = log_eff + (eff_coeff[i] * ((1 / energy) ** i))
            i = i + 1
        eff = np.exp(log_eff)
    else:
        raise ValueError('The selected eff_fit is not valid')
        eff = 0

    return eff
