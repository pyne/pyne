"""
    This module is for gamma spectrometry analysis
    November 2013

"""

# import useful modules
import numpy as np
# import local modules
import spectAnalysis


class gammaSpec(spectAnalysis.phSpectra):
    """ gamma Spec class includes gamma specific functionality

        inherits from PHspectra class in spectAnalysis

    """

    def __init__(self):
        """ define gamma spectra specific variables """    
        spectAnalysis.phSpectra.__init__(self)
        self.real_time = 0
        self.live_time = 0
        self.dead_time = 0
        self.det_ID = ""
        self.det_descp = ""
        self.start_date = ""
        self.start_time = ""
        self.calib_E_fit = []
        self.calib_FWHM_fit = []
        self.Eres = 0
        self.file_name = ""

    def calc_dead_time(self):
        """Calculate the dead time. """
        self.dead_time = self.real_time - self.live_time

    def calc_ebins(self):
        """ Calculate the energy value for each channel."""
        for item in self.channels:
            energy = self.calib_E_fit[0] + (self.calib_E_fit[1]*item)
            + (self.calib_E_fit[2]*item*item)
            self.Ebin.append(energy)
        self.Ebin = np.array(self.Ebin)

    def print_parameters(self):
        """Print all the gamma spectrum parameters"""
        print "Debug print of all header variables"
        print "The real time is:" + str(self.real_time)
        print "The live time is:" + str(self.live_time)
        print "The dead time is:" + str(self.dead_time)
        print "Detector ID:" + self.det_ID
        print "Detector description:" + self.det_descp
        print "Start date:" + self.start_date
        print "Start time:" + self.start_time
        print "Start channel number:" + str(self.start_chan_num)
        print "Number of channels:" + str(self.num_channels)
        print "Energy calibration fit:" + str(self.calib_E_fit)
        print "FWHM calibration fit:" + str(self.calib_FWHM_fit)
        print "Spectrum:" + str(self.counts)
        print "File name:" + self.file_name


def readspefile(spec_file_path):
    """ Read spe file

    Reads a gamma vision .spe file and assigns all the useful
    parts to the correct variables.
    arguments: a file path to read in
    returns: gamma spectra object

    """

    print "Starting to read file:" + spec_file_path
    spectra = gammaSpec()
    spectra.file_name = spec_file_path
   
    spec_file = open(spec_file_path, "r")
    full_file_text = spec_file.read()
    file_split = full_file_text.split("\n")
    spec_file.close()
    inspec = False

    for item in file_split:
        line = item.split(":")
        # it processes the spectrum into 2 lists 1 for channel numbers
        # the other for counts
        if(inspec):
            if (len(line) > 1):
                spectra.channels.append(int(line[0]))
                temp = line[1].strip()
                spectra.counts.append(float(temp))

        if(line[0] == "Spectrum name"):
            spectra.spec_name = line[1]
        elif(line[0] == "Detector ID"):
            spectra.det_ID = line[1].strip()
        elif(line[0] == "Detector description"):
            spectra.det_descp = line[1].strip()
        elif(line[0] == "Real Time"):
            spectra.real_time = float(line[1])
        elif(line[0] == "Live Time"):
            spectra.live_time = float(line[1])
        elif(line[0] == "Acquisition start date"):
            spectra.start_date = line[1].strip()
        elif(line[0] == "Acquisition start time"):
            spectra.start_time = line[1].strip() + ":" + line[2] + ":"+line[3]
        elif(line[0] == "Starting channel number"):
            spectra.start_chan_num = int(line[1].strip())
        elif(line[0] == "Number of channels"):
            spectra.num_channels = int(line[1].strip())
        elif(line[0] == "Energy Fit"):
            temp = line[1].strip()
            temp = temp.split(" ")
            spectra.calib_E_fit.append(float(temp[0]))
            spectra.calib_E_fit.append(float(temp[2]))
            spectra.calib_E_fit.append(float(temp[4]))
        elif(line[0] == "FWHM Fit"):
            temp = line[1].strip()
            temp = temp.split(" ")
            spectra.calib_FWHM_fit.append(float(temp[0]))
            spectra.calib_FWHM_fit.append(float(temp[2]))
            spectra.calib_FWHM_fit.append(float(temp[4]))
        elif(line[0] == "SPECTRUM"):
            inspec = True

    print "Finished file read"
    #convert to numpy arrays
    spectra.counts = np.array(spectra.counts)
    spectra.channels = np.array(spectra.channels)
    # calculate additional parameters
    spectra.calc_dead_time()
    spectra.calc_ebins()
    return spectra
