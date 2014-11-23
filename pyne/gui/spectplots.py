"""Plotting routines for spectrometry modules"""
from warnings import warn
from pyne.utils import QAWarning

import matplotlib.pyplot as plt

warn(__name__ + " is not yet QA compliant.", QAWarning)


def plot_spectrum(spect):
    """Create a standard spectrum pulse height plot"""
    plt.plot(spect.ebin, spect.counts)
    plt.xlim(xmin=10)
    plt.ylim(ymin=1)
    plt.xlabel("Energy (keV)")
    plt.ylabel("Counts")
    plt.yscale("log")
    plt.title(spect.file_name)
    plt.show()


def plot_peak(spect, energy, spread=4):
    """Create a filled plot of a region of the spectra around a given energy
    value.
    """
    plt.plot(spect.ebin, spect.counts)
    plt.ylim(ymin=1)
    plt.xlabel("Energy (keV)")
    plt.ylabel("Counts")
    plt.xlim(xmin=energy-spread, xmax=energy+spread)
    plt.fill(spect.ebin, spect.counts, "g")
    plt.title(spect.file_name + " " + str(energy) + "keV")
    plt.show()
