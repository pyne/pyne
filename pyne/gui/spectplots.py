""" """

import matplotlib.pyplot as plt


def plot_spectrum(spect):
    """ Create a standard spectrum pulse height plot"""
    plt.plot(spect.ebin, spect.counts)
    plt.xlim(xmin=10)
    plt.ylim(ymin=1)
    plt.xlabel("Energy (KeV)")
    plt.ylabel("Counts")
    plt.yscale("log")
    plt.title(spect.file_name)
    plt.show()


def plot_peak(spec1, energy, spread=4):
    """
    Create a filled plot of a region of the spectra around a given energy
    value.
    """
    plt.plot(spec1.ebin, spec1.counts)
    plt.ylim(ymin=1)
    plt.xlabel("Energy (KeV)")
    plt.ylabel("Counts")
    plt.xlim(xmin=energy-spread, xmax=energy+spread)
    plt.fill(spec1.ebin, spec1.counts, "g")
    plt.title(spec1.file_name + " " + str(energy))
    plt.show()