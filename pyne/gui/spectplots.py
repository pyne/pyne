""" """

import matplotlib.pyplot as plt


def plot_spectra(spect):
    """ Create a standard spectrum pulse height plot"""
    plt.plot(spect.Ebin, spect.counts)
    plt.xlim(xmin=10)
    plt.ylim(ymin=1)
    plt.xlabel("Energy (KeV)")
    plt.ylabel("Counts")
    plt.yscale("log")
    plt.show()


def plot_peak(spec1, energy, spread=4):
    """
    Create a filled plot of a region of the spectra around a given energy
    value.

    Inputs: spec1=spectrum object from which to plot
            energy=the energy to centre the plot on
            spread (default=4KeV)= half the energy range to plot i.e the
            number of KeV above and below the energy value to plot

    Returns: nothing
    Displays: matplotlib filled plot

    """
    plt.plot(spec1.Ebin, spec1.counts)
    plt.ylim(ymin=1)
    plt.xlabel("Energy (KeV)")
    plt.ylabel("Counts")
    plt.xlim(xmin=energy-spread, xmax=energy+spread)
    plt.fill(spec1.Ebin, spec1.counts, "g")
    plt.show()