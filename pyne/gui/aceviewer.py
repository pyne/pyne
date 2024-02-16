#!/usr/bin/env python

"""This module is for plotting ACE-format cross sections based on the ACE
module.

.. moduleauthor:: Paul Romano <paul.k.romano@gmail.com>
"""
import sys
from bisect import bisect_right
from pyne.utils import QA_warn

import numpy as np
import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import (
    NavigationToolbar2QTAgg as NavigationToolbar,
)
from matplotlib.figure import Figure

# For the time being we're using PyQt4 since matplotlib-support for PySide is
# very recent
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from pyne import ace

QA_warn(__name__)

matplotlib.use("Qt4Agg")


class AceViewer(QMainWindow):
    def __init__(self, parent=None):
        super(AceViewer, self).__init__(parent)

        # Create GUI elements
        self._create_gui()

        # Initial data structures
        self.tables = []

        self.populate_reactions()

    def _create_gui(self):
        # Set title of window
        self.setWindowTitle("ACE Data Viewer")

        # Create widgets
        self.main = QWidget()
        self.setCentralWidget(self.main)

        # Create reaction list view
        self.reactionTree = MyTreeWidget()
        self.reactionTree.setColumnCount(1)
        self.reactionTree.setHeaderLabels(["Nuclides/Reactions"])
        self.reactionTree.setMinimumWidth(200)
        self.reactionTree.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.reactionTree.setContextMenuPolicy(Qt.DefaultContextMenu)

        # Create canvas
        self.fig = Figure(dpi=72, facecolor=(1, 1, 1), edgecolor=(0, 0, 0))
        self.canvas = FigureCanvas(self.fig)

        # Create toolbar
        self.mplToolbar = NavigationToolbar(self.canvas, self.main)

        # Create layout for canvas and toolbar
        drawLayout = QVBoxLayout()
        drawLayout.addWidget(self.canvas)
        drawLayout.addWidget(self.mplToolbar)

        layout = QHBoxLayout()
        layout.addWidget(self.reactionTree)
        layout.addLayout(drawLayout)

        # Set layout
        self.main.setLayout(layout)

        # Create menu bar
        self.menubar = QMenuBar(self)
        self.menuFile = QMenu("&File", self.menubar)
        self.menuHelp = QMenu("&Help", self.menubar)
        self.menubar.addActions([self.menuFile.menuAction()])
        self.setMenuBar(self.menubar)

        # File menu
        self.actionOpen = QAction("&Open Library...", self)
        self.actionOpenPartial = QAction("Open &Partial Library...", self)
        self.actionExit = QAction("E&xit", self)
        self.menuFile.addActions(
            [self.actionOpen, self.actionOpenPartial, self.actionExit]
        )

        # Actions
        self.connect(self.actionOpen, SIGNAL("triggered()"), self.open_library)
        self.connect(
            self.actionOpenPartial, SIGNAL("triggered()"), self.open_partial_library
        )
        self.connect(self.actionExit, SIGNAL("triggered()"), self.close)
        self.connect(
            self.reactionTree, SIGNAL("itemSelectionChanged()"), self.draw_plot
        )

    def open_library(self):
        """Select and open an ACE file and store data in memory."""

        filename = QFileDialog.getOpenFileName(
            self, "Load ACE Library", "./", "ACE Libraries (*)"
        )

        try:
            if filename:
                # Parse ACE library
                lib = ace.Library(filename)
                lib.read()

                # Append tables into self.tables object
                for table in lib.tables.values():
                    self.tables.append(table)
        except:
            pass

        # Sort tables based on name
        self.tables.sort(key=lambda table: table.name)

        # Reset combo box
        self.populate_reactions()

    def open_partial_library(self):
        """Select and open an ACE file and store data in memory."""

        filename = QFileDialog.getOpenFileName(
            self, "Load ACE Library", "./", "ACE Libraries (*)"
        )

        table_names, completed = QInputDialog.getText(
            self, "Nuclides", "Enter nuclides:"
        )
        if completed:
            table_names = str(table_names).split()
        else:
            return

        try:
            if filename:
                # Parse ACE library
                lib = ace.Library(filename)
                lib.read(table_names)

                # Append tables into self.tables object
                for table in lib.tables.values():
                    self.tables.append(table)
        except:
            pass

        # Sort tables based on name
        self.tables.sort(key=lambda table: table.name)

        # Reset combo box
        self.populate_reactions()

    def populate_reactions(self):
        self.reactionTree.clear()

        for table in self.tables:
            # Add top-level item
            tableItem = QTreeWidgetItem(self.reactionTree, [table.name])
            tableItem.setData(0, Qt.UserRole, table)

            # Add item for total reaction
            item = QTreeWidgetItem(tableItem, ["(n,total)"])

            for reaction in table:
                # Add sub-item
                try:
                    reactionName = ace.reaction_names[reaction.MT]
                except:
                    reactionName = "MT = {0}".format(reaction.MT)
                item = QTreeWidgetItem(tableItem, [reactionName])
                item.setData(0, Qt.UserRole, reaction)

        self.draw_plot()

    def draw_plot(self):
        # Clears the current figure
        self.fig.clear()

        # Get all selected reactions
        items = self.reactionTree.selectedItems()

        if len(items) > 0:
            # Create instance of Axes on the Figure
            self.axes = self.fig.add_subplot(111)

            for item in items:
                # Get Reaction object stored in QTreeWidgetItem
                reaction = item.data(0, Qt.UserRole).toPyObject()

                # Handle total reaction separately
                if item.text(0) == "(n,total)":
                    # Get NeutronTable object
                    table = item.parent().data(0, Qt.UserRole).toPyObject()

                    # Plot total cross section
                    self.axes.loglog(table.energy, table.sigma_t)
                    continue

                # Make sure that the data stored in QTreeWidgetItem is actually
                # a Reaction instance
                if not isinstance(reaction, ace.Reaction):
                    continue

                # Get reference to NeutronTable containing Reaction
                table = reaction.table

                # Plot reaction cross section
                self.axes.loglog(table.energy[reaction.IE :], reaction.sigma)

            # Customize plot
            self.axes.grid(True)
            self.axes.set_xlabel("Energy (MeV)")
            self.axes.set_ylabel("Cross section (barns)")

            # Display plot on FigureCanvas
            self.canvas.draw()


class MyTreeWidget(QTreeWidget):
    def __init__(self):
        super(MyTreeWidget, self).__init__()

    def mousePressEvent(self, event):
        if event.button() == Qt.RightButton:
            self.temporaryPos = self.viewport().mapFromGlobal(event.globalPos())
            menu = QMenu()
            menu.addAction("Plot Angle Distribution", self.create_angle)
            menu.addAction("Plot Angle Distribution (Polar)", self.create_angle_polar)
            menu.addAction("Plot Energy Distribution", self.create_energy)
            menu.exec_(event.globalPos())
        else:
            super(MyTreeWidget, self).mousePressEvent(event)

    def create_angle(self, polar=False):
        # Get QTreeWidgetItem that was right-clicked
        item = self.itemAt(self.temporaryPos)

        # Make sure there's an actual item
        if not item:
            return

        # Get reaction data from the QTreeWidgetItem
        reaction = item.data(0, Qt.UserRole).toPyObject()

        if isinstance(reaction, ace.Reaction):
            # Ask user for incoming energies
            energies, completed = QInputDialog.getText(
                self.parent(),
                "Energies",
                "Enter incoming energies at "
                "which to plot angular distribution (MeV):",
            )

            # Plot angular distribution
            if completed:
                energies = map(float, str(energies).split())
                anglePlot = DistributionPlot(
                    reaction, energies, polar=polar, parent=self.parent()
                )
                anglePlot.show()

    def create_angle_polar(self):
        self.create_angle(polar=True)

    def create_energy(self):
        # Get QTreeWidgetItem that was right-clicked
        item = self.itemAt(self.temporaryPos)

        # Make sure there's an actual item
        if not item:
            return

        # Get reaction data from the QTreeWidgetItem
        reaction = item.data(0, Qt.UserRole).toPyObject()

        if isinstance(reaction, ace.Reaction):
            # Ask user for incoming energies
            energies, completed = QInputDialog.getText(
                self.parent(),
                "Energies",
                "Enter incoming energies at "
                "which to plot energy distribution (MeV):",
            )

            # Plot angular distribution
            if completed:
                energies = map(float, str(energies).split())
                anglePlot = DistributionPlot(
                    reaction, energies, energydist=True, parent=self.parent()
                )
                anglePlot.show()


class DistributionPlot(QMainWindow):
    def __init__(self, reaction, energies, polar=False, energydist=False, parent=None):
        super(DistributionPlot, self).__init__(parent)

        # Initialize data
        self.reaction = reaction
        self.energies = energies
        self.polar = polar
        self.energydist = energydist

        # Draw widgets
        self._create_gui()

        # Draw plot
        self._draw_plot()

    def _create_gui(self):
        # Set title of window
        if self.energydist:
            self.setWindowTitle("Energy Distribution")
        else:
            self.setWindowTitle("Angle Distribution")

        # Create widgets
        self.main = QWidget()
        self.setCentralWidget(self.main)

        # Create canvas
        self.fig = Figure(dpi=72, facecolor=(1, 1, 1), edgecolor=(0, 0, 0))
        self.canvas = FigureCanvas(self.fig)

        # Create toolbar
        self.mplToolbar = NavigationToolbar(self.canvas, self.main)

        # Create layout for canvas and toolbar
        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        layout.addWidget(self.mplToolbar)

        # Set layout
        self.main.setLayout(layout)

    def _draw_plot(self):
        # Create instance of Axes on the Figure
        if self.polar:
            self.axes = self.fig.add_subplot(111, polar=True)
        else:
            self.axes = self.fig.add_subplot(111)

        # Loop over each incoming energy
        for E in self.energies:
            if self.energydist:
                self._plot_energy_dist(E)
                self.axes.grid(True)
            else:
                if self.polar:
                    self._plot_angle_dist_polar(E)
                else:
                    self._plot_angle_dist(E)
                    self.axes.grid(True)
        self.axes.legend()

    def _plot_angle_dist(self, E_in):

        # determine index for incoming energy
        index = bisect_right(self.reaction.ang_energy_in, E_in)

        # plot distribution
        self.axes.plot(
            self.reaction.ang_cos[index],
            self.reaction.ang_pdf[index],
            label="E = {0} MeV".format(E_in),
        )

    def _plot_angle_dist_polar(self, E_in):
        """
        Plots the secondary angle distribution for this reaction at a
        given incoming energy of the particle.
        """

        # determine index for incoming energy
        index = bisect_right(self.reaction.ang_energy_in, E_in)

        # Find angles and probabilities (cos from 0 to pi)
        angles = np.arccos(self.reaction.ang_cos[index])[::-1]
        pdf = self.reaction.ang_pdf[index][::-1]

        theta = np.linspace(0, np.pi, 100)
        r = np.interp(theta, angles, pdf)

        theta = np.concatenate((theta, theta + np.pi))
        r = np.concatenate((r, r[::-1]))

        # plot angle distribution
        self.axes.plot(theta, r, label="E = {0} MeV".format(E_in))

    def _plot_energy_dist(self, E_in):
        """
        Plots the secondary energy distribution for this reaction at a
        given incoming energy if data are available.
        """

        try:
            # determine index for incoming energy
            index = bisect_right(self.reaction.e_dist_energy_in, E_in)

            # plot energy distribution
            self.axes.semilogx(
                self.reaction.e_dist_energy_out[index],
                self.reaction.e_dist_pdf[index],
                label="E = {0} MeV".format(E_in),
            )
        except:
            QMessageBox.warning(
                self,
                "Unsupported Energy Distribution",
                "Energy distribution for law {0} not yet supported!".format(
                    self.reaction.e_dist_law
                ),
            )


def main():
    app = QApplication(sys.argv)
    window = AceViewer()
    window.show()
    app.exec_()
    sys.exit()


if __name__ == "__main__":
    main()
