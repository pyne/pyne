#!/usr/bin/env python

"""This module is for plotting ACE-format cross sections based on the ACE
module.

.. moduleauthor:: Paul Romano <paul.k.romano@gmail.com>
"""

import sys

import matplotlib
matplotlib.use('Qt4Agg')
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# For the time being we're using 
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from pyne import ace

class AceViewer(QMainWindow):

    def __init__(self, parent=None):
        super(AceViewer, self).__init__(parent)

        # Create GUI elements
        self._CreateGui()
        
        # Initial data structures
        self.tables = []

        self.populateReactions()
        
    def _CreateGui(self):
        # Set title of window
        self.setWindowTitle("ACE Data Viewer")

        # Create widgets
        self.main = QWidget()
        self.setCentralWidget(self.main)

        # Create reaction list view
        self.reactionTree = QTreeWidget()
        self.reactionTree.setColumnCount(1)
        self.reactionTree.setHeaderLabels(["Nuclides/Reactions"])
        self.reactionTree.setMinimumWidth(200)
        self.reactionTree.setSelectionMode(QAbstractItemView.ExtendedSelection)

        # Create canvas
        self.fig = Figure(figsize=(400,200), dpi=72, facecolor=(1,1,1), edgecolor=(0,0,0))
        self.canvas = FigureCanvas(self.fig)

        layout = QHBoxLayout()
        layout.addWidget(self.reactionTree)
        layout.addWidget(self.canvas)

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
        self.actionOpenPartial = QAction("&Open Partial Library...", self)
        self.actionExit = QAction("E&xit", self)
        self.menuFile.addActions([self.actionOpen, self.actionOpenPartial,
                                  self.actionExit])

        # Actions
        self.connect(self.actionOpen, SIGNAL("triggered()"), self.openLibrary)
        self.connect(self.actionOpenPartial, SIGNAL("triggered()"),
                     self.openPartialLibrary)
        self.connect(self.actionExit, SIGNAL("triggered()"), self.close)
        self.connect(self.reactionTree, SIGNAL("itemSelectionChanged()"),
                     self.drawPlot)

    def openLibrary(self): 
        """Select and open an ACE file and store data in memory."""

        filename = QFileDialog.getOpenFileName(self, "Load ACE Library", "./",
                                               "ACE Libraries (*)")

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
        self.populateReactions()

    def openPartialLibrary(self):
        """Select and open an ACE file and store data in memory."""

        filename = QFileDialog.getOpenFileName(self, "Load ACE Library", "./",
                                               "ACE Libraries (*)")

        table_names, completed = QInputDialog.getText(
            self, "Nuclides", "Enter nuclides:")
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
        self.populateReactions()

    def populateReactions(self):
        self.reactionTree.clear()

        for table in self.tables:
            # Add top-level item
            tableItem = QTreeWidgetItem(self.reactionTree, [table.name])
            tableItem.setData(0, Qt.UserRole, table)

            for reaction in table:
                # Add sub-item
                try:
                    reactionName = ace.reaction_names[reaction.MT]
                except:
                    reactionName = "MT = {0}".format(reaction.MT)
                item = QTreeWidgetItem(tableItem, [reactionName])
                item.setData(0, Qt.UserRole, reaction)

        self.drawPlot()

    def drawPlot(self):
        self.fig.clear()

        items = self.reactionTree.selectedItems()

        if len(items) > 0:
            self.axes = self.fig.add_subplot(111)

            for item in items:
                reaction = item.data(0, Qt.UserRole).toPyObject()

                if not isinstance(reaction, ace.Reaction):
                    continue

                # Get table
                table = reaction.table
                
                if reaction.MT == 1:
                    self.axes.loglog(table.energy, reaction.sigma)
                elif reaction.MT == 27:
                    self.axes.loglog(table.energy, reaction.sigma)
                else:
                    self.axes.loglog(table.energy[reaction.IE-1:], reaction.sigma)
                
            self.axes.grid(True)
            self.axes.set_xlabel('Energy (MeV)')
            self.axes.set_ylabel('Cross section (barns)')
            self.canvas.draw()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = AceViewer()
    window.show()
    app.exec_()
    sys.exit()
