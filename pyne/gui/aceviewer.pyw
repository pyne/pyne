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

        self.populateTables()
        
    def _CreateGui(self):
        # Set title of window
        self.setWindowTitle("ACE Data Viewer")

        # Create widgets
        self.main = QWidget()
        self.setCentralWidget(self.main)

        tableLabel = QLabel("ACE Table:")
        self.tableCombo = QComboBox()
        tableLayout = QHBoxLayout()
        tableLayout.addWidget(tableLabel)
        tableLayout.addWidget(self.tableCombo)

        # Create reaction list view
        self.reactionList = QListWidget()
        self.reactionList.setSelectionMode(QAbstractItemView.ExtendedSelection)

        # Create canvas
        self.fig = Figure(figsize=(400,400), dpi=72, facecolor=(1,1,1), edgecolor=(0,0,0))
        self.canvas = FigureCanvas(self.fig)

        bottomLayout = QHBoxLayout()
        bottomLayout.addWidget(self.reactionList)
        bottomLayout.addWidget(self.canvas)

        line = QFrame(self.main)
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)

        # Set layout
        layout = QVBoxLayout()
        layout.addLayout(tableLayout)
        layout.addWidget(line)
        layout.addLayout(bottomLayout)
        self.main.setLayout(layout)

        # Create menu bar
        self.menubar = QMenuBar(self)
        self.menuFile = QMenu("&File", self.menubar)
        self.menuHelp = QMenu("&Help", self.menubar)
        self.menubar.addActions([self.menuFile.menuAction()])
        self.setMenuBar(self.menubar)

        # File menu
        self.actionOpen = QAction("&Open ACE Library...", self)
        self.actionExit = QAction("E&xit", self)
        self.menuFile.addActions([self.actionOpen, self.actionExit])

        # Actions
        self.connect(self.actionOpen, SIGNAL("triggered()"), self.openLibrary)
        self.connect(self.actionExit, SIGNAL("triggered()"), self.close)
        self.connect(self.tableCombo, SIGNAL("currentIndexChanged(int)"),
                     self.populateReactions)
        self.connect(self.reactionList, SIGNAL("itemSelectionChanged()"),
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
        self.populateTables()

    def populateTables(self):
        self.tableCombo.clear()
        for table in self.tables:
            self.tableCombo.addItem(table.name)
        if self.tables:
            self.tableCombo.setEnabled(True)
        else:
            self.tableCombo.setDisabled(True)

    def populateReactions(self):
        self.reactionList.clear()

        table = self._currentTable()
        if table:
            for reaction in table:
                self.reactionList.addItem(ace.reaction_names[reaction.MT])

        self.drawPlot()
            

    def drawPlot(self):
        self.fig.clear()

        items = self.reactionList.selectedItems()

        if items: 
            self.axes = self.fig.add_subplot(111)
            table = self._currentTable()

            for item in self.reactionList.selectedItems():
                index = self.reactionList.row(item)
                reaction = table.reactions[index]
                
                if reaction.MT == 1:
                    self.axes.loglog(table.energy, reaction.sigma)
                elif reaction.MT == 27:
                    self.axes.loglog(table.energy, reaction.sigma)
                else:
                    self.axes.loglog(table.energy[reaction.IE-1:], reaction.sigma)
                self.axes.grid(True)

        self.canvas.draw()

    def _currentTable(self):
        indexTable = self.tableCombo.currentIndex()
        return self.tables[indexTable]

    def _currentReaction(self):
        indexReaction = self.reactionList.currentRow()
        print(indexReaction)
        if indexReaction >= 0:
            return self._currentTable().reactions[indexReaction]
        else:
            return None


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = AceViewer()
    window.show()
    app.exec_()
    sys.exit()
