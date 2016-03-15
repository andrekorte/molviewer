# -*- coding: utf-8 -*-
'''

'''
import sys
from jcamp import JCAMP_reader
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore
from pyqtgraph.Point import Point
import pyqtgraph.opengl as gl
import numpy as np
from random import randint
from moldata import spectrumcolors
from spectrum import Spectrum

import pyqtgraph.parametertree.parameterTypes as pTypes
from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType


class SpectrumViewer(QtGui.QMainWindow):
        '''The application main window.'''

        def __init__(self, parent=None):
            # We're not calling the super method so we can span an new instance in
            # the new function.
            # However, when the first created window is closed, all windoes close
            # super(SpectrumViewer, self).__init__()
            QtGui.QMainWindow.__init__(self, parent)
            self.spectra = []
            self.plots = []
            self.initUI()

        def initToolBar(self):
            '''The main toolbar.'''
            # File actions
            self.newAction = QtGui.QAction(QtGui.QIcon('icons/new.png'), '&New', self)
            self.newAction.setStatusTip('New molecule viewer')
            self.newAction.setShortcut('Ctrl+N')
            self.newAction.triggered.connect(self.new)

            self.openAction = QtGui.QAction(QtGui.QIcon('icons/open.png'), '&Open', self)
            self.openAction.setStatusTip('Read molecule from file')
            self.openAction.setShortcut('Ctrl+O')
            self.openAction.triggered.connect(self.open)

            self.saveAction = QtGui.QAction(QtGui.QIcon('icons/save.png'), '&Save', self)
            self.saveAction.setStatusTip('Write molecule to file')
            self.saveAction.setShortcut('Ctrl+S')
            self.saveAction.triggered.connect(self.save)

            self.saveimageAction = QtGui.QAction(QtGui.QIcon('icons/picture-edit.png'), 'Save &image as PNG', self)
            self.saveimageAction.setStatusTip('Save image as PNG')
            self.saveimageAction.setShortcut('Ctrl+I')
            self.saveimageAction.triggered.connect(self.saveimage)

            self.exitAction = QtGui.QAction(QtGui.QIcon('icons/exit.png'), '&Exit', self)
            self.exitAction.setShortcut('Ctrl+Q')
            self.exitAction.setStatusTip('Exit application')
            self.exitAction.triggered.connect(QtGui.qApp.quit)

            # View actions
            self.fullScreenAction = QtGui.QAction(QtGui.QIcon('icons/fullscreen.png'), '&Full Screen', self)
            self.fullScreenAction.setStatusTip('Toggle full screen')
            self.fullScreenAction.setShortcut('F11')
            self.fullScreenAction.triggered.connect(self.fullscreen)

            self.toggleSidebarAction = QtGui.QAction('&Toggle sidebar', self)
            self.toggleSidebarAction.setStatusTip('Toggle sidebar')
            self.toggleSidebarAction.triggered.connect(self.togglesidebar)

            # Spectrum actions
            self.addspectrumAction = QtGui.QAction(QtGui.QIcon('icons/addspectrum.png'), '&Add spectrum', self)
            self.addspectrumAction.setStatusTip('Add a new spectrum')
            self.addspectrumAction.setShortcut('Ctrl + A')
            self.addspectrumAction.triggered.connect(self.open)

            self.removespectrumAction = QtGui.QAction(QtGui.QIcon('icons/removespectrum.png'), '&Remove spectrum', self)
            self.removespectrumAction.setStatusTip('Remove a spectrum')
            self.removespectrumAction.setShortcut('Ctrl + R')
            self.removespectrumAction.triggered.connect(self.removespectrum)


            # Toolbars
            self.filetoolbar = self.addToolBar('File')
            self.filetoolbar.addAction(self.newAction)
            self.filetoolbar.addAction(self.openAction)
            self.filetoolbar.addAction(self.saveAction)
            self.filetoolbar.addAction(self.saveimageAction)

            self.viewtoolbar = self.addToolBar('View')
            self.viewtoolbar.addAction(self.fullScreenAction)

            self.spectrumtoolbar = self.addToolBar('Spectrum')
            self.spectrumtoolbar.addAction(self.addspectrumAction)
            self.spectrumtoolbar.addAction(self.removespectrumAction)

        def initMenuBar(self):
            '''The main menubar.'''
            menubar = self.menuBar()
            menubar.setNativeMenuBar(True)

            # Create menus
            file = menubar.addMenu('&File')
            edit = menubar.addMenu('&Edit')
            view = menubar.addMenu('&View')
            spectrum = menubar.addMenu('&Spectrum')
            about = menubar.addMenu('&About')

            # File menu
            file.addAction(self.newAction)
            file.addAction(self.openAction)
            file.addAction(self.saveAction)
            file.addAction(self.saveimageAction)
            file.addSeparator()
            file.addAction(self.exitAction)

            # View menu
            view.addAction(self.fullScreenAction)
            view.addAction(self.toggleSidebarAction)

            # Spectrum menu
            spectrum.addAction(self.addspectrumAction)
            spectrum.addAction(self.removespectrumAction)

        def initParameterTree(self):
            '''The sidebar.'''
            self.parametertree = ParameterTree()
            self.params = []

        def initUI(self):
            '''Initializes the main window widgets.'''
            self.initToolBar()
            self.initMenuBar()
            self.initParameterTree()

            self.plotwidget = pg.PlotWidget()
            self.plot = self.plotwidget.plotItem
            self.plot.setMenuEnabled(False)
            #cross hair
            self.vLine = pg.InfiniteLine(angle=90, movable=False)
            self.hLine = pg.InfiniteLine(angle=0, movable=False)
            self.plotwidget.addItem(self.vLine, ignoreBounds=True)
            self.plotwidget.addItem(self.hLine, ignoreBounds=True)
            self.vb = self.plot.vb
            self.plot.scene().sigMouseMoved.connect(self.mouseMoved)

            widget = QtGui.QWidget()
            hbox = QtGui.QHBoxLayout(widget)
            hbox.addWidget(self.parametertree)
            hbox.addWidget(self.plotwidget, 1)

            self.setCentralWidget(widget)
            self.resize(1200, 600)
            self.setWindowTitle('Molviewer :: Spectrum main window')
            self.show()

        def mouseMoved(self, evt):
            '''Detect mouse movement.'''
            pos = evt
            if self.plot.sceneBoundingRect().contains(pos):
                mousePoint = self.vb.mapSceneToView(pos)
                self.vLine.setPos(mousePoint.x())
                self.hLine.setPos(mousePoint.y())

        def new(self):
            '''Spawn new main window.'''
            spawn = SpectrumViewer(self)
            spawn.setWindowTitle('Molviewer :: Spectrum main window (child)')
            spawn.show()

        def open(self):
            '''Open a new spectrum using the file dialog.'''
            filenames = QtGui.QFileDialog.getOpenFileNames(self, 'Open Files', '', '*.dx;; *.jdx')
            if filenames:
                for f in filenames:
                    spectrum = Spectrum(f)
                    self.addspectrum(spectrum)

        def save(self):
            pass

        def saveimage(self):
            pass

        def draw(self,spectrum):
            '''Draw a new spectrum.

            :param spectrum: The spectrum to draw.
            :type spectrum: spectrum.Spectrum
            :param color: The color in which the spectrum is drawn.

            '''
            p = self.plotwidget.plot(spectrum.x, spectrum.y, pen=spectrum.color)
            if p not in self.plots:
                self.plots.append(p)
            if spectrum.datatype.lower() == 'infrared spectrum':
                # Plot higher wavenumbers to the left hand side
                self.plotwidget.invertX(True)
                self.plotwidget.setLabel('bottom', '1/cm')
            else:
                self.plotwidget.setLabel('bottom', spectrum.xunits)

        def fullscreen(self):
            '''Toggle full screen window.'''
            if self.windowState() & QtCore.Qt.WindowFullScreen:
                self.showNormal()
            else:
                self.showFullScreen()

        def togglesidebar(self):
            if self.parametertree.isVisible():
                self.parametertree.hide()
            else:
                self.parametertree.show()

        def addspectrum(self, spectrum):
            '''Add a new spectrum.

            :param spectrum: The new spectrum.
            :type spectrum: spectrum.Spectrum

            '''
            self.spectra.append(spectrum)
            ind = self.spectra.index(spectrum)
            spectrum.color = spectrumcolors[ind]
            params = [
                    {'name': spectrum.title, 'type': 'group', 'children': [
                    {'name': 'Number', 'type': 'int', 'value': self.spectra.index(spectrum)},
                    {'name': 'Color', 'type': 'color', 'value': spectrum.color},
                    {'name': 'Visible', 'type': 'bool', 'value': True},
                    {'name': 'Action Parameter', 'type': 'action'},
                    {'name': 'Spectrum', 'type': 'str', 'value': spectrum, 'visible': False}]
                    }]

            p = Parameter.create(name='Files', type='group', children=params)
            p.sigTreeStateChanged.connect(self.change)
            self.parametertree.addParameters(p, showTop=False)
            self.params.append(p)
            self.draw(spectrum)

        def removespectrum(self):
            '''Remove a spectrum from display.'''
            num, ok = QtGui.QInputDialog.getText(self, 'Remove a spectrum', 'Enter spectrum number:')
            if ok and num:
                num = int(num)
                self.plotwidget.removeItem(self.plots[num])
                self.plots.pop(num)
                self.spectra.pop(num)
                self.params.pop(num)
                self.parametertree.clear()
                for p in self.params:
                    self.parametertree.addParameters(p, showTop=False)

        def change(self, param, changes):
            for param, change, data in changes:
                spec = param.parent().child('Spectrum').value()
                index = self.spectra.index(spec)

                if isinstance(param.value(), QtGui.QColor):
                    self.plots[index].setPen(param.value())
                    self.spectra[index].color = param.value()

                if param.name() == 'Visible':
                    if data == False:
                        self.plots[index].clear()
                    elif data == True:
                        self.plots[index].setData(self.spectra[index].x, self.spectra[index].y)


def main():
    '''The main function

    Args: None
    Returns: None
    '''
    app = QtGui.QApplication(sys.argv)
    app.setWindowIcon(QtGui.QIcon('icons/flask-256.png'))
    spec = SpectrumViewer()
    s = Spectrum('spectrum.dx')
    spec.addspectrum(s)
    #spec.draw()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
