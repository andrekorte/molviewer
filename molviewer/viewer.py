# -*- coding: utf-8 -*-
'''
A simple molecular viewer.

'''
import sys
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph.opengl as gl
from pyqtgraph.parametertree import ParameterTree, Parameter
import numpy as np
from molio import get_mol_by_name, read_xyz, write_xyz
from molfuncs import angle_between_vectors, distance_matrix, moleculeinfo


class MolViewer(QtGui.QMainWindow):
    '''The application main window.'''

    def __init__(self, parent=None):
        # We're not calling the super method so we can span an new instance in
        # the new function.
        # However, when the first created window is closed, all windoes close
        # super(MolViewer, self).__init__()
        QtGui.QMainWindow.__init__(self, parent)
        self.mol = None
        self.initUI()

    def initToolBar(self):
        '''The main toolbar.'''
        # File actions
        self.newAction = QtGui.QAction(QtGui.QIcon('../icons/new.png'), '&New', self)
        self.newAction.setStatusTip('New molecule viewer')
        self.newAction.setShortcut('Ctrl+N')
        self.newAction.triggered.connect(self.new)

        self.openAction = QtGui.QAction(QtGui.QIcon('../icons/open.png'), '&Open', self)
        self.openAction.setStatusTip('Read molecule from file')
        self.openAction.setShortcut('Ctrl+O')
        self.openAction.triggered.connect(self.open)

        self.saveAction = QtGui.QAction(QtGui.QIcon('../icons/save.png'), '&Save', self)
        self.saveAction.setStatusTip('Write molecule to file')
        self.saveAction.setShortcut('Ctrl+S')
        self.saveAction.triggered.connect(self.save)

        self.saveasAction = QtGui.QAction(QtGui.QIcon.fromTheme('document-save-as'), '&Save as...', self)
        self.saveasAction.setStatusTip('Write molecule to file')
        self.saveasAction.setShortcut('Ctrl+Shift+S')
        self.saveasAction.triggered.connect(self.saveas)

        self.saveimageAction = QtGui.QAction(QtGui.QIcon('../icons/picture-edit.png'), 'Save &image as PNG', self)
        self.saveimageAction.setStatusTip('Save image as PNG')
        self.saveimageAction.setShortcut('Ctrl+I')
        self.saveimageAction.triggered.connect(self.saveimage)

        self.exitAction = QtGui.QAction(QtGui.QIcon('../icons/exit.png'), '&Exit', self)
        self.exitAction.setShortcut('Ctrl+Q')
        self.exitAction.setStatusTip('Exit application')
        self.exitAction.triggered.connect(QtGui.qApp.quit)

        # Edit actions
        self.clearAction = QtGui.QAction(QtGui.QIcon('../icons/clear.png'), 'Clea&r', self)
        self.clearAction.setShortcut('Ctrl+R')
        self.clearAction.setStatusTip('Remove molecule from viewer')
        self.clearAction.triggered.connect(self.clear)

        self.preferencesAction = QtGui.QAction(QtGui.QIcon('../icons/preferences.png'), '&Options', self)
        self.preferencesAction.setStatusTip('Edit Preferences')
        # self.preferencesAction.setShortcut()
        self.preferencesAction.triggered.connect(self.editpreferences)

        # View actions
        self.fullScreenAction = QtGui.QAction(QtGui.QIcon('../icons/fullscreen.png'), '&Full Screen', self)
        self.fullScreenAction.setStatusTip('Toggle full screen')
        self.fullScreenAction.setShortcut('F11')
        self.fullScreenAction.triggered.connect(self.fullscreen)

        self.setwindowsizeAction = QtGui.QAction(QtGui.QIcon('../icons/fullscreen.png'), '&Set Geometry', self)
        self.setwindowsizeAction.setStatusTip('Set height and width of main window')
        self.setwindowsizeAction.setShortcut('F12')
        self.setwindowsizeAction.triggered.connect(self.setwindowsize)

        # Molecule action
        self.getbynameAction = QtGui.QAction(QtGui.QIcon('../icons/find.png'), '&Get molecule by name or CAS number', self)
        self.getbynameAction.setStatusTip('Get a molecule by name or CAS number')
        self.getbynameAction.setShortcut('Ctrl+G')
        self.getbynameAction.triggered.connect(self.getbyname)

        self.distmatAction = QtGui.QAction(QtGui.QIcon('../icons/table.png'), '&Distance matrix', self)
        self.distmatAction.setShortcut('Ctrl+D')
        self.distmatAction.triggered.connect(self.distmatWindow)

        self.infoAction = QtGui.QAction(QtGui.QIcon('../icons/info.png'), '&Information', self)
        self.infoAction.setStatusTip('Show basic information about molecule')
        self.infoAction.setShortcut('Ctrl+I')
        self.infoAction.triggered.connect(self.info)

        # Spectrum action
        self.showspectrumAction = QtGui.QAction(QtGui.QIcon('../icons/spectrum.png'), '&Spectrum', self)
        self.showspectrumAction.setStatusTip('Show spectrum')
        self.showspectrumAction.triggered.connect(self.showspectrum)

        # About action
        self.aboutAction = QtGui.QAction(QtGui.QIcon('../icons/spectrum.png'), '&About', self)
        self.aboutAction.setStatusTip('Show the about dialog')
        self.aboutAction.triggered.connect(self.about)

        # Toolbars
        self.filetoolbar = self.addToolBar('File')
        self.filetoolbar.addAction(self.newAction)
        self.filetoolbar.addAction(self.openAction)
        self.filetoolbar.addAction(self.saveAction)
        # self.filetoolbar.addAction(self.saveasAction)
        self.filetoolbar.addAction(self.saveimageAction)

        self.edittoolbar = self.addToolBar('Edit')
        self.edittoolbar.addAction(self.preferencesAction)
        self.edittoolbar.addAction(self.clearAction)

        self.viewtoolbar = self.addToolBar('View')
        self.viewtoolbar.addAction(self.fullScreenAction)

        self.moleculetoolbar = self.addToolBar('Molecule')
        self.moleculetoolbar.addAction(self.getbynameAction)
        self.moleculetoolbar.addAction(self.distmatAction)
        self.moleculetoolbar.addAction(self.infoAction)

        self.spectrumtoolbar = self.addToolBar('Spectrum')
        self.spectrumtoolbar.addAction(self.showspectrumAction)

    def initMenuBar(self):
        '''The main menubar.'''
        menubar = self.menuBar()
        menubar.setNativeMenuBar(True)

        # Create menus
        file = menubar.addMenu('&File')
        edit = menubar.addMenu('&Edit')
        view = menubar.addMenu('&View')
        molmenu = menubar.addMenu('&Molecule')
        spectrummenu = menubar.addMenu('&Spectrum')
        aboutmenu = menubar.addMenu('&About')

        # File menu
        file.addAction(self.newAction)
        file.addAction(self.openAction)
        file.addAction(self.saveAction)
        file.addAction(self.saveasAction)
        file.addAction(self.saveimageAction)
        file.addSeparator()
        file.addAction(self.exitAction)

        # Edit menu
        edit.addAction(self.clearAction)
        edit.addSeparator()
        edit.addAction(self.preferencesAction)

        # View menu
        view.addAction(self.fullScreenAction)
        view.addAction(self.setwindowsizeAction)

        # Molecule menu
        molmenu.addAction(self.getbynameAction)
        molmenu.addSeparator()
        molmenu.addAction(self.distmatAction)
        molmenu.addAction(self.infoAction)

        # Spectrum menu
        spectrummenu.addAction(self.showspectrumAction)

        # About menu
        aboutmenu.addAction(self.aboutAction)

    def initUI(self):
        '''Initializes the main window widgets.'''
        self.initToolBar()
        self.initMenuBar()

        self.parametertree = ParameterTree()
        self.viewwidget = gl.GLViewWidget()

        widget = QtGui.QWidget()
        hbox = QtGui.QHBoxLayout(widget)
        hbox.addWidget(self.parametertree)
        hbox.addWidget(self.viewwidget, 1)

        self.setCentralWidget(widget)
        self.resize(1200, 600)
        self.setWindowTitle('MolViewer :: Main window')
        self.drawgrid()
        self.show()

    def distmatWindow(self):
        '''Shows the distance matrix window.'''
        if self.mol is None:
            QtGui.QMessageBox.information(self, 'Information', 'No molecule loaded')
        self.matrixtable = TableWindow(parent=self, mol=self.mol, kind='distance')
        self.matrixtable.show()

    def _open(self, mol):
        '''Draws a new molecule

        :param molecule: The molecule for which to draw the atoms.
        :type molecule: molecule.Molecule

        '''
        self.mol = mol
        self.clearview()
        self.drawmolecule(mol)
        self.drawgrid()
        self.setWindowTitle('MolViewer :: Main window :: {:s}'.format(self.mol.filename))

    def open(self):
        '''Open a new molecule using the file dialog.'''
        filename = QtGui.QFileDialog.getOpenFileName(self, 'Open File', './xyz', '*.*;; *.xyz;; *.mol;; *.sdf')
        if filename:
            mol = read_xyz(filename)
            self._open(mol)

    def new(self):
        '''Spawns a new main window.'''
        spawn = MolViewer(self)
        spawn.show()

    def saveas(self):
        '''Save a molecule to file with a new filename.'''
        if self.mol is None:
            QtGui.QMessageBox.information(self, 'Information', 'No molecule to saves')
        else:
            filename = QtGui.QFileDialog.getSaveFileName(self, 'Save file as', './molecules/xyz', '*.xyz')
            if filename:
                write_xyz(self.mol, filename)

    def save(self):
        '''Save the molecule to file.'''
        if self.mol and not self.mol.filename:
            self.saveas()
        elif self.mol and self.mol.filename:
            write_xyz(self.mol, self.mol.filename)

    def saveimage(self):
        '''Saves the viewwidget as a png image.'''
        filename = QtGui.QFileDialog.getSaveFileName(self, 'Save image as PNG', '', '*.png')
        img = self.viewwidget.readQImage()
        # print('Writing image {:s}'.format(filename), end=' ')
        img.save(filename, 'png', 100)
        print('Done')

    def fullscreen(self):
        '''Toggle full screen window.'''
        if self.windowState() & QtCore.Qt.WindowFullScreen:
            self.showNormal()
        else:
            self.showFullScreen()

    def setwindowsize(self):
        '''Set the window to specified size.'''
        width, ok = QtGui.QInputDialog.getText(self, 'Set width', 'Enter window width in pixels:')
        height, ok = QtGui.QInputDialog.getText(self, 'Set height', 'Enter window height in pixels:')
        self.resize(int(width), int(height))

    def editpreferences(self):
        '''Opens the preferences dialog.'''
        dlg = OptionsDialog(self)
        dlg.exec_()

    def getbyname(self):
        '''Displays the input dialog to find a molecule by its identifier.'''
        text, ok = QtGui.QInputDialog.getText(
            self, 'Get molecule by name or CAS number', 'Enter molecule name or CAS number:')
        if ok:
            if isinstance(text, str) and len(text) > 0:
                mol = get_mol_by_name(text)
                if mol:
                    self._open(mol)
                else:
                    QtGui.QMessageBox.information(self, 'Information', 'Molecule not found')

    def drawmolecule(self, mol):
        '''Creates the polygons to draw the molecule to the viewwidget.

        :param molecule: The molecule for which to draw the atoms.
        :type molecule: molecule.Molecule

        '''
        self.viewwidget.setCameraPosition(distance=10)
        self.drawatoms(mol)
        self.drawbonds(mol)
        self.setparameter(mol)

    def drawatoms(self, mol):
        '''Draw the atoms.

        :param molecule: The molecule for which to draw the atoms.
        :type molecule: molecule.Molecule

        '''
        scale = 0.25
        rows = 16
        cols = 16
        shader = 'shaded'
        smooth = True

        for atom in mol.atoms:
            radius = atom.radius() * scale
            pos = atom.r
            md = gl.MeshData.sphere(rows=rows, cols=cols, radius=radius)
            m = gl.GLMeshItem(meshdata=md, smooth=smooth, shader=shader)
            colors = np.ones((md.faceCount(), 4), dtype=float)
            rgb = atom.color()
            colors[:, 0] = rgb[0]  # r
            colors[:, 1] = rgb[1]  # g
            colors[:, 2] = rgb[2]  # b
            colors[:, 3] = 1.0     # a
            md.setFaceColors(colors)
            m.translate(pos[0], pos[1], pos[2])
            self.viewwidget.addItem(m)

    def drawbonds(self, mol):
        '''Draw the bonds.

        :param molecule: The molecule for which to show the bonds.
        :type molecule: molecule.Molecule

        '''
        rows = 16
        cols = 16
        shader = 'shaded'
        smooth = True
        bond_radius = [0.1, 0.1]
        draw_edges = False
        edge_color = (0, 0, 0, 1)
        z = np.asarray([0, 0, 1])
        for bond in mol.bonds():
            vec = bond.p1 - bond.p2
            angle = angle_between_vectors(vec, z) * 180 / np.pi
            axis = np.cross(z, vec)
            md = gl.MeshData.cylinder(rows=rows, cols=cols, radius=bond_radius, length=bond.length())
            m = gl.GLMeshItem(meshdata=md, smooth=smooth, drawEdges=draw_edges, edgeColor=edge_color, shader=shader)
            colors = np.ones((md.faceCount(), 4), dtype=float)
            rgb = bond.atoms[0].color()
            colors[:, 0] = rgb[0]
            colors[:, 1] = rgb[1]
            colors[:, 2] = rgb[2]
            colors[:, 3] = 1.0
            md.setFaceColors(colors)
            m.rotate(angle, axis[0], axis[1], axis[2])
            m.translate(bond.p2[0], bond.p2[1], bond.p2[2])
            self.viewwidget.addItem(m)

    def drawgrid(self):
        '''Show a grid in x,y,z planes.'''
        # gx = gl.GLGridItem()
        # gx.rotate(90, 0, 1, 0)
        # w.addItem(gx)
        # gy = gl.GLGridItem()
        # gy.rotate(90, 1, 0, 0)
        # w.addItem(gy)
        gz = gl.GLGridItem()
        self.viewwidget.addItem(gz)

    def clearview(self):
        '''Removes all items from ViewWidget.'''
        # Not all items are removed in one loop. Why?
        while len(self.viewwidget.items) > 0:
            for i in self.viewwidget.items:
                self.viewwidget.removeItem(i)

    def clear(self):
        '''Clears the viewwidget and deletes the molecule.'''
        self.clearview()
        self.mol = None

    def info(self):
        '''Shows a info box with basic infoemation about the molecule.'''
        if self.mol is None:
            QtGui.QMessageBox.information(self, 'Information', 'No molecule loaded')
        self.infotree = TreeWindow(parent=self, mol=self.mol)
        self.infotree.show()

    def showspectrum(self):
        win = SpectrumWindow(self)
        win.show()

    def about(self):
        dlg = AboutDialog(self)
        dlg.exec_()

    def setparameter(self, mol):
            params = [
                {'name': mol.name, 'type': 'group', 'children': [
                    {'name': 'File', 'type': 'str', 'value': mol.filename, 'readonly': True},
                    {'name': 'Formula', 'type': 'str', 'value': mol.stoichiometry(), 'readonly': True},
                    {'name': 'Mass', 'type': 'float', 'value': mol.mass(), 'readonly': True},
                    {'name': 'Electrons', 'type': 'int', 'value': mol.nel(), 'readonly': True},
                    {'name': 'Charge', 'type': 'float', 'value': mol.charge, 'readonly': True}]
                }]

            p = Parameter.create(name='Molecule', type='group', children=params)
            self.parametertree.setParameters(p, showTop=True)


class TableWindow(QtGui.QDialog):
    '''A window with a table widget to show molecule data.

    :param molecule: The molecule for which to show the data.
    :type molecule: molecule.Molecule

    '''
    def __init__(self, mol, parent, **kwargs):
        super(TableWindow, self).__init__(parent)
        kind = kwargs.get('kind', None)
        if kind == 'distance':
            self.setWindowTitle('Molviewer :: Distance matrix :: {:s}'.format(mol.filename))
            data = distance_matrix(mol)
            table = pg.TableWidget(sortable=False)
            table.setData(data)
        # elif kind == 'info':
            # self.setWindowTitle('MolViewer :: Info table :: {:s}'.format(mol.filename))
            # data = moleculeinfo(mol)
            # table = pg.TableWidget(sortable=False)
            # table.setData(data)

        closebutton = QtGui.QPushButton("&Close")
        closebutton.setFocus()
        closebutton.clicked.connect(self.close)

        buttonbox = QtGui.QDialogButtonBox()
        buttonbox.addButton(closebutton, QtGui.QDialogButtonBox.ActionRole)

        verticalLayout = QtGui.QVBoxLayout(self)
        verticalLayout.addWidget(table)
        verticalLayout.addWidget(buttonbox)


class TreeWindow(QtGui.QDialog):
    '''A window to show molecule information in a tree view.

    :param molecule: The molecule for which to show the data.
    :type molecule: molecule.Molecule

    '''
    def __init__(self, mol, parent):
        super(TreeWindow, self).__init__(parent)
        self.setWindowTitle('MolViewer :: Info table :: {:s}'.format(mol.filename))

        data = moleculeinfo(mol)
        tree = pg.DataTreeWidget()
        tree.setData(data, hideRoot=True)

        closebutton = QtGui.QPushButton("&Close")
        closebutton.setFocus()
        closebutton.clicked.connect(self.close)

        buttonbox = QtGui.QDialogButtonBox()
        buttonbox.addButton(closebutton, QtGui.QDialogButtonBox.ActionRole)

        verticalLayout = QtGui.QVBoxLayout(self)
        verticalLayout.addWidget(tree)
        verticalLayout.addWidget(buttonbox)


class OptionsDialog(QtGui.QDialog):
    def __init__(self, parent):
        super(OptionsDialog, self).__init__(parent)
        self.setWindowTitle('Options')


class AboutDialog(QtGui.QDialog):
    def __init__(self, parent):
        super(AboutDialog, self).__init__(parent)
        self.setWindowTitle('About')


class SpectrumWindow(QtGui.QMainWindow):
    def __init__(self, parent):
        super(SpectrumWindow, self).__init__(parent)
        self.initUI()

    def initUI(self):
        self.setWindowTitle('MolViewer :: Spectrum window')
        self.resize(800, 300)

        cw = QtGui.QWidget()
        self.setCentralWidget(cw)
        layout = QtGui.QVBoxLayout()
        cw.setLayout(layout)

        plot = pg.PlotWidget()

        closebutton = QtGui.QPushButton("&Close")
        closebutton.setFocus()
        closebutton.clicked.connect(self.close)

        buttonbox = QtGui.QDialogButtonBox()
        buttonbox.addButton(closebutton, QtGui.QDialogButtonBox.ActionRole)

        layout.addWidget(plot)
        layout.addWidget(buttonbox)


def main():
    '''The main function

    Args: None
    Returns: None
    '''
    app = QtGui.QApplication(sys.argv)
    app.setWindowIcon(QtGui.QIcon('../icons/flask-256.png'))
    viewer = MolViewer()
    mol = read_xyz('../molecules/xyz/methane.xyz')
    viewer._open(mol)
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
