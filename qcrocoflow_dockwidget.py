# -*- coding: utf-8 -*-
"""
/***************************************************************************
 qcrocoflowDockWidget
                                 A QGIS plugin
 A QGIS plugin to manage CROCO projectsqcrocoflow
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                             -------------------
        begin                : 2022-12-03
        git sha              : $Format:%H$
        copyright            : (C) 2022 by Emmanuel Poizot
        email                : emmanuel.poizot@lecnam.net
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

import os

from qgis.PyQt import QtGui, QtWidgets, uic
from qgis.PyQt.QtCore import pyqtSignal

from qgis.PyQt.QtWidgets import QDialog, QMessageBox, QFileDialog, QGridLayout, QComboBox, QMenuBar, QAction

from .qcrocoflow_croco2qgis import qcrocoflowCROCO2QGIS

FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'qcrocoflow_dockwidget_base.ui'))

class qcrocoflowDockWidget(QtWidgets.QDockWidget, FORM_CLASS):

    closingPlugin = pyqtSignal()

    def __init__(self, _iface, parent=None):
        """Constructor."""
        super(qcrocoflowDockWidget, self).__init__(parent)
        # Set up the user interface from Designer.
        # After setupUI you can access any designer object by doing
        # self.<objectname>, and you can use autoconnect slots - see
        # http://doc.qt.io/qt-5/designer-using-a-ui-file.html
        # #widgets-and-dialogs-with-auto-connect

        # Get pointer on QGIS interface
        self.iface  = _iface
        # Construction of menu bar entries and relative actions
        self.myQMenuBar = QMenuBar(self)
        ### Project ################################
        projectMenu = self.myQMenuBar.addMenu('Project')
        newProjectAction = QAction('New', self)
        newProjectAction.triggered.connect(self.NewProject)
        projectMenu.addAction(newProjectAction)
        # ---
        openProjectAction = QAction('Open', self)
        openProjectAction.triggered.connect(self.OpenProject)
        projectMenu.addAction(openProjectAction)
        # ---
        saveProjectAction = QtWidgets.QAction('Save', self)
        saveProjectAction.triggered.connect(self.SaveProject)
        projectMenu.addAction(saveProjectAction)
        # ---
        projectMenu.addSeparator()
        # ---
        closeProjectAction = QAction('Close', self)
        closeProjectAction.triggered.connect(self.CloseProject)
        projectMenu.addAction(closeProjectAction)
        ### Grid ################################
        gridMenu = self.myQMenuBar.addMenu('Grid')
        newGridAction = QAction('New', self)
        newGridAction.triggered.connect(self.NewGrid)
        gridMenu.addAction(newGridAction)
        # ---
        openGridAction = QAction('Open', self)
        openGridAction.triggered.connect(self.OpenGrid)
        gridMenu.addAction(openGridAction)
        # ---
        gridMenu.addSeparator()
        # ---
        initialConditionAction = QAction('IC', self)
        initialConditionAction.triggered.connect(self.printHello)
        gridMenu.addAction(initialConditionAction)
        # ---
        openBoundaryConditionAction = QAction('OBC', self)
        openBoundaryConditionAction.triggered.connect(self.printHello)
        gridMenu.addAction(openBoundaryConditionAction)
        ### Sediment ################################
        sedimentMenu = self.myQMenuBar.addMenu('Sediment')
        newSedimentAction = QAction('New', self)
        newSedimentAction.triggered.connect(self.NewSediment)
        sedimentMenu.addAction(newSedimentAction)
        # ---
        openSedimentAction = QAction('Open', self)
        openSedimentAction.triggered.connect(self.OpenSediment)
        sedimentMenu.addAction(openSedimentAction)
        # ---
        gridMenu.addSeparator()
        # ---

        # TODO: define other menu entries from herre

        #self.setupUi(self)
        # Variables
        self.projectOpened = False
        self.projectName = None
        self.projectDirectory = os.path.expanduser("~user") # Default value for project directory

    def printHello(self) -> None:
        print('Hello')

    def NewProject(self) -> None:
    # Create a new empty QCrocoFlow project
        if self.projectOpened and self.projectName:
            QMessageBox.warning(self, "Project file", f"The QCrocoFlow project {self.projectName} is currently opened.\nClose it before create a new one")
            return
        dialog = QFileDialog(self)
        dialog.setFileMode(QFileDialog.AnyFile)
        dialog.setNameFilter("QCF project (*.qcf *.QCF)")  # This is the good one
        dialog.setWindowTitle("Create a new QCrocoFlow project")
        dialog.setViewMode(QFileDialog.Detail)
        if (dialog.exec()):
            selectedFileName = dialog.selectedFiles()[0]  # Get the fisrt element of the returned list
            if os.path.isfile(selectedFileName):
                ans = QMessageBox.information(self, "Project file exist", f"Do you want to overwrite {selectedFileName} ?", \
                            buttons=QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel, \
                             defaultButton = QMessageBox.StandardButton.Cancel)
                if ans == QMessageBox.StandardButton.Cancel:
                    return
            self.projectName = os.path.basename(selectedFileName)
            self.projectDirectory = os.path.dirname(selectedFileName)
            # TODO: manage creation of empty project file here
            self.projectOpened = True
        else:
            QMessageBox.information(self, "Project file", "No QCrocoFlow project file created")
        return

    def OpenProject(self) -> None:
    # Open an existing QCrocoFlow  project
        if self.projectOpened and self.projectName:
            QMessageBox.warning(self, "Project file", f"The QCrocoFlow project {self.projectName} is currently opened.\nClose it before open a new one")
            return
        dialog = QFileDialog(self)
        dialog.setFileMode(QFileDialog.ExistingFile)
        dialog.setNameFilter("QCF project (*.qcf *.QCF)") # This is the good one
        dialog.setNameFilter("QCF project (*.txt *.TXT)") # Just for dev purpose -> TODO: eliminate asap
        dialog.setWindowTitle("Open QCrocoFlow existing project file")
        dialog.setViewMode(QFileDialog.Detail)
        if (dialog.exec()):
            selectedFileName = dialog.selectedFiles()[0] # Get the fisrt element of the returned list
            self.projectName = os.path.basename(selectedFileName)
            self.projectDirectory = os.path.dirname(selectedFileName)
            self.projectOpened = True
        else:
            QMessageBox.information(self, "Project file", "No QCrocoFlow project file selected")
        return

    def SaveProject(self) -> None:
        QMessageBox.information(self, "Project file", "Not implemented yet")
        return

    def CloseProject(self) -> None:
    # Manage close de project: empty runtime variables, etc.
        if not self.projectOpened:
            QMessageBox.information(self, "Project file", "No QCrocoFlow project yet opened")
            return
        ans = QMessageBox.question(self, "Project file", f"Do you really want to close {self.projectName} project ?", \
                             buttons=QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No, \
                             defaultButton = QMessageBox.StandardButton.No)
        if ans == QMessageBox.StandardButton.Yes:
            # TODO: manage cleanup variables here
            self.projectName = None
            self.projectDirectory = os.path.expanduser("~user")
            self.projectOpened = False
        return

    def NewGrid(self) -> None:
        QMessageBox.information(self, "Project file", "Not implemented yet")
        return
    def OpenGrid(self) -> None:
        grid = qcrocoflowCROCO2QGIS(self.projectDirectory, self)
        grid.show()
        return

    def NewSediment(self) -> None:
        QMessageBox.information(self, "Project file", "Not implemented yet")
        return

    def OpenSediment(self) -> None:
        QMessageBox.information(self, "Project file", "Not implemented yet")
        return

    def closeEvent(self, event)-> None:
        self.closingPlugin.emit()
        event.accept()
