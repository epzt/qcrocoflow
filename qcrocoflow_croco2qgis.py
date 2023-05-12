# -*- coding: utf-8 -*-
"""
/***************************************************************************
 qcrocoflowDockWidget
                                 A QGIS plugin
 Plugin to manage CROCO intial and result files
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                             -------------------
        begin                : 2022-04-12
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
import netCDF4 as nc
import difflib
import numpy as np
from datetime import datetime, timedelta

from qgis.PyQt import QtGui, QtWidgets, uic, Qt
#from qgis.PyQt.QtWidgets import QDialog
from qgis.PyQt.QtCore import pyqtSignal

from PyQt5.QtWidgets import QDialog, QMessageBox, QFileDialog, QVBoxLayout, QComboBox
from .qcrocoflow_config import *
from .qcrocoflow_timereference import qcrocoflowTIMEREFERENCE

FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'qcrocoflow_croco2qgis_dialog.ui'), resource_suffix='') # Change 13/05/2022
class qcrocoflowCROCO2QGIS(QDialog, FORM_CLASS):
    def __init__(self, currentDir=None, parent=None):
    # Constructor of qcrocoflowCROCO2QGIS
        super(qcrocoflowCROCO2QGIS, self).__init__(parent)
        # Final Setup
        self.setupUi(self)
        # Variables initialisations
        self.innetCDFFileName = None
        self.currentDirectory = currentDir
        self.dVars = dict()
        self.dCoords = dict()

        self.innetCDFFileNameButton.clicked.connect(self.SelectNetCDFInFile)
        self.innetCDFStartDateDate.dateTimeChanged.connect(self.StartDateTimeChanged)
        self.innetCDFEndDateDate.dateTimeChanged.connect(self.EndDateTimeChanged)
        self.lonVariableNameComboBox.currentTextChanged.connect(self.LongitudeVariableCoordsChanged)
        self.latVariableNameComboBox.currentTextChanged.connect(self.LatitudeVariableCoordsChanged)

    def GetVariablesWithDim(self, _dataset, _dim) -> list:
    # Return a list of variables width _dim dimensions """
        retList = list()
        for v in _dataset.variables:
            if _dataset[v].ndim == _dim:
                retList.append(v)
        return retList

    #def IsVariablesWithTime(self, _var, _dimtimename) -> list:
    # Return True or False if _var has _dimtimename in its dimensions """
    #    return _dimtimename in _var.dimensions

    def IsExistingDimension(self, _dataset, _dim, _ratio=DEFAULTRATIO) -> tuple:
    # Search for time dimmension in the netCDF file """
        dims = [d for d in _dataset.dimensions]
        maxMatch = 0.0
        bestdName = None
        dSize = 0
        for _, dname in enumerate(dims):
            match = difflib.SequenceMatcher(None, dname, _dim)
            ratio = match.ratio()
            if ratio > maxMatch and ratio > _ratio:
                maxMatch = ratio
                dSize = _dataset.dimensions[dname].size
                bestdName = dname
        return (maxMatch > 0, bestdName, dSize)

    def StartDateTimeChanged(self, _date) -> None:
        if _date > self.innetCDFEndDateDate.dateTime():
            return
        else:
            self.innetCDFStartDateDate.setDateTime(_date)
        return

    def EndDateTimeChanged(self, _date) -> None:
        if _date < self.innetCDFStartDateDate.dateTime():
            return
        else:
            self.innetCDFEndDateDate.setDateTime(_date)
        return

    def DeleteLayout(self, _layout):
        if _layout is not None:
            while _layout.count():
                item = _layout.takeAt(0)
                widget = item.widget()
                if widget is not None:
                    widget.deleteLater()
                else:
                    self.DeleteLayout(item.layout())

    def UpdateVariablesComboBox(self, _ldicts: list) -> None:
    # Create combobox and fill them according variable lists
        # Fill the coordinates comboboxes with 2D variables
        self.lonVariableNameComboBox.addItems(_ldicts[0]["vars"])
        self.latVariableNameComboBox.addItems(_ldicts[0]["vars"])
        if self.latVariableNameComboBox.findText('lat_rho') < 0 or self.lonVariableNameComboBox.findText('lon_rho') < 0:
            QMessageBox.information(self, "Coordinates warning", "Select proper coordinate variables.")
        else:
            self.lonVariableNameComboBox.setCurrentIndex(self.lonVariableNameComboBox.findText('lon_rho'))
            self.latVariableNameComboBox.setCurrentIndex(self.latVariableNameComboBox.findText('lat_rho'))
            self.dCoords['lat'] = "lat_rho" # Default values for coordinates
            self.dCoords['lon'] = "lon_rho"
        # Delete eventually previous comboboxes
        self.DeleteLayout(self.variablesGroupBox.layout())
        if self.variablesGroupBox.layout() is not None:
            gbox = self.variablesGroupBox.layout()
        else:
            gbox = QVBoxLayout()
        cbs = list()
        for d in _ldicts:
            if len(d["vars"]) == 0:
                continue
            cbs.append(QComboBox(self))
            cbs[-1].setAccessibleName(str(d["dim"]))
            cbs[-1].addItem("(var{}D)".format(d["dim"]))
            cbs[-1].addItems(d["vars"])
            cbs[-1].currentTextChanged.connect(self.VariableChanged)
            gbox.addWidget(cbs[-1])
        gbox.addStretch(10)
        self.variablesGroupBox.setLayout(gbox)
        return

    def LongitudeVariableCoordsChanged(self, _str: str) -> None:
    # Create/update the dictionary of coordinates selected
        self.dCoords["lon"] = _str
    def LatitudeVariableCoordsChanged(self, _str: str) -> None:
    # Create/update the dictionary of coordinates selected
        self.dCoords["lat"] = _str
    def VariableChanged(self, _str: str) -> None:
    # Create/update the dictionary of variables selected for each dimension of the netcdf file
        self.dVars[self.sender().accessibleName()] = _str

    def GetTimeReference(self, _minT, _maxT) -> tuple:
    # Ask user for reference time of the netCDF file
        dlg = qcrocoflowTIMEREFERENCE(_minT, _maxT, self)
        dlg.exec()
        return dlg.GetMinRefrenceTime(), dlg.GetMaxRefrenceTime()

    def UpdateStartEndDate(self, _min, _max) -> None:
    # Update calendar widget according time """
        self.innetCDFStartDateDate.setDateTime(_min)
        self.innetCDFStartDateDate.setDateRange(_min, _max)
        self.innetCDFEndDateDate.setDateTime(_max)
        self.innetCDFEndDateDate.setDateRange(_min, _max)
        return

    def GetInformationFromNetCDFFile(self, _file: str) -> list:
    # Create/update widgets of the dialog according variables found"""
        try:
            dataset = nc.Dataset(_file, 'r')
        except:
            QMessageBox.warning(self, "Error", f"Can not open/access file {self.innetCDFFileName}")
            return []
        # look if a time variable is present in the netCDF dataset
        present, tName, nbrTimeSteps = self.IsExistingDimension(dataset, 'time')
        if present:
            minTime = np.nanmin(dataset[tName][:])
            maxTime = np.nanmax(dataset[tName][:])
            minDate = datetime.utcfromtimestamp(minTime)
            maxDate = datetime.utcfromtimestamp(maxTime)
            minDate, maxDate = self.GetTimeReference(minDate, maxDate)
            self.UpdateStartEndDate(minDate, maxDate)
        else:
            QMessageBox.information(self, "Information", f"{self.innetCDFFileName} as no time dimension")
            minDate = datetime.fromtimestamp(0)
            maxDate = datetime.fromtimestamp(0)
        # Normaly 2D variables have no time dimension
        dVar2D = {"dim":2, "vars":self.GetVariablesWithDim(dataset, 2)}
        # 3D variables may or not have time dimension
        dVar3D = {"dim":3, "vars":self.GetVariablesWithDim(dataset, 3)}
        # Normaly all 4D variables have a time dimension
        dVar4D = {"dim":4, "vars":self.GetVariablesWithDim(dataset, 4)}
        return [dVar2D, dVar3D, dVar4D]

    def SelectNetCDFInFile(self) -> tuple:
    # Allows the user to choose a netCDF file to read
        dialog = QFileDialog(self)
        dialog.setFileMode(QFileDialog.ExistingFile)
        dialog.setNameFilter("netCDF (*.nc *.NC)")  # Just for dev purpose -> TODO: eliminate asap
        dialog.setWindowTitle("Open CROCO result netCDF file")
        dialog.setDirectory(self.currentDirectory)
        dialog.setViewMode(QFileDialog.Detail)
        if not dialog.exec():
            return
        selectedFileName = dialog.selectedFiles()[0] # Get the fisrt element of the returned list
        if not selectedFileName:
            return
        # Updating variables according last choice
        self.innetCDFFileName = os.path.basename(selectedFileName)
        self.currentDirectory = os.path.dirname(selectedFileName)
        # Updating of linedit on the dialog
        if self.innetCDFFileName == "":
            self.innetCDFFileNameLineEdit.setText("")
            return
        self.innetCDFFileNameLineEdit.setText(self.innetCDFFileName)
        # Get dates (if present) and dictionary of variables present in the netCDF file
        ldVars = self.GetInformationFromNetCDFFile(selectedFileName)
        self.UpdateVariablesComboBox(ldVars)
        return

    def GetNetCDFFileName(self) -> tuple:
        if self.currentDirectory is not None and self.innetCDFFileName is not None:
            return os.path.join(self.currentDirectory, self.innetCDFFileName)
        else:
            return ""

    def GetDictOfVars(self) -> dict:
        return self.dVars

    def GetDictOfCoords(self) -> dict:
        return self.dCoords
