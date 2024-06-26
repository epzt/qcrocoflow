# Standard imports
"""
/***************************************************************************
 qcrocoflow_croco_grid_dialog
                                 Aa QGIS plugin
 qcrocoflow_crocogrid
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                             -------------------
        begin                : 2024-05-25
        git sha              : $Format:%H$
        copyright            : by Emmanuel Poizot
        email                : emmanuel.poizot@lecnam.net
 ***************************************************************************/
This module provides a function to create a regular grid from the coordinates
minimum and maximum latitude and longitude provided.
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

from netCDF4 import Dataset
import numpy as np
from datetime import date
import os
# QT5 imports
from PyQt5 import uic
from PyQt5.QtWidgets import QDialog, QCheckBox, QGridLayout, QStackedWidget, QVBoxLayout, QWidget, QMessageBox, QFileDialog
from PyQt5.QtCore import Qt
# Specific imports
import copernicusmarine

FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'qcrocoflow_ini_copernicus_dialog.ui'))
class qcrocoflow_Copernicus(QDialog, FORM_CLASS):
    def __init__(self, _argDict, parent=None):
        super(qcrocoflow_Copernicus, self).__init__(parent)
        self.setupUi(self)
        # Initialisation of Copernicus products list
        productNames = ["cmems_mod_ibi_phy_anfc_0.027deg-2D_PT15M-i",
                        "cmems_mod_ibi_phy_anfc_0.027deg-2D_PT1H-m",
                        "cmems_mod_ibi_phy_anfc_0.027deg-3D_PT1H-m",
                        "cmems_mod_ibi_phy_anfc_0.027deg-3D_P1D-m",
                        "cmems_mod_ibi_phy_anfc_0.027deg-3D_P1M-m"]
        self.productListCopernicusComboBox.insertItems(0,productNames)
        self.productListCopernicusComboBox.setCurrentIndex(3)
        procductVariables = [["time","longitude","latitude","zos","uo","vo"],
                             ["time","longitude","latitude","thetao","uo","vo","ubar","vbar","zos","mlotst"],
                             ["time","longitude","latitude","depth","thetao","so","uo","vo"],
                             ["time","longitude","latitude","depth","thetao","so","uo","vo","zos","42ercat","mlotst"],
                             ["time","longitude","latitude","42ercat","depth","mlotst","so","thetao","uo","vo","zos"]]
        # to manage checboxes variables
        self.outputFileName = None
        self.outputDirName = None
        self.stackedWidget = QStackedWidget()
        self.productPageList = list()
        self.productPageList.append(self.CreateCheckBoxVariable(procductVariables[0]))
        self.stackedWidget.addWidget(self.productPageList[-1])
        self.productPageList.append(self.CreateCheckBoxVariable(procductVariables[1]))
        self.stackedWidget.addWidget(self.productPageList[-1])
        self.productPageList.append(self.CreateCheckBoxVariable(procductVariables[2]))
        self.stackedWidget.addWidget(self.productPageList[-1])
        self.productPageList.append(self.CreateCheckBoxVariable(procductVariables[3]))
        self.stackedWidget.addWidget(self.productPageList[-1])
        self.productPageList.append(self.CreateCheckBoxVariable(procductVariables[4]))
        self.stackedWidget.addWidget(self.productPageList[-1])
        layout = QVBoxLayout()
        layout.addWidget(self.stackedWidget)
        self.variableCheckboxGroupBox.setLayout(layout)
        # Connections
        self.productListCopernicusComboBox.currentIndexChanged.connect(self.SwitchPage)
        self.downloadCopernicusFilePushButton.clicked.connect(self.GetCopernicusDataFile)
        self.outputDirCopernicusPushButton.clicked.connect(self.SelectOutputCopernicusDirName)
        self.outputCopernicusFileLineEdit.textChanged.connect(self.SetOutputCopernicusFileName)

    def SetOutputCopernicusFileName(self, _name):
        self.outputFileName = os.path.basename(self.outputCopernicusFileLineEdit.text())
        file_name, file_extension = os.path.splitext(self.outputFileName)
        if file_name == "" and file_extension != ".nc":
            return
        self.outputFileName = file_name + ".nc"
        self.outputCopernicusFileLineEdit.setText(self.outputFileName)
        self.outputCopernicusFileLineEdit.setCursorPosition(len(self.outputFileName)-3)
    def SelectOutputCopernicusDirName(self):
        self.outputDirName = QFileDialog.getExistingDirectory(self, "Select a directory", os.path.expanduser("~user"), QFileDialog.ShowDirsOnly | QFileDialog.DontResolveSymlinks)
        self.outputDirNameLabel.setText(self.outputDirName)

    def CreateCheckBoxVariable(self,_variableList):
        page = QWidget()
        checboxList = list()
        i = j = 0
        layout = QGridLayout()
        for v in _variableList:
            checboxList.append(QCheckBox(v))
            if v == "latitude" or v == "longitude" or v == "time":
                checboxList[-1].setChecked(True)  # Lat, lon , time checked by default
                checboxList[-1].setEnabled(False)
            layout.addWidget(checboxList[-1], j, i)
            i += 1
            if i >= 3:
                j += 1
                i = 0
        page.setLayout(layout)
        return page

    def SwitchPage(self):
        self.stackedWidget.setCurrentIndex(self.productListCopernicusComboBox.currentIndex())

    def GetCopernicusDataFile(self):
        thePage = self.stackedWidget.currentWidget()
        variableList = list()
        for i in range(thePage.layout().count()):
            item = thePage.layout().itemAt(i)
            if item.widget().isChecked():
                variableList.append(item.widget().text())
        copernicusmarine.subset(
            dataset_id=self.productListCopernicusComboBox.currentText(),
            variables=variableList,
            minimum_longitude=self.minLongitudeCopernicusLineEdit.text(),
            maximum_longitude=self.maxLongitudeCopernicusLineEdit.text(),
            minimum_latitude=self.minLatitudeCopernicusLineEdit.text(),
            maximum_latitude=self.maxLatitudeCopernicusLineEdit.text(),
            start_datetime=self.startDateCopernicusDateEdit.date().toString(Qt.ISODate),
            end_datetime=self.endDateCopernicusDateEdit.date().toString(Qt.ISODate),
            minimum_depth=self.minDepthCopernicusDoubleSpinBox.value(),
            maximum_depth=self.maxDepthCopernicusDoubleSpinBox.value(),
            output_filename=self.outputFileName,
            output_directory=self.outputDirName
            )
class qcrocoflow_CreateIniFile:
    def __init__(self, _inifilename = None, _gridfilename = None):
        self.inifilename = _inifilename
        self.gridfilename = _gridfilename
        self.thetas = 0.0
        self.thetab = 0.0
        self.hc = 200.0
        self.N = 0
        self.Vtransform = 2

    def GetInifileName(self):
        return self.inifilename
    def SetIniFileName(self, _nn):
        self.inifilename = _nn

    def GetGridfileName(self):
        return self.gridfilename
    def SetGridFileName(self, _nn):
        self.gridfilenamefilename = _nn

    def GetThetaS(self):
        return self.thetas
    def SetThetaS(self, _nv):
        assert(0 <= _nv <= 10)
        self.thetas = _nv

    def GetThetaB(self):
        return self.thetab
    def SetThetaB(self, _nv):
        assert(0 <= _nv <= 4)
        self.thetab = _nv

    def GetHc(self):
        return self.hc
    def SetHc(self, _nv):
        assert(0 <= _nv <= 12000) # no depth higher than 12000m
        self.hc = _nv

    def GetN(self):
        return self.N
    def SetN(self, _nv):
        self.N = _nv

    def GetVtransform(self):
        return self.Vtransform
    def SetVtransform(self, _nv):
        assert(1 <= _nv <= 2) # no depth higher than 12000m
        self.Vtransform = _nv

    def csf(self, sc, theta_s, theta_b):
        if theta_s > 0:
            csrf = (1 - np.cosh(sc * theta_s)) / (np.cosh(theta_s) - 1)
        else:
            csrf = -np.power(sc,2)
        if theta_b > 0:
            h = (np.exp(theta_b * csrf)-1) / (1 - np.exp(-theta_b))
        else:
            h = csrf
        return h

    def scoordinate(self, theta_s, theta_b, N, hc = None, vtransform = 2):
        if vtransform == 0:
            vtransform = 2  # Default new value

        # Set S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
        sc_r = np.zeros(N)
        Cs_r = np.zeros(N)
        sc_w = np.zeros(N+1)
        Cs_w = np.zeros(N+1)

        if vtransform == 2:
            sc_r = (np.arange(0, N) - N - 0.5) / N
            Cs_r = self.csf(sc_r, theta_s, theta_b)
            sc_w[0] = -1.0
            sc_w[N] = 0
            Cs_w[0] = -1.0
            Cs_w[N] =  0
            sc_w[1:N] = (np.arange(1, N) - N) / N
            Cs_w = self.csf(sc_w, theta_s, theta_b)
        else:
            cff1 = 1 / np.sinh(theta_s)
            cff2 = 0.5 / np.tanh(0.5 * theta_s)
            sc_w = (np.arange(0, N) - N) / N
            Cs_w = (1. - theta_b) * cff1 * np.sinh(theta_s * sc_w) + theta_b * (cff2 * np.tanh(theta_s * (sc_w + 0.5)) - 0.5)
            sc_r = (np.arange(0, N) - N - 0.5) / N
            Cs_r = (1. - theta_b) * cff1 * np.sinh(theta_s * sc_r) + theta_b * (cff2 * np.tanh(theta_s * (sc_r + 0.5)) - 0.5)

        return (sc_r, Cs_r, sc_w, Cs_w)

    def create_inifile(self, inifile, gridfile, title, theta_s, theta_b, hc, N, time, vtransform=2, vstreching=4):
        print(" Creating the file : {}".format(inifile))
        if vtransform == 0:
           print(" NO valid value for VTRANSFORM")
           print(" USE TRANSFORM default value vtransform = 2")
           vtransform = 2

        #  Read the grid file
        nc = Dataset(gridfile,"r")
        h = nc.variables["h"][:]
        mask = nc.variables["mask_rho"][:]
        nc.close()

        hmin = min(min(h(mask==1)))
        if vtransform ==1:
            if hc > hmin:
                print(f"hc ({hc}) m) > hmin ({hmin} m)")
        Mp = np.zeros(h)
        Lp = np.zeros(h)
        L = Lp-1
        M = Mp-1
        Np = N+1

        # Create the initial file
        type = "INITIAL file"
        history = "CROCO"
        nc = Dataset(inifile,"w")

        #  Create dimensions
        nc.createDimension("xi_u",L)
        nc.createDimension("xi_v",Lp)
        nc.createDimension("xi_rho",Lp)
        nc.createDimension("eta_u",Mp)
        nc.createDimension("eta_v",M)
        nc.createDimension("eta_rho",Mp)
        nc.createDimension("s_rho",N)
        nc.createDimension("s_w",Np)
        nc.createDimension("tracer",2)
        nc.createDimension("time",None)
        nc.createDimension("one",1)

        #  Create variables
        spherical = nc.createVariable("spherical","str","one")
        Vtransform = nc.createVariable("Vtransform","i4", "one")
        Vstretching = nc.createVariable("Vstretching","i4","one")
        tstart = nc.createVariable("tstart","f8","one")
        tend = nc.createVariable("tend","f8","one")
        theta_s = nc.createVariable("theta_s","f8","one")
        theta_b = nc.createVariable("theta_b","f8","one")
        Tcline = nc.createVariable("Tcline","f8","one")
        hc = nc.createVariable("hc","f8","one")
        sc_r = nc.createVariable("sc_r","f8","s_rho")
        Cs_r = nc.createVariable("Cs_r","f8","s_rho")
        ocean_time = nc.createVariable("ocean_time","f8","time")
        scrum_time = nc.createVariable("scrum_time","f8","time")
        u = nc.createVariable("u","f8","time","s_rho","eta_u","xi_u")
        v = nc.createVariable("v","f8","time","s_rho","eta_v","xi_v")
        ubar = nc.createVariable("ubar","f8","time","eta_u","xi_u")
        vbar = nc.createVariable("vbar","f8","time","eta_v","xi_v")
        zeta = nc.createVariable("zeta","f8","time","eta_rho","xi_rho")
        temp = nc.createVariable("temp","f8","time","s_rho","eta_rho","xi_rho")
        salt = nc.createVariable("salt","f8","time","s_rho","eta_rho","xi_rho")

        #  Create attributes
        Vtransform.long_name = "vertical terrain-following transformation equation"
        Vstretching.long_name = "vertical terrain-following stretching function"

        tstart.long_name = "start processing day"
        tstart.units = "day"

        tend.long_name = "end processing day"
        tend.units = "day"

        theta_s.long_name = "S-coordinate surface control parameter"
        theta_s.units = "nondimensional"

        theta_b.long_name = "S-coordinate bottom control parameter"
        theta_b.units = "nondimensional"

        Tcline.long_name = "S-coordinate surface/bottom layer width"
        Tcline.units = "meter"

        hc.long_name = "S-coordinate parameter, critical depth"
        hc.units = "meter"

        sc_r.long_name = "S-coordinate at RHO-points"
        sc_r.units = "nondimensional"
        sc_r.valid_min = -1
        sc_r.valid_max = 0

        Cs_r.long_name = "S-coordinate stretching curves at RHO-points"
        Cs_r.units = "nondimensional"
        Cs_r.valid_min = -1
        Cs_r.valid_max = 0

        ocean_time.long_name = "time since initialization"
        ocean_time.units = "second"

        scrum_time.long_name = "time since initialization"
        scrum_time.units = "second"

        u.long_name = "u-momentum component"
        u.units = "meter second-1"

        v.long_name = "v-momentum component"
        v.units = "meter second-1"

        ubar.long_name = "vertically integrated u-momentum component"
        ubar.units = "meter second-1"

        vbar.long_name = "vertically integrated v-momentum component"
        vbar.units = "meter second-1"

        zeta.long_name = "free-surface"
        zeta.units = "meter"

        temp.long_name = "potential temperature"
        temp.units = "Celsius"

        salt.long_name = "salinity"
        salt.units = "PSU"

        # Create global attributes
        nc.title = f"{title}"
        nc.date = f"{date.today()}"
        nc.clim_file = f"{inifile}"
        nc.grd_file = f"{gridfile}"
        nc.type = f"{type}"
        nc.history = f"{history}"

        # Compute S coordinates
        sc_r, Cs_r, sc_w, Cs_w = self.scoordinate(theta_s, theta_b, N, hc, vtransform)

        # Write variables
        nc['spherical'][:] = 'T'
        nc['Vtransform'][:] = vtransform
        nc['Vstretching'][:] = vstreching
        nc['tstart'][:] = time
        nc['tend'][:] = time
        nc['theta_s'][:] = theta_s
        nc['theta_b'][:] = theta_b
        nc['Tcline'][:] = hc
        nc['hc'][:] = hc
        nc['sc_r'][:] = sc_r
        nc['Cs_r'][:] = Cs_r
        nc['scrum_time'][0] = time*24*3600
        nc['ocean_time'][0] = time*24*3600
        nc['u'][:] = 0
        nc['v'][:] = 0
        nc['zeta'][:] = 0
        nc['ubar'][:] = 0
        nc['vbar'][:] = 0
        nc['temp'][:] = 0
        nc['salt'][:] = 0

        # Synchronize on disk
        nc.close()
        return


