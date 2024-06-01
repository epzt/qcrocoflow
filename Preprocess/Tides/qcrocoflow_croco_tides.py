"""
/***************************************************************************
 qcrocoflow_croco_grid_dialog
                                 Aa QGIS plugin
 qcrocoflow_crocogrid
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                             -------------------
        begin                : 2024-05-25
        git sha              : $Format:%H$
        copyright            : by Jonathan nejmann & Emmanuel Poizot
        emails               : jonathan.nejmann@hotmail.fr
                               emmanuel.poizot@lecnam.net
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

from datetime import datetime, timedelta
from .qcrocoflow_croco_TPXO import qcrocoflow_croco_TPXO
from PyQt5.QtWidgets import QMessageBox
import os

class qcrocoflow_qcroco_tides():
    def __init__(self, _parent=None):
        self.parentApp = _parent
    @classmethod
    def make_tides_script(self, TPXO_path, grdname_full_path, ROMSnames, simudate, date_min_datetime, grdname, directory):
        QMessageBox.warning(None, "Warning",
                            "Warning: Please note, the generation of forcing files may take some time. Do not quit QGIS, and wait for the loading to finish. You will be notified once the process is complete.")
        # Define initial simulation date
        t0 = (date_min_datetime - datetime(1, 1, 1)).days + 367

        # Define simulation length in days
        lengthSim = simudate  # Approximate length of model run in days. TODO: replace with difference between simulation dates
        grdnamefullpath = grdname_full_path.replace("/", "\\")
        fnGrid = r'{}'.format(grdnamefullpath)
        # Define paths to grid and output files

        fnOut = os.path.join(directory, grdname + "_tides.nc")
        print(fnOut)
        # Define list of harmonique
        ROMSnames= ROMSnames


        #ROMSnames = ['m2', 'n2']
        # Define path to TPXO files
        TPXO_path = TPXO_path

        #Run
        qcrocoflow_croco_TPXO.TPXO2CROCO(t0, ROMSnames, fnGrid, fnOut, lengthSim, TPXO_path, grdname)