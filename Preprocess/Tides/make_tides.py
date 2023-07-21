from datetime import datetime, timedelta
from .TPXO2CROCO_Script import TPXO2CROCO
from PyQt5.QtWidgets import QMessageBox
import os
def make_tides_script(TPXO_path, grdname_full_path, ROMSnames, simudate, date_min_datetime, grdname, directory):
    QMessageBox.warning(None, "Warning",
                        "Warning: Please note, the generation of forcing files may take some time. Do not quit QGIS, and wait for the loading to finish. You will be notified once the process is complete.")
    # Define initial simulation date
    t0 = (date_min_datetime - datetime(1, 1, 1)).days + 367

    # Define simulation length in days
    lengthSim = simudate  # Approximate length of model run in days. TODO: replace with difference between simulation dates
    #grdnamefullpath = grdname_full_path.replace("/", "\\")

    fnGrid = r'{}'.format(grdname_full_path)
    # Define paths to grid and output files

    fnOut = os.path.join(directory, grdname + "_tides.nc")
    print(fnOut)
    # Define list of harmonique
    ROMSnames= ROMSnames


    #ROMSnames = ['m2', 'n2']
    # Define path to TPXO files
    TPXO_path = TPXO_path

    #Run
    TPXO2CROCO(t0, ROMSnames, fnGrid, fnOut, lengthSim, TPXO_path, grdname)
