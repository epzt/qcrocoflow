from qgis.analysis import QgsZonalStatistics
from qgis.core import QgsProject, QgsField, QgsCoordinateReferenceSystem
from qgis import processing
import os
import numpy as np

def add_topo_file2(topofile, numColumns, numLines):
    """
    This function interpolates the topography data into the CROCO grid.

    Args:
        topofile (str): topography file.
        numColumns (int): x.shape.
        numLines (int): y.shape.

    Returns:
        h (np.array): Array with the reshaped topography data.
    """

    # Get raster and grid layers from the project
    basename = os.path.basename(topofile)
    raster_layer = QgsProject.instance().mapLayersByName(basename)[0]
    grid_layer = QgsProject.instance().mapLayersByName('grid')[0]

    # Reproject grid layer to WGS84
    crs_wgs84 = QgsCoordinateReferenceSystem("EPSG:4326")
    parameters = {'INPUT': grid_layer, 'TARGET_CRS': crs_wgs84, 'OUTPUT': 'memory:'}
    result = processing.run('native:reprojectlayer', parameters)

    # Check if the re-projection was successful
    if 'OUTPUT' not in result:
        raise Exception('Grid layer re-projection failed')

    grid_layer = result['OUTPUT']

    # Calculate zonal statistics
    zoneStat = QgsZonalStatistics(grid_layer, raster_layer, '', 1, QgsZonalStatistics.Mean)
    zoneStat.calculateStatistics(None)

    grid_layer.startEditing()

    # Rename the 'mean' field to 'topo'
    for field in grid_layer.fields():
        if field.name() == 'mean':
            grid_layer.renameAttribute(grid_layer.fields().indexFromName('mean'), 'topo')

    # Get the field indices
    field_index_x = grid_layer.fields().indexFromName("lower_left_x")
    field_index_y = grid_layer.fields().indexFromName("lower_left_y")
    field_index_topo = grid_layer.fields().indexFromName("topo")

    # Initialize numpy arrays
    features_count = grid_layer.featureCount()
    lat_rho, lon_rho, topo = np.zeros(features_count), np.zeros(features_count), np.zeros(features_count)

    # Update the fields with the coordinates and topography values
    for i, feature in enumerate(grid_layer.getFeatures()):
        bbox = feature.geometry().boundingBox()
        lower_left_x, lower_left_y = bbox.xMinimum(), bbox.yMinimum()

        topo[i] = 100 if feature['topo'] is None else feature['topo']
        lat_rho[i], lon_rho[i] = lower_left_y, lower_left_x

        grid_layer.changeAttributeValue(feature.id(), field_index_x, lower_left_x)
        grid_layer.changeAttributeValue(feature.id(), field_index_y, lower_left_y)
        grid_layer.changeAttributeValue(feature.id(), field_index_topo, topo[i])

    grid_layer.commitChanges()

    # Generate the grid and reshape the topography array
    h = np.reshape(topo, (numLines[0], numColumns[0]))

    return h
