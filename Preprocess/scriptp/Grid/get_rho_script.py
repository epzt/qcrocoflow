from qgis.core import QgsProject
import numpy as np

def get_rho(dl):
    """
    This function calculates the horizontal (x) and vertical (y) distances of rho grid points.

    Parameters:
        dl (float): The grid resolution.

    Returns:
        dx_rho (np.array): The horizontal distances of rho grid points.
        dy_rho (np.array): The vertical distances of rho grid points.
    """

    # Access the 'grid' layer from the current QGIS project instance
    grid_layer = QgsProject.instance().mapLayersByName('grid')[0]

    # Initialize numpy arrays for x and y coordinates with the size of the grid_layer
    features_count = grid_layer.featureCount()
    x_rho, y_rho = np.zeros(features_count), np.zeros(features_count)

    # Iterate over the features (grid points) in the grid layer
    for i, feature in enumerate(grid_layer.getFeatures()):
        # Get the bounding box of the feature (grid point)
        bbox = feature.geometry().boundingBox()

        # Extract the lower left x and y coordinates of the bounding box
        lower_left_x = bbox.xMinimum()
        lower_left_y = bbox.yMinimum()

        # Store the lower left coordinates in the corresponding numpy arrays
        x_rho[i] = lower_left_x
        y_rho[i] = lower_left_y

    # Calculate the distances from the first point
    dx_rho = x_rho - x_rho[0]
    dy_rho = y_rho - y_rho[0]

    # Calculate the number of points along each axis
    num_x_points = int((x_rho.max() - x_rho.min()) / dl) + 1
    num_y_points = int((y_rho.max() - y_rho.min()) / dl) + 1

    # Reshape the arrays to match the grid dimensions
    dx_rho = dx_rho.reshape((num_y_points, num_x_points))
    dy_rho = dy_rho.reshape((num_y_points, num_x_points))

    return dx_rho, dy_rho
