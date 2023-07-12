import datetime
from netCDF4 import Dataset


def add_global_attributes(filename, year, month, day, product='unknown'):
    """
    Add global attributes to a NetCDF file.

    Parameters:
    filename (str): The path to the NetCDF file.
    year (int): The origin year.
    month (int): The origin month.
    day (int): The origin day.
    product (str, optional): The product name. Defaults to 'unknown'.

    Returns:
    None
    """

    # Ensure the input file exists
    try:
        nc = Dataset(filename, 'a')
    except FileNotFoundError:
        print(f"File {filename} not found.")
        return

    # Assign default time values for hour, minute, and second
    hour = 0
    minute = 0
    second = 0

    # Create a datetime object from the input arguments
    origin_date = datetime.datetime(year, month, day, hour, minute, second)

    # Add global attributes
    nc.origin_date = origin_date.strftime('%d-%b-%Y %H:%M:%S')
    nc.product = product

    # Close the NetCDF file
    nc.close()
