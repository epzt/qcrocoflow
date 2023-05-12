import numpy as np
from netCDF4 import Dataset
from math import pi
from scipy.interpolate import griddata

def add_topo_file(grdname, topofile):
    # Lire la grille CROCO
    grd = Dataset(grdname, 'r')
    lon = grd.variables['lon_rho'][:]
    lat = grd.variables['lat_rho'][:]
    pm = grd.variables['pm'][:]
    pn = grd.variables['pn'][:]
    grd.close()

    # Obtenir la résolution moyenne de CROCO
    dx = np.mean(1 / pm)
    dy = np.mean(1 / pn)
    dx_croco = np.mean([dx, dy])
    print(f"   Résolution CROCO : {dx_croco / 1000:.3f} km")

    dl = max([1, 2 * (dx_croco / (60 * 1852))])
    lonmin = np.min(lon) - dl
    lonmax = np.max(lon) + dl
    latmin = np.min(lat) - dl
    latmax = np.max(lat) + dl

    # Ouvrir le fichier topo
    nc = Dataset(topofile, 'r')
    tlon = nc.variables['lon'][:]
    tlat = nc.variables['lat'][:]

    # Obtenir une sous-grille
    j = np.where((tlat >= latmin) & (tlat <= latmax))
    i1 = np.where((tlon - 360 >= lonmin) & (tlon - 360 <= lonmax))
    i2 = np.where((tlon >= lonmin) & (tlon <= lonmax))
    i3 = np.where((tlon + 360 >= lonmin) & (tlon + 360 <= lonmax))
    x = np.concatenate((tlon[i1] - 360, tlon[i2], tlon[i3] + 360))
    y = tlat[j]

    topo = []
    if i2:
        topo = -nc.variables['topo'][j[0], i2[0]]
    else:
        topo = []

    if i1 and i1[0].size > 0:
        X, Y = np.meshgrid(tlon[i1] - 360, tlat[j])
        topo1 = -nc.variables['topo'][Y.ravel().astype(int), X.ravel().astype(int)].reshape(Y.shape)
        topo = np.hstack((topo1, topo))

    if i3 and i3[0].size > 0:
        X, Y = np.meshgrid(tlon[i3] + 360, tlat[j])
        topo3 = -nc.variables['topo'][Y.ravel().astype(int), X.ravel().astype(int)].reshape(Y.shape)
        topo = np.hstack((topo, topo3))

    nc.close()
    # Obtenir la résolution moyenne de TOPO
    R = 6367442.76
    deg2rad = np.pi / 180

    dg = np.mean(x[1:] - x[:-1])
    dphi = y[1:] - y[:-1]

    dy = R * deg2rad * dphi
    dx = R * deg2rad * dg * np.cos(deg2rad * y[:-1])

    dx_topo = np.mean([dx, dy])
    print("Taille de topo:", topo.shape)
    print(f"   Résolution des données topographiques : {dx_topo / 1000:.3f} km")

    # Dégrader la résolution TOPO
    n = 0
    while dx_croco > dx_topo:
        n += 1

        x = 0.5 * (x[1:] + x[:-1])
        x = x[::2]

        y = 0.5 * (y[1:] + y[:-1])
        y = y[::2]

        topo = 0.25 * (topo[1:, :-1] + topo[1:, 1:] + topo[:-1, :-1] + topo[:-1, 1:])
        topo = topo[::2, ::2]

        dg = np.mean(x[1:] - x[:-1])
        dphi = y[1:] - y[:-1]

        dy = R * deg2rad * dphi
        dx = R * deg2rad * dg * np.cos(deg2rad * y[:-1])

        dx_topo = np.mean([dx, dy])

    print("Taille de dx:", dx.shape)
    print("Taille de dy:", dy.shape)
    print("Taille de topo:", topo.shape)
    
    print(f"   Résolution topographique réduite de moitié {n} fois")
    print(f"   Nouvelle résolution topographique : {dx_topo / 1000:.3f} km")




    # Ajuster les dimensions de x et y pour correspondre à celles de lon et lat
    x, y = np.meshgrid(x, y)

    # Interpoler le topo
    h = griddata((x.ravel(), y.ravel()), topo.ravel(), (lon, lat), method='cubic')

    return h


  

