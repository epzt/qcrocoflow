import netCDF4 as nc
import numpy as np
import processing

# Variables globales
'''S_RHO = 10
S_W = 11
L = 50
M = 50

HC = 1000.0
THETA_S = 0.0
THETA_B = 0.0
TCLINE = 0.0'''
CELLSIZEMETER = 50

# Fichiers
GRDNAME = '/home/epoizot/recherche/dev/QCrocoflow/DEVTESTS/croco_grid_test.nc'
BATHYGRDNAME = 'lid_tati_zero_hydro'

def check_epsg():
        sysref_epsg = iface.mapCanvas().mapSettings().destinationCrs().authid()
        if sysref_epsg != "EPSG:4326":
            print( "EPSG error", f"The reference system of project is {sysref_epsg}.\nChange it to WGS84.")
            return False
        return True
        
class RectangleMapTool(QgsMapTool):
    rectangleDrawn = pyqtSignal(QgsRectangle)
    def __init__(self, canvas):
        super().__init__(canvas)
        self.canvas = canvas
        self.rubberBand = QgsRubberBand(self.canvas, Qgis.GeometryType.Polygon)  # Crée un polygone pour le rectangle
        self.rubberBand.setColor(QColor(0, 255, 0, 100))  # vert semi-transparent
        self.rubberBand.setWidth(2)
        self.startPoint = None
        self.isDrawing = False
        self.currentMapTool = canvas.mapTool() # Save current map tool
        self.rect = QgsRectangle()

    def canvasPressEvent(self, event):
        # Enregistre le point de départ et commence le dessin
        self.startPoint = self.toMapCoordinates(event.pos())
        self.isDrawing = True
        self.rubberBand.reset(Qgis.GeometryType.Polygon)  # Réinitialise le polygone du rectangle

    def canvasMoveEvent(self, event):
        # Met à jour le rectangle pendant que l'utilisateur dessine
        if self.isDrawing:
            endPoint = self.toMapCoordinates(event.pos())
            self.updateRubberBand(self.startPoint, endPoint)

    def canvasReleaseEvent(self, event):
        # Finalise le rectangle lorsque l'utilisateur relâche la souris
        endPoint = self.toMapCoordinates(event.pos())
        self.rect = QgsRectangle(self.startPoint, endPoint)
        #print(f"Rectangle dessiné : {rect.toString()}")
        self.isDrawing = False
        self.rectangleDrawn.emit(self.rect)
        self.canvas.setMapTool(self.currentMapTool)
        
    def updateRubberBand(self, startPoint, endPoint):
        # Dessine un rectangle entre startPoint et endPoint
        self.rect = QgsRectangle(startPoint, endPoint)
        points = [
            QgsPointXY(self.rect.xMinimum(), self.rect.yMinimum()),  # Coin inférieur gauche
            QgsPointXY(self.rect.xMaximum(), self.rect.yMinimum()),  # Coin inférieur droit
            QgsPointXY(self.rect.xMaximum(), self.rect.yMaximum()),  # Coin supérieur droit
            QgsPointXY(self.rect.xMinimum(), self.rect.yMaximum()),  # Coin supérieur gauche
            QgsPointXY(self.rect.xMinimum(), self.rect.yMinimum())   # Retour au début pour fermer le polygone
        ]
        self.rubberBand.setToGeometry(QgsGeometry.fromPolygonXY([points]), None)

    def deactivate(self):
        # Nettoie le rectangle lors de la désactivation
        self.rubberBand.reset(Qgis.GeometryType.Polygon)
        super().deactivate()

def pick_rectangle():
    loop = QEventLoop()
    def on_rectangle_drawn(rect):
        #print(f"Rectangle finalisé : {rect.toString()}")
        x_min = rect.xMinimum()
        y_min = rect.yMinimum()
        x_max = rect.xMaximum()
        y_max = rect.yMaximum()
        loop.quit()  # Quitter la boucle d'attente lorsque le rectangle est dessiné

    # Initialiser et activer l'outil
    canvas = iface.mapCanvas()
    rect_tool = RectangleMapTool(canvas)
    rect_tool.rectangleDrawn.connect(on_rectangle_drawn)  # Connecter le signal
    canvas.setMapTool(rect_tool)

    # Entrer dans la boucle d'attente
    loop.exec_()
    return rect_tool.rect

class grid():
    def __init__(_nbline, _nbcol, _nblevel=1):
        self.nbline = _nbline
        self.nbcol = _nbcol
        self.nblevel = _nblevel
        
    def grid_init(L, M, S_RHO):
        Mp = M + 1
        Lp = L + 1
        S_W = S_RHO + 1
        nw = nc.Dataset(GRDNAME, 'w', format='NETCDF4')
        # Create dimensions
        xi_u = nw.createDimension('xi_u', L)
        eta_u = nw.createDimension('eta_u', Mp)
        xi_v = nw.createDimension('xi_v', Lp)
        eta_v = nw.createDimension('eta_v', M)
        xi_rho = nw.createDimension('xi_rho', Lp)
        eta_rho = nw.createDimension('eta_rho', Mp)
        xi_psi = nw.createDimension('xi_psi', L)
        eta_psi = nw.createDimension('eta_psi', M)
        one = nw.createDimension('one', 1)
        s_rho = nw.createDimension('s_rho', S_RHO)
        s_w = nw.createDimension('s_w', S_W)
        time = nw.createDimension('time', None)

        lon_rho = nw.createVariable("lonr","f8", (eta_rho, xi_rho))
        lon_rho.units = "degrees_east"
        lat_rho = nw.createVariable("latr","f8", (eta_rho, xi_rho))
        lat_rho.units = "degrees_north"
        lon_u = nw.createVariable("lonu","f8", (eta_u, xi_u))
        lon_u.units = "degrees_east"
        lat_u = nw.createVariable("latu","f8", (eta_u, xi_u))
        lat_u.units = "degrees_north"
        lon_v = nw.createVariable("lonv","f8", (eta_v, xi_v))
        lon_v.units = "degrees_east"
        lat_v = nw.createVariable("latv","f8", (eta_v, xi_v))
        lat_v.units = "degrees_north"
        lon_psi = nw.createVariable("lonp","f8", (eta_psi, xi_psi))
        lon_psi.units = "degrees_east"
        lat_psi = nw.createVariable("latp","f8", (eta_psi, xi_psi))
        lat_psi.units = "degrees_north"

        h = nw.createVariable("h", "f8", (eta_rho, xi_rho))
        h.units = "meters"
        h.long_name = "Bathymetry"
        f = nw.createVariable("f", "f8", (eta_rho, xi_rho))
        f.units = "second-1"
        f.long_name = "Coriolis parameter"
        pm = nw.createVariable("pm", "f8", (eta_rho, xi_rho))
        pm.units = "meter-1"
        pm.long_name = "inverse grid spacing in XI"
        pn = nw.createVariable("pn", "f8", (eta_rho, xi_rho))
        pn.units = "meter-1"
        pn.long_name = "inverse grid spacing in ETA"
        dmde = nw.createVariable("dmde", "f8", (eta_rho, xi_rho))
        dmde.long_name = "" # "ETA-derivative of inverse metric factor m, d(1/M)/d(ETA)"
        dndx = nw.createVariable("dndx", "f8", (eta_rho, xi_rho))
        dndx.long_name = "" # XI-derivative of inverse metric factor n, d(1/N)/d(XI)"
        angle = nw.createVariable("angler", "f8", (eta_rho, xi_rho))
        angle.units = "radians"
        angle.long_name = "angle between XI axis and east"
        mask_rho = nw.createVariable("rmask", "i4", (eta_rho, xi_rho))
        mask_rho.flag_values = "0, 1"
        mask_rho.flag_meanings = "land, water"
        mask_u = nw.createVariable("umask", "i4", (eta_u, xi_u))
        mask_u.flag_values = "0, 1"
        mask_u.flag_meanings = "land, water"
        mask_v = nw.createVariable("vmask", "i4", (eta_v, xi_v))
        mask_v.flag_values = "0, 1"
        mask_v.flag_meanings = "land, water"
        mask_psi = nw.createVariable("pmask", "i4",(eta_psi, xi_psi))
        mask_pdi.flag_values = "0, 1"
        mask_psi.flag_meanings = "land, water"
        s_rho = nw.createVariable("s_rho", "f8", s_rho)
        s_rho.long_name = "S-coordinate at RHO-points"
        s_rho.positive = "up"
        s_rho.standard_name = "ocean_s_coordinate"
        s_w = nw.createVariable("s_w", "f8", s_w)
        s_w.long_name = "S-coordinate at W-points"
        s_w.positive = "up"
        s_w.standard_name = "ocean_s_coordinate"
        hc  = nw.createVariable("hc", "f8", one)
        hc.long_name = "S-coordinate parameter, critical depth"
        hc.units = "meter"
        theta_s  = nw.createVariable("theta_s", "f8", one)
        theta_s.long_name = "S-coordinate surface control parameter"
        theta_b = nw.createVariable("theta_b", "f8", one)
        theta_b.long_name = "S-coordinate bottom control parameter"
        Tcline = nw.createVariable("Tcline", "f8", one)
        Tcline.long_name = "S-coordinate surface/bottom layer width"
        Tcline.units = "meter"
        zeta = nw.createVariable("zeta", "f8", (time, eta_rho, xi_rho))
        zeta.units = "meter"
        zeta.long_name = "sea surface height"
        time = nw.createVariable("time", "f8", time)
        time.units = "days since 0001-01-01 00:00:00"
        time.calendar = "gregorian"

        nw.variables['hc'][:] = HC
        nw.variables['s_rho'][:] = S_RHO
        nw.variables['s_w'][:] = S_W
        nw.variables['theta_s'][:] = 0.0
        nw.variables['theta_b'][:] = 0.0
        nw.variables['Tcline'][:] = 0.0
        
        nw.close()
    
    # _rho, _u, _v, _psi : list of respective 2D arrays of coordinates
    def grid_fill_coords(_rho, _u, _v, _psi):
        nw = nc.Dataset(GRDNAME, 'a', format='NETCDF4')
        nw['lat_rho'][:] = _rho[0]
        nw['lon_rho'][:] = _rho[1]
        nw['lat_u'][:] = _u[0]
        nw['lon_u'][:] = _u[1]
        nw['lat_v'][:] = _v[0]
        nw['lon_v'][:] = _v[1]
        nw['lat_psi'][:] = _psi[0]
        nw['lon_psi'][:] = _psi[1]
        nw.close()
        
    def grid_fill_h():
        return
    
    def grid_fill_f():
        return 
    
    def grid_fill_pmpn():
        return
        
    def grid_fill_dmdedndx():
        return

    def grid_fill_mask():
        return
        
        
def reproject_layer(input_layer, target_crs = 'EPSG:4326'):
    output = processing.run('native:reprojectlayer', {
                    'INPUT': input_layer,
                    'TARGET_CRS': target_crs,
                    'OUTPUT' : 'TEMPORARY_OUTPUT'})
    return output['OUTPUT']

def create_regular_grid(_rect, _cell_size):
    x_min = _rect.xMinimum()
    y_min = _rect.yMinimum()
    x_max = _rect.xMaximum()
    y_max = _rect.yMaximum()
    
    dx = x_max - x_min
    dy = y_max - y_min
    
    nb_max_col = int(np.floor(dx / _cell_size))
    nb_max_lin = int(np.floor(dy / _cell_size))
    
    lon_coords = np.linspace(x_min, x_max, nb_max_col)
    lat_coords = np.linspace(y_min, y_max, nb_max_lin)
    
    layer = QgsVectorLayer("Point?crs=EPSG:4326", "Grid_Points", "memory")
    provider = layer.dataProvider()
    provider.addAttributes([QgsField("id", QVariant.Int)])
    layer.updateFields()
    features = list()
    point_id = 0
    for lat in lat_coords:
        for lon in lon_coords:
            point = QgsPointXY(lon, lat)
            feature = QgsFeature()
            feature.setGeometry(QgsGeometry.fromPointXY(point))
            feature.setAttributes([point_id])
            features.append(feature)
            point_id += 1
            features.append(feature)
    provider.addFeatures(features)
    return layer
    
def create_regular_grids(_rect, _cell_size, _delta):
    x_min = _rect.xMinimum()
    y_min = _rect.yMinimum()
    x_max = _rect.xMaximum()
    y_max = _rect.yMaximum()
    # RHO grid
    reg_grid_point = processing.run('qgis:regularpoints', {'EXTENT': f"{x_min}, {x_max}, {y_min}, {y_max}",
                                         'SPACING' : float(_cell_size),
                                         'INSET' : 0.0,
                                         'RANDOMIZE' : False,
                                         'IS_SPACING' : True,
                                         'CRS' : TARGETCRS,
                                         'OUTPUT' : 'TEMPORARY_OUTPUT'
                                        })
    tmpRLayer = reg_grid_point['OUTPUT']
    tmpRLayer.setName('RHO')
    # U grid
    reg_grid_point = processing.run('qgis:regularpoints', {'EXTENT': f"{x_min - _delta}, {x_max + _delta}, {y_min}, {y_max}",
                                         'SPACING' : float(_cell_size),
                                         'INSET' : 0.0,
                                         'RANDOMIZE' : False,
                                         'IS_SPACING' : True,
                                         'CRS' : TARGETCRS,
                                         'OUTPUT' : 'TEMPORARY_OUTPUT'
                                        })
    tmpULayer = reg_grid_point['OUTPUT']
    tmpULayer.setName('U')
    # V grid
    reg_grid_point = processing.run('qgis:regularpoints', {'EXTENT': f"{x_min}, {x_max}, {y_min - _delta}, {y_max + _delta}",
                                         'SPACING' : float(_cell_size),
                                         'INSET' : 0.0,
                                         'RANDOMIZE' : False,
                                         'IS_SPACING' : True,
                                         'CRS' : TARGETCRS,
                                         'OUTPUT' : 'TEMPORARY_OUTPUT'
                                        })
    tmpVLayer = reg_grid_point['OUTPUT']
    tmpVLayer.setName('V')
    # PSI grid
    reg_grid_point = processing.run('qgis:regularpoints', {'EXTENT': f"{x_min + _delta}, {x_max - _delta}, {y_min + _delta}, {y_max - _delta}",
                                         'SPACING' : float(_cell_size),
                                         'INSET' : 0.0,
                                         'RANDOMIZE' : False,
                                         'IS_SPACING' : True,
                                         'CRS' : TARGETCRS,
                                         'OUTPUT' : 'TEMPORARY_OUTPUT'
                                        })
    tmpPLayer = reg_grid_point['OUTPUT']
    tmpPLayer.setName('PSI')
    return tmpRLayer, tmpULayer, tmpVLayer, tmpPLayer

'''grid = processing.run('native:creategrid', {'CRS' : TARGETCRS,
                            'EXTENT' : f"{x_min}, {x_max}, {y_min}, {y_max}",
                            'HOVERLAY' : 0,
                            'HSPACING' : float(cell_size),
                            'OUTPUT' : 'TEMPORARY_OUTPUT',
                            'TYPE' : 1,
                            'VOVERLAY' : 0,
                            'VSPACING' : float(cell_size)
                            })
tmpVLayer = reproject_layer(grid['OUTPUT'])
tmpVLayer.setName('GRID')
QgsProject.instance().addMapLayer(tmpVLayer)
'''
def get_h_from_raster(_vecteur_layer, _raster_layer):
    res = processing.run('native:rastersampling', {'INPUT': _vecteur_layer,
                                        'RASTERCOPY': _raster_layer,
                                        'COLUMN_PREFIX': "DEPTH_",
                                        'OUTPUT': 'TEMPORARY_OUTPUT'
                                        })
    return res['OUTPUT']

# Compte the Corilis factor f
def coriolis_factor(_lat):
    # _lat must be in degrees
    return 2.0 * 7.2921e-05 * np.sin(np.radians(_lat))
    
def getK0e(_targetcrs, _pt):
    # Définir un CRS spécifique (par exemple EPSG:32631 pour UTM Zone 31N)
    #crs = QgsCoordinateReferenceSystem("EPSG:32631")  # Remplace par le code EPSG de ta projection

    # Extraire le facteur d'échelle au méridien central (k0)
    k0 = _targetcrs.factors(QgsPoint(_pt.x(), _pt.y())).arealScale()

    # Obtenir l'ellipsoïde associé au CRS
    ellipsoid_acronym = _targetcrs.ellipsoidAcronym()
    
    # Calculer l'excentricité (e) de l'ellipsoïde
    semi_major = QgsEllipsoidUtils.EllipsoidParameters().semiMajor  # Demi-grand axe
    semi_minor = QgsEllipsoidUtils.EllipsoidParameters().semiMinor # Demi-petit axe

    # Formule de l'excentricité e = sqrt(1 - (b^2 / a^2))
    excentricite = np.sqrt(1 - (np.power(semi_minor,2) / np.power(semi_major,2)))
    return k0, excentricite


#######################################################################################################
#
#######################################################################################################

#x_min, y_min = -1.3280, 49.5117   # Coordonnées géographiques (long/lat) pour le coin inférieur gauche
#x_max, y_max = -1.2130, 49.5956    # Coordonnées géographiques (long/lat) pour le coin supérieur droit
cell_size = CELLSIZEMETER * 8.9993e-6  # Cell size en degre decimal
delta = 0.5 * cell_size

TARGETCRS = QgsCoordinateReferenceSystem('EPSG:4326') # QgsCoordinateReferenceSystem('EPSG:3857')
SOURCECRS = QgsProject.instance().crs() # QgsCoordinateReferenceSystem('EPSG:4326')
tr = QgsCoordinateTransform(SOURCECRS, TARGETCRS, QgsProject.instance())
if not tr.isValid():
    print("Unable to manage geodetic datum of current project.")
    exit()

# Get rectangle drawed by mouse
grid_rectangle = pick_rectangle()
geo_rectangle = tr.transformBoundingBox(grid_rectangle, Qgis.TransformDirection.Forward, handle180Crossover=False)

# Initial bounding box using given coordinates
#geo_rectangle = QgsRectangle()
#geo_rectangle.set(x_min, y_min, x_max, y_max)

# Create new bounding box using a projected geodetic system to deal with metric coordinates
# metric_rectangle = tr.transformBoundingBox(geo_rectangle, Qgis.TransformDirection.Forward, handle180Crossover=False)

# Modify the metric bounding box to have rounded number of rows and columns
geo_rectangle = geo_rectangle.snappedToGrid(cell_size)
# Refresh min and max coordinates

#rho_grid = create_regular_grid(geo_rectangle, cell_size)
#QgsProject.instance().addMapLayer(rho_grid)

rho_grid, u_grid, v_grid, psi_grid = create_regular_grids(geo_rectangle, cell_size, delta)
QgsProject.instance().addMapLayer(rho_grid)
QgsProject.instance().addMapLayer(u_grid)
QgsProject.instance().addMapLayer(v_grid)
QgsProject.instance().addMapLayer(psi_grid)

da = QgsDistanceArea()
ctex = QgsCoordinateTransformContext()
da.setSourceCrs(SOURCECRS, ctex)

x_min = geo_rectangle.xMinimum()
y_min = geo_rectangle.yMinimum()
x_max = geo_rectangle.xMaximum()
y_max = geo_rectangle.yMaximum()
    
ncol = int((x_max - x_min) / cell_size) + 1
nlin = int((y_max - y_min) / cell_size) + 1

print(f"Nbre col: {ncol}\nNbre lig: {nlin}")

'''selectid = [1,2,ncol+1]
rho_grid.select(selectid)
pts = []
for feat in rho_grid.selectedFeatures():
    pts.append(feat.geometry().asPoint())
    
print(f"Distance geographique (Lon): {da.measureLine(pts[0],pts[1])}")
print(f"Distance geographique (Lat): {da.measureLine(pts[0],pts[2])}")
print(f"Distance metrique (X): {da.convertLengthMeasurement(da.measureLine(pts[0],pts[1]), Qgis.DistanceUnit.DistanceMeters)}")
print(f"Distance metrique (Y): {da.convertLengthMeasurement(da.measureLine(pts[0],pts[2]), Qgis.DistanceUnit.DistanceMeters)}")'''


# Create a vector layer witch contains an attribute 'DEPTH_' filled by values from a raster layer
tmp_h_layer = get_h_from_raster(rho_grid, QgsProject.instance().mapLayersByName(BATHYGRDNAME)[0])
tmp_h_layer.setName('H')
QgsProject.instance().addMapLayer(tmp_h_layer)

selectid = [1,2,ncol+1]
rho_grid.select(selectid)
for featr in rho_grid.selectedFeatures():
    pt=featr.geometry().asPoint()
    print("f={} s-1".format(coriolis_factor(pt.y())))
    
psi_grid.select(selectid)
list_pt = list()
for featp in  psi_grid.selectedFeatures():
    pt=featp.geometry().asPoint()
    k0, e = getK0e(TARGETCRS, pt)
    print(f"K0: {k0} - Exen.: {e}")



    
    
