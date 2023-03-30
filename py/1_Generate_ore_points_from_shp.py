#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 14:19:23 2020

@author: geolog
"""

import gdal
import ogr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
import os #file system tools

from scipy.sparse import csr_matrix #https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html
from scipy.ndimage.morphology import binary_dilation
from skimage.morphology import dilation
from skimage.morphology import disk
from skimage.draw import circle

from mygdal_functions0_9 import *
from configuration import *

#some solutions were taken from
##https://pcjericks.github.io/py-gdalogr-cookbook/vector_layers.html 



ColMinInd=0; RowMinInd=0; #because we work on already cropped pictures

#2 Чтение extent AOI shp
try:
    aoi_ext=get_shp_extent(shpfilepath); #x_min, x_max, y_min, y_max;
    print(aoi_ext)
except:
    print("Can not read shp AOI file.")
    
#3 Чтение Points shp
try:
    pnts_ext=get_shp_extent(points_shp_filepath);
    print(pnts_ext)
except:
    print("Can not read shp Points file.")

#3.5 Получение разрешения пикселя в метрах из точечного слоя
try:
    for file in os.listdir(dir_cropped_path):         #exclude 3band tif files
        if file.lower().endswith("."+fileext.lower()): 
           gdal_object = gdal.Open(os.path.join(dir_cropped_path,file)) #as new gdal_object was created, no more ColMinInd,RowMinInd
           bands1=gdal_object.GetRasterBand(1).ReadAsArray() ; 
           
           #pixel_x_size,pixel_y_size=GetResolutionMeters(gdal_object);  #неверное разрешени
           pixel_x_size,pixel_y_size=gdal_object.GetGeoTransform()[1],gdal_object.GetGeoTransform()[5];
           print('x res=%.4f, y res=%.4f'%(pixel_x_size,pixel_y_size));
           break;
           
except(FileNotFoundError):
        print("Input image folder doesn\'t exist...");  
       
   

    
#4 Чтение объектов из точечного слоя
driver = ogr.GetDriverByName('ESRI Shapefile')

dataSource = driver.Open(points_shp_filepath, 0) # 0 means read-only. 1 means writeable.
#количество объектов в слое
layer = dataSource.GetLayer()
projection=layer.GetSpatialRef();
featureCount = layer.GetFeatureCount()
print ("Number of features in " + str(featureCount))

#5Create points list inside shp AOI
points_inside=list(); 
for feature in layer:
     #shp_ext 
    geom=feature.GetGeometryRef();
    mx,my=geom.GetX(), geom.GetY()  #coord in map units
    print("point on the map:")
    print(mx,my);
    if (mx<=aoi_ext[1]) and (mx>=aoi_ext[0]) and (my>=aoi_ext[2]) and (my<=aoi_ext[3]):
        points_inside.append((mx,my));
#convert global coordinates to local
lpntx,lpnty=list(),list(); #points list with local coordinates
for point in points_inside:
    lpntx.append(point[0]-aoi_ext[0]); 
    lpnty.append(aoi_ext[3]-point[1]);  #meters into pixels
    #if Landsat coordinates grow from top to bottom aoi_ext[3]-point[1], if otherwise point[1] - aoi_ext[2]
    
    

#graphic output
plt.plot(lpntx,lpnty,'o')
plt.show()

#finding resolution
x_res = int((aoi_ext[1] - aoi_ext[0]) / abs(pixel_x_size));
y_res = int((aoi_ext[3] - aoi_ext[2]) / abs(pixel_y_size));

#6 Points raster map 
#making sparse matrix

pntx=np.int32(np.array(lpntx)/abs(pixel_x_size));
pnty=np.int32(np.array(lpnty)/abs(pixel_y_size));
pntval=np.ones(pnty.size);

#обязательно shape=(y_res-1, x_res-1)!!!!!!!!!!!!!!
#obj_raster=csr_matrix((pntval, (pnty, pntx)), shape=(y_res-1, x_res-1)).toarray(); 
obj_raster=np.zeros([y_res-1, x_res-1]);

for r,c in zip(pnty,pntx):
    try:
        obj_raster=draw_circle(obj_raster,r,c,radius);
    except IndexError:
        print('Index error, possibly marginal position of point?');

#7 Output into Geotiff
#data conversion
obj_raster[np.isnan(obj_raster)==True]=65536
obj_raster=np.uint16(obj_raster)


#saveGeoTiff(obj_raster,raster_path,gdal_object,ColMinInd,RowMinInd);    
driver = gdal.GetDriverByName("GTiff")
outdata = driver.Create(raster_path, x_res, y_res, 1, gdal.GDT_UInt16)
outdata.SetGeoTransform((aoi_ext[0], pixel_x_size, 0, aoi_ext[3], 0, pixel_y_size));
outdata.SetProjection(projection.ExportToWkt())##sets same projection as input
outdata.GetRasterBand(1).WriteArray(obj_raster);


outdata.GetRasterBand(1).SetNoDataValue(0);
outdata.FlushCache();

#8 Image output
plt.imshow(obj_raster)
plt.colorbar()
plt.show()

plt.imshow(bands1)
plt.colorbar()
plt.show()