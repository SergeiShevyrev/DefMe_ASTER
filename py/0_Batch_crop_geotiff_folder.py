# DefMe scripts set

#some methods are from https://pcjericks.github.io/py-gdalogr-cookbook/vector_layers.html
#this script has been composed and written by Sergei L Shevirev http://lefa.geologov.net

import gdal,ogr #OpenGIS Simple Features Reference Implementation
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import time;
import copy;

from skimage.transform import resize #function for SRTM resize according to Landsat resolution
import elevation
import richdem as rd

from mygdal_functions0_9 import *
from configuration import *

"""
libraries version required:
    gdal 2.3.3 from conda install
    scikit-learn 0.23
    richdem 0.3.4 from pip install 
    elevation 1.1.3  pip install
"""
time_start=time.time()

#topocorrection flag
is_topocorrection=True; #defaut on 

#files for processing, input and output directory
# are taking from configuration.py

#no edits beyong that line ->

file_for_crop=[];
metadata=[];

try:
    for file in os.listdir(os.path.join('..',out_dir_name)):
        #file=file.lower();
        if file.lower().endswith("."+fileext.lower()):
            file_for_crop.append(file);
            print(file+" was added to crop queue.");
        if file.lower().endswith("."+metafileext.lower()):
            print(file+" was defined as metadata file.");
            try:
                metadata=parse_MTL(os.path.join('..',out_dir_name,file));
            except:
                print('Unable to parse metadata from {}'.format(os.path.join('..',out_dir_name,file)));
except(FileNotFoundError):
        print("Input image folder doesn\'t exist...");

SunElevation=metadata['SUN_ELEVATION'];      
SunAzimuth=metadata['SUN_AZIMUTH']; 
SolarZenith=90-SunElevation;

#STEP 0. Prepare for the topocorrection
try:
    shp_extent=get_shp_extent(shpfilepath);
except:
    print("Can not read shp AOI file.")
            
#crop dem file
if is_topocorrection==True:
    print("Perform cropping of srtm");
    
    #read DEM geotiff
    srtm_gdal_object = gdal.Open(drm_filepath)
    srtm_band = srtm_gdal_object.GetRasterBand(1)
    srtm_band_array = srtm_band.ReadAsArray() 
    
    #get spatial resolution
    srtm_gt=srtm_gdal_object.GetGeoTransform()
    srtm_xsize = srtm_gdal_object.RasterXSize
    srtm_ysize = srtm_gdal_object.RasterYSize #x and y raster size in pixels
    srtm_ext=GetExtent(srtm_gt,srtm_ysize,srtm_xsize); 
    #resolution in meters
    srtm_dpx=(srtm_ext[3][0]-srtm_ext[0][0])/srtm_xsize
    srtm_dpy=(srtm_ext[0][1]-srtm_ext[2][1])/srtm_ysize
    
    if check_shp_inside_raster(srtm_ext,shp_extent):
#        sampleSrtmImage,ColMinIndSRTM,RowMinIndSRTM =crop_by_shp(shp_extent,srtm_ext,\
#                                                    srtm_dpx,srtm_dpy,srtm_band_array); 
        srtm_band = rd.LoadGDAL(drm_filepath);

        slope = rd.TerrainAttribute(srtm_band, attrib='slope_degrees')
        aspect = rd.TerrainAttribute(srtm_band, attrib='aspect')
        
        if os.path.exists(os.path.join("..","tmp"))==False:
            os.mkdir(os.path.join("..","tmp"));
        
        rd.SaveGDAL(os.path.join("..","tmp","aspectInitialRes.tif"), aspect);
        rd.SaveGDAL(os.path.join("..","tmp","SlopeInitialRes.tif"), slope);
    else:
        print("AOI shp file" +shpfilepath + "is not inside of DEM"+drm_filepath+". Exiting.");
        input('Press Enter for exit...')
        exit;    

    #reopening SRTM products
    #read srtm products
    aspect_gdal_object = gdal.Open(os.path.join("..","tmp","aspectInitialRes.tif"))       #aspect
    aspect_band = aspect_gdal_object.GetRasterBand(1)
    aspect_band_array = aspect_band.ReadAsArray() 
    
    slope_gdal_object = gdal.Open(os.path.join("..","tmp","SlopeInitialRes.tif"))        #slope
    slope_band = slope_gdal_object.GetRasterBand(1)
    slope_band_array = slope_band.ReadAsArray() 
    
    #get PRODUCTS spatial resolution
    srtm_gt,srtm_xsize,srtm_ysize,srtm_ext,srtm_dpx,srtm_dpy=getGeotiffParams(aspect_gdal_object);
        
    
    #check if SRTM products inside of SHP AOI ad crop it
    if check_shp_inside_raster(srtm_ext,shp_extent):
        #do image crop
        aspect_cropped,ColMinInd,RowMinInd =crop_by_shp(shp_extent,srtm_ext,srtm_dpx,srtm_dpy,aspect_band_array)
        slope_cropped,ColMinInd,RowMinInd =crop_by_shp(shp_extent,srtm_ext,srtm_dpx,srtm_dpy,slope_band_array)
        
    else:
        print("SRTM is outside of the AOI, exiting...")
        exit();

was_corrected=False; #flag to check if resolution and scale were corrected to landsat8        
#STEP 1. CROP geotiffs one by one with AOI shape file
print("Step. 1. Starting geotiff crop operation...")        
for myfile in file_for_crop:
    
   
    #read geotiff
    gdal_object = gdal.Open(os.path.join('..',out_dir_name,myfile))
    band = gdal_object.GetRasterBand(1)
    band_array = band.ReadAsArray() 
    
    #get spatial resolution
    #do image crop
    gt,xsize,ysize,ext,dpx,dpy=getGeotiffParams(gdal_object)

    #check shp posiiton inside of tiff
    if check_shp_inside_raster(ext,shp_extent):
        #do image crop
        sampleImage,ColMinInd,RowMinInd =crop_by_shp(shp_extent,ext,dpx,dpy,band_array)
        
    else:
        print("AOI shp file" +shpfilepath + "is not inside of tiff"+myfile+". Exiting.");
        input('Press Enter for exit...')
        exit;
    
    #PARSING OF MTL FILE 
    #atmospheric correction 
    print('Aster files have been corrected radiometrically...');
        
                
    #topocorrection
    if is_topocorrection==True: #topocorrection flag    
        if was_corrected==False:
            
            #коррекция aspect по aster
            #adjust srtm resolution to aster
            #[hlc,wlc]=np.shape(sampleImage);
            [hlc,wlc]=np.shape(sampleImage);
            aspect_band_cropped=resize(aspect_cropped,(hlc,wlc),preserve_range=True,mode="wrap") #it works with scikit-image resize
    
            #коррекция slope по aster
            slope_band_cropped=resize(slope_cropped,(hlc,wlc),preserve_range=True,mode="wrap") #it works with scikit-image resize
            
            
            Cos_i=np.cos(np.deg2rad(slope_band_cropped))*np.cos(np.deg2rad(SolarZenith))+\
            np.sin(np.deg2rad(slope_band_cropped))*np.sin(np.deg2rad(SolarZenith))*\
            np.cos(np.deg2rad(SunAzimuth-aspect_band_cropped));
            
            #ЭТОТ РАСЧЕТ КОС I РАССМАТРИВАЕТ ВСЕ СКЛОНЫ КАК ОСВЕЩЕННЫЕ ПОД ПРЯМЫМ УГЛОМ!                                                                            
            #Cos_i=np.cos(np.deg2rad(SolarZenith-slope_band_cropped));
                        
            #Do SCS+C correction anyway
           
            #(b,a)=np.polyfit(Cos_i.ravel(),sampleImage.ravel(),1);
            (b,a)=np.polyfit(Cos_i.ravel(),sampleImage.ravel(),1);
            C=a/b;                                                                   
            was_corrected=True; #switch the flag to true                                                                                        
        
        print("Performing topographic correction.. Please, wait..")
        #Sun-Canopy-Sensor Correction (SCS)+C
        
        sampleImage_correct=np.float64(sampleImage*\
                ((np.cos(np.deg2rad(SolarZenith))*np.cos(np.deg2rad(slope_band_cropped))+C)\
                 /(C+Cos_i)));
        
        
        #band_array=np.float32(L_lambda); 
        pic_show(sampleImage,"aster band initial ");
        hist_show(sampleImage);
        pic_show(sampleImage_correct,"aster band  corrected to topography");
        hist_show(sampleImage_correct);
        #band_array_out=copy.copy(sampleImage);   #reflectance                    
    else: #no topocorrection
        print("No topocorrection was selected")
        sampleImage_correct=sampleImage.copy();
        
        
    #drop image to the disk
    print("drop image to the disk")
    outfilename=os.path.join(dir_cropped_path,"crop_"+myfile.lower());
    if not os.path.isdir(dir_cropped_path):
        os.makedirs(dir_cropped_path) #create output directory if none
    try:
        saveGeoTiff(sampleImage_correct,outfilename,gdal_object,ColMinInd,RowMinInd) #save topocorrected Landsat crop
    except:
        print("Can not write on a disk... and/or error(s) in saveGeoTiff function")
          
        
#STEP 2. COMPUTE pseudocolor RGB stacks and satellite indexes
"""
Autodetect BANDs for default names, if the user has not specified names (for now, default names),
skip index or RGB stack if we don't find BAND NUMBER
"""
print("Step. 2. Getting names of the cropped files...")        
#getting names of the cropped files, aquire band names
file_for_processing=[];
try:
    for file in os.listdir(dir_cropped_path): #We collect files from the folder with cropped images
        file=file.lower();
        if file.endswith("."+fileext.lower()):
            file_for_processing.append(file);
            print(file+" was added to crop queue.");
except(FileNotFoundError):
        print("Input image folder doesn\'t exist...");

bands={};  #dictionary storing band names 
for myfile in file_for_processing:
    for N in range(1,15):
        #populating bands dictionary
        if band_number_inname.replace('%n%',str(N),1) in myfile:
            try:
                if '3n' in myfile:
                    suffix='n';
                elif '3b' in myfile:
                    suffix='b';
                else:
                    suffix='';
                gdal_object = gdal.Open(os.path.join(dir_cropped_path,myfile)) #as new gdal_object was created, no more ColMinInd,RowMinInd
                bands['band'+str(N)+suffix]=gdal_object.GetRasterBand(1).ReadAsArray() ;
            except:
                print("Error! Can not read cropped bands!")
print("Bands dictionary output:")
print(bands) 

#create RGB stacks:
#truecolor
try:
    truecolorRGB=image_stack(bands['band3n'],bands['band2'],bands['band1'],do_norm8=1,do_show=1)  
except:
    truecolorRGB=image_stack(bands['band3b'],bands['band2'],bands['band1'],do_norm8=1,do_show=1)  
#Комбинация 8-3-2. Изображение близкое к естественным цветам, позволяет анализировать состояние атмосферы и дым. Здоровая растительность выглядит ярко зеленой, ярко розовые участки детектируют открытую почву, коричневые и оранжевые тона характерны для разреженной растительности.
try:
    b832RGB=image_stack(bands['band8'],bands['band3b'],bands['band2'],do_norm8=1,do_show=1)
except:
    b832RGB=image_stack(bands['band8'],bands['band3n'],bands['band2'],do_norm8=1,do_show=1)
#Комбинация RGB: 4/6, 4/7, 3/1 (Fig. 7) shows good results for lithologic discrimination  I. Di Tommaso, N. Rubinstein / Ore Geology Reviews 32 (2007) 275–290
try:
    b46_47_31RGB=image_stack(bands['band4']/bands['band6'],bands['band4']/bands['band7'],bands['band3n']/bands['band1'],do_norm8=1,do_show=1) 
except:
    b46_47_31RGB=image_stack(bands['band4']/bands['band6'],bands['band4']/bands['band7'],bands['band3b']/bands['band1'],do_norm8=1,do_show=1) 
#color composition RGB: 461 shows the Infiernillo alteration halo enhanced in two different concentric color zones, with the external magenta halo due to band 6 (Al–OH) absorption and the central yellow halo due to band 1 (Fe-oxides) absorption. I. Di Tommaso, N. Rubinstein / Ore Geology Reviews 32 (2007) 275–290
b461RGB=image_stack(bands['band4'],bands['band6'],bands['band1'],do_norm8=1,do_show=1)

#RGB: 4/5, 4/6, 4/7 (Fig. 6), the area in white shows a response of band 5 and band 6 (Al–OH) and band 7 (Fe–OH). I. Di Tommaso, N. Rubinstein / Ore Geology Reviews 32 (2007) 275–290
b45_46_47RGB=image_stack(bands['band4']/bands['band5'],bands['band4']/bands['band6'],bands['band4']/bands['band7'],do_norm8=1,do_show=1)    

#Комбинация 3-2-1. Растительность - красная 
try:
    b321RGB=image_stack(bands['band3n'],bands['band2'],bands['band1'],do_norm8=1,do_show=1)
except:
    b321RGB=image_stack(bands['band3b'],bands['band2'],bands['band1'],do_norm8=1,do_show=1)

#create indexes
try:
    NDVIC=((bands['band3n']-bands['band2'])/(bands['band3n']+bands['band2']))*\
        (1-(bands['band6']-bands['band6'].min())/(bands['band6'].max()-bands['band6'].min())); 
    #NDVI=(bands['band3n']-bands['band2'])/(bands['band3n']+bands['band2']) #NDVI
except:
    NDVIC=((bands['band3b']-bands['band2'])/(bands['band3b']+bands['band2']))*\
        (1-(bands['band6']-bands['band6'].min())/(bands['band6'].max()-bands['band6'].min())); 
    #NDVI=(bands['band3b']-bands['band2'])/(bands['band3b']+bands['band2']) #NDVI

#GENERAL OUTPUT
#print("Prepare to show PCA images")

#later incorporate path into functions
if not os.path.isdir(dir_products_path):
            os.makedirs(dir_products_path) #create output products directory if none
            

#COMPUTE Landsat and PCA stat for the CROSTA METHOD

stat_bands_save=os.path.join(dir_products_path,"bands_stat.xls");
cor_bands_save=os.path.join(dir_products_path,"bands_cor_stat.xls");
#cov_bands_pca_save=os.path.join(dir_products_path,"bands_pca_cov_stat.xls");

print("Saving band stat to {}".format(stat_bands_save));
save_landsat_bands_stat(bands,stat_bands_save);

print("Saving bands mutual correlation to {}".format(cor_bands_save));
save_landsat_mutual_cor(bands,cor_bands_save);

#save RGB's and index to the disk
print("Saving products on a disk")
if not os.path.isdir(dir_products_path):
    os.makedirs(dir_products_path) #create output directory if none
try:
    print("Saving RGBs...")
    ColMinInd=0; RowMinInd=0; #because we work on already cropped pictures
    saveGeoTiff(truecolorRGB,os.path.join(dir_products_path,"truecolorRGB.tif"),gdal_object,ColMinInd,RowMinInd);     
    saveGeoTiff(b832RGB,os.path.join(dir_products_path,"b832RGB"+".tif"),gdal_object,ColMinInd,RowMinInd);
    saveGeoTiff(b321RGB,os.path.join(dir_products_path,"b321RGB"+".tif"),gdal_object,ColMinInd,RowMinInd);
    saveGeoTiff(b46_47_31RGB,os.path.join(dir_products_path,"b46_47_31RGB"+".tif"),gdal_object,ColMinInd,RowMinInd);
    saveGeoTiff(b45_46_47RGB,os.path.join(dir_products_path,"b45_46_47RGB"+".tif"),gdal_object,ColMinInd,RowMinInd);
     #Aydal pseudocolor:
    saveGeoTiff(b461RGB,os.path.join(dir_products_path,"b461RGB"+".tif"),gdal_object,ColMinInd,RowMinInd);
    print("Saving Indexes...")
    saveGeoTiff(NDVIC,os.path.join(dir_products_path,"NDVIC"+".tif"),gdal_object,ColMinInd,RowMinInd);
    
    print("Products data were saved.")
except:
    print("Can not write PRODUCTS on a disk... and/or error(s) in saveGeoTiff function")

print("Operations were finished. It took {} sec".format(time.time()-time_start))