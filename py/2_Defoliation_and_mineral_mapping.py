#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 15:57:43 2020

Defoliation script with mineral principal components computing.

Считается на основе спектральных каналов, прошедших топографическую коррекцию.
Вычисляются Direct Principal Components, на основе которых рассчитывается модель
оконтуривания месторождений
 
Вычисление компонентов и их нагрузок отсюда https://scentellegher.github.io/machine-learning/2020/01/27/pca-loadings-sklearn.html

"""

import gdal,ogr #OpenGIS Simple Features Reference Implementation
import numpy as np
from mygdal_functions0_9 import *
import matplotlib.pyplot as plt
import os
import pandas as pd
import time;
from sklearn.preprocessing import scale
from sklearn import decomposition
import copy;
from sklearn.cluster import KMeans
from matplotlib.colors import LinearSegmentedColormap
from sklearn.model_selection import train_test_split

from mygdal_functions0_9 import *
from configuration import *

#1 Settings

#files for processing, input and output directory
start_time=time.time();

cov_ratios_dpca_save_iron=os.path.join(dir_products_path,"ratios_dpca_cov_stat_iron.xls");
cov_ratios_dpca_save_clay=os.path.join(dir_products_path,"ratios_dpca_cov_stat_clay.xls");

loadings_dpca_save_iron=os.path.join(dir_products_path,"loadings_stat_iron.xls");
loadings_dpca_save_clay=os.path.join(dir_products_path,"ratios_dpca_cov_stat_clay.xls");
loadings_filemask='loadings_DPCA_{}.xls';
loadings_filemask2='loadings_DPCA_2_{}.xls';
variance_filemask='variance_DPCA_{}.xls';

band_number_inname='_b%N%.' #%N% - for band number e.g. LC81050292016143LGN00_B6.TIF NOT A CASE SENSITIVE
band_number_inname=band_number_inname.lower();
ColMinInd=0; RowMinInd=0; #because we work on already cropped pictures

#list of bands to be used in computations
bands=[2,3,4,5,6,7];

DO_KMEANS=False; #if this flag set to FALSE image will be neglected (*-1), if TRUE
                #k-means classififed
DO_NORM_NEG_NORM=True; #нормализует данные 0 - 1, шкалирует 1 - 0 в случае 
                        #отрицательных нагрузок loadings 
NDVI_CLASS_MASK=-1; #ndvi class to be excluded from computations -1 if unused


#bands to be used (Carranza, Hale, 2002)
#Alteration mineral to map	Band ratio images input to DPC (Landsat 8 OLI)

#Quartz	3/4, 7/1
#muscovite	 3/4, 6/7
#kaolinite	5/1, 7/4
#Chlorite	3/4, 7/5
#hematite 	5/4, 7/1
#limonite 3/4, 4/2


#A stabilized vegetation index (StVI)  (Ninomiya, Y. and Fu., B., 2005)

#2 data processing
files_for_processing=[];

try:
    for file in os.listdir(dir_cropped_path):         #exclude 3band tif files
        #file=file.lower();
        if file.lower().endswith("."+fileext.lower())  \
        and not file.lower().startswith("b") and not file.lower().startswith("true")\
        and not file.lower().startswith("cum"): 
            
            files_for_processing.append(file);
            print(file+" was added to data collecting queue.");
except(FileNotFoundError):
        print("Input image folder doesn\'t exist...");

#создание словаря каналов bands и загрузка туда файлов
bands={};  #dictionary storing band names 
for myfile in files_for_processing:
    for N in range(1,15):
        #populating bands dictionary
        #if band_number_inname.replace('%n%',str(N),1) in myfile:
         if str(N) in myfile.split('_')[-1]:   
            try:
                gdal_object = gdal.Open(os.path.join(dir_cropped_path,myfile)) #as new gdal_object was created, no more ColMinInd,RowMinInd
                bands['band'+str(N)]=gdal_object.GetRasterBand(1).ReadAsArray() ;
            except:
                print("Error! Can not read cropped bands!")


#Compute J.Carranza band rations and their PCs
#mineral band ratios were updated according to Prof. Carranza letter of 2021/02/05

m,n=bands['band4'].shape;

#Vegetation Normalized Difference Vegetation Index C
#https://www.indexdatabase.de/db/i-single.php?id=377 
NDVIC=((bands['band3']-bands['band2'])/(bands['band3']+bands['band2']))*\
    (1-(bands['band6']-bands['band6'].min())/(bands['band6'].max()-bands['band6'].min())); 


#Enhanced Vegetation Index 2 https://www.indexdatabase.de/db/s-single.php?id=9
EVI2=2.4*(bands['band3']-bands['band2'])/(bands['band3']+bands['band2']+1);

#Vegetation index... pickup other indeces from NVDI vs MIneral HERE ratios https://www.indexdatabase.de/db/s-single.php?id=9
#NDVIC=

#muscovite https://www.indexdatabase.de/db/i-single.php?id=49
#RMSC=bands['band3']/bands['band8']; #corrected to make vegs index >1
#RMSC=bands['band6']/bands['band5']; #- wrong loadings
#RMSC=bands['band6']/bands['band7'];
RMSC=bands['band4']/bands['band5']; 
#RMSC=bands['band5']/bands['band9'];
VI0=bands['band3']/bands['band4'];

#kaolinite 
#RKL=bands['band5']/bands['band7']; #corrected to make vegs index >1 
RKL=bands['band7']/bands['band9'];


#Carbonate/Chlorite/Epidote https://www.indexdatabase.de/db/s-single.php?id=9 
#RCCE=(bands['band7']+bands['band9'])/bands['band8'];      
RE=bands['band5']/bands['band8'];
#RE=bands['band3']/bands['band5']; #in fact that is chlorite
#RE=bands['band5']/bands['band2']; 
#RE=bands['band8']/bands['band9']; 
#RE=bands['band4']/bands['band8']; 

#quartz
RQ_Ref=bands['band11']**2/(bands['band10']*bands['band12']); #minerals Ration Quartz Mineral          
#RQ=bands['band1']/bands['band9']; # - положительная зависимость для PDC, слабая  модель, до 0.702 
#RQ=bands['band4']/bands['band1']; #- слабая отрицательная связь
RQ=bands['band1']/bands['band8']; #нейтральная или слабая положительная связь
#RQ=bands['band11']/bands['band12']; #- один знак с вегетацией в компоненте
#RQ=bands['band10']/bands['band11'];
#RQ=bands['band13']/bands['band12'];
#RQ=bands['band10']/bands['band13']; #дает обратную корреляцию в модели
#RQ=bands['band13']/bands['band12']; # - это отношение и bands['band14']/bands['band12'] для раст дает обратную корреляцию в модели
#RQ=RQ_Ref.copy();
#VI1=bands['band14']/bands['band12'];
VI2=bands['band3']/bands['band2']; #- слабая отрицательная связь
 
#hematite, Fe2+ https://www.indexdatabase.de/db/i-single.php?id=18 
#RHEM=(bands['band5']/bands['band3'])+(bands['band1']/bands['band2']); #minerals          
#RHEM=bands['band2']/bands['band9'];
RHEM=bands['band3']/bands['band9'];

#limonite, Fe3+ https://www.indexdatabase.de/db/i-single.php?id=19 
#RLIM=bands['band3']/bands['band7']; #changed to make vegetation >1 -0.000079
#RLIM=bands['band4']/bands['band1']
#RLIM=bands['band3']/bands['band6'] # - small value for important mineral
#RLIM=bands['band1']/bands['band5']
#RLIM=bands['band5']/bands['band2']
#RLIM=bands['band2']/bands['band7']
RLIM=bands['band5']/bands['band9'] # - 3%, ok? 
#RLIM=bands['band4']/bands['band9']


#nvdi classes using K-means
ndvi_classes=kmeans_sort(NDVIC,do_reverse=0,n_clusters=total_number_ndvi_classes,elbow=1);
             
#show k-maens classes ndvi


colors = [(1, 0, 0),(0, 1, 0),(0, 0.8, 0.4),(0, 0.7, 0.1),(0, 0.5, 0.2),(0, 0, 1)]  # R -> G -> B
n_bins = total_number_ndvi_classes  # Discretizes the interpolation into bins
cmap_name = 'my_list'
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins);

#tell the colorbar to tick at integers
im=plt.imshow(ndvi_classes,cmap=cm);
plt.colorbar(im,ticks=np.arange(np.min(ndvi_classes),np.max(ndvi_classes)+1))

plt.savefig('ndvi_classes_area.png',dpi=300);
plt.savefig('ndvi_classes_area.svg',dpi=300);
plt.show();

#save ndvi classes
saveGeoTiff(ndvi_classes,ndvi_file_path,gdal_object,ColMinInd,RowMinInd);

"""
NDVI classes (empirical description)
3 - water and ice
0 - open soil
1 - rare vegetation
2 - dense vegetation
"""


#assessment of channel standart deviation, selecting channels with the higher deviation
minerals={}
vegetation={}


minerals.update({'quartz':RQ});
vegetation.update({'quartz':VI2});

minerals.update({'muscovite':RMSC});
vegetation.update({'muscovite':VI0});

minerals.update({'kaolinite':RKL});
vegetation.update({'kaolinite':VI0});

minerals.update({'chlorite':RE});
vegetation.update({'chlorite':VI0});

minerals.update({'hematite':RHEM});
vegetation.update({'hematite':VI0});

minerals.update({'limonite':RLIM});
vegetation.update({'limonite':VI0});

#compute leaves/minerals DPC's
#compute DPCA for minerals vs vegetation 5/4
#1 flatten image matrixes
print("Started to compute DPCA for the iron oxides and clay minerals...")
print("Flatten image matrix...")

mineral_vegetation_flat={};  #compute flattened matrix for all minerals
for key in [*minerals]:
    mineral_vegetation_flat.update({key:mat4pca((minerals[key],vegetation[key]))});     

#2 compute PCA for minerals
pca_models_minerals={};
pca_dpca_minerals={};

"""
mineral_vegetation_flat[key] - минерал и его вегетация в виде таблицы с двумя столбцами

NDVI_CLASS_MASK=0

mask_index=(ndvi_classes.flatten()!=NDVI_CLASS_MASK); #булевы индексы "не-вода"
table4pca=mineral_vegetation_flat[key][mask_index] #выбираем в таблицу для PCA только классы суши

"""

#selecting/masking non-water NDVI pixels 
for key in [*minerals]:
    pca_models_minerals.update({key:decomposition.PCA(n_components=2)});
    #updates of 2021/02/26
    mask_index=(ndvi_classes.flatten()!=NDVI_CLASS_MASK); #булевы индексы "не-вода"
    blank_table_pca=np.zeros(mineral_vegetation_flat[key].shape);
    table4pca=mineral_vegetation_flat[key][mask_index]
    found_pca=pca_models_minerals[key].fit_transform(table4pca); #table with pca inside    
    blank_table_pca[mask_index]=found_pca; #place computed DPCA into blank table, mask stays zero
    pca_dpca_minerals.update({key:blank_table_pca});
    

#3 loadings for minerals/vegetation (eigenvectors/eigenvalues)
                                            
                                           
loadings_minerals={};
for key in [*minerals]:
    #split 30% of data set for computation loadings
    #mineral_veg_train, mineral_veg_test = train_test_split(mineral_vegetation_flat[key],\
    #                                test_size=0.66, random_state=42);
    loadings = get_comp_loadings2(pca_models_minerals[key],\
        features=[key,'vegetation'],columns=['DPC1', 'DPC2']);
    print(loadings);                              
    loadings_minerals.update({key:loadings}); #    

#В ЭТОЙ ВЕРСИИ ИСПРАВЛЕНО ЗНАЧЕНИЕ МНОЖИТЕЛЕЙ ДЛЯ negation 
#4 DPCA for minerals
DPCA_minerals={};
for key in [*minerals]:
    pca_image=get_pca_image(pca_dpca_minerals[key],m,n,n_components=2);
    pca_image_classes={};
    for pckey in pca_image: #classify PCA images
        if loadings_minerals[key][('DPC'+pckey)][key]<0:
            do_reverse=1; #должно быть обусловлено знаками DPCA ЕСЛИ НАГРУЗКА ОТРИЦАТЕЛЬНА, МНОЖИМ на -1
            multiplier=-1;
        else:
            do_reverse=0;
            multiplier=1; #NO negation
        
        
        if DO_KMEANS==False:
            #no K-means classification!!!!!!!!!! 
            if DO_NORM_NEG_NORM==False:
                img_class=pca_image[pckey]*multiplier;
            else:
                #print('multiplier={}'.format(multiplier));
                if multiplier==1:
                    img_class=(pca_image[pckey]-np.min(pca_image[pckey]))/\
                        (np.max(pca_image[pckey])-np.min(pca_image[pckey]));
                    #print('img_class={}'.format(img_class));
                else: #multiplier==-1
                    img_class=(np.max(pca_image[pckey])-pca_image[pckey])/\
                        (np.max(pca_image[pckey])-np.min(pca_image[pckey]));
                    #print('img_class={}'.format(img_class));
        else:
            #uncomment line below if you wish k-means classification
            img_class=kmeans_sort(pca_image[pckey],do_reverse=do_reverse,n_clusters=10);
        
        #img_class=mean_shift_sort(pca_image[pckey],do_reverse=do_reverse);
        pca_image_classes.update({pckey:img_class});
        #pca_image_classes.update({pckey:get_image_classes(pca_image[pckey])});
    DPCA_minerals.update({key:pca_image_classes});    


#5 save loadings
for key in [*minerals]:
    loadings_minerals[key].to_excel(os.path.join(dir_products_path,\
                     loadings_filemask.format(key)),index=True);
    
                    
#6 save DPCs into geotiff and matplotlib image
for key in [*minerals]:
    for n in DPCA_minerals[key]: #save normalize values of DPCA matrix
        #img=normalize_matrix(DPCA_minerals[key][n]);
        #apply negation 
        img=DPCA_minerals[key][n];
        
        #if loadings_minerals[key][key]['DPC'+n]<0 and loadings_minerals[key]['vegetation']['DPC'+n]>0:
        #    print('negation for mineral %s, %s'%(key,'DPC'+n));
            #img=np.ones(img.shape)*np.max(img)-img;
        #    img=np.ones(img.shape)-img;
        
        #DPCA_minerals[key][n]=img; #replace minerals
        
        saveGeoTiff(img,os.path.join(dir_products_path,"DPC{}_{}_".format(n,key)+"_"+".tif"),gdal_object,ColMinInd,RowMinInd);        

saveGeoTiff(NDVIC,os.path.join(dir_products_path,"NDVIC"+"_"+".tif"),gdal_object,ColMinInd,RowMinInd); 
saveGeoTiff(EVI2,os.path.join(dir_products_path,"EVI2"+"_"+".tif"),gdal_object,ColMinInd,RowMinInd); 

plt.figure();
plt.subplot(231);
plt.imshow(DPCA_minerals['limonite']['1'],cmap='gray');
plt.title('\'limonite\' image');
plt.axis('off');
plt.subplot(232);
plt.imshow(DPCA_minerals['quartz']['1'],cmap='gray');
plt.title('\'quartz\' image');
plt.axis('off');
plt.subplot(233);
plt.imshow(DPCA_minerals['hematite']['1'],cmap='gray');
plt.title('\'hematite\' image');
plt.axis('off');
plt.subplot(234);
plt.imshow(DPCA_minerals['muscovite']['1'],cmap='gray');
plt.title('\'muscovite\' image');
plt.axis('off');
plt.subplot(235);
plt.imshow(DPCA_minerals['kaolinite']['1'],cmap='gray');
plt.title('\'kaolinite\' image');
plt.axis('off');
plt.subplot(236);
plt.imshow(DPCA_minerals['chlorite']['1'],cmap='gray');
plt.title('\'chlorite\' image');
plt.axis('off');
plt.savefig('mineral_maps_DPC1.png',dpi=300);
plt.savefig('mineral_maps_DPC1.svg',dpi=300);
plt.show();

###

plt.figure();
plt.subplot(231);
plt.imshow(DPCA_minerals['limonite']['2'],cmap='gray');
plt.title('\'limonite\' image');
plt.axis('off');
plt.subplot(232);
plt.imshow(DPCA_minerals['quartz']['2'],cmap='gray');
plt.title('\'quartz\' image');
plt.axis('off');
plt.subplot(233);
plt.imshow(DPCA_minerals['hematite']['2'],cmap='gray');
plt.title('\'hematite\' image');
plt.axis('off');
plt.subplot(234);
plt.imshow(DPCA_minerals['muscovite']['2'],cmap='gray');
plt.title('\'muscovite\' image');
plt.axis('off');
plt.subplot(235);
plt.imshow(DPCA_minerals['kaolinite']['2'],cmap='gray');
plt.title('\'kaolinite\' image');
plt.axis('off');
plt.subplot(236);
plt.imshow(DPCA_minerals['chlorite']['2'],cmap='gray');
plt.title('\'chlorite\' image');
plt.axis('off');
plt.savefig('mineral_maps_DPC2.png',dpi=300);
plt.savefig('mineral_maps_DPC2.svg',dpi=300);
plt.show();

plt.figure();
plt.imshow(RQ_Ref,cmap='gray');
plt.title('RQ as per reference');
plt.colorbar();
plt.show();


plt.figure();
plt.imshow(EVI2,cmap='gray');
plt.title('EVI2 as per reference');
plt.colorbar();
plt.show();

plt.figure();
plt.imshow(NDVIC,cmap='gray');
plt.title('NDVIC as per reference');
plt.colorbar();
plt.show();

plt.figure();
plt.imshow(bands['band2'],cmap='gray');
plt.title('BAND2');
plt.colorbar();
plt.show();

print('computations took {} s'.format(time.time()-start_time));
