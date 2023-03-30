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
from scipy.sparse import csr_matrix 
import pickle

#one class SVM
from sklearn import svm, metrics
#selecting training and testing set
from sklearn.model_selection import train_test_split
import random
from sklearn.inspection import plot_partial_dependence, partial_dependence

#RF
#importing Random Forest model
from sklearn.ensemble import RandomForestRegressor
from sklearn.datasets import make_regression

#importing classes for validation
from sklearn.metrics import r2_score,mean_squared_error
from sklearn.preprocessing import StandardScaler  #https://scikit-learn.org/stable/modules/neural_networks_supervised.html

from mygdal_functions0_9 import *
from configuration import *


#1 Исходные параметры
filename_model=os.path.join(filedir_model,model_maxent_name);


file_model_data_name='model_data_{}{}'.format(pathrowfolder_new,datefolder_new);
file_model_data_name_path=os.path.join(dir_products_path_new,file_model_data_name);
file_model_data_name_path_out=os.path.join(dir_products_path_new,(file_model_data_name));


#getting resolution from file    
for file in os.listdir(dir_cropped_path_new):         #exclude 3band tif files
        if file.lower().endswith("."+fileext.lower()):
            gdal_object = gdal.Open(os.path.join(dir_cropped_path_new,file)); #as new gdal_object was created, no more ColMinInd,RowMinInd
            break;                          

row,col=gdal_object.GetRasterBand(1).ReadAsArray().shape;  


#2 Чтение
#read saved model
with open(filename_model, 'rb') as f:
    clf = pickle.load(f)

#read saved data
with open(file_model_data_name_path, 'rb') as f:
    model_data_df = pickle.load(f)

#convert datatable to dictionary
column_names=[*model_data_df];


#select DPC1 only indexes
dpc1_columns=[];
for col_name in [*model_data_df]:
    if col_name.lower().startswith('dpc'):
        
        #select only specific minerals

        if  col_name == 'DPC2_goethite' or col_name == 'DPC2_quartz' \
            or col_name == 'DPC2_alunite' or col_name == 'DPC2_illite':
                dpc1_columns.append(col_name);

#arrange columns in the necessary order into model_np_array
dpc1_columns=sorted(dpc1_columns, reverse=False);

#
all_data_np_array=np.array(np.array(model_data_df['id']).reshape(-1,1));
all_data_np_array=np.hstack([all_data_np_array,np.array(model_data_df['row']).reshape(-1,1),\
                             np.array(model_data_df['col']).reshape(-1,1)]);

for colname in dpc1_columns:
    col_np=np.array(model_data_df[colname]);
    all_data_np_array=np.hstack([all_data_np_array,col_np.reshape(-1,1)]);

#
#APPLYING THE MODEL TO 100% STRINGS AND OUTPUT
mean_all= all_data_np_array[:,3:].mean(axis=0);
std_all = all_data_np_array[:,3:].std(axis=0);
all_data_std = (all_data_np_array[:,3:] - mean_all) / std_all;                    
pred = clf.decision_function(all_data_std) #prediction

#select 1/1000 of data for responce plots
all_data_std_train,all_data_std_test=train_test_split(all_data_std, test_size=0.01); #,random_state=42);  #


#normalize output
pred_norm=(pred-pred.min())/(pred.max()-pred.min());    
    
    

#save geotiff file of predicted values
#render predicted values into raster
predicted_raster=np.empty([row,col],dtype=np.float64);
predicted_raster[:]=np.nan;

for i in range(0,len(all_data_np_array[:,0]),1):
    r,c=all_data_np_array[i,1],all_data_np_array[i,2];
    predicted_raster[int(r),int(c)]=pred_norm[i];



#save geotiff with normalized predicted data
saveGeoTiff(predicted_raster,os.path.join(dir_products_path_new,\
    "predicted_oneClassSVM_NEW{}_".format('')+pathrowfolder_new+"_"+datefolder_new+".tif"),\
    gdal_object,0,0);
    

    #save dataframe
model_data_df=model_data_df.assign(Prob=pred);  
    
#save model data into file
with open(file_model_data_name_path_out, 'wb') as f:
    pickle.dump(model_data_df, f);
    




#create partial dependence plots
feature_names=dpc1_columns;
features = np.arange(0,len(feature_names),1)

fig, ax = plt.subplots(figsize=(12, 6))
#fig, ax = plt.subplots(figsize=(2, 2))
part_dep=plot_partial_dependence(clf, all_data_std_test,features,feature_names,n_cols=3, grid_resolution=50,ax=ax)
fig = plt.gcf()
plt.setp(ax.get_yticklabels(), fontsize=3);
plt.setp(ax.get_xticklabels(), fontsize=3);
fig.suptitle('Response curves for testing (Iturup) dataset')
fig.subplots_adjust(hspace=0.7)
#fig.text(0.06, 0.5, 'common ylabel', ha='center', va='center', rotation='vertical')
fig.savefig('Partial_dependence_testing.svg',dpi=300);
fig.savefig('Partial_dependence_testing.png',dpi=300);  
