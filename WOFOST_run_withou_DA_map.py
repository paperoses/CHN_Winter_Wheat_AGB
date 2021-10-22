# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 21:15:15 2021
空跑结果生成tif
@author: Administrator
"""
from osgeo import gdal
import glob
import numpy as np
import matplotlib.pyplot as plt

# for year in range(2014,2016,1):
    # print(year)
g = gdal.Open('D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/11省分区范围/tif_01dg_clip.tif')
geo_transform=g.GetGeoTransform()
rows=g.RasterYSize
cols=g.RasterXSize

# for stage in [[2001,2002,2003,2004,2005],[2006,2007,2008,2009,2010],[2011,2012,2013,2014,2015]][:]:
for stage in [[2001,2002,2003],[2004,2005,2006],[2007,2008,2009],[2010,2011,2012],[2013,2014,2015]][2:]:
# for stage in [[2001],[2002],[2003],[2004],[2005],[2006],[2007],[2008],[2009],[2010],[2011],[2012],[2013],[2014],[2015]][:]:    
    
    add=-1
    for item in ['yield','AGB']:
        add+=1
        fs=glob.glob( "D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/空跑/TEST%d-%d_V2-3*.csv"%(stage[0],stage[-1]) )
    
        # data = np.empty((0,13), int)
        data = np.empty((0,3+2*len(stage)), int)
        for f in fs:   
            try:                                
                x = np.loadtxt(f,delimiter=',',dtype=int)
            except:
                f_lines=open(f).readlines()
                x=[]
                for f_line in f_lines:
                    # if len(f_line)==13:
                    if len(f_line)==9:    
                        continue
                    else:
                        l=f_line.split(',')
                        l2=[int(a) for a in l]
                        x.append(l2)
            if len(x)!=0:
                data= np.vstack((data,x)) 
        for year in stage:
            print(year,item)
            image = np.zeros((cols, rows), int)
            for line in data:
                row = line[1]
                col = line[0]
                value=line[3+add+2*(year-stage[0])]
                '''异常值处理筛选'''
                if (value>1000 and value<8500 and item=='yield') or (value>2000 and value<18000 and item=='AGB') :
                    # image[col-1, row-1] = value
                    image[col, row] = value
            image=np.ma.masked_where(image==0, image)
            plt.imshow(image.transpose(), interpolation='nearest', vmin=1000, vmax=10000, cmap=plt.cm.jet)
            plt.show()
        
            driver = gdal.GetDriverByName ( "GTiff" ) # Get a handler to a driver
            dataset_y = driver.Create ( 'D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/空跑/open_loop_%s_%d_V2-3.tif'%(item,year),cols, rows,1, gdal.GDT_Int16, options=['COMPRESS=LZW'] )
            dataset_y.SetGeoTransform ( geo_transform)
            dataset_y.SetProjection ( g.GetProjectionRef() )
            dataset_y.GetRasterBand(1).SetNoDataValue(0)
            dataset_y.GetRasterBand(1).WriteArray(image.transpose())
            dataset_y.FlushCache()
            del dataset_y
        

del g