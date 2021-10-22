# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 15:33:10 2021

@author: Administrator
"""
import geopandas as gpd
from osgeo import ogr
import numpy as np
import pandas as pd
from rasterstats import zonal_stats
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error, r2_score
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.patches import PathPatch
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import gridspec

import os
# os.environ['PROJ_LIB'] = 'C:/ProgramData/Anaconda3/envs/geo_env/Library/share/proj'
os.environ['PROJ_LIB'] = 'C:/Users/HH/anaconda3/envs/geo_env/Library/share/proj' #win server512-2

import warnings
warnings.filterwarnings("ignore")
# 注册所有的驱动
ogr.RegisterAll()

polygons='县域产量2001-2015_研究区裁剪_filter_手动再检查_加均值_补空值.shp'
HI_data=pd.read_excel('各省HI.xls')
HI_data.index=HI_data.省名
HI_dct=HI_data.to_dict()['HI']

data = gpd.read_file(polygons,encoding ='UTF-8')
# data.plot()
# plt.show()#简单展示


# In[]
# year=2001
for year in range(2007,2016):
    print(year)
    data_new=data[['PAC','NAME','Prefecture','ProvinceNa','geometry','%d_yield'%year]]
    data_new['PAC']=data_new['PAC'].astype(int)
    
    # #即df.insert(添加列位置索引序号，添加列名，数值，是否允许列名重复)
    # data_new.insert(6,'err',0,allow_duplicates=False)
    
    #读取shp特定字段
    ds=ogr.Open(polygons, update=True) #打开矢量
    oLayer = ds.GetLayerByIndex(0)#获得图层，一般就是第一个
    nf=oLayer.GetFeatureCount(0)#获取Feature个数
    pac0=[0]*nf
    
    
    
    yields0=[0]*nf#县域统计
    for i in range(nf):
        oFeature = oLayer.GetNextFeature()
        pac0[i]=int(oFeature['PAC'])
        yields0[i]=oFeature['%d_yield'%year]
    # del ds
    #分区统计
    yield_raster='D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/空跑/open_loop_yield_%d_V2-3.tif'%year
    agb_raster='D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/空跑/open_loop_AGB_%d_V2-3.tif'%year
       
    yield_result=zonal_stats(polygons, yield_raster, stats="mean")
    yield_mean0=np.array([x['mean'] for x in yield_result])

    agb_result=zonal_stats(polygons, agb_raster, stats="mean")
    agb_mean0=np.array([x['mean'] for x in agb_result])
    
    
    #仅保留未参与标定的县
    yield_mean,agb_mean,pac,yields=[],[],[],[]
    used_pac=np.load('used_pac.npy')
    ip=-1
    for p in pac0:
        ip+=1
        if p in used_pac:
            continue
        yield_mean.append(yield_mean0[ip])
        agb_mean.append(agb_mean0[ip])
        yields.append(yields0[ip])
        pac.append(p)
        
    yield_zonal=dict(zip(pac,yield_mean))
    agb_zonal=dict(zip(pac,agb_mean))         
            
    temp=[]
    for p in data_new['PAC']:
        if p in used_pac:
            temp.append(None)
        else:
            temp.append(yield_zonal[p])
    data_new['sim']=temp#[yield_zonal[p] for p in data_new['PAC']]
    # data_new['yield_err']=(data_new['sim']-data_new['%d_yield'%year]*0.875)/(data_new['%d_yield'%year]*0.875)
    data_new['yield_err']=data_new['sim']-data_new['%d_yield'%year]*0.875
    agbs,agb_err=[],[]
    # yield_err=[]
    for n in range(nf):
        p=data_new['PAC'][n]
       
        proName=data_new['ProvinceNa'][n]
        if proName in HI_data.省名 and p not in used_pac:
            try:
                HI=HI_dct[proName]
                obs_agb=data_new['%d_yield'%year][n]*0.875/HI
                sim_agb=agb_zonal[p]
                if obs_agb:
                    # agb_err.append( (sim_agb-obs_agb)/obs_agb)
                    agb_err.append(sim_agb-obs_agb)
                else:
                    agb_err.append(np.NaN)
                agbs.append(obs_agb)
            except:
                agbs.append(np.NaN)
                agb_err.append(np.NaN)
        else:
            # agbs.append(np.NaN)
            agb_err.append(np.NaN)
    
    data_new['AGB_err']= agb_err   
 
    
    data_new.to_file('D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/空跑/%d_err_V2-3.shp'%year,encoding ='UTF-8')
    sim_data,obs_data={'yield':[],'AGB':[]},{'yield':[],'AGB':[]}
    for i in range(len(yield_mean)):
        '''异常值处理筛选'''
        if yields[i]>500 and yield_mean[i]!=None and yields[i]<8000 and \
            yield_mean[i]<8000 and yield_mean[i]>500 and agbs[i]==agbs[i]:
            sim_data['yield'].append(yield_mean[i])
            obs_data['yield'].append(yields[i]*0.875)
            sim_data['AGB'].append(agb_mean[i])
            obs_data['AGB'].append(agbs[i])        
    
    for item in ['yield','AGB']:   
        sim,obs=sim_data[item],obs_data[item]
        sim=pd.Series(sim)
        obs=pd.Series(obs)
        # fig, ax = plt.subplots(1, 2, figsize=(15,5))
        fig = plt.figure(figsize=(15,5.2))
        spec = gridspec.GridSpec(ncols=3, nrows=1,
                             width_ratios=[3,0.12, 2])#中间一列仅为了留出空隙
        ax0 = fig.add_subplot(spec[0])
        ax1 = fig.add_subplot(spec[2])
        ax=[ax0,ax1]
        m=Basemap(rsphere=(6378137.00,6356752.31414),llcrnrlon=79, llcrnrlat=14.5, \
              urcrnrlon=146.5, urcrnrlat=53, projection='aea', \
              lat_1=25, lat_2=47, lon_0=105, ax = ax[0])
      
        m.readshapefile('D:/hh/矢量边界/中国9段线边界','China',drawbounds=True,zorder=1,linewidth=0.2)
        m.readshapefile('D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/空跑/%d_err_V2-3'%year,'%s_err'%item,drawbounds=False,zorder=1,linewidth=0.2)    
        cmap=plt.cm.jet#.reversed()
        # colors = ['red','orangered',"darkorange", "gold", "greenyellow", "#2A8000",'forestgreen'][::-1]
        # nodes = [0.0, 0.15, 0.3, 0.45, 0.6, 0.8, 1.0]
        # cmap = LinearSegmentedColormap.from_list("mycmap", list(zip(nodes, colors)))
        
        # norm=plt.Normalize(-0.5,0.5)
        kv=2000 if item=='yield' else 4000
        norm=plt.Normalize(-kv,kv)
        
        parallels = np.arange(20.,90,15.)
        m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10, zorder=1, linewidth=0.2) # 绘制纬线
        
        meridians = np.arange(70.,140.,15.)
        m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10, zorder=1, linewidth=0.2) # 绘制经线
        
        patches   = []
        if item=='yield':
            zips=zip(m.yield_err_info, m.yield_err)
        else:
            zips=zip(m.AGB_err_info, m.AGB_err)
        for info, shape in zips:
          try:
            # if info['%s_err'%item]>=-0.5 and info['%s_err'%item]<=0.5:
            #   color=cmap(norm(info['%s_err'%item]))
            # elif info['%s_err'%item]<-0.5:
            #   color=cmap(norm(-0.5))
            # else:
            #   color=cmap(norm(0.5)) 
            
            if info['%s_err'%item]>=-kv and info['%s_err'%item]<=kv:
              color=cmap(norm(info['%s_err'%item]))
            elif info['%s_err'%item]<-kv:
              color=cmap(norm(-kv))
            else:
              color=cmap(norm(kv)) 
            patches.append( Polygon(np.array(shape), True, color=color) )
          except:
            pass
        
        pc = PatchCollection(patches, match_original=True, edgecolor=None, zorder=2)
        ax[0].add_collection(pc)
        
        ######局部小地图：九段线区域########
        axins = zoomed_inset_axes(ax[0], 0.53, loc = 4)
        axins.set_xlim(107, 123.5)
        axins.set_ylim(2, 26)
        
        map2 = Basemap(rsphere=(6378137.00,6356752.31414),llcrnrlon = 106.1, llcrnrlat = 3, urcrnrlon = 122.4, \
                         urcrnrlat = 24,projection='aea', lat_1=25, lat_2=47, lon_0=105, ax = axins)                      
        shpfile = 'D:/hh/矢量边界/中国9段线边界'             
        map2.readshapefile(shpfile, 'China',linewidth=0.2)              
        mark_inset(ax[0], axins, loc1=2, loc2=4, fc = "none", ec = "none")
        #####################################    
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        # sm.set_array(colvals)
        cb=fig.colorbar(sm, ax=ax[0],fraction=0.032,extend='both')
        font={'size':14}
        # cb.set_label('Relative error',fontdict=font) #设置colorbar的标签字体及其大小
        cb.set_label('Error',fontdict=font) #设置colorbar的标签字体及其大小
        plt.text(0.5,0.92,'%d-%s'%(year,item), fontsize=15, transform=ax[0].transAxes ,horizontalalignment='center')
        # plt.show() 
        scale=2 if item=='AGB' else 1 
        # ax[1].set_title("yield (with zero water)" )
        ax[1].set_xlim(xmin=100,xmax=8000*scale)
        ax[1].set_ylim(ymin=100,ymax=8000*scale)
        ax[1].set_xlabel("county statistics",fontdict=font)
        ax[1].set_ylabel("zonal simulation",fontdict=font)
        for i in range(len(obs)):
            r_err=(sim[i]-obs[i])#/obs[i]
            if r_err>=-kv and r_err<=kv:
                color=cmap(norm(r_err))
            elif r_err<-kv:
              color=cmap(norm(-kv))
            else:
                color=cmap(norm(kv))            
            ax[1].scatter(obs[i:i+1],sim[i:i+1],color=color, marker='v', edgecolors=color)
        ax[1].plot([0,20000],[0,20000])
        rmse=np.sqrt(mean_squared_error(obs,sim))
        r2=r2_score(obs,sim)
        r=obs.corr(sim)
        plt.text(0.2,0.85,'$\mathrm{R}^{2}$=%.2f'%r2,horizontalalignment='center',verticalalignment='center',transform=ax[1].transAxes)
        plt.text(0.2,0.8,'$\mathrm{RMSE}$=%d'%rmse,horizontalalignment='center',verticalalignment='center',transform=ax[1].transAxes)
        fig.savefig('D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/空跑/validation/%d_%s_V2-3.png'%(year,item),dpi=600,bbox_inches='tight')