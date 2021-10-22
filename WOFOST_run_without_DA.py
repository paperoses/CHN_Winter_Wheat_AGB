# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 11:36:32 2021
新疆参数的差异在于TBASE叶片生长的最低温设为了-10
@author: Administrator
"""
# import gdal
import os
import time
import datetime
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
import geopandas as gpd
import scipy.stats as ss
import csv
from multiprocessing import Pool #导入进程池

from osgeo import gdal
import matplotlib.colors as mcolors
from pcse.util import Afgen
from pcse.models import Wofost71_PP,Wofost71_WLP_FD
from pcse.base import ParameterProvider
from pcse.fileinput import YAMLAgroManagementReader, CABOFileReader
from pcse.fileinput import YAMLCropDataProvider,CABOWeatherDataProvider
from pcse.util import WOFOST71SiteDataProvider
from netCDF4 import Dataset,date2index


# In[]
def year_num(d,lst):
    '''主要针对冬小麦出苗在前一年，lst第一项是出苗的doy'''
    dd=2 if d==lst[0] else 1
    return dd
    
class creatAgro():    
    def __init__(self, code, run_years, calendar_doys, outputfile,level, SM, irrigation):
        self.code = code
        if type(run_years)==int:
            run_years=[run_years]
        self.run_years=run_years
        self.calendar_doys=calendar_doys
        self.outputfile = outputfile 
        self.level=level
        self.SM=SM
        self.irrigation=irrigation
    
    def __call__(self):
#         year = self.year_begin
        n=0
#         while year <= self.year_end: 
        for year in self.run_years:
            em,flw,mt=[datetime.datetime(year-year_num(d,self.calendar_doys[n]),12,31)+datetime.timedelta(d) for d in self.calendar_doys[n]]
            file=os.path.join(self.outputfile,'amgt','%d_%d_V2-3_L%d.amgt'%(self.code,year,self.level))
            # if os.path.exists(file):
            #     continue
            f=open(file,'w+')
            f.write('Version: 1.0.0\n'
                    +'AgroManagement:\n- %s:\n'%em.strftime("%Y-%m-%d")
                    +'    CropCalendar:\n'
                    +'        crop_name: wheat\n'
                    +'        variety_name: Winter_wheat_105\n'
                    +'        crop_start_date: %s\n'%em.strftime("%Y-%m-%d")
                    +'        crop_start_type: emergence\n'
                    +'        crop_end_date: %s\n'%mt.strftime("%Y-%m-%d")
                    +'        crop_end_type: earliest\n'
                    +'        max_duration: 330\n'
                    +'    TimedEvents:\n'
                    +'    StateEvents:\n'
                    # +'    -   event_signal: irrigate\n'
                    # +'        event_state: DVS\n'
                    # +'        zero_condition: rising\n'
                    # +'        name:  Soil moisture driven irrigation scheduling\n'
                    # +'        comment: All irrigation amounts in cm\n'
                    # +'        events_table:\n'
                    # +'        - 0.3: {amount:  %.1f, efficiency: 0.8}\n'%self.irrigation1                     
                    # +'        - 1.2: {amount:  %.1f, efficiency: 0.8}\n'%self.irrigation2
                    +'    -   event_signal: irrigate\n'
                    +'        event_state: SM\n'
                    +'        zero_condition: falling\n'
                    +'        name:  Soil moisture driven irrigation scheduling\n'
                    +'        comment: All irrigation amounts in cm\n'
                    +'        events_table:\n'
                    +'        - %.4f: {amount:  %.3f, efficiency: 0.7}\n'%(self.SM,self.irrigation)                    
                    +'- %s:'%mt.strftime("%Y-%m-%d")
                    )
            f.close()
#             year+=1
            n+=1
        return file#os.path.join(self.outputfile,'amgt','%d_????_V2-3_L1.amgt'%(self.code))

class ModelRerunner4(object):
    """recaliobrate irragation"""
    
    parameters = ["SLATB", "SPAN", "AMAXTB", "TMPFTB", "TMNFTB", "CVO", 'FLTB', 'FSTB', 'FOTB']
    data_dir = './'
    # data_dir = 'H:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定'
    
    def __init__(self, code, params, wdp, years, tsum=None, DTSM=0, level='PP',level2=1,x=None,calendar_doys=None):
        self.code = code
        self.params = params
        self.wdp = wdp
        # self.agro = agro
        self.years = years
        self.tsum=tsum
        self.DTSM=DTSM
        self.level = level    
        self.level2=level2
        self.par_values=x
        self.calendar_doys=calendar_doys
        self.VB,self.VS=self.params['VERNBASE'],self.params['VERNSAT']
        
    def __call__(self, irragation, return_df=False):
        par_values=self.par_values
        out=[]
        new_pars=new_pars=[[],0,[],[],[],0,[],[],[]]
            
        
        ny=len(self.years)

        my_year=list(self.years)
        irr=np.array(irragation).reshape([2,ny]).T
        irr_dct=dict(zip(my_year,irr))
        
        for year in self.years:
            agromanagement_file =creatAgro(self.code, year, self.calendar_doys, self.data_dir, self.level2, irr_dct[year][0],irr_dct[year][1])()
            agro=YAMLAgroManagementReader(agromanagement_file)
            # Clear any existing overrides
            self.params.clear_override()
            # Set overrides for the new parameter values
            # for parname, value in zip(['VERNBASE','VERNSAT'], [self.VB,self.VS]):
            #     self.params.set_override(parname, value) 
            SLATB=[0.0, par_values[0], 0.5, par_values[0], 2.0, par_values[0]]
            AMAXTB=[0.0, par_values[2], 1.0, par_values[2], 1.3, par_values[2], 2.0, par_values[2]]
            TMPFTB=[0.0+par_values[4], 0.01, 10.0, 0.6, 15.0, 1.0, par_values[3], 1.0, 35.0, 0.0]
            TMNFTB=[0.0+par_values[4], 0.0, 3.0+par_values[4], 1.0]
            k_dvs = par_values[6]
            FLTB  = [0.00,    0.65,  
                     0.10,    0.65,  
                     0.25,    0.70,
                     0.50,    0.50,
                     0.646,   0.30,
                     0.95,    0.00,
                     1.00,    0.00,
                     2.00,    0.00]
            FLTBs=Afgen(FLTB)
            
            FSTB  = [0.00,    0.35,
                     0.10,    0.35, 
                     0.25,    0.30,
                     0.50,    0.50,
                     0.646,   0.70,
                     k_dvs,   1-FLTBs(k_dvs),
                     1.00,    0.00,
                     2.00,    0.00]
            FSTBs=Afgen(FSTB)
            
            FOTB  = [0.00,    0.00,
                     0.10,    0.00,  
                     0.25,    0.00,
                     0.50,    0.00,
                     0.646,   0.00,
                     k_dvs,   0.00,
                     0.95,    1-FSTBs(0.95),
                     1.00,    1.00,
                     2.00,    1.00]   
            # print(len(new_pars),len(par_values))
            new_pars[1],new_pars[5]=par_values[1],par_values[5]
            new_pars[0],new_pars[2],new_pars[3],new_pars[4],new_pars[6],new_pars[7],new_pars[8]= \
            SLATB,AMAXTB,TMPFTB,TMNFTB,FLTB,FSTB,FOTB
            for parname, value in zip(self.parameters, new_pars):
                self.params.set_override(parname, value)
            if self.tsum!=None:
                self.params.set_override('TSUM1', self.tsum[0])
                self.params.set_override('TSUM2', self.tsum[1])
            if self.DTSM!=0:
                self.params.set_override('DTSMTB', [self.DTSM,0,30.0,30.0,45.0,30.0])            
            # Run the model with given parameter values
            if self.level=='PP':
                wofost = Wofost71_PP(self.params, self.wdp, agro)
            elif self.level=='WLP_FD':
                wofost = Wofost71_WLP_FD(self.params, self.wdp, agro)
            else:
                msg = "Wrong level name!"
                raise RuntimeError(msg)
            wofost.run_till_terminate()
            df = pd.DataFrame(wofost.get_output())
            df.index = pd.to_datetime(df.day)
            if return_df:
                out.append(df) 
            else:
                out.append([df.TWSO.max(),df.TAGP.max(),df.LAI.max()])          
        return out
    
def mycallback(x):
    print(x[:2])
    csv_write.writerow(x)  



def sayHi(num):
     
    data_dir='D:\hh\@进行中的工作\ESSD_Wheat_Multi-vars\分区标定'   
    # data_dir='H:\hh\@进行中的工作\ESSD_Wheat_Multi-vars\分区标定' 
    
    x,y,lon,lat,df_v,par_lst,level2=num[0],num[1],num[2],num[3],num[4],num[5],num[6]
    sited,cropd,soild=par_lst
   
    code=df_v['区站号']
    if code in replace_code.keys():
        code=replace_code[code]
    
    p_file='cropdata/calibrated//%d-%d_%d_V2-3_L%d.txt'%(run_years[0],run_years[-1],code,level2)
    if os.path.exists(p_file):
        pars=np.loadtxt(p_file)
        print('Calibrated crop pars %d-%d_%d_V2-3_L%d.txt 已读取'%(run_years[0],run_years[-1],code,level2))
    else:
        print('no calibrated crop pars %d-%d_%d_V2-3_L%d.txt'%(run_years[0],run_years[-1],code,level2))
        pars=[829,699,0.0019,22,32,22,1,0.7,0.66]
    site="%3.1f_%2.1f"%(np.floor(lon*10)/10.,np.floor(lat*10)/10.)#这里要和气象生成一致，前面代码注意修改
    calendar_doys=[]
    for year in run_years:
        calendar_doys.append(list(map(int,[df_v['%d_Em'%year],df_v['%d_Hg'%year]+5,df_v['%d_Mt'%year]])))

    weatherfile = os.path.join(data_dir, 'weather', site)
    wdp = CABOWeatherDataProvider(weatherfile) 
    
    
    #update soil data accroding He Liang's data
    soild.update({'SMW': df_v['SMW'],
                  'SMFCF': df_v['SMFCF'],
                  'SM0': df_v['SM0'],
                  'CRAIRC': df_v['CRAIRC']})
    # Package everyting into a single parameter object
    params = ParameterProvider(cropdata=cropd, sitedata=sited, soildata=soild)
    
    # tsum1,tsum2,SLATB, SPAN, AMAXTB, TMPFTB, TMNFTB, CVO, shift_dvs=pars
    tsum1,tsum2,DTSM,SLATB, SPAN, AMAXTB, TMPFTB, TMNFTB, CVO, shift_dvs=pars
    '''modelrunner2'''
    # model=ModelRerunner2(code, params, wdp, run_years, tsum=[tsum1,tsum2], level='WLP_FD')  
    # dfs=model([SLATB, SPAN, AMAXTB, TMPFTB, TMNFTB, CVO, shift_dvs],return_df=True)
    '''modelrunner4'''
    model=ModelRerunner4(code, params, wdp, run_years, tsum=[tsum1,tsum2], DTSM=DTSM,\
                          level='WLP_FD',level2=level2,x=[SLATB, SPAN, AMAXTB, TMPFTB, TMNFTB, CVO, shift_dvs], \
                              calendar_doys=calendar_doys)  
    

    i_file='amgt/calibrated//%d-%d_%d_V2-3_L%d.txt'%(run_years[0],run_years[-1],code,level2)
    if os.path.exists(i_file):
        irragation=np.loadtxt(i_file)
        print('Calibrated irragation 已读取')
    else:
        print('no calibrated irragation')
        # irragation=[0.2,4]
    
    dfs=model(irragation,return_df=True)
    
    print('%d'%code)
    df_dct=dict(zip(run_years,dfs))
    out=[str(x),str(y),str(code)]
    for year in run_years:
        out.append(df_dct[year].TWSO.max()) 
        out.append(df_dct[year].TAGP.max()) 
    return out



replace_code = {56768: 56571,
                56651: 56571,
                56598: 57508,
                56778: 56571,
                57517: 57411,
                57625: 57508,
                57432: 57326,
                57512: 57411,
                52876: 52983,
                52868: 52983,
                53914: 57002,
                53817: 53930} 

# run_years=[2001,2002,2003]
# run_years=[2004,2005,2006]
run_years=[2007,2008,2009]  
run_years=[2010,2011,2012]
run_years=[2013,2014,2015]  

# run_years=[2001]
# run_years=[2002]
# run_years=[2003]
# run_years=[2004]
# run_years=[2005]
# run_years=[2006]
# run_years=[2007]
# run_years=[2008]
# run_years=[2009]
# run_years=[2010]
# run_years=[2011]
# run_years=[2012]
# run_years=[2013]
# run_years=[2014]
# run_years=[2015]

# In[]
if __name__ == '__main__':
    
    data_dir='D:\hh\@进行中的工作\ESSD_Wheat_Multi-vars\分区标定'
    # data_dir='H:\hh\@进行中的工作\ESSD_Wheat_Multi-vars\分区标定'
    # Some site parameters
    sited = WOFOST71SiteDataProvider(
            IFUNRN = 0,#Indicates whether non-infiltrating fraction of rain is a function of storm size (1)
            SSMAX  = 0.000000, #Maximum depth of water that can be stored on the soil surface [cm]
            WAV    = 20.000000, #Initial amount of water in total soil profile [cm]
            NOTINF = 0.000000,#Maximum fraction of rain not-infiltrating into the soil [0-1], default 0.
            SSI    = 0.000000, #Initial depth of water stored on the surface [cm]
            SMLIM  = 0.500000,# Initial maximum moisture content in initial rooting depth zone [0-1], default 0.4
            CO2    = 360.)#Atmospheric CO2 level (ppm), default 360.
    
    cropfile = os.path.join(data_dir, 'cropdata')
    
    

    soilfile = os.path.join(data_dir, 'soildata', 'zhengzhou.soil')
    soild = CABOFileReader(soilfile)     
    # par_lst=[sited,cropd,soild]#土壤数据还需要更新！！！
    
    ds=gdal.Open('D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/11省分区范围/tif_01dg_clip.tif')
    data=ds.ReadAsArray()
    gt=ds.GetGeoTransform()
    
    gdf = gpd.read_file(data_dir+'/站点生育期整理_3期标定_join_yield_xie_and_soil_he_临县补齐.shp',encoding ='UTF-8')
    
    ds_mask=gdal.Open('D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/Mask/mask01/union_mask_01deg.tif')
    mask=ds_mask.ReadAsArray()
    gt_mask=ds_mask.GetGeoTransform()
    ds_mean_yield=gdal.Open('D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/县域_mean_yield.tif')
    mean_yield=ds_mean_yield.ReadAsArray()
    gt_yield=ds_mean_yield.GetGeoTransform()
    yield_class_dct=np.load('final_dct.npy',allow_pickle=True).item()
    nums=[]
    for x in range(ds.RasterXSize):
        for y in range(ds.RasterYSize):
            lon=gt[0] + x * gt[1]
            lat=gt[3] + y * gt[5]
            #用一个mask去掩膜，但分辨率不太一样，近似10倍
            if lon<gt_mask[0] or  lon>gt_mask[0]+ds_mask.RasterXSize*gt_mask[1] \
                or lat>gt_mask[3] or lat<gt_mask[3]+ds_mask.RasterYSize*gt_mask[5]:
                    continue
            xx=int((lon-gt_mask[0])/gt_mask[1]) 
            yy=int((lat-gt_mask[3])/gt_mask[5])
            xx=5 if xx<5 else xx
            yy=5 if yy<5 else yy
            data_mask=mask[yy-5:yy+5,xx-5:xx+5]
            
            #判断是否在省界内
            v=data[y,x]
            #判断是否存在气象
            site="%3.1f_%2.1f.008"%(np.floor(lon*10)/10.,np.floor(lat*10)/10.)
            weatherfile = os.path.join(data_dir, 'weather', site)
            exist=os.path.exists(weatherfile)
            
            if v==255 or np.sum(data_mask)==0 or (not exist):#255是空，<20是新疆一带or v<180
                continue
            df_v=gdf.loc[v]
            
            cropd = YAMLCropDataProvider(cropfile)
            if v <18:
                cropd.set_active_crop('wheat_xj', "Winter_wheat_105")
                print('这是新疆')
            else:
                cropd.set_active_crop('wheat', "Winter_wheat_105")            
            
            #更新参数
            cropd=dict(cropd)
            VS=2.4665*lat - 35.751 #根据山东的结果回归拟合
            if VS<0:
                VS=0.1
            VB=VS/5
            cropd.update({'VERNBASE': VB,
                          'VERNSAT': VS,
                          'DLO': 15})
            par_lst=[sited,cropd,soild]
            
            #判断yield_level
            ix=int((lon-gt_yield[0])/gt_yield[1]) 
            iy=int((lat-gt_yield[3])/gt_yield[5])
            pac_yield=mean_yield[iy,ix]
            code=int(df_v['区站号'])#最新的用了研究区裁剪，这里偏多，下一步用try筛选
            if code in replace_code.keys():
                code=replace_code[code]

            try:
                yields=yield_class_dct[code]
            except:
                continue
            if len(yields)==1:
                level=1
            elif pac_yield>20000:#空值
                level=2
            else:
                if pac_yield<yields[1]:
                    level=1
                elif pac_yield >yields[3]:
                    level=3
                else:
                    level=2
            nums.append([x,y,lon,lat,df_v,par_lst,level])
    # nums=np.load('D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/空跑/nums.npy')
    e1 = datetime.datetime.now()
    # with open('D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/空跑/TEST%d-%d.csv'%(run_years[0],run_years[-1]), 'w') as csv_file:
    csv_file=open('D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/空跑/TEST%d-%d_V2-3.csv'%(run_years[0],run_years[-1]), 'w')
    csv_write = csv.writer(csv_file)
    p = Pool(6)
    
    count=1
    for i in nums:
        # p.apply_async(sayHi, (i,),callback=mycallback)
        last=p.apply_async(func=sayHi, args=(i,),callback=mycallback)
        count=count+1
        print(count,len(p._cache))
        
        if len(p._cache) > 500:
            print("waiting for cache to clear...")
            last.wait()

            p.close()
            p.join()
            p = Pool(6)
    e2 = datetime.datetime.now()
    print((e2-e1)) 
    
    time.sleep( 600 )#似乎不这样会出现文件关闭还没运行完的情况
    csv_file.close() 
        
        