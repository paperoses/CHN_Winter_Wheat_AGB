#!/usr/bin/env python
# coding: utf-8

'''
针对新疆低温，DTSM01调到负值来校准
物候，TMPFTB，和TMNFTB校准生长量
此外改下Tbase

v2-3: 修改了kdvs优化方式，不变动叶片
v2-1：增加分区内分级
V4版结果:三年标定的筛选了单产数据
    单年标定的在上述基础上补充了部分空值，基本采用前后空值（前六年联合标定效果不好）
'''
# In[1]:


import os
import sys
import csv
import nlopt
import time
import datetime
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
import geopandas as gpd
import scipy.stats as ss

import matplotlib.colors as mcolors
from pathlib import Path
from pcse.util import Afgen
from pcse.models import Wofost71_PP,Wofost71_WLP_FD
from pcse.base import ParameterProvider
from pcse.fileinput import YAMLAgroManagementReader, CABOFileReader
from pcse.fileinput import YAMLCropDataProvider#,CABOWeatherDataProvider
sys.path.append('D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/ModifiedClass')
#修改了取消不同年份气象文件经纬度海拔一致性检查
from cabo_weather import CABOWeatherDataProvider
from pcse.util import WOFOST71SiteDataProvider
from netCDF4 import Dataset,date2index
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
from multiprocessing import Pool #导入进程池

# In[]:
#生成路径

# fps=['cropdata/calibrated', 'amgt/calibrated', '标定png', 'weather']
# for fp in fps:
#     Path(fp).mkdir(parents=True, exist_ok=True)


# In[2]:
# 气象文件自动生成模块
ranges = {"LAT": (-90., 90.),
          "LON": (-180., 180.),
          "ELEV": (-300, 6000),
          "IRRAD": (0., 40e6),
          "TMIN": (-50., 60.),
          "TMAX": (-50., 60.),
          "VAP": (0.06, 199.3),  # hPa, computed as sat. vapour pressure at -50, 60 Celsius
          "RAIN": (0, 25),
          "E0": (0., 2.5),
          "ES0": (0., 2.5),
          "ET0": (0., 2.5),
          "WIND": (0., 100.),
          "SNOWDEPTH": (0., 250.),
          "TEMP": (-50., 60.),
          "TMINRA": (-50., 60.)}
class creatWeather(object):
    
    def __init__(self, mylon, mylat, year_begin, year_end, outputfile):
        self.mylon = mylon
        self.mylat = mylat
        self.year_begin = year_begin
        self.year_end = year_end
        self.outputfile = outputfile
        self.inputfile='//server512-1/ITPCAS_NC_10_10_grid_1979-2018/six_in_one'
        
    def  __call__(self):
        station_number=1
        
        site="%3.1f_%2.1f"%(np.floor(self.mylon*10)/10.,np.floor(self.mylat*10)/10.)
        c1=-0.18
        c2=-0.55
        elev_nc=Dataset(os.path.join(self.inputfile,"elev_ITPCAS-CMFD_V0106_B-01_010deg.nc"))
        elev=elev_nc.variables["elev"][:]

        x_lon=int(self.mylon)
        y_lat=int(self.mylat)
        parnames = {"rad","tmin","tmax","hum","wind","prec"}
        createVar = locals()

        dataset=Dataset(os.path.join(self.inputfile,"%03d_%03d.nc"%(x_lon,y_lat)))
        for par in parnames:    
            createVar[par]=dataset.variables[par][:]
#         print(createVar["rad"])
        #计算所在1度*1度格网上的经纬度索引值（用于变量）
        x=int(round((self.mylon-x_lon-0.05)/0.1+0.0000000001))
        y=int(round((self.mylat-y_lat-0.05)/0.1+0.0000000001))
        #计算所在整幅影像上的经纬度索引值（用于高程）
        x_e=int(round((self.mylon-70.05)/0.1+0.0000000001))
        y_e=int(round((self.mylat-15.05)/0.1+0.0000000001))
        alt=elev[0,y_e,x_e]

        times=dataset.variables["time"]

        year = self.year_begin-1
        while year <= self.year_end:
            t = datetime.datetime.strptime('%d-01-01 10:30:00'%(year), '%Y-%m-%d %H:%M:%S')
            index=date2index(t,times)

            f=open(os.path.join(self.outputfile,'weather',site+".%s"%(str(year)[-3:])),"w+")
#             print (year,"begins for site",site)
            f.write("*------------------------------------------------------------*"+"\n"
                        +'*'+"%12s"%("Country: ")+"China"+"\n"
                        +'*'+"%12s"%("Station: ")+site+"\n"
                        +'*'+"%12s"%("Year: ")+"%d"%(year)+"\n"
                        +'*'+"%12s"%("Origin: ")+"ITPCAS China Meteorological Forcing Dataset"+"\n"
                        +'*'+"%12s"%("Author: ")+"CAU, HuangHai"+"\n"
                        +'*'+"%12s"%("Longitude: ")+"%.2f"%(self.mylon)+" E"+"\n"
                        +'*'+"%12s"%("Latitude: ")+"%.2f"%(self.mylat)+" N"+"\n"
                        +'*'+"%12s"%("Elevation: ")+"%.2f"%(alt)+" m"+"\n"
                        +'*'+"%12s"%("Columns: ")+"\n"
                        +'*'+"%12s"%("======== ")+"\n"
                        +'*'+"  station number"+"\n"
                        +'*'+"  year"+"\n"
                        +'*'+"  day"+"\n"
                        +'*'+"  irradiation (kJ·m-2·d-1)"+"\n"
                        +'*'+"  minimum temperature (degrees Celsius)"+"\n"
                        +'*'+"  maximum temperature (degrees Celsius)"+"\n"
                        +'*'+"  vapour pressure (kPa)"+"\n" 
                        +'*'+"  mean wind speed (m·s-1)"+"\n" 
                        +'*'+"  precipitation (mm·d-1)"+"\n" 
                        +'**'+" WCCDESCRIPTION="+site+", China"+"\n" 
                        +'**'+" WCCFORMAT=2"+"\n" 
                        +'**'+" WCCYEARNR="+"%d"%(year)+"\n" 
                        +"*------------------------------------------------------------*"+"\n"
                        +"%.2f  %.2f  %.2f  %.2f  %.2f\n"%(self.mylon, self.mylat, alt, c1, c2)
                        )
            if (year % 4) == 0 and (year % 100) != 0 or (year % 400) == 0:
                #print year,u"闰年"
                numdays = 366
            else:
                numdays = 365
            for d in range(0,numdays):
                if createVar["prec"][index+d,y,x]<0:
                   createVar["prec"][index+d,y,x]=0 
                if createVar["prec"][index+d,y,x]>250:#模型中降水不超过250mm
                   createVar["prec"][index+d,y,x] =200
                if createVar["hum"][index+d,y,x]<0.06:#模型中VAP不超过0.06-199.3
                   createVar["hum"][index+d,y,x] =0.1
                if createVar["tmin"][index+d,y,x]<-50:#模型中tmin不超过-50
                   createVar["tmin"][index+d,y,x] =-49
                if createVar["tmax"][index+d,y,x]>60:#模型中tmax不超过60
                   createVar["tmax"][index+d,y,x] =59
                if createVar["wind"][index+d,y,x]>100:#模型中wind不超过100
                   createVar["wind"][index+d,y,x] =99    
                try:
                    f.write("%d"%(station_number)+"\t"+"%d"%(year)+"\t"+"%3d"%(d+1)+"\t"
                                +"%5d"%(round(createVar["rad"][index+d,y,x]))+"\t"
                                +"%5.1f"%(round(createVar["tmin"][index+d,y,x]*10)/10)+"\t"
                                +"%5.1f"%(round(createVar["tmax"][index+d,y,x]*10)/10)+"\t"
                                +"%5.3f"%(round(createVar["hum"][index+d,y,x]*1000)/1000)+"\t"
                                +"%4.1f"%(round(createVar["wind"][index+d,y,x]*10)/10)+"\t"
                                +"%4.1f"%(round(createVar["prec"][index+d,y,x]*10)/10)+"\n")
                except:#masked value
                    return -99
            f.close()
            year += 1
        dataset.close()   
        elev_nc.close()
        return os.path.join(self.outputfile,'weather',site) 





# In[77]:
# 管理文件生成模块

def year_num(d,lst):
    '''主要针对冬小麦出苗在前一年，lst第一项是出苗的doy'''
    dd=2 if d==lst[0] else 1
    return dd
    
class creatAgro():    
    def __init__(self, code, year, calendar_doy, outputfile, level2, SM=0.15, irrigation=5.0):
        self.code = code
        self.year=year
        self.calendar_doy=calendar_doy
        self.outputfile = outputfile 
        self.level2=level2
        self.SM=SM
        self.irrigation=irrigation
    
    def __call__(self):

        em,flw,mt=[datetime.datetime(self.year-year_num(d,self.calendar_doy),12,31)+datetime.timedelta(d) for d in self.calendar_doy]
        file=os.path.join(self.outputfile,'amgt','%d_%d_V2-3_L%d.amgt'%(self.code,self.year,self.level2))
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

        return file
    




# In[105]:
# 模型运行以及代价函数模块

class ModelRerunner1(object):
    """for DVS
    """
    parameters = ["TSUM1", "TSUM2"]
    data_dir = './'
    
    def __init__(self, code, params, wdp, years,calendar_doys,level2):
        self.code = code
        self.params = params
        self.wdp = wdp
        self.years = years
        self.calendar_doys=calendar_doys
        self.level2=level2
        self.VB,self.VS=self.params['VERNBASE'],self.params['VERNSAT']

        # self.out=[]
        
    def __call__(self, par_values):
        out=[]
        # Check if correct number of parameter values were provided
        if len(par_values) != len(self.parameters):
            msg = "Optimizing %d parameters, but only %d values were provided!" % (len(self.parameters), len(par_values))
            raise RuntimeError(msg)
        for year in self.years:
            calendar_doy=self.calendar_doys[self.years.index(year)]
            agromanagement_file =creatAgro(self.code, year, calendar_doy, self.data_dir, self.level2)()# os.path.join(self.data_dir, 'amgt', '%d_%d_V2-3_L0.amgt'%(self.code,year))
            agro=YAMLAgroManagementReader(agromanagement_file)
            # Clear any existing overrides
            self.params.clear_override()
            
            # for parname, value in zip(['VERNBASE','VERNSAT'], [self.VB,self.VS]):
            #     self.params.set_override(parname, value)            
            # Set overrides for the new parameter values
            for parname, value in zip(self.parameters, par_values):
                self.params.set_override(parname, value)
            # print('VB:%d,VS:%d'%(self.params['VERNBASE'],self.params['VERNSAT']))
            # print('SMW:%f'%self.params['SMW'])
            # Run the model with given parameter values
            wofost = Wofost71_PP(self.params, self.wdp, agro)
            wofost.run_till_terminate()
            df = pd.DataFrame(wofost.get_output())
            df.index = pd.to_datetime(df.day)
            ix1 = (df.DVS >0.99) & (df.DVS <1.01)
            #模型结束还没开花，近似一天DVS+0.02
            if len(df.loc[ix1].day)==0:
                ix1=(df.DVS ==df.DVS.max())
                doy_flw=df.loc[ix1].day.max().timetuple().tm_yday+int((1-df.DVS.max())/0.02)
            else:
                doy_flw=df.loc[ix1].day.max().timetuple().tm_yday 
            ix2=(df.DVS ==df.DVS.max())
            doy_mat=df.loc[ix2].day.max().timetuple().tm_yday
            #模型结束还未成熟
            if df.DVS.max()<1.98:
               doy_mat+= int((2-df.DVS.max())/0.02)
            out.append([doy_flw,doy_mat])
        return out

class ModelRerunner2(object):
    """for DVS 2
    """
    parameters = ["TSUM1", "TSUM2", 'DTSM']
    data_dir = './'
    
    def __init__(self, code, params, wdp, years,level2):
        self.code = code
        self.params = params
        self.wdp = wdp
        self.years = years
        self.level2=level2
        self.VB,self.VS=self.params['VERNBASE'],self.params['VERNSAT']

        # self.out=[]
        
    def __call__(self, par_values):
        out=[]
        # Check if correct number of parameter values were provided
        if len(par_values) != len(self.parameters):
            msg = "Optimizing %d parameters, but only %d values were provided!" % (len(self.parameters), len(par_values))
            raise RuntimeError(msg)
        for year in self.years:
            #因为在第一次作物物候标定时已经生成了文件，此处直接读取
            agromanagement_file = os.path.join(self.data_dir, 'amgt', '%d_%d_V2-3_L%d.amgt'%(self.code,year,self.level2))
            agro=YAMLAgroManagementReader(agromanagement_file)
            # Clear any existing overrides
            self.params.clear_override()
            
            # for parname, value in zip(['VERNBASE','VERNSAT'], [self.VB,self.VS]):
            #     self.params.set_override(parname, value)            
            # Set overrides for the new parameter values
            for parname, value in zip(self.parameters, par_values):
                if parname=='DTSM':
                    self.params.set_override('DTSMTB', [value,0,30.0,30.0,45.0,30.0])
                else:
                    self.params.set_override(parname, value)
            # print('VB:%d,VS:%d'%(self.params['VERNBASE'],self.params['VERNSAT']))
            # print('SMW:%f'%self.params['SMW'])
            # Run the model with given parameter values
            wofost = Wofost71_PP(self.params, self.wdp, agro)
            wofost.run_till_terminate()
            df = pd.DataFrame(wofost.get_output())
            df.index = pd.to_datetime(df.day)
            ix1 = (df.DVS >0.99) & (df.DVS <1.01)
            #模型结束还没开花，近似一天DVS+0.02
            if len(df.loc[ix1].day)==0:
                ix1=(df.DVS ==df.DVS.max())
                doy_flw=df.loc[ix1].day.max().timetuple().tm_yday+int((1-df.DVS.max())/0.02)
            else:
                doy_flw=df.loc[ix1].day.max().timetuple().tm_yday 
            ix2=(df.DVS ==df.DVS.max())
            doy_mat=df.loc[ix2].day.max().timetuple().tm_yday
            #模型结束还未成熟
            if df.DVS.max()<1.98:
               doy_mat+= int((2-df.DVS.max())/0.02)
            out.append([doy_flw,doy_mat])
        return out


class ObjectiveFunctionCalculator1(object):
    """Computes the objective function for yield or DVS
    """
    
    def __init__(self, code, params, wdp, years, observations,calendar_doys,level2,dvs_2rd=False):
        if dvs_2rd:
            self.modelrerunner = ModelRerunner2(code, params, wdp, years,level2)
        else:
            self.modelrerunner = ModelRerunner1(code, params, wdp, years,calendar_doys,level2)
        self.observations = observations
        self.n_calls = 0
       
    def __call__(self, par_values, grad=None):
        """Runs the model and computes the objective function for given par_values.
        
        The input parameter 'grad' must be defined in the function call, but is only
        required for optimization methods where analytical gradients can be computed.
        """
        self.n_calls += 1
        print(".", end="")
        # Run the model and collect output
        self.simulations = self.modelrerunner(par_values)            
        obj_func=np.mean(np.sqrt((np.array(self.simulations)-np.array(self.observations))**2))
        # print(self.simulations,'\n',self.observations,obj_func,'\n')
        return obj_func

class ModelRerunner3(object):
    """标定作物参数
    """
    parameters = ["SLATB", "SPAN", "AMAXTB", "TMPFTB", "TMNFTB", "CVO", 'FLTB', 'FSTB', 'FOTB']
    data_dir = './'
    
    def __init__(self, code, params, wdp, years, tsum=None, DTSM=0, level='PP',level2=1):
        self.code = code
        self.params = params
        self.wdp = wdp
        # self.agro = agro
        self.years = years
        self.tsum=tsum
        self.DTSM=DTSM
        self.level = level
        self.level2=level2
        self.VB,self.VS=self.params['VERNBASE'],self.params['VERNSAT']
        
    def __call__(self, par_values,return_df=False):
        out=[]
        new_pars=[[],0,[],[],[],0,[],[],[]]
        # Check if correct number of parameter values were provided
        if len(par_values) != len(self.parameters)-2:#'FLTB', 'FSTB', 'FOTB'中优化同一个参数
            msg = "Optimizing %d parameters, but only %d values were provided!" % (len(self.parameters), len(par_values))
            raise RuntimeError(msg)
        for year in self.years:
            #此时不需要修改灌溉，直接读就行，前面生成过
            agromanagement_file = os.path.join(self.data_dir, 'amgt', '%d_%d_V2-3_L%d.amgt'%(self.code,year,self.level2))
            agro=YAMLAgroManagementReader(agromanagement_file)    
            # Clear any existing overrides
            self.params.clear_override()
            # for parname, value in zip(['VERNBASE','VERNSAT'], [self.VB,self.VS]):
            #     self.params.set_override(parname, value) 
            # Set overrides for the new parameter values
            SLATB=[0.0, par_values[0], 0.5, par_values[0], 2.0, par_values[0]]
            AMAXTB=[0.0, par_values[2], 1.0, par_values[2], 1.3, par_values[2], 2.0, par_values[2]]
            TMPFTB=[0.0+par_values[4], 0.01, 10.0, 0.6, 15.0, 1.0, par_values[3], 1.0, 35.0, 0.0]#同步将TMNFTB低温阈值放进来
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
                # print('LAI max:%.2f'%df.LAI.max())
                
                # if abs(df.TWSO.max()-3537)<1 and abs(df.TAGP.max()-6198)<1:
                #     df.to_csv('D:/temp/2007.csv',index=False)
                # elif abs(df.TWSO.max()-2325)<1 and abs(df.TAGP.max()-4191)<1:
                #     df.to_csv('D:/temp/2008.csv',index=False)
                # elif abs(df.TWSO.max()-1826)<1 and abs(df.TAGP.max()-3260)<1:
                #     df.to_csv('D:/temp/2009.csv',index=False)                 
        return out

class ModelRerunner4(object):
    """recaliobrate irragation"""
    
    parameters = ["SLATB", "SPAN", "AMAXTB", "TMPFTB", "TMNFTB", "CVO", 'FLTB', 'FSTB', 'FOTB']
    data_dir = './'
    
    def __init__(self, code, params, wdp, years, tsum=None, DTSM=0, level='PP',level2=1,x=None,calendar_doys=None):
        self.code = code
        self.params = params
        self.wdp = wdp
        # self.agro = agro
        self.years = years
        self.tsum=tsum
        self.DTSM=DTSM
        self.level = level  
        self.level2 = level2          
        self.par_values=x
        self.calendar_doys=calendar_doys
        self.VB,self.VS=self.params['VERNBASE'],self.params['VERNSAT']#可以取消了，通过外层写入了
        
    def __call__(self, irragation, return_df=False):
        par_values=self.par_values
        out=[]
        new_pars=new_pars=[[],0,[],[],[],0,[],[],[]]

        ny=len(self.years)

        my_year=list(self.years)
        irr=np.array(irragation).reshape([2,ny]).T
        irr_dct=dict(zip(my_year,irr))
        
        for year in self.years:
            calendar_doy=self.calendar_doys[self.years.index(year)]
            agromanagement_file =creatAgro(self.code, year, calendar_doy, self.data_dir,self.level2, irr_dct[year][0],irr_dct[year][1])()
            agro=YAMLAgroManagementReader(agromanagement_file.replace('????','%d'%year))    
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
    
class ObjectiveFunctionCalculator2(object):
    """TWSO & AGB
    当calibrate_irragation=True时，calendar_doys必须给定
    """
    
    def __init__(self, code, params, wdp, years, observations,HI, calendar_doys=None, tsum=None, DTSM=0,\
                 recalibrate_tsum=False,level='PP',level2=1,param_names=None,calibrate_irragation=False,x=None):
        if calibrate_irragation:
            self.modelrerunner = ModelRerunner4(code,  params, wdp, years, tsum, DTSM,level,level2,x,calendar_doys)
        else:
            self.modelrerunner = ModelRerunner3(code, params, wdp, years, tsum, DTSM, level, level2)
        self.observations = np.array(observations)
        self.HI=HI
        self.n_calls = 0
       
    def __call__(self, par_values, grad=None):
        """Runs the model and computes the objective function for given par_values.
        The costfunction take the HI into consideration by introducing a scale 
        factor alpha for costfunction.
        
        The input parameter 'grad' must be defined in the function call, but is only
        required for optimization methods where analytical gradients can be computed.
        """
        self.n_calls += 1
        print(".", end="")
        # Run the model and collect output
        self.simulations = np.array(self.modelrerunner(par_values)) 
        twso=list(map(int,self.simulations[:,0]))    
        tagp=list(map(int,self.simulations[:,1]))
        lai=list(map(int,self.simulations[:,2]))
             
        sigma_twso = self.observations*0.1
        sigma_tagp = self.observations/self.HI*0.1
        sigma_lai = self.observations*0+1.5
        
        #产量的代价函数
        obj_func1=-ss.multivariate_normal.logpdf(np.array(twso), \
                  mean=self.observations, cov=np.diag(sigma_twso**2))
        #通过收获指数转换产量到地上生物量
        obs_agb=list(map(int,self.observations/self.HI))
        #地上生物量的代价函数
        obj_func2=-ss.multivariate_normal.logpdf(np.array(tagp), \
                  mean=obs_agb, cov=np.diag(sigma_tagp**2)) 
        #叶面积指数最大值的代价函数   
        obj_func3=-ss.multivariate_normal.logpdf(np.array(lai), \
                  mean=self.observations*0+6.5, cov=np.diag(sigma_lai**2))  
            
        obj_func=obj_func1+obj_func2+obj_func3
        
        # print(self.n_calls,'\n',twso,tagp,'\n',list(self.observations),obs_agb,obj_func,'\n')

        return obj_func



# In[115]:
# 标定模块
class  CalibrateManage(object):
    data_dir = 'D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/'
    '''calibration TSUM and other parameters '''
    def __init__(self, code,params, wdp, level='PP',level2=1,key_irr=0.25):
        self.code=code
        self.params=params
        self.wdp=wdp
        self.level=level
        self.level2=level2
        self.key_irr=key_irr
    def __call__(self,run_years,dvs_obs,yield_obs,HI,calendar_doys,tsum1=610,tsum2=900, \
                SLATB=0.0020, SPAN=25, AMAXTB=35, TMPFTB=28, TMNFTB=-2, CVO=0.71, shift_dvs=0.75):
        '''\n对物候进行优化，如果对应level2，3，直接读取level1结果''' 
        if self.level2 >1:
            print('读取物候标定结果')
            x=np.loadtxt(self.data_dir+'/cropdata/calibrated/%d-%d_%d_V2-3_L1.txt'%(run_years[0],run_years[-1],self.code))
            tsum_result={'TSUM1':x[0], 'TSUM2': x[1], 'DTSM':x[2],\
                         'optimum_value':-1, \
                         'result_code':-1, \
                         'ncalls':-1}
            dvs_2rd=True
            DTSM=x[2]
            for year in run_years:
                calendar_doy=calendar_doys[run_years.index(year)]
                _ = creatAgro(self.code, year, calendar_doy, self.data_dir, self.level2)()
            
        else:
            print('标定物候')
            dvs_2rd=False
            TSUM1_range = [150, 1050]
            TSUM2_range = [600, 1850]
            stepsize1 = 50.
            stepsize2 = 50.
            objfunc_calculator = ObjectiveFunctionCalculator1(self.code,self.params, \
                                 self.wdp, run_years, dvs_obs,calendar_doys,level2=self.level2)
            # Start optimizer with the SUBPLEX algorithm for two parameters
            opt = nlopt.opt(nlopt.LN_SBPLX, 2)
            opt.set_min_objective(objfunc_calculator)
            opt.set_lower_bounds([TSUM1_range[0], TSUM2_range[0]])
            opt.set_upper_bounds([TSUM1_range[1], TSUM2_range[1]])
            opt.set_initial_step([stepsize1, stepsize2])
            # Maximum number of evaluations allowed
            opt.set_maxeval(100)
            # Relative tolerance for convergence
            opt.set_ftol_rel(0.001)        
            # Start the optimization with the first guess
            firstguess = [tsum1, tsum2]
            x = opt.optimize(firstguess)
            tsum_result={'TSUM1':x[0], 'TSUM2': x[1], \
                         'optimum_value':opt.last_optimum_value(), \
                         'result_code':opt.last_optimize_result(), \
                         'ncalls':objfunc_calculator.n_calls} 
            # print("\noptimum at TSUM1: %s, TSUM2: %s" % (x[0], x[1]))
            # print("minimum value = %f d/year/dvs"%opt.last_optimum_value())
            # print("result code = %d"%opt.last_optimize_result())
            # #=3表示/ *由于达到了ftol_rel或ftol_abs（上方），优化已停止。* /
            # print("With %d function calls" % objfunc_calculator.n_calls)
            
            '''如果误差大于阈值，则加入DTSM优化'''
            if opt.last_optimum_value()>1.0*len(run_years):
                dvs_2rd=True
                print('重新标定物候')
                TSUM1_range = [150, 1050]
                TSUM2_range = [600, 1850]
                DTSM_range=[-30,10]
                stepsize1 = 50.
                stepsize2 = 50.
                stepsize3 = 5
                objfunc_calculator = ObjectiveFunctionCalculator1(self.code,self.params, \
                                     self.wdp, run_years, dvs_obs,calendar_doys,level2=self.level2,dvs_2rd=True)
                # Start optimizer with the SUBPLEX algorithm for two parameters
                opt = nlopt.opt(nlopt.LN_SBPLX, 3)
                opt.set_min_objective(objfunc_calculator)
                opt.set_lower_bounds([TSUM1_range[0], TSUM2_range[0],DTSM_range[0]])
                opt.set_upper_bounds([TSUM1_range[1], TSUM2_range[1],DTSM_range[1]])
                opt.set_initial_step([stepsize1, stepsize2, stepsize3])
                # Maximum number of evaluations allowed
                opt.set_maxeval(200)
                # Relative tolerance for convergence
                opt.set_ftol_rel(0.001)        
                # Start the optimization with the first guess
                firstguess = [tsum1, tsum2, -15]
                x = opt.optimize(firstguess)
                tsum_result={'TSUM1':x[0], 'TSUM2': x[1], 'DTSM':x[2],\
                             'optimum_value':opt.last_optimum_value(), \
                             'result_code':opt.last_optimize_result(), \
                             'ncalls':objfunc_calculator.n_calls} 
                # print("\noptimum at TSUM1: %.1f, TSUM2: %.1f, DTSM: %.1f" % (x[0], x[1], x[2]))
                # print("minimum value = %f d/year/dvs"%opt.last_optimum_value())
                # print("result code = %d"%opt.last_optimize_result())
                # #=3表示/ *由于达到了ftol_rel或ftol_abs（上方），优化已停止。* /
                # print("With %d function calls" % objfunc_calculator.n_calls) 
            
        

        
        #other params
        SLATB_range = [0.0019, 0.0025]
        SPAN_range = [20, 38]
        AMAXTB_range = [25, 42]
        TMPFTB_04_range = [20, 30]
        TMNFTB_shift_range = [-20, 1]
        CVO_range = [0.55, 0.8]
        shift_dvs_range=[0.65,0.94]
        
        # if self.code>53986:
        #     SLATB_range = [0.0020, 0.0025]
        #     SPAN_range = [22, 38]
        #     AMAXTB_range = [35, 42]
        #     TMPFTB_04_range = [20, 30]
        #     TMNFTB_shift_range = [-5, 1]
        #     CVO_range = [0.55, 0.8]
        #     shift_dvs_range=[0.7,0.999]            
        
        
        stepsize1 = 0.0001
        stepsize2 = 2
        stepsize3 = 0.5
        stepsize4 = 1
        stepsize5 = 0.5
        stepsize6 = 0.005
        stepsize7 = 0.05
 
        '''对作物参数进行优化'''
        print('\n标定作物参数')
        try:
            DTSM=tsum_result['DTSM']
        except:
            DTSM=0
        objfunc_calculator = ObjectiveFunctionCalculator2(self.code,self.params, \
                             self.wdp, run_years, yield_obs, HI,\
                             tsum=[tsum_result['TSUM1'],tsum_result['TSUM2']], \
                             DTSM=DTSM,level=self.level,level2=self.level2)
        # Start optimizer with the SUBPLEX algorithm for 6 parameters
        opt = nlopt.opt(nlopt.LN_SBPLX, 7)
        opt.set_min_objective(objfunc_calculator)
        opt.set_lower_bounds([SLATB_range[0], SPAN_range[0], AMAXTB_range[0], \
                              TMPFTB_04_range[0], TMNFTB_shift_range[0], \
                              CVO_range[0],shift_dvs_range[0]])
        opt.set_upper_bounds([SLATB_range[1], SPAN_range[1], AMAXTB_range[1], \
                              TMPFTB_04_range[1], TMNFTB_shift_range[1], \
                              CVO_range[1],shift_dvs_range[1]])
        opt.set_initial_step([stepsize1, stepsize2,stepsize3, stepsize4, \
                              stepsize5, stepsize6, stepsize7])
        opt.set_maxeval(300)
        # Relative tolerance for convergence
        opt.set_ftol_rel(0.005)
        firstguess = [SLATB, SPAN, AMAXTB, TMPFTB, TMNFTB, CVO, shift_dvs]
        x = opt.optimize(firstguess)
        params_result={'SLA':x[0], 'SPAN':x[1], 'AMAX':x[2], \
                        'TMPFTB_04':x[3], 'TMNFTB_shift':x[4], \
                        'CVO':x[5], 'shift_dvs':x[6], \
                        'optimum_value':opt.last_optimum_value(), \
                        'result_code':opt.last_optimize_result(), \
                        'ncalls':objfunc_calculator.n_calls}  
            
            
        # print("\noptimum at \nSLATB: %.5f, SPAN: %.1f, AMAXTB: %.1f, TMPFTB: %.1f, \nTMNFTB: %.1f, CVO: %.2f, shift_dvs: %.2f" \
        #       % (x[0], x[1], x[2], x[3], x[4], x[5], x[6]))
        # print("minimum value = %f "% opt.last_optimum_value())
        # print("result code = %d"% opt.last_optimize_result())
        # print("With %i function calls" % objfunc_calculator.n_calls) 
        
        #保留目前的结果
        # params_result_b = params_result
        last_b = opt.last_optimum_value()
        x_b = x
            
            
        '''对灌溉进行优化'''
        print('\n优化灌溉量')
        objfunc_calculator = ObjectiveFunctionCalculator2(self.code,self.params, \
                             self.wdp, run_years, yield_obs, HI, calendar_doys, \
                             tsum=[tsum_result['TSUM1'],tsum_result['TSUM2']], \
                             DTSM=DTSM,level=self.level,level2=self.level2, calibrate_irragation=True,x=x)
        # Start optimizer with the SUBPLEX algorithm for 10 parameters
        ny=len(run_years)
        opt = nlopt.opt(nlopt.LN_SBPLX, ny*2)
        opt.set_min_objective(objfunc_calculator)
        opt.set_lower_bounds(np.zeros(ny*2))
        #这个最大值设置很关键，太大了会导致卡死在最大值，SM不经历由大到小的下降过程；太小了没有足够优化空间，对萎蔫含水量大的土壤会出问题
        #这里根据SMW和SMFCF设置了一个估计值
        opt.set_upper_bounds(list(np.zeros(ny)+self.key_irr)+list(np.zeros(ny)+16))
        opt.set_initial_step(list(np.zeros(ny)+0.1)+list(np.zeros(ny)+2))          
        opt.set_maxeval(300)
        # Relative tolerance for convergence
        opt.set_ftol_rel(0.005)
        firstguess = list(np.zeros(ny)+(self.key_irr-0.02))+list(np.zeros(ny)+6)
        x1 = opt.optimize(firstguess)
        #这里针对x1再生成一次管理文件，因为最后一次优化结果不一定最优
        y=0
        for year in run_years:
            _ =creatAgro(self.code, year, calendar_doys[y], self.data_dir, self.level2, SM=x[y], irrigation=x[y+3])()
            y+=1
        

        # print("minimum value = %f "% opt.last_optimum_value())
        # print("result code = %d"% opt.last_optimize_result())
        # print("With %i function calls" % objfunc_calculator.n_calls)             
        
        params_result={'SLA':x[0], 'SPAN':x[1], 'AMAX':x[2], \
                        'TMPFTB_04':x[3], 'TMNFTB_shift':x[4], \
                        'CVO':x[5], 'shift_dvs':x[6],'irragation':x1}             
            
        '''参数二次优化'''
        print('\n再次标定作物参数')
        objfunc_calculator = ObjectiveFunctionCalculator2(self.code,self.params, \
                              self.wdp, run_years, yield_obs, HI,\
                              tsum=[tsum_result['TSUM1'],tsum_result['TSUM2']], \
                              DTSM=DTSM,level=self.level,level2=self.level2)
        # Start optimizer with the SUBPLEX algorithm for 6 parameters
        opt = nlopt.opt(nlopt.LN_SBPLX, 7)
        opt.set_min_objective(objfunc_calculator)
        opt.set_lower_bounds([SLATB_range[0], SPAN_range[0], AMAXTB_range[0], \
                              TMPFTB_04_range[0], TMNFTB_shift_range[0], \
                              CVO_range[0],shift_dvs_range[0]])
        opt.set_upper_bounds([SLATB_range[1], SPAN_range[1], AMAXTB_range[1], \
                              TMPFTB_04_range[1], TMNFTB_shift_range[1], \
                              CVO_range[1],shift_dvs_range[1]])
        opt.set_initial_step([stepsize1, stepsize2,stepsize3, stepsize4, \
                              stepsize5, stepsize6, stepsize7])
        opt.set_maxeval(300)
        # Relative tolerance for convergence
        opt.set_ftol_rel(0.005)
        firstguess = x
        x = opt.optimize(firstguess)           

        
        if opt.last_optimum_value()>last_b:
            print('使用未标定灌溉结果')
            x=x_b
            x1=list(np.zeros(ny)+0.15)+list(np.zeros(ny)+5)#这里要和管理参数生成的默认参数设置契合
            
        params_result={'SLA':x[0], 'SPAN':x[1], 'AMAX':x[2], \
                        'TMPFTB_04':x[3], 'TMNFTB_shift':x[4], \
                        'CVO':x[5], 'shift_dvs':x[6],'irragation':x1} 
            
        # print("\noptimum at \nSLATB: %.5f, SPAN: %.1f, AMAXTB: %.1f, TMPFTB: %.1f, \nTMNFTB: %.1f, CVO: %.2f, shift_dvs: %.2f" \
        #       % (x[0], x[1], x[2], x[3], x[4], x[5], x[6]))
        # print("minimum value = %.2f --> %.2f "% (last_b, opt.last_optimum_value()))
        # print("result code = %d"% opt.last_optimize_result())
        # print("With %i function calls" % objfunc_calculator.n_calls)
        
        if dvs_2rd:
            return {'TSUM':tsum_result, 'DTSM':DTSM,'params':params_result,'min_v':np.min([last_b,opt.last_optimum_value()])}
        else:
            return {'TSUM':tsum_result, 'params':params_result,'min_v':np.min([last_b,opt.last_optimum_value()])}
    

# In[56]:
# ## 生育期和产量数据
# 通过Arcgis spatial join

data = gpd.read_file('站点生育期整理_3期标定_join_yield_xie_and_soil_he_临县补齐_研究区.shp',encoding ='UTF-8')
yield_class_dct=np.load('final_dct.npy',allow_pickle=True).item()
pac_yield_data=gpd.read_file('县域产量2001-2015_研究区裁剪_filter_手动再检查_加均值_补空值.shp',encoding ='UTF-8')
pac_yield_data['PAC']=pac_yield_data['PAC'].astype(int)
pac_yield_data=pac_yield_data.set_index ( ['PAC']).T.to_dict()
# print(data.crs)  # 查看数据对应的投影信息
# print(data.head())  # 查看前5行数据
# data.plot()
# plt.show()#简单展示




# In[10]:
# ## 读取各省冬小麦收获指数
# Xie, G., D. Han and X. Wang, Harvest index and residue factor of cereal crops in China. Journal of China agricultural university, 2011. 16(1): p. 1-8.    
# 省份缺值用全国平均

HI_data=pd.read_excel('各省HI.xls')
HI_data.index=HI_data.省名
HI_dct=HI_data.to_dict()['HI']
# VB_dct=HI_data.to_dict()['VERNBASE']
# VS_dct=HI_data.to_dict()['VERNSAT']


# In[]:
    
def mycallback(x):
    # print(x[:2])
    csv_write.writerow(x)  
    
def sayHi(num):
    data_dir = './'
    i=num[0]
    run_years=num[1]
    # jdct=num[2]

    # print('\n%d'%i)
    code=int(data.loc[i,'区站号'])
    print('\n区站号：%d'%code)
    mylon=data.loc[i,'经度']
    mylat=data.loc[i,'纬度']     
    site="%3.1f_%2.1f"%(np.floor(mylon*10)/10.,np.floor(mylat*10)/10.)#这里要和气象生成一致，前面代码注意修改
    HI=HI_dct[data.loc[i,'ProvinceNa']]
    VS=2.4665*mylat - 35.751 #根据山东的结果回归拟合
    if VS<0:
        VS=0.1
    VB=VS/5
    
    # VB=VB_dct[data.loc[i,'ProvinceNa']]
    # VS=VS_dct[data.loc[i,'ProvinceNa']]
    calendar_doys=[]
    county_yields=[]
    new_run_years=[]
    cw=False#creat weather
    
    station_code=data.loc[i,'区站号']#这张表实际少了一小部分
    yield_class=yield_class_dct[station_code]
    calibrate_result=[code]#主要为了记录最小代价函数值
    for yy in run_years:
        calibrate_result.append(yy)
    for level2 in [1,2,3]:#每个分区内三个等级
        # #############
        # min_v=jdct[code,run_years[0]][2*level2-1]
        # if min_v<70:
        #     calibrate_result.append(level2)
        #     calibrate_result.append(min_v)
        #     print('小于70，跳过再标定')
        #     continue
    
    
        # ###################
        # #临时修补
        # if os.path.exists('D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/cropdata/calibrated/%d-%d_%d_V2-3_L%d.txt'%(run_years[0],run_years[-1],code,level2)):
        #     print('Skiping %d-%d_%d_V2-3_L%d.txt'%(run_years[0],run_years[-1],code,level2) )
        #     continue
        # ###################
    
    
        calendar_doys=[]
        county_yields=[]
        new_run_years=[]
        cw=False#creat weather    
        if len(yield_class)==1 and level2==1:
            pac=yield_class[0]
            yields=pac_yield_data[pac]
            for year in run_years:    
                if yields['%d_yield'%year]==0:
                    continue
                new_run_years.append(year)
                calendar_doys.append(list(map(int,[data.loc[i,'%d_Em'%year],data.loc[i,'%d_Hg'%year]+5,data.loc[i,'%d_Mt'%year]])))
                county_yields.append(int(yields['%d_yield'%year]))
                if not os.path.exists(os.path.join(data_dir,'weather/%s.%s'%(site,str(year)[-3:]))):
                    cw=True
        elif len(yield_class)==5:
            pac=yield_class[2*level2-2]#0,2,4
            yields=pac_yield_data[pac]
            for year in run_years:    
                if yields['%d_yield'%year]==0:
                    continue
                new_run_years.append(year)
                calendar_doys.append(list(map(int,[data.loc[i,'%d_Em'%year],data.loc[i,'%d_Hg'%year]+5,data.loc[i,'%d_Mt'%year]])))
                county_yields.append(int(yields['%d_yield'%year]))
                if not os.path.exists(os.path.join(data_dir,'weather/%s.%s'%(site,str(year)[-3:]))):
                    cw=True
            
        else:
            print('  ')
            continue
 
        year_begin,year_end = new_run_years[0],new_run_years[-1]
        
        #only check the last year for simplication
        
        if cw or (not os.path.exists(os.path.join(data_dir,'weather/%s.%s'%(site,str(year_begin-1)[-3:])))):
            weaCreat=creatWeather(mylon, mylat, year_begin, year_end, data_dir) 
            new_weather=weaCreat()
        
            
        weatherfile = os.path.join(data_dir, 'weather', site)
        # print(weatherfile)
        wdp = CABOWeatherDataProvider(weatherfile)
        cropfile = os.path.join(data_dir, 'cropdata')
        cropd = YAMLCropDataProvider(cropfile)
        if i <18:
            cropd.set_active_crop('wheat_xj', "Winter_wheat_105")
            print('这是新疆')
        else:
            cropd.set_active_crop('wheat', "Winter_wheat_105")
        
        #对春化参数永久固定
        cropd=dict(cropd)
        cropd.update({'VERNBASE': VB,
                      'VERNSAT': VS,
                      'DLO': 15})    
        
        soilfile = os.path.join(data_dir, 'soildata', 'zhengzhou.soil')
        soild = CABOFileReader(soilfile)
        #update soil data accroding He Liang's data
        # print('SMW0:%f'%data.loc[i,'SMW'])
        SMW=data.loc[i,'SMW']
        SMFCF=data.loc[i,'SMFCF']
        soild.update({'SMW': SMW,
                      'SMFCF': SMFCF,
                      'SM0': data.loc[i,'SM0'],
                      'CRAIRC': data.loc[i,'CRAIRC']})
        # Some site parameters
        sited = WOFOST71SiteDataProvider(
                IFUNRN = 0,#Indicates whether non-infiltrating fraction of rain is a function of storm size (1)
                SSMAX  = 0.000000, #Maximum depth of water that can be stored on the soil surface [cm]
                WAV    = 20.000000, #Initial amount of water in total soil profile [cm]
                NOTINF = 0.000000,#Maximum fraction of rain not-infiltrating into the soil [0-1], default 0.
                SSI    = 0.000000, #Initial depth of water stored on the surface [cm]
                #设的大一点，避免出现第二天水就降到很低
                SMLIM  = 0.500000,# Initial maximum moisture content in initial rooting depth zone [0-1], default 0.4
                CO2    = 360.)#Atmospheric CO2 level (ppm), default 360.
    
        # Package everyting into a single parameter object
        
        params = ParameterProvider(cropdata=cropd, sitedata=sited, soildata=soild)
    
        # # Here we define the agromanagement 
        # agromanagement_file = os.path.join(data_dir, 'data/amgt', 'Zhengzhou_2013.amgt')
        # agro = YAMLAgroManagementReader(agromanagement_file)
    
        '''Calibration'''
        dvs_obs=[dd[1:] for dd in calendar_doys] 
        yield_obs=list(map(int,np.array(county_yields)*0.875))
        # print('county_yields:',county_yields)
        #计算一个灌水优化的临界高值,经验设定
        key_irr=SMW+9*(SMFCF-SMW)/10
        calibrate_manage=CalibrateManage(code,params,wdp,'WLP_FD',level2,key_irr)

    
        '''run with calibrated result'''
        
        results=calibrate_manage(new_run_years,dvs_obs,yield_obs,HI,calendar_doys)
        # print('\nrun with calibrated result')
        tsum1,tsum2=results['TSUM']['TSUM1'],results['TSUM']['TSUM2']
        dvs_2rd=False
        try:
            DTSM=results['DTSM']
            dvs_2rd=True
        except:
            DTSM=0
        SLATB, SPAN, AMAXTB, TMPFTB, TMNFTB, CVO, shift_dvs=results['params']['SLA'], \
                results['params']['SPAN'],results['params']['AMAX'], \
                results['params']['TMPFTB_04'],results['params']['TMNFTB_shift'], \
                results['params']['CVO'], results['params']['shift_dvs']
        irragation=results['params']['irragation']
        min_v=results['min_v']
        # print('\ttsum1: %d, tsum2: %d'%(tsum1,tsum2))
        if dvs_2rd:
            print('\tDTSM: %.2f'%DTSM)
        # print('\tirragation:',irragation)
        # print('%.5f,%.1f,%.1f,%.1f,%.1f,%.2f,%.2f'%(SLATB, SPAN, AMAXTB, TMPFTB, TMNFTB, CVO, shift_dvs))
        model=ModelRerunner4(code, params, wdp, new_run_years, tsum=[tsum1,tsum2], \
                             DTSM=DTSM,level='WLP_FD',level2=level2,x=[SLATB, SPAN, AMAXTB, TMPFTB, \
                            TMNFTB, CVO, shift_dvs],calendar_doys=calendar_doys)
        dfs=model(irragation,return_df=True)
    
    
        with open('cropdata/calibrated/%d-%d_%d_V2-3_L%d.txt'%(run_years[0],run_years[-1],code,level2),'w+') as ff:
            for item in [tsum1,tsum2,DTSM,SLATB, SPAN, AMAXTB, TMPFTB, TMNFTB, CVO, shift_dvs]:
                ff.write('%f\n'%item)
        with open('amgt/calibrated/%d-%d_%d_V2-3_L%d.txt'%(run_years[0],run_years[-1],code,level2),'w+') as ff:
            for item in irragation:
                ff.write('%f\n'%item)        
    
    
        day_max=np.max([len(df) for df in dfs])
        fig, axes = plt.subplots(figsize=(11,6))
        axes2=axes.twinx() # 创建第二个坐标轴
        colors=list(mcolors.TABLEAU_COLORS.keys()) #颜色变化
        n_y=0
        for year in new_run_years:
            color=mcolors.TABLEAU_COLORS[colors[n_y]]
            # axes.plot((dfs[n_y].index-dfs[n_y].index[0]).days, dfs[n_y].LAI,'-',color=color)#,label="LAI_sim_%d"%year
            # #just for legend
            # axes2.plot((dfs[n_y].index-dfs[n_y].index[0]).days[0], dfs[n_y].LAI[0],'-',color=color,label="LAI_sim_%d"%year)
            # ######
            # axes2.plot((dfs[n_y].index-dfs[n_y].index[0]).days, dfs[n_y].TWSO,'--',color=color,label="Yield_sim_%d"%year)
            # temp=(dfs[n_y].index-dfs[n_y].index[0]).days[-1]-dfs[n_y].index[-1].dayofyear
            # axes2.plot(temp+dvs_obs[n_y][-1], yield_obs[n_y],'*',color=color, label="Yield_obs_%d"%year)
            # axes2.text(10,np.mean(yield_obs)-400*n_y,'%d HI = %.2f (%d / %d)'% (year,dfs[n_y].TWSO.max()/dfs[n_y].TAGP.max(),  \
            #             dfs[n_y].TWSO.max(),dfs[n_y].TAGP.max()), bbox=dict(facecolor='white', alpha=0), zorder=2)
            
            axes.plot(dfs[n_y].index-datetime.timedelta(365*n_y), dfs[n_y].LAI,'-',color=color)#,label="LAI_sim_%d"%year
            #just for legend
            axes2.plot(dfs[n_y].index[0]-datetime.timedelta(365*n_y), dfs[n_y].LAI[0],'-',color=color,label="LAI_sim_%d"%year)
            ######
            axes2.plot(dfs[n_y].index-datetime.timedelta(365*n_y), dfs[n_y].TWSO,'--',color=color,label="Yield_sim_%d"%year)
            temp=dfs[n_y].index[-1]-datetime.timedelta(dfs[n_y].index[-1].dayofyear+365*n_y)
            axes2.plot(temp+datetime.timedelta(dvs_obs[n_y][-1]), yield_obs[n_y],'*',color=color, label="Yield_obs_%d"%year)            
            axes2.text(temp-datetime.timedelta(50),np.mean(yield_obs)-400*n_y,'%d HI = %.2f (%d / %d)'% (year,dfs[n_y].TWSO.max()/dfs[n_y].TAGP.max(),  \
                        dfs[n_y].TWSO.max(),dfs[n_y].TAGP.max()), bbox=dict(facecolor='white', alpha=0), zorder=2)
            #仅显示月份，隐去年份，因为实际是三年，目前是归到了第一年
            # axes2.xaxis.set_major_locator(mdates.DayLocator(interval=30))
            axes2.xaxis.set_major_formatter(DateFormatter("%m"))
            n_y+=1
            
        fig.legend(loc=9,ncol=len(new_run_years),framealpha=0)
        print('Saving figs:%d-%d_%d_V2-3_L%d.png'%(run_years[0],run_years[-1],code,level2))
        fig.savefig('标定png/%d-%d_%d_V2-3_L%d.png'%(run_years[0],run_years[-1],code,level2),bbox_inches='tight')
        calibrate_result.append(level2)
        calibrate_result.append(min_v)
        
    return calibrate_result
## 主程序-

# In[ ]:
#这里是为了对部分分区重新标定
##################
rd=pd.read_csv('D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/标定记录/add.csv',header=None)
rds=np.array(rd)
a=rds[:,:2]
a1=[tuple(aa) for aa in a]
b=rds[:,4:]#np.max(rds[:,4:],axis=1)
jdct=dict(zip(a1,b))
# #################

if __name__=='__main__':

    nums=[]
    for i in range(177):#202
    # for i in [12]:#202    
        # for years in [[2001,2002,2003,2004,2005],[2006,2007,2008,2009,2010],[2011,2012,2013,2014,2015]]:
        for years in [[2001,2002,2003],[2004,2005,2006],[2007,2008,2009],[2010,2011,2012],[2013,2014,2015]][2:]:
        # for years in [[2001],[2002],[2003],[2004],[2005],[2006],[2007],[2008],[2009],[2010],[2011],[2012],[2013],[2014],[2015]]:
            # for years in [[2004,2005,2006],[2010,2011,2012]]:    
        # for years in [[2001,2002,2003],[2007,2008,2009],[2013,2014,2015]]:    
            # if i<50 and years==[2001,2002,2003,2004,2005]:
            #     continue
            nums.append([i,years])
    e1 = datetime.datetime.now()
    # with open('D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/空跑/TEST%d-%d.csv'%(run_years[0],run_years[-1]), 'w') as csv_file:
    csv_file=open(('标定记录/标定记录_%s'% \
                    (e1)).replace(':','-').replace(' ','-').replace('.','_')+'.csv', 'w')
    csv_write = csv.writer(csv_file)
    p = Pool(14)

    count=1
    for i in nums:
        # ######
        code=data.loc[i[0],'区站号']
        year=i[1][0]
        # try:
        #     maxv= jdct[(code,year)]  
        # except:
        #     maxv=0
        # if maxv>70:#
        if True:
        # #########
            last=p.apply_async(func=sayHi, args=(i,),callback=mycallback)
            count=count+1
            print(count,len(p._cache))
            
            if len(p._cache) > 80:
                print("waiting for cache to clear...")
                last.wait()
    
                p.close()
                p.join()
                p = Pool(14)
    e2 = datetime.datetime.now()
    print((e2-e1)) 
    
    # # time.sleep( 6000 )#似乎不这样会出现文件关闭还没运行完的情况
    # # csv_file.close() 
    
    
    
