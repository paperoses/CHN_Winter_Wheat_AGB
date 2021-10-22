#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import csv
import glob
import nlopt
import socket
import datetime
import logging
import logging.handlers
import pandas as pd
import numpy as np
import geopandas as gpd
# import multiprocessing
from multiprocessing import Pool #导入进程池
from osgeo import gdal

# from functools import lru_cache
from pcse.util import Afgen
from pcse.models import Wofost71_PP,Wofost71_WLP_FD
from pcse.base import ParameterProvider
from pcse.fileinput import YAMLAgroManagementReader, CABOFileReader
from pcse.fileinput import YAMLCropDataProvider,CABOWeatherDataProvider
from pcse.util import WOFOST71SiteDataProvider

# from modules import ModelRerunner

from netCDF4 import Dataset,date2index
from imp import reload 
import warnings
warnings.filterwarnings("ignore")


logger = logging.getLogger(__name__)
# logger.setLevel(level = logging.INFO)
logger.setLevel(level = logging.ERROR)
timestr=datetime.datetime.strftime(datetime.datetime.now(), '%Y%m%d')#%H%M%S

handler = logging.handlers.RotatingFileHandler("G:/结果/log/log_%s.log"%timestr, maxBytes=1024000, backupCount=100)
  
handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

# # 模型运行，更新SPAN，TDWI，TSUM1三参数

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
        try:
            dataset=Dataset(os.path.join(self.inputfile,"%03d_%03d.nc"%(x_lon,y_lat)))
        except:
            return -99
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
    

def year_num(d,lst):
    '''主要针对冬小麦出苗在前一年，lst第一项是出苗的doy'''
    dd=2 if d==lst[0] else 1
    return dd
    
class creatAgro():    
    def __init__(self, code, year, calendar_doy, outputfile,level, SM, irrigation):
        self.code = code
        self.year=year
        self.calendar_doy=calendar_doy
        self.outputfile = outputfile 
        self.level=level
        self.SM=SM
        self.irrigation=irrigation
    
    def __call__(self):

        em,flw,mt=[datetime.datetime(self.year-year_num(d,self.calendar_doy),12,31)+datetime.timedelta(d) for d in self.calendar_doy]
        file=os.path.join(self.outputfile,'amgt','%d_%d_V2-3_L%d.amgt'%(self.code,self.year,self.level))
        if os.path.exists(file):
            return file
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

class ModelRerunner(object):
    """标定作物参数
    """
    parameters = ["IDEM","SPAN", "TSUM1"]
    # if  socket.gethostname()=='Server512-2':
    #     data_dir = 'D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/'
    # elif  socket.gethostname()=='Server512-1':
    #     data_dir = 'D:/_Huanghai/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/'        
    # elif socket.gethostname()=='HH-PC':
    #     data_dir = 'D:/a进行中工作/同化AGB产品/分区标定/'
    # else:
    #     print('Please define data_dir')
    
    
    def __init__(self, code, params, calendar_doys, irr,wdp, year, level='PP', level2=1):
        self.code = code
        self.params = params
        self.calendar_doys = calendar_doys
        self.irr = irr       
        
        self.wdp = wdp
        # self.agro = agro
        self.year = year
        self.level = level
        self.level2 = level2

    

    def __call__(self, par_values0,return_df=True,use_cache=True,_cache={}):
        par_values=[]
        par_values.append(round(par_values0[0]))#EM
        par_values.append(round(par_values0[1]))#Span
        par_values.append(round(par_values0[2])*10)#TSUM1 10以内的变化忽略不计，为了便于优化，优化时除以了10
        lon= np.floor(self.wdp.longitude*10)
        lat=np.floor(self.wdp.latitude*10)#扩大了10倍，仅用于缓存
        par_cache=[self.code,self.level2,lon,lat,par_values[0],par_values[1],par_values[2]]
        # print(par_values,self.year,par_cache)
        # print(par_values)
        
        # Check if correct number of parameter values were provided
        if len(par_values) != len(self.parameters):
            msg = "Optimizing %d parameters, but only %d values were provided!" % (len(self.parameters), len(par_values))
            raise RuntimeError(msg)
        if use_cache:
            try:  
                print('n_cache:',len(_cache),)
                return _cache[tuple(map(round,par_cache))]
            except KeyError:
                logger.info("caches: %d %d %d %d %d %d %d %d"%(len(_cache.keys()),self.code,self.level2,lon,lat,par_values[0],par_values[1],par_values[2]))
                pass    
        data_dir='D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/'
        agromanagement_file =creatAgro(self.code, self.year, self.calendar_doys, data_dir, self.level2, self.irr[0],self.irr[1])()
        agro=YAMLAgroManagementReader(agromanagement_file) 
        date_em=list(agro[0].keys())[0]
        idem=par_values[0]
        if idem:
            t_dict=agro[0][date_em]['CropCalendar']
            t_dict.update({'crop_start_date':date_em+datetime.timedelta(idem)})
            date_end=t_dict['crop_end_date']
            t_dict.update({'crop_end_date':date_end+datetime.timedelta(5)})#由于是earliest end，此处相当于放宽成熟区间
            agro[0][date_em].update({'CropCalendar':t_dict})
            agro[0][date_em+datetime.timedelta(idem)]=agro[0].pop(date_em)#替换出苗日期键值
            agro[1][date_end+datetime.timedelta(5)]=agro[1].pop(date_end)#替换end日期键值
        
        
        # Clear any existing overrides
        self.params.clear_override()
        # Set overrides for the new parameter values
        for parname, value in zip(self.parameters[1:], par_values[1:]):
            self.params.set_override(parname, value)
        # print('model_runner',self.params._cropdata['TSUM2'])
        # Run the model with given parameter values
        if self.level=='PP':
            wofost = Wofost71_PP(self.params, self.wdp, agro)
        elif self.level=='WLP_FD':
            # print(self.wdp.longitude,self.wdp.latitude)
            wofost = Wofost71_WLP_FD(self.params, self.wdp, agro)
        else:
            msg = "Wrong level name!"
            raise RuntimeError(msg)
        wofost.run_till_terminate()
        df = pd.DataFrame(wofost.get_output())
        df.index = pd.to_datetime(df.day)
        if return_df:
            out=df
        else:
            out=[df.TWSO.max(),df.TAGP.max()]  
        # if len(_cache)>30:
        #     _cache={}  #不起作用？？？
        #     # print('???',len(_cache))
        if use_cache:
            _cache[tuple(par_cache)] = out
        
        return out


# # 代价函数

# In[3]:


class CostFunstion(object):
    '''代价函数'''
    def __init__(self, code, params,calendar_doys, irr, wdp, year, obs_doy, RS_data, level='WLP_FD', level2=1,obj_type=1):
        self.code = code
        self.params = params
        self.calendar_doys = calendar_doys
        self.irr = irr           
        self.wdp = wdp
        # self.agro = agro
        self.year = year
        self.obs_doy = obs_doy
        self.RS_data = RS_data
        self.level = level  
        self.level2 = level2  
        self.n_calls = 0
        self.type=obj_type
        
    def __call__(self,par_values, grad=None):
        self.n_calls += 1
        if self.RS_data.max()==0:
            #对于nlopt，代价函数的返回值一定要是numpy.float64
            print('no RS data')
            return -1.0
#         par_values=tuple(map(round,(par_values)))
        # print('cost',self.params._cropdata['TSUM2'])
        m=ModelRerunner(self.code, self.params, self.calendar_doys, self.irr,self.wdp, self.year, self.level,self.level2)

        df = m(par_values,return_df=True)
        diff=(datetime.date(self.year-1,12,31)-df['day'][0]).days
        # print('diff',diff)
        last_day=df.index[-1].day_of_year
        #根据结束日期调整同化结束观测点
        my_obs_doy=self.obs_doy[self.obs_doy<last_day]
        ns=len(my_obs_doy)
        # print(dir(df.index[-1]))

        sim=np.array(df.iloc[my_obs_doy+diff].LAI)
        obs=self.RS_data[:ns]
        
        #如果出现观测后期一直不变或者模拟值为nan，则去掉后期值
        stop1,stop2=len(obs),len(obs)
        for i in range(2,len(obs)):
            if sim[i]!=sim[i]:
                stop1=i
                break
        for i in range(2,len(obs)):        
            if np.min(obs[i:])==np.max(obs[i:]):
                stop2=i+1
                break         
        stop=np.min([stop1,stop2])
        
        sim,obs=sim[:stop],obs[:stop]
        
        # print(sim,'',obs,'',stop)
        
        #如果模拟值是恒定小值，一般是气象参数问题，生长发育停止
        if sim.max()==sim.min():
            obj_func= 99.99

        #如果obs是恒定值
        if obs.max()==obs.min():
            obj_func= 88.88
                
        elif self.type==1:
            #建立归一化的代价函数
            norm_sim=(sim-sim.min())/(sim.max()-sim.min())
            norm_obs=(obs-obs.min())/(obs.max()-obs.min())
            obj_func=np.mean(np.sqrt((np.array(norm_sim)-np.array(norm_obs))**2))
        
        else:
            #建立向量夹角最小化的代价函数
            x,y=sim,obs
            
            # 分别计算两个向量的模：
            l_x=np.sqrt(x.dot(x))
            l_y=np.sqrt(y.dot(y))
            # print('向量的模=',l_x,l_y)
            
            # 计算两个向量的点积
            dian=x.dot(y)
            # print('向量的点积=',dian)
            
            # 计算夹角的cos值：
            cos_=dian/(l_x*l_y)
            # print('夹角的cos值=',cos_)
            
            # 求得夹角（弧度制）：
            # obj_func=np.arccos(cos_)
            
            obj_func=-cos_#单调递减
            
        
        # print('obj:',obj_func)
        
        return obj_func


# # 一个点的Test

# In[4]:

# if  socket.gethostname()=='Server512-2':
#     data_dir = 'D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/'
# elif  socket.gethostname()=='Server512-1':
#     data_dir = 'D:/_Huanghai/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/'        
# elif socket.gethostname()=='HH-PC':
#     data_dir = 'D:/a进行中工作/同化AGB产品/分区标定/'
# else:
#     print('Please define data_dir')
    
# code=53956
# year=2012
# mylon,mylat=110.83,35.4
# site='%3.1f_%2.1f'%(mylon,mylat)
# weatherfile = os.path.join(data_dir, 'weather', site)
# wdp = CABOWeatherDataProvider(weatherfile)
# cropfile = os.path.join(data_dir, 'cropdata')
# cropd = YAMLCropDataProvider(cropfile)
# cropd.set_active_crop('wheat', "Winter_wheat_105")
# soilfile = os.path.join(data_dir, 'soildata', 'zhengzhou.soil')
# soild = CABOFileReader(soilfile)
# sited = WOFOST71SiteDataProvider(
#             IFUNRN = 0,#Indicates whether non-infiltrating fraction of rain is a function of storm size (1)
#             SSMAX  = 0.000000, #Maximum depth of water that can be stored on the soil surface [cm]
#             WAV    = 20.000000, #Initial amount of water in total soil profile [cm]
#             NOTINF = 0.000000,#Maximum fraction of rain not-infiltrating into the soil [0-1], default 0.
#             SSI    = 0.000000, #Initial depth of water stored on the surface [cm]
#             SMLIM  = 0.400000,# Initial maximum moisture content in initial rooting depth zone [0-1], default 0.4
#             CO2    = 360.)#Atmospheric CO2 level (ppm), default 360.

# # Package everyting into a single parameter object
# params = ParameterProvider(cropdata=cropd, sitedata=sited, soildata=soild)
# m=ModelRerunner(code, params, wdp, year, level='WLP_FD')
# par_values1=(28,100,600)
# df1 = m(par_values1,return_df=True)
# print(df1)
# par_values2=(25,200,800)
# df2 = m(par_values2,return_df=True)
# print(df2)


# In[ ]:


# dir(ModelRerunner)


# # 虚拟观测

# In[5]:


# obs_doy=np.array([73,81,89,97,105,113,121,129,137,145])
# diff1=(datetime.date(df1['day'][0].year,12,31)-df1['day'][0]).days
# diff2=(datetime.date(df2['day'][0].year,12,31)-df2['day'][0]).days
# RS_data1=df1.iloc[obs_doy+diff1].LAI
# RS_data2=df2.iloc[obs_doy+diff2].LAI

# Cost1=CostFunstion(code, params, wdp, year, obs_doy, RS_data1)
# Cost2=CostFunstion(code, params, wdp, year, obs_doy, RS_data2)
# Cost1(par_values1),Cost2(par_values2)


# In[11]:
def sayHi(point):
    
    xx,yy=point
    lon,lat=gt[0]+xx*gt[1],gt[3]+yy*gt[5]  
    # print('******',[xx,yy],'******')
    RS_data=data[:,yy,xx]
    '''此处还需要对RS LAI进行筛选判断'''
    if RS_data.max()==0:
        print('zero LAI')
        return list(point)+list([0,0,0])
        
    #根据坐标匹配参数
    code=data1[yy,xx]
    
    

    #判断yield_level
    pac_yield=mean_yield[yy,xx]#空值对应65533,这里似乎有问题，这些县是筛选过的！    
    # level2=1
    if code in replace_code.keys():
        code0=code+0
        code=replace_code[code]
        print('code %d replaced by %d'%(code0,code))
    try:
        yields=yield_class_dct[code]
    except:
        print('No yields info')
        return list(point)+list([0,0,0])
    
    if len(yields)==1:
        level2=1
    elif pac_yield>20000:#空值
        level2=2
    else:
        if pac_yield<yields[1]:
            level2=1
        elif pac_yield >yields[3]:
            level2=3
        else:
            level2=2    
    
    
    
    site="%3.1f_%2.1f"%(np.floor(lon*10)/10.,np.floor(lat*10)/10.)#这里要和气象生成一致，前面代码注意修改
    # site='81.3_43.9'
    weatherfile = os.path.join(data_dir, 'weather', site)
    try:
        wdp = CABOWeatherDataProvider(weatherfile)
        if wdp.lastyear-wdp.firstyear<16:
            0/0
    except:
        re=creatWeather(lon,lat,2000,2016,'D:\hh\@进行中的工作\ESSD_Wheat_Multi-vars\分区标定')()
        if re==-99:
            for i1 in np.arange(-0.1,0.11,0.1):
                for i2 in np.arange(-0.1,0.11,0.1):
                    site="%3.1f_%2.1f"%(np.floor(lon*10)/10.+i1,np.floor(lat*10)/10.+i2)
                    weatherfile = os.path.join(data_dir, 'weather', site)
                    try:
                        wdp = CABOWeatherDataProvider(weatherfile)
                    except:
                        pass                    
                    if wdp.lastyear-wdp.firstyear>=16:
                        # wdp = CABOWeatherDataProvider(weatherfile)
                        break
        wdp = CABOWeatherDataProvider(weatherfile)
    cropfile = os.path.join(data_dir, 'cropdata')
    cropd = YAMLCropDataProvider(cropfile)
    # if int(code) in [51133,51346,51358,51368,51431,51436,51436,51628,51633, \
    #                  51642,51644,51709,51716,51777,51810,51811,51828,51855,51931]:
    #     cropd.set_active_crop('wheat_xj', "Winter_wheat_105")
    #     print('这是新疆')
    # else:
    #     cropd.set_active_crop('wheat', "Winter_wheat_105")
    # cropd.set_active_crop('wheat', "Winter_wheat_105")
    soilfile = os.path.join(data_dir, 'soildata', 'zhengzhou.soil')
    soild = CABOFileReader(soilfile)
    
    
    sited = WOFOST71SiteDataProvider(
                IFUNRN = 0,#Indicates whether non-infiltrating fraction of rain is a function of storm size (1)
                SSMAX  = 0.000000, #Maximum depth of water that can be stored on the soil surface [cm]
                WAV    = 20.000000, #Initial amount of water in total soil profile [cm]
                NOTINF = 0.000000,#Maximum fraction of rain not-infiltrating into the soil [0-1], default 0.
                SSI    = 0.000000, #Initial depth of water stored on the surface [cm]
                SMLIM  = 0.500000,# Initial maximum moisture content in initial rooting depth zone [0-1], default 0.4
                CO2    = 360.)#Atmospheric CO2 level (ppm), default 360.
    

    #读取已经标定的参数
    index=conde_index_dct[int(code)]
    calendar_doys=list(map(int,[gdata.loc[index,'%d_Em'%year],gdata.loc[index,'%d_Hg'%year]+5,gdata.loc[index,'%d_Mt'%year]]))
    irr_file='D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/amgt/calibrated/%d-%d_%d_V2-3_L%d.txt'%(run_years[0],run_years[-1],code,level2)
    
    ny=len(run_years)
    irragation=np.loadtxt(irr_file)
    irr=np.array(irragation).reshape([2,ny]).T
    irr_dct=dict(zip(run_years,irr))
    
    p_file='D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/cropdata/calibrated//%d-%d_%d_V2-3_L%d.txt'%(run_years[0],run_years[-1],code,level2)
    par_values=np.loadtxt(p_file)
    # tsum1,tsum2,DTSM,SLATB, SPAN, AMAXTB, TMPFTB, TMNFTB, CVO, shift_dvs=par_values
    tsum1,tsum2,dtsm,span,cvo=par_values[0],par_values[1],par_values[2],par_values[4],par_values[8]
    parameters = ["SLATB", "SPAN", "AMAXTB", "TMPFTB", "TMNFTB", "CVO", 'FLTB', 'FSTB', 'FOTB',"TSUM1", "TSUM2",'DTSM' ]
    new_pars=[[],0,[],[],[],0,[],[],[],tsum1,tsum2,dtsm]
    SLATB=[0.0, par_values[3], 0.5, par_values[3], 2.0, par_values[3]]
    AMAXTB=[0.0, par_values[5], 1.0, par_values[5], 1.3, par_values[5], 2.0, par_values[5]]
    TMPFTB=[0.0+par_values[7], 0.01, 10.0, 0.6, 15.0, 1.0, par_values[6], 1.0, 35.0, 0.0]#同步将TMNFTB低温阈值放进来
    TMNFTB=[0.0+par_values[7], 0.0, 3.0+par_values[7], 1.0]
    k_dvs = par_values[9]
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
    new_pars[1],new_pars[5]=span,cvo
    new_pars[0],new_pars[2],new_pars[3],new_pars[4],new_pars[6],new_pars[7],new_pars[8]= \
    SLATB,AMAXTB,TMPFTB,TMNFTB,FLTB,FSTB,FOTB
    par_dct=dict(zip(parameters,new_pars))
    # cropd.update(par_dct)
    if index <18:
        cropd.set_active_crop('wheat_xj', "Winter_wheat_105")
        print('这是新疆')
    else:
        cropd.set_active_crop('wheat', "Winter_wheat_105")    
    cropd1=dict(cropd)
    cropd1.update(par_dct)
    VS=2.4665*lat - 35.751 #根据山东的结果回归拟合
    if VS<0:
        VS=0.1
    VB=VS/5
    cropd1.update({'VERNBASE': VB,
                  'VERNSAT': VS,
                  'DLO': 15})    
    #update soil data accroding He Liang's data 
    soild.update({'SMW':gdata.loc[index,'SMW'],
                  'SMFCF': gdata.loc[index,'SMFCF'],
                  'SM0': gdata.loc[index,'SM0'],
                  'CRAIRC': gdata.loc[index,'CRAIRC']})    
    # Package everyting into a single parameter object
    params = ParameterProvider(cropdata=cropd1, sitedata=sited, soildata=soild)
    # print('SPAN:%d, TSUM1:%d'%(params._cropdata['SPAN'],params._cropdata['TSUM1']))
    
    if RS_data.max():
        # 同化：利用遥感数据优化参数
        EM_range=[-3,3]
        SPAN_range=[round(span)-1,round(span)+1]
        TSUM1_range=[round(tsum1/10)-2,round(tsum1/10)+2]

        stepsize1=1
        stepsize2=1
        stepsize3=1
        objfunc_calculator = CostFunstion(code, params, calendar_doys, irr_dct[year], wdp, year, obs_doy, RS_data, level2=level2, obj_type=1)
        # Start optimizer with the SUBPLEX algorithm for 3 parameters
        opt = nlopt.opt(nlopt.LN_SBPLX, 3)
        opt.set_min_objective(objfunc_calculator)
        opt.set_lower_bounds([EM_range[0],SPAN_range[0], TSUM1_range[0]])
        opt.set_upper_bounds([EM_range[1],SPAN_range[1], TSUM1_range[1]])
        opt.set_initial_step([stepsize1, stepsize2, stepsize3])
        # Maximum number of evaluations allowed
        opt.set_maxeval(50)
        # Relative tolerance for convergence
        opt.set_ftol_rel(0.01)        
        # Start the optimization with the first guess
        firstguess = [0,round(span),round(tsum1/10)]
        x = opt.optimize(firstguess)
        x=np.array(list(map(round,x)))
        x1=x*np.array([1,1,10])
        print(x1,'optimum_value:%.3f'%opt.last_optimum_value(),'stop code:',opt.last_optimize_result(),'n_calls',objfunc_calculator.n_calls)
        m=ModelRerunner(code,params, calendar_doys, irr_dct[year], wdp, year, level='WLP_FD', level2=level2)
        # best = m(x,return_df=False)
        # # print('best:',best)
        # if len(best)!=2:
        #     best=[best.TWSO.max(),best.TAGP.max()] 
        # print('best result:',best)
    else:
        print('直接用原始标定后的参数')
        # m=ModelRerunner(code,params, wdp, year, level='WLP_FD')
        # print('No data assimilation')
        # x=[0,span,tsum1]
        return list(point)+list([0,span,tsum1])
    
    #模型运行结果
    df=m(x,return_df=True)
    #如果结果差的离谱，多数是因为气象差距太大，用站点标定结果替代
    if df.TWSO.max()<1000:
        print('TWSO max=%d, bad result, %d_L%d replaced!'%(df.TWSO.max(),code,level2))
        lon,lat=gdata.loc[index,'经度'],gdata.loc[index,'纬度']
        # print('lon=%.2f, lat=%.2f'%(lon,lat))
        site="%3.1f_%2.1f"%(np.floor(lon*10)/10.,np.floor(lat*10)/10.)#这里要和气象生成一致，前面代码注意修改
        weatherfile = os.path.join(data_dir, 'weather', site)
        try:
            wdp = CABOWeatherDataProvider(weatherfile)
            if wdp.lastyear-wdp.firstyear<16:
                0/0
        except:
            re=creatWeather(lon,lat,2000,2016,'D:\hh\@进行中的工作\ESSD_Wheat_Multi-vars\分区标定')()
            if re==-99:
                for i1 in np.arange(-0.1,0.11,0.1):
                    for i2 in np.arange(-0.1,0.11,0.1):
                        site="%3.1f_%2.1f"%(np.floor(lon*10)/10.+i1,np.floor(lat*10)/10.+i2)
                        weatherfile = os.path.join(data_dir, 'weather', site)
                        try:
                            wdp = CABOWeatherDataProvider(weatherfile)
                        except:
                            pass                    
                        if wdp.lastyear-wdp.firstyear>=16:
                            # wdp = CABOWeatherDataProvider(weatherfile)
                            break
        wdp = CABOWeatherDataProvider(weatherfile)      
        m=ModelRerunner(code,params, calendar_doys, irr_dct[year], wdp, year, level='WLP_FD', level2=level2)
        # print(calendar_doys,irr_dct[year],span,tsum1)
        df=m([0,span,tsum1/10.],return_df=True,use_cache=False)#记住优化时TSUM1是缩小了的
    # print('df:',df)
    lai=df.LAI
    twso=df.TWSO
    tagp=df.TAGP
    dvs=df.DVS
    doy=df.day[0].timetuple().tm_yday
    np.savez('G:/结果/V2-3/%d/%d_%d_V2-3.npz'%(year,point[0],point[1]),lai=lai,twso=twso,tagp=tagp,dvs=dvs,doy=doy)
    return list(point)+list(map(round,x))


def mycallback(x):
    # print(x)
    try:
        csv_write.writerow(x)
    except:
        pass

# sayHi(RS_data1),sayHi(RS_data2)


# In[ ]:
    
gdata = gpd.read_file('D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/站点生育期整理_3期标定_join_yield_xie_and_soil_he_临县补齐_研究区.shp',encoding ='UTF-8')
conde_index_dct=dict(zip(gdata['区站号'].astype(int),gdata.index)) 

# if  socket.gethostname()=='Server512-2':
#     my_dir = 'D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/'
# elif  socket.gethostname()=='Server512-1':
#     my_dir = 'D:/_Huanghai/@进行中的工作/ESSD_Wheat_Multi-vars/'        
# elif socket.gethostname()=='HH-PC':
#     my_dir = 'D:/a进行中工作/同化AGB产品/'
# else:
#     print('Please define data_dir')
    
my_dir = 'D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/'   

data_dir=my_dir+'分区标定/'

obs_doy=np.array([73,81,89,97,105,113,121,129,137,145])

year=2013
for years in [[2001,2002,2003],[2004,2005,2006],[2007,2008,2009],[2010,2011,2012],[2013,2014,2015]]:
    if year in years:
        run_years=years
        continue
e1 = datetime.datetime.now()
ds=gdal.Open(my_dir+'MOD_LAI_2000-2015/BNU_LAI/masked_tif/BNU_MODIS_LAI_%d.tif'%year)
data=ds.ReadAsArray()[:11,:,:]
gt=ds.GetGeoTransform()

ds=gdal.Open('D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/区站号_01deg.tif')
data1=ds.ReadAsArray()
gt1=ds.GetGeoTransform() #分辨率严格为0.01°   

# codeList=np.load(my_dir+'分区标定/11省分区范围/CodeList.npy')

ds_mean_yield=gdal.Open('D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/县域_mean_yield_01deg.tif')
mean_yield=ds_mean_yield.ReadAsArray()
gt_yield=ds_mean_yield.GetGeoTransform()
yield_class_dct=np.load('D:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/分区标定/final_dct.npy',allow_pickle=True).item()
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

mask=np.sum(data,axis=0)
print('\t')


if __name__=='__main__':
    n_cores=14

    # try:
    #     csv_file = open('G:/结果/V2-3/DA_Pars_V2-3_%d.csv'%year, 'w')
    # except:
    #     csv_file = open('E:/hh/@进行中的工作/ESSD_Wheat_Multi-vars/结果/V2-3/DA_Pars_V2-3_%d.csv'%year, 'w')
    
    # csv_write = csv.writer(csv_file)
    # # multiprocessing.freeze_support()
    # p = Pool(18)

    # lst0=[]
    # for lon in np.arange(79.4,122.3,0.1):
    #     for lat in np.arange(47.1,21.7,-0.1):
    #         lst0.append([round(lon*10)/10,round(lat*10)/10])
    # points=[]
    # for l in lst0:
    #     for i in range(10):
    #         for j in range(10):
    #             lon,lat=round((l[0]+i/100.)*100)/100,round((l[1]+j/100.)*100)/100
    #             ix=round((lon-gt[0])/gt[1]) 
    #             iy=round((lat-gt[3])/gt[5])
    #             #print(lon,lat,ix,iy)
    #             if ix>=0 and iy>=0 and iy<2538 and ix<4277 and mask[iy,ix]>0:
    #                 if os.path.exists('G:/结果/V2-3/%d/%d_%d_V2-3.npz'%(year,ix,iy)):
    #                     continue
    #                 points.append([ix,iy])

    # for point in points:
    #     p.apply_async(func=sayHi, args=(point,),callback=mycallback)
    # p.close()
    # p.join()

    # csv_file.close()
    # e2 = datetime.datetime.now()
    # print(e2-e1)
    
    lst0=[]
    for lon in np.arange(79.4,122.3,0.1):
        for lat in np.arange(47.1,21.7,-0.1):
            lst0.append([round(lon*10)/10,round(lat*10)/10])
    points=[]
    for l in lst0:
        for i in range(10):
            for j in range(10):
                lon,lat=round((l[0]+i/100.)*100)/100,round((l[1]+j/100.)*100)/100
                ix=round((lon-gt[0])/gt[1]) 
                iy=round((lat-gt[3])/gt[5])
                #print(lon,lat,ix,iy)
                if ix>=0 and iy>=0 and iy<2538 and ix<4277 and mask[iy,ix]>0:
                    if os.path.exists('G:/结果/V2-3/%d/%d_%d_V2-3.npz'%(year,ix,iy)):
                        continue
                    points.append([ix,iy]) 
                    
    files=glob.glob('G:/结果/V2-3/DA_Pars_V2-3_%d*.csv'%year)
    csv_file = open('G:/结果/V2-3/DA_Pars_V2-3_%d-%d.csv'%(year,len(files)+1), 'w')
    csv_write = csv.writer(csv_file)
    
    s_file='G:/结果/V2-3/start/%d.txt'%year
    if not os.path.exists(s_file):
        f=open(s_file,'w')
        f.write('%d'%1)
        f.close()
        start=1
    else:
        start=int(np.loadtxt(s_file))+1
        f=open(s_file,'w')
        f.write('%d'%(start))
        f.close()
        
    total=len(points)
    s=int((start-1)*total/n_cores)
    
    n=s
    for point in points[::-1][s:]:#逆序是因为西侧区域小麦稀疏，缓存利用率低，速度慢
        n+=1
          
        print("%d / %d"%(n,total))
        ix,iy=point
        if os.path.exists('G:/结果/V2-3/%d/%d_%d_V2-3.npz'%(year,ix,iy)):
            continue     
        try:
            mycallback(sayHi(point))
        except:
            print('wrong point: %s'%point)
    csv_file.close()
    
    