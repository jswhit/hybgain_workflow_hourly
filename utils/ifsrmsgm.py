from netCDF4 import Dataset
import numpy as np
import os,sys,dateutils
from bisect import bisect

def getmean(diff,coslats):
    meancoslats = coslats.mean()
    mean = np.empty(diff.shape[0],diff.dtype)
    for k in range(diff.shape[0]):
        mean[k] = (coslats*diff[k]).mean()/meancoslats
    return mean

def getmean2d(diff,coslats):
    meancoslats = coslats.mean()
    return (coslats*diff).mean()/meancoslats

date1 = sys.argv[1]
date2 = sys.argv[2]
expt = sys.argv[3]
label = sys.argv[4]

ifsanldir = '/work2/noaa/gsienkf/whitaker/ifsanal'
expbasedir = '/work2/noaa/gsienkf/whitaker'
hr=6
prefix='sfg2'

levtop=125; levbot=925 # plot limits

dates = dateutils.daterange(date1,date2,6)

# save time series to netcdf file

exptdir = os.path.join(expbasedir+'/'+expt,dates[0])
ufsfcst = Dataset(os.path.join(exptdir,'%s_%s_fhr%02i_ensmean'%(prefix,dates[0],hr)))
nlats = len(ufsfcst.dimensions['grid_yt'])
nlons = len(ufsfcst.dimensions['grid_xt'])
nlevs = len(ufsfcst.dimensions['pfull'])
LEVS=nlevs; RES=nlons//4
ufsfcst.close()
outfile = 'ifsdiff_%s_globalmean' % label
ncout = Dataset(outfile+'.nc','w')
plevs = ncout.createDimension('plevs',nlevs)
lats = ncout.createDimension('lat',nlats)
timed = ncout.createDimension('time',len(dates))
plevsm = ncout.createVariable('plevs_mid',np.float32,'plevs')
plevsd = ncout.createVariable('plevs_down',np.float32,'plevs')
plevsu = ncout.createVariable('plevs_up',np.float32,'plevs')
plevsu.units = 'hPa'
plevsd.units = 'hPa'
plevsm.units = 'hPa'
#latitudes = ncout.createVariable('lat',np.float32,'lat')
#latitudes.units = 'degrees'
times = ncout.createVariable('time',np.float64,'time')
times.units = 'hours since 01-01-01'
winderrsq = ncout.createVariable('winderr',np.float32, ('time','plevs'))
tmperrsq = ncout.createVariable('temperr',np.float32, ('time','plevs'))
tmpbias = ncout.createVariable('tempbias',np.float32, ('time','plevs'))
qerrsq = ncout.createVariable('qerr',np.float32, ('time','plevs'))
qbias = ncout.createVariable('qbias',np.float32, ('time','plevs'))
ncout.exptname = expt
ncout.fhour = hr

ncount=0
for date in dates:
    # 6hourly
    datev = date
    exptdir = os.path.join(expbasedir+'/'+expt,datev)
    ufsfcst = os.path.join(exptdir,'%s_%s_fhr%02i_ensmean'%(prefix,datev,hr))
    # hourly (dates list should and ends 6-h earlier)
    #datev = dateutils.dateshift(date,6)
    #exptdir = os.path.join(expbasedir+'/'+expt,date)
    #ufsfcst = os.path.join(exptdir,'%s_%s_fhr%02i_ensmean'%(prefix,date,hr))
    print(ufsfcst)
    ncufs = Dataset(ufsfcst)
    ifsanl = os.path.join(ifsanldir,'C%s_ifsanl_%s.nc'% (RES,datev))
    print(ifsanl)
    ncifs = Dataset(ifsanl)
    if not ncount:
        lats = ncufs['lat'][:]
        coslats = np.cos(np.radians(lats[:,0]))
        coslats = coslats[:,np.newaxis]
        ak = ncufs.ak; bk = ncufs.bk # hybrid vertical coordinates
        ncout.ak=ak; ncout.bk=bk
        plevs = 0.01*(ak[::-1] + bk[::-1]*1.e5) # hPa
        plevs_up = plevs[1:]; plevs_dn = plevs[:-1] 
        plevs_mid = 0.5*(plevs[1:]+plevs[:-1]) # mid-levels
        plevsm[:]=plevs_mid
        plevsd[:]=plevs_dn
        plevsu[:]=plevs_up
    plevsml = plevsm[:].tolist()
    plevsml.sort()
    nlevbot = nlevs-bisect(plevsml,levbot)+1
    nlevtop = nlevs-bisect(plevsml,levtop)
           
    tmpdiff = (ncufs['tmp'][:].squeeze()-ncifs['tmp'][:].squeeze())[::-1,...]
    udiff = (ncufs['ugrd'][:].squeeze()-ncifs['ugrd'][:].squeeze())[::-1,...]
    vdiff = (ncufs['vgrd'][:].squeeze()-ncifs['vgrd'][:].squeeze())[::-1,...]
    #for k in range(127):
    #    #print(k,udiff[k].min(),udiff[k].max(),vdiff[k].min(),vdiff[k].max())
    #    print(k,getmean2d(np.sqrt(udiff[k]**2+vdiff[k]**2),coslats[k]))
    #raise SystemExit
    # convert humidity to g/kg
    qdiff = (1000.*(ncufs['spfh'][:].squeeze()-ncifs['spfh'][:].squeeze()))[::-1]
    ncufs.close()

    tmperrsq[ncount] = np.sqrt(getmean(tmpdiff**2,coslats))
    tmpbias[ncount] = getmean(tmpdiff,coslats)
    qerrsq[ncount] = np.sqrt(getmean(qdiff**2,coslats))
    qbias[ncount] = getmean(qdiff,coslats)
    #winderrsq[ncount] = np.sqrt(getmean(udiff**2 + vdiff**2,coslats)
    winderrsq[ncount] = getmean(np.sqrt(udiff**2 + vdiff**2),coslats)
    times[ncount] = dateutils.datetohrs(date)

    print(ncount,date,winderrsq[ncount,nlevbot:nlevtop].mean(),tmperrsq[ncount,nlevbot:nlevtop].mean())
    ncount+=1

ncout.close()
