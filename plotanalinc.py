from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from dateutils import daterange
expt='C192_hybcov_6hourly_iau'
expt='C192_hybcov_hourly'
date1='2021090500'; date2='2021090500'
date1='2021090421'; date2='2021090503'
dates = daterange(date1,date2,1)
for date in dates:
    filename_fg = '/work2/noaa/gsienkf/whitaker/%s/%s/sfg_%s_fhr01_ensmean' % (expt,date,date)
    filename_an = '/work2/noaa/gsienkf/whitaker/%s/%s/sanl_%s_fhr01_ensmean' % (expt,date,date)
    #filename_fg = '/work2/noaa/gsienkf/whitaker/%s/%s/sfg_%s_fhr06_ensmean' % (expt,date,date)
    #filename_an = '/work2/noaa/gsienkf/whitaker/%s/%s/sanl_%s_fhr06_ensmean' % (expt,date,date)
    
    nc = Dataset(filename_fg)
    lats = nc['grid_yt'][:]
    lons = nc['grid_xt'][:]
    lons, lats = np.meshgrid(lons, lats)
    varname = 'spfh'
    varname = 'tmp'
    varname = 'vgrd'
    #varname = 'pressfc'
    fg = nc[varname][:].squeeze()
    nc.close()
    
    nc = Dataset(filename_an)
    an = nc[varname][:].squeeze()
    nc.close()
    inc = an - fg
    print('RMS',np.sqrt((inc**2).mean()))
    
    print(inc.shape, inc.min(), inc.max())
    clevs = np.arange(-3.0,3.0001,0.25)
    #clevs = np.arange(-0.25,0.251,0.0251)
    #clevs = np.arange(-7.5e-4,7.5001e-4,0.75e-4)
    if varname == 'pressfc':
        inc = 0.01*inc
        print(inc.min(), inc.max())
        incplot = inc
        nlev = 0
    else:
        nlev = 68
        for k in range(inc.shape[0]):
            print(k,inc[k].min(),inc[k].max())
        incplot = inc[nlev]
    
    plt.figure()
    plt.contourf(lons,lats,incplot,clevs,cmap=plt.cm.bwr,extend='both')
    plt.colorbar(shrink=0.6)
    plt.title('%s level %s increment %s' % (varname, nlev, date))
    plt.gca().set_aspect('equal')
    plt.savefig('%s_level%s_inc1_%s.png' % (varname, nlev, date))
#plt.show()
