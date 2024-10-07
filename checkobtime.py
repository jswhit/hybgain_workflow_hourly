from netCDF4 import Dataset
import sys
diagfile = sys.argv[1]
nc = Dataset(diagfile)
try:
    time = nc['Time'][:]
except:
    time = nc['Obs_Time'][:]
print(time.min(), time.max())
