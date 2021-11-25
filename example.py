from TinySatTracker import SatTracker
from datetime import datetime
from dateutil import rrule

# NOAA-18 NORAD ID
norad_id = 28654
tr = SatTracker(norad_id)

# MIPT LK geodetic coordinates
lk_lat = 55.930195540685304
lk_lon = 37.5182925079981
lk_alt = 193

# Get azimuth and elevation angles at specific date
az, al = tr.getSatLookAngles([lk_lat, lk_lon, lk_alt], [2021, 11, 25, 14, 5, 0])
print(az, al)

start_date = datetime(2021, 11, 25, 13, 30, 0)
end_date = datetime(2021, 11, 25, 16, 0, 0)

# Find and plot satellite passes in a given time range
intercepts = tr.findIntercepts([lk_lat, lk_lon, lk_alt], start_date, end_date,
        rule=rrule.SECONDLY, plot=True)

# Max elevation points are also stored
for inter in intercepts:
    print(inter['max_alt'])
