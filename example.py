from TinySatTracker import SatTracker
from datetime import datetime
from dateutil import rrule

norad_id = 28654
tr = SatTracker(norad_id)

lk_lat = 55.930195540685304
lk_lon = 37.5182925079981
lk_alt = 193

az, al = tr.getSatLookAngles([lk_lat, lk_lon, lk_alt], [2021, 11, 25, 14, 5, 0])
print(az, al)

start_date = datetime(2021, 11, 25, 13, 30, 0)
end_date = datetime(2021, 11, 25, 16, 0, 0)

intercepts = tr.findIntercepts([lk_lat, lk_lon, lk_alt], start_date, end_date,
        rule=rrule.SECONDLY, plot=True)

for inter in intercepts:
    print(inter['max_alt'])
