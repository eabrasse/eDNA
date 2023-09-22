"""
generate bash file of tracker commands to perform hourly particle releases
"""
from datetime import datetime, timedelta
import pytz

Pacific = pytz.timezone("PST8PDT")
utc = pytz.utc

first_sample_utc = Pacific.localize(datetime(2023,6,5,10,30)).astimezone(utc)
first_release = first_sample_utc-timedelta(days=2)
last_sample_utc = Pacific.localize(datetime(2023,6,5,15,30)).astimezone(utc)

fname = 'June5_hc_bangor_perimeter_point_rels.sh'
f = open(fname, 'w')

f.write('#!/bin/bash\ncd /data2/pmr4/eab32/LO/tracker/\n')

logcount=0
for i in range(first_release.day,last_sample_utc.day+1):
    ytd = '2023.06.0{}'.format(i)
    for sh in range(24):
        mydt = utc.localize(datetime(2023,6,i,sh))
        if mydt<last_sample_utc:
            if mydt>first_release:
                logstring = str(logcount).zfill(2)
                my_str = f'python tracker.py -gtx hc11_v01_uu0k -exp hc_bangor_perimeter_point -3d True -dtt 2 -clb True -d {ytd} -sh {str(sh)} -sub_tag {ytd} -sph 4 > /data2/pmr4/eab32/LO_output/track_logs/log{logstring}.txt &\n'
                f.write(my_str)
                logcount+=1
        else:
            break
f.close()