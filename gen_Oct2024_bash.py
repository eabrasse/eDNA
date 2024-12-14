"""
generate bash file of tracker commands to perform hourly particle releases
"""
from datetime import datetime, timedelta
import pytz
import numpy as np

Pacific = pytz.timezone("PST8PDT")
utc = pytz.utc

filecount=0
runcount=0
fname = 'hcdolphOct16.2024.sh'
f = open(fname, 'w')
f.write('#!/bin/bash\ncd /data2/pmr4/eab32/LO/tracker2/\n')



for sh in range(19,24):

    logstring = str(runcount).zfill(2)
    my_str = f'python tracker.py -gtx hc11_v01_uu0k -exp hc_dolph2 -3d True -dtt 1 -clb True -d 2024.10.16 -sh {str(sh)} -sph 12 -sub_tag 2024.10.16_sh{str(sh)} > /data2/pmr4/eab32/LO_output/track_logs/hcdolphOct2024_log{logstring}.txt &\n'
    f.write(my_str)
    runcount+=1

for sh in range(0,4):
    logstring = str(runcount).zfill(2)
    my_str = f'python tracker.py -gtx hc11_v01_uu0k -exp hc_dolph2 -3d True -dtt 1 -clb True -d 2024.10.17 -sh {str(sh)} -sph 12 -sub_tag 2024.10.17_sh{str(sh)} > /data2/pmr4/eab32/LO_output/track_logs/hcdolphOct2024_log{logstring}.txt &\n'
    f.write(my_str)
    runcount+=1

f.close()