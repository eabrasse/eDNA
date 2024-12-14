"""
generate bash file of tracker commands to perform hourly particle releases
"""
from datetime import datetime, timedelta
import pytz
import numpy as np

Pacific = pytz.timezone("PST8PDT")
utc = pytz.utc


ytd_list = ['2024.06.11','2024.06.12','2024.06.13','2024.06.14']
dtt_dict = {'2024.06.11':4,'2024.06.12':3,'2024.06.13':2,'2024.06.14':1}

filecount=0
runcount=0
fname = 'hcdolphJune2024_{}.sh'.format(str(filecount).zfill(2))
f = open(fname, 'w')
f.write('#!/bin/bash\ncd /data2/pmr4/eab32/LO/tracker2/\n')

for ytd in ytd_list:
    dtt = dtt_dict[ytd]

    for sh in range(24):
        if np.mod(runcount,10)==0: #only ten runs per file. Start a new file after 10.
            if runcount>0:
                f.close()
                filecount+=1
            
                fname = 'hcdolphJune2024_{}.sh'.format(str(filecount).zfill(2))
                f = open(fname, 'w')
                f.write('#!/bin/bash\ncd /data2/pmr4/eab32/LO/tracker2/\n')
        
        mydt = utc.localize(datetime(int(ytd[:4]),int(ytd[5:7]),int(ytd[8:10]),sh))
        logstring = str(runcount).zfill(2)
        my_str = f'python tracker.py -gtx hc11_v01_uu0k -exp hc_dolph -3d True -dtt {dtt} -clb True -d {ytd} -sh {str(sh)} -sub_tag {ytd}_sh{str(sh)} > /data2/pmr4/eab32/LO_output/track_logs/hcdolphJune2024_log{logstring}.txt &\n'
        f.write(my_str)
        runcount+=1

f.close()