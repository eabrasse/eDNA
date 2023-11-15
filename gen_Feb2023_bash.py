"""
generate bash file of tracker commands to perform hourly particle releases
"""
from datetime import datetime, timedelta
import pytz
import numpy as np

Pacific = pytz.timezone("PST8PDT")
utc = pytz.utc


ytd_list = ['2023.01.31','2023.02.01','2023.02.02','2023.02.03']
dtt_dict = {'2023.01.31':4,'2023.02.01':3,'2023.02.02':2,'2023.02.03':1}

filecount=0
runcount=0
fname = 'hcdolph_3d_sampling_{}.sh'.format(str(filecount).zfill(2))
f = open(fname, 'w')
f.write('#!/bin/bash\ncd /data2/pmr4/eab32/LO/tracker/\n')

for ytd in ytd_list:
    dtt = dtt_dict[ytd]

    for sh in range(24):
        if np.mod(runcount,10)==0: #only ten runs per file. Start a new file after 10.
            if runcount>0:
                f.close()
                filecount+=1
            
                fname = 'hcdolph_3d_sampling_{}.sh'.format(str(filecount).zfill(2))
                f = open(fname, 'w')
                f.write('#!/bin/bash\ncd /data2/pmr4/eab32/LO/tracker/\n')
        
        mydt = utc.localize(datetime(int(ytd[:4]),int(ytd[5:7]),int(ytd[8:10]),sh))
        logstring = str(runcount).zfill(2)
        my_str = f'python tracker.py -gtx hc11_v01_uu0k -exp hc_dolph -3d True -dtt {dtt} -clb True -d {ytd} -sh {str(sh)} -sub_tag {ytd} > /data2/pmr4/eab32/LO_output/track_logs/log{logstring}.txt &\n'
        f.write(my_str)
        runcount+=1

f.close()