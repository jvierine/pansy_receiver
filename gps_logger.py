#!/usr/bin/env python3

import os
import time

tnow=time.time()

logfname="logs/gpslog-%f.log"%(tnow)

os.environ['LD_LIBRARY_PATH']="/home/radar/.syowa/lib"
lfo=open(logfname,"w")
while True:
    os.system("rm -f /tmp/gpslock")
    cmd='./test_gps 2>/dev/null |grep "GPS lock" 1> /tmp/gpslock'
    os.system(cmd)
    f=open("/tmp/gpslock")
    l=f.readlines()[0].strip()
    tnow=time.time()
    lfo.write("%f %s\n"%(tnow,l))
    lfo.flush()
    time.sleep(1)

lfo.close()
