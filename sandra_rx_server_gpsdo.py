import os
import sys
import subprocess
import time

import digital_rf as drf
import pansy_config as pc
subdirectory_cadence_seconds = 3600
file_cadence_seconds = 60
samples_per_second_numerator = 1000000
samples_per_second_denominator = 1
file_name = "lock"
omddir=pc.gpslock_metadata_dir
os.system("mkdir -p %s"%(omddir))

dmw = drf.DigitalMetadataWriter(
    omddir,
    subdirectory_cadence_seconds,
    file_cadence_seconds,
    samples_per_second_numerator,
    samples_per_second_denominator,
    file_name,
)

GPS_LOCKED_SUBSTRINGS = [
    "GPS lock status: locked",
    #"GPSDO detected: true",
]

PATH_PYTHON_INTERPRETER = "/home/radar/.sandra/envs/sandra-2.7/bin/python2.7"
PATH_SANDRA_CLIENT = "/home/radar/.sandra/envs/sandra-2.7/bin/sandra-rx-client.py"

CMD_SANDRA  = PATH_PYTHON_INTERPRETER+" "
CMD_SANDRA += PATH_SANDRA_CLIENT

CMD_SANDRA_STOP = CMD_SANDRA+" stop"
CMD_SANDRA_START = CMD_SANDRA+" start"

# when was gps lock last seen.
require_gps=True
when_last_started=time.time()
when_last_locked=time.time()
holdover=0

while True:

    gps_status = subprocess.run(
        ["/bin/bash","/home/radar/src/git/pansy_receiver/test_gps.sh",],
        capture_output=True,
        text=True,
    )
    gps_status = str(gps_status.stdout)
    
    is_gps_locked = all([
        GPS_LOCKED_SUBSTRING in gps_status
        for GPS_LOCKED_SUBSTRING in GPS_LOCKED_SUBSTRINGS
    ])

    acq_status = subprocess.run(
        ["tail","-1","/home/radar/.sandra/envs/sandra-2.7/var/log/rx.events",],
        capture_output=True,
        text=True,
    )
    acq_status = str(acq_status.stdout)
    acq_status = acq_status.split(" ")[-1]
    acq_status = acq_status.split("\n")[0]

    if acq_status=="stop":
        is_acq_active = False
    else:
        is_acq_active = True

    time_now=time.time()
    if require_gps:
        if is_gps_locked:
            when_last_locked=time_now
            if is_acq_active:
                print("GPS locked {}, ACQ active {}, ACQ running!".format(is_gps_locked,is_acq_active),flush=True)
                if (time_now-when_last_started)>24*3600:
                    # if over 24 hours since restart and we have a lock, then restart now when most likely there is 
                    # not going to be a wait
                    os.system(CMD_SANDRA_STOP)

            else:
                print("GPS locked {}, ACQ active {}, ACQ start!".format(is_gps_locked,is_acq_active),flush=True)
                when_last_started = time_now
                os.system(CMD_SANDRA_START)
        else:
            if is_acq_active:
                # switch to no gps mode if we have lost the lock for more than 24 hours
                holdover = time_now - when_last_locked
                if holdover > 24*3600:
                    print("GPS locked {}, ACQ active {}, ACQ stop!".format(is_gps_locked,is_acq_active),flush=True)
                    os.system(CMD_SANDRA_STOP)
                else:
                    print("GPS locked {}, ACQ active {}, holdover {} s".format(is_gps_locked,is_acq_active,holdover),flush=True)
            else:
                print("GPS locked {}, ACQ active {}, ACQ wait!".format(is_gps_locked,is_acq_active),flush=True)
        data_dict={"holdover":holdover,"acq":is_acq_active,"lock":is_gps_locked}
        dmw.write(int(time_now*1000000),data_dict)
    else:
        # if we don't require GPS, then restart once a day to keep some level of alignment with system clock
        if is_acq_active:
            print("GPS locked {}, ACQ active {}, ACQ running!".format(is_gps_locked,is_acq_active),flush=True)
            if (time_now-when_last_started)>24*3600:
                os.system(CMD_SANDRA_STOP)
        else:
            when_last_started=time_now
            print("GPS locked {}, ACQ active {}, ACQ start!".format(is_gps_locked,is_acq_active),flush=True)
            os.system(CMD_SANDRA_START)

    time.sleep(1.0)
