import os
import sys
import subprocess
import time

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
when_last_locked=-1
# when did we restart last
when_last_stopped=-1
require_gps=False
sandra_control=True
when_last_started=time.time()

while True:

    gps_status = subprocess.run(
        ["/bin/bash","/home/radar/src/git/pansy_receiver/test_gps.sh",],
        capture_output=True,
        text=True,
    )
    gps_status = str(gps_status.stdout)
    #print(gps_status)
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

    if require_gps:
        if is_gps_locked:
            when_last_locked=time.time()
            if is_acq_active:
                print("GPS locked {}, ACQ active {}, ACQ running!".format(is_gps_locked,is_acq_active),flush=True)
            else:
                when_last_started=time.time()
                print("GPS locked {}, ACQ active {}, ACQ start!".format(is_gps_locked,is_acq_active),flush=True)
                if sandra_control:
                    os.system(CMD_SANDRA_START)
        else:
            time_now=time.time()
            holdover_time = time_now-when_last_locked
            time_since_stopped=time_now-when_last_stopped
            # 
            # There is a loss of lock. We stop recording if:
            # 1. loss of gpslock for more than two hours or
            # 2. there has been more than 24 hours since restart
            if is_acq_active and ((holdover_time > 2*3600) or (time_since_stopped > 24*3600)):
                when_last_stopped=time.time()
                print("GPS locked {}, ACQ active {}, ACQ stop!".format(is_gps_locked,is_acq_active),flush=True)
                if sandra_control:
                    os.system(CMD_SANDRA_STOP)
            else:
                print("GPS locked {}, ACQ active {}, ACQ wait!".format(is_gps_locked,is_acq_active),flush=True)
    else:
        # if we don't require GPS, then restart once a day to keep some level of alignment with system clock
        if is_acq_active:
            time_now=time.time()
            print("GPS locked {}, ACQ active {}, ACQ running!".format(is_gps_locked,is_acq_active),flush=True)
            if (time_now-when_last_started)>24*3600:
                os.system(CMD_SANDRA_STOP)
        else:
            when_last_started=time.time()
            print("GPS locked {}, ACQ active {}, ACQ start!".format(is_gps_locked,is_acq_active),flush=True)
            os.system(CMD_SANDRA_START)

    # 
    time.sleep(1.0)
