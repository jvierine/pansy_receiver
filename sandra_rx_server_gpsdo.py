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

    if is_gps_locked:
        if is_acq_active:
            print("GPS locked {}, ACQ active {}, ACQ running!".format(is_gps_locked,is_acq_active),flush=True)
        else:
            print("GPS locked {}, ACQ active {}, ACQ start!".format(is_gps_locked,is_acq_active),flush=True)
            os.system(CMD_SANDRA_START)
    else:
        if is_acq_active:
            print("GPS locked {}, ACQ active {}, ACQ stop!".format(is_gps_locked,is_acq_active),flush=True)
            os.system(CMD_SANDRA_STOP)
        else:
            print("GPS locked {}, ACQ active {}, ACQ wait!".format(is_gps_locked,is_acq_active),flush=True)

    time.sleep(1.0)
