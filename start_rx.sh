#!/usr/bin/bash
#
# start a ringbuffer
#
# activate conda environment
source /home/radar/.miniforge3/bin/activate 

# todo, move this to startup
sudo sysctl -w net.core.rmem_max=100000000

RUNDIR=/home/radar/src/pansy_receiver/
DDIR=/media/buffer/rf

mkdir -p $RUNDIR/logs

# delete old data from ram disk
rm -Rf $DDIR
mkdir -p $DDIR

# add ntp?

# setup ringbuffer
echo "Ringbuffer"
drf ringbuffer -z 127000MB $DDIR -p 2 > $RUNDIR/logs/ringbuffer.log 2>&1 &

FREQ=47.5e6
RATE=1e6
#
#thor.py -m 192.168.11.2 -d "A:A A:B" -c ch0,ch1 -f 47.5e6 -r 1e6 /media/buffer/rf
echo "192.168.12.2"
thor.py -m 192.168.12.2 -d "A:A A:B" -c ch2,ch3 -f $FREQ -r $RATE /media/buffer/rf > $RUNDIR/logs/thor.12.2.log 2>&1 &
echo "192.168.12.2"
thor.py -m 192.168.13.2 -d "A:A A:B" -c ch4,ch5 -f $FREQ -r $RATE /media/buffer/rf > $RUNDIR/logs/thor.13.2.log 2>&1 &
echo "192.168.12.2"
thor.py -m 192.168.14.2 -d "A:A A:B" -c ch6,ch7 -f $FREQ -r $RATE /media/buffer/rf > $RUNDIR/logs/thor.14.2.log 2>&1 &

echo "done"
echo "to stop service, do something"
