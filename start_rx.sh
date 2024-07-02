#!/usr/bin/bash
#
# start a ringbuffer
#
# activate conda environment
#source /home/radar/.miniforge3/bin/activate syowa-meteor-radar
#source /home/radar/.miniforge3/bin/activate 

# todo, move this to startup
#sudo sysctl -w net.core.rmem_max=100000000
ENVPYTHON=/home/radar/.miniforge3/envs/syowa-meteor-radar/bin/python3.8
ENVDIR=/home/radar/.miniforge3/envs/syowa-meteor-radar/bin
RUNDIR=/home/radar/src/pansy_receiver
DDIR=/media/buffer/rf

mkdir -p $RUNDIR/logs

# delete old data from ram disk
rm -Rf $DDIR
mkdir -p $DDIR

# add ntp?

# setup ringbuffer
echo "Ringbuffer on RAM disk"
$ENVPYTHON $ENVDIR/drf ringbuffer -z 100000MB $DDIR -p 60 > $RUNDIR/logs/ringbuffer_ramdisk.log 2>&1 &
echo "Mirror - RAM disk to archive"
$ENVPYTHON $ENVDIR/drf mirror cp /media/buffer/rf/ /media/archive/rf/ > $RUNDIR/logs/mirror.log 2>&1 &
echo "Ringbuffer on archive"
$ENVPYTHON $ENVDIR/drf ringbuffer -z 15TB /media/archive/rf/ -p 60 > $RUNDIR/logs/ringbuffer_archive.log 2>&1 &

FREQ=47.5e6
RATE=1e6
#

echo "192.168.11.2"
$ENVPYTHON $ENVDIR/thor.py -m 192.168.11.2 -d "A:A A:B" -c ch0,ch1 -f $FREQ -r $RATE /media/buffer/rf > $RUNDIR/logs/thor.11.2.log 2>&1 &
echo "192.168.12.2"
$ENVPYTHON $ENVDIR/thor.py -m 192.168.12.2 -d "A:A A:B" -c ch2,ch3 -f $FREQ -r $RATE /media/buffer/rf > $RUNDIR/logs/thor.12.2.log 2>&1 &
echo "192.168.12.2"
$ENVPYTHON $ENVDIR/thor.py -m 192.168.13.2 -d "A:A A:B" -c ch4,ch5 -f $FREQ -r $RATE /media/buffer/rf > $RUNDIR/logs/thor.13.2.log 2>&1 &
echo "192.168.12.2"
$ENVPYTHON $ENVDIR/thor.py -m 192.168.14.2 -d "A:A A:B" -c ch6,ch7 -f $FREQ -r $RATE /media/buffer/rf > $RUNDIR/logs/thor.14.2.log 2>&1 &

echo "done"
echo "to stop service, do something"
