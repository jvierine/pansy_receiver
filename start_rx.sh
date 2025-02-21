#!/usr/bin/bash
#
# start a ringbuffer
#
# activate conda environment
#source /home/radar/.miniforge3/bin/activate syowa-meteor-radar
#source /home/radar/.miniforge3/bin/activate 

# todo, move this to startup
#sudo sysctl -w net.core.rmem_max=100000000
ENVPYTHON=/home/radar/.sandra/envs/sandra-2.7/bin/python2.7
#/home/radar/.miniforge3/envs/syowa-meteor-radar/bin/python3.8
ENVDIR=/home/radar/.sandra/envs/sandra-2.7/bin
#/home/radar/.miniforge3/envs/syowa-meteor-radar/bin
RUNDIR=/home/radar/src/git/pansy_receiver
DDIR=/media/buffer/

mkdir -p $RUNDIR/logs

# delete old data from ram disk
#rm -Rf $DDIR
#mkdir -p $DDIR

# add ntp?

# setup ringbuffer
echo "Ringbuffer on RAM disk"
$ENVPYTHON $ENVDIR/drf ringbuffer -z 30GB /media/buffer -p 10 > $RUNDIR/logs/ringbuffer_ramdisk.log 2>&1 &
echo "Mirror - RAM disk to archive"
$ENVPYTHON $ENVDIR/drf mirror cp /media/buffer/ /media/archive/ > $RUNDIR/logs/mirror.log 2>&1 &
echo "Ringbuffer on archive"
$ENVPYTHON $ENVDIR/drf ringbuffer -z 30000GB /media/archive/ -p 60 > $RUNDIR/logs/ringbuffer_archive.log 2>&1 &

FREQ=47.0e6
RATE=1000000
#

echo "192.168.11.2"
$ENVPYTHON $ENVDIR/thor.py -m 192.168.11.2 -d "A:A A:B" -c ch000,ch001 -f $FREQ -r $RATE --clock_source external --time_source external /media/buffer > $RUNDIR/logs/thor.11.2.log 2>&1 &
echo "192.168.12.2"
$ENVPYTHON $ENVDIR/thor.py -m 192.168.12.2 -d "A:A A:B" -c ch002,ch003 -f $FREQ -r $RATE --clock_source external --time_source external  /media/buffer > $RUNDIR/logs/thor.12.2.log 2>&1 &
echo "192.168.12.2"
$ENVPYTHON $ENVDIR/thor.py -m 192.168.13.2 -d "A:A A:B" -c ch004,ch005 -f $FREQ -r $RATE --clock_source external --time_source external  /media/buffer > $RUNDIR/logs/thor.13.2.log 2>&1 &
echo "192.168.12.2"
$ENVPYTHON $ENVDIR/thor.py -m 192.168.14.2 -d "A:A A:B" -c ch006,ch007 -f $FREQ -r $RATE --clock_source external --time_source external  /media/buffer > $RUNDIR/logs/thor.14.2.log 2>&1 &

echo "done"
echo "to stop service, do something"
