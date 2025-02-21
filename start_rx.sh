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
nohup $ENVPYTHON $ENVDIR/drf ringbuffer -z 30GB /media/buffer -p 10 > $RUNDIR/logs/ringbuffer_ramdisk.log 2>&1 &
echo "Mirror - RAM disk to archive"
nohup $ENVPYTHON $ENVDIR/drf mirror cp /media/buffer/ /media/archive/ > $RUNDIR/logs/mirror.log 2>&1 &
echo "Ringbuffer on archive"
nohup $ENVPYTHON $ENVDIR/drf ringbuffer -z 30000GB /media/archive/ -p 60 > $RUNDIR/logs/ringbuffer_archive.log 2>&1 &

echo "done"
echo "to stop service, do something"
