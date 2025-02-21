#!/usr/bin/bash
ENVPYTHON=/home/radar/.sandra/envs/sandra-2.7/bin/python2.7
ENVDIR=/home/radar/.sandra/envs/sandra-2.7/bin
RUNDIR=/home/radar/src/git/pansy_receiver
DDIR=/media/buffer/

mkdir -p $RUNDIR/logs

FREQ=47.0e6
RATE=1000000

$ENVPYTHON $ENVDIR/thor.py -m 192.168.11.2 -d "A:A A:B" -c ch000,ch001 -f $FREQ -r $RATE --clock_source external --time_source external /media/buffer > $RUNDIR/logs/thor.11.2.log 2>&1 &
$ENVPYTHON $ENVDIR/thor.py -m 192.168.12.2 -d "A:A A:B" -c ch002,ch003 -f $FREQ -r $RATE --clock_source external --time_source external  /media/buffer > $RUNDIR/logs/thor.12.2.log 2>&1 &
$ENVPYTHON $ENVDIR/thor.py -m 192.168.13.2 -d "A:A A:B" -c ch004,ch005 -f $FREQ -r $RATE --clock_source external --time_source external  /media/buffer > $RUNDIR/logs/thor.13.2.log 2>&1 &
$ENVPYTHON $ENVDIR/thor.py -m 192.168.14.2 -d "A:A A:B" -c ch006,ch007 -f $FREQ -r $RATE --clock_source external --time_source external  /media/buffer > $RUNDIR/logs/thor.14.2.log 2>&1 &
