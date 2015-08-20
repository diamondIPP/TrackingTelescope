#!/bin/sh

if [ "$1" = "" ]; then
	echo "Forget to specify root file! --> exiting"
	exit
fi

if [ "$2" != "" ]; then
	mode=$2
else
	mode=0
fi

if [ "$3" != "" ]; then
    telescope=$3
else
    telescope=9
fi

file=`pwd`"/"$1

cd ~/software/trackingTelescope
./TrackingTelescope $file $mode $telescope
cd -

