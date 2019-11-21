#!/bin/bash

if [ "$(whoami)" != "root" ]; then
  echo -e "This is the auto start script.  Don't run this."
  exit -1
fi

if [ "$2" != "auto-start" ]; then
  echo -e "This is the auto start script. Don't run this!"
  exit -1
fi

# enable multicast and add route for lcm out the top
#sudo ifconfig enxa0cec80e3ced multicast
#sudo route add -net 224.0.0.0 netmask 240.0.0.0 dev enxa0cec80e3ced
python ./setup_network_mc.py

# configure libraries
sudo LD_LIBRARY_PATH=. ldconfig
#sudo LD_LIBRARY_PATH=. ldd ./robot
sudo LD_LIBRARY_PATH=. $1 m r f
