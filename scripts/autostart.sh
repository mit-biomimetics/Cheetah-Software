#/bin/bash

# this script should be copied to /home/user/ on mini cheetahs for autostarting.
# comment out the following line to enable auto-start

#exit 1

cd /home/user/robot-software/build
./run_mc_2_autostart.sh ./mit_ctrl auto-start
