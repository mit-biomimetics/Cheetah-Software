#!/bin/bash
if [ $# -eq 0 ]; then
    echo "Usage: config_network_lcm.sh -I [interface]"
    echo "or config_network_lcm.sh [computer]"
    echo "interface: network interface to configure"
    echo "computer: use stored interface for computer"
    echo "current computers:"
    echo " name       interface    description"
    echo " ----       --------     -----------"
    echo " thinkpad   enp0s25      "
    echo " cynergy    enxa0cec808fb18"
    echo " xps-wifi   wlp2s0"
    echo " mc         enp1s0"
    echo " mc-usb     enxa0cec80424d3"
    echo " mc-top     enxa0cec80e3ced"
fi

if [ "$1" == "-I" ]; then
    sudo ifconfig $2 multicast
    sudo route add -net 224.0.0.0 netmask 240.0.0.0 dev $2
    #sudo route add -net 10.0.0.0 netmask 255.0.0.0 dev $2
fi

if [ "$1" == "thinkpad" ]; then
    sudo ifconfig enp0s25 multicast
    sudo route add -net 224.0.0.0 netmask 240.0.0.0 dev enp0s25
fi

if [ "$1" == "cynergy" ]; then
    sudo ifconfig eenxa0cec808fb18 multicast
    sudo route add -net 224.0.0.0 netmask 240.0.0.0 dev enxa0cec808fb18
fi

if [ "$1" == "xps-wifi" ]; then
    sudo ifconfig wlp2s0 multicast
    sudo route add -net 224.0.0.0 netmask 240.0.0.0 dev wlp2s0
fi

if [ "$1" == "mc" ]; then
    sudo ifconfig enp1s0 multicast
    sudo route add -net 224.0.0.0 netmask 240.0.0.0 dev enp1s0
fi

if [ "$1" == "mc-usb" ]; then
    sudo ifconfig enxa0cec80424d3 multicast
    sudo route add -net 224.0.0.0 netmask 240.0.0.0 dev enxa0cec80424d3
fi

if [ "$1" == "mc-top" ]; then
    sudo ifconfig enxa0cec80e3ced multicast
    sudo route add -net 224.0.0.0 netmask 240.0.0.0 dev enxa0cec80e3ced
fi

if [ "$1" == "dhkim" ]; then
    sudo ifconfig enx70886b885732 multicast
    sudo route add -net 224.0.0.0 netmask 240.0.0.0 dev enx70886b885732
fi


if [ "$1" == "ben-vm" ]; then
    sudo ifconfig enxa0cec80424d3 multicast
    sudo route add -net 224.0.0.0 netmask 240.0.0.0 dev enxa0cec80e3ced
fi

if [ "$1" == "dhmac" ]; then
    sudo ifconfig enx70886b887f40 multicast
    sudo route add -net 224.0.0.0 netmask 240.0.0.0 dev enx70886b887f40
fi
