#!/usr/bin/env python3

from os.path import expanduser
import netifaces
import sys
import subprocess


def get_saved_interface_name():
    home = expanduser("~")
    name = ""
    try:
        with open(home + "/.cheetah_network.txt"):
            name = f.read().split()[0]
    except:
        name = ""
    return name


def get_likely_iface():
    ifs = netifaces.interfaces()
    print("Found {} interfaces:".format(len(ifs)))

    if_to_addrs = {}

    for i in ifs:
        if_to_addrs[i] = []
        if netifaces.AF_INET in netifaces.ifaddresses(i).keys():
            for ad in netifaces.ifaddresses(i)[netifaces.AF_INET]:
                if_to_addrs[i].append(ad['addr'])


    for i in range(len(ifs)):
        print("  [{}] : {} : {}".format(i, ifs[i], if_to_addrs[ifs[i]]))

    found_10_ip = 0
    selected_if = ""


    for i in ifs:
        match_string = "10.0.0."
        if len(if_to_addrs[i]) > 0 and if_to_addrs[i][0][:len(match_string)] == match_string:
            found_10_ip = found_10_ip + 1
            selected_if = i

    if found_10_ip == 0:
        print("None of the network adapters look correct.  Make sure you have set a 10.0.0.x static ip!")
        return ""

    elif found_10_ip == 1:
        print("The adapter {} seems correct".format(selected_if))
        return selected_if

    else:
        print("Found {} possible adapters, giving up".format(found_10_ip))
        return ""
    


def main():
    name = get_saved_interface_name()
    if not name:
        print("Didn't find saved interface, searching...")
        name = get_likely_iface()
        if not name:
            sys.exit("Failed to find network adapter name")
    else:
        print("found saved interface {}".format(name))

    print("setup for interface {}".format(name))
    subprocess.call(['sudo', 'ifconfig', name, 'multicast'])
    subprocess.call(['sudo', 'route', 'add', '-net', '224.0.0.0', 'netmask', '240.0.0.0', 'dev', name])


    return 0


main()
