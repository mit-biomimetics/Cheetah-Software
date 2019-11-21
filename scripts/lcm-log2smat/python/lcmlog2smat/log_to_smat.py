#!/usr/bin/python
#
# Converts a LCM log to a "structured" format that is easier to work with
# external tools such as Matlab or Python. The set of messages on a given
# channel can be represented as a structure preserving the original lcm message
# structure.
# 
# lcm-log2smat is based on libbot2 script bot-log2mat.
# Modified by G.Troni

import os
import sys
import binascii
import types
import numpy
import re
import getopt

# Python 2.7 and 3 
if sys.version_info >= (3, 0):
    import pickle
    xrange  = range
    long    = int
    unicode = int
else:
    import cPickle as pickle

# check which version for mio location
if sys.version_info < (2, 6):
    import scipy.io.mio
else:
    import scipy.io.matlab.mio



from lcm import EventLog
from .scan_for_lcmtypes import *

longOpts = ["help", "print", "pickle", "format", "separator", "channelsToProcess", "ignore", "outfile", "lcm_packages"]

def usage():
    pname, sname = os.path.split(sys.argv[0])
    sys.stderr.write("usage: % s %s < filename > \n" % (sname, str(longOpts)))
    print( """
    -h --help                 print this message
    -p --print                Output log data to stdout instead of to .mat or .pkl
    -k --pickle               Output log data to python pickle .pkl instead of to .mat
    -f --format               print the data format to stderr
    -s --seperator=sep        print data with separator [sep] instead of default to ["" ""]
    -c --channelsToProcess=chan        Parse channelsToProcess that match Python regex [chan] defaults to [".*"]
    -i --ignore=chan          Ignore channelsToProcess that match Python regex [chan]
                              ignores take precedence over includes!
    -o --outfile=ofname       output data to [ofname] instead of default [filename.mat or filename.pkd or stdout]
    -l --lcmtype_pkgs=pkgs    load python modules from comma seperated list of packages [pkgs] defaults to ["botlcm"]
    -v                        Verbose

    """ )
    sys.exit()


def msg_getfields (lcm_msg):
    return lcm_msg.__slots__


def msg_getconstants (lcm_msg):
    # Get full list of valid attributes
    fulllist = dir(lcm_msg)
    # Get constants
    constantslist = [x for x in fulllist if not(x[0]=='_') 
                    if not(x=='decode')
                    if not(x=='encode')
                    if x not in msg_getfields (lcm_msg)]
    return constantslist


def msg_to_dict (data, e_channel, msg, statusMsg, verbose=False, lcm_timestamp=-1):

    # Initializing channel
    if e_channel not in data:
        data[e_channel] = dict()

        # Iterate each constant of the LCM message
        constants = msg_getconstants (msg)
        for i in xrange(len(constants)):
            myValue = None
            myValue = eval('msg.' + constants[i])
            data[e_channel][constants[i][:31]] = myValue

    # Get lcm fields and constants
    fields = msg_getfields (msg)

    # Iterate each field of the LCM message
    for i in xrange(len(fields)):
        myValue = None
        myValue = eval(' msg.' + fields[i])
        if (isinstance(myValue,int)     or
            isinstance(myValue,long)    or 
            isinstance(myValue,float)   or 
            isinstance(myValue,tuple)   or 
            isinstance(myValue,unicode) or 
            isinstance(myValue,str)):
            try:
                data[e_channel][fields[i][:31]].append(myValue)
            except KeyError as AttributeError:
                data[e_channel][fields[i][:31]] = [(myValue)]

        elif (hasattr(myValue,'__slots__')):
            submsg = eval('msg.' + fields[i])
            msg_to_dict (data[e_channel], fields[i][:31], submsg, statusMsg, verbose)

        else:
            if verbose:
                statusMsg = deleteStatusMsg(statusMsg)
                sys.stderr.write("ignoring field %s from channel %s. \n" %(fields[i], e_channel))
            continue

    # Add extra field with lcm_timestamp
    if lcm_timestamp > 0:
        try:
            data[e_channel]['lcm_timestamp'].append(lcm_timestamp)
        except KeyError as AttributeError:
            data[e_channel]['lcm_timestamp'] = [(lcm_timestamp)]


def deleteStatusMsg(statMsg):
    if statMsg:
        sys.stderr.write("\r")
        sys.stderr.write(" " * (len(statMsg)))
        sys.stderr.write("\r")
    return ""



def parse_and_save (args, opts={}):

    if isinstance(args, (list, tuple)):
        fname = args[0]
    elif isinstance(args, str):
        fname = args
    else:
        usage()


    #default options
    lcm_packages = [ "botlcm"]

    printFname = "stdout"
    printFile = sys.stdout
    verbose = False
    printOutput = False
    savePickle  = False
    printFormat = False
    channelsToIgnore = ""
    checkIgnore = False
    channelsToProcess = ".*"
    separator = ' '
    for o, a in opts:
        if o == "-v":
            verbose = True
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-p", "--print"):
            printOutput = True
        elif o in ("-k", "--pickle"):
            savePickle = True
        elif o in ("-f", "--format"):
            printFormat = True
        elif o in ("-s", "--separator="):
            separator = a
        elif o in ("-o", "--outfile="):
            outFname = a
            printFname = a
        elif o in ("-c", "--channelsToProcess="):
            channelsToProcess = a
        elif o in ("-i", "--ignore="):
            channelsToIgnore = a
            checkIgnore = True
        elif o in ("-l", "--lcm_packages="):
            lcm_packages = a.split(",")
        else:
            assert False, "unhandled option"

    try:
        outFname
    except NameError:
        outDir, outFname = os.path.split(os.path.abspath(fname))
       
        if savePickle:
            outFname = outDir + "/" + outFname + ".pkl"
        else:
            outFname = outFname.replace(".", "_")
            outFname = outFname.replace("-", "_")
            outFname = outDir + "/" + outFname + ".mat"



    fullPathName = os.path.abspath(outFname)
    dirname = os.path.dirname(fullPathName)
    outBaseName = ".".join(os.path.basename(outFname).split(".")[0:-1])
    fullBaseName = dirname + "/" + outBaseName

    data = {}

    print("Searching for LCM types...")
    type_db = make_lcmtype_dictionary()

    channelsToProcess = re.compile(channelsToProcess)
    channelsToIgnore = re.compile(channelsToIgnore)
    log = EventLog(fname, "r")

    if printOutput:
        sys.stderr.write("opened % s, printing output to %s \n" % (fname, printFname))
        if printFname == "stdout":
            printFile = sys.stdout
        else:
            printFile = open(printFname, "w")
    else:
        sys.stderr.write("opened % s, outputing to % s\n" % (fname, outFname))

    ignored_channels = []
    msgCount = 0
    statusMsg = ""
    startTime = 0


    # Iterate LCM log file
    for e in log:
        if msgCount == 0: 
            startTime = e.timestamp
            
        if e.channel in ignored_channels:
            continue
        if ((checkIgnore and channelsToIgnore.match(e.channel) and len(channelsToIgnore.match(e.channel).group())==len(e.channel)) \
             or (not channelsToProcess.match(e.channel))):
            if verbose:
                statusMsg = deleteStatusMsg(statusMsg)
                sys.stderr.write("ignoring channel %s\n" % e.channel)
            ignored_channels.append(e.channel)
            continue

        packed_fingerprint = e.data[:8]
        lcmtype = type_db.get(packed_fingerprint, None)
        if not lcmtype:
            if verbose:
                statusMsg = deleteStatusMsg(statusMsg)
                sys.stderr.write("ignoring channel %s -not a known LCM type\n" % e.channel)
            ignored_channels.append(e.channel)
            continue
        try:
            msg = lcmtype.decode(e.data)
        except:
            statusMsg = deleteStatusMsg(statusMsg)
            sys.stderr.write("error: couldn't decode msg on channel %s\n" % e.channel)
            continue
        
        msgCount = msgCount + 1
        if (msgCount % 5000) == 0:
            statusMsg = deleteStatusMsg(statusMsg)
            statusMsg = "read % d messages, % d %% done" % (msgCount, log.tell() / float(log.size())*100)
            sys.stderr.write(statusMsg)
            sys.stderr.flush()
        
        msg_to_dict (data, e.channel, msg, statusMsg, verbose, (e.timestamp - startTime) / 1e6)



    deleteStatusMsg(statusMsg)
    if not printOutput:
                
        sys.stderr.write("loaded all %d messages, saving to % s\n" % (msgCount, outFname))

        if savePickle:

            # Pickle the list/dictonary using the highest protocol available.
            output = open(outFname, 'wb')
            pickle.dump(data, output, -1)
            output.close()
        else:
            # Matlab format using scipy
            if sys.version_info < (2, 6):
                scipy.io.mio.savemat(outFname, data)
            else:
                scipy.io.matlab.mio.savemat(outFname, data, oned_as='row')



if __name__ == "__main__":
    # Parse command line arguments
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "hpktvfs:c:i:o:l:", longOpts)
    except getopt.GetoptError as err:
        # print help information and exit:
        print (str(err)) # will print something like "option -a not recognized"
        usage()
    if len(args) != 1:
        usage()

    # Run main parser
    parse_and_save (args,opts)
