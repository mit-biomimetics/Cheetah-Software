##lcm-log2smat
-----

  `lcm-log2smat` is a python utility to convert Lightweight Communications and
  Marshalling (LCM) logfiles to Matlab (.mat) or Python Pickle (.pkl) files.

  The script converts a LCM logfile to a "structured" format that is easier to
  work with external tools such as Matlab or Python. The set of messages on a
  given channel is saved as a structure preserving the original LCM message
  structure. It works with nested messages, it is flexible and independent to
  the LCM message type definitions and it can be used as a script or imported
  as a Python module.


###Requirements
-----

  `lcm-log2smat` has been tested only on the GNU/Linux operating systems.
  But it should work in other operating systems running Python.

  The requirements are:
  * [Python](https://www.python.org/)
  * [LCM](https://lcm-proj.github.io/)  
  * [CMake](http://www.cmake.org/) (only required for installation) 


###Install
-----

  `lcm-log2smat` adheres to the Pods core policy, and makes use of the Pods
  build tools. For more information about Pods, see:
  https://sourceforge.net/p/pods/

  You can install `lcm-log2smat` to your system by running:

    $ sudo make BUILD_PREFIX=/usr/local

  For a local installation, you can build `lcm-log2smat` inside the source
  directory. To do this, run:

    $ make
  
  In this case, the executables, headers, libraries, etc. will be installed
  according to the Pods core policy. If no other "build" directory is found,
  it will default to: lcm-log2smat/build/


###Uninstall
-----

  If compiled and installed to your system, run:

    $ sudo make clean

  Note: this can only be done from the lcm-log2smat source that was used for
  installation


###Examples
-----

  To convert a LCM log file (data.lcm) to a Matlab file (data.mat)

    $ lcm-log2smat data.lcm -o data.mat


  To convert a LCM log file (data.lcm) to Python Pickle file (data.pkd)

    $ lcm-log2smat data.lcm -k  -o data.pkd


###Author
-----

  `lcm-log2smat` is based on [libbot2](https://code.google.com/p/libbot/)
  script `bot-log2mat`.

  Modified and extended by Giancarlo.Troni "at" GMail


###Licenses
-----

(libbot2)

  libbot2 is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published 
  by the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  libbot2 is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.


(lcm-log2smat)

  lcm-log2smat is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published 
  by the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  lcm-log2smat is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.


  You should have received a copy of the GNU Lesser General Public License
  along with lcm-log2smat.  If not, see <http://www.gnu.org/licenses/>.
