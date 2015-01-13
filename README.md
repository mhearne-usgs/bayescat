Introduction
------------

bayescat is a program which uses the Lawrence Livermore National
Laboratory (LLNL) BayesLoc software
(https://missions.llnl.gov/nonproliferation/nuclear-explosion-monitoring/bayesloc)
and the NEIC Comprehensive Catalog (ComCat) to allow the user to build
a local re-located earthquake catalog using NEIC events with Phase
Data as input.



Installation and Dependencies
-----------------------------

This package depends on:
 * numpy, the fundamental package for scientific computing with Python. <a href="http://www.numpy.org/">http://www.numpy.org/</a>  
 * matplotlib, a Python 2D plotting library which produces publication quality figures. <a href="<a href="http://matplotlib.org/index.html">http://matplotlib.org/index.html</a>
 * obspy, a Seismological data processing package.  <a href="<a href="https://github.com/obspy/obspy/wiki">https://github.com/obspy/obspy/wiki</a>
 * libcomcat, an NEIC package for retrieving data from ComCat.
 * neicio, an NEIC package with a module for running external commands.
 * neicutil, an NEIC package with a module for manipulating text representations of numbers.


The best way to install numpy and matplotlib is to use one of the Python distributions described here:

<a href="http://www.scipy.org/install.html">http://www.scipy.org/install.html</a>

The Anaconda distribution has been successfully tested with bayescat.

Most of those distributions should include <em>pip</em>, a command line tool for installing and 
managing Python packages.  You will use pip to install the other dependencies and libcomcat itself.  
 
You may need to open a new terminal window to ensure that the newly installed versions of python and pip
are in your path.

To install obspy:

pip install obspy

To install neicio:

pip install git+git://github.com/usgs/neicio.git

To install neicutil:

pip install git+git://github.com/usgs/neicutil.git

To install libcomcat:

pip install git+git://github.com/usgs/libcomcat.git

To install this package:

pip install git+git://github.com/mhearne-usgs/bayescat.git

Uninstalling and Updating
-------------------------

To uninstall:

pip uninstall bayescat

To update:

pip install git+git://github.com/mhearne-usgs/bayescat.git

Configuration
-------------
This software requires the user to have previously downloaded and compiled BayesLoc following
the instructions provided by LLNL.  Then the user must create the following directory structure:

~/bayesloc/bin/bayesloc (the executable)
~/bayesloc/events/
~/bayesloc/ttimes/ (put here the ak135.* files from ftp://ftpext.cr.usgs.gov/pub/cr/co/golden/benz/ak135.tar)

Usage
-----
<pre>
usage: bayescat.py [-h] [-i ID] [-r RADIUS] [-s] [-c] [-f]
                   [-d [DELETE [DELETE ...]]]

Use BayesLoc to help automate the creation of a relocated earthquake catalog.

optional arguments:
  -h, --help            show this help message and exit
  -i ID, --id ID        The NEIC event ID to relocate (default: None)
  -r RADIUS, --radius RADIUS
                        The radial search distance. (default: 15)
  -s, --stats           Display database statistics (number of events,
                        stations, etc.) (default: False)
  -c, --count           Display number of events in search radius. (default:
                        False)
  -f, --force           Force relocation of an already relocated event.
                        (default: False)
  -d [DELETE [DELETE ...]], --delete [DELETE [DELETE ...]]
                        Delete event(s) from database. (default: None)
</pre>

