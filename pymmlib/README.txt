mmLib-0.7: Macromolecular Library (mmLib)
http://pymmlib.sourceforget.net/
-------------------------------------------------------------------------------

DESCRIPTION::
The Python Macromolecular Library (mmLib) is a collection of Python 
modules the examination and manipulation of macromolecular structures,
and the files which describe them.

Included is a complete PDB v2.2, and mmCIF parser.  The parsers can 
read and write to both formats.  Although several Python PDB parsers
already exist, this is the only one supporting all PDB v2.2 records.

Also included is a new data structure PDB and mmCIF files can be loaded
into for manipulation and examination.  For more information, see
the website or some of the example programs included.

REQUIREMENTS::
  * Python           >= 2.2.1 (the iterators use yield, True/False)
    NOTE: I use Python 2.3, and it is about 20% faster than 2.2.1!
  * Numerical Python >= 22.0  (seems to work fine on 20.0)

Additional requirements for GUI App/Viewer viewer.py, the mmCIF editor
only uses PyGTK.  It should "just work" on Redhat 9.0 systems.
  * PyOpenGL         >= 2.0.0.44
  * gtk              >= 2.0/2.2
  * PyGtk            >= 1.99.16
  * GtkGLExt         >= 0.7.1
  * PyGtkGLExt       >= 0.0.2

To hijack CCP4's monomer library:
  * CCP4             >= 4.2.1

INSTALL::
  PLEASE NOTE:
    Many Linux distributions split the Python distribution into
    separate packages.  Most Linux distributions use the Redhat
    Package Manager (.rpm) format, and have a python-*.rpm 
    and python-dev-*.rpm/python-devel-*.rpm packages for the 
    complete distribution.  However, I've noticed some 
    Linux distributions are missing importent Python libraries
    even when both of these packages are installed, causing the
    mmLib installer to fail.  

    For this reason, you may have to compile and install Python from 
    source.  This may sound difficult, but in my experience it 
    is usually easier than finding the python-devel-*.rpm for
    your distribution.  I reccomend compiling and installing Python-2.3.3
    from http://www.python.org.  

    Here are some simple instructions if you want to do this:

    # (get Python-2.3.3.tgz from http://www.python.org, save to your home dir)
    # tar xzf Python-2.3.3.tgz
    # cd Python-2.3.3
    # ./configure --prefix=/usr/local/python-2.3.3
    # make
    # su (root passwd)
    # make install
    # exit
    # export PATH=/usr/local/python-2.3.3/bin:$PATH
    
    # (get Numeric Python from http://numpy.sourceforge.net/)
    # tar xzf Numeric-23.1.tar.gz
    # cd Numeric-23.1
    # python setup.py build
    # su (root passwd)
    # python setup.py install
    
    # tar xzf pymmlib-0.7.tar.gz
    # cd pymmlib-0.7
    # su (root passwd)
    # python setup.py install
    
    Now you have a working installation of Python-2.3.3, Numeric, and mmLib.

Python includes its own system for installing new modules called distutils.
To use the distutils, one needs to write a setup.py file which acts much
like a Makefile for Python.  The installer does not check for any of the
other required Python modules.  It will install, but the programs will not
run.  There are links to all the required software on the mmLib website.

If you want to install mmLib into your Python distribution's standard
library, run the following command as root:

# python setup.py install

If you want to use mmLib without installing it as root, you can set the
PYTHONPATH environment variable.  For example, to unpack mmLib and run the
example programs and applications:

# tar xzf pymmlib-0.7.tar.gz
# export PYTHONPATH=/home/jpaint/pymmlib-0.7
# cd pymmlib-0.7/tests
# python mmlib_test.py /home/jpaint/myfile.[pdb|cif]

PROBLEMS/BUG REPORTS/CONTACT::
You can contact us through our SourceForge site, or email us directly at:
    * Ethan Merritt <merritt@u.washington.edu>
    * Jay Painter <jpaint@u.washington.edu>

