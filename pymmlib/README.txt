mmLib-0.5: Macromolecular Library (mmLib)
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
  * Scientific Python>= 2.4.1 (only used for the Vector class)

Additional requirements for GUI App/Viewer mmView.py, the mmCIF editor
only uses PyGTK.  It should "just work" on Redhat 9.0 systems.
  * PyOpenGL         >= 2.0.0.44
  * gtk              >= 2.0/2.2
  * PyGtk            >= 1.99.16
  * GtkGLExt         >= 0.7.1
  * PyGtkGLExt       >= 0.0.2

To hijack CCP4's monomer library:
  * CCP4             >= 4.2.1

INSTALL::
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

# tar xzf pymmlib-0.5.tar.gz
# export PYTHONPATH=/home/jpaint/pymmlib-0.5
# cd pymmlib-0.5/tests
# python mmlib_test.py /home/jpaint/myfile.pdb

PROBLEMS/BUG REPORTS/CONTACT::
You can contact us through our SourceForge site, or email us directly at:
    * Ethan Merritt <merritt@u.washington.edu>
    * Jay Painter <jpaint@u.washington.edu>

