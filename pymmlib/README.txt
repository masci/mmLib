PyMMLib-0.3pre1: Macromolecular Library (mmLib)
http://pymmlib.sourceforget.net/
-------------------------------------------------------------------------------

Requirements:
  * Python              2.2.1 (the iterators use yield, True/False)
  * Numerical Python    22.0  (seems to work fine on 20.0)
  * Scientific Python   2.4.1 (only used for the Vector class)

Additional requirements for GUI App/Viewer mmView.py:
  * PyOpenGL            2.0.0.44
  * gtk                 2.0/2.2
  * PyGtk               1.99.16
  * GtkGLExt            0.7.1
  * PyGtkGLExt          0.0.2

To hijack CCP4's monomer library:
  * CCP4             >= 4.2.1


The Python Macromolecular Library (PyMMLib) is a collection of Python 
modules the examination and manipulation of macromolecular structures,
and the files which describe them.

Included is a complete PDB v2.2, and mmCIF parser.  The parsers can 
read and write to both formats.  Although several Python PDB parsers
already exist, this is the only one supporting all PDB v2.2 records.

Also included is a new data structure PDB and mmCIF files can be loaded
into for manipulation and examination.  For more information, see
the website or some of the example programs included.

