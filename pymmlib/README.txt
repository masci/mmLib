mmLib-0.9.5: Python Macromolecular Library (mmLib)
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

Additional requirements for the mmCIF Editor and mmLib Molecular Viewer:
  * PyOpenGL         >= 2.0.0.44 (http://pyopengl.sourceforge.net/)
  * gtk              >= 2.0/2.2  (http://www.gtk.org/)
  * PyGtk            >= 1.99.16  (http://www.pygtk.org/)
  * GtkGLExt         >= 0.7.1    (http://gtkglext.sourceforge.net)
  * PyGtkGLExt       >= 0.0.2    (http://gtkglext.sourceforge.net)

To hijack CCP4's monomer library:
  * CCP4             >= 4.2.1

HARDWARE::
If you want to use the GLViewer.py module, or use the mmLib Molecular 
Viewer, you will need OpenGL hardware acceleration provided by your video
card and driver.  I have tested this on a NVIDIA 5200FX using RedHat 9.0
and the latest NVIDIA-6106 driver.

INSTALLATION::
Read INSTALL.txt

PROBLEMS/BUG REPORTS/CONTACT::
You can contact us through our SourceForge site, or email us directly at:
    * Ethan Merritt <merritt@u.washington.edu>
    * Jay Painter <jpaint@u.washington.edu>

