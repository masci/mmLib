#!/usr/bin/env python

## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys


## try to import the distutils package
try:
    from distutils.core import setup, Extension
    from distutils.command.install_data import install_data
except ImportError:
    DISTUTILS_FOUND = False
else:
    DISTUTILS_FOUND = True
    

class package_install_data(install_data):
    def run(self):
        ## need to change self.install_dir to the actual library dir
        install_cmd = self.get_finalized_command('install')
        self.install_dir = getattr(install_cmd, 'install_lib')
        return install_data.run(self)


def assemble_paths_list():
    """Returns a list of paths to search for libraries and headers.
    """
    PATHS = [
        ("/usr/lib",       "/usr/include"),
        ("/usr/X11/lib",   "/usr/X11/include"),
        ("/usr/X11R6/lib", "/usr/X11R6/include"),
        ("/usr/local/lib", "/usr/local/include"),
        ("/opt/lib",       "/opt/include") ]

    ## add path from environment var LD_LIBRARY_PATH
    try:
        ld_lib_path = os.environ["LD_LIBRARY_PATH"]
    except KeyError:
        pass
    else:
        for path in ld_lib_path.split(":"):
            path = path.strip()
            (dir, first) = os.path.split(path)
            include_path = os.path.join(dir, "include")

            if dir not in PATHS:
                PATHS.append((path, include_path))

    ## add paths from /etc/ld.so.conf
    try:
        ld_so_conf = open("/etc/ld.so.conf", "r").readlines()
    except IOError:
        return PATHS

    for path in ld_so_conf:
        path = path.strip()
        (dir, first) = os.path.split(path)
        include_path = os.path.join(dir, "include")

        if dir not in PATHS:
            PATHS.append((path, include_path))

    return PATHS


def find_lib_paths(library, include):
    """Find the location of the library.
    """
    import glob

    PATHS = assemble_paths_list()

    shared_lib = "lib%s.so" % (library)
    static_lib = "lib%s.a" % (library)

    found_lib_path = None
    found_inc_path = None
    found_library  = None
        
    for lib_path, inc_path in PATHS:
        lib_glob    = os.path.join(lib_path, shared_lib)
        shared_libs = glob.glob(lib_glob)

        if len(shared_libs)==0:
            print lib_glob
            continue

        print lib_glob + "*"

        found_lib_path = lib_path
        found_library  = shared_libs[0]
        
        inc_check = os.path.join(inc_path, include)

	print "Checking Library/Include: %s %s" % (found_library, inc_check)

        if not os.path.isfile(inc_check):
            continue

        found_inc_path = inc_path
        break
        
    return found_lib_path, found_inc_path, found_library


def library_data():
    """Install mmLib/Data/Monomer library.
    """
    ## start with the mmLib data files
    inst_list = [
        (os.path.join("mmLib", "Data"), [
        os.path.join(os.curdir, "mmLib", "Data", "elements.cif"),
        os.path.join(os.curdir, "mmLib", "Data", "monomers.cif") ])
    ]
  
    ## add all the monomer mmCIF files 
    mon_dir = os.path.join(os.curdir, "mmLib", "Data", "Monomers")
    
    for dir1 in os.listdir(mon_dir):
        dir2 = os.path.join(mon_dir, dir1)

        inst_dir = os.path.join("mmLib", "Data", "Monomers", dir1)

        file_list = []
        inst_list.append((inst_dir, file_list)) 
        
        for fil in os.listdir(dir2):
            file_list.append(os.path.join(dir2, fil))
            
    return inst_list


def pdbmodule_extension():
    """Add the PDB Accelerator module.
    """
    print "PDB C Module (pdbmodule)"
    ext = Extension(
        "mmLib.pdbmodule", 
        ["src/pdbmodule.c"],
        include_dirs = [],
        library_dirs = [],
        libraries    = [])
    
    return ext


def glaccel_extension():
    """Add the OpenGL Accelorator module.
    """
    print "OpenGL C Module (glaccel)"
    
    glaccel_libs = [
        {"library":    "m",
         "header":     "math.h",
         "desc":       "Math Library"},

        {"library":    "X11",
         "header":     "X11/X.h",
         "desc":       "X11 Client Library"},

        {"library":    "Xext",
         "header":     "X11/extensions/Xext.h",
         "desc":       "X11 Extensions Library"},

    #    {"library":    "Xi",
    #     "header":     "X11/extensions/XIE.h",
    #     "desc":       "X11 Extended Input Library"},

        {"library":    "Xmu",
         "header":     "X11/Xmu/Xmu.h",
         "desc":       "X11 Misc Library"},
        
        {"library":    "GL",
         "header":     "GL/gl.h",
         "desc":       "OpenGL Library"},

        {"library":    "GLU",
         "header":     "GL/glu.h",
         "desc":       "OpenGLU Library"},
        
        {"library":    "glut",
         "header":     "GL/glut.h",
         "desc":       "OpenGL Utility Library (or FreeGLUT)"}
        ]

    include_dirs = []
    library_dirs = []
    libraries    = []

    for lib_dict in glaccel_libs:

        lib  = lib_dict["library"]
        inc  = lib_dict["header"]
        desc = lib_dict["desc"]

        lib_path, inc_path, library_path = find_lib_paths(lib, inc)

        print "  Searching For: %s" % (desc)

        if lib_path==None and inc_path==None:
            print "    ERROR: Library Not Found: %s" % (lib)
            return None
        elif lib_path!=None and inc_path==None:
            print "    ERROR: Header Files Not Found for Library: %s" % (lib)
            return None

        print "    Found Library: %s" % (library_path)

        if inc_path not in include_dirs:
            include_dirs.append(inc_path)

        if lib_path not in library_dirs:
            library_dirs.append(lib_path)

        libraries.append(lib)

    ext = Extension(
        "mmLib.glaccel", 
        ["src/glaccel.c"],
        include_dirs = include_dirs,
        library_dirs = library_dirs,
        libraries    = libraries)

    return ext


def extension_list():
    """Assemble the list of C extensions which need to be compiled.
    """
    print "="*79
    print "Checking Libraries and Headers for C Compiled Modules"
    print "-"*79

    ext_list = []

    pdbmodule = pdbmodule_extension()
    if pdbmodule!=None:
        print "  SUCCESS:  Module will be built."
        ext_list.append(pdbmodule)
    else:
        print "  FAILURE: Module will NOT be built."

    print "-"*79

    glaccel = glaccel_extension()
    if glaccel!=None:
        print "  SUCCESS: Module will be built."
        ext_list.append(glaccel)
    else:
        print "  FAILURE: Module will NOT be built."

    print "="*79
    
    return ext_list


def run_setup():
    """Invoke the Python Distutils setup function.
    """
    s0 = setup(
        cmdclass = {'install_data': package_install_data},
        
        name         = "pymmlib",
        version      = "0.9.5",
        author       = "Jay Painter",
        author_email = "jpaint@u.washington.edu",
        url          = "http://pymmlib.sourceforge.net/",
        packages     = ["mmLib", "mmLib/Extensions"],
        ext_modules  = extension_list(),
        data_files   = library_data()
        )


def make_doc():
    """This is a special function to generate the documentation with
    Epidoc.  It once used happydoc.
    """
    #os.system("happydoc -d doc/happy -t 'PyMMLib Documentation' --no-comments --no-private-names mmLib")

    os.system(
        'epydoc --html --output doc/api_reference '\
        '--name "mmLib Documentation" mmLib mmLib/Extensions') 


def check_deps():
    """
    Checks for all required dependancies.
    XXX: This is only checking for the Python modules.
    """

    print "="*79
    print "Running Depancy Checks"
    print "-"*79

    ## check Python version
    (major, minor, mminor, junk1, junk2) = sys.version_info
    ver_string = "%d.%d.%d" % (major, minor, mminor)

    too_old = False
    if major<2:
        too_old = True
    elif major==2:
        if minor<2:
            too_old = True
        elif minor==2:
            if mminor<1:
                too_old = True
            
    if too_old==True:
        print "ERROR: Python %s too old version >= 2.2.1 required." % (
            ver_string)
    else:
        print "OK:    Python %s Found." % (ver_string)

    ## check Python distutils
    if DISTUTILS_FOUND==False:
        print "ERROR: Python Distutils not found.  You may need to install"
        print "       the python-devel package for your Linux distribution."
    else:
        print "OK:    Python Distutils found."

    ## check Python Numeric
    try:
        import Numeric
    except ImportError:
        print "ERROR: Numeric Python not found."
    else:
        print "OK:    Numeric Python found."

    ## check PyOpenGL
    try:
        import OpenGL.GL
    except ImportError:
        print "ERROR: OpenGL.GL not found."
    else:
        print "OK:    OpenGL.GL found."

    try:
        import OpenGL.GLU
    except ImportError:
        print "ERROR: OpenGL.GLU not found."
    else:
        print "OK:    OpenGL.GLU found."

    try:
        import OpenGL.GLUT
    except ImportError:
        print "ERROR: OpenGL.GLUT not found.  PyOpenGL may need to be"
        print "       rebuilt after installing GLUT or FreeGLUT."
    else:
        print "OK:    OpenGL.GLUT found."

    ## check PyGTK >= 2.0
    try:
        import pygtk
        pygtk.require("2.0")
    except (ImportError, AssertionError):
        print "ERROR: PyGTK not found.  PyGTK 2.0 required."
    else:
        print "OK:    PyGTK found."

    ## check of PyGTKGLExt
    try:
        import gtk.gtkgl
    except ImportError:
        print "ERROR: PyGtkGLExt not found."
    else:
        print "OK:    PyGtkGLExt found."

    print "="*79


def usage():
    """Print setup.py usage.
    """
    print "SYNOPSIS"
    print "    python setup.py <command>"
    print
    print "DESCRIPTION"
    print "    mmLib build/install program based on the Python"
    print "    distutils package.  This setup.py script has several"
    print "    features not found in most Python seutp.py scrips."
    print
    print "COMMANDS"
    print "    build"
    print "        Comples/builds mmLib C modules and preps Python"
    print "        files for installation."
    print "    install"
    print "        Installs all mmLib files into the Python library"
    print "        directory."
    print "    checkdeps"
    print "        Checks for all dependant Python modules used by"
    print "        the mmLib library and associated applications."
    print "        This produces a nice report to help you figure out"
    print "        what dependancies need to be installed."
    print "        See INSTALL.txt for more details."
    print "    doc"
    print "        Build the mmLib developers documentation using the"
    print "        Epidoc program."
    print
    

if __name__ == "__main__":

    print
    print "PYTHON MACROMOLECULAR LIBRARY -- SETUP PROGRAM"
    print

    if len(sys.argv)==1:
        usage()

    elif sys.argv[1] == "doc":
        make_doc()

    elif sys.argv[1] == "checkdeps":
        check_deps()

    else:
        if DISTUTILS_FOUND==True:
            run_setup()
        else:
            print """\
            ERROR: Python Distuils Not Found.  You may have to install
            the python-devel package on some Linux distructions.  See
            INSTALL.txt for details.
            """
