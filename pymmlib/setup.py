#!/usr/bin/env python
## Copyright 2002-2010 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys
import copy

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
        ("/usr/lib64",     "/usr/include"),
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
    import sys
    if sys.platform == "darwin":
        shared_ext = "dylib"
        static_ext = "a"
    elif sys.platform == "win32":
        shared_ext = "dll"
        static_ext = "lib"
    else:
        shared_ext = "so"
        static_ext = "a"

    PATHS = assemble_paths_list()

    shared_lib = "lib%s.%s" % (library, shared_ext)
    static_lib = "lib%s.%s" % (library, static_ext)

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


def library_data(opts):
    """Install mmLib/Data/Monomer library.
    """
    ## start with the mmLib data files
    inst_list = [
        (os.path.join("mmLib", "Data"), [
        os.path.join(os.curdir, "mmLib", "Data", "elements.cif"),
        os.path.join(os.curdir, "mmLib", "Data", "monomers.cif") ])
    ]

    if opts["zip"]:
        inst_list.append((os.path.join("mmLib", "Data"),
            [ os.path.join(os.curdir, "mmLib", "Data", "Monomers.zip") ]))
    else:
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

        if lib_path is None and inc_path is None:
            print "    ERROR: Library Not Found: %s" % (lib)
            return None
        elif lib_path is not None and inc_path is None:
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


def extension_list(opts):
    """Assemble the list of C extensions which need to be compiled.
    """
    print "="*79
    print "Checking Libraries and Headers for C Compiled Modules"
    print "-"*79

    ext_list = []

    if not opts["pdb"]:
        print "  USER_OPTION: PDB module will NOT be built."
    else:
        pdbmodule = pdbmodule_extension()
        if pdbmodule is not None:
            print "  SUCCESS: PDB Module will be built."
            ext_list.append(pdbmodule)
        else:
            print "  FAILURE: PDB Module will NOT be built."

    print "-"*79

    if not opts["opengl"]:
        print "  USER OPTION: OpenGL module will NOT be built."
    else:
        glaccel = glaccel_extension()
        if glaccel is not None:
            print "  SUCCESS: OpenGL Module will be built."
            ext_list.append(glaccel)
        else:
            print "  FAILURE: OpenGL Module will NOT be built."

    print "="*79

    return ext_list


def run_setup(opts):
    """Invoke the Python Distutils setup function.
    """
    s0 = setup(
        cmdclass     = {'install_data': package_install_data},
        name         = "pymmlib",
        version      = "1.2.1",
        author       = "Jay Painter",
        author_email = "jpaint@u.washington.edu",
        url          = "http://pymmlib.sourceforge.net/",
        packages     = ["mmLib"],
        ext_modules  = extension_list(opts),
        data_files   = library_data(opts)
        )


def make_doc(opts):
    """This is a special function to generate the documentation with
    Epydoc.
    """
    ## Requires: http://epydoc.sourceforge.net/
    os.system('epydoc --html --output doc/api_reference --name "Python Macromolecular Library" mmLib/*.py') 

def check_deps_numeric(opts):
    ## check NumPy
    print "Checking for NumPy..."
    try:
        import numpy
    except ImportError:
        print "ERROR: NumPy not found."
    else:
        print "OK:    NumPy found."

def check_deps_opengl(opts):
    if not opts["opengl"]:
        return

    ## check PyOpenGL
    print "Checking for Python OpenGL Bindings..."
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
        print "ERROR: OpenGL.GLUT not found. PyOpenGL may need to be"
        print "       rebuilt after installing GLUT or FreeGLUT."
    else:
        print "OK:    OpenGL.GLUT found."

    ## check PyGTK >= 2.0
    print "Checking for Python GTK+ Bindings..."
    try:
        import pygtk
        pygtk.require("2.0")
    except (ImportError, AssertionError):
        print "ERROR: PyGTK not found. PyGTK 2.0 required."
    else:
        print "OK:    PyGTK found."

    ## check of PyGTKGLExt
    print "Checking for Python GTK+ OpenGL Widget Bindings..."
    try:
        import gtk.gtkgl
    except ImportError:
        print "ERROR: PyGtkGLExt not found."
    else:
        print "OK:    PyGtkGLExt found."

def check_deps(opts):
    """Checks for all required dependancies.
    XXX: This is only checking for the Python modules.
    """

    print "="*79
    print "Running Python Depancy Checks"
    print "-"*79

    ## check Python version
    print "Checking for Python Interpreter..."

    (major, minor, mminor, junk1, junk2) = sys.version_info
    ver_string = "%d.%d.%d" % (major, minor, mminor)

    too_old = False
    if major < 2:
        too_old = True
    elif major == 2:
        if minor < 4:
            too_old = True
        elif minor == 0:
            if mminor < 1:
                too_old = True

    if too_old is True:
        print "ERROR: Python %s too old version >= 2.4.0 required." % (
            ver_string)
    else:
        print "OK:    Python %s Found." % (ver_string)

    ## check Python distutils
    if DISTUTILS_FOUND is False:
        print "ERROR: Python Distutils not found. You may need to install"
        print "       the python-devel package for your Linux distribution."
    else:
        print "OK:    Python Distutils found."

    check_deps_numeric(opts)
    check_deps_opengl(opts)

    print "=" * 79


def buildlib(opts):
    """Download and construct the mmLib monomer library.
    """
    import urllib2

    LIB_FILE = os.path.join("mmLib", "Data", "Monomers.zip")
    LIB_PATH = os.path.join("mmLib", "Data", "Monomers")
    TMP_PATH = "components.cif"
    URL      = "ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif"

    print "[BUILDLIB] downloading %s" % (URL)

    fil = open(TMP_PATH, "w")

    furl = urllib2.urlopen(URL)
    meta = furl.info()
    size = None
    if len(meta.getheaders("Content-Length")):
        size = int(meta.getheaders("Content-Length")[0])/1000

    block_sz = 8192
    downloaded = 0
    while True:
        buffer = furl.read(block_sz)
        if not buffer: break
        downloaded += len(buffer)
        fil.write(buffer)
        print "\r[BUILDLIB] progress: %d KB out of %s    " % (downloaded/1000, 
            '--' if not size else str(size)),  
    fil.close()
    print ''
    print "[BUILDLIB] constructing library from %s" % (TMP_PATH)

    import mmLib.mmCIF

    if opts["zip"]:
        import zipfile
        import cStringIO
        zf = zipfile.ZipFile(LIB_FILE, "w")

    cif_file = mmLib.mmCIF.mmCIFFile()
    cif_file.load_file(TMP_PATH)

    if not os.path.isdir(LIB_PATH):
        os.mkdir(LIB_PATH)

    while len(cif_file) > 0:
        cif_data = cif_file[0]
        cif_file.remove(cif_data)
        cf = mmLib.mmCIF.mmCIFFile()
        cf.append(cif_data)

        if opts["zip"]:
            print "[BUILDLIB] writing %s" % (cif_data.name)
            sf = cStringIO.StringIO()
            cf.save_file(sf)
            zf.writestr(cif_data.name, sf.getvalue())
            sf.close()
        else:
            mkdir_path = os.path.join(LIB_PATH, cif_data.name[0])
            if not os.path.isdir(mkdir_path):
                os.mkdir(mkdir_path)
            save_path = os.path.join(mkdir_path, "%s.cif" % (cif_data.name))
            print "[BUILDLIB] writing %s" % (save_path)
            cf.save_file(save_path)

    if opts["zip"]:
        zf.close()


def check_pymmlib_options():
    import sys

    ## Changed zip=True to zip=False, from Jay's suggestion.
    ## Christoph Champ, 2008-02-14
    opt_defaults = {
        "pdb": True,
        "opengl": True,
        "zip": False,
        }

    opts = {}
    for opt, default in opt_defaults.iteritems():
        opts[opt] = default
        flag = "--with-" + opt
        if flag in sys.argv:
            sys.argv.remove(flag)
            opts[opt] = True
        flag = "--without-" + opt
        if flag in sys.argv:
            sys.argv.remove(flag)
            opts[opt] = False

    return opts


def usage():
    """Print setup.py usage.
    """
    print "SYNOPSIS"
    print "    python setup.py <command>"
    print
    print "DESCRIPTION"
    print "    mmLib build/install program based on the Python"
    print "    distutils package. This setup.py script has several"
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
    print "    buildlib"
    print "        Download the RCSB monomer library from"
    print "        http://pdb.rutgers.edu/public-component-erf.cif and"
    print "        use it to build the mmLib monomer library in"
    print "        mmLib/Data/Monomers"


def main():
    print
    print "PYTHON MACROMOLECULAR LIBRARY -- SETUP PROGRAM"
    print

    opts = check_pymmlib_options()

    if len(sys.argv) == 1:
        usage()

    elif sys.argv[1] == "doc":
        make_doc(opts)

    elif sys.argv[1] == "checkdeps":
        check_deps(opts)

    elif sys.argv[1] == "buildlib":
        buildlib(opts)

    else:
        if DISTUTILS_FOUND is True:
            run_setup(opts)
        else:
            print """\
            ERROR: Python Distuils Not Found. You may have to install
            the python-devel package on some Linux distructions.  See
            INSTALL.txt for details.
            """

if __name__ == "__main__":
    main()
