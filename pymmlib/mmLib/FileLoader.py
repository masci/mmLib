## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

"""This module provides loading and saving of Structure classes from and to
various supported file formats.  The mmCIF and PDB file formats are currently
supported."""

import os
import string
from   mmTypes          import *
from   StructureBuilder import PDBStructureBuilder, mmCIFStructureBuilder
from   PDB              import StructurePDBFileBuilder


def decode_format(path, format):
    """Returns the 3-letter MIME code for the file format."""
    ## check for compressed file
    ext = format
    if type(path) == StringType and format == "":
        (base, ext) = os.path.splitext(path)
        ext = string.lower(ext)
    elif type(path) != StringType and format == "":
        ## default format
        format = ".pdb"

    ## remove the compression extention if it is there
    if ext in ['.z', '.gz', '.bz2']:
        (base, ext) = os.path.splitext(base)
        ext = string.lower(ext)

    ## select adaptor based on file extention
    if ext == '.cif':
        return "CIF"
    return "PDB"


def LoadStructure(fil,
                  library          = None,
                  format           = "",
                  build_properties = ()):
    """Loads a mmCIF file(.cif) or PDB file(.pdb) into a mmPython
    Structure class and returns it.  The first argument is either a
    path string or a file object opened for reading."""
    format = decode_format(fil, format)
    if   format == "PDB":
        return PDBStructureBuilder(
            fil              = fil,
            library          = library,
            build_properties = build_properties).structure
    elif format == "CIF":
        return mmCIFStructureBuilder(
            fil              = fil,
            library          = library,
            build_properties = build_properties).structure
    raise FileLoaderError, "Unsupported file format %s" % (str(fil))


def SaveStructure(fil,
                  structure,
                  format     = ""):
    """Saves a mmPython Structure class into a supported file type."""
    format = decode_format(fil, format)
    if   format == "PDB":
        return StructurePDBFileBuilder(structure).pdb_file.saveFile(fil)
    elif format == "CIF":
        pass
    raise FileLoaderError, "Unsupported file format %s" % (str(fil))


### <TESTING>
if __name__ == "__main__":
    import sys
    struct = LoadStructure(sys.argv[1], build_properties=("polymers","bonds"))
    SaveStructure(sys.stdout, struct)
### </TESTING>
