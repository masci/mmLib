## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Load and save mmLib.Structure objects from/to mmLib supported formats.
The mmCIF and PDB file formats are currently supported.
"""
import os
from mmTypes import *
from  StructureBuilder import PDBStructureBuilder, mmCIFStructureBuilder
from mmCIF import mmCIFFile, mmCIFFileBuilder
from PDB import PDBFile, PDBFileBuilder


def decode_format(path, format):
    """Returns the 3-letter MIME code for the file format."""
    if format: return format

    ## check/remove compressed file extention
    if path and type(path) == StringType:
        (base, ext) = os.path.splitext(path)
        if ext.lower() in ['.z', '.gz', '.bz2']:
            path = base

        (base, ext) = os.path.splitext(path)
        ext = ext.lower()

        if ext == ".cif": return "CIF"

    return "PDB"


def LoadStructure(fil,
                  library = None,
                  format = "",
                  build_properties = (),
                  struct = None):

    """Loads a mmCIF file(.cif) or PDB file(.pdb) into a 
    Structure class and returns it.  The first argument is either a
    path string or a file object opened for reading.
    """
    format = decode_format(fil, format)

    if format == "PDB":
        return PDBStructureBuilder(
            fil              = fil,
            library          = library,
            build_properties = build_properties,
            struct           = struct).struct

    elif format == "CIF":
        return mmCIFStructureBuilder(
            fil              = fil,
            library          = library,
            build_properties = build_properties,
            struct           = struct).struct

    raise FileLoaderError, "Unsupported file format %s" % (str(fil))


def SaveStructure(fil, struct, format = ""):
    """Saves a Structure object into a supported file type.
    """
    format = decode_format(fil, format)

    if format == "PDB":
        pdb_file = PDBFile()
        PDBFileBuilder(struct, pdb_file)
        pdb_file.save_file(fil)
        return

    elif format == "CIF":
        cif_file = mmCIFFile()
        mmCIFFileBuilder(struct, cif_file)
        cif_file.save_file(fil)
        return

    raise FileLoaderError, "Unsupported file format %s" % (str(fil))


### <TESTING>
if __name__ == "__main__":
    import sys
    struct = LoadStructure(sys.argv[1], build_properties=("polymers","bonds"))
    SaveStructure(sys.stdout, struct, "PDB")
### </TESTING>
