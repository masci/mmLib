## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Load and save mmLib.Structure objects from/to mmLib supported formats.
The mmCIF and PDB file formats are currently supported.
"""
import os

from mmTypes      import *
from mmCIF        import mmCIFFile
from mmCIFBuilder import mmCIFStructureBuilder, mmCIFFileBuilder
from PDB          import PDBFile
from PDBBuilder   import PDBStructureBuilder, PDBFileBuilder


def decode_format(path):
    """Returns the 3-letter MIME code for the file format.
    """
    ## check/remove compressed file extention
    if path and type(path) == StringType:
        (base, ext) = os.path.splitext(path)
        if ext.lower() in ['.z', '.gz', '.bz2']:
            path = base

        (base, ext) = os.path.splitext(path)
        ext = ext.lower()

        if ext == ".cif":
            return "CIF"

    return "PDB"


def LoadStructure(**args):
    """Loads a mmCIF file(.cif) or PDB file(.pdb) into a 
    Structure class and returns it.  The function takes 5 named arguments,
    one is required:

    fil = <file object or path; required>
    format = <PDB|CIF; defaults to PDB>
    library = <mmLib.Library class; defaults to mmLib.Library>
    struct = <mmLib.Structure object to build on; defaults to createing new>
    build_properties = <tuple of strings for the StructureBuilder>
    """
    try:
        fil = args["fil"]
    except KeyError:
        raise TypeError,"LoadStructure(fil=) argument required"

    update_cb        = args.get("update_cb")
    format           = args.get("format") or decode_format(fil)
    library          = args.get("library")
    struct           = args.get("struct") or args.get("structure")
    build_properties = args.get("build_properties", ())

    if format == "PDB":
        return PDBStructureBuilder(
            fil              = fil,
            update_cb        = update_cb,
            library          = library,
            build_properties = build_properties,
            struct           = struct).struct

    elif format == "CIF":
        return mmCIFStructureBuilder(
            fil              = fil,
            update_cb        = update_cb,
            library          = library,
            build_properties = build_properties,
            struct           = struct).struct

    raise FileLoaderError, "Unsupported file format %s" % (str(fil))


def SaveStructure(**args):
    """Saves a Structure object into a supported file type.
    fil = <file object or path; required>
    struct = <mmLib.Structure object to save; required>
    format = <PDB|CIF; defaults to PDB>
    """
    try:
        fil = args["fil"]
    except KeyError:
        raise TypeError,"LoadStructure(fil=) argument required"

    try:
        struct = args["struct"]
    except KeyError:
        try:
            struct = args["structure"]
        except KeyError:
            raise TypeError,"LoadStructure(struct=) argument required"

    format = args.get("format") or decode_format(fil)

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
    struct = LoadStructure(fil = sys.argv[1])
    SaveStructure(fil=sys.stdout, struct=struct)
### </TESTING>
