## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Load and save mmLib.Structure objects from/to mmLib supported formats.
The mmCIF and PDB file formats are currently supported.
"""
import os
import types

from mmCIF        import mmCIFFile
from mmCIFBuilder import mmCIFStructureBuilder, mmCIFFileBuilder
from PDB          import PDBFile
from PDBBuilder   import PDBStructureBuilder, PDBFileBuilder


def OpenFile(path, mode):
    """Right now this only supports opening GZip'ed files, in the future
    it might be extended for URLs.
    """
    ## if path is not a string, assume it is a file object and
    ## return it
    
    if isinstance(path, str):
        base, ext = os.path.splitext(path)
        if ext == ".gz":
            import gzip
            return gzip.open(path, mode)
        return open(path, mode)

    return path


def decode_format(path):
    """Returns the 3-letter MIME code for the file format.
    """
    ## check/remove compressed file extention
    if path and isinstance(path, str):
        base, ext = os.path.splitext(path)
        if ext.lower() in ['.z', '.gz', '.bz2']:
            path = base

        base, ext = os.path.splitext(path)
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
    struct = <mmLib.Structure object to build on; defaults to createing new>
    build_properties = <tuple of strings for the StructureBuilder>
    """
    try:
        fil = args["fil"]
    except KeyError:
        raise TypeError,"LoadStructure(fil=) argument required"

    if isinstance(fil, str):
        fileobj = OpenFile(fil, "r")
    else:
        fileobj = fil

    update_cb        = args.get("update_cb")
    format           = args.get("format") or decode_format(fil)
    struct           = args.get("struct") or args.get("structure")
    build_properties = args.get("build_properties", ())

    if format == "PDB":
        return PDBStructureBuilder(
            fil              = fileobj,
            update_cb        = update_cb,
            build_properties = build_properties,
            struct           = struct).struct

    elif format == "CIF":
        return mmCIFStructureBuilder(
            fil              = fileobj,
            update_cb        = update_cb,
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

    if isinstance(fil, str):
        fileobj = OpenFile(fil, "w")
    else:
        fileobj = fil
        
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
        pdb_file.save_file(fileobj)
        return

    elif format == "CIF":
        cif_file = mmCIFFile()
        mmCIFFileBuilder(struct, cif_file)
        cif_file.save_file(fileobj)
        return

    raise FileLoaderError, "Unsupported file format %s" % (str(fil))


### <TESTING>
def test_module():
    import sys
    struct = LoadStructure(fil = sys.argv[1])
    SaveStructure(fil=sys.stdout, struct=struct)

if __name__ == "__main__":
    test_module()
### </TESTING>
