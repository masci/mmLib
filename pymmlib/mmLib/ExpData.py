## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Classes which hold information about experimental method and refinement.
"""


class ExpData(dict):
    """Stroage class for experimental information from the structure.  This
    class inherits from the Python dictionary class, and exists mostly
    for the purpose of documenting the common attributes and providing
    for future expansion.

    key                 definition
    --------------------------------------------------------------------------
    id                  (string) PDB ID
    date                (string) date of original deposition
    title               (string) 
    author              (string)
    exp_method          (string) experimental method
    cif_file            (string) the mmCIFFile object used
    pdb_file            (string) the PDBFile object used
    """
    def __init__(self):
        dict.__init__(self)

        self["id"] = "XXX"
        self["date"] = "unknown"
        self["title"] = "none"
        self["author"] = "unknown"
        self["exp_method"] = "unknown"
