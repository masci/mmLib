## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Classes which hold information about experimental method and refinement.
"""


class StructureXRayData(dict):
    """Data describing the XRay experiment

    program             (string)
    resolution_high     (float)
    resolution_low      (float)
    data_cutoff         (float) SIGMA(F)
    data_cutoff_high    (float) ABS(F)
    data_cutoof_low     (float) ABS(F)
    completeness        (float) working + test %
    reflections         (int)   
    cross_val_method    (string) cross validation method
    R_free_sel_details  (string) R free selection details
    R                   (float)
    free_R              (float)
    free_R_set_size     (float) percent reflections in Rfree
    free_R_set_count    (int)   number of free reflections
    free_R_error        (float)

    ## need shell data here ##

    num_atoms_protein   (int)
    num_atoms_nucleic_acid (int)
    num_atoms_ligand    (int)
    num_atoms_solvent   (int)

    B_iso_Wilson        (float)
    B_iso_mean          (float)
    B11                 (float)
    B22                 (float)
    B33                 (float)
    B12                 (float)
    B13                 (float)
    B23                 (float)
    """


class StructureData(dict):
    """Contains non-structural data associated with a structure.
    Attributes:

    id                  (string) PDB ID
    keywords            (string) 
    date_original       (string) date of original deposition
    title               (string) 
    author              (string)
    exp_method          (string) experimental method
    

    pdbx_keywords       (string) PDB keywords
    pdbx_description    (string) PDB molecule description
    pdbx_fragment       (string) PDB fragment keywords
    pdbx_ec             (string) PDB ec keywords
    pdbx_mutation       (string) PDB mutation keywords
    """


class StructureSource(dict):
    """
    Data about the source of the structure.

    src_fragment        (string) gene fragment
    src_organism        (string) organism name
    src_organism_common (string) common organism name
    src_genus           (string) genus of source organism
    src_species         (string) species of source organism
    src_strain          (string) strain of source orgaanism
    src_variant         (string)
    src_cell_line       (string)
    src_atcc            (string)
    src_organ           (string)
    src_tissue          (string)
    src_cell            (string)
    src_organelle       (string)
    src_secretion       (string)
    src_cellular_loc    (string) cellular location
    src_plasmid         (string)
    src_gene            (string)
    src_details         (string)

    Source Expression System:
    src_es              (string)
    src_es_common       (string)
    src_es_genus        (string)
    src_es_species      (string)
    src_es_strain       (string)
    src_es_variant      (string)
    src_es_cell_line    (string)
    src_es_atcc_num     (string)
    src_es_organ        (string)
    src_es_tissue       (string)
    src_es_cell         (string)
    src_es_organelle    (string)
    src_es_cellular_loc (string)
    src_es_vector_type  (string)
    src_es_vector       (string)
    src_es_plasmid      (string)
    src_es_gene
    """

    
class PDBData(object):
    """
    id               (string) PDB ID 
    date             (string) original deposition date
    keywords         (string) structure keywords
    pdbx_keywords    (string) PDB keywords
    title            (stirng) title
    """
    def __init__(self):
        self.id            = ""
        self.date          = ""
        self.keywords      = ""
        self.pdbx_keywords = ""
        self.title         = ""


class XRayData(object):
    """
    R_fact           (float)  R factor
    free_R_fact      (float)  free R factor
    res_high         (float)  highest resolution of structure data
    res_low          (float)  lowest resolution of structure
    """
    def __init__(self):
        self.R_fact      = None
        self.free_R_fact = None
        self.res_high    = None
        self.res_low     = None


class NMRData(object):
    pass
