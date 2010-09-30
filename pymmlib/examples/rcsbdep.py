#!/usr/bin/env python
## Copyright 2002-2010 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import sys
import getopt
import copy
import popen2

from mmLib import mmCIF


##
## RCSB Deposition Required mmCIF Section/Subsections
##
AUDIT_AUTHOR = {
    "table": "audit_author",
    "single_row": False,
    "columns": [{"name": "name"}]
    }

AUDIT_CONTACT_AUTHOR = {
    "table": "audit_contact_author",
    "single_row": False,
    "columns": [{"name": "name"},
                {"name": "email"},
                {"name": "address"},
                {"name": "phone"},
                {"name": "fax"}]
    }

NDB_DATABASE_STATUS = {
    "table": "ndb_database_status",
    "single_row": True,
    "columns": [{"name": "entry_id"},
                {"name":        "dep_release_code_coordinates",
                 "enum_values": ["HOLD FOR PUBLICATION", "RELEASE NOW"] },
                {"name":        "dep_release_code_struct_fact",
                 "enum_values": ["HOLD FOR PUBLICATION", "RELEASE NOW"] },
                {"name":        "dep_release_code_sequence",
                 "enum_values": ["HOLD FOR PUBLICATION", "RELEASE NOW"] } ]
    }

STRUCT = {
    "table": "struct",
    "single_row": True,
    "columns": [{"name": "entry_id"},
                {"name": "title"},
                {"name": "ndb_details"}]
    }

RCSB_DATABASE_RELATED = {
    "table": "rcsb_database_related",
    "single_row": False,
    "columns": [{"name": "db_name"},
                {"name": "db_id"},
                {"name": "details"}]
    }

CITATION_AUTHOR = {
    "table": "citation_author",
    "single_row": False,
    "columns": [{"name": "citation_id", "default": "primary"},
                {"name": "name"}]
    }

CITATION = {
    "table": "citation",
    "single_row": False,
    "columns": [{"name": "journal_abbrev"},
                {"name": "title"},
                {"name": "year"},
                {"name": "journal_volume"},
                {"name": "page_first"},
                {"name": "page_last"}]
    }

ENTITY = {
    "table": "entity",
    "single_row": True,
    "columns": [{"name": "id"},
                {"name": "ndb_description"}]
    }

ENTITY_KEYWORDS = {
    "table": "entity_keywords",
    "single_row": True,
    "columns": [{"name": "entity_id"},
                {"name": "ndb_fragment"},
                {"name": "ndb_mutation"},
                {"name": "ndb_ec"}]
    }

ENTITY_POLY = {
    "table": "entity_poly",
    "single_row": False,
    "columns": [{"name": "ndb_chain_id"},
                {"name": "ndb_seq_one_letter_code"}]
    }

ENTITY_SRC_GEN = {
    "table": "entity_src_gen",
    "single_row": True,
    "columns": [{"name": "entry_id"},
                {"name": "ndb_gene_src_scientific_name"},
                {"name": "ndb_gene_src_gene"},
                {"name": "ndb_host_org_scientific_name"},
                {"name": "ndb_host_org_strain"},
                {"name": "ndb_host_org_vector_type"},
                {"name": "plasmid_name"},
                {"name": "gene_src_details"}],
    "comment": "Source table for genetically manipulated sources."
    }

ENTITY_SRC_NAT = {
    "table": "entity_src_nat",
    "single_row": True,
    "columns": [{"name": "entry_id"},
                {"name": "ndb_organism_scientific"},
                {"name": "strain"},
                {"name": "details"}],
    "comment": "Source table for natural sources."
    }

RCSB_ENTITY_SRC_SYN = {
    "table": "rcsb_entity_src_syn",
    "single_row": True,
    "columns": [{"name": "entry_id"},
                {"name": "details"}],
    "comment": "Source table for synthetic sources."
    }

STRUCT_KEYWORDS = {
    "table": "struct_keywords",
    "single_row": True,
    "columns": [{"name": "entry_id"},
                {"name": "text"}]
    }

STRUCT_BIOL = {
    "table": "struct_biol",
    "single_row": False,
    "columns": [{"name": "id"},
                {"name": "details"}],
    "comment": "Symmetry operators for the construction of the biological unit"
    }

EXPTL_CRYSTAL_GROW = {
    "table": "exptl_crystal_grow",
    "single_row": True,
    "columns": [{"name": "entry_id"},
                {"name": "method"},
                {"name": "pH"},
                {"name": "temp"},
                {"name": "rcsb_details"}]
    }

EXPTL_CRYSTAL = {
    "table": "exptl_crystal",
    "single_row": True,
    "columns": [{"name": "density_percent_sol"},
                {"name": "density_Matthews"}]
    }

CELL = {
    "table": "cell",
    "single_row": True,
    "columns": [{"name": "entry_id"},
                {"name": "length_a"},
                {"name": "length_b"},
                {"name": "length_c"},
                {"name": "angle_alpha"},
                {"name": "angle_beta"},
                {"name": "angle_gamma"}]
    }

SYMMETRY = {
    "table": "symmetry",
    "single_row": True,
    "columns": [{"name": "entry_id"},
                {"name": "space_group_name_H-M"}]
    }

EXPTL = {
    "table": "exptl",
    "single_row": True,
    "columns": [{"name": "entry_id"},
                {"name": "crystals_number", "default": "1"}]
    }

DIFFRN = {
    "table": "diffrn",
    "single_row": True,
    "columns": [{"name": "id"},
                {"name": "ambient_temp"}]
    }

DIFFRN_SOURCE = {
    "table": "diffrn_source",
    "single_row": True,
    "columns": [{"name": "diffrn_id"},
                {"name": "source"},
                {"name": "type"},
                {"name": "rcsb_wavelength_list"}]
    }

DIFFRN_DETECTOR = {
    "table": "diffrn_detector",
    "single_row": True,
    "columns": [{"name": "diffrn_id"},
                {"name": "detector"},
                {"name": "type"},
                {"name": "ndb_collection_date"},
                {"name": "details"}]
    }

DIFFRN_RADIATION = {
    "table": "diffrn_radiation",
    "single_row": True,
    "columns": [{"name": "diffrn_id"},
                {"name": "rcsb_diffrn_protocol",
                 "enum_values": ["SINGLE WAVELENGTH", "LAUE", "MAD"]},
                {"name": "monochromator"}]
    }

REFLNS = {
    "table": "reflns",
    "single_row": True,
    "columns": [{"name": "observed_criterion_sigma_F"},
                {"name": "observed_criterion_sigma_I"},
                {"name": "d_resolution_high"},
                {"name": "d_resolution_low"},
                {"name": "ndb_Rmerge_I_obs"},
                {"name": "ndb_Rsym_value"},
                {"name": "ndb_netI_over_av_sigmaI"},
                {"name": "B_iso_Wilson_estimate"},
                {"name": "ndb_redundancy"}]
    }

REFLNS_SHELL = {
    "table": "reflns_shell",
    "single_row": False,
    "columns": [{"name": "d_res_high"},
                {"name": "d_res_low"},
                {"name": "percent_possible_all"},
                {"name": "Rmerge_I_obs"},
                {"name": "meanI_over_sigI_obs"},
                {"name": "ndb_Rsym_value"},
                {"name": "ndb_redundancy"},
                {"name": "number_unique_all"}]
    }

REFINE = {
    "table": "refine",
    "single_row": True,
    "columns": [{"name": "entry_id"},
                {"name": "ls_d_res_high"},
                {"name": "ls_d_res_low"},
                {"name": "ndb_ls_sigma_F"},
                {"name": "ndb_ls_sigma_I"},
                {"name": "ls_number_reflns_all"},
                {"name": "ls_number_reflns_obs"},
                {"name": "ls_number_reflns_R_free"},
                {"name": "ls_R_factor_all"},
                {"name": "ls_R_factor_obs"},
                {"name": "ls_R_factor_R_work"},
                {"name": "ls_R_factor_R_free"},
                {"name": "ndb_method_to_determine_struct"},
                {"name": "ndb_starting_model"},
                {"name": "ndb_ls_cross_valid_method"},
                {"name": "ndb_R_Free_selection_details"},
                {"name": "ndb_stereochemistry_target_values"},
                {"name": "ndb_isotropic_thermal_model"},
                {"name": "B_iso_mean"},
                {"name": "aniso_B[1][1]"},
                {"name": "aniso_B[1][2]"},
                {"name": "aniso_B[1][3]"},
                {"name": "aniso_B[2][2]"},
                {"name": "aniso_B[2][3]"},
                {"name": "aniso_B[3][3]"},
                {"name": "details"}]
    }

REFINE_LS_SHELL = {
    "table": "refine_ls_shell",
    "single_row": True,
    "columns": [{"name": "d_res_high"},
                {"name": "d_res_low"},
                {"name": "number_reflns_obs"},
                {"name": "number_reflns_R_free"},
                {"name": "R_factor_R_work"},
                {"name": "R_factor_R_free"},
                {"name": "R_factor_R_free_error"},
                {"name": "percent_reflns_obs"}]
    }

REFINE_LS_RESTR = {
    "table": "refine_ls_restr",
    "single_row": True,
    "columns": [{"name": "type"},
                {"name": "dev_ideal"}]
    }

REFINE_ANALYZE = {
    "table": "refine_analyze",
    "single_row": True,
    "columns": [{"name": "Luzzati_coordinate_error_obs"},
                {"name": "Luzzati_sigma_a_obs"},
                {"name": "Luzzati_d_res_low_obs"},
                {"name": "Luzzati_coordinate_error_free"},
                {"name": "Luzzati_sigma_a_free"}]
    }

SOFTWARE = {
    "table": "software",
    "single_row": False,
    "columns": [{"name": "entry_id"}]
    }

COMPUTING = {
    "table": "computing",
    "single_row": True,
    "columns": [{"name": "entry_id"},
                {"name": "data_collection"},
                {"name": "data_reduction"},
                {"name": "structure_solution"},
                {"name": "structure_refinement"}]
    }

DEPOSITION_TABLES = [
    AUDIT_AUTHOR,
    AUDIT_CONTACT_AUTHOR,
    NDB_DATABASE_STATUS,
    STRUCT,
    STRUCT_KEYWORDS,
    STRUCT_BIOL,
    RCSB_DATABASE_RELATED,
    CITATION_AUTHOR,
    CITATION,
    ENTITY,
    ENTITY_KEYWORDS,
    ENTITY_POLY,
    ENTITY_SRC_GEN,
    ENTITY_SRC_NAT,
    RCSB_ENTITY_SRC_SYN,
    EXPTL_CRYSTAL_GROW,
    CELL,
    SYMMETRY,
    EXPTL,
    DIFFRN_SOURCE,
    DIFFRN_DETECTOR,
    DIFFRN,
    DIFFRN_RADIATION,
    REFLNS,
    REFLNS_SHELL,
    REFINE,
    REFINE_LS_SHELL,
    REFINE_LS_RESTR,
    REFINE_ANALYZE,
    SOFTWARE,
    COMPUTING,
    ]

def rcsb_get_deposition_table(table_name):
    for table_desc in DEPOSITION_TABLES:
        if table_desc["table"] == table_name:
            return table_desc
    return None

def table_desc_columns(table_desc):
    columns = []
    for col_desc in table_desc["columns"]:
        columns.append(col_desc["name"])
    return columns

def rcsb_report_missing(cif_data):
    """Reports any missing tables/rows.
    """
    for table_desc in DEPOSITION_TABLES:
        table_name = table_desc["table"]

        cif_table = cif_data.get_table(table_name)

        ## table does not exist
        if cif_table == None:
            print "[ERROR] missing RCSB required table: %s" % (
                table_name)

        ## check table for missing columns
        else:
            for col_desc in table_desc["columns"]:
                col = col_desc["name"]

                for row in cif_table:
                    if not row.has_key(col):
                        print "[ERROR] missing RCSB required column %s.%s" % (
                            table_name, col)


##
## USAGE INFO
##
USAGE = """
rcsbdep.py  - RCSB Deposition Preperation
usage: rcsbdep.py [-f cif path] [-f cif path] ... <mmCIF OUTPUT FILE>

    -e <entry id>  Specify a new entry ID for the structure.
    -f <cif path>  Input mmCIF File
"""

def usage():
    print USAGE
    raise SystemExit

def vet_merge_cif_table_entry_id(cif_table, entry_id):
    """Gives a cif_table a proper vetting before it is merged.
    """
    if cif_table.name == "entry":
        return

    for col in cif_table.columns:
        if col == "entry_id":
            for cif_row in cif_table:
                if cif_row.has_key("entry_id"):
                    cif_row["entry_id"] = entry_id

def add_entry_id(cif_table, entry_id):
    """If there is an entry_id column in the table, ensure all rows
    have an entry_id field.
    """
    if cif_table.has_column("entry_id"):
        for row in cif_table:
            row["entry_id"] = entry_id

def merge_cif_table_multi(dest_table, src_table):
    """Merge the rows from src_table into dest_table and add any
    missing column names to the dest table.
    """
    for row in src_table:
        dest_table.append(copy.deepcopy(row))

        for col in row.keys():
            if not dest_table.has_column(col):
                dest_table.append_column(col)

def null_val(val):
    """Return True if the value is a mmCIF NULL charactor: ? or .
    """
    if val == "?" or val == ".":
        return True
    return False

def merge_cif_table_single(dest_table, src_table):
    """Merge the row from src_table into the row of dest_table and
    add any missing column names to dest_table.
    """
    try:
        src_row = src_table[0]
    except IndexError:
        return

    if len(dest_table) == 0:
        dest_row = dest_table.new_row()
    else:
        dest_row = dest_table[0]

    for key, val in src_row.items():

        ## If the source and destination rows both have values for the same 
        ## field, check for NULL values in the source row so the dest row 
        ## does not have a good value overwriten by a NULL value.
        if dest_row.has_key(key):

            ## if the destination value is null, merge the value without 
            ## complaint
            dval = dest_row[key]

            if null_val(dval):
                dest_row[key] = val

            else:
                ## do not merge a null source value into a non-null
                ## destination value
                if null_val(val):
                    pass

                ## if source and destination values are non-null, then warn 
                ## the user but merge anyway
                else:
                    print "[WARNING] merge overwrite: %s.%s = %s -> %s" % (
                        dest_table.name, key, str(dest_row[key]), str(val))
                    dest_row[key] = val

    for col in src_table.columns:
        if not dest_table.has_column(col):
            dest_table.append_column(col)

def default_table_merge(dest_data, src_table):
    """Merges the src_table into the dest_data mmCIF block. The merging is 
    controlled by optional table definitions.
    """
    ## get a description of this table if it is defined
    table_desc = rcsb_get_deposition_table(src_table.name)

    ## retrieve the destination table if it exists
    dest_table = dest_data.get_table(src_table.name)

    ## mergin a new mmCIF table
    if dest_table == None:

        ## undefined tables are copied in directly
        if table_desc == None:
            dest_data.append(copy.deepcopy(src_table))

        ## defined tables are handled more carefully
        else:
            dest_table = dest_data.new_table(src_table.name)
            dest_table.set_columns(table_desc_columns(table_desc))
            merge_cif_table_multi(dest_table, src_table)

    ## merging rows into an existing table
    else:
        ## undefined tables are merged as multi-row tables
        if table_desc == None:
            merge_cif_table_multi(dest_table, src_table)

        ## defined tables are merged by the defined method
        else:
            if table_desc["single_row"] == True:
                merge_cif_table_single(dest_table, src_table)
            else:
                merge_cif_table_multi(dest_table, src_table)



def cif_merge(merge_paths, output_path, entry_id):
    """Merge the mmCIF files in merge_paths into a mmCIF file which will be 
    written to the output_path. This merging is configured for a 
    crystallographic structure, and the entry_ids of all tables are
    overridden with the given entry_id.
    """
    output_cif_file = mmCIF.mmCIFFile()
    output_cif_data = output_cif_file.new_data(entry_id)

    ## add the entry table to the output mmCIFFile
    entry = output_cif_data.new_table("entry", ["id"])
    row = entry.new_row()
    row["id"] = entry_id

    ## merge in mmCIF files
    for mpath in merge_paths:
        merge_cif_file = mmCIF.mmCIFFile()
        merge_cif_file.load_file(mpath)

        for src_data in merge_cif_file:
            for src_table in src_data:

                ## don't merge the entry table
                if src_table.name == "entry":
                    continue

                vet_merge_cif_table_entry_id(src_table, entry_id)
                default_table_merge(output_cif_data, src_table)

    ## post-process merged mmCIF file

    ## make sure the entry_id fields are filled in
    for cif_table in output_cif_data:
        add_entry_id(cif_table, entry_id)

    ## report missing items
    rcsb_report_missing(output_cif_data)

    output_cif_file.save_file(output_path)

def main():
    ## parse options
    (opts, args) = getopt.getopt(sys.argv[1:], "h?f:e:")

    entry_id    = "XXXX"
    merge_paths = []
    output_path = None

    for (opt, item) in opts:
        if opt == "-h" or opt == "-?":
            usage()
        elif opt == "-f":
            merge_paths.append(item)
        elif opt == "-e":
            entry_id = item

    if len(args) != 1:
        print "[ERROR] Output mmCIF path required."
        usage()

    if len(merge_paths) == 0:
        print "[ERROR] Input mmCIF paths required."
        usage()

    output_path = args[0]

    cif_merge(merge_paths, output_path, entry_id)


if __name__ == "__main__":
    main()
