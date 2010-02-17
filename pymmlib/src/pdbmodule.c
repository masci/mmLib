/* Copyright 2002-2010 by PyMMLib Development Group (see AUTHORS file)
 * This code is part of the PyMMLib distribution and governed by
 * its license.  Please see the LICENSE file that should have been
 * included as part of this package.
 *
 * pdbmodule.c - PDB parser/accelorator for mmLib
 *
 */
#include "Python.h"
#include <stdio.h>

#define MAX_LINE 82

static PyObject *PDBModuleErr = NULL;

typedef enum {
  PDB_FIELD_TYPE_STRING,
  PDB_FIELD_TYPE_INTEGER,
  PDB_FIELD_TYPE_FLOAT_2,
  PDB_FIELD_TYPE_FLOAT_3,
  PDB_FIELD_TYPE_FLOAT_4,
  PDB_FIELD_TYPE_FLOAT_5,
  PDB_FIELD_TYPE_FLOAT_6
} PDBFieldType;

typedef enum {
  PDB_FIELD_RJUST,
  PDB_FIELD_LJUST
} PDBFieldAlign;

struct PDBFieldDef {
  char          *name;
  int            istart;
  int            iend;
  PDBFieldType   type;
  PDBFieldAlign  align;
};

static struct PDBFieldDef keywds_fields[] = {
  {"continuation", 9, 10, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"keywds",      11, 70, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
};

static struct PDBFieldDef seqres_fields[] = {
  {"serNum",       9, 10, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"chainID",     12, 12, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"numRes",      14, 17, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"resName1",    20, 22, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName2",    24, 26, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName3",    28, 30, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName4",    32, 34, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName5",    36, 38, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName6",    40, 42, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName7",    44, 46, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName8",    48, 50, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName9",    52, 54, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName10",   56, 58, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName11",   60, 62, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName12",   64, 66, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName13",   68, 70, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
};

static struct PDBFieldDef link_fields[] = {
  {"name1",       13, 16, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"altLoc1",     17, 17, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName1",    18, 20, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainID1",    22, 22, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resSeq1",     23, 26, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCode1",      27, 27, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"name2",       43, 46, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"altLoc2",     47, 47, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName2",    48, 50, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainID2",    52, 52, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resSeq2",     53, 56, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCode2",      57, 57, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"sym1",        60, 65, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"sym2",        67, 72, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
};

static struct PDBFieldDef origx1_fields[] = {
  {"o[n][1]",     11, 20, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"o[n][2]",     21, 30, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"o[n][3]",     31, 40, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"t[n]",        46, 55, PDB_FIELD_TYPE_FLOAT_5, PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_FLOAT_5, PDB_FIELD_LJUST},
};

static struct PDBFieldDef endmdl_fields[] = {
  {NULL,           0,  0, PDB_FIELD_TYPE_FLOAT_5, PDB_FIELD_LJUST},
};

static struct PDBFieldDef hetsyn_fields[] = {
  {"continuation", 9, 10, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {"hetID",       12, 14, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"hetSynonyms", 16, 70, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
};

static struct PDBFieldDef sigatm_fields[] = {
  {"serial",       7, 11, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"name",        13, 16, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {"altLoc",      17, 17, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName",     18, 20, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainID",     22, 22, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resSeq",      23, 26, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCode",       27, 27, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"sigX",        31, 38, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"sigY",        39, 46, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"sigZ",        47, 54, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"sigOccupancy",55, 60, PDB_FIELD_TYPE_FLOAT_2, PDB_FIELD_LJUST},
  {"sigTempFactor",61, 66, PDB_FIELD_TYPE_FLOAT_2, PDB_FIELD_LJUST},
  {"segID",       73, 76, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"element",     77, 78, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"charge",      79, 80, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
};

static struct PDBFieldDef mtrix2_fields[] = {
  {"serial",       8, 10, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"s[n][1]",     11, 20, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"s[n][2]",     21, 30, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"s[n][3]",     31, 40, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"v[n]",        46, 55, PDB_FIELD_TYPE_FLOAT_5, PDB_FIELD_LJUST},
  {"iGiven",      60, 60, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
};

static struct PDBFieldDef mtrix3_fields[] = {
  {"serial",       8, 10, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"s[n][1]",     11, 20, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"s[n][2]",     21, 30, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"s[n][3]",     31, 40, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"v[n]",        46, 55, PDB_FIELD_TYPE_FLOAT_5, PDB_FIELD_LJUST},
  {"iGiven",      60, 60, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
};

static struct PDBFieldDef helix_fields[] = {
  {"serNum",       8, 10, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"helixID",     12, 14, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"initResName", 16, 18, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"initChainID", 20, 20, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"initSeqNum",  22, 25, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"initICode",   26, 26, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"endResName",  28, 30, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"endChainID",  32, 32, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"endSeqNum",   34, 37, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"endICode",    38, 38, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"helixClass",  39, 40, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"comment",     41, 70, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {"length",      72, 76, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
};

static struct PDBFieldDef mtrix1_fields[] = {
  {"serial",       8, 10, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"s[n][1]",     11, 20, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"s[n][2]",     21, 30, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"s[n][3]",     31, 40, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"v[n]",        46, 55, PDB_FIELD_TYPE_FLOAT_5, PDB_FIELD_LJUST},
  {"iGiven",      60, 60, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
};

static struct PDBFieldDef turn_fields[] = {
  {"seq",          8, 10, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"turnID",      12, 14, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"initResName", 16, 18, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"initChainID", 20, 20, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"initSeqNum",  21, 24, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"initICode",   25, 25, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"endResName",  27, 29, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"endChainID",  31, 31, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"endSeqNum",   32, 35, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"endICode",    36, 36, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"comment",     41, 70, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
};

static struct PDBFieldDef anisou_fields[] = {
  {"serial",       7, 11, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"name",        13, 16, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {"altLoc",      17, 17, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName",     18, 20, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainID",     22, 22, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resSeq",      23, 26, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCode",       27, 27, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"u[0][0]",     29, 35, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"u[1][1]",     36, 42, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"u[2][2]",     43, 49, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"u[0][1]",     50, 56, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"u[0][2]",     57, 63, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"u[1][2]",     64, 70, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"segID",       73, 76, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"element",     77, 78, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"charge",      79, 80, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
};

static struct PDBFieldDef title_fields[] = {
  {"continuation", 9, 10, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"title",       11, 70, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
};

static struct PDBFieldDef sltbrg_fields[] = {
  {"name1",       13, 16, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"altLoc1",     17, 17, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName1",    18, 20, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainID1",    22, 22, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resSeq1",     23, 26, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCode1",      27, 27, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"name2",       43, 46, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"altLoc2",     47, 47, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName2",    48, 50, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainID2",    52, 52, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resSeq2",     53, 56, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCode2",      57, 57, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"sym1",        60, 65, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"sym2",        67, 72, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
};

static struct PDBFieldDef remark_fields[] = {
  {"remarkNum",    8, 10, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"text",        12, 70, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
};

static struct PDBFieldDef end_fields[] = {
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
};

static struct PDBFieldDef sheet_fields[] = {
  {"strand",       8, 10, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"sheetID",     12, 14, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"numStrands",  15, 16, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"initResName", 18, 20, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"initChainID", 22, 22, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"initSeqNum",  23, 26, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"initICode",   27, 27, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"endResName",  29, 31, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"endChainID",  33, 33, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"endSeqNum",   34, 37, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"endICode",    38, 38, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"sense",       39, 40, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"curAtom",     42, 45, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"curResName",  46, 48, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"curChainID",  50, 50, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"curResSeq",   51, 54, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"curICode",    55, 55, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"prevAtom",    57, 60, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"prevResName", 61, 63, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"prevChainID", 65, 65, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"prevResSeq",  66, 69, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"prevICode",   70, 70, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
};

static struct PDBFieldDef scale2_fields[] = {
  {"s[n][1]",     11, 20, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"s[n][2]",     21, 30, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"s[n][3]",     31, 40, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"u[n]",        46, 55, PDB_FIELD_TYPE_FLOAT_5, PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_FLOAT_5, PDB_FIELD_LJUST},
};

static struct PDBFieldDef author_fields[] = {
  {"continuation", 9, 10, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"authorList",  11, 70, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
};

static struct PDBFieldDef cryst1_fields[] = {
  {"a",            7, 15, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"b",           16, 24, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"c",           25, 33, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"alpha",       34, 40, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"beta",        41, 47, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"gamma",       48, 54, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"sgroup",      56, 66, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"z",           67, 70, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
};

static struct PDBFieldDef siguij_fields[] = {
  {"serial",       7, 11, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"name",        13, 16, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {"altLoc",      17, 17, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName",     18, 20, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainID",     22, 22, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resSeq",      23, 26, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCode",       27, 27, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"sig[1][1]",   29, 35, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"sig[2][2]",   36, 42, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"sig[3][3]",   43, 49, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"sig[1][2]",   50, 56, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"sig[1][3]",   57, 63, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"sig[2][3]",   64, 70, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"segID",       73, 76, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"element",     77, 78, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"charge",      79, 80, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
};

static struct PDBFieldDef cryst3_fields[] = {
  {"a",            7, 15, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"b",           16, 24, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"c",           25, 33, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"alpha",       34, 40, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"beta",        41, 47, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"gamma",       48, 54, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"sgroup",      56, 66, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"z",           67, 70, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
};

static struct PDBFieldDef cryst2_fields[] = {
  {"a",            7, 15, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"b",           16, 24, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"c",           25, 33, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"alpha",       34, 40, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"beta",        41, 47, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"gamma",       48, 54, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"sgroup",      56, 66, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"z",           67, 70, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
};

static struct PDBFieldDef cispep_fields[] = {
  {"serial",       8, 10, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"resName1",    12, 14, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainID1",    16, 16, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"seqNum1",     18, 21, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCode1",      22, 22, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName2",    26, 28, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainID2",    30, 30, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"seqNum2",     32, 35, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCode2",      36, 36, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"modNum",      44, 46, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"measure",     54, 59, PDB_FIELD_TYPE_FLOAT_2, PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_FLOAT_2, PDB_FIELD_LJUST},
};

static struct PDBFieldDef atom_fields[] = {
  {"serial",       7, 11, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"name",        13, 16, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {"altLoc",      17, 17, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName",     18, 20, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainID",     22, 22, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resSeq",      23, 26, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCode",       27, 27, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"x",           31, 38, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"y",           39, 46, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"z",           47, 54, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"occupancy",   55, 60, PDB_FIELD_TYPE_FLOAT_2, PDB_FIELD_LJUST},
  {"tempFactor",  61, 66, PDB_FIELD_TYPE_FLOAT_2, PDB_FIELD_LJUST},
  {"segID",       73, 76, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"element",     77, 78, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"charge",      79, 80, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
};

static struct PDBFieldDef origx3_fields[] = {
  {"o[n][1]",     11, 20, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"o[n][2]",     21, 30, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"o[n][3]",     31, 40, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"t[n]",        46, 55, PDB_FIELD_TYPE_FLOAT_5, PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_FLOAT_5, PDB_FIELD_LJUST},
};

static struct PDBFieldDef origx2_fields[] = {
  {"o[n][1]",     11, 20, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"o[n][2]",     21, 30, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"o[n][3]",     31, 40, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"t[n]",        46, 55, PDB_FIELD_TYPE_FLOAT_5, PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_FLOAT_5, PDB_FIELD_LJUST},
};

static struct PDBFieldDef modres_fields[] = {
  {"idCode",       8, 11, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName",     13, 15, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainID",     17, 17, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"seqNum",      19, 22, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCode",       23, 23, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"stdRes",      25, 27, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"comment",     30, 70, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
};

static struct PDBFieldDef source_fields[] = {
  {"continuation", 9, 10, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"srcName",     11, 70, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
};

static struct PDBFieldDef formul_fields[] = {
  {"compNum",      9, 10, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"hetID",       13, 15, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"continuation",17, 18, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"asterisk",    19, 19, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"text",        20, 70, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
};

static struct PDBFieldDef master_fields[] = {
  {"numRemark",   11, 15, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"O",           16, 20, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"numHet",      21, 25, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"numHelix",    26, 30, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"numSheet",    31, 35, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"numTurn",     36, 40, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"numSite",     41, 45, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"numXForm",    46, 50, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"numCoord",    51, 55, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"numTer",      56, 60, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"numConect",   61, 65, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"numSeq",      66, 70, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
};

static struct PDBFieldDef caveat_fields[] = {
  {"continuation", 9, 10, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"idCode",      12, 15, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"comment",     20, 70, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
};

static struct PDBFieldDef het_fields[] = {
  {"hetID",        8, 10, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainID",     13, 13, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"seqNum",      14, 17, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCode",       18, 18, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"numHetAtoms", 21, 25, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"text",        31, 70, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
};

static struct PDBFieldDef compnd_fields[] = {
  {"continuation", 9, 10, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"compound",    11, 70, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
};

static struct PDBFieldDef model_fields[] = {
  {"serial",      11, 14, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
};

static struct PDBFieldDef revdat_fields[] = {
  {"modNum",       8, 10, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"continuation",11, 12, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"modDate",     14, 22, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"modID",       24, 28, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"modType",     32, 32, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"record1",     40, 45, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {"record2",     47, 52, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {"record3",     54, 59, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {"record4",     61, 66, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
};

static struct PDBFieldDef ssbond_fields[] = {
  {"serNum",       8, 10, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"resName1",    12, 14, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainID1",    16, 16, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"seqNum1",     18, 21, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCode1",      22, 22, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName2",    26, 28, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainID2",    30, 30, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"seqNum2",     32, 35, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCode2",      36, 36, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"sym1",        60, 65, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"sym2",        67, 72, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
};

static struct PDBFieldDef obslte_fields[] = {
  {"continuation", 9, 10, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"repDate",     12, 20, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"idCode",      22, 25, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"rIdCode1",    32, 35, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"rIdCode2",    37, 40, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"rIdCode3",    42, 45, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"rIdCode4",    47, 50, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"rIdCode5",    52, 55, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"rIdCode6",    57, 60, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"rIdCode7",    62, 65, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"rIdCode8",    67, 70, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
};

static struct PDBFieldDef conect_fields[] = {
  {"serial",       7, 11, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"serialBond1", 12, 16, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"serialBond2", 17, 21, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"serialBond3", 22, 26, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"serialBond4", 27, 31, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"serialHydBond1",32, 36, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"serialHydBond2",37, 41, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"serialSaltBond1",42, 46, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"serialHydBond3",47, 51, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"serialHydBond4",52, 56, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"serialSaltBond2",57, 61, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
};

static struct PDBFieldDef jrnl_fields[] = {
  {"text",        13, 70, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
};

static struct PDBFieldDef hetnam_fields[] = {
  {"continuation", 9, 10, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {"hetID",       12, 14, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"text",        16, 70, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
};

static struct PDBFieldDef sprsde_fields[] = {
  {"continuation", 9, 10, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"sprsdeDate",  12, 20, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"idCode",      22, 25, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"sIdCode1",    32, 35, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"sIdCode2",    37, 40, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"sIdCode3",    42, 45, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"sIdCode4",    47, 50, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"sIdCode5",    52, 55, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"sIdCode6",    57, 60, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"sIdCode7",    62, 65, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"sIdCode8",    67, 70, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
};

static struct PDBFieldDef header_fields[] = {
  {"classification",11, 50, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"depDate",     51, 59, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"idCode",      63, 66, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
};

static struct PDBFieldDef site_fields[] = {
  {"seqNum",       8, 10, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"siteID",      12, 14, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"numRes",      16, 17, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"resName1",    19, 21, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainID1",    23, 23, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"seq1",        24, 27, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCode1",      28, 28, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName2",    30, 32, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainID2",    34, 34, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"seq2",        35, 38, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCode2",      39, 39, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName3",    41, 43, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainID3",    45, 45, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"seq3",        46, 49, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCode3",      50, 50, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName4",    52, 54, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainID4",    56, 56, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"seq4",        57, 60, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCode4",      61, 61, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
};

static struct PDBFieldDef scale1_fields[] = {
  {"s[n][1]",     11, 20, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"s[n][2]",     21, 30, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"s[n][3]",     31, 40, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"u[n]",        46, 55, PDB_FIELD_TYPE_FLOAT_5, PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_FLOAT_5, PDB_FIELD_LJUST},
};

static struct PDBFieldDef hydbnd_fields[] = {
  {"name1",       13, 16, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"altLoc1",     17, 17, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName1",    18, 20, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainID1",    22, 22, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resSeq1",     23, 27, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCode1",      28, 28, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"nameH",       30, 33, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"altLocH",     34, 34, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainH",      36, 36, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resSeqH",     37, 41, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCodeH",      42, 42, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"name2",       44, 47, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"altLoc2",     48, 48, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName2",    49, 51, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainID2",    53, 53, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resSeq2",     54, 58, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCode2",      59, 59, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"sym1",        60, 65, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"sym2",        67, 72, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
};

static struct PDBFieldDef ter_fields[] = {
  {"serial",       7, 11, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"resName",     18, 20, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainID",     22, 22, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resSeq",      23, 26, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCode",       27, 27, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
};

static struct PDBFieldDef dbref_fields[] = {
  {"idCode",       8, 11, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chain_ID",    13, 13, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"seqBegin",    15, 18, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"insertBegin", 19, 19, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"seqEnd",      21, 24, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"insertEnd",   25, 25, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"database",    27, 32, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {"dbAccession", 34, 41, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {"dbIdCode",    43, 54, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {"dbseqBegin",  56, 60, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"idbnsBeg",    61, 61, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"dbseqEnd",    63, 67, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"dbinsEnd",    68, 68, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
};

static struct PDBFieldDef expdta_fields[] = {
  {"continuation", 9, 10, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"technique",   11, 70, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
};

static struct PDBFieldDef scale3_fields[] = {
  {"s[n][1]",     11, 20, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"s[n][2]",     21, 30, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"s[n][3]",     31, 40, PDB_FIELD_TYPE_FLOAT_6, PDB_FIELD_LJUST},
  {"u[n]",        46, 55, PDB_FIELD_TYPE_FLOAT_5, PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_FLOAT_5, PDB_FIELD_LJUST},
};

static struct PDBFieldDef hetatm_fields[] = {
  {"serial",       7, 11, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"name",        13, 16, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {"altLoc",      17, 17, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName",     18, 20, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainID",     22, 22, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resSeq",      23, 26, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCode",       27, 27, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"x",           31, 38, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"y",           39, 46, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"z",           47, 54, PDB_FIELD_TYPE_FLOAT_3, PDB_FIELD_LJUST},
  {"occupancy",   55, 60, PDB_FIELD_TYPE_FLOAT_2, PDB_FIELD_LJUST},
  {"tempFactor",  61, 66, PDB_FIELD_TYPE_FLOAT_2, PDB_FIELD_LJUST},
  {"segID",       73, 76, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"element",     77, 78, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"charge",      79, 80, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
};

static struct PDBFieldDef seqadv_fields[] = {
  {"idCode",       8, 11, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"resName",     13, 15, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"chainID",     17, 17, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"seqNum",      19, 22, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"iCode",       23, 23, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"database",    25, 28, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {"dbIDCode",    30, 38, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {"dbRes",       40, 42, PDB_FIELD_TYPE_STRING,  PDB_FIELD_LJUST},
  {"dbSeq",       44, 48, PDB_FIELD_TYPE_INTEGER, PDB_FIELD_LJUST},
  {"convlict",    40, 70, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
  {NULL,           0,  0, PDB_FIELD_TYPE_STRING,  PDB_FIELD_RJUST},
};

static struct PDBRecordDef {
  char               *name;
  struct PDBFieldDef *fields;
} pdb_record_defs[] = {
  {"ATOM  ", atom_fields},
  {"ANISOU", anisou_fields},
  {"REMARK", remark_fields},
  {"KEYWDS", keywds_fields},
  {"SEQRES", seqres_fields},
  {"LINK  ", link_fields},
  {"ORIGX1", origx1_fields},
  {"ENDMDL", endmdl_fields},
  {"HETSYN", hetsyn_fields},
  {"SIGATM", sigatm_fields},
  {"MTRIX2", mtrix2_fields},
  {"MTRIX3", mtrix3_fields},
  {"HELIX ", helix_fields},
  {"MTRIX1", mtrix1_fields},
  {"TURN  ", turn_fields},
  {"TITLE ", title_fields},
  {"SLTBRG", sltbrg_fields},
  {"END   ", end_fields},
  {"SHEET ", sheet_fields},
  {"SCALE2", scale2_fields},
  {"AUTHOR", author_fields},
  {"CRYST1", cryst1_fields},
  {"SIGUIJ", siguij_fields},
  {"CRYST3", cryst3_fields},
  {"CRYST2", cryst2_fields},
  {"CISPEP", cispep_fields},
  {"ORIGX3", origx3_fields},
  {"ORIGX2", origx2_fields},
  {"MODRES", modres_fields},
  {"SOURCE", source_fields},
  {"FORMUL", formul_fields},
  {"MASTER", master_fields},
  {"CAVEAT", caveat_fields},
  {"HET   ", het_fields},
  {"COMPND", compnd_fields},
  {"MODEL ", model_fields},
  {"REVDAT", revdat_fields},
  {"SSBOND", ssbond_fields},
  {"OBSLTE", obslte_fields},
  {"CONECT", conect_fields},
  {"JRNL  ", jrnl_fields},
  {"HETNAM", hetnam_fields},
  {"SPRSDE", sprsde_fields},
  {"HEADER", header_fields},
  {"SITE  ", site_fields},
  {"SCALE1", scale1_fields},
  {"HYDBND", hydbnd_fields},
  {"TER   ", ter_fields},
  {"DBREF ", dbref_fields},
  {"EXPDTA", expdta_fields},
  {"SCALE3", scale3_fields},
  {"HETATM", hetatm_fields},
  {"SEQADV", seqadv_fields},
  {NULL,     NULL}
};

static PyObject *
pdb_read(PyObject *self, PyObject *args)
{
  int       i;
  int       n;
  int       irec;
  int       ifie;
  int       ibuf;
  char     *pdb_path;
  char      in_buff[MAX_LINE];
  char      field[MAX_LINE];
  FILE     *pdb_fil;
  PyObject *py_pdb_list = NULL;
  PyObject *py_rec_dict = NULL;
  PyObject *py_strx     = NULL;
  PyObject *py_intx     = NULL;
  PyObject *py_floatx   = NULL;


  if (!PyArg_ParseTuple(args, "s", &pdb_path))
    return NULL;
  
  /* open input/output files */
  if ((pdb_fil = fopen(pdb_path, "r")) == NULL) {
    Py_INCREF(Py_None);
    return Py_None;
  }

  py_pdb_list = PyList_New(0);

  while (fgets(in_buff, MAX_LINE, pdb_fil) != NULL) {
    
    /* strip newlines */
    ibuf = strnlen(in_buff, MAX_LINE) - 1;
    while (in_buff[ibuf] == '\n' && ibuf >= 0) {
      in_buff[ibuf] = '\0';
      ibuf--;
    }
    if (ibuf < 0) {
      continue;
    }

    if (ibuf < 5) {
      n = ibuf + 1;
    } else {
      n = 6;
    }
    
    /* find the index in  */
    for (irec = 0; pdb_record_defs[irec].name != NULL; irec++) {      
      if (strncmp(pdb_record_defs[irec].name, in_buff, n) == 0) {
	break;
      }
    }
    if (pdb_record_defs[irec].name == NULL) {
      continue;
    }

    /* create the new python dictionary for this record */
    py_rec_dict = PyDict_New();
    PyList_Append(py_pdb_list, py_rec_dict);
  
    /* set the record name */
    py_strx = PyString_FromString(pdb_record_defs[irec].name);
    PyDict_SetItemString(py_rec_dict, "RECORD", py_strx);
    Py_XDECREF(py_strx);
    
    for (ifie = 0;
	 pdb_record_defs[irec].fields[ifie].name != NULL;
	 ifie++) {
	  
      /* beginning index of the field */
      i = pdb_record_defs[irec].fields[ifie].istart - 1;
      /* length of the field */
      n = pdb_record_defs[irec].fields[ifie].iend   - i;
      
      /* stop reading at the end of the buffer */
      if (i > ibuf) {
	break;
      }      

      /* reduce the length of the field if the field definition
       * extends beyond the end of the in_buff
       */
      if (i+n-1 > ibuf) {
	n = ibuf - i + 1;
      }

      switch (pdb_record_defs[irec].fields[ifie].type) {

      case PDB_FIELD_TYPE_STRING:
	strncpy(field, &in_buff[i], n);
	field[n] = '\0';

	/* strip right hand whitespace */
	i = strnlen(field, MAX_LINE) - 1;
	while (field[i] == ' ' && i >= 0) {
	  field[i] = '\0';
	  i--;
	}

	/* no blanks */
	if (i >= 0) {
	  py_strx = PyString_FromString(field);

	  PyDict_SetItemString(
	    py_rec_dict,
	    pdb_record_defs[irec].fields[ifie].name,
	    py_strx);

	  Py_XDECREF(py_strx);
	}
	break;

      case PDB_FIELD_TYPE_INTEGER:
	strncpy(field, &in_buff[i], n);
	field[n] = '\0';

	py_intx = PyInt_FromString(field, NULL, 10);
	if (py_intx != NULL) {
	  PyDict_SetItemString(
	    py_rec_dict,
	    pdb_record_defs[irec].fields[ifie].name,
	    py_intx);

	  Py_XDECREF(py_intx);
	} else {
	  PyErr_Clear();
	}
	break;

      case PDB_FIELD_TYPE_FLOAT_2:
      case PDB_FIELD_TYPE_FLOAT_3:
      case PDB_FIELD_TYPE_FLOAT_4:
      case PDB_FIELD_TYPE_FLOAT_5:
      case PDB_FIELD_TYPE_FLOAT_6:
	py_strx   = PyString_FromStringAndSize(&in_buff[i], n);
	py_floatx = PyFloat_FromString(py_strx, NULL);
	Py_XDECREF(py_strx);

	if (py_floatx != NULL) {
	  PyDict_SetItemString(
	    py_rec_dict,
	    pdb_record_defs[irec].fields[ifie].name,
	    py_floatx);

	  Py_XDECREF(py_floatx);
	} else {
	  PyErr_Clear();
	}
	break;
      } /* switch */
    } /* for (fields) */

    Py_XDECREF(py_rec_dict);
  } /* while (records) */

  /* close file */
  fclose(pdb_fil);
  return py_pdb_list;
}


static PyMethodDef PDBModuleMethods[] = {
  {"read",
   pdb_read,
   METH_VARARGS,
   "Reads a PDB file and returns a list."},
  
  {NULL, NULL, 0, NULL}
};


DL_EXPORT(void)
initpdbmodule(void)
{
  PyObject *m;
  
  m = Py_InitModule("pdbmodule", PDBModuleMethods);
  
  PDBModuleErr = PyErr_NewException("pdbmodule.error", NULL, NULL);
  Py_INCREF(PDBModuleErr);
  PyModule_AddObject(m, "error", PDBModuleErr);
}
