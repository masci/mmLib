#!/usr/bin/env python
## Copyright 2002-2010 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

## Python
import os
import sys
import popen2
import gc

## pymmlib
import kid
import test_util
from mmLib import mmCIF, FileIO

NOTFOUND = "########"


class mmCIFRowKid(object):
    def __init__(self, cif_row = None):
        self.cif_row = cif_row

    def __getattr__(self, name):
        try:
            return self.cif_row[name]
        except (KeyError, TypeError):
            return NOTFOUND

    def __getitem__(self, x):
        try:
            return self.cif_row[x]
        except (KeyError, TypeError):
            return NOTFOUND


class mmCIFTableKid(object):
    def __init__(self, cif_table = None):
        self.cif_table = cif_table

    def __getattr__(self, name):
        try:
            return self.cif_table[name]
        except (KeyError, TypeError):
            return NOTFOUND

    def __getitem__(self, x):
        try:
            return self.cif_table[x]
        except (KeyError, TypeError):
            return NOTFOUND

    def __iter__(self):
        if self.cif_table is not None:
            for row in self.cif_table:
                yield mmCIFRowKid(row)


class mmCIFDataKid(object):
    def __init__(self, cif_data = None):
        self.cif_data = cif_data

    def __getattr__(self, name):
        try:
            return mmCIFTableKid(self.cif_data[name])
        except (KeyError, TypeError):
            return mmCIFTableKid()

    def __getitem__(self, x):
        try:
            return self.cif_data[x]
        except (KeyError, TypeError):
            return mmCIFTableKid()


def cif2html(cif_path, html_path):
    fileobj = FileIO.OpenFile(cif_path, "r")

    cif_file = mmCIF.mmCIFFile()
    print "loading..."
    cif_file.load_file(fileobj)
    print "converting..."
    cif_data = cif_file[0]

    c2h_template = kid.Template(file = "cif2html.kid")
    c2h_template.cif = mmCIFDataKid(cif_data)
    c2h_template.write(html_path)


def path_split2(base_path, path):
    full_path, filename = os.path.split(path)
    return full_path[len(base_path) + 1:], filename


def main():
    try:
        cif_root = sys.argv[1]
        html_root = sys.argv[2]
    except IndexError:
        sys.stderr.write("cif2html.py <cif_dir> <html_dir>\n\n")
        raise SystemExit

    for cif_path in test_util.walk_cif(cif_root):
        relative_path, filename = path_split2(cif_root, cif_path)

        html_dir = os.path.join(html_root, relative_path)
        if not os.path.isdir(html_dir):
            os.makedirs(html_dir)
        struct_id = filename[:4]
        html_path = os.path.join(html_dir, "%s.xhtml" % (struct_id))

        if os.path.isfile(cif_path):
            if not os.path.isfile(html_path):
                print "%s -> %s" % (cif_path, html_path)
                cif2html(cif_path, html_path)
                gc.collect()
            else:
                print "*%s -> %s" % (cif_path, html_path)


if __name__ == "__main__":
    main()
