## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

NucleicAcidNames = [
    "A", "C", "G", "T", "U", "I", "+A", "+C", "+G", "+T", "+U", "+I"]


class NucleicAcid:
    """Empty definition class for building a Nucleic Acid library."""

    def __init__(self):
        self.name = ""
        self.full_name = ""
        self.one_letter_name = ""


    def __str__(self):
        return "NucleicAcid=%s" % (self.name)



NucleicAcidMap = {
    }
