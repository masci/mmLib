#!/home/tlsmd/local/bin/python
# -*- Mode: Python -*-
import sys, os
PYMMLIB_ROOT = "/home/tlsmd/pymmlib"
TLSMD_ROOT = "/home/tlsmd/tlsmd"
MMOV_DIR = os.path.join(TLSMD_ROOT, "molmovdb")

sys.path.insert(0, PYMMLIB_ROOT)
sys.path.insert(0, os.path.join(TLSMD_ROOT, "bin"))
sys.path.insert(0, MMOV_DIR)
os.chdir(MMOV_DIR)

import molmovdb_cgi
molmovdb_cgi.main()
sys.exit(0)
