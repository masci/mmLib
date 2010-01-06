## TLS Motion Determination (TLSMD)
## Copyright 2002-2010 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

## Python modules
import time
import datetime
import random
import string
import re
import subprocess

## TLSMD
import console
import conf

_STIME = 0.0

def start_timing():
    global _STIME
    _STIME = time.time()

def end_timing():
    global _STIME
    tm = time.time() - _STIME
    return "Computation Time: %5.2f sec" % (tm)

def timestamp():
    """Generate a timestamp at the moment this function is called. This is the
    main/global timestamp generator for all TLSMD scripts.
    """
    ## TODO: Allow for either "-7" or "-10". 2009-03-25
    ## Also, allow for "Y-M-D" format, etc.
    return datetime.datetime.fromtimestamp(time.time()).isoformat(' ')[:-7]

def start_time():
    ## returns format e.g., "26 May 2009"
    return time.strftime("%d %b %Y", time.localtime(conf.globalconf.start_time))

def begin_chain_timing(chain_id):
    console.stdoutln("BEGIN TIMING CHAIN %s %f" % (chain_id, time.time()))

def end_chain_timing(chain_id):
    console.stdoutln("END TIMING CHAIN %s %f" % (chain_id, time.time()))

def rgb_f2i(rgb):
    """Transforms the float 0.0-1.0 RGB color values to
    integer 0-255 RGB values.
    """
    r, g, b = rgb
    ri = int(255.0 * r)
    gi = int(255.0 * g)
    bi = int(255.0 * b)
    return (ri, gi, bi)

def rgb_f2s(rgbf):
    rgbs = "#%2x%2x%2x" % rgb_f2i(rgbf)
    rgbs = rgbs.replace(" ", "0")
    return rgbs

def generate_security_code(code_length = 8):
    """Generates a random 8-character string
    """
    random.seed()
    codelist = list(5 * string.ascii_letters)
    random.shuffle(codelist)
    code = "".join(random.sample(codelist, code_length))
    return code

def parse_chains(chains):
    """Parses a given chains string returned from the MySQL DB
    """
    ## Turns "A:100:1:aa" -> "A", "100", "1", "aa"
    ## I.e.: chain_id = A, num_residues = 100,
    ##       selected = 1/True, type = aa/amino acids
    ## NOTE:
    ##   - aa = amino acid
    ##   - na = nucleic acid
    ##   - ot = other
    chid     = re.sub(r'([A-Za-z0-9]):([0-9]{1,}):([01]):([naot]{2});?', '\\1;', chains).rstrip(";")
    length   = re.sub(r'([A-Za-z0-9]):([0-9]{1,}):([01]):([naot]{2});?', '\\2;', chains).rstrip(";")
    selected = re.sub(r'([A-Za-z0-9]):([0-9]{1,}):([01]):([naot]{2});?', '\\3;', chains).rstrip(";")
    type     = re.sub(r'([A-Za-z0-9]):([0-9]{1,}):([01]):([naot]{2});?', '\\4;', chains).rstrip(";")

    return chid, length, selected, type

def parse_molauto(infile, outfile):
    """Parses the molauto output to force each chain to have its own unique
    colour.
    """
    file = open(outfile, "w")
    for line in open(infile).readlines():
        file.write("%s" % line)
        if(re.match(r'^  set segments', line)):
            file.write("  set segments 10;\n")
            file.write("  set planecolour hsb 0.6667 1 1;")
        elif(re.match(r'^  set planecolour', line)):
            colour = line
        elif(re.match(r'^  .* from ', line)):
            chain1 = chain2 = line
            chain1 = re.sub(r'^  .* from ([A-Z])[0-9]{1,} to ([A-Z])[0-9].*$', '\\1', line).strip()
            chain2 = re.sub(r'^  .* from ([A-Z])[0-9]{1,} to ([A-Z])[0-9].*$', '\\2', line).strip()
            file.write("%s" % line)
            if(chain1 != chain2):
                file.write("%s" % colour)
        else:
            file.write("%s" % line)
    file.close()
    return

def run_subprocess(cmdlist):
    """General call to 'subprocess' for a given command list.
    """
    proc = subprocess.Popen(cmdlist,
                            shell = True,
                            stdin = subprocess.PIPE,
                            stdout = subprocess.PIPE,
                            stderr = subprocess.PIPE,
                            close_fds = True,
                            bufsize = 32768)

def render_struct(job_dir):
    """Generate struct.png via molauto/parse_molauto/molscript/render
    """
    ## cmd: molauto smallAB.pdb|parse_molauto.pl|molscript -r |\
    ##      render -bg white -size 200x200 -png mymol.png
    cmdlist = ["%s %s/struct.pdb | %s | %s -r | %s -bg white -size %s -png %s/struct.png 1>&2" % (
              conf.MOLAUTO, job_dir, conf.PARSE_MOLAUTO_SCRIPT,
              conf.MOLSCRIPT, conf.RENDER,
              conf.RENDER_SIZE, job_dir)]
    run_subprocess(cmdlist)

def extract_raw_backbone(infile, outfile):
    """Extracts all backbone atoms (amino acid or nucleic acid)
    for use in tlsanim2r3d
    """
    file = open(outfile, "w")
    for line in open(infile).readlines():
        if(re.match(r'^ATOM........ CA .*', line) or \
           re.match(r'^ATOM........ P  .*', line) or \
           re.match(r'^ATOM........ C[543]\'.*', line) or \
           re.match(r'^ATOM........ O[53]\'.*', line)):
            chain_id = line[21:22]
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())

            file.write("1 0 %s 0 0 %.3f %.3f %.3f\n" % (chain_id, x, y, z))
    file.close()
    return

def generate_raw_grey_struct(job_dir):
    """Generate 'raw' input for tlsanim2r3d, but only for the non-animated
    sections, for _all_ chains.
    """
    ## cmd: ./extract_raw_chains.pl <smallAB.pdb |./tlsanim2r3d - >ANALYSIS/struct.r3d
    extract_raw_backbone("%s/struct.pdb" % job_dir, "%s/struct.raw" % job_dir)

    ## cmd: ./tlsanim2r3d < struct.raw >ANALYSIS/struct.r3d
    cmdlist = ["%s < %s/struct.raw > %s/struct.r3d 2>/dev/null" % (
              conf.TLSANIM2R3D, job_dir, job_dir)]
    run_subprocess(cmdlist)

def generate_bases_r3d(job_dir, chain_id):
    """Generate 'raw' input for tlsanim2r3d, but only for the non-animated
    sections, for _all_ chains.
    """
    ## cmd: grep '^ATOM.................B.*' | rings3d -bases >>bases.r3d
    cmdlist = ["grep '^ATOM.................%s.*' %s/struct.pdb | %s -bases >>%s/bases.r3d 2>/dev/null" % (
              chain_id, job_dir, conf.RINGS3D, job_dir)]
    run_subprocess(cmdlist)

def generate_sugars_r3d(job_dir, chain_id):
    """Generate 'raw' input for tlsanim2r3d, but only for the non-animated
    sections, for _all_ chains.
    """
    ## cmd: grep '^ATOM.................B.*' | rings3d -ribose >>sugars.r3d
    cmdlist = ["grep '^ATOM.................%s.*' %s/struct.pdb | %s -ribose >>%s/sugars.r3d 2>/dev/null" % (
              chain_id, job_dir, conf.RINGS3D, job_dir)]
    run_subprocess(cmdlist)
