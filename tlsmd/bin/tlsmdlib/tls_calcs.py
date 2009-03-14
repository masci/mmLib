## TLS Minimized Domains (TLSMD)
## Copyright 2002-2009 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

## Python
import math
import numpy
import console
import re ## to force residue numbers to be integers

## pymmlib
from mmLib import Constants, TLS, FileIO


def create_fractional_residue_number(res_num):
    """Converts insertion residues to fractional residue numbers.
       E.g., "5A" -> "5.0"
       This is so gnuplot can handle x-axis number values.
    """
    ## TODO: Figure out a better way to handle insertion residues, 2008-12-03

    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" ## used for index position

    if re.sub(r'.?([A-Za-z]?)', '\\1', res_num) != '':
        ## seg_start has insertion residues (e.g., "5A")
        fraction = alphabet.index(re.sub(r'.?([A-Za-z]?)', '\\1', res_num.upper()))

        ## remove any non-integer chars
        res_num = re.sub(r'[^0-9]', '', res_num)

        ## create fractional residue number (e.g., "5.0")
        res_num = res_num + "." + str(fraction)

    return res_num

def tlsdict2tensors(tlsdict):
    """Convert the result dictionaries returned by tlsmdmodule to NumPy
    tensors.
    """
    origin = numpy.array([tlsdict["x"], tlsdict["y"], tlsdict["z"]], float)

    T = numpy.array(
        [ [tlsdict["t11"], tlsdict["t12"], tlsdict["t13"]],
          [tlsdict["t12"], tlsdict["t22"], tlsdict["t23"]],
          [tlsdict["t13"], tlsdict["t23"], tlsdict["t33"]] ], float)
    
    L = numpy.array(
        [ [tlsdict["l11"], tlsdict["l12"], tlsdict["l13"]],
          [tlsdict["l12"], tlsdict["l22"], tlsdict["l23"]],
          [tlsdict["l13"], tlsdict["l23"], tlsdict["l33"]] ], float)
    
    s11, s22, s33 = TLS.calc_s11_s22_s33(tlsdict["s2211"], tlsdict["s1133"]) 
        
    S = numpy.array(
        [ [       s11, tlsdict["s12"], tlsdict["s13"]],
          [tlsdict["s21"],        s22, tlsdict["s23"]],
          [tlsdict["s31"], tlsdict["s32"],       s33] ], float)

    return T, L, S, origin

def isotlsdict2tensors(itlsdict):
    """Convert the result dictionaries returned by tlsmdmodule to NumPy
    tensors.
    """
    origin = numpy.array([itlsdict["x"], itlsdict["y"], itlsdict["z"]], float)
    
    IT = itlsdict["it"]
        
    IL = numpy.array(
        [ [itlsdict["il11"], itlsdict["il12"], itlsdict["il13"]],
          [itlsdict["il12"], itlsdict["il22"], itlsdict["il23"]],
          [itlsdict["il13"], itlsdict["il23"], itlsdict["il33"]] ], float)
    
    IS = numpy.array([itlsdict["is1"], itlsdict["is2"], itlsdict["is3"]], float)

    return IT, IL, IS, origin 

class ITLSParameters(object):
    """Not used yet.
    """
    def __init__(self, arg):
        if isinstance(arg, ITLSParameters):
            self.IT = arg.IT.copy()
            self.IL = arg.IL.copy()
            self.IS = arg.IS.copy()
            self.IO = arg.IO.copy()
        elif isinstance(arg, dict):
            self.IT, self.IL, self.IS, self.IO = isotlsdict2tensors(arg)
        else:
            raise ValueError

class TLSParameters(object):
    """Not used yet.
    """
    def __init__(self, arg):
        if isinstance(arg, TLSParameters):
            self.T = arg.T.copy()
            self.L = arg.L.copy()
            self.S = arg.S.copy()
            self.O = arg.O.copy()
        elif isinstance(arg, dict):
            self.T, self.L, self.S, self.O = tlsdict2tensors(arg)
        else:
            raise ValueError

def calc_rmsd_tls_biso(tls_group):
    """Calculate the RMSD of the tls_group using the isotropic TLS model.
    """
    T = tls_group.itls_T
    L = tls_group.itls_L
    S = tls_group.itls_S
    O = tls_group.origin
    
    msd_sum = 0.0
    
    for atm, uiso_tls in TLS.iter_itls_uiso(iter(tls_group), T, L, S, O):
        msd_sum += (Constants.U2B*uiso_tls - atm.temp_factor)**2
        
    if len(tls_group)>0:
        msd = msd_sum / len(tls_group)
        rmsd = math.sqrt(msd)
    else:
        rmsd = 0.0

    return rmsd

def calc_mean_biso_obs(chain):
    """Calculates the mean B value per residue in the chain (as observed in the input structure).
    """
    num_res = chain.count_fragments()
    biso = numpy.zeros(num_res, float)

    for frag in chain.iter_fragments():
        n = 0
        b_sum_obs = 0.0

        for atm in frag.iter_all_atoms():
            if atm.include == False:
                continue

            n += 1
            b_sum_obs += atm.temp_factor

        if n > 0:
            biso[frag.ifrag] = b_sum_obs / n

    return biso


def calc_mean_biso_tls(chain, cpartition):
    """Calculated the mean B value per residue in the chain
    as calculated in the chain optimization.
    """
    num_res = chain.count_fragments()
    biso = numpy.zeros(num_res, float)

    for i, tls in enumerate(cpartition.iter_tls_segments()):
        tls_group = tls.tls_group

        T = tls_group.itls_T
        L = tls_group.itls_L
        S = tls_group.itls_S
        O = tls_group.origin

        for frag in tls.iter_fragments():
            n = 0
            b_sum_tls = 0.0

            for atm in frag.iter_all_atoms():
                if atm.include is False:
                    continue

                n += 1
                b_sum_tls += Constants.U2B * TLS.calc_itls_uiso(T, L, S, atm.position - O)

            if n > 0:
                biso[frag.ifrag] = b_sum_tls / n

    return biso


def calc_residue_mean_rmsd(chain, cpartition):
    num_tls = cpartition.num_tls_segments()
    num_res = chain.count_fragments()

    cmtx = numpy.zeros((num_tls, num_res), float)

    i_ntls = 0
    for i, tls in enumerate(cpartition.iter_tls_segments()):
        tls_group = tls.tls_group

        T = tls_group.itls_T # float(3)
        L = tls_group.itls_L # array(3,3)
        S = tls_group.itls_S # array(3): S[0], S[1], S[2]
        O = tls_group.origin # array(3)

        ##======================================================================
        ##<TLS_FILE>
        ##TLS 3
        ##CHAIN D
        ##RANGE 47-103
        ##ORIGIN -3.806092 0.302217 8.665245
        ##T 0.147191
        ##L 0.001070 0.000058 0.000790 0.000007 0.000229 -0.000206
        ##S -0.004597 0.001594 -0.002114
        ##
        ## note: S1 == S21-S12; S2 == S13-S31; S3 == S32-S23
        ## u_tls = T + (
        ##         L[0,0]*(zz+yy) + L[1,1]*(xx+zz) + L[2,2]*(xx+yy)
        ##         - 2.0*L[0,1]*x*y - 2.0*L[0,2]*x*z - 2.0*L[1,2]*y*z
        ##         + 2.0*S[0]*z + 2.0*S[1]*y + 2.0*S[2]*x) / 3.0
        ##
        O1, O2, O3 = O
        ## tls = "A:1-10" || "G:5A-14" || "G:-1--10"
        chain_id  = re.sub(r'^([A-Za-z0-9]):(-?[A-Za-z0-9]{1,})-(-?[A-Za-z0-9]{1,})', '\\1', str(tls))
        seg_start = re.sub(r'^([A-Za-z0-9]):(-?[A-Za-z0-9]{1,})-(-?[A-Za-z0-9]{1,})', '\\2', str(tls))
        seg_end   = re.sub(r'^([A-Za-z0-9]):(-?[A-Za-z0-9]{1,})-(-?[A-Za-z0-9]{1,})', '\\3', str(tls))

        seg_start = create_fractional_residue_number(seg_start)
        seg_end   = create_fractional_residue_number(seg_end)

        filename = "CHAIN%s_NTLS%s-%s.tls_model" % (chain_id, num_tls, i+1)
        tls_file = open(filename, "w")
        console.stdoutln("TLS_MODEL: Saving %s" % filename) ## LOGFILE
        tls_file.write("TLS %s SEG %s\nCHAIN %s\nRANGE %s-%s\nORIGIN %f %f %f\n" % (
                       num_tls, i_ntls + 1, chain_id, seg_start, seg_end,
                       O1, O2, O3))
        tls_file.write("T %f\nL %f %f %f %f %f %f\nS %f %f %f\n" % (
                       T, L[0,0], L[1,1], L[2,2], L[0,1], L[0,2], L[1,2], 
                       S[0], S[1], S[2]))
        tls_file.close()
        ##</TLS_FILE>
        ##======================================================================

        filename = "CHAIN%s_NTLS%s-%s.b_obs" % (chain_id, num_tls, i+1)
        b_obs_file = open(filename, "w")
        console.stdoutln("B_OBS: Saving %s" % filename) ## LOGFILE

        for j, frag in enumerate(chain):
            ## NOTE: j = res_num, frag = Res(ALA,23,A)

            ## calculate a atom-normalized rmsd deviation for each residue
            num_atoms = 0
            msd_sum = 0.0

            for atm in frag.iter_all_atoms():
                if atm.include == False:
                    continue

                num_atoms += 1

                ##==============================================================
                ##<DEBUG>
                ## Equation:
                ##u_tls = T + (
                ##        L[0,0]*(zz+yy) + L[1,1]*(xx+zz) + L[2,2]*(xx+yy)
                ##        - 2.0*L[0,1]*x*y - 2.0*L[0,2]*x*z - 2.0*L[1,2]*y*z
                ##        + 2.0*S[0]*z + 2.0*S[1]*y + 2.0*S[2]*x) / 3.0
                ##    0   1   2
                ## 0 0,0 0,1 0,2
                ## 1 1,0 1,1 1,2
                ## 2 2,0 2,1 2,2
                ##</DEBUG>

                b_iso_tls = Constants.U2B * TLS.calc_itls_uiso(T, L, S, atm.position - O)
                delta = atm.temp_factor - b_iso_tls
                msd_sum += delta**2

                ##==============================================================
                ##<B_OBS_DATA>
                ## frag = "Res(VAL,39,A)" || "Res(G,5A,A)" || "Res(GLY,-1,A)"
                x, y, z = atm.position

                res_name = re.sub(r'^Res\(([A-Za-z]{1,3}),([A-Za-z0-9-]{1,}),[A-Za-z0-9]\)', '\\1', str(frag)).upper()
                res_num  = re.sub(r'^Res\(([A-Za-z]{1,3}),([A-Za-z0-9-]{1,}),[A-Za-z0-9]\)', '\\2', str(frag)).upper()
                res_num  = create_fractional_residue_number(res_num)

                if num_tls > 0 and (
                   (float(res_num) >= float(seg_start)) and 
                   (float(res_num) <= float(seg_end)) ):
                    b_obs_file.write("TLS:%s %s,%s,%s %.3f %.3f %.3f %.2f %.2f %.2f\n" % (
                        num_tls, res_name, res_num, chain_id, x, y, z,
                        atm.occupancy, atm.temp_factor, b_iso_tls))
                ##</B_OBS_DATA>
                ##==============================================================

            if num_atoms > 0:
                msd = msd_sum / num_atoms
                rmsd = math.sqrt(msd)

                ## set the cross prediction matrix
                cmtx[i,j] = rmsd

        b_obs_file.close()

    return cmtx
     
def refmac5_prep(xyzin, tlsin_list, xyzout, tlsout):
    """Use TLS model + Uiso for each atom.  Output xyzout with the
    residual Uiso only.
    """
    ## load structure
    struct = FileIO.LoadStructure(fil = xyzin)

    ## load and construct TLS groups
    tls_group_list = []
    tls_file = TLS.TLSFile()
    tls_file.set_file_format(TLS.TLSFileFormatTLSOUT())
    tls_file_format = TLS.TLSFileFormatTLSOUT()
    for tlsin in tlsin_list:
        tls_desc_list = tls_file_format.load(open(tlsin, "r"))
        for tls_desc in tls_desc_list:
            tls_file.tls_desc_list.append(tls_desc)
            tls_group = tls_desc.construct_tls_group_with_atoms(struct)
            tls_group.tls_desc = tls_desc
            tls_group_list.append(tls_group)

    ## set the extra Uiso for each atom
    for tls_group in tls_group_list:

        ## minimal/maximal amount of Uiso which has to be added
        ## to the group's atoms to to make Uiso == Uiso_tls
        min_Uiso = 0.0
        max_Uiso = 0.0

        for atm, Utls in tls_group.iter_atm_Utls():
            tls_tf = numpy.trace(Utls) / 3.0
            ref_tf = numpy.trace(atm.get_U()) / 3.0
                
            if ref_tf > tls_tf:
                max_Uiso = max(ref_tf - tls_tf, max_Uiso)
            else:
                min_Uiso = max(tls_tf - ref_tf, min_Uiso)

        ## reduce the TLS group T tensor by min_Uiso so that
        ## a PDB file can be written out where all atoms
        ## Uiso == Uiso_tls

        ## we must rotate the T tensor to its primary axes before
        ## subtracting min_Uiso magnitude from it
        (T_eval, TR) = numpy.linalg.eig(tls_group.T)
        T = numpy.dot(TR, numpy.dot(tls_group.T, numpy.transpose(TR)))

	# FIXME: allclose(some_array, some_scalar)
	# Christoph: The next three lines appear to be the problem (2007-10-04)
        #assert numpy.allclose(T[0,1], 0.0)
        #assert numpy.allclose(T[0,2], 0.0)
        #assert numpy.allclose(T[1,2], 0.0)

        T[0,0] = T[0,0] - min_Uiso
        T[1,1] = T[1,1] - min_Uiso
        T[2,2] = T[2,2] - min_Uiso

        ## now take some of the smallest principal component of T and
        ## move it into the individual atomic temperature factors
        min_T    = min(T[0,0], min(T[1,1], T[2,2]))
        sub_T    = min_T * 0.50
        add_Uiso = min_T - sub_T
        
        T[0,0] = T[0,0] - sub_T
        T[1,1] = T[1,1] - sub_T
        T[2,2] = T[2,2] - sub_T
        
        ## rotate T back to original orientation
        tls_group.T = numpy.dot(
            numpy.transpose(TR),
            numpy.dot(T, TR))

        ## reset the TLS tensor values in the TLSDesc object so they can be saved
        tls_group.tls_desc.set_tls_group(tls_group)

        ## set atm.temp_factor
        for atm, Utls in tls_group.iter_atm_Utls():
            tls_tf = numpy.trace(Utls) / 3.0
            ref_tf = numpy.trace(atm.get_U()) / 3.0

            if ref_tf > tls_tf:
                atm.temp_factor = ((add_Uiso) + ref_tf - tls_tf)*Constants.U2B
                atm.U = None
            else:
                atm.temp_factor = (add_Uiso) * Constants.U2B
                atm.U = None

    FileIO.SaveStructure(fil = xyzout, struct = struct)
    tls_file.save(open(tlsout, "w"))

def phenix_prep(xyzin, phenix_tlsin_list, phenix_tlsout):
    """PHENIX input file. Tells 'phenix.refine' what the TLS groups are.
       Use TLS model + Uiso for each atom.  Output xyzout with the
       residual Uiso only.
    """
    ## load structure
    struct = FileIO.LoadStructure(fil = xyzin)

    ## load and construct TLS groups
    tls_group_list = []
    tls_file = TLS.TLSFile()
    tls_file.set_file_format(TLS.TLSFileFormatPHENIX())
    tls_file_format = TLS.TLSFileFormatPHENIX()
    for tlsin in phenix_tlsin_list:
        tls_desc_list = tls_file_format.load(open(tlsin, "r"))
        for tls_desc in tls_desc_list:
            tls_file.tls_desc_list.append(tls_desc)
            tls_group = tls_desc.construct_tls_group_with_atoms(struct)
            tls_group.tls_desc = tls_desc
            tls_group_list.append(tls_group)

    ## set the extra Uiso for each atom
    for tls_group in tls_group_list:

        ## minimal/maximal amount of Uiso which has to be added
        ## to the group's atoms to to make Uiso == Uiso_tls
        min_Uiso = 0.0
        max_Uiso = 0.0

        for atm, Utls in tls_group.iter_atm_Utls():
            tls_tf = numpy.trace(Utls) / 3.0
            ref_tf = numpy.trace(atm.get_U()) / 3.0
                
            if ref_tf > tls_tf:
                max_Uiso = max(ref_tf - tls_tf, max_Uiso)
            else:
                min_Uiso = max(tls_tf - ref_tf, min_Uiso)

        ## reduce the TLS group T tensor by min_Uiso so that
        ## a PDB file can be written out where all atoms
        ## Uiso == Uiso_tls

        ## we must rotate the T tensor to its primary axes before
        ## subtracting min_Uiso magnitude from it
        (T_eval, TR) = numpy.linalg.eig(tls_group.T)
        T = numpy.dot(TR, numpy.dot(tls_group.T, numpy.transpose(TR)))

	# FIXME: allclose(some_array, some_scalar)
	# Christoph: The next three lines appear to be the problem (2007-10-04)
        #assert numpy.allclose(T[0,1], 0.0)
        #assert numpy.allclose(T[0,2], 0.0)
        #assert numpy.allclose(T[1,2], 0.0)

        T[0,0] = T[0,0] - min_Uiso
        T[1,1] = T[1,1] - min_Uiso
        T[2,2] = T[2,2] - min_Uiso

        ## now take some of the smallest principal component of T and
        ## move it into the individual atomic temperature factors
        min_T    = min(T[0,0], min(T[1,1], T[2,2]))
        sub_T    = min_T * 0.50
        add_Uiso = min_T - sub_T
        
        T[0,0] = T[0,0] - sub_T
        T[1,1] = T[1,1] - sub_T
        T[2,2] = T[2,2] - sub_T
        
        ## rotate T back to original orientation
        tls_group.T = numpy.dot(
            numpy.transpose(TR),
            numpy.dot(T, TR))

        ## reset the TLS tensor values in the TLSDesc object so they can be saved
        tls_group.tls_desc.set_tls_group(tls_group)

        ## set atm.temp_factor
        for atm, Utls in tls_group.iter_atm_Utls():
            tls_tf = numpy.trace(Utls) / 3.0
            ref_tf = numpy.trace(atm.get_U()) / 3.0

            if ref_tf > tls_tf:
                atm.temp_factor = ((add_Uiso) + ref_tf - tls_tf)*Constants.U2B
                atm.U = None
            else:
                atm.temp_factor = (add_Uiso) * Constants.U2B
                atm.U = None

    tls_file.save(open(phenix_tlsout, "w"))
