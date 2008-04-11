## TLS Minimized Domains (TLSMD)
## Copyright 2002-2008 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import math
import numpy

from mmLib import Constants, TLS, FileIO

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
                biso[frag.ifrag] =  b_sum_tls / n

    return biso


def calc_residue_mean_rmsd(chain, cpartition):
    num_tls = cpartition.num_tls_segments()
    num_res = chain.count_fragments()

    cmtx = numpy.zeros((num_tls, num_res), float)

    for i, tls in enumerate(cpartition.iter_tls_segments()):
        tls_group = tls.tls_group

        T = tls_group.itls_T
        L = tls_group.itls_L
        S = tls_group.itls_S
        O = tls_group.origin

        for j, frag in enumerate(chain):
            ## calculate a atom-normalized rmsd deviation for each residue
            num_atoms = 0
            msd_sum = 0.0

            for atm in frag.iter_all_atoms():
                if atm.include == False:
                    continue

                num_atoms += 1
                b_iso_tls = Constants.U2B * TLS.calc_itls_uiso(T, L, S, atm.position - O)
                delta = atm.temp_factor - b_iso_tls
                msd_sum += delta**2

            if num_atoms > 0:
                msd = msd_sum / num_atoms
                rmsd = math.sqrt(msd)

                ## set the cross prediction matrix
                cmtx[i,j] = rmsd

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

	# Christoph DEBUG: allclose(some_array, some_scalar)
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
    
    # EAM DEBUG
    #return "success"

## PHENIX input file. Tells 'phenix.refine' what the TLS groups are. Christoph Champ, 2007-11-06
#def phenix_prep(xyzin, tlsin_list, xyzout, tlsout):
def phenix_prep(xyzin, phenix_tlsin_list, phenix_tlsout):
    """Use TLS model + Uiso for each atom.  Output xyzout with the
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

	# Christoph DEBUG: allclose(some_array, some_scalar)
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

    #FileIO.SaveStructure(fil = xyzout, struct = struct)
    tls_file.save(open(phenix_tlsout, "w"))
    
    # EAM DEBUG
    #return "success"
