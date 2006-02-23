## TLS Minimized Domains (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import math
import numpy

from mmLib import Constants, TLS, FileIO


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
            if atm.include == False: continue

            n += 1
            b_sum_obs += atm.temp_factor

        if n>0:
            biso[frag.ichain] = b_sum_obs / n

    return biso


def calc_mean_biso_tls(chain, cpartition):
    """Calculated the mean B value per residue in the chain as calculated in the chain optimization.
    """
    num_res = chain.count_fragments()
    biso = numpy.zeros(num_res, float)

    for i, tls in cpartition.enumerate_tls_segments():
        tls_group = tls.tls_group
        segment   = tls.segment

        T = tls_group.itls_T
        L = tls_group.itls_L
        S = tls_group.itls_S
        O = tls_group.origin

        for frag in segment.iter_fragments():
            n = 0
            b_sum_tls = 0.0

            for atm in frag.iter_all_atoms():
                if atm.include == False:
                    continue

                n += 1
                b_sum_tls += Constants.U2B * TLS.calc_itls_uiso(T, L, S, atm.position - O)

            if n > 0:
                biso[frag.ichain] =  b_sum_tls / n

    return biso


def calc_cross_prediction_matrix_rmsd(chain, cpartition):
    num_tls = cpartition.num_tls_segments()
    num_res = chain.count_fragments()

    cmtx = numpy.zeros((num_tls, num_res), float)

    for i, tls in cpartition.enumerate_tls_segments():
        tls_group = tls.tls_group

        T = tls_group.itls_T
        L = tls_group.itls_L
        S = tls_group.itls_S
        O = tls_group.origin

        for j in range(len(chain)):
            frag = chain[j]

            ## calculate a atom-normalized rmsd deviation for each residue
            n = 0
            delta2 = 0.0

            for atm in frag.iter_all_atoms():
                if atm.include == False:
                    continue

                n += 1
                b_iso_tls = Constants.U2B * TLS.calc_itls_uiso(T, L, S, atm.position - O)
                delta = atm.temp_factor - b_iso_tls
                delta2 += delta**2

            if n>0:
                msd = delta2 / n
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

        n         = 0
        sum_diff2 = 0.0

        for atm, Utls in tls_group.iter_atm_Utls():
            for aatm in atm.iter_alt_loc():
                tls_tf = numpy.trace(Utls)/3.0
                ref_tf = numpy.trace(aatm.get_U())/3.0

                n += 1
                sum_diff2 += (tls_tf - ref_tf)**2
                
                if ref_tf > tls_tf:
                    max_Uiso = max(ref_tf - tls_tf, max_Uiso)
                else:
                    min_Uiso = max(tls_tf - ref_tf, min_Uiso)

        msd = sum_diff2 / n
        rmsd = math.sqrt(msd)

        ## report the percentage of atoms with Uiso within the RMSD
        ntotal = 0
        nrmsd  = 0
        
        for atm, Utls in tls_group.iter_atm_Utls():
            for aatm in atm.iter_alt_loc():
                tls_tf = numpy.trace(Utls)/3.0
                ref_tf = numpy.trace(aatm.get_U())/3.0

                ntotal += 1
                deviation = math.sqrt((tls_tf - ref_tf)**2)
                
                if deviation <= rmsd:
                    nrmsd += 1

        ## reduce the TLS group T tensor by min_Uiso so that
        ## a PDB file can be written out where all atoms
        ## Uiso == Uiso_tls

        ## we must rotate the T tensor to its primary axes before
        ## subtracting min_Uiso magnitude from it
        (T_eval, TR) = numpy.linalg.eigenvectors(tls_group.T)
        T = numpy.matrixmultiply(TR, numpy.matrixmultiply(tls_group.T, numpy.transpose(TR)))

        assert numpy.allclose(T[0,1], 0.0)
        assert numpy.allclose(T[0,2], 0.0)
        assert numpy.allclose(T[1,2], 0.0)

        T[0,0] = T[0,0] - min_Uiso
        T[1,1] = T[1,1] - min_Uiso
        T[2,2] = T[2,2] - min_Uiso

        ## now take half of the smallest principal component of T and
        ## move it into the individual atomic temperature factors

        min_T    = min(T[0,0], min(T[1,1], T[2,2]))
        sub_T    = min_T * 0.80
        add_Uiso = min_T - sub_T
        
        T[0,0] = T[0,0] - sub_T
        T[1,1] = T[1,1] - sub_T
        T[2,2] = T[2,2] - sub_T
        
        ## rotate T back to original orientation
        tls_group.T = numpy.matrixmultiply(numpy.transpose(TR), numpy.matrixmultiply(T, TR))

        ## reset the TLS tensor values in the TLSDesc object so they can be saved
        tls_group.tls_desc.set_tls_group(tls_group)
        
        ## set atm.temp_factor
        for atm, Utls in tls_group.iter_atm_Utls():
            for aatm in atm.iter_alt_loc():
                tls_tf = numpy.trace(Utls)/3.0
                ref_tf = numpy.trace(aatm.get_U())/3.0
                
                if ref_tf > tls_tf:
                    aatm.temp_factor = ((add_Uiso) + ref_tf - tls_tf)*Constants.U2B
                    aatm.U = None
                else:
                    aatm.temp_factor = (add_Uiso) * Constants.U2B
                    aatm.U = None

    FileIO.SaveStructure(fil = xyzout, struct = struct)
    tls_file.save(open(tlsout, "w"))
    
