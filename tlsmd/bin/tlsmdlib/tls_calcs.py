## TLS Minimized Domains (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import popen2
import string

## mmLib
from mmLib.Colors         import *
## tlsmdlib
from misc                 import *


def calc_rmsd_tls_biso(tls_group):
    """Calculate the RMSD of the tls_group using the isotropic TLS model.
    """
    T = tls_group.itls_T
    L = tls_group.itls_L
    S = tls_group.itls_S
    O = tls_group.origin

    msd_sum = 0.0
    
    for atm, uiso_tls in iter_itls_uiso(iter(tls_group), T, L, S, O):
        msd_sum += (U2B*uiso_tls - atm.temp_factor)**2
        
    if len(tls_group)>0:
        msd = msd_sum / len(tls_group)
        rmsd = math.sqrt(msd)
    else:
        rmsd = 0.0

    return rmsd
     

def calc_mean_biso_obs(chainopt):
    """Calculates the mean B value per residue in the chain (as observed in the input structure).
    """
    chain = chainopt["chain"]
    num_res = chain.count_fragments()
    biso = zeros(num_res, Float)

    for frag in chain.iter_fragments():
        n = 0
        b_sum_obs = 0.0

        for atm in frag.iter_all_atoms():
            if atm.include==False: continue

            n += 1
            b_sum_obs += atm.temp_factor

        if n>0:
            biso[frag.ichain] = b_sum_obs / n

    return biso


def calc_mean_biso_tls(chainopt, tlsopt):
    """Calculated the mean B value per residue in the chain as calculated in the chain optimization.
    """
    chain = chainopt["chain"]
    num_res = chain.count_fragments()
    biso = zeros(num_res, Float)

    for i in range(len(tlsopt.tls_list)):
        tls       = tlsopt.tls_list[i]
        tls_group = tls["tls_group"]
        segment   = tls["segment"]

        T = tls_group.itls_T
        L = tls_group.itls_L
        S = tls_group.itls_S
        O = tls_group.origin

        for frag in segment.iter_fragments():
            n = 0
            b_sum_tls = 0.0

            for atm in frag.iter_all_atoms():
                if atm.include==False: continue

                n += 1
                b_sum_tls += U2B * calc_itls_uiso(T, L, S, atm.position - O)

            if n>0:
                biso[frag.ichain] =  b_sum_tls / n

    return biso


def calc_cross_prediction_matrix_rmsd(chainopt, tlsopt):
    chain = chainopt["chain"]

    num_tls = len(tlsopt.tls_list)
    num_res = chain.count_fragments()

    cmtx = zeros((num_tls, num_res), Float)

    for i in range(len(tlsopt.tls_list)):
        tls = tlsopt.tls_list[i]
        tls_group = tls["tls_group"]

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
                if atm.include==False: continue

                n += 1
                b_iso_tls = U2B * calc_itls_uiso(T, L, S, atm.position - O)
                delta = atm.temp_factor - b_iso_tls
                delta2 += delta**2

            if n>0:
                msd = delta2 / n
                rmsd = math.sqrt(msd)

                ## set the cross prediction matrix
                cmtx[i,j] = rmsd

    return cmtx
     

