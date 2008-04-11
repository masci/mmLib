## TLS Motion Determination (TLSMD)
## Copyright 2002-2008 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import math
import numpy

from mmLib import Constants, TLS

import console
import atom_selection
import tls_calcs
import tlsmdmodule


def IsotropicFitSegmentOutlierRejection(chain, frag_id1, frag_id2):
    segment = chain[frag_id1:frag_id2]
    atoms = list(segment.iter_all_atoms())
    
    orig_num_atoms = len(atoms)
    rejected = 0
    
    while True:
        tls_analyzer = tlsmdmodule.TLSModelAnalyzer()
        xlist = atom_selection.chain_to_xmlrpc_list(iter(atoms))
        tls_analyzer.set_xmlrpc_chain(xlist)

        tlsdict = tls_analyzer.isotropic_fit_segment(frag_id1, frag_id2)
        IT, IL, IS, IOrigin = tls_calcs.isotlsdict2tensors(tlsdict)

        num_atoms = 0
        msd_sum = 0.0
        atm_deltab = []

        for atm, uiso in TLS.iter_itls_uiso(iter(atoms), IT, IL, IS, IOrigin):
            num_atoms += 1;
            deltab = atm.temp_factor - (Constants.U2B * uiso)
            msd_sum += deltab**2
            atm_deltab.append((deltab, atm))

        sigma = math.sqrt((msd_sum / num_atoms))
        sigma2 = 2.0 * sigma

        outliers = 0
        for deltab, atm in atm_deltab:
            if abs(deltab) > sigma2:
                atoms.remove(atm)
                outliers += 1

        rejected += outliers

        if outliers == 0 or (num_atoms - outliers) < 10:
            console.stdoutln("SEGMENT %s-%s %d->%d" % (frag_id1, frag_id2, orig_num_atoms, orig_num_atoms - rejected))
            return IT, IL, IS, IOrigin


def IsotropicADPDataSmoother(chain, num_smooth = 1):
    """Experimental data smoothing of temperature factors
    """
    console.endln()
    console.stdoutln("SMOOTHING CHAIN %s ADPs" % (chain.chain_id))
    conesole.kvformat("SMOOTH WINDOW", 2 * num_smooth + 1)
    
    num_frags = len(chain)

    smooth_uiso = dict()
    ifrag_start = num_smooth
    ifrag_end = num_frags - num_smooth - 1

    for ifrag in xrange(ifrag_start, ifrag_end + 1):
        smooth_frag = chain[ifrag]
        frag1 = chain[ifrag - num_smooth]
        frag2 = chain[ifrag + num_smooth]

        IT, IL, IS, IOrigin = IsotropicFitSegmentOutlierRejection(
            chain, frag1.fragment_id, frag2.fragment_id)

        for atm, uiso in TLS.iter_itls_uiso(smooth_frag.iter_all_atoms(), IT, IL, IS, IOrigin):
            smooth_uiso[atm] = uiso
        
        if ifrag == ifrag_start:
            for i in range(ifrag_start):
                smooth_frag = chain[i]
                for atm, uiso in TLS.iter_itls_uiso(smooth_frag.iter_all_atoms(), IT, IL, IS, IOrigin):
                    smooth_uiso[atm] = uiso
        elif ifrag == ifrag_end:
            for i in range(ifrag_end + 1, num_frags):
                smooth_frag = chain[i]
                for atm, uiso in TLS.iter_itls_uiso(smooth_frag.iter_all_atoms(), IT, IL, IS, IOrigin):
                    smooth_uiso[atm] = uiso

    for atm, uiso in smooth_uiso.iteritems():
        atm.temp_factor = Constants.U2B * uiso
        atm.U = numpy.identity(3, float) * uiso

