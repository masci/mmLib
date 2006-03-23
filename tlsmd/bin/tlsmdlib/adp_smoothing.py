## TLS Motion Determination (TLSMD)
## Copyright 2002-2006 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import numpy

from mmLib import Constants, TLS

import atom_selection
import tls_calcs
import tlsmdmodule


def IsotropicADPDataSmoother(chain, num_smooth = 1):
    """Experimental data smoothing of temperature factors
    """
    print "SMOOTHING CHAIN %s ADPs" % (chain.chain_id)
    print "SMOOTH WINDOW......................: %d" % (2 * num_smooth + 1)
    
    tls_analyzer = tlsmdmodule.TLSModelAnalyzer()
    xlist = atom_selection.chain_to_xmlrpc_list(chain)
    tls_analyzer.set_xmlrpc_chain(xlist)

    num_frags = len(chain)

    smooth_uiso = dict()
    ifrag_start = num_smooth
    ifrag_end = num_frags - num_smooth - 1

    for ifrag in xrange(ifrag_start, ifrag_end + 1):
        smooth_frag = chain[ifrag]
        frag1 = chain[ifrag - num_smooth]
        frag2 = chain[ifrag + num_smooth]

        tlsdict = tls_analyzer.isotropic_fit_segment(frag1.fragment_id, frag2.fragment_id)
        IT, IL, IS, IOrigin = tls_calcs.isotlsdict2tensors(tlsdict)

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

