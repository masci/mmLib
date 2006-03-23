## TLS Motion Determination (TLSMD)
## Copyright 2002-2006 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import numpy

from mmLib import Constants

import conf
import const

def calc_include_atom(atm, reject_messages = False):
    """Filter out atoms from the model which will cause problems or
    cont contribute to the TLS analysis.
    """
    if atm.position == None:
        return False

    if atm.occupancy < 0.1:
        if reject_messages == True:
            print "calc_include_atom(%s): rejected because of low occupancy" % (atm)
	return False
    
    if numpy.trace(atm.get_U()) <= const.TSMALL:
        if reject_messages == True:
            print "calc_include_atom(%s): rejected because of small Uiso magnitude " % (atm)
        return False

    elif conf.globalconf.include_atoms == "MAINCHAIN":
        if atm.name not in const.MAINCHAIN_ATOMS:
            if reject_messages == True:
                print "calc_include_atom(%s): rejected non-mainchain atom" % (atm)
            return False
    
    return True

def calc_atom_weight(atm):
    """Weight the least-squares fit according to this function.
    """
    assert atm.occupancy >= 0.0 and atm.occupancy <= 1.0
    return atm.occupancy

def chain_to_xmlrpc_list(chain):
    """Converts the Atoms of a Chain/Segment to a list of dictionaries 
    for transfer over xmlrpc.  Only the information required to fit
    TLS groups to the Atoms is set in the Atom dictionaries to
    reduce traffic over the xmlrpc calls.
    """
    xmlrpc_chain = []

    for atm in chain.iter_all_atoms():
        if atm.include is False:
            continue

        atm_desc = {}
        xmlrpc_chain.append(atm_desc)

        atm_desc["name"] = atm.name
        atm_desc["frag_id"] = atm.fragment_id

        frag = atm.get_fragment()
        atm_desc["ifrag"] = chain.index(frag)

        atm_desc["x"] = atm.position[0]
        atm_desc["y"] = atm.position[1]
        atm_desc["z"] = atm.position[2]

        atm_desc["u_iso"] = Constants.B2U * atm.temp_factor

        U = atm.get_U()
        atm_desc["u11"] = U[0,0]
        atm_desc["u22"] = U[1,1]
        atm_desc["u33"] = U[2,2]
        atm_desc["u12"] = U[0,1]
        atm_desc["u13"] = U[0,2]
        atm_desc["u23"] = U[1,2]

        ## calculate weight
        atm_desc["weight"] = calc_atom_weight(atm)

    return xmlrpc_chain
