## TLS Motion Determination (TLSMD)
## Copyright 2002-2009 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import numpy
import itertools

from mmLib import Constants, Library

import conf
import const
import console

def iter_mainchain_atoms(atom_iter):
    filter = lambda atm: const.MAINCHAIN_ATOMS.has_key(atm.res_name)
    return itertools.ifilter(filter, atom_iter)

def calc_include_atom(atm, reject_messages = False):
    """Filter out atoms from the model which will cause problems or
    cont contribute to the TLS analysis.
    """
    if atm.position == None:
        return False

    if atm.occupancy < 0.1:
        if reject_messages == True:
            console.stdoutln("calc_include_atom(%s): rejected because of low occupancy" % (atm))
	return False

    if atm.occupancy > 1.0:
        atm.occupancy = 1.0
        console.stdoutln("calc_include_atom(%s): atom occupancy greator than 1.0; truncating" % (atm))
    
    if numpy.trace(atm.get_U()) <= const.TSMALL:
        if reject_messages == True:
            console.stdoutln("calc_include_atom(%s): rejected because of small Uiso magnitude " % (atm))
        return False

    elif conf.globalconf.include_atoms == "MAINCHAIN":
        if const.MAINCHAIN_ATOM_DICT.has_key(atm.name) is False:
            if reject_messages == True:
                console.stdoutln("calc_include_atom(%s): rejected non-mainchain atom" % (atm))
            return False
    
    return True

def calc_atom_weight(atm):
    """Weight the least-squares fit according to this function.
    """
    return atm.occupancy

def chain_to_xmlrpc_list(atom_iter):
    """Converts the Atoms of a Chain/Segment to a list of dictionaries 
    for transfer over xmlrpc.  Only the information required to fit
    TLS groups to the Atoms is set in the Atom dictionaries to
    reduce traffic over the xmlrpc calls.
    """
    xmlrpc_chain = []

    atom_iter = itertools.ifilter(lambda atm: atm.include, atom_iter)
    for atm in atom_iter:
        atm_desc = {}
        xmlrpc_chain.append(atm_desc)

        atm_desc["name"] = atm.name
        atm_desc["frag_id"] = atm.fragment_id

        atm_desc["ifrag"] = 0

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
