## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

from   mmTypes              import *
from   AtomMath             import *



class AminoAcidResidueAlgorithms(object):
    def calcMainchainBondLength(self):
        """Calculates the main chain bond lengths: (N-CA, CA-C, C-O, CA-CB,
        CA-(next)N).  The result is returned as a 5-tuple in that order.  Bond
        lengths involving missing atoms are returned as None in the tuple."""
        aN  = self.getAtom('N')
        aCA = self.getAtom('CA')
        aC  = self.getAtom('C')
        aO  = self.getAtom('O')
        aCB = self.getAtom('CB')

        try:
            naN = self.getResidue(1).getAtom('N')
        except AttributeError:
            naN = None
     
        N_CA  = calculateDistance(aN, aCA)
        CA_C  = calculateDistance(aCA, aC)
        C_O   = calculateDistance(aC, aO)
        C_nN  = calculateDistance(aC, naN)
        CA_CB = calculateDistance(aCA, aCB)
        return (N_CA, CA_C, C_O, CA_CB, C_nN)

    def calcMainchainBondAngle(self, conf_id = None):
        """Calculates main chain bond angles (N-CA-C, N-CA-CB, CB-CA-C,
        CA-C-O, CA-C-(next)N, C-(next residue)N-(next residue)CA) and
        returnst the result as a 6-tuple in that order.  Angles involving
        missing atoms are returned as None in the tuple."""
        aN       = self.getAtom('N')
        aCA      = self.getAtom('CA')
        aC       = self.getAtom('C')
        aO       = self.getAtom('O')
        aCB      = self.getAtom('CB')

        naN      = None
        naCA     = None
        next_res = self.getResidue(1)
        if next_res:
            naN  = next_res.getAtom('N')
            naCA = next_res.getAtom('CA')

        N_CA_C   = calculateAngle(aN, aCA, aC)
        CA_C_O   = calculateAngle(aCA, aC, aO)
        N_CA_CB  = calculateAngle(aN, aCA, aCB)
        CB_CA_C  = calculateAngle(aCB, aCA, aC)
        CA_C_nN  = calculateAngle(aCA, aC, naN)
        C_nN_nCA = calculateAngle(aC, naN, naCA)

        return (N_CA_C, N_CA_CB, CB_CA_C, CA_C_O, CA_C_nN, C_nN_nCA) 

    def calcTorsionPsi(self):
        """Calculates the Psi torsion angle of the amino acid.  Raises a
        CTerminal exception if called on a C-terminal residue which does
        not have a Psi torsion angle."""
        next_res = self.getResidue(1)
        if not next_res:
            return None

        aN  = self.getAtom('N')
        aCA = self.getAtom('CA')
        aC  = self.getAtom('C')
        naN = next_res.getAtom('N')
        return calculateTorsionAngle(aN, aCA, aC, naN)

    def calcTorsionPhi(self):
        """Calculates the Phi torsion angle of the amino acid.  Raises a
        NTerminal exception if called on a N-terminal residue which does
        not have a Phi torsion angle."""
        prev_res = self.getResidue(-1)
        if not prev_res:
            return None

        paC = prev_res.getAtom('C')
        aN  = self.getAtom('N')
        aCA = self.getAtom('CA')
        aC  = self.getAtom('C')
        return calculateTorsionAngle(paC, aN, aCA, aC)

    def calcTorsionOmega(self):
        """Calculates the Omega torsion angle of the amino acid. Raises a
        CTerminal exception if called on a C-terminal residue which does
        not have a Omega torsion angle."""
        next_res = self.getResidue(1)
        if not next_res:
            return None

        aCA  = self.getAtom('CA')
        aC   = self.getAtom('C')
        naN  = next_res.getAtom('N')
        naCA = next_res.getAtom('CA')
        return calculateTorsionAngle(aCA, aC, naN, naCA)

    def isCis(self):
        """Returns true if this is a CIS amino acid, otherwise returns false.
        It uses calcTorsionOmega."""
        omega = self.calcTorsionOmega()
        return abs(omega) > (math.pi / 2.0)

    def calcPuckerTorsion(self, conf_id = None):
        """Calculates the Pucker torsion of a ring system.  Returns None
        for Amino Acids which do not have Pucker torsion angles."""
        mon = self.structure().library[self.res_name]
        if not mon.pucker_definition:
            return None

        a1 = self.getAtom(mon.pucker_definition[0])
        a2 = self.getAtom(mon.pucker_definition[1])
        a3 = self.getAtom(mon.pucker_definition[2])
        a4 = self.getAtom(mon.pucker_definition[3])
        return calculateTorsionAngle(a1, a2, a3, a4)

    def calcTorsionChi1(self):
        mon = self.structure().library[self.res_name]
        if not mon.chi1_definition:
            return None
        
        a1 = self.getAtom(mon.chi1_definition[0])
        a2 = self.getAtom(mon.chi1_definition[1])
        a3 = self.getAtom(mon.chi1_definition[2])
        a4 = self.getAtom(mon.chi1_definition[3])
        return calculateTorsionAngle(a1, a2, a3, a4)

    def calcTorsionChi2(self):
        mon = self.structure().library[self.res_name]
        if not mon.chi2_definition:
            return None
        
        a1 = self.getAtom(mon.chi2_definition[0])
        a2 = self.getAtom(mon.chi2_definition[1])
        a3 = self.getAtom(mon.chi2_definition[2])
        a4 = self.getAtom(mon.chi2_definition[3])
        return calculateTorsionAngle(a1, a2, a3, a4)

    def calcTorsionChi3(self):
        mon = self.structure().library[self.res_name]
        if not mon.chi3_definition:
            return None
        
        a1 = self.getAtom(mon.chi3_definition[0])
        a2 = self.getAtom(mon.chi3_definition[1])
        a3 = self.getAtom(mon.chi3_definition[2])
        a4 = self.getAtom(mon.chi3_definition[3])
        return calculateTorsionAngle(a1, a2, a3, a4)

    def calcTorsionChi4(self):
        mon = self.structure().library[self.res_name]
        if not mon.chi4_definition:
            return None
        
        a1 = self.getAtom(mon.chi4_definition[0])
        a2 = self.getAtom(mon.chi4_definition[1])
        a3 = self.getAtom(mon.chi4_definition[2])
        a4 = self.getAtom(mon.chi4_definition[3])
        return calculateTorsionAngle(a1, a2, a3, a4)

    def calcTorsionChi(self):
        """Calculates CHI side-chain torsion angles according to the
        amino acid specific definitions in the AminoAcids library.
        Returns the 4-tuple (CHI1, CHI2, CHI3, CHI4).  Angles involving
        missing atoms, or angles which do not exist for the amino acid
        are returned as None in the tuple."""
        chi1 = self.calcTorsionChi1()
        chi2 = self.calcTorsionChi2()
        chi3 = self.calcTorsionChi3()
        chi4 = self.calcTorsionChi4()
        return (chi1, chi2, chi3, chi4)



class AtomAlgorithms(object):
    def calcAnisotropy(self):
        """Calculates the ansitropy of that atom."""
        ## no Anisotropic values, we have a spherical atom
        if not self.U: return 1.0

        ## build Numeric Python (NumPy) matrix
        m = array([[ self.U[0], self.U[3], self.U[4] ],
                   [ self.U[3], self.U[1], self.U[5] ],
                   [ self.U[4], self.U[5], self.U[2] ]])

        evals = eigenvalues(m)
        ansotropy = min(evals) / max(evals)
        return ansotropy
        
    def iterAtomsByDistance(self):
        """Iterates all atoms in the Structure object from the closest to the
        farthest.  Yields the 2-tuple (dist, atm)."""
        list = []
        for atm in self.structure().iterAtoms():
            list.append((calculateDistance(self, atm), atm))
        list.sort()

        return iter(list)



class BondAlgorithms(object):
    def calcLength(self):
        """Calculate the length of the bond."""
        atm1 = self.atom1()
        atm2 = self.atom2()
        return calculateDistance(atm1, atm2)
    
