from __future__ import generators

## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import fpformat

from   Numeric              import array
from   LinearAlgebra        import eigenvalues
from   Scientific.Geometry  import Vector, Tensor
from   Composite            import *
from   AtomMath             import *

from   Elements             import ElementNames, ElementMap
from   AminoAcids           import AminoAcidNames, AminoAcidMap
from   NucleicAcids         import NucleicAcidNames, NucleicAcidMap


## exceptions
MissingAtom  = "MissingAtom"
MissingChain = "MissingChain"



class AtomContainer(Composite):
    """Any type of structure container sould be subclassed from AtomContainer.
    This includes atoms."""
    def __init__(self):
        Composite.__init__(self)


    def __str__(self):
        return "AtomContainer"

    def structureIterator(self):
        """Iterates all Structure classes rooted at self."""
        for x in self:
            if isinstance(x, Structure): yield x

    def conformationIterator(self):
        """Iterates all Conformation classes rooted at self."""
        for x in self:
            if isinstance(x, Conformation): yield x

    def moleculeIterator(self):
        """Iterates all Molecule classes rooted at self."""
        for x in self:
            if isinstance(x, Molecule): yield x

    def chainIterator(self):
        """Iterates all Chain classes rooted at self."""
        for x in self:
            if isinstance(x, Chain): yield x

    def polypeptideChainIterator(self):
        """Iterates all PolyeptideChain classes rooted at self."""
        for x in self:
            if isinstance(x, PolypeptideChain): yield x
                
    def residueIterator(self):
        """Iterates all Residue classes rooted at self."""
        for x in self:
            if isinstance(x, Residue): yield x
                
    def aminoAcidResidueIterator(self):
        """Iterates all AminoAcideResidue classes rooted at self."""
        for x in self:
            if isinstance(x, AminoAcidResidue): yield x

    def nucleicAcidResidueIterator(self):
        """Iterates all NucleicAcideResidue classes rooted at self."""
        for x in self:
            if isinstance(x, NucleicAcidResidue): yield x

    def atomIterator(self, conf_id = None):
        """Iterates all Atom classes rooted at self."""
        if conf_id:
            def block_func(x):
                if not isinstance(x, Conformation): return 0
                x_conf_id = x.getID()
                return x_conf_id != conf_id

            for x in self.blockingIterator(block_func):
                if isinstance(x, Atom): yield x
            
        else:
            for x in self:
                if isinstance(x, Atom): yield x
        
    def bondIterator(self, conf_id = None):
        """Special case of atomIterator which iterates all bonds in
        atoms rooted at self."""
        visited_bonds = []
        for atm in self.atomIterator(conf_id):
            for bond in atm.bond_list:
                if bond in visited_bonds:
                    continue
                visited_bonds.append(bond)
                yield bond

    def getStructure(self):
        """Returns self's closest parent of class System."""
        for parent in self.getParentList():
            if isinstance(parent, Structure):
                return parent
        return None

    def getConformation(self):
        """Returns self's closest parent of class System."""
        for parent in self.getParentList():
            if isinstance(parent, Conformation):
                return parent
        return None
   
    def getMolecule(self):
        """Returns self's closest parent of class Molecule."""
        for parent in self.getParentList():
            if isinstance(parent, Molecule):
                return parent
        return None

    def getChain(self):
        """Returns self's closest parent of class Chain."""
        for parent in self.getParentList():
            if isinstance(parent, Chain):
                return parent
        return None

    def getPolypeptideChain(self):
        """Returns self's closest parent of class PolypeptideChain."""
        for parent in self.getParentList():
            if isinstance(parent, PolypeptideChain):
                return parent
        return None
    
    def getResidue(self):
        """Returns self's closest parent of class Residue."""
        for parent in self.getParentList():
            if isinstance(parent, Residue):
                return parent
        return None
    
    def getAminoAcidResidue(self):
        """Returns self's closest parent of class Residue."""
        for parent in self.getParentList():
            if isinstance(parent, AminoAcidResidue):
                return parent
        return None

    def getNucleicAcidResidue(self):
        """Returns self's closest parent of class NucleicAcidResidue."""
        for parent in self.getParentList():
            if isinstance(parent, NucleicAcidResidue):
                return parent
        return None

    def getConformations(self):
        """Iterates the tree starting at self and returns a sorted unique
        list of all the conformation IDs."""
        conf_list = []
        for conf in self.conformationIterator():
            conf_id = conf.getID()
            if conf_id and conf_id not in conf_list:
                conf_list.append(conf_id)
        conf_list.sort()
        return conf_list



class Conformation(AtomContainer):
    """Container class descirbing multiple conformations of molecules,
    atoms, and chains."""

    def __init__(self):
        AtomContainer.__init__(self)

        self.conformation_id = ""

    def __str__(self):
        return "Conformation::%s" % (self.conformation_id)

    def setID(self, cid):
        self.conformation_id = str(cid)

    def getID(self):
        return self.conformation_id



class Structure(AtomContainer):
    """Container for a Structure of molecules."""
    def __init__(self):
        AtomContainer.__init__(self)

        self.title                = ""
        self.id                   = ""
        self.key_words            = ""
        self.orig_deposition_date = ""


    def __str__(self):
        return "Structure::%s" % (str(self.name))

    def setTitle(self, title):
        """Set the title of the structure."""
        self.title = str(title)

    def getTitle(self):
        """Return the title of the structure."""
        return self.title

    def setID(self, id):
        """Set the ID of the structure, comes from PDB::HEADER.id_code or
        mmCIF::_entity.id."""
        self.id = str(id)

    def getID(self):
        """Get the ID of the structure, comes from PDB::HEADER.id_code or
        mmCIF::_entity.id."""
        return self.id

    def setKeywords(self, key_words):
        self.key_words = str(key_words)

    def getKeywords(self):
        return self.key_words

    def setOriginalDepositionDate(self, date):
        self.orig_deposition_date = str(date)

    def getOriginalDepositionDate(self):
        return self.orig_deposition_date

    def setAuthorList(self, list):
        self.author_list = list

    def getAuthorList(self):
        return self.author_list

    def setDescription(self, desc):
        """Description from PDB::COMPND or mmCIF::entity.pdbx_description."""
        self.description = str(desc)

    def getDescription(self):
        return self.description

    def setExperimentalMethod(self, experimental_method):
        self.experimental_method = str(experimental_method)
    
    def getExperimentalMethod(self):
        return self.experimental_method



class Molecule(AtomContainer):
    """Container for a molecule."""
    def __init__(self):
        AtomContainer.__init__(self)

        self.name          = ""
        self.molecule_type = ""
        self.polymer_type  = ""


    def __str__(self):
        tstr = "Molecule"

        if self.molecule_type != None:
            tstr += "::" + str(self.molecule_type)
        
        if self.name != None:
            tstr += "::" + str(self.name)
        
        return tstr

    def setName(self, name):
        self.name = str(name)

    def getName(self):
        return self.name

    def setType(self, mtype):
        """Sets the type of molecule."""
        self.molecule_type = str(mtype)

    def getType(self):
        """Returns the molecule type."""
        return self.molecule_type

    def setPolymerType(self, ptype):
        """Sets the polymer type of the molecule."""
        self.polymer_type = str(ptype)

    def getPolymerType(self):
        """Returns the polymer type of the molecule."""
        return self.polymer_type

    def isProtein(self):
        """Returns true if the molecule is a protein, otherwise
        returns false."""
        protein_types = ["polypeptide", "polypeptide(L)"]

        if self.molecule_type and self.molecule_type != "polymer":
            return 0

        if self.polymer_type:
            return self.polymer_type in protein_types

        ## simple algorithim to figure out if this is a protein
        num_aa_residues = 0
        for res in self.aminoAcidResidueIterator():
            num_aa_residues += 1

        if num_aa_residues > 0:
            self.molecule_type = "polymer"
            self.polymer_type  = "polypeptide"
            return 1

        return 0
            


class Solvent(AtomContainer):
    """Container for solvent molecules."""
    def __init__(self):
        AtomContainer.__init__(self)

    def __str__(self):
        return "Solvent"



## NOTE: move the constructPeptide bonds here, rename it 
##       constructLinkageBonds or somthing like that

class Chain(AtomContainer):
    """Conatiner for residues of the same type."""
    def __init__(self):
        AtomContainer.__init__(self)

        self.chain_id = ""


    def __str__(self):
        return "Chain::%s" % (self.chain_id)

    def setID(self, id):
        """Set the chain ID.  This should be a 1-char string.  It is
        converted to upper case set as lower case."""
        self.chain_id = str(id).upper()

    def getID(self):
        """Returns the chain ID."""
        return self.chain_id

    def getResidueBySequenceNumber(self, seq_num):
        """Returns the residue of the given sequence ID."""
        for res in self.residueIterator():
            if res.getSequenceNumber() == seq_num:
                return res
        return None



class NucleicAcidChain(Chain):
    def __init__(self):
        Chain.__init__(self)

    def __str__(self):
        return "NucleicAcidChain::%s" % (self.chain_id)



class PolypeptideChain(Chain):
    def __init__(self):
        Chain.__init__(self)

    def __str__(self):
        return "PolypeptideChain::%s" % (self.chain_id)

    def constructPeptideBonds(self):
        """Constructs all inter-residue bonds in the Polypeptide chain."""
        for aa_res in self.aminoAcidResidueIterator():

            next_aa_res = aa_res.getNextResidue()
            if next_aa_res == None:
                continue

            ## make a unique list of the conformations found in both residues
            conf_list = aa_res.getConformations()
            for conf_id in next_aa_res.getConformations():
                if conf_id not in conf_list:
                    conf_list.append(conf_id)

            if not conf_list:
                conf_list.append(None)

            ## iterate through the conformation list constructing bonds
            for conf_id in conf_list:
                try:
                    aC = aa_res.getAtomByLabel("C", conf_id)
                except MissingAtom, label:
                    continue

                try:
                    naN = next_aa_res.getAtomByLabel("N", conf_id)
                except MissingAtom, label:
                    continue
            
                aC.createBond(naN)



class Residue(AtomContainer):
    """Abstract class representing a monomer residue."""
    def __init__(self):
        AtomContainer.__init__(self)

        self.name         = ""
        self.sequence_num = 0
        self.sequence_id  = ""
        self.icode        = ""


    def __str__(self):
        return "Residue::%s::%s" % (str(self.name), str(self.sequence_id))

    def getName(self):
        """Returns the name of the Residue."""
        return self.name

    def setName(self, name):
        """Sets the name of the Residue."""
        self.name = str(name).upper()

    def setSequenceNumber(self, num):
        self.sequence_num = int(num)

    def getSequenceNumber(self):
        return self.sequence_num

    def getSequenceID(self):
        """Returns the sequence ID of the Residue.  This is the (resSeq, iCdoe)
        tuple."""
        return self.sequence_id

    def setSequenceID(self, seq_id):
        """Sets the sequence ID of  the Residue.  This is the resSeq+iCode
        string."""
        self.sequence_id = str(seq_id)

    def getAtomByLabel(self, label, conf_id = None):
        """Searches for and returns the atom in the residue with the
        argument label.  A ResidueError exception is raised if the atom
        is not found."""
        for atom in self.atomIterator(conf_id):
            if atom.getAtomLabel() == label:
                return atom
        raise MissingAtom, label

    def getAtom(self, label, conf_id = None):
        """Searches for and returns the atom in the residue with the
        argument label.  Returns None if the atom is not found."""
        try:
            return self.getAtomByLabel(label, conf_id)
        except MissingAtom:
            return None

    def getPrevResidue(self):
        """Returns the previous Residue in the chain."""
        prev_res = None

        chain = self.getChain()
        if not chain:
            raise MissingChain, str(self)
        
        for res in chain.residueIterator():
            if  res == self: return prev_res
            prev_res = res

        return None

    def getNextResidue(self):
        """Returns the next Residue in the chain."""
        nflag = 0
        
        chain = self.getChain()
        if not chain:
            raise MissingChain, str(self)

        for res in chain.residueIterator():
            if   res == self: nflag = 1
            elif nflag:       return res

        return None



class NucleicAcidResidue(Residue):
    """Container and methods for nucleic acid residues."""
    def __str__(self):
        return "NucleicAcidResidue::%s::%s" % (
            str(self.name), str(self.sequence_id))



class AminoAcidResidue(Residue):
    """Container and methods for amino acid residues."""
    def __init__(self):
        Residue.__init__(self)
        self.aa = None


    def __str__(self):
        return "AminoAcidResidue::%s::%s" % (
            str(self.name), str(self.sequence_id))

    def setName(self, name):
        """Set the name of the amino acid residue.  The name is checked
        aginst the names in the AminoAcids library and a ValueError exception
        is raised if the name isn't defined."""
        Residue.setName(self, name)
        try:
            self.aa = AminoAcidMap[self.name]
        except KeyError:
            print "unknown amino acid type=%s" % (self.name)
            self.aa = AminoAcidMap["UNK"]

    def isNTerminal(self):
        """Returns true if the residue is the N-Terminal, otherwise
        returns false."""
        if self.sequence_num == 1:
            return 1
        return 0
    
    def isCTerminal(self):
        """Returns true if the residue is the C-Terminal, otherwise
        returns false."""
        if not self.getNextResidue():
            return 1
        return 0

    def constructBonds(self):
        """Constructs all bonds defined in the AminoAcid.py dictionary
        for this monomer."""
        conf_list = self.getConformations()
        if not conf_list:
            conf_list.append(None)

        for (label1, label2) in self.aa.bond_list:
            for conf_id in conf_list:
                atm1 = self.getAtom(label1, conf_id)
                atm2 = self.getAtom(label2, conf_id)
                if atm1 and atm2:
                    atm1.createBond(atm2)

    def calcMainchainBondLength(self, conf_id = None):
        """Calculates the main chain bond lengths: (N-CA, CA-C, C-O, CA-CB,
        CA-(next)N).  The result is returned as a 5-tuple in that order.  Bond
        lengths involving missing atoms are returned as None in the tuple."""
        aN = self.getAtom('N', conf_id)
        aCA = self.getAtom('CA', conf_id)
        aC = self.getAtom('C', conf_id)
        aO = self.getAtom('O', conf_id)
        aCB = self.getAtom('CB', conf_id)

        naN = None
        next_res = self.getNextResidue()
        if next_res:
            naN = next_res.getAtom('N', conf_id)
     
        N_CA_distance = calculateDistance(aN, aCA)
        CA_C_distance = calculateDistance(aCA, aC)
        C_O_distance = calculateDistance(aC, aO)
        
        if naN:
            C_nN_distance = calculateDistance(aC, naN)
        else:
            C_nN_distance = None

        if aCB:
            CA_CB_distance = calculateDistance(aCA, aCB)
        else:
            CA_CB_distance = None

        return (N_CA_distance, CA_C_distance, C_O_distance,
                CA_CB_distance, C_nN_distance)

    def calcMainchainBondAngle(self, conf_id = None):
        """Calculates main chain bond angles (N-CA-C, N-CA-CB, CB-CA-C,
        CA-C-O, CA-C-(next)N, C-(next residue)N-(next residue)CA) and
        returnst the result as a 6-tuple in that order.  Angles involving
        missing atoms are returned as None in the tuple."""
        aN = self.getAtomByLabel('N', conf_id)
        aCA = self.getAtomByLabel('CA', conf_id)
        aC = self.getAtomByLabel('C', conf_id)
        aO = self.getAtomByLabel('O', conf_id)
        aCB = self.getAtom('CB', conf_id)

        naN = None
        naCA = None
        next_res = self.getNextResidue()
        if next_res:
            naN = next_res.getAtomByLabel('N', conf_id)
            naCA = next_res.getAtomByLabel('CA', conf_id)

        N_CA_C_angle = calculateAngle(aN, aCA, aC)
        CA_C_O_angle = calculateAngle(aCA, aC, aO)
        
        if aCB:
            N_CA_CB_angle = calculateAngle(aN, aCA, aCB)
            CB_CA_C_angle = calculateAngle(aCB, aCA, aC)
        else:
            N_CA_CB_angle = None
            CB_CA_C_angle = None

        if naN:
            CA_C_nN_angle = calculateAngle(aCA, aC, naN)
        else:
            CA_C_nN_angle = None

        if naCA:
            C_nN_nCA_angle = calculateAngle(aC, naN, naCA)
        else:
            C_nN_nCA_angle = None

        return (N_CA_C_angle, N_CA_CB_angle, CB_CA_C_angle,
                CA_C_O_angle, CA_C_nN_angle, C_nN_nCA_angle) 
        
    def calcTorsionPsi(self, conf_id = None):
        """Calculates the Psi torsion angle of the amino acid.  Raises a
        CTerminal exception if called on a C-terminal residue which does
        not have a Psi torsion angle."""
        next_res = self.getNextResidue()
        if not next_res:
            return None

        aN = self.getAtomByLabel('N', conf_id)
        aCA = self.getAtomByLabel('CA', conf_id)
        aC = self.getAtomByLabel('C', conf_id)
        naN = next_res.getAtomByLabel('N', conf_id)
        return calculateTorsionAngle(aN, aCA, aC, naN)

    def calcTorsionPhi(self, conf_id = None):
        """Calculates the Phi torsion angle of the amino acid.  Raises a
        NTerminal exception if called on a N-terminal residue which does
        not have a Phi torsion angle."""
        prev_res = self.getPrevResidue()
        if not prev_res:
            return None

        paC = prev_res.getAtomByLabel('C', conf_id)
        aN = self.getAtomByLabel('N', conf_id)
        aCA = self.getAtomByLabel('CA', conf_id)
        aC = self.getAtomByLabel('C', conf_id)
        return calculateTorsionAngle(paC, aN, aCA, aC)

    def calcTorsionOmega(self, conf_id = None):
        """Calculates the Omega torsion angle of the amino acid. Raises a
        CTerminal exception if called on a C-terminal residue which does
        not have a Omega torsion angle."""
        next_res = self.getNextResidue()
        if not next_res:
            return None

        aCA = self.getAtomByLabel('CA', conf_id)
        aC = self.getAtomByLabel('C', conf_id)
        naN = next_res.getAtomByLabel('N', conf_id)
        naCA = next_res.getAtomByLabel('CA', conf_id)
        return calculateTorsionAngle(aCA, aC, naN, naCA)

    def isProline(self):
        """Returns true if the amino acid is Proline, otherwise false.
        This is for convience."""
        return self.name == "PRO"
    
    def isCis(self, conf_id = None):
        """Returns true if this is a CIS amino acid, otherwise returns false.
        It uses calcTorsionOmega."""
        omega = self.calcTorsionOmega(conf_id)
        if abs(omega) > math.pi / 2.0:
            return 1
        return 0

    def calcPuckerTorsion(self, conf_id = None):
        """Calculates the Pucker torsion of a ring system.  Returns None
        for Amino Acids which do not have Pucker torsion angles."""
        if not self.aa.pucker_definition:
            return None

        a1 = self.getAtomByLabel(self.aa.pucker_definition[0], conf_id)
        a2 = self.getAtomByLabel(self.aa.pucker_definition[1], conf_id)
        a3 = self.getAtomByLabel(self.aa.pucker_definition[2], conf_id)
        a4 = self.getAtomByLabel(self.aa.pucker_definition[3], conf_id)
        return calculateTorsionAngle(a1, a2, a3, a4)

    def calcTorsionChi1(self, conf_id = None):
        if not self.aa.chi1_definition:
            return None
        
        a1 = self.getAtomByLabel(self.aa.chi1_definition[0], conf_id)
        a2 = self.getAtomByLabel(self.aa.chi1_definition[1], conf_id)
        a3 = self.getAtomByLabel(self.aa.chi1_definition[2], conf_id)
        a4 = self.getAtomByLabel(self.aa.chi1_definition[3], conf_id)
        return calculateTorsionAngle(a1, a2, a3, a4)

    def calcTorsionChi2(self, conf_id = None):
        if not self.aa.chi2_definition:
            return None
        
        a1 = self.getAtomByLabel(self.aa.chi2_definition[0], conf_id)
        a2 = self.getAtomByLabel(self.aa.chi2_definition[1], conf_id)
        a3 = self.getAtomByLabel(self.aa.chi2_definition[2], conf_id)
        a4 = self.getAtomByLabel(self.aa.chi2_definition[3], conf_id)
        return calculateTorsionAngle(a1, a2, a3, a4)

    def calcTorsionChi3(self, conf_id = None):
        if not self.aa.chi1_definition:
            return None
        
        a1 = self.getAtomByLabel(self.aa.chi3_definition[0], conf_id)
        a2 = self.getAtomByLabel(self.aa.chi3_definition[1], conf_id)
        a3 = self.getAtomByLabel(self.aa.chi3_definition[2], conf_id)
        a4 = self.getAtomByLabel(self.aa.chi3_definition[3], conf_id)
        return calculateTorsionAngle(a1, a2, a3, a4)

    def calcTorsionChi4(self, conf_id = None):
        if not self.aa.chi1_definition:
            return None
        
        a1 = self.getAtomByLabel(self.aa.chi4_definition[0], conf_id)
        a2 = self.getAtomByLabel(self.aa.chi4_definition[1], conf_id)
        a3 = self.getAtomByLabel(self.aa.chi4_definition[2], conf_id)
        a4 = self.getAtomByLabel(self.aa.chi4_definition[3], conf_id)
        return calculateTorsionAngle(a1, a2, a3, a4)

    def calcTorsionChi(self, conf_id = None):
        """Calculates CHI side-chain torsion angles according to the
        amino acid specific definitions in the AminoAcids library.
        Returns the 4-tuple (CHI1, CHI2, CHI3, CHI4).  Angles involving
        missing atoms, or angles which do not exist for the amino acid
        are returned as None in the tuple."""
        try:
            chi1 = self.calcTorsionChi1(conf_id)
        except MissingAtom:
            chi1 = None

        try:
            chi2 = self.calcTorsionChi2(conf_id)
        except MissingAtom:
            chi2 = None
            
        try:
            chi3 = self.calcTorsionChi3(conf_id)
        except MissingAtom:
            chi3 = None
            
        try:
            chi4 = self.calcTorsionChi4(conf_id)
        except MissingAtom:
            chi4 = None

        return (chi1, chi2, chi3, chi4)


class Atom(AtomContainer):
    """Individual atom and properties."""

    def __init__(self):
        AtomContainer.__init__(self)

        self.element      = ""
        self.el           = None
        self.charge       = 0.0
        self.position     = Vector(0.0, 0.0, 0.0)
        self.occupancy    = 1.0
        self.temp_factor  = 1.0
        self.U            = None
        self.radius       = 0.0
        self.atom_label   = ""
        self.velocity     = None
        self.force        = None
        self.bond_list    = []


    def __str__(self):
        return "Atom::%s" % (str(self.atom_label))

    def setElement(self, element):
        try:
            self.el = ElementMap[element]
        except KeyError:
            raise SystemExit, "Unknown element %s" % (element)
        self.element = self.el.symbol
        
    def getElement(self):
        return self.element
    
    def setCharge(self, charge):
        self.charge = float(charge)

    def getCharge(self):
        return self.charge

    def setOccupancy(self, occupancy):
        self.occupancy = float(occupancy)

    def getOccupancy(self):
        return self.occupancy

    def setTemperatureFactor(self, tf):
        self.temp_factor = float(tf)

    def getTemperatureFactor(self):
        return self.temp_factor

    def setU(self, U):
        self.U = U

    def getU(self):
        return self.U

    def setPosition(self, position):
        self.position = position

    def getPosition(self):
        return self.position

    def setRadius(self, radius):
        self.radius = float(radius)

    def setAtomLabel(self, atom_label):
        self.atom_label = str(atom_label)

    def getAtomLabel(self):
        return str(self.atom_label).upper()

    def setVelocity(self, velocity):
        self.velocity = velocity

    def setForce(self, force):
        self.force = force

    def getForce(self):
        return self.force

    def countBonds(self):
        return len(self.bond_list)

    def getBondList(self):
        return self.bond_list

    def getBond(self, atom):
        """Returns the bond connection self to the argument atom."""
        for bond in self.bond_list:
            if bond.isBondOf(atom):
                return bond
        return None

    def createBond(self, atom):
        """Create a new bond between self and the argument atom."""
        ## check if the bond already exists
        bond = self.getBond(atom)
        if bond:
            return bond
        
        ## create new bond
        bond = Bond()
        bond.setFirstAtom(self)
        bond.setSecondAtom(atom)
        self.bond_list.append(bond)
        atom.bond_list.append(bond)
        return bond
                              
    def destroyBond(self, atom):
        """Destroy bond to atom."""
        bond = self.getBond(atom)
        if bond:
            self.bond_list.remove(bond)
            atom.bond_list.remove(bond)
            del bond

    def destroyBonds(self):
        """Destory all bonds."""
        for bond in self.bond_list:
            self.bond_list.remove(bond)
            del bond

    def hasBond(self, bond):
        """Returns true if the atom is involved in the bond."""
        return bond in self.bond_list

    def isBound(self):
        """Returns true if the atom is bound to any other atom."""
        return len(self.bond_list) > 0

    def calcAnisotropy(self):
        """Calculates the ansitropy of that atom."""
        ## no Anisotropic values, we have a spherical atom
        if self.U == None:
            return 1.0

        ## build Numeric Python (NumPy) matrix
        m = array([[ self.U[0], self.U[3], self.U[4] ],
                   [ self.U[3], self.U[1], self.U[5] ],
                   [ self.U[4], self.U[5], self.U[2] ]])

        evals = eigenvalues(m)
        ansotropy = min(evals) / max(evals)
        return ansotropy
        


class Bond:
    """Represents a bond between two atoms."""
    def __init__(self):
        self.first_atom   = None
        self.second_atom  = None
        self.bond_type    = "covalent"
        self.distance     = 0.0
        self.distance_esd = 0.0


    def __str__(self):
        a1 = "X"
        a2 = "X"
        conn = "-"

        if self.first_atom: a1 = self.first_atom.getAtomLabel() or "X"
        if self.second_atom: a2 = self.second_atom.getAtomLabel() or "X"

        if self.distance and self.distance_esd:
            d = fpformat.fix(self.distance, 3)
            d_esd = fpformat.fix(self.distance_esd, 3)
            conn = "-(%s, %s)-" % (d, d_esd)
        
        return a1 + conn + a2
    
    def setFirstAtom(self, atom):
        self.first_atom = atom

    def getFirstAtom(self):
        return self.first_atom

    def setSecondAtom(self, atom):
        self.second_atom = atom

    def getSecondAtom(self):
        return self.second_atom

    def getAtoms(self):
        """Return a tuple of the atoms in the bond."""
        return (self.first_atom, self.second_atom)

    def setType(self, btype):
        """Sets the type of the bond.  This must be one of the text
        strings in Bond.TYPES."""
        self.bond_type = btype  

    def getType(self):
        return self.bond_type

    def getPartner(self, atom):
        """Return the partner (bound) atom in the bond."""
        if atom == self.first_atom:
            return self.second_atom
        elif atom == self.second_atom:
            return self.first_atom
        raise StructureError, 'Argument atom not in bond'

    def isBondOf(self, atom):
        """Return false if the atom is not a member of the bond."""
        return atom == self.first_atom or atom == self.second_atom

    def setDistance(self, dist, dist_esd):
        """Set a distance for the bond."""
        self.distance = float(dist)
        self.distance_esd = float(dist_esd)

    def getDistance(self):
        return (self.distance, self.distance_esd)

    def calcLength(self):
        """Calculate the bond length bond length."""
        return calculateDistance(self.first_atom, self.second_atom)
    
