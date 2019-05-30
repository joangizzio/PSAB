#!/usr/bin/env python2
from __future__ import with_statement
from scipy import io, empty, array, rec, dtype
import sys, os, cPickle

"""
ATOM      8  C  BASP A   2      76.283  57.361  60.309  0.50 84.80           C  
ATOM  44035  H2  TIP3 9165     -42.898 -14.686  13.233  1.00  0.00      WT8
ATOM   3949  OXT ARG E 118       1.647   7.647  53.839  1.00 68.00           O
AAAAAAIIIII AAAA RRRRCNNNN    XXXXXXXXYYYYYYYYZZZZZZZZOOOOOOBBBBBB      SSSSEECC
                    ?     ?                                             ????
12345678901234567890123456789012345678901234567890123456789012345678901234567890
0        1         2         3         4         5         6         7         8

#I don't quite follow this table, see ? above.
#note: This is in base-1 indexing, but my arrays are base-0

COLUMNS      DATA TYPE        FIELD      DEFINITION
------------------------------------------------------
 1 -  6      Record name      "ATOM    "
 7 - 11      Integer          serial     Atom serial number.
13 - 16      Atom             name       Atom name.
17           Character        altLoc     Alternate location indicator.
18 - 20      Residue name     resName    Residue name.
22           Character        chainID    Chain identifier.
23 - 26      Integer          resSeq     Residue sequence number.
27           AChar            iCode      Code for insertion of residues.
31 - 38      Real(8.3)        x          Orthogonal coordinates for X in Angs.
39 - 46      Real(8.3)        y          Orthogonal coordinates for Y in Angs.
47 - 54      Real(8.3)        z          Orthogonal coordinates for Z in Angs.
55 - 60      Real(6.2)        occupancy  Occupancy.
61 - 66      Real(6.2)        tempFactor Temperature factor.
77 - 78      LString(2)       element    Element symbol, right-justified.
79 - 80      LString(2)       charge     Charge on the atom.

the atom name field (4 chars) is further subdivided: 1st two chars are element
(right justified, so " C  " is carbon, "CA  " is calcium. Next char is 'distance indicator',
named by greek letters (ABGDEZH), so " CA " is an alpha carbon. Last char is branch number,
so eg " HA3" is the third hydrogen attatched to an atom. Hydrogens are named specially:
The 'distance indicator' is labelled by which carbon the H is attached to, and if there are
many H attached to an atom, they are given names "1H", "2H" etc as the elemtn name, so eg
"1HA3" would be H #1 attached to alpha carbon, 1st branch (If that is even possible).
"""

residueCode = { "GLY": "G", "PRO": "P", "ALA": "A", "VAL": "V", "LEU": "L", 
                "ILE": "I", "MET": "M", "CYS": "C", "PHE": "F", "TYR": "Y", 
                "TRP": "W", "HIS": "H", "LYS": "K", "ARG": "R", "GLN": "Q", 
                "ASN": "N", "GLU": "E", "ASP": "D", "SER": "S", "THR": "T"} 

#to speed up loading time, if a directory 'pickedPDBs' exists the data will
#be stored in pickled form in that directory, so it can be loaded faster the next time.
def loadPDB(filename):    
    #in principle more careful error checking is needed but good enough for now
    if(os.path.exists('./pickledPDBs')):
        pickledName = os.path.join('./pickledPDBs', filename + '.pkl')
        
        if( os.path.exists(pickledName) ):
            with open(pickledName, 'rb') as f:
                pdb = cPickle.load(f)
            return pdb
        else:
            print "creating pickle file"
            pdb = loadTextPDB(filename)
            with open(pickledName, 'wb') as f:
                cPickle.dump(pdb, f, 2)
            return pdb        
    else:
        return loadTextPDB(filename)
        
                    
#loads a pdb from text. Returns a recarray, with entries as in pdbtype
def loadTextPDB(filename):
    with open(filename) as f:
        data = f.readlines()
    
    pdbtype = [ ('atomId',    int), 
                ('atom',      ('S4', [('element', 'S2'), ('distance', 'c'), ('branch', 'c')])), 
                ('altLoc',    'S1'),
                ('resName',   'S4'),
                ('chain',     'S1'),
                ('resID',     'S5'),
                ('resNum',    int),
                ('iCode',     'S1'),
                ('coords',    float,3),
                ('occupancy', float),
                ('beta',      float),
                ('segment',   'S4'),
                ('element',   'S2'),
                ('charge',    'S2'),
                ('weight',    float) ]

    atomLines = [line for line in data if line[0:6] == 'ATOM  ']
    pdb = empty(len(atomLines), dtype(pdbtype))

    for n,line in enumerate(atomLines):
        pdb[n] = ( int(line[6:11]),
                   line[12:16],
                   line[16],
                   line[17:21].strip(),
                   line[21],
                   line[22:27],
                   int(line[22:26]),
                   line[26],
                    
                   ( float(line[30:38]),
                     float(line[38:46]),
                     float(line[46:54]) ),
                    
                   float(line[54:60]),
                   float(line[60:67]),
                   line[72:76].strip(),
                   line[76:78].strip(),
                   line[78:80].strip(),
                   0)

    pdb = rec.array(pdb)

    weights = [(' C', 12.01), (' H', 1.01), (' O', 16.00), (' P', 30.97), 
               (' N', 14.01), ('Na', 22.99), ('Cl', 35.45), (' S', 28.09)]
    for at,w in weights:
        pdb.weight[pdb.element == at] = w
    
    return pdb

#rough preliminary version. Ignores TER etc.
def writePDB(atoms, filename, renumber=False):
    f = open(filename, 'wt')
    writeAtoms(atoms, f)
    f.close()

def writeAtoms(atoms, f):
    for a in atoms:
        f.write("ATOM  %5d %4s%1s%3s%1s%5s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n" % (
                a.atomId, a.atom, a.altLoc, a.resName, a.chain, a.resID, 
                a.coords[0], a.coords[1], a.coords[2], a.occupancy, a.beta, a.element, a.charge))
    
if __name__ == '__main__':
    loadPDB(sys.argv[1])



