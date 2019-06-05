import sys, os
from subprocess import call
from numpy import *

"""Calls Byungkook Lee's Seed Extension program for a list of PDBs and their "names"; generates (N*(N-1)/2) fasta files for N PDBs and renames the files"""

pdbnames = genfromtxt(str(sys.argv[1]), dtype=str)[1:,:] 
SEpath = str(sys.argv[2]) #the se program needs to be in this directory
pdbpath = str(sys.argv[3]) #directory of PDB files which are superimposed
pdbprefix = 'theseus_' #structural superposition programs usually put a prefix on the output PDB files


def SE(pdbpairs, namepairs=[], SEpath=SEpath, pdbpath=os.path.join(pdbpath,pdbprefix)):
    if len(namepairs)==0:
        namepairs = pdbpairs
    for pdb, name in zip(pdbpairs, namepairs):
        call('{}se -fasta {}{}.pdb {}{}.pdb'.format(SEpath, pdbpath, pdb[0], pdbpath, pdb[1]).split())
        call('mv {}{}-{}{}.fasta {}{}-{}.fasta'.format(pdbprefix, pdb[0], pdbprefix, pdb[1], SEpath, name[0].replace('/', '_'), name[1].replace('/','_')).split())


names = pdbnames[:,0]
pdbs = pdbnames[:,1]

pdbpairs = array([(pdbs[i], pdbs[j]) for i in range(len(pdbs)-1) for j in range(i+1, len(pdbs))])

namepairs = array([(names[i], names[j]) for i in range(len(names)-1) for j in range(i+1, len(names))])


SE(pdbpairs, namepairs=namepairs)
