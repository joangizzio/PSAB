import sys, os
from subprocess import call
from numpy import *

pdbnames = genfromtxt(str(sys.argv[1]), dtype=str)[1:,:]
SEpath = str(sys.argv[2])
pdbpath = str(sys.argv[3])
pdbprefix = 'theseus_'


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
