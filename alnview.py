#!/usr/bin/env python
from pymol.cgo import *
import draw_links
from pymol import cmd
import psab
import sys,re
from scipy import array, zeros
from numpy import *
from pdbParse import loadPDB

PSAdir='/se2.02'
pdbpath='/Theseus'
L=232
alnpath='/msa'


def show(name1, name2, alnpath=alnpath, PSAdir=PSAdir, pdbpath=pdbpath, L=232, seq1color='ruby', seq2color='smudge', benchmarkcolor='yellow', querycolor='grey40'):
    alnstats = psab.alnCM(name1, name2, PSAdir=PSAdir, L=L, alnpath_=alnpath, pdbpath=pdbpath)
    pdbid1 = alnstats['pdbpair'][0]
    chain1 = pdbid1[-1]
    pdbid2 = alnstats['pdbpair'][1]
    chain2 = pdbid2[-1]

    alnmap = load(os.path.join(alnpath, 'alnMap.npz'))

    colors = {}
    colors['gold'] = [1.0, 0.81, 0.14]
    colors['yellow'] = [1.0, 1.0, 0.0]
    colors['red']   = [1.0,0.0,0.0]
    colors['ruby'] = [0.6,0.2,0.2]
    colors['green'] = [0.0,1.0,0.0]
    colors['smudge'] = [0.55, 0.7, 0.4]
    colors['black'] = [0.0,0.0,0.0]
    colors['grey50'] = [0.5,0.5,0.5]
    colors['grey40'] = [0.4,0.4,0.4]
    colors['blue'] = [0.0, 0.0, 1.0]

    objname1 = name1+'-'+pdbid1
    objname2 = name2+'-'+pdbid2

    cmd.load(os.path.join(pdbpath, 'theseus_'+pdbid1+'.pdb'), objname1)
    cmd.load(os.path.join(pdbpath, 'theseus_'+pdbid2+'.pdb'), objname2)
    cmd.hide('all')
    cmd.show('cartoon')
    cmd.color(seq1color, objname1)
    cmd.color(seq2color, objname2)

    pdb1 = loadPDB(os.path.join(pdbpath, 'theseus_'+pdbid1+'.pdb'))
    pdb1 = pdb1[pdb1.chain == chain1.encode('ascii')]

    pdb2 = loadPDB(os.path.join(pdbpath, 'theseus_'+pdbid2+'.pdb'))
    pdb2 = pdb2[pdb2.chain == chain2.encode('ascii')]

    CA1 = (pdb1.atom == b' CA ')
    CA2 = (pdb2.atom == b' CA ')

    alnmap1 = alnmap[pdbid1]
    alnmap2 = alnmap[pdbid2]

    goldmatches1 = alnstats['goldinds'][:,0][alnstats['goldmatches'][:,0]]
    goldmatches2 = alnstats['goldinds'][:,1][alnstats['goldmatches'][:,1]]


    pairs = ((i,j) for i in range(L-1) for j in range(i+1,L))
    p = 1.0
    assert(len(alnstats['posCM'])==L)
    disagreements = alnstats['posCM']=='F'
    F1 = alnmap1[disagreements]
    F2 = alnmap2[disagreements]
    for n,(a,b) in enumerate(zip(F1,F2)):
        print a,b
        if (a == '-') and (b != '-'):
            cmd.color('black', objname2+' and resi '+b)
            continue

        if (b == '-') and (a != '-'):
            cmd.color('black', objname1+' and resi '+a)
            continue

        if (a in goldmatches1):
            cmd.color('red', objname1+' and resi '+a)
        
        if (b in goldmatches2):
            cmd.color('green', objname2+' and resi '+b)
     
        if (a in goldmatches1) or (b in goldmatches2):  
            if a in goldmatches1:
                b_gold = goldmatches2[goldmatches1==a][0]
                draw_links.draw_links(objname1+' and resi {} and n. CA'.format(a), objname2+' and resi {} and n. CA'.format(b_gold), color=colors[benchmarkcolor], object_name='benchmark', radius=0.1)
            if b in goldmatches2:
                a_gold = goldmatches1[goldmatches2==b][0]
                draw_links.draw_links(objname2+' and resi {} and n. CA'.format(b), objname1+' and resi {} and n. CA'.format(a_gold), color=colors[benchmarkcolor], object_name='benchmark', radius=0.1)
  
        if (a != '-') and (b != '-'):
             draw_links.draw_links(objname1+' and resi {} and n. CA'.format(a), objname2+' and resi {} and n. CA'.format(b), color=colors[querycolor], object_name='query', radius=0.1)
            
    cmd.zoom()
    cmd.orient()
