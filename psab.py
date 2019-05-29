#Pairwise Sequence Alignment Benchmarking
import os, glob
from numpy import *
import pdbParse

L = 232
pdbpath = '../Theseus/split_peng/'
alnpath = '../../../../Kinase_Dataset/structure/k232_04-01-2019/'
pdbprefix = 'theseus_'
SEpath = '/se2.02'
#MLSpath = '../../../../Bromodomain_Dataset/structure/Bromodomain80/Theseus_useqs/'
#alnpath = '../../../../Bromodomain_Dataset/structure/Bromodomain80_v2'


def findall(cond):
    return argwhere(cond).flatten()

def findfirst(cond):
    res = argwhere(cond).flatten()
    return res[0] if len(res) > 0 else None

def MSAparse(alnpath=alnpath, L=L): 

    with open(os.path.join(alnpath, 'pdbseqs_full_ID')) as f: #pdbseqs_full_ID
        MSA = dict([l.split() for l in f.readlines()])

    for name in MSA.keys():
        MSA[name.upper()] = MSA.pop(name)

    with open(os.path.join(alnpath, 'pdbseqIndices')) as f:
        dat = [l.split() for l in f.readlines()]
        indmaps = dict([(d[0].upper(), array([int(x) for x in d[1:]])) for d in dat])

    for name in indmaps.keys():
        indmaps[name.upper()] = indmaps.pop(name)

    pdb2alnmap = load(os.path.join(alnpath,'alnMap.npz'))

    return {'dots': MSA, 'pdbinds': indmaps, 'alnPDBinds': pdb2alnmap}


def getpairs(names, mode='Lc2'):
    if mode=='Lc2':
        pairs = array([(names[a], names[b]) for a in range(len(names)-1) for b in range(a+1, len(names))])
    if mode=='LxL':
        pairs = array([(names[a], names[b]) for a in range(len(names)) for b in range(len(names))])
    return pairs

def missingres(pdbs, L=None, alnpath=None, msaparse=None):
    validids = []
    Nmissing = []
    if (alnpath != None) and (L != None):
        msaparse = MSAparse(alnpath=alnpath, L=L)
    for name in pdbs:
        try:
            inds = msaparse['pdbinds'][name]
            seqstart = argmax(inds>-1)
            seqend = -1 * argmax(inds[::-1]>-1) - 1
            Nmissing.append(sum(inds[seqstart:seqend]==-1))
            validids.append(name)
        except:
            print('ERROR: '+name)
    validids = array(validids)
    Nmissing = array(Nmissing)

    return validids, Nmissing
    
def genes_missingres(pdbs, genes, L=None, alnpath=None, msaparse=None):
    assert pdbs.shape == genes.shape, 'PDB ids must map directly to gene names'
    geneprofile = {}
    for gene in unique(genes):
        structs = pdbs[genes==gene]
        geneprofile[gene] = missingres(structs, L=L, alnpath=alnpath, msaparse=msaparse) #returns (PDB ids, number of missing residues)
    return geneprofile
 
def goodpdbs(geneprofile, Nmissingthresh=0):
    chosenpdbs = []
    chosengenes = []
    Nmissing = []
    skippedgenes = []
    for gene in geneprofile.keys():
        if any(geneprofile[gene][1] <= Nmissingthresh):
            chosenpdbs.append(geneprofile[gene][0][geneprofile[gene][1]<=Nmissingthresh])
            chosengenes.append(gene)
            Nmissing.append(geneprofile[gene][1][geneprofile[gene][1]<=Nmissingthresh])
        else:
            skippedgenes.append(gene)
    return {'goodpdbs': chosenpdbs, 'goodgenes': array(chosengenes), 'Nmissing': Nmissing, 'skippedgenes': array(skippedgenes)}


def hamming(seqpairs, namepairs=None, lopass=L+1, L=L, alnpath=None, PSAdir=None): #unnormalized
    hamhits = []
    pdbhits = []
    genepairs = []
    seqidhits = []
    if alnpath != None:
        msaparse = MSAparse(alnpath=alnpath, L=L) 
        MSA = msaparse['dots']

    for n,(a,b) in enumerate(seqpairs):
        if PSAdir != None:
            psadat = getPSA(a,b, PSAdir=PSAdir)
            M1 = psadat[0][:,0]
            M2 = psadat[0][:,1]
            assert(len(M1) == len(M2))
            M1upper = empty(len(M1), dtype=bool)
            M2upper = empty(len(M1), dtype=bool)

            for n,(A,B) in enumerate(zip(M1,M2)):
                M1upper[n], M2upper[n] = (A.isupper(), B.isupper())

            seq_a = M1[M1upper&M2upper]
            seq_b = M2[M1upper&M2upper]
            L = len(seq_a)
            assert(len(seq_a)==len(seq_b))

       
        if alnpath != None:
            m1 = array(list(MSA[a]))
            m2 = array(list(MSA[b]))
            dotpairs = []
            assert(len(m1) == len(m2))
            for a,b in zip(m1,m2):
                dotpairs.append((a==b) and (a=='.'))
            dotpairs = array(dotpairs).astype(bool) #insertions in other parts of the MSA

            M1 = m1[~dotpairs]
            M2 = m2[~dotpairs]  #two sequences from the MSA aligned together, with insertions

            M1upper = empty(len(M1), dtype=bool)
            M2upper = empty(len(M1), dtype=bool)

            for n,(A,B) in enumerate(zip(M1,M2)):
                M1upper[n], M2upper[n] = (A.isupper(), B.isupper())

            seq_a = M1[M1upper&M2upper]
            seq_b = M2[M1upper&M2upper]

        mismatch = (seq_a!=seq_b).astype(int)
        hamming = sum(mismatch)
        seqid = 1-(hamming/float(L))
        if hamming < lopass:
            hamhits.append(hamming)
            pdbhits.append([a,b])
            seqidhits.append(seqid)
            if all(namepairs!=None):
                namehits.append(namepairs[n])

    print len(hamhits)

    return array(hamhits), array(pdbhits), array(seqidhits), array(namepairs)

def hamfilt(superensemble, targets, lopass=L+1, L=175):
    seqpairs = array([(a,b) for a in superensemble for b in targets])
    hamhits, pdbhits, genepairs = hamming(seqpairs, lopass=lopass, L=L, rids=rids, rgenes=rgenes, seqs=seqs, seqids=seqids)
    pdboptions = unique(pdbhits[:,0])
    iterpass = 0
    while True:
        iterpass += 1
        PDBpairinds = array([(a,b) for a in range(len(pdboptions)-1) for b in range(a+1, len(pdboptions))])
        hamhits, pdbhits, genepairs = hamming(pdboptions[PDBpairinds], lopass=lopass, L=L, rids=rids, rgenes=rgenes, seqs=seqs, seqids=seqids)
        genefailscore = []
        for name in unique(genepairs[:,0]):
            genefailscore.append(sum(hamhits[genepairs[:,0]==name]>=lopass)/float(sum(genepairs[:,0]==name)))
        genefailscore = array(genefailscore)
        if all(genefailscore==0):
            print 'Done'
            break
        else:
            print 'Pass '+str(iterpass)
        genesrt = argsort(genefailscore)[::-1]
        print genefailscore[genesrt]
        assert(len(genesrt) == len(unique(genepairs[:,0])))
        trimgenes = unique(genepairs[:,0])[genesrt][1:]
        pdboptions = unique(pdbhits[:,0][isin(genepairs[:,0], trimgenes)])
        
    return array(unique(pdbhits[:,0])), hamhits, pdbhits, genepairs




def SE(pdbpairs, namepairs=[], SEpath=SEpath, pdbpath=os.path.join(pdbpath,pdbprefix)):
    if len(namepairs)==0:
        namepairs = pdbpairs
    for pdb, name in zip(pdbpairs, namepairs):
        call('{}/se -fasta {}{}.pdb {}{}.pdb'.format(SEpath, pdbpath, pdb[0], pdbpath, pdb[1]).split())
        call('mv {}/{}{}-{}{}.fasta {}/{}-{}.fasta'.format(SEpath, pdbprefix, pdb[0], pdbprefix, pdb[1], SEpath, name[0].replace('/', '_'), name[1].replace('/','_')).split())

def get_resids(pdbdata, pdbsplit=True, pdbdata_chain=None):
    p = pdbdata
    if pdbsplit:
        chains = list(set(p.chain))
        chain = chains[0]
    else:
        chain = pdbdata_chain
    c = p[(p.chain == chain)]
    CA = c[c.atom == ' CA ']
    CB = c[(c.atom == ' CB ') | ((c.atom == ' CA ') & (c.resName == 'GLY'))] #for CB calculation, count CA as CB for gly

	#remove unknown residues and remove duplicate conformers
	#Important that this is done in exactly the same way as in 
	#getPDBseq.py so that "seen" residue indexes match

    fidsA = CA['resID'] #work with full resID to account for inserts
    filtA = (CA.resName != 'UNK') & (r_[True, fidsA[1:] != fidsA[:-1]])

    fidsB = CB['resID'] #work with full resID to account for inserts
    filtB = (CB.resName != 'UNK') & (r_[True, fidsB[1:] != fidsB[:-1]])

    resids = CA['resID'][filtA]

    nres = len(resids)
    dat = {'chain':c, 'CA':CA, 'CB':CB, 'fidsA':fidsA, 'filtA':filtA, 'fidsB':fidsB, 'filtB':filtB, 'resids':resids, 'nres': nres}
    return dat

def alnseq(msaparse, pdb):
    seq = msaparse['dots'][pdb]
    seq = array(list(seq))
    alnmsk = array([(seq[n].isupper() or seq[n]=='-') for n in range(len(seq))])
    aln = seq[alnmsk]
    return aln

def getalnseqs(msaparse=None, namelist=None, pdblist=['2F4J_A', '4LRM_A']):
    seqs = []
    validids = []
    validnames = []
    verbosenames = True
    if all(namelist==None):
        verbosenames = False
    for n,pdb in enumerate(pdblist):
        try:
            seqs.append(alnseq(msaparse, pdb))
            validids.append(pdb)
            if verbosenames == True:
                validnames.append(namelist[n])
        except:
            print 'Cannot align '+pdb
            continue
    dat = {'seqs':array(seqs), 'seqIDs':array(validids)}
    if verbosenames == True:
        dat['seqnames'] = array(namelist)
    return dat

def getPSA(seqAname, seqBname, PSAdir='../../../Alignment/Kinase/dfgin/se2.02/', headerlength=15):
    try:
        seqAind = 0 #indicates which sequence comes first in the PSA
        seqBind = 1
        fn = os.path.join(PSAdir, seqAname+'-'+seqBname+'.fasta')
        f = open(fn)
    except:
        seqBind = 0 #indicates which sequence comes first in the PSA
        seqAind = 1
        fn = os.path.join(PSAdir, seqBname+'-'+seqAname+'.fasta')
        f = open(fn)  

    lines = f.read()
    f.close()
    charlist = array(list(lines))
    charlist = charlist[charlist!='\n']
    seqAstart = findall(charlist=='>')[seqAind] + headerlength
    seqBstart = findall(charlist=='>')[seqBind] + headerlength
    pdbA = charlist[seqAstart-headerlength+9 : seqAstart]
    pdbB = charlist[seqBstart-headerlength+9 : seqBstart]

    if seqAind == 0:
        seqA = charlist[seqAstart : seqBstart-headerlength] #seqAname comes first in the PSA, so need to stop at the header of seqB
        seqB = charlist[seqBstart :]
        pdb1 = pdbA
        pdb2 = pdbB
        seq1 = seqA
        seq2 = seqB
    elif seqAind == 1:
        seqB = charlist[seqBstart : seqAstart-headerlength] #seqBname comes first in the PSA, so need to stop at the header of seqA
        seqA = charlist[seqAstart :]
        pdb1 = pdbB
        pdb2 = pdbA
        seq1 = seqB
        seq2 = seqA
    
    assert len(seqA) == len(seqB), 'Sequences in '+fn+' are not the same length'
    matches = zeros(len(seq1), dtype=bool)
    for n,(k,l) in enumerate(zip(seq1,seq2)): #doesn't matter if seqA or seqB is listed first, here
        matches[n] = k.isupper() and l.isupper()
    return array(zip(seqA, seqB)), (pdbA.tostring(), pdbB.tostring()), matches

def alnbenchmark(seqAname, seqBname, PSAdir='../../../Alignment/Kinase/dfgin/se2.02/', alnpath_=None, pdbpath=MLSpath, headerlength=15, L=L, MSA_=None, indmaps_=None, pdb2alnmap_=None):
    if alnpath_ == None:
#        print 'MSA from '+alnpath
        msaparse = MSAparse(alnpath=alnpath, L=L)
        MSA = msaparse['dots']
        indmaps = msaparse['pdbinds']
        pdb2alnmap = msaparse['alnPDBinds']

    elif all(array([MSA_, indmaps_, pdb2alnmap_]) == None):
#        print 'MSA from '+alnpath_
        msaparse = MSAparse(alnpath=alnpath_, L=L)
        MSA = msaparse['dots']
        indmaps = msaparse['pdbinds']
        pdb2alnmap = msaparse['alnPDBinds']

    elif all(array([MSA_, indmaps_, pdb2alnmap_]) != None):
        MSA = MSA_
        indmaps = indmaps_
        pdb2alnmap = pdb2alnmap_

    else:
        raise 'Invalid arguments'

    PSA = getPSA(seqAname, seqBname, PSAdir=PSAdir, headerlength=headerlength)
    pdb1 = PSA[1][0] #seqAname's pdb ID
    pdb2 = PSA[1][1]
    P1 = PSA[0][:,0]
    P2 = PSA[0][:,1]
    assert (len(pdb2alnmap[pdb1]) == L) and (len(pdb2alnmap[pdb2]) == L), str(pdb1)+','+str(pdb2)

######## Query MSA pre-processing ########   

    m1 = array(list(MSA[pdb1]))
    pdbinds1 = indmaps[pdb1]
    m2 = array(list(MSA[pdb2]))
    pdbinds2 = indmaps[pdb2]
    dotpairs = []
    assert(len(m1) == len(m2))
    for a,b in zip(m1,m2):
        dotpairs.append((a==b) and (a=='.'))
    dotpairs = array(dotpairs).astype(bool) #insertions in other parts of the MSA

    M1 = m1[~dotpairs]
    M2 = m2[~dotpairs]  #two sequences from the MSA aligned together, with insertions
    
    M1alpha = empty(len(M1), dtype=bool)
    M2alpha = empty(len(M1), dtype=bool)

    M1upper = empty(len(M1), dtype=bool)
    M2upper = empty(len(M1), dtype=bool)

    M1lower = empty(len(M1), dtype=bool)
    M2lower = empty(len(M1), dtype=bool)

    M1gaps = empty(len(M1), dtype=bool)
    M2gaps = empty(len(M1), dtype=bool)

    M1dots = empty(len(M1), dtype=bool)
    M2dots = empty(len(M1), dtype=bool)


    for n,(a,b) in enumerate(zip(M1,M2)):
#        print a,b
        M1alpha[n], M2alpha[n] = (a.isalpha(), b.isalpha())
        M1upper[n], M2upper[n] = (a.isupper(), b.isupper())
        M1lower[n], M2lower[n] = (a.islower(), b.islower())
        M1gaps[n], M2gaps[n] = (a=='-', b=='-')
        M1dots[n], M2dots[n] = (a=='.', b=='.')
        
    M1pdbinds = empty(len(M1), dtype='S5')
    M2pdbinds = empty(len(M1), dtype='S5')

    M1pdbinds[:] = '.'
    M2pdbinds[:] = '.'

    M1pdbinds[M1alpha] = pdbinds1
    M2pdbinds[M2alpha] = pdbinds2

    M1pdbinds[M1gaps] = '-'
    M2pdbinds[M2gaps] = '-'


######## Pairwise alignment "gold-standard" pre-processing #########

    P1alpha = empty(len(P1), dtype=bool)
    P2alpha = empty(len(P2), dtype=bool)

    P1upper = empty(len(P1), dtype=bool)
    P2upper = empty(len(P1), dtype=bool)

    P1lower = empty(len(P1), dtype=bool)
    P2lower = empty(len(P1), dtype=bool)

    P1gaps = empty(len(P1), dtype=bool)
    P2gaps = empty(len(P1), dtype=bool)

    for n,(a,b) in enumerate(zip(P1,P2)):
        P1alpha[n], P2alpha[n] = (a.isalpha(), b.isalpha())
        P1upper[n], P2upper[n] = (a.isupper(), b.isupper())
        P1lower[n], P2lower[n] = (a.islower(), b.islower())
        P1gaps[n], P2gaps[n] = (a=='-', b=='-')

    P1pdbinds = empty(len(P1), dtype='S5')
    P2pdbinds = empty(len(P1), dtype='S5')
 
    P1pdbinds[:] = '-'
    P2pdbinds[:] = '-'

    P1pdbinds[P1alpha] = pdbinds1[pdbinds1>-1]
    P2pdbinds[P2alpha] = pdbinds2[pdbinds2>-1]



#################  Benchmarking  ####################

    rids1 = get_resids(pdbParse.loadPDB(os.path.join(pdbpath, 'theseus_'+pdb1+'.pdb')))['resids']
    rids2 = get_resids(pdbParse.loadPDB(os.path.join(pdbpath, 'theseus_'+pdb2+'.pdb')))['resids']

    rids1 = rids1.astype(int)
    rids2 = rids2.astype(int)

    P1resinds = zeros_like(P1pdbinds, dtype='S4')
    P1resinds[P1alpha] = rids1

    P2resinds = zeros_like(P2pdbinds, dtype='S4')
    P2resinds[P2alpha] = rids2



# "Go back to the Shadow! You cannot pass. - Gandalf"
    assert(sum(logical_or(logical_or(M1upper, M2upper), M1gaps & M2gaps)) == L)
    assert(all((P1pdbinds=='-')==(P1=='-')))
    assert(all((P2pdbinds=='-')==(P2=='-')))
    assert(all((M1pdbinds=='-')==(M1=='-')))
    assert(all((M1pdbinds=='.')==(M1=='.')))
    assert(all((M2pdbinds=='-')==(M2=='-')))
    assert(all((M2pdbinds=='.')==(M2=='.')))
    assert(sum(pdbinds1>-1) == len(rids1))
    assert(sum(pdbinds2>-1) == len(rids2))


# Limit our analysis to the stretch of sequence spanned by both alignments

    start = findall(logical_or(M1upper, M1gaps))[0]
    end = findall(logical_or(M1upper, M1gaps))[-1]
# Generate a sequence of T, F, U spanned by the query and benchmark alignment
    CM = zeros(len(M1), dtype='S1')
    
    for n,(a,b) in enumerate(zip(M1pdbinds, M2pdbinds)):
        if (n < start) or (n > end):
            continue

        if M1upper[n] and M2upper[n]:
            m = findall(P1pdbinds==a)
            if (P1pdbinds[m] == a) and (P2pdbinds[m] == b):
                CM[n] = 'T'
            
            elif (P1upper[m] and P2upper[m]): #05/14/19
                CM[n] = 'F'

        elif any([M1gaps[n], M1lower[n], M1dots[n], M2gaps[n], M2dots[n]]):
            CM[n] = 'F'
        
    posCM = CM[logical_or(logical_or(M1upper, M2upper), M1gaps & M2gaps)]
    posCM = CM[logical_or(M1upper, M1gaps)]
    assert len(posCM) == L, posCM

    pdbpair = (pdb1, pdb2)    
    
    return {'CM':CM, 'pdbpair': pdbpair, 'posCM': posCM, 'goldinds': vstack([P1resinds, P2resinds]).T, 'goldmatches': vstack([P1upper, P2upper]).T}
    
def alnStats(pdbpath=os.path.join(pdbpath, pdbprefix), alnpath_=None, PSAdir=SEpath, headerlength=1+len(pdbprefix)+6, L=L, MSA_=None, indmaps_=None, pdb2alnmap_=None):
    files = glob.glob(os.path.join(PSAdir, '*.fasta'))
    T = empty(len(files), dtype=int)
    F = empty(len(files), dtype=int)
    U = empty(len(files), dtype=int)
    posCM = empty((len(files),L), dtype='S2')
    pdbpairs = empty((len(files),2), dtype='S6')
    namepairs = empty((len(files),2), dtype='S6')
    poshits = zeros(L, dtype=int)
    count = 0
    for i,filename in enumerate(files):
        fasta = os.path.split(filename)[-1]
        seq1name, seq2name = fasta[:-6].split('-')
        namepairs[i,:] = (seq1name, seq2name)
        try:
            dat = alnbenchmark(seq1name, seq2name, PSAdir=PSAdir, alnpath_=alnpath_, pdbpath=pdbpath, headerlength=headerlength, L=L, MSA_=MSA_, indmaps_=indmaps_, pdb2alnmap_=pdb2alnmap_)
            CM = dat['CM']
            T[i] = sum(CM=='T')
            F[i] = sum(CM=='F')
            U[i] = sum(CM=='')
            posCM[i,:] = dat['posCM']
            pdbpairs[i,:] = dat['pdbpair']
            count += 1
        except:
            T[i] = -1
            F[i] = -1
            U[i] = -1
            posCM[i,:] = ''
            pdbpairs[i] = ''
            continue

    print str(count)+' benchmark alignments'
    return {'T': T, 'F': F, 'posCM': posCM, 'pdbpairs': pdbpairs, 'namepairs': namepairs, 'benchmarks': count}


def report(alnstats, alnpath=None, PSAdir=SEpath):
    #def hamming(seqpairs, lopass=L+1, L=L, rids=rIDs, alnpath=None, PSAdir=None): #unnormalized

    nanmsk = alnstats['T']!=-1

    seqdat = hamming(alnstats['namepairs'][nanmsk,:],lopass=1000, alnpath=alnpath, PSAdir=PSAdir1)
    hamming = seqdat[0]
    validids = seqdat[1]
    seqidentity = seqdat[2]

    T = alnstats['T'][nanmsk]
    F = alnstats['F'][nanmsk]
    t = sum(alnstats['posCM'][nanmsk,:]=='T', axis=1)
    f = sum(alnstats['posCM'][nanmsk,:]=='F', axis=1)
    u = sum(alnstats['posCM'][nanmsk,:]=='', axis=1)
    alnu = sum(alnstats['posCM'][nanmsk,:]=='', axis=0)
    alnt = sum(alnstats['posCM'][nanmsk, :]=='T', axis=0)
    alnf = sum(alnstats['posCM'][nanmsk, :]=='F', axis=0)
    ACC = T.astype(float)/(T+F)
    acc = t.astype(float)/(t+f)
    alnacc = alnt.astype(float)/(alnt+alnf)
        
    stats = {'T':T, 't':t, 'F':F, 'f': f, 'alnt':alnt, 'alnf': alnf, 'u':u, 'alnu': alnu, 'ACC': ACC, 'acc': acc, 'alnacc': alnacc, 'seqidentity': seqidentity, 'validids': validids}
 
    return stats


