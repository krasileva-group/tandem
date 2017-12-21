#!/usr/bin/env python
import os
import sys
import csv
import argparse

from collections import Counter, namedtuple

GFeature = namedtuple('GFeature', 'seqname source feature start end score strand frame attribute'.split(' '))
Tandem = namedtuple('Tandem', 'seqname nlr1 start1 end1 strand1 parent1 nlr2 start2 end2 strand2 parent2 type delta'.split(' '))

VERBOSE = False

def getAttributeDict(attribute):
    return dict(item.split('=') for item in attribute.split(';'))
def getID(feature):
    return feature.attribute.get('transcript_id', feature.attribute.get('Name', 'NA'))

def readGFF(fn):
    with open(fn) as fi:
        for row in csv.reader(fi, delimiter='\t'):
            if not row[0].startswith('#'):
                row[3:5] = list(map(int, row[3:5]))
                row[8] = getAttributeDict(row[8])
                yield GFeature(*row)

def extractFeature(fn, ftype='mRNA'):
    for feature in readGFF(fn):
        if feature.feature == ftype:
            yield feature

def tandemCheck(set1, set2, geneset, maxdist=15000, maxoverlap=2000):
    """ checks two sets of GFF-records (transcript-clusters) for Tandem-architecture """

    def makeTandem(nlr1, nlr2, ttype, delta):
        """ takes two (NLR) GFF-records and returns Tandem-namedtuples """
        def getParent(attribute):
            """ returns parent information for GFF (mRNA) records """
            return attribute.get('Parent', attribute.get('gene_id', 'NA'))

        return Tandem(nlr1.seqname, getID(nlr1), nlr1.start, nlr1.end, nlr1.strand, getParent(nlr1.attribute), getID(nlr2), nlr2.start, nlr2.end, nlr2.strand, getParent(nlr2.attribute), ttype, delta)

    def mergeTranscripts(transcripts):
        """ 
        helper function to create dummy representative transcripts from transcript clusters, created to circumvent S. italica mess
        """
        start = min(t.start for t in transcripts)
        end = max(t.end for t in transcripts)
        attrDict = transcripts[0].attribute
        attrDict['Name'] = attrDict['transcript_id'] = ':'.join(map(getID, transcripts))
        return [GFeature(transcripts[0].seqname, transcripts[0].source, transcripts[0].feature, start, end, '.', transcripts[0].strand, '.', attrDict)]

    """ test which/how many transcripts in the cluster are of interest according to specified gene set """
    ids1 = list(map(getID, set1))
    ids2 = list(map(getID, set2))
    if VERBOSE:
        print('SET1:', ids1)
        print('SET2:', ids2)
    targets1 = list(set(ids1).intersection(geneset))
    targets2 = list(set(ids2).intersection(geneset))

    """ both clusters have to have at least one transcript in the gene set """
    if targets1 and targets2:
        if VERBOSE:
            print('T1:', set1, sep='\n')
            print('T2:', set2, sep='\n')
            print('------')
        """ S. italica mess: multiple representative transcripts: use dummy information, was corrected in input later on, can probably be disabled now  """
        if len(targets1) > 1:
            set1 = mergeTranscripts(set1)
        if len(targets2) > 1:
            set2 = mergeTranscripts(set2)
        if len(targets1) > 1 or len(targets2) > 1:
            if VERBOSE:
                print('MERGED T1:', set1, sep='\n')
                print('MERGED T2:', set2, sep='\n')
                print('------')

        """ pull out the representatives according to gene set """
        nlr1 = [item for item in set1 if getID(item) == targets1[0]][0]
        nlr2 = [item for item in set2 if getID(item) == targets2[0]][0]
        # nlr2 = set1[0], set2[0]

        """ assert that both NLRs are on the same sequence (chr, scaff, ...) """
        if nlr1.seqname == nlr2.seqname:
            """ tandem-architecture: minus-transcript is upstream of plus-transcript """
            if nlr1.strand == '-' and nlr2.strand == '+':
                """ remove pairs with full overlaps """
                if nlr1.start <= nlr2.start <= nlr2.end <= nlr1.end:
                    return tuple() # nlr1, nlr2, 'DEVOURED'
                if nlr2.start <= nlr1.start <= nlr1.end <= nlr2.end:
                    return tuple() # nlr1, nlr2, 'DEVOURED'
                """ overlap """
                if nlr1.start <= nlr2.start <= nlr1.end <= nlr2.end:
                    d = nlr1.end - nlr2.start
                    if d <= maxoverlap:
                        return makeTandem(nlr1, nlr2, 'OVL-TANDEM', d)
                    return tuple() # nlr1, nlr2, 'OVERLAP TOO BIG'
                d = nlr2.start - nlr1.end
                """ non-overlapping tandem """
                if d <= maxdist:
                    return makeTandem(nlr1, nlr2, 'TANDEM', d)
                return tuple() # nlr1, nlr2, 'TOO FAR APART'
    return tuple()

def tandemScan(fn, targets, maxdist=15000, maxoverlap=2000):
    """ scans mRNA features from ORDERED! (by start position) GFF for tandems """
    features = extractFeature(fn)
    targetset = set(targets)

    try:
        transcripts1 = [next(features)]
    except:
        return None

    """ method: 
         - cluster all transcripts from same gene/locus G1 together
         - cluster all transcripts from next gene/locus G2 together
         - check for tandem between cluster representatives (according to specified gene 'whitelist')
         - pop G1, G2 -> G!
         - build next cluster -> G2, then check for tandem ... repeat until done
         This was necessary, because original dataset contained conflicting transcripts, now probably can be done on representative transcripts.
    """
    gene_locus1 = transcripts1[0].attribute.get('Parent', 'NA')
    gene_locus2, transcripts2 = '', list()

    for transcript in features:
        current_locus = transcript.attribute.get('Parent', 'NA')
        if current_locus == gene_locus1:
            transcripts1.append(transcript)
        else:
            if gene_locus2:
                if VERBOSE:
                    print(gene_locus2)
                    print(*(t.attribute['ID'] for t in transcripts2), sep=',')
                    print('------')

                tandem = tandemCheck(transcripts2, transcripts1, targetset, maxdist=maxdist, maxoverlap=maxoverlap)
                if tandem:
                    ann1, ann2 = targets.get(tandem.nlr1, TargetGene(tandem.nlr1, 'NA', 'NA', 'NA')), targets.get(tandem.nlr2, TargetGene(tandem.nlr1, 'NA', 'NA', 'NA'))
                    print(*(tandem), *ann1[1:], *ann2[1:], sep='\t')

            gene_locus2, transcripts2 = gene_locus1, transcripts1
            gene_locus1, transcripts1 = current_locus, [transcript]

    if VERBOSE:
        print(gene_locus2)
        print('x', *(t.attribute['ID'] for t in transcripts2), sep=',')
        print('------')
        print(gene_locus1)
        print('x', *(t.attribute['ID'] for t in transcripts1), sep=',')
        print('------')

    tandem = tandemCheck(transcripts2, transcripts1, targetset, maxdist=maxdist, maxoverlap=maxoverlap)
    if tandem:
        ann1, ann2 = targets.get(tandem.nlr1, TargetGene(tandem.nlr1, 'NA', 'NA', 'NA')), targets.get(tandem.nlr2, TargetGene(tandem.nlr1, 'NA', 'NA', 'NA'))
        print(*(tandem), *ann1[1:], *ann2[1:], sep='\t')

TargetGene = namedtuple('TargetGene', 'id gclass gsubclass hasid'.split(' '))

def readTargets(fn):
    """
     read gene white list, only gene-pairs with both partners in this list are tested for tandem-ness (e.g. all NLRs as whitelist)
    """
    targets = dict()
    with open(fn) as fi:
        for row in csv.reader(fi, delimiter='\t'):
            # print(row)
            targets[row[0]] = TargetGene(row[0], 'TREE' if row[1] != 'TRUNC' else 'TRUNC', row[1], 'NLRID' if len(row) == 3 and row[2] == 'NLRID' else 'NLR')
    return targets





if __name__ == '__main__':

    ap = argparse.ArgumentParser()
    ap.add_argument('transcriptome', type=str, help='A sorted(!) GFF file (best performance with only mRNA records.)')
    ap.add_argument('target_genes', type=str, help='A list of transcript identifiers to be considered in the scan.')
    ap.add_argument('--max-distance', '-d', type=int, default=15000)
    ap.add_argument('--max-overlap', '-l', type=int, default=2000)
    args = ap.parse_args()

    assert os.path.exists(args.transcriptome)
    assert os.path.exists(args.target_genes)
    assert 0 < args.max_distance
    assert 0 < args.max_overlap

    # with open(args.target_genes) as fi:
    #    targets = set(line.strip() for line in fi)
    targets = readTargets(args.target_genes)

    tandemScan(args.transcriptome, targets, maxdist=args.max_distance, maxoverlap=args.max_overlap)
    pass
