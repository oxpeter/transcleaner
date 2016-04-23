#! /usr/bin/env python
"""
A wrapper to identify, extract, align and phylogenise a target gene. The primary purpose
of this module is for a user to confirm orthology of their gene, not to understand in
detail the phylogenetic relationships between all closely related proteins.
"""

import argparse
import datetime
import itertools
import os
import re
import tempfile

import matplotlib.pyplot as plt
import networkx
from networkx.algorithms.components.connected import connected_components

from ortholotree import config
from orthomods import internal

import blast_breakdown as bb

############################################################################

class Transcript(object):
    """
    This class holds all the relevant details about a transcript
    """
    def __init__(self, trinity_id, transdecoder_id, blastline=None, seq=None):
        self.trinity_id      = trinity_id
        self.td_id           = transdecoder_id

        if blastline:
            cols = blastline.split()
            self.blast      = cols[1]
            self.start      = int(cols[6]) # the trinity pep alignment start position
            self.end        = int(cols[7]) # the trinity pep alignment end position
            self.bstart     = int(cols[8]) # the target protein start position
            self.bend       = int(cols[9]) # the target protein end position
        else:
            self.blast  = None
            self.start  = None
            self.end    = None
            self.bstart = None
            self.bend   = None

        self.seq             = seq

        self.geneid          = "_".join(trinity_id.split("_")[:-1])
        self.iso_id          = trinity_id

    def __str__(self):
        if self.seq:
            s1 = "...".join([self.seq[:5], self.seq[-5:]])
        else:
            s1 = "---"

        if self.blast:
            b1 = self.blast
        else:
            b1 = "---"

        if self.start and self.end:
            pos = "%s-%s" % (self.start, self.end)
        else:
            pos = ""

        if self.bstart and self.bend:
            bpos = "%s-%s" % (self.bstart, self.bend)
        else:
            bpos = ""

        return "%s [%s %s %s] (%s)" % (self.td_id, pos, b1, bpos, s1 )


    def __repr__(self):
        return "%r %r %r %r %r" % (self.trinity_id,
                                    self.td_id,
                                    self.blast,
                                    self.start,
                                    self.end)

    def __len__(self):
        if self.seq:
            return len(self.seq)
        else:
            return 0

class Gene_family(object):
    """
    This class defines a set of transcripts grouped into a "gene" based on several
    different criteria. Primarily this class allows fast indexing of all the
    transcript members, criteria members, and merging of two gene_families.
    """
    _ids = itertools.count(0) # to give a unique id to each instance

    def __init__(self, transcript_list):
        self.id = self._ids.next()

        self.members = { t:True for t in transcript_list }
        self.memberid = { t.trinity_id:True for t in transcript_list }
        self.blast_hits = { t.blast:True for t in transcript_list if t.blast is not None }
        self.sequences = { t.seq:True for t in transcript_list if t.seq is not None }
        self.geneids = { t.geneid:True for t in transcript_list }

    def __str__(self):
        return "%d members (%d w blast)" % (len(self.members), len(self.blast_hits))

    def __repr__(self):
        return "%r" % (self.members)

    def __len__(self):
        return len(self.members)

    def __contains__(self, item):
        if item in self.members:
            return True
        else:
            return False

    def __iter__(self):
        for t in self.members:
            yield t

    def __comp__(self, other):
        if self.members == other.members:
            return True
        else:
            return False

    def add_transcript(self, transcript):
        self.members[transcript] = True
        self.memberid[transcript.trinity_id] = True
        self.geneids[transcript.geneid] = True
        if transcript.blast:
            self.blast_hits[transcript.blast] = True
        if transcript.seq:
            self.sequences[transcript.seq] = True

    def merge(self, gf):
        for t in gf.members:
            self.add_transcript(t)

############################################################################

def define_arguments():
    parser = argparse.ArgumentParser(description=
            "A module to identify and remove transcript duplicates from fasta files")
    ### input options ###
    # logging options:
    parser.add_argument("-q", "--quiet", action='store_true',default=False,
                        help="print fewer messages and output details")
    parser.add_argument("-o", "--output", type=str, default='remove_duplicates',
                        help="specify the filename to save results to")
    parser.add_argument("-d", "--directory", type=str,
                        help="specify the directory to save results to")
    parser.add_argument("-D", "--display_on", action='store_true',default=False,
                        help="display graph results (eg for p value calculation)")


    # data file options:
    parser.add_argument("-t", "--transdecoder", type=str,
                        help="""Transdecoder peptide fasta file
                        (eg Trinity.fasta.transdecoder.pep) to analyse""")
    parser.add_argument("-f", "--fasta", type=str,
                        help="""fasta file (eg Trinity.fasta) to analyse""")
    parser.add_argument("-g", "--gff", type=str,
                        help="""gff file to enable new gff to be constructed from
                        combined gene families.""")
    parser.add_argument("-b", "--blast", type=str,
                        help="""blast results file (outfmt 6) to assign genes to
                        transcripts""")
    parser.add_argument("-S", "--stats", action='store_true',
                        help="""perform basic stats on the fasta file (N50 etc)""")
    parser.add_argument("-L", "--longest", action='store_true',
                        help="""extract the longest transcript from each gene, and
                        save a fasta file of the peptide sequences.""")

    return parser

def parse_defline(line):
    tdid = line.split()[0][1:]
    trinityid = tdid.split("|")[0]
    return tdid, trinityid

def get_full_blast_idx(blastfile):
    idx = {}
    handle = open(blastfile, 'rb')
    for line in handle:
        idx[line.split()[0]] = line.strip()

    return idx

def report(gene_families, num=5, handle=None, verbalise=lambda *x:None):

    verbalise("G", "%d transcripts in %d gene families" % (sum(len(gf) for gf in gene_families),
                                                len(gene_families) ))
    verbalise("G",
        "%d gene families with >1 member" % (sum(1 for gf in gene_families if len(gf) > 1)))
    verbalise("G",
        "%d gene families with >1 blast result" % (sum(1 for gf in gene_families if len(gf.blast_hits) > 1)))

    print ""

    # sanity check output:
    for gcount, gf in enumerate(gene_families.keys()[:num]):
        if handle:
            handle.write("%s\n" % gf)
        verbalise("M", gf)
        for t in gf:
            if handle:
                handle.write("%d %s\n" % (gcount, t))
            verbalise("B", t)

def make_gff(outhandle, gff, gene_families, gf_idx):
    no_family = {}
    no_transcript = {}
    handle = open(gff, 'rb')
    for line in handle:
        cols = line.split()
        if len(cols) >= 8:
            fields = cols[:8]
            atts = " ".join(cols[8:])
            if fields[2] == 'exon': # exon includes five_prime_utr, cds, and three_prime_utr
                ts_s = re.search('ID=([\w\|\.]+).exon\d+;', atts)
                if ts_s:
                    ts = ts_s.group(1)
                    if ts in gf_idx:
                        gf = gf_idx[ts]
                        newatts = ['gene_id "SGF%d"; transcript_id "%s";\n' % (gf.id, ts)]
                        newline = '\t'.join(fields + newatts)
                        outhandle.write(newline)
                    else:
                        no_family[ts_s.group(1)] = line.strip()
                else:
                    no_transcript[line.strip()] = line.strip()

    if len(no_family) > 0:
        verbalise("M", "%d transcripts found in gff that are not in the gf_idx" % len(no_family))
        verbalise("R", "\n".join(str(k)+" "+str(v) for k,v in no_family.items()[:5]))
    if len(no_transcript)>0:
        verbalise("M", "%d transcripts from gff not matched by regex profile" % len(no_transcript))



def progress(counter, total, t0):
    counter += 1
    if counter % int(1.0 * total / 10000) == 0:
        t1 = datetime.datetime.now()
        diff = t1-t0
        if diff.seconds < 2:
            left = 9999999999
        else:
            rate = 1.0 *  diff.seconds / counter
            left =  (total * rate) - diff.seconds

        try:
            remaining = datetime.timedelta(seconds=left)
        except OverflowError:
            remaining = datetime.timedelta(days=999999999)

        print '\r%d of %d complete (%.2f%%) Time remaining: %s       \r' % (counter,
                total, 100.0*counter/total, remaining),

    return counter

def to_graph(l):
    """
    This solution was modified from a stack overflow answer:
    http://stackoverflow.com/questions/4842613/merge-lists-that-share-common-elements
    """
    G = networkx.Graph()
    for gf in l:
        # each sublist is a bunch of nodes
        G.add_nodes_from(gf.geneids.keys())
        # it also imlies a number of edges:
        G.add_edges_from(to_edges(gf.geneids.keys()))
    return G

def to_edges(l):
    """
        treat `l` as a Graph and returns it's edges
        to_edges(['a','b','c','d']) -> [(a,b), (b,c),(c,d)]
    """
    it = iter(l)
    last = next(it)

    for current in it:
        yield last, current
        last = current

def fasta_stats(fastafile, display=False):
    # get basic stats on the fasta file:
    genelens = sorted( [ len(seq[1]) for seq in internal.parsefasta(fastafile)] , reverse=True)
    totallen = sum(genelens)

    print "Number of transcripts: %d" % len(genelens)
    print "Total length of transcripts: %d" % totallen
    print "longest transcript: %d" % genelens[0]

    if display:
        plt.hist(genelens)
        plt.title("Histogram of gene lengths (bp)")
        plt.show()
    prevl = 0
    lensum = 0
    trancount = 0
    for l in genelens:
        if lensum + l >= (totallen / 2):
            print "N50: %d (%d transcripts)" % (prevl, trancount)
            break
        lensum += l
        prevl = l
        trancount += 1

if __name__ == '__main__':
    dbpaths = config.import_paths()

    parser = define_arguments()
    args = parser.parse_args()

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output)

    temp_dir = tempfile.mkdtemp()
    os.rmdir(temp_dir)  # dir must be empty!

    # if requested, report the basic stats about the supplied fasta file:
    if args.stats and args.fasta:
        fasta_stats(args.fasta, args.display_on)

    #############################################################################
    ### extract sequence information from each gene in supplied transdecoder file
    transcript_dic = {}
    gene_families = {}
    seq_families = {}
    geneid_idx = {}
    gf_idx = {}


    full_blast = get_full_blast_idx(args.blast)
    verbalise("Y", "Created full blast index with %d entries" % (len(full_blast)))

    for defline, seq in internal.parsefasta(args.transdecoder):
        # get trinity and transdecoder gene ids:
        tdid, trinityid = parse_defline(defline)

        # get any blast results:
        if tdid in full_blast:
            blastline = full_blast[tdid]
        else:
            blastline = None

        # create new transcript instance
        newtranscript = Transcript(trinity_id=trinityid,
                                    transdecoder_id=tdid,
                                    blastline=blastline,
                                    seq=seq)

        # file based on sequence similarity
        if seq in seq_families:
            seq_families[seq].append(newtranscript)
        else:
            seq_families[seq] = [newtranscript]

        # create index of geneids for faster merging:
        if newtranscript.geneid in geneid_idx:
            geneid_idx[newtranscript.geneid].append(newtranscript)
        else:
            geneid_idx[newtranscript.geneid] = [newtranscript]

    verbalise("Y", "%d trinity gene ids\n\n" % (len(geneid_idx)))

    ############################################################
    ### assemble gene_families starting with sequence similarity
    for ts in seq_families.values():
        newfamily = Gene_family(ts)

        # files for indexing
        gene_families[newfamily] = True
        for t in ts:
            if t in gf_idx:
                gf_idx[t].append(newfamily)
            else:
                gf_idx[t] = [newfamily]
    verbalise("Y", "Results after merging based on sequence identity:")
    report(gene_families, 5, verbalise=verbalise)

    #############################################################################
    # merge gene families based on Trinity gene groups (solution from SO):
    G = to_graph(gene_families)

    # collect new gene families:
    trinity_pool = {}
    gf_idx = {} # reset the gene family index to the new groups
    count = 0

    for group in connected_components(G):
        count += 1
        ts = []
        for gid in group:
            ts += geneid_idx[gid]
        newgf = Gene_family(ts)

        # set index to find gf for given transcript:
        for t in newgf:
            gf_idx[t.td_id] = newgf

        trinity_pool[newgf] = True

    # report results:
    verbalise("Y", "\n\nResults after merging based on Trinity gene assignment:")
    report(trinity_pool, 5, verbalise=verbalise)

    ##############################
    # write gene families to file:
    if args.gff:
        outhandle = open(logfile[:-3] + "families.gtf", 'w')
        make_gff(outhandle, args.gff, trinity_pool, gf_idx)
        outhandle.close()

    if args.display_on:
        plt.hist([len(gf) for gf in trinity_pool], bins=range(12))
        plt.title("Number of genes in each gene family")
        plt.ylabel("Number of families")
        plt.xlabel("Number of trinity genes")
        plt.savefig(logfile[:-3] + "family_size.pdf", format='pdf')
        plt.show()


    #############################################################################
    # get longest peptide for each gene (not necessarily the longest transcript)
    if args.longest:







