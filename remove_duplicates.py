#! /usr/bin/env python
"""
A wrapper to identify, extract, align and phylogenise a target gene. The primary purpose
of this module is for a user to confirm orthology of their gene, not to understand in
detail the phylogenetic relationships between all closely related proteins.
"""

import argparse
import datetime
import os
import re
import tempfile

import matplotlib.pyplot as plt

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
    def __init__(self, transcript_list):
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
    parser.add_argument("-b", "--blast", type=str,
                        help="""blast results file (outfmt 6) to assign genes to
                        transcripts""")

    return parser

def parse_transcript(defline, blast_idx=None):
    """
    assumes a Trinity-styled + TransDecoder transcript ID,
    and extracts the gene and isoform id tags, plus positional information.
    """
    # get columns with relevant information:
    cols = defline.split()
    base = cols[0]
    posn = cols[-1]

    # search for Trinity gene and isoform ID:
    base_s = re.search( '>(\w+)\|?', base)

    # get TransDecoder peptide positions from Trinity transcript:
    tdpos_s = re.search( ':(\d+)-(\d+)', posn)

    # get TransDecoder gene id and blast result:
    tdgid = base[1:]
    if blast_idx and tdgid in blast_idx:
        blast = blast_idx[tdgid]
    else:
        blast = None

    if base_s and tdpos_s:
        # gene and isoform
        el = base_s.group(1).split('_')
        gene = "_".join(el[:-1])
        isoform = base

        # TD positions
        start = int(tdpos_s.group(1))
        end   = int(tdpos_s.group(2))

        return gene, isoform, start, end, blast

    else:
        return None, None, None, None, None

def get_blast_idx(blastfile, level='transdecoder'):
    taxalist = bb.get_taxa_groups(ants=True)
    blast_out = {}
    handle = open(args.blast, 'rb')
    for line in handle:
        if level == 'gene':
            trinity = "_".join(line.split()[0].split("|")[0].split("_")[:-1])
        elif level == 'transdecoder':
            trinity = line.split()[0]
        else:
            trinity = line.split()[0].split("|")[0]
        species, blast = bb.parse_gene(line, taxalist)
        if trinity in blast_out:
            blast_out[trinity].append(blast)
        else:
            blast_out[trinity] = [blast]

    return blast_out

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
    verbalise("Y",
        "%d gene families with >1 member" % (sum(1 for gf in gene_families if len(gf) > 1)))
    verbalise("Y",
        "%d gene families with >1 blast result" % (sum(1 for gf in gene_families if len(gf.blast_hits) > 1)))

    print ""

    # sanity check output:
    for gf in gene_families.keys()[:num]:
        if handle:
            handle.write("%s\n" % gf)
        verbalise("M", gf)
        for t in gf:
            if handle:
                handle.write("%s\n" % t)
            verbalise("B", t)

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

if __name__ == '__main__':
    dbpaths = config.import_paths()

    parser = define_arguments()
    args = parser.parse_args()

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output)

    temp_dir = tempfile.mkdtemp()
    os.rmdir(temp_dir)  # dir must be empty!

    # collate genes with identical sequences:
    seqdic_isoform = {}###
    seqdic_gene = {}###
    seqdic_plain = {}###
    gene_isoforms = {}###
    all_genes = {}###
    blast_genes = {}###
    blast_results = {}###

    transcript_dic = {}
    gene_families = {}
    seq_families = {}
    geneid_idx = {}
    gf_idx = {}

    blast_idx = get_blast_idx(args.blast, level='transdecoder')
    verbalise("Y", "Created blast index with %d entries" % (len(blast_idx)))

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


    # assemble gene_families starting with sequence similarity
    for ts in seq_families.values():
        newfamily = Gene_family(ts)

        # files for indexing
        gene_families[newfamily] = True
        for t in ts:
            if t in gf_idx:
                gf_idx[t].append(newfamily)
            else:
                gf_idx[t] = [newfamily]

    report(gene_families, 5, verbalise=verbalise)
    verbalise("Y", "%d trinity gene ids found" % (len(geneid_idx)))



    # merge families based on trinity gene name:
    trinity_pool = {} # index on members to removed duplicate genefamilies
    verbalise("B", "\nMerging datasets")
    total = len(gene_families)
    counter = 0
    t0 = datetime.datetime.now()
    for gf in gene_families.keys():
        for gid in gf.geneids.keys():
            if gid in geneid_idx:
                for t in geneid_idx[gid]:
                    if t in gf_idx:
                        for mgf in gf_idx[t]:
                            if mgf == gf:
                                pass
                            else:
                                gf.merge(mgf)

        trinity_pool[tuple(gf.members)] = gf
        counter = progress(counter, total, t0)

    trinity_dic = { v:True for v in trinity_pool.values()  }

    report(trinity_dic, 5, verbalise=verbalise)

    # write gene families to file:
    outhandle = open(logfile[:-3] + "families.out", 'w')
    report(trinity_dic, None, outhandle)
    outhandle.close()

    plt.hist([len(gf) for gf in trinity_dic], bins=range(12))
    plt.title("Number of genes in each gene family")
    plt.show()





