#! /usr/bin/env python
"""
A wrapper to identify, extract, align and phylogenise a target gene. The primary purpose
of this module is for a user to confirm orthology of their gene, not to understand in
detail the phylogenetic relationships between all closely related proteins.
"""

import argparse
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
        if self.start and self.end:
            return "%s [%s-%s]" % (self.td_id, self.start, self.end)
        else:
            return "%s" % (self.td_id)

    def __repr__(self):
        return "%r %r %r %r %r" % (self.trinity_id,
                                    self.td_id,
                                    self.blast,
                                    self.start,
                                    self.end)


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

    def add_transcript(self, transcript):
        self.members[transcript] = True
        self.memberid[transcript.trinity_id] = True
        if transcript.blast:
            self.blast_hits[transcript.blast] = True
        if transcript.seq:
            self.sequences[transcript.seq] = True

    def merge(self, gf):
        for t in gf.members:
            self.add_transcript(t)
        del gf



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
    blast_results = {}##3
    transcript_dic = {}

    blast_idx = get_blast_idx(args.blast, level='transdecoder')
    verbalise("Y", "Created blast index with %d entries" % (len(blast_idx)))

    full_blast = get_full_blast_idx(args.blast)
    verbalise("Y", "Created full blast index with %d entries" % (len(full_blast)))

    for defline, seq in internal.parsefasta(args.transdecoder):
        gene, isoform, start, end, blast = parse_transcript(defline, blast_idx)
        #newtranscript = Transcript(trinity_id=, transdecoder_id=, blastline=, seq=)

        if seq in seqdic_gene:
            seqdic_gene[seq].append((gene, start, end, blast))
            seqdic_isoform[seq].append(isoform)
            seqdic_plain[seq].append(gene)
        else:
            seqdic_gene[seq] = [(gene, start, end, blast),]
            seqdic_isoform[seq] = [isoform,]
            seqdic_plain[seq] = [gene,]

        if gene in gene_isoforms:
            gene_isoforms[gene].append(isoform)
            all_genes[gene].append(seq)
        else:
            gene_isoforms[gene] = [isoform,]
            all_genes[gene] = [seq,]

        if blast and tuple(blast) in blast_results:
            blast_results[tuple(blast)].append(gene)
        elif blast:
            blast_results[tuple(blast)] = [gene]

        if blast and (gene,tuple(blast)) in blast_genes:
            blast_genes[(gene,tuple(blast))].append(seq)
        elif blast:
            blast_genes[(gene,tuple(blast))] = [seq,]

    verbalise("G", "%d unique sequences found" % (len(seqdic_gene)))
    verbalise("G", "%d unique trinity gene ids found" % (len(all_genes)))
    verbalise("Y",
        "%d blast results found, %d map to more than 1 gene" % (len(blast_results),
                                    sum( 1 for i in blast_results.values() if len(i) > 1)))

    #verbalise("B", blast_results.keys())
    genesets = []
    unclassified_genes = { k:1 for k in all_genes.keys()[:] }

    for gene in all_genes.keys()[:]:
        newset = [gene]
        for seq in all_genes[gene]:
            for g2 in seqdic_plain[seq]:
                newset.append(g2)

        genesets.append(tuple(set(newset)))

    genesets = list(set(genesets))

    verbalise("M",
        "%d sets with more than one trinity gene" % (sum([ 1 for gl in seqdic_plain.values() if len(set(gl)) > 1])))

    # report some stats:
    verbalise("G", "%d gene families created" % (len(genesets)))
    verbalise("C", "\n".join([ str(g) for g in genesets[:5]]))
    print ""



    plt.hist([len(set(i)) for i in genesets], bins=range(12))
    plt.title("Number of genes in each gene family")
    plt.show()

    plt.hist([len(set(i)) for i in seqdic_gene.values()],
                    bins=range(12), facecolor='b', alpha=0.7)
    plt.hist([ len(set(i)) for i in seqdic_isoform.values()],
                    bins=range(12), facecolor='g', alpha=0.7)
    plt.show()





