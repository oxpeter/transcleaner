#! /usr/bin/env python
"""
A wrapper to identify, extract, align and phylogenise a target gene. The primary purpose
of this module is for a user to confirm orthology of their gene, not to understand in
detail the phylogenetic relationships between all closely related proteins.
"""

import argparse
import os
import tempfile

from ortholotree import config
from orthomods import internal, external

############################################################################

def define_arguments():
    parser = argparse.ArgumentParser(description=
            "A module to identify and remove transcript duplicates from fasta files")
    ### input options ###
    # logging options:
    parser.add_argument("-q", "--quiet", action='store_true',default=False,
                        help="print fewer messages and output details")
    parser.add_argument("-o", "--output", type=str, default='genematch.out',
                        help="specify the filename to save results to")
    parser.add_argument("-d", "--directory", type=str,
                        help="specify the directory to save results to")
    parser.add_argument("-D", "--display_on", action='store_true',default=False,
                        help="display graph results (eg for p value calculation)")


    # data file options:
    parser.add_argument("-f", "--fasta", type=str,
                        help="Fasta file to analyse")

    return parser

if __name__ == '__main__':
    dbpaths = config.import_paths()

    parser = define_arguments()
    args = parser.parse_args()

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output)

    temp_dir = tempfile.mkdtemp()
    os.rmdir(temp_dir)  # dir must be empty!

