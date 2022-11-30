#!/usr/bin/env python

import sys
import argparse
from run import runCistromeGO

def prepareParser():
    """
    :return: argparse object
    """
    description = "Cistrome-GO --- Enrichment analysis of ChIP-seq peaks"
    epilog = "For command line options , type: %(prog)s -h"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    version = "1.0.0"
    parser.add_argument("-v", "--version", action="version", version="%(prog)s " + version)
    parser.add_argument('-g', '--genome', dest='assembly', choices=['mm9', 'hg19', 'mm10', 'hg38'], required=True, help='genome')
    parser.add_argument('-b', '--bed', dest='bed', type=str, required=True, help='Peak bed file path')
    parser.add_argument('-pn', '--peaknumber', dest='peakn', type=str, required=False, default="10000", help='Peak number to use. Default: 10000.')
    parser.add_argument('-d', '--decay', dest='decay', type=str, required=False, default='auto', help='Half decay distance to calculate RP score. Default: auto.')
    parser.add_argument('-n', '--name', dest='name', type=str, required=False, default='cistromego', help='Prefix of output file. Default: cistromego')
    parser.add_argument('-o', '--output', dest='output', type=str, required=False, default="./", help='Output directory. Default: ./')
    parser.add_argument('-dg', '--dego', dest='dego', help='upgenes, downgenes, allgenes for gene ontology. Default: allgenes.', default="allgenes")
    parser.add_argument('-e', '--expr', dest='expr', required=False, help='Expression file path')
    parser.add_argument('-ei', '--exprinfo', dest='exprinfo', type=str, required=False, default='1,3,7',
                        help='Expression file format setting. Column index of geneID, logFoldChange, FDR in differential expression file. Default: 1,3,7.')
    parser.add_argument('-max', '--maxgenenumber', dest='max_gene_number', type=int, required=False, default=2000,
                        help='maximum gene number in the GO or kegg terms. Default: 2000.')
    parser.add_argument('-min', '--mingenenumber', dest='min_gene_number', type=int, required=False, default=10,
                        help='minimum gene number in the GO or kegg terms. Default: 10.')
    parser.add_argument('-logfc', '--logfccut', dest='logfc_cut', type=float, required=False, default=1.0,
                        help='logFoldChange of differential expression genes. Default: 1.0.')
    parser.add_argument('-fdr', '--fdrcut', dest='fdr_cut', type=float, required=False, default=0.05,
                        help='FDR cutoff of differential expression genes. Default: 0.05. ')
    args = parser.parse_args()
    return args


def main():
    parser = prepareParser()
    runCistromeGO(parser)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrunpt me! ;-) Bye!\n")
        sys.exit(0)