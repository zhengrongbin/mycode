import os,sys
import pandas as pd
import numpy as np
import argparse
from argparse import RawDescriptionHelpFormatter


## === gene annotation file
global mm10_gene_ann_path; mm10_gene_ann_path = '/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/gencode.vM27.annotation.gene_annotation.csv'
global hg38_gene_ann_path; hg38_gene_ann_path = '/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/hg38/gencode.v38.annotation.gene_annotation.csv'

class FriendlyArgumentParser(argparse.ArgumentParser):
    """
    Override argparse to show full length help information
    """

    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(1)


def parse_args(args=None):
    """
    parse input arguments
    """
    parser = FriendlyArgumentParser(description=__doc__)
    parser.add_argument('-V', '--version', action='version',
                        version="%(prog)s (code version 1.0, db version 1.0)")

    sub_parsers = parser.add_subparsers(
        help="sub-command help", dest="sub_command")

    parser_ann_gene = sub_parsers.add_parser("ann_gene", help="annotate gene symbol",
                                          description="annotate gene symbol using ensembl ID")
    parser_merge_mat = sub_parsers.add_parser("merge_matrix", help="merge matrix",
                                         description="merge expression matrix by giving a path list")
    parser_merge_stat = sub_parsers.add_parser("merge_stat", help="merge stat",
                                         description="merge statstic matrix by giving a path list")

    # parser_corr arguments
    parser_ann_gene.add_argument("-c", "--count", dest="count_file", required=False, type=str,
                             help="path of STAR output read count")
    parser_ann_gene.add_argument("-t", "--tpm", dest="tpm_file", required=False, type=str,
                             help="path of RSEM output gene rsem")
    parser_ann_gene.add_argument("-s", "--species", dest="species", required=True, type=str,
                             help="should be one of [hg38, mm10]")
    parser_ann_gene.add_argument("-d", "--dir", dest="dir_path", required=False, type=str,
                             help="new directory for saving the result, default back to the original path")


    # evaluate the function of given genes using GSEA
    parser_merge_mat.add_argument('-i', '--inFile', dest='inFile_path', required=True,
                            help="file path, the files contains a column for STAR output count file path, should be annotated to gene symbol using ann_gene sub-command")
    parser_merge_mat.add_argument("-ft", "--ftype", dest="file_type", required=False,
                            help="should be one of [count, tpm, fpkm], indicating you are merging count matrix or tpm/fpkm matrix")
    parser_merge_mat.add_argument("-d", "--dir", dest="dir_path", required=False, type=str,
                             help="new directory for saving the result, default back to the original path")

    # evaluate the function of given genes using GSEA
    parser_merge_stat.add_argument('-i', '--inFile', dest='inFile_path', required=True,
                            help="file path, the files contains a column for STAR output statstic file path")
    parser_merge_stat.add_argument("-d", "--dir", dest="dir_path", required=False, type=str,
                             help="new directory for saving the result, default back to the original path")

    return (parser.parse_args(args), parser)


def _annotate_gene_name_(count_file, tpm_file, species, dir_path):
    """
    annotate gene name using esembl gene id
    """
    if species == 'hg38':
        gene_ann_path = hg38_gene_ann_path
    elif species == 'mm10':
        gene_ann_path = mm10_gene_ann_path
    else:
        print('species problem')
        sys.exit(0)
    ## read gene annotation
    gene_ann = pd.read_csv(gene_ann_path)
    gene_ann.index = gene_ann.gene_id.tolist()

    # count_file = '/lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/RNA-seq/adipocyte_CL_rep1/result/adipocyte_CL_rep1.gene_count.txt'
    # tpm_file = '/lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/RNA-seq/adipocyte_CL_rep1/result/adipocyte_CL_rep1.gene.rsem.txt'

    ## read count
    if count_file:
        count = pd.read_csv(count_file, sep = '\t', header = None)
        count = count[~count[0].isin(['N_unmapped', 'N_multimapping', 'N_noFeature', 'N_ambiguous'])]
        count = pd.merge(count, gene_ann[['gene_name']], left_on = 0, right_index = True)
        if dir_path:
            count.to_csv(os.path.join(dir_path, os.path.basename(count_file)), sep = '\t', header = None, index = None)
        else:
            count.to_csv(count_file, sep = '\t', header = None, index = None)
    ## TPM
    if tpm_file:
        TPM = pd.read_csv(tpm_file, sep = '\t')
        TPM = pd.merge(TPM, gene_ann[['gene_name']], left_on = 'gene_id', right_index = True)
        if dir_path:
            TPM.to_csv(os.path.join(dir_path, os.path.basename(tpm_file)), sep = '\t', index = None)
        else:
            TPM.to_csv(tpm_file, sep = '\t', index = None)

def _merge_matrix_(inFile_path, file_type, dir_path):
    """
    merge count, tpm, and fpkm result
    """
    path_set = dict()
    for path in [x.rstrip().split('\t')[0] for x in open(inFile_path)]:
        label = os.path.basename(path).replace(".gene.rsem.txt", '').replace('.isoforms.rsem.txt', '').replace('.gene_count.txt', '')
        path_set[label] = path
    ### iterate each file
    res = pd.DataFrame()
    if file_type == 'count':
        for label in path_set.keys():
            path = path_set[label]
            cont = pd.read_csv(path, sep = '\t', header = None)
            cont = cont.groupby(4)[1].max() ## unique gene
            res = pd.concat([res, pd.DataFrame(cont.tolist(), index = cont.index.tolist(), columns = [label])], axis = 1)
    elif file_type == 'tpm':
        for label in path_set.keys():
            path = path_set[label]
            cont = pd.read_csv(path, sep = '\t')
            cont = cont.groupby('gene_name')['TPM'].max() ## unique gene
            res = pd.concat([res, pd.DataFrame(cont.tolist(), index = cont.index.tolist(), columns = [label])], axis = 1)
    elif file_type == 'fpkm':
        for label in path_set.keys():
            path = path_set[label]
            cont = pd.read_csv(path, sep = '\t')
            cont = cont.groupby('gene_name')['FPKM'].max() ## unique gene
            res = pd.concat([res, pd.DataFrame(cont.tolist(), index = cont.index.tolist(), columns = [label])], axis = 1)
    else:
        print('data type problem')
        sys.exit(0)
    if dir_path:
        res.to_csv(os.path.join(dir_path, file_type+'_matrix.csv'))
    else:
        res.to_csv(file_type+'_matrix.csv')

def _star_stat_(inFile_path, dir_path):
    """
    merge stat into matrix
    """
    path_set = dict()
    for path in [x.rstrip().split('\t')[0] for x in open(inFile_path)]:
        label = os.path.basename(path).replace("_STAR.stat", '')
        path_set[label] = path
    # for path in [x.rstrip().split('\t')[0] for x in open(inFile_path)]:
    #     label = path.split('/')[0]
    #     path_set[label] = path
    # ## 
    res = pd.DataFrame()
    for label in path_set:
        path = path_set[label]
        cont = [x.strip().split('\t') for x in open(path).readlines()[5:] if not x.rstrip().endswith(':')] 
        cont = pd.DataFrame(cont)
        cont[0] = cont[0].str.rstrip(' |')
        cont.index = cont[0].tolist()
        cont.columns = ['term', label]
        res = pd.concat([res, cont[[label]]], axis = 1) 
    if dir_path:
        res.to_csv(os.path.join(dir_path, 'STAR_mapping_stat.csv'))
    else:
        res.to_csv('STAR_mapping_stat.csv')



def main(args=None):
    args, parser = parse_args(args)
    if not args.sub_command:
        sys.exit(1)
    if args.sub_command == 'ann_gene':
        _annotate_gene_name_(count_file=args.count_file, tpm_file=args.tpm_file, species=args.species, dir_path=args.dir_path)
    if args.sub_command == 'merge_matrix':
        _merge_matrix_(inFile_path=args.inFile_path, file_type=args.file_type, dir_path=args.dir_path)
    if args.sub_command == 'merge_stat':
        _star_stat_(inFile_path=args.inFile_path, dir_path=args.dir_path)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr("User interrupt:) \n")











