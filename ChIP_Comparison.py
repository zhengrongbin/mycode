import os,sys
import collections
import pandas as pd
import numpy as np
import subprocess
import argparse
from optparse import OptionParser

def get_parameter():
    # peaks_files = sys.argv[1] ## a set of peak files, comma separated
    # bigwig_files = sys.argv[2] ## a set of bigwig files, comma separated, correspond with peaks_files
    # bam_files = sys.argv[3]

    parser = argparse.ArgumentParser(description="""ChIP-seq downstream analysis""")

    parser.add_argument('-p', dest='peak_files', type=str, required=True, 
        help='path of peak files, comma separated')
    parser.add_argument('-lp', dest='peakfile_label', type=str, required=False, 
        help='label for peak files, comma separated, should be same with ordering of peak files given by -p')

    parser.add_argument('-bw', dest='bw_files', type=str, required=False, 
        help='path of bigwig files, comma separated')
    parser.add_argument('-lbw', dest='bwfile_label', type=str, required=False, 
        help='label for bigwig files, comma separated, should be same with ordering of bigwig files given by -p')

    parser.add_argument('-bm', dest='bam_files', type=str, required=False, 
        help='path of BAM files, dedup BAM recommanded, comma separated')
    parser.add_argument('-lbm', dest='bamfile_label', type=str, required=False, 
        help='label for bam files, comma separated, should be same with ordering of bam files given by -p')

    parser.add_argument('-d', dest='diff_design', type=str, required=False, 
        help='path of design matrix to conduct differential peak analysis, the sample name should be same as what given by -lbm')
    parser.add_argument('-k', dest='kmeans', type=str, required=False, 
        help='k clusters for signal heatmap, default is 3')

    parser.add_argument('-f5', dest='fold5', action = 'store_true',
        help='by giving this parameters, only 5 fold peaks were used in all analyses')

    parser.add_argument('-O', dest='output_dir', type=str, required=False, 
        help='the path of output derectory, default current working path')
    parser.add_argument('-n', dest='prefix', type=str, required=False, 
        help='prefix of file names, default NA')

    args = parser.parse_args()
    pk_files = args.peak_files
    pk_files = pk_files.split(',')
    pfile_label = args.peakfile_label
    pfile_label = pfile_label.split(',')
    bw_files = args.bw_files
    bw_files = bw_files.split(',')
    bwfile_label = args.bwfile_label
    bwfile_label = bwfile_label.split(',')
    bm_files = args.bam_files
    bm_files = bm_files.split(',')
    bamfile_label = args.bamfile_label
    bamfile_label = bamfile_label.split(',')
    diffdesign = args.diff_design
    kmeans = args.kmeans
    f5 = args.fold5
    output_dir = args.output_dir if args.output_dir is not None else "./"
    prefix = args.prefix+'_' if args.prefix is not None else ""
    parameters = {'peak_files':pk_files, 'peakfile_label':pfile_label,
                    'bw_files': bw_files, 'bw_label':bwfile_label,
                    'bam_files': bm_files, 'bam_label': bamfile_label,
                    'design_matrix':diffdesign, 'kmeans': kmeans, 'fold5': f5,
                    'output_dir': output_dir, 'prefix': prefix}
    return(parameters)

# if len(pk_files) != len(bw_files) or len(pk_files) != len(bm_files):
#     print('bigwig, bam, and peaks files with diff number')
#     sys.exit(0)
def _peakStat(parameters):
    ## read peak file and stat peak numbers
    print('++ stat peak number')
    pstats = {}
    pk_files = parameters['peak_files']
    for p in pk_files:
        pstats[p] = {}
        ## total
        total = subprocess.getoutput("cat "+p+"|grep -w -E 'chr[0-9]|chrX|chrY'|wc -l")
        total = int(total)
        pstats[p]['TotalPeak'] = total
        ## 5 fold
        fold5 = subprocess.getoutput("awk '{if($7>5){print}}' "+p+"|grep -w -E 'chr[0-9]|chrX|chrY'|wc -l")
        fold5 = int(fold5)
        pstats[p]['5FoldPeak'] = fold5
        ## 10 fold peaks
        fold10 = subprocess.getoutput("awk '{if($7>10){print}}' "+p+"|grep -w -E 'chr[0-9]|chrX|chrY'|wc -l")
        fold10 = int(fold10)
        pstats[p]['10FoldPeak'] = fold10
        ## 20 fold peaks
        fold20 = subprocess.getoutput("awk '{if($7>20){print}}' "+p+"|grep -w -E 'chr[0-9]|chrX|chrY'|wc -l")
        fold20 = int(fold20)
        pstats[p]['20FoldPeak'] = fold20

    pstats_df = pd.DataFrame(pstats).T
    pstats_df.to_csv(os.path.join(parameters['output_dir'], parameters['prefix']+'peak_stats.txt'), sep = '\t')
    return(parameters)

def _getUnionPeak(parameters):
    ## get union peaks and dedup
    print('++ union peaks')
    pk_files = parameters['peak_files']
    f5 = parameters['fold5']
    uPeakOut = os.path.join(parameters['output_dir'], parameters['prefix']+"union_dedup_peaks.bed")
    parameters['unionPeak'] = uPeakOut
    if f5:
        print('++++ 5 fold peak used')
        cmd1 = "cat %s |awk '{if($7>5){print}}' | bedtools sort |bedtools merge |grep -w -E 'chr[0-9]|chrX|chrY'> %s"%(' '.join(pk_files), parameters['unionPeak'])
    else:
        cmd1 = "cat %s | bedtools sort |bedtools merge |grep -w -E 'chr[0-9]|chrX|chrY'> %s"%(' '.join(pk_files), parameters['unionPeak'])
    os.system(cmd1)
    return(parameters)

def _bwHeatmap(parameters):
    ## make bigwig heatmap
    print('++ bw heatmap')
    bw_files = parameters['bw_files']
    parameters['bwSignalOutput'] = os.path.join(parameters['output_dir'], parameters['prefix']+'union_dedup_peaks_3k.mat.gz')
    if not os.path.exists(parameters['bwSignalOutput']):
        cmd3 = "computeMatrix reference-point -S %s -R %s -a 3000 -b 3000 --outFileName %s --numberOfProcessors 8"%(' '.join(bw_files),
            parameters['unionPeak'], parameters['bwSignalOutput'])
        os.system(cmd3)
    ## plot heatmap
    cmd4 = "plotHeatmap -m %s -out %s --colorMap RdBu --refPointLabel 'center' --whatToShow 'heatmap and colorbar' --kmeans %s --samplesLabel %s"%(parameters['bwSignalOutput'],
        os.path.join(parameters['output_dir'], parameters['prefix']+"union_dedup_peaks_3k.png"),
        parameters['kmeans'], ' '.join(parameters['bw_label']))
    os.system(cmd4)
    return(parameters)

def _cleanStat(stat):
    stat_dict = {}
    for x in stat.split('\n'):
        x = x.split(' ')
        n = int(x[0])
        label = ' '.join(x[3:])
        stat_dict[label] = n
    return(stat_dict)

def _bamCount(parameters):
    ## make count matrix with bam files
    print('++ bam count')
    ## build index
    bm_files = parameters['bam_files']
    bam_label = parameters['bam_label']
    for b in bm_files:
        if os.path.exists(os.path.join(b, '.bai')):
            continue
        os.system('samtools index %s'%b)
          
    parameters['bamCountOutput'] = os.path.join(parameters['output_dir'], parameters['prefix']+'union_dedup_peaks_count_mat.txt')
    
    ## bam stat
    parameters['BamMapStat'] = os.path.join(parameters['output_dir'], parameters['prefix']+'mapping_stat.txt')
    if not parameters['BamMapStat']:
        print('++ bam mapping stat')
        stats = {}
        for b in bm_files:
            ## get stat
            stat = subprocess.getoutput("samtools flagstat %s"%b)
            stats[bam_label[bm_files.index(b)]] = _cleanStat(stat)
      
        stats = pd.DataFrame(stats).T
        stats.to_csv(parameters['BamMapStat'], sep = '\t')
    ## extract count matrix
    if not os.pathe.exists(parameters['bamCountOutput']):
        # return (parameters)
        print('++ bam read count')
        cmd2 = "bedtools multicov -bams %s -bed %s > %s"%(' '.join(bm_files), parameters['unionPeak'], parameters['bamCountOutput'])
        os.system(cmd2)
        ## add column names
        bcount = pd.read_csv(parameters['bamCountOutput'], sep = '\t')
        bcount.columns = ['chr', 'start', 'end']+bam_label
        bcount.to_csv(parameters['bamCountOutput'], sep = '\t')

    return(parameters)

def _diffPeakAnalysis(parameters):
    cmd = 'Rscript diffPeak.R -c ./allPeak/Par_ELF3_ER_union_dedup_peaks_count_mat.txt -m design.csv -s ./allPeak/Par_ELF3_ER_mapping_stat.txt -e hsa -p Par_ELF3_ER -o ./allPeak'

def main():
    parameters = get_parameter()
    print(parameters)
    parameters = _peakStat(parameters)
    parameters = _getUnionPeak(parameters)
    parameters = _bwHeatmap(parameters)
    parameters = _bamCount(parameters)
    print(parameters)

main()







