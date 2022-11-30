import os,sys
import collections
import pandas as pd
import numpy as np
import subprocess
import argparse
from optparse import OptionParser
import gzip
import json
import pybedtools
from adjustText import adjust_text
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams.update(plt.rcParamsDefault)
rc={"axes.labelsize": 16, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "figure.titleweight":"bold", #"font.size":14,
    "figure.figsize":(5.5,4.2), "font.weight":"regular", "legend.fontsize":10,
    'axes.labelpad':8, 'figure.dpi':300}
plt.rcParams.update(**rc)

def get_parameter():
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

    parser.add_argument('-dd', dest='diff_design_deseq', type=str, required=False, 
        help='path of design matrix to conduct differential peak analysis for DESeq2, the sample name should be same as what given by -lbm')
    parser.add_argument('-db', dest='diff_design_bedtools', type=str, required=False, 
        help='path of design matrix to conduct differential peak analysis for Bedtools, the sample name should be same as what given by -lp')

    parser.add_argument('-dm', dest='diff_method', type=str, required=False, 
        help='method for diff peaks, can be deseq2 or bedtools, or deseq2,bedtools')

    parser.add_argument('-k', dest='kmeans', type=str, required=False, 
        help='k clusters for signal heatmap, default is 3')

    parser.add_argument('-s', dest='species', type=str, required=False, 
        help='can be hg38, mm10')

    parser.add_argument('-f5', dest='fold5', action = 'store_true',
        help='by giving this parameters, only 5 fold peaks were used in all analyses')

    parser.add_argument('-O', dest='output_dir', type=str, required=False, 
        help='the path of output derectory, default current working path')
    parser.add_argument('-n', dest='prefix', type=str, required=False, 
        help='prefix of file names, default NA')

    parser.add_argument('-c', dest='colmap', type=str, required=False, 
        help='color map for plotHeatmap, default will be RdBu, can be any in matplotlib colmap, such as Reds_r, Blue_r')
    parser.add_argument('-decay', dest='decay_dist', type=str, required=False, 
        help="decay distance to calculate RP score for peaks, should be one of ['auto', '0.5k', '1k', '10k', '50k', '100k', '500k'], default will be 'auto'")
    parser.add_argument('-sf', dest='size_factor', type=str, required=False, 
        help="should be true or false, true to take total mapped read count as size factor in DESeq2, false for not.")


    args = parser.parse_args()
    pk_files = args.peak_files
    pk_files = pk_files.split(',')
    pfile_label = args.peakfile_label
    pfile_label = pfile_label.split(',')
    bw_files = args.bw_files
    bw_files = bw_files.split(',') if bw_files is not None or bw_files is not False  else None  
    bwfile_label = args.bwfile_label
    bwfile_label = bwfile_label.split(',') if bwfile_label is not None or bwfile_label is not False  else None 
    bm_files = args.bam_files
    bm_files = bm_files.split(',') if bm_files is not None or bm_files is not False  else None
    bamfile_label = args.bamfile_label
    bamfile_label = bamfile_label.split(',') if bamfile_label is not None or bamfile_label is not False  else None
    diffdesign = args.diff_design_deseq
    diffdesignBedtools = args.diff_design_bedtools
    diffmethod = args.diff_method
    diffmethod = diffmethod.split(',') if diffmethod is not None or diffmethod is not False else None
    kmeans = args.kmeans
    f5 = args.fold5
    species = args.species
    colMap = args.colmap
    colmap = 'RdBu' if colMap is None or colMap is False else colMap
    decay = args.decay_dist
    decay = 'auto' if decay is None or decay is False else decay
    sizef = args.size_factor 
    sizef = 'false' if sizef is None or sizef is False else sizef

    output_dir = args.output_dir if args.output_dir is not None or args.output_dir is not False else "./"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    prefix = args.prefix+'_' if args.prefix is not None or args.prefix is not False else ""
    parameters = {'peak_files':pk_files, 'peakfile_label':pfile_label,
                    'bw_files': bw_files, 'bw_label':bwfile_label,
                    'bam_files': bm_files, 'bam_label': bamfile_label,
                    'design_matrix':diffdesign, 'kmeans': kmeans, 'fold5': f5,
                    'output_dir': output_dir, 'prefix': prefix, 'species':species,
                    'diffmethod':diffmethod, 'diffdesignBedtools':diffdesignBedtools,
                    'colmap':colmap, 'decay':decay, 'sizefactor': sizef}
    return(parameters)

# if len(pk_files) != len(bw_files) or len(pk_files) != len(bm_files):
#     print('bigwig, bam, and peaks files with diff number')
#     sys.exit(0)

def _plot_peakstat(pstats_df, parameters):
    data = pstats_df.T.values

    columns = pstats_df.index
    rows = pstats_df.columns

    values = np.arange(0, 2500, 500)
    value_increment = 1000

    # Get some pastel shades for the colors
    color = {'TotalPeak':'lightgrey', '5FoldPeak':'#87CEFA',
                         '10FoldPeak':'#1E90FF', '20FoldPeak':'#4169E1'}
    colors = list(color.values())
    n_rows = len(data)

    index = np.arange(len(columns)) + 0.3
    bar_width = 0.4

    # Initialize the vertical-offset for the stacked bar chart.
    y_offset = np.zeros(len(columns))
    colors = colors[::-1]

    # Plot bars and create text labels for the table
    cell_text = []
    for row in range(n_rows):
        plt.bar(index, data[row], bar_width, bottom=y_offset, color=colors[row])
        y_offset = y_offset + data[row]
        cell_text.append([int(x) for x in y_offset])
    # Reverse colors and text labels to display the last value at the top.
    # cell_text.reverse()

    # Add a table at the bottom of the axes
    the_table = plt.table(cellText=cell_text,
                          rowLabels=rows,
                          rowColours=colors,
                          colLabels=columns,
                          loc='bottom')

    # Adjust layout to make room for the table:
    plt.subplots_adjust(left=0.2, bottom=0.2)

    plt.ylabel("# of peaks")
    plt.xticks([])
    plt.tight_layout()
    plt.savefig(os.path.join(parameters['output_dir'], parameters['prefix']+'peak_stats.pdf'))
    plt.close()

def _peakStat(parameters):
    ## read peak file and stat peak numbers
    print('++ stat peak number')
    pstats = {}
    pk_files = parameters['peak_files']
    pk_labels = parameters['peakfile_label']
    for p in pk_files:
        l = pk_labels[pk_files.index(p)]
        ## ceas
        pstats[l] = {}
        ## total
        total = subprocess.getoutput("cat "+p+"|grep -w -E 'chr[0-9]*|chrX|chrY'|awk '{if($9>2){print}}' |wc -l")
        total = int(total)
        pstats[l]['TotalPeak'] = total
        ## 5 fold
        fold5 = subprocess.getoutput("awk '{if(($7>5)&&($9>2)){print}}' "+p+"|grep -w -E 'chr[0-9]*|chrX|chrY'|wc -l")
        fold5 = int(fold5)
        pstats[l]['5FoldPeak'] = fold5
        ## 10 fold peaks
        fold10 = subprocess.getoutput("awk '{if(($7>10)&&($9>2)){print}}' "+p+"|grep -w -E 'chr[0-9]*|chrX|chrY'|wc -l")
        fold10 = int(fold10)
        pstats[l]['10FoldPeak'] = fold10
        ## 20 fold peaks
        fold20 = subprocess.getoutput("awk '{if(($7>20)&&($9>2)){print}}' "+p+"|grep -w -E 'chr[0-9]*|chrX|chrY'|wc -l")
        fold20 = int(fold20)
        pstats[l]['20FoldPeak'] = fold20

    pstats_df = pd.DataFrame(pstats).T
    pstats_df.to_csv(os.path.join(parameters['output_dir'], parameters['prefix']+'peak_stats.txt'), sep = '\t')
    ## plot
    pstats_df['TotalPeak'] = pstats_df['TotalPeak'] - pstats_df['5FoldPeak']
    pstats_df['5FoldPeak'] = pstats_df['5FoldPeak'] - pstats_df['10FoldPeak']
    pstats_df['10FoldPeak'] = pstats_df['10FoldPeak'] - pstats_df['20FoldPeak']
    pstats_df.index = pstats_df.index.str.replace('.narrowPeak', '')
    pstats_df = pstats_df[['20FoldPeak','10FoldPeak','5FoldPeak','TotalPeak']]
    _plot_peakstat(pstats_df, parameters)
    ## output good peaks
    for i in range(len(pk_labels)):
        p_path = pk_files[i]
        pl = pk_labels[i]
        if parameters['fold5']:
            cmd3 = "awk '{if(($7>5)&&($9>2)){print}}' %s |grep -w -E 'chr[0-9]*|chrX|chrY' > %s"%(
                p_path, os.path.join(parameters['output_dir'], pl+'.GoodPeak.bed'),
                )
        else:
            cmd3 = "awk '{if($9>2){print}}' %s |grep -w -E 'chr[0-9]*|chrX|chrY' > %s"%(
                p_path, os.path.join(parameters['output_dir'], pl+'.GoodPeak.bed')
                )
        os.system(cmd3)
        peakfile = os.path.join(parameters['output_dir'], pl+'.GoodPeak.bed')
        _ceas(peakfile, parameters)        
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
        cmd1 = "cat %s |awk '{if(($7>5)&&($9>2)){print}}' | bedtools sort |bedtools merge |grep -w -E 'chr[0-9]*|chrX|chrY'> %s"%(' '.join(pk_files), parameters['unionPeak'])
    else:
        cmd1 = "cat %s |awk '{if($9>2){print}}'| bedtools sort |bedtools merge |grep -w -E 'chr[0-9]*|chrX|chrY'> %s"%(' '.join(pk_files), parameters['unionPeak'])
    os.system(cmd1)
    return(parameters)

def _bwHeatmap(parameters):
    ## make bigwig heatmap
    print('++ bw heatmap')
    bw_files = parameters['bw_files']
    parameters['bwSignalOutput'] = os.path.join(parameters['output_dir'], parameters['prefix']+'union_dedup_peaks_3k.mat.gz')
    ## define peak center bed for heatmap
    pmat = pd.read_csv(parameters['unionPeak'], header = None, sep = '\t')
    pmat[1] = pmat[1] + ((pmat[2]-pmat[1])/2).astype('int')
    pmat[2] = pmat[1] + 1
    pmat.to_csv(parameters['unionPeak']+'.center.bed', header = None, index = None, sep = '\t')

    if not os.path.exists(parameters['bwSignalOutput']):
        cmd3 = "computeMatrix reference-point -S %s -R %s -a 3000 -b 3000 --outFileName %s --numberOfProcessors 8"%(' '.join(bw_files),
            parameters['unionPeak']+'.center.bed', parameters['bwSignalOutput'])
        os.system(cmd3)
    if not os.path.exists(os.path.join(parameters['output_dir'], parameters['prefix']+'union_dedup_peaks_3k.png')):
        ## plot heatmap
        cmd4 = "plotHeatmap -m %s -out %s --colorMap %s --refPointLabel 'center' --xAxisLabel 'distance (bp)' --numberOfProcessors 8 --yAxisLabel 'peaks' --whatToShow 'heatmap and colorbar' --kmeans %s --samplesLabel %s --outFileSortedRegions %s"%(
            parameters['bwSignalOutput'], 
            os.path.join(parameters['output_dir'], parameters['prefix']+'union_dedup_peaks_3k.png'),
            parameters['colmap'],
            parameters['kmeans'], 
            ' '.join(parameters['bw_label']), 
            os.path.join(parameters['output_dir'], parameters['prefix']+'union_dedup_peaks_3k.order.txt'))
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
        if not os.path.exists(b+'.bai'):
            print('+++ index bam')
            os.system('samtools index %s'%b)
    parameters['bamCountOutput'] = os.path.join(parameters['output_dir'], parameters['prefix']+'union_dedup_peaks_count_mat.txt')
    ## bam stat
    parameters['BamMapStat'] = os.path.join(parameters['output_dir'], parameters['prefix']+'mapping_stat.txt')
    if not os.path.exists(parameters['BamMapStat']):
        print('+++ bam stat')
        stats = {}
        for b in bm_files:
            ## get stat
            stat = subprocess.getoutput("samtools flagstat %s"%b)
            stats[bam_label[bm_files.index(b)]] = _cleanStat(stat)
        stats = pd.DataFrame(stats).T
        stats.to_csv(parameters['BamMapStat'], sep = '\t')
    ## extract count matrix
    if not os.path.exists(parameters['bamCountOutput']):
        # return (parameters)
        print('+++ bam count')
        cmd2 = "bedtools multicov -bams %s -bed %s > %s"%(' '.join(bm_files), parameters['unionPeak'], parameters['bamCountOutput'])
        os.system(cmd2)
        ## add column names
        bcount = pd.read_csv(parameters['bamCountOutput'], sep = '\t')
        bcount.columns = ['chr', 'start', 'end']+bam_label
        bcount.to_csv(parameters['bamCountOutput'], sep = '\t', index = None)
    return(parameters)


def _bwHeatmap_from_mat(peakbed, prefix, parameters, denovel=False):
    if denovel is False:
        peaks = []
        for x in open(peakbed, 'r'):
            x = x.rstrip().split('\t')
            peaks.append(x[0]+':'+x[1]+'-'+x[2])
        ## get subset of mat
        matpath = parameters['bwSignalOutput']
        regionLines = []
        for line in gzip.open(matpath, 'r'):
            lineutf8 = line.decode('utf8')
            if lineutf8.startswith('@'):
                titleLine = lineutf8
                # out.write(str.encode(line))
            else:
                if lineutf8.split('\t')[3] in peaks:
                    regionLines.append(lineutf8)
        i = titleLine.index('group_boundaries')
        j = titleLine.index('],"sample_labels')
        titleLine = titleLine.replace(titleLine[i:j], 'group_boundaries":[0,'+str(len(regionLines)))
        out = gzip.open('%s_peaks.mat.gz'%prefix, 'wb')
        out.write(str.encode(titleLine))
        for l in regionLines:
            out.write(str.encode(l))
        out.close()
        ## start to plot bwheatmap
        cmd = "plotHeatmap -m %s -out %s --colorMap %s --refPointLabel 'center' --xAxisLabel 'distance (bp)' --numberOfProcessors 8 --yAxisLabel 'peaks' --whatToShow 'heatmap and colorbar' --samplesLabel %s --outFileSortedRegions %s"%(
            '%s_peaks.mat.gz'%prefix, 
            prefix+'_peaks.png',
            parameters['colmap'],
            ' '.join(parameters['bw_label']), prefix+'.order.txt')
        os.system(cmd)
    else:
        bw_files = parameters['bw_files']
        if not os.path.exists('%s_peaks.mat.gz'%prefix):
            ## define peak center bed for heatmap
            pmat = pd.read_csv(peakbed, header = None, sep = '\t')
            pmat[1] = pmat[1] + ((pmat[2]-pmat[1])/2).astype('int')
            pmat[2] = pmat[1] + 1
            pmat.to_csv(peakbed+'.center.bed', header = None, index = None, sep = '\t')

            cmd1 = "computeMatrix reference-point -S %s -R %s -a 3000 -b 3000 --outFileName %s --numberOfProcessors 8"%(' '.join(bw_files),
                peakbed+'.center.bed', '%s_peaks.mat.gz'%prefix)
            os.system(cmd1)
        cmd = "plotHeatmap -m %s -out %s --colorMap %s --refPointLabel 'center' --xAxisLabel 'distance (bp)' --numberOfProcessors 8 --yAxisLabel 'peaks' --whatToShow 'heatmap and colorbar' --samplesLabel %s --outFileSortedRegions %s"%(
            '%s_peaks.mat.gz'%prefix, 
            prefix+'_peaks.png',
            parameters['colmap'],
            ' '.join(parameters['bw_label']), prefix+'.order.txt')
        os.system(cmd)


    # bw_files = parameters['bw_files']
    # mat_path = '%s_peaks.mat.gz'%prefix
    # if not os.path.exists(mat_path):
    #     cmd3 = "/lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/envs/chips/bin/computeMatrix reference-point -S %s -R %s -a 3000 -b 3000 --outFileName %s --numberOfProcessors 8"%(' '.join(bw_files),
    #         peakbed, '%s_peaks.mat.gz'%prefix)
    #     os.system(cmd3)
    # if not os.path.exists(prefix+'_peaks.png'):
    #     ## plot heatmap
    #     cmd4 = "/lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/envs/chips/bin/plotHeatmap -m %s -out %s --colorMap RdBu --refPointLabel 'center' --whatToShow 'heatmap and colorbar' --samplesLabel %s --outFileSortedRegions %s"%('%s_peaks.mat.gz'%prefix, prefix+'_peaks.png',
    #     ' '.join(parameters['bw_label']), prefix+'.order.txt')
    #     os.system(cmd4)

def _diff_deseq(parameters):
    if not os.path.exists(parameters['output_dir']+'/DESeq2'):
        os.mkdir(parameters['output_dir']+'/DESeq2')
    cmat_path = '{outdir}/{pref}union_dedup_peaks_count_mat.txt'.format(
        outdir = parameters['output_dir'], pref = parameters['prefix'])
    mapstat_path = '{outdir}/{pref}mapping_stat.txt '.format(
        outdir = parameters['output_dir'], pref = parameters['prefix'])
    # cmat_path_for_deseq = '{outdir}/{pref}union_dedup_peaks_count_mat.txt'.format(
    #     outdir = parameters['output_dir']+'/DESeq2', pref = parameters['prefix'])
    sizefactor = parameters['sizefactor']
    if sizefactor == 'true':
        cmd = 'Rscript diffPeak.R -c {cmat} -m {design} -s {mapstat} -p {pref} -o {outdir}'.format(
            cmat= cmat_path,
            mapstat = mapstat_path,
            outdir = parameters['output_dir']+'/DESeq2',
            pref = parameters['prefix'],
            design = parameters['design_matrix']
            )
    else:
        cmd = 'Rscript diffPeak.R -c {cmat} -m {design} -s {mapstat} -p {pref} -o {outdir}'.format(
            cmat= cmat_path,
            outdir = parameters['output_dir']+'/DESeq2',
            pref = parameters['prefix'],
            design = parameters['design_matrix']
            )
    os.system(cmd)



def venn_mpl(a, b, colors=None, 
             outfn="out.pdf", labels=None,
             dpi=300, figsize = (5, 4)):
    """
    *a*, *b*, and *c* are filenames to BED-like files.
    *colors* is a list of matplotlib colors for the Venn diagram circles.
    *outfn* is the resulting output file.  This is passed directly to
    fig.savefig(), so you can supply extensions of .png, .pdf, or whatever your
    matplotlib installation supports.
    *labels* is a list of labels to use for each of the files; by default the
    labels are ['a','b','c']
    
    *dpi* is the dpi setting passed to matplotlib savefig
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Circle
    except ImportError:
        sys.stderr.write(
            "matplotlib is required to make a Venn diagram with %s\n"
            % os.path.basename(sys.argv[0])
        )
        sys.exit(1)

    a = pybedtools.BedTool(a)
    b = pybedtools.BedTool(b)

    if colors is None:
        colors = ["r", "b"] #if c is not None else ["r", "b"]

    radius = 1.0
    center = 0.0
    offset = radius / 2
    s = sum([a.count(), b.count()]) #if c is None else sum([a.count(), b.count(), c.count()])
    
    if labels is None:
        labels = ["a", "b"] #if c is not None else ["a", "b"]
    aradius = radius
    bradius = aradius * b.count() / a.count()
    ab = (a + b).count()
    
    Oa = ab * aradius / a.count()
    Ob = ab * bradius / b.count()
    
    aoffset = aradius - Oa
    boffset = bradius - Ob 
    
    circle_a = Circle(
        xy=(center - aoffset, center),
        radius=aradius,
        edgecolor=colors[0],
        label=labels[0],
    )
    
    circle_b = Circle(
        xy=(center + boffset, center),
        radius=bradius,
        edgecolor=colors[1],
        label=labels[1],
    )
    
    fig = plt.figure(facecolor="w", figsize = figsize)
    ax = fig.add_subplot(111)

    for circle in (circle_a, circle_b):
        circle.set_facecolor("none")
        circle.set_linewidth(3)
        ax.add_patch(circle)


    ax.axis("tight")
    ax.axis("equal")
    ax.set_axis_off()

    kwargs = dict(horizontalalignment="center")

    
    atextset = aradius - 0.5*Oa
    btextset = bradius - 0.5*Ob 
    
    # Unique to A
    t1 = ax.text(center - atextset, center, str((a - b).count()), **kwargs)

    # Unique to B
    t2 = ax.text(center + btextset, center, str((b - a).count()), **kwargs)

    t3 = ax.text(
            center, center, str((a + b).count()), **kwargs
        )
    adjust_text([t1, t2, t3], arrowprops=dict(arrowstyle="-", lw=0.5), save_steps = False, **kwargs)
    
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 8)
    plt.tight_layout()
    fig.savefig(outfn, dpi=dpi)
    plt.close(fig)

def _diff_bedtools(parameters):
    peakbeds = parameters['peak_files']
    peak_labels = parameters['peakfile_label']

    design = pd.read_csv(parameters['diffdesignBedtools'], index_col = 0)
    comprisons = design.index.tolist()
    for comp in comprisons:
        ccomp = design.loc[comp].dropna()
        p1 = ccomp[ccomp == 1].index.tolist()
        if len(p1) != 1:
            print('design problem for bedtools, multiple 1 or no 1')
        else:
            p1 = p1[0]
        p2 = ccomp[ccomp == 0].index.tolist()
        if len(p2) != 1:
            print('design problem for bedtools, multiple 1 or no 1')
        else:
            p2 = p2[0]
        # p1_path = peakbeds[peak_labels.index(p1)]
        # p2_path = peakbeds[peak_labels.index(p2)]
        # ## select peaks, p1
        # if parameters['fold5']:
        #     cmd3 = "awk '{if(($7>5)&&($9>2)){print}}' %s |grep -w -E 'chr[0-9]*|chrX|chrY' > %s"%(
        #         p1_path, os.path.join(parameters['output_dir'], p1+'.GoodPeak.bed'),
        #         )
        # else:
        #     cmd3 = "awk '{if($9>2){print}}' %s |grep -w -E 'chr[0-9]*|chrX|chrY' > %s"%(
        #         p1_path, os.path.join(parameters['output_dir'], p1+'.GoodPeak.bed')
        #         )
        # os.system(cmd3)
        # ## p2
        # if parameters['fold5']:
        #     cmd4 = "awk '{if(($7>5)&&($9>2)){print}}' %s |grep -w -E 'chr[0-9]*|chrX|chrY' > %s"%(
        #         p2_path, os.path.join(parameters['output_dir'], p2+'.GoodPeak.bed'),
        #         )
        # else:
        #     cmd4 = "awk '{if($9>2){print}}' %s |grep -w -E 'chr[0-9]*|chrX|chrY' > %s"%(
        #         p2_path, os.path.join(parameters['output_dir'], p2+'.GoodPeak.bed'),
        #         )
        # os.system(cmd4)
        ## replace p1 and p2
        p1_path = os.path.join(parameters['output_dir'], p1+'.GoodPeak.bed')
        p2_path = os.path.join(parameters['output_dir'], p2+'.GoodPeak.bed')
        ## venn
        venn_mpl(a=p1_path, b=p2_path, colors=['b', 'g'], 
            outfn=os.path.join(parameters['output_dir']+'/bedtools', parameters['prefix']+comp+'_venn.pdf'), 
            labels=[p1, p2], dpi=300, figsize = (5, 4))

        ## diff by bedtools
        cmd1 = "bedtools intersect -a %s -b %s -v | sort|uniq > %s"%(
            p1_path, p2_path,
            os.path.join(parameters['output_dir']+'/bedtools', parameters['prefix']+comp+'_UPpeak.bed'))
        cmd2 = "bedtools intersect -a %s -b %s -v | sort|uniq > %s"%(
            p2_path, p1_path,
            os.path.join(parameters['output_dir']+'/bedtools', parameters['prefix']+comp+'_DOWNpeak.bed'))
        cmd5="bedtools intersect -a %s -b %s -wa | sort|uniq > %s"%(
            p2_path, p1_path,
            os.path.join(parameters['output_dir']+'/bedtools', parameters['prefix']+comp+'_SHAREDpeak.bed'))
        os.system(cmd1)
        os.system(cmd2)
        os.system(cmd5)

def _diffPeakAnalysis(parameters):
    peakPathCollect = {'DESeq2':{}, 'bedtools':{}}
    if 'deseq2' in parameters['diffmethod']:
        print('+++ diff by DESeq2')
        if not os.path.exists(os.path.join(parameters['output_dir']+'/DESeq2')):
            os.mkdir(os.path.join(parameters['output_dir']+'/DESeq2'))
        _diff_deseq(parameters)
        
        print('+++ plot heatmap and Cistrome-GO for diff by DESeq2')
        design = pd.read_csv(parameters['design_matrix'], index_col = 0)
        comprisons = design.index.tolist()
        for comp in comprisons:
            ## up and down peaks
            peakpath = {"%s_UPpeak"%comp:os.path.join(parameters['output_dir']+'/DESeq2', parameters['prefix']+comp+'_UPpeak.bed'),
                        "%s_DOWNpeak"%comp:os.path.join(parameters['output_dir']+'/DESeq2', parameters['prefix']+comp+'_DOWNpeak.bed'),
                        "%s_SHAREDpeak"%comp:os.path.join(parameters['output_dir']+'/DESeq2', parameters['prefix']+comp+'_SHAREDpeak.bed')}
            peakPathCollect['DESeq2'][comp] = peakpath
            for p in peakpath:
                peakbed = peakpath[p]
                linen = int(subprocess.getoutput("wc -l %s|cut -f1 -d ' '"%peakbed))
                if linen > 50:
                    ## ceas
                    _ceas(peakbed, parameters)
                    ## heatmap
                    print('+++ heatmap for %s'%p)
                    if not os.path.exists(os.path.join(parameters['output_dir']+'/DESeq2', parameters['prefix']+p)+'_peaks.png'):
                        _bwHeatmap_from_mat(peakbed, os.path.join(parameters['output_dir']+'/DESeq2', parameters['prefix']+p), parameters, True)
                    ## cistrome-go
                    print('+++ Cistrome-GO for %s'%p)
                    go_folder = '%s_go'%os.path.join(parameters['output_dir']+'/DESeq2', parameters['prefix']+p)
                    cmd1='mkdir -p %s'%go_folder
                    os.system(cmd1)
                    prefix_new = os.path.join(go_folder, p)
                    cmd = 'python cistromego_py3/cistromego.py -g %s -b %s -pn all -n %s -d %s'%(parameters['species'], peakbed, prefix_new, parameters['decay'])
                    os.system(cmd)
                else:
                    print('+++ skip heatmap for %s due to low peak number %s'%(p, linen))
    if 'bedtools' in parameters['diffmethod']:
        print('+++ diff by Bedtools')
        if not os.path.exists(os.path.join(parameters['output_dir']+'/bedtools')):
            os.mkdir(os.path.join(parameters['output_dir']+'/bedtools'))
        _diff_bedtools(parameters)
        ## ceas for good peaks
        peakbeds = parameters['peak_files']
        peak_labels = parameters['peakfile_label']
        design = pd.read_csv(parameters['diffdesignBedtools'], index_col = 0)
        comprisons = design.index.tolist()
        # for cond in design.columns.tolist():
        #     peakfile = os.path.join(parameters['output_dir'], cond+'.GoodPeak.bed')
        #     _ceas(peakfile, parameters)

        print('+++ plot heatmap and Cistrome-GO for diff by Bedtools')
        for comp in comprisons:
            ## up and down peaks
            peakpath = {"%s_UPpeak"%comp:os.path.join(parameters['output_dir']+'/bedtools', parameters['prefix']+comp+'_UPpeak.bed'),
                        "%s_DOWNpeak"%comp:os.path.join(parameters['output_dir']+'/bedtools', parameters['prefix']+comp+'_DOWNpeak.bed'),
                        "%s_SHAREDpeak"%comp:os.path.join(parameters['output_dir']+'/bedtools', parameters['prefix']+comp+'_SHAREDpeak.bed')}
            peakPathCollect['bedtools'][comp] = peakpath
            for p in peakpath:
                peakbed = peakpath[p]
                linen = int(subprocess.getoutput("wc -l %s|cut -f1 -d ' '"%peakbed))
                if linen > 50:
                    ## ceas
                    _ceas(peakbed, parameters)
                    ## heatmap
                    print('+++ heatmap for %s'%p)
                    if not os.path.exists(os.path.join(parameters['output_dir']+'/bedtools', parameters['prefix']+p)+'_peaks.png'):
                        if p.endswith('_SHAREDpeak'):
                            _bwHeatmap_from_mat(peakbed, os.path.join(parameters['output_dir']+'/bedtools', parameters['prefix']+p), parameters, True)
                        else:
                            _bwHeatmap_from_mat(peakbed, os.path.join(parameters['output_dir']+'/bedtools', parameters['prefix']+p), parameters, True)
                    ## cistrome-go
                    print('+++ Cistrome-GO for %s'%p)
                    go_folder = '%s_go'%os.path.join(parameters['output_dir']+'/bedtools', parameters['prefix']+p)
                    cmd1='mkdir -p %s'%go_folder
                    os.system(cmd1)
                    prefix_new = os.path.join(go_folder, p)
                    cmd = 'python cistromego_py3/cistromego.py -g %s -b %s -pn all -n %s -d %s'%(parameters['species'], peakbed, prefix_new, parameters['decay'])
                    os.system(cmd)
                else:
                    print('+++ skip heatmap for %s due to low peak number %s'%(p, linen))
def _ceas(peakbed, parameters):
    species = parameters['species']
    if species == 'hg38':
        genetable = 'ceas/geneTable/hg38.refGene'
    elif species == 'mm10':
        genetable = 'ceas/geneTable/mm10.refGene'
    else:
        print('species error')
        sys.exit(0)
    outdir = peakbed.replace(os.path.basename(peakbed), '')
    basename = os.path.basename(peakbed)
    cmd = 'python ceas/ceas_bedAnnotate.v2.py -g %s -b %s -o %s -n %s'%(genetable, peakbed, outdir, basename)
    os.system(cmd)

def _motif_homer(parameters):
    ## diff peak path
    print("++ motif homer")
    if 'deseq2' in parameters['diffmethod']:
        design = pd.read_csv(parameters['design_matrix'], index_col = 0)
        comprisons = design.index.tolist()
        for comp in comprisons:
            ## up and down peaks
            peakpath = {"%s_UPpeak"%comp:os.path.join(parameters['output_dir']+'/DESeq2', parameters['prefix']+comp+'_UPpeak.bed'),
                        "%s_DOWNpeak"%comp:os.path.join(parameters['output_dir']+'/DESeq2', parameters['prefix']+comp+'_DOWNpeak.bed'),
                        "%s_SHAREDpeak"%comp:os.path.join(parameters['output_dir']+'/DESeq2', parameters['prefix']+comp+'_SHAREDpeak.bed')}
            for p in peakpath:
                print('++++ %s'%p)
                peakbed = peakpath[p]
                linen = int(subprocess.getoutput("wc -l %s|cut -f1 -d ' '"%peakbed))
                if linen > 100:
                    if not os.path.exists(peakbed+'.motif'):
                        cmd = 'findMotifsGenome.pl {peak} {species} {output} -size 600 -mask'.format(
                            peak = peakbed,
                            species = parameters['species'],
                            output = peakbed+'.motif'
                            )
                        print(cmd)
                        os.system(cmd)
    if 'bedtools' in parameters['diffmethod']:
        design = pd.read_csv(parameters['diffdesignBedtools'], index_col = 0)
        comprisons = design.index.tolist()
        for comp in comprisons:
            ## up and down peaks
            peakpath = {"%s_UPpeak"%comp:os.path.join(parameters['output_dir']+'/bedtools', parameters['prefix']+comp+'_UPpeak.bed'),
                        "%s_DOWNpeak"%comp:os.path.join(parameters['output_dir']+'/bedtools', parameters['prefix']+comp+'_DOWNpeak.bed'),
                        "%s_SHAREDpeak"%comp:os.path.join(parameters['output_dir']+'/bedtools', parameters['prefix']+comp+'_SHAREDpeak.bed')}
            for p in peakpath:
                print('++++ %s'%p)
                peakbed = peakpath[p]
                linen = int(subprocess.getoutput("wc -l %s|cut -f1 -d ' '"%peakbed))
                if linen > 100:
                    if not os.path.exists(peakbed+'.motif'):
                        cmd = 'findMotifsGenome.pl {peak} {species} {output} -size 600 -mask -seqlogo'.format(
                            peak = peakbed,
                            species = parameters['species'],
                            output = peakbed+'.motif'
                            )
                        print(cmd)
                        os.system(cmd)
        
def main():
    parameters = get_parameter()
    print(parameters)
    parameters = _peakStat(parameters)
    parameters = _getUnionPeak(parameters)
    if parameters['bw_files'] is not None: 
        parameters = _bwHeatmap(parameters)
    parameters = _bamCount(parameters)
#    parameters = {'peak_files': ['Parental_K27_peaks.narrowPeak', 'TWF20_K27_peaks.narrowPeak', 'TWF30_K27_peaks.narrowPeak'], 'peakfile_label': ['Parental_K27', 'TWF20_K27', 'TWF30_K27'], 'bw_files': ['Parental_K27.rep1_treat_pileup.bw', 'Parental_K27.rep2_treat_pileup.bw', 'TWF20_K27.rep1_treat_pileup.bw', 'TWF20_K27.rep2_treat_pileup.bw', 'TWF30_K27.rep1_treat_pileup.bw', 'TWF30_K27.rep2_treat_pileup.bw'], 'bw_label': ['Parental_K27_1', 'Parental_K27_2', 'TWF20_K27_1', 'TWF20_K27_2', 'TWF30_K27_1', 'TWF30_K27_2'], 'bam_files': ['Parental_K27_1_unique.sorted.bam.dedup.bam', 'Parental_K27_2_unique.sorted.bam.dedup.bam', 'TWF20_K27_1_unique.sorted.bam.dedup.bam', 'TWF20_K27_2_unique.sorted.bam.dedup.bam', 'TWF30_K27_1_unique.sorted.bam.dedup.bam', 'TWF30_K27_2_unique.sorted.bam.dedup.bam'], 'bam_label': ['Parental_K27_1', 'Parental_K27_2', 'TWF20_K27_1', 'TWF20_K27_2', 'TWF30_K27_1', 'TWF30_K27_2'], 'design_matrix': 'design.csv', 'kmeans': '7', 'fold5': False, 'output_dir': './allPeak_new', 'prefix': 'K27_Par_TWF20_TWF30_', 'species': 'hg38', 'diffmethod': ['deseq2', 'bedtools'], 'diffdesignBedtools': 'design_bed.csv', 'unionPeak': './allPeak_new/K27_Par_TWF20_TWF30_union_dedup_peaks.bed', 'bwSignalOutput': './allPeak_new/K27_Par_TWF20_TWF30_union_dedup_peaks_3k.mat.gz', 'bamCountOutput': './allPeak_new/K27_Par_TWF20_TWF30_union_dedup_peaks_count_mat.txt', 'BamMapStat': './allPeak_new/K27_Par_TWF20_TWF30_mapping_stat.txt'}
    _diffPeakAnalysis(parameters)
    print(parameters)
    _motif_homer(parameters)
    ## record parameters
    out = open('%s/parameters.json'%parameters['output_dir'], 'w')
    json.dump(parameters, out)
    out.close()
    print(parameters)
    print('++++ Finished ++++')

main()



