import xlmhg,os
import numpy as np
from mne.stats import fdr_correction
from corelib import *
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

plt.rcParams.update(plt.rcParamsDefault)
rc={"axes.labelsize": 16, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "figure.titleweight":"bold", #"font.size":14,
    "figure.figsize":(5.5,4.2), "font.weight":"regular", "legend.fontsize":10,
    'axes.labelpad':8, 'figure.dpi':300}
plt.rcParams.update(**rc)

def _plot_(res_go, title, prefix, top_n):
    res_go = pd.DataFrame(res_go)
    res_go.columns = ['ID', 'term', 'N', 'p', 'padj']
    res_go = res_go.sort_values('padj')
    res_go = res_go[res_go['padj']<0.25]
    if res_go.shape[0] > 0:
        res_go = res_go.head(top_n)#[d['padj']<0.01]
        res_go['-log10padj'] = -np.log10(res_go['padj'])
        fig, ax = plt.subplots(figsize = (12, 5))
        sns.barplot(data = res_go, x = '-log10padj', y = 'term', color = 'darkorange')
        sns.despine()
        ax.set_title(title)
        ax.set_ylabel('')
        plt.tight_layout()
        fig.savefig(prefix+'_top_terms.png')
        plt.close()
        

### Change the input of function to a dict {gene: rp/rpt}
def geneIndices(gene_score_dict, assembly, rp=True, up=2000, down=10):
    go_terms = ["kegg", "cc", "mf", "bp"]
    genes_go_indices = {}
    for go_term in go_terms:
        Info("Indexing genes of %s..." % go_term)
        genes_go_indices[go_term] = {}
        go_symbol = readJson("%s_%s_symbol_dedup.json" % (assembly[:2], go_term))
        genes_go_indices[go_term]["maxL"] = len(set(go_symbol["genes"]) & set(gene_score_dict["genes"]))
        genes_go = list(set(go_symbol["genes"]) & set(gene_score_dict["scores"].keys()))
        genes_score_list = [gene_score_dict["scores"][i] for i in genes_go]
        genes_score_rank = calculateRank2(genes_score_list, rp)
        genes_go_indices[go_term]["N"] = max(genes_score_rank)
        genes_score_rank_dict = {}
        for igene in range(len(genes_go)):
            genes_score_rank_dict[genes_go[igene]] = genes_score_rank[igene]
        genes_go_indices[go_term]["indices"] = {}
        for term in list(go_symbol["ontology"].keys()):
            term_len = len(go_symbol["ontology"][term])
            if term_len >= up or term_len <= down:
                continue
            genes_go_term_index = [genes_score_rank_dict[i] for i in go_symbol["ontology"][term]]
            genes_go_indices[go_term]["indices"][term] = sorted(set(genes_go_term_index))
    return genes_go_indices

### Function for gene annotation
def annotIndices(genes_score_dict, assembly, rp, up, down, prefix):
    genes_go_indices = geneIndices(genes_score_dict, assembly, rp, up, down)
    gokegg_id_link = readJson("gokegg_id_link.json")
    go_terms = ["kegg", "cc", "mf", "bp"]
    for go_type in go_terms:
        Info("Enrichment analysis of %s..." % go_type)
        gene_dict = genes_go_indices[go_type]
        genes_indices = list(gene_dict["indices"].values())
        genes_indices_keys = list(gene_dict["indices"].keys())
        N = gene_dict["N"]
        L = min(2000, gene_dict["maxL"])
        def getmHG(genes_indices, N=N, X=1, L=L):
            # genes_indices = sorted(set(genes_indices))
            go_indices = np.uint16(genes_indices)
            res = xlmhg.get_xlmhg_test_result(N, go_indices, X=X, L=L)
            try:
                res_list = ["%s (%s,%s,%s,%s)"%(round(res.fold_enrichment, 3), res.N, res.cutoff, res.K, res.k), res.pval]
            except:
                res_list = ["0 (%s,%s,%s,%s)"%(res.N, res.cutoff, res.K, res.k), res.pval]
            # res_list = [res.N, res.cutoff, res.K, res.k, res.fold_enrichment, res.pval]
            return res_list
        res = list(map(getmHG, genes_indices))
        pvals = [x[-1] for x in res]
        fdr = list(fdr_correction(pvals)[1])
        res_go = [[genes_indices_keys[i],
                    gokegg_id_link[genes_indices_keys[i]].split(">")[1].split("<")[0]] +
                    res[i] + [ fdr[i] ] for i in range(len(res))]
        res_go.sort(key=lambda x: x[-1])
        _plot_(res_go=res_go, title=os.path.basename(prefix)+'_'+go_type, prefix = prefix+'_'+go_type, top_n=20)
        res_go_str = "\n".join(["\t".join(map(str, x)) for x in res_go])
        f = open("%s_%s.txt" % (prefix, go_type), "w")
        f.write('ID\tterm\tN\tp\tpadj\n')
        f.write(res_go_str)
        f.close()
