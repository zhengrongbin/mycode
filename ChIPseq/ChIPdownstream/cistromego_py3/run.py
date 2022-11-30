from rpscore import calcRPscore
from enrichment import annotIndices
from opt_validator import optValidate
import json
import pickle as pk
import pandas as pd
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
from scipy.stats import spearmanr, pearsonr
import seaborn as sns


def _plot_AUC_(x=[], y=[], label=[], score = [], test = '', pdf = False):
    fig, ax = plt.subplots()
    if test == 'AUROC':
        plt.plot([0, 1], [0, 1], color="black", lw=.5, linestyle="--")
    for i in range(len(x)):
        plt.plot(x[i], y[i], 
            lw=1,
            label="%s (AUC=%.2f)" % (label[i], score[i]))
    ax.set_xlabel('FPR') if test == 'AUROC' else ax.set_xlabel('Recall')
    ax.set_ylabel('TPR') if test == 'AUROC' else ax.set_xlabel('Precision')
    ax.legend()
    plt.tight_layout()
    pdf.savefig(fig) if pdf is not False else None
    plt.close()


def runCistromeGO(args):
    args = optValidate(args)
    res_rp_table, symbol_rp_dict, decay = calcRPscore(args.bed, args.peakn, args.assembly, args.decay)
    res_rp_table = pd.DataFrame(res_rp_table['data'])
    res_rp_table.columns = ['symbol', 'cordinate', 'peakn', 'RP', 'normRP', 'rank_by_normRP']
    res_rp_table['decay'] = decay
    res_rp_table = res_rp_table.sort_values('rank_by_normRP')
    # json.dump(symbol_rp_dict, open("%s_rp.json"%args.prefix, "w"))
    if not args.expr:
        res = annotIndices(symbol_rp_dict, args.assembly, True, args.max_gene_number, args.min_gene_number, args.prefix)
    else:
        from expr_combine import calcRPT
        symbol_rpt_dict, df, res_rate_dict = calcRPT(args.expr, args.exprinfo, symbol_rp_dict, args.assembly, args.dego, args.logfc_cut, args.fdr_cut)
        res_rp_table = pd.merge(res_rp_table, df[['Log2FC', 'FDR', 'rankproduct']], left_on = 'symbol', right_index = True)
        # print(res_rp_table.head())
        ## plot AUC
        pdf = PdfPages("%s_auc.pdf"%args.prefix)
        _plot_AUC_(x=[res_rate_dict['up']['FPR'], res_rate_dict['down']['FPR']],
                   y=[res_rate_dict['up']['TPR'], res_rate_dict['down']['TPR']], 
                   label=['Up-regulation', 'Down-regulation'],
                    score = [res_rate_dict['up']['auroc'], res_rate_dict['down']['auroc']],
                    test = 'AUROC',
                    pdf = pdf)
        _plot_AUC_(x=[res_rate_dict['up']['recall'], res_rate_dict['down']['recall']],
                   y=[res_rate_dict['up']['precision'], res_rate_dict['down']['precision']], 
                   label=['Up-regulation', 'Down-regulation'],
                    score = [res_rate_dict['up']['auprc'], res_rate_dict['down']['auprc']],
                    test = 'AUPRC',
                    pdf = pdf)
        pdf.close()
        ## enrichment
        res = annotIndices(symbol_rpt_dict, args.assembly, False, args.max_gene_number, args.min_gene_number, args.prefix)

    if 'rankproduct' in res_rp_table.columns.tolist():
        res_rp_table['rank_by_rankproduct'] = res_rp_table['rankproduct'].rank(ascending = True)
        res_rp_table = res_rp_table.sort_values('rank_by_rankproduct')
        
        ## plot a scatter plot to show correlation
        plot_df = res_rp_table[res_rp_table['FDR']<0.25]
        fig, ax = plt.subplots(figsize = (4.5, 4.2))
        # ax.scatter(np.sqrt(plot_df['RP']), plot_df['Log2FC'], s = 10, alpha = .6, edgecolor = 'none') 
        # ax.set_title('Pearson R=%.4f, Spearmanr R=%.4f'%(pearsonr(np.sqrt(plot_df['RP']), plot_df['Log2FC'])[0],
        #     spearmanr(np.sqrt(plot_df['RP']), plot_df['Log2FC'])[0]))
        ax.scatter(plot_df['normRP'], plot_df['Log2FC'], s = 18, alpha = .4, edgecolor = 'none') 
        # sns.regplot(data = plot_df, x = 'normRP', y = 'Log2FC', 
        #             scatter_kws={'s': 20, 'alpha': .4, "edgecolor": 'none'},
        #             line_kws = {'color': 'grey'}, ci = None) 
        # ax.set_title('Pearson R=%.4f, Spearmanr R=%.4f'%(pearsonr(plot_df['normRP'], plot_df['Log2FC'])[0],
        #     spearmanr(plot_df['normRP'], plot_df['Log2FC'])[0]))
        ax.set_title('Pearson R=%.4f'%(pearsonr(plot_df['normRP'], plot_df['Log2FC'])[0]))
        plt.axhline(y = 0, color = 'black', linestyle = '--')
        ax.set_xlabel('Regulatory Potential')
        ax.set_ylabel('Log2FC')
        sns.despine()
        plt.tight_layout()
        fig.savefig("%s_sig_corr.pdf"%args.prefix)
        plt.close()

        fig, ax = plt.subplots(figsize = (4.5, 4.2))
        ax.scatter(res_rp_table['normRP'], res_rp_table['Log2FC'], s = 18, alpha = .4, edgecolor = 'none') 
        # sns.regplot(data = res_rp_table, x = 'normRP', y = 'Log2FC', 
        #             scatter_kws={'s': 20, 'alpha': .4, "edgecolor": 'none'},
        #             line_kws = {'color': 'grey'}, ci = None) 
        # ax.set_title('Pearson R=%.4f, Spearmanr R=%.4f'%(pearsonr(res_rp_table['normRP'], res_rp_table['Log2FC'])[0],
        #     spearmanr(np.sqrt(res_rp_table['RP']), res_rp_table['Log2FC'])[0]))
        ax.set_title('Pearson R=%.4f'%(pearsonr(res_rp_table['normRP'], res_rp_table['Log2FC'])[0]))
        plt.axhline(y = 0, color = 'black', linestyle = '--')
        ax.set_xlabel('Regulatory Potential')
        ax.set_ylabel('Log2FC')
        sns.despine()
        plt.tight_layout()
        fig.savefig("%s_all_corr.pdf"%args.prefix)
        plt.close()

    res_rp_table.to_csv("%s_rp.txt"%args.prefix, index = None, sep = '\t')
    return res




