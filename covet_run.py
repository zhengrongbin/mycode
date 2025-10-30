import os,sys
import numpy as np
import pandas as pd
import collections
import pickle as pk
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import colors as mcolors
import seaborn as sns
# from adjustText import adjust_text
from scipy.stats import pearsonr, spearmanr, chisquare, ttest_ind, ranksums, wilcoxon, fisher_exact, mannwhitneyu
from scipy import stats
import scanpy as sc
# import squidpy as sq
from scipy.cluster import hierarchy as sch
from sklearn.cluster import KMeans
import scenvi
import sklearn.neighbors
import sklearn
import scipy.sparse
import umap.umap_ as umap
from fa2 import ForceAtlas2

os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"]="0" # Change to -1 if you want to use CPU!


def force_directed_layout(affinity_matrix, cell_names=None, verbose=True, iterations=500, device='cpu'):
    """" Function to compute force directed layout from the affinity_matrix
    :param affinity_matrix: Sparse matrix representing affinities between cells
    :param cell_names: pandas Series object with cell names
    :param verbose: Verbosity for force directed layout computation
    :param iterations: Number of iterations used by ForceAtlas
    :return: Pandas data frame representing the force directed layout
    """

    init_coords = np.random.random((affinity_matrix.shape[0], 2))

    if device == 'cpu':
        forceatlas2 = ForceAtlas2(
            # Behavior alternatives
            outboundAttractionDistribution=False,
            linLogMode=False,
            adjustSizes=False,
            edgeWeightInfluence=1.0,
            # Performance
            jitterTolerance=1.0,
            barnesHutOptimize=True,
            barnesHutTheta=1.2,
            multiThreaded=False,
            # Tuning
            scalingRatio=2.0,
            strongGravityMode=False,
            gravity=1.0,
            # Log
            verbose=verbose)

        positions = forceatlas2.forceatlas2(
            affinity_matrix, pos=init_coords, iterations=iterations)
        positions = np.array(positions)


    positions = pd.DataFrame(positions,
                             index=np.arange(affinity_matrix.shape[0]), columns=['x', 'y'])
    return positions

def run_diffusion_maps(data_df, n_components=10, knn=30, alpha=0):
    """Run Diffusion maps using the adaptive anisotropic kernel
    :param data_df: PCA projections of the data or adjacency matrix
    :param n_components: Number of diffusion components
    :param knn: Number of nearest neighbors for graph construction
    :param alpha: Normalization parameter for the diffusion operator
    :return: Diffusion components, corresponding eigen values and the diffusion operator
    """

    # Determine the kernel
    N = data_df.shape[0]

    if(type(data_df).__module__ == np.__name__):
        data_df = pd.DataFrame(data_df)

    if not scipy.sparse.issparse(data_df):
        print("Determing nearest neighbor graph...")
        temp = sc.AnnData(data_df.values)
        sc.pp.neighbors(temp, n_pcs=0, n_neighbors=knn)
        kNN = temp.obsp['distances']

        # Adaptive k
        adaptive_k = int(np.floor(knn / 3))
        adaptive_std = np.zeros(N)

        for i in np.arange(len(adaptive_std)):
            adaptive_std[i] = np.sort(kNN.data[kNN.indptr[i] : kNN.indptr[i + 1]])[
                adaptive_k - 1
            ]

        # Kernel
        x, y, dists = scipy.sparse.find(kNN)

        # X, y specific stds
        dists = dists / adaptive_std[x]
        W = scipy.sparse.csr_matrix((np.exp(-dists), (x, y)), shape=[N, N])

        # Diffusion components
        kernel = W + W.T
    else:
        kernel = data_df

    # Markov
    D = np.ravel(kernel.sum(axis=1))

    if alpha > 0:
        # L_alpha
        D[D != 0] = D[D != 0] ** (-alpha)
        mat = scipy.sparse.csr_matrix((D, (range(N), range(N))), shape=[N, N])
        kernel = mat.dot(kernel).dot(mat)
        D = np.ravel(kernel.sum(axis=1))

    D[D != 0] = 1 / D[D != 0]
    T = scipy.sparse.csr_matrix((D, (range(N), range(N))), shape=[N, N]).dot(kernel)
    # Eigen value dcomposition
    D, V = scipy.sparse.linalg.eigs(T, n_components, tol=1e-4, maxiter=1000)
    D = np.real(D)
    V = np.real(V)
    inds = np.argsort(D)[::-1]
    D = D[inds]
    V = V[:, inds]

    # Normalize
    for i in range(V.shape[1]):
        V[:, i] = V[:, i] / np.linalg.norm(V[:, i])

    # Create are results dictionary
    res = {"T": T, "EigenVectors": V, "EigenValues": D}
    res["EigenVectors"] = pd.DataFrame(res["EigenVectors"])
    if not scipy.sparse.issparse(data_df):
        res["EigenVectors"].index = data_df.index
    res["EigenValues"] = pd.Series(res["EigenValues"])
    res["kernel"] = kernel

    return res


def FDL(data, k = 30):


    nbrs = sklearn.neighbors.NearestNeighbors(n_neighbors=int(k), metric='euclidean',
                               n_jobs=5).fit(data)
    kNN = nbrs.kneighbors_graph(data, mode='distance')
    # Adaptive k

    adaptive_k = int(np.floor(k / 3))
    nbrs = sklearn.neighbors.NearestNeighbors(n_neighbors=int(adaptive_k),
                           metric='euclidean', n_jobs=5).fit(data)
    adaptive_std = nbrs.kneighbors_graph(data, mode='distance').max(axis=1)
    adaptive_std = np.ravel(adaptive_std.todense())
    # Kernel
    x, y, dists = scipy.sparse.find(kNN)
    # X, y specific stds
    dists = dists / adaptive_std[x]
    N = data.shape[0]
    W = scipy.sparse.csr_matrix((np.exp(-dists), (x, y)), shape=[N, N])
    # Diffusion components
    kernel = W + W.T
    layout = force_directed_layout(kernel)
    return(layout)

def flatten(arr):
    return(np.reshape(arr, [arr.shape[0], -1]))



st_data = sc.read_h5ad('./label_transfer/adata_vgz_combined_process_volume_norm_label_transfer_mouse_tabulaMuscle_updated_celllabel_updated.h5ad')
st_data = st_data[st_data.obs['SampleName'].isin(['BAT_Cold7_1', 'BAT_Cold7_2', 'BAT_Cold7_3'])]

res  = scenvi.compute_covet(st_data)
st_data.obsm['COVET'], st_data.obsm['COVET_SQRT'], st_data.uns['CovGenes'] = res
out = open('COVET_cold7.pk', 'wb')
pk.dump(res, out)
out.close()

data = flatten(st_data.obsm['COVET_SQRT'])
leiden_ann = sc.AnnData(X = data)
sc.pp.neighbors(leiden_ann, n_neighbors = 30)
sc.tl.leiden(leiden_ann, resolution=0.2)
out = open('COVET_leiden_cluster.pk', 'wb')
pk.dump(leiden_ann.obs, out)
out.close()


FDL_COVET = np.asarray(FDL(flatten(st_data.obsm['COVET_SQRT']), k = 10))
FDL_COVET_df = pd.DataFrame(FDL_COVET, index = st_data.obs_names, columns = ['COVET_FDL1', 'COVET_FDL2'])
FDL_COVET_df['SampleName'] = st_data.obs['SampleName']
FDL_COVET_df['cell_label'] = st_data.obs['cell_label_updated']
FDL_COVET_df['X'] = st_data.obsm['spatial'][:,0]
FDL_COVET_df['Y'] = st_data.obsm['spatial'][:,1]
FDL_COVET_df['leiden'] = leiden_ann.obs['leiden'].tolist()
out = open('FDL_COVET_df_cold7.pk', 'wb')
pk.dump(FDL_COVET_df, out)
out.close()

fit = umap.UMAP(
    n_neighbors = 10,
    min_dist = 0.8,
    n_components = 2,
)

UMAP_COVET = fit.fit_transform(flatten(st_data.obsm['COVET_SQRT']))
UMAP_COVET_df = pd.DataFrame(UMAP_COVET, index = st_data.obs_names, columns = ['COVET_UMAP1', 'COVET_UMAP2'])
UMAP_COVET_df['cell_label'] = st_data.obs['cell_label_updated']
UMAP_COVET_df['X'] = st_data.obsm['spatial'][:,0]
UMAP_COVET_df['Y'] = st_data.obsm['spatial'][:,1]
UMAP_COVET_df['leiden'] = leiden_ann.obs['leiden'].tolist()

out = open('UMAP_COVET_df_cold7.pk', 'wb')
pk.dump(UMAP_COVET_df, out)
out.close()

