## below code is to run the ASW and NMI evaluation metrics
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import time
from datetime import timedelta
import scanpy as sc
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
print(sc.logging.print_versions())
import os
dirname = os.getcwd()
print(dirname)
from sklearn.metrics import silhouette_score
import random
from scipy import sparse
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import silhouette_score, normalized_mutual_info_score, silhouette_samples

def silhouette_coeff_ASW(adata, method_use='raw',save_dir='', task_name='', percent_extract=0.8):
    random.seed(0)
    asw_fscore = []
    asw_bn = []
    asw_bn_sub = []
    asw_ctn = [] 
    iters = []
    for i in range(20):
        iters.append('iteration_'+str(i+1))
        rand_cidx = np.random.choice(adata.obs_names, size=int(len(adata.obs_names) * percent_extract), replace=False)
        adata_ext = adata[rand_cidx,:]
        asw_batch = silhouette_score(adata_ext.X, adata_ext.obs['batch'])
        asw_celltype = silhouette_score(adata_ext.X, adata_ext.obs['cell_type'])
        min_val = -1
        max_val = 1
        asw_batch_norm = (asw_batch - min_val) / (max_val - min_val)
        asw_celltype_norm = (asw_celltype - min_val) / (max_val - min_val)
        
        fscoreASW = (2 * (1 - asw_batch_norm)*(asw_celltype_norm))/(1 - asw_batch_norm + asw_celltype_norm)
        asw_fscore.append(fscoreASW)
        asw_bn.append(asw_batch_norm)
        asw_bn_sub.append(1-asw_batch_norm)
        asw_ctn.append(asw_celltype_norm)

    df = pd.DataFrame({'asw_batch_norm':asw_bn, 'asw_batch_norm_sub': asw_bn_sub,
                       'asw_celltype_norm': asw_ctn, 'fscore':asw_fscore,
                       'method_use':np.repeat(method_use, len(asw_fscore))})
    df_mean=df.mean()
    df=df.append(df_mean,ignore_index=True)
    df['method_use'][20]='mean'
    df.to_csv(save_dir + task_name + "_"+method_use+'_ASW_metric.csv')
    print('Save output of pca in: ',save_dir)
    return df


def createAnnData(data_dir,myDatafn,sep_symbol):
    myData = pd.read_table(os.path.join(data_dir, myDatafn),header=0, index_col=0,sep=sep_symbol)
    myData.index=myData.index.astype(str)
    bex = ['batch','Batch','Batchlb','batchlb','BATCH']
    ib = np.isin(myData.keys(), bex)
    cex = ['celltype','CellType','cell_type','Cell_Type','ct']
    ict = np.isin(myData.keys(), cex)
    adata = sc.AnnData(myData.values[:,:-2])
    adata.obs_names = myData.index
    adata.obs['batch'] = myData.values[:, np.where(ib)[0][0]]  # factor function in R
    adata.obs['cell_type'] = myData.values[:, np.where(ict)[0][0]]
    print(adata)
    return adata

## below function is to run the ASW scores for the species specific cell types
def silhouette_coeff_ASW_species_specific_celltype(adata, celltype_names=[],method_use='raw',save_dir='', task_name='', percent_extract=0.8):
    random.seed(0)
    asw_fscore = []
    asw_bn = []
    asw_bn_sub = []
    asw_ctn = []
    iters = []
    for celltype_name in celltype_names:
        i=0
        iters.append('iteration_'+str(i+1))
        #rand_cidx = np.random.choice(adata.obs_names, size=int(len(adata.obs_names) * percent_extract), replace=False)
        #rand_cidx = adata.obs[adata.obs['cell_type']==celltype_name].index
        # Get the boolean mask for the desired cell type
        mask = adata.obs['cell_type'] == celltype_name

        # Get the positional indices where the mask is True
        rand_cidx = list(adata.obs.index.get_indexer(mask[mask].index))
        #import pdb;pdb.set_trace()
        #asw_batch = silhouette_samples(adata_ext.X, adata_ext.obs['batch'])
        asw_celltype = silhouette_samples(adata.X, adata.obs['cell_type'])
        min_val = -1
        max_val = 1
        #asw_batch_norm = (asw_batch - min_val) / (max_val - min_val)
        asw_celltype_norm = (asw_celltype - min_val) / (max_val - min_val)

        #fscoreASW = (2 * (1 - asw_batch_norm)*(asw_celltype_norm))/(1 - asw_batch_norm + asw_celltype_norm)
        #asw_fscore.append(fscoreASW)
        #asw_bn.append(asw_batch_norm)
        #asw_bn_sub.append(1-asw_batch_norm)
        asw_ctn.append(np.mean(asw_celltype_norm[rand_cidx]))
    
    #import pdb;pdb.set_trace()
    df = pd.DataFrame({'asw_celltype_norm': asw_ctn,
                       'method_use':np.repeat(method_use, len(asw_ctn))})
    df_mean=df.mean()
    df=df.append(df_mean,ignore_index=True)
    df.to_csv(save_dir + task_name + "_"+method_use+'_ASW_metric_species_sepcific_celltype.csv')
    print('Save output of pca in: ',save_dir)
    return df

########################### nmi 
def opt_louvain(adata, label_key, cluster_key, function=None, resolutions=None,
                inplace=True, plot=False, verbose=True, **kwargs):
    """
       This Louvain Clustering method was taken from scIB:
       Title: scIB
       Authors: Malte Luecken,
                Maren Buettner,
                Daniel Strobl,
                Michaela Mueller
       Date: 4th October 2020
       Code version: 0.2.0
       Availability: https://github.com/theislab/scib/blob/master/scIB/clustering.py
    params:
        label_key: name of column in adata.obs containing biological labels to be
            optimised against
        cluster_key: name of column to be added to adata.obs during clustering.
            Will be overwritten if exists and `force=True`
        function: function that computes the cost to be optimised over. Must take as
            arguments (adata, group1, group2, **kwargs) and returns a number for maximising
        resolutions: list if resolutions to be optimised over. If `resolutions=None`,
            default resolutions of 20 values ranging between 0.1 and 2 will be used
    returns:
        res_max: resolution of maximum score
        score_max: maximum score
        score_all: `pd.DataFrame` containing all scores at resolutions. Can be used to plot the score profile.
        clustering: only if `inplace=False`, return cluster assignment as `pd.Series`
        plot: if `plot=True` plot the score profile over resolution
    """
    adata = remove_sparsity(adata)

    if resolutions is None:
        n = 20
        resolutions = [2 * x / n for x in range(1, n + 1)]

    score_max = 0
    res_max = resolutions[0]
    clustering = None
    score_all = []

    # maren's edit - recompute neighbors if not existing
    try:
        adata.uns['neighbors']
    except KeyError:
        if verbose:
            print('computing neigbours for opt_cluster')
        sc.pp.neighbors(adata)

    for res in resolutions:
        sc.tl.louvain(adata, resolution=res, key_added=cluster_key)
        score = function(adata, label_key, cluster_key, **kwargs)
        score_all.append(score)
        if score_max < score:
            score_max = score
            res_max = res
            clustering = adata.obs[cluster_key]
        del adata.obs[cluster_key]

    if verbose:
        print(f'optimised clustering against {label_key}')
        print(f'optimal cluster resolution: {res_max}')
        print(f'optimal score: {score_max}')

    score_all = pd.DataFrame(zip(resolutions, score_all), columns=('resolution', 'score'))
    if plot:
        # score vs. resolution profile
        sns.lineplot(data=score_all, x='resolution', y='score').set_title('Optimal cluster resolution profile')
        plt.show()

    if inplace:
        adata.obs[cluster_key] = clustering
        return res_max, score_max, score_all
    else:
        return res_max, score_max, score_all, clustering

def nmi_helper(adata, group1, group2, method="arithmetic"):
    """
       This NMI function was taken from scIB:
       Title: scIB
       Authors: Malte Luecken,
                Maren Buettner,
                Daniel Strobl,
                Michaela Mueller
       Date: 4th October 2020
       Code version: 0.2.0
       Availability: https://github.com/theislab/scib/blob/master/scIB/metrics.py
       Normalized mutual information NMI based on 2 different cluster assignments `group1` and `group2`
       params:
        adata: Anndata object
        group1: column name of `adata.obs` or group assignment
        group2: column name of `adata.obs` or group assignment
        method: NMI implementation
            'max': scikit method with `average_method='max'`
            'min': scikit method with `average_method='min'`
            'geometric': scikit method with `average_method='geometric'`
            'arithmetic': scikit method with `average_method='arithmetic'`
       return:
        normalized mutual information (NMI)
    """
    adata = remove_sparsity(adata)

    if isinstance(group1, str):
        group1 = adata.obs[group1].tolist()
    elif isinstance(group1, pd.Series):
        group1 = group1.tolist()

    labels = adata.obs[group2].values
    labels_encoded = LabelEncoder().fit_transform(labels)
    group2 = labels_encoded

    if len(group1) != len(group2):
        raise ValueError(f'different lengths in group1 ({len(group1)}) and group2 ({len(group2)})')

    # choose method
    if method in ['max', 'min', 'geometric', 'arithmetic']:
        nmi_value = normalized_mutual_info_score(group1, group2, average_method=method)
    else:
        raise ValueError(f"Method {method} not valid")

    return nmi_value

def remove_sparsity(adata):
    """
        If ``adata.X`` is a sparse matrix, this will convert it in to normal matrix.
        Parameters
        ----------
        adata: :class:`~anndata.AnnData`
            Annotated data matrix.
        Returns
        -------
        adata: :class:`~anndata.AnnData`
            Annotated dataset.
    """
    if sparse.issparse(adata.X):
        new_adata = sc.AnnData(X=adata.X.A, obs=adata.obs.copy(deep=True), var=adata.var.copy(deep=True))
        return new_adata

    return adata


def nmi(adata, label_key, verbose=False, nmi_method='arithmetic'):
    cluster_key = 'cluster'
    opt_louvain(adata, label_key=label_key, cluster_key=cluster_key, function=nmi_helper,
                plot=False, verbose=verbose, inplace=True)
    #import pdb; pdb.set_trace()
    print('NMI...')
    nmi_score = nmi_helper(adata, group1=cluster_key, group2=label_key, method=nmi_method)

    return nmi_score
