## below code is to run the LISI metric for the knn output methods (samap and bbknn)
import pandas as pd
import numpy as np
import scanpy as sc
import os
import sys
import scipy.sparse
from scipy.io import mmwrite
import itertools
import logging
import multiprocessing as mp
import pathlib
import subprocess
import tempfile
import warnings

def Hbeta(D_row, beta):
    """
    Helper function for simpson index computation
    """
    P = np.exp(-D_row * beta)
    sumP = np.nansum(P)
    if sumP == 0:
        H = 0
        P = np.zeros(len(D_row))
    else:
        H = np.log(sumP) + beta * np.nansum(D_row * P) / sumP
        P /= sumP
    return H, P
    
def convert_to_one_hot(vector, num_classes=None):
    """
    Converts an input 1-D vector of integers into an output 2-D array of one-hot vectors,
    where an i'th input value of j will set a '1' in the i'th row, j'th column of the
    output array.
    Example:
    .. code-block:: python
        v = np.array((1, 0, 4))
        one_hot_v = convertToOneHot(v)
        print(one_hot_v)
    .. code-block::
        [[0 1 0 0 0]
         [1 0 0 0 0]
         [0 0 0 0 1]]
    """

    # assert isinstance(vector, np.ndarray)
    # assert len(vector) > 0

    if num_classes is None:
        num_classes = np.max(vector) + 1
    # else:
    #    assert num_classes > 0
    #    assert num_classes >= np.max(vector)

    result = np.zeros(shape=(len(vector), num_classes))
    result[np.arange(len(vector)), vector] = 1
    return result.astype(int)

def compute_simpson_index(
    D=None, knn_idx=None, batch_labels=None, n_batches=None, perplexity=15, tol=1e-5
):
    """
    Simpson index of batch labels subset by group.
    :param D: distance matrix ``n_cells x n_nearest_neighbors``
    :param knn_idx: index of ``n_nearest_neighbors`` of each cell
    :param batch_labels: a vector of length n_cells with batch info
    :param n_batches: number of unique batch labels
    :param perplexity: effective neighborhood size
    :param tol: a tolerance for testing effective neighborhood size
    :returns: the simpson index for the neighborhood of each cell
    """
    n = D.shape[0]
    P = np.zeros(D.shape[1])
    simpson = np.zeros(n)
    logU = np.log(perplexity)

    # loop over all cells
    for i in np.arange(0, n, 1):
        beta = 1
        # negative infinity
        betamin = -np.inf
        # positive infinity
        betamax = np.inf
        # get active row of D
        D_act = D[i, :]
        H, P = Hbeta(D_act, beta)
        Hdiff = H - logU
        tries = 0
        # first get neighbor probabilities
        while np.logical_and(np.abs(Hdiff) > tol, tries < 50):
            if Hdiff > 0:
                betamin = beta
                if betamax == np.inf:
                    beta *= 2
                else:
                    beta = (beta + betamax) / 2
            else:
                betamax = beta
                if betamin == -np.inf:
                    beta /= 2
                else:
                    beta = (beta + betamin) / 2

            H, P = Hbeta(D_act, beta)
            Hdiff = H - logU
            tries += 1

        if H == 0:
            simpson[i] = -1
            continue

            # then compute Simpson's Index
        non_nan_knn = knn_idx[i][np.invert(np.isnan(knn_idx[i]))].astype("int")
        batch = batch_labels[non_nan_knn]
        # convertToOneHot omits all nan entries.
        # Therefore, we run into errors in np.matmul.
        if len(batch) == len(P):
            B = convert_to_one_hot(batch, n_batches)
            sumP = np.matmul(P, B)  # sum P per batch
            simpson[i] = np.dot(sumP, sumP)  # sum squares
        else:  # assign worst possible score
            simpson[i] = 1

    return simpson

def evaluate_LISI_knn(adata,batch_key,label_key,perplexity=40,verbose=False):
    # get the knn index matrix
    dist_mat = scipy.sparse.find(adata.obsp["connectivities"])
    n_nn = adata.uns["neighbors"]["params"]["n_neighbors"] - 1

    # initialise index and fill it with NaN values
    nn_index = np.empty(shape=(adata.obsp["connectivities"].shape[0], n_nn))
    nn_index[:] = np.nan
    nn_dists = np.empty(shape=(adata.obsp["connectivities"].shape[0], n_nn))
    nn_dists[:] = np.NaN
    index_out = []
    for cell_id in np.arange(np.min(dist_mat[0]), np.max(dist_mat[0]) + 1):
        get_idx = dist_mat[0] == cell_id
        num_idx = get_idx.sum()
        # in case that get_idx contains more than n_nn neighbours, cut away the outlying ones
        fin_idx = np.min([num_idx, n_nn])
        nn_index[cell_id, :fin_idx] = dist_mat[1][get_idx][
            np.argsort(dist_mat[2][get_idx])
        ][:fin_idx]
        nn_dists[cell_id, :fin_idx] = np.sort(dist_mat[2][get_idx])[:fin_idx]
        if num_idx < n_nn:
            index_out.append(cell_id)

    out_cells = len(index_out)

    if out_cells > 0:
        if verbose:
            print(f"{out_cells} had less than {n_nn} neighbors.")

    if perplexity is None:
        # use LISI default
        perplexity = np.floor(nn_index.shape[1] / 3)

    # run LISI in python
    if verbose:
        print("importing knn-graph")

    batch = adata.obs[batch_key].cat.codes.values
    n_batches = len(np.unique(adata.obs[batch_key]))
    label = adata.obs[label_key].cat.codes.values
    n_labels = len(np.unique(adata.obs[label_key]))

    if verbose:
        print("running LISI")

    simpson_estimate_batch = compute_simpson_index(
        D=nn_dists,
        knn_idx=nn_index,
        batch_labels=batch,
        n_batches=n_batches,
        perplexity=perplexity,
    )
    simpson_estimate_label = compute_simpson_index(
        D=nn_dists,
        knn_idx=nn_index,
        batch_labels=label,
        n_batches=n_labels,
        perplexity=perplexity,
    )
    simpson_est_batch = 1 / simpson_estimate_batch
    simpson_est_label = 1 / simpson_estimate_label
    # extract results
    d = {batch_key: simpson_est_batch, label_key: simpson_est_label}
    lisi_estimate = pd.DataFrame(data=d, index=np.arange(0, len(simpson_est_label)))

    return lisi_estimate

