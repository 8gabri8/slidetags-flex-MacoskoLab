"""LOCAL TESTS ONLY: a tiny stand-in for scanpy (same function names, simple sklearn maths), used
when the real package is not installed. The Docker image has the real scanpy."""
import types

import numpy as np
import scipy.sparse
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA


def _dense(x):
    return x.toarray() if scipy.sparse.issparse(x) else np.asarray(x)


def normalize_total(data, target_sum):
    x = _dense(data.X)
    data.X = x / np.maximum(x.sum(1, keepdims=True), 1) * target_sum


def log1p(data):
    data.X = np.log1p(_dense(data.X))


def highly_variable_genes(data, n_top_genes, flavor=None):
    variance = _dense(data.X).var(0)
    data.var["highly_variable"] = variance >= np.sort(variance)[::-1][n_top_genes - 1]


def pca(data, n_comps, random_state=0):
    data.obsm["X_pca"] = PCA(n_comps, random_state=random_state).fit_transform(_dense(data.X))


def neighbors(data, n_neighbors, random_state=0):
    pass


def umap(data, random_state=0):
    data.obsm["X_umap"] = data.obsm["X_pca"][:, :2]


def leiden(data, resolution=1.0, random_state=0):
    k = min(8, max(2, data.n_obs // 50))
    data.obs["leiden"] = KMeans(k, n_init=4, random_state=random_state).fit_predict(data.obsm["X_pca"][:, :10]).astype(str)


pp = types.SimpleNamespace(normalize_total=normalize_total, log1p=log1p, highly_variable_genes=highly_variable_genes,
                           pca=pca, neighbors=neighbors)
tl = types.SimpleNamespace(umap=umap, leiden=leiden)
