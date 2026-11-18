"""LOCAL TESTS ONLY: a tiny stand-in for anndata, used when the real package is not installed.
The Docker image has the real anndata/scanpy (checked at build time by merge_sections.py --self-test)."""
import numpy as np
import pandas as pd


class AnnData:
    def __init__(self, X, var=None):
        self.X = X
        self.obs = pd.DataFrame(index=[str(i) for i in range(X.shape[0])])
        self.var = var if var is not None else pd.DataFrame(index=[str(i) for i in range(X.shape[1])])
        self.obsm = {}

    n_obs = property(lambda self: self.X.shape[0])
    n_vars = property(lambda self: self.X.shape[1])

    def __getitem__(self, key):
        rows, columns = key
        columns = np.asarray(columns)
        return AnnData(self.X[:, columns], self.var[columns].copy())

    def copy(self):
        return self
