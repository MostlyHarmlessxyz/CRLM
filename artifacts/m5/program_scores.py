"""Positive-expression rank scores on a supplied common measurement axis."""

import numpy as np
import pandas as pd
from scipy import sparse


def rank_scores(matrix, symbols, programs, *, max_rank=1500):
    """Score fixed, already measurement-eligible symbol sets on cells × genes.

    Positive values use descending midranks, including all measured genes
    when computing ranks. Measured zeros have zero weight. The caller supplies
    exact measured program members; absence from an axis is not an expression
    zero. Returns one column per program in input order and one row per cell.
    """
    symbols = pd.Index(symbols)
    if not symbols.is_unique or len(symbols) != matrix.shape[1] or max_rank < 1:
        raise ValueError("Require a unique gene axis matching the matrix and a positive rank cap")
    names = list(programs)
    genes_in_sets, sets, weights = [], [], []
    for number, name in enumerate(names):
        members = list(programs[name])
        indices = symbols.get_indexer(members)
        if not members or len(set(members)) != len(members) or (indices < 0).any():
            raise ValueError(f"{name}: supply unique, nonempty measured program members")
        genes_in_sets.extend(indices)
        sets.extend([number] * len(indices))
        weights.extend([1 / len(indices)] * len(indices))
    design = sparse.csr_matrix((weights, (genes_in_sets, sets)), shape=(len(symbols), len(names)))
    relevant = np.diff(design.indptr) > 0
    scores = np.zeros((matrix.shape[0], len(names)), dtype=np.float64)
    for start in range(0, matrix.shape[0], 512):
        chunk = sparse.csr_matrix(matrix[start:start + 512], copy=True)
        chunk.sum_duplicates()
        chunk.eliminate_zeros()
        if not np.isfinite(chunk.data).all() or (chunk.data < 0).any():
            raise ValueError("Rank inputs must have finite nonnegative expression")
        ranks = np.zeros(chunk.nnz, dtype=np.float64)
        for a, b in zip(chunk.indptr[:-1], chunk.indptr[1:]):
            values = chunk.data[a:b]
            keep = relevant[chunk.indices[a:b]]
            if not keep.any():
                continue
            ordered = np.sort(values)
            target = values[keep]
            left = np.searchsorted(ordered, target, side="left")
            right = np.searchsorted(ordered, target, side="right")
            descending_midrank = len(ordered) + 0.5 - 0.5 * (left + right)
            ranks[a:b][keep] = np.maximum(0, 1 - (descending_midrank - 1) / max_rank)
        chunk.data = ranks
        chunk.eliminate_zeros()
        scores[start:start + chunk.shape[0]] = (chunk @ design).toarray()
    return pd.DataFrame(scores, columns=names)
