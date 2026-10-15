"""Nearest-neighbour search (replacement for MATLAB knnsearch, CPU or torch GPU)."""
import numpy as np

try:
    import torch

    _DEVICE = "cuda" if torch.cuda.is_available() else None
except ImportError:  # pragma: no cover
    torch, _DEVICE = None, None


def knnsearch(X, Y, chunk=2048):
    """For every row of Y return the index of its nearest row in X."""
    if Y.shape[0] == 0:
        return np.zeros(0, dtype=np.int64)
    if _DEVICE is not None:
        Xt = torch.as_tensor(X, device=_DEVICE)
        Yt = torch.as_tensor(Y, device=_DEVICE)
        out = [torch.cdist(Yt[s:s + chunk], Xt).argmin(dim=1)
               for s in range(0, Yt.shape[0], chunk)]
        return torch.cat(out).cpu().numpy()

    X = np.ascontiguousarray(X)
    xx = np.einsum("ij,ij->i", X, X)
    out = np.empty(Y.shape[0], dtype=np.int64)
    for s in range(0, Y.shape[0], chunk):
        Yc = Y[s:s + chunk]
        # ||x||^2 - 2 x.y  (||y||^2 is constant per query)
        out[s:s + chunk] = np.argmin(xx[None, :] - 2.0 * (Yc @ X.T), axis=1)
    return out
