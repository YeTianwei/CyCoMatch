"""Pairwise stage: landmark filtering (LMD_Filter.m), MWP reconstruction
(MWP_reconstruct.m), fusion (add_fmap.m) and fmap -> pmap conversion
(fmap2pmap_group.m)."""
import numpy as np

from .knn import knnsearch


def find_k(D):
    """findK.m: number of significant singular values."""
    d = D / np.mean(D[:100])
    for s in range(len(d) - 10):
        if np.sum(np.abs(d[s:s + 11])) < 0.1:
            return s + 1
    return len(d)


def _svd_orth(H):
    U, s, Vt = np.linalg.svd(H)
    return U @ Vt, s


def _distortion(Di, Dj, T, area, matlab_compat, chunk=256):
    """Landmark metric distortion for every source vertex of the map T: i -> j."""
    n = len(T)
    R_max = Dj.max()
    a = area[:n, None]
    out = np.empty(n)
    for s in range(0, n, chunk):
        cols = np.arange(s, min(s + chunk, n))
        ds = Di[:n, cols].astype(np.float64)
        if matlab_compat:
            # The original LMD_Filter.m indexed the target distances with the *source*
            # vertex ids: Dist{j}(idx_source, T(k)).
            dt = Dj[:n][:, T[cols]].astype(np.float64)
        else:
            # d_j(T(m), T(k)), the metric distortion as described in the paper
            dt = Dj[np.ix_(T, T[cols])].astype(np.float64)
        dt[dt == 0] = R_max
        valid = ds != 0
        maxd = ds.max(axis=0)
        num = (np.abs(ds - dt) / maxd * a * valid).sum(axis=0)
        den = (a * valid).sum(axis=0)
        out[cols] = num / den
    return out


def lmd_filter(shapes, T0, threshold, matlab_compat=False):
    """Returns the refined maps, the spectral dimension per pair and C1."""
    n = len(shapes)
    T = [[None] * n for _ in range(n)]
    dim = np.zeros((n, n), dtype=int)
    C1 = [[None] * n for _ in range(n)]
    for i in range(n):
        Si = shapes[i]
        for j in range(n):
            Sj = shapes[j]
            nv = min(Si.nv, Sj.nv)
            Tij = np.array(T0[i][j], copy=True)
            dist = _distortion(Si.Gamma, Sj.Gamma, Tij[:nv], Sj.area_v, matlab_compat)

            L = np.flatnonzero(dist < threshold)
            NL = np.setdiff1d(np.arange(nv), L)
            NL_corr = np.setdiff1d(np.arange(nv), Tij[L])

            H = Si.evecs[L].T @ Sj.evecs[Tij[L]]
            C1[i][j], sv = _svd_orth(H)
            d = find_k(sv)
            dim[i, j] = d

            C_fmap, _ = _svd_orth(Si.evecs[L, :d].T @ Sj.evecs[Tij[L], :d])
            if len(NL) and len(NL_corr):
                nn = knnsearch(Sj.evecs[NL_corr, :d] @ C_fmap.T, Si.evecs[NL, :d])
                Tij[NL] = NL_corr[nn]
            T[i][j] = Tij
    return T, dim, C1


def _meyer_kernel(x, kind):
    l1, l2, l3 = 2 / 3, 4 / 3, 8 / 3
    v = lambda t: t ** 4 * (35 - 84 * t + 70 * t ** 2 - 20 * t ** 3)
    r = np.zeros_like(x)
    if kind == "sf":
        r[x < l1] = 1
        m = (x >= l1) & (x < l2)
        r[m] = np.cos(np.pi / 2 * v(np.abs(x[m]) / l1 - 1))
    else:
        m = (x >= l1) & (x < l2)
        r[m] = np.sin(np.pi / 2 * v(np.abs(x[m]) / l1 - 1))
        m = (x >= l2) & (x < l3)
        r[m] = np.cos(np.pi / 2 * v(np.abs(x[m]) / l2 - 1))
    return r


def meyer_filterbank(lmax, Nf, x):
    """gsp_design_meyer(lmax, Nf) evaluated at x; returns an (Nf, len(x)) array."""
    t = (4 / (3 * lmax)) * 2.0 ** np.arange(Nf - 2, -1, -1)
    rows = [_meyer_kernel(t[0] * x, "sf")]
    rows += [_meyer_kernel(t[k] * x, "wavelet") for k in range(Nf - 1)]
    return np.stack(rows)


def mwp_reconstruct(shapes, T0, Nf, num_iters, dim):
    n = len(shapes)
    Cout = [[None] * n for _ in range(n)]
    pinv_cache = {}
    for i in range(n):
        for j in range(n):
            if i == j:
                continue  # the diagonal is replaced by C1 in add_fmap
            Si, Sj = shapes[i], shapes[j]
            d = dim[i, j]
            ev_i, ev_j = Si.evals[:d], Sj.evals[:d]
            g = meyer_filterbank(Si.evals[d - 1], Nf, ev_i)
            ref_g = meyer_filterbank(Sj.evals[d - 1], Nf, ev_j)
            # sum_s diag(g_s) C diag(ref_g_s) == C .* (g' * ref_g)
            Wf = g.T @ ref_g
            if (i, d) not in pinv_cache:
                pinv_cache[(i, d)] = np.linalg.pinv(Si.evecs[:, :d])
            P = pinv_cache[(i, d)]

            T = T0[i][j]
            for _ in range(num_iters):
                C = (P @ Sj.evecs[T, :d]) * Wf
                T = knnsearch(Sj.evecs[:, :d] @ C.T, Si.evecs[:, :d])
            Cout[i][j] = C
    return Cout


def add_fmap(C1, C2, alpha):
    n = len(C1)
    C = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                C[i][j] = C1[i][j]
                continue
            k = C1[i][j].shape[0]
            r, c = C2[i][j].shape
            pad = np.zeros((k, k))
            pad[:r, :c] = C2[i][j]
            C[i][j] = alpha * C1[i][j] + (1 - alpha) * pad
    return C


def fmap2pmap_group(shapes, C):
    n = len(C)
    return [[knnsearch(shapes[j].evecs @ C[i][j].T, shapes[i].evecs) for j in range(n)]
            for i in range(n)]
