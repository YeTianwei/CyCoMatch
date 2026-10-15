"""Cycle-consistency basis generation (port of external/cycle_basis_generator.m,
cycle_weight_opt_admm.m and cycle_sampling.m, from NeurIPS19_Cycle)."""
import numpy as np


def adj_mat_2_edge(A):
    """Undirected edges (i < j) in MATLAB `find` (column-major) order."""
    A = np.asarray(A)
    cols, rows = np.nonzero(A.T)
    keep = rows < cols
    edges = np.stack([rows[keep], cols[keep]], axis=0)
    edge_mat = -np.ones(A.shape, dtype=int)
    for e in range(edges.shape[1]):
        edge_mat[edges[0, e], edges[1, e]] = e
        edge_mat[edges[1, e], edges[0, e]] = e
    return edge_mat, edges


def _bfs_tree(nbrs, root):
    n = len(nbrs)
    parent = np.full(n, -2)  # -2: unvisited, -1: root
    parent[root] = -1
    visited = [root]
    start, end, depth = 0, 1, 0
    while start < end:
        depth += 1
        for vid in visited[start:end]:
            for nid in nbrs[vid]:
                if parent[nid] == -2:
                    parent[nid] = vid
                    visited.append(nid)
        start, end = end, len(visited)
    return parent, depth


def _path_to_root(parent, v):
    path = [v]
    while parent[v] != -1:
        v = parent[v]
        path.append(v)
    return path


def _cycles_from_tree(edges, parent):
    cycles = []
    for v1, v2 in edges.T:
        if parent[v1] != v2 and parent[v2] != v1:
            p1 = _path_to_root(parent, v1)
            p2 = _path_to_root(parent, v2)
            cycles.append(p1 + p2[::-1][1:])
    return cycles


def _jacobi_matrix(cycles, edge_mat, num_e):
    J = np.zeros((num_e, len(cycles)))
    for c, cyc in enumerate(cycles):
        prev = [cyc[-1]] + cyc[:-1]
        for a, b in zip(cyc, prev):
            J[edge_mat[a, b], c] += 1 if a < b else -1
    return J


def _sdp_cone_proj(X):
    v, u = np.linalg.eigh((X + X.T) / 2)
    pos = v > 0
    return (u[:, pos] * v[pos]) @ u[:, pos].T


def _cycle_weight_opt_admm(J, cons_ids, eps2, mu_init=1.0, mu_rho=1.006, num_iters=600):
    # remove the null space of the edge indicator space
    A = J @ J.T
    v, u = np.linalg.eigh((A + A.T) / 2)
    J = u[:, v > 1e-5].T @ J

    dim, dim_w = J.shape
    w = np.ones(dim_w) / dim_w
    X1 = (J * w[None, :]) @ J.T
    X2 = X1.copy()
    ev = np.linalg.eigvalsh((X1 + X1.T) / 2)
    s1, s2 = ev.min(), ev.max()
    y = 0.0
    Y1 = np.zeros((dim, dim))
    Y2 = np.zeros((dim, dim))
    z = np.zeros(dim_w)
    mu = mu_init

    dd = dim * dim
    dids = np.flatnonzero(np.eye(dim).reshape(-1, order="F"))
    J2 = np.kron(J, np.ones((dim, 1))) * np.kron(np.ones((dim, 1)), J)
    L = np.eye(dd)
    L[np.ix_(dids, dids)] += 1
    L = np.kron(np.eye(2), L)
    JJ = np.vstack([J2, -J2])
    L = L + JJ @ JJ.T
    invL = np.linalg.inv(L)

    vec = lambda M: M.reshape(-1, order="F")
    for _ in range(num_iters):
        S1 = _sdp_cone_proj(Y1 - mu * X1)
        S2 = _sdp_cone_proj(Y2 - mu * X2)

        fixed = z - y + mu * w
        b1 = np.zeros(dd)
        b2 = np.zeros(dd)
        b1[dids] += 1 + mu * s1
        b2[dids] += 1 - mu * s2
        b1 += vec(S1 + mu * X1) - J2 @ fixed
        b2 += vec(S2 + mu * X2) + J2 @ fixed
        x = invL @ np.concatenate([b1, b2])
        Y1 = x[:dd].reshape(dim, dim, order="F")
        Y2 = x[dd:].reshape(dim, dim, order="F")

        fixed2 = np.sum(J * ((Y1 - Y2) @ J), axis=0)
        z = -fixed2 + y - mu * w
        z[cons_ids] += eps2 * mu
        z = np.maximum(0, z)

        y = (np.sum(z + fixed2 + mu * w) - mu) / dim_w

        s1 = s1 + (1 - np.trace(Y1)) / mu
        s2 = s2 + (np.trace(Y2) - 1) / mu
        X1 = X1 + (S1 - Y1) / mu
        X2 = X2 + (S2 - Y2) / mu
        w = w + (fixed2 + z - y) / mu
        mu *= mu_rho
    return w


def _sub_sampling(weights, num_cycles, rng):
    for p in np.arange(100, -1, -1) / 100.0:
        tp = np.power(weights / weights.max(), p)
        if tp.sum() >= num_cycles:
            prob = tp * num_cycles / tp.sum()
            adjusted = np.power(weights, 1 - p)
            ids = np.flatnonzero(rng.random(len(tp)) < prob)
            new_w = adjusted[ids]
            return ids, new_w * weights.sum() / new_w.sum()
    return np.arange(len(weights)), weights


def cycle_basis_generator(A, num_cycles, rng=None):
    """Returns a list of dicts {'cycle': [v0, v1, ...], 'weight': w}."""
    rng = np.random.default_rng(0) if rng is None else rng
    A = np.asarray(A)
    nv = A.shape[0]
    nbrs = [np.flatnonzero(A[v]) for v in range(nv)]
    trees = [_bfs_tree(nbrs, v) for v in range(nv)]
    depths = np.array([t[1] for t in trees])
    edge_mat, edges = adj_mat_2_edge(A)
    root = int(np.argmin(depths))

    num_e = edges.shape[1]
    per_root = num_e - nv + 1
    cand, cons_ids = [], None
    for v in range(nv):
        if v == root:
            cons_ids = np.arange(len(cand), len(cand) + per_root)
        cand += _cycles_from_tree(edges, trees[v][0])

    J = _jacobi_matrix(cand, edge_mat, num_e)
    eps2 = 1 / 4 / per_root
    w = _cycle_weight_opt_admm(J, cons_ids, eps2)

    cycles = [{"cycle": cand[c], "weight": w[c]} for c in cons_ids]
    remaining = np.setdiff1d(np.arange(len(cand)), cons_ids)
    ids, new_w = _sub_sampling(w[remaining], num_cycles - len(cons_ids), rng)
    cycles += [{"cycle": cand[remaining[k]], "weight": nw} for k, nw in zip(ids, new_w)]
    return cycles
