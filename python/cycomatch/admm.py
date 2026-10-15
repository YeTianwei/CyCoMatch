"""Cycle-consistency refinement of the functional maps (admm_optimization.m)."""
import numpy as np

from .cycles import adj_mat_2_edge


def _directed_edges(G):
    _, edges = adj_mat_2_edge(G)
    return [(int(i), int(j)) for i, j in edges.T]


def _cycle_terms(cycles):
    """For every directed pair (i, j) the (z, weight) of the last cycle writing Z2{i,j}.

    The MATLAB loop overwrites Z2{i,j} for each cycle containing i -> j, so only
    the last one counts; precomputing it avoids redundant 500x500 products.
    """
    last = {}
    for c in cycles:
        cyc, L = c["cycle"], len(c["cycle"])
        for k in range(L):
            last[(cyc[k], cyc[(k + 1) % L])] = (cyc[(k + 2) % L], c["weight"])
    return last, sum(c["weight"] for c in cycles)


def admm_optimization(C0, beta, mu, G, cycles, rho, max_iter, tol, eps=1e-6, verbose=False):
    n_sh = len(C0)
    k = C0[0][0].shape[0]
    I = np.eye(k)
    edges = _directed_edges(G)
    last, sum_w = _cycle_terms(cycles)

    C = [row[:] for row in C0]
    Z = [row[:] for row in C0]
    U = [[np.zeros_like(c) for c in row] for row in C0]

    for it in range(max_iter):
        C_old, Z_old = C, Z

        # C-update
        C = [row[:] for row in C_old]
        for i, j in edges:
            for a, b in ((i, j), (j, i)):
                C[a][b] = (C0[a][b] + rho * (Z[a][b] - U[a][b])) / (1 + rho)

        # Z-update
        Z = [row[:] for row in Z_old]
        for i, j in edges:
            for a, b in ((i, j), (j, i)):
                Z1 = C[a][b] + U[a][b]
                if (a, b) in last:
                    z, w = last[(a, b)]
                    Z2 = beta * w * (Z_old[a][b] @ Z_old[b][z] - I)
                else:
                    Z2 = 0
                Zab = Z_old[a][b]
                Z3 = mu * Zab @ (Zab.T @ Zab - I)
                if (a, b) == (i, j):
                    Z3 = Z3 + eps * I
                Z[a][b] = (Z1 + Z2 + Z3) / (sum_w + beta + mu)

        # dual update
        U = [row[:] for row in U]
        for i, j in edges:
            for a, b in ((i, j), (j, i)):
                U[a][b] = U[a][b] + C[a][b] - Z[a][b]

        dC = sum(np.linalg.norm(C[a][b] - C_old[a][b]) for a in range(n_sh) for b in range(n_sh))
        dZ = sum(np.linalg.norm(Z[a][b] - Z_old[a][b]) for a in range(n_sh) for b in range(n_sh))
        if verbose:
            print(f"  ADMM iter {it + 1}: dC={dC:.3e} dZ={dZ:.3e}")
        if dC < tol and dZ < tol:
            break
    return C
