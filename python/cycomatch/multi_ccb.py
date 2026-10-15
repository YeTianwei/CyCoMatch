"""Main pipeline (Multi_CCB.m)."""
import numpy as np

from .admm import admm_optimization
from .cycles import cycle_basis_generator
from .fmaps import add_fmap, fmap2pmap_group, lmd_filter, mwp_reconstruct


def multi_ccb(shapes, T0, G, thresholds=(0.12, 0.1), Nf=6, num_iters=7, alpha=0.3,
              beta=1.0, mu=1.0, rho=1e-3, max_iter=100, tol=1e-4,
              matlab_compat=False, seed=0, verbose=True):
    n = len(shapes)
    num_cycles = n * (n - 1) * (n - 2) // 2
    cycles = cycle_basis_generator(G, num_cycles, np.random.default_rng(seed))

    C_opt, T = None, T0
    for it, thr in enumerate(thresholds):
        log = print if verbose else (lambda *a, **k: None)
        log(f"Iteration {it + 1}")
        log("Identify landmarks...")
        TMap, dim, C1 = lmd_filter(shapes, T, thr, matlab_compat)
        log("Reconstruct functional maps...")
        C2 = mwp_reconstruct(shapes, TMap, Nf, num_iters, dim)
        C = add_fmap(C1, C2, alpha)
        log("Cycle Consistency Refinement...")
        C_opt = admm_optimization(C, beta, mu, G, cycles, rho, max_iter, tol)
        T = fmap2pmap_group(shapes, C_opt)
    return C_opt, T
